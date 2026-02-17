!!TODO_ACC: remove `V_ts` from transfers...

MODULE icedyn_rhg_bbm
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_bbm  ***
   !!   Sea-Ice dynamics : rheology Britle Maxwell X
   !!======================================================================
   !! History :
   !!            4.2  !  2022     (L. Brodeau) `BBM`
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_bbm : computes ice velocities from BBM rheology
   !!   rhg_bbm_rst     : read/write BBM fields in ice restart
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE sbc_oce , ONLY : nn_fsbc
   USE oss_nnq , ONLY : ssh_m, ln_ice_embd
   USE sbc_ice , ONLY : utau_ice, vtau_ice, snwice_mass_b
   USE par_ice
   USE ice            ! sea-ice: ice variables
   USE icevar,   ONLY : ice_var_sshdyn
   USE bdy ,     ONLY : ln_bdy
   USE bdyice
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#if defined _OPENACC
   USE lbclnk_gpu
#endif
   USE prtctl         ! Print control
   !
   USE timing

   !USE remap_weno   , ONLY: rmpT2U_wn5s, rmpT2V_wn5s
   USE remap_classic, ONLY: rmpT2F, rmpT2U, rmpT2V, rmpU2V, rmpV2U, do_rmpT2F

   USE ice_util

   USE icedyn_rhg_tools

   USE icedyn_rhg_vel

   !USE iceistate , ONLY : ln_iceini


   !USE io_ezcdf , ONLY : DUMP_FIELD ; !#LOLOdebug

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_bbm_init  ! called by icedyn_rhg.F90
   PUBLIC   ice_dyn_rhg_bbm       ! called by icedyn_rhg.F90
   PUBLIC   rhg_bbm_rst           ! called by icedyn_rhg.F90

   REAL(wp), SAVE :: rk0  ! factor to stiffness matrix => 1._wp / ( 1._wp - rnup*rnup)
   REAL(wp), SAVE :: rk11, rk22, rk12, rk33 ! elements of stiffness matrix
   REAL(wp), SAVE :: rlambda0, rsqrt_E0 ! Constant part of Eq.28
   REAL(wp), SAVE :: rdtbbm, r1_dtbbm !: small time step (time splitting) [s] and its inverse [1/s]
   !$acc declare create( rk0, rk11, rk22, rk12, rk33, rlambda0, rsqrt_E0, rdtbbm, r1_dtbbm )


   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   kmsk01x, kmsk01y                ! dummy arrays
   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   kmsk00x, kmsk00y                ! mask for ice presence
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ztmp1, ztmp2, ztmp3, ztmp4, zht, zhf, zxpCt, zxpCf
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zPmax_t, zPmax_f, zxpCtbet, zxpCfbet
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zSclH_t, zSclH_f
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   zmU_dt, zmV_dt                  ! (ice-snow_mass / dt) on U/V points
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zgrdSH                          ! surface pressure gradient at U/V points




   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icedyn_rhg_bbm.F90 13646 2020-10-20 15:33:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE ice_dyn_rhg_bbm( kt, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_bbm  ***
      !!                             BBM-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  the BBM rheology of Olason et al., 2022.
      !!
      !! ** Inputs  : - wind forcing (stress), oceanic currents
      !!                ice total volume (vt_i) per unit area
      !!                snow total volume (vt_s) per unit area
      !!
      !! ** Action  : - compute u_ice, v_ice : the components of the
      !!                sea-ice velocity vector
      !!              - compute delta_i, shear_i, divu_i, which are inputs
      !!                of the ice thickness distribution
      !!
      !! ** Steps   : 0) compute mask at F point
      !!              1) Compute ice snow mass, ice strength
      !!              2) Compute wind, oceanic stresses, mass terms and
      !!                 coriolis terms of the momentum equation
      !!              3) Solve the momentum equation (iterative procedure)
      !!              4) Recompute delta, shear and divergence
      !!                 (which are inputs of the ITD) & store stress
      !!                 for the next time step
      !!              5) Diagnostics including charge ellipse
      !!
      !! ** Notes   :
      !!
      !!
      !!
      !!
      !! References : Olason et al., 2022 #fixme
      !!              Hunke and Dukowicz, JPO97
      !!              Bouillon et al., Ocean Modelling 2009
      !!              Bouillon et al., Ocean Modelling 2013
      !!              Kimmritz et al., Ocean Modelling 2016 & 2017
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in ) :: kt                                    ! time step
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pstress1_i, pstress2_i, pstress12_i   !
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pshear_i  , pdivu_i   , pdelta_i      !
      !!-------------------------------------------------------------------
      !
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   jter         ! local integers
      !
      REAL(wp) ::   zmassU, zmassV                       ! ice/snow mass and volume
      !
      REAL(wp) :: zr, zr1, zr2, zr3, zmsk, zzt, zzf

      REAL(wp) :: zravrg

      LOGICAL, PARAMETER :: l_apply_d_healing = .TRUE. !LOLO!

      !CHARACTER(len=64) :: cf_tmp ; !#LOLOdebug

      !!-------------------------------------------------------------------
      !$acc data present( pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i, u_ice, v_ice, uVice, vUice, SIGMAt, SIGMAf, dmdt, dmdf )

      IF( ln_timing )   CALL timing_start('ice_dyn_rhg_bbm')

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_rhg_bbm: BBM sea-ice rheology'

      !------------------------------------------------------------------------------!
      ! 1) define some variables and initialize arrays
      !------------------------------------------------------------------------------!

      !! Ice thickness @T & @F we are going to work with:
      IF( l_use_v_for_h ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zht(ji,jj) = MAX( vt_i(ji,jj) , 0._wp)
            END DO
         END DO
         !$acc end parallel loop
         CALL do_rmpT2F( zht, zhf, lconserv=.TRUE. )
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zhf(ji,jj) = MAX( zhf(ji,jj) , 0._wp)
            END DO
         END DO
         !$acc end parallel loop
# if ! defined _OPENACC
         CALL lbc_lnk( 'icedyn_rhg_bbm', zhf,'F',1._wp )
# endif
         !
      ELSE
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zht(ji,jj) = hm_i  (ji,jj)
               zhf(ji,jj) = hm_i_f(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF !IF( l_use_v_for_h )
      !


      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      ! `zxpCt` & `zxpCf` are used as temporary arrays for SSH@T & SSH@F, respectively !
      CALL ice_var_sshdyn( ssh_m, snwice_mass_b, zxpCt )

      CALL do_rmpT2F( zxpCt, zxpCf )  ! `zxpCf` -> SSH at F-points

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            kmsk00x(ji,jj)  = 1
            kmsk00y(ji,jj)  = 1
            kmsk01x(ji,jj)  = 1
            kmsk01y(ji,jj)  = 1
            zgrdSH(ji,jj,1) = 0._wp
            zgrdSH(ji,jj,2) = 0._wp
            zgrdSH(ji,jj,3) = 0._wp
            zgrdSH(ji,jj,4) = 0._wp
            !
            ! Ice/snow mass (kg/m^2):
            zr  = ( rhos*vt_s(ji  ,jj  ) + rhoi*vt_i(ji  ,jj  ) ) * e1e2t(ji  ,jj  )  ! mass @T [kg]
            zr1 = ( rhos*vt_s(ji+1,jj  ) + rhoi*vt_i(ji+1,jj  ) ) * e1e2t(ji+1,jj  )
            zr2 = ( rhos*vt_s(ji  ,jj+1) + rhoi*vt_i(ji  ,jj+1) ) * e1e2t(ji  ,jj+1)
            ztmp3(ji,jj) = 0.5_wp*( zr + zr1 )*r1_e1e2u(ji,jj)*umask(ji,jj,1) ! ztmp3 -> mass @ U [kg/m^2]
            ztmp4(ji,jj) = 0.5_wp*( zr + zr2 )*r1_e1e2v(ji,jj)*vmask(ji,jj,1) ! ztmp4 -> mass @ V [kg/m^2]
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2) present(umask,vmask)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zmassU = ztmp3(ji,jj)
            zmassV = ztmp4(ji,jj)

            ! m/dt
            zmU_dt(ji,jj)   = zmassU * r1_dtbbm
            zmV_dt(ji,jj)   = zmassV * r1_dtbbm

            ! Surface pressure gradient (- m*g*GRAD(ssh)) at U-V points
            zgrdSH(ji,jj,1) = - zmassU * grav * ( zxpCt(ji+1,jj) - zxpCt(ji,jj)   ) * r1_e1u(ji,jj)
            zgrdSH(ji,jj,2) = - zmassV * grav * ( zxpCt(ji,jj+1) - zxpCt(ji,jj)   ) * r1_e2v(ji,jj)
            zgrdSH(ji,jj,3) = - zmassV * grav * ( zxpCf(ji,jj)   - zxpCf(ji-1,jj) ) * r1_e1v(ji,jj)  ! `zxpCf` is `zxpCt` interpolated  @F !
            zgrdSH(ji,jj,4) = - zmassU * grav * ( zxpCf(ji,jj)   - zxpCf(ji,jj-1) ) * r1_e2u(ji,jj)  ! `zxpCf` is `zxpCt` interpolated  @F !

            ! masks
            IF( zmassU<=0._wp .OR. umask(ji,jj,1)<0.1_wp ) kmsk00x(ji,jj) = 0
            IF( zmassV<=0._wp .OR. vmask(ji,jj,1)<0.1_wp ) kmsk00y(ji,jj) = 0

            ! switches
            IF( zmassU <= rmass_min .AND. au_i(ji,jj) <= rconc_min )  kmsk01x(ji,jj) = 0
            IF( zmassV <= rmass_min .AND. av_i(ji,jj) <= rconc_min )  kmsk01y(ji,jj) = 0

         END DO
      ENDDO
      !$acc end parallel loop


      ! --- Healing of damage with time [Eq.30, Olason et al., 2022]
      !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF( l_apply_d_healing ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               zmsk = xmskt(ji,jj)
               zr1 = 1._wp / MAX(at_i(ji,jj),epsi06)
               zr3 = rcnd_i*vt_s(ji,jj)*zr1 / MAX( rcnd_s*vt_i(ji,jj)*zr1, epsi06 ) * zmsk   ! => `C` of the
               IF(ln_icethd) THEN ! Thermo is on, normal stuff
                  ztmp1(ji,jj) = (t_bo(ji,jj) - tm_su(ji,jj)) / (1._wp + zr3 ) * zmsk ! temp. difference between bottom and surface
               ELSE
                  !! Thermo is off, yet we want som refreezing!
                  ztmp1(ji,jj) = (-1.8_wp + 25._wp) / (1._wp + zr3 ) * zmsk ! faked temp. difference between bottom and surface
               END IF
               !
               ztmp2(ji,jj) = rdt_ice * MIN( ztmp1(ji,jj) / rn_kth , 1._wp/rdt_ice )  ! dt * 1/T_relax => `1-d` increment @ T
            END DO
         END DO
         !$acc end parallel loop
         CALL do_rmpT2F( ztmp1, ztmp3 )
         CALL do_rmpT2F( ztmp2, ztmp4 )
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ! Apply healing on damage:
               IF( ztmp1(ji,jj) > 0._wp )  dmdt(ji,jj) = dmdt(ji,jj) + ztmp2(ji,jj)
               IF( ztmp3(ji,jj) > 0._wp )  dmdf(ji,jj) = dmdf(ji,jj) + ztmp4(ji,jj)
            END DO
         END DO
         !
      ENDIF !IF( l_apply_d_healing )

      CALL cap_1md( at_i, af_i, dmdt, dmdf ) ! Capping for both  post-healing ("post-advection" is done in `icedyn_adv`):


      zravrg = 1._wp/REAL(nbbm) ! going to average (set to 0 before accumulating during the `nbbm` sub time steps)

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            u_ice(ji,jj) = 0._wp
            v_ice(ji,jj) = 0._wp
            uVice(ji,jj) = 0._wp
            vUice(ji,jj) = 0._wp
            !
            ! We do not want to do the following stuff `nbbm` times below because they remain unchanged:
            zzt = EXP( rn_C0*(1._wp - at_i(ji,jj)) )
            zzf = EXP( rn_C0*(1._wp - af_i(ji,jj)) )
            zxpCt(ji,jj) = zzt
            zxpCf(ji,jj) = zzf
            zxpCtbet(ji,jj) = zzt ** nn_btrlx
            zxpCfbet(ji,jj) = zzf ** nn_btrlx
            !
            zPmax_t(ji,jj)  =  -rn_P0 * zzt * zht(ji,jj) ** 2.5_wp ! `2.5` because working with vertically integrated sigmas => `h^3/2 * h`
            zPmax_f(ji,jj)  =  -rn_P0 * zzf * zhf(ji,jj) ** 2.5_wp !   "               "               "             "
            !
            zSclH_t(ji,jj)  = SQRT( rn_l_ref / REAL(res_grd_loc_t(ji,jj),wp) ) * zht(ji,jj)  ! required for Mohr-Coulomn test (! multiply with `h` "  " ")
            zSclH_f(ji,jj)  = SQRT( rn_l_ref / REAL(res_grd_loc_f(ji,jj),wp) ) * zhf(ji,jj)  ! required for Mohr-Coulomn test (! multiply with `h` "  " ")
         ENDDO
      ENDDO
      !$acc end parallel loop


      !$acc loop seq                                  ! ==================== !
      DO jter = 1 , nbbm                           !    loop over jter    !
         !                                            ! ==================== !

         ! ---  Updates the components of the vertically-integrated internal stress tensor and the damage in both T- & F-centric worlds ---
         !           => based on previously computed ice velocities...

         CALL update_sigma_d( kt, jter, rdtbbm, V_ts, at_i, af_i, zxpCt, zxpCf, zht, zhf, SIGMAt, SIGMAf, dmdt, dmdf )
         !           => `sigmas` & `d` ok on whole `Nis0-1:Nie0+1,Njs0-1:Nje0+1` ! (provided `V_ts` was fully lbclinked!)

         ! --- Computation of ice velocity --- ! (`Nis0:Nie0,Njs0:Nje0`)
         CALL update_uv_euler_si( jter, rdtbbm, au_i, av_i, zmU_dt, zmV_dt, SIGMAt, SIGMAf, zgrdSH, V_oce, utau_ice, vtau_ice, &
            &                     kmsk01x, kmsk01y, kmsk00x, kmsk00y,  V_ts )

         !CALL update_uv_rk3( jter, rdtbbm, au_i, av_i, zmU_dt, zmV_dt, SIGMAt, SIGMAf, zgrdSH, V_oce, Tau_ai, kmsk01x, kmsk01y, kmsk00x, kmsk00y,  V_ts )
         !
#if defined _OPENACC
         IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', V_ts )
         IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', V_ts )
#else
         CALL lbc_lnk( 'icedyn_rhg_bbm',  V_ts(:,:,1),'U',-1._wp, V_ts(:,:,2),'V',-1._wp,  V_ts(:,:,3),'V',-1._wp, V_ts(:,:,4),'U',-1._wp )
#endif

         IF( ln_bdy ) THEN
            CALL bdy_ice_dyn( 'U', V_ts(:,:,1) )
            CALL bdy_ice_dyn( 'V', V_ts(:,:,2) )
            CALL bdy_ice_dyn( 'V', V_ts(:,:,3), l_FcVel=.TRUE. )
            CALL bdy_ice_dyn( 'U', V_ts(:,:,4), l_FcVel=.TRUE. )
            CALL bdy_ice_dmg( kt, jter, V_ts )
         ENDIF

         ! Average the velocity (to be used for advection at the big time step):
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               u_ice(ji,jj) = u_ice(ji,jj) + zravrg*V_ts(ji,jj,1)
               v_ice(ji,jj) = v_ice(ji,jj) + zravrg*V_ts(ji,jj,2)
               uVice(ji,jj) = uVice(ji,jj) + zravrg*V_ts(ji,jj,3)
               vUice(ji,jj) = vUice(ji,jj) + zravrg*V_ts(ji,jj,4)
            ENDDO
         ENDDO
         !$acc end parallel loop
         !
         !                                             ! ==================== !
      END DO !DO jter = 1 , nbbm                       !  end loop over jter  !
      !                                                ! ==================== !


      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!
      CALL strain_rate_dsd( 'T', u_ice, v_ice, uVice, vUice, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e1t2, e2t2, xmskt, &
         &                       pdivu_i, pshear_i, pdelta_i )
# if ! defined _OPENACC
      CALL lbc_lnk( 'icedyn_rhg_bbm', pdivu_i,'T',1._wp, pshear_i,'T',1._wp, pdelta_i,'T',1._wp )
# endif

      !! sigma_1, sigma_2 & sigma_12 in Pa.m:
      !$acc parallel loop collapse(2) present( pstress1_i, pstress2_i, pstress12_i )
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pstress1_i (ji,jj) = ( SIGMAt(ji,jj,1) + SIGMAt(ji,jj,2) ) ! @T
            pstress2_i (ji,jj) = ( SIGMAt(ji,jj,1) - SIGMAt(ji,jj,2) ) ! @T
            pstress12_i(ji,jj) =          SIGMAt(ji,jj,3)              ! @F
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_dyn_rhg_bbm')
      !
   END SUBROUTINE ice_dyn_rhg_bbm


   SUBROUTINE ice_dyn_rhg_bbm_init( )
      !!-------------------------------------------------------------------
      !! Called into `ice_dyn_rhg_init()@icedyn_rhg.F90`
      !!-------------------------------------------------------------------
      INTEGER  ::   ierr
      INTEGER  ::   icycle
      REAL(wp) ::   ztmp, zdx_m, zce, zdts
      REAL(wp), DIMENSION(jpi,jpj) :: zt1, zt2
      !!-------------------------------------------------------------------
      IF( lwp ) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) '    *** Initialization of BBM rheology (ice_dyn_rhg_bbm_init) ***'
      ENDIF

      !! Stiffness matrix
      !! ****************
      rk0  = 1._wp / ( 1._wp - rnup*rnup)
      rk11 = rk0
      rk12 = rk0 * rnup
      rk22 = rk0
      rk33 = rk0 * (1._wp - rnup)

      rlambda0 = rn_eta0 / rn_E0     !: Viscosity / Elasticity of undamaged ice (aka relaxation time) [s]
      IF( lwp ) WRITE(numout,*) '  * Viscous relaxation time scale => rlambda0 =', REAL(rlambda0,4), ' [s]'

      rsqrt_E0      = SQRT( rn_E0 )

      ! Find the smallest `dx` of the whole WET domain:
      zt1(:,:) = REAL( res_grd_loc_t(:,:) , 4 ) ! SQRT(dx*dy)
      zt1(:,:) = MERGE( zt1(:,:) , 1.E12_wp , (xmskt(:,:) > 0.9_wp) )  ! => stupidly big value over continents...
      zdx_m = MINVAL( zt1 )                        ! min of dx local
      CALL mpp_min( 'ice_dyn_rhg_bbm_init', zdx_m) ! min of dx over the whole domain

      ! time step adjusted automatically
      ! ********************************
      ! we look at 2) propagation speed of elastic waves (zce), 2) the time step to solve them (zdts)
      !            3) the number of iterations needed (icycle) to go from t to t+dt
      !            4) the number of iterations over which u_ice is averaged (nflt)
      !       then 5) the number of total iterations of the rheology
      zce = rsqrt_E0 / rsqrt_nu_rhoi ! propagation speed of shearing elastic waves based on the mean `dx`
      zdts = 0.5_wp*(zdx_m/zce)      ! largest possible small time-step to consider... #LOLO: Clem uses `0.9` rather than `0.5`
      icycle = CEILING( rDt_ice / zdts )        ! number of iterations (local to init, cf below for nbbm)
      icycle = icycle + MOD(icycle, 2)          ! + make it an odd number
      rdtbbm   = rDt_ice / REAL( icycle, wp )   ! small time step used here
      r1_dtbbm = 1._wp / rdtbbm
      !
      nflt = FLOOR( rn_bbm_flt * icycle ) ! number of iterations over which u_ice is averaged
      nflt = nflt + MOD( nflt, 2)         ! + make it an odd number
      nbbm = icycle + nflt/2              ! total number of iterations

      CALL cross_nudging_init()

      IF( lwp ) THEN
         WRITE(numout,*) '  * Big time step (advection & thermo)  => rdt_ice  =', rdt_ice, ' [s]'
         WRITE(numout,*) '  * Min `dx` of wet computational domain = ', REAL(zdx_m/1000._wp,4), ' [km]'
         WRITE(numout,*) '     ==> propagation speed of shearing elastic waves =>',  zce, '[m/s]'
         WRITE(numout,*) '     ==> time-step requirement to resolve these waves: dt = 0.5*(dx/c_e) =', zdts, ' [s]'
         WRITE(numout,*) '     ==> implies a `nbbm` =', INT(nbbm,2)
         WRITE(numout,*) '     ==> implies a `nflt` =', INT(nflt,2)
         WRITE(numout,*) '     ==> implies a small time step (rheology) => rdtbbm =', rdtbbm,  ' [s]'
         IF(ln_x_MC_test) WRITE(numout,*) '  * Will perform only 1 Mohr-Coulomb test, at mid-point between T & F points!'
         IF(l_CN) THEN
            WRITE(numout,*) '  * About cross-nudging:'
            WRITE(numout,*) '      - CN parameter (gamma) => rn_crndg =', rn_crndg,' [-]'
         ENDIF
         WRITE(numout,*) '  * (scaled) Compression threshod => N_lim =',REAL(rn_Nref*SQRT( rn_l_ref/zdx_m ),4),' [Pa]'
         WRITE(numout,*) ''
      ENDIF

      IF(rdtbbm>zdts) CALL ctl_warn( 'ice_dyn_rhg_bbm_init: `nbbm` is probably to small' )

      IF( lwp ) THEN
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) ''
      ENDIF

      ALLOCATE( kmsk01x(jpi,jpj), kmsk01y(jpi,jpj),  kmsk00x(jpi,jpj), kmsk00y(jpi,jpj), zmU_dt(jpi,jpj), zmV_dt(jpi,jpj), &
         &        ztmp1(jpi,jpj), ztmp2(jpi,jpj),      ztmp3(jpi,jpj), ztmp4(jpi,jpj), zht(jpi,jpj), zhf(jpi,jpj),         &
         &        zxpCt(jpi,jpj), zxpCf(jpi,jpj),   zgrdSH(jpi,jpj,4), zPmax_t(jpi,jpj), zPmax_f(jpi,jpj),                 &
         &     zxpCtbet(jpi,jpj), zxpCfbet(jpi,jpj), zSclH_t(jpi,jpj), zSclH_f(jpi,jpj),      STAT = ierr )

      CALL mpp_sum( 'ice_dyn_rhg_bbm_init', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ice_dyn_rhg_bbm_init : unable to allocate work arrays')

# if defined _OPENACC
      !$acc update device ( nbbm, nflt, rk0, rk11, rk22, rk12, rk33, rlambda0, rsqrt_E0, rdtbbm, r1_dtbbm )
      PRINT *, ' * info GPU: ice_dyn_rhg_bbm_init() => adding work arrays to memory!'
      !$acc enter data copyin( kmsk01x, kmsk01y, kmsk00x, kmsk00y, zmU_dt, zmV_dt, ztmp1, ztmp2, ztmp3, ztmp4, zht, zhf )
      PRINT *, '    ==> kmsk01x, kmsk01y, kmsk00x, kmsk00y, zmU_dt, zmV_dt, ztmp1, ztmp2, ztmp3, ztmp4, zht, zhf'
      !$acc enter data copyin( zxpCt, zxpCf, zgrdSH, zPmax_t, zPmax_f, zxpCtbet, zxpCfbet, zSclH_t, zSclH_f )
      PRINT *, '    ==> zxpCt, zxpCf, zgrdSH, zPmax_t, zPmax_f, zxpCtbet, zxpCfbet, zSclH_t, zSclH_f'
# endif

   END SUBROUTINE ice_dyn_rhg_bbm_init


   SUBROUTINE update_sigma_d( kt, kts, pdt, pV4, pAt, pAf, pxpCt, pxpCf, pht, phf,  psgmt, psgmf, p1mdt, p1mdf )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE UPDATE_SIGMA_D  ***
      !! ** Purpose :
      !!
      !! ** Method  :
      !!
      !! ** Note    : Called at the sub-time-stepping level!
      !!
      !! ** Author : L. Brodeau, 2022
      !!----------------------------------------------------------------------
      INTEGER,                        INTENT(in)    :: kt, kts        ! # of current big and small/sub-time step
      REAL(wp),                       INTENT(in)    :: pdt             ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in)    :: pV4             ! the 4 components of sea-ice velocity
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pAt, pAf        ! Ice concentration @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf    ! EXP( rn_C0*(1 - A) ) @T & @F !lili
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pht, phf        ! Ice thickness @T & @F
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt           ! T-centric vertically-integrated stress tensor [N/m^2*m]
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf           ! F-centric vertically-integrated stress tensor [N/m^2*m]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf    ! `1 - ice damage` @T & @F
      !!----------------------------------------------------------------------
      REAL(wp)                       :: zfc, zml, zmsk
      REAL(wp)                       :: zh, zE, zeta, zL, zang, zc0, zmul
      REAL(wp)                       :: zE1, zE2
      REAL(wp)                       :: ze11t, ze22t, ze12t, ze11f, ze22f, ze12f
      REAL(wp)                       :: zxpC, z1md, zsigI, zsigII, zPmax, zPtld
      REAL(wp)                       :: z1_zsigI
      INTEGER :: ji, jj, km, kp
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_sigma_d')
      !$acc data present( pV4, pxpCt, pxpCf, pht, phf, psgmt, psgmf, p1mdt, p1mdf )

      IF( (.NOT. l_CN).AND.(kt==nit000).AND.(lwp) ) WRITE(numout,*) ' *** MIND: no cross-nudging will be applied between stress tensors!'

      !LOLOfix:
      !IF( l_CN ) THEN
      kp = 0
      km = 0
      !ELSE
      !   kp = MIN( nn_hls, 2 )
      !   km = kp - 1
      !ENDIF
      !LOLOfix.

      !$acc parallel loop collapse(2)
      DO jj=Njs0-km, Nje0+kp
         DO ji=Nis0-km, Nie0+kp

            zh = pht(ji,jj)

            ! --- Strain rate tensors ---
            !     *******************
#           include "icedyn_rhg_bbm_strn_t.h90"
            !        => uses `ji-1` & `jj-1`
            ! ==> ze11t, ze22t, ze12t

            ! --- E, Lambda & multiplicator (uses `sigmas`!!!) ---
            !     ********************************************
#           include "icedyn_rhg_bbm_elm_t.h90"
            ! ==> zEt, zLt, zmult

            ! --- Predictor estimate of stress tensor at k+1 ---
            !     ******************************************
            zml = zmul * xmskt(ji,jj)
            zfc =  zh * zE * pdt
            psgmt(ji,jj,1) = zml * ( zfc * ( rk11*ze11t + rk12*ze22t) + psgmt(ji,jj,1) )
            psgmt(ji,jj,2) = zml * ( zfc * ( rk12*ze11t + rk22*ze22t) + psgmt(ji,jj,2) )
            psgmf(ji,jj,3) = zml * ( zfc *        rk33 * ze12t        + psgmf(ji,jj,3) )
            !
         END DO !DO ji=Nis0-1, Nie0+2
      END DO !DO jj=Njs0-1, Nje0+2
      !$acc end parallel loop


      !$acc parallel loop collapse(2)
      DO jj=Njs0-kp, Nje0+km
         DO ji=Nis0-kp, Nie0+km

            zh = phf(ji,jj)

            ! --- Strain rate tensors ---
            !     *******************
#           include "icedyn_rhg_bbm_strn_f.h90"
            !        => uses `ji+1` & `jj+1`
            ! ==> ze11f, ze22f, ze12f

            ! --- E, Lambda & multiplicator (uses `sigmas`!!!) ---
            !     ********************************************
#           include "icedyn_rhg_bbm_elm_f.h90"
            ! ==> zEf, zLf, zmulf

            ! --- Predictor estimate of stress tensor at k+1 ---
            !     ******************************************
            zml = zmul * xmskf(ji,jj)
            zfc = zh * zE * pdt
            psgmf(ji,jj,1) = zml * ( zfc * ( rk11*ze11f + rk12*ze22f) + psgmf(ji,jj,1) )
            psgmf(ji,jj,2) = zml * ( zfc * ( rk12*ze11f + rk22*ze22f) + psgmf(ji,jj,2) )
            psgmt(ji,jj,3) = zml * ( zfc *        rk33 * ze12f        + psgmt(ji,jj,3) )
            !
         END DO !DO ji=Nis0-kp, Nie0+km
      END DO !DO jj=Njs0-kp, Nje0+km
      !$acc end parallel loop



      !! --- lbc-linking of updated stress tensors ---
#if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', psgmt, psgmf )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', psgmt, psgmf )
#else
      CALL lbc_lnk( 'icedyn_rhg_bbm', psgmt(:,:,1),'T',1._wp, psgmt(:,:,2),'T', 1._wp, psgmt(:,:,3),'F', 1._wp, &
         &                            psgmf(:,:,1),'F',1._wp, psgmf(:,:,2),'F', 1._wp, psgmf(:,:,3),'T', 1._wp  )
#endif

      IF( l_CN ) THEN
         CALL apply_cn_gpu( kts, psgmt, psgmf ) !  `Nis0-1:Nie0+1,Njs0-1:Nje0+1`
         !CALL apply_cn_trd( kts, psgmt, psgmf ) !  `Nis0-1:Nie0+1,Njs0-1:Nje0+1`
         !CALL apply_cn_wn7( kts, psgmt, psgmf ) !  `Nis0-1:Nie0+1,Njs0-1:Nje0+1`
         !      !!          !        => @T uses `ji+1` & `jj+1`  result => ok on whole `Nis0:Nie0,Njs0:Nje0`
         !      !!          !        => @F uses `ji-1` & `jj-1`  result => ok on whole      "          "
      ENDIF

      IF( ln_x_MC_test ) THEN
         CALL mohr_coulomb_dmg_mp( pdt, pxpCt, pxpCf, zSclH_t, zSclH_f, p1mdt, p1mdf, psgmt, psgmf )
      ELSE
         CALL mohr_coulomb_dmg(    pdt, pxpCt, pxpCf, zSclH_t, zSclH_f, p1mdt, p1mdf, psgmt, psgmf )
      ENDIF

      CALL clean_small_a_all( pAt, pAf,  p1mdt, p1mdf, psgmt, psgmf )

      !! --- lbc-linking of updated stress tensors and `1-d` ---
#if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', psgmt, psgmf, p1mdt, p1mdf )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', psgmt, psgmf, p1mdt, p1mdf )
#else
      CALL lbc_lnk( 'icedyn_rhg_bbm', psgmt(:,:,1),'T',1._wp, psgmt(:,:,2),'T', 1._wp, psgmt(:,:,3),'F', 1._wp, p1mdt,'T',1._wp, &
         &                            psgmf(:,:,1),'F',1._wp, psgmf(:,:,2),'F', 1._wp, psgmf(:,:,3),'T', 1._wp, p1mdf,'F',1._wp  )
#endif

      !$acc end data
      IF( ln_timing )   CALL timing_stop('update_sigma_d')
      !
   END SUBROUTINE update_sigma_d



   SUBROUTINE rhg_bbm_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_bbm_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter            ! local integer
      INTEGER  ::   id01, id02
      INTEGER  ::   id1, id2, id3, id4, id5, id6, id7, id8
      INTEGER  ::   id11, id12, id13, id14
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id01 = iom_varid( numrir, 'dmdt' , ldstop = .FALSE. )
            id02 = iom_varid( numrir, 'dmdf' , ldstop = .FALSE. )
            !
            id1 = iom_varid( numrir, 'sgm11t' , ldstop = .FALSE. )
            id2 = iom_varid( numrir, 'sgm22t' , ldstop = .FALSE. )
            id3 = iom_varid( numrir, 'sgm12f' , ldstop = .FALSE. )
            id4 = iom_varid( numrir, 'sgm11f' , ldstop = .FALSE. )
            id5 = iom_varid( numrir, 'sgm22f' , ldstop = .FALSE. )
            id6 = iom_varid( numrir, 'sgm12t' , ldstop = .FALSE. )
            !
            id7 = iom_varid( numrir, 'uVice' , ldstop = .FALSE. )
            id8 = iom_varid( numrir, 'vUice' , ldstop = .FALSE. )
            !
            id11 = iom_varid( numrir, 'Uu_sub' , ldstop = .FALSE. )
            id12 = iom_varid( numrir, 'Uv_sub' , ldstop = .FALSE. )
            id13 = iom_varid( numrir, 'Vv_sub' , ldstop = .FALSE. )
            id14 = iom_varid( numrir, 'Vu_sub' , ldstop = .FALSE. )

            IF( MIN( id01, id02 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'dmdt' , dmdt , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'dmdf' , dmdf , cd_type = 'F' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without rheology, set damage @T and @F to 0'
               dmdt(:,:) = 1._wp
               dmdf(:,:) = 1._wp
            ENDIF

            IF( MIN( id1, id2, id3, id4, id5, id6 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'sgm11t', SIGMAt(:,:,1), cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22t', SIGMAt(:,:,2), cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12f', SIGMAt(:,:,3), cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm11f', SIGMAf(:,:,1), cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22f', SIGMAf(:,:,2), cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12t', SIGMAf(:,:,3), cd_type = 'T' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>> did not find components of stress tensors in restart file => set to 0'
               SIGMAt(:,:,:) =  0._wp
               SIGMAf(:,:,:) =  0._wp
            ENDIF

            IF( MIN( id7, id8 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'uVice' , uVice , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'vUice' , vUice , cd_type = 'U', psgn = -1._wp )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without BBM rheology, interpolate F-centric velocities'
               uVice(:,:) = rmpU2V( u_ice )
               vUice(:,:) = rmpV2U( v_ice )
               CALL lbc_lnk( 'rhg_bbm_rst',  uVice,'V',-1._wp, vUice,'U',-1._wp )
            ENDIF

            IF( MIN( id11, id12, id13, id14 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'Uu_sub' , V_ts(:,:,1) , cd_type = 'U', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Uv_sub' , V_ts(:,:,3) , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Vv_sub' , V_ts(:,:,2) , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Vu_sub' , V_ts(:,:,4) , cd_type = 'U', psgn = -1._wp )
            ELSE
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without BBM rheology, fill sub-ts velocities'
               V_ts(:,:,1) = u_ice(:,:)
               V_ts(:,:,3) = uVice(:,:)
               V_ts(:,:,2) = v_ice(:,:)
               V_ts(:,:,4) = vUice(:,:)
            ENDIF
            !
         ELSE                                   !* Start from rest
            !
            IF(lwp) WRITE(numout,*)
            !
            IF(.NOT. ln_iceini) THEN
               IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set T- and F- centric damage to 0'
               dmdt(:,:) = 1._wp
               dmdf(:,:) = 1._wp
            ELSE
               IF(lwp) WRITE(numout,*) '   ==>>>   damage@T taken from `sn_dmg@namini` file => damage@F interpolated!'
               dmdt(:,:) = MIN( MAX(         dmdt(:,:)                , r_dmd_min ) , 1._wp )
               dmdf(:,:) = MIN( MAX( rmpT2F( dmdt,  lconserv=.TRUE. ) , r_dmd_min ) , 1._wp )
            ENDIF
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set T-centric stresses to 0'
            SIGMAt(:,:,:) = 0._wp
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set F-centric stresses to 0'
            SIGMAf(:,:,:) = 0._wp
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set F-centric velocities to 0'
            uVice(:,:)  = 0._wp
            vUice(:,:)  = 0._wp
            V_ts(:,:,:) = 0._wp
            !
         ENDIF
         !
         ! Update onto the GPU:
         !$acc update device( dmdt, dmdf, u_ice, v_ice, uVice, vUice, V_ts, SIGMAt, SIGMAf )

         !
         !
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         ! Update onto the HOST:
         !$acc update self( dmdt, dmdf, u_ice, v_ice, uVice, vUice, V_ts, SIGMAt, SIGMAf )
         !
         IF(lwp) WRITE(numout,*) '---- rhg-rst ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         CALL iom_rstput( iter, nitrst, numriw, 'dmdt' , dmdt  )
         CALL iom_rstput( iter, nitrst, numriw, 'dmdf' , dmdf  )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sgm11t' , SIGMAt(:,:,1) )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm22t' , SIGMAt(:,:,2) )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm12f' , SIGMAt(:,:,3) )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm11f' , SIGMAf(:,:,1) )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm22f' , SIGMAf(:,:,2) )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm12t' , SIGMAf(:,:,3) )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'uVice' , uVice )
         CALL iom_rstput( iter, nitrst, numriw, 'vUice' , vUice )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'Uu_sub' , V_ts(:,:,1) )
         CALL iom_rstput( iter, nitrst, numriw, 'Uv_sub' , V_ts(:,:,3) )
         CALL iom_rstput( iter, nitrst, numriw, 'Vv_sub' , V_ts(:,:,2) )
         CALL iom_rstput( iter, nitrst, numriw, 'Vu_sub' , V_ts(:,:,4) )
         !
      ENDIF
      !
   END SUBROUTINE rhg_bbm_rst

   !!==============================================================================
END MODULE icedyn_rhg_bbm
