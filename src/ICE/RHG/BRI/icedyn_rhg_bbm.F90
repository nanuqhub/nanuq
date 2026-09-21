MODULE icedyn_rhg_bbm
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_bbm  ***
   !!   Sea-Ice dynamics : rheology Britle Maxwell X
   !!======================================================================
   !! History :
   !!            4.2  !  2022     (L. Brodeau)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_bbm : computes ice velocities from BBM rheology
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE oss_nnq , ONLY : ssh_m, ln_ice_embd
   USE sbc_ice , ONLY : taux_ai_t, tauy_ai_t, snwice_mass_b
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
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu
#else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#endif
   USE prtctl         ! Print control
   !
   USE timing

   !USE remap_weno   , ONLY: rmpT2U_wn5s, rmpT2V_wn5s
   USE remap_classic, ONLY: rmpT2F, rmpT2U, rmpT2V, do_rmpT2F

   USE ice_util

   USE icedyn_rhg_tools, ONLY: sigmaII_sclr, strain_rate_dsd
   USE icedyn_rhg_bri

   USE icedyn_rhg_vel

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_bbm_init  ! called by icedyn_rhg.F90
   PUBLIC   ice_dyn_rhg_bbm       ! called by icedyn_rhg.F90

   REAL(wp), SAVE :: rk0  ! factor to stiffness matrix => 1._wp / ( 1._wp - rnup*rnup)
   REAL(wp), SAVE :: rk11, rk22, rk12, rk33 ! elements of stiffness matrix
   REAL(wp), SAVE :: rsqrt_E0 ! Constant part of Eq.28
   REAL(wp), SAVE :: rdtbri, r1_dtbri !: small time step (time splitting) [s] and its inverse [1/s]
   !$acc declare create( rk0, rk11, rk22, rk12, rk33, rsqrt_E0, rdtbri, r1_dtbri )

   !! Work arrays specific to BBM:
   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   kmsk01x, kmsk01y     ! dummy arrays
   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   kmsk00x, kmsk00y     ! mask for ice presence
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xxpCt, xxpCf
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xPmax_t, xPmax_f
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xScHt, xScHf     ! factors to scale cohesion and `Nlim` to local grid resolution
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   xxmU, xxmV           ! `ln_bri_rk3=F` => `ice-snow_mass / dt` on U/V points
   !                                                                         ! `ln_bri_rk3=T` => `1/ice-snow_mass`
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xgrdH               ! surface pressure gradient at U/V points

   LOGICAL, PARAMETER :: l_apply_d_healing = .TRUE.

   CHARACTER(len=15), PARAMETER :: crtnm = 'ice_dyn_rhg_bbm'

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE ice_dyn_rhg_bbm( kt, pshear_i, pdivu_i, pdelta_i )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_bbm  ***
      !!                             BBM-E-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  the BBM rheology of Òlason et al., 2022.
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
      !! References : Brodeau et al., 2024, GMD
      !!              Olason et al., 2022
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in ) :: kt                                    ! time step
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pshear_i, pdivu_i, pdelta_i      !
      !!-------------------------------------------------------------------
      INTEGER  ::  ji, jj
      INTEGER  ::  jter
      REAL(wp) ::  zmassU, zmassV   ! ice/snow mass and volume
      REAL(wp) ::  zr, zr1, zr2, zr3, zmsk, zzt, zzf, zht, zhf, zravrg
      !!-------------------------------------------------------------------
      !$acc data present( pshear_i,pdivu_i,pdelta_i,u_ice,v_ice,uVice,vUice,SIGMAt,SIGMAf,dmdt,dmdf,SI1t,SI2t,SI1f,SI2f )

      IF( ln_timing )   CALL timing_start(crtnm)

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- '//crtnm//': BBM sea-ice rheology'

      !------------------------------------------------------------------------------!
      ! 1) define some variables and initialize arrays
      !------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      ! `xtmp1` & `xtmp2` are used as temporary arrays for SSH@T & SSH@F, respectively !
      CALL ice_var_sshdyn( ssh_m, snwice_mass_b, xtmp1 )

      CALL do_rmpT2F( xtmp1, xtmp2 )  ! `xtmp2` -> SSH at F-points

      !$acc parallel loop collapse(2) present(umask,vmask,xtmp1,xtmp2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ! Ice/snow mass (kg/m^2):
            zr  = ( rhos*vt_s(ji  ,jj  ) + rhoi*vt_i(ji  ,jj  ) ) * e1e2t(ji  ,jj  )  ! mass @T [kg/m^2] (`vt_*` is in m!)
            zr1 = ( rhos*vt_s(ji+1,jj  ) + rhoi*vt_i(ji+1,jj  ) ) * e1e2t(ji+1,jj  )
            zr2 = ( rhos*vt_s(ji  ,jj+1) + rhoi*vt_i(ji  ,jj+1) ) * e1e2t(ji  ,jj+1)
            zmassU = 0.5_wp*( zr + zr1 )*r1_e1e2u(ji,jj)*umask(ji,jj,1) ! -> mass @ U [kg/m^2]
            zmassV = 0.5_wp*( zr + zr2 )*r1_e1e2v(ji,jj)*vmask(ji,jj,1) ! -> mass @ V [kg/m^2]

            IF( ln_bri_rk3 ) THEN
               ! 1/m !
               xxmU(ji,jj) = 1._wp / MAX( zmassU, epsi20 )
               xxmV(ji,jj) = 1._wp / MAX( zmassV, epsi20 )
            ELSE
               ! m/dt !
               xxmU(ji,jj)   = zmassU * r1_dtbri  ! [kg/m^2/s]
               xxmV(ji,jj)   = zmassV * r1_dtbri  !      "
            ENDIF

            ! Surface pressure gradient (- m*g*GRAD(ssh)) at U-V points
            !  ==> not multiplied by the mass here, do it if required in the momentum equation (Euler, not RK3)
            !      because `m*d[u]/dt = ... -m.g.grad[SSH]` => canceled depending of the form...
            xgrdH(ji,jj,1) = - grav * ( xtmp1(ji+1,jj) - xtmp1(ji,jj)   ) * r1_e1u(ji,jj)
            xgrdH(ji,jj,2) = - grav * ( xtmp1(ji,jj+1) - xtmp1(ji,jj)   ) * r1_e2v(ji,jj)
            xgrdH(ji,jj,3) = - grav * ( xtmp2(ji,jj)   - xtmp2(ji-1,jj) ) * r1_e1v(ji,jj)  ! `xtmp2` is `xtmp1` interpolated  @F !
            xgrdH(ji,jj,4) = - grav * ( xtmp2(ji,jj)   - xtmp2(ji,jj-1) ) * r1_e2u(ji,jj)  ! `xtmp2` is `xtmp1` interpolated  @F !

            ! masks
            kmsk00x(ji,jj) = MERGE( 0 , 1 , zmassU<=0._wp .OR. umask(ji,jj,1)<0.1_wp )
            kmsk00y(ji,jj) = MERGE( 0 , 1 , zmassV<=0._wp .OR. vmask(ji,jj,1)<0.1_wp )

            ! switches
            kmsk01x(ji,jj) = MERGE( 0 , 1 , zmassU <= rMmin_vel .AND. au_i(ji,jj) <= rAmin_vel )
            kmsk01y(ji,jj) = MERGE( 0 , 1 , zmassV <= rMmin_vel .AND. av_i(ji,jj) <= rAmin_vel )

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
               ! * zr => temp. difference between bottom and surface:
               zr = MERGE( t_bo(ji,jj) - tm_su(ji,jj) , -1.8_wp + 25._wp ,  ln_icethd ) ! use a "fake" one when no thermo...
               zr1 = 1._wp / MAX(at_i(ji,jj),epsi06)
               zr3 = rcnd_i*vt_s(ji,jj)*zr1 / MAX( rcnd_s*vt_i(ji,jj)*zr1, epsi06 ) * zmsk   ! => `C` of the
               xtmp1(ji,jj) = zr / (1._wp + zr3 ) * zmsk
               xtmp2(ji,jj) = rDt_ice * MIN( xtmp1(ji,jj) / rn_kth , 1._wp/rDt_ice )  ! dt * 1/T_relax => `1-d` increment @ T
            END DO
         END DO
         !$acc end parallel loop
         CALL do_rmpT2F( xtmp1, xtmp3 )
         CALL do_rmpT2F( xtmp2, xtmp4 )
         !$acc parallel loop collapse(2) present(xtmp3,xtmp4)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               ! Apply healing on damage:
               zr1 = MERGE( xtmp2(ji,jj) , 0._wp ,  xtmp1(ji,jj) > 0._wp )
               zr3 = MERGE( xtmp4(ji,jj) , 0._wp ,  xtmp3(ji,jj) > 0._wp )
               dmdt(ji,jj) = dmdt(ji,jj) + zr1
               dmdf(ji,jj) = dmdf(ji,jj) + zr3
            END DO
         END DO
         !
         CALL cap_1md( at_i, af_i, dmdt, dmdf ) ! Capping for both  post-healing ("post-advection" is done in `icedyn_adv`):
         !
         !! => so that we avoid a `lbc_lnk` prior of cross-nudging into `update_sigma_d()`:
         CALL lbc_lnk( crtnm, dmdt,'T',1._wp,  dmdf,'F',1._wp  )
         !
      ENDIF !IF( l_apply_d_healing )


      zravrg = 1._wp/REAL(nbrttl) ! going to average (set to 0 before accumulating during the `nbrttl` sub time steps)


      !! We do not want to update the following arrays`nbrttl` times in the insane upcoming
      !! `jter=1, nbrttl` loop below; because they remain unchanged since they only depend on
      !! stuff that does not change during the `jter=1, nbrttl` loop (such as `A` for example).
      !!
      !! Mind: stepping 1 point inside the halo is mandatory here for:
      !!       `xxpCt` & `xPmax_t` as wells as their F-point counterparts
      !!        so that we can avoid calling `nbrttl` times a `lbc_lnk` prior to calling
      !!        the cross-nudging into `update_sigma_d()`...

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ! => `nn_hls` required if cross-nudging is used...
            !
            zht = xht(ji,jj)
            zhf = xhf(ji,jj)
            !
            zzt = EXP( rn_C0*(1._wp - at_i(ji,jj)) )
            zzf = EXP( rn_C0*(1._wp - af_i(ji,jj)) )
            !
            xxpCt(ji,jj) = zzt
            xxpCf(ji,jj) = zzf
            !
            xPmax_t(ji,jj)  = -rn_P0 * zzt * zht**2.5_wp ! `2.5` because vert. integrated sigmas => `h^3/2 * h`
            xPmax_f(ji,jj)  = -rn_P0 * zzf * zhf**2.5_wp !   "               "               "             "
            !
            xScHt(ji,jj)  = ( rn_l_ref / REAL(res_grd_loc_t(ji,jj),wp) )**rn_pow_scl_res * zht  ! to scale cohesion / Δx in Mohr-Coulomn test (! multiply with `h` "  " ")
            xScHf(ji,jj)  = ( rn_l_ref / REAL(res_grd_loc_f(ji,jj),wp) )**rn_pow_scl_res * zhf  !       "        "       "        "       "        "       "        "
            !
            ! Interpolate components of air-ice wind stress vector at relevant points:
            zr1 = (2._wp-umask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji+1,jj))
            zr2 = (2._wp-vmask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji,jj+1))
            xtmp1(ji,jj) = 0.5_wp*(taux_ai_t(ji,jj) + taux_ai_t(ji+1,jj)) * zr1 ! x-component of air-ice wind stress at T-point to U-point [kg m^-1 s^-2]
            xtmp2(ji,jj) = 0.5_wp*(tauy_ai_t(ji,jj) + tauy_ai_t(ji,jj+1)) * zr2 ! y-component of air-ice wind stress at T-point to V-point [kg m^-1 s^-2]
            xtmp3(ji,jj) = 0.5_wp*(taux_ai_t(ji,jj) + taux_ai_t(ji,jj+1)) * zr2 ! x-component of air-ice wind stress at T-point to V-point [kg m^-1 s^-2]
            xtmp4(ji,jj) = 0.5_wp*(tauy_ai_t(ji,jj) + tauy_ai_t(ji+1,jj)) * zr1 ! y-component of air-ice wind stress at T-point to U-point [kg m^-1 s^-2]
            !
         ENDDO
      ENDDO
      !$acc end parallel loop

      !LOLOdebug:
      !CALL lbc_lnk( crtnm, xPmax_t(:,:),'T',1._wp, xScHt(:,:),'T',1._wp, xxpCt(:,:),'T',1._wp, &
      !   &                 xPmax_f(:,:),'F',1._wp, xScHf(:,:),'F',1._wp, xxpCf(:,:),'F',1._wp )
      !LOLOdebug.

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            u_ice(ji,jj) = 0._wp
            v_ice(ji,jj) = 0._wp
            uVice(ji,jj) = 0._wp
            vUice(ji,jj) = 0._wp
         ENDDO
      ENDDO
      !$acc end parallel loop

      !$acc loop seq                         ! ==================== !
      DO jter = 1 , nbrttl                   !    loop over jter    !
         !                                   ! ==================== !

         ! ---  Updates the components of the vertically-integrated internal stress tensor and the damage in both T- & F-centric worlds ---
         !           => based on previously computed ice velocities...
         CALL update_sigma_d( kt, jter, rdtbri, V_ts, at_i, af_i, xxpCt, xxpCf, xScHt, xScHf, xPmax_t, xPmax_f, &
            &                                   xht, xhf, SIGMAt, SIGMAf, dmdt, dmdf, SI1t, SI2t, SI1f, SI2f )
         !           => `sigmas` & `d` ok on whole `Nis0-1:Nie0+1,Njs0-1:Nje0+1` ! (provided `V_ts` was fully lbclinked!)

         ! --- Computation of ice velocity --- ! (`Nis0:Nie0,Njs0:Nje0`)
         IF( ln_bri_rk3 ) THEN
            ! ==> Update velocities using IMPLICIT RK3 scheme "Lobatto IIIA" (implicit in terms of ice velocity for the bottom ice-water drag)
            CALL update_uv_rk3( jter, rdtbri, au_i, av_i, xxmU, xxmV, SIGMAt, SIGMAf, xgrdH, V_oce, xtmp1, xtmp2, xtmp3, xtmp4, &
               &                              kmsk01x, kmsk01y, kmsk00x, kmsk00y,  V_ts )
         ELSE
            ! ==> Update velocities using 1st-order Euler scheme in time (implicit in terms of ice velocity for the bottom ice-water drag)
            CALL update_uv_eul( jter, rdtbri, au_i, av_i, xxmU, xxmV, SIGMAt, SIGMAf, xgrdH, V_oce, xtmp1, xtmp2, xtmp3, xtmp4, &
               &                              kmsk01x, kmsk01y, kmsk00x, kmsk00y,  V_ts )
         ENDIF
         !
         !
#if defined _OPENACC || defined _OPENMP
         CALL lbc_lnk_gpu( crtnm, V_ts )
#else
         CALL lbc_lnk(     crtnm, V_ts(:,:,1),'U',-1._wp, V_ts(:,:,2),'V',-1._wp,  V_ts(:,:,3),'V',-1._wp, V_ts(:,:,4),'U',-1._wp )
#endif

         IF( ln_bdy ) THEN
            CALL bdy_ice_dyn( 'U', V_ts(:,:,1) )
            CALL bdy_ice_dyn( 'V', V_ts(:,:,2) )
            CALL bdy_ice_dyn( 'V', V_ts(:,:,3), l_FcVel=.TRUE. )
            CALL bdy_ice_dyn( 'U', V_ts(:,:,4), l_FcVel=.TRUE. )
            CALL bdy_ice_dmg( kt, jter, V_ts )
         ENDIF

         !! Average the velocity (to be used for advection at the big time step):
         !!  ==> MIND: this is the ice velocity at time level `kt+1/2` (not `kt+1`)
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
      END DO !DO jter = 1 , nbrttl                       !  end loop over jter  !
      !                                                ! ==================== !


      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!
      CALL strain_rate_dsd( 'T', u_ice, v_ice, uVice, vUice, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e1t2, e2t2, xmskt, &
         &                       pdivu_i, pshear_i, pdelta_i )
#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( crtnm, pdivu_i, pshear_i, pdelta_i )
#else
      CALL lbc_lnk(     crtnm, pdivu_i,'T',1._wp, pshear_i,'T',1._wp, pdelta_i,'T',1._wp )
#endif

      !$acc end data
      IF( ln_timing )   CALL timing_stop(crtnm)
      !
   END SUBROUTINE ice_dyn_rhg_bbm


   SUBROUTINE ice_dyn_rhg_bbm_init( )
      !!-------------------------------------------------------------------
      !! Called into `ice_dyn_rhg_init()@icedyn_rhg.F90`
      !!-------------------------------------------------------------------
      INTEGER  ::   ierr
      REAL(wp) ::   ztmp, zdx_m, zce, zdts
      REAL(wp), DIMENSION(jpi,jpj) :: zt1, zt2
      !!-------------------------------------------------------------------
      IF( lwp ) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) '    *** Initialization of BBM rheology ('//crtnm//'_init) ***'
      ENDIF

      !! Stiffness matrix
      !! ****************
      rk0  = 1._wp / ( 1._wp - rnup*rnup)
      rk11 = rk0
      rk12 = rk0 * rnup
      rk22 = rk0
      rk33 = rk0 * (1._wp - rnup)

      rsqrt_E0      = SQRT( rn_E0 )

      ! Find the smallest `Δx` of the whole WET domain:
      zt1(:,:) = REAL( res_grd_loc_t(:,:) , 4 ) ! SQRT(Δx*Δy)
      zt1(:,:) = MERGE( zt1(:,:) , 1.E12_wp , (xmskt(:,:) > 0.9_wp) )  ! => stupidly big value over continents...
      zdx_m = MINVAL( zt1 )                        ! min of Δx local
      CALL mpp_min( crtnm//'_init', zdx_m) ! min of dx over the whole domain

      ! time step adjusted automatically
      ! ********************************
      ! we look at the propagation speed of elastic waves (zce), 2) the time step to solve them (zdts)
      !        and the the number of iterations needed (`nbrttl`) to go from t to t+dt
      zce = rsqrt_E0 / rsqrt_nu_rhoi            ! propagation speed of shearing elastic waves based on the mean `Δx`
      zdts = 0.5_wp*(zdx_m/zce)                 ! largest theoretical small time-step to consider to resolve these waves
      zdts = rn_cfl_coeff * zdts                ! with user input
      nbrttl = INT( CEILING( rDt_ice / zdts ) ) ! number of iterations (local to init, cf below for nbrttl)
      nbrttl = nbrttl + MOD(nbrttl, 2)          ! + make it an odd number

      rdtbri   = rDt_ice/REAL(nbrttl,wp)
      r1_dtbri = 1._wp / rdtbri

      CALL cross_nudging_init()

      IF( lwp ) THEN
         WRITE(numout,*) '  * Big time step (advection & thermo)  => rDt_ice  =', REAL(rDt_ice,4), ' [s]'
         WRITE(numout,*) '  * Min `dx` of wet computational domain = ', REAL(zdx_m/1000._wp,4), ' [km]'
         WRITE(numout,*) '     ==> propagation speed of shearing elastic waves => c_e =', REAL(zce,4), '[m/s]'
         WRITE(numout,*) '     ==> time-step requirement to resolve these waves: dt = 0.5*(dx/c_e) =', REAL(0.5_wp*(zdx_m/zce),4), ' [s]'
         WRITE(numout,*) '     ==> adaptation of `dt` based on pick of `rn_cfl_coeff`  :        dt =', REAL(zdts,4), ' [s]'
         WRITE(numout,*) '         => implies a `nbrttl` =', INT(nbrttl,2)
         WRITE(numout,*) '         => implies a sub time step (rheology) => rdtbri =', REAL(rdtbri,4),  ' [s]'
         IF(ln_MCx_test) WRITE(numout,*) '  * Will perform only 1 Mohr-Coulomb test, at mid-point between T & F points!'
         IF(l_CN) THEN
            WRITE(numout,*) '  * About cross-nudging:'
            WRITE(numout,*) '      - CN parameter (gamma) => rn_crndg =', REAL(rn_crndg,4),' [-]'
         ENDIF
         WRITE(numout,*) '  * (scaled) Compression threshod => N_lim =',REAL(rn_N_ref*SQRT( rn_l_ref/zdx_m ),4),' [Pa]'
         WRITE(numout,*) ''
      ENDIF

      IF(rdtbri>zdts) CALL ctl_stop( crtnm//'_init: `nbrttl` is probably to small' )

      IF( lwp ) THEN
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) ''
      ENDIF

      ALLOCATE( kmsk01x(jpi,jpj), kmsk01y(jpi,jpj),  kmsk00x(jpi,jpj), kmsk00y(jpi,jpj),  &
         &      xxmU(jpi,jpj), xxmV(jpi,jpj),      xxpCt(jpi,jpj),   xxpCf(jpi,jpj),  &
         &      xgrdH(jpi,jpj,4), xScHt(jpi,jpj), xScHf(jpi,jpj),                    &
         &      xPmax_t(jpi,jpj), xPmax_f(jpi,jpj),                                       &
         &      STAT = ierr )
      !
      kmsk01x(:,:) = 0 ;  kmsk01y(:,:) = 0 ; kmsk00x(:,:) = 0 ; kmsk00y(:,:) = 0
      xgrdH(:,:,:) = 0._wp

      CALL mpp_sum( crtnm//'_init', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', crtnm//'_init : unable to allocate work arrays')

#if defined _OPENACC || defined _OPENMP
      !$acc update device ( nbrttl, rk0, rk11, rk22, rk12, rk33, rsqrt_E0, rdtbri, r1_dtbri )
      PRINT *, ' * info GPU: '//crtnm//'_init() => adding work arrays to memory!'
      !$acc enter data copyin( kmsk01x, kmsk01y, kmsk00x, kmsk00y, xxmU, xxmV )
      PRINT *, '    ==> kmsk01x, kmsk01y, kmsk00x, kmsk00y, xxmU, xxmV'
      !$acc enter data copyin( xxpCt, xxpCf, xgrdH, xPmax_t, xPmax_f, xScHt, xScHf )
      PRINT *, '    ==> xxpCt, xxpCf, xgrdH, xPmax_t, xPmax_f, xScHt, xScHf'
#endif

   END SUBROUTINE ice_dyn_rhg_bbm_init


   SUBROUTINE update_sigma_d( kt, kts, pdt, pV4, pAt, pAf, pxpCt, pxpCf, pScHt, pScHf, pPmax_t, pPmax_f, &
      &                                     pht, phf, psgmt, psgmf, p1mdt, p1mdf, pSI1t, pSI2t, pSI1f, pSI2f )
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
      INTEGER,                        INTENT(in)    :: kt, kts          ! # of current big and small/sub-time step
      REAL(wp),                       INTENT(in)    :: pdt              ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in)    :: pV4              ! the 4 components of sea-ice velocity
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pAt, pAf         ! Ice concentration @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf     ! EXP( rn_C0*(1 - A) ) @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pScHt, pScHf ! `Scaling factor * h` (w.r.t. local `Δx`) for the cohesion and `Nlim` @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pPmax_t, pPmax_f ! `Pmax * h`  @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pht, phf         ! Ice thickness @T & @F
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt            ! T-centric vertically-integrated stress tensor [N/m^2*m]
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf            ! F-centric vertically-integrated stress tensor [N/m^2*m]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf     ! `1 - ice damage` @T & @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: pSI1t, pSI2t, pSI1f, pSI2f ! 1st & 2nd invariants of vertically-integrated stress tensor [N/m^2*m]
      !!----------------------------------------------------------------------
      REAL(wp) :: zfc, zml, zdum
      REAL(wp) :: zs11, zs22, zs12
      REAL(wp) :: zh, zE, zeta, zL, zang, zc0, zmul
      REAL(wp) :: zE1, zE2
      REAL(wp) :: ze11t, ze22t, ze12t, ze11f, ze22f, ze12f
      REAL(wp) :: zxpC, z1md, zsigI, zsigII, zPmax, zPtld
      REAL(wp) :: z1_sigI, zecc, zCohe
      REAL(wp) :: zrr
      !!----------------------------------------------------------------------
      REAL(wp)               :: zxpCt, zScHt, zPmax_t, zht, z1mdt, zmskt
      REAL(wp)               :: zxpCf, zScHf, zPmax_f, zhf, z1mdf, zmskf
      REAL(wp)               :: zSI1t, zSI2t, zSI1f, zSI2f
      REAL(wp), DIMENSION(6) :: zSR ! T- & F- centric strain-rate tensor components at time `t=k` [s^-1] // [1,2,6]@T [4,5,3]@F
      REAL(wp), DIMENSION(3) :: zsgmt, zsgmf
      REAL(wp), DIMENSION(6) :: zK
      REAL(wp), DIMENSION(2) :: zG
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jp
      INTEGER  :: khep ! go `khep` points into the halo
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_sigma_d')

      IF( (.NOT. l_CN).AND.(kt==nit000).AND.(lwp) ) WRITE(numout,*) ' *** MIND: no cross-nudging will be applied between stress tensors!'

      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, psgmt(:,:,1), 'Sigma11t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, psgmt(:,:,2), 'Sigma22t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, psgmt(:,:,3), 'Sigma12f' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, psgmf(:,:,1), 'Sigma11f' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, psgmf(:,:,2), 'Sigma22f' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, psgmf(:,:,3), 'Sigma12t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, pSI1t, 'sigmaI_t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, pSI2t, 'sigmaII_t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, pSI1f, 'sigmaI_f' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, pSI2f, 'sigmaII_f' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'U', 1._wp, pV4(:,:,1), 'u_ice' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'V', 1._wp, pV4(:,:,2), 'v_ice' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'U', 1._wp, pV4(:,:,3), 'uVice' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'V', 1._wp, pV4(:,:,4), 'vUice' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, pht,      'pht' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, pxpCt,    'pxpCt' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, p1mdt,    'p1mdt' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, psgmt,    'psgmt' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'T', 1._wp, pPmax_t,  'pPmax_t' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, phf,      'phf' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, pxpCf,    'pxpCf' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, p1mdf,    'p1mdf' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, psgmf,    'psgmf' )
      !CALL debug_test_if_lbclnked( 'update_sigma_d', 'F', 1._wp, pPmax_f,  'pPmax_f' )

      !**************************************
      IF( ln_bri_rk3 ) THEN
         ! =================================================================
         ! Update of stress tensors using implicit RK3 scheme "Lobatto IIIA"
         ! =================================================================

         !$acc data present( pV4,pxpCt,pxpCf,pScHt,pScHf,pPmax_t,pPmax_f,pht,phf,psgmt,psgmf,p1mdt,p1mdf ) create( zSR,zsgmt,zsgmf,zK,zG )

         !$acc parallel loop collapse(2) private( zSR,zsgmt,zsgmf,zK,zG )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1

               zmskt = xmskt(ji,jj) ; zmskf = xmskf(ji,jj)

               !! *** T-centric grid ***
#              include "icedyn_rhg_bri_strn_t.h90"
               !        => uses `ji-1` & `jj-1`
               ! ==> ze11t, ze22t, ze12t

               !! *** F-centric grid ***

#              include "icedyn_rhg_bri_strn_f.h90"
               !        => uses `ji+1` & `jj+1`
               ! ==> ze11f, ze22f, ze12f

               zSR(1) = ze11t ; zSR(4) = ze11f
               zSR(2) = ze22t ; zSR(5) = ze22f
               zSR(6) = ze12t ; zSR(3) = ze12f

               ! Now time for the RK3 stuff..
               zxpCt    =   pxpCt(ji,jj)    ;   zxpCf    =   pxpCf(ji,jj)
               zScHt  = pScHt(ji,jj)    ;   zScHf  = pScHf(ji,jj)
               zPmax_t  = pPmax_t(ji,jj)    ;   zPmax_f  = pPmax_f(ji,jj)
               zht      =     pht(ji,jj)    ;   zhf      =     phf(ji,jj)
               z1mdt    =   p1mdt(ji,jj)    ;   z1mdf    =   p1mdf(ji,jj)
               zmskt    =   xmskt(ji,jj)    ;   zmskf    =   xmskf(ji,jj)
               !$acc loop seq
               DO jp=1, 3
                  ! Because of a `ld.lld: error: undefined symbol: _FortranAAssign` bug!
                  zsgmt(jp) = psgmt(ji,jj,jp)
                  zsgmf(jp) = psgmf(ji,jj,jp)
               END DO
               zSI1t    =   pSI1t(ji,jj)    ;   zSI1f    =   pSI1f(ji,jj)
               zSI2t    =   pSI2t(ji,jj)    ;   zSI2f    =   pSI2f(ji,jj)

               ! IMPLICIT RK3
               CALL RHS_sigmah_impl( pdt, zSR, zxpCt, zxpCf, zScHt, zScHf, zPmax_t, zPmax_f, &
                  &                  zht, zhf, zsgmt, zsgmf, zSI1t, zSI2t, zSI1f, zSI2f, z1mdt, z1mdf, &
                  &                  zmskt, zmskf,  zK, zG )

               ! stresses @ T-points:
               psgmt(ji,jj,1) = newS_rk3_imp( pdt, zsgmt(1), zK(1), zG(1) )
               psgmt(ji,jj,2) = newS_rk3_imp( pdt, zsgmt(2), zK(2), zG(1) )
               psgmf(ji,jj,3) = newS_rk3_imp( pdt, zsgmf(3), zK(6), zG(1) )
               ! stresses @ F-points:
               psgmf(ji,jj,1) = newS_rk3_imp( pdt, zsgmf(1), zK(4), zG(2) )
               psgmf(ji,jj,2) = newS_rk3_imp( pdt, zsgmf(2), zK(5), zG(2) )
               psgmt(ji,jj,3) = newS_rk3_imp( pdt, zsgmt(3), zK(3), zG(2) )

            END DO !ji=Nis0-1, Nie0+1
         END DO !jj=Njs0-1, Nje0+1

         !$acc end data
         !**************************************
      ELSE
         ! =======================================================
         ! Update of stress tensors using implicit Euler 1st order
         ! =======================================================

         !$acc data present( pV4,pxpCt,pxpCf,pScHt,pScHf,pPmax_t,pPmax_f,pht,phf,psgmt,psgmf,p1mdt,p1mdf )

         khep = MERGE( 1 , 0 , ln_MCx_test )

         !! *** T-centric grid ***
         !! ~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-khep, Nje0+1
            DO ji=Nis0-khep, Nie0+1

               zmskt = xmskt(ji,jj)
               zh    =   pht(ji,jj)

               ! --- Strain rate tensors ---
               !     *******************
#              include "icedyn_rhg_bri_strn_t.h90"
               !        => uses `ji-1` & `jj-1`
               ! ==> ze11t, ze22t, ze12t

               ! --- E, Lambda & multiplicator (uses `sigmas`!!!) ---
               !     ********************************************
#              include "icedyn_rhg_bbm_elm_t.h90"
               ! ==> zEt, zLt, zmult

               ! --- Predictor estimate of stress tensor at k+1 ---
               !     ******************************************
               zml = zmul * zmskt
               zfc =  zh * zE * pdt
               psgmt(ji,jj,1) = zml * ( zfc * ( rk11*ze11t + rk12*ze22t) + psgmt(ji,jj,1) )
               psgmt(ji,jj,2) = zml * ( zfc * ( rk12*ze11t + rk22*ze22t) + psgmt(ji,jj,2) )
               psgmf(ji,jj,3) = zml * ( zfc *        rk33 * ze12t        + psgmf(ji,jj,3) )

            END DO !DO ji=Nis0, Nie0
         END DO !DO jj=Njs0, Nje0
         !$acc end parallel loop

         !! *** F-centric grid ***
         !! ~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+khep
            DO ji=Nis0-1, Nie0+khep

               zmskf = xmskf(ji,jj)
               zh    =   phf(ji,jj)

               ! --- Strain rate tensors ---
               !     *******************
#              include "icedyn_rhg_bri_strn_f.h90"
               !        => uses `ji+1` & `jj+1`
               ! ==> ze11f, ze22f, ze12f

               ! --- E, Lambda & multiplicator (uses `sigmas`!!!) ---
               !     ********************************************
#              include "icedyn_rhg_bbm_elm_f.h90"
               ! ==> zEf, zLf, zmulf

               ! --- Predictor estimate of stress tensor at k+1 ---
               !     ******************************************
               zml = zmul * zmskf
               zfc = zh * zE * pdt
               psgmf(ji,jj,1) = zml * ( zfc * ( rk11*ze11f + rk12*ze22f) + psgmf(ji,jj,1) )
               psgmf(ji,jj,2) = zml * ( zfc * ( rk12*ze11f + rk22*ze22f) + psgmf(ji,jj,2) )
               psgmt(ji,jj,3) = zml * ( zfc *        rk33 * ze12f        + psgmt(ji,jj,3) )

            END DO !DO ji=Nis0, Nie0
         END DO !DO jj=Njs0, Nje0
         !$acc end parallel loop

         !$acc end data
         !**************************************
      ENDIF !IF( ln_bri_rk3 )
      !**************************************

      IF( l_CN ) THEN
         !LOLOdebug:
         !CALL lbc_lnk( crtnm, psgmt(:,:,1),'T',1._wp, psgmt(:,:,2),'T',1._wp, psgmt(:,:,3),'F',1._wp, &
         !   &                 psgmf(:,:,1),'F',1._wp, psgmf(:,:,2),'F',1._wp, psgmf(:,:,3),'T',1._wp  )
         !LOLOdebug.
         !              ******** CROSS NUDGING *********
         CALL apply_CN( kts, psgmt, psgmf ) !  `Nis0-1:Nie0+1,Njs0-1:Nje0+1`
         !      !!          !        => @T uses `ji+1` & `jj+1`  result => ok on whole `Nis0:Nie0,Njs0:Nje0`
         !      !!          !        => @F uses `ji-1` & `jj-1`  result => ok on whole      "          "
      ENDIF

      ! ******** MOHR-COULOMB FAILURE TEST & UPDATE of damage and stress tensors  *********
      IF( ln_MCx_test ) THEN
         IF( kt==nit000 .AND. kts==1 ) CALL MC_ud_d_s_mwp_init()
         CALL MC_ud_d_s_mwp( pdt, pxpCt, pxpCf, pScHt, pScHf, p1mdt, p1mdf, psgmt, psgmf )
      ELSE
         CALL MC_ud_d_s(     pdt, pxpCt, pxpCf, pScHt, pScHf, p1mdt, p1mdf, psgmt, psgmf )
      ENDIF
      ! ***********************************************************************************

      CALL clean_small_a_all( pAt, pAf,  p1mdt, p1mdf, psgmt, psgmf )

      !! --- lbc-linking of updated stress tensors and `1-damage` ---
#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( crtnm, psgmt, psgmf, p1mdt, p1mdf )
#else
      CALL lbc_lnk(     crtnm, psgmt(:,:,1),'T',1._wp, psgmt(:,:,2),'T',1._wp, psgmt(:,:,3),'F',1._wp, p1mdt,'T',1._wp, &
         &                     psgmf(:,:,1),'F',1._wp, psgmf(:,:,2),'F',1._wp, psgmf(:,:,3),'T',1._wp, p1mdf,'F',1._wp  )
#endif

      !! Now that stress tensors are `lbc_lnk`ed, we can update invariants:
      !$acc data present( pSI1t, pSI2t, pSI1f, pSI2f )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            zs11 = psgmt(ji,jj,1) ; zs22 = psgmt(ji,jj,2) ; zs12 = psgmf(ji,jj,3)
            zrr          = 0.5_wp * ( zs11 - zs22 )
            pSI1t(ji,jj) = 0.5_wp * ( zs11 + zs22 )
            pSI2t(ji,jj) = SQRT( zrr*zrr + zs12*zs12 )
            !
            zs11 = psgmf(ji,jj,1) ; zs22 = psgmf(ji,jj,2) ; zs12 = psgmt(ji,jj,3)
            zrr          = 0.5_wp * ( zs11 - zs22 )
            pSI1f(ji,jj) = 0.5_wp * ( zs11 + zs22 )
            pSI2f(ji,jj) = SQRT( zrr*zrr + zs12*zs12 )
         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !$acc end parallel loop
      !$acc end data

      IF( ln_timing )   CALL timing_stop('update_sigma_d')
      !
   END SUBROUTINE update_sigma_d
















   !! RK3 implicit for Sigma update !!



   FUNCTION newS_rk3_imp( pdt, pSh, pK, pG )
      !!----------------------------------------------------------------------------------------------
      !$acc routine seq
      !!----------------------------------------------------------------------------------------------
      !!   Implicit RK3 operator for internal stress tensor update
      !!     => using the coefficients of "Lobatto IIIA"
      !!----------------------------------------------------------------------------------------------
      REAL(wp)             :: newS_rk3_imp
      !!----------------------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pdt   ! small time-step [s]
      REAL(wp), INTENT(in) :: pSh   ! a `σ*h` at time k
      REAL(wp), INTENT(in) :: pK    ! `E*h*|K|:ε` i.e. all the content of the RHS exluding the `-σ*h/λ*(1+P~)` term
      REAL(wp), INTENT(in) :: pG    ! `1/λ*(1+P~)` term
      !!----------------------------------------------------------------------------------------------
      REAL(wp) :: zdum
      !!----------------------------------------------------------------------------------------------
      zdum = pdt*pG
      !
      newS_rk3_imp =  pSh  +  pdt * ( 12._wp*(pK - pG*pSh) ) / ( 12._wp + 6._wp*zdum + zdum*zdum )
      !
   END FUNCTION newS_rk3_imp


   SUBROUTINE RHS_sigmah_impl( pdt, pSR, pxpCt, pxpCf, pSclHt, pSclHf, pPmaxt, pPmaxf, pht, phf, &
      &                             psgmt, psgmf, pSI1t, pSI2t, pSI1f, pSI2f, p1mdt, p1mdf, pmskt, pmskf,  pK, pG )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE RHS_SIGMA  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!               implicitness is in the `σ*h` of term `-σ*h/λ*(1+P~)`, eq.(33) of Olason et al.
      !!
      !! ** Note    : Called at the sub-time-stepping level!
      !!
      !! ** Author : L. Brodeau, 2026
      !!----------------------------------------------------------------------
      !$acc routine seq
      !!----------------------------------------------------------------------
      REAL(wp),               INTENT(in)  :: pdt              ! (small) time-step [s]
      REAL(wp), DIMENSION(6), INTENT(in)  :: pSR              ! T- & F- centric strain-rate tensor components at time `t=k` [s^-1]
      REAL(wp),               INTENT(in)  :: pxpCt, pxpCf     ! EXP( rn_C0*(1 - A) ) @T & @F
      REAL(wp),               INTENT(in)  :: pSclHt, pSclHf ! `Scaling factor * h` (w.r.t. local `Δx`) for the cohesion and `Nlim` @T & @F
      REAL(wp),               INTENT(in)  :: pPmaxt, pPmaxf ! `Pmax * h` @T & @F (<0)
      REAL(wp),               INTENT(in)  :: pht, phf         ! Ice thickness @T & @F
      REAL(wp), DIMENSION(3), INTENT(in)  :: psgmt, psgmf     ! T- & F-centric vert.-int. stress tensor components at time `t=k` [N/m^2*m]
      REAL(wp),               INTENT(in)  :: pSI1t, pSI2t, pSI1f, pSI2f  ! 1st & 2nd invariant of vert.-int. stress tensor comp. [N/m^2*m]
      REAL(wp),               INTENT(in)  :: p1mdt, p1mdf     ! `1 - ice damage` @T & @F
      REAL(wp),               INTENT(in)  :: pmskt, pmskf     ! land sea mask @T & @F
      REAL(wp), DIMENSION(6), INTENT(out) :: pK               ! terms `E*h*|K|:ε`  => 1 for each stress tensor component
      REAL(wp), DIMENSION(2), INTENT(out) :: pG               ! terms `1/λ*(1+P~)` => 1 for each grid  !#LOLOfixme: shoul be of size 2 !!!
      !!----------------------------------------------------------------------
      REAL(wp) :: zfc, zml, zmsk
      REAL(wp) :: zE, zeta, zL, zc0
      REAL(wp) :: ze11, ze22, ze12
      REAL(wp) :: zxpC, z1md, zsigI, zsigII, zPmax, zSclH, zPtld
      REAL(wp) :: z1_sigI, zecc, zCohe
      !!----------------------------------------------------------------------
      !! ~~~~~~~~~~~~~~~~~~~~~~
      !! *** T-centric grid ***
      !! ~~~~~~~~~~~~~~~~~~~~~~
      zmsk = pmskt
      ! --- Internal stress tensor components @T ---
      zsigI  = pSI1t
      zsigII = pSI2t
      ! --- `E`, `Lambda` & `P~` (uses `sigmas`!!!) ---
      zxpC  = pxpCt
      z1md  = p1mdt
      zPmax = pPmaxt ! `P_max * h` (<0)
      zSclH = pSclHt
#     include "icedyn_rhg_bbm_elmexpl.h90"
      ! --- Strain rate tensor components @T ---
      ze11 = pSR(1)
      ze22 = pSR(2)
      ze12 = pSR(6)
      ! --- RHS of Eq.(20) of Olason et al.  ---
      zfc   =  pht * zE
      pK(1) = zfc * ( rk11*ze11 + rk12*ze22) * zmsk
      pK(2) = zfc * ( rk12*ze11 + rk22*ze22) * zmsk
      pK(6) = zfc *         rk33 * ze12      * zmsk
      pG(1) = (1._wp + zPtld) / zL * zmsk

      !! ~~~~~~~~~~~~~~~~~~~~~~
      !! *** F-centric grid ***
      !! ~~~~~~~~~~~~~~~~~~~~~~
      zmsk = pmskf
      ! --- Internal stress tensor components @T ---
      zsigI  = pSI1f
      zsigII = pSI2f
      ! --- `E`, `Lambda` & `P~` (uses `sigmas`!!!) ---
      zxpC  = pxpCf
      z1md  = p1mdf
      zPmax = pPmaxf ! `P_max * h` (<0)
      zSclH = pSclHf
#     include "icedyn_rhg_bbm_elmexpl.h90"
      ! --- Strain rate tensor components @F ---
      ze11 = pSR(4)
      ze22 = pSR(5)
      ze12 = pSR(3)
      ! --- RHS of Eq.(20) of Olason et al.  ---
      zfc   = phf * zE
      pK(4) = zfc * ( rk11*ze11 + rk12*ze22) * zmsk
      pK(5) = zfc * ( rk12*ze11 + rk22*ze22) * zmsk
      pK(3) = zfc *         rk33 * ze12      * zmsk
      pG(2) = (1._wp + zPtld) / zL * zmsk
      !
   END SUBROUTINE RHS_sigmah_impl


   !!==============================================================================
END MODULE icedyn_rhg_bbm
