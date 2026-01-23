MODULE icedyn_rhg_evp
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_evp  ***
   !!   Sea-Ice dynamics : rheology Elasto-Viscous-Plastic
   !!======================================================================
   !! History :   -   !  2007-03  (M.A. Morales Maqueda, S. Bouillon) Original code
   !!            3.0  !  2008-03  (M. Vancoppenolle) adaptation to new model
   !!             -   !  2008-11  (M. Vancoppenolle, S. Bouillon, Y. Aksenov) add surface tilt in ice rheolohy
   !!            3.3  !  2009-05  (G.Garric)    addition of the evp case
   !!            3.4  !  2011-01  (A. Porter)   dynamical allocation
   !!            3.5  !  2012-08  (R. Benshila) AGRIF
   !!            3.6  !  2016-06  (C. Rousset)  Rewriting + landfast ice + mEVP (Bouillon 2013)
   !!            3.7  !  2017     (C. Rousset)  add aEVP (Kimmritz 2016-2017)
   !!            4.0  !  2018     (many people) SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_evp : computes ice velocities from EVP rheology
   !!   rhg_evp_rst     : read/write EVP fields in ice restart
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE sbc_oce , ONLY : nn_fsbc
   USE oss_nnq , ONLY : ssh_m
   USE sbc_ice , ONLY : utau_ice, vtau_ice, snwice_mass_b
   USE par_ice
   USE ice            ! sea-ice: ice variables
   USE icevar         ! ice_var_sshdyn
   USE icedyn_rdgrft, ONLY : ice_strength
   USE bdy , ONLY : ln_bdy
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

   USE remap_classic,  ONLY: rmpT2F, rmpU2V, rmpV2U, do_rmpT2F

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_evp_init
   PUBLIC   ice_dyn_rhg_evp   ! called by icedyn_rhg.F90
   PUBLIC   rhg_evp_rst       ! called by icedyn_rhg.F90

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zdelta, zp_delt                 ! delta and P/delta at T points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zbeta                           ! beta coef from Kimmritz 2017
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zdt_m                           ! (dt / ice-snow_mass) on T points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zmU_dt, zmV_dt                  ! (ice-snow_mass / dt) on U/V points
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zht, zhf
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zds                             ! shear
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zten_i, zshear                  ! tension, shear
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zs1, zs2, zs12                  ! stress tensor components
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zsshdyn                         ! array used for the calculation of ice surface slope:
   !                                                                           !    ocean surface (ssh_m) if ice is not embedded
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zfU  , zfV                      ! internal stresses
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   zspgU, zspgV                    ! surface pressure gradient at U/V points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   ztaux_bi, ztauy_bi              ! ice-OceanBottom stress at U-V points (landfast)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)  ::   ztaux_base, ztauy_base          ! ice-bottom stress at U-V points (landfast)
   !
   REAL(wp),   ALLOCATABLE, SAVE, DIMENSION(:,:)  :: zmsk
   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: kmsk01x, kmsk01y                ! dummy arrays
   INTEGER(1), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: kmsk00x, kmsk00y                ! mask for ice presence

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icedyn_rhg_evp.F90 15550 2021-11-28 20:02:31Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg_evp( kt, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_evp  ***
      !!                             EVP-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  a non-linear elasto-viscous-plastic (EVP) law including shear
      !!  strength and a bulk rheology (Hunke and Dukowicz, 2002).
      !!
      !!  The points in the C-grid look like this, dear reader
      !!
      !!                              (ji,jj)
      !!                                 |
      !!                                 |
      !!                      (ji-1,jj)  |  (ji,jj)
      !!                             ---------
      !!                            |         |
      !!                            | (ji,jj) |------(ji,jj)
      !!                            |         |
      !!                             ---------
      !!                     (ji-1,jj-1)     (ji,jj-1)
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
      !! ** Notes   : This aEVP! Based on the nice work of Kimmritz et al. (2016 & 2017)
      !!              (i.e. changing alpha and beta parameters).
      !!              This is an upgraded version of mEVP from Bouillon et al. 2013
      !!              (i.e. more stable and better convergence)
      !!
      !! References : Hunke and Dukowicz, JPO97
      !!              Bouillon et al., Ocean Modelling 2009
      !!              Bouillon et al., Ocean Modelling 2013
      !!              Kimmritz et al., Ocean Modelling 2016 & 2017
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in   ) ::   kt                                    ! time step
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pstress1_i, pstress2_i, pstress12_i   !
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pshear_i  , pdivu_i   , pdelta_i      !
      !!
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   jter         ! local integers
      !
      REAL(wp) ::   zswitch, zmask
      REAL(wp) ::   zrhoco                                              ! rho0 * rn_Cd_io
      REAL(wp) ::   ztau_ai                                             ! ice-atm. stress at U-V points
      REAL(wp) ::   zdtevp, z1_dtevp                                    ! time step for subcycling
      REAL(wp) ::   ecc2, z1_ecc2                                       ! square of yield ellipse eccenticity
      REAL(wp) ::   zalph1, z1_alph1, zalph2, z1_alph2                  ! alpha coef from Bouillon 2009 or Kimmritz 2017
      REAl(wp) ::   zbetau, zbetav
      REAL(wp) ::   zm1, zm2, zm3, zmassU, zmassV, zvU, zvV             ! ice/snow mass and volume
      REAL(wp) ::   zp_delf, zds2, zdt, zdt2, zdiv, zdiv2               ! temporary scalars
      REAL(wp) ::   zTauO, zTauB, zRHS, zvel                            ! temporary scalars
      REAL(wp) ::   ztaux_oi, ztauy_oi                                  ! ice-ocean stress at U-V points
      REAL(wp) ::   zkt                                                 ! isotropic tensile strength for landfast ice
      REAL(wp) ::   zvCr                                                ! critical ice volume above which ice is landfast
      REAL(wp) ::   zcorio, zA, zUi, zVi, zM_dt, zUo, zVo, zt1, zt2
      REAL(wp) ::   zfac_x, zfac_y
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_dyn_rhg_evp')
      !$acc data present( pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i, u_ice, v_ice, uVice, vUice, SIGMAt )

      IF( kt == nit000 .AND. lwp ) WRITE(numout,*) '-- ice_dyn_rhg_evp: EVP sea-ice rheology'

      ! for diagnostics and convergence tests
      !$acc parallel loop collapse(2) present( at_i, zmsk )
      DO jj = Njs0-1, Nje0+1
         DO ji = Nis0-1, Nie0+1
            zmsk  (ji,jj) = MERGE( 1._wp,  0._wp,  at_i(ji,jj) >= epsi10 )   ! 1 if ice    , 0 if no ice
         END DO
      END DO
      !$acc end parallel loop

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
               zds(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
# if ! defined _OPENACC
         CALL lbc_lnk( 'icedyn_rhg_bbm', zhf,'F',1._wp )
# endif
      ELSE
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zht(ji,jj) = hm_i  (ji,jj)
               zhf(ji,jj) = hm_i_f(ji,jj)
               zds(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF !IF( l_use_v_for_h )




      zrhoco = rho0 * rn_Cd_io

      ! ecc2: square of yield ellipse eccenticrity
      ecc2    = rn_ecc * rn_ecc
      z1_ecc2 = 1._wp / ecc2

      ! alpha parameters (Bouillon 2009)
      zdtevp   = rDt_ice
      ! zalpha parameters set later on adaptatively
      z1_dtevp = 1._wp / zdtevp

      ! Initialise stress tensor
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            zs1 (ji,jj) = pstress1_i (ji,jj)
            zs2 (ji,jj) = pstress2_i (ji,jj)
            zs12(ji,jj) = pstress12_i(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      ! Ice strength
      CALL ice_strength

      ! landfast param from Lemieux(2016): add isotropic tensile strength (following Konig Beatty and Holland, 2010)
      zkt = MERGE( rn_lf_tensile,  0._wp,  ln_landfast_L16 )
      !
      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      CALL ice_var_sshdyn( ssh_m, snwice_mass_b, zsshdyn )

      !$acc parallel loop collapse(2) present( zdt_m )
      DO jj = Njs0-1, Nje0+1
         DO ji = Nis0-1, Nie0+1
            zm1          = ( rhos * vt_s(ji,jj) + rhoi * vt_i(ji,jj) )  ! Ice/snow mass at U-V points
            zdt_m(ji,jj) = zdtevp / MAX( zm1, rmass_min )               ! dt/m at T points (for alpha and beta coefficients)
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2) present( zmU_dt, zmV_dt, zspgU, zspgV, kmsk00x, kmsk00y, kmsk01x, kmsk01y )
      DO jj = Njs0, Nje0
         DO ji = Nis0, Nie0

            ! Ice/snow mass at U-V points
            zm1 = ( rhos * vt_s(ji  ,jj  ) + rhoi * vt_i(ji  ,jj  ) )
            zm2 = ( rhos * vt_s(ji+1,jj  ) + rhoi * vt_i(ji+1,jj  ) )
            zm3 = ( rhos * vt_s(ji  ,jj+1) + rhoi * vt_i(ji  ,jj+1) )
            zmassU = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm2 * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
            zmassV = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm3 * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)

            ! m/dt
            zmU_dt(ji,jj)    = zmassU * z1_dtevp
            zmV_dt(ji,jj)    = zmassV * z1_dtevp

            ! Surface pressure gradient (- m*g*GRAD(ssh)) at U-V points
            zspgU(ji,jj) = - zmassU * grav * ( zsshdyn(ji+1,jj) - zsshdyn(ji,jj) ) * r1_e1u(ji,jj)
            zspgV(ji,jj) = - zmassV * grav * ( zsshdyn(ji,jj+1) - zsshdyn(ji,jj) ) * r1_e2v(ji,jj)

            ! masks
            !kmsk00x(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassU ) )
            !kmsk00y(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassV ) )
            kmsk00x(ji,jj) = MERGE( 1,  0,  zmassU > 0._wp )   ! 0 if no ice
            kmsk00y(ji,jj) = MERGE( 1,  0,  zmassV > 0._wp )   ! 0 if no ice

            ! switches
            kmsk01x(ji,jj) = MERGE( 0,  1,  zmassU <= rmass_min .AND. au_i(ji,jj) <= rconc_min )
            kmsk01y(ji,jj) = MERGE( 0,  1,  zmassV <= rmass_min .AND. av_i(ji,jj) <= rconc_min )

         END DO
      END DO
      !$acc end parallel loop

      !
      !                                  !== Landfast ice parameterization ==!
      !
      IF( ln_landfast_L16 ) THEN         !-- Lemieux 2016
# if defined _OPENACC
         CALL ctl_stop('STOP', 'ice_dyn_rhg_evp: add the OpenACC stuff for `ln_landfast_L16` !')
# endif
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               ! ice thickness at U-V points
               zvU = 0.5_wp * ( vt_i(ji,jj) * e1e2t(ji,jj) + vt_i(ji+1,jj) * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
               zvV = 0.5_wp * ( vt_i(ji,jj) * e1e2t(ji,jj) + vt_i(ji,jj+1) * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)
               ! ice-bottom stress at U points
               !LOLO:zvCr = au_i(ji,jj) * rn_lf_depfra * hu(ji,jj,Kmm) * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
               zvCr = au_i(ji,jj) * rn_lf_depfra * 1._wp * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
               ztaux_base(ji,jj) = - rn_lf_bfr * MAX( 0._wp, zvU - zvCr ) * EXP( -rn_crhg * ( 1._wp - au_i(ji,jj) ) )
               ! ice-bottom stress at V points
               !LOLO:zvCr = av_i(ji,jj) * rn_lf_depfra * hv(ji,jj,Kmm) * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
               zvCr = av_i(ji,jj) * rn_lf_depfra * 1._wp * ( 1._wp - icb_mask(ji,jj) ) ! if grounded icebergs are read: ocean depth = 0
               ztauy_base(ji,jj) = - rn_lf_bfr * MAX( 0._wp, zvV - zvCr ) * EXP( -rn_crhg * ( 1._wp - av_i(ji,jj) ) )
               ! ice_bottom stress at T points
               !LOLO:zvCr = at_i(ji,jj) * rn_lf_depfra * ht(ji,jj) * ( 1._wp - icb_mask(ji,jj) )    ! if grounded icebergs are read: ocean depth = 0
               zvCr = at_i(ji,jj) * rn_lf_depfra * 1._wp * ( 1._wp - icb_mask(ji,jj) )    ! if grounded icebergs are read: ocean depth = 0
               tau_icebfr(ji,jj) = - rn_lf_bfr * MAX( 0._wp, vt_i(ji,jj) - zvCr ) * EXP( -rn_crhg * ( 1._wp - at_i(ji,jj) ) )
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_rhg_evp', tau_icebfr(:,:), 'T', 1.0_wp )
         !
      ELSE
         !                   !-- no landfast
         !$acc parallel loop collapse(2) present( ztaux_base, ztauy_base )
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               ztaux_base(ji,jj) = 0._wp
               ztauy_base(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop

      ENDIF !IF( ln_landfast_L16 )

      !------------------------------------------------------------------------------!
      ! 3) Solution of the momentum equation, iterative procedure
      !------------------------------------------------------------------------------!
      !
      !$acc loop seq                                  ! ==================== !
      DO jter = 1 , nn_nevp                           !    loop over jter    !
         !                                            ! ==================== !
         l_full_nf_update = jter == nn_nevp   ! false: disable full North fold update (performances) for iter = 1 to nn_nevp-1

         ! --- divergence, tension & shear (Appendix B of Hunke & Dukowicz, 2002) --- !
         !$acc parallel loop collapse(2) present( zds, r1_e1u, e1f2, r1_e2v, e2f2, r1_e1e2f )
         DO jj = Njs0-1, Nje0
            DO ji = Nis0-1, Nie0
               ! shear at F points
               zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f2(ji,jj)   &
                  &         + ( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f2(ji,jj)   &
                  &         ) * r1_e1e2f(ji,jj) * fmask(ji,jj,1)
            END DO
         END DO
         !$acc end parallel loop

         !$acc parallel loop collapse(2) present( e1e2f, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e2t2, e1t2, zdelta, zp_delt, strength, zmsk )
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0

               ! shear**2 at T points (doc eq. A16)
               zds2 = ( zds(ji,jj  ) * zds(ji,jj  ) * e1e2f(ji,jj  ) + zds(ji-1,jj  ) * zds(ji-1,jj  ) * e1e2f(ji-1,jj  )  &
                  &   + zds(ji,jj-1) * zds(ji,jj-1) * e1e2f(ji,jj-1) + zds(ji-1,jj-1) * zds(ji-1,jj-1) * e1e2f(ji-1,jj-1)  &
                  &   ) * 0.25_wp * r1_e1e2t(ji,jj)

               ! divergence at T points
               zdiv  = ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
                  &    + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
                  &    ) * r1_e1e2t(ji,jj)
               zdiv2 = zdiv * zdiv

               ! tension at T points
               zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t2(ji,jj)   &
                  &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t2(ji,jj)   &
                  &   ) * r1_e1e2t(ji,jj)
               zdt2 = zdt * zdt

               ! delta at T points
               zdelta(ji,jj) = SQRT( zdiv2 + ( zdt2 + zds2 ) * z1_ecc2 ) * zmsk(ji,jj)        ! zmsk is for reducing cpu

               ! P/delta at T points
               zp_delt(ji,jj) = strength(ji,jj) / ( zdelta(ji,jj) + rn_creepl ) * zmsk(ji,jj) ! zmsk is for reducing cpu

            END DO
         END DO
         !$acc end parallel loop

# if ! defined _OPENACC
         CALL lbc_lnk( 'icedyn_rhg_evp', zdelta, 'T', 1.0_wp, zp_delt, 'T', 1.0_wp )
# endif

         !$acc parallel loop collapse(2) present( zp_delt, zdt_m, zs1, zs2 )
         DO jj = Njs0, Nje0+1
            DO ji = Nis0, Nie0+1   ! loop ends at jpi,jpj so that no lbc_lnk are needed for zs1 and zs2

               ! divergence at T points (duplication to avoid communications)
               ! (brackets added to fix the order of floating point operations for halo 1 - halo 2 compatibility)
               zdiv  = ( (e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj))   &
                  &    + (e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1))   &
                  &    ) * r1_e1e2t(ji,jj)

               ! tension at T points (duplication to avoid communications)
               zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t2(ji,jj)   &
                  &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t2(ji,jj)   &
                  &   ) * r1_e1e2t(ji,jj)

               ! alpha for aEVP
               !   gamma = 0.5*P/(delta+creepl) * (c*pi)**2/Area * dt/m
               !   alpha = beta = sqrt(4*gamma)
               zalph1   = MAX( 50._wp, rpi * SQRT( 0.5_wp * zp_delt(ji,jj) * r1_e1e2t(ji,jj) * zdt_m(ji,jj) ) )
               z1_alph1 = 1._wp / ( zalph1 + 1._wp )
               zalph2   = zalph1
               z1_alph2 = z1_alph1
               ! explicit:
               ! z1_alph1 = 1._wp / zalph1
               ! z1_alph2 = 1._wp / zalph1
               ! zalph1 = zalph1 - 1._wp
               ! zalph2 = zalph1

               ! stress at T points (zkt/=0 if landfast)
               zs1(ji,jj) = ( zs1(ji,jj)*zalph1 + zp_delt(ji,jj) * ( zdiv*(1._wp + zkt) - zdelta(ji,jj)*(1._wp - zkt) ) ) &
                  &         * z1_alph1 * zmsk(ji,jj) ! zmsk is for reducing cpu
               zs2(ji,jj) = ( zs2(ji,jj)*zalph2 + zp_delt(ji,jj) * ( zdt * z1_ecc2 * (1._wp + zkt) ) ) &
                  &         * z1_alph2 * zmsk(ji,jj) ! zmsk is for reducing cpu

            END DO
         END DO
         !$acc end parallel loop

         ! Save beta at T-points for further computations
         !$acc parallel loop collapse(2) present( zbeta )
         DO jj = Njs0-1, Nje0+1
            DO ji = Nis0-1, Nie0+1
               zbeta(ji,jj) = MAX( 50._wp, rpi * SQRT( 0.5_wp * zp_delt(ji,jj) * r1_e1e2t(ji,jj) * zdt_m(ji,jj) ) )
            END DO
         END DO
         !$acc end parallel loop

         !$acc parallel loop collapse(2) present( zs12 )
         DO jj = Njs0-1, Nje0
            DO ji = Nis0-1, Nie0

               ! alpha for aEVP
               zalph2   = MAX( zbeta(ji,jj), zbeta(ji+1,jj), zbeta(ji,jj+1), zbeta(ji+1,jj+1) )
               z1_alph2 = 1._wp / ( zalph2 + 1._wp )
               ! explicit:
               ! z1_alph2 = 1._wp / zalph2
               ! zalph2 = zalph2 - 1._wp

               ! P/delta at F points
               ! (brackets added to fix the order of floating point operations for halo 1 - halo 2 compatibility)
               zp_delf = 0.25_wp * ( (zp_delt(ji,jj) + zp_delt(ji+1,jj)) + (zp_delt(ji,jj+1) + zp_delt(ji+1,jj+1)) )

               ! stress at F points (zkt/=0 if landfast)
               zs12(ji,jj)= ( zs12(ji,jj) * zalph2 + zp_delf * ( zds(ji,jj) * z1_ecc2 * (1._wp + zkt) ) * 0.5_wp ) &
                  &         * z1_alph2

            END DO
         END DO
         !$acc end parallel loop

         ! --- Ice internal stresses (Appendix C of Hunke and Dukowicz, 2002) --- !
         ! (brackets added to fix the order of floating point operations for halo 1 - halo 2 compatibility)
         !$acc parallel loop collapse(2) present( zfU, zfV )
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               !                   !--- U points
               zfU(ji,jj) = 0.5_wp * ( (( zs1(ji+1,jj) - zs1(ji,jj) ) * e2u(ji,jj)                                             &
                  &                  + ( zs2(ji+1,jj) * e2t(ji+1,jj) * e2t(ji+1,jj) - zs2(ji,jj) * e2t2(ji,jj)    &
                  &                    ) * r1_e2u(ji,jj))                                                                      &
                  &                  + ( zs12(ji,jj) * e1f2(ji,jj) - zs12(ji,jj-1) * e1f(ji,jj-1) * e1f(ji,jj-1)  &
                  &                    ) * 2._wp * r1_e1u(ji,jj)                                                              &
                  &                  ) * r1_e1e2u(ji,jj)
               !
               !                !--- V points
               zfV(ji,jj) = 0.5_wp * ( (( zs1(ji,jj+1) - zs1(ji,jj) ) * e1v(ji,jj)                                             &
                  &                  - ( zs2(ji,jj+1) * e1t(ji,jj+1) * e1t(ji,jj+1) - zs2(ji,jj) * e1t2(ji,jj)    &
                  &                    ) * r1_e1v(ji,jj))                                                                      &
                  &                  + ( zs12(ji,jj) * e2f2(ji,jj) - zs12(ji-1,jj) * e2f(ji-1,jj) * e2f(ji-1,jj)  &
                  &                    ) * 2._wp * r1_e2v(ji,jj)                                                              &
                  &                  ) * r1_e1e2v(ji,jj)
               !
            END DO
         END DO
         !$acc end parallel loop

         ! --- Computation of ice velocity --- !
         !  Bouillon et al. 2013 (eq 47-48) => unstable unless alpha, beta vary as in Kimmritz 2016 & 2017
         !  Bouillon et al. 2009 (eq 34-35) => stable

         IF( ln_landfast_L16 ) THEN

            IF( MOD(jter,2) == 0 ) THEN ! even iterations
               !
#include    "icedyn_rhg_evp_updtV_landfast.h90"
               !
#include    "icedyn_rhg_evp_updtU_landfast.h90"
               !
            ELSE ! odd iterations
               !
#include    "icedyn_rhg_evp_updtU_landfast.h90"
               !
#include    "icedyn_rhg_evp_updtV_landfast.h90"
               !
            ENDIF

         ELSE
            
            IF( MOD(jter,2) == 0 ) THEN ! even iterations
               !
#include    "icedyn_rhg_evp_updtV.h90"
               !
#include    "icedyn_rhg_evp_updtU.h90"
               !
            ELSE ! odd iterations
               !
#include    "icedyn_rhg_evp_updtU.h90"
               !
#include    "icedyn_rhg_evp_updtV.h90"
               !
            ENDIF

         ENDIF !IF( ln_landfast_L16 )

         IF( ln_bdy )   CALL bdy_ice_dyn( 'U', u_ice )
         IF( ln_bdy )   CALL bdy_ice_dyn( 'V', v_ice )
         !
         !                                                ! ==================== !
      END DO !DO jter = 1 , nn_nevp                       !  end loop over jter  !
      !                                                   ! ==================== !

      IF(iom_use('fUu')) THEN
         !$acc update self ( zfU )
         CALL iom_put( 'fUu' , zfU )
      ENDIF
      IF(iom_use('fVv')) THEN
         !$acc update self ( zfV )
         CALL iom_put( 'fVv' , zfV )
      ENDIF
      IF(iom_use('beta_evp')) THEN
         !$acc update self ( zbeta )
         CALL iom_put( 'beta_evp' , zbeta )
      ENDIF

      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!
      !$acc parallel loop collapse(2)
      DO jj = Njs0-1, Nje0
         DO ji = Nis0-1, Nie0
            ! shear at F points
            zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f2(ji,jj)   &
               &         + ( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f2(ji,jj)   &
               &         ) * r1_e1e2f(ji,jj) * fmask(ji,jj,1)
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj = Njs0, Nje0
         DO ji = Nis0, Nie0   ! no vector loop
            ! tension**2 at T points
            zdt  = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t2(ji,jj)   &
               &   - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t2(ji,jj)   &
               &   ) * r1_e1e2t(ji,jj)
            zdt2 = zdt * zdt

            zten_i(ji,jj) = zdt

            ! shear**2 at T points (doc eq. A16)
            zds2 = ( zds(ji,jj  ) * zds(ji,jj  ) * e1e2f(ji,jj  ) + zds(ji-1,jj  ) * zds(ji-1,jj  ) * e1e2f(ji-1,jj  )  &
               &   + zds(ji,jj-1) * zds(ji,jj-1) * e1e2f(ji,jj-1) + zds(ji-1,jj-1) * zds(ji-1,jj-1) * e1e2f(ji-1,jj-1)  &
               &   ) * 0.25_wp * r1_e1e2t(ji,jj)

            ! maximum shear rate at T points (includes tension, output only)
            pshear_i(ji,jj) = SQRT( zdt2 + zds2 ) * zmsk(ji,jj) !

            ! shear at T-points
            zshear(ji,jj)   = SQRT( zds2 ) * zmsk(ji,jj)

            ! divergence at T points
            pdivu_i(ji,jj) = ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
               &             + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
               &             ) * r1_e1e2t(ji,jj) * zmsk(ji,jj)

            ! delta at T points
            zdelta(ji,jj)   = SQRT( pdivu_i(ji,jj) * pdivu_i(ji,jj) + ( zdt2 + zds2 ) * z1_ecc2 ) * zmsk(ji,jj) ! delta

            ! delta* at T points (pdelta_i)
            zswitch         = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zdelta(ji,jj) ) ) ! 0 if delta=0
            pdelta_i(ji,jj) = zdelta(ji,jj) + rn_creepl * zswitch
            ! it seems that deformation used for advection and mech redistribution is delta*
            ! MV in principle adding creep limit is a regularization for viscosity not for delta
            ! delta_star should not (in my view) be used in a replacement for delta
         END DO
      END DO
      !$acc end parallel loop

# if ! defined _OPENACC
      CALL lbc_lnk( 'icedyn_rhg_evp', pshear_i, 'T', 1._wp, pdivu_i, 'T', 1._wp, pdelta_i, 'T', 1._wp, zten_i, 'T', 1._wp, &
         &                            zshear  , 'T', 1._wp, zdelta , 'T', 1._wp, zs1     , 'T', 1._wp, zs2   , 'T', 1._wp, zs12, 'F', 1._wp )
      !CALL lbc_lnk( 'icedyn_rhg_evp', zs1, 'T', 1._wp, zs2,'T', 1._wp, zs12, 'F', 1._wp )
# endif

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !! Vertically-integrated components of the T-centric stress tensor (in Pa.m):
            SIGMAt(ji,jj,1) = 0.5_wp*( zs1(ji,jj) + zs2(ji,jj) ) * xmskt(ji,jj)
            SIGMAt(ji,jj,2) = 0.5_wp*( zs1(ji,jj) - zs2(ji,jj) ) * xmskt(ji,jj)
            SIGMAt(ji,jj,3) =                 zs12(ji,jj)        * xmskf(ji,jj)
            ! --- Store the stress tensor for the next time step --- !
            pstress1_i (ji,jj) =  zs1(ji,jj)
            pstress2_i (ji,jj) =  zs2(ji,jj)
            pstress12_i(ji,jj) = zs12(ji,jj)
         END DO
      END DO
      !$acc end parallel loop



      ! 5) diagnostics
      !------------------------------------------------------------------------------!

      ! --- ice/atm. & ice-ocean/bottom(landfast) stresses ---
      IF( ln_landfast_L16 ) THEN
         IF( iom_use('utau_bi') .OR. iom_use('vtau_bi') ) THEN
            CALL lbc_lnk( 'icedyn_rhg_evp', ztaux_bi, 'U', -1._wp, ztauy_bi, 'V', -1._wp )
            CALL iom_put( 'utau_bi' , ztaux_bi * zmsk )
            CALL iom_put( 'vtau_bi' , ztauy_bi * zmsk )
         ENDIF
      ENDIF

      ! --- divergence, shear and strength --- !
      !IF( iom_use('icediv') )   CALL iom_put( 'icediv' , pdivu_i  * zmsk )   ! divergence
      !IF( iom_use('iceshe') )   CALL iom_put( 'iceshe' , pshear_i * zmsk )   ! shear
      IF( iom_use('icestr') ) THEN
         !$acc update self ( strength )
         CALL iom_put( 'icestr' , strength ) !* zmsk )   ! strength
      ENDIF
      IF( iom_use('icedlt') ) THEN
         !$acc update self ( zdelta )
         CALL iom_put( 'icedlt' , zdelta ) !  * zmsk )   ! delta
      ENDIF
      ! --- total deformation of velocity field @T:
      !IF( iom_use('icedeft') ) CALL iom_put( 'icedeft', SQRT( pshear_i*pshear_i + pdivu_i*pdivu_i ) )

      ! --- Stress tensor invariants (SIMIP diags) --- !
      !IF( iom_use('normstr') .OR. iom_use('sheastr') ) THEN
      !   !
      !   ALLOCATE( zsig_I(jpi,jpj) , zsig_II(jpi,jpj) )
      !   !
      !   DO jj = Njs0-1, Nje0+1
      !      DO ji = Nis0-1, Nie0+1
      !         ! Ice stresses
      !         ! sigma1, sigma2, sigma12 are some recombination of the stresses (HD MWR002, Bouillon et al., OM2013)
      !         ! not to be confused with stress tensor components, stress invariants, or stress principal components
      !         zfac             =   strength(ji,jj) / ( zdelta(ji,jj) + rn_creepl )          ! viscosity
      !         zsig1            =   zfac * ( pdivu_i(ji,jj) - zdelta(ji,jj) )
      !         zsig2            =   zfac * z1_ecc2 * zten_i(ji,jj)
      !         zsig12           =   zfac * z1_ecc2 * zshear(ji,jj) * 0.5_wp
      !         ! Stress invariants (sigma_I, sigma_II, Coon 1974, Feltham 2008)
      !         zsig_I (ji,jj)   =   0.5_wp * zsig1
      !         zsig_II(ji,jj)   =   0.5_wp * SQRT ( zsig2 * zsig2 + 4._wp * zsig12 * zsig12 )
      !      END DO
      !   END DO
      !   !
      !   IF( iom_use('normstr') )   CALL iom_put( 'normstr', zsig_I (:,:) * zmsk(:,:) ) ! Normal stress
      !   IF( iom_use('sheastr') )   CALL iom_put( 'sheastr', zsig_II(:,:) * zmsk(:,:) ) ! Maximum shear stress
      !   DEALLOCATE ( zsig_I, zsig_II )
      !ENDIF

      ! --- Normalized stress tensor principal components --- !
      ! This are used to plot the normalized yield curve, see Lemieux & Dupont, 2020
      ! Recommendation 1 : we use ice strength, not replacement pressure
      ! Recommendation 2 : for EVP, no need to use viscosities at last iteration (stress is properly iterated)
      !IF( iom_use('sig1_pnorm') .OR. iom_use('sig2_pnorm') ) THEN
      !   !
      !   ALLOCATE( zsig1_p(jpi,jpj) , zsig2_p(jpi,jpj) , zsig_I(jpi,jpj) , zsig_II(jpi,jpj) )
      !   !
      !   DO jj = Njs0-1, Nje0+1
      !      DO ji = Nis0-1, Nie0+1
      !
      !         ! For EVP solvers, ice stresses at current iterates can be used
      !         !                        following Lemieux & Dupont (2020)
      !         zfac             =   strength(ji,jj) / ( zdelta(ji,jj) + rn_creepl )
      !         zsig1            =   zfac * ( pdivu_i(ji,jj) - zdelta(ji,jj) )
      !         zsig2            =   zfac * z1_ecc2 * zten_i(ji,jj)
      !         zsig12           =   zfac * z1_ecc2 * zshear(ji,jj) * 0.5_wp
      !
      !         ! Stress invariants (sigma_I, sigma_II, Coon 1974, Feltham 2008), T-point
      !         zsig_I(ji,jj)    =   0.5_wp * zsig1                                         ! normal stress
      !         zsig_II(ji,jj)   =   0.5_wp * SQRT ( zsig2 * zsig2 + 4._wp * zsig12 * zsig12 ) ! max shear stress
      !
      !         ! Normalized  principal stresses (used to display the ellipse)
      !         z1_strength      =   1._wp / MAX( 1._wp, strength(ji,jj) )
      !         zsig1_p(ji,jj)   =   ( zsig_I(ji,jj) + zsig_II(ji,jj) ) * z1_strength
      !         zsig2_p(ji,jj)   =   ( zsig_I(ji,jj) - zsig_II(ji,jj) ) * z1_strength
      !END DO
      !   END DO
      !   !
      !   CALL iom_put( 'sig1_pnorm' , zsig1_p * zmsk )
      !   CALL iom_put( 'sig2_pnorm' , zsig2_p * zmsk )
      !
      !   DEALLOCATE( zsig1_p , zsig2_p , zsig_I, zsig_II )
      !
      !ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_dyn_rhg_evp')
      !
   END SUBROUTINE ice_dyn_rhg_evp


   SUBROUTINE rhg_evp_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_evp_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter            ! local integer
      INTEGER  ::   id1, id2, id3   ! local integers
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id1 = iom_varid( numrir, 'stress1_i' , ldstop = .FALSE. )
            id2 = iom_varid( numrir, 'stress2_i' , ldstop = .FALSE. )
            id3 = iom_varid( numrir, 'stress12_i', ldstop = .FALSE. )
            !
            IF( MIN( id1, id2, id3 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'stress1_i' , stress1_i , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'stress2_i' , stress2_i , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'stress12_i', stress12_i, cd_type = 'F' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without rheology, set stresses to 0'
               stress1_i (:,:) = 0._wp
               stress2_i (:,:) = 0._wp
               stress12_i(:,:) = 0._wp
            ENDIF
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set stresses to 0'
            stress1_i (:,:) = 0._wp
            stress2_i (:,:) = 0._wp
            stress12_i(:,:) = 0._wp
         ENDIF
         !
         !$acc update device ( stress1_i, stress2_i, stress12_i )
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- rhg-rst ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         !$acc update self ( stress1_i, stress2_i, stress12_i )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'stress1_i' , stress1_i )
         CALL iom_rstput( iter, nitrst, numriw, 'stress2_i' , stress2_i )
         CALL iom_rstput( iter, nitrst, numriw, 'stress12_i', stress12_i )
         !
      ENDIF
      !
   END SUBROUTINE rhg_evp_rst


   SUBROUTINE ice_dyn_rhg_evp_init( )
      !!-------------------------------------------------------------------
      !! Called into `ice_dyn_rhg_init()@icedyn_rhg.F90`
      !!-------------------------------------------------------------------
      INTEGER, DIMENSION(3) ::   ierr
      INTEGER :: k_alloc
      !!-------------------------------------------------------------------
      IF( lwp ) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) '    *** Initialization of EVP rheology (ice_dyn_rhg_evp_init) ***'
      ENDIF

      ALLOCATE( zdelta(jpi,jpj), zp_delt(jpi,jpj), zbeta(jpi,jpj), zdt_m(jpi,jpj), zmU_dt(jpi,jpj), zmV_dt(jpi,jpj), &
         &      zht(jpi,jpj), zhf(jpi,jpj), zds(jpi,jpj), zten_i(jpi,jpj), zshear(jpi,jpj), zs1(jpi,jpj), zs2(jpi,jpj), zs12(jpi,jpj), &
         &      zsshdyn(jpi,jpj), zfU(jpi,jpj), zfV(jpi,jpj), zspgU(jpi,jpj), zspgV(jpi,jpj),   &
         &      zmsk(jpi,jpj), kmsk01x(jpi,jpj), kmsk01y(jpi,jpj), kmsk00x(jpi,jpj), kmsk00y(jpi,jpj), &
         &      ztaux_base(jpi,jpj), ztauy_base(jpi,jpj), ztaux_bi(jpi,jpj), ztauy_bi(jpi,jpj), &
         &      STAT = ierr(1) )
      !
      zdelta(:,:) = 0._wp ;  zp_delt(:,:) = 0._wp ;  zbeta(:,:) = 0._wp ;  zdt_m(:,:) = 0._wp ;  zmU_dt(:,:) = 0._wp ;  zmV_dt(:,:) = 0._wp
      zht(:,:) = 0._wp ;  zhf(:,:) = 0._wp ;  zds(:,:) = 0._wp ;  zten_i(:,:) = 0._wp ;  zshear(:,:) = 0._wp ;  zs1(:,:) = 0._wp ;  zs2(:,:) = 0._wp ;  zs12(:,:) = 0._wp
      zsshdyn(:,:) = 0._wp ;  zfU(:,:) = 0._wp ;  zfV(:,:) = 0._wp ;  zspgU(:,:) = 0._wp ;  zspgV(:,:) = 0._wp
      zmsk(:,:) = 0._wp
      ztaux_base(:,:) = 0._wp ;  ztauy_base(:,:) = 0._wp ;  ztaux_bi(:,:) = 0._wp ;  ztauy_bi(:,:) = 0._wp

# if defined _OPENACC
      PRINT *, ' * info GPU: ice_dyn_rhg_evp_init() => adding aEVP work arrays to memory!'
      !$acc enter data copyin( zdelta, zp_delt, zbeta, zdt_m, zmU_dt, zmV_dt, zht, zhf, zds )
      PRINT *, '    ==> zdelta, zp_delt, zbeta, zdt_m, zmU_dt, zmV_dt, zht, zhf, zds'
      !$acc enter data copyin( zten_i, zshear, zs1, zs2, zs12, zsshdyn, zfU, zfV, zspgU, zspgV )
      PRINT *, '    ==> zten_i, zshear, zs1, zs2, zs12, zsshdyn, zfU, zfV, zspgU, zspgV'
      !$acc enter data copyin( zmsk, kmsk01x, kmsk01y, kmsk00x, kmsk00y, ztaux_base, ztauy_base, ztaux_bi, ztauy_bi )
      PRINT *, '    ==> zmsk, kmsk01x, kmsk01y, kmsk00x, kmsk00y, ztaux_base, ztauy_base, ztaux_bi, ztauy_bi'
# endif

      !      IF( ln_landfast_L16 ) THEN
      !         PRINT *, 'LOLO_EVP => allocating work landfast arrays!'
      !         ALLOCATE( ,  STAT = ierr(2)  )
      !# if defined _OPENACC
      !         !$acc enter data copyin( ztaux_bi, ztauy_bi )
      !         PRINT *, '    ==> ztaux_bi, ztauy_bi'
      !# endif
      !      ENDIF

      k_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( 'ice_dyn_rhg_evp_init', k_alloc )
      IF( k_alloc > 0 ) CALL ctl_stop('STOP', 'ice_dyn_rhg_evp_init: unable to allocate work arrays for EVP')

   END SUBROUTINE ice_dyn_rhg_evp_init


   !!==============================================================================
END MODULE icedyn_rhg_evp
