MODULE sbc_oce
   !!======================================================================
   !!                       ***  MODULE  sbc_oce  ***
   !! Surface module :  provides air-sea fluxes over liquid water
   !!======================================================================
   !!
   !! History :  3.0  ! 2006-06  (G. Madec)  Original code
   !!             -   ! 2008-08  (G. Madec)  namsbc moved from sbcmod
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!            3.3  ! 2010-10  (J. Chanut, C. Bricaud)  add the surface pressure forcing
   !!            4.0  ! 2012-05  (C. Rousset) add attenuation coef for use in ice model
   !!            4.0  ! 2016-06  (L. Brodeau) new unified bulk routine (based on AeroBulk)
   !!            4.0  ! 2019-03  (F. Lemarié, G. Samson) add compatibility with ABL mode
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) modified wave parameters in namelist
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_oce_alloc : allocation of sbc arrays
   !!   sbc_tau2wnd   : wind speed estimated from wind stress
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_oce_alloc   ! routine called in sbcmod.F90
   PUBLIC   sbc_tau2wnd     ! routine called in several sbc modules

   !!----------------------------------------------------------------------
   !!           Namelist for the Ocean Surface Boundary Condition
   !!----------------------------------------------------------------------
   !                                   !!* namsbc namelist *
   LOGICAL , PUBLIC ::   ln_flx         !: flux      formulation
   LOGICAL , PUBLIC ::   ln_blk         !: bulk formulation
   LOGICAL , PUBLIC ::   ln_abl         !: Atmospheric boundary layer model
   LOGICAL , PUBLIC ::   ln_wave        !: wave in the system (forced or coupled)
   !
   !LOLO: coupling with atmosphere disabled for now
   !#if defined key_oasis3
   !   LOGICAL , PUBLIC ::   lk_oasis_atm = .TRUE.  !: OASIS used
   !#else
   LOGICAL , PUBLIC ::   lk_oasis_atm = .FALSE. !: OASIS unused
   !#endif
   !
   LOGICAL , PUBLIC ::   ln_cpl_atm     !: ocean/sea-ice - atmosphere coupled formulation
   LOGICAL , PUBLIC ::   ln_dm2dc       !: Daily mean to Diurnal Cycle short wave (qsr)
   !LOGICAL , PUBLIC ::   ln_icebergs    !: Icebergs
   !
   INTEGER , PUBLIC ::   nn_lsm         !: Number of iteration if seaoverland is applied
   !
   !                                   !!* namsbc_cpl namelist *
   INTEGER , PUBLIC ::   nn_cats_cpl    !: Number of sea ice categories over which the coupling is carried out
   !
   !!----------------------------------------------------------------------
   !!           switch definition (improve readability)
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_flx = 2        !: flux                          formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_blk = 3        !: bulk                          formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_abl = 4        !: Atmospheric boundary layer    formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_cpl_atm = 5        !: Pure ocean-atmosphere Coupled formulation
   !
   !!----------------------------------------------------------------------
   !!              Ocean Surface Boundary Condition fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC ::  ncpl_qsr_freq = 0        !: qsr coupling frequency per days from atmosphere (used by top)
   !
   !!                                   !!   now    ! before   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau   , utau_b   !: sea surface i-stress (ocean referential)     [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vtau   , vtau_b   !: sea surface j-stress (ocean referential)     [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau_icb, vtau_icb !: sea surface (i,j)-stress used by icebergs   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   taum              !: module of sea surface stress (at T-point)    [N/m2]
   !! wndm is used compute surface gases exchanges in ice-free ocean or leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wndm              !: wind speed module at T-point (=|U10m-Uoce|)  [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rhoa              !: air density at "rn_zu" m above the sea       [kg/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr               !: sea heat flux:     solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns    , qns_b    !: sea heat flux: non solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr_tot           !: total     solar heat flux (over sea and ice) [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns_tot           !: total non solar heat flux (over sea and ice) [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp    , emp_b    !: freshwater budget: volume flux               [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx    , sfx_b    !: salt flux                                    [PSS.kg/m2/s]
   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp_tot           !: total E-P over ocean and ice                 [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fmmflx            !: freshwater budget: freezing/melting          [Kg/m2/s]
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fr_i              !: ice fraction = 1 - lead fraction      (between 0 to 1)

   !!---------------------------------------------------------------------
   !! ABL Vertical Domain size
   !!---------------------------------------------------------------------
   INTEGER , PUBLIC            ::   jpka   = 2     !: ABL number of vertical levels (default definition)
   INTEGER , PUBLIC            ::   jpkam1 = 1     !: jpka-1
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   ght_abl, ghw_abl          !: ABL geopotential height (needed for iom)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_abl, e3w_abl          !: ABL vertical scale factors (needed for iom)

   !!----------------------------------------------------------------------
   !!                     Sea Surface Mean fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC                     ::   nn_fsbc   !: frequency of sbc computation (as well as sea-ice model)

   !!----------------------------------------------------------------------
   !!                     Surface atmospheric fields
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_theta, fatm_q, fatm_slp, fatm_wnd, fatm_u, fatm_v
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_prcp, fatm_snow !: total & snow precipitation    [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_dqsw, fatm_dqlw !: downwelling short- and long-wave radiation [W/m2]


   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: sbc_oce.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_oce_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION sbc_oce_alloc  ***
      !!---------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!---------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( utau(jpi,jpj) , utau_b(jpi,jpj) , taum(jpi,jpj) ,     &
         &      vtau(jpi,jpj) , vtau_b(jpi,jpj) , wndm(jpi,jpj) , rhoa(jpi,jpj) , STAT=ierr(1) )
      !
      ALLOCATE( qns_tot(jpi,jpj) , qns  (jpi,jpj) , qns_b(jpi,jpj),        &
         &      qsr_tot(jpi,jpj) , qsr  (jpi,jpj) ,                        &
         &      emp    (jpi,jpj) , emp_b(jpi,jpj) ,                        &
         &      sfx    (jpi,jpj) , sfx_b(jpi,jpj), fmmflx(jpi,jpj), STAT=ierr(2) )
      !  , emp_tot(jpi,jpj)
      !
      ALLOCATE(  fr_i(jpi,jpj), STAT=ierr(3) )
      !
      ALLOCATE( fatm_theta(jpi,jpj), fatm_q(jpi,jpj), fatm_slp(jpi,jpj), fatm_wnd(jpi,jpj), &
         &      fatm_prcp(jpi,jpj),  fatm_snow(jpi,jpj), fatm_u(jpi,jpj), fatm_v(jpi,jpj),  &
         &      fatm_dqsw(jpi,jpj),  fatm_dqlw(jpi,jpj),  STAT=ierr(4) ) !#LB
      !
      sbc_oce_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbc_oce', sbc_oce_alloc )
      IF( sbc_oce_alloc > 0 )   CALL ctl_warn('sbc_oce_alloc: allocation of arrays failed')
      !
# if defined _OPENACC
      PRINT *, ' * info GPU: sbc_oce_alloc() => adding SBC-related arrays to memory!'
      PRINT *, '            => qns, qsr, emp, sfx, fmmflx, qns_tot, qsr_tot, qns_b, emp_b, sfx_b'
      !$acc enter data copyin( qns, qsr, emp, sfx, fmmflx, qns_tot, qsr_tot, qns_b, emp_b, sfx_b )
      PRINT *, '            => qsr, wndm, taum, fr_i, rhoa, utau, vtau, utau_b, vtau_b'
      !$acc enter data copyin( qsr, wndm, taum, fr_i, rhoa, utau, vtau, utau_b, vtau_b )
      PRINT *, '            => fatm_theta, fatm_q, fatm_slp, fatm_wnd, fatm_u, fatm_v, fatm_prcp, fatm_snow, fatm_dqsw, fatm_dqlw'
      !$acc enter data copyin( fatm_theta, fatm_q, fatm_slp, fatm_wnd, fatm_u, fatm_v, fatm_prcp, fatm_snow, fatm_dqsw, fatm_dqlw )
# endif
      !
   END FUNCTION sbc_oce_alloc


   SUBROUTINE sbc_tau2wnd
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_tau2wnd  ***
      !!
      !! ** Purpose : Estimation of wind speed as a function of wind stress
      !!
      !! ** Method  : |tau|=rhoa*Cd*|U|^2
      !!---------------------------------------------------------------------
      USE dom_oce         ! ocean space and time domain
      USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, ztau, zcoef ! temporary variables
      INTEGER  ::   ji, jj                ! dummy indices
      !!---------------------------------------------------------------------
      zcoef = 0.5 / ( zrhoa * zcdrag )
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            ztx = utau(ji-1,jj  ) + utau(ji,jj)
            zty = vtau(ji  ,jj-1) + vtau(ji,jj)
            ztau = SQRT( ztx * ztx + zty * zty )
            wndm(ji,jj) = SQRT ( ztau * zcoef ) * xmskt(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( 'sbc_oce', wndm(:,:) , 'T', 1.0_wp )
      !
   END SUBROUTINE sbc_tau2wnd

   !!======================================================================
END MODULE sbc_oce
