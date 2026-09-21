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
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_oce_alloc   ! routine called in sbcmod.F90

   !!----------------------------------------------------------------------
   !!           Namelist for the Ocean Surface Boundary Condition
   !!----------------------------------------------------------------------
   !                                   !!* namsbc namelist *
   LOGICAL , PUBLIC ::   ln_flx         !: flux      formulation
   LOGICAL , PUBLIC ::   ln_blk         !: bulk formulation
   LOGICAL , PUBLIC ::   ln_abl         !: Atmospheric boundary layer model
   !
   CHARACTER(len=1) , PUBLIC ::  sn_loc_vct_tau  !: C-grid point location of surface (air-sea or ice-sea) stress to be computed and/or transmitted to OASIS
   !                                             !:   * 'T'   => both `utau` & `vtau` are located at T-point
   !                                             !:   * 'C'   => follow 'C-grid' convention => `utau` @ U-point & `vtau` @ V-point
   INTEGER(1), PUBLIC :: k_tau_air_at_T
   !$acc declare create( k_tau_air_at_T )
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
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau   , utau_b   !: sea surface i-stress (ocean referential) @ T [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vtau   , vtau_b   !: sea surface j-stress (ocean referential) @ T [N/m2]
   !! IMPORTANT: about `utau` & `vtau`:
   !!            regardless of the value of `sn_loc_vct_tau`, `utau` and `vtau` computed for ice-free ocean (IOMed as
   !!            `utau_oce` & `vtau_oce`) ARE ALWAYS DEFINED AT T POINTS WHEN EXITING `SBC()`!!!
   !!            (because prescribed atmospheric wind is generally read at T-points as well).
   !!   However, later on, when updated for transmission to the ocean component (or simple diagnostic in standalone mode),
   !!   taking into account the potential contribution of sea ice, in `ice_update_tau@iceupdate.F90`, they are interpolated
   !!   to U- & V-point if `sn_loc_vct_tau==C` or kept at T-points if `sn_loc_vct_tau==T`.
   !!   => for example, when coupled to NEMOv4.2.2, which expects `utau` and `vtau` received from OASIS to be @U & @V,
   !!      one should USE `sn_loc_vct_tau==C` !
   !!
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau_icb, vtau_icb !: sea surface (i,j)-stress used by icebergs   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   taum              !: module of sea surface stress (at T-point)    [N/m2]
   !! wndm is used compute surface gases exchanges in ice-free ocean or leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wndm              !: wind speed module at T-point (=|U10m-Uoce|)  [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rhoa              !: air density at "rn_zu" m above the sea       [kg/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr    , qsr_b    !: sea heat flux:     solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns    , qns_b    !: sea heat flux: non solar                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp    , emp_b    !: freshwater budget: volume flux (>0 => LOSS for liquid ocean)      [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx    , sfx_b    !: salt flux                                                         [PSS.kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fmmflx            !: freshwater budget: freezing/melting (>0 => LOSS for liquid ocean) [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   t_air_zu          !: air potential temperature (adjusted) at wind height (10m) [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   q_air_zu          !: air specific humidity (adjusted)     at wind height (10m) [kg/kg]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns_oce           !: non solar heat flux over ice-free ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr_oce           !: non solar heat flux over ice-free ocean              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qemp_oce          !: heat flux of precip and evap over ice-free ocean     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp_oce           !: evap - precip over ice-free ocean                 [kg/m2/s]

   !!---------------------------------------------------------------------
   !! ABL Vertical Domain size
   !!---------------------------------------------------------------------
   INTEGER , PUBLIC            ::   jpka   = 2     !: ABL number of vertical levels (default definition)
   INTEGER , PUBLIC            ::   jpkam1 = 1     !: jpka-1
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   ght_abl, ghw_abl          !: ABL geopotential height (needed for iom)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_abl, e3w_abl          !: ABL vertical scale factors (needed for iom)

   !!----------------------------------------------------------------------
   !!                     Surface atmospheric fields
   !!----------------------------------------------------------------------
   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_theta, fatm_q, fatm_slp, fatm_wnd, fatm_u, fatm_v
   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_prcp, fatm_snow !: total & snow precipitation    [Kg/m2/s]
   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fatm_dqsw, fatm_dqlw !: downwelling short- and long-wave radiation [W/m2]


   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! $Id: sbc_oce.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_oce_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION sbc_oce_alloc  ***
      !!---------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!---------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( utau(jpi,jpj) , utau_b(jpi,jpj) , taum(jpi,jpj) ,     &
         &      vtau(jpi,jpj) , vtau_b(jpi,jpj) , wndm(jpi,jpj) ,     &
         &      rhoa(jpi,jpj) , t_air_zu(jpi,jpj), q_air_zu(jpi,jpj),  STAT=ierr(1) )
      !
      ALLOCATE( qns(jpi,jpj) , qns_b(jpi,jpj) , qsr(jpi,jpj) , qsr_b(jpi,jpj),  &
         &      emp    (jpi,jpj) , emp_b(jpi,jpj) ,                        &
         &      sfx    (jpi,jpj) , sfx_b(jpi,jpj), fmmflx(jpi,jpj), STAT=ierr(2) )

      ALLOCATE( qns_oce(jpi,jpj), qsr_oce(jpi,jpj), qemp_oce(jpi,jpj), emp_oce(jpi,jpj), STAT=ierr(3))

      sbc_oce_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbc_oce', sbc_oce_alloc )
      IF( sbc_oce_alloc > 0 )   CALL ctl_warn('sbc_oce_alloc: allocation of arrays failed')

      utau(:,:)=0._wp ; utau_b(:,:)=0._wp ; taum(:,:)=0._wp
      vtau(:,:)=0._wp ; vtau_b(:,:)=0._wp ; wndm(:,:)=0._wp
      rhoa(:,:)=0._wp ; t_air_zu(:,:)=0._wp ; q_air_zu(:,:)=0._wp
      !
      qns(:,:)=0._wp ; qns_b(:,:)=0._wp
      qsr(:,:)=0._wp ; qsr_b(:,:)=0._wp
      emp(:,:)=0._wp ; emp_b(:,:)=0._wp
      sfx(:,:)=0._wp ; sfx_b(:,:)=0._wp ; fmmflx(:,:)=0._wp
      !
      qns_oce(:,:)=0._wp ; qsr_oce(:,:)=0._wp ; qemp_oce(:,:)=0._wp ; emp_oce(:,:)=0._wp
      !
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: sbc_oce_alloc() => adding SBC-related arrays to memory!'
      PRINT *, '            => qns, qsr, emp, sfx, fmmflx, qns_b, qsr_b, emp_b, sfx_b'
      !$acc enter data copyin( qns, qsr, emp, sfx, fmmflx, qns_b, qsr_b, emp_b, sfx_b )
      PRINT *, '            => qsr, wndm, taum, rhoa, utau, vtau, utau_b, vtau_b, t_air_zu, q_air_zu'
      !$acc enter data copyin( qsr, wndm, taum, rhoa, utau, vtau, utau_b, vtau_b, t_air_zu, q_air_zu )
      PRINT *, '            => qns_oce, qsr_oce, qemp_oce, emp_oce'
      !$acc enter data copyin( qns_oce, qsr_oce, qemp_oce, emp_oce )
#endif
      !
   END FUNCTION sbc_oce_alloc

   !!======================================================================
END MODULE sbc_oce
