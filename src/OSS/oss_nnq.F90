MODULE oss_nnq
   !!======================================================================
   !!                       ***  MODULE  oss_nnq  ***
   !!                       Ocean Surface State arrays
   !!                          => OSS for NANUQ
   !!======================================================================
   !! History :  new for NANUQ
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   oss_nnq_alloc : allocation of OSS arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss_nnq_alloc   ! routine called in ossmod.F90
   PUBLIC   sbcblk_cs_wl_alloc, sbcblk_cs_wl_dealloc

   !!----------------------------------------------------------------------
   !!           Namelist for the Ocean Surface Boundary Condition
   !!----------------------------------------------------------------------
   !                                   !!* namoss namelist *
#if defined key_oasis3
   LOGICAL , PUBLIC ::   lk_oasis_oce = .TRUE.  !: OASIS used    !#lolo?
#else
   LOGICAL , PUBLIC ::   lk_oasis_oce = .FALSE. !: OASIS unused  !#lolo?
#endif
   LOGICAL , PUBLIC ::   ln_prs_oce     !: prescribed surface state of the ocean => standalone formulation
   LOGICAL , PUBLIC ::   ln_cpl_oce     !: ice-ocean coupled formulation
   !
   LOGICAL , PUBLIC ::   ln_ice_embd    !: flag for levitating/embedding sea-ice in the ocean
   !                                             !: =F levitating ice (no presure effect) with mass and salt exchanges
   !                                             !: =T embedded sea-ice (pressure effect + mass and salt exchanges)
   !
   REAL(wp), PARAMETER, PUBLIC :: rdwl0  = 3.    !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)
   !!----------------------------------------------------------------------
   !!           switch definition (improve readability)
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_prs_oce = 1    !: prescribed surface state of the ocean is used as OSS for NANUK
   INTEGER , PUBLIC, PARAMETER ::   jp_cpl_oce = 2    !: coupling with an ocean model is used as OSS for NANUK
   !
   !!----------------------------------------------------------------------
   !!                     Sea Surface Mean fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC                     ::   nn_foss   !: frequency of oss computation
   !$acc declare create( ln_prs_oce, ln_cpl_oce, ln_ice_embd, nn_foss )

   !! Fields read into netCDF file(s) or received from OASIS:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssu_m     !: prescribed or received (coupled) surface sea i-current (U-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssv_m     !: prescribed or received (coupled) surface sea j-current (V-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sst_m     !: prescribed or received (coupled) surface sea temperature (bulk)  [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sss_m     !: prescribed or received (coupled) surface sea salinity            [psu]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssh_m     !: prescribed or received (coupled) sea surface height                [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e3t_m     !: prescribed or received (coupled) sea surface layer thickness       [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   frq_m     !: prescribed or received (coupled) fraction of solar net radiation absorbed in the 1st T level [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   mld_m     !: prescribed mixed layer depth (standalone mode)                     [m]
   ! Following is only for idealized test-cases for 2D advection featuring advection of scalars @ F-points:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssu_v_m   !: prescribed or received (coupled) surface sea i-current (V-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssv_u_m   !: prescribed or received (coupled) surface sea j-current (U-point) [m/s]

   !! The skin temperature in NANUQ, `ssst` should be:
   !! - (A) the temperature of the water directly in contact with the very bottom of the sea-ice layer if sea-ice is present => `t_bo` !?
   !!   or
   !! - (B) the temperature of the water directly in contact with the air (when no sea-ice) if no sea-ice is present
   !!       => see the SBC bulk routines (ECMWF, COARE)...
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssst       !: surface skin temperature of liquid ocean [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sst_s      !: slab bulk surface temperature            [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sss_s      !: slab bulk surface salinity               []
   !                                                                   !:  => basically `sst_m` & `sss_m` corrected when `ln_slab_sst=T`
   !
   !! Cool-skin
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   dT_cs      !: cool-skin temperature anomaly  [K]
   !                                                                   ! => temperature difference between air-sea interface (z=0)
   !                                                                   !    and right below viscous layer (z=delta)
   !! Warm-layer related parameters:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   dT_wl      !: warm-layer temperature anomaly [K]
   !                                                                   ! => difference between "almost surface (right below
   !                                                                   !    viscous layer, z=delta) and depth of bulk SST (z=gdept_1d(1))
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Hz_wl      !: warm-layer thickness           [m]
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: Qnt_ac !: time integral / accumulated heat stored by the warm layer
   !                                                      !         Qxdt => [J/m^2] (reset to zero every midnight)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: Tau_ac !: time integral / accumulated momentum

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: oss_nnq.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION oss_nnq_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION oss_nnq_alloc  ***
      !!---------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: ierr
      !!---------------------------------------------------------------------
      ierr = 0
      !
      IF(lwp) WRITE(numout,*) '  *** oss_nnq_alloc => allocating ssx_m arrays !!!'
      !
      ALLOCATE( ssu_m(jpi,jpj) , sst_m(jpi,jpj) , frq_m(jpi,jpj) , e3t_m(jpi,jpj) ,  &
         &      ssv_m(jpi,jpj) , sss_m(jpi,jpj) , ssh_m(jpi,jpj) ,                   &
         &      ssst(jpi,jpj)  , sst_s(jpi,jpj) , sss_s(jpi,jpj) ,  STAT=ierr(1) )
      !
      oss_nnq_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( 'oss_nnq_alloc', oss_nnq_alloc )
      IF( oss_nnq_alloc > 0 ) CALL ctl_warn('oss_nnq_alloc: allocation of arrays failed')
      !
# if defined _OPENACC
      PRINT *, ' * info GPU: oss_nnq_alloc() => adding OSS arrays to memory!'
      PRINT *, '             => ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m, ssst, sst_s, sss_s'
      !$acc enter data copyin( ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m, ssst, sst_s, sss_s )
# endif
      !
   END FUNCTION oss_nnq_alloc


   INTEGER FUNCTION sbcblk_cs_wl_alloc(l_use_cs, l_use_wl, l_coare)
      !!---------------------------------------------------------------------
      !! INPUT :
      !! -------
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL, INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
      LOGICAL, OPTIONAL, INTENT(in) ::   l_coare  ! we use a COARE algorithm...
      !
      INTEGER, DIMENSION(3) :: ierr
      LOGICAL :: lcoare = .false.
      !!---------------------------------------------------------------------
      IF( PRESENT(l_coare) ) lcoare = l_coare
      ierr(:) = 0
      !
      IF( l_use_wl ) THEN
         ALLOCATE ( dT_wl(jpi,jpj), Hz_wl(jpi,jpj), STAT=ierr(1) )
         dT_wl(:,:)  = 0._wp
         Hz_wl(:,:)  = rdwl0 ! (rdwl0, constant, = 3m is default for Zeng & Beljaars)
         IF( lcoare ) THEN
            ALLOCATE ( Qnt_ac(jpi,jpj), Tau_ac(jpi,jpj), STAT=ierr(2) )
            Qnt_ac(:,:)  = 0._wp
            Tau_ac(:,:)  = 0._wp
         ENDIF
      ENDIF
      IF( l_use_cs ) THEN
         ALLOCATE ( dT_cs(jpi,jpj), STAT=ierr(3) )
         dT_cs(:,:) = -0.25_wp  ! First guess of skin correction
      ENDIF

      sbcblk_cs_wl_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( 'sbcblk_cs_wl_alloc', sbcblk_cs_wl_alloc )
      IF( sbcblk_cs_wl_alloc > 0 ) CALL ctl_warn('sbcblk_cs_wl_alloc: allocation of arrays failed')

   END FUNCTION sbcblk_cs_wl_alloc


   SUBROUTINE sbcblk_cs_wl_dealloc( l_use_wl )
      !!---------------------------------------------------------------------
      !! INPUT :
      !! -------
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      LOGICAL , INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
      INTEGER :: ierr
      !!---------------------------------------------------------------------
      IF( l_use_wl ) THEN
         ierr = 0
         DEALLOCATE ( dT_wl, Hz_wl, STAT=ierr )
         IF( ierr > 0 ) CALL ctl_stop( 'sbcblk_cs_wl_dealloc => deallocation of dT_wl & Hz_wl failed!' )
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   DEALLOCATE ( dT_cs, STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( 'sbcblk_cs_wl_dealloc => deallocation of dT_cs failed!' )
      !ENDIF
   END SUBROUTINE sbcblk_cs_wl_dealloc



   !!======================================================================
END MODULE oss_nnq
