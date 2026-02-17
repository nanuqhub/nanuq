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


   !!--------------------------------------------------
   !! *** Skin arrays (always used) ***
   !!--------------------------------------------------
   !! The skin temperature in NANUQ, `ssst` should be:
   !! - (A) the temperature of the water directly in contact with the very bottom of the sea-ice layer if sea-ice is present => `t_bo` !?
   !!   or
   !! - (B) the temperature of the water directly in contact with the air (when no sea-ice) if no sea-ice is present
   !!       => see the SBC bulk routines (ECMWF, COARE)...
   !!
   !!  * Arrays `ssst` => contains the skin SST if "cool-skin or/and warm-layer" params are used
   !!                  => contains a simple copy of the bulk SST otherwize
   !!                  => `ssst` remains in degree Celsius !!!
   !!
   !!      ==> same for `sssq`, which is `q_sat(ssst)` in kg/kg
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssst       !: surface skin temperature of liquid ocean [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sssq       !: q_sat(ssst)                              [kg/kq]


   !!--------------------------------------------------
   !! *** SlabOcean stuff ***
   !!--------------------------------------------------
   !! Here we are dealing with a modification of the prescribed bulk SST and SSS
   !! based on consideration of the net heat & freshwater fluxes received by the
   !! ocean and "absorbed" by a layer corresponding to the prescribed MLD !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sst_s      !: slab bulk surface temperature            [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sss_s      !: slab bulk surface salinity               []
   !                                                                   !:  => basically `sst_m` & `sss_m` corrected when `ln_slab_sst=T`

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
      IF(lwp) WRITE(numout,*) '  *** oss_nnq_alloc => allocating `ss*` arrays !!!'
      !
      ALLOCATE( ssu_m(jpi,jpj) , sst_m(jpi,jpj) , frq_m(jpi,jpj) , e3t_m(jpi,jpj) ,  &
         &      ssv_m(jpi,jpj) , sss_m(jpi,jpj) , ssh_m(jpi,jpj) ,                   &
         &      ssst(jpi,jpj)  , sssq(jpi,jpj)  ,                                    &
         &      sst_s(jpi,jpj) , sss_s(jpi,jpj) ,  STAT=ierr(1) )
      !
      oss_nnq_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( 'oss_nnq_alloc', oss_nnq_alloc )
      IF( oss_nnq_alloc > 0 ) CALL ctl_warn('oss_nnq_alloc: allocation of arrays failed')
      !
# if defined _OPENACC
      PRINT *, ' * info GPU: oss_nnq_alloc() => adding OSS arrays to memory!'
      PRINT *, '             => ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m, ssst, sssq, sst_s, sss_s'
      !$acc enter data copyin( ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m, ssst, sssq, sst_s, sss_s )
# endif
      !
   END FUNCTION oss_nnq_alloc


   !!======================================================================
END MODULE oss_nnq
