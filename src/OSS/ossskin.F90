MODULE ossskin
   !!======================================================================
   !!                       ***  MODULE  ossskin  ***
   !!   Consideration of Sea Surface Temperature Skin: cool skin and/or warm layer
   !!           in the computation of turbulent fluxes with bulk formulas
   !!
   !!   MIND: ==> ONLY relevant when ECMWF or COAREX bulk algos are used !
   !!             when these are not used ssst == sst_m & sssq = q_sat(sst_m) !
   !!
   !!======================================================================
   !! History :  new for NANUQ => Laurent Brodeau 2026
   !!
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   !PUBLIC  oss_skin_init

   PUBLIC  oss_skin_alloc
   PUBLIC  oss_skin_dealloc


   !!--------------------------------------------------
   !! *** Skin temperature stuff ***
   !!--------------------------------------------------
   !!
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
   !! `ssst` & `ssss` declared and allocated in `oss_nnq.F90`...
   !!
   REAL(wp), PARAMETER, PUBLIC :: rdwl0  = 3.    !: Depth scale [m] of warm layer, "d" in Eq.11 (Zeng & Beljaars 2005)

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
   !! $Id: ossskin.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !
   !# if defined _OPENACC
   !      PRINT *, ' * info GPU: ossskin_alloc() => adding OSS arrays to memory!'
   !      PRINT *, '             => ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m, sst_s, sss_s'
   !      !$acc enter data copyin( ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m,  sst_s, sss_s )
   !# endif


   INTEGER FUNCTION oss_skin_alloc( l_use_cs, l_use_wl, l_coare )
      !!---------------------------------------------------------------------
      !! INPUT :
      !! -------
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!---------------------------------------------------------------------
      LOGICAL, OPTIONAL, INTENT(in) ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL, OPTIONAL, INTENT(in) ::   l_use_wl ! use the warm-layer parameterization
      LOGICAL, OPTIONAL, INTENT(in) ::   l_coare  ! we use a COARE algorithm...
      !!---------------------------------------------------------------------
      INTEGER, DIMENSION(4) :: ierr
      LOGICAL :: llcs, llwl, llcoare
      !!---------------------------------------------------------------------
      llcs = .FALSE. ; llwl = .FALSE. ; llcoare = .FALSE.

      IF( PRESENT( l_use_cs) ) llcs    = l_use_cs
      IF( PRESENT( l_use_wl) ) llwl    = l_use_wl
      IF( PRESENT( l_coare ) ) llcoare = l_coare


      ierr(:) = 0

      !      ALLOCATE( ssst(jpi,jpj)  , sssq(jpi,jpj)  ,  STAT=ierr(1)  )


      IF( llwl ) THEN
         ALLOCATE ( dT_wl(jpi,jpj), Hz_wl(jpi,jpj), STAT=ierr(2) )
         dT_wl(:,:)  = 0._wp
         Hz_wl(:,:)  = rdwl0 ! (rdwl0, constant, = 3m is default for Zeng & Beljaars)
         IF( llcoare ) THEN
            ALLOCATE ( Qnt_ac(jpi,jpj), Tau_ac(jpi,jpj), STAT=ierr(3) )
            Qnt_ac(:,:)  = 0._wp
            Tau_ac(:,:)  = 0._wp
         ENDIF
      ENDIF
      IF( llcs ) THEN
         ALLOCATE ( dT_cs(jpi,jpj), STAT=ierr(4) )
         dT_cs(:,:) = -0.25_wp  ! First guess of skin correction
      ENDIF


      !# if defined _OPENACC
      !      PRINT *, ' * info GPU: oss_skin_alloc() => adding skin SST arrays to memory!'
      !      PRINT *, '             => ssst, sssq'
      !      !$acc enter data copyin( ssst, sssq )
      !# endif


      oss_skin_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( 'oss_skin_alloc', oss_skin_alloc )
      !IF( oss_skin_alloc > 0 ) CALL ctl_warn('oss_skin_alloc: allocation of arrays failed')

   END FUNCTION oss_skin_alloc


   SUBROUTINE oss_skin_dealloc( l_use_wl )
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
         IF( ierr > 0 ) CALL ctl_stop( 'oss_skin_dealloc => deallocation of dT_wl & Hz_wl failed!' )
      ENDIF
      !IF( l_use_cs ) THEN
      !   ierr = 0
      !   DEALLOCATE ( dT_cs, STAT=ierr )
      !   IF( ierr > 0 ) CALL ctl_stop( 'oss_skin_dealloc => deallocation of dT_cs failed!' )
      !ENDIF
   END SUBROUTINE oss_skin_dealloc


   !!======================================================================
END MODULE ossskin
