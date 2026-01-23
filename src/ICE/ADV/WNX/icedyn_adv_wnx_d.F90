MODULE icedyn_adv_wnx_d
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_wnx_d   ***
   !!   sea-ice : advection => WENOX scheme
   !!======================================================================
   !! History :  NANUQ 0.1 !  2025-02 Laurent Brodeau, original code
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE par_ice, ONLY: rDt_ice
   USE iom            ! I/O manager library
   !
   USE icedyn_adv_wnx_adv, ONLY : wenoX_rk3
   !
# if defined _OPENACC
   USE lbclnk_gpu
# else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
# endif
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv_wnx_d   ! called by icedyn_adv

   !!----------------------------------------------------------------------
   !! NANUQ_beta
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ice_dyn_adv_wnx_d( kt, cgt, pe1e2, p1_e1e2, pmsk, pUx, pVx, pmlbc,  p1md )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_wnx_d  **
      !!
      !! ** BRODEAU, 2025 => Advect `damage` at T- & F- points
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!                but also their quadratic forms. The method preserves
      !!                tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in   ) :: kt         ! current time step
      CHARACTER(len=1),             INTENT(in   ) :: cgt        ! 'T' or 'F' points ???
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pe1e2, p1_e1e2 ! at the relevant point (T or F) ! ALREADY IN KM^2 !!!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pmsk       ! at the relevant point (T or F)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pUx     ! ice i-velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pVx     ! ice j-velocity
      INTEGER(1), DIMENSION(jpi,jpj,nn_hls,4), INTENT(in ) ::   pmlbc      ! masks for solid LBCs
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: p1md       ! (1 - damage) @ T
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj
      REAL(wp) :: zdt, zr
      CHARACTER(len=64) :: cstr
      !
      CHARACTER(len=19), PARAMETER :: crtnm = 'ice_dyn_adv_wnx_d'
      REAL(wp), PARAMETER :: rr_scl_fct = 1.E-6_wp
      REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start(crtnm)
      !$acc data present( pe1e2, p1_e1e2, pmsk, pUx, pVx, pmlbc, p1md )

      cstr = 'WENOX advection scheme for BBM (damage only)'
      IF( kt == nit000 .AND. lwp ) WRITE(numout,*) '-- '//crtnm//': '//TRIM(cstr)//' at '//cgt//'-points'

      zdt = rDt_ice

      !LOLO: as oposed to `ice_dyn_adv_wnx`, here input arrays to advect seem to be `lbc_lnk`ed
      !      => so no `lbc_lnk`ing of `p1md` !!!
      ! ==> REALLY ???
#if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', p1md )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', p1md )
#else
      CALL lbc_lnk( crtnm, p1md,cgt,1._wp )
#endif


      !! ADVECT !

      CALL wenoX_rk3( kt, cgt, zdt, pe1e2, p1_e1e2, pUx, pVx, pmlbc, p1md )
      
      !$acc end data
      IF( ln_timing )   CALL timing_stop(crtnm)
      !
   END SUBROUTINE ice_dyn_adv_wnx_d

   !!======================================================================
END MODULE icedyn_adv_wnx_d
