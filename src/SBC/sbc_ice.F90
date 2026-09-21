MODULE sbc_ice
   !!======================================================================
   !!                 ***  MODULE  sbc_ice  ***
   !! Surface module - NANUQ: parameters & variables defined in memory
   !!======================================================================
   !! History :  3.0   !  2006-08  (G. Madec)        Surface module
   !!            3.2   !  2009-06  (S. Masson)       merge with ice_oce
   !!            3.3.1 !  2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4   !  2011-11  (C. Harris)       CICE added as an option
   !!            4.0   !  2018     (many people)     SI3 compatibility
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE par_ice, ONLY : jpl
   USE lib_mpp, ONLY : mpp_sum, ctl_warn
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ice_alloc   ! called in sbcmod.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qns_ice        !: non solar heat flux over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice        !: solar heat flux over ice                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qlw_ice        !: IR (longwave) heat flux over ice              [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qla_ice        !: latent heat flux over ice                     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsb_ice        !: sensible heat flux over ice                   [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqns_ice       !: non solar heat flux over ice (LW+SEN+LA)      [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   alb_ice        !: ice albedo                                       [-]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qml_ice        !: heat available for snow / ice surface melting     [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qcn_ice        !: heat conduction flux in the layer below surface   [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qtr_ice_top    !: solar flux transmitted below the ice surface      [W/m2]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   taux_ai_t       !: atmos-ice u-stress. T-pts                  [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   tauy_ai_t       !: atmos-ice v-stress. T-pts                  [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice        !: sublimation - precip over sea ice          [kg/m2/s]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt            !: category topmelt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   botmelt            !: category botmelt

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   evap_ice       !: sublimation of ice to the atmo (<0 when ice losing FW to atmo) [kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qemp_ice       !: heat flux of precip and evap over ice       [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qprec_ice      !: enthalpy of precip over ice                 [J/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   CD_ice, CE_ice, CH_ice !: bulk transfer coefficients over sea-ice [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   t_air_zu_i     !: air pot. temperature (adjusted) at wind height (10m)  over sea-ice [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   q_air_zu_i     !: air specific humidity (adjusted) at wind height (10m) over sea-ice [kg/kg]

   !! arrays relating to embedding ice in the ocean
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass        !: mass of snow and ice at current  ice time step   [Kg/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass_b      !: mass of snow and ice at previous ice time step   [Kg/m2]

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_ice_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION sbc_ice_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(5), ii
      !!----------------------------------------------------------------------
      ierr(:) = 0
      ii = 0

      ii = ii + 1
      ALLOCATE( snwice_mass(jpi,jpj) , snwice_mass_b(jpi,jpj) ,   STAT=ierr(ii) )

      ii = ii + 1
      ALLOCATE( taux_ai_t(jpi,jpj) , tauy_ai_t(jpi,jpj) , STAT= ierr(ii) )
      !   &      rCdU_ice(jpi,jpj)                     , STAT= ierr(ii) )


      ii = ii + 1
      ALLOCATE( t_air_zu_i(jpi,jpj)    ,  q_air_zu_i(jpi,jpj),        &
         &      CD_ice(jpi,jpj) , CE_ice(jpi,jpj) , CH_ice(jpi,jpj) ,                         &
         &      qns_ice (jpi,jpj,jpl) , qsr_ice  (jpi,jpj,jpl) ,     qlw_ice(jpi,jpj,jpl) ,   &
         &      qla_ice (jpi,jpj,jpl) , qsb_ice  (jpi,jpj,jpl) ,      &
         &      dqns_ice(jpi,jpj,jpl) , alb_ice  (jpi,jpj,jpl) ,   &
         &      qml_ice (jpi,jpj,jpl) , qcn_ice  (jpi,jpj,jpl) , qtr_ice_top(jpi,jpj,jpl) ,   &
         &      evap_ice(jpi,jpj,jpl) , qprec_ice(jpi,jpj)     ,   &
         &      qemp_ice(jpi,jpj)     , emp_ice (jpi,jpj) , STAT= ierr(ii) )
      !
      sbc_ice_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbc_ice', sbc_ice_alloc )
      IF( sbc_ice_alloc > 0 )   CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
      !
      t_air_zu_i(:,:)=0._wp    ; q_air_zu_i(:,:)=0._wp
      CD_ice(:,:)=1._wp   ;  CE_ice(:,:)=1._wp   ;  CH_ice(:,:)=1._wp   ;
      qns_ice (:,:,:)=0._wp ; qsr_ice  (:,:,:)=0._wp ;     qlw_ice(:,:,:)=0._wp
      qla_ice (:,:,:)=0._wp ; qsb_ice  (:,:,:)=0._wp
      dqns_ice(:,:,:)=0._wp ; alb_ice  (:,:,:)=0._wp
      qml_ice (:,:,:)=0._wp ; qcn_ice  (:,:,:)=0._wp ; qtr_ice_top(:,:,:)=0._wp
      evap_ice(:,:,:)=0._wp ; qprec_ice(:,:)=0._wp
      qemp_ice(:,:)=0._wp   ; emp_ice  (:,:)=0._wp
      !
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: sbc_ice_alloc() => adding SBC-related arrays to memory!'
      PRINT *, '            => taux_ai_t, tauy_ai_t, t_air_zu_i, q_air_zu_i, CD_ice, CE_ice, CH_ice'
      !$acc enter data copyin( taux_ai_t, tauy_ai_t, t_air_zu_i, q_air_zu_i, CD_ice, CE_ice, CH_ice )
      PRINT *, '            => qns_ice, qsr_ice, qlw_ice, qla_ice, qsb_ice, qsb_ice, dqns_ice, alb_ice'
      !$acc enter data copyin( qns_ice, qsr_ice, qlw_ice, qla_ice, qsb_ice, qsb_ice, dqns_ice, alb_ice )
      PRINT *, '            => qml_ice, qcn_ice, qtr_ice_top, evap_ice, qprec_ice'
      !$acc enter data copyin( qml_ice, qcn_ice, qtr_ice_top, evap_ice, qprec_ice )
      PRINT *, '            => qemp_ice, emp_ice, snwice_mass, snwice_mass_b'
      !$acc enter data copyin( qemp_ice, emp_ice, snwice_mass, snwice_mass_b )
#endif

   END FUNCTION sbc_ice_alloc


   !!======================================================================
END MODULE sbc_ice
