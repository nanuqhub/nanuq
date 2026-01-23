MODULE icethd_zdf
   !!======================================================================
   !!                       ***  MODULE icethd_zdf ***
   !!   sea-ice: master routine for vertical heat diffusion in sea ice
   !!======================================================================
   !! History :  4.0  !  2018     (C. Rousset)      Original code SI3
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  ice_thd_zdf      : select the appropriate routine for vertical heat diffusion calculation
   !!  ice_thd_zdf_BL99 : heat diffusion from Bitz and Lipscomb 1999
   !!  ice_thd_zdf_init : initialization
   !!----------------------------------------------------------------------
   USE par_ice
   USE par_oce
   USE par_kind , ONLY : wp
   USE phycst   , ONLY : rcnd_s
   USE icethd_zdf_BL99 ! sea-ice: vertical diffusion (Bitz and Lipscomb, 1999)
   !
   USE in_out_manager , ONLY : numnam_ice_ref, numnam_ice_cfg, numout, numoni, lwp, lwm, ln_timing ! I/O manager
   USE lib_mpp        , ONLY : ctl_stop, ctl_warn, ctl_nam                              ! MPP library

   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_zdf        ! called by icethd
   PUBLIC   ice_thd_zdf_init   ! called by icestp

   INTEGER ::   nice_zdf       ! Choice of the type of vertical heat diffusion formulation
   !                                 ! associated indices:
   INTEGER, PARAMETER ::   np_BL99 = 1   ! Bitz and Lipscomb (1999)

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_zdf( kl, ll_ice_present )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_zdf  ***
      !!
      !! ** Purpose :   select the appropriate routine for the computation
      !!              of vertical diffusion
      !!-------------------------------------------------------------------
      INTEGER,                     INTENT(in) :: kl
      LOGICAL, DIMENSION(jpi,jpj), INTENT(in) :: ll_ice_present
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_thd_zdf')
      !
      SELECT CASE ( nice_zdf )      ! Choose the vertical heat diffusion solver
         !
         !                          !-------------!
      CASE( np_BL99 )               ! BL99 solver !
         !                          !-------------!
         IF( .NOT.ln_cndflx ) THEN                           ! No conduction flux ==> default option
            CALL ice_thd_zdf_BL99( kl, np_cnd_OFF, ll_ice_present )
         ELSEIF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN   ! Conduction flux as surface boundary condition ==> Met Office default option
            CALL ice_thd_zdf_BL99( kl, np_cnd_ON, ll_ice_present  )
         ELSEIF( ln_cndflx .AND.      ln_cndemulate ) THEN   ! Conduction flux is emulated
            CALL ice_thd_zdf_BL99( kl, np_cnd_EMU, ll_ice_present )
            CALL ice_thd_zdf_BL99( kl, np_cnd_ON, ll_ice_present  )
         ENDIF
         !
      END SELECT
      !
      IF( ln_timing )   CALL timing_stop('ice_thd_zdf')
      !
   END SUBROUTINE ice_thd_zdf


   SUBROUTINE ice_thd_zdf_init
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_zdf_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice thermodynamics
      !!
      !! ** Method  :   Read the namthd_zdf namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_zdf
      !!-------------------------------------------------------------------
      INTEGER  ::   ios, ioptio   ! Local integer
      !!
      NAMELIST/namthd_zdf/ ln_zdf_BL99, ln_cndi_U64, ln_cndi_P07, rn_cnd_s, &
         &                 rn_kappa_i, rn_kappa_s, rn_kappa_smlt, rn_kappa_sdry, ln_zdf_chkcvg
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_zdf)
      READ_NML_CFG(numnam_ice,namthd_zdf)
      IF(lwm) WRITE( numoni, namthd_zdf )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_zdf_init: Ice vertical heat diffusion'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_zdf:'
         WRITE(numout,*) '      Bitz and Lipscomb (1999) formulation                      ln_zdf_BL99   = ', ln_zdf_BL99
         WRITE(numout,*) '      thermal conductivity in the ice (Untersteiner 1964)       ln_cndi_U64   = ', ln_cndi_U64
         WRITE(numout,*) '      thermal conductivity in the ice (Pringle et al 2007)      ln_cndi_P07   = ', ln_cndi_P07
         WRITE(numout,*) '      thermal conductivity in the snow                          rn_cnd_s      = ', rn_cnd_s
         WRITE(numout,*) '      extinction radiation parameter in sea ice                 rn_kappa_i    = ', rn_kappa_i
         WRITE(numout,*) '      extinction radiation parameter in snw      (nn_qtrice=0)  rn_kappa_s    = ', rn_kappa_s
         WRITE(numout,*) '      extinction radiation parameter in melt snw (nn_qtrice=1)  rn_kappa_smlt = ', rn_kappa_smlt
         WRITE(numout,*) '      extinction radiation parameter in dry  snw (nn_qtrice=1)  rn_kappa_sdry = ', rn_kappa_sdry
         WRITE(numout,*) '      check convergence of heat diffusion scheme                ln_zdf_chkcvg = ', ln_zdf_chkcvg
      ENDIF
      !
      rcnd_s = rn_cnd_s ! to be "name compliant" with ice thermal conductivity
      !
      IF ( ( ln_cndi_U64 .AND. ln_cndi_P07 ) .OR. ( .NOT. ln_cndi_U64 .AND. .NOT. ln_cndi_P07 ) ) THEN
         CALL ctl_stop( 'ice_thd_zdf_init: choose 1 and only 1 formulation for thermal conduction (ln_cndi_U64 or ln_cndi_P07)' )
      ENDIF
      !                             !== set the choice of ice vertical thermodynamic formulation ==!
      ioptio = 0
      IF( ln_zdf_BL99 ) THEN
         ioptio = ioptio + 1   ;   nice_zdf = np_BL99
      ENDIF   ! BL99 thermodynamics (linear liquidus + constant thermal properties)
      !      IF( .NOT. ln_zdf_BL99 ) ioptio = ioptio + 1  ! Prevents model from crashing if ln_zdf_BL99 is false
      !!    IF( ln_zdf_XXXX ) THEN   ;   ioptio = ioptio + 1   ;   nice_zdf = np_XXXX   ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_thd_init: one and only one ice thermo option has to be defined ' )
      !
      !$acc update device( ln_zdf_BL99, ln_cndi_U64, ln_cndi_P07, rn_cnd_s, rcnd_s, rn_kappa_i, rn_kappa_s, rn_kappa_smlt, rn_kappa_sdry, ln_zdf_chkcvg )
      !
   END SUBROUTINE ice_thd_zdf_init

   !!======================================================================
END MODULE icethd_zdf
