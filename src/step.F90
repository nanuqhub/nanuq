MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!                    version for standalone surface scheme
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             .   !    .
   !!             .   !    .
   !!   NEMO     3.5  !  2012-03  (S. Alderson)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   stp             : OCE system time-stepping
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean space and time domain variables
   USE daymod           ! calendar                         (day     routine)
   USE icestp,  ONLY : ice_stp  !lulu    ! surface boundary condition: fields
   USE sbcmod,  ONLY : nsbc, sbc, sbc_write ! surface boundary condition       (sbc     routine)
   USE oss_nnq, ONLY : lk_oasis_oce, ln_cpl_oce     ! bottom boundary condition: fields
   USE ossmod , ONLY : oss       ! oceans surface state (oss     routine)
   USE osscpl , ONLY : oss_cpl_snd        ! surface boundary condition: coupled interface
   !
   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE bdy   , ONLY: ln_bdy
   USE bdydta           ! mandatory for sea-ice
   USE stpctl           ! time stepping control            (stp_ctl routine)
   !
   USE in_out_manager   ! I/O manager
   USE prtctl           ! Print control                    (prt_ctl routine)
   USE iom              !
   USE lbclnk           !
   USE timing           ! Timing
   !
   USE par_ice, ONLY : ln_dynADV2D
   !
   USE xios

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nanuqgcm.F90

   !!----------------------------------------------------------------------
   !! time level indices
   !!----------------------------------------------------------------------
   INTEGER, PUBLIC :: Nbb, Nnn, Naa, Nrhs          !! used by nanuq_init
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: step.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of SBC (surface boundary)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Outputs and diagnostics
      !!----------------------------------------------------------------------

      IF( kstp == nit000 )   CALL iom_init( cxios_context ) ! iom_put initialization (must be done after nanuq_init for AGRIF+XIOS+OASIS)
      CALL iom_setkt( kstp - nit000 + 1, cxios_context )   ! tell iom we are at time step kstp
      IF((kstp == nitrst) .AND. lwxios) THEN
         CALL iom_swap(      cw_ocerst_cxt          )
         CALL iom_init_closedef(cw_ocerst_cxt)
         CALL iom_setkt( kstp - nit000 + 1,      cw_ocerst_cxt          )
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )             ! Calendar (day was already called at nit000 in day_init)

      IF(((kstp + nn_fsbc - 1) == nitrst) .AND. lwxios) THEN
         CALL iom_swap(      cw_icerst_cxt          )
         CALL iom_init_closedef(cw_icerst_cxt)
         CALL iom_setkt( kstp - nit000 + 1,      cw_icerst_cxt          )
      ENDIF

      IF( ln_bdy ) CALL bdy_dta( kstp, Nnn )  ! update sea-ice data at open boundaries


      !! Ocean Surface State:
      CALL oss( kstp )       ! obtain an Ocean Surface State (prescribed OR via OASIS coupling to ocean model)
      ! Arrays that have been update by `oss()` (same for standalone and coupled mode) have been put on the GPU, namely:
      ! *** ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m ***


      !IF( .NOT. ln_dynADV2D ) THEN
      !! Surface Boundary Condition for the liquid ocean:
      CALL sbc( kstp )
      ! Arrays that have been update by `sbc()` have been put on the GPU, namely:
      ! *** emp, qsr, qns, qns_oce, qsr_oce, fatm_prcp, fatm_snow, wndm, utau, vtau, taum, rhoa ***
      !!    => not `sfx` & `fmmflx` because they are given a value in ICE model...
      !!    => not `theta_zu` & `q_zu` I guess..
      !ENDIF

      CALL ice_stp( kstp, nsbc )  ! Sea-ice model

      CALL sbc_write( kstp )      ! Outputs and control print for fluxes that serve as SBCs to the liquid ocean component

      CALL dia_wri( kstp,      Nnn )        ! Outputs

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CALL stp_ctl( kstp )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! File manipulation at the end of the first time step
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nit000   ) THEN
         CALL iom_close( numror )                          ! close input  ocean restart file
         IF( lrxios )     CALL iom_context_finalize(      cr_ocerst_cxt      )
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( lk_oasis_oce .AND. nstop == 0 ) THEN
         CALL oss_cpl_snd( kstp )   ! coupled mode: send the SBC for liquid ocean to the ocean component !
      END IF
      IF( kstp == nitrst ) THEN
         IF(.NOT.lwxios) THEN
            CALL iom_close( numrow )
         ELSE
            CALL iom_context_finalize( cw_ocerst_cxt )
            iom_file(numrow)%nfid       = 0
            numrow = 0
         ENDIF
      ENDIF
      IF( kstp == nitend .OR. nstop > 0 ) THEN
         CALL iom_context_finalize( cxios_context ) ! needed for XIOS+AGRIF
      ENDIF
      !
      IF( ln_timing .AND.  kstp == nit000  )   CALL timing_reset
      !
   END SUBROUTINE stp

   !!======================================================================
END MODULE step
