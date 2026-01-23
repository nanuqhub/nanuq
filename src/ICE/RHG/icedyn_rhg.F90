MODULE icedyn_rhg
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg  ***
   !!   Sea-Ice dynamics : master routine for rheology
   !!======================================================================
   !! history :  4.0  !  2018     (C. Rousset)      Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!    ice_dyn_rhg      : computes ice velocities
   !!    ice_dyn_rhg_init : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE par_ice
   USE ice            ! sea-ice: variables
   USE icedyn_rhg_evp ! sea-ice: EVP rheology
   USE icedyn_rhg_bbm ! sea-ice: BBM rheology !#bbm
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg        ! called by icestp.F90
   PUBLIC   ice_dyn_rhg_init   ! called by icestp.F90



   INTEGER ::              nice_rhg   ! choice of the type of rheology
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_rhgEVP = 1   ! EVP rheology
   INTEGER, PARAMETER ::   np_rhgBBM = 2   ! BBM rheology

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icedyn_rhg.F90 14072 2020-12-04 07:48:38Z laurent $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn_rhg  ***
      !!
      !! ** Purpose :   compute ice velocity
      !!
      !! ** Action  : comupte - ice velocity (u_ice, v_ice)
      !!                      - 3 components of the stress tensor (stress1_i, stress2_i, stress12_i)
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      !!--------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_rhg')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn_rhg: sea-ice rheology'
         WRITE(numout,*)'~~~~~~~~~~~'
      ENDIF
      !
      !--------------!
      !== Rheology ==!
      !--------------!
      SELECT CASE( nice_rhg )
         !                                !------------------------!
      CASE( np_rhgEVP )                ! Elasto-Viscous-Plastic !
         !                             !------------------------!
         CALL ice_dyn_rhg_evp( kt, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )
         !
         !                             !-------------------------!
      CASE( np_rhgBBM )                ! Brittle Bingham Maxwell !
         !                             !-------------------------!
         CALL ice_dyn_rhg_bbm( kt, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )
         !
      END SELECT
      !%acc update self ( stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )

      IF( lrst_ice ) THEN
         IF( ln_rhg_EVP )   CALL rhg_evp_rst( 'WRITE', kt )
         IF( ln_rhg_BBM )   CALL rhg_bbm_rst( 'WRITE', kt )
      ENDIF
      !
      ! controls
      IF( sn_cfctl%l_prtctl ) &
         &                 CALL ice_prt3D   ('icedyn_rhg')                                                             ! prints
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      IF( ln_timing    )   CALL timing_stop ('icedyn_rhg')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_rhg


   SUBROUTINE ice_dyn_rhg_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_rhg_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn_rhg namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_rhg
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!  ln_x_MC_test  = .true.          ! Perform only 1 Mohr-Coulomb test, at mid-point between T & F points (implicit smoothing)
      NAMELIST/namdyn_rhg/  ln_idealized, &
         &                  ln_rhg_EVP, rn_creepl, rn_ecc, nn_nevp, rn_relast, nn_rhg_chkcvg,  &  !-- evp
         &                  ln_rhg_BBM, rn_Nref, rn_E0, rn_eta0, rn_P0, rn_kth, nn_nbbm,                 &
         &                  nn_d_adv, ln_adv_d_pra, ln_adv_d_wnx, ln_x_MC_test, rn_crndg,  &  !-- bbm
         &                  rn_dmg_max, rn_C0, nn_alrlx, nn_btrlx, rn_c_ref, rn_l_ref                       !-- bbm
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn_rhg)
      !901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_rhg in reference namelist' )
      READ_NML_CFG(numnam_ice,namdyn_rhg)
      !902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namdyn_rhg in configuration namelist' )
      IF(lwm) WRITE ( numoni, namdyn_rhg )
      !

      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_rhg_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist : namdyn_rhg:'
         IF(ln_idealized) WRITE(numout,*) '   Disregarding Coriolis & SSH terms in the momentum eq.    ln_idealized = ', ln_idealized
         IF( ln_rhg_EVP ) THEN
            WRITE(numout,*) '      rheology EVP (icedyn_rhg_evp) (adaptive, aka aEVP)   ln_rhg_EVP    = ', ln_rhg_EVP
            WRITE(numout,*) '         creep limit                                       rn_creepl     = ', rn_creepl ! also used by vp
            WRITE(numout,*) '         eccentricity of the elliptical yield curve        rn_ecc        = ', rn_ecc    ! also used by vp
            WRITE(numout,*) '         number of iterations for subcycling               nn_nevp       = ', nn_nevp
            WRITE(numout,*) '         ratio of elastic timescale over ice time step     rn_relast     = ', rn_relast
            WRITE(numout,*) '         check convergence of rheology                     nn_rhg_chkcvg = ', nn_rhg_chkcvg
         ENDIF
         IF( ln_rhg_BBM ) THEN
            IF(.NOT.ln_damage) CALL ctl_stop( 'ice_dyn_rhg_init: BBM rheology => set `ln_damage=.true` in `nampar`' )
            WRITE(numout,*) '    rheology BBM (icedyn_rhg_bbm)                          ln_rhg_BBM    = ', ln_rhg_BBM !#bbm
            !IF(ln_MEB) WRITE(numout,*) '         will use the MEB rheology variant rather than pure BBM!'
            WRITE(numout,*) '         max. compressive stress at the ref. scale [Pa]    rn_Nref       = ', rn_Nref  !#bbm
            WRITE(numout,*) '         elasticity of undamaged ice [Pa]                  rn_E0         = ', rn_E0  !#bbm
            WRITE(numout,*) '         viscosity of undamaged ice  [Pa.s]                rn_eta0       = ', rn_eta0  !#bbm
            WRITE(numout,*) '         compression factor "P" at play in "P_max"         rn_P0         = ', rn_P0  !#bbm
            WRITE(numout,*) '         healing constant for damage                       rn_kth        = ', rn_kth  !#bbm
            WRITE(numout,*) '         number of iterations for subcycling               nn_nbbm       = ', nn_nbbm !#bbm
            WRITE(numout,*) '         advection of damage and stresses @T & @F          nn_d_adv  = ', nn_d_adv !#bbm
            IF( nn_d_adv==0 ) WRITE(numout,*) '           => no advection at all!' !#bbm
            IF( nn_d_adv==1 ) WRITE(numout,*) '           => advection of damage only' !#bbm
            IF( nn_d_adv >1 ) WRITE(numout,*) '           => advection of damage + stress tensors' !#bbm
            IF( nn_d_adv==3 ) WRITE(numout,*) '             ==> add "lower-convected" term in tensor advection' !#bbm
            IF( nn_d_adv==4 ) WRITE(numout,*) '             ==> add "upper-convected" term in tensor advection' !#bbm
            IF( nn_d_adv >4 ) CALL ctl_stop( 'ice_dyn_rhg_init: valid choices for `nn_d_adv` span 0 to 4' )
            IF( nn_d_adv >0 ) THEN
               IF( ln_adv_d_pra ) WRITE(numout,*) '         => will use Prather advection scheme'
               IF( ln_adv_d_wnx ) WRITE(numout,*) '         => will use WENO advection scheme'
            ENDIF
            !
            WRITE(numout,*) '         perform Mohr-Coulomb test at mid-point            ln_x_MC_test  = ', ln_x_MC_test
            WRITE(numout,*) '         cross-nudging coeff. for stress tensor            rn_crndg      = ', rn_crndg !#bbm
            !IF(rn_crndg>0._wp) THEN
            !   WRITE(numout,*) '      => boost the CN at the coastline?             ln_boost_CN_coast = ', ln_boost_CN_coast
            !   IF(ln_boost_CN_coast) WRITE(numout,*) '                   ==>          rn_max_CN_coast = ', rn_max_CN_coast
            !   WRITE(numout,*) '      => boost the CN where damage is high?      ln_boost_CN_high_dmg = ', ln_boost_CN_high_dmg
            !   IF(ln_boost_CN_high_dmg) WRITE(numout,*) '                   ==>         rn_max_CN_dmg = ', rn_max_CN_dmg
            !ENDIF
            WRITE(numout,*) '         ceiling value to cap max damage with              rn_dmg_max    = ', rn_dmg_max !#bbm
            WRITE(numout,*) '         compaction paramater (coeff. of exponential)      rn_C0         = ', rn_C0 !#bbm
            WRITE(numout,*) '         `alpha` of viscosity (dep. on `A` and `h`)        nn_alrlx      = ', nn_alrlx !#bbm
            WRITE(numout,*) '         `alpha`-like of Olason/Boutin in `lamdda`         nn_btrlx      = ', nn_btrlx !#bbm
            WRITE(numout,*) '         ice cohesion value at the lab scale               rn_c_ref      = ', rn_c_ref !#bbm
            WRITE(numout,*) '         scaling parameter for cohesion                    rn_l_ref      = ', rn_l_ref  !#bbm
         END IF


         IF( ln_rhg_EVP ) THEN
            IF    ( nn_rhg_chkcvg == 0 ) THEN   ;   WRITE(numout,*) '         no check cvg'
            ELSEIF( nn_rhg_chkcvg == 1 ) THEN   ;   WRITE(numout,*) '         check cvg at the main time step'
            ELSEIF( nn_rhg_chkcvg == 2 ) THEN   ;   WRITE(numout,*) '         check cvg at both main and rheology time steps'
            ENDIF
         ENDIF

         WRITE(numout,*) ''
      ENDIF!IF(lwp)

      !! Min value for `1-d`:
      r_dmd_min = 1._wp - rn_dmg_max
      !WRITE(numout,*) '         `r_dmd_min` set to = ', r_dmd_min, ' on proc #', narea

      !$acc update device( rn_dmg_max, r_dmd_min, rn_C0, nn_alrlx, nn_btrlx, rn_c_ref, rn_l_ref )

      ridlzd = 1._wp
      IF(ln_idealized) THEN
         ridlzd = 0._wp
         IF( lwp ) WRITE(numout,*) ' * Disregarding Coriolis and SSH terms in momentum eq.:',ln_idealized,'=> ridlzd =',INT(ridlzd,1)
      ENDIF

      !                             !== set the choice of ice advection ==!
      ioptio = 0
      IF( ln_rhg_EVP ) THEN
         ioptio = ioptio + 1   ;   nice_rhg = np_rhgEVP
      ENDIF
      IF( ln_rhg_BBM ) THEN
         ioptio = ioptio + 1   ;   nice_rhg = np_rhgBBM
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_dyn_rhg_init: choose one and only one ice rheology' )

      !LOLOrm: IF( ln_rhg_EVP  )   CALL rhg_evp_rst( 'READ' )  !* read or initialize all required files

      IF( ln_rhg_EVP  ) THEN
         CALL ice_dyn_rhg_evp_init() !* allocation of EVP-specific arrays (LB: needs to be done before reading restarts...)
         CALL rhg_evp_rst( 'READ' )  !* read or initialize all required files
         !$acc update device( ln_rhg_EVP, rn_creepl, rn_ecc, nn_nevp, rn_relast )
      END IF

      IF( ln_rhg_BBM  ) THEN
         CALL ice_dyn_rhg_bbm_init() !* allocation of BBM-specific arrays (LB: needs to be done before reading restarts...)
         CALL rhg_bbm_rst( 'READ' )  !* read or initialize all required files
         !$acc update device( ln_rhg_BBM, nn_nbbm, rn_Nref, rn_P0, rn_E0, rn_eta0, rn_kth, nn_d_adv, rn_crndg )
      END IF

   END SUBROUTINE ice_dyn_rhg_init

   !!======================================================================
END MODULE icedyn_rhg
