MODULE icecor
   !!======================================================================
   !!                     ***  MODULE  icecor  ***
   !!   sea-ice: Corrections on sea-ice variables at the end of the time step
   !!======================================================================
   !! History :  3.0  !  2006-04  (M. Vancoppenolle) Original code
   !!            3.5  !  2014-06  (C. Rousset)       Complete rewriting/cleaning
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!    ice_cor      : corrections on sea-ice variables
   !!----------------------------------------------------------------------
   USE par_ice
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variable
   USE iceitd  , ONLY : ice_itd_reb
   USE icevar  , ONLY : ice_var_zapsmall
   USE icectl         ! sea-ice: control prints
   USE oss_nnq , ONLY : sss_s

   USE in_out_manager ! I/O manager
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_cor   ! called by icestp.F90

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_cor( kt, kn )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE ice_cor  ***
      !!
      !! ** Purpose :   Computes corrections on sea-ice global variables at
      !!              the end of the dynamics (kn=1) and thermodynamics (kn=2)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      INTEGER, INTENT(in) ::   kn    ! 1 = after dyn ; 2 = after thermo
      !
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zsal, zdum, zrhoi_dt
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_cor')
      !$acc data present( a_i, h_i, v_i, sv_i, sss_s, sfx_res, szv_i )

      zrhoi_dt = rhoi * r1_Dt_ice

      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'ice_cor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icediachk )   CALL ice_cons2D  (0, 'ice_cor',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp .AND. kn == 2 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_cor:  correct sea ice variables if out of bounds '
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !                             !-----------------------------------------------------
      !                             !  ice thickness must exceed himin (for temp. diff.) !
      !                             !-----------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl
               IF( a_i(ji,jj,jl) >= epsi20 ) THEN
                  h_i(ji,jj,jl) = v_i(ji,jj,jl) / a_i(ji,jj,jl)
               ELSE
                  h_i(ji,jj,jl) = 0._wp
               ENDIF
               !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               !   IF( h_i(ji,jj,jl) < rn_himin )  a_ip(ji,jj,jl) = a_ip(ji,jj,jl) * h_i(ji,jj,jl) / rn_himin
               !ENDIF
               IF( h_i(ji,jj,jl) < rn_himin )     a_i (ji,jj,jl) = a_i (ji,jj,jl) * h_i(ji,jj,jl) / rn_himin
            END DO
         END DO
      END DO
      !$acc end parallel loop
      !
      !                             !-----------------------------------------------------
      !                             !  ice concentration should not exceed amax          !
      !                             !-----------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
            !$acc loop seq
            DO jl = 1, jpl
               IF( at_i(ji,jj) > rn_amax )   a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / at_i(ji,jj)
            END DO
            !
         END DO
      END DO
      !$acc end parallel loop

      !                             !-----------------------------------------------------
      !                             !  Rebin categories with thickness out of bounds     !
      !                             !-----------------------------------------------------
      IF( jpl > 1 )   CALL ice_itd_reb( kt )
      !
      !                             !-----------------------------------------------------
      !                             !  salinity must stay in bounds [Simin,Simax]        !
      !                             !-----------------------------------------------------
      IF ( nn_icesal == 2 ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !$acc loop seq
               DO jl = 1, jpl
                  zsal = sv_i(ji,jj,jl)
                  sv_i(ji,jj,jl) = MIN( MAX( rn_simin*v_i(ji,jj,jl) , zsal ) , rn_sinew*sss_s(ji,jj)*v_i(ji,jj,jl)  )
                  IF( kn /= 0 ) & ! no ice-ocean exchanges if kn=0 (for bdy for instance) otherwise conservation diags will fail
                     &   sfx_res(ji,jj) = sfx_res(ji,jj) - ( sv_i(ji,jj,jl) - zsal ) * zrhoi_dt   ! associated salt flux
               END DO
            END DO
         END DO
         !$acc end parallel loop
      ELSEIF ( nn_icesal == 4 ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !$acc loop seq
               DO jl = 1, jpl
                  !$acc loop seq
                  DO jk=1, nlay_i
                     zsal = szv_i(ji,jj,jk,jl)
                     zdum = v_i(ji,jj,jl) * r1_nlay_i
                     szv_i(ji,jj,jk,jl) = MIN( MAX( rn_simin * zdum , zsal ) , rn_sinew * sss_s(ji,jj) * zdum )
                     IF( kn /= 0 ) & ! no ice-ocean exchanges if kn=0 (for bdy for instance) otherwise conservation diags will fail
                        &   sfx_res(ji,jj) = sfx_res(ji,jj) - ( szv_i(ji,jj,jk,jl) - zsal ) * zrhoi_dt   ! associated salt flux
                  END DO
               END DO
            END DO
         END DO
         !$acc end parallel loop
      ENDIF
      !
      IF( kn /= 0 ) THEN   ! no zapsmall if kn=0 (for bdy for instance) because we do not want ice-ocean exchanges (wfx,sfx,hfx)
         !                                                              otherwise conservation diags will fail
         !                          !-----------------------------------------------------
         CALL ice_var_zapsmall      !  Zap small values                                  !
         !                          !-----------------------------------------------------
      ENDIF
      !
      ! controls
      !IF( sn_cfctl%l_prtctl ) &
      !   &                 CALL ice_prt3D   ('ice_cor')                                                             ! prints
      !IF( ln_icectl .AND. kn == 2 ) &
      !   &                 CALL ice_prt     ( kt, iiceprt, jiceprt, 2, ' - Final state - ' )                       ! prints
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'ice_cor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'ice_cor',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      !$acc end data
      IF( ln_timing )   CALL timing_stop ('ice_cor')
      !
   END SUBROUTINE ice_cor


   !!======================================================================
END MODULE icecor
