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
   !! NANUQ 1.0.0, Brodeau (2026)
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
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zsal, zdum, zrhoi_dt, zA, zz1, zz2
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_cor')
      !$acc data present( a_i, h_i, v_i, sss_s, sfx_res, szv_i )

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
               zA = a_i(ji,jj,jl)
               h_i(ji,jj,jl) = MERGE( v_i(ji,jj,jl) / MAX( zA, epsi20 )  ,  0._wp  ,  zA >= epsi20 )
               !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               !   IF( h_i(ji,jj,jl) < rn_himin )  a_ip(ji,jj,jl) = a_ip(ji,jj,jl) * h_i(ji,jj,jl) / rn_himin
               !ENDIF
               a_i(ji,jj,jl) = MERGE( zA * h_i(ji,jj,jl) / rn_himin  ,  zA  ,  h_i(ji,jj,jl) < rn_himin )
            END DO
         END DO
      END DO
      !$acc end parallel loop

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
               zA = a_i(ji,jj,jl)
               a_i(ji,jj,jl) = MERGE( zA * rn_amax / MAX( at_i(ji,jj), epsi20 )  ,  zA  ,  at_i(ji,jj) > rn_amax )
               !IF( at_i(ji,jj) > rn_amax )   a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / at_i(ji,jj)
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
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            zz1 = rn_sinew * sss_s(ji,jj)
            !$acc loop seq
            DO jl = 1, jpl
               zdum = v_i(ji,jj,jl) * r1_nlay_i
               zz2 = zz1 * zdum
               zdum = rn_simin * zdum
               !$acc loop seq
               DO jk=1, nlay_i
                  zsal = szv_i(ji,jj,jk,jl)
                  !szv_i(ji,jj,jk,jl) = MIN( MAX( rn_simin * zdum , zsal ) , rn_sinew * sss_s(ji,jj) * zdum )
                  szv_i(ji,jj,jk,jl) = MIN( MAX( zdum , zsal ) , zz2 )
                  ! no ice-ocean exchanges if kn=0 (for bdy for instance) otherwise conservation diags will fail
                  sfx_res(ji,jj) = sfx_res(ji,jj) - MERGE( 0._wp,  ( szv_i(ji,jj,jk,jl) - zsal ) * zrhoi_dt,  kn==0 )
               END DO
            END DO
         END DO
      END DO
      !$acc end parallel loop

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
