MODULE icedyn_adv_util
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_util   ***
   !!   sea-ice : advection
   !!======================================================================
   !! History :       !  2008-03  (M. Vancoppenolle) original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE phycst        ! physical constant
   USE dom_oce        ! ocean domain
   USE par_ice, ONLY : jpl, nlay_i, nlay_s, epsi10, epsi20, ln_pnd_LEV, ln_pnd_TOPO
   USE ice            ! sea-ice variables
   USE lib_mpp        ! MPP library
   USE par_ice, ONLY : rconc_min

   IMPLICIT NONE
   PRIVATE

   PUBLIC   Hbig
   PUBLIC   Hsnow
   PUBLIC   icemax3D
   PUBLIC   icemax3D_cnt
   PUBLIC   icemax4D
   PUBLIC   icemax4D_cnt

   INTERFACE Hbig
      MODULE PROCEDURE Hbig_dyn_thd_pnd, Hbig_dyn_thd_s0, Hbig_dyn_thd_sz, Hbig_dyn
   END INTERFACE Hbig

   INTERFACE Hsnow
      MODULE PROCEDURE Hsnow_dyn_thd_pnd, Hsnow_dyn_thd
   END INTERFACE Hsnow

   !!----------------------------------------------------------------------
   !! NANUQ_beta
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE Hbig_dyn_thd_pnd( pdt, phi_max, phs_max, phip_max, psi_max, pes_max, pei_max, &
      &                              pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i, pe_s, pe_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hbig_dyn_thd_pnd  ***
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice or snow
      !!
      !! ** Method  : 1- check whether ice thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by adapting ice concentration
      !!              2- check whether snow thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by sending the excess in the ocean
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdt                                   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(in   ) ::   phi_max, phs_max, phip_max, psi_max   ! max ice thick from surrounding 9-pts
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(in   ) ::   pes_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(in   ) ::   pei_max
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   z1_dt, zhip, zhi, zhs, zsi, zes, zei, zfra
      !!-------------------------------------------------------------------
      !$acc data present( phi_max, phs_max, phip_max, psi_max, pes_max, pei_max, pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i, pe_s, pe_i, wfx_res, hfx_res, sfx_res, hi_max )

      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(3)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl = 1, jpl
               !
               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !
                  !                               ! -- check h_ip -- !
                  ! if h_ip is larger than the surrounding 9 pts => reduce h_ip and increase a_ip
                  IF( ln_pnd_LEV .OR. ln_pnd_TOPO .AND. pv_ip(ji,jj,jl) > 0._wp ) THEN
                     zhip = pv_ip(ji,jj,jl) / MAX( epsi20, pa_ip(ji,jj,jl) )
                     IF( zhip > phip_max(ji,jj,jl) .AND. pa_ip(ji,jj,jl) < 0.15 ) THEN
                        pa_ip(ji,jj,jl) = pv_ip(ji,jj,jl) / phip_max(ji,jj,jl)
                     ENDIF
                  ENDIF
                  !                               ! -- check h_i -- !
                  ! if h_i is larger than the surrounding 9 pts => reduce h_i and increase a_i
                  zhi = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)
                  !LOLO: IF( zhi > phi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                  IF( (phi_max(ji,jj,jl)>epsi10).AND.(zhi > phi_max(ji,jj,jl)).AND.(pa_i(ji,jj,jl) < 0.15) ) THEN !LOLO!
                     pa_i(ji,jj,jl) = pv_i(ji,jj,jl) / MIN( phi_max(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
                  ENDIF
                  !
                  !                               ! -- check h_s -- !
                  ! if h_s is larger than the surrounding 9 pts => put the snow excess in the ocean
                  zhs = pv_s(ji,jj,jl) / pa_i(ji,jj,jl)
                  IF( pv_s(ji,jj,jl) > 0._wp .AND. zhs > phs_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                     zfra = phs_max(ji,jj,jl) / MAX( zhs, epsi20 )
                     !
                     wfx_res(ji,jj) = wfx_res(ji,jj) + ( pv_s(ji,jj,jl) - pa_i(ji,jj,jl) * phs_max(ji,jj,jl) ) * rhos * z1_dt
                     hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                     !
                     pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
                     pv_s(ji,jj,jl)          = pa_i(ji,jj,jl) * phs_max(ji,jj,jl)
                  ENDIF
                  !
                  !                               ! -- check s_i -- !
                  ! if s_i is larger than the surrounding 9 pts => put salt excess in the ocean
                  zsi = psv_i(ji,jj,jl) / pv_i(ji,jj,jl)
                  IF( zsi > psi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                     zfra = psi_max(ji,jj,jl) / zsi
                     sfx_res(ji,jj) = sfx_res(ji,jj) + psv_i(ji,jj,jl) * ( 1._wp - zfra ) * rhoi * z1_dt
                     psv_i(ji,jj,jl) = psv_i(ji,jj,jl) * zfra
                  ENDIF
                  !
               ENDIF

               !                                           ! -- check e_s/v_s -- !
               !$acc loop seq
               DO jk = 1,  nlay_i
                  IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                     ! if e_i/v_i is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zei = pe_i(ji,jj,jk,jl) / pv_i(ji,jj,jl)
                     IF( zei > pei_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                        zfra = pei_max(ji,jj,jk,jl) / zei
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_i(ji,jj,jk,jl) = pe_i(ji,jj,jk,jl) * zfra
                     ENDIF
                  ENDIF
               END DO
               !$acc loop seq
               DO jk = 1,  nlay_s
                  IF ( pv_s(ji,jj,jl) > 0._wp ) THEN
                     ! if e_s/v_s is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zes = pe_s(ji,jj,jk,jl) / pv_s(ji,jj,jl)
                     IF( zes > pes_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                        zfra = pes_max(ji,jj,jk,jl) / zes
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_s(ji,jj,jk,jl) = pe_s(ji,jj,jk,jl) * zfra
                     ENDIF
                  ENDIF
               END DO

            END DO
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data

   END SUBROUTINE Hbig_dyn_thd_pnd


   SUBROUTINE Hbig_dyn_thd_s0( pdt, phi_max, phs_max, psi_max, pes_max, pei_max, pv_i, pv_s, pa_i, psv_i, pe_s, pe_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hbig_dyn_thd_s0  ***
      !!
      !!     ===> For Salt Content array of shape (jpi,jpj,jpl) !!!
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice or snow
      !!
      !! ** Method  : 1- check whether ice thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by adapting ice concentration
      !!              2- check whether snow thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by sending the excess in the ocean
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp)                    ,            INTENT(in   ) ::   pdt                                   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in   ) ::   phi_max, phs_max, psi_max   ! max ice thick from surrounding 9-pts
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(in   ) ::   pes_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(in   ) ::   pei_max
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(inout) ::   pv_i, pv_s, pa_i, psv_i
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   z1_dt, zhi, zhs, zsi, zes, zei, zfra
      !!-------------------------------------------------------------------
      !$acc data present( phi_max, phs_max, psi_max, pes_max, pei_max, pv_i, pv_s, pa_i, psv_i, pe_s, pe_i, wfx_res, hfx_res, sfx_res, hi_max )

      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl

               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !
                  !                               ! -- check h_i -- !
                  ! if h_i is larger than the surrounding 9 pts => reduce h_i and increase a_i
                  zhi = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)
                  !LOLO: IF( zhi > phi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                  IF( (phi_max(ji,jj,jl)>epsi10).AND.(zhi > phi_max(ji,jj,jl)).AND.(pa_i(ji,jj,jl) < 0.15) ) THEN !LOLO!
                     pa_i(ji,jj,jl) = pv_i(ji,jj,jl) / MIN( phi_max(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
                  ENDIF
                  !
                  !                               ! -- check h_s -- !
                  ! if h_s is larger than the surrounding 9 pts => put the snow excess in the ocean
                  zhs = pv_s(ji,jj,jl) / pa_i(ji,jj,jl)
                  IF( pv_s(ji,jj,jl) > 0._wp .AND. zhs > phs_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                     zfra = phs_max(ji,jj,jl) / MAX( zhs, epsi20 )
                     !
                     wfx_res(ji,jj) = wfx_res(ji,jj) + ( pv_s(ji,jj,jl) - pa_i(ji,jj,jl) * phs_max(ji,jj,jl) ) * rhos * z1_dt
                     hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                     !
                     pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
                     pv_s(ji,jj,jl)          = pa_i(ji,jj,jl) * phs_max(ji,jj,jl)
                  ENDIF
                  !
                  !                               ! -- check s_i -- !
                  ! if s_i is larger than the surrounding 9 pts => put salt excess in the ocean
                  zsi = psv_i(ji,jj,jl) / pv_i(ji,jj,jl)
                  IF( zsi > psi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                     zfra = psi_max(ji,jj,jl) / zsi
                     sfx_res(ji,jj) = sfx_res(ji,jj) + psv_i(ji,jj,jl) * ( 1._wp - zfra ) * rhoi * z1_dt
                     psv_i(ji,jj,jl) = psv_i(ji,jj,jl) * zfra
                  ENDIF
                  !
               ENDIF

               !                                           ! -- check e_s/v_s -- !
               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 1,  nlay_i
                     ! if e_i/v_i is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zei = pe_i(ji,jj,jk,jl) / pv_i(ji,jj,jl)
                     !                               ! -- check e_i -- !
                     IF( zei > pei_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                        zfra = pei_max(ji,jj,jk,jl) / zei
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_i(ji,jj,jk,jl) = pe_i(ji,jj,jk,jl) * zfra
                     ENDIF
                  END DO
               ENDIF
               !
               IF ( pv_s(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 1,  nlay_s
                     ! if e_s/v_s is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zes = pe_s(ji,jj,jk,jl) / pv_s(ji,jj,jl)
                     IF( zes > pes_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                        zfra = pes_max(ji,jj,jk,jl) / zes
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_s(ji,jj,jk,jl) = pe_s(ji,jj,jk,jl) * zfra
                     ENDIF
                  END DO
               ENDIF

            END DO !DO jl = 1, jpl
         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop

      !$acc end data

   END SUBROUTINE Hbig_dyn_thd_s0

   SUBROUTINE Hbig_dyn_thd_sz( pdt, phi_max, phs_max, pszi_max, pes_max, pei_max, pv_i, pv_s, pa_i, pszv_i, pe_s, pe_i )
      !!-------------------------------------------------------------------
      !!  => MODEL!
      !!                  ***  ROUTINE Hbig_dyn_thd_sz  ***
      !!
      !!     ===> For Salt Content array of shape (jpi,jpj,nlay_i,jpl) !!!
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice or snow
      !!
      !! ** Method  : 1- check whether ice thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by adapting ice concentration
      !!              2- check whether snow thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by sending the excess in the ocean
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp)                    ,            INTENT(in   ) ::   pdt                                   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in   ) ::   phi_max, phs_max   ! max ice thick from surrounding 9-pts
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(in   ) ::   pszi_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(in   ) ::   pes_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(in   ) ::   pei_max
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in   ) ::   pv_i
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(inout) ::   pv_s, pa_i
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   z1_dt, zA, zV, zhi, zhs, zsi, zes, zei, zfra, zdum
      !!-------------------------------------------------------------------
      !$acc data present( phi_max, phs_max, pszi_max, pes_max, pei_max, pv_i, pv_s, pa_i, pszv_i, pe_s, pe_i, wfx_res, hfx_res, sfx_res, hi_max )

      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl
               
               zV = pv_i(ji,jj,jl)
               zA = pa_i(ji,jj,jl)
               
               IF ( zV > 0._wp .AND. zA > 0._wp ) THEN
                  !
                  !                               ! -- check h_i -- !
                  ! if h_i is larger than the surrounding 9 pts => reduce h_i and increase a_i
                  zhi = zV / zA
                  !LOLO: IF( zhi > phi_max(ji,jj,jl) .AND. zA < 0.15 ) THEN
                  IF( (phi_max(ji,jj,jl)>epsi10).AND.(zhi > phi_max(ji,jj,jl)).AND.(zA < 0.15) ) THEN !LOLO!
                     zA = zV / MIN( phi_max(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
                  ENDIF
                  !
                  !                               ! -- check h_s -- !
                  ! if h_s is larger than the surrounding 9 pts => put the snow excess in the ocean
                  zhs = pv_s(ji,jj,jl) / zA
                  IF( pv_s(ji,jj,jl) > 0._wp .AND. zhs > phs_max(ji,jj,jl) .AND. zA < 0.15 ) THEN
                     zfra = phs_max(ji,jj,jl) / MAX( zhs, epsi20 )
                     !
                     wfx_res(ji,jj) = wfx_res(ji,jj) + ( pv_s(ji,jj,jl) - zA * phs_max(ji,jj,jl) ) * rhos * z1_dt
                     hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                     !
                     pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
                     pv_s(ji,jj,jl)          = zA * phs_max(ji,jj,jl)
                  ENDIF
                  !
               ENDIF
               !
               !                                           ! -- check s_i & e_s/v_s -- !
               IF ( zV > 0._wp .AND. zA > 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 1,  nlay_i
                     zdum = 1._wp / zV
                     !
                     ! if szv_i/v_i is larger than the surrounding 9 pts => put the salt excess in the ocean
                     ! if   e_i/v_i is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zsi = pszv_i(ji,jj,jk,jl) * zdum
                     zei =   pe_i(ji,jj,jk,jl) * zdum
                     !                               ! -- check e_i -- !
                     IF( zei > pei_max(ji,jj,jk,jl) .AND. zA < 0.15 ) THEN
                        zfra = pei_max(ji,jj,jk,jl) / zei
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_i(ji,jj,jk,jl) = pe_i(ji,jj,jk,jl) * zfra
                     ENDIF
                     !                               ! -- check s_i -- !
                     IF( zsi > pszi_max(ji,jj,jk,jl) .AND. zA < 0.15 ) THEN
                        zfra = pszi_max(ji,jj,jk,jl) / zsi
                        sfx_res(ji,jj) = sfx_res(ji,jj) + pszv_i(ji,jj,jk,jl) * ( 1._wp - zfra ) * rhoi * z1_dt
                        pszv_i(ji,jj,jk,jl) = pszv_i(ji,jj,jk,jl) * zfra
                     ENDIF
                     !
                  END DO
               ENDIF
               !
               IF ( pv_s(ji,jj,jl) > 0._wp .AND. zA > 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 1,  nlay_s
                     ! if e_s/v_s is larger than the surrounding 9 pts => put the heat excess in the ocean
                     zes = pe_s(ji,jj,jk,jl) / pv_s(ji,jj,jl)
                     IF( zes > pes_max(ji,jj,jk,jl) .AND. zA < 0.15 ) THEN
                        zfra = pes_max(ji,jj,jk,jl) / zes
                        hfx_res(ji,jj) = hfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                        pe_s(ji,jj,jk,jl) = pe_s(ji,jj,jk,jl) * zfra
                     ENDIF
                  END DO
               ENDIF

               pa_i(ji,jj,jl) = zA
               
            END DO !DO jl = 1, jpl
         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop

      !$acc end data

   END SUBROUTINE Hbig_dyn_thd_sz




   SUBROUTINE Hbig_dyn( pdt, phi_max, phs_max, pv_i, pv_s, pa_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hbig_dyn  ***
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice or snow
      !!
      !! ** Method  : 1- check whether ice thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by adapting ice concentration
      !!              2- check whether snow thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by sending the excess in the ocean
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdt                                   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(in   ) ::   phi_max, phs_max   ! max ice thick from surrounding 9-pts
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_i, pv_s, pa_i
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   z1_dt, zhi, zhs
      !!-------------------------------------------------------------------
      !$acc data present( phi_max, phs_max, pv_i, pv_s, pa_i, hi_max )
      !
      z1_dt = 1._wp / pdt
      !
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl
               !
               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !
                  !                               ! -- check h_i -- !
                  ! if h_i is larger than the surrounding 9 pts => reduce h_i and increase a_i
                  zhi = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)
                  !LOLO:orig:IF( zhi > phi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                  IF( (phi_max(ji,jj,jl)>epsi10).AND.(zhi > phi_max(ji,jj,jl)).AND.(pa_i(ji,jj,jl) < 0.15) ) THEN !LOLO!
                     pa_i(ji,jj,jl) = pv_i(ji,jj,jl) / MIN( phi_max(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
                  ENDIF
                  !
                  !                               ! -- check h_s -- !
                  ! if h_s is larger than the surrounding 9 pts => put the snow excess in the ocean
                  zhs = pv_s(ji,jj,jl) / pa_i(ji,jj,jl)
                  IF( pv_s(ji,jj,jl) > 0._wp .AND. zhs > phs_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
                     !wfx_res(ji,jj) = wfx_res(ji,jj) + ( pv_s(ji,jj,jl) - pa_i(ji,jj,jl) * phs_max(ji,jj,jl) ) * rhos * z1_dt
                     pv_s(ji,jj,jl)          = pa_i(ji,jj,jl) * phs_max(ji,jj,jl)
                  ENDIF
                  !
               ENDIF
               !
            END DO
            !
         END DO
      END DO
      !$acc end parallel loop
      !
      !$acc end data
      !
   END SUBROUTINE Hbig_dyn
   ! zfra

   SUBROUTINE Hsnow_dyn_thd_pnd( pdt, pv_i, pv_s, pa_i, pa_ip, pe_s )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hsnow_dyn_thd_pnd  ***
      !!
      !! ** Purpose : 1- Check snow load after advection
      !!              2- Correct pond concentration to avoid a_ip > a_i
      !!
      !! ** Method :  If snow load makes snow-ice interface to deplet below the ocean surface
      !!              then put the snow excess in the ocean
      !!
      !! ** Notes :   This correction is crucial because of the call to routine icecor afterwards
      !!              which imposes a mini of ice thick. (rn_himin). This imposed mini can artificially
      !!              make the snow very thick (if concentration decreases drastically)
      !!              This behavior has been seen in Ultimate-Macho and supposedly it can also be true for Prather
      !!-------------------------------------------------------------------
      REAL(wp)                    ,            INTENT(in   ) ::   pdt   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in)    ::   pv_i
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(inout) ::   pv_s
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in)    ::   pa_i
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(inout) ::   pa_ip
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      REAL(wp) ::   z1_dt, zvs_excess, zfra
      !!-------------------------------------------------------------------
      !$acc data present( pv_i, pv_s, pa_i, pa_ip, pe_s )
      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(3)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl = 1, jpl
               !
               ! -- check snow load -- !
               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !
                  zvs_excess = MAX( 0._wp, pv_s(ji,jj,jl) - pv_i(ji,jj,jl) * (rho0-rhoi) * r1_rhos )
                  !
                  IF( zvs_excess > 0._wp ) THEN   ! snow-ice interface deplets below the ocean surface
                     ! put snow excess in the ocean
                     zfra = ( pv_s(ji,jj,jl) - zvs_excess ) / MAX( pv_s(ji,jj,jl), epsi20 )
                     wfx_res(ji,jj) = wfx_res(ji,jj) + zvs_excess * rhos * z1_dt
                     hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                     ! correct snow volume and heat content
                     pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
                     pv_s(ji,jj,jl)          = pv_s(ji,jj,jl) - zvs_excess
                  ENDIF
                  !
               ENDIF
               !
               !-- correct pond concentration to avoid a_ip > a_i -- !
               IF( pa_ip(ji,jj,jl) > pa_i(ji,jj,jl) )   pa_ip(ji,jj,jl) = pa_i(ji,jj,jl)
               !
            END DO
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data

   END SUBROUTINE Hsnow_dyn_thd_pnd

   SUBROUTINE Hsnow_dyn_thd( pdt, pv_i, pv_s, pa_i, pe_s )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hsnow_dyn_thd  ***
      !!
      !! ** Purpose : 1- Check snow load after advection
      !!              2- Correct pond concentration to avoid a_ip > a_i
      !!
      !! ** Method :  If snow load makes snow-ice interface to deplet below the ocean surface
      !!              then put the snow excess in the ocean
      !!
      !! ** Notes :   This correction is crucial because of the call to routine icecor afterwards
      !!              which imposes a mini of ice thick. (rn_himin). This imposed mini can artificially
      !!              make the snow very thick (if concentration decreases drastically)
      !!              This behavior has been seen in Ultimate-Macho and supposedly it can also be true for Prather
      !!-------------------------------------------------------------------
      REAL(wp)                    ,            INTENT(in)    :: pdt   ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in)    :: pv_i
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(inout) :: pv_s
      REAL(wp), DIMENSION(jpi,jpj,jpl)  ,      INTENT(in)    :: pa_i !LOLO: not used
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) :: pe_s
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      REAL(wp) ::   z1_dt, zvs_excess, zfra
      !!-------------------------------------------------------------------
      !$acc data present( pv_i, pv_s, pa_i, pe_s )
      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(3)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl = 1, jpl
               !
               ! -- check snow load -- !
               IF ( pv_i(ji,jj,jl) > 0._wp .AND. pa_i(ji,jj,jl) > 0._wp ) THEN
                  !
                  zvs_excess = MAX( 0._wp, pv_s(ji,jj,jl) - pv_i(ji,jj,jl) * (rho0-rhoi) * r1_rhos )
                  !
                  IF( zvs_excess > 0._wp ) THEN   ! snow-ice interface deplets below the ocean surface
                     ! put snow excess in the ocean
                     zfra = ( pv_s(ji,jj,jl) - zvs_excess ) / MAX( pv_s(ji,jj,jl), epsi20 )
                     wfx_res(ji,jj) = wfx_res(ji,jj) + zvs_excess * rhos * z1_dt
                     hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
                     ! correct snow volume and heat content
                     pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
                     pv_s(ji,jj,jl)          = pv_s(ji,jj,jl) - zvs_excess
                  ENDIF
                  !
               ENDIF
               !
            END DO
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data

   END SUBROUTINE Hsnow_dyn_thd




   SUBROUTINE icemax3D( pice , pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax3D ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in ) ::   pice   ! input
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(out) ::   pmax   ! output
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      !!----------------------------------------------------------------------
      ! basic version: get the max of epsi20 + 9 neighbours

      !$acc data present( pice, pmax )

      !$acc parallel loop collapse(3)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            DO jl = 1, jpl
               pmax(ji,jj,jl) = MAX( epsi20, pice(ji-1,jj-1,jl), pice(ji,jj-1,jl), pice(ji+1,jj-1,jl),   &
                  &                          pice(ji-1,jj  ,jl), pice(ji,jj  ,jl), pice(ji+1,jj  ,jl),   &
                  &                          pice(ji-1,jj+1,jl), pice(ji,jj+1,jl), pice(ji+1,jj+1,jl) )
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      !   ==> screws up on GPU !!!!

      !$acc end data
   END SUBROUTINE icemax3D

   SUBROUTINE icemax3D_cnt( pcv, pv, pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax3D_cnt ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in ) ::   pcv    ! content of X
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in ) ::   pv     ! volume of X
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(out) ::   pmax   !
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      REAL(wp) ::   zv
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zc
      !!----------------------------------------------------------------------
      !$acc data present( pcv, pv, pmax ) create( zc )

      ! basic version: get the max of epcvi20 + 9 neighbours

      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl
               zv = MAX( pv(ji,jj,jl), epsi20 )
               zc(ji,jj,jl) = MERGE( pcv(ji,jj,jl) / zv, 0._wp,  zv >= epsi10 )
            END DO
         END DO
      END DO
      !$acc end parallel loop


      !$acc parallel loop collapse(3)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            DO jl = 1, jpl
               pmax(ji,jj,jl) = MAX( epsi20, zc(ji-1,jj-1,jl), zc(ji,jj-1,jl), zc(ji+1,jj-1,jl),   &
                  &                          zc(ji-1,jj  ,jl), zc(ji,jj  ,jl), zc(ji+1,jj  ,jl),   &
                  &                          zc(ji-1,jj+1,jl), zc(ji,jj+1,jl), zc(ji+1,jj+1,jl) )
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      !   ==> screws up on GPU !!!!

      !$acc end data
   END SUBROUTINE icemax3D_cnt


   SUBROUTINE icemax4D( pice , pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax4D ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in ) ::   pice   ! input
      REAL(wp), DIMENSION(:,:,:,:), INTENT(out) ::   pmax   ! output
      !
      INTEGER  ::   jlay, ji, jj, jk, jl   ! dummy loop indices
      !!----------------------------------------------------------------------
      jlay = SIZE( pice , 3 )   ! size of input arrays
      !
      !$acc data present( pice, pmax )

      ! basic version: get the max of epsi20 + 9 neighbours

      !$acc parallel loop collapse(4)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            DO jl = 1, jpl
               DO jk = 1, jlay
                  !
                  pmax(ji,jj,jk,jl) = MAX( epsi20, pice(ji-1,jj-1,jk,jl), pice(ji,jj-1,jk,jl), pice(ji+1,jj-1,jk,jl),   &
                     &                             pice(ji-1,jj  ,jk,jl), pice(ji,jj  ,jk,jl), pice(ji+1,jj  ,jk,jl),   &
                     &                             pice(ji-1,jj+1,jk,jl), pice(ji,jj+1,jk,jl), pice(ji+1,jj+1,jk,jl) )
                  !
               END DO
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      !   ==> screws up on GPU !!!!

      !$acc end data
   END SUBROUTINE icemax4D


   SUBROUTINE icemax4D_cnt( klay, pe, pv, pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax4D_cnt ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      INTEGER,                               INTENT(in ) ::   klay   ! number of layers
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(in ) ::   pe     ! enthalpy of X
      REAL(wp), DIMENSION(jpi,jpj,     jpl), INTENT(in ) ::   pv     ! volume   of X
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(out) ::   pmax   ! output
      !
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) :: zv
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl) ::   ze     ! enthalpy of X / volume of x
      !!----------------------------------------------------------------------
      !$acc data present( pe, pv, pmax ) create( ze )

      !jlay = SIZE( pe , 3 )   ! size of input arrays
      !
      ! basic version: get the max of epsi20 + 9 neighbours

      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl
               !
               zv = MAX( pv(ji,jj,jl), epsi20 )
               !$acc loop seq
               DO jk = 1, klay                  
                  ze(ji,jj,jk,jl) = MERGE( pe(ji,jj,jk,jl) / zv, 0._wp,  zv >= epsi10 )
               END DO
               !
            END DO
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(4)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            DO jl = 1, jpl
               DO jk = 1, klay
                  !
                  pmax(ji,jj,jk,jl) = MAX( epsi20, ze(ji-1,jj-1,jk,jl), ze(ji,jj-1,jk,jl), ze(ji+1,jj-1,jk,jl),   &
                     &                             ze(ji-1,jj  ,jk,jl), ze(ji,jj  ,jk,jl), ze(ji+1,jj  ,jk,jl),   &
                     &                             ze(ji-1,jj+1,jk,jl), ze(ji,jj+1,jk,jl), ze(ji+1,jj+1,jk,jl) )
                  !
               END DO
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      !   ==> screws up on GPU !!!!
      
      !$acc end data
   END SUBROUTINE icemax4D_cnt




   !!======================================================================
END MODULE icedyn_adv_util
