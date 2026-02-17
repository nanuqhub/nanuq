MODULE icevar
   !!======================================================================
   !!                       ***  MODULE icevar ***
   !!   sea-ice:  series of functions to transform or compute ice variables
   !!======================================================================
   !! History :   -   !  2006-01  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!
   !!                 There are three sets of variables
   !!                 VGLO : global variables of the model
   !!                        - v_i (jpi,jpj,jpl)
   !!                        - v_s (jpi,jpj,jpl)
   !!                        - a_i (jpi,jpj,jpl)
   !!                        - t_s (jpi,jpj,jpl)
   !!                        - e_i (jpi,jpj,nlay_i,jpl)
   !!                        - e_s (jpi,jpj,nlay_s,jpl)
   !!                        - sv_i(jpi,jpj,jpl)
   !!                        - oa_i(jpi,jpj,jpl)
   !!                 VEQV : equivalent variables sometimes used in the model
   !!                        - h_i(jpi,jpj,jpl)
   !!                        - h_s(jpi,jpj,jpl)
   !!                        - t_i(jpi,jpj,nlay_i,jpl)
   !!                        ...
   !!                 VAGG : aggregate variables, averaged/summed over all
   !!                        thickness categories
   !!                        - vt_i(jpi,jpj)
   !!                        - vt_s(jpi,jpj)
   !!                        - at_i(jpi,jpj)
   !!                        - st_i(jpi,jpj)
   !!                        - et_s(jpi,jpj)  total snow heat content
   !!                        - et_i(jpi,jpj)  total ice thermal content
   !!                        - sm_i(jpi,jpj)  mean ice salinity
   !!                        - tm_i(jpi,jpj)  mean ice temperature
   !!                        - tm_s(jpi,jpj)  mean snw temperature
   !!----------------------------------------------------------------------
   !!   ice_var_agg       : integrate variables over layers and categories
   !!   ice_var_glo2eqv   : transform from VGLO to VEQV
   !!   ice_var_salprof   : salinity profile in the ice
   !!   ice_var_zapsmall  : remove very small area and volume
   !!   ice_var_zapneg    : remove negative ice fields
   !!   ice_var_roundoff  : remove negative values arising from roundoff erros
   !!   ice_var_brine     : brine volume
   !!   ice_var_enthalpy  : compute ice and snow enthalpies from temperature
   !!   ice_var_sshdyn    : compute equivalent ssh in lead
   !!   ice_var_itd       : convert N-cat to M-cat
   !!   ice_var_snwfra    : fraction of ice covered by snow
   !!   ice_var_snwblow   : distribute snow fall between ice and ocean
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constants (ocean directory)
   USE ice            ! sea-ice: variables
   USE oss_nnq , ONLY : ln_ice_embd, sss_s, sst_s
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   USE remap_classic, ONLY: rmpT2F, do_rmpT2U, do_rmpT2V, rmpT2U, rmpT2V, do_rmpT2F
   USE remap_weno,    ONLY: rmpT2F_A_h_wn5s

   USE lbclnk

   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_var_agg
   PUBLIC   ice_var_agg_gpu
   PUBLIC   ice_var_agg_adv2d_gpu

   PUBLIC   ice_var_glo2eqv
   PUBLIC   ice_var_glo2eqv_gpu

   PUBLIC   ice_var_salprof
   PUBLIC   ice_var_salprof_gpu

   PUBLIC   ice_var_itd_1cMc_2d

   PUBLIC   ice_var_zapsmall
   PUBLIC   ice_var_zapsmall_dyn
   PUBLIC   ice_var_zapneg
   PUBLIC   ice_var_roundoff
   PUBLIC   ice_var_brine
   !PUBLIC   ice_var_enthalpy ! ==> inlined where needed
   PUBLIC   ice_var_vremap
   !PUBLIC   snw_ent
   PUBLIC   ice_var_sshdyn
   PUBLIC   ice_var_itd
   PUBLIC   ice_var_snwfra
   PUBLIC   ice_var_snwfra_sclr
   PUBLIC   ice_var_snwblow
   PUBLIC   ice_var_hpiling
   PUBLIC   ice_var_cap_at


   PUBLIC   l_is_it_a_nan
   PUBLIC   l_is_it_a_inf
   PUBLIC   test4inf
   PUBLIC   test4nan

   INTERFACE ice_var_zapneg
      MODULE PROCEDURE ice_var_zapneg_dyn_thd_pnd, ice_var_zapneg_dyn_thd, ice_var_zapneg_dyn
   END INTERFACE ice_var_zapneg

   INTERFACE ice_var_roundoff
      MODULE PROCEDURE ice_var_roundoff_dyn_thd_pnd, ice_var_roundoff_dyn_thd
   END INTERFACE ice_var_roundoff

   INTERFACE ice_var_itd
      MODULE PROCEDURE ice_var_itd_1c1c, ice_var_itd_Nc1c, ice_var_itd_1cMc, ice_var_itd_NcMc
   END INTERFACE ice_var_itd

   INTERFACE ice_var_snwfra
      MODULE PROCEDURE ice_var_snwfra_1d, ice_var_snwfra_2d, ice_var_snwfra_3d
   END INTERFACE ice_var_snwfra

   INTERFACE ice_var_snwblow
      MODULE PROCEDURE ice_var_snwblow_1d, ice_var_snwblow_2d
   END INTERFACE ice_var_snwblow

   INTERFACE test4inf
      MODULE PROCEDURE test4inf_2d, test4inf_3d
   END INTERFACE test4inf


   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

# include "icevar_agg.h90"


# include "icevar_agg_gpu.h90"


   SUBROUTINE ice_var_agg_adv2d_gpu( )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_agg_adv2d_gpu  ***
      !!
      !! ** Purpose :   aggregates ice-thickness-category variables to
      !!                all-ice variables, i.e. it turns VGLO into VAGG
      !!
      !!   => only used when `ln_dynADV2D==T` with `jpl==1` !!!
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   z1_at_i
      !!-------------------------------------------------------------------
      !$acc data present( af_i, a_i, at_i, ato_i, hm_i, hm_i_f, hm_s, kmsk_ice_f, kmsk_ice_t, v_i, v_s, vt_i, vt_s )
      !    au_i, av_i,  kmsk_ice_u, kmsk_ice_v, umask, vmask

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            vt_i(ji,jj)  =  v_i(ji,jj,1)
            vt_s(ji,jj)  =  v_s(ji,jj,1)
            st_i(ji,jj)  = sv_i(ji,jj,1)
            at_i(ji,jj)  =  a_i(ji,jj,1)
            !
            ato_i(ji,jj) = 1._wp - at_i(ji,jj)         ! open water fraction
            !
            z1_at_i = 0._wp
            IF( at_i(ji,jj) > epsi20 ) z1_at_i = 1._wp / at_i(ji,jj)
            hm_i(ji,jj) = vt_i(ji,jj) * z1_at_i
            hm_s(ji,jj) = vt_s(ji,jj) * z1_at_i
            !
         END DO
      END DO
      !$acc end parallel loop

# if ! defined _OPENACC
      CALL lbc_lnk( 'ice_var_agg_adv2d_gpu',  at_i,'T',1._wp, hm_i,'T',1._wp ) !LOLOfixme: `at_i` or/and `hm_i` LBC_LNKed!
      !                                                                        ! => shows up when `ln_use_weno_rmp` !
      !                                                                        ! => scary, so doing it here...
# endif

      !! Ice concentration and thickness @F,U,V:
      IF( ln_use_weno_rmp ) THEN
         CALL rmpT2F_A_h_wn5s( at_i, hm_i, af_i, hm_i_f )
      ELSE
         CALL do_rmpT2F( at_i, af_i,    lconserv=.TRUE. )
         CALL do_rmpT2F( hm_i, hm_i_f,  lconserv=.TRUE. )
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               af_i(ji,jj)   = MIN( MAX(   af_i(ji,jj), 0._wp ), rn_amax )
               hm_i_f(ji,jj) =      MAX( hm_i_f(ji,jj), 0._wp )
            END DO
         END DO
         !$acc end parallel loop
      ENDIF
      !
      !CALL do_rmpT2U( at_i, au_i,  lconserv=.TRUE. )
      !CALL do_rmpT2V( at_i, av_i,  lconserv=.TRUE. )
!!$acc parallel loop collapse(2)
      !DO jj=Njs0-nn_hls, Nje0+nn_hls
      !   DO ji=Nis0-nn_hls, Nie0+nn_hls
      !      au_i(ji,jj) = MIN( MAX( au_i(ji,jj) , 0._wp ) , rn_amax )
      !      av_i(ji,jj) = MIN( MAX( av_i(ji,jj) , 0._wp ) , rn_amax )
      !   END DO
      !END DO
!!$acc end parallel loop

# if !defined _OPENACC
      !CALL lbc_lnk( 'ice_var_agg_adv2d_gpu',  af_i,'F',1._wp,  hm_i_f,'F',1._wp,  au_i,'U',1._wp,  av_i,'V',1._wp )
      CALL lbc_lnk( 'ice_var_agg_adv2d_gpu',  af_i,'F',1._wp,  hm_i_f,'F',1._wp )
# endif

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls

            ! For diagnostics:
            kmsk_ice_t(ji,jj) = xmskt(ji,jj)
            kmsk_ice_f(ji,jj) = xmskf(ji,jj)
            IF( at_i(ji,jj) < rAmin_fld ) kmsk_ice_t(ji,jj) = 0
            IF( af_i(ji,jj) < rAmin_fld ) kmsk_ice_f(ji,jj) = 0

            ! For rheology and diags:
            !kmsk_ice_u(ji,jj) = INT(umask(ji,jj,1))
            !kmsk_ice_v(ji,jj) = INT(vmask(ji,jj,1))
            !IF( au_i(ji,jj) < rAmin_fld ) kmsk_ice_u(ji,jj) = 0
            !IF( av_i(ji,jj) < rAmin_fld ) kmsk_ice_v(ji,jj) = 0

         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      !
   END SUBROUTINE ice_var_agg_adv2d_gpu


# include "icevar_glo2eqv.h90"

# include "icevar_glo2eqv_gpu.h90"

# include "icevar_salprof.h90"

# include "icevar_salprof_gpu.h90"


   SUBROUTINE ice_var_zapsmall
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapsmall ***
      !!
      !! ** Purpose :   Remove too small sea ice areas and correct fluxes
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   zsmall
      !!-------------------------------------------------------------------
      !$acc data present( a_i, at_i, e_i, e_s, hfx_res, h_i, h_s, oa_i, sfx_res, sv_i, t_i, t_s, t_su, v_i, v_s, vt_i, wfx_pnd, wfx_res, sz_i, szv_i )
      !
      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl
               IF( a_i(ji,jj,jl) > epsi10 ) THEN
                  h_i(ji,jj,jl) = v_i(ji,jj,jl) / a_i(ji,jj,jl)
               ELSE
                  h_i(ji,jj,jl) = 0._wp
               ENDIF
            END DO
         END DO
      END DO
      !$acc end parallel loop
      !
      !-----------------------------------------------------------------
      ! Zap ice volume, add salt to ocean
      !-----------------------------------------------------------------
      IF( nn_icesal == 4 ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !$acc loop seq
               DO jl = 1, jpl
                  !$acc loop seq
                  DO jk=1, nlay_i
                     zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
                     IF( zsmall < epsi10 ) THEN
                        ! update exchanges with ocean
                        sfx_res(ji,jj)  = sfx_res(ji,jj) + szv_i(ji,jj,jk,jl) * rhoi * r1_Dt_ice
                        szv_i(ji,jj,jk,jl) = 0._wp
                        sz_i (ji,jj,jk,jl) = rn_simin
                     ENDIF
                  END DO
               END DO
               !
            END DO
         ENDDO
         !$acc end parallel loop
      ELSE
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !$acc loop seq
               DO jl = 1, jpl
                  zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
                  IF( zsmall < epsi10 ) THEN
                     ! update exchanges with ocean
                     sfx_res(ji,jj)  = sfx_res(ji,jj) + sv_i(ji,jj,jl)   * rhoi * r1_Dt_ice
                     sv_i(ji,jj,jl) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl
               zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
               IF( zsmall < epsi10 ) THEN
                  !$acc loop seq
                  DO jk=1, nlay_i
                     ! update exchanges with ocean
                     hfx_res(ji,jj)   = hfx_res(ji,jj) - e_i(ji,jj,jk,jl) * r1_Dt_ice ! W.m-2 <0
                     e_i(ji,jj,jk,jl) = 0._wp
                     t_i(ji,jj,jk,jl) = rt0
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_s
                     ! update exchanges with ocean
                     hfx_res(ji,jj)   = hfx_res(ji,jj) - e_s(ji,jj,jk,jl) * r1_Dt_ice ! W.m-2 <0
                     e_s(ji,jj,jk,jl) = 0._wp
                     t_s(ji,jj,jk,jl) = rt0
                  END DO
               ENDIF
            END DO !DO jl = 1, jpl
            !
         END DO
      ENDDO
      !$acc end parallel loop
      !
      !-----------------------------------------------------------------
      ! zap ice and snow volume, add water to ocean
      !-----------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl               !
               zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
               IF( zsmall < epsi10 ) THEN
                  ! update exchanges with ocean
                  wfx_res(ji,jj)  = wfx_res(ji,jj) + v_i (ji,jj,jl)   * rhoi * r1_Dt_ice
                  wfx_res(ji,jj)  = wfx_res(ji,jj) + v_s (ji,jj,jl)   * rhos * r1_Dt_ice
                  !wfx_res(ji,jj)  = wfx_res(ji,jj) + ( v_ip(ji,jj,jl)+v_il(ji,jj,jl) ) * rhow * r1_Dt_ice
                  !
                  a_i  (ji,jj,jl) = 0._wp
                  v_i  (ji,jj,jl) = 0._wp
                  v_s  (ji,jj,jl) = 0._wp
                  t_su (ji,jj,jl) = sst_s(ji,jj) + rt0
                  oa_i (ji,jj,jl) = 0._wp
                  !
                  h_i (ji,jj,jl) = 0._wp
                  h_s (ji,jj,jl) = 0._wp
                  !
                  !a_ip (ji,jj,jl) = 0._wp
                  !v_ip (ji,jj,jl) = 0._wp
                  !v_il (ji,jj,jl) = 0._wp
                  !h_ip (ji,jj,jl) = 0._wp
                  !h_il (ji,jj,jl) = 0._wp
               ENDIF
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! to be sure that at_i is the sum of a_i(jl)
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            at_i(ji,jj) = 0._wp
            vt_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
               vt_i(ji,jj) = vt_i(ji,jj) + v_i(ji,jj,jl)
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! open water = 1 if at_i=0
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF( at_i(ji,jj) == 0._wp )   ato_i(ji,jj) = 1._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE ice_var_zapsmall

   SUBROUTINE ice_var_zapsmall_dyn
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapsmall_dyn ***
      !!
      !! ** Purpose :   Remove too small sea ice areas and correct fluxes
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   zsmall
      !!-------------------------------------------------------------------
      !$acc data present( a_i, at_i, h_i, h_s, oa_i, v_i, v_s, vt_i )
      !
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl
               IF( a_i(ji,jj,jl) > epsi10 ) THEN
                  h_i(ji,jj,jl) = v_i(ji,jj,jl) / a_i(ji,jj,jl)
               ELSE
                  h_i(ji,jj,jl) = 0._wp
               ENDIF
            END DO
         END DO
      END DO
      !$acc end parallel loop

      !-----------------------------------------------------------------
      ! zap ice and snow volume, add water to ocean
      !-----------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !$acc loop seq
            DO jl = 1, jpl               !
               zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
               IF( zsmall < epsi10 ) THEN
                  a_i (ji,jj,jl) = 0._wp
                  v_i (ji,jj,jl) = 0._wp
                  v_s (ji,jj,jl) = 0._wp
                  oa_i(ji,jj,jl) = 0._wp
                  h_i (ji,jj,jl) = 0._wp
                  h_s (ji,jj,jl) = 0._wp
               ENDIF
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! to be sure that at_i is the sum of a_i(jl)
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            at_i(ji,jj) = 0._wp
            vt_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
               vt_i(ji,jj) = vt_i(ji,jj) + v_i(ji,jj,jl)
            END DO
         END DO
      END DO
      !$acc end parallel loop

      ! open water = 1 if at_i=0
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( at_i(ji,jj) == 0._wp )   ato_i(ji,jj) = 1._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE ice_var_zapsmall_dyn


   SUBROUTINE ice_var_zapneg_dyn_thd_pnd( pdt, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i ) !, pszv_i )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapneg_dyn_thd_pnd ***
      !!
      !! ** Purpose :   Remove negative sea ice fields and correct fluxes
      !!-------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_il      ! melt pond lid volume
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s  ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i  ! ice heat content
      !REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i     ! ice salt content
      !
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   z1_dt
      REAL(wp), DIMENSION(jpi,jpj) ::   zwfx_res, zhfx_res, zsfx_res ! needed since loop is not (0,0,0,0)
      !!-------------------------------------------------------------------
      !%acc data present( pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, hfx_res, sfx_res, wfx_res, wfx_pnd )
      !! --> , pszv_i
      !
      !%acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            zwfx_res(ji,jj) = 0._wp
            zhfx_res(ji,jj) = 0._wp
            zsfx_res(ji,jj) = 0._wp
         END DO
      END DO
      !%acc end parallel loop

      z1_dt = 1._wp / pdt
      !
      ! make sure a_i=0 where v_i<=0
      WHERE( pv_i(:,:,:) <= 0._wp )   pa_i(:,:,:) = 0._wp

      !--------------------------------------
      ! zap ice salt and send it to the ocean
      !--------------------------------------
      IF( nn_icesal == 4 ) THEN
         PRINT *, 'STOP! ice_var_zapneg_dyn_thd_pnd@icevar.F90 ==> re-add the `nn_icesal == 4` !!!'
         STOP
         !DO jl = 1, jpl
         !   DO jj=Njs0-1, Nje0+1
         !      DO ji=Nis0-1, Nie0+1
         !         DO jk=1, nlay_i
         !            IF( pszv_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
         !               zsfx_res(ji,jj)     = zsfx_res(ji,jj) + pszv_i(ji,jj,jk,jl) * rhoi * z1_dt
         !               pszv_i(ji,jj,jk,jl) = 0._wp
         !            ENDIF
         !         END DO
         !      END DO
         !   END DO
         !ENDDO
      ELSE
         DO jl = 1, jpl
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  IF( psv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                     zsfx_res(ji,jj)    = zsfx_res(ji,jj) + psv_i(ji,jj,jl) * rhoi * z1_dt
                     psv_i   (ji,jj,jl) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !
      !----------------------------------------
      ! zap ice energy and send it to the ocean
      !----------------------------------------
      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               DO jk=1, nlay_i
                  IF( pe_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                     zhfx_res(ji,jj)   = zhfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * z1_dt ! W.m-2 >0
                     pe_i(ji,jj,jk,jl) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
         !
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               DO jk=1, nlay_s
                  IF( pe_s(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_s(ji,jj,jl) <= 0._wp ) THEN
                     zhfx_res(ji,jj)   = zhfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * z1_dt ! W.m-2 <0
                     pe_s(ji,jj,jk,jl) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
      ENDDO
      !
      !--------------------------------------------
      ! zap ice and snow volume, add water to ocean
      !--------------------------------------------
      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               IF( pv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_i (ji,jj,jl) * rhoi * z1_dt
                  pv_i    (ji,jj,jl) = 0._wp
               ENDIF
               IF( pv_s(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_s (ji,jj,jl) * rhos * z1_dt
                  pv_s    (ji,jj,jl) = 0._wp
               ENDIF
               IF( pv_ip(ji,jj,jl) < 0._wp .OR. pv_il(ji,jj,jl) < 0._wp .OR. pa_ip(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_il(ji,jj,jl) * rhow * z1_dt
                  pv_il   (ji,jj,jl) = 0._wp
               ENDIF
               IF( pv_ip(ji,jj,jl) < 0._wp .OR. pa_ip(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_ip(ji,jj,jl) * rhow * z1_dt
                  pv_ip   (ji,jj,jl) = 0._wp
               ENDIF
            END DO
         END DO
      END DO
      !%acc end parallel loop
      !
      ! record residual fluxes
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            wfx_res(ji,jj) = wfx_res(ji,jj) + zwfx_res(ji,jj)
            hfx_res(ji,jj) = hfx_res(ji,jj) + zhfx_res(ji,jj)
            sfx_res(ji,jj) = sfx_res(ji,jj) + zsfx_res(ji,jj)
         END DO
      END DO
      !
      WHERE( poa_i (:,:,:) < 0._wp )   poa_i (:,:,:) = 0._wp
      WHERE( pa_i  (:,:,:) < 0._wp )   pa_i  (:,:,:) = 0._wp
      WHERE( pa_ip (:,:,:) < 0._wp )   pa_ip (:,:,:) = 0._wp
      !
      !%acc end data
   END SUBROUTINE ice_var_zapneg_dyn_thd_pnd

   SUBROUTINE ice_var_zapneg_dyn_thd( pdt, pv_i, pv_s, psv_i, poa_i, pa_i, pe_s, pe_i,  pszv_i )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapneg_dyn_thd ***
      !!
      !! ** Purpose :   Remove negative sea ice fields and correct fluxes
      !!-------------------------------------------------------------------
      REAL(wp)                               , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i       ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), OPTIONAL, INTENT(inout) ::   pszv_i     ! ice salt content
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   z1_dt
      REAL(wp) ::   zwfx_res, zhfx_res, zsfx_res ! needed since loop is not (0,0,0,0)
      LOGICAL  ::   l_do_szv_i
      !!-------------------------------------------------------------------
      !$acc data present( pv_i, pv_s, psv_i, poa_i, pa_i, pe_s, pe_i, hfx_res, sfx_res, wfx_res )
      !
      l_do_szv_i = PRESENT( pszv_i )

      z1_dt = 1._wp / pdt

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            zwfx_res = 0._wp
            zhfx_res = 0._wp
            zsfx_res = 0._wp

            !$acc loop seq
            DO jl=1, jpl

               ! make sure a_i=0 where v_i<=0
               IF( pv_i(ji,jj,jl) <= 0._wp )   pa_i(ji,jj,jl) = 0._wp

               !--------------------------------------
               ! zap ice salt and send it to the ocean
               !--------------------------------------
               !IF( nn_icesal == 4 ) THEN
               IF( l_do_szv_i ) THEN
                  !!  ==> it obviously implies that `nn_icesal == 4` !
                  !$acc loop seq
                  DO jk=1, nlay_i
                     IF( pszv_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                        zsfx_res     = zsfx_res + pszv_i(ji,jj,jk,jl) * rhoi * z1_dt
                        pszv_i(ji,jj,jk,jl) = 0._wp
                     ENDIF
                  END DO
               ELSE
                  IF( psv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                     zsfx_res    = zsfx_res + psv_i(ji,jj,jl) * rhoi * z1_dt
                     psv_i   (ji,jj,jl) = 0._wp
                  ENDIF
               ENDIF

               !----------------------------------------
               ! zap ice energy and send it to the ocean
               !----------------------------------------
               !$acc loop seq
               DO jk=1, nlay_i
                  IF( pe_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                     zhfx_res   = zhfx_res - pe_i(ji,jj,jk,jl) * z1_dt ! W.m-2 >0
                     pe_i(ji,jj,jk,jl) = 0._wp
                  ENDIF
               END DO
               !$acc loop seq
               DO jk=1, nlay_s
                  IF( pe_s(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_s(ji,jj,jl) <= 0._wp ) THEN
                     zhfx_res   = zhfx_res - pe_s(ji,jj,jk,jl) * z1_dt ! W.m-2 <0
                     pe_s(ji,jj,jk,jl) = 0._wp
                  ENDIF
               END DO

               !--------------------------------------------
               ! zap ice and snow volume, add water to ocean
               !--------------------------------------------
               IF( pv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res    = zwfx_res + pv_i (ji,jj,jl) * rhoi * z1_dt
                  pv_i    (ji,jj,jl) = 0._wp
               ENDIF
               IF( pv_s(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  zwfx_res    = zwfx_res + pv_s (ji,jj,jl) * rhos * z1_dt
                  pv_s    (ji,jj,jl) = 0._wp
               ENDIF

               IF( poa_i (ji,jj,jl) < 0._wp )   poa_i (ji,jj,jl) = 0._wp
               IF( pa_i  (ji,jj,jl) < 0._wp )   pa_i  (ji,jj,jl) = 0._wp

            END DO !DO jl=1, jpl

            ! record residual fluxes
            wfx_res(ji,jj) = wfx_res(ji,jj) + zwfx_res
            hfx_res(ji,jj) = hfx_res(ji,jj) + zhfx_res
            sfx_res(ji,jj) = sfx_res(ji,jj) + zsfx_res

         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE ice_var_zapneg_dyn_thd


   SUBROUTINE ice_var_zapneg_dyn( pdt, pv_i, pv_s, pa_i )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapneg_dyn ***
      !!
      !! ** Purpose :   Remove negative sea ice fields and correct fluxes
      !!-------------------------------------------------------------------
      REAL(wp)                               , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i       ! ice concentration
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl  ! dummy loop indices
      !!-------------------------------------------------------------------
      IF( ln_timing )  CALL timing_start('ice_var_zapneg_dyn')
      !$acc data present( pv_i, pv_s, pa_i )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            !$acc loop seq
            DO jl=1, jpl

               ! make sure a_i=0 where v_i<=0
               IF( pv_i(ji,jj,jl) <= 0._wp )   pa_i(ji,jj,jl) = 0._wp

               !--------------------------------------------
               ! zap ice and snow volume, add water to ocean
               !--------------------------------------------
               IF( pv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  pv_i(ji,jj,jl) = 0._wp
               ENDIF
               IF( pv_s(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
                  pv_s(ji,jj,jl) = 0._wp
               ENDIF

               IF( pa_i  (ji,jj,jl) < 0._wp )   pa_i  (ji,jj,jl) = 0._wp

            END DO !DO jl=1, jpl

         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )  CALL timing_stop('ice_var_zapneg_dyn')
      !
   END SUBROUTINE ice_var_zapneg_dyn

   SUBROUTINE ice_var_roundoff_dyn_thd_pnd( pa_i, pv_i, pv_s, psv_i, poa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i, ll_ice_present )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_roundoff_dyn_thd_pnd ***
      !!
      !! ** Purpose :   Remove negative sea ice values arising from roundoff errors
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_il      ! melt pond lid volume
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i       ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i     ! ice salt content
      LOGICAL,  DIMENSION(jpi,jpj),            INTENT(in)    :: ll_ice_present
      !!-------------------------------------------------------------------
      INTEGER :: ji, jj, jk, jl
      !!-------------------------------------------------------------------
      !$acc data present( pa_i, pv_i, pv_s, psv_i, poa_i, pe_s, pe_i, pszv_i, pa_ip, pv_ip, pv_il, ll_ice_present )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF (ll_ice_present(ji,jj)) THEN
               !
               !$acc loop seq
               DO jl=1, jpl
                  pa_i (ji,jj,jl) = MAX( pa_i(ji,jj,jl), 0._wp)
                  pv_i (ji,jj,jl) = MAX( pv_i(ji,jj,jl), 0._wp)
                  pv_s (ji,jj,jl) = MAX( pv_s(ji,jj,jl), 0._wp)
                  poa_i(ji,jj,jl) = MAX(poa_i(ji,jj,jl), 0._wp)
                  IF( nn_icesal /= 4 ) psv_i(ji,jj,jl) = MAX(psv_i(ji,jj,jl), 0._wp)
                  !$acc loop seq
                  DO jk=1, nlay_i
                     pe_i  (ji,jj,jk,jl) = MAX(  pe_i(ji,jj,jk,jl), 0._wp)
                     IF( nn_icesal == 4 ) pszv_i(ji,jj,jk,jl) = MAX(pszv_i(ji,jj,jk,jl), 0._wp)
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_s
                     pe_s(ji,jj,jk,jl) = MAX(pe_s(ji,jj,jk,jl), 0._wp)
                  END DO
               END DO
               !
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               IF (ll_ice_present(ji,jj)) THEN
                  !$acc loop seq
                  DO jl=1, jpl
                     pa_ip(ji,jj,jl) = MAX(pa_ip(ji,jj,jl), 0._wp)
                     pv_ip(ji,jj,jl) = MAX(pv_ip(ji,jj,jl), 0._wp)
                  END DO
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
         IF( ln_pnd_lids ) THEN
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  IF (ll_ice_present(ji,jj)) THEN
                     !$acc loop seq
                     DO jl=1, jpl
                        pv_il(ji,jj,jl) = MAX(pv_il(ji,jj,jl), 0._wp)
                     END DO
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop
         ENDIF
      ENDIF

      !$acc end data
   END SUBROUTINE ice_var_roundoff_dyn_thd_pnd

   SUBROUTINE ice_var_roundoff_dyn_thd( pa_i, pv_i, pv_s, psv_i, poa_i, pe_s, pe_i, pszv_i, ll_ice_present )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_roundoff_dyn_thd ***
      !!
      !! ** Purpose :   Remove negative sea ice values arising from roundoff errors
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i       ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i     ! ice salt content
      LOGICAL,  DIMENSION(jpi,jpj),            INTENT(in)    ::   ll_ice_present
      !!-------------------------------------------------------------------
      INTEGER :: ji, jj, jk, jl
      !!-------------------------------------------------------------------
      !$acc data present( pa_i, pv_i, pv_s, psv_i, poa_i, pe_s, pe_i, pszv_i, ll_ice_present )
      !!
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF (ll_ice_present(ji,jj)) THEN
               !
               !$acc loop seq
               DO jl=1, jpl
                  pa_i (ji,jj,jl) = MAX( pa_i(ji,jj,jl), 0._wp)
                  pv_i (ji,jj,jl) = MAX( pv_i(ji,jj,jl), 0._wp)
                  pv_s (ji,jj,jl) = MAX( pv_s(ji,jj,jl), 0._wp)
                  poa_i(ji,jj,jl) = MAX(poa_i(ji,jj,jl), 0._wp)
                  IF( nn_icesal /= 4 ) psv_i(ji,jj,jl) = MAX(psv_i(ji,jj,jl), 0._wp)
                  !$acc loop seq
                  DO jk=1, nlay_i
                     pe_i  (ji,jj,jk,jl) = MAX(  pe_i(ji,jj,jk,jl), 0._wp)
                     IF( nn_icesal == 4 ) pszv_i(ji,jj,jk,jl) = MAX(pszv_i(ji,jj,jk,jl), 0._wp)
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_s
                     pe_s(ji,jj,jk,jl) = MAX(pe_s(ji,jj,jk,jl), 0._wp)
                  END DO
               END DO
               !
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE ice_var_roundoff_dyn_thd


   SUBROUTINE ice_var_brine
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_brine ***
      !!
      !! ** Purpose :   computes brine volume fraction (%)
      !!                         and salinity of the brine in sea ice
      !!
      !! ** Method  : e = - 0.054 * S (ppt) / T (C)
      !!
      !! References : Vancoppenolle et al., JGR, 2007
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zt1, zt2, zt3, zs_br
      !!-------------------------------------------------------------------
      !
      v_ibr(:,:,:) = 0._wp
      DO jl = 1, jpl
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               DO jk=1, nlay_i
                  ! brine salinity
                  zt1 = t_i(ji,jj,jk,jl) - rt0
                  zt2 = zt1 * zt1
                  zt3 = zt2 * zt1
                  IF    ( nn_liquidus == 1 ) THEN
                     zs_br = - zt1 / rTmlt                                      ! --- Linear liquidus
                  ELSEIF( nn_liquidus == 2 ) THEN
                     zs_br = -18.7_wp * zt1 - 0.519_wp * zt2 - 0.00535_wp * zt3 ! --- 3rd order liquidus, VC19
                  ELSEIF( nn_liquidus == 3 ) THEN
                     zs_br = -17.6_wp * zt1 - 0.389_wp * zt2 - 0.00362_wp * zt3 ! --- Weast 71 liquidus in RJW14
                  ENDIF
                  ! brine volume fraction
                  IF( zt1 < - epsi10 )   v_ibr(ji,jj,jl) = v_ibr(ji,jj,jl) + r1_nlay_i * sz_i(ji,jj,jk,jl) / zs_br
               END DO
            END DO
         END DO
      ENDDO
      !
      ! mean brine volume fraction
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF( vt_i(ji,jj) > epsi20 ) THEN
               vm_ibr(ji,jj) = SUM( v_ibr(ji,jj,:) * v_i(ji,jj,:) ) / vt_i(ji,jj)
            ELSE
               vm_ibr(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE ice_var_brine

   SUBROUTINE ice_var_enthalpy(jl_cat, ll_ice_present)
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_enthalpy ***
      !!
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3) from temperature
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER,                     INTENT(in) :: jl_cat
      LOGICAL, DIMENSION(jpi,jpj), INTENT(in) :: ll_ice_present
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztmelts  ! local scalar
      !!-------------------------------------------------------------------
      !$acc data present( ll_ice_present, sz_i, t_i, e_i, t_s, e_s )
      !%acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF(ll_ice_present(ji,jj)) THEN
               !%acc loop seq
               DO jk = 1, nlay_i             ! Sea ice energy of melting
                  ztmelts       = - rTmlt  * sz_i(ji,jj,jk,jl_cat)
                  t_i(ji,jj,jk,jl_cat) = MIN( t_i(ji,jj,jk,jl_cat), ztmelts + rt0 ) ! Force t_i_1d to be lower than melting point => likely conservation issue
                  !   (sometimes zdf scheme produces abnormally high temperatures)
                  e_i(ji,jj,jk,jl_cat) = rhoi * ( rcpi  * ( ztmelts - ( t_i(ji,jj,jk,jl_cat) - rt0 ) )           &
                     &                   + rLfus * ( 1._wp - ztmelts / ( t_i(ji,jj,jk,jl_cat) - rt0 ) )   &
                     &                   - rcp   * ztmelts )
               END DO
            ENDIF
         END DO
      END DO
      !%acc end parallel loop

      !%acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF(ll_ice_present(ji,jj)) THEN
               !%acc loop seq
               DO jk = 1, nlay_s             ! Snow energy of melting
                  e_s(ji,jj,jk,jl_cat) = rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl_cat) ) + rLfus )
               END DO
            ENDIF
         END DO
      END DO
      !%acc end parallel loop
      !$acc end data
   END SUBROUTINE ice_var_enthalpy

   SUBROUTINE ice_var_vremap( ph_old, pts_old, pts_i )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_var_vremap  ***
      !!
      !! ** Purpose :
      !!           This routine computes new vertical grids in the ice,
      !!           and consistently redistributes temperatures and salinities
      !!           Redistribution is made so as to ensure energy/salt conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/salinities/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_i+2)                 ------------- cum1(nlay_i)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(0:nlay_i+1), INTENT(in)    ::   ph_old, pts_old  ! old tickness (m), enthlapy (J.m-2) or salt (m.g/kg)
      REAL(wp), DIMENSION(1:nlay_i)  , INTENT(inout) ::   pts_i            ! new enthlapies (J.m-3, remapped) or salt (g/kg)
      !!-------------------------------------------------------------------
      INTEGER  ::   jk0, jk1   !  old/new layer indices
      !
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zts_cum0, zh_cum0   ! old cumulative enthlapies/salinities and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zts_cum1, zh_cum1   ! new cumulative enthlapies/salinities and layers interfaces
      REAL(wp)                        ::   zhnew               ! new layers thicknesses
      !!-------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy/salt * thickness and layers interfaces
      !-------------------------------------------------------------------------------
      zts_cum0(0) = 0._wp
      zh_cum0 (0) = 0._wp
      DO jk0 = 1, nlay_i+2
         zts_cum0(jk0) = zts_cum0(jk0-1) + pts_old(jk0-1)
         zh_cum0 (jk0) = zh_cum0 (jk0-1) + ph_old(jk0-1)
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      zhnew = SUM( ph_old(0:nlay_i+1) ) * r1_nlay_i

      ! new layers interfaces
      zh_cum1(0) = 0._wp
      DO jk1 = 1, nlay_i
         zh_cum1(jk1) = zh_cum1(jk1-1) + zhnew
      END DO

      zts_cum1(0:nlay_i) = 0._wp
      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_i+2
         DO jk1 = 1, nlay_i-1
            IF( zh_cum1(jk1) <= zh_cum0(jk0) .AND. zh_cum1(jk1) > zh_cum0(jk0-1) )   THEN
               zts_cum1(jk1) = ( zts_cum0(jk0-1) * ( zh_cum0(jk0) - zh_cum1(jk1  ) ) +  &
                  &              zts_cum0(jk0  ) * ( zh_cum1(jk1) - zh_cum0(jk0-1) ) )  &
                  &            / ( zh_cum0(jk0) - zh_cum0(jk0-1) )
            ENDIF
         END DO
      END DO
      ! to ensure that total heat/salt content is strictly conserved, set:
      zts_cum1(nlay_i) = zts_cum0(nlay_i+2)

      ! new enthalpies/salinities
      DO jk1 = 1, nlay_i
         pts_i(jk1) = MAX( 0._wp, zts_cum1(jk1) - zts_cum1(jk1-1) ) / MAX( zhnew, epsi20 ) ! max for roundoff error
      END DO

   END SUBROUTINE ice_var_vremap

   FUNCTION snw_ent( ph_old, pe_old )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE snw_ent  ***
      !! ** Purpose :
      !!           This routine computes new vertical grids in the snow,
      !!           and consistently redistributes temperatures.
      !!           Redistribution is made so as to ensure to energy conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_s+1)                 ------------- cum1(nlay_s)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(0:nlay_s), INTENT(in) ::   ph_old             ! old thicknesses (m)
      REAL(wp), DIMENSION(0:nlay_s), INTENT(in) ::   pe_old             ! old enthlapies (J.m-3)
      REAL(wp), DIMENSION(1:nlay_s)             ::   snw_ent            ! new enthlapies (J.m-3, remapped)
      !
      INTEGER  :: ji         !  dummy loop indices
      INTEGER  :: jk0, jk1   !  old/new layer indices
      !
      REAL(wp), DIMENSION(0:nlay_s+1) ::   zeh_cum0, zh_cum0   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(0:nlay_s)   ::   zeh_cum1, zh_cum1   ! new cumulative enthlapies and layers interfaces
      REAL(wp)                        ::   zhnew               ! new layers thicknesses
      !!-------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy * thickness and layers interfaces
      !--------------------------------------------------------------------------
      zeh_cum0(0) = 0._wp
      zh_cum0 (0) = 0._wp
      DO jk0 = 1, nlay_s+1
         zeh_cum0(jk0) = zeh_cum0(jk0-1) + pe_old(jk0-1) * ph_old(jk0-1)
         zh_cum0 (jk0) = zh_cum0 (jk0-1) + ph_old(jk0-1)
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      zhnew = SUM( ph_old(0:nlay_s) ) * r1_nlay_s

      ! new layers interfaces
      zh_cum1(0) = 0._wp
      DO jk1 = 1, nlay_s
         zh_cum1(jk1) = zh_cum1(jk1-1) + zhnew
      END DO

      zeh_cum1(0:nlay_s) = 0._wp
      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_s+1
         DO jk1 = 1, nlay_s-1
            IF( zh_cum1(jk1) <= zh_cum0(jk0) .AND. zh_cum1(jk1) > zh_cum0(jk0-1) )   THEN
               zeh_cum1(jk1) = ( zeh_cum0(jk0-1) * ( zh_cum0(jk0) - zh_cum1(jk1  ) ) +  &
                  &              zeh_cum0(jk0  ) * ( zh_cum1(jk1) - zh_cum0(jk0-1) ) )  &
                  &            / ( zh_cum0(jk0) - zh_cum0(jk0-1) )
            ENDIF
         END DO
      END DO
      ! to ensure that total heat content is strictly conserved, set:
      zeh_cum1(nlay_s) = zeh_cum0(nlay_s+1)

      ! new enthalpies
      DO jk1 = 1, nlay_s
         snw_ent(jk1) = MAX( 0._wp, zeh_cum1(jk1) - zeh_cum1(jk1-1) ) / MAX( zhnew, epsi20 ) ! max for roundoff error
      END DO

   END FUNCTION snw_ent



   SUBROUTINE ice_var_sshdyn( pssh, psnwice_mass_b, psshdyn )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_sshdyn  ***
      !!
      !! ** Purpose :  compute the equivalent ssh in lead when sea ice is embedded
      !!
      !! ** Method  :  ssh_lead = ssh + (Mice + Msnow) / rho0
      !!
      !! ** Reference : Jean-Michel Campin, John Marshall, David Ferreira,
      !!                Sea ice-ocean coupling using a rescaled vertical coordinate z*,
      !!                Ocean Modelling, Volume 24, Issues 1-2, 2008
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pssh            !: ssh [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: psnwice_mass_b  !: mass of snow and ice at previous ice time step [Kg/m2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: psshdyn         ! [m]
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj
      !!---------------------------------------------------------------------
      !$acc data present( pssh, psnwice_mass_b, psshdyn )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ! compute ice load used to define the equivalent ssh in lead
            psshdyn(ji,jj) = MERGE( pssh(ji,jj) + psnwice_mass_b(ji,jj)*r1_rho0,  pssh(ji,jj),  ln_ice_embd )
            !
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE ice_var_sshdyn


   !!-------------------------------------------------------------------
   !!                ***  INTERFACE ice_var_itd   ***
   !!
   !! ** Purpose :  converting N-cat ice to jpl ice categories
   !!-------------------------------------------------------------------
   SUBROUTINE ice_var_itd_1c1c( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !! ** Purpose :  converting 1-cat ice to 1 ice category
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !!-------------------------------------------------------------------
      ! == thickness and concentration == !
      ph_i(:) = phti(:)
      ph_s(:) = phts(:)
      pa_i(:) = pati(:)
      !
      ! == temperature and salinity and ponds == !
      pt_i (:) = ptmi (:)
      pt_s (:) = ptms (:)
      pt_su(:) = ptmsu(:)
      ps_i (:) = psmi (:)
      pa_ip(:) = patip(:)
      ph_ip(:) = phtip(:)
      ph_il(:) = phtil(:)

   END SUBROUTINE ice_var_itd_1c1c

   SUBROUTINE ice_var_itd_Nc1c( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !! ** Purpose :  converting N-cat ice to 1 ice category
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:)  , INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:)  , INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   z1_ai, z1_vi, z1_vs
      !
      INTEGER ::   idim
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti, 1 )
      !
      ! == thickness and concentration == !
      ALLOCATE( z1_ai(idim), z1_vi(idim), z1_vs(idim) )
      !
      pa_i(:) = SUM( pati(:,:), dim=2 )

      WHERE( ( pa_i(:) ) /= 0._wp )
         z1_ai(:) = 1._wp / pa_i(:)
      ELSEWHERE
         z1_ai(:) = 0._wp
      END WHERE

      ph_i(:) = SUM( phti(:,:) * pati(:,:), dim=2 ) * z1_ai(:)
      ph_s(:) = SUM( phts(:,:) * pati(:,:), dim=2 ) * z1_ai(:)
      !
      ! == temperature and salinity == !
      WHERE( ( pa_i(:) * ph_i(:) ) /= 0._wp )
         z1_vi(:) = 1._wp / ( pa_i(:) * ph_i(:) )
      ELSEWHERE
         z1_vi(:) = 0._wp
      END WHERE
      WHERE( ( pa_i(:) * ph_s(:) ) /= 0._wp )
         z1_vs(:) = 1._wp / ( pa_i(:) * ph_s(:) )
      ELSEWHERE
         z1_vs(:) = 0._wp
      END WHERE
      pt_i (:) = SUM( ptmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
      pt_s (:) = SUM( ptms (:,:) * pati(:,:) * phts(:,:), dim=2 ) * z1_vs(:)
      pt_su(:) = SUM( ptmsu(:,:) * pati(:,:)            , dim=2 ) * z1_ai(:)
      ps_i (:) = SUM( psmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)

      ! == ponds == !
      pa_ip(:) = SUM( patip(:,:), dim=2 )
      WHERE( pa_ip(:) /= 0._wp )
         ph_ip(:) = SUM( phtip(:,:) * patip(:,:), dim=2 ) / pa_ip(:)
         ph_il(:) = SUM( phtil(:,:) * patip(:,:), dim=2 ) / pa_ip(:)
      ELSEWHERE
         ph_ip(:) = 0._wp
         ph_il(:) = 0._wp
      END WHERE
      !
      DEALLOCATE( z1_ai, z1_vi, z1_vs )
      !
   END SUBROUTINE ice_var_itd_Nc1c

   SUBROUTINE ice_var_itd_1cMc( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :  converting 1-cat ice to jpl ice categories
      !!
      !!
      !! ** Method:   ice thickness distribution follows a gamma function from Abraham et al. (2015)
      !!              it has the property of conserving total concentration and volume
      !!
      !!
      !! ** Arguments : phti: 1-cat ice thickness
      !!                phts: 1-cat snow depth
      !!                pati: 1-cat ice concentration
      !!
      !! ** Output    : jpl-cat
      !!
      !!  Abraham, C., Steiner, N., Monahan, A. and Michel, C., 2015.
      !!               Effects of subgrid‐scale snow thickness variability on radiative transfer in sea ice.
      !!               Journal of Geophysical Research: Oceans, 120(8), pp.5597-5614
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:),   INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:)  , INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zfra, z1_hti
      INTEGER  ::   ji, jk, jl
      INTEGER  ::   idim
      REAL(wp) ::   zv, zdh
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti , 1 )
      !
      ph_i(1:idim,1:jpl) = 0._wp
      ph_s(1:idim,1:jpl) = 0._wp
      pa_i(1:idim,1:jpl) = 0._wp
      !
      ALLOCATE( z1_hti(idim) )
      WHERE( phti(:) /= 0._wp )
         z1_hti(:) = 1._wp / phti(:)
      ELSEWHERE
         z1_hti(:) = 0._wp
      END WHERE
      !
      ! == thickness and concentration == !
      ! for categories 1:jpl-1, integrate the gamma function from hi_max(jl-1) to hi_max(jl)
      DO jl = 1, jpl-1
         DO ji = 1, idim
            !
            IF( phti(ji) > 0._wp ) THEN
               ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
               pa_i(ji,jl) = pati(ji) * z1_hti(ji) * (  ( phti(ji) + 2.*hi_max(jl-1) ) * EXP( -2.*hi_max(jl-1)*z1_hti(ji) ) &
                  &                                   - ( phti(ji) + 2.*hi_max(jl  ) ) * EXP( -2.*hi_max(jl  )*z1_hti(ji) ) )
               !
               ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
               zv = pati(ji) * z1_hti(ji) * (  ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jl-1) + 2.*hi_max(jl-1)*hi_max(jl-1) ) &
                  &                            * EXP( -2.*hi_max(jl-1)*z1_hti(ji) ) &
                  &                          - ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jl) + 2.*hi_max(jl)*hi_max(jl) ) &
                  &                            * EXP(-2.*hi_max(jl)*z1_hti(ji)) )
               ! thickness
               IF( pa_i(ji,jl) > epsi06 ) THEN
                  ph_i(ji,jl) = zv / pa_i(ji,jl)
               ELSE
                  ph_i(ji,jl) = 0.
                  pa_i(ji,jl) = 0.
               ENDIF
            ENDIF
            !
         ENDDO
      ENDDO
      !
      ! for the last category (jpl), integrate the gamma function from hi_max(jpl-1) to infinity
      DO ji = 1, idim
         !
         IF( phti(ji) > 0._wp ) THEN
            ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jpl-1) to infinity
            pa_i(ji,jpl) = pati(ji) * z1_hti(ji) * ( phti(ji) + 2.*hi_max(jpl-1) ) * EXP( -2.*hi_max(jpl-1)*z1_hti(ji) )

            ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jpl-1) to infinity
            zv = pati(ji) * z1_hti(ji) * ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jpl-1) + 2.*hi_max(jpl-1)*hi_max(jpl-1) ) &
               &                         * EXP( -2.*hi_max(jpl-1)*z1_hti(ji) )
            ! thickness
            IF( pa_i(ji,jpl) > epsi06 ) THEN
               ph_i(ji,jpl) = zv / pa_i(ji,jpl)
            else
               ph_i(ji,jpl) = 0.
               pa_i(ji,jpl) = 0.
            ENDIF
         ENDIF
         !
      ENDDO
      !
      ! Add Snow in each category where pa_i is not 0
      DO jl = 1, jpl
         DO ji = 1, idim
            IF( pa_i(ji,jl) > 0._wp ) THEN
               ph_s(ji,jl) = ph_i(ji,jl) * phts(ji) * z1_hti(ji)
               ! In case snow load is in excess that would lead to transformation from snow to ice
               ! Then, transfer the snow excess into the ice (different from icethd_dh)
               zdh = MAX( 0._wp, ( rhos * ph_s(ji,jl) + ( rhoi - rho0 ) * ph_i(ji,jl) ) * r1_rho0 )
               ! recompute h_i, h_s avoiding out of bounds values
               ph_i(ji,jl) = MIN( hi_max(jl), ph_i(ji,jl) + zdh )
               ph_s(ji,jl) = MAX( 0._wp, ph_s(ji,jl) - zdh * rhoi * r1_rhos )
            ENDIF
         END DO
      END DO
      !
      DEALLOCATE( z1_hti )
      !
      ! == temperature and salinity == !
      DO jl = 1, jpl
         pt_i (:,jl) = ptmi (:)
         pt_s (:,jl) = ptms (:)
         pt_su(:,jl) = ptmsu(:)
         ps_i (:,jl) = psmi (:)
      END DO
      !
      ! == ponds == !
      ALLOCATE( zfra(idim) )
      ! keep the same pond fraction atip/ati for each category
      WHERE( pati(:) /= 0._wp )
         zfra(:) = patip(:) / pati(:)
      ELSEWHERE
         zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         pa_ip(:,jl) = zfra(:) * pa_i(:,jl)
      END DO
      ! keep the same v_ip/v_i ratio for each category
      WHERE( ( phti(:) * pati(:) ) /= 0._wp )
         zfra(:) = ( phtip(:) * patip(:) ) / ( phti(:) * pati(:) )
      ELSEWHERE
         zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,jl) /= 0._wp )
            ph_ip(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
         ELSEWHERE
            ph_ip(:,jl) = 0._wp
         END WHERE
      END DO
      ! keep the same v_il/v_i ratio for each category
      WHERE( ( phti(:) * pati(:) ) /= 0._wp )
         zfra(:) = ( phtil(:) * patip(:) ) / ( phti(:) * pati(:) )
      ELSEWHERE
         zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,jl) /= 0._wp )
            ph_il(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
         ELSEWHERE
            ph_il(:,jl) = 0._wp
         END WHERE
      END DO
      DEALLOCATE( zfra )
      !
   END SUBROUTINE ice_var_itd_1cMc

   SUBROUTINE ice_var_itd_NcMc( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :  converting N-cat ice to jpl ice categories
      !!
      !!                  ice thickness distribution follows a gaussian law
      !!               around the concentration of the most likely ice thickness
      !!                           (similar as iceistate.F90)
      !!
      !! ** Method:   Iterative procedure
      !!
      !!               1) Fill ice cat that correspond to input thicknesses
      !!                  Find the lowest(jlmin) and highest(jlmax) cat that are filled
      !!
      !!               2) Expand the filling to the cat jlmin-1 and jlmax+1
      !!                   by removing 25% ice area from jlmin and jlmax (resp.)
      !!
      !!               3) Expand the filling to the empty cat between jlmin and jlmax
      !!                   by a) removing 25% ice area from the lower cat (ascendant loop jlmin=>jlmax)
      !!                      b) removing 25% ice area from the higher cat (descendant loop jlmax=>jlmin)
      !!
      !! ** Arguments : phti: N-cat ice thickness
      !!                phts: N-cat snow depth
      !!                pati: N-cat ice concentration
      !!
      !! ** Output    : jpl-cat
      !!
      !!  (Example of application: BDY forcings when inputs have N-cat /= jpl)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      INTEGER , ALLOCATABLE, DIMENSION(:,:) ::   jlfil, jlfil2
      INTEGER , ALLOCATABLE, DIMENSION(:)   ::   jlmax, jlmin
      REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   z1_ai, z1_vi, z1_vs, ztmp, zfra
      !
      REAL(wp), PARAMETER ::   ztrans = 0.25_wp
      INTEGER  ::   ji, jl, jl1, jl2
      INTEGER  ::   idim, icat
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti, 1 )
      icat = SIZE( phti, 2 )
      !
      ! == thickness and concentration == !
      !                                 ! ---------------------- !
      IF( icat == jpl ) THEN            ! input cat = output cat !
         !                              ! ---------------------- !
         ph_i(:,:) = phti(:,:)
         ph_s(:,:) = phts(:,:)
         pa_i(:,:) = pati(:,:)
         !
         ! == temperature and salinity and ponds == !
         pt_i (:,:) = ptmi (:,:)
         pt_s (:,:) = ptms (:,:)
         pt_su(:,:) = ptmsu(:,:)
         ps_i (:,:) = psmi (:,:)
         pa_ip(:,:) = patip(:,:)
         ph_ip(:,:) = phtip(:,:)
         ph_il(:,:) = phtil(:,:)
         !                              ! ---------------------- !
      ELSEIF( icat == 1 ) THEN          ! input cat = 1          !
         !                              ! ---------------------- !
         CALL  ice_var_itd_1cMc( phti(:,1), phts(:,1), pati (:,1), &
            &                    ph_i(:,:), ph_s(:,:), pa_i (:,:), &
            &                    ptmi(:,1), ptms(:,1), ptmsu(:,1), psmi(:,1), patip(:,1), phtip(:,1), phtil(:,1), &
            &                    pt_i(:,:), pt_s(:,:), pt_su(:,:), ps_i(:,:), pa_ip(:,:), ph_ip(:,:), ph_il(:,:)  )
         !                              ! ---------------------- !
      ELSEIF( jpl == 1 ) THEN           ! output cat = 1         !
         !                              ! ---------------------- !
         CALL  ice_var_itd_Nc1c( phti(:,:), phts(:,:), pati (:,:), &
            &                    ph_i(:,1), ph_s(:,1), pa_i (:,1), &
            &                    ptmi(:,:), ptms(:,:), ptmsu(:,:), psmi(:,:), patip(:,:), phtip(:,:), phtil(:,:), &
            &                    pt_i(:,1), pt_s(:,1), pt_su(:,1), ps_i(:,1), pa_ip(:,1), ph_ip(:,1), ph_il(:,1)  )
         !                              ! ----------------------- !
      ELSE                              ! input cat /= output cat !
         !                              ! ----------------------- !

         ALLOCATE( jlfil(idim,jpl), jlfil2(idim,jpl) )       ! allocate arrays
         ALLOCATE( jlmin(idim), jlmax(idim) )

         ! --- initialize output fields to 0 --- !
         ph_i(1:idim,1:jpl) = 0._wp
         ph_s(1:idim,1:jpl) = 0._wp
         pa_i(1:idim,1:jpl) = 0._wp
         !
         ! --- fill the categories --- !
         !     find where cat-input = cat-output and fill cat-output fields
         jlmax(:) = 0
         jlmin(:) = 999
         jlfil(:,:) = 0
         DO jl1 = 1, jpl
            DO jl2 = 1, icat
               DO ji = 1, idim
                  IF( hi_max(jl1-1) <= phti(ji,jl2) .AND. hi_max(jl1) > phti(ji,jl2) ) THEN
                     ! fill the right category
                     ph_i(ji,jl1) = phti(ji,jl2)
                     ph_s(ji,jl1) = phts(ji,jl2)
                     pa_i(ji,jl1) = pati(ji,jl2)
                     ! record categories that are filled
                     jlmax(ji) = MAX( jlmax(ji), jl1 )
                     jlmin(ji) = MIN( jlmin(ji), jl1 )
                     jlfil(ji,jl1) = jl1
                  ENDIF
               END DO
            END DO
         END DO
         !
         ! --- fill the gaps between categories --- !
         !     transfer from categories filled at the previous step to the empty ones in between
         DO ji = 1, idim
            jl1 = jlmin(ji)
            jl2 = jlmax(ji)
            IF( jl1 > 1 ) THEN
               ! fill the lower cat (jl1-1)
               pa_i(ji,jl1-1) = ztrans * pa_i(ji,jl1)
               ph_i(ji,jl1-1) = hi_mean(jl1-1)
               ! remove from cat jl1
               pa_i(ji,jl1  ) = ( 1._wp - ztrans ) * pa_i(ji,jl1)
            ENDIF
            IF( jl2 < jpl ) THEN
               ! fill the upper cat (jl2+1)
               pa_i(ji,jl2+1) = ztrans * pa_i(ji,jl2)
               ph_i(ji,jl2+1) = hi_mean(jl2+1)
               ! remove from cat jl2
               pa_i(ji,jl2  ) = ( 1._wp - ztrans ) * pa_i(ji,jl2)
            ENDIF
         END DO
         !
         jlfil2(:,:) = jlfil(:,:)
         ! fill categories from low to high
         DO jl = 2, jpl-1
            DO ji = 1, idim
               IF( jlfil(ji,jl-1) /= 0 .AND. jlfil(ji,jl) == 0 ) THEN
                  ! fill high
                  pa_i(ji,jl) = ztrans * pa_i(ji,jl-1)
                  ph_i(ji,jl) = hi_mean(jl)
                  jlfil(ji,jl) = jl
                  ! remove low
                  pa_i(ji,jl-1) = ( 1._wp - ztrans ) * pa_i(ji,jl-1)
               ENDIF
            END DO
         END DO
         !
         ! fill categories from high to low
         DO jl = jpl-1, 2, -1
            DO ji = 1, idim
               IF( jlfil2(ji,jl+1) /= 0 .AND. jlfil2(ji,jl) == 0 ) THEN
                  ! fill low
                  pa_i(ji,jl) = pa_i(ji,jl) + ztrans * pa_i(ji,jl+1)
                  ph_i(ji,jl) = hi_mean(jl)
                  jlfil2(ji,jl) = jl
                  ! remove high
                  pa_i(ji,jl+1) = ( 1._wp - ztrans ) * pa_i(ji,jl+1)
               ENDIF
            END DO
         END DO
         !
         DEALLOCATE( jlfil, jlfil2 )      ! deallocate arrays
         DEALLOCATE( jlmin, jlmax )
         !
         ! == temperature and salinity == !
         !
         ALLOCATE( z1_ai(idim), z1_vi(idim), z1_vs(idim), ztmp(idim) )
         !
         WHERE( SUM( pa_i(:,:), dim=2 ) /= 0._wp )
            z1_ai(:) = 1._wp / SUM( pa_i(:,:), dim=2 )
         ELSEWHERE
            z1_ai(:) = 0._wp
         END WHERE
         WHERE( SUM( pa_i(:,:) * ph_i(:,:), dim=2 ) /= 0._wp )
            z1_vi(:) = 1._wp / SUM( pa_i(:,:) * ph_i(:,:), dim=2 )
         ELSEWHERE
            z1_vi(:) = 0._wp
         END WHERE
         WHERE( SUM( pa_i(:,:) * ph_s(:,:), dim=2 ) /= 0._wp )
            z1_vs(:) = 1._wp / SUM( pa_i(:,:) * ph_s(:,:), dim=2 )
         ELSEWHERE
            z1_vs(:) = 0._wp
         END WHERE
         !
         ! fill all the categories with the same value
         ztmp(:) = SUM( ptmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
         DO jl = 1, jpl
            pt_i (:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( ptms (:,:) * pati(:,:) * phts(:,:), dim=2 ) * z1_vs(:)
         DO jl = 1, jpl
            pt_s (:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( ptmsu(:,:) * pati(:,:)            , dim=2 ) * z1_ai(:)
         DO jl = 1, jpl
            pt_su(:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( psmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
         DO jl = 1, jpl
            ps_i (:,jl) = ztmp(:)
         END DO
         !
         DEALLOCATE( z1_ai, z1_vi, z1_vs, ztmp )
         !
         ! == ponds == !
         ALLOCATE( zfra(idim) )
         ! keep the same pond fraction atip/ati for each category
         WHERE( SUM( pati(:,:), dim=2 ) /= 0._wp )
            zfra(:) = SUM( patip(:,:), dim=2 ) / SUM( pati(:,:), dim=2 )
         ELSEWHERE
            zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            pa_ip(:,jl) = zfra(:) * pa_i(:,jl)
         END DO
         ! keep the same v_ip/v_i ratio for each category
         WHERE( SUM( phti(:,:) * pati(:,:), dim=2 ) /= 0._wp )
            zfra(:) = SUM( phtip(:,:) * patip(:,:), dim=2 ) / SUM( phti(:,:) * pati(:,:), dim=2 )
         ELSEWHERE
            zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            WHERE( pa_ip(:,jl) /= 0._wp )
               ph_ip(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
            ELSEWHERE
               ph_ip(:,jl) = 0._wp
            END WHERE
         END DO
         ! keep the same v_il/v_i ratio for each category
         WHERE( SUM( phti(:,:) * pati(:,:), dim=2 ) /= 0._wp )
            zfra(:) = SUM( phtil(:,:) * patip(:,:), dim=2 ) / SUM( phti(:,:) * pati(:,:), dim=2 )
         ELSEWHERE
            zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            WHERE( pa_ip(:,jl) /= 0._wp )
               ph_il(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
            ELSEWHERE
               ph_il(:,jl) = 0._wp
            END WHERE
         END DO
         DEALLOCATE( zfra )
         !
      ENDIF
      !
   END SUBROUTINE ice_var_itd_NcMc


   SUBROUTINE ice_var_itd_1cMc_2d( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !!      ==> 2D version of `ice_var_itd_1cMc` for 2D thermo (with )
      !!
      !! ** Purpose :  converting 1-cat ice to jpl ice categories
      !!
      !!
      !! ** Method:   ice thickness distribution follows a gamma function from Abraham et al. (2015)
      !!              it has the property of conserving total concentration and volume
      !!
      !!
      !! ** Arguments : phti: 1-cat ice thickness
      !!                phts: 1-cat snow depth
      !!                pati: 1-cat ice concentration
      !!
      !! ** Output    : jpl-cat
      !!
      !!  Abraham, C., Steiner, N., Monahan, A. and Michel, C., 2015.
      !!               Effects of subgrid‐scale snow thickness variability on radiative transfer in sea ice.
      !!               Journal of Geophysical Research: Oceans, 120(8), pp.5597-5614
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(jpi,jpj)  ,   INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      REAL(wp), DIMENSION(jpi,jpj) ::  z1_hti, zfra
      INTEGER  ::   ji, jj, jk, jl
      REAL(wp) ::   zv, zdh
      !!-------------------------------------------------------------------
      !
      ph_i(:,:,:) = 0._wp
      ph_s(:,:,:) = 0._wp
      pa_i(:,:,:) = 0._wp
      !
      z1_hti(:,:) = 0._wp
      WHERE( phti(:,:) /= 0._wp ) z1_hti(:,:) = 1._wp / phti(:,:)
      !
      ! == thickness and concentration == !
      ! for categories 1:jpl-1, integrate the gamma function from hi_max(jl-1) to hi_max(jl)
      DO jl = 1, jpl-1
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !
               IF( phti(ji,jj) > 0._wp ) THEN
                  ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
                  pa_i(ji,jj,jl) = pati(ji,jj) * z1_hti(ji,jj) * (  ( phti(ji,jj) + 2.*hi_max(jl-1) ) * EXP( -2.*hi_max(jl-1)*z1_hti(ji,jj) ) &
                     &                                   - ( phti(ji,jj) + 2.*hi_max(jl  ) ) * EXP( -2.*hi_max(jl  )*z1_hti(ji,jj) ) )
                  !
                  ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
                  zv = pati(ji,jj) * z1_hti(ji,jj) * (  ( phti(ji,jj)*phti(ji,jj) + 2.*phti(ji,jj)*hi_max(jl-1) + 2.*hi_max(jl-1)*hi_max(jl-1) ) &
                     &                            * EXP( -2.*hi_max(jl-1)*z1_hti(ji,jj) ) &
                     &                          - ( phti(ji,jj)*phti(ji,jj) + 2.*phti(ji,jj)*hi_max(jl) + 2.*hi_max(jl)*hi_max(jl) ) &
                     &                            * EXP(-2.*hi_max(jl)*z1_hti(ji,jj)) )
                  ! thickness
                  IF( pa_i(ji,jj,jl) > epsi06 ) THEN
                     ph_i(ji,jj,jl) = zv / pa_i(ji,jj,jl)
                  ELSE
                     ph_i(ji,jj,jl) = 0.
                     pa_i(ji,jj,jl) = 0.
                  ENDIF
               ENDIF
               !
            ENDDO
         ENDDO
      ENDDO !DO jl = 1, jpl-1
      !
      ! for the last category (jpl), integrate the gamma function from hi_max(jpl-1) to infinity
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !
            IF( phti(ji,jj) > 0._wp ) THEN
               ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jpl-1) to infinity
               pa_i(ji,jj,jpl) = pati(ji,jj) * z1_hti(ji,jj) * ( phti(ji,jj) + 2.*hi_max(jpl-1) ) * EXP( -2.*hi_max(jpl-1)*z1_hti(ji,jj) )

               ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jpl-1) to infinity
               zv = pati(ji,jj) * z1_hti(ji,jj) * ( phti(ji,jj)*phti(ji,jj) + 2.*phti(ji,jj)*hi_max(jpl-1) + 2.*hi_max(jpl-1)*hi_max(jpl-1) ) &
                  &                         * EXP( -2.*hi_max(jpl-1)*z1_hti(ji,jj) )
               ! thickness
               IF( pa_i(ji,jj,jpl) > epsi06 ) THEN
                  ph_i(ji,jj,jpl) = zv / pa_i(ji,jj,jpl)
               else
                  ph_i(ji,jj,jpl) = 0.
                  pa_i(ji,jj,jpl) = 0.
               ENDIF
            ENDIF
            !
         ENDDO
      ENDDO
      !
      ! Add Snow in each category where pa_i is not 0
      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1

               IF( pa_i(ji,jj,jl) > 0._wp ) THEN
                  ph_s(ji,jj,jl) = ph_i(ji,jj,jl) * phts(ji,jj) * z1_hti(ji,jj)
                  ! In case snow load is in excess that would lead to transformation from snow to ice
                  ! Then, transfer the snow excess into the ice (different from icethd_dh)
                  zdh = MAX( 0._wp, ( rhos * ph_s(ji,jj,jl) + ( rhoi - rho0 ) * ph_i(ji,jj,jl) ) * r1_rho0 )
                  ! recompute h_i, h_s avoiding out of bounds values
                  ph_i(ji,jj,jl) = MIN( hi_max(jl), ph_i(ji,jj,jl) + zdh )
                  ph_s(ji,jj,jl) = MAX( 0._wp, ph_s(ji,jj,jl) - zdh * rhoi * r1_rhos )
               ENDIF
            END DO
         END DO
      END DO
      !
      ! == temperature and salinity == !
      DO jl = 1, jpl
         pt_i (:,:,jl) = ptmi (:,:)
         pt_s (:,:,jl) = ptms (:,:)
         pt_su(:,:,jl) = ptmsu(:,:)
         ps_i (:,:,jl) = psmi (:,:)
      END DO
      !
      ! == ponds == !
      ! keep the same pond fraction atip/ati for each category
      WHERE( pati(:,:) /= 0._wp )
         zfra(:,:) = patip(:,:) / pati(:,:)
      ELSEWHERE
         zfra(:,:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         pa_ip(:,:,jl) = zfra(:,:) * pa_i(:,:,jl)
      END DO
      ! keep the same v_ip/v_i ratio for each category
      WHERE( ( phti(:,:) * pati(:,:) ) /= 0._wp )
         zfra(:,:) = ( phtip(:,:) * patip(:,:) ) / ( phti(:,:) * pati(:,:) )
      ELSEWHERE
         zfra(:,:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,:,jl) /= 0._wp )
            ph_ip(:,:,jl) = zfra(:,:) * ( ph_i(:,:,jl) * pa_i(:,:,jl) ) / pa_ip(:,:,jl)
         ELSEWHERE
            ph_ip(:,:,jl) = 0._wp
         END WHERE
      END DO
      ! keep the same v_il/v_i ratio for each category
      WHERE( ( phti(:,:) * pati(:,:) ) /= 0._wp )
         zfra(:,:) = ( phtil(:,:) * patip(:,:) ) / ( phti(:,:) * pati(:,:) )
      ELSEWHERE
         zfra(:,:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,:,jl) /= 0._wp )
            ph_il(:,:,jl) = zfra(:,:) * ( ph_i(:,:,jl) * pa_i(:,:,jl) ) / pa_ip(:,:,jl)
         ELSEWHERE
            ph_il(:,:,jl) = 0._wp
         END WHERE
      END DO
      !
   END SUBROUTINE ice_var_itd_1cMc_2d



   !!-------------------------------------------------------------------
   !! INTERFACE ice_var_snwfra
   !!
   !! ** Purpose :  fraction of ice covered by snow
   !!
   !! ** Method  :  In absence of proper snow model on top of sea ice,
   !!               we argue that snow does not cover the whole ice because
   !!               of wind blowing...
   !!
   !! ** Arguments : ph_s: snow thickness
   !!
   !! ** Output    : pa_s_fra: fraction of ice covered by snow
   !!
   !!-------------------------------------------------------------------
   SUBROUTINE ice_var_snwfra_3d( ph_s, pa_s_fra )
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ph_s        ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pa_s_fra    ! ice fraction covered by snow
      !IF(     nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
      !   pa_s_fra = 1._wp - EXP( -0.2_wp * rhos * ph_s )
      !ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
      pa_s_fra = ph_s / ( ph_s + 0.02_wp )
      !ENDIF
   END SUBROUTINE ice_var_snwfra_3d

   SUBROUTINE ice_var_snwfra_2d( ph_s, pa_s_fra )
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   ph_s        ! snow thickness
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pa_s_fra    ! ice fraction covered by snow
      !IF(     nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
      !   pa_s_fra = 1._wp - EXP( -0.2_wp * rhos * ph_s )
      !ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
      pa_s_fra = ph_s / ( ph_s + 0.02_wp )
      !ENDIF
   END SUBROUTINE ice_var_snwfra_2d

   SUBROUTINE ice_var_snwfra_1d( ph_s, pa_s_fra )
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   ph_s        ! snow thickness
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pa_s_fra    ! ice fraction covered by snow
      !IF(     nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
      !   pa_s_fra = 1._wp - EXP( -0.2_wp * rhos * ph_s )
      !ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
      pa_s_fra = ph_s / ( ph_s + 0.02_wp )
      !ENDIF
   END SUBROUTINE ice_var_snwfra_1d

   SUBROUTINE ice_var_snwfra_sclr( ph_s, pa_s_fra )
      !$acc routine
      REAL(wp), INTENT(in   ) :: ph_s        ! snow thickness
      REAL(wp), INTENT(  out) :: pa_s_fra    ! ice fraction covered by snow
      !IF(     nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
      !   pa_s_fra = 1._wp - EXP( -0.2_wp * rhos * ph_s )
      !ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
      pa_s_fra = ph_s / ( ph_s + 0.02_wp )
      !ENDIF
   END SUBROUTINE ice_var_snwfra_sclr

   !!--------------------------------------------------------------------------
   !! INTERFACE ice_var_snwblow
   !!
   !! ** Purpose :   Compute distribution of precip over the ice
   !!
   !!                Snow accumulation in one thermodynamic time step
   !!                snowfall is partitionned between leads and ice.
   !!                If snow fall was uniform, a fraction (1-at_i) would fall into leads
   !!                but because of the winds, more snow falls on leads than on sea ice
   !!                and a greater fraction (1-at_i)^beta of the total mass of snow
   !!                (beta < 1) falls in leads.
   !!                In reality, beta depends on wind speed,
   !!                and should decrease with increasing wind speed but here, it is
   !!                considered as a constant. an average value is 0.66
   !!--------------------------------------------------------------------------
   !!gm  I think it can be usefull to set this as a FUNCTION, not a SUBROUTINE....
   SUBROUTINE ice_var_snwblow_2d( pin, pout )
      REAL(wp), DIMENSION(:,:), INTENT(in   ) :: pin   ! previous fraction lead ( 1. - a_i_b )
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pout
      pout = ( 1._wp - MAX(pin,0._wp)**rn_snwblow )
   END SUBROUTINE ice_var_snwblow_2d

   SUBROUTINE ice_var_snwblow_1d( pin, pout )
      REAL(wp), DIMENSION(:), INTENT(in   ) :: pin
      REAL(wp), DIMENSION(:), INTENT(inout) :: pout
      pout = ( 1._wp - MAX(pin,0._wp)**rn_snwblow )
   END SUBROUTINE ice_var_snwblow_1d



   SUBROUTINE ice_var_hpiling
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_var_hpiling  ***
      !!
      !! ** Purpose : Simple conservative piling comparable with 1-cat models
      !!
      !! ** Method  : pile-up ice when no ridging/rafting
      !!
      !! ** input   : a_i
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      IF( ln_timing )  CALL timing_start('ice_var_hpiling')

      !$acc parallel loop collapse(2) present( at_i, a_i )
      DO jj = Njs0-nn_hls, Nje0+nn_hls
         DO ji = Nis0-nn_hls, Nie0+nn_hls
            !
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
            !$acc loop seq
            DO jl = 1, jpl
               IF( at_i(ji,jj) > epsi20 ) THEN
                  a_i(ji,jj,jl) = a_i(ji,jj,jl) * (  1._wp + MIN( rn_amax - at_i(ji,jj) , 0._wp ) / at_i(ji,jj)  )
               ENDIF
            END DO
            !
         END DO
      END DO
      !$acc end parallel loop

      IF( ln_timing )  CALL timing_stop('ice_var_hpiling')
      !
   END SUBROUTINE ice_var_hpiling

   SUBROUTINE ice_var_cap_at
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_var_cap_at  ***
      !!
      !! ** Purpose : Correct `a_i` so that `at_i` cannot overshoot `rn_amax`
      !!
      !! ** input   : a_i
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      !$acc parallel loop collapse(2) present( at_i, a_i )
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
            !$acc loop seq
            DO jl=1, jpl
               IF( at_i(ji,jj) > rn_amax )   a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / at_i(ji,jj)
            END DO
         END DO
      END DO
      !$acc end parallel loop
   END SUBROUTINE ice_var_cap_at


   SUBROUTINE test4inf_2d( cstr, px )
      !!-------------------------------------------------------------------
      CHARACTER(len=22),            INTENT(in) :: cstr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, icpt!, jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      !$acc data present( px )
      icpt = 0
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF( px(ji,jj) > HUGE(px(ji,jj)) ) THEN
               !PRINT *, ' *** Infinite value for at ji,jj=',ji,jj, cstr
               icpt = icpt + 1
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      IF( icpt > 0 ) THEN
         PRINT *, ' *** Infinite value for: ', cstr, icpt
         STOP
      ENDIF
      !$acc end data
   END SUBROUTINE test4inf_2d

   SUBROUTINE test4inf_3d( cstr, px )
      !!-------------------------------------------------------------------
      CHARACTER(len=22),                INTENT(in) :: cstr
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) :: px
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, jl, icpt         ! dummy loop indices
      !!-------------------------------------------------------------------
      !$acc data present( px )
      icpt = 0
      !$acc parallel loop collapse(3)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl=1, jpl
               IF( px(ji,jj,jl) > HUGE(px(ji,jj,jl)) ) THEN
                  !PRINT *, ' *** Infinite value for at ji,jj=',ji,jj, cstr
                  icpt = icpt + 1
               ENDIF
            END DO
         END DO
      END DO
      !$acc end parallel loop

      IF( icpt > 0 ) THEN
         PRINT *, ' *** Infinite value for: ', cstr, icpt
         STOP
      ENDIF
      !$acc end data
   END SUBROUTINE test4inf_3d

   SUBROUTINE test4nan( cstr, px )
      !!-------------------------------------------------------------------
      CHARACTER(len=22),            INTENT(in) :: cstr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, icpt!, jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      !$acc data present( px )
      icpt = 0
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF( px(ji,jj) /= px(ji,jj) ) THEN
               !PRINT *, ' *** Infinite value for at ji,jj=',ji,jj, cstr
               icpt = icpt + 1
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      IF( icpt > 0 ) THEN
         PRINT *, ' *** NaN value for: ', cstr, icpt
         STOP
      ENDIF
      !$acc end data
   END SUBROUTINE test4nan


   FUNCTION l_is_it_a_nan( px )
      !!-------------------------------------------------------------------
      LOGICAL              :: l_is_it_a_nan
      REAL(wp), INTENT(in) :: px
      !!-------------------------------------------------------------------
      l_is_it_a_nan = ( px /= px )
   END FUNCTION l_is_it_a_nan


   FUNCTION l_is_it_a_inf( px )
      !!-------------------------------------------------------------------
      LOGICAL              :: l_is_it_a_inf
      REAL(wp), INTENT(in) :: px
      !!-------------------------------------------------------------------
      l_is_it_a_inf = ( px > HUGE(px) )
   END FUNCTION l_is_it_a_inf

   !!======================================================================
END MODULE icevar
