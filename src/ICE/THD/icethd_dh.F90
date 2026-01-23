MODULE icethd_dh
   !!======================================================================
   !!                       ***  MODULE icethd_dh ***
   !!   seaice : thermodynamic growth and melt
   !!======================================================================
   !! History :       !  2003-05  (M. Vancoppenolle) Original code in 1D
   !!                 !  2005-06  (M. Vancoppenolle) 3D version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd_dh        : vertical sea-ice growth and melt
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE par_ice
   USE ice            ! sea-ice: variables

   USE sbc_oce , ONLY : fatm_snow
   USE sbc_ice , ONLY : qns_ice, qtr_ice_top, qsr_ice, qml_ice, qprec_ice, evap_ice
   USE oss_nnq , ONLY : sst_s, sss_s, frq_m

   USE icethd_sal     ! sea-ice: salinity profiles
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_dh        ! called by ice_thd

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_dh(jl_cat, ll_ice_present)
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_dh  ***
      !!
      !! ** Purpose :   compute ice and snow thickness changes due to growth/melting
      !!
      !! ** Method  :   Ice/Snow surface melting arises from imbalance in surface fluxes
      !!                Bottom accretion/ablation arises from flux budget
      !!                Snow thickness can increase by precipitation and decrease by sublimation
      !!                If snow load excesses Archmiede limit, snow-ice is formed by
      !!                the flooding of sea-water in the snow
      !!
      !!                - Compute available flux of heat for surface ablation
      !!                - Compute snow and sea ice enthalpies
      !!                - Surface ablation and sublimation
      !!                - Bottom accretion/ablation
      !!                - Snow ice formation
      !!
      !! ** Note     :  h=max(0,h+dh) are often used to ensure positivity of h.
      !!                very small negative values can occur otherwise (e.g. -1.e-20)
      !!
      !! References : Bitz and Lipscomb, 1999, J. Geophys. Res.
      !!              Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646
      !!              Vancoppenolle, Fichefet and Bitz, 2005, Geophys. Res. Let.
      !!              Vancoppenolle et al.,2009, Ocean Modelling
      !!------------------------------------------------------------------
      INTEGER,                     INTENT(in)    :: jl_cat        ! ice-category we are working with
      LOGICAL, DIMENSION(jpi,jpj), INTENT(inout) :: ll_ice_present
      !!------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl ! dummy loop indices
      !
      REAL(wp) ::   ztmelts      ! local scalar
      REAL(wp) ::   zdum
      REAL(wp) ::   zt_i_new     ! bottom formation temperature
      REAL(wp) ::   z1_rho       ! 1/(rhos+rho0-rhoi)
      !
      REAL(wp) ::   zQm          ! enthalpy exchanged with the ocean (J/m2), >0 towards the ocean
      REAL(wp) ::   zEi          ! specific enthalpy of sea ice (J/kg)
      REAL(wp) ::   zEw          ! specific enthalpy of exchanged water (J/kg)
      REAL(wp) ::   zdE          ! specific enthalpy difference (J/kg)
      REAL(wp) ::   zfmdt        ! exchange mass flux x time step (J/m2), >0 towards the ocean
      REAL(wp) ::   zevap_rema   ! remaining mass flux from sublimation        (kg.m-2)
      REAL(wp) ::   zdeltah, zs_i_new, zds, zs_sni
      REAL(wp) ::   zswitch_sal
      !
      REAL(wp) ::   zq_top      ! heat for surface ablation                   (J.m-2)
      REAL(wp) ::   zq_bot      ! heat for bottom ablation                    (J.m-2)
      REAL(wp) ::   zf_tt       ! Heat budget to determine melting or freezing(W.m-2)
      REAL(wp) ::   zsnw        ! distribution of snow after wind blowing
      !
      INTEGER , DIMENSION(nlay_i)     ::   icount    ! number of layers vanishing by melting
      REAL(wp), DIMENSION(nlay_i)     ::   zs_i      ! ice salinity
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zh_i      ! ice layer thickness (m)
      REAL(wp), DIMENSION(0:nlay_s  ) ::   zh_s      ! snw layer thickness (m)
      REAL(wp), DIMENSION(0:nlay_s  ) ::   ze_s      ! snw layer enthalpy (J.m-3)
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zh_i_o    ! old thickness
      REAL(wp), DIMENSION(0:nlay_i+1) ::   ze_i_o    ! old enthalpy
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zs_i_o    ! old salt content
      !
      ! For "manual inlining" of routine `ice_var_vremap` and `snw_ent`:
      INTEGER  ::   jk0, jk1   !  old/new layer indices
      REAL(wp) ::   zhnew      ! new layers thicknesses
      ! For ice enthalpy and salt content:
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zxi_cum0, zhi_cum0   ! old cumulative enthlapies/salinities and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zxi_cum1, zhi_cum1   ! new cumulative enthlapies/salinities and layers interfaces
      ! For snow enthalpy:
      REAL(wp), DIMENSION(0:nlay_s+1) ::   zes_cum0, zhs_cum0   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(0:nlay_s)   ::   zes_cum1, zhs_cum1   ! new cumulative enthlapies and layers interfaces
      !!------------------------------------------------------------------
      IF( ln_timing    )   CALL timing_start('ice_thd_dh')

      !$acc data present( a_i,at_i,dh_i_bog,dh_i_bom,dh_i_itm,dh_i_sub,dh_i_sum_2d,dh_s_itm,dh_snowice,dh_s_sum_2d,e_i,e_s,evap_ice,fatm_snow,fhld,frq_m )
      !$acc data present( hfx_bog,hfx_bom,hfx_res,hfx_snw,hfx_spr,hfx_sub,hfx_sum,hfx_thd,h_i,h_i,h_s,ll_ice_present,qcn_ice_top )
      !$acc data present( qml_ice,qprec_ice,qsb_ice_bot,qsr_ice,qtr_ice_bot,qtr_ice_top,sfx_bog,sfx_bom,sfx_bri,sfx_res,sfx_sni,sfx_sub,sfx_sum,s_i )
      !$acc data present( sss_s,sst_s,sz_i,t_bo,t_i,t_s,t_su,wfx_bog,wfx_bom,wfx_err_sub,wfx_ice_sub,wfx_res,wfx_sni,wfx_snw_sni,wfx_snw_sub,wfx_snw_sum,wfx_spr,wfx_sum )
      !KEEP: fhld, qsb_ice_bot

      !$acc data create( icount,zs_i,zh_i,zh_s,ze_s,zh_i_o,ze_i_o,zs_i_o, zxi_cum0,zhi_cum0,zxi_cum1,zhi_cum1,zes_cum0,zhs_cum0,zes_cum1,zhs_cum1 )

      ! Discriminate between time varying salinity and constant
      SELECT CASE( nn_icesal )                  ! varying salinity or not
      CASE( 1 , 3 )   ;   zswitch_sal = 0._wp   ! prescribed salinity profile
      CASE( 2 , 4 )   ;   zswitch_sal = 1._wp   ! varying salinity profile
      END SELECT
      !
      ! for snw-ice formation
      z1_rho = 1._wp / ( rhos+rho0-rhoi )


      !                       ! ==================== !
      !                       ! Start main loop here !
      !                       ! ==================== !

      !$acc parallel loop collapse(2) private(icount,zs_i,zh_i,zh_s,ze_s,zh_i_o,ze_i_o,zs_i_o, zxi_cum0,zhi_cum0,zxi_cum1,zhi_cum1,zes_cum0,zhs_cum0,zes_cum1,zhs_cum1)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !
            IF( ll_ice_present(ji,jj) ) THEN
               !                       ! ============================================== !
               !                       ! Available heat for surface and bottom ablation !
               !                       ! ============================================== !
               IF( .NOT.ln_cndflx .OR. ln_cndemulate ) THEN
                  IF( t_su(ji,jj,jl_cat) >= rt0 ) THEN
                     qml_ice(ji,jj,jl_cat) = qns_ice(ji,jj,jl_cat) + qsr_ice(ji,jj,jl_cat) - qtr_ice_top(ji,jj,jl_cat) - qcn_ice_top(ji,jj,jl_cat)
                  ELSE
                     qml_ice(ji,jj,jl_cat) = 0._wp
                  ENDIF
               ENDIF
               !
               zq_top = MAX( 0._wp, qml_ice(ji,jj,jl_cat) * rDt_ice )
               zf_tt  = qcn_ice_bot(ji,jj,jl_cat) + qsb_ice_bot(ji,jj) + fhld(ji,jj) + qtr_ice_bot(ji,jj,jl_cat) * frq_m(ji,jj)
               zq_bot = MAX( 0._wp, zf_tt * rDt_ice )
               !
               ! initialize salinity
               IF( nn_icesal == 4 ) THEN
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     zs_i(jk) = sz_i(ji,jj,jk,jl_cat)  ! use layer salinity
                  END DO
               ELSE
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     zs_i(jk) = s_i(ji,jj,jl_cat)      !     bulk salinity otherwise (for conservation purpose)  !BUG reported to la Rousette! (used to be s_i(ji,jj,:))
                  END DO
               ENDIF
               !
               ! initialize ice layer thicknesses and enthalpies
               !$acc loop seq
               DO jk = 0, nlay_i+1
                  zs_i_o(jk) = 0._wp
                  ze_i_o(jk) = 0._wp
                  zh_i_o(jk) = 0._wp
                  zh_i    (jk) = 0._wp
               END DO
               !$acc loop seq
               DO jk = 1, nlay_i
                  zs_i_o(jk) = h_i(ji,jj,jl_cat) * r1_nlay_i * zs_i  (jk)
                  ze_i_o(jk) = h_i(ji,jj,jl_cat) * r1_nlay_i * e_i(ji,jj,jk,jl_cat)
                  zh_i_o(jk) = h_i(ji,jj,jl_cat) * r1_nlay_i
                  zh_i    (jk) = h_i(ji,jj,jl_cat) * r1_nlay_i
               END DO
               !
               ! initialize snw layer thicknesses and enthalpies
               !$acc loop seq
               DO jk = 0, nlay_s
                  zh_s(jk) = 0._wp
                  ze_s(jk) = 0._wp
               END DO
               !$acc loop seq
               DO jk = 1, nlay_s
                  zh_s(jk) = h_s(ji,jj,   jl_cat) * r1_nlay_s
                  ze_s(jk) = e_s(ji,jj,jk,jl_cat)
               END DO

               !                       ! ============ !
               !                       !     Snow     !
               !                       ! ============ !
               !
               ! Internal melting
               ! ----------------
               ! IF snow temperature is above freezing point, THEN snow melts (should not happen but sometimes it does)
               !$acc loop seq
               DO jk = 1, nlay_s
                  IF( t_s(ji,jj,jk,jl_cat) > rt0 ) THEN
                     hfx_res(ji,jj)     = hfx_res(ji,jj) - ze_s(jk) * zh_s(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! heat flux to the ocean [W.m-2], < 0
                     wfx_snw_sum(ji,jj) = wfx_snw_sum(ji,jj) + rhos * zh_s(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! mass flux
                     ! updates
                     dh_s_itm(ji,jj) =             dh_s_itm(ji,jj) - zh_s(jk)
                     h_s(ji,jj,jl_cat) = MAX( 0._wp, h_s(ji,jj,jl_cat) - zh_s(jk) )
                     zh_s(jk) = 0._wp
                     ze_s(jk) = 0._wp
                  ENDIF
               END DO

               ! Snow precipitation
               !-------------------
               IF( fatm_snow(ji,jj) > 0._wp ) THEN
                  zsnw =   1._wp - MAX(1._wp - at_i(ji,jj),0._wp )**rn_snwblow   ! snow distribution over ice after wind blowing   !LOLO inlining!
                  zh_s(0) = zsnw * fatm_snow(ji,jj) * rDt_ice * r1_rhos / at_i(ji,jj)   ! thickness of precip
                  ze_s(0) = MAX( 0._wp, - qprec_ice(ji,jj) )                              ! enthalpy of the precip (>0, J.m-3)
                  !
                  hfx_spr(ji,jj) = hfx_spr(ji,jj) + ze_s(0) * zh_s(0) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! heat flux from snow precip (>0, W.m-2)
                  wfx_spr(ji,jj) = wfx_spr(ji,jj) - rhos    * zh_s(0) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! mass flux, <0
                  !
                  ! update thickness
                  h_s(ji,jj,jl_cat) = h_s(ji,jj,jl_cat) + zh_s(0)
               ENDIF

               ! Snow sublimation and deposition
               !-----------------
               ! if qla_ice is >=0 (upwards), heat goes to the atmosphere, therefore snow sublimates
               ! else                       , there is snow deposition
               !    comment: not counted in mass/heat exchange in iceupdate.F90 since this is an exchange with atm. (not ocean)
               zdeltah    = MAX( - evap_ice(ji,jj,jl_cat) * r1_rhos * rDt_ice, - h_s(ji,jj,jl_cat) )   ! amount of snw that sublimates (<0) or deposition (>0)
               zevap_rema =        evap_ice(ji,jj,jl_cat)           * rDt_ice + zdeltah * rhos  ! remaining evap in kg.m-2 (used for ice sublimation later on)
               IF( zdeltah > 0._wp .AND. ze_s(0) == 0._wp ) THEN   ! if snow deposition and no snow precip, then estimate ze_s(0) with t_su
                  ze_s(0) = rhos * ( rLfus - rcpi * ( t_su(ji,jj,jl_cat) - rt0 ) )
               ENDIF
               !$acc loop seq
               DO jk = 0, nlay_s
                  zdum = MAX( -zh_s(jk), zdeltah ) ! snow layer thickness that sublimates (<0) or deposits (>0)
                  !
                  hfx_sub(ji,jj) = hfx_sub(ji,jj) + ze_s(jk) * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice  ! Heat flux of snw that sublimates/deposits [W.m-2], <0 or >0
                  wfx_snw_sub(ji,jj) = wfx_snw_sub(ji,jj) - rhos     * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice  ! Mass flux by sublimation or deposition
                  ! update thickness
                  h_s(ji,jj,jl_cat) = MAX( 0._wp , h_s(ji,jj,jl_cat) + zdum )
                  zh_s(jk) = MAX( 0._wp , zh_s(jk) + zdum )
                  ! update sublimation left (if any)
                  zdeltah = MIN( zdeltah - zdum, 0._wp )
               END DO
               !

               ! Snow melting
               ! ------------
               ! If heat still available (zq_top > 0)
               ! then all snw precip has been melted and we need to melt more snow
               !$acc loop seq
               DO jk = 0, nlay_s
                  IF( zh_s(jk) > 0._wp .AND. zq_top > 0._wp ) THEN
                     !
                     zdum = - zq_top / MAX( ze_s(jk), epsi20 )   ! thickness change
                     zdum = MAX( zdum , - zh_s(jk) )                 ! bound melting

                     hfx_snw(ji,jj) = hfx_snw(ji,jj) - ze_s(jk) * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! heat used to melt snow(W.m-2, >0)
                     wfx_snw_sum(ji,jj) = wfx_snw_sum(ji,jj) - rhos     * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! snow melting only = water into the ocean

                     ! updates available heat + thickness
                     dh_s_sum_2d(ji,jj,jl_cat) =              dh_s_sum_2d(ji,jj,jl_cat) + zdum
                     zq_top       = MAX( 0._wp , zq_top       + zdum * ze_s(jk) )
                     h_s(ji,jj,jl_cat) = MAX( 0._wp , h_s(ji,jj,jl_cat) + zdum )
                     zh_s    (jk) = MAX( 0._wp , zh_s    (jk) + zdum )
                     !
                  ENDIF
               END DO

               !
               !                       ! ============ !
               !                       !     Ice      !
               !                       ! ============ !

               ! Surface ice melting
               !--------------------
               !$acc loop seq
               DO jk = 1, nlay_i
                  ztmelts = - rTmlt * sz_i(ji,jj,jk,jl_cat)   ! Melting point of layer k [C]

                  IF( t_i(ji,jj,jk,jl_cat) >= (ztmelts+rt0) ) THEN   !-- Internal melting

                     zEi            = - e_i(ji,jj,jk,jl_cat) * r1_rhoi             ! Specific enthalpy of layer k [J/kg, <0]
                     zdE            =   0._wp                               ! Specific enthalpy difference (J/kg, <0)
                     !                                                          set up at 0 since no energy is needed to melt water...(it is already melted)
                     zdum           = MIN( 0._wp , - zh_i(jk) )             ! internal melting occurs when the internal temperature is above freezing
                     !                                                          this should normally not happen, but sometimes, heat diffusion leads to this
                     zfmdt          = - zdum * rhoi                         ! Recompute mass flux [kg/m2, >0]
                     !
                     dh_i_itm(ji,jj)   = dh_i_itm(ji,jj) + zdum                   ! Cumulate internal melting
                     !
                     hfx_res(ji,jj) = hfx_res(ji,jj) + zEi  * zfmdt           * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Heat flux to the ocean [W.m-2], <0
                     !                                                                                        ice enthalpy zEi is "sent" to the ocean
                     wfx_res(ji,jj) = wfx_res(ji,jj) - rhoi * zdum            * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Mass flux
                     sfx_res(ji,jj) = sfx_res(ji,jj) - rhoi * zdum * zs_i(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Salt flux
                     !
                  ELSE                                        !-- Surface melting

                     zEi            = - e_i(ji,jj,jk,jl_cat) * r1_rhoi             ! Specific enthalpy of layer k [J/kg, <0]
                     zEw            =    rcp * ztmelts                      ! Specific enthalpy of resulting meltwater [J/kg, <0]
                     zdE            =    zEi - zEw                          ! Specific enthalpy difference < 0

                     zfmdt          = - zq_top / zdE                        ! Mass flux to the ocean [kg/m2, >0]

                     zdum           = - zfmdt * r1_rhoi                     ! Melt of layer jk [m, <0]

                     zdum           = MIN( 0._wp , MAX( zdum , - zh_i(jk) ) )    ! Melt of layer jk cannot exceed the layer thickness [m, <0]

                     zq_top         = MAX( 0._wp , zq_top - zdum * rhoi * zdE ) ! update available heat

                     dh_i_sum_2d(ji,jj,jl_cat)   = dh_i_sum_2d(ji,jj,jl_cat) + zdum                   ! Cumulate surface melt

                     zfmdt          = - rhoi * zdum                         ! Recompute mass flux [kg/m2, >0]

                     zQm            = zfmdt * zEw                           ! Energy of the melt water sent to the ocean [J/m2, <0]

                     hfx_thd(ji,jj) = hfx_thd(ji,jj) + zEw  * zfmdt           * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Heat flux [W.m-2], < 0
                     hfx_sum(ji,jj) = hfx_sum(ji,jj) - zdE  * zfmdt           * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Heat flux used in this process [W.m-2], > 0
                     wfx_sum(ji,jj) = wfx_sum(ji,jj) - rhoi * zdum            * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Mass flux
                     sfx_sum(ji,jj) = sfx_sum(ji,jj) - rhoi * zdum * zs_i(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice    ! Salt flux >0
                     !
                  ENDIF
                  ! update thickness
                  zh_i  (jk) = MAX( 0._wp, zh_i  (jk) + zdum )
                  h_i(ji,jj,jl_cat) = MAX( 0._wp, h_i(ji,jj,jl_cat) + zdum )
                  !
                  ! update heat content (J.m-2), salt content and layer thickness
                  zs_i_o(jk) = zs_i_o(jk) + zdum * zs_i(jk)
                  ze_i_o(jk) = ze_i_o(jk) + zdum * e_i(ji,jj,jk,jl_cat)
                  zh_i_o(jk) = zh_i_o(jk) + zdum
                  !
                  !
                  ! Ice sublimation
                  ! ---------------
                  zdum               = MAX( - zh_i(jk) , - zevap_rema * r1_rhoi )
                  !
                  hfx_sub(ji,jj)     = hfx_sub(ji,jj) + e_i(ji,jj,jk,jl_cat) * zdum  * a_i(ji,jj,jl_cat) * r1_Dt_ice ! Heat flux [W.m-2], < 0
                  wfx_ice_sub(ji,jj) = wfx_ice_sub(ji,jj) - rhoi * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice ! Mass flux > 0
                  sfx_sub(ji,jj)     = sfx_sub(ji,jj)     - rhoi * zdum * zs_i(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice ! Salt flux >0
                  !                                                             clem: flux is sent to the ocean for simplicity
                  !                                                                   but salt should remain in the ice except
                  !                                                                   if all ice is melted. => must be corrected
                  ! update remaining mass flux and thickness
                  zevap_rema   = zevap_rema + zdum * rhoi
                  zh_i  (jk)   = MAX( 0._wp, zh_i  (jk) + zdum )
                  h_i(ji,jj,jl_cat)   = MAX( 0._wp, h_i(ji,jj,jl_cat) + zdum )
                  dh_i_sub(ji,jj) = dh_i_sub(ji,jj) + zdum

                  ! update heat content (J.m-2), salt content and layer thickness
                  zs_i_o(jk) = zs_i_o(jk) + zdum * zs_i(jk)
                  ze_i_o(jk) = ze_i_o(jk) + zdum * e_i(ji,jj,jk,jl_cat)
                  zh_i_o(jk) = zh_i_o(jk) + zdum

                  ! record which layers have disappeared (for bottom melting)
                  !    => icount=0 : no layer has vanished
                  !    => icount=5 : 5 layers have vanished
                  IF( zh_i(jk) > 0._wp ) THEN
                     icount(jk) = 0
                  ELSE
                     icount(jk) = 1
                  ENDIF

               END DO

               ! remaining "potential" evap is sent to ocean
               wfx_err_sub(ji,jj) = wfx_err_sub(ji,jj) - zevap_rema * a_i(ji,jj,jl_cat) * r1_Dt_ice  ! <=0 (net evap for the ocean in kg.m-2.s-1)


               ! Ice Basal growth
               !------------------
               ! Basal growth is driven by heat imbalance at the ice-ocean interface,
               ! between the inner conductive flux  (qcn_ice_bot), from the open water heat flux
               ! (fhld) and the sensible ice-ocean flux (qsb_ice_bot).
               ! qcn_ice_bot is positive downwards. qsb_ice_bot and fhld are positive to the ice
               !
               zs_i_new = 0._wp
               !
               IF(  zf_tt < 0._wp  ) THEN

                  zs_i_new       = zswitch_sal * rn_sinew * sss_s(ji,jj) + ( 1. - zswitch_sal ) * zs_i(1)          ! New ice salinity

                  ztmelts        = - rTmlt * zs_i_new                                                            ! New ice melting point (C)

                  zt_i_new       = zswitch_sal * t_bo(ji,jj) + ( 1. - zswitch_sal) * t_i(ji,jj, nlay_i,jl_cat)

                  zEi            = rcpi * ( zt_i_new - (ztmelts+rt0) ) &                                         ! Specific enthalpy of forming ice (J/kg, <0)
                     &             - rLfus * ( 1.0 - ztmelts / ( MIN( zt_i_new - rt0, -epsi10 ) ) ) + rcp * ztmelts

                  zEw            = rcp  * ( t_bo(ji,jj) - rt0 )                                                  ! Specific enthalpy of seawater (J/kg, < 0)

                  zdE            = zEi - zEw                                                                     ! Specific enthalpy difference (J/kg, <0)

                  dh_i_bog(ji,jj)   = rDt_ice * MAX( 0._wp , zf_tt / ( zdE * rhoi ) )

                  ! Contribution to Energy and Salt Fluxes
                  zfmdt = - rhoi * dh_i_bog(ji,jj)                                                                  ! Mass flux x time step (kg/m2, < 0)

                  hfx_thd(ji,jj) = hfx_thd(ji,jj) + zEw  * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], >0
                  hfx_bog(ji,jj) = hfx_bog(ji,jj) - zdE  * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Heat flux used in this process [W.m-2], <0
                  wfx_bog(ji,jj) = wfx_bog(ji,jj) - rhoi * dh_i_bog(ji,jj) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Mass flux, <0
                  sfx_bog(ji,jj) = sfx_bog(ji,jj) - rhoi * dh_i_bog(ji,jj) * zs_i_new * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Salt flux, <0

                  ! update thickness
                  zh_i(nlay_i+1) = zh_i(nlay_i+1) + dh_i_bog(ji,jj)
                  h_i(ji,jj,jl_cat)     = h_i(ji,jj,jl_cat) + dh_i_bog(ji,jj)

                  ! update heat content (J.m-2), salt content and layer thickness
                  zs_i_o(nlay_i+1) = zs_i_o(nlay_i+1) + dh_i_bog(ji,jj) * zs_i_new
                  ze_i_o(nlay_i+1) = ze_i_o(nlay_i+1) + dh_i_bog(ji,jj) * (-zEi * rhoi)
                  zh_i_o(nlay_i+1) = zh_i_o(nlay_i+1) + dh_i_bog(ji,jj)

               ENDIF

               ! Ice Basal melt
               !---------------
               !$acc loop seq
               DO jk = nlay_i, 1, -1
                  IF(  zf_tt  >  0._wp  .AND. jk > icount(jk) ) THEN   ! do not calculate where layer has already disappeared by surface melting

                     ztmelts = - rTmlt * sz_i(ji,jj,jk,jl_cat)  ! Melting point of layer jk (C)

                     IF( t_i(ji,jj,jk,jl_cat) >= (ztmelts+rt0) ) THEN   !-- Internal melting

                        zEi            = - e_i(ji,jj,jk,jl_cat) * r1_rhoi     ! Specific enthalpy of melting ice (J/kg, <0)
                        zdE            = 0._wp                         ! Specific enthalpy difference   (J/kg, <0)
                        !                                                  set up at 0 since no energy is needed to melt water...(it is already melted)
                        zdum           = MIN( 0._wp , - zh_i(jk) )  ! internal melting occurs when the internal temperature is above freezing
                        !                                                  this should normally not happen, but sometimes, heat diffusion leads to this
                        dh_i_itm (ji,jj)  = dh_i_itm(ji,jj) + zdum
                        !
                        zfmdt          = - zdum * rhoi                 ! Mass flux x time step > 0
                        !
                        hfx_res(ji,jj) = hfx_res(ji,jj) + zEi  * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], <0
                        !                                                                                       ice enthalpy zEi is "sent" to the ocean
                        wfx_res(ji,jj) = wfx_res(ji,jj) - rhoi * zdum * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Mass flux
                        sfx_res(ji,jj) = sfx_res(ji,jj) - rhoi * zdum * zs_i(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Salt flux
                        !
                     ELSE                                        !-- Basal melting

                        zEi            = - e_i(ji,jj,jk,jl_cat) * r1_rhoi                ! Specific enthalpy of melting ice (J/kg, <0)
                        zEw            = rcp * ztmelts                                   ! Specific enthalpy of meltwater (J/kg, <0)
                        zdE            = zEi - zEw                                       ! Specific enthalpy difference   (J/kg, <0)

                        zfmdt          = - zq_bot / zdE                                  ! Mass flux x time step (kg/m2, >0)

                        zdum           = - zfmdt * r1_rhoi                               ! Gross thickness change

                        zdum           = MIN( 0._wp , MAX( zdum, - zh_i(jk) ) )          ! bound thickness change

                        zq_bot         = MAX( 0._wp , zq_bot - zdum * rhoi * zdE )       ! update available heat. MAX is necessary for roundup errors

                        dh_i_bom(ji,jj)   = dh_i_bom(ji,jj) + zdum                             ! Update basal melt

                        zfmdt          = - zdum * rhoi                                   ! Mass flux x time step > 0

                        zQm            = zfmdt * zEw                                     ! Heat exchanged with ocean

                        hfx_thd(ji,jj) = hfx_thd(ji,jj) + zEw  * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Heat flux to the ocean [W.m-2], <0
                        hfx_bom(ji,jj) = hfx_bom(ji,jj) - zdE  * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Heat used in this process [W.m-2], >0
                        wfx_bom(ji,jj) = wfx_bom(ji,jj) - rhoi * zdum  * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Mass flux
                        sfx_bom(ji,jj) = sfx_bom(ji,jj) - rhoi * zdum * zs_i(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice   ! Salt flux
                        !
                     ENDIF
                     ! update thickness
                     zh_i  (jk) = MAX( 0._wp, zh_i  (jk) + zdum )
                     h_i(ji,jj,jl_cat) = MAX( 0._wp, h_i(ji,jj,jl_cat) + zdum )
                     !
                     ! update heat content (J.m-2), salt content and layer thickness
                     zs_i_o(jk) = zs_i_o(jk) + zdum * zs_i(jk)
                     ze_i_o(jk) = ze_i_o(jk) + zdum * e_i(ji,jj,jk,jl_cat)
                     zh_i_o(jk) = zh_i_o(jk) + zdum
                  ENDIF
               END DO

               ! Remove snow if ice has melted entirely
               ! --------------------------------------
               IF( h_i(ji,jj,jl_cat) == 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 0, nlay_s
                     ! mass & energy loss to the ocean
                     hfx_res(ji,jj) = hfx_res(ji,jj) - ze_s(jk) * zh_s(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice  ! heat flux to the ocean [W.m-2], < 0
                     wfx_res(ji,jj) = wfx_res(ji,jj) + rhos     * zh_s(jk) * a_i(ji,jj,jl_cat) * r1_Dt_ice  ! mass flux

                     ! update thickness and energy
                     h_s(ji,jj,jl_cat) = 0._wp
                     ze_s(jk) = 0._wp
                     zh_s(jk) = 0._wp
                  END DO
               ENDIF


               ! Snow-Ice formation
               ! ------------------
               ! When snow load exceeds Archimede's limit, snow-ice interface goes down under sea-level,
               ! flooding of seawater transforms snow into ice. Thickness that is transformed is dh_snowice (positive for the ice)
               !
               dh_snowice(ji,jj) = MAX( 0._wp , ( rhos * h_s(ji,jj,jl_cat) + (rhoi-rho0) * h_i(ji,jj,jl_cat) ) * z1_rho )

               h_i(ji,jj,jl_cat)    = h_i(ji,jj,jl_cat) + dh_snowice(ji,jj)
               h_s(ji,jj,jl_cat)    = h_s(ji,jj,jl_cat) - dh_snowice(ji,jj)

               ! Contribution to energy flux to the ocean [J/m2], >0 (if sst<0)
               zfmdt          = ( rhos - rhoi ) * dh_snowice(ji,jj)    ! <0
               zEw            = rcp * sst_s(ji,jj)
               zQm            = zfmdt * zEw

               hfx_thd(ji,jj) = hfx_thd(ji,jj) + zEw * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice ! Heat flux
               sfx_sni(ji,jj) = sfx_sni(ji,jj) + sss_s(ji,jj) * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice ! Salt flux

               ! Case constant salinity in time: virtual salt flux to keep salinity constant
               IF( nn_icesal == 1 .OR. nn_icesal == 3 )  THEN
                  sfx_bri(ji,jj) = sfx_bri(ji,jj) - sss_s(ji,jj) * zfmdt * a_i(ji,jj,jl_cat) * r1_Dt_ice  &  ! put back sss_s     into the ocean
                     & - zs_i(1) * dh_snowice(ji,jj) * rhoi * a_i(ji,jj,jl_cat) * r1_Dt_ice     ! and get  rn_icesal from the ocean
               ENDIF

               ! Mass flux: All snow is thrown in the ocean, and seawater is taken to replace the volume
               wfx_sni(ji,jj) = wfx_sni(ji,jj) - dh_snowice(ji,jj) * rhoi * a_i(ji,jj,jl_cat) * r1_Dt_ice
               wfx_snw_sni(ji,jj) = wfx_snw_sni(ji,jj) + dh_snowice(ji,jj) * rhos * a_i(ji,jj,jl_cat) * r1_Dt_ice

               ! update thickness
               zh_i(0) = zh_i(0) + dh_snowice(ji,jj)
               zdeltah =           dh_snowice(ji,jj)

               ! update heat content (J.m-2), salt content and layer thickness
               zs_i_o(0) = zs_i_o(0) - zfmdt * sss_s(ji,jj) * r1_rhoi      ! clem: s(0) could be > rn_sinew*sss
               zh_i_o(0) = zh_i_o(0) + dh_snowice(ji,jj)
               ze_i_o(0) = ze_i_o(0) + zfmdt * zEw           ! 1st part (sea water enthalpy)

               !$acc loop seq
               DO jk = nlay_s, 0, -1   ! flooding of snow starts from the base
                  zdum        = MIN( zdeltah, zh_s(jk) )         ! amount of snw that floods, > 0
                  zh_s(jk)    = MAX( 0._wp, zh_s(jk) - zdum )    ! remove some snow thickness
                  ze_i_o(0) = ze_i_o(0) + zdum * ze_s(jk)    ! 2nd part (snow enthalpy)
                  ! update dh_snowice
                  zdeltah     = MAX( 0._wp, zdeltah - zdum )
               END DO
               !
               !
               ! Remapping of snw enthalpy on a regular grid
               !--------------------------------------------
               !e_s(ji,jj,:,jl_cat) = snw_ent( zh_s, ze_s )   !lili: manual inlining instead of call to `snw_ent()`:
               !  ==> inlining, better for GPU...
# include      "ice_var_snw_ent.h90"
               !
               ! recalculate t_s from e_s
               IF( h_s(ji,jj,jl_cat) > 0._wp ) THEN
                  !$acc loop seq
                  DO jk = 1, nlay_s
                     t_s(ji,jj,jk,jl_cat) = rt0 + ( - e_s(ji,jj,jk,jl_cat) * r1_rhos * r1_rcpi + rLfus * r1_rcpi )
                  END DO
               ELSE
                  !$acc loop seq
                  DO jk = 1, nlay_s
                     t_s(ji,jj,jk,jl_cat) = rt0
                  END DO
               ENDIF


               ! Remapping of ice enthalpy on a regular grid
               !--------------------------------------------
               !CALL ice_var_vremap( zh_i_o, ze_i_o, e_i(ji,jj,:,jl_cat) ) ! manual inlining of `ice_var_vremap`:
               !  ==> inlining, better for GPU...
               jl = jl_cat
# include      "ice_var_vremap_e.h90"

               ! Remapping of ice salt on a regular grid
               !----------------------------------------
               IF( nn_icesal == 4 ) THEN
                  !CALL ice_var_vremap( zh_i_o, zs_i_o, sz_i(ji,jj,:,jl_cat) ) ! manual inlining of `ice_var_vremap`:
                  !  ==> inlining, better for GPU...
# include         "ice_var_vremap_s.h90"
                  !$acc loop seq
                  DO jk1 = 1, nlay_i
                     sz_i(ji,jj,jk1,jl_cat) = MAX( 0._wp, zxi_cum1(jk1) - zxi_cum1(jk1-1) ) * zdum ! max for roundoff error
                  END DO
               ENDIF
               !
               IF( nn_icesal == 2 )   THEN ! Update ice salinity from snow-ice and bottom growth
                  zs_sni = sss_s(ji,jj) * ( rhoi - rhos ) * r1_rhoi                                       ! salinity of snow ice
                  zds    =       ( zs_sni   - s_i(ji,jj,jl_cat) ) * dh_snowice(ji,jj) / MAX( epsi10, h_i(ji,jj,jl_cat) ) ! snow-ice
                  zds    = zds + ( zs_i_new - s_i(ji,jj,jl_cat) ) * dh_i_bog(ji,jj) / MAX( epsi10, h_i(ji,jj,jl_cat) ) ! bottom growth
                  !
                  s_i(ji,jj,jl_cat) = s_i(ji,jj,jl_cat) + zds
               ENDIF


            ENDIF ! ll_ice_present
            !
         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !                       ! ================== !
      !                       ! End main loop here !
      !                       ! ================== !


      ! --- ensure that a_i = 0 & h_s = 0 where h_i = 0 ---
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF( ll_ice_present(ji,jj) .AND. h_i(ji,jj,jl_cat) <= 0._wp) THEN
               a_i(ji,jj,jl_cat)  = 0._wp
               h_i(ji,jj,jl_cat)  = 0._wp
               h_s(ji,jj,jl_cat)  = 0._wp
               t_su(ji,jj,jl_cat) = 0._wp
               ll_ice_present(ji,jj) = .FALSE.
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      !$acc end data
      !$acc end data
      !$acc end data
      !$acc end data
      IF( ln_timing    )   CALL timing_stop('ice_thd_dh')

   END SUBROUTINE ice_thd_dh

   !!======================================================================
END MODULE icethd_dh
