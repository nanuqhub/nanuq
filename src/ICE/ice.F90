MODULE ice
   !!======================================================================
   !!                        ***  MODULE  ice  ***
   !!   sea-ice:  ice variables defined in memory
   !!======================================================================
   !! History :  3.0  !  2008-03  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_alloc   ! called by icestp.F90

   !!======================================================================
   !!                                                                     |
   !!              I C E   S T A T E   V A R I A B L E S                  |
   !!                                                                     |
   !! Introduction :                                                      |
   !! --------------                                                      |
   !! Every ice-covered grid cell is characterized by a series of state   |
   !! variables. To account for unresolved spatial variability in ice     |
   !! thickness, the ice cover in divided in ice thickness categories.    |
   !!                                                                     |
   !! Sea ice state variables depend on the ice thickness category        |
   !!                                                                     |
   !! Those variables are divided into two groups                         |
   !! * Extensive (or global) variables.                                  |
   !!   These are the variables that are transported by all means         |
   !! * Intensive (or equivalent) variables.                              |
   !!   These are the variables that are either physically more           |
   !!   meaningful and/or used in ice thermodynamics                      |
   !!                                                                     |
   !! List of ice state variables :                                       |
   !! -----------------------------                                       |
   !!                                                                     |
   !!-------------|-------------|---------------------------------|-------|
   !!   name in   |   name in   |              meaning            | units |
   !! 2D routines | 1D routines |                                 |       |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Dynamical variables (prognostic)                    *** |
   !! ******************************************************************* |
   !!                                                                     |
   !! u_ice       |      -      |    ice velocity in i-direction  | m/s   |
   !! v_ice       |      -      |    ice velocity in j-direction  | m/s   |
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Category dependent state variables (prognostic)     *** |
   !! ******************************************************************* |
   !!                                                                     |
   !! ** Global variables                                                 |
   !!-------------|-------------|---------------------------------|-------|
   !! a_i         |   a_i_1d    |    Ice concentration            |       |
   !! v_i         |      -      |    Ice volume per unit area     | m     |
   !! v_s         |      -      |    Snow volume per unit area    | m     |
   !! sv_i        |      -      |    Sea ice salt content (3D)    | g/kg.m|
   !! szv_i       |      -      |    Sea ice salt content (4D)    | g/kg.m|
   !! oa_i        |      -      |    Sea ice areal age content    | s     |
   !! e_i         |             |    Ice enthalpy                 | J/m2  |
   !!             |    e_i_1d   |    Ice enthalpy per unit vol.   | J/m3  |
   !! e_s         |             |    Snow enthalpy                | J/m2  |
   !!             |    e_s_1d   |    Snow enthalpy per unit vol.  | J/m3  |
   !! a_ip        |      -      |    Ice pond concentration       |       |
   !! v_ip        |      -      |    Ice pond volume per unit area| m     |
   !! v_il        |    v_il_1d  |    Ice pond lid volume per area | m     |
   !!                                                                     |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! ** Equivalent variables                                             |
   !!-------------|-------------|---------------------------------|-------|
   !!                                                                     |
   !! h_i         | h_i_1d      |    Ice thickness                | m     |
   !! h_s         ! h_s_1d      |    Snow depth                   | m     |
   !! s_i         ! s_i_1d      |    Sea ice bulk salinity        | g/kg  |
   !! sz_i        ! sz_i_1d     |    Sea ice salinity profile     | g/kg  |
   !! o_i         !      -      |    Sea ice Age                  | s     |
   !! dmdt        !      -      |    1-Sea ice Damage @T            !       !
   !! dmdf        !      -      |    1-Sea ice Damage @F            !       !
   !! t_i         ! t_i_1d      |    Sea ice temperature          | K     |
   !! t_s         ! t_s_1d      |    Snow temperature             | K     |
   !! t_su        ! t_su_1d     |    Sea ice surface temperature  | K     |
   !! h_ip        | h_ip_1d     |    Ice pond thickness           | m     |
   !! h_il        | h_il_1d     |    Ice pond lid thickness       | m     |
   !!                                                                     |
   !! notes: the ice model only sees a bulk (i.e., vertically averaged)   |
   !!        salinity, except in thermodynamic computations, for which    |
   !!        the salinity profile is computed as a function of bulk       |
   !!        salinity                                                     |
   !!                                                                     |
   !!        the sea ice surface temperature is not associated to any     |
   !!        heat content. Therefore, it is not a state variable and      |
   !!        does not have to be advected. Nevertheless, it has to be     |
   !!        computed to determine whether the ice is melting or not      |
   !!                                                                     |
   !! ******************************************************************* |
   !! ***         Category-summed state variables (diagnostic)        *** |
   !! ******************************************************************* |
   !! at_i        | at_i_1d     |    Total ice concentration      |       |
   !! vt_i        |      -      |    Total ice vol. per unit area | m     |
   !! vt_s        |      -      |    Total snow vol. per unit ar. | m     |
   !! st_i        |      -      |    Total Sea ice salt content   | g/kg.m|
   !! sm_i        |      -      |    Mean sea ice salinity        | g/kg  |
   !! tm_i        |      -      |    Mean sea ice temperature     | K     |
   !! tm_s        |      -      |    Mean snow    temperature     | K     |
   !! et_i        |      -      |    Total ice enthalpy           | J/m2  |
   !! et_s        |      -      |    Total snow enthalpy          | J/m2  |
   !! v_ibr       |      -      |    relative brine volume        | ???   |
   !! at_ip       |      -      |    Total ice pond concentration |       |
   !! at_ip_eff   !      -      !    Effective pond concentration |       |
   !! hm_ip       |      -      |    Mean ice pond depth          | m     |
   !! vt_ip       |      -      |    Total ice pond vol. per unit area| m |
   !! hm_il       |      -      |    Mean ice pond lid depth      | m     |
   !! vt_il       |      -      |    Total ice pond lid vol. per area | m |
   !!=====================================================================

   !!----------------------------------------------------------------------
   !! * Share Module variables
   !!----------------------------------------------------------------------

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: weno_lw_t_x, weno_lw_t_y ! WENOX linear  weights for curvilinear T-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: weno_ow_t_x, weno_ow_t_y ! WENOX optimal weights for curvilinear T-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: weno_lw_f_x, weno_lw_f_y ! WENOX linear  weights for curvilinear F-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: weno_ow_f_x, weno_ow_f_y ! WENOX optimal weights for curvilinear F-grid

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: weno_s_lw_t_x, weno_s_lw_t_y ! WENOX linear  weights for curvilinear T-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: weno_s_ow_t_x, weno_s_ow_t_y ! WENOX optimal weights for curvilinear T-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: weno_s_lw_f_x, weno_s_lw_f_y ! WENOX linear  weights for curvilinear F-grid
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: weno_s_ow_f_x, weno_s_ow_f_y ! WENOX optimal weights for curvilinear F-grid

   !                                     !!** define arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   V_oce           !: surface ocean velocity used in ice dynamics
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ht_i_new        !: ice collection thickness accreted in leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fraz_frac       !: fraction of frazil ice accreted at the ice bottom
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   strength        !: ice strength
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   stress1_i, stress2_i, stress12_i   !: 1st, 2nd & diagonal stress tensor element
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   delta_i         !: ice rheology elta factor (Flato & Hibler 95) [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   divu_i          !: Divergence of the velocity field             [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   shear_i         !: Shear of the velocity field                  [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rdg_conv
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   t_bo            !: Sea-Ice bottom temperature [Kelvin]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qlead           !: heat balance of the lead (or of the open ocean)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qsb_ice_bot     !: net downward heat flux from the ice to the ocean
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fhld            !: heat flux from the lead used for bottom melting

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_snw         !: mass flux from snow-ocean mass exchange             [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_snw_sni     !: mass flux from snow ice growth component of wfx_snw [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_snw_sum     !: mass flux from surface melt component of wfx_snw    [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_pnd         !: mass flux from melt pond-ocean mass exchange        [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_spr         !: mass flux from snow precipitation on ice            [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_sub         !: mass flux from sublimation of snow/ice              [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_snw_sub     !: mass flux from snow sublimation                     [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_ice_sub     !: mass flux from ice sublimation                      [kg.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_snw_dyn     !: mass flux from dynamical component of wfx_snw       [kg.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_ice         !: mass flux from ice-ocean mass exchange                   [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_sni         !: mass flux from snow ice growth component of wfx_ice      [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_opw         !: mass flux from lateral ice growth component of wfx_ice   [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_bog         !: mass flux from bottom ice growth component of wfx_ice    [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_dyn         !: mass flux from dynamical ice growth component of wfx_ice [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_bom         !: mass flux from bottom melt component of wfx_ice          [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_sum         !: mass flux from surface melt component of wfx_ice         [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_lam         !: mass flux from lateral melt component of wfx_ice         [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_res         !: mass flux from residual component of wfx_ice             [kg.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   wfx_err_sub     !: mass flux error after sublimation                        [kg.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_bog         !: salt flux due to ice bottom growth                   [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_bom         !: salt flux due to ice bottom melt                     [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_lam         !: salt flux due to ice lateral melt                    [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_sum         !: salt flux due to ice surface melt                    [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_sni         !: salt flux due to snow-ice growth                     [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_opw         !: salt flux due to growth in open water                [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_bri         !: salt flux due to brine rejection                     [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_dyn         !: salt flux due to porous ridged ice formation         [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_res         !: salt flux due to correction on ice thick. (residual) [g/kg.kg.m-2.s-1 => g.m-2.s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_sub         !: salt flux due to ice sublimation                     [g/kg.kg.m-2.s-1 => g.m-2.s-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_bog         !: total heat flux causing bottom ice growth           [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_bom         !: total heat flux causing bottom ice melt             [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_sum         !: total heat flux causing surface ice melt            [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_opw         !: total heat flux causing open water ice formation    [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_dif         !: total heat flux causing Temp change in the ice      [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_snw         !: heat flux for snow melt                             [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_err_dif     !: heat flux remaining due to change in non-solar flux [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qt_atm_oi       !: heat flux at the interface atm-[oce+ice]            [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qt_oce_ai       !: heat flux at the interface oce-[atm+ice]            [W.m-2]

   ! heat flux associated with ice-atmosphere mass exchange
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_sub         !: heat flux for sublimation            [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_spr         !: heat flux of the snow precipitation  [W.m-2]

   ! heat flux associated with ice-ocean mass exchange
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_thd         !: ice-ocean heat flux from thermo processes (icethd_dh) [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_dyn         !: ice-ocean heat flux from ridging                      [W.m-2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfx_res         !: heat flux due to correction on ice thick. (residual)  [W.m-2]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qtr_ice_bot     !: transmitted solar radiation under ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   t1_ice          !: temperature of the first layer          (ln_cndflx=T) [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   cnd_ice         !: effective conductivity of the 1st layer (ln_cndflx=T) [W.m-2.K-1]

   !!----------------------------------------------------------------------
   !! * Ice global state variables
   !!----------------------------------------------------------------------
   !! Variables defined for each ice category
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   h_i           !: Ice thickness                           (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_i           !: Ice fractional areas (concentration)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_i           !: Ice volume per unit area                (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_s           !: Snow volume per unit area               (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   h_s           !: Snow thickness                          (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   t_su          !: Sea-Ice Surface Temperature             (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   s_i           !: Sea-Ice Bulk salinity                   (g/kg)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sv_i          !: Sea-Ice Bulk salinity * volume per area (g/kg.m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   o_i           !: Sea-Ice Age                             (s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   oa_i          !: Sea-Ice Age times ice area              (s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_ibr         !: brine volume

   !! Variables summed over all categories, or associated to all the ice in a single grid cell
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   u_ice, v_ice  !: components of the ice velocity                          (m/s)
   !#bbm:
   INTEGER(1), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   kmsk_ice_t, kmsk_ice_f, kmsk_ice_u, kmsk_ice_v !: mask for presence of ice
   REAL(wp),   PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   uVice, vUice  !: components of the ice velocity in F-centric formalism   (m/s)
   REAL(wp),   PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   V_ts ! all the "time-splitting" sea-ice velocities u,v@U,V & u,v@V,U
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: dmdt, dmdf    !: `1-d` 1 - ice damage / #bbm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: SIGMAt !: T-centric vertically-integrated internal stress tensor => (s11@t, s22@t, s12@f) (N/m2*m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: SIGMAf !: F-centric vertically-integrated internal stress tensor => (s11@f, s22@f, s12@t) (N/m2*m)
   !#bbm.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   vt_i , vt_s   !: ice and snow total volume per unit area                 (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   st_i          !: Total ice salinity content                              (g/kg.m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   at_i          !: ice total fractional area (ice concentration)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   af_i          !: ice total fractional area @F (interpolated)             !!clembbm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   au_i          !: ice total fractional area @U (interpolated)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   av_i          !: ice total fractional area @V (interpolated)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   ato_i         !: =1-at_i ; total open water fractional area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   et_i , et_s   !: ice and snow total heat content                         (J/m2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   tm_i          !: mean ice temperature over all categories                (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   tm_s          !: mean snw temperature over all categories                (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   vm_ibr        !: brine volume averaged over all categories
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sm_i          !: mean sea ice salinity averaged over all categories      (g/kg)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   tm_su         !: mean surface temperature over all categories            (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hm_i          !: mean ice  thickness over all categories                 (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hm_s          !: mean snow thickness over all categories                 (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hm_i_f        !: mean ice  thickness @ F                                 (m)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   om_i          !: mean ice age over all categories                        (s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   tau_icebfr    !: ice friction on ocean bottom (landfast param activated)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   icb_mask      !: mask of grounded icebergs if landfast [0-1]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   t_s           !: Snow temperatures     [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_s           !: Snow enthalpy         [J/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   t_i           !: ice temperatures      [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_i           !: ice enthalpy          [J/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sz_i          !: ice salinity          [g/kg]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   szv_i         !: ice salinity content  [g/kg]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_ip          !: melt pond concentration
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_ip          !: melt pond volume per grid cell area      [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_ip_frac     !: melt pond fraction (a_ip/a_i)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_ip_eff      !: melt pond effective fraction (not covered up by lid) (a_ip/a_i)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   h_ip          !: melt pond depth                          [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_il          !: melt pond lid volume                     [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   h_il          !: melt pond lid thickness                  [m]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   at_ip         !: total melt pond concentration
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hm_ip         !: mean melt pond depth                     [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   vt_ip         !: total melt pond volume per gridcell area [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   hm_il         !: mean melt pond lid depth                     [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   vt_il         !: total melt pond lid volume per gridcell area [m]

   !!----------------------------------------------------------------------
   !! * Global variables at before time step
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_s_b, v_i_b, h_s_b, h_i_b !: snow and ice volumes/thickness
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   v_ip_b, v_il_b             !: ponds and lids volumes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_i_b, sv_i_b              !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_s_b                      !: snow heat content
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   e_i_b, szv_i_b             !: ice temperatures and salt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   u_ice_b, v_ice_b           !: ice velocity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   at_i_b                     !: ice concentration (total)

   !!----------------------------------------------------------------------
   !! * Ice thickness distribution variables
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   hi_max            !: Boundary of ice thickness categories in thickness space
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   hi_mean           !: Mean ice thickness in catgories
   !
   !!----------------------------------------------------------------------
   !! * Ice diagnostics
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_trp_vi       !: transport of ice volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_trp_vs       !: transport of snw volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_trp_ei       !: transport of ice enthalpy [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_trp_es       !: transport of snw enthalpy [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_trp_sv       !: transport of salt content
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_heat         !: snw/ice heat content variation   [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_sice         !: ice salt content variation   []
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_vice         !: ice volume variation   [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_vsnw         !: snw volume variation   [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_aice         !: ice conc.  variation   [s-1]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_vpnd         !: pond volume variation  [m/s]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_adv_mass     !: advection of mass (kg/m2/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_adv_salt     !: advection of salt (g/m2/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_adv_heat     !: advection of heat (W/m2)
   !
   !!----------------------------------------------------------------------
   !! * Ice conservation
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_v            !: conservation of ice volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_s            !: conservation of ice salt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_t            !: conservation of ice heat
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_fv           !: conservation of ice volume
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_fs           !: conservation of ice salt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   diag_ft           !: conservation of ice heat
   !
   !!----------------------------------------------------------------------
   !! * Ice salinity checks in drainage and flushing
   !! * Note these arrays are allocated in icethd.F90 (icethd_salchk), depending on whether ln_sal_chk is enabled
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   zsneg_drain , zsneg_flush
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   zcfl_drain  , zcfl_flush
   !! * For convergence tests; ztice_cvgerr and ztice_cvgstp are allocated in icethd.F90 depending on whether ln_zdf_chkcvg is enabled
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ztice_cvgerr, ztice_cvgstp
   !
   !!----------------------------------------------------------------------
   !! * SIMIP extra diagnostics
   !!----------------------------------------------------------------------
   ! Extra sea ice diagnostics to address the data request
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   t_si            !: Temperature at Snow-ice interface (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   tm_si           !: mean temperature at the snow-ice interface (K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qcn_ice_bot     !: Bottom  conduction flux (W/m2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qcn_ice_top     !: Surface conduction flux (W/m2)
   !
   !#thd2d:
   !!----------------------------------------------------------------------
   !! * 2D version of a few arrays thar are (used to be) declared and
   !!   allocated in `ice1d.F90` for 1d thermo...
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: s_i_new, dh_i_itm, dh_s_itm, dh_s_tot, dh_i_bom, &
      &                                                   dh_i_sub, dh_i_bog, dh_snowice
   ! meltwater arrays to save for melt ponds (mv - could be grouped in a single meltwater volume array)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dh_i_sum_2d, dh_s_sum_2d

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: eh_i_old    !: ice heat content (q*h, J.m-2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: h_i_old     !: ice thickness layer (m)

   !! Transports for advection:
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)   ::   sudy_u, svdx_v  ! transports (u*dy, v*dx) for T-centric mesh (m^2/s)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)   ::   sudy_v, svdx_u  ! transports (u*dy, v*dx) for F-centric mesh (m^2/s)

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION ice_alloc()
      !!-----------------------------------------------------------------
      !!               *** Routine ice_alloc ***
      !!-----------------------------------------------------------------
      INTEGER :: ice_alloc
      !
      INTEGER :: ierr(34), ii
      !!-----------------------------------------------------------------
      ierr(:) = 0

      ii = 1
      ALLOCATE( ht_i_new (jpi,jpj) , fraz_frac (jpi,jpj) ,  &
         &      strength (jpi,jpj) , stress1_i(jpi,jpj) , stress2_i(jpi,jpj) , stress12_i(jpi,jpj), &
         &      delta_i  (jpi,jpj) , divu_i   (jpi,jpj) , shear_i  (jpi,jpj) , STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 2D arrays to memory'
      PRINT *, '            => ht_i_new, fraz_frac, strength, stress1_i, stress2_i, stress12_i, delta_i, divu_i, shear_i'
      !$acc enter data copyin( ht_i_new, fraz_frac, strength, stress1_i, stress2_i, stress12_i, delta_i, divu_i, shear_i  )
# endif

      !LOLOfixme: make it inside a `ln_icethd` flag here?
      ii = ii + 1
      ALLOCATE( t_bo       (jpi,jpj) , wfx_snw_sni(jpi,jpj) ,                                                &
         &      wfx_snw    (jpi,jpj) , wfx_snw_dyn(jpi,jpj) , wfx_snw_sum(jpi,jpj) , wfx_snw_sub(jpi,jpj) ,  &
         &      wfx_ice    (jpi,jpj) , wfx_sub    (jpi,jpj) , wfx_ice_sub(jpi,jpj) , wfx_lam    (jpi,jpj) ,  &
         &      wfx_pnd    (jpi,jpj) ,                                                                       &
         &      wfx_bog    (jpi,jpj) , wfx_dyn   (jpi,jpj) , wfx_bom(jpi,jpj) , wfx_sum(jpi,jpj) ,           &
         &      wfx_res    (jpi,jpj) , wfx_sni   (jpi,jpj) , wfx_opw(jpi,jpj) , wfx_spr(jpi,jpj) ,           &
         &      qsb_ice_bot(jpi,jpj) , qlead     (jpi,jpj) ,                                                 &
         &      sfx_res    (jpi,jpj) , sfx_bri   (jpi,jpj) , sfx_dyn(jpi,jpj) , sfx_sub(jpi,jpj) , sfx_lam(jpi,jpj) ,  &
         &      sfx_bog    (jpi,jpj) , sfx_bom   (jpi,jpj) , sfx_sum(jpi,jpj) , sfx_sni(jpi,jpj) , sfx_opw(jpi,jpj) ,  &
         &      hfx_res    (jpi,jpj) , hfx_snw   (jpi,jpj) , hfx_sub(jpi,jpj) ,                        &
         &      qt_atm_oi  (jpi,jpj) , qt_oce_ai (jpi,jpj) , fhld   (jpi,jpj) ,                        &
         &      hfx_sum    (jpi,jpj) , hfx_bom   (jpi,jpj) , hfx_bog(jpi,jpj) , hfx_dif(jpi,jpj) ,     &
         &      hfx_opw    (jpi,jpj) , hfx_thd   (jpi,jpj) , hfx_dyn(jpi,jpj) , hfx_spr(jpi,jpj) ,     &
         &      hfx_err_dif(jpi,jpj) , wfx_err_sub(jpi,jpj)                   , STAT=ierr(ii) )

      t_bo(:,:)=0._wp ; wfx_snw_sni(:,:)=0._wp
      wfx_snw(:,:)=0._wp ; wfx_snw_dyn(:,:)=0._wp ; wfx_snw_sum(:,:)=0._wp ; wfx_snw_sub(:,:)=0._wp
      wfx_ice(:,:)=0._wp ; wfx_sub(:,:)=0._wp ; wfx_ice_sub(:,:)=0._wp ; wfx_lam(:,:)=0._wp
      wfx_pnd(:,:)=0._wp
      wfx_bog(:,:)=0._wp ; wfx_dyn(:,:)=0._wp ; wfx_bom(:,:)=0._wp ; wfx_sum(:,:)=0._wp
      wfx_res(:,:)=0._wp ; wfx_sni(:,:)=0._wp ; wfx_opw(:,:)=0._wp ; wfx_spr(:,:)=0._wp
      qsb_ice_bot(:,:)=0._wp ; qlead(:,:)=0._wp
      sfx_res(:,:)=0._wp ; sfx_bri(:,:)=0._wp ; sfx_dyn(:,:)=0._wp ; sfx_sub(:,:)=0._wp ; sfx_lam(:,:)=0._wp
      sfx_bog(:,:)=0._wp ; sfx_bom(:,:)=0._wp ; sfx_sum(:,:)=0._wp ; sfx_sni(:,:)=0._wp ; sfx_opw(:,:)=0._wp
      hfx_res(:,:)=0._wp ; hfx_snw(:,:)=0._wp ; hfx_sub(:,:)=0._wp
      qt_atm_oi(:,:)=0._wp ; qt_oce_ai(:,:)=0._wp ; fhld(:,:)=0._wp
      hfx_sum(:,:)=0._wp ; hfx_bom(:,:)=0._wp ; hfx_bog(:,:)=0._wp ; hfx_dif(:,:)=0._wp
      hfx_opw(:,:)=0._wp ; hfx_thd(:,:)=0._wp ; hfx_dyn(:,:)=0._wp ; hfx_spr(:,:)=0._wp
      hfx_err_dif(:,:)=0._wp ; wfx_err_sub(:,:)=0._wp


# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 2D thermo-only arrays to memory'
      PRINT *, '            => t_bo, wfx_snw_sni, wfx_snw, wfx_snw_dyn, wfx_snw_sum, wfx_snw_sub, wfx_ice, wfx_sub, wfx_ice_sub'
      !$acc enter data copyin( t_bo, wfx_snw_sni, wfx_snw, wfx_snw_dyn, wfx_snw_sum, wfx_snw_sub, wfx_ice, wfx_sub, wfx_ice_sub )
      PRINT *, '            => wfx_lam, wfx_pnd, wfx_bog, wfx_dyn, wfx_bom, wfx_sum, wfx_res, wfx_sni, wfx_opw, wfx_spr'
      !$acc enter data copyin( wfx_lam, wfx_pnd, wfx_bog, wfx_dyn, wfx_bom, wfx_sum, wfx_res, wfx_sni, wfx_opw, wfx_spr )
      PRINT *, '            => qsb_ice_bot, qlead, sfx_res, sfx_bri, sfx_dyn, sfx_sub, sfx_lam, sfx_bog, sfx_bom, sfx_sum'
      !$acc enter data copyin( qsb_ice_bot, qlead, sfx_res, sfx_bri, sfx_dyn, sfx_sub, sfx_lam, sfx_bog, sfx_bom, sfx_sum )
      PRINT *, '            => sfx_sni, sfx_opw, hfx_res, hfx_snw, hfx_sub, qt_atm_oi, qt_oce_ai, fhld, hfx_sum, hfx_bom'
      !$acc enter data copyin( sfx_sni, sfx_opw, hfx_res, hfx_snw, hfx_sub, qt_atm_oi, qt_oce_ai, fhld, hfx_sum, hfx_bom )
      PRINT *, '            => hfx_bog, hfx_dif, hfx_opw, hfx_thd, hfx_dyn, hfx_spr, hfx_err_dif, wfx_err_sub'
      !$acc enter data copyin( hfx_bog, hfx_dif, hfx_opw, hfx_thd, hfx_dyn, hfx_spr, hfx_err_dif, wfx_err_sub )
# endif
      !LOLOfixme.

      ! * Ice global state variables
      ii = ii + 1
      ALLOCATE( qtr_ice_bot(jpi,jpj,jpl) , cnd_ice(jpi,jpj,jpl) , t1_ice(jpi,jpj,jpl) ,  &
         &      h_i        (jpi,jpj,jpl) , a_i    (jpi,jpj,jpl) , v_i   (jpi,jpj,jpl) ,  &
         &      v_s        (jpi,jpj,jpl) , h_s    (jpi,jpj,jpl) , t_su  (jpi,jpj,jpl) ,  &
         &      s_i        (jpi,jpj,jpl) , sv_i   (jpi,jpj,jpl) , o_i   (jpi,jpj,jpl) ,  &
         &      oa_i       (jpi,jpj,jpl) , v_ibr   (jpi,jpj,jpl) , STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 3D arrays to memory'
      PRINT *, '            => qtr_ice_bot, cnd_ice, t1_ice, h_i, a_i, v_i, v_s, h_s, t_su, s_i, sv_i, o_i, oa_i'
      !$acc enter data copyin( qtr_ice_bot, cnd_ice, t1_ice, h_i, a_i, v_i, v_s, h_s, t_su, s_i, sv_i, o_i, oa_i )
      !!   --> ignored for now: v_ibr
# endif

      ii = ii + 1
      ALLOCATE( u_ice(jpi,jpj) , v_ice(jpi,jpj) , uVice(jpi,jpj), vUice(jpi,jpj), SIGMAt(jpi,jpj,3), &
         &      vt_i (jpi,jpj) , vt_s (jpi,jpj) , st_i(jpi,jpj) , at_i(jpi,jpj) , ato_i(jpi,jpj) ,   &
         &      et_i (jpi,jpj) , et_s (jpi,jpj) , tm_i(jpi,jpj) , tm_s(jpi,jpj) ,  &
         &      sm_i (jpi,jpj) , tm_su(jpi,jpj) , hm_i(jpi,jpj) , hm_s(jpi,jpj) ,  &
         &      om_i (jpi,jpj) , vm_ibr(jpi,jpj) , tau_icebfr(jpi,jpj), icb_mask(jpi,jpj), &
         &      af_i(jpi,jpj),  au_i(jpi,jpj),  av_i(jpi,jpj),   STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 2D arrays to memory'
      PRINT *, '            => u_ice, v_ice, uVice, vUice, SIGMAt'
      !$acc enter data copyin( u_ice, v_ice, uVice, vUice, SIGMAt )
      PRINT *, '            => vt_i, vt_s, st_i, at_i, ato_i, et_i, et_s, tm_i, tm_s, sm_i, tm_su, hm_i, hm_s, om_i, af_i, au_i, av_i'
      !$acc enter data copyin( vt_i, vt_s, st_i, at_i, ato_i, et_i, et_s, tm_i, tm_s, sm_i, tm_su, hm_i, hm_s, om_i, af_i, au_i, av_i )
      !! --> ignored for now: vm_ibr, tau_icebfr, icb_mask
# endif

      ii = ii + 1
      ALLOCATE( t_s(jpi,jpj,nlay_s,jpl) , e_s(jpi,jpj,nlay_s,jpl) , STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 4D snow arrays to memory'
      PRINT *, '            => t_s, e_s'
      !$acc enter data copyin( t_s, e_s )
# endif

      ii = ii + 1
      ALLOCATE( t_i(jpi,jpj,nlay_i,jpl) , e_i(jpi,jpj,nlay_i,jpl) , szv_i(jpi,jpj,nlay_i,jpl) , sz_i(jpi,jpj,nlay_i,jpl) , STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 4D ice arrays to memory'
      PRINT *, '            => t_i, e_i, szv_i, sz_i'
      !$acc enter data copyin( t_i, e_i, szv_i, sz_i )
# endif

      ii = ii + 1
      ALLOCATE( a_ip(jpi,jpj,jpl) , v_ip(jpi,jpj,jpl) , a_ip_frac(jpi,jpj,jpl) , h_ip(jpi,jpj,jpl),  &
         &      v_il(jpi,jpj,jpl) , h_il(jpi,jpj,jpl) , a_ip_eff (jpi,jpj,jpl) , STAT = ierr(ii) )

      ii = ii + 1
      ALLOCATE( at_ip(jpi,jpj) , hm_ip(jpi,jpj) , vt_ip(jpi,jpj) , hm_il(jpi,jpj) , vt_il(jpi,jpj) , STAT = ierr(ii) )

      ! * Old values of global variables
      ii = ii + 1
      ALLOCATE( v_s_b (jpi,jpj,jpl) , v_i_b (jpi,jpj,jpl) , h_s_b(jpi,jpj,jpl)        , h_i_b(jpi,jpj,jpl),         &
         &      v_ip_b(jpi,jpj,jpl) , v_il_b(jpi,jpj,jpl) , szv_i_b (jpi,jpj,nlay_i,jpl),                           &
         &      a_i_b (jpi,jpj,jpl) , sv_i_b(jpi,jpj,jpl) , e_i_b(jpi,jpj,nlay_i,jpl) , e_s_b(jpi,jpj,nlay_s,jpl) , &
         &      STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 3D/4D arrays to memory'
      PRINT *, '            => v_s_b, v_i_b, h_s_b, h_i_b, szv_i_b, a_i_b, sv_i_b, e_i_b, e_s_b'
      !$acc enter data copyin( v_s_b, v_i_b, h_s_b, h_i_b, szv_i_b, a_i_b, sv_i_b, e_i_b, e_s_b )
      !!   --> ignored for now: v_ip_b, v_il_b
# endif

      ii = ii + 1
      ALLOCATE( u_ice_b(jpi,jpj) , v_ice_b(jpi,jpj) , at_i_b(jpi,jpj) , STAT=ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding 2D arrays to memory'
      PRINT *, '            => u_ice_b, v_ice_b, at_i_b'
      !$acc enter data copyin( u_ice_b, v_ice_b, at_i_b )
# endif


      ! * Ice thickness distribution variables
      ii = ii + 1
      ALLOCATE( hi_max(0:jpl), hi_mean(jpl),  STAT=ierr(ii) )

      ! * Ice diagnostics
      ii = ii + 1
      ALLOCATE( diag_trp_vi(jpi,jpj) , diag_trp_vs (jpi,jpj) , diag_trp_ei(jpi,jpj),                      &
         &      diag_trp_es(jpi,jpj) , diag_trp_sv (jpi,jpj) , diag_heat  (jpi,jpj),                      &
         &      diag_sice  (jpi,jpj) , diag_vice   (jpi,jpj) , diag_vsnw  (jpi,jpj), diag_aice(jpi,jpj), diag_vpnd(jpi,jpj),  &
         &      diag_adv_mass(jpi,jpj), diag_adv_salt(jpi,jpj), diag_adv_heat(jpi,jpj), STAT=ierr(ii) )

      ! * Ice conservation
      ii = ii + 1
      ALLOCATE( diag_v (jpi,jpj) , diag_s (jpi,jpj) , diag_t (jpi,jpj),   &
         &      diag_fv(jpi,jpj) , diag_fs(jpi,jpj) , diag_ft(jpi,jpj), STAT=ierr(ii) )

      ! * SIMIP diagnostics
      ii = ii + 1
      ALLOCATE( t_si(jpi,jpj,jpl) , tm_si(jpi,jpj) , qcn_ice_bot(jpi,jpj,jpl) , qcn_ice_top(jpi,jpj,jpl) , STAT = ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding arrays to memory'
      PRINT *, '            => t_si, tm_si, qcn_ice_bot, qcn_ice_top'
      !$acc enter data copyin( t_si, tm_si, qcn_ice_bot, qcn_ice_top )
# endif

      ii = ii + 1
      ALLOCATE( kmsk_ice_t(jpi,jpj), kmsk_ice_f(jpi,jpj), kmsk_ice_u(jpi,jpj), kmsk_ice_v(jpi,jpj), &
         &           V_oce(jpi,jpj,4), hm_i_f(jpi,jpj),  STAT = ierr(ii) ) ! Uu_oce, Vv_oce, Uv_oce & Vu_oce
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding arrays to memory'
      PRINT *, '            => kmsk_ice_t, kmsk_ice_f, kmsk_ice_u, kmsk_ice_v, V_oce, hm_i_f'
      !$acc enter data copyin( kmsk_ice_t, kmsk_ice_f, kmsk_ice_u, kmsk_ice_v, V_oce, hm_i_f )
# endif

      ! * damage tracers, and components of internal stress tensor at both T- and F-points:
      IF( ln_damage ) THEN
         ii = ii + 1
         ALLOCATE( V_ts(jpi,jpj,4), dmdt(jpi,jpj), dmdf(jpi,jpj), SIGMAf(jpi,jpj,3),    STAT = ierr(ii) )
# if defined _OPENACC
         PRINT *, ' * info GPU: ice_alloc() => adding brittle-reology-specific 2D arrays to memory'
         PRINT *, '            => dmdt, dmdf, SIGMAf, V_ts'
         !$acc enter data copyin( dmdt, dmdf, SIGMAf, V_ts )
# endif
      END IF

      !LOLOfixme: make it inside a `ln_icethd` flag here?
      ii = ii + 1
      ALLOCATE( s_i_new(jpi,jpj),  dh_i_itm(jpi,jpj), dh_s_tot(jpi,jpj),    &
         &      dh_i_bom(jpi,jpj), dh_s_itm(jpi,jpj), &
         &      dh_i_sub(jpi,jpj), dh_i_bog(jpi,jpj), dh_snowice(jpi,jpj), &
         &      eh_i_old(jpi,jpj,0:nlay_i+1), h_i_old(jpi,jpj,0:nlay_i+1), &
         &      dh_i_sum_2d(jpi,jpj,jpl), dh_s_sum_2d(jpi,jpj,jpl),  STAT = ierr(ii) )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding thermo-only arrays to memory'
      PRINT *, '            => s_i_new, dh_i_itm, dh_s_tot, dh_i_bom, dh_s_itm, dh_i_sub, dh_i_bog, dh_snowice, eh_i_old, h_i_old, dh_i_sum_2d, dh_s_sum_2d'
      !$acc enter data copyin( s_i_new, dh_i_itm, dh_s_tot, dh_i_bom, dh_s_itm, dh_i_sub, dh_i_bog, dh_snowice, eh_i_old, h_i_old, dh_i_sum_2d, dh_s_sum_2d )
# endif
      !LOLO.

      ii = ii + 1
      ALLOCATE( sudy_u(jpi,jpj), svdx_v(jpi,jpj),   STAT = ierr(2) )
      sudy_u(:,:) = 0._wp ;  svdx_v(:,:) = 0._wp
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_alloc() => adding advection transport arrays'
      PRINT *, '            => sudy_u, svdx_v'
      !$acc enter data copyin( sudy_u, svdx_v )
# endif

      IF( ln_damage ) THEN
         ii = ii + 1
         ALLOCATE( sudy_v(jpi,jpj), svdx_u(jpi,jpj),   STAT = ierr(2) )
         sudy_v(:,:) = 0._wp ;  svdx_u(:,:) = 0._wp
# if defined _OPENACC
         PRINT *, ' * info GPU: ice_alloc() => adding advection transport arrays'
         PRINT *, '            => sudy_v, svdx_u'
         !$acc enter data copyin( sudy_v, svdx_u )
# endif
      ENDIF


      ice_alloc = MAXVAL( ierr(:) )
      IF( ice_alloc /= 0 )   CALL ctl_stop( 'STOP', 'ice_alloc: failed to allocate arrays.' )
      !

   END FUNCTION ice_alloc

   !!======================================================================
END MODULE ice
