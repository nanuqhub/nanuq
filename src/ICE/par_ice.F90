MODULE par_ice
   !!======================================================================
   !!                        ***  par_ice  ***
   !! Ice :   set the ice parameters
   !!======================================================================
   !! History :  4.x  !  2023     (Rousset)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   LOGICAL,  PARAMETER, PUBLIC :: l_use_v_for_h = .FALSE.   ! `h` or `v` in rheology for `h` as in vertically integrated stresses ? #LB: we should not! But `h` is sloppy by essence in SI3...
   REAL(wp), PARAMETER, PUBLIC :: rconc_min  = 0.001_wp      ! ice concentration below which ice velocity becomes very small
   REAL(wp), PARAMETER, PUBLIC :: rmass_min  = 1._wp         ! ice mass (kg/m2)  below which ice velocity becomes very small
   REAL(wp), PARAMETER, PUBLIC :: rAmin_fld  = 0.01_wp       ! ice concentration below which we hide relevant fields in output netCDF files...
   REAL(wp), PARAMETER, PUBLIC :: rAmin_dmg  = 0.1_wp        ! ice concentration below which ice damage & stresses becomes a nonsense and are forced to 0
   REAL(wp), PARAMETER, PUBLIC :: rcloud_fra = 0.81_wp       ! prescribed cloud fraction
   !$acc declare create( rconc_min, rmass_min, rAmin_fld, rAmin_dmg, rcloud_fra )


   !!----------------------------------------------------------------------
   !!                   shared namelist parameters
   !!----------------------------------------------------------------------
   !                                     !!** ice-generic parameters namelist (nampar) **
   INTEGER           , PUBLIC ::   jpl              !: number of ice  categories
   INTEGER           , PUBLIC ::   nlay_i           !: number of ice  layers
   INTEGER           , PUBLIC ::   nlay_s           !: number of snow layers
   LOGICAL           , PUBLIC ::   ln_virtual_itd   !: virtual ITD mono-category parameterization (T) or not (F)
   LOGICAL           , PUBLIC ::   ln_icedyn        !: flag for ice dynamics (T) or not (F)
   LOGICAL           , PUBLIC ::   ln_icethd        !: flag for ice thermo   (T) or not (F)
   REAL(wp)          , PUBLIC ::   rn_amax          !: maximum ice concentration
   CHARACTER(len=256), PUBLIC ::   cn_icerst_in     !: suffix of ice restart name (input)
   CHARACTER(len=256), PUBLIC ::   cn_icerst_out    !: suffix of ice restart name (output)
   CHARACTER(len=256), PUBLIC ::   cn_icerst_indir  !: ice restart input directory
   CHARACTER(len=256), PUBLIC ::   cn_icerst_outdir !: ice restart output directory
   LOGICAL           , PUBLIC ::   ln_damage        !: existence of the "ice damage" tracer in the code? => (T) if `ln_rhg_BBM==.true.`
   !$acc declare create( jpl, nlay_i, nlay_s, rn_amax )

   !                                     !! ** namelist (namini) **
   LOGICAL, PUBLIC  ::   ln_iceini        !: Ice initialization or not
   INTEGER, PUBLIC  ::   nn_iceini_file   !: Ice initialization:
   !                                      !        0 = Initialise sea ice based on SSTs
   !                                      !        1 = Initialise sea ice from single category netcdf file
   !                                      !        2 = Initialise sea ice from multi category restart file

   !                                     !!** ice-itd namelist (namitd) **
   LOGICAL , PUBLIC                   ::   ln_cat_hfn   ! ice categories are defined by function like rn_himean**(-0.05)
   REAL(wp), PUBLIC                   ::   rn_himean    ! mean thickness of the domain
   LOGICAL , PUBLIC                   ::   ln_cat_usr   ! ice categories are defined by rn_catbnd
   REAL(wp), PUBLIC, DIMENSION(0:100) ::   rn_catbnd    ! ice categories bounds
   REAL(wp), PUBLIC                   ::   rn_himin         !: minimum ice thickness
   REAL(wp), PUBLIC                   ::   rn_himax     ! maximum ice thickness allowed
   !$acc declare create( ln_cat_hfn, rn_himean, ln_cat_usr, rn_catbnd, rn_himin, rn_himax )

   !                                     !!** ice-dynamics namelist (namdyn) **
   REAL(wp), PUBLIC ::   rn_ishlat        !: lateral boundary condition for sea-ice
   LOGICAL , PUBLIC ::   ln_landfast_L16  !: landfast ice parameterizationfrom lemieux2016
   REAL(wp), PUBLIC ::   rn_lf_depfra     !:    fraction of ocean depth that ice must reach to initiate landfast ice
   REAL(wp), PUBLIC ::   rn_lf_bfr        !:    maximum bottom stress per unit area of contact (lemieux2016) or per unit volume (home)
   REAL(wp), PUBLIC ::   rn_lf_relax      !:    relaxation time scale (s-1) to reach static friction
   REAL(wp), PUBLIC ::   rn_lf_tensile    !:    isotropic tensile strength
   !
   ! ** namelist (namdyn) **
   LOGICAL, PUBLIC, SAVE :: ln_dynALL     ! full ice dynamics                      (rheology + advection + ridging/rafting + correction)
   LOGICAL, PUBLIC, SAVE :: ln_dynRHGADV  ! no ridge/raft & no corrections         (rheology + advection)
   LOGICAL, PUBLIC, SAVE :: ln_dynADV1D   ! only advection in 1D w ice convergence (test case from Schar & Smolarkiewicz 1996)
   LOGICAL, PUBLIC, SAVE :: ln_dynADV2D   ! only advection in 2D w prescribed vel.
   !                                      !   => uses ocean velocities provided in netCDF SSX forcing as sea-ice velocities
   LOGICAL, PUBLIC, SAVE :: ln_pureADV2D  !   `T` => only pure advection, do not apply any king of corrections on advected fields
   !
   !           !!** ice-ridging/rafting namelist (namdyn_rdgrft) **
   REAL(wp), PUBLIC ::   rn_delta_ecc     !: eccentricity of the elliptical yield curve to use when computing `delta`
   !                                      !:  => overwritten with `rn_ecc` (&namdyn_rhg) when VP rheologies are used
   !                                      !:  `rn_delta_ecc = 1` implies that `delta` = `total deformation`
   LOGICAL , PUBLIC ::   ln_str_H79       !: ice strength parameterization: Hibler 79 (can be used in rheology)                       
   REAL(wp), PUBLIC ::   rn_crhg          !:    determines changes in ice strength                                                    
   REAL(wp), PUBLIC ::   rn_pstar         !:    determines ice strength, Hibler 79                                                    
   LOGICAL , PUBLIC ::   ln_str_R75       !: ice strength parameterization: Rothrock 75                                               
   REAL(wp), PUBLIC ::   rn_pe_rdg        !:    coef accounting for frictional dissipation                                            
   LOGICAL , PUBLIC ::   ln_str_CST       !: ice strength parameterization: Constant                                                  
   REAL(wp), PUBLIC ::   rn_str           !:    constant value of ice strength
   !$acc declare create( rn_delta_ecc, ln_str_H79, rn_crhg, rn_pstar, ln_str_R75, rn_pe_rdg, ln_str_CST, rn_str )
   !
   !                                     !!** ice-rheology namelist (namdyn_rhg) **
   LOGICAL,  PUBLIC ::   ln_idealized     !: if set to true: Coriolis and SSH tilt terms will be disregarded in the momentum eq.
   REAL(wp), PUBLIC ::   ridlzd
   ! -- evp
   LOGICAL , PUBLIC ::   ln_rhg_EVP       ! EVP rheology switch, used for rdgrft and rheology
   REAL(wp), PUBLIC ::   rn_creepl        !: creep limit (has to be low enough, circa 10-9 m/s, depending on rheology)
   REAL(wp), PUBLIC ::   rn_ecc           !: eccentricity of the elliptical yield curve
   INTEGER , PUBLIC ::   nn_nevp          !: number of iterations for subcycling
   REAL(wp), PUBLIC ::   rn_relast        !: ratio => telast/rDt_ice (1/3 or 1/9 depending on nb of subcycling nevp)
   INTEGER , PUBLIC ::   nn_rhg_chkcvg    !: check ice rheology convergence
   !$acc declare create( ln_rhg_EVP, rn_creepl, rn_ecc, nn_nevp, rn_relast )
   !
   ! -- bbm
   LOGICAL , PUBLIC ::   ln_rhg_BBM       ! BBM rheology #bbm   
   INTEGER , PUBLIC ::   nbbm, nflt       !: number of iterations for subcycling !#bbm
   REAL(wp), PUBLIC ::   rn_Nref          !: Maximum compressive stress at the reference scale [Pa] / neXtSIM => `compr_strength`
   REAL(wp), PUBLIC ::   rn_P0            !: Compression factor "P" at play in P_max, in Eq.8 of [Olason al.2022]
   REAL(wp), PUBLIC ::   rn_E0            !: Elasticity of undamaged ice [Pa]
   REAL(wp), PUBLIC ::   rn_eta0          !: Viscosity of Undamaged ice [Pa.s]
   REAL(wp), PUBLIC ::   rn_kth           !: healing constant [Eq.30 of Olason et al.,2022]
   REAL(wp), PUBLIC ::   rn_bbm_flt       !: filtering window for ice velocities 
   INTEGER , PUBLIC ::   nn_d_adv         !: advection of damage and stress tensor @T and @F
   LOGICAL,  PUBLIC ::   ln_adv_d_pra     !: use Prather advection scheme to advect damage (and stress components)
   LOGICAL,  PUBLIC ::   ln_adv_d_wnx     !: use  WENO advection scheme to advect damage (and stress components)
   LOGICAL,  PUBLIC ::   ln_x_MC_test     !: Perform only 1 Mohr-Coulomb test, at mid-point between T & F points (implicit smoothing)
   REAL(wp), PUBLIC ::   rn_crndg         !: cross-nudging coefficient ... #bbm
   !$acc declare create( ln_rhg_BBM, nbbm, nflt, rn_Nref, rn_P0, rn_E0, rn_eta0, rn_kth, rn_bbm_flt, nn_d_adv, rn_crndg )
   !
   REAL(wp), PUBLIC ::   rn_dmg_max, r_dmd_min !: max value possible for capping damage (~1-eps)
   REAL(wp), PUBLIC ::   rn_C0            !: compaction parameter "C" (Hibler's exponential)                                         =  -20.
   INTEGER,  PUBLIC ::   nn_alrlx         !: `alpha` inherent to viscosity, Eq.10 of [Olason al.2022], used by Dansereau             =   5.
   INTEGER,  PUBLIC ::   nn_btrlx         !: `beta`: sort of an `alpha` to go in the exp[] of `lambda`, used by Olason   Boutin      =   5.
   REAL(wp), PUBLIC ::   rn_c_ref         !: Cohesion value at the lab scale                                                         = 2.E6
   REAL(wp), PUBLIC ::   rn_l_ref         !: scaling paramater for cohesion, `l_ref` in [Eq.30 of Olason et al.,2022]
   !$acc declare create( rn_dmg_max, r_dmd_min, rn_C0, nn_alrlx, nn_btrlx, rn_c_ref, rn_l_ref )
   !
   !                                      !!** ice-advection namelist (namdyn_adv) **
   LOGICAL, PUBLIC :: ln_adv_Pra       !: Prather        advection scheme
   LOGICAL, PUBLIC :: ln_adv_UMx       !: Ultimate-Macho advection scheme
   LOGICAL, PUBLIC :: ln_adv_WNx       !: WENO advection scheme
   INTEGER, PUBLIC :: nn_UMx       ! order of the UMx advection scheme
   INTEGER, PUBLIC :: nn_WNx       ! order of the WENO advection scheme
   INTEGER, PUBLIC :: kp_weno      ! = `(nn_WNx+1)/2`
   CHARACTER(len=512), PUBLIC :: cn_weno_wght
   !$acc declare create( nn_UMx, nn_WNx, kp_weno, cn_weno_wght )


   !! Remapping usinng centered WENO:
   !                                      !!** remapping namelist (namrmp) **
   LOGICAL,            PUBLIC :: ln_use_weno_rmp
   CHARACTER(len=256), PUBLIC :: cn_rmp_weno_wght
   


   !                                     !!** ice-surface boundary conditions namelist (namsbc) **
   REAL(wp), PUBLIC ::   rn_snwblow       !: coef. for partitioning of snowfall between leads and sea ice
   ! -- icethd_zdf and icealb -- !
   INTEGER , PUBLIC ::   nn_snwfra        !: calculate the fraction of ice covered by snow
   !                                      !   = 1  fraction = 1-exp(-0.2*rhos*hsnw) [MetO formulation]
   !                                      !   = 2  fraction = hsnw / (hsnw+0.02)    [CICE formulation]
   !$acc declare create( nn_snwfra, rn_snwblow )

   ! -- icesbc -- !
   REAL(wp), PUBLIC ::   rn_Cd_io         !: drag coefficient for oceanic stress
   LOGICAL , PUBLIC ::   ln_drgice_imp    !: use implicit ice-ocean drag
   !$acc declare create( rn_Cd_io, ln_drgice_imp )

   ! -- icethd_zdf -- !
   LOGICAL , PUBLIC ::   ln_cndflx        !: use conduction flux as surface boundary condition (instead of qsr and qns)
   LOGICAL , PUBLIC ::   ln_cndemulate    !: emulate conduction flux (if not provided)
   INTEGER , PUBLIC ::   nn_qtrice        !: Solar flux transmitted thru the surface scattering layer
   !$acc declare create( ln_cndflx, ln_cndemulate, nn_qtrice )

   !                                      ! Conduction flux as surface forcing or not
   INTEGER, PUBLIC, PARAMETER ::   np_cnd_OFF = 0  !: no forcing from conduction flux (ice thermodynamics forced via qsr and qns)
   INTEGER, PUBLIC, PARAMETER ::   np_cnd_ON  = 1  !: forcing from conduction flux (SM0L) (compute qcn and qsr_tr via sbcblk.F90 or *ccpl.F90)
   INTEGER, PUBLIC, PARAMETER ::   np_cnd_EMU = 2  !: emulate conduction flux via icethd_zdf.F90 (BL99) (1st round compute qcn and qsr_tr, 2nd round use it)
   !                                      !   = 0  Grenfell and Maykut 1977 (depends on cloudiness and is 0 when there is snow)
   !                                      !   = 1  Lebrun 2019 (equals 0.3 anytime with different melting/dry snw conductivities)

   !                                     !!** namelist (namthd) **
   LOGICAL , PUBLIC ::   ln_icedH         ! activate ice thickness change from growing/melting (T) or not (F)
   LOGICAL , PUBLIC ::   ln_icedA         ! activate lateral melting param. (T) or not (F)
   LOGICAL , PUBLIC ::   ln_icedO         ! activate ice growth in open-water (T) or not (F)
   LOGICAL , PUBLIC ::   ln_leadhfx       ! heat in the leads is used to melt sea-ice before warming the ocean
   !
   !                                     !!** ice-vertical diffusion namelist (namthd_zdf) **
   LOGICAL , PUBLIC ::   ln_zdf_BL99      !: heat diffusion follows Bitz and Lipscomb (1999)
   LOGICAL , PUBLIC ::   ln_cndi_U64      !: thermal conductivity: Untersteiner (1964)
   LOGICAL , PUBLIC ::   ln_cndi_P07      !: thermal conductivity: Pringle et al (2007)
   REAL(wp), PUBLIC ::   rn_cnd_s         !: thermal conductivity of the snow [W/m/K]
   REAL(wp), PUBLIC ::   rn_kappa_i       !: coef. for the extinction of radiation in sea ice, Grenfell et al. (2006) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_s       !: coef. for the extinction of radiation in snw (nn_qtrice=0) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_smlt    !: coef. for the extinction of radiation in melt snw (nn_qtrice=1) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_sdry    !: coef. for the extinction of radiation in dry  snw (nn_qtrice=1) [1/m]
   LOGICAL , PUBLIC ::   ln_zdf_chkcvg    !: check convergence of heat diffusion scheme
   !$acc declare create( ln_zdf_BL99, ln_cndi_U64, ln_cndi_P07, rn_cnd_s, rn_kappa_i, rn_kappa_s, rn_kappa_smlt, rn_kappa_sdry, ln_zdf_chkcvg )

   !                                     !!** ice-salinity namelist (namthd_sal) **
   INTEGER , PUBLIC ::   nn_icesal        !: salinity configuration used in the model
   !                                      ! 1 - constant salinity in both space and time
   !                                      ! 2 - prognostic salinity (s(z,t))
   !                                      ! 3 - salinity profile, constant in time
   !                                      ! 4 - flushing and gravity drainage
   REAL(wp), PUBLIC ::   rn_icesal        !: bulk salinity (ppt) in case of constant salinity
   REAL(wp), PUBLIC ::   rn_sinew         !: fraction of sss that is kept in new ice
   REAL(wp), PUBLIC ::   rn_simin         !: minimum ice salinity [PSU]
   LOGICAL , PUBLIC ::   ln_sal_chk       !: sanity checks for salt drainage and flushing
   INTEGER , PUBLIC ::   nn_liquidus      !: formulation of liquidus
   !                                        1 = linear liquidus
   !                                        2 = Vancopenolle et al (2019) formulation
   !                                        3 = Weast formulation (used in RJW2014)
   REAL(wp), PUBLIC ::   rn_sal_gd        !: restoring salinity for gravity drainage [PSU]
   REAL(wp), PUBLIC ::   rn_time_gd       !: restoring time constant for gravity drainage (= 20 days) [s]
   REAL(wp), PUBLIC ::   rn_sal_fl        !: restoring salinity for flushing [PSU]
   REAL(wp), PUBLIC ::   rn_time_fl       !: restoring time constant for gravity drainage (= 10 days) [s]
   INTEGER , PUBLIC ::   nn_sal_scheme    !: convection scheme
   LOGICAL , PUBLIC ::   ln_flushing      !: activate flushing
   LOGICAL , PUBLIC ::   ln_drainage      !: activate gravity drainage
   INTEGER , PUBLIC ::   nn_drainage      !: number of subcycles for gravity drainage
   INTEGER , PUBLIC ::   nn_flushing      !: number of subcycles for flushing
   REAL(wp), PUBLIC ::   rn_flushrate     !: rate of flushing (fraction of melt water used for flushing)
   REAL(wp), PUBLIC ::   rn_alpha_CW      !: Brine flow for CW1988
   REAL(wp), PUBLIC ::   rn_alpha_RJW     !: Brine flow for RJW2014
   REAL(wp), PUBLIC ::   rn_alpha_GN      !: Brine flow for GN2013 (kg/m3/s)
   REAL(wp), PUBLIC ::   rn_Rc_RJW        !: critical Rayleigh number for RJW
   REAL(wp), PUBLIC ::   rn_Rc_GN         !:                          for GN
   REAL(wp), PUBLIC ::   rn_sal_himin     !: min ice thickness for gravity drainage and flushing calculation
   REAL(wp), PUBLIC ::   rn_vbrc          !: critical brines volume above which flushing can occur
   !$acc declare create( nn_icesal, rn_icesal, rn_sinew, rn_simin, ln_sal_chk, nn_liquidus, rn_sal_gd, rn_time_gd, rn_sal_fl, rn_time_fl, nn_sal_scheme, ln_flushing )
   !$acc declare create( ln_drainage, nn_drainage, nn_flushing, rn_flushrate, rn_alpha_CW, rn_alpha_RJW, rn_alpha_GN, rn_Rc_RJW, rn_Rc_GN, rn_sal_himin, rn_vbrc )

   !                                     !!** ice-ponds namelist (namthd_pnd)
   LOGICAL , PUBLIC ::   ln_pnd           !: Melt ponds (T) or not (F)
   LOGICAL , PUBLIC ::   ln_pnd_TOPO      !: Topographic Melt ponds scheme (Flocco et al 2007, 2010)
   LOGICAL , PUBLIC ::   ln_pnd_LEV       !: Simple melt pond scheme
   REAL(wp), PUBLIC ::   rn_apnd_min      !: Minimum fraction of melt water contributing to ponds
   REAL(wp), PUBLIC ::   rn_apnd_max      !: Maximum fraction of melt water contributing to ponds
   REAL(wp), PUBLIC ::   rn_pnd_flush     !: Pond flushing efficiency (tuning parameter)
   LOGICAL , PUBLIC ::   ln_pnd_CST       !: Melt ponds scheme with constant fraction and depth
   REAL(wp), PUBLIC ::   rn_apnd          !: prescribed pond fraction (0<rn_apnd<1)
   REAL(wp), PUBLIC ::   rn_hpnd          !: prescribed pond depth    (0<rn_hpnd<1)
   LOGICAL,  PUBLIC ::   ln_pnd_lids      !: Allow ponds to have frozen lids
   LOGICAL,  PUBLIC ::   ln_pnd_rain      !: rain added to melt ponds
   LOGICAL , PUBLIC ::   ln_pnd_alb       !: melt ponds affect albedo
   REAL(wp), PUBLIC ::   rn_pnd_hl_min    !: pond lid thickness below which full pond area used in albedo calculation
   REAL(wp), PUBLIC ::   rn_pnd_hl_max    !: pond lid thickness above which ponds disappear from albedo calculation
   !
   !                                     !!** ice-diagnostics namelist (namdia) **
   LOGICAL , PUBLIC ::   ln_icediachk     !: flag for ice diag (T) or not (F)
   REAL(wp), PUBLIC ::   rn_icechk_cel    !: rate of ice spuriously gained/lost (at any gridcell)
   REAL(wp), PUBLIC ::   rn_icechk_glo    !: rate of ice spuriously gained/lost (globally)
   LOGICAL , PUBLIC ::   ln_icediahsb     !: flag for ice diag (T) or not (F)
   LOGICAL , PUBLIC ::   ln_icectl        !: flag for sea-ice points output (T) or not (F)
   INTEGER , PUBLIC ::   iiceprt          !: debug i-point
   INTEGER , PUBLIC ::   jiceprt          !: debug j-point

   !!----------------------------------------------------------------------
   !!                   shared other parameters
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC ::   kt_ice           !: iteration number
   REAL(wp), PUBLIC ::   rDt_ice          !: ice time step
   REAL(wp), PUBLIC ::   r1_Dt_ice        !: = 1. / rDt_ice
   REAL(wp), PUBLIC ::   r1_nlay_i        !: 1 / nlay_i
   REAL(wp), PUBLIC ::   r1_nlay_s        !: 1 / nlay_s
   REAL(wp), PUBLIC ::   rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft   !: conservation diagnostics
   REAL(wp), PUBLIC, PARAMETER ::   epsi06 = 1.e-06_wp  !: small number
   REAL(wp), PUBLIC, PARAMETER ::   epsi10 = 1.e-10_wp  !: small number
   REAL(wp), PUBLIC, PARAMETER ::   epsi20 = 1.e-20_wp  !: small number
   !$acc declare create( rDt_ice, r1_Dt_ice, r1_nlay_i, r1_nlay_s, epsi06, epsi10, epsi20 )

   !!======================================================================
END MODULE par_ice
