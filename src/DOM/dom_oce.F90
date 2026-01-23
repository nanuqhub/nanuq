MODULE dom_oce
   !!======================================================================
   !!                       ***  MODULE dom_oce  ***
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!======================================================================
   !! History :  1.0  ! 2005-10  (A. Beckmann, G. Madec)  reactivate s-coordinate
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012     (S. Mocavero, I. Epicoco) Add arrays associated
   !!                             to the optimization of BDY communications
   !!            3.7  ! 2015-11  (G. Madec) introduce surface and scale factor ratio
   !!             -   ! 2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!            4.1  ! 2019-08  (A. Coward, D. Storkey) rename prognostic variables in preparation for new time scheme.
   !!            4.x  ! 2020-02  (G. Madec, S. Techene) introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Root    : dummy function used when lk_agrif=F
   !!   Agrif_Fixed   : dummy function used when lk_agrif=F
   !!   Agrif_CFixed  : dummy function used when lk_agrif=F
   !!   dom_oce_alloc : dynamical allocation of dom_oce arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PUBLIC             ! allows the acces to par_oce when dom_oce is used (exception to coding rules)

   PUBLIC dom_oce_alloc  ! Called from nanuqgcm.F90

   INTEGER, PUBLIC, PARAMETER :: iverbose = 2

   INTEGER, PUBLIC, SAVE :: jperio    ! => as read into `domain_cfg.nc` file...
#if defined _OPENACC
   LOGICAL, PUBLIC, SAVE :: l_1c1g    ! => we use 1 CPU core together with 1 GPU
   LOGICAL, PUBLIC, SAVE :: l_need4ll ! => the ocean domain is not entirely landlocked by solid boundaries
   !                                  !    ==> `lbc_linking` must be used, even when no MPP !
#endif
   
   !!----------------------------------------------------------------------
   !! time & space domain namelist
   !! ----------------------------
   !                                   !!* Namelist namdom : time & space domain *
   LOGICAL , PUBLIC ::   ln_meshmask    !: =T  create a mesh-mask file (mesh_mask.nc)
   REAL(wp), PUBLIC ::   rn_Dt          !: time step for the dynamics and tracer

   !! Free surface parameters
   !! =======================
   LOGICAL , PUBLIC :: ln_dynspg_exp    !: Explicit free surface flag
   LOGICAL , PUBLIC :: ln_dynspg_ts     !: Split-Explicit free surface flag

   !! Time splitting parameters
   !! =========================
   LOGICAL,  PUBLIC :: ln_bt_fw         !: Forward integration of barotropic sub-stepping
   LOGICAL,  PUBLIC :: ln_bt_av         !: Time averaging of barotropic variables
   LOGICAL,  PUBLIC :: ln_bt_auto       !: Set number of barotropic iterations automatically
   INTEGER,  PUBLIC :: nn_bt_flt        !: Filter choice
   INTEGER,  PUBLIC :: nn_e          !: Number of barotropic iterations during one baroclinic step (rn_Dt)
   REAL(wp), PUBLIC :: rn_bt_cmax       !: Maximum allowed courant number (used if ln_bt_auto=T)
   REAL(wp), PUBLIC :: rn_bt_alpha      !: Time stepping diffusion parameter


   !                                   !!! associated variables
   REAL(wp), PUBLIC ::   rDt, r1_Dt     !: Current model timestep and reciprocal

   !!----------------------------------------------------------------------
   !! space domain parameters
   !!----------------------------------------------------------------------
   LOGICAL         , PUBLIC ::   l_Iperio, l_Jperio   ! i- j-periodicity
   LOGICAL         , PUBLIC ::   l_NFold              ! North Pole folding
   CHARACTER(len=1), PUBLIC ::   c_NFtype             ! type of North pole Folding: T or F point

   ! Tiling namelist
   LOGICAL, PUBLIC ::   ln_tile
   INTEGER         ::   nn_ltile_i, nn_ltile_j

   ! Domain tiling
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ntsi_a       !: start of internal part of tile domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ntsj_a       !
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ntei_a       !: end of internal part of tile domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ntej_a       !
   LOGICAL, PUBLIC                                  ::   l_istiled    ! whether tiling is currently active or not

   !                             !: domain MPP decomposition parameters
   INTEGER              , PUBLIC ::   nimpp, njmpp     !: i- & j-indexes for mpp-subdomain left bottom
   INTEGER              , PUBLIC ::   narea            !: number for local area (starting at 1) = MPI rank + 1
   INTEGER,               PUBLIC ::   nidom      !: IOIPSL things...

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig        !: local ==> global domain, including halos (jpiglo), i-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg        !: local ==> global domain, including halos (jpjglo), j-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig0       !: local ==> global domain, excluding halos (Ni0glo), i-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg0       !: local ==> global domain, excluding halos (Nj0glo), j-index
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mi0, mi1   !: global, including halos (jpiglo) ==> local domain i-index
   !                                                                !:    (mi0=1 and mi1=0 if global index not in local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   mj0, mj1   !: global, including halos (jpjglo) ==> local domain j-index
   !                                                                !:    (mj0=1 and mj1=0 if global index not in local domain)
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nfimpp, nfproc, nfjpi

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   glamt , glamu, glamv , glamf    !: longitude at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   gphit , gphiu, gphiv , gphif    !: latitude  at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: r1_e1t, r1_e2t!: t-point horizontal scale factors    [m]
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e1t, e2t!: t-point horizontal scale factors    [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e2u, r1_e1u, r1_e2u!: horizontal scale factors at u-point [m]
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e1u!: horizontal scale factors at u-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e1v, r1_e1v, r1_e2v!: horizontal scale factors at v-point [m]
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e2v!: horizontal scale factors at v-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: r1_e1f, r1_e2f!: horizontal scale factors at f-point [m]
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)  :: e1f, e2f!: horizontal scale factors at f-point [m]
   !
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2t , r1_e1e2t                !: associated metrics at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2u , r1_e1e2u , e2_e1u       !: associated metrics at u-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2v , r1_e1e2v , e1_e2v       !: associated metrics at v-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1e2f , r1_e1e2f                !: associated metrics at f-point
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1t2, e2t2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e1f2, e2f2
   !
   REAL(dp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   res_grd_loc_t, res_grd_loc_f    !: local mean size of the mesh
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ff_u, ff_v                !: Coriolis factor at U- & V-points  [1/s]

   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! ---------------------------------------------------------------------
   LOGICAL, PUBLIC ::   ln_zco       !: z-coordinate - full step
   LOGICAL, PUBLIC ::   ln_zps       !: z-coordinate - partial step
   LOGICAL, PUBLIC ::   ln_sco       !: s-coordinate or hybrid z-s coordinate
   LOGICAL, PUBLIC ::   ln_isfcav    !: presence of ISF
   !                                                        ! time-dependent heights of ocean water column   (domvvl)
   INTEGER, PUBLIC ::   nla10              !: deepest    W level Above  ~10m (nlb10 - 1)
   INTEGER, PUBLIC ::   nlb10              !: shallowest W level Bellow ~10m (nla10 + 1)

   !!----------------------------------------------------------------------
   !! masks, top and bottom ocean point position
   !! ---------------------------------------------------------------------
   REAL(wp),   PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tmask_i                  !: interior (excluding halos+duplicated points) domain T-point mask
   REAL(wp),   PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:), TARGET :: tmask, umask, vmask, fmask  !: land/ocean mask at T-, U-, V- and F-pts
   REAL(wp),   PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: xmskt, xmskf, xmsku     ! for ice dynamics

   INTEGER(1), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: klbct, klbcf, klbcu ! masks for weno5 solid latera BCs, last axis: E,N,W,S
   !INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: kmEt, kmWt
   !INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: kmEf, kmWf
   !INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: kmNt, kmSt
   !INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: kmNf, kmSf


   
   !!----------------------------------------------------------------------
   !! calendar variables
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC ::   nyear         !: current year
   INTEGER , PUBLIC ::   nmonth        !: current month
   INTEGER , PUBLIC ::   nday          !: current day of the month
   INTEGER , PUBLIC ::   nhour         !: current hour
   INTEGER , PUBLIC ::   nminute       !: current minute
   INTEGER , PUBLIC ::   ndastp        !: time step date in yyyymmdd format
   INTEGER , PUBLIC ::   nday_year     !: current day counted from jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_year     !: seconds between 00h jan 1st of the current  year and half of the current time step
   INTEGER , PUBLIC ::   nsec_month    !: seconds between 00h 1st day of the current month and half of the current time step
   INTEGER , PUBLIC ::   nsec_monday   !: seconds between 00h         of the last Monday   and half of the current time step
   INTEGER , PUBLIC ::   nsec_day      !: seconds between 00h         of the current   day and half of the current time step
   REAL(dp), PUBLIC ::   fjulday       !: current julian day
   REAL(dp), PUBLIC ::   fjulstartyear !: first day of the current year in julian days
   REAL(wp), PUBLIC ::   adatrj        !: number of elapsed days since the begining of the whole simulation
   !                                   !: (cumulative duration of previous runs that may have used different time-step size)
   INTEGER , PUBLIC, DIMENSION(  0: 2) ::   nyear_len     !: length in days of the previous/current/next year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_len    !: length in days of the months of the current year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_beg    !: second since Jan 1st 0h of the current year and the half of the months
   INTEGER , PUBLIC                  ::   nsec1jan000     !: second since Jan 1st 0h of nit000 year and Jan 1st 0h the current year
   INTEGER , PUBLIC                  ::   nsec000_1jan000   !: second since Jan 1st 0h of nit000 year and nit000
   INTEGER , PUBLIC                  ::   nsecend_1jan000   !: second since Jan 1st 0h of nit000 year and nitend

   !!----------------------------------------------------------------------
   !! variable defined here to avoid circular dependencies...
   !! ---------------------------------------------------------------------
   INTEGER, PUBLIC ::   nbasin         ! number of basin to be considered in diaprt (glo, atl, pac, ind, ipc)

   !!----------------------------------------------------------------------
   !! agrif domain
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .FALSE.   !: agrif flag

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: dom_oce.F90 15556 2021-11-29 15:23:06Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   LOGICAL FUNCTION Agrif_Root()
      Agrif_Root = .TRUE.
   END FUNCTION Agrif_Root

   INTEGER FUNCTION Agrif_Fixed()
      Agrif_Fixed = 0
   END FUNCTION Agrif_Fixed

   CHARACTER(len=3) FUNCTION Agrif_CFixed()
      Agrif_CFixed = '0'
   END FUNCTION Agrif_CFixed

   INTEGER FUNCTION dom_oce_alloc()
      !!----------------------------------------------------------------------
      INTEGER                ::   ii
      INTEGER, DIMENSION(30) :: ierr
      !!----------------------------------------------------------------------
      ii = 0   ;   ierr(:) = 0
      !
      ii = ii+1
      ALLOCATE( glamt(jpi,jpj) ,    glamu(jpi,jpj) ,  glamv(jpi,jpj) ,  glamf(jpi,jpj) ,     &
         &      gphit(jpi,jpj) ,    gphiu(jpi,jpj) ,  gphiv(jpi,jpj) ,  gphif(jpi,jpj) ,     &
         &       e1t (jpi,jpj) ,     e2t (jpi,jpj) , r1_e1t(jpi,jpj) , r1_e2t(jpi,jpj) ,     &
         &       e1u (jpi,jpj) ,     e2u (jpi,jpj) , r1_e1u(jpi,jpj) , r1_e2u(jpi,jpj) ,     &
         &       e1v (jpi,jpj) ,     e2v (jpi,jpj) , r1_e1v(jpi,jpj) , r1_e2v(jpi,jpj) ,     &
         &       e1f (jpi,jpj) ,     e2f (jpi,jpj) , r1_e1f(jpi,jpj) , r1_e2f(jpi,jpj) ,     &
         &      e1e2t(jpi,jpj) , r1_e1e2t(jpi,jpj)                                     ,     &
         &      e1e2u(jpi,jpj) , r1_e1e2u(jpi,jpj) , e2_e1u(jpi,jpj)                   ,     &
         &      e1e2v(jpi,jpj) , r1_e1e2v(jpi,jpj) , e1_e2v(jpi,jpj)                   ,     &
         &      e1e2f(jpi,jpj) , r1_e1e2f(jpi,jpj) , ff_u(jpi,jpj),    ff_v(jpi,jpj)   ,     &
         &      res_grd_loc_t(jpi,jpj) ,  res_grd_loc_f(jpi,jpj),    STAT=ierr(ii) )
      !
      ii = ii+1
      ALLOCATE( e1t2(jpi,jpj), e2t2(jpi,jpj), xmskt(jpi,jpj), &
         &      e1f2(jpi,jpj), e2f2(jpi,jpj), xmskf(jpi,jpj), &
         &      xmsku(jpi,jpj),   STAT=ierr(ii) )
      !
      ii = ii+1
      ALLOCATE(    klbct(jpi,jpj,nn_hls,4), klbcf(jpi,jpj,nn_hls,4), klbcu(jpi,jpj,nn_hls,4),  STAT = ierr(ii) )
      !
      !
      PRINT *, 'LOLO `dom_oce_alloc@dom_oce.F90` allocating mask arrays with jpk =', jpk, ', proc #', narea
      !
      ii = ii+1
      ALLOCATE( tmask(jpi,jpj,jpk), umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), fmask(jpi,jpj,jpk), tmask_i(jpi,jpj), STAT=ierr(ii) )
      !
      dom_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION dom_oce_alloc
   
   !!======================================================================
END MODULE dom_oce
