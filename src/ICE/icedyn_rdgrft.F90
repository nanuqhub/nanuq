MODULE icedyn_rdgrft
   !!======================================================================
   !!                       ***  MODULE icedyn_rdgrft ***
   !!    sea-ice : Mechanical impact on ice thickness distribution
   !!======================================================================
   !! History :       !  2006-02  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rdgrft       : ridging/rafting of sea ice
   !!   ice_dyn_rdgrft_init  : initialization of ridging/rafting of sea ice
   !!   ice_strength         : ice strength calculation
   !!----------------------------------------------------------------------
   USE par_ice
   USE phycst         ! physical constants (ocean directory)
   USE oss_nnq , ONLY : sss_s, sst_s   ! surface boundary condition: ocean fields
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_roundoff
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rdgrft        ! called by icestp
   PUBLIC   ice_dyn_rdgrft_init   ! called by icedyn
   PUBLIC   ice_strength          ! called by icedyn_rhg_evp

   INTEGER ::              nice_str   ! choice of the type of strength
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_strh79     = 1   ! Hibler 79
   INTEGER, PARAMETER ::   np_strr75     = 2   ! Rothrock 75
   INTEGER, PARAMETER ::   np_strcst     = 3   ! Constant value

   ! Variables shared among ridging subroutines
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   closing_net     ! net rate at which area is removed    (1/s)
   !                                                               ! (ridging ice area - area of new ridges) / dt
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   opning          ! rate of opening due to divergence/shear
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   closing_gross   ! rate at which area removed, not counting area of new ridges
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   apartf          ! participation function; fraction of ridging/closing associated w/ category n
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrmin           ! minimum ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrmax           ! maximum ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrexp           ! e-folding ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hraft           ! thickness of rafted ice
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hi_hrdg         ! thickness of ridging ice / mean ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   aridge          ! participating ice ridging
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   araft           ! participating ice rafting
   !
   ! For ridging diagnostics
   LOGICAL                                         ::   ll_diag_rdg     ! activate ridging diagnostics or not
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airdg1          ! ridging ice area loss
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airft1          ! rafting ice area loss
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airdg2          ! new ridged ice area gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airft2          ! new rafted ice area gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   opning_2d       ! lead opening rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   dairdg1dt       ! ridging ice area loss rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   dairft1dt       ! rafting ice area loss rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   dairdg2dt       ! new ridged ice area gain rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   dairft2dt       ! new rafted ice area gain rate diagnostic
   !
   REAL(wp), PARAMETER ::   hrdg_hi_min = 1.1_wp    ! min ridge thickness multiplier: min(hrdg/hi)
   REAL(wp), PARAMETER ::   hi_hrft     = 0.5_wp    ! rafting multiplier: (hi/hraft)
   !$acc declare create( np_strh79, np_strr75, np_strcst, hrdg_hi_min, hi_hrft )
   !
   ! ** namelist (namdyn_rdgrft) **
   LOGICAL  ::   ln_str_smooth    ! ice strength spatial smoothing
   LOGICAL  ::   ln_distf_lin     ! redistribution of ridged ice: linear (Hibler 1980)
   LOGICAL  ::   ln_distf_exp     ! redistribution of ridged ice: exponential (Lipscomb et al 2017)
   REAL(wp) ::   rn_murdg         !    gives e-folding scale of ridged ice (m^.5)
   REAL(wp) ::   rn_csrdg         ! fraction of shearing energy contributing to ridging
   LOGICAL  ::   ln_partf_lin     ! participation function linear (Thorndike et al. (1975))
   REAL(wp) ::   rn_gstar         !    fractional area of young ice contributing to ridging
   LOGICAL  ::   ln_partf_exp     ! participation function exponential (Lipscomb et al. (2007))
   REAL(wp) ::   rn_astar         !    equivalent of G* for an exponential participation function
   LOGICAL  ::   ln_ridging       ! ridging of ice or not
   REAL(wp) ::   rn_hstar         !    thickness that determines the maximal thickness of ridged ice
   REAL(wp) ::   rn_porordg       !    initial porosity of ridges (0.3 regular value)
   REAL(wp) ::   rn_fsnwrdg       !    fractional snow loss to the ocean during ridging
   REAL(wp) ::   rn_fpndrdg       !    fractional pond loss to the ocean during ridging
   LOGICAL  ::   ln_rafting       ! rafting of ice or not
   REAL(wp) ::   rn_hraft         !    threshold thickness (m) for rafting / ridging
   REAL(wp) ::   rn_craft         !    coefficient for smoothness of the hyperbolic tangent in rafting
   REAL(wp) ::   rn_fsnwrft       !    fractional snow loss to the ocean during rafting
   REAL(wp) ::   rn_fpndrft       !    fractional pond loss to the ocean during rafting
   !
   !$acc declare create( nice_str, ln_str_smooth, ln_distf_lin, ln_distf_exp, rn_murdg, rn_csrdg, ln_partf_lin, rn_gstar, ln_partf_exp, rn_astar )
   !$acc declare create( ln_ridging, rn_hstar, rn_porordg, rn_fsnwrdg, rn_fpndrdg, ln_rafting, rn_hraft, rn_craft, rn_fsnwrft, rn_fpndrft )



   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ice_dyn_rdgrft_alloc()
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_dyn_rdgrft_alloc ***
      !!-------------------------------------------------------------------
      IF( ll_diag_rdg ) THEN
         ALLOCATE( closing_net(jpi,jpj) , opning(jpi,jpj)    , closing_gross(jpi,jpj),                   &
            &      apartf(jpi,jpj,0:jpl), hrmin (jpi,jpj,jpl), hraft(jpi,jpj,jpl)    , aridge(jpi,jpj,jpl), &
            &      hrmax (jpi,jpj,jpl)  , hrexp (jpi,jpj,jpl), hi_hrdg(jpi,jpj,jpl)  , araft(jpi,jpj,jpl) , &
            &      airdg1(jpi,jpj)      , airft1(jpi,jpj)    , airdg2(jpi,jpj)       , airft2(jpi,jpj)    , &
                                ! diagnostics
            &      opning_2d(jpi,jpj) , dairdg1dt(jpi,jpj) , dairft1dt(jpi,jpj) , dairdg2dt(jpi,jpj) , dairft2dt(jpi,jpj) , &
                                !
            &      STAT=ice_dyn_rdgrft_alloc )
      ELSE
         ALLOCATE( closing_net(jpi,jpj) , opning(jpi,jpj)    , closing_gross(jpi,jpj),                   &
            &      apartf(jpi,jpj,0:jpl), hrmin (jpi,jpj,jpl), hraft(jpi,jpj,jpl)    , aridge(jpi,jpj,jpl), &
            &      hrmax (jpi,jpj,jpl)  , hrexp (jpi,jpj,jpl), hi_hrdg(jpi,jpj,jpl)  , araft(jpi,jpj,jpl) , &
            &      airdg1(jpi,jpj)      , airft1(jpi,jpj)    , airdg2(jpi,jpj)       , airft2(jpi,jpj)    , &
            &      STAT=ice_dyn_rdgrft_alloc )
# if defined _OPENACC
         PRINT *, ' * info GPU: icedyn_rdgrft() => adding ridging/rafting arrays to memory!'
         PRINT *, '            => closing_net, opning, closing_gross, apartf, hrmin, hraft, aridge'
         PRINT *, '            => hrmax, hrexp, hi_hrdg, araft, airdg1, airft1, airdg2, airft2'
         !$acc enter data copyin( closing_net, opning, closing_gross, apartf, hrmin, hraft, aridge, hrmax, hrexp, hi_hrdg, araft, airdg1, airft1, airdg2, airft2 )
# endif
      ENDIF

      CALL mpp_sum ( 'icedyn_rdgrft', ice_dyn_rdgrft_alloc )
      IF( ice_dyn_rdgrft_alloc /= 0 )   CALL ctl_stop( 'STOP',  'ice_dyn_rdgrft_alloc: failed to allocate arrays'  )
      !
   END FUNCTION ice_dyn_rdgrft_alloc


   SUBROUTINE ice_dyn_rdgrft( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_dyn_rdgrft ***
      !!
      !! ** Purpose :   computes the mechanical redistribution of ice thickness
      !!
      !! ** Method  :   Steps :
      !!       0) Identify grid cells with ice
      !!       1) Calculate closing rate, divergence and opening
      !!       2) Identify grid cells with ridging
      !!       3) Start ridging iterations
      !!          - prep = ridged and rafted ice + closing_gross
      !!          - shift = move ice from one category to another
      !!
      !! ** Details
      !!    step1: The net rate of closing is due to convergence and shear, based on Flato and Hibler (1995).
      !!           The energy dissipation rate is equal to the net closing rate times the ice strength.
      !!
      !!    step3: The gross closing rate is equal to the first two terms (open
      !!           water closing and thin ice ridging) without the third term
      !!           (thick, newly ridged ice).
      !!
      !! References :   Flato, G. M., and W. D. Hibler III, 1995, JGR, 100, 18,611-18,626.
      !!                Hibler, W. D. III, 1980, MWR, 108, 1943-1973, 1980.
      !!                Rothrock, D. A., 1975: JGR, 80, 4514-4519.
      !!                Thorndike et al., 1975, JGR, 80, 4501-4513.
      !!                Bitz et al., JGR, 2001
      !!                Amundrud and Melling, JGR 2005
      !!                Babko et al., JGR 2002
      !!
      !!     This routine is based on CICE code and authors William H. Lipscomb,
      !!     and Elizabeth C. Hunke, LANL are gratefully acknowledged
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl                 ! dummy loop index
      INTEGER  ::   iter, iterate_ridging      ! local integer
      INTEGER  ::   k_np_ice, k_np_rdg         ! n. of points of domain with sea-ice
      REAL(wp) ::   zfac, zsum                 ! local scalar
      !
      LOGICAL,  DIMENSION(jpi,jpj) ::   ll_ridge ! To determine whether ridging is required or not within local domain
      LOGICAL,  DIMENSION(jpi,jpj) ::   ll_ice_present
      REAL(wp), DIMENSION(jpi,jpj) ::   zdivu    ! divu_i
      !!
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk      ! Temporary array for ice presence mask  ! used only when ll_diag_rdg !
      !
      INTEGER, PARAMETER ::   jp_itermax = 20
      !!-------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_rdgrft')
      !$acc data present( a_i,at_i,ato_i,e_i,e_s,v_i,v_s,s_i,sv_i,szv_i,divu_i,delta_i,hi_max,oa_i,hfx_dyn,sfx_bri,sfx_dyn,wfx_dyn,wfx_snw_dyn )
      !$acc data present( sst_s, sss_s, airdg1,airdg2,airft1,airft2,apartf,opning,closing_gross,closing_net )
      !$acc data create( ll_ridge, ll_ice_present, zdivu )

      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icediachk )   CALL ice_cons2D  (0, 'icedyn_rdgrft',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'ice_dyn_rdgrft: ice ridging and rafting'
         IF(lwp) WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF
      !PRINT *, ' *** LOLO entering `ice_dyn_rdgrft` kt =', kt, ' (ll_diag_rdg =',ll_diag_rdg,')'

      ! Initialise ridging diagnostics if required
      IF( ll_diag_rdg ) THEN
# if defined _OPENACC
         CALL ctl_stop( 'STOP',  'icedyn_rdgrft: `ll_diag_rdg=T` not supported yet on GPU!'  )
# endif
         opning_2d(:,:) = 0.0_wp
         dairdg1dt(:,:) = 0.0_wp ; dairft1dt(:,:) = 0.0_wp
         dairdg2dt(:,:) = 0.0_wp ; dairft2dt(:,:) = 0.0_wp
         zmsk(:,:) = MERGE( 1._wp, 0._wp, at_i(:,:) >= epsi10 ) ! 1 if ice, 0 if no ice
      ENDIF

      !--------------------------------
      ! 0) Identify grid cells with ice
      !--------------------------------

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
            ! Initialise logical array for ice points
            ll_ice_present(ji,jj) = .FALSE.
            ! Initialise logical array for ice points where ridging occurs
            ll_ridge(ji,jj) = .FALSE.
         END DO
      END DO
      !$acc end parallel loop


      k_np_ice = 0
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF( at_i(ji,jj) > epsi10 ) THEN
               k_np_ice = k_np_ice + 1
               ll_ice_present(ji,jj) = .TRUE.
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !--------------------------------------------------------
      ! 1) Dynamical inputs (closing rate, divergence, opening)
      !--------------------------------------------------------
      IF( k_np_ice > 0 ) THEN

         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               IF(ll_ice_present(ji,jj)) THEN
                  zdivu(ji,jj) = divu_i(ji,jj)
                  ! closing_net = rate at which open water area is removed + ice area removed by ridging
                  !                                                        - ice area added in new ridges
                  closing_net(ji,jj) = rn_csrdg * 0.5_wp * ( delta_i(ji,jj) - ABS( zdivu(ji,jj) ) ) - MIN( zdivu(ji,jj), 0._wp )
                  !
                  IF( zdivu(ji,jj) < 0._wp )   closing_net(ji,jj) = MAX( closing_net(ji,jj), -zdivu(ji,jj) )   ! make sure the closing rate is large enough
                  !                                                                                ! to give asum = 1.0 after ridging
                  ! Opening rate (non-negative) that will give asum = 1.0 after ridging.
                  opning(ji,jj) = closing_net(ji,jj) + zdivu(ji,jj)
               ENDIF
            END DO
         END DO
         !$acc end parallel loop


         !------------------------------------
         ! 2) Identify grid cells with ridging
         !------------------------------------
         CALL rdgrft_prep( a_i, v_i, ato_i, ll_ice_present, &                                                   ! <<== in
            &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, pclosing_gross=closing_gross, popning=opning, & ! ==>> out
            &              pclosing_net=closing_net )                                                           ! <<== in

         k_np_rdg = 0
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               IF(ll_ice_present(ji,jj))  THEN
                  zsum = 0._wp
                  !$acc loop seq
                  DO jl = 1, jpl
                     zsum = zsum + apartf(ji,jj,jl)
                  END DO
                  IF( zsum > 0._wp .AND. closing_gross(ji,jj) > 0._wp ) THEN
                     k_np_rdg = k_np_rdg + 1
                     ll_ridge(ji,jj)= .TRUE. ! Set to true only if ridging is required
                  ENDIF
               ENDIF
            END DO
         END DO
         !$acc end parallel loop

      ENDIF

      !-----------------
      ! 3) Start ridging
      !-----------------
      IF( k_np_rdg > 0 ) THEN

         iter            = 1
         iterate_ridging = 1
         !                                                        !----------------------!
         DO WHILE( iterate_ridging > 0 .AND. iter < jp_itermax )  !  ridging iterations  !
            !                                                     !----------------------!
            ! Calculate participation function (apartf)
            !       and transfer      function
            !       and closing_gross (+correction on opening)

            CALL rdgrft_prep( a_i, v_i, ato_i, ll_ridge, &                       ! <<== in
               &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, pclosing_gross=closing_gross, popning=opning, & ! ==>> out
               &              pclosing_net=closing_net )                                                           ! <<== in

            ! Redistribute area, volume, and energy between categories
            CALL rdgrft_shift( apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, closing_gross, opning, & ! <<== in
               &               sst_s, sss_s, ll_ridge )                                       ! <<== in


            ! Do we keep on iterating?
            !-------------------------
            ! Check whether a_i + ato_i = 0
            ! If not, because the closing and opening rates were reduced above, ridge again with new rates

            iterate_ridging = 0
            !$acc parallel loop collapse(2)
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF(ll_ridge(ji,jj) ) THEN
                     zsum = 0._wp
                     !$acc loop seq
                     DO jl = 1, jpl
                        zsum = zsum + a_i(ji,jj,jl)
                     END DO
                     !
                     zfac = 1._wp - ( ato_i(ji,jj) + zsum )
                     IF( ABS( zfac ) < epsi10 ) THEN
                        closing_net(ji,jj) = 0._wp
                        opning     (ji,jj) = 0._wp
                        ato_i      (ji,jj) = MAX( 0._wp, 1._wp - zsum )
                     ELSE
                        iterate_ridging  = iterate_ridging + 1
                        zdivu      (ji,jj) = zfac * r1_Dt_ice
                        closing_net(ji,jj) = MAX( 0._wp, -zdivu(ji,jj) )
                        opning     (ji,jj) = MAX( 0._wp,  zdivu(ji,jj) )
                     ENDIF
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop
            !
            iter = iter + 1
            IF( iter  >  jp_itermax )    CALL ctl_stop( 'STOP',  'icedyn_rdgrft: non-converging ridging scheme'  )

         END DO !DO WHILE( iterate_ridging > 0 .AND. iter < jp_itermax )

      ENDIF !IF( k_np_rdg > 0 )

      IF( ll_diag_rdg ) THEN
         !         ! --- Ridging diagnostics --- !
         ! Update diagnostic arrays
         opning_2d(:,:) = opning(:,:)
         dairdg1dt(:,:) = airdg1(:,:) * r1_Dt_ice
         dairft1dt(:,:) = airft1(:,:) * r1_Dt_ice
         dairdg2dt(:,:) = airdg2(:,:) * r1_Dt_ice
         dairft2dt(:,:) = airft2(:,:) * r1_Dt_ice
         CALL iom_put( 'lead_open', opning_2d(:,:) * zmsk(:,:) )  ! Lead area opening rate
         CALL iom_put( 'rdg_loss',  dairdg1dt(:,:) * zmsk(:,:) )  ! Ridging ice area loss rate
         CALL iom_put( 'rft_loss',  dairft1dt(:,:) * zmsk(:,:) )  ! Rafting ice area loss rate
         CALL iom_put( 'rdg_gain',  dairdg2dt(:,:) * zmsk(:,:) )  ! New ridged ice area gain rate
         CALL iom_put( 'rft_gain',  dairft2dt(:,:) * zmsk(:,:) )  ! New rafted ice area gain rate
      ENDIF

      ! clem: those fields must be updated on the halos: ato_i, a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_i, e_s, szv_i

      ! clem: I think we can comment this line but I am not sure it does not change results

      ! controls
      !IF( sn_cfctl%l_prtctl )   CALL ice_prt3D('icedyn_rdgrft')                                                           ! prints
      !IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn rdgrft - ')                             ! prints
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'icedyn_rdgrft',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      !$acc end data
      !$acc end data
      !$acc end data
      !PRINT *, ' *** LOLO exiting `ice_dyn_rdgrft` kt =', kt ;      PRINT *, ''
      IF( ln_timing    )   CALL timing_stop ('icedyn_rdgrft')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_rdgrft


   SUBROUTINE rdgrft_prep( pa_i, pv_i, pato_i, ll_ice_present, &                                                          ! <<== in
      &                    papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp, pclosing_gross, popning, & ! ==>> out
      &                    pclosing_net )                                                                                 ! <<== in
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE rdgrft_prep ***
      !!
      !! ** Purpose :   preparation for ridging calculations
      !!
      !! ** Method  :   Compute the thickness distribution of the ice and open water
      !!                participating in ridging and of the resulting ridges.
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl),   INTENT(in   ) ::   pa_i, pv_i
      REAL(wp), DIMENSION(jpi,jpj),       INTENT(in   ) ::   pato_i
      LOGICAL,  DIMENSION(jpi,jpj),       INTENT(in)    ::   ll_ice_present
      REAL(wp), DIMENSION(jpi,jpj,0:jpl), INTENT(  out) ::   papartf
      REAL(wp), DIMENSION(jpi,jpj,jpl),   INTENT(  out) ::   paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp
      REAL(wp), DIMENSION(jpi,jpj),       INTENT(inout), OPTIONAL ::   pclosing_gross, popning
      REAL(wp), DIMENSION(jpi,jpj),       INTENT(in   ), OPTIONAL ::   pclosing_net
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jm              ! dummy loop indices
      REAL(wp) ::   z1_gstar, z1_astar, zhmean, zfac, zAt, zAl, z1_hi_hrft   ! local scalar
      LOGICAL  ::   l_do_closing_net
      REAL(wp) ::   zasum, z1_asum              ! sum of a_i+ato_i and inverse
      REAL(wp) ::   zaksum                      ! normalisation factor
      REAL(wp), DIMENSION(jpl)    ::   zhi      ! ice thickness
      REAL(wp), DIMENSION(-1:jpl) ::   zGsum    ! zGsum(n) = sum of areas in categories 0 to n
      !--------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('rdgrft_prep')
      !$acc data present( ll_ice_present, papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp ) create( zhi, zGsum )

      l_do_closing_net = ( PRESENT(pclosing_gross) .AND. PRESENT(popning) .AND. PRESENT(pclosing_net) )

      z1_gstar = 1._wp / rn_gstar
      z1_astar = 1._wp / rn_astar
      z1_hi_hrft = 1._wp / hi_hrft

      !$acc parallel loop collapse(2) private( zhi, zGsum )
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF(ll_ice_present(ji,jj)) THEN

               ! Ice thickness needed for rafting
               ! In single precision there were floating point invalids due a sqrt of zhi which happens to have negative values
               ! To solve that an extra check about the value of pv_i was added.
               ! Although adding this condition is safe, the double definition (one for single other for double) has been kept to preserve the results of the sette test.
               zAt = 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  zAt     = zAt + pa_i(ji,jj,jl)
                  zhi(jl) = pv_i(ji,jj,jl) / MAX( pa_i(ji,jj,jl), epsi10 )
               END DO

               ! 1) Participation function (apartf): a(h) = b(h).g(h)
               !-----------------------------------------------------------------
               ! Compute the participation function = total area lost due to ridging/closing
               ! This is analogous to
               !   a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
               !   assuming b(h) = (2/Gstar) * (1 - G(h)/Gstar).
               !
               ! apartf = integrating b(h)g(h) between the category boundaries
               ! apartf is always >= 0 and SUM(apartf(0:jpl))=1
               !-----------------------------------------------------------------
               !
               ! Compute total area of ice plus open water.
               ! This is in general not equal to one because of divergence during transport
               zasum = pato_i(ji,jj) + zAt
               z1_asum = 1._wp / MAX( zasum, epsi10 )
               !
               ! Compute cumulative thickness distribution function
               ! Compute the cumulative thickness distribution function zGsum,
               ! where zGsum(n) is the fractional area in categories 0 to n.
               ! initial value (in h = 0) = open water area
               zGsum(-1) = 0._wp
               zGsum(0 ) = pato_i(ji,jj) * z1_asum
               !$acc loop seq
               DO jl = 1, jpl
                  zAl = 0._wp
                  !$acc loop seq
                  DO jm = 1, jl
                     zAl = zAl + pa_i(ji,jj,jm)  ! sum(1:jl) is correct (and not jpl)
                  END DO
                  zGsum(jl) = ( pato_i(ji,jj) + zAl ) * z1_asum
               END DO

               IF( ln_partf_lin ) THEN          !--- Linear formulation (Thorndike et al., 1975)
                  !$acc loop seq
                  DO jl = 0, jpl
                     IF    ( zGsum(jl)   < rn_gstar ) THEN
                        papartf(ji,jj,jl) = z1_gstar * ( zGsum(jl) - zGsum(jl-1) ) * &
                           &                       ( 2._wp - ( zGsum(jl-1) + zGsum(jl) ) * z1_gstar )
                     ELSEIF( zGsum(jl-1) < rn_gstar ) THEN
                        papartf(ji,jj,jl) = z1_gstar * ( rn_gstar     - zGsum(jl-1) ) *  &
                           &                       ( 2._wp - ( zGsum(jl-1) + rn_gstar     ) * z1_gstar )
                     ELSE
                        papartf(ji,jj,jl) = 0._wp
                     ENDIF
                  END DO
                  !
               ELSEIF( ln_partf_exp ) THEN      !--- Exponential, more stable formulation (Lipscomb et al, 2007)
                  !$acc loop seq
                  DO jl = 0, jpl
                     zfac = 1._wp / ( 1._wp - EXP(-z1_astar) )
                     zGsum(-1) = EXP( -zGsum(-1) * z1_astar ) * zfac
                     zGsum(jl) = EXP( -zGsum(jl) * z1_astar ) * zfac
                     papartf(ji,jj,jl) = zGsum(jl-1) - zGsum(jl)
                  END DO
                  !
               ENDIF

               !                                !--- Ridging and rafting participation concentrations
               IF( ln_rafting .AND. ln_ridging ) THEN             !- ridging & rafting
                  !$acc loop seq
                  DO jl = 1, jpl
                     paridge(ji,jj,jl) = ( 1._wp + TANH ( rn_craft * ( zhi(jl) - rn_hraft ) ) ) * 0.5_wp * papartf(ji,jj,jl)
                     paraft (ji,jj,jl) = papartf(ji,jj,jl) - paridge(ji,jj,jl)
                  END DO
                  !
               ELSEIF( ln_ridging .AND. .NOT. ln_rafting ) THEN   !- ridging alone
                  !$acc loop seq
                  DO jl = 1, jpl
                     paridge(ji,jj,jl) = papartf(ji,jj,jl)
                     paraft (ji,jj,jl) = 0._wp
                  END DO
                  !
               ELSEIF( ln_rafting .AND. .NOT. ln_ridging ) THEN   !- rafting alone
                  !$acc loop seq
                  DO jl = 1, jpl
                     paridge(ji,jj,jl) = 0._wp
                     paraft (ji,jj,jl) = papartf(ji,jj,jl)
                  END DO
                  !
               ELSE                                               !- no ridging & no rafting
                  !$acc loop seq
                  DO jl = 1, jpl
                     paridge(ji,jj,jl) = 0._wp
                     paraft (ji,jj,jl) = 0._wp
                  END DO
                  !
               ENDIF

               ! 2) Transfer function
               !-----------------------------------------------------------------
               ! If assuming ridged ice is uniformly distributed between hrmin and
               ! hrmax (ln_distf_lin):
               !
               ! Compute max and min ridged ice thickness for each ridging category.
               !
               ! This parameterization is a modified version of Hibler (1980).
               ! The mean ridging thickness, zhmean, is proportional to hi^(0.5)
               !  and for very thick ridging ice must be >= hrdg_hi_min*hi
               !
               ! The minimum ridging thickness, hrmin, is equal to 2*hi
               !  (i.e., rafting) and for very thick ridging ice is
               !  constrained by hrmin <= (zhmean + hi)/2.
               !
               ! The maximum ridging thickness, hrmax, is determined by zhmean and hrmin.
               !
               ! These modifications have the effect of reducing the ice strength
               ! (relative to the Hibler formulation) when very thick ice is ridging.
               !
               !-----------------------------------------------------------------
               ! If assuming ridged ice ITD is a negative exponential
               ! (ln_distf_exp) and following CICE implementation:
               !
               !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin
               !
               ! where hrmin is the minimum thickness of ridging ice and
               ! hrexp is the e-folding thickness.
               !
               ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
               ! That is, the minimum ridge thickness results from rafting,
               !  unless the ice is thicker than maxraft.
               !
               ! Also, assume that hrexp = mu_rdg*sqrt(hi).
               ! The parameter mu_rdg is tuned to give e-folding scales mostly
               !  in the range 2-4 m as observed by upward-looking sonar.
               !
               ! Values of mu_rdg in the right column give ice strengths
               !  roughly equal to values of Hstar in the left column
               !  (within ~10 kN/m for typical ITDs):
               !
               !   Hstar     mu_rdg
               !
               !     25        3.0
               !     50        4.0
               !     75        5.0
               !    100        6.0
               !
               ! zaksum = net area removed/ total area participating
               ! where total area participating = area of ice that ridges
               !         net area removed = total area participating - area of new ridges
               !-----------------------------------------------------------------
               zaksum = papartf(ji,jj,0)
               !$acc loop seq
               DO jl = 1, jpl
                  !
                  IF( papartf(ji,jj,jl) > 0._wp ) THEN
                     zhmean        = MAX( SQRT( rn_hstar * zhi(jl) ), zhi(jl) * hrdg_hi_min )
                     phrmin(ji,jj,jl) = MIN( 2._wp * zhi(jl), 0.5_wp * ( zhmean + zhi(jl) ) )
                     phraft(ji,jj,jl) = zhi(jl) * z1_hi_hrft
                     !
                     IF( ln_distf_lin ) THEN
                        phrmax  (ji,jj,jl) = 2._wp * zhmean - phrmin(ji,jj,jl)
                        phi_hrdg(ji,jj,jl) = zhi(jl) / MAX( zhmean, epsi20 )
                     ELSEIF( ln_distf_exp ) THEN
                        phrexp  (ji,jj,jl) = rn_murdg * SQRT( zhi(jl) )
                        phi_hrdg(ji,jj,jl) = zhi(jl) / MAX( epsi20, phrmin(ji,jj,jl) + phrexp(ji,jj,jl) )
                     ENDIF
                     !
                     ! Normalization factor : zaksum, ensures mass conservation
                     zaksum = zaksum + paridge(ji,jj,jl) * ( 1._wp - phi_hrdg(ji,jj,jl) )    &
                        &                    + paraft (ji,jj,jl) * ( 1._wp - hi_hrft )
                  ELSE
                     phrmin  (ji,jj,jl) = 0._wp
                     phrmax  (ji,jj,jl) = 0._wp
                     phrexp  (ji,jj,jl) = 0._wp
                     phraft  (ji,jj,jl) = 0._wp
                     phi_hrdg(ji,jj,jl) = 1._wp
                  ENDIF
                  !
               END DO

               IF( l_do_closing_net ) THEN
                  !
                  ! 3) closing_gross
                  !-----------------
                  ! Based on the ITD of ridging and ridged ice, convert the net closing rate to a gross closing rate.
                  ! NOTE: 0 < aksum <= 1
                  IF(zaksum > epsi10 ) THEN
                     pclosing_gross(ji,jj) = pclosing_net(ji,jj) / zaksum
                  ELSE
                     pclosing_gross(ji,jj) = 0._wp
                  ENDIF

                  ! correction to closing rate if excessive ice removal
                  !----------------------------------------------------
                  ! Reduce the closing rate if more than 100% of any ice category would be removed
                  ! Reduce the opening rate in proportion
                  !$acc loop seq
                  DO jl = 1, jpl
                     zfac = papartf(ji,jj,jl) * pclosing_gross(ji,jj) * rDt_ice
                     IF( zfac > pa_i(ji,jj,jl) .AND. papartf(ji,jj,jl) /= 0._wp ) THEN
                        pclosing_gross(ji,jj) = pa_i(ji,jj,jl) / papartf(ji,jj,jl) * r1_Dt_ice
                     ENDIF
                  END DO

                  ! 4) correction to opening if excessive open water removal
                  !---------------------------------------------------------
                  ! Reduce the closing rate if more than 100% of the open water would be removed
                  ! Reduce the opening rate in proportion
                  zfac = pato_i(ji,jj) + ( popning(ji,jj) - papartf(ji,jj,0) * pclosing_gross(ji,jj) ) * rDt_ice
                  IF( zfac < 0._wp ) THEN           ! would lead to negative ato_i
                     popning(ji,jj) = papartf(ji,jj,0) * pclosing_gross(ji,jj) - pato_i(ji,jj) * r1_Dt_ice
                  ELSEIF( zfac > zasum ) THEN   ! would lead to ato_i > asum
                     popning(ji,jj) = papartf(ji,jj,0) * pclosing_gross(ji,jj) + ( zasum - pato_i(ji,jj) ) * r1_Dt_ice
                  ENDIF

               ENDIF !IF( l_do_closing_net ) THEN

            ENDIF !IF(ll_ice_present(ji,jj))
         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('rdgrft_prep')
   END SUBROUTINE rdgrft_prep


   SUBROUTINE rdgrft_shift( papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp, pclosing_gross, popning, & ! <<== in
      &                     psst, psss, ll_ice_present )                                                                   ! <<== in
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE rdgrft_shift ***
      !!
      !! ** Purpose :   shift ridging ice among thickness categories of ice thickness
      !!
      !! ** Method  :   Remove area, volume, and energy from each ridging category
      !!                and add to thicker ice categories.
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,0:jpl)     , INTENT(in) ::   papartf
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in) ::   paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in) ::   pclosing_gross, popning
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in) ::   psst, psss
      LOGICAL,  DIMENSION(jpi,jpj),            INTENT(in) ::   ll_ice_present
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jl1, jl2, jk       ! dummy loop indices
      REAL(wp) ::   hL, hR, farea              ! left and right limits of integration and new area going to jl2
      REAL(wp) ::   expL, expR                 ! exponentials involving hL, hR
      REAL(wp) ::   vsw                        ! vol of water trapped into ridges
      REAL(wp) ::   afrdg, afrft               ! fraction of category area ridged/rafted
      REAL(wp) ::   zoirdg, zvirdg, zvsrdg !, zaprdg, zvprdg, zvlrdg  ! area etc of new ridges
      REAL(wp) ::   zoirft, zvirft, zvsrft !, zaprft, zvprft, vzlrft  ! area etc of rafted ice
      !
      REAL(wp) ::   ersw             ! enthalpy of water trapped into ridges
      REAL(wp) ::   zswitch, fvol    ! new ridge volume going to jl2
      REAL(wp) ::   z1_ai            ! 1 / a
      REAL(wp) ::   zvti             ! sum(v_i)
      !
      REAL(wp), DIMENSION(nlay_i) ::   zsirdg, zsirft
      REAL(wp), DIMENSION(nlay_s) ::   zesrft     ! snow energy of rafting ice
      REAL(wp), DIMENSION(nlay_i) ::   zeirft     ! ice  energy of rafting ice
      REAL(wp), DIMENSION(nlay_s) ::   zesrdg     ! enth*volume of new ridges
      REAL(wp), DIMENSION(nlay_i) ::   zeirdg     ! enth*volume of new ridges
      !
      INTEGER ::   itest_rdg, itest_rft   ! test for conservation
      LOGICAL ::   ll_shift         ! logical for doing calculation or not
      !!-------------------------------------------------------------------
      IF( ln_timing )  CALL timing_start('rdgrft_shift')
      !$acc data present( papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp, pclosing_gross, popning, psst, psss, ll_ice_present )
      !$acc data  create( zsirdg, zsirft, zesrft, zeirft, zesrdg, zeirdg )

      !$acc parallel loop collapse(2) private( zsirdg, zsirft, zesrft, zeirft, zesrdg, zeirdg )
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ll_shift = .FALSE.
            zvti = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               zvti = zvti + v_i(ji,jj,jl)   ! total ice volume
            END DO


            IF(ll_ice_present(ji,jj)) THEN

               ! 1) Change in open water area due to closing and opening
               !--------------------------------------------------------
               ato_i(ji,jj) = MAX( 0._wp, ato_i(ji,jj) + ( popning(ji,jj) - papartf(ji,jj,0) * pclosing_gross(ji,jj) ) * rDt_ice )
               !$acc loop seq
               DO jl1 = 1, jpl

                  ! set logical to true when ridging
                  IF( papartf(ji,jj,jl1) > 0._wp .AND. pclosing_gross(ji,jj) > 0._wp ) THEN
                     ll_shift = .TRUE.
                  ELSE
                     ll_shift = .FALSE.
                  ENDIF

                  IF( ll_shift ) THEN   ! only if ice is ridging

                     IF( a_i(ji,jj,jl1) > epsi10 ) THEN
                        z1_ai = 1._wp / a_i(ji,jj,jl1)
                     ELSE
                        z1_ai = 0._wp
                     ENDIF

                     ! area of ridging / rafting ice (airdg1) and of new ridge (airdg2)
                     airdg1(ji,jj) = paridge(ji,jj,jl1) * pclosing_gross(ji,jj) * rDt_ice
                     airft1(ji,jj) = paraft (ji,jj,jl1) * pclosing_gross(ji,jj) * rDt_ice

                     airdg2(ji,jj) = airdg1(ji,jj) * phi_hrdg(ji,jj,jl1)
                     airft2(ji,jj) = airft1(ji,jj) *  hi_hrft

                     ! ridging /rafting fractions
                     afrdg = airdg1(ji,jj) * z1_ai
                     afrft = airft1(ji,jj) * z1_ai

                     ! volume and enthalpy (J/m2, >0) of seawater trapped into ridges
                     IF    ( zvti <= 10. ) THEN
                        vsw = v_i(ji,jj,jl1) * afrdg * rn_porordg                                           ! v <= 10m then porosity = rn_porordg
                     ELSEIF( zvti >= 20. ) THEN
                        vsw = 0._wp                                                                         ! v >= 20m then porosity = 0
                     ELSE
                        vsw = v_i(ji,jj,jl1) * afrdg * rn_porordg * MAX( 0._wp, 2._wp - 0.1_wp * zvti ) ! v > 10m and v < 20m then porosity = linear transition to 0
                     ENDIF
                     ersw = -rhoi * vsw * rcp * psst(ji,jj)   ! clem: if sst>0, then ersw <0 (is that possible?)

                     ! volume etc of ridging / rafting ice and new ridges (vi, vs, sm, oi, es, ei)
                     zvirdg = v_i(ji,jj,jl1)    * afrdg + vsw
                     zvsrdg = v_s(ji,jj,jl1)    * afrdg
                     zoirdg = oa_i(ji,jj,jl1)    * afrdg * phi_hrdg(ji,jj,jl1)

                     zvirft = v_i(ji,jj,jl1)    * afrft
                     zvsrft = v_s(ji,jj,jl1)    * afrft
                     zoirft = oa_i(ji,jj,jl1)    * afrft * hi_hrft

                     !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   zaprdg = a_ip(ji,jj,jl1)    * afrdg * phi_hrdg(ji,jj,jl1)
                     !   zvprdg = v_ip(ji,jj,jl1)    * afrdg
                     !   zaprft = a_ip(ji,jj,jl1)    * afrft *  hi_hrft
                     !   zvprft = v_ip(ji,jj,jl1)    * afrft
                     !   IF( ln_pnd_lids ) THEN
                     !      vlrdg(ji,jj) = v_il(ji,jj,jl1) * afrdg
                     !      vlrft(ji,jj) = v_il(ji,jj,jl1) * afrft
                     !   ENDIF
                     !ENDIF
                     !$acc loop seq
                     DO jk=1, nlay_s
                        zesrdg(jk) = e_s(ji,jj,jk,jl1) * afrdg
                        zesrft(jk) = e_s(ji,jj,jk,jl1) * afrft
                     END DO
                     !$acc loop seq
                     DO jk=1, nlay_i
                        zeirdg(jk) = e_i(ji,jj,jk,jl1) * afrdg + ersw * r1_nlay_i
                        zeirft(jk) = e_i(ji,jj,jk,jl1) * afrft
                     END DO

                     IF( nn_icesal == 4 ) THEN
                        !$acc loop seq
                        DO jk=1, nlay_i
                           zsirdg(jk) = szv_i(ji,jj,jk,jl1) * afrdg + vsw * psss(ji,jj) * r1_nlay_i
                           zsirft(jk) = szv_i(ji,jj,jk,jl1) * afrft
                        END DO
                     ELSE
                        zsirdg(1) = sv_i(ji, jj, jl1) * afrdg + vsw * psss(ji,jj)
                        zsirft(1) = sv_i(ji, jj, jl1) * afrft
                     ENDIF

                     ! Ice-ocean exchanges associated with ice porosity
                     wfx_dyn(ji,jj) = wfx_dyn(ji,jj) - vsw * rhoi * r1_Dt_ice   ! increase in ice volume due to seawater frozen in voids
                     sfx_dyn(ji,jj) = sfx_dyn(ji,jj) - vsw * psss(ji,jj) * rhoi * r1_Dt_ice
                     hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + ersw * r1_Dt_ice          ! > 0 [W.m-2]

                     ! Put the snow and pond lost by ridging into the ocean
                     !  Note that esrdg > 0; the ocean must cool to melt snow. If the ocean temp = Tf already, new ice must grow.
                     wfx_snw_dyn(ji,jj) = wfx_snw_dyn(ji,jj) + ( rhos * zvsrdg * ( 1._wp - rn_fsnwrdg )   &   ! fresh water source for ocean
                        &                                      + rhos * zvsrft * ( 1._wp - rn_fsnwrft ) ) * r1_Dt_ice
                     !$acc loop seq
                     DO jk = 1, nlay_s
                        hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + ( - zesrdg(jk) * ( 1._wp - rn_fsnwrdg )   &          ! heat sink for ocean (<0, W.m-2)
                           &                                - zesrft(jk) * ( 1._wp - rn_fsnwrft ) ) * r1_Dt_ice
                     END DO

                     !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   wfx_pnd(ji,jj)    = wfx_pnd(ji,jj)   + ( rhow * zvprdg * ( 1._wp - rn_fpndrdg )   &   ! fresh water source for ocean
                     !      &                                   + rhow * zvprft * ( 1._wp - rn_fpndrft ) ) * r1_Dt_ice
                     !   IF( ln_pnd_lids ) THEN
                     !      wfx_pnd(ji,jj) = wfx_pnd(ji,jj)   + ( rhow * vlrdg(ji,jj) * ( 1._wp - rn_fpndrdg )   &   ! fresh water source for ocean
                     !         &                                + rhow * vlrft(ji,jj) * ( 1._wp - rn_fpndrft ) ) * r1_Dt_ice
                     !   ENDIF
                     !ENDIF

                     ! virtual salt flux to keep salinity constant
                     IF( nn_icesal == 1 .OR. nn_icesal == 3 )  THEN
                        zsirdg(1)    = zsirdg(1)    - ( psss(ji,jj) - s_i(ji,jj,jl1) ) * vsw                      ! ridge salinity = s_i
                        sfx_bri(ji,jj) = sfx_bri(ji,jj) + ( psss(ji,jj) - s_i(ji,jj,jl1) ) * vsw * rhoi * r1_Dt_ice   ! put back sss_s into the ocean
                        !                                                                                          ! and get  s_i  from the ocean
                     ENDIF

                     ! Remove area, volume of new ridge to each category jl1
                     !------------------------------------------------------
                     a_i (ji,jj,jl1) = a_i (ji,jj,jl1) - airdg1(ji,jj) - airft1(ji,jj)
                     v_i (ji,jj,jl1) = v_i (ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                     v_s (ji,jj,jl1) = v_s (ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                     oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                     !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   a_ip(ji,jj,jl1)  = a_ip(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     !   v_ip(ji,jj,jl1)  = v_ip(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     !   IF( ln_pnd_lids ) v_il(ji,jj,jl1) = v_il(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     !ENDIF
                     !$acc loop seq
                     DO jk=1, nlay_s
                        e_s(ji,jj,jk,jl1) = e_s(ji,jj,jk,jl1)   * ( 1._wp - afrdg - afrft )
                     END DO
                     !$acc loop seq
                     DO jk=1, nlay_i
                        e_i(ji,jj,jk,jl1) = e_i(ji,jj,jk,jl1)   * ( 1._wp - afrdg - afrft )
                     END DO
                     !
                     IF( nn_icesal == 4 ) THEN
                        !$acc loop seq
                        DO jk=1, nlay_i
                           szv_i(ji,jj,jk,jl1) = szv_i(ji,jj,jk,jl1) * ( 1._wp - afrdg - afrft )
                        END DO
                     ELSE
                        sv_i(ji,jj,jl1) = sv_i(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     ENDIF

                  ENDIF

                  ! 2) compute categories in which ice is added (jl2)
                  !--------------------------------------------------
                  itest_rdg = 0
                  itest_rft = 0
                  !$acc loop seq
                  DO jl2  = 1, jpl

                     !IF( ll_ice_present(ji,jj) .AND. ll_shift ) THEN
                     IF( ll_shift ) THEN

                        ! Compute the fraction of ridged ice area and volume going to thickness category jl2
                        IF( ln_distf_lin ) THEN ! Hibler (1980) linear formulation
                           !
                           IF( phrmin(ji,jj,jl1) <= hi_max(jl2) .AND. phrmax(ji,jj,jl1) > hi_max(jl2-1) ) THEN
                              hL = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                              hR = MIN( phrmax(ji,jj,jl1), hi_max(jl2)   )
                              farea = ( hR      - hL      ) / ( phrmax(ji,jj,jl1)                  - phrmin(ji,jj,jl1)                  )
                              fvol  = ( hR * hR - hL * hL ) / ( phrmax(ji,jj,jl1) * phrmax(ji,jj,jl1) - phrmin(ji,jj,jl1) * phrmin(ji,jj,jl1) )
                              !
                              itest_rdg = 1   ! test for conservation
                           ELSE
                              farea = 0._wp
                              fvol  = 0._wp
                           ENDIF
                           !
                        ELSEIF( ln_distf_exp ) THEN ! Lipscomb et al. (2007) exponential formulation
                           !
                           IF( jl2 < jpl ) THEN
                              !
                              IF( phrmin(ji,jj,jl1) <= hi_max(jl2) ) THEN
                                 hL    = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                                 hR    = hi_max(jl2)
                                 expL  = EXP( -( hL - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                                 expR  = EXP( -( hR - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                                 farea = expL - expR
                                 fvol  = ( ( hL + phrexp(ji,jj,jl1) ) * expL  &
                                    &    - ( hR + phrexp(ji,jj,jl1) ) * expR ) / MAX( epsi20, phrmin(ji,jj,jl1) + phrexp(ji,jj,jl1) )
                              ELSE
                                 farea = 0._wp
                                 fvol  = 0._wp
                              ENDIF
                              !
                           ELSE             ! jl2 = jpl
                              !
                              hL    = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                              expL  = EXP(-( hL - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                              farea = expL
                              fvol  = ( hL + phrexp(ji,jj,jl1) ) * expL / MAX( epsi20, phrmin(ji,jj,jl1) + phrexp(ji,jj,jl1) )
                              !
                           ENDIF            ! jl2 < jpl
                           !
                           itest_rdg = 1   ! test for conservation => clem: I am not sure about that
                           !
                        ENDIF             ! ridge redistribution

                        ! Compute the fraction of rafted ice area and volume going to thickness category jl2
                        IF( phraft(ji,jj,jl1) <= hi_max(jl2) .AND. phraft(ji,jj,jl1) >  hi_max(jl2-1) ) THEN
                           zswitch = 1._wp
                           !
                           itest_rft = 1   ! test for conservation
                        ELSE
                           zswitch = 0._wp
                        ENDIF
                        !
                        ! Patch to ensure perfect conservation if ice thickness goes mad
                        ! Sometimes thickness is larger than hi_max(jpl) because of advection scheme (for very small areas)
                        ! Then ice volume is removed from one category but the ridging/rafting scheme
                        ! does not know where to move it, leading to a conservation issue.
                        IF( itest_rdg == 0 .AND. jl2 == jpl ) THEN
                           farea = 1._wp
                           fvol = 1._wp
                        ENDIF
                        IF( itest_rft == 0 .AND. jl2 == jpl ) zswitch = 1._wp
                        !
                        ! Add area, volume of new ridge to category jl2
                        !----------------------------------------------
                        a_i (ji,jj,jl2) = a_i (ji,jj,jl2) + ( airdg2(ji,jj)              * farea + airft2(ji,jj)              * zswitch )
                        oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ( zoirdg              * farea + zoirft              * zswitch )
                        v_i (ji,jj,jl2) = v_i (ji,jj,jl2) + ( zvirdg              * fvol  + zvirft              * zswitch )
                        v_s (ji,jj,jl2) = v_s (ji,jj,jl2) + ( zvsrdg * rn_fsnwrdg * fvol  + zvsrft * rn_fsnwrft * zswitch )
                        !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                        !   v_ip (ji,jj,jl2) = v_ip(ji,jj,jl2)    + ( zvprdg * rn_fpndrdg * fvol  + zvprft * rn_fpndrft * zswitch )
                        !   a_ip (ji,jj,jl2) = a_ip(ji,jj,jl2)    + ( zaprdg * rn_fpndrdg * farea + zaprft * rn_fpndrft * zswitch )
                        !   IF( ln_pnd_lids ) THEN
                        !      v_il (ji,jj,jl2) = v_il(ji,jj,jl2) + ( vlrdg(ji,jj) * rn_fpndrdg * fvol  + vlrft(ji,jj) * rn_fpndrft * zswitch )
                        !   ENDIF
                        !ENDIF
                        !$acc loop seq
                        DO jk=1, nlay_s
                           e_s(ji,jj,jk,jl2) = e_s(ji,jj,jk,jl2) + ( zesrdg(jk) * rn_fsnwrdg * fvol + zesrft(jk) * rn_fsnwrft * zswitch )
                        END DO
                        !$acc loop seq
                        DO jk=1, nlay_i
                           e_i  (ji,jj,jk,jl2) = e_i  (ji,jj,jk,jl2) + ( zeirdg(jk)              * fvol + zeirft(jk)              * zswitch )
                        END DO
                        IF( nn_icesal == 4 ) THEN
                           !$acc loop seq
                           DO jk=1, nlay_i
                              szv_i(ji,jj,jk,jl2) = szv_i(ji,jj,jk,jl2) + ( zsirdg(jk) * fvol + zsirft(jk) * zswitch )
                           END DO
                        ELSE
                           sv_i (ji,jj,  jl2) = sv_i (ji, jj, jl2) + ( zsirdg(1) * fvol + zsirft(1) * zswitch )
                        ENDIF
                        !
                     ENDIF


                  END DO ! jl2
                  !
               END DO ! jl1


            END IF !IF(ll_ice_present(ji,jj))
         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !$acc end parallel loop

      ! roundoff errors
      !----------------
      ! In case ridging/rafting lead to very small negative values (sometimes it happens)
      !CALL ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, ll_ice_present)
      CALL  ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i,                   e_s, e_i, szv_i, ll_ice_present)

      !$acc end data
      !$acc end data
      IF( ln_timing )  CALL timing_stop('rdgrft_shift')
      !
   END SUBROUTINE rdgrft_shift


   SUBROUTINE ice_strength
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ice_strength ***
      !!
      !! ** Purpose :   computes ice strength used in dynamics routines of ice thickness
      !!
      !! ** Method  :   Compute the strength of the ice pack, defined as the energy (J m-2)
      !!              dissipated per unit area removed from the ice pack under compression,
      !!              and assumed proportional to the change in potential energy caused
      !!              by ridging. Note that ice strength using Hibler's formulation must be
      !!              smoothed.
      !!----------------------------------------------------------------------
      INTEGER             ::   ji, jj, jl  ! dummy loop indices
      INTEGER             ::   kice_p      ! n. of points with ice
      REAL(wp)            ::   z1_3        ! local scalars
      REAL(wp)            ::   zaksum      ! normalisation factor
      REAL(wp), DIMENSION(jpi,jpj)   ::   zistr         ! temporary array used here
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   za_i_cap       ! local capped ice concentration
      !!
      LOGICAL, DIMENSION(jpi,jpj) :: ll_ice_present
      REAL(wp)            ::   zhi, zcp, zv
      REAL(wp)            ::   h2rdg                     ! mean value of h^2 for new ridge
      !
      REAL(wp), PARAMETER ::   zmax_strength = 200.e3_wp ! Max strength for R75 formulation.
      !                                                  ! Richter-Menge and Elder (1998) estimate maximum
      !                                                  ! in Beaufort Sea in wintertime of the order 150 kN/m.
      !
      ! Coon et al. (2007) state that 20 kN/m is ~10% of the maximum compressive strength of isotropic ice, giving max strength of 200 kN/m.
      
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_strength')
      !$acc data present( strength, apartf, araft, aridge, hi_hrdg, hi_hrdg, hraft, hrmin, hrmax, hrexp ) create( ll_ice_present, zistr )

      !LOLOrm:
      ! The 2 following sanity tests should not belong here, it is a consequence of `ice_strength()` being always used by EVP, and using
      ! EVP does not necessarilly mean that the user as selected `ln_dynALL=T`, so for example if EVP is used and `ln_dynRHGADV=T`
      ! (which implies `ln_dynALL=F`), all namelist parameters and allocated arrays needed here in `ice_strength()` are not set/allocated,
      ! because `namdyn_rdgrft` is technically not read...
      ! But I added `IF( ln_dynALL .OR. ln_rhg_EVP )` for calling `ice_dyn_rdgrft_init` in `ice_dyn_init` so we should be fine...
      !IF( nice_str <= 0                      ) CALL ctl_stop( 'ice_strength@icedyn_rdgrft.F90: `nice_str <=0 ` ! Should be set!' )
      !IF( ln_str_R75 .AND. (.NOT. ln_dynALL) ) CALL ctl_stop( 'ice_strength@icedyn_rdgrft.F90: `ln_str_R75=T` & `ln_dynALL=F` !' )
      !LOLOrm.

      kice_p = 0
      ! `at_i` needed for strength
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
            !
            IF( at_i(ji,jj) > epsi10) THEN
               kice_p = kice_p + 1
               ll_ice_present(ji,jj) = .TRUE.
            ELSE
               ll_ice_present(ji,jj) = .FALSE.
            ENDIF
            !
            zistr(ji,jj) = 0._wp
            !
         END DO
      END DO
      !$acc end parallel loop



      SELECT CASE( nice_str )          !--- Set which ice strength is chosen

      CASE ( np_strr75 )           !== Rothrock(1975)'s method ==!
         !$acc data present( apartf, araft, aridge, hi_hrdg, hi_hrdg, hraft, hrmin, hrmax, hrexp ) create( za_i_cap )

         ! this should be defined once for all at the 1st time step
         zcp = 0.5_wp * grav * (rho0-rhoi) * rhoi * r1_rho0   ! proport const for PE
         !
         ! Initialise local capped a_i to zero
         ! Note that if 0 < a_i < epsi10, can end up with zhi=0 but apartf>0 rdgrft_prep; za_i_cap avoids this here
         !$acc parallel loop collapse(3)
         DO jl = 1, jpl
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  za_i_cap(ji,jj,jl) = 0._wp
               END DO
            END DO
         END DO
         !$acc end parallel loop


         IF( kice_p > 0 ) THEN

            ! Cap a_i to avoid zhi in rdgrft_prep going below minimum
            !$acc parallel loop collapse(3)
            DO jl = 1, jpl
               DO jj=Njs0, Nje0
                  DO ji=Nis0, Nie0
                     za_i_cap(ji,jj,jl) = a_i(ji,jj,jl)
                  END DO
               END DO
            END DO
            !$acc end parallel loop

            !$acc parallel loop collapse(2)
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF(ll_ice_present(ji,jj)) THEN
                     !$acc loop seq
                     DO jl = 1, jpl
                        zv = v_i(ji,jj,jl) / a_i(ji,jj,jl)
                        IF( zv < rn_himin )  za_i_cap(ji,jj,jl) = a_i(ji,jj,jl) * zv / rn_himin
                     END DO
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop

            CALL rdgrft_prep( za_i_cap, v_i, ato_i, ll_ice_present, &                ! <<== in
               &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp )             ! ==>> out
            
            z1_3 = 1._wp / 3._wp
            !$acc parallel loop collapse(2)
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF(ll_ice_present(ji,jj)) THEN
                     !
                     zaksum = apartf(ji,jj,0) !clem: aksum should be defined in the header => local to module
                     !$acc loop seq
                     DO jl = 1, jpl
                        IF( apartf(ji,jj,jl) > 0._wp ) THEN
                           zaksum = zaksum + aridge(ji,jj,jl) * ( 1._wp - hi_hrdg(ji,jj,jl) ) &
                              &                          + araft (ji,jj,jl) * ( 1._wp - hi_hrft )
                        ENDIF
                     END DO
                     !
                     !$acc loop seq
                     DO jl = 1, jpl
                        !
                        IF( apartf(ji,jj,jl) > 0._wp ) THEN
                           !
                           IF( ln_distf_lin ) THEN       ! Uniform redistribution of ridged ice
                              h2rdg = z1_3 * ( hrmax(ji,jj,jl) * hrmax(ji,jj,jl) +     & ! (a**3-b**3)/(a-b) = a*a+ab+b*b
                                 &             hrmin(ji,jj,jl) * hrmin(ji,jj,jl) +     &
                                 &             hrmax(ji,jj,jl) * hrmin(ji,jj,jl) )
                              !
                           ELSEIF( ln_distf_exp ) THEN   ! Exponential redistribution of ridged ice
                              h2rdg =          hrmin(ji,jj,jl) * hrmin(ji,jj,jl)   &
                                 &   + 2._wp * hrmin(ji,jj,jl) * hrexp(ji,jj,jl)   &
                                 &   + 2._wp * hrexp(ji,jj,jl) * hrexp(ji,jj,jl)
                           ENDIF
                           !
                           zhi = MERGE( v_i(ji,jj,jl) / MAX( a_i(ji,jj,jl), epsi10 )  ,  0._wp  ,  a_i(ji,jj,jl) > epsi10 )
                                                      
                           ! Make sure ice thickness is not below the minimum
                           ! Do not adjust concentration as don't want strength routine to be able to do this
                           zhi = MAX( zhi, rn_himin )

                           zistr(ji,jj) = zistr(ji,jj) - apartf(ji,jj,jl) * zhi * zhi                  ! PE loss
                           zistr(ji,jj) = zistr(ji,jj) + 2._wp * araft(ji,jj,jl) * zhi * zhi           ! PE gain (rafting)
                           zistr(ji,jj) = zistr(ji,jj) + aridge(ji,jj,jl) * h2rdg *  hi_hrdg(ji,jj,jl)    ! PE gain (ridging)
                           
                        ENDIF
                        !
                     END DO !DO jl = 1, jpl
                     !
                     zistr(ji,jj) = MIN( rn_pe_rdg * zcp * zistr(ji,jj) / zaksum , zmax_strength )   ! Enforce a maximum for R75 strength
                     !
                  ENDIF !IF(ll_ice_present(ji,jj))
               END DO !DO ji=Nis0, Nie0
            END DO !DO jj=Njs0, Nje0
            !$acc end parallel loop
            
         ENDIF !IF( kice_p > 0 )
         !
# if ! defined _OPENACC
         CALL lbc_lnk( 'icedyn_rdgrft', zistr, 'T', 1.0_wp ) ! this call could be removed if calculations were done on the full domain
         !                                                   ! but we decided it is more efficient this way
# endif
         !$acc end data

         
      CASE ( np_strh79 )           !== Hibler(1979)'s method ==!
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zv = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zv = zv + v_i(ji,jj,jl)
               END DO
               zistr(ji,jj) = MERGE( rn_pstar * zv * EXP( -rn_crhg * ( 1._wp - at_i(ji,jj) ) ),  0._wp,  ll_ice_present(ji,jj) )
            END DO
         END DO
         !$acc end parallel loop
         !
      CASE ( np_strcst )           !== Constant strength ==!
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zistr(ji,jj) = MERGE( rn_str,  0._wp,  ll_ice_present(ji,jj) )
            END DO
         END DO
         !$acc end parallel loop
         !
      END SELECT


      IF( ln_str_smooth ) THEN         !--- Spatial smoothing
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               IF( ll_ice_present(ji,jj) ) THEN
                  strength(ji,jj) = ( 4._wp * zistr(ji,jj)              &
                     &                    + ( ( zistr(ji-1,jj) * xmskt(ji-1,jj) + zistr(ji+1,jj) * xmskt(ji+1,jj) ) &
                     &                      + ( zistr(ji,jj-1) * xmskt(ji,jj-1) + zistr(ji,jj+1) * xmskt(ji,jj+1) ) ) &
                     &            ) / ( 4._wp + xmskt(ji-1,jj) + xmskt(ji+1,jj) + xmskt(ji,jj-1) + xmskt(ji,jj+1) )
               ELSE
                  strength(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
# if ! defined _OPENACC
         CALL lbc_lnk( 'icedyn_rdgrft', strength, 'T', 1._wp )
# endif
         !
      ELSE
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               strength(ji,jj) = zistr(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF !IF( ln_str_smooth )
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_strength')
   END SUBROUTINE ice_strength


   SUBROUTINE ice_dyn_rdgrft_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_rdgrft_init ***
      !!
      !! ** Purpose :   Physical constants and parameters linked
      !!                to the mechanical ice redistribution
      !!
      !! ** Method  :   Read the namdyn_rdgrft namelist
      !!                and check the parameters values
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_rdgrft
      !!-------------------------------------------------------------------
      INTEGER :: ios, ioptio                ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_rdgrft/ rn_delta_ecc, ln_str_H79, rn_pstar, rn_crhg, ln_str_R75, rn_pe_rdg, ln_str_CST, rn_str, ln_str_smooth, &
         &                    ln_distf_lin, ln_distf_exp, rn_murdg, rn_csrdg,            &
         &                    ln_partf_lin, rn_gstar, ln_partf_exp, rn_astar,            &
         &                    ln_ridging, rn_hstar, rn_porordg, rn_fsnwrdg, rn_fpndrdg,  &
         &                    ln_rafting, rn_hraft, rn_craft  , rn_fsnwrft, rn_fpndrft
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn_rdgrft)
      READ_NML_CFG(numnam_ice,namdyn_rdgrft)
      IF(lwm) WRITE ( numoni, namdyn_rdgrft )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_rdgrft_init: ice parameters for ridging/rafting '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_rdgrft:'
         WRITE(numout,*) '      eccentricity of yield curve to compute `delta`           rn_delta_ecc = ', rn_delta_ecc
         WRITE(numout,*) '      ice strength parameterization Hibler (1979)              ln_str_H79   = ', ln_str_H79
         WRITE(numout,*) '            1st bulk-rheology parameter                        rn_pstar     = ', rn_pstar
         WRITE(numout,*) '            2nd bulk-rhelogy parameter                         rn_crhg      = ', rn_crhg
         WRITE(numout,*) '      ice strength parameterization Rothrock (1975)            ln_str_R75   = ', ln_str_R75
         WRITE(numout,*) '            coef accounting for frictional dissipation         rn_pe_rdg    = ', rn_pe_rdg
         WRITE(numout,*) '      ice strength parameterization Constant                   ln_str_CST   = ', ln_str_CST
         WRITE(numout,*) '            ice strength value                                 rn_str       = ', rn_str
         WRITE(numout,*) '      spatial smoothing of the strength                        ln_str_smooth= ', ln_str_smooth
         WRITE(numout,*) '      redistribution of ridged ice: linear (Hibler 1980)       ln_distf_lin = ', ln_distf_lin
         WRITE(numout,*) '      redistribution of ridged ice: exponential(Lipscomb 2017) ln_distf_exp = ', ln_distf_exp
         WRITE(numout,*) '            e-folding scale of ridged ice                      rn_murdg     = ', rn_murdg
         WRITE(numout,*) '      Fraction of shear energy contributing to ridging         rn_csrdg     = ', rn_csrdg
         WRITE(numout,*) '      linear ridging participation function                    ln_partf_lin = ', ln_partf_lin
         WRITE(numout,*) '            Fraction of ice coverage contributing to ridging   rn_gstar     = ', rn_gstar
         WRITE(numout,*) '      Exponential ridging participation function               ln_partf_exp = ', ln_partf_exp
         WRITE(numout,*) '            Equivalent to G* for an exponential function       rn_astar     = ', rn_astar
         WRITE(numout,*) '      Ridging of ice sheets or not                             ln_ridging   = ', ln_ridging
         WRITE(numout,*) '            max ridged ice thickness                           rn_hstar     = ', rn_hstar
         WRITE(numout,*) '            Initial porosity of ridges                         rn_porordg   = ', rn_porordg
         WRITE(numout,*) '            Fraction of snow volume conserved during ridging   rn_fsnwrdg   = ', rn_fsnwrdg
         WRITE(numout,*) '            Fraction of pond volume conserved during ridging   rn_fpndrdg   = ', rn_fpndrdg
         WRITE(numout,*) '      Rafting of ice sheets or not                             ln_rafting   = ', ln_rafting
         WRITE(numout,*) '            Parmeter thickness (threshold between ridge-raft)  rn_hraft     = ', rn_hraft
         WRITE(numout,*) '            Rafting hyperbolic tangent coefficient             rn_craft     = ', rn_craft
         WRITE(numout,*) '            Fraction of snow volume conserved during rafting   rn_fsnwrft   = ', rn_fsnwrft
         WRITE(numout,*) '            Fraction of pond volume conserved during rafting   rn_fpndrft   = ', rn_fpndrft
      ENDIF
      !
      ioptio = 0
      IF( ln_str_H79    ) THEN
         ioptio = ioptio + 1
         nice_str = np_strh79
      ENDIF
      IF( ln_str_R75    ) THEN
         ioptio = ioptio + 1
         nice_str = np_strr75
      ENDIF
      IF( ln_str_CST    ) THEN
         ioptio = ioptio + 1
         nice_str = np_strcst
      ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: one and only one ice strength option has to be defined ' )
      !
      ioptio = 0
      IF( ln_distf_lin ) THEN
         ioptio = ioptio + 1
      ENDIF
      IF( ln_distf_exp ) THEN
         ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: choose one and only one redistribution function (ln_distf_lin or ln_distf_exp)' )
      !
      ioptio = 0
      IF( ln_partf_lin ) THEN
         ioptio = ioptio + 1
      ENDIF
      IF( ln_partf_exp ) THEN
         ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: choose one and only one participation function (ln_partf_lin or ln_partf_exp)' )
      !
      IF( .NOT. ln_icethd ) THEN
         rn_porordg = 0._wp
         rn_fsnwrdg = 1._wp ; rn_fsnwrft = 1._wp
         rn_fpndrdg = 1._wp ; rn_fpndrft = 1._wp
         IF( lwp ) THEN
            WRITE(numout,*) '      ==> only ice dynamics is activated, thus some parameters must be changed'
            WRITE(numout,*) '            rn_porordg   = ', rn_porordg
            WRITE(numout,*) '            rn_fsnwrdg   = ', rn_fsnwrdg
            WRITE(numout,*) '            rn_fpndrdg   = ', rn_fpndrdg
            WRITE(numout,*) '            rn_fsnwrft   = ', rn_fsnwrft
            WRITE(numout,*) '            rn_fpndrft   = ', rn_fpndrft
         ENDIF
      ENDIF
      !
      ! diagnostics
      !IF( iom_use('lead_open') .OR. iom_use('rdg_loss') .OR. iom_use('rft_loss') .OR. &
      !   &                          iom_use('rdg_gain') .OR. iom_use('rft_gain') ) THEN
      !   ll_diag_rdg = .TRUE.
      !ELSE
      ll_diag_rdg = .FALSE.
      !ENDIF
      !                              ! allocate arrays
      IF( ln_dynALL .OR. (ln_rhg_EVP .AND. ln_str_R75) ) THEN
         !  => LB: because one might want to use the rheology with EVP and `ice_strength` using "Rothrock_75" without necessarily
         !         using the "ridging/rafting", that the case when `ln_rhg_EVP=T` & `ln_dynRHGADV=T` !!!
         !        `ice_dyn_rdgrft_init` is called `IF( ln_dynALL .OR. ln_rhg_EVP )` in `ice_dyn_init@icedyn.F90`...
         IF( ice_dyn_rdgrft_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'ice_dyn_rdgrft_init: unable to allocate arrays' )
      ENDIF

# if defined _OPENACC
      ! Those declared into `par_ice.F90`:
      !$acc update device( rn_delta_ecc, ln_str_H79, rn_crhg, rn_pstar, ln_str_R75, rn_pe_rdg, ln_str_CST, rn_str )
      ! The rest:
      !$acc update device( nice_str, ln_str_smooth, ln_distf_lin, ln_distf_exp, rn_murdg, rn_csrdg, ln_partf_lin, rn_gstar, ln_partf_exp, rn_astar )
      !$acc update device( ln_ridging, rn_hstar, rn_porordg, rn_fsnwrdg, rn_fpndrdg, ln_rafting, rn_hraft, rn_craft, rn_fsnwrft, rn_fpndrft )
# endif

   END SUBROUTINE ice_dyn_rdgrft_init

   !!======================================================================
END MODULE icedyn_rdgrft
