MODULE icethd_zdf_BL99
   !!======================================================================
   !!                       ***  MODULE icethd_zdf_BL99 ***
   !!   sea-ice: vertical heat diffusion in sea ice (computation of temperatures)
   !!======================================================================
   !! History :       !  2003-02  (M. Vancoppenolle) original 1D code
   !!                 !  2005-06  (M. Vancoppenolle) 3d version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  ice_thd_zdf_BL99 : vertical diffusion computation
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_ice        ! SI3 parameters
   USE par_kind, ONLY : wp
   USE phycst
   USE ice
   USE sbc_ice , ONLY : qns_ice, dqns_ice, qcn_ice, qtr_ice_top, qsr_ice

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_zdf_BL99   ! called by icethd_zdf

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_zdf_BL99( jl_cat, k_cnd, ll_ice_present )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_zdf_BL99  ***
      !!
      !! ** Purpose : computes the time evolution of snow and sea-ice temperature
      !!              profiles, using the original Bitz and Lipscomb (1999) algorithm
      !!
      !! ** Method  : solves the heat equation diffusion with a Neumann boundary
      !!              condition at the surface and a Dirichlet one at the bottom.
      !!              Solar radiation is partially absorbed into the ice.
      !!              The specific heat and thermal conductivities depend on ice
      !!              salinity and temperature to take into account brine pocket
      !!              melting. The numerical scheme is an iterative Crank-Nicolson
      !!              on a non-uniform multilayer grid in the ice and snow system.
      !!
      !!           The successive steps of this routine are
      !!           1.  initialization of ice-snow layers thicknesses
      !!           2.  Internal absorbed and transmitted radiation
      !!           Then iterative procedure begins
      !!           3.  Thermal conductivity
      !!           4.  Kappa factors
      !!           5.  specific heat in the ice
      !!           6.  eta factors
      !!           7.  surface flux computation
      !!           8.  tridiagonal system terms
      !!           9.  solving the tridiagonal system with Gauss elimination
      !!           Iterative procedure ends according to a criterion on evolution
      !!           of temperature
      !!           10. Fluxes at the interfaces
      !!
      !! ** Inputs / Ouputs : (global commons)
      !!           surface temperature              : t_su
      !!           ice/snow temperatures            : t_i, t_s
      !!           ice salinities                   : sz_i
      !!           number of layers in the ice/snow : nlay_i, nlay_s
      !!           total ice/snow thickness         : h_i, h_s
      !!-------------------------------------------------------------------
      INTEGER,                     INTENT(in) ::   jl_cat        ! ice-category we are working with...
      INTEGER,                     INTENT(in) ::   k_cnd     ! conduction flux (off, on, emulated)
      LOGICAL, DIMENSION(jpi,jpj), INTENT(in) ::   ll_ice_present
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj, jk                 ! spatial loop index
      INTEGER ::   jm                         ! current reference number of equation
      INTEGER ::   iconv                      ! number of iterations in iterative procedure
      INTEGER ::   iconv_max = 50             ! max number of iterations in iterative procedure
      INTEGER ::   k_T_converged              ! `1` when T has converged (per grid point)
      !
      REAL(wp) ::   zg1s      =  2._wp        ! for the tridiagonal system
      REAL(wp) ::   zg1       =  2._wp        !
      REAL(wp) ::   zgamma    =  18009._wp    ! for specific heat
      REAL(wp) ::   zbeta     =  0.117_wp     ! for thermal conductivity (could be 0.13)
      REAL(wp) ::   zkimin    =  0.10_wp      ! minimum ice thermal conductivity
      REAL(wp) ::   ztsu_err  =  1.e-5_wp     ! range around which t_su is considered at 0C
      REAL(wp) ::   zdti_bnd  =  1.e-4_wp     ! maximal authorized error on temperature
      REAL(wp) ::   zhs_ssl   =  0.03_wp      ! surface scattering layer in the snow
      REAL(wp) ::   zhi_ssl   =  0.10_wp      ! surface scattering layer in the ice
      REAL(wp) ::   zh_min    =  1.e-3_wp     ! minimum ice/snow thickness for conduction
      REAL(wp) ::   ztmelts                   ! ice melting temperature
      REAL(wp) ::   zdti_max                  ! current maximal error on temperature
      REAL(wp) ::   zcpi                      ! Ice specific heat
      REAL(wp) ::   zhfx_err, zdq             ! diag errors on heat
      REAL(wp) ::   zfac
      !
      REAL(wp) ::   za_s_fra    ! ice fraction covered by snow
      REAL(wp) ::   zraext_s     ! extinction coefficient of radiation in the snow
      REAL(wp) ::   zghe        ! G(he), th. conduct enhancement factor, mono-cat
      REAL(wp) ::   z1_h_i, z1_h_s
      REAL(wp) ::   zisnow_comb  ! snow presence for met-office
      !
      REAL(wp) ::   ztsub        ! surface temperature at previous iteration
      REAL(wp) ::   zh_i         ! ice layer thickness
      REAL(wp) ::   zh_s         ! snow layer thickness
      REAL(wp) ::   zqns_ice_b   ! solar radiation absorbed at the surface
      REAL(wp) ::   zdqns_ice_b  ! derivative of the surface flux function
      !
      REAL(wp) ::   zkappa_comb ! Combined snow and ice surface conductivity
      REAL(wp) ::   zq_ini      ! diag errors on heat
      REAL(wp) ::   zisnow      ! snow presence (1) or not (0)
      !
      INTEGER  ::   jm_min    ! reference number of top equation
      INTEGER  ::   jm_max    ! reference number of bottom equation
      REAL(wp) ::   zfnet     ! surface flux function
      !
      REAL(wp), DIMENSION(0:nlay_s) ::   zradtr_s    ! Radiation transmited through the snow
      REAL(wp), DIMENSION(0:nlay_s) ::   zradab_s    ! Radiation absorbed in the snow
      REAL(wp), DIMENSION(0:nlay_i) ::   zradtr_i    ! Radiation transmitted through the ice
      REAL(wp), DIMENSION(0:nlay_i) ::   zradab_i    ! Radiation absorbed in the ice
      REAL(wp), DIMENSION(0:nlay_i) ::   ztcond_i    ! Ice thermal conductivity
      REAL(wp), DIMENSION(0:nlay_i) ::   ztcond_i_cp ! copy
      REAL(wp), DIMENSION(nlay_i)   ::   ztiold      ! Old temperature in the ice
      REAL(wp), DIMENSION(nlay_s)   ::   ztsold      ! Old temperature in the snow
      REAL(wp), DIMENSION(0:nlay_i) ::   zkappa_i    ! Kappa factor in the ice
      REAL(wp), DIMENSION(0:nlay_i) ::   zeta_i      ! Eta factor in the ice
      REAL(wp), DIMENSION(0:nlay_s) ::   zkappa_s    ! Kappa factor in the snow
      REAL(wp), DIMENSION(0:nlay_s) ::   zeta_s      ! Eta factor in the snow
      REAL(wp), DIMENSION(nlay_i)   ::   ztib        ! Temporary temperature in the ice to check the convergence
      REAL(wp), DIMENSION(nlay_s)   ::   ztsb        ! Temporary temperature in the snow to check the convergence
      !
      REAL(wp), DIMENSION(nlay_i+nlay_s+1)   ::   zindterm    ! 'Ind'ependent term
      REAL(wp), DIMENSION(nlay_i+nlay_s+1)   ::   zindtbis    ! Temporary 'ind'ependent term
      REAL(wp), DIMENSION(nlay_i+nlay_s+1)   ::   zdiagbis    ! Temporary 'dia'gonal term
      REAL(wp), DIMENSION(nlay_i+nlay_s+1,3) ::   ztrid       ! Tridiagonal system terms
      !
      ! zradtr_s,zradab_s,zradtr_i,zradab_i,ztcond_i,ztcond_i_cp,ztiold,ztsold,zkappa_i,zeta_i,zkappa_s,zeta_s,ztib,ztsb,zindterm,zindtbis,zdiagbis,ztrid
      ! Mono-category
      REAL(wp) ::   zepsilon, zthres  ! determines thres. above which computation of G(h) is done
      REAL(wp) ::   zhe               ! dummy factor
      REAL(wp) ::   zcnd_i            ! mean sea ice thermal conductivity
      REAL(wp) ::   z1_hi_ssl, zt_su, zA, zdum, zsum_i, zsum_s
      !!------------------------------------------------------------------
      !$acc data present( cnd_ice,dqns_ice,hfx_dif,hfx_err_dif,h_i,h_s,ll_ice_present,qcn_ice,qcn_ice_bot,qcn_ice_top,qns_ice,qsr_ice,qtr_ice_bot,qtr_ice_top,sz_i,t1_ice )
      !$acc data create( zradtr_s,zradab_s,zradtr_i,zradab_i,ztcond_i,ztcond_i_cp,ztiold,ztsold,zkappa_i,zeta_i,zkappa_s,zeta_s,ztib,ztsb,zindterm,zindtbis,zdiagbis,ztrid )

      IF( ln_virtual_itd ) THEN
         zepsilon = 0.1_wp
         zthres   = zepsilon * 0.5_wp * EXP(1._wp)
      ENDIF

      z1_hi_ssl = 1._wp / zhi_ssl      
      
      !$acc parallel loop collapse(2) private(zradtr_s,zradab_s,zradtr_i,zradab_i,ztcond_i,ztcond_i_cp,ztiold,ztsold,zkappa_i,zeta_i,zkappa_s,zeta_s,ztib,ztsb,zindterm,zindtbis,zdiagbis,ztrid)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            IF( ll_ice_present(ji,jj) ) THEN

               zt_su = t_su(ji,jj,jl_cat)
               
               zsum_i = 0._wp ; zsum_s = 0._wp
               !$acc loop seq
               DO jk = 1, nlay_s
                  zsum_s = zsum_s + e_s(ji,jj,jk,jl_cat)
               END DO
               !$acc loop seq
               DO jk = 1, nlay_i
                  zsum_i = zsum_i + e_i(ji,jj,jk,jl_cat)
               END DO
               zq_ini = ( zsum_i * h_i(ji,jj,jl_cat) * r1_nlay_i + zsum_s * h_s(ji,jj,jl_cat) * r1_nlay_s )

               !------------------
               ! 1) Initialization
               !------------------
               !
               ! thicknesses
               zh_i = h_i(ji,jj,jl_cat)
               zh_i = MERGE( MAX( zh_min , zh_i )*r1_nlay_i, 0._wp, zh_i>0._wp )
               zh_s = h_s(ji,jj,jl_cat)
               zisnow = MERGE( 1._wp,                 0._wp, zh_s>0._wp )
               zh_s = MERGE( MAX( zh_min , zh_s )*r1_nlay_s, 0._wp, zh_s>0._wp )

               ! clem: we should apply correction on snow thickness to take into account snow fraction
               !       it must be a distribution, so it is a bit complicated
               !
               ! Store initial temperatures and non solar heat fluxes
               !IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
               ztsub              =      zt_su        ! surface temperature at iteration n-1
               zt_su = MIN( ztsub, rt0 - ztsu_err )   ! required to leave the choice between melting or not
               zdqns_ice_b        =  dqns_ice(ji,jj,jl_cat)        ! derivative of incoming nonsolar flux
               zqns_ice_b         =   qns_ice(ji,jj,jl_cat)        ! store previous qns_ice_1d value
               !ENDIF
               !
               !$acc loop seq
               DO jk = 1, nlay_s
                  ztsold(jk) = t_s(ji,jj,jk,jl_cat)   ! Old snow temperature
               END DO
               !$acc loop seq
               DO jk = 1, nlay_i
                  ztiold(jk) = t_i(ji,jj,jk,jl_cat)   ! Old ice temperature
               END DO


               !-------------
               ! 2) Radiation
               !-------------
               ! --- Transmission/absorption of solar radiation in the ice --- !
               ! extinction radiation in the snow
               zraext_s = MERGE( rn_kappa_sdry ,  rn_kappa_smlt ,  zt_su < rt0 )
               zraext_s = MERGE( rn_kappa_s    ,  zraext_s      ,           nn_qtrice == 0  )
               !
               zradtr_s(0) = qtr_ice_top(ji,jj,jl_cat)
               !$acc loop seq
               DO jk = 1, nlay_s
                  zradtr_s(jk) = zradtr_s(0) * EXP( - zraext_s * MAX( 0._wp, zh_s * REAL(jk) - zhs_ssl ) ) ! radiation transmitted below the layer-th snow layer
                  zradab_s(jk) = zradtr_s(jk-1) - zradtr_s(jk) ! radiation absorbed by the layer-th snow layer
               END DO

               ! Manual inlining of `ice_var_snwfra`:
               zdum = h_s(ji,jj,jl_cat)
               za_s_fra = zdum / ( zdum + 0.02_wp ) ! snow cover depends on hsnow (CICE style)
               !
               zradtr_i(0) = zradtr_s(nlay_s) * za_s_fra + qtr_ice_top(ji,jj,jl_cat) * ( 1._wp - za_s_fra )
               !$acc loop seq
               DO jk = 1, nlay_i
                  zdum = zh_i * REAL(jk)
                  !                             ! radiation transmitted below the layer-th ice layer
                  zradtr_i(jk) =           za_s_fra   * zradtr_s(nlay_s)                       &   ! part covered by snow
                     &                                       * EXP( - rn_kappa_i * MAX( 0._wp, zdum - zh_min  ) ) &
                     &            + ( 1._wp - za_s_fra ) * qtr_ice_top(ji,jj,jl_cat)                        &   ! part snow free
                     &                                       * EXP( - rn_kappa_i * MAX( 0._wp, zdum - zhi_ssl ) )
                  !                             ! radiation absorbed by the layer-th ice layer
                  zradab_i(jk) = zradtr_i(jk-1) - zradtr_i(jk)
               END DO

               qtr_ice_bot(ji,jj,jl_cat) = zradtr_i(nlay_i)   ! record radiation transmitted below the ice

            ENDIF ! IF( ll_ice_present(ji,jj) )



            !************************************************************
            !                      Iteration block
            !************************************************************

            k_T_converged = MERGE( 0, 1,  ll_ice_present(ji,jj) )  ! => 1, aka "converged" where no ice

            ! Convergence calculated until all sub-domain grid points have converged
            ! Calculations keep going for all grid points until sub-domain convergence (vectorisation optimisation)
            ! but values are not taken into account (results independant of MPI partitioning)
            !
            iconv = 0          ! number of iterations
            !                                                                                  !============================!
            DO WHILE ( (k_T_converged < 1).AND.( iconv < iconv_max ) )   ! Iterative procedure begins !
               !                                                                               !============================!
               iconv = iconv + 1

               IF( ll_ice_present(ji,jj) ) THEN

                  ! thicknesses
                  zh_i = h_i(ji,jj,jl_cat)
                  zh_i = MERGE( MAX( zh_min , zh_i )*r1_nlay_i, 0._wp, zh_i>0._wp )
                  zh_s = h_s(ji,jj,jl_cat)
                  zh_s = MERGE( MAX( zh_min , zh_s )*r1_nlay_s, 0._wp, zh_s>0._wp )

                  !--------------------------------
                  ! 3) Sea ice thermal conductivity
                  !--------------------------------
                  IF( ln_cndi_U64 ) THEN         !-- Untersteiner (1964) formula: k = k0 + beta.S/T
                     ztcond_i_cp(0) = rcnd_i + zbeta * sz_i(ji,jj,1,jl_cat)      / MIN( -epsi10, t_i(ji,jj,1,jl_cat) - rt0 )
                     !$acc loop seq
                     DO jk = 1, nlay_i-1
                        ztcond_i_cp(jk) = rcnd_i + zbeta * 0.5_wp * ( sz_i(ji,jj,jk,jl_cat) + sz_i(ji,jj,jk+1,jl_cat) ) /  &
                           &              MIN( -epsi10, 0.5_wp * (  t_i(ji,jj,jk,jl_cat) +  t_i(ji,jj,jk+1,jl_cat) ) - rt0 )
                     END DO
                     ztcond_i_cp(nlay_i) = rcnd_i + zbeta * sz_i(ji,jj,nlay_i,jl_cat) / MIN( -epsi10, t_bo(ji,jj)  - rt0 )
                     !
                  ELSEIF( ln_cndi_P07 ) THEN     !-- Pringle et al formula: k = k0 + beta1.S/T - beta2.T
                     ztcond_i_cp(0) = rcnd_i + 0.09_wp  *  sz_i(ji,jj,1,jl_cat)      / MIN( -epsi10, t_i(ji,jj,1,jl_cat) - rt0 )  &
                        &                            - 0.011_wp * ( t_i(ji,jj,1,jl_cat) - rt0 )
                     !$acc loop seq
                     DO jk = 1, nlay_i-1
                        zdum = 0.5_wp * (  t_i(ji,jj,jk,jl_cat) +  t_i(ji,jj,jk+1,jl_cat) ) - rt0
                        ztcond_i_cp(jk) = rcnd_i + 0.09_wp*0.5_wp * ( sz_i(ji,jj,jk,jl_cat) + sz_i(ji,jj,jk+1,jl_cat) ) &
                           &               / MIN(-epsi10,zdum) - 0.011_wp*zdum
                     END DO
                     ztcond_i_cp(nlay_i) = rcnd_i + 0.09_wp  *  sz_i(ji,jj,nlay_i,jl_cat) / MIN( -epsi10, t_bo(ji,jj)  - rt0 )  &
                        &                            - 0.011_wp * ( t_bo(ji,jj) - rt0 )
                     !
                  ENDIF
                  !
                  ! Variable used after iterations
                  ! Value must be frozen after convergence for MPP independance reason
                  IF( k_T_converged==0 ) THEN
                     !$acc loop seq
                     DO jk = 0, nlay_i
                        ztcond_i(jk) = MAX( zkimin, ztcond_i_cp(jk) )
                     END DO
                  ENDIF

                  !--- G(he) : enhancement of thermal conductivity in mono-category case
                  ! Computation of effective thermal conductivity G(h)
                  ! Used in mono-category case only to simulate an ITD implicitly
                  ! Fichefet and Morales Maqueda, JGR 1997
                  zghe = 1._wp
                  !
                  !LOLO CHECK:
                  IF( ln_virtual_itd ) THEN
                     !
                     zdum = 0._wp
                     !$acc loop seq
                     DO jk = 0, nlay_i
                        zdum = zdum + ztcond_i(jk)  ! => SUM( ztcond_i(:) )
                     END DO
                     zcnd_i = zdum / REAL( nlay_i+1, wp )                            ! Mean sea ice thermal conductivity
                     zhe = ( rcnd_s * h_i(ji,jj,jl_cat) + zcnd_i * h_s(ji,jj,jl_cat) ) / ( rcnd_s + zcnd_i )        ! Effective thickness he (zhe)
                     !IF( zhe >=  zthres )  zghe = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( 2._wp * zhe / zepsilon ) ) )   ! G(he)
                     zghe = MERGE( MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( 2._wp * zhe / zepsilon ) ) ) ,  zghe ,  zhe >=  zthres )
                     !
                  ENDIF
                  !LOLO CHECK.

                  IF( k_T_converged==0 ) THEN

                     !-----------------
                     ! 4) kappa factors
                     !-----------------
                     z1_h_s      = MERGE( 1._wp / MAX(zh_s, epsi10),  0._wp,  h_s(ji,jj,jl_cat)>0._wp  )
                     z1_h_i      = MERGE( 1._wp / MAX(zh_i, epsi10),  0._wp,  h_i(ji,jj,jl_cat)>0._wp  )
                     zisnow_comb = MERGE( h_s(ji,jj,jl_cat) / zh_min, 1._wp,  h_s(ji,jj,jl_cat)<zh_min )

                     !--- Snow
                     ! Variable used after iterations
                     ! Value must be frozen after convergence for MPP independance reason
                     !$acc loop seq
                     DO jk = 0, nlay_s-1
                        zkappa_s(jk) = zghe * rcnd_s * z1_h_s
                     END DO
                     zfac = 0.5_wp * (  ztcond_i(0) * zh_s + rcnd_s * zh_i )
                     IF( zfac<=0._wp) THEN
                        PRINT *, 'LOLO: `zfac<=0.` !!! `icethd_zdf_bl99.F90` (did not think it could be negative)'
                        STOP
                     ENDIF
                     !zkappa_s(nlay_s) = zisnow * zghe * rcnd_s * ztcond_i(0) / zfac   ! Snow-ice interface
                     zkappa_s(nlay_s) = zisnow * zghe * rcnd_s * ztcond_i(0) / MAX( zfac, epsi10 ) !LOLO
                     !
                     !--- Ice
                     ! Variable used after iterations
                     ! Value must be frozen after convergence for MPP independance reason
                     !$acc loop seq
                     DO jk = 0, nlay_i
                        zkappa_i(jk) = zghe * ztcond_i(jk) * z1_h_i
                     END DO
                     ! Calculate combined surface snow and ice conductivity to pass through the coupler (met-office)
                     zkappa_comb = zisnow_comb * zkappa_s(0) + ( 1._wp - zisnow_comb ) * zkappa_i(0)
                     ! If there is snow then use the same snow-ice interface conductivity for the top layer of ice
                     zkappa_i(0) = MERGE( zkappa_s(nlay_s) ,  zkappa_i(0) ,  h_s(ji,jj,jl_cat) > 0._wp )  ! Snow-ice interface


                     IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
#                       include "icethd_zdf_bl99_cnd_OFF.h90"
                     ELSEIF( k_cnd == np_cnd_ON ) THEN
#                       include "icethd_zdf_bl99_cnd_ON.h90"
                     ENDIF !IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU )
                     IF( zdti_max < zdti_bnd )   k_T_converged = 1


                  ENDIF !IF( k_T_converged==0 )

               ENDIF !IF( ll_ice_present(ji,jj) )

            END DO !DO WHILE ( ANY(k_T_converged(Nis0:Nie0,Njs0:Nje0)<1) .AND. (iconv < iconv_max) )

            !************************************************************
            !                End of iteration block
            !************************************************************

            IF( ll_ice_present(ji,jj) ) THEN

               zA = a_i(ji,jj,jl_cat)

               !-----------------------------
               ! 10) Fluxes at the interfaces
               !-----------------------------
               zdum = 1._wp - zisnow
               !
               IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
                  ! --- ice conduction fluxes (positive downward)
                  qcn_ice_bot(ji,jj,jl_cat) = - zkappa_i(nlay_i) * zg1 * ( t_bo(ji,jj ) - t_i(ji,jj,nlay_i,jl_cat) )                ! bottom
                  qcn_ice_top(ji,jj,jl_cat) = - zisnow   * zkappa_s(0) * zg1s * ( t_s(ji,jj,1,jl_cat) - zt_su ) & ! surface
                     &                 - zdum * zkappa_i(0) * zg1  * ( t_i(ji,jj,1,jl_cat) - zt_su )
                  !
                  ! --- Diagnose the heat loss due to changing non-solar / conduction flux
                  hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) - ( qns_ice(ji,jj,jl_cat) - zqns_ice_b ) * zA
                  !
               ELSEIF( k_cnd == np_cnd_ON ) THEN
                  ! --- ice conduction fluxes (positive downward)
                  qcn_ice_bot(ji,jj,jl_cat) = - zkappa_i(nlay_i) * zg1 * ( t_bo(ji,jj ) - t_i(ji,jj,nlay_i,jl_cat) ) ! bottom
                  qcn_ice_top(ji,jj,jl_cat) = qcn_ice(ji,jj,jl_cat)                                                        ! surface
                  !
                  ! --- surface ice temperature
                  IF( ln_cndemulate ) THEN
                     zt_su = ( qcn_ice_top(ji,jj,jl_cat) + zisnow * zkappa_s(0) * zg1s * t_s(ji,jj,1,jl_cat) + &
                        &                                zdum * zkappa_i(0) * zg1  * t_i(ji,jj,1,jl_cat) ) &
                        &          / MAX( epsi10, zisnow * zkappa_s(0) * zg1s + zdum * zkappa_i(0) * zg1 )
                     zt_su = MAX( MIN( zt_su, rt0 ), rt0 - 100._wp )  ! cap t_su
                  ENDIF
                  !
               ENDIF
               !
               ! --- Diagnose the heat loss due to non-fully converged temperature solution (should not be larger than 10-4 W-m2)
               !
               IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_ON ) THEN
                  !
                  !CALL ice_var_enthalpy(jl_cat, ll_ice_present)  ! ==> manual inlining:
                  !$acc loop seq
                  DO jk = 1, nlay_i             ! Sea ice energy of melting
                     ztmelts       = - rTmlt  * sz_i(ji,jj,jk,jl_cat)
                     t_i(ji,jj,jk,jl_cat) = MIN( t_i(ji,jj,jk,jl_cat), ztmelts + rt0 ) ! Force t_i_1d to be lower than melting point => likely conservation issue
                     !   (sometimes zdf scheme produces abnormally high temperatures)
                     zdum = t_i(ji,jj,jk,jl_cat) - rt0
                     IF( zdum==0._wp) THEN
                        PRINT *, 'LOLO: `zdum==0` !!! `icethd_zdf_bl99.F90`'
                        STOP
                     ENDIF
                     e_i(ji,jj,jk,jl_cat) = rhoi*( rcpi*( ztmelts - ( zdum ) ) + rLfus*( 1._wp - ztmelts /  zdum  ) - rcp*ztmelts )
                  END DO
                  !$acc loop seq
                  DO jk = 1, nlay_s             ! Snow energy of melting
                     e_s(ji,jj,jk,jl_cat) = rhos*( rcpi*( rt0 - t_s(ji,jj,jk,jl_cat) ) + rLfus )
                  END DO

                  ! zhfx_err = correction on the diagnosed heat flux due to non-convergence of the algorithm used to solve heat equation
                  zsum_i = 0._wp ; zsum_s = 0._wp
                  !$acc loop seq
                  DO jk = 1, nlay_s
                     zsum_s = zsum_s + e_s(ji,jj,jk,jl_cat)
                  END DO
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     zsum_i = zsum_i + e_i(ji,jj,jk,jl_cat)
                  END DO
                  zdq = - zq_ini + ( zsum_i * h_i(ji,jj,jl_cat) * r1_nlay_i + zsum_s * h_s(ji,jj,jl_cat) * r1_nlay_s )

                  zdum   = zradtr_i(nlay_i) - qcn_ice_bot(ji,jj,jl_cat) + zdq * r1_Dt_ice
                  zsum_i = qcn_ice_top(ji,jj,jl_cat) + qtr_ice_top(ji,jj,jl_cat) - zdum
                  IF( k_cnd == np_cnd_OFF ) THEN
                     zhfx_err = MERGE( (     qns_ice(ji,jj,jl_cat) +     qsr_ice(ji,jj,jl_cat) - zdum ) ,  & ! case T_su < 0degC
                        &                                      zsum_i                                   ,  & ! case T_su = 0degC
                        &                         zt_su < rt0       ) * zA
                  ELSEIF( k_cnd == np_cnd_ON ) THEN
                     zhfx_err =  zsum_i * zA
                  ENDIF
                  
                  ! total heat sink to be sent to the ocean
                  hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) + zhfx_err
                  !
                  ! hfx_dif = Heat flux diagnostic of sensible heat used to warm/cool ice in W.m-2
                  hfx_dif(ji,jj) = hfx_dif(ji,jj) - zdq * r1_Dt_ice * zA
                  !
                  !
               ENDIF


               !--------------------------------------------------------------------
               ! 11) reset inner snow and ice temperatures, update conduction fluxes
               !--------------------------------------------------------------------
               ! effective conductivity and 1st layer temperature (needed by Met Office)
               ! this is a conductivity at mid-layer, hence the factor 2
               cnd_ice(ji,jj,jl_cat) = 2._wp * MERGE( zkappa_comb ,  ztcond_i(0) * z1_hi_ssl ,  h_i(ji,jj,jl_cat) >= zhi_ssl ) ! 2nd: cnd_ice is capped by: cond_i/zhi_ssl
               t1_ice(ji,jj,jl_cat) = zisnow * t_s(ji,jj,1,jl_cat) + ( 1._wp - zisnow ) * t_i(ji,jj,1,jl_cat)
               !
               IF( k_cnd == np_cnd_EMU ) THEN
                  ! Restore temperatures to their initial values
                  !$acc loop seq
                  DO jk = 1, nlay_s
                     t_s(ji,jj,jk,jl_cat) = ztsold(jk)
                  END DO
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     t_i(ji,jj,jk,jl_cat) = ztiold(jk)
                  END DO
                  qcn_ice(ji,jj,jl_cat) = qcn_ice_top(ji,jj,jl_cat)
               ENDIF
               !
               ! --- SIMIP diagnostics (Snow-ice interfacial temperature)
               IF( h_s(ji,jj,jl_cat) >= zhs_ssl ) THEN
                  zdum   = h_i(ji,jj,jl_cat) * r1_nlay_i
                  zsum_s = ztcond_i(1) * h_s(ji,jj,jl_cat) * r1_nlay_s
                  IF( (rcnd_s*zdum + zsum_s)==0._wp) THEN
                     PRINT *, 'LOLO: `(rcnd_s*zdum + zsum_s)==0` !!! `icethd_zdf_bl99.F90`'
                     STOP
                  ENDIF
                  t_si(ji,jj,jl_cat) = ( rcnd_s*zdum*t_s(ji,jj,nlay_s,jl_cat) + zsum_s*t_i(ji,jj,1,jl_cat) ) / ( rcnd_s*zdum + zsum_s )
               ELSE
                  t_si(ji,jj,jl_cat) = zt_su
               ENDIF
               
               t_su(ji,jj,jl_cat) = zt_su
               
            ENDIF !IF( ll_ice_present(ji,jj) )

         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !$acc end parallel loop

      !$acc end data
      !$acc end data
   END SUBROUTINE ice_thd_zdf_BL99

   !!======================================================================
END MODULE icethd_zdf_BL99
