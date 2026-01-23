MODULE icethd_sal
   !!======================================================================
   !!                       ***  MODULE icethd_sal ***
   !!   sea-ice : computation of salinity variations in the ice
   !!======================================================================
   !! History :   -   !  2003-05  (M. Vancoppenolle) original code 1-D
   !!            3.0  !  2005-12  (M. Vancoppenolle) adapted to the 3-D version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!---------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd_sal      : salinity variations in the ice
   !!   ice_thd_sal_init : initialization
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_oce
   USE par_kind, ONLY : wp
   USE phycst
   USE oss_nnq , ONLY : sss_s
   USE ice
   USE icevar  , ONLY : ice_var_salprof_gpu
   !
   USE in_out_manager , ONLY : numnam_ice_ref, numnam_ice_cfg, numout, numoni, lwp, lwm, ln_timing
   USE lib_mpp        , ONLY : ctl_stop, ctl_warn, ctl_nam                               ! MPP library
   USE iom            , ONLY : iom_miss_val                                              ! I/O manager library
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_sal        ! called by icethd
   PUBLIC   ice_thd_sal_init   ! called by ice_init

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_sal(jl_cat, ll_ice_present)
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_sal  ***
      !!
      !! ** Purpose :   computes new salinities in the ice
      !!
      !! ** Method  :  4 possibilities
      !!               -> nn_icesal = 1 -> Sice = cst    [ice salinity constant in both time & space]
      !!               -> nn_icesal = 2 -> Sice = S(z,t) [Vancoppenolle et al. 2005]
      !!               -> nn_icesal = 3 -> Sice = S(z)   [multiyear ice]
      !!               -> nn_icesal = 4 -> Sice = S(z,t) [Gravity Drainage and Flushing parameterizations]
      !!
      !! ** Case 4 details :
      !!
      !!    For both gravity drainage and flushing, brines are calculated depending on ice temperature (liquidus formulation):
      !!       Sbr = - T / mu                                [linear liquidus]   ( nn_liquidus == 1 )
      !!       Sbr = -18.7 * T - 0.519 * T2 - 0.00535 * T3   [VC2019]            ( nn_liquidus == 2 )
      !!       Sbr = -17.6 * T - 0.389 * T2 - 0.00362 * T3   [Weast]             ( nn_liquidus == 3 )
      !!
      !! ****************
      !! Gravity Drainage
      !! ****************
      !!
      !! we want to solve this equation:
      !! ==============================
      !!    dS/dt = -w dSbr/dz
      !!
      !!    with S   = sea ice salinity
      !!         Sbr = brine salinity
      !!         w   = upwelling Darcy velocity of the return flow (i.e. vertical velocity of the brines, positive downward => >0)
      !!
      !!    discrete form is solved using upward scheme (such as in CICE):
      !!    (S(t+dt)-S(t))/dt = -w(k) * (Sbr(k+1)-Sbr(k))/dz
      !!
      !! 3 schemes are proposed based on the paper from Thomas et al. (2020):
      !! ======================
      !!   0 |  ----------------------------------- surface
      !!     |
      !!     |  ----------------------------------- zc
      !!   z |
      !!     |      Ra > Rac => brine convection
      !!     |
      !!   h |  ------------------------------------ bottom
      !!     v
      !!
      !!     Ra = cp_br * g * beta * (Sbr(z) - Sw) * perm * (h-z) / (cnd_br*visc)   [RWJ2014 formulation]
      !!
      !!        with Ra     : Rayleigh number
      !!             cp_br  : brine heat capacity         (J/m3/K)
      !!             g      : gravity                     (m/s2)
      !!             beta   : saline density coefficient  (g/kg)-1
      !!             cnd_br : brine thermal conductivity  (W/m/K)
      !!             visc   : brine kinematic viscosity   (m2/s)
      !!             Sw     : ocean salinity              (g/kg)
      !!             zc     : critical depth below which convection occurs (m)
      !!             h      : total ice thickness         (m)
      !!             perm   : effective permeability      (m2)
      !!                      = 3.e-8    * (S/Sbr)^3      [SI3]     ( hard coded )
      !!                      = 1.995e-8 * (S/Sbr)^3.1    [Freitag]
      !!                      = 1.0e-8   * (S/Sbr)^3.     [RJW2014]
      !!
      !!    1) === Reese Jones & Worster 2014 (refer to as RJW2014) ===
      !!
      !!        w(z) = - alpha_rjw * cnd_br / cp_br * max(Ra(z)-Rac) * (z-zc)/(h-zc)^2
      !!           with alpha_rjw : intensity parameter
      !!
      !!    2) === Griewank & Notz 2013 (refer to as GN2013) ===
      !!
      !!        w(k) = - alpha_gn/rho * sum( (Ra(kk)-Rac) * dz(kk), [from kk=1 to k] )
      !!           with rho      : brine density       (kg/m3)
      !!                alpha_gn : intensity parameter (kg/m3/s)
      !!
      !!    3) === Cox and Weeks 1988 (refer to as CW1988) ===
      !!
      !!        w(k) = - alpha_cw * 0.0589_wp * MAX( 0._wp, zv_br(z)/rn_vbrc - 1._wp )
      !!           with alpha_cw : intensity parameter
      !!                rn_vbrc  : critical brines volume (for permeability)
      !!
      !! ********
      !! Flushing
      !! ********
      !!
      !! we want to solve this equation:
      !! ==============================
      !!    dS/dt = -w dSbr/dz
      !!
      !!    with Sbr = brine salinity
      !!         w   = upwelling velocity (i.e. vertical velocity of the brines, negative upward => < 0)
      !!
      !!         w   = Fmass / rhob            if v_br > v_brc (= 5%)
      !!             = 0                       otherwise
      !!
      !!           with Fmass = -Flush * rhoi * dh / dt : mass flux (kg/m2/s, >0 since dh<0)
      !!                rhob = rhow * ( 1 + c*Sbr )     : brine density
      !!                v_br = S / Sbr                  : brine volume fraction
      !!                rhoi                            : ice density
      !!                rhow                            : fresh water density (kg/m3)
      !!                c                               : empirical coef (0.8e-3 ‰-1)
      !!           tuning parameters:
      !!                Flush                           : fraction of melt water allowed to percolate thru the ice (30%)
      !!                v_brc                           : critical brine volume above which there is flushing (5%)
      !!
      !!    discrete form is solved using upward scheme (such as in CICE):
      !!    (S(t+dt)-S(t))/dt = -w(k) * (Sbr(k-1)-Sbr(k))/dz
      !!
      !!
      !! ** References
      !! Thomas, M., Vancoppenolle, M., France, J. L., Sturges, W. T., Bakker, D. C. E., Kaiser, J., & von Glasow, R. (2020).
      !!             Tracer measurements in growing sea ice support convective gravity drainage parameterizations.
      !              Journal of Geophysical Research: Oceans, 125, e2019JC015791. https://doi.org/10. 1029/2019JC015791
      !!---------------------------------------------------------------------
      INTEGER,                     INTENT(in) :: jl_cat
      LOGICAL, DIMENSION(jpi,jpj), INTENT(in) :: ll_ice_present
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jk2            ! dummy loop indices
      REAL(wp) ::   z1_time_gd, z1_time_fl
      !
      ! for gravity drainage and flushing
      INTEGER  ::   iter, jc
      REAL(wp) ::   z1_h_i, zhmelt, zc, zcfl, zperm, ztmp, zRae, zdt, zt1, zt2, zt3
      REAL(wp) ::   zv_brmin, zs_brmax, z1_cp_br, z1_cnd_br, z1_visc, z1_c2
      !
      REAL(wp), DIMENSION(nlay_i)     ::   z_mid
      REAL(wp), DIMENSION(nlay_i+1)   ::   z_edge
      REAL(wp), DIMENSION(nlay_i)     ::   zds, zv_br, zRa, zperm_eff, zw_br, zmsk
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zs_br
      !
      ! permeability
      REAL(wp), PARAMETER ::   np_perm_eff = 2     ! 1 = vertical minimum
      !                                              2 = harmonic mean
      ! Rayleigh
      REAL(wp), PARAMETER ::   zcp_br    = 4.e6    ! heat capacity of brine (J/m3)
      REAL(wp), PARAMETER ::   zbeta     = 7.5e-4  ! saline density coefficient (g/kg)-1
      REAL(wp), PARAMETER ::   zcnd_br   = 0.523   ! thermal conductivity of brine W/m/K
      REAL(wp), PARAMETER ::   zvisc     = 1.8e-6  ! Kinematic viscosity of brine
      ! GN scheme constant
      REAL(wp), PARAMETER ::   zrhob_GN  = 1020.   ! Brine density (kg/m3)
      !
      ! for sanity checks
      REAL(wp) ::   zs_min, zcfl_max, zhi, zSi, zdum

      !!clem test
      INTEGER  ::   ndrainage, nflushing
      REAL(wp) ::  zcfl_test
      LOGICAL  ::  lDo
      !!---------------------------------------------------------------------
      IF( ln_timing    )   CALL timing_start('ice_thd_sal')
      !$acc data present( a_i, dh_i_sum_2d, dh_s_sum_2d, h_i, ll_ice_present, sfx_bri, sfx_res, s_i, sss_s, sz_i, t_i, t_su, t_bo )

      SELECT CASE ( nn_icesal )
         !
         !            !---------------------------------------------!
      CASE( 2 )       !  time varying salinity with linear profile  !
         !            !---------------------------------------------!         
         z1_time_gd = MERGE( rDt_ice / rn_time_gd ,  0._wp ,  ln_drainage )
         z1_time_fl = MERGE( rDt_ice / rn_time_fl ,  0._wp ,  ln_flushing )
         !
         !$acc parallel loop collapse(2) create( zds ) private(zds)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1

               zhi = h_i(ji,jj,jl_cat)
               lDo = ( ll_ice_present(ji,jj) .AND. ( zhi > rn_sal_himin ) )
               
               IF( lDo ) THEN
                  zSi = s_i(ji,jj,jl_cat)
                  !
                  ! --- Update ice salinity from brine drainage and flushing --- !
                  IF( t_su(ji,jj,jl_cat) >= rt0 ) THEN             ! flushing (summer time)
                     zds(1) = - MAX( zSi - rn_sal_fl , 0._wp ) * z1_time_fl
                  ELSEIF( t_su(ji,jj,jl_cat) <= t_bo(ji,jj) ) THEN ! gravity drainage
                     zds(1) = - MAX( zSi - rn_sal_gd , 0._wp ) * z1_time_gd
                  ELSE
                     zds(1) = 0._wp
                  ENDIF
                  ! update salinity
                  zSi = zSi + zds(1)
                  zdum = rhoi * a_i(ji,jj,jl_cat) * zhi * r1_Dt_ice
                  ! salt flux
                  sfx_bri(ji,jj) = sfx_bri(ji,jj) - zdum * zds(1)
                  !
                  ! --- salinity must stay inbounds --- !
                  IF( ln_drainage .OR. ln_flushing ) THEN
                     zds(1) =          MAX( 0._wp, rn_simin              - zSi ) ! > 0 if s_i < simin
                     zds(1) = zds(1) + MIN( 0._wp, rn_sinew*sss_s(ji,jj) - zSi ) ! < 0 if s_i > simax
                     ! update salinity
                     s_i(ji,jj,jl_cat) = zSi + zds(1)
                     ! salt flux
                     sfx_res(ji,jj) = sfx_res(ji,jj) - zdum * zds(1)
                  ENDIF
                  !
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
         !
         ! Salinity profile (gives sz_i)
         CALL ice_var_salprof_gpu

         !             !----------------------------------------!
      CASE( 3 )        ! constant salinity with a fixed profile ! (Schwarzacher (1959) multiyear salinity profile (mean = 2.30)
         !             !----------------------------------------!
         CALL ice_var_salprof_gpu

         !
         !             !--------------------------------!
      CASE( 4 )        ! Gravity Drainage and Flushing  !
         !             !--------------------------------!
         !$acc data create( z_mid, z_edge, zds, zv_br, zRa, zperm_eff, zw_br, zmsk, zs_br )
         ! ==============
         ! Initialization
         ! ==============
         !$acc parallel loop
         DO jk = 1, nlay_i
            z_mid (jk) = ( REAL( jk, wp ) - 0.5_wp ) * r1_nlay_i
         END DO
         !$acc end parallel loop
         !$acc parallel loop
         DO jk = 1, nlay_i+1
            z_edge(jk) = ( REAL( jk, wp ) - 1.0_wp ) * r1_nlay_i
         END DO
         !$acc end parallel loop

         ! ================
         ! Gravity Drainage
         ! ================
         IF( ln_drainage ) THEN
            !
            z1_cp_br  = 1._wp / zcp_br
            z1_cnd_br = 1._wp / zcnd_br
            z1_visc   = 1._wp / zvisc
            !LOLO: `z_mid, z_edge` must not be private as they are the same for all `ji,jj` points...
            !$acc parallel loop collapse(2) private( zds, zv_br, zRa, zperm_eff, zw_br, zmsk, zs_br )
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF( ll_ice_present(ji,jj) ) THEN
                     !IF(ji==kdi .AND. jj==kdj) PRINT *, '(0) z_mid =', REAL(z_mid(1),4), REAL(z_mid(2),4)
                     !IF(ji==kdi .AND. jj==kdj) PRINT *, '(0) z_edge =', REAL(z_edge(1),4), REAL(z_edge(2),4), REAL(z_edge(3),4)
                     ! ice thickness ( we do not want to do anything for salt when ice is thinner than the minimum allowed )
                     z1_h_i = 1._wp / MAX( epsi10, h_i(ji,jj,jl_cat) * r1_nlay_i )
                     !
                     ! surface melting (m)
                     zhmelt = dh_s_sum_2d(ji,jj,jl_cat) + dh_i_sum_2d(ji,jj,jl_cat) ! =0 if no melt, <0 otherwise
                     !
                     ! brine salinity
                     zs_brmax = 0._wp
                     !$acc loop seq
                     DO jk = 1, nlay_i
                        zt1 = t_i(ji,jj,jk,jl_cat) - rt0
                        zt2 = zt1 * zt1
                        zt3 = zt2 * zt1
                        zs_br(jk) = -18.7_wp * zt1 - 0.519_wp * zt2 - 0.00535_wp * zt3
                        zs_brmax = MAX( zs_brmax, zs_br(jk) )
                     END DO
                     zs_br(0) = 0._wp               ! brine salinity at the interfaces
                     zs_br(nlay_i+1) = sss_s(ji,jj)
                     !
                     zcfl_max = 0._wp
                     zs_min   = 0._wp
                     IF( h_i(ji,jj,jl_cat) >= rn_sal_himin .AND. zhmelt >= 0._wp .AND. zs_brmax > sss_s(ji,jj) ) THEN
                        !                               ! during melting season, salt flux can turn upward
                        !
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           zmsk(jk) = MERGE( 1._wp, 0._wp,  t_i(ji,jj,jk,jl_cat)-rt0 < -epsi06 )
                        END DO
                        !
                        ! Compute CFL
                        ! ===========
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           zv_br(jk) = zmsk(jk) * sz_i(ji,jj,jk,jl_cat) / zs_br(jk)
                        END DO

                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(1) zmsk =', REAL(zmsk(1),4), REAL(zmsk(2),4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(1) zv_br =', REAL(zv_br(1),4), REAL(zv_br(2),4)

                        ! Effective permeability
                        ! ----------------------
                        IF( np_perm_eff == 1 ) THEN ! Minimum
                           STOP'LOLO: `np_perm_eff == 1` not ready!'
                           DO jk = 1, nlay_i
                              ztmp = MINVAL( zv_br(jk:nlay_i) )
                              zperm_eff(jk) = 3.e-8_wp * ztmp * ztmp * ztmp
                           END DO
                        ELSEIF( np_perm_eff == 2 ) THEN ! Harmonic Mean
                           !$acc loop seq
                           DO jk = 1, nlay_i
                              ztmp = 0._wp
                              !$acc loop seq
                              DO jk2 = jk, nlay_i
                                 zperm = 3.e-8_wp * zv_br(jk2)*zv_br(jk2)*zv_br(jk2)
                                 ztmp = ztmp + 1._wp / MAX( zperm, epsi10 ) !LOLOfix
                              END DO
                              zperm_eff(jk) = REAL( nlay_i-jk+1, wp ) / ztmp
                           END DO
                        END IF
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(2) zperm_eff =', REAL(zperm_eff(1),4), REAL(zperm_eff(2),4)

                        ! Rayleigh number
                        ! ---------------
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           ! Ra = cp_br * g * beta * (Sbr(z) - Sw) * perm * (h-z) / (cnd_br*visc)   [RWJ2014 formulation]
                           zRa(jk) = zcp_br * grav * zbeta * MAX( 0._wp, zs_br(jk) - sss_s(ji,jj)) * zperm_eff(jk) &
                              &                            *  h_i(ji,jj,jl_cat) * ( 1._wp - z_mid(jk) ) * ( z1_cnd_br * z1_visc )
                        END DO
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) h_i(ji,jj,jl_cat) =', REAL(h_i(ji,jj,jl_cat), 4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) A =', REAL(zcp_br * grav * zbeta, 4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) B =', REAL(zs_br(1) - sss_s(ji,jj), 4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) C =', REAL( zperm_eff(1)*  h_i(ji,jj,jl_cat) * ( 1._wp - z_mid(1) ), 4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) D =', REAL(z1_cnd_br * z1_visc, 4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) zRa =', REAL(zRa(1),4), REAL(zRa(2), 4)


                        ! Vertical velocity
                        ! -----------------
                        IF( nn_sal_scheme == 1 ) THEN     ! *** RJW 2014 ***
                           !
                           ! if Ra is everywhere < Rc : no convection => Rae =0
                           ! else                     : convection until zc => Rae /= 0
                           jc = nlay_i+1
                           zc = 0._wp
                           zRae = 0._wp
                           !$acc loop seq
                           DO jk = nlay_i,1,-1
                              IF( zRa(jk) >= rn_Rc_RJW ) THEN
                                 zRae = MAX( zRae , zRa(jk) - rn_Rc_RJW )
                                 zc = z_edge(jk)
                                 jc = jk
                              END IF
                           END DO
                           z1_c2 = 1._wp / ( ( 1._wp - zc )*( 1._wp - zc ) )

                           !$acc loop seq
                           DO jk = 1, nlay_i
                              zw_br(jk) = 0._wp
                           END DO
                           !$acc loop seq
                           DO jk = jc, nlay_i
                              zw_br(jk) = - rn_alpha_RJW * zRae * ( zcnd_br * z1_cp_br ) * ( z_mid(jk) - zc ) * z1_h_i * z1_c2
                           END DO
                           !
                        ELSEIF( nn_sal_scheme == 2 ) THEN    ! *** GN 2013 ***
                           STOP'LOLO: `nn_sal_scheme == 2` not ready!'
                           !zRae = 0._wp
!!$acc loop seq
                           !DO jk = 1, nlay_i
                           !   zRae = zRae + MAX( zRa(jk) - rn_Rc_GN, 0._wp )
                           !   zw_br(jk) = - rn_alpha_GN / zrhob_GN * zRae * h_i(ji,jj,jl_cat) * r1_nlay_i
                           !END DO
                           !
                        ELSEIF( nn_sal_scheme == 3 ) THEN    ! *** CW 1988 ***
                           STOP'LOLO: `nn_sal_scheme == 3` not ready!'
                           !DO jk = 1, nlay_i
                           !   zw_br(jk) = - rn_alpha_CW  * 0.0589_wp * MAX( 0._wp, zv_br(jk)/rn_vbrc - 1._wp )
                           !END DO
                           !
                        END IF

                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(3) zw_br =', REAL(zw_br(1),4), REAL(zw_br(2),4)


                        ! CFL
                        ! ---
                        zcfl_test = 0._wp
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           zcfl = zw_br(jk) * rDt_ice * z1_h_i
                           zcfl_test = MAX( zcfl_test, ABS(zcfl) )
                        END DO
                        !
                        ndrainage = MIN( nn_drainage, CEILING(zcfl_test*2._wp) )
                        zdt   = rDt_ice / REAL( MAX( 1, ndrainage ) )

                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(4) zcfl_testr =', REAL(zcfl_test,4)
                        !IF(ji==kdi .AND. jj==kdj) PRINT *, '(4) ndrainage =', REAL(ndrainage,4)


                        ! Iteration Loop
                        ! ===============
                        ! ndrainage is larger than 1 only if CFL>0.5
                        !$acc loop seq
                        DO iter = 1, ndrainage
                           !
                           IF( iter > 1 ) THEN
                              !$acc loop seq
                              DO jk = 1, nlay_i
                                 zv_br(jk) = zmsk(jk) * sz_i(ji,jj,jk,jl_cat) / zs_br(jk)
                              END DO

                              ! Effective permeability
                              ! ----------------------
                              IF( np_perm_eff == 1 ) THEN ! Minimum
                                 STOP'LOLO: `np_perm_eff == 1` not ready!'
                                 DO jk = 1, nlay_i
                                    ztmp = MINVAL( zv_br(jk:nlay_i) )
                                    zperm_eff(jk) = 3.e-8_wp * ztmp * ztmp * ztmp
                                 END DO
                              ELSEIF( np_perm_eff == 2 ) THEN ! Harmonic Mean
                                 !$acc loop seq
                                 DO jk = 1, nlay_i
                                    ztmp = 0._wp
                                    !$acc loop seq
                                    DO jk2 = jk, nlay_i
                                       zperm = 3.e-8_wp * zv_br(jk2)*zv_br(jk2)*zv_br(jk2)
                                       ztmp = ztmp + 1._wp / MAX( zperm, epsi10 ) !LOLOfix
                                    END DO
                                    zperm_eff(jk) = REAL( nlay_i-jk+1, wp ) / ztmp
                                 END DO
                              END IF

                              ! Rayleigh number
                              ! ---------------
                              !$acc loop seq
                              DO jk = 1, nlay_i
                                 ! Ra = cp_br * g * beta * (Sbr(z) - Sw) * perm * (h-z) / (cnd_br*visc)   [RWJ2014 formulation]
                                 zRa(jk) = zcp_br * grav * zbeta * MAX( 0._wp, zs_br(jk) - sss_s(ji,jj)) * zperm_eff(jk) &
                                    &                            *  h_i(ji,jj,jl_cat) * ( 1._wp - z_mid(jk) ) * ( z1_cnd_br * z1_visc )
                              END DO

                              ! Vertical velocity
                              ! -----------------
                              IF( nn_sal_scheme == 1 ) THEN     ! *** RJW 2014 ***
                                 !
                                 ! if Ra is everywhere < Rc : no convection => Rae =0
                                 ! else                     : convection until zc => Rae /= 0
                                 jc = nlay_i+1
                                 zc = 0._wp
                                 zRae = 0._wp
                                 !$acc loop seq
                                 DO jk = nlay_i,1,-1
                                    IF( zRa(jk) >= rn_Rc_RJW ) THEN
                                       zRae = MAX( zRae , zRa(jk) - rn_Rc_RJW )
                                       zc = z_edge(jk)
                                       jc = jk
                                    END IF
                                 END DO
                                 z1_c2 = 1._wp / ( ( 1._wp - zc )*( 1._wp - zc ) )

                                 !$acc loop seq
                                 DO jk = 1, nlay_i
                                    zw_br(jk) = 0._wp
                                 END DO
                                 !$acc loop seq
                                 DO jk = jc, nlay_i
                                    zw_br(jk) = - rn_alpha_RJW * zRae * ( zcnd_br * z1_cp_br ) * ( z_mid(jk) - zc ) * z1_h_i * z1_c2
                                 END DO

                              ELSEIF( nn_sal_scheme == 2 ) THEN    ! *** GN 2013 ***
                                 STOP'LOLO: `nn_sal_scheme == 2` not ready!'
                                 !zRae = 0._wp
                                 !DO jk = 1, nlay_i
                                 !   zRae = zRae + MAX( zRa(jk) - rn_Rc_GN, 0._wp )
                                 !   zw_br(jk) = - rn_alpha_GN / zrhob_GN * zRae * h_i(ji,jj,jl_cat) * r1_nlay_i
                                 !END DO
                                 !
                              ELSEIF( nn_sal_scheme == 3 ) THEN    ! *** CW 1988 ***
                                 STOP'LOLO: `nn_sal_scheme == 3` not ready!'
                                 !DO jk = 1, nlay_i
                                 !   zw_br(jk) = - rn_alpha_CW  * 0.0589_wp * MAX( 0._wp, zv_br(jk)/rn_vbrc - 1._wp )
                                 !END DO
                                 !
                              END IF
                           ENDIF

                           ! Salinity
                           ! --------
                           ! upstream scheme as in CICE: ds = -w(k)*dt/dz * ( s_br(k+1) - s_br(k) ), w > 0 or < 0
                           !$acc loop seq
                           DO jk = 1, nlay_i
                              !
                              zcfl = zw_br(jk) * zdt * z1_h_i
                              zds(jk) = - zcfl * ( zs_br(jk+1) - zs_br(jk) )
                              !
                              zcfl_max = MAX( zcfl_max, ABS(zcfl) )
                              zs_min   = MIN( zs_min , sz_i(ji,jj,jk,jl_cat) + zds(jk) ) ! record what salinity would be without the trick below
                              !
                              !!clem trick
                              zds(jk) = MAX( zds(jk), -sz_i(ji,jj,jk,jl_cat)+rn_simin )
                              !
                              ! new salinity
                              sz_i(ji,jj,jk,jl_cat) = sz_i(ji,jj,jk,jl_cat) + zds(jk)
                              !
                              ! salt flux
                              sfx_bri(ji,jj) = sfx_bri(ji,jj) - rhoi * a_i(ji,jj,jl_cat) * h_i(ji,jj,jl_cat) * r1_nlay_i * zds(jk) * r1_Dt_ice ! r1_Dt_ice is ok
                              !
                           END DO

                        END DO
                     ENDIF

                     ! sanity check
                     !IF( ln_sal_chk) THEN
                     !   zcfl_drain(ji,jj,jl_cat) = zcfl_max
                     !   IF( zs_min < 0._wp )  THEN
                     !      zsneg_drain(ji,jj,jl_cat) = zs_min
                     !   ELSE
                     !      zsneg_drain(ji,jj,jl_cat) = 0._wp
                     !   ENDIF
                     !ENDIF

                  ENDIF !IF( ll_ice_present(ji,jj) )
               END DO !DO ji=Nis0, Nie0
            END DO !DO jj=Njs0, Nje0
            !$acc end parallel loop

         ENDIF ! ln_drainage


         ! ========
         ! Flushing
         ! ========
         IF( ln_flushing ) THEN
            !
            !$acc parallel loop collapse(2) private( zds, zv_br, zw_br, zmsk, zs_br )
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF( ll_ice_present(ji,jj) ) THEN
                     ! ice thickness ( we do not want to do anything for salt when ice is thinner than the minimum allowed )
                     z1_h_i = 1._wp / MAX( h_i(ji,jj,jl_cat) * r1_nlay_i, epsi10 )
                     !
                     ! surface melting (m)
                     zhmelt = dh_s_sum_2d(ji,jj,jl_cat) + dh_i_sum_2d(ji,jj,jl_cat) ! =0 if no melt, <0 otherwise
                     !
                     zcfl_max = 0._wp
                     zs_min   = 0._wp
                     IF( h_i(ji,jj,jl_cat) >= rn_sal_himin .AND. zhmelt < 0._wp ) THEN     ! Flushing if  surface melting

                        ! brine salinity
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           zt1 = t_i(ji,jj,jk,jl_cat) - rt0
                           zt2 = zt1 * zt1
                           zt3 = zt2 * zt1
                           zs_br(jk) = -18.7_wp * zt1 - 0.519_wp * zt2 - 0.00535_wp * zt3
                        END DO
                        !               ! brine salinity at the interfaces
                        zs_br(0) = 0._wp
                        zs_br(nlay_i+1) = sss_s(ji,jj)

                        !$acc loop seq
                        DO jk = 1, nlay_i
                           IF( (t_i(ji,jj,jk,jl_cat)-rt0) <  - epsi06 ) THEN
                              zmsk(jk) = 1._wp
                           ELSE
                              zmsk(jk) = 0._wp
                           ENDIF
                        END DO


                        ! Compute CFL
                        ! ===========
                        zcfl_test = 0._wp
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           zw_br(jk) = -rn_flushrate * ( dh_i_sum_2d(ji,jj,jl_cat)*rhoi + dh_s_sum_2d(ji,jj,jl_cat)*rhos )  &
                              &                 / ( rhow * ( 1._wp + 0.8e-3_wp * zs_br(jk) ) ) * r1_Dt_ice ! r1_Dt_ice is ok
                           !                        can be replaced by rhow but in theory rhow should be rho_br = (rho0*(1+c*S_br)), with c = 0.8e-3
                           zcfl = zw_br(jk) * rDt_ice * z1_h_i
                           zcfl_test = MAX( zcfl_test, ABS(zcfl) )
                        END DO
                        nflushing = MIN( nn_flushing, CEILING(zcfl_test*2._wp) )
                        zdt   = rDt_ice / REAL( MAX( 1, nflushing ) )

                        !
                        ! iteration Loop
                        ! ==============
                        ! nflushing is larger than 1 only if CFL>0.5
                        !$acc loop seq
                        DO iter = 1, nflushing
                           !
                           zv_brmin = 1.e+20_wp
                           !$acc loop seq
                           DO jk = 1, nlay_i
                              zv_br(jk) = zmsk(jk) * sz_i(ji,jj,jk,jl_cat) / zs_br(jk)
                              zv_brmin  = MIN( zv_brmin, zv_br(jk) )
                           END DO
                           !
                           IF( zv_brmin >= rn_vbrc ) THEN            ! Flushing if brine volume fraction exceeds a certain treshold
                              ! .AND. t_i_1d(ji,1) >= t_i_1d(ji,nlay_i)      !          and surface temperature is warmer than bottom temperature
                              !
                              !$acc loop seq
                              DO jk = 1, nlay_i
                                 ! Vertical velocity
                                 ! -----------------
                                 zw_br(jk) = -rn_flushrate * ( dh_i_sum_2d(ji,jj,jl_cat)*rhoi + dh_s_sum_2d(ji,jj,jl_cat)*rhos )  &
                                    &                 / ( rhow * ( 1._wp + 0.8e-3_wp * zs_br(jk) ) ) * r1_Dt_ice ! r1_Dt_ice is ok
                                 !                        can be replaced by rhow but in theory rhow should be rho_br = (rho0*(1+c*S_br)), with c = 0.8e-3
                                 ! Salinity
                                 ! --------
                                 ! upstream scheme as in CICE: ds = -w*dt/dz * ( s_br(k) - s_br(k-1) ), w > 0
                                 ! zcfl = w*dt/dz
                                 zcfl = zw_br(jk) * zdt * z1_h_i
                                 !
                                 zcfl_max = MAX( zcfl_max, ABS(zcfl) )
                                 !
                                 zds(jk) = - zcfl * ( zs_br(jk) - zs_br(jk-1) )
                                 !
                                 zs_min  = MIN( zs_min , sz_i(ji,jj,jk,jl_cat) + zds(jk) )   ! record what salinity would be without the trick below
                                 !
                                 !!clem trick
                                 zds(jk) = MAX( MIN( 0._wp, zds(jk) ), -sz_i(ji,jj,jk,jl_cat)+rn_simin )
                                 !            ! min to block flushing when temperature profile is not ok
                                 !
                                 ! new salinity
                                 sz_i(ji,jj,jk,jl_cat) = sz_i(ji,jj,jk,jl_cat) + zds(jk)
                                 !
                                 ! salt flux
                                 sfx_bri(ji,jj) = sfx_bri(ji,jj) - rhoi * a_i(ji,jj,jl_cat) * h_i(ji,jj,jl_cat) * r1_nlay_i * zds(jk) * r1_Dt_ice ! r1_Dt_ice is ok

                              END DO
                              !
                           ENDIF

                        END DO
                     ENDIF

                     ! sanity check
                     !IF( ln_sal_chk ) THEN
                     !   zcfl_flush(ji,jj,jl_cat) = zcfl_max
                     !   IF( zs_min < 0._wp ) THEN
                     !      zsneg_flush(ji,jj,jl_cat) = zs_min
                     !   ELSE
                     !      zsneg_flush(ji,jj,jl_cat) = 0._wp
                     !   ENDIF
                     !ENDIF

                  ENDIF !IF( ll_ice_present(ji,jj) )
               END DO !DO ji=Nis0, Nie0
            END DO !DO jj=Njs0, Nje0
            !$acc end parallel loop

         ENDIF

         ! --- salinity must stay inbounds --- !
         IF( ln_drainage .OR. ln_flushing ) THEN
            !$acc parallel loop collapse(2) private( zds )
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF( ll_ice_present(ji,jj) ) THEN
                     !
                     !$acc loop seq
                     DO jk = 1, nlay_i
                        zds(jk) =           MAX( 0._wp, rn_simin                - sz_i(ji,jj,jk,jl_cat) ) ! > 0 if s_i < simin
                        zds(jk) = zds(jk) + MIN( 0._wp, rn_sinew * sss_s(ji,jj) - sz_i(ji,jj,jk,jl_cat) ) ! < 0 if s_i > simax
                        ! update salinity
                        sz_i(ji,jj,jk,jl_cat) = sz_i(ji,jj,jk,jl_cat) + zds(jk)
                        ! salt flux
                        sfx_res(ji,jj) = sfx_res(ji,jj) - rhoi * a_i(ji,jj,jl_cat) * h_i(ji,jj,jl_cat) * r1_nlay_i * zds(jk) * r1_Dt_ice
                     END DO
                     !
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop
         ENDIF

         !$acc end data
         
      END SELECT

      !$acc end data
      IF( ln_timing    )   CALL timing_stop('ice_thd_sal')
      !
   END SUBROUTINE ice_thd_sal

   SUBROUTINE ice_thd_sal_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_thd_sal_init  ***
      !!
      !! ** Purpose :   initialization of ice salinity parameters
      !!
      !! ** Method  :   Read the namthd_sal namelist and check the parameter
      !!                values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_sal
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer
      !!
      NAMELIST/namthd_sal/ nn_icesal, ln_flushing, ln_drainage, rn_sinew, rn_simin, &
         &                 rn_icesal, rn_sal_gd, rn_time_gd, rn_sal_fl, rn_time_fl, &
         &                 rn_sal_himin, nn_liquidus, nn_drainage, nn_flushing, rn_flushrate, rn_vbrc, &
         &                 nn_sal_scheme, rn_alpha_RJW, rn_Rc_RJW, rn_alpha_GN, rn_Rc_GN, rn_alpha_CW, ln_sal_chk
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_sal)
      READ_NML_CFG(numnam_ice,namthd_sal)
      IF(lwm) WRITE ( numoni, namthd_sal )
      !
      IF(lwp) THEN                           ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_sal_init : Ice parameters for salinity '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_sal:'
         WRITE(numout,*) '      switch for salinity                                     nn_icesal     = ', nn_icesal
         WRITE(numout,*) '      activate flushing                                       ln_flushing   = ', ln_flushing
         WRITE(numout,*) '      activate gravity drainage                               ln_drainage   = ', ln_drainage
         WRITE(numout,*) '      New ice salinity (fraction of sss)                      rn_sinew      = ', rn_sinew
         WRITE(numout,*) '      Minimum tolerated ice salinity                          rn_simin      = ', rn_simin
         ! -- nn_icesal=1 -- !
         WRITE(numout,*) '      bulk salinity value if nn_icesal = 1                    rn_icesal     = ', rn_icesal
         ! -- nn_icesal=2 -- !
         WRITE(numout,*) '      restoring salinity for gravity drainage                 rn_sal_gd     = ', rn_sal_gd
         WRITE(numout,*) '      restoring time for for gravity drainage                 rn_time_gd    = ', rn_time_gd
         WRITE(numout,*) '      restoring salinity for flushing                         rn_sal_fl     = ', rn_sal_fl
         WRITE(numout,*) '      restoring time for flushing                             rn_time_fl    = ', rn_time_fl
         ! -- nn_icesal=4 -- !
         WRITE(numout,*) '      min ice thickness for drainage and flushing             rn_sal_himin  = ', rn_sal_himin
         WRITE(numout,*) '      liquidous formulation (1=linear, 2=VC2019, 3=Weast71)   nn_liquidus   = ', nn_liquidus
         WRITE(numout,*) '      number of subcycles for gravity drainage                nn_drainage   = ', nn_drainage
         WRITE(numout,*) '      number of subcycles for flushing                        nn_flushing   = ', nn_flushing
         WRITE(numout,*) '      fraction of melt water used for flushing                rn_flushrate  = ', rn_flushrate
         WRITE(numout,*) '      critical brines volume above which flushing can occur   rn_vbrc       = ', rn_vbrc
         WRITE(numout,*) '      convection scheme (1=RJW2014, 2=GN2013, 3=CW88)         nn_sal_scheme = ', nn_sal_scheme
         WRITE(numout,*) '      brine flow for RJW2014 scheme                           rn_alpha_RJW  = ', rn_alpha_RJW
         WRITE(numout,*) '      critical Rayleigh number for RJW2014 scheme             rn_Rc_RJW     = ', rn_Rc_RJW
         WRITE(numout,*) '      brine flow for GN2013 scheme (kg/m3/s)                  rn_alpha_GN   = ', rn_alpha_GN
         WRITE(numout,*) '      critical Rayleigh number for GN2013 scheme              rn_Rc_GN      = ', rn_Rc_GN
         WRITE(numout,*) '      brine flow for CW1988 scheme                            rn_alpha_CW   = ', rn_alpha_CW
         WRITE(numout,*) '      sanity checks (output diags)                            ln_sal_chk    = ', ln_sal_chk
      ENDIF
      !
      IF( nn_icesal /= 4 )   ln_sal_chk=.FALSE. ! option only valid for nn_icesal = 4
      !
      !$acc update device( nn_icesal, rn_icesal, rn_sinew, rn_simin, ln_sal_chk, nn_liquidus, rn_sal_gd, rn_time_gd, rn_sal_fl, rn_time_fl, nn_sal_scheme, ln_flushing )
      !$acc update device( ln_drainage, nn_drainage, nn_flushing, rn_flushrate, rn_alpha_CW, rn_alpha_RJW, rn_alpha_GN, rn_Rc_RJW, rn_Rc_GN, rn_sal_himin, rn_vbrc )
      !
   END SUBROUTINE ice_thd_sal_init

   !!======================================================================
END MODULE icethd_sal
