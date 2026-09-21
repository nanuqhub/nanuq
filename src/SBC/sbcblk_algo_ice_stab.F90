MODULE sbcblk_algo_ice_stab
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ice_stab  ***
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!   What is the "STAB" algorithm ?
   !!    => Given a constant value for the neutral coefficients C_D_N, C_E_N and C_H_N over sea-ice
   !!       it will commpute the C_D, C_E and C_H consistent with the near-surface atmospheric stability
   !!    ==> which is already a better approach than using a constant value for C_D, C_E and C_H as in
   !!        the default (and simplest) option in NEMO for instance.
   !!    ==> the only room for improvement is the pick of the stability functions: we use those of Andreas 2005 for now...
   !!
   !! What it computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H over sea-ice
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ub (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_ice_stab maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, Summer 2026
   !!
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE par_oce,  ONLY: jpi, jpj, Nis0, Nie0, Njs0, Nje0, nn_hls, ntsi, ntsj, ntei, ntej
   USE lib_mpp,  ONLY: ctl_stop         ! distribued memory computing library
   USE phycst          ! physical constants
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
   USE par_ice,         ONLY: epsi10, epsi20
   USE in_out_manager , ONLY: ln_timing, lwp
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_stab

   INTEGER , PARAMETER ::   nb_iter = 8   ! number of itterations, pretty big because usually stable over sea-ice
   !                                      !   => and "stable" means that more itterations may be required

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_stab( zt, zu, pts_i, pt_zt, pqs_i, pq_zt, pU_zu, &
      &                      pCdN, pChN, pCeN,                          &
      &                      pCd_i, pCh_i, pCe_i, pt_zu_i, pq_zu_i ) !, &
      !&                      xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_stab  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to:
      !!   Andreas, E.L., Jordan, R.E. & Makshtas, A.P. Parameterizing turbulent exchange over sea ice: the ice station weddell results.
      !!   Boundary-Layer Meteorology 114, 439–460 (2005). https://doi.org/10.1007/s10546-004-1414-7
      !!
      !!           If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!           Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  pts_i  : surface temperature of sea-ice                         [K]
      !!    *  pt_zt : potential air temperature at zt                         [K]
      !!    *  pqs_i  : saturation specific humidity at temp. pts_i over ice    [kg/kg]
      !!    *  pq_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  pU_zu : scalar wind speed at zu                                 [m/s]
      !!    * pCdN     : neutral-stability drag coefficient
      !!    * pChN     : neutral-stability sensible heat coefficient
      !!    * pCeN     : neutral-stability evaporation coefficient
      !!
      !! OUTPUT :
      !! --------
      !!    *  pCd_i   : drag coefficient over sea-ice
      !!    *  pCh_i   : sensible heat coefficient over sea-ice
      !!    *  pCe_i   : sublimation coefficient over sea-ice
      !!    *  pt_zu_i : pot. air temp. adjusted at zu over sea-ice             [K]
      !!    *  pq_zu_i : spec. hum. of air adjusted at zu over sea-ice          [kg/kg]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0     : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star : return u* the friction velocity                    [m/s]
      !!    * xL      : return the Obukhov length                          [m]
      !!    * xUN10   : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, July 2023 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in )                     :: zt    ! height for pt_zt and pq_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for pU_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: pts_i  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: pt_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: pqs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: pq_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: pU_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in )                     :: pCdN, pChN, pCeN  ! neutral coefficients             [-]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: pCd_i  ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: pCh_i  ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: pCe_i  ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: pt_zu_i ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: pq_zu_i ! spec. humidity adjusted at zu           [kg/kg]
      !!----------------------------------------------------------------------------------
      INTEGER  :: ji, jj, jit
      REAL(wp) :: zUbzu, z1_vkarmn
      REAL(wp) :: zdum, zdvm, z1_L, zsqrtCD
      REAL(wp) :: zdt_zu, zdq_zu
      REAL(wp) :: zus, zts, zqs
      REAL(wp) :: zzeta_u, zzeta_t           ! stability parameter at height zu
      REAL(wp) :: zCD, zCE, zCH, zt_0, zq_0, zt_zu, zq_zu, zt_zt, zq_zt
      REAL(wp) :: zsqrtCDN, z1_sqrtCDN, zlog1, zlog2
      LOGICAL  :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_stab@sbcblk_algo_ice_stab.f90'
      !!----------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('turb_ice_stab')
      !$acc data present( pts_i, pt_zt, pqs_i, pq_zt, pU_zu, pCdN, pChN, pCeN, pCd_i, pCh_i, pCe_i, pt_zu_i, pq_zu_i )

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      ! Optimization:
      zlog1    = LOG(zt/zu)
      zlog2    = LOG(zu/10._wp)

      z1_vkarmn  = 1._wp / vkarmn

      zsqrtCDN   = SQRT( pCdN )
      z1_sqrtCDN = 1._wp / zsqrtCDN

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            !! Theta & q at the air-sea interface:
            zt_0 = pts_i(ji,jj)
            zq_0 = pqs_i(ji,jj)

            !! Scalar wind speed cannot be below 0.2 m/s
            zUbzu = MAX( pU_zu(ji,jj), wspd_thrshld_ice )

            zt_zt = pt_zt(ji,jj)
            zq_zt = pq_zt(ji,jj)

            !! First guess of temperature and humidity at height zu:
            zt_zu = MAX( zt_zt , 100._wp )   ! who knows what's given on masked-continental regions...
            zq_zu = MAX( zq_zt , epsi10  )   !               "

            !! First guess of transfer coefficients:
            zCD = pCdN
            zCH = pChN
            zCE = pCeN

            zsqrtCD  = SQRT(zCD)             ! == u*/zUbzu

            !#DEBUG:
#if defined key_verbose
            IF( lwp .AND. (ji==4).AND.(jj==4) ) THEN
               zdt_zu = zt_zu - zt_0
               PRINT *, ' BEFORE LOOP:'
               PRINT *, '  zt_zu =', zt_zu
               PRINT *, '  zt_0 =', zt_0
               PRINT *, '  zdt_zu =', zdt_zu
               PRINT *, ''
            ENDIF
#endif
            !#DEBUG.

            !! ITERATION BLOCK

            !$acc loop seq
            DO jit = 1, nb_iter

               !! Air-Ice differences:
               zdt_zu = zt_zu - zt_0   !  zdt_zu = SIGN( MAX(ABS(zdt_zu),1.E-6_wp), zdt_zu )   !RM 2nd part!!!
               zdq_zu = zq_zu - zq_0   !  zdq_zu = SIGN( MAX(ABS(zdq_zu),1.E-9_wp), zdq_zu )

               !! Now we can get the urbulent scales:
               zus = zsqrtCD * zUbzu
               zdum  = 1._wp / MAX( zsqrtCD, epsi10 )   ! == zUbzu/u*
               zts =      zCH  * zdt_zu * zdum
               zqs =      zCE  * zdq_zu * zdum

               !!Inverse of Obukov length (1/L) :
               z1_L = One_on_L(zt_zu, zq_zu, zus, zts, zqs)  ! 1/L == 1/[Obukhov length]
               z1_L = SIGN( MIN(ABS(z1_L),200._wp), z1_L ) ! (prevents FPE from stupid values from masked region later on...)

               !! Stability parameters "zeta" :
               zzeta_u = zu*z1_L
               zzeta_t = zt*z1_L
               zzeta_u = SIGN( MIN(ABS(zzeta_u),50.0_wp), zzeta_u )
               zzeta_t = SIGN( MIN(ABS(zzeta_t),50.0_wp), zzeta_t )
               
               !! Update C_D:
               zdum = 1._wp + zsqrtCDN*z1_vkarmn*(zlog2 - psi_m_ice(zzeta_u))
               zCD  = MAX( pCdN / ( zdum*zdum ) , Cx_min )
               !zCD  = MIN( MAX( pCdN / ( zdum*zdum ), Cx_min ) , 1.9E-3_wp ) ! capped version

               zsqrtCD  = SQRT(zCD)             ! == u*/zUbzu

               !! Update C_H and C_E
               zdum = ( zlog2 - psi_h_ice(zzeta_u) ) * z1_vkarmn * z1_sqrtCDN
               zdvm = zsqrtCD * z1_sqrtCDN
               zCH  = MAX( pChN*zdvm / ( 1._wp + pChN*zdum ) , Cx_min )
               zCE  = MAX( pCeN*zdvm / ( 1._wp + pCeN*zdum ) , Cx_min )
               !zCH  = MIN( MAX( pChN*zdvm / ( 1._wp + pChN*zdum ) , Cx_min ) , 1.999E-3_wp ) ! capped version
               !zCE  = MIN( MAX( pCeN*zdvm / ( 1._wp + pCeN*zdum ) , Cx_min ) , 1.999E-3_wp ) ! capped version

               !! Re-updating temperature and humidity at zu if zt /= zu :
               zdum  = z1_vkarmn * ( psi_h_ice(zzeta_u) - psi_h_ice(zzeta_t) + zlog1 )
               zt_zu = MERGE( zt_zu  ,             zt_zt - zts*zdum   ,  l_zt_equal_zu )
               zq_zu = MERGE( zq_zu  ,  MAX(0._wp, zq_zt - zqs*zdum)  ,  l_zt_equal_zu )

            END DO !DO jit = 1, nb_iter

            pCd_i(ji,jj) = zCD
            pCh_i(ji,jj) = zCH
            pCe_i(ji,jj) = zCE

            pt_zu_i(ji,jj) = zt_zu
            pq_zu_i(ji,jj) = zq_zu

            !#DEBUG:
#if defined key_verbose
            IF( lwp .AND. (ji==4).AND.(jj==4) ) THEN
               zdt_zu = pt_zu_i(ji,jj) - pts_i(ji,jj)
               PRINT *, ' AFTER LOOP:'
               PRINT *, '  pt_zu_i =', pt_zu_i(ji,jj)
               PRINT *, '  pts_i =', pts_i(ji,jj)
               PRINT *, '  zdt_zu =', zdt_zu
               PRINT *, '  pCd_i =', pCd_i(ji,jj)
               PRINT *, ''
            ENDIF
#endif
            !#DEBUG.


         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('turb_ice_stab')

   END SUBROUTINE turb_ice_stab




   FUNCTION psi_m_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      !$acc routine
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: psi_m_ice
      REAL(wp), INTENT(in) :: pzeta
      !!----------------------------------------------------------------------------------
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta
      !
      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u = LOG( 0.5_wp*(1._wp + zx*zx) ) + 2._wp*LOG( 0.5_wp*(1._wp + zx) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi  ! Eq.(30) Jordan et al. 1999:

      ! Stable stratification:
      !zpsi_s = 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35_wp*zta) + 10.7_wp  ! Eq.(33) Jordan et al. 1999
      zpsi_s = 1._wp + 6.5_wp*zta * ABS(1._wp + zta)**0.3333333_wp / ( 1.3_wp + zta )     ! Eq.(9.a) Grachev et al. 2007

      !! Combine:
      psi_m_ice = MERGE( zpsi_u  ,  -1._wp*zpsi_s  ,  zta < 0._wp )   ! Unstable <= zta < 0)  |  Stable <= zta > 0)
      !
   END FUNCTION psi_m_ice


   FUNCTION psi_h_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for
      !!             temperature and humidity
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      !$acc routine
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: psi_h_ice
      REAL(wp), INTENT(in) :: pzeta
      !!----------------------------------------------------------------------------------
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s
      !!----------------------------------------------------------------------------------
      zta = pzeta

      ! Unstable stratification:
      zx = ABS(1._wp - 16._wp*zta)**0.25_wp              !  (16 here, not 15!)

      zpsi_u =   2._wp*LOG( 0.5_wp * (1._wp + zx*zx) )  ! Eq.(31) Jordan et al. 1999

      ! Stable stratification (identical to Psi_m!):
      !zpsi_s = 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35_wp*zta) + 10.7_wp       ! Eq.(33) Jordan et al. 1999
      zpsi_s = 1._wp + 5._wp*zta*( 1._wp + zta ) / ( 1._wp + 3._wp*zta + zta*zta )       ! Eq.(9.b) Grachev et al. 2007

      !! Combine:
      psi_h_ice = MERGE( zpsi_u  ,  -1._wp*zpsi_s  ,  zta < 0._wp )   ! Unstable <= zta < 0)  |  Stable <= zta > 0)

   END FUNCTION psi_h_ice

   !!======================================================================
END MODULE sbcblk_algo_ice_stab
