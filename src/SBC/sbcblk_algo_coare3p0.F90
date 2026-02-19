MODULE sbcblk_algo_coare3p0
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_coare3p0  ***
   !!
   !!       Computes turbulent components of surface fluxes
   !!         according the formulation/param. of Fairall et al, 2003
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ubzu
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_coare3p0 maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk)
   !!
   !!   When using AeroBulk to produce scientific work, please acknowledge with the following citation:
   !!
   !!   Brodeau, L., B. Barnier, S. Gulev, and C. Woods, 2016: Climatologically
   !!   significant effects of some approximations in the bulk parameterizations of
   !!   turbulent air-sea fluxes. J. Phys. Oceanogr., doi:10.1175/JPO-D-16-0169.1.
   !!
   !!
   !!            Author: Laurent Brodeau, 2016
   !!
   !!=====================================================================
   !! History :  4.0  ! 2016-02  (L.Brodeau)   Original code
   !!            4.2  ! 2020-12  (L. Brodeau) Introduction of various air-ice bulk parameterizations + improvements
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_coare3p0  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE phycst          ! physical constants
   USE lib_mpp,        ONLY: ctl_stop         ! distribued memory computing library
   USE in_out_manager, ONLY: nit000, nitend, ln_timing  ! I/O manager
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
   USE sbcblk_coare_util, ONLY : first_guess_coare
   USE sbcblk_skin_coare ! cool-skin/warm layer scheme
   USE oss_nnq , ONLY : e3t_m
   USE ossskin , ONLY : oss_skin_alloc, oss_skin_dealloc, dT_cs, dT_wl, Hz_wl

   USE timing         ! Timing

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TURB_COARE3P0

   !! COARE own values for given constants:
   REAL(wp), PARAMETER :: zi0   = 600._wp     ! scale height of the atmospheric boundary layer...
   REAL(wp), PARAMETER :: Beta0 =  1.25_wp    ! gustiness parameter
   REAL(wp), PARAMETER :: zeta_abs_max = 50._wp

CONTAINS

   SUBROUTINE turb_coare3p0( kt, zt, zu, pSST, pT_s, pt_zt, pq_s, pq_zt, pU_zu, l_use_cs, l_use_wl, &
      &                       pCd, pCh, pCe, pt_zu, pq_zu, pUbzu,                                   &
      &                       nb_iter,                                                              & ! optional input
      &                       pQsw, prad_lw, pslp )                                                   ! opt. cool-skin & warm-layer
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_coare3p0  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Fairall et al. (2003)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to pT_s
      !!                if the net shortwave flux at the surface (pQsw), the downwelling longwave
      !!                radiative fluxes at the surface (prad_lw), and the sea-leve pressure (pslp)
      !!                are provided as (optional) arguments!
      !!
      !! INPUT :
      !! -------
      !!    *  kt    : current time step (starts at 1)
      !!    *  zt    : height for temperature and spec. hum. of air            [m]
      !!    *  zu    : height for wind speed (usually 10m)                     [m]
      !!    *  pSST  : bulk SST                                                [deg.C]
      !!    *  pt_zt : potential air temperature at zt                         [K]
      !!    *  pq_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  pU_zu : scalar wind speed at zu                                 [m/s]
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  pT_s  : surface skin temperature                               [deg.C]
      !!    *  pq_s  : saturation specific humidity at temp. pT_s             [kg/kg]
      !!  ==> these 2 are identical `pSST` & `q_sat(pSST)` when CSWL not used !!!
      !!
      !! OPTIONAL INPUT:
      !! ---------------
      !!    *  pQsw    : net solar flux (after albedo) at the surface (>0)     [W/m^2]
      !!    *  prad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  pslp    : sea-level pressure                                    [Pa]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * pdT_cs  : SST increment "dT" for cool-skin correction           [K]
      !!    * pdT_wl  : SST increment "dT" for warm-layer correction          [K]
      !!    * pHz_wl  : thickness of warm-layer                               [m]
      !!
      !! OUTPUT :
      !! --------
      !!    *  pCd     : drag coefficient
      !!    *  pCh     : sensible heat coefficient
      !!    *  pCe     : evaporation coefficient
      !!    *  pt_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  pq_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  pUbzu   : bulk wind speed at zu                                 [m/s]
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in   )                     ::   kt       ! current time step
      REAL(wp), INTENT(in   )                     ::   zt       ! height for pt_zt and pq_zt                [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for pU_zu                          [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pSST     ! BULK SST                              [deg.C]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   pT_s     ! SKIN SST                              [deg.C]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pt_zt    ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   pq_s     ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pq_zt    ! specific air humidity at zt           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pU_zu    ! relative wind module at zu              [m/s]
      LOGICAL , INTENT(in   )                     ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in   )                     ::   l_use_wl ! use the warm-layer parameterization
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pCd      ! transfer coefficient for momentum         [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pCh      ! transfer coefficient for sensible heat    [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pCe      ! transfert coefficient for evaporation     [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pt_zu    ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pq_zu    ! spec. humidity adjusted at zu         [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pUbzu    ! bulk wind speed at zu                   [m/s]
      !
      INTEGER , INTENT(in   ), OPTIONAL           :: nb_iter    ! number of iterations
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   pQsw      !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   prad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   pslp      !             [Pa]
      !!----------------------------------------------------------------------------------
      LOGICAL  :: l_skin
      INTEGER  :: nbit
      INTEGER  :: ji, jj, jit
      REAL(wp) :: zdum, zm_ztzu                ! => `1.` if `zu /= zt`, `0.` otherwize
      REAL(wp) :: zdt, zdq, zus, zus2, zUzu, zts, zqs, zNu_a, z1oL, zdT_cs
      REAL(wp) :: zUn10, zgust2, zz0, zz0t, zzta_u, zzta_t
      REAL(wp) :: zlog_10, zlog_zu, zlog_zt, zlog_ztu, zlog_z0, zlog_z0t
      REAL(wp) :: zQns, zQlat, zTau
      REAL(wp) :: zSST, zT_s, zq_s, zubzu, zt_zt, zq_zt, zt_zu, zq_zu
      REAL(wp) :: ztmp0, ztmp1
      !
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_coare3p0@sbcblk_algo_coare3p0'
      !!----------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('turb_coare3p0')
      !$acc data present( pSST, pT_s, pt_zt, pq_s, pq_zt, pU_zu, pCd, pCh, pCe, pt_zu, pq_zu, pUbzu )


      !IF( kt == nit000 ) CALL oss_skin_alloc(l_use_cs, l_use_wl)

      nbit = nb_iter0
      IF( PRESENT(nb_iter) ) nbit = nb_iter

      zm_ztzu = MERGE( 0._wp, 1._wp,  ABS(zu - zt) < 0.01_wp )  ! => `1.` if `zu /= zt`, `0.` otherwize

      !! Initializations for cool skin and warm layer:
      IF( l_use_cs .AND. (.NOT.(PRESENT(pQsw) .AND. PRESENT(prad_lw) .AND. PRESENT(pslp))) ) &
         & CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , ' use of `l_use_cs` (cool-skin) not supporterted yet! FIXME!!!' )
      !&   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide pQsw, prad_lw & pslp to use cool-skin param!' )

      IF( l_use_wl .AND. (.NOT.(PRESENT(pQsw) .AND. PRESENT(prad_lw) .AND. PRESENT(pslp))) ) &
         CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , ' use of `l_use_wl` (warm-layer) not supporterted yet! FIXME!!!' )
      !&   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide pQsw, prad_lw & pslp to use warm-layer param!' )


      l_skin = l_use_cs .OR. l_use_wl

      !! Constants:
      zlog_10  = LOG(10._wp)
      zlog_zt  = LOG(zt)
      zlog_zu  = LOG(zu)
      zlog_ztu = LOG(zt/zu)

      !$acc parallel loop collapse(2)
      DO jj = Njs0, Nje0
         DO ji = Nis0, Nie0

            zSST  = pSST(ji,jj) + rt0   ! => to Kelvin
            zT_s  = pT_s(ji,jj) + rt0   ! => to Kelvin
            zt_zt = pt_zt(ji,jj)
            zq_s  = pq_s(ji,jj)
            zq_zt = pq_zt(ji,jj)
            zUzu  = pU_zu(ji,jj)

            CALL first_guess_coare( zt, zu, zT_s, zt_zt, zq_s, zq_zt, zUzu, &
               &                    charn_coare3p0(zUzu),  zus, zts, zqs, &
               &                    zt_zu, zq_zu, zUbzu,  pz0=zz0 )

            zlog_z0 = LOG(zz0)
            znu_a   = visc_air(zt_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = zt_zu - zT_s ;   zdt = SIGN( MAX(ABS(zdt),1.E-09_wp), zdt )
            zdq = zq_zu - zq_s ;   zdq = SIGN( MAX(ABS(zdq),1.E-12_wp), zdq )


            !! ITERATION BLOCK
            !$acc loop seq
            DO jit = 1, nb_iter

               zus2    = zus*zus   ! u*^2

               !!Inverse of Obukov length (1/L) :
               z1oL = One_on_L(zt_zu, zq_zu, zus, zts, zqs)  ! 1/L == 1/[Obukhov length]
               z1oL = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! 1/L (prevents FPE from stupid values from masked region later on...)

               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Fairall et al. 2003, Eq.8):
               zgust2 = Beta0*Beta0*zus2*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution, zgust2 == Ug^2
               !!   ! Only true when unstable (L<0) => when z1oL < 0 => explains "-" before zi0
               zUbzu = MAX(SQRT(zUzu*zUzu + zgust2), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents pUbzu to be 0 in stable case when pU_zu=0.

               !! Stability parameters:
               zzta_u = zu*z1oL
               zzta_u = SIGN( MIN(ABS(zzta_u),zeta_abs_max), zzta_u )
               zzta_t = zt*z1oL
               zzta_t = SIGN( MIN(ABS(zzta_t),zeta_abs_max), zzta_t )

               !! Roughness lengthes z0, z0t (z0q = z0t):
               zUn10 = zus/vkarmn*(zlog_10 - zlog_z0)         ! Neutral wind speed at 10m
               zz0    = charn_coare3p0(zUn10)*zus2/grav + 0.11_wp*znu_a/zus ! Roughness length (eq.6)
               zz0     = MIN( MAX(ABS(zz0), 1.E-9) , 1._wp )  ! (prevents FPE from stupid values from masked region later on)
               zlog_z0 = LOG(zz0)

               ztmp1 = ( znu_a / (zz0*zus) )**0.6_wp         ! (1./Re_r)^0.6 (Re_r: roughness Reynolds number) COARE 3.0 - specific!
               zz0t   = MIN( 1.1E-4_wp , 5.5E-5_wp*ztmp1 )   ! Scalar roughness for temp. and q (eq.28) #LB: some use 1.15 not 1.1 !!!
               zz0t   = MIN( MAX(ABS(zz0t), 1.E-9) , 1._wp ) ! (prevents FPE from stupid values from masked region later on)
               zlog_z0t = LOG(zz0t)

               !! Turbulent scales at zu:
               ztmp0   = psi_h_coare(zzta_u)
               ztmp1   = vkarmn/(zlog_zu - zlog_z0t - ztmp0) ! #LB: in ztmp0, some use psi_h_coare(zzta_t) rather than psi_h_coare(zzta_t) ???

               zts = zdt*ztmp1
               zqs = zdq*ztmp1
               zus = MAX( zUbzu*vkarmn/(zlog_zu - zlog_z0 - psi_m_coare(zzta_u)) , 1.E-9 )

               !! Adjusting temperature and humidity at zu if required by `zm_ztzu`:
               ztmp1 = zlog_zt - zlog_zu + ztmp0 - psi_h_coare(zzta_t)
               zt_zu = zt_zt - zm_ztzu*zts/vkarmn*ztmp1
               zq_zu = zq_zt - zm_ztzu*zqs/vkarmn*ztmp1

               !IF( l_use_cs ) THEN
               !   !! Cool-skin contribution
               !   CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, zt_zu, zq_zu, zus, zts, zqs, &
               !      &                   zUzu, zUbzu, pslp(ji,jj), prad_lw(ji,jj), zQns, zTau, Qlat=zQlat )
               !
               !   CALL CS_COARE( pQsw(ji,jj), zQns, zus, zSST, zQlat, zdT_cs )
               !   IF( PRESENT(pdT_cs) ) pdT_cs(ji,jj) = zdT_cs
               !   zT_s = zSST + zdT_cs
               !   IF( l_use_wl ) zT_s = zT_s + dT_wl(ji,jj)
               !   zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), pslp(ji,jj))
               !ENDIF

               !IF( l_use_wl ) THEN
               !   !! Warm-layer contribution
               !   CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, zt_zu, zq_zu, zus, zts, zqs, &
               !      &                   zUzu, zUbzu, pslp(ji,jj), prad_lw(ji,jj),  zQns, zTau)
               !   !! In WL_COARE or , Tau_ac and Qnt_ac must be updated at the final itteration step => add a flag to do this!
               !   CALL WL_COARE( ji, jj, pQsw(ji,jj), zQns, zTau, zSST, plong(ji,jj), isecday_utc, MOD(nb_iter,jit) )
               !   !! Updating pT_s and pq_s !!!
               !   zT_s = zSST + dT_wl(ji,jj)
               !   IF( l_use_cs ) zT_s = zT_s + zdT_cs
               !   zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), pslp(ji,jj))
               !ENDIF

               zdt = zt_zu - zT_s ;  zdt = SIGN( MAX(ABS(zdt),1.E-09_wp), zdt )
               zdq = zq_zu - zq_s ;  zdq = SIGN( MAX(ABS(zdq),1.E-12_wp), zdq )

            END DO !DO jit = 1, nb_iter

            !! Update arrays that are returned by the routine:
            IF( l_skin ) THEN
               pT_s(ji,jj)  = zT_s - rt0  ! back to deg.C
               pq_s(ji,jj)  = zq_s
            ENDIF
            pt_zu(ji,jj) = zt_zu
            pq_zu(ji,jj) = zq_zu
            pUbzu(ji,jj) = zUbzu


            ! compute transfer coefficients at zu :
            ztmp0 = zus/zUbzu
            pCd(ji,jj) = MAX( ztmp0*ztmp0   , Cx_min )
            pCh(ji,jj) = MAX( ztmp0*zts/zdt , Cx_min )
            pCe(ji,jj) = MAX( ztmp0*zqs/zdq , Cx_min )

         END DO !DO ji = Nis0, Nie0
      END DO !DO jj = Njs0, Nje0
      !$acc end parallel loop

      !$acc end data

      !IF( kt == nitend ) CALL oss_skin_dealloc( l_use_wl )

      IF( ln_timing )   CALL timing_stop('turb_coare3p0')

   END SUBROUTINE turb_coare3p0


   !!===============================================================================================
   FUNCTION charn_coare3p0( pwnd )
      !$acc routine
      !!-------------------------------------------------------------------
      !! Compute the Charnock parameter as a function of the wind speed
      !!
      !! (Fairall et al., 2003 p.577-578)
      !!
      !! Wind below 10 m/s :  alfa = 0.011
      !! Wind between 10 and 18 m/s : linear increase from 0.011 to 0.018
      !! Wind greater than 18 m/s :  alfa = 0.018
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk  (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pwnd   ! wind speed
      REAL(wp) :: charn_coare3p0
      !
      REAL(wp) :: zw, zgt10, zgt18
      !!-------------------------------------------------------------------
      zw = pwnd   ! wind speed
      !!
      !! Charnock's constant, increases with the wind :
      zgt10 = 0.5_wp + SIGN(0.5_wp,(zw - 10._wp))  ! If zw<10. --> 0, else --> 1
      zgt18 = 0.5_wp + SIGN(0.5_wp,(zw - 18._wp))  ! If zw<18. --> 0, else --> 1
      !
      charn_coare3p0 =  (1. - zgt10)*0.011    &    ! wind is lower than 10 m/s
         &              + zgt10*((1. - zgt18)*(0.011 + (0.018 - 0.011) &
         &              *(zw - 10.)/(18. - 10.)) + zgt18*( 0.018 ) )    ! Hare et al. (1999)
      !!
   END FUNCTION charn_coare3p0



   !!===============================================================================================
   FUNCTION psi_m_coare( pzeta )
      !$acc routine
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!             COARE 3.0, Fairall et al. 2003
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!       Stability function for wind speed and scalars matching Kansas and free
      !!       convection forms with weighting f convective form, follows Fairall et
      !!       al (1996) with profile constants from Grachev et al (2000) BLM stable
      !!       form from Beljaars and Holtslag (1991)
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp) :: psi_m_coare
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------------------------
      zphi_m = ABS(1. - 15.*pzeta)**.25    !!Kansas unstable
      !
      zpsi_k = 2.*LOG((1. + zphi_m)/2.) + LOG((1. + zphi_m*zphi_m)/2.)   &
         & - 2.*ATAN(zphi_m) + 0.5*rpi
      !
      zphi_c = ABS(1. - 10.15*pzeta)**.3333                   !!Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &     - 1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp, 0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_m_coare = (1. - zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) & ! (pzeta < 0)
         &           -   zstb  * ( 1. + 1.*pzeta     &                ! (pzeta > 0)
         &                          + 0.6667*(pzeta - 14.28)/EXP(zc) + 8.525 )  !     "
      !!
   END FUNCTION psi_m_coare

   !!===============================================================================================


   !!===============================================================================================
   FUNCTION psi_h_coare( pzeta )
      !$acc routine
      !!---------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !! COARE 3.0, Fairall et al. 2003
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! Stability function for wind speed and scalars matching Kansas and free
      !! convection forms with weighting f convective form, follows Fairall et
      !! al (1996) with profile constants from Grachev et al (2000) BLM stable
      !! form from Beljaars and Holtslag (1991)
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------
      REAL(wp) :: psi_h_coare
      REAL(wp), INTENT(in) :: pzeta
      !!
      REAL(wp) :: zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstb
      !!----------------------------------------------------------------
      zphi_h = (ABS(1. - 15.*pzeta))**.5  !! Kansas unstable   (zphi_h = zphi_m**2 when unstable, zphi_m when stable)
      !
      zpsi_k = 2.*LOG((1. + zphi_h)/2.)
      !
      zphi_c = (ABS(1. - 34.15*pzeta))**.3333   !! Convective
      !
      zpsi_c = 1.5*LOG((1. + zphi_c + zphi_c*zphi_c)/3.) &
         &    -1.7320508*ATAN((1. + 2.*zphi_c)/1.7320508) + 1.813799447
      !
      zf = pzeta*pzeta
      zf = zf/(1. + zf)
      zc = MIN(50._wp,0.35_wp*pzeta)
      zstb = 0.5 + SIGN(0.5_wp, pzeta)
      !
      psi_h_coare = (1.-zstb) * ( (1. - zf)*zpsi_k + zf*zpsi_c ) &
         &                  -zstb  * ( (ABS(1. + 2.*pzeta/3.))**1.5     &
         &                            + .6667*(pzeta - 14.28)/EXP(zc) + 8.525 )
      !!
   END FUNCTION psi_h_coare

   !!===============================================================================================

   !!======================================================================
END MODULE sbcblk_algo_coare3p0
