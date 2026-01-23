MODULE sbcblk_algo_ecmwf
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ecmwf  ***
   !!
   !!       Computes turbulent components of surface fluxes
   !!         according the formulation/param. of IFS of ECMWF (cycle 40r1)
   !!         based on IFS doc (avaible online on the ECMWF's website)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ubzu
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!
   !!       Routine turb_ecmwf maintained and developed in AeroBulk
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
   !! History :  3.6  !  2016-02  (L.Brodeau) successor of old turb_ncar of former sbcblk_core.F90
   !!            4.2  !  2020-12  (L. Brodeau) Introduction of various air-ice bulk parameterizations + improvements
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ecmwf  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE phycst          ! physical constants
   USE lib_mpp,        ONLY: ctl_stop         ! distribued memory computing library
   USE in_out_manager, ONLY: nit000, nitend, ln_timing  ! I/O manager
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
   USE sbcblk_coare_util, ONLY : first_guess_coare
   USE sbcblk_skin_ecmwf ! cool-skin/warm layer scheme !LB
   USE oss_nnq , ONLY : sbcblk_cs_wl_alloc, sbcblk_cs_wl_dealloc, e3t_m, dT_cs, dT_wl, Hz_wl

   USE timing         ! Timing

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TURB_ECMWF

   !! ECMWF own values for given constants, taken form IFS documentation...
   REAL(wp), PARAMETER, PUBLIC :: charn0_ecmwf = 0.018_wp    ! Charnock constant (pretty high value here !!!
   !                                          !    =>  Usually 0.011 for moderate winds)
   REAL(wp), PARAMETER ::   zi0     = 1000.   ! scale height of the atmospheric boundary layer...1
   REAL(wp), PARAMETER ::   Beta0    = 1.     ! gustiness parameter ( = 1.25 in COAREv3)
   REAL(wp), PARAMETER ::   alpha_M = 0.11    ! For roughness length (smooth surface term)
   REAL(wp), PARAMETER ::   alpha_H = 0.40    ! (Chapter 3, p.34, IFS doc Cy31r1)
   REAL(wp), PARAMETER ::   alpha_Q = 0.62    !


CONTAINS

   SUBROUTINE turb_ecmwf( kt, zt, zu, T_s, t_zt, q_s, q_zt, U_zu, l_use_cs, l_use_wl, &
      &                      Cd, Ch, Ce, t_zu, q_zu, Ubzu,                            &
      &                      nb_iter,                                         & ! optional input
      &                      Qsw, rad_lw, slp )                                 ! optionals for warm-layer only
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ecmwf  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to IFS doc. (cycle 45r1)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !!                Applies the cool-skin warm-layer correction of the SST to T_s
      !!                if the net shortwave flux at the surface (Qsw), the downwelling longwave
      !!                radiative fluxes at the surface (rad_lw), and the sea-leve pressure (slp)
      !!                are provided as (optional) arguments!
      !!
      !! INPUT :
      !! -------
      !!    *  kt   : current time step (starts at 1)
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * l_use_cs : use the cool-skin parameterization
      !!    * l_use_wl : use the warm-layer parameterization
      !!
      !! INPUT/OUTPUT:
      !! -------------
      !!    *  T_s  : always "bulk SST" as input                              [K]
      !!              -> unchanged "bulk SST" as output if CSWL not used      [K]
      !!              -> skin temperature as output if CSWL used              [K]
      !!
      !!    *  q_s  : SSQ aka saturation specific humidity at temp. T_s       [kg/kg]
      !!              -> doesn't need to be given a value if skin temp computed (in case l_use_cs=True or l_use_wl=True)
      !!              -> MUST be given the correct value if not computing skint temp. (in case l_use_cs=False or l_use_wl=False)
      !!
      !! OPTIONAL INPUT:
      !! ---------------
      !!    *  Qsw    : net solar flux (after albedo) at the surface (>0)     [W/m^2]
      !!    *  rad_lw : downwelling longwave radiation at the surface  (>0)   [W/m^2]
      !!    *  slp    : sea-level pressure                                    [Pa]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * pdT_cs  : SST increment "dT" for cool-skin correction           [K]
      !!    * pdT_wl  : SST increment "dT" for warm-layer correction          [K]
      !!    * pHz_wl  : thickness of warm-layer                               [m]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ubzu   : bulk wind speed at zu                                 [m/s]
      !!
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      INTEGER,  INTENT(in   )                     ::   kt       ! current time step
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                  [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                           [m]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   T_s      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   q_s      ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu              [m/s]
      LOGICAL , INTENT(in   )                     ::   l_use_cs ! use the cool-skin parameterization
      LOGICAL , INTENT(in   )                     ::   l_use_wl ! use the warm-layer parameterization
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat    [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation     [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu         [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ubzu     ! bulk wind speed at zu                   [m/s]
      INTEGER , INTENT(in   ), OPTIONAL           :: nb_iter    ! number of iterations
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   Qsw      !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   rad_lw   !             [W/m^2]
      REAL(wp), INTENT(in   ), OPTIONAL, DIMENSION(jpi,jpj) ::   slp      !             [Pa]
      !!----------------------------------------------------------------------------------
      INTEGER  :: ji, jj, nbit, jit
      REAL(wp) :: zm_ztzu                   ! => `1.` if `zu /= zt`, `0.` otherwize
      REAL(wp), DIMENSION(jpi,jpj) :: zTsblk  ! to back up the initial bulk SST
      REAL(wp) :: zdt, zdq, zus, zus2, zUzu, zts, zqs, zNu_a, z1oL, zdT_cs
      REAL(wp) :: zRib, zpsi_m_u, zpsi_h_u, zpsi_h_t
      REAL(wp) :: zz0, zz0t, zz0q, zpsi_h_z0t, zpsi_h_z0q, zpsi_m_z0, zzeta_u, zzeta_t
      REAL(wp) :: zlog_10, zlog_ztu, zlog_z0, zlog_zu, zlog_z0t, zlog_z0q
      REAL(wp) :: zFm, zFh, zFq, zQns
      REAL(wp) :: zSST, zT_s, zq_s, zubzu, zt_zt, zq_zt
      REAL(wp) :: ztmp0, ztmp1, zqzu, ztzu
      !
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ecmwf@sbcblk_algo_ecmwf.F90'
      !!----------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('turb_ecmwf')
      !$acc data present( T_s, t_zt, q_s, q_zt, U_zu, Cd, Ch, Ce, t_zu, q_zu, Ubzu ) create(zTsblk)

      
      !IF( kt == nit000 ) THEN sbcblk_cs_wl_alloc(l_use_cs, l_use_wl)


      nbit = nb_iter0
      IF( PRESENT(nb_iter) ) nbit = nb_iter

      zm_ztzu = 0.5_wp + SIGN( 0.5_wp, ABS(zu-zt) - 0.01_wp ) ! => `1.` if `zu /= zt`, `0.` otherwize

      !! Initializations for cool skin and warm layer:
      IF( l_use_cs .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         & CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , ' use of `l_use_cs` (cool-skin) not supporterted yet! FIXME!!!' )
         !&   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use cool-skin param!' )

      IF( l_use_wl .AND. (.NOT.(PRESENT(Qsw) .AND. PRESENT(rad_lw) .AND. PRESENT(slp))) ) &
         CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , ' use of `l_use_wl` (warm-layer) not supporterted yet! FIXME!!!' )
         !&   CALL ctl_stop( '['//TRIM(crtnm)//'] => ' , 'you need to provide Qsw, rad_lw & slp to use warm-layer param!' )

      !! Constants:
      zlog_10  = LOG(10._wp)
      zlog_zu  = LOG(zu)
      zlog_ztu = LOG(zt/zu)

      ! Backing up the bulk SST (`T_s` contains the bulk SST as input)
      !$acc parallel loop collapse(2)
      DO jj = Njs0, Nje0
         DO ji = Nis0, Nie0
            zTsblk(ji,jj) = T_s(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      IF( l_use_cs ) THEN         
         !*acc parallel loop collapse(2) present(T_s, q_s, slp)
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               zT_s       = T_s(ji,jj) - 0.25_wp   ! First guess of cool-skin correction
               T_s(ji,jj) = zT_s
               q_s(ji,jj) = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), slp(ji,jj)) ! decrease `q_s` accordingly
            END DO
         END DO
         !*acc end parallel loop
      ENDIF

      !$acc parallel loop collapse(2)
      DO jj = Njs0, Nje0
         DO ji = Nis0, Nie0

            zSST  = zTsblk(ji,jj)
            zT_s  = T_s(ji,jj)
            zt_zt = t_zt(ji,jj)
            zq_s  = q_s(ji,jj)
            zq_zt = q_zt(ji,jj)
            zUzu  = U_zu(ji,jj)

            CALL first_guess_coare( zt, zu, zT_s, zt_zt, zq_s, zq_zt, zUzu, &
               &                    charn0_ecmwf,  zus, zts, zqs, &
               &                    ztzu, zqzu, zubzu,  pz0=zz0 )

            zlog_z0 = LOG(zz0)
            znu_a   = visc_air(zt_zt) ! Air viscosity (m^2/s) at zt given from temperature in (K)

            !! Pot. temp. difference (and we don't want it to be 0!)
            zdt = ztzu - zT_s ;   zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
            zdq = zqzu - zq_s ;   zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )

            !! First guess of inverse of Obukov length (1/L) :
            z1oL    = One_on_L( ztzu, zqzu, zus, zts, zqs )
            zzeta_u = zu*z1oL
            zzeta_t = zt*z1oL

            zz0t    = MIN( MAX(ABS(  1._wp / ( 0.1_wp*EXP(vkarmn/(0.00115/( vkarmn/(zlog_10-zlog_z0) ))) )   ), 1.E-9) , 1._wp )
            zlog_z0t = LOG(zz0t)

            !! Functions such as  u* = Ubzu*vkarmn/zFm
            zFm = zlog_zu - zlog_z0  - psi_m_ecmwf_sclr(zzeta_u) + psi_m_ecmwf_sclr( zz0*z1oL)
            zpsi_h_u = psi_h_ecmwf_sclr(zzeta_u)
            zFh = zlog_zu - zlog_z0t - zpsi_h_u + psi_h_ecmwf_sclr(zz0t*z1oL)

            !! ITERATION BLOCK
            !$acc loop seq
            DO jit = 1, nb_iter

               !! Bulk Richardson Number at z=zu (Eq. 3.25)
               zRib = Ri_bulk( zu, zT_s, ztzu, zq_s, zqzu, zubzu ) ! Bulk Richardson Number (BRN)

               !! New estimate of the inverse of the Obukhon length (z1oL == zeta/zu) :
               z1oL = zRib*zFm*zFm/zFh / zu     ! From Eq. 3.23, Chap.3.2.3, IFS doc - Cy40r1
               !! Note: it is slightly different that the L we would get with the usual
               z1oL   = SIGN( MIN(ABS(z1oL),200._wp), z1oL ) ! (prevent FPE from stupid values from masked region later on...)

               zzeta_u  = zu*z1oL
               zpsi_m_u = psi_m_ecmwf_sclr(zzeta_u)
               zpsi_h_u = psi_h_ecmwf_sclr(zzeta_u)

               zzeta_t  = zt*z1oL
               zpsi_h_t = psi_h_ecmwf_sclr(zzeta_t)

               !! Update zFm with new z1oL:
               zFm = zlog_zu -zlog_z0 - zpsi_m_u + psi_m_ecmwf_sclr(zz0*z1oL) ! LB: should be "zu+z0" rather than "zu" alone, but z0 is tiny wrt zu!

               !! Need to update roughness lengthes:
               zus = zubzu*vkarmn/zFm
               zus2  = zus*zus
               ztmp0  = znu_a/zus
               zz0     = MIN( ABS( alpha_M*ztmp0 + charn0_ecmwf*zus2/grav ) , 0.001_wp)
               zz0t    = MIN( ABS( alpha_H*ztmp0                           ) , 0.001_wp)   ! eq.3.26, Chap.3, p.34, IFS doc - Cy31r1
               zz0q    = MIN( ABS( alpha_Q*ztmp0                           ) , 0.001_wp)

               zlog_z0  = LOG(zz0 )
               zlog_z0t = LOG(zz0t)
               zlog_z0q = LOG(zz0q)

               zpsi_m_z0   = psi_m_ecmwf_sclr(zz0 *z1oL)  ! LB: should be "zu+z0" rather than "zu" alone, but z0 is tiny wrt zu!
               zpsi_h_z0t  = psi_h_ecmwf_sclr(zz0t*z1oL)  !              "                           "
               zpsi_h_z0q  = psi_h_ecmwf_sclr(zz0q*z1oL)  !              "                           "


               !! Update wind at zu with convection-related wind gustiness in unstable conditions (Chap. 3.2, IFS doc - Cy40r1, Eq.3.17 and Eq.3.18 + Eq.3.8)
               ztmp0 = Beta0*Beta0*zus2*(MAX(-zi0*z1oL/vkarmn,0._wp))**(2._wp/3._wp) ! square of wind gustiness contribution  (combining Eq. 3.8 and 3.18, hap.3, IFS doc - Cy31r1)
               !!   ! Only true when unstable (L<0) => when zRib < 0 => explains "-" before zi0
               zubzu = MAX(SQRT(U_zu(ji,jj)*U_zu(ji,jj) + ztmp0), 0.2_wp)        ! include gustiness in bulk wind speed
               ! => 0.2 prevents Ubzu to be 0 in stable case when U_zu=0.


               !! Shifting temperature and humidity at zu if required by `zm_ztzu`:
               !IF( .NOT. l_zt_equal_zu ) THEN
               ztmp0  = zpsi_h_u - zpsi_h_z0t
               ztmp1  = vkarmn/(zlog_zu - zlog_z0t - ztmp0)
               zts = zdt*ztmp1
               ztmp1  = zlog_ztu + ztmp0 - zpsi_h_t + zpsi_h_z0t
               ztzu   = zt_zt - zm_ztzu*zts/vkarmn*ztmp1

               zpsi_h_z0q  = zpsi_h_z0q
               ztmp0  = zpsi_h_u - zpsi_h_z0q
               ztmp1  = vkarmn/(zlog_zu - zlog_z0q - ztmp0)
               zqs = zdq*ztmp1
               ztmp1  = zlog_ztu + ztmp0 - zpsi_h_t + zpsi_h_z0q
               zqzu   = MAX( zq_zt - zm_ztzu*zqs/vkarmn*ztmp1 , 0._wp )
               !ENDIF

               !! Updating because of updated z0 and z0t and new z1oL...
               zFm = zlog_zu - zlog_z0  - zpsi_m_u + zpsi_m_z0
               zFh = zlog_zu - zlog_z0t - zpsi_h_u + zpsi_h_z0t

               !IF( l_use_cs ) THEN
               !   !! Cool-skin contribution
               !   CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, ztzu, zqzu, zus, zts, zqs, &
               !      &                   U_zu(ji,jj), zubzu, slp(ji,jj), rad_lw(ji,jj), zQns, ztmp0 )  ! Tau -> ztmp0
               !
               !   CALL CS_ECMWF( Qsw(ji,jj), zQns, zus, zSST, zdT_cs )
               !   !IF( PRESENT(pdT_cs) ) pdT_cs(ji,jj) = zdT_cs
               !   zT_s = zSST + zdT_cs
               !   IF( l_use_wl ) zT_s = zT_s + dT_wl(ji,jj)
               !   zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), slp(ji,jj))
               !ENDIF

               !IF( l_use_wl ) THEN
               !   !! Warm-layer contribution
               !   CALL UPDATE_QNSOL_TAU( zu, zT_s, zq_s, ztzu, zqzu, zus, zts, zqs, U_zu(ji,jj), zubzu, &
               !      &                   slp(ji,jj), rad_lw(ji,jj), zQns, ztmp0)  ! Tau -> ztmp0
               !   !IF((ji==10).AND.(jj==10)) PRINT *, 'LOLO: sbcblk_algo_ecmwf.F90 => depth SST = ', 0.5*e3t_m(ji,jj)
               !   !CALL WL_ECMWF( ji, jj, Qsw(ji,jj), zQns, zus, zsst(ji,jj), 0.5*e3t_m(ji,jj) )
               !   CALL WL_ECMWF( ji, jj, Qsw(ji,jj), zQns, zus, zSST, 0.5*e3t_m(ji,jj) )
               !   !! Updating T_s and q_s !!!
               !   zT_s = zSST + dT_wl(ji,jj)
               !   IF( l_use_cs ) zT_s = zT_s + zdT_cs
               !   zq_s = rdct_qsat_salt*q_sat(MAX(zT_s, 200._wp), slp(ji,jj))
               !ENDIF

               !IF( l_use_cs .OR. l_use_wl .OR. (.NOT. l_zt_equal_zu) ) THEN
               zdt = ztzu - zT_s ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
               zdq = zqzu - zq_s ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
               !ENDIF

            END DO !DO jit = 1, nb_iter

            ! compute transfer coefficients at zu :
            zFq = zlog_zu - zlog_z0q - zpsi_h_u + zpsi_h_z0q

            Cd(ji,jj)   = MAX( vkarmn2/(zFm*zFm) , Cx_min )
            Ch(ji,jj)   = MAX( vkarmn2/(zFm*zFh) , Cx_min )
            Ce(ji,jj)   = MAX( vkarmn2/(zFm*zFq) , Cx_min )
            !T_s(ji,jj) = zT_s
            !q_s(ji,jj) = zq_s
            t_zu(ji,jj) = ztzu
            q_zu(ji,jj) = zqzu
            Ubzu(ji,jj) = zubzu

         END DO !DO ji = Nis0, Nie0
      END DO !DO jj = Njs0, Nje0
      !$acc end parallel loop

      !$acc end data

      !IF( kt == nitend ) CALL sbcblk_cs_wl_dealloc( l_use_wl )

      IF( ln_timing )   CALL timing_stop('turb_ecmwf')

   END SUBROUTINE turb_ecmwf


   !!===============================================================================================
   FUNCTION psi_m_ecmwf_sclr( pzeta )
      !!--------------------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!--------------------------------------------------------------------------------------------
      !$acc routine
      REAL(wp), INTENT(in) :: pzeta
      REAL(wp)             :: psi_m_ecmwf_sclr
      !!
      REAL(wp) :: zta, zx2, zx, ztmp, zpsi_unst, zpsi_stab, zstab, zc
      !!--------------------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp

      zta = pzeta
      CALL cap_zeta( zta )

      ! *** Unstable (Paulson 1970)    [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
      zx  = SQRT(zx2)                        ! (1 - 16z)^0.25
      ztmp = 1._wp + zx
      zpsi_unst = LOG( 0.125_wp*ztmp*ztmp*(1._wp + zx2) ) - 2._wp*ATAN( zx ) + 0.5_wp*rpi

      ! *** Stable                   [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
      zpsi_stab = -2._wp/3._wp*(zta - zc)*EXP(-0.35_wp*zta) &
         &        - zta - 2._wp/3._wp*zc
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_m_ecmwf_sclr =         zstab    * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
      !
   END FUNCTION psi_m_ecmwf_sclr
   !!===============================================================================================


   !!===============================================================================================
   FUNCTION psi_h_ecmwf_sclr( pzeta )
      !!--------------------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!     ECMWF / as in IFS cy31r1 documentation, available online
      !!     at ecmwf.int
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!--------------------------------------------------------------------------------------------
      !$acc routine
      REAL(wp), INTENT(in) :: pzeta
      REAL(wp)             :: psi_h_ecmwf_sclr
      !!
      REAL(wp) ::  zta, zx2, zpsi_unst, zpsi_stab, zstab, zc
      !!--------------------------------------------------------------------------------------------
      zc = 5._wp/0.35_wp

      zta = pzeta
      CALL cap_zeta( zta )

      ! *** Unstable (Paulson 1970)   [eq.3.20, Chap.3, p.33, IFS doc - Cy31r1] :
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
      zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
      !
      ! *** Stable [eq.3.22, Chap.3, p.33, IFS doc - Cy31r1] :
      zpsi_stab = -2._wp/3._wp*(zta - zc)*EXP(-0.35_wp*zta) &
         &       - ABS(1._wp + 2._wp/3._wp*zta)**1.5_wp - 2._wp/3._wp*zc + 1._wp
      !! LB: added ABS() to avoid NaN values when unstable, which contaminates the unstable solution...
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_h_ecmwf_sclr =        zstab     * zpsi_stab   &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst      ! (zta < 0) Unstable
      !
   END FUNCTION psi_h_ecmwf_sclr
   !!===============================================================================================

   SUBROUTINE cap_zeta( pzeta )
      !!--------------------------------------------------------------------------------------------
      !$acc routine
      REAL(wp), INTENT(inout) :: pzeta
      REAL(wp) ::  zta
      !!--------------------------------------------------------------------------------------------
      !#fixme: Jean Bildot & Sam Hatfield @ ECMWF, complain that
      !        `EXP(-0.35_wp*zta)` later blows up in single precision when unstable with big `zta`
      !#fixme: LB suggests:
      zta = MAX( pzeta , -50._wp ) ! => regions where `zeta<-50.` are given value -50 (still unrealistic but numerically safe?)
      !                            !  ==> prevents numerical problems such as overflows...
      zta = MIN(  zta ,   5._wp )  !`zeta` plateaus at 5 in very stable conditions (L>0 and small!), inherent to ECMWF algo!
      !
      pzeta = zta
   END SUBROUTINE cap_zeta
   !!======================================================================

END MODULE sbcblk_algo_ecmwf
