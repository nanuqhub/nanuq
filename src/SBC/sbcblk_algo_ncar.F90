MODULE sbcblk_algo_ncar
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ncar  ***
   !!
   !!       Computes turbulent components of surface fluxes
   !!         according to Large & Yeager (2004,2008)
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m Ubzu
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
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
   !!   turb_ncar  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE in_out_manager, ONLY: ln_timing
   USE phycst          ! physical constants
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR   ! called by sbcblk.F90

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu,   &
      &                        Cd, Ch, Ce, t_zu, q_zu, Ubzu,   &
      &                        nb_iter )
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ncar  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  sst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!
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
      !! OPTIONAL INPUT:
      !! ---------------
      !!    * nb_iter : nb. of itterations (default is `nb_iter0`)
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                  [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                           [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature            [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity         [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu              [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat    [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation     [-]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu             [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu         [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ubzu     ! bulk wind speed at zu                   [m/s]
      INTEGER , INTENT(in   ), OPTIONAL           :: nb_iter    ! number of iterations
      !!----------------------------------------------------------------------------------
      INTEGER  :: ji, jj, nbit, jit
      REAL(wp) :: zm_ztzu                   ! => `1.` if `zu /= zt`, `0.` otherwize
      REAL(wp) :: zstab, zCdN, zCeN, zChN   ! 10m neutral latent/sensible coefficient
      REAL(wp) :: zsqrt_Cd, zsqrt_CdN       ! square root of Cd_n10
      REAL(wp) :: zeta_u, zeta_t            ! stability parameter at height zu and zt
      REAL(wp) :: zlog1, zlog2, ztmp, ztmp2
      REAL(wp) :: zdt, zdq, zus, zts, zqs, z1oL, zpsi_m, zUn10
      !!----------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('turb_ncar')

      nbit = nb_iter0
      IF( PRESENT(nb_iter) ) nbit = nb_iter

      zm_ztzu = 0.5_wp + SIGN( 0.5_wp, ABS(zu-zt) - 0.01_wp ) ! => `1.` if `zu /= zt`, `0.` otherwize

      Ubzu = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! ij-independant constants:
      zlog1 = LOG(zt/zu)
      zlog2 = LOG(zu/10._wp)

      DO jj = Njs0-1, Nje0+1
         DO ji = Nis0-1, Nie0+1

            !! First guess of stability:
            zstab = 0.5_wp + SIGN( 0.5_wp , virt_temp(t_zt(ji,jj), q_zt(ji,jj)) - virt_temp(sst(ji,jj), ssq(ji,jj)) )

            !! Neutral coefficients at 10m:
            zCdN = cd_n10_ncar( Ubzu(ji,jj) )
            zsqrt_CdN = SQRT( zCdN )

            !! Initializing transf. coeff. with their first guess neutral equivalents :
            Cd(ji,jj) = zCdN
            Ce(ji,jj) = ce_n10_ncar( zsqrt_CdN )
            Ch(ji,jj) = ch_n10_ncar( zsqrt_CdN , zstab )   ! zstab is stability (1/0)
            zsqrt_Cd = zsqrt_CdN

            !! Initializing values at z_u with z_t values:
            t_zu(ji,jj) = MAX( t_zt(ji,jj) ,  180._wp )   ! who knows what's given on masked-continental regions...
            q_zu(ji,jj) = MAX( q_zt(ji,jj) , 1.e-6_wp )   !               "

            !! ITERATION BLOCK
            DO jit = 1, nbit
               !
               zdt = t_zu(ji,jj) - sst(ji,jj)   ! Updating air/sea differences
               zdq = q_zu(ji,jj) - ssq(ji,jj)

               ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
               zus = zsqrt_Cd*Ubzu(ji,jj)      ! u*
               zts = Ch(ji,jj)/zsqrt_Cd*zdt    ! theta*
               zqs = Ce(ji,jj)/zsqrt_Cd*zdq    ! q*

               ! Estimate the inverse of Obukov length (1/L) at height zu:
               z1oL = One_on_L( t_zu(ji,jj), q_zu(ji,jj), zus, zts, zqs )

               !! Stability parameters :
               zeta_u   = zu*z1oL
               zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

               !! Shifting temperature and humidity at zu if required by `zm_ztzu` (L&Y 2004 Eq.9b-9c)
               zeta_t = zt*z1oL ! zeta_t !
               zeta_t = SIGN( MIN(ABS(zeta_t),10._wp), zeta_t )
               ztmp = zlog1 + psi_h_ncar(zeta_u) - psi_h_ncar(zeta_t)
               t_zu(ji,jj) =       t_zt(ji,jj) - zm_ztzu*zts/vkarmn*ztmp
               q_zu(ji,jj) = MAX(  q_zt(ji,jj) - zm_ztzu*zqs/vkarmn*ztmp , 0._wp )

               ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 Eq. 9a)...
               !   In very rare low-wind conditions, the old way of estimating the
               !   neutral wind speed at 10m leads to a negative value that causes the code
               !   to crash. To prevent this a threshold of 0.25m/s is imposed.
               zpsi_m = psi_m_ncar(zeta_u)
               zUn10 = MAX( 0.25_wp , UN10_from_CD(zu, Ubzu(ji,jj), Cd(ji,jj), ppsi=zpsi_m) )
               zCdN = cd_n10_ncar(zUn10)
               zsqrt_CdN = SQRT(zCdN)

               !! Update of transfer coefficients:
               !! C_D
               ztmp  = 1._wp + zsqrt_CdN/vkarmn*(zlog2 - zpsi_m)   ! L&Y 2004 Eq. (10a) (zpsi_m == psi_m(zeta_u))
               Cd(ji,jj)     = MAX( zCdN / ( ztmp*ztmp ), Cx_min )
               !! C_H and C_E
               zsqrt_Cd = SQRT( Cd(ji,jj) )
               ztmp = ( zlog2 - psi_h_ncar(zeta_u) ) / vkarmn / zsqrt_CdN
               ztmp2 = zsqrt_Cd / zsqrt_CdN

               zstab = 0.5_wp + SIGN(0.5_wp,zeta_u)                                ! update stability
               zChN  = 1.e-3_wp * zsqrt_CdN*(18._wp*zstab + 32.7_wp*(1._wp - zstab))  ! L&Y 2004 eq. (6c-6d)
               zCeN  = 1.e-3_wp * (34.6_wp * zsqrt_CdN)                             ! L&Y 2004 eq. (6b)

               Ch(ji,jj)    = MAX( zChN*ztmp2 / ( 1._wp + zChN*ztmp ) , Cx_min ) ! L&Y 2004 eq. (10b)
               Ce(ji,jj)    = MAX( zCeN*ztmp2 / ( 1._wp + zCeN*ztmp ) , Cx_min ) ! L&Y 2004 eq. (10c)

            END DO !DO jit = 1, nb_iter

         END DO
      END DO

      IF( ln_timing )   CALL timing_stop('turb_ncar')

   END SUBROUTINE turb_ncar


   !!===============================================================================================
   FUNCTION cd_n10_ncar( pw10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008, Eq. (11)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp)             :: cd_n10_ncar
      !!
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      zw  = pw10
      zw6 = zw*zw*zw
      zw6 = zw6*zw6
      !
      ! When wind speed > 33 m/s => Cyclone conditions => special treatment
      zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
      !
      cd_n10_ncar = 1.e-3_wp * ( &
         &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
         &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
      !
      cd_n10_ncar = MAX( cd_n10_ncar, Cx_min )
      !
   END FUNCTION cd_n10_ncar
   !!===============================================================================================

   !!===============================================================================================
   FUNCTION ch_n10_ncar( psqrtcdn10 , pstab )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (12)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp), INTENT(in) :: pstab      ! stable ABL => 1 / unstable ABL => 0
      REAL(wp)             :: ch_n10_ncar
      !!----------------------------------------------------------------------------------
      IF( (pstab < -0.00001).OR.(pstab >  1.00001) ) THEN
         PRINT *, 'ERROR: ch_n10_ncar@mod_blk_ncar.f90: pstab ='
         PRINT *, pstab
         STOP
      END IF
      ch_n10_ncar = MAX( 1.e-3_wp * psqrtcdn10*( 18._wp*pstab + 32.7_wp*(1._wp - pstab) )  , Cx_min )   ! Eq. (9) & (12) Large & Yeager, 2008
   END FUNCTION ch_n10_ncar
   !!===============================================================================================

   !!===============================================================================================
   FUNCTION ce_n10_ncar( psqrtcdn10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (13)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp)             :: ce_n10_ncar
      !!----------------------------------------------------------------------------------
      ce_n10_ncar = MAX( 1.e-3_wp * ( 34.6_wp * psqrtcdn10 ) , Cx_min )
   END FUNCTION ce_n10_ncar
   !!===============================================================================================


   !!===============================================================================================
   FUNCTION psi_m_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pzeta
      REAL(wp)             :: psi_m_ncar
      !!
      REAL(wp) :: zta, zx2, zx, zpsi_unst, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      zta = pzeta
      !
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 - 16z)^0.5
      zx2 = MAX( zx2 , 1._wp )
      zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
      zpsi_unst = 2._wp*LOG( (1._wp + zx )*0.5_wp )   &
         &            + LOG( (1._wp + zx2)*0.5_wp )   &
         &          - 2._wp*ATAN(zx) + rpi*0.5_wp
      !
      zpsi_stab = -5._wp*zta
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_m_ncar =          zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
   END FUNCTION psi_m_ncar
   !!===============================================================================================

   !!===============================================================================================
   FUNCTION psi_h_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pzeta
      REAL(wp)             :: psi_h_ncar
      !!
      REAL(wp) :: zta, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      zta = pzeta
      !
      zx2 = SQRT( ABS(1._wp - 16._wp*zta) )  ! (1 -16z)^0.5
      zx2 = MAX( zx2 , 1._wp )
      zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
      !
      zpsi_stab = -5._wp*zta
      !
      zstab = 0.5_wp + SIGN(0.5_wp, zta) ! zta > 0 => zstab = 1
      !
      psi_h_ncar =          zstab  * zpsi_stab &  ! (zta > 0) Stable
         &              + (1._wp - zstab) * zpsi_unst    ! (zta < 0) Unstable
      !
   END FUNCTION psi_h_ncar
   !!===============================================================================================

   !!======================================================================
END MODULE sbcblk_algo_ncar
