MODULE icedyn_rhg_tools
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_tools  ***
   !!   Sea-Ice dynamics : misc. rheology tools...
   !!======================================================================
   !! History : L. Brodeau, 2024
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE dom_oce        ! Ocean domain
   USE par_ice, ONLY: rn_delta_ecc, rn_creepl
   USE lib_mpp,        ONLY: ctl_stop
   !USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE in_out_manager, ONLY: ln_timing, lwp
   USE timing

   IMPLICIT NONE

   PRIVATE

   !INTERFACE strain_rate
   !   MODULE PROCEDURE strain_rate_Cgrid, strain_rate_Egrid
   !END INTERFACE strain_rate

   PUBLIC sigmaII_sclr
   PUBLIC sigmaII_full
   PUBLIC update_invariants_Egrid

   PUBLIC strain_rate_all
   PUBLIC strain_rate_dsd
   PUBLIC strain_rate_min

   PUBLIC low_conc_canceler
   PUBLIC cancel_low_conc
   PUBLIC fdamp_low_conc

   PUBLIC div_stress_tensor

   PUBLIC vel_div_t
   PUBLIC vel_ten_t
   PUBLIC vel_shear_f
   PUBLIC vel_maxshr_t
   !PUBLIC vel_delta_t



   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
CONTAINS


   FUNCTION sigmaII_sclr( ps11, ps22, ps12 )
      !!------------------------------------------------------------------------------------
      !$acc routine
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor => same units as input stresses
      !!------------------------------------------------------------------------------------
      REAL(wp)             :: sigmaII_sclr
      REAL(wp), INTENT(in) :: ps11, ps22, ps12  ! sigma_11, sigma_22, sigma_12
      !!
      REAL(wp) :: ztmp
      !!------------------------------------------------------------------------------------
      ztmp  = 0.5_wp * (ps11 - ps22)
      sigmaII_sclr = SQRT( ztmp*ztmp + ps12*ps12 )
      !!
   END FUNCTION sigmaII_sclr

   SUBROUTINE sigmaII_full( pSt, pSf, pSII )
      !!------------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor => same units as input stresses
      !!------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(in)  :: pSt, pSf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(out) :: pSII
      !!
      REAL(wp) :: ztmp, zs
      INTEGER  :: ji, jj
      !!------------------------------------------------------------------------------------
      !$acc data present( pSt, pSf, pSII )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ztmp        = 0.5_wp * ( pSt(ji,jj,1) - pSt(ji,jj,2) )
            zs          = pSf(ji,jj,3)
            pSII(ji,jj) = SQRT( ztmp*ztmp + zs*zs )
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE sigmaII_full


   SUBROUTINE update_invariants_Egrid( pSt, pSf,  pSIt, pSIIt, pSIf, pSIIf )
      !!------------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor => same units as input stresses
      !!------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(in)  :: pSt, pSf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(out) :: pSIt, pSIIt, pSIf, pSIIf
      !!------------------------------------------------------------------------------------
      REAL(wp) :: ztmp, zs11, zs22, zs12
      INTEGER  :: ji, jj
      !!------------------------------------------------------------------------------------
      !%acc data present( pSt, pSf, pSI, pSII )
      !%acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            ! T-points:
            zs11 = pSt(ji,jj,1)
            zs22 = pSt(ji,jj,2)
            zs12 = pSf(ji,jj,3)
            !
            pSIt(ji,jj)  = 0.5_wp * ( zs11 + zs22 )
            ztmp         = 0.5_wp * ( zs11 - zs22 )
            pSIIt(ji,jj) = SQRT( ztmp*ztmp + zs12*zs12 )
            !
            ! F-points:
            zs11 = pSf(ji,jj,1)
            zs22 = pSf(ji,jj,2)
            zs12 = pSt(ji,jj,3)
            !
            pSIf(ji,jj)  = 0.5_wp * ( zs11 + zs22 )
            ztmp         = 0.5_wp * ( zs11 - zs22 )
            pSIIf(ji,jj) = SQRT( ztmp*ztmp + zs12*zs12 )
            !
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
   END SUBROUTINE update_invariants_Egrid





   SUBROUTINE strain_rate_all( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                               pe11, pe22, pe12, pdudy, pdvdx, pdiv, pmshr, pdelta )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the 3 elements of the strain rate tensor, e11, e22 & e12, at either T- or F-points
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx` and `pmshr` !)
      !!----------------------------------------------------------------------------------------------------------
      CHARACTER(len=1),                   INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pe11, pe22, pe12          ! e11, e22 & e12 @ `cgt` points                   [1/s]
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pdudy, pdvdx, pdiv, pmshr  ! @ `cgt` points                [1/s]
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pdelta
      !!----------------------------------------------------------------------------------------------------------
      LOGICAL  :: l_r_e11, l_r_e22, l_r_e12, l_r_dudy, l_r_dvdx, l_r_div, l_r_mshr, l_r_dlt
      REAL(wp) :: zE1, zE2, zS1, zS2, zSHR, zdlt, z1_e1e2, zzf, ze2e2, ze1e1, zmask, zswitch, z1_ecc2
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!----------------------------------------------------------------------------------------------------------
      !!
      l_r_e11  = PRESENT(  pe11  )
      l_r_e22  = PRESENT(  pe22  )
      l_r_e12  = PRESENT(  pe12  )
      l_r_dudy = PRESENT(  pdudy )
      l_r_dvdx = PRESENT(  pdvdx )
      l_r_div  = PRESENT(  pdiv  )
      l_r_mshr = PRESENT(  pmshr )
      l_r_dlt  = PRESENT( pdelta )

      IF( l_r_dlt ) THEN
#if defined key_verbose
         IF(lwp) PRINT *, 'LOLO [strain_rate_all@icedyn_rhg_tools.F90]: for `delta` => using ecc =', REAL(rn_delta_ecc)
#endif
         z1_ecc2 = 1._wp / ( rn_delta_ecc*rn_delta_ecc )
      ENDIF

      ! Prevent the occurence of NaN on the halos:
      IF( l_r_e11  )  pe11(:,:)   = 0._wp
      IF( l_r_e22  )  pe22(:,:)   = 0._wp
      IF( l_r_e12  )  pe12(:,:)   = 0._wp
      IF( l_r_dudy )  pdudy(:,:)  = 0._wp
      IF( l_r_dvdx )  pdvdx(:,:)  = 0._wp
      IF( l_r_div  )  pdiv(:,:)   = 0._wp
      IF( l_r_mshr )  pmshr(:,:)  = 0._wp
      IF( l_r_dlt  )  pdelta(:,:) = 0._wp

      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_all(): unknown grid-point type: '//cgt//'!')
      ENDIF

      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            IF( l_r_div .OR. l_r_e11 .OR. l_r_e22 .OR. l_r_mshr .OR. l_r_dlt ) THEN
               !! Divergence at cgt-points, `dU/dx + dV/dy` :
               zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
                  &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
                  &  ) * z1_e1e2
            ENDIF
            IF( l_r_e11 .OR. l_r_e22 .OR. l_r_mshr .OR. l_r_dlt  ) THEN
               !! Tension at cgt-points, `dU/dx - dV/dy` :
               zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
                  &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
                  &  ) * z1_e1e2
            ENDIF
            IF( l_r_e11 ) pe11(ji,jj) = 0.5_wp * ( zE1 + zE2 )
            IF( l_r_e22 ) pe22(ji,jj) = 0.5_wp * ( zE1 - zE2 )
            IF( l_r_div )  pdiv(ji,jj)  = zE1

            IF( l_r_e12 .OR. l_r_dudy .OR. l_r_dvdx .OR. l_r_mshr .OR. l_r_dlt ) THEN
               zzf = z1_e1e2 * zmask
               zS1 = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf    ! du/dy
               zS2 = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf    ! dv/dx
            ENDIF

            IF( l_r_e12 .OR. l_r_mshr .OR. l_r_dlt ) zSHR =  zS1 + zS2   ! shearing strain rate == 2*eps12 !

            IF( l_r_e12 ) pe12(ji,jj) = 0.5_wp * zSHR      ! pe12 == eps12 = 1/2 `shearing strain rate` !

            IF( l_r_dudy )  pdudy(ji,jj) = zS1
            IF( l_r_dvdx )  pdvdx(ji,jj) = zS2
            IF( l_r_mshr )  pmshr(ji,jj) = SQRT( zE2*zE2 + zSHR*zSHR )  ! Maximum shear: == SQRT( T^2 + S^2 )

            IF( l_r_dlt ) THEN
               zdlt = SQRT( zE1*zE1 + ( zE2*zE2 + zSHR*zSHR ) * z1_ecc2 ) ! `Delta` with z1_ecc2 => 0.25   ! Hunke & Dukowicz, 2002, Eq.5
               zswitch       = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zdlt ) ) ! 0 if delta=0
               pdelta(ji,jj) = zdlt + rn_creepl * zswitch
            ENDIF

         END DO
      END DO

   END SUBROUTINE strain_rate_all





   SUBROUTINE strain_rate_dsd( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                                 pdiv, pmshr, pdelta )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the divergence, maximum shear & delta
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx` and `pmshr` !)
      !!----------------------------------------------------------------------------------------------------------
      CHARACTER(len=1),             INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdiv, pmshr, pdelta  ! @ `cgt` points                [1/s]
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp) :: zE1, zE2, zS1, zS2, zSHR, zdlt, z1_e1e2, zzf, ze2e2, ze1e1, zmask, zsw, z1_ecc2
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!----------------------------------------------------------------------------------------------------------
      !$acc data present( pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, pdiv, pmshr, pdelta )
      z1_ecc2 = 1._wp / ( rn_delta_ecc*rn_delta_ecc )

      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_dsd(): unknown grid-point type: '//cgt//'!')
      ENDIF

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pdiv(ji,jj)   = 0._wp
            pmshr(ji,jj)  = 0._wp
            pdelta(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            !! Divergence at cgt-points, `dU/dx + dV/dy` :
            zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
               &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
               &  ) * z1_e1e2

            !! Tension at cgt-points, `dU/dx - dV/dy` :
            zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
               &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
               &  ) * z1_e1e2

            pdiv(ji,jj)  = zE1

            zzf = z1_e1e2 * zmask
            zS1 = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf    ! du/dy
            zS2 = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf    ! dv/dx
            zSHR =  zS1 + zS2   ! shearing strain rate == 2*eps12 !
            pmshr(ji,jj) = SQRT( zE2*zE2 + zSHR*zSHR )  ! Maximum shear: == SQRT( T^2 + S^2 )

            zdlt = SQRT( zE1*zE1 + ( zE2*zE2 + zSHR*zSHR ) * z1_ecc2 ) ! `Delta` with z1_ecc2 => 0.25   ! Hunke & Dukowicz, 2002, Eq.5
            zsw       = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zdlt ) ) ! 0 if delta=0
            pdelta(ji,jj) = zdlt + rn_creepl * zsw

         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE strain_rate_dsd




   SUBROUTINE strain_rate_min( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                               pe11, pe22, pdiv, pdudy, pdvdx )
      !!
      !! Computes the 3 elements of the strain rate tensor, e11, e22 & e12, at either T- or F-points
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx`!)
      !!
      CHARACTER(len=1),         INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pe11, pe22, pdiv
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdudy, pdvdx
      !!
      REAL(wp) :: zE1, zE2, z1_e1e2, zzf, ze2e2, ze1e1, zmask
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!
      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_min(): unknown grid-point type: '//cgt//'!')
      ENDIF

      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            !! Divergence at cgt-points, `dU/dx + dV/dy` :
            zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
               &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
               &  ) * z1_e1e2
            !! Tension at cgt-points, `dU/dx - dV/dy` :
            zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
               &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
               &  ) * z1_e1e2
            !!
            pe11(ji,jj) = 0.5_wp * ( zE1 + zE2 )
            pe22(ji,jj) = 0.5_wp * ( zE1 - zE2 )
            pdiv(ji,jj)  = zE1

            !! 2 * shear at cgt-points, `dU/dy + dV/dx` :
            zzf = z1_e1e2 * zmask
            pdudy(ji,jj) = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf
            pdvdx(ji,jj) = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf
            !
         END DO
      END DO

   END SUBROUTINE strain_rate_min





   SUBROUTINE low_conc_canceler( pA, pcncl )
      !!------------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------------
      !! Create an array intended to be used (trough multiplication) to gradually cancel
      !! a given array field at low ice concentration.
      !!
      !!
      !! It is mainly used to cancel the components of the divergence of the `h*SIGMA` tensors
      !! as spatial derivatives of `h*SIGMA` tend to become a nonsense at low ice resolution
      !!
      !! Here is the "gnuplot-read" equation of the function we use
      !!  ``` plot 0.51 * ( 1. + 20.*(x-0.3) / sqrt( 1 + (20.*(x-0.3))**2 ) ) - 0.008 ```
      !!
      !! `x` being the ice concentration
      !!
      !! => looks like a smooth step function that is 0 at `x=0` and reaches 1 at about `x=0.5`
      !!
      !!------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pA    ! ice concentration at point "X" [0:1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pcncl ! correction factor [0:1]
      !!------------------------------------------------------------------------------------
      REAL(wp) :: zx, zc
      INTEGER  :: ji, jj
      !!------------------------------------------------------------------------------------
      !$acc data present( pA, pcncl )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            zx = 20._wp * (pA(ji,jj) - 0.3_wp)
            zc = 0.51_wp * ( 1._wp + zx / SQRT(1._wp + zx*zx) ) - 0.008_wp
            pcncl(ji,jj) = MIN( MAX( zc , 0._wp ) , 1._wp )
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE low_conc_canceler


   SUBROUTINE cancel_low_conc(  pA, pF )
      !!------------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------------
      !! Create an array intended to be used (trough multiplication) to gradually cancel
      !! a given array field at low ice concentration.
      !!
      !!
      !! It is mainly used to cancel the components of the divergence of the `h*SIGMA` tensors
      !! as spatial derivatives of `h*SIGMA` tend to become a nonsense at low ice resolution
      !!
      !! Here is the "gnuplot-read" equation of the function we use
      !!  ```plot 0.505 * ( 1. + 25.*(x-0.5) / sqrt( 1 + (25.*(x-0.5))**2 ) ) - 0.003`
      !!
      !! `x` being the ice concentration
      !!
      !! => looks like a smooth step function that is 0 at `x=0` and reaches 1 at about `x=0.7`
      !!
      !!------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pA    ! ice concentration at point "X" [0:1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pF    ! field to correct
      !!------------------------------------------------------------------------------------
      REAL(wp) :: zx, zc, zF
      INTEGER  :: ji, jj
      !!------------------------------------------------------------------------------------
      !$acc data present( pA, pF )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            zF = pF(ji,jj)
            zx = 25._wp * (pA(ji,jj) - 0.5_wp)
            zc = 0.505_wp * ( 1._wp + zx / SQRT(1._wp + zx*zx) ) - 0.003_wp
            zF = MIN( MAX( zc , 0._wp ) , 1._wp ) * zF
            pF(ji,jj) = zF
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE cancel_low_conc


   FUNCTION fdamp_low_conc(  pA )
      !!------------------------------------------------------------------------------------
      !$acc routine seq
      !!------------------------------------------------------------------------------------
      !! Create an array intended to be used (trough multiplication) to gradually cancel
      !! a given array field at low ice concentration.
      !!
      !!
      !! It is mainly used to cancel the components of the divergence of the `h*SIGMA` tensors
      !! as spatial derivatives of `h*SIGMA` tend to become a nonsense at low ice resolution
      !!
      !! Here is the "gnuplot-read" equation of the function we use
      !! OLD:  ```plot 0.505 * ( 1. + 25.*(x-0.5 ) / sqrt( 1 + (25.*(x-0.5 ))**2 ) ) - 0.003`
      !! NEW:  ```plot 0.505 * ( 1. + 40.*(x-0.15) / sqrt( 1 + (40.*(x-0.15))**2 ) ) - 0.005`
      !!
      !! `x` being the ice concentration
      !!
      !! => looks like a smooth step function that is 0 at `x=0` and reaches 1 at about `x=0.7`
      !!
      !!------------------------------------------------------------------------------------
      REAL(wp)             :: fdamp_low_conc
      REAL(wp), INTENT(in) :: pA    ! ice concentration at point "X" [0:1]
      !!------------------------------------------------------------------------------------
      REAL(wp) :: zx, zc
      !!------------------------------------------------------------------------------------
      zx = 40._wp * (pA - 0.15_wp)
      zc = 0.505_wp * ( 1._wp + zx / SQRT(1._wp + zx*zx) ) - 0.005_wp
      fdamp_low_conc = MIN( MAX( zc , 0._wp ) , 1._wp )
      !!
   END FUNCTION fdamp_low_conc







   SUBROUTINE div_stress_tensor( cgt, phc, phx, pe1e1, pe2e2,  pe1e1_e, pe2e2_e,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y,  &
      &                               ps11c, ps22c, ps12x,  pdSx, pdSy,  pm0 )
      !!----------------------------------------------------------------------------------------------
      !! Computes the vector (pdSx,pdSy) = divergence of the h-integrated internal stress tensor
      !!
      !!   depending on the grid: T-centric grid => cgt='T' or F-centric grid => cgt='F'
      !!
      !! INPUT:                                               |     cgt=='T'   |    cgt=='F'    |
      !!   * ps11c, ps22c: sigma11, sigma22           =>  ! @ point T[i,j] | @ point F[i,j] |
      !!   * ps12x       :       sigma12                =>  ! @ point F[i,j] | @ point T[i,j] |
      !!
      !! RETURNS:                                             |     cgt=='T'   |    cgt=='F'    |
      !!   * pdSx: x-component of the div of the tensor =>  | @ point U[i,j] | @ point V[i,j] |
      !!   * pdSy: y-component of the div of the tensor =>  | @ point V[i,j] | @ point U[i,j] |
      !!
      !!----------------------------------------------------------------------------------------------
      CHARACTER(len=1),             INTENT(in)  :: cgt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: phc, phx   ! ice thickness at center and corner point [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pe1e1, pe2e2, pe1e1_e, pe2e2_e
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e1e2x, p1_e1e2y
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: ps11c, ps22c, ps12x ! components of stress tensors on T- or F-centric grids x h !!!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdSx, pdSy      ! x,y components of the divergence of the tensor
      INTEGER,        OPTIONAL,     INTENT(in)  :: pm0
      !!
      INTEGER  :: ip, im, jp, jm, ji, jj, m0
      !!--------------------------------------------------------------------------------------------
      IF( ln_timing ) CALL timing_start('div_stress_tensor')
      m0 = 0
      IF( PRESENT(pm0) ) m0 = pm0

      IF ( cgt == 'T' ) THEN
         ip =  1
         im =  0
         jp =  0
         jm = -1
      ELSEIF ( cgt == 'F' ) THEN
         ip =  0
         im = -1
         jp =  1
         jm =  0
      ELSE
         CALL ctl_stop( 'STOP', 'div_stress_tensor(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      !pdSx(:,:) = 0._wp
      !pdSy(:,:) = 0._wp
      !
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !                   !--- ds11/dx + ds12/dy
            pdSx(ji,jj) = ( ( ps11c(ji+ip,jj)*phc(ji+ip,jj)*pe2e2(ji+ip,jj)   - ps11c(ji+im,jj)*phc(ji+im,jj)*pe2e2(ji+im,jj)   ) * p1_e2x(ji,jj) &
               &          + ( ps12x(ji,jj+jp)*phx(ji,jj+jp)*pe1e1_e(ji,jj+jp) - ps12x(ji,jj+jm)*phx(ji,jj+jm)*pe1e1_e(ji,jj+jm) ) * p1_e1x(ji,jj) &
               &                 ) * p1_e1e2x(ji,jj)
            !                   !--- ds22/dy + ds12/dx
            pdSy(ji,jj) = ( ( ps22c(ji,jj-jm)*phc(ji,jj-jm)*pe1e1(ji,jj-jm)   - ps22c(ji,jj-jp)*phc(ji,jj-jp)*pe1e1(ji,jj-jp)   ) * p1_e1y(ji,jj) &
               &          + ( ps12x(ji-im,jj)*phx(ji-im,jj)*pe2e2_e(ji-im,jj) - ps12x(ji-ip,jj)*phx(ji-ip,jj)*pe2e2_e(ji-ip,jj) ) * p1_e2y(ji,jj) &
               &               ) * p1_e1e2y(ji,jj)
            !
         END DO
      END DO
      !
      IF( ln_timing )   CALL timing_stop('div_stress_tensor')
      !
   END SUBROUTINE div_stress_tensor



   SUBROUTINE vel_div_t( pU, pV, p1_e1e2t, pe2u, pe1v, pmskt, pdivt )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) divergence of sea-ice velocity vector => T-point
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t         ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2u, pe1v       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt           ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdivt   ! divergence @ T points               [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Divergence at T-points, `dU/dx + dV/dy` :
            pdivt(ji,jj) = (   pe2u(ji,jj)*pU(ji,jj) - pe2u(ji-1,jj)*pU(ji-1,jj) &
               &             + pe1v(ji,jj)*pV(ji,jj) - pe1v(ji,jj-1)*pV(ji,jj-1) &
               &            )           * p1_e1e2t(ji,jj) * pmskt(ji,jj)

         END DO
      END DO

   END SUBROUTINE vel_div_t


   SUBROUTINE vel_ten_t( pU, pV, p1_e1e2t, p1_e2u, p1_e1v, pe1e1t, pe2e2t, pmskt, ptent )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) tension of sea-ice velocity vector => T-point
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t         ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2u, p1_e1v   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1t, pe2e2t
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt           ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: ptent   ! tension @ T points               [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Tension at T-points, `dU/dx - dV/dy` :
            ptent(ji,jj) = (   ( pU(ji,jj)*p1_e2u(ji,jj) - pU(ji-1,jj)*p1_e2u(ji-1,jj) ) * pe2e2t(ji,jj) &
               &             - ( pV(ji,jj)*p1_e1v(ji,jj) - pV(ji,jj-1)*p1_e1v(ji,jj-1) ) * pe1e1t(ji,jj) &
               &            )           * p1_e1e2t(ji,jj) * pmskt(ji,jj)

         END DO
      END DO

   END SUBROUTINE vel_ten_t


   SUBROUTINE vel_shear_f( pU, pV, p1_e1e2f, p1_e1u, p1_e2v, pe1e1f, pe2e2f, pmskf, pe12f )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) shear of sea-ice velocity vector => F-point
      !!
      !! Note: the mask must be the actual `fmask` that takes into account the slip/no-slip conditions
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV            ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2f         !     1/(e1f*e2f)               [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1u, p1_e2v
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1f, pe2e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskf         ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pe12f  !  e12 @ F points            [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      REAL(wp) :: zzf, zS1, zS2
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            ! Shear at F points:
            zzf = p1_e1e2f(ji,jj) * pmskf(ji,jj)
            zS1 = ( pU(ji,jj+1) * p1_e1u(ji,jj+1) - pU(ji,jj) * p1_e1u(ji,jj) ) * pe1e1f(ji,jj) * zzf
            zS2 = ( pV(ji+1,jj) * p1_e2v(ji+1,jj) - pV(ji,jj) * p1_e2v(ji,jj) ) * pe2e2f(ji,jj) * zzf
            pe12f(ji,jj) = 0.5_wp * ( zS1 + zS2 )    ! eps12 =  1/2 `shearing strain rate` !
            !
         END DO
      END DO

   END SUBROUTINE vel_shear_f


   SUBROUTINE vel_maxshr_t( pU, pV, p1_e1e2t, p1_e1e2f, p1_e1u, p1_e2v, p1_e2u, p1_e1v, pe1e1t, pe2e2t, pe1e1f, pe2e2f, &
      &                     pe1e2f, pmskt, pmskf, pms )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) tension of sea-ice velocity vector => T-point
      !!
      !! Note: the mask@F must be the actual `fmask` that takes into account the slip/no-slip conditions
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t, p1_e1e2f       ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1u, p1_e2v
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2u, p1_e1v   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1t, pe2e2t
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1f, pe2e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt, pmskf         ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pms   ! maximum shear @ T points               [1/s]
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zSHR
      REAL(wp) :: zzf, zten, zs2
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------

      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ! 2*Shear at F points:
            zzf = p1_e1e2f(ji,jj) * pmskf(ji,jj)
            zSHR(ji,jj) = ( ( pU(ji,jj+1) * p1_e1u(ji,jj+1) - pU(ji,jj) * p1_e1u(ji,jj) ) * pe1e1f(ji,jj) * zzf  &
               &         + ( pV(ji+1,jj) * p1_e2v(ji+1,jj) - pV(ji,jj) * p1_e2v(ji,jj) ) * pe2e2f(ji,jj) * zzf  &
               &         ) * p1_e1e2f(ji,jj) * pmskf(ji,jj)

         END DO
      END DO

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Tension at T-points, `dU/dx - dV/dy` :
            zten = (   ( pU(ji,jj)*p1_e2u(ji,jj) - pU(ji-1,jj)*p1_e2u(ji-1,jj) ) * pe2e2t(ji,jj) &
               &     - ( pV(ji,jj)*p1_e1v(ji,jj) - pV(ji,jj-1)*p1_e1v(ji,jj-1) ) * pe1e1t(ji,jj) &
               &    )           * p1_e1e2t(ji,jj)

            !! Shear**2 at T points (doc eq. A16)
            zs2 =  ( zSHR(ji,jj  ) * zSHR(ji,jj  ) * pe1e2f(ji,jj  ) + zSHR(ji-1,jj  ) * zSHR(ji-1,jj  ) * pe1e2f(ji-1,jj  )  &
               &   + zSHR(ji,jj-1) * zSHR(ji,jj-1) * pe1e2f(ji,jj-1) + zSHR(ji-1,jj-1) * zSHR(ji-1,jj-1) * pe1e2f(ji-1,jj-1)  &
               &   ) * 0.25_wp * p1_e1e2t(ji,jj)

            !! Maximum shear rate at T points
            pms(ji,jj) = SQRT( zten*zten + zs2 ) * pmskt(ji,jj)

         END DO
      END DO
   END SUBROUTINE vel_maxshr_t

   !!==============================================================================
END MODULE icedyn_rhg_tools
