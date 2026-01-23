MODULE ice_util
   !!======================================================================
   !!                     ***  MODULE  ice_util  ***
   !!   Sea-Ice dynamics : master routine for rheology
   !!======================================================================
   !! history :  4.0  !  2018     (C. Rousset)      Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!    ...
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE par_ice
   USE ice            ! sea-ice: variables
   USE lib_mpp
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   INTERFACE cap_1md
      MODULE PROCEDURE cap_1md_1g, cap_1md_2g
   END INTERFACE cap_1md

   PUBLIC   round
   !PUBLIC   pos_or_0
   PUBLIC   cap_1md    ! cap `1-d`
   PUBLIC   smooth5p
   PUBLIC   smooth9p

   PUBLIC   smoothCrossTF

   PUBLIC clean_small_a_all

   PUBLIC trace_mean_array_dbg

   REAL(wp), PARAMETER :: rtol_dmg = 0.1_wp   ! tolerance for damage overshoot (above/below 1/0)

   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION round( pr, kn)
      !!--------------------------------------
      !$acc routine
      !!--------------------------------------
      REAL(wp)             :: round
      REAL(wp), INTENT(in) :: pr
      INTEGER , INTENT(in) :: kn
      !!--------------------------------------
      REAL(wp) :: zr
      !!--------------------------------------
      zr = 10._wp**kn
      round = REAL( ANINT( pr * zr ) / zr , wp )
      !!
   END FUNCTION round

   !PURE REAL(wp) FUNCTION pos_or_0( pb, px )
   !   !!--------------------------------------
   !   !$acc routine vector
   !   !!--------------------------------------
   !   !! Returns `px` if `pb  > 0`                                                                                                                            
   !   !! Returns  `0` if `pb <= 0`
   !   !!--------------------------------------
   !   REAL(wp), INTENT(in) :: pb, px
   !   !!--------------------------------------
   !   !IF( pb <= 0._wp ) THEN
   !   !   pos_or_0 = 0._wp
   !   !ELSE
   !   !   pos_or_0 = px
   !   !ENDIF
   !   pos_or_0 = px * ( 1._wp + REAL( pb==0. , wp) )
   !   !      
   !END FUNCTION pos_or_0

   
   SUBROUTINE cap_1md_1g( pA, p1md )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cap_1md  ***
      !! ** Purpose :   Constrain ice damage to sound values
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pA        ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: p1md      ! damage
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------
      !$acc data present( pA, p1md )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            p1md(ji,jj) = MAX( p1md(ji,jj), r_dmd_min )
            p1md(ji,jj) = MIN( p1md(ji,jj),   1._wp   )
            IF( pA(ji,jj) < rAmin_dmg ) p1md(ji,jj) = 1._wp ! => damage is forced to where A<rAmin_dmg
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE cap_1md_1g

   SUBROUTINE cap_1md_2g( pAt, pAf, p1md_t, p1md_f )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cap_1md  ***
      !! ** Purpose :   Constrain ice damage to sound values
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pAt,    pAf     ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: p1md_t, p1md_f  ! damage
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------
      !$acc data present( pAt, p1md_t, pAf, p1md_f )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            p1md_t(ji,jj) = MAX( p1md_t(ji,jj), r_dmd_min )
            p1md_t(ji,jj) = MIN( p1md_t(ji,jj),   1._wp   )
            IF( pAt(ji,jj) < rAmin_dmg ) p1md_t(ji,jj) = 1._wp ! => damage is forced to where A<rAmin_dmg
            !
            p1md_f(ji,jj) = MAX( p1md_f(ji,jj), r_dmd_min )
            p1md_f(ji,jj) = MIN( p1md_f(ji,jj),   1._wp   )
            IF( pAf(ji,jj) < rAmin_dmg ) p1md_f(ji,jj) = 1._wp ! => damage is forced to where A<rAmin_dmg
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE cap_1md_2g


   SUBROUTINE clean_small_a_all( pAt, pAf,  p1mdt, p1mdf,  pSt, pSf )
      !!-------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pAt, pAf     ! ice concentration @T and @F
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf ! ice damage @T and @F
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSt          ! T-centric Sigmas [Pa]
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSf          ! F-centric Sigmas [Pa]
      !!-------------------------------------------------------------------------------------------
      INTEGER :: ji, jj
      !!-------------------------------------------------------------------------------------------
      !$acc data present( pAt, pAf, p1mdt, p1mdf, pSt, pSf )
      !!
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            IF( pAt(ji,jj) < rAmin_dmg ) THEN
               p1mdt(ji,jj) = 1._wp ! => damage=0 where almost no sea-ice left...
               pSt(ji,jj,1) = 0._wp
               pSt(ji,jj,2) = 0._wp
               pSf(ji,jj,3) = 0._wp
            ENDIF
            IF( pAf(ji,jj) < rAmin_dmg ) THEN
               p1mdf(ji,jj) = 1._wp ! => damage=0 max where almost no sea-ice left...
               pSf(ji,jj,1) = 0._wp
               pSf(ji,jj,2) = 0._wp
               pSt(ji,jj,3) = 0._wp
            ENDIF
            !
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE clean_small_a_all





   FUNCTION smooth5p( cgt, px, pwij,  lbcl )
      REAL(wp), DIMENSION(jpi,jpj) :: smooth5p
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt       ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4
      REAL(wp) :: zs, zaa, zbb, zfc
      LOGICAL  :: l_bcl
      INTEGER  :: kh
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      smooth5p(:,:) = 0._wp
      !
      IF    ( cgt=='T') THEN
         zma(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSEIF( cgt=='F') THEN
         zma(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smooth5p(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      kh = nn_hls - 1
      !
      DO jj=Njs0-kh, Nje0+kh
         DO ji=Nis0-kh, Nie0+kh
            it1 = ji+1 ; jt1 = jj
            it2 = ji   ; jt2 = jj+1
            it3 = ji-1 ; jt3 = jj
            it4 = ji   ; jt4 = jj-1
            !!
            zm0 = zma( ji,jj )
            zm1 = zma(it1,jt1)
            zm2 = zma(it2,jt2)
            zm3 = zma(it3,jt3)
            zm4 = zma(it4,jt4)
            !!
            zt0 = px( ji,jj )*zm0
            zt1 = px(it1,jt1)*zm1
            zt2 = px(it2,jt2)*zm2
            zt3 = px(it3,jt3)*zm3
            zt4 = px(it4,jt4)*zm4
            !!
            zs     =        MAX( zaa*zm0  +  zbb*( zm1 + zm2 + zm3 + zm4 ) , 1.E-12_wp ) ! sum of wheights
            !
            smooth5p(ji,jj) = (  zaa*zt0  +  zbb*( zt1 + zt2 + zt3 + zt4 ) ) / zs
            !
         END DO
      END DO
      !
      IF( cgt=='T') THEN
         smooth5p(:,:) = smooth5p(:,:)*xmskt(:,:)
      ELSE
         smooth5p(:,:) = smooth5p(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smooth5p@icedyn_rhg_bbm', smooth5p, cgt, 1._wp )
      !
   END FUNCTION smooth5p


   FUNCTION smooth9p( cgt, px, pwij, lbcl )
      REAL(wp), DIMENSION(jpi,jpj)             :: smooth9p
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt       ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4
      INTEGER  :: it5, jt5, it6, jt6, it7, jt7, it8, jt8
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4, zm5, zm6, zm7, zm8
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8
      REAL(wp) :: zs, zaa, zbb, zd
      LOGICAL  :: l_bcl
      INTEGER  :: kh
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      zd  = 0.7071067811865475_wp ! 1/sqrt(2)
      !
      IF    ( cgt=='T') THEN
         zma(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSEIF( cgt=='F') THEN
         zma(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smooth9p(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      smooth9p(:,:) = 0._wp
      !
      kh = nn_hls - 1
      !
      DO jj=Njs0-kh, Nje0+kh
         DO ji=Nis0-kh, Nie0+kh
            it1 = ji+1 ; jt1 = jj
            it2 = ji   ; jt2 = jj+1
            it3 = ji-1 ; jt3 = jj
            it4 = ji   ; jt4 = jj-1
            !!
            it5 = ji+1 ; jt5 = jj+1
            it6 = ji-1 ; jt6 = jj+1
            it7 = ji-1 ; jt7 = jj-1
            it8 = ji+1 ; jt8 = jj-1
            !!
            zm0 = zma( ji,jj )
            zm1 = zma(it1,jt1)
            zm2 = zma(it2,jt2)
            zm3 = zma(it3,jt3)
            zm4 = zma(it4,jt4)
            !!
            zm5 = zma(it5,jt5) * zd
            zm6 = zma(it6,jt6) * zd
            zm7 = zma(it7,jt7) * zd
            zm8 = zma(it8,jt8) * zd
            !!
            zt0 = px( ji,jj )*zm0
            zt1 = px(it1,jt1)*zm1 ; zt5 = px(it5,jt5)*zm5
            zt2 = px(it2,jt2)*zm2 ; zt6 = px(it6,jt6)*zm6
            zt3 = px(it3,jt3)*zm3 ; zt7 = px(it7,jt7)*zm7
            zt4 = px(it4,jt4)*zm4 ; zt8 = px(it8,jt8)*zm8
            !
            zs     =        MAX( zaa * zm0  +  zbb * ( zm1 + zm2 + zm3 + zm4 + zm5 + zm6 + zm7 + zm8 ) , 1.E-12_wp ) ! sum of wheights
            !
            smooth9p(ji,jj) = (  zaa * zt0  +  zbb * ( zt1 + zt2 + zt3 + zt4 + zt5 + zt6 + zt7 + zt8 ) ) / zs
            !
         END DO
      END DO
      !
      IF( cgt=='T') THEN
         smooth9p(:,:) = smooth9p(:,:)*xmskt(:,:)
      ELSE
         smooth9p(:,:) = smooth9p(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smooth9p@icedyn_rhg_bbm', smooth9p, cgt, 1._wp )
      !
   END FUNCTION smooth9p



   FUNCTION smoothCrossTF( cgt, px, pxE, pwij,  lbcl )
      REAL(wp), DIMENSION(jpi,jpj) :: smoothCrossTF
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt   ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px    ! field on `cgt` grid
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxE   ! counterpart field !
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma0, zmaE
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4, idxp
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4
      REAL(wp) :: zs, zaa, zbb, zfc
      LOGICAL  :: l_bcl
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      smoothCrossTF(:,:) = 0._wp
      !
      IF    ( cgt=='T') THEN
         idxp = 0
         zma0(:,:) = xmskt(:,:)*e1e2t(:,:)
         zmaE(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSEIF( cgt=='F') THEN
         idxp = 1
         zma0(:,:) = xmskf(:,:)*e1e2f(:,:)
         zmaE(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smoothCrossTF(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            it1 = ji+idxp   ; jt1 = jj+idxp
            it2 = ji-1+idxp ; jt2 = jj+idxp
            it3 = ji-1+idxp ; jt3 = jj-1+idxp
            it4 = ji+idxp   ; jt4 = jj-1+idxp
            !IF on F-grid:
            !it1 = ji+1 ; jt1 = jj+1
            !it2 = ji   ; jt2 = jj+1
            !it3 = ji   ; jt3 = jj
            !it4 = ji+1 ; jt4 = jj
            !!
            zm0 = zma0( ji,jj )
            zm1 = zmaE(it1,jt1)
            zm2 = zmaE(it2,jt2)
            zm3 = zmaE(it3,jt3)
            zm4 = zmaE(it4,jt4)
            !!
            zt0 =  px( ji,jj )*zm0
            zt1 = pxE(it1,jt1)*zm1
            zt2 = pxE(it2,jt2)*zm2
            zt3 = pxE(it3,jt3)*zm3
            zt4 = pxE(it4,jt4)*zm4
            !!
            zs     =        MAX( zaa*zm0  +  zbb*( zm1 + zm2 + zm3 + zm4 ) , 1.E-12_wp ) ! sum of wheights
            !
            smoothCrossTF(ji,jj) = (  zaa*zt0  +  zbb*( zt1 + zt2 + zt3 + zt4 ) ) / zs
            !
         END DO
      END DO
      !
      IF( cgt=='T') THEN
         smoothCrossTF(:,:) = smoothCrossTF(:,:)*xmskt(:,:)
      ELSE
         smoothCrossTF(:,:) = smoothCrossTF(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smoothCrossTF@icedyn_rhg_bbm', smoothCrossTF, cgt, 1._wp )
      !
   END FUNCTION smoothCrossTF




   SUBROUTINE div_stress_tensor_v2( cgt,  pe2x, pe1y, pe1e1, pe2e2,  pe1e1_e, pe2e2_e,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y,  &
      &                               ps11h, ps22h, ps12h,  pdivSx, pdivSy )
      !! => use the same discretization approach as what's being done in EVP, using `sigma_1` and `sigma_2` rather than `sigma_11` and `sigma_22`...
      !!----------------------------------------------------------------------------------------------
      !! Computes the vector (pdivSx,pdivSy) = divergence of the h-integrated internal stress tensor
      !!
      !!   depending on the grid: T-centric grid => cgt='T' or F-centric grid => cgt='F'
      !!
      !! INPUT:                                               |     cgt=='T'   |    cgt=='F'    |
      !!   * ps11h, ps22h: sigma11*h, sigma22*h           =>  ! @ point T[i,j] | @ point F[i,j] |
      !!   * ps12h       :       sigma12*h                =>  ! @ point F[i,j] | @ point T[i,j] |
      !!
      !! RETURNS:                                             |     cgt=='T'   |    cgt=='F'    |
      !!   * pdivSx: x-component of the div of the tensor =>  | @ point U[i,j] | @ point V[i,j] |
      !!   * pdivSy: y-component of the div of the tensor =>  | @ point V[i,j] | @ point U[i,j] |
      !!
      !!----------------------------------------------------------------------------------------------
      CHARACTER(len=1),         INTENT(in)  :: cgt
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2x, pe1y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2, pe1e1_e, pe2e2_e
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2x, p1_e1e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11h, ps22h, ps12h ! components of stress tensors on T- or F-centric grids x h !!!
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdivSx, pdivSy      ! x,y components of the divergence of the tensor
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zs1, zs2, zs3
      INTEGER  :: ip, im, jp, jm, ji, jj
      !!--------------------------------------------------------------------------------------------
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
         CALL ctl_stop( 'STOP', 'div_stress_tensor_v2(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      pdivSx(:,:) = 0._wp
      pdivSy(:,:) = 0._wp
      !
      zs1(:,:) = ps11h(:,:) + ps22h(:,:) ! h * sigma_1
      zs2(:,:) = ps11h(:,:) - ps22h(:,:) ! h * sigma_2
      zs3(:,:) = 2._wp * ps12h(:,:)      ! 2 * h * sigma_12
      !
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !                   !--- U points
            pdivSx(ji,jj) = 0.5_wp * ( (( zs1(ji+ip,jj)                   - zs1(ji+im,jj)                  ) *   pe2x(ji,jj)  &
               &                      + ( zs2(ji+ip,jj)*pe2e2(ji+ip,jj)   - zs2(ji+im,jj)*pe2e2(ji+im,jj)  ) * p1_e2x(ji,jj)) &
               &                      + ( zs3(ji,jj+jp)*pe1e1_e(ji,jj+jp) - zs3(ji,jj+jm)*pe1e1_e(ji,jj+jm)) * p1_e1x(ji,jj)  &
               &                      ) * p1_e1e2x(ji,jj)
            !
            !                !--- V points
            pdivSy(ji,jj) = 0.5_wp * ( (( zs1(ji,jj-jm)                   - zs1(ji,jj-jp)                  ) *   pe1y(ji,jj)  &
               &                      - ( zs2(ji,jj-jm)*pe1e1(ji,jj-jm)   - zs2(ji,jj-jp)*pe1e1(ji,jj-jp)  ) * p1_e1y(ji,jj)) &
               &                      + ( zs3(ji-im,jj)*pe2e2_e(ji-im,jj) - zs3(ji-ip,jj)*pe2e2_e(ji-ip,jj)) * p1_e2y(ji,jj)  &
               &                      ) * p1_e1e2y(ji,jj)
            !
         END DO
      END DO
      !
   END SUBROUTINE div_stress_tensor_v2





   SUBROUTINE trace_mean_array_dbg( cname, kt, kts, pmsk, pX )
      !!-------------------------------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in)     :: cname
      INTEGER,                    INTENT(in)     :: kt, kts
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in) :: pmsk  ! mask to apply
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in) :: pX  ! array to trace
      !!-------------------------------------------------------------------------------------------
      REAL(wp) :: zs, zw
      INTEGER :: ji, jj
      !!-------------------------------------------------------------------------------------------
      !
      zw = 0._wp
      !
      !$acc parallel loop collapse(2) present( pmsk, pX )
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            zs = zs + pmsk(ji,jj)
            zw = zw + pmsk(ji,jj)*pX(ji,jj)
            !
         END DO
      END DO
      !$acc end parallel loop
      WRITE(*,'(" *** Mean val of ",a," at kt, kts =",i3.3,", ",i3.3," =>",f)') TRIM(cname), kt, kts, REAL(zw/zs,4)
      !!
      !!
   END SUBROUTINE trace_mean_array_dbg




   !!======================================================================
END MODULE ice_util
