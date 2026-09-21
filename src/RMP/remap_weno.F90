MODULE remap_weno

   USE par_kind, ONLY: wp

   USE lib_mpp,  ONLY: ctl_stop

   USE in_out_manager, ONLY: numout, lwp, numoni

   USE iom            ! I/O library

   USE dom_oce, ONLY : xmskt, xmskf

   USE par_ice, ONLY : ln_damage, ln_use_weno_rmp, cn_rmp_weno_wght, epsi10, epsi20, rn_amax

   USE ice,     ONLY : weno_s_lw_t_x, weno_s_lw_t_y, weno_s_ow_t_x, weno_s_ow_t_y,   &
      &                weno_s_lw_f_x, weno_s_lw_f_y, weno_s_ow_f_x, weno_s_ow_f_y

   USE timing

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: remap_weno_init

   PUBLIC :: weno5_ISx_G
   PUBLIC :: weno7_ISx_G

   PUBLIC :: weno5_sym_ISx_G

   !PUBLIC :: rmpT2U_wn5s
   !PUBLIC :: rmpT2V_wn5s
   !PUBLIC :: rmpT2F_wn5s

   PUBLIC :: rmpT2F_A_h_wn5s   ! interpolate both `A` & `h` from T- to F-points...
   PUBLIC :: rmpT2UV_wn5s      ! interpolate a scalar defined at T to both U- and V-points...
   !PUBLIC :: rmpT2F_A_h_v2_wn5s   ! interpolate both `A` & `h` from T- to F-points...


   INTEGER, PARAMETER, PUBLIC :: ko_wrmp = 5   ! order of WENO scheme for remapping

   INTEGER, PARAMETER, PUBLIC :: kp_wrmp = (ko_wrmp + 1)/2   ! order of each of the `kp_wrmp+1` stencils used

   LOGICAL, PARAMETER :: ll_SI_wght = .TRUE. ! Weight the optimal weights with the smoothness indicators

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   FUNCTION rmpT2U_wn5s( pxt,  lbcl )
      !!--------------------------------------------------------------
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2U_wn5s
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zxu
      !!--------------------------------------------------------------
      CALL ctl_stop( 'STOP', 'rmpT2U_wn5s: => make me GPU compliant !!!' )
      rmpT2U_wn5s(:,:) = 0._wp
      !
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pxt, zxu )
      !
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2U_wn5s@remap_weno.F90', zxu,'U',1._wp )
      END IF
      !
      rmpT2U_wn5s(:,:) = zxu(:,:)
      !
   END FUNCTION rmpT2U_wn5s



   SUBROUTINE rmpT2UV_wn5s( pxt, pxu, pxv,  vminmax, lbcl )
      !!--------------------------------------------------------------
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in)  :: pxt
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(out) :: pxu
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(out) :: pxv
      REAL(wp), DIMENSION(2), OPTIONAL, INTENT(in)  :: vminmax
      LOGICAL,                OPTIONAL, INTENT(in)  :: lbcl
      !!--------------------------------------------------------------
      INTEGER  :: ji, jj
      LOGICAL  :: l_min_max
      REAL(wp) :: zmin, zmax
      !!--------------------------------------------------------------
      !$acc data present( pxt, pxu, pxv )
      l_min_max = PRESENT( vminmax )

      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pxt, pxu )
      CALL rmp_wn5s_y( klbct, e2t, weno_s_lw_t_y, weno_s_ow_t_y, pxt, pxv )
      
      IF( l_min_max ) THEN
         zmin = vminmax(1) ; zmax = vminmax(2)
         !$acc parallel loop collapse(2) copyin(zmin,zmax)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               pxu(ji,jj) = MIN( MAX( pxu(ji,jj), zmin), zmax)
               pxv(ji,jj) = MIN( MAX( pxv(ji,jj), zmin), zmax)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF
      !
#if ! defined _OPENACC || defined _OPENMP
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2UV_wn5ss@remap_weno', pxu,'U',1._wp, pxv,'V',1._wp ) ! Needed!
      END IF
#endif
      !
      !$acc end data
   END SUBROUTINE rmpT2UV_wn5s




















   FUNCTION rmpT2V_wn5s( pxt,  lbcl )
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2V_wn5s
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zxv
      !!--------------------------------------------------------------
      CALL ctl_stop( 'STOP', 'rmpT2V_wn5s: => make me GPU compliant !!!' )
      rmpT2V_wn5s(:,:) = 0._wp
      !
      CALL rmp_wn5s_y( klbct, e2t, weno_s_lw_t_y, weno_s_ow_t_y, pxt, zxv )
      !
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2V_wn5s@remap_weno.F90', zxv,'V',1._wp )
      END IF
      !
      rmpT2V_wn5s(:,:) = zxv(:,:)
      !
   END FUNCTION rmpT2V_wn5s



   FUNCTION rmpT2F_wn5s( pxt,  lbcl )
      !!---------------------------------------------------------------------------
      !! WENO symetric 5th order generalized to orthogonal curvilinear coordinates
      !!  ==> used for high-order interpolation of fields from T to F (or F to T) points
      !!---------------------------------------------------------------------------
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2F_wn5s
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zxu, zxf
      !!--------------------------------------------------------------
      CALL ctl_stop( 'STOP', 'rmpT2F_wn5s: => make me GPU compliant !!!' )
      rmpT2F_wn5s(:,:) = 0._wp
      !
      ! First interp `pxt` from T[i] to U[i] points:
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pxt, zxu )
      CALL lbc_lnk( 'rmpT2F_wn5s@remap_weno.F90', zxu,'U',1._wp )

      ! Now interp `zxu` from U[i] to F[i] points:
      CALL rmp_wn5s_y( klbcu, e2u, weno_s_lw_t_y, weno_s_ow_t_y, zxu, zxf ) !LOLO: need weights for U point ??
      !       "vu depuis `U` le `F` est l'equivalent du `V` vu depuis `T`", hence the `t` here...

      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F_wn5s@remap_weno.F90', zxf,'F',1._wp )
      END IF

      rmpT2F_wn5s(:,:) = zxf(:,:) * xmskf(:,:)

   END FUNCTION rmpT2F_wn5s




   SUBROUTINE rmpT2F_A_h_wn5s( pAt, pht, pAf, phf,  lbcl )
      !!--------------------------------------------------------------
      !! ==> pAt & pht have to be LBC_LNKed !!!!
      !!
      !! zAu => xtmp1
      !! zhu => xtmp2
      !! zAf => xtmp3
      !! zhf => xtmp4
      !!
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pAt, pht
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pAf, phf
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl
      !!--------------------------------------------------------------
      INTEGER  :: ji, jj
      !!--------------------------------------------------------------
      !IF( ln_timing ) CALL timing_start('rmpT2F_A_h_wn5s')
      !$acc data present( pAt, pht, pAf, phf, xtmp1, xtmp2, xtmp3, xtmp4 )

      ! First interp from T[i] to U[i] points:
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pAt, xtmp1 )
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pht, xtmp2 )
#if ! defined _OPENACC || defined _OPENMP
      CALL lbc_lnk( 'rmpT2F_A_h_wn5s@remap_weno', xtmp1,'U',1._wp,  xtmp2,'U',1._wp )
#endif
      !
      ! Now interp from U[i] to F[i] points:
      CALL rmp_wn5s_y( klbcu, e2u, weno_s_lw_t_y, weno_s_ow_t_y, xtmp1, xtmp3 ) !LOLO: need weights for U point ??
      CALL rmp_wn5s_y( klbcu, e2u, weno_s_lw_t_y, weno_s_ow_t_y, xtmp2, xtmp4 ) !LOLO: need weights for U point ??
      !
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            pAf(ji,jj) = MIN( MAX( xtmp3(ji,jj), 0._wp), rn_amax) * xmskf(ji,jj)
            phf(ji,jj) =      MAX( xtmp4(ji,jj), 0._wp)           * xmskf(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      !
#if ! defined _OPENACC || defined _OPENMP
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F_A_h_wn5s@remap_weno', pAf,'F',1._wp, phf,'F',1._wp ) ! Needed!
      END IF
#endif
      !
      !$acc end data
      !IF( ln_timing ) CALL timing_stop('rmpT2F_A_h_wn5s')
      !
   END SUBROUTINE rmpT2F_A_h_wn5s


   SUBROUTINE rmpT2F_A_h_v2_wn5s( pAt, pVt, pAf, phf,  lbcl )
      !!--------------------------------------------------------------
      !! ==> pAt & pVt have to be LBC_LNKed !!!!
      !!
      !! zAu => xtmp1
      !! zVu => xtmp2
      !! zAf => xtmp3
      !! zVf => xtmp4
      !!
      !!--------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pAt, pVt  ! Ice concentration & volume at T-points [-], [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pAf, phf
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl
      !!--------------------------------------------------------------
      INTEGER  :: ji, jj
      REAL(wp) :: zAf
      !!--------------------------------------------------------------
      !$acc data present( pAt, pVt, pAf, phf, xtmp1, xtmp2, xtmp3, xtmp4 )

      ! First interp from T[i] to U[i] points:
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pAt, xtmp1 )
      CALL rmp_wn5s_x( klbct, e1t, weno_s_lw_t_x, weno_s_ow_t_x, pVt, xtmp2 )
#if ! defined _OPENACC || defined _OPENMP
      CALL lbc_lnk( 'rmpT2F_A_h_v2_wn5s@remap_weno', xtmp1,'U',1._wp,  xtmp2,'U',1._wp )
#endif
      !
      ! Now interp from U[i] to F[i] points:
      CALL rmp_wn5s_y( klbcu, e2u, weno_s_lw_t_y, weno_s_ow_t_y, xtmp1, xtmp3 ) !LOLO: need weights for U point ??
      CALL rmp_wn5s_y( klbcu, e2u, weno_s_lw_t_y, weno_s_ow_t_y, xtmp2, xtmp4 ) !LOLO: need weights for U point ??
      !
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            zAf = MIN( MAX( xtmp3(ji,jj), 0._wp), rn_amax) * xmskf(ji,jj)
            phf(ji,jj) = MERGE( MAX(xtmp4(ji,jj),0._wp) / MAX(zAf,epsi20)  ,  0._wp  , zAf>epsi20 )
            pAf(ji,jj) = zAf
         END DO
      END DO
      !$acc end parallel loop

#if ! defined _OPENACC || defined _OPENMP
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F_A_h_v2_wn5s@remap_weno', pAf,'F',1._wp, phf,'F',1._wp ) ! Needed!
      END IF
#endif

      !$acc end data
   END SUBROUTINE rmpT2F_A_h_v2_wn5s!










#include "remap_wn5_sym_x.h90"

#include "remap_wn5_sym_y.h90"




   FUNCTION weno5_ISx_G( pF, pFm, pFp )
      !!=======================================================================
      !!
      !! Computes the "smoothness indicators" for WENO 5 in the generalized
      !! case of a curvilinear grid, i.e. Cartesian coordinates with a variying
      !! mesh size `h` (aka `dx` or `dy`) !
      !!
      !! Author: Laurent Brodeau, February 2025
      !!
      !! => I used symbolic math software to do all the algebra & development
      !!    that leads to the following expression of
      !!    IS[l] (smoothness indicators) for WENO5
      !!
      !!=======================================================================
      !%acc routine
      !!=======================================================================
      REAL(wp)                :: weno5_ISx_G ! Smoothness indicator
      REAL(wp),    INTENT(in) :: pF          ! scalar field
      REAL(wp),    INTENT(in) :: pFm, pFp    ! WENO `p`th order (and not `W`th order!) interp of F at i-1/2 & i+1/2 ()
      !!=======================================================================
      REAL(wp) :: zF2, zFm2, zFp2, zS
      !!=======================================================================
      zF2  = pF*pF
      zFm2 = pFm*pFm
      zFp2 = pFp*pFp
      zS   = pF*(pFm + pFp)
      !
      weno5_ISx_G = 4._wp*( 39._wp*zF2 + 10._wp*zFm2 + 10._wp*zFp2 + 19._wp*pFm*pFp - 39._wp*zS )
      !
   END FUNCTION weno5_ISx_G


   FUNCTION weno7_ISx_G( kl, prh_m, prh_p, pFkm1, pFk, pFkp1, pFm, pFp )
      !!=======================================================================
      !!
      !! Computes the "smoothness indicators" for WENO 7 in the generalized
      !! case of a curvilinear grid, i.e. Cartesian coordinates with a variying
      !! mesh size `h` (aka `dx` or `dy`) !
      !!
      !! Author: Laurent Brodeau, February 2025
      !!
      !! => I used symbolic math software to do all the algebra & development
      !!    that leads to the following expression of
      !!    IS[l] (smoothness indicators) for WENO7
      !!
      !!=======================================================================
      !$acc routine
      !!=======================================================================
      REAL(wp)              :: weno7_ISx_G  ! Smoothness indicator for `l=kl`
      INTEGER,  INTENT(in)  :: kl           ! `l` stencil bias: kl=[0:3] for WENO7
      REAL(wp), INTENT(in)  :: prh_m, prh_p ! `h` local scale factor of mesh size
      !                                     ! ==> prh_m = e1(i-1)/e1(i), prh_p = e1(i+1)/e1(i) or e2(j-1)/e2(j), prh_p = e2(j+1)/e2(j)
      REAL(wp), INTENT(in)  :: pFkm1, pFk, pFkp1 ! scalar field          (k-1,k,k+1)
      REAL(wp), INTENT(in)  :: pFm, pFp     ! WENO interp of F at i-1/2 & i+1/2 (are also functions of `h`)
      !!=======================================================================
      REAL(wp) :: zQi2, zqip2, zqim2, zqimp, zqq, zA, zB, zC, zD, zE, zfct, zX, zY, zZ
      REAL(wp) :: zQiN, zqim, zqip, zeh, ze2, zdif, zj
      !!=======================================================================

      zQi2  = pFk*pFk
      zqimp = pFm + pFp
      zqq   = pFm*pFp
      zA    = 4294._wp*pFk*zqimp
      zB    =   39._wp*pFk*zqimp
      zC    = 1562._wp*zqimp

      IF( kl< 2) THEN
         ! => `pFkp1` & `prh_p` not used
         zQiN = pFkm1
         zqim = pFm
         zqip = pFp
         zeh  = prh_m
         !
      ELSE
         ! => `pFkm1` & `prh_m` not used
         zQiN = pFkp1
         zqim = pFp
         zqip = pFm
         zeh  = prh_p
         !
      ENDIF

      zqip2 = zqip*zqip
      zqim2 = zqim*zqim
      zD    = 4294._wp*zQi2 + 1081._wp*zqim2 + 2132._wp*zqq + 1081._wp*zqip2 - zA
      zE    =   39._wp*zQi2 +   10._wp*zqim2 +   19._wp*zqq +   10._wp*zqip2 - zB
      ze2   = zeh*zeh
      zdif  = zqim - zQiN
      zj    = 1._wp + zeh

      zfct = 4._wp / ( 5._wp * ze2 * zj*zj*zj*zj )

      zX =   ze2*( pFk*(-12691._wp*zqim + 3124._wp*zQiN - 4881._wp*zqip) + 4781._wp*zqq &
         & + 831._wp*zqip2 + 7224._wp*zQi2 - zC*zQiN + 4736._wp*zqim2 )

      zY =   2._wp*ze2*zeh*(5076._wp*zQi2 - 4295._wp*pFk*zqip - 5857._wp*pFk*zqim &
         & + 1662._wp*zqim2 + 2533._wp*zqq + 881._wp*zqip2)

      zZ = ze2*ze2*( zD + 20._wp*zeh*zE  + 5._wp*ze2*zE )

      weno7_ISx_G = zfct * ( 781._wp*zdif*zdif + 1562._wp*zeh*zdif*(2._wp*zqim + zqip -3._wp*pFk) + zX + zY + zZ )

   END FUNCTION weno7_ISx_G


   FUNCTION weno5_sym_ISx_G( kl, prh_m, prh_p, pFkm1, pFk, pFkp1, pFp )
      !!=======================================================================
      !!
      !! Computes the "smoothness indicators" for the symetric WENO 5
      !! in the generalized case of a curvilinear grid, i.e. Cartesian
      !! coordinates WITH a variying mesh size `h` (aka `dx` or `dy`) !
      !!
      !! Author: Laurent Brodeau, March 2025
      !!
      !! => I used symbolic math software to do all the algebra & development
      !!    that leads to the following expression of
      !!    IS[l] (smoothness indicators) for WENO5-sym
      !!
      !!=======================================================================
      !$acc routine
      !!=======================================================================
      REAL(wp)            :: weno5_sym_ISx_G   ! Smoothness indicator for stencil # `kl`
      INTEGER, INTENT(in) :: kl                ! which of the 4 stencils kl=[1:4]
      REAL(wp),INTENT(in) :: prh_m, prh_p      ! `h` local scale factor of mesh size
      !                                        ! ==> prh_m = e1(i-1)/e1(i), prh_p = e1(i+1)/e1(i) or e2(j-1)/e2(j), prh_p = e2(j+1)/e2(j)
      REAL(8), INTENT(in) :: pFkm1, pFk, pFkp1 ! scalar field   (k-1,k,k+1)
      REAL(8), INTENT(in) :: pFp               ! interp of F at i+1/2 by stencil # `kl`
      !!=======================================================================
      INTEGER :: ki
      REAL(8) :: em, ep, qip
      REAL(8) :: ep2, ep3, ep4, zA, zA2, zB, zC, zD, zE, em2
      !!=======================================================================

      !IF( (kl>4).OR.(kl<1) ) CALL ctl_stop( 'STOP', 'weno5_sym_ISx_G: only `1 <= kl <= 4` allowed' )

      ki = 2

      !Qi  =  pF(ki)

      qip = pFp

      IF( kl < 3 ) THEN

         em   = prh_m
         !Qim1 = pF(ki-1)

         em2 = em*em

         zA = 1._wp+ em
         zB = 2._wp+ em

         zC =          -zB      *pFk +    pFkm1 +        zA *qip
         zC = SIGN( 1._wp, zC ) * MAX( ABS(zC) ,epsi20 )  ! zC without the 0 singularity
         zD =       (3._wp- em2)  *pFk - 2._wp*pFkm1 + (em2 - 1.)*qip
         zE = (3._wp+ 3._wp*em + em2)*pFk -    pFkm1 -      zA*zB*qip

         weno5_sym_ISx_G = ( 4._wp*( 81._wp*zC*zC - (zD*zD*zD + zE*zE*zE)/zC ) ) / (9._wp*zA*zA*zA*zA)

      ELSE

         ep   = prh_p
         !Qip1 = pF(ki+1)

         ep2 = ep**2
         ep3 = ep2*ep
         ep4 = ep3*ep

         zA  = pFk - qip
         zA2 = zA*zA
         zB  = qip - pFkp1

         zC  = 12._wp*pFk - 13._wp*qip + pFkp1

         zD  = ep2 + 2._wp*ep3 + ep4

         weno5_sym_ISx_G = ( 4._wp*( 3._wp*ep3*zA2 + ep4*zA2 - 21._wp*ep*zA*zB + 10._wp*zB*zB + ep2*zA*zC ) ) / zD

      ENDIF

   END FUNCTION weno5_sym_ISx_G








   SUBROUTINE read_weno_sym_w( kpi, kpj, kpp, kwo, cfilew, pLWtx, pOWtx, pLWty, pOWty, &
      &                                                    pLWfx, pOWfx, pLWfy, pOWfy  )
      !!================================================================================================
      !!================================================================================================
      INTEGER,                                            INTENT(in)  :: kpi, kpj, kpp, kwo
      CHARACTER(len=*),                                   INTENT(in)  :: cfilew
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo+1),           INTENT(out) :: pLWtx
      REAL(8),    DIMENSION(kpi,kpj,kpp),                 INTENT(out) :: pOWtx
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo+1),           INTENT(out) :: pLWty
      REAL(8),    DIMENSION(kpi,kpj,kpp),                 INTENT(out) :: pOWty
      !
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo+1), OPTIONAL, INTENT(out) :: pLWfx
      REAL(8),    DIMENSION(kpi,kpj,kpp),       OPTIONAL, INTENT(out) :: pOWfx
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo+1), OPTIONAL, INTENT(out) :: pLWfy
      REAL(8),    DIMENSION(kpi,kpj,kpp),       OPTIONAL, INTENT(out) :: pOWfy
      !!================================================================================================
      !!
      CHARACTER(len=13), PARAMETER :: crtn = 'read_weno_sym_w'
      !!
      INTEGER :: inum, js
      CHARACTER(len=64) :: cvar
      CHARACTER(len=13) :: cvlw    ! ex: `w5s_lw_t_01_x`
      CHARACTER(len=1)  :: cwo, cp
      LOGICAL           :: lReadFp ! read weights at F-points?
      !!================================================================================================
      lReadFp = ( PRESENT(pLWfx) .AND. PRESENT(pOWfx) .AND. PRESENT(pLWfy) .AND. PRESENT(pOWfy) )

      WRITE(cwo,'(i1)') kwo

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' * Will read linear and optimal weights for symetric WENO'//cwo
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   => in file ', TRIM(cfilew)
      ENDIF

      CALL iom_open( cfilew, inum )

      cp = 't'

      pLWtx = 0._wp
      pOWtx = 0._wp
      pLWty = 0._wp
      pOWty = 0._wp

      cvar = 'weno'//cwo//'_ow_'//cp//'_x'
      IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
      CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvar), pOWtx(:,:,:) )

      cvar = 'weno'//cwo//'_ow_'//cp//'_y'
      IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
      CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvar), pOWty(:,:,:) )

      DO js = 1, kwo+1

         WRITE(cvlw,'("w",i1,"s_lw_",a1,"_",i2.2,"_x")') kwo, cp, js   ! Name of linear weight array to read
         IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
         CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvlw), pLWtx(:,:,:,js) )

         WRITE(cvlw,'("w",i1,"s_lw_",a1,"_",i2.2,"_y")') kwo, cp, js   ! Name of linear weight array to read
         IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
         CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvlw), pLWty(:,:,:,js) )

      END DO

      IF( lReadFp ) THEN

         cp = 'f'

         pLWfx = 0._wp
         pOWfx = 0._wp
         pLWfy = 0._wp
         pOWfy = 0._wp

         cvar = 'weno'//cwo//'_ow_'//cp//'_x'
         IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
         CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvar), pOWfx(:,:,:) )

         cvar = 'weno'//cwo//'_ow_'//cp//'_y'
         IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
         CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvar), pOWfy(:,:,:) )

         DO js = 1, kwo+1

            WRITE(cvlw,'("w",i1,"s_lw_",a1,"_",i2.2,"_x")') kwo, cp, js   ! Name of linear weight array to read
            IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
            CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvlw), pLWfx(:,:,:,js) )

            WRITE(cvlw,'("w",i1,"s_lw_",a1,"_",i2.2,"_y")') kwo, cp, js   ! Name of linear weight array to read
            IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
            CALL iom_get( 'read_weno_sym_w', inum, jpdom_global, TRIM(cvlw), pLWfy(:,:,:,js) )

         END DO

      ENDIF !IF( lReadFp )

      CALL iom_close( inum )

      IF(lwp) WRITE(numout,*) ' * Linear and Optimal weights for symetric WENO'//cwo//' read !'
      IF(lwp) WRITE(numout,*) ''

   END SUBROUTINE read_weno_sym_w

   SUBROUTINE remap_weno_init()
      !!================================================================================================
      !!  * `kwo` is the order of the cuurent WENO scheme: 5 for now...
      !!  * `p` is the order of each of the `p` `p`-point stencils used
      !!  *        `w = 2*p -1`
      !!
      !!  => WENO5-symetric: `kwo=5`, `p=3` ==> 5th order accuracy, combines 4 3rd order 3-point stencils.
      !!================================================================================================
      INTEGER :: nS, nP, ierr
      INTEGER :: ios
      !!
      NAMELIST/namrmp/ ln_use_weno_rmp, cn_rmp_weno_wght
      !!================================================================================================
      !
      READ_NML_REF(numnam_ice,namrmp)
      READ_NML_CFG(numnam_ice,namrmp)
      IF(lwm) WRITE( numoni, namrmp )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'remap_weno_init: parameters remapping '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namrmp:'
         WRITE(numout,*) '      use WENO-sym p2p remapping                   ln_use_weno_rmp  = ', ln_use_weno_rmp
         WRITE(numout,*) '      file containing WENO-symetric weights        cn_rmp_weno_wght = ', TRIM(cn_rmp_weno_wght)
         WRITE(numout,*)
      ENDIF

      IF( ln_use_weno_rmp ) THEN

         IF( nn_hls<3 ) CALL ctl_stop( 'STOP', 'remap_weno_init: use a `nn_hls` of at least 3 for WENO5-sym interpolation.' )

         nS = kp_wrmp + 1
         nP = ko_wrmp + 1

         IF(lwp) THEN                     ! control print
            WRITE(numout,*)
            WRITE(numout,*) 'remap_weno_init: initialization of WENO symetric scheme for p2p interpolation'
            WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~'
            WRITE(numout,*) '      * order of the scheme: `ko_wrmp` =', ko_wrmp
            WRITE(numout,*) '                    ==>             `p` =', kp_wrmp
            WRITE(numout,*) '      * allocating linear and optimal weight arrays'
         ENDIF

         !! Allocation of linear and optimal weights (curvilinear grid!):
         ALLOCATE( weno_s_lw_t_x(jpi,jpj,nS,nP), weno_s_ow_t_x(jpi,jpj,nS),  &
            &      weno_s_lw_t_y(jpi,jpj,nS,nP), weno_s_ow_t_y(jpi,jpj,nS), STAT=ierr )
         IF( ierr/=0 ) CALL ctl_stop( 'STOP', 'remap_weno_init: failed to allocate WENO T-point arrays.' )
         IF( ln_damage ) THEN
            ALLOCATE( weno_s_lw_f_x(jpi,jpj,nS,nP), weno_s_ow_f_x(jpi,jpj,nS),  &
               &      weno_s_lw_f_y(jpi,jpj,nS,nP), weno_s_ow_f_y(jpi,jpj,nS), STAT=ierr )
            IF( ierr/=0 ) CALL ctl_stop( 'STOP', 'remap_weno_init: failed to allocate WENO F-point arrays.' )
         ENDIF
         IF(lwp) THEN
            WRITE(numout,*) '                    ==>             done!'
            WRITE(numout,*) '      * now filling them with weights found in file '//TRIM(cn_rmp_weno_wght)
         ENDIF

         IF( ln_damage ) THEN
            CALL read_weno_sym_w( jpi, jpj, nS, ko_wrmp, cn_rmp_weno_wght, weno_s_lw_t_x, weno_s_ow_t_x, weno_s_lw_t_y, weno_s_ow_t_y, &
               &                                                           weno_s_lw_f_x, weno_s_ow_f_x, weno_s_lw_f_y, weno_s_ow_f_y )
            CALL lbc_lnk('remap_weno_init', weno_s_lw_t_x,'T',1._wp, weno_s_lw_t_y,'T',1._wp, weno_s_lw_f_x,'F',1._wp, weno_s_lw_f_y,'F',1._wp )
            CALL lbc_lnk('remap_weno_init', weno_s_ow_t_x,'T',1._wp, weno_s_ow_t_y,'T',1._wp, weno_s_ow_f_x,'F',1._wp, weno_s_ow_f_y,'F',1._wp )
         ELSE
            CALL read_weno_sym_w( jpi, jpj, nS, ko_wrmp, cn_rmp_weno_wght, weno_s_lw_t_x, weno_s_ow_t_x, weno_s_lw_t_y, weno_s_ow_t_y )
            CALL lbc_lnk('remap_weno_init', weno_s_lw_t_x,'T',1._wp, weno_s_lw_t_y,'T',1._wp )
            CALL lbc_lnk('remap_weno_init', weno_s_ow_t_x,'T',1._wp, weno_s_ow_t_y,'T',1._wp )
         ENDIF

         IF(lwp) THEN
            WRITE(numout,*) '                    ==>             done!'
            WRITE(numout,'("      Symetric WENO",i1," initialization done for remapping :D")') ko_wrmp
            WRITE(numout,*) ''
         ENDIF

#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: remap_weno_init() => adding the 8 `weno_s_*w_*_*` arrays to memory!'
         !$acc enter data copyin(weno_s_lw_t_x,weno_s_lw_t_y,weno_s_ow_t_x,weno_s_ow_t_y,weno_s_lw_f_x,weno_s_lw_f_y,weno_s_ow_f_x,weno_s_ow_f_y)
#endif

         !! Some debug control:
         !IF(lwp) THEN
         !   WRITE(numout,*) ''
         !   WRITE(numout,*) ' *** WENO weights at T-points ***'
         !   CALL PRINT_LINEAR_W(  ko_wrmp, weno_s_lw_t_x(jpi/2,jpj/2,:,:),  'for X, at a given T-point', fmult=6._wp )
         !   CALL PRINT_LINEAR_W(  ko_wrmp, weno_s_lw_t_y(jpi/2,jpj/2,:,:),  'for Y, at a given T-point', fmult=6._wp )
         !   CALL PRINT_OPTIMAL_W( ko_wrmp, weno_s_ow_t_x(jpi/2,jpj/2,:),    'for X, at a given T-point' )
         !   CALL PRINT_OPTIMAL_W( ko_wrmp, weno_s_ow_t_y(jpi/2,jpj/2,:),    'for Y, at a given T-point' )
         !   WRITE(numout,*) ''
         !   IF( ln_damage ) THEN
         !      WRITE(numout,*) ''
         !      WRITE(numout,*) ' *** WENO weights at F-points ***'
         !      CALL PRINT_LINEAR_W(  ko_wrmp, weno_s_lw_f_x(jpi/2,jpj/2,:,:),  'for X, at a given F-point', fmult=6._wp )
         !      CALL PRINT_LINEAR_W(  ko_wrmp, weno_s_lw_f_y(jpi/2,jpj/2,:,:),  'for Y, at a given F-point', fmult=6._wp )
         !      CALL PRINT_OPTIMAL_W( ko_wrmp, weno_s_ow_f_x(jpi/2,jpj/2,:),    'for X, at a given F-point' )
         !      CALL PRINT_OPTIMAL_W( ko_wrmp, weno_s_ow_f_y(jpi/2,jpj/2,:),    'for Y, at a given F-point' )
         !      WRITE(numout,*) ''
         !   ENDIF
         !ENDIF

      ENDIF !IF( ln_use_weno_rmp )

   END SUBROUTINE remap_weno_init


END MODULE remap_weno
