MODULE remap_classic
   !!======================================================================
   !!                     ***  MODULE  remap_classic  ***
   !!======================================================================
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!    ...
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE par_oce,  ONLY: jpi, jpj, Nis0, Nie0, Njs0, Nje0, nn_hls
   USE dom_oce
   USE lib_mpp,  ONLY: ctl_stop

   USE in_out_manager, ONLY: numout, lwp, numoni

   USE lbclnk         ! lateral boundary conditions (or mpp links)
#if defined _OPENACC
   USE lbclnk_gpu
#endif

   IMPLICIT NONE

   PRIVATE

   INTERFACE rmpT2F
      MODULE PROCEDURE rmpT2F_default, rmpT2F_spcmsk !, rmpT2Fnm
   END INTERFACE rmpT2F

   INTERFACE rmpF2T
      MODULE PROCEDURE rmpF2T_default !, rmpF2Tnm
   END INTERFACE rmpF2T

   INTERFACE rmpU2V
      MODULE PROCEDURE rmpU2V_default !, rmpU2Vnm
   END INTERFACE rmpU2V

   INTERFACE rmpV2U
      MODULE PROCEDURE rmpV2U_default !, rmpV2Unm
   END INTERFACE rmpV2U

   PUBLIC   rmpT2F
   PUBLIC   rmpT2U
   PUBLIC   rmpT2V
   PUBLIC   rmpF2T
   PUBLIC   rmpU2V
   PUBLIC   rmpV2U

   !! Better for GPU:
   PUBLIC do_rmpT2F
   PUBLIC do_rmpF2T

   PUBLIC do_rmpU2V
   PUBLIC do_rmpV2U

   PUBLIC do_rmpT2U
   PUBLIC do_rmpT2V

   PUBLIC do_Voce

   !!----------------------------------------------------------------------

CONTAINS


   FUNCTION rmpT2F_default( pxt,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2F_default
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpT2F_default(:,:) = 0._wp
      !!
      !%acc data copyin( pxt, xmskt, xmskf, e1e2t, r1_e1e2f) create(rmpT2F_default)
      !%acc parallel loop
      DO jj=Njs0-1, Nje0
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj+1
            i4 = ji+1 ; j4 = jj+1
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zt3 = pxt(i3,j3)*xmskt(i3,j3)
            zt4 = pxt(i4,j4)*xmskt(i4,j4)
            zfc = xmskf(ji,jj)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zt3 = zt3 * e1e2t(i3,j3)
               zt4 = zt4 * e1e2t(i4,j4)
               zfc = zfc * r1_e1e2f(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2) + xmskt(i3,j3) + xmskt(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpT2F_default(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F_default@icedyn_rhg_bbm', rmpT2F_default, 'F', 1._wp )
      END IF
      !!
   END FUNCTION rmpT2F_default


   FUNCTION rmpT2F_spcmsk( pxt, pmt, pmf, lbcl, lconserv )
      !! => provide required LSMs as arguments
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2F_spcmsk
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt, pmt, pmf
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpT2F_spcmsk(:,:) = 0._wp
      !!
      PRINT *, 'LOLO: `rmpT2F_spcmsk`: Njs0-1, Nje0, Nis0-1, Nie0=',Njs0-1, Nje0, Nis0-1, Nie0

      !%acc data copyin( pxt, pmt, pmf, e1e2t, r1_e1e2f) create(rmpT2F_spcmsk)
      !%acc parallel loop
      DO jj=Njs0-1, Nje0
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj+1
            i4 = ji+1 ; j4 = jj+1
            !!
            zt1 = pxt(i1,j1)*pmt(i1,j1)
            zt2 = pxt(i2,j2)*pmt(i2,j2)
            zt3 = pxt(i3,j3)*pmt(i3,j3)
            zt4 = pxt(i4,j4)*pmt(i4,j4)
            zfc = pmf(ji,jj)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zt3 = zt3 * e1e2t(i3,j3)
               zt4 = zt4 * e1e2t(i4,j4)
               zfc = zfc * r1_e1e2f(ji,jj)
            END IF
            !!
            zm = pmt(i1,j1) + pmt(i2,j2) + pmt(i3,j3) + pmt(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpT2F_spcmsk(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F_spcmsk@icedyn_rhg_bbm', rmpT2F_spcmsk, 'F', 1._wp )
      END IF
      !!
   END FUNCTION rmpT2F_spcmsk

   FUNCTION rmpT2Fnm( kh, pxt )
      !! The most simple version possible...
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2Fnm
      INTEGER,                      INTENT(in) :: kh  ! halo reach...
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      !!
      INTEGER  :: km, kp, ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4
      !
      !rmpT2Fnm(:,:) = 0._wp
      !
      km = MIN(nn_hls  ,kh)
      kp = MIN(nn_hls-1,kh)

      !%acc data copyin( pxt, xmskt, xmskf, e1e2t, r1_e1e2f) create(rmpT2Fnm)
      !%acc parallel loop
      DO jj=Njs0-km, Nje0+kp
         DO ji=Nis0-km, Nie0+kp
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj+1
            i4 = ji+1 ; j4 = jj+1
            !!
            zt1 = pxt(i1,j1)
            zt2 = pxt(i2,j2)
            zt3 = pxt(i3,j3)
            zt4 = pxt(i4,j4)
            !!
            rmpT2Fnm(ji,jj) = 0.25_wp * ( zt1 + zt2 + zt3 + zt4 )
            !!
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
      !!
   END FUNCTION rmpT2Fnm

   FUNCTION rmpF2Tnm( kh, pxt )
      !! The most simple version possible...
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpF2Tnm
      INTEGER,                      INTENT(in) :: kh  ! halo reach...
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      !!
      INTEGER  :: km, kp, ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4
      !
      !rmpF2Tnm(:,:) = 0._wp
      !
      km = MIN(nn_hls-1,kh)
      kp = MIN(nn_hls  ,kh)
      !
      !%acc data copyin( pxt, xmskt, xmskf, e1e2t, r1_e1e2f) create(rmpF2Tnm)
      !%acc parallel loop
      DO jj=Njs0-km, Nje0+kp
         DO ji=Nis0-km, Nie0+kp
            !!
            i1 = ji   ; j1 = jj
            i2 = ji-1 ; j2 = jj
            i3 = ji-1 ; j3 = jj-1
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pxt(i1,j1)
            zt2 = pxt(i2,j2)
            zt3 = pxt(i3,j3)
            zt4 = pxt(i4,j4)
            !!
            rmpF2Tnm(ji,jj) = 0.25_wp * ( zt1 + zt2 + zt3 + zt4 )
            !!
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
      !!
   END FUNCTION rmpF2Tnm



   FUNCTION rmpT2U( pxt,  pm0, lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2U
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      INTEGER,            OPTIONAL, INTENT(in) :: pm0
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, m0
      REAL(wp) :: zt1, zt2, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv

      lcnsrv = .FALSE.
      m0 = 0
      IF( PRESENT( pm0      ) ) m0     = pm0
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpT2U(:,:) = 0._wp
      !!
      DO jj=Njs0-1-m0, Nje0+1+m0
         DO ji=Nis0-1-m0, Nie0+m0
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zfc = umask(ji,jj,1)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zfc = zfc * r1_e1e2u(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpT2U(ji,jj) = ( zt1 + zt2 ) * zfc / zs
            !!
         END DO
      END DO
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2U@icedyn_rhg_bbm', rmpT2U, 'U', 1._wp )
      END IF
      !!
   END FUNCTION rmpT2U


   FUNCTION rmpT2V( pxt,  pm0, lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2V
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      INTEGER,            OPTIONAL, INTENT(in) :: pm0
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, m0
      REAL(wp) :: zt1, zt2, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv

      m0 = 0
      IF( PRESENT( pm0      ) ) m0     = pm0
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpT2V(:,:) = 0._wp
      !!
      DO jj=Njs0-1-m0, Nje0+m0
         DO ji=Nis0-1-m0, Nie0+1+m0
            !!
            i1 = ji ; j1 = jj
            i2 = ji ; j2 = jj+1
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zfc = vmask(ji,jj,1)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zfc = zfc * r1_e1e2v(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpT2V(ji,jj) = ( zt1 + zt2 ) * zfc / zs
            !!
         END DO
      END DO
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2V@icedyn_rhg_bbm', rmpT2V, 'V', 1._wp )
      END IF
      !!
   END FUNCTION rmpT2V







   FUNCTION rmpF2T_default( pxf, lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpF2T_default
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxf
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zf1, zf2, zf3, zf4, zs, zm, zz, zfc
      LOGICAL  :: lcnsrv
      !!
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpF2T_default(:,:) = 0._wp

      !%acc data copyin( pxf, xmskf, xmskt, e1e2f, r1_e1e2t) create(rmpF2T_default)
      !%acc parallel loop
      DO jj=Njs0, Nje0+1
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji-1 ; j2 = jj
            i3 = ji-1 ; j3 = jj-1
            i4 = ji   ; j4 = jj-1
            !!
            zf1 = pxf(i1,j1)*xmskf(i1,j1)
            zf2 = pxf(i2,j2)*xmskf(i2,j2)
            zf3 = pxf(i3,j3)*xmskf(i3,j3)
            zf4 = pxf(i4,j4)*xmskf(i4,j4)
            zfc = xmskt(ji,jj)
            IF( lcnsrv ) THEN
               zf1 = zf1 * e1e2f(i1,j1)
               zf2 = zf2 * e1e2f(i2,j2)
               zf3 = zf3 * e1e2f(i3,j3)
               zf4 = zf4 * e1e2f(i4,j4)
               zfc = zfc * r1_e1e2t(ji,jj)
            END IF
            !!
            zm = xmskf(i1,j1) + xmskf(i2,j2) + xmskf(i3,j3) + xmskf(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet F-point, `0` otherwize
            zfc = zfc * zz
            !
            zs = MAX( zm , 1.E-12_wp )
            !!
            rmpF2T_default(ji,jj) = ( zf1 + zf2 + zf3 + zf4 ) * zfc / zs
         END DO
      END DO
      !%acc end parallel loop
      !%acc end data
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpF2T_default@icedyn_rhg_bbm', rmpF2T_default, 'T', 1._wp )
      END IF
      !!
   END FUNCTION rmpF2T_default


   FUNCTION rmpU2V_default( pxu,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpU2V_default
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxu
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpU2V_default(:,:) = 0._wp
      !!
      DO jj=Njs0-1, Nje0
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji   ; j2 = jj+1
            i3 = ji-1 ; j3 = jj+1
            i4 = ji-1 ; j4 = jj
            !!
            zt1 = pxu(i1,j1)*umask(i1,j1,1)
            zt2 = pxu(i2,j2)*umask(i2,j2,1)
            zt3 = pxu(i3,j3)*umask(i3,j3,1)
            zt4 = pxu(i4,j4)*umask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2u(i1,j1)
               zt2 = zt2 * e1e2u(i2,j2)
               zt3 = zt3 * e1e2u(i3,j3)
               zt4 = zt4 * e1e2u(i4,j4)
               zfc = zfc * r1_e1e2v(ji,jj)
            END IF
            !!
            zm = umask(i1,j1,1) + umask(i2,j2,1) + umask(i3,j3,1) + umask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpU2V_default(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpU2V_default@icedyn_rhg_bbm', rmpU2V_default, 'V', 1._wp )
      END IF
      !!
   END FUNCTION rmpU2V_default


   SUBROUTINE do_rmpU2V( pxu,  pxv,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxu
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxv
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !
      !$acc parallel loop collapse(2) present(pxu, pxv, umask, e1e2u, r1_e1e2v)
      DO jj=Njs0-1, Nje0
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji   ; j2 = jj+1
            i3 = ji-1 ; j3 = jj+1
            i4 = ji-1 ; j4 = jj
            !!
            zt1 = pxu(i1,j1)*umask(i1,j1,1)
            zt2 = pxu(i2,j2)*umask(i2,j2,1)
            zt3 = pxu(i3,j3)*umask(i3,j3,1)
            zt4 = pxu(i4,j4)*umask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2u(i1,j1)
               zt2 = zt2 * e1e2u(i2,j2)
               zt3 = zt3 * e1e2u(i3,j3)
               zt4 = zt4 * e1e2u(i4,j4)
               zfc = zfc * r1_e1e2v(ji,jj)
            END IF
            !!
            zm = umask(i1,j1,1) + umask(i2,j2,1) + umask(i3,j3,1) + umask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pxv(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'do_rmpU2V@icedyn_rhg_bbm', pxv, 'V', 1._wp )
      END IF
      !!
   END SUBROUTINE do_rmpU2V


   FUNCTION rmpV2U_default( pxv,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpV2U_default
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxv
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpV2U_default(:,:) = 0._wp
      !!
      DO jj=Njs0, Nje0+1
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji+1 ; j1 = jj-1
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pxv(i1,j1)*vmask(i1,j1,1)
            zt2 = pxv(i2,j2)*vmask(i2,j2,1)
            zt3 = pxv(i3,j3)*vmask(i3,j3,1)
            zt4 = pxv(i4,j4)*vmask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2v(i1,j1)
               zt2 = zt2 * e1e2v(i2,j2)
               zt3 = zt3 * e1e2v(i3,j3)
               zt4 = zt4 * e1e2v(i4,j4)
               zfc = zfc * r1_e1e2u(ji,jj)
            END IF
            !!
            zm = vmask(i1,j1,1) + vmask(i2,j2,1) + vmask(i3,j3,1) + vmask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpV2U_default(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpV2U_default@icedyn_rhg_bbm', rmpV2U_default, 'U', 1._wp )
      END IF
      !!
   END FUNCTION rmpV2U_default

   SUBROUTINE do_rmpV2U( pxv, pxu,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxv
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxu
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !
      !$acc parallel loop collapse(2) present(pxv, pxu, vmask, e1e2v, r1_e1e2u)
      DO jj=Njs0, Nje0+1
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji+1 ; j1 = jj-1
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pxv(i1,j1)*vmask(i1,j1,1)
            zt2 = pxv(i2,j2)*vmask(i2,j2,1)
            zt3 = pxv(i3,j3)*vmask(i3,j3,1)
            zt4 = pxv(i4,j4)*vmask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2v(i1,j1)
               zt2 = zt2 * e1e2v(i2,j2)
               zt3 = zt3 * e1e2v(i3,j3)
               zt4 = zt4 * e1e2v(i4,j4)
               zfc = zfc * r1_e1e2u(ji,jj)
            END IF
            !!
            zm = vmask(i1,j1,1) + vmask(i2,j2,1) + vmask(i3,j3,1) + vmask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pxu(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop
      !
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'do_rmpV2U@icedyn_rhg_bbm', pxu, 'U', 1._wp )
      END IF
      !
   END SUBROUTINE do_rmpV2U




   FUNCTION rmpU2Vnm( pxu )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpU2Vnm
      !INTEGER,                      INTENT(in) :: kh  ! halo reach...
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxu
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4
      !
      !kp = MIN(nn_hls  ,kh)
      !km = MIN(nn_hls-1,kh)
      !
      !DO jj=Njs0-kp, Nje0+km
      !   DO ji=Nis0+km, Nie0+kp
      DO jj=Njs0-1, Nje0
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji   ; j2 = jj+1
            i3 = ji-1 ; j3 = jj+1
            i4 = ji-1 ; j4 = jj
            !!
            zt1 = pxu(i1,j1)
            zt2 = pxu(i2,j2)
            zt3 = pxu(i3,j3)
            zt4 = pxu(i4,j4)
            !!
            rmpU2Vnm(ji,jj) = 0.25_wp * ( zt1 + zt2 + zt3 + zt4 )
            !!
         END DO
      END DO
      !!
   END FUNCTION rmpU2Vnm

   FUNCTION rmpV2Unm( pxv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpV2Unm
      !INTEGER,                      INTENT(in) :: kh  ! halo reach...
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4
      !
      !kp = MIN(nn_hls  ,kh)
      !km = MIN(nn_hls-1,kh)
      !
      !DO jj=Njs0-km, Nje0+kp
      !   DO ji=Nis0-kp, Nie0+km
      DO jj=Njs0, Nje0+1
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji+1 ; j1 = jj-1
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pxv(i1,j1)
            zt2 = pxv(i2,j2)
            zt3 = pxv(i3,j3)
            zt4 = pxv(i4,j4)
            !!
            rmpV2Unm(ji,jj) = 0.25_wp * ( zt1 + zt2 + zt3 + zt4 )
            !!
         END DO
      END DO
      !!
   END FUNCTION rmpV2Unm



   SUBROUTINE do_Voce( pum,  pvm,  pVoce )
      !!---------------------------------------------------------
      !! => prepare `V_oce(:,:,1:4) out of `ssu_m` & `ssv_m` !
      !!
      REAL(wp), DIMENSION(jpi,jpj),    INTENT(in) :: pum, pvm
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(out) :: pVoce
      !!---------------------------------------------------------
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      !!---------------------------------------------------------
      !$acc data present( pum, pvm, pVoce, umask, vmask )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pVoce(ji,jj,1) = pum(ji,jj)
            pVoce(ji,jj,2) = pvm(ji,jj)
            pVoce(ji,jj,3) = 0._wp
            pVoce(ji,jj,4) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !! U 2 V:
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji   ; j2 = jj+1
            i3 = ji-1 ; j3 = jj+1
            i4 = ji-1 ; j4 = jj
            !!
            zt1 = pum(i1,j1)*umask(i1,j1,1)
            zt2 = pum(i2,j2)*umask(i2,j2,1)
            zt3 = pum(i3,j3)*umask(i3,j3,1)
            zt4 = pum(i4,j4)*umask(i4,j4,1)
            zfc = 1._wp
            !!
            zm = umask(i1,j1,1) + umask(i2,j2,1) + umask(i3,j3,1) + umask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pVoce(ji,jj,3) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop

      !! V 2 U:
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0+1
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji+1 ; j1 = jj-1
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pvm(i1,j1)*vmask(i1,j1,1)
            zt2 = pvm(i2,j2)*vmask(i2,j2,1)
            zt3 = pvm(i3,j3)*vmask(i3,j3,1)
            zt4 = pvm(i4,j4)*vmask(i4,j4,1)
            zfc = 1._wp
            !!
            zm = vmask(i1,j1,1) + vmask(i2,j2,1) + vmask(i3,j3,1) + vmask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pVoce(ji,jj,4) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop
      !
# if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'do_Voce', pVoce )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'do_Voce', pVoce )
# else
      CALL lbc_lnk( 'do_Voce', pVoce(:,:,1),'U',-1._wp, pVoce(:,:,2),'V',-1._wp, pVoce(:,:,3),'V',-1._wp, pVoce(:,:,4),'U',-1._wp )
# endif
      !
      !$acc end data
      !
   END SUBROUTINE do_Voce





   SUBROUTINE do_rmpT2F( pxt, pxf,  lbcl, lconserv )
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxf
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl, lconserv
      !!-------------------------------------------------------------------
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!-------------------------------------------------------------------
      !$acc data present( pxt, pxf, xmskt, xmskf, e1e2t, r1_e1e2f )
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pxf(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0
         DO ji=Nis0-1, Nie0
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj+1
            i4 = ji+1 ; j4 = jj+1
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zt3 = pxt(i3,j3)*xmskt(i3,j3)
            zt4 = pxt(i4,j4)*xmskt(i4,j4)
            zfc = xmskf(ji,jj)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zt3 = zt3 * e1e2t(i3,j3)
               zt4 = zt4 * e1e2t(i4,j4)
               zfc = zfc * r1_e1e2f(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2) + xmskt(i3,j3) + xmskt(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pxf(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop

# if ! defined _OPENACC
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'do_rmpT2F@icedyn_rhg_bbm', pxf, 'F', 1._wp )
      END IF
# endif
      !$acc end data
      !!
   END SUBROUTINE do_rmpT2F

   SUBROUTINE do_rmpF2T( pxf, pxt, lbcl, lconserv )
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl, lconserv
      !!-------------------------------------------------------------------
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zf1, zf2, zf3, zf4, zs, zm, zz, zfc
      LOGICAL  :: lcnsrv
      !!-------------------------------------------------------------------
      !$acc data present( pxf, pxt, xmskt, xmskf, e1e2f, r1_e1e2t )
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pxt(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0+1
         DO ji=Nis0, Nie0+1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji-1 ; j2 = jj
            i3 = ji-1 ; j3 = jj-1
            i4 = ji   ; j4 = jj-1
            !!
            zf1 = pxf(i1,j1)*xmskf(i1,j1)
            zf2 = pxf(i2,j2)*xmskf(i2,j2)
            zf3 = pxf(i3,j3)*xmskf(i3,j3)
            zf4 = pxf(i4,j4)*xmskf(i4,j4)
            zfc = xmskt(ji,jj)
            IF( lcnsrv ) THEN
               zf1 = zf1 * e1e2f(i1,j1)
               zf2 = zf2 * e1e2f(i2,j2)
               zf3 = zf3 * e1e2f(i3,j3)
               zf4 = zf4 * e1e2f(i4,j4)
               zfc = zfc * r1_e1e2t(ji,jj)
            END IF
            !!
            zm = xmskf(i1,j1) + xmskf(i2,j2) + xmskf(i3,j3) + xmskf(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet F-point, `0` otherwize
            zfc = zfc * zz
            !
            zs = MAX( zm , 1.E-12_wp )
            !!
            pxt(ji,jj) = ( zf1 + zf2 + zf3 + zf4 ) * zfc / zs
         END DO
      END DO
      !$acc end parallel loop
      !!
# if ! defined _OPENACC
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'pxt@icedyn_rhg_bbm', pxt, 'T', 1._wp )
      END IF
# endif
      !$acc end data
      !!
   END SUBROUTINE do_rmpF2T




   SUBROUTINE do_rmpT2U( pxt, pxu,  lbcl, lconserv )
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxu
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl, lconserv
      !!-------------------------------------------------------------------
      INTEGER  :: ji, jj, i1, j1, i2, j2
      REAL(wp) :: zt1, zt2, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!-------------------------------------------------------------------
      !$acc data present( umask, xmskt, e1e2t, r1_e1e2u ) pcopyin( pxt ) pcopyout( pxu )
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pxu(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls-1
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zfc = umask(ji,jj,1)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zfc = zfc * r1_e1e2u(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pxu(ji,jj) = ( zt1 + zt2 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop
      !!
# if ! defined _OPENACC
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'do_rmpT2U', pxu, 'U', 1._wp )
      END IF
# endif
      !$acc end data
      !!
   END SUBROUTINE do_rmpT2U


   SUBROUTINE do_rmpT2V( pxt, pxv,  lbcl, lconserv )
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pxt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pxv
      LOGICAL,            OPTIONAL, INTENT(in)  :: lbcl, lconserv
      !!-------------------------------------------------------------------
      INTEGER  :: ji, jj, i1, j1, i2, j2
      REAL(wp) :: zt1, zt2, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv
      !!-------------------------------------------------------------------
      !$acc data present( vmask, xmskt, e1e2t, r1_e1e2v ) pcopyin( pxt ) pcopyout( pxv )
      lcnsrv = .FALSE.
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pxv(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls-1
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !!
            i1 = ji ; j1 = jj
            i2 = ji ; j2 = jj+1
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zfc = vmask(ji,jj,1)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zfc = zfc * r1_e1e2v(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            pxv(ji,jj) = ( zt1 + zt2 ) * zfc / zs
            !!
         END DO
      END DO
      !$acc end parallel loop
      !!
# if ! defined _OPENACC
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'do_rmpT2V', pxv, 'V', 1._wp )
      END IF
# endif
      !$acc end data
      !!
   END SUBROUTINE do_rmpT2V



   !!======================================================================
END MODULE remap_classic
