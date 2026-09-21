MODULE lib_fortran
   !!======================================================================
   !!                       ***  MODULE  lib_fortran  ***
   !! Fortran utilities:  includes some low levels fortran functionality
   !!======================================================================
   !! History :  3.2  !  2010-05  (M. Dunphy, R. Benshila)  Original code
   !!            3.4  !  2013-06  (C. Rousset)  add glob_min, glob_max
   !!                                           + 3d dim. of input is fexible (jpk, jpl...)
   !!            4.0  !  2016-06  (T. Lovato)  double precision global sum by default
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   local_2Dmin/max/sum : generic interface for 2D local min/max/sum
   !!                         masked with tmask_i, for 2D or 3D arrays
   !!                         it works only for T points
   !!   local_3Dmin/max/sum : generic interface for 3D local min/max/sum
   !!                         masked with tmask, for 3D or 4D arrays
   !!                         it works only for T points
   !!   glob_2Dmin/max/sum  : generic interface for 2D global min/max/sum
   !!                         masked with tmask_i, for 2D or 3D arrays
   !!                         it works only for T points
   !!                         calls local_2Dmin/max/sum + mpp_min/max/sum
   !!   glob_3Dmin/max/sum  : generic interface for 3D global min/max/sum
   !!                         masked with tmask, for 3D or 4D arrays
   !!                         it works only for T points
   !!                         calls local_3Dmin/max/sum + mpp_min/max/sum
   !!   SIGN        : generic interface for SIGN to overwrite f95 behaviour
   !!                 of intrinsinc sign function
   !!----------------------------------------------------------------------
   USE par_oce         ! Ocean parameter
   USE dom_oce         ! ocean domain
   USE domutl, ONLY : is_tile
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions

   IMPLICIT NONE
   PRIVATE

   PUBLIC   local_2Dmin, local_2Dmax, local_2Dsum
   PUBLIC    glob_2Dmin,  glob_2Dmax,  glob_2Dsum
   PUBLIC   local_3Dmin, local_3Dmax, local_3Dsum
   PUBLIC    glob_3Dmin,  glob_3Dmax,  glob_3Dsum
   PUBLIC   sum3x3        ! used in trcrad, do a sum over 3x3 boxes
#if defined key_nosignedzero
   PUBLIC SIGN
#endif

   INTERFACE local_2Dmin
      MODULE PROCEDURE local_2Dmin_2Din_sp, local_2Dmin_2Din_dp
      MODULE PROCEDURE local_2Dmin_3Din_sp, local_2Dmin_3Din_dp
   END INTERFACE local_2Dmin
   INTERFACE local_3Dmin
      MODULE PROCEDURE local_3Dmin_3Din_sp, local_3Dmin_3Din_dp
      MODULE PROCEDURE local_3Dmin_4Din_sp, local_3Dmin_4Din_dp
   END INTERFACE local_3Dmin
   INTERFACE glob_2Dmin
      MODULE PROCEDURE glob_2Dmin_2Din_sp, glob_2Dmin_2Din_dp
      MODULE PROCEDURE glob_2Dmin_3Din_sp, glob_2Dmin_3Din_dp
   END INTERFACE glob_2Dmin
   INTERFACE glob_3Dmin
      MODULE PROCEDURE glob_3Dmin_3Din_sp, glob_3Dmin_3Din_dp
      MODULE PROCEDURE glob_3Dmin_4Din_sp, glob_3Dmin_4Din_dp
   END INTERFACE glob_3Dmin

   INTERFACE local_2Dmax
      MODULE PROCEDURE local_2Dmax_2Din_sp, local_2Dmax_2Din_dp
      MODULE PROCEDURE local_2Dmax_3Din_sp, local_2Dmax_3Din_dp
   END INTERFACE local_2Dmax
   INTERFACE local_3Dmax
      MODULE PROCEDURE local_3Dmax_3Din_sp, local_3Dmax_3Din_dp
      MODULE PROCEDURE local_3Dmax_4Din_sp, local_3Dmax_4Din_dp
   END INTERFACE local_3Dmax
   INTERFACE glob_2Dmax
      MODULE PROCEDURE glob_2Dmax_2Din_sp, glob_2Dmax_2Din_dp
      MODULE PROCEDURE glob_2Dmax_3Din_sp, glob_2Dmax_3Din_dp
   END INTERFACE glob_2Dmax
   INTERFACE glob_3Dmax
      MODULE PROCEDURE glob_3Dmax_3Din_sp, glob_3Dmax_3Din_dp
      MODULE PROCEDURE glob_3Dmax_4Din_sp, glob_3Dmax_4Din_dp
   END INTERFACE glob_3Dmax

   INTERFACE local_2Dsum
      MODULE PROCEDURE local_2Dsum_2Din_sp, local_2Dsum_2Din_dp
      MODULE PROCEDURE local_2Dsum_3Din_sp, local_2Dsum_3Din_dp
   END INTERFACE local_2Dsum
   INTERFACE local_3Dsum
      MODULE PROCEDURE local_3Dsum_3Din_sp, local_3Dsum_3Din_dp
      MODULE PROCEDURE local_3Dsum_4Din_sp, local_3Dsum_4Din_dp
   END INTERFACE local_3Dsum
   INTERFACE glob_2Dsum
      MODULE PROCEDURE glob_2Dsum_2Din_sp, glob_2Dsum_2Din_dp
      MODULE PROCEDURE glob_2Dsum_3Din_sp, glob_2Dsum_3Din_dp
   END INTERFACE glob_2Dsum
   INTERFACE glob_3Dsum
      MODULE PROCEDURE glob_3Dsum_3Din_sp, glob_3Dsum_3Din_dp
      MODULE PROCEDURE glob_3Dsum_4Din_sp, glob_3Dsum_4Din_dp
   END INTERFACE glob_3Dsum

   INTERFACE sum3x3
      MODULE PROCEDURE sum3x3_2d, sum3x3_3d
   END INTERFACE sum3x3

#if defined key_nosignedzero
   INTERFACE SIGN
      MODULE PROCEDURE SIGN_POSZERO
   END INTERFACE SIGN
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   /**/
#  define LOCAL_GLOBAL
   /**/
#   define OPERMIN
#    define REALSP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#        include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef XDOPER
#    undef  REALSP_TYPE
#    define REALDP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#    undef  REALDP_TYPE
#   undef OPERMIN
   /**/
#   define OPERMAX
#    define REALSP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#        include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef XDOPER
#    undef  REALSP_TYPE
#    define REALDP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#    undef  REALDP_TYPE
#   undef OPERMAX
   /**/
#   define OPERSUM
#    define REALSP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#        include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef XDOPER
#    undef  REALSP_TYPE
#    define REALDP_TYPE
#     define XDOPER 2
#      define XDIN 2
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#     define XDOPER 3
#      define XDIN 3
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#      define XDIN 4
#       include "lib_fortran_generic.h90"
#      undef  XDIN
#     undef  XDOPER
#    undef  REALDP_TYPE
#   undef OPERSUM
   /**/
#  undef  LOCAL_GLOBAL


   SUBROUTINE sum3x3_2d( p2d )
      !!-----------------------------------------------------------------------
      !!                  ***  routine sum3x3_2d  ***
      !!
      !! ** Purpose : sum over 3x3 boxes
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION (:,:), INTENT(inout) ::   p2d
      !
      INTEGER ::   ji, ji2, jj, jj2     ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( SIZE(p2d,1) /= jpi ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_2d, the first dimension is not equal to jpi' )
      IF( SIZE(p2d,2) /= jpj ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_2d, the second dimension is not equal to jpj' )
      !
      ! work over the whole domain (guarantees all internal cells are set when nn_hls=2)
      !
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF( MOD(mig(ji,nn_hls), 3) == MOD(nn_hls, 3) .AND.   &              ! 1st bottom left corner always at (Nis0-1, Njs0-1)
               & MOD(mjg(jj,nn_hls), 3) == MOD(nn_hls, 3)         ) THEN         ! bottom left corner of a 3x3 box
               ji2 = MIN(mig(ji,nn_hls)+2, jpiglo) - nimpp + 1                  ! right position of the box
               jj2 = MIN(mjg(jj,nn_hls)+2, jpjglo) - njmpp + 1                  ! upper position of the box
               IF( ji2 <= jpi .AND. jj2 <= jpj ) THEN                           ! the box is fully included in the local mpi domain
                  p2d(ji:ji2,jj:jj2) = SUM(p2d(ji:ji2,jj:jj2))
               ENDIF
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'lib_fortran', p2d, 'T', 1.0_wp )

   END SUBROUTINE sum3x3_2d

   SUBROUTINE sum3x3_3d( p3d )
      !!-----------------------------------------------------------------------
      !!                  ***  routine sum3x3_3d  ***
      !!
      !! ** Purpose : sum over 3x3 boxes
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   p3d
      !
      INTEGER ::   ji, ji2, jj, jj2, jn     ! dummy loop indices
      INTEGER ::   ipn                      ! Third dimension size
      !!----------------------------------------------------------------------
      !
      IF( SIZE(p3d,1) /= jpi ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_3d, the first dimension is not equal to jpi' )
      IF( SIZE(p3d,2) /= jpj ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_3d, the second dimension is not equal to jpj' )
      ipn = SIZE(p3d,3)
      !
      DO jn = 1, ipn
         !
         ! work over the whole domain (guarantees all internal cells are set when nn_hls=2)
         !
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               IF( MOD(mig(ji,nn_hls), 3) == MOD(nn_hls, 3) .AND.   &              ! 1st bottom left corner always at (Nis0-1, Njs0-1)
                  & MOD(mjg(jj,nn_hls), 3) == MOD(nn_hls, 3)         ) THEN         ! bottom left corner of a 3x3 box
                  ji2 = MIN(mig(ji,nn_hls)+2, jpiglo) - nimpp + 1                  ! right position of the box
                  jj2 = MIN(mjg(jj,nn_hls)+2, jpjglo) - njmpp + 1                  ! upper position of the box
                  IF( ji2 <= jpi .AND. jj2 <= jpj ) THEN                           ! the box is fully included in the local mpi domain
                     p3d(ji:ji2,jj:jj2,jn) = SUM(p3d(ji:ji2,jj:jj2,jn))
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'lib_fortran', p3d, 'T', 1.0_wp )

   END SUBROUTINE sum3x3_3d


#if defined key_nosignedzero
   !!----------------------------------------------------------------------
   !!   'key_nosignedzero'                                         F90 SIGN
   !!----------------------------------------------------------------------

   ELEMENTAL FUNCTION SIGN_POSZERO( a, b )
      !!-----------------------------------------------------------------------
      !!                    ***  FUNCTION SIGN_POSZERO  ***
      !!
      !! ** Purpose : substitute for the intrinsic SIGN function; it ignores
      !!              the sign of signed zeros that may be present as the
      !!              value of argument b; its dummy-argument names are exempt
      !!              from compliance with the NEMO coding convention in order
      !!              to enable the conformity between its interface and the
      !!              interface of the intrinsic SIGN function
      !!
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   a, b           ! input
      REAL(wp)             ::   SIGN_POSZERO   ! result
      !!-----------------------------------------------------------------------
      IF( b >= 0._wp ) THEN
         SIGN_POSZERO =  ABS(a)
      ELSE
         SIGN_POSZERO = -ABS(a)
      ENDIF
   END FUNCTION SIGN_POSZERO

#endif

   !!======================================================================
END MODULE lib_fortran
