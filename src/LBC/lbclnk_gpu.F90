MODULE lbclnk_gpu
   !!*******************************************************************************************************************
   !! INTENDED TO BE USED ONLY WHEN NO MPI DOMAIN DECOMPOSITION IS USED == WHEN NANUQ IS RUNNING ON 1 CORE (+ 1 GPU ) !!
   !!  + they do not support North-Pole folding of ORCA* grids yet...
   !!    => so vector components are just treated as other fields...
   !!*******************************************************************************************************************

   USE dom_oce        ! ocean space and time domain
   USE lib_mpp        ! distributed memory computing library
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk_EW_gpu
      MODULE PROCEDURE lbc_lnk_EW_1f2d_r8, lbc_lnk_EW_2f2d_r8, lbc_lnk_EW_3f2d_r8, lbc_lnk_EW_4f2d_r8, lbc_lnk_EW_1f3d4_r8, lbc_lnk_EW_2f3d3_r8, lbc_lnk_EW_2f3d3_2f2d_r8
   END INTERFACE lbc_lnk_EW_gpu

   INTERFACE lbc_lnk_NS_gpu
      MODULE PROCEDURE lbc_lnk_NS_1f2d_r8, lbc_lnk_NS_2f2d_r8, lbc_lnk_NS_3f2d_r8, lbc_lnk_NS_4f2d_r8, lbc_lnk_NS_1f3d4_r8, lbc_lnk_NS_2f3d3_r8, lbc_lnk_NS_2f3d3_2f2d_r8
   END INTERFACE lbc_lnk_NS_gpu

   PUBLIC   lbc_lnk_EW_gpu
   PUBLIC   lbc_lnk_NS_gpu

CONTAINS

   !!####################################################################################################################
   !!     Example for `nn_hls=3`:
   !!
   !!    ji:        1                                                                                            jpi
   !!               |xxxx|xxxx|xxxx|----|----|----|----|----|----|----   ...   |----|----|----|----|xxxx|xxxx|xxxx|
   !!               1    2    3    4    5    6    7                                               N-3  N-2  N-1   N
   !!                             Nis0                                                           Nie0
   !!
   !!                             jpi
   !! ... |----|----|xxxx|xxxx|xxxx|
   !!              N-3  N-2  N-1   N
   !!              Nie0
   !!####################################################################################################################


   !!##################################################################################################
   !!            * 1 2D array (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_1f2d_r8( cdname, pf2d )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d
      !!----------------------------------------------------------------------
      INTEGER    :: jj, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO jj = 1, jpj
            pf2d(    jh  ,jj) = pf2d(jpi-jc  ,jj)
            pf2d(jpi-jh+1,jj) = pf2d(    jc+1,jj)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_1f2d_r8
   !
   SUBROUTINE lbc_lnk_NS_1f2d_r8( cdname, pf2d )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO ji = 1, jpi
            pf2d(ji,    jh  ) = pf2d(ji,jpj-jc  )
            pf2d(ji,jpj-jh+1) = pf2d(ji,    jc+1)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_1f2d_r8
   !!##################################################################################################

   !!##################################################################################################
   !!            * 2 2D arrays (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_2f2d_r8( cdname, pf2d1, pf2d2 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2
      !!----------------------------------------------------------------------
      INTEGER    :: jj, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO jj = 1, jpj
            pf2d1(    jh  ,jj) = pf2d1(jpi-jc  ,jj)
            pf2d1(jpi-jh+1,jj) = pf2d1(    jc+1,jj)
            pf2d2(    jh  ,jj) = pf2d2(jpi-jc  ,jj)
            pf2d2(jpi-jh+1,jj) = pf2d2(    jc+1,jj)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_2f2d_r8
   !
   SUBROUTINE lbc_lnk_NS_2f2d_r8( cdname, pf2d1, pf2d2 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO ji = 1, jpi
            pf2d1(ji,    jh  ) = pf2d1(ji,jpj-jc  )
            pf2d1(ji,jpj-jh+1) = pf2d1(ji,    jc+1)
            pf2d2(ji,    jh  ) = pf2d2(ji,jpj-jc  )
            pf2d2(ji,jpj-jh+1) = pf2d2(ji,    jc+1)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_2f2d_r8
   !!##################################################################################################

   !!##################################################################################################
   !!            * 3 2D arrays (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_3f2d_r8( cdname, pf2d1, pf2d2, pf2d3 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2, pf2d3
      !!----------------------------------------------------------------------
      INTEGER    :: jj, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2, pf2d3)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO jj = 1, jpj
            pf2d1(    jh  ,jj) = pf2d1(jpi-jc  ,jj)
            pf2d1(jpi-jh+1,jj) = pf2d1(    jc+1,jj)
            pf2d2(    jh  ,jj) = pf2d2(jpi-jc  ,jj)
            pf2d2(jpi-jh+1,jj) = pf2d2(    jc+1,jj)
            pf2d3(    jh  ,jj) = pf2d3(jpi-jc  ,jj)
            pf2d3(jpi-jh+1,jj) = pf2d3(    jc+1,jj)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_3f2d_r8
   !
   SUBROUTINE lbc_lnk_NS_3f2d_r8( cdname, pf2d1, pf2d2, pf2d3 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2, pf2d3
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2, pf2d3)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO ji = 1, jpi
            pf2d1(ji,    jh  ) = pf2d1(ji,jpj-jc  )
            pf2d1(ji,jpj-jh+1) = pf2d1(ji,    jc+1)
            pf2d2(ji,    jh  ) = pf2d2(ji,jpj-jc  )
            pf2d2(ji,jpj-jh+1) = pf2d2(ji,    jc+1)
            pf2d3(ji,    jh  ) = pf2d3(ji,jpj-jc  )
            pf2d3(ji,jpj-jh+1) = pf2d3(ji,    jc+1)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_3f2d_r8
   !!##################################################################################################

   !!##################################################################################################
   !!            * 4 2D arrays (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_4f2d_r8( cdname, pf2d1, pf2d2, pf2d3, pf2d4 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2, pf2d3, pf2d4
      !!----------------------------------------------------------------------
      INTEGER    :: jj, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2, pf2d3, pf2d4)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO jj = 1, jpj
            pf2d1(    jh  ,jj) = pf2d1(jpi-jc  ,jj)
            pf2d1(jpi-jh+1,jj) = pf2d1(    jc+1,jj)
            pf2d2(    jh  ,jj) = pf2d2(jpi-jc  ,jj)
            pf2d2(jpi-jh+1,jj) = pf2d2(    jc+1,jj)
            pf2d3(    jh  ,jj) = pf2d3(jpi-jc  ,jj)
            pf2d3(jpi-jh+1,jj) = pf2d3(    jc+1,jj)
            pf2d4(    jh  ,jj) = pf2d4(jpi-jc  ,jj)
            pf2d4(jpi-jh+1,jj) = pf2d4(    jc+1,jj)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_4f2d_r8
   !
   SUBROUTINE lbc_lnk_NS_4f2d_r8( cdname, pf2d1, pf2d2, pf2d3, pf2d4 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,        INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf2d1, pf2d2, pf2d3, pf2d4
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pf2d1, pf2d2, pf2d3, pf2d4)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop
         DO ji = 1, jpi
            pf2d1(ji,    jh  ) = pf2d1(ji,jpj-jc  )
            pf2d1(ji,jpj-jh+1) = pf2d1(ji,    jc+1)
            pf2d2(ji,    jh  ) = pf2d2(ji,jpj-jc  )
            pf2d2(ji,jpj-jh+1) = pf2d2(ji,    jc+1)
            pf2d3(ji,    jh  ) = pf2d3(ji,jpj-jc  )
            pf2d3(ji,jpj-jh+1) = pf2d3(ji,    jc+1)
            pf2d4(ji,    jh  ) = pf2d4(ji,jpj-jc  )
            pf2d4(ji,jpj-jh+1) = pf2d4(ji,    jc+1)
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_4f2d_r8
   !!##################################################################################################



   !!##################################################################################################
   !!            * 1 3D array (3rd dim len 4) (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_1f3d4_r8( cdname, pv4 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(inout) :: pv4
      !!----------------------------------------------------------------------
      INTEGER  :: jj, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pv4)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 4
            DO jj = 1, jpj
               pv4(    jh  ,jj,jk) = pv4(jpi-jc  ,jj,jk)
               pv4(jpi-jh+1,jj,jk) = pv4(    jc+1,jj,jk)
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_1f3d4_r8

   SUBROUTINE lbc_lnk_NS_1f3d4_r8( cdname, pv4 )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(inout) :: pv4
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pv4)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 4
            DO ji = 1, jpi
               pv4(ji,    jh  ,jk) = pv4(ji,jpj-jc  ,jk)
               pv4(ji,jpj-jh+1,jk) = pv4(ji,    jc+1,jk)
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_1f3d4_r8
   !!##################################################################################################


   !!##################################################################################################
   !!            * 2 3D arrays (3rd dim len 3) (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_2f3d3_r8( cdname, pSt, psf )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSt, psf
      !!----------------------------------------------------------------------
      INTEGER  :: jj, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pSt, psf)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 3
            DO jj = 1, jpj
               pSt(    jh  ,jj,jk) = pSt(jpi-jc  ,jj,jk)
               pSt(jpi-jh+1,jj,jk) = pSt(    jc+1,jj,jk)
               psf(    jh  ,jj,jk) = psf(jpi-jc  ,jj,jk)
               psf(jpi-jh+1,jj,jk) = psf(    jc+1,jj,jk)
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_2f3d3_r8

   SUBROUTINE lbc_lnk_NS_2f3d3_r8( cdname, pSt, pSf )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSt, pSf
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pSt, pSf)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 3
            DO ji = 1, jpi
               pSt(ji,    jh  ,jk) = pSt(ji,jpj-jc  ,jk)
               pSt(ji,jpj-jh+1,jk) = pSt(ji,    jc+1,jk)
               pSf(ji,    jh  ,jk) = pSf(ji,jpj-jc  ,jk)
               pSf(ji,jpj-jh+1,jk) = pSf(ji,    jc+1,jk)
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_2f3d3_r8
   !!##################################################################################################


   !!##################################################################################################
   !!            * 2 3D arrays (3rd dim len 3) (real 8)
   !!##################################################################################################
   SUBROUTINE lbc_lnk_EW_2f3d3_2f2d_r8( cdname, pSt, pSf, pdt, pdf )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSt, pSf
      REAL(wp), DIMENSION(jpi,jpj)  , INTENT(inout) :: pdt, pdf
      !!----------------------------------------------------------------------
      INTEGER  :: jj, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pSt, pSf, pdt, pdf)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 3
            DO jj = 1, jpj
               pSt(    jh  ,jj,jk) = pSt(jpi-jc  ,jj,jk)
               pSt(jpi-jh+1,jj,jk) = pSt(    jc+1,jj,jk)
               pSf(    jh  ,jj,jk) = pSf(jpi-jc  ,jj,jk)
               pSf(jpi-jh+1,jj,jk) = pSf(    jc+1,jj,jk)
               pdt(    jh  ,jj   ) = pdt(jpi-jc  ,jj   )
               pdt(jpi-jh+1,jj   ) = pdt(    jc+1,jj   )
               pdf(    jh  ,jj   ) = pdf(jpi-jc  ,jj   )
               pdf(jpi-jh+1,jj   ) = pdf(    jc+1,jj   )
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_EW_2f3d3_2f2d_r8

   SUBROUTINE lbc_lnk_NS_2f3d3_2f2d_r8( cdname, pSt, pSf, pdt, pdf )
      !!---------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(len=*)     ,          INTENT(in   ) :: cdname    ! name of the calling subroutine
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pSt, pSf
      REAL(wp), DIMENSION(jpi,jpj)  , INTENT(inout) :: pdt, pdf
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jk, jh, jc
      !!----------------------------------------------------------------------
      !$acc data present(pSt, pSf, pdt, pdf)
      !$acc loop seq
      DO jh = 1, nn_hls
         jc = nn_hls-jh+1
         !$acc parallel loop collapse(2)
         DO jk = 1, 3
            DO ji = 1, jpi
               pSt(ji,    jh  ,jk) = pSt(ji,jpj-jc  ,jk)
               pSt(ji,jpj-jh+1,jk) = pSt(ji,    jc+1,jk)
               pSf(ji,    jh  ,jk) = pSf(ji,jpj-jc  ,jk)
               pSf(ji,jpj-jh+1,jk) = pSf(ji,    jc+1,jk)
               pdt(ji,    jh     ) = pdt(ji,jpj-jc     )
               pdt(ji,jpj-jh+1   ) = pdt(ji,    jc+1   )
               pdf(ji,    jh     ) = pdf(ji,jpj-jc     )
               pdf(ji,jpj-jh+1   ) = pdf(ji,    jc+1   )
            END DO
         END DO
         !$acc end parallel loop
      END DO
      !$acc end data
   END SUBROUTINE lbc_lnk_NS_2f3d3_2f2d_r8
   !!##################################################################################################


































   SUBROUTINE lbc_nfd_gpu_2d( npolj, pt2d, cd_type, psgn)
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd_gpu_2d  ***
      !!
      !! ** Purpose :   2D lateral boundary condition : North fold treatment
      !!       without processor exchanges.
      !!
      !! ** Method  :
      !!
      !! ** Action  :   pt2d with updated values along the north fold
      !!----------------------------------------------------------------------
      INTEGER                 , INTENT(in   ) ::   npolj     !: north fold mark (0, 3 or 4)
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! define the nature of ptab array grid-points
      !                                                      ! = T , U , V , F , W points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! control of the sign change
      !                                                      !   = -1. , the sign is changed if north fold boundary
      !                                                      !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d      ! 2D array on which the boundary condition is applied
      !INTEGER , OPTIONAL      , INTENT(in   ) ::   pr2dj     !OPTIONAL! number of additional halos
      !
      INTEGER  ::   ji, jl, ipr2dj
      INTEGER  ::   ijt, iju, ijpj, ijpjm1
      !!----------------------------------------------------------------------
      INTEGER :: jpi, jpj, jpiglo, jpjglo
      !!----------------------------------------------------------------------

      jpi = SIZE(pt2d,1) ; jpiglo = jpi
      jpj = SIZE(pt2d,2) ; jpjglo = jpj


      !SELECT CASE ( jpni )
      !CASE ( 1 )     ;   ijpj = nlcj      ! 1 proc only  along the i-direction
      !CASE DEFAULT   ;   ijpj = 4         ! several proc along the i-direction
      !END SELECT
      !
      !IF( PRESENT(pr2dj) ) THEN           ! use of additional halos
      !   ipr2dj = pr2dj
      !   IF( jpni > 1 )   ijpj = ijpj + ipr2dj
      !ELSE
      !   ipr2dj = 0
      !ENDIF
      !

      !! 1 proc (jpni=1):
      ijpj   = jpj
      ipr2dj = 0


      ijpjm1 = ijpj-1

      !! ORCA2 => npolj => 4
      !! ORCA1 => npolj => 6

      SELECT CASE ( npolj )
         !
      CASE ( 3, 4 )                       ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_type )
            !
         CASE ( 'T' , 'W' )                               ! T- , W-points
            DO jl = 0, ipr2dj
               DO ji = 2, jpiglo
                  ijt=jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            pt2d(1,ijpj)   = psgn * pt2d(3,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo
               ijt=jpiglo-ji+2
               pt2d(ji,ijpj-1) = psgn * pt2d(ijt,ijpj-1)
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj  ) = psgn * pt2d(    2   ,ijpj-2)
            pt2d(jpiglo,ijpj  ) = psgn * pt2d(jpiglo-1,ijpj-2)
            pt2d(1     ,ijpj-1) = psgn * pt2d(jpiglo  ,ijpj-1)
            DO ji = jpiglo/2, jpiglo-1
               iju = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'V' )                                     ! V-point
            DO jl = -1, ipr2dj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-3-jl)
               END DO
            END DO
            pt2d( 1 ,ijpj)   = psgn * pt2d( 3 ,ijpj-3)
         CASE ( 'F' )                                     ! F-point
            DO jl = -1, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-3-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj)   = psgn * pt2d(    2   ,ijpj-3)
            pt2d(jpiglo,ijpj)   = psgn * pt2d(jpiglo-1,ijpj-3)
            pt2d(jpiglo,ijpj-1) = psgn * pt2d(jpiglo-1,ijpj-2)
            pt2d(   1  ,ijpj-1) = psgn * pt2d(    2   ,ijpj-2)
         CASE ( 'I' )                                     ! ice U-V point (I-point)
            DO jl = 0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'J' )                                     ! first ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                     ! second ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE ( 5, 6 )                        ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'W' )                               ! T-, W-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-1)
         CASE ( 'V' )                                     ! V-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(ijt,ijpjm1)
            END DO
         CASE ( 'F' )                               ! F-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo-1
               iju = jpiglo-ji
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'I' )                                  ! ice U-V point (I-point)
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= 0.5 * ( pt2d(ji,ijpj-1-jl) + psgn * pt2d(ijt,ijpj-1-jl) )
               END DO
            END DO
         CASE ( 'J' )                                  ! first ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ji,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                  ! second ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE DEFAULT                           ! *  closed : the code probably never go through
         !
         SELECT CASE ( cd_type)
         CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'F' )                                   ! F-point
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'I' )                                   ! ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'J' )                                   ! first ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'K' )                                   ! second ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         END SELECT
         !
      END SELECT
      !
   END SUBROUTINE lbc_nfd_gpu_2d





   !!======================================================================
END MODULE lbclnk_gpu
