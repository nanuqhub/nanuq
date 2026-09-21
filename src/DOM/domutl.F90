MODULE domutl
   !!======================================================================
   !!                       ***  MODULE  domutl  ***
   !! Grid utilities:
   !!======================================================================
   !! History : 4.2  !  2020-04  (S. Masson)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_ngb       : find the closest grid point from a given lon/lat position
   !!   dom_uniq      : identify unique point of a grid (TUVF)
   !!----------------------------------------------------------------------
   !
   USE dom_oce        ! ocean space and time domain
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! for mppsum

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbnd_ij
      MODULE PROCEDURE arr_lbnd_2d_i,  arr_lbnd_3d_i,  arr_lbnd_4d_i,  arr_lbnd_5d_i,  &
         &             arr_lbnd_2d_sp, arr_lbnd_3d_sp, arr_lbnd_4d_sp, arr_lbnd_5d_sp, &
         &             arr_lbnd_2d_dp, arr_lbnd_3d_dp, arr_lbnd_4d_dp, arr_lbnd_5d_dp
   END INTERFACE lbnd_ij

   PUBLIC dom_ngb    ! routine called in iom.F90 module
   PUBLIC dom_uniq   ! Called by dommsk and domwri
   PUBLIC is_tile
   PUBLIC in_hdom
   PUBLIC lbnd_ij
   PUBLIC arr_hls

   !! * Substitutions
   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_ngb( plon, plat, kii, kjj, cdgrid, kkk )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_ngb  ***
      !!
      !! ** Purpose :   find the closest grid point from a given lon/lat position
      !!
      !! ** Method  :   look for minimum distance in cylindrical projection
      !!                -> not good if located at too high latitude...
      !!----------------------------------------------------------------------
      REAL(wp)        , INTENT(in   ) ::   plon, plat   ! longitude,latitude of the point
      INTEGER         , INTENT(  out) ::   kii, kjj     ! i-,j-index of the closes grid point
      INTEGER         , INTENT(in   ), OPTIONAL :: kkk  ! k-index of the mask level used
      CHARACTER(len=1), INTENT(in   ) ::   cdgrid       ! grid name 'T', 'U', 'V', 'W'
      !
      INTEGER :: ik         ! working level
      INTEGER , DIMENSION(2)      ::   iloc
      REAL(wp)                    ::   zlon, zmini
      REAL(wp), DIMENSION(Nis0:Nie0,Njs0:Nje0) :: zglam, zgphi, zdist
      LOGICAL , DIMENSION(Nis0:Nie0,Njs0:Nje0) :: llmsk
      !!--------------------------------------------------------------------
      !
      ik = 1
      IF( PRESENT(kkk) ) ik=kkk
      !
      SELECT CASE( cdgrid )
      CASE( 'U' ) ; zglam(Nis0:Nie0,Njs0:Nje0) = glamu(Nis0:Nie0,Njs0:Nje0) ; zgphi(Nis0:Nie0,Njs0:Nje0) = gphiu(Nis0:Nie0,Njs0:Nje0) ; llmsk(Nis0:Nie0,Njs0:Nje0) = tmask_i(Nis0:Nie0,Njs0:Nje0) * umask(Nis0:Nie0,Njs0:Nje0,ik) == 1._wp
      CASE( 'V' ) ; zglam(Nis0:Nie0,Njs0:Nje0) = glamv(Nis0:Nie0,Njs0:Nje0) ; zgphi(Nis0:Nie0,Njs0:Nje0) = gphiv(Nis0:Nie0,Njs0:Nje0) ; llmsk(Nis0:Nie0,Njs0:Nje0) = tmask_i(Nis0:Nie0,Njs0:Nje0) * vmask(Nis0:Nie0,Njs0:Nje0,ik) == 1._wp
      CASE( 'F' ) ; zglam(Nis0:Nie0,Njs0:Nje0) = glamf(Nis0:Nie0,Njs0:Nje0) ; zgphi(Nis0:Nie0,Njs0:Nje0) = gphif(Nis0:Nie0,Njs0:Nje0) ; llmsk(Nis0:Nie0,Njs0:Nje0) = tmask_i(Nis0:Nie0,Njs0:Nje0) * fmask(Nis0:Nie0,Njs0:Nje0,ik) == 1._wp
      CASE DEFAULT; zglam(Nis0:Nie0,Njs0:Nje0) = glamt(Nis0:Nie0,Njs0:Nje0) ; zgphi(Nis0:Nie0,Njs0:Nje0) = gphit(Nis0:Nie0,Njs0:Nje0) ; llmsk(Nis0:Nie0,Njs0:Nje0) = tmask_i(Nis0:Nie0,Njs0:Nje0) * tmask(Nis0:Nie0,Njs0:Nje0,ik) == 1._wp
      END SELECT
      !
      zlon       = MOD( plon       + 720., 360. )                                     ! plon between    0 and 360
      zglam(:,:) = MOD( zglam(:,:) + 720., 360. )                                     ! glam between    0 and 360
      IF( zlon > 270. )   zlon = zlon - 360.                                          ! zlon between  -90 and 270
      IF( zlon <  90. )   WHERE( zglam(:,:) > 180. ) zglam(:,:) = zglam(:,:) - 360.   ! glam between -180 and 180
      zglam(:,:) = zglam(:,:) - zlon
      !
      zgphi(:,:) = zgphi(:,:) - plat
      zdist(:,:) = zglam(:,:) * zglam(:,:) + zgphi(:,:) * zgphi(:,:)
      !
      CALL mpp_minloc( 'domngb', zdist(:,:), llmsk, zmini, iloc, ldhalo = .TRUE. )
      kii = iloc(1)
      kjj = iloc(2)
      !
   END SUBROUTINE dom_ngb


   SUBROUTINE dom_uniq( puniq, cdgrd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_uniq  ***
      !!
      !! ** Purpose :   identify unique point of a grid (TUVF)
      !!
      !! ** Method  :   1) aplly lbc_lnk on an array with different values for each element
      !!                2) check which elements have been changed
      !!----------------------------------------------------------------------
      CHARACTER(len=1)        , INTENT(in   ) ::   cdgrd   !
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   puniq   !
      !
      REAL(wp)                       ::  zshift   ! shift value link to the process number
      INTEGER                        ::  ji       ! dummy loop indices
      LOGICAL , DIMENSION(jpi,jpj,1) ::   lluniq  ! store whether each point is unique or not
      REAL(wp), DIMENSION(jpi,jpj  ) ::   ztstref
      !!----------------------------------------------------------------------
      !
      ! build an array with different values for each element
      ! in mpp: make sure that these values are different even between process
      ! -> apply a shift value according to the process number
      zshift = jpimax * jpjmax * ( narea - 1 )
      ztstref(:,:) = RESHAPE( (/ (zshift + REAL(ji,wp), ji = 1, jpi*jpj) /), (/ jpi, jpj /) )
      !
      puniq(:,:) = ztstref(:,:)                       ! default definition
      CALL lbc_lnk( 'domwri', puniq, cdgrd, 1._wp )   ! apply boundary conditions
      lluniq(:,:,1) = puniq(:,:) == ztstref(:,:)      ! check which values have not been changed
      !
      puniq(:,:) = REAL( COUNT( lluniq(:,:,:), dim = 3 ), wp )
      !
   END SUBROUTINE dom_uniq


   FUNCTION arr_lbnd( kasi, kasj, ldtile ) RESULT(kbnd)
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION arr_lbnd  ***
      !!
      !! ** Purpose :   Determine the lower bounds of an array's internal area from its i-j shape
      !!
      !! ** Method  :   1) Identify whether the array shape corresponds to a tile or MPI domain (is_tile)
      !!                2) Compare the array shape to the internal area shape to get the halo size
      !!                3) Determine the lower bounds of the internal area
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           :: kasi, kasj                 ! Size of i, j dimensions
      LOGICAL, INTENT(in), OPTIONAL :: ldtile                     ! Is the array on the tile (T) or MPI (F) domain?
      !
      INTEGER, DIMENSION(2) :: kbnd                               ! Lower bounds of i, j dimensions
      INTEGER, DIMENSION(2) :: ihls                               ! Halo size along i, j dimensions
      LOGICAL :: lltile                                           ! Is the array on the tile (T) or MPI (F) domain?
      !!----------------------------------------------------------------------
      ! Is this array on the tile or MPI domain?
      IF( PRESENT(ldtile) ) THEN
         lltile = ldtile
      ELSE
         lltile = is_tile(kasi, kasj)
      ENDIF

      ! Halo size
      ihls = arr_hls(kasi, kasj, ldtile=lltile, ldsize=.TRUE.)

      ! Lower bounds
      IF( lltile ) THEN                      ! Tile domain
         kbnd(1) = ntsi - ihls(1)
         kbnd(2) = ntsj - ihls(2)
      ELSE                                   ! MPI domain
         kbnd(1) = Nis0 - ihls(1)
         kbnd(2) = Njs0 - ihls(2)
      ENDIF
   END FUNCTION arr_lbnd


   FUNCTION arr_hls( kasi, kasj, ldtile, ldsize ) RESULT(khls)
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION arr_hls  ***
      !!
      !! ** Purpose :   Determine the halo size or number of halo points in an array from its i-j shape
      !!
      !! ** Method  :   Compare the array shape to that of the MPI domain internal area
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           :: kasi, kasj                 ! Size of i, j dimensions
      LOGICAL, INTENT(in), OPTIONAL :: ldtile                     ! Is the array on the tile (T) or MPI (F) domain?
      LOGICAL, INTENT(in), OPTIONAL :: ldsize                     ! Return the halo size (T) or total number of halo points (F)
      !
      INTEGER, DIMENSION(2) :: khls                               ! No. halo points or halo size along i, j dimensions
      LOGICAL               :: lltile, llsize
      !!----------------------------------------------------------------------
      llsize = .FALSE.
      IF( PRESENT(ldsize) ) llsize = ldsize

      IF( PRESENT(ldtile) ) THEN                                  ! Is this array on the tile or MPI domain?
         lltile = ldtile
      ELSE
         lltile = is_tile(kasi, kasj)
      ENDIF

      ! Number of halo points
      IF( lltile ) THEN                                           ! Tile domain (smaller than MPI domain internal area)
         khls(1) = kasi - (ntei - ntsi + 1)
         khls(2) = kasj - (ntej - ntsj + 1)
      ELSE                                                        ! MPI domain
         khls(1) = kasi - Ni_0
         khls(2) = kasj - Nj_0
      ENDIF

      ! Halo size
      IF( llsize ) THEN
         IF( MOD(khls(1), 2) == 1 .OR. MOD(khls(2), 2) == 1 ) THEN
            WRITE(ctmp1,*) 'Cannot determine halo size of an array with an odd number of i-j halo points:'
            WRITE(ctmp2,*) khls(1), 'x', khls(2)
            CALL ctl_stop('STOP', ctmp1, ctmp2)                         ! Stop before out of bounds errors
         ELSE
            khls(1) = khls(1) / 2
            khls(2) = khls(2) / 2
         ENDIF
      ENDIF
   END FUNCTION arr_hls


   LOGICAL FUNCTION check_hdom( kasi, kasj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION check_hdom  ***
      !!
      !! ** Purpose :   Determine if an array's shape corresponds to the horizontal domain
      !!
      !! ** Method  :   The array must be at least as large as the internal domain,
      !!                but have no more than 2 * nn_hls halo points in total.
      !!                This must be true for either the tile or MPI domain.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kasi, kasj                           ! Size of i, j dimensions
      INTEGER :: imini, imaxi, iminj, imaxj
      LOGICAL :: lltile, llmpi                                    ! Whether array is on the tile / MPI domain
      !!----------------------------------------------------------------------
      imini = ntei - ntsi + 1 ; iminj = ntej - ntsj + 1
      imaxi = imini + 2 * nn_hls ; imaxj = iminj + 2 * nn_hls

      lltile = (imini <= kasi .AND. kasi <= imaxi) .AND. (iminj <= kasj .AND. kasj <= imaxj)
      llmpi = (Ni_0 <= kasi .AND. kasi <= jpi) .AND. (Nj_0 <= kasj .AND. kasj <= jpj)
      check_hdom = lltile .OR. llmpi
   END FUNCTION check_hdom


   LOGICAL FUNCTION is_tile( kasi, kasj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION is_tile  ***
      !!
      !! ** Purpose :   Determine whether an array's shape corresponds to the tile or MPI domain
      !!
      !! ** Method  :   Compare the array shape to that of the MPI domain internal area.
      !!                The maximum tile size is limited (in domain.F90) to ensure this is not ambiguous.
      !!                If the array shape is not valid for either the tile or MPI domain, raise an error.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kasi, kasj                           ! Size of i, j dimensions
      !!----------------------------------------------------------------------
      ! kas[ij] must represent a valid shape for either the tile or MPI domain
      IF( .NOT. check_hdom(kasi, kasj) ) THEN
         WRITE(ctmp1,*) 'is_tile: the following i-j shape does not represent an array on the horizontal domain:'
         WRITE(ctmp2,*) kasi, 'x', kasj
         CALL ctl_stop('STOP', ctmp1, ctmp2)                         ! Stop before out of bounds errors
      ENDIF

      is_tile = l_istiled .AND. (kasi < Ni_0 .OR. kasj < Nj_0)
   END FUNCTION is_tile


   ELEMENTAL FUNCTION in_hdom( ki, kj, khls, cddom )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION in_hdom  ***
      !!
      !! ** Purpose :  Check whether a given index is within the current tile or MPI domain
      !!
      !! ** Method  :  Returns true if point (ki,kj) is within:
      !!
      !!                   * Nis0:Nie0, Njs0:Nje0  (cddom='mpi', MPI domain)
      !!                   * ntsi:ntei, ntsj:ntej  (cddom='tile', tile domain)
      !!
      !!               `khls` halo points will be added to these bounds if specified.
      !!----------------------------------------------------------------------
      LOGICAL                                :: in_hdom
      INTEGER,          INTENT(in)           :: ki, kj            ! Indices of point to check
      INTEGER,          INTENT(in), OPTIONAL :: khls              ! Number of halo points to include in domain
      CHARACTER(len=*), INTENT(in), OPTIONAL :: cddom             ! Domain to check ('tile'|'mpi')
      !
      INTEGER :: ihls, isi, iei, isj, iej
      CHARACTER(len=10) :: cldom
      !!----------------------------------------------------------------------
      cldom = 'tile'
      IF( PRESENT(cddom) ) cldom = cddom
      ihls = 0
      IF( PRESENT(khls) ) ihls = khls

      SELECT CASE( TRIM(cldom) )
      CASE( 'mpi' )        ! The MPI case is redundant but included for completeness
         isi = Nis0 - ihls ; isj = Njs0 - ihls
         iei = Nie0 + ihls ; iej = Nje0 + ihls
      CASE DEFAULT
         isi = ntsi - ihls ; isj = ntsj - ihls
         iei = ntei + ihls ; iej = ntej + ihls
      END SELECT
      in_hdom = (ki >= isi) .AND. (ki <= iei) .AND. (kj >= isj) .AND. (kj <= iej)
   END FUNCTION in_hdom


   FUNCTION arr_lbnd_2d_i( pt ) RESULT(kbnd)
      INTEGER, DIMENSION(:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_2d_i


   FUNCTION arr_lbnd_2d_sp( pt ) RESULT(kbnd)
      REAL(sp), DIMENSION(:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_2d_sp


   FUNCTION arr_lbnd_2d_dp( pt ) RESULT(kbnd)
      REAL(dp), DIMENSION(:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_2d_dp


   FUNCTION arr_lbnd_3d_i( pt ) RESULT(kbnd)
      INTEGER, DIMENSION(:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_3d_i


   FUNCTION arr_lbnd_3d_sp( pt ) RESULT(kbnd)
      REAL(sp), DIMENSION(:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_3d_sp


   FUNCTION arr_lbnd_3d_dp( pt ) RESULT(kbnd)
      REAL(dp), DIMENSION(:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_3d_dp


   FUNCTION arr_lbnd_4d_i( pt ) RESULT(kbnd)
      INTEGER, DIMENSION(:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_4d_i


   FUNCTION arr_lbnd_4d_sp( pt ) RESULT(kbnd)
      REAL(sp), DIMENSION(:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_4d_sp


   FUNCTION arr_lbnd_4d_dp( pt ) RESULT(kbnd)
      REAL(dp), DIMENSION(:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_4d_dp


   FUNCTION arr_lbnd_5d_i( pt ) RESULT(kbnd)
      INTEGER, DIMENSION(:,:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_5d_i


   FUNCTION arr_lbnd_5d_sp( pt ) RESULT(kbnd)
      REAL(sp), DIMENSION(:,:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_5d_sp


   FUNCTION arr_lbnd_5d_dp( pt ) RESULT(kbnd)
      REAL(dp), DIMENSION(:,:,:,:,:), INTENT(in) ::   pt
      INTEGER, DIMENSION(2) :: kbnd

      kbnd = arr_lbnd(SIZE(pt, 1), SIZE(pt, 2))
   END FUNCTION arr_lbnd_5d_dp
   !!======================================================================
END MODULE domutl
