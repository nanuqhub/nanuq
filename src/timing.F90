MODULE timing
   !!========================================================================
   !!                     ***  MODULE  timing  ***
   !!========================================================================
   !! History : 4.0  ! 2001-05  (R. Benshila)
   !!           4.x  ! 2023-05  (G. Irrmann, S. Masson)
   !!------------------------------------------------------------------------

   !!------------------------------------------------------------------------
   !!   timing_open     : open timing_nanuq.output file
   !!   timing_start    : start Timerg
   !!   timing_stop     : stop  Timer
   !!   timing_finalize : compute stats and write timing_nanuq.output
   !!------------------------------------------------------------------------
   USE par_kind, ONLY: dp, i4
   USE par_oce , ONLY: ntile, jpni, jpnj, jpnij
   USE dom_oce , ONLY: l_istiled, narea, nimpi, njmpi, Agrif_Root, Agrif_CFixed, rn_Dt
   !! WARNING: we cannot use lib_mpp because of circular dependencies

   USE netcdf          ! NetCDF library
#if ! defined key_mpi_off
   USE MPI
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   timing_start, timing_stop, timing_open      ! called in each routine to time

   INTEGER, PARAMETER :: nszmax = 1000000   ! 10 * 100 * 1000 -> equivalent of ten (33x33x100) 3D arrays
   INTEGER, PARAMETER :: ntest = 2
   INTEGER, DIMENSION(ntest) :: n1st = 0
   INTEGER, DIMENSION(ntest) :: nend = 0
   LOGICAL, DIMENSION(ntest) :: luse = .FALSE.

   ! Variables for fine grain timing
   TYPE timer
      CHARACTER(LEN=32)  :: cname, cleanm
      INTEGER(8)                   :: n8start, n8childsum
      INTEGER(8), DIMENSION(ntest) :: n8tnet, n8tfull
      INTEGER   , DIMENSION(ntest) :: ncalls
      INTEGER  :: nrecurrent
      INTEGER  :: nwrt, nchunk, ncid
      LOGICAL  :: ldone
      REAL(dp), DIMENSION(:), ALLOCATABLE :: tkeep
      REAL(dp) ::  tnet , tnetavg,  tnetmin,  tnetmax,  tnetblc   ! net  time average, min, max and load balance (max-min)
      REAL(dp) :: tfull, tfullavg, tfullmin, tfullmax, tfullblc   ! full time average, min, max and load balance (max-min)
      TYPE(timer), POINTER :: s_next, s_prev, s_parent
   END TYPE timer

   TYPE(timer), POINTER :: s_timer_root => NULL()
   TYPE(timer), POINTER :: s_timer      => NULL()

   INTEGER    :: numtime   = -1      ! logical unit for timing
   INTEGER    :: nrunid
   INTEGER    :: ncall_clock

   INTEGER    :: nmpicom   !: we cannot use mpi_comm_oce as we cannot use lib_mpp
   INTEGER(8) :: n8start000
   REAL(dp)   :: secondclock

   INTEGER, PARAMETER :: jp_cname = 1
   INTEGER, PARAMETER :: jp_tnet  = 2   ! local net time (of the current MPI process)
   INTEGER, PARAMETER :: jp_tavg  = 3   ! mean net time (among MPI processes)
   INTEGER, PARAMETER :: jp_tmin  = 4   ! min net time (among MPI processes)
   INTEGER, PARAMETER :: jp_tmax  = 5   ! max net time (among MPI processes)
   INTEGER, PARAMETER :: jp_tblc  = 6   ! load balance between MPI processes (max-min)

   INTEGER, PARAMETER :: jpmaxline = 50   ! max number of line to be printed
   INTEGER, PARAMETER :: nbequal   = 70   ! number of "=" sign in separation print lines

   INTEGER, DIMENSION(8)           :: nvalues
   CHARACTER(LEN= 8), DIMENSION(2) :: cdate
   CHARACTER(LEN=10), DIMENSION(2) :: ctime
   CHARACTER(LEN= 5)               :: czone

   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE timing_open( ldwp, kmpicom, cdname )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_open  ***
      !! ** Purpose :   Open timing output file.
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(in) :: ldwp
      INTEGER                   , INTENT(in) :: kmpicom
      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cdname
      !
      CHARACTER(len=32) ::   cln
      CHARACTER(LEN=10) ::   clfmt
      INTEGER           ::   idg
      !!----------------------------------------------------------------------

      IF( .NOT. Agrif_Root() ) RETURN

      nmpicom = kmpicom

      IF( PRESENT(cdname) ) THEN
         cln = cdname
      ELSE
         cln = 'timing_nanuq.output'
      ENDIF

      IF( ldwp ) THEN
         ! we cannot use ctl_open as we cannot use lib_mpp
         OPEN(NEWUNIT = numtime, FILE = TRIM(cln)//add_xxx(), STATUS = "REPLACE", ACTION = "write")
         WRITE(numtime,*)
         WRITE(numtime,*) '      CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC'
         WRITE(numtime,*) '                             NEMO team'
         WRITE(numtime,*) '                  Ocean General Circulation Model'
         WRITE(numtime,*) '                        version 5.0  (2024) '
         WRITE(numtime,*)
         WRITE(numtime,*) '                        Timing Informations '
      ENDIF
      !
   END SUBROUTINE timing_open


   SUBROUTINE timing_start( cdinfo, kt, kt000, ktend, kfsbc, keepnc )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_start  ***
      !! ** Purpose :   collect execution time
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) :: cdinfo
      INTEGER, OPTIONAL, INTENT(in) :: kt, kt000, ktend, kfsbc, keepnc
      !
      CHARACTER(LEN=32)    :: clinfo
      TYPE(timer), POINTER :: s_wrk
      INTEGER(8)           :: i8rate
      INTEGER              :: ji
      !!----------------------------------------------------------------------
      !
      IF( .NOT. luse(1) )   luse(1) = .TRUE.   ! switch on the timing for the entire simulation
      !
      IF( PRESENT(kt000) .AND. n1st(1) == 0 ) THEN   ! definition of n1st and nend according to kt000, ktend and kfsbc
         n1st(1) = kt000
         nend(1) = ktend
         n1st(2) = kt000 + 3 * kfsbc   ! exclude open(write) communication report at kt000 + 1(2) * kfsbc
         nend(2) = ktend - 2 * kfsbc   ! exclude open(write) restarts at ktend - 2(1) * kfsbc + 1
      ENDIF
      IF( PRESENT(kt) ) THEN
         DO ji = 2, ntest
            luse(ji) = kt >= n1st(ji) .AND. kt <= nend(ji)
         END DO
      ENDIF
      !
      clinfo = cdinfo
      IF( .NOT. Agrif_Root() )   clinfo = TRIM(Agrif_CFixed())//'_'//clinfo

      IF( .NOT. ASSOCIATED(s_timer_root) ) THEN          ! this is the first call to timing_start
         s_timer_root => def_newlink( clinfo, keepnc )   ! define the root link
         CALL SYSTEM_CLOCK( COUNT_RATE = i8rate )        ! define rateclock
         secondclock = 1._dp / REAL(i8rate, dp)
         CALL DATE_AND_TIME( cdate(1), ctime(1), czone, nvalues )
         ncall_clock = 0
      ENDIF

      IF( ASSOCIATED(s_timer) ) THEN                     ! at the first call, s_timer is not yet associated
         IF( s_timer%cname == clinfo ) THEN              ! we are timing a recurrent routine (a routine directly calling itself)
            s_timer%nrecurrent = s_timer%nrecurrent + 1
            RETURN
         ENDIF
      ENDIF

      ! store s_timer chain link in s_wrk
      s_wrk => s_timer
      !
      ! make s_timer pointing toward the chain link corresponding to clinfo
      s_timer => find_link( clinfo, s_timer_root, keepnc )

      ! we must take care if we do a timing inside another timing...
      ! if s_wrk did not finish is timing, this means that s_timer in part of s_wrk.
      ! in this case we link s_timer to s_wrk
      s_timer%s_parent => NULL()             ! default not the part of another timing
      IF( ASSOCIATED(s_wrk) ) THEN
         IF( .NOT. s_wrk%ldone )   s_timer%s_parent => s_wrk
      ENDIF

      ! initialisation
      s_timer%ldone = .FALSE.                ! we are just starting the timing (not done)
      s_timer%n8childsum = 0_8               ! not yet any my children count

      ! clock time collection
      CALL SYSTEM_CLOCK( COUNT = s_timer%n8start )   ;   ncall_clock = ncall_clock + 1
      IF( ncall_clock == 1 )   n8start000 =  s_timer%n8start   ! keep track of the first call

   END SUBROUTINE timing_start


   SUBROUTINE timing_stop( cdinfo, kt, ld_finalize )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_stop  ***
      !! ** Purpose :   stop timing window
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdinfo
      INTEGER, OPTIONAL, INTENT(in) ::   kt
      LOGICAL, OPTIONAL, INTENT(in) ::   ld_finalize
      !
      CHARACTER(LEN=32)  :: clinfo
      INTEGER(8) :: i8end, i8full, i8net
      INTEGER    :: ji
      INTEGER    :: icode
      INTEGER    :: ivid
      LOGICAL    :: ll_finalize
      !!----------------------------------------------------------------------
      !
      IF( s_timer%nrecurrent > 0 ) THEN
         s_timer%nrecurrent = s_timer%nrecurrent - 1   ! we are timing a recurrent routine (a routine directly calling itself)
         RETURN
      ENDIF
      !
      clinfo = cdinfo
      IF( .NOT. Agrif_Root() )   clinfo = TRIM(Agrif_CFixed())//'_'//clinfo

      IF( s_timer%cname /= clinfo )   CALL local_stop('try to stop '//TRIM(clinfo)//' but we point toward '//TRIM(s_timer%cname))

      IF( PRESENT(ld_finalize) ) THEN
         ll_finalize = ld_finalize
      ELSE
         ll_finalize = .FALSE.
      ENDIF

      IF( .NOT. l_istiled .OR. ntile == 1 ) THEN
         DO ji = 1, ntest
            IF( luse(ji) )   s_timer%ncalls(ji)  = s_timer%ncalls(ji) + 1      ! All tiles count as one iteration
         END DO
      ENDIF

      ! clock time collection
      CALL SYSTEM_CLOCK( COUNT = i8end )   ;   ncall_clock = ncall_clock + 1
      i8full = i8end - s_timer%n8start                                         ! count between  timing_start and timing_stop
      DO ji = 1, ntest
         IF( luse(ji) )   s_timer%n8tfull(ji) = s_timer%n8tfull(ji) + i8full   ! cumulate my full time
      END DO

      ! store timing
      IF( ALLOCATED(s_timer%tkeep) ) THEN
         s_timer%nwrt = s_timer%nwrt + 1
         s_timer%tkeep(s_timer%nwrt) = REAL(i8full,dp) * secondclock
         IF( s_timer%nwrt == SIZE(s_timer%tkeep) ) THEN
            IF( s_timer%ncid /= -1 ) THEN   ! write timing chunk to NetCDF
               CALL nf90chk( NF90_INQ_VARID(s_timer%ncid, 'timing_'//TRIM(s_timer%cleanm), ivid) )
               CALL nf90chk( NF90_PUT_VAR(  s_timer%ncid, ivid, s_timer%tkeep,   &
                  &                         (/s_timer%nchunk * s_timer%nwrt + 1/), (/s_timer%nwrt/) ) )
            ENDIF
            s_timer%nwrt   = 0
            s_timer%nchunk = s_timer%nchunk + 1
         ENDIF
      ENDIF

      ! time diagnostics
      i8net = i8full - s_timer%n8childsum                                      ! don't take into account my cildren count
      DO ji = 1, ntest
         IF( luse(ji) )   s_timer%n8tnet(ji) = s_timer%n8tnet(ji) + i8net      ! cumulate my net time
      END DO
      s_timer%ldone = .TRUE.                                                   ! I am done with this counting

      ! we come back to the parent
      s_timer => s_timer%s_parent
      IF( ASSOCIATED(s_timer) )   s_timer%n8childsum = s_timer%n8childsum + i8full   ! add myself to the children of my parents

      IF( ll_finalize ) CALL timing_finalize( s_timer_root )

   END SUBROUTINE timing_stop


   SUBROUTINE timing_finalize( sd_root )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_finalize ***
      !! ** Purpose :  compute average time
      !!               write timing output file
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      !
      TYPE(timer), POINTER :: s_wrk
      REAL(dp) :: zmytime, zmysum, zval, zavgtime, zavgsum, zavgextra, zmin, zmax, zblc, zldsum, zsypd
      REAL(dp), DIMENSION(:), ALLOCATABLE :: zalltime
      INTEGER :: inb
      INTEGER :: jpnbtest = 100
      INTEGER :: ji
      INTEGER :: icode, istatus, icmdstat
      INTEGER, DIMENSION(:), ALLOCATABLE :: iallmpi
      INTEGER(8) :: i8start, i8end
      CHARACTER(len=128) :: clfmt, cline
      LOGICAL :: ll_avg, llwrt, ll_execmd
      !!----------------------------------------------------------------------

      llwrt = numtime /= -1

      ! check that all timing are done (all childrend count are properly reported)
      s_wrk => sd_root
      DO WHILE( ASSOCIATED(s_wrk) )
         IF( .NOT. s_wrk%ldone )   CALL local_stop( TRIM(s_wrk%cname)//'did not finish its timing' )
         s_wrk => s_wrk%s_next
      END DO

      ! get elapsed time between now and the first call to timing start
      CALL SYSTEM_CLOCK( COUNT = i8end )
      zmytime = REAL(i8end - n8start000, dp) * secondclock

      ! End of timings on date & time
      CALL DATE_AND_TIME( cdate(2), ctime(2), czone, nvalues )
      !
      IF( llwrt ) THEN
         CALL write_header('Timing report:')
         clfmt='(1X,"Timing started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
         WRITE(numtime, clfmt) &
            &       cdate(1)(7:8), cdate(1)(5:6), cdate(1)(1:4),   &
            &       ctime(1)(1:2), ctime(1)(3:4), ctime(1)(5:6),   &
            &       czone(1:3),    czone(4:5)
         clfmt='(1X,"Timing   ended on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
         WRITE(numtime, clfmt) &
            &       cdate(2)(7:8), cdate(2)(5:6), cdate(2)(1:4),   &
            &       ctime(2)(1:2), ctime(2)(3:4), ctime(2)(5:6),   &
            &       czone(1:3),    czone(4:5)
      ENDIF

#if ! defined key_mpi_off
      ! Compute the number of routines
      inb = nblinks( sd_root )
      ALLOCATE(iallmpi(jpnij))
      CALL MPI_ALLGATHER(    inb , 1, MPI_INTEGER,   &
         &                iallmpi, 1, MPI_INTEGER, nmpicom, icode)
      IF( SUM( iallmpi ) /= inb*jpnij ) THEN
         IF( llwrt ) THEN
            CALL write_separator()
            WRITE(numtime,*) '        ===> W A R N I N G: '
            WRITE(numtime,*) ' Some CPU have different number of routines instrumented for timing'
            WRITE(numtime,*) ' No detailed report on averaged timing can be provided'
            WRITE(numtime,*) ' The following detailed report only deals with the current processor'
         ENDIF
         ll_avg = .FALSE.
      ELSE
         ll_avg = jpnij > 1
      ENDIF
      DEALLOCATE(iallmpi)
#else
      ll_avg = .FALSE.
#endif

      ! get elapsed time avg/min/max
      CALL mpi_avgminmax( zmytime, zavgtime, zmin, zmax, ll_avg )
      IF( llwrt ) THEN
         CALL write_header('Elapsed Time mesured by timing (s):')
         clfmt = "(a,f0.6,'s')"
         IF( ll_avg .and. narea == 1 ) THEN
            WRITE(numtime,clfmt) '       avg over all MPI processes = ', zavgtime
            WRITE(numtime,clfmt) '       min over all MPI processes = ', zmin
            WRITE(numtime,clfmt) '       max over all MPI processes = ', zmax
            WRITE(numtime,clfmt) '              local MPI process   = ', zmytime
         ELSE
            WRITE(numtime,clfmt) '                    local process = ', zmytime
         ENDIF
      ENDIF

      ! add an aveluation of the timing itself...
      CALL SYSTEM_CLOCK( COUNT = i8start )
      DO ji = 1, jpnbtest
         CALL SYSTEM_CLOCK( COUNT = i8end )
      ENDDO
      zval = REAL( (i8end-i8start) * ncall_clock, dp) * secondclock / REAL(jpnbtest, dp)
      CALL mpi_avgminmax( zval, zavgextra, zmin, zmax, ll_avg )

      IF( llwrt ) THEN
         CALL write_header('Evaluation of the extra coast due to the timing itself (% of avg elapsed):')
         WRITE(numtime,*) '   Number calls to SYSTEM_CLOCK = ', ncall_clock
         WRITE(numtime,'(a,i3,a)') '    Avg Estimation over ', jpnbtest,' tests'
         clfmt ="(a,f0.6,'s (',f0.3,'%)')"
         IF( ll_avg .and. narea == 1 ) THEN
            WRITE(numtime,clfmt) '       avg over all MPI processes = ', zavgextra, zavgextra / zavgtime * 100._dp
            WRITE(numtime,clfmt) '       min over all MPI processes = ', zmin, zmin / zavgtime * 100._dp
            WRITE(numtime,clfmt) '       max over all MPI processes = ', zmax, zmax / zavgtime * 100._dp
            WRITE(numtime,clfmt) '              local MPI process   = ', zval, zval / zavgtime * 100._dp
         ELSEIF( narea > 1 ) THEN
            WRITE(numtime,clfmt) '                    local process = ', zval, zval / zavgtime * 100._dp
         ELSE
            WRITE(numtime,clfmt) '                    local process = ', zval, 100._dp
         ENDIF
      ENDIF

      ! write all kept timing time series in NetCDF files
      call write_allts_nc( sd_root )

      ! write gnuplot script to get statistics and plot on timing time series (e.g. step)
      IF( narea == 1 )   CALL write_gnuplotscript()   ! write gnuplot script
      IF( llwrt ) THEN
         CALL write_header('Info: gnuplots statistics and plots')
         CALL EXECUTE_COMMAND_LINE('which gnuplot &> /dev/null', EXITSTAT = istatus, CMDSTAT = icmdstat)   ! is gnuplot found?
         ll_execmd = istatus == 0 .AND. icmdstat == 0
         IF( ll_execmd ) THEN                         ! execute gnuplot script
            WRITE(numtime,*) '        Shell script used to create png plots: timing_nanuq_gnuplot.sh'
         ELSE
            WRITE(numtime,*) '        gnuplot not found but we created the shell script timing_nanuq_gnuplot.sh'
            WRITE(numtime,*) '        that can be used later to create gnuplot stats and analyses'
         ENDIF
         WRITE(numtime,*) '        Usage:'
         WRITE(numtime,*) '           ./timing_nanuq_gnuplot.sh'
         WRITE(numtime,*) '           ./timing_nanuq_gnuplot.sh --png    # save plots in png'
         WRITE(numtime,*) '           ./timing_nanuq_gnuplot.sh --help   # see all options'
      ELSE
         ll_execmd = .FALSE.
      ENDIF

      ! write performance for each time window and each link of the chain
      DO ji = 1, ntest   ! for each time window

         IF( n1st(ji) == 0 .OR. nend(ji) < n1st(ji) )   CYCLE   ! no timing taken

         ! get avg/min/max for each time window and each link of the chain, write 'timing_mpitasks.nc'
         s_wrk => sd_root
         DO WHILE( ASSOCIATED(s_wrk) )
            s_wrk%tnet  = REAL(s_wrk%n8tnet( ji), dp) * secondclock
            s_wrk%tfull = REAL(s_wrk%n8tfull(ji), dp) * secondclock
            s_wrk => s_wrk%s_next
         END DO
         IF( ll_avg                 )   CALL mpi_t_avgminmax(   sd_root, ji )
         IF( ll_avg .AND. jpnij > 1 )   CALL write_tsum_allmpi( sd_root, ji )

         zmysum = 0._dp
         zldsum = 0._dp
         s_wrk => sd_root
         DO WHILE ( ASSOCIATED(s_wrk) )
            zmysum = zmysum + s_wrk%tnet

            IF( ll_avg .AND. narea == 1 )   zldsum = zldsum + s_wrk%tnetblc
            s_wrk => s_wrk%s_next
         END DO
         CALL mpi_avgminmax( zmysum, zavgsum, zmin, zmax, ll_avg )

         IF( llwrt ) THEN

            IF( ji == 1 ) THEN
               cline = 'Performance statistics for the ENTIRE simulation'
            ELSE
               WRITE(cline, "('Performance statistics between time step ',i0,' and ',i0)") n1st(ji), nend(ji)
            ENDIF
            CALL write_bigheader(cline)

            ! Print model speed
            IF( ji == 2 ) THEN

               zsypd  = 0._dp
               ! convert time-to-solution (inner time steps)
               ! to speed in simulated years per day (SYPD)
               IF( zmax > 0. ) THEN
                  zsypd = ( nend(ji) - n1st(ji) + 1 ) * rn_Dt / ( zmax * 365._dp )
                  WRITE(numtime,*) ' '
                  WRITE(numtime,'(a,f0.2)') '            Model speed (SYPD) : ', zsypd
                  WRITE(numtime,*) ' '
               ELSE
                  WRITE(numtime,*) ' '
                  WRITE(numtime,*) ' unable to calculate SYPD '
                  WRITE(numtime,*) ' '
               ENDIF

            ENDIF

            ! execute gnuplot script to get statistics and plot on timing time series (e.g. step) for the given time window
            IF( ll_execmd )   CALL gnuplot_statplot( sd_root, ji )
            IF( ll_avg .AND. narea == 1 ) THEN
               CALL timer_write( 'Timing : AVG values over all MPI processes:', sd_root, ji, jp_tavg, zavgsum, zavgextra )
               CALL timer_write( 'Timing : MAX load unbalance over all MPI processes (max-min) :', sd_root, ji, jp_tblc, zldsum, 0._dp )
               CALL timer_write( 'Timing : MIN values over all MPI processes:', sd_root, ji, jp_tmin, zavgsum, zavgextra )
               CALL timer_write( 'Timing : MAX values over all MPI processes:', sd_root, ji, jp_tmax, zavgsum, zavgextra )
            ENDIF
            CALL timer_write( 'Timing : values for local MPI process:', sd_root, ji, jp_tnet, zmysum, zavgextra )
         ENDIF

      END DO

      IF( llwrt ) THEN
         CALL write_separator()
         CLOSE(numtime)
      ENDIF
      !
   END SUBROUTINE timing_finalize


   SUBROUTINE timer_write( cdinfo, sd_root, kn, kswitch, ptimetot, ptextra )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wcurrent_info ***
      !! ** Purpose :  compute and write timing output file
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      INTEGER             , INTENT(in   ) :: kn           ! time window index
      INTEGER,              INTENT(in   ) :: kswitch      ! select the variable to be printed
      REAL(dp),             INTENT(in   ) :: ptimetot     ! total elapsed time
      REAL(dp),             INTENT(in   ) :: ptextra      ! extrat total elapsed time from the calls to system_clock
      !
      TYPE(timer), POINTER :: s_wrk
      REAL(dp) :: ztnet, ztfull, zpcent
      INTEGER  :: icnt
      INTEGER  :: idgnet, idgfull, idgpct
      INTEGER  :: itot0, itot1, itot2, itot3
      INTEGER  :: iblk1, iblk2
      CHARACTER(len= 32) :: cltitle0, cltitle1, cltitle2, cltitle3
      CHARACTER(len=128) :: clfmt, clflt1, clflt2
      LOGICAL :: llwarning
      !!----------------------------------------------------------------------

      IF( ptimetot == 0._dp )   RETURN

      llwarning = .FALSE.
      zpcent = 100._dp / ptimetot

      CALL sort_chain( sd_root, kswitch )   ! reorder the current list by decreassing values of kswitch

      IF(     kswitch == jp_tnet ) THEN
         ztnet = sd_root%tnet
      ELSEIF( kswitch == jp_tavg ) THEN
         ztnet = sd_root%tnetavg
      ELSEIF( kswitch == jp_tmin ) THEN
         ztnet = sd_root%tnetmin
      ELSEIF( kswitch == jp_tmax ) THEN
         ztnet = sd_root%tnetmax
      ELSEIF( kswitch == jp_tblc ) THEN
         ztnet = sd_root%tnetblc
      ENDIF
      idgnet  = INT(LOG10(MAX(1._dp,ztnet   ))) + 1 + 7   ! how many digits to we need to write? add 7 difgits for '.6'
      idgfull = INT(LOG10(MAX(1._dp,ptimetot))) + 1 + 7   ! how many digits to we need to write? add 7 difgits for '.6'
      idgpct  =                                   3 + 4   ! add 4 digits for '.3'

      cltitle0 = 'Section'
      cltitle1 = 'Net elapsed Time'
      cltitle2 = 'Full elapsed Time'
      cltitle3 = 'Frequency'

      itot0 = MAX(LEN_TRIM(cltitle0),               12)
      itot1 = MAX(LEN_TRIM(cltitle1), idgnet +idgpct+6)   ;   iblk1 = itot1 - (idgnet +idgpct+6) + 1   ! add 6 for 's ( %)'
      itot2 = MAX(LEN_TRIM(cltitle1), idgfull+idgpct+6)   ;   iblk2 = itot2 - (idgfull+idgpct+6) + 1   ! add 6 for 's ( %)'
      itot3 = MAX(LEN_TRIM(cltitle3),                9)

      ! write current info
      CALL write_header(cdinfo)
      WRITE(clfmt, '(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
         &         '(1x,a', itot0, ",' | ',a", itot1, ",' | ',a", itot2, ",' | ',a", itot3, ')'
      WRITE(numtime,clfmt) cltitle0, cltitle1, cltitle2, cltitle3

      WRITE(clfmt, "('(1x,',i2,a,i2,a,i2,a,i2,a)")   &
         &         itot0+1, "('-'),'|',", itot1+2, "('-'),'|',", itot2+2, "('-'),'|',", itot3+1,"('-'))"
      WRITE(numtime,clfmt)

      WRITE(clflt1, "('f',i2.2,'.6,',a,',f',i2.2,'.3,',a)") idgnet , "'s ('", idgpct, "' %) |',"    ! "fx.6,'s (',fx.3,' %) |',"
      WRITE(clflt2, "('f',i2.2,'.6,',a,',f',i2.2,'.3,',a)") idgfull, "'s ('", idgpct, "' %) |',i"   ! "fx.6,'s (',fx.3,' %) | ,i'"
      WRITE(clfmt, "('(1x,a',i2.2,a,i2,'x,',a,i2,'x,',a,i2.2,')')")   &
         &         itot0, ",' |',", iblk1, TRIM(clflt1), iblk2, TRIM(clflt2), itot3+1

      icnt = 0
      s_wrk => sd_root
      DO WHILE ( ASSOCIATED(s_wrk) .AND. icnt < jpmaxline )   ! we print only the first imaxline lines
         IF(     kswitch == jp_tnet ) THEN
            ztnet = s_wrk%tnet
            ztfull = s_wrk%tfull
         ELSEIF( kswitch == jp_tavg ) THEN
            ztnet = s_wrk%tnetavg
            ztfull = s_wrk%tfullavg
         ELSEIF( kswitch == jp_tmin ) THEN
            ztnet = s_wrk%tnetmin
            ztfull = s_wrk%tfullmin
         ELSEIF( kswitch == jp_tmax ) THEN
            ztnet = s_wrk%tnetmax
            ztfull = s_wrk%tfullmax
         ELSEIF( kswitch == jp_tblc ) THEN
            ztnet = s_wrk%tnetblc
            ztfull = s_wrk%tfullblc
         ENDIF
         IF( ztnet == 0._dp )   RETURN
         IF( ztnet < ptextra .AND. .NOT. llwarning ) THEN
            WRITE(numtime,"('WARNING: timings below are smaller than the estimation of the timing itself: ',f0.6,'s')") ptextra
            llwarning = .TRUE.
         ENDIF
         WRITE(numtime,clfmt) s_wrk%cname, ztnet, ztnet*zpcent, ztfull, ztfull*zpcent, s_wrk%ncalls(kn)
         s_wrk => s_wrk%s_next
         icnt = icnt + 1
      END DO

      IF(  ASSOCIATED(s_wrk) )   WRITE(numtime,*) '...'   ! show that there is still more lines that could have been printed
      !
   END SUBROUTINE timer_write


   FUNCTION find_link( cdinfo, sd_root, keepnc )   RESULT( ptr )
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION find_link  ***
      !! ** Purpose :   find the link named cdinfo in the chain link starting with the link sd_root
      !!----------------------------------------------------------------------
      CHARACTER(len=*)    , INTENT(in   ) :: cdinfo
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      INTEGER, OPTIONAL,    INTENT(in   ) :: keepnc
      !
      TYPE(timer), POINTER :: ptr, s_wrk
      !!----------------------------------------------------------------------

      ! case of already existing area (typically inside a loop)
      ptr => sd_root
      DO WHILE( ASSOCIATED(ptr) )
         IF( ptr%cname == cdinfo ) THEN
            IF( .NOT. ptr%ldone )   &
               CALL local_stop('Indirect recurrence (A->B->..->A) with '//TRIM(s_timer%cname)//'?? This is not coded...')
            RETURN   ! cdinfo is already in the chain
         ENDIF
         s_wrk => ptr                         ! store ptr in s_wrk
         ptr   => ptr%s_next
      END DO

      ! cdinfo not found, we reach the end of the chain list -> ptr is NULL()
      ptr => s_wrk   ! go back to the last chain link.

      ! we are at the end of the chain -> add a new chain link
      ptr%s_next => def_newlink( cdinfo, keepnc )   ! define a new chain link and link it to the end of the chain
      ptr%s_next%s_prev => ptr                          ! link the new chain link to the current chain link
      ptr => ptr%s_next                                 ! move the current chain link to this new chain link
      !
   END FUNCTION find_link


   FUNCTION def_newlink( cdinfo, keepnc )   RESULT( ptr )
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION def_newlink  ***
      !! ** Purpose :   add a new link to the chain link
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) :: cdinfo
      INTEGER, OPTIONAL, INTENT(in) :: keepnc
      !
      TYPE(timer), POINTER :: ptr
      INTEGER :: idtime, ivid, ioldMode
      INTEGER :: ji
      INTEGER :: ikeepnc
      CHARACTER(LEN=1 ) :: cl1
      CHARACTER(LEN=64) :: clname
      !!----------------------------------------------------------------------

      IF( PRESENT(keepnc) ) THEN
         ikeepnc = keepnc
      ELSE
         ikeepnc = 0
      ENDIF

      ALLOCATE(ptr)   ! allocate memory space associated with ptr
      ! default required definitions
      ptr%cname      = cdinfo
      ptr%cleanm     = ''
      DO ji = 1, LEN_TRIM(cdinfo)
         cl1 = cdinfo(ji:ji)
         IF( ('0' <= cl1 .AND. cl1 <= '9') .OR. ('A' <= cl1 .AND. cl1 <= 'Z') .OR. ('a' <= cl1 .AND. cl1 <= 'z') ) THEN
            ptr%cleanm(ji:ji) = cl1   ! alphanumeric characters: keep it
         ELSE
            ptr%cleanm(ji:ji) = '_'   ! non-alphanumeric characters: replace it by '_'
         ENDIF
      END DO
      ptr%n8tnet(:)  = 0_8
      ptr%n8tfull(:) = 0_8
      ptr%ncalls(:)  = 0
      ptr%nrecurrent = 0
      ptr%ldone      = .TRUE.
      ptr%ncid       = -1
      ptr%s_parent   => NULL()
      ptr%s_prev     => NULL()
      ptr%s_next     => NULL()
      IF( ikeepnc > 0 ) THEN
         IF( numtime /= -1 ) THEN   ! create nc file only if we also opened the timing_nanuq.output file
            clname = 'timing_'//TRIM(ptr%cleanm)//TRIM(add_xxx())//'.nc'
            CALL nf90chk( NF90_CREATE( TRIM(clname), IOR( NF90_64BIT_OFFSET, NF90_CLOBBER ), ptr%ncid ) )
            CALL nf90chk( NF90_SET_FILL( ptr%ncid, NF90_NOFILL, ioldMode) )
            CALL nf90chk( NF90_DEF_DIM(  ptr%ncid, 'kt', NF90_UNLIMITED, idtime ) )
            CALL nf90chk( NF90_DEF_VAR(  ptr%ncid, 'timing_'//TRIM(ptr%cleanm), NF90_DOUBLE, (/ idtime /), ivid ) )
            CALL nf90chk( NF90_PUT_ATT(  ptr%ncid, ivid,  'name', 'elapse time spent on each call to '//TRIM(ptr%cleanm) ) )
            CALL nf90chk( NF90_PUT_ATT(  ptr%ncid, ivid, 'units', 's' ) )
            CALL nf90chk( NF90_ENDDEF(   ptr%ncid ) )
         ENDIF
         ALLOCATE( ptr%tkeep(ikeepnc) )
         ptr%nwrt   = 0
         ptr%nchunk = 0
      ENDIF

   END FUNCTION def_newlink


   SUBROUTINE switch_links(sd_root, sd_current)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE switch_links  ***
      !! ** Purpose :   switch current and current%s_next chain links
      !!                switch_links is called only if current%s_next is associated
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      TYPE(timer), POINTER, INTENT(inout) :: sd_current   ! current chain link
      !!----------------------------------------------------------------------
      !
      ! initial setup: current%prev <-> current <-> current%next <-> current%next%next
      ! --------------
      !
      !                      current
      !                         2
      !        current%prev 1       4 current%next%next
      !                         3
      !                    current%next
      !
      ! after switch: current%prev%prev <-> current%prev <-> current <-> current%next
      ! ------------
      !                      current
      !                         3
      !   current%prev%prev 1       4 current%next
      !                         2
      !                   current*prev
      !
      ! link current%s_prev and current%s_next
      IF( ASSOCIATED( sd_current%s_prev ) ) THEN
         sd_current%s_prev%s_next => sd_current%s_next   ! link forward current%s_prev to current%s_next
         sd_current%s_next%s_prev => sd_current%s_prev   ! link backward current%s_next to current%s_prev
      ELSE
         sd_current%s_next%s_prev => NULL()              ! current%next becomes the new root of the chain
         sd_root => sd_current%s_next
      ENDIF

      sd_current%s_prev => sd_current%s_next             ! link backward current to its new prev (former next)

      ! cannot yet modify sd_current%s_prev%s_next (which is also sd_current%s_next%s_next) as
      IF( ASSOCIATED( sd_current%s_next%s_next ) ) THEN
         sd_current%s_next => sd_current%s_next%s_next   ! link forward current to its new next (former next%next)
         sd_current%s_next%s_prev => sd_current          ! link backward new current next to current
      ELSE
         sd_current%s_next => NULL()   ! current is now at the end of the chain
      ENDIF
      ! we can now update sd_current%s_prev%s_next
      sd_current%s_prev%s_next => sd_current   ! link new current prev to current
      !
   END SUBROUTINE switch_links


   FUNCTION nblinks( sd_root ) RESULT(inb)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE nblinks  ***
      !! ** Purpose :   return the number of links in the chain link
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(in) :: sd_root      ! root chain link of the chain
      !
      TYPE(timer), POINTER :: s_wrk
      INTEGER :: inb
      !!----------------------------------------------------------------------
      inb = 0
      s_wrk => sd_root
      DO WHILE( ASSOCIATED(s_wrk) )
         inb = inb + 1
         s_wrk => s_wrk%s_next
      END DO
   END FUNCTION nblinks


   SUBROUTINE sort_chain( sd_root, kswitch )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE sort_chain  ***
      !! ** Purpose :   sort the chain link
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      INTEGER,              INTENT(in   ) :: kswitch      ! select the sorting variable
      !
      TYPE(timer), POINTER :: s_wrk
      LOGICAL              :: ll_ord, ll_test
      !!----------------------------------------------------------------------

      ll_ord = .FALSE.
      DO WHILE ( .NOT. ll_ord )
         ll_ord = .TRUE.
         s_wrk => sd_root
         DO WHILE ( ASSOCIATED( s_wrk%s_next ) )
            IF(     kswitch == jp_cname ) THEN
               ll_test = s_wrk%cname    < s_wrk%s_next%cname
            ELSEIF( kswitch == jp_tnet  ) THEN
               ll_test = s_wrk%tnet     < s_wrk%s_next%tnet
            ELSEIF( kswitch == jp_tavg  ) THEN
               ll_test = s_wrk%tnetavg  < s_wrk%s_next%tnetavg
            ELSEIF( kswitch == jp_tmin  ) THEN
               ll_test = s_wrk%tnetmin  < s_wrk%s_next%tnetmin
            ELSEIF( kswitch == jp_tmax  ) THEN
               ll_test = s_wrk%tnetmax  < s_wrk%s_next%tnetmax
            ELSEIF( kswitch == jp_tblc  ) THEN
               ll_test = s_wrk%tnetblc  < s_wrk%s_next%tnetblc
            ENDIF
            IF( ll_test ) THEN
               CALL switch_links(sd_root, s_wrk)
               ll_ord = .FALSE.
               CYCLE
            ENDIF
            IF( ASSOCIATED(s_wrk%s_next) ) s_wrk => s_wrk%s_next
         END DO
      END DO

   END SUBROUTINE sort_chain


   SUBROUTINE mpi_avgminmax( pval, pavg, pmin, pmax, ld_avg )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpi_avgminmax  ***
      !! ** Purpose :   get average, min, max over all MPI processes
      !!----------------------------------------------------------------------
      REAL(dp), INTENT(in   ) ::   pval              ! local value (of the current MPI process)
      REAL(dp), INTENT(  out) ::   pavg, pmin, pmax  ! mean, min, max and load balance (among MPI processes)
      LOGICAL , INTENT(in   ) ::   ld_avg
      !
      REAL(dp), DIMENSION(:  ), ALLOCATABLE ::  zallmpi
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  z2d
      INTEGER :: icode
      !!----------------------------------------------------------------------

      IF(ld_avg) THEN
#if ! defined key_mpi_off
         ALLOCATE(zallmpi(jpnij))
         CALL MPI_ALLGATHER(    pval, 1, MPI_DOUBLE_PRECISION,   &
            &                zallmpi, 1, MPI_DOUBLE_PRECISION, nmpicom, icode)
#endif
         pavg = SUM(    zallmpi ) / REAL(jpnij, dp)
         pmin = MINVAL( zallmpi )
         pmax = MAXVAL( zallmpi )
         DEALLOCATE(zallmpi)
      ELSE
         pavg = pval
         pmin = pval
         pmax = pval
      ENDIF

   END SUBROUTINE mpi_avgminmax


   SUBROUTINE mpi_t_avgminmax( sd_root, kn )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE   ***
      !! ** Purpose :
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_root      ! root chain link of the chain
      INTEGER             , INTENT(in   ) :: kn           ! time window index
      !
      TYPE(timer), POINTER :: s_wrk
      REAL(dp)                                ::  zjpnij_r
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::  zallmpi
      REAL(dp), DIMENSION(:,:  ), ALLOCATABLE ::  zlocal
      INTEGER :: ji, jn, ii
      INTEGER :: inb, ichksz, inmax
      INTEGER :: ins, ine, insz, inszold
      INTEGER :: icode, ierr
      !!----------------------------------------------------------------------
      zjpnij_r = 1._dp / REAL(jpnij, dp)

      ! reorder the chain list according to cname to make sure that each MPI process has the chain links in the same order
      CALL sort_chain( sd_root, jp_cname )

      inb = nblinks( sd_root )

      ! cut the work into chunks to avoid to allocate huge arrays when jpnij is very big...
      ichksz = nszmax / (2*jpnij)
      inmax = (inb-1) / ichksz + 1
      ins = 1 - ichksz
      inszold = -1
      DO jn = 1, inmax

         ins  = ins + ichksz
         ine  = MIN( inb, ins + ichksz - 1 )
         insz = ine - ins + 1

         IF( insz /= inszold ) THEN
            IF( ALLOCATED(zlocal) ) DEALLOCATE(zlocal)
            ALLOCATE( zlocal(insz,2) )
         ENDIF
         s_wrk => sd_root
         DO ji = 1, inb
            IF( ji >= ins ) THEN
               ii = ji - ins + 1
               zlocal(ii,1) = s_wrk%tnet
               zlocal(ii,2) = s_wrk%tfull
               IF( ji == ine )   EXIT
            ENDIF
            s_wrk => s_wrk%s_next
         END DO

         IF( insz /= inszold ) THEN
            IF( ALLOCATED(zallmpi) )   DEALLOCATE(zallmpi)
            IF( narea == 1 ) THEN
               ALLOCATE( zallmpi(insz,2,jpnij), STAT=ierr )
            ELSE
               ALLOCATE( zallmpi(   1,1,    1), STAT=ierr )   ! not used, allocate less memory
            ENDIF
            IF( ierr /= 0 )   CALL local_stop( 'cannot allocate zallmpi in mpi_t_avgminmax' )
         ENDIF
#if ! defined key_mpi_off
         CALL MPI_GATHER(  zlocal, insz*2, MPI_DOUBLE_PRECISION,   &
            &             zallmpi, insz*2, MPI_DOUBLE_PRECISION, 0, nmpicom, icode )
#endif
         IF( narea == 1 ) THEN
            s_wrk => sd_root
            DO ji = 1, inb
               IF( ji >= ins ) THEN
                  ii = ji - ins + 1
                  s_wrk%tnetavg  = SUM(    zallmpi(ii,1,:) ) * zjpnij_r
                  s_wrk%tnetmin  = MINVAL( zallmpi(ii,1,:) )
                  s_wrk%tnetmax  = MAXVAL( zallmpi(ii,1,:) )
                  s_wrk%tfullavg = SUM(    zallmpi(ii,2,:) ) * zjpnij_r
                  s_wrk%tfullmin = MINVAL( zallmpi(ii,2,:) )
                  s_wrk%tfullmax = MAXVAL( zallmpi(ii,2,:) )
                  s_wrk%tfullblc = s_wrk%tfullmax - s_wrk%tfullmin
                  s_wrk%tnetblc  = s_wrk%tnetmax  - s_wrk%tnetmin
                  IF( ji == ine )   EXIT
               ENDIF
               s_wrk => s_wrk%s_next
            END DO
         ENDIF

         inszold = insz
      END DO   ! jn, chunks

      DEALLOCATE(zlocal, zallmpi)

   END SUBROUTINE mpi_t_avgminmax


   SUBROUTINE write_tsum_allmpi( sd_root, kn )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE   ***
      !! ** Purpose :
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER      , INTENT(inout) :: sd_root      ! root chain link of the chain
      INTEGER                   , INTENT(in   ) :: kn           ! time window index
      !
      TYPE(timer), POINTER :: s_wrk
      INTEGER :: ioldMode, ierr, icode
      INTEGER :: ji, jn, jr, ii, ij, iii, inb
      INTEGER :: icuti, icutj, intfu, incid, ivid
      INTEGER :: ichksz, inmax, ins, ine, insz, inszold
      INTEGER,  DIMENSION(2  )                :: ilocal
      INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iallmpi
      INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: irank2d
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zallmpi
      REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: zlocal
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: z3d
      CHARACTER(LEN=9 ) :: cl1st, clend
      CHARACTER(LEN=64) :: clname
      !!----------------------------------------------------------------------
      inb = nblinks( sd_root )
      !
      IF( narea == 1 ) THEN
         WRITE(cl1st, '(i9)') n1st(kn)   ;   cl1st = ADJUSTL(cl1st)
         WRITE(clend, '(i9)') nend(kn)   ;   clend = ADJUSTL(clend)
         clname = 'timing_tsum_allmpi_t'//TRIM(cl1st)//'_t'//TRIM(clend)//'.nc'
         CALL nf90chk( NF90_CREATE( clname, IOR( NF90_64BIT_OFFSET, NF90_CLOBBER ), incid ) )
         CALL nf90chk( NF90_SET_FILL( incid, NF90_NOFILL, ioldMode ) )
         CALL nf90chk( NF90_DEF_DIM( incid, 'mpicut_i', jpni, icuti ) )
         CALL nf90chk( NF90_DEF_DIM( incid, 'mpicut_j', jpnj, icutj ) )
         CALL nf90chk( NF90_DEF_DIM( incid, 'net_full',    2, intfu ) )
         !
         CALL nf90chk( NF90_DEF_VAR( incid, 'mpirank', NF90_INT, (/ icuti, icutj /), ivid ) )
         CALL nf90chk( NF90_PUT_ATT( incid, ivid,       'name', 'MPI rank' ) )
         CALL nf90chk( NF90_PUT_ATT( incid, ivid, '_FillValue', -1_i4 ) )
         !
         CALL sort_chain( sd_root, jp_tavg )   ! sort accoring to tnetavg
         s_wrk => sd_root
         DO WHILE( ASSOCIATED(s_wrk) )
            IF( s_wrk%tfullmax > 0._dp ) THEN   ! can happen if not present in the timing time window
               CALL nf90chk( NF90_DEF_VAR( incid, TRIM(s_wrk%cleanm), NF90_DOUBLE, (/ icuti, icutj, intfu /), ivid ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,        'name', 'deviation from MPI average of the total elapse time spent on '//TRIM(s_wrk%cleanm) ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,  '_FillValue', 0._dp ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,       'units', '%' ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid, 'time_window', (/ n1st(kn), nend(kn) /) ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,     'mpi_avg', (/    0._dp,    0._dp /) ) )
            ENDIF
            s_wrk => s_wrk%s_next
         END DO
         !
         CALL nf90chk( NF90_ENDDEF(incid) )
      ENDIF

      ! reorder the chain list according to cname to make sure that each MPI process has the chain links in the same order
      CALL sort_chain( sd_root, jp_cname )

      ! fill mpi rank 2D array
      ilocal = (/ nimpi, njmpi /)
      IF( narea == 1 ) THEN
         ALLOCATE( iallmpi(2,jpnij), STAT = ierr )   ! can be huge if jpnij is big...
      ELSE
         ALLOCATE( iallmpi(1,    1), STAT = ierr )   ! not used, allocate less memory
      ENDIF
      IF( ierr /= 0 )   CALL local_stop( 'cannot allocate iallmpi in write_tsum_allmpi' )
#if ! defined key_mpi_off
      CALL MPI_GATHER(  ilocal, 2, MPI_INTEGER,   &
         &             iallmpi, 2, MPI_INTEGER, 0, nmpicom, icode )
#endif
      IF( narea == 1 ) THEN
         ALLOCATE( irank2d(jpni,jpnj), STAT = ierr )
         irank2d(:,:) = -1
         DO jr = 1, jpnij
            ii = iallmpi(1,jr)
            ij = iallmpi(2,jr)
            irank2d(ii,ij) = jr - 1
         END DO
         CALL nf90chk( NF90_INQ_VARID(incid, 'mpirank', ivid) )
         CALL nf90chk( NF90_PUT_VAR(  incid, ivid, irank2d ) )
         DEALLOCATE(irank2d)   ! free memory as soon as possible as this array can be very big...
      ENDIF

      IF( narea == 1 ) THEN
         ALLOCATE( z3d(jpni,jpnj,2) )
         z3d(:,:,:) = 0._dp
      ENDIF

      ! cut the work into chunks to avoid to allocate huge arrays when jpnij is very big...
      ichksz = nszmax / (2*jpnij)
      inmax = (inb-1) / ichksz + 1
      ins = 1 - ichksz
      inszold = -1
      DO jn = 1, inmax

         ins  = ins + ichksz
         ine  = MIN( inb, ins + ichksz - 1 )
         insz = ine - ins + 1

         IF( insz /= inszold ) THEN
            IF( ALLOCATED(zlocal) ) DEALLOCATE(zlocal)
            ALLOCATE( zlocal(insz,2) )
         ENDIF
         s_wrk => sd_root
         DO ji = 1, inb
            IF( ji >= ins ) THEN
               ii = ji - ins + 1
               zlocal(ii,1) = s_wrk%tnet
               zlocal(ii,2) = s_wrk%tfull
               IF( ji == ine )   EXIT
            ENDIF
            s_wrk => s_wrk%s_next
         END DO

         IF( insz /= inszold ) THEN
            IF( ALLOCATED(zallmpi) )   DEALLOCATE(zallmpi)
            IF( narea == 1 ) THEN
               ALLOCATE( zallmpi(insz,2,jpnij), STAT=ierr )
            ELSE
               ALLOCATE( zallmpi(   1,1,    1), STAT=ierr )   ! not used, allocate less memory
            ENDIF
            IF( ierr /= 0 )   CALL local_stop( 'cannot allocate zallmpi in write_tsum_allmpi' )
         ENDIF
#if ! defined key_mpi_off
         CALL MPI_GATHER(  zlocal, insz*2, MPI_DOUBLE_PRECISION,   &
            &             zallmpi, insz*2, MPI_DOUBLE_PRECISION, 0, nmpicom, icode )
#endif
         IF( narea == 1 ) THEN
            s_wrk => sd_root
            DO ji = 1, inb
               IF( ji >= ins ) THEN
                  IF( s_wrk%tfullmax > 0._dp ) THEN   ! can not happen if not present in the timing time window
                     iii = ji - ins + 1
                     DO jr = 1, jpnij
                        ii = iallmpi(1,jr)
                        ij = iallmpi(2,jr)
                        z3d(ii,ij,1) = ( zallmpi(iii,1,jr) / s_wrk%tnetavg  - 1._dp ) * 100._dp
                        z3d(ii,ij,2) = ( zallmpi(iii,2,jr) / s_wrk%tfullavg - 1._dp ) * 100._dp
                     END DO
                     CALL nf90chk( NF90_INQ_VARID(incid, TRIM(s_wrk%cleanm), ivid) )
                     CALL nf90chk( NF90_PUT_VAR(  incid, ivid, z3d ) )
                     CALL nf90chk( NF90_PUT_ATT(  incid, ivid, 'mpi_avg', (/ s_wrk%tnetavg, s_wrk%tfullavg /) ) )
                  ENDIF
                  IF( ji == ine )   EXIT
               ENDIF
               s_wrk => s_wrk%s_next
            END DO
            IF( ine == inb )   CALL nf90chk( NF90_CLOSE(incid) )
         ENDIF

         inszold = insz
      END DO   ! jn, chunks

      IF( narea == 1 )     DEALLOCATE(z3d)
      DEALLOCATE(iallmpi, zlocal, zallmpi)

   END SUBROUTINE write_tsum_allmpi


   SUBROUTINE write_allts_nc( sd_root )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE   ***
      !! ** Purpose :
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER      , INTENT(inout) :: sd_root      ! root chain link of the chain
      !
      TYPE(timer), POINTER :: s_wrk
      INTEGER :: ioldMode, ierr, icode
      INTEGER :: ji, ii, ij, jn, jr
      INTEGER :: isz, ishft, iks, ike
      INTEGER :: idtime, icuti, icutj, intfu, incid, ivid
      INTEGER :: ichksz, inmax, ins, ine, insz, inszold
      INTEGER,  DIMENSION(2  )                :: ilocal
      INTEGER,  DIMENSION(3  )                :: ichklocal
      INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iallmpi, ichkmpi
      INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: irank2d
      REAL(dp), DIMENSION(:    ), ALLOCATABLE :: zavg
      REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: zallmpi
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: z3d
      CHARACTER(LEN=64) :: clname
      !!----------------------------------------------------------------------
      inszold = -1

      ! reorder the chain list according to cname to make sure that each MPI process has the chain links in the same order
      CALL sort_chain( sd_root, jp_cname )
      !
      s_wrk => sd_root
      DO WHILE( ASSOCIATED(s_wrk) )

         IF( s_wrk%ncid /= -1 ) THEN
            IF( s_wrk%nwrt > 0 ) THEN   ! remaining data to be written
               CALL nf90chk( NF90_INQ_VARID(s_wrk%ncid, 'timing_'//TRIM(s_wrk%cleanm), ivid) )
               CALL nf90chk( NF90_PUT_VAR(  s_wrk%ncid, ivid, s_wrk%tkeep(1:s_wrk%nwrt),  &
                  &                         (/ s_wrk%nchunk * SIZE(s_wrk%tkeep) + 1 /), (/s_wrk%nwrt/) ) )
            ENDIF
            CALL nf90chk( NF90_CLOSE( s_wrk%ncid ) )   ! close NetCDF timing time series
         ENDIF

         IF( ALLOCATED(s_wrk%tkeep) .AND. jpnij > 1 ) THEN   ! create a NetCDF with the last part of all MPI timing time series

            isz = SIZE(s_wrk%tkeep)
            ike = s_wrk%nchunk * isz + s_wrk%nwrt
            ! use a maximum of available data
            IF( s_wrk%nwrt > 0 .AND. s_wrk%nchunk > 0 ) THEN
               ishft = isz - s_wrk%nwrt   ! > 0
               s_wrk%tkeep = CSHIFT(s_wrk%tkeep, -ishft )   ! shift is < 0
               s_wrk%nwrt = isz
            ELSEIF( s_wrk%nwrt == 0 .AND. s_wrk%nchunk > 0 ) THEN
               s_wrk%nwrt = isz
            ELSEIF( s_wrk%nwrt == 0 .AND. s_wrk%nchunk == 0 ) THEN
               EXIT
            ENDIF
            iks = ike - s_wrk%nwrt + 1

            ! some (useless) checks
            IF( narea == 1 ) THEN
               ALLOCATE( ichkmpi(3,jpnij), STAT = ierr )   ! can be huge if jpnij is big...
            ELSE
               ALLOCATE( ichkmpi(1,    1), STAT = ierr )   ! not used, allocate less memory
            ENDIF
            IF( ierr /= 0 )   CALL local_stop( 'cannot allocate ichkmpi in write_allts_nc' )
            ichklocal = (/ iks, ike, s_wrk%nwrt /)
#if ! defined key_mpi_off
            CALL MPI_GATHER(  ichklocal, 3, MPI_INTEGER,   &
               &                ichkmpi, 3, MPI_INTEGER, 0, nmpicom, icode )
#endif
            IF( narea == 1 ) THEN
               IF( COUNT( SUM( ichkmpi, dim=2 ) /= ichklocal*jpnij ) /= 0 )   &
                  CALL local_stop( TRIM(s_wrk%cleanm)//' the number of measurements is not the same on all processes' )
            ENDIF
            DEALLOCATE(ichkmpi)

            IF( narea == 1 ) THEN
               clname = 'timing_ts_allmpi_'//TRIM(s_wrk%cleanm)//'.nc'
               CALL nf90chk( NF90_CREATE( clname, IOR( NF90_64BIT_OFFSET, NF90_CLOBBER ), incid ) )
               CALL nf90chk( NF90_SET_FILL( incid, NF90_NOFILL, ioldMode) )
               CALL nf90chk( NF90_DEF_DIM( incid, 'kt', NF90_UNLIMITED, idtime ) )
               CALL nf90chk( NF90_DEF_DIM( incid, 'mpicut_i', jpni, icuti ) )
               CALL nf90chk( NF90_DEF_DIM( incid, 'mpicut_j', jpnj, icutj ) )
               !
               CALL nf90chk( NF90_DEF_VAR( incid,      'kt', NF90_INT, (/ idtime       /), ivid ) )
               CALL nf90chk( NF90_DEF_VAR( incid, 'mpirank', NF90_INT, (/ icuti, icutj /), ivid ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,       'name', 'MPI rank' ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid, '_FillValue', -1_i4 ) )
               !
               CALL nf90chk( NF90_DEF_VAR( incid, 'timing_mpiavg_'//TRIM(s_wrk%cleanm), NF90_DOUBLE, (/ idtime /), ivid ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,       'name', 'MPI averaged elapse of the time spent on each call to '//TRIM(s_wrk%cleanm) ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,      'units', 's' ) )
               !
               CALL nf90chk( NF90_DEF_VAR( incid, 'timing_mpipc_'//TRIM(s_wrk%cleanm), NF90_DOUBLE, (/ icuti, icutj, idtime /), ivid ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,       'name', 'deviation from MPI average on each call to '//TRIM(s_wrk%cleanm) ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid, '_FillValue', 0._dp ) )
               CALL nf90chk( NF90_PUT_ATT( incid, ivid,      'units', '%' ) )
               CALL nf90chk( NF90_ENDDEF(incid) )

               CALL nf90chk( NF90_INQ_VARID(incid, 'kt', ivid) )
               CALL nf90chk( NF90_PUT_VAR(  incid, ivid, (/ (ji, ji=iks,ike) /) ) )
            ENDIF

            ! fill mpi rank 2D array
            IF( .NOT. ALLOCATED(iallmpi) ) THEN
               ilocal=(/ nimpi, njmpi /)
               IF( narea == 1 ) THEN
                  ALLOCATE( iallmpi(2,jpnij), STAT = ierr )   ! can be huge if jpnij is big...
               ELSE
                  ALLOCATE( iallmpi(1,    1), STAT = ierr )   ! not used, allocate less memory
               ENDIF
               IF( ierr /= 0 )   CALL local_stop( 'cannot allocate iallmpi in write_allts_nc' )
#if ! defined key_mpi_off
               CALL MPI_GATHER(  ilocal, 2, MPI_INTEGER,   &
                  &             iallmpi, 2, MPI_INTEGER, 0, nmpicom, icode )
#endif
            ENDIF
            IF( narea == 1 ) THEN
               IF( .NOT. ALLOCATED(irank2d) ) THEN
                  ALLOCATE( irank2d(jpni,jpnj), STAT = ierr )
                  irank2d(:,:) = -1
                  DO jr = 1, jpnij
                     ii = iallmpi(1,jr)
                     ij = iallmpi(2,jr)
                     irank2d(ii,ij) = jr - 1
                  END DO
               ENDIF
               CALL nf90chk( NF90_INQ_VARID(incid, 'mpirank', ivid) )
               CALL nf90chk( NF90_PUT_VAR(  incid, ivid, irank2d ) )
            ENDIF

            ! cut the work into chunks to avoid to allocate huge arrays when jpnij is very big...
            ichksz = nszmax / jpnij
            inmax = (s_wrk%nwrt-1) / ichksz + 1
            ins = 1 - ichksz
            DO jn = 1, inmax

               ins  = ins + ichksz
               ine  = MIN( s_wrk%nwrt, ins + ichksz - 1 )
               insz = ine - ins + 1

               IF( insz /= inszold ) THEN
                  IF( ALLOCATED(zallmpi) )   DEALLOCATE(zallmpi)
                  IF( narea == 1 ) THEN
                     ALLOCATE( zallmpi(insz,jpnij), STAT=ierr )
                  ELSE
                     ALLOCATE( zallmpi(   1,    1), STAT=ierr )   ! not used, allocate less memory
                  ENDIF
                  IF( ierr /= 0 )   CALL local_stop( 'cannot allocate zallmpi in write_allts_nc' )
               ENDIF
#if ! defined key_mpi_off
               CALL MPI_GATHER(  s_wrk%tkeep(ins:ine), insz, MPI_DOUBLE_PRECISION,   &
                  &                           zallmpi, insz, MPI_DOUBLE_PRECISION, 0, nmpicom, icode )
#endif
               IF( narea == 1 ) THEN
                  IF( insz /= inszold ) THEN
                     IF( ALLOCATED(z3d) )   DEALLOCATE(z3d, zavg)
                     ALLOCATE( z3d(jpni,jpnj,insz), zavg(insz) )
                     z3d(:,:,:) = 0._dp
                  ENDIF
                  zavg(:) = SUM( zallmpi(:,:), dim = 2 ) / REAL(jpnij, dp)
                  CALL nf90chk( NF90_INQ_VARID(incid, 'timing_mpiavg_'//TRIM(s_wrk%cleanm), ivid) )
                  CALL nf90chk( NF90_PUT_VAR(  incid, ivid, zavg, (/ins/), (/insz/) ) )
                  DO jr = 1, jpnij
                     ii = iallmpi(1,jr)
                     ij = iallmpi(2,jr)
                     z3d(ii,ij,:) = ( zallmpi(:,jr) / zavg(:) - 1._dp ) * 100._dp
                  END DO
                  CALL nf90chk( NF90_INQ_VARID(incid, 'timing_mpipc_'//TRIM(s_wrk%cleanm), ivid) )
                  CALL nf90chk( NF90_PUT_VAR(  incid, ivid, z3d, (/1,1,ins/), (/jpni,jpnj,insz/) ) )
                  IF( ine == s_wrk%nwrt ) CALL nf90chk( NF90_CLOSE(incid) )
               ENDIF

               inszold = insz
            END DO   ! jn, chunks

         ENDIF   ! ALLOCATED(s_wrk%tkeep) .AND. jpnij > 1
         s_wrk => s_wrk%s_next
      END DO   ! WHILE( ASSOCIATED(s_wrk) )

      IF( ALLOCATED(irank2d) )   DEALLOCATE(irank2d, z3d, zavg)
      IF( ALLOCATED(iallmpi) )   DEALLOCATE(iallmpi, zallmpi)

   END SUBROUTINE write_allts_nc


   SUBROUTINE write_header(cdname)
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
      !
      CHARACTER(len=64 ) :: clfmt
      CHARACTER(LEN=128) :: cldash
      !!----------------------------------------------------------------------
      !
      CALL write_separator()
      WRITE(numtime,'(a)') ' '//TRIM(cdname)
      WRITE(clfmt, "('(x,',i2,'(''-''))')") LEN_TRIM(cdname)
      WRITE(cldash, clfmt)
      WRITE(numtime,'(a)') TRIM(cldash)
      WRITE(numtime,*)
      !
   END SUBROUTINE write_header


   SUBROUTINE write_bigheader(cdname)
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
      !
      INTEGER :: inb1, inb2
      CHARACTER(len=64 ) :: clfmt
      CHARACTER(LEN=128) :: clequal, cldash, cline
      !!----------------------------------------------------------------------
      !
      WRITE(clfmt, "('(3x,',i2,'(''-''))')") nbequal
      WRITE(cldash , clfmt)
      WRITE(clfmt, "('(3x,',i2,'(''=''))')") nbequal
      WRITE(clequal, clfmt)

      inb1 = ( nbequal - LEN_TRIM(cdname) - 2 + 1 ) / 2
      inb2 = ( nbequal - LEN_TRIM(cdname) - 2     ) / 2
      WRITE(clfmt, "('(3x,''|'',', i2,'('' ''),a,', i2,'('' ''),''|'')')") inb1, inb2
      WRITE(cline,clfmt) TRIM(cdname)

      WRITE(numtime,*)
      WRITE(numtime,*)
      WRITE(numtime,'(a)') TRIM(cldash)
      WRITE(numtime,'(a)') TRIM(clequal)
      WRITE(numtime,'(a)') TRIM(cline)
      WRITE(numtime,'(a)') TRIM(clequal)
      WRITE(numtime,'(a)') TRIM(cldash)
      WRITE(numtime,*)
      !
   END SUBROUTINE write_bigheader


   SUBROUTINE write_separator()
      !!----------------------------------------------------------------------
      CHARACTER(len=64 ) ::   clfmt
      CHARACTER(LEN=128) :: clequal
      !!----------------------------------------------------------------------
      !
      WRITE(clfmt, "('(3x,',i2,'(''=''))')") nbequal
      WRITE(clequal, clfmt)

      WRITE(numtime,*)
      WRITE(numtime,'(a)') TRIM(clequal)
      WRITE(numtime,*)
      WRITE(numtime,*)
      !
   END SUBROUTINE write_separator


   FUNCTION add_xxx()   RESULT( cldgt )
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION add_xxx  ***
      !! ** Purpose :   add '_xxx' if narea > 1
      !!----------------------------------------------------------------------
      INTEGER :: idg
      CHARACTER(LEN=10) :: clfmt, cldgt
      !!----------------------------------------------------------------------

      IF( narea > 1 ) THEN   ! add _xxxx
         idg = MAX( INT(LOG10(REAL(MAX(1,jpnij)))) + 1, 4 )
         WRITE(clfmt, "('( a, i', i1, '.', i1, ')')") idg, idg   ! '(a,ix.x)'
         WRITE(cldgt, clfmt) '_', narea
      ELSE
         cldgt = ''
      ENDIF

   END FUNCTION add_xxx


   SUBROUTINE local_stop( cdname )
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION local_stop  ***
      !! ** Purpose :   stop the model
      !!                we cannot USE ctl_stop as we cannot USE lib_mpp
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
      INTEGER :: icode
      !!----------------------------------------------------------------------
      !
      WRITE(*,'(a,i0)') ' STOP from timing: '//TRIM(cdname)//', MPI rank: ', narea-1
#if ! defined key_mpi_off
      CALL mpi_abort( MPI_COMM_WORLD, 123, icode )
#else
      STOP 123
#endif
   END SUBROUTINE local_stop


   SUBROUTINE nf90chk( kstatus )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE nf90chk  ***
      !!
      !! ** Purpose :   check nf90 errors
      !!--------------------------------------------------------------------
      INTEGER,          INTENT(in) :: kstatus
      !---------------------------------------------------------------------
      IF(kstatus /= nf90_NOERR)   CALL local_stop( NF90_STRERROR(kstatus) )
   END SUBROUTINE nf90chk


   SUBROUTINE gnuplot_statplot( sd_root, kn )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE gnuplot_statplot  ***
      !! ** Purpose :   Perform gnuplot statplot/plots
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(in) :: sd_root      ! root chain link of the chain
      INTEGER             , INTENT(in) :: kn           ! time window index
      !
      TYPE(timer), POINTER :: s_wrk
      INTEGER :: icmdstat
      CHARACTER(LEN=512) :: clcmd
      CHARACTER(LEN=256) :: clfile
      CHARACTER(LEN=9  ) :: cl1st, clend, clrnk
      !!----------------------------------------------------------------------

      s_wrk => sd_root
      DO WHILE( ASSOCIATED(s_wrk) )

         IF( s_wrk%ncid /= -1 ) THEN

            CALL write_header('Stats report on the timing of '//TRIM(s_wrk%cleanm)//' (in microsecond):')

            INQUIRE(unit = numtime, NAME = clfile)   ! get numtime file name (i.e. timing_nanuq.output)
            CLOSE(numtime)                           ! close it as EXECUTE_COMMAND_LINE will write in this file

            ! write gnoplot stats in clfile (i.e. timing_nanuq.output)
            WRITE(clrnk, '(i9)') narea-1    ;   clrnk = ADJUSTL(clrnk)
            WRITE(cl1st, '(i9)') n1st(kn)   ;   cl1st = ADJUSTL(cl1st)
            WRITE(clend, '(i9)') nend(kn)   ;   clend = ADJUSTL(clend)

            clcmd = './timing_nanuq_gnuplot.sh --vname timing_'//TRIM(s_wrk%cleanm)//' --rnk '//TRIM(clrnk)//   &
               &    ' --png --t1st '//TRIM(cl1st)//' --tend '//TRIM(clend)//' 2>> '//TRIM(clfile)
            CALL EXECUTE_COMMAND_LINE(TRIM(clcmd), CMDSTAT = icmdstat)

            OPEN(NEWUNIT = numtime, FILE = TRIM(clfile), POSITION = "append", STATUS = "old", ACTION = "write")   ! reopen it
            IF( icmdstat > 0 )   RETURN   ! EXECUTE_COMMAND_LINE not properly working

         ENDIF
         s_wrk => s_wrk%s_next
      END DO

   END SUBROUTINE gnuplot_statplot


   SUBROUTINE write_gnuplotscript( )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE write_gnuplotscript  ***
      !! ** Purpose :   Write gnuplot statplot/plots
      !!----------------------------------------------------------------------
      INTEGER :: inumsh, icmdstat
      !!----------------------------------------------------------------------

      OPEN(NEWUNIT = inumsh, FILE = 'timing_nanuq_gnuplot.sh', POSITION = 'rewind', STATUS = 'replace', ACTION = 'write')
      WRITE(inumsh,'(a)') '#!/bin/bash'
      WRITE(inumsh,'(a)') 'set -u'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'vname=timing_step'
      WRITE(inumsh,'(a)') 'doplot=1'
      WRITE(inumsh,'(a)') 'dopng=0'
      WRITE(inumsh,'(a)') 't1st=1    ; tend=0'
      WRITE(inumsh,'(a)') 'vminup=-1 ; vmaxup=-1'
      WRITE(inumsh,'(a)') 'vmindn=-1 ; vmaxdn=-1'
      WRITE(inumsh,'(a)') 'dinput="not defined"'
      WRITE(inumsh,'(a)') 'finput="not defined"'
      WRITE(inumsh,'(a)') 'rnk=0'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'while [ $# -gt 0 ]'
      WRITE(inumsh,'(a)') 'do'
      WRITE(inumsh,'(a)') "    case $( echo $1 | tr '[:upper:]' '[:lower:]' ) in"
      WRITE(inumsh,'(a)') '     -h|--help)'
      WRITE(inumsh,'(a)') '         echo "   Purpose:"'
      WRITE(inumsh,'(a)') '         echo "       Print statistics and plot of NEMO stats time series in micro-second"'
      WRITE(inumsh,'(a)') '         echo "   Usage:"'
      WRITE(inumsh,'(a)') '         echo "       ./timing_nanuq_gnuplot.sh"'
      WRITE(inumsh,'(a)') '         echo "   Options:"'
      WRITE(inumsh,'(a)') '         echo "       --png               # save plots in png"'
      WRITE(inumsh,'(a)') '         echo "       --noplot            # print only the statistics (no plot)"'
      WRITE(inumsh,'(a)') '         echo "       --vname vname       # NetCDF variable name (timing_step by default)"'
      WRITE(inumsh,'(a)') '         echo "       --dir dir_name      # dir_name: directory where is located (./ by default)"'
      WRITE(inumsh,'(a)') '         echo "       --finput file_name  # file_name: NetCDF input file (timing_step.nc by default)"'
      WRITE(inumsh,'(a)') '         echo "       --rnk               # MPI rank who wrote the NetCDF input file (0 by default)"'
      WRITE(inumsh,'(a)') '         echo "       --t1st time1        # time1: first time step to select ($t1st by default)"'
      WRITE(inumsh,'(a)') '         echo "       --tend time2        # time2: last time step to select (last by default)"'
      WRITE(inumsh,'(a)') '         echo "       --rgup vmin vmax    # [vmin, vmax]: upper plot range ([0, max] by default)"'
      WRITE(inumsh,'(a)') '         echo "       --rgdn vmin vmax    # [vmin, vmax]: lower plot range"'
      WRITE(inumsh,'(a)') '         exit 0 ;;'
      WRITE(inumsh,'(a)') '     --png)     dopng=1 ;;'
      WRITE(inumsh,'(a)') '     --noplot)  doplot=0 ;;'
      WRITE(inumsh,'(a)') '     --vname)   vname=${2}  ; shift ;;'
      WRITE(inumsh,'(a)') '     --dir)     dinput=${2} ; shift ;;'
      WRITE(inumsh,'(a)') '     --finput)  finput=${2} ; shift ;;'
      WRITE(inumsh,'(a)') '     --rnk)     rnk=${2}    ; shift ;;'
      WRITE(inumsh,'(a)') '     --t1st)    t1st=${2}   ; shift ;;'
      WRITE(inumsh,'(a)') '     --tend)    tend=${2}   ; shift ;;'
      WRITE(inumsh,'(a)') '     --rgup)    vminup=${2} ; shift ; vmaxup=${2} ; shift ;;'
      WRITE(inumsh,'(a)') '     --rgdn)    vmindn=${2} ; shift ; vmaxdn=${2} ; shift ;;'
      WRITE(inumsh,'(a)') '    esac'
      WRITE(inumsh,'(a)') '    shift'
      WRITE(inumsh,'(a)') 'done'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') '[ "$finput" = "not defined" ] && finput=${vname}.nc'
      WRITE(inumsh,'(a)') 'if [ "$dinput" != "not defined" ] && finput=$dinput/$finput'
      WRITE(inumsh,'(a)') 'then'
      WRITE(inumsh,'(a)') '   [ ! -d $dinput ] && echo "ERROR: $dinput is not a directory" && exit 1'
      WRITE(inumsh,'(a)') '   finput=$dinput/$finput'
      WRITE(inumsh,'(a)') 'fi'
      WRITE(inumsh,'(a)') 'if [ $rnk -gt 0 ]'
      WRITE(inumsh,'(a)') 'then'
      WRITE(inumsh,'(a)') '   finput=$( ls -1 ${finput/.nc/_*.nc} | grep -e "_0*${rnk}.nc$" )'
      WRITE(inumsh,'(a)') '   rnk=$( echo $finput | sed -e "s/.*_/_/" -e "s/\.nc//" )'
      WRITE(inumsh,'(a)') 'else'
      WRITE(inumsh,'(a)') '   rnk=""'
      WRITE(inumsh,'(a)') 'fi'
      WRITE(inumsh,'(a)') '[ ! -f $finput ] && echo "ERROR: $finput not found" && exit 1'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'ok=$( ncdump -h $finput | grep -c "double *${vname}(" )'
      WRITE(inumsh,'(a)') '[ $ok -eq 0 ] && echo "ERROR: $vname not found in $finput" && exit 1'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') '[ $tend -eq 0 ] && tend=$( ncdump -h $finput | grep UNLIMITED | sed -e "s/[^0-9]*//g" )'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'which gnuplot &> /dev/null'
      WRITE(inumsh,'(a)') '[ $? -ne 0 ] && echo "$( basename $0 ) uses gnuplot which was not found" && exit 1'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'ncdump -f f -v $vname $finput | sed -n -e "/\/\/ *${vname}(${t1st})/,/\/\/ *${vname}(${tend})/p" | sed -e "s/[;,] .*//" -e "s/.* /scale=6 ; 1000000. * /" | bc > timing_gnuplot.$$'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'gnuplot -persist << EOF'
      WRITE(inumsh,'(a)') 'dopng = $dopng'
      WRITE(inumsh,'(a)') 'doplot = $doplot'
      WRITE(inumsh,'(a)') 'vminup = $vminup ; vmaxup = $vmaxup'
      WRITE(inumsh,'(a)') 'vmindn = $vmindn ; vmaxdn = $vmaxdn'
      WRITE(inumsh,'(a)') 'if ( dopng == 0 ) {'
      WRITE(inumsh,'(a)') '   set terminal x11 size 800,800'
      WRITE(inumsh,'(a)') '} else {'
      WRITE(inumsh,'(a)') '   set terminal png size 800,800'
      WRITE(inumsh,'(a)') '   set output "${finput/${rnk}.nc/_t${t1st}_t${tend}${rnk}.png}"'
      WRITE(inumsh,'(a)') '}'
      WRITE(inumsh,'(a)') 'stats "timing_gnuplot.$$" name "ST"'
      WRITE(inumsh,'(a)') 'if ( doplot == 1 ) {'
      WRITE(inumsh,'(a)') '   mn = ST_mean ; md = ST_median ; std = ST_stddev'
      WRITE(inumsh,'(a)') '   if ( vminup == -1 ) { ; vminup = ST_min ; }'
      WRITE(inumsh,'(a)') '   if ( vmaxup == -1 ) { ; vmaxup = ST_max ; }'
      WRITE(inumsh,'(a)') '   lo = ST_lo_quartile ; up = ST_up_quartile ; iqr = up-lo'
      WRITE(inumsh,'(a)') '   set xrange [0:ST_records]'
      WRITE(inumsh,'(a)') '   set xlabel "${vname//_/ } number"'
      WRITE(inumsh,'(a)') '   set ylabel " ${vname//_/ } elapsed time (microsecond)"'
      WRITE(inumsh,'(a)') '   set multiplot layout 2,1'
      WRITE(inumsh,'(a)') '   set yrange [vminup:vmaxup]'
      WRITE(inumsh,'(a)') '   set title sprintf("FULL RANGE: ${vname//_/ } (microsecond), mean = %f, median = %f", mn, md)'
      WRITE(inumsh,'(a)') '   plot "timing_gnuplot.$$" notitle, mn title "Mean" lw 2, md title "Median" lw 2, mn-std title "Mean-StdDev" lw 2, mn+std title "Mean+StdDev" lw 2'
      WRITE(inumsh,'(a)') '   if ( vmindn == -1 ) {'
      WRITE(inumsh,'(a)') '      vmindn = lo - 1.5*iqr'
      WRITE(inumsh,'(a)') '      if ( vmindn < 0      ) { ; vmindn = 0 ; }'
      WRITE(inumsh,'(a)') '   }'
      WRITE(inumsh,'(a)') '   if ( vmaxdn == -1 ) {'
      WRITE(inumsh,'(a)') '      vmaxdn = up + 1.5*iqr'
      WRITE(inumsh,'(a)') '      if ( vmaxdn > ST_max ) { ; vmaxdn = ST_max ; }'
      WRITE(inumsh,'(a)') '   }'
      WRITE(inumsh,'(a)') '   set yrange [vmindn:vmaxdn]'
      WRITE(inumsh,'(a)') '   set title sprintf("ZOOM: ${vname//_/ } (microsecond), mean = %f, median = %f", mn, md)'
      WRITE(inumsh,'(a)') '   plot "timing_gnuplot.$$" notitle, mn title "Mean" lw 2, md title "Median" lw 2, lo title "1st Quatile" lw 2, up title "3rd Quartile" lw 2'
      WRITE(inumsh,'(a)') '}'
      WRITE(inumsh,'(a)') 'EOF'
      WRITE(inumsh,'(a)') '#'
      WRITE(inumsh,'(a)') 'rm -f timing_gnuplot.$$ &>/dev/null'
      CLOSE(inumsh)

      CALL EXECUTE_COMMAND_LINE('chmod u+x timing_nanuq_gnuplot.sh', CMDSTAT = icmdstat)   ! must use icmdstat to avoid stop if error

   END SUBROUTINE write_gnuplotscript




   !!=====================================================================
END MODULE timing
