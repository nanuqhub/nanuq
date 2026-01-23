MODULE icedyn_adv_pra_d
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_pra_d   ***
   !!   sea-ice : advection => Prather scheme
   !!
   !!        => advects ice damage  @T or @F -points
   !!
   !!======================================================================
   !! History :       !  2008-03  (M. Vancoppenolle) original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_pra_d : advection of sea ice using Prather scheme
   !!   adv_x_2d, adv_y_2d    : Prather scheme applied in i- and j-direction, resp.
   !!   adv_pra_d_init    : initialisation of the Prather scheme
   !!   adv_pra_d_rst     : read/write Prather field in ice restart file, or initialized to zero
   !!----------------------------------------------------------------------
   USE par_ice, ONLY: rDt_ice, ln_adv_pra
   USE iom            ! I/O manager library
   !
   USE icedyn_adv_pra_adv
   !
# if defined _OPENACC
   USE lbclnk_gpu
# else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
# endif
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv_pra_d   ! called by icedyn_adv
   PUBLIC   adv_pra_d_init      ! called by icedyn_adv

   ! Moments for advection
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: sx1md_t, sy1md_t, sxx1md_t, syy1md_t, sxy1md_t ! moments for `1-damage` @T
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: sx1md_f, sy1md_f, sxx1md_f, syy1md_f, sxy1md_f ! moments for `1-damage` @F

   !!----------------------------------------------------------------------
   !! NANUQ_beta
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ice_dyn_adv_pra_d( kt, cgt, pe1e2, p1_e1e2, pmsk, pUx, pVx,  p1md )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_pra_d  **
      !!
      !! ** BRODEAU, 2022 => Advect `damage` at T- & F- points
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!                but also their quadratic forms. The method preserves
      !!                tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in   ) :: kt         ! current time step
      CHARACTER(len=1),             INTENT(in   ) :: cgt        ! 'T' or 'F' points ???
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pe1e2, p1_e1e2 ! at the relevant point (T or F) ! ALREADY IN KM^2 !!!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pmsk       ! at the relevant point (T or F)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pUx     ! ice i-velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pVx     ! ice j-velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: p1md       ! (1 - damage) @ T
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj
      REAL(wp) :: zdt, z1_dt, zd
      CHARACTER(len=64) :: cstr
      !
      CHARACTER(len=17), PARAMETER :: crtnm = 'ice_dyn_adv_pra_d'
      REAL(wp), PARAMETER :: rr_scl_fct = 1.E-6_wp
      REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start(crtnm)
      !$acc data present( pe1e2, p1_e1e2, pmsk, pUx, pVx, p1md )

      cstr = 'Prather advection scheme for BBM (damage only)'
      IF( kt == nit000 .AND. lwp ) WRITE(numout,*) '-- '//crtnm//': '//TRIM(cstr)//' at '//cgt//'-points'

      zdt = rDt_ice
      z1_dt = 1._wp / zdt


      !! ADVECT !

      IF( cgt=='T' ) THEN

         CALL adv_pra_2d( kt, zdt, pUx, pVx, pe1e2, p1_e1e2, pmsk,  p1md, sx1md_t, sxx1md_t, sy1md_t, syy1md_t, sxy1md_t ) !--- ice damage @ T

         ! --- Lateral boundary conditions --- !
         !     caution: for gradients (sx and sy) the sign changes
         !# if defined _OPENACC
         !         PRINT *, ' *** Skipping `lbc_lnk` for DAMAGE@T Prather !, kt =',kt
         !# else
# if ! defined _OPENACC
         CALL lbc_lnk( crtnm, p1md,cgt,1._wp,     sx1md_t,cgt,-1._wp, sy1md_t,cgt,-1._wp,  &
            &                 sxx1md_t,cgt,1._wp, syy1md_t,cgt,1._wp, sxy1md_t,cgt,1._wp  )
# endif
         !
      ELSEIF( cgt=='F' ) THEN

         CALL adv_pra_2d( kt, zdt, pUx, pVx, pe1e2, p1_e1e2, pmsk,  p1md, sx1md_f, sxx1md_f, sy1md_f, syy1md_f, sxy1md_f ) !--- ice damage @ F

         !                                                                  !--------------------------------------------!
         ! --- Lateral boundary conditions --- !
         !     caution: for gradients (sx and sy) the sign changes
# if ! defined _OPENACC
         CALL lbc_lnk( crtnm, p1md,cgt,1._wp,    sx1md_f,cgt,-1._wp, sy1md_f,cgt,-1._wp,  &
            &                 sxx1md_f,cgt,1._wp, syy1md_f,cgt,1._wp, sxy1md_f,cgt,1._wp  )
# endif
         !
         !
         !
      ELSE
         CALL ctl_stop('STOP', 'ice_dyn_adv_pra_d: wrong value for grid type `cgt`')
         !
      ENDIF !IF( cgt=='T' )

      IF( lrst_ice )   CALL adv_pra_d_rst( 'WRITE', kt )   !* write Prather fields in the restart file

      !$acc end data
      IF( ln_timing )   CALL timing_stop(crtnm)
      !
   END SUBROUTINE ice_dyn_adv_pra_d



   SUBROUTINE adv_pra_d_init( )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE adv_pra_d_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for Prather advection
      !!-------------------------------------------------------------------
      INTEGER, DIMENSION(3) ::   ierr
      INTEGER :: k_alloc
      CHARACTER(len=16), PARAMETER :: crtnm = 'adv_pra_d_init'
      !!-------------------------------------------------------------------
      ierr(:) = 0
      !                             !* allocate prather fields
      ALLOCATE( sx1md_t(jpi,jpj), sy1md_t(jpi,jpj), sxx1md_t(jpi,jpj), syy1md_t(jpi,jpj), sxy1md_t(jpi,jpj), &
         &      sx1md_f(jpi,jpj), sy1md_f(jpi,jpj), sxx1md_f(jpi,jpj), syy1md_f(jpi,jpj), sxy1md_f(jpi,jpj), &
         &      STAT = ierr(1) )
      sx1md_t(:,:) = 0._wp ;  sy1md_t(:,:) = 0._wp ;  sxx1md_t(:,:) = 0._wp ;  syy1md_t(:,:) = 0._wp ;  sxy1md_t(:,:) = 0._wp
      sx1md_f(:,:) = 0._wp ;  sy1md_f(:,:) = 0._wp ;  sxx1md_f(:,:) = 0._wp ;  syy1md_f(:,:) = 0._wp ;  sxy1md_f(:,:) = 0._wp
      !$acc enter data copyin( sx1md_t, sy1md_t, sxx1md_t, syy1md_t, sxy1md_t )
      !$acc enter data copyin( sx1md_f, sy1md_f, sxx1md_f, syy1md_f, sxy1md_f )
      !
      ALLOCATE( sa2d(jpi,jpj),    STAT = ierr(2) )
      sa2d(:,:) = 0._wp
      !$acc enter data copyin( sa2d )


      ! If Prather is not used to advect generic fields, then we must allocate the following arrays
      !! => because `adv_pra_init` is not doing it...
      IF( .NOT. ln_adv_Pra ) THEN
         ALLOCATE( zfld(jpi,jpj), zf0(jpi,jpj), zbet(jpi,jpj), zfm(jpi,jpj), zfx(jpi,jpj), zfy(jpi,jpj),  &
            &      zfxx(jpi,jpj), zfyy(jpi,jpj), zfxy(jpi,jpj), zpm(jpi,jpj), zpx(jpi,jpj), zpy(jpi,jpj), &
            &      zpxx(jpi,jpj), zpyy(jpi,jpj), zpxy(jpi,jpj), zalg(jpi,jpj), zalg1(jpi,jpj), zalg1q(jpi,jpj), STAT = ierr(3) )
         zfld(:,:) = 0._wp; zf0(:,:) = 0._wp; zbet(:,:) = 0._wp; zfm(:,:) = 0._wp; zfx(:,:) = 0._wp; zfy(:,:) = 0._wp
         zfxx(:,:) = 0._wp; zfyy(:,:) = 0._wp; zfxy(:,:) = 0._wp; zpm(:,:) = 0._wp; zpx(:,:) = 0._wp; zpy(:,:) = 0._wp
         zpxx(:,:) = 0._wp; zpyy(:,:) = 0._wp; zpxy(:,:) = 0._wp; zalg(:,:) = 0._wp; zalg1(:,:) = 0._wp; zalg1q(:,:) = 0._wp
# if defined _OPENACC
         PRINT *, ' * info GPU: adv_pra_d_init() => adding Prather advection workspace arrays to memory'
         PRINT *, '            => zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy'
         PRINT *, '            => zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q'
         !$acc enter data copyin( zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )
# endif
      ENDIF


      k_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( crtnm, k_alloc )
      IF( k_alloc > 0 ) CALL ctl_stop('STOP', crtnm//' : unable to allocate 1md array for Prather advection scheme')
      !
      CALL adv_pra_d_rst( 'READ' )    !* read or initialize all required files
      !
   END SUBROUTINE adv_pra_d_init


   SUBROUTINE adv_pra_d_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE adv_pra_d_rst  ***
      !!
      !! ** Purpose :   Read or write file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER ::   iter     ! local integer
      INTEGER ::   id1      ! local integer
      !!----------------------------------------------------------------------
      !                                      !==========================!
      IF( TRIM(cdrw) == 'READ' ) THEN        !==  Read or initialize  ==!
         !                                   !==========================!
         !
         IF( ln_rstart ) THEN
            id1 = iom_varid( numrir, 'sx1md_t', ldstop = .FALSE. )    ! file exist: id1>0
         ELSE
            id1 = 0                                                  ! no restart: id1=0
         ENDIF
         !
         IF( id1 > 0 ) THEN                     !**  Read the restart file  **!
            !
            !                                                        ! ice damage !#bbm
            CALL iom_get( numrir, jpdom_auto, 'sx1md_t' ,  sx1md_t ,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sy1md_t' ,  sy1md_t ,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxx1md_t', sxx1md_t )
            CALL iom_get( numrir, jpdom_auto, 'syy1md_t', syy1md_t )
            CALL iom_get( numrir, jpdom_auto, 'sxy1md_t', sxy1md_t )
            !
            CALL iom_get( numrir, jpdom_auto, 'sx1md_f' ,  sx1md_f ,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sy1md_f' ,  sy1md_f ,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxx1md_f', sxx1md_f )
            CALL iom_get( numrir, jpdom_auto, 'syy1md_f', syy1md_f )
            CALL iom_get( numrir, jpdom_auto, 'sxy1md_f', sxy1md_f )
            !
         ELSE                                   !**  start rheology from rest  **!
            !
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest OR previous run without Prather, set moments to 0'
            !
            sx1md_t = 0._wp   ;   sy1md_t = 0._wp   ;   sxx1md_t = 0._wp   ;   syy1md_t = 0._wp   ;   sxy1md_t = 0._wp      ! ice damage
            sx1md_f = 0._wp   ;   sy1md_f = 0._wp   ;   sxx1md_f = 0._wp   ;   syy1md_f = 0._wp   ;   sxy1md_f = 0._wp      ! ice damage
            !
            ! ==> no need to send to GPU here, it's done via `$acc enter data` right after the call to present routine for 'READ' mode...
            !
         ENDIF !IF( id1 > 0 )
         !                                   !=====================================!
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   !==  write in the ice restart file  ==!
         !                                   !=====================================!
         IF(lwp) WRITE(numout,*) '----  ice-adv-rst-d  ----'
         iter = kt                           ! ice restarts are written at kt == nitrst
         !
         !$acc update self( sx1md_t, sy1md_t, sxx1md_t, syy1md_t, sxy1md_t )
         !$acc update self( sx1md_f, sy1md_f, sxx1md_f, syy1md_f, sxy1md_f )
         !
         ! In case Prather scheme is used for advection, write second order moments
         ! ------------------------------------------------------------------------
         CALL iom_rstput( iter, nitrst, numriw, 'sx1md_t' , sx1md_t  )
         CALL iom_rstput( iter, nitrst, numriw, 'sy1md_t' , sy1md_t  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxx1md_t', sxx1md_t )
         CALL iom_rstput( iter, nitrst, numriw, 'syy1md_t', syy1md_t )
         CALL iom_rstput( iter, nitrst, numriw, 'sxy1md_t', sxy1md_t )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sx1md_f' , sx1md_f  )
         CALL iom_rstput( iter, nitrst, numriw, 'sy1md_f' , sy1md_f  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxx1md_f', sxx1md_f )
         CALL iom_rstput( iter, nitrst, numriw, 'syy1md_f', syy1md_f )
         CALL iom_rstput( iter, nitrst, numriw, 'sxy1md_f', sxy1md_f )
         !
      ENDIF
      !
   END SUBROUTINE adv_pra_d_rst


   !!======================================================================
END MODULE icedyn_adv_pra_d
