MODULE icedyn
   !!======================================================================
   !!                     ***  MODULE  icedyn  ***
   !!   Sea-Ice dynamics : master routine for sea ice dynamics
   !!======================================================================
   !! history :  4.0  ! 2018  (C. Rousset)  original code SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn       : dynamics of sea ice
   !!   ice_dyn_init  : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE par_ice
   USE ice            ! sea-ice: variables
   USE icedyn_rhg     ! sea-ice: rheology
   USE icedyn_adv     ! sea-ice: advection
   USE icedyn_rdgrft  ! sea-ice: ridging/rafting
   USE icecor         ! sea-ice: corrections
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing
   USE fldread        ! read input fields

   USE oss_nnq  , ONLY: ln_cpl_oce, ssu_m, ssv_m, ssu_v_m, ssv_u_m
   USE ossprs   , ONLY: ln_ssv_Fgrid

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn        ! called by icestp.F90
   PUBLIC   ice_dyn_init   ! called by icestp.F90

   INTEGER ::              nice_dyn   ! choice of the type of dynamics
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_dynALL     = 1   ! full ice dynamics               (rheology + advection + ridging/rafting + correction)
   INTEGER, PARAMETER ::   np_dynRHGADV  = 2   ! pure dynamics                   (rheology + advection)
   INTEGER, PARAMETER ::   np_dynADV1D   = 3   ! only advection 1D - test case from Schar & Smolarkiewicz 1996
   INTEGER, PARAMETER ::   np_dynADV2D   = 4   ! only advection 2D w prescribed vel.
   !                                           !   => uses ocean velocities provided in netCDF SSX forcing as sea-ice velocities

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_icbmsk   ! structure of input grounded icebergs mask (file informations, fields read)

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icedyn.F90 14997 2021-06-16 06:43:57Z smasson $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn  ***
      !!
      !! ** Purpose :   this routine manages sea ice dynamics
      !!
      !! ** Action : - calculation of friction in case of landfast ice
      !!             - call ice_dyn_rhg    = rheology
      !!             - call ice_dyn_adv    = advection
      !!             - call ice_dyn_rdgrft = ridging/rafting
      !!             - call ice_cor        = corrections if fields are out of bounds
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      !!
      INTEGER  ::   ji, jj, jl        ! dummy loop indices
      LOGICAL  ::   lAp
      !!--------------------------------------------------------------------
      !
      ! controls
      IF( ln_timing )   CALL timing_start('ice_dyn')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn: sea-ice dynamics'
         WRITE(numout,*)'~~~~~~~'
      ENDIF
      !
      ! retrieve thickness from volume for landfast param. and UMx advection scheme

      !$acc parallel loop collapse(3)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl=1, jpl
               !
               lAp = ( a_i(ji,jj,jl) >= epsi20 )
               h_i(ji,jj,jl) = MERGE( v_i(ji,jj,jl) / MAX( a_i(ji,jj,jl), epsi20 )  ,  0._wp  , lAp )
               h_s(ji,jj,jl) = MERGE( v_s(ji,jj,jl) / MAX( a_i(ji,jj,jl), epsi20 )  ,  0._wp  , lAp )

               !
               !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               !   IF( a_ip(ji,jj,jl) >= epsi20 ) THEN
               !      h_ip(ji,jj,jl) = v_ip(ji,jj,jl) / a_ip(ji,jj,jl)
               !      h_il(ji,jj,jl) = v_il(ji,jj,jl) / a_ip(ji,jj,jl)
               !   ELSE
               !      h_ip(ji,jj,jl) = 0._wp
               !      h_il(ji,jj,jl) = 0._wp
               !   ENDIF
               !ENDIF
               !
            END DO
         END DO
      END DO
      !$acc end parallel loop

      IF( ln_landfast_L16 ) THEN
         CALL fld_read( kt, 1, sf_icbmsk )
         icb_mask(:,:) = sf_icbmsk(1)%fnow(:,:,1)
      ENDIF

      SELECT CASE( nice_dyn )          !-- Set which dynamics is running

      CASE ( np_dynALL )           !==  all dynamical processes  ==!
         !
         CALL ice_dyn_rhg   ( kt )                                     ! -- rheology
         CALL ice_dyn_adv   ( kt )                                     ! -- advection of ice
         CALL ice_dyn_rdgrft( kt )                                     ! -- ridging/rafting
         CALL ice_cor       ( kt , 1 )                                 ! -- Corrections

      CASE ( np_dynRHGADV  )       !==  no ridge/raft & no corrections ==!
         !
         CALL ice_dyn_rhg( kt )                                     ! -- rheology
         !PRINT *, ' LOLOdyn1' ; CALL test4inf( ' a_i@icedyn 3 ice_dyn  ', a_i )

         CALL ice_dyn_adv( kt )                                     ! -- advection of ice
         !PRINT *, ' LOLOdyn2' ;  CALL test4inf( ' a_i@icedyn 4 ice_dyn  ', a_i )

         CALL ice_var_hpiling()                                     ! -- simple pile-up (replaces ridging/rafting)

         IF( ln_icethd ) THEN
            CALL ice_var_zapsmall()                                    ! -- zap small areas
            !    --> a_ip,h_ip,v_ip,h_il,v_il
         ELSE
            CALL ice_var_zapsmall_dyn()
         ENDIF

         !
         !CALL ice_var_cap_at()  !LOLO: not needed when `hpiling` & `zapsmall` on    ! -- correct `a_i` so that `at_i` remains below `rn_amax`...




      CASE ( np_dynADV1D )         !==  pure advection ==!   (1D)
         !
         CALL ctl_stop( 'ice_dyn: FIXME! finish adapting `dynADV1D` as done in `dynADV2D` !' )
         ! --- monotonicity test from Schar & Smolarkiewicz 1996 --- !
         ! CFL = 0.5 at a distance from the bound of 1/6 of the basin length
         ! Then for dx = 2m and dt = 1s => rn_uice = u (1/6th) = 1m/s
         !DO jj=Njs0-1, Nje0+1
         !   DO ji=Nis0-1, Nie0+1
         !      zcoefu = ( REAL(jpiglo+1)*0.5_wp - REAL(ji+nimpp-1) ) / ( REAL(jpiglo+1)*0.5_wp - 1._wp )
         !      zcoefv = ( REAL(jpjglo+1)*0.5_wp - REAL(jj+njmpp-1) ) / ( REAL(jpjglo+1)*0.5_wp - 1._wp )
         !      u_ice(ji,jj) = rn_uice * 1.5_wp * SIGN( 1.0_wp, zcoefu ) * ABS( zcoefu ) * umask(ji,jj,1)
         !      v_ice(ji,jj) = rn_vice * 1.5_wp * SIGN( 1.0_wp, zcoefv ) * ABS( zcoefv ) * vmask(ji,jj,1)
         !   END DO
         !END DO
         ! ---
         !CALL ice_dyn_adv   ( kt )                                          ! -- advection of ice


      CASE ( np_dynADV2D )         !==  pure advection ==!   (2D w prescribed velocities)
         !
         !  => uses ocean velocities provided in netCDF SSX forcing as sea-ice velocities:
         !
         !$acc parallel loop collapse(2) present( u_ice, v_ice, uVice, vUice, ssu_m, ssv_m )
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               u_ice(ji,jj) = ssu_m(ji,jj)
               v_ice(ji,jj) = ssv_m(ji,jj)
               IF( ln_ssv_Fgrid ) THEN
                  uVice(ji,jj) = ssu_v_m(ji,jj)
                  vUice(ji,jj) = ssv_u_m(ji,jj)
               ELSE
                  uVice(ji,jj) = 0._wp
                  vUice(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
         !$acc end parallel loop

         CALL ice_dyn_adv( kt )                                     ! -- advection of ice

         IF( .NOT. ln_pureADV2D ) THEN
            CALL ice_var_hpiling()                                  ! -- simple pile-up (replaces ridging/rafting)
            CALL ice_var_zapsmall_dyn()
         ENDIF

      END SELECT

# if ! defined _OPENACC
      ! --- Lateral boundary conditions --- !
      !     caution: t_su update needed from itd_reb
      !              plus, one needs ldfull=T to deal with the NorthFold in case of Prather advection
      IF( ln_icethd ) THEN
         !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         !   CALL lbc_lnk( 'icedyn', a_i , 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp, oa_i, 'T', 1._wp, &
         !      &                    t_su, 'T', 1._wp, a_ip, 'T', 1._wp, v_ip, 'T', 1._wp, v_il, 'T', 1._wp, ldfull = .TRUE. )
         !ELSE
         CALL lbc_lnk( 'icedyn', a_i ,'T',1._wp, v_i ,'T',1._wp, v_s ,'T',1._wp, sv_i,'T',1._wp, oa_i,'T',1._wp, &
            &                    t_su,'T',1._wp )
         !ENDIF
         CALL lbc_lnk( 'icedyn', e_i,'T',1._wp, e_s,'T',1._wp, szv_i,'T',1._wp )
      ELSE
         CALL lbc_lnk( 'icedyn', a_i ,'T',1._wp, v_i ,'T',1._wp, v_s ,'T',1._wp )
      ENDIF
# endif

      IF( ln_timing )   CALL timing_stop ('ice_dyn')
      !
   END SUBROUTINE ice_dyn





   SUBROUTINE ice_dyn_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio, ierror   ! Local integer output status for namelist read
      INTEGER ::   ji, jj
      !
      CHARACTER(len=256) ::   cn_dir     ! Root directory for location of ice files
      TYPE(FLD_N)        ::   sn_icbmsk  ! informations about the grounded icebergs field to be read
      !!
      NAMELIST/namdyn/ ln_dynALL, ln_dynRHGADV, ln_dynADV1D, ln_dynADV2D, ln_pureADV2D, rn_ishlat, &
         &             ln_landfast_L16, rn_lf_depfra, rn_lf_bfr, rn_lf_relax, rn_lf_tensile,       &
         &             sn_icbmsk, cn_dir
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn in reference namelist' )
      READ_NML_CFG(numnam_ice,namdyn)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn in configuration namelist' )
      IF(lwm) WRITE( numoni, namdyn )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn:'
         WRITE(numout,*) '      Full ice dynamics      (rhg + adv + ridge/raft + corr) ln_dynALL       = ', ln_dynALL
         WRITE(numout,*) '      No ridge/raft & No cor (rhg + adv)                     ln_dynRHGADV    = ', ln_dynRHGADV
         WRITE(numout,*) '      Advection 1D only      (Schar & Smolarkiewicz 1996)    ln_dynADV1D     = ', ln_dynADV1D
         WRITE(numout,*) '      Advection 2D only (with oce current used as ice vel)   ln_dynADV2D     = ', ln_dynADV2D
         WRITE(numout,*) '               ==> skip ALL post-advection corrections ?        ln_pureADV2D = ', ln_pureADV2D
         WRITE(numout,*) '      lateral boundary condition for sea ice dynamics        rn_ishlat       = ', rn_ishlat
         WRITE(numout,*) '      Landfast: param from Lemieux 2016                      ln_landfast_L16 = ', ln_landfast_L16
         WRITE(numout,*) '         fraction of ocean depth that ice must reach         rn_lf_depfra    = ', rn_lf_depfra
         WRITE(numout,*) '         maximum bottom stress per unit area of contact      rn_lf_bfr       = ', rn_lf_bfr
         WRITE(numout,*) '         relax time scale (s-1) to reach static friction     rn_lf_relax     = ', rn_lf_relax
         WRITE(numout,*) '         isotropic tensile strength                          rn_lf_tensile   = ', rn_lf_tensile
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_dynALL .AND. (.NOT. ln_icethd) ) THEN
         ln_dynALL    = .FALSE.
         ln_dynRHGADV = .TRUE.
         CALL ctl_warn( 'ice_dyn_init: forcing `ln_dynALL=F` & `ln_dynRHGADV=T` because thermo is OFF!' )
      ENDIF
      !
      !                             !== set the choice of ice dynamics ==!
      ioptio = 0
      !      !--- full dynamics                               (rheology + advection + ridging/rafting + correction)
      IF( ln_dynALL    ) THEN
         ioptio = ioptio + 1
         nice_dyn = np_dynALL
      ENDIF
      !      !--- dynamics without ridging/rafting and corr   (rheology + advection)
      IF( ln_dynRHGADV ) THEN
         ioptio = ioptio + 1
         nice_dyn = np_dynRHGADV
      ENDIF
      !      !--- advection 1D only - test case from Schar & Smolarkiewicz 1996
      IF( ln_dynADV1D  ) THEN
         ioptio = ioptio + 1
         nice_dyn = np_dynADV1D
      ENDIF
      !      !--- advection 2D only with prescribed ice velocities (from namelist)
      IF( ln_dynADV2D  ) THEN
         IF( jpl /= 1   ) CALL ctl_stop( 'ice_dyn_init: only 1 ice category for 2D advection test-case please !' )
         IF( ln_icethd  ) CALL ctl_stop( 'ice_dyn_init: `ln_icethd` cannot be used with `ln_dynADV2D=T` ! ' )
         IF( ln_cpl_oce ) CALL ctl_stop( 'ice_dyn_init: `ln_dynADV2D` cannot be set to true in coupled mode! ' )
         ioptio = ioptio + 1
         nice_dyn = np_dynADV2D
      ELSE
         ln_pureADV2D = .FALSE.  ! force to false as it makes no sense otherwize
      ENDIF
      !
      IF( ioptio /= 1 )    CALL ctl_stop( 'ice_dyn_init: one and only one ice dynamics option has to be defined ' )


      !                                      !--- Lateral boundary conditions
      IF    (      rn_ishlat == 0.                ) THEN
         IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  free-slip'
      ELSEIF(      rn_ishlat == 2.                ) THEN
         IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  no-slip'
      ELSEIF( 0. < rn_ishlat .AND. rn_ishlat < 2. ) THEN
         IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  partial-slip'
      ELSEIF( 2. < rn_ishlat                      ) THEN
         IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  strong-slip'
      ENDIF
      !
      !!      The vorticity mask (fmask) is deduced from tmask taking
      !!      into account the choice of lateral boundary condition (rn_ishlat) :
      !!         rn_ishlat = 0, free slip  (no shear along the coast)
      !!         rn_ishlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_ishlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_ishlat, strong slip        | in the lateral boundary layer
      !
      fmask(:,:,:) = 0._wp
      !
      IF( rn_ishlat == 0._wp ) THEN
         IF(lwp )PRINT *, ' * IMPORTANT: `icedyn.F90` building sea-ice `fmask` (shear) with rn_ishlat=0 !'
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               fmask(ji,jj,1) = xmskt(ji,jj) * xmskt(ji+1,jj) * xmskt(ji,jj+1) * xmskt(ji+1,jj+1)
            END DO
         END DO
         !
      ELSE
         IF(lwp )PRINT *, ' * IMPORTANT: `icedyn.F90` building sea-ice `fmask` (shear) with rn_ishlat/=0 !'
         DO jj = Njs0, Nje0
            DO ji = Nis0, Nie0
               fmask(ji,jj,1) = xmskt(ji,jj) * xmskt(ji+1,jj) * xmskt(ji,jj+1) * xmskt(ji+1,jj+1)
               ! Lateral boundary conditions on velocity (modify fmask)
               IF( fmask(ji,jj,1) == 0._wp ) THEN
                  fmask(ji,jj,1) = rn_ishlat * MIN( 1._wp , MAX( umask(ji,jj,1), umask(ji,jj+1,1), &
                     &                                          vmask(ji,jj,1), vmask(ji+1,jj,1) ) )
               ENDIF
            END DO
         END DO
         !
      ENDIF
      CALL lbc_lnk( 'icedyn_rhg_evp', fmask,'F',1._wp )
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_dyn_init() => adding `fmask` arrays to memory!'
      !$acc enter data copyin( fmask )
# endif

      !                                      !--- Landfast ice
      IF( .NOT.ln_landfast_L16 )   tau_icebfr(:,:) = 0._wp
      !
      !                                      !--- allocate and fill structure for grounded icebergs mask
      IF( ln_landfast_L16 ) THEN
         ALLOCATE( sf_icbmsk(1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'ice_dyn_init: unable to allocate sf_icbmsk structure' ) ; RETURN
         ENDIF
         !
         CALL fld_fill( sf_icbmsk, (/ sn_icbmsk /), cn_dir, 'ice_dyn_init',   &
            &                                               'landfast ice is a function of read grounded icebergs', 'icedyn' )
         !
         ALLOCATE( sf_icbmsk(1)%fnow(jpi,jpj,1) )
         IF( sf_icbmsk(1)%ln_tint )   ALLOCATE( sf_icbmsk(1)%fdta(jpi,jpj,1,2) )
         IF( TRIM(sf_icbmsk(1)%clrootname) == 'NOT USED' ) sf_icbmsk(1)%fnow(:,:,1) = 0._wp   ! not used field  (set to 0)
      ELSE
         icb_mask(:,:) = 0._wp
      ENDIF


# if defined _OPENACC
      PRINT *, ' * info GPU: ice_dyn_init() => adding ice velocities arrays to memory!'
      !$acc enter data copyin( u_ice, v_ice, uVice, vUice )
      PRINT *, '    ==> u_ice, v_ice, uVice, vUice'
# endif

      IF( .NOT. ln_dynADV2D ) THEN

         CALL ice_dyn_rhg_init             ! set ice rheology parameters

         !                                      !--- other init
         IF(.NOT. ln_rhg_EVP) rn_creepl = 2.E-9_wp
         rn_delta_ecc = 2._wp ! default in case `ice_dyn_rdgrft_init` not called (`rn_delta_ecc` always needed by `strain_rate_dsd`!)
         IF( ln_dynALL .OR. ln_rhg_EVP ) THEN   ! `ln_rhg_EVP` => because we need params initialized by `ice_dyn_rdgrft_init` for `ice_strength` in EVP !
            CALL ice_dyn_rdgrft_init          ! set ice ridging/rafting parameters
         ENDIF

         IF( ln_rhg_EVP ) THEN
            rn_delta_ecc = rn_ecc
            IF(lwp) WRITE(numout,*) '  *** since EVP used forcing `rn_delta_ecc` to `rn_ecc`! => rn_delta_ecc=',rn_delta_ecc
         ENDIF

         !$acc update device ( rn_creepl, rn_delta_ecc )

      ENDIF


      CALL ice_dyn_adv_init             ! set ice advection parameters


   END SUBROUTINE ice_dyn_init

   !!======================================================================
END MODULE icedyn
