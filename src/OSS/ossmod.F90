MODULE ossmod
   !!======================================================================
   !!                       ***  MODULE  ossmod  ***
   !!                      Ocean Surface State module
   !!
   !!======================================================================
   !! History :  0.1  ! 2024-09  (L. Brodeau) starting from `sbcmod.F90` of NEMO 4.2.2
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   oss_init      : read namoss namelist
   !!   oss           :
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbcblk,  ONLY : ln_skin_cs, ln_skin_wl, jp_tair, jp_humi, jp_wndi, jp_wndj, jp_dqsw, jp_dqlw, jp_prcp, jp_snow, sf, ll_skin
   USE sbc_oce
   !
   USE oss_nnq        ! Surface boundary condition: ocean fields
   USE ossskin, ONLY : oss_skin_alloc
   USE ossprs         ! surface boundary condition: sea-surface mean variables
   USE par_ice
   USE ice
   !
   USE remap_classic, ONLY : do_rmpVecT2UV, do_Voce
   !
   USE osscpl         ! surface boundary condition: coupled formulation
   USE cpl_oasis3     ! OASIS routines for coupling
   USE bdy   , ONLY: ln_bdy
   !
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu
#else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#endif
   !
   USE eosbn2, ONLY: eos10_fzp_2d, eos10_fzp_2d_gpu
   !
   USE prtctl         ! Print control                    (prt_ctl routine)
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss             ! routine called by step.F90
   PUBLIC   oss_init        ! routine called by opa.F90
   PUBLIC   oss_blk_ssx     ! routine called by step.F90
   PUBLIC   oss_flx_write   ! deals with all the surface fluxes for the liquid ocean

   INTEGER ::   noss   ! type of surface boundary condition (deduced from namsbc informations)

   !! * Substitutions
#  include "single_precision_substitute.h90"
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! $Id: ossmod.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE oss_init()
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE oss_init ***
      !!
      !! ** Purpose :   Initialisation of the ice surface boundary computation
      !!
      !! ** Method  :   Read the namoss namelist and set derived parameters
      !!                Call init routines for all other OSS modules that have one
      !!
      !! ** Action  : - read namoss parameters
      !!              - noss: type of oss
      !!----------------------------------------------------------------------
      INTEGER ::   ios, icpt                         ! local integer
      !!
      NAMELIST/namoss/ rn_Cd_io, ln_drgice_imp, nn_foss, ln_prs_oce, ln_cpl_oce, ln_ice_embd
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'oss_init : surface boundary condition setting'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      !                       !**  read Surface Module namelist
      READ_NML_REF(numnam,namoss)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namoss in reference namelist' )
      READ_NML_CFG(numnam,namoss)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namoss in configuration namelist' )
      IF(lwm) WRITE( numond, namoss )
      !
      IF(lwp) THEN                  !* Control print
         WRITE(numout,*) '   Namelist namoss (partly overwritten with CPP key setting)'
         WRITE(numout,*) '      drag coefficient for oceanic stress           rn_Cd_io  = ', rn_Cd_io
         WRITE(numout,*) '      implicit ice-ocean drag                ln_drgice_imp  =', ln_drgice_imp
         WRITE(numout,*) '      frequency update of oss (and ice)             nn_foss = ', nn_foss
         WRITE(numout,*) '      Type of coupling (Ocean/Ice/Atmosphere) : '
         WRITE(numout,*) '         prescribed surface ocean state          ln_prs_oce = ', ln_prs_oce
         WRITE(numout,*) '         ice-ocean coupled formulation           ln_cpl_oce = ', ln_cpl_oce
         WRITE(numout,*) '         OASIS coupling                        lk_oasis_oce = ', lk_oasis_oce
         WRITE(numout,*) '         ice embedded into ocean              ln_ice_embd   = ', ln_ice_embd
         !WRITE(numout,*) '      Misc. options of oss : '
         !WRITE(numout,*) '         nb of iterations if land-sea-mask applied  nn_lsm        = ', nn_lsm
      ENDIF
      !
      !IF( MOD( rday , rn_Dt ) /= 0. )   CALL ctl_stop( 'the time step must devide the number of second of in a day' )
      !IF( MOD( rday , 2.  ) /= 0. )   CALL ctl_stop( 'the number of second of in a day must be an even number'    )
      !IF( MOD( rn_Dt  , 2.  ) /= 0. )   CALL ctl_stop( 'the time step (in second) must be an even number'           )
      !
      !                       !**  check option consistency


      !! Coupled to an ocean model or standalone ?
      IF(       ln_prs_oce .AND.      ln_cpl_oce ) CALL ctl_stop( 'oss_init : cannot have `ln_prs_oce=T` with `ln_cpl_oce=T`' )
      IF((.NOT. ln_prs_oce).AND.(.NOT.ln_cpl_oce)) CALL ctl_stop( 'oss_init : pick `ln_prs_oce=T` OR `ln_cpl_oce=T`' )
      IF( ln_cpl_oce ) THEN
         IF(lwp) WRITE(numout,*) '   ==>>>  NANUQ is coupled to an ocean model!'
         IF( .NOT.lk_oasis_oce )   CALL ctl_stop( 'oss_init : coupled mode but key_oasis3 disabled' )
      ELSE
         IF(lwp) WRITE(numout,*) '   ==>>>  NANUQ will use a prescribed surface state of the ocean!'
      ENDIF

      !
      !                       !**  allocate and set required variables

      !                             !* allocate oss arrays
      IF( oss_nnq_alloc() /= 0 )   CALL ctl_stop( 'oss_init : unable to allocate `oss_nnq` arrays' )

      !                             !* allocate skin arrays
      IF(lwp) PRINT *, ' *** Calling `oss_skin_alloc` with ln_skin_cs, ln_skin_wl =', ln_skin_cs, ln_skin_wl
      IF( oss_skin_alloc( l_use_cs=ln_skin_cs, l_use_wl=ln_skin_wl ) /= 0 )  CALL ctl_stop( 'oss_init : unable to allocate `oss_skin` arrays' )
      !
      !
      !                             !* Choice of the Surface Boudary Condition
      !                             (set noss)
      !
      noss = 1                  ! prescribed sea-surface state (default)
      IF( ln_cpl_oce ) noss = 2 ! coupled to ocean model
      !
      IF(lwp) THEN                     !- print the choice of surface flux formulation
         WRITE(numout,*) '  *** As a surface boundary condition NANUQ will use:'
         WRITE(numout,*)
         SELECT CASE( noss )
         CASE( jp_prs_oce ) ;   WRITE(numout,*) '   ==>>>   prescribed sea-surface state'
         CASE( jp_cpl_oce ) ;   WRITE(numout,*) '   ==>>>   coupled to ocean model'
         END SELECT
      ENDIF
      !
      !                             !* OASIS initialization
      !
      !!IF( ln_cpl )   CALL sbc_cpl_init( nn_ice )   ! Must be done before: (1) first time step
      !!                                             !                      (2) the use of nn_fsbc
      IF( lk_oasis_oce )  THEN
         CALL oss_cpl_init( 2 )    ! Must be done before: (1) first time step
         !                         !                     (2) the use of nn_foss
         !
         CALL cpl_enddef           ! terminate coupling initialization LOLO: must be here !!! Before call to `cpl_freq` !!!
      ENDIF




      !     nn_foss initialization if OCE-SAS coupling via OASIS
      !     SAS time-step has to be declared in OASIS (mandatory) -> nn_foss has to be modified accordingly
      IF( ln_cpl_oce ) THEN
#if defined key_verbose
         IF(lwp) PRINT *, ' * LOLO calling `cpl_freq` with midcpl=',midcpl
#endif
         ! try: I_OTaux1
         nn_foss = cpl_freq('I_SFLX',midcpl) / NINT(rn_Dt)
#if defined key_verbose
         IF(lwp) PRINT *, ' * LOLO DONE with `cpl_freq` ! nn_foss =', nn_foss
#endif
         !
         IF(lwp)THEN
            WRITE(numout,*)
            WRITE(numout,*)"   NANUQ to an ocean model via OASIS : nn_foss re-defined from OASIS namcouple ", nn_foss
            WRITE(numout,*)
         ENDIF
      ENDIF
      !
      !                             !* check consistency between model timeline and nn_foss
      IF( ln_rst_list .OR. nn_stock /= -1 ) THEN   ! we will do restart files
         IF( MOD( nitend - nit000 + 1, nn_foss) /= 0 ) THEN
            WRITE(ctmp1,*) 'oss_init : experiment length (', nitend - nit000 + 1, ') is NOT a multiple of nn_foss (', nn_foss, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
         IF( .NOT. ln_rst_list .AND. MOD( nn_stock, nn_foss) /= 0 ) THEN   ! we don't use nn_stock if ln_rst_list
            WRITE(ctmp1,*) 'oss_init : nn_stock (', nn_stock, ') is NOT a multiple of nn_foss (', nn_foss, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
      ENDIF
      !
      IF( MOD( rday, REAL(nn_foss, wp) * rn_Dt ) /= 0 )   &
         &  CALL ctl_warn( 'oss_init : nn_foss is NOT a multiple of the number of time steps in a day' )


      CALL oss_prs_init() ! Prescribed sea surface state fields initialization
      !                   ! or initial state in coupled mode

      !$acc update device ( rn_Cd_io, ln_drgice_imp, nn_foss, ln_cpl_oce, ln_ice_embd )

      IF( .NOT. ln_rstart ) THEN
         IF(lwp) WRITE(numout,*)"  'oss_init()' => `ssst`, `sssq`, `sst_s` & `sss_s` initialized with `sst_m`, `sst_m` & `sss_m`!"
         sst_s(:,:) = sst_m(:,:) ! slab bulk SST
         sss_s(:,:) = sss_m(:,:) ! slab bulk SSS
         !
         ssst(:,:)  = sst_m(:,:) ! water skin temperature
         !sssq(:,:)  = rdct_qsat_salt * q_sat( ssst(:,:)+rt0, sst_m(:,:)*0._wp+101000. )   ! (kg/kg)
         sssq(:,:) = 0._wp
         !$acc update device ( ssst, sssq, sst_s, sss_s )
      ELSE
         !
         IF( ln_cpl_oce .OR. (ln_prs_oce .AND. (.NOT. ln_slab_sst) ) ) THEN
            sst_s(:,:) = sst_m(:,:) ! slab bulk SST
            sss_s(:,:) = sss_m(:,:) ! slab bulk SSS
            ssst(:,:)  = sst_m(:,:) ! water skin temperature
            sssq(:,:) = 0._wp
            !$acc update device ( ssst, sssq, sst_s, sss_s )
         ELSE
            CALL ctl_stop( 'oss_init : FIXME!!! => add restart capability for `sst_s` & `sss_s` ! (-->ossmod.F90)' )
         ENDIF
         !
         IF( ll_skin ) CALL ctl_stop( 'oss_init : FIXME!!! => add restart capability for `ssst` (skin) ! (-->ossmod.F90)' )
         !
      ENDIF

   END SUBROUTINE oss_init


   SUBROUTINE oss( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE oss  ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!                condition (momentum, heat and freshwater fluxes)
      !!
      !! ** Method  :   blah blah  to be written ?????????
      !!                CAUTION : never mask the surface stress field (tke oss)
      !!
      !! ** Action  : - set the ocean surface boundary condition at before and now
      !!                time step, i.e.
      !!
      !! zus -> xtmp1
      !! zvs -> xtmp2
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      INTEGER  ::   jj, ji
      REAL(wp) ::   zthscl        ! wd  tanh scale
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('oss')
      !
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      !
      !                                            !==  oss formulation  ==!
      IF(ln_cpl_oce) THEN
         CALL oss_cpl_rcv( kt, nn_foss, 2 ) ! => coupled to an ocean model: NANUQ receives sea surface state fields via OASIS
      ELSE
         CALL oss_prs_rcv( kt )             ! => prescribed sea surface state read into netCDF files
      END IF
      !! Either way, we have updated the following fields:
      !!  sst_m, sss_m, ssh_m, ssu_m, ssv_m, ( mld_m, frq_m, e3t_m )



      !#lolo: I think that is the correct place to output the SSX field for the ocean we just read or received:
      !!  No update from GPU to CPU required because they have just been read by the CPU...
      IF( iom_use('ssh_m') ) CALL iom_put( 'ssh_m', ssh_m )
      IF( iom_use('e3t_m') ) CALL iom_put( 'e3t_m', e3t_m )
      IF( iom_use('frq_m') ) CALL iom_put( 'frq_m', frq_m )
      IF( iom_use('frz_m') ) THEN
         CALL eos10_fzp_2d( sss_m(:,:), xtmp1(:,:) )
         CALL iom_put( 'frz_m', xtmp1 * xmskt )
         xtmp1(:,:) = 0._wp
      ENDIF

#if defined _OPENACC || defined _OPENMP
      ! ==> Sending all the `*_m` arrays to GPU memory
      !$acc update device ( ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m )
      IF( ln_ssv_Fgrid ) THEN
         !PRINT *, ' *** oss@ossmod.F90 => updating `ssu_v_m` & `ssv_u_m` on the GPU!'
         !$acc update device ( ssu_v_m, ssv_u_m )
      ENDIF
      IF( ln_slab_sst ) THEN
         IF( iom_use('mld_m') ) CALL iom_put( 'mld_m', mld_m )
         !$acc update device ( mld_m )
      ENDIF
#else
      IF( ln_slab_sst ) THEN
         IF( iom_use('mld_m') ) CALL iom_put( 'mld_m', mld_m )
      ENDIF
#endif

      ! If prescribed velocities provided at T- rather than U,V- points, interpolate from T to U & T to V
      IF( ln_ssv_T ) THEN
#if defined key_verbose
         IF(lwp) PRINT *, ' * [oss@ossmod.F90] => interpolates `SSU,SSV` from T to U,V!', kt
#endif
         !$acc data present( xtmp1, xtmp2 )
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               xtmp1(ji,jj) = ssu_m(ji,jj)
               xtmp2(ji,jj) = ssv_m(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
         !$acc end data
         CALL do_rmpVecT2UV( xtmp1, xtmp2, ssu_m, ssv_m )
         !
#if defined _OPENACC || defined _OPENMP
         IF( iom_use('ssu_m') ) THEN
            !$acc update self ( ssu_m )
         ENDIF
         IF( iom_use('ssv_m') ) THEN
            !$acc update self ( ssv_m )
         ENDIF
#endif
      ENDIF !IF( ln_ssv_T )

      IF( iom_use('ssu_m') ) CALL iom_put( 'ssu_m', ssu_m )
      IF( iom_use('ssv_m') ) CALL iom_put( 'ssv_m', ssv_m )

      IF( ln_ssv_Fgrid ) THEN
         IF( iom_use('ssu_v_m') ) CALL iom_put( 'ssu_v_m', ssu_v_m )
         IF( iom_use('ssv_u_m') ) CALL iom_put( 'ssv_u_m', ssv_u_m )
         !CALL ctl_stop( 'oss: WORRIED that `ssu_v_m` & `ssv_u_m` not sent to V_oce(:,:,3) & V_oce(:,:,4) ?' )
      ENDIF


      ! -- mean surface ocean current for the E-grid
      !    => 4 compoments (x@U, y@V, x@V, y@U) stored into array `V_oce`

      IF( ln_ssv_Fgrid ) THEN
         !! We have read Uv & Vu in the netCDF file!
         !$acc data present( ssu_m, ssv_m, ssu_v_m, ssv_u_m, V_oce )
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               V_oce(ji,jj,1) = ssu_m(ji,jj)
               V_oce(ji,jj,2) = ssv_m(ji,jj)
               V_oce(ji,jj,3) = ssu_v_m(ji,jj)
               V_oce(ji,jj,4) = ssv_u_m(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
#if defined _OPENACC || defined _OPENMP
         CALL lbc_lnk_gpu( '', V_oce )
#else
         CALL lbc_lnk(     '', V_oce(:,:,1),'U',-1._wp, V_oce(:,:,2),'V',-1._wp, V_oce(:,:,3),'V',-1._wp, V_oce(:,:,4),'U',-1._wp )
#endif
         !$acc end data
         !
      ELSE
         !!
         !! We have to create Uv & Vu (interpolation)
         CALL do_Voce( ssu_m, ssv_m,  V_oce )
         !
      ENDIF

      IF( ln_timing )   CALL timing_stop('oss')

   END SUBROUTINE oss



   !SUBROUTINE oss_blk_ssx( kt, psst_m, psss_m, pmld_m, pqsr_b, pqns_b, pemp_b,   psst_s, psss_s, pt_bo )
   SUBROUTINE oss_blk_ssx( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE oss_blk_ssx  ***
      !!
      !! ** Purpose :   - prevent prescribed SST to be colder than freezing point
      !!                - apply the SLAB ocean correction if required
      !!
      !! ** Method  :
      !!
      !! ** Action  :   - prevent prescribed SST to be colder than freezing point
      !!                - apply the SLAB ocean correction if required
      !!
      !!                       => uses:       sst_m, sss_m, mld_m, qsr_b, qns_b, emp_b
      !!                       => may update: sst_m (freezing point consistency)
      !!                       => updates:    sst_s, sss_s, t_bo
      !!
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in)    ::   kt
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psst_m  ! bulk SST as prescribed (in netCDF) or received from OASIS [deg.C]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   psss_m  ! bulk SSS as prescribed (in netCDF) or received from OASIS [-]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmld_m  ! MLD as prescribed (in netCDF) (irrelevant when coupled)   [m]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pqsr_b, pqns_b  ! heat fluxes  [W/m^2]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pemp_b  ! E-P freshwater flux [kg/m^2/s]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   psst_s  ! adjusted bulk SST to use [deg.C]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   psss_s  ! adjusted bulk SSS to use [-]
      !REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pt_bo   ! freezing temperature based on SSS [deg.C]
      !!----------------------------------------------------------------------
      INTEGER  ::   jj, ji          ! dummy loop argument
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('oss_blk_ssx')
      !$acc data present ( sst_m, sss_m, mld_m, qsr_b, qns_b, emp_b, sst_s, sss_s, t_bo )

      CALL eos10_fzp_2d_gpu( sss_m, t_bo )   ! -- freezing temperature based on salinity [C]

      IF( ln_prs_oce ) THEN
         !! => we use a prescribed surface state of the ocean read into netCDF files,
         !!    sometimes `sst_m` is not consistent with the equation of state and can
         !!    be colder than the freezing point (true at least in GLORYS4...)
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               sst_m(ji,jj) = MAX( sst_m(ji,jj) , t_bo(ji,jj) ) ! prescribed/observed SST, cannot be colder than freezing-point temperature.
            END DO
         END DO
         !$acc end parallel loop

         !! Update `sst_s`,  `sss_s` & `t_bo` based on a simplistic slab-ocean model approach:
         IF( ln_slab_sst ) THEN
            CALL oss_prs_slab( kt, rDt_ice, sst_m, sss_m, mld_m, qsr_b, qns_b, emp_b, sst_s, sss_s, t_bo )
         ENDIF
         !
      ENDIF !IF( ln_prs_oce .AND. ln_icethd )


      !! Now for coupled oce-ice setup or a standalone run without the SLAB, `sst_s, sss_s` default to `sst_m, sss_m`
      !!  ==> this is important because surface bulk fluxes and sea-ice thermo are going to use `sst_s, sss_s` as BULK SS* !!!
      IF( ln_cpl_oce .OR. (ln_prs_oce .AND. (.NOT. ln_slab_sst) ) ) THEN
         !! ==> `sst_s` & `sss_s` default to `sst_m` & `sst_s`:
#if defined key_verbose
         IF(lwp) PRINT *, ' * `sst_s & sss_s` forced to (prescribed) `sst_m & sss_m` ! kt =', kt
#endif
         !$acc parallel loop collapse(2) present( sst_m, sss_m, sst_s, sss_s )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               sst_s(ji,jj) = sst_m(ji,jj)
               sss_s(ji,jj) = sss_m(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      IF( .NOT. ll_skin ) THEN
         !! No cool-skin/warm-layer param. in use, `ssst` defaults to `sst_s`
#if defined key_verbose
         IF(lwp) PRINT *, 'LOLO: oss_blk_ssx@ossmod.F90: ssst <- sst_s, kt=', kt
#endif
         !$acc parallel loop collapse(2) present( sst_m, sss_m, sst_s, sss_s )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ssst(ji,jj) = sst_s(ji,jj)
            END DO
         END DO
         !$acc end parallel loop

      ENDIF



#if defined _TRDBG
      !$acc update self( at_i, t_bo, oa_i, sst_m, sst_s, sss_m, sss_s )
      CALL TRDBG( 'oss_blk_ssx',  't_bo', t_bo )
      CALL TRDBG( 'oss_blk_ssx',  'om_i', SUM( oa_i(:,:,:), dim=3 ) / MAX(at_i(:,:),epsi20) )
      CALL TRDBG( 'oss_blk_ssx', 'sst_m', sst_m )
      CALL TRDBG( 'oss_blk_ssx', 'sst_s', sst_s )
      CALL TRDBG( 'oss_blk_ssx', 'sst_m', sss_m )
      CALL TRDBG( 'oss_blk_ssx', 'sss_s', sss_s )
#endif


      IF( iom_use('sst_m') .OR. iom_use('dsst_s') ) THEN
         !$acc update self(sst_m)
         IF( iom_use('sst_m') ) CALL iom_put ( 'sst_m' ,  sst_m * xmskt )
      ENDIF
      IF( iom_use('sss_m') .OR. iom_use('dsss_s') ) THEN
         !$acc update self(sss_m)
         IF( iom_use('sss_m') ) CALL iom_put ( 'sss_m' ,  sss_m * xmskt )
      ENDIF


      ! Actual bulk SST & SSS:
      IF( iom_use( 'sst_s') .OR. iom_use('dsst_s') .OR. iom_use('dt_skin') ) THEN
         !$acc update self(sst_s)
         IF( iom_use( 'sst_s') ) CALL iom_put (  'sst_s' ,  sst_s * xmskt )
      ENDIF
      IF( iom_use( 'sss_s') .OR. iom_use('dsss_s') ) THEN
         !$acc update self(sss_s)
         IF( iom_use( 'sss_s') ) CALL iom_put (  'sss_s' ,  sss_s * xmskt )
      ENDIF

      ! Deviation of actual SST & SSS from `sst_m` & `sss_m` due to use of slab-ocean:
      IF( iom_use('dsss_s') ) CALL iom_put ( 'dsss_s' , (sss_s - sss_m)*xmskt )
      IF( iom_use('dsst_s') ) CALL iom_put ( 'dsst_s' , (sst_s - sst_m)*xmskt )
      ! Skin temperature used in surface heat flux estimate (if relevant):

      ! Deviation of `Skin SST` from actual bulk SST due to use of cool-skin/warm-layer:
      IF( iom_use('ssst').OR.iom_use('dt_skin') ) THEN
         !$acc update self ( ssst )
         IF( iom_use('ssst') ) CALL iom_put ( 'ssst'    ,  ssst         *xmskt )
      ENDIF

      IF( iom_use('dt_skin') ) CALL iom_put ( 'dt_skin' , (ssst - sst_s)*xmskt )

      !$acc end data
      IF( ln_timing )   CALL timing_stop('oss_blk_ssx')

   END SUBROUTINE oss_blk_ssx




   SUBROUTINE oss_flx_write( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE oss_flx_write  ***
      !!
      !!       Deals with all the surface fluxes for the liquid ocean
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('oss_flx_write')
      !$acc data present( sfx, fmmflx, emp, qns, qsr, taum, wndm, qsr_oce, qns_oce, qemp_oce )

      !                                                ! ---------------------------------------- !
      !                                                !        Outputs and control print         !
      !                                                ! ---------------------------------------- !

      !! LB => when coupled to ocean component via OASIS I rather save the exacts same fields (as saved here)
      !!       right at the place where they are sent to OASIS, so technically in `oss_cpl_snd()@osscpl.F90` !

      IF( iom_use('saltflx' ) ) THEN
         !$acc update self( sfx )
         CALL iom_put( 'saltflx', sfx * xmskt         ) !#LOLOfixme  ! downward salt flux (includes virtual salt flux beneath ice in linear free surface case)
      ENDIF
      IF( iom_use('fmmflx' ) ) THEN
         !$acc update self( fmmflx )
         CALL iom_put( 'fmmflx' , fmmflx * xmskt ) ! Freezing-melting upward water flux for liquid ocean (`fmmflx>0` => LOSS for the liquid ocean)
      ENDIF
      IF( iom_use('emp' ) ) THEN
         !$acc update self( emp )
         CALL iom_put( 'emp'    , emp * xmskt )    ! E-P (net upward freshwater, (`emp>0` => LOSS for the liquid ocean)
      ENDIF

      !! In the following, a positive flux means a gain for the liquid ocean...
      IF( iom_use('qt') ) THEN
         !$acc update self( qns, qsr )
         CALL iom_put( 'qt'  , (qns + qsr) * xmskt    )       ! total heat flux
      ENDIF
      IF( iom_use('qns' ) ) THEN
         !$acc update self( qns )
         CALL iom_put( 'qns'    , qns ) ! * xmskt )       ! non-solar heat flux
      ENDIF
      IF( iom_use('qsr' ) ) THEN
         !$acc update self( qsr )
         CALL iom_put( 'qsr'    , qsr * xmskt )       ! solar heat flux
      ENDIF
      IF( iom_use('taum' ) ) THEN
         !$acc update self( taum )
         CALL iom_put( 'taum'   ,  taum * xmskt )       ! wind stress module
      ENDIF
      !
      IF( sn_loc_vct_tau=='T' ) THEN
         IF( iom_use('utau_t') ) THEN
            !$acc update self( utau )
            CALL iom_put( 'utau_t',  utau * xmskt )
         ENDIF
         IF( iom_use('vtau_t') ) THEN
            !$acc update self( vtau )
            CALL iom_put( 'vtau_t',  vtau * xmskt )
         ENDIF
      ENDIF
      IF( sn_loc_vct_tau=='C' ) THEN
         IF( iom_use('utau_u') ) THEN
            !$acc update self( utau )
            CALL iom_put( 'utau_u',  utau * umask(:,:,1) )
         ENDIF
         IF( iom_use('vtau_v') ) THEN
            !$acc update self( vtau )
            CALL iom_put( 'vtau_v',  vtau * vmask(:,:,1) )
         ENDIF
      ENDIF
      !
      IF( iom_use('wspd' ) ) THEN
         !$acc update self( wndm )
         CALL iom_put( 'wspd',     wndm * xmskt ) ! wind speed  module over free ocean or leads in presence of sea-ice
      ENDIF
      !

      ! Atmospheric forcing input fields for bulk forcing (for debugging purposes)
      ! **************************************************************************
      IF( iom_use('theta_zt') ) CALL iom_put('theta_zt', (sf(jp_tair)%fnow(:,:,1) -rt0) * xmskt ) ! absolute temperature at z=zt
      IF( iom_use('q_zt')     ) THEN
         !$acc update self( sf(jp_humi)%fnow(:,:,1) )
         CALL                        iom_put('q_zt',      sf(jp_humi)%fnow(:,:,1)       * xmskt ) ! specific humidity    at z=zt
      ENDIF

      IF( iom_use("wnd_x") )   CALL iom_put( "wnd_x",     sf(jp_wndi)%fnow(:,:,1)*xmskt(:,:) )
      IF( iom_use("wnd_y") )   CALL iom_put( "wnd_y",     sf(jp_wndj)%fnow(:,:,1)*xmskt(:,:) )

      IF( iom_use("snowpre" ) ) CALL iom_put( "snowpre",  sf(jp_snow)%fnow(:,:,1)*xmskt(:,:) )    ! output solid precipitation [kg/m2/s]
      IF( iom_use("precip" ) )  CALL iom_put( "precip" ,  sf(jp_prcp)%fnow(:,:,1)*xmskt(:,:) )    ! output total precipitation [kg/m2/s]

      IF( iom_use('theta_zu') .OR. iom_use('dtheta_zu_zt')  ) THEN
         !$acc update self(t_air_zu)
         IF( iom_use('theta_zu') )     CALL iom_put('theta_zu', (t_air_zu -rt0) * xmskt) ! potential temperature at z=zu
         IF( iom_use('dtheta_zu_zt') ) THEN
            !! IF(ln_abl), we are comparing to what was read in the forcing field, not what has been updated by the ABL (because `t_air_zu` forced to that!)
            !$acc update self( sf(jp_tair)%fnow(:,:,1) )
            CALL iom_put('dtheta_zu_zt', (t_air_zu(:,:) - sf(jp_tair)%fnow(:,:,1)) * xmskt) ! diff potential temperature between z=zu and z=zt
         ENDIF
      ENDIF
      IF( iom_use('q_zu')     .OR. iom_use('dq_zu_zt')      ) THEN
         !$acc update self(q_air_zu)
         IF( iom_use('q_zu') )     CALL iom_put('q_zu',      q_air_zu       * xmskt) ! specific humidity       '
         IF( iom_use('dq_zu_zt') ) THEN
            !! IF(ln_abl), we are comparing to what was read in the forcing field, not what has been updated by the ABL (because `q_air_zu` forced to that!)
            !$acc update self( sf(jp_humi)%fnow(:,:,1) )
            CALL iom_put('dq_zu_zt',  (q_air_zu(:,:) - sf(jp_humi)%fnow(:,:,1)) * xmskt) ! diff specific humidity between z=zu and z=zt
         ENDIF
      ENDIF

      IF( iom_use("rad_sw" ) ) CALL iom_put( "rad_sw",  sf(jp_dqsw)%fnow(:,:,1)*xmskt(:,:) ) ! Downwelling shortwave radiative flux available at sea level
      IF( iom_use("rad_lw" ) ) CALL iom_put( "rad_lw",  sf(jp_dqlw)%fnow(:,:,1)*xmskt(:,:) ) ! Downwelling  longwave radiative flux available at sea level
      ! **************************************************************************


      !ENDIF !IF( .NOT. ln_cpl_oce )

      !! Arrays that exist because of presence of sea-ice
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF( iom_use('Cd_oce') ) THEN
         !$acc update self( CD_oce )
         CALL iom_put( 'Cd_oce'    , CD_oce * REAL(1-kmsk_ice_t,wp) * xmskt )
      ENDIF
      IF( iom_use('Ch_oce') ) THEN
         !$acc update self( CH_oce )
         CALL iom_put( 'Ch_oce'    , CH_oce * REAL(1-kmsk_ice_t,wp) * xmskt )
      ENDIF
      IF( iom_use('Ce_oce') ) THEN
         !$acc update self( CE_oce )
         CALL iom_put( 'Ce_oce'    , CE_oce * REAL(1-kmsk_ice_t,wp) * xmskt )
      ENDIF
      IF( iom_use('qsr_oce') ) THEN
         !$acc update self( qsr_oce )
         CALL iom_put( 'qsr_oce'    , qsr_oce * REAL(1-kmsk_ice_t,wp) * xmskt ) ! solar flux availabe at open-ocean surface (computed as if no sea-ice was present)
      ENDIF
      IF( iom_use('qns_oce') ) THEN
         !$acc update self( qns_oce )
         CALL iom_put( 'qns_oce'    , qns_oce * REAL(1-kmsk_ice_t,wp) * xmskt ) ! non-solar flux at ocean surface (qlat+qsen+qlw) (computed as if no sea-ice was present)
      ENDIF
      IF( iom_use('qemp_oce') ) THEN
         !$acc update self( qemp_oce )
         CALL iom_put( 'qemp_oce'   , qemp_oce * REAL(1-kmsk_ice_t,wp) * xmskt ) ! Downward Heat Flux from E-P over open-ocean
      ENDIF
      IF( iom_use('qt_oce') ) THEN
         !$acc update self( qsr_oce, qns_oce, qemp_oce )
         IF( iom_use('qt_oce') ) CALL iom_put( 'qt_oce', (qsr_oce + qns_oce + qemp_oce) * REAL(1-kmsk_ice_t,wp) * xmskt )
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN     ! print mean trends (used for debugging)
         CALL prt_ctl(tab2d_1=CASTDP(qns)    , clinfo1=' qns      - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=CASTDP(qsr)    , clinfo1=' qsr      - : ', mask1=tmask )
         CALL prt_ctl(tab3d_1=CASTDP(tmask)  , clinfo1=' tmask    - : ', mask1=tmask, kdim=1 )
         CALL prt_ctl(tab2d_1=CASTDP(utau)   , clinfo1=' utau     - : ', mask1=umask,        &
            &         tab2d_2=CASTDP(vtau)   , clinfo2=' vtau     - : ', mask2=vmask )
      ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('oss_flx_write')
      !
   END SUBROUTINE oss_flx_write



   !!======================================================================
END MODULE ossmod
