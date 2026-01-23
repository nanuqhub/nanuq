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
   USE oss_nnq        ! Surface boundary condition: ocean fields
   USE ossprs         ! surface boundary condition: sea-surface mean variables
   USE par_ice
   USE ice
   USE osscpl         ! surface boundary condition: coupled formulation
   USE cpl_oasis3     ! OASIS routines for coupling
   USE bdy   , ONLY: ln_bdy
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE eosbn2, ONLY: eos_fzp_2d
   !
   USE prtctl         ! Print control                    (prt_ctl routine)
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss        ! routine called by step.F90
   PUBLIC   oss_init   ! routine called by opa.F90

   INTEGER ::   noss   ! type of surface boundary condition (deduced from namsbc informations)

   !! * Substitutions
#  include "single_precision_substitute.h90"
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
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
      NAMELIST/namoss/ rn_Cd_io, ln_drgice_imp, nn_foss, ln_prs_oce, ln_cpl_oce, ln_ice_embd, ln_ssv_Fgrid
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
         WRITE(numout,*) '      Read SSU@V & SSV@U in prescribed OSS     ln_ssv_Fgrid = ', ln_ssv_Fgrid
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
      !
      !                             !* allocate oss arrays
      IF( oss_nnq_alloc() /= 0 )   CALL ctl_stop( 'oss_init : unable to allocate oss_nnq arrays' )
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
      IF( lk_oasis_oce )   CALL oss_cpl_init( 2 )   ! Must be done before: (1) first time step
      !                                              !                     (2) the use of nn_foss
      !     nn_foss initialization if OCE-SAS coupling via OASIS
      !     SAS time-step has to be declared in OASIS (mandatory) -> nn_foss has to be modified accordingly
      IF( ln_cpl_oce ) THEN
         nn_foss = cpl_freq('I_SFLX') / NINT(rn_Dt)
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
      
      !$acc update device ( rn_Cd_io, ln_drgice_imp, nn_foss, ln_cpl_oce, ln_ssv_Fgrid, ln_ice_embd )


      IF( .NOT. ln_rstart ) THEN
         IF(lwp) WRITE(numout,*)"  'oss_init()' => `ssst`, `sst_s` & `sss_s` initialized with `sst_m`, `sst_m` & `sss_m`!"
         ssst(:,:)  = sst_m(:,:) ! water skin temperature
         sst_s(:,:) = sst_m(:,:) ! slab bulk SST
         sss_s(:,:) = sss_m(:,:) ! slab bulk SSS
         !$acc update device ( ssst, sst_s, sss_s )

      ELSE
         CALL ctl_stop( 'oss_init : FIXME!!! => add restart capability for `ssst`, `sst_s` & `sss_s` ! (-->ossmod.F90)' )         
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
      !!              - updte the ice fraction : fr_i
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER  ::   jj, ji          ! dummy loop argument
      !
      REAL(wp) ::     zthscl        ! wd  tanh scale
      REAL(wp), DIMENSION(jpi,jpj) ::  zwdht, zwght  ! wd dep over wd limit, wgt
      REAL(wp), DIMENSION(jpi,jpj) ::  zus, zvs      ! temporary arrays

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('oss')
      !
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      !                                            !==  oss formulation  ==!
      IF(ln_cpl_oce) THEN
         CALL oss_cpl_rcv( kt, nn_foss, 2 )   ! NANUQ to ocean coupling: NANUQ receiving fields from the ocean model
      ELSE
         CALL oss_prs_rcv( kt )  ! prescribed mean ocean sea surface variables (sst_m, sss_m, ssu_m, ssv_m, etc)
      END IF

      !#lolo: I think that is the correct place to output the SSX field for the ocean we just read or received:
      CALL iom_put( 'sst_m', sst_m )
      CALL iom_put( 'sss_m', sss_m )
      CALL iom_put( 'ssh_m', ssh_m )
      IF( ln_slab_sst .AND. iom_use('mld_m') ) CALL iom_put( 'mld_m', mld_m )
      IF( iom_use('e3t_m') ) CALL iom_put( 'e3t_m', e3t_m )
      IF( iom_use('frq_m') ) CALL iom_put( 'frq_m', frq_m )
      !
      IF( iom_use('frz_m') ) THEN
         CALL eos_fzp_2d( sss_m(:,:), zus(:,:) )
         CALL iom_put( 'frz_m', zus * xmskt )
         zus(:,:) = 0._wp
      ENDIF

# if defined _OPENACC
      !$acc update device ( ssu_m, ssv_m, ssh_m, sst_m, sss_m, frq_m, e3t_m )
      IF( ln_ssv_Fgrid ) THEN
         !$acc update device ( ssu_v_m, ssv_u_m )
      ENDIF
# endif

      ! If prescribed velocities provided at T- rather than U,V- points, interpolate from T to U & T to V
      IF( ln_ssV_T ) THEN
         IF(lwp) PRINT *, ' * [oss@ossmod.F90] => interpolates `SSU,SSV` from T to U,V!', kt
         !$acc data create(zus, zvs) present( ssu_m, ssv_m, umask(:,:,1), vmask(:,:,1), xmskt )
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zus(ji,jj) = ssu_m(ji,jj)
               zvs(ji,jj) = ssv_m(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls-1
            DO ji=Nis0-nn_hls, Nie0+nn_hls-1
               ssu_m(ji,jj) = 0.5_wp*(zus(ji,jj) + zus(ji+1,jj) )* (2._wp-umask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji+1,jj))
               ssv_m(ji,jj) = 0.5_wp*(zvs(ji,jj) + zvs(ji,jj+1)) * (2._wp-vmask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji,jj+1))
            END DO
         END DO
         !$acc end parallel loop
         IF( iom_use('ssu_m') ) THEN
            !$acc update self ( ssu_m )
         ENDIF
         IF( iom_use('ssv_m') ) THEN
            !$acc update self ( ssv_m )
         ENDIF
         !$acc end data
         !
      ENDIF !IF( ln_ssV_T )

      CALL iom_put( 'ssu_m', ssu_m )
      CALL iom_put( 'ssv_m', ssv_m )

      IF( ln_timing )   CALL timing_stop('oss')

   END SUBROUTINE oss

   !!======================================================================
END MODULE ossmod
