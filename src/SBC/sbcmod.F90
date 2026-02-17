MODULE sbcmod
   !!======================================================================
   !!                       ***  MODULE  sbcmod  ***
   !! Surface module :  provide to the ocean its surface boundary condition
   !!======================================================================
   !! History :  3.0  ! 2006-07  (G. Madec)  Original code
   !!            3.1  ! 2008-08  (S. Masson, A. Caubel, E. Maisonnave, G. Madec) coupled interface
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.3  ! 2010-10  (S. Masson)  add diurnal cycle
   !!            3.3  ! 2010-09  (D. Storkey) add ice boundary conditions (BDY)
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!             -   ! 2010-10  (J. Chanut, C. Bricaud, G. Madec)  add the surface pressure forcing
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!            3.5  ! 2012-11  (A. Coward, G. Madec) Rethink of heat, mass and salt surface fluxes
   !!            3.6  ! 2014-11  (P. Mathiot, C. Harris) add ice shelves melting
   !!            4.0  ! 2016-06  (L. Brodeau) new general bulk formulation
   !!            4.0  ! 2019-03  (F. Lemarié & G. Samson)  add ABL compatibility (ln_abl=TRUE)
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) modified wave forcing and coupling
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_init      : read namsbc namelist
   !!   sbc           : surface ocean momentum, heat and freshwater boundary conditions
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY : rn_Dt
   USE phycst,  ONLY : rday         ! physical constants
   USE sbc_oce        ! Surface boundary condition: ocean fields
   !
   USE sbcdcy, ONLY : nday_qsr              ! surface boundary condition: diurnal cycle
   USE sbcflx, ONLY : sbc_flx               ! surface boundary condition: flux formulation
   USE sbcblk, ONLY : sbc_blk_init, sbc_blk ! surface boundary condition: bulk formulation
   USE sbcabl, ONLY : sbc_abl_init, sbc_abl   ! atmospheric boundary layer
   !
   USE cpl_oasis3, ONLY : cpl_freq     ! OASIS routines for coupling
   !
# if ! defined _OPENACC   
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
# endif
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp, ONLY : ncom_fsbc, ctl_nam, ctl_stop, ctl_warn
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc        ! routine called by step.F90
   PUBLIC   sbc_init   ! routine called by opa.F90

   INTEGER, PUBLIC ::   nsbc  ! type of surface boundary condition (deduced from namsbc informations)

   !! * Substitutions
#  include "single_precision_substitute.h90"
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: sbcmod.F90 15372 2021-10-14 15:47:24Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_init()
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_init ***
      !!
      !! ** Purpose :   Initialisation of the ocean surface boundary computation
      !!
      !! ** Method  :   Read the namsbc namelist and set derived parameters
      !!                Call init routines for all other SBC modules that have one
      !!
      !! ** Action  : - read namsbc parameters
      !!              - nsbc: type of sbc
      !!----------------------------------------------------------------------
      INTEGER ::   ios, icpt                         ! local integer
      !!
      NAMELIST/namsbc/ nn_fsbc, ln_flx, ln_blk, ln_abl,                      &
         &             ln_cpl_atm, nn_lsm
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_init : surface boundary condition setting'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      !                       !**  read Surface Module namelist
      READ_NML_REF(numnam,namsbc)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc in reference namelist' )
      READ_NML_CFG(numnam,namsbc)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc in configuration namelist' )
      IF(lwm) WRITE( numond, namsbc )
      !
#if ! defined key_mpi_off
      ncom_fsbc = nn_fsbc    ! make nn_fsbc available for lib_mpp
#endif
      !
      IF(lwp) THEN                  !* Control print
         WRITE(numout,*) '   Namelist namsbc (partly overwritten with CPP key setting)'
         WRITE(numout,*) '      frequency update of sbc (and ice)             nn_fsbc       = ', nn_fsbc
         WRITE(numout,*) '      Type of air-sea fluxes : '
         WRITE(numout,*) '         flux         formulation                   ln_flx        = ', ln_flx
         WRITE(numout,*) '         bulk         formulation                   ln_blk        = ', ln_blk
         WRITE(numout,*) '         ABL          formulation                   ln_abl        = ', ln_abl
         WRITE(numout,*) '      Type of coupling (Ocean/Ice/Atmosphere) : '
         WRITE(numout,*) '         ocean-atmosphere coupled formulation       ln_cpl_atm        = ', ln_cpl_atm
         WRITE(numout,*) '         OASIS coupling                             lk_oasis_atm      = ', lk_oasis_atm
         WRITE(numout,*) '      Misc. options of sbc : '
         WRITE(numout,*) '         nb of iterations if land-sea-mask applied  nn_lsm        = ', nn_lsm
      ENDIF
      !
      IF( MOD( rday , rn_Dt ) /= 0. )   CALL ctl_stop( 'the time step must devide the number of second of in a day' )
      IF( MOD( rday , 2.  ) /= 0. )   CALL ctl_stop( 'the number of second of in a day must be an even number'    )
      IF( MOD( rn_Dt  , 2.  ) /= 0. )   CALL ctl_stop( 'the time step (in second) must be an even number'           )
      !
      !                       !**  check option consistency
      !
      IF(lwp) WRITE(numout,*)

      IF( ln_cpl_atm ) THEN
         CALL ctl_stop( 'sbc : LOLO TODO add the atmo coupling capability!!!')
         IF(lwp) WRITE(numout,*) '   ==>>>  NANUQ is coupled to an atmosphere model!'
         IF( .NOT.lk_oasis_atm )   CALL ctl_stop( 'sbc_init : coupled mode but key_oasis3 disabled' )
      ELSE
         IF(lwp) WRITE(numout,*) '   ==>>>  NANUQ is not coupled to an atmosphere model.'
      ENDIF
      !
      !                             !* sea-ice
      IF( .NOT.( ln_blk .OR. ln_cpl_atm .OR. ln_abl ) )   &
         &                   CALL ctl_stop( 'sbc_init : SI3 sea-ice model requires ln_blk or ln_cpl_atm or ln_abl = T' )
      !
      !                       !**  allocate and set required variables
      !
      !                             !* allocate sbc arrays
      IF( sbc_oce_alloc() /= 0 )   CALL ctl_stop( 'sbc_init : unable to allocate sbc_oce arrays' )
      !
      sfx   (:,:) = 0._wp           !* salt flux due to freezing/melting
      fmmflx(:,:) = 0._wp           !* freezing minus melting flux
      taum  (:,:) = 0._wp           !* wind stress module (needed in GLS in case of reduced restart)

      !                          ! Choice of the Surface Boudary Condition (set nsbc)
      nday_qsr = -1   ! allow initialization at the 1st call !LB: now warm-layer of COARE* calls "sbc_dcy_param" of sbcdcy.F90!
      IF( ln_dm2dc ) THEN           !* daily mean to diurnal cycle
         CALL ctl_stop( 'FIXME!!! Reactivate `ln_dm2dc`!!! See blk_oce_1()@sbcblk.F90 for GPU compliance' )
         !LB:nday_qsr = -1   ! allow initialization at the 1st call
         IF( .NOT.( ln_flx .OR. ln_blk .OR. ln_abl ) )   &
            &   CALL ctl_stop( 'qsr diurnal cycle from daily values requires flux, bulk or abl formulation' )
      ENDIF
      !                             !* Choice of the Surface Boudary Condition
      !                             (set nsbc)
      !
      icpt = 0
      !
      IF( ln_flx ) THEN
         nsbc = jp_flx
         icpt = icpt + 1
      ENDIF       ! flux                 formulation
      IF( ln_blk ) THEN
         nsbc = jp_blk
         icpt = icpt + 1
      ENDIF       ! bulk                 formulation
      IF( ln_abl ) THEN
         nsbc = jp_abl
         icpt = icpt + 1
      ENDIF       ! ABL                  formulation
      IF( ln_cpl_atm ) THEN
         nsbc = jp_cpl_atm
         icpt = icpt + 1
      ENDIF       ! Pure Coupled         formulation
      !
      IF( icpt /= 1 )    CALL ctl_stop( 'sbc_init : choose ONE and only ONE sbc option' )
      !
      IF(lwp) THEN                     !- print the choice of surface flux formulation
         WRITE(numout,*)
         SELECT CASE( nsbc )
         CASE( jp_flx )   ;   WRITE(numout,*) '   ==>>>   flux formulation'
         CASE( jp_blk )   ;   WRITE(numout,*) '   ==>>>   bulk formulation'
         CASE( jp_abl )   ;   WRITE(numout,*) '   ==>>>   ABL  formulation'
         CASE( jp_cpl_atm )   ;   WRITE(numout,*) '   ==>>>   pure coupled formulation'
         END SELECT
      ENDIF
      !
      !                             !* OASIS initialization
      !
      !IF( lk_oasis_atm )   CALL sbc_cpl_init( 2 )   ! Must be done before: (1) first time step
      !                                          !                      (2) the use of nn_fsbc
      !     nn_fsbc initialization if NANUQ-ATMO coupling via OASIS
      !     NANUQ time-step has to be declared in OASIS (mandatory) -> nn_fsbc has to be modified accordingly
      IF( ln_cpl_atm ) THEN !LOLO!!!
         nn_fsbc = cpl_freq('I_SFLX') / NINT(rn_Dt)
         !
         IF(lwp)THEN
            WRITE(numout,*)
            WRITE(numout,*)" * NANUQ is coupled to an atmospheric model via OASIS : nn_fsbc re-defined from OASIS namcouple ", nn_fsbc
            WRITE(numout,*)
         ENDIF
      ENDIF
      !
      !                             !* check consistency between model timeline and nn_fsbc
      IF( ln_rst_list .OR. nn_stock /= -1 ) THEN   ! we will do restart files
         IF( MOD( nitend - nit000 + 1, nn_fsbc) /= 0 ) THEN
            WRITE(ctmp1,*) 'sbc_init : experiment length (', nitend - nit000 + 1, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
         IF( .NOT. ln_rst_list .AND. MOD( nn_stock, nn_fsbc) /= 0 ) THEN   ! we don't use nn_stock if ln_rst_list
            WRITE(ctmp1,*) 'sbc_init : nn_stock (', nn_stock, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
      ENDIF
      !
      IF( MOD( rday, REAL(nn_fsbc, wp) * rn_Dt ) /= 0 )   &
         &  CALL ctl_warn( 'sbc_init : nn_fsbc is NOT a multiple of the number of time steps in a day' )
      !
      IF( ln_dm2dc .AND. NINT(rday) / ( nn_fsbc * NINT(rn_Dt) ) < 8  )   &
         &   CALL ctl_warn( 'sbc_init : diurnal cycle for qsr: the sampling of the diurnal cycle is too small...' )
      !
      !                       !**  associated modules : initialization
      !
      IF( ln_blk      )   CALL sbc_blk_init              ! bulk formulae initialization
      !
      IF( ln_abl      )   CALL sbc_abl_init              ! Atmospheric Boundary Layer (ABL)
      !
   END SUBROUTINE sbc_init


   SUBROUTINE sbc( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc  ***
      !!
      !! ** Purpose :   compute the liquid ocean surface boundary conditions,
      !!                namely the fluxes of momentum, heat and freshwater at
      !!                the top of the liquid ocean interface
      !!
      !! ** Method  :   blah blah  to be written ?????????
      !!                CAUTION : never mask the surface stress field (tke sbc)
      !!
      !! ** Action  : - set the ocean surface boundary condition at before and now
      !!                time step, i.e.
      !!                utau_b, vtau_b, qns_b, qsr_b, emp_n, sfx_b
      !!                utau  , vtau  , qns  , qsr  , emp  , sfx
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      INTEGER  ::   jj, ji          ! dummy loop argument
      REAL(wp) ::     zthscl        ! wd  tanh scale
      REAL(wp), DIMENSION(jpi,jpj) ::  zwdht, zwght  ! wd dep over wd limit, wgt
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('sbc')
      !$acc data present( utau, vtau, qns, emp, sfx, utau_b, vtau_b, qns_b, qsr_b, emp_b, sfx_b )

      !                                            ! ---------------------------------------- !
      IF( kt > nit000 ) THEN                       !          Swap of forcing fields          !
         !$acc parallel loop collapse(2)           ! ---------------------------------------- !
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               utau_b(ji,jj) = utau(ji,jj)                       ! Swap the ocean forcing fields
               vtau_b(ji,jj) = vtau(ji,jj)                       ! (except at nit000 where before fields
               qns_b (ji,jj) = qns (ji,jj)                       !  are set at the end of the routine)
               qsr_b (ji,jj) = qsr (ji,jj)                       !  are set at the end of the routine)               
               emp_b (ji,jj) = emp (ji,jj)
               sfx_b (ji,jj) = sfx (ji,jj)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      !                                            !==  sbc formulation  ==!

      SELECT CASE( nsbc )               ! Compute ocean surface boundary condition
         !                              ! (i.e. utau,vtau, qns, qsr, emp)      !LOLOfixme: NOT `sfx` !!! `sfx` should come from ICE model !!!
      CASE( jp_flx )
# if defined _OPENACC
         CALL ctl_stop( 'sbc : prescribed surface fluxes method `jp_flx` not done yet for GPU! (only bulk for now)')
# endif
         CALL sbc_flx( kt )                        ! flux formulation
         !
         !
      CASE( jp_blk )
         CALL sbc_blk( kt )                        ! bulk formulation for the ocean
         !
      CASE( jp_abl )
# if defined _OPENACC
         CALL ctl_stop( 'sbc : ABL coupling method `jp_abl` not done yet for GPU! (only bulk for now)')
# endif
         CALL sbc_abl( kt )                        ! ABL  formulation for the ocean
         !
         !
      CASE( jp_cpl_atm )                           ! Coupling to an atmosphere GCM
         CALL ctl_stop( 'sbc : CREATE the `sbc_cpl_rcv` routine for atmo GCM coupling!!!')
         !CALL sbc_cpl_rcv   ( kt, nn_fsbc, 2 )  ! pure coupled formulation
         !
         !
      END SELECT

      !! ==> has just updated the following arrays:
      !!       emp, qsr, qns, qns_oce, qsr_oce, fatm_prcp, fatm_snow, wndm, utau, vtau, taum, theta_zu, q_zu, rhoa
      !!     Note: not `sfx` and not `fmmflx` !!!

# if ! defined _OPENACC
      CALL lbc_lnk( 'sbcmod', utau,'T',-1._wp, vtau,'T',-1._wp, emp,'T',1._wp, qns,'T',1._wp, taum,'T',1._wp )
# endif


      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         !LOLOfixme:
         !IF( ln_rstart .AND. .NOT.l_1st_euler ) THEN            !* Restart: read in restart file
         !   IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields read in the restart file'
         !   CALL iom_get( numror, jpdom_auto, 'utau_b', utau_b )   ! i-stress
         !   CALL iom_get( numror, jpdom_auto, 'vtau_b', vtau_b )   ! j-stress
         !   CALL iom_get( numror, jpdom_auto,  'qns_b',  qns_b )   ! non solar heat flux
         !   CALL iom_get( numror, jpdom_auto,  'qsr_b',  qsr_b )   !     solar heat flux
         !   CALL iom_get( numror, jpdom_auto,  'emp_b',  emp_b )   ! freshwater flux
         !   ! To ensure restart capability with 3.3x/3.4 restart files    !! to be removed in v3.6
         !   IF( iom_varid( numror, 'sfx_b', ldstop = .FALSE. ) > 0 ) THEN
         !      CALL iom_get( numror, jpdom_auto, 'sfx_b', sfx_b )  ! before salt flux (T-point)
         !   ELSE
         !      sfx_b (:,:) = sfx(:,:)
         !   ENDIF
         !ELSE                                                   !* no restart: set from nit000 values
         IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields set to nit000'
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               utau_b(ji,jj) = utau(ji,jj)
               vtau_b(ji,jj) = vtau(ji,jj)
               qns_b (ji,jj) = qns (ji,jj)
               qsr_b (ji,jj) = qsr (ji,jj)               
               emp_b (ji,jj) = emp (ji,jj)
               sfx_b (ji,jj) = sfx (ji,jj)  !#LOLOfixme: it should be provided by ICE model, it is not updated by the 3 `sbc_*` routines above!!!
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('sbc')
      !
   END SUBROUTINE sbc


   !!======================================================================
END MODULE sbcmod
