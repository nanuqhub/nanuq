MODULE nanuqgcm
   !!======================================================================
   !!                       ***  MODULE nanuqgcm   ***
   !! StandAlone Surface module : surface fluxes + sea-ice + iceberg floats + ABL
   !!======================================================================
   !! History :  3.6  ! 2011-11  (S. Alderson, G. Madec) original code
   !!             -   ! 2013-06  (I. Epicoco, S. Mocavero, CMCC) northcomms: setup avoiding MPI communication
   !!             -   ! 2014-12  (G. Madec) remove KPP scheme and cross-land advection (cla)
   !!            4.0  ! 2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   nanuq_gcm      : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
   !!   nanuq_init     : initialization of the NANUQ system
   !!   nanuq_ctl      : initialisation of the contol print
   !!   nanuq_closefile: close remaining open files
   !!   nanuq_alloc    : dynamical allocation
   !!----------------------------------------------------------------------
   USE sbcmod         ! surface boundary condition       (sbc     routine)
   USE ossmod         ! surface boundary condition       (sbc     routine)
   USE oss_nnq !lolo, ONLY : lk_oasis_oce        ! surface boundary condition
   USE icestp,  ONLY : ice_init         ! surface boundary condition: SI3 sea-ice model
   USE phycst         ! physical constant                  (par_cst routine)
   USE domain         ! domain initialization   (dom_init & dom_cfg routines)
   USE daymod         ! calendar
   USE step           ! NANUQ time-stepping                 (stp     routine)
   USE cpl_oasis3     !
   USE bdyini         ! open boundary cond. setting       (bdy_init routine). mandatory for sea-ice
   USE bdydta         ! open boundary cond. setting   (bdy_dta_init routine). mandatory for sea-ice
   !
   USE eosbn2, ONLY : eos_init ! equation of state of sea-water
   !
   USE stpctl         ! time stepping control            (stp_ctl routine)
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE in_out_manager ! I/O manager
   USE iom            !
   USE lib_mpp        ! distributed memory computing
   USE mppini         ! shared/distributed memory setting (mpp_init routine)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   USE lbclnk
   USE timing          ! Timing
   USE xios            ! I/O server

   USE halo_mng

#if defined _OPENACC
   USE dom_oce, ONLY : l_1c1g, l_need4ll, klbct, klbcf, klbcu
   USE par_ice, ONLY : jpl, nlay_i, nlay_s, ln_icethd
   !USE ice
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   nanuq_gcm    ! called by model.F90

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

#if ! defined key_mpi_off
   ! need MPI_Wtime
   INCLUDE 'mpif.h'
#endif

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: nanuqgcm.F90 15267 2021-09-17 09:04:34Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nanuq_gcm
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_gcm  ***
      !!
      !! ** Purpose :   NANUQ solves the primitive equations on an orthogonal
      !!              curvilinear mesh on the sphere.
      !!
      !! ** Method  : - model general initialization
      !!              - launch the time-stepping (stp routine)
      !!              - finalize the run by closing files and communications
      !!
      !! References : Madec, Delecluse, Imbard, and Levy, 1997:  internal report, IPSL.
      !!              Madec, 2008, internal report, IPSL.
      !!----------------------------------------------------------------------
      INTEGER ::   istp   ! time step index
      REAL(wp)::   zstptiming   ! elapsed time for 1 time step
      !!----------------------------------------------------------------------

      !                            !-----------------------!
      CALL nanuq_init              !==  Initialisations  ==!
      !                            !-----------------------!

      ! check that all process are still there... If some process have an error,
      ! they will never enter in step and other processes will wait until the end of the cpu time!
      CALL mpp_max( 'nanuqgcm', nstop )

      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA

      !                            !-----------------------!
      !                            !==   time stepping   ==!
      !                            !-----------------------!

#if defined _OPENACC
      l_1c1g = (jpnij==1)  ! => we use 1 CPU core together with 1 GPU
      !
      IF(lwp) WRITE(numout,*) ''
      IF(l_1c1g) THEN
         WRITE(numout,*) ' ********************** OpenACC *************************'
         WRITE(numout,*) ' *** NANUQ will use 1 CPU core and offload on 1 GPU ! ***'
         ! * jperio= 0, landlocked
         ! * jperio= 1, CYCLIC east-west
         ! * jperio= 2, equatorial symmetric (i.e. CYCLIC north-south)
         ! * jperio= 3, north fold WITH T-point pivot
         ! * jperio= 4, CYCLIC east-west and north fold WITH T-point pivot
         ! * jperio= 5, north fold WITH F-point pivot
         ! * jperio= 6, CYCLIC east-west and north fold WITH F-point pivot
         ! * jperio= 7, CYCLIC east-west and north-south
         l_need4ll = (jperio>0)
         !
         IF( l_need4ll ) THEN
            WRITE(numout,'("      => since jperio = ",i1,", `lbc_linking` is required!")') jperio
            WRITE(numout,*) '         * l_Iperio =',l_Iperio
            WRITE(numout,*) '         * l_Jperio =',l_Jperio
         ELSE
            WRITE(numout,'("      =>  domain is landlocked (jperio = ",i1,") => No `lbc_linking` required!")') jperio
         ENDIF
         WRITE(numout,*) ' ********************************************************'
      ELSE
         CALL ctl_stop( 'STOP', 'nanuqgcm : with GPU, no MPP domain decomposition allowed for now!')
      ENDIF
      IF(lwp) WRITE(numout,*) ''
#endif

      !                                               !== set the model time-step  ==!
      !
      istp = nit000
      !
      !
      DO WHILE( istp <= nitend .AND. nstop == 0 )

         ncom_stp = istp
         IF( ln_timing ) THEN
            zstptiming = MPI_Wtime()
            IF ( istp == ( nit000 + 1 ) ) elapsed_time = zstptiming
            IF ( istp ==         nitend ) elapsed_time = zstptiming - elapsed_time
         ENDIF

         CALL stp( istp )
         istp = istp + 1

         IF( lwp .AND. ln_timing )   WRITE(numtime,*) 'timing step ', istp-1, ' : ', MPI_Wtime() - zstptiming

      END DO
      !
      !
      !*acc end data
      !*acc end data
      !*acc end data
      !*acc end data
      !*acc end data
      !********************************************************************************************************************
      !*acc end data
      !********************************************************************************************************************

      !
      !                            !------------------------!
      !                            !==  finalize the run  ==!
      !                            !------------------------!
      IF(lwp) WRITE(numout,cform_aaa)        ! Flag AAAAAAA
      !
      IF( nstop /= 0 .AND. lwp ) THEN        ! error print
         WRITE(ctmp1,*) '   ==>>>   nanuq_gcm: a total of ', nstop, ' errors have been found'
         IF( ngrdstop > 0 ) THEN
            WRITE(ctmp9,'(i2)') ngrdstop
            WRITE(ctmp2,*) '           E R R O R detected in grid '//TRIM(ctmp9)
            WRITE(ctmp3,*) '           Look for "E R R O R" messages in all existing '//TRIM(ctmp9)//'_ocean_output* files'
            CALL ctl_stop( ' ', ctmp1, ' ', ctmp2, ' ', ctmp3 )
         ELSE
            WRITE(ctmp2,*) '           Look for "E R R O R" messages in all existing ocean_output* files'
            CALL ctl_stop( ' ', ctmp1, ' ', ctmp2 )
         ENDIF
      ENDIF
      !
      IF( ln_timing )   CALL timing_finalize
      !
      CALL nanuq_closefile
      !
      CALL xios_finalize  ! end mpp communications with xios
      IF( lk_oasis_oce ) CALL cpl_finalize   ! end coupling and mpp communications with OASIS

      !
      IF(lwm) THEN
         IF( nstop == 0 ) THEN   ;   STOP 0
         ELSE                    ;   STOP 123
         ENDIF
      ENDIF
      !
   END SUBROUTINE nanuq_gcm


   SUBROUTINE nanuq_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_init  ***
      !!
      !! ** Purpose :   initialization of the NANUQ GCM
      !!----------------------------------------------------------------------
      INTEGER ::   ios, ilocal_comm   ! local integers
      !!
      NAMELIST/namctl/ sn_cfctl, ln_timing, ln_diacfl,                                &
         &             nn_isplt,  nn_jsplt,  nn_ictls, nn_ictle, nn_jctls, nn_jctle
      NAMELIST/namcfg/ ln_read_cfg, cn_domcfg, ln_write_cfg, cn_domcfg_out, ln_use_jattr
      !!----------------------------------------------------------------------
      !
      !PRINT *, 'LOLO: entering `nanuq_init()`!, NAREA = ', narea
      !
      cxios_context = 'nanuq'
      !
      nn_hls = 1
      !
      !
      !                             !-------------------------------------------------!
      !                             !     set communicator & select the local rank    !
      !                             !  must be done as soon as possible to get narea  !
      !                             !-------------------------------------------------!
      !
      IF( lk_oasis_oce ) THEN
         CALL cpl_init( "nanuq", ilocal_comm )                                  ! nanuq local communicator given by oasis
         CALL xios_initialize( "not used",local_comm=ilocal_comm )            ! send nanuq communicator to xios
      ELSE
         CALL xios_initialize( "for_xios_mpi_id",return_comm=ilocal_comm )    ! nanuq local communicator given by xios
      ENDIF
      CALL mpp_start( ilocal_comm )
      !
      narea = mpprank + 1                                   ! mpprank: the rank of proc (0 --> mppsize -1 )
      lwm = (narea == 1)                ! control of output namelists
      !
      !                             !---------------------------------------------------------------!
      !                             ! Open output files, reference and configuration namelist files !
      !                             !---------------------------------------------------------------!
      !
      ! open nanuq.output as soon as possible to get all output prints (including errors messages)
      IF( lwm )   CALL ctl_opn(      numout,            'nanuq.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )

      ! open reference and configuration namelist files
      CALL load_nml( numnam_ref,        'namelist_dom_ref',                                           -1, lwm )
      CALL load_nml( numnam_cfg,        'namelist_dom_cfg',                                           -1, lwm )
      IF( lwm )   CALL ctl_opn(      numond, 'output.namelist_dom', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      !ENDIF
      ! open /dev/null file to be able to supress output write easily
      CALL ctl_opn(     numnul,               '/dev/null', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      !
      !                             !--------------------!
      !                             ! Open listing units !  -> need sn_cfctl from namctl to define lwp
      !                             !--------------------!
      !
      READ_NML_REF(numnam,namctl)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namctl in reference namelist' )
      READ_NML_CFG(numnam,namctl)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namctl in configuration namelist' )
      !
      ! finalize the definition of namctl variables
      IF( narea < sn_cfctl%procmin .OR. narea > sn_cfctl%procmax .OR. MOD( narea - sn_cfctl%procmin, sn_cfctl%procincr ) /= 0 )   &
         &   CALL nanuq_set_cfctl( sn_cfctl, .FALSE. )
      !
      lwp = (narea == 1) .OR. sn_cfctl%l_oceout    ! control of all listing output print
      !
      IF(lwp) THEN                      ! open listing units
         !
         IF( .NOT. lwm ) THEN           ! alreay opened for narea == 1
            IF(lk_oasis_oce) THEN; CALL ctl_opn( numout,   'nanuq.output','REPLACE','FORMATTED','SEQUENTIAL',-1,-1, .FALSE., narea )
            ELSE                 ; CALL ctl_opn( numout, 'nanuq.output','REPLACE','FORMATTED','SEQUENTIAL',-1,-1, .FALSE., narea )
            ENDIF
         ENDIF
         !
         WRITE(numout,*)
         WRITE(numout,*) '                     N A N U Q '
         WRITE(numout,*) '          Sea-Ice General Circulation Model'
         WRITE(numout,*) '            A standalone fork of NEMO/SI3'
         WRITE(numout,*) '                version 0.1  (2024) '
         WRITE(numout,*)
         WRITE(numout,*) " "
         WRITE(numout,*) '     .-""-.          ( )-"```"-( )          .-""-.   '
         WRITE(numout,*) "    / O O  \          /         \          /  O O \  "
         WRITE(numout,*) "    |O .-.  \        /   0 _ 0   \        /  .-. O|  "
         WRITE(numout,*) "    \ (   )  '.    _|     (_)     |     .'  (   ) /  "
         WRITE(numout,*) "     '.`-'     '-./ |             |`\.-'     '-'.'   "
         WRITE(numout,*) "       \         |  \   \     /   /  |         /     "
         WRITE(numout,*) "        \        \   '.  '._.'  .'   /        /      "
         WRITE(numout,*) "         \        '.   `'-----'`   .'        /       "
         WRITE(numout,*) "          \   .'    '-._        .-'\   '.   /        "
         WRITE(numout,*) "           |/`          `'''''')    )    `\|         "
         WRITE(numout,*) "           /                  (    (      ,\         "
         WRITE(numout,*) "          ;                    \    '-..-'/ ;        "
         WRITE(numout,*) "          |                     '.       /  |        "
         WRITE(numout,*) "          |                       `'---'`   |        "
         WRITE(numout,*) "          ;                                 ;        "
         WRITE(numout,*) "           \                               /         "
         WRITE(numout,*) "            `.                           .'          "
         WRITE(numout,*) "              '-._                   _.-'            "
         WRITE(numout,*) "        jgs    __/`'  '  - - -  ' '`` \__            "
         WRITE(numout,*) "             /`            /^\           `\          "
         WRITE(numout,*) "             \(          .'   '.         )/          "
         WRITE(numout,*) "              '.(__(__.-'       '.__)__).'           "
         WRITE(numout,*)
         WRITE(numout,*)
         !
         WRITE(numout,cform_aaa)                                        ! Flag AAAAAAA
         !
      ENDIF
      !
      IF(lwm) WRITE( numond, namctl )
      !
      !                             !------------------------------------!
      !                             !  Set global domain size parameters !
      !                             !------------------------------------!
      !
      READ_NML_REF(numnam,namcfg)
      !903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcfg in reference namelist' )
      READ_NML_CFG(numnam,namcfg)
      !904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcfg in configuration namelist' )
      !
      CALL domain_cfg ( cn_cfg, nn_cfg, Ni0glo, Nj0glo, jpkglo, l_Iperio, l_Jperio, l_NFold, c_NFtype )
      !
      IF(lwm)   WRITE( numond, namcfg )
      !
      !                             !-----------------------------------------!
      !                             ! mpp parameters and domain decomposition !
      !                             !-----------------------------------------!
      CALL mpp_init
      !$acc update device( nn_hls, Nis0, Nie0, Njs0, Nje0 )

#if defined key_loop_fusion
      IF( nn_hls == 1 ) THEN
         CALL ctl_stop( 'STOP', 'nanuqgcm : Loop fusion can be used only with extra-halo' )
      ENDIF
      CALL ctl_warn( 'nanuq_init', 'you use key_loop_fusion, this may significantly slow down NANUQ performances' )
#endif

      CALL halo_mng_init()
      ! Now we know the dimensions of the grid and numout has been set: we can allocate arrays
      CALL nanuq_alloc()

      ! Initialise time level indices
      Nbb = 1; Nnn = 2; Naa = 3; Nrhs = Naa

      !                             !-------------------------------!
      !                             !  NANUQ general initialization  !
      !                             !-------------------------------!

      CALL nanuq_ctl                          ! Control prints
      !
      !                                      ! General initialization
      IF( ln_timing    )   CALL timing_init ( 'timing_nanuq.output' )
      IF( ln_timing    )   CALL timing_start( 'nanuq_init')

      CALL phy_cst         ! Physical constants

      CALL eos_init        ! Equation of seawater

      CALL dom_init()      ! Domain

      IF( sn_cfctl%l_prtctl )   &
         &                 CALL prt_ctl_init        ! Print control

      !IF( ln_rstart )      CALL rst_read_open
      CALL day_init()        ! model calendar (using both namelist and restart infos)

      !                                      ! external forcing
      !
      !IF(lwp) WRITE(numout,*) 'LOLO: calling `sbc_init` from `nanuq_init` of nanuqgcm.F90 !'
      CALL sbc_init()  ! Forcings : surface module
      !IF(lwp) WRITE(numout,*) 'LOLO: exiting `sbc_init` from `nanuq_init` of nanuqgcm.F90 !'
      !
      !IF(lwp) WRITE(numout,*) 'LOLO: calling `oss_init` from `nanuq_init` of nanuqgcm.F90 !'
      CALL oss_init()  ! Forcings : bottom module
      !IF(lwp) WRITE(numout,*) 'LOLO: exiting `oss_init` from `nanuq_init` of nanuqgcm.F90 !'
      !
      !
      !IF(lwp) WRITE(numout,*) 'LOLO: calling `ice_init` from `nanuq_init` of nanuqgcm.F90 !'
      CALL ice_init()         ! ICE initialization
      !IF(lwp) WRITE(numout,*) 'LOLO: exiting `ice_init` from `nanuq_init` of nanuqgcm.F90 !'

      ! ==> clem: open boundaries init. is mandatory for sea-ice because ice BDY is not decoupled from
      !           the environment of ocean BDY. Therefore bdy is called by NANUQ
      !           This is not clean and should be changed in the future.
      CALL bdy_init()

      IF(lwp) WRITE(numout,cform_aaa)           ! Flag AAAAAAA
      !
      IF( ln_timing    )   CALL timing_stop( 'nanuq_init')
      !
   END SUBROUTINE nanuq_init


   SUBROUTINE nanuq_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_ctl  ***
      !!
      !! ** Purpose :   control print setting
      !!
      !! ** Method  : - print namctl and namcfg information and check some consistencies
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nanuq_ctl: Control prints'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '                              sn_cfctl%l_runstat = ', sn_cfctl%l_runstat
         WRITE(numout,*) '                              sn_cfctl%l_trcstat = ', sn_cfctl%l_trcstat
         WRITE(numout,*) '                              sn_cfctl%l_oceout  = ', sn_cfctl%l_oceout
         WRITE(numout,*) '                              sn_cfctl%l_layout  = ', sn_cfctl%l_layout
         WRITE(numout,*) '                              sn_cfctl%l_prtctl  = ', sn_cfctl%l_prtctl
         WRITE(numout,*) '                              sn_cfctl%l_prttrc  = ', sn_cfctl%l_prttrc
         WRITE(numout,*) '                              sn_cfctl%l_oasout  = ', sn_cfctl%l_oasout
         WRITE(numout,*) '                              sn_cfctl%procmin   = ', sn_cfctl%procmin
         WRITE(numout,*) '                              sn_cfctl%procmax   = ', sn_cfctl%procmax
         WRITE(numout,*) '                              sn_cfctl%procincr  = ', sn_cfctl%procincr
         WRITE(numout,*) '                              sn_cfctl%ptimincr  = ', sn_cfctl%ptimincr
         WRITE(numout,*) '      timing by routine               ln_timing  = ', ln_timing
         WRITE(numout,*) '      CFL diagnostics                 ln_diacfl  = ', ln_diacfl
      ENDIF
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      read domain configuration file                ln_read_cfg      = ', ln_read_cfg
         WRITE(numout,*) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
         WRITE(numout,*) '      create a configuration definition file        ln_write_cfg     = ', ln_write_cfg
         WRITE(numout,*) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
         WRITE(numout,*) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
      ENDIF
      !
      !LOLOfixme:
      !IF( 1._wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nanuq_ctl: The intrinsec SIGN function follows f2003 standard.',  &
      !   &                                                'Compile with key_nosignedzero enabled:',   &
      !   &                                                '--> add -Dkey_nosignedzero to the definition of %CPP in your arch file' )
      !LOLOfixme.
      !
      !
   END SUBROUTINE nanuq_ctl


   SUBROUTINE nanuq_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!----------------------------------------------------------------------
      !
      IF( lk_mpp )   CALL mppsync
      !
      CALL iom_close                                 ! close all input/output files managed by iom_*
      !
      IF( numstp          /= -1 )   CLOSE( numstp     )   ! time-step file
      IF( numrun          /= -1 )   CLOSE( numrun     )   ! run statistics file
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond     )   ! oce output namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni     )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice )   ! ice variables (temp. evolution)
      IF( numout          /=  6 )   CLOSE( numout     )   ! standard model output file
      !
      numout = 6                                     ! redefine numout in case it is used after this point...
      !
   END SUBROUTINE nanuq_closefile


   SUBROUTINE nanuq_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OCE modules
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE diawri , ONLY : dia_wri_alloc
      USE dom_oce, ONLY : dom_oce_alloc
      USE bdy    , ONLY : ln_bdy, bdy_alloc
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      ierr =        dia_wri_alloc()
      ierr = ierr + dom_oce_alloc()          ! ocean domain
      ierr = ierr + bdy_alloc()          ! bdy masks (incl. initialization)
      !
      CALL mpp_sum( 'nanuqgcm', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nanuq_alloc: unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE nanuq_alloc


   SUBROUTINE nanuq_set_cfctl(sn_cfctl, setto )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nanuq_set_cfctl  ***
      !!
      !! ** Purpose :   Set elements of the output control structure to setto.
      !!
      !! ** Method  :   Note this routine can be used to switch on/off some
      !!                types of output for selected areas.
      !!----------------------------------------------------------------------
      TYPE(sn_ctl), INTENT(inout) :: sn_cfctl
      LOGICAL     , INTENT(in   ) :: setto
      !!----------------------------------------------------------------------
      sn_cfctl%l_runstat = setto
      sn_cfctl%l_trcstat = setto
      sn_cfctl%l_oceout  = setto
      sn_cfctl%l_layout  = setto
      sn_cfctl%l_prtctl  = setto
      sn_cfctl%l_prttrc  = setto
      sn_cfctl%l_oasout  = setto
   END SUBROUTINE nanuq_set_cfctl

   !!======================================================================
END MODULE nanuqgcm
