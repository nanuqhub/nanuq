!! What NANUQ receives from the ocean model
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! * O_SSTSST   1          jpr_sst
! * O_SSSal    2          jpr_sss
! * O_OCurx1   3          jpr_ssu
! * O_OCury1   4          jpr_ssv
! * O_SSHght   5          jpr_ssh
! * O_E3T1st   6          jpr_e3t
! * O_FraQsr   7          jpr_frq
!   => 7 fields

!! What NANUQ sends to the ocean model:
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! * I_OTaux1   1          jps_otx1
! * I_OTauy1   2          jps_oty1
! * I_QnsOce   3          jps_qnsoce
! * I_QsrOce   4          jps_qsroce
! * IOEvaMPr   5          jps_oemp
! * I_SFLX     6          jps_sflx
! * I_TauMod   7          jps_taum
! * IIceFrc    8          jps_fice2
!   => 8 fields


MODULE osscpl
   !!======================================================================
   !!                       ***  MODULE  osscpl  ***
   !! Bottom Boundary Condition
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!            3.4  ! 2011_11  (C. Harris) more flexibility + multi-category fields
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi)  wave coupling updates
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namoss_cpl      : coupled formulation namlist
   !!   oss_cpl_init    : initialisation of the coupled exchanges
   !!   oss_cpl_rcv     : receive fields from the ocean over the ocean (ocean only)
   !!   oss_cpl_snd     : send     fields to the ocean
   !!----------------------------------------------------------------------
   !USE dom_oce, ONLY: rn_Dt, narea, idbg, jdbg !LOLOdbg
   USE dom_oce, ONLY: rn_Dt
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE ice    , ONLY: at_i ! sea-ice fraction
   USE oss_nnq        ! Surface boundary condition: ocean fields
   USE cpl_oasis3     ! OASIS3 coupling
   !USE eosbn2
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library

#if defined key_oasis3
   USE mod_oasis, ONLY : OASIS_Sent, OASIS_ToRest, OASIS_SentOut, OASIS_ToRestOut
#endif

   USE iom            ! IOM library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss_cpl_init      ! routine called by sbcmod.F90
   PUBLIC   oss_cpl_rcv       ! routine called by icestp.F90
   PUBLIC   oss_cpl_snd       ! routine called by step.F90


   !! Received:
   INTEGER, PARAMETER ::   jpr_sst = 1   ! ocean temperature
   INTEGER, PARAMETER ::   jpr_sss = 2   ! ocean salinity
   INTEGER, PARAMETER ::   jpr_ssu = 3   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jpr_ssv = 4   !
   INTEGER, PARAMETER ::   jpr_ssh = 5   ! sea surface height
   INTEGER, PARAMETER ::   jpr_e3t = 6   ! first T level thickness
   INTEGER, PARAMETER ::   jpr_frq = 7   ! fraction of solar net radiation absorbed in the first ocean level
   !
   INTEGER, PARAMETER ::   jprcv   = 7   ! total number of fields received

   !! Sent:
   INTEGER, PARAMETER ::   jps_qsroce = 1   ! Qsr above the ocean
   INTEGER, PARAMETER ::   jps_qnsoce = 2   ! Qns above the ocean
   INTEGER, PARAMETER ::   jps_oemp   = 3   ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jps_sflx   = 4   ! salt flux
   INTEGER, PARAMETER ::   jps_otx1   = 5   ! 2 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jps_oty1   = 6   !
   INTEGER, PARAMETER ::   jps_taum   = 7   ! wind stress module
   INTEGER, PARAMETER ::   jps_fice2  = 8   ! ice fraction sent to OCE (by NANUQ when doing NANUQ-OCE coupling)
   !
   INTEGER, PARAMETER ::   jpsnd      = 8   ! total number of fields sent


#if ! defined key_oasis3
   ! Dummy variables to enable compilation when oasis3 is not being used
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1
#endif

   !                                  !!** namelist namoss_cpl **
   TYPE ::   FLD_C                     !
      CHARACTER(len = 32) ::   cldes      ! desciption of the coupling strategy
      CHARACTER(len = 32) ::   clcat      ! multiple ice categories strategy
      CHARACTER(len = 32) ::   clvref     ! reference of vector ('spherical' or 'cartesian')
      CHARACTER(len = 32) ::   clvor      ! orientation of vector fields ('eastward-northward' or 'local grid')
      CHARACTER(len = 32) ::   clvgrd     ! grids on which is located the vector fields
   END TYPE FLD_C
   !
   INTEGER     ::   nn_cplmodel           ! Maximum number of models to/from which NANUQ is potentialy sending/receiving data
   !LOGICAL     ::   ln_scale_ice_flux     !  use ice fluxes that are already "ice weighted" ( i.e. multiplied ice concentration)

   TYPE ::   DYNARR
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   z3
   END TYPE DYNARR

   TYPE( DYNARR ), SAVE, DIMENSION(jprcv) ::   frcv                ! all fields recieved from the ocean

   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i_last_couple !: Ice fractional area at last coupling time

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   nrcvinfo           ! OASIS info argument

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: xcplmask

   !! Substitution
#  include "single_precision_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: osscpl.F90 15551 2021-11-28 20:19:36Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION oss_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION oss_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(1)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( nrcvinfo(jprcv), xcplmask(jpi,jpj,0:nn_cplmodel), STAT=ierr(1) )
      !
      oss_cpl_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'osscpl', oss_cpl_alloc )
      IF( oss_cpl_alloc > 0 )   CALL ctl_warn('oss_cpl_alloc: allocation of arrays failed')
      !
   END FUNCTION oss_cpl_alloc




   SUBROUTINE oss_cpl_init( k_ice )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE oss_cpl_init  ***
      !!
      !! ** Purpose :   Initialisation of send and received information from
      !!                the atmospheric component
      !!
      !! ** Method  : * Read namoss_cpl namelist
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_ice   ! ice management in the sbc (=0/1/2/3)
      !
      INTEGER ::   jn          ! dummy loop index
      INTEGER ::   ios, inum   ! Local integer
      REAL(wp), DIMENSION(jpi,jpj) ::   zacs, zaos
      !!

      !! #lolo: What actually needs to be defined:
      nn_cplmodel = 1


      !
      !                                   Allocate osscpl arrays
      !                                   **********************
      IF( oss_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'oss_cpl_alloc : unable to allocate arrays' )


      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress

      ! for each field: define the OASIS name                              (srcv(:)%clname)
      !                 define receive or not from the namelist parameters (srcv(:)%laction)
      !                 define the north fold type of lbc                  (srcv(:)%nsgn)

      ! default definitions of srcv
      srcv(:)%laction = .FALSE.   ;   srcv(:)%clgrid = 'T'   ;   srcv(:)%nsgn = 1.   ;   srcv(:)%nct = 1


      !                                                      ! ------------------------------------- !
      !                                                      !   OASIS coupling - received by NANUQ  !
      !                                                      ! ------------------------------------- !
      srcv(jpr_sst)%clname = 'I_SSTSST'
      srcv(jpr_sss)%clname = 'I_SSSal'
      srcv(jpr_ssu)%clname = 'I_OCurx1'
      srcv(jpr_ssv)%clname = 'I_OCury1'
      srcv(jpr_ssh)%clname = 'I_SSHght'
      srcv(jpr_e3t)%clname = 'I_E3T1st'
      srcv(jpr_frq)%clname = 'I_FraQsr'
      !
      srcv(:)%laction = .FALSE.   ! force default definition in case of OceanModel <-> NANUQ coupling
      srcv(:)%clgrid  = 'T'       ! force default definition in case of OceanModel <-> NANUQ coupling
      srcv(:)%nsgn    = 1.        ! force default definition in case of OceanModel <-> NANUQ coupling
      !
      srcv( (/ jpr_sst, jpr_sss, jpr_ssh, jpr_frq, jpr_ssu, jpr_ssv, jpr_e3t /) )%laction = .TRUE.
      !
      srcv(jpr_ssu)%clgrid = 'U'        ! oce components given at U-point
      srcv(jpr_ssv)%clgrid = 'V'        !           and           V-point
      ! Vectors: change of sign at north fold ONLY if on the local grid
      srcv(jpr_ssu:jpr_ssv)%nsgn = -1.
      ! Change first letter to couple with ocean if already coupled OCE
      ! this is nedeed as each variable name used in the namcouple must be unique:
      ! for example O_Runoff received by OCE from NANUQ and therefore S_Runoff received by NANUQ from the Ocean
      DO jn = 1, jprcv
         IF( srcv(jn)%clname(1:1) == "O" ) srcv(jn)%clname = "S"//srcv(jn)%clname(2:LEN(srcv(jn)%clname))
      END DO
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'               Special conditions for NANUQ-OCE coupling  '
         WRITE(numout,*)'               NANUQ component  '
         WRITE(numout,*)
         WRITE(numout,*)'  received (7) fields from OCE component '
         WRITE(numout,*)'               sea surface temperature (Celsius) '
         WRITE(numout,*)'               sea surface salinity '
         WRITE(numout,*)'               surface currents '
         WRITE(numout,*)'               sea surface height '
         WRITE(numout,*)'               thickness of first ocean T level '
         WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
         WRITE(numout,*)
      ENDIF


      ! =================================================== !
      ! Allocate all parts of frcv used for received fields !
      ! =================================================== !
      DO jn = 1, jprcv
         IF( srcv(jn)%laction ) THEN
            !PRINT *, ' LOLO: allocating `frcv%z3`, 3rd dim =', srcv(jn)%nct, 'for jn=',jn
            ALLOCATE( frcv(jn)%z3(jpi,jpj,srcv(jn)%nct) )
         ENDIF
      END DO

      !LULU



      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (ssnd(:)%clname)
      !                 define send or not from the namelist parameters (ssnd(:)%laction)
      !                 define the north fold type of lbc               (ssnd(:)%nsgn)

      ! default definitions of nsnd
      ssnd(:)%laction = .FALSE.   ;   ssnd(:)%clgrid = 'T'   ;   ssnd(:)%nsgn = 1.  ; ssnd(:)%nct = 1
      !
      !
      !
      !                                                      ! -------------------------------- !
      !                                                      !   OASIS coupling - sent by NANUQ !
      !                                                      ! -------------------------------- !
      ssnd(jps_sflx  )%clname = 'I_SFLX'
      ssnd(jps_qsroce)%clname = 'I_QsrOce'
      ssnd(jps_qnsoce)%clname = 'I_QnsOce'
      ssnd(jps_oemp  )%clname = 'IOEvaMPr'
      ssnd(jps_otx1  )%clname = 'I_OTaux1'
      ssnd(jps_oty1  )%clname = 'I_OTauy1'
      ssnd(jps_taum  )%clname = 'I_TauMod'
      ssnd(jps_fice2 )%clname = 'IIceFrc'
      !
      IF( .NOT. ln_cpl_oce ) ssnd(:)%laction = .FALSE.   ! force default definition in case of OceanModel <-> NANUQ coupling
      ssnd( (/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/) )%laction = .TRUE.
      !
      ! Change first letter to couple with ocean if already coupled with sea-ice
      ! this is nedeed as each variable name used in the namcouple must be unique:
      ! for example O_SSTSST sent by OCE to NANUQ and therefore S_SSTSST sent by NANUQ to the Ocean
      DO jn = 1, jpsnd
         IF( ssnd(jn)%clname(1:1) == "O" ) ssnd(jn)%clname = "S"//ssnd(jn)%clname(2:LEN(ssnd(jn)%clname))
      END DO
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'  sent (8) fields to OCE component '
         WRITE(numout,*)'                  ice cover '
         WRITE(numout,*)'                  oce only EMP  '
         WRITE(numout,*)'                  salt flux  '
         WRITE(numout,*)'                  mixed oce-ice solar flux  '
         WRITE(numout,*)'                  mixed oce-ice non solar flux  '
         WRITE(numout,*)'                  wind stress U,V components'
         WRITE(numout,*)'                  wind stress module'
      ENDIF

      !
      ! ================================ !
      !   initialisation of the coupler  !
      ! ================================ !
      !
      nn_cplmodel = 1   !#lolo: fix this the day couling with AGCM is considered!!!
      CALL cpl_define(jprcv, jpsnd, nn_cplmodel)
      !
      xcplmask(:,:,:) = 1.
      xcplmask(:,:,0) = 1. - SUM( xcplmask(:,:,1:nn_cplmodel), dim = 3 )
      !
   END SUBROUTINE oss_cpl_init



   SUBROUTINE oss_cpl_rcv( kt, k_fsbc, k_ice )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE oss_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in frcv array). called at each time step.
      !!                OASIS controls if there is something do receive or not. nrcvinfo contains the info
      !!                to know if the field was really received or not
      !!
      !!              --> If ocean stress was really received:
      !!
      !!                  - transform the received ocean stress vector from the received
      !!                 referential and grid into an atmosphere-ocean stress in
      !!                 the (i,j) ocean referencial and at the ocean velocity point.
      !!                    The received stress are :
      !!                     - defined by 3 components (if cartesian coordinate)
      !!                            or by 2 components (if spherical)
      !!                     - oriented along geographical   coordinate (if eastward-northward)
      !!                            or  along the local grid coordinate (if local grid)
      !!                     - given at U- and V-point, resp.   if received on 2 grids
      !!                            or at T-point               if received on 1 grid
      !!                    Therefore and if necessary, they are successively
      !!                  processed in order to obtain them
      !!                     first  as  2 components on the sphere
      !!                     second as  2 components oriented along the local grid
      !!                     third  as  2 components on the U,V grid
      !!
      !!              -->
      !!
      !!              - In 'ocean only' case, non solar and solar ocean heat fluxes
      !!             and total ocean freshwater fluxes
      !!
      !! ** Method  :   receive all fields from the atmosphere and transform
      !!              them into ocean surface boundary condition fields
      !!
      !! ** Action  :   update  utau, vtau   ocean stress at U,V grid
      !!                        taum         wind stress module at T-point
      !!                        wndm         wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!                        qns          non solar heat fluxes including emp heat content    (ocean only case)
      !!                                     and the latent heat flux of solid precip. melting
      !!                        qsr          solar ocean heat fluxes   (ocean only case)
      !!                        emp          upward mass flux [evap. - precip. (- runoffs) (- calving)] (ocean only case)
      !!----------------------------------------------------------------------
      !USE zdf_oce,  ONLY :   ln_zdfswm
      !
      INTEGER, INTENT(in) ::   kt          ! ocean model time step index
      INTEGER, INTENT(in) ::   k_fsbc      ! frequency of sbc (-> ice model) computation
      INTEGER, INTENT(in) ::   k_ice       ! ice management in the sbc (=0/1/2/3)
      !!
      LOGICAL  ::   llnewtx, llnewtau      ! update wind stress components and module??
      INTEGER  ::   ji, jj, jn             ! dummy loop indices
      INTEGER  ::   isec                   ! number of seconds since nit000 (assuming rdt did not change since nit000)
      REAL(wp) ::   zcumulneg, zcumulpos   ! temporary scalars
      REAL(wp) ::   zcoef                  ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3        ! drag coefficient
      REAL(wp) ::   zzx, zzy               ! temporary variables
      REAL(wp) ::   r1_grau                ! = 1.e0 / (grav * rho0)
      !REAL(wp), DIMENSION(jpi,jpj) ::   ztx, zty, zmsk, zemp, zqns, zqsr
      !!----------------------------------------------------------------------
      !
      !                                                      ! ======================
      !                                                      ! Receive all the fields
      !                                                      ! ======================
      isec = ( kt - nit000 ) * NINT( rn_Dt )                    ! date of exchanges
      DO jn = 1, jprcv                                          ! received fields sent by the ocean
         IF( srcv(jn)%laction )  THEN
            CALL cpl_rcv( jn, isec, frcv(jn)%z3, xcplmask(:,:,1:nn_cplmodel), nrcvinfo(jn) )
         ENDIF
      END DO
      !
      !
      !                                                      ! ================== !
      !    1                                                 !        SST         !
      !                                                      ! ================== !
      IF( srcv(jpr_sst)%laction ) THEN                      ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         sst_m(:,:) = frcv(jpr_sst)%z3(:,:,1)
         !IF( srcv(jpr_sss)%laction .AND. l_useCT ) THEN    ! make sure that sst_m is the potential temperature
         !   sst_m(:,:) = eos_pt_from_ct( sst_m(:,:), sss_m(:,:) )
         !ENDIF
         !$acc update device (sst_m)
      ENDIF
      !
      !                                                      ! ================== !
      !    2                                                 !        SSS         !
      !                                                      ! ================== !
      IF( srcv(jpr_sss)%laction ) THEN                      ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         sss_m(:,:) = frcv(jpr_sss)%z3(:,:,1)
         !$acc update device (sss_m)
      ENDIF
      !                                                      ! ================== !
      !                                                      !        SSH         !
      !    3                                                 ! ================== !
      IF( srcv(jpr_ssh )%laction ) THEN                      ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         ssh_m(:,:) = frcv(jpr_ssh )%z3(:,:,1)
         !$acc update device (ssh_m)
      ENDIF
      !                                                      ! ================== !
      !                                                      !  surface currents  !
      !    4, 5                                              ! ================== !
      IF( srcv(jpr_ssu)%laction ) THEN                      ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         ssu_m(:,:) = frcv(jpr_ssu)%z3(:,:,1)
         !$acc update device (ssu_m)
      ENDIF
      IF( srcv(jpr_ssv)%laction ) THEN
         ssv_m(:,:) = frcv(jpr_ssv)%z3(:,:,1)
         !$acc update device (ssv_m)
      ENDIF
      !                                                      ! ======================== !
      !                                                      !  first T level thickness !
      !   6                                                  ! ======================== !
      IF( srcv(jpr_e3t )%laction ) THEN                   ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         e3t_m(:,:) = frcv(jpr_e3t )%z3(:,:,1)
         !$acc update device (e3t_m)
      ENDIF
      !                                                      ! ================================ !
      !                                                      !  fraction of solar net radiation !
      !    7                                                 ! ================================ !
      IF( srcv(jpr_frq)%laction ) THEN                    ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         frq_m(:,:) = frcv(jpr_frq)%z3(:,:,1)
         !$acc update device (frq_m)
      ENDIF
      !
   END SUBROUTINE oss_cpl_rcv



   SUBROUTINE oss_cpl_snd( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE oss_cpl_snd  ***
      !!
      !! ** Purpose :   provide the ocean-ice informations to the atmosphere
      !!
      !! ** Method  :   send to the atmosphere through a call to cpl_snd
      !!              all the needed fields (as defined in oss_cpl_init)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      !
      INTEGER ::   ji, jj, jl   ! dummy loop indices
      INTEGER ::   isec, info   ! local integer
      !!----------------------------------------------------------------------

      isec = ( kt - nit000 ) * NINT( rn_Dt )        ! date of exchanges
      info = OASIS_idle

      !  Fields sent by NANUQ to OCE when OASIS coupling

      ! i-component of surface stress (aka flux of momentum):
      IF( ssnd(jps_otx1  )%laction ) THEN
         !$acc update self (utau)
         CALL cpl_snd( jps_otx1  , isec, RESHAPE ( utau, (/jpi,jpj,1/) ), info ) ! 1
      ENDIF

      ! j-component of surface stress (aka flux of momentum):
      IF( ssnd(jps_oty1  )%laction ) THEN
         !$acc update self (vtau)
         CALL cpl_snd( jps_oty1  , isec, RESHAPE ( vtau, (/jpi,jpj,1/) ), info ) ! 2
      ENDIF

      ! Non-solar heat flux:
      IF( ssnd(jps_qnsoce)%laction ) THEN
         !$acc update self (qns)
         CALL cpl_snd( jps_qnsoce, isec, RESHAPE ( qns , (/jpi,jpj,1/) ), info ) ! 3
      ENDIF

      ! Solar heat flux:
      IF( ssnd(jps_qsroce)%laction ) THEN
         !$acc update self (qsr)
         CALL cpl_snd( jps_qsroce, isec, RESHAPE ( qsr , (/jpi,jpj,1/) ), info ) ! 4
      ENDIF

      ! Net freshwater flux (excluding continental runoffs):
      IF( ssnd(jps_oemp  )%laction ) THEN
         !$acc update self (emp)
         CALL cpl_snd( jps_oemp  , isec, RESHAPE ( emp , (/jpi,jpj,1/) ), info ) ! 5
      ENDIF

      ! Salt flux:
      IF( ssnd(jps_sflx  )%laction ) THEN
         !$acc update self (sfx)
         CALL cpl_snd( jps_sflx  , isec, RESHAPE ( sfx , (/jpi,jpj,1/) ), info ) ! 6
      ENDIF

      ! Modulus of surface stress:
      IF( ssnd(jps_taum  )%laction ) THEN
         !$acc update self (taum)
         CALL cpl_snd( jps_taum  , isec, RESHAPE ( taum, (/jpi,jpj,1/) ), info ) ! 7
      ENDIF

      ! Sea-ice concentration:
      IF( ssnd(jps_fice2 )%laction ) THEN
         !$acc update self (at_i)
         CALL cpl_snd( jps_fice2,  isec, RESHAPE ( at_i, (/jpi,jpj,1/) ), info ) ! 8
      ENDIF


      !! The following is certainly a "doublon" with what's done in `sbc()@sbcmod.F90`, but I want it
      !! to be right here to be sure about what fluxes are passed to OASIS...

      IF( kt == nit000 ) PRINT *, 'LOLO `oss_cpl_snd()@osscpl.F90` => COUPLED: potentially IOMing the liquid ocean SBC fluxes !, kt =', kt

      IF( iom_use('utau_oa3' ) ) THEN
         !$acc update self( utau )
         CALL iom_put( "utau_oa3"   ,  utau )
      ENDIF

      IF( iom_use('vtau_oa3' ) ) THEN
         !$acc update self( vtau )
         CALL iom_put( "vtau_oa3"   ,  vtau )
      ENDIF

      IF( iom_use('qns_oa3' ) ) THEN
         !$acc update self( qns )
         CALL iom_put( "qns_oa3"    ,  qns  )
      ENDIF

      IF( iom_use('qsr_oa3' ) ) THEN
         !$acc update self( qsr )
         CALL iom_put( "qsr_oa3"    ,  qsr  )
      ENDIF

      IF( iom_use('emp_oa3' ) ) THEN
         !$acc update self( emp )
         CALL iom_put( "emp_oa3"    ,  emp  )
      ENDIF

      IF( iom_use('sfx_oa3' ) ) THEN
         !$acc update self( sfx )
         CALL iom_put( "sfx_oa3",  sfx  )
      ENDIF

      IF( iom_use('taum_oa3' ) ) THEN
         !$acc update self( taum )
         CALL iom_put( "taum_oa3",  taum  )
      ENDIF

      IF( iom_use('fr_i_oa3' ) ) THEN
         !$acc update self( at_i )
         CALL iom_put( "fr_i_oa3", at_i )
      ENDIF


   END SUBROUTINE oss_cpl_snd

   !!======================================================================
END MODULE osscpl
