# define smyrcv srcv(midcpl)%fld  /* alias to have a shorter name */
# define smysnd ssnd(midcpl)%fld  /* alias to have a shorter name */

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

   !! What NANUQ receives from the ocean model
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! * O_SSTSST   1          jpr_sst
   ! * O_SSSal    2          jpr_sss
   ! * O_OCurx1   3          jpr_ssu --  @ U-points
   ! * O_OCury1   4          jpr_ssv --  @ V-points
   ! * O_SSHght   5          jpr_ssh
   ! * O_E3T1st   6          jpr_e3t
   ! * O_FraQsr   7          jpr_frq
   !   => 7 fields

   !! What NANUQ sends to the ocean model
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! * I_OTaux1   1          jps_otx1 -- IF(sn_loc_vct_tau=='C') => @ U-points | IF(sn_loc_vct_tau=='T') => @ T-points
   ! * I_OTauy1   2          jps_oty1 -- IF(sn_loc_vct_tau=='C') => @ V-points | IF(sn_loc_vct_tau=='T') => @ T-points
   ! * I_QnsOce   3          jps_qnsoce
   ! * I_QsrOce   4          jps_qsroce
   ! * IOEvaMPr   5          jps_oemp
   ! * I_SFLX     6          jps_sflx
   ! * I_TauMod   7          jps_taum
   ! * IIceFrc    8          jps_fice2
   !   => 8 fields



   !!----------------------------------------------------------------------
   !!   namoss_cpl      : coupled formulation namlist
   !!   oss_cpl_init    : initialisation of the coupled exchanges
   !!   oss_cpl_rcv     : receive surface state fields from the ocean
   !!   oss_cpl_snd     : send surface fluxes for the liqui ocean
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
   USE iom            ! NetCDF library
   USE lib_mpp        ! distribued memory computing library

#if defined key_oasis3
   USE mod_oasis, ONLY : OASIS_Sent, OASIS_ToRest, OASIS_SentOut, OASIS_ToRestOut
#endif


   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss_cpl_init      ! routine called by sbcmod.F90
   PUBLIC   oss_cpl_rcv       ! routine called by icestp.F90
   PUBLIC   oss_cpl_snd       ! routine called by step.F90
   PUBLIC   oss_cpl_alloc     ! routine called in sbcice_cice.F90


   !! Received:
   INTEGER, PARAMETER ::   jpr_sst = 1   ! ocean temperature
   INTEGER, PARAMETER ::   jpr_sss = 2   ! ocean salinity
   INTEGER, PARAMETER ::   jpr_ssu = 3   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jpr_ssv = 4   !
   INTEGER, PARAMETER ::   jpr_ssh = 5   ! sea surface height
   INTEGER, PARAMETER ::   jpr_e3t = 6   ! first T level thickness
   INTEGER, PARAMETER ::   jpr_frq = 7   ! fraction of solar net radiation absorbed in the first ocean level

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

   INTEGER, PARAMETER ::   jpsnd      = 8   ! total number of fields sent


#if ! defined key_oasis3
   ! Dummy variables to enable compilation when oasis3 is not being used
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1
#endif

   !                                  !!** namelist namoss_cpl **
   TYPE ::  FLD_C                     !
      CHARACTER(len = 32) ::   cldes      ! desciption of the coupling strategy
      CHARACTER(len = 32) ::   clcat      ! multiple ice categories strategy
      CHARACTER(len = 32) ::   clvref     ! reference of vector ('spherical' or 'cartesian')
      CHARACTER(len = 32) ::   clvor      ! orientation of vector fields ('eastward-northward' or 'local grid')
      CHARACTER(len = 32) ::   clvgrd     ! grids on which is located the vector fields
   ENDTYPE FLD_C
   !
   INTEGER     ::   nn_cplmodel=1          ! Maximum number of models to/from which NANUQ is potentialy sending/receiving data
   !LOGICAL     ::   ln_scale_ice_flux     !  use ice fluxes that are already "ice weighted" ( i.e. multiplied ice concentration)

   TYPE ::  DYNARR
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   z3
   ENDTYPE DYNARR

   TYPE( DYNARR ), SAVE, DIMENSION(jprcv) ::   frcv                ! all fields recieved from the ocean

   !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i_last_couple !: Ice fractional area at last coupling time

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   nrcvinfo           ! OASIS info argument

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: xcplmask

   !! Substitution
#  include "single_precision_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION oss_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION oss_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( nrcvinfo(jprcv), STAT=ierr(1) )

      ALLOCATE( xcplmask(Nis0:Nie0,Njs0:Nje0,1,0:nn_cplmodel) , STAT=ierr(2) )
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
      REAL(wp), DIMENSION(Nis0:Nie0,Njs0:Nje0) :: zacs, zaos
      !!

      !!---------------------------------------------------------------------


      !! #lolo: What actually needs to be defined:
      nn_cplmodel = 1


      !
      !                                   ! allocate osscpl arrays
      IF( oss_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'oss_cpl_alloc : unable to allocate arrays' )


      IF( ln_cpl_oce_croco .AND. sn_loc_vct_tau=='T' ) &
         &   CALL ctl_stop( 'STOP', 'oss_cpl_init : you cannot have "sn_loc_vct_tau=T" with "ln_cpl_oce_croco=.true."' )
      !!    => because CROCO expects surface stress components to be defined at U- and V-points, not at T-points !!!



      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress

      ! for each field: define the OASIS name                              (smyrcv(:)%clname)
      !                 define receive or not from the namelist parameters (smyrcv(:)%laction)
      !                 define the north fold type of lbc                  (smyrcv(:)%nsgn)

      ! default definitions of srcv
      ALLOCATE( smyrcv(jprcv) )
      smyrcv(:)%laction = .FALSE.   ;   smyrcv(:)%clgrid = 'T'   ;   smyrcv(:)%nsgn = 1.
      smyrcv(:)%nct     = 1         ;   smyrcv(:)%nlvl   = 1     ;   smyrcv(:)%ncplmodel = nn_cplmodel


      !                                                      ! ------------------------------------- !
      !                                                      !   OASIS coupling - received by NANUQ  !
      !                                                      ! ------------------------------------- !
      smyrcv(jpr_sst)%clname = 'I_SSTSST'
      smyrcv(jpr_sss)%clname = 'I_SSSal'
      smyrcv(jpr_ssu)%clname = 'I_OCurx1'
      smyrcv(jpr_ssv)%clname = 'I_OCury1'
      smyrcv(jpr_ssh)%clname = 'I_SSHght'
      smyrcv(jpr_e3t)%clname = 'I_E3T1st'
      smyrcv(jpr_frq)%clname = 'I_FraQsr'
      !
      smyrcv(:)%laction = .FALSE.   ! force default definition in case of OceanModel <-> NANUQ coupling
      smyrcv(:)%clgrid  = 'T'       ! force default definition in case of OceanModel <-> NANUQ coupling
      smyrcv(:)%nsgn    = 1.        ! force default definition in case of OceanModel <-> NANUQ coupling
      !
      smyrcv( (/ jpr_sst, jpr_sss, jpr_ssh, jpr_frq, jpr_ssu, jpr_ssv, jpr_e3t /) )%laction = .TRUE.
      !
      smyrcv(jpr_ssu)%clgrid = 'U'        ! oce components given at U-point
      smyrcv(jpr_ssv)%clgrid = 'V'        !           and           V-point
      ! Vectors: change of sign at north fold ONLY if on the local grid
      smyrcv(jpr_ssu:jpr_ssv)%nsgn = -1.
      ! Change first letter to couple with ocean if already coupled OCE
      ! this is nedeed as each variable name used in the namcouple must be unique:
      ! for example O_Runoff received by OCE from NANUQ and therefore S_Runoff received by NANUQ from the Ocean
      DO jn = 1, jprcv
         IF( smyrcv(jn)%clname(1:1) == "O" ) smyrcv(jn)%clname = "S"//smyrcv(jn)%clname(2:LEN(smyrcv(jn)%clname))
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

      ! ===================================================== !
      ! Allocate all parts of smyrcv used for received fields !
      ! ===================================================== !
      DO jn = 1, jprcv
         IF( smyrcv(jn)%laction ) ALLOCATE( smyrcv(jn)%z3(Nis0:Nie0,Njs0:Nje0,smyrcv(jn)%nct) )
      END DO



      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (smysnd(:)%clname)
      !                 define send or not from the namelist parameters (smysnd(:)%laction)
      !                 define the north fold type of lbc               (smysnd(:)%nsgn)

      ! default definitions of nsnd
      ALLOCATE( smysnd(jpsnd) )
      smysnd(:)%laction = .FALSE.   ;   smysnd(:)%clgrid = 'T'   ;   smysnd(:)%nsgn      = 1.
      smysnd(:)%nct     = 1         ;   smysnd(:)%nlvl   = 1     ;   smysnd(:)%ncplmodel = nn_cplmodel


      !                                                      ! -------------------------------- !
      !                                                      !   OASIS coupling - sent by NANUQ !
      !                                                      ! -------------------------------- !
      smysnd(jps_sflx  )%clname = 'I_SFLX'
      smysnd(jps_qsroce)%clname = 'I_QsrOce'
      smysnd(jps_qnsoce)%clname = 'I_QnsOce'
      smysnd(jps_oemp  )%clname = 'IOEvaMPr'
      smysnd(jps_otx1  )%clname = 'I_OTaux1'
      smysnd(jps_oty1  )%clname = 'I_OTauy1'
      smysnd(jps_taum  )%clname = 'I_TauMod'
      smysnd(jps_fice2 )%clname = 'IIceFrc'
      !
      IF( .NOT. ln_cpl_oce ) smysnd(:)%laction = .FALSE.   ! force default definition in case of OceanModel <-> NANUQ coupling
      smysnd( (/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/) )%laction = .TRUE.


      IF( sn_loc_vct_tau == 'T' ) THEN
         ! ==> we send `utau` & `vtau` both defined at T-points !
         smysnd(jps_otx1)%clgrid = 'T'        ! oce components given at T-point
         smysnd(jps_oty1)%clgrid = 'T'
      ELSE
         ! ==> we send `utau` & `vtau` defined at U- & V-points, respectively!
         smysnd(jps_otx1)%clgrid = 'U'        ! oce components given at T-point
         smysnd(jps_oty1)%clgrid = 'V'
      ENDIF

      ! Change first letter to couple with ocean if already coupled with sea-ice
      ! this is nedeed as each variable name used in the namcouple must be unique:
      ! for example O_SSTSST sent by OCE to NANUQ and therefore S_SSTSST sent by NANUQ to the Ocean
      DO jn = 1, jpsnd
         IF( smysnd(jn)%clname(1:1) == "O" ) smysnd(jn)%clname = "S"//smysnd(jn)%clname(2:LEN(smysnd(jn)%clname))
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





      ! =================================== !
      !   define variables for the coupler  !
      ! =================================== !
      CALL cpl_vardef( midcpl )

      !IF(ln_usecplmask) THEN
      !   CALL iom_open( 'cplmask', inum )
      !   CALL iom_get( inum, jpdom_unknown, 'cplmask', xcplmask(Nis0:Nie0,Njs0:Nje0,1,1:nn_cplmodel), &
      !      &          kstart = (/ mig(Nis0,0),mjg(Njs0,0),1 /), kcount = (/ Ni_0,Nj_0,nn_cplmodel /) )
      !   CALL iom_close( inum )
      !   xcplmask(Nis0:Nie0,Njs0:Nje0,1,0) = 1. - SUM( xcplmask(Nis0:Nie0,Njs0:Nje0,1,1:nn_cplmodel), dim = 3 )
      !ELSE
      xcplmask(Nis0:Nie0,Njs0:Nje0,1,:) = 1.
      xcplmask(Nis0:Nie0,Njs0:Nje0,1,0) = 1. - SUM( xcplmask(Nis0:Nie0,Njs0:Nje0,1,1:nn_cplmodel), dim = 3 )
      !ENDIF
      !
   END SUBROUTINE oss_cpl_init








   SUBROUTINE oss_cpl_rcv( kt, k_foss, k_ice )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE oss_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in srcv(midcpl)%fld%z3 array). called at each time step.
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
      !! ** Action  :   update  utau, vtau   ocean stress
      !!                      (@ U- and V-points, respectively, if `sn_loc_vct_tau=='C'`)
      !!                      (@ T-points, if `sn_loc_vct_tau=='T'`)
      !!
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
      INTEGER, INTENT(in) ::   k_foss      ! frequency of sbc (-> ice model) computation
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
      !REAL(wp), DIMENSION(Nis0:Nie0,Njs0:Nje0) :: ztx, zty, zmsk, zemp, zqns, zqsr
      !!----------------------------------------------------------------------

      !                                                      ! ======================
      !                                                      ! Receive all the fields
      !                                                      ! ======================

      isec = ( kt - nit000 ) * NINT( rn_Dt )                 ! date of exchanges

      DO jn = 1, jprcv                                       ! receive fields sent by the ocean model
         IF( smyrcv(jn)%laction ) CALL cpl_rcv( midcpl, jn, isec, smyrcv(jn)%z3, nrcvinfo(jn), xcplmask(Nis0:Nie0,Njs0:Nje0,1:1,1:nn_cplmodel) )
      END DO

      !                                                      ! ================== !
      !    1                                                 !        SST         !
      !                                                      ! ================== !
      IF( smyrcv(jpr_sst)%laction ) THEN                     ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         sst_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_sst)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      !                                                      ! ================== !
      !    2                                                 !        SSS         !
      !                                                      ! ================== !
      IF( smyrcv(jpr_sss)%laction ) THEN                     ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         sss_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_sss)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      !                                                      ! ================== !
      !                                                      !        SSH         !
      !    3                                                 ! ================== !
      IF( smyrcv(jpr_ssh )%laction ) THEN                    ! received by sas in case of opa <-> sas coupling
         ssh_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_ssh )%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      !                                                      ! ================== !
      !                                                      !  surface currents  !
      !    4, 5                                              ! ================== !
      IF( smyrcv(jpr_ssu)%laction ) THEN                     ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         ssu_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_ssu)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      IF( smyrcv(jpr_ssv)%laction ) THEN
         ssv_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_ssv)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      !                                                      ! ======================== !
      !                                                      !  first T level thickness !
      !   6                                                  ! ======================== !
      IF( smyrcv(jpr_e3t)%laction ) THEN                     ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         e3t_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_e3t)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF
      !                                                      ! ================================ !
      !                                                      !  fraction of solar net radiation !
      !    7                                                 ! ================================ !
      IF( smyrcv(jpr_frq)%laction ) THEN                     ! received by NANUQ in case of OceanModel <-> NANUQ coupling
         frq_m(Nis0:Nie0,Njs0:Nje0) = smyrcv(jpr_frq)%z3(Nis0:Nie0,Njs0:Nje0,1)
      ENDIF

#if defined _TRDBG
      CALL test4nan('`ssh_m` <= `oss_cpl_rcv`', ssh_m(Nis0:Nie0,Njs0:Nje0), kkt=kt, lStop=.FALSE. )
      CALL test4nan('`sss_m` <= `oss_cpl_rcv`', sss_m(Nis0:Nie0,Njs0:Nje0), kkt=kt, lStop=.FALSE. )
      CALL test4nan('`sst_m` <= `oss_cpl_rcv`', sst_m(Nis0:Nie0,Njs0:Nje0), kkt=kt, lStop=.FALSE. )
      CALL test4nan('`ssu_m` <= `oss_cpl_rcv`', ssu_m(Nis0:Nie0,Njs0:Nje0), kkt=kt, lStop=.FALSE. )
      CALL test4nan('`ssv_m` <= `oss_cpl_rcv`', ssv_m(Nis0:Nie0,Njs0:Nje0), kkt=kt, lStop=.FALSE. )
#endif

      CALL lbc_lnk( 'osscpl', sst_m,'T', 1._wp, sss_m,'T', 1._wp, ssh_m,'T',1._wp, e3t_m,'T',1._wp, frq_m,'T',1._wp, &
         &                    ssu_m,'U',-1._wp, ssv_m,'V',-1._wp,   ldfull = .TRUE. )

      !! GPU ==> they are all updated onto device later on in: `oss()@ossmod.F90` !

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
      IF( smysnd(jps_otx1  )%laction ) THEN
         !
#if defined key_verbose
         IF(lwp) THEN
            IF( sn_loc_vct_tau == 'T' ) THEN
               PRINT *, ' *** LOLO: `oss_cpl_snd` sending `utau` @ T-point to ocean model, kt =',kt
            ELSE
               PRINT *, ' *** LOLO: `oss_cpl_snd` sending `utau` @ U-point to ocean model, kt =',kt
            ENDIF
         ENDIF
#endif
         !
         !$acc update self ( utau )
         IF( ln_cpl_oce_croco ) THEN
            CALL cpl_snd( midcpl, jps_otx1 , isec, RESHAPE ( utau(Nis0-1:Nie0-1,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 1
         ELSE
            !! Ocean component below is on a normal C-grid:
            CALL cpl_snd( midcpl, jps_otx1 , isec, RESHAPE ( utau(Nis0:Nie0,Njs0:Nje0),     (/Ni_0,Nj_0,1/) ), info ) ! 1
         ENDIF
         !
         IF( (sn_loc_vct_tau=='T').AND.(iom_use('utau_oa3_t')) )  CALL iom_put( 'utau_oa3_t', utau )
         IF( (sn_loc_vct_tau=='C').AND.(iom_use('utau_oa3_u')) )  CALL iom_put( 'utau_oa3_u', utau )
         !
      ENDIF

      ! j-component of surface stress (aka flux of momentum):
      IF( smysnd(jps_oty1  )%laction ) THEN
         !
#if defined key_verbose
         IF(lwp) THEN
            IF( sn_loc_vct_tau == 'T' ) THEN
               PRINT *, ' *** LOLO: `oss_cpl_snd` sending `vtau` @ T-point to ocean model, kt =',kt
            ELSE
               PRINT *, ' *** LOLO: `oss_cpl_snd` sending `vtau` @ V-point to ocean model, kt =',kt
            ENDIF
         ENDIF
#endif
         !
         !$acc update self ( vtau )
         IF( ln_cpl_oce_croco ) THEN
            CALL cpl_snd( midcpl, jps_oty1 , isec, RESHAPE ( vtau(Nis0:Nie0,Njs0-1:Nje0-1), (/Ni_0,Nj_0,1/) ), info ) ! 2
         ELSE
            !! Ocean component below is on a normal C-grid:
            CALL cpl_snd( midcpl, jps_oty1 , isec, RESHAPE ( vtau(Nis0:Nie0,Njs0:Nje0),     (/Ni_0,Nj_0,1/) ), info ) ! 2
         ENDIF
         !
         IF( (sn_loc_vct_tau=='T').AND.(iom_use('vtau_oa3_t')) )  CALL iom_put( 'vtau_oa3_t', vtau )
         IF( (sn_loc_vct_tau=='C').AND.(iom_use('vtau_oa3_v')) )  CALL iom_put( 'vtau_oa3_v', vtau )
         !
      ENDIF

      ! Non-solar heat flux:
      IF( smysnd(jps_qnsoce)%laction ) THEN
         !$acc update self ( qns(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_qnsoce, isec, RESHAPE ( qns(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 3
         !
         IF( iom_use('qns_oa3' ) )  CALL iom_put( "qns_oa3"    ,  qns  )
      ENDIF

      ! Solar heat flux:
      IF( smysnd(jps_qsroce)%laction ) THEN
         !$acc update self ( qsr(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_qsroce, isec, RESHAPE ( qsr(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 4
         !
         IF( iom_use('qsr_oa3' ) )  CALL iom_put( "qsr_oa3"    ,  qsr  )
      ENDIF

      ! Net freshwater flux (excluding continental runoffs):
      IF( smysnd(jps_oemp  )%laction ) THEN
         !$acc update self ( emp(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_oemp  , isec, RESHAPE ( emp(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 5
         !
         IF( iom_use('emp_oa3' ) )  CALL iom_put( "emp_oa3"    ,  emp  )
      ENDIF

      ! Salt flux:
      IF( smysnd(jps_sflx  )%laction ) THEN
         !$acc update self ( sfx(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_sflx  , isec, RESHAPE ( sfx(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 6
         !
         IF( iom_use('sfx_oa3' ) )  CALL iom_put( "sfx_oa3",  sfx  )
      ENDIF

      ! Modulus of surface stress:
      IF( smysnd(jps_taum  )%laction ) THEN
         !$acc update self ( taum(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_taum  , isec, RESHAPE ( taum(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 7
         !
         IF( iom_use('taum_oa3' ) )  CALL iom_put( "taum_oa3",  taum  )
      ENDIF

      ! Sea-ice concentration:
      IF( smysnd(jps_fice2 )%laction ) THEN
         !$acc update self ( at_i(Nis0:Nie0,Njs0:Nje0) )
         CALL cpl_snd( midcpl, jps_fice2,  isec, RESHAPE ( at_i(Nis0:Nie0,Njs0:Nje0), (/Ni_0,Nj_0,1/) ), info ) ! 8
         !
         IF( iom_use('fr_i_oa3' ) )  CALL iom_put( "fr_i_oa3", at_i )
      ENDIF

   END SUBROUTINE oss_cpl_snd

   !!======================================================================
END MODULE osscpl

#undef smyrcv
#undef smysnd
