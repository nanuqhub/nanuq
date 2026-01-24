MODULE ossprs
   !!======================================================================
   !!                       ***  MODULE  ossprs  ***
   !!
   !!     Prescribed ocean surface state for standalone simulations
   !!
   !!======================================================================
   !!
   !!----------------------------------------------------------------------
   !!   oss_prs_init  : initialization, namelist read, and SAVEs control
   !!   oss_prs_rcv       : Interpolation of the fields
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain: variables
   USE oss_nnq        ! surface module: variables
   USE phycst         ! physical constants
   USE eosbn2 , ONLY : t_eos10_fzp_scl, eos10_fzp_2d
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp        ! distributed memory computing library
   USE prtctl         ! print control
   USE fldread        ! read input fields

   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oss_prs_init   ! called by sbc_init
   PUBLIC   oss_prs_rcv    ! called by sbc
   PUBLIC   oss_prs_slab   ! called by ice_stp

   LOGICAL, PUBLIC    ::   ln_ssxread    ! Ice initialisation: =T read a file ; =F anaytical initilaistion
   LOGICAL            ::   ln_3d_uve     ! specify whether input velocity data is 3D
   LOGICAL,  PUBLIC   ::   ln_read_e3t   ! specify whether we must read `e3t` or not
   LOGICAL,  PUBLIC   ::   ln_read_frq   ! specify whether we must read `frq` or not
   REAL(wp), PUBLIC   ::   rn_e3t_0, rn_frq_0
   REAL(wp), PUBLIC   ::   rn_mld_0      ! global value of mixed layer depth to fall back on when `ln_slab_sst=T` & `sn_mld='NOT USED'`
   LOGICAL,  PUBLIC   ::   ln_ssv_T      ! prescribed sea surface velocities are provided on T grid points (default is U,V!)
   LOGICAL,  PUBLIC   ::   ln_ssv_Fgrid  ! also read SSU @ V-points & SSV @ U-points in prescribed ocean surface state (IF ln_ssxread )
   LOGICAL,  PUBLIC   ::   ln_slab_sst   ! For standalone mode only: if ln_slab_sst=T => correct the prescribed SST based on a slab model approach
   REAL(wp), PUBLIC   ::   rn_ncs_sst    ! nudging coefficient for SLAB bulk SST correction

   !$acc declare create( rn_e3t_0, rn_frq_0, rn_mld_0, ln_ssv_T, ln_ssv_Fgrid, rn_ncs_sst )

   CHARACTER(len=100) ::   cn_dir        ! Root directory for location of ssm files

   LOGICAL     ::   l_initdone = .false.
   INTEGER     ::   nfld_3d
   INTEGER     ::   nfld_2d

   INTEGER     ::   jf_tem         ! index of temperature
   INTEGER     ::   jf_sal         ! index of salinity
   INTEGER     ::   jf_usp         ! index of u velocity component
   INTEGER     ::   jf_vsp         ! index of v velocity component
   INTEGER     ::   jf_ssh         ! index of sea surface height
   INTEGER     ::   jf_mld         ! index of mixed layer depth
   INTEGER     ::   jf_e3t         ! index of first T level thickness
   INTEGER     ::   jf_frq         ! index of fraction of qsr absorbed in the 1st T level

   LOGICAL     ::   lk_read_mld

   ! Following is only for idealized test-cases for 2D advection featuring advection of scalars @ F-points:
   INTEGER     ::   jf_usf         ! index of u velocity component @ V-points
   INTEGER     ::   jf_vsf         ! index of v velocity component @ U-points

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_3d  ! structure of input fields (file information, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_2d  ! structure of input fields (file information, fields read)

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE oss_prs_rcv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE oss_prs_rcv  ***
      !!
      !! ** Purpose :  Prepares dynamics and physics fields from a NEMO run
      !!               for an off-line simulation using surface processes only
      !!
      !! ** Method : calculates the position of data
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'oss_prs_rcv')

      IF ( ln_ssxread ) THEN
         IF( nfld_3d > 0 ) CALL fld_read( kt, 1, sf_ssm_3d )      !==   read data at kt time step   ==!
         IF( nfld_2d > 0 ) CALL fld_read( kt, 1, sf_ssm_2d )      !==   read data at kt time step   ==!
         !
         IF( TRIM(sf_ssm_2d(jf_usp)%clrootname)=='NOT USED' )  sf_ssm_2d(jf_usp)%fnow(:,:,1) = 0._wp
         IF( TRIM(sf_ssm_2d(jf_vsp)%clrootname)=='NOT USED' )  sf_ssm_2d(jf_vsp)%fnow(:,:,1) = 0._wp
         ssu_m(:,:) = sf_ssm_2d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
         ssv_m(:,:) = sf_ssm_2d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity
         !
         IF( TRIM(sf_ssm_2d(jf_sal)%clrootname)=='NOT USED' )                 sf_ssm_2d(jf_sal)%fnow(:,:,1) = 35._wp
         IF( TRIM(sf_ssm_2d(jf_tem)%clrootname)=='NOT USED' ) CALL eos10_fzp_2d(sf_ssm_2d(jf_sal)%fnow(:,:,1), sf_ssm_2d(jf_tem)%fnow(:,:,1))
         IF( TRIM(sf_ssm_2d(jf_ssh)%clrootname)=='NOT USED' )                 sf_ssm_2d(jf_ssh)%fnow(:,:,1) = 0._wp
         sst_m(:,:) = sf_ssm_2d(jf_tem)%fnow(:,:,1) * xmskt(:,:)    ! temperature
         sss_m(:,:) = sf_ssm_2d(jf_sal)%fnow(:,:,1) * xmskt(:,:)    ! salinity
         ssh_m(:,:) = sf_ssm_2d(jf_ssh)%fnow(:,:,1) * xmskt(:,:)    ! sea surface height
         !
         IF(ln_slab_sst) THEN
            IF( TRIM(sf_ssm_2d(jf_mld)%clrootname)=='NOT USED' )   sf_ssm_2d(jf_mld)%fnow(:,:,1) = rn_mld_0
            IF( lk_read_mld ) THEN
               mld_m(:,:) = sf_ssm_2d(jf_mld)%fnow(:,:,1) * xmskt(:,:)    ! mixed layer depth
            ELSE
               mld_m(:,:) = rn_mld_0
            ENDIF
         ENDIF
         !
         IF( ln_read_e3t ) THEN
            e3t_m(:,:) = sf_ssm_2d(jf_e3t)%fnow(:,:,1)
         ELSE
            e3t_m(:,:) = rn_e3t_0
         ENDIF
         IF( ln_read_frq ) THEN
            frq_m(:,:) = sf_ssm_2d(jf_frq)%fnow(:,:,1) * xmskt(:,:) ! solar penetration
         ELSE
            frq_m(:,:) = rn_frq_0
         ENDIF
         !
         IF( ln_ssv_Fgrid ) THEN
            IF( TRIM(sf_ssm_2d(jf_usf)%clrootname)=='NOT USED' )  sf_ssm_2d(jf_usf)%fnow(:,:,1) = 0._wp
            IF( TRIM(sf_ssm_2d(jf_vsf)%clrootname)=='NOT USED' )  sf_ssm_2d(jf_vsf)%fnow(:,:,1) = 0._wp
            ssu_v_m(:,:) = sf_ssm_2d(jf_usf)%fnow(:,:,1) * vmask(:,:,1)    ! u-velocity
            ssv_u_m(:,:) = sf_ssm_2d(jf_vsf)%fnow(:,:,1) * umask(:,:,1)    ! v-velocity
         ENDIF
         !
      ELSE
         ! When forcing not read in netCDF file(s), fall back on the following:
         sss_m(:,:) = 35._wp                             ! =35. to obtain a physical value for the freezing point
         CALL eos10_fzp_2d( sss_m(:,:), sst_m(:,:) )     ! sst_m is set at the freezing point
         ssu_m(:,:) = 0._wp
         ssv_m(:,:) = 0._wp
         ssh_m(:,:) = 0._wp
         IF( ln_slab_sst ) mld_m(:,:) = rn_mld_0
         e3t_m(:,:) = rn_e3t_0
         frq_m(:,:) = rn_frq_0
         IF( ln_ssv_Fgrid ) THEN
            ssu_v_m(:,:) = 0._wp
            ssv_u_m(:,:) = 0._wp
         ENDIF
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN            ! print control
         CALL prt_ctl(tab2d_1=sst_m, clinfo1=' sst_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=sss_m, clinfo1=' sss_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=ssu_m, clinfo1=' ssu_m   - : ', mask1=umask   )
         CALL prt_ctl(tab2d_1=ssv_m, clinfo1=' ssv_m   - : ', mask1=vmask   )
         IF( lk_read_mld ) CALL prt_ctl(tab2d_1=mld_m, clinfo1=' mld_m   - : ', mask1=tmask   )
         IF( ln_read_frq ) CALL prt_ctl(tab2d_1=frq_m, clinfo1=' frq_m   - : ', mask1=tmask   )
      ENDIF

      IF( ln_timing )   CALL timing_stop( 'oss_prs_rcv')
      !
   END SUBROUTINE oss_prs_rcv







   SUBROUTINE oss_prs_slab( kt, pdt, pA, psst_m, psss_m, pmld_m, pqsr, pqns, pemp, psst_s, psss_s, pt_bo )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE oss_prs_slab  ***
      !!
      !! ** Purpose : correct the prescribed SST & SSS following a simplistic
      !!              slab ocean.
      !!
      !!
      !! ** Method : uses the net heat & freshwater flux into the liquid ocean
      !!             and the prescribed MLD...
      !!             with a nudging approach not to deviate to much from
      !!             prescribed SST & SSS...
      !!
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in)    ::   kt   ! ocean time-step
      REAL(wp),                     INTENT(in)    :: pdt  ! time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pA   ! sea-ice concentration      
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: psst_m, psss_m, pmld_m
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pqsr   ! solar heat flux to (with liquid ocean albedo decrease considered)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: pqns, pemp
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: psst_s, psss_s, pt_bo
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      LOGICAL  ::   ldebg
      REAL(wp) :: zFWkg, zQjoules, zinc, zsss_p, zsst_p, zsss_n, zsst_n, znc, z1_mld, zdum
      !REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      !REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start( 'oss_prs_slab')
      !$acc data present( psst_m, psss_m, pmld_m, pqns, pemp, psst_s, psss_s, pt_bo )

      ! ==> we use a prescribed SST read in netCDF files
      !     and we are correcting this SST in a simple SLAB ocean fashion:
      !
      ! A/ temperature increment for the whole mixed layer based on the flux received by the liquid
      !    ocean during the previous time step
      
      IF( ln_rstart ) THEN
         PRINT *, 'FIX `sst_s` & `sss_s` for restarts!!!'; STOP
      ENDIF

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            !ldebg = (narea==3).AND.(ji==17).AND.(jj==44)
            !ldebg = .true.

            zsss_p = psss_m(ji,jj) ! prescribed/observed SSS as read in netCDF file
            zsst_p = psst_m(ji,jj) ! prescribed/observed SST as read in netCDF file and potentially corrected / `eos_fzp`...

            !IF( ldebg ) THEN
            !   PRINT *, '---'
            !   !PRINT *, '  * glamt   = ', REAL(glamt(ji,jj),4)
            !   !PRINT *, '  * gphit   = ', REAL(gphit(ji,jj),4)
            !   !PRINT *, '  * qns   = ', REAL(qns(ji,jj),4)
            !   PRINT *, ' *-emp_b        = ', REAL(-pemp(ji,jj),4)
            !   PRINT *, ' * qns_b        = ', REAL(pqns(ji,jj),4)
            !   PRINT *, ' * mld_m        = ', REAL(pmld_m(ji,jj),4)
            !   PRINT *, ' * sss_m        = ', REAL(psss_m(ji,jj),4)
            !   PRINT *, ' * FPT ==> t_bo = ', REAL( pt_bo(ji,jj),4)
            !   PRINT *, ' * sst_m        = ', REAL(psst_m(ji,jj),4)
            !   PRINT *, ' * zsst_p       = ', REAL(zsst_p,4)
            !   PRINT *, '---'
            !ENDIF

            z1_mld  = 1._wp / MAX( pmld_m(ji,jj), 0.1_wp )

            zFWkg    = -pemp(ji,jj) * pdt            ! Mass of freshwater received/lost per surface area during `pdt` (kg/m2)
            zinc =  zFWkg/rhow * z1_mld * psss_m(ji,jj) !  zFWkg/rhow ~ volume / m2 = m ; zFWkg/rhow * z1_mld ~ '-'
            zsss_n  = psss_s(ji,jj) + zinc              !  expected new salinity in the MLD...
            !
            !IF( ldebg ) THEN
            !   PRINT *, '---'
            !   PRINT *, '  ==> zFWkg = ', REAL(zFWkg,4),'kg/m^2, during',REAL(pdt/3600._wp,4),'hours'
            !   PRINT *, '            ===> equivalent to', REAL(zFWkg/rhow*1000._wp,4),'mm of water'
            !   PRINT *, '  ==> zds_inc = ', REAL(zinc,4),'psu'
            !   PRINT *, '  * zsss_n = ', REAL(zsss_n,4)
            !ENDIF
            
            !IF((ji==40).AND.(jj==40)) THEN
            !   PRINT *, ''
            !   PRINT *, 'LOLO: solar heat flux available and actually counted for SLAB warming with A=',REAL(pA(ji,jj),4)
            !   PRINT *, 'LOLO: qsr, qsr_counted:', REAL(pqsr(ji,jj),4), REAL((1._wp -pA(ji,jj))*pqsr(ji,jj),4)
            !   PRINT *, ''
            !ENDIF

            
            zQjoules = ( pqns(ji,jj) +  (1._wp -pA(ji,jj))*pqsr(ji,jj) ) * pdt   ! Energy received/lost per surface area during `pdt` (J/m2) !#LOLOfixme: consider solar penetration for qsr !?
            zinc =  zQjoules * r1_rho0_rcp * z1_mld   !  rho0_rcp ~ J/K/m3 => rho0_rcp*mld ~ J/K/m2 => zinc = (J/m2) / (J/K/m2) ==> K !
            zsst_n  = psst_s(ji,jj) + zinc                    !  expected new temperature in the MLD...
            !CALL eos_fzp_0d( zsss_n, zdum )
            zdum = t_eos10_fzp_scl( zsss_n )
            zsst_n = MAX( zsst_n , zdum )  ! slab SST cannot be colder than slab freezing-point temperature.


            ! We want the nudging coefficient to be `1` where open ocean and `rn_ncs_sst` where ice fraction = 1:
            !znc = (1._wp - at_i(ji,jj)) + at_i(ji,jj)*rn_ncs_sst
            znc = rn_ncs_sst

            !IF( ldebg ) THEN
            !   PRINT *, '---'
            !   PRINT *, '  ==> zQjoules = ', REAL(zQjoules,4),'J/m^2, during',REAL(pdt/3600._wp,4),'hours'
            !   PRINT *, '  ==> zdt_inc = ', REAL(zinc,4),'K'
            !   PRINT *, '  * zsst_n = ', REAL(zsst_n,4)
            !   PRINT *, ''
            !   PRINT *, '  * NDG coeff, znc = ', REAL(znc,4)
            !ENDIF

            psss_s(ji,jj) = zsss_n - znc * ( zsss_n - zsss_p ) ! nudging towards the prescribed SSS
            psst_s(ji,jj) = zsst_n - znc * ( zsst_n - zsst_p ) ! nudging towards the prescribed SST
            
            !CALL eos_fzp_0d( psss_s(ji,jj), pt_bo(ji,jj) )
            pt_bo(ji,jj) = t_eos10_fzp_scl( psss_s(ji,jj) )
            psst_s(ji,jj) = MAX( psst_s(ji,jj) , pt_bo(ji,jj) )  ! slab SST cannot be colder than slab freezing-point temperature.


            ! DEBUG:
            !CALL eos_fzp_0d( psss_s(ji,jj), zsst_n )
            !IF(psst_s(ji,jj)<zsst_n) THEN
            !   PRINT *, 'LOLO FUCKUP: `sst_s(ji,jj)<zsst_n`!!! icestp.F90'
            !   STOP
            !ENDIF
            !DEBUG.

            !IF( ldebg ) THEN
            !   PRINT *, ' *  NEW t_bo = ', REAL( pt_bo(ji,jj),4)
            !   PRINT *, '  * sst_s - sst_m =', REAL(psst_s(ji,jj) - psst_m(ji,jj),4)
            !   PRINT *, '  * sss_s - sss_m =', REAL(psss_s(ji,jj) - psss_m(ji,jj),4)
            !   PRINT *, '---'
            !   PRINT *, ''
            !ENDIF

         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop( 'oss_prs_slab')
      !
   END SUBROUTINE oss_prs_slab






   SUBROUTINE oss_prs_init( )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE oss_prs_init  ***
      !!
      !! ** Purpose :   Initialisation of sea surface mean data
      !!----------------------------------------------------------------------
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: inum, idv, idimv, jpm, kmld, ke3t, kfrq  ! local integer
      INTEGER  :: ios                                ! Local integer output status for namelist read
      !!
      CHARACTER(len=100)                     ::  cn_dir       ! Root directory for location of core files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_3d       ! array of namelist information on the fields to read
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_2d       ! array of namelist information on the fields to read
      TYPE(FLD_N) ::   sn_tem, sn_sal                     ! information about the fields to be read
      TYPE(FLD_N) ::   sn_usp, sn_vsp
      TYPE(FLD_N) ::   sn_usf, sn_vsf
      TYPE(FLD_N) ::   sn_ssh, sn_mld
      TYPE(FLD_N) ::   sn_e3t, sn_frq
      !!
      TYPE(FLD_N) ::   sn_ifr, sn_tic, sn_ial
      !!
      NAMELIST/namoss_ssx/ ln_ssxread, ln_3d_uve, ln_read_e3t, rn_e3t_0, ln_read_frq, rn_frq_0, &
         &                 ln_ssv_T, ln_ssv_Fgrid, ln_slab_sst, rn_mld_0, rn_ncs_sst,           &
         &                 cn_dir, sn_usp, sn_vsp, sn_tem, sn_sal, sn_ssh, sn_mld,              &
         &                 sn_e3t, sn_frq, sn_usf, sn_vsf, sn_ifr, sn_tic, sn_ial
      !!----------------------------------------------------------------------
      !
      IF( ln_rstart .AND. ln_cpl_oce ) RETURN
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'oss_prs_init : sea surface mean data initialisation '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      READ_NML_REF(numnam,namoss_ssx)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namoss_ssx in reference namelist' )
      READ_NML_CFG(numnam,namoss_ssx)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namoss_ssx in configuration namelist' )
      IF(lwm) WRITE ( numond, namoss_ssx )
      !
      IF(lwp) THEN                              ! Control print
         WRITE(numout,*) '   Namelist namoss_ssx'
         WRITE(numout,*) '      Initialisation using an input file                                 ln_ssxread   = ', ln_ssxread
         WRITE(numout,*) '      Are we supplying a 3D u,v and e3 field                             ln_3d_uve    = ', ln_3d_uve
         WRITE(numout,*) '      Are we reading e3t (depth of 1st ocean level                   )   ln_read_e3t  = ', ln_read_e3t
         WRITE(numout,*) '                => prescribed value of `e3t` to fall back on                rn_e3t_0  = ', rn_e3t_0
         WRITE(numout,*) '      Are we reading frq (fraction of qsr absorbed in the 1st T level)   ln_read_frq  = ', ln_read_frq
         WRITE(numout,*) '                => prescribed value of `frq` to fall back on                rn_frq_0  = ', rn_frq_0
         WRITE(numout,*) '      Prescribed surface velocities provided @T rather than @U,V  points    ln_ssv_T  = ', ln_ssv_T
         WRITE(numout,*) '      Read SSU@V & SSV@U in prescribed OSS                               ln_ssv_Fgrid = ', ln_ssv_Fgrid
         WRITE(numout,*) '      Correction of the prescibed SST based on a "slab model" approach    ln_slab_sst = ', ln_slab_sst
         WRITE(numout,*) '                => prescribed value of `mld` to fall back on                rn_mld_0  = ', rn_mld_0
         WRITE(numout,*) '                => nudging coefficient for SLAB bulk SST correction       rn_ncs_sst  = ', rn_ncs_sst

         WRITE(numout,*) '                         cn_dir =', cn_dir
         IF( ln_ssxread ) THEN
            WRITE(numout,*) '           *      sn_usp =', sn_usp
            WRITE(numout,*) '           *      sn_vsp =', sn_vsp
            WRITE(numout,*) '           *      sn_tem =', sn_tem
            WRITE(numout,*) '           *      sn_sal =', sn_sal
            WRITE(numout,*) '           *      sn_ssh =', sn_ssh
            IF(ln_slab_sst) WRITE(numout,*) '           *      sn_mld =', sn_mld
            IF(ln_read_e3t) WRITE(numout,*) '           *      sn_e3t =', sn_e3t
            IF(ln_read_frq) WRITE(numout,*) '           *      sn_frq =', sn_frq

            IF( ln_ssv_Fgrid ) THEN
               WRITE(numout,*) '           *      sn_usf =', sn_usf
               WRITE(numout,*) '           *      sn_vsf =', sn_vsf
            ENDIF
         ENDIF
      ENDIF

      lk_read_mld = ( ln_slab_sst .AND. (TRIM(sn_mld%clname)/='NOT USED') )
      IF(lwp) WRITE(numout,*) '           *      lk_read_mld =', lk_read_mld

      !! switch off stuff that isn't sensible with a standalone module
      !! note that we need oss_prs_rcv called first in sbc
      !
      IF( ln_ssxread ) THEN                       ! store namelist information in an array
         !

         !! Allocating non-mandatory arrays:
         IF( ln_slab_sst ) THEN
            IF(lwp) WRITE(numout,*) '  *** Allocating `mld_m` array !!!'
            ALLOCATE( mld_m(jpi,jpj) ,  STAT=ierr )
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate `mld_m`' )   ;   RETURN
            ENDIF
# if defined _OPENACC
            PRINT *, ' * info GPU: oss_prs_init() => adding `mld_m` array to memory!'
            !$acc enter data copyin( mld_m )
# endif
         ENDIF
         IF( ln_ssv_Fgrid ) THEN
            IF(lwp) WRITE(numout,*) '  *** Allocating `ssu_v_m` & `ssv_u_m` arrays !!!'
            ALLOCATE( ssu_v_m(jpi,jpj) , ssv_u_m(jpi,jpj),  STAT=ierr )
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate `ssu_v_m` & `ssv_u_m`' )   ;   RETURN
            ENDIF
# if defined _OPENACC
            PRINT *, ' * info GPU: oss_prs_init() => adding `ssu_v_m` & `ssv_u_m` arrays to memory!'
            !$acc enter data copyin( ssu_v_m, ssv_u_m )
# endif
         ENDIF


         !! following code is a bit messy, but distinguishes between when u,v are 3d arrays and
         !! when we have other 3d arrays that we need to read in
         !! so if a new field is added i.e. jf_new, just give it the next integer in sequence
         !! for the corresponding dimension (currently if ln_3d_uve is true, 4 for 2d and 3 for 3d,
         !! alternatively if ln_3d_uve is false, 6 for 2d and 1 for 3d), reset nfld_3d, nfld_2d,
         !! and the rest of the logic should still work
         !
         kmld = COUNT( (/lk_read_mld/) )
         ke3t = COUNT( (/ln_read_e3t/) )
         kfrq = COUNT( (/ln_read_frq/) )

         jf_tem = 1 ; jf_sal = 2 ; jf_ssh = 3   ! default 2D fields index
         !
         IF( ln_3d_uve ) THEN
            jf_usp = 1   ;   jf_vsp = 2   ;  jf_mld = 3  ; jf_e3t = 4  ; jf_frq = 5     ! define 3D fields index
            nfld_3d  = 2 + 0                                 ! number of 3D fields to read
            nfld_2d  = 3 + ke3t + kfrq                       ! number of 2D fields to read
         ELSE
            jf_usp = 4
            jf_vsp = 5
            jf_mld = jf_vsp + kmld
            jf_e3t = jf_mld + ke3t
            jf_frq = jf_e3t + kfrq
            jf_usf = jf_frq+1
            jf_vsf = jf_frq+2
            !
            nfld_3d  = 0                                     ! no 3D fields to read
            nfld_2d  = 3 + 2 + kmld + ke3t + kfrq + 2*COUNT( (/ln_ssv_Fgrid/) )     ! number of 2D fields to read
            IF(.NOT. lk_read_mld)  jf_mld=-1  ! safety! => will cause error
            IF(.NOT. ln_read_e3t)  jf_e3t=-1  ! safety! => will cause error
            IF(.NOT. ln_read_frq)  jf_frq=-1  ! safety! => will cause error
         ENDIF
         !
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( slf_3d(nfld_3d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate slf 3d structure' )   ;   RETURN
            ENDIF
            slf_3d(jf_usp) = sn_usp
            slf_3d(jf_vsp) = sn_vsp
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( slf_2d(nfld_2d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate slf 2d structure' )   ;   RETURN
            ENDIF
            slf_2d(jf_tem) = sn_tem ; slf_2d(jf_sal) = sn_sal ; slf_2d(jf_ssh) = sn_ssh

            IF( lk_read_mld )   slf_2d(jf_mld) = sn_mld
            IF( ln_read_e3t )   slf_2d(jf_e3t) = sn_e3t
            IF( ln_read_frq )   slf_2d(jf_frq) = sn_frq
            IF( .NOT. ln_3d_uve ) THEN
               slf_2d(jf_usp) = sn_usp ; slf_2d(jf_vsp) = sn_vsp
            ENDIF
            IF( ln_ssv_Fgrid ) THEN
               slf_2d(jf_usf) = sn_usf ; slf_2d(jf_vsf) = sn_vsf
            ENDIF
         ENDIF
         !
         ierr1 = 0    ! default definition if slf_?d(ifpr)%ln_tint = .false.
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( sf_ssm_3d(nfld_3d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate sf structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_3d
               ALLOCATE( sf_ssm_3d(ifpr)%fnow(jpi,jpj,jpk)    , STAT=ierr0 )
               IF( slf_3d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_3d(ifpr)%fdta(jpi,jpj,jpk,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'oss_prs_init : unable to allocate sf_ssm_3d array structure' )   ;   RETURN
               ENDIF
            END DO
            !                                         ! fill sf with slf_i and control print
            CALL fld_fill( sf_ssm_3d, slf_3d, cn_dir, 'oss_prs_init', '3D Data in file', 'namoss_ssx' )
            sf_ssm_3d(jf_usp)%cltype = 'U'   ;   sf_ssm_3d(jf_usp)%zsgn = -1._wp
            sf_ssm_3d(jf_vsp)%cltype = 'V'   ;   sf_ssm_3d(jf_vsp)%zsgn = -1._wp
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( sf_ssm_2d(nfld_2d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'oss_prs_init: unable to allocate sf 2d structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_2d
               ALLOCATE( sf_ssm_2d(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
               IF( slf_2d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_2d(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'oss_prs_init : unable to allocate sf_ssm_2d array structure' )   ;   RETURN
               ENDIF
            END DO
            !
            CALL fld_fill( sf_ssm_2d, slf_2d, cn_dir, 'oss_prs_init', '2D Data in file', 'namoss_ssx' )
            IF( .NOT. ln_3d_uve ) THEN
               IF( ln_ssv_T ) THEN
                  sf_ssm_2d(jf_usp)%cltype = 'T'   ;   sf_ssm_2d(jf_usp)%zsgn = -1._wp
                  sf_ssm_2d(jf_vsp)%cltype = 'T'   ;   sf_ssm_2d(jf_vsp)%zsgn = -1._wp
                  IF( ln_ssv_Fgrid ) THEN
                     CALL ctl_stop( 'oss_prs_init: `ln_ssv_Fgrid=T` prohibited when `ln_ssv_T=T` !!!' )   ;   RETURN
                  ENDIF
               ELSE
                  sf_ssm_2d(jf_usp)%cltype = 'U'   ;   sf_ssm_2d(jf_usp)%zsgn = -1._wp
                  sf_ssm_2d(jf_vsp)%cltype = 'V'   ;   sf_ssm_2d(jf_vsp)%zsgn = -1._wp
                  IF( ln_ssv_Fgrid ) THEN
                     sf_ssm_2d(jf_usf)%cltype = 'V'   ;   sf_ssm_2d(jf_usf)%zsgn = -1._wp
                     sf_ssm_2d(jf_vsf)%cltype = 'U'   ;   sf_ssm_2d(jf_vsf)%zsgn = -1._wp
                  ENDIF
               ENDIF !IF( ln_ssv_T )
            ENDIF !IF( .NOT. ln_3d_uve )
         ENDIF !IF( nfld_2d > 0 )
         !
         IF( nfld_3d > 0 )   DEALLOCATE( slf_3d, STAT=ierr )
         IF( nfld_2d > 0 )   DEALLOCATE( slf_2d, STAT=ierr )
         !
      ENDIF !IF( ln_ssxread )
      !
      CALL oss_prs_rcv( nit000 )   ! need to define ssx_m arrays used in iceistate
      l_initdone = .TRUE.
      !
      !$acc update device ( rn_e3t_0, rn_frq_0, rn_mld_0, ln_ssv_T, ln_ssv_Fgrid, rn_ncs_sst )

   END SUBROUTINE oss_prs_init

   !!======================================================================
END MODULE ossprs
