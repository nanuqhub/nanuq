MODULE sbcblk
   !!======================================================================
   !!                       ***  MODULE  sbcblk  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!                         Aerodynamic Bulk Formulas
   !!                        SUCCESSOR OF "sbcblk_core"
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original CORE code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier)  improved CORE bulk and its user interface
   !!            3.0  !  2006-06  (G. Madec)  sbc rewritting
   !!             -   !  2006-12  (L. Brodeau)  Original code for turb_core
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!            3.4  !  2011-11  (C. Harris)  Fill arrays required by CICE
   !!            3.7  !  2014-06  (L. Brodeau)  simplification and optimization of CORE bulk
   !!            4.0  !  2016-06  (L. Brodeau)  sbcblk_core becomes sbcblk and is not restricted to the CORE algorithm anymore
   !!                 !                        ==> based on AeroBulk (https://github.com/brodeau/aerobulk/)
   !!            4.0  !  2016-10  (G. Madec)  introduce a sbc_blk_init routine
   !!            4.0  !  2016-10  (M. Vancoppenolle)  Introduce conduction flux emulator (M. Vancoppenolle)
   !!            4.0  !  2019-03  (F. Lemarié & G. Samson)  add ABL compatibility (ln_abl=TRUE)
   !!            4.2  !  2020-12  (L. Brodeau) Introduction of various air-ice bulk parameterizations + improvements
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_init  : initialisation of the chosen bulk formulation as ocean surface boundary condition
   !!   sbc_blk       : bulk formulation as ocean surface boundary condition
   !!   blk_oce_1     : computes pieces of momentum, heat and freshwater fluxes over ocean for ABL model  (ln_abl=TRUE)
   !!   blk_oce_2     : finalizes momentum, heat and freshwater fluxes computation over ocean after the ABL step  (ln_abl=TRUE)
   !!             sea-ice case only :
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn   : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE ossprs  , ONLY : ln_slab_sst ! #DEBUG
   USE oss_nnq , ONLY : ssh_m, ssh_m, ssu_m, ssv_m, ssst, sssq, sst_s
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE lib_fortran    ! to use key_nosignedzero and glob_2Dsum
   !
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE par_ice , ONLY : jpl, nn_qtrice, rcloud_fra
   USE ice     , ONLY : at_i, a_i_b, at_i_b, hfx_err_dif  !#LOLOfixme => WHY USE `a_i_b` & `at_i_b` and not `a_i` & `at_i` ????
   USE icevar         ! for CALL ice_var_snwblow
   !
   USE sbcblk_algo_ice_an05
   USE sbcblk_algo_ice_lu12
   USE sbcblk_algo_ice_lg15
   !
   USE sbcblk_algo_ncar     ! => turb_ncar     : NCAR - (formerly known as CORE, Large & Yeager, 2009)
   USE sbcblk_algo_coare3p0 ! => turb_coare3p0 : COAREv3.0 (Fairall et al. 2003)
   !USE sbcblk_algo_coare3p6 ! => turb_coare3p6 : COAREv3.6 (Fairall et al. 2018 + Edson et al. 2013)
   USE sbcblk_algo_ecmwf    ! => turb_ecmwf    : ECMWF (IFS cycle 45r1)
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE timing         ! Timing

   USE sbc_phy        ! Catalog of functions for physical/meteorological parameters in the marine boundary layer

   USE par_ice , ONLY : rn_snwblow

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_init  ! called in sbcmod
   PUBLIC   sbc_blk       ! called in sbcmod
   PUBLIC   blk_oce_1     ! called in sbcabl
   PUBLIC   blk_oce_2     ! called in sbcabl

   PUBLIC   blk_ice_1     ! routine called in icesbc
   PUBLIC   blk_ice_2     ! routine called in icesbc
   PUBLIC   blk_ice_qcn   ! routine called in icesbc

   INTEGER , PUBLIC, PARAMETER ::   jp_wndi  =  1   ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_wndj  =  2   ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_tair  =  3   ! index of 10m air temperature             (Kelvin)
   INTEGER , PUBLIC, PARAMETER ::   jp_humi  =  4   ! index of specific humidity               (kg/kg)
   INTEGER , PUBLIC, PARAMETER ::   jp_qsr   =  5   ! index of solar heat                      (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_qlw   =  6   ! index of Long wave                       (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_prec  =  7   ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_snow  =  8   ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_slp   =  9   ! index of sea level pressure              (Pa)
   INTEGER , PUBLIC, PARAMETER ::   jp_uoatm = 10   ! index of surface current (i-component)
   !                                                !          seen by the atmospheric forcing (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_voatm = 11   ! index of surface current (j-component)
   !                                                !          seen by the atmospheric forcing (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_cc    = 12   ! index of cloud cover                     (-)      range:0-1
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgi  = 13   ! index of ABL geostrophic wind or hpg (i-component) (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgj  = 14   ! index of ABL geostrophic wind or hpg (j-component) (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jpfld    = 14   ! maximum number of files to read

   ! Warning: keep this structure allocatable for Agrif...
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input atmospheric fields (file informations, fields read)

   !                           !!* Namelist namsbc_blk : bulk parameters
   LOGICAL  ::   ln_NCAR        ! "NCAR"      algorithm   (Large and Yeager 2008)
   LOGICAL  ::   ln_COARE_3p0   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   !LOGICAL  ::   ln_COARE_3p6   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   LOGICAL  ::   ln_ECMWF       ! "ECMWF"     algorithm   (IFS cycle 45r1)
   !LOGICAL  ::   ln_ANDREAS     ! "ANDREAS"   algorithm   (Andreas et al. 2015)
   !
   !#LB:
   LOGICAL  ::   ln_Cx_ice_cst             ! use constant air-ice bulk transfer coefficients (value given in namelist's rn_Cd_i, rn_Ce_i & rn_Ch_i)
   !LOGICAL  ::   ln_Cx_ice_EASY            ! air-ice bulk transfer coefficients based on Andreas et al., 2005
   REAL(wp) ::   rn_Cd_i, rn_Ce_i, rn_Ch_i ! values for  "    "
   LOGICAL  ::   ln_Cx_ice_AN05            ! air-ice bulk transfer coefficients based on Andreas et al., 2005
   LOGICAL  ::   ln_Cx_ice_LU12            ! air-ice bulk transfer coefficients based on Lupkes et al., 2012
   LOGICAL  ::   ln_Cx_ice_LG15            ! air-ice bulk transfer coefficients based on Lupkes & Gryanik, 2015
   !#LB.
   !
   LOGICAL  ::   ln_crt_fbk     ! Add surface current feedback to the wind stress computation  (Renault et al. 2020)
   REAL(wp) ::   rn_stau_a      ! Alpha and Beta coefficients of Renault et al. 2020, eq. 10: Stau = Alpha * Wnd + Beta
   REAL(wp) ::   rn_stau_b      !
   !
   REAL(wp)         ::   rn_zqt    ! z(q,t) : height of humidity and temperature measurements
   REAL(wp)         ::   rn_zu     ! z(u)   : height of wind measurements
   !$acc declare create( rn_zqt, rn_zu )
   !
   INTEGER          :: nn_iter_algo   !  Number of iterations in bulk param. algo ("stable ABL + weak wind" requires more)

   LOGICAL, PUBLIC  ::   ln_skin_cs     ! use the cool-skin (only available in ECMWF and COARE algorithms) !LB
   LOGICAL, PUBLIC  ::   ln_skin_wl     ! use the warm-layer parameterization (only available in ECMWF and COARE algorithms) !LB
   LOGICAL, PUBLIC  ::   ll_skin        ! `ll_skin = ln_skin_cs .OR. ln_skin_wl`
   LOGICAL  ::   ln_humi_sph    ! humidity read in files ("sn_humi") is specific humidity [kg/kg] if .true. !LB
   LOGICAL  ::   ln_humi_dpt    ! humidity read in files ("sn_humi") is dew-point temperature [K] if .true. !LB
   LOGICAL  ::   ln_humi_rlh    ! humidity read in files ("sn_humi") is relative humidity     [%] if .true. !LB
   LOGICAL  ::   ln_tair_pot    ! temperature read in files ("sn_tair") is already potential temperature (not absolute)
   !
   INTEGER  ::   nhumi          ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER :: np_humi_sph = 1
   INTEGER, PARAMETER :: np_humi_dpt = 2
   INTEGER, PARAMETER :: np_humi_rlh = 3

   INTEGER  ::   nblk           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_NCAR      = 1   ! "NCAR" algorithm        (Large and Yeager 2008)
   INTEGER, PARAMETER ::   np_COARE_3p0 = 2   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   !INTEGER, PARAMETER ::   np_COARE_3p6 = 3   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   INTEGER, PARAMETER ::   np_ECMWF     = 4   ! "ECMWF" algorithm       (IFS cycle 45r1)
   !#LB:
   ! Same, over sea-ice:
   INTEGER  ::   nblk_ice           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_ice_cst  = 1   ! constant transfer coefficients
   INTEGER, PARAMETER ::   np_ice_an05 = 3   ! Andreas et al., 2005
   INTEGER, PARAMETER ::   np_ice_lu12 = 4   ! Lupkes el al., 2012
   INTEGER, PARAMETER ::   np_ice_lg15 = 5   ! Lupkes & Gryanik, 2015
   !LB.



   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_blk_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_init  ***
      !!
      !! ** Purpose :   choose and initialize a bulk formulae formulation
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   jfpr                  ! dummy loop indice and argument
      INTEGER  ::   ios, ierror, ioptio   ! Local integer
      !!
      CHARACTER(len=100)            ::   cn_dir                ! Root directory for location of atmospheric forcing files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                 ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj , sn_humi, sn_qsr      ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair , sn_prec, sn_snow     !       "                        "
      TYPE(FLD_N) ::   sn_slp , sn_uoatm, sn_voatm             !       "                        "
      TYPE(FLD_N) ::   sn_cc, sn_hpgi, sn_hpgj                 !       "                        "
      INTEGER     ::   ipka                                    ! number of levels in the atmospheric variable
      NAMELIST/namsbc_blk/ ln_NCAR, ln_COARE_3p0, ln_ECMWF, &   ! bulk algorithm    | ln_COARE_3p6, 
         &                 rn_zqt, rn_zu, nn_iter_algo, ln_skin_cs, ln_skin_wl,       &
         &                 ln_crt_fbk, rn_stau_a, rn_stau_b,                          &   ! current feedback
         &                 ln_humi_sph, ln_humi_dpt, ln_humi_rlh, ln_tair_pot,        &
         &                 ln_Cx_ice_cst, rn_Cd_i, rn_Ce_i, rn_Ch_i,  &
         &                 ln_Cx_ice_AN05, ln_Cx_ice_LU12, ln_Cx_ice_LG15,            &
         &                 cn_dir,                                                    &
         &                 sn_wndi, sn_wndj, sn_qsr, sn_qlw ,                         &   ! input fields
         &                 sn_tair, sn_humi, sn_prec, sn_snow, sn_slp,                &
         &                 sn_uoatm, sn_voatm, sn_cc, sn_hpgi, sn_hpgj

      ! cool-skin / warm-layer !LB
      !!---------------------------------------------------------------------
      !
      !                             !** read bulk namelist
      READ_NML_REF(numnam,namsbc_blk)
      READ_NML_CFG(numnam,namsbc_blk)
      IF(lwm) WRITE( numond, namsbc_blk )
      !
      !                             !** initialization of the chosen bulk formulae (+ check)
      !                                   !* select the bulk chosen in the namelist and check the choice
      ioptio = 0
      IF( ln_NCAR      ) THEN
         nblk =  np_NCAR        ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_COARE_3p0 ) THEN
         nblk =  np_COARE_3p0   ;   ioptio = ioptio + 1
      ENDIF
      !IF( ln_COARE_3p6 ) THEN
      !   nblk =  np_COARE_3p6   ;   ioptio = ioptio + 1
      !ENDIF
      IF( ln_ECMWF     ) THEN
         nblk =  np_ECMWF       ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one bulk algorithm' )


      ll_skin = ln_skin_cs .OR. ln_skin_wl

      IF(ll_skin) PRINT *, 'LOLO [sbcblk.F90] => ll_skin =', ll_skin

      !                             !** initialization of the cool-skin / warm-layer parametrization
      IF( ll_skin ) THEN
         !! Some namelist sanity tests:
         IF( ln_NCAR )      &
            & CALL ctl_stop( 'sbc_blk_init: Cool-skin/warm-layer param. not compatible with NCAR algorithm' )
         !IF( nn_fsbc /= 1 ) &
         !   & CALL ctl_stop( 'sbc_blk_init: Please set "nn_fsbc" to 1 when using cool-skin/warm-layer param.')
      ENDIF

      IF( ln_skin_wl ) THEN
         !! Check if the frequency of downwelling solar flux input makes sense and if ln_dm2dc=T if it is daily!
         IF( (sn_qsr%freqh  < 0.).OR.(sn_qsr%freqh  > 24.) ) &
            & CALL ctl_stop( 'sbc_blk_init: Warm-layer param. (ln_skin_wl) not compatible with freq. of solar flux > daily' )
         IF( (sn_qsr%freqh == 24.).AND.(.NOT. ln_dm2dc) ) &
            & CALL ctl_stop( 'sbc_blk_init: Please set ln_dm2dc=T for warm-layer param. (ln_skin_wl) to work properly' )
      ENDIF

      ioptio = 0
      IF( ln_humi_sph ) THEN
         nhumi =  np_humi_sph    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_dpt ) THEN
         nhumi =  np_humi_dpt    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_rlh ) THEN
         nhumi =  np_humi_rlh    ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one type of air humidity' )
      !
      IF( ln_dm2dc ) THEN                 !* check: diurnal cycle on Qsr
         IF( sn_qsr%freqh /= 24. )   CALL ctl_stop( 'sbc_blk_init: ln_dm2dc=T only with daily short-wave input' )
         IF( sn_qsr%ln_tint ) THEN
            CALL ctl_warn( 'sbc_blk_init: ln_dm2dc=T daily qsr time interpolation done by sbcdcy module',   &
               &           '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
      ENDIF

      ioptio = 0
      IF( ln_Cx_ice_cst ) THEN
         nblk_ice =  np_ice_cst     ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_Cx_ice_LU12 ) THEN
         nblk_ice =  np_ice_lu12    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_Cx_ice_LG15 ) THEN
         nblk_ice =  np_ice_lg15   ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one ice-atm bulk algorithm' )


      !                                   !* set the bulk structure
      !                                      !- store namelist information in an array
      !
      slf_i(jp_wndi ) = sn_wndi    ;   slf_i(jp_wndj ) = sn_wndj
      slf_i(jp_qsr  ) = sn_qsr     ;   slf_i(jp_qlw  ) = sn_qlw
      slf_i(jp_tair ) = sn_tair    ;   slf_i(jp_humi ) = sn_humi
      slf_i(jp_prec ) = sn_prec    ;   slf_i(jp_snow ) = sn_snow
      slf_i(jp_slp  ) = sn_slp     ;   slf_i(jp_cc   ) = sn_cc
      slf_i(jp_uoatm) = sn_uoatm   ;   slf_i(jp_voatm) = sn_voatm
      slf_i(jp_hpgi ) = sn_hpgi    ;   slf_i(jp_hpgj ) = sn_hpgj
      !
      IF( .NOT. ln_abl ) THEN   ! force to not use jp_hpgi and jp_hpgj, should already be done in namelist_* but we never know...
         slf_i(jp_hpgi)%clname = 'NOT USED'
         slf_i(jp_hpgj)%clname = 'NOT USED'
      ENDIF
      !
      !                                      !- allocate the bulk structure
      ALLOCATE( sf(jpfld), STAT=ierror )
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_init: unable to allocate sf structure' )
      !
      !                                      !- fill the bulk structure with namelist informations
      CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_init', 'surface boundary condition -- bulk formulae', 'namsbc_blk' )
      sf(jp_wndi )%zsgn = -1._wp   ;   sf(jp_wndj )%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      sf(jp_uoatm)%zsgn = -1._wp   ;   sf(jp_voatm)%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      sf(jp_hpgi )%zsgn = -1._wp   ;   sf(jp_hpgj )%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
      !
      DO jfpr= 1, jpfld
         !
         IF(   ln_abl    .AND.                                                      &
            &    ( jfpr == jp_wndi .OR. jfpr == jp_wndj .OR. jfpr == jp_humi .OR.   &
            &      jfpr == jp_hpgi .OR. jfpr == jp_hpgj .OR. jfpr == jp_tair     )  ) THEN
            ipka = jpka   ! ABL: some fields are 3D input
         ELSE
            ipka = 1
         ENDIF
         !
         ALLOCATE( sf(jfpr)%fnow(jpi,jpj,ipka) )
         !
         IF( TRIM(sf(jfpr)%clrootname) == 'NOT USED' ) THEN    !--  not used field  --!   (only now allocated and set to default)
            IF(     jfpr == jp_slp ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 101325._wp   ! use standard pressure in Pa
            ELSEIF( jfpr == jp_prec .OR. jfpr == jp_snow .OR. jfpr == jp_uoatm .OR. jfpr == jp_voatm ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 0._wp        ! no precip or no snow or no surface currents
            ELSEIF( jfpr == jp_wndi .OR. jfpr == jp_wndj ) THEN
               sf(jfpr)%fnow(:,:,1:ipka) = 0._wp
            ELSEIF( jfpr == jp_hpgi .OR. jfpr == jp_hpgj ) THEN
               IF( .NOT. ln_abl ) THEN
                  DEALLOCATE( sf(jfpr)%fnow )   ! deallocate as not used in this case
               ELSE
                  sf(jfpr)%fnow(:,:,1:ipka) = 0._wp
               ENDIF
            ELSEIF( jfpr == jp_cc  ) THEN
               sf(jp_cc)%fnow(:,:,1:ipka) = pp_cldf
            ELSE
               WRITE(ctmp1,*) 'sbc_blk_init: no default value defined for field number', jfpr
               CALL ctl_stop( ctmp1 )
            ENDIF
         ELSE                                                  !-- used field  --!
            IF( sf(jfpr)%ln_tint )   ALLOCATE( sf(jfpr)%fdta(jpi,jpj,ipka,2) )   ! allocate array for temporal interpolation
            !
            IF( sf(jfpr)%freqh > 0. .AND. MOD( NINT(3600. * sf(jfpr)%freqh), nn_fsbc * NINT(rn_Dt) ) /= 0 )   &
               &  CALL ctl_warn( 'sbc_blk_init: sbcmod timestep rn_Dt*nn_fsbc is NOT a submultiple of atmospheric forcing frequency.',   &
               &                 '               This is not ideal. You should consider changing either rn_Dt or nn_fsbc value...' )
         ENDIF
      END DO
      !
      IF( ln_abl ) THEN       ! ABL: read 3D fields for wind, temperature, humidity and pressure gradient
         rn_zqt = ght_abl(2)          ! set the bulk altitude to ABL first level
         rn_zu  = ght_abl(2)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ABL formulation: overwrite rn_zqt & rn_zu with ABL first level altitude'
      ENDIF
      !
      !
      IF(lwp) THEN                     !** Control print
         !
         WRITE(numout,*)                  !* namelist
         WRITE(numout,*) '   Namelist namsbc_blk (other than data information):'
         WRITE(numout,*) '      "NCAR"      algorithm   (Large and Yeager 2008)      ln_NCAR      = ', ln_NCAR
         WRITE(numout,*) '      "COARE 3.0" algorithm   (Fairall et al. 2003)       ln_COARE_3p0 = ', ln_COARE_3p0
         !WRITE(numout,*) '      "COARE 3.6" algorithm (Fairall 2018 + Edson al 2013) ln_COARE_3p6 = ', ln_COARE_3p6
         WRITE(numout,*) '      "ECMWF"     algorithm   (IFS cycle 45r1)             ln_ECMWF     = ', ln_ECMWF
         WRITE(numout,*) '      Air temperature and humidity reference height (m)   rn_zqt       = ', rn_zqt
         WRITE(numout,*) '      Wind vector reference height (m)                    rn_zu        = ', rn_zu
         WRITE(numout,*) '         (form absolute (=0) to relative winds(=1))'
         WRITE(numout,*) '      use surface current feedback on wind stress         ln_crt_fbk   = ', ln_crt_fbk
         IF(ln_crt_fbk) THEN
            WRITE(numout,*) '         Renault et al. 2020, eq. 10: Stau = Alpha * Wnd + Beta'
            WRITE(numout,*) '            Alpha                                         rn_stau_a    = ', rn_stau_a
            WRITE(numout,*) '            Beta                                          rn_stau_b    = ', rn_stau_b
         ENDIF
         !
         WRITE(numout,*)
         SELECT CASE( nblk )              !* Print the choice of bulk algorithm
         CASE( np_NCAR      )   ;   WRITE(numout,*) '   ==>>>   "NCAR" algorithm        (Large and Yeager 2008)'
         CASE( np_COARE_3p0 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.0" algorithm   (Fairall et al. 2003)'
         !CASE( np_COARE_3p6 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.6" algorithm (Fairall 2018+Edson et al. 2013)'
         CASE( np_ECMWF     )   ;   WRITE(numout,*) '   ==>>>   "ECMWF" algorithm       (IFS cycle 45r1)'
         END SELECT
         !
         WRITE(numout,*)
         WRITE(numout,*) '      use cool-skin  parameterization (SSST)  ln_skin_cs  = ', ln_skin_cs
         WRITE(numout,*) '      use warm-layer parameterization (SSST)  ln_skin_wl  = ', ln_skin_wl
         !
         WRITE(numout,*)
         SELECT CASE( nhumi )              !* Print the choice of air humidity
         CASE( np_humi_sph )   ;   WRITE(numout,*) '   ==>>>   air humidity is SPECIFIC HUMIDITY     [kg/kg]'
         CASE( np_humi_dpt )   ;   WRITE(numout,*) '   ==>>>   air humidity is DEW-POINT TEMPERATURE [K]'
         CASE( np_humi_rlh )   ;   WRITE(numout,*) '   ==>>>   air humidity is RELATIVE HUMIDITY     [%]'
         END SELECT
         !
         WRITE(numout,*)
         WRITE(numout,*) '      use constant ice-atm bulk transfer coeff.           ln_Cx_ice_cst  = ', ln_Cx_ice_cst
         WRITE(numout,*) '      use ice-atm bulk coeff. from Lupkes et al., 2012    ln_Cx_ice_LU12 = ', ln_Cx_ice_LU12
         WRITE(numout,*) '      use ice-atm bulk coeff. from Lupkes & Gryanik, 2015 ln_Cx_ice_LG15 = ', ln_Cx_ice_LG15
         WRITE(numout,*)
         SELECT CASE( nblk_ice )              !* Print the choice of bulk algorithm
         CASE( np_ice_cst  )
            WRITE(numout,*) '   ==>>>   Constant bulk transfer coefficients over sea-ice:'
            WRITE(numout,*) '      => Cd_ice, Ce_ice, Ch_ice =', REAL(rn_Cd_i,4), REAL(rn_Ce_i,4), REAL(rn_Ch_i,4)
            IF( (rn_Cd_i<0._wp).OR.(rn_Cd_i>1.E-2_wp).OR.(rn_Ce_i<0._wp).OR.(rn_Ce_i>1.E-2_wp).OR.(rn_Ch_i<0._wp).OR.(rn_Ch_i>1.E-2_wp) ) &
               & CALL ctl_stop( 'Be realistic in your pick of Cd_ice, Ce_ice & Ch_ice ! (0 < Cx < 1.E-2)')
         CASE( np_ice_an05 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Andreas et al, 2005'
         CASE( np_ice_lu12 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Lupkes et al, 2012'
         CASE( np_ice_lg15 )   ;   WRITE(numout,*) '   ==>>> bulk algo over ice: Lupkes & Gryanik, 2015'
         END SELECT
         !#LB.
         !
      ENDIF
      !$acc update device( rn_zqt, rn_zu )
      !
   END SUBROUTINE sbc_blk_init


   SUBROUTINE sbc_blk( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!              (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  :
      !!              (1) READ each fluxes in NetCDF files:
      !!      the wind velocity (i-component) at z=rn_zu  (m/s) at T-point
      !!      the wind velocity (j-component) at z=rn_zu  (m/s) at T-point
      !!      the specific humidity           at z=rn_zqt (kg/kg)
      !!      the air temperature             at z=rn_zqt (Kelvin)
      !!      the solar heat                              (W/m2)
      !!      the Long wave                               (W/m2)
      !!      the total precipitation (rain+snow)         (Kg/m2/s)
      !!      the snow (solid precipitation)              (kg/m2/s)
      !!      ABL dynamical forcing (i/j-components of either hpg or geostrophic winds)
      !!              (2) CALL blk_oce_1 and blk_oce_2
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress at T-point
      !!              - taum        wind stress module at T-point
      !!              - wndm        wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-solar and solar heat fluxes
      !!              - emp         upward mass flux (evapo. - precip.)
      !!              - sfx         salt flux due to freezing/melting (non-zero only if ice is present)
      !!
      !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !! ==> updates the following arrays:
      !!       emp, qsr, qns, qns_oce, qsr_oce, wndm, utau, vtau, taum, rhoa
      !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!                   Brodeau et al. Ocean Modelling 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zcd_du, zsen, zlat, zevap, zpre
      REAL(wp) :: ztst, zpa
      LOGICAL  :: llerr
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('sbc_blk')  !lili
      !
      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step

      ! Sanity/consistence test on humidity at first time step to detect potential screw-up:
      IF( kt == nit000 ) THEN
         ! mean humidity over ocean on proc
         ztst = glob_sum( 'sbcblk', sf(jp_humi)%fnow(:,:,1) * e1e2t(:,:) * xmskt(:,:) ) / glob_sum( 'sbcblk', e1e2t(:,:) * xmskt(:,:) )
         llerr = .FALSE.
         SELECT CASE( nhumi )
         CASE( np_humi_sph ) ! specific humidity => expect: 0. <= something < 0.065 [kg/kg] (0.061 is saturation at 45degC !!!)
            IF( (ztst <   0._wp) .OR. (ztst > 0.065_wp) )   llerr = .TRUE.
         CASE( np_humi_dpt ) ! dew-point temperature => expect: 110. <= something < 320. [K]
            IF( (ztst < 110._wp) .OR. (ztst >  320._wp) )   llerr = .TRUE.
         CASE( np_humi_rlh ) ! relative humidity => expect: 0. <= something < 100. [%]
            IF( (ztst <   0._wp) .OR. (ztst >  100._wp) )   llerr = .TRUE.
         END SELECT
         IF(llerr) THEN
            WRITE(ctmp1,'("   Error on mean humidity value: ",f10.5)') ztst
            CALL ctl_stop( 'STOP', ctmp1, 'Something is wrong with air humidity!!!', &
               &   ' ==> check the unit in your input files'       , &
               &   ' ==> check consistence of namelist choice: specific? relative? dew-point?', &
               &   ' ==> ln_humi_sph -> [kg/kg] | ln_humi_rlh -> [%] | ln_humi_dpt -> [K] !!!' )
         ENDIF
         IF(lwp) THEN
            WRITE(numout,*) ''
            WRITE(numout,*) ' Global mean humidity at kt = nit000: ', ztst
            WRITE(numout,*) ' === Sanity/consistence test on air humidity sucessfuly passed! ==='
            WRITE(numout,*) ''
         ENDIF
      ENDIF   !IF( kt == nit000 )
      !                                            ! compute the surface ocean fluxes using bulk formulea


      !! CPU STUFF:
      fatm_slp(:,:)   = sf(jp_slp )%fnow(:,:,1)
      fatm_theta(:,:) = sf(jp_tair)%fnow(:,:,1)
      fatm_q(:,:)     = sf(jp_humi)%fnow(:,:,1)
      fatm_u(:,:)     = sf(jp_wndi)%fnow(:,:,1)
      fatm_v(:,:)     = sf(jp_wndj)%fnow(:,:,1)
      fatm_prcp(:,:)  = sf(jp_prec)%fnow(:,:,1)
      fatm_snow(:,:)  = sf(jp_snow)%fnow(:,:,1)
      fatm_dqsw(:,:)  = sf(jp_qsr )%fnow(:,:,1)
      fatm_dqlw(:,:)  = sf(jp_qlw )%fnow(:,:,1)
      !$acc update device( fatm_slp, fatm_theta, fatm_q, fatm_u, fatm_v, fatm_prcp, fatm_snow, fatm_dqsw, fatm_dqlw )

      IF( iom_use('snowpre') ) CALL iom_put( 'snowpre', fatm_prcp )                  ! Snow precipitation
      IF( iom_use('precip' ) ) CALL iom_put( 'precip' , fatm_snow )                  ! Total precipitation


      !$acc data create( zsen, zlat, zevap ) present( qsr, wndm, utau, vtau, ssst, taum, rhoa, emp, qns, qns_oce, qsr_oce )

      ! Specific humidity of air at z=rn_zqt
      SELECT CASE( nhumi )
      CASE( np_humi_dpt )
         IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => computing q_air out of dew-point and P !'
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               fatm_q(ji,jj) = q_sat( fatm_q(ji,jj), fatm_slp(ji,jj) )
            ENDDO
         ENDDO
         !$acc end parallel loop
      CASE( np_humi_rlh )
         IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => computing q_air out of RH, t_air and slp !' !LBrm
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               fatm_q(ji,jj) = q_air_rh( 0.01_wp*fatm_q(ji,jj), fatm_theta(ji,jj), fatm_slp(ji,jj) ) !#LB: 0.01 => RH is % percent in file
            ENDDO
         ENDDO
         !$acc end parallel loop
      END SELECT

      ! Potential temperature of air at z=rn_zqt (most reanalysis products provide absolute temp., not potential temp.)
      IF( .NOT. ln_tair_pot ) THEN
         ! temperature read into file is ABSOLUTE temperature (that's the case for ECMWF products for example...)
         IF((kt==nit000).AND.lwp) WRITE(numout,*) ' *** sbc_blk() => air temperature converted from ABSOLUTE to POTENTIAL!'
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zpa = pres_temp( fatm_q(ji,jj), fatm_slp(ji,jj), rn_zqt, pta=fatm_theta(ji,jj) )
               fatm_theta(ji,jj) = theta_exner( fatm_theta(ji,jj), zpa )
            ENDDO
         ENDDO
         !$acc end parallel loop
         !
      ENDIF


      !DEBUG/SANITY:
      !IF( .NOT. ln_slab_sst ) THEN
      !   IF( SUM(ABS(sst_m-sst_s)) > 0.001_wp ) THEN
      !      PRINT *, 'ERROR: SLAB not used, yet (SSS_s /= sst_m) !!! kt=',kt
      !      PRINT *, ' ==> `sst_s-sst_m` ='
      !      PRINT *, sst_s-sst_m
      !      PRINT *, ''
      !      STOP
      !   ENDIF
      !ENDIF
      !DEBUG/SANITY.

      
      !DEBUG:
      !ji=22 ; jj=120
      !PRINT *, 'BULK INPUT:'
      !PRINT *, ' * SST   =', sst_m(ji,jj)
      !PRINT *, ' * t_air =', fatm_theta(ji,jj)-rt0  ! => C
      !PRINT *, ' * q_air =', fatm_q(ji,jj)*1000._wp ! => g/kg
      !PRINT *, ' * slp   =', fatm_slp(ji,jj)/100.   ! => hPa
      !PRINT *, ' * wspd  =', SQRT( fatm_u(ji,jj)*fatm_u(ji,jj) + fatm_v(ji,jj)*fatm_v(ji,jj) )
      !ENDIF

      !! REMINDER:
      !! - at this stage, `sst_s` is the SLAB bulk SST if the slab ocean is used, otherwize it is the bulk SST (`sst_m`) !
      

      CALL blk_oce_1( kt, fatm_u, fatm_v, fatm_theta, fatm_q, fatm_slp,       &   !   <<= in
         &                sst_s, ssu_m, ssv_m,                                &   !   <<= in
         &                sf(jp_uoatm)%fnow(:,:,1), sf(jp_voatm)%fnow(:,:,1), &   !   <<= in
         &                fatm_dqsw, fatm_dqlw,                               &   !   <<= in (wl/cs)
         &                ssst, sssq, zsen, zlat, zevap, qsr, wndm, taum, utau, vtau ) ! =>> in/out or out
      !&                 zcd_du, zsen, zlat, zevap )                    !   =>> out
      !!
      !!    ==> updates: qsr, wndm, utau, vtau, ssst, taum, rhoa, zsen, zlat, zevap

      !ji=22 ; jj=120
      !PRINT *, 'BULK OUTPUT:'
      !PRINT *, ' * Qsen   =', zsen(ji,jj)
      !PRINT *, ' * Qlat   =', zlat(ji,jj)
      !PRINT *, ' * taum   =', taum(ji,jj)*1000. ! => mN/m^2
      !PRINT *, ' * ssst   =', ssst(ji,jj)
      !PRINT *, ' * sssq   =', sssq(ji,jj)
      !PRINT *, ''; PRINT *, ''


      CALL blk_oce_2( fatm_theta, fatm_dqlw, fatm_prcp, fatm_snow, ssst, zsen, zlat, zevap ) !   <<= in
      !!   ==> updates: emp, qns, qns_oce, qsr_oce

      !$acc end data
      IF( ln_timing )   CALL timing_stop('sbc_blk')
      !
   END SUBROUTINE sbc_blk


   SUBROUTINE blk_oce_1( kt, pwndi, pwndj, ptair, pqair, pslp, psst, pu, pv,     &  ! in
      &                      puatm, pvatm, pdqsr, pdqlw,                         &  ! in
      &                      pssst, psssq,                                        &  ! in/out
      &                      psen, plat, pevap, pqsr, pwndm, ptaum, putau, pvtau  )  ! out
      !&                             pcd_du, psen, plat, pevap ) ! out
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_1  ***
      !!
      !! ** Purpose :   if ln_blk=T, computes surface momentum, heat and freshwater fluxes
      !!                if ln_abl=T, computes Cd x |U|, Ch x |U|, Ce x |U| for ABL integration
      !!
      !! ** Method  :   bulk formulae using atmospheric fields from :
      !!                if ln_blk=T, atmospheric fields read in sbc_read
      !!                if ln_abl=T, the ABL model at previous time-step
      !!
      !! ** Outputs : - psssq    : surface humidity used to compute latent heat flux (kg/kg)
      !!              - pcd_du  : Cd x |dU| at T-points  (m/s)
      !!              - psen    : sensible heat flux (W/m^2)
      !!              - plat    : latent heat flux   (W/m^2)
      !!              - pevap    : evaporation        (mm/s) #lolo
      !!              - pqsr    : net shortwave radiation available for the ocean (after albedo) (W/m^2)
      !!              - pwndm   : Wind speed module at T-point          (m/s)
      !!              - ptaum   : module of air-sea wind stress at T-point (N/m^2)
      !!              - putau   : i-component of the stress at T-point     (N/m2)
      !!              - pvtau   : j-component of the stress at T-point     (N/m2)
      !!
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in   )                     ::   kt     ! time step index
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pwndi  ! atmospheric wind at T-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pwndj  ! atmospheric wind at T-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pqair  ! specific humidity at T-points            [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ptair  ! potential temperature at T-points        [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pslp   ! sea-level pressure                       [Pa]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   psst   ! (always) BULK SST                        [Celsius]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pu     ! surface current at U-point (i-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pv     ! surface current at V-point (j-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   puatm  ! surface current seen by the atm at T-point (i-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pvatm  ! surface current seen by the atm at T-point (j-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pdqsr  ! downwelling solar (shortwave) radiation at surface [W/m^2]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   pdqlw  ! downwelling longwave radiation at surface [W/m^2]
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj) ::   pssst   ! SKIN SST (or BULK if CS & WL not used)  [Celsius]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   psssq   ! SKIN SSQ (or BULK if CS & WL not used)  [Celsius]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   psen
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   plat
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pevap
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pqsr   !NEW
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pwndm  !NEW
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   ptaum
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   putau  !NEW
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   pvtau  !NEW
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zmsk, z1_alb, zztmp  ! local variable
      REAL(wp) ::   zstmax, zstau, zpatm
      !REAL(wp), DIMENSION(jpi,jpj) ::   ztau_i, ztau_j    ! wind stress components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   ztheta_zu, zq_zu, zU_zu   ! bulk wind speed at height zu  [m/s]
      REAL(wp), DIMENSION(jpi,jpj) ::   zcd_oce           ! momentum transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zch_oce           ! sensible heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zce_oce           ! latent   heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   ztabs       ! air pressure [Pa] & absolute temperature [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zztmp1, zztmp2
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('blk_oce_1')
      !
      !$acc data present( pwndi,pwndj,ptair,pqair,pslp,psst,pu,pv,pdqsr,pdqlw, pssst,psssq,psen,plat,pevap,ptaum,pqsr,pwndm ) create(zcd_oce,zch_oce,zce_oce,ztheta_zu,zq_zu,zU_zu)
      !

      z1_alb = 1. - albo

      ! Required so that no funny FP stuff occurs on the halos...
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pwndm(ji,jj) = 0._wp !lili
            psen(ji,jj)  = 0._wp
            plat(ji,jj)  = 0._wp
            ptaum(ji,jj) = 0._wp
            pevap(ji,jj)  = 0._wp
            !
            pqsr(ji,jj)  = 0._wp
         END DO
      END DO
      !$acc end parallel loop


      !! Probably not needed, but:
      IF( .NOT. ll_skin ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               pssst(ji,jj) = psst(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF





      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ! local scalars ( place there for vector optimisation purposes)
            !                           ! Temporary conversion from Celcius to Kelvin (and set minimum value far above 0 K)
            !pssst(ji,jj) = psst(ji,jj) + rt0  ! by default: skin temperature = "bulk SST" (will remain this way if NCAR algorithm used!)

            ! sea surface potential temperature [K]
            !zsspt(ji,jj) = theta_exner( pssst(ji,jj), pslp(ji,jj) )

            ! ----------------------------------------------------------------------------- !
            !      0   Wind components and module at T-point relative to the moving ocean   !
            ! ----------------------------------------------------------------------------- !

            ! ... components ( U10m - U_oce ) at T-point (unmasked)
            ! ... scalar wind module at T-point (not masked)
            pwndm(ji,jj) = SQRT(  pwndi(ji,jj) * pwndi(ji,jj) + pwndj(ji,jj) * pwndj(ji,jj)  )


            ! ----------------------------------------------------------------------------- !
            !      I   Solar FLUX                                                           !
            ! ----------------------------------------------------------------------------- !

            ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave

            !IF( ln_dm2dc ) THEN !#LOLOfixme
            !   pqsr(:,:) = z1_alb * sbc_dcy( pdqsr(:,:) ) * xmskt(:,:)
            !ELSE
            pqsr(ji,jj) = z1_alb * pdqsr(ji,jj) * xmskt(ji,jj)
            !ENDIF

            ! ----------------------------------------------------------------------------- !
            !     II   Turbulent FLUXES                                                     !
            ! ----------------------------------------------------------------------------- !

            ! specific humidity at SKIN SST (if `ll_skin`) or BULK SST (if `.NOT. ll_skin`):
            psssq(ji,jj) = rdct_qsat_salt * q_sat( pssst(ji,jj)+rt0 , pslp(ji,jj) )


         END DO !DO ji=Nis0, Nie0
      END DO !DO jj=Njs0, Nje0
      !$acc end parallel loop

      ! transfer coefficients (Cd, Ch, Ce at T-point, and more)
      SELECT CASE( nblk )   ! user-selected bulk parameterization
         !
      CASE( np_NCAR      )
         CALL turb_ncar    (     rn_zqt, rn_zu, psst, ptair, psssq, pqair, pwndm,   &
            &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu ,  &
            &                nb_iter=nn_iter_algo )
         !
      CASE( np_COARE_3p0 )
         CALL turb_coare3p0( kt, rn_zqt, rn_zu, psst, pssst, ptair, psssq, pqair, pwndm,  &
            &                ln_skin_cs, ln_skin_wl,                             &
            &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu, &
            &                nb_iter=nn_iter_algo ) !LOLOfixme: no WL/CS for now...
         !
         !CASE( np_COARE_3p6 )
         !   CALL turb_coare3p6( kt, rn_zqt, rn_zu, psst, ptair, psssq, pqair, pwndm, &
         !      &                ln_skin_cs, ln_skin_wl,                             &
         !      &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu, &
         !      &                nb_iter=nn_iter_algo,                               &
         !      &                Qsw=pqsr(:,:), rad_lw=pdqlw(:,:), slp=pslp(:,:) )
         !
      CASE( np_ECMWF     )
         CALL turb_ecmwf   ( kt, rn_zqt, rn_zu, psst, pssst, ptair, psssq, pqair, pwndm,  &
            &                ln_skin_cs, ln_skin_wl,                             &
            &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu, &
            &                nb_iter=nn_iter_algo ) !LOLOfixme: no WL/CS for now...
         !
         !CALL turb_ecmwf   ( kt, rn_zqt, rn_zu, psst, ptair, psssq, pqair, pwndm, &
         !   &                ln_skin_cs, ln_skin_wl,                            &
         !   &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu,  &
         !   &                nb_iter=nn_iter_algo,                              &
         !   &                Qsw=pqsr(:,:), rad_lw=pdqlw(:,:), slp=pslp(:,:) )
         !
         !CASE( np_ANDREAS   )
         !   CALL turb_andreas (     rn_zqt, rn_zu, psst, ptair, psssq, pqair, pwndm, &
         !      &                zcd_oce, zch_oce, zce_oce, ztheta_zu, zq_zu, zU_zu , &
         !      &                nb_iter=nn_iter_algo   )
         !
      CASE DEFAULT
         CALL ctl_stop( 'STOP', 'sbc_oce: non-existing bulk parameterizaton selected' )
      END SELECT


      !DEBUG:
      !ji=22 ; jj=120
      !PRINT *, 'BULK OUTPUT:'
      !PRINT *, ' * C_D   =', zcd_oce(ji,jj) * 1000.
      !PRINT *, ' * C_E   =', zce_oce(ji,jj) * 1000.
      !PRINT *, ' * C_H   =', zch_oce(ji,jj) * 1000.
      !PRINT *, ' * t_zu  =', ztheta_zu(ji,jj) -rt0
      !PRINT *, ' * q_zu  =', zq_zu(ji,jj) * 1000.
      !PRINT *, ''
      !DEBUG.

      !IF( iom_use('Cd_oce') )   CALL iom_put("Cd_oce",   zcd_oce * xmskt(:,:))
      !IF( iom_use('Ce_oce') )   CALL iom_put("Ce_oce",   zce_oce * xmskt(:,:))
      !IF( iom_use('Ch_oce') )   CALL iom_put("Ch_oce",   zch_oce * xmskt(:,:))
!!! LB: mainly here for debugging purpose:
      !IF( iom_use('theta_zt') ) CALL iom_put("theta_zt", (ptair-rt0) * xmskt(:,:)) ! potential temperature at z=zt
      !IF( iom_use('q_zt') )     CALL iom_put("q_zt",     pqair       * xmskt(:,:)) ! specific humidity       "
      !IF( iom_use('theta_zu') ) CALL iom_put("theta_zu", (ztheta_zu -rt0) * xmskt(:,:)) ! potential temperature at z=zu
      !IF( iom_use('q_zu') )     CALL iom_put("q_zu",     zq_zu        * xmskt(:,:)) ! specific humidity       "
      !IF( iom_use('ssq') )      CALL iom_put("ssq",      zssq        * xmskt(:,:)) ! saturation specific humidity at z=0
      !IF( iom_use('wspd_blk') ) CALL iom_put("wspd_blk", zU_zu       * xmskt(:,:)) ! bulk wind speed at z=zu

      !LOLOfixme:
      !IF( ll_skin ) THEN
      !   !! In the presence of sea-ice we forget about the cool-skin/warm-layer update of zsspt, zssq & pssst:
      !   WHERE( at_i(:,:) > 0.001_wp )
      !      ! sea-ice present, we forget about the update, using what we backed up before call to turb_*()
      !      zsspt(:,:) = zztmp1(:,:)
      !      zssq(:,:)  = zztmp2(:,:)
      !   END WHERE
      !   ! apply potential temperature increment to abolute SST
      !   pssst(:,:) = pssst(:,:) + ( zsspt(:,:) - zztmp1(:,:) )
      !ENDIF
      !LOLOfixme.

      !  Turbulent fluxes over ocean  => BULK_FORMULA @ sbc_phy.F90
      ! -------------------------------------------------------------

      !IF( ln_abl ) THEN         !==  ABL formulation  ==!   multiplication by rho_air and turbulent fluxes computation done in ablstp
      !
      !   DO jj=Njs0, Nje0
      !      DO ji=Nis0, Nie0
      !         zztmp = zU_zu(ji,jj)
      !         pwndm(ji,jj)   = zztmp                   ! Store zU_zu in pwndm to compute ustar2 in ablmod
      !         pcd_du(ji,jj) = zztmp * zcd_oce(ji,jj)
      !         psen(ji,jj)   = zztmp * zch_oce(ji,jj)
      !         pevap(ji,jj)   = zztmp * zce_oce(ji,jj)
      !         zpre(ji,jj)   = pres_temp( pqair(ji,jj), pslp(ji,jj), rn_zu, ptpot=ptair(ji,jj), pta=ztabs(ji,jj) )
      !         rhoa(ji,jj)   = rho_air( ztabs(ji,jj), pqair(ji,jj), zpre(ji,jj) )
      !      END DO
      !   END DO
      !
      !ELSE                      !==  BLK formulation  ==!   turbulent fluxes computation


      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zmsk = xmskt(ji,jj)
            
            !zpatm       = pres_temp( zq_zu(ji,jj), pslp(ji,jj), rn_zu, ptpot=ztheta_zu(ji,jj), pta=ztabs(ji,jj) )
            !rhoa(ji,jj) = rho_air( ztabs(ji,jj), zq_zu(ji,jj), zpatm ) !LOLObug: no need as returned by `bulk_formula_sclr` !!!

            CALL bulk_formula_sclr( rn_zu, pssst(ji,jj), psssq(ji,jj), ztheta_zu(ji,jj), zq_zu(ji,jj), &
               &                    zcd_oce(ji,jj), zch_oce(ji,jj), zce_oce(ji,jj),               &
               &                    pwndm(ji,jj), zU_zu(ji,jj), pslp(ji,jj),                       &
               &                    pTau=ptaum(ji,jj), pQsen=psen(ji,jj), pQlat=plat(ji,jj),       &
               &                    pEvap=pevap(ji,jj), prhoa=rhoa(ji,jj) )

            !DEBUG:
            !IF( ji==22 .AND. jj==120 ) THEN
            !   PRINT *, 'BULK_FORMULA OUTPUT:'
            !   PRINT *, ' * Rho_air =', rhoa(ji,jj)
            !   PRINT *, ' * E       =', pevap(ji,jj)
            !   PRINT *, ' * Qsens   =', psen(ji,jj)
            !   PRINT *, ' * Qlat    =', plat(ji,jj)
            !   PRINT *, ''
            !ENDIF
            !DEBUG.
            
            psen(ji,jj)  = psen(ji,jj)  * zmsk
            plat(ji,jj)  = plat(ji,jj)  * zmsk
            ptaum(ji,jj) = ptaum(ji,jj) * zmsk
            pevap(ji,jj)  = pevap(ji,jj)  * zmsk
            rhoa(ji,jj)  = rhoa(ji,jj)  * zmsk
            !
         END DO
      END DO
      !$acc end parallel loop

      !! `utau` & `vtau` at T-points:
      !! ----------------------------
      !$acc parallel loop collapse(2) present(putau,pvtau)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            zztmp = ptaum(ji,jj)/MAX(pwndm(ji,jj),1.E-12_wp)*xmskt(ji,jj)
            putau(ji,jj) =  zztmp*pwndi(ji,jj)
            pvtau(ji,jj) =  zztmp*pwndj(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      ! Saving open-ocean air-water fluxes:
      ! -----------------------------------
      IF( iom_use('taum_oce' ) ) THEN
         !$acc update self( ptaum )
         CALL iom_put( "taum_oce",   ptaum )   ! output wind stress module
      ENDIF
      IF( iom_use('rho_air' ) ) THEN
         !$acc update self( rhoa )
         CALL iom_put( "rho_air"  ,  rhoa )     ! output air density [kg/m^3]
      ENDIF
      IF( iom_use('evap_oce' ) ) THEN
         !$acc update self( pevap )
         CALL iom_put( "evap_oce" ,  pevap )     ! evaporation
      ENDIF
      IF( iom_use('qsb_oce' ) ) THEN
         !$acc update self( psen )
         CALL iom_put( "qsb_oce"  ,  psen )     ! output downward sensible heat over the ocean
      ENDIF
      IF( iom_use('qla_oce' ) ) THEN
         !$acc update self( plat )
         CALL iom_put( "qla_oce"  ,  plat )     ! output downward latent   heat over the ocean
      ENDIF
      IF( iom_use('utau_oce' ) ) THEN
         !$acc update self( putau )
         CALL iom_put( "utau_oce" ,  putau )     ! output downward latent   heat over the ocean
      ENDIF
      IF( iom_use('vtau_oce' ) ) THEN
         !$acc update self( pvtau )
         CALL iom_put( "vtau_oce" ,  pvtau )     ! output downward latent   heat over the ocean
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl( tab2d_1=psssq, clinfo1=' blk_oce_1: sssq    : ', mask1=tmask )
         CALL prt_ctl( tab2d_1=pwndm, clinfo1=' blk_oce_1: pwndm   : ', mask1=tmask )
         CALL prt_ctl( tab2d_1=putau, clinfo1=' blk_oce_1: putau   : ', mask1=umask,   &
            &          tab2d_2=pvtau, clinfo2='            pvtau   : ', mask2=vmask )
         CALL prt_ctl( tab2d_1=zcd_oce, clinfo1=' blk_oce_1: Cd    : ', mask1=tmask )
      ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('blk_oce_1')
      !
   END SUBROUTINE blk_oce_1
   !qsr, wndm, utau, vtau

   SUBROUTINE blk_oce_2( ptair, pdqlw, pprec, psnow, pssst, psen, plat, pevap )   ! <<= in
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_2  ***
      !!
      !! ** Purpose :   finalize the momentum, heat and freshwater fluxes computation
      !!                at the ocean surface at each time step knowing Cd, Ch, Ce and
      !!                atmospheric variables (from ABL or external data)
      !!
      !! ** Outputs : - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!---------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   ptair   ! potential temperature of air #LB: confirm!
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   pdqlw   ! downwelling longwave radiation at surface [W/m^2]
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   pprec
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   psnow
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   pssst   ! SKIN surface temperature   [Celsius]
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   psen
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   plat
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::   pevap
      !REAL(wp), INTENT(in), OPTIONAL, DIMENSION(jpi,jpj) ::   qlwn   ! net longwave radiation at surface (MFS only) [W/m^2]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp,zz1,zz2,zz3, zmsk    ! local variable
      REAL(wp), DIMENSION(jpi,jpj) ::   zqlw              ! net long wave radiative heat flux
      REAL(wp)                     ::   zcptrain, zcptsnw, zcptn ! Heat content per unit mass (J/kg)
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('blk_oce_2')
      !
      !$acc data present(ptair, pdqlw, pprec, psnow, pssst, psen, plat, pevap, emp, qns, qns_oce, qsr_oce) create( zqlw )


      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zmsk = xmskt(ji,jj)

            ! Heat content per unit mass (J/kg)
            zcptrain = (     ptair(ji,jj)        - rt0 ) * rcp  * zmsk
            zcptsnw = ( MIN( ptair(ji,jj), rt0 ) - rt0 ) * rcpi * zmsk
            zcptn =          pssst(ji,jj)                 * rcp  * zmsk

            ! ----------------------------------------------------------------------------- !
            !     III    Net longwave radiative FLUX                                        !
            ! ----------------------------------------------------------------------------- !
            !! #LB: now moved after Turbulent fluxes because must use the skin temperature rather than bulk SST
            !! (pssst is skin temperature if ln_skin_cs==.TRUE. .OR. ln_skin_wl==.TRUE., bulk SST otherwize)
            zqlw(ji,jj) = qlw_net( pdqlw(ji,jj), pssst(ji,jj)+rt0 ) * zmsk

            ! ----------------------------------------------------------------------------- !
            !     IV    Total FLUXES                                                       !
            ! ----------------------------------------------------------------------------- !
            !
            emp(ji,jj) = ( pevap(ji,jj) - pprec(ji,jj) ) * zmsk             ! mass flux (evap. - precip.)
            !
            qns(ji,jj) = (   zqlw(ji,jj) + psen(ji,jj) + plat(ji,jj)    &   ! Downward Non Solar
               &           - psnow(ji,jj) *  rLfus                      &   ! remove latent melting heat for solid precip
               &           - pevap(ji,jj) * zcptn                       &   ! remove evap heat content at SST
               &           + ( pprec(ji,jj) - psnow(ji,jj) ) * zcptrain &   ! add liquid precip heat content at Tair
               &           + psnow(ji,jj) * zcptsnw                     &   ! add solid  precip heat content at min(Tair,Tsnow)
               &          ) * zmsk
            !
            qns_oce(ji,jj) = ( zqlw(ji,jj) + psen(ji,jj) + plat(ji,jj) ) * zmsk ! non solar without emp (only needed by SI3)
            !
            qsr_oce(ji,jj) =    qsr(ji,jj)
            !
         END DO
      END DO
      !$acc end parallel loop

      IF( iom_use('qlw_oce' ) ) THEN
         !$acc update self( zqlw )
         CALL iom_put( "qlw_oce"  ,  zqlw )
      ENDIF
      IF( iom_use('snowpre' ) ) THEN
         !$acc update self( psnow )
         CALL iom_put( 'snowpre',   psnow )    ! output solid precipitation [kg/m2/s]
      ENDIF
      IF( iom_use('precip' ) ) THEN
         !$acc update self( pprec )
         CALL iom_put( 'precip' ,   pprec )    ! output total precipitation [kg/m2/s]
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl(tab2d_1=zqlw , clinfo1=' blk_oce_2: zqlw  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=psen , clinfo1=' blk_oce_2: psen  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=plat , clinfo1=' blk_oce_2: plat  : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=qns  , clinfo1=' blk_oce_2: qns   : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=emp  , clinfo1=' blk_oce_2: emp   : ', mask1=tmask )
      ENDIF
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('blk_oce_2')
      !
   END SUBROUTINE blk_oce_2


   !!----------------------------------------------------------------------
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!----------------------------------------------------------------------

   SUBROUTINE blk_ice_1( pwndi, pwndj, ptair, pqair, pslp, puice, pvice, ptsui,   & ! inputs
      &                  pCHi, pCEi, ptheta_zu_i, pq_zu_i,                        & ! outputs
      &                  putaui, pvtaui, pseni, pevapi, pssqi, pcd_dui )             ! optional outputs
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_1  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum using bulk formulation
      !!                formulea, ice variables and read atmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   pwndi   ! atmospheric wind at T-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   pwndj   ! atmospheric wind at T-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   ptair   ! atmospheric potential temperature at T-point [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   pqair   ! atmospheric specific humidity at T-point [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   pslp    ! sea-level pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   puice   ! sea-ice velocity on I or C grid [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   pvice   ! "
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::   ptsui   ! sea-ice surface temperature [K]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCEi    ! evap/sublim. transfer coefficient  [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   ptheta_zu_i ! air temperature adjusted at heaight `zu` [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pq_zu_i      ! air spec. hum. adjusted at heaight `zu` [kg/kg]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   putaui  ! if ln_blk air-ice wind stress i-component at T-point [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   pvtaui  ! if ln_blk air-ice wind stress j-component at T-point [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   pseni   ! if ln_abl
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   pevapi   ! if ln_abl
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   pssqi   ! if ln_abl
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL ::   pcd_dui ! if ln_abl
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zootm_su                      ! sea-ice surface mean temperature
      REAL(wp) ::   zztmp1, zztmp2                ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) :: zCDi, ztmp1, ztmp3 ! temporary array
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('blk_ice_1')
      ! UPDATES: wndm_ice     , pCHi,pCEi,ptheta_zu_i,pq_zu_i
      !$acc data present( pwndi,pwndj,ptair,pqair,pslp,puice,pvice,ptsui,pCHi,pCEi,ptheta_zu_i,pq_zu_i,putaui,pvtaui,wndm_ice,rhoa ) create( ztmp1, ztmp3, zCDi )
      !
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ! ------------------------------------------------------------ !
            !    Wind module relative to the moving ice ( U10m - U_ice )   !
            ! ------------------------------------------------------------ !
            wndm_ice(ji,jj) = SQRT( pwndi(ji,jj) * pwndi(ji,jj) + pwndj(ji,jj) * pwndj(ji,jj) )

            IF( nblk_ice>1 ) ztmp1(ji,jj) = q_sat( ptsui(ji,jj), pslp(ji,jj), l_ice=.TRUE. ) ! temporary array for SSQ over ice

            IF( nblk_ice==np_ice_cst ) THEN
               ! Constant bulk transfer coefficients over sea-ice:
               zCDi(ji,jj) = rn_Cd_i
               pCHi(ji,jj) = rn_Ch_i
               pCEi(ji,jj) = rn_Ce_i
               ! no height adjustment, keeping zt values:
               ptheta_zu_i(ji,jj) = ptair(ji,jj)
               pq_zu_i(ji,jj)     = pqair(ji,jj)
            ENDIF

         END DO
      END DO
      !$acc end parallel loop

      ! sea-ice <-> atmosphere bulk transfer coefficients
      SELECT CASE( nblk_ice )

      CASE( np_ice_an05 )  ! calculate new drag from Lupkes(2015) equations
# if defined _OPENACC
         CALL ctl_stop( 'blk_ice_1: ADAPT `turb_ice_an05` for GPU before using it!')
# endif
         CALL turb_ice_an05( rn_zqt, rn_zu, ptsui, ptair, ztmp1, pqair, wndm_ice,       &
            &                      zCDi, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
         !!
      CASE( np_ice_lu12 )
# if defined _OPENACC
         CALL ctl_stop( 'blk_ice_1: ADAPT `turb_ice_lu12` for GPU before using it!')
# endif
         CALL turb_ice_lu12( rn_zqt, rn_zu, ptsui, ptair, ztmp1, pqair, wndm_ice, at_i, &
            &                      zCDi, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
         !!
      CASE( np_ice_lg15 )  ! calculate new drag from Lupkes(2015) equations
# if defined _OPENACC
         CALL ctl_stop( 'blk_ice_1: ADAPT `turb_ice_lg15` for GPU before using it!')
# endif
         CALL turb_ice_lg15( rn_zqt, rn_zu, ptsui, ptair, ztmp1, pqair, wndm_ice, at_i, &
            &                      zCDi, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
         !!
      END SELECT

      IF( iom_use('Cd_ice') ) THEN
         !$acc update self( zCDi )
         CALL iom_put("Cd_ice", zCDi*1000._wp*xmskt(:,:))
      ENDIF
      IF( iom_use('Ce_ice') ) THEN
         !$acc update self( pCEi )
         CALL iom_put("Ce_ice", pCEi*1000._wp*xmskt(:,:))
      ENDIF
      IF( iom_use('Ch_ice') ) THEN
         !$acc update self( pCHi )
         CALL iom_put("Ch_ice", pCHi*1000._wp*xmskt(:,:))
      ENDIF


      IF( ln_blk ) THEN
         ! ---------------------------------------------------------------- !
         !    Wind stress relative to nonmoving ice at T-points ( U10m )    !
         ! ---------------------------------------------------------------- !
         ! supress moving ice in wind stress computation as we don't know how to do it properly...
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zztmp1        = rhoa(ji,jj) * zCDi(ji,jj) * wndm_ice(ji,jj)
               putaui(ji,jj) = zztmp1 * pwndi(ji,jj)
               pvtaui(ji,jj) = zztmp1 * pwndj(ji,jj)
            END DO
         END DO
         !$acc end parallel loop

         IF(sn_cfctl%l_prtctl)  CALL prt_ctl( tab2d_1=putaui, clinfo1=' blk_ice: putaui : ', mask1=umask   &
            &                               , tab2d_2=pvtaui, clinfo2='          pvtaui : ', mask2=vmask )

      ELSE

         ! ln_abl
         !$acc parallel loop collapse(2) present( pcd_dui, pseni, pevapi, pssqi )
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               pcd_dui(ji,jj) = wndm_ice(ji,jj) * zCDi(ji,jj)
               pseni  (ji,jj) = wndm_ice(ji,jj) * pCHi(ji,jj)
               pevapi (ji,jj) = wndm_ice(ji,jj) * pCEi(ji,jj)
               !
               pssqi(ji,jj) = q_sat( ptsui(ji,jj), pslp(ji,jj), l_ice=.TRUE. ) ; ! more accurate way to obtain ssq !
            END DO
         END DO
         !$acc end parallel loop

      ENDIF ! ln_blk  / ln_abl

      IF(sn_cfctl%l_prtctl)  CALL prt_ctl(tab2d_1=wndm_ice, clinfo1=' blk_ice: wndm_ice : ', mask1=tmask )

      !$acc end data
      IF( ln_timing )   CALL timing_stop('blk_ice_1')
   END SUBROUTINE blk_ice_1


   SUBROUTINE blk_ice_2( ptsu, phs, phi, palb, ptair, pslp, pdqlw, pprec, psnow, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_2  ***
      !!
      !! ** Purpose :   provide the heat and mass fluxes at air-ice interface
      !!
      !! ** Method  :   compute heat and freshwater exchanged
      !!                between atmosphere and sea-ice using bulk formulation
      !!                formulea, ice variables and read atmmospheric fields.
      !!
      !! caution : the net upward water flux has with mm/day unit
      !!
      !! UPDATED ARRAYS:
      !! `devap_ice,dqla_ice,dqns_ice,emp_ice,emp_oce,evap_ice,qemp_ice,qemp_oce,
      !!  qevap_ice,qla_ice,qsb_ice,qns_ice,qprec_ice,qsr_ice,qlw_ice,qtr_ice_top`
      !!
      !!
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   ptsu   ! sea ice surface temperature [K]
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   phs    ! snow thickness
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   phi    ! ice thickness
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   palb   ! ice albedo (all skies)
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   ptair  ! potential temperature of air #LB: okay ???
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pslp
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pdqlw
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pprec
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   psnow
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pCEi    ! evap/sublim. transfer coefficient  [-]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   ptheta_zu_i ! air temperature adjusted at heaight `zu` [K]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pq_zu_i      ! air spec. hum. adjusted at heaight `zu` [kg/kg]
      !!
      INTEGER  ::   ji, jj, jl               ! dummy loop indices
      REAL(wp) ::   zst, zst3, zsq, zsipt, zmsk, zsum, zhi    ! local variable
      REAL(wp) ::   zcoef_dqlw, zcoef_dqla, zsnw, zsnwprc   !   -      -
      REAL(wp) ::   zztmp, zzblk, zztmp1, z1_rLsub, zA_b, z1mA_b   !   -      -
      !REAL(wp), DIMENSION(jpi,jpj,jpl), ALLOCATABLE ::   zmsk   ! temporary mask for prt_ctl
      REAL(wp)  ::   z_dqlw        ! long wave heat sensitivity over ice
      REAL(wp)  ::   z_dqsb        ! sensible  heat sensitivity over ice
      REAL(wp)  ::   zevap         ! evaporation and snw distribution after wind blowing (SI3)
      REAL(wp)  ::   ztri
      REAL(wp)  ::   zcptrain, zcptsnw, zcptn ! Heat content per unit mass (J/kg)
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('blk_ice_2')
      !$acc data present( ptsu,phs,phi,palb,ptair,pslp,pdqlw,pprec,psnow,pCHi,pCEi,ptheta_zu_i,pq_zu_i,sst_s,qsr_ice,qlw_ice,qsr,rhoa,wndm_ice,qla_ice,qsb_ice,dqla_ice,qns_ice,dqns_ice,qsr_oce )

      zcoef_dqlw = 4._wp * emiss_i * stefan             ! local scalars
      zztmp = 1. / ( 1. - albo )
      z1_rLsub = 1._wp / rLsub

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zmsk = xmskt(ji,jj)

            !$acc loop seq
            DO jl = 1, jpl                        !  Loop over ice categories  !

               zst   = ptsu(ji,jj,jl)                                ! surface temperature of sea-ice [K]
               zsq   = q_sat( zst, pslp(ji,jj), l_ice=.TRUE. )       ! surface saturation specific humidity when ice present
               zsipt = theta_exner( zst, pslp(ji,jj) )               ! potential sea-ice surface temperature [K]

               ! ----------------------------!
               !      I   Radiative FLUXES   !
               ! ----------------------------!
               ! Short Wave (sw)
               qsr_ice(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)

               ! Long  Wave (lw)
               zst3 = zst * zst * zst
               qlw_ice(ji,jj,jl)   = emiss_i * ( pdqlw(ji,jj) - stefan * zst * zst3 ) * zmsk
               ! lw sensitivity
               z_dqlw  = zcoef_dqlw * zst3

               ! ----------------------------!
               !     II    Turbulent FLUXES  !
               ! ----------------------------!

               ! ... turbulent heat fluxes with pCHi recalculated in blk_ice_1

               ! Common term in bulk F. equations...
               zzblk = rhoa(ji,jj) * wndm_ice(ji,jj)

               ! Sensible Heat
               zztmp1 = zzblk * rCp_air * pCHi(ji,jj)
               qsb_ice(ji,jj,jl) = zztmp1 * (zsipt - ptheta_zu_i(ji,jj))
               z_dqsb = zztmp1                        ! ==> Qsens sensitivity (Dqsb_ice/Dtn_ice)

               ! Latent Heat
               zztmp1 = zzblk * rLsub * pCEi(ji,jj)
               qla_ice(ji,jj,jl) = MAX( zztmp1 * (zsq - pq_zu_i(ji,jj)) , 0._wp )   ! #LB: only sublimation (and not condensation) ???
               dqla_ice(ji,jj,jl) =  zztmp1 * dq_sat_dt_ice( zst, pslp(ji,jj) )
               ! if negative `qla_ice` allowed`:
               !dqla_ice(ji,jj,jl) = (0.5_wp + SIGN(0.5_wp,qla_ice(ji,jj,jl))) * zztmp1 * dq_sat_dt_ice(zst, pslp(ji,jj))
               !!IF(qla_ice(ji,jj,jl)>0._wp) dqla_ice(ji,jj,jl) = zztmp1*dq_sat_dt_ice(zst, pslp(ji,jj))  ! ==> Qlat sensitivity  (dQlat/dT)
               !                                                                          !#LB: dq_sat_dt_ice() in "sbc_phy.F90"
               !#LB: without this unjustified "condensation sensure":
               !qla_ice( ji,jj,jl) = zztmp1 * (zsq - pq_zu_i(ji,jj))
               !dqla_ice(ji,jj,jl) = zztmp1 * dq_sat_dt_ice(zst, pslp(ji,jj)) ! ==> Qlat sensitivity  (dQlat/dT)


               ! ----------------------------!
               !     III    Total FLUXES     !
               ! ----------------------------!
               ! Downward Non Solar flux
               qns_ice (ji,jj,jl) =     qlw_ice(ji,jj,jl) - qsb_ice(ji,jj,jl) - qla_ice(ji,jj,jl)
               ! Total non solar heat flux sensitivity for ice
               dqns_ice(ji,jj,jl) = - ( z_dqlw + z_dqsb + dqla_ice(ji,jj,jl) ) !#LB: correct signs ????

            END DO !DO jl = 1, jpl


            zsnwprc = fatm_snow(ji,jj)
            zA_b    = at_i_b(ji,jj)
            z1mA_b  = 1._wp - zA_b

            ! Heat content per unit mass (J/kg)
            zcptrain = (      ptair(ji,jj)        - rt0 ) * rcp  * zmsk
            zcptsnw  = ( MIN( ptair(ji,jj), rt0 ) - rt0 ) * rcpi * zmsk
            zcptn    =        sst_s(ji,jj)                * rcp  * zmsk

            ! --- evaporation --- !
            zevap = emp(ji,jj) + fatm_prcp(ji,jj)   ! evaporation over ocean   !#LOLOfixme WTF??? WHY precip here ???? remove and use evap already computed by bulks!
            !$acc loop seq
            DO jl = 1, jpl                        !  Loop over ice categories  !
               evap_ice (ji,jj,jl) =  qla_ice(ji,jj,jl) * z1_rLsub    ! sublimation
               devap_ice(ji,jj,jl) = dqla_ice(ji,jj,jl) * z1_rLsub    ! d(sublimation)/dT
            END DO

            ! --- evaporation minus precipitation --- !
            zsnw = 1._wp - MAX(1._wp-zA_b , 0._wp)**rn_snwblow  ! used to be a call to `ice_var_snwblow()`

            emp_oce(ji,jj) = z1mA_b * zevap - ( fatm_prcp(ji,jj) - zsnwprc ) - zsnwprc * (1._wp - zsnw )
            zsum = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               zsum = zsum + a_i_b(ji,jj,jl)*evap_ice(ji,jj,jl)
            END DO
            emp_ice(ji,jj) = zsum - zsnwprc * zsnw
            !emp_tot(ji,jj) = emp_oce(ji,jj) + emp_ice(ji,jj)

            ! --- heat flux associated with emp --- !
            qemp_oce(ji,jj) = ( - z1mA_b * zevap * zcptn                             & ! evap at sst !#LOLOfixme => use the real deal !!!
               &                + ( fatm_prcp(ji,jj) - zsnwprc )   *   zcptrain      & ! liquid precip at Tair
               &                +   zsnwprc * ( 1._wp - zsnw ) * ( zcptsnw - rLfus ) &  ! solid precip at min(Tair,Tsnow)
               &                      ) * zmsk
            
            qemp_ice(ji,jj) = ( zsnwprc * zsnw * ( zcptsnw - rLfus ) ) * zmsk   ! solid precip (only)

            ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
            qprec_ice(ji,jj) = rhos * ( zcptsnw - rLfus )

            ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) ---
            !$acc loop seq
            DO jl = 1, jpl
               qevap_ice(ji,jj,jl) = 0._wp ! should be -evap_ice(:,:,jl)*( ( Tice - rt0 ) * rcpi * xmskt(:,:) )
            END DO                         ! But we do not have Tice => consider it at 0degC => evap=0


            ! --- shortwave radiation transmitted thru the surface scattering layer (W/m2) --- !
            IF( nn_qtrice == 0 ) THEN
               ! formulation derived from Grenfell and Maykut (1977), where transmission rate
               !    1) depends on cloudiness
               !    2) is 0 when there is any snow
               !    3) tends to 1 for thin ice
               ztri = 0.18_wp * ( 1._wp - rcloud_fra ) + 0.35_wp * rcloud_fra  ! surface transmission when hi>10cm
               !$acc loop seq
               DO jl = 1, jpl
                  zhi = phi(ji,jj,jl)
                  IF    ( phs(ji,jj,jl) <= 0._wp .AND. zhi <  0.1_wp ) THEN    ! linear decrease from hi=0 to 10cm
                     qtr_ice_top(ji,jj,jl) = qsr_ice(ji,jj,jl) * ( ztri + ( 1._wp - ztri ) * ( 1._wp - zhi * 10._wp ) )
                  ELSEIF( phs(ji,jj,jl) <= 0._wp .AND. zhi >= 0.1_wp ) THEN    ! constant (ztri) when hi>10cm
                     qtr_ice_top(ji,jj,jl) = qsr_ice(ji,jj,jl) * ztri
                  ELSE                                                         ! zero when hs>0
                     qtr_ice_top(ji,jj,jl) = 0._wp
                  ENDIF
               ENDDO
            ELSEIF( nn_qtrice == 1 ) THEN
               ! formulation is derived from the thesis of M. Lebrun (2019).
               !    It represents the best fit using several sets of observations
               !    It comes with snow conductivities adapted to freezing/melting conditions (see icethd_zdf_bl99.F90)
               !$acc loop seq
               DO jl = 1, jpl
                  qtr_ice_top(ji,jj,jl) = 0.3_wp * qsr_ice(ji,jj,jl)
               END DO
               !
            ENDIF !IF( nn_qtrice == 0 )



         END DO
      END DO
      !$acc end parallel loop

      IF( iom_use('qsr_ice') ) THEN
         !$acc update self( qsr_ice )
         CALL iom_put( 'qsr_ice', SUM(   qsr_ice * a_i_b, dim=3 ) )
      ENDIF
      IF( iom_use('qlw_ice') ) THEN
         !$acc update self( qlw_ice )
         CALL iom_put( 'qlw_ice', SUM(   qlw_ice * a_i_b, dim=3 ) )
      ENDIF
      IF( iom_use('qla_ice') ) THEN
         !$acc update self( qla_ice )
         CALL iom_put( 'qla_ice', SUM( - qla_ice * a_i_b, dim=3 ) ) !#LB: sign consistent with what's done for ocean
      ENDIF
      IF( iom_use('qsb_ice') ) THEN
         !$acc update self( qsb_ice )
         CALL iom_put( 'qsb_ice', SUM( - qsb_ice * a_i_b, dim=3 ) ) !#LB: sign consistent with what's done for ocean
      ENDIF
      IF( iom_use('qns_ice') ) THEN
         !$acc update self( qns_ice )
         CALL iom_put( 'qns_ice', SUM(   qns_ice * a_i_b, dim=3 ) )
      ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('blk_ice_2')
   END SUBROUTINE blk_ice_2




   SUBROUTINE blk_ice_qcn( ld_virtual_itd, ptsu, ptb, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_qcn  ***
      !!
      !! ** Purpose :   Compute surface temperature and snow/ice conduction flux
      !!                to force sea ice / snow thermodynamics
      !!                in the case conduction flux is emulated
      !!
      !! ** Method  :   compute surface energy balance assuming neglecting heat storage
      !!                following the 0-layer Semtner (1976) approach
      !!
      !! ** Outputs : - ptsu    : sea-ice / snow surface temperature (K)
      !!              - qcn_ice : surface inner conduction flux (W/m2)
      !!
      !!---------------------------------------------------------------------
      LOGICAL                   , INTENT(in   ) ::   ld_virtual_itd  ! single-category option
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptsu            ! sea ice / snow surface temperature
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   ptb             ! sea ice base temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phs             ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phi             ! sea ice thickness
      !
      INTEGER , PARAMETER ::   nit = 10                  ! number of iterations
      REAL(wp), PARAMETER ::   zepsilon = 0.1_wp         ! characteristic thickness for enhanced conduction
      !
      INTEGER  ::   ji, jj, jl           ! dummy loop indices
      INTEGER  ::   iter                 ! local integer
      REAL(wp) ::   zfac, zfac2, zfac3   ! local scalars
      REAL(wp) ::   zkeff_h, ztsu, ztsu0 !
      REAL(wp) ::   zqc, zqnet           !
      REAL(wp) ::   zhe, zqa0            !
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zgfac   ! enhanced conduction factor
      !!---------------------------------------------------------------------

      ! -------------------------------------!
      !      I   Enhanced conduction factor  !
      ! -------------------------------------!
      ! Emulates the enhancement of conduction by unresolved thin ice (ld_virtual_itd = T)
      ! Fichefet and Morales Maqueda, JGR 1997
      !
      zgfac(:,:,:) = 1._wp

      IF( ld_virtual_itd ) THEN
         !
         zfac  = 1._wp /  ( rcnd_s + rcnd_i )
         zfac2 = EXP(1._wp) * 0.5_wp * zepsilon
         zfac3 = 2._wp / zepsilon
         !
         DO jl = 1, jpl
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  zhe = ( rcnd_s * phi(ji,jj,jl) + rcnd_i * phs(ji,jj,jl) ) * zfac                            ! Effective thickness
                  IF( zhe >=  zfac2 )   zgfac(ji,jj,jl) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( zhe * zfac3 ) ) ) ! Enhanced conduction factor
               END DO
            END DO
         END DO
         !
      ENDIF

      ! -------------------------------------------------------------!
      !      II   Surface temperature and conduction flux            !
      ! -------------------------------------------------------------!
      !
      zfac = rcnd_i * rcnd_s
      !
      DO jl = 1, jpl
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               !
               zkeff_h = zfac * zgfac(ji,jj,jl) / &                                    ! Effective conductivity of the snow-ice system divided by thickness
                  &      ( rcnd_i * phs(ji,jj,jl) + rcnd_s * MAX( 0.01, phi(ji,jj,jl) ) )
               ztsu    = ptsu(ji,jj,jl)                                                ! Store current iteration temperature
               ztsu0   = ptsu(ji,jj,jl)                                                ! Store initial surface temperature
               zqa0    = qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) ! Net initial atmospheric heat flux
               !
               DO iter = 1, nit     ! --- Iterative loop
                  zqc   = zkeff_h * ( ztsu - ptb(ji,jj) )                              ! Conduction heat flux through snow-ice system (>0 downwards)
                  zqnet = zqa0 + dqns_ice(ji,jj,jl) * ( ztsu - ptsu(ji,jj,jl) ) - zqc  ! Surface energy budget
                  ztsu  = ztsu - zqnet / ( dqns_ice(ji,jj,jl) - zkeff_h )              ! Temperature update
               END DO
               !
               ptsu   (ji,jj,jl) = MIN( rt0, ztsu )
               qcn_ice(ji,jj,jl) = zkeff_h * ( ptsu(ji,jj,jl) - ptb(ji,jj) )
               qns_ice(ji,jj,jl) = qns_ice(ji,jj,jl) + dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 )
               qml_ice(ji,jj,jl) = ( qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) - qcn_ice(ji,jj,jl) )  &
                  &   * MAX( 0._wp , SIGN( 1._wp, ptsu(ji,jj,jl) - rt0 ) )

               ! --- Diagnose the heat loss due to changing non-solar flux (as in icethd_zdf_bl99) --- !
               hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) - ( dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 ) ) * a_i_b(ji,jj,jl)

            END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE blk_ice_qcn

   !!======================================================================
END MODULE sbcblk

! LocalWords:  sst
