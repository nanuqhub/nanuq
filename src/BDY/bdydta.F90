MODULE bdydta
   !!======================================================================
   !!                       ***  MODULE bdydta  ***
   !! Open boundary data : read the data for the unstructured open boundaries.
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2012-01  (C. Rousset) add ice boundary conditions for sea ice
   !!            4.0  !  2018     (C. Rousset) SI3 compatibility
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    bdy_dta      : read external data along open boundaries from file
   !!    bdy_dta_init : initialise arrays etc for reading of external data
   !!----------------------------------------------------------------------
   !USE sbc_oce !lolo        ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE bdy        ! ocean open boundary conditions

   USE par_ice, ONLY : r1_nlay_i, r1_nlay_s, ln_pnd, ln_pnd_lids

   USE ice            ! sea-ice variables

   USE icevar         ! redistribute ice input into categories
   !
   USE lib_mpp, ONLY: ctl_stop, ctl_nam
   USE fldread        ! read input fields
   USE iom            ! IOM library
   USE in_out_manager ! I/O logical units
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta          ! routine called by step.F90 and dynspg_ts.F90
   PUBLIC   bdy_dta_init     ! routine called by nanuqgcm.F90

   INTEGER , PARAMETER ::   jpbdyfld  = 11    ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_bdya_i = 1
   INTEGER , PARAMETER ::   jp_bdyh_i = 2
   INTEGER , PARAMETER ::   jp_bdyh_s = 3
   INTEGER , PARAMETER ::   jp_bdyt_i = 4
   INTEGER , PARAMETER ::   jp_bdyt_s = 5
   INTEGER , PARAMETER ::   jp_bdytsu = 6
   INTEGER , PARAMETER ::   jp_bdys_i = 7
   INTEGER , PARAMETER ::   jp_bdydmg = 8
   INTEGER , PARAMETER ::   jp_bdyaip = 9
   INTEGER , PARAMETER ::   jp_bdyhip = 10
   INTEGER , PARAMETER ::   jp_bdyhil = 11

   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:,:), TARGET ::   bf   ! structure of input fields (file informations, fields read)

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: bdydta.F90 15368 2021-10-14 08:25:34Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dta( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta  ***
      !!
      !! ** Purpose :   Update external data for open boundary conditions
      !!
      !! ** Method  :   Use fldread.F90
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           ::   kt           ! ocean time-step index
      INTEGER, INTENT(in)           ::   Kmm          ! ocean time level index
      !
      INTEGER ::  jbdy, jfld, jstart, jend, ib, jl    ! dummy loop indices
      INTEGER ::  ii, ij, ik, igrd, ipl               ! local integers
      TYPE(OBC_DATA)         , POINTER ::   dta_alias        ! short cut
      TYPE(FLD), DIMENSION(:), POINTER ::   bf_alias
      !!---------------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_dta')
      !
      ! Initialise data arrays once for all from initial conditions where required
      !---------------------------------------------------------------------------
      IF( kt == nit000 ) THEN

         ! Calculate depth-mean currents
         !-----------------------------

         DO jbdy = 1, nb_bdy
            !
            !IF( nn_dyn2d_dta(jbdy) == 0 ) THEN
            !   IF( dta_bdy(jbdy)%lneed_ssh ) THEN
            !      igrd = 1
            !      DO ib = 1, idx_bdy(jbdy)%nblenrim(igrd)   ! ssh is allocated and used only on the rim
            !         ii = idx_bdy(jbdy)%nbi(ib,igrd)
            !         ij = idx_bdy(jbdy)%nbj(ib,igrd)
            !         dta_bdy(jbdy)%ssh(ib) = ssh(ii,ij,Kmm) * tmask(ii,ij,1)
            !      END DO
            !   ENDIF
            !   IF( ASSOCIATED(dta_bdy(jbdy)%u2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
            !      igrd = 2
            !      DO ib = 1, SIZE(dta_bdy(jbdy)%u2d)      ! u2d is used either over the whole bdy or only on the rim
            !         ii = idx_bdy(jbdy)%nbi(ib,igrd)
            !         ij = idx_bdy(jbdy)%nbj(ib,igrd)
            !         dta_bdy(jbdy)%u2d(ib) = uu_b(ii,ij,Kmm) * umask(ii,ij,1)
            !      END DO
            !   ENDIF
            !   IF( ASSOCIATED(dta_bdy(jbdy)%v2d) ) THEN   ! no SIZE with a unassociated pointer. v2d and u2d can differ on subdomain
            !      igrd = 3
            !      DO ib = 1, SIZE(dta_bdy(jbdy)%v2d)      ! v2d is used either over the whole bdy or only on the rim
            !         ii = idx_bdy(jbdy)%nbi(ib,igrd)
            !         ij = idx_bdy(jbdy)%nbj(ib,igrd)
            !         dta_bdy(jbdy)%v2d(ib) = vv_b(ii,ij,Kmm) * vmask(ii,ij,1)
            !      END DO
            !   ENDIF
            !ENDIF
            !

            !IF( nn_tra_dta(jbdy) == 0 ) THEN
            !   IF( dta_bdy(jbdy)%lneed_tra ) THEN
            !      igrd = 1
            !      DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
            !         DO ik = 1, jpkm1
            !            ii = idx_bdy(jbdy)%nbi(ib,igrd)
            !            ij = idx_bdy(jbdy)%nbj(ib,igrd)
            !            dta_bdy(jbdy)%tem(ib,ik) = ts(ii,ij,ik,jp_tem,Kmm) * tmask(ii,ij,ik)
            !            dta_bdy(jbdy)%sal(ib,ik) = ts(ii,ij,ik,jp_sal,Kmm) * tmask(ii,ij,ik)
            !         END DO
            !      END DO
            !   ENDIF
            !ENDIF

            IF( nn_ice_dta(jbdy) == 0 ) THEN    ! set ice to initial values
               IF( dta_bdy(jbdy)%lneed_ice ) THEN
                  igrd = 1
                  DO jl = 1, jpl
                     DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%a_i(ib,jl) =  a_i(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%h_i(ib,jl) =  h_i(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%h_s(ib,jl) =  h_s(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%t_i(ib,jl) =  SUM(t_i(ii,ij,:,jl)) * r1_nlay_i * tmask(ii,ij,1)
                        dta_bdy(jbdy)%t_s(ib,jl) =  SUM(t_s(ii,ij,:,jl)) * r1_nlay_s * tmask(ii,ij,1)
                        dta_bdy(jbdy)%tsu(ib,jl) =  t_su(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%s_i(ib,jl) =  s_i(ii,ij,jl) * tmask(ii,ij,1)
                        ! melt ponds
                        dta_bdy(jbdy)%aip(ib,jl) =  a_ip(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%hip(ib,jl) =  h_ip(ii,ij,jl) * tmask(ii,ij,1)
                        dta_bdy(jbdy)%hil(ib,jl) =  h_il(ii,ij,jl) * tmask(ii,ij,1)
                     END DO
                  END DO
               ENDIF
            ENDIF ! IF( nn_ice_dta(jbdy) == 0 )

            !! Damage:
            IF( nn_dmg_dta(jbdy) == 0 ) THEN    ! set damage to initial values
               IF( dta_bdy(jbdy)%lneed_dmg ) THEN
                  igrd = 1
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%dmg(ib) = (1._wp - dmdt(ii,ij)) * tmask(ii,ij,1)
                  END DO
               ENDIF
            ENDIF ! IF( nn_dmg_dta(jbdy) == 0 )




         END DO ! jbdy
         !
      ENDIF ! kt == nit000

      ! update external data from files
      !--------------------------------

      DO jbdy = 1, nb_bdy

         dta_alias => dta_bdy(jbdy)
         bf_alias  => bf(:,jbdy)

         ! read/update all bdy data
         ! ------------------------
         ! BDY: use pt_offset=0.5 as applied at the end of the step and fldread is referenced at the middle of the step
         CALL fld_read( kt, 1, bf_alias, pt_offset = 0.5_wp, Kmm = Kmm )
         ! apply some corrections in some specific cases...
         ! --------------------------------------------------

         IF( dta_alias%lneed_ice .AND. idx_bdy(jbdy)%nblen(1) > 0 ) THEN
            ! fill temperature and salinity arrays
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyt_i)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyt_s)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' )   bf_alias(jp_bdytsu)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf_alias(jp_bdys_i)%clrootname) == 'NOT USED' )   bf_alias(jp_bdys_i)%fnow(:,1,:) = rice_sal(jbdy)
            IF( TRIM(bf_alias(jp_bdydmg)%clrootname) == 'NOT USED' )   bf_alias(jp_bdydmg)%fnow(:,1,:) = rice_dmg(jbdy)
            !
            IF( TRIM(bf_alias(jp_bdyaip)%clrootname) == 'NOT USED' )   &               ! rice_apnd is the pond fraction
               &   bf_alias(jp_bdyaip)%fnow(:,1,:) = rice_apnd(jbdy) * bf_alias(jp_bdya_i)%fnow(:,1,:)   ! ( a_ip = rice_apnd*a_i )
            IF( TRIM(bf_alias(jp_bdyhip)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyhip)%fnow(:,1,:) = rice_hpnd(jbdy)
            IF( TRIM(bf_alias(jp_bdyhil)%clrootname) == 'NOT USED' )   bf_alias(jp_bdyhil)%fnow(:,1,:) = rice_hlid(jbdy)

            ! if T_i is read and not T_su, set T_su = T_i
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdytsu)%fnow(:,1,:) = bf_alias(jp_bdyt_i)%fnow(:,1,:)
            ! if T_s is read and not T_su, set T_su = T_s
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdytsu)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdytsu)%fnow(:,1,:) = bf_alias(jp_bdyt_s)%fnow(:,1,:)
            ! if T_i is read and not T_s, set T_s = T_i
            IF( TRIM(bf_alias(jp_bdyt_i)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_s)%fnow(:,1,:) = bf_alias(jp_bdyt_i)%fnow(:,1,:)
            ! if T_su is read and not T_s, set T_s = T_su
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_s)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_s)%fnow(:,1,:) = bf_alias(jp_bdytsu)%fnow(:,1,:)
            ! if T_su is read and not T_i, set T_i = (T_su + T_freeze)/2
            IF( TRIM(bf_alias(jp_bdytsu)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_i)%fnow(:,1,:) = 0.5_wp * ( bf_alias(jp_bdytsu)%fnow(:,1,:) + 271.15 )
            ! if T_s is read and not T_i, set T_i = (T_s + T_freeze)/2
            IF( TRIM(bf_alias(jp_bdyt_s)%clrootname) /= 'NOT USED' .AND. TRIM(bf_alias(jp_bdyt_i)%clrootname) == 'NOT USED' ) &
               &   bf_alias(jp_bdyt_i)%fnow(:,1,:) = 0.5_wp * ( bf_alias(jp_bdyt_s)%fnow(:,1,:) + 271.15 )

            ! make sure ponds = 0 if no ponds scheme
            IF ( .NOT.ln_pnd ) THEN
               bf_alias(jp_bdyaip)%fnow(:,1,:) = 0._wp
               bf_alias(jp_bdyhip)%fnow(:,1,:) = 0._wp
               bf_alias(jp_bdyhil)%fnow(:,1,:) = 0._wp
            ENDIF
            IF ( .NOT.ln_pnd_lids ) THEN
               bf_alias(jp_bdyhil)%fnow(:,1,:) = 0._wp
            ENDIF

            ! convert N-cat fields (input) into jpl-cat (output)
            ipl = SIZE(bf_alias(jp_bdya_i)%fnow, 3)
            IF( ipl /= jpl ) THEN      ! ice: convert N-cat fields (input) into jpl-cat (output)
               CALL ice_var_itd( bf_alias(jp_bdyh_i)%fnow(:,1,:), bf_alias(jp_bdyh_s)%fnow(:,1,:), bf_alias(jp_bdya_i)%fnow(:,1,:), & ! in
                  &              dta_alias%h_i                  , dta_alias%h_s                  , dta_alias%a_i                  , & ! out
                  &              bf_alias(jp_bdyt_i)%fnow(:,1,:), bf_alias(jp_bdyt_s)%fnow(:,1,:), &                                  ! in (optional)
                  &              bf_alias(jp_bdytsu)%fnow(:,1,:), bf_alias(jp_bdys_i)%fnow(:,1,:), &                                  ! in     -
                  &              bf_alias(jp_bdyaip)%fnow(:,1,:), bf_alias(jp_bdyhip)%fnow(:,1,:), bf_alias(jp_bdyhil)%fnow(:,1,:), & ! in     -
                  &              dta_alias%t_i                  , dta_alias%t_s                  , &                                  ! out    -
                  &              dta_alias%tsu                  , dta_alias%s_i                  , &                                  ! out    -
                  &              dta_alias%aip                  , dta_alias%hip                  , dta_alias%hil )                    ! out    -
            ENDIF
         ENDIF
      END DO  ! jbdy

      !
      IF( ln_timing )   CALL timing_stop('bdy_dta')
      !
   END SUBROUTINE bdy_dta


   SUBROUTINE bdy_dta_init
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta_init  ***
      !!
      !! ** Purpose :   Initialise arrays for reading of external data
      !!                for open boundary conditions
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   jbdy, jfld    ! Local integers
      INTEGER ::   ierror, ios     !
      INTEGER(1) :: iread
      !
      INTEGER ::   nbdy_rdstart, nbdy_loc
      CHARACTER(LEN=50)                      ::   cerrmsg       ! error string
      CHARACTER(len=3)                       ::   cl3           !
      CHARACTER(len=100)                     ::   cn_dir        ! Root directory for location of data files
      REAL(wp)                               ::   rn_ice_tem, rn_ice_sal, rn_ice_age, rn_ice_dmg, rn_ice_apnd, rn_ice_hpnd, rn_ice_hlid
      INTEGER                                ::   ipk,ipl       !
      INTEGER                                ::   idvar         ! variable ID
      INTEGER                                ::   indims        ! number of dimensions of the variable
      INTEGER                                ::   iszdim        ! number of dimensions of the variable
      INTEGER, DIMENSION(4)                  ::   i4dimsz       ! size of variable dimensions
      INTEGER                                ::   igrd          ! index for grid type (1,2,3 = T,U,V)
      LOGICAL                                ::   lluld         ! is the variable using the unlimited dimension
      LOGICAL                                ::   llneed        !
      LOGICAL                                ::   llread        !
      LOGICAL                                ::   llfullbdy     !
      TYPE(FLD_N), DIMENSION(1), TARGET  ::   bn_a_i, bn_h_i, bn_h_s, bn_t_i, bn_t_s, bn_tsu, bn_s_i, bn_dmg, bn_aip, bn_hip, bn_hil
      TYPE(FLD_N), DIMENSION(:), POINTER ::   bn_alias                        ! must be an array to be used with fld_fill
      TYPE(FLD  ), DIMENSION(:), POINTER ::   bf_alias
      !
      NAMELIST/nambdy_dta/ cn_dir, &
         & bn_a_i, bn_h_i, bn_h_s, bn_t_i, bn_t_s, bn_tsu, bn_s_i, bn_dmg, bn_aip, bn_hip, bn_hil, &
         & rn_ice_tem, rn_ice_sal, rn_ice_age, rn_ice_dmg, rn_ice_apnd, rn_ice_hpnd, rn_ice_hlid
      !!---------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'bdy_dta_ini : initialization of data at the open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) ''

      ALLOCATE( bf(jpbdyfld,nb_bdy), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'bdy_dta: unable to allocate bf structure' )   ;   RETURN
      ENDIF
      bf(:,:)%clrootname = 'NOT USED'   ! default definition used as a flag in fld_read to do nothing.
      bf(:,:)%lzint      = .FALSE.      ! default definition
      bf(:,:)%ltotvel    = .FALSE.      ! default definition

      ! Read namelists
      ! --------------
      nbdy_rdstart = 1
      DO jbdy = 1, nb_bdy

         WRITE(ctmp1, '(a,i2)') 'BDY number ', jbdy
         WRITE(ctmp2, '(a,i2)') 'block nambdy_dta number ', jbdy

         ! There is only one nambdy_dta block in namelist_ref -> use it for each bdy so we read from the beginning
         READ_NML_REF(numnam,nambdy_dta)

         !   by-pass nambdy_dta reading if no input data used in this bdy
         IF(       ( dta_bdy(jbdy)%lneed_ice .AND. nn_ice_dta(jbdy) == 1 )   &
            & .OR. ( dta_bdy(jbdy)%lneed_dmg .AND. nn_dmg_dta(jbdy) == 1 )   )   THEN
            !
            ! Need to support possibility of reading more than one
            ! nambdy_dta from the namelist_cfg internal file.
            ! Do this by finding the jbdy'th occurence of nambdy_dta in the
            ! character buffer as the starting point.
            !
            nbdy_loc = INDEX( numnam_cfg( nbdy_rdstart: ), 'nambdy_dta' )
            IF( nbdy_loc .GT. 0 ) THEN
               nbdy_rdstart = nbdy_rdstart + nbdy_loc
            ELSE
               WRITE(cerrmsg,'(A,I4,A)') 'Error: entry number ',jbdy,' of nambdy_dta not found'
               ios = -1
               CALL ctl_nam ( ios , cerrmsg )
            ENDIF
            READ(numnam_cfg( MAX( 1, nbdy_rdstart - 2 ): ), nambdy_dta)
            IF(lwm) WRITE( numond, nambdy_dta )
         ENDIF

         ! get the number of ice categories in bdy data file (use a_i information to do this)
         ipl = jpl   ! default definition
         IF( dta_bdy(jbdy)%lneed_ice ) THEN    ! if we need ice bdy data
            IF( nn_ice_dta(jbdy) == 1 ) THEN   ! if we get ice bdy data from netcdf file
               CALL fld_fill(  bf(jp_bdya_i,jbdy:jbdy), bn_a_i, cn_dir, 'bdy_dta', 'a_i'//' '//ctmp1, ctmp2 )   ! use namelist info
               CALL fld_def( bf(jp_bdya_i,jbdy) )
               CALL iom_open( bf(jp_bdya_i,jbdy)%clname, bf(jp_bdya_i,jbdy)%num )
               idvar = iom_varid( bf(jp_bdya_i,jbdy)%num, bf(jp_bdya_i,jbdy)%clvar, kndims=indims, kdimsz=i4dimsz, lduld=lluld )
               IF( indims == 4 .OR. ( indims == 3 .AND. .NOT. lluld ) ) THEN   ;   ipl = i4dimsz(3)   ! xylt or xyl
               ELSE                                                            ;   ipl = 1            ! xy or xyt
               ENDIF
               CALL iom_close( bf(jp_bdya_i,jbdy)%num )
               bf(jp_bdya_i,jbdy)%clrootname = 'NOT USED'   ! reset to default value as this subdomain may not need to read this bdy
            ENDIF
         ENDIF
         IF(lwp) WRITE(numout,*) ' *** LOLO: n. of ice categories deduced from NC ice bdys: ipl=', ipl !lolorm
         
         IF( .NOT.ln_pnd ) THEN
            rn_ice_apnd = 0. ; rn_ice_hpnd = 0. ; rn_ice_hlid = 0.
            CALL ctl_warn( 'rn_ice_apnd & rn_ice_hpnd = 0 & rn_ice_hlid = 0 when no ponds' )
         ENDIF
         IF( .NOT.ln_pnd_lids ) THEN
            rn_ice_hlid = 0.
         ENDIF

         ! temp, salt, age and ponds of incoming ice
         rice_tem (jbdy) = rn_ice_tem
         rice_sal (jbdy) = rn_ice_sal
         rice_age (jbdy) = rn_ice_age
         rice_dmg (jbdy) = rn_ice_dmg
         rice_apnd(jbdy) = rn_ice_apnd
         rice_hpnd(jbdy) = rn_ice_hpnd
         rice_hlid(jbdy) = rn_ice_hlid


         DO jfld = 1, jpbdyfld

            iread = 1
            
            ! =====================
            !          ice
            ! =====================
            IF(  jfld == jp_bdya_i .OR. jfld == jp_bdyh_i .OR. jfld == jp_bdyh_s .OR. &
               & jfld == jp_bdyt_i .OR. jfld == jp_bdyt_s .OR. jfld == jp_bdytsu .OR. &
               & jfld == jp_bdys_i .OR. jfld == jp_bdyaip .OR. &
               & jfld == jp_bdyhip .OR. jfld == jp_bdyhil ) THEN
               igrd = 1                                                    ! T point
               ipk = ipl                                                   ! jpl-cat data
               llneed = dta_bdy(jbdy)%lneed_ice                            ! ice will be needed
               llread = nn_ice_dta(jbdy) == 1                              ! get data from NetCDF file
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
            ENDIF
            IF(  jfld == jp_bdydmg ) THEN
               igrd = 1                                                    ! T point
               ipk =  1             !lolo                                  ! jpl-cat data
               llneed = dta_bdy(jbdy)%lneed_dmg                            ! dmg will be needed
               llread = nn_dmg_dta(jbdy) == 1                              ! get data from NetCDF file
               iszdim = idx_bdy(jbdy)%nblen(igrd)                          ! length of this bdy on this MPI processus
            ENDIF
            !
            !
            IF( jfld == jp_bdya_i ) THEN
               cl3 = 'a_i'
               bf_alias => bf(jp_bdya_i,jbdy:jbdy)                         ! alias for a_i structure of bdy number jbdy
               bn_alias => bn_a_i                                          ! alias for a_i structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyh_i ) THEN
               cl3 = 'h_i'
               bf_alias => bf(jp_bdyh_i,jbdy:jbdy)                         ! alias for h_i structure of bdy number jbdy
               bn_alias => bn_h_i                                          ! alias for h_i structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyh_s ) THEN
               cl3 = 'h_s'
               bf_alias => bf(jp_bdyh_s,jbdy:jbdy)                         ! alias for h_s structure of bdy number jbdy
               bn_alias => bn_h_s                                          ! alias for h_s structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyt_i ) THEN
               cl3 = 't_i'
               bf_alias => bf(jp_bdyt_i,jbdy:jbdy)                         ! alias for t_i structure of bdy number jbdy
               bn_alias => bn_t_i                                          ! alias for t_i structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyt_s ) THEN
               cl3 = 't_s'
               bf_alias => bf(jp_bdyt_s,jbdy:jbdy)                         ! alias for t_s structure of bdy number jbdy
               bn_alias => bn_t_s                                          ! alias for t_s structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdytsu ) THEN
               cl3 = 'tsu'
               bf_alias => bf(jp_bdytsu,jbdy:jbdy)                         ! alias for tsu structure of bdy number jbdy
               bn_alias => bn_tsu                                          ! alias for tsu structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdys_i ) THEN
               cl3 = 's_i'
               bf_alias => bf(jp_bdys_i,jbdy:jbdy)                         ! alias for s_i structure of bdy number jbdy
               bn_alias => bn_s_i                                          ! alias for s_i structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdydmg ) THEN
               cl3 = 'dmg'
               bf_alias => bf(jp_bdydmg,jbdy:jbdy)                         ! alias for dmg structure of bdy number jbdy
               bn_alias => bn_dmg                                          ! alias for dmg structure of nambdy_dta
               IF( TRIM(cn_dmg(jbdy))/='frs' ) iread = 0
            ENDIF
            IF( jfld == jp_bdyaip ) THEN
               cl3 = 'aip'
               bf_alias => bf(jp_bdyaip,jbdy:jbdy)                         ! alias for aip structure of bdy number jbdy
               bn_alias => bn_aip                                          ! alias for aip structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyhip ) THEN
               cl3 = 'hip'
               bf_alias => bf(jp_bdyhip,jbdy:jbdy)                         ! alias for hip structure of bdy number jbdy
               bn_alias => bn_hip                                          ! alias for hip structure of nambdy_dta
            ENDIF
            IF( jfld == jp_bdyhil ) THEN
               cl3 = 'hil'
               bf_alias => bf(jp_bdyhil,jbdy:jbdy)                         ! alias for hil structure of bdy number jbdy
               bn_alias => bn_hil                                          ! alias for hil structure of nambdy_dta
            ENDIF

            IF( llneed .AND. iszdim > 0 .AND. (iread==1) ) THEN            ! dta_bdy(jbdy)%xxx will be needed
               !                                                           !   -> must be associated with an allocated target
               ALLOCATE( bf_alias(1)%fnow( iszdim, 1, ipk ) )              ! allocate the target
               !
               IF( llread  ) THEN                                           ! get data from NetCDF file
                  CALL fld_fill( bf_alias, bn_alias, cn_dir, 'bdy_dta', cl3//' '//ctmp1, ctmp2 )   ! use namelist info
                  IF( bf_alias(1)%ln_tint ) ALLOCATE( bf_alias(1)%fdta( iszdim, 1, ipk, 2 ) )
                  bf_alias(1)%imap    => idx_bdy(jbdy)%nbmap(1:iszdim,igrd)   ! associate the mapping used for this bdy
                  bf_alias(1)%igrd    = igrd                                  ! used only for vertical integration of 3D arrays
                  bf_alias(1)%ibdy    = jbdy                                  !  "    "    "     "          "      "  "    "
                  bf_alias(1)%ltotvel = .TRUE.   !LOLO don't need             ! T if u3d is full velocity
                  bf_alias(1)%lzint   = .FALSE.  !LOLO don't need             ! T if it requires a vertical interpolation
               ENDIF

               ! associate the pointer and get rid of the dimensions with a size equal to 1
               !IF( jfld == jp_bdyssh )        dta_bdy(jbdy)%ssh => bf_alias(1)%fnow(:,1,1)
               !IF( jfld == jp_bdyu2d )        dta_bdy(jbdy)%u2d => bf_alias(1)%fnow(:,1,1)
               !IF( jfld == jp_bdyv2d )        dta_bdy(jbdy)%v2d => bf_alias(1)%fnow(:,1,1)
               !IF( jfld == jp_bdyu3d )        dta_bdy(jbdy)%u3d => bf_alias(1)%fnow(:,1,:)
               !IF( jfld == jp_bdyv3d )        dta_bdy(jbdy)%v3d => bf_alias(1)%fnow(:,1,:)
               !IF( jfld == jp_bdytem )        dta_bdy(jbdy)%tem => bf_alias(1)%fnow(:,1,:)
               !IF( jfld == jp_bdysal )        dta_bdy(jbdy)%sal => bf_alias(1)%fnow(:,1,:)
               IF( jfld == jp_bdydmg )        dta_bdy(jbdy)%dmg => bf_alias(1)%fnow(:,1,1)

               IF( jfld == jp_bdya_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%a_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%a_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyh_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%h_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%h_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyh_s ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%h_s => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%h_s(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyt_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%t_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%t_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyt_s ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%t_s => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%t_s(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdytsu ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%tsu => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%tsu(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdys_i ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%s_i => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%s_i(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyaip ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%aip => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%aip(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyhip ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%hip => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%hip(iszdim,jpl) )
                  ENDIF
               ENDIF
               IF( jfld == jp_bdyhil ) THEN
                  IF( ipk == jpl ) THEN   ;   dta_bdy(jbdy)%hil => bf_alias(1)%fnow(:,1,:)
                  ELSE                    ;   ALLOCATE( dta_bdy(jbdy)%hil(iszdim,jpl) )
                  ENDIF
               ENDIF
               
            ENDIF !IF( llneed .AND. iszdim > 0 .AND. (iread==1) )

         END DO   ! jpbdyfld
         !
      END DO ! jbdy
      !
   END SUBROUTINE bdy_dta_init

   !!==============================================================================
END MODULE bdydta
