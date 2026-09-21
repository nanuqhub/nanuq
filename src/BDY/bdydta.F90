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

   USE par_ice, ONLY : nlay_i, nlay_s, r1_nlay_i, r1_nlay_s, ln_pnd, ln_pnd_lids
   USE ice     , ONLY : a_i, h_i, h_s, t_i, t_s, t_su, sz_i, dmdt  !, a_ip, h_ip, h_il
   USE icevar  , ONLY : ice_var_itd
   !
   USE lib_mpp , ONLY : ctl_stop, ctl_nam
   USE fldread        ! read input fields
   USE iom            ! IOM library
   USE in_out_manager ! I/O logical units
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta          ! routine called by step.F90
   PUBLIC   bdy_dta_init     ! routine called by nanuqgcm.F90

   TYPE :: VFLD_N
      TYPE(FLD_N), DIMENSION(1) :: XX
   END TYPE VFLD_N

   INTEGER , PARAMETER ::   jpbdyfld  = 11    ! maximum number of files to read
   !!
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

   CHARACTER(len=3), DIMENSION(jpbdyfld) :: vnames = (/ 'a_i', 'h_i', 'h_s', 't_i', 't_s', 'tsu', 's_i', 'dmg', 'aip', 'hip', 'hil' /)


   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: bf   ! structure of input fields (file informations, fields read)

   ! Flags used to indicate boundary-data arrays that have been separated from the corresponding input-data arrays and need to be
   ! reset at each time step
   LOGICAL, DIMENSION(:,:), ALLOCATABLE ::   l_bdydta_reset

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! NEMO/OCE 5.1.a, NEMO Consortium (2026)
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
      INTEGER, INTENT(in)              ::   kt                                  ! ocean time-step index
      INTEGER, INTENT(in)              ::   Kmm                                 ! ocean time level index
      !
      INTEGER ::  jbdy, jfld, jstart, jend, ib, jl, jk
      INTEGER ::  ii, ij, ik, igrd, ipl
      REAL(wp) :: zmsk, ztim_k, zsim_k, ztsm_k
      !!---------------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_dta')
      !$acc data present( idx_bdy, dta_bdy )


      ! Initialise data arrays once for all from initial conditions where required
      !---------------------------------------------------------------------------
      IF( kt == nit000 ) THEN

         !$acc loop seq
         DO jbdy = 1, nb_bdy

            IF( nn_ice_dta(jbdy) == 0 ) THEN    ! set BDY ice values to initial state values
               IF( dta_bdy(jbdy)%lneed_ice ) THEN
                  igrd = 1
                  !$acc parallel loop collapse(2) present( idx_bdy(jbdy)%nblen, idx_bdy(jbdy)%nbi, idx_bdy(jbdy)%nbj, dta_bdy(jbdy)%a_i, dta_bdy(jbdy)%h_i )
                  DO jl = 1, jpl
                     DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        zmsk = xmskt(ii,ij)
                        dta_bdy(jbdy)%a_i(ib,jl) =  a_i (ii,ij,jl) * zmsk
                        dta_bdy(jbdy)%h_i(ib,jl) =  h_i (ii,ij,jl) * zmsk
                        dta_bdy(jbdy)%h_s(ib,jl) =  h_s (ii,ij,jl) * zmsk
                        dta_bdy(jbdy)%tsu(ib,jl) =  t_su(ii,ij,jl) * zmsk
                        !
                        ztim_k = 0 ; zsim_k = 0 ; ztsm_k = 0
                        !$acc loop seq
                        DO jk =1, nlay_i
                           ztim_k = ztim_k +  t_i(ii,ij,jk,jl)
                           zsim_k = zsim_k + sz_i(ii,ij,jk,jl)
                        END DO
                        !$acc loop seq
                        DO jk =1, nlay_s
                           ztsm_k = ztsm_k +  t_s(ii,ij,jk,jl)
                        END DO
                        dta_bdy(jbdy)%t_i(ib,jl) =  ztim_k * r1_nlay_i * zmsk
                        dta_bdy(jbdy)%s_i(ib,jl) =  zsim_k * r1_nlay_i * zmsk
                        dta_bdy(jbdy)%t_s(ib,jl) =  ztsm_k * r1_nlay_s * zmsk
                        ! melt ponds
                        !dta_bdy(jbdy)%aip(ib,jl) =  a_ip(ii,ij,jl) * zmsk
                        !dta_bdy(jbdy)%hip(ib,jl) =  h_ip(ii,ij,jl) * zmsk
                        !dta_bdy(jbdy)%hil(ib,jl) =  h_il(ii,ij,jl) * zmsk
                     END DO
                  END DO
                  !$acc end parallel loop
               ENDIF
            ENDIF

            !! Damage:
            IF( nn_dmg_dta(jbdy) == 0 ) THEN    ! set damage to initial values
               IF( dta_bdy(jbdy)%lneed_dmg ) THEN
                  igrd = 1
                  !$acc parallel loop present( dmdt, dta_bdy(jbdy)%dmg )
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%dmg(ib) = (1._wp - dmdt(ii,ij)) * xmskt(ii,ij)
                  END DO
                  !$acc end parallel loop
               ENDIF
            ENDIF ! IF( nn_dmg_dta(jbdy) == 0 )


         END DO ! jbdy
         !
      ENDIF ! kt == nit000

      ! update external data from files
      !--------------------------------
      !$acc loop seq
      DO jbdy = 1, nb_bdy

         ! read/update all bdy data
         ! ------------------------
         ! BDY: use pt_offset=0.5 as applied at the end of the step and fldread is referenced at the middle of the step
         CALL fld_read( kt, bf(:,jbdy), pt_offset = 0.5_wp, Kmm = Kmm )

         ! apply some corrections in some specific cases...
         ! --------------------------------------------------

         IF( dta_bdy(jbdy)%lneed_ice .AND. idx_bdy(jbdy)%nblen(1) > 0 ) THEN
            ! fill temperature and salinity arrays
            IF( TRIM(bf(jp_bdyt_i,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdyt_i,jbdy)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf(jp_bdyt_s,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdyt_s,jbdy)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf(jp_bdytsu,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdytsu,jbdy)%fnow(:,1,:) = rice_tem(jbdy)
            IF( TRIM(bf(jp_bdys_i,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdys_i,jbdy)%fnow(:,1,:) = rice_sal(jbdy)
            IF( TRIM(bf(jp_bdydmg,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdydmg,jbdy)%fnow(:,1,:) = rice_dmg(jbdy)
            !
            !IF( TRIM(bf(jp_bdyaip,jbdy)%clrootname) == 'NOT_USED' )   &               ! rice_apnd is the pond fraction
            !   &   bf(jp_bdyaip,jbdy)%fnow(:,1,:) = rice_apnd(jbdy) * bf(jp_bdya_i,jbdy)%fnow(:,1,:)   ! ( a_ip = rice_apnd*a_i )
            !IF( TRIM(bf(jp_bdyhip,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdyhip,jbdy)%fnow(:,1,:) = rice_hpnd(jbdy)
            !IF( TRIM(bf(jp_bdyhil,jbdy)%clrootname) == 'NOT_USED' )   bf(jp_bdyhil,jbdy)%fnow(:,1,:) = rice_hlid(jbdy)

            ! if T_i is read and not T_su, set T_su = T_i
            IF( TRIM(bf(jp_bdyt_i,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdytsu,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdytsu,jbdy)%fnow(:,1,:) = bf(jp_bdyt_i,jbdy)%fnow(:,1,:)
            ! if T_s is read and not T_su, set T_su = T_s
            IF( TRIM(bf(jp_bdyt_s,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdytsu,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdytsu,jbdy)%fnow(:,1,:) = bf(jp_bdyt_s,jbdy)%fnow(:,1,:)
            ! if T_i is read and not T_s, set T_s = T_i
            IF( TRIM(bf(jp_bdyt_i,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdyt_s,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdyt_s,jbdy)%fnow(:,1,:) = bf(jp_bdyt_i,jbdy)%fnow(:,1,:)
            ! if T_su is read and not T_s, set T_s = T_su
            IF( TRIM(bf(jp_bdytsu,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdyt_s,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdyt_s,jbdy)%fnow(:,1,:) = bf(jp_bdytsu,jbdy)%fnow(:,1,:)
            ! if T_su is read and not T_i, set T_i = (T_su + T_freeze)/2
            IF( TRIM(bf(jp_bdytsu,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdyt_i,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdyt_i,jbdy)%fnow(:,1,:) = 0.5_wp * ( bf(jp_bdytsu,jbdy)%fnow(:,1,:) + 271.15 )
            ! if T_s is read and not T_i, set T_i = (T_s + T_freeze)/2
            IF( TRIM(bf(jp_bdyt_s,jbdy)%clrootname) /= 'NOT_USED' .AND. TRIM(bf(jp_bdyt_i,jbdy)%clrootname) == 'NOT_USED' ) &
               &   bf(jp_bdyt_i,jbdy)%fnow(:,1,:) = 0.5_wp * ( bf(jp_bdyt_s,jbdy)%fnow(:,1,:) + 271.15 )

            ! make sure ponds = 0 if no ponds scheme
            !IF( .NOT.ln_pnd ) THEN
            !   bf(jp_bdyaip,jbdy)%fnow(:,1,:) = 0._wp
            !   bf(jp_bdyhip,jbdy)%fnow(:,1,:) = 0._wp
            !   bf(jp_bdyhil,jbdy)%fnow(:,1,:) = 0._wp
            !ENDIF
            !IF( .NOT.ln_pnd_lids ) THEN
            !   bf(jp_bdyhil,jbdy)%fnow(:,1,:) = 0._wp
            !ENDIF

#if defined _OPENACC || defined _OPENMP
            ! ==> updating read+corrected data into GPU's memory:
#if defined key_verbose
            PRINT *, '*LOLO [bdy_dta()]: updating freshly read `bf(jp_bdyXXX,jbdy)%fnow` onto GPU [1], jbdy, kt=',jbdy,kt
#endif
            !$acc update device( bf(jp_bdya_i,jbdy)%fnow, bf(jp_bdyh_i,jbdy)%fnow, bf(jp_bdyh_s,jbdy)%fnow, bf(jp_bdyt_i,jbdy)%fnow )
            !$acc update device( bf(jp_bdyt_s,jbdy)%fnow, bf(jp_bdytsu,jbdy)%fnow, bf(jp_bdys_i,jbdy)%fnow, bf(jp_bdydmg,jbdy)%fnow )
#endif

            ! convert N-cat fields (input) into jpl-cat (output)
            ipl = SIZE(bf(jp_bdya_i,jbdy)%fnow, 3)
            IF( ipl /= jpl ) THEN      ! ice: convert N-cat fields (input) into jpl-cat (output)
               CALL ice_var_itd( bf(jp_bdyh_i,jbdy)%fnow(:,1,:), bf(jp_bdyh_s,jbdy)%fnow(:,1,:), bf(jp_bdya_i,jbdy)%fnow(:,1,:), & ! in
                  &              dta_bdy(jbdy)%h_i                  , dta_bdy(jbdy)%h_s                  , dta_bdy(jbdy)%a_i   , & ! out
                  &              bf(jp_bdyt_i,jbdy)%fnow(:,1,:), bf(jp_bdyt_s,jbdy)%fnow(:,1,:), &                                 ! in (optional)
                  &              bf(jp_bdytsu,jbdy)%fnow(:,1,:), bf(jp_bdys_i,jbdy)%fnow(:,1,:), &                                 ! in     -
                  &              bf(jp_bdyaip,jbdy)%fnow(:,1,:), bf(jp_bdyhip,jbdy)%fnow(:,1,:), bf(jp_bdyhil,jbdy)%fnow(:,1,:), & ! in     -
                  &              dta_bdy(jbdy)%t_i                  , dta_bdy(jbdy)%t_s                  , &                       ! out    -
                  &              dta_bdy(jbdy)%tsu                  , dta_bdy(jbdy)%s_i                 ) ! , &                    ! out    -
               !&              dta_bdy(jbdy)%aip                  , dta_bdy(jbdy)%hip                  , dta_bdy(jbdy)%hil )       ! out    -
            ENDIF

         ELSE

#if defined _OPENACC || defined _OPENMP
            ! ==> updating read+corrected data into GPU's memory:
#if defined key_verbose
            PRINT *, '*LOLO [bdy_dta()]: updating freshly read `bf(jp_bdyXXX)%fnow` onto GPU [2], kt=',kt
#endif
            !$acc update device( bf(jp_bdya_i,jbdy)%fnow, bf(jp_bdyh_i,jbdy)%fnow, bf(jp_bdyh_s,jbdy)%fnow, bf(jp_bdyt_i,jbdy)%fnow )
            !$acc update device( bf(jp_bdyt_s,jbdy)%fnow, bf(jp_bdytsu,jbdy)%fnow, bf(jp_bdys_i,jbdy)%fnow, bf(jp_bdydmg,jbdy)%fnow )
#endif

         ENDIF !IF( dta_bdy(jbdy)%lneed_ice .AND. idx_bdy(jbdy)%nblen(1) > 0 )

      END DO  ! jbdy

      !LB [GPU]: we cannot have a single `!$acc update device( bf(jp_bdyt_s)%fnow, ...)` block because `ice_var_itd()` works on the GPU!

      !$acc end data
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
      CHARACTER(LEN=50)           ::   cerrmsg       ! error string
      CHARACTER(len=3)            ::   cl3           !
      CHARACTER(len=100)          ::   cn_dir        ! Root directory for location of data files
      REAL(wp)                    ::   rn_ice_tem, rn_ice_sal, rn_ice_age, rn_ice_dmg, rn_ice_apnd, rn_ice_hpnd, rn_ice_hlid
      INTEGER                     ::   ipk,ipl       !
      INTEGER                     ::   idvar         ! variable ID
      INTEGER                     ::   indims        ! number of dimensions of the variable
      INTEGER                     ::   iszdim        ! number of dimensions of the variable
      INTEGER, DIMENSION(4)       ::   i4dimsz       ! size of variable dimensions
      INTEGER                     ::   igrd          ! index for grid type (1,2,3 = T,U,V)
      LOGICAL                     ::   lluld         ! is the variable using the unlimited dimension
      LOGICAL                     ::   llneed        !
      LOGICAL                     ::   llread        !
      LOGICAL                     ::   llfullbdy     !
      TYPE(FLD_N), DIMENSION(1)   ::   bn_a_i, bn_h_i, bn_h_s, bn_t_i, bn_t_s, bn_tsu, bn_s_i, bn_dmg, bn_aip, bn_hip, bn_hil
      !
      TYPE(VFLD_N), DIMENSION(jpbdyfld) ::  vbn
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
      bf(:,:)%clrootname = 'NOT_USED'   ! default definition used as a flag in fld_read to do nothing.
      bf(:,:)%lzint      = .FALSE.      ! default definition
      bf(:,:)%ltotvel    = .FALSE.      ! default definition

#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: bdy_dta_init() => adding `bf(:)` derived type array to memory'
      PRINT *, '            => bf'
      !$acc enter data copyin( bf )
      PRINT *, '   => will add `bf(:)%fnow(:,:,:)` arrays 1 by 1...'
#endif

      ! Prepare flags that indicate the presence of boundary-data arrays that have been separated from the corresponding input-data
      ! arrays
      ALLOCATE( l_bdydta_reset(jpbdyfld,nb_bdy), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'bdy_dta: memory-allocation failure' )   ;   RETURN
      ENDIF
      l_bdydta_reset(:,:) = .FALSE.


#if defined _OPENACC || defined _OPENMP
      PRINT *, ''
      PRINT *, ' * info GPU: bdy_dta_init() => adding derived type array `dta_bdy` to memory'
      !$acc enter data copyin(dta_bdy)
#endif

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
            READ( numnam_cfg( MAX( 1, nbdy_rdstart - 2 ): ), nambdy_dta, IOSTAT=ios )
            CALL ctl_nam( ios, 'nambdy_dta (numnam_cfg)', .FALSE.)
            IF(lwm) WRITE( numond, nambdy_dta )
         ENDIF

         ! get the number of ice categories in bdy data file (use a_i information to do this)
         ipl = jpl   ! default definition
         IF( dta_bdy(jbdy)%lneed_ice ) THEN    ! if we need ice bdy data
            IF( nn_ice_dta(jbdy) == 1 ) THEN   ! if we get ice bdy data from netcdf file
               CALL fld_fill(  bf(jp_bdya_i,jbdy:jbdy), bn_a_i, cn_dir, 'bdy_dta', vnames(jp_bdya_i)//' '//ctmp1, ctmp2 )   ! use namelist info
               CALL fld_def( bf(jp_bdya_i,jbdy) )
               CALL iom_open( bf(jp_bdya_i,jbdy)%clname, bf(jp_bdya_i,jbdy)%num )
               idvar = iom_varid( 'bdy_dta_init', bf(jp_bdya_i,jbdy)%num, bf(jp_bdya_i,jbdy)%clvar, kndims=indims, kdimsz=i4dimsz, lduld=lluld )
               IF( indims == 4 .OR. ( indims == 3 .AND. .NOT. lluld ) ) THEN
                  ipl = i4dimsz(3) ! xylt or xyl
               ELSE
                  ipl = 1 ! xy or xyt
               ENDIF
               CALL iom_close( bf(jp_bdya_i,jbdy)%num )
               bf(jp_bdya_i,jbdy)%clrootname = 'NOT_USED'   ! reset to default value as this subdomain may not need to read this bdy
            ENDIF
         ENDIF
#if defined key_verbose
         IF(lwp) WRITE(numout,*) ' *** LOLO: n. of ice categories deduced from NC ice bdys: ipl=', ipl !lolorm
#endif

         IF( .NOT.ln_pnd ) THEN
            rn_ice_apnd = 0. ; rn_ice_hpnd = 0. ; rn_ice_hlid = 0.
            CALL ctl_warn( 'rn_ice_apnd & rn_ice_hpnd = 0 & rn_ice_hlid = 0 when no ponds' )
         ENDIF
         IF( .NOT.ln_pnd_lids ) THEN
            rn_ice_hlid = 0.
         ENDIF

         vbn( 1)%XX = bn_a_i
         vbn( 2)%XX = bn_h_i
         vbn( 3)%XX = bn_h_s
         vbn( 4)%XX = bn_t_i
         vbn( 5)%XX = bn_t_s
         vbn( 6)%XX = bn_tsu
         vbn( 7)%XX = bn_s_i
         vbn( 8)%XX = bn_dmg
         vbn( 9)%XX = bn_aip
         vbn(10)%XX = bn_hip
         vbn(11)%XX = bn_hil

         ! temp, salt, age and ponds of incoming ice
         rice_tem (jbdy) = rn_ice_tem
         rice_sal (jbdy) = rn_ice_sal
         rice_age (jbdy) = rn_ice_age
         rice_dmg (jbdy) = rn_ice_dmg
         !rice_apnd(jbdy) = rn_ice_apnd
         !rice_hpnd(jbdy) = rn_ice_hpnd
         !rice_hlid(jbdy) = rn_ice_hlid


         DO jfld = 1, jpbdyfld

            iread = 1

            ! =====================
            !          ice
            ! =====================
            IF(  jfld == jp_bdya_i .OR. jfld == jp_bdyh_i .OR. jfld == jp_bdyh_s .OR. &
               & jfld == jp_bdyt_i .OR. jfld == jp_bdyt_s .OR. jfld == jp_bdytsu .OR. &
               & jfld == jp_bdys_i .OR. jfld == jp_bdyaip .OR. jfld == jp_bdyhip .OR. jfld == jp_bdyhil ) THEN
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

            cl3 = vnames(jfld)

            IF( jfld == jp_bdydmg ) THEN
               IF( TRIM(cn_dmg(jbdy))/='frs' ) iread = 0
            ENDIF

            IF( llneed .AND. iszdim > 0 .AND. (iread==1) ) THEN            ! dta_bdy(jbdy)%xxx will be needed
               ALLOCATE( bf(jfld,jbdy)%fnow( iszdim, 1, ipk ) )
#if defined _OPENACC || defined _OPENMP
               PRINT *, '            => bf(jfld,jbdy)%fnow for jfld,jbdy =', jfld,jbdy
               !$acc enter data copyin( bf(jfld,jbdy)%fnow )
#endif
               !
               IF( llread ) THEN                                           ! get data from NetCDF file
                  CALL fld_fill( bf(jfld,jbdy:jbdy), vbn(jfld)%XX, cn_dir, 'bdy_dta', cl3//' '//ctmp1, ctmp2 )   ! use namelist info
                  IF( bf(jfld,jbdy)%ln_tint ) ALLOCATE( bf(jfld,jbdy)%fdta( iszdim, 1, ipk, 2 ) )
                  ALLOCATE( bf(jfld,jbdy)%imap(iszdim) )
                  bf(jfld,jbdy)%imap    = idx_bdy(jbdy)%nbmap(1:iszdim,igrd)   ! associate the mapping used for this bdy
                  bf(jfld,jbdy)%igrd    = igrd                                  ! used only for vertical integration of 3D arrays
                  bf(jfld,jbdy)%ibdy    = jbdy                                  !  "    "    "     "          "      "  "    "
                  bf(jfld,jbdy)%ltotvel = .TRUE.   !LOLO don't need             ! T if u3d is full velocity
                  bf(jfld,jbdy)%lzint   = .FALSE.  !LOLO don't need             ! T if it requires a vertical interpolation
               ENDIF

               ! associate the pointer and get rid of the dimensions with a size equal to 1

               IF( jfld == jp_bdydmg ) THEN
#if defined key_verbose
                  PRINT *, '*LOLO [bdy_dta_init()] => allocating and filling `dta_bdy(jbdy)%dmg`, jbdy =',jbdy
#endif
                  ALLOCATE( dta_bdy(jbdy)%dmg(iszdim) )
                  dta_bdy(jbdy)%dmg = bf(jfld,jbdy)%fnow(:,1,1)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%dmg(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%dmg)
#endif
               ENDIF

               IF( jfld == jp_bdya_i ) THEN
                  ALLOCATE( dta_bdy(jbdy)%a_i(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%a_i = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%a_i(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%a_i)
#endif
               ENDIF

               IF( jfld == jp_bdyh_i ) THEN
                  ALLOCATE( dta_bdy(jbdy)%h_i(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%h_i = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%h_i(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%h_i)
#endif
               ENDIF

               IF( jfld == jp_bdyh_s ) THEN
                  ALLOCATE( dta_bdy(jbdy)%h_s(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%h_s = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%h_s(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%h_s)
#endif
               ENDIF

               IF( jfld == jp_bdyt_i ) THEN
                  ALLOCATE( dta_bdy(jbdy)%t_i(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%t_i = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%t_i(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%t_i)
#endif
               ENDIF

               IF( jfld == jp_bdyt_s ) THEN
                  ALLOCATE( dta_bdy(jbdy)%t_s(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%t_s = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%t_s(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%t_s)
#endif
               ENDIF

               IF( jfld == jp_bdytsu ) THEN
                  ALLOCATE( dta_bdy(jbdy)%tsu(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%tsu = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%tsu(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%tsu)
#endif
               ENDIF

               IF( jfld == jp_bdys_i ) THEN
                  ALLOCATE( dta_bdy(jbdy)%s_i(iszdim,jpl) )
                  IF( ipk == jpl )  dta_bdy(jbdy)%s_i = bf(jfld,jbdy)%fnow(:,1,:)
#if defined _OPENACC || defined _OPENMP
                  PRINT *, '            => dta_bdy(jbdy)%s_i(:,:), jbdy=',jbdy
                  !$acc enter data copyin(dta_bdy(jbdy)%s_i)
#endif
               ENDIF

               !IF( jfld == jp_bdyaip ) THEN
               !   IF( ipk == jpl ) THEN
               !      dta_bdy(jbdy)%aip = bf(jfld,jbdy)%fnow(:,1,:)
               !   ELSE
               !      ALLOCATE( dta_bdy(jbdy)%aip(iszdim,jpl) )
               !   ENDIF
               !ENDIF
               !IF( jfld == jp_bdyhip ) THEN
               !   IF( ipk == jpl ) THEN
               !      dta_bdy(jbdy)%hip = bf(jfld,jbdy)%fnow(:,1,:)
               !   ELSE
               !      ALLOCATE( dta_bdy(jbdy)%hip(iszdim,jpl) )
               !   ENDIF
               !ENDIF
               !IF( jfld == jp_bdyhil ) THEN
               !   IF( ipk == jpl ) THEN
               !      dta_bdy(jbdy)%hil = bf(jfld,jbdy)%fnow(:,1,:)
               !   ELSE
               !      ALLOCATE( dta_bdy(jbdy)%hil(iszdim,jpl) )
               !   ENDIF
               !ENDIF

            ENDIF !IF( llneed .AND. iszdim > 0 .AND. (iread==1) )

         END DO   ! jpbdyfld
         !
      END DO ! jbdy

#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: bdy_dta_init() => adding `rice_*` 1D arrays to memory'
      PRINT *, '            => rice_tem, rice_sal, rice_dmg, rice_age'
      !$acc enter data copyin( rice_tem, rice_sal, rice_dmg, rice_age )
#endif

   END SUBROUTINE bdy_dta_init

   !!==============================================================================
END MODULE bdydta
