MODULE bdyice
   !!======================================================================
   !!                       ***  MODULE  bdyice  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (SI3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2012-01 (C. Rousset)  add new sea ice model
   !!             4.0  !  2018    (C. Rousset)  SI3 compatibility
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   bdy_ice        : Application of open boundaries to ice
   !!   bdy_ice_frs    : Application of Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice             ! sea-ice: variables
   USE icevar          ! sea-ice: operations
   USE icecor          ! sea-ice: corrections
   USE icectl          ! sea-ice: control prints
   USE phycst          ! physical constant
   USE par_oce         ! ocean parameters
   USE dom_oce         ! ocean space and time domain variables
   USE oss_nnq         ! Bottom boundary condition: ocean fields
   USE bdy             ! Lateral open boundary conditions
   !
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! write to numout file
   USE lib_mpp         ! distributed memory computing
   USE lib_fortran     ! to use key_nosignedzero
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_ice     ! routine called in sbcmod
   PUBLIC   bdy_ice_dyn ! routine called in icedyn_rhg_evp

   PUBLIC   bdy_ice_dmg


   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: bdyice.F90 15368 2021-10-14 08:25:34Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_ice( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ice  ***
      !!
      !! ** Purpose : Apply open boundary conditions for sea ice
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      !
      INTEGER ::   jbdy, ir                             ! BDY set index, rim index
      INTEGER ::   ibeg, iend                           ! length of rim to be treated (rim 0 or rim 1)
      LOGICAL ::   llrim0                               ! indicate if rim 0 is treated
      LOGICAL, DIMENSION(8)  :: llsend1, llrecv1        ! indicate how communications are to be carried out
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('bdy_ice_thd')   ! timing
      !
      CALL ice_var_glo2eqv_gpu(2)
      !

      llsend1(:) = .FALSE.   ;   llrecv1(:) = .FALSE.

      DO ir = 1, 0, -1   ! treat rim 1 before rim 0

         llrim0 = ( ir == 0 )

         DO jbdy = 1, nb_bdy
            !
            SELECT CASE( cn_ice(jbdy) )
            CASE('none')
               CYCLE
            CASE('frs' )
               IF(lwp.AND.(ir==1)) WRITE(numout,*) ' * LOLO [bdy_ice()]: calling `bdy_ice_frs()` for kt:',kt
               CALL bdy_ice_frs( idx_bdy(jbdy), dta_bdy(jbdy), kt, jbdy, llrim0 )
               !lili
            CASE DEFAULT
               CALL ctl_stop( 'bdy_ice : unrecognised `'//TRIM(cn_ice(jbdy))//'` option for open boundaries for ice fields' )
            END SELECT
            !
         END DO
         !
         ! Update bdy points
         IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
         IF( nn_hls == 1 ) THEN
            llsend1(:) = .false.   ;   llrecv1(:) = .false.
         END IF
         DO jbdy = 1, nb_bdy
            IF( cn_ice(jbdy) == 'frs' ) THEN
               llsend1(:) = llsend1(:) .OR. lsend_bdyint(jbdy,1,:,ir)   ! possibly every direction, T points
               llrecv1(:) = llrecv1(:) .OR. lrecv_bdyint(jbdy,1,:,ir)   ! possibly every direction, T points
            END IF
         END DO   ! jbdy
         IF( ANY(llsend1) .OR. ANY(llrecv1) ) THEN   ! if need to send/recv in at least one direction
            ! exchange 3d arrays
            CALL lbc_lnk('bdyice', a_i , 'T', 1._wp, h_i , 'T', 1._wp, h_s , 'T', 1._wp, oa_i, 'T', 1._wp                   &
               &                 , s_i , 'T', 1._wp, t_su, 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp &
               &                 , a_ip, 'T', 1._wp, v_ip, 'T', 1._wp, v_il, 'T', 1._wp                                     &
               &                 , kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
            ! exchange 4d arrays :   third dimension = 1   and then   third dimension = jpk
            CALL lbc_lnk('bdyice', t_s , 'T', 1._wp, e_s , 'T', 1._wp, kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
            CALL lbc_lnk('bdyice', t_i , 'T', 1._wp, e_i , 'T', 1._wp, kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
         END IF
      END DO   ! ir
      !
      CALL ice_cor( kt , 0 )      ! -- In case categories are out of bounds, do a remapping
      !                           !    i.e. inputs have not the same ice thickness distribution (set by rn_himean)
      !                           !         than the regional simulation
      !CALL ice_var_agg(1)
      CALL ice_var_agg_gpu(1)
      !
      ! controls
      IF( ln_icectl )   CALL ice_prt     ( kt, iiceprt, jiceprt, 1, ' - ice thermo bdy - ' )   ! prints
      IF( ln_timing )   CALL timing_stop ('bdy_ice_thd')                                       ! timing
      !
   END SUBROUTINE bdy_ice


   SUBROUTINE bdy_ice_frs( idx, dta, kt, jbdy, llrim0 )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_frs  ***
      !!
      !! ** Purpose : Apply the Flow Relaxation Scheme for sea-ice fields
      !!
      !! Reference : Engedahl H., 1995: Use of the flow relaxation scheme in a three-
      !!             dimensional baroclinic ocean model with realistic topography. Tellus, 365-382.
      !!------------------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx      ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta      ! OBC external data
      INTEGER,         INTENT(in) ::   kt       ! main time-step counter
      INTEGER,         INTENT(in) ::   jbdy     ! BDY set index
      LOGICAL,         INTENT(in) ::   llrim0   ! indicate if rim 0 is treated
      !
      INTEGER  ::   jpbound            ! 0 = incoming ice
      !                                ! 1 = outgoing ice
      INTEGER  ::   ibeg, iend         ! length of rim to be treated (rim 0 or rim 1)
      INTEGER  ::   i_bdy, jgrd        ! dummy loop indices
      INTEGER  ::   ji, jj, jk, jl, ib, jb
      REAL(wp) ::   zmsk, zwgt, zwgt1        ! local scalar
      REAL(wp) ::   ztmelts, zdh
      REAL(wp), POINTER  :: flagu, flagv              ! short cuts
      !!------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('bdy_ice_frs')
      !$acc data present( a_i, h_i, h_s, t_i, t_s, t_su, s_i, sz_i, oa_i, hi_max ) copyin( idx, dta )

      jgrd = 1      ! Everything is at T-points here
      IF( llrim0 ) THEN
         ibeg = 1                       ;   iend = idx%nblenrim0(jgrd)
      ELSE
         ibeg = idx%nblenrim0(jgrd)+1   ;   iend = idx%nblenrim(jgrd)
      END IF
      !
      !$acc loop seq
      DO jl = 1, jpl
         !$acc loop seq
         DO i_bdy = ibeg, iend
            ji    = idx%nbi(i_bdy,jgrd)
            jj    = idx%nbj(i_bdy,jgrd)
            zmsk  = xmskt(ji,jj)
            zwgt  = idx%nbw(i_bdy,jgrd)
            zwgt1 = 1.e0 - idx%nbw(i_bdy,jgrd)
            a_i (ji,jj,  jl) = ( a_i (ji,jj,  jl) * zwgt1 + dta%a_i(i_bdy,jl) * zwgt ) * zmsk  ! Ice  concentration
            h_i (ji,jj,  jl) = ( h_i (ji,jj,  jl) * zwgt1 + dta%h_i(i_bdy,jl) * zwgt ) * zmsk  ! Ice  depth
            h_s (ji,jj,  jl) = ( h_s (ji,jj,  jl) * zwgt1 + dta%h_s(i_bdy,jl) * zwgt ) * zmsk  ! Snow depth
            t_su(ji,jj,  jl) =                              dta%tsu(i_bdy,jl)          * zmsk  ! Surf temperature
            s_i (ji,jj,  jl) =                              dta%s_i(i_bdy,jl)          * zmsk  ! Ice  salinity
            !a_ip(ji,jj,  jl) = ( a_ip(ji,jj,  jl) * zwgt1 + dta%aip(i_bdy,jl) * zwgt ) * zmsk  ! Ice  pond concentration
            !h_ip(ji,jj,  jl) = ( h_ip(ji,jj,  jl) * zwgt1 + dta%hip(i_bdy,jl) * zwgt ) * zmsk  ! Ice  pond depth
            !h_il(ji,jj,  jl) = ( h_il(ji,jj,  jl) * zwgt1 + dta%hil(i_bdy,jl) * zwgt ) * zmsk  ! Ice  pond lid depth
            !
            !$acc loop seq
            DO jk = 1, nlay_i
               sz_i(ji,jj,jk,jl) = s_i(ji,jj,jl)
               t_i (ji,jj,jk,jl) =                          dta%t_i(i_bdy,jl)          * zmsk  ! Ice  temperature
            END DO
            !$acc loop seq
            DO jk = 1, nlay_s
               t_s (ji,jj,jk,jl) =                          dta%t_s(i_bdy,jl)          * zmsk  ! Snow temperature
            END DO
            !lulu

            ! make sure ponds = 0 if no ponds scheme
            !IF( .NOT.ln_pnd ) THEN
            !   a_ip(ji,jj,jl) = 0._wp
            !   h_ip(ji,jj,jl) = 0._wp
            !   h_il(ji,jj,jl) = 0._wp
            !ENDIF

            !IF( .NOT.ln_pnd_lids ) THEN
            !   h_il(ji,jj,jl) = 0._wp
            !ENDIF
            !
            ! -----------------
            ! Pathological case
            ! -----------------
            ! In case a) snow load would be in excess or b) ice is coming into a warmer environment that would lead to
            ! very large transformation from snow to ice (see icethd_dh.F90)

            ! Then, a) transfer the snow excess into the ice (different from icethd_dh)
            zdh = MAX( 0._wp, ( rhos * h_s(ji,jj,jl) + ( rhoi - rho0 ) * h_i(ji,jj,jl) ) * r1_rho0 )
            ! Or, b) transfer all the snow into ice (if incoming ice is likely to melt as it comes into a warmer environment)
            !zdh = MAX( 0._wp, h_s(ji,jj,jl) * rhos / rhoi )

            ! recompute h_i, h_s
            h_i(ji,jj,jl) = MIN( hi_max(jl), h_i(ji,jj,jl) + zdh )
            h_s(ji,jj,jl) = MAX( 0._wp, h_s(ji,jj,jl) - zdh * rhoi / rhos )
            !
         ENDDO
      ENDDO

      !$acc loop seq
      DO jl = 1, jpl
         !$acc loop seq
         DO i_bdy = ibeg, iend

            !IF(lwp .AND. (jl==1) .AND. (i_bdy==ibeg)) WRITE(numout,*) ' * LOLO doing `bdy_ice` stuff with potential bug for kt:',kt

            ji = idx%nbi(i_bdy,jgrd)
            jj = idx%nbj(i_bdy,jgrd)
            flagu => idx%flagu(i_bdy,jgrd)
            flagv => idx%flagv(i_bdy,jgrd)
            ! condition on ice thickness depends on the ice velocity
            ! if velocity is outward (strictly), then ice thickness, volume... must be equal to adjacent values
            jpbound = 0   ;   ib = ji   ;   jb = jj
            !
            IF( flagu ==  1. )   THEN
               IF( ji+1 > jpi )   CYCLE
               IF( u_ice(ji  ,jj  ) < 0. )   jpbound = 1 ; ib = ji+1
            END IF
            IF( flagu == -1. )   THEN
               IF( ji-1 < 1   )   CYCLE
               IF( u_ice(ji-1,jj  ) < 0. )   jpbound = 1 ; ib = ji-1   !LOLO: why not `u > 0.` ??? bug ???
            END IF
            IF( flagv ==  1. )   THEN
               IF( jj+1 > jpj )   CYCLE
               IF( v_ice(ji  ,jj  ) < 0. )   jpbound = 1 ; jb = jj+1
            END IF
            IF( flagv == -1. )   THEN
               IF( jj-1 < 1   )   CYCLE
               IF( v_ice(ji  ,jj-1) < 0. )   jpbound = 1 ; jb = jj-1   !LOLO: why not `v > 0.` ??? bug ???
            END IF
            !
            IF( nn_ice_dta(jbdy) == 0 )   jpbound = 0 ; ib = ji ; jb = jj   ! case ice boundaries = initial conditions
            !                                                               !      do not make state variables dependent on velocity
            !
            IF( a_i(ib,jb,jl) > 0._wp ) THEN   ! there is ice at the boundary
               !
               a_i (ji,jj,  jl) = a_i (ib,jb,  jl)
               h_i (ji,jj,  jl) = h_i (ib,jb,  jl)
               h_s (ji,jj,  jl) = h_s (ib,jb,  jl)
               t_su(ji,jj,  jl) = t_su(ib,jb,  jl)
               s_i (ji,jj,  jl) = s_i (ib,jb,  jl)
               !a_ip(ji,jj,  jl) = a_ip(ib,jb,  jl)
               !h_ip(ji,jj,  jl) = h_ip(ib,jb,  jl)
               !h_il(ji,jj,  jl) = h_il(ib,jb,  jl)
               !
               !$acc loop seq
               DO jk = 1, nlay_i
                  sz_i(ji,jj,jk,jl) = sz_i(ib,jb,jk,jl)
                  t_i (ji,jj,jk,jl) = t_i (ib,jb,jk,jl)
               END DO
               !$acc loop seq
               DO jk = 1, nlay_s
                  t_s (ji,jj,jk,jl) = t_s (ib,jb,jk,jl)
               END DO
               !
               ! ice age
               IF    ( jpbound == 0 ) THEN  ! velocity is inward
                  oa_i(ji,jj,jl) = rice_age(jbdy) * a_i(ji,jj,jl)
               ELSEIF( jpbound == 1 ) THEN  ! velocity is outward
                  oa_i(ji,jj,jl) = oa_i(ib,jb,jl)
               ENDIF
               !
               IF( nn_icesal == 1 ) THEN     ! if constant salinity
                  s_i (ji,jj  ,jl) = rn_icesal
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     sz_i(ji,jj,jk,jl) = rn_icesal
                  END DO
               ENDIF
               !
               ! global fields
               v_i (ji,jj,jl) = h_i(ji,jj,jl) * a_i(ji,jj,jl)                       ! volume ice
               v_s (ji,jj,jl) = h_s(ji,jj,jl) * a_i(ji,jj,jl)                       ! volume snw
               sv_i(ji,jj,jl) = MIN( s_i(ji,jj,jl) , sss_m(ji,jj) ) * v_i(ji,jj,jl) ! salt content
               !$acc loop seq
               DO jk = 1, nlay_s
                  t_s(ji,jj,jk,jl) = MIN( t_s(ji,jj,jk,jl), -0.15_wp + rt0 )           ! Force t_s to be lower than -0.15deg (arbitrary) => likely conservation issue
                  !                                                                    !       otherwise instant melting can occur
                  e_s(ji,jj,jk,jl) = rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl) ) + rLfus )   ! enthalpy in J/m3
                  e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) * r1_nlay_s           ! enthalpy in J/m2
               END DO
               t_su(ji,jj,jl) = MIN( t_su(ji,jj,jl), -0.15_wp + rt0 )                  ! Force t_su to be lower than -0.15deg (arbitrary)
               !$acc loop seq
               DO jk = 1, nlay_i
                  ztmelts          = - rTmlt  * sz_i(ji,jj,jk,jl)             ! Melting temperature in C
                  t_i(ji,jj,jk,jl) = MIN( t_i(ji,jj,jk,jl), (ztmelts-0.15_wp) + rt0 )  ! Force t_i to be lower than melting point (-0.15) => likely conservation issue
                  !
                  e_i(ji,jj,jk,jl) = rhoi * ( rcpi  * ( ztmelts - ( t_i(ji,jj,jk,jl) - rt0 ) )           &   ! enthalpy in J/m3
                     &                      + rLfus * ( 1._wp - ztmelts / ( t_i(ji,jj,jk,jl) - rt0 ) )   &
                     &                      - rcp   *   ztmelts )
                  e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i                            ! enthalpy in J/m2
               END DO
               !
               ! melt ponds
               !v_ip(ji,jj,jl) = h_ip(ji,jj,jl) * a_ip(ji,jj,jl)
               !v_il(ji,jj,jl) = h_il(ji,jj,jl) * a_ip(ji,jj,jl)
               !
            ELSE   ! no ice at the boundary
               !
               a_i(ji,jj,  jl) = 0._wp
               h_i(ji,jj,  jl) = 0._wp
               h_s(ji,jj,  jl) = 0._wp
               oa_i(ji,jj,  jl) = 0._wp
               t_su(ji,jj,  jl) = rt0
               !$acc loop seq
               DO jk = 1, nlay_s
                  t_s (ji,jj,jk,jl) = rt0
               END DO
               !$acc loop seq
               DO jk = 1, nlay_i
                  t_i (ji,jj,jk,jl) = rt0
               END DO

               !a_ip(ji,jj,jl) = 0._wp
               !h_ip(ji,jj,jl) = 0._wp
               !h_il(ji,jj,jl) = 0._wp

               IF( nn_icesal == 1 ) THEN     ! if constant salinity
                  s_i (ji,jj  ,jl) = rn_icesal
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     sz_i(ji,jj,jk,jl) = rn_icesal
                  END DO
               ELSE                          ! if variable salinity
                  s_i (ji,jj,jl)   = rn_simin
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     sz_i(ji,jj,jk,jl) = rn_simin
                  END DO
               ENDIF
               !
               ! global fields
               v_i (ji,jj,  jl) = 0._wp
               v_s (ji,jj,  jl) = 0._wp
               sv_i(ji,jj,  jl) = 0._wp
               !$acc loop seq
               DO jk = 1, nlay_s
                  e_s(ji,jj,jk,jl) = 0._wp
               END DO
               !$acc loop seq
               DO jk = 1, nlay_i
                  e_i(ji,jj,jk,jl) = 0._wp
               END DO
               !v_ip(ji,jj,  jl) = 0._wp
               !v_il(ji,jj,  jl) = 0._wp

            ENDIF

         END DO
         !
      END DO ! jl

      !$acc end data
      IF( ln_timing )   CALL timing_stop('bdy_ice_frs')
      !
   END SUBROUTINE bdy_ice_frs


   SUBROUTINE bdy_ice_dyn( cd_type, pvel,  l_FcVel )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_dyn  ***
      !!
      !! ** Purpose : Apply dynamics boundary conditions for sea-ice.
      !!
      !! ** Method :  if this adjacent grid point is not ice free, then u_ice and v_ice take its value
      !!              if                          is     ice free, then u_ice and v_ice are unchanged by BDY
      !!                                                           they keep values calculated in rheology
      !!
      !!------------------------------------------------------------------------------
      CHARACTER(len=1),             INTENT(in)    ::   cd_type   ! nature of velocity grid-points (U or V)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pvel      ! ice velocity component [m/s]
      LOGICAL,     OPTIONAL,        INTENT(in)    ::   l_FcVel   ! true if dealing with a F-centric velocitie component (u@V or v@U)
      !
      INTEGER               ::   i_bdy, jgrd       ! dummy loop indices
      INTEGER               ::   ji, jj            ! local scalar
      INTEGER               ::   jbdy, ir          ! BDY set index, rim index
      INTEGER               ::   ibeg, iend        ! length of rim to be treated (rim 0 or rim 1)
      INTEGER               ::   ixU, ixV          ! indices to access the relevant ocean velocity component (V_oce(:,:,1->4))
      INTEGER, DIMENSION(3) ::   idir3
      REAL(wp)              ::   zmsk1, zmsk2, zflag
      LOGICAL               :: lFcV
      LOGICAL, DIMENSION(8) :: llsend2, llrecv2, llsend3, llrecv3  ! indicate how communications are to be carried out
      !!------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('bdy_ice_dyn')
      !
      lFcV = .FALSE.
      IF(PRESENT(l_FcVel)) lFcV = l_FcVel
      !
      !! Index to access the relevant ocean velocity component (V_oce(:,:,1->4)):
      ixU = 1 ; ixV = 2
      IF( lFcV ) THEN
         !PRINT *, 'LOLO: bdy_ice_dyn for F-centric vel:'//TRIM(cd_type)
         ixU = 4 ; ixV = 3
      ENDIF
      !
      llsend2(:) = .false.   ;   llrecv2(:) = .false.
      llsend3(:) = .false.   ;   llrecv3(:) = .false.
      DO ir = 1, 0, -1
         DO jbdy = 1, nb_bdy
            !
            SELECT CASE( cn_ice(jbdy) )
               !
            CASE('none')
               CYCLE
               !
            CASE('frs')
               !
               IF( nn_ice_dta(jbdy) == 0 ) CYCLE            ! case ice boundaries = initial conditions
               !                                            !      do not change ice velocity (it is only computed by rheology)
               SELECT CASE ( cd_type )
                  !
               CASE ( 'U' )
                  jgrd = 2      ! u velocity
                  IF( ir == 0 ) THEN
                     ibeg = 1                                 ;   iend = idx_bdy(jbdy)%nblenrim0(jgrd)
                  ELSE
                     ibeg = idx_bdy(jbdy)%nblenrim0(jgrd)+1   ;   iend = idx_bdy(jbdy)%nblenrim(jgrd)
                  END IF
                  DO i_bdy = ibeg, iend
                     ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
                     jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
                     zflag = idx_bdy(jbdy)%flagu(i_bdy,jgrd)
                     !     i-1  i   i    |  !        i  i i+1 |  !          i  i i+1 |
                     !      >  ice  >    |  !        o  > ice |  !          o  >  o  |
                     ! => set at u_ice(i-1) !  => set to O       !  => unchanged
                     IF( zflag == -1. .AND. ji > 1 .AND. ji < jpi )   THEN
                        IF    ( vt_i(ji  ,jj) > 0. )   THEN
                           pvel(ji,jj) = pvel(ji-1,jj)
                        ELSEIF( vt_i(ji+1,jj) > 0. )   THEN
                           pvel(ji,jj) = V_oce(ji,jj,ixU)
                        END IF
                     END IF
                     ! |    i  i+1 i+1        !  |  i   i i+1        !  | i  i i+1
                     ! |    >  ice  >         !  | ice  >  o         !  | o  >  o
                     ! => set at u_ice(i+1)   !     => set to O      !     =>  unchanged
                     IF( zflag ==  1. .AND. ji+1 < jpi+1 )   THEN
                        IF    ( vt_i(ji+1,jj) > 0. )   THEN
                           pvel(ji,jj) = pvel(ji+1,jj)
                        ELSEIF( vt_i(ji  ,jj) > 0. )   THEN
                           pvel(ji,jj) = V_oce(ji,jj,ixU)
                        END IF
                     END IF
                     !
                     IF( zflag ==  0. )   pvel(ji,jj) = 0._wp   ! u_ice = 0  if north/south bdy
                     !
                  END DO
                  !
               CASE ( 'V' )
                  jgrd = 3      ! v velocity
                  IF( ir == 0 ) THEN
                     ibeg = 1                                 ;   iend = idx_bdy(jbdy)%nblenrim0(jgrd)
                  ELSE
                     ibeg = idx_bdy(jbdy)%nblenrim0(jgrd)+1   ;   iend = idx_bdy(jbdy)%nblenrim(jgrd)
                  END IF
                  DO i_bdy = ibeg, iend
                     ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
                     jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
                     zflag = idx_bdy(jbdy)%flagv(i_bdy,jgrd)
                     !                         !      ice   (jj+1)       !       o    (jj+1)
                     !       ^    (jj  )       !       ^    (jj  )       !       ^    (jj  )
                     !      ice   (jj  )       !       o    (jj  )       !       o    (jj  )
                     !       ^    (jj-1)       !                         !
                     ! => set to v_ice(jj-1)   !  =>   set to 0          !   => unchanged
                     IF( zflag == -1. .AND. jj > 1 .AND. jj < jpj )   THEN
                        IF    ( vt_i(ji,jj  ) > 0. )   THEN
                           pvel(ji,jj) = pvel(ji,jj-1)
                        ELSEIF( vt_i(ji,jj+1) > 0. )   THEN
                           pvel(ji,jj) = V_oce(ji,jj,ixV)
                        END IF
                     END IF
                     !       ^    (jj+1)       !                         !
                     !      ice   (jj+1)       !       o    (jj+1)       !       o    (jj+1)
                     !       ^    (jj  )       !       ^    (jj  )       !       ^    (jj  )
                     !   ________________      !  ____ice___(jj  )_      !  _____o____(jj  )
                     ! => set to v_ice(jj+1)   !    => set to 0          !    => unchanged
                     IF( zflag ==  1. .AND. jj < jpj )   THEN
                        IF    ( vt_i(ji,jj+1) > 0. )   THEN
                           pvel(ji,jj) = pvel(ji,jj+1)
                        ELSEIF( vt_i(ji,jj  ) > 0. )   THEN
                           pvel(ji,jj) = V_oce(ji,jj,ixV)
                        END IF
                     END IF
                     !
                     IF( zflag ==  0. )   pvel(ji,jj) = 0._wp   ! v_ice = 0  if west/east bdy
                     !
                  END DO
                  !
               END SELECT !SELECT CASE ( cd_type )
               !
            CASE DEFAULT
               CALL ctl_stop( 'bdy_ice_dyn : unrecognised `'//TRIM(cn_ice(jbdy))//'` option for open boundaries for ice fields' )
            END SELECT !SELECT CASE( cn_ice(jbdy) )
            !
         END DO !DO jbdy = 1, nb_bdy
         !
         SELECT CASE ( cd_type )
         CASE ( 'U' )
            IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
            IF( nn_hls == 1 ) THEN
               llsend2(:) = .false.   ;   llrecv2(:) = .false.
            END IF
            DO jbdy = 1, nb_bdy
               IF( cn_ice(jbdy) == 'frs' .AND. nn_ice_dta(jbdy) /= 0 ) THEN
                  llsend2(  :  ) = llsend2(  :  ) .OR. lsend_bdyint(jbdy,2,  :  ,ir)   ! possibly every direction, U points
                  idir3 = (/ jpwe, jpsw, jpnw /)
                  llsend2(idir3) = llsend2(idir3) .OR. lsend_bdyext(jbdy,2,idir3,ir)   ! nei might search point towards its ea bdy
                  llrecv2(  :  ) = llrecv2(  :  ) .OR. lrecv_bdyint(jbdy,2,  :  ,ir)   ! possibly every direction, U points
                  idir3 = (/ jpea, jpse, jpne /)
                  llrecv2(idir3) = llrecv2(idir3) .OR. lrecv_bdyext(jbdy,2,idir3,ir)   ! might search point towards east bdy
               END IF
            END DO
            IF( ANY(llsend2) .OR. ANY(llrecv2) ) THEN   ! if need to send/recv in at least one direction
               CALL lbc_lnk( 'bdyice', pvel, 'U', -1.0_wp, kfillmode=jpfillnothing ,lsend=llsend2, lrecv=llrecv2 )
            END IF
         CASE ( 'V' )
            IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
            IF( nn_hls == 1 ) THEN
               llsend3(:) = .false.   ;   llrecv3(:) = .false.
            END IF
            DO jbdy = 1, nb_bdy
               IF( cn_ice(jbdy) == 'frs' .AND. nn_ice_dta(jbdy) /= 0 ) THEN
                  llsend3(  :  ) = llsend3(  :  ) .OR. lsend_bdyint(jbdy,3,  :  ,ir)   ! possibly every direction, V points
                  idir3 = (/ jpso, jpsw, jpse /)
                  llsend3(idir3) = llsend3(idir3) .OR. lsend_bdyext(jbdy,3,idir3,ir)   ! nei might search point towards its no bdy
                  llrecv3(  :  ) = llrecv3(  :  ) .OR. lrecv_bdyint(jbdy,3,  :  ,ir)   ! possibly every direction, V points
                  idir3 = (/ jpno, jpnw, jpne /)
                  llrecv3(idir3) = llrecv3(idir3) .OR. lrecv_bdyext(jbdy,3,idir3,ir)   ! might search point towards north bdy
               END IF
            END DO
            IF( ANY(llsend3) .OR. ANY(llrecv3) ) THEN   ! if need to send/recv in at least one direction
               CALL lbc_lnk( 'bdyice', pvel, 'V', -1.0_wp, kfillmode=jpfillnothing ,lsend=llsend3, lrecv=llrecv3 )
            END IF
         END SELECT
      END DO   ! DO ir = 1, 0, -1
      !
      IF( ln_timing )   CALL timing_stop('bdy_ice_dyn')
      !
   END SUBROUTINE bdy_ice_dyn









   !lili


   SUBROUTINE bdy_ice_dmg( kt, kts, pV )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ice_dmg  ***
      !!
      !! ** Purpose : Apply open boundary conditions for damage
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      INTEGER, INTENT(in) ::   kts  ! small time-step counter
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in) :: pV        ! the 4 sea-ice velocity components
      !
      INTEGER ::   jbdy, ir                             ! BDY set index, rim index
      INTEGER ::   ibeg, iend                           ! length of rim to be treated (rim 0 or rim 1)
      LOGICAL ::   llrim0                               ! indicate if rim 0 is treated
      LOGICAL, DIMENSION(8)  :: llsend1, llrecv1        ! indicate how communications are to be carried out
      !!----------------------------------------------------------------------
      ! controls
      IF( ln_timing )   CALL timing_start('bdy_ice_dmg')   ! timing

      llsend1(:) = .false.   ;   llrecv1(:) = .false.
      DO ir = 1, 0, -1   ! treat rim 1 before rim 0
         IF( ir == 0 ) THEN
            llrim0 = .TRUE.
         ELSE
            llrim0 = .FALSE.
         END IF
         DO jbdy = 1, nb_bdy
            !
            SELECT CASE( cn_dmg(jbdy) )
            CASE('none')
               IF(lwp.AND.(ir==1).AND.(kts==10)) WRITE(numout,*) ' * LOLO [bdy_ice_dmg()]: calling NOTHING!!! for kt, kts:',kt,kts
               CYCLE
            CASE('frs' )
               IF(lwp.AND.(ir==1).AND.(kts==10)) WRITE(numout,*) ' * LOLO [bdy_ice_dmg()]: calling `bdy_ice_dmg_frs()` for kt, kts:',kt,kts
               CALL bdy_ice_dmg_frs( idx_bdy(jbdy), dta_bdy(jbdy), kt, kts, jbdy, llrim0, pV )
            CASE('neumann' )
               IF(lwp.AND.(ir==1).AND.(kts==10)) WRITE(numout,*) ' * LOLO [bdy_ice_dmg()]: calling `bdy_dmg_neumann()` for kt, kts:',kt,kts
               CALL bdy_dmg_neumann( idx_bdy(jbdy), dta_bdy(jbdy), kt, kts, ir, jbdy )
            CASE DEFAULT
               CALL ctl_stop( 'bdy_ice_dmg : unrecognised `'//TRIM(cn_dmg(jbdy))//'` option for open boundaries for damage field' )
            END SELECT
            !
         END DO
         !
         ! Update bdy points
         IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
         IF( nn_hls == 1 ) THEN
            llsend1(:) = .false.   ;   llrecv1(:) = .false.
         END IF
         DO jbdy = 1, nb_bdy
            IF( cn_dmg(jbdy) == 'frs' ) THEN
               llsend1(:) = llsend1(:) .OR. lsend_bdyint(jbdy,1,:,ir)   ! possibly every direction, T points
               llrecv1(:) = llrecv1(:) .OR. lrecv_bdyint(jbdy,1,:,ir)   ! possibly every direction, T points
            END IF

            !LOLOfixme: what to do for Neumann !?

         END DO   ! jbdy
         IF( ANY(llsend1) .OR. ANY(llrecv1) ) THEN   ! if need to send/recv in at least one direction
            ! exchange 3d arrays
            CALL lbc_lnk('bdy_ice_dmg', dmdt,'T',1._wp, kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
            !LOLO: add dmg @ F !!!
         END IF
      END DO   ! ir
      !
      IF( ln_timing )   CALL timing_stop ('bdy_ice_dmg')                                       ! timing
      !
   END SUBROUTINE bdy_ice_dmg


   SUBROUTINE bdy_ice_dmg_frs( idx, dta, kt, kts, jbdy, llrim0, pV )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_dmg_frs  ***
      !!
      !! ** Purpose : Apply the Flow Relaxation Scheme on sea-ice damage
      !!
      !! Reference : Engedahl H., 1995: Use of the flow relaxation scheme in a three-
      !!             dimensional baroclinic ocean model with realistic topography. Tellus, 365-382.
      !!------------------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx      ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta      ! OBC external data
      INTEGER,         INTENT(in) ::   kt       ! main time-step counter
      INTEGER,         INTENT(in) ::   kts      ! small time-step counter
      INTEGER,         INTENT(in) ::   jbdy     ! BDY set index
      LOGICAL,         INTENT(in) ::   llrim0   ! indicate if rim 0 is treated
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in) :: pV        ! the 4 sea-ice velocity components
      !
      INTEGER  ::   jpbound            ! 0 = incoming ice
      !                                ! 1 = outgoing ice
      INTEGER  ::   ibeg, iend         ! length of rim to be treated (rim 0 or rim 1)
      INTEGER  ::   i_bdy, jgrd        ! dummy loop indices
      INTEGER  ::   ji, jj, jk, ibt, jbt
      REAL(wp) ::   zwgt, zwgt1        ! local scalar
      REAL(wp), POINTER  :: flagu, flagv ! shortcuts
      !!------------------------------------------------------------------------------
      !
      IF(lwp .AND. (kts==10)) WRITE(numout,*) ' * doing `bdy_ice_dmg_frs` for kt, kts',kt,kts

      jgrd = 1      ! Everything is at T-points here
      IF( llrim0 ) THEN
         ibeg = 1
         iend = idx%nblenrim0(jgrd)
      ELSE
         ibeg = idx%nblenrim0(jgrd)+1
         iend = idx%nblenrim(jgrd)
      END IF
      !
      DO i_bdy = ibeg, iend
         ji    = idx%nbi(i_bdy,jgrd)
         jj    = idx%nbj(i_bdy,jgrd)
         zwgt  = idx%nbw(i_bdy,jgrd)
         zwgt1 = 1.e0 - idx%nbw(i_bdy,jgrd)
         !
         !h_s(ji,jj,  jl) = ( h_s(ji,jj,  jl) * zwgt1 + dta%h_s(i_bdy,jl) * zwgt ) * xmskt(ji,jj)  ! Snow depth
         dmdt(ji,jj) = MIN( MAX( dmdt(ji,jj) * zwgt1 + (1._wp - dta%dmg(i_bdy)) * zwgt , r_dmd_min ), 1._wp ) * xmskt(ji,jj)  ! damage, lolo !
         !
      ENDDO

      DO i_bdy = ibeg, iend
         ji = idx%nbi(i_bdy,jgrd)
         jj = idx%nbj(i_bdy,jgrd)
         flagu => idx%flagu(i_bdy,jgrd)
         flagv => idx%flagv(i_bdy,jgrd)
         ! condition on ice damage depends on the ice velocity
         ! if velocity is outward (strictly), then ice damage must be equal to adjacent values
         jpbound = 0   ;   ibt = ji   ;   jbt = jj
         !
         IF( flagu ==  1. )   THEN
            IF( ji+1 > jpi )   CYCLE
            IF( pV(ji  ,jj  ,1) < 0. )   jpbound = 1 ; ibt = ji+1
         END IF
         IF( flagu == -1. )   THEN
            IF( ji-1 < 1   )   CYCLE
            IF( pV(ji-1,jj  ,1) < 0. )   jpbound = 1 ; ibt = ji-1   !LOLO: why not `u > 0.` ??? bug ???
         END IF
         IF( flagv ==  1. )   THEN
            IF( jj+1 > jpj )   CYCLE
            IF( pV(ji  ,jj  ,2) < 0. )   jpbound = 1 ; jbt = jj+1
         END IF
         IF( flagv == -1. )   THEN
            IF( jj-1 < 1   )   CYCLE
            IF( pV(ji  ,jj-1,2) < 0. )   jpbound = 1 ; jbt = jj-1
         END IF
         !
         IF( nn_dmg_dta(jbdy) == 0 )   jpbound = 0 ; ibt = ji ; jbt = jj   ! case ice boundaries = initial conditions
         !                                                               !      do not make state variables dependent on velocity
         !lolo
         IF( at_i(ibt,jbt) > 0._wp ) THEN   ! there is ice at the boundary
            !
            dmdt(ji,jj) = dmdt(ibt,jbt)
            !
         ELSE   ! no ice at the boundary
            !
            dmdt(ji,jj) = 0._wp
            !lolo add dmdf !!!
            !
         ENDIF

      END DO
      !
   END SUBROUTINE bdy_ice_dmg_frs




   SUBROUTINE bdy_dmg_neumann( idx, dta, kt, kts, ir, jbdy )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dmg_neumann  ***
      !!
      !! ** Purpose :
      !!
      !! Reference :
      !!------------------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx      ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta      ! OBC external data
      INTEGER,         INTENT(in) ::   kt       ! main time-step counter
      INTEGER,         INTENT(in) ::   kts      ! small time-step counter
      INTEGER,         INTENT(in) ::   ir, jbdy     ! BDY set index
      !
      INTEGER  ::   jpbound            ! 0 = incoming ice
      !                                ! 1 = outgoing ice
      INTEGER  ::   ibeg, iend         ! length of rim to be treated (rim 0 or rim 1)
      INTEGER  ::   i_bdy, jgrd, jgrdU, jgrdV        ! dummy loop indices
      INTEGER  ::   ji, jj, jk, ibt, jbt
      !INTEGER  ::   ir, jbdy
      REAL(wp) ::   zwgt, zwgt1        ! local scalar
      REAL(wp) ::   zflagu, zflagv
      LOGICAL  :: ldebug=.FALSE.
      !!------------------------------------------------------------------------------
      !

      jgrd  = 1      ! Everything is at T-points here
      jgrd = 1
      jgrdV = 1


      IF(lwp .AND. (kts==10)) WRITE(numout,*) ' * doing `bdy_dmg_neumann` for kt, kts, jbdy,ir',kt,kts,jbdy, ir
      !
      !
      !!IF( nn_dmg_dta(jbdy) == 0 ) CYCLE            ! case ice boundaries = initial conditions
      !!                                            !      do not change ice velocity (it is only computed by rheology)


      IF( ir == 0 ) THEN
         ibeg = 1                                 ;   iend = idx_bdy(jbdy)%nblenrim0(jgrd)
      ELSE
         ibeg = idx_bdy(jbdy)%nblenrim0(jgrd)+1   ;   iend = idx_bdy(jbdy)%nblenrim(jgrd)
      END IF


      DO i_bdy = ibeg, iend


         !! Eastern/Western BDYs:
         !! ~~~~~~~~~~~~~~~~~~~~~~
         ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
         jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
         zflagu = idx_bdy(jbdy)%flagu(i_bdy,jgrd)

         !     i-1  i   i    |  !        i  i i+1 |  !          i  i i+1 |
         !      >  ice  >    |  !        o  > ice |  !          o  >  o  |
         ! => set at u_ice(i-1) !  => set to O       !  => unchanged
         IF( zflagu == -1. .AND. ji > 1 .AND. ji < jpi )   THEN
            IF    ( vt_i(ji  ,jj) > 0. )   THEN
               !pvel(ji,jj) = pvel(ji-1,jj)
               IF(ldebug) PRINT *, '* LOLO: dmdt at ji=',INT(ji-nn_hls,2),' takes value from ji=',INT(ji-nn_hls-1,2)
               dmdt(ji,jj) = dmdt(ji-1,jj)
               dmdf(ji,jj) = dmdt(ji-1,jj) !lolo?
            ELSEIF( vt_i(ji+1,jj) > 0. )   THEN
               !pvel(ji,jj) = V_oce(ji,jj,ixU)
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes `1.` value'
               dmdt(ji,jj) = 1._wp !lolo?
               dmdf(ji,jj) = 1._wp !lolo?
            END IF
         END IF
         ! |    i  i+1 i+1        !  |  i   i i+1        !  | i  i i+1
         ! |    >  ice  >         !  | ice  >  o         !  | o  >  o
         ! => set at u_ice(i+1)   !     => set to O      !     =>  unchanged
         IF( zflagu ==  1. .AND. ji+1 < jpi+1 )   THEN
            IF    ( vt_i(ji+1,jj) > 0. )   THEN
               !pvel(ji,jj) = pvel(ji+1,jj)
               IF(ldebug) PRINT *, '* LOLO: dmdt at ji=',INT(ji-nn_hls,2),' takes value from ji=',INT(ji-nn_hls+1,2)
               dmdt(ji,jj) = dmdt(ji+1,jj)
               dmdf(ji,jj) = dmdt(ji+1,jj)  !lolo?
            ELSEIF( vt_i(ji  ,jj) > 0. )   THEN
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes `1.` value'
               dmdt(ji,jj) = 1._wp !lolo?
               dmdf(ji,jj) = 1._wp !lolo?
            END IF
         END IF
         !
         !IF( zflagu ==  0. )   pvel(ji,jj) = 0._wp   ! u_ice = 0  if north/south bdy


         !! Northern/Southern BDYs:
         !! ~~~~~~~~~~~~~~~~~~~~~~
         ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
         jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
         zflagv = idx_bdy(jbdy)%flagv(i_bdy,jgrd)
         !                         !      ice   (jj+1)       !       o    (jj+1)
         !       ^    (jj  )       !       ^    (jj  )       !       ^    (jj  )
         !      ice   (jj  )       !       o    (jj  )       !       o    (jj  )
         !       ^    (jj-1)       !                         !
         ! => set to v_ice(jj-1)   !  =>   set to 0          !   => unchanged
         IF( zflagv == -1. .AND. jj > 1 .AND. jj < jpj )   THEN
            IF    ( vt_i(ji,jj  ) > 0. )   THEN
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes value from jj=',INT(jj-nn_hls-1,2)
               dmdt(ji,jj) = dmdt(ji,jj-1)
               dmdf(ji,jj) = dmdt(ji,jj-1)
            ELSEIF( vt_i(ji,jj+1) > 0. )   THEN
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes `1.` value'
               dmdt(ji,jj) = 1._wp
               dmdf(ji,jj) = 1._wp
            END IF
         END IF
         !       ^    (jj+1)       !                         !
         !      ice   (jj+1)       !       o    (jj+1)       !       o    (jj+1)
         !       ^    (jj  )       !       ^    (jj  )       !       ^    (jj  )
         !   ________________      !  ____ice___(jj  )_      !  _____o____(jj  )
         ! => set to v_ice(jj+1)   !    => set to 0          !    => unchanged
         IF( zflagv ==  1. .AND. jj < jpj )   THEN
            IF    ( vt_i(ji,jj+1) > 0. )   THEN
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes value from jj=',INT(jj-nn_hls+1,2)
               dmdt(ji,jj) = dmdt(ji,jj+1)
               dmdf(ji,jj) = dmdt(ji,jj+1)
            ELSEIF( vt_i(ji,jj  ) > 0. )   THEN
               IF(ldebug) PRINT *, '* LOLO: dmdt at jj=',INT(jj-nn_hls,2),' takes `1.` value'
               dmdt(ji,jj) = 1._wp
               dmdf(ji,jj) = 1._wp
            END IF
         END IF
         !
         !IF( zflagv ==  0. )   pvel(ji,jj) = 0._wp   ! v_ice = 0  if west/east bdy
         !
      END DO !DO i_bdy = ibeg, iend
      !
      !

   END SUBROUTINE bdy_dmg_neumann


   !!=================================================================================
END MODULE bdyice
