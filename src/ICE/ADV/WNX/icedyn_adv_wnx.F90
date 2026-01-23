MODULE icedyn_adv_wnx
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_wnx   ***
   !!   sea-ice : advection => WENOX scheme
   !!======================================================================
   !! History :  NANUQ 0.1 !  2025-02 Laurent Brodeau, original code
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice
   USE icevar,   ONLY : ice_var_zapneg
   !USE icedyn_adv_util
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
# if defined _OPENACC
   USE lbclnk_gpu
# else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
# endif
   USE timing         ! Timing

   USE icedyn_adv_wnx_adv

   IMPLICIT NONE
   PRIVATE

   PUBLIC   adv_wnx_init      ! called by icedyn_adv
   PUBLIC   ice_dyn_adv_wnx   ! called by icedyn_adv

   ! Work arrays (that should remain once for all in the memory of the GPU)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sati1

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_wnx( kt, pUu, pVv, pmlbc, ph_i, ph_s, ph_ip,  &
      &                            pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_wnx  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  : WENOX in space and RK3 in time
      !!
      !! Reference:  ...
      !!----------------------------------------------------------------------
      INTEGER                                , INTENT(in   ) ::   kt     ! time step
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pUu    ! ice i-velocity * dy in km^2/s
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pVv    ! ice j-velocity * dx in km^2/s
      INTEGER(1), DIMENSION(jpi,jpj,nn_hls,4), INTENT(in   ) ::   pmlbc  ! masks for solid LBCs
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_i   ! ice thickness
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_s   ! snw thickness
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_ip  ! ice pond thickness
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(inout) ::   pato_i ! open water area
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i   ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s   ! snw volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   psv_i  ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   poa_i  ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i   ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_ip  ! melt pond concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_ip  ! melt pond volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_il  ! melt pond lid volume
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s   ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i   ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i ! ice salt content
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl          ! dummy loop indices
      REAL(wp) ::   zdt, zati2
      !
      CHARACTER(len=15), PARAMETER :: crtnm = 'ice_dyn_adv_wnx'
      CHARACTER(len=1),  PARAMETER :: cgt    = 'T'
      REAL(wp), PARAMETER :: rr_scl_fct = 1.E-6_wp
      REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start(crtnm)
      !$acc data present( pUu, pVv, pmlbc, ph_i, ph_s, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pe_s, pe_i, pszv_i, e1e2t, r1_e1e2t )

      IF( (kt == nit000) .AND. lwp )   WRITE(numout,*) '-- '//crtnm//': WenoX advection scheme'

      zdt = rDt_ice

      ! --- transport --- !
      !      + record at_i before advection (for open water)
      !$acc parallel loop collapse(2) present( sati1 )
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            sati1(ji,jj) = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               sati1(ji,jj) = sati1(ji,jj) + pa_i(ji,jj,jl)
            END DO
            !
         END DO
      END DO
      !$acc end parallel loop

      !IF( ln_icediachk ) THEN
      !   ! diagnostics
      !   zdiag_adv_mass(:,:) =   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
      !      &                  + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow
      !   zdiag_adv_salt(:,:) =   SUM( psv_i(:,:,:) , dim=3 ) * rhoi
      !   zdiag_adv_heat(:,:) = - SUM(SUM( pe_i(:,:,1:nlay_i,:) , dim=4 ), dim=3 ) &
      !      &                  - SUM(SUM( pe_s(:,:,1:nlay_s,:) , dim=4 ), dim=3 )
      !ENDIF

      !! --- ADVECT ---

      ! --- transported fields --- !
      DO jl = 1, jpl

         !== Ice area ==
         CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pa_i(:,:,jl),  lSmesh=.FALSE. )

         !== Ice volume ==
         CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pv_i(:,:,jl),  lSmesh=.FALSE. )

         !== Snow volume ==
         CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pv_s(:,:,jl),  lSmesh=.FALSE. )


         IF( ln_icethd ) THEN

            !== Ice age ==
            CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, poa_i(:,:,jl) )

            !== Salt content ==
            IF( nn_icesal == 4 ) THEN
               DO jk = 1, nlay_i
                  CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pszv_i(:,:,jk,jl),  lSmesh=.FALSE. )
               END DO
            ELSE
               CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, psv_i(:,:,jl),  lSmesh=.FALSE. )
            ENDIF

            !== Ice heat content ==
            DO jk = 1, nlay_i
               CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pe_i(:,:,jk,jl),  lSmesh=.FALSE. )
            END DO

            !== Snow heat content ==
            DO jk = 1, nlay_s
               CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pe_s(:,:,jk,jl),  lSmesh=.FALSE. )
            END DO

            !== melt ponds ==!
            !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   ! concentration
            !   CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pa_ip(:,:,jl),  lSmesh=.FALSE. )
            !   ! volume
            !   CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pv_ip(:,:,jl),  lSmesh=.FALSE. )
            !   ! lid
            !   IF ( ln_pnd_lids ) THEN
            !      CALL wenoX_rk3( kt, cgt, zdt, e1e2t, r1_e1e2t, pUu, pVv, pmlbc, pv_il(:,:,jl),  lSmesh=.FALSE. )
            !   ENDIF
            !ENDIF

         ENDIF !IF( ln_icethd )

      END DO !DO jl = 1, jpl

      !--- derive open water from ice concentration
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !
            zati2 = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zati2 = zati2 + pa_i(ji,jj,jl)
            END DO
            !
            pato_i(ji,jj) = pato_i(ji,jj) - ( zati2 - sati1(ji,jj) ) &
               &                          - ( pUu(ji,jj) - pUu(ji-1,jj) + pVv(ji,jj) - pVv(ji,jj-1) ) * r1_scl_fct * r1_e1e2t(ji,jj) * zdt
         END DO
      END DO
      !$acc end parallel loop
# if ! defined _OPENACC
      CALL lbc_lnk( crtnm, pato_i,cgt,1._wp )
# endif

      !IF( ln_icediachk ) THEN
      !   ! --- diagnostics --- !
      !   diag_adv_mass(:,:) = diag_adv_mass(:,:) + (   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
      !      &                                        + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow &
      !      &                                        - zdiag_adv_mass(:,:) ) * z1_dt
      !   diag_adv_salt(:,:) = diag_adv_salt(:,:) + (   SUM( psv_i(:,:,:) , dim=3 ) * rhoi &
      !      &                                        - zdiag_adv_salt(:,:) ) * z1_dt
      !   diag_adv_heat(:,:) = diag_adv_heat(:,:) + ( - SUM(SUM( pe_i(:,:,1:nlay_i,:) , dim=4 ), dim=3 ) &
      !      &                                        - SUM(SUM( pe_s(:,:,1:nlay_s,:) , dim=4 ), dim=3 ) &
      !      &                                        - zdiag_adv_heat(:,:) ) * z1_dt
      !ENDIF

      IF( .NOT. ln_pureADV2D ) THEN
         ! --- Ensure non-negative fields --- !
         !     Remove negative values (conservation is ensured)
         !     (because advected fields are not perfectly bounded and tiny negative values can occur, e.g. -1.e-20)
         IF( ln_icethd ) THEN
            !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   CALL ice_var_zapneg( zdt, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i )
            !ELSE
            IF( nn_icesal == 4 ) THEN
               CALL ice_var_zapneg( zdt, pv_i, pv_s, psv_i, poa_i, pa_i, pe_s, pe_i,  pszv_i=pszv_i )
            ELSE
               CALL ice_var_zapneg( zdt, pv_i, pv_s, psv_i, poa_i, pa_i, pe_s, pe_i )
            ENDIF
            !ENDIF
         ELSE
            CALL ice_var_zapneg( zdt, pv_i, pv_s,               pa_i )
         ENDIF

      ENDIF !IF( .NOT. ln_pureADV2D )

      !$acc end data
      IF( ln_timing )   CALL timing_stop(crtnm)
      !
   END SUBROUTINE ice_dyn_adv_wnx


   SUBROUTINE adv_wnx_init()
      !!-------------------------------------------------------------------
      !!  * `kwo` is the order of the cuurent WENO scheme: 5, 7, etc
      !!  * `p` is the order of each of the `p` `p`-point stencils used
      !!  *        `w = 2*p -1`
      !!
      !!  => WENO5: `kwo=5`, `p=3` ==> 5th order accuracy, combination of 3 3rd order 3-point stencils.
      !!  => WENO7: `kwo=7`, `p=4` ==> 7th order accuracy, combination of 4 4th order 5-point stencils.
      !!
      !!  Hence the array for linear weights have 4 dimmensions : [Ni x Nj x 2p x w ]
      !!
      !!-------------------------------------------------------------------
      INTEGER, DIMENSION(5) ::   ierr
      INTEGER :: k2p, k_alloc
      CHARACTER(len=12), PARAMETER :: crtnm = 'adv_wnx_init'
      !!-------------------------------------------------------------------
      ierr(:) = 0

      kp_weno = (nn_WNx + 1)/2 ! `p` => order of each of the `p` `p`th order stencil used in the `kp_weno`th order scheme

      !$acc update device ( kp_weno )

      k2p = 2*kp_weno

      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'adv_wnx_init: initialization of WENO advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '      * order of the scheme: `nn_WNx` =', nn_WNx
         WRITE(numout,*) '                    ==>             `p` =', kp_weno
         WRITE(numout,*) '      * allocating linear and optimal weight arrays'
      ENDIF

      ALLOCATE( weno_lw_t_x(jpi,jpj,k2p,nn_WNx), weno_ow_t_x(jpi,jpj,k2p),  &
         &      weno_lw_t_y(jpi,jpj,k2p,nn_WNx), weno_ow_t_y(jpi,jpj,k2p), STAT=ierr(1) )
      IF( ln_damage ) THEN
         ALLOCATE( weno_lw_f_x(jpi,jpj,k2p,nn_WNx), weno_ow_f_x(jpi,jpj,k2p),  &
            &      weno_lw_f_y(jpi,jpj,k2p,nn_WNx), weno_ow_f_y(jpi,jpj,k2p), STAT=ierr(2) )
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*) '                    ==>             done!'
         WRITE(numout,*) '      * now filling them with weights found in file '//TRIM(cn_weno_wght)
      ENDIF

      IF( ln_damage ) THEN
         CALL read_weno_w( jpi, jpj, k2p, nn_WNx, cn_weno_wght, weno_lw_t_x, weno_ow_t_x, weno_lw_t_y, weno_ow_t_y, &
            &                                                   weno_lw_f_x, weno_ow_f_x, weno_lw_f_y, weno_ow_f_y )
         !
         CALL lbc_lnk('adv_wnx_init', weno_lw_t_x,'T',1._wp, weno_lw_t_y,'T',1._wp, weno_lw_f_x,'F',1._wp, weno_lw_f_y,'F',1._wp )
         CALL lbc_lnk('adv_wnx_init', weno_ow_t_x,'T',1._wp, weno_ow_t_y,'T',1._wp, weno_ow_f_x,'F',1._wp, weno_ow_f_y,'F',1._wp )
         !
      ELSE
         CALL read_weno_w( jpi, jpj, k2p, nn_WNx, cn_weno_wght, weno_lw_t_x, weno_ow_t_x, weno_lw_t_y, weno_ow_t_y )
         !
         CALL lbc_lnk('adv_wnx_init', weno_lw_t_x,'T',1._wp, weno_lw_t_y,'T',1._wp )
         CALL lbc_lnk('adv_wnx_init', weno_ow_t_x,'T',1._wp, weno_ow_t_y,'T',1._wp )
      ENDIF

# if defined _OPENACC
      PRINT *, ' * info GPU: adv_wnx_init() => adding linear and optimal weight arrays @T to memory!'
      !$acc enter data copyin( weno_lw_t_x, weno_ow_t_x, weno_lw_t_y, weno_ow_t_y )
      IF( ln_damage ) THEN
         PRINT *, ' * info GPU: adv_wnx_init() => adding linear and optimal weight arrays @F to memory!'
         !$acc enter data copyin( weno_lw_f_x, weno_ow_f_x, weno_lw_f_y, weno_ow_f_y )
      ENDIF
# endif


      IF(lwp) THEN
         WRITE(numout,*) '                    ==>             done!'
         WRITE(numout,'("      WENO",i1," initialization done :D")') nn_WNx
         WRITE(numout,*) ''
      ENDIF

      !! Some debug control:
      IF(lwp) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) ' *** WENO weights at T-points ***'
         CALL print_linear_w(  nn_WNx, weno_lw_t_x(jpi/2,jpj/2,:,:),  'for X, at a given T-point', fmult=6._wp )
         CALL print_linear_w(  nn_WNx, weno_lw_t_y(jpi/2,jpj/2,:,:),  'for Y, at a given T-point', fmult=6._wp )
         CALL print_optimal_w( nn_WNx, weno_ow_t_x(jpi/2,jpj/2,:),    'for X, at a given T-point' )
         CALL print_optimal_w( nn_WNx, weno_ow_t_y(jpi/2,jpj/2,:),    'for Y, at a given T-point' )
         WRITE(numout,*) ''
         IF( ln_damage ) THEN
            WRITE(numout,*) ''
            WRITE(numout,*) ' *** WENO weights at F-points ***'
            CALL print_linear_w(  nn_WNx, weno_lw_f_x(jpi/2,jpj/2,:,:),  'for X, at a given F-point', fmult=6._wp )
            CALL print_linear_w(  nn_WNx, weno_lw_f_y(jpi/2,jpj/2,:,:),  'for Y, at a given F-point', fmult=6._wp )
            CALL print_optimal_w( nn_WNx, weno_ow_f_x(jpi/2,jpj/2,:),    'for X, at a given F-point' )
            CALL print_optimal_w( nn_WNx, weno_ow_f_y(jpi/2,jpj/2,:),    'for Y, at a given F-point' )
            WRITE(numout,*) ''
         ENDIF
      ENDIF

      ALLOCATE( sati1(jpi,jpj), &
         &      zfs1(jpi,jpj), zfs2(jpi,jpj), zfs3(jpi,jpj), zfs4(jpi,jpj),   &
         &      ztrk1(jpi,jpj), ztrk2(jpi,jpj), zoper(jpi,jpj),    STAT = ierr(3) )
      !
      sati1(:,:) = 0._wp
      zfs1(:,:) = 0._wp; zfs2(:,:) = 0._wp; zfs3(:,:) = 0._wp; zfs4(:,:) = 0._wp
      ztrk1(:,:) = 0._wp; ztrk2(:,:) = 0._wp; zoper(:,:) = 0._wp
# if defined _OPENACC
      PRINT *, ' * info GPU: adv_wnx_init() => adding WENO advection workspace arrays to memory'
      PRINT *, '            => sati1, zfs1, zfs2, zfs3, zfs4, ztrk1, ztrk2, zoper'
      !$acc enter data copyin( sati1, zfs1, zfs2, zfs3, zfs4, ztrk1, ztrk2, zoper )
# endif
      
      k_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( crtnm, k_alloc )
      IF( k_alloc > 0 ) CALL ctl_stop('STOP', crtnm//': unable to allocate ice arrays for Prather advection scheme')

   END SUBROUTINE adv_wnx_init

   !!======================================================================
END MODULE icedyn_adv_wnx
