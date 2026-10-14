MODULE icedyn_adv_pra
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_pra   ***
   !!   sea-ice : advection => Prather scheme
   !!======================================================================
   !! History :       !  2008-03  (M. Vancoppenolle) original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_pra : advection of sea ice using Prather scheme
   !!   adv_x_3d, adv_y_3d    : Prather scheme applied in i- and j-direction, resp.
   !!   adv_pra_init    : initialisation of the Prather scheme
   !!   adv_pra_rst     : read/write Prather field in ice restart file, or initialized to zero
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice
   USE icevar,   ONLY : ice_var_zapneg
   USE icedyn_adv_util
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu
#else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#endif
   USE timing         ! Timing

   USE icedyn_adv_pra_adv

   IMPLICIT NONE
   PRIVATE

   PUBLIC   adv_pra_init      ! called by icedyn_adv
   PUBLIC   ice_dyn_adv_pra   ! called by icedyn_adv

   ! Moments for advection
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxice, syice, sxxice, syyice, sxyice   ! ice volume
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxsnw, sysnw, sxxsnw, syysnw, sxysnw   ! snow volume
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxa  , sya  , sxxa  , syya  , sxya     ! ice concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxsz, sysz, sxxsz, syysz, sxysz        ! ice salinity (layers)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxage, syage, sxxage, syyage, sxyage   ! ice age
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxes , syes , sxxes , syyes , sxyes    ! snow layers heat content
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxei , syei , sxxei , syyei , sxyei    ! ice layers heat content
   !REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxap , syap , sxxap , syyap , sxyap    ! melt pond fraction
   !REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxvp , syvp , sxxvp , syyvp , sxyvp    ! melt pond volume
   !REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxvl , syvl , sxxvl , syyvl , sxyvl    ! melt pond lid volume
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg   ! ridged-ice fraction

   ! Work arrays (that should remain once for all in the memory of the GPU)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sati1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   shi_max, shs_max
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sei_max, ses_max, sszi_max

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !SUBROUTINE ice_dyn_adv_pra( kt, pmsk, pUu, pVv, ph_i, ph_s, ph_ip,  &
   !   &                            pato_i, pv_i, pv_s, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i )
   SUBROUTINE ice_dyn_adv_pra( kt, pmsk, pUu, pVv, ph_i, ph_s, pato_i, pv_i, pv_s, poa_i, pa_i, pe_s, pe_i, pszv_i, prdgc )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_pra  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!                but also their quadratic forms. The method preserves
      !!                tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!----------------------------------------------------------------------
      INTEGER                                , INTENT(in   ) ::   kt     ! time step
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pmsk   ! LSM at the relevant point (T)
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pUu    ! ice i-velocity * dy in km^2/s
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pVv    ! ice j-velocity * dx in km^2/s
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_i   ! ice thickness
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_s   ! snw thickness
      !REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(in   ) ::   ph_ip  ! ice pond thickness
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(inout) ::   pato_i ! open water area
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_i   ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_s   ! snw volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   poa_i  ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_i   ! ice concentration
      !REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pa_ip  ! melt pond concentration
      !REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_ip  ! melt pond volume
      !REAL(wp), DIMENSION(jpi,jpj,jpl)       , INTENT(inout) ::   pv_il  ! melt pond lid volume
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s   ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i   ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i ! ice salt content
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(inout) ::   prdgc  ! ridged-ice fraction
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl          ! dummy loop indices
      REAL(wp) ::   zdt, zati2
      REAL(wp) ::   z1_dt, zds
      REAL(wp) ::   zvt                     ! total ice volume (sum over the categories) of a box
      !REAL(wp) ::   zcfl, zcflmax          ! Courant number of a box, and its maximum over the domain
      !CHARACTER(len=8) :: cwcfl             ! `zcflmax` as a string (for the warning)
      !
      CHARACTER(len=15), PARAMETER :: crtnm = 'ice_dyn_adv_pra'
      CHARACTER(len=1),  PARAMETER :: cgt    = 'T'
      REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp ! because using km rather than m...
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start(crtnm)
      !$acc data present(pmsk,pUu,pVv,ph_i,ph_s,pato_i,pv_i,pv_s,poa_i,pa_i,pe_s,pe_i,pszv_i,e1e2t,r1_e1e2t,shi_max,shs_max,sszi_max,sei_max,ses_max)

      IF( (kt == nit000) .AND. lwp )   WRITE(numout,*) '-- '//crtnm//': Prather advection scheme'

      IF( .NOT. ln_pureADV2D ) THEN
         ! --- Record max of the surrounding 9-pts (for call Hbig) --- !
         ! ice + snow thickness
         CALL icemax3D(     ph_i ,       shi_max )
         CALL icemax3D(     ph_s ,       shs_max )
         !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) CALL icemax3D( ph_ip, zhip_max)
         IF( ln_icethd ) THEN
            !! Salt content & Enthalpies
            CALL icemax4D_cnt( nlay_i, pszv_i, pv_i, sszi_max )
            CALL icemax4D_cnt( nlay_i,      pe_i, pv_i,  sei_max )
            CALL icemax4D_cnt( nlay_s,      pe_s, pv_s,  ses_max )

            ! CPU build only (lbc_lnk_gpu is used otherwise):
#if ! ( defined _OPENACC || defined _OPENMP )
            IF( nn_hls <= 2 ) THEN   ! (maybe even `<= 1` would do it actually, but not tested yet...)
               CALL lbc_lnk( crtnm, shi_max,cgt,1._wp, shs_max,cgt,1._wp )   ! 3D
               CALL lbc_lnk( crtnm, sei_max,cgt,1._wp )   ! 4D / nlay_i
               CALL lbc_lnk( crtnm, ses_max,cgt,1._wp )   ! 4D / nlay_s
               CALL lbc_lnk( crtnm, sszi_max,cgt,1._wp )  ! 4D / nlay_i
            ELSE
               CALL lbc_lnk( crtnm, shi_max,cgt,1._wp, shs_max,cgt,1._wp ) ! 3D
            ENDIF
#endif
         ENDIF !IF( ln_icethd )

      ENDIF !IF( .NOT. ln_pureADV2D )

      zdt = rDt_ice
      z1_dt = 1._wp / zdt

      !! --- Courant number check --- !
      !! The Prather scheme moves the fraction `|U|*dt/V` of a box through each of its faces, so the total fraction leaving a box
      !! (through its 2 faces in x, and in y) must not exceed 1, otherwise volumes (and moments) can become negative. Nothing else
      !! checks this, so we compute here the largest such fraction over the domain and send a warning if it is above 1.
      !! `pUu` & `pVv` are transports in km^2/s and `r1_e1e2t*r1_scl_fct` is 1/(box area) in km^-2; `r1_e1e2t` is 0 over land.
      !zcflmax = 0._wp
!!$acc parallel loop collapse(2) reduction(max:zcflmax)
      !DO jj=Njs0, Nje0
      !   DO ji=Nis0, Nie0
      !      ! The x and y passes are done one after the other, so the x and y fractions are tested separately (not added up):
      !      zcfl = MAX( ( MAX( 0._wp, pUu(ji,jj) ) + MAX( 0._wp, -pUu(ji-1,jj) ) ) ,   &
      !         &        ( MAX( 0._wp, pVv(ji,jj) ) + MAX( 0._wp, -pVv(ji,jj-1) ) ) ) * zdt * r1_e1e2t(ji,jj) * r1_scl_fct
      !      zcflmax = MAX( zcflmax, zcfl )
      !   END DO
      !END DO
!!$acc end parallel loop
      !CALL mpp_max( crtnm, zcflmax )   ! max over the whole domain (all the MPI processes)
      !IF( zcflmax > 1._wp ) THEN
      !   WRITE(cwcfl,'(F8.2)') zcflmax
      !   CALL ctl_warn( crtnm//': the Courant number of the Prather advection is above 1 (max = '//TRIM(ADJUSTL(cwcfl))//')', &
      !      &           'The volumes and moments of some boxes can become negative: reduce the time step or the ice velocity' )
      !ENDIF

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

      IF( ln_rdgtrc ) THEN
         ! The ridged-ice fraction `rdgc` is the volume of ridged ice divided by the total volume of ice, so it is an intensive
         ! field weighted by the ice volume: advecting it directly would mix the boxes with a wrong (box area) weight, and would
         ! change it in a divergent flow. We advect instead the volume of ridged ice `Vr = rdgc * vt_i` (a quantity per unit area,
         ! so the standard advection is the right one for it) and we get `rdgc = Vr / vt_i` back at the end of this routine,
         ! with the ice volume after advection. This must be done HERE, before `v_i` is advected below.
         ! (`prdgc` then holds `Vr`, and `sxrdg`, `syrdg`, ... the moments of `Vr`)
         !$acc parallel loop collapse(2) present( prdgc, pv_i )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               zvt = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zvt = zvt + pv_i(ji,jj,jl)
               END DO
               prdgc(ji,jj) = prdgc(ji,jj) * zvt
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      !IF( ln_icediachk ) THEN
      !   ! diagnostics
      !   zdiag_adv_mass(:,:) =   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
      !      &                  + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow
      !   zdiag_adv_heat(:,:) = - SUM(SUM( pe_i(:,:,1:nlay_i,:) , dim=4 ), dim=3 ) &
      !      &                  - SUM(SUM( pe_s(:,:,1:nlay_s,:) , dim=4 ), dim=3 )
      !ENDIF

      !! --- ADVECT ---

      ! --- transported fields --- !
      CALL adv_pra_3d(    kt, zdt, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pa_i,  sxa  , sxxa  , sya  , syya  , sxya   ) !--- ice concentration
      CALL adv_pra_3d(    kt, zdt, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pv_i,  sxice, sxxice, syice, syyice, sxyice ) !--- ice volume
      CALL adv_pra_3d(    kt, zdt, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pv_s,  sxsnw, sxxsnw, sysnw, syysnw, sxysnw ) !--- snow volume

      IF( ln_icethd .AND. ln_age) THEN
         CALL adv_pra_3d( kt, zdt, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  poa_i, sxage, sxxage, syage, syyage, sxyage ) !--- ice age
      ENDIF

      IF( ln_icethd ) THEN
         CALL adv_pra_4d( kt, zdt, nlay_i, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pszv_i, sxsz, sxxsz, sysz, syysz, sxysz ) !--- ice salinity
         CALL adv_pra_4d( kt, zdt, nlay_s,    pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pe_s,   sxes, sxxes, syes, syyes, sxyes ) !--- snow heat content
         CALL adv_pra_4d( kt, zdt, nlay_i,    pUu, pVv, e1e2t, r1_e1e2t, pmsk,  pe_i,   sxei, sxxei, syei, syyei, sxyei ) !--- ice  heat content
      ENDIF

      !== Ridged ice fraction ==
      IF( ln_rdgtrc ) THEN
         CALL adv_pra_2d( kt, zdt, pUu, pVv, e1e2t, r1_e1e2t, pmsk,  prdgc, sxrdg, sxxrdg, syrdg, syyrdg, sxyrdg )
      ENDIF


      !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
      !   CALL adv_x_3d( zdt, pUu, 1._wp, e1e2t, r1_e1e2t, pmsk,  z0ap, sxap, sxxap, syap, syyap, sxyap )    !--- melt pond fraction
      !   CALL adv_y_3d( zdt, pVv, 0._wp, e1e2t, r1_e1e2t, pmsk,  z0ap, sxap, sxxap, syap, syyap, sxyap )
      !   CALL adv_x_3d( zdt, pUu, 1._wp, e1e2t, r1_e1e2t, pmsk,  z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )    !--- melt pond volume
      !   CALL adv_y_3d( zdt, pVv, 0._wp, e1e2t, r1_e1e2t, pmsk,  z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )
      !   IF ( ln_pnd_lids ) THEN
      !      CALL adv_x_3d( zdt, pUu, 1._wp, e1e2t, r1_e1e2t, pmsk,  z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl ) !--- melt pond lid volume
      !      CALL adv_y_3d( zdt, pVv, 0._wp, e1e2t, r1_e1e2t, pmsk,  z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl )
      !   ENDIF
      !ENDIF


      ! --- Lateral boundary conditions --- !
      !        caution: for gradients (sx and sy) the sign changes
#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( crtnm,  pv_i,  sxice,  syice,  pv_s,   sxsnw,   sysnw,  pa_i,    sxa,    sya    )
      CALL lbc_lnk_gpu( crtnm,  sxxice, syyice, sxyice, sxxsnw,  syysnw,  sxysnw, sxxa,   syya,   sxya  )
      !
      IF( ln_icethd ) THEN
         ! Same fields as in the CPU branch below. The GPU build runs on 1 GPU (no MPI), so `lbc_lnk_gpu` only takes care of the
         ! East-West & North-South periodicity (no sign change is needed): without this, the periodic halos of these fields
         ! would be stale at the next time step, and the fluxes across the periodic seam would be wrong.
         ! One call per tracer: the field itself + its 5 moments, with `nlay_i` or `nlay_s` layers (4D arrays):
         CALL lbc_lnk_gpu( crtnm, nlay_i, pszv_i, sxsz, sysz, sxxsz, syysz, sxysz )   ! ice salinity
         CALL lbc_lnk_gpu( crtnm, nlay_i, pe_i,   sxei, syei, sxxei, syyei, sxyei )   ! ice enthalpy
         CALL lbc_lnk_gpu( crtnm, nlay_s, pe_s,   sxes, syes, sxxes, syyes, sxyes )   ! snow enthalpy
         IF( ln_age ) THEN
            ! Ice age (3D arrays): the field itself + its 5 moments:
            CALL lbc_lnk_gpu( crtnm, poa_i, sxage, syage, sxxage, syyage, sxyage )
         ENDIF
      ENDIF
#else
      CALL lbc_lnk( crtnm, pv_i,cgt,1._wp,  sxice,cgt,-1._wp,  syice,cgt,-1._wp,  & ! ice volume
         &               sxxice,cgt,1._wp, syyice,cgt, 1._wp, sxyice,cgt, 1._wp,  &
         &                 pv_s,cgt,1._wp,   sxsnw,cgt,-1._wp,   sysnw,cgt,-1._wp,  & ! snw volume
         &                sxxsnw,cgt,1._wp,  syysnw,cgt, 1._wp,  sxysnw,cgt, 1._wp,  &
         &                 pa_i,cgt,1._wp,    sxa,cgt,-1._wp,    sya,cgt,-1._wp,  & ! ice concentration
         &                 sxxa,cgt,1._wp,   syya,cgt, 1._wp,   sxya,cgt, 1._wp  )
      !
      IF( ln_icethd ) THEN
         CALL lbc_lnk( crtnm,  pszv_i,cgt,1._wp,  sxsz,cgt,-1._wp,  sysz,cgt,-1._wp,   & ! ice salinity 3D+nlay_i
            &                 sxxsz,cgt,1._wp, syysz,cgt, 1._wp, sxysz,cgt, 1._wp )
         CALL lbc_lnk( crtnm,  pe_i,cgt,1._wp,  sxei,cgt,-1._wp,  syei,cgt,-1._wp,     & ! ice enthalpy 3D+nlay_i
            &                 sxxei,cgt,1._wp, syyei,cgt, 1._wp, sxyei,cgt, 1._wp )
         CALL lbc_lnk( crtnm,  pe_s,cgt,1._wp,  sxes,cgt,-1._wp,  syes,cgt,-1._wp,     & ! snw enthalpy 3D+nlay_s
            &                 sxxes,cgt,1._wp, syyes,cgt, 1._wp, sxyes,cgt, 1._wp )
         IF( ln_age ) THEN
            CALL lbc_lnk( crtnm,  poa_i,cgt,1._wp,  sxage,cgt,-1._wp,  syage,cgt,-1._wp,  & ! ice age 3D
               &                 sxxage,cgt,1._wp, syyage,cgt, 1._wp, sxyage,cgt, 1._wp  )
         ENDIF
      ENDIF
      !
      !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
      !   IF( ln_pnd_lids ) THEN
      !      CALL lbc_lnk( crtnm, z0ap,cgt,1._wp, sxap,cgt,-1._wp, syap,cgt,-1._wp  & ! melt pond fraction
      !         &                          , sxxap,cgt,1._wp, syyap,cgt,1._wp, sxyap,cgt,1._wp  &
      !         &                          , z0vp,cgt,1._wp, sxvp,cgt,-1._wp, syvp,cgt,-1._wp  & ! melt pond volume
      !         &                          , sxxvp,cgt,1._wp, syyvp,cgt,1._wp, sxyvp,cgt,1._wp  &
      !         &                          , z0vl,cgt,1._wp, sxvl,cgt,-1._wp, syvl,cgt,-1._wp  & ! melt pond lid volume
      !         &                          , sxxvl,cgt,1._wp, syyvl,cgt,1._wp, sxyvl,cgt,1._wp  )
      !   ELSE
      !      CALL lbc_lnk( crtnm, z0ap,cgt,1._wp, sxap,cgt,-1._wp, syap,cgt,-1._wp  & ! melt pond fraction
      !         &                          , sxxap,cgt,1._wp, syyap,cgt,1._wp, sxyap,cgt,1._wp  &
      !         &                          , z0vp,cgt,1._wp, sxvp,cgt,-1._wp, syvp,cgt,-1._wp  & ! melt pond volume
      !         &                          , sxxvp,cgt,1._wp, syyvp,cgt,1._wp, sxyvp,cgt,1._wp  )
      !   ENDIF
      !ENDIF
#endif


      !--- derive open water from ice concentration
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
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

      IF( ln_rdgtrc ) THEN
         ! Open water and ridged-ice fraction: lateral boundary conditions.
         ! The 5 Prather moments of the ridged-ice fraction must be exchanged as well (as is done for all the other advected
         ! tracers): the advection routines do not update all the halo points, so without this the moments in the halos are stale
         ! and the first pass of the next time step would use them to compute fluxes into the first interior points.
         ! Sign change for the gradients (`sx` and `sy`), as for the other tracers (cf. the `lbc_lnk` calls above).
#if defined _OPENACC || defined _OPENMP
         ! `pato_i`, `prdgc` and its 5 moments in a single call (7 fields):
         CALL lbc_lnk_gpu( crtnm, pato_i, prdgc, sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg )
#else
         CALL lbc_lnk(     crtnm, pato_i,cgt,1._wp, prdgc,cgt,1._wp,                                  &
            &                     sxrdg,cgt,-1._wp, syrdg,cgt,-1._wp,                                 &
            &                     sxxrdg,cgt,1._wp, syyrdg,cgt,1._wp, sxyrdg,cgt,1._wp )
#endif
      ELSE
#if defined _OPENACC || defined _OPENMP
         CALL lbc_lnk_gpu( crtnm, pato_i )
#else
         CALL lbc_lnk(     crtnm, pato_i,cgt,1._wp )
#endif
      ENDIF
      !IF( ln_icediachk ) THEN
      !   ! --- diagnostics --- !
      !   diag_adv_mass(:,:) = diag_adv_mass(:,:) + (   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
      !      &                                        + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow &
      !      &                                        - zdiag_adv_mass(:,:) ) * z1_dt
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
            !   CALL ice_var_zapneg( zdt, pv_i, pv_s, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i )
            !ELSE
            CALL ice_var_zapneg( zdt, pv_i, pv_s, poa_i, pa_i, pe_s, pe_i, pszv_i, prdgc )
            !ENDIF
         ELSE
            CALL ice_var_zapneg( zdt, pv_i, pv_s,        pa_i,                     prdgc )
         ENDIF
         !
         ! --- Make sure ice thickness is not too big --- !
         !     (because ice thickness can be too large where ice concentration is very small)
         IF( ln_icethd ) THEN
            !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   CALL Hbig( zdt, shi_max, shs_max, zhip_max, ssi_max, ses_max, sei_max, pv_i, pv_s, pa_i, pa_ip, pv_ip, pe_s, pe_i )
            !ELSE
            CALL Hbig( zdt, shi_max, shs_max, sszi_max, ses_max, sei_max, pv_i, pv_s, pa_i, pszv_i, pe_s, pe_i )
            !ENDIF
         ELSE
            CALL Hbig( zdt, shi_max, shs_max,                             pv_i, pv_s, pa_i )
         ENDIF

         ! --- Ensure snow load is not too big --- !
         !   ==> included into `Hbig` !!!

      ENDIF !IF( .NOT. ln_pureADV2D )

      IF( ln_rdgtrc ) THEN
         ! Back from the advected volume of ridged ice `Vr` (held by `prdgc` since before the advection) to the ridged-ice fraction:
         ! `rdgc = Vr / vt_i`, with the total ice volume `vt_i` AFTER advection (and after the corrections of `v_i` just above, so
         ! that `rdgc` is consistent with the final `v_i`). The fraction is bounded by [0,1], and set to 0 where there is
         ! (almost) no ice. The whole domain, halos included, is converted: the exchange above filled the halos with `Vr`
         ! (and `v_i`), so they must be converted too, otherwise `rdgc` would be `Vr` in the halos.
         !$acc parallel loop collapse(2) present( prdgc, pv_i )
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               zvt = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zvt = zvt + pv_i(ji,jj,jl)
               END DO
               prdgc(ji,jj) = MERGE( MIN( MAX( prdgc(ji,jj) / MAX( zvt, epsi10 ) , 0._wp ) , 1._wp ) , 0._wp , zvt > epsi10 )
            END DO
         END DO
         !$acc end parallel loop
      ENDIF


      IF( lrst_ice )   CALL adv_pra_rst( 'WRITE', kt )   !* write Prather fields in the restart file !LOLOfixme!

      !$acc end data
      IF( ln_timing )   CALL timing_stop(crtnm)
      !
   END SUBROUTINE ice_dyn_adv_pra


   SUBROUTINE adv_pra_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE adv_pra_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for Prather advection
      !!-------------------------------------------------------------------
      INTEGER, DIMENSION(15) ::  ierr
      INTEGER :: k_alloc
      CHARACTER(len=12), PARAMETER :: crtnm = 'adv_pra_init'
      !!-------------------------------------------------------------------
      ierr(:) = 0
      !                             !* allocate prather fields
      ALLOCATE( sxice(jpi,jpj,jpl) , syice(jpi,jpj,jpl) , sxxice(jpi,jpj,jpl) , syyice(jpi,jpj,jpl) , sxyice(jpi,jpj,jpl) ,   &
         &      sxsnw(jpi,jpj,jpl) , sysnw(jpi,jpj,jpl) , sxxsnw(jpi,jpj,jpl) , syysnw(jpi,jpj,jpl) , sxysnw(jpi,jpj,jpl) ,   &
         &      sxa  (jpi,jpj,jpl) , sya  (jpi,jpj,jpl) , sxxa  (jpi,jpj,jpl) , syya  (jpi,jpj,jpl) , sxya  (jpi,jpj,jpl) ,   &
         &      STAT=ierr(1) )
      sxice(:,:,:) = 0._wp ; syice(:,:,:) = 0._wp ; sxxice(:,:,:) = 0._wp ; syyice(:,:,:) = 0._wp ; sxyice(:,:,:) = 0._wp
      sxsnw(:,:,:) = 0._wp ; sysnw(:,:,:) = 0._wp ; sxxsnw(:,:,:) = 0._wp ; syysnw(:,:,:) = 0._wp ; sxysnw(:,:,:) = 0._wp
      sxa(:,:,:) = 0._wp   ;   sya(:,:,:) = 0._wp ;   sxxa(:,:,:) = 0._wp ;   syya(:,:,:) = 0._wp ;   sxya(:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
      PRINT *, '            => sxice, syice, sxxice, syyice, sxyice, sxsnw, sysnw, sxxsnw, syysnw, sxysnw, sxa, sya, sxxa, syya, sxya'
      !$acc enter data copyin( sxice, syice, sxxice, syyice, sxyice, sxsnw, sysnw, sxxsnw, syysnw, sxysnw, sxa, sya, sxxa, syya, sxya )
#endif
      !
      ALLOCATE( sa3d(jpi,jpj,jpl), sati1(jpi,jpj),   STAT=ierr(2) )
      sa3d(:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
      PRINT *, '            => sa3d, sati1'
      !$acc enter data copyin( sa3d, sati1 )
#endif
      !
      IF( ln_rdgtrc .AND. .NOT. ALLOCATED(sa2d) ) THEN
         ALLOCATE( sa2d(jpi,jpj),   STAT=ierr(3) )
         sa2d(:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding array to memory'
         PRINT *, '            => sa2d'
         !$acc enter data copyin( sa2d )
#endif
      ENDIF
      !
      IF( .NOT. ln_pureADV2D ) THEN
         ALLOCATE( shi_max(jpi,jpj,jpl), shs_max(jpi,jpj,jpl),    STAT=ierr(4) )
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
         PRINT *, '            => shi_max, shs_max'
         !$acc enter data copyin( shi_max, shs_max )
#endif
      ENDIF
      !
      IF( ln_icethd ) THEN
         ALLOCATE( sxes (jpi,jpj,nlay_s,jpl) ,  syes (jpi,jpj,nlay_s,jpl) , sxxes(jpi,jpj,nlay_s,jpl) , &
            &      syyes(jpi,jpj,nlay_s,jpl) ,  sxyes(jpi,jpj,nlay_s,jpl)                             , &
            &      sxei  (jpi,jpj,nlay_i,jpl), syei  (jpi,jpj,nlay_i,jpl) , sxxei (jpi,jpj,nlay_i,jpl), &
            &      syyei (jpi,jpj,nlay_i,jpl), sxyei (jpi,jpj,nlay_i,jpl) ,       STAT=ierr(5) )
         sxes(:,:,:,:) = 0._wp ; syes(:,:,:,:) = 0._wp ; sxxes(:,:,:,:) = 0._wp ; syyes(:,:,:,:) = 0._wp ; sxyes(:,:,:,:) = 0._wp
         sxei(:,:,:,:) = 0._wp ; syei(:,:,:,:) = 0._wp ; sxxei(:,:,:,:) = 0._wp ; syyei(:,:,:,:) = 0._wp ; sxyei(:,:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
         PRINT *, '            => sxes, syes, sxxes, syyes, sxyes, sxei, syei, sxxei, syyei, sxyei'
         !$acc enter data copyin( sxes, syes, sxxes, syyes, sxyes, sxei, syei, sxxei, syyei, sxyei )
#endif

         ALLOCATE( sxsz  (jpi,jpj,nlay_i,jpl), sysz  (jpi,jpj,nlay_i,jpl) , sxxsz (jpi,jpj,nlay_i,jpl), &
            &      syysz (jpi,jpj,nlay_i,jpl), sxysz (jpi,jpj,nlay_i,jpl) ,    STAT=ierr(6) )
         sxsz(:,:,:,:) = 0._wp ; sysz(:,:,:,:) = 0._wp ; sxxsz(:,:,:,:) = 0._wp ; syysz(:,:,:,:) = 0._wp ; sxysz(:,:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
         PRINT *, '            => sxsz, sysz, sxxsz, syysz, sxysz'
         !$acc enter data copyin( sxsz, sysz, sxxsz, syysz, sxysz )
#endif

         ALLOCATE( sa4d(jpi,jpj,nlay_i,jpl),    STAT=ierr(7) )
         sa4d(:,:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
         PRINT *, '            => sa4d'
         !$acc enter data copyin( sa4d )
#endif
         IF( .NOT. ln_pureADV2D ) THEN
            ALLOCATE( sei_max(jpi,jpj,nlay_i,jpl), ses_max(jpi,jpj,nlay_s,jpl), sszi_max(jpi,jpj,nlay_i,jpl), STAT=ierr(8) )
#if defined _OPENACC || defined _OPENMP
            PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
            PRINT *, '            => sei_max, ses_max, sszi_max'
            !$acc enter data copyin( sei_max, ses_max, sszi_max )
#endif
         ENDIF

         IF( ln_age ) THEN
            ALLOCATE( sxage(jpi,jpj,jpl), syage(jpi,jpj,jpl), sxxage(jpi,jpj,jpl), syyage(jpi,jpj,jpl), sxyage(jpi,jpj,jpl),  STAT=ierr(9) )
            sxage(:,:,:) = 0._wp ; syage(:,:,:) = 0._wp ; sxxage(:,:,:) = 0._wp ; syyage(:,:,:) = 0._wp ; sxyage(:,:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
            PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
            PRINT *, '            => sxage, syage, sxxage, syyage, sxyage'
            !$acc enter data copyin( sxage, syage, sxxage, syyage, sxyage )
#endif
         ENDIF

      ENDIF !IF( ln_icethd )


      IF( ln_rdgtrc ) THEN
         ALLOCATE( sxrdg(jpi,jpj), syrdg(jpi,jpj), sxxrdg (jpi,jpj), syyrdg(jpi,jpj), sxyrdg(jpi,jpj),  STAT=ierr(10) )
         sxrdg(:,:) = 0._wp ; syrdg(:,:) = 0._wp ; sxxrdg(:,:) = 0._wp ; syyrdg(:,:) = 0._wp ; sxyrdg(:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
         PRINT *, ' * info GPU: adv_pra_init() => adding arrays to memory'
         PRINT *, '            => sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg'
         !$acc enter data copyin( sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg )
#endif
      ENDIF


      !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
      !   ALLOCATE( sxap (jpi,jpj,jpl) , syap (jpi,jpj,jpl) , sxxap (jpi,jpj,jpl) , syyap (jpi,jpj,jpl) , sxyap (jpi,jpj,jpl), &
      !      &      sxvp (jpi,jpj,jpl) , syvp (jpi,jpj,jpl) , sxxvp (jpi,jpj,jpl) , syyvp (jpi,jpj,jpl) , sxyvp (jpi,jpj,jpl), &
      !      &      STAT=ierr(10) )
      !   IF ( ln_pnd_lids ) THEN
      !      ALLOCATE( sxvl (jpi,jpj,jpl) , syvl (jpi,jpj,jpl) , sxxvl (jpi,jpj,jpl) , syyvl (jpi,jpj,jpl) , sxyvl (jpi,jpj,jpl), &
      !         &      STAT=ierr(11) )
      !   ENDIF
      !ENDIF
      !
      ALLOCATE( zfld(jpi,jpj), zf0(jpi,jpj), zbet(jpi,jpj), zfm(jpi,jpj), zfx(jpi,jpj), zfy(jpi,jpj),  &
         &      zfxx(jpi,jpj), zfyy(jpi,jpj), zfxy(jpi,jpj), zpm(jpi,jpj), zpx(jpi,jpj), zpy(jpi,jpj), &
         &      zpxx(jpi,jpj), zpyy(jpi,jpj), zpxy(jpi,jpj), zalg(jpi,jpj), zalg1(jpi,jpj), zalg1q(jpi,jpj), STAT=ierr(12) )
      zfld(:,:) = 0._wp; zf0(:,:) = 0._wp; zbet(:,:) = 0._wp; zfm(:,:) = 0._wp; zfx(:,:) = 0._wp; zfy(:,:) = 0._wp
      zfxx(:,:) = 0._wp; zfyy(:,:) = 0._wp; zfxy(:,:) = 0._wp; zpm(:,:) = 0._wp; zpx(:,:) = 0._wp; zpy(:,:) = 0._wp
      zpxx(:,:) = 0._wp; zpyy(:,:) = 0._wp; zpxy(:,:) = 0._wp; zalg(:,:) = 0._wp; zalg1(:,:) = 0._wp; zalg1q(:,:) = 0._wp
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: adv_pra_init() => adding Prather advection workspace arrays to memory'
      PRINT *, '            => zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy'
      PRINT *, '            => zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q'
      !$acc enter data copyin( zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )
#endif
      !
      k_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( crtnm, k_alloc )
      IF( k_alloc > 0 ) CALL ctl_stop('STOP', crtnm//': unable to allocate ice arrays for Prather advection scheme')
      !
      CALL adv_pra_rst( 'READ' )    !* read or initialize all required files

   END SUBROUTINE adv_pra_init


   SUBROUTINE adv_pra_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE adv_pra_rst  ***
      !!
      !! ** Purpose :   Read or write file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER ::   jk, jl
      INTEGER ::   iter
      INTEGER ::   id1, id0, idx
      LOGICAL ::   ll_lay   ! the number of snow & ice layers of the restart matches `nlay_s` & `nlay_i`
      CHARACTER(len=25) ::   znam
      CHARACTER(len=2)  ::   zchar, zchar1
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      !                                      !==========================!
      IF( TRIM(cdrw) == 'READ' ) THEN        !==  Read or initialize  ==!
         !                                   !==========================!
         !
         id1 = 0                                                                              ! no restart: id1=0
         IF( ln_rstart ) id1 = iom_varid( 'adv_pra_rst', numrir, 'sxice' , ldstop = .FALSE. ) ! file exist: id1>0
         !
         ! check size of the input fields
         id0=0
         DO jk = 1, 99
            WRITE(zchar1,'(I2.2)') jk
            znam = 'sxes'//'_l'//zchar1
            IF( iom_varid( 'adv_pra_rst', numrir, znam , ldstop = .FALSE. ) > 0 )   id0 = id0+1
            znam = 'sxei'//'_l'//zchar1
            IF( iom_varid( 'adv_pra_rst', numrir, znam , ldstop = .FALSE. ) > 0 )   id0 = id0+1
         END DO

         ! The mandatory fields (ice & snow volume, ice concentration) are read as soon as a restart exists. The layer-dependent
         ! fields (salinity, snow & ice enthalpy), which only exist if the run that wrote the restart had `ln_icethd=T`, are
         ! checked separately below (`ll_lay`): a mismatch there must not reset the moments that can be read correctly.
         IF( id1 > 0 ) THEN   !**  Read the restart file  **!
            !
            ! ---------------------- !
            ! == mandatory fields == !
            ! ---------------------- !
            !                                                        ! ice thickness
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxice' , sxice , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syice' , syice , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxice', sxxice )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyice', syyice )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyice', sxyice )
            !                                                        ! snow thickness
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxsnw' , sxsnw , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sysnw' , sysnw , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxsnw', sxxsnw )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syysnw', syysnw )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxysnw', sxysnw )
            !                                                        ! ice concentration
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxa'   , sxa   , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sya'   , sya   , psgn = -1._wp )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxa'  , sxxa   )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syya'  , syya   )
            CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxya'  , sxya   )

            ! The Prather moments of the ridged-ice fraction are now those of the VOLUME of ridged ice (cf. `ice_dyn_adv_pra`),
            ! hence the new names `s*vrdg`: the moments of an older restart (moments of the fraction) are not read, but set to 0.
            IF( ln_rdgtrc ) THEN                                    ! ridged-ice fraction
               idx = iom_varid( 'adv_pra_rst',numrir, 'sxvrdg', ldstop = .FALSE. )
               IF( idx > 0 ) THEN      ! field exists
                  CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxvrdg'  , sxrdg   , psgn = -1._wp )
                  CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syvrdg'  , syrdg   , psgn = -1._wp )
                  CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxvrdg' , sxxrdg   )
                  CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyvrdg' , syyrdg   )
                  CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyvrdg' , sxyrdg   )
               ELSE
                  IF(lwp) WRITE(numout,*)
                  IF(lwp) WRITE(numout,*) '   ==>>> did not find the "ridged fraction"-related Prather adv. fields in the restart => set to 0'
                  sxrdg(:,:) = 0._wp ; syrdg(:,:) = 0._wp ; sxxrdg(:,:) = 0._wp ; syyrdg(:,:) = 0._wp ; sxyrdg(:,:) = 0._wp
               ENDIF
            ENDIF
            !
            IF( ln_icethd ) THEN
               !
               ! Do the snow & ice layers of the restart match the ones of this run? (`id0` counts the layers found for the
               ! snow and ice enthalpy moments; it is 0 if the restart was written by a run without thermodynamics)
               ll_lay = ( id0 == (nlay_s+nlay_i) )
               IF( .NOT. ll_lay ) THEN
                  ! No match: only the layer-dependent moments are set to 0 (all the other moments are read normally):
                  IF(lwp) WRITE(numout,*)
                  IF(lwp) WRITE(numout,*) '   ==>>> layers of the restart do not match `nlay_i`/`nlay_s` (or no thermo in previous run)' &
                     &                    //' => salinity & enthalpy Prather adv. fields set to 0'
                  sxsz = 0._wp   ;   sysz = 0._wp   ;   sxxsz = 0._wp   ;   syysz = 0._wp   ;   sxysz = 0._wp
                  sxes = 0._wp   ;   syes = 0._wp   ;   sxxes = 0._wp   ;   syyes = 0._wp   ;   sxyes = 0._wp
                  sxei = 0._wp   ;   syei = 0._wp   ;   sxxei = 0._wp   ;   syyei = 0._wp   ;   sxyei = 0._wp
               ENDIF
               !                                                        ! ice salinity
               IF( ll_lay ) THEN
                  DO jk = 1, nlay_i
                     WRITE(zchar1,'(I2.2)') jk
                     znam = 'sxsz'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sxsz (:,:,jk,:) = z3d(:,:,:)
                     znam = 'sysz'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sysz (:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxxsz'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxxsz(:,:,jk,:) = z3d(:,:,:)
                     znam = 'syysz'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   syysz(:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxysz'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxysz(:,:,jk,:) = z3d(:,:,:)
                  END DO
               ENDIF !IF( ll_lay )

               IF( ln_age ) THEN                                        ! ice age
                  idx = iom_varid( 'adv_pra_rst',numrir, 'sxage' , ldstop = .FALSE. )
                  IF( idx > 0 ) THEN      ! field exists
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxage' , sxage , psgn = -1._wp )
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syage' , syage , psgn = -1._wp )
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxage', sxxage )
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyage', syyage )
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyage', sxyage )
                  ELSE
                     IF(lwp) WRITE(numout,*)
                     IF(lwp) WRITE(numout,*) '   ==>>> did not find the "age"-related Prather adv. fields in the restart => set to 0'
                     sxage(:,:,:) = 0._wp ; syage(:,:,:) = 0._wp ; sxxage(:,:,:) = 0._wp ; syyage(:,:,:) = 0._wp ; sxyage(:,:,:) = 0._wp
                  ENDIF
               ENDIF
               !                                                        ! snow layers heat content
               IF( ll_lay ) THEN
                  DO jk = 1, nlay_s
                     WRITE(zchar1,'(I2.2)') jk
                     znam = 'sxes'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sxes (:,:,jk,:) = z3d(:,:,:)
                     znam = 'syes'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   syes (:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxxes'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxxes(:,:,jk,:) = z3d(:,:,:)
                     znam = 'syyes'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   syyes(:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxyes'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxyes(:,:,jk,:) = z3d(:,:,:)
                  END DO
                  !                                                        ! ice layers heat content
                  DO jk = 1, nlay_i
                     WRITE(zchar1,'(I2.2)') jk
                     znam = 'sxei'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sxei (:,:,jk,:) = z3d(:,:,:)
                     znam = 'syei'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   syei (:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxxei'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxxei(:,:,jk,:) = z3d(:,:,:)
                     znam = 'syyei'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   syyei(:,:,jk,:) = z3d(:,:,:)
                     znam = 'sxyei'//'_l'//zchar1
                     CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, znam , z3d )   ;   sxyei(:,:,jk,:) = z3d(:,:,:)
                  END DO
               ENDIF !IF( ll_lay ) (this `IF` was opened before the snow layers)
               !
            ENDIF
            !
            !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN                                    ! melt pond fraction
            !   IF( iom_varid( 'adv_pra_rst', numrir, 'sxap', ldstop = .FALSE. ) > 0 ) THEN
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxap' , sxap , psgn = -1._wp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syap' , syap , psgn = -1._wp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxap', sxxap )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyap', syyap )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyap', sxyap )
            !      !                                                     ! melt pond volume
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxvp' , sxvp , psgn = -1._wp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syvp' , syvp , psgn = -1._wp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxvp', sxxvp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyvp', syyvp )
            !      CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyvp', sxyvp )
            !   ELSE
            !      sxap = 0._wp ;   syap = 0._wp    ;   sxxap = 0._wp    ;   syyap = 0._wp    ;   sxyap = 0._wp   ! melt pond fraction
            !      sxvp = 0._wp ;   syvp = 0._wp    ;   sxxvp = 0._wp    ;   syyvp = 0._wp    ;   sxyvp = 0._wp   ! melt pond volume
            !   ENDIF
            !   !
            !   IF ( ln_pnd_lids ) THEN                               ! melt pond lid volume
            !      IF( iom_varid( 'adv_pra_rst', numrir, 'sxvl', ldstop = .FALSE. ) > 0 ) THEN
            !         CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxvl' , sxvl , psgn = -1._wp )
            !         CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syvl' , syvl , psgn = -1._wp )
            !         CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxxvl', sxxvl )
            !         CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'syyvl', syyvl )
            !         CALL iom_get( 'adv_pra_rst', numrir, jpdom_auto, 'sxyvl', sxyvl )
            !      ELSE
            !         sxvl = 0._wp; syvl = 0._wp    ;   sxxvl = 0._wp    ;   syyvl = 0._wp    ;   sxyvl = 0._wp   ! melt pond lid volume
            !      ENDIF
            !   ENDIF
            !ENDIF
            !
         ELSE                                   !**  start rheology from rest  **!
            !
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest OR previous run without Prather, set moments to 0'
            !
            sxice = 0._wp   ;   syice = 0._wp   ;   sxxice = 0._wp   ;   syyice = 0._wp   ;   sxyice = 0._wp      ! ice thickness
            sxsnw = 0._wp   ;   sysnw = 0._wp   ;   sxxsnw = 0._wp   ;   syysnw = 0._wp   ;   sxysnw = 0._wp      ! snow thickness
            sxa   = 0._wp   ;   sya   = 0._wp   ;   sxxa   = 0._wp   ;   syya   = 0._wp   ;   sxya   = 0._wp      ! ice concentration
            IF( ln_rdgtrc ) THEN
               sxrdg = 0._wp   ;   syrdg = 0._wp   ;   sxxrdg = 0._wp   ;   syyrdg = 0._wp   ;   sxyrdg = 0._wp
            ENDIF
            IF( ln_icethd ) THEN
               sxsz = 0._wp   ;   sysz = 0._wp   ;   sxxsz = 0._wp   ;   syysz = 0._wp   ;   sxysz = 0._wp      ! ice salinity
               sxes  = 0._wp   ;   syes  = 0._wp   ;   sxxes  = 0._wp   ;   syyes  = 0._wp   ;   sxyes  = 0._wp      ! snow layers heat content
               sxei   = 0._wp   ;   syei   = 0._wp   ;   sxxei   = 0._wp   ;   syyei   = 0._wp   ;   sxyei   = 0._wp      ! ice layers heat content
               IF( ln_age ) THEN
                  sxage = 0._wp   ;   syage = 0._wp   ;   sxxage = 0._wp   ;   syyage = 0._wp   ;   sxyage = 0._wp      ! ice age
               ENDIF
            ENDIF
            !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   sxap = 0._wp ;   syap = 0._wp    ;   sxxap = 0._wp    ;   syyap = 0._wp    ;   sxyap = 0._wp       ! melt pond fraction
            !   sxvp = 0._wp ;   syvp = 0._wp    ;   sxxvp = 0._wp    ;   syyvp = 0._wp    ;   sxyvp = 0._wp       ! melt pond volume
            !   IF ( ln_pnd_lids ) THEN
            !      sxvl = 0._wp; syvl = 0._wp    ;   sxxvl = 0._wp    ;   syyvl = 0._wp    ;   sxyvl = 0._wp       ! melt pond lid volume
            !   ENDIF
            !ENDIF
         ENDIF !IF( id1 > 0 )
         !
#if defined _OPENACC || defined _OPENMP
         !$acc update device( sxice, syice, sxxice, syyice, sxyice, sxsnw, sysnw, sxxsnw, syysnw, sxysnw, sxa, sya, sxxa, syya, sxya )
         IF( ln_rdgtrc ) THEN
            !$acc update device( sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg )
         ENDIF
         IF( ln_icethd) THEN
            !$acc update device( sxes, syes, sxxes, syyes, sxyes, sxei, syei, sxxei, syyei, sxyei )
            !$acc update device( sxsz, sysz, sxxsz, syysz, sxysz )
            IF( ln_age ) THEN
               !$acc update device( sxage, syage, sxxage, syyage, sxyage )
            ENDIF
         ENDIF
#endif
         !
         !                                   !=====================================!
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   !==  write in the ice restart file  ==!
         !                                   !=====================================!
         IF(lwp) WRITE(numout,*) '----  ice-adv-rst  ----'
         iter = kt              ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         !
         ! Update CPU with GPU's content:
         !$acc update self( sxice, syice, sxxice, syyice, sxyice )
         !$acc update self( sxsnw, sysnw, sxxsnw, syysnw, sxysnw )
         !$acc update self( sxa  , sya  , sxxa  , syya  , sxya   )
         !LOLO:PND !LOLOfixme
         !0acc update self( sxap , syap , sxxap , syyap , sxyap  )
         !0acc update self( sxvp , syvp , sxxvp , syyvp , sxyvp  )
         !0acc update self( sxvl , syvl , sxxvl , syyvl , sxyvl  )
         !
         ! In case Prather scheme is used for advection, write second order moments
         ! ------------------------------------------------------------------------
         !
         ! ---------------------- !
         ! == mandatory fields == !
         ! ---------------------- !
         !                                                           ! ice thickness
         CALL iom_rstput( kt, nitrst, numriw, 'sxice' , sxice  )
         CALL iom_rstput( kt, nitrst, numriw, 'syice' , syice  )
         CALL iom_rstput( kt, nitrst, numriw, 'sxxice', sxxice )
         CALL iom_rstput( kt, nitrst, numriw, 'syyice', syyice )
         CALL iom_rstput( kt, nitrst, numriw, 'sxyice', sxyice )
         !                                                           ! snow thickness
         CALL iom_rstput( kt, nitrst, numriw, 'sxsnw' , sxsnw  )
         CALL iom_rstput( kt, nitrst, numriw, 'sysnw' , sysnw  )
         CALL iom_rstput( kt, nitrst, numriw, 'sxxsnw', sxxsnw )
         CALL iom_rstput( kt, nitrst, numriw, 'syysnw', syysnw )
         CALL iom_rstput( kt, nitrst, numriw, 'sxysnw', sxysnw )
         !                                                           ! ice concentration
         CALL iom_rstput( kt, nitrst, numriw, 'sxa'   , sxa    )
         CALL iom_rstput( kt, nitrst, numriw, 'sya'   , sya    )
         CALL iom_rstput( kt, nitrst, numriw, 'sxxa'  , sxxa   )
         CALL iom_rstput( kt, nitrst, numriw, 'syya'  , syya   )
         CALL iom_rstput( kt, nitrst, numriw, 'sxya'  , sxya   )
         !                                                           ! snow layers heat content
         IF( ln_icethd ) THEN
            !                                                           ! ice salinity
            !$acc update self( sxsz, sysz, sxxsz, syysz, sxysz )
            DO jk = 1, nlay_i
               WRITE(zchar1,'(I2.2)') jk
               znam = 'sxsz'//'_l'//zchar1   ;   z3d(:,:,:) = sxsz (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sysz'//'_l'//zchar1   ;   z3d(:,:,:) = sysz (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxxsz'//'_l'//zchar1  ;   z3d(:,:,:) = sxxsz(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'syysz'//'_l'//zchar1  ;   z3d(:,:,:) = syysz(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxysz'//'_l'//zchar1  ;   z3d(:,:,:) = sxysz(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
            END DO
            !
            !$acc update self( sxei  , syei  , sxxei  , syyei  , sxyei   )
            DO jk = 1, nlay_i
               WRITE(zchar1,'(I2.2)') jk
               znam = 'sxei'//'_l'//zchar1   ;   z3d(:,:,:) = sxei (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'syei'//'_l'//zchar1   ;   z3d(:,:,:) = syei (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxxei'//'_l'//zchar1  ;   z3d(:,:,:) = sxxei(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'syyei'//'_l'//zchar1  ;   z3d(:,:,:) = syyei(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxyei'//'_l'//zchar1  ;   z3d(:,:,:) = sxyei(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
            END DO
            !
            !$acc update self( sxes , syes , sxxes , syyes , sxyes  )
            DO jk = 1, nlay_s
               WRITE(zchar1,'(I2.2)') jk
               znam = 'sxes'//'_l'//zchar1  ;   z3d(:,:,:) = sxes (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'syes'//'_l'//zchar1  ;   z3d(:,:,:) = syes (:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxxes'//'_l'//zchar1 ;   z3d(:,:,:) = sxxes(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'syyes'//'_l'//zchar1 ;   z3d(:,:,:) = syyes(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
               znam = 'sxyes'//'_l'//zchar1 ;   z3d(:,:,:) = sxyes(:,:,jk,:)
               CALL iom_rstput( kt, nitrst, numriw, znam , z3d )
            END DO
            !
            IF( ln_age ) THEN
               !$acc update self( sxage, syage, sxxage, syyage, sxyage )
               CALL iom_rstput( kt, nitrst, numriw, 'sxage' , sxage  )
               CALL iom_rstput( kt, nitrst, numriw, 'syage' , syage  )
               CALL iom_rstput( kt, nitrst, numriw, 'sxxage', sxxage )
               CALL iom_rstput( kt, nitrst, numriw, 'syyage', syyage )
               CALL iom_rstput( kt, nitrst, numriw, 'sxyage', sxyage )
            ENDIF
            !
         ENDIF !IF( ln_icethd )
         !
         IF( ln_rdgtrc ) THEN    !  ridged-ice fraction
            !$acc update self( sxrdg, syrdg, sxxrdg, syyrdg, sxyrdg )
            CALL iom_rstput( kt, nitrst, numriw, 'sxvrdg' ,   sxrdg(:,:) )
            CALL iom_rstput( kt, nitrst, numriw, 'syvrdg'  ,  syrdg(:,:) )
            CALL iom_rstput( kt, nitrst, numriw, 'sxxvrdg' , sxxrdg(:,:) )
            CALL iom_rstput( kt, nitrst, numriw, 'syyvrdg' , syyrdg(:,:) )
            CALL iom_rstput( kt, nitrst, numriw, 'sxyvrdg' , sxyrdg(:,:) )
         ENDIF
         !
         !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN                       !  melt pond fraction
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxap' , sxap  )
         !   CALL iom_rstput( kt, nitrst, numriw, 'syap' , syap  )
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxxap', sxxap )
         !   CALL iom_rstput( kt, nitrst, numriw, 'syyap', syyap )
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxyap', sxyap )
         !   !                                                        ! melt pond volume
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxvp' , sxvp  )
         !   CALL iom_rstput( kt, nitrst, numriw, 'syvp' , syvp  )
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxxvp', sxxvp )
         !   CALL iom_rstput( kt, nitrst, numriw, 'syyvp', syyvp )
         !   CALL iom_rstput( kt, nitrst, numriw, 'sxyvp', sxyvp )
         !   !
         !   IF ( ln_pnd_lids ) THEN                                  ! melt pond lid volume
         !      CALL iom_rstput( kt, nitrst, numriw, 'sxvl' , sxvl  )
         !      CALL iom_rstput( kt, nitrst, numriw, 'syvl' , syvl  )
         !      CALL iom_rstput( kt, nitrst, numriw, 'sxxvl', sxxvl )
         !      CALL iom_rstput( kt, nitrst, numriw, 'syyvl', syyvl )
         !      CALL iom_rstput( kt, nitrst, numriw, 'sxyvl', sxyvl )
         !   ENDIF
         !ENDIF
         !
      ENDIF
      !
   END SUBROUTINE adv_pra_rst

   !!======================================================================
END MODULE icedyn_adv_pra
