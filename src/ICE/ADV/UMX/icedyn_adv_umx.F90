MODULE icedyn_adv_umx
   !!==============================================================================
   !!                       ***  MODULE  icedyn_adv_umx  ***
   !! sea-ice : advection using the ULTIMATE-MACHO scheme
   !!==============================================================================
   !! History :  3.6  !  2014-11  (C. Rousset, G. Madec)  Original code
   !!            4.0  !  2018     (many people)           SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_umx   : update the tracer fields
   !!   ultimate_x(_y)    : compute a tracer value at velocity points using ULTIMATE scheme at various orders
   !!   macho             : compute the fluxes
   !!   nonosc_ice        : limit the fluxes using a non-oscillatory algorithm
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice            ! sea-ice variables
   USE icevar,   ONLY : ice_var_zapneg
   USE icedyn_adv_util
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! to use sign with key_nosignedzero
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   adv_umx_init
   PUBLIC   ice_dyn_adv_umx   ! called by icedyn_adv.F90
   !
   INTEGER, PARAMETER ::   np_advS = 2         ! advection for S and T:    dVS/dt = -div(      uVS     ) => np_advS = 1
   !                                                                    or dVS/dt = -div( uA * uHS / u ) => np_advS = 2
   !                                                                    or dVS/dt = -div( uV * uS  / u ) => np_advS = 3
   INTEGER, PARAMETER ::   np_limiter = 1      ! limiter: 1 = nonosc
   !                                                      2 = superbee
   !                                                      3 = h3
   LOGICAL            ::   ll_upsxy  = .TRUE.   ! alternate directions for upstream
   LOGICAL            ::   ll_hoxy   = .TRUE.   ! alternate directions for high order
   LOGICAL            ::   ll_neg    = .TRUE.   ! if T interpolated at u/v points is negative or v_i < 1.e-6
   !                                                 then interpolate T at u/v points using the upstream scheme
   LOGICAL            ::   ll_prelim = .FALSE.  ! prelimiter from: Zalesak(1979) eq. 14 => not well defined in 2D
   !
   REAL(wp)           ::   r1_6   = 1._wp /   6._wp   ! =1/6
   REAL(wp)           ::   r1_120 = 1._wp / 120._wp   ! =1/120

   ! Work array (that should remain once for all in the memory of the GPU)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   sati1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   shi_max, shs_max, ssi_max
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sei_max, sszi_max
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   ses_max
   !
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) ::   imsk_small, jmsk_small

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_umx( kn_umx, kt, pUu, pVv, ph_i, ph_s, ph_ip,  &
      &                        pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_umx  ***
      !!
      !! **  Purpose :   Compute the now trend due to total advection of
      !!                 tracers and add it to the general trend of tracer equations
      !!                 using an "Ultimate-Macho" scheme
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74.
      !!----------------------------------------------------------------------
      INTEGER                                , INTENT(in   ) ::   kn_umx ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                                , INTENT(in   ) ::   kt     ! time step
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pUu    ! ice i-velocity
      REAL(wp), DIMENSION(jpi,jpj)           , INTENT(in   ) ::   pVv    ! ice j-velocity
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
      INTEGER  ::   ji, jj, jk, jl, jm      ! dummy loop indices
      INTEGER  ::   ndim                    ! number of variables to advect
      REAL(wp) ::   zdt, z1_dt, zvi_cen
      !
      CHARACTER(len=19), PARAMETER :: crtnm = 'ice_dyn_adv_umx'
      REAL(wp) ::   zati2
      REAL(wp), DIMENSION(jpi,jpj)            ::   zcu_box, zcv_box
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   zu_cat, zv_cat
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   zua_ho, zva_ho, zua_ups, zva_ups
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   zuap_ho, zvap_ho, zuap_ups, zvap_ups
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   z1_ai
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ze_i
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ze_s
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   z1_aip
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zs_i
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zsz_i
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:)       ::   zamsk                   ! 1 if advection of concentration, 0 if advection of other tracers
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zvar, zhvar
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zuv_ho, zvv_ho, zuv_ups, zvv_ups, z1_vi, z1_vs
      !! diagnostics
      REAL(wp), DIMENSION(jpi,jpj) ::   zdiag_adv_mass, zdiag_adv_salt, zdiag_adv_heat
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
      !
      ndim = nlay_s + 2*nlay_i + 5 ! max number of tracers to advect at the same time
      !
      ! --- Allocate arrays --- !
      ALLOCATE( zvar(jpi,jpj,jpl,ndim), zhvar(jpi,jpj,jpl,ndim), zamsk(ndim) )
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO )   ALLOCATE( z1_aip(jpi,jpj,jpl) )
      IF( np_advS == 3 )   ALLOCATE( zuv_ho(jpi,jpj,jpl), zvv_ho(jpi,jpj,jpl), zuv_ups(jpi,jpj,jpl), zvv_ups(jpi,jpj,jpl), &
         &                           z1_vi (jpi,jpj,jpl), z1_vs (jpi,jpj,jpl) )

      IF( (kt == nit000) .AND. lwp )   WRITE(numout,*) '-- '//crtnm//': Ultimate-Macho advection scheme'

      IF( .NOT. ln_pureADV2D ) THEN
         ! --- Record max of the surrounding 9-pts (for call Hbig) --- !
         ! thickness and salinity
         CALL icemax3D(     ph_i ,       shi_max )
         CALL icemax3D(     ph_s ,       shs_max )
         !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) CALL icemax3D( ph_ip, zhip_max)

         IF( ln_icethd ) THEN
            ! Salt content
            IF( nn_icesal == 4 ) THEN
               CALL icemax4D_cnt( nlay_i, pszv_i, pv_i, sszi_max )
            ELSE
               CALL icemax3D_cnt( psv_i,          pv_i,  ssi_max )
            ENDIF
            ! Enthalpies
            CALL icemax4D_cnt( nlay_i, pe_i,   pv_i,  sei_max )
            CALL icemax4D_cnt( nlay_s, pe_s,   pv_s,  ses_max )
         ENDIF !IF( ln_icethd )

         !! LOLO: with my wide halos, the following might not be needed:
         IF( ln_icethd ) THEN
            IF( nn_icesal == 4 ) THEN
               CALL lbc_lnk( crtnm, shi_max,'T',1._wp, shs_max,'T',1._wp )   ! 3D
               CALL lbc_lnk( crtnm, sei_max,'T',1._wp, sszi_max,'T',1._wp )  ! 4D / nlay_i
            ELSE
               CALL lbc_lnk( crtnm, shi_max,'T',1._wp, shs_max,'T',1._wp,  ssi_max,'T',1._wp )   ! 3D
               CALL lbc_lnk( crtnm, sei_max,'T',1._wp )                      ! 4D / nlay_i
            ENDIF
            CALL lbc_lnk(    crtnm, ses_max,'T',1._wp )                      ! 4D / nlay_s
         ELSE
            CALL lbc_lnk( crtnm, shi_max,'T',1._wp, shs_max,'T',1._wp ) ! 3D
         ENDIF

      ENDIF !IF( .NOT. ln_pureADV2D )

      zdt = rDt_ice
      z1_dt = 1._wp / zdt

      ! --- transport --- !
      !      + record at_i before advection (for open water)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            sati1(ji,jj) = 0._wp
            !%acc loop seq
            DO jl=1, jpl
               sati1(ji,jj) = sati1(ji,jj) + pa_i(ji,jj,jl)
            END DO
            !
         END DO
      END DO

      !
      ! setup transport for each ice cat
      DO jl = 1, jpl
         zu_cat(:,:,jl) = pUu(:,:)
         zv_cat(:,:,jl) = pVv(:,:)
      END DO
      !
      ! --- define velocity for advection: u*grad(H) --- !
      DO jj=Njs0-2, Nje0+2
         DO ji=Nis0-1, Nie0+2
            IF    ( pUu(ji,jj) * pUu(ji-1,jj) <= 0._wp ) THEN
               zcu_box(ji,jj) = 0._wp
            ELSEIF( pUu(ji,jj)                   >  0._wp ) THEN
               zcu_box(ji,jj) = pUu(ji-1,jj)
            ELSE
               zcu_box(ji,jj) = pUu(ji  ,jj)
            ENDIF
         END DO
      END DO
      DO jj=Njs0-1, Nje0+2
         DO ji=Nis0-2, Nie0+2
            IF    ( pVv(ji,jj) * pVv(ji,jj-1) <= 0._wp ) THEN
               zcv_box(ji,jj) = 0._wp
            ELSEIF( pVv(ji,jj)                   >  0._wp ) THEN
               zcv_box(ji,jj) = pVv(ji,jj-1)
            ELSE
               zcv_box(ji,jj) = pVv(ji,jj  )
            ENDIF
         END DO
      END DO


      !---------------!
      !== advection ==!
      !---------------!

      ! inverse of A and Ap
      WHERE( pa_i(:,:,:) >= epsi20 )
         z1_ai(:,:,:) = 1._wp / pa_i(:,:,:)
      ELSEWHERE
         z1_ai(:,:,:) = 0.
      END WHERE
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         WHERE( pa_ip(:,:,:) >= epsi20 )
            z1_aip(:,:,:) = 1._wp / pa_ip(:,:,:)
         ELSEWHERE
            z1_aip(:,:,:) = 0.
         END WHERE
      ENDIF
      !
      ! setup a mask where advection will be upstream
      IF( ll_neg ) THEN
         IF( .NOT. ALLOCATED(imsk_small) )   ALLOCATE( imsk_small(jpi,jpj,jpl) )
         IF( .NOT. ALLOCATED(jmsk_small) )   ALLOCATE( jmsk_small(jpi,jpj,jpl) )
         DO jl = 1, jpl
            DO jj=Njs0-2, Nje0+2
               DO ji=Nis0-2, Nie0+1
                  zvi_cen = 0.5_wp * ( pv_i(ji+1,jj,jl) + pv_i(ji,jj,jl) )
                  IF( zvi_cen < epsi06) THEN
                     imsk_small(ji,jj,jl) = 0
                  ELSE
                     imsk_small(ji,jj,jl) = 1
                  ENDIF
               END DO
            END DO
            DO jj=Njs0-2, Nje0+1
               DO ji=Nis0-2, Nie0+2
                  zvi_cen = 0.5_wp * ( pv_i(ji,jj+1,jl) + pv_i(ji,jj,jl) )
                  IF( zvi_cen < epsi06) THEN
                     jmsk_small(ji,jj,jl) = 0
                  ELSE
                     jmsk_small(ji,jj,jl) = 1
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !
      ! diagnostics
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            zdiag_adv_mass(ji,jj) =   SUM( pv_i (ji,jj,:) ) * rhoi + SUM( pv_s (ji,jj,:) ) * rhos &
               &                    + SUM( pv_ip(ji,jj,:) ) * rhow + SUM( pv_il(ji,jj,:) ) * rhow
            zdiag_adv_heat(ji,jj) = - SUM( SUM( pe_i(ji,jj,1:nlay_i,:), dim=2 ) ) - SUM( SUM( pe_s(ji,jj,1:nlay_s,:), dim=2 ) )
         END DO
      END DO
      IF( nn_icesal == 4 ) THEN
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zdiag_adv_salt(ji,jj) = SUM( SUM( pszv_i(ji,jj,:,:), dim=2 ) ) * rhoi
            END DO
         END DO
      ELSE
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zdiag_adv_salt(ji,jj) = SUM( psv_i(ji,jj,:) ) * rhoi
            END DO
         END DO
      ENDIF
      !
      ! record at_i before advection (for open water)
      sati1(:,:) = SUM( pa_i(:,:,:), dim=3 )
      !
      ! ----------------------- !
      ! ==> start advection <== !
      ! ----------------------- !
      !
      zamsk(1) = 1._wp
      !
      !== Ice area ==!
      zvar(:,:,:,1) = pa_i(:,:,:)
      CALL adv_umx( zamsk(1:1), kn_umx, kt, zdt, pUu, pVv, zu_cat , zv_cat , zcu_box, zcv_box, &
         &                                           zvar(:,:,:,1:1), zvar(:,:,:,1:1), zua_ups, zva_ups, zua_ho, zva_ho )
      pa_i(:,:,:) = zvar(:,:,:,1)

      !== Ice age ==!
      zvar(:,:,:,1) = poa_i(:,:,:)
      CALL adv_umx( zamsk(1:1), kn_umx, kt, zdt, pUu , pVv , zu_cat, zv_cat, zcu_box, zcv_box, &
         &                                           zvar(:,:,:,1:1), zvar(:,:,:,1:1) )
      poa_i(:,:,:) = zvar(:,:,:,1)

      !== melt ponds area ==!
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         zvar(:,:,:,1) = pa_ip(:,:,:)
         CALL adv_umx( zamsk(1:1), kn_umx, kt, zdt, pUu , pVv , zu_cat , zv_cat , zcu_box, zcv_box,  &
            &                                           zvar(:,:,:,1:1), zvar(:,:,:,1:1), zuap_ups, zvap_ups, zuap_ho, zvap_ho )
         pa_ip(:,:,:) = zvar(:,:,:,1)
      ENDIF
      !
      !                             ! --------------------------------- !
      IF( np_advS == 1 ) THEN       ! -- advection form: -div( uVS ) -- !
         !                          ! --------------------------------- !

         !== Ice volume ==!
         jm = 1    ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_i(:,:,:)
         zhvar(:,:,:,jm) = pv_i(:,:,:) * z1_ai(:,:,:)
         !== Snw volume ==!
         jm = jm+1 ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_s(:,:,:)
         zhvar(:,:,:,jm) = pv_s(:,:,:) * z1_ai(:,:,:)
         !== Ice heat content ==!
         DO jk = 1, nlay_i
            jm = jm+1 ; zamsk(jm) = 1._wp
            zvar (:,:,:,jm) = pe_i(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_i(:,:,jk,:)
         ENDDO
         !== Snw heat content ==!
         DO jk = 1, nlay_s
            jm = jm+1 ; zamsk(jm) = 1._wp
            zvar (:,:,:,jm) = pe_s(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_s(:,:,jk,:)
         ENDDO
         !== Ice salt content ==!
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; zamsk(jm) = 1._wp
               zvar (:,:,:,jm) = pszv_i(:,:,jk,:)
               zhvar(:,:,:,jm) = pszv_i(:,:,jk,:)
            ENDDO
         ELSE
            jm = jm+1 ; zamsk(jm) = 1._wp
            zvar (:,:,:,jm) = psv_i(:,:,:)
            zhvar(:,:,:,jm) = psv_i(:,:,:)
         ENDIF

         !== advection ==!
         CALL adv_umx( zamsk(1:jm), kn_umx, kt, zdt, pUu, pVv, zua_ho, zva_ho, zcu_box, zcv_box, &
            &                                            zhvar(:,:,:,1:jm), zvar(:,:,:,1:jm), zua_ups, zva_ups )
         !

         !== Recover quantities ==!
         jm = 1       ;    pv_i  (:,:,:)    = zvar (:,:,:,jm)
         jm = jm+1    ;    pv_s  (:,:,:)    = zvar (:,:,:,jm)
         DO jk = 1, nlay_i
            jm = jm+1 ;    pe_i  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         DO jk = 1, nlay_s
            jm = jm+1 ;    pe_s  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; pszv_i(:,:,jk,:) = zvar (:,:,:,jm)
            ENDDO
         ELSE
            jm = jm+1 ;    psv_i (:,:,:)    = zvar (:,:,:,jm)
         ENDIF

         !
         !                          ! ------------------------------------------ !
      ELSEIF( np_advS == 2 ) THEN   ! -- advection form: -div( uA * uHS / u ) -- !
         !                          ! ------------------------------------------ !

         !== Ice volume ==!
         jm = 1    ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_i(:,:,:)
         zhvar(:,:,:,jm) = pv_i(:,:,:) * z1_ai(:,:,:)
         !== Snw volume ==!
         jm = jm+1 ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_s(:,:,:)
         zhvar(:,:,:,jm) = pv_s(:,:,:) * z1_ai(:,:,:)
         !== Ice heat content ==!
         DO jk = 1, nlay_i
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = pe_i(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_i(:,:,jk,:) * z1_ai(:,:,:)
         ENDDO
         !== Snw heat content ==!
         DO jk = 1, nlay_s
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = pe_s(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_s(:,:,jk,:) * z1_ai(:,:,:)
         ENDDO
         !== Ice salt content ==!
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; zamsk(jm) = 0._wp
               zvar (:,:,:,jm) = pszv_i(:,:,jk,:)
               zhvar(:,:,:,jm) = pszv_i(:,:,jk,:) * z1_ai(:,:,:)
            ENDDO
         ELSE
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = psv_i(:,:,:)
            zhvar(:,:,:,jm) = psv_i(:,:,:) * z1_ai(:,:,:)
         ENDIF

         !== advection ==!
         CALL adv_umx( zamsk(1:jm), kn_umx, kt, zdt, pUu, pVv, zua_ho, zva_ho, zcu_box, zcv_box, &
            &                                            zhvar(:,:,:,1:jm), zvar(:,:,:,1:jm), zua_ups, zva_ups )
         !
         !== Recover quantities ==!
         jm = 1       ;    pv_i  (:,:,:)    = zvar (:,:,:,jm)
         jm = jm+1    ;    pv_s  (:,:,:)    = zvar (:,:,:,jm)
         DO jk = 1, nlay_i
            jm = jm+1 ;    pe_i  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         DO jk = 1, nlay_s
            jm = jm+1 ;    pe_s  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; pszv_i(:,:,jk,:) = zvar (:,:,:,jm)
            ENDDO
         ELSE
            jm = jm+1 ;    psv_i (:,:,:)    = zvar (:,:,:,jm)
         ENDIF

         !                          ! ----------------------------------------- !
      ELSEIF( np_advS == 3 ) THEN   ! -- advection form: -div( uV * uS / u ) -- !
         !                          ! ----------------------------------------- !
         !
         ! inverse of Vi
         WHERE( pv_i(:,:,:) >= epsi20 )
            z1_vi(:,:,:) = 1._wp / pv_i(:,:,:)
         ELSEWHERE
            z1_vi(:,:,:) = 0.
         END WHERE
         ! inverse of Vs
         WHERE( pv_s(:,:,:) >= epsi20 )
            z1_vs(:,:,:) = 1._wp / pv_s(:,:,:)
         ELSEWHERE
            z1_vs(:,:,:) = 0.
         END WHERE
         !
         ! It is important to first calculate the ice fields and then the snow fields (because we use the same arrays)
         !
         !== Ice volume ==!
         jm = 1 ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_i(:,:,:)
         zhvar(:,:,:,jm) = pv_i(:,:,:) * z1_ai(:,:,:)
         zuv_ups = zua_ups
         zvv_ups = zva_ups
         CALL adv_umx( zamsk(1:1), kn_umx, kt, zdt, pUu, pVv, zua_ho, zva_ho, zcu_box, zcv_box, &
            &                                           zhvar(:,:,:,1:1), zvar(:,:,:,1:1), zuv_ups, zvv_ups, zuv_ho, zvv_ho )
         !== Ice heat content ==!
         DO jk = 1, nlay_i
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = pe_i(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_i(:,:,jk,:) * z1_vi(:,:,:)
         ENDDO
         !== Ice salt content ==!
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; zamsk(jm) = 0._wp
               zvar (:,:,:,jm) = pszv_i(:,:,jk,:)
               zhvar(:,:,:,jm) = pszv_i(:,:,jk,:) * z1_vi(:,:,:)
            ENDDO
         ELSE
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = psv_i(:,:,:)
            zhvar(:,:,:,jm) = psv_i(:,:,:) * z1_vi(:,:,:)
         ENDIF
         CALL adv_umx( zamsk(2:jm), kn_umx, kt, zdt, pUu, pVv, zuv_ho, zvv_ho, zcu_box, zcv_box, &
            &                                            zhvar(:,:,:,2:jm), zvar(:,:,:,2:jm), zuv_ups, zvv_ups )
         !
         !== Recover quantities ==!
         jm = 1       ;    pv_i  (:,:,:)    = zvar (:,:,:,jm)
         DO jk = 1, nlay_i
            jm = jm+1 ;    pe_i  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         IF( nn_icesal == 4 ) THEN
            DO jk = 1, nlay_i
               jm = jm+1 ; pszv_i(:,:,jk,:) = zvar (:,:,:,jm)
            ENDDO
         ELSE
            jm = jm+1 ;    psv_i (:,:,:)    = zvar (:,:,:,jm)
         ENDIF

         !== Snw volume ==!
         jm = 1 ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_s(:,:,:)
         zhvar(:,:,:,jm) = pv_s(:,:,:) * z1_ai(:,:,:)
         zuv_ups = zua_ups
         zvv_ups = zva_ups
         CALL adv_umx( zamsk(1:1), kn_umx, kt, zdt, pUu , pVv, zua_ho , zva_ho , zcu_box, zcv_box, &
            &                                           zhvar(:,:,:,1:1), zvar(:,:,:,1:1), zuv_ups, zvv_ups, zuv_ho , zvv_ho )
         !== Snw heat content ==!
         DO jk = 1, nlay_s
            jm = jm+1 ; zamsk(jm) = 0._wp
            zvar (:,:,:,jm) = pe_s(:,:,jk,:)
            zhvar(:,:,:,jm) = pe_s(:,:,jk,:) * z1_vs(:,:,:)
         ENDDO
         CALL adv_umx( zamsk(2:jm), kn_umx, kt, zdt, pUu, pVv, zuv_ho, zvv_ho, zcu_box, zcv_box, &
            &                                            zhvar(:,:,:,2:jm), zvar(:,:,:,2:jm), zuv_ups, zvv_ups )
         !
         !== Recover quantities ==!
         jm = 1       ;    pv_s  (:,:,:)    = zvar (:,:,:,jm)
         DO jk = 1, nlay_s
            jm = jm+1 ;    pe_s  (:,:,jk,:) = zvar (:,:,:,jm)
         ENDDO
         !
         !
      ENDIF
      !
      !
      !== melt ponds ==!
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN

         !== pond volume ==!
         jm = 1    ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_ip(:,:,:)
         zhvar(:,:,:,jm) = pv_ip(:,:,:) * z1_aip(:,:,:)
         !== lid volume ==!
         jm = jm+1 ; zamsk(jm) = 0._wp
         zvar (:,:,:,jm) = pv_il(:,:,:)
         zhvar(:,:,:,jm) = pv_il(:,:,:) * z1_aip(:,:,:)
         !
         !== advection ==!
         CALL adv_umx( zamsk(1:jm), kn_umx, kt, zdt, pUu, pVv, zuap_ho, zvap_ho, zcu_box, zcv_box, &
            &                                            zhvar(:,:,:,1:jm), zvar(:,:,:,1:jm), zuap_ups, zvap_ups )

         !== Recover quantities ==!
         jm = 1       ;    pv_ip  (:,:,:)    = zvar (:,:,:,jm)
         jm = jm+1    ;    pv_il  (:,:,:)    = zvar (:,:,:,jm)
      ENDIF

      ! --- diagnostics --- !
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            diag_adv_mass(ji,jj) = diag_adv_mass(ji,jj) + (   SUM( pv_i (ji,jj,:) ) * rhoi + SUM( pv_s (ji,jj,:) ) * rhos &
               &                                            + SUM( pv_ip(ji,jj,:) ) * rhow + SUM( pv_il(ji,jj,:) ) * rhow &
               &                                          - zdiag_adv_mass(ji,jj) ) * r1_Dt_ice
            diag_adv_heat(ji,jj) = diag_adv_heat(ji,jj) + ( - SUM(SUM( pe_i(ji,jj,1:nlay_i,:) , dim=2 ) ) &
               &                                            - SUM(SUM( pe_s(ji,jj,1:nlay_s,:) , dim=2 ) ) &
               &                                          - zdiag_adv_heat(ji,jj) ) * r1_Dt_ice
         END DO
      END DO
      IF( nn_icesal == 4 ) THEN
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               diag_adv_salt(ji,jj) = diag_adv_salt(ji,jj) + ( SUM( SUM( pszv_i(ji,jj,:,:), dim=2 ) ) * rhoi &
                  &                                          - zdiag_adv_salt(ji,jj) ) * r1_Dt_ice
            END DO
         END DO
      ELSE
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               diag_adv_salt(ji,jj) = diag_adv_salt(ji,jj) + ( SUM( psv_i(ji,jj,:) ) * rhoi &
                  &                                          - zdiag_adv_salt(ji,jj) ) * r1_Dt_ice
            END DO
         END DO
      ENDIF

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
         !
         ! --- Make sure ice thickness is not too big --- !
         !     (because ice thickness can be too large where ice concentration is very small)
         IF( ln_icethd ) THEN
            !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   CALL Hbig( zdt, shi_max, shs_max, zhip_max, ssi_max, ses_max, sei_max, pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i, pe_s, pe_i )
            !ELSE
            IF( nn_icesal == 4 ) THEN
               CALL Hbig( zdt, shi_max, shs_max, sszi_max, ses_max, sei_max, pv_i, pv_s, pa_i, pszv_i, pe_s, pe_i )
            ELSE
               CALL Hbig( zdt, shi_max, shs_max,  ssi_max, ses_max, sei_max, pv_i, pv_s, pa_i,  psv_i, pe_s, pe_i )
            ENDIF
            !ENDIF
         ELSE
            CALL    Hbig( zdt, shi_max, shs_max,                            pv_i, pv_s, pa_i )
         ENDIF

         ! --- Ensure snow load is not too big --- !
         IF( ln_icethd ) THEN
            !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            !   CALL Hsnow( zdt, pv_i, pv_s, pa_i, pa_ip, pe_s )
            !ELSE
            CALL Hsnow( zdt, pv_i, pv_s, pa_i,        pe_s )
            !ENDIF
         ENDIF

      ENDIF !IF( .NOT. ln_pureADV2D )
      !
      !== open water area ==!
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            zati2 = SUM( pa_i(ji,jj,:) )
            pato_i(ji,jj) = MAX( 0._wp, pato_i(ji,jj) - ( zati2 - sati1(ji,jj) )            &   ! derive open water from ice concentration
               &                                      - (   ( pUu(ji,jj) - pUu(ji-1,jj) ) &   ! ad () for NP repro
               &                                          + ( pVv(ji,jj) - pVv(ji,jj-1) ) ) * r1_e1e2t(ji,jj) * zdt )
         END DO
      END DO
      ! note: no need of lbc_lnk for open water (never used in the halos)
      !
      !
      !
      ! --- Deallocate arrays --- !
      DEALLOCATE( zvar, zhvar, zamsk )
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) DEALLOCATE( z1_aip )
      IF( np_advS == 3 )                DEALLOCATE( zuv_ho, zvv_ho, zuv_ups, zvv_ups, z1_vi, z1_vs )
      !
   END SUBROUTINE ice_dyn_adv_umx


   SUBROUTINE adv_umx_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE adv_umx_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for Prather advection
      !!-------------------------------------------------------------------
      INTEGER, DIMENSION(11) ::   ierr
      INTEGER :: k_alloc
      CHARACTER(len=4), PARAMETER :: crtnm = 'adv_umx_init'
      !!-------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( sati1(jpi,jpj),   STAT = ierr(1) )
      !%acc enter data copyin( sati1 )
      IF( .NOT. ln_pureADV2D ) THEN
         ALLOCATE( shi_max(jpi,jpj,jpl), shs_max(jpi,jpj,jpl),    STAT = ierr(2) )
         !%acc enter data copyin( shi_max, shs_max )
      ENDIF
      !
      IF( ln_icethd ) THEN
         IF( .NOT. ln_pureADV2D ) THEN
            IF( nn_icesal == 4 ) THEN
               ALLOCATE( sei_max(jpi,jpj,nlay_i,jpl), ses_max(jpi,jpj,nlay_s,jpl), sszi_max(jpi,jpj,nlay_i,jpl), STAT = ierr(3) )
               !%acc enter data copyin( sei_max, ses_max, sszi_max )
            ELSE
               ALLOCATE( sei_max(jpi,jpj,nlay_i,jpl), ses_max(jpi,jpj,nlay_s,jpl), ssi_max(jpi,jpj,jpl), STAT = ierr(4) )
               !%acc enter data copyin( sei_max, ses_max, ssi_max )
            ENDIF
         ENDIF
      ENDIF
      !
      k_alloc = MAXVAL( ierr(:) )
      CALL mpp_sum ( crtnm, k_alloc )
      IF( k_alloc > 0 ) CALL ctl_stop('STOP', crtnm//': unable to allocate ice arrays for Prather advection scheme')
      !
   END SUBROUTINE adv_umx_init


   SUBROUTINE adv_umx( pamsk, kn_umx, kt, pdt, pu, pv, puc, pvc, pubox, pvbox,  &
      &                                            pt, ptc, pua_ups, pva_ups, pua_ho, pva_ho )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE adv_umx  ***
      !!
      !! **  Purpose :   Compute the now trend due to total advection of
      !!                 tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   - calculate upstream fluxes and upstream solution for tracers V/A(=H) etc
      !!                 - calculate tracer H at u and v points (Ultimate)
      !!                 - calculate the high order fluxes using alterning directions (Macho)
      !!                 - apply a limiter on the fluxes (nonosc_ice)
      !!                 - convert this tracer flux to a "volume" flux (uH -> uV)
      !!                 - apply a limiter a second time on the volumes fluxes (nonosc_ice)
      !!                 - calculate the high order solution for V
      !!
      !! ** Action : solve 3 equations => a) dA/dt  = -div(uA)
      !!                                  b) dV/dt  = -div(uV)  using dH/dt = -u.grad(H)
      !!                                  c) dVS/dt = -div(uVS) using either dHS/dt = -u.grad(HS) or dS/dt = -u.grad(S)
      !!
      !!             in eq. b), - fluxes uH are evaluated (with UMx) and limited with nonosc_ice. This step is necessary to get a good H.
      !!                        - then we convert this flux to a "volume" flux this way => uH * uA / u
      !!                             where uA is the flux from eq. a)
      !!                             this "volume" flux is also limited with nonosc_ice (otherwise overshoots can occur)
      !!                        - at last we estimate dV/dt = -div(uH * uA / u)
      !!
      !!             in eq. c), one can solve the equation for  S (ln_advS=T), then dVS/dt = -div(uV * uS  / u)
      !!                                                or for HS (ln_advS=F), then dVS/dt = -div(uA * uHS / u)
      !!
      !! ** Note : - this method can lead to tiny negative V (-1.e-20) => set it to 0 while conserving mass etc.
      !!           - At the ice edge, Ultimate scheme can lead to:
      !!                              1) negative interpolated tracers at u-v points
      !!                              2) non-zero interpolated tracers at u-v points eventhough there is no ice and velocity is outward
      !!                              Solution for 1): apply an upstream scheme when it occurs. A better solution would be to degrade the order of
      !!                                               the scheme automatically by applying a mask of the ice cover inside Ultimate (not done).
      !!                              Solution for 2): we set it to 0 in this case
      !!           - Eventhough 1D tests give very good results (typically the one from Schar & Smolarkiewiecz), the 2D is less good.
      !!             Large values of H can appear for very small ice concentration, and when it does it messes the things up since we
      !!             work on H (and not V). It is partly related to the multi-category approach
      !!             Therefore, after advection we limit the thickness to the largest value of the 9-points around (only if ice
      !!             concentration is small). We also limit S and T.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)          , INTENT(in   )           ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   )           ::   kn_umx           ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                         , INTENT(in   )           ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   )           ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   )           ::   pu   , pv        ! 2 ice velocity components => u*e2
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   )           ::   puc  , pvc       ! 2 ice velocity components => u*e2 or u*a*e2u
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   )           ::   pubox, pvbox     ! upstream velocity
      REAL(wp), DIMENSION(:,:,:,:)    , INTENT(inout)           ::   pt               ! tracer field
      REAL(wp), DIMENSION(:,:,:,:)    , INTENT(inout)           ::   ptc              ! tracer content field
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout), OPTIONAL ::   pua_ups, pva_ups ! upstream u*a fluxes
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out), OPTIONAL ::   pua_ho, pva_ho   ! high order u*a fluxes
      !
      INTEGER  ::   ji, jj, jl, jm   ! dummy loop indices
      INTEGER  ::   ndim             ! number of variables to advect
      REAL(wp) ::   ztra             ! local scalar
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zfu_ho , zfv_ho
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zfu_ups, zfv_ups, zt_ups
      !!----------------------------------------------------------------------
      !
      ndim = SIZE( ptc, dim=4 )
      !
      ALLOCATE( zfu_ho (jpi,jpj,jpl,ndim), zfv_ho (jpi,jpj,jpl,ndim), &
         &      zfu_ups(jpi,jpj,jpl,ndim), zfv_ups(jpi,jpj,jpl,ndim), zt_ups(jpi,jpj,jpl,ndim)  )
      !
      ! Upstream (_ups) fluxes
      ! -----------------------
      DO jm = 1, ndim
         CALL upstream( pamsk(jm), kt, pdt, pt(:,:,:,jm), pu(:,:), pv(:,:), &                        ! <<= in
            &                                   zt_ups(:,:,:,jm), zfu_ups(:,:,:,jm), zfv_ups(:,:,:,jm) ) ! =>> out ( gives zt_ups(1,1,1,1), zfu_ups(2,1,1,1) & zfv_ups(1,1,2,1) )
      ENDDO
      !
      ! High order (_ho) fluxes
      ! -----------------------
      SELECT CASE( kn_umx )
         !
      CASE ( 20 )                          !== centered second order ==!
         DO jm = 1, ndim
            CALL cen2( pamsk(jm), kt, pdt, pt(:,:,:,jm), pu(:,:), pv(:,:), &  ! <<= in
               &                         zfu_ups(:,:,:,jm), zfv_ups(:,:,:,jm), &  ! <<= in  (upstream)
               &                         zfu_ho (:,:,:,jm), zfv_ho (:,:,:,jm)  )  ! =>> out (high order) ( gives zfu_ho(2,1,1,1) & zfv_ho(1,1,2,1) )
         ENDDO
      CASE ( 1:5 )                         !== 1st to 5th order ULTIMATE-MACHO scheme ==!
         CALL macho( pamsk(:), kn_umx, kt, pdt, pt(:,:,:,:), pu(:,:), pv(:,:), pubox(:,:), pvbox(:,:), &  ! <<= in
            &                                                          zfu_ups(:,:,:,:), zfv_ups(:,:,:,:), &  ! <<= in  (upstream)
            &                                                          zfu_ho (:,:,:,:), zfv_ho (:,:,:,:)  )  ! =>> out (high order) ( gives zfu_ho(2,1,1,1) & zfv_ho(1,1,2,1) )
      END SELECT
      !
      ! Flux limiter
      ! ------------
      IF( np_limiter == 1 ) THEN
         CALL lbc_lnk( 'icedyn_adv_umx', zt_ups, 'T', 1.0_wp )   ! nonosc needs zt_ups over the whole domain
         DO jm = 1, ndim
            CALL nonosc_ice( pamsk(jm), pdt, pu, pv, pt(:,:,:,jm), zt_ups(:,:,:,jm), zfu_ups(:,:,:,jm), zfv_ups(:,:,:,jm), &
               &                                                                     zfu_ho (:,:,:,jm), zfv_ho (:,:,:,jm)  ) ! gives zfu_ho(1,0,0,0) & zfv_ho(0,0,1,0)
         ENDDO
         !
      ENDIF
      !              --ho    --ho
      ! new fluxes = u*H  *  u*a / u
      ! ----------------------------
      DO jm = 1, ndim
         IF( pamsk(jm) == 0._wp ) THEN

            DO jl = 1, jpl
               DO jj=Njs0-0, Nje0+0
                  DO ji=Nis0-1, Nie0+0
                     IF( ABS( pu(ji,jj) ) > epsi10 ) THEN
                        zfu_ho(ji,jj,jl,jm) = zfu_ho(ji,jj,jl,jm) * puc(ji,jj,jl) / pu(ji,jj)
                     ELSE
                        zfu_ho(ji,jj,jl,jm) = 0._wp
                     ENDIF
                  END DO
               END DO
               DO jj=Njs0-1, Nje0+0
                  DO ji=Nis0-0, Nie0+0
                     IF( ABS( pv(ji,jj) ) > epsi10 ) THEN
                        zfv_ho(ji,jj,jl,jm) = zfv_ho(ji,jj,jl,jm) * pvc(ji,jj,jl) / pv(ji,jj)
                     ELSE
                        zfv_ho(ji,jj,jl,jm) = 0._wp
                     ENDIF
                  END DO
               END DO
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+1
                     IF( ABS( pu(ji,jj) ) > epsi10 ) THEN
                        zfu_ups(ji,jj,jl,jm) = zfu_ups(ji,jj,jl,jm) * pua_ups(ji,jj,jl) / pu(ji,jj)
                     ELSE
                        zfu_ups(ji,jj,jl,jm) = 0._wp
                     ENDIF
                  END DO
               END DO
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-1, Nie0+1
                     IF( ABS( pv(ji,jj) ) > epsi10 ) THEN
                        zfv_ups(ji,jj,jl,jm) = zfv_ups(ji,jj,jl,jm) * pva_ups(ji,jj,jl) / pv(ji,jj)
                     ELSE
                        zfv_ups(ji,jj,jl,jm) = 0._wp
                     ENDIF
                  END DO
               END DO

               ! the new "volume" fluxes must also be "flux corrected" => we calculate the upstream solution and apply a limiter again
               DO jj=Njs0, Nje0
                  DO ji=Nis0, Nie0
                     ztra = - (  ( zfu_ups(ji,jj,jl,jm) - zfu_ups(ji-1,jj,jl,jm) ) &   ! add () for NP repro
                        &      + ( zfv_ups(ji,jj,jl,jm) - zfv_ups(ji,jj-1,jl,jm) ) )
                     !
                     zt_ups(ji,jj,jl,jm) = ( ptc(ji,jj,jl,jm) + ztra * r1_e1e2t(ji,jj) * pdt ) * xmskt(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
      ENDDO
      ! lbc needed for nonosc
      IF( MINVAL( pamsk(:) ) == 0._wp .AND. np_limiter == 1 )  &
         &                   CALL lbc_lnk( 'icedyn_adv_umx', zt_ups, 'T', 1.0_wp, zfu_ho, 'U', -1.0_wp, zfv_ho, 'V', -1.0_wp  )
      ! flux limiter
      DO jm = 1, ndim
         IF( pamsk(jm) == 0._wp ) THEN
            IF    ( np_limiter == 1 ) THEN
               CALL nonosc_ice( 1._wp, pdt, pu, pv, ptc(:,:,:,jm), zt_ups (:,:,:,jm), zfu_ups(:,:,:,jm), zfv_ups(:,:,:,jm), &
                  &                                                                   zfu_ho (:,:,:,jm), zfv_ho (:,:,:,jm)  ) ! gives zfu_ho(1,0,0,0) & zfv_ho(0,0,1,0)
            ELSEIF( np_limiter == 2 .OR. np_limiter == 3 ) THEN
               CALL limiter_x( pdt, pu, ptc(:,:,:,jm), zfu_ups(:,:,:,jm), zfu_ho(:,:,:,jm) )
               CALL limiter_y( pdt, pv, ptc(:,:,:,jm), zfv_ups(:,:,:,jm), zfv_ho(:,:,:,jm) )

            ENDIF
         ENDIF
      ENDDO
      !
      !                                   --ho    --ups
      ! in case of advection of A: output u*a and u*a
      ! -----------------------------------------------
      IF( PRESENT( pua_ho ) ) THEN
         DO jl = 1, jpl
            DO jj=Njs0-0, Nje0+0
               DO ji=Nis0-1, Nie0+0
                  pua_ho (ji,jj,jl) = zfu_ho (ji,jj,jl,1)
               END DO
            END DO
            DO jj=Njs0-1, Nje0+0
               DO ji=Nis0-0, Nie0+0
                  pva_ho (ji,jj,jl) = zfv_ho (ji,jj,jl,1)
               END DO
            END DO
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-2, Nie0+1
                  pua_ups(ji,jj,jl) = zfu_ups(ji,jj,jl,1)
               END DO
            END DO
            DO jj=Njs0-2, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  pva_ups(ji,jj,jl) = zfv_ups(ji,jj,jl,1)
               END DO
            END DO
         END DO
      ENDIF
      !
      ! final trend with corrected fluxes
      ! ---------------------------------
      DO jm = 1, ndim
         DO jl = 1, jpl
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  ztra = - (  ( zfu_ho(ji,jj,jl,jm) - zfu_ho(ji-1,jj,jl,jm) ) &   ! add () for NP repro
                     &      + ( zfv_ho(ji,jj,jl,jm) - zfv_ho(ji,jj-1,jl,jm) ) )
                  !
                  ptc(ji,jj,jl,jm) = ( ptc(ji,jj,jl,jm) + ztra * r1_e1e2t(ji,jj) * pdt ) * xmskt(ji,jj)
               END DO
            END DO
         END DO
      ENDDO
      !
      DEALLOCATE( zfu_ho, zfv_ho, zfu_ups, zfv_ups, zt_ups )
      !
   END SUBROUTINE adv_umx


   SUBROUTINE upstream( pamsk, kt, pdt, pt, pu, pv, pt_ups, pfu_ups, pfv_ups )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE upstream  ***
      !!
      !! **  Purpose :   compute the upstream fluxes and upstream guess of tracer
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pt_ups           ! upstream guess of tracer
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ups, pfv_ups ! upstream fluxes
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   ztra          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zpt
      !!----------------------------------------------------------------------

      IF( .NOT. ll_upsxy ) THEN         !** no alternate directions **!
         !
         DO jl = 1, jpl
            DO jj=Njs0-2, Nje0+1
               DO ji=Nis0-2, Nie0+1
                  pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                  pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
               END DO
            END DO
         END DO
         !
      ELSE                              !** alternate directions **!
         !
         IF( MOD( (kt - 1) , 2 ) == 0 ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj=Njs0-2, Nje0+2
                  DO ji=Nis0-2, Nie0+1
                     pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pu(ji,jj), 0._wp ) * pt(ji+1,jj,jl)
                  END DO
               END DO
               !                        !-- first guess of tracer from u-flux
               DO jj=Njs0-2, Nje0+2
                  DO ji=Nis0-1, Nie0+1
                     ztra = - ( pfu_ups(ji,jj,jl) - pfu_ups(ji-1,jj,jl) )              &
                        &   + ( pu     (ji,jj   ) - pu     (ji-1,jj   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * xmskt(ji,jj)
                  END DO
               END DO
               !                        !-- flux in y-direction
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-1, Nie0+1
                     pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( pv(ji,jj), 0._wp ) * zpt(ji,jj+1)
                  END DO
               END DO
            END DO
            !
         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            DO jl = 1, jpl              !-- flux in y-direction
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-2, Nie0+2
                     pfv_ups(ji,jj,jl) = MAX( pv(ji,jj), 0._wp ) * pt(ji,jj,jl) + MIN( pv(ji,jj), 0._wp ) * pt(ji,jj+1,jl)
                  END DO
               END DO
               !                        !-- first guess of tracer from v-flux
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+2
                     ztra = - ( pfv_ups(ji,jj,jl) - pfv_ups(ji,jj-1,jl) )  &
                        &   + ( pv     (ji,jj   ) - pv     (ji,jj-1   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * xmskt(ji,jj)
                  END DO
               END DO
               !                        !-- flux in x-direction
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+1
                     pfu_ups(ji,jj,jl) = MAX( pu(ji,jj), 0._wp ) * zpt(ji,jj) + MIN( pu(ji,jj), 0._wp ) * zpt(ji+1,jj)
                  END DO
               END DO
            END DO
            !
         ENDIF

      ENDIF
      !
      DO jl = 1, jpl                    !-- after tracer with upstream scheme
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ztra = - (   ( pfu_ups(ji,jj,jl) - pfu_ups(ji-1,jj  ,jl) )   &   ! add () for NP repro
                  &       + ( pfv_ups(ji,jj,jl) - pfv_ups(ji  ,jj-1,jl) ) ) &
                  &   + (   ( pu     (ji,jj   ) - pu     (ji-1,jj     ) )   &
                  &       + ( pv     (ji,jj   ) - pv     (ji  ,jj-1   ) ) ) * pt(ji,jj,jl) * (1.-pamsk)
               !
               pt_ups(ji,jj,jl) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * xmskt(ji,jj)
            END DO
         END DO
      END DO

   END SUBROUTINE upstream


   SUBROUTINE cen2( pamsk, kt, pdt, pt, pu, pv, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE cen2  ***
      !!
      !! **  Purpose :   compute the high order fluxes using a centered
      !!                 second order scheme
      !!----------------------------------------------------------------------
      REAL(wp)                        , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                         , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                        , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:  )      , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:,:)      , INTENT(in   ) ::   pfu_ups, pfv_ups ! upstream fluxes
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(  out) ::   pfu_ho, pfv_ho   ! high order fluxes
      !
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp) ::   ztra          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zpt
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ll_hoxy ) THEN           !** no alternate directions **!
         !
         DO jl = 1, jpl
            DO jj=Njs0-2, Nje0+2
               DO ji=Nis0-2, Nie0+1
                  pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj  ,jl) )
               END DO
            END DO
            DO jj=Njs0-2, Nje0+1
               DO ji=Nis0-2, Nie0+2
                  pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji  ,jj+1,jl) )
               END DO
            END DO
         END DO
         !
         IF( np_limiter == 2 .OR. np_limiter == 3 ) THEN
            CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
            CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
         ENDIF
         !
      ELSE                              !** alternate directions **!
         !
         IF( MOD( (kt - 1) , 2 ) == 0 ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !
            DO jl = 1, jpl              !-- flux in x-direction
               DO jj=Njs0-2, Nje0+2
                  DO ji=Nis0-2, Nie0+1
                     pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( pt(ji,jj,jl) + pt(ji+1,jj,jl) )
                  END DO
               END DO
            END DO
            IF( np_limiter == 2 .OR. np_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )

            DO jl = 1, jpl              !-- first guess of tracer from u-flux
               DO jj=Njs0-2, Nje0+2
                  DO ji=Nis0-1, Nie0+1
                     ztra = - ( pfu_ho(ji,jj,jl) - pfu_ho(ji-1,jj,jl) )              &
                        &   + ( pu    (ji,jj   ) - pu    (ji-1,jj   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * xmskt(ji,jj)
                  END DO
               END DO

               !                        !-- flux in y-direction
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-1, Nie0+1
                     pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( zpt(ji,jj) + zpt(ji,jj+1) )
                  END DO
               END DO
            END DO
            IF( np_limiter == 2 .OR. np_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )

         ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
            !
            DO jl = 1, jpl              !-- flux in y-direction
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-2, Nie0+2
                     pfv_ho(ji,jj,jl) = 0.5_wp * pv(ji,jj) * ( pt(ji,jj,jl) + pt(ji,jj+1,jl) )
                  END DO
               END DO
            END DO
            IF( np_limiter == 2 .OR. np_limiter == 3 )   CALL limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
            !
            DO jl = 1, jpl              !-- first guess of tracer from v-flux
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+2
                     ztra = - ( pfv_ho(ji,jj,jl) - pfv_ho(ji,jj-1,jl) )  &
                        &   + ( pv    (ji,jj   ) - pv    (ji,jj-1   ) ) * pt(ji,jj,jl) * (1.-pamsk)
                     !
                     zpt(ji,jj) = ( pt(ji,jj,jl) + ztra * pdt * r1_e1e2t(ji,jj) ) * xmskt(ji,jj)
                  END DO
               END DO
               !
               !                        !-- flux in x-direction
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+1
                     pfu_ho(ji,jj,jl) = 0.5_wp * pu(ji,jj) * ( zpt(ji,jj) + zpt(ji+1,jj) )
                  END DO
               END DO
            END DO
            IF( np_limiter == 2 .OR. np_limiter == 3 )   CALL limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )

         ENDIF

      ENDIF
      !
      !
   END SUBROUTINE cen2


   SUBROUTINE macho( pamsk, kn_umx, kt, pdt, pt, pu, pv, pubox, pvbox, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE macho  ***
      !!
      !! **  Purpose :   compute the high order fluxes using Ultimate-Macho scheme
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)      , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   ) ::   kn_umx           ! order of the scheme (1-5=UM or 20=CEN2)
      INTEGER                     , INTENT(in   ) ::   kt               ! number of iteration
      REAL(wp)                    , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pt               ! tracer fields
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pu, pv           ! 2 ice velocity components
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pubox, pvbox     ! upstream velocity
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pfu_ups, pfv_ups ! upstream fluxes
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pfu_ho, pfv_ho   ! high order fluxes (only out)
      !
      INTEGER  ::   ji, jj, jl, jm    ! dummy loop indices
      INTEGER  ::   ndim              ! number of variables to advect
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zt_u, zt_v, zpt
      !!----------------------------------------------------------------------
      ndim = SIZE( pt, dim=4 )
      !
      ALLOCATE( zt_u(jpi,jpj,jpl,ndim), zt_v(jpi,jpj,jpl,ndim), zpt(jpi,jpj,jpl,ndim) )

      !
      IF( MOD( (kt - 1) , 2 ) == 0 ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
         !
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         CALL ultimate_x( 2, pamsk, kn_umx, pdt, pt, pu, zt_u, pfu_ho )
         !
         !                                                        !--  limiter in x --!
         DO jm = 1, ndim
            IF( np_limiter == 2 .OR. np_limiter == 3 ) CALL limiter_x( pdt, pu, pt(:,:,:,jm), pfu_ups(:,:,:,jm), pfu_ho(:,:,:,jm) )
         END DO
         !                                                        !--  advective form update in zpt  --!
         DO jm = 1, ndim
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+2
                  DO ji=Nis0-1, Nie0+1
                     zpt(ji,jj,jl,jm) = (pt(ji,jj,jl,jm) - ( pubox(ji,jj   ) * ( zt_u(ji,jj,jl,jm) - zt_u(ji-1,jj,jl,jm) )       &
                        &                                                    * r1_e1t(ji,jj)                                     &
                        &                                  + pt(ji,jj,jl,jm) * ( pu  (ji,jj)       - pu  (ji-1,jj) ) * pamsk(jm) &
                        &                                                    * r1_e1e2t(ji,jj)                                   &
                        &                                  ) * pdt) * xmskt(ji,jj)
                  END DO
               END DO
            END DO
         END DO
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_y( 1, pamsk, kn_umx, pdt, zpt, pv, zt_v, pfv_ho )
         ELSE
            CALL ultimate_y( 1, pamsk, kn_umx, pdt,  pt, pv, zt_v, pfv_ho )
         ENDIF
         !                                                        !--  limiter in y --!
         DO jm = 1, ndim
            IF( np_limiter == 2 .OR. np_limiter == 3 ) CALL limiter_y( pdt, pv, pt(:,:,:,jm), pfv_ups(:,:,:,jm), pfv_ho(:,:,:,jm) )
         END DO
         !
         !
      ELSE                                                               !==  even ice time step:  adv_y then adv_x  ==!
         !
         !                                                        !--  ultimate interpolation of pt at v-point  --!
         CALL ultimate_y( 2, pamsk, kn_umx, pdt, pt, pv, zt_v, pfv_ho )
         !
         !                                                        !--  limiter in y --!
         DO jm = 1, ndim
            IF( np_limiter == 2 .OR. np_limiter == 3 ) CALL limiter_y( pdt, pv, pt(:,:,:,jm), pfv_ups(:,:,:,jm), pfv_ho(:,:,:,jm) )
         END DO
         !                                                        !--  advective form update in zpt  --!
         DO jm = 1, ndim
            DO jl = 1, jpl
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-2, Nie0+2
                     zpt(ji,jj,jl,jm) = (pt(ji,jj,jl,jm) - ( pvbox(ji,jj   ) * ( zt_v(ji,jj,jl,jm) - zt_v(ji,jj-1,jl,jm) )       &
                        &                                                    * r1_e2t(ji,jj)                                     &
                        &                                  + pt(ji,jj,jl,jm) * ( pv  (ji,jj)       - pv  (ji,jj-1) ) * pamsk(jm) &
                        &                                                    * r1_e1e2t(ji,jj)                                   &
                        &                                  ) * pdt) * xmskt(ji,jj)
                  END DO
               END DO
            END DO
         END DO
         !                                                        !--  ultimate interpolation of pt at u-point  --!
         IF( ll_hoxy ) THEN
            CALL ultimate_x( 1, pamsk, kn_umx, pdt, zpt, pu, zt_u, pfu_ho )
         ELSE
            CALL ultimate_x( 1, pamsk, kn_umx, pdt,  pt, pu, zt_u, pfu_ho )
         ENDIF
         !                                                        !--  limiter in x --!
         DO jm = 1, ndim
            IF( np_limiter == 2 .OR. np_limiter == 3 ) CALL limiter_x( pdt, pu, pt(:,:,:,jm), pfu_ups(:,:,:,jm), pfu_ho(:,:,:,jm) )
         END DO
         !
      ENDIF

      DEALLOCATE( zt_u, zt_v, zpt )

      !
   END SUBROUTINE macho


   SUBROUTINE ultimate_x( kloop, pamsk, kn_umx, pdt, pt, pu, pt_u, pfu_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_x  ***
      !!
      !! **  Purpose :   compute tracer at u-points
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74.
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kloop     ! either 0 or nn_hls depending on the order of the call
      REAL(wp), DIMENSION(:)      , INTENT(in   ) ::   pamsk     ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pu        ! ice i-velocity component
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pt_u      ! tracer at u-point (only out)
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pfu_ho    ! high order flux   (only out)
      !
      INTEGER  ::   ji, jj, jl, jm   ! dummy loop indices
      INTEGER  ::   ndim             ! number of variables to advect
      REAL(wp) ::   zcu, zdx2, zdx4        !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::   ztu1, ztu3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ztu2, ztu4
      !!----------------------------------------------------------------------
      ndim = SIZE( pt, dim=4 )
      !
      IF( kn_umx >= 3 )    ALLOCATE( ztu1(jpi,jpj), ztu2(jpi,jpj,jpl,ndim) )
      IF( kn_umx == 5 )    ALLOCATE( ztu3(jpi,jpj), ztu4(jpi,jpj,jpl,ndim) )
      !
      DO jm = 1, ndim
         !                                                     !--  Laplacian in i-direction  --!
         IF( kn_umx >= 3 ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     ztu1(ji,jj) = ( pt(ji+1,jj,jl,jm) - pt(ji,jj,jl,jm) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
                  END DO
               END DO
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-1, Nie0+1
                     ztu2(ji,jj,jl,jm) = ( ztu1(ji,jj) - ztu1(ji-1,jj) ) * r1_e1t(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
         !                                                     !--  BiLaplacian in i-direction  --!
         IF( kn_umx == 5 ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-1, Nie0+0
                     ztu3(ji,jj) = ( ztu2(ji+1,jj,jl,jm) - ztu2(ji,jj,jl,jm) ) * r1_e1u(ji,jj) * umask(ji,jj,1)
                  END DO
               END DO
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-0, Nie0+0
                     ztu4(ji,jj,jl,jm) = ( ztu3(ji,jj) - ztu3(ji-1,jj) ) * r1_e1t(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
         !
      ENDDO
      ! lbc only needed for some orders
      IF    ( kn_umx == 3 .OR. kn_umx == 4 ) THEN
         CALL lbc_lnk( 'icedyn_adv_umx', ztu2, 'T', 1.0_wp )
      ELSEIF( kn_umx == 5 )                  THEN
         CALL lbc_lnk( 'icedyn_adv_umx', ztu2, 'T', 1.0_wp, ztu4, 'T', 1.0_wp )
      ENDIF
      !
      !
      DO jm = 1, ndim
         !
         SELECT CASE ( kn_umx )
            !
         CASE( 1 )                                                   !==  1st order central TIM  ==! (Eq. 21)
            !
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * (                          ( pt(ji+1,jj,jl,jm) + pt(ji,jj,jl,jm) ) &
                        &                                        - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl,jm) - pt(ji,jj,jl,jm) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 2 )                                                   !==  2nd order central TIM  ==! (Eq. 23)
            !
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                     pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * (                    ( pt(ji+1,jj,jl,jm) + pt(ji,jj,jl,jm) )  &
                        &                                                     - zcu   * ( pt(ji+1,jj,jl,jm) - pt(ji,jj,jl,jm) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 3 )                                                   !==  3rd order central TIM  ==! (Eq. 24)
            !
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                     zdx2 = e1u(ji,jj) * e1u(ji,jj)
                     !!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                     pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * ( (                  ( pt  (ji+1,jj,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcu   * ( pt  (ji+1,jj,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6 * zdx2 * ( zcu*zcu - 1._wp ) *      (                  ( ztu2(ji+1,jj,jl,jm) + ztu2(ji,jj,jl,jm) )   &
                        &                                        - SIGN( 1._wp, zcu ) * ( ztu2(ji+1,jj,jl,jm) - ztu2(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 4 )                                                   !==  4th order central TIM  ==! (Eq. 27)
            !
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                     zdx2 = e1u(ji,jj) * e1u(ji,jj)
                     !!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                     pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * ( (                  ( pt  (ji+1,jj,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcu   * ( pt  (ji+1,jj,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6 * zdx2 * ( zcu*zcu - 1._wp ) *      (                  ( ztu2(ji+1,jj,jl,jm) + ztu2(ji,jj,jl,jm) )   &
                        &                                            - 0.5_wp * zcu   * ( ztu2(ji+1,jj,jl,jm) - ztu2(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 5 )                                                   !==  5th order central TIM  ==! (Eq. 29)
            !
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     zcu  = pu(ji,jj) * r1_e2u(ji,jj) * pdt * r1_e1u(ji,jj)
                     zdx2 = e1u(ji,jj) * e1u(ji,jj)
                     !!rachid          zdx2 = e1u(ji,jj) * e1t(ji,jj)
                     zdx4 = zdx2 * zdx2
                     pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * ( (                  ( pt  (ji+1,jj,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcu   * ( pt  (ji+1,jj,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6   * zdx2 * ( zcu*zcu - 1._wp ) *    (                  ( ztu2(ji+1,jj,jl,jm) + ztu2(ji,jj,jl,jm) )   &
                        &                                            - 0.5_wp * zcu   * ( ztu2(ji+1,jj,jl,jm) - ztu2(ji,jj,jl,jm) ) ) &
                        & + r1_120 * zdx4 * ( zcu*zcu - 1._wp ) * ( zcu*zcu - 4._wp ) * ((ztu4(ji+1,jj,jl,jm) + ztu4(ji,jj,jl,jm) )   &
                        &                                        - SIGN( 1._wp, zcu ) * ( ztu4(ji+1,jj,jl,jm) - ztu4(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         END SELECT
         !
         !
         ! if pt at u-point is negative then use the upstream value
         !    this should not be necessary if a proper sea-ice mask is set in Ultimate
         !    to degrade the order of the scheme when necessary (for ex. at the ice edge)
         IF( ll_neg ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-kloop, Nje0+kloop
                  DO ji=Nis0-2, Nie0+1
                     IF( pt_u(ji,jj,jl,jm) < 0._wp .OR. ( imsk_small(ji,jj,jl) == 0 .AND. pamsk(jm) == 0. ) ) THEN
                        pt_u(ji,jj,jl,jm) = 0.5_wp * umask(ji,jj,1) * (                      ( pt(ji+1,jj,jl,jm) + pt(ji,jj,jl,jm) ) &
                           &                                    - SIGN( 1._wp, pu(ji,jj) ) * ( pt(ji+1,jj,jl,jm) - pt(ji,jj,jl,jm) ) )
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
         !                                                     !-- High order flux in i-direction  --!
         DO jl = 1, jpl
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-2, Nie0+1
                  pfu_ho(ji,jj,jl,jm) = pu(ji,jj) * pt_u(ji,jj,jl,jm)
               END DO
            END DO
         END DO
         !
      ENDDO

      IF( kn_umx >= 3 )    DEALLOCATE( ztu1, ztu2 )
      IF( kn_umx == 5 )    DEALLOCATE( ztu3, ztu4 )
      !
   END SUBROUTINE ultimate_x


   SUBROUTINE ultimate_y( kloop, pamsk, kn_umx, pdt, pt, pv, pt_v, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ultimate_y  ***
      !!
      !! **  Purpose :   compute tracer at v-points
      !!
      !! **  Method  :   ...
      !!
      !! Reference : Leonard, B.P., 1991, Comput. Methods Appl. Mech. Eng., 88, 17-74.
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kloop     ! either 0 or nn_hls depending on the order of the call
      REAL(wp), DIMENSION(:)      , INTENT(in   ) ::   pamsk     ! advection of concentration (1) or other tracers (0)
      INTEGER                     , INTENT(in   ) ::   kn_umx    ! order of the scheme (1-5=UM or 20=CEN2)
      REAL(wp)                    , INTENT(in   ) ::   pdt       ! tracer time-step
      REAL(wp), DIMENSION(:,:  )  , INTENT(in   ) ::   pv        ! ice j-velocity component
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pt        ! tracer fields
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pt_v      ! tracer at v-point (only out)
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pfv_ho    ! high order flux   (only out)
      !
      INTEGER  ::   ji, jj, jl, jm   ! dummy loop indices
      INTEGER  ::   ndim             ! number of variables to advect
      REAL(wp) ::   zcv, zdy2, zdy4    !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::   ztv1, ztv3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ztv2, ztv4
      !!----------------------------------------------------------------------
      ndim = SIZE( pt, dim=4 )
      !
      IF( kn_umx >= 3 )    ALLOCATE( ztv1(jpi,jpj), ztv2(jpi,jpj,jpl,ndim) )
      IF( kn_umx == 5 )    ALLOCATE( ztv3(jpi,jpj), ztv4(jpi,jpj,jpl,ndim) )
      !
      DO jm = 1, ndim
         !                                                     !--  Laplacian in j-direction  --!
         IF( kn_umx >= 3 ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     ztv1(ji,jj) = ( pt(ji,jj+1,jl,jm) - pt(ji,jj,jl,jm) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
                  END DO
               END DO
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     ztv2(ji,jj,jl,jm) = ( ztv1(ji,jj) - ztv1(ji,jj-1) ) * r1_e2t(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
         !                                                     !--  BiLaplacian in j-direction  --!
         IF( kn_umx == 5 ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-1, Nje0+0
                  DO ji=Nis0-kloop, Nie0+kloop
                     ztv3(ji,jj) = ( ztv2(ji,jj+1,jl,jm) - ztv2(ji,jj,jl,jm) ) * r1_e2v(ji,jj) * vmask(ji,jj,1)
                  END DO
               END DO
               DO jj=Njs0-0, Nje0+0
                  DO ji=Nis0-kloop, Nie0+kloop
                     ztv4(ji,jj,jl,jm) = ( ztv3(ji,jj) - ztv3(ji,jj-1) ) * r1_e2t(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
         !
      ENDDO
      ! lbc only needed for some orders
      IF    ( kn_umx == 3 .OR. kn_umx == 4 ) THEN
         CALL lbc_lnk( 'icedyn_adv_umx', ztv2, 'T', 1.0_wp )
      ELSEIF( kn_umx == 5 )                  THEN
         CALL lbc_lnk( 'icedyn_adv_umx', ztv2, 'T', 1.0_wp, ztv4, 'T', 1.0_wp )
      ENDIF
      !
      !
      DO jm = 1, ndim
         SELECT CASE ( kn_umx )
            !
         CASE( 1 )                                                !==  1st order central TIM  ==! (Eq. 21)
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * (                          ( pt(ji,jj+1,jl,jm) + pt(ji,jj,jl,jm) ) &
                        &                                        - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl,jm) - pt(ji,jj,jl,jm) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 2 )                                                !==  2nd order central TIM  ==! (Eq. 23)
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                     pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * (                    ( pt(ji,jj+1,jl,jm) + pt(ji,jj,jl,jm) ) &
                        &                                                     - zcv   * ( pt(ji,jj+1,jl,jm) - pt(ji,jj,jl,jm) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 3 )                                                !==  3rd order central TIM  ==! (Eq. 24)
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                     zdy2 = e2v(ji,jj) * e2v(ji,jj)
                     !!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                     pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * ( (                  ( pt  (ji,jj+1,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcv   * ( pt  (ji,jj+1,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                       ( ztv2(ji,jj+1,jl,jm) + ztv2(ji,jj,jl,jm) )   &
                        &                                        - SIGN( 1._wp, zcv ) * ( ztv2(ji,jj+1,jl,jm) - ztv2(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 4 )                                                !==  4th order central TIM  ==! (Eq. 27)
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                     zdy2 = e2v(ji,jj) * e2v(ji,jj)
                     !!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                     pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * ( (                  ( pt  (ji,jj+1,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcv   * ( pt  (ji,jj+1,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6 * zdy2 * ( zcv*zcv - 1._wp ) * (                       ( ztv2(ji,jj+1,jl,jm) + ztv2(ji,jj,jl,jm) )   &
                        &                                            - 0.5_wp * zcv   * ( ztv2(ji,jj+1,jl,jm) - ztv2(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         CASE( 5 )                                                !==  5th order central TIM  ==! (Eq. 29)
            !
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     zcv  = pv(ji,jj) * r1_e1v(ji,jj) * pdt * r1_e2v(ji,jj)
                     zdy2 = e2v(ji,jj) * e2v(ji,jj)
                     !!rachid          zdy2 = e2v(ji,jj) * e2t(ji,jj)
                     zdy4 = zdy2 * zdy2
                     pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * ( (                  ( pt  (ji,jj+1,jl,jm) + pt  (ji,jj,jl,jm) )   &
                        &                                                     - zcv   * ( pt  (ji,jj+1,jl,jm) - pt  (ji,jj,jl,jm) ) ) &
                        & + r1_6   * zdy2 * ( zcv*zcv - 1._wp ) *    (                  ( ztv2(ji,jj+1,jl,jm) + ztv2(ji,jj,jl,jm) )   &
                        &                                            - 0.5_wp * zcv   * ( ztv2(ji,jj+1,jl,jm) - ztv2(ji,jj,jl,jm) ) ) &
                        & + r1_120 * zdy4 * ( zcv*zcv - 1._wp ) * ( zcv*zcv - 4._wp ) * ((ztv4(ji,jj+1,jl,jm) + ztv4(ji,jj,jl,jm) )   &
                        &                                        - SIGN( 1._wp, zcv ) * ( ztv4(ji,jj+1,jl,jm) - ztv4(ji,jj,jl,jm) ) ) )
                  END DO
               END DO
            END DO
            !
         END SELECT
         !
         ! if pt at v-point is negative then use the upstream value
         !    this should not be necessary if a proper sea-ice mask is set in Ultimate
         !    to degrade the order of the scheme when necessary (for ex. at the ice edge)
         IF( ll_neg ) THEN
            DO jl = 1, jpl
               DO jj=Njs0-2, Nje0+1
                  DO ji=Nis0-kloop, Nie0+kloop
                     IF( pt_v(ji,jj,jl,jm) < 0._wp .OR. ( jmsk_small(ji,jj,jl) == 0 .AND. pamsk(jm) == 0. ) ) THEN
                        pt_v(ji,jj,jl,jm) = 0.5_wp * vmask(ji,jj,1) * (                      ( pt(ji,jj+1,jl,jm) + pt(ji,jj,jl,jm) ) &
                           &                                    - SIGN( 1._wp, pv(ji,jj) ) * ( pt(ji,jj+1,jl,jm) - pt(ji,jj,jl,jm) ) )
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
         !                                                     !-- High order flux in j-direction  --!
         DO jl = 1, jpl
            DO jj=Njs0-2, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  pfv_ho(ji,jj,jl,jm) = pv(ji,jj) * pt_v(ji,jj,jl,jm)
               END DO
            END DO
         END DO
         !
      ENDDO
      !
      IF( kn_umx >= 3 )    DEALLOCATE( ztv1, ztv2 )
      IF( kn_umx == 5 )    DEALLOCATE( ztv3, ztv4 )
      !
   END SUBROUTINE ultimate_y


   SUBROUTINE nonosc_ice( pamsk, pdt, pu, pv, pt, pt_ups, pfu_ups, pfv_ups, pfu_ho, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc_ice  ***
      !!
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream
      !!       scheme and the before field by a non-oscillatory algorithm
      !!
      !! **  Method  :   ...
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   ) ::   pamsk            ! advection of concentration (1) or other tracers (0)
      REAL(wp)                   , INTENT(in   ) ::   pdt              ! tracer time-step
      REAL(wp), DIMENSION (:,:)  , INTENT(in   ) ::   pu               ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:)  , INTENT(in   ) ::   pv               ! ice j-velocity => v*e1
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pt, pt_ups       ! before field & upstream guess of after field
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pfv_ups, pfu_ups ! upstream flux
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_ho, pfu_ho   ! monotonic flux
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      REAL(wp) ::   zpos, zneg, zbig, zup, zdo, z1_dt              ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zcoef, zzt       !   -      -
      REAL(wp), DIMENSION(jpi,jpj)     ::   zbup, zbdo
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zbetup, zbetdo
      !!----------------------------------------------------------------------
      zbig = 1.e+20_wp   ! works ok with simple/double precison

      !
      ! antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-2, Nie0+1
               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) - pfu_ups(ji,jj,jl)
            END DO
         END DO
         DO jj=Njs0-2, Nje0+1
            DO ji=Nis0-1, Nie0+1
               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) - pfv_ups(ji,jj,jl)
            END DO
         END DO
      END DO

      ! extreme case where pfu_ho has to be zero
      ! ----------------------------------------
      !                                    pfu_ho
      !                           *         --->
      !                        .xOx.      .xOx.  *   .xOx.        .xOx.
      !                        .xOx.      .xOx.      .xOx.    *   .xOx.
      !                        .xOx.      .xOx.      .xOx.        .xOx.    *
      !            t_ups :       i-1     i       i+1       i+2
      IF( ll_prelim ) THEN

         DO jl = 1, jpl
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  IF( pfu_ho(ji,jj,jl) * ( pt_ups(ji+1,jj  ,jl) - pt_ups(ji,jj,jl) ) <= 0._wp .AND.  &
                     & pfv_ho(ji,jj,jl) * ( pt_ups(ji  ,jj+1,jl) - pt_ups(ji,jj,jl) ) <= 0._wp ) THEN
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( pt_ups(ji+2,jj  ,jl) - pt_ups(ji+1,jj  ,jl) ) <= 0._wp .AND.  &
                        & pfv_ho(ji,jj,jl) * ( pt_ups(ji  ,jj+2,jl) - pt_ups(ji  ,jj+1,jl) ) <= 0._wp ) THEN
                        pfu_ho(ji,jj,jl)=0._wp
                        pfv_ho(ji,jj,jl)=0._wp
                     ENDIF
                     !
                     IF(  pfu_ho(ji,jj,jl) * ( pt_ups(ji,jj,jl) - pt_ups(ji-1,jj  ,jl) ) <= 0._wp .AND.  &
                        & pfv_ho(ji,jj,jl) * ( pt_ups(ji,jj,jl) - pt_ups(ji  ,jj-1,jl) ) <= 0._wp ) THEN
                        pfu_ho(ji,jj,jl)=0._wp
                        pfv_ho(ji,jj,jl)=0._wp
                     ENDIF
                     !
                  ENDIF
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'icedyn_adv_umx', pfu_ho, 'U', -1.0_wp, pfv_ho, 'V', -1.0_wp )   ! lateral boundary cond.

      ENDIF

      ! Search local extrema
      ! --------------------
      ! max/min of pt & pt_ups with large negative/positive value (-/+zbig) outside ice cover
      z1_dt = 1._wp / pdt

      DO jl = 1, jpl

         DO jj=Njs0-2, Nje0+2
            DO ji=Nis0-2, Nie0+2
               IF    ( pt(ji,jj,jl) <= 0._wp .AND. pt_ups(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = -zbig
                  zbdo(ji,jj) =  zbig
               ELSEIF( pt(ji,jj,jl) <= 0._wp .AND. pt_ups(ji,jj,jl) > 0._wp ) THEN
                  zbup(ji,jj) = pt_ups(ji,jj,jl)
                  zbdo(ji,jj) = pt_ups(ji,jj,jl)
               ELSEIF( pt(ji,jj,jl) > 0._wp .AND. pt_ups(ji,jj,jl) <= 0._wp ) THEN
                  zbup(ji,jj) = pt(ji,jj,jl)
                  zbdo(ji,jj) = pt(ji,jj,jl)
               ELSE
                  zbup(ji,jj) = MAX( pt(ji,jj,jl) , pt_ups(ji,jj,jl) )
                  zbdo(ji,jj) = MIN( pt(ji,jj,jl) , pt_ups(ji,jj,jl) )
               ENDIF
            END DO
         END DO

         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !
               zup  = MAX( zbup(ji,jj), zbup(ji-1,jj), zbup(ji+1,jj), zbup(ji,jj-1), zbup(ji,jj+1) )  ! search max/min in neighbourhood
               zdo  = MIN( zbdo(ji,jj), zbdo(ji-1,jj), zbdo(ji+1,jj), zbdo(ji,jj-1), zbdo(ji,jj+1) )
               !
               zpos = MAX( 0._wp, pfu_ho(ji-1,jj  ,jl) ) - MIN( 0._wp, pfu_ho(ji  ,jj  ,jl) ) &  ! positive/negative part of the flux
                  & + MAX( 0._wp, pfv_ho(ji  ,jj-1,jl) ) - MIN( 0._wp, pfv_ho(ji  ,jj  ,jl) )
               zneg = MAX( 0._wp, pfu_ho(ji  ,jj  ,jl) ) - MIN( 0._wp, pfu_ho(ji-1,jj  ,jl) ) &
                  & + MAX( 0._wp, pfv_ho(ji  ,jj  ,jl) ) - MIN( 0._wp, pfv_ho(ji  ,jj-1,jl) )
               !
               zpos = zpos - (  pt(ji,jj,jl) * MIN( 0., pu(ji,jj) - pu(ji-1,jj) )   &
                  &           + pt(ji,jj,jl) * MIN( 0., pv(ji,jj) - pv(ji,jj-1) ) ) * ( 1. - pamsk )
               zneg = zneg + (  pt(ji,jj,jl) * MAX( 0., pu(ji,jj) - pu(ji-1,jj) )   &
                  &           + pt(ji,jj,jl) * MAX( 0., pv(ji,jj) - pv(ji,jj-1) ) ) * ( 1. - pamsk )
               !
               !                                  ! up & down beta terms
               ! clem: zbetup and zbetdo must be 0 for zpos>1.e-10 & zneg>1.e-10 (do not put 0 instead of 1.e-10 !!!)
               IF( zpos > epsi10 ) THEN
                  zbetup(ji,jj,jl) = MAX( 0._wp, zup - pt_ups(ji,jj,jl) ) / zpos * e1e2t(ji,jj) * z1_dt
               ELSE
                  zbetup(ji,jj,jl) = 0._wp ! zbig
               ENDIF
               !
               IF( zneg > epsi10 ) THEN
                  zbetdo(ji,jj,jl) = MAX( 0._wp, pt_ups(ji,jj,jl) - zdo ) / zneg * e1e2t(ji,jj) * z1_dt
               ELSE
                  zbetdo(ji,jj,jl) = 0._wp ! zbig
               ENDIF
               !
               ! if all the points are outside ice cover
               IF( zup == -zbig )   zbetup(ji,jj,jl) = 0._wp ! zbig
               IF( zdo ==  zbig )   zbetdo(ji,jj,jl) = 0._wp ! zbig
               !
            END DO
         END DO
      END DO

      ! monotonic flux in the y direction
      ! ---------------------------------
      DO jl = 1, jpl
         DO jj=Njs0-0, Nje0+0
            DO ji=Nis0-1, Nie0+0
               zau = MIN( 1._wp , zbetdo(ji,jj,jl) , zbetup(ji+1,jj,jl) )
               zbu = MIN( 1._wp , zbetup(ji,jj,jl) , zbetdo(ji+1,jj,jl) )
               zcu = 0.5_wp + SIGN( 0.5_wp , pfu_ho(ji,jj,jl) )
               !
               zcoef = ( zcu * zau + ( 1._wp - zcu ) * zbu )
               !
               pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) * zcoef + pfu_ups(ji,jj,jl)
               !
            END DO
         END DO

         DO jj=Njs0-1, Nje0+0
            DO ji=Nis0-0, Nie0+0
               zav = MIN( 1._wp , zbetdo(ji,jj,jl) , zbetup(ji,jj+1,jl) )
               zbv = MIN( 1._wp , zbetup(ji,jj,jl) , zbetdo(ji,jj+1,jl) )
               zcv = 0.5_wp + SIGN( 0.5_wp , pfv_ho(ji,jj,jl) )
               !
               zcoef = ( zcv * zav + ( 1._wp - zcv ) * zbv )
               !
               pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) * zcoef + pfv_ups(ji,jj,jl)
               !
            END DO
         END DO

      END DO
      !

   END SUBROUTINE nonosc_ice


   SUBROUTINE limiter_x( pdt, pu, pt, pfu_ups, pfu_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_x  ***
      !!
      !! **  Purpose :   compute flux limiter
      !!----------------------------------------------------------------------
      REAL(wp)                  , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION(:,:  ), INTENT(in   ) ::   pu           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pt           ! ice tracer
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pfu_ups      ! upstream flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pfu_ho       ! high order flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, uCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj) ::   zslpx       ! tracer slopes
      !!----------------------------------------------------------------------
      !

      DO jl = 1, jpl

         DO jj=Njs0-0, Nje0+0
            DO ji=Nis0-2, Nie0+1
               zslpx(ji,jj) = ( pt(ji+1,jj,jl) - pt(ji,jj,jl) ) * umask(ji,jj,1)
            END DO
         END DO

         DO jj=Njs0-0, Nje0+0
            DO ji=Nis0-1, Nie0+0
               uCFL = pdt * ABS( pu(ji,jj) ) * r1_e1e2t(ji,jj)

               Rjm = zslpx(ji-1,jj)
               Rj  = zslpx(ji  ,jj)
               Rjp = zslpx(ji+1,jj)

               IF( np_limiter == 3 ) THEN

                  IF( pu(ji,jj) > 0. ) THEN
                     Rr = Rjm
                  ELSE
                     Rr = Rjp
                  ENDIF

                  zh3 = pfu_ho(ji,jj,jl) - pfu_ups(ji,jj,jl)
                  IF( Rj > 0. ) THEN
                     zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(pu(ji,jj)),  &
                        &        MIN( 2. * Rr * 0.5 * ABS(pu(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(pu(ji,jj)) ) ) ) )
                  ELSE
                     zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(pu(ji,jj)),  &
                        &        MIN(-2. * Rr * 0.5 * ABS(pu(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(pu(ji,jj)) ) ) ) )
                  ENDIF
                  pfu_ho(ji,jj,jl) = pfu_ups(ji,jj,jl) + zlimiter

               ELSEIF( np_limiter == 2 ) THEN
                  IF( Rj /= 0. ) THEN
                     IF( pu(ji,jj) > 0. ) THEN
                        Cr = Rjm / Rj
                     ELSE
                        Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
                  ENDIF

                  ! -- superbee --
                  zpsi = MAX( 0., MAX( MIN(1.,2.*Cr), MIN(2.,Cr) ) )
                  ! -- van albada 2 --
                  !!zpsi = 2.*Cr / (Cr*Cr+1.)
                  ! -- sweby (with beta=1) --
                  !!zpsi = MAX( 0., MAX( MIN(1.,1.*Cr), MIN(1.,Cr) ) )
                  ! -- van Leer --
                  !!zpsi = ( Cr + ABS(Cr) ) / ( 1. + ABS(Cr) )
                  ! -- ospre --
                  !!zpsi = 1.5 * ( Cr*Cr + Cr ) / ( Cr*Cr + Cr + 1. )
                  ! -- koren --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( (1.+2*Cr)/3., 2. ) ) )
                  ! -- charm --
                  !IF( Cr > 0. ) THEN   ;   zpsi = Cr * (3.*Cr + 1.) / ( (Cr + 1.) * (Cr + 1.) )
                  !ELSE                 ;   zpsi = 0.
                  !ENDIF
                  ! -- van albada 1 --
                  !!zpsi = (Cr*Cr + Cr) / (Cr*Cr +1)
                  ! -- smart --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, 4. ) ) )
                  ! -- umist --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, MIN(0.75+0.25*Cr, 2. ) ) ) )

                  ! high order flux corrected by the limiter
                  pfu_ho(ji,jj,jl) = pfu_ho(ji,jj,jl) - ABS( pu(ji,jj) ) * ( (1.-zpsi) + uCFL*zpsi ) * Rj * 0.5

               ENDIF
            END DO
         END DO
      END DO
      !
   END SUBROUTINE limiter_x


   SUBROUTINE limiter_y( pdt, pv, pt, pfv_ups, pfv_ho )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE limiter_y  ***
      !!
      !! **  Purpose :   compute flux limiter
      !!----------------------------------------------------------------------
      REAL(wp)                   , INTENT(in   ) ::   pdt          ! tracer time-step
      REAL(wp), DIMENSION (:,:  ), INTENT(in   ) ::   pv           ! ice i-velocity => u*e2
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pt           ! ice tracer
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pfv_ups      ! upstream flux
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   pfv_ho       ! high order flux
      !
      REAL(wp) ::   Cr, Rjm, Rj, Rjp, vCFL, zpsi, zh3, zlimiter, Rr
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      REAL(wp), DIMENSION (jpi,jpj) ::   zslpy       ! tracer slopes
      !!----------------------------------------------------------------------
      !

      DO jl = 1, jpl

         DO jj=Njs0-2, Nje0+1
            DO ji=Nis0-0, Nie0+0
               zslpy(ji,jj) = ( pt(ji,jj+1,jl) - pt(ji,jj,jl) ) * vmask(ji,jj,1)
            END DO
         END DO

         DO jj=Njs0-1, Nje0+0
            DO ji=Nis0-0, Nie0+0
               vCFL = pdt * ABS( pv(ji,jj) ) * r1_e1e2t(ji,jj)

               Rjm = zslpy(ji,jj-1)
               Rj  = zslpy(ji,jj  )
               Rjp = zslpy(ji,jj+1)

               IF( np_limiter == 3 ) THEN

                  IF( pv(ji,jj) > 0. ) THEN
                     Rr = Rjm
                  ELSE
                     Rr = Rjp
                  ENDIF

                  zh3 = pfv_ho(ji,jj,jl) - pfv_ups(ji,jj,jl)
                  IF( Rj > 0. ) THEN
                     zlimiter =  MAX( 0., MIN( zh3, MAX(-Rr * 0.5 * ABS(pv(ji,jj)),  &
                        &        MIN( 2. * Rr * 0.5 * ABS(pv(ji,jj)),  zh3,  1.5 * Rj * 0.5 * ABS(pv(ji,jj)) ) ) ) )
                  ELSE
                     zlimiter = -MAX( 0., MIN(-zh3, MAX( Rr * 0.5 * ABS(pv(ji,jj)),  &
                        &        MIN(-2. * Rr * 0.5 * ABS(pv(ji,jj)), -zh3, -1.5 * Rj * 0.5 * ABS(pv(ji,jj)) ) ) ) )
                  ENDIF
                  pfv_ho(ji,jj,jl) = pfv_ups(ji,jj,jl) + zlimiter

               ELSEIF( np_limiter == 2 ) THEN

                  IF( Rj /= 0. ) THEN
                     IF( pv(ji,jj) > 0. ) THEN
                        Cr = Rjm / Rj
                     ELSE
                        Cr = Rjp / Rj
                     ENDIF
                  ELSE
                     Cr = 0.
                  ENDIF

                  ! -- superbee --
                  zpsi = MAX( 0., MAX( MIN(1.,2.*Cr), MIN(2.,Cr) ) )
                  ! -- van albada 2 --
                  !!zpsi = 2.*Cr / (Cr*Cr+1.)
                  ! -- sweby (with beta=1) --
                  !!zpsi = MAX( 0., MAX( MIN(1.,1.*Cr), MIN(1.,Cr) ) )
                  ! -- van Leer --
                  !!zpsi = ( Cr + ABS(Cr) ) / ( 1. + ABS(Cr) )
                  ! -- ospre --
                  !!zpsi = 1.5 * ( Cr*Cr + Cr ) / ( Cr*Cr + Cr + 1. )
                  ! -- koren --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( (1.+2*Cr)/3., 2. ) ) )
                  ! -- charm --
                  !IF( Cr > 0. ) THEN   ;   zpsi = Cr * (3.*Cr + 1.) / ( (Cr + 1.) * (Cr + 1.) )
                  !ELSE                 ;   zpsi = 0.
                  !ENDIF
                  ! -- van albada 1 --
                  !!zpsi = (Cr*Cr + Cr) / (Cr*Cr +1)
                  ! -- smart --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, 4. ) ) )
                  ! -- umist --
                  !!zpsi = MAX( 0., MIN( 2.*Cr, MIN( 0.25+0.75*Cr, MIN(0.75+0.25*Cr, 2. ) ) ) )

                  ! high order flux corrected by the limiter
                  pfv_ho(ji,jj,jl) = pfv_ho(ji,jj,jl) - ABS( pv(ji,jj) ) * ( (1.-zpsi) + vCFL*zpsi ) * Rj * 0.5

               ENDIF
            END DO
         END DO
      END DO
      !
   END SUBROUTINE limiter_y



   !!======================================================================
END MODULE icedyn_adv_umx
