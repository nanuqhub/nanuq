MODULE icethd
   !!======================================================================
   !!                  ***  MODULE icethd   ***
   !!   sea-ice : master routine for thermodynamics
   !!======================================================================
   !! History :  1.0  !  2000-01  (M.A. Morales Maqueda, H. Goosse, T. Fichefet) original code 1D
   !!            4.0  !  2018     (many people)       SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd       : thermodynamics of sea ice
   !!   ice_thd_init  : initialisation of sea-ice thermodynamics
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constants
   USE par_oce
   USE ice            ! sea-ice: variables
   USE sbc_oce , ONLY : utau, vtau, fatm_snow, ln_cpl_atm
   !USE oss_nnq , ONLY : sss_s, e3t_m, ssu_m, ssv_m, frq_m, ln_cpl_oce
   USE sbc_ice , ONLY : qsr_oce, qns_oce, qemp_oce, qsr_ice, qns_ice, dqns_ice, evap_ice, qprec_ice, qevap_ice, &
      &                 qml_ice, qcn_ice, qtr_ice_top
   USE icethd_zdf     ! sea-ice: vertical heat diffusion
   USE icethd_dh      ! sea-ice: ice-snow growth and melt
   USE icethd_da      ! sea-ice: lateral melting
   USE icethd_sal     ! sea-ice: salinity
   USE icethd_do      ! sea-ice: growth in open water
   !LOLOreadd:USE icethd_pnd     ! sea-ice: melt ponds
   USE iceitd  , ONLY : ice_itd_rem
   USE icecor         ! sea-ice: corrections
   USE icectl         ! sea-ice: control print
   !
   USE in_out_manager ! I/O manager
   USE iom            , ONLY : iom_miss_val, iom_put       ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd         ! called by limstp module
   PUBLIC   ice_thd_init    ! called by ice_init

   LOGICAL , ALLOCATABLE, DIMENSION(:,:,:) ::   llmsk
   CHARACTER(LEN=50)      ::   clname="cfl_icesalt.ascii"    ! ascii filename
   INTEGER , DIMENSION(3) ::   iloc
   REAL(wp)               ::   zcfl_drain_max, zcfl_flush_max
   INTEGER                ::   numcfl                        ! outfile unit


   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd  ***
      !!
      !! ** Purpose : This routine manages ice thermodynamics
      !!
      !! ** Action : - computation of oceanic sensible heat flux at the ice base
      !!                              energy budget in the leads
      !!                              net fluxes on top of ice and of ocean
      !!             - selection of grid cells with ice
      !!                - call ice_thd_zdf  for vertical heat diffusion
      !!                - call ice_thd_dh   for vertical ice growth and melt
      !!                - call ice_thd_pnd  for melt ponds
      !!                - call ice_thd_temp to  retrieve temperature from ice enthalpy
      !!                - call ice_thd_sal  for ice desalination
      !!                - call ice_thd_temp to  retrieve temperature from ice enthalpy
      !!                - call ice_thd_mono for extra lateral ice melt if active virtual thickness distribution
      !!                - call ice_thd_da   for lateral ice melt
      !!             - back to the geographic grid
      !!                - call ice_thd_rem  for remapping thickness distribution
      !!                - call ice_thd_do   for ice growth in leads
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      !
      LOGICAL , DIMENSION(jpi,jpj) ::   ll_ice_present
      INTEGER ::   ji, jj, jk, jl   ! dummy loop indices
      INTEGER :: k_np_ice           ! n. of points of domain with sea-ice for given category
      INTEGER :: k_np_mlt           ! n. of points of domain where complete ice-melting occured
      !!-------------------------------------------------------------------

      ! controls
      IF( ln_timing    )   CALL timing_start('icethd')                                                             ! timing
      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icediachk )   CALL ice_cons2D  (0, 'icethd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd: sea-ice thermodynamics'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! convergence tests
      !IF( ln_zdf_chkcvg ) THEN
      !   ALLOCATE( ztice_cvgerr(jpi,jpj,jpl) , ztice_cvgstp(jpi,jpj,jpl) )
      !   ztice_cvgerr = 0._wp ; ztice_cvgstp = 0._wp
      !ENDIF
      !
      !IF( ln_sal_chk )   CALL ice_thd_salchk( kt, 1 )
      !
      !-------------------------------------------------------------------------------------------!
      ! Thermodynamic computation (only on grid points covered by ice) => loop over ice categories
      !-------------------------------------------------------------------------------------------!
      !
      CALL ice_thd_frazil             !--- frazil ice: collection thickness (ht_i_new) & fraction of frazil (fraz_frac)
      !%acc update zelf( ht_i_new, fraz_frac )

      !$acc data create( ll_ice_present )

      DO jl = 1, jpl

         k_np_ice = 0
         !$acc parallel loop collapse(2) reduction(+:k_np_ice)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               IF ( a_i(ji,jj,jl) > epsi10 ) THEN
                  ll_ice_present(ji,jj) = .TRUE. ! select ice covered grid points
                  k_np_ice = k_np_ice + 1
               ELSE
                  ll_ice_present(ji,jj) = .FALSE.
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
         !%acc update zelf( ll_ice_present )

         IF ( k_np_ice > 0 ) THEN

            CALL ice_thd_unit_convert( jl, 1 )            ! --- & Change units of e_i, e_s from J/m2 to J/m3 --- !
            !=> UPDATES: e_i, e_s
            !%acc update zelf( e_i, e_s )

            !$acc parallel loop collapse(2)
            DO jj=Njs0-nn_hls, Nje0+nn_hls
               DO ji=Nis0-nn_hls, Nie0+nn_hls
                  ! --- some init --- !  (important to have them here)
                  dh_s_tot(ji,jj) = 0._wp
                  dh_i_bom(ji,jj) = 0._wp ; dh_i_itm(ji,jj) = 0._wp
                  dh_i_sub(ji,jj) = 0._wp ; dh_i_bog(ji,jj) = 0._wp
                  dh_snowice(ji,jj) = 0._wp ; dh_s_itm(ji,jj) = 0._wp
               END DO
            END DO
            !$acc end parallel loop

            CALL ice_thd_zdf( jl, ll_ice_present )                  ! --- Ice-Snow temperature --- !
            !=> UPDATES: cnd_ice, hfx_dif, hfx_err_dif, h_s, qcn_ice, qcn_ice_bot, qcn_ice_top, qns_ice, qtr_ice_bot, t1_ice, t_i, t_s, t_si, t_su, e_i, e_s

            ! ll_ice_present can be updated by thd_dh, due to melting
            IF( ln_icedH ) CALL ice_thd_dh(jl, ll_ice_present)                   ! --- Growing/Melting --- !
            !=> UPDATES:  h_i, h_s, ll_ice_present,a_i, e_s,  qml_ice, s_i, t_s, t_su
            !             dh_i_bog,dh_i_bom,dh_i_itm,dh_i_sub,dh_i_sum_2d,dh_s_itm,dh_snowice,dh_s_sum_2d
            !             sfx_bog,sfx_bom,sfx_bri,sfx_res,sfx_sni,sfx_sub,sfx_sum,hfx_bog,hfx_bom
            !             hfx_res,hfx_snw,hfx_spr,hfx_sub,hfx_sum,hfx_thd wfx_bog,wfx_bom,wfx_err_sub
            !             wfx_ice_sub,wfx_res,wfx_sni,wfx_snw_sni,wfx_snw_sub,wfx_snw_sum,wfx_spr,wfx_sum
            
            CALL ice_thd_temp(jl, ll_ice_present)  !LOLO maybe unnecessary!???    ! --- Temperature update --- !
            !

            !=> USED: a_i,dh_i_sum_2d,dh_s_sum_2d,h_i,ll_ice_present,sfx_bri,sfx_res,s_i,sss_s,sz_i,t_i,t_su
            CALL ice_thd_sal(jl, ll_ice_present)                  ! --- Ice salinity --- !
            !=> UPDATES: sfx_bri,sfx_res,s_i,sz_i
            
            CALL ice_thd_temp(jl, ll_ice_present)                 ! --- Temperature update --- !



            !lolo:IF( ln_icedH .AND. ln_virtual_itd ) &
            !lolo:&              CALL ice_thd_mono(jl, ll_ice_present)                 ! --- Extra lateral melting if virtual_itd --- !
            !
            ! ll_ice_present can be updated by thd_da
            IF( ln_icedA )    CALL ice_thd_da(jl, ll_ice_present)                   ! --- Lateral melting --- !
            
            !%acc update device( h_i, a_i, h_s, s_i, o_i, e_i, e_s )

            !
            CALL ice_thd_unit_convert( jl, 2 )            ! --- Change units of e_i, e_s from J/m3 to J/m2 --- !
            !%acc update zelf( e_i, e_s, szv_i, oa_i, sv_i, v_s, v_i )
            !
         ENDIF ! ll_ice_present
         !
         !
      END DO ! jl loop
      !
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'icethd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      !lolo:IF ( ln_pnd .AND. ln_icedH ) &
      !lolo:   &                    CALL ice_thd_pnd                      ! --- Melt ponds --- !
      
      IF( jpl > 1  ) THEN
         !=> USED:    a_i, a_i_b, e_i, e_s, h_i, h_i_b, oa_i, sv_i, szv_i, t_su, v_i, v_s
         !%acc update zelf( a_i, a_i_b, e_i, e_s, h_i, h_i_b, oa_i, sv_i, szv_i, t_su, v_i, v_s )
         CALL ice_itd_rem( kt )                ! --- Transport ice between thickness categories --- !
         !=> UPDATES: a_i, e_i, e_s, h_i, h_i_b, oa_i, sv_i, szv_i, t_su, v_i, v_s
         !%acc update zelf(   a_i, e_i, e_s, h_i, h_i_b, oa_i, sv_i, szv_i, t_su, v_i, v_s )
      ENDIF
      !
      IF( ln_icedO )   CALL ice_thd_do                       ! --- Frazil ice growth in leads --- !
      !

      !$acc end data

      CALL ice_cor( kt , 2 )                ! --- Corrections --- !

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            ! --- Ice natural aging incrementation
            !$acc loop seq
            DO jl=1, jpl
               oa_i(ji,jj,jl) = oa_i(ji,jj,jl) + a_i(ji,jj,jl) * rDt_ice
            END DO
         END DO
      END DO
      !$acc end parallel loop

# if ! defined _OPENACC
      !                                                             ! --- LBC for the halos --- !
      CALL lbc_lnk( 'icethd', a_i , 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp, oa_i, 'T', 1._wp, t_su, 'T', 1._wp )
      ! a_ip, 'T', 1._wp, v_ip, 'T', 1._wp, v_il, 'T', 1._wp
      CALL lbc_lnk( 'icethd', e_i , 'T', 1._wp, e_s , 'T', 1._wp, szv_i , 'T', 1._wp )
# endif
      !

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !#halo: ok because `a_i` has been lbc_lnked !
            at_i(ji,jj) = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
            END DO
         END DO
      END DO
      !$acc end parallel loop
      !%acc update zelf( at_i )

      ! --- Ice velocity corrections
      k_np_mlt = 0
      !$acc parallel loop collapse(2) reduction(+:k_np_mlt) present( u_ice, v_ice )
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF( at_i(ji,jj) == 0._wp ) THEN   ! if ice has melted
               k_np_mlt = k_np_mlt + 1
               IF( at_i(ji+1,jj) == 0._wp )   u_ice(ji  ,jj) = 0._wp   ! right side
               IF( at_i(ji-1,jj) == 0._wp )   u_ice(ji-1,jj) = 0._wp   ! left side
               IF( at_i(ji,jj+1) == 0._wp )   v_ice(ji,jj  ) = 0._wp   ! upper side
               IF( at_i(ji,jj-1) == 0._wp )   v_ice(ji,jj-1) = 0._wp   ! bottom side
            ENDIF
         END DO
      END DO
      !$acc end parallel loop
      !
# if ! defined _OPENACC      
      CALL mpp_sum( 'ice_thd', k_np_mlt)    ! very important !!!
      IF( k_np_mlt > 0 )  CALL lbc_lnk( 'icethd', u_ice,'U',-1._wp, v_ice,'V',-1._wp )
# endif      

      ! convergence tests
      IF( ln_zdf_chkcvg ) THEN
         CALL iom_put( 'tice_cvgerr', ztice_cvgerr ) ; DEALLOCATE( ztice_cvgerr )
         CALL iom_put( 'tice_cvgstp', ztice_cvgstp ) ; DEALLOCATE( ztice_cvgstp )
      ENDIF
      !
      ! sanity checks for salt drainage and flushing
      IF( ln_sal_chk )   CALL ice_thd_salchk( kt, 2 )

      ! controls
      IF( ln_icectl )   CALL ice_prt    (kt, iiceprt, jiceprt, 1, ' - ice thermodyn. - ') ! prints
      IF( sn_cfctl%l_prtctl )   &
         &               CALL ice_prt3D  ('icethd')                                        ! prints

      !%acc update zelf( e_i, e_s, szv_i, oa_i, sv_i, v_s, v_i, t_su )

      IF( ln_timing )   CALL timing_stop('icethd')                                        ! timing
      !
   END SUBROUTINE ice_thd


   SUBROUTINE ice_thd_temp(jl_cat, ll_ice_present)
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_temp ***
      !!
      !! ** Purpose :   Computes sea ice temperature (Kelvin) from enthalpy
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER,                     INTENT(in) :: jl_cat
      LOGICAL, DIMENSION(jpi,jpj), INTENT(in) :: ll_ice_present
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztmelts, zbbb, zccc  ! local scalar
      !!-------------------------------------------------------------------
      IF( ln_timing    )   CALL timing_start('ice_thd_temp')
      ! Recover ice temperature
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !
            IF ( ll_ice_present(ji,jj) ) THEN
               !$acc loop seq
               DO jk = 1, nlay_i
                  IF( h_i(ji,jj,jl_cat) > 0._wp ) THEN
                     ztmelts       = -rTmlt * sz_i(ji,jj,jk,jl_cat)
                     ! Conversion q(S,T) -> T (second order equation)
                     zbbb          = ( rcp - rcpi ) * ztmelts + e_i(ji,jj,jk,jl_cat) * r1_rhoi - rLfus
                     zccc          = SQRT( MAX( zbbb * zbbb - 4._wp * rcpi * rLfus * ztmelts, 0._wp ) )
                     t_i(ji,jj,jk,jl_cat) = rt0 - ( zbbb + zccc ) * 0.5_wp * r1_rcpi
                  ELSE
                     t_i(ji,jj,jk,jl_cat) = rt0
                  ENDIF
               END DO
            ENDIF ! ll_ice_present
            !
         END DO
      END DO
      !$acc end parallel loop
      !
      IF( ln_timing    )   CALL timing_stop('ice_thd_temp')
      !
   END SUBROUTINE ice_thd_temp


   SUBROUTINE ice_thd_mono(jl_cat, ll_ice_present)
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_mono ***
      !!
      !! ** Purpose :   Lateral melting in case virtual_itd
      !!                          ( dA = A/2h dh )
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: jl_cat
      LOGICAL, DIMENSION(jpi,jpj) :: ll_ice_present
      INTEGER  ::   ji,jj              ! dummy loop indices
      REAL(wp) ::   zhi_bef            ! ice thickness before thermo
      REAL(wp) ::   zdh_mel, zda_mel   ! net melting
      REAL(wp) ::   zvi, zvs           ! ice/snow volumes
      !!-----------------------------------------------------------------------
      !%acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !
            IF (ll_ice_present(ji,jj)) THEN
               zdh_mel = MIN( 0._wp, dh_i_itm(ji,jj) + dh_i_sum_2d(ji,jj,jl_cat) + dh_i_bom(ji,jj) + dh_snowice(ji,jj) + dh_i_sub(ji,jj) )
               IF( zdh_mel < 0._wp .AND. a_i(ji,jj,jl_cat) > 0._wp )  THEN
                  zvi          = a_i(ji,jj,jl_cat) * h_i(ji,jj,jl_cat)
                  zvs          = a_i(ji,jj,jl_cat) * h_s(ji,jj,jl_cat)
                  ! lateral melting = concentration change
                  zhi_bef     = h_i(ji,jj,jl_cat) - zdh_mel
                  zda_mel     = MAX( -a_i(ji,jj,jl_cat) , a_i(ji,jj,jl_cat) * zdh_mel / ( 2._wp * MAX( zhi_bef, epsi20 ) ) )
                  a_i(ji,jj,jl_cat)  = MAX( epsi20, a_i(ji,jj,jl_cat) + zda_mel )
                  ! adjust thickness
                  h_i(ji,jj,jl_cat) = zvi / a_i(ji,jj,jl_cat)
                  h_s(ji,jj,jl_cat) = zvs / a_i(ji,jj,jl_cat)
                  ! retrieve total concentration
                  at_i(ji,jj) = a_i(ji,jj,jl_cat)
               ENDIF
            ENDIF
            !
         END DO
      END DO
      !%acc end parallel loop
      !
   END SUBROUTINE ice_thd_mono

   SUBROUTINE ice_thd_salchk( kt, kn )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_salchk ***
      !!
      !! ** Purpose :   checking salt drainage and flushing
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, kn
      !
      !INTEGER ::   jk   ! dummy loop indices
      !!-----------------------------------------------------------------------

      ! sanity checks for salt flushing and drainage
      IF( kn == 1 ) THEN
         !
         ALLOCATE( llmsk(jpi,jpj,jpl) )
         ALLOCATE( zcfl_flush(jpi,jpj,jpl) , zcfl_drain(jpi,jpj,jpl), zsneg_flush(jpi,jpj,jpl) , zsneg_drain(jpi,jpj,jpl) )
         !
         zcfl_flush = 0._wp ; zcfl_drain = 0._wp
         zsneg_flush = 0._wp ; zsneg_drain = 0._wp
         !
      ELSEIF( kn == 2 ) THEN
         !
         CALL iom_put( 'sice_flush_neg', zsneg_flush ) ; DEALLOCATE( zsneg_flush )
         CALL iom_put( 'sice_drain_neg', zsneg_drain ) ; DEALLOCATE( zsneg_drain )
         !
         CALL iom_put( 'cfl_flush', zcfl_flush )
         CALL iom_put( 'cfl_drain', zcfl_drain )

         !                    ! calculate maximum values and locations
         llmsk(Nis0:Nie0,Njs0:Nje0,:) = h_i(Nis0:Nie0,Njs0:Nje0,:) > rn_himin        ! define only where h > 0.10m
         CALL mpp_maxloc( 'icethd', zcfl_drain, llmsk, zcfl_drain_max, iloc )
         CALL mpp_maxloc( 'icethd', zcfl_flush, llmsk, zcfl_flush_max, iloc )

         IF( lwp ) THEN       ! write out to file
            WRITE(numcfl,FMT='(2x,i6,3x,a10,4x,f8.4,1x,i4,1x,i4,1x,i4)') kt, 'Max Cdrain', zcfl_drain_max, iloc(1), iloc(2), iloc(3)
            WRITE(numcfl,FMT='(11x,     a10,4x,f8.4,1x,i4,1x,i4,1x,i4)')     'Max Cflush', zcfl_flush_max, iloc(1), iloc(2), iloc(3)
         ENDIF
         DEALLOCATE( zcfl_flush, zcfl_drain )
         DEALLOCATE( llmsk )

         IF( kt == nitend .AND. lwp )   CLOSE( numcfl )
      ENDIF

   END SUBROUTINE ice_thd_salchk

   SUBROUTINE ice_thd_unit_convert( kl, kn )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_unit_convert ***
      !!
      !! ** Purpose :   to handle the conversion of units from J/m2 to J/m3 and reverse
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kl   ! index of the ice category
      INTEGER, INTENT(in) ::   kn   ! 1= from J/m2 to J/m3   ;   2= from J/m3 to J/m2
      !!-----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      SELECT CASE( kn )
         !                 !-------------------------!
      CASE( 1 )            !==  from J/m2 to J/m3  ==!
         !                 !-------------------------!
         ! --- Change units of e_i, e_s from J/m2 to J/m3 --- !
         ! Here we make sure that we don't divide by very small, but physically
         ! meaningless, products of sea ice thicknesses/snow depths and sea ice
         ! concentration
         !
         !$acc parallel loop collapse(2) present( h_i, a_i, e_i, h_s, e_s )
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               !$acc loop seq
               DO jk = 1, nlay_i
                  IF( (h_i(ji,jj,kl) * a_i(ji,jj,kl)) > epsi20 ) THEN
                     e_i(ji,jj,jk,kl) = e_i(ji,jj,jk,kl) / (h_i(ji,jj,kl) * a_i(ji,jj,kl)) * nlay_i
                  ELSE
                     e_i(ji,jj,jk,kl) = 0._wp
                  ENDIF
               END DO
               !$acc loop seq
               DO jk = 1, nlay_s
                  IF( (h_s(ji,jj,kl) * a_i(ji,jj,kl)) > epsi20 ) THEN
                     e_s(ji,jj,jk,kl) = e_s(ji,jj,jk,kl) / (h_s(ji,jj,kl) * a_i(ji,jj,kl)) * nlay_s
                  ELSE
                     e_s(ji,jj,jk,kl) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
         !$acc end parallel loop
         !
         !                 !-------------------------!
      CASE( 2 )            !==  from J/m3 to J/m2  ==!
         !                 !-------------------------!
         !$acc parallel loop collapse(2) present( h_i, a_i, e_i, h_s, e_s, v_i, v_s, sv_i, oa_i, s_i, o_i, szv_i, sz_i )
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               ! Change thickness to volume (replaces routine ice_var_eqv2glo)
               v_i(ji,jj,kl)   = h_i(ji,jj,kl)   * a_i(ji,jj,kl)
               v_s(ji,jj,kl)   = h_s(ji,jj,kl)   * a_i(ji,jj,kl)
               sv_i(ji,jj,kl)  = s_i(ji,jj,kl)   * v_i(ji,jj,kl)
               oa_i(ji,jj,kl)  = o_i(ji,jj,kl)   * a_i(ji,jj,kl)

               ! --- Change units of e_i, e_s from J/m3 to J/m2 --- !
               !$acc loop seq
               DO jk = 1, nlay_i
                  e_i(ji,jj,jk,kl)   =  e_i(ji,jj,jk,kl) * h_i(ji,jj,kl) * a_i(ji,jj,kl) * r1_nlay_i
                  szv_i(ji,jj,jk,kl) = sz_i(ji,jj,jk,kl) * v_i(ji,jj,kl)                 * r1_nlay_i
               END DO
               !$acc loop seq
               DO jk = 1, nlay_s
                  e_s(ji,jj,jk,kl) = e_s(ji,jj,jk,kl) * h_s(ji,jj,kl) * a_i(ji,jj,kl) * r1_nlay_s
               END DO
               !
            END DO
         END DO
         !$acc end parallel loop

      END SELECT
      !
   END SUBROUTINE ice_thd_unit_convert


   SUBROUTINE ice_thd_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice thermodynamics
      !!
      !! ** Method  :   Read the namthd namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer output status for namelist read
      !!
      NAMELIST/namthd/ ln_icedH, ln_icedA, ln_icedO, ln_leadhfx
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd)
      READ_NML_CFG(numnam_ice,namthd)
      IF(lwm) WRITE( numoni, namthd )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_init: Ice Thermodynamics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd:'
         WRITE(numout,*) '      activate ice thick change from top/bot (T) or not (F)                ln_icedH   = ', ln_icedH
         WRITE(numout,*) '      activate lateral melting (T) or not (F)                              ln_icedA   = ', ln_icedA
         WRITE(numout,*) '      activate ice growth in open-water (T) or not (F)                     ln_icedO   = ', ln_icedO
         WRITE(numout,*) '      heat in the leads is used to melt sea-ice before warming the ocean   ln_leadhfx = ', ln_leadhfx
      ENDIF
      !
      CALL ice_thd_zdf_init   ! set ice heat diffusion parameters
      IF( ln_icedA )   CALL ice_thd_da_init    ! set ice lateral melting parameters
      IF( ln_icedO )   CALL ice_thd_do_init    ! set ice growth in open water parameters
      CALL ice_thd_sal_init   ! set ice salinity parameters
      !LOLOaddme:CALL ice_thd_pnd_init   ! set melt ponds parameters
      !
      IF( ln_sal_chk ) THEN
         ! create output ascii file
         CALL ctl_opn( numcfl, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
         WRITE(numcfl,*) 'Timestep  Direction   Max C     i    j    k'
         WRITE(numcfl,*) '*******************************************'
      ENDIF
   END SUBROUTINE ice_thd_init

   !!======================================================================
END MODULE icethd
