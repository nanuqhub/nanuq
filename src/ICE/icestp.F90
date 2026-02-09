MODULE icestp
   !!======================================================================
   !!                       ***  MODULE  icestp  ***
   !! sea ice : Master routine for all the sea ice model
   !!=====================================================================
   !!
   !!    The Nanuq sea ice model is a fork of the sea ice model SI3 (Sea Ice
   !!    modelling Integrated Initiative), aka Sea Ice cube:
   !!    which is originally based on LIM3, developed in Louvain-la-Neuve by:
   !!       * Martin Vancoppenolle (UCL-ASTR, Belgium)
   !!       * Sylvain Bouillon (UCL-ASTR, Belgium)
   !!       * Miguel Angel Morales Maqueda (NOC-L, UK)
   !!      thanks to valuable earlier work by
   !!       * Thierry Fichefet
   !!       * Hugues Goosse
   !!      thanks also to the following persons who contributed
   !!       * Gurvan Madec, Claude Talandier, Christian Ethe (LOCEAN, France)
   !!       * Xavier Fettweis (UCL-ASTR), Ralph Timmermann (AWI, Germany)
   !!       * Bill Lipscomb (LANL), Cecilia Bitz (UWa) and Elisabeth Hunke (LANL), USA.
   !!
   !! SI3, and hence Nanuq, has been made possible by a handful of persons who met as working group
   !!      (from France, Belgium, UK and Italy)
   !!    * Clement Rousset, Martin Vancoppenolle & Gurvan Madec (LOCEAN, France)
   !!    * Matthieu Chevalier & David Salas (Meteo France, France)
   !!    * Gilles Garric (Mercator Ocean, France)
   !!    * Thierry Fichefet & Francois Massonnet (UCL, Belgium)
   !!    * Ed Blockley & Jeff Ridley (Met Office, UK)
   !!    * Danny Feltham & David Schroeder (CPOM, UK)
   !!    * Yevgeny Aksenov (NOC, UK)
   !!    * Paul Holland (BAS, UK)
   !!    * Dorotea Iovino (CMCC, Italy)
   !!======================================================================
   !! History :  4.0  !  2018     (C. Rousset)      Original code SI3
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_stp       : sea-ice model time-stepping and update ocean SBC over ice-covered area
   !!   ice_init      : initialize sea-ice
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   !
   USE par_ice
   USE ice
   !
   USE phycst         ! Define parameters for the routines
   USE eosbn2  , ONLY : eos10_fzp_2d_gpu
   USE oss_nnq , ONLY : ssu_m, ssv_m, sss_m, sss_s, sst_m, sst_s, ln_prs_oce, ln_cpl_oce, mld_m
   USE ossprs  , ONLY : ln_slab_sst, oss_prs_slab
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice   fields
   !
   USE icesbc         ! sea-ice: Surface boundary conditions
   USE icebbc         ! sea-ice: Bottom boundary conditions
   USE remap_weno
   USE icedyn         ! sea-ice: dynamics
   USE icethd         ! sea-ice: thermodynamics
   USE iceupdate      ! sea-ice: sea surface boundary condition update
   USE icedia         ! sea-ice: budget diagnostics
   USE icewri         ! sea-ice: outputs
   USE icerst         ! sea-ice: restarts
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control
   USE iceistate      ! sea-ice: initial state
   USE iceitd         ! sea-ice: remapping thickness distribution
   USE icealb         ! sea-ice: albedo
   !
   USE remap_classic, ONLY : rmpT2F, do_Voce
   !
   USE bdy , ONLY : ln_bdy   ! flag for bdy
   USE bdyice         ! unstructured open boundary data for sea-ice
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_stp    ! called by sbcmod.F90
   PUBLIC   ice_init   ! called by sbcmod.F90

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icestp.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_stp( kt, ksbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_stp  ***
      !!
      !! ** Purpose :   sea-ice model time-stepping and update ocean surface
      !!              boundary condition over ice-covered area
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine
      !!              - call the ice advection/diffusion routine
      !!              - call the ice thermodynamics routine
      !!              - call the routine that computes mass and
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      INTEGER, INTENT(in) ::   ksbc     ! flux formulation (user defined, bulk, or Pure Coupled)
      !!---------------------------------------------------------------------
      LOGICAL :: l_normal, l_do_diags
      INTEGER :: ji, jj, jl   ! dummy loop index
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_stp')

      kt_ice = kt                              ! -- Ice model time step

      l_normal   = ( .NOT. ln_dynADV2D )       ! i.e. it's not a super-idealized simulation

      l_do_diags = ( ln_icediachk .OR. iom_use('hfxdhc') )


      IF( l_normal ) THEN

         !#LOLO: mv to OSS !?
         ! -- mean surface ocean current
         CALL do_Voce( ssu_m, ssv_m,  V_oce )
         !#LOLO.

         CALL eos10_fzp_2d_gpu( sss_m, t_bo )   ! -- freezing temperature based on salinity [C]



         IF( ln_prs_oce ) THEN
            !! => we use a prescribed surface state of the ocean read into netCDF files,
            !!    sometimes `sst_m` is not consistent with the equation of state and can
            !!    be colder than the freezing point (true at least in GLORYS4...)
            !$acc parallel loop collapse(2) present( sst_m, t_bo )
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  sst_m(ji,jj) = MAX( sst_m(ji,jj) , t_bo(ji,jj) ) ! prescribed/observed SST, cannot be colder than freezing-point temperature.
               END DO
            END DO
            !$acc end parallel loop

            !! Update `sst_s`,  `sss_s` & `t_bo` based on a simplistic slab-ocean model approach:
            IF( ln_slab_sst )  CALL oss_prs_slab( kt, rDt_ice, at_i, sst_m, sss_m, mld_m, qsr_tot, qns_b, emp_b, sst_s, sss_s, t_bo )
            !
         ENDIF !IF( ln_prs_oce .AND. ln_icethd )

         IF( ln_cpl_oce .OR. (ln_prs_oce .AND. (.NOT. ln_slab_sst) ) ) THEN
            !! ==> `sst_s` & `sss_s` default to `sst_m` & `sst_s`:
            !IF(lwp) PRINT *, ' * `sst_s & sss_s` forced to (prescribed) `sst_m & sss_m` ! kt =', kt
            !$acc parallel loop collapse(2) present( sst_m, sss_m, sst_s, sss_s )
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  sst_s(ji,jj) = sst_m(ji,jj)
                  sss_s(ji,jj) = sss_m(ji,jj)
               END DO
            END DO
            !$acc end parallel loop
         ENDIF


         !! * `t_bo` to Kelvin & `=rt0` over land:
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               t_bo(ji,jj) = ( t_bo(ji,jj) + rt0 ) * xmskt(ji,jj) + rt0 * ( 1._wp - xmskt(ji,jj) )
            END DO
         END DO
         !$acc end parallel loop

      ENDIF !IF( l_normal )


      IF( iom_use( 'sss_s') ) CALL iom_put (  'sss_s' ,  sss_s         *xmskt )
      IF( iom_use( 'sst_s') ) CALL iom_put (  'sst_s' ,  sst_s         *xmskt )
      IF( iom_use('dsss_s') ) CALL iom_put ( 'dsss_s' , (sss_s - sss_m)*xmskt )
      IF( iom_use('dsst_s') ) CALL iom_put ( 'dsst_s' , (sst_s - sst_m)*xmskt )



      CALL store_fields()             ! Store now ice values


      !------------------------------------------------!
      ! --- Dynamical coupling with the atmosphere --- !
      !------------------------------------------------!
      ! It provides the following fields used in sea ice model:
      !    utau_ice, vtau_ice   = surface ice stress at T-points [N/m2]
      !------------------------------------------------!
      !------------------------------------------------------!
      ! --- Thermodynamical coupling with the atmosphere --- !
      !------------------------------------------------------!
      ! It provides the following fields used in sea ice model:
      !    emp_oce , emp_ice    = E-P over ocean and sea ice                    [Kg/m2/s]
      !    fatm_snow            = solid precipitation                           [Kg/m2/s]
      !    evap_ice             = sublimation                                   [Kg/m2/s]
      !    qsr_tot , qns_tot    = solar & non solar heat flux (total)           [W/m2]
      !    qsr_ice , qns_ice    = solar & non solar heat flux over ice          [W/m2]
      !    dqns_ice             = non solar  heat sensistivity                  [W/m2]
      !    qemp_oce, qemp_ice,  = sensible heat (associated with evap & precip) [W/m2]
      !    qprec_ice, qevap_ice
      !------------------------------------------------------!
      IF( l_normal ) THEN
         CALL ice_sbc( kt, ksbc, utau_ice, vtau_ice )
      ENDIF

      !PRINT *, ' * LOLO1 max At=', MAXVAL(at_i(:,:)*xmskt(:,:)), MAX_VAR_GPU( at_i, xmskt )

      !-------------------------------------!
      ! --- ice dynamics and advection  --- !
      !-------------------------------------!
      IF(ln_icethd) THEN
         CALL diag_set0_gpu()              ! set diag of mass, heat and salt fluxes to 0
         !%acc update zelf( sfx, sfx_bri, sfx_sni, sfx_bog, sfx_bom, sfx_res, sfx_lam, sfx_opw, sfx_dyn, sfx_sum, sfx_sub )
         !%acc update zelf( wfx_snw, wfx_ice, wfx_sni, wfx_opw, wfx_bog, wfx_dyn, wfx_bom, wfx_sum, wfx_res, wfx_sub, wfx_spr, wfx_lam, wfx_snw_dyn, wfx_snw_sum, wfx_ice_sub, wfx_snw_sni, wfx_pnd )
         !%acc update zelf( hfx_thd, hfx_snw, hfx_bog, hfx_bom, hfx_res, hfx_spr, hfx_err_dif, wfx_err_sub, hfx_opw, hfx_dyn, hfx_sum, hfx_sub, hfx_dif, qsb_ice_bot, fhld )
         !%acc update zelf( dh_s_sum_2d, dh_i_sum_2d, qml_ice, qtr_ice_bot, qcn_ice, cnd_ice, qcn_ice_top, qcn_ice_bot, t_si )
      ENDIF

      CALL ice_rst_opn( kt )        ! Open Ice restart file (if necessary)

      !PRINT *, ' * LOLO2 max At=', MAXVAL(at_i(:,:)*xmskt(:,:)), MAX_VAR_GPU( at_i, xmskt )
      !CALL test4inf( ' a_i@icestp b ice_dyn  ', a_i )


      IF( ln_icedyn ) CALL ice_dyn( kt )       ! -- Ice dynamics


      !IF( l_do_diags ) CALL diag_trends(1)         ! record dyn trends  ! NOT GPU PORTED YET !!!!


      !                          !==  lateral boundary conditions  ==!
      IF( ln_icethd .AND. ln_bdy ) THEN
# if defined _OPENACC
         CALL ctl_stop('STOP', 'ice_stp() : finish GPU for `ln_bdy=T` !!!')
# endif
         CALL bdy_ice( kt )            ! -- bdy ice thermo
      ENDIF


      IF( l_normal ) THEN
         !CALL test4inf( ' a_i@icestp b glo2eqv  ', a_i )
         !                          !==  previous lead fraction and ice volume for flux calculations
         CALL ice_var_glo2eqv_gpu(2)   !2       ! h_i and h_s for ice albedo calculation

         !CALL test4inf( ' a_i@icestp p glo2eqv  ', a_i )
         CALL ice_var_agg_gpu(1)         ! at_i for coupling
         ! => UPDATES: af_i,at_i,ato_i,au_i,av_i,et_i,et_s,hm_i,hm_i_f,hm_s,kmsk_ice_f,kmsk_ice_t,kmsk_ice_u,kmsk_ice_v,om_i,sm_i,st_i,tm_i,tm_s,tm_si,tm_su,vt_i,vt_s
         !         --> at_ip,hm_ip,vt_ip,hm_il,vt_il

      ELSE
         !! => we are the idealized 2D advection test-case:
         CALL ice_var_agg_adv2d_gpu()

      ENDIF !IF( l_normal )


      CALL store_fields()             ! Store now ice values


      !----------------------------!
      ! --- ice thermodynamics --- !
      !----------------------------!
      IF( ln_icethd ) THEN

         CALL ice_thd( kt )            ! -- Ice thermodynamics
         !%acc update zelf( e_i, e_s, szv_i, oa_i, sv_i, v_s, v_i, t_su )

         IF( l_do_diags ) CALL diag_trends( 2 )         ! record thermo trends

         CALL ice_var_glo2eqv_gpu(2)  ! 2!    ! necessary calls (at least for coupling)

         CALL ice_var_agg_gpu( 2 )     ! necessary calls (at least for coupling)
         ! ==> UPDATES: af_i,at_i,ato_i,au_i,av_i,et_i,et_s,hm_i,hm_i_f,hm_s,kmsk_ice_f,kmsk_ice_t,kmsk_ice_u,kmsk_ice_v,om_i,sm_i,st_i,tm_i,tm_s,tm_si,tm_su,vt_i,vt_s )

         CALL ice_update_flx( kt )     ! -- Update ocean surface mass, heat and salt fluxes (calls ice_alb() !!!)


      ELSE
         ! --> thermo is turned off !

         CALL ice_var_glo2eqv_gpu(2)          ! necessary calls (at least for coupling)

         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               tm_i(ji,jj) = rt0               ! prevents the controls (`nanuq.stat`) to blow when thermo turned off...
            END DO
         END DO
         !$acc end parallel loop

         !%acc update zelf( tm_i )

      ENDIF !IF( ln_icethd )


# if defined _OPENACC
      ! UPDATE TO CPU !
      ! We do it here once for all (needed in controls, icewri...), updated in `adv` in `ice_dyn`
      ! & potentially corrected in `ice_thd`, hence the present position:
      !$acc update self( u_ice, v_ice, vt_i, tm_i )
# endif

      !IF( ln_icediahsb )      CALL ice_dia( kt )            ! -- Diagnostics outputs                    ! NOT GPU PORTED YET !!!!
      !
      !IF( ln_icediachk )      CALL ice_drift_wri( kt )      ! -- Diagnostics outputs for conservation   ! NOT GPU PORTED YET !!!
      !
      IF( ln_dynADV2D ) THEN
         CALL ice_wri_adv( kt )        ! -- Ice outputs (minimal set for advection-only tests)
      ELSE
         CALL ice_wri( kt )            ! -- Ice outputs
      ENDIF
      !
      IF( lrst_ice )           CALL ice_rst_write( kt )      ! -- Ice restart file
      !
      IF( ln_icectl )          CALL ice_ctl( kt )            ! -- Control checks
      !

      !-------------------------!
      ! --- Ocean time step --- !
      !-------------------------!

      IF( l_normal ) THEN
         CALL ice_update_tau( kt )         ! -- update surface ocean stresses   !LOLO_NANUQ
         !%acc update zelf( taum, utau, vtau )
      ENDIF

      !!gm   remark, the ocean-ice stress is not saved in ice diag call above .....  find a solution!!!
      !
      IF( ln_timing )   CALL timing_stop('ice_stp')
      !
   END SUBROUTINE ice_stp


   SUBROUTINE ice_init()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init  ***
      !!
      !! ** purpose :   Initialize sea-ice parameters
      !!----------------------------------------------------------------------
      INTEGER ::   ierr
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'Sea Ice Model: NANUQ'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_init: Arrays allocation & Initialization of all routines & init state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      !                                ! Load the reference and configuration namelist files and open namelist output file
      CALL load_nml( numnam_ice_ref, 'namelist_ice_ref',    numout, lwm )
      CALL load_nml( numnam_ice_cfg, 'namelist_ice_cfg',    numout, lwm )
      IF(lwm) CALL ctl_opn( numoni , 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )
      !
      CALL par_init                ! set some ice run parameters
      !
      !
      !                                ! Allocate the ice arrays (sbc_ice already allocated in sbc_init)
      ierr =        ice_alloc        ()      ! ice variables
      ierr = ierr + sbc_ice_alloc    ()      ! surface boundary conditions


      CALL mpp_sum( 'ice_init', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ice_init : unable to allocate ice arrays')
      !
      IF(ln_icethd) CALL diag_set0()                 ! set diag of mass, heat and salt fluxes to 0: needed for Agrif child grids
      !
      CALL ice_itd_init                ! ice thickness distribution initialization
      !
      IF( ln_icethd ) CALL ice_thd_init     ! set ice thermodynics parameters (clem: important to call it first for melt ponds)
      !
      !IF(lwp) WRITE(numout,*) 'LOLO: calling `ice_sbc_init` from `ice_init` of icestp.F90 !'
      CALL ice_sbc_init                ! set ice-ocean and ice-atm. coupling parameters
      !IF(lwp) WRITE(numout,*) 'LOLO: exiting `ice_sbc_init` from `ice_init` of icestp.F90 !'
      !
      ! There are no `ice_oss_init` needed for now!!!
      !WRITE(numout,*) 'LOLO: calling `ice_oss_init` from icestp.F90 !'
      !CALL ice_oss_init                ! set ice-ocean and ice-atm. coupling parameters
      !
      CALL ice_istate_init             ! Initial sea-ice state
      IF ( ln_rstart .OR. nn_iceini_file == 2 ) THEN
         CALL ice_rst_read()         ! start from a restart file
      ELSE
         !IF(lwp) WRITE(numout,*) 'LOLO: calling `ice_istate` from `ice_init` of icestp.F90 !'
         CALL ice_istate( nit000 )   ! start from rest or read a file
         !IF(lwp) WRITE(numout,*) 'LOLO: exiting `ice_istate` from `ice_init` of icestp.F90 !'
      ENDIF
      CALL ice_var_glo2eqv(1)
      CALL ice_var_agg(1)
      !
      CALL remap_weno_init
      !
      CALL ice_dyn_init                ! set ice dynamics parameters
      !
      CALL ice_update_init             ! ice surface boundary condition
      !
      CALL ice_alb_init                ! ice surface albedo
      !
      CALL ice_dia_init                ! initialization for diags
      !
      CALL ice_drift_init              ! initialization for diags of conservation
      !
      fr_i  (:,:)   = at_i(:,:)        ! initialisation of sea-ice fraction
      tn_ice(:,:,:) = t_su(:,:,:)      ! initialisation of surface temp for coupled simu
      !
      IF( ln_rstart )  THEN
         CALL iom_close( numrir )  ! close input ice restart file
         IF(lrxios) CALL iom_context_finalize(      cr_icerst_cxt         )
      ENDIF
      !
      IF( ln_dynADV2D ) tm_i(:,:) = 271._wp   ! prevents the controls to blow up in ADV2D test-case...
      !

# if defined _OPENACC
      !! Time to update initialized arrays into GPU's memory before moving on !
      PRINT *, ''
      PRINT *, ' * info GPU: ice_init() => updating all processed arrays following initial state into GPU !'
      !$acc update device ( t_bo, cnd_ice, tn_ice, t1_ice, e_i, t_i, szv_i, sz_i, e_s, t_s, a_i, v_i, v_s, sv_i, oa_i, h_i, h_s, s_i, o_i, t_su, u_ice, v_ice, uVice, vUice, snwice_mass, snwice_mass_b )
      !$acc update device ( vt_i, vt_s, at_i, af_i, fr_i, tm_i, hm_i, ato_i, et_i, et_s, tm_su, st_i, hm_i_f, au_i, av_i, kmsk_ice_t, kmsk_ice_f, kmsk_ice_u, kmsk_ice_v, SIGMAt )
      IF( ln_damage ) THEN
         !$acc update device ( dmdt, dmdf, SIGMAf )
      ENDIF
      PRINT *, ''
# endif

   END SUBROUTINE ice_init


   SUBROUTINE par_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE par_init ***
      !!
      !! ** Purpose :   Definition generic parameters for ice model
      !!
      !! ** Method  :   Read namelist and check the parameter
      !!                values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampar
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer
      !!
      NAMELIST/nampar/ jpl, nlay_i, nlay_s, ln_virtual_itd, ln_icedyn, ln_icethd, rn_amax,  &
         &             cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir, ln_damage
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,nampar)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampar in reference namelist' )
      READ_NML_CFG(numnam_ice,nampar)
      !902   IF( ios > 0 )   CALL ctl_nam ( ios , 'nampar in configuration namelist' )
      IF(lwm) WRITE( numoni, nampar )
      !
      IF( .NOT. ln_icethd ) THEN
         jpl    = 1
         nlay_i = 1
         nlay_s = 1
         IF(lwp) WRITE(numout,*) '   `jpl` forced to `1` as thermodynamics is turned off!'
         IF(lwp) WRITE(numout,*) '   `nlay_i` & `nlay_s` forced to `1` "      "'
      ENDIF
      IF ( jpl > 1 .AND. ln_virtual_itd ) THEN
         ln_virtual_itd = .FALSE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   `ln_virtual_itd` forced to false as jpl>1, no need with multiple categories to emulate them'
      ENDIF
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   par_init: ice parameters shared among all the routines'
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      Namelist nampar: '
         WRITE(numout,*) '         number of ice  categories                           jpl       = ', jpl
         WRITE(numout,*) '         number of ice  layers                               nlay_i    = ', nlay_i
         WRITE(numout,*) '         number of snow layers                               nlay_s    = ', nlay_s
         WRITE(numout,*) '         virtual ITD param for jpl=1 (T) or not (F)     ln_virtual_itd = ', ln_virtual_itd
         WRITE(numout,*) '         Ice dynamics       (T) or not (F)                   ln_icedyn = ', ln_icedyn
         !!
         WRITE(numout,*) '         Ice thermodynamics (T) or not (F)                   ln_icethd = ', ln_icethd
         WRITE(numout,*) '         maximum ice concentration                                     = ', rn_amax
         WRITE(numout,*) '         use of "ice damage" tracer (for brittle rheologies) ln_damage = ', ln_damage
      ENDIF
      !                                        !--- change max ice concentration for roundoff errors
      rn_amax = MIN( rn_amax, 1._wp - epsi10 )

      !                                        !--- check consistency
      !
      IF( ln_cpl_atm .AND. nn_cats_cpl /= 1 .AND. nn_cats_cpl /= jpl ) THEN
         CALL ctl_stop( 'STOP', 'par_init: in coupled mode, nn_cats_cpl should be either 1 or jpl' )
      ENDIF
      !
      rDt_ice   = REAL(nn_fsbc) * rn_Dt          !--- sea-ice timestep and its inverse
      r1_Dt_ice = 1._wp / rDt_ice
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      ice timestep rDt_ice = nn_fsbc*rn_Dt = ', rDt_ice
      !
      r1_nlay_i = 1._wp / REAL( nlay_i, wp )   !--- inverse of nlay_i and nlay_s
      r1_nlay_s = 1._wp / REAL( nlay_s, wp )
      !

      !$acc update device( rn_amax, jpl, rDt_ice, r1_Dt_ice, nlay_i, nlay_s, r1_nlay_i, r1_nlay_s )
   END SUBROUTINE par_init


   SUBROUTINE store_fields()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE store_fields()  ***
      !!
      !! ** purpose :  store ice variables at "before" time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl      ! dummy loop index
      REAL(wp) ::   zA
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('store_fields')
      !$acc parallel loop collapse(3) present(a_i_b, v_i_b, v_s_b, h_i_b, h_s_b)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            DO jl=1, jpl
               !
               zA = a_i(ji,jj,jl)
               a_i_b(ji,jj,jl)  =  zA                ! ice area
               v_i_b(ji,jj,jl)  =  v_i(ji,jj,jl)     ! ice volume
               v_s_b(ji,jj,jl)  =  v_s(ji,jj,jl)     ! snow volume
               !
               h_i_b(ji,jj,jl) = 0._wp
               h_s_b(ji,jj,jl) = 0._wp
               IF( zA >= epsi20 ) THEN
                  h_i_b(ji,jj,jl) = v_i_b(ji,jj,jl) / zA   ! ice thickness
                  h_s_b(ji,jj,jl) = v_s_b(ji,jj,jl) / zA   ! snw thickness
               ENDIF
               !
            END DO
         END DO
      END DO
      !$acc end parallel loop

      IF( ln_icethd ) THEN
         !$acc parallel loop collapse(3) present(sv_i_b, e_s_b, e_i_b)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               DO jl=1, jpl
                  !
                  sv_i_b(ji,jj,jl) = sv_i(ji,jj,jl)     ! salt content
                  !$acc loop seq
                  DO jk=1, nlay_s
                     e_s_b(ji,jj,jk,jl) = e_s(ji,jj,jk,jl)   ! snow thermal energy
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_i
                     e_i_b(ji,jj,jk,jl) = e_i(ji,jj,jk,jl)   ! ice thermal energy
                  END DO
                  !
               END DO
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
      !   !%acc update device( v_ip, v_il )
      !   !%acc parallel loop collapse(3) present(v_ip_b,v_il_b, v_ip,v_il)
      !   DO jj=Njs0-nn_hls, Nje0+nn_hls
      !      DO ji=Nis0-nn_hls, Nie0+nn_hls
      !         DO jl=1, jpl
      !            v_ip_b(ji,jj,jl)  = v_ip(ji,jj,jl)     ! pond volume
      !            v_il_b(ji,jj,jl)  = v_il(ji,jj,jl)     ! pond lid volume
      !         END DO
      !      END DO
      !   END DO
      !   !%acc end parallel loop
      !   !%acc update zelf (v_ip_b,v_il_b)
      !ENDIF

      ! Ice velocities & total concentration
      !$acc parallel loop collapse(2) present(at_i_b,u_ice,v_ice,u_ice_b,v_ice_b)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            !
            u_ice_b(ji,jj) = u_ice(ji,jj)
            v_ice_b(ji,jj) = v_ice(ji,jj)
            !
            zA = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zA = zA + a_i_b(ji,jj,jl)
            END DO
            at_i_b(ji,jj) = zA
            !
         END DO
      END DO
      !$acc end parallel loop

      IF( ln_timing )   CALL timing_stop('store_fields')
   END SUBROUTINE store_fields


   SUBROUTINE diag_set0()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_set0()  ***
      !!
      !! ** purpose :  set ice-ocean and ice-atm. fluxes to zeros at the beggining
      !!               of the time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl      ! dummy loop index
      !!----------------------------------------------------------------------

      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ! needed for (at least) diag_adv_mass -> to be removed
            sfx    (ji,jj) = 0._wp   ;
            sfx_bri(ji,jj) = 0._wp   ;   sfx_lam(ji,jj) = 0._wp
            sfx_sni(ji,jj) = 0._wp   ;   sfx_opw(ji,jj) = 0._wp
            sfx_bog(ji,jj) = 0._wp   ;   sfx_dyn(ji,jj) = 0._wp
            sfx_bom(ji,jj) = 0._wp   ;   sfx_sum(ji,jj) = 0._wp
            sfx_res(ji,jj) = 0._wp   ;   sfx_sub(ji,jj) = 0._wp
            !
            wfx_snw(ji,jj) = 0._wp   ;   wfx_ice(ji,jj) = 0._wp
            wfx_sni(ji,jj) = 0._wp   ;   wfx_opw(ji,jj) = 0._wp
            wfx_bog(ji,jj) = 0._wp   ;   wfx_dyn(ji,jj) = 0._wp
            wfx_bom(ji,jj) = 0._wp   ;   wfx_sum(ji,jj) = 0._wp
            wfx_res(ji,jj) = 0._wp   ;   wfx_sub(ji,jj) = 0._wp
            wfx_spr(ji,jj) = 0._wp   ;   wfx_lam(ji,jj) = 0._wp
            wfx_snw_dyn(ji,jj) = 0._wp ; wfx_snw_sum(ji,jj) = 0._wp
            wfx_snw_sub(ji,jj) = 0._wp ; wfx_ice_sub(ji,jj) = 0._wp
            wfx_snw_sni(ji,jj) = 0._wp
            wfx_pnd(ji,jj) = 0._wp

            hfx_thd(ji,jj) = 0._wp   ;
            hfx_snw(ji,jj) = 0._wp   ;   hfx_opw(ji,jj) = 0._wp
            hfx_bog(ji,jj) = 0._wp   ;   hfx_dyn(ji,jj) = 0._wp
            hfx_bom(ji,jj) = 0._wp   ;   hfx_sum(ji,jj) = 0._wp
            hfx_res(ji,jj) = 0._wp   ;   hfx_sub(ji,jj) = 0._wp
            hfx_spr(ji,jj) = 0._wp   ;   hfx_dif(ji,jj) = 0._wp
            hfx_err_dif(ji,jj) = 0._wp
            wfx_err_sub(ji,jj) = 0._wp
            !
            diag_heat(ji,jj) = 0._wp ;   diag_sice(ji,jj) = 0._wp
            diag_vice(ji,jj) = 0._wp ;   diag_vsnw(ji,jj) = 0._wp
            diag_aice(ji,jj) = 0._wp ;   diag_vpnd(ji,jj) = 0._wp

            tau_icebfr (ji,jj) = 0._wp   ! landfast ice param only (clem: important to keep the init here)
            qsb_ice_bot(ji,jj) = 0._wp   ! (needed if ln_icethd=F)

            fhld(ji,jj) = 0._wp   ! needed if ln_icethd=F

            ! for control checks (ln_icediachk)
            diag_trp_vi(ji,jj) = 0._wp   ;   diag_trp_vs(ji,jj) = 0._wp
            diag_trp_ei(ji,jj) = 0._wp   ;   diag_trp_es(ji,jj) = 0._wp
            diag_trp_sv(ji,jj) = 0._wp
            !
            diag_adv_mass(ji,jj) = 0._wp
            diag_adv_salt(ji,jj) = 0._wp
            diag_adv_heat(ji,jj) = 0._wp
         END DO
      END DO

      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ! SIMIP diagnostics
               t_si       (ji,jj,jl) = rt0     ! temp at the ice-snow interface
               qcn_ice_bot(ji,jj,jl) = 0._wp
               qcn_ice_top(ji,jj,jl) = 0._wp   ! conductive fluxes
               cnd_ice    (ji,jj,jl) = 0._wp   ! effective conductivity at the top of ice/snow (ln_cndflx=T)
               qcn_ice    (ji,jj,jl) = 0._wp   ! conductive flux (ln_cndflx=T & ln_cndemule=T)
               qtr_ice_bot(ji,jj,jl) = 0._wp   ! part of solar radiation transmitted through the ice needed at least for outputs
               qml_ice    (ji,jj,jl) = 0._wp   ! surface melt heat flux
               ! Melt pond surface melt diagnostics (mv - more efficient: grouped into one water volume flux)
               dh_i_sum_2d(ji,jj,jl) = 0._wp
               dh_s_sum_2d(ji,jj,jl) = 0._wp
            END DO
         END DO
      ENDDO

   END SUBROUTINE diag_set0



   SUBROUTINE diag_set0_gpu()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_set0_gpu()  ***
      !!
      !! ** purpose :  set ice-ocean and ice-atm. fluxes to zeros at the beggining
      !!               of the time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl      ! dummy loop index
      !!----------------------------------------------------------------------
      !$acc data present( sfx, sfx_bri, sfx_sni, sfx_bog, sfx_bom, sfx_res, sfx_lam, sfx_opw, sfx_dyn, sfx_sum, sfx_sub )
      !$acc data present( wfx_snw, wfx_ice, wfx_sni, wfx_opw, wfx_bog, wfx_dyn, wfx_bom, wfx_sum, wfx_res, wfx_sub, wfx_spr, wfx_lam, wfx_snw_dyn, wfx_snw_sum, wfx_ice_sub, wfx_snw_sni, wfx_pnd )
      !$acc data present( hfx_thd, hfx_snw, hfx_bog, hfx_bom, hfx_res, hfx_spr, hfx_err_dif, wfx_err_sub, hfx_opw, hfx_dyn, hfx_sum, hfx_sub, hfx_dif, qsb_ice_bot, fhld )
      !$acc data present( dh_s_sum_2d, dh_i_sum_2d, qml_ice, qtr_ice_bot, qcn_ice, cnd_ice, qcn_ice_top, qcn_ice_bot, t_si )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ! needed for (at least) diag_adv_mass -> to be removed
            sfx    (ji,jj) = 0._wp   ;
            sfx_bri(ji,jj) = 0._wp   ;   sfx_lam(ji,jj) = 0._wp
            sfx_sni(ji,jj) = 0._wp   ;   sfx_opw(ji,jj) = 0._wp
            sfx_bog(ji,jj) = 0._wp   ;   sfx_dyn(ji,jj) = 0._wp
            sfx_bom(ji,jj) = 0._wp   ;   sfx_sum(ji,jj) = 0._wp
            sfx_res(ji,jj) = 0._wp   ;   sfx_sub(ji,jj) = 0._wp
            !
            wfx_snw(ji,jj) = 0._wp   ;   wfx_ice(ji,jj) = 0._wp
            wfx_sni(ji,jj) = 0._wp   ;   wfx_opw(ji,jj) = 0._wp
            wfx_bog(ji,jj) = 0._wp   ;   wfx_dyn(ji,jj) = 0._wp
            wfx_bom(ji,jj) = 0._wp   ;   wfx_sum(ji,jj) = 0._wp
            wfx_res(ji,jj) = 0._wp   ;   wfx_sub(ji,jj) = 0._wp
            wfx_spr(ji,jj) = 0._wp   ;   wfx_lam(ji,jj) = 0._wp
            wfx_snw_dyn(ji,jj) = 0._wp ; wfx_snw_sum(ji,jj) = 0._wp
            wfx_snw_sub(ji,jj) = 0._wp ; wfx_ice_sub(ji,jj) = 0._wp
            wfx_snw_sni(ji,jj) = 0._wp
            wfx_pnd(ji,jj) = 0._wp

            hfx_thd(ji,jj) = 0._wp   ;
            hfx_snw(ji,jj) = 0._wp   ;   hfx_opw(ji,jj) = 0._wp
            hfx_bog(ji,jj) = 0._wp   ;   hfx_dyn(ji,jj) = 0._wp
            hfx_bom(ji,jj) = 0._wp   ;   hfx_sum(ji,jj) = 0._wp
            hfx_res(ji,jj) = 0._wp   ;   hfx_sub(ji,jj) = 0._wp
            hfx_spr(ji,jj) = 0._wp   ;   hfx_dif(ji,jj) = 0._wp
            hfx_err_dif(ji,jj) = 0._wp
            wfx_err_sub(ji,jj) = 0._wp
            !
            !diag_heat(ji,jj) = 0._wp ;   diag_sice(ji,jj) = 0._wp
            !diag_vice(ji,jj) = 0._wp ;   diag_vsnw(ji,jj) = 0._wp
            !diag_aice(ji,jj) = 0._wp ;   diag_vpnd(ji,jj) = 0._wp

            !tau_icebfr (ji,jj) = 0._wp   ! landfast ice param only (clem: important to keep the init here)
            qsb_ice_bot(ji,jj) = 0._wp   ! (needed if ln_icethd=F)

            fhld(ji,jj) = 0._wp   ! needed if ln_icethd=F

            ! for control checks (ln_icediachk)
            !diag_trp_vi(ji,jj) = 0._wp   ;   diag_trp_vs(ji,jj) = 0._wp
            !diag_trp_ei(ji,jj) = 0._wp   ;   diag_trp_es(ji,jj) = 0._wp
            !diag_trp_sv(ji,jj) = 0._wp
            !
            !diag_adv_mass(ji,jj) = 0._wp
            !diag_adv_salt(ji,jj) = 0._wp
            !diag_adv_heat(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(3)
      DO jl = 1, jpl
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ! SIMIP diagnostics
               t_si       (ji,jj,jl) = rt0     ! temp at the ice-snow interface
               qcn_ice_bot(ji,jj,jl) = 0._wp
               qcn_ice_top(ji,jj,jl) = 0._wp   ! conductive fluxes
               cnd_ice    (ji,jj,jl) = 0._wp   ! effective conductivity at the top of ice/snow (ln_cndflx=T)
               qcn_ice    (ji,jj,jl) = 0._wp   ! conductive flux (ln_cndflx=T & ln_cndemule=T)
               qtr_ice_bot(ji,jj,jl) = 0._wp   ! part of solar radiation transmitted through the ice needed at least for outputs
               qml_ice    (ji,jj,jl) = 0._wp   ! surface melt heat flux
               ! Melt pond surface melt diagnostics (mv - more efficient: grouped into one water volume flux)
               dh_i_sum_2d(ji,jj,jl) = 0._wp
               dh_s_sum_2d(ji,jj,jl) = 0._wp
            END DO
         END DO
      ENDDO
      !$acc end parallel loop

      !$acc end data
      !$acc end data
      !$acc end data
      !$acc end data

   END SUBROUTINE diag_set0_gpu






   SUBROUTINE diag_trends( kn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_trends  ***
      !!
      !! ** purpose : diagnostics of the trends. Used for conservation purposes
      !!              and outputs
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kn    ! 1 = after dyn ; 2 = after thermo
      !!----------------------------------------------------------------------
      !
      ! --- trends of heat, salt, mass (used for conservation controls)
      IF( ln_icediachk .OR. iom_use('hfxdhc') ) THEN
         !
         diag_heat(:,:) = diag_heat(:,:) &
            &             - SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice &
            &             - SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
         diag_sice(:,:) = diag_sice(:,:) &
            &             + SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
         diag_vice(:,:) = diag_vice(:,:) &
            &             + SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
         diag_vsnw(:,:) = diag_vsnw(:,:) &
            &             + SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhos
         diag_vpnd(:,:) = diag_vpnd(:,:) &
            &             + SUM(     v_ip + v_il          - v_ip_b - v_il_b                , dim=3 ) * r1_Dt_ice * rhow
         !
         IF( kn == 2 )    CALL iom_put ( 'hfxdhc' , diag_heat )   ! output of heat trend
         !
      ENDIF
      !
      ! --- trends of concentration (used for simip outputs)
      IF( iom_use('afxdyn') .OR. iom_use('afxthd') .OR. iom_use('afxtot') ) THEN
         !
         diag_aice(:,:) = diag_aice(:,:) + SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice
         !
         IF( kn == 1 )   CALL iom_put( 'afxdyn' , diag_aice )                                           ! dyn trend
         IF( kn == 2 )   CALL iom_put( 'afxthd' , SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice ) ! thermo trend
         IF( kn == 2 )   CALL iom_put( 'afxtot' , diag_aice )                                           ! total trend
         !
      ENDIF
      !
   END SUBROUTINE diag_trends


   FUNCTION MAX_VAR_GPU( pv, pmsk )
      !----------------------------------------------------
      REAL(wp) :: MAX_VAR_GPU
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pv
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pmsk
      !----------------------------------------------------
      INTEGER :: ji, jj
      REAL(wp) :: zmax
      !----------------------------------------------------
      !$acc data present( pv, pmsk )
      zmax = -1.E-20_wp
      !$acc loop seq
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         !$acc loop seq
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF( pmsk(ji,jj) > 0.5_wp ) THEN
               zmax = MAX( zmax, pv(ji,jj) )
            ENDIF
         END DO
      END DO
      !
      !$acc end data
      MAX_VAR_GPU = zmax
      !!
   END FUNCTION MAX_VAR_GPU


   !!======================================================================
END MODULE icestp
