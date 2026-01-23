MODULE icesbc
   !!======================================================================
   !!                       ***  MODULE  icesbc  ***
   !! Sea-Ice :   air-ice sbc fields
   !!=====================================================================
   !! History :  4.0  !  2017-08  (C. Rousset)       Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE par_ice
   USE ice            ! sea-ice: variables
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE oss_nnq , ONLY : sst_s, ssh_m, ssh_m, ssu_m, ssv_m, e3t_m, frq_m
   USE sbc_ice        ! Surface boundary condition: ice   fields
   USE sbcblk         ! Surface boundary condition: bulk
   USE icealb         ! sea-ice: albedo
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing
   USE fldread        !!GS: needed by agrif

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_sbc       ! called by icestp.F90
   PUBLIC ice_sbc_init  ! called by icestp.F90

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icesbc.F90 15388 2021-10-17 11:33:47Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_sbc( kt, ksbc, ptaux_ai, ptauy_ai )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc  ***
      !!  ==> Wrapper for ice_sbc_tau + ice_sbc_flx
      !!-------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kt                   ! ocean time step
      INTEGER                     , INTENT(in   ) ::   ksbc                 ! type of sbc flux
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptaux_ai, ptauy_ai   ! air-ice stress at T-point  [N/m2]
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) ::   zCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj) ::   zCEi    ! evap/sublim. transfer coefficient  [-]
      REAL(wp), DIMENSION(jpi,jpj) ::   ztheta_zu_i ! air temperature adjusted at heaight `zu` [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zq_zu_i      ! air spec. hum. adjusted at heaight `zu` [kg/kg]
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc')
      !$acc data present( ptaux_ai, ptauy_ai ) create( zCHi, zCEi, ztheta_zu_i, zq_zu_i )
      !
      !
      CALL ice_sbc_tau( kt, ksbc, ptaux_ai, ptauy_ai, zCHi, zCEi, ztheta_zu_i, zq_zu_i )
      !%acc update self( wndm_ice, ptaux_ai, ptauy_ai, zCHi, zCEi, ztheta_zu_i, zq_zu_i )
      !
      CALL ice_sbc_flx( kt, ksbc, zCHi, zCEi, ztheta_zu_i, zq_zu_i )
      !%acc update self( alb_ice,qsr_ice,qla_ice,dqla_ice,qns_ice,dqns_ice,evap_ice,devap_ice,emp_oce,emp_ice,qemp_oce,qemp_ice,qns_tot,qsr_tot,qprec_ice,qevap_ice,qtr_ice_top )
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_sbc')
      !
   END SUBROUTINE ice_sbc



   SUBROUTINE ice_sbc_tau( kt, ksbc, ptaux_ai, ptauy_ai, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc_tau  ***
      !!
      !! ** Purpose : provide surface boundary condition for sea ice (momentum)
      !!
      !! ** Action  : It provides the following fields:
      !!              ptaux_ai, ptauy_ai : surface ice stress (U- & V-points) [N/m2]
      !!-------------------------------------------------------------------
      INTEGER                     , INTENT(in ) ::   kt                   ! ocean time step
      INTEGER                     , INTENT(in ) ::   ksbc                 ! type of sbc flux
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   ptaux_ai, ptauy_ai   ! air-ice stress at T-points  [N/m2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCEi    ! evap/sublim. transfer coefficient  [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   ptheta_zu_i ! air temperature adjusted at heaight `zu` [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pq_zu_i      ! air spec. hum. adjusted at heaight `zu` [kg/kg]
      !!
      INTEGER  ::   ji, jj                 ! dummy loop index
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_ai, ztauy_ai
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc_tau')
      !$acc data present( ptaux_ai, ptauy_ai, pCHi, pCEi, ptheta_zu_i, pq_zu_i, fatm_u, fatm_v, fatm_theta, fatm_q, fatm_slp, u_ice, v_ice, tm_su )
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_sbc_tau: Surface boundary condition for sea ice (momentum)'
         WRITE(numout,*)'~~~~~~~~~~~~~~~'
      ENDIF
      !
      SELECT CASE( ksbc )
      CASE( jp_blk     )
         CALL blk_ice_1( fatm_u, fatm_v, fatm_theta, fatm_q,  &   ! #LB: known from "sbc_oce" module...
            &            fatm_slp, u_ice, v_ice, tm_su    ,   &   ! inputs
            &            pCHi, pCEi, ptheta_zu_i, pq_zu_i,    &   ! outputs
            &            putaui = ptaux_ai, pvtaui = ptauy_ai            )       ! outputs
         !        CASE( jp_abl     )    ptaux_ai & ptauy_ai are computed in ablmod
         !         CASE( jp_cpl_atm )   ;    CALL sbc_cpl_ice_tau( ptaux_ai , ptauy_ai )   ! Coupled      formulation LOLO: coupled with atmosphere, not ocean
      END SELECT
      
# if ! defined _OPENACC
      CALL lbc_lnk( 'ice_sbc_tau', ptaux_ai,'T',-1._wp, ptauy_ai,'T',-1._wp )
# endif
      
      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_sbc_tau')
      !
   END SUBROUTINE ice_sbc_tau


   SUBROUTINE ice_sbc_flx( kt, ksbc, pCHi, pCEi, ptheta_zu_i, pq_zu_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc_flx  ***
      !!
      !! ** Purpose : provide surface boundary condition for sea ice (flux)
      !!
      !! ** Action  : It provides the following fields used in sea ice model:
      !!                emp_oce , emp_ice                        = E-P over ocean and sea ice                    [Kg/m2/s]
      !!                fatm_snow                                = solid precipitation                           [Kg/m2/s]
      !!                evap_ice                                 = sublimation                                   [Kg/m2/s]
      !!                qsr_tot , qns_tot                        = solar & non solar heat flux (total)           [W/m2]
      !!                qsr_ice , qns_ice                        = solar & non solar heat flux over ice          [W/m2]
      !!                dqns_ice                                 = non solar  heat sensistivity                  [W/m2]
      !!                qemp_oce, qemp_ice, qprec_ice, qevap_ice = sensible heat (associated with evap & precip) [W/m2]
      !!            + these fields
      !!                qsb_ice_bot                              = sensible heat at the ice bottom               [W/m2]
      !!                fhld, qlead                              = heat budget in the leads                      [W/m2]
      !!            + some fields that are not used outside this module:
      !!                qla_ice                                  = latent heat flux over ice                     [W/m2]
      !!                dqla_ice                                 = latent heat sensistivity                      [W/m2]
      !!                fatm_prcp                                = total  precipitation                          [Kg/m2/s]
      !!                alb_ice                                  = albedo above sea ice
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time step
      INTEGER, INTENT(in) ::   ksbc   ! flux formulation (user defined, bulk or Pure Coupled)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pCEi    ! evap/sublim. transfer coefficient  [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptheta_zu_i ! air temperature adjusted at heaight `zu` [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pq_zu_i      ! air spec. hum. adjusted at heaight `zu` [kg/kg]
      !!--------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc_flx')
      !$acc data present( pCHi, pCEi, ptheta_zu_i, pq_zu_i, t_su, h_i, h_s, alb_ice, fatm_theta, fatm_slp, fatm_dqlw, fatm_prcp, fatm_snow )

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_sbc_flx: Surface boundary condition for sea ice (flux)'
         WRITE(numout,*)'~~~~~~~~~~~~~~~'
      ENDIF
      !                     !== ice albedo ==!
      IF( ln_pnd_alb ) THEN
         CALL ice_alb_pnd( t_su, h_i, h_s, a_ip_eff, h_ip, alb_ice )
      ELSE
         CALL ice_alb(     t_su, h_i, h_s,                 alb_ice )
      ENDIF
# if defined _TRDBG
      !$acc update self( alb_ice )
      CALL TRDBG( 'ice_sbc_flx', 'alb_ice', alb_ice )
# endif

      !
      !SELECT CASE( ksbc )   !== fluxes over sea ice ==!
      !CASE( jp_blk, jp_abl )      !--- bulk formulation & ABL formulation

      CALL blk_ice_2( t_su, h_s, h_i, alb_ice, fatm_theta,       &   ! #LB: known from "sbc_oce" module...
         &            fatm_slp, fatm_dqlw, fatm_prcp, fatm_snow, &
         &            pCHi, pCEi, ptheta_zu_i, pq_zu_i )

      !                        !    compute conduction flux and surface temperature (as in Jules surface module)
      IF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN
# if defined _OPENACC
         CALL ctl_stop( 'ice_sbc_flx: routine `blk_ice_qcn` not ported yet to GPU!!! (`ln_cndflx .AND. .NOT.ln_cndemulate`)!' )
# endif
         CALL blk_ice_qcn( ln_virtual_itd, t_su, t_bo, h_s, h_i )
      ENDIF
      !
      !CASE ( jp_cpl_atm )         !--- coupled formulation to atmosphere
      !                            CALL sbc_cpl_ice_flx( kt, picefr=at_i_b, palbi=alb_ice, psst=sst_s, pist=t_su, phs=h_s, phi=h_i )
      !END SELECT
      !                     !== some fluxes at the ice-ocean interface and in the leads
      !
      IF( ln_icethd )  THEN
         CALL ice_flx_other()
         !%acc update self( fhld, qsb_ice_bot, qlead )
      ENDIF
      !
      !$acc end data
      IF( ln_timing )  CALL timing_stop('ice_sbc_flx')
      !
   END SUBROUTINE ice_sbc_flx


   SUBROUTINE ice_flx_dist( ptn_ice, palb_ice, pqns_ice, pqsr_ice, pdqn_ice, pevap_ice, pdevap_ice, k_flxdist )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_flx_dist  ***
      !!
      !! ** Purpose :   update the ice surface boundary condition by averaging
      !!              and/or redistributing fluxes on ice categories
      !!
      !! ** Method  :   average then redistribute
      !!
      !! ** Action  :   depends on k_flxdist
      !!                = -1  Do nothing (needs N(cat) fluxes)
      !!                =  0  Average N(cat) fluxes then apply the average over the N(cat) ice
      !!                =  1  Average N(cat) fluxes then redistribute over the N(cat) ice
      !!                                                 using T-ice and albedo sensitivity
      !!                =  2  Redistribute a single flux over categories
      !!-------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   k_flxdist  ! redistributor
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptn_ice    ! ice surface temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   palb_ice   ! ice albedo
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pqns_ice   ! non solar flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pqsr_ice   ! net solar flux
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdqn_ice   ! non solar flux sensitivity
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pevap_ice  ! sublimation
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdevap_ice ! sublimation sensitivity
      !
      INTEGER  ::   jl      ! dummy loop index
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z1_at_i   ! inverse of concentration
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z_qsr_m   ! Mean solar heat flux over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z_qns_m   ! Mean non solar heat flux over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z_evap_m  ! Mean sublimation over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z_dqn_m   ! Mean d(qns)/dT over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z_devap_m ! Mean d(evap)/dT over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zalb_m    ! Mean albedo over all categories
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   ztem_m    ! Mean temperature over all categories
      !!----------------------------------------------------------------------
      !
      WHERE ( at_i (:,:) > 0._wp )
         z1_at_i(:,:) = 1._wp / at_i (:,:)
      ELSEWHERE
         z1_at_i(:,:) = 0._wp
      END WHERE

      SELECT CASE( k_flxdist )       !==  averaged on all ice categories  ==!
         !
      CASE( 0 , 1 )
         !
         ALLOCATE( z_qns_m(jpi,jpj), z_qsr_m(jpi,jpj), z_dqn_m(jpi,jpj), z_evap_m(jpi,jpj), z_devap_m(jpi,jpj) )
         !
         z_qns_m  (:,:) = SUM( a_i(:,:,:) * pqns_ice  (:,:,:) , dim=3 ) * z1_at_i(:,:)
         z_qsr_m  (:,:) = SUM( a_i(:,:,:) * pqsr_ice  (:,:,:) , dim=3 ) * z1_at_i(:,:)
         z_dqn_m  (:,:) = SUM( a_i(:,:,:) * pdqn_ice  (:,:,:) , dim=3 ) * z1_at_i(:,:)
         z_evap_m (:,:) = SUM( a_i(:,:,:) * pevap_ice (:,:,:) , dim=3 ) * z1_at_i(:,:)
         z_devap_m(:,:) = SUM( a_i(:,:,:) * pdevap_ice(:,:,:) , dim=3 ) * z1_at_i(:,:)
         DO jl = 1, jpl
            pqns_ice  (:,:,jl) = z_qns_m (:,:)
            pqsr_ice  (:,:,jl) = z_qsr_m (:,:)
            pdqn_ice  (:,:,jl) = z_dqn_m  (:,:)
            pevap_ice (:,:,jl) = z_evap_m(:,:)
            pdevap_ice(:,:,jl) = z_devap_m(:,:)
         END DO
         !
         DEALLOCATE( z_qns_m, z_qsr_m, z_dqn_m, z_evap_m, z_devap_m )
         !
      END SELECT
      !
      SELECT CASE( k_flxdist )       !==  redistribution on all ice categories  ==!
         !
      CASE( 1 , 2 )
         !
         ALLOCATE( zalb_m(jpi,jpj), ztem_m(jpi,jpj) )
         !
         zalb_m(:,:) = SUM( a_i(:,:,:) * palb_ice(:,:,:) , dim=3 ) * z1_at_i(:,:)
         ztem_m(:,:) = SUM( a_i(:,:,:) * ptn_ice (:,:,:) , dim=3 ) * z1_at_i(:,:)
         DO jl = 1, jpl
            pqns_ice (:,:,jl) = pqns_ice (:,:,jl) + pdqn_ice  (:,:,jl) * ( ptn_ice(:,:,jl) - ztem_m(:,:) )
            pevap_ice(:,:,jl) = pevap_ice(:,:,jl) + pdevap_ice(:,:,jl) * ( ptn_ice(:,:,jl) - ztem_m(:,:) )
            pqsr_ice (:,:,jl) = pqsr_ice (:,:,jl) * ( 1._wp - palb_ice(:,:,jl) ) / ( 1._wp - zalb_m(:,:) )
         END DO
         !
         DEALLOCATE( zalb_m, ztem_m )
         !
      END SELECT
      !
   END SUBROUTINE ice_flx_dist


   SUBROUTINE ice_flx_other
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_flx_other ***
      !!
      !! ** Purpose :   prepare necessary fields for thermo calculations
      !!
      !! ** Inputs  :   u_ice, v_ice, ssu_m, ssv_m, utau, vtau
      !!                frq_m, qsr_oce, qns_oce, qemp_oce, e3t_m, sst_s
      !! ** Outputs :   qsb_ice_bot, fhld, qlead
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj             ! dummy loop indices
      REAL(wp) ::   zswitch
      REAL(wp) ::   zfric_u, zqld, zqfr, zqfr_neg, zqfr_pos, zu_io, zv_io, zu_iom1, zv_iom1
      REAL(wp), PARAMETER ::   zfric_umin = 0._wp       ! lower bound for the friction velocity (cice value=5.e-04)
      REAL(wp), PARAMETER ::   zch        = 0.0057_wp   ! heat transfer coefficient
      REAL(wp), DIMENSION(jpi,jpj) ::  zfric, zvel      ! ice-ocean velocity (m/s) and frictional velocity (m2/s2)
      !!-----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_flx_other')
      !$acc data present( ssu_m, ssv_m, u_ice, v_ice, qsb_ice_bot, fhld, qlead ) create( zfric, zvel )
      !
      ! computation of friction velocity at T points
      IF( ln_icedyn ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zu_io   = u_ice(ji  ,jj  ) - ssu_m(ji  ,jj  )
               zu_iom1 = u_ice(ji-1,jj  ) - ssu_m(ji-1,jj  )
               zv_io   = v_ice(ji  ,jj  ) - ssv_m(ji  ,jj  )
               zv_iom1 = v_ice(ji  ,jj-1) - ssv_m(ji  ,jj-1)

               zfric(ji,jj) = rn_Cd_io * ( 0.5_wp * ( zu_io*zu_io + zu_iom1*zu_iom1 + zv_io*zv_io + zv_iom1*zv_iom1 ) ) * xmskt(ji,jj)
               zvel (ji,jj) = 0.5_wp * SQRT( ( u_ice(ji-1,jj  ) + u_ice(ji,jj) ) * ( u_ice(ji-1,jj  ) + u_ice(ji,jj) ) + &
                  &                          ( v_ice(ji  ,jj-1) + v_ice(ji,jj) ) * ( v_ice(ji  ,jj-1) + v_ice(ji,jj) ) )
            END DO
         END DO
         !$acc end parallel loop
         !
      ELSE      !  if no ice dynamics => transfer directly the atmospheric stress to the ocean
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zfric(ji,jj) = r1_rho0 * SQRT( 0.5_wp *  &
                  &                         (  utau(ji,jj) * utau(ji,jj) + utau(ji-1,jj) * utau(ji-1,jj)   &
                  &                          + vtau(ji,jj) * vtau(ji,jj) + vtau(ji,jj-1) * vtau(ji,jj-1) ) ) * xmskt(ji,jj)
               zvel(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF


# if ! defined _OPENACC
      CALL lbc_lnk( 'icesbc', zfric, 'T',  1.0_wp, zvel, 'T', 1.0_wp )
# endif

      !--------------------------------------------------------------------!
      ! Partial computation of forcing for the thermodynamic sea ice model
      !--------------------------------------------------------------------!
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0! needed for qlead
            !
            zswitch  = xmskt(ji,jj) * MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - epsi10 ) ) ! 0 if no ice
            !
            ! --- Energy received in the lead from atm-oce exchanges, zqld is defined everywhere (J.m-2) --- !
            zqld =  xmskt(ji,jj) * rDt_ice *  &
               &    ( ( 1._wp - at_i_b(ji,jj) ) * qsr_oce(ji,jj) * frq_m(ji,jj) +  &
               &      ( 1._wp - at_i_b(ji,jj) ) * qns_oce(ji,jj) + qemp_oce(ji,jj) )

            ! --- Energy needed to bring ocean surface layer until its freezing, zqfr is defined everywhere (J.m-2) --- !
            !     (mostly<0 but >0 if supercooling)
            zqfr     = rho0 * rcp * e3t_m(ji,jj) * ( t_bo(ji,jj) - ( sst_s(ji,jj) + rt0 ) ) * xmskt(ji,jj)  ! both < 0 (t_bo < sst) and > 0 (t_bo > sst)
            zqfr_neg = MIN( zqfr , 0._wp )                                                                    ! only < 0
            zqfr_pos = MAX( zqfr , 0._wp )                                                                    ! only > 0

            ! --- Sensible ocean-to-ice heat flux (W/m2) --- !
            !     (mostly>0 but <0 if supercooling)
            zfric_u            = MAX( SQRT( zfric(ji,jj) ), zfric_umin )
            qsb_ice_bot(ji,jj) = zswitch * rho0 * rcp * zch * zfric_u * ( ( sst_s(ji,jj) + rt0 ) - t_bo(ji,jj) )

            ! upper bound for qsb_ice_bot: the heat retrieved from the ocean must be smaller than the heat necessary to reach
            !                              the freezing point, so that we do not have SST < T_freeze
            !                              This implies: qsb_ice_bot(ji,jj) * at_i(ji,jj) * rtdice <= - zqfr_neg
            !                              The following formulation is ok for both normal conditions and supercooling
            qsb_ice_bot(ji,jj) = zswitch * MIN( qsb_ice_bot(ji,jj), - zqfr_neg * r1_Dt_ice / MAX( at_i(ji,jj), epsi10 ) )

            ! If conditions are always supercooled (such as at the mouth of ice-shelves), then ice grows continuously
            ! ==> stop ice formation by artificially setting up the turbulent fluxes to 0 when volume > 20m (arbitrary)
            IF( ( t_bo(ji,jj) - ( sst_s(ji,jj) + rt0 ) ) > 0._wp .AND. vt_i(ji,jj) >= 20._wp ) THEN
               zqfr               = 0._wp
               zqfr_pos           = 0._wp
               qsb_ice_bot(ji,jj) = 0._wp
            ENDIF

            ! --- Energy Budget of the leads (qlead, J.m-2) --- !
            !     qlead is the energy received from the atm. in the leads.
            !     If warming (zqld >= 0), then the energy in the leads is used to melt ice (bottom melting) => fhld  (W/m2)
            !     If cooling (zqld <  0), then the energy in the leads is used to grow ice in open water    => qlead (J.m-2)
            IF( ( zqld - zqfr ) < 0._wp .OR. at_i(ji,jj) < epsi10 ) THEN
               fhld (ji,jj) = 0._wp
               ! upper bound for qlead: qlead should be equal to zqld
               !                        but before using this heat for ice formation, we suppose that the ocean cools down till the freezing point.
               !                        The energy for this cooling down is zqfr and freezing point is reached if zqfr = zqld
               !                        so the max heat that can be pulled out of the ocean is zqld - zqfr
               !                        The following formulation is ok for both normal conditions and supercooling
               qlead(ji,jj) = MIN( 0._wp , zqld - zqfr )
            ELSE
               ! upper bound for fhld: fhld should be equal to zqld
               !                        but we have to make sure that this heat will not make the sst drop below the freezing point
               !                        so the max heat that can be pulled out of the ocean is zqld - zqfr_pos
               !                        The following formulation is ok for both normal conditions and supercooling
               fhld (ji,jj) = zswitch * MAX( 0._wp, ( zqld - zqfr_pos ) * r1_Dt_ice / MAX( at_i(ji,jj), epsi10 ) )  ! divided by at_i since this is (re)multiplied by a_i in icethd_dh.F90
               qlead(ji,jj) = 0._wp
            ENDIF
            !
            ! If ice is landfast and ice concentration reaches its max
            ! => stop ice formation in open water
            IF(  zvel(ji,jj) <= 5.e-04_wp .AND. at_i(ji,jj) >= rn_amax-epsi06 )   qlead(ji,jj) = 0._wp
            !
            ! If the grid cell is almost fully covered by ice (no leads)
            ! => stop ice formation in open water
            IF( at_i(ji,jj) >= (1._wp - epsi10) )   qlead(ji,jj) = 0._wp
            !
            ! If ln_leadhfx is false
            ! => do not use energy of the leads to melt sea-ice
            IF( .NOT.ln_leadhfx )   fhld(ji,jj) = 0._wp
            !
         END DO
      END DO
      !$acc end parallel loop

      ! In case we bypass open-water ice formation
      IF( .NOT. ln_icedO ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               qlead(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      ! In case we bypass growing/melting from top and bottom
      IF( .NOT. ln_icedH ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               qsb_ice_bot(ji,jj) = 0._wp
               fhld       (ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_flx_other')
      !
   END SUBROUTINE ice_flx_other


   SUBROUTINE ice_sbc_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc_init  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice dynamics
      !!
      !! ** Method  :   Read the namsbc namelist and check the ice-dynamic
      !!              parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namsbc
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer
      !!
      NAMELIST/namsbc/ rn_snwblow, ln_cndflx, ln_cndemulate, nn_qtrice
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namsbc)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc in reference namelist' )
      READ_NML_CFG(numnam_ice,namsbc)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc in configuration namelist' )
      IF(lwm) WRITE( numoni, namsbc )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_sbc_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsbc:'
         WRITE(numout,*) '      coefficient for ice-lead partition of snowfall            rn_snwblow    = ', rn_snwblow
         WRITE(numout,*) '      Use conduction flux as surface condition                  ln_cndflx     = ', ln_cndflx
         WRITE(numout,*) '         emulate conduction flux                                ln_cndemulate = ', ln_cndemulate
         WRITE(numout,*) '      solar flux transmitted thru the surface scattering layer  nn_qtrice     = ', nn_qtrice
         WRITE(numout,*) '         = 0  Grenfell and Maykut 1977'
         WRITE(numout,*) '         = 1  Lebrun 2019'
         WRITE(numout,*) '   SI3: use per-category fluxes'
      ENDIF
      !
      !IF( (nn_snwfra<1).OR.(nn_snwfra>2) ) CALL ctl_stop( 'ice_sbc_init: `nn_snwfra` can omly be `1` or `2`!' )
      !
      !$acc update device( rn_snwblow, nn_qtrice )
   END SUBROUTINE ice_sbc_init

   !!======================================================================
END MODULE icesbc
