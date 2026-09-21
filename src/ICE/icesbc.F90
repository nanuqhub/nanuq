MODULE icesbc
   !!======================================================================
   !!                       ***  MODULE  icesbc  ***
   !! Sea-Ice :   air-ice sbc fields
   !!=====================================================================
   !! History :  4.0  !  2017-08  (C. Rousset)       Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE dom_oce , ONLY : xmskt, xtmp1, xtmp2
   USE phycst
   USE par_ice
   USE ice            ! sea-ice: variables
   USE sbc_oce , ONLY : sn_loc_vct_tau, jp_blk, jp_abl, utau, vtau, qemp_oce, qsr_oce, qns_oce, sfx, t_air_zu, q_air_zu
   USE oss_nnq , ONLY : sst_s, ssh_m, ssh_m, ssu_m, ssv_m, e3t_m, frq_m
   USE sbc_ice        ! Surface boundary condition: ice   fields
   USE sbcblk         ! Surface boundary condition: bulk
   !
   USE abl     , ONLY : nt_n, u_abl, v_abl, tq_abl
   USE par_abl , ONLY : jp_ta, jp_qa
   USE icealb         ! sea-ice: albedo
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp , ONLY : ctl_nam
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu
#else
   USE lbclnk
#endif

   USE iom
   USE timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_sbc       ! called by icestp.F90
   PUBLIC ice_sbc_init  ! called by icestp.F90
   PUBLIC ice_sbc_wri   ! called bu icestp.F90

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
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
      !$acc data present( ptaux_ai, ptauy_ai, CH_ice, CE_ice )
      !
      CALL ice_sbc_tau( kt, ksbc, ptaux_ai, ptauy_ai, CH_ice, CE_ice )
      !
      CALL ice_sbc_flx( kt, ksbc,                     CH_ice, CE_ice )
      !
      !$acc end data
   END SUBROUTINE ice_sbc



   SUBROUTINE ice_sbc_tau( kt, ksbc, ptaux_ai, ptauy_ai, pCHi, pCEi )
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
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCHi                 ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   pCEi                 ! evap/sublim. transfer coefficient  [-]
      !!
      INTEGER  ::   ji, jj                 ! dummy loop index
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc_tau')
      !%acc data present( ptaux_ai, ptauy_ai, pCHi, pCEi, t_air_zu, q_air_zu, sf(jp_wndi)%fnow(:,:,1), sf(jp_wndj)%fnow(:,:,1), sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1), sf(jp_mslp)%fnow(:,:,1), u_ice, v_ice, tm_su )
      !$acc data present(ptaux_ai,ptauy_ai,pCHi,pCEi,t_air_zu,q_air_zu,sf,u_ice,v_ice,tm_su)
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_sbc_tau: Surface boundary condition for sea ice (momentum)'
         WRITE(numout,*)'~~~~~~~~~~~~~~~'
      ENDIF
      !
      SELECT CASE( ksbc )
         !
      CASE( jp_blk )
         !
         CALL blk_ice_1( sf(jp_wndi)%fnow(:,:,1), sf(jp_wndj)%fnow(:,:,1), sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1), &  ! #LB: `sf(jp_tair)` & `sf(jp_humi)`
            &            sf(jp_mslp)%fnow(:,:,1), u_ice, v_ice, tm_su    ,                              &   ! inputs            !     have been updated at 10m in blk_oce_1
            &            pCHi=pCHi, pCEi=pCEi, putaui=ptaux_ai, pvtaui=ptauy_ai            )                ! outputs
         !
         !
         !
      CASE( jp_abl ) !
         !! This call of `blk_ice_1` is solely done to obtain `CHi` & `CEi` (even if right now they are constant and that's therefore useless)!
         !!  => `blk_ice_1` is also called into `sbc_abl()@sbcabl.F90` to obtaine Qsens, Evap & CDi*U !!
         !!  `ptaux_ai` & `ptauy_ai` are computed in ablmod
         CALL blk_ice_1(  u_abl(:,:,2,nt_n)      ,  v_abl(:,:,2,nt_n) ,      &  !   <<= in
            &            tq_abl(:,:,2,nt_n,jp_ta), tq_abl(:,:,2,nt_n,jp_qa), &  !   <<= in
            &            sf(jp_mslp)%fnow(:,:,1), u_ice, v_ice, tm_su    ,   &  !   <<= in
            &            pCHi=pCHi, pCEi=pCEi )                                 !   <<= out

#if defined key_verbose
         IF(lwp) PRINT *, '* LOLO: `t_air_zu, q_air_zu` overwriten with `tq_abl` @ level=2!'
#endif
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               t_air_zu(ji,jj) = tq_abl(ji,jj,2,nt_n,jp_ta)
               q_air_zu(ji,jj) = tq_abl(ji,jj,2,nt_n,jp_qa)
            END DO
         END DO
         !$acc end parallel loop
         !
#if defined _ABLDBG
         CALL TRDBG( 'ice_sbc_tau: `blk_ice_1:out`',  'pCHi, pCEi', pCHi, pCEi )
#endif
         !
         !
         !
         ! CASE( jp_cpl_atm )   ;    CALL sbc_cpl_ice_tau( ptaux_ai , ptauy_ai )   ! Coupled      formulation LOLO: coupled with atmosphere, not ocean
         !
         !
         !
      END SELECT

#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( 'ice_sbc_tau', ptaux_ai, ptauy_ai )
#else
      CALL lbc_lnk(     'ice_sbc_tau', ptaux_ai,'T',-1._wp, ptauy_ai,'T',-1._wp )
#endif

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_sbc_tau')
      !
   END SUBROUTINE ice_sbc_tau



   SUBROUTINE ice_sbc_flx( kt, ksbc, pCHi, pCEi )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_sbc_flx  ***
      !!
      !! ** Purpose : provide surface boundary condition for sea ice (flux)
      !!
      !! ** Action  : It provides the following fields used in sea ice model:
      !!                emp_oce , emp_ice                        = E-P over ocean and sea ice                    [Kg/m2/s]
      !!                sf(jp_snow)%fnow(:,:,1)                                = solid precipitation                           [Kg/m2/s]
      !!                evap_ice                                 = sublimation (<0 when ice losing FW to atmo)   [Kg/m2/s]
      !!                qsr_ice , qns_ice                        = solar & non solar heat flux over ice          [W/m2]
      !!                dqns_ice                                 = non solar  heat sensistivity                  [W/m2]
      !!                qemp_oce, qemp_ice, qprec_ice            = sensible heat (associated with evap & precip) [W/m2]
      !!            + these fields
      !!                qsb_ice_bot                              = sensible heat at the ice bottom               [W/m2]
      !!                fhld, qlead                              = heat budget in the leads                      [W/m2]
      !!            + some fields that are not used outside this module:
      !!                qla_ice                                  = latent heat flux over ice                     [W/m2]
      !!                sf(jp_prcp)%fnow(:,:,1)                                = total  precipitation                          [Kg/m2/s]
      !!                alb_ice                                  = albedo above sea ice
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time step
      INTEGER, INTENT(in) ::   ksbc   ! flux formulation (user defined, bulk or Pure Coupled)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pCHi    ! sensible heat transfer coefficient [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pCEi    ! evap/sublim. transfer coefficient  [-]
      !!--------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc_flx')
      !%acc data present( pCHi, pCEi, t_su, h_i, h_s, alb_ice, sf(jp_tair)%fnow(:,:,1),sf(jp_humi)%fnow(:,:,1),sf(jp_mslp)%fnow(:,:,1),sf(jp_dqlw)%fnow(:,:,1),sf(jp_prcp)%fnow(:,:,1),sf(jp_snow)%fnow(:,:,1) )
      !$acc data present(pCHi,pCEi,t_su,h_i,h_s,alb_ice,sf)

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_sbc_flx: Surface boundary condition for sea ice (flux)'
         WRITE(numout,*)'~~~~~~~~~~~~~~~'
      ENDIF
      !                     !== ice albedo ==!
      !IF( ln_pnd_alb ) THEN
      !   CALL ice_alb_pnd( t_su, h_i, h_s, a_ip_eff, h_ip, alb_ice )
      !ELSE
      CALL ice_alb(     t_su, h_i, h_s,                 alb_ice )
      !ENDIF

      !
      !SELECT CASE( ksbc )   !== fluxes over sea ice ==!
      !CASE( jp_blk, jp_abl )      !--- bulk formulation & ABL formulation

#if defined _ABLDBG
      CALL TRDBG( 'ice_sbc_flx@icesbc.F90: `blk_ice_2:in`',  't_su, h_s, h_i, alb_ice', t_su, h_s, h_i, alb_ice )
      CALL TRDBG( 'ice_sbc_flx@icesbc.F90: `blk_ice_2:in`',  'sf(jp_tair)%fnow(:,:,1), sf(jp_mslp)%fnow(:,:,1), sf(jp_dqlw)%fnow(:,:,1)', sf(jp_tair)%fnow(:,:,1), sf(jp_mslp)%fnow(:,:,1), sf(jp_dqlw)%fnow(:,:,1) )
      CALL TRDBG( 'ice_sbc_flx@icesbc.F90: `blk_ice_2:in`',  'sf(jp_prcp)%fnow(:,:,1), sf(jp_snow)%fnow(:,:,1)', sf(jp_prcp)%fnow(:,:,1), sf(jp_snow)%fnow(:,:,1) )
      CALL TRDBG( 'ice_sbc_flx@icesbc.F90: `blk_ice_2:in`',  'pCHi, pCEi', pCHi, pCEi )
#endif

      CALL blk_ice_2( t_su, h_s, h_i, alb_ice, sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1),       &
         &            sf(jp_mslp)%fnow(:,:,1), sf(jp_dqlw)%fnow(:,:,1), sf(jp_prcp)%fnow(:,:,1), sf(jp_snow)%fnow(:,:,1), &
         &            pCHi, pCEi )

      !                        !    compute conduction flux and surface temperature (as in Jules surface module)
      IF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN
#if defined _OPENACC || defined _OPENMP
         CALL ctl_stop( 'ice_sbc_flx: routine `blk_ice_qcn` not ported yet to GPU!!! (`ln_cndflx .AND. .NOT.ln_cndemulate`)!' )
#endif
         CALL blk_ice_qcn( ln_virtual_itd, t_su, t_bo, h_s, h_i )
      ENDIF

      IF( ln_icethd )  CALL ice_flx_other()

      !$acc end data
      IF( ln_timing )  CALL timing_stop('ice_sbc_flx')
      !
   END SUBROUTINE ice_sbc_flx



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
      REAL(wp) ::   zfric_u, zqld, zqfr, zqfr_neg, zqfr_pos, zu_io, zv_io, zu_iom1, zv_iom1, zz1, zz2, zmsk
      REAL(wp), PARAMETER ::   zfric_umin = 0._wp       ! lower bound for the friction velocity (cice value=5.e-04)
      REAL(wp), PARAMETER ::   zch        = 0.0057_wp   ! heat transfer coefficient
      !REAL(wp), DIMENSION(jpi,jpj) ::  zfric, zvel      ! ice-ocean velocity (m/s) and frictional velocity (m2/s2)
      ! xtmp1 -> zfric
      ! xtmp2 -> zvel
      !!-----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_flx_other')
      !$acc data present( ssu_m, ssv_m, u_ice, v_ice, qsb_ice_bot, fhld, qlead, xtmp1, xtmp2 )
      !
      ! computation of friction velocity at T points
      IF( ln_icedyn ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               zmsk = xmskt(ji,jj)
               zu_io   = u_ice(ji  ,jj  ) - ssu_m(ji  ,jj  )
               zu_iom1 = u_ice(ji-1,jj  ) - ssu_m(ji-1,jj  )
               zv_io   = v_ice(ji  ,jj  ) - ssv_m(ji  ,jj  )
               zv_iom1 = v_ice(ji  ,jj-1) - ssv_m(ji  ,jj-1)
               xtmp1(ji,jj) = rn_Cd_io * ( 0.5_wp * ( zu_io*zu_io + zu_iom1*zu_iom1 + zv_io*zv_io + zv_iom1*zv_iom1 ) ) * zmsk
               !
               zz1 = u_ice(ji-1,jj  ) + u_ice(ji,jj)
               zz2 = v_ice(ji  ,jj-1) + v_ice(ji,jj)
               xtmp2(ji,jj) = 0.5_wp * SQRT( zz1*zz1 + zz2*zz2 ) * zmsk
            END DO
         END DO
         !$acc end parallel loop
         !
      ELSE      !  if no ice dynamics => transfer directly the atmospheric stress to the ocean
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               xtmp1(ji,jj) = r1_rho0 * SQRT( 0.5_wp *  &
                  &                         (  utau(ji,jj) * utau(ji,jj) + utau(ji-1,jj) * utau(ji-1,jj)   &
                  &                          + vtau(ji,jj) * vtau(ji,jj) + vtau(ji,jj-1) * vtau(ji,jj-1) ) ) * xmskt(ji,jj)
               xtmp2(ji,jj) = 0._wp
            END DO
         END DO
         !$acc end parallel loop
         !
      ENDIF


#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( 'icesbc', xtmp1, xtmp2 )
#else
      CALL lbc_lnk(     'icesbc', xtmp1, 'T',  1.0_wp, xtmp2, 'T', 1.0_wp )
#endif

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
            zfric_u            = MAX( SQRT( xtmp1(ji,jj) ), zfric_umin )
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
            IF(  xtmp2(ji,jj) <= 5.e-04_wp .AND. at_i(ji,jj) >= rn_amax-epsi06 )   qlead(ji,jj) = 0._wp
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
      READ_NML_CFG(numnam_ice,namsbc)
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


   SUBROUTINE ice_sbc_wri( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_sbc_wri ***
      !!
      !! ** Purpose :   save all possible fluxes involved in sea-ice
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_sbc_wri')


      !*****************************
      !  A  Sea-ice to Air fluxes
      !*****************************

      ! A.1 Bulk approach
      ! ~~~~~~~~~~~~~~~~~
      !IF( iom_use('dt_ice_air_zu') ) THEN
      !   !$acc update self( tm_su, t_air_zu )
      !   CALL iom_put("dt_ice_air_zu", (tm_su-t_air_zu)*xmskt * REAL(kmsk_ice_t,wp))
      !ENDIF
      !IF( iom_use('dt_ice_air_zt') ) THEN
      !   !$acc update self( tm_su, sf(jp_tair)%fnow(:,:,1) )
      !   CALL iom_put("dt_ice_air_zt", (tm_su-sf(jp_tair)%fnow(:,:,1))*xmskt * REAL(kmsk_ice_t,wp))
      !ENDIF
      !  ==> still into `blk_ice_2@sbcblk.F90`
      !
      IF( iom_use('Cd_ice') ) THEN
         !$acc update self( CD_ice )
         CALL iom_put("Cd_ice", CD_ice*1000._wp*xmskt * REAL(kmsk_ice_t,wp))
      ENDIF
      IF( iom_use('Ce_ice') ) THEN
         !$acc update self( CE_ice )
         CALL iom_put("Ce_ice", CE_ice*1000._wp*xmskt* REAL(kmsk_ice_t,wp))
      ENDIF
      IF( iom_use('Ch_ice') ) THEN
         !$acc update self( CH_ice )
         CALL iom_put("Ch_ice", CH_ice*1000._wp*xmskt* REAL(kmsk_ice_t,wp))
      ENDIF
      !
      ! A.2 Sea-ice to Air HEAT fluxes (W/m^2)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF( iom_use('qsr_ice') ) THEN
         !$acc update self( qsr_ice )
         CALL iom_put( 'qsr_ice', SUM(   qsr_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use('qlw_ice') ) THEN
         !$acc update self( qlw_ice )
         CALL iom_put( 'qlw_ice', SUM(   qlw_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use('qla_ice') ) THEN
         !$acc update self( qla_ice )
         CALL iom_put( 'qla_ice', SUM( - qla_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) ) !#LB: sign consistent with what's done for ocean
      ENDIF
      IF( iom_use('qsb_ice') ) THEN
         !$acc update self( qsb_ice )
         CALL iom_put( 'qsb_ice', SUM( - qsb_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) ) !#LB: sign consistent with what's done for ocean
      ENDIF
      IF( iom_use('qns_ice') ) THEN
         !$acc update self( qns_ice )
         CALL iom_put( 'qns_ice', SUM(   qns_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use("qemp_ice") ) THEN
         !$acc update self( qemp_ice )
         CALL iom_put( "qemp_ice"   , (  qemp_ice                                            ) * REAL(kmsk_ice_t,wp) ) ! Downward Heat Flux from E-P over ice
      ENDIF
      IF( iom_use("qt_ice") ) THEN
         !$acc update self( qns_ice, qsr_ice, qemp_ice )
         CALL iom_put( "qt_ice" , (SUM( ( qns_ice + qsr_ice ) * a_i_b, dim=3 )     + qemp_ice) * REAL(kmsk_ice_t,wp) )
      ENDIF
      !
      ! A.3 Sea-ice to Air FRESHWATER fluxes (kg/m2/s)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF( iom_use("evap_ice") ) THEN
         !$acc update self( evap_ice )
         CALL iom_put( "evap_ice", SUM( evap_ice * a_i_b, dim=3 )* REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use("emp_ice") ) THEN
         !$acc update self( emp_ice )
         CALL iom_put( "emp_ice" , emp_ice  )   ! emp over ice   (taking into account the snow blown away from the ice)
      ENDIF



      ! A.4 Sea-ice MOMENTUM fluxes (top & bottom) (N/m2)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! --- ice-atm. stress: No `A` involved !!!
      IF( iom_use('taux_ai_t') ) THEN
         !$acc update self(taux_ai_t)
         CALL iom_put( 'taux_ai_t'   , taux_ai_t * REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use('tauy_ai_t') ) THEN
         !$acc update self(tauy_ai_t)
         CALL iom_put( 'tauy_ai_t'   , tauy_ai_t * REAL(kmsk_ice_t,wp) )
      ENDIF
      IF( iom_use('taum_ai') ) THEN
         IF(.NOT. iom_use('taux_ai_t') ) THEN
            !$acc update self(taux_ai_t)
         ENDIF
         IF(.NOT. iom_use('tauy_ai_t') ) THEN
            !$acc update self(tauy_ai_t)
         ENDIF
         CALL iom_put( 'taum_ai' , SQRT(taux_ai_t*taux_ai_t + tauy_ai_t*tauy_ai_t)* REAL(kmsk_ice_t,wp) )
      END IF
      !
      ! --- oce-ice stress under sea-ice (with sign as seen from ice): No `A` involved !!!
      IF( iom_use('taux_oi_u') ) THEN
         IF( sn_loc_vct_tau /= 'C' ) CALL ctl_stop( 'STOP', 'ice_sbc_wri : we have "sn_loc_vct_tau/=C", cannot save "taux_oi_u"' )
         !$acc update self( taux_oi_u )
         CALL iom_put( 'taux_oi_u' , taux_oi_u * REAL(kmsk_ice_u, wp) )
      ENDIF
      IF( iom_use('tauy_oi_v') ) THEN
         IF( sn_loc_vct_tau /= 'C' ) CALL ctl_stop( 'STOP', 'ice_sbc_wri : we have "sn_loc_vct_tau/=C", cannot save "tauy_oi_v"' )
         !$acc update self( tauy_oi_v )
         CALL iom_put( 'tauy_oi_v' , tauy_oi_v * REAL(kmsk_ice_v, wp) )
      ENDIF
      !
      !IF( iom_use('taum_oi') ) THEN
      !   IF( ln_damage ) THEN
      !      ztmp1(:,:) = V_oce(:,:,1) - u_ice(:,:)  ! dU @U
      !      ztmp2(:,:) = V_oce(:,:,4) - vUice(:,:)  ! dV @U
      !      ztmp3(:,:) = rho0*rn_Cd_io * SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 ) * ztmp1(:,:) * REAL(kmsk_ice_u(:,:),wp)
      !      IF( iom_use('taux_oi_u') ) CALL iom_put( 'taux_oi_u' , ztmp3 )
      !   ELSE
      !      CALL ctl_stop( 'STOP', 'ice_wri: FIXME! Add `taux_oi_u` for EVP' )
      !   ENDIF
      !   ztmp1(2:jpi,:) = 0.5_wp * ( ztmp3(2:jpi,:) + ztmp3(1:jpi-1,:) )
      !   ztmp2(:,2:jpj) = 0.5_wp * ( ztmp4(:,2:jpj) + ztmp4(:,1:jpj-1) )
      !   CALL iom_put( 'taum_oi' , SQRT(ztmp1*ztmp1 + ztmp2*ztmp2) )
      !ENDIF



      !! Heat fluxes from the system "ocean+ice" to the atmosphere
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF( iom_use("qns2atm") ) THEN
         !! Non-solar heat flux to the atmosphere
         !!  ==> same as `zqns_tot` in `ice_update_flx@iceupdate.F90` => `-zqns_tot` because "to the atmosphere"!
         !$acc update self(at_i_b, a_i_b, qns_ice, qns_oce, qemp_ice, qemp_oce)
         CALL iom_put( "qns2atm" , -1._wp * ( (1._wp-at_i_b(:,:)) * qns_oce(:,:)   +  SUM( a_i_b(:,:,:)*qns_ice(:,:,:), dim=3 )  +  qemp_ice(:,:) + qemp_oce(:,:) ) * xmskt(:,:) )
      ENDIF

      IF( iom_use("qsr2atm") ) THEN
         !! Solar heat flux to the atmosphere (from ocean+ice)
         !!  ==> same as `zqsr_tot` in `ice_update_flx@iceupdate.F90` => `-zqsr_tot` because "to the atmosphere"!
         !$acc update self(at_i_b, a_i_b, qsr_ice, qsr_oce, qemp_ice, qemp_oce)
         CALL iom_put( "qsr2atm" , -1._wp * ( (1._wp-at_i_b(:,:)) * qsr_oce(:,:)   +  SUM( a_i_b(:,:,:)*qsr_ice(:,:,:), dim=3 ) ) * xmskt(:,:) )
      ENDIF
      
      IF( iom_use("qt2atm") ) THEN
         !! Net heat flux to the atmosphere (this is `-qt_atm_oi`)
         !$acc update self( qt_atm_oi )
         CALL iom_put( "qt2atm" , -1._wp * qt_atm_oi(:,:) * xmskt(:,:) )
      ENDIF

      

      ! --- salt fluxes [kg/m2/s] --- !
      !                           ! sfxice =  sfxbog + sfxbom + sfxsum + sfxsni + sfxopw + sfxres + sfxdyn + sfxbri + sfxsub + sfxlam
      IF( iom_use("sfxice"  ) )   CALL iom_put( "sfxice", sfx     * 1.e-03 )   ! salt flux from total ice growth/melt
      IF( iom_use("sfxbog"  ) )   CALL iom_put( "sfxbog", sfx_bog * 1.e-03 )   ! salt flux from bottom growth
      IF( iom_use("sfxbom"  ) )   CALL iom_put( "sfxbom", sfx_bom * 1.e-03 )   ! salt flux from bottom melting
      IF( iom_use("sfxsum"  ) )   CALL iom_put( "sfxsum", sfx_sum * 1.e-03 )   ! salt flux from surface melting
      IF( iom_use("sfxlam"  ) )   CALL iom_put( "sfxlam", sfx_lam * 1.e-03 )   ! salt flux from lateral melting
      IF( iom_use("sfxsni"  ) )   CALL iom_put( "sfxsni", sfx_sni * 1.e-03 )   ! salt flux from snow ice formation
      IF( iom_use("sfxopw"  ) )   CALL iom_put( "sfxopw", sfx_opw * 1.e-03 )   ! salt flux from open water formation
      IF( iom_use("sfxdyn"  ) )   CALL iom_put( "sfxdyn", sfx_dyn * 1.e-03 )   ! salt flux from ridging rafting
      IF( iom_use("sfxbri"  ) )   CALL iom_put( "sfxbri", sfx_bri * 1.e-03 )   ! salt flux from brines
      IF( iom_use("sfxres"  ) )   CALL iom_put( "sfxres", sfx_res * 1.e-03 )   ! salt flux from undiagnosed processes
      IF( iom_use("sfxsub"  ) )   CALL iom_put( "sfxsub", sfx_sub * 1.e-03 )   ! salt flux from sublimation

      ! --- mass fluxes [kg/m2/s] --- !
      !CALL iom_put( "emp_oce", emp_oce )   ! emp over ocean (taking into account the snow blown away from the ice)


      !                           ! vfxice = vfxbog + vfxbom + vfxsum + vfxsni + vfxopw + vfxdyn + vfxres + vfxlam + vfxpnd
      CALL iom_put( "vfxice"    , wfx_ice     )   ! mass flux from total ice growth/melt
      CALL iom_put( "vfxbog"    , wfx_bog     )   ! mass flux from bottom growth
      CALL iom_put( "vfxbom"    , wfx_bom     )   ! mass flux from bottom melt
      CALL iom_put( "vfxsum"    , wfx_sum     )   ! mass flux from surface melt
      CALL iom_put( "vfxlam"    , wfx_lam     )   ! mass flux from lateral melt
      CALL iom_put( "vfxsni"    , wfx_sni     )   ! mass flux from snow-ice formation
      CALL iom_put( "vfxopw"    , wfx_opw     )   ! mass flux from growth in open water
      CALL iom_put( "vfxdyn"    , wfx_dyn     )   ! mass flux from dynamics (ridging)
      CALL iom_put( "vfxres"    , wfx_res     )   ! mass flux from undiagnosed processes
      CALL iom_put( "vfxpnd"    , wfx_pnd     )   ! mass flux from melt ponds
      CALL iom_put( "vfxsub"    , wfx_ice_sub )   ! mass flux from ice sublimation (ice-atm.)
      CALL iom_put( "vfxsub_err", wfx_err_sub )   ! "excess" of sublimation sent to ocean

      !                            ! vfxsnw = vfxsnw_sni + vfxsnw_dyn + vfxsnw_sum
      CALL iom_put( "vfxsnw"     , wfx_snw     )   ! mass flux from total snow growth/melt
      CALL iom_put( "vfxsnw_sum" , wfx_snw_sum )   ! mass flux from snow melt at the surface
      CALL iom_put( "vfxsnw_sni" , wfx_snw_sni )   ! mass flux from snow melt during snow-ice formation
      CALL iom_put( "vfxsnw_dyn" , wfx_snw_dyn )   ! mass flux from dynamics (ridging)
      CALL iom_put( "vfxsnw_sub" , wfx_snw_sub )   ! mass flux from snow sublimation (ice-atm.)
      CALL iom_put( "vfxsnw_pre" , wfx_spr     )   ! snow precip

      IF( iom_use("qtr_ice_bot") ) CALL iom_put( "qtr_ice_bot", (SUM( qtr_ice_bot * a_i_b, dim=3 ) ) * REAL(kmsk_ice_t,wp) ) !     solar flux transmitted thru ice
      IF( iom_use("qtr_ice_top") ) CALL iom_put( "qtr_ice_top", (SUM( qtr_ice_top * a_i_b, dim=3 ) ) * REAL(kmsk_ice_t,wp) ) !     solar flux transmitted thru ice surface
      IF( iom_use("qt_oce_ai"  ) ) CALL iom_put( "qt_oce_ai"  , qt_oce_ai                            * REAL(kmsk_ice_t,wp) ) ! total heat flux at the ocean   surface: interface oce-(ice+atm)
      IF( iom_use("qt_atm_oi"  ) ) CALL iom_put( "qt_atm_oi"  , qt_atm_oi                            * REAL(kmsk_ice_t,wp) ) ! total heat flux at the oce-ice surface: interface atm-(ice+oce)


      ! heat fluxes from ice transformations
      !                            ! hfxdhc = hfxbog + hfxbom + hfxsum + hfxopw + hfxdif + hfxsnw - ( hfxthd + hfxdyn + hfxres + hfxsub + hfxspr )
      CALL iom_put ("hfxbog"     , hfx_bog     )   ! heat flux used for ice bottom growth
      CALL iom_put ("hfxbom"     , hfx_bom     )   ! heat flux used for ice bottom melt
      CALL iom_put ("hfxsum"     , hfx_sum     )   ! heat flux used for ice surface melt
      CALL iom_put ("hfxopw"     , hfx_opw     )   ! heat flux used for ice formation in open water
      CALL iom_put ("hfxdif"     , hfx_dif     )   ! heat flux used for ice temperature change
      CALL iom_put ("hfxsnw"     , hfx_snw     )   ! heat flux used for snow melt
      CALL iom_put ("hfxerr"     , hfx_err_dif )   ! heat flux error after heat diffusion

      ! heat fluxes associated with mass exchange (freeze/melt/precip...)
      CALL iom_put ("hfxthd"     , hfx_thd     )   !
      CALL iom_put ("hfxdyn"     , hfx_dyn     )   !
      CALL iom_put ("hfxres"     , hfx_res     )   !
      CALL iom_put ("hfxsub"     , hfx_sub     )   !
      CALL iom_put ("hfxspr"     , hfx_spr     )   ! Heat flux from snow precip heat content

      ! other heat fluxes
      IF( iom_use("hfxsensib"  ) )   CALL iom_put( "hfxsensib"  ,      qsb_ice_bot * at_i_b         )   ! Sensible oceanic heat flux
      IF( iom_use("hfxcndbot"  ) )   CALL iom_put( "hfxcndbot"  , SUM( qcn_ice_bot * a_i_b, dim=3 ) )   ! Bottom conduction flux
      IF( iom_use("hfxcndtop"  ) )   CALL iom_put( "hfxcndtop"  , SUM( qcn_ice_top * a_i_b, dim=3 ) )   ! Surface conduction flux
      IF( iom_use("hfxmelt"    ) )   CALL iom_put( "hfxmelt"    , SUM( qml_ice     * a_i_b, dim=3 ) )   ! Surface melt flux
      IF( iom_use("hfxldmelt"  ) )   CALL iom_put( "hfxldmelt"  ,      fhld        * at_i_b         )   ! Heat in lead for ice melting
      IF( iom_use("hfxldgrow"  ) )   CALL iom_put( "hfxldgrow"  ,      qlead       * r1_Dt_ice      )   ! Heat in lead for ice growth

      IF( ln_timing         )   CALL timing_stop   ('ice_sbc_wri')                                       ! timing
      !
   END SUBROUTINE ice_sbc_wri


   !!======================================================================
END MODULE icesbc
