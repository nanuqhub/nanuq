MODULE iceupdate
   !!======================================================================
   !!                       ***  MODULE iceupdate   ***
   !!  Sea-ice :   computation of the flux at the sea ice/ocean interface
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_update_init  : initialisation
   !!   ice_update_flx   : updates mass, heat and salt fluxes at the ocean surface
   !!   ice_update_tau   : update i- and j-stresses, and its modulus at the ocean surface
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean domain
   USE par_ice
   USE ice            ! sea-ice: variables
   USE sbc_ice        ! Surface boundary condition: ice   fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE oss_nnq , ONLY : frq_m
   USE icealb         ! sea-ice: albedo parameters
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu     ! lateral boundary conditions (or mpp links)
#else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#endif
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_update_init   ! called by ice_init
   PUBLIC   ice_update_flx    ! called by ice_stp
   PUBLIC   ice_update_tau    ! called by ice_stp

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! $Id: iceupdate.F90 15385 2021-10-15 13:52:48Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_update_flx( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_update_flx ***
      !!
      !! ** Purpose :   Update the surface ocean boundary condition for heat
      !!                salt and mass over areas where sea-ice is non-zero
      !!
      !! ** Action  : - computes the heat and freshwater/salt fluxes
      !!                at the ice-ocean interface.
      !!              - Update the ocean sbc
      !!
      !! ** Outputs : - qsr     : sea heat flux:     solar
      !!              - qns     : sea heat flux: non solar
      !!              - emp     : freshwater budget: volume flux
      !!              - sfx     : salt flux
      !!              - t_su    : sea-ice surface temperature
      !!              - alb_ice : sea-ice albedo (recomputed only for coupled mode)
      !!
      !! References : Goosse, H. et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90.
      !!              Tartinville et al. 2001 Ocean Modelling, 3, 95-108.
      !!              These refs are now obsolete since everything has been revised
      !!              The ref should be Rousset et al., 2015
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   zqsr             ! New solar flux received by the ocean
      REAL(wp) ::   zA_b, z1mA_b, zsum, zqns_tot, zqsr_tot
      !REAL(wp), DIMENSION(jpi,jpj) ::   z2d    ! 2D workspace for IOM stuff
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_update_flx')
      !$acc data present( a_i_b, alb_ice, at_i_b )

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_flx: update fluxes (mass, salt and heat) at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

#if defined _OPENACC || defined _OPENMP
      IF( ln_cndflx )  CALL ctl_stop( 'STOP', 'ice_update_flx : adapt option `ln_cndflx` for GPU!')
#endif

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zA_b    = at_i_b(ji,jj)
            z1mA_b  = 1._wp - zA_b

            ! Net heat flux on top of the ice-ocean (W.m-2)
            !----------------------------------------------
            ! MIND => when the system "ice-ocean" is receiving heat from the atmosphere => `flux > 0` !
            
            zsum = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               zsum = zsum + a_i_b(ji,jj,jl) * qns_ice(ji,jj,jl)
            END DO
            zqns_tot = z1mA_b * qns_oce(ji,jj)   +   zsum   +   qemp_ice(ji,jj) + qemp_oce(ji,jj) ! saved as `-qns2atm` by XIOS (`ice_sbc_wri@icesbc.F90`)

            zsum = 0._wp
            !$acc loop seq
            DO jl = 1, jpl
               zsum = zsum + a_i_b(ji,jj,jl) * qsr_ice(ji,jj,jl)
            END DO
            zqsr_tot = z1mA_b * qsr_oce(ji,jj)   +   zsum              ! saved as `-qsr2atm` by XIOS (`ice_sbc_wri@icesbc.F90`)

            !IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
            !   qt_atm_oi(ji,jj) = z1mA_b * ( qns_oce(ji,jj) + qsr_oce(ji,jj) ) + qemp_oce(ji,jj) + &
            !      &             SUM( a_i_b(ji,jj,1:jpl) * ( qcn_ice(ji,jj,1:jpl) + qml_ice(ji,jj,1:jpl) + qtr_ice_top(ji,jj,1:jpl) ), dim=3 ) + qemp_ice(ji,jj)
            !ELSE                   ! ice-atm interface = solar and non-solar fluxes
            qt_atm_oi(ji,jj) = zqns_tot + zqsr_tot     ! saved as `-qt2atm` by XIOS (`ice_sbc_wri@icesbc.F90`)
            !ENDIF

            ! --- case we bypass ice thermodynamics --- !
            IF( .NOT. ln_icethd ) THEN   ! we suppose ice is impermeable => ocean is isolated from atmosphere
#if defined key_verbose
               IF(lwp) PRINT *, 'LOLO: `ice_update_flx@iceupdate.F90` => `qt_oce_ai` bypasses thermo!, kt =', kt
#endif
               qt_atm_oi(ji,jj)   = z1mA_b * ( qns_oce(ji,jj) + qsr_oce(ji,jj) ) + qemp_oce(ji,jj)
               qt_oce_ai(ji,jj)   = z1mA_b *   qns_oce(ji,jj)                    + qemp_oce(ji,jj)
               emp_ice  (ji,jj)   = 0._wp
               qemp_ice (ji,jj)   = 0._wp
            ENDIF


            ! Solar heat flux reaching the ocean (max) = zqsr (W.m-2)
            !---------------------------------------------------
            !IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
            !   zqsr = z1mA_b * qsr_oce(ji,jj) + SUM( a_i_b (ji,jj,:) * qtr_ice_bot(ji,jj,:) )
            !ELSE                   ! ice-atm interface = solar and non-solar fluxes
            zsum = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zsum = zsum + a_i_b(ji,jj,jl) * ( qsr_ice(ji,jj,jl) - qtr_ice_bot(ji,jj,jl) )
            END DO
            zqsr = zqsr_tot - zsum
            !ENDIF

            ! Total heat flux reaching the ocean = qt_oce_ai (W.m-2)
            !---------------------------------------------------
            IF( ln_icethd ) THEN
               qt_oce_ai(ji,jj) = qt_atm_oi(ji,jj) - hfx_sum(ji,jj) - hfx_bom(ji,jj) - hfx_bog(ji,jj) &
                  &                                - hfx_dif(ji,jj) - hfx_opw(ji,jj) - hfx_snw(ji,jj) &
                  &                                + hfx_thd(ji,jj) + hfx_dyn(ji,jj) + hfx_res(ji,jj) &
                  &                                + hfx_sub(ji,jj) + hfx_spr(ji,jj)
               !
            ENDIF

            ! New qsr and qns used to compute the oceanic heat flux at the next time step
            !----------------------------------------------------------------------------
            ! if warming and some ice remains, then we suppose that the whole solar flux has been consumed to melt the ice
            ! else ( cooling or no ice left ), then we suppose that     no    solar flux has been consumed
            !
            IF( fhld(ji,jj) > 0._wp .AND. at_i(ji,jj) > 0._wp ) THEN   !-- warming and some ice remains
               zsum = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zsum = zsum + a_i_b(ji,jj,jl) * qtr_ice_bot(ji,jj,jl)
               END DO
               !                                        solar flux transmitted thru the 1st level of the ocean (i.e. not used by sea-ice)
               qsr(ji,jj) = z1mA_b * qsr_oce(ji,jj) * ( 1._wp - frq_m(ji,jj) ) &
                                !                                   + solar flux transmitted thru ice and the 1st ocean level (also not used by sea-ice)
                  &             + zsum * ( 1._wp - frq_m(ji,jj) )
               !
            ELSE                                                       !-- cooling or no ice left
               qsr(ji,jj) = zqsr
            ENDIF
            !
            ! the non-solar is simply derived from the solar flux
            !IF(lwp .AND. (ji==10 .AND. jj==10)) PRXNT *, 'LOLO: `ice_update_flx@iceupdate.F90`: qns = qt_oce_ai - qsr ,  kt =', kt
            qns(ji,jj) = qt_oce_ai(ji,jj) - qsr(ji,jj)

            ! Mass flux at the atm. surface
            !-----------------------------------
            wfx_sub(ji,jj) = wfx_snw_sub(ji,jj) + wfx_ice_sub(ji,jj)

            ! Mass flux at the ocean surface
            !------------------------------------
            ! ice-ocean  mass flux
            wfx_ice(ji,jj) = wfx_bog(ji,jj) + wfx_bom(ji,jj) + wfx_sum(ji,jj) + wfx_sni(ji,jj)   &
               &           + wfx_opw(ji,jj) + wfx_dyn(ji,jj) + wfx_res(ji,jj) + wfx_lam(ji,jj)

            ! snw-ocean mass flux
            wfx_snw(ji,jj) = wfx_snw_sni(ji,jj) + wfx_snw_dyn(ji,jj) + wfx_snw_sum(ji,jj)

            ! total mass flux at the ocean/ice interface
            !IF( ln_pnd ) THEN
            !   fmmflx(ji,jj) =                - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_pnd(ji,jj) - wfx_err_sub(ji,jj) ! ice-ocean mass flux saved at least for biogeochemical model
            !   emp   (ji,jj) = emp_oce(ji,jj) - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_pnd(ji,jj) - wfx_err_sub(ji,jj) ! atm-ocean + ice-ocean mass flux
            !ELSE
            fmmflx(ji,jj) = - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_err_sub(ji,jj)  ! ice-ocean mass flux: `fmmflx>0` => LOSS for the liquid ocean
            emp   (ji,jj) = emp_oce(ji,jj) + fmmflx(ji,jj)              ! atm-ocean + ice-ocean mass flux:    `emp>0` => LOSS for the liquid ocean
            !ENDIF

            ! Salt flux at the ocean surface
            !------------------------------------------
            sfx(ji,jj) = sfx_bog(ji,jj) + sfx_bom(ji,jj) + sfx_sum(ji,jj) + sfx_sni(ji,jj) + sfx_opw(ji,jj)   &
               &       + sfx_res(ji,jj) + sfx_dyn(ji,jj) + sfx_bri(ji,jj) + sfx_sub(ji,jj) + sfx_lam(ji,jj)

            ! Mass of snow and ice per unit area
            !----------------------------------------
            snwice_mass_b(ji,jj) = snwice_mass(ji,jj)       ! save mass from the previous ice time step
            !                                               ! new mass per unit area
            !IF( ln_pnd ) THEN
            !   snwice_mass  (ji,jj) = xmskt(ji,jj) * ( rhos * vt_s(ji,jj) + rhoi * vt_i(ji,jj) + rhow * (vt_ip(ji,jj) + vt_il(ji,jj)) )
            !ELSE
            snwice_mass  (ji,jj) = xmskt(ji,jj) * ( rhos * vt_s(ji,jj) + rhoi * vt_i(ji,jj) )
            !ENDIF
            !                                               ! time evolution of snow+ice mass

         END DO
      END DO
      !$acc end parallel loop

#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( 'ice_update_flx', emp )             ! IMPORTANT
#else
      CALL lbc_lnk(     'ice_update_flx', emp,'T',1._wp )   ! IMPORTANT
#endif

      ! Snow/ice albedo (only if sent to coupler, useless in forced mode)
      !------------------------------------------------------------------
      !IF( ln_pnd_alb ) THEN
      !   CALL ice_alb_pnd( t_su, h_i, h_s, a_ip_eff, h_ip, alb_ice ) ! ice albedo
      !ELSE
      CALL ice_alb(     t_su, h_i, h_s,                  alb_ice ) ! ice albedo
      !ENDIF

      !
      IF( lrst_ice ) THEN                       !* write snwice_mass fields in the restart file
         CALL update_rst( 'WRITE', kt )
      ENDIF

      ! controls
      !---------
      !IF( ln_icediachk      )   CALL ice_cons_final('ice_update_flx')                                       ! conservation
      IF( ln_icectl         )   CALL ice_prt       (kt, iiceprt, jiceprt, 3, 'Final state ice_update') ! prints
      IF( sn_cfctl%l_prtctl )   CALL ice_prt3D     ('ice_update_flx')                                       ! prints

      !$acc end data
      IF( ln_timing         )   CALL timing_stop   ('ice_update_flx')                                       ! timing
      !
   END SUBROUTINE ice_update_flx


   SUBROUTINE ice_update_tau( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_update_tau ***
      !!
      !! ** Purpose : Update the ocean surface stresses due to the ice
      !!
      !! ** Action  : * at each ice time step (every time step):
      !!                - compute the modulus of ice-ocean relative velocity
      !!                  (*rho*Cd) at T-point (C-grid) or I-point (B-grid)
      !!                      ztmod_io = rhoco * | U_ice-U_oce |
      !!                - update the modulus of stress at ocean surface
      !!                      taum = (1-a) * taum + a * ztmod_io * | U_ice-U_oce |
      !!              * at each ocean time step (every kt):
      !!                  compute linearized ice-ocean stresses as
      !!                      Utau = ztmod_io * | U_ice - U_oce |
      !!                using instantaneous current ocean velocity (usually before)
      !!
      !!    NB: - ice-ocean rotation angle no more allowed
      !!        - here we make an approximation: taum is only computed every ice time step
      !!          This avoids mutiple average to pass from T -> U,V grids and next from U,V grids
      !!          to T grid. taum is used in TKE and GLS, which should not be too sensitive to this approximaton...
      !!
      !! ** Outputs : - utau, vtau : surface ocean i- and j-stress updated WITH ice-ocean fluxes
      !!                               (@ U- and V-points, respectively, if `sn_loc_vct_tau=='C'`)
      !!                               (@ T-points, if `sn_loc_vct_tau=='T'`)
      !!              - taum       : modulus of the surface ocean stress (T-point) updated with ice-ocean fluxes
      !!---------------------------------------------------------------------
      INTEGER ,                     INTENT(in) ::   kt               ! ocean time-step index
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   za_tot_u, ztaux_ai_u, zu_t, zmodt   ! local scalar
      REAL(wp) ::   za_tot_v, ztauy_ai_v, zv_t, zrhoco  !   -      -
      REAL(wp) ::   ztaux_oi_u, ztauy_oi_v, ztmod_io
      REAL(wp) ::   za_tot_t, zA, zB, ztaux, ztauy
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_update_tau')
      !$acc data present( utau, vtau, u_ice, v_ice, uVice, vUice, V_oce, taum, taux_oi_u, tauy_oi_v, at_i, au_i, av_i )

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_tau: update stress at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

      zrhoco = rho0 * rn_Cd_io

      !nb = MAX( nn_hls-1, 0 )

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ! 1/ Update modulus of stress received by the liquid ocean @T (weight contributions of air-sea & ice-sea stresses / A )
            zA = at_i(ji,jj)
            !
            zu_t = u_ice(ji,jj) + u_ice(ji-1,jj) - V_oce(ji,jj,1) - V_oce(ji-1,jj,1)  ! 2*(U_ice-U_oce) at T-point
            zv_t = v_ice(ji,jj) + v_ice(ji,jj-1) - V_oce(ji,jj,2) - V_oce(ji,jj-1,2)  ! 2*(V_ice-V_oce) at T-point
            zmodt =  0.25_wp * ( zu_t*zu_t + zv_t*zv_t  )   ! |U_ice-U_oce|^2
            ! Making stress modulus received by the liquid ocean @T ice-aware:
            taum(ji,jj) = (1._wp - zA) * taum(ji,jj) + zA * zrhoco * zmodt

            ! 2/ Update stress vector components received by the liquid ocean @U,V (weight contributions of air-sea & ice-sea stresses / Au,Av )
            zA = u_ice(ji,jj) - V_oce(ji,jj,1)  ! U_ice - U_oce @ U
            zB = vUice(ji,jj) - V_oce(ji,jj,4)  ! V_ice - V_oce @ U
            ztmod_io = SQRT( zA*zA + zB*zB )    ! modulus of ice - oce @ U
            ztaux_oi_u = zrhoco * ztmod_io * zA
            taux_oi_u(ji,jj) = -1._wp * ztaux_oi_u ! reverse sign because `taux_oi_u` is what is felt by the ice, not the ocean !
            !
            zA = uVice(ji,jj) - V_oce(ji,jj,3) ! U_ice - U_oce @ V
            zB = v_ice(ji,jj) - V_oce(ji,jj,2) ! V_ice - V_oce @ V
            ztmod_io = SQRT( zA*zA + zB*zB )   ! modulus of ice - oce @ U
            ztauy_oi_v = zrhoco * ztmod_io * zB
            tauy_oi_v(ji,jj) = -1._wp * ztauy_oi_v ! reverse sign because `taux_oi_u` is what is felt by the ice, not the ocean !
            !
            ! Making stresses received by the liquid ocean @U,V ice-aware:
            IF( k_tau_air_at_T == 1 ) THEN
               ! => both `utau` & `vtau` have to be located at T-point
               !    ==> ice-sea stress needs to be interpolated at T-points (`utau` & `vtau` already @T)
               ztaux = 0.5_wp * ( taux_oi_u(ji,jj) + taux_oi_u(ji-1,jj) ) ! => ztaux is `taux_oi_u` interp @T
               ztauy = 0.5_wp * ( tauy_oi_v(ji,jj) + tauy_oi_v(ji,jj-1) ) ! => ztauy is `tauy_oi_v` interp @T
               utau(ji,jj) = (1._wp - zA) * utau(ji,jj)  +  zA * ztaux ! @T
               vtau(ji,jj) = (1._wp - zA) * vtau(ji,jj)  +  zA * ztauy ! @T
               !
            ELSE
               ! => follow 'C-grid' convention => `utau` @ U-point & `vtau` @ V-point
               za_tot_u = au_i(ji,jj)
               za_tot_v = av_i(ji,jj)
               !    ==> air-sea stress needs to be interpolated at U,V-points
               ztaux = 0.5_wp*(utau(ji,jj) + utau(ji+1,jj)) * (2._wp-umask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji+1,jj)) ! => ztaux is utau interp @U
               ztauy = 0.5_wp*(vtau(ji,jj) + vtau(ji,jj+1)) * (2._wp-vmask(ji,jj,1))*MAX(xmskt(ji,jj),xmskt(ji,jj+1)) ! => ztauy is vtau interp @V
               !
               utau(ji,jj) = (1._wp - za_tot_u) * ztaux  +  za_tot_u * ztaux_oi_u  ! @U
               vtau(ji,jj) = (1._wp - za_tot_v) * ztauy  +  za_tot_v * ztauy_oi_v  ! @V
            ENDIF
         END DO
      END DO
      !$acc end parallel loop


#if defined _OPENACC || defined _OPENMP
      !#LOLOfixme: really ok for periodic LBC to do the same shit for T, U or V like here????
      CALL lbc_lnk_gpu( 'ice_update_tau', taum, utau, vtau, taux_oi_u, tauy_oi_v )            ! lateral boundary condition
#else
      IF( k_tau_air_at_T == 1 ) THEN
#if defined key_verbose
         IF(lwp) PRINT *, ' *** LOLO:ice_update_tau => `utau, vtau` updated at T-points!!!'
#endif
         CALL lbc_lnk(     'ice_update_tau', taum,'T',1._wp, utau,'T',-1._wp,      vtau,'T',-1._wp, &
            &                                           taux_oi_u,'U',-1._wp, tauy_oi_v,'V',-1._wp  )   ! lateral boundary condition
      ELSE
#if defined key_verbose
         IF(lwp) PRINT *, ' *** LOLO:ice_update_tau => `utau, vtau` updated at U- & V-points!!!'
#endif
         CALL lbc_lnk(     'ice_update_tau', taum,'T',1._wp, utau,'U',-1._wp,      vtau,'V',-1._wp, &
            &                                           taux_oi_u,'U',-1._wp, tauy_oi_v,'V',-1._wp  )   ! lateral boundary condition
      ENDIF
#endif

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_update_tau')
      !
   END SUBROUTINE ice_update_tau


   SUBROUTINE ice_update_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_update_init  ***
      !!
      !! ** Purpose :   allocate ice-ocean stress fields and read restarts
      !!                containing the snow & ice mass
      !!
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   zcoefu, zcoefv, zcoeff   ! local scalar
      !!-------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_update_init: ice-ocean stress init'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      !
      CALL update_rst( 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE ice_update_init


   SUBROUTINE update_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_evp_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! 'READ'/'WRITE' flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter   ! local integer
      INTEGER  ::   id1    ! local integer
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id1 = iom_varid( 'update_rst', numrir, 'snwice_mass' , ldstop = .FALSE. )
            !
            IF( id1 > 0 ) THEN                       ! fields exist
               CALL iom_get( 'update_rst', numrir, jpdom_auto, 'snwice_mass'  , snwice_mass   )
               CALL iom_get( 'update_rst', numrir, jpdom_auto, 'snwice_mass_b', snwice_mass_b )
            ELSE                                     ! start from rest
               IF(lwp) WRITE(numout,*) '   ==>>   previous run without snow-ice mass output then set it'
               snwice_mass  (:,:) = xmskt(:,:) * ( rhos * vt_s(:,:) + rhoi * vt_i(:,:) ) ! &
               !&  + rhow * (vt_ip(:,:) + vt_il(:,:))  )
               snwice_mass_b(:,:) = snwice_mass(:,:)
            ENDIF
         ELSE                                   !* Start from rest
            !JC: I think this is useless with what is now done in ice_istate
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest: set the snow-ice mass'
            !IF( ln_pnd ) THEN
            !   snwice_mass(:,:) = xmskt(:,:) * ( rhos*vt_s(:,:) + rhoi*vt_i(:,:) + rhow*(vt_ip(:,:) + vt_il(:,:))  )
            !ELSE
            snwice_mass(:,:) = xmskt(:,:) * ( rhos*vt_s(:,:) + rhoi*vt_i(:,:) )
            !ENDIF
            snwice_mass_b(:,:) = snwice_mass(:,:)
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- update-rst ----'
         iter = kt             ! ice restarts are written at kt == nitrst
         !
         !$acc update self( snwice_mass, snwice_mass_b )
         CALL iom_rstput( iter, nitrst, numriw, 'snwice_mass'  , snwice_mass   )
         CALL iom_rstput( iter, nitrst, numriw, 'snwice_mass_b', snwice_mass_b )
         !
      ENDIF
      !
   END SUBROUTINE update_rst

   !!======================================================================
END MODULE iceupdate
