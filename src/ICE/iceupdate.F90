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
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_update_init   ! called by ice_init
   PUBLIC   ice_update_flx    ! called by ice_stp
   PUBLIC   ice_update_tau    ! called by ice_stp

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
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
      !!              - fr_i    : ice fraction
      !!              - tn_ice  : sea-ice surface temperature
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
      REAL(wp) ::   zsum
      !REAL(wp), DIMENSION(jpi,jpj) ::   z2d    ! 2D workspace for IOM stuff
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_update_flx')
      !$acc data present( a_i_b, alb_ice, at_i, at_i_b, emp, emp_ice, emp_oce, fhld, fmmflx, fr_i, frq_m, hfx_bog, hfx_bom, hfx_dif, hfx_dyn, hfx_opw, hfx_res  )
      !$acc data present( hfx_snw, hfx_spr, hfx_sub, hfx_sum, hfx_thd, h_i, h_s, qemp_ice, qemp_oce, qevap_ice, qns, qns_oce, qns_tot, qsr, qsr_oce, qsr_tot, qt_atm_oi    )
      !$acc data present( qt_oce_ai, qtr_ice_bot, sfx, sfx_bog, sfx_bom, sfx_bri, sfx_dyn, sfx_lam, sfx_opw, sfx_res, sfx_sni, sfx_sub, sfx_sum )
      !$acc data present( snwice_mass_b, tn_ice, t_su, vt_i, vt_s, wfx_bog, wfx_bom, wfx_dyn, wfx_err_sub, wfx_ice, wfx_ice_sub, wfx_lam, wfx_opw, wfx_pnd, wfx_res        )
      !$acc data present( wfx_sni, wfx_snw, wfx_snw_dyn, wfx_snw_sni, wfx_snw_sub, wfx_snw_sum, wfx_sub, wfx_sum )

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_flx: update fluxes (mass, salt and heat) at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

# if defined _OPENACC
      IF( ln_cndflx )  CALL ctl_stop( 'STOP', 'ice_update_flx : adapt option `ln_cndflx` for GPU!')
# endif

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            ! Net heat flux on top of the ice-ocean (W.m-2)
            !----------------------------------------------
            !IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
            !   qt_atm_oi(ji,jj) = ( 1._wp - at_i_b(ji,jj) ) * ( qns_oce(ji,jj) + qsr_oce(ji,jj) ) + qemp_oce(ji,jj) + &
            !      &             SUM( a_i_b(ji,jj,1:jpl) * ( qcn_ice(ji,jj,1:jpl) + qml_ice(ji,jj,1:jpl) + qtr_ice_top(ji,jj,1:jpl) ), dim=3 ) + qemp_ice(ji,jj)
            !ELSE                   ! ice-atm interface = solar and non-solar fluxes
            qt_atm_oi(ji,jj) = qns_tot(ji,jj) + qsr_tot(ji,jj)
            !ENDIF

            ! --- case we bypass ice thermodynamics --- !
            IF( .NOT. ln_icethd ) THEN   ! we suppose ice is impermeable => ocean is isolated from atmosphere
               IF(lwp) PRINT *, 'LOLO: `ice_update_flx@iceupdate.F90`:  qt_oce_ai bypass thermo! ,  kt =', kt
               qt_atm_oi(ji,jj)   = ( 1._wp - at_i_b(ji,jj) ) * ( qns_oce(ji,jj) + qsr_oce(ji,jj) ) + qemp_oce(ji,jj)
               qt_oce_ai(ji,jj)   = ( 1._wp - at_i_b(ji,jj) ) *   qns_oce(ji,jj)                    + qemp_oce(ji,jj)
               emp_ice  (ji,jj)   = 0._wp
               qemp_ice (ji,jj)   = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  qevap_ice(ji,jj,jl) = 0._wp
               END DO
            ENDIF


            ! Solar heat flux reaching the ocean (max) = zqsr (W.m-2)
            !---------------------------------------------------
            !IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
            !   zqsr = ( 1._wp - at_i_b(ji,jj) ) * qsr_oce(ji,jj) + SUM( a_i_b (ji,jj,:) * qtr_ice_bot(ji,jj,:) )
            !ELSE                   ! ice-atm interface = solar and non-solar fluxes
            zsum = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zsum = zsum + a_i_b(ji,jj,jl) * ( qsr_ice(ji,jj,jl) - qtr_ice_bot(ji,jj,jl) )
            END DO
            zqsr = qsr_tot(ji,jj) - zsum
            !ENDIF

            ! Total heat flux reaching the ocean = qt_oce_ai (W.m-2)
            !---------------------------------------------------
            IF( ln_icethd ) THEN
               !
               zsum = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zsum = zsum + qevap_ice(ji,jj,jl) * a_i_b(ji,jj,jl)
               END DO

               qt_oce_ai(ji,jj) = qt_atm_oi(ji,jj) - hfx_sum(ji,jj) - hfx_bom(ji,jj) - hfx_bog(ji,jj) &
                  &                                - hfx_dif(ji,jj) - hfx_opw(ji,jj) - hfx_snw(ji,jj) &
                  &                                + hfx_thd(ji,jj) + hfx_dyn(ji,jj) + hfx_res(ji,jj) &
                  &                                + hfx_sub(ji,jj) - zsum + hfx_spr(ji,jj)
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
               qsr(ji,jj) = ( 1._wp - at_i_b(ji,jj) ) * qsr_oce(ji,jj) * ( 1._wp - frq_m(ji,jj) ) &
                                !                                   + solar flux transmitted thru ice and the 1st ocean level (also not used by sea-ice)
                  &             + zsum * ( 1._wp - frq_m(ji,jj) )
               !
            ELSE                                                       !-- cooling or no ice left
               qsr(ji,jj) = zqsr
            ENDIF
            !
            ! the non-solar is simply derived from the solar flux
            !IF(lwp .AND. (ji==10 .AND. jj==10)) PRINT *, 'LOLO: `ice_update_flx@iceupdate.F90`: qns = qt_oce_ai - qsr ,  kt =', kt
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
            fmmflx(ji,jj) = - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_err_sub(ji,jj) ! ice-ocean mass flux saved at least for biogeochemical model
            emp   (ji,jj) = emp_oce(ji,jj) + fmmflx(ji,jj)   ! atm-ocean + ice-ocean mass flux
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

            ! Storing the transmitted variables
            !----------------------------------
            fr_i(ji,jj)   = at_i(ji,jj)             ! Sea-ice fraction
            !$acc loop seq
            DO jl=1, jpl
               tn_ice(ji,jj,jl) = t_su(ji,jj,jl)           ! Ice surface temperature
            END DO

         END DO
      END DO
      !$acc end parallel loop

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
      !
      ! output all fluxes
      !------------------
      !
      ! --- salt fluxes [kg/m2/s] --- !
      !                           ! sfxice =  sfxbog + sfxbom + sfxsum + sfxsni + sfxopw + sfxres + sfxdyn + sfxbri + sfxsub + sfxlam
      IF( iom_use('sfxice'  ) )   CALL iom_put( 'sfxice', sfx     * 1.e-03 )   ! salt flux from total ice growth/melt
      IF( iom_use('sfxbog'  ) )   CALL iom_put( 'sfxbog', sfx_bog * 1.e-03 )   ! salt flux from bottom growth
      IF( iom_use('sfxbom'  ) )   CALL iom_put( 'sfxbom', sfx_bom * 1.e-03 )   ! salt flux from bottom melting
      IF( iom_use('sfxsum'  ) )   CALL iom_put( 'sfxsum', sfx_sum * 1.e-03 )   ! salt flux from surface melting
      IF( iom_use('sfxlam'  ) )   CALL iom_put( 'sfxlam', sfx_lam * 1.e-03 )   ! salt flux from lateral melting
      IF( iom_use('sfxsni'  ) )   CALL iom_put( 'sfxsni', sfx_sni * 1.e-03 )   ! salt flux from snow ice formation
      IF( iom_use('sfxopw'  ) )   CALL iom_put( 'sfxopw', sfx_opw * 1.e-03 )   ! salt flux from open water formation
      IF( iom_use('sfxdyn'  ) )   CALL iom_put( 'sfxdyn', sfx_dyn * 1.e-03 )   ! salt flux from ridging rafting
      IF( iom_use('sfxbri'  ) )   CALL iom_put( 'sfxbri', sfx_bri * 1.e-03 )   ! salt flux from brines
      IF( iom_use('sfxres'  ) )   CALL iom_put( 'sfxres', sfx_res * 1.e-03 )   ! salt flux from undiagnosed processes
      IF( iom_use('sfxsub'  ) )   CALL iom_put( 'sfxsub', sfx_sub * 1.e-03 )   ! salt flux from sublimation

      ! --- mass fluxes [kg/m2/s] --- !
      CALL iom_put( 'emp_oce', emp_oce )   ! emp over ocean (taking into account the snow blown away from the ice)
      CALL iom_put( 'emp_ice', emp_ice )   ! emp over ice   (taking into account the snow blown away from the ice)

      !                           ! vfxice = vfxbog + vfxbom + vfxsum + vfxsni + vfxopw + vfxdyn + vfxres + vfxlam + vfxpnd
      CALL iom_put( 'vfxice'    , wfx_ice     )   ! mass flux from total ice growth/melt
      CALL iom_put( 'vfxbog'    , wfx_bog     )   ! mass flux from bottom growth
      CALL iom_put( 'vfxbom'    , wfx_bom     )   ! mass flux from bottom melt
      CALL iom_put( 'vfxsum'    , wfx_sum     )   ! mass flux from surface melt
      CALL iom_put( 'vfxlam'    , wfx_lam     )   ! mass flux from lateral melt
      CALL iom_put( 'vfxsni'    , wfx_sni     )   ! mass flux from snow-ice formation
      CALL iom_put( 'vfxopw'    , wfx_opw     )   ! mass flux from growth in open water
      CALL iom_put( 'vfxdyn'    , wfx_dyn     )   ! mass flux from dynamics (ridging)
      CALL iom_put( 'vfxres'    , wfx_res     )   ! mass flux from undiagnosed processes
      CALL iom_put( 'vfxpnd'    , wfx_pnd     )   ! mass flux from melt ponds
      CALL iom_put( 'vfxsub'    , wfx_ice_sub )   ! mass flux from ice sublimation (ice-atm.)
      CALL iom_put( 'vfxsub_err', wfx_err_sub )   ! "excess" of sublimation sent to ocean

      !IF ( iom_use( 'vfxthin' ) ) THEN   ! mass flux from ice growth in open water + thin ice (<20cm) => comparable to observations
      !   WHERE( hm_i(:,:) < 0.2 .AND. hm_i(:,:) > 0. )
      !      z2d = wfx_bog
      !   ELSEWHERE
      !      z2d = 0._wp
      !   END WHERE
      !   CALL iom_put( 'vfxthin', wfx_opw + z2d )
      !ENDIF

      !                            ! vfxsnw = vfxsnw_sni + vfxsnw_dyn + vfxsnw_sum
      CALL iom_put( 'vfxsnw'     , wfx_snw     )   ! mass flux from total snow growth/melt
      CALL iom_put( 'vfxsnw_sum' , wfx_snw_sum )   ! mass flux from snow melt at the surface
      CALL iom_put( 'vfxsnw_sni' , wfx_snw_sni )   ! mass flux from snow melt during snow-ice formation
      CALL iom_put( 'vfxsnw_dyn' , wfx_snw_dyn )   ! mass flux from dynamics (ridging)
      CALL iom_put( 'vfxsnw_sub' , wfx_snw_sub )   ! mass flux from snow sublimation (ice-atm.)
      CALL iom_put( 'vfxsnw_pre' , wfx_spr     )   ! snow precip

      ! --- heat fluxes [W/m2] --- !
      !                              ! qt_atm_oi - qt_oce_ai = hfxdhc - ( dihctrp + dshctrp )

      !IF( iom_use('qsr_oce_si') .OR. iom_use('qns_oce_si') .OR. iom_use('qemp_oce_si') .OR. iom_use('qns_atmo') ) THEN
      !LB:
      !      !z2d(:,:) = xmskt(:,:)
      !
      !   !! LB: '_si' means we just keep regions where there is sea-ice (field will be 0 where A=0) !
      !   WHERE( at_i_b(:,:) <= 0.01_wp ) z2d(:,:) = 0._wp
      !   IF( iom_use('qemp_oce_si') ) CALL iom_put( 'qemp_oce_si', qemp_oce                      * z2d ) ! Downward Heat Flux from E-P over ocean
      !   IF( iom_use('qsr_oce_si' ) ) CALL iom_put( 'qsr_oce_si' , qsr_oce  * ( 1._wp - at_i_b ) * z2d ) !     solar flux at ocean surface
      !   IF( iom_use('qns_oce_si' ) ) CALL iom_put( 'qns_oce_si' , qns_oce  * ( 1._wp - at_i_b ) * z2d ) ! non-solar flux at ocean surface !LOLO: add `qemp_oce` ??? don't get it...
      !   IF( iom_use('qns_atmo'   ) ) CALL iom_put( 'qns_atmo'   , ( -SUM( qns_ice * a_i_b, dim=3 ) - qns_oce*( 1._wp - at_i_b ) ) * z2d ) ! Non solar heat flux to the atmosphere
      !   !!
      !   z2d(:,:) = xmskt(:,:)
      !ENDIF
      !LB.

      IF( iom_use('qsr_oce'    ) ) CALL iom_put( 'qsr_oce'    , (qsr_oce * ( 1._wp - at_i_b )                              ) * xmskt ) !     solar flux at ocean surface
      IF( iom_use('qns_oce'    ) ) CALL iom_put( 'qns_oce'    , (qns_oce * ( 1._wp - at_i_b ) + qemp_oce                   ) * xmskt ) ! non-solar flux at ocean surface
      !IF( iom_use('qsr_ice'    ) ) CALL iom_put( 'qsr_ice'    , (SUM( qsr_ice * a_i_b, dim=3 )                             ) * xmskt ) !     solar flux at ice surface ! --> moved to `blk_ice_2` of `sbcblk.F90`
      !#LOLOFixme: the true `qns_ice` is saved in `blk_ice_2` of `sbcblk.F90`, change the name for this one:
      !IF( iom_use('qns_ice'    ) ) CALL iom_put( 'qns_ice'    , (SUM( qns_ice * a_i_b, dim=3 ) + qemp_ice                  ) * xmskt ) ! non-solar flux at ice surface
      !#LOLOfixme.
      IF( iom_use('qtr_ice_bot') ) CALL iom_put( 'qtr_ice_bot', (SUM( qtr_ice_bot * a_i_b, dim=3 )                         ) * xmskt ) !     solar flux transmitted thru ice
      IF( iom_use('qtr_ice_top') ) CALL iom_put( 'qtr_ice_top', (SUM( qtr_ice_top * a_i_b, dim=3 )                         ) * xmskt ) !     solar flux transmitted thru ice surface
      IF( iom_use('qt_oce'     ) ) CALL iom_put( 'qt_oce'     , (     ( qsr_oce + qns_oce ) * ( 1._wp - at_i_b ) + qemp_oce) * xmskt )
      IF( iom_use('qt_ice'     ) ) CALL iom_put( 'qt_ice'     , (SUM( ( qns_ice + qsr_ice ) * a_i_b, dim=3 )     + qemp_ice) * xmskt )
      IF( iom_use('qt_oce_ai'  ) ) CALL iom_put( 'qt_oce_ai'  , qt_oce_ai                                                    * xmskt ) ! total heat flux at the ocean   surface: interface oce-(ice+atm)
      IF( iom_use('qt_atm_oi'  ) ) CALL iom_put( 'qt_atm_oi'  , qt_atm_oi                                                    * xmskt ) ! total heat flux at the oce-ice surface: interface atm-(ice+oce)
      IF( iom_use('qemp_oce'   ) ) CALL iom_put( 'qemp_oce'   , (qemp_oce                                                  ) * xmskt ) ! Downward Heat Flux from E-P over ocean
      IF( iom_use('qemp_ice'   ) ) CALL iom_put( 'qemp_ice'   , (qemp_ice                                                  ) * xmskt ) ! Downward Heat Flux from E-P over ice

      ! heat fluxes from ice transformations
      !                            ! hfxdhc = hfxbog + hfxbom + hfxsum + hfxopw + hfxdif + hfxsnw - ( hfxthd + hfxdyn + hfxres + hfxsub + hfxspr )
      CALL iom_put ('hfxbog'     , hfx_bog     )   ! heat flux used for ice bottom growth
      CALL iom_put ('hfxbom'     , hfx_bom     )   ! heat flux used for ice bottom melt
      CALL iom_put ('hfxsum'     , hfx_sum     )   ! heat flux used for ice surface melt
      CALL iom_put ('hfxopw'     , hfx_opw     )   ! heat flux used for ice formation in open water
      CALL iom_put ('hfxdif'     , hfx_dif     )   ! heat flux used for ice temperature change
      CALL iom_put ('hfxsnw'     , hfx_snw     )   ! heat flux used for snow melt
      CALL iom_put ('hfxerr'     , hfx_err_dif )   ! heat flux error after heat diffusion

      ! heat fluxes associated with mass exchange (freeze/melt/precip...)
      CALL iom_put ('hfxthd'     , hfx_thd     )   !
      CALL iom_put ('hfxdyn'     , hfx_dyn     )   !
      CALL iom_put ('hfxres'     , hfx_res     )   !
      CALL iom_put ('hfxsub'     , hfx_sub     )   !
      CALL iom_put ('hfxspr'     , hfx_spr     )   ! Heat flux from snow precip heat content

      ! other heat fluxes
      IF( iom_use('hfxsensib'  ) )   CALL iom_put( 'hfxsensib'  ,      qsb_ice_bot * at_i_b         )   ! Sensible oceanic heat flux
      IF( iom_use('hfxcndbot'  ) )   CALL iom_put( 'hfxcndbot'  , SUM( qcn_ice_bot * a_i_b, dim=3 ) )   ! Bottom conduction flux
      IF( iom_use('hfxcndtop'  ) )   CALL iom_put( 'hfxcndtop'  , SUM( qcn_ice_top * a_i_b, dim=3 ) )   ! Surface conduction flux
      IF( iom_use('hfxmelt'    ) )   CALL iom_put( 'hfxmelt'    , SUM( qml_ice     * a_i_b, dim=3 ) )   ! Surface melt flux
      IF( iom_use('hfxldmelt'  ) )   CALL iom_put( 'hfxldmelt'  ,      fhld        * at_i_b         )   ! Heat in lead for ice melting
      IF( iom_use('hfxldgrow'  ) )   CALL iom_put( 'hfxldgrow'  ,      qlead       * r1_Dt_ice      )   ! Heat in lead for ice growth

      ! controls
      !---------
      IF( ln_icediachk      )   CALL ice_cons_final('ice_update_flx')                                       ! conservation
      IF( ln_icectl         )   CALL ice_prt       (kt, iiceprt, jiceprt, 3, 'Final state ice_update') ! prints
      IF( sn_cfctl%l_prtctl )   CALL ice_prt3D     ('ice_update_flx')                                       ! prints


      !$acc end data
      !$acc end data
      !$acc end data
      !$acc end data
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
      !! ** Action  : * at each ice time step (every nn_fsbc time step):
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
      !! ** Outputs : - utau, vtau   : surface ocean i- and j-stress (u- & v-pts) updated with ice-ocean fluxes
      !!              - taum         : modulus of the surface ocean stress (T-point) updated with ice-ocean fluxes
      !!---------------------------------------------------------------------
      INTEGER ,                     INTENT(in) ::   kt               ! ocean time-step index
      !
      INTEGER  ::   ji, jj, nb   ! dummy loop indices
      REAL(wp) ::   za_tot_u, ztaux_ai_u, zu_t, zmodt   ! local scalar
      REAL(wp) ::   za_tot_v, ztauy_ai_v, zv_t, zrhoco  !   -      -
      REAL(wp) ::   ztaux_oi_u, ztauy_oi_v, ztmod_io
      REAL(wp) ::   zflagi, zA, zB                  !   -      -
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_update_tau')
      !$acc data present( utau, vtau, u_ice, v_ice, uVice, vUice, V_oce, taum, at_i, au_i, av_i )

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_tau: update stress at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

      zrhoco = rho0 * rn_Cd_io

      zflagi = MERGE( 0._wp   ,   1._wp  ,   ln_drgice_imp )

      nb = MAX( nn_hls-1, 0 )
      
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nb, Nje0+nb
         DO ji=Nis0-nb, Nie0+nb
            !                                               ! 2*(U_ice-U_oce) at T-point
            zu_t = u_ice(ji,jj) + u_ice(ji-1,jj) - V_oce(ji,jj,1) - V_oce(ji-1,jj,1)
            zv_t = v_ice(ji,jj) + v_ice(ji,jj-1) - V_oce(ji,jj,2) - V_oce(ji,jj-1,2)
            !                                              ! |U_ice-U_oce|^2
            zmodt =  0.25_wp * (  zu_t * zu_t + zv_t * zv_t  )
            !                                               ! update the ocean stress modulus
            taum(ji,jj) = ( 1._wp - at_i(ji,jj) ) * taum(ji,jj) + at_i(ji,jj) * zrhoco * zmodt
            ztmod_io = zrhoco * SQRT( zmodt )          ! rhoco * |U_ice-U_oce| at T-point
            !
            !
            ztaux_oi_u = utau(ji,jj)
            ztauy_oi_v = vtau(ji,jj)

            !#LOLOreview!
            ! ice area at u and v-points
            !za_tot_u  = ( at_i(ji,jj) * xmskt(ji,jj) + at_i (ji+1,jj    ) * xmskt(ji+1,jj  ) )  &
            !   &     / MAX( 1.0_wp , xmskt(ji,jj) + xmskt(ji+1,jj  ) )
            !za_tot_v  = ( at_i(ji,jj) * xmskt(ji,jj) + at_i (ji  ,jj+1  ) * xmskt(ji  ,jj+1) )  &
            !   &     / MAX( 1.0_wp , xmskt(ji,jj) + xmskt(ji  ,jj+1) )
            za_tot_u = au_i(ji,jj)
            za_tot_v = av_i(ji,jj)
            !                                                   ! linearized quadratic drag formulation
            !ztaux_ai_u   = 0.5_wp * ( ztmod_io + tmod_io(ji+1,jj) ) * ( u_ice(ji,jj) - zflagi * pu_oce(ji,jj) )
            !ztauy_ai_v   = 0.5_wp * ( ztmod_io + tmod_io(ji,jj+1) ) * ( v_ice(ji,jj) - zflagi * pv_oce(ji,jj) )
            !
            zA = u_ice(ji,jj) - V_oce(ji,jj,1) ! u_ice - u_oce @ U
            zB = vUice(ji,jj) - V_oce(ji,jj,4) ! v_ice - v_oce @ U
            ztmod_io = SQRT( zA*zA + zB*zB )   ! modulus of ice - oce @ U 
            ztaux_ai_u = zrhoco * ztmod_io * ( u_ice(ji,jj) - zflagi * V_oce(ji,jj,1) )

            zA = uVice(ji,jj) - V_oce(ji,jj,3) ! u_ice - u_oce @ V
            zB = v_ice(ji,jj) - V_oce(ji,jj,2) ! v_ice - v_oce @ V
            ztmod_io = SQRT( zA*zA + zB*zB )   ! modulus of ice - oce @ U 
            ztauy_ai_v = zrhoco * ztmod_io * ( v_ice(ji,jj) - zflagi * V_oce(ji,jj,2) )
            
            !#LOLOreview.
            !                                                   ! stresses at the ocean surface
            utau(ji,jj) = ( 1._wp - za_tot_u ) * ztaux_oi_u + za_tot_u * ztaux_ai_u
            vtau(ji,jj) = ( 1._wp - za_tot_v ) * ztauy_oi_v + za_tot_v * ztauy_ai_v
            !
         END DO
      END DO
      !$acc end parallel loop

# if ! defined _OPENACC
      CALL lbc_lnk( 'ice_update_tau', taum,'T',1._wp, utau,'U',-1._wp, vtau,'V',-1._wp )   ! lateral boundary condition
# endif

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
            id1 = iom_varid( numrir, 'snwice_mass' , ldstop = .FALSE. )
            !
            IF( id1 > 0 ) THEN                       ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'snwice_mass'  , snwice_mass   )
               CALL iom_get( numrir, jpdom_auto, 'snwice_mass_b', snwice_mass_b )
            ELSE                                     ! start from rest
               IF(lwp) WRITE(numout,*) '   ==>>   previous run without snow-ice mass output then set it'
               snwice_mass  (:,:) = xmskt(:,:) * ( rhos * vt_s(:,:) + rhoi * vt_i(:,:) &
                  &  + rhow * (vt_ip(:,:) + vt_il(:,:))  )
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
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
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
