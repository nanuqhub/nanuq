MODULE icewri
   !!======================================================================
   !!                     ***  MODULE  icewri  ***
   !!   sea-ice : output ice variables
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_wri       : write of the diagnostics variables in ouput file
   !!   ice_wri_state : write for initial state or/and abandon
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! domain: ocean
   USE sbc_oce        ! surf. boundary cond.: ocean
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE par_ice
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   USE icealb  , ONLY : rn_alb_oce
   USE oss_nnq , ONLY : sst_m, e3t_m, frq_m
   !
   USE ioipsl         !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   USE remap_classic, ONLY : rmpT2F, rmpF2T, rmpT2U, rmpT2V, do_rmpF2T
   USE remap_weno,    ONLY : rmpT2F_wn5s !rmpT2U_wn5s, rmpT2V_wn5s

   USE icedyn_rhg_tools, ONLY : strain_rate_all, sigmaII_full, P_max_diag, P_tilde_diag, Elast_diag, Visco_diag, Lambda_diag, vel_div_t, vel_shear_f, vel_maxshr_t

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_wri        ! called by ice_stp
   PUBLIC ice_wri_adv    ! called by ice_stp (minimal set for advection-only tests)
   PUBLIC ice_wri_state  ! called by dia_wri_state
   

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icewri.F90 15388 2021-10-17 11:33:47Z clem $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_wri( kt )
      !!-------------------------------------------------------------------
      !!  This routine ouputs some (most?) of the sea ice fields
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time-step
      !
      INTEGER  ::   ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) ::   z2da, z2db, zrho1, zrho2
      REAL(wp) ::   zmiss_val       ! missing value retrieved from xios
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d                            ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zmskt00l, zmsktsnl               ! cat masks
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::   zfast, zalb, zmsktalb      ! 2D workspace
      REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::   z1_h_t, z1_h_f
      REAL(wp), DIMENSION(jpi,jpj) :: zimskt, zimskf
      REAL(wp), DIMENSION(jpi,jpj) :: ztmp1, ztmp2, ztmp3, ztmp4
      REAL(wp)                     :: zt1, zt2
      !
      ! Global ice diagnostics (SIMIP)
      REAL(wp) ::   zdiag_area_nh, zdiag_extt_nh, zdiag_volu_nh   ! area, extent, volume
      REAL(wp) ::   zdiag_area_sh, zdiag_extt_sh, zdiag_volu_sh
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_wri')
      !$acc data create( z2d )


      ! Ice-extent masks at T- & F-points:
      !$acc update self( kmsk_ice_t, kmsk_ice_f )
      zimskt(:,:) = REAL(kmsk_ice_t,wp)
      zimskf(:,:) = REAL(kmsk_ice_f,wp)

      !-----------------
      ! Standard outputs
      !-----------------
      zrho1 = ( rho0 - rhoi ) * r1_rho0 ; zrho2 = rhos * r1_rho0
      ! masks
      CALL iom_put( 'icemask'  , zimskt )   ! ice mask 0%
      !
      ! general fields
      
      IF( iom_use('ff_u' ) )   CALL iom_put( 'ff_u', ff_u )
      IF( iom_use('ff_v' ) )   CALL iom_put( 'ff_v', ff_v )

      
      IF( iom_use('icemass' ) )   CALL iom_put( 'icemass', vt_i * rhoi * zimskt )       ! Ice mass per cell area
      IF( iom_use('snwmass' ) )   CALL iom_put( 'snwmass', vt_s * rhos * zimskt )       ! Snow mass per cell area
      IF( iom_use('iceconc' ) ) THEN
         !$acc update self( at_i )
         CALL iom_put( 'iceconc',   at_i ) ! * zimskt )   ! ice concentration
      ENDIF
      IF( iom_use('iceconc-f') )  THEN
         !$acc update self( af_i )
         CALL iom_put( 'iceconc-f', af_i ) ! * zimskf )   ! ice concentration
      ENDIF
      IF( iom_use('icevolu' ) ) THEN
         CALL iom_put( 'icevolu', vt_i * zimskt )       ! ice volume = mean ice thickness over the cell
      ENDIF
      IF( iom_use('icethic' ) ) THEN
         !$acc update self( hm_i )
         CALL iom_put( 'icethic', hm_i * zimskt )   ! ice thickness
      ENDIF
      IF( iom_use('icethic-f' ) ) THEN
         !$acc update self( hm_i_f )
         CALL iom_put( 'icethic-f', hm_i_f * zimskf )   ! ice thickness
      ENDIF
      
      IF( iom_use('snwthic' ) ) THEN
         !$acc update self( hm_s )
         CALL iom_put( 'snwthic', hm_s      * zimskt )       ! snw thickness
      ENDIF
      IF( iom_use('snwvolu' ) ) THEN
         !$acc update self( vt_s )
         CALL iom_put( 'snwvolu', vt_s        * zimskt )       ! snow volume
      ENDIF
      IF( ln_damage ) THEN
         IF( iom_use('icedmgt') ) THEN
            !$acc update self( dmdt )
            CALL iom_put( 'icedmgt', (1._wp - dmdt) * zimskt )          ! ice damage @T !
         ENDIF
         IF( iom_use('icedmgf') ) THEN
            !$acc update self( dmdf )
            CALL iom_put( 'icedmgf', (1._wp - dmdf) * zimskf )          ! ice damage @T !
         ENDIF
      END IF


      IF( iom_use('icevolu-u') ) THEN
         z2d(:,:) = MAX( rmpT2U( vt_i, lbcl=.TRUE., lconserv=.TRUE. ) , 0._wp )
         CALL iom_put( 'icevolu-u', z2d )
      ENDIF
      IF( iom_use('icevolu-v') ) THEN
         z2d(:,:) = MAX( rmpT2V( vt_i, lbcl=.TRUE., lconserv=.TRUE. ) , 0._wp )
         CALL iom_put( 'icevolu-v', z2d )
      ENDIF
      IF( iom_use('icevolu-f') ) THEN
         IF(ln_use_weno_rmp) THEN
            z2d(:,:) = MAX( rmpT2F_wn5s( vt_i, lbcl=.TRUE. )                  , 0._wp )
         ELSE
            z2d(:,:) = MAX( rmpT2F(      vt_i, lbcl=.TRUE., lconserv=.TRUE. ) , 0._wp )
         ENDIF
         CALL iom_put( 'icevolu-f', z2d * zimskf )
      ENDIF

      IF( iom_use('iceconc-u') ) CALL iom_put( 'iceconc-u', au_i )
      IF( iom_use('iceconc-v') ) CALL iom_put( 'iceconc-v', av_i )


      ! momentum
      IF( iom_use('uice'    ) ) CALL iom_put( 'uice'   , u_ice ) !* zimskt   )                     ! ice velocity u
      IF( iom_use('vice'    ) ) CALL iom_put( 'vice'   , v_ice ) !* zimskt   )                     ! ice velocity v

      IF( iom_use('uice_v'  ) ) THEN
         !$acc update self( uVice )
         CALL iom_put( 'uice_v' , uVice   )
      ENDIF
      IF( iom_use('vice_u'  ) ) THEN
         !$acc update self( vUice )
         CALL iom_put( 'vice_u' , vUice   )
      ENDIF
      !
      IF( iom_use('icevel') .OR. iom_use('fasticepres') ) THEN                                                              ! module of ice velocity & fast ice
         ALLOCATE( zfast(jpi,jpj) )
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               z2da  = u_ice(ji,jj) + u_ice(ji-1,jj)
               z2db  = v_ice(ji,jj) + v_ice(ji,jj-1)
               z2d(ji,jj) = 0.5_wp * SQRT( z2da * z2da + z2db * z2db )
            END DO
         END DO
         CALL lbc_lnk( 'icewri', z2d, 'T', 1.0_wp )
         CALL iom_put( 'icevel', z2d )

         WHERE( z2d(:,:) < 5.e-04_wp .AND. zimskt(:,:) == 1._wp )
            zfast(:,:) = 1._wp                                      ! record presence of fast ice
         ELSEWHERE
            zfast(:,:) = 0._wp
         END WHERE
         CALL iom_put( 'fasticepres', zfast )
         DEALLOCATE( zfast )
      ENDIF

      IF( ln_damage ) THEN
         IF( iom_use('icevelf') ) THEN                                                              ! module of ice velocity @ F points
            DO jj=Njs0, Nje0
               DO ji=Nis0, Nie0
                  z2da  = uVice(ji+1,jj) + uVice(ji,jj)
                  z2db  = vUice(ji,jj+1) + vUice(ji,jj)
                  z2d(ji,jj) = 0.5_wp * SQRT( z2da * z2da + z2db * z2db )
               END DO
            END DO
            CALL lbc_lnk( 'icewri',  z2d, 'F', 1.0_wp )
            CALL iom_put( 'icevelf', z2d * zimskf )
         ENDIF
      ENDIF




      !! ***********
      !!  Rheology
      !! ***********


      !! ~~~~~~~~~~~~~~~~~~~~
      !!  Internal stresses
      !! ~~~~~~~~~~~~~~~~~~~~      
      IF( iom_use('sheastr') .OR. iom_use('ice_sig11f') .OR. iom_use('ice_sig22f') .OR. iom_use('ice_sig12t') .OR. &
         & iom_use('normstrf') .OR. iom_use('ice_sig11') .OR. iom_use('ice_sig22') .OR. iom_use('ice_sig12') .OR. iom_use('normstr') ) THEN

         ALLOCATE( z1_h_t(jpi,jpj), z1_h_f(jpi,jpj) )

         ! Since internal stresses are living in the code as vertically-integrated (`sigma*h` in [N/m^2*m] rather than sigma` in [N/m^2]
         !! => we need the multiplicator to fall back on actuall stresses for the output files
         IF( l_use_v_for_h ) THEN
            z1_h_t(:,:) = MAX( vt_i(:,:) , 0._wp)
            z1_h_f(:,:) = MAX( rmpT2F( z1_h_t,  lconserv=.TRUE. ) , 0._wp )
            CALL lbc_lnk( 'icedyn_rhg_bbm', z1_h_f,'F',1._wp )
         ELSE
            !$acc update self( hm_i, hm_i_f )
            z1_h_t(:,:) = hm_i(:,:)
            z1_h_f(:,:) = hm_i_f(:,:)
         ENDIF
         z1_h_t(:,:) = zimskt(:,:) / MAX( z1_h_t(:,:), epsi20 )
         z1_h_f(:,:) = zimskf(:,:) / MAX( z1_h_f(:,:), epsi20 )


         IF( iom_use('sheastr') ) THEN
            IF( ln_damage ) THEN
               !! We don't need to interpolate S12 from the corner to center !
               CALL sigmaII_full( SIGMAt, SIGMAf, z2d )
            ELSE
               !! We need to interpolate....
               CALL do_rmpF2T( SIGMAt(:,:,3), z2d  ) !LOLOfixme: ACC `present` might not recognize `SIGMAt(:,:,3)` !
               !$acc parallel loop collapse(2)
               DO jj=Njs0-1, Nje0+1
                  DO ji=Nis0-1, Nie0+1
                     zt1         = 0.5_wp * ( SIGMAt(ji,jj,1) - SIGMAt(ji,jj,2) )
                     zt2         = z2d(ji,jj)
                     z2d(ji,jj) = SQRT( zt1*zt1 + zt2*zt2 )
                  END DO
               END DO
               !$acc end parallel loop
            ENDIF
            !$acc update self( z2d )
            CALL iom_put( 'sheastr',  z2d(:,:)*z1_h_t(:,:) )
         ENDIF
         !! Same at F-points when E-grid is used:
         !IF( ln_damage ) THEN
         !IF( iom_use('sheastrf')) CALL iom_put( 'sheastrf', sigmaII( SIGMAf(:,:,1), SIGMAf(:,:,2), SIGMAt(:,:,3) )*z1_h_f(:,:) )
         !ENDIF

         !! Saving stress tensor components in the proper units! i.e. Pa :
         IF( iom_use('ice_sig11') .OR. iom_use('ice_sig22') .OR. iom_use('ice_sig12') .OR. iom_use('normstr') ) THEN
            !$acc update self( SIGMAt )
            IF( iom_use('ice_sig11')  ) CALL iom_put( 'ice_sig11' ,  SIGMAt(:,:,1)*z1_h_t(:,:) )
            IF( iom_use('ice_sig22')  ) CALL iom_put( 'ice_sig22' ,  SIGMAt(:,:,2)*z1_h_t(:,:) )
            IF( iom_use('ice_sig12')  ) CALL iom_put( 'ice_sig12' ,  SIGMAt(:,:,3)*z1_h_f(:,:) )
            IF( iom_use('normstr')    ) CALL iom_put( 'normstr', 0.5_wp*(SIGMAt(:,:,1)+SIGMAt(:,:,2))*z1_h_t(:,:) ) ! Sigma_I
         ENDIF

         IF( ln_damage ) THEN
            !! --> following fields only exist on the E-grid...

            IF( iom_use('ice_sig11f') .OR. iom_use('ice_sig22f') .OR. iom_use('ice_sig12t') .OR. iom_use('normstrf') ) THEN
               !$acc update self( SIGMAf )
               IF( iom_use('ice_sig11f') ) CALL iom_put( 'ice_sig11f',  SIGMAf(:,:,1)*z1_h_f(:,:) )
               IF( iom_use('ice_sig22f') ) CALL iom_put( 'ice_sig22f',  SIGMAf(:,:,2)*z1_h_f(:,:) )
               IF( iom_use('ice_sig12t') ) CALL iom_put( 'ice_sig12t',  SIGMAf(:,:,3)*z1_h_t(:,:) )
               IF( iom_use('normstrf')   ) CALL iom_put( 'normstrf', 0.5_wp*(SIGMAf(:,:,1)+SIGMAf(:,:,2))*z1_h_f(:,:)) ! Sigma_I
            ENDIF

         ENDIF !IF( ln_damage )
         
         DEALLOCATE( z1_h_t, z1_h_f )
         
      END IF !IF( iom_use('sheastr') .OR. iom_use('ice_sig11f') .OR. iom_use('ice_sig22f') .OR. iom_use('ice_sig12t') ....


      !IF( iom_use('pmaxt') ) THEN
      !   z2d(:,:) = P_max_diag(   at_i, zht, SIGMAt(:,:,1), SIGMAt(:,:,2), SIGMAf(:,:,3) ) !ok to be vert.-int.
      !   CALL iom_put( 'pmaxt',  z2d(:,:)*zimskt )
      !ENDIF
      !IF( iom_use('ptildet') ) THEN
      !   z2d(:,:) = P_tilde_diag( at_i, zht, SIGMAt(:,:,1), SIGMAt(:,:,2), SIGMAf(:,:,3) ) !ok to be vert.-int.
      !   CALL iom_put( 'ptildet',  z2d(:,:)*zimskt )
      !ENDIF
      !
      !IF( iom_use('pmaxf') ) THEN
      !   z2d(:,:) = P_max_diag(   af_i, zhf, SIGMAf(:,:,1), SIGMAf(:,:,2), SIGMAt(:,:,3) ) !ok to be vert.-int.
      !   CALL iom_put( 'pmaxf',  z2d(:,:)*zimskf )
      !ENDIF
      !IF( iom_use('ptildef') ) THEN
      !   z2d(:,:) = P_tilde_diag( af_i, zhf, SIGMAf(:,:,1), SIGMAf(:,:,2), SIGMAt(:,:,3) ) !ok to be vert.-int.
      !   CALL iom_put( 'ptildef',  z2d(:,:)*zimskf )
      !ENDIF
      !
      !IF( iom_use('zelat') ) THEN
      !   z2d(:,:) = Elast_diag( at_i, dmdt )
      !   CALL iom_put( 'zelat',  z2d(:,:)*zimskt )
      !ENDIF
      !
      !IF( iom_use('zetat') ) THEN
      !   z2d(:,:) = Visco_diag( at_i, dmdt )
      !   CALL iom_put( 'zetat',  z2d(:,:)*zimskt )
      !ENDIF
      !
      !IF( iom_use('zlambt') ) THEN
      !   z2d(:,:) = Lambda_diag( at_i, dmdt, rdt_ice/REAL(nn_nbbm,wp) )
      !   CALL iom_put( 'zlambt',  z2d(:,:)*zimskt )
      !ENDIF
      !



      !! ~~~~~~~~~~~~~~~~~~~~
      !! Strain rate and co |
      !! ~~~~~~~~~~~~~~~~~~~~

      !! Strain rate of ice velocity stuff:
      IF( iom_use('icedivt') .OR. iom_use('iceshet') .OR. iom_use('icedeft') .OR. iom_use('icedelt') ) THEN
         !$acc update self( divu_i, shear_i, delta_i )
         ! --- divergence of velocity field @T:
         IF( iom_use('icedivt') )  CALL iom_put( 'icedivt' ,  divu_i*zimskt )
         ! --- MAXIMUM shear of velocity field @T:
         IF( iom_use('iceshet') )  CALL iom_put( 'iceshet' , shear_i*zimskt )
         !! --- total deformation of velocity field @T:
         IF( iom_use('icedeft') )  CALL iom_put( 'icedeft', SQRT( shear_i*shear_i + divu_i*divu_i )*zimskt )
         !! --- Delta @T:
         IF( iom_use('icedelt') )  CALL iom_put( 'icedelt' , delta_i*zimskt )
      ENDIF


      !!  ==> output of those that are computed by rheology:
      IF( ln_damage ) THEN

         IF( iom_use('e11t').OR.iom_use('e22t').OR.iom_use('e12t').OR.iom_use('iceshrt') ) THEN
            CALL strain_rate_all( 'T', u_ice, v_ice, uVice, vUice, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e1t2, e2t2, xmskt, &
               &                       pe11=ztmp1, pe22=ztmp2, pe12=ztmp3 )
            !LOLOfixme: should lbc_lnk these fields !!!
            IF( iom_use('e11t') ) CALL iom_put( 'e11t' , ztmp1*zimskt )
            IF( iom_use('e22t') ) CALL iom_put( 'e22t' , ztmp2*zimskt )
            IF( iom_use('e12t') ) CALL iom_put( 'e12t' , ztmp3*zimskt )
            ! --- shear of velocity field @T:
            IF( iom_use('iceshrt') )  CALL iom_put( 'iceshrt' , 2._wp*ztmp3*zimskt )
         ENDIF

         IF( iom_use('e11f').OR.iom_use('e22f').OR.iom_use('e12f').OR.iom_use('icedivf').OR.iom_use('iceshrf').OR.iom_use('iceshef').OR.iom_use('icedeff') ) THEN
            CALL strain_rate_all( 'F', uVice, vUice, u_ice, v_ice, r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, e1f2, e2f2, fmask(:,:,1), &
               &              pe11=ztmp1, pe22=ztmp2, pe12=ztmp3, pdiv=ztmp4, pmshr=z2d )
            !LOLOfixme: should lbc_lnk these fields !!!
            IF( iom_use('e11f') ) CALL iom_put( 'e11f' , ztmp1*zimskf )
            IF( iom_use('e22f') ) CALL iom_put( 'e22f' , ztmp2*zimskf )
            IF( iom_use('e12f') ) CALL iom_put( 'e12f' , ztmp3*zimskf )
            ! --- divergence of velocity field @F:
            IF( iom_use('icedivf') ) CALL iom_put( 'icedivf' , ztmp4*zimskf )
            ! --- shear of velocity field @F:
            IF( iom_use('iceshrf') ) CALL iom_put( 'iceshrf' , 2._wp*ztmp3*zimskf )
            ! --- MAXIMUM shear of velocity field @F:
            IF( iom_use('iceshef') ) CALL iom_put( 'iceshef' , z2d*zimskf )
            ! --- total deformation of velocity field @F:
            IF( iom_use('icedeff') ) CALL iom_put( 'icedeff', SQRT( z2d*z2d + ztmp4*ztmp4 )*zimskf )
         ENDIF

      ELSE
         !IF( iom_use('icedivt') .OR. iom_use('icedeft') ) THEN
         !   CALL vel_div_t( u_ice, v_ice, r1_e1e2t, e2u, e1v, xmskt, ztmp4 )
         !   CALL iom_put( 'icedivt' , ztmp4*zimskt )
         !ENDIF
         IF( iom_use('iceshrf') ) THEN
            CALL vel_shear_f( u_ice, v_ice, r1_e1e2f, r1_e1u, r1_e2v, e1f2, e2f2, fmask(:,:,1), ztmp3 )
            CALL iom_put( 'iceshrf' , ztmp3*zimskf )
         ENDIF
         !IF( iom_use('iceshet') .OR. iom_use('icedeft') ) THEN
         !   CALL vel_maxshr_t( u_ice, v_ice, r1_e1e2t, r1_e1e2f, r1_e1u, r1_e2v, r1_e2u, r1_e1v, e1t2, e2t2, e1f2, e2f2, &
         !      &                     e1e2f, xmskt, fmask(:,:,1), z2d )
         !   CALL iom_put( 'iceshet' , z2d*zimskt )
         !ENDIF
         !IF( iom_use('icedeft') ) CALL iom_put( 'icedeft', SQRT( z2d*z2d + ztmp4*ztmp4 )*zimskt )

      ENDIF !IF( ln_damage )


      !------------------------------------------------------------------------------!
      ! 5) diagnostics
      !------------------------------------------------------------------------------!

      ! --- vorticity of velocity field @F:
      IF( iom_use('icevorf') ) THEN
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               ztmp3(ji,jj) = (   ( v_ice(ji+1,jj)*r1_e2v(ji+1,jj) - v_ice(ji,jj)*r1_e2v(ji,jj) ) * e2f2(ji,jj) &
                  &             - ( u_ice(ji,jj+1)*r1_e1u(ji,jj+1) - u_ice(ji,jj)*r1_e1u(ji,jj) ) * e1f2(ji,jj) &
                  &           ) * r1_e1e2f(ji,jj) * fmask(ji,jj,1)  !#fixme: sure about `fmask` here ?
            END DO
         ENDDO
         CALL iom_put( 'icevorf' , ztmp3*zimskf )
      ENDIF
      ! --- vorticity of velocity field @T:
      IF( iom_use('icevort') ) THEN
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0
               ztmp3(ji,jj) = (   ( vUice(ji,jj)*r1_e2u(ji,jj) - vUice(ji-1,jj)*r1_e2u(ji-1,jj) ) * e2t2(ji,jj) &
                  &             - ( uVice(ji,jj)*r1_e1v(ji,jj) - uVice(ji,jj-1)*r1_e1v(ji,jj-1) ) * e1t2(ji,jj) &
                  &           ) * r1_e1e2t(ji,jj) * xmskt(ji,jj)
            END DO
         ENDDO
         CALL iom_put( 'icevort' , ztmp3*zimskt )
      ENDIF


      ! --- ice-atm. stress: No `A` involved !!!
      IF( iom_use('taux_ai')   ) CALL iom_put( 'taux_ai'   , utau_ice * zimskt )
      IF( iom_use('tauy_ai')   ) CALL iom_put( 'tauy_ai'   , vtau_ice * zimskt )

      IF( iom_use('taum_ai') ) THEN
         CALL iom_put( 'taum_ai' , SQRT(utau_ice*utau_ice + vtau_ice*vtau_ice)* zimskt )
      END IF

      ! --- oce-ice stress under sea-ice (with sign as seen from ice): No `A` involved !!!
      IF( iom_use('taux_oi_u') .OR. iom_use('taum_oi') ) THEN
         IF( ln_damage ) THEN
            ztmp1(:,:) = V_oce(:,:,1) - u_ice(:,:)  ! dU @U
            ztmp2(:,:) = V_oce(:,:,4) - vUice(:,:)  ! dV @U
            ztmp3(:,:) = rho0*rn_Cd_io * SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 ) * ztmp1(:,:) * REAL(kmsk_ice_u(:,:),wp)
            IF( iom_use('taux_oi_u') ) CALL iom_put( 'taux_oi_u' , ztmp3 )
         ELSE
            CALL ctl_stop( 'STOP', 'ice_wri: FIXME! Add `taux_oi_u` for EVP' )
         ENDIF
      ENDIF
      IF( iom_use('tauy_oi_v') .OR. iom_use('taum_oi') ) THEN
         IF( ln_damage ) THEN
            ztmp1(:,:) = V_oce(:,:,2) - v_ice(:,:)  ! dV @V
            ztmp2(:,:) = V_oce(:,:,3) - uVice(:,:)  ! dU @V
            ztmp4(:,:) = rho0*rn_Cd_io * SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 ) * ztmp1(:,:) * REAL(kmsk_ice_v(:,:),wp)
            IF( iom_use('tauy_oi_v') ) CALL iom_put( 'tauy_oi_v' , ztmp4 )
         ELSE
            CALL ctl_stop( 'STOP', 'ice_wri: FIXME! Add `tauy_oi_v` for EVP' )
         ENDIF
      ENDIF
      IF( iom_use('taum_oi') ) THEN
         ztmp1(2:jpi,:) = 0.5_wp * ( ztmp3(2:jpi,:) + ztmp3(1:jpi-1,:) )
         ztmp2(:,2:jpj) = 0.5_wp * ( ztmp4(:,2:jpj) + ztmp4(:,1:jpj-1) )
         CALL iom_put( 'taum_oi' , SQRT(ztmp1*ztmp1 + ztmp2*ztmp2) )
      END IF



      IF( ln_icethd ) THEN
         !-------------------------------------------------------------------------------------------------

         ! tresholds for outputs
         ALLOCATE( zmskt00l(jpi,jpj,jpl), zmsktsnl(jpi,jpj,jpl) )
         !
         DO jl = 1, jpl
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  zmskt00l(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , a_i(ji,jj,jl) - epsi06 ) )
                  zmsktsnl(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , v_s(ji,jj,jl) - epsi06 ) )
               END DO
            END DO
         END DO

         ! get missing value from xml
         ! CAUSES CRASHES, because `zmiss_val=1.E+020
         !CALL iom_miss_val( 'icetemp', zmiss_val )
         zmiss_val = -9999._wp


         ! brine volume
         CALL ice_var_brine

         IF( iom_use('icebrv'  ) ) THEN
            !$acc update self( vm_ibr )
            CALL iom_put( 'icebrv' , vm_ibr* 100. * zimskt )       ! brine volume
         ENDIF
         
         IF( iom_use('iceage') ) THEN
            !$acc update self( om_i )
            CALL iom_put( 'iceage' , om_i / rday * zimskt )
            !CALL iom_put( 'iceage' , om_i / rday * zimskt + zmiss_val * ( 1._wp - zimskt ) )          ! ice age            
         ENDIF

         IF( iom_use('icehnew' ) )   CALL iom_put( 'icehnew', ht_i_new             )                                           ! new ice thickness formed in the leads
         IF( iom_use('icefrb'  ) ) THEN                                                                                        ! Ice freeboard
            z2d(:,:) = ( zrho1 * hm_i(:,:) - zrho2 * hm_s(:,:) )
            WHERE( z2d < 0._wp )   z2d = 0._wp
            CALL iom_put( 'icefrb' , z2d * zimskt         )
         ENDIF
         ! melt ponds
         IF( iom_use('iceapnd' ) )   CALL iom_put( 'iceapnd', at_ip  * zimskt      )                                           ! melt pond total fraction
         IF( iom_use('icehpnd' ) )   CALL iom_put( 'icehpnd', hm_ip  * zimskt      )                                           ! melt pond depth
         IF( iom_use('icevpnd' ) )   CALL iom_put( 'icevpnd', vt_ip  * zimskt      )                                           ! melt pond total volume per unit area
         IF( iom_use('icehlid' ) )   CALL iom_put( 'icehlid', hm_il  * zimskt      )                                           ! melt pond lid depth
         IF( iom_use('icevlid' ) )   CALL iom_put( 'icevlid', vt_il  * zimskt      )                                           ! melt pond lid total volume per unit area
         ! salt
         IF( iom_use('icesalt' ) ) THEN
            !$acc update self( sm_i )
            CALL iom_put( 'icesalt', sm_i                 * zimskt ) !+ zmiss_val * ( 1._wp - zimskt ) ) ! mean ice salinity
         ENDIF
         IF( iom_use('icesalm' ) )   CALL iom_put( 'icesalm', st_i * rhoi * 1.0e-3 * zimskt )                                  ! Mass of salt in sea ice per cell area

         ! heat         
         IF( iom_use('icetemp' ) ) CALL iom_put( 'icetemp', ( tm_i  - rt0 ) ) !* zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! ice mean temperature

         IF( iom_use('e3t_m'   ) )   CALL iom_put( 'e3t_m',  e3t_m(:,:) )
         IF( iom_use('frq_m'   ) )   CALL iom_put( 'frq_m',  frq_m(:,:) )

         
         IF( iom_use('snwtemp' ) )   CALL iom_put( 'snwtemp', ( tm_s  - rt0 ) * zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! snw mean temperature
         IF( iom_use('icettop' ) )   CALL iom_put( 'icettop', ( tm_su - rt0 ) * zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! temperature at the ice surface
         IF( iom_use('icetbot' ) )   CALL iom_put( 'icetbot', ( t_bo  - rt0 ) * zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! temperature at the ice bottom
         IF( iom_use('dtsrfbt' ) )   CALL iom_put( 'dtsrfbt', ( t_bo-(sst_m+rt0))*zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! `t_bo - sst`
         IF( iom_use('icetsni' ) )   CALL iom_put( 'icetsni', ( tm_si - rt0 ) * zimskt + zmiss_val * ( 1._wp - zimskt ) )      ! temperature at the snow-ice interface
         IF( iom_use('icehc'   ) )   CALL iom_put( 'icehc'  ,  -et_i          * zimskt )                                       ! ice heat content
         IF( iom_use('snwhc'   ) )   CALL iom_put( 'snwhc'  ,  -et_s          * zimskt )                                       ! snow heat content


         IF( iom_use('icealb') .OR. iom_use('albedo') ) THEN                                                                   ! ice albedo and surface albedo
            ALLOCATE( zalb(jpi,jpj), zmsktalb(jpi,jpj) )
            ! ice albedo
            WHERE( at_i_b < 1.e-03 )
               zmsktalb(:,:) = 0._wp
               zalb   (:,:) = rn_alb_oce
            ELSEWHERE
               zmsktalb(:,:) = 1._wp
               zalb   (:,:) = SUM( alb_ice * a_i_b, dim=3 ) / at_i_b
            END WHERE
            CALL iom_put( 'icealb' , zalb * zmsktalb + zmiss_val * ( 1._wp - zmsktalb ) )
            ! ice+ocean albedo
            zalb(:,:) = SUM( alb_ice * a_i_b, dim=3 ) + rn_alb_oce * ( 1._wp - at_i_b )
            CALL iom_put( 'albedo' , zalb )
            DEALLOCATE( zalb, zmsktalb )
         ENDIF

         ! --- category-dependent fields --- !
         IF( iom_use('icemask_cat' ) )   CALL iom_put( 'icemask_cat' ,                  zmskt00l                                   ) ! ice mask 0%
         IF( iom_use('iceconc_cat' ) )   CALL iom_put( 'iceconc_cat' , a_i            * zmskt00l                                   ) ! area for categories
         IF( iom_use('icethic_cat' ) )   CALL iom_put( 'icethic_cat' , h_i            * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! thickness for categories
         IF( iom_use('snwthic_cat' ) )   CALL iom_put( 'snwthic_cat' , h_s            * zmsktsnl + zmiss_val * ( 1._wp - zmsktsnl ) ) ! snow depth for categories
         IF( iom_use('icesalt_cat' ) )   CALL iom_put( 'icesalt_cat' , s_i            * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! salinity for categories
         IF( iom_use('iceage_cat'  ) )   CALL iom_put( 'iceage_cat'  , o_i / rday     * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! ice age
         IF( iom_use('icetemp_cat' ) )   CALL iom_put( 'icetemp_cat' , ( SUM( t_i, dim=3 ) * r1_nlay_i - rt0 ) &
            &                                                                         * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! ice temperature
         IF( iom_use('snwtemp_cat' ) )   CALL iom_put( 'snwtemp_cat' , ( SUM( t_s, dim=3 ) * r1_nlay_s - rt0 ) &
            &                                                                         * zmsktsnl + zmiss_val * ( 1._wp - zmsktsnl ) ) ! snow temperature
         IF( iom_use('icettop_cat' ) )   CALL iom_put( 'icettop_cat' , ( t_su - rt0 ) * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! surface temperature
         !IF( iom_use('icebrv_cat'  ) )   CALL iom_put( 'icebrv_cat'  ,   v_ibr * 100.  * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! brine volume
         !IF( iom_use('iceapnd_cat' ) )   CALL iom_put( 'iceapnd_cat' ,   a_ip         * zmskt00l                                   ) ! melt pond frac for categories
         !IF( iom_use('icevpnd_cat' ) )   CALL iom_put( 'icevpnd_cat' ,   v_ip         * zmskt00l                                   ) ! melt pond volume for categories
         !IF( iom_use('icehpnd_cat' ) )   CALL iom_put( 'icehpnd_cat' ,   h_ip         * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! melt pond thickness for categories
         !IF( iom_use('icehlid_cat' ) )   CALL iom_put( 'icehlid_cat' ,   h_il         * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! melt pond lid thickness for categories
         !IF( iom_use('iceafpnd_cat') )   CALL iom_put( 'iceafpnd_cat',   a_ip_frac    * zmskt00l                                   ) ! melt pond frac per ice area for categories
         !IF( iom_use('iceaepnd_cat') )   CALL iom_put( 'iceaepnd_cat',   a_ip_eff     * zmskt00l                                   ) ! melt pond effective frac for categories
         !IF( iom_use('icealb_cat'  ) )   CALL iom_put( 'icealb_cat'  ,   alb_ice      * zmskt00l + zmiss_val * ( 1._wp - zmskt00l ) ) ! ice albedo for categories

         !------------------
         ! Add-ons for SIMIP
         !------------------
         ! trends
         IF( iom_use('dmithd') )   CALL iom_put( 'dmithd', - wfx_bog - wfx_bom - wfx_sum - wfx_sni - wfx_opw - wfx_lam - wfx_res ) ! Sea-ice mass change from thermodynamics
         IF( iom_use('dmidyn') )   CALL iom_put( 'dmidyn', - wfx_dyn + rhoi * diag_trp_vi                                        ) ! Sea-ice mass change from dynamics(kg/m2/s)
         IF( iom_use('dmiopw') )   CALL iom_put( 'dmiopw', - wfx_opw                                                             ) ! Sea-ice mass change through growth in open water
         IF( iom_use('dmibog') )   CALL iom_put( 'dmibog', - wfx_bog                                                             ) ! Sea-ice mass change through basal growth
         IF( iom_use('dmisni') )   CALL iom_put( 'dmisni', - wfx_sni                                                             ) ! Sea-ice mass change through snow-to-ice conversion
         IF( iom_use('dmisum') )   CALL iom_put( 'dmisum', - wfx_sum                                                             ) ! Sea-ice mass change through surface melting
         IF( iom_use('dmibom') )   CALL iom_put( 'dmibom', - wfx_bom                                                             ) ! Sea-ice mass change through bottom melting
         IF( iom_use('dmilam') )   CALL iom_put( 'dmilam', - wfx_lam                                                             ) ! Sea-ice mass change through lateral melting
         IF( iom_use('dmtsub') )   CALL iom_put( 'dmtsub', - wfx_sub                                                             ) ! Sea-ice mass change through evaporation and sublimation
         IF( iom_use('dmssub') )   CALL iom_put( 'dmssub', - wfx_snw_sub                                                         ) ! Snow mass change through sublimation
         IF( iom_use('dmisub') )   CALL iom_put( 'dmisub', - wfx_ice_sub                                                         ) ! Sea-ice mass change through sublimation
         IF( iom_use('dmsspr') )   CALL iom_put( 'dmsspr', - wfx_spr                                                             ) ! Snow mass change through snow fall
         IF( iom_use('dmsssi') )   CALL iom_put( 'dmsssi',   wfx_sni*rhos*r1_rhoi                                                ) ! Snow mass change through snow-to-ice conversion
         IF( iom_use('dmsmel') )   CALL iom_put( 'dmsmel', - wfx_snw_sum                                                         ) ! Snow mass change through melt
         IF( iom_use('dmsdyn') )   CALL iom_put( 'dmsdyn', - wfx_snw_dyn + rhos * diag_trp_vs                                    ) ! Snow mass change through dynamics(kg/m2/s)

         ! Global ice diagnostics
         IF(  iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt') .OR. &
            & iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt') ) THEN
            !
            z2d(:,:) = MERGE( 1._wp, 0._wp,  ff_u(:,:) > 0._wp )
            !
            IF( iom_use('NH_icearea') )   zdiag_area_nh = glob_sum( 'icewri', at_i *           z2d   * e1e2t * 1.e-12 )
            IF( iom_use('NH_icevolu') )   zdiag_volu_nh = glob_sum( 'icewri', vt_i *           z2d   * e1e2t * 1.e-12 )
            IF( iom_use('NH_iceextt') )   zdiag_extt_nh = glob_sum( 'icewri',                  z2d   * e1e2t * 1.e-12 * zimskt )
            !
            IF( iom_use('SH_icearea') )   zdiag_area_sh = glob_sum( 'icewri', at_i * ( 1._wp - z2d ) * e1e2t * 1.e-12 )
            IF( iom_use('SH_icevolu') )   zdiag_volu_sh = glob_sum( 'icewri', vt_i * ( 1._wp - z2d ) * e1e2t * 1.e-12 )
            IF( iom_use('SH_iceextt') )   zdiag_extt_sh = glob_sum( 'icewri',        ( 1._wp - z2d ) * e1e2t * 1.e-12 * zimskt )
            !
            CALL iom_put( 'NH_icearea' , zdiag_area_nh )
            CALL iom_put( 'NH_icevolu' , zdiag_volu_nh )
            CALL iom_put( 'NH_iceextt' , zdiag_extt_nh )
            CALL iom_put( 'SH_icearea' , zdiag_area_sh )
            CALL iom_put( 'SH_icevolu' , zdiag_volu_sh )
            CALL iom_put( 'SH_iceextt' , zdiag_extt_sh )
            !
         ENDIF

         DEALLOCATE( zmskt00l, zmsktsnl )
         !
         !-------------------------------------------------------------------------------------------------
      ENDIF !IF( ln_icethd )



      !$acc end data
      IF( ln_timing )  CALL timing_stop('ice_wri')
      !
   END SUBROUTINE ice_wri






   SUBROUTINE ice_wri_adv( kt )
      !!-------------------------------------------------------------------
      !!  This routine ouputs `dynADV2D`-related fields
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time-step
      REAL(wp), DIMENSION(jpi,jpj) :: zimskt, zimskf
      !!-------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_wri_adv')
      
      ! Ice-extent masks at T- & F-points:
      !$acc update self( kmsk_ice_t, kmsk_ice_f )
      zimskt(:,:) = REAL(kmsk_ice_t,wp)
      zimskf(:,:) = REAL(kmsk_ice_f,wp)

      !-----------------
      ! Standard outputs
      !-----------------

      ! masks
      !CALL iom_put( 'icemask'  , zimskt )   ! ice mask 0%
      
      IF( iom_use('iceconc' ) ) THEN
         !$acc update self( at_i )
         CALL iom_put( 'iceconc',   at_i )    ! ice concentration
      ENDIF
      IF( iom_use('iceconc-f') )  THEN
         !$acc update self( af_i )
         CALL iom_put( 'iceconc-f', af_i )    ! ice concentration
      ENDIF
      IF( iom_use('icevolu' ) ) THEN
         CALL iom_put( 'icevolu', vt_i * zimskt )       ! ice volume = mean ice thickness over the cell
      ENDIF
      IF( iom_use('icethic' ) ) THEN
         !$acc update self( hm_i )
         CALL iom_put( 'icethic', hm_i * zimskt )   ! ice thickness
      ENDIF
      
      IF( ln_damage ) THEN
         IF( iom_use('icedmgt') ) THEN
            !$acc update self( dmdt )
            CALL iom_put( 'icedmgt', (1._wp - dmdt) * zimskt )          ! ice damage @T !
         ENDIF
         IF( iom_use('icedmgf') ) THEN
            !$acc update self( dmdf )
            CALL iom_put( 'icedmgf', (1._wp - dmdf) * zimskf )          ! ice damage @T !
         ENDIF
      END IF
      
      ! momentum
      IF( iom_use('uice'    ) ) THEN
         CALL iom_put( 'uice'   , u_ice ) !* zimskt   )                                            ! ice velocity u
      ENDIF
      IF( iom_use('vice'    ) ) THEN
         CALL iom_put( 'vice'   , v_ice ) !* zimskt   )                                            ! ice velocity v
      ENDIF
      IF( iom_use('uice_v'  ) ) THEN
         !$acc update self( uVice )
         CALL iom_put( 'uice_v' , uVice   )
      ENDIF
      IF( iom_use('vice_u'  ) ) THEN
         !$acc update self( vUice )
         CALL iom_put( 'vice_u' , vUice   )
      ENDIF
      
      IF( ln_timing )  CALL timing_stop('ice_wri_adv')
      !
   END SUBROUTINE ice_wri_adv












   

   SUBROUTINE ice_wri_state( kid )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wri_state  ***
      !!
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :   4.0  !  2013-06  (C. Rousset)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kid
      !!----------------------------------------------------------------------
      !
      !! The file is open in dia_wri_state (ocean routine)


      !$acc update self( hm_i, at_i, tm_i, u_ice, v_ice, qsr, qns, fatm_snow, sm_i, vt_i, h_i, a_i, s_i, h_s )

      CALL iom_rstput( 0, 0, kid, 'sithic', hm_i         )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconc', at_i         )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sitemp', tm_i - rt0   )   ! Ice temperature
      CALL iom_rstput( 0, 0, kid, 'sivelu', u_ice        )   ! i-Ice speed
      CALL iom_rstput( 0, 0, kid, 'sivelv', v_ice        )   ! j-Ice speed
      CALL iom_rstput( 0, 0, kid, 'sisflx', qsr          )   ! Solar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'sinflx', qns          )   ! NonSolar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'snwpre', fatm_snow    )   ! Snow precipitation
      CALL iom_rstput( 0, 0, kid, 'sisali', sm_i         )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'sivolu', vt_i         )   ! Ice volume
      !CALL iom_rstput( 0, 0, kid, 'sidive', divu_i*1.0e8 )   ! Ice divergence
      CALL iom_rstput( 0, 0, kid, 'si_amp', at_ip        )   ! Melt pond fraction
      CALL iom_rstput( 0, 0, kid, 'si_vmp', vt_ip        )   ! Melt pond volume
      CALL iom_rstput( 0, 0, kid, 'sithicat', h_i        )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconcat', a_i        )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sisalcat', s_i        )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'snthicat', h_s        )   ! Snw thickness
      IF(ln_damage) THEN
         !$acc update self( dmdt, dmdf, uVice, vUice )
         CALL iom_rstput( 0, 0, kid, 'dmd-t',    dmdt )
         CALL iom_rstput( 0, 0, kid, 'dmd-f',    dmdf )
         CALL iom_rstput( 0, 0, kid, 'sivelu_v', uVice        )
         CALL iom_rstput( 0, 0, kid, 'sivelv_u', vUice        )
      ENDIF

   END SUBROUTINE ice_wri_state

   !!======================================================================
END MODULE icewri
