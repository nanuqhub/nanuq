MODULE icedyn_adv
   !!======================================================================
   !!                       ***  MODULE icedyn_adv   ***
   !!   sea-ice: advection
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv   : advection of sea ice variables
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice            ! sea-ice: variables
   !
   USE ice_util, ONLY : cap_1md
   !
   USE icedyn_adv_pra ! sea-ice: advection scheme (Prather)
   USE icedyn_adv_umx ! sea-ice: advection scheme (ultimate-macho)
   USE icedyn_adv_wnx ! sea-ice: advection scheme WENO[5/7] => needs nn_hls=2 for WENO5, nn_hls=3 for WENO7 !!!
   !
   USE icedyn_adv_pra_d ! sea-ice: advection scheme (Prather) for damage @T or @F
   USE icedyn_adv_wnx_d ! sea-ice: advection scheme WENO[5/7] for damage @T or @F
   !
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library

   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv        ! called by icestp
   PUBLIC   ice_dyn_adv_init   ! called by icedyn

   !! Advection of non-brittle fields:
   INTEGER ::              nice_adv   ! choice of the type of advection scheme
   !                                       ! associated indices:
   INTEGER, PARAMETER ::   np_advPRA = 1   ! Prather scheme
   INTEGER, PARAMETER ::   np_advUMx = 2   ! Ultimate-Macho scheme
   INTEGER, PARAMETER ::   np_advWNx = 3   ! WENO scheme
   !
   !! Advection of damage (and stress tensor components):
   INTEGER ::              nice_d_adv   ! choice of the type of advection scheme
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_d_advPRA = 1   ! Prather scheme
   INTEGER, PARAMETER ::   np_d_advWNx = 2   ! WENO scheme

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icedyn_adv.F90 13472 2020-09-16 13:05:19Z smasson $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_dyn_adv ***
      !!
      !! ** purpose : advection of sea ice
      !!
      !! ** method  : One can choose between
      !!     a) an Ultimate-Macho scheme (with order defined by nn_UMx) => ln_adv_UMx
      !!     b) and a second order Prather scheme => ln_adv_Pra
      !!
      !! ** action :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zdudx, zdvdy, zdudy, zdvdx, zdiv
      REAL(wp), DIMENSION(jpi,jpj) :: zs11, zs22, zs12
      REAL(wp), DIMENSION(jpi,jpj) :: zs11_ci, zs22_ci, zs12_ci ! increments for extra upper- or lower-convected terms...
      REAL(wp)                     :: zm1t, zm1f,zm2t, zm2f, zm, zk
      REAL(wp)                     :: z_skk_ao, z_s12_ao
      INTEGER :: ji, jj
      LOGICAL :: lwpv, lwpw
      CHARACTER(len=1) :: cgt
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('icedyn_adv')
      !$acc data present( sudy_u, svdx_v, u_ice, v_ice, e2u, e1v, sudy_v, svdx_u, uVice, vUice, e2v, e1u )

      lwpv = ( lwp .AND. (iverbose>0) )
      lwpw = ( lwp .AND. (iverbose>1) )

      ! controls
      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv: sea-ice advection'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      !---------------!
      !== Advection ==!
      !---------------!

      ! Compute ice velocity transports:
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            sudy_u(ji,jj) = u_ice(ji,jj) * e2u(ji,jj) * 1.E-6_wp    ! `1.E-6_wp` because we work with km^2 rather than m^2 !!!
            svdx_v(ji,jj) = v_ice(ji,jj) * e1v(ji,jj) * 1.E-6_wp    !       "                     "
         END DO
      END DO
      !$acc end parallel loop

      IF( iom_use('u_ice_trsprt') ) THEN
         !$acc update self( sudy_u )
         CALL iom_put( 'u_ice_trsprt' , sudy_u )
      ENDIF
      IF( iom_use('v_ice_trsprt') ) THEN
         !$acc update self( svdx_v )
         CALL iom_put( 'v_ice_trsprt' , svdx_v )
      ENDIF


      IF(lwpv) WRITE(numout,*) ''
      SELECT CASE( nice_adv )
         !                             !-----------------------!
      CASE( np_advUMx )                ! ULTIMATE-MACHO scheme !
         !                             !-----------------------!
         IF(lwpv) WRITE(numout,'("  *** advects GENERIC fields @T with UMX, order = ",i1," kt=",i)') nn_UMx, kt
         CALL ice_dyn_adv_umx( nn_UMx, kt, sudy_u, svdx_v, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i, szv_i )
         !
         !                             !-----------------------!
      CASE( np_advPRA )                ! PRATHER scheme        !
         !                             !-----------------------!
         IF(lwpv) WRITE(numout,'("  *** advects GENERIC fields @T with Prather, kt=",i)') kt
         !
         CALL ice_dyn_adv_pra(         kt, xmskt, sudy_u, svdx_v, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i, szv_i )
         !
         !                             !-----------------------!
      CASE( np_advWNx )                ! WENOX scheme          !
         !                             !-----------------------!
         IF(lwpv) WRITE(numout,'("  *** advects GENERIC fields @T with WENO",i1,", kt=",i)') nn_WNx, kt
         CALL ice_dyn_adv_wnx(         kt, sudy_u, svdx_v, klbct, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i, szv_i )
         !
      END SELECT


      IF( ln_damage .AND. nn_d_adv >= 1 ) THEN

         !! Advect at T points with u@U, v@V velocities
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         cgt = 'T'
         !
         SELECT CASE( nice_d_adv )
         CASE( np_d_advPRA )
            IF(lwpv) WRITE(numout,'("  *** advects only damage @'//cgt//' with Prather, kt=",i)') kt
            CALL ice_dyn_adv_pra_d( kt, cgt, e1e2t, r1_e1e2t, xmskt, sudy_u, svdx_v,         dmdt )
            !
         CASE( np_d_advWNx )
            IF(lwpv) WRITE(numout,'("  *** advects only damage @'//cgt//' with WENO",i1,", kt=",i)') nn_WNx, kt
            CALL ice_dyn_adv_wnx_d( kt, cgt, e1e2t, r1_e1e2t, xmskt, sudy_u, svdx_v, klbct,  dmdt )
            !
         END SELECT

         !! Advect at F points with u@V, v@U velocities
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         cgt = 'F'
         !
         ! Compute transports:
         ! (MIND: in the F-centric mesh, to have U pt at the RHS of the center point & V above, need to use `i+1` & `j+1`, respectively)
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls-1
            DO ji=Nis0-nn_hls, Nie0+nn_hls-1
               sudy_v(ji,jj) = uVice(ji+1,jj) * e2v(ji+1,jj) * 1.E-6_wp   ! `1.E-6_wp` because we work with km^2 rather than m^2 !!!
               svdx_u(ji,jj) = vUice(ji,jj+1) * e1u(ji,jj+1) * 1.E-6_wp   !       "                     "
            END DO
         END DO
         !$acc end parallel loop

         IF( iom_use('uVice_trsprt') ) THEN
            !$acc update self( sudy_v )
            CALL iom_put( 'uVice_trsprt' , sudy_v )
         ENDIF
         IF( iom_use('vUice_trsprt') ) THEN
            !$acc update self( svdx_u )
            CALL iom_put( 'vUice_trsprt' , svdx_u )
         ENDIF

         SELECT CASE( nice_d_adv )
         CASE( np_d_advPRA )
            IF(lwpv) WRITE(numout,'("  *** advects only damage @'//cgt//' with Prather, kt=",i)') kt
            CALL ice_dyn_adv_pra_d( kt, cgt, e1e2f, r1_e1e2f, xmskf, sudy_v, svdx_u,         dmdf )
            !
         CASE( np_d_advWNx )
            IF(lwpv) WRITE(numout,'("  *** advects only damage @'//cgt//' with WENO",i1,", kt=",i)') nn_WNx, kt
            CALL ice_dyn_adv_wnx_d( kt, cgt, e1e2f, r1_e1e2f, xmskf, sudy_v, svdx_u, klbcf,  dmdf )
            !
            !CASE( np_d_advUMx )
            !   CALL ctl_stop( 'ice_dyn_adv: UMX advection not ready yet to advect at F-points!')
            !
         END SELECT

         IF( .NOT. ln_pureADV2D )  CALL cap_1md( at_i, af_i, dmdt, dmdf )   ! Mind that they are `at_i` & `af_i` prior to advection...

      END IF !IF( ln_damage .AND. nn_d_adv >= 1 ) THEN

      !------------
      ! diagnostics
      !------------
      !diag_trp_ei(:,:) = SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice
      !diag_trp_es(:,:) = SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
      !diag_trp_sv(:,:) = SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice
      !diag_trp_vi(:,:) = SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice
      !diag_trp_vs(:,:) = SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice
      !IF( iom_use('icemtrp') )   CALL iom_put( 'icemtrp' ,  diag_trp_vi * rhoi          )   ! ice mass transport
      !IF( iom_use('snwmtrp') )   CALL iom_put( 'snwmtrp' ,  diag_trp_vs * rhos          )   ! snw mass transport
      !IF( iom_use('salmtrp') )   CALL iom_put( 'salmtrp' ,  diag_trp_sv * rhoi * 1.e-03 )   ! salt mass transport (kg/m2/s)
      !IF( iom_use('dihctrp') )   CALL iom_put( 'dihctrp' , -diag_trp_ei                 )   ! advected ice heat content (W/m2)
      !IF( iom_use('dshctrp') )   CALL iom_put( 'dshctrp' , -diag_trp_es                 )   ! advected snw heat content (W/m2)

      ! controls
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn & trp - ')                           ! prints

      !$acc end data
      IF( ln_timing)   CALL timing_stop ('icedyn_adv')
      !
   END SUBROUTINE ice_dyn_adv


   SUBROUTINE ice_dyn_adv_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_init  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice
      !!                dynamics
      !!
      !! ** Method  :   Read the namdyn_adv namelist and check the ice-dynamic
      !!                parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_adv
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_adv/ ln_adv_Pra, ln_adv_UMx, nn_UMx, ln_adv_WNx, nn_WNx, cn_weno_wght
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn_adv)
      !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_adv in reference namelist' )
      READ_NML_CFG(numnam_ice,namdyn_adv)
      !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_adv in configuration namelist' )
      IF(lwm) WRITE( numoni, namdyn_adv )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_adv:'
         WRITE(numout,*) '      type of advection scheme (Prather)             ln_adv_Pra = ', ln_adv_Pra
         WRITE(numout,*) '      type of advection scheme (Ulimate-Macho)       ln_adv_UMx = ', ln_adv_UMx
         WRITE(numout,*) '         order of the Ultimate-Macho scheme          nn_UMx     = ', nn_UMx
         WRITE(numout,*) '      type of advection scheme (WENO)                ln_adv_WNx = ', ln_adv_WNx
         WRITE(numout,*) '         order of the WENO scheme                    nn_WNx     = ', nn_WNx
         WRITE(numout,*) '         file in which WENO weights are read       cn_weno_wght = ', TRIM(cn_weno_wght)
      ENDIF
      !
      !                             !== set the choice of ice advection ==!
      ioptio = 0
      IF( ln_adv_Pra ) THEN
         ioptio = ioptio + 1
         nice_adv = np_advPRA
      ENDIF
      IF( ln_adv_UMx ) THEN
         ioptio = ioptio + 1
         nice_adv = np_advUMx
# if defined _OPENACC
         CALL ctl_stop( 'ice_dyn_adv_init: UMX advection has not been adapted yet for GPU!' )
# endif
      ENDIF
      IF( ln_adv_WNx ) THEN
         IF(     nn_WNx == 5 ) THEN
            IF( nn_hls<2 ) CALL ctl_stop( 'ice_dyn_adv_init: with WENO5 you need `nn_hls>=2`' )
         ELSEIF( nn_WNx == 7 ) THEN
            IF( nn_hls<3 ) CALL ctl_stop( 'ice_dyn_adv_init: with WENO7 you need `nn_hls>=3`' )
         ELSE
            CALL ctl_stop( 'ice_dyn_adv_init: for WENO, only order 5 or 7 are available for now' )
         ENDIF
         ioptio = ioptio + 1
         nice_adv = np_advWNx
      ENDIF
      !
      IF( ioptio < 1 ) CALL ctl_stop( 'ice_dyn_adv_init: pick one advection scheme for tracers! (ln_adv_Pra/ln_adv_UMx/ln_adv_WNx)' )
      IF( ioptio > 1 ) CALL ctl_stop( 'ice_dyn_adv_init: pick ONLY one advection scheme for tracers! (ln_adv_Pra/ln_adv_UMx/ln_adv_WNx)' )
      !
      IF( ln_adv_Pra ) CALL adv_pra_init()  !* read or initialize all required files
      !
      IF( ln_adv_UMx ) CALL adv_umx_init()
      !
      IF( ln_adv_WNx .OR. ln_adv_d_wnx ) CALL adv_wnx_init()
      !
      IF( ln_damage ) THEN
         IF(nn_d_adv > 1) CALL ctl_stop( 'ice_dyn_adv_init: Advection of stress tensors not supported anymore as it made no sense (nn_d_adv)' )

         IF( ln_dynADV2D ) THEN
            !! Idealized test-cases to test sea-ice advection scheme (on both T- and F-points...)
            IF( nn_d_adv/=1 ) THEN
               ! This is an idealized advection test-case and `ln_damage=T`, this means that we want to test the advection of damage
               ! at both T- and F-points !
               nn_d_adv = 1
               CALL ctl_warn( 'ice_dyn_adv_init: forcing `nn_d_adv` to 1! Damage at T and F to be advected!' )
            ENDIF
            IF( ln_adv_WNx )  ln_adv_d_wnx = .TRUE.
            IF( ln_adv_Pra )  ln_adv_d_pra = .TRUE.
         ENDIF

         IF(ln_adv_d_pra) THEN
            nice_d_adv = np_d_advPRA
            CALL adv_pra_d_init()
         ENDIF

         IF( (ln_adv_d_wnx).AND.(nn_d_adv>0) ) THEN
            IF( .NOT. ln_adv_WNx ) CALL ctl_stop( 'ice_dyn_adv_init: with `ln_adv_d_wnx=T`, use `ln_adv_WNx=T`!' )
            nice_d_adv = np_d_advWNx
         ENDIF

      END IF

      !$acc update device ( nn_UMx, nn_WNx, cn_weno_wght )

   END SUBROUTINE ice_dyn_adv_init

   !SUBROUTINE lower_convected_inc( pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk, ps11, ps22, ps12,  pinc11, pinc22, pinc12 )
   !   REAL(wp), DIMENSION(:,:), INTENT(in)  :: pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk
   !   REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11, ps22, ps12 ! tensor components before any form of advection
   !   REAL(wp), DIMENSION(:,:), INTENT(out) :: pinc11, pinc22, pinc12  ! lower-convected contribution increment
   !   !! Lower-convected version (sign is inversed because moved from LHS to RHS):
   !   pinc11(:,:) = -2._wp*rDt_ice*( pdudx(:,:)*ps11(:,:) + pdvdx(:,:)*ps12(:,:) )*zmsk(:,:)
   !   pinc22(:,:) = -2._wp*rDt_ice*( pdvdy(:,:)*ps22(:,:) + pdudy(:,:)*ps12(:,:) )*zmsk(:,:)
   !   pinc12(:,:) = -rDt_ice*( pdiv(:,:)*ps12(:,:) + pdudy(:,:)*ps11(:,:) + pdvdx(:,:)*ps22(:,:) )*zmsk(:,:)
   !END SUBROUTINE lower_convected_inc

   !SUBROUTINE upper_convected_inc( pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk, ps11, ps22, ps12,  pinc11, pinc22, pinc12 )
   !   REAL(wp), DIMENSION(:,:), INTENT(in)  :: pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk
   !   REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11, ps22, ps12 ! tensor components before any form of advection
   !   REAL(wp), DIMENSION(:,:), INTENT(out) :: pinc11, pinc22, pinc12  ! lower-convected contribution increment
   !   !! Upper-convected version (sign is inversed because moved from LHS to RHS):
   !   pinc11(:,:) =   2._wp*rDt_ice*( pdudx(:,:)*ps11(:,:) + pdudy(:,:)*ps12(:,:) )*zmsk(:,:)
   !   pinc22(:,:) =   2._wp*rDt_ice*( pdvdy(:,:)*ps22(:,:) + pdvdx(:,:)*ps12(:,:) )*zmsk(:,:)
   !   pinc12(:,:) =   rDt_ice*( pdiv(:,:)*ps12(:,:) + pdudy(:,:)*ps22(:,:) + pdvdx(:,:)*ps11(:,:) )*zmsk(:,:)
   !END SUBROUTINE upper_convected_inc

   !!======================================================================
END MODULE icedyn_adv
