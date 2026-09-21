MODULE icedyn_rhg_dummy
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_dummy  ***
   !!   Sea-Ice dynamics : rheology Britle Maxwell X
   !!======================================================================
   !! History :
   !!            4.2  !  2022     (L. Brodeau) `BBM`
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_dummy : computes ice velocities from BBM rheology
   !!----------------------------------------------------------------------
   USE dom_oce        ! Ocean domain
   USE ice,      ONLY : u_ice, v_ice, uVice, vUice
   USE oss_nnq,  ONLY : ssu_m, ssv_m
   USE icevar,   ONLY : ice_var_sshdyn
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary conditions (or mpp links)
#if defined _OPENACC || defined _OPENMP
   USE lbclnk_gpu
#endif
   !USE prtctl         ! Print control
   !
   USE timing

   !USE ice_util

   USE icedyn_rhg_tools, ONLY : strain_rate_dsd

   USE remap_classic, ONLY : do_rmpU2V, do_rmpV2U


   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_dummy       ! called by icedyn_rhg.F90



   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! $Id: icedyn_rhg_dummy.F90 13646 2020-10-20 15:33:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE ice_dyn_rhg_dummy( kt, pshear_i, pdivu_i, pdelta_i )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_dummy  ***
      !!                             BBM-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  the BBM rheology of Olason et al., 2022.
      !!
      !! ** Inputs  : - wind forcing (stress), oceanic currents
      !!                ice total volume (vt_i) per unit area
      !!                snow total volume (vt_s) per unit area
      !!
      !! ** Action  : - compute u_ice, v_ice : the components of the
      !!                sea-ice velocity vector
      !!              - compute delta_i, shear_i, divu_i, which are inputs
      !!                of the ice thickness distribution
      !!
      !! ** Steps   : 0) compute mask at F point
      !!              1) Compute ice snow mass, ice strength
      !!              2) Compute wind, oceanic stresses, mass terms and
      !!                 coriolis terms of the momentum equation
      !!              3) Solve the momentum equation (iterative procedure)
      !!              4) Recompute delta, shear and divergence
      !!                 (which are inputs of the ITD) & store stress
      !!                 for the next time step
      !!              5) Diagnostics including charge ellipse
      !!
      !! ** Notes   :
      !!
      !!
      !!
      !!
      !! References : Brodeau et al., 2024, GMD
      !!              Olason et al., 2022
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in ) :: kt                                    ! time step
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pshear_i, pdivu_i, pdelta_i      !
      !!-------------------------------------------------------------------
      !
      INTEGER ::   ji, jj       ! dummy loop indices

      !!-------------------------------------------------------------------
      !$acc data present( pshear_i, pdivu_i, pdelta_i, u_ice, v_ice, uVice, vUice, ssu_m, ssv_m )

      IF( ln_timing )   CALL timing_start('ice_dyn_rhg_dummy')

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_rhg_dummy: NO sea-ice rheology!'

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            u_ice(ji,jj) = ssu_m(ji,jj)
            v_ice(ji,jj) = ssv_m(ji,jj)
         END DO
      END DO
      !$acc end parallel loop
      CALL do_rmpU2V( u_ice,  uVice )
      CALL do_rmpV2U( v_ice,  vUice )
      CALL lbc_lnk( 'icedyn_rhg_dummy', u_ice,'U',-1._wp, v_ice,'V',-1._wp, uVice,'V',-1._wp, vUice,'U',-1._wp )

      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!
      CALL strain_rate_dsd( 'T', u_ice, v_ice, uVice, vUice, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e1t2, e2t2, xmskt, &
         &                       pdivu_i, pshear_i, pdelta_i )
#if defined _OPENACC || defined _OPENMP
      CALL lbc_lnk_gpu( 'icedyn_rhg_dummy', pdivu_i, pshear_i, pdelta_i )
#else
      CALL lbc_lnk(     'icedyn_rhg_dummy', pdivu_i,'T',1._wp, pshear_i,'T',1._wp, pdelta_i,'T',1._wp )
#endif

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_dyn_rhg_dummy')
      !
   END SUBROUTINE ice_dyn_rhg_dummy

   !!==============================================================================
END MODULE icedyn_rhg_dummy
