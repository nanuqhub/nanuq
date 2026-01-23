MODULE icebbc
   !!======================================================================
   !!                       ***  MODULE  icebbc  ***
   !! Sea-Ice :   liquid-ocean/ice bbc fields
   !!=====================================================================
   !! History :  L. Brodeau, 2024/12
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE phycst,  ONLY : rho0
   USE dom_oce, ONLY : umask, vmask
   USE par_ice, ONLY : rn_Cd_io
   USE ice            ! sea-ice: variables
   USE remap_classic, ONLY : rmpU2V, rmpV2U
   USE in_out_manager ! I/O manager
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_bbc_tau   ! called by icestp.F90

   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_bbc_tau( kt, ptx_oi_u, pty_oi_v,  ptx_oi_v, pty_oi_u )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_bbc_tau  ***
      !!
      !! ** Purpose : provide bottom boundary condition for sea ice (momentum)
      !!
      !! ** Action  : It provides the following fields:
      !!              ptx_oi_u, pty_oi_v : bottom ice stress (U- & V-points) [N/m2]
      !!-------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kt                   ! ocean time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptx_oi_u, pty_oi_v   ! liquid-ocean/ice stress   [N/m2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out), OPTIONAL :: ptx_oi_v, pty_oi_u   ! liquid-ocean/ice stress   [N/m2] !#bbm
      !!
      INTEGER  ::   ji, jj                 ! dummy loop index
      REAL(wp) :: zrhoco, zTauO, zt1, zt2, zUi, zUo, zVi, zVo
      LOGICAL :: lEgrid
      !!-------------------------------------------------------------------
      !$acc data present( V_oce(:,:,1:4), u_ice, v_ice, uVice, vUice, ptx_oi_u, pty_oi_v )
      
      lEgrid = ( PRESENT(ptx_oi_v) .AND. PRESENT(pty_oi_u) )
      IF( ln_timing )   CALL timing_start('ice_bbc_tau')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_bbc_tau: Bottom boundary condition for sea ice (momentum)'
         WRITE(numout,*)'~~~~~~~~~~~~~~~'
      ENDIF
      !
      zrhoco   = rho0 * rn_Cd_io

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            ! These 2 are computed both by BBM and EVP:
            !vUice(ji,jj) = 0.25_wp * ( (v_ice(ji,jj) + v_ice(ji,jj-1)) + (v_ice(ji+1,jj) + v_ice(ji+1,jj-1)) ) * umask(ji,jj,1)
            !uVice(ji,jj) = 0.25_wp * ( (u_ice(ji,jj) + u_ice(ji-1,jj)) + (u_ice(ji,jj+1) + u_ice(ji-1,jj+1)) ) * vmask(ji,jj,1)

            !! X-component of stress at U-points:
            zUi = u_ice(ji,jj)
            zUo = V_oce(ji,jj,1)
            zt1 = zUi - zUo
            zt2 = vUice(ji,jj) - V_oce(ji,jj,4)
            zTauO = zrhoco * SQRT( zt1*zt1 + zt2*zt2 )
            ptx_oi_u = zTauO * ( zUo - zUi )

            !! Y-component of stress at V-points:
            zVi = v_ice(ji,jj)
            zVo = V_oce(ji,jj,2)
            zt1 = zVi - zVo
            zt2 = uVice(ji,jj) - V_oce(ji,jj,3)
            zTauO = zrhoco * SQRT( zt1*zt1 + zt2*zt2 )
            pty_oi_v = zTauO * ( zVo - zVi )

         ENDDO
      ENDDO
      !$acc end parallel loop

      IF( lEgrid ) THEN
         !$acc parallel loop collapse(2) present( ptx_oi_v, pty_oi_u )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               !! X-component of stress at V-points:
               zUi  = uVice(ji,jj)
               zUo  = V_oce(ji,jj,3)
               zt1   = zUi - zUo
               zt2   = v_ice(ji,jj) - V_oce(ji,jj,2)
               zTauO = zrhoco * SQRT( zt1*zt1 + zt2*zt2 )
               ptx_oi_v = zTauO * ( zUo - zUi )

               !! Y-component of stress at U-points:
               zVi  = vUice(ji,jj)
               zVo  = V_oce(ji,jj,4)
               zt1 = zVi - zVo
               zt2 = u_ice(ji,jj) - V_oce(ji,jj,1)
               zTauO = zrhoco * SQRT( zt1*zt1 + zt2*zt2 )
               pty_oi_u = zTauO * ( zVo - zVi )

            ENDDO
         ENDDO
         !$acc end parallel loop
      ENDIF
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_bbc_tau')
      !
   END SUBROUTINE ice_bbc_tau

   !!======================================================================
END MODULE icebbc
