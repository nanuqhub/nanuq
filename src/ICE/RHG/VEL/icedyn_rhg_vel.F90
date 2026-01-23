MODULE icedyn_rhg_vel
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_vel  ***
   !!   Sea-Ice dynamics : rheology Britle tools...
   !!======================================================================
   !! History :
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE sbc_oce , ONLY : nn_fsbc
   USE par_ice
   USE ice            ! => taux_ai_v, tauy_ai_u are there
   USE lib_mpp,  ONLY: ctl_stop
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE ice_util
   USE in_out_manager, ONLY : ln_timing
   USE timing

   IMPLICIT NONE

   PRIVATE

   !INTERFACE strain_rate
   !   MODULE PROCEDURE strain_rate_Cgrid, strain_rate_Egrid
   !END INTERFACE strain_rate

   PUBLIC FEEO_velocities
   PUBLIC update_uv_euler_si
   PUBLIC update_uv_rk3

   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE update_uv_euler_si( kts, pdt, pAu, pAv, pmu_dt, pmv_dt, pSt, pSf, pgrdSH, pV_oce, &
      &                                     putau_i, pvtau_i, pm01x, pm01y, pm00x, pm00y, pV )
      !!----------------------------------------------------------------------------------------------
      !!   Semi-Implicit Euler update operator for ice velocities update
      !!----------------------------------------------------------------------------------------------
      INTEGER,                          INTENT(in)    :: kts       ! current small time step
      REAL(wp),                         INTENT(in)    :: pdt       ! (small) time-step [s]
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pAu, pAv  ! Ice concentration
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pmu_dt, pmv_dt ! Mass / dt @U and @V
      REAL(wp),   DIMENSION(jpi,jpj,3), INTENT(in)    :: pSt, pSf  ! the 3 components the VERTICALLY-INTEGRATED stress tensor [N/m^2*m] !
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pgrdSH    ! gradient of SSH
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pV_oce    ! the 4 ocean velocity components
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: putau_i, pvtau_i ! air-ice windstress components at T-points [N/m^2]
      INTEGER(1), DIMENSION(jpi,jpj),   INTENT(in)    :: pm01x, pm01y, pm00x, pm00y
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(inout) :: pV        ! the 4 sea-ice velocity components
      !!
      REAL(wp) :: zA, zUi, zVi, zUo, zVo, zM_dt, zmsk, ztau_ai
      REAL(wp) :: zTauO, zcorio, zt1, zt2, zRHS
      REAL(wp) :: zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_uv_euler_si')
      !$acc data present( ff_u,ff_v,pAu,pAv,pmu_dt,pmv_dt,pSt,pSf,pgrdSH,pV_oce,putau_i,pvtau_i,pm01x,pm01y,pm00x,pm00y,pV,umask,vmask )
      
      IF( MOD(kts,2) == 0 ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0

               !! Divergence of the 2 vertically-integrated stress tensors
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_divs_t.h90"

#              include "icedyn_rhg_vel_divs_f.h90"

               !! Update `Vv_sub` at V-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtVv.h90"
               !!
               !! Update `Vu_sub` at U-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtVu.h90"
               !!
               !! Update `Uu_sub` at U-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtUu.h90"
               !!
               !! Update `Uv_sub` at V-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtUv.h90"

            END DO
         END DO
         !$acc end parallel loop

      ELSE

         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0

               !! Divergence of the 2 vertically-integrated stress tensors
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_divs_t.h90"

#              include "icedyn_rhg_vel_divs_f.h90"

               !! Update `Uu_sub` at U-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtUu.h90"
               !!
               !! Update `Uv_sub` at V-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtUv.h90"
               !!
               !! Update `Vv_sub` at V-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtVv.h90"
               !!
               !! Update `Vu_sub` at U-points
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#              include "icedyn_rhg_vel_updtVu.h90"

            END DO
         END DO
         !$acc end parallel loop

      ENDIF

      !$acc end data
      
      IF( ln_timing )   CALL timing_stop('update_uv_euler_si')
   END SUBROUTINE update_uv_euler_si





   SUBROUTINE update_uv_rk3( kts, pdt, pAu, pAv, pmu_dt, pmv_dt, pSt, pSf, pgrdSH, pV_oce, pTau_ai, pm01x, pm01y, pm00x, pm00y, pV )
      !! ACHTUNG WE HAVE AND NEED `pmu, pmv` and not `pAu_dt, pAv_dt` !!!!
      !!----------------------------------------------------------------------------------------------
      !!   Semi-Implicit Euler update operator for ice velocities update
      !!----------------------------------------------------------------------------------------------
      INTEGER,                        INTENT(in)    :: kts       ! current small time step
      REAL(wp),                       INTENT(in)    :: pdt       ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pAu, pAv  ! Ice concentration
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pmu_dt, pmv_dt ! Mass / dt @U and @V
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(in)    :: pSt, pSf    ! the 3 components the stress tensor
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in)    :: pgrdSH    ! gradient of SSH
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in)    :: pV_oce    ! the 4 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(in)    :: pTau_ai   ! the 4 air-ice wind stress components
      INTEGER(1), DIMENSION(jpi,jpj), INTENT(in)    :: pm01x, pm01y, pm00x, pm00y
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(inout) :: pV        ! the 4 sea-ice velocity components
      !!
      REAL(wp) :: zA, zUi, zVi, zUo, zVo, zM, zmsk
      REAL(wp), DIMENSION(4) :: zRHS
      REAL(wp) :: zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f
      REAL(wp) :: zms0x, zmsax, zmsbx, zms0y, zmsay, zmsby, zrhoco, zAu, zAv, zMu, zMv, zcorio_u, zCorio_v
      REAL(wp), DIMENSION(4) :: zVn, zV0, zV_old
      !REAL(wp), DIMENSION(jpi,jpj,4) :: z3d !debug
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_uv_rk3')

      !! Divergence of the 2 stress tensors
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !%acc parallel loop collapse(2) present(ff_u,ff_v,pAu,pAv,pmu_dt,pmv_dt,pSt,pSf,pgrdSH,pV_oce,pTau_ai,pm01x,pm01y,pm00x,pm00y,pV)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !
#           include "icedyn_rhg_vel_divs_t.h90"
            !
#           include "icedyn_rhg_vel_divs_f.h90"
            !
            zrhoco = rho0*rn_Cd_io
            zAu = pAu(ji,jj)
            zAv = pAv(ji,jj)
            zMu = pmu_dt(ji,jj)*pdt
            zMv = pmv_dt(ji,jj)*pdt
            zCorio_u = ff_u(ji,jj)
            zCorio_v = ff_v(ji,jj)
            !
            zms0x = pm00x(ji,jj)
            zmsax = pm01x(ji,jj)
            zmsbx = 1._wp - zmsax
            zms0y = pm00y(ji,jj)
            zmsay = pm01y(ji,jj)
            zmsby = 1._wp - zmsay
            !
            zV_old(:) = pV(ji,jj,:)

            ! Now time for the RK3 stuff..
            !! # 1
            zRHS(:) = FEEO_velocities( ridlzd, zrhoco, zAu, zAv, zMu, zMv, zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f, &
               &                      pgrdSH(ji,jj,:), pV_oce(ji,jj,:), pTau_ai(ji,jj,:), zV_old(:), zCorio_u, zCorio_v )
            zV0(:) = zV_old(:) + pdt * zRHS(:)
            !! # 2
            zRHS(:) = FEEO_velocities( ridlzd, zrhoco, zAu, zAv, zMu, zMv, zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f, &
               &                      pgrdSH(ji,jj,:), pV_oce(ji,jj,:), pTau_ai(ji,jj,:),   zV0(:) , zCorio_u, zCorio_v )
            zVn(:) = 0.75_wp*zV_old(:)     +     0.25_wp*( zV0(:) + pdt * zRHS(:) )
            !! # 2
            zRHS(:) = FEEO_velocities( ridlzd, zrhoco, zAu, zAv, zMu, zMv, zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f, &
               &                      pgrdSH(ji,jj,:), pV_oce(ji,jj,:), pTau_ai(ji,jj,:),   zVn(:) , zCorio_u, zCorio_v )
            zV0(:) = 1._wp/3._wp*zV_old(:) + 2._wp/3._wp*( zVn(:) + pdt * zRHS(:) )

            !! Apply masking:
            pV(ji,jj,1) = ( zmsax * zV0(1)  +  zmsbx * pV_oce(ji,jj,1) * 0.01_wp ) * zms0x
            pV(ji,jj,2) = ( zmsay * zV0(2)  +  zmsby * pV_oce(ji,jj,2) * 0.01_wp ) * zms0y
            pV(ji,jj,3) = ( zmsay * zV0(3)  +  zmsby * pV_oce(ji,jj,3) * 0.01_wp ) * zms0y
            pV(ji,jj,4) = ( zmsax * zV0(4)  +  zmsbx * pV_oce(ji,jj,4) * 0.01_wp ) * zms0x
            !
         END DO
      END DO
      !%acc end parallel loop

      IF( ln_timing )   CALL timing_stop('update_uv_rk3')
   END SUBROUTINE update_uv_rk3



   FUNCTION FEEO_velocities( pidlzd, prhoco, pAu, pAv, pmu, pmv, pdivSx_t, pdivSy_t, pdivSx_f, pdivSy_f, pgrdSH, pV_oce, pTau_ai, pV, zcor_u, zcor_v )
      !!----------------------------------------------------------------------------------------------
      !%acc routine
      !!----------------------------------------------------------------------------------------------
      !!   Fully-Explicit Euler Operator for ice velocities update
      !!----------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(4)             :: FEEO_velocities
      REAL(wp),               INTENT(in) :: pidlzd      ! idealized (0.) case, or not (1.)
      REAL(wp),               INTENT(in) :: prhoco      ! \rho_oce * CD_oce
      REAL(wp),               INTENT(in) :: pAu, pAv    ! Ice concentration
      REAL(wp),               INTENT(in) :: pmu, pmv
      REAL(wp),               INTENT(in) :: pdivSx_t, pdivSy_t, pdivSx_f, pdivSy_f ! 4 comp. of the divergence of the stress tensor
      REAL(wp), DIMENSION(4), INTENT(in) :: pgrdSH   ! gradient of SSH
      REAL(wp), DIMENSION(4), INTENT(in) :: pV_oce    ! the 4 ocean velocity components
      REAL(wp), DIMENSION(4), INTENT(in) :: pTau_ai   ! the 4 air-ice wind stress components
      REAL(wp), DIMENSION(4), INTENT(in) :: pV       ! the 4 sea-ice velocity components
      REAL(wp),               INTENT(in) :: zcor_u, zcor_v
      !!
      REAL(wp) :: zAu, zAv, zUu, zVu, zUu_o, zVu_o, zUv, zVv, zUv_o, zVv_o, z1_Mu, z1_Mv
      REAL(wp) :: zTauO, zcorio, zt1, zt2, zRHS
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------
      !
      z1_Mu  = 1._wp / MAX(pmu,epsi20) ! 1/mass (ice+snow) per surface unit at u-points
      zUu    = pV(1)    ! U at u-points!
      zUu_o  = pV_oce(1)   ! U ocean current at u-points!
      zVu    = pV(4)   ! V at u-points!
      zVu_o  = pV_oce(4)   ! V ocean current at u-points!
      !
      z1_Mv  = 1._wp / MAX(pmv,epsi20) ! 1/mass (ice+snow) per surface unit at v-points
      zVv    = pV(2)    ! V at v-points!
      zVv_o  = pV_oce(2)   ! V ocean current at v-points!
      zUv    = pV(3)   ! V at u-points!
      zUv_o  = pV_oce(3)   ! U ocean current at v-points!

      !====================================q==========================================================
      zt1   = zUu_o - zUu
      zt2   = zVu_o - zVu
      zTauO = prhoco * SQRT( zt1*zt1 + zt2*zt2 ) * zt1
      !
      zCorio =   pmu * zcor_u * zVu ! Coriolis at Uu-points (no interpolation needed)
      !
      !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
      zRHS = pdivSx_t + pAu*( pTau_ai(1) + zTauO ) + pidlzd*(zCorio + pgrdSH(1)) ! [Pa] = [kg.m-1.s-2]
      !
      FEEO_velocities(1) = zRHS * z1_Mu
      !==============================================================================================

      !==============================================================================================
      zt1   = zVv_o - zVv
      zt2   = zUv_o - zUv
      zTauO = prhoco * SQRT( zt1*zt1 + zt2*zt2 ) * zt1
      !
      zCorio = - pmv * zcor_v * zUv ! Coriolis at Vv-points (no interpolation needed)
      !
      !                 !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
      zRHS = pdivSy_t + pAv*( pTau_ai(2) + zTauO ) + pidlzd*(zCorio + pgrdSH(2)) ! [Pa] = [kg.m-1.s-2]
      !
      FEEO_velocities(2) = zRHS * z1_Mv
      !==============================================================================================

      !==============================================================================================
      zt1   = zUv_o - zUv
      zt2   = zVv_o - zVv
      zTauO = prhoco * SQRT( zt1*zt1 + zt2*zt2 ) * zt1
      !
      zCorio =   pmv * zcor_v * zVv ! Coriolis at Uv-points (no interpolation needed)
      !
      !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
      zRHS = pdivSx_f + pAv*( pTau_ai(3) + zTauO ) + pidlzd*(zCorio + pgrdSH(3)) ! [Pa] = [kg.m-1.s-2]
      !
      FEEO_velocities(3) = zRHS * z1_Mv
      !==============================================================================================

      !==============================================================================================
      zt1   = zVu_o - zVu
      zt2   = zUu_o - zUu
      zTauO = prhoco * SQRT( zt1*zt1 + zt2*zt2 ) * zt1
      !
      zCorio = - pmu * zcor_u * zUu ! Coriolis at Vu-points (no interpolation needed)
      !
      !--- Sum of external forces (explicit solution) = F + tau_ia + Coriolis + spg + tau_io
      zRHS = pdivSy_f + pAu*( pTau_ai(4) + zTauO ) + pidlzd*(zCorio + pgrdSH(4)) ! [Pa] = [kg.m-1.s-2]
      !
      FEEO_velocities(4) = zRHS * z1_Mu
      !==============================================================================================

      !!
   END FUNCTION FEEO_velocities

   !!==============================================================================
END MODULE icedyn_rhg_vel
