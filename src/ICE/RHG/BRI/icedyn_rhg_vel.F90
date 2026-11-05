MODULE icedyn_rhg_vel
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_vel  ***
   !!
   !!   Sea-Ice dynamics : ice momentum equations and ice velocity components
   !!                      on the E-grid for Brittle Rheologies
   !!
   !!======================================================================
   !! History : Brodeau et al. 2024, GMD
   !!----------------------------------------------------------------------
   !!
   USE phycst, ONLY : rho0
   USE dom_oce        ! Ocean domain
   USE par_ice
   USE in_out_manager, ONLY : ln_timing

#if defined key_dbgvel
   USE ice, ONLY: divSx_t, divSy_t, divSx_f, divSy_f
#endif

   USE timing

   IMPLICIT NONE

   PRIVATE

   PUBLIC update_uv_eul
   PUBLIC update_uv_rk3

   !! Coefficients of the SDIRK3 method of Alexander (1977), used in `icedyn_rhg_vel_sdirk3.h90`:
   !!   => 3 stages, order 3, L-stable, stiffly accurate (last stage = new solution)
   !!   => Butcher tableau:
   !!            γ      |    γ       0      0
   !!         (1+γ)/2   | (1-γ)/2    γ      0
   !!            1      |    b1      b2     γ
   !!         ----------+-----------------------
   !!                   |    b1      b2     γ
   !!   => γ is the root of `γ^3 - 3γ^2 + 3/2γ - 1/6 = 0` that lies in ]1/6,1/2[:
   !!         γ = 1 - sqrt(2)*cos( π/3 + acos(2*sqrt(2)/3)/3 ) ≈ 0.435866521508459
   !!   => b1 = -(6γ^2 - 16γ + 1)/4 ≈  1.208496649
   !!      b2 =  (6γ^2 - 20γ + 5)/4 ≈ -0.644363171
   !!      and b1 + b2 = 1 - γ
   REAL(wp), PARAMETER :: rgam_sdk3 = 0.43586652150845899941601945_wp                                    ! γ (diagonal coefficient)
   REAL(wp), PARAMETER :: ra21_sdk3 = 0.5_wp*(1._wp - rgam_sdk3)                                         ! a21 = (1-γ)/2
   REAL(wp), PARAMETER :: rc2_sdk3  = 0.5_wp*(1._wp + rgam_sdk3)                                         ! c2  = (1+γ)/2 = a21 + γ
   REAL(wp), PARAMETER :: rb1_sdk3  = -0.25_wp*(6._wp*rgam_sdk3*rgam_sdk3 - 16._wp*rgam_sdk3 + 1._wp)    ! b1
   REAL(wp), PARAMETER :: rb2_sdk3  =  0.25_wp*(6._wp*rgam_sdk3*rgam_sdk3 - 20._wp*rgam_sdk3 + 5._wp)    ! b2

   !! Work array of `update_uv_rk3`: the 4 components of `A` (all the terms of RHS/mass except the ice-ocean stress)
   !!   => filled by the 1st loop of `update_uv_rk3` (stress divergence + RHS), read by the 2nd loop (SDIRK3 step)
   !!   => allocated at the first call of `update_uv_rk3`
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: xA_sdk3   ! [m.s-2]

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.1, Brodeau (2026)
   !! -- a standalone fork of NEMO/SI3 version 5 ---
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE update_uv_eul( kts, pdt, pAu, pAv, pmu_dt, pmv_dt, pSt, pSf, pgrdH, pVoce, ptaux_i_u, ptauy_i_v, ptaux_i_v, ptauy_i_u, &
      &                                     pm01x, pm01y, pm00x, pm00y, pucncl, pvcncl, pV )
      !!----------------------------------------------------------------------------------------------
      !!   Semi-Implicit Euler update operator for ice velocities update
      !!----------------------------------------------------------------------------------------------
      INTEGER,                          INTENT(in)    :: kts       ! current small time step
      REAL(wp),                         INTENT(in)    :: pdt       ! (small) time-step [s]
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pAu, pAv  ! Ice concentration @U and @V
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pmu_dt, pmv_dt ! Mass / dt @U and @V    [kg.m-2.s-1]
      REAL(wp),   DIMENSION(jpi,jpj,3), INTENT(in)    :: pSt, pSf  ! the 3 components the VERTICALLY-INTEGRATED stress tensor [N/m^2*m] !
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pgrdH    ! gradient of SSH * m
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pVoce    ! the 4 ocean velocity components
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: ptaux_i_u, ptauy_i_v, ptaux_i_v, ptauy_i_u ! air-ice windstress components at relevant points [N/m^2]
      INTEGER(1), DIMENSION(jpi,jpj),   INTENT(in)    :: pm01x, pm01y, pm00x, pm00y
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pucncl, pvcncl ! cancelers for the div. of the stress tensors where tiny ice concentration (@U,V points)
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(inout) :: pV        ! the 4 sea-ice velocity components
      !!----------------------------------------------------------------------------------------------
      REAL(wp) :: zS1, zS2, zS12
      REAL(wp) :: zA, zUi, zVi, zUo, zVo, zM_dt, zmsk, ztau_ai
      REAL(wp) :: zTauO, zcorio, zt1, zt2, zRHS
      REAL(wp) :: zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f
      REAL(wp) :: zmltU, zmltV
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_uv_eul')
      !$acc data present( ff_u,ff_v,pAu,pAv,pmu_dt,pmv_dt,pSt,pSf,pgrdH,pVoce,ptaux_i_u,ptauy_i_v,ptaux_i_v,ptauy_i_u,pm01x,pm01y,pm00x,pm00y,pV )

      IF( MOD(kts,2) == 0 ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0

               !! Divergence of the 2 vertically-integrated stress tensors
               !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! makes no sense to have huge gradient where the 2 involved points have a tiny concencentration of ice (A<~0.1)
               zmltU = r1_e1e2u(ji,jj) * pucncl(ji,jj)
               zmltV = r1_e1e2v(ji,jj) * pvcncl(ji,jj)

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
               ! makes no sense to have huge gradient where the 2 involved points have a tiny concencentration of ice (A<~0.1)
               zmltU = r1_e1e2u(ji,jj) * pucncl(ji,jj)
               zmltV = r1_e1e2v(ji,jj) * pvcncl(ji,jj)

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

      ENDIF !IF( MOD(kts,2) == 0 )

      !$acc end data
      IF( ln_timing )   CALL timing_stop('update_uv_eul')

   END SUBROUTINE update_uv_eul



   !! ************ !!
   !! RK3 IMPLICIT !!
   !! ************ !!

   SUBROUTINE update_uv_rk3( kts, pdt, pAu, pAv, p1_Mu, p1_Mv, pSt, pSf, pgrdH, pVoce,  &
      &                                     ptaux_i_u, ptauy_i_v, ptaux_i_v, ptauy_i_u, &
      &                                     pm01x, pm01y, pm00x, pm00y, pucncl, pvcncl, pVice )
      !!----------------------------------------------------------------------------------------------
      !!   Implicit RK3 method for ice velocities update
      !!     => SDIRK3 method of Alexander (1977), see `icedyn_rhg_vel_sdirk3.h90`
      !!
      !!   2 loops on purpose (for performance on CPU):
      !!     1. stress divergence + all the terms `A` of RHS/mass except the ice-ocean stress => `xA_sdk3`
      !!     2. SDIRK3 step + masking, with scalars only and no CALL, so that the CPU can work on several
      !!        grid points at the same time (the 3 stages form a long chain of dependent SQRT and divisions)
      !!----------------------------------------------------------------------------------------------
      INTEGER,                          INTENT(in)    :: kts       ! current small time step
      REAL(wp),                         INTENT(in)    :: pdt       ! (small) time-step [s]
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pAu, pAv  ! Ice concentration
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: p1_Mu, p1_Mv ! `1/mass` @U and @V                 [m2.kg-1]
      REAL(wp),   DIMENSION(jpi,jpj,3), INTENT(in)    :: pSt, pSf    ! the 3 components the stress tensor
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pgrdH    ! gradient of SSH
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(in)    :: pVoce    ! the 4 ocean velocity components
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: ptaux_i_u, ptauy_i_v, ptaux_i_v, ptauy_i_u ! air-ice windstress components at relevant points [N/m^2]
      INTEGER(1), DIMENSION(jpi,jpj),   INTENT(in)    :: pm01x, pm01y, pm00x, pm00y
      REAL(wp),   DIMENSION(jpi,jpj),   INTENT(in)    :: pucncl, pvcncl ! cancelers for the div. of the stress tensors where tiny ice concentration (@U,V points)
      REAL(wp),   DIMENSION(jpi,jpj,4), INTENT(inout) :: pVice          ! the 4 sea-ice velocity components
      !!----------------------------------------------------------------------------------------------
      REAL(wp) :: zS1, zS2, zS12
      REAL(wp) :: zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f
      REAL(wp) :: zms0x, zmsax, zmsbx, zms0y, zmsay, zmsby, zrhoco, zAu, zAv, z1_Mu, z1_Mv, zcorio_u, zCorio_v
      INTEGER  :: ji, jj, jp
      REAL(wp), DIMENSION(4) :: zTau_ai      ! air-ice windstress components at relevant points [N/m^2]
      REAL(wp), DIMENSION(4) :: zgrdH
      REAL(wp), DIMENSION(4) :: zVice        ! ice velocities at time k (for the Coriolis term)                                     [m.s-1]
      REAL(wp), DIMENSION(4) :: zA           ! all the terms of the RHS/mass except the ice-ocean stress (at time k)                [m.s-2]
      REAL(wp) :: zmltU, zmltV
      !! Scalars of the SDIRK3 step (loop 2 and `icedyn_rhg_vel_sdirk3.h90`):
      REAL(wp) :: zui, zvi, zuo, zvo, zAx, zAy, zBp, zun, zvn     ! inputs / outputs of the step for one 2D pair
      REAL(wp) :: zDx, zDy, zk, zSx, zSy, zSn, zfct               ! slip, `γ*Δt*B`, known part `S` of a stage, `|S|`, closed-form factor
      REAL(wp) :: zW1x, zW1y, zn1, zW2x, zW2y, zn2, zW3x, zW3y    ! `W = U_o - Y` of stages 1,2,3 and norms `|W|` of stages 1,2
      REAL(wp) :: zUu_n, zVu_n, zUv_n, zVv_n                      ! new velocities before masking: U-pair (1),(4) and V-pair (3),(2)
      !!----------------------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('update_uv_rk3')

      !! Allocate the work array `xA_sdk3` at the first call (kept for all the next calls)
      IF( .NOT. ALLOCATED(xA_sdk3) ) THEN
         ALLOCATE( xA_sdk3(jpi,jpj,4) )
         !$acc enter data create( xA_sdk3 )
      ENDIF

      !$acc data present( pAu,pAv,p1_Mu,p1_Mv,pSt,pSf,pgrdH,pVoce,ptaux_i_u,ptauy_i_v,ptaux_i_v,ptauy_i_u,pm01x,pm01y,pm00x,pm00y,pVice,xA_sdk3 ) create( zTau_ai,zgrdH,zVice,zA )

      !! `rho_o*CD_o`: constant, needed in loop 2 to compute `B = A_ice*rho_o*CD_o/mass`
      zrhoco = rho0*rn_Cd_io

      !! LOOP 1: divergence of the 2 stress tensors + the terms `A` of RHS/mass => `xA_sdk3`
      !! ~~~~~~
      !!   This loop reads neighbour points of the stress tensors but never reads `pVice` at other points,
      !!   so splitting it from loop 2 (which writes `pVice`) does not change the results.
      !$acc parallel loop collapse(2) private(zTau_ai,zgrdH,zVice,zA)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zAu = pAu(ji,jj)
            zAv = pAv(ji,jj)

            !! Divergence of the 2 vertically-integrated stress tensors
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! makes no sense to have huge gradient where the 2 involved points have a tiny concencentration of ice (A<~0.1)
            zmltU = r1_e1e2u(ji,jj) * pucncl(ji,jj)
            zmltV = r1_e1e2v(ji,jj) * pvcncl(ji,jj)

#           include "icedyn_rhg_vel_divs_t.h90"

#           include "icedyn_rhg_vel_divs_f.h90"

            !! RHS of the momentum equation
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            z1_Mu    = p1_Mu(ji,jj)
            z1_Mv    = p1_Mv(ji,jj)
            zCorio_u = ff_u(ji,jj)
            zCorio_v = ff_v(ji,jj)
            !$acc loop seq
            DO jp=1, 4
               ! Because of a `ld.lld: error: undefined symbol: _FortranAAssign` bug!
               zVice(jp) = pVice(ji,jj,jp)
               zgrdH(jp) = pgrdH(ji,jj,jp)
            END DO

            !BECAUSE `_FortranAAssign` bug (ROCm): zTau_ai(:) = (/ ptaux_i_u(ji,jj), ptauy_i_v(ji,jj), ptaux_i_v(ji,jj), ptauy_i_u(ji,jj) /)
            zTau_ai(1) = ptaux_i_u(ji,jj)
            zTau_ai(2) = ptauy_i_v(ji,jj)
            zTau_ai(3) = ptaux_i_v(ji,jj)
            zTau_ai(4) = ptauy_i_u(ji,jj)

            !! Term `A` of the momentum equation `dU/dt = A + B*|U_o - U|*(U_o - U)`:
            !!   all the terms of RHS/mass except the ice-ocean stress, at time k (Coriolis included)
            !!   (`B` is computed directly in loop 2, cheaper than storing it)
            CALL RHS_vel_sdirk3( ridlzd, zAu, zAv, z1_Mu, z1_Mv, zdivSx_t, zdivSy_t, zdivSx_f, zdivSy_f, &
               &                 zgrdH, zTau_ai, zVice, zCorio_u, zCorio_v,  zA )

            !! Store `A` for loop 2
            !$acc loop seq
            DO jp=1, 4
               xA_sdk3(ji,jj,jp) = zA(jp)
            END DO

         END DO
      END DO
      !$acc end parallel loop

      !! LOOP 2: SDIRK3 step for the 2 velocity pairs + masking
      !! ~~~~~~
      !!   Scalars only and no CALL (the SDIRK3 step is in `icedyn_rhg_vel_sdirk3.h90`), so that the CPU can
      !!   overlap the long chain of SQRT and divisions of several grid points, and the compiler can vectorize it.
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! U-points pair: (1) u@U and (4) v@U, coupled by the norm `|U_o - U|` at U-points
            !! `B = A_ice*rho_o*CD_o/mass` at U-points
            zui = pVice(ji,jj,1)
            zvi = pVice(ji,jj,4)
            zuo = pVoce(ji,jj,1)
            zvo = pVoce(ji,jj,4)
            zAx = xA_sdk3(ji,jj,1)
            zAy = xA_sdk3(ji,jj,4)
            zBp = zrhoco * pAu(ji,jj) * p1_Mu(ji,jj)
#           include "icedyn_rhg_vel_sdirk3.h90"
            zUu_n = zun
            zVu_n = zvn

            !! V-points pair: (3) u@V and (2) v@V, coupled by the norm `|U_o - U|` at V-points
            !! `B = A_ice*rho_o*CD_o/mass` at V-points
            zui = pVice(ji,jj,3)
            zvi = pVice(ji,jj,2)
            zuo = pVoce(ji,jj,3)
            zvo = pVoce(ji,jj,2)
            zAx = xA_sdk3(ji,jj,3)
            zAy = xA_sdk3(ji,jj,2)
            zBp = zrhoco * pAv(ji,jj) * p1_Mv(ji,jj)
#           include "icedyn_rhg_vel_sdirk3.h90"
            zUv_n = zun
            zVv_n = zvn

            !! Apply masking and store the new velocities (same masking as before)
            zms0x = pm00x(ji,jj)
            zmsax = pm01x(ji,jj)
            zmsbx = 1._wp - zmsax
            zms0y = pm00y(ji,jj)
            zmsay = pm01y(ji,jj)
            zmsby = 1._wp - zmsay
            pVice(ji,jj,1) = ( zmsax * zUu_n  +  zmsbx * pVoce(ji,jj,1) * 0.01_wp ) * zms0x
            pVice(ji,jj,2) = ( zmsay * zVv_n  +  zmsby * pVoce(ji,jj,2) * 0.01_wp ) * zms0y
            pVice(ji,jj,3) = ( zmsay * zUv_n  +  zmsby * pVoce(ji,jj,3) * 0.01_wp ) * zms0y
            pVice(ji,jj,4) = ( zmsax * zVu_n  +  zmsbx * pVoce(ji,jj,4) * 0.01_wp ) * zms0x

         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('update_uv_rk3')
   END SUBROUTINE update_uv_rk3


   SUBROUTINE RHS_vel_sdirk3( pidlzd, pAu, pAv, p1_Mu, p1_Mv, pdivSx_t, pdivSy_t, pdivSx_f, pdivSy_f, &
      &                       pgrdH, pTau_ai, pV, pcor_u, pcor_v,  pA )
      !!----------------------------------------------------------------------------------------------
      !$acc routine seq
      !!----------------------------------------------------------------------------------------------
      !!   Terms of the momentum equation, in the form used by the SDIRK3 step (`icedyn_rhg_vel_sdirk3.h90`):
      !!
      !!      `dU/dt = A + B*|U_o - U|*(U_o - U)`
      !!
      !!   * `A = 1/mass * [ div(h*sigma) + A_ice*Tau_air ] + Coriolis(U[k]) + (-g*grad[SSH])`   [m.s-2]
      !!        => all the terms except the ice-ocean stress
      !!        => the Coriolis term uses the ice velocities at time k, so `A` is constant in terms of `U[k+1]`
      !!   * `B = A_ice*rho_o*CD_o/mass`   [m-1]
      !!        => coefficient of the ice-ocean stress term, NOT computed here:
      !!           it is computed directly in loop 2 of `update_uv_rk3` (cheaper than storing it)
      !!
      !!   Only `A` is computed here. `mass` is in [kg.m-2], `A_ice` is the ice concentration.
      !!   The 4 components are: (1) u@U, (2) v@V, (3) u@V, (4) v@U
      !!----------------------------------------------------------------------------------------------
      REAL(wp),               INTENT(in)  :: pidlzd      ! idealized (0.) case, or not (1.)
      REAL(wp),               INTENT(in)  :: pAu, pAv    ! Ice concentration @U and @V
      REAL(wp),               INTENT(in)  :: p1_Mu, p1_Mv ! `1/mass` @U and @V                                [m2.kg-1]
      REAL(wp),               INTENT(in)  :: pdivSx_t, pdivSy_t, pdivSx_f, pdivSy_f ! 4 comp. of the divergence of the stress tensor
      REAL(wp), DIMENSION(4), INTENT(in)  :: pgrdH       ! `-g*grad[SSH]`                                     [m.s-2]
      REAL(wp), DIMENSION(4), INTENT(in)  :: pTau_ai     ! the 4 air-ice wind stress components               [N.m-2]
      REAL(wp), DIMENSION(4), INTENT(in)  :: pV          ! the 4 sea-ice velocity components at time k        [m.s-1]
      REAL(wp),               INTENT(in)  :: pcor_u, pcor_v ! Coriolis parameter @U and @V                    [s-1]
      REAL(wp), DIMENSION(4), INTENT(out) :: pA          ! all the terms except the ice-ocean stress          [m.s-2]
      !!----------------------------------------------------------------------------------------------
      REAL(wp) :: zUu, zVv, zUv, zVu
      !!----------------------------------------------------------------------------------------------

      !! Ice velocities at time k, needed for the Coriolis term
      zUu = pV(1)   ! U at U-points
      zVv = pV(2)   ! V at V-points
      zUv = pV(3)   ! U at V-points
      zVu = pV(4)   ! V at U-points

      !! (1) u@U: stress divergence + wind stress, then Coriolis (uses V@U, no interpolation needed) + SSH gradient
      pA(1) = p1_Mu * ( pdivSx_t + pAu*pTau_ai(1) )  +  pidlzd*(   pcor_u*zVu + pgrdH(1) )
      !! (2) v@V: stress divergence + wind stress, then Coriolis (uses U@V, no interpolation needed) + SSH gradient
      pA(2) = p1_Mv * ( pdivSy_t + pAv*pTau_ai(2) )  +  pidlzd*( - pcor_v*zUv + pgrdH(2) )
      !! (3) u@V: stress divergence + wind stress, then Coriolis (uses V@V, no interpolation needed) + SSH gradient
      pA(3) = p1_Mv * ( pdivSx_f + pAv*pTau_ai(3) )  +  pidlzd*(   pcor_v*zVv + pgrdH(3) )
      !! (4) v@U: stress divergence + wind stress, then Coriolis (uses U@U, no interpolation needed) + SSH gradient
      pA(4) = p1_Mu * ( pdivSy_f + pAu*pTau_ai(4) )  +  pidlzd*( - pcor_u*zUu + pgrdH(4) )
      !
   END SUBROUTINE RHS_vel_sdirk3




   !!==============================================================================
END MODULE icedyn_rhg_vel
