MODULE icedyn_rhg_bri
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_bri  ***
   !!   Sea-Ice dynamics : Brittle rheology misc. functions and routines
   !!======================================================================
   !! History : L. Brodeau, 2026
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE phycst,  ONLY: rsqrt_nu_rhoi, rmuMC
   USE dom_oce
   USE par_ice
   USE in_out_manager,   ONLY: nit000, ln_timing !, lwp
   USE timing
   USE lib_mpp, ONLY : ctl_stop

   IMPLICIT NONE

   PRIVATE

   PUBLIC MC_ud_d_s
   PUBLIC MC_ud_d_s_mwp_init
   PUBLIC MC_ud_d_s_mwp

   PUBLIC cross_nudging_init
   PUBLIC apply_CN

   LOGICAL,  PUBLIC, SAVE :: l_CN        !: whether cross nudging is used ?
   REAL(wp), PUBLIC, SAVE :: rCNC_eff    !: effective cross-nudging coefficient [-]
   !$acc declare create( l_CN, rCNC_eff )

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: udMCxt, udMCxf

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE MC_ud_d_s( pdt, pxpCt, pxpCf, pScHt, pScHf, p1mdt, p1mdf, psgmt, psgmf )
      !!======================================================================
      !! --- Mohr-Coulomb test and britle update of damage and stress tensors if necessary ---
      !!======================================================================
      !INTEGER ,                       INTENT(in)    :: kt   ! current advective time-step
      REAL(wp),                       INTENT(in)    :: pdt            ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pScHt, pScHf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt ! vert.-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf ! vert.-integrated F-centric stress tensor (mind that `s12` is @T)
      !!======================================================================
      REAL(wp) :: z1md, zs11, zs22, zs12
      REAL(wp) :: zrr, zinc, zCohe, zTd, zsigI, zsigII
      INTEGER  :: ji, jj
      !!======================================================================
      IF( ln_timing ) CALL timing_start('MC_ud_d_s')
      !$acc data present( pxpCt,pxpCf,pScHt,pScHf,p1mdt,p1mdf,psgmt,psgmf,res_grd_loc_t,res_grd_loc_f )
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
#           include "icedyn_rhg_bri_mc_t.h90"
#           include "icedyn_rhg_bri_mc_f.h90"
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
      IF( ln_timing ) CALL timing_stop('MC_ud_d_s')
   END SUBROUTINE MC_ud_d_s


   SUBROUTINE MC_ud_d_s_mwp( pdt, pxpCt, pxpCf, pScHt, pScHf, p1mdt, p1mdf, psgmt, psgmf )
      !!======================================================================
      !! Version using a common mid-point MC test (@ points located midway between T and F points)
      !!======================================================================
      !INTEGER ,                       INTENT(in)    :: kt   ! current advective time-step
      REAL(wp),                       INTENT(in)    :: pdt            ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pScHt, pScHf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt ! vert.-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf ! vert.-integrated F-centric stress tensor (mind that `s12` is @T)
      !!======================================================================
      REAL(wp) :: z1md, zs11, zs22, zs12, zdx, zrr, zE, zsqrtE, zinc
      REAL(wp) :: zTd_t, zI1_t, zI2_t, zSc_t
      REAL(wp) :: zTd_f, zI1_f, zI2_f, zSc_f
      REAL(wp) :: zinc_se0, zinc_ne0, zinc_nw0, zinc_sw0, zinc_sex, zinc_nex, zinc_nwx
      REAL(wp) :: zTd, zCh, zI1, zI2
      INTEGER  :: ji, jj, jit, jjt, jif, jjf
      !!======================================================================
      IF( ln_timing ) CALL timing_start('MC_ud_d_s_mwp')
      !$acc data present( pxpCt,pxpCf,pScHt,pScHf,p1mdt,p1mdf,psgmt,psgmf,udMCxt,udMCxf,res_grd_loc_t,res_grd_loc_f )

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !!------------------------------------------------------------------------------------------
            !! For T & F points of cell [ji,jj] we need to know everything at the following 7 points:
            !!
            !! pSE0 -> midway between T[i,j] & F[i  ,j-1]
            !! pNE0 -> midway between T[i,j] & F[i  ,j  ]
            !! pNW0 -> midway between T[i,j] & F[i-1,j  ]
            !! pSW0 -> midway between T[i,j] & F[i-1,j-1]
            !!
            !! pSEx -> midway between T[i+1,j  ] & F[i,j]
            !! pNEx -> midway between T[i+1,j+1] & F[i,j]
            !! pNWx -> midway between T[i  ,j+1] & F[i,j]
            !!        pSWx would be == pNE0
            !!
            !! For upcomming `ji` increment, solution at the new `pNW0` doesn't need to be recomputed
            !!   ==> it can inherit values of previous `pSEx` !
            !!------------------------------------------------------------------------------------------

            !! For pSE0, pNE0, pNW0, pSW0 we use the T-value remains the same
            zTd_t = TimeDamagePropag( res_grd_loc_t(ji,jj), p1mdt(ji,jj), pxpCt(ji,jj) )
            zSc_t = pScHt(ji,jj)
            zs11  = psgmt(ji,jj,1) ; zs22 = psgmt(ji,jj,2) ; zs12 = psgmf(ji,jj,3)
            zrr   = 0.5_wp * (zs11 - zs22)
            zI1_t = 0.5_wp * (zs11 + zs22)
            zI2_t = SQRT( zrr*zrr + zs12*zs12 )

            !! For pSEx, pNEx, pNWx       we use the F-value remains the same
            zTd_f = TimeDamagePropag( res_grd_loc_f(ji,jj), p1mdf(ji,jj), pxpCf(ji,jj) )
            zSc_f = pScHf(ji,jj)
            zs11  = psgmf(ji,jj,1) ; zs22 = psgmf(ji,jj,2) ; zs12 = psgmt(ji,jj,3)
            zrr   = 0.5_wp * (zs11 - zs22)
            zI1_f = 0.5_wp * (zs11 + zs22)
            zI2_f = SQRT( zrr*zrr + zs12*zs12 )

            !! Characteristic time, stress tensor invariants and cohesion @ pSE0:
            jif = ji ; jjf = jj-1
#           include "icedyn_rhg_bri_mcX_t.h90"
            zinc_se0 = mc_incrmt( zI1, zI2, zCh, zTd )

            !! Characteristic time, stress tensor invariants and cohesion @ pNE0:
            zTd = 0.5_wp*( zTd_t + zTd_f )
            zCh = 0.5_wp*( zSc_t + zSc_f ) * rn_c_ref
            zI1 = 0.5_wp*( zI1_t + zI1_f )
            zI2 = 0.5_wp*( zI2_t + zI2_f )
            zinc_ne0 = MC_incrmt( zI1, zI2, zCh, zTd )

            !! Characteristic time, stress tensor invariants and cohesion @ pNW0:
#if defined _OPENACC || defined _OPENMP
            jif = ji-1 ; jjf = jj
#           include "icedyn_rhg_bri_mcX_t.h90"
            zinc_nw0 = mc_incrmt( zI1, zI2, zCh, zTd )
#else
            IF( ji>Nis0 ) THEN
               !! The new `zinc_nw0` is the previous `zinc_sex` (valid only for cpu!)
               zinc_nw0 = zinc_sex
            ELSE
               jif = ji-1 ; jjf = jj
#              include "icedyn_rhg_bri_mcX_t.h90"
               zinc_nw0 = mc_incrmt( zI1, zI2, zCh, zTd )
            ENDIF
#endif

            !! Characteristic time, stress tensor invariants and cohesion @ pSW0:
            jif = ji-1 ; jjf = jj-1
#           include "icedyn_rhg_bri_mcX_t.h90"
            zinc_sw0 = mc_incrmt( zI1, zI2, zCh, zTd )

            !! Characteristic time, stress tensor invariants and cohesion @ pSEx:
            jit = ji+1 ; jjt = jj
#           include "icedyn_rhg_bri_mcX_f.h90"
            zinc_sex = mc_incrmt( zI1, zI2, zCh, zTd )

            !! Characteristic time, stress tensor invariants and cohesion @ pNEx:
            jit = ji+1 ; jjt = jj+1
#           include "icedyn_rhg_bri_mcX_f.h90"
            zinc_nex = mc_incrmt( zI1, zI2, zCh, zTd )

            !! Characteristic time, stress tensor invariants and cohesion @ pNWx:
            jit = ji ; jjt = jj+1
#           include "icedyn_rhg_bri_mcX_f.h90"
            zinc_nwx = mc_incrmt( zI1, zI2, zCh, zTd )

            !! Increment at both points:
            udMCxt(ji,jj) = 0.25_wp * ( zinc_se0 + zinc_ne0 + zinc_nw0 + zinc_sw0 ) * pdt
            udMCxf(ji,jj) = 0.25_wp * ( zinc_sex + zinc_nex + zinc_nwx + zinc_ne0 ) * pdt

         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! MC adjustment at T-points
            zinc = udMCxt(ji,jj)
            !
            z1md = p1mdt(ji,jj)
            zs11 = psgmt(ji,jj,1) ; zs22 = psgmt(ji,jj,2) ; zs12 = psgmf(ji,jj,3)
            !
            p1mdt(ji,jj) = MIN( MAX( z1md - z1md * zinc , r_dmd_min ) , 1._wp )
            zs11         =           zs11 - zs11 * zinc
            zs22         =           zs22 - zs22 * zinc
            zs12         =           zs12 - zs12 * zinc
            !
            psgmt(ji,jj,1) = zs11 ; psgmt(ji,jj,2) = zs22 ; psgmf(ji,jj,3) = zs12

            !! MC adjustment at F-points
            zinc = udMCxf(ji,jj)
            !
            z1md = p1mdf(ji,jj)
            zs11 = psgmf(ji,jj,1) ; zs22 = psgmf(ji,jj,2) ; zs12 = psgmt(ji,jj,3)
            !
            p1mdf(ji,jj) = MIN( MAX( z1md - z1md * zinc , r_dmd_min ) , 1._wp )
            zs11         =           zs11 - zs11 * zinc
            zs22         =           zs22 - zs22 * zinc
            zs12         =           zs12 - zs12 * zinc
            !
            psgmf(ji,jj,1) = zs11 ; psgmf(ji,jj,2) = zs22 ; psgmt(ji,jj,3) = zs12

         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing ) CALL timing_stop('MC_ud_d_s_mwp')
   END SUBROUTINE MC_ud_d_s_mwp

   SUBROUTINE MC_ud_d_s_mwp_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE MC_ud_d_s_mwp_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for "mid-way-point"
      !!                MC test and update...
      !!-------------------------------------------------------------------
      INTEGER :: ierr
      !!-------------------------------------------------------------------
      ierr = 0
      ALLOCATE( udMCxt(jpi,jpj) , udMCxf(jpi,jpj) , STAT=ierr )
      udMCxt(:,:)=0._wp ; udMCxf(:,:)=0._wp
      !
#if defined _OPENACC || defined _OPENMP
      PRINT *, ' * info GPU: MC_ud_d_s_mwp_init() => adding arrays to memory'
      PRINT *, '            => udMCxt, udMCxf'
      !$acc enter data copyin( udMCxt, udMCxf )
      PRINT *, ''
#endif
      !
      IF( ierr/=0 ) CALL ctl_stop('STOP', 'MC_ud_d_s_mwp_init: unable to allocate `udMCx*` arrays')
      !
   END SUBROUTINE MC_ud_d_s_mwp_init

   FUNCTION TimeDamagePropag( pdx, p1md, pxpC )
      !!---------------------------------------------------------------------
      !! Characteristic time for the propagation of damage / elastic shear waves
      !!---------------------------------------------------------------------
      !$acc routine seq
      !!---------------------------------------------------------------------
      REAL(wp) :: TimeDamagePropag  !: [s]
      !!----------------------------------------------------------------------
      REAL(4),  INTENT(in) :: pdx  ! local resolution (dx) of the grid [m]
      REAL(wp), INTENT(in) :: p1md      ! 1-damage
      REAL(wp), INTENT(in) :: pxpC      ! Hibler's exponential function / A
      !!----------------------------------------------------------------------
      REAL(wp) :: zE, zsqrtE, zdx
      !!----------------------------------------------------------------------
      zE     = rn_E0 * p1md * pxpC   ! Elasticity of ice
      zsqrtE = SQRT(MAX(zE,epsi10))
      zdx    = REAL( pdx , wp )
      !
      TimeDamagePropag = MAX( zdx * rsqrt_nu_rhoi / zsqrtE , epsi10 )
      !
   END FUNCTION TimeDamagePropag



   SUBROUTINE cross_nudging_init( )
      !!-------------------------------------------------------------------
      !! Called into `ice_dyn_rhg_init()@icedyn_rhg.F90`
      !!-------------------------------------------------------------------
      INTEGER  :: ierror
      REAL(wp) :: zr
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zt1, zt2, zt3, zt4
      INTEGER :: jm
      !!-------------------------------------------------------------------
      l_CN = ( rn_crndg > 0._wp )
      rCNC_eff = rn_crndg / REAL( nbrttl, wp )
      !$acc update device(l_CN, rCNC_eff )
   END SUBROUTINE cross_nudging_init

   SUBROUTINE apply_CN( kts, pS_t, pS_f )
      !!========================================================================================================================
      INTEGER,                        INTENT(in)    :: kts  ! current small time step
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_t ! vert.-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_f ! vert.-integrated F-centric stress tensor (mind that `s12` is @T)
      !!========================================================================================================================
      REAL(wp) :: zms, zrc, zml, zr1, zr2, zr3, zr4, zs11x, zs22x, zs12x
      INTEGER  :: ji, jj, i2, i3, i4, j2, j3, j4
      INTEGER  :: khep ! go `khep` points into the halo
      !!========================================================================================================================
      !$acc data present( pS_t, pS_f )
      IF( ln_timing )   CALL timing_start('apply_CN')

      khep = MERGE( 1 , 0 , ln_MCx_test )

      IF( MOD(kts,2) == 0 ) THEN
         !! Correction of T-centric stress tensor components
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-khep, Nje0+khep
            DO ji=Nis0-khep, Nie0+khep
#              include "icedyn_rhg_bri_cn_t.h90"
            END DO
         END DO
         !$acc end parallel loop
      ELSE
         !! Correction of F-centric stress tensor components
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-khep, Nje0+khep
            DO ji=Nis0-khep, Nie0+khep
#              include "icedyn_rhg_bri_cn_f.h90"
            END DO
         END DO
         !$acc end parallel loop
      END IF
      !$acc end data
      IF( ln_timing )   CALL timing_stop('apply_CN')
      !
   END SUBROUTINE apply_CN



   !FUNCTION d_crit( pcohe, pNlim, pE, pdx, pSGM )
   !   !!----------------------------------------------------------------------
   !   !! Fully Explicit Euler Operator for damage update
   !   !!----------------------------------------------------------------------
   !   REAL(wp), DIMENSION(jpi,jpj)                :: d_crit
   !   REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)  :: pcohe  ! cohesion
   !   REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)  :: pNlim  ! N
   !   REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)  :: pE     ! Elasticity of damaged ice
   !   REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)  :: pdx    ! Local grid resolution [m]
   !   REAL(wp), DIMENSION(jpi,jpj,3), INTENT(in)  :: pSGM   ! Stress tensor components
   !   !!
   !   REAL(wp) :: zsigI, zsigII, zMC
   !   REAL(wp) :: zsqrtE, zTd, zc0, z1_zsigI, z1_zMC, zNlim, ztmp
   !   INTEGER  :: ji, jj
   !   !!----------------------------------------------------------------------
   !   DO jj=Njs0, Nje0
   !      DO ji=Nis0, Nie0
   !
   !         zNlim = pNlim(ji,jj)
   !
   !         zsqrtE = SQRT(MAX(pE(ji,jj),epsi06))                               ! `sqrt(E)` (damaged ice)...
   !         zTd    = MAX( pdx(ji,jj) * rsqrt_nu_rhoi / zsqrtE , epsi06 )       ! characteristic time for damage [s] |  (we shall divide by it)...
   !
   !         zsigI  = 0.5_wp * (pSGM(ji,jj,1) + pSGM(ji,jj,2))
   !         ztmp   =           pSGM(ji,jj,1) - pSGM(ji,jj,2)
   !         zsigII = SQRT( 0.25_wp*ztmp*ztmp +  pSGM(ji,jj,3)*pSGM(ji,jj,3) )
   !
   !         z1_zsigI = SIGN( 1._wp , zsigI ) / MAX( ABS(zsigI), epsi20 )   ! 1/SigI without the SigI=0 singularity...
   !
   !         zMC = zsigII + rmuMC*zsigI                             ! Mohr-Coulomb  [Eq.29.2]
   !         z1_zMC = SIGN( 1._wp , zMC ) / MAX( ABS(zMC), epsi20 )   ! 1/MC without the MC=0 singularity...
   !
   !         zc0 = 0.5_wp + SIGN( 0.5_wp , zsigI + zNlim       )   ! if zsigI<-Nlim => zc0=0 ; zc0=1 otherwize
   !
   !         d_crit(ji,jj) = zc0 * pcohe(ji,jj) * z1_zMC  +  (zc0-1._wp) * zNlim * z1_zsigI   ! `zc0-1` because we need `-Nlim`
   !
   !      END DO
   !   ENDDO
   !END FUNCTION d_crit



   FUNCTION Visco_sclr( pexpC, p1md )
      !***************************************************************************************
      !   Returns `eta`, the viscosity of sea-ice [N/m^2.s]
      !***************************************************************************************
      !*acc routine
      !***************************************************************************************
      REAL(wp)                       :: Visco_sclr ! [s]
      REAL(wp),           INTENT(in) :: pexpC       ! `EXP[ rn_C0*(1 - pA) ) ]` with `rn_C0=-20`
      REAL(wp),           INTENT(in) :: p1md        ! `1-damage`   [:]
      !***************************************************************************************
      REAL(wp) :: zr1, zr2
      !***************************************************************************************
      ! Viscosity [Pa.s]:
      !    *** MEB (Dansereau et al., 2016):
      !       * V = V0 * (1 - d)**a * exp[  -C*(1-A)]  (viscosity)
      !    *** BBM (Olason et al. 2022) [Eq.10/Eq.9]:
      !       * V = V0 * (1 - d)**a * exp[b*-C*(1-A)]    (with b=a in Olason et al. 2022)
      zr1 = p1md*pexpC
      zr2 = zr1*zr1
      Visco_sclr = rn_eta0 * zr2 * zr2 * zr1  ! viscosity [Pa.s]
      !
   END FUNCTION Visco_sclr


   FUNCTION Lambda_sclr( pexpC, p1md, pE, pdt )
      !***************************************************************************************
      !   Returns `Lambda`, the " viscous relaxation time" [s]
      !***************************************************************************************
      !*acc routine
      !***************************************************************************************
      REAL(wp)                       :: Lambda_sclr ! [s]
      REAL(wp),           INTENT(in) :: pexpC       ! `EXP[ rn_C0*(1 - pA) ) ]` with `rn_C0=-20`
      REAL(wp),           INTENT(in) :: p1md        ! `1-damage`   [:]
      REAL(wp),           INTENT(in) :: pE          ! elasticity           [N/m^2]
      REAL(wp),           INTENT(in) :: pdt         ! small time step used [s]
      !***************************************************************************************
      REAL(wp) :: zeta
      !***************************************************************************************
      !
      zeta = Visco_sclr( pexpC, p1md ) ! viscosity [Pa.s]
      !
      Lambda_sclr = MAX( zeta / MAX( pE, epsi20 ) , pdt )
      !
   END FUNCTION Lambda_sclr



   FUNCTION MC_incrmt( pSI, pSII, pC, pTd )
      !!---------------------------------------------------------------------
      !! Mohr-Coulomb-test-based increment to update damage and internal stresses
      !!---------------------------------------------------------------------
      !$acc routine seq
      !!---------------------------------------------------------------------
      REAL(wp) :: MC_incrmt
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pSI, pSII ! 1st and second invariant of vert.-integrated (or not) stress tensor [N/m^2] or [N/m^2*m]
      REAL(wp), INTENT(in) :: pC        ! Cohesion [N/m^2] or [N/m^2*m]
      REAL(wp), INTENT(in) :: pTd       ! Characteristic time of propagation of damage [s]
      !!----------------------------------------------------------------------
      REAL(wp) :: zN, z1_zsigI, zMC, z1_zMC, zdcrit
      !!----------------------------------------------------------------------
      zN = pC * r_c_to_N    ! `zN == -Nlim` [N/m^2] or [N/m^2*m]
      z1_zsigI = SIGN( 1._wp , pSI ) / MAX( ABS(pSI), epsi10 )
      !
      zMC      = pSII + rmuMC*pSI
      z1_zMC   = SIGN( 1._wp , zMC ) / MAX( ABS(zMC), epsi10 )
      !
      zdcrit   = MERGE( zN * z1_zsigI  ,  pC * z1_zMC  ,  pSI < zN )
      !
      MC_incrmt = MERGE( (1._wp - zdcrit) / pTd  ,  0._wp  ,  (zdcrit>0._wp).AND.(zdcrit<1._wp) )
      ! Comprehensive version:
      !zdcrit = 9999._wp
      !IF( pSI < -zN ) THEN
      !   zdcrit = -zN / pSI
      !ELSEIF( ABS(zMC) > epsi10 ) THEN
      !   zdcrit = pC / zMC
      !ENDIF
      !MC_incrmt = 0._wp
      !IF( (zdcrit>0._wp).AND.(zdcrit<1._wp) ) MC_incrmt = (1._wp - zdcrit) / pTd
      !
   END FUNCTION MC_incrmt

   !!==============================================================================
END MODULE icedyn_rhg_bri
