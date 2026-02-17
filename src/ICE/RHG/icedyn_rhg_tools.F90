MODULE icedyn_rhg_tools
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_tools  ***
   !!   Sea-Ice dynamics : rheology Britle tools...
   !!======================================================================
   !! History :
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE par_ice
   USE ice            ! => taux_ai_v, tauy_ai_u are there
   USE lib_mpp,        ONLY: ctl_stop
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE remap_classic,  ONLY: rmpF2T, rmpT2F
   USE in_out_manager, ONLY: ln_timing
   USE timing

   !USE icevar, ONLY : test4nan

   IMPLICIT NONE

   PRIVATE

   !INTERFACE strain_rate
   !   MODULE PROCEDURE strain_rate_Cgrid, strain_rate_Egrid
   !END INTERFACE strain_rate

   INTERFACE P_max_sclr
      MODULE PROCEDURE P_max_bbm_sclr, P_max_cos_sclr
   END INTERFACE P_max_sclr

   PUBLIC sigmaII_sclr
   PUBLIC sigmaII_full
   PUBLIC strain_rate_all
   PUBLIC strain_rate_dsd
   PUBLIC strain_rate_min
   PUBLIC vel_div_t
   PUBLIC vel_ten_t
   PUBLIC vel_shear_f
   PUBLIC vel_maxshr_t
   !PUBLIC vel_delta_t

   PUBLIC div_stress_tensor
   !
   PUBLIC mohr_coulomb_dmg
   PUBLIC mohr_coulomb_dmg_mp
   !
   PUBLIC cross_nudging_init
   PUBLIC apply_cn_trd
   !PUBLIC apply_cn_wn5s
   PUBLIC apply_cn_gpu
   !
   !PUBLIC d_crit
   PUBLIC Visco_sclr
   PUBLIC Lambda_sclr
   !
   PUBLIC P_max_sclr
   !
   PUBLIC P_tilde_sclr
   !
   PUBLIC mc_incrmt
   !
   PUBLIC Elast_diag
   PUBLIC Visco_diag
   PUBLIC Lambda_diag
   PUBLIC P_tilde_diag
   PUBLIC P_max_diag



   !REAL(wp), DIMENSION(:,:), ALLOCATABLE, PUBLIC, SAVE :: xtcoast, xfcoast ! to prevent doing cross-nudging at the coast!
   !REAL(wp), DIMENSION(:,:), ALLOCATABLE, PUBLIC, SAVE :: xCNt, xCNf   ! cross nudging coefficients (space-dependant)

   LOGICAL,  PUBLIC, SAVE :: l_CN        !: whether cross nudging is used ?
   !LOGICAL,  PUBLIC, SAVE :: l_CN_is_2d  !: whether cross nudging coefficient is a 2D array, not a scalar
   REAL(wp), PUBLIC, SAVE :: rCNC_eff    !: effective cross-nudging coefficient [-]
   !$acc declare create( l_CN, rCNC_eff )

   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
CONTAINS


   FUNCTION sigmaII_sclr( ps11, ps22, ps12 )
      !!------------------------------------------------------------------------------------
      !$acc routine
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor => same units as input stresses
      !!------------------------------------------------------------------------------------
      REAL(wp)             :: sigmaII_sclr
      REAL(wp), INTENT(in) :: ps11, ps22, ps12  ! sigma_11, sigma_22, sigma_12
      !!
      REAL(wp) :: ztmp
      !!------------------------------------------------------------------------------------
      ztmp  = 0.5_wp * (ps11 - ps22)
      sigmaII_sclr = SQRT( ztmp*ztmp + ps12*ps12 )
      !!
   END FUNCTION sigmaII_sclr

   SUBROUTINE sigmaII_full( pSt, pSf, pSII )
      !!------------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor => same units as input stresses
      !!------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(in)  :: pSt, pSf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(out) :: pSII
      !!
      REAL(wp) :: ztmp, zs
      INTEGER  :: ji, jj
      !!------------------------------------------------------------------------------------
      !$acc data present( pSt, pSf, pSII )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ztmp        = 0.5_wp * ( pSt(ji,jj,1) - pSt(ji,jj,2) )
            zs          = pSf(ji,jj,3)
            pSII(ji,jj) = SQRT( ztmp*ztmp + zs*zs )
         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE sigmaII_full




   SUBROUTINE mohr_coulomb_dmg( pdt, pxpCt, pxpCf, pSclH_t, pSclH_f, p1mdt, p1mdf, psgmt, psgmf )
      !!======================================================================
      REAL(wp),                       INTENT(in)    :: pdt            ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pSclH_t, pSclH_f
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt          ! vertically-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf          ! vertically-integrated F-centric stress tensor (mind that `s12` is @T)
      !!======================================================================
      REAL(wp) :: zs11t, zs22t, zs12t, zs11f, zs22f, zs12f, zE
      REAL(wp) :: zdx, zrr, zCohe, zNlim, zsqrtE, zTd, zsigI, zsigII
      INTEGER  :: ji, jj
      !!======================================================================
      IF( ln_timing ) CALL timing_start('mohr_coulomb_dmg')
      !$acc data present( pxpCt, pxpCf, pSclH_t, pSclH_f, p1mdt, p1mdf, psgmt, psgmf, res_grd_loc_t, res_grd_loc_f )

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            zs11t = psgmt(ji,jj,1)
            zs22t = psgmt(ji,jj,2)
            zs12f = psgmt(ji,jj,3)
            zs11f = psgmf(ji,jj,1)
            zs22f = psgmf(ji,jj,2)
            zs12t = psgmf(ji,jj,3)

            !! --- Mohr-Coulomb test and britle update if necessary ---
#           include "icedyn_rhg_tools_mc_t.h90"

#           include "icedyn_rhg_tools_mc_f.h90"
            !
            ! ==> pdmgt(:,:), psgmt(:,:,3) | pdmgf(:,:), psgmf(:,:,3)
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing ) CALL timing_stop('mohr_coulomb_dmg')
   END SUBROUTINE mohr_coulomb_dmg



   SUBROUTINE mohr_coulomb_dmg_mp( pdt, pxpCt, pxpCf, pSclH_t, pSclH_f, p1mdt, p1mdf, psgmt, psgmf )
      !!======================================================================
      !! Version using a common mid-point MC test (@ points located midway between T and F points)
      !!======================================================================
      REAL(wp),                       INTENT(in)    :: pdt            ! (small) time-step [s]
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pxpCt, pxpCf
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)    :: pSclH_t, pSclH_f
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(inout) :: p1mdt, p1mdf
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmt          ! vertically-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: psgmf          ! vertically-integrated F-centric stress tensor (mind that `s12` is @T)
      !!======================================================================
      REAL(wp), DIMENSION(jpi,jpj)   :: zTd_t, zCh_t, zNl_t, zSI_h_t, zSII_h_t
      REAL(wp), DIMENSION(jpi,jpj)   :: zTd_f, zCh_f, zNl_f, zSI_h_f, zSII_h_f
      REAL(wp), DIMENSION(jpi,jpj,4) :: zTd_mp, zCh_mp, zNl_mp, zSI_h_mp, zSII_h_mp !! Contains the 4 mid-points `T-F mean value` of the T-centric cell ji,jj
      !                                                                             !! NE=1, NW=2, SW=3, SE=4
      REAL(wp) :: zs11, zs22, zs12, zdx, zrr, zE, zsqrtE, zinc
      REAL(wp) :: zTd_ne, zTd_nw, zTd_sw, zTd_se, zinc_ne, zinc_nw, zinc_sw, zinc_se
      INTEGER  :: ji, jj
      !!======================================================================
      IF( ln_timing ) CALL timing_start('mohr_coulomb_dmg_mp')
      !$acc data create(zTd_t,zCh_t,zNl_t,zSI_h_t,zSII_h_t,zTd_f,zCh_f,zNl_f,zSI_h_f,zSII_h_f,zTd_mp,zCh_mp,zNl_mp,zSI_h_mp,zSII_h_mp) present(pxpCt,pxpCf,pSclH_t,pSclH_f,p1mdt,p1mdf,psgmt,psgmf)

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0+1
         DO ji=Nis0, Nie0+1
#           include "icedyn_rhg_tools_mp_bmc_t.h90"
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0
         DO ji=Nis0-1, Nie0
#           include "icedyn_rhg_tools_mp_bmc_f.h90"
         END DO
      END DO
      !$acc end parallel loop

      CALL build_T_F_mp_val( zSI_h_t,  zSI_h_f,  zSI_h_mp )
      CALL build_T_F_mp_val( zSII_h_t, zSII_h_f, zSII_h_mp )
      CALL build_T_F_mp_val( zTd_t,    zTd_f,    zTd_mp )
      CALL build_T_F_mp_val( zCh_t,    zCh_f,    zCh_mp )
      CALL build_T_F_mp_val( zNl_t,    zNl_f,    zNl_mp )

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
#           include "icedyn_rhg_tools_mp_mc_t.h90"
#           include "icedyn_rhg_tools_mp_mc_f.h90"
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data

      IF( ln_timing ) CALL timing_stop('mohr_coulomb_dmg_mp')
   END SUBROUTINE mohr_coulomb_dmg_mp


   SUBROUTINE build_T_F_mp_val( pXt, pXf, pXmp )
      !!-------------------------------------------------------------------
      !! Contains the 4 mid-T-F-point mean value for the T-centric cell ji,jj
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj),   INTENT(in)  :: pXt, pXf
      REAL(wp), DIMENSION(jpi,jpj,4), INTENT(out) :: pXmp
      !!-------------------------------------------------------------------
      REAL(wp) :: zxc
      INTEGER  :: ji, jj
      !!-------------------------------------------------------------------
      !LOLOfixme: decrease the halo use (check what is actually needed in `icedyn_rhg_tools_mp_mc_*.h90` !!!)
      !$acc parallel loop collapse(2) present(pXt, pXf, pXmp)
      DO jj=Njs0-(nn_hls-1), Nje0+nn_hls
         DO ji=Nis0-(nn_hls-1), Nie0+nn_hls
            zxc = pXt(ji,jj)
            pXmp(ji,jj,1) = 0.5_wp * ( zxc + pXf(ji  ,jj  ) )   !   ne
            pXmp(ji,jj,2) = 0.5_wp * ( zxc + pXf(ji-1,jj  ) )   !   nw
            pXmp(ji,jj,3) = 0.5_wp * ( zxc + pXf(ji-1,jj-1) )   !   sw
            pXmp(ji,jj,4) = 0.5_wp * ( zxc + pXf(ji  ,jj-1) )   !   se
         END DO
      END DO
      !$acc end parallel loop
      !
   END SUBROUTINE build_T_F_mp_val



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
      rCNC_eff = rn_crndg / REAL( nbbm, wp )
      !$acc update device(l_CN, rCNC_eff )
   END SUBROUTINE cross_nudging_init


   SUBROUTINE cross_nudging_trd( cgt, ps11x, ps22x, ps12x,   ps11, ps22, ps12 )
      !!
      CHARACTER(len=1),         INTENT(in)    :: cgt                  ! 'T' or 'F' -centric grid ?
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: ps11x, ps22x, ps12x  ! 3 stresses reference on grid [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: ps11,  ps22,  ps12   ! 3 stresses on grid to correct [Pa]
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zIs11x, zIs22x, zIs12x
      INTEGER :: i1,i2, j1,j2
      !
      i1=Nis0 ; i2=Nie0
      j1=Njs0 ; j2=Nje0
      !
      IF( cgt == 'T' ) THEN
         !
         zIs11x(:,:) = rmpF2T( ps11x )
         zIs22x(:,:) = rmpF2T( ps22x )
         zIs12x(:,:) = rmpT2F( ps12x )
         !
         ps11(i1:i2,j1:j2) = ps11(i1:i2,j1:j2) - rCNC_eff*( ps11(i1:i2,j1:j2) - zIs11x(i1:i2,j1:j2) )
         ps22(i1:i2,j1:j2) = ps22(i1:i2,j1:j2) - rCNC_eff*( ps22(i1:i2,j1:j2) - zIs22x(i1:i2,j1:j2) )
         ps12(i1:i2,j1:j2) = ps12(i1:i2,j1:j2) - rCNC_eff*( ps12(i1:i2,j1:j2) - zIs12x(i1:i2,j1:j2) )
         !
      ELSEIF( cgt == 'F' ) THEN
         !
         zIs11x(:,:) = rmpT2F( ps11x )
         zIs22x(:,:) = rmpT2F( ps22x )
         zIs12x(:,:) = rmpF2T( ps12x )
         !
         ps11(i1:i2,j1:j2) = ps11(i1:i2,j1:j2) - rCNC_eff*( ps11(i1:i2,j1:j2) - zIs11x(i1:i2,j1:j2) )
         ps22(i1:i2,j1:j2) = ps22(i1:i2,j1:j2) - rCNC_eff*( ps22(i1:i2,j1:j2) - zIs22x(i1:i2,j1:j2) )
         ps12(i1:i2,j1:j2) = ps12(i1:i2,j1:j2) - rCNC_eff*( ps12(i1:i2,j1:j2) - zIs12x(i1:i2,j1:j2) )
         !
      ELSE
         CALL ctl_stop( 'STOP', 'cross_nudging_trd() => wrong type of grid-point: '//cgt )
      ENDIF
      !!
   END SUBROUTINE cross_nudging_trd


   !SUBROUTINE cross_nudging_wn5s( cgt, ps11x, ps22x, ps12x,   ps11, ps22, ps12 )
   !   !!
   !   CHARACTER(len=1),         INTENT(in)    :: cgt                  ! 'T' or 'F' -centric grid ?
   !   REAL(wp), DIMENSION(jpi,jpj), INTENT(in)    :: ps11x, ps22x, ps12x  ! 3 stresses reference on grid [Pa]
   !   REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: ps11,  ps22,  ps12   ! 3 stresses on grid to correct [Pa]
   !   !!
   !   REAL(wp), DIMENSION(jpi,jpj) :: zIs11x, zIs22x, zIs12x
   !   INTEGER :: i1,i2, j1,j2
   !   !
   !   i1=Nis0 ; i2=Nie0
   !   j1=Njs0 ; j2=Nje0
   !   !
   !   IF( cgt == 'T' ) THEN
   !      !
   !      zIs11x(:,:) = rmpF2T_wn5s( ps11x )
   !      zIs22x(:,:) = rmpF2T_wn5s( ps22x )
   !      zIs12x(:,:) = rmpT2F_wn5s( ps12x )
   !      !
   !      ps11(i1:i2,j1:j2) = ps11(i1:i2,j1:j2) - rCNC_eff*( ps11(i1:i2,j1:j2) - zIs11x(i1:i2,j1:j2) )
   !      ps22(i1:i2,j1:j2) = ps22(i1:i2,j1:j2) - rCNC_eff*( ps22(i1:i2,j1:j2) - zIs22x(i1:i2,j1:j2) )
   !      ps12(i1:i2,j1:j2) = ps12(i1:i2,j1:j2) - rCNC_eff*( ps12(i1:i2,j1:j2) - zIs12x(i1:i2,j1:j2) )
   !      !
   !   ELSEIF( cgt == 'F' ) THEN
   !      !
   !      zIs11x(:,:) = rmpT2F_wn5s( ps11x )
   !      zIs22x(:,:) = rmpT2F_wn5s( ps22x )
   !      zIs12x(:,:) = rmpF2T_wn5s( ps12x )
   !      !
   !      ps11(i1:i2,j1:j2) = ps11(i1:i2,j1:j2) - rCNC_eff*( ps11(i1:i2,j1:j2) - zIs11x(i1:i2,j1:j2) )
   !      ps22(i1:i2,j1:j2) = ps22(i1:i2,j1:j2) - rCNC_eff*( ps22(i1:i2,j1:j2) - zIs22x(i1:i2,j1:j2) )
   !      ps12(i1:i2,j1:j2) = ps12(i1:i2,j1:j2) - rCNC_eff*( ps12(i1:i2,j1:j2) - zIs12x(i1:i2,j1:j2) )
   !      !
   !   ELSE
   !      CALL ctl_stop( 'STOP', 'cross_nudging_wn5s() => wrong type of grid-point: '//cgt )
   !   ENDIF
   !   !!
   !END SUBROUTINE cross_nudging_wn5s


   SUBROUTINE apply_cn_trd( kts, pS_t, pS_f )
      !!--------------------------------------------------------------------------------
      INTEGER,                        INTENT(in)    :: kts       ! current small time step
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_t      ! vertically-integrated T-centric stress tensor (mind that `s12` is @F)
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_f      ! vertically-integrated F-centric stress tensor (mind that `s12` is @T)
      !!--------------------------------------------------------------------------------
      !
      IF( MOD(kts,2) == 0 ) THEN
         !! Correction of T-centric stress tensor components:
         CALL cross_nudging_trd( 'T', pS_f(:,:,1), pS_f(:,:,2), pS_f(:,:,3),   pS_t(:,:,1), pS_t(:,:,2), pS_t(:,:,3) )
      ELSE
         !! Correction of F-centric stress tensor components:
         CALL cross_nudging_trd( 'F', pS_t(:,:,1), pS_t(:,:,2), pS_t(:,:,3),   pS_f(:,:,1), pS_f(:,:,2), pS_f(:,:,3) )
      END IF
      !
   END SUBROUTINE apply_cn_trd


   !SUBROUTINE apply_cn_wn5s( kts, pS_t, pS_f )
   !   !!--------------------------------------------------------------------------------
   !   INTEGER,                        INTENT(in)    :: kts       ! current small time step
   !   REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_t      ! vertically-integrated T-centric stress tensor (mind that `s12` is @F)
   !   REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_f      ! vertically-integrated F-centric stress tensor (mind that `s12` is @T)
   !   !!--------------------------------------------------------------------------------
   !   !
   !   IF( MOD(kts,2) == 0 ) THEN
   !      !! Correction of T-centric stress tensor components:
   !      CALL cross_nudging_wn5s( 'T', pS_f(:,:,1), pS_f(:,:,2), pS_f(:,:,3),   pS_t(:,:,1), pS_t(:,:,2), pS_t(:,:,3) )
   !   ELSE
   !      !! Correction of F-centric stress tensor components:
   !      CALL cross_nudging_wn5s( 'F', pS_t(:,:,1), pS_t(:,:,2), pS_t(:,:,3),   pS_f(:,:,1), pS_f(:,:,2), pS_f(:,:,3) )
   !   END IF
   !   !
   !END SUBROUTINE apply_cn_wn5s





   SUBROUTINE strain_rate_all( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                               pe11, pe22, pe12, pdudy, pdvdx, pdiv, pmshr, pdelta )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the 3 elements of the strain rate tensor, e11, e22 & e12, at either T- or F-points
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx` and `pmshr` !)
      !!----------------------------------------------------------------------------------------------------------
      CHARACTER(len=1),                   INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pe11, pe22, pe12          ! e11, e22 & e12 @ `cgt` points                   [1/s]
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pdudy, pdvdx, pdiv, pmshr  ! @ `cgt` points                [1/s]
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pdelta
      !!----------------------------------------------------------------------------------------------------------
      LOGICAL  :: l_r_e11, l_r_e22, l_r_e12, l_r_dudy, l_r_dvdx, l_r_div, l_r_mshr, l_r_dlt
      REAL(wp) :: zE1, zE2, zS1, zS2, zSHR, zdlt, z1_e1e2, zzf, ze2e2, ze1e1, zmask, zswitch, z1_ecc2
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!----------------------------------------------------------------------------------------------------------
      !!
      l_r_e11  = PRESENT(  pe11  )
      l_r_e22  = PRESENT(  pe22  )
      l_r_e12  = PRESENT(  pe12  )
      l_r_dudy = PRESENT(  pdudy )
      l_r_dvdx = PRESENT(  pdvdx )
      l_r_div  = PRESENT(  pdiv  )
      l_r_mshr = PRESENT(  pmshr )
      l_r_dlt  = PRESENT( pdelta )

      IF( l_r_dlt ) THEN
         PRINT *, 'LOLO [strain_rate_all@icedyn_rhg_tools.F90]: for `delta` => using ecc =', REAL(rn_delta_ecc)
         z1_ecc2 = 1._wp / ( rn_delta_ecc*rn_delta_ecc )
      ENDIF

      ! Prevent the occurence of NaN on the halos:
      IF( l_r_e11  )  pe11(:,:)   = 0._wp 
      IF( l_r_e22  )  pe22(:,:)   = 0._wp 
      IF( l_r_e12  )  pe12(:,:)   = 0._wp 
      IF( l_r_dudy )  pdudy(:,:)  = 0._wp 
      IF( l_r_dvdx )  pdvdx(:,:)  = 0._wp 
      IF( l_r_div  )  pdiv(:,:)   = 0._wp 
      IF( l_r_mshr )  pmshr(:,:)  = 0._wp 
      IF( l_r_dlt  )  pdelta(:,:) = 0._wp 

      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_all(): unknown grid-point type: '//cgt//'!')
      ENDIF

      !*acc data present(pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask) copyout(pe11, pe22, pe12, pdudy, pdvdx, pdiv, pmshr, pdelta)

      !*acc parallel loop collapse(2)
      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            IF( l_r_div .OR. l_r_e11 .OR. l_r_e22 .OR. l_r_mshr .OR. l_r_dlt ) THEN
               !! Divergence at cgt-points, `dU/dx + dV/dy` :
               zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
                  &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
                  &  ) * z1_e1e2
            ENDIF
            IF( l_r_e11 .OR. l_r_e22 .OR. l_r_mshr .OR. l_r_dlt  ) THEN
               !! Tension at cgt-points, `dU/dx - dV/dy` :
               zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
                  &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
                  &  ) * z1_e1e2
            ENDIF
            IF( l_r_e11 ) pe11(ji,jj) = 0.5_wp * ( zE1 + zE2 )
            IF( l_r_e22 ) pe22(ji,jj) = 0.5_wp * ( zE1 - zE2 )
            IF( l_r_div )  pdiv(ji,jj)  = zE1

            IF( l_r_e12 .OR. l_r_dudy .OR. l_r_dvdx .OR. l_r_mshr .OR. l_r_dlt ) THEN
               zzf = z1_e1e2 * zmask
               zS1 = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf    ! du/dy
               zS2 = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf    ! dv/dx
            ENDIF

            IF( l_r_e12 .OR. l_r_mshr .OR. l_r_dlt ) zSHR =  zS1 + zS2   ! shearing strain rate == 2*eps12 !

            IF( l_r_e12 ) pe12(ji,jj) = 0.5_wp * zSHR      ! pe12 == eps12 = 1/2 `shearing strain rate` !

            IF( l_r_dudy )  pdudy(ji,jj) = zS1
            IF( l_r_dvdx )  pdvdx(ji,jj) = zS2
            IF( l_r_mshr )  pmshr(ji,jj) = SQRT( zE2*zE2 + zSHR*zSHR )  ! Maximum shear: == SQRT( T^2 + S^2 )

            IF( l_r_dlt ) THEN
               zdlt = SQRT( zE1*zE1 + ( zE2*zE2 + zSHR*zSHR ) * z1_ecc2 ) ! `Delta` with z1_ecc2 => 0.25   ! Hunke & Dukowicz, 2002, Eq.5
               zswitch       = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zdlt ) ) ! 0 if delta=0
               pdelta(ji,jj) = zdlt + rn_creepl * zswitch
            ENDIF

         END DO
      END DO
      !*acc end parallel loop
      !*acc end data

   END SUBROUTINE strain_rate_all





   SUBROUTINE strain_rate_dsd( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                                 pdiv, pmshr, pdelta )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the divergence, maximum shear & delta
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx` and `pmshr` !)
      !!----------------------------------------------------------------------------------------------------------
      CHARACTER(len=1),                   INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:),           INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:),           INTENT(out) :: pdiv, pmshr, pdelta  ! @ `cgt` points                [1/s]
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp) :: zE1, zE2, zS1, zS2, zSHR, zdlt, z1_e1e2, zzf, ze2e2, ze1e1, zmask, zsw, z1_ecc2
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!----------------------------------------------------------------------------------------------------------
      !$acc data present( pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, pdiv, pmshr, pdelta )
      z1_ecc2 = 1._wp / ( rn_delta_ecc*rn_delta_ecc )

      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_dsd(): unknown grid-point type: '//cgt//'!')
      ENDIF

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pdiv(ji,jj) = 0._wp
            pmshr(ji,jj) = 0._wp
            pdelta(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      !$acc parallel loop collapse(2)
      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            !! Divergence at cgt-points, `dU/dx + dV/dy` :
            zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
               &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
               &  ) * z1_e1e2

            !! Tension at cgt-points, `dU/dx - dV/dy` :
            zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
               &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
               &  ) * z1_e1e2

            pdiv(ji,jj)  = zE1

            zzf = z1_e1e2 * zmask
            zS1 = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf    ! du/dy
            zS2 = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf    ! dv/dx
            zSHR =  zS1 + zS2   ! shearing strain rate == 2*eps12 !
            pmshr(ji,jj) = SQRT( zE2*zE2 + zSHR*zSHR )  ! Maximum shear: == SQRT( T^2 + S^2 )

            zdlt = SQRT( zE1*zE1 + ( zE2*zE2 + zSHR*zSHR ) * z1_ecc2 ) ! `Delta` with z1_ecc2 => 0.25   ! Hunke & Dukowicz, 2002, Eq.5
            zsw       = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zdlt ) ) ! 0 if delta=0
            pdelta(ji,jj) = zdlt + rn_creepl * zsw

         END DO
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE strain_rate_dsd




   SUBROUTINE strain_rate_min( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, pmask, &
      &                               pe11, pe22, pdiv, pdudy, pdvdx )
      !!
      !! Computes the 3 elements of the strain rate tensor, e11, e22 & e12, at either T- or F-points
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pe12`, `pdudy`, `pdvdx`!)
      !!
      CHARACTER(len=1),         INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pe11, pe22, pdiv
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdudy, pdvdx
      !!
      REAL(wp) :: zE1, zE2, z1_e1e2, zzf, ze2e2, ze1e1, zmask
      INTEGER  :: kq, ip, im, jp, jm, ji, jj, k1, k2
      !!
      kq = MAX( nn_hls-1, 0 )
      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  kq
         k2 =  1 + kq
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1 + kq
         k2 =  kq
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate_min(): unknown grid-point type: '//cgt//'!')
      ENDIF


      !*acc data present(pU,pV,pUd,pVd,p1_e1e2,pe2X,pe1Y,p1_e2X,p1_e1Y,pe1e1,pe2e2,pmask) copyout(pe11,pe22,pdiv,pdudy,pdvdx)

      !*acc parallel loop collapse(2)
      DO jj=Njs0-k1, Nje0+k2
         DO ji=Nis0-k1, Nie0+k2

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            !! Divergence at cgt-points, `dU/dx + dV/dy` :
            zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
               &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
               &  ) * z1_e1e2
            !! Tension at cgt-points, `dU/dx - dV/dy` :
            zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
               &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
               &  ) * z1_e1e2
            !!
            pe11(ji,jj) = 0.5_wp * ( zE1 + zE2 )
            pe22(ji,jj) = 0.5_wp * ( zE1 - zE2 )
            pdiv(ji,jj)  = zE1

            !! 2 * shear at cgt-points, `dU/dy + dV/dx` :
            zzf = z1_e1e2 * zmask
            pdudy(ji,jj) = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf
            pdvdx(ji,jj) = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf
            !
         END DO
      END DO
      !*acc end parallel loop
      !*acc end data

   END SUBROUTINE strain_rate_min








   SUBROUTINE div_stress_tensor( cgt, phc, phx, pe1e1, pe2e2,  pe1e1_e, pe2e2_e,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y,  &
      &                               ps11c, ps22c, ps12x,  pdivSx, pdivSy,  pm0 )
      !!----------------------------------------------------------------------------------------------
      !! Computes the vector (pdivSx,pdivSy) = divergence of the h-integrated internal stress tensor
      !!
      !!   depending on the grid: T-centric grid => cgt='T' or F-centric grid => cgt='F'
      !!
      !! INPUT:                                               |     cgt=='T'   |    cgt=='F'    |
      !!   * ps11c, ps22c: sigma11, sigma22           =>  ! @ point T[i,j] | @ point F[i,j] |
      !!   * ps12x       :       sigma12                =>  ! @ point F[i,j] | @ point T[i,j] |
      !!
      !! RETURNS:                                             |     cgt=='T'   |    cgt=='F'    |
      !!   * pdivSx: x-component of the div of the tensor =>  | @ point U[i,j] | @ point V[i,j] |
      !!   * pdivSy: y-component of the div of the tensor =>  | @ point V[i,j] | @ point U[i,j] |
      !!
      !!----------------------------------------------------------------------------------------------
      CHARACTER(len=1),         INTENT(in)  :: cgt
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: phc, phx   ! ice thickness at center and corner point [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pe1e1, pe2e2, pe1e1_e, pe2e2_e
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e1e2x, p1_e1e2y
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: ps11c, ps22c, ps12x ! components of stress tensors on T- or F-centric grids x h !!!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdivSx, pdivSy      ! x,y components of the divergence of the tensor
      INTEGER,        OPTIONAL, INTENT(in)  :: pm0
      !!
      INTEGER  :: ip, im, jp, jm, ji, jj, m0
      !!--------------------------------------------------------------------------------------------
      IF( ln_timing ) CALL timing_start('div_stress_tensor')
      m0 = 0
      IF( PRESENT(pm0) ) m0 = pm0

      IF ( cgt == 'T' ) THEN
         ip =  1
         im =  0
         jp =  0
         jm = -1
      ELSEIF ( cgt == 'F' ) THEN
         ip =  0
         im = -1
         jp =  1
         jm =  0
      ELSE
         CALL ctl_stop( 'STOP', 'div_stress_tensor(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      !pdivSx(:,:) = 0._wp
      !pdivSy(:,:) = 0._wp
      !
      !*acc data copyin(ip, im, jp, jm, phc, phx, ps11c, ps22c, ps12x, p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y)  copyout(pdivSx, pdivSy)
      !*acc parallel loop
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            !                   !--- ds11/dx + ds12/dy
            pdivSx(ji,jj) = ( ( ps11c(ji+ip,jj)*phc(ji+ip,jj)*pe2e2(ji+ip,jj)   - ps11c(ji+im,jj)*phc(ji+im,jj)*pe2e2(ji+im,jj)   ) * p1_e2x(ji,jj) &
               &            + ( ps12x(ji,jj+jp)*phx(ji,jj+jp)*pe1e1_e(ji,jj+jp) - ps12x(ji,jj+jm)*phx(ji,jj+jm)*pe1e1_e(ji,jj+jm) ) * p1_e1x(ji,jj) &
               &                 ) * p1_e1e2x(ji,jj)
            !                   !--- ds22/dy + ds12/dx
            pdivSy(ji,jj) = ( ( ps22c(ji,jj-jm)*phc(ji,jj-jm)*pe1e1(ji,jj-jm)   - ps22c(ji,jj-jp)*phc(ji,jj-jp)*pe1e1(ji,jj-jp)   ) * p1_e1y(ji,jj) &
               &            + ( ps12x(ji-im,jj)*phx(ji-im,jj)*pe2e2_e(ji-im,jj) - ps12x(ji-ip,jj)*phx(ji-ip,jj)*pe2e2_e(ji-ip,jj) ) * p1_e2y(ji,jj) &
               &                 ) * p1_e1e2y(ji,jj)
            !
         END DO
      END DO
      !*acc end parallel loop
      !*acc end data
      !
      IF( ln_timing )   CALL timing_stop('div_stress_tensor')
      !
   END SUBROUTINE div_stress_tensor




   SUBROUTINE apply_cn_gpu( kts, pS_t, pS_f )
      !!========================================================================================================================
      INTEGER,                        INTENT(in)    :: kts       ! current small time step
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_t      ! MIND that pS_t is the stress tensor with all the stress components @T,
      REAL(wp), DIMENSION(jpi,jpj,3), INTENT(inout) :: pS_f      !  it is not the same thing as the "T-centric" stress tensor (s12@F)!
      !
      REAL(wp) :: zms, zrc, zml, zr1, zr2, zr3, zr4, zs11x, zs22x, zs12x
      INTEGER  :: ji, jj, i2, i3, i4, j2, j3, j4, kp
      !!========================================================================================================================
      !$acc data present( pS_t, pS_f )
      IF( ln_timing )   CALL timing_start('apply_cn_gpu')
      kp = MIN( nn_hls-1 , 1 )

      IF( MOD(kts,2) == 0 ) THEN
         !! Correction of T-centric stress tensor components
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-kp, Nje0+kp
            DO ji=Nis0-kp, Nie0+kp
               !!
#              include "icedyn_rhg_tools_cn_t.h90"
               !!
            END DO
         END DO
         !$acc end parallel loop
      ELSE
         !!
         !! Correction of F-centric stress tensor components
         !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !$acc parallel loop collapse(2)
         DO jj=Njs0-kp, Nje0+kp
            DO ji=Nis0-kp, Nie0+kp
               !!
#              include "icedyn_rhg_tools_cn_f.h90"
               !!
            END DO
         END DO
         !$acc end parallel loop
      END IF
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('apply_cn_gpu')
   END SUBROUTINE apply_cn_gpu



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
   !   !*acc data copyin(rsqrt_nu_rhoi, epsi06, epsi20, pcohe, pNlim, pE, pdx) copyout(d_crit) present(pSGM)
   !   !*acc parallel loop collapse(2)
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
   !   !*acc end parallel loop
   !   !*acc end data
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
      !REAL(wp) :: zE, zeta
      !***************************************************************************************
      ! Viscosity [Pa.s]:
      !    *** MEB (Dansereau et al., 2016):
      !       * V = V0 * (1 - d)**a * exp[-C*(1-A)]  (viscosity)
      !    *** BBM (Olason et al. 2022) [Eq.10/Eq.9]:
      !       *    V = V0 * (1 - d)**a * exp[b*-C*(1-A)]    (with b=a in Olason et al. 2022)
      Visco_sclr = rn_eta0 * p1md**nn_alrlx * pexpC**nn_btrlx  ! viscosity [Pa.s]
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



   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !! `P_max` interface
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   FUNCTION P_max_bbm_sclr( pexpC, ph )
      !***************************************************************************************
      ! This function is to be used with vertically-integrated stresses in [N/m^2*m]
      !    => hence the `h**2.5` in place of the `h**1.5`
      !    => returns `-P_max*h`, [Eq.8] of Olason et al. 2022
      !
      !***************************************************************************************
      !*acc routine
      !***************************************************************************************
      REAL(wp)                       :: P_max_bbm_sclr
      REAL(wp),           INTENT(in) :: pexpC       ! `EXP[ rn_C0*(1 - pA) ) ]` with `rn_C0=-20`
      REAL(wp),           INTENT(in) :: ph          ! Ice thickness            [m]
      !***************************************************************************************
      !
      !P_max_bbm_sclr = -rn_P0 * ph**1.5_wp * pexpC  ! `-P_max` (for sigI<0)
      P_max_bbm_sclr = -rn_P0 * ph**2.5_wp * pexpC   ! `-P_max` (for sigI<0)  !#LOLOsh `Pmax` must be in [N/m^2*m]
      !
   END FUNCTION P_max_bbm_sclr

   FUNCTION P_max_cos_sclr( pexpC, ph, p1_SgmI, pSgmII )
      !***************************************************************************************
      ! This function is to be used with vertically-integrated stresses in [N/m^2*m]
      !    => hence the `h**2.5` in place of the `h**1.5`
      !    => returns `-P_max*h`, [Eq.8] of Olason et al. 2022
      !***************************************************************************************
      !*acc routine
      !***************************************************************************************
      REAL(wp)                       :: P_max_cos_sclr
      REAL(wp),           INTENT(in) :: pexpC       ! `EXP[ rn_C0*(1 - pA) ) ]` with `rn_C0=-20`
      REAL(wp),           INTENT(in) :: ph          ! Ice thickness            [m]
      REAL(wp),           INTENT(in) :: p1_SgmI     ! `1/sigma_I`              [m^2/N/m] vertically-integrated stress !
      REAL(wp),           INTENT(in) :: pSgmII      ! `sigma_II`               [N/m^2*m] vertically-integrated stress !
      !***************************************************************************************
      REAL(wp) :: zang
      !***************************************************************************************
      !
      zang  = ATAN( pSgmII * p1_SgmI )
      !
      !!P_max_cos_sclr = -rn_P0 * ph**1.5_wp * pexpC  * COS( zang )   ! `-P_max` (for sigI<0)
      !
      P_max_cos_sclr = -rn_P0 * ph**2.5_wp * pexpC  * COS( zang )   ! `-P_max` (for sigI<0)  !#LOLOsh `Pmax` must be in [N/m^2*m]
      !
      !P_max_cos_sclr = -rn_P0 * ph**2.5_wp          * COS( zang )   ! `-P_max` (for sigI<0)  !#LOLOsh `Pmax` must be in [N/m^2*m]
      !
   END FUNCTION P_max_cos_sclr

   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   FUNCTION P_tilde_sclr( pPmax, pSgmI, p1_SgmI )
      !***************************************************************************************
      ! Expects vertically-integrated stress !!!
      !    => returns `-P_max*h`, [Eq.8] of Olason et al. 2022
      !
      !***************************************************************************************
      !*acc routine
      !***************************************************************************************
      REAL(wp)                       :: P_tilde_sclr
      REAL(wp),           INTENT(in) :: pPmax        ! vertically-integrated `P_max`   [N/m^2*m]
      REAL(wp),           INTENT(in) :: pSgmI        ! vertically-integrated `sigma_I` [N/m^2*m]
      REAL(wp),           INTENT(in) :: p1_SgmI      ! `1/pSgmI`
      !***************************************************************************************
      REAL(wp) :: zc0
      !***************************************************************************************
      !
      zc0          = 0.5_wp + SIGN( 0.5_wp, -pSgmI-epsi20 ) ! => if sigI<-epsi20 => zc0=1 else: zc0=0
      P_tilde_sclr = -zc0 * MIN( pPmax*p1_SgmI , 1._wp )
      !
   END FUNCTION P_tilde_sclr



   FUNCTION mc_incrmt( pSI, pSII, pN, pC, pTd )
      !$acc routine
      !!---------------------------------------------------------------------
      REAL(wp) :: mc_incrmt
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pSI, pSII ! 1st and second invariant of vertically-integrated (or not) stress tensor [N/m^2] or [N/m^2*m]
      REAL(wp), INTENT(in) :: pN        ! `Nlim` [N/m^2] or [N/m^2*m]
      REAL(wp), INTENT(in) :: pC        ! Cohesion [N/m^2] or [N/m^2*m]
      REAL(wp), INTENT(in) :: pTd       ! Characteristic time of propagation of damage [s]
      !!----------------------------------------------------------------------
      REAL(wp) :: z1_zsigI, zMC, z1_zMC, zc0, zc1, zdcrit
      !!----------------------------------------------------------------------
      z1_zsigI = SIGN( 1._wp , pSI ) / MAX( ABS(pSI), epsi20 )
      zMC = pSII + rmuMC*pSI
      z1_zMC   = SIGN( 1._wp , zMC ) / MAX( ABS(zMC), epsi20 )
      zc0 = 0.5_wp + SIGN( 0.5_wp , pSI + pN )     ! => `zc0=0` if `pSI < -pN` (`pSI` is negative!)
      zdcrit = zc0 * pC * z1_zMC  +  (zc0-1._wp) * pN * z1_zsigI
      zc0 = 0.5_wp + SIGN( 0.5_wp , zdcrit-epsi20 )        ! => `zc0=1` if `dcrit>0`
      zc1 = 0.5_wp + SIGN( 0.5_wp , 1._wp-zdcrit-epsi20 )  ! => `zc1=1` if `dcrit<1`
      !
      mc_incrmt = zc0*zc1 * (1._wp - zdcrit) / pTd
      !
      ! Comprehensive version:
      !zdcrit = 9999._wp
      !IF( pSI < -pN ) THEN
      !   zdcrit = -pN / pSI
      !ELSEIF( ABS(zMC) > epsi10 ) THEN
      !   zdcrit = pC / zMC
      !ENDIF
      !mc_incrmt = 0._wp
      !IF( (zdcrit>0._wp).AND.(zdcrit<1._wp) ) mc_incrmt = (1._wp - zdcrit) / pTd
      !
   END FUNCTION mc_incrmt



   FUNCTION Elast_diag( pA, p1md )
      !!----------------------------------------------------------------------
      !! Returns the elasticity of (damaged) sea-ice [N/m^2]
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: Elast_diag       ! [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pA               ! sea-ice concentration [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: p1md             ! sea-ice damage [-]
      !!----------------------------------------------------------------------
      REAL(wp) :: zxpC, z1md
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      !*acc data pcopyin(pA, p1md) copyout(Elast_diag)
      !*acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            z1md = p1md(ji,jj)
            zxpC = EXP( rn_C0*(1._wp - pA(ji,jj)) )  ! `expC` [Eq.8]
            Elast_diag(ji,jj) = rn_E0 * z1md * zxpC                 !  `E = E0 * (1 - d) * exp[-C*(1-A)]`
         END DO
      END DO
      !*acc end parallel loop
      !*acc end data
   END FUNCTION Elast_diag


   FUNCTION Visco_diag( pA, p1md )
      !!----------------------------------------------------------------------
      !! Returns the viscosity of (damaged) sea-ice [Pa.s]
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: Visco_diag       ! [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pA               ! sea-ice concentration [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: p1md               ! sea-ice damage [-]
      !!----------------------------------------------------------------------
      REAL(wp) :: zxpC, z1md
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            z1md = p1md(ji,jj)
            zxpC = EXP( rn_C0*(1._wp - pA(ji,jj)) )
            Visco_diag(ji,jj) = Visco_sclr( zxpC, z1md )
         END DO
      END DO
   END FUNCTION Visco_diag


   FUNCTION Lambda_diag( pA, p1md, pdt )
      !!----------------------------------------------------------------------
      !! Returns the viscosity of (damaged) sea-ice [Pa.s]
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: Lambda_diag       ! [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pA               ! sea-ice concentration [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: p1md               ! sea-ice damage [-]
      REAL(wp),                     INTENT(in) :: pdt              ! small time step used [s]
      !!----------------------------------------------------------------------
      REAL(wp) :: zxpC, z1md, zE
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            z1md = p1md(ji,jj)
            zxpC = EXP( rn_C0*(1._wp - pA(ji,jj)) )
            zE   = rn_E0 * z1md * zxpC     ! elasticity [N/m^2]
            Lambda_diag(ji,jj) = Lambda_sclr( zxpC, z1md, zE, pdt )
         END DO
      END DO
   END FUNCTION Lambda_diag




   FUNCTION P_max_diag( pA, ph, ps11, ps22, ps12 )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: P_max_diag       ! [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pA               ! Ice concentration
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ph               ! Ice thickness
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ps11, ps22, ps12 ! Vertically-integrated stress tensor components at given point! [N/m^2*m]
      !!----------------------------------------------------------------------
      REAL(wp) :: zxpC, zsigI, zsigII, z1_zsigI
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            zxpC   = EXP( rn_C0*(1._wp - pA(ji,jj)) )  ! `expC` [Eq.8]
            zsigI  = 0.5_wp * ( ps11(ji,jj) + ps22(ji,jj) ) ! sigI: normal stress aka first invariant
            zsigII = sigmaII_sclr( ps11(ji,jj), ps22(ji,jj), ps12(ji,jj) )
            z1_zsigI = SIGN( 1._wp , zsigI ) / MAX( ABS(zsigI), epsi20 )   ! 1/SigI without the SigI=0 singularity...
            !
            P_max_diag(ji,jj) = -1._wp * P_max_sclr( zxpC, ph(ji,jj), z1_zsigI, zsigII ) !COS `P_max` with the right sign !
            !P_max_diag(ji,jj) = -1._wp * P_max_sclr( zxpC, ph(ji,jj) ) !BBM `P_max` with the right sign !
            !
         END DO
      END DO
   END FUNCTION P_max_diag

   FUNCTION P_tilde_diag( pA, ph, ps11, ps22, ps12 )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: P_tilde_diag     !   [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pA               ! Ice concentration [-]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ph               ! Ice thickness [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ps11, ps22, ps12 ! Vertically-integrated stress tensor components at given point! [N/m^2*m]
      !!----------------------------------------------------------------------
      REAL(wp) :: zxpC, zsigI, zsigII, zPmax, zc0, z1_zsigI
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            zxpC     = EXP( rn_C0*(1._wp - pA(ji,jj)) )  ! `expC` [Eq.8]
            zsigI    = 0.5_wp * ( ps11(ji,jj) + ps22(ji,jj) ) ! sigI: normal stress aka first invariant
            zsigII   = sigmaII_sclr( ps11(ji,jj), ps22(ji,jj), ps12(ji,jj) )
            z1_zsigI = SIGN( 1._wp , zsigI ) / MAX( ABS(zsigI), epsi20 )   ! 1/SigI without the SigI=0 singularity...
            zPmax    = P_max_sclr( zxpC, ph(ji,jj), z1_zsigI, zsigII ) !COS      ! `-P_max` (for sigI<0)
            !zPmax    = P_max_sclr( zxpC, ph(ji,jj) ) !BBM      ! `-P_max` (for sigI<0)
            zc0      = 0.5_wp + SIGN( 0.5_wp, -zsigI-epsi20 )           ! => if sigI<-epsi20 => zc0=1 else: zc0=0
            P_tilde_diag(ji,jj) = -1._wp * zc0 * MIN( zPmax*z1_zsigI , 1._wp )  ! => P~ with the right sign !
         END DO
      END DO
   END FUNCTION P_tilde_diag




   SUBROUTINE vel_div_t( pU, pV, p1_e1e2t, pe2u, pe1v, pmskt, pdivt )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) divergence of sea-ice velocity vector => T-point
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t         ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2u, pe1v       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt           ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdivt   ! divergence @ T points               [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Divergence at T-points, `dU/dx + dV/dy` :
            pdivt(ji,jj) = (   pe2u(ji,jj)*pU(ji,jj) - pe2u(ji-1,jj)*pU(ji-1,jj) &
               &             + pe1v(ji,jj)*pV(ji,jj) - pe1v(ji,jj-1)*pV(ji,jj-1) &
               &            )           * p1_e1e2t(ji,jj) * pmskt(ji,jj)

         END DO
      END DO

   END SUBROUTINE vel_div_t


   SUBROUTINE vel_ten_t( pU, pV, p1_e1e2t, p1_e2u, p1_e1v, pe1e1t, pe2e2t, pmskt, ptent )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) tension of sea-ice velocity vector => T-point
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t         ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2u, p1_e1v   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1t, pe2e2t
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt           ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: ptent   ! tension @ T points               [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Tension at T-points, `dU/dx - dV/dy` :
            ptent(ji,jj) = (   ( pU(ji,jj)*p1_e2u(ji,jj) - pU(ji-1,jj)*p1_e2u(ji-1,jj) ) * pe2e2t(ji,jj) &
               &             - ( pV(ji,jj)*p1_e1v(ji,jj) - pV(ji,jj-1)*p1_e1v(ji,jj-1) ) * pe1e1t(ji,jj) &
               &            )           * p1_e1e2t(ji,jj) * pmskt(ji,jj)

         END DO
      END DO

   END SUBROUTINE vel_ten_t


   SUBROUTINE vel_shear_f( pU, pV, p1_e1e2f, p1_e1u, p1_e2v, pe1e1f, pe2e2f, pmskf, pe12f )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) shear of sea-ice velocity vector => F-point
      !!
      !! Note: the mask must be the actual `fmask` that takes into account the slip/no-slip conditions
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV            ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2f         !     1/(e1f*e2f)               [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1u, p1_e2v
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1f, pe2e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskf         ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pe12f  !  e12 @ F points            [1/s]
      !LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      !!----------------------------------------------------------------------------------------------------------
      !LOGICAL  :: l_b_lnk
      REAL(wp) :: zzf, zS1, zS2
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------
      !IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            ! Shear at F points:
            zzf = p1_e1e2f(ji,jj) * pmskf(ji,jj)
            zS1 = ( u_ice(ji,jj+1) * p1_e1u(ji,jj+1) - u_ice(ji,jj) * p1_e1u(ji,jj) ) * pe1e1f(ji,jj) * zzf
            zS2 = ( v_ice(ji+1,jj) * p1_e2v(ji+1,jj) - v_ice(ji,jj) * p1_e2v(ji,jj) ) * pe2e2f(ji,jj) * zzf
            pe12f(ji,jj) = 0.5_wp * ( zS1 + zS2 )    ! eps12 =  1/2 `shearing strain rate` !
            !
         END DO
      END DO

   END SUBROUTINE vel_shear_f


   SUBROUTINE vel_maxshr_t( pU, pV, p1_e1e2t, p1_e1e2f, p1_e1u, p1_e2v, p1_e2u, p1_e1v, pe1e1t, pe2e2t, pe1e1f, pe2e2f, &
      &                     pe1e2f, pmskt, pmskf, pms )
      !!----------------------------------------------------------------------------------------------------------
      !! Computes the (strain-rate) tension of sea-ice velocity vector => T-point
      !!
      !! Note: the mask@F must be the actual `fmask` that takes into account the slip/no-slip conditions
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV           ! u@U & v@V                                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2t, p1_e1e2f       ! 1/(e1t*e2t)              [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1u, p1_e2v
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2u, p1_e1v   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1t, pe2e2t
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1f, pe2e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e2f
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmskt, pmskf         ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pms   ! maximum shear @ T points               [1/s]
      !!----------------------------------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zSHR
      REAL(wp) :: zzf, zten, zs2
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------------------------------

      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ! 2*Shear at F points:
            zzf = p1_e1e2f(ji,jj) * pmskf(ji,jj)
            zSHR(ji,jj) = ( ( pU(ji,jj+1) * p1_e1u(ji,jj+1) - pU(ji,jj) * p1_e1u(ji,jj) ) * pe1e1f(ji,jj) * zzf  &
               &         + ( pV(ji+1,jj) * p1_e2v(ji+1,jj) - pV(ji,jj) * p1_e2v(ji,jj) ) * pe2e2f(ji,jj) * zzf  &
               &         ) * p1_e1e2f(ji,jj) * pmskf(ji,jj)

         END DO
      END DO

      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! Tension at T-points, `dU/dx - dV/dy` :
            zten = (   ( pU(ji,jj)*p1_e2u(ji,jj) - pU(ji-1,jj)*p1_e2u(ji-1,jj) ) * pe2e2t(ji,jj) &
               &     - ( pV(ji,jj)*p1_e1v(ji,jj) - pV(ji,jj-1)*p1_e1v(ji,jj-1) ) * pe1e1t(ji,jj) &
               &    )           * p1_e1e2t(ji,jj)

            !! Shear**2 at T points (doc eq. A16)
            zs2 =  ( zSHR(ji,jj  ) * zSHR(ji,jj  ) * pe1e2f(ji,jj  ) + zSHR(ji-1,jj  ) * zSHR(ji-1,jj  ) * pe1e2f(ji-1,jj  )  &
               &   + zSHR(ji,jj-1) * zSHR(ji,jj-1) * pe1e2f(ji,jj-1) + zSHR(ji-1,jj-1) * zSHR(ji-1,jj-1) * pe1e2f(ji-1,jj-1)  &
               &   ) * 0.25_wp * p1_e1e2t(ji,jj)

            !! Maximum shear rate at T points
            pms(ji,jj) = SQRT( zten*zten + zs2 ) * pmskt(ji,jj)

         END DO
      END DO
   END SUBROUTINE vel_maxshr_t




   !!==============================================================================
END MODULE icedyn_rhg_tools
