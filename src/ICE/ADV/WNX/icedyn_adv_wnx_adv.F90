MODULE icedyn_adv_wnx_adv

   !!---------------------------------------------------------------------
   !! WENO order 5 or 7 generalized to orthogonal curvilinear coordinates
   !!---------------------------------------------------------------------

   USE par_kind, ONLY: wp
   USE par_oce,  ONLY: jpi, jpj, Nis0, Njs0, Nie0, Nje0
   USE lib_mpp,   ONLY: ctl_stop
   USE in_out_manager, ONLY: numout, lwp
   USE iom            ! I/O library
   !
   USE par_ice, ONLY : nn_WNx, kp_weno, epsi20
   USE ice,     ONLY : weno_lw_t_x, weno_lw_t_y, weno_ow_t_x, weno_ow_t_y, &
      &                weno_lw_f_x, weno_lw_f_y, weno_ow_f_x, weno_ow_f_y
   USE remap_weno, ONLY : weno5_ISx_G, weno7_ISx_G
   !
# if defined _OPENACC
   USE lbclnk_gpu
# else
   USE lbclnk         ! lateral boundary conditions (or mpp links)
# endif
   !
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC  wenoX_rk3

   PUBLIC  read_weno_w, print_linear_w, print_optimal_w

   REAL(wp), PARAMETER :: &
      &                      r13_12 = 13._wp/12._wp, &
      &                      r3_4   =  3._wp/4._wp,  &
      &                      r1_4   =  1._wp/4._wp,  &
      &                      r1_3   =  1._wp/3._wp,  &
      &                      r2_3   =  2._wp/3._wp
   !$acc declare create( r13_12, r3_4, r1_4, r1_3, r2_3 )


   ! Work array (that should remain once for all in the memory of the GPU)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)     ::   zfs1, zfs2, zfs3, zfs4
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)     ::   ztrk1, ztrk2, zoper

   !!----------------------------------------------------------------------
   !! NANUQ_beta
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS



#include "icedyn_adv_wn5_x.h90"

#include "icedyn_adv_wn5_y.h90"

#include "icedyn_adv_wn7_x.h90"

#include "icedyn_adv_wn7_y.h90"








   SUBROUTINE wenoX_rk3( kt, cgt, pdt, pe1e2, p1_e1e2, pu, pv, pmlbc,  pf,  psf,  lSmesh )
      !!----------------------------------------------------------------------
      !! Time discretization operator for WENOX
      !! Strong Stability Preserving RK3 (Jiang & Shu, 1996)
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in) :: kt
      CHARACTER(len=1)            , INTENT(in) :: cgt  ! grid-point location of tracer to advect ('T'/'F')
      REAL(wp),                     INTENT(in) :: pdt  ! advection time step [s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pe1e2   ! `e1e2` at the relevant point (T or F)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: p1_e1e2 ! `1/e1e2 * msk` at the relevant point (T or F)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pu   ! sea-ice x-velocity @ U-points (or V-points)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pv   ! sea-ice y-velocity @ V-points (or U-points)
      INTEGER(1), DIMENSION(jpi,jpj,nn_hls,4), INTENT(in) ::   pmlbc      ! masks for solid LBCs
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pf   ! field to advect
      REAL(wp),          OPTIONAL,  INTENT(in)    :: psf        ! Scale factor to scale the advected field
      LOGICAL,           OPTIONAL,  INTENT(in)    :: lSmesh     ! Advect `pf*dx*dy` rather than `pf`
      !
      LOGICAL  :: l_SF=.FALSE., l_Sm=.FALSE.
      REAL(wp) :: zsf, z1_sf
      INTEGER :: ji, jj
      !!----------------------------------------------------------------------
      CALL timing_start('wenoX_rk3')
      !$acc data present( pe1e2, p1_e1e2, pu, pv, pmlbc, pf, ztrk1, ztrk2, zoper )

      l_SF = ( PRESENT(psf) )
      IF( PRESENT(lSmesh) ) l_Sm = lSmesh

      ! `pf` must be fully `lbc_lnk`ed !!!

      IF(l_SF) THEN
         zsf = psf
         z1_sf = 1._wp / zsf
         !
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               pf(ji,jj) = zsf * pf(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      IF( l_Sm ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               pf(ji,jj)  = pf(ji,jj) * pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      ! # 1
      CALL wnx_spc_op( cgt, p1_e1e2, pu, pv, pmlbc, pf, zoper )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ztrk1(ji,jj) = pf(ji,jj)  - pdt * zoper(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      ! # 2
      CALL wnx_spc_op( cgt, p1_e1e2, pu, pv, pmlbc, ztrk1, zoper )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            ztrk2(ji,jj) = r3_4*pf(ji,jj) + r1_4 * ( ztrk1(ji,jj) - pdt * zoper(ji,jj) )
         END DO
      END DO
      !$acc end parallel loop

      ! # 3 !
      CALL wnx_spc_op( cgt, p1_e1e2, pu, pv, pmlbc, ztrk2, zoper )
      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            pf(ji,jj)  = r1_3*pf(ji,jj) + r2_3 * ( ztrk2(ji,jj) - pdt * zoper(ji,jj) )
         END DO
      END DO
      !$acc end parallel loop

      IF( l_Sm ) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               pf(ji,jj)  = pf(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp  ! needs to be in km
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      IF(l_SF) THEN
         !$acc parallel loop collapse(2)
         DO jj=Njs0-nn_hls, Nje0+nn_hls
            DO ji=Nis0-nn_hls, Nie0+nn_hls
               pf(ji,jj) = z1_sf * pf(ji,jj)
            END DO
         END DO
         !$acc end parallel loop
      ENDIF

      !$acc end data
      CALL timing_stop('wenoX_rk3')
      !
   END SUBROUTINE wenoX_rk3


   SUBROUTINE wnx_spc_op( cgt, p1_e1e2, pu, pv, pmlbc, pf, pop )
      !!----------------------------------------------------------------------
      !! Spatial discretization operator for WENOX
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in)  :: cgt         ! grid-point location of tracer to advect ('T'/'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: p1_e1e2     ! `1/e1e2 * msk` at the relevant point (T or F)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pu          ! sea-ice velocity @ U/V-points * dy
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pv          ! sea-ice velocity @ U/V-points * dx
      INTEGER(1), DIMENSION(jpi,jpj,nn_hls,4), INTENT(in) :: pmlbc  ! masks for solid LBCs
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pf          ! field to advect  @ T/F-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pop         ! field to advect  @ T/F-points
      !!----------------------------------------------------------------------
      REAL(wp)  :: zup, zvp
      CHARACTER(len=1) :: cgtu, cgtv
      INTEGER  :: ji, jj, jk
      !!----------------------------------------------------------------------
      CALL timing_start('wnx_spc_op')
      !$acc data present( p1_e1e2, pu, pv, pmlbc, pf, pop, zfs1, zfs2, zfs3, zfs4 )

      IF( cgt=='T' )  THEN
         cgtu = 'U' ; cgtv = 'V'
      ELSEIF( cgt=='F' )  THEN
         cgtu = 'V' ; cgtv = 'U'
      ENDIF

      ! `pf` & `pop` are present, because `created` in `wenoX_rk3` !

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            zfs1(ji,jj) = 0._wp
            zfs2(ji,jj) = 0._wp
            zfs3(ji,jj) = 0._wp
            zfs4(ji,jj) = 0._wp
         END DO
      END DO
      !$acc end parallel loop

      IF(     nn_WNx == 5 ) THEN
         STOP'LOLO: FIXME WENO5 F weights!!!!'
         ! Upwind (for U,V > 0):
         CALL wn5_x(  1, pmlbc, pf, zfs1 ) ! zf@T[i] => zfuwx@U[i] !!!
         CALL wn5_y(  1, pmlbc, pf, zfs2 ) ! zf@T[j] => zfuwy@V[j] !!!
         ! Downwind (for U,V < 0):
         CALL wn5_x( -1, pmlbc, pf, zfs3 ) ! zf@T[i] => zfuwx@U[i-1] !!!
         CALL wn5_y( -1, pmlbc, pf, zfs4 ) ! zf@T[j] => zfdwy@V[j-1] !!!
         !
      ELSEIF( nn_WNx == 7 ) THEN
         IF( cgt=='T' ) THEN
            ! Upwind (for U,V > 0):
            CALL wn7_x(  1, pmlbc, e1t, weno_lw_t_x, weno_ow_t_x, pf, zfs1 ) ! zf@T[i] => zfuwx@U[i] !!!
            CALL wn7_y(  1, pmlbc, e2t, weno_lw_t_y, weno_ow_t_y, pf, zfs2 ) ! zf@T[j] => zfuwy@V[j] !!!
            ! Downwind (for U,V < 0):
            CALL wn7_x( -1, pmlbc, e1t, weno_lw_t_x, weno_ow_t_x, pf, zfs3 ) ! zf@T[i] => zfuwx@U[i-1] !!!
            CALL wn7_y( -1, pmlbc, e2t, weno_lw_t_y, weno_ow_t_y, pf, zfs4 ) ! zf@T[j] => zfdwy@V[j-1] !!!
            !
         ELSEIF( cgt=='F' ) THEN
            ! Upwind (for U,V > 0):
            CALL wn7_x(  1, pmlbc, e1f, weno_lw_f_x, weno_ow_f_x, pf, zfs1 )
            CALL wn7_y(  1, pmlbc, e2f, weno_lw_f_y, weno_ow_f_y, pf, zfs2 )
            ! Downwind (for U,V < 0):
            CALL wn7_x( -1, pmlbc, e1f, weno_lw_f_x, weno_ow_f_x, pf, zfs3 )
            CALL wn7_y( -1, pmlbc, e2f, weno_lw_f_y, weno_ow_f_y, pf, zfs4 )
         ENDIF
      ENDIF

      !! Compute transport of `f` on velocity points
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            zfs1(ji,jj) = zfs1(ji,jj) * pu(ji,jj)   ! transport of `f` @U[i] when U>0
            zfs2(ji,jj) = zfs2(ji,jj) * pv(ji,jj)   ! transport of `f` @V[i] when V>0
            !
            zfs3(ji,jj) = zfs3(ji,jj) * pu(ji-1,jj) ! transport of `f` @U[i-1] when U<0
            zfs4(ji,jj) = zfs4(ji,jj) * pv(ji,jj-1) ! transport of `f` @V[j-1] when V<0
         END DO
      END DO
      !$acc end parallel loop

      !! Linking required here!
      !!  => because `zfs` used +-2 stencils, so can't do better than `Nis0:Nie0,Njs0:Nje0` for them.
#if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', zfs1, zfs2, zfs3, zfs4 )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', zfs1, zfs2, zfs3, zfs4 )
#else
      CALL lbc_lnk( 'wnx_spc_op', zfs1,cgtu,-1._wp, zfs2,cgtv,-1._wp, zfs3,cgtu,-1._wp, zfs4,cgtv,-1._wp )
#endif

      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0

            !! We need to know the sign of velocity components at center grid points (0.5 factor not needed):
            zup = 0.5_wp + SIGN( 0.5_wp, pu(ji,jj) + pu(ji-1,jj) )
            zvp = 0.5_wp + SIGN( 0.5_wp, pv(ji,jj) + pv(ji,jj-1) )
            !
            pop(ji,jj) = (      zup     * ( zfs1(ji  ,jj  ) - zfs1(ji-1,jj) )   &
               &          + (1._wp-zup) * ( zfs3(ji+1,jj  ) - zfs3(ji  ,jj) )   &
               &          +     zvp     * ( zfs2(ji  ,jj  ) - zfs2(ji,jj-1) )   &
               &          + (1._wp-zvp) * ( zfs4(ji  ,jj+1) - zfs4(ji,jj  ) ) ) &
               &                  * p1_e1e2(ji,jj) * 1.E6_wp  ! needs to be in km    ! => does either the /dx or /dy, because u=u*dy, v=v*dx (curvilinear coordinates)

         END DO
      END DO
      !$acc end parallel loop

      !! Linking required here!
      !!  => needed for the next round in RK3 (`pf` has to be complete!)
#if defined _OPENACC
      IF( l_Iperio ) CALL lbc_lnk_EW_gpu( 'icedyn_rhg_bbm', pop )
      IF( l_Jperio ) CALL lbc_lnk_NS_gpu( 'icedyn_rhg_bbm', pop )
#else
      CALL lbc_lnk( 'wnx_spc_op', pop,cgt,1._wp )
#endif

      !$acc end data
      CALL timing_stop('wnx_spc_op')
      !
   END SUBROUTINE wnx_spc_op


   SUBROUTINE read_weno_w( kpi, kpj, kpp, kwo, cnm_weno_wght, pLWtx, pOWtx, pLWty, pOWty, &
      &                                                       pLWfx, pOWfx, pLWfy, pOWfy  )
      !!================================================================================================
      !! For each point of the 2D grid we are going to read a [2p x w] matrix (`2*p` rows, `w` columns)
      !!
      !!  * `w` is the order of the cuurent WENO scheme, 5, 7, etc
      !!  * `p` is the order of each of the `p` `p`-point stencils used
      !!  *        `w = 2*p -1`
      !!
      !!  => WENO5: `w=5`, `p=3` ==> 5th order accuracy, combination of 3 3rd order 3-point stencils.
      !!  => WENO7: `w=7`, `p=4` ==> 7th order accuracy, combination of 4 4th order 5-point stencils.
      !!
      !!  Hence the array for linear weights have 4 dimmensions : [Ni x Nj x 2p x w ]
      !!
      !!   Double precision is strongly advised here !!!
      !!
      !!================================================================================================
      INTEGER,                                          INTENT(in)  :: kpi, kpj, kpp, kwo
      CHARACTER(len=*),                                 INTENT(in)  :: cnm_weno_wght
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo),           INTENT(out) :: pLWtx
      REAL(8),    DIMENSION(kpi,kpj,kpp),               INTENT(out) :: pOWtx
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo),           INTENT(out) :: pLWty
      REAL(8),    DIMENSION(kpi,kpj,kpp),               INTENT(out) :: pOWty
      !
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo), OPTIONAL, INTENT(out) :: pLWfx
      REAL(8),    DIMENSION(kpi,kpj,kpp),     OPTIONAL, INTENT(out) :: pOWfx
      REAL(8),    DIMENSION(kpi,kpj,kpp,kwo), OPTIONAL, INTENT(out) :: pLWfy
      REAL(8),    DIMENSION(kpi,kpj,kpp),     OPTIONAL, INTENT(out) :: pOWfy
      !!================================================================================================
      !!
      CHARACTER(len=13), PARAMETER :: crtn = 'read_weno_w'
      !!
      INTEGER :: inum, js
      CHARACTER(len=64) :: cvar
      CHARACTER(len=12) :: cvlw    ! ex: `w5_lw_t_01_x`
      CHARACTER(len=1)  :: cwo, cp
      LOGICAL           :: lReadFp ! read weights at F-points?
      !!================================================================================================
      lReadFp = ( PRESENT(pLWfx) .AND. PRESENT(pOWfx) .AND. PRESENT(pLWfy) .AND. PRESENT(pOWfy) )

      WRITE(cwo,'(i1)') kwo

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' * Will read linear and optimal weights for WENO'//cwo
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   => in file ', TRIM(cnm_weno_wght)
      ENDIF


      CALL iom_open( cnm_weno_wght, inum )

      cp = 't'

      pLWtx = 0._wp
      pOWtx = 0._wp
      pLWty = 0._wp
      pOWty = 0._wp

      cvar = 'weno'//cwo//'_ow_'//cp//'_x' !
      IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
      CALL iom_get ( inum, jpdom_global, TRIM(cvar), pOWtx(:,:,:) )

      cvar = 'weno'//cwo//'_ow_'//cp//'_y'
      IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
      CALL iom_get ( inum, jpdom_global, TRIM(cvar), pOWty(:,:,:) )

      DO js = 1, kwo

         WRITE(cvlw,'("w",i1,"_lw_",a1,"_",i2.2,"_x")') kwo, cp, js   ! Name of linear weight array to read
         IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
         CALL iom_get ( inum, jpdom_global, TRIM(cvlw), pLWtx(:,:,:,js) )

         WRITE(cvlw,'("w",i1,"_lw_",a1,"_",i2.2,"_y")') kwo, cp, js   ! Name of linear weight array to read
         IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
         CALL iom_get ( inum, jpdom_global, TRIM(cvlw), pLWty(:,:,:,js) )

      END DO

      IF( lReadFp ) THEN

         cp = 'f'

         pLWfx = 0._wp
         pOWfx = 0._wp
         pLWfy = 0._wp
         pOWfy = 0._wp

         cvar = 'weno'//cwo//'_ow_'//cp//'_x'
         IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
         CALL iom_get ( inum, jpdom_global, TRIM(cvar), pOWfx(:,:,:) )

         cvar = 'weno'//cwo//'_ow_'//cp//'_y' !
         IF(lwp) WRITE(numout,*) '     --- reading optimal w. 3D array ',TRIM(cvar)
         CALL iom_get ( inum, jpdom_global, TRIM(cvar), pOWfy(:,:,:) )

         DO js = 1, kwo

            WRITE(cvlw,'("w",i1,"_lw_",a1,"_",i2.2,"_x")') kwo, cp, js   ! Name of linear weight array to read
            IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw)
            CALL iom_get ( inum, jpdom_global, TRIM(cvlw), pLWfx(:,:,:,js) )

            WRITE(cvlw,'("w",i1,"_lw_",a1,"_",i2.2,"_y")') kwo, cp, js ! Name of linear weight array to read
            IF(lwp) WRITE(numout,*) '     --- reading linear w. 3D array ',TRIM(cvlw) !
            CALL iom_get ( inum, jpdom_global, TRIM(cvlw), pLWfy(:,:,:,js) )

         END DO

      ENDIF !IF( lReadFp )

      CALL iom_close( inum )

      IF(lwp) WRITE(numout,*) ' * Linear and Optimal weights for WENO'//cwo//' read !'
      IF(lwp) WRITE(numout,*) ''

   END SUBROUTINE read_weno_w


   SUBROUTINE print_linear_w( kWo, pLWt,  cinfo, fmult )
      !
      !   Shows a fancy & comprehensive representatio of the linear weights
      !
      INTEGER,                    INTENT(in) :: kWo        ! order of WENO scheme
      REAL(8), DIMENSION(:,:),    INTENT(in) :: pLWt        ! matrix of linear weights
      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cinfo ! info you want to add to the blabla
      REAL(8)         , OPTIONAL, INTENT(in) :: fmult      ! field multiplicator to express values as intergers
      !
      CHARACTER(len=64)  :: cstrl, cstrr, cstr
      CHARACTER(len=128) :: cbla=''
      INTEGER :: kp, n1, n2, jl, js
      LOGICAL :: lAsInt
      !
      IF(PRESENT(cinfo)) cbla = TRIM(cinfo)
      lAsInt = PRESENT(fmult)
      !
      kp = (kWo + 1)/2 ! order of each stencil
      !
      n1 = SIZE(pLWt,1)
      n2 = SIZE(pLWt,2)
      !
      IF( (n1/=2*kp).OR.(n2/=kWo) ) THEN
         WRITE(numout,'(" * ERROR [print_linear_w] * : `pLWt` has not the expected shape of ",i2,"x",i2," !")') n1,n2
         WRITE(numout,*)
         STOP
      ENDIF
      !
      cstrl  = '';  cstrr  = ''
      DO js=1,kp-1
         WRITE(cstr,'(" F[i-",i1,"]")') js
         cstrl = TRIM(cstr)//TRIM(cstrl)
         WRITE(cstr,'(" F[i+",i1,"]")') kp-js
         cstrr = TRIM(cstr)//TRIM(cstrr)
      END DO
      cstr = '  '//TRIM(cstrl)//'  F[i]'//TRIM(cstrr)
      !
      WRITE(numout,*)
      WRITE(numout,*) '  *** Matrix of the linear weights '//TRIM(cbla)//':'
      !               -1      5      2      0      0     (l =    0 )
      WRITE(numout,*) cstr
      DO jl = 1, 2*kp ! stencil bias x 2
         IF(lAsInt) THEN
            WRITE(numout,*) ( INT(REAL(pLWt(jl,js)*fmult,4),2), js=1,kWo ) , '  (l=', INT(MOD(jl-1,kp),1),') SUM=', REAL(SUM(pLWt(jl,:)),4)
         ELSE
            WRITE(numout,*) (     REAL(pLWt(jl,js)      ,4),    js=1,kWo ) , '  (l=', INT(MOD(jl-1,kp),1),') SUM=', REAL(SUM(pLWt(jl,:)),4)
         ENDIF
      END DO
      WRITE(numout,*)
      IF(lAsInt) WRITE(numout,'("      (this matrix was multiplied by ",i2,")")') INT(fmult)
      WRITE(numout,*)
      !
   END SUBROUTINE print_linear_w

   SUBROUTINE print_optimal_w( kWo, pOWt,  cinfo, fmult )
      !
      !   Shows a fancy & comprehensive representatio of the optimal weights
      !
      INTEGER,                    INTENT(in) :: kWo        ! order of WENO scheme
      REAL(8), DIMENSION(:),      INTENT(in) :: pOWt        ! matrix of linear weights
      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cinfo      ! info you want to add to the blabla
      REAL(8)         , OPTIONAL, INTENT(in) :: fmult      ! field multiplicator to express values as intergers
      !
      CHARACTER(len=64)  :: cstrl, cstrr, cstr
      CHARACTER(len=128) :: cbla=''
      INTEGER :: kp, n1, jl
      LOGICAL :: lAsInt
      !
      IF(PRESENT(cinfo)) cbla = TRIM(cinfo)
      lAsInt = PRESENT(fmult)
      !
      kp = (kWo + 1)/2 ! order of each stencil
      !
      n1 = SIZE(pOWt,1)
      !
      IF( n1/=2*kp ) THEN
         WRITE(numout,'(" * ERROR [print_optimal_w] * : `pOWt` has not the expected shape of ",i2," !")') n1
         WRITE(numout,*)
         STOP
      ENDIF
      !
      WRITE(numout,*)
      WRITE(numout,*) '  *** Matrix of the optimal weights '//TRIM(cbla)//':'
      !               -1      5      2      0      0     (l =    0 )
      !WRITE(numout,*) cstr

      DO jl = 1, kp ! stencil bias x 2
         IF(lAsInt) THEN
            !                            3      1
            IF(jl==1) WRITE(numout,*) '   i-1/2  i+1/2'
            WRITE(numout,*) INT(REAL(pOWt(jl)*fmult,4),2), INT(REAL(pOWt(jl+kp)*fmult,4),2), '    (l =', INT(MOD(jl-1,kp),1),')'
         ELSE
            IF(jl==1) WRITE(numout,*) '         i-1/2          i+1/2'
            WRITE(numout,*) '    ', REAL(pOWt(jl),4), REAL(pOWt(jl+kp),4) , '    (l =', INT(MOD(jl-1,kp),1),')'
         ENDIF
      END DO
      WRITE(numout,*) 'SUM=',REAL(SUM(pOWt(:kp)),4),REAL(SUM(pOWt(kp+1:)),4)
      WRITE(numout,*)
      IF(lAsInt) WRITE(numout,'("      (this matrix was multiplied by ",i2,")")') INT(fmult)
      WRITE(numout,*)
      !
   END SUBROUTINE print_optimal_w

END MODULE icedyn_adv_wnx_adv
