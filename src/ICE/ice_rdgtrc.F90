MODULE ice_rdgtrc
   !!======================================================================
   !!                       ***  MODULE  ice_rdgtrc  ***
   !!
   !!   Update the `ridge concentration` tracer  based on the change in
   !!   sea-ice volume after advection in regions of negative divergence
   !!   aka convergence.
   !!
   !!=====================================================================
   !! History :  NANUQ 1.0  !  2026-05  (L. Brodeau)       Split ice and ocean albedos
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_dyn_rdgtrc        : albedo for ice (clear and overcast skies)
   !!   ice_dyn_rdgtrc_init   : initialisation of albedo computation
   !!----------------------------------------------------------------------
   !USE phycst         ! physical constants
   USE dom_oce, ONLY: xmskt, xtmp1, xtmp2
   USE par_ice, ONLY: jpl, epsi10, epsi20
   !USE ice,
   !USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   !USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   !PUBLIC   ice_dyn_rdgtrc_init   ! called in icestp
   PUBLIC   ice_check_mean
   PUBLIC   ice_dyn_rdgtrc        ! called in icesbc.F90 and iceupdate.F90
   PUBLIC   ice_thd_rdgtrc

   REAL(wp), PARAMETER :: rA0 = 0.95_wp  ! minimum ice concentration required for ridging to occur...

   !! * Substitutions
   !#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 1.0.0, Brodeau (2026)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS



   SUBROUTINE ice_check_mean( kt, pa_i, pv_i, pat_i, pvt_i )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE ice_check_mean  ***
      !!
      !! ** Purpose :   Simply make sure that the current `pat_i` & `pvt_i` contain
      !!                what is expected
      !! References :
      !!
      !!----------------------------------------------------------------------
      INTEGER,                          INTENT(in) ::   kt       ! current time step
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   pa_i
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in) ::   pv_i
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pat_i
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in) ::   pvt_i
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl
      REAL(wp) ::   zAt, zVt
      INTEGER(1), DIMENSION(jpi,jpj) :: kfu
      !!---------------------------------------s------------------------------
      IF( ln_timing )   CALL timing_start('ice_check_mean')
      !$acc data present( pa_i, pv_i, pat_i, pvt_i, xtmp1, xtmp2 ) copyout( kfu )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls

            kfu(ji,jj) = 0

            zAt = 0._wp
            zVt = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zAt = zAt + pa_i(ji,jj,jl)
               zVt = zVt + pv_i(ji,jj,jl)
            END DO

            kfu(ji,jj) = kfu(ji,jj) + MERGE( 1 , 0 , ( pat_i(ji,jj) /= zAt ).OR.( pvt_i(ji,jj) /= zVt ) )

            !IF( pat_i(ji,jj) /= zAt ) CALL ctl_stop( 'ice_check_mean: mean `a_i` does not agree with `at_i` !')
            !IF( pvt_i(ji,jj) /= zVt ) CALL ctl_stop( 'ice_check_mean: mean `v_i` does not agree with `vt_i` !')

         END DO
      END DO
      !$acc end parallel loop

      !$acc end data

      IF( ANY( kfu > 0 ) ) CALL ctl_stop( 'ice_check_mean: mean `v_i` or `a_i` not consistent with `vt_i` or `at_i`, respectively !')

      IF( ln_timing )   CALL timing_stop('ice_check_mean')

   END SUBROUTINE ice_check_mean



   SUBROUTINE ice_dyn_rdgtrc( kt, pdiv, pa_i, pat_i_b, pv_i, pvt_i_b, prdgc )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn_rdgtrc  ***
      !!
      !! ** Purpose :   Update the ridged-ice concentration tracer after the advection.
      !!
      !! ** Method  :   based on Olason et al 2026 (neXtSIM v2)
      !!
      !!                If `R[k]` is the fraction of the tolal volume of ice at time `k`,
      !!                namely `V[k]`, that is ridged ice, then we have the following relation:
      !!
      !!                R[k+1] = 1 + (R[k]-1) * V[k]/V[k+1]
      !!
      !!
      !! ** Note    :
      !!
      !!
      !! References :
      !!
      !!
      !!----------------------------------------------------------------------
      INTEGER,                          INTENT(in   ) ::   kt       ! current time step
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pdiv     !  divergence (@ T-points) of sea-ice velocity [s^-1]
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in   ) ::   pa_i     !  "now" volume of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pat_i_b   !  "before" mean concentration of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in   ) ::   pv_i     !  "now" volume of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pvt_i_b   !  "before" mean volume of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(inout) ::   prdgc    !  ridge concentration tracer
      !
      INTEGER  ::  ji, jj, jl
      REAL(wp) ::  zAt_b, zAt_n, zVt_b, zVt_n, zdum   !, zu_t, zv_t
      LOGICAL  ::  lDo
      !!---------------------------------------s------------------------------
      IF( ln_timing )   CALL timing_start('ice_dyn_rdgtrc')
      !$acc data present( pdiv, pa_i, pat_i_b, pv_i, pvt_i_b, prdgc )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-(nn_hls-1), Nje0+(nn_hls-1)
         DO ji=Nis0-(nn_hls-1), Nie0+(nn_hls-1)

            prdgc(ji,jj) = MIN( prdgc(ji,jj) , 0.9999_wp ) ! correct potential overshoots created by the advection scheme used

            zAt_b = pat_i_b(ji,jj)
            zVt_b = pvt_i_b(ji,jj)

            !! Only bother if there is actually sea-ice & where advective velocity flow is convergent:
            IF( (zAt_b > rA0).AND.(pdiv(ji,jj) < -epsi20) ) THEN

               !zu_t = 0.5_wp*( pu(ji,jj) + pu(ji-1,jj) )  ! X-velocity @ T
               !zv_t = 0.5_wp*( pv(ji,jj) + pv(ji,jj-1) )  ! Y-velocity @ T

               !! Ridging can only occur when there is ice around or against a coastline:
               lDo =     ((pat_i_b(ji+1,jj)>rA0).OR.(xmskt(ji+1,jj)<0.1_wp)).AND.((pat_i_b(ji,jj+1)>rA0).OR.(xmskt(ji,jj+1)<0.1_wp)) &
                  & .AND.((pat_i_b(ji-1,jj)>rA0).OR.(xmskt(ji-1,jj)<0.1_wp)).AND.((pat_i_b(ji,jj-1)>rA0).OR.(xmskt(ji,jj-1)<0.1_wp))

               zAt_n = 0._wp
               zVt_n = 0._wp
               !$acc loop seq
               DO jl=1, jpl
                  zAt_n = zAt_n + pa_i(ji,jj,jl)
                  zVt_n = zVt_n + pv_i(ji,jj,jl)
               END DO

               zdum = 1._wp / MAX( zVt_n, epsi10 )   ! => 1/V[k+1]

               !  R[k+1]    =   1   + (   R[k]]     - 1    ) *              V[k]    / V[k+1]
               prdgc(ji,jj) = 1._wp + (prdgc(ji,jj) - 1._wp) * MERGE( MIN( ( zVt_b ) * zdum , 1._wp )  , 1._wp  , lDo ) ! Eq.(105) Olason et al., 2026

               prdgc(ji,jj) = MIN( prdgc(ji,jj) , 0.999_wp )

            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_dyn_rdgtrc')

   END SUBROUTINE ice_dyn_rdgtrc


   SUBROUTINE ice_thd_rdgtrc( kt, pv_i, pvt_i_b, prdgc )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE ice_thd_rdgtrc  ***
      !!
      !! ** Purpose :   Update the ridged-ice concentration tracer after thermodynamics
      !!
      !! ** Method  :   
      !!
      !! ** Note    :
      !!
      !! References :
      !!
      !!----------------------------------------------------------------------
      INTEGER,                          INTENT(in   ) ::   kt       ! current time step
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(in   ) ::   pv_i     !  "now" volume of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pvt_i_b  !  "before" mean volume of sea-ice [m]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(inout) ::   prdgc    !  ridge concentration tracer
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl
      REAL(wp) ::   zVt_b, zVt_n, z1_evol
      !!---------------------------------------s------------------------------
      IF( ln_timing )   CALL timing_start('ice_thd_rdgtrc')
      !$acc data present( pv_i, pvt_i_b, prdgc )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0-nn_hls, Nie0+nn_hls

            zVt_b = pvt_i_b(ji,jj)

            zVt_n = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zVt_n = zVt_n + pv_i(ji,jj,jl)
            END DO

            z1_evol = zVt_b / MAX( zVt_n, epsi20 )

            !                       Growth (zVt_n>zVt_b => z1_evol < 1) | Same melting for everyone |
            prdgc(ji,jj) = MERGE(        prdgc(ji,jj) * z1_evol         ,  prdgc(ji,jj)             ,  z1_evol < 1._wp  )

            !IF(  (z1_evol < 1._wp).AND.(zVt_b>0.5_wp)  ) THEN
            !   PRXNT *, 'LOLO: ridged ice volume decreased by: ', z1_evol, ji,jj
            !ENDIF
            
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_thd_rdgtrc')

   END SUBROUTINE ice_thd_rdgtrc







   !SUBROUTINE ice_dyn_rdgtrc_init
   !   !!----------------------------------------------------------------------
   !   !!                 ***  ROUTINE alb_init  ***
   !   !!
   !   !! ** Purpose :   initializations for the albedo parameters
   !   !!
   !   !! ** Method  :   Read the namelist namalb
   !   !!----------------------------------------------------------------------
   !   INTEGER ::   ios   ! Local integer output status for namelist read
   !   !!
   !   NAMELIST/namalb/ rn_alb_sdry, rn_alb_smlt, rn_alb_idry, rn_alb_imlt, rn_alb_dpnd, rn_alb_hpiv
   !   !!----------------------------------------------------------------------
   !   !
   !   READ_NML_REF(numnam_ice,namalb)
   !   !901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namalb in reference namelist' )
   !   READ_NML_CFG(numnam_ice,namalb)
   !   !902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namalb in configuration namelist' )
   !   IF(lwm) WRITE( numoni, namalb )
   !   !
   !   IF(lwp) THEN                      ! Control print
   !      WRITE(numout,*)
   !      WRITE(numout,*) 'ice_dyn_rdgtrc_init: set albedo parameters'
   !      WRITE(numout,*) '~~~~~~~~~~~~'
   !      WRITE(numout,*) '   Namelist namalb:'
   !      WRITE(numout,*) '      albedo of dry snow                   rn_alb_sdry = ', rn_alb_sdry
   !      WRITE(numout,*) '      albedo of melting snow               rn_alb_smlt = ', rn_alb_smlt
   !      WRITE(numout,*) '      albedo of dry ice                    rn_alb_idry = ', rn_alb_idry
   !      WRITE(numout,*) '      albedo of bare puddled ice           rn_alb_imlt = ', rn_alb_imlt
   !      WRITE(numout,*) '      albedo of ponded ice                 rn_alb_dpnd = ', rn_alb_dpnd
   !      WRITE(numout,*) '      pivotal ice thickness (m)            rn_alb_hpiv = ', rn_alb_hpiv
   !   ENDIF
   !   !
   !   !$acc update device( rn_alb_sdry, rn_alb_smlt, rn_alb_idry, rn_alb_imlt, rn_alb_dpnd, rn_alb_hpiv )
   !   !
   !END SUBROUTINE ice_dyn_rdgtrc_init

   !!======================================================================
END MODULE ice_rdgtrc
