MODULE icethd_da
   !!======================================================================
   !!                       ***  MODULE icethd_da ***
   !!   sea-ice : lateral melting
   !!======================================================================
   !! History :  3.7  !  2016-03  (C. Rousset)       Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!---------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd_da      : sea ice lateral melting
   !!   ice_thd_da_init : sea ice lateral melting initialization
   !!----------------------------------------------------------------------
   USE par_kind, ONLY : wp
   USE par_ice        ! SI3 parameters
   USE par_oce
   USE phycst  , ONLY : rpi, rt0, rhoi, rhos
   USE oss_nnq , ONLY : sst_s
   USE ice
   !
   USE in_out_manager , ONLY : numnam_ice_ref, numnam_ice_cfg, numout, numoni, lwp, lwm, ln_timing  ! I/O manager
   USE lib_mpp        , ONLY : ctl_stop, ctl_warn, ctl_nam                               ! MPP library
   !
   USE timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_da        ! called by icethd.F90
   PUBLIC   ice_thd_da_init   ! called by icestp.F90

   !                      !!** namelist (namthd_da) **
   REAL(wp) ::   rn_beta   ! coef. beta for lateral melting param.
   REAL(wp) ::   rn_dmin   ! minimum floe diameter for lateral melting param.
   !$acc declare create( rn_beta, rn_dmin )

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2025)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_da(jl_cat, ll_ice_present)
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_da  ***
      !!
      !! ** Purpose :   computes sea ice lateral melting
      !!
      !! ** Method  :   dA/dt = - P * W   [s-1]
      !!                   W = melting velocity [m.s-1]
      !!                   P = perimeter of ice-ocean lateral interface normalized by grid cell area [m.m-2]
      !!
      !!                   W = m1 * (Tw -Tf)**m2                    --- originally from Josberger 1979 ---
      !!                      (Tw - Tf) = elevation of water temp above freezing
      !!                      m1 and m2 = (1.6e-6 , 1.36) best fit from field experiment near the coast of Prince Patrick Island
      !!                                                                                           (Perovich 1983) => static ice
      !!                      m1 and m2 = (3.0e-6 , 1.36) best fit from MIZEX 84 experiment
      !!                                                                                (Maykut and Perovich 1987) => moving ice
      !!
      !!                   P = N * pi * D                           --- from Rothrock and Thorndike 1984 ---
      !!                      D = mean floe caliper diameter
      !!                      N = number of floes = ice area / floe area(average) = A / (Cs * D**2)
      !!                         A = ice concentration
      !!                         Cs = deviation from a square (square:Cs=1 ; circle:Cs=pi/4 ; floe:Cs=0.66)
      !!
      !!                   D = Dmin * ( Astar / (Astar-A) )**beta   --- from Lupkes et al., 2012 (eq. 26-27) ---
      !!
      !!                      Astar = 1 / ( 1 - (Dmin/Dmax)**(1/beta) )
      !!                      Dmin = minimum floe diameter (recommended to be 8m +- 20%)
      !!                      Dmax = maximum floe diameter (recommended to be 300m,
      !!                                                    but it does not impact melting much except for Dmax<100m)
      !!                      beta = 1.0 +-20% (recommended value)
      !!                           = 0.3 best fit for western Fram Strait and Antarctica
      !!                           = 1.4 best fit for eastern Fram Strait
      !!
      !! ** Tunable parameters  :   We propose to tune the lateral melting via 2 parameters
      !!                               Dmin [6-10m]   => 6  vs 8m = +40% melting at the peak (A~0.5)
      !!                                                 10 vs 8m = -20% melting
      !!                               beta [0.8-1.2] => decrease = more melt and melt peaks toward higher concentration
      !!                                                                  (A~0.5 for beta=1 ; A~0.8 for beta=0.2)
      !!                                                 0.3 = best fit for western Fram Strait and Antarctica
      !!                                                 1.4 = best fit for eastern Fram Strait
      !!
      !! ** Note   :   Former and more simple formulations for floe diameters can be found in Mai (1995),
      !!               Birnbaum and Lupkes (2002), Lupkes and Birnbaum (2005). They are reviewed in Lupkes et al 2012
      !!               A simpler implementation for CICE can be found in Bitz et al (2001) and Tsamados et al (2015)
      !!
      !! ** References
      !!    Bitz, C. M., Holland, M. M., Weaver, A. J., & Eby, M. (2001).
      !!              Simulating the ice‐thickness distribution in a coupled climate model.
      !!              Journal of Geophysical Research: Oceans, 106(C2), 2441-2463.
      !!    Josberger, E. G. (1979).
      !!              Laminar and turbulent boundary layers adjacent to melting vertical ice walls in salt water
      !!              (No. SCIENTIFIC-16). WASHINGTON UNIV SEATTLE DEPT OF ATMOSPHERIC SCIENCES.
      !!    Lüpkes, C., Gryanik, V. M., Hartmann, J., & Andreas, E. L. (2012).
      !!              A parametrization, based on sea ice morphology, of the neutral atmospheric drag coefficients
      !!              for weather prediction and climate models.
      !!              Journal of Geophysical Research: Atmospheres, 117(D13).
      !!    Maykut, G. A., & Perovich, D. K. (1987).
      !!              The role of shortwave radiation in the summer decay of a sea ice cover.
      !!              Journal of Geophysical Research: Oceans, 92(C7), 7032-7044.
      !!    Perovich, D. K. (1983).
      !!              On the summer decay of a sea ice cover. (Doctoral dissertation, University of Washington).
      !!    Rothrock, D. A., & Thorndike, A. S. (1984).
      !!              Measuring the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 89(C4), 6477-6486.
      !!    Tsamados, M., Feltham, D., Petty, A., Schroeder, D., & Flocco, D. (2015).
      !!              Processes controlling surface, bottom and lateral melt of Arctic sea ice in a state of the art sea ice model.
      !!              Phil. Trans. R. Soc. A, 373(2052), 20140167.
      !!---------------------------------------------------------------------
      INTEGER,                     INTENT(IN)    :: jl_cat
      LOGICAL, DIMENSION(jpi,jpj), INTENT(inout) :: ll_ice_present
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk ! dummy loop indices
      REAL(wp)            ::   zastar, zdfloe, zperi, zwlat, zda, zda_tot, zt1, zt2
      REAL(wp), PARAMETER ::   zdmax = 300._wp
      REAL(wp), PARAMETER ::   zcs   = 0.66_wp
      REAL(wp), PARAMETER ::   zm1   = 3.e-6_wp
      REAL(wp), PARAMETER ::   zm2   = 1.36_wp
      REAL(wp)            ::   zsum_s_i, zsum_e_i, zsum_e_s      ! SUM along `nlay_i` or `nlay_s`...
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('ice_thd_da')
      !$acc data present( a_i, at_i, e_i, e_s, hfx_thd, h_i, h_s, ll_ice_present, rDt_ice, at_i, sfx_lam, sst_s, t_bo, wfx_lam )
      
      zastar = 1._wp / ( 1._wp - (rn_dmin / zdmax)**(1._wp/rn_beta) )
      
      !$acc parallel loop collapse(2)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            IF(ll_ice_present(ji,jj)) THEN
               !
               zsum_s_i = 0._wp ; zsum_e_i = 0._wp               
               !$acc loop seq
               DO jk=1, nlay_i
                  zsum_e_i = zsum_e_i +  e_i(ji,jj,jk,jl_cat)
                  IF( nn_icesal == 4 ) THEN
                     zsum_s_i = zsum_s_i + sz_i(ji,jj,jk,jl_cat)  ! use layer salinity if nn_icesal=4
                  ELSE
                     zsum_s_i = zsum_s_i + s_i (ji,jj,  jl_cat)  !     bulk salinity otherwise (for conservation purpose)
                  ENDIF
               END DO
               !
               zsum_e_s = 0._wp               
               !$acc loop seq
               DO jk=1, nlay_s
                  zsum_e_s =    zsum_e_s +  e_s(ji,jj,jk,jl_cat)
               END DO
               !
               ! --- Calculate reduction of total sea ice concentration --- !
               zdfloe = rn_dmin * ( zastar / ( zastar - at_i(ji,jj) ) )**rn_beta           ! Mean floe caliper diameter [m]
               !
               zperi  = at_i(ji,jj) * rpi / ( zcs * zdfloe )                               ! Mean perimeter of the floe [m.m-2]
               !                                                                           !    = N*pi*D = (A/cs*D^2)*pi*D
               zwlat  = zm1 * ( MAX( 0._wp, sst_s(ji,jj) - ( t_bo(ji,jj) - rt0 ) ) )**zm2  ! Melt speed rate [m/s]
               !
               zda_tot = MIN( zwlat * zperi * rDt_ice, at_i(ji,jj) )                     ! sea ice concentration decrease (>0)

               ! --- Distribute reduction among ice categories and calculate associated ice-ocean fluxes --- !
               ! decrease of concentration for the category jl
               ! each category contributes to melting in proportion to its concentration
               zda = MIN( a_i(ji,jj,jl_cat), zda_tot * a_i(ji,jj,jl_cat) / at_i(ji,jj) )

               zt1 = zda * r1_Dt_ice
               zt2 = h_i(ji,jj,jl_cat) * r1_nlay_i
               
               ! Contribution to salt flux
               sfx_lam(ji,jj) = sfx_lam(ji,jj) + rhoi * zt1 * zt2 * zsum_s_i

               ! Contribution to heat flux into the ocean [W.m-2], (<0)
               hfx_thd(ji,jj) = hfx_thd(ji,jj) - zt1 * ( zt2 * zsum_e_i + h_s(ji,jj,jl_cat)*r1_nlay_s * zsum_e_s )
               
               ! Contribution to mass flux
               wfx_lam(ji,jj) = wfx_lam(ji,jj) + zt1 * ( rhoi * h_i(ji,jj,jl_cat) + rhos * h_s(ji,jj,jl_cat) )

               ! new concentration
               a_i(ji,jj,jl_cat) = a_i(ji,jj,jl_cat) - zda

               ! ensure that h_i = 0 where a_i = 0
               IF( a_i(ji,jj,jl_cat) == 0._wp ) THEN
                  h_i(ji,jj,jl_cat) = 0._wp
                  h_s(ji,jj,jl_cat) = 0._wp
                  ll_ice_present(ji,jj) = .false.
               ENDIF
               
            ENDIF
         END DO
      END DO
      !$acc end parallel loop
      !
      !$acc end data
      IF( ln_timing )   CALL timing_stop('ice_thd_da')
   END SUBROUTINE ice_thd_da


   SUBROUTINE ice_thd_da_init
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_da_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice thermodynamics
      !!
      !! ** Method  :   Read the namthd_da namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_da
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer
      !!
      NAMELIST/namthd_da/ rn_beta, rn_dmin
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_da)
      READ_NML_CFG(numnam_ice,namthd_da)
      IF(lwm) WRITE( numoni, namthd_da )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_da_init: Ice lateral melting'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_da:'
         WRITE(numout,*) '      Coef. beta for lateral melting param.               rn_beta = ', rn_beta
         WRITE(numout,*) '      Minimum floe diameter for lateral melting param.    rn_dmin = ', rn_dmin
      ENDIF
      !
      !$acc update device( rn_beta, rn_dmin )
      !
   END SUBROUTINE ice_thd_da_init

   !!======================================================================
END MODULE icethd_da
