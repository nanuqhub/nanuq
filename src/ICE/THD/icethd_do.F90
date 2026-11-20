MODULE icethd_do
   !!======================================================================
   !!                       ***  MODULE icethd_do   ***
   !!   sea-ice: sea ice growth in the leads (open water)
   !!======================================================================
   !! History :       !  2005-12  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd_do        : ice growth in open water (=lateral accretion of ice)
   !!   ice_thd_do_init   : initialization
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_oce
   USE dom_oce , ONLY : xmskt, umask, vmask !, smask0
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE oss_nnq , ONLY : sss_s
   USE sbc_ice , ONLY : taux_ai_t, tauy_ai_t
   USE icectl         ! sea-ice: conservation
   USE icethd_sal     ! sea-ice: salinity profiles

   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_do        ! called by ice_thd
   PUBLIC   ice_thd_do_init   ! called by ice_stp
   !
   !                             !!** namelist (namthd_do) **
   REAL(wp) ::   rn_hinew         !  thickness for new ice formation (m)
   LOGICAL  ::   ln_frazil        !  use of frazil ice collection as function of wind (T) or not (F)
   REAL(wp) ::   rn_maxfraz       !  maximum portion of frazil ice collecting at the ice bottom
   REAL(wp) ::   rn_vfraz         !  threshold drift speed for collection of bottom frazil ice
   REAL(wp) ::   rn_Cfraz         !  squeezing coefficient for collection of bottom frazil ice
   !$acc declare create( rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz )

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 1.0.2, Brodeau (2026)
   !! -- a standalone fork of NEMO/SI3 version 5 ---
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_do
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_thd_do  ***
      !!
      !! ** Purpose : Computation of the evolution of the ice thickness and
      !!              concentration as a function of the heat balance in the leads
      !!
      !! ** Method  : Ice is formed in the open water when ocean looses heat
      !!              (heat budget of open water is negative) following
      !!
      !!       (dA/dt)acc = F[ (1-A)/(1-a) ] * [ Bl / (Li*h0) ]
      !!          where - h0 is the thickness of ice created in the lead
      !!                - a is a minimum fraction for leads
      !!                - F is a monotonic non-increasing function defined as:
      !!                  F(X)=( 1 - X**exld )**(1.0/exld)
      !!                - exld is the exponent closure rate (=2 default val.)
      !!
      !! ** Action : - Adjustment of snow and ice thicknesses and heat
      !!                content in brine pockets
      !!             - Updating ice internal temperature
      !!             - Computation of variation of ice volume and mass
      !!             - Computation of a_i after lateral accretion and
      !!               update h_s, h_i
      !!
      !! ** Involves : qlead,a_i,at_i,e_i,v_i,szv_i,hi_max,t_bo,hfx_thd,hfx_opw,wfx_opw,sfx_opw,sss_s,u_ice,v_ice,taux_ai_t,tauy_ai_t
      !! ** Updates :  a_i,e_i,hfx_opw,hfx_thd,sfx_opw,szv_i,v_i,wfx_opw
      !!
      !!------------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl
      !
      REAL(wp) ::   ztmelts
      REAL(wp) ::   zdE
      REAL(wp) ::   zQm          ! enthalpy exchanged with the ocean (J/m2, >0 towards ocean)
      REAL(wp) ::   zEi          ! sea ice specific enthalpy (J/kg)
      REAL(wp) ::   zEw          ! seawater specific enthalpy (J/kg)
      REAL(wp) ::   zfmdt        ! mass flux x time step (kg/m2, >0 towards ocean)
      !
      INTEGER  ::   jcat        ! indexes of categories where new ice grows
      INTEGER  ::   npti
      !
      REAL(wp) ::   zAt, zdum
      REAL(wp) ::   zv_newfra
      REAL(wp) ::   zv_newice   ! volume of accreted ice
      REAL(wp) ::   za_newice   ! fractional area of accreted ice
      REAL(wp) ::   ze_newice   ! heat content of accreted ice
      REAL(wp) ::   zo_newice   ! age of accreted ice
      REAL(wp) ::   zdv_res     ! residual volume in case of excessive heat budget
      REAL(wp) ::   zda_res     ! residual area in case of excessive heat budget
      REAL(wp) ::   zv_frazb    ! accretion of frazil ice at the ice bottom
      REAL(wp) ::   zs_newice   ! salinity of accreted ice
      LOGICAL  ::   lDo
      !
      ! Frazil ice collection thickness and fraction, at the current grid point.
      ! They used to be permanent 2D module arrays (`ht_i_new`, `fraz_frac`) computed in a separate routine
      ! (`ice_thd_frazil`, called at the start of `ice_thd`). Their inputs (qlead, wind stress, ice velocity, at_i)
      ! do not change between that call and `ice_thd_do`, so they are now computed here, inside the main kernel,
      ! as local scalars: no permanent memory, and no extra kernel.
      REAL(wp) ::   zht_i_new   ! thickness of the new ice formed in open water (collection thickness) (m)
      REAL(wp) ::   zfraz_frac  ! fraction of frazil ice accreted at the ice bottom (-)
      ! Local scalars and parameters for the frazil computation (formerly in `ice_thd_frazil`):
      INTEGER  ::   iter
      REAL(wp) ::   zfr_dum, ztmp, zhti, zvfrx, zvgx, ztaux, zf, ztenagm, zvfry, zvgy, ztauy, zvrel2, zfp, ztwogp
      REAL(wp), PARAMETER ::   zcai    = 1.4e-3_wp                       ! ice-air drag (clem: should be dependent on coupling/forcing used)
      REAL(wp), PARAMETER ::   zhicrit = 0.04_wp                         ! frazil ice thickness
      REAL(wp), PARAMETER ::   zsqcd   = 1.0_wp / SQRT( 1.3_wp * zcai )  ! 1/SQRT(airdensity*drag)
      REAL(wp), PARAMETER ::   zgamafr = 0.03_wp
      !
      REAL(wp), DIMENSION(jpl)        ::   zv_b    ! old volume of ice in category jl
      REAL(wp), DIMENSION(jpl)        ::   za_b    ! old area of ice in category jl
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zh_i_o, ze_i_o, zs_i_o
      !
      ! For "manual inlining" of routine `ice_var_vremap`:
      INTEGER  ::   jk0, jk1   !  old/new layer indices
      REAL(wp) ::   zhnew      ! new layers thicknesses
      ! For ice enthalpy and salt content:
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zxi_cum0, zhi_cum0   ! old cumulative enthlapies/salinities and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zxi_cum1, zhi_cum1   ! new cumulative enthlapies/salinities and layers interfaces
      !!-----------------------------------------------------------------------!
      IF( ln_timing    )   CALL timing_start('icethd_do')
      !$acc data present( qlead,hi_max,t_bo,hfx_thd,hfx_opw,wfx_opw,sfx_opw,sss_s,u_ice,v_ice,taux_ai_t,tauy_ai_t ) create( zv_b,za_b,zh_i_o,ze_i_o,zs_i_o, zxi_cum0,zhi_cum0,zxi_cum1,zhi_cum1 )

      ! Reduced gravity, needed for the frazil collection thickness (only if `ln_frazil=T`)
      ztwogp  = 2._wp * rho0 / ( grav * 0.3_wp * ( rho0 - rhoi ) )

      !IF( ln_icediachk )   CALL ice_cons_hsm( 0, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft )
      !IF( ln_icediachk )   CALL ice_cons2D  ( 0, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft )

      ! Identify grid points where new ice forms
      npti = 0
      !$acc parallel loop collapse(2) reduction(+:npti)
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            npti = npti + MERGE( 1, 0,  qlead(ji,jj) < 0._wp )
         END DO
      END DO
      !$acc end parallel loop

      IF( npti > 0 ) THEN

         !$acc parallel loop collapse(2)  private( zv_b,za_b,zh_i_o,ze_i_o,zs_i_o, zxi_cum0,zhi_cum0,zxi_cum1,zhi_cum1 )
         DO jj=Njs0, Nje0
            DO ji=Nis0, Nie0

               !------------------------------------------------------------------------------!
               ! 1) Compute thickness, salinity, enthalpy, age, area and volume of new ice
               !------------------------------------------------------------------------------!
               ! it occurs if cooling
               zAt = 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  zAt = zAt + a_i(ji,jj,jl)
               END DO

               IF(qlead(ji,jj) < 0._wp) THEN
                  ! Convert units for ice internal energy and salt content
                  !$acc loop seq
                  DO jl = 1, jpl
                     zdum = REAL( nlay_i, wp ) / MAX( v_i(ji,jj,jl), epsi20 )
                     lDo = v_i(ji,jj,jl) > 0._wp
                     !$acc loop seq
                     DO jk = 1, nlay_i
                        e_i  (ji,jj,jk,jl) = MERGE(   e_i(ji,jj,jk,jl) * zdum , 0._wp , lDo )
                     END DO
                     !$acc loop seq
                     DO jk = 1, nlay_i
                        szv_i(ji,jj,jk,jl) = MERGE( szv_i(ji,jj,jk,jl) * zdum , 0._wp , lDo )
                     END DO
                  END DO
               ENDIF

               ! --- Salinity of new ice --- !
               ! Sice = F(z,t) [Griewank and Notz 2013 ; Rees Jones and Worster 2014]
               zs_newice = MERGE( rn_sinew * sss_s(ji,jj) , zs_newice , qlead(ji,jj) < 0._wp )


               !                       ! ==================== !
               !                       ! Start main loop here !
               !                       ! ==================== !
               IF( qlead(ji,jj) < 0._wp ) THEN ! qlead is the heat budget in the first ocean level. Only grow ice when it is negative

                  !---------------------------------------------------------!
                  ! Collection thickness of ice formed in leads and polynyas
                  !---------------------------------------------------------!
                  ! (formerly done in routine `ice_thd_frazil`, same computation)
                  ! zht_i_new is the thickness of new ice formed in open water
                  ! zht_i_new can be either prescribed (ln_frazil=F) or computed (ln_frazil=T)
                  ! Frazil ice forms in open water, is transported by wind, accumulates at the edge of the consolidated ice edge
                  ! where it forms aggregates of a specific thickness called collection thickness.
                  ! Default values (prescribed new ice thickness, no frazil accreted at the ice bottom):
                  zfraz_frac = 0._wp
                  zht_i_new  = rn_hinew
                  !
                  IF( ln_frazil ) THEN
                     ! -- Wind stress -- !
                     ztaux = taux_ai_t(ji,jj) * xmskt(ji,jj)
                     ztauy = tauy_ai_t(ji,jj) * xmskt(ji,jj)
                     ! Square root of wind stress
                     ztenagm = SQRT( SQRT( ztaux * ztaux + ztauy * ztauy ) )

                     ! -- Frazil ice velocity -- !
                     lDo = ztenagm >= epsi10
                     zfr_dum = zgamafr * zsqcd / MAX( ztenagm, epsi10 )
                     zvfrx = MERGE( zfr_dum * ztaux , 0._wp , lDo )
                     zvfry = MERGE( zfr_dum * ztauy , 0._wp , lDo )
                     ! -- Pack ice velocity -- !
                     zvgx = ( u_ice(ji-1,jj  ) * umask(ji-1,jj  ,1)  + u_ice(ji,jj) * umask(ji,jj,1) ) * 0.5_wp
                     zvgy = ( v_ice(ji  ,jj-1) * vmask(ji  ,jj-1,1)  + v_ice(ji,jj) * vmask(ji,jj,1) ) * 0.5_wp

                     ! -- Relative frazil/pack ice velocity & fraction of frazil ice-- !
                     lDo = at_i(ji,jj) >= epsi10
                     zvrel2     = MERGE( MAX( (zvfrx - zvgx)*(zvfrx - zvgx) + (zvfry - zvgy)*(zvfry - zvgy), 0.15_wp*0.15_wp ) , 0._wp , lDo )
                     zfraz_frac = MERGE( ( TANH( rn_Cfraz * ( SQRT(zvrel2) - rn_vfraz ) ) + 1._wp ) * 0.5_wp * rn_maxfraz      , 0._wp , lDo )

                     ! -- new ice thickness (iterative loop) -- !
                     zhti = zhicrit +   ( zhicrit + 0.1_wp )    &
                        &                / ( ( zhicrit + 0.1_wp ) * ( zhicrit + 0.1_wp ) -  zhicrit * zhicrit ) * ztwogp * zvrel2
                     iter = 1
                     DO WHILE ( iter < 20 )
                        zfr_dum = zhti - zhicrit
                        ztmp = zhicrit * ztwogp * zvrel2
                        zf  = zfr_dum * ( zhti * zhti - zhicrit * zhicrit ) - zhti * ztmp
                        zfp = zfr_dum * ( 3.0_wp * zhti + zhicrit ) -  ztmp
                        zhti = zhti - zf / MAX( zfp, epsi20 )
                        iter = iter + 1
                     END DO
                     !
                     ! bound zht_i_new (though I don't see why it should be necessary)
                     zht_i_new = MAX( 0.01_wp, MIN( zhti, hi_max(jpl) ) )
                  ENDIF

                  ! Keep old ice areas and volume in memory
                  !$acc loop seq
                  DO jl = 1, jpl
                     zv_b(jl) = v_i(ji,jj,jl)
                     za_b(jl) = a_i(ji,jj,jl)
                  ENDDO

                  ! --- Heat content of new ice --- !
                  ! We assume that new ice is formed at the seawater freezing point
                  ztmelts   = - rTmlt * zs_newice                  ! Melting point (C)
                  ! Enthalpy of the new ice (J/m3). The latent part uses the solid fraction (1 - Tmelt/T), bounded at 0
                  ! (as in SI3) to stay physical if the ice melting point were ever colder than the ocean freezing point
                  ze_newice =   rhoi * (  rcpi  * ( ztmelts - ( t_bo(ji,jj) - rt0 ) )                                   &
                     &                  + rLfus * MAX( 0._wp, 1.0 - ztmelts / MIN( t_bo(ji,jj) - rt0, -epsi10 ) )       &
                     &                  - rcp   *         ztmelts )

                  ! --- Age of new ice --- !
                  zo_newice = 0._wp

                  ! --- Volume of new ice --- !
                  zEi           = - ze_newice * r1_rhoi                  ! specific enthalpy of forming ice [J/kg]

                  zEw           = rcp * ( t_bo(ji,jj) - rt0 )            ! specific enthalpy of seawater at t_bo_1d [J/kg]
                  ! clem: we suppose we are already at the freezing point (condition qlead<0 is satisfyied)

                  zdE           = zEi - zEw                              ! specific enthalpy difference [J/kg]

                  zfmdt         = - qlead(ji,jj) / zdE                   ! Fm.dt [kg/m2] (<0)
                  ! clem: we use qlead instead of zqld (icethd) because we suppose we are at the freezing point
                  zv_newice     = - zfmdt * r1_rhoi

                  zQm           = zfmdt * zEw                            ! heat to the ocean >0 associated with mass flux

                  ! Contribution to heat flux to the ocean [W.m-2], >0
                  hfx_thd(ji,jj) = hfx_thd(ji,jj) + zfmdt * zEw * r1_Dt_ice
                  ! Total heat flux used in this process [W.m-2]
                  hfx_opw(ji,jj) = hfx_opw(ji,jj) - zfmdt * zdE * r1_Dt_ice
                  ! mass flux
                  wfx_opw(ji,jj) = wfx_opw(ji,jj) - zv_newice * rhoi * r1_Dt_ice
                  ! salt flux
                  sfx_opw(ji,jj) = sfx_opw(ji,jj) - zv_newice * rhoi * zs_newice * r1_Dt_ice

                  ! A fraction zfraz_frac of frazil ice is accreted at the ice bottom
                  lDo = zAt > 0._wp
                  zv_frazb  = MERGE(           zfraz_frac   * zv_newice , 0._wp     , lDo )
                  zv_newice = MERGE( ( 1._wp - zfraz_frac ) * zv_newice , zv_newice , lDo )
                  ! --- Area of new ice --- !
                  za_newice = zv_newice / zht_i_new

                  ! --- Redistribute new ice area and volume into ice categories --- !

                  ! --- lateral ice growth --- !
                  ! If lateral ice growth gives an ice concentration > amax, then
                  ! we keep the excessive volume in memory and attribute it later to bottom accretion
                  IF( za_newice >  MAX( 0._wp, rn_amax - zAt ) ) THEN ! max is for roundoff error
                     zda_res   = za_newice - MAX( 0._wp, rn_amax - zAt )
                     zdv_res   = zda_res * zht_i_new
                     za_newice = MAX( 0._wp, za_newice - zda_res )
                     zv_newice = MAX( 0._wp, zv_newice - zdv_res )
                  ELSE
                     zda_res = 0._wp
                     zdv_res = 0._wp
                  ENDIF

                  ! find which category to fill
                  zAt = 0._wp
                  !$acc loop seq
                  DO jl = 1, jpl
                     IF( zht_i_new > hi_max(jl-1) .AND. zht_i_new <= hi_max(jl) ) THEN
                        a_i(ji,jj,jl) = a_i(ji,jj,jl) + za_newice
                        v_i(ji,jj,jl) = v_i(ji,jj,jl) + zv_newice
                        jcat = jl
                     ENDIF
                     zAt = zAt + a_i(ji,jj,jl)
                  END DO

                  ! Heat content
                  jl   = jcat                                       ! categroy in which new ice is put
                  zdum = 1._wp / MAX( v_i(ji,jj,jl), epsi20 )
                  lDo = za_b(jl) > 0._wp
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     e_i  (ji,jj,jk,jl) = MERGE( ( ze_newice * zv_newice +   e_i(ji,jj,jk,jl) * zv_b(jl) ) * zdum , ze_newice , lDo )
                  END DO
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     szv_i(ji,jj,jk,jl) = MERGE( ( zs_newice * zv_newice + szv_i(ji,jj,jk,jl) * zv_b(jl) ) * zdum , zs_newice , lDo )
                  END DO

                  ! --- bottom ice growth + ice enthalpy remapping --- !
                  !$acc loop seq
                  DO jl = 1, jpl
                     ! for remapping
                     !$acc loop seq
                     DO jk = 0, nlay_i+1
                        zh_i_o(jk) = 0._wp
                        ze_i_o(jk) = 0._wp
                        zs_i_o(jk) = 0._wp
                     END DO
                     !$acc loop seq
                     DO jk = 1, nlay_i
                        zdum = v_i(ji,jj,jl) * r1_nlay_i
                        zh_i_o(jk) =                      zdum
                        ze_i_o(jk) =   e_i(ji,jj,jk,jl) * zdum
                        zs_i_o(jk) = szv_i(ji,jj,jk,jl) * zdum
                     END DO

                     ! new volumes including lateral/bottom accretion + residual
                     lDo = zAt >= epsi20
                     zv_newfra     = MERGE( ( zdv_res + zv_frazb ) * a_i(ji,jj,jl) / MAX( zAt , epsi20 ) , 0._wp , lDo )
                     a_i(ji,jj,jl) = MERGE(  a_i(ji,jj,jl)                                               , 0._wp , lDo )
                     v_i(ji,jj,jl) = v_i(ji,jj,jl) + zv_newfra
                     ! for remapping
                     zh_i_o(nlay_i+1) = zv_newfra
                     ze_i_o(nlay_i+1) = ze_newice * zv_newfra
                     zs_i_o(nlay_i+1) = zs_newice * zv_newfra

                     ! --- Ice enthalpy remapping --- !
                     !CALL ice_var_vremap( zh_i_o, ze_i_o,   e_i(ji,jj,:,jl) )
                     !  ==> inlining, better for GPU...
# include            "ice_var_vremap_e.h90"

                     ! --- Ice salt content remapping --- !
                     !CALL ice_var_vremap( zh_i_o, zs_i_o, szv_i(ji,jj,:,jl) )
# include            "ice_var_vremap_s.h90"
                     !$acc loop seq
                     DO jk1 = 1, nlay_i
                        szv_i(ji,jj,jk1,jl) = MAX( 0._wp, zxi_cum1(jk1) - zxi_cum1(jk1-1) ) * zdum ! max for roundoff error
                     END DO

                  END DO !DO jl = 1, jpl

               ENDIF ! qlead < 0

               !                       ! ================== !
               !                       ! End main loop here !
               !                       ! ================== !
               !
               ! Change units for e_i/szv_i
               lDo = qlead(ji,jj) < 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  zdum = v_i(ji,jj,jl) * r1_nlay_i
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     e_i  (ji,jj,jk,jl) = MERGE(   e_i(ji,jj,jk,jl) * zdum ,   e_i(ji,jj,jk,jl) , lDo )
                  END DO
                  !$acc loop seq
                  DO jk = 1, nlay_i
                     szv_i(ji,jj,jk,jl) = MERGE( szv_i(ji,jj,jk,jl) * zdum , szv_i(ji,jj,jk,jl) , lDo )
                  END DO
               END DO

            END DO !DO ji=Nis0, Nie0
         END DO !DO jj=Njs0, Nje0

      ENDIF !IF( npti > 0 )

      ! the following fields need to be updated on the halos (done in icethd): a_i, v_i, szv_i, e_i
      !
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)

      !$acc end data
      IF( ln_timing    )   CALL timing_stop ('icethd_do')
      !
   END SUBROUTINE ice_thd_do


   SUBROUTINE ice_thd_do_init
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_do_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice growth in the leads
      !!
      !! ** Method  :   Read the namthd_do namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_do
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer
      !!
      NAMELIST/namthd_do/ rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_do)
      READ_NML_CFG(numnam_ice,namthd_do)
      IF(lwm) WRITE( numoni, namthd_do )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_do_init: Ice growth in open water'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_do:'
         WRITE(numout,*) '      ice thickness for lateral accretion                       rn_hinew   = ', rn_hinew
         WRITE(numout,*) '      Frazil ice thickness as a function of wind or not         ln_frazil  = ', ln_frazil
         WRITE(numout,*) '      Maximum proportion of frazil ice collecting at bottom     rn_maxfraz = ', rn_maxfraz
         WRITE(numout,*) '      Threshold relative drift speed for collection of frazil   rn_vfraz   = ', rn_vfraz
         WRITE(numout,*) '      Squeezing coefficient for collection of frazil            rn_Cfraz   = ', rn_Cfraz
      ENDIF
      !
      IF( rn_hinew < rn_himin )   CALL ctl_stop( 'ice_thd_do_init : rn_hinew should be >= rn_himin' )
      !
      !$acc update device( rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz )
      !
   END SUBROUTINE ice_thd_do_init

   !!======================================================================
END MODULE icethd_do
