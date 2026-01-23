MODULE iceitd
   !!======================================================================
   !!                       ***  MODULE iceitd ***
   !!   sea-ice : ice thickness distribution
   !!======================================================================
   !! History :  3.0  !  2005-12  (M. Vancoppenolle) original code (based on CICE)
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_itd_rem   : redistribute ice thicknesses after thermo growth and melt
   !!   itd_glinear   : build g(h) satisfying area and volume constraints
   !!   itd_shiftice  : shift ice across category boundaries, conserving everything
   !!   ice_itd_reb   : rebin ice thicknesses into bounded categories
   !!   ice_itd_init  : read ice thicknesses mean and min from namelist
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst,  ONLY : rt0
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_roundoff
   USE icectl         ! sea-ice: conservation tests
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_itd_init  ! called in icestp
   PUBLIC   ice_itd_rem   ! called in icethd
   PUBLIC   ice_itd_reb   ! called in icecor

   INTEGER            ::   nice_catbnd     ! choice of the type of ice category function
   !                                       ! associated indices:
   INTEGER, PARAMETER ::   np_cathfn = 1   ! categories defined by a function
   INTEGER, PARAMETER ::   np_catusr = 2   ! categories defined by the user
   !$acc declare create( nice_catbnd, np_cathfn, np_catusr )
   !
   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_itd_rem( kt )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_rem ***
      !!
      !! ** Purpose :   computes the redistribution of ice thickness
      !!                after thermodynamic growth of ice thickness
      !!
      !! ** Method  :   Linear remapping
      !!
      !! References :   W.H. Lipscomb, JGR 2001
      !!------------------------------------------------------------------
      INTEGER , INTENT(in) ::   kt      ! Ocean time step
      !!------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jcat !, npti ! dummy loop index
      REAL(wp) ::   zx1, zwk1, zdh0, zetamin, zdamax   ! local scalars
      REAL(wp) ::   zx2, zwk2, zda0, zetamax           !   -      -
      REAL(wp) ::   zx3
      REAL(wp) ::   zslope          ! used to compute local thermodynamic "speeds"
      !
      LOGICAL ::   lptidx
      INTEGER , DIMENSION(jpl-1) ::   jdonor          ! donor category index
      REAL(wp), DIMENSION(jpl)   ::   zdhice          ! ice thickness increment
      REAL(wp), DIMENSION(jpl)   ::   zg0, zg1          ! coefficients for fitting the line of the ITD
      REAL(wp), DIMENSION(jpl)   ::   zhL, zhR          ! left and right boundary for the ITD for each thickness
      REAL(wp), DIMENSION(jpl-1) ::   zdaice, zdvice  ! local increment of ice area and volume
      REAL(wp) ::   zhb0, zhb1      ! category boundaries for thinnes categories
      REAL(wp), DIMENSION(0:jpl) ::   zhbnew          ! new boundaries of ice categories
      REAL(wp) :: zAt  ! mean sea-ice concentration
      !
      !
      INTEGER  ::   jk, jl2, jl1           ! local integers
      REAL(wp) ::   zworka, zworkv, ztrans ! ice/snow transferred
      REAL(wp) ::   ztmp, za1, za2             ! workspace
      REAL(wp), DIMENSION(jpl) ::   zaTsfn           !  -    -
      !
      !!------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceitd_rem')
      !$acc data present( a_i, a_i_b, e_i, e_s, h_i, h_i_b, oa_i, sv_i, szv_i, t_su, v_i, v_s, hi_max, hi_mean )
      !$acc data create( jdonor, zdhice, zg0, zg1, zhL, zhR, zdaice, zdvice, zhbnew, zaTsfn )

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_itd_rem: remapping ice thickness distribution'

      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (0, 'iceitd_rem',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)

      !$acc parallel loop collapse(2) private( jdonor, zdhice, zg0, zg1, zhL, zhR, zdaice, zdvice, zhbnew, zaTsfn )
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1

            !$acc loop seq
            DO jl = 0, jpl
               zhbnew(jl) = 0._wp  !LOLO !!! Initialize it with zeros!
            END DO
            
            !-----------------------------------------------------------------------------------------------
            !  1) Identify grid cells with ice
            !-----------------------------------------------------------------------------------------------
            zAt   = 0._wp
            !$acc loop seq
            DO jl=1, jpl
               zAt = zAt + a_i(ji,jj,jl)
            END DO
            lptidx = ( zAt > epsi10 )


            IF( lptidx ) THEN
               !-----------------------------------------------------------------------------------------------
               !  2) Compute new category boundaries
               !-----------------------------------------------------------------------------------------------
               !
               ! Compute thickness change in each ice category
               !$acc loop seq
               DO jl = 1, jpl
                  zdhice(jl) = 0._wp
                  zhbnew(jl) = 0._wp
                  IF( a_i(ji,jj,jl) > epsi10 )   zdhice(jl) = h_i(ji,jj,jl) - h_i_b(ji,jj,jl)
               END DO
               !
               ! --- New boundaries for category 1:jpl-1 --- !
               !$acc loop seq
               DO jl = 1, jpl - 1
                  ! --- New boundary: Hn* = Hn + Fn*dt --- !
                  !     Fn*dt = ( fn + (fn+1 - fn)/(hn+1 - hn) * (Hn - hn) ) * dt = zdhice + zslope * (Hmax - h_i_b)
                  za1 = a_i_b(ji,jj,jl)  ;  za2 = a_i_b(ji,jj,jl+1)
                  IF    ( za1 >  epsi10 .AND. za2 >  epsi10 ) THEN   ! a(jl+1) & a(jl) /= 0
                     zslope     = ( zdhice(jl+1) - zdhice(jl) ) / ( h_i_b(ji,jj,jl+1) - h_i_b(ji,jj,jl) )
                     zhbnew(jl) = hi_max(jl) + zdhice(jl) + zslope * ( hi_max(jl) - h_i_b(ji,jj,jl) )
                  ELSEIF( za1 >  epsi10 .AND. za2 <= epsi10 ) THEN   ! a(jl+1)=0 => Hn* = Hn + fn*dt
                     zhbnew(jl) = hi_max(jl) + zdhice(jl)
                  ELSEIF( za1 <= epsi10 .AND. za2 >  epsi10 ) THEN   ! a(jl)=0 => Hn* = Hn + fn+1*dt
                     zhbnew(jl) = hi_max(jl) + zdhice(jl+1)
                  ELSE                                                                       ! a(jl+1) & a(jl) = 0
                     zhbnew(jl) = hi_max(jl)
                  ENDIF
                  !
                  ! --- 2 conditions for remapping --- !
                  ! 1) hn(t+1)+espi < Hn* < hn+1(t+1)-epsi
                  !    Note: hn(t+1) must not be too close to either HR or HL otherwise a division by nearly 0 is possible
                  !          in itd_glinear in the case (HR-HL) = 3(Hice - HL) or = 3(HR - Hice)
                  IF( a_i(ji,jj,jl  ) > epsi10 .AND. h_i(ji,jj,jl  ) > ( zhbnew(jl) - epsi10 ) )   lptidx = .false.
                  IF( a_i(ji,jj,jl+1) > epsi10 .AND. h_i(ji,jj,jl+1) < ( zhbnew(jl) + epsi10 ) )   lptidx = .false.
                  !
                  ! 2) Hn-1 < Hn* < Hn+1
                  IF( zhbnew(jl) < hi_max(jl-1) )   lptidx = .false.
                  IF( zhbnew(jl) > hi_max(jl+1) )   lptidx = .false.
                  !
               END DO !DO jl = 1, jpl - 1
               !
               !
               ! --- New boundaries for category jpl --- !
               IF( a_i(ji,jj,jpl) > epsi10 ) THEN
                  zhbnew(jpl) = MAX( hi_max(jpl-1), 3._wp * h_i(ji,jj,jpl) - 2._wp * zhbnew(jpl-1) )
               ELSE
                  zhbnew(jpl) = hi_max(jpl)
               ENDIF
               !
               ! --- 1 additional condition for remapping (1st category) --- !
               ! H0+epsi < h1(t) < H1-epsi
               !    h1(t) must not be too close to either HR or HL otherwise a division by nearly 0 is possible
               !    in itd_glinear in the case (HR-HL) = 3(Hice - HL) or = 3(HR - Hice)
               IF( h_i_b(ji,jj,1) < ( hi_max(0) + epsi10 ) )   lptidx = .false.
               IF( h_i_b(ji,jj,1) > ( hi_max(1) - epsi10 ) )   lptidx = .false.

            ENDIF !IF( lptidx )


            IF( lptidx ) THEN

               !-----------------------------------------------------------------------------------------------
               !  3) Compute g(h)
               !-----------------------------------------------------------------------------------------------

               zhb0 = hi_max(0)   ;   zhb1 = hi_max(1)
               !$acc loop seq
               DO jl=1, jpl
                  zg0(jl) = 0._wp       ;   zg1(jl) = 0._wp
                  zhL(jl) = 0._wp       ;   zhR(jl) = 0._wp
               END DO

               ! Area lost due to melting of thin ice

               !$acc loop seq
               DO jl = 1, jpl
                  !
                  IF( jl == 1 ) THEN
                     IF( a_i(ji,jj,jl) > epsi10 ) THEN
                        !
                        ! --- g(h) for category 1 --- !
                        CALL itd_glinear_sclr( zhb0, zhb1, h_i_b(ji,jj,1), a_i(ji,jj,1),  &   ! in
                           &                   zg0(1), zg1(1), zhL(1), zhR(1)   )   ! out
                        !
                        zdh0 =  h_i(ji,jj,1) - h_i_b(ji,jj,1)
                        IF( zdh0 < 0.0 ) THEN      ! remove area from category 1
                           zdh0 = MIN( -zdh0, hi_max(1) )
                           !Integrate g(1) from 0 to dh0 to estimate area melted
                           zetamax = MIN( zdh0, zhR(1) ) - zhL(1)
                           !
                           IF( zetamax > 0.0 ) THEN
                              zx1    = zetamax
                              zx2    = 0.5 * zetamax * zetamax
                              zda0   = zg1(1) * zx2 + zg0(1) * zx1                ! ice area removed
                              zdamax = a_i(ji,jj,1) * (1.0 - h_i(ji,jj,1) / h_i_b(ji,jj,1) ) ! Constrain new thickness <= h_i
                              zda0   = MIN( zda0, zdamax )                            ! ice area lost due to melting of thin ice (zdamax > 0)
                              ! Remove area, conserving volume
                              h_i(ji,jj,1) = h_i(ji,jj,1) * a_i(ji,jj,1) / ( a_i(ji,jj,1) - zda0 )
                              a_i(ji,jj,1) = a_i(ji,jj,1) - zda0
                              v_i(ji,jj,1) = a_i(ji,jj,1) * h_i(ji,jj,1) ! useless ?
                           ENDIF
                           !
                        ELSE ! if ice accretion zdh0 > 0
                           ! zhbnew was 0, and is shifted to the right to account for thin ice growth in openwater (F0 = f1)
                           zhbnew(0) = MIN( zdh0, hi_max(1) )
                        ENDIF
                        !
                     ENDIF
                     !
                  ENDIF !IF( jl == 1 )

                  ! --- g(h) for each thickness category --- !
                  CALL itd_glinear_sclr( zhbnew(jl-1), zhbnew(jl), h_i(ji,jj,jl), a_i(ji,jj,jl),  &   ! in
                     &                   zg0(jl), zg1(jl), zhL(jl), zhR(jl)   )   ! out

                  !
               END DO !DO jl = 1, jpl


               !-----------------------------------------------------------------------------------------------
               !  4) Compute area and volume to be shifted across each boundary (Eq. 18)
               !-----------------------------------------------------------------------------------------------
               !
               !$acc loop seq
               DO jl = 1, jpl - 1

                  ! left and right integration limits in eta space
                  IF (zhbnew(jl) > hi_max(jl)) THEN ! Hn* > Hn => transfer from jl to jl+1
                     zetamin = MAX( hi_max(jl)   , zhL(jl) ) - zhL(jl)   ! hi_max(jl) - zhL
                     zetamax = MIN( zhbnew(jl), zhR(jl) ) - zhL(jl)   ! zhR - zhL
                     jdonor(jl) = jl
                  ELSE                                 ! Hn* <= Hn => transfer from jl+1 to jl
                     zetamin = 0.0
                     zetamax = MIN( hi_max(jl), zhR(jl+1) ) - zhL(jl+1)  ! hi_max(jl) - zhL
                     jdonor(jl) = jl + 1
                  ENDIF
                  zetamax = MAX( zetamax, zetamin ) ! no transfer if etamax < etamin
                  !
                  zx1  = zetamax - zetamin
                  zwk1 = zetamin * zetamin
                  zwk2 = zetamax * zetamax
                  zx2  = 0.5 * ( zwk2 - zwk1 )
                  zwk1 = zwk1 * zetamin
                  zwk2 = zwk2 * zetamax
                  zx3  = 1.0 / 3.0 * ( zwk2 - zwk1 )
                  jcat = jdonor(jl)
                  zdaice(jl) = zg1(jcat)*zx2 + zg0(jcat)*zx1
                  zdvice(jl) = zg1(jcat)*zx3 + zg0(jcat)*zx2 + zdaice(jl)*zhL(jcat)
                  !
               END DO


               !----------------------------------------------------------------------------------------------
               ! 5) Shift ice between categories
               !----------------------------------------------------------------------------------------------

               ! Inlining of `itd_shiftice_gpu`:
               !#################################
               !
               !----------------------------------------------------------------------------------------------
               ! 5.1) Define a variable equal to a_i*T_su
               !----------------------------------------------------------------------------------------------
               !$acc loop seq
               DO jl = 1, jpl
                  zaTsfn(jl) = a_i(ji,jj,jl) * t_su(ji,jj,jl)
               END DO
               !-------------------------------------------------------------------------------
               ! 5.2) Transfer volume and energy between categories
               !-------------------------------------------------------------------------------
               !$acc loop seq
               DO jl = 1, jpl - 1
                  jl1 = jdonor(jl)
                  IF( jl1 > 0 ) THEN
                     IF ( jl1 == jl  ) THEN
                        jl2 = jl1+1
                     ELSE
                        jl2 = jl
                     ENDIF
                     IF( v_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworkv = zdvice(jl) / v_i(ji,jj,jl1)
                     ELSE
                        zworkv = 0._wp
                     ENDIF
                     IF( a_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworka = zdaice(jl) / a_i(ji,jj,jl1)
                     ELSE
                        zworka = 0._wp
                     ENDIF
                     !
                     a_i(ji,jj,jl1) = a_i(ji,jj,jl1) - zdaice(jl)       ! Ice areas
                     a_i(ji,jj,jl2) = a_i(ji,jj,jl2) + zdaice(jl)
                     !
                     v_i(ji,jj,jl1) = v_i(ji,jj,jl1) - zdvice(jl)       ! Ice volumes
                     v_i(ji,jj,jl2) = v_i(ji,jj,jl2) + zdvice(jl)
                     !
                     ztrans         = v_s(ji,jj,jl1) * zworkv              ! Snow volumes
                     v_s(ji,jj,jl1) = v_s(ji,jj,jl1) - ztrans
                     v_s(ji,jj,jl2) = v_s(ji,jj,jl2) + ztrans
                     !
                     ztrans          = oa_i(ji,jj,jl1) * zworka            ! Ice age
                     oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1) - ztrans
                     oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ztrans
                     !
                     ztrans          = zaTsfn(jl1) * zworka             ! Surface temperature
                     zaTsfn(jl1)  = zaTsfn(jl1) - ztrans
                     zaTsfn(jl2)  = zaTsfn(jl2) + ztrans
                     !
                     ! Screws up ACC if uncommented (because not declared on GPU):
                     !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   ztrans          = a_ip(ji,jj,jl1) * zworka         ! Pond fraction
                     !   a_ip(ji,jj,jl1) = a_ip(ji,jj,jl1) - ztrans
                     !   a_ip(ji,jj,jl2) = a_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   ztrans          = v_ip(ji,jj,jl1) * zworkv         ! Pond volume
                     !   v_ip(ji,jj,jl1) = v_ip(ji,jj,jl1) - ztrans
                     !   v_ip(ji,jj,jl2) = v_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   IF ( ln_pnd_lids ) THEN                            ! Pond lid volume
                     !      ztrans          = v_il(ji,jj,jl1) * zworkv
                     !      v_il(ji,jj,jl1) = v_il(ji,jj,jl1) - ztrans
                     !      v_il(ji,jj,jl2) = v_il(ji,jj,jl2) + ztrans
                     !   ENDIF
                     !ENDIF
                     !
                     !$acc loop seq
                     DO jk = 1, nlay_s                                     ! Snow heat content
                        ztrans            = e_s(ji,jj,jk,jl1) * zworkv
                        e_s(ji,jj,jk,jl1) = e_s(ji,jj,jk,jl1) - ztrans
                        e_s(ji,jj,jk,jl2) = e_s(ji,jj,jk,jl2) + ztrans
                     END DO
                     !$acc loop seq
                     DO jk = 1, nlay_i                                     ! Ice heat content
                        ztrans            = e_i(ji,jj,jk,jl1) * zworkv
                        e_i(ji,jj,jk,jl1) = e_i(ji,jj,jk,jl1) - ztrans
                        e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + ztrans
                     END DO
                     !                                                     ! Ice salinity
                     IF( nn_icesal == 4 ) THEN
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           ztrans              = szv_i(ji,jj,jk,jl1) * zworkv
                           szv_i(ji,jj,jk,jl1) = szv_i(ji,jj,jk,jl1) - ztrans
                           szv_i(ji,jj,jk,jl2) = szv_i(ji,jj,jk,jl2) + ztrans
                        END DO
                     ELSE
                        ztrans          = sv_i(ji,jj,jl1) * zworkv
                        sv_i(ji,jj,jl1) = sv_i(ji,jj,jl1) - ztrans
                        sv_i(ji,jj,jl2) = sv_i(ji,jj,jl2) + ztrans
                     ENDIF
                     !
                  ENDIF !IF( jl1 > 0 )
                  !
               END DO !DO jl = 1, jpl - 1

               !-------------------
               ! 5.3) roundoff errors
               !-------------------
               !$acc loop seq
               DO jl = 1, jpl
                  ! clem: The transfer between one category to another can lead to very small negative values (-1.e-20)
                  !       because of truncation error ( i.e. 1. - 1. /= 0 )
                  a_i(ji,jj,jl) = MAX(a_i(ji,jj,jl), 0._wp)
                  v_i(ji,jj,jl) = MAX(v_i(ji,jj,jl), 0._wp)
                  v_s(ji,jj,jl) = MAX(v_s(ji,jj,jl), 0._wp)
                  oa_i(ji,jj,jl) = MAX(oa_i(ji,jj,jl), 0._wp)
                  !$acc loop seq
                  DO jk=1, nlay_i
                     e_i(ji,jj,jk,jl) = MAX(e_i(ji,jj,jk,jl), 0._wp)
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_s
                     e_s(ji,jj,jk,jl) = MAX(e_s(ji,jj,jk,jl), 0._wp)
                  END DO
                  !
                  ! Screws up ACC if uncommented (because not declared on GPU):
                  !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                  !   a_ip(ji,jj,jl) = MAX(a_ip(ji,jj,jl), 0._wp)
                  !   v_ip(ji,jj,jl) = MAX(v_ip(ji,jj,jl), 0._wp)
                  !   IF( ln_pnd_lids ) THEN
                  !      v_il(ji,jj,jl) = MAX(v_il(ji,jj,jl), 0._wp)
                  !   ENDIF
                  !ENDIF
                  !
                  IF( nn_icesal == 4 ) THEN
                     !$acc loop seq
                     DO jk=1, nlay_i
                        szv_i(ji,jj,jk,jl) = MAX(szv_i(ji,jj,jk,jl), 0._wp)
                     END DO
                  ELSE
                     sv_i(ji,jj,jl) = MAX(sv_i(ji,jj,jl), 0._wp)
                  ENDIF
               END DO

               ! at_i must be <= rn_amax
               ztmp = 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  ztmp = ztmp + a_i(ji,jj,jl)
               END DO

               IF ( ztmp > rn_amax ) THEN
                  !$acc loop seq
                  DO jl = 1, jpl
                     a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / ztmp
                  END DO
               ENDIF

               !-------------------------------------------------------------------------------
               ! 5.4) Update ice thickness and temperature
               !-------------------------------------------------------------------------------
               !$acc loop seq
               DO jl = 1, jpl
                  IF ( a_i(ji,jj,jl) >= epsi20 ) THEN
                     h_i (ji,jj,jl)  =  v_i(ji,jj,jl) / a_i(ji,jj,jl)
                     t_su(ji,jj,jl)  =  zaTsfn(jl) / a_i(ji,jj,jl)
                  ELSE
                     h_i (ji,jj,jl)  = 0._wp
                     t_su(ji,jj,jl)  = rt0
                  ENDIF
               END DO !DO jl = 1, jpl

               ! #### end inlining of `itd_shiftice_gpu` ######



            ENDIF !IF( lptidx )

         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+
      !$acc end parallel loop

      ! the following fields need to be updated in the halos (done afterwards):
      ! a_i, v_i, v_s, sv_i, oa_i, h_i, a_ip, v_ip, v_il, t_su, e_i, e_s
      !
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'iceitd_rem',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      !$acc end data
      !$acc end data
      IF( ln_timing    )   CALL timing_stop ('iceitd_rem')
      !
   END SUBROUTINE ice_itd_rem


   SUBROUTINE itd_glinear( iptidx, HbL, Hbr, phice, paice, pg0, pg1, phL, phR )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_glinear ***
      !!
      !! ** Purpose :   build g(h) satisfying area and volume constraints (Eq. 6 and 9)
      !!
      !! ** Method  :   g(h) is linear and written as: g(eta) = zg1(eta) + zg0
      !!                with eta = h - HL
      !!------------------------------------------------------------------
      LOGICAL,  DIMENSION(jpi,jpj), INTENT(in   ) ::   iptidx
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   HbL, HbR      ! left and right category boundaries
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   phice, paice  ! ice thickness and concentration
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pg0, pg1      ! coefficients in linear equation for g(eta)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   phL, phR      ! min and max value of range over which g(h) > 0
      !
      INTEGER  ::   ji,jj        ! horizontal indices
      REAL(wp) ::   z1_3 , z2_3  ! 1/3 , 2/3
      REAL(wp) ::   zh13         ! HbL + 1/3 * (HbR - HbL)
      REAL(wp) ::   zh23         ! HbL + 2/3 * (HbR - HbL)
      REAL(wp) ::   zdhr         ! 1 / (zhR - zhL)
      REAL(wp) ::   zwk1, zwk2   ! temporary variables
      !!------------------------------------------------------------------
      !
      z1_3 = 1._wp / 3._wp
      z2_3 = 2._wp / 3._wp
      !
      !%acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            !
            IF (iptidx(ji,jj) ) THEN

               IF( paice(ji,jj) > epsi10  .AND. phice(ji,jj) > epsi10 )  THEN
                  !
                  ! Initialize zhL and zhR
                  phL(ji,jj) = HbL(ji,jj)
                  phR(ji,jj) = HbR(ji,jj)
                  !
                  ! Change zhL or zhR if hice falls outside central third of range,
                  ! so that hice is in the central third of the range [HL HR]
                  zh13 = z1_3 * ( 2._wp * phL(ji,jj) +         phR(ji,jj) )
                  zh23 = z1_3 * (         phL(ji,jj) + 2._wp * phR(ji,jj) )
                  !
                  IF    ( phice(ji,jj) < zh13 ) THEN
                     phR(ji,jj) = 3._wp * phice(ji,jj) - 2._wp * phL(ji,jj) ! move HR to the left
                  ELSEIF( phice(ji,jj) > zh23 ) THEN
                     phL(ji,jj) = 3._wp * phice(ji,jj) - 2._wp * phR(ji,jj) ! move HL to the right
                  ENDIF
                  !
                  ! Compute coefficients of g(eta) = zg0 + zg1*eta
                  IF( phR(ji,jj) > phL(ji,jj) ) THEN
                     zdhr = 1._wp / (phR(ji,jj) - phL(ji,jj))
                  ELSE
                     zdhr = 0._wp ! if zhR=zhL=hice => no remapping
                  ENDIF
                  !!zdhr = 1._wp / (phR(ii) - phL(ii))
                  zwk1 = 6._wp * paice(ji,jj) * zdhr
                  zwk2 = ( phice(ji,jj) - phL(ji,jj) ) * zdhr
                  pg0(ji,jj) = zwk1 * ( z2_3 - zwk2 )                    ! Eq. 14
                  pg1(ji,jj) = 2._wp * zdhr * zwk1 * ( zwk2 - 0.5_wp )   ! Eq. 14
                  !
               ELSE  ! remap_flag = .false. or a_i < epsi10
                  phL(ji,jj) = 0._wp
                  phR(ji,jj) = 0._wp
                  pg0(ji,jj) = 0._wp
                  pg1(ji,jj) = 0._wp
               ENDIF
               !
            ENDIF ! iptidx
         END DO
      END DO
      !%acc end parallel loop
      !
   END SUBROUTINE itd_glinear

   SUBROUTINE itd_glinear_sclr( pHbL, pHbR, phice, paice, pg0, pg1, phL, phR )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_glinear_sclr ***
      !!
      !! ** Purpose :   build g(h) satisfying area and volume constraints (Eq. 6 and 9)
      !!
      !! ** Method  :   g(h) is linear and written as: g(eta) = zg1(eta) + zg0
      !!                with eta = h - HL
      !!------------------------------------------------------------------
      !%acc routine
      !!------------------------------------------------------------------
      REAL(wp), INTENT(in   ) ::   pHbL, pHbR      ! left and right category boundaries
      REAL(wp), INTENT(in   ) ::   phice, paice  ! ice thickness and concentration
      REAL(wp), INTENT(inout) ::   pg0, pg1      ! coefficients in linear equation for g(eta)
      REAL(wp), INTENT(inout) ::   phL, phR      ! min and max value of range over which g(h) > 0
      !
      INTEGER  ::   ji,jj        ! horizontal indices
      REAL(wp) ::   z1_3 , z2_3  ! 1/3 , 2/3
      REAL(wp) ::   zh13         ! pHbL + 1/3 * (pHbR - pHbL)
      REAL(wp) ::   zh23         ! pHbL + 2/3 * (pHbR - pHbL)
      REAL(wp) ::   zdhr         ! 1 / (zhR - zhL)
      REAL(wp) ::   zwk1, zwk2   ! temporary variables
      !!------------------------------------------------------------------
      !
      z1_3 = 1._wp / 3._wp
      z2_3 = 2._wp / 3._wp
      !
      IF( paice > epsi10  .AND. phice > epsi10 )  THEN
         !
         ! Initialize zhL and zhR
         phL = pHbL
         phR = pHbR
         !
         ! Change zhL or zhR if hice falls outside central third of range,
         ! so that hice is in the central third of the range [HL HR]
         zh13 = z1_3 * ( 2._wp * phL +         phR )
         zh23 = z1_3 * (         phL + 2._wp * phR )
         !
         IF    ( phice < zh13 ) THEN
            phR = 3._wp * phice - 2._wp * phL ! move HR to the left
         ELSEIF( phice > zh23 ) THEN
            phL = 3._wp * phice - 2._wp * phR ! move HL to the right
         ENDIF
         !
         ! Compute coefficients of g(eta) = zg0 + zg1*eta
         IF( phR > phL ) THEN
            zdhr = 1._wp / (phR - phL)
         ELSE
            zdhr = 0._wp ! if zhR=zhL=hice => no remapping
         ENDIF
         !!zdhr = 1._wp / (phR(ii) - phL(ii))
         zwk1 = 6._wp * paice * zdhr
         zwk2 = ( phice - phL ) * zdhr
         pg0 = zwk1 * ( z2_3 - zwk2 )                    ! Eq. 14
         pg1 = 2._wp * zdhr * zwk1 * ( zwk2 - 0.5_wp )   ! Eq. 14
         !
      ELSE  ! remap_flag = .false. or a_i < epsi10
         phL = 0._wp
         phR = 0._wp
         pg0 = 0._wp
         pg1 = 0._wp
      ENDIF
      !
      !
   END SUBROUTINE itd_glinear_sclr




   SUBROUTINE itd_shiftice( lptidx, kdonor, pdaice, pdvice )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_shiftice ***
      !!
      !! ** Purpose :   shift ice across category boundaries, conserving everything
      !!              ( area, volume, energy, age*vol, and mass of salt )
      !!------------------------------------------------------------------
      LOGICAL , DIMENSION(jpi,jpj),       INTENT(in) ::   lptidx   !
      INTEGER , DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   kdonor   ! donor category index
      REAL(wp), DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   pdaice   ! ice area transferred across boundary
      REAL(wp), DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   pdvice   ! ice volume transferred across boundary
      !
      INTEGER  ::   ji, jj, jl, jk         ! dummy loop indices
      INTEGER  ::   jl2, jl1           ! local integers
      REAL(wp) ::   zworka, zworkv, ztrans, zAt ! ice/snow transferred
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zaTsfn           !  -    -
      !!------------------------------------------------------------------
      !$acc data create( zaTsfn ) present( lptidx, kdonor, pdaice, pdvice, a_i, t_su, a_i, v_i, oa_i, v_s, e_s, e_i, szv_i )

      !----------------------------------------------------------------------------------------------
      ! 1) Define a variable equal to a_i*T_su
      !----------------------------------------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( lptidx(ji,jj) ) THEN
               !$acc loop seq
               DO jl = 1, jpl
                  zaTsfn(ji,jj,jl) = a_i(ji,jj,jl) * t_su(ji,jj,jl)
               END DO
            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !-------------------------------------------------------------------------------
      ! 2) Transfer volume and energy between categories
      !-------------------------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( lptidx(ji,jj) ) THEN

               !$acc loop seq
               DO jl = 1, jpl - 1
                  !
                  jl1 = kdonor(ji,jj,jl)
                  !
                  IF( jl1 > 0 ) THEN
                     !
                     IF ( jl1 == jl  ) THEN
                        jl2 = jl1+1
                     ELSE
                        jl2 = jl
                     ENDIF
                     !
                     IF( v_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworkv = pdvice(ji,jj,jl) / v_i(ji,jj,jl1)
                     ELSE
                        zworkv = 0._wp
                     ENDIF
                     IF( a_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworka = pdaice(ji,jj,jl) / a_i(ji,jj,jl1)
                     ELSE
                        zworka = 0._wp
                     ENDIF
                     !
                     a_i(ji,jj,jl1) = a_i(ji,jj,jl1) - pdaice(ji,jj,jl)       ! Ice areas
                     a_i(ji,jj,jl2) = a_i(ji,jj,jl2) + pdaice(ji,jj,jl)
                     !
                     v_i(ji,jj,jl1) = v_i(ji,jj,jl1) - pdvice(ji,jj,jl)       ! Ice volumes
                     v_i(ji,jj,jl2) = v_i(ji,jj,jl2) + pdvice(ji,jj,jl)
                     !
                     ztrans         = v_s(ji,jj,jl1) * zworkv              ! Snow volumes
                     v_s(ji,jj,jl1) = v_s(ji,jj,jl1) - ztrans
                     v_s(ji,jj,jl2) = v_s(ji,jj,jl2) + ztrans
                     !
                     ztrans          = oa_i(ji,jj,jl1) * zworka            ! Ice age
                     oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1) - ztrans
                     oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ztrans
                     !
                     ztrans          = zaTsfn(ji,jj,jl1) * zworka             ! Surface temperature
                     zaTsfn(ji,jj,jl1)  = zaTsfn(ji,jj,jl1) - ztrans
                     zaTsfn(ji,jj,jl2)  = zaTsfn(ji,jj,jl2) + ztrans
                     !
                     !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   ztrans          = a_ip(ji,jj,jl1) * zworka         ! Pond fraction
                     !   a_ip(ji,jj,jl1) = a_ip(ji,jj,jl1) - ztrans
                     !   a_ip(ji,jj,jl2) = a_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   ztrans          = v_ip(ji,jj,jl1) * zworkv         ! Pond volume
                     !   v_ip(ji,jj,jl1) = v_ip(ji,jj,jl1) - ztrans
                     !   v_ip(ji,jj,jl2) = v_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   IF ( ln_pnd_lids ) THEN                            ! Pond lid volume
                     !      ztrans          = v_il(ji,jj,jl1) * zworkv
                     !      v_il(ji,jj,jl1) = v_il(ji,jj,jl1) - ztrans
                     !      v_il(ji,jj,jl2) = v_il(ji,jj,jl2) + ztrans
                     !   ENDIF
                     !ENDIF
                     !
                     !$acc loop seq
                     DO jk = 1, nlay_s                                     ! Snow heat content
                        ztrans            = e_s(ji,jj,jk,jl1) * zworkv
                        e_s(ji,jj,jk,jl1) = e_s(ji,jj,jk,jl1) - ztrans
                        e_s(ji,jj,jk,jl2) = e_s(ji,jj,jk,jl2) + ztrans
                     END DO
                     !$acc loop seq
                     DO jk = 1, nlay_i                                     ! Ice heat content
                        ztrans            = e_i(ji,jj,jk,jl1) * zworkv
                        e_i(ji,jj,jk,jl1) = e_i(ji,jj,jk,jl1) - ztrans
                        e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + ztrans
                     END DO
                     !                                                     ! Ice salinity
                     IF( nn_icesal == 4 ) THEN
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           ztrans              = szv_i(ji,jj,jk,jl1) * zworkv
                           szv_i(ji,jj,jk,jl1) = szv_i(ji,jj,jk,jl1) - ztrans
                           szv_i(ji,jj,jk,jl2) = szv_i(ji,jj,jk,jl2) + ztrans
                        END DO
                     ELSE
                        ztrans          = sv_i(ji,jj,jl1) * zworkv
                        sv_i(ji,jj,jl1) = sv_i(ji,jj,jl1) - ztrans
                        sv_i(ji,jj,jl2) = sv_i(ji,jj,jl2) + ztrans
                     ENDIF
                     !
                  ENDIF   ! jl1 >0
                  !
               END DO !DO jl = 1, jpl - 1

            ENDIF ! lptidx true
         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop

      !-------------------
      ! 3) roundoff errors
      !-------------------
      ! clem: The transfer between one category to another can lead to very small negative values (-1.e-20)
      !       because of truncation error ( i.e. 1. - 1. /= 0 )

      !CALL ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, lptidx )
      CALL  ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i,                   e_s, e_i, szv_i, lptidx )

      ! at_i must be <= rn_amax
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( lptidx(ji,jj) ) THEN
               !
               zAt = 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  zAt = zAt + a_i(ji,jj,jl)
               END DO

               IF ( zAt > rn_amax ) THEN
                  !$acc loop seq
                  DO jl  = 1, jpl
                     a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / zAt
                  END DO
               ENDIF

            ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !-------------------------------------------------------------------------------
      ! 4) Update ice thickness and temperature
      !-------------------------------------------------------------------------------
      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( lptidx(ji,jj) ) THEN
               !$acc loop seq
               DO jl = 1, jpl
                  IF ( a_i(ji,jj,jl) >= epsi20 ) THEN
                     h_i (ji,jj,jl)  =  v_i(ji,jj,jl) / a_i(ji,jj,jl)
                     t_su(ji,jj,jl)  =  zaTsfn(ji,jj,jl) / a_i(ji,jj,jl)
                  ELSE
                     h_i (ji,jj,jl)  = 0._wp
                     t_su(ji,jj,jl)  = rt0
                  ENDIF
               END DO
            ENDIF
         END DO
      END DO
      !$acc end parallel loop
      !
      !$acc end data
      !
   END SUBROUTINE itd_shiftice


   SUBROUTINE itd_shiftice_gpu( lptidx, kdonor, pdaice, pdvice )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_shiftice_gpu ***
      !!
      !! ** Purpose :   shift ice across category boundaries, conserving everything
      !!              ( area, volume, energy, age*vol, and mass of salt )
      !!------------------------------------------------------------------
      LOGICAL , DIMENSION(jpi,jpj),       INTENT(in) ::   lptidx   !
      INTEGER , DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   kdonor   ! donor category index
      REAL(wp), DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   pdaice   ! ice area transferred across boundary
      REAL(wp), DIMENSION(jpi,jpj,jpl-1), INTENT(in) ::   pdvice   ! ice volume transferred across boundary
      !
      INTEGER  ::   ji, jj, jl, jk         ! dummy loop indices
      INTEGER  ::   jl2, jl1           ! local integers
      REAL(wp) ::   zworka, zworkv, ztrans, zAt ! ice/snow transferred
      REAL(wp), DIMENSION(jpl) ::   zaTsfn           !  -    -
      !!------------------------------------------------------------------
      !$acc data present( lptidx, kdonor, pdaice, pdvice, a_i, t_su, a_i, v_i, oa_i, v_s, e_s, e_i, szv_i ) create( zaTsfn )

      !$acc parallel loop collapse(2) private( zaTsfn )
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            IF( lptidx(ji,jj) ) THEN

               !----------------------------------------------------------------------------------------------
               ! 1) Define a variable equal to a_i*T_su
               !----------------------------------------------------------------------------------------------

               !$acc loop seq
               DO jl = 1, jpl
                  zaTsfn(jl) = a_i(ji,jj,jl) * t_su(ji,jj,jl)
               END DO

               !-------------------------------------------------------------------------------
               ! 2) Transfer volume and energy between categories
               !-------------------------------------------------------------------------------
               !$acc loop seq
               DO jl = 1, jpl - 1
                  !
                  jl1 = kdonor(ji,jj,jl)
                  !
                  IF( jl1 > 0 ) THEN
                     !
                     IF ( jl1 == jl  ) THEN
                        jl2 = jl1+1
                     ELSE
                        jl2 = jl
                     ENDIF
                     !
                     IF( v_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworkv = pdvice(ji,jj,jl) / v_i(ji,jj,jl1)
                     ELSE
                        zworkv = 0._wp
                     ENDIF
                     IF( a_i(ji,jj,jl1) >= epsi10 ) THEN
                        zworka = pdaice(ji,jj,jl) / a_i(ji,jj,jl1)
                     ELSE
                        zworka = 0._wp
                     ENDIF
                     !
                     a_i(ji,jj,jl1) = a_i(ji,jj,jl1) - pdaice(ji,jj,jl)       ! Ice areas
                     a_i(ji,jj,jl2) = a_i(ji,jj,jl2) + pdaice(ji,jj,jl)
                     !
                     v_i(ji,jj,jl1) = v_i(ji,jj,jl1) - pdvice(ji,jj,jl)       ! Ice volumes
                     v_i(ji,jj,jl2) = v_i(ji,jj,jl2) + pdvice(ji,jj,jl)
                     !
                     ztrans         = v_s(ji,jj,jl1) * zworkv              ! Snow volumes
                     v_s(ji,jj,jl1) = v_s(ji,jj,jl1) - ztrans
                     v_s(ji,jj,jl2) = v_s(ji,jj,jl2) + ztrans
                     !
                     ztrans          = oa_i(ji,jj,jl1) * zworka            ! Ice age
                     oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1) - ztrans
                     oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ztrans
                     !
                     ztrans             = zaTsfn(jl1) * zworka             ! Surface temperature
                     zaTsfn(jl1)  = zaTsfn(jl1) - ztrans
                     zaTsfn(jl2)  = zaTsfn(jl2) + ztrans
                     !
                     !IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     !   ztrans          = a_ip(ji,jj,jl1) * zworka         ! Pond fraction
                     !   a_ip(ji,jj,jl1) = a_ip(ji,jj,jl1) - ztrans
                     !   a_ip(ji,jj,jl2) = a_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   ztrans          = v_ip(ji,jj,jl1) * zworkv         ! Pond volume
                     !   v_ip(ji,jj,jl1) = v_ip(ji,jj,jl1) - ztrans
                     !   v_ip(ji,jj,jl2) = v_ip(ji,jj,jl2) + ztrans
                     !   !
                     !   IF ( ln_pnd_lids ) THEN                            ! Pond lid volume
                     !      ztrans          = v_il(ji,jj,jl1) * zworkv
                     !      v_il(ji,jj,jl1) = v_il(ji,jj,jl1) - ztrans
                     !      v_il(ji,jj,jl2) = v_il(ji,jj,jl2) + ztrans
                     !   ENDIF
                     !ENDIF
                     !
                     !$acc loop seq
                     DO jk = 1, nlay_s                                     ! Snow heat content
                        ztrans            = e_s(ji,jj,jk,jl1) * zworkv
                        e_s(ji,jj,jk,jl1) = e_s(ji,jj,jk,jl1) - ztrans
                        e_s(ji,jj,jk,jl2) = e_s(ji,jj,jk,jl2) + ztrans
                     END DO
                     !$acc loop seq
                     DO jk = 1, nlay_i                                     ! Ice heat content
                        ztrans            = e_i(ji,jj,jk,jl1) * zworkv
                        e_i(ji,jj,jk,jl1) = e_i(ji,jj,jk,jl1) - ztrans
                        e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + ztrans
                     END DO
                     !                                                     ! Ice salinity
                     IF( nn_icesal == 4 ) THEN
                        !$acc loop seq
                        DO jk = 1, nlay_i
                           ztrans              = szv_i(ji,jj,jk,jl1) * zworkv
                           szv_i(ji,jj,jk,jl1) = szv_i(ji,jj,jk,jl1) - ztrans
                           szv_i(ji,jj,jk,jl2) = szv_i(ji,jj,jk,jl2) + ztrans
                        END DO
                     ELSE
                        ztrans          = sv_i(ji,jj,jl1) * zworkv
                        sv_i(ji,jj,jl1) = sv_i(ji,jj,jl1) - ztrans
                        sv_i(ji,jj,jl2) = sv_i(ji,jj,jl2) + ztrans
                     ENDIF
                     !
                  ENDIF !IF( jl1 > 0 )
                  !
               END DO !DO jl = 1, jpl - 1

               !-------------------
               ! 3) roundoff errors
               !-------------------
               ! ==> inlining of `ice_var_roundoff`
               !$acc loop seq
               DO jl = 1, jpl
                  ! clem: The transfer between one category to another can lead to very small negative values (-1.e-20)
                  !       because of truncation error ( i.e. 1. - 1. /= 0 )
                  !CALL ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, lptidx )
                  a_i(ji,jj,jl) = MAX(a_i(ji,jj,jl), 0._wp)
                  v_i(ji,jj,jl) = MAX(v_i(ji,jj,jl), 0._wp)
                  v_s(ji,jj,jl) = MAX(v_s(ji,jj,jl), 0._wp)
                  oa_i(ji,jj,jl) = MAX(oa_i(ji,jj,jl), 0._wp)
                  !$acc loop seq
                  DO jk=1, nlay_i
                     e_i(ji,jj,jk,jl) = MAX(e_i(ji,jj,jk,jl), 0._wp)
                  END DO
                  !$acc loop seq
                  DO jk=1, nlay_s
                     e_s(ji,jj,jk,jl) = MAX(e_s(ji,jj,jk,jl), 0._wp)
                  END DO
                  !
                  !IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                  !   a_ip(ji,jj,jl) = MAX(a_ip(ji,jj,jl), 0._wp)
                  !   v_ip(ji,jj,jl) = MAX(v_ip(ji,jj,jl), 0._wp)
                  !   IF( ln_pnd_lids ) THEN
                  !      v_il(ji,jj,jl) = MAX(v_il(ji,jj,jl), 0._wp)
                  !   ENDIF
                  !ENDIF
                  !
                  IF( nn_icesal == 4 ) THEN
                     !$acc loop seq
                     DO jk=1, nlay_i
                        szv_i(ji,jj,jk,jl) = MAX(szv_i(ji,jj,jk,jl), 0._wp)
                     END DO
                  ELSE
                     sv_i(ji,jj,jl) = MAX(sv_i(ji,jj,jl), 0._wp)
                  ENDIF
               END DO
               !! ----- END roundoff error ------

               ! at_i must be <= rn_amax
               !
               zAt = 0._wp
               !$acc loop seq
               DO jl = 1, jpl
                  zAt = zAt + a_i(ji,jj,jl)
               END DO

               IF ( zAt > rn_amax ) THEN
                  !$acc loop seq
                  DO jl  = 1, jpl
                     a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax / zAt
                  END DO
               ENDIF

               !-------------------------------------------------------------------------------
               ! 4) Update ice thickness and temperature
               !-------------------------------------------------------------------------------
               !$acc loop seq
               DO jl = 1, jpl
                  IF ( a_i(ji,jj,jl) >= epsi20 ) THEN
                     h_i (ji,jj,jl) =  v_i(ji,jj,jl) / a_i(ji,jj,jl)
                     t_su(ji,jj,jl) =  zaTsfn(jl)    / a_i(ji,jj,jl)
                  ELSE
                     h_i (ji,jj,jl) = 0._wp
                     t_su(ji,jj,jl) = rt0
                  ENDIF

               END DO !DO jl = 1, jpl

            ENDIF !IF( lptidx(ji,jj) )
         END DO !DO ji=Nis0-1, Nie0+1
      END DO !DO jj=Njs0-1, Nje0+1
      !$acc end parallel loop
      !
      !$acc end data
      !
   END SUBROUTINE itd_shiftice_gpu



   SUBROUTINE ice_itd_reb( kt )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_reb ***
      !!
      !! ** Purpose : rebin - rebins thicknesses into defined categories
      !!
      !! ** Method  : If a category thickness is out of bounds, shift part (for down to top)
      !!              or entire (for top to down) area, volume, and energy
      !!              to the neighboring category
      !!------------------------------------------------------------------
      INTEGER , INTENT (in) ::   kt      ! current time step
      !!------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl !, npti   ! dummy loop indices
      INTEGER  ::   ksatisfy
      REAL(wp) ::   zA
      LOGICAL , DIMENSION(jpi,jpj)       ::   lptidx
      INTEGER , DIMENSION(jpi,jpj,jpl-1) ::   jdonor           ! donor category index
      REAL(wp), DIMENSION(jpi,jpj,jpl-1) ::   zdaice, zdvice   ! ice area and volume transferred
      !!------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceitd_reb')
      !$acc data present( a_i, v_i, hi_max, hi_mean ) create( lptidx, jdonor, zdaice, zdvice )
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_itd_reb: rebining ice thickness distribution'
      !
      !IF( ln_icediachk )   CALL ice_cons_hsm(0, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (0, 'iceitd_reb',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)

      !$acc parallel loop collapse(3)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            DO jl = 1, jpl-1
               jdonor(ji,jj,jl) = 0
               zdaice(ji,jj,jl) = 0._wp
               zdvice(ji,jj,jl) = 0._wp
            END DO
         END DO
      END DO
      !$acc end parallel loop

      !$acc loop seq          !---------------------------------------
      DO jl = 1, jpl-1        ! identify thicknesses that are too big
         !                    !---------------------------------------
         ksatisfy = 0
         !$acc parallel loop collapse(2) reduction(+:ksatisfy)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               IF( a_i(ji,jj,jl) > 0._wp .AND. v_i(ji,jj,jl) > (a_i(ji,jj,jl) * hi_max(jl)) ) THEN
                  lptidx(ji,jj) = .TRUE.
                  ksatisfy = ksatisfy + 1
               ELSE
                  lptidx(ji,jj) = .FALSE.
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
         !
         IF( ksatisfy > 0 ) THEN
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  IF( lptidx(ji,jj) ) THEN
                     zA = a_i(ji,jj,jl)
                     jdonor(ji,jj,jl)  = jl
                     ! how much of a_i you send in cat sup is somewhat arbitrary
                     ! these are from CICE => transfer everything
                     ! these are from LLN => transfer only half of the category
                     zdaice(ji,jj,jl) =                          0.5_wp  * zA
                     zdvice(ji,jj,jl) = v_i(ji,jj,jl) - (1._wp - 0.5_wp) * zA * hi_mean(jl)
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop
            !
            !CALL itd_shiftice( lptidx, jdonor, zdaice, zdvice )  ! Shift jl=>jl+1
            CALL itd_shiftice_gpu( lptidx, jdonor, zdaice, zdvice )  ! Shift jl=>jl+1

            ! Reset shift parameters
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  jdonor(ji,jj,jl) = 0
                  zdaice(ji,jj,jl) = 0._wp
                  zdvice(ji,jj,jl) = 0._wp
               END DO
            END DO
            !$acc end parallel loop

         ENDIF ! IF( ksatisfy > 0 )
         !
      END DO ! DO jl = 1, jpl-1



      !$acc loop seq          !-----------------------------------------
      DO jl = jpl-1, 1, -1    ! Identify thicknesses that are too small
         !                    !-----------------------------------------
         ksatisfy = 0
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               IF( a_i(ji,jj,jl+1) > 0._wp .AND. v_i(ji,jj,jl+1) <= (a_i(ji,jj,jl+1) * hi_max(jl)) ) THEN
                  lptidx(ji,jj) = .TRUE.
                  ksatisfy = ksatisfy + 1
               ELSE
                  lptidx(ji,jj) = .FALSE.
               ENDIF
            END DO
         END DO
         !$acc end parallel loop
         !
         IF( ksatisfy > 0 ) THEN
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  IF( lptidx(ji,jj) ) THEN
                     jdonor(ji,jj,jl) = jl + 1
                     zdaice(ji,jj,jl) = a_i(ji,jj,jl+1)
                     zdvice(ji,jj,jl) = v_i(ji,jj,jl+1)
                  ENDIF
               END DO
            END DO
            !$acc end parallel loop
            !
            !CALL itd_shiftice( lptidx, jdonor, zdaice, zdvice )  ! Shift jl+1=>jl
            CALL itd_shiftice_gpu( lptidx, jdonor, zdaice, zdvice )  ! Shift jl+1=>jl

            ! Reset shift parameters
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  jdonor(ji,jj,jl) = 0
                  zdaice(ji,jj,jl) = 0._wp
                  zdvice(ji,jj,jl) = 0._wp
               END DO
            END DO
            !$acc end parallel loop

         ENDIF ! IF( ksatisfy > 0 )
         !
      END DO !DO jl = jpl-1, 1, -1
      !
      ! clem: those fields must be updated on the halos: a_i, v_i, v_s, sv_i, oa_i, h_i, t_su, a_ip, v_ip, v_il, e_i, e_s
      !       note: ice_itd_reb is called in icedyn
      !             and in icethd (but once the arrays are already updated on the boundaries)
      !
      !IF( ln_icediachk )   CALL ice_cons_hsm(1, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      !IF( ln_icediachk )   CALL ice_cons2D  (1, 'iceitd_reb',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      !$acc end data
      IF( ln_timing    )   CALL timing_stop ('iceitd_reb')
      !
   END SUBROUTINE ice_itd_reb


   SUBROUTINE ice_itd_init
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_init ***
      !!
      !! ** Purpose :   Initializes the ice thickness distribution
      !! ** Method  :   ...
      !! ** input   :   Namelist namitd
      !!-------------------------------------------------------------------
      INTEGER  ::   jl            ! dummy loop index
      INTEGER  ::   ios, ioptio   ! Local integer output status for namelist read
      REAL(wp) ::   zhmax, znum, zden, zalpha   !   -      -
      !
      NAMELIST/namitd/ ln_cat_hfn, rn_himean, ln_cat_usr, rn_catbnd, rn_himin, rn_himax
      !!------------------------------------------------------------------
      !
      rn_catbnd(:) =  0._wp ! Circumvent possible initialization by compiler
      ! to prevent from errors when writing output
      READ_NML_REF(numnam_ice,namitd)
      READ_NML_CFG(numnam_ice,namitd)
      IF(lwm) WRITE( numoni, namitd )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_itd_init: Initialization of ice cat distribution '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namitd: '
         WRITE(numout,*) '      Ice categories are defined by a function of rn_himean**(-0.05)    ln_cat_hfn = ', ln_cat_hfn
         WRITE(numout,*) '         mean ice thickness in the domain                               rn_himean  = ', rn_himean
         WRITE(numout,*) '      Ice categories are defined by rn_catbnd                           ln_cat_usr = ', ln_cat_usr
         WRITE(numout,*) '      minimum ice thickness allowed                                     rn_himin   = ', rn_himin
         WRITE(numout,*) '      maximum ice thickness allowed                                     rn_himax   = ', rn_himax
      ENDIF
      !
      !-----------------------------------!
      !  Thickness categories boundaries  !
      !-----------------------------------!
      !                             !== set the choice of ice categories ==!
      ioptio = 0
      IF( ln_cat_hfn ) THEN
         ioptio = ioptio + 1
         nice_catbnd = np_cathfn
      ENDIF
      IF( ln_cat_usr ) THEN
         ioptio = ioptio + 1
         nice_catbnd = np_catusr
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_itd_init: choose one and only one ice categories boundaries' )
      !
      SELECT CASE( nice_catbnd )
         !                             !------------------------!
      CASE( np_cathfn )                ! h^(-alpha) function
         !                             !------------------------!
         zalpha = 0.05_wp
         zhmax  = 3._wp * rn_himean
         hi_max(0) = 0._wp
         DO jl = 1, jpl
            znum = jpl * ( zhmax+1 )**zalpha
            zden = REAL( jpl-jl , wp ) * ( zhmax + 1._wp )**zalpha + REAL( jl , wp )
            hi_max(jl) = ( znum / zden )**(1./zalpha) - 1
         END DO
         !                             !------------------------!
      CASE( np_catusr )                ! user defined
         !                             !------------------------!
         DO jl = 0, jpl
            hi_max(jl) = rn_catbnd(jl)
         END DO
         !
      END SELECT
      !
      DO jl = 1, jpl                ! mean thickness by category
         hi_mean(jl) = ( hi_max(jl) + hi_max(jl-1) ) * 0.5_wp
      END DO
      !
      hi_max(jpl) = rn_himax        ! set to a big value to ensure that all ice is thinner than hi_max(jpl)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   ===>>>   resulting thickness category boundaries :'
      IF(lwp) WRITE(numout,*) '            hi_max(:)= ', hi_max(0:jpl)
      !
      IF( hi_max(1) < rn_himin )   CALL ctl_stop('ice_itd_init: the upper bound of the 1st category must be bigger than rn_himin')
      !
      !$acc update device( ln_cat_hfn, rn_himean, ln_cat_usr, rn_catbnd, rn_himin, rn_himax, nice_catbnd )
      !
# if defined _OPENACC
      PRINT *, ' * info GPU: ice_itd_init() => adding `hi_max` & `hi_mean` arrays to memory!'
      !$acc enter data copyin( hi_max(0:jpl), hi_mean(1:jpl) )
# endif
      !
   END SUBROUTINE ice_itd_init

   !!======================================================================
END MODULE iceitd
