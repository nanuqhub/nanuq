MODULE icedyn_adv_pra_adv
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_pra_adv   ***
   !!   sea-ice : advection => Prather scheme
   !!======================================================================
   !!----------------------------------------------------------------------
   USE par_oce,  ONLY: jpi, jpj, Nis0, Nie0, Njs0, Nje0, nn_hls
   USE par_kind, ONLY: wp
   USE par_ice,  ONLY: jpl, ln_damage, epsi20

   IMPLICIT NONE
   PRIVATE

   PUBLIC   adv_pra_2d
   PUBLIC   adv_pra_3d
   PUBLIC   adv_pra_4d

   ! Work array (that should remain once for all in the memory of the GPU)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)     ::   sa2d
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:,:)   ::   sa3d
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:,:,:) ::   sa4d

   ! 2D workspace arrays for advection routines   
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   zfld, zf0, zbet
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   zfm, zfx, zfy, zfxx, zfyy, zfxy
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   zpm, zpx, zpy, zpxx, zpyy, zpxy
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   zalg, zalg1, zalg1q
   
   !!----------------------------------------------------------------------
   !! NANUQ_beta
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE adv_pra_2d( kt, pdt, pU, pV, pe1e2, p1_e1e2, pmsk,  pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      INTEGER,                      INTENT(in   ) ::   kt                 ! time step
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step length [s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      IF( MOD( (kt - 1) , 2 ) == 0 ) THEN                           !==  odd ice time step:  adv_x_2d then adv_y_2d  ==!
         !                                                          !--------------------------------------------!
         CALL adv_x_2d( pdt, pU, 1._wp, pe1e2, p1_e1e2, pmsk,   sa2d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_y_2d( pdt, pV, 0._wp, pe1e2, p1_e1e2, pmsk,   sa2d, pF, psx, psxx, psy, psyy, psxy )
         !                                                          !--------------------------------------------!
      ELSE                                                          !== even ice time step:  adv_y_2d then adv_x_2d  ==!
         !                                                          !--------------------------------------------!
         CALL adv_y_2d( pdt, pV, 1._wp, pe1e2, p1_e1e2, pmsk,   sa2d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_x_2d( pdt, pU, 0._wp, pe1e2, p1_e1e2, pmsk,   sa2d, pF, psx, psxx, psy, psyy, psxy )
         !
      ENDIF
   END SUBROUTINE adv_pra_2d


   SUBROUTINE adv_pra_3d( kt, pdt, pU, pV, pe1e2, p1_e1e2, pmsk,  pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      INTEGER,                          INTENT(in   ) ::   kt                 ! time step
      REAL(wp),                         INTENT(in   ) ::   pdt                ! time step length [s]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      IF( MOD( (kt - 1) , 2 ) == 0 ) THEN  !==  odd ice time step:  adv_x_3d then adv_y_3d  ==!
         CALL adv_x_3d( pdt, pU, 1._wp, pe1e2, p1_e1e2, pmsk,   sa3d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_y_3d( pdt, pV, 0._wp, pe1e2, p1_e1e2, pmsk,   sa3d, pF, psx, psxx, psy, psyy, psxy )
      ELSE                                 !== even ice time step:  adv_y_3d then adv_x_3d  ==!
         CALL adv_y_3d( pdt, pV, 1._wp, pe1e2, p1_e1e2, pmsk,   sa3d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_x_3d( pdt, pU, 0._wp, pe1e2, p1_e1e2, pmsk,   sa3d, pF, psx, psxx, psy, psyy, psxy )
      ENDIF
   END SUBROUTINE adv_pra_3d


   SUBROUTINE adv_pra_4d( kt, pdt, klay, pU, pV, pe1e2, p1_e1e2, pmsk,  pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      INTEGER,                               INTENT(in   ) ::   kt                 ! time step
      REAL(wp),                              INTENT(in   ) ::   pdt                ! time step length [s]
      INTEGER,                               INTENT(in   ) ::   klay               ! number of layers
      REAL(wp), DIMENSION(jpi,jpj),          INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj),          INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp), DIMENSION(jpi,jpj),          INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj),          INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      IF( MOD( (kt - 1) , 2 ) == 0 ) THEN  !==  odd ice time step:  adv_x_4d then adv_y_4d  ==!
         CALL adv_x_4d( pdt, klay, pU, 1._wp, pe1e2, p1_e1e2, pmsk,   sa4d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_y_4d( pdt, klay, pV, 0._wp, pe1e2, p1_e1e2, pmsk,   sa4d, pF, psx, psxx, psy, psyy, psxy )
      ELSE                                 !== even ice time step:  adv_y_4d then adv_x_4d  ==!
         CALL adv_y_4d( pdt, klay, pV, 1._wp, pe1e2, p1_e1e2, pmsk,   sa4d, pF, psx, psxx, psy, psyy, psxy )
         CALL adv_x_4d( pdt, klay, pU, 0._wp, pe1e2, p1_e1e2, pmsk,   sa4d, pF, psx, psxx, psy, psyy, psxy )
      ENDIF
   END SUBROUTINE adv_pra_4d









   SUBROUTINE adv_x_2d( pdt, pU, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_x_2d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on x axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj                               ! dummy loop indices
      INTEGER  ::   kj0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zU, zswitch
      !---------------------------------------------------------------------
      !$acc data present( pU, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
            zfld(ji,jj) = pF  (ji,jj) * ze1e2
            zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj) , epsi20 )
            zpx (ji,jj) = psx (ji,jj)
            zpxx(ji,jj) = psxx(ji,jj)
            zpy (ji,jj) = psy (ji,jj)
            zpxy(ji,jj) = psxy(ji,jj)
            zpyy(ji,jj) = psyy(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

# include "icedyn_adv_pra_adv_x.h90"

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            pF  (ji,jj) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km
            psm (ji,jj) = zpm (ji,jj)
            psx (ji,jj) = zpx (ji,jj)
            psxx(ji,jj) = zpxx(ji,jj)
            psy (ji,jj) = zpy (ji,jj)
            psxy(ji,jj) = zpxy(ji,jj)
            psyy(ji,jj) = zpyy(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE adv_x_2d

   SUBROUTINE adv_y_2d( pdt, pV, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_y_2d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on y axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj                               ! dummy loop indices
      INTEGER  ::   ki0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zV, zswitch
      !---------------------------------------------------------------------
      !$acc data present( pV, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
            zfld(ji,jj) = pF  (ji,jj) * ze1e2
            zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj) , epsi20 )
            zpx (ji,jj) = psx (ji,jj)
            zpxx(ji,jj) = psxx(ji,jj)
            zpy (ji,jj) = psy (ji,jj)
            zpxy(ji,jj) = psxy(ji,jj)
            zpyy(ji,jj) = psyy(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

# include "icedyn_adv_pra_adv_y.h90"

      !$acc parallel loop collapse(2)
      DO jj=Njs0-1, Nje0+1
         DO ji=Nis0-1, Nie0+1
            pF  (ji,jj) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km 
            psm (ji,jj) = zpm (ji,jj)
            psx (ji,jj) = zpx (ji,jj)
            psxx(ji,jj) = zpxx(ji,jj)
            psy (ji,jj) = zpy (ji,jj)
            psxy(ji,jj) = zpxy(ji,jj)
            psyy(ji,jj) = zpyy(ji,jj)
         END DO
      END DO
      !$acc end parallel loop

      !$acc end data
   END SUBROUTINE adv_y_2d



   SUBROUTINE adv_x_3d( pdt, pU, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_x_3d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on x axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl                           ! dummy loop indices
      INTEGER  ::   kj0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zU, zswitch
      !---------------------------------------------------------------------
      !$acc data present(pU, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc loop seq
      DO jl = 1, jpl   ! loop on categories
         !
         ! Limitation of moments.
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
               zfld(ji,jj) = pF  (ji,jj,jl) * ze1e2
               zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj,jl) , epsi20 )
               zpx (ji,jj) = psx (ji,jj,jl)
               zpxx(ji,jj) = psxx(ji,jj,jl)
               zpy (ji,jj) = psy (ji,jj,jl)
               zpxy(ji,jj) = psxy(ji,jj,jl)
               zpyy(ji,jj) = psyy(ji,jj,jl)
            END DO
         END DO
         !$acc end parallel loop

#        include "icedyn_adv_pra_adv_x.h90"

         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               pF  (ji,jj,jl) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km
               psm (ji,jj,jl) = zpm (ji,jj)
               psx (ji,jj,jl) = zpx (ji,jj)
               psxx(ji,jj,jl) = zpxx(ji,jj)
               psy (ji,jj,jl) = zpy (ji,jj)
               psxy(ji,jj,jl) = zpxy(ji,jj)
               psyy(ji,jj,jl) = zpyy(ji,jj)
            END DO
         END DO
         !$acc end parallel loop

      END DO !DO jl = 1, jpl

      !$acc end data
   END SUBROUTINE adv_x_3d

   SUBROUTINE adv_y_3d( pdt, pV, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_y_3d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on y axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl                           ! dummy loop indices
      INTEGER  ::   ki0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zV, zswitch
      !---------------------------------------------------------------------
      !$acc data present( pV, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc loop seq
      DO jl = 1, jpl   ! loop on categories
         !
         ! Limitation of moments.
         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
               zfld(ji,jj) = pF  (ji,jj,jl) * ze1e2
               zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj,jl) , epsi20 )
               zpx (ji,jj) = psx (ji,jj,jl)
               zpxx(ji,jj) = psxx(ji,jj,jl)
               zpy (ji,jj) = psy (ji,jj,jl)
               zpxy(ji,jj) = psxy(ji,jj,jl)
               zpyy(ji,jj) = psyy(ji,jj,jl)
            END DO
         END DO
         !$acc end parallel loop

#        include "icedyn_adv_pra_adv_y.h90"

         !$acc parallel loop collapse(2)
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               pF  (ji,jj,jl) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km
               psm (ji,jj,jl) = zpm (ji,jj)
               psx (ji,jj,jl) = zpx (ji,jj)
               psxx(ji,jj,jl) = zpxx(ji,jj)
               psy (ji,jj,jl) = zpy (ji,jj)
               psxy(ji,jj,jl) = zpxy(ji,jj)
               psyy(ji,jj,jl) = zpyy(ji,jj)
            END DO
         END DO
         !$acc end parallel loop

      END DO !DO jl = 1, jpl

      !$acc end data
   END SUBROUTINE adv_y_3d

   SUBROUTINE adv_x_4d( pdt, klay, pU, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_x_4d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on x axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      INTEGER ,                     INTENT(in   ) ::   klay               ! number of layers
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pU                 ! i-direction ice velocity at U-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk                       ! dummy loop indices
      INTEGER  ::   kj0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zU, zswitch
      !---------------------------------------------------------------------
      !$acc data present( pU, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc loop seq
      DO jk = 1, klay   ! loop on layers
         !
         !$acc loop seq
         DO jl = 1, jpl   ! loop on categories
            !
            ! Limitation of moments.
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
                  zfld(ji,jj) = pF (ji,jj,jk,jl) * ze1e2
                  zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj,jk,jl) , epsi20 )
                  zpx (ji,jj) = psx (ji,jj,jk,jl)
                  zpxx(ji,jj) = psxx(ji,jj,jk,jl)
                  zpy (ji,jj) = psy (ji,jj,jk,jl)
                  zpxy(ji,jj) = psxy(ji,jj,jk,jl)
                  zpyy(ji,jj) = psyy(ji,jj,jk,jl)
               END DO
            END DO
            !$acc end parallel loop

#           include "icedyn_adv_pra_adv_x.h90"

            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  pF  (ji,jj,jk,jl) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km
                  psm (ji,jj,jk,jl) = zpm (ji,jj)
                  psx (ji,jj,jk,jl) = zpx (ji,jj)
                  psxx(ji,jj,jk,jl) = zpxx(ji,jj)
                  psy (ji,jj,jk,jl) = zpy (ji,jj)
                  psxy(ji,jj,jk,jl) = zpxy(ji,jj)
                  psyy(ji,jj,jk,jl) = zpyy(ji,jj)
               END DO
            END DO
            !$acc end parallel loop

         END DO !DO jl = 1, jpl

      END DO !DO jk = 1, klay

      !$acc end data
   END SUBROUTINE adv_x_4d

   SUBROUTINE adv_y_4d( pdt, klay, pV, pcrh, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_y_4d  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on y axis
      !!---------------------------------------------------------------------
      REAL(wp),                     INTENT(in   ) ::   pdt                ! time step
      INTEGER ,                     INTENT(in   ) ::   klay               ! number of layers
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pV                 ! j-direction ice velocity at V-point [m/s]
      REAL(wp),                     INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pe1e2, p1_e1e2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pmsk
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   pF                 ! field to be advected
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(jpi,jpj,klay,jpl), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk                       ! dummy loop indices
      INTEGER  ::   ki0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zV, zswitch
      !---------------------------------------------------------------------
      !$acc data present( pV, pe1e2, p1_e1e2, pmsk, psm, pF, psx, psxx, psy, psyy, psxy, zfld, zf0, zbet, zfm, zfx, zfy, zfxx, zfyy, zfxy, zpm, zpx, zpy, zpxx, zpyy, zpxy, zalg, zalg1, zalg1q )

      !$acc loop seq
      DO jk = 1, klay   ! loop on layers

         !$acc loop seq
         DO jl = 1, jpl   ! loop on categories

            ! Limitation of moments.
            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  ze1e2 = pe1e2(ji,jj) * 1.E-6_wp  ! needs to be in km
                  zfld(ji,jj) = pF (ji,jj,jk,jl) * ze1e2
                  zpm (ji,jj) = MAX( pcrh * ze1e2 + ( 1._wp - pcrh ) * psm(ji,jj,jk,jl) , epsi20 )
                  zpx (ji,jj) = psx (ji,jj,jk,jl)
                  zpxx(ji,jj) = psxx(ji,jj,jk,jl)
                  zpy (ji,jj) = psy (ji,jj,jk,jl)
                  zpxy(ji,jj) = psxy(ji,jj,jk,jl)
                  zpyy(ji,jj) = psyy(ji,jj,jk,jl)
               END DO
            END DO
            !$acc end parallel loop

#           include "icedyn_adv_pra_adv_y.h90"

            !$acc parallel loop collapse(2)
            DO jj=Njs0-1, Nje0+1
               DO ji=Nis0-1, Nie0+1
                  pF  (ji,jj,jk,jl) = zfld(ji,jj) * p1_e1e2(ji,jj) * 1.E6_wp   ! needs to be in km
                  psm (ji,jj,jk,jl) = zpm (ji,jj)
                  psx (ji,jj,jk,jl) = zpx (ji,jj)
                  psxx(ji,jj,jk,jl) = zpxx(ji,jj)
                  psy (ji,jj,jk,jl) = zpy (ji,jj)
                  psxy(ji,jj,jk,jl) = zpxy(ji,jj)
                  psyy(ji,jj,jk,jl) = zpyy(ji,jj)
               END DO
            END DO
            !$acc end parallel loop

         END DO !DO jl = 1, jpl

      END DO !DO jk = 1, klay

      !$acc end data
   END SUBROUTINE adv_y_4d

END MODULE icedyn_adv_pra_adv
