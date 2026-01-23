MODULE icethd_ent
   !!======================================================================
   !!                       ***  MODULE icethd_ent   ***
   !!   sea-ice: redistribution of enthalpy in the ice on the new vertical grid
   !!                       after vertical growth/melt
   !!======================================================================
   !! History :       !  2003-05  (M. Vancoppenolle) Original code in 1D
   !!                 !  2005-07  (M. Vancoppenolle) 3D version
   !!            3.6  !  2014-05  (C. Rousset)       New version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ice_thd_ent   : ice redistribution of enthalpy
   !!----------------------------------------------------------------------
   USE dom_oce        ! domain variables
   USE domain         !
   USE phycst         ! physical constants
   USE par_ice
   USE ice            ! sea-ice: variables
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   USE timing         ! Timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_ent         ! called by icethd and icethd_do
   PUBLIC   ice_thd_ent_scl     ! called by icethd and icethd_do

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: icethd_ent.F90 14778 2021-05-03 08:58:22Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_ent( qnew_2d )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_thd_ent  ***
      !!
      !! ** Purpose :
      !!           This routine computes new vertical grids in the ice,
      !!           and consistently redistributes temperatures.
      !!           Redistribution is made so as to ensure to energy conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_i+2)                 ------------- cum1(nlay_i)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   qnew_2d        ! new enthlapies (J.m-3, remapped)
      !
      INTEGER  :: ji         !  dummy loop indices
      !#thd2d:
      INTEGER  ::   jj, i1,i2, j1,j2
      !#thd2d.
      INTEGER  :: jk0, jk1   !  old/new layer indices
      REAL(wp) ::   zswitch
      !
      !#thd2d:
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zeh_cum0_2d, zh_cum0_2d   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zeh_cum1_2d, zh_cum1_2d   ! new cumulative enthlapies and layers interfaces
      REAL(wp)                        ::   zhnew_2d                  ! new layers thicknesses
      !#thd2d.
      !!-------------------------------------------------------------------
      IF( ln_timing    )   CALL timing_start('ice_thd_ent')

      !--------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy * thickness and layers interfaces
      !--------------------------------------------------------------------------
      !#thd2d:
      i1 = Nis0-1 ; i2 = Nie0+1
      j1 = Njs0-1 ; j2 = Nje0+1

      DO jj=j1, j2
         DO ji=i1, i2

            zeh_cum0_2d(0) = 0._wp
            zh_cum0_2d(0)  = 0._wp

            DO jk0 = 1, nlay_i+2
               zeh_cum0_2d(jk0) = zeh_cum0_2d(jk0-1) + eh_i_old(ji,jj,jk0-1)
               zh_cum0_2d(jk0)  = zh_cum0_2d(jk0-1) + h_i_old(ji,jj,jk0-1)
            END DO

            !------------------------------------
            !  2) Interpolation on the new layers
            !------------------------------------
            ! new layer thickesses
            zhnew_2d = SUM( h_i_old(ji,jj,0:nlay_i+1) ) * r1_nlay_i

            ! new layers interfaces
            zh_cum1_2d(0) = 0._wp

            DO jk1 = 1, nlay_i
               zh_cum1_2d(jk1) = zh_cum1_2d(jk1-1) + zhnew_2d
            END DO

            zeh_cum1_2d(0:nlay_i) = 0._wp

            ! new cumulative q*h => linear interpolation
            DO jk0 = 1, nlay_i+2
               DO jk1 = 1, nlay_i-1
                  IF( zh_cum1_2d(jk1) <= zh_cum0_2d(jk0) .AND. zh_cum1_2d(jk1) > zh_cum0_2d(jk0-1) ) THEN
                     zeh_cum1_2d(jk1) = ( zeh_cum0_2d(jk0-1) * ( zh_cum0_2d(jk0) - zh_cum1_2d(jk1  ) ) +  &
                        &                 zeh_cum0_2d(jk0  ) * ( zh_cum1_2d(jk1) - zh_cum0_2d(jk0-1) ) )  &
                        &             / ( zh_cum0_2d(jk0) - zh_cum0_2d(jk0-1) )
                  ENDIF
               END DO
            END DO
            ! to ensure that total heat content is strictly conserved, set:
            zeh_cum1_2d(nlay_i) = zeh_cum0_2d(nlay_i+2)

            ! new enthalpies
            DO jk1 = 1, nlay_i
               zswitch      = MAX( 0._wp , SIGN( 1._wp , zhnew_2d - epsi20 ) )
               qnew_2d(ji,jj,jk1) = zswitch * MAX( 0._wp, zeh_cum1_2d(jk1) - zeh_cum1_2d(jk1-1) ) / MAX( zhnew_2d, epsi20 ) ! max for roundoff error
            END DO



         END DO
      END DO


      ! --- diag error on heat remapping --- !
      ! comment: if input h_i_old and eh_i_old are already multiplied by a_i (as in icethd_do),
      ! then we should not (* a_i) again but not important since this is just to check that remap error is ~0
      !DO ji = 1, npti
      !   hfx_err_rem_1d(ji) = hfx_err_rem_1d(ji) + a_i_1d(ji) * r1_Dt_ice *  &
      !      &               ( SUM( qnew(ji,1:nlay_i) ) * zhnew(ji) - SUM( eh_i_old(ji,0:nlay_i+1) ) )
      !END DO

      IF( ln_timing    )   CALL timing_stop('ice_thd_ent')

   END SUBROUTINE ice_thd_ent



   SUBROUTINE ice_thd_ent_scl( ph_i_old, peh_i_old, qnew )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_thd_ent_scl  ***
      !!
      !! ** Purpose :
      !!           This routine computes new vertical grids in the ice,
      !!           and consistently redistributes temperatures.
      !!           Redistribution is made so as to ensure to energy conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_i+2)                 ------------- cum1(nlay_i)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(0:nlay_i+1), INTENT(in)    ::   ph_i_old, peh_i_old
      REAL(wp), DIMENSION(nlay_i),     INTENT(inout) ::   qnew        ! new enthlapies (J.m-3, remapped)
      !!-------------------------------------------------------------------
      INTEGER  :: jk0, jk1   !  old/new layer indices
      REAL(wp) ::   zswitch
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zeh_cum0_2d, zh_cum0_2d   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zeh_cum1_2d, zh_cum1_2d   ! new cumulative enthlapies and layers interfaces
      REAL(wp)                        ::   zhnew_2d                  ! new layers thicknesses
      !!-------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy * thickness and layers interfaces
      !--------------------------------------------------------------------------
      zeh_cum0_2d(0) = 0._wp
      zh_cum0_2d(0)  = 0._wp

      DO jk0 = 1, nlay_i+2
         zeh_cum0_2d(jk0) = zeh_cum0_2d(jk0-1) + peh_i_old(jk0-1)
         zh_cum0_2d(jk0)  = zh_cum0_2d(jk0-1) + ph_i_old(jk0-1)
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      zhnew_2d = SUM( ph_i_old(0:nlay_i+1) ) * r1_nlay_i

      ! new layers interfaces
      zh_cum1_2d(0) = 0._wp

      DO jk1 = 1, nlay_i
         zh_cum1_2d(jk1) = zh_cum1_2d(jk1-1) + zhnew_2d
      END DO

      zeh_cum1_2d(0:nlay_i) = 0._wp

      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_i+2
         DO jk1 = 1, nlay_i-1
            IF( zh_cum1_2d(jk1) <= zh_cum0_2d(jk0) .AND. zh_cum1_2d(jk1) > zh_cum0_2d(jk0-1) ) THEN
               zeh_cum1_2d(jk1) = ( zeh_cum0_2d(jk0-1) * ( zh_cum0_2d(jk0) - zh_cum1_2d(jk1  ) ) +  &
                  &                 zeh_cum0_2d(jk0  ) * ( zh_cum1_2d(jk1) - zh_cum0_2d(jk0-1) ) )  &
                  &             / ( zh_cum0_2d(jk0) - zh_cum0_2d(jk0-1) )
            ENDIF
         END DO
      END DO
      ! to ensure that total heat content is strictly conserved, set:
      zeh_cum1_2d(nlay_i) = zeh_cum0_2d(nlay_i+2)

      ! new enthalpies
      DO jk1 = 1, nlay_i
         zswitch      = MAX( 0._wp , SIGN( 1._wp , zhnew_2d - epsi20 ) )
         qnew(jk1) = zswitch * MAX( 0._wp, zeh_cum1_2d(jk1) - zeh_cum1_2d(jk1-1) ) / MAX( zhnew_2d, epsi20 ) ! max for roundoff error
      END DO

   END SUBROUTINE ice_thd_ent_scl


   !!======================================================================
END MODULE icethd_ent
