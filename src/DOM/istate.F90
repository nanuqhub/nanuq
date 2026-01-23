MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec) Merge of dtatem and dtasal & suppression of tb,tn/sb,sn 
   !!            3.7  !  2016-04  (S. Flavoni) introduce user defined initial state 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   !USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain 
   USE daymod         ! calendar
   !USE dtatsd         ! data temperature and salinity   (dta_tsd routine)
   !USE dtauvd         ! data: U & V current             (dta_uvd routine)
   !USE domvvl          ! varying vertical mesh
   !USE wet_dry         ! wetting and drying (needed for wad_istate)
   !USE usrdef_istate   ! User defined initial state
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE lbclnk         ! lateal boundary condition / mpp exchanges
   !USE restart         ! restart

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by nanuqgcm.F90

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: istate.F90 15052 2021-06-24 14:39:14Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::  Kbb, Kmm, Kaa   ! ocean time level indices
      !
      INTEGER ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_init : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      CALL dta_tsd_init                 ! Initialisation of T & S input data

         IF( ln_rstart ) THEN                    ! Restart from a file
            !                                    ! -------------------
            !CALL rst_read( Kbb, Kmm )            ! Read the restart file !#lolo: commented because ocean restarts!?
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !
         ELSE                                    ! Start from rest
            !                                    ! ---------------
            numror = 0                           ! define numror = 0 -> no restart file to read
            !
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !                                    ! Initialization of ocean to zero
            !
         ENDIF 
      ! 
   END SUBROUTINE istate_init

   !!======================================================================
END MODULE istate
