MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of kbat and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.6  ! 2015-05  (P. Mathiot) ISF: add wmask,wumask and wvmask
   !!            4.0  ! 2016-06  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!            4.x  ! 2020-02  (G. Madec, S. Techene) introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk       : compute land/ocean mask
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE domutl, ONLY : dom_uniq
   USE bdy        ! open boundary
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! Massively Parallel Processing library
   !
   !USE io_ezcdf , ONLY : DUMP_FIELD ; !#LOLOdebug

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk    ! routine called by inidom.F90

   LOGICAL, PARAMETER :: ldebug = .FALSE.

   !! * Substitutions
#  include "read_nml_substitute.h90"

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: dommsk.F90 15556 2021-11-29 15:23:06Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE dom_msk( pbathy )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) points.
      !!
      !! ** Method  :   The ocean/land mask  at t-point is deduced from ko_top
      !!      and ko_bot, the indices of the fist and last ocean t-levels which
      !!      are either defined in usrdef_lsm or read in lsm_read.
      !!                The velocity masks (umask, vmask, wmask, wumask, wvmask)
      !!      are deduced from a product of the two neighboring tmask.
      !!
      !! ** Action :   tmask, umask, vmask, wmask, wumask, wvmask : land/ocean mask
      !!                         at t-, u-, v- w, wu-, and wv-points (=0. or 1.)
      !!               tmask_i : ssmask * ( excludes halo+duplicated points (NP folding) )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) :: pbathy
      !
      REAL(wp),   DIMENSION(jpi,jpj) :: z2dt, z2df, z2du
      INTEGER(1), DIMENSION(jpi,jpj) :: i2dt, i2df, i2du
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   ios, inum
      !
      CHARACTER(len=64) :: cf_tmp ; !#LOLOdebug
      !!
      NAMELIST/nambdy/ ln_bdy ,nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice, nn_ice_dta, cn_dmg, nn_dmg_dta, nn_rimwidth
      !!---------------------------------------------------------------------

      !  Ocean/land mask at t-point  (computed from ko_top and ko_bot)
      ! ----------------------------
      !
      tmask(:,:,:) = 0._wp
      WHERE( pbathy(:,:) > 0.01_wp ) tmask(:,:,1) = 1._wp

      ! Mask corrections for bdy (read in mppini2)
      READ_NML_REF(numnam,nambdy)
      !903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist' )
      READ_NML_CFG(numnam,nambdy)
      !904   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist' )
      ! ------------------------
      IF ( ln_bdy .AND. ln_mask_file ) THEN
         CALL iom_open( cn_mask_file, inum )
         CALL iom_get ( inum, jpdom_global, 'bdy_msk', bdytmask(:,:) )
         CALL iom_close( inum )
         DO jj=Njs0-1, Nje0+1
            DO ji=Nis0-1, Nie0+1
               tmask(ji,jj,1) = tmask(ji,jj,1) * bdytmask(ji,jj)
            END DO
         ENDDO
      ENDIF

      !  Ocean/land mask at u-, v-, and f-points   (computed from tmask)
      ! ----------------------------------------
      DO jj=Njs0, Nje0
         DO ji=Nis0, Nie0
            umask(ji,jj,1) = tmask(ji,jj  ,1) * tmask(ji+1,jj  ,1)
            vmask(ji,jj,1) = tmask(ji,jj  ,1) * tmask(ji  ,jj+1,1)
         END DO
      ENDDO
      !
      CALL lbc_lnk( 'dommsk', umask, 'U', 1.0_wp, vmask, 'V', 1.0_wp )

      ! Interior domain mask  (used for global sum) : 2D ocean mask x (halo+duplicated points) mask
      ! --------------------
      !
      CALL dom_uniq( tmask_i, 'T' )
      tmask_i(:,:) = MAXVAL( tmask(:,:,:), DIM=3 ) * tmask_i(:,:) ; ! ssmask = MAXVAL( tmask(:,:,:), DIM=3 )

      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------
      ! ==> this is now done in `icedyn.F90` because `rn_ishlat` has yet to be known!


      !! Masks at T- & F-points for the E-grid-based brittle rheology approach
      !! *********************************************************************
      !!  => for F-points it's tricky, because at a solid boundary, the F-point is at the interface between
      !!     ocean and land, and as such, it can be considered both as ocean and land depending on the variable
      !!     of interest. So far, given variables that depend on it, it is considered as a wet point !!!
      !!     ==>  `xmskf==1` wherever the F-point is in contact with water (=> includes the land-water interface)
      !!  As opposed to that the mask `fmask` (constructed in `icedyn.F90`) is the one used to compute the velocity
      !!  shear, based on `rn_ishlat` the type of expected boundaries at the coast (slip / no-slip).
      !!
      !! Examples:
      !!    * set ice concentration at F-points over land => USE `xmskf` !
      !!    * compute the shear of ice velocity vector => USE `fmask` !

      xmskt(:,:) = tmask(:,:,1)
      xmskf(:,:) = 0._wp
      i2dt(:,:) = INT( tmask(:,:,1), 1 )
      xmskf(Nis0:Nie0,Njs0:Nje0) = REAL( MIN( i2dt(Nis0:Nie0,Njs0:Nje0) + i2dt(Nis0+1:Nie0+1,Njs0+1:Nje0+1) &
         &                       + i2dt(Nis0+1:Nie0+1,Njs0:Nje0) + i2dt(Nis0:Nie0,Njs0+1:Nje0+1), 1 ), wp )

      !! Same spirit as for `xmskt`: `xmsku` is the mask at U-points for tracers that are not velocities and can therefore
      !! have a value other than 0. at the solid interface...
      xmsku(:,:) = 0._wp
      xmsku(Nis0:Nie0,Njs0:Nje0) = REAL( MIN( i2dt(Nis0:Nie0,Njs0:Nje0) + i2dt(Nis0+1:Nie0+1,Njs0:Nje0), 1 ), wp )

      CALL lbc_lnk( 'dommsk',  xmskt,'T',1._wp, xmskf,'F',1._wp, xmsku,'U',1._wp ) !LOLOlbclnk, probably not needed for `xmskt`

      !IF(ldebug) THEN
      !   WRITE(cf_tmp,'("xmskt_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(xmskt(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'xmskt' )
      !   WRITE(cf_tmp,'("xmskf_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(xmskf(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'xmskf' )
      !   WRITE(cf_tmp,'("xmsku_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(xmsku(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'xmsku' )
      !ENDIF

      r1_e1e2t(:,:) = r1_e1e2t(:,:) * xmskt(:,:)
      r1_e1e2f(:,:) = r1_e1e2f(:,:) * xmskf(:,:)



      !! Mask for solid lateral boundary conditions
      !! ------------------------------------------

      i2dt(:,:) = 0
      i2dt(:,:) = 1 - INT(xmskt(:,:),1)
      i2df(:,:) = 0
      i2df(:,:) = 1 - INT(xmskf(:,:),1)
      i2du(:,:) = 0
      i2du(:,:) = 1 - INT(xmsku(:,:),1)

      !! #LOLOfixme: add `ln_damage` test or equivalent to prevent working with the F-arrays if they are not needed!

      klbct(:,:,:,:) = 0
      klbcf(:,:,:,:) = 0
      klbcu(:,:,:,:) = 0

      !! *** Eastern solid boundary ***
      z2dt(:,:) = 0._wp ; z2df(:,:) = 0._wp ; z2du(:,:) = 0._wp
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nie0, Nis0, -1
            IF( (i2dt(ji,jj)==0).AND.(i2dt(ji+1,jj)==1) )                                  z2dt(ji,jj) = 1.
            IF( (i2df(ji,jj)==0).AND.(i2df(ji+1,jj)==1) )                                  z2df(ji,jj) = 1.
            IF( (i2du(ji,jj)==0).AND.(i2du(ji+1,jj)==1) )                                  z2du(ji,jj) = 1.
            IF( nn_hls>1 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji+2,jj)==1).AND.(    i2dt(ji+1     ,jj) ==0) ) z2dt(ji,jj) = 2.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji+2,jj)==1).AND.(    i2df(ji+1     ,jj) ==0) ) z2df(ji,jj) = 2.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji+2,jj)==1).AND.(    i2du(ji+1     ,jj) ==0) ) z2du(ji,jj) = 2.
            ENDIF
            IF( nn_hls>2 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji+3,jj)==1).AND.(SUM(i2dt(ji+1:ji+2,jj))==0) ) z2dt(ji,jj) = 3.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji+3,jj)==1).AND.(SUM(i2df(ji+1:ji+2,jj))==0) ) z2df(ji,jj) = 3.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji+3,jj)==1).AND.(SUM(i2du(ji+1:ji+2,jj))==0) ) z2du(ji,jj) = 3.
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'dommsk', z2dt,'T',1._wp, z2df,'F',1._wp, z2du,'U',1._wp )
      WHERE(    INT(z2dt(:,:),1) == 1 ) klbct(:,:,1,1) = 1 ! wet T-point at i & the T-point at i+1 is a land point
      WHERE(    INT(z2df(:,:),1) == 1 ) klbcf(:,:,1,1) = 1 ! wet F-point at i & the F-point at i+1 is a land point
      WHERE(    INT(z2du(:,:),1) == 1 ) klbcu(:,:,1,1) = 1
      IF( nn_hls > 1 ) THEN
         WHERE( INT(z2dt(:,:),1) == 2 ) klbct(:,:,2,1) = 1 ! wet T-point at i and i+1 & the T-point at i+2 is a land point
         WHERE( INT(z2df(:,:),1) == 2 ) klbcf(:,:,2,1) = 1 ! wet F-point at i and i+1 & the F-point at i+2 is a land point
         WHERE( INT(z2du(:,:),1) == 2 ) klbcu(:,:,2,1) = 1
      ENDIF
      IF( nn_hls > 2 ) THEN
         WHERE( INT(z2dt(:,:),1) == 3 ) klbct(:,:,3,1) = 1 ! wet T-point at i, i+1, and i+2 & the T-point at i+3 is a land point
         WHERE( INT(z2df(:,:),1) == 3 ) klbcf(:,:,3,1) = 1 ! wet F-point at i, i+1, and i+2 & the F-point at i+3 is a land point
         WHERE( INT(z2du(:,:),1) == 3 ) klbcu(:,:,3,1) = 1
      ENDIF
      !
      !IF(ldebug) THEN
      !   WRITE(cf_tmp,'("kmEt_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2dt(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmEt' )
      !   WRITE(cf_tmp,'("kmEf_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2df(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmEf' )
      !   WRITE(cf_tmp,'("kmEu_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2du(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmEu' )
      !ENDIF

      !! *** Northern solid boundary ***
      z2dt(:,:) = 0._wp ; z2df(:,:) = 0._wp ; z2du(:,:) = 0._wp
      DO jj=Nje0, Njs0, -1
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF(    (i2dt(ji,jj)==0).AND.(i2dt(ji,jj+1)==1) )                                  z2dt(ji,jj) = 1.
            IF(    (i2df(ji,jj)==0).AND.(i2df(ji,jj+1)==1) )                                  z2df(ji,jj) = 1.
            IF(    (i2du(ji,jj)==0).AND.(i2du(ji,jj+1)==1) )                                  z2du(ji,jj) = 1.
            IF( nn_hls>1 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji,jj+2)==1).AND.(    i2dt(ji,jj+1)      ==0) ) z2dt(ji,jj) = 2.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji,jj+2)==1).AND.(    i2df(ji,jj+1)      ==0) ) z2df(ji,jj) = 2.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji,jj+2)==1).AND.(    i2du(ji,jj+1)      ==0) ) z2du(ji,jj) = 2.
            ENDIF
            IF( nn_hls>2 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji,jj+3)==1).AND.(SUM(i2dt(ji,jj+1:jj+2))==0) ) z2dt(ji,jj) = 3.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji,jj+3)==1).AND.(SUM(i2df(ji,jj+1:jj+2))==0) ) z2df(ji,jj) = 3.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji,jj+3)==1).AND.(SUM(i2du(ji,jj+1:jj+2))==0) ) z2du(ji,jj) = 3.
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'dommsk', z2dt,'T',1._wp, z2df,'F',1._wp, z2du,'U',1._wp )
      WHERE(    INT(z2dt(:,:),1) == 1 ) klbct(:,:,1,2) = 1 ! wet T-point at j & the T-point at j+1 is a land point
      WHERE(    INT(z2df(:,:),1) == 1 ) klbcf(:,:,1,2) = 1 ! wet F-point at j & the F-point at j+1 is a land point
      WHERE(    INT(z2du(:,:),1) == 1 ) klbcu(:,:,1,2) = 1 ! wet F-point at j & the F-point at j+1 is a land point
      IF( nn_hls > 1 ) THEN
         WHERE( INT(z2dt(:,:),1) == 2 ) klbct(:,:,2,2) = 1 ! wet T-point at j and j+1 & the T-point at j+2 is a land point
         WHERE( INT(z2df(:,:),1) == 2 ) klbcf(:,:,2,2) = 1 ! wet F-point at j and j+1 & the F-point at j+2 is a land point
         WHERE( INT(z2du(:,:),1) == 2 ) klbcu(:,:,2,2) = 1 ! wet F-point at j and j+1 & the F-point at j+2 is a land point
      ENDIF
      IF( nn_hls > 2 ) THEN
         WHERE( INT(z2dt(:,:),1) == 3 ) klbct(:,:,3,2) = 1 ! wet T-point at j, j+1, and j+2 & the T-point at j+3 is a land point
         WHERE( INT(z2df(:,:),1) == 3 ) klbcf(:,:,3,2) = 1 ! wet F-point at j, j+1, and j+2 & the F-point at j+3 is a land point
         WHERE( INT(z2du(:,:),1) == 3 ) klbcu(:,:,3,2) = 1 ! wet F-point at j, j+1, and j+2 & the F-point at j+3 is a land point
      ENDIF
      !
      !IF(ldebug) THEN
      !   WRITE(cf_tmp,'("kmNt_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2dt(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmNt' )
      !   WRITE(cf_tmp,'("kmNf_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2df(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmNf' )
      !   WRITE(cf_tmp,'("kmNu_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2du(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmNu' )
      !ENDIF

      !! *** Western solid boundary ***
      z2dt(:,:) = 0._wp ; z2df(:,:) = 0._wp ; z2du(:,:) = 0._wp
      DO jj=Njs0-nn_hls, Nje0+nn_hls
         DO ji=Nis0, Nie0
            IF(    (i2dt(ji,jj)==0).AND.(i2dt(ji-1,jj)==1) )                                  z2dt(ji,jj) = 1.
            IF(    (i2df(ji,jj)==0).AND.(i2df(ji-1,jj)==1) )                                  z2df(ji,jj) = 1.
            IF(    (i2du(ji,jj)==0).AND.(i2du(ji-1,jj)==1) )                                  z2du(ji,jj) = 1.
            IF( nn_hls>1 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji-2,jj)==1).AND.(    i2dt(     ji-1,jj) ==0) ) z2dt(ji,jj) = 2.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji-2,jj)==1).AND.(    i2df(     ji-1,jj) ==0) ) z2df(ji,jj) = 2.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji-2,jj)==1).AND.(    i2du(     ji-1,jj) ==0) ) z2du(ji,jj) = 2.
            ENDIF
            IF( nn_hls>2 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji-3,jj)==1).AND.(SUM(i2dt(ji-2:ji-1,jj))==0) ) z2dt(ji,jj) = 3.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji-3,jj)==1).AND.(SUM(i2df(ji-2:ji-1,jj))==0) ) z2df(ji,jj) = 3.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji-3,jj)==1).AND.(SUM(i2du(ji-2:ji-1,jj))==0) ) z2du(ji,jj) = 3.
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'dommsk', z2dt,'T',1._wp, z2df,'F',1._wp, z2du,'U',1._wp )
      WHERE(    INT(z2dt(:,:),1) == 1 ) klbct(:,:,1,3) = 1 ! wet T-point at i & the T-point at i-1 is a land point
      WHERE(    INT(z2df(:,:),1) == 1 ) klbcf(:,:,1,3) = 1 ! wet F-point at i & the F-point at i-1 is a land point
      WHERE(    INT(z2du(:,:),1) == 1 ) klbcu(:,:,1,3) = 1
      IF( nn_hls > 1 ) THEN
         WHERE( INT(z2dt(:,:),1) == 2 ) klbct(:,:,2,3) = 1 ! wet T-point at i and i-1 & the T-point at i-2 is a land point
         WHERE( INT(z2df(:,:),1) == 2 ) klbcf(:,:,2,3) = 1 ! wet F-point at i and i-1 & the F-point at i-2 is a land point
         WHERE( INT(z2du(:,:),1) == 2 ) klbcu(:,:,2,3) = 1
      ENDIF
      IF( nn_hls > 2 ) THEN
         WHERE( INT(z2dt(:,:),1) == 3 ) klbct(:,:,3,3) = 1 ! wet T-point at i, i-1, and i-2 & the T-point at i-3 is a land point
         WHERE( INT(z2df(:,:),1) == 3 ) klbcf(:,:,3,3) = 1 ! wet F-point at i, i-1, and i-2 & the F-point at i-3 is a land point
         WHERE( INT(z2du(:,:),1) == 3 ) klbcu(:,:,3,3) = 1
      ENDIF
      !
      !IF(ldebug) THEN
      !   WRITE(cf_tmp,'("kmWt_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2dt(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmWt' )
      !   WRITE(cf_tmp,'("kmWf_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2df(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmWf' )
      !   WRITE(cf_tmp,'("kmWu_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2du(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmWu' )
      !ENDIF

      !! *** Southern solid boundary ***
      z2dt(:,:) = 0._wp ; z2df(:,:) = 0._wp ; z2du(:,:) = 0._wp
      DO jj=Njs0, Nje0
         DO ji=Nis0-nn_hls, Nie0+nn_hls
            IF(    (i2dt(ji,jj)==0).AND.(i2dt(ji,jj-1)==1) )                                  z2dt(ji,jj) = 1.
            IF(    (i2df(ji,jj)==0).AND.(i2df(ji,jj-1)==1) )                                  z2df(ji,jj) = 1.
            IF(    (i2du(ji,jj)==0).AND.(i2du(ji,jj-1)==1) )                                  z2du(ji,jj) = 1.
            IF( nn_hls>1 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji,jj-2)==1).AND.(    i2dt(ji,     jj-1) ==0) ) z2dt(ji,jj) = 2.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji,jj-2)==1).AND.(    i2df(ji,     jj-1) ==0) ) z2df(ji,jj) = 2.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji,jj-2)==1).AND.(    i2du(ji,     jj-1) ==0) ) z2du(ji,jj) = 2.
            ENDIF
            IF( nn_hls>2 ) THEN
               IF( (i2dt(ji,jj)==0).AND.(i2dt(ji,jj-3)==1).AND.(SUM(i2dt(ji,jj-2:jj-1))==0) ) z2dt(ji,jj) = 3.
               IF( (i2df(ji,jj)==0).AND.(i2df(ji,jj-3)==1).AND.(SUM(i2df(ji,jj-2:jj-1))==0) ) z2df(ji,jj) = 3.
               IF( (i2du(ji,jj)==0).AND.(i2du(ji,jj-3)==1).AND.(SUM(i2du(ji,jj-2:jj-1))==0) ) z2du(ji,jj) = 3.
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'dommsk', z2dt,'T',1._wp, z2df,'F',1._wp, z2du,'U',1._wp )
      WHERE(    INT(z2dt(:,:),1) == 1 ) klbct(:,:,1,4) = 1 ! wet T-point at j & the T-point at j-1 is a land point
      WHERE(    INT(z2df(:,:),1) == 1 ) klbcf(:,:,1,4) = 1 ! wet F-point at j & the F-point at j-1 is a land point
      WHERE(    INT(z2du(:,:),1) == 1 ) klbcu(:,:,1,4) = 1
      IF( nn_hls > 1 ) THEN
         WHERE( INT(z2dt(:,:),1) == 2 ) klbct(:,:,2,4) = 1 ! wet T-point at j and j-1 & the T-point at j-2 is a land point
         WHERE( INT(z2df(:,:),1) == 2 ) klbcf(:,:,2,4) = 1 ! wet F-point at j and j-1 & the F-point at j-2 is a land point
         WHERE( INT(z2du(:,:),1) == 2 ) klbcu(:,:,2,4) = 1
      ENDIF
      IF( nn_hls > 2 ) THEN
         WHERE( INT(z2dt(:,:),1) == 3 ) klbct(:,:,3,4) = 1 ! wet T-point at j, j-1, and j-2 & the T-point at j-3 is a land point
         WHERE( INT(z2df(:,:),1) == 3 ) klbcf(:,:,3,4) = 1 ! wet F-point at j, j-1, and j-2 & the F-point at j-3 is a land point
         WHERE( INT(z2du(:,:),1) == 3 ) klbcu(:,:,3,4) = 1
      ENDIF
      !
      !IF(ldebug) THEN
      !   WRITE(cf_tmp,'("kmSt_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2dt(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmSt' )
      !   WRITE(cf_tmp,'("kmSf_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2df(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmSf' )
      !   WRITE(cf_tmp,'("kmSu_",i4.4,".tmp")') narea
      !   CALL DUMP_FIELD( REAL(z2du(Nis0:Nie0,Njs0:Nje0),4), TRIM(cf_tmp), 'kmSu' )
      !ENDIF

# if defined _OPENACC
      PRINT *, ' * info GPU: dom_msk() => adding `umask,vmask,xmskt,xmskf,klbct,klbcf,klbcu` arrays to memory!'
      !$acc enter data copyin( umask,vmask, xmskt,xmskf, klbct,klbcf,klbcu )
# endif

   END SUBROUTINE dom_msk




   !!======================================================================
END MODULE dommsk
