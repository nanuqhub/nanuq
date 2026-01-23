MODULE dombth
   !!==============================================================================
   !!                       ***  MODULE dombth   ***
   !! Ocean domain : definition of the vertical coordinate system
   !!==============================================================================
   !! History :  OPA  ! 1995-12  (G. Madec)  Original code : s vertical coordinate
   !!                 ! 1997-07  (G. Madec)  lbc_lnk call
   !!                 ! 1997-04  (J.-O. Beismann) 
   !!            8.5  ! 2002-09  (A. Bozec, G. Madec)  F90: Free form and module
   !!             -   ! 2002-09  (A. de Miranda)  rigid-lid + islands
   !!  NEMO      1.0  ! 2003-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-10  (A. Beckmann)  modifications for hybrid s-ccordinates & new stretching function
   !!            2.0  ! 2006-04  (R. Benshila, G. Madec)  add zgr_zco
   !!            3.0  ! 2008-06  (G. Madec)  insertion of domzgr_zps.h90 & conding style
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2012-08  (J. Siddorn) added Siddorn and Furner stretching function
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  modify C1D case  
   !!            3.6  ! 2014-11  (P. Mathiot and C. Harris) add ice shelf capabilitye  
   !!            3.?  ! 2015-11  (H. Liu) Modifications for Wetting/Drying
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_bth       : read or set the ocean vertical coordinate system
   !!   bth_read      : read the vertical information in the domain configuration file
   !!---------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_bth   ! called by dom_init.F90

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: dombth.F90 15556 2021-11-29 15:23:06Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS       

   SUBROUTINE bathy_read( pbathy )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE bathy_read  ***
      !!
      !! ** Purpose :   Read the vertical information in the domain configuration file
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pbathy
      !
      INTEGER  ::   inum
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   bathy_read : read the bathymetry in ', TRIM( cn_domcfg ), ' file'
         WRITE(numout,*) '   ~~~~~~~~'
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      !                          !* ocean top and bottom level
      CALL iom_get( inum, jpdom_global, 'bathy_metry'    , z2d   )
      pbathy(:,:) = NINT( z2d(:,:) )
      !
      CALL iom_close( inum )
      !
   END SUBROUTINE bathy_read


   SUBROUTINE dom_bth( pbathy )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_bth  ***
      !!                   
      !! ** Purpose :   set the depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  : - reference 1D vertical coordinate (gdep._1d, e3._1d)
      !!              - read/set ocean depth and ocean levels (bathy, mbathy)
      !!              - vertical coordinate (gdep., e3.) depending on the 
      !!                coordinate chosen :
      !!                   ln_zco=T   z-coordinate   
      !!                   ln_zps=T   z-coordinate with partial steps
      !!                   ln_zco=T   s-coordinate 
      !!
      !! ** Action  :   define gdep., e3., mbathy and bathy
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pbathy
      !
      INTEGER  ::   ji,jj,jk            ! dummy loop index
      INTEGER  ::   ikt, ikb            ! top/bot index
      INTEGER  ::   ioptio, ibat, ios   ! local integer
      INTEGER  ::   is_mbkuvf           ! ==0 if mbku, mbkv, mbkf to be computed
      !REAL(wp) ::   zrefdep             ! depth of the reference level (~10m)
      REAL(wp), DIMENSION(jpi,jpj  ) ::   zmsk
      REAL(wp), DIMENSION(jpi,jpj,2) ::   ztopbot
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_bth : preparing bth'
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !
      CALL bathy_read( pbathy )
      !
      ! the following is mandatory
      ! make sure that closed boundaries are correctly defined in pbathy that will be used to compute all mask arrays
      !
      zmsk(:,:) = 1._wp                                       ! default: no closed boundaries
      IF( .NOT. l_Iperio ) THEN                                    ! E-W closed:
         zmsk(  mi0(     1+nn_hls):mi1(     1+nn_hls),:) = 0._wp   ! first column of inner global domain at 0
         zmsk(  mi0(jpiglo-nn_hls):mi1(jpiglo-nn_hls),:) = 0._wp   ! last  column of inner global domain at 0 
      ENDIF
      IF( .NOT. l_Jperio ) THEN                                    ! S closed:
         zmsk(:,mj0(     1+nn_hls):mj1(     1+nn_hls)  ) = 0._wp   ! first   line of inner global domain at 0
      ENDIF
      IF( .NOT. ( l_Jperio .OR. l_NFold ) ) THEN                   ! N closed:
         zmsk(:,mj0(jpjglo-nn_hls):mj1(jpjglo-nn_hls)  ) = 0._wp   ! last    line of inner global domain at 0
      ENDIF
      pbathy(:,:) = pbathy(:,:) * NINT( zmsk(:,:) )
      !
      IF( lwp ) WRITE(numout,*) ' MIN val pbathy   ', MINVAL(   pbathy(:,:) ), ' MAX ', MAXVAL( pbathy(:,:) )
      !
   END SUBROUTINE dom_bth
      
   !!======================================================================
END MODULE dombth
