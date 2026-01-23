MODULE domwri
   !!======================================================================
   !!                       ***  MODULE domwri  ***
   !! Ocean initialization : write the ocean domain mesh file(s)
   !!======================================================================
   !! History :  OPA  ! 1997-02  (G. Madec)  Original code
   !!            8.1  ! 1999-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90 and several file
   !!            3.0  ! 2008-01  (S. Masson)  add dom_uniq
   !!            4.0  ! 2016-01  (G. Madec)  simplified mesh_mask.nc file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_wri        : create and write mesh and mask file(s)
   !!----------------------------------------------------------------------
   !
   USE dom_oce         ! ocean space and time domain
   USE domutl, ONLY : dom_uniq
   USE phycst, ONLY : rsmall
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_wri              ! routine called by inidom.F90

   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: domwri.F90 15033 2021-06-21 10:24:45Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_wri
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_wri  ***
      !!
      !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
      !!      ocean domain informations (mesh and mask arrays). This (these)
      !!      file(s) is (are) used for visualisation (SAXO software) and
      !!      diagnostic computation.
      !!
      !! ** Method  :   create a file with all domain related arrays
      !!
      !! ** output file :   meshmask.nc  : domain size, horizontal grid-point position,
      !!                                   masks, depth and vertical scale factors
      !!----------------------------------------------------------------------
      INTEGER           ::   inum    ! temprary units for 'mesh_mask.nc' file
      CHARACTER(len=21) ::   clnam   ! filename (mesh and mask informations)
      INTEGER           ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj)     ::   zprt, zprw     ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zdepu, zdepv   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dom_wri : create NetCDF mesh and mask information file(s)'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      clnam = 'mesh_mask'  ! filename (mesh and mask informations)

      !                                  ! ============================
      !                                  !  create 'mesh_mask.nc' file
      !                                  ! ============================
      CALL iom_open( TRIM(clnam), inum, ldwrt = .TRUE. )
      !                                                         ! Configuration specificities
      CALL iom_putatt( inum,  'CfgName', TRIM(cn_cfg) )
      CALL iom_putatt( inum, 'CfgIndex',      nn_cfg  )
      !                                                         ! lateral boundary of the global domain
      CALL iom_putatt( inum,   'Iperio', COUNT( (/l_Iperio/) ) )
      CALL iom_putatt( inum,   'Jperio', COUNT( (/l_Jperio/) ) )
      CALL iom_putatt( inum,    'NFold', COUNT( (/l_NFold /) ) )
      CALL iom_putatt( inum,   'NFtype',          c_NFtype     )
      !                                                         ! type of vertical coordinate
      IF(ln_zco)   CALL iom_putatt( inum, 'VertCoord', 'zco' )
      IF(ln_zps)   CALL iom_putatt( inum, 'VertCoord', 'zps' )
      IF(ln_sco)   CALL iom_putatt( inum, 'VertCoord', 'sco' )
      !                                                         ! ocean cavities under iceshelves
      CALL iom_putatt( inum,   'IsfCav', COUNT( (/ln_isfcav/) ) )
      !                                                         ! masks
      CALL iom_rstput( 0, 0, inum, 'tmask', tmask, ktype = jp_i1 )     !    ! land-sea mask
      CALL iom_rstput( 0, 0, inum, 'umask', umask, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum, 'vmask', vmask, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum, 'fmask', fmask, ktype = jp_i1 )

      CALL iom_rstput( 0, 0, inum, 'glamt', glamt, ktype = jp_r8 )     !    ! latitude
      CALL iom_rstput( 0, 0, inum, 'glamu', glamu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamv', glamv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamf', glamf, ktype = jp_r8 )

      CALL iom_rstput( 0, 0, inum, 'gphit', gphit, ktype = jp_r8 )     !    ! longitude
      CALL iom_rstput( 0, 0, inum, 'gphiu', gphiu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphiv', gphiv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphif', gphif, ktype = jp_r8 )

      CALL iom_rstput( 0, 0, inum, 'e1t', e1t, ktype = jp_r8 )         !    ! e1 scale factors
      CALL iom_rstput( 0, 0, inum, 'e1u', e1u, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1v', e1v, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1f', e1f, ktype = jp_r8 )

      CALL iom_rstput( 0, 0, inum, 'e2t', e2t, ktype = jp_r8 )         !    ! e2 scale factors
      CALL iom_rstput( 0, 0, inum, 'e2u', e2u, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2v', e2v, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2f', e2f, ktype = jp_r8 )

      !                                     ! ============================
      CALL iom_close( inum )                !        close the files
      !                                     ! ============================
   END SUBROUTINE dom_wri



   !!======================================================================
END MODULE domwri
