MODULE phycst
   !!======================================================================
   !!                    ***  MODULE  phycst  ***
   !!     Definition of of both ocean and ice parameters used in the code
   !!=====================================================================
   !! History :   OPA  !  1990-10  (C. Levy - G. Madec)  Original code
   !!             8.1  !  1991-11  (G. Madec, M. Imbard)  cosmetic changes
   !!   NEMO      1.0  !  2002-08  (G. Madec, C. Ethe)  F90, add ice constants
   !!              -   !  2006-08  (G. Madec)  style
   !!             3.2  !  2006-08  (S. Masson, G. Madec)  suppress useless variables + style
   !!             3.4  !  2011-11  (C. Harris)  minor changes for CICE constants
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   phy_cst  : define and print physical constant and domain parameters
   !!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   phy_cst     ! routine called by inipar.F90

   REAL(wp), PARAMETER, PUBLIC ::   rpi      = 3.141592653589793_wp             !: pi
   REAL(wp), PARAMETER, PUBLIC ::   rad      = 3.141592653589793_wp / 180._wp   !: conversion from degre into radian
   REAL(wp), PARAMETER, PUBLIC ::   rsmall   = 0.5 * EPSILON( 1.e0 )            !: smallest real computer value

   REAL(wp), PARAMETER, PUBLIC ::   rday     = 24.*60.*60.      !: day                                [s]
   REAL(wp), PUBLIC ::   rsiyea                      !: sideral year                       [s]
   REAL(wp), PUBLIC ::   rsiday                      !: sideral day                        [s]
   REAL(wp), PARAMETER, PUBLIC ::   raamo    =  12._wp          !: number of months in one year
   REAL(wp), PARAMETER, PUBLIC ::   rjjhh    =  24._wp          !: number of hours in one day
   REAL(wp), PARAMETER, PUBLIC ::   rhhmm    =  60._wp          !: number of minutes in one hour
   REAL(wp), PARAMETER, PUBLIC ::   rmmss    =  60._wp          !: number of seconds in one minute
   REAL(wp), PUBLIC ::   omega                       !: earth rotation parameter           [s-1]
   REAL(wp), PARAMETER, PUBLIC ::   ra       = 6371229._wp      !: earth radius                       [m]
   REAL(wp), PARAMETER, PUBLIC ::   grav     = 9.80665_wp       !: gravity                            [m/s2]
   REAL(wp), PARAMETER, PUBLIC ::   rt0      = 273.15_wp        !: freezing point of fresh water [Kelvin]
   !$acc declare create( rsiyea, rsiday, omega, grav, rt0 )
   
   REAL(wp), PUBLIC ::   rho0                        !: volumic mass of reference     [kg/m3]
   REAL(wp), PUBLIC ::   r1_rho0                     !: = 1. / rho0                   [m3/kg]
   REAL(wp), PUBLIC ::   rcp                         !: ocean specific heat           [J/Kelvin/kg]
   REAL(wp), PUBLIC ::   r1_rcp                      !: = 1. / rcp                    [Kelvin/J]
   REAL(wp), PUBLIC ::   rho0_rcp                    !: = rho0 * rcp                  [J/K/m3)]
   REAL(wp), PUBLIC ::   r1_rho0_rcp                 !: = 1. / ( rho0 * rcp )
   !$acc declare create(  rho0, r1_rho0, rcp, r1_rcp, rho0_rcp, r1_rho0_rcp )

   REAL(wp), PARAMETER, PUBLIC ::   emic     =    0.97_wp       !: emissivity of snow or ice (not used?)

   REAL(wp), PARAMETER, PUBLIC ::   sice     =    6.0_wp        !: salinity of ice (for pisces)          [psu]
   REAL(wp), PARAMETER, PUBLIC ::   soce     =   34.7_wp        !: salinity of sea (for pisces and isf)  [psu]
   REAL(wp), PARAMETER, PUBLIC ::   rLevap   =    2.5e+6_wp     !: latent heat of evaporation (water)
   REAL(wp), PARAMETER, PUBLIC ::   vkarmn   =    0.4_wp        !: von Karman constant
   REAL(wp), PARAMETER, PUBLIC ::   vkarmn2  =    0.4_wp*0.4_wp !: square of von Karman constant
   REAL(wp), PARAMETER, PUBLIC ::   stefan   =    5.67e-8_wp    !: Stefan-Boltzmann constant
   !$acc declare create( emic, sice, soce, rLevap, vkarmn, vkarmn2, stefan )
   
   REAL(wp), PARAMETER, PUBLIC ::   rhos     =  330._wp         !: volumic mass of snow                                  [kg/m3]
   REAL(wp), PARAMETER, PUBLIC ::   rhoi     =  917._wp         !: volumic mass of sea ice                               [kg/m3]
   REAL(wp), PARAMETER, PUBLIC ::   rhow     = 1000._wp         !: volumic mass of freshwater in melt ponds              [kg/m3]
   REAL(wp), PARAMETER, PUBLIC ::   rcnd_i   = 2.034396_wp      !: thermal conductivity of fresh ice                     [W/m/K]
   REAL(wp),            PUBLIC ::   rcnd_s   = 0.31_wp          !: thermal conductivity of snow (0.31 W/m/K, Maykut and Untersteiner, 1971) [W/m/K]
   REAL(wp), PARAMETER, PUBLIC ::   rcpi     = 2067.0_wp        !: specific heat of fresh ice                            [J/kg/K]
   REAL(wp), PARAMETER, PUBLIC ::   rLsub    =    2.834e+6_wp   !: pure ice latent heat of sublimation                   [J/kg]
   REAL(wp), PARAMETER, PUBLIC ::   rLfus    =    0.334e+6_wp   !: latent heat of fusion of fresh ice                    [J/kg]
   REAL(wp), PARAMETER, PUBLIC ::   rTmlt    =    0.054_wp      !: decrease of seawater meltpoint with salinity
   !$acc declare create( rhos, rhoi, rhow, rcnd_i, rcnd_s, rcpi, rLsub, rLfus, rTmlt )

   REAL(wp), PUBLIC ::   r1_rhoi                     !: 1 / rhoi
   REAL(wp), PUBLIC ::   r1_rhos                     !: 1 / rhos
   REAL(wp), PUBLIC ::   r1_rcpi                     !: 1 / rcpi

   REAL(wp), PARAMETER, PUBLIC :: rnup   = 1._wp/3._wp  !: Poisson's ratio
   REAL(wp), PARAMETER, PUBLIC :: rmuMC  = 0.7_wp       !: slope of Mohr-Coulomb enveloppe
   REAL(wp), PARAMETER, PUBLIC :: rsqrt_nu_rhoi = SQRT( 2._wp*(1._wp + rnup)*rhoi )
   
   REAL(wp), PUBLIC, PARAMETER :: ref_tau_max = 10._wp ! Wind stress [N/m2]
   !$acc declare create( r1_rhoi, r1_rhos, r1_rcpi, rnup, rmuMC, rsqrt_nu_rhoi, ref_tau_max )
   
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! $Id: phycst.F90 14072 2020-12-04 07:48:38Z laurent $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE phy_cst
      !!----------------------------------------------------------------------
      !!                       ***  ROUTINE phy_cst  ***
      !!
      !! ** Purpose :   set and print the constants
      !!----------------------------------------------------------------------

      rsiyea = 365.25_wp * rday * 2._wp * rpi / 6.283076_wp
      rsiday = rday / ( 1._wp + rday / rsiyea )
      omega  = 2._wp * rpi / rsiday

      r1_rhoi = 1._wp / rhoi
      r1_rhos = 1._wp / rhos
      r1_rcpi = 1._wp / rcpi

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'phy_cst : initialization of ocean parameters and constants'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '      mathematical constant                 rpi = ', rpi
         WRITE(numout,*) '      day                                rday   = ', rday,   ' s'
         WRITE(numout,*) '      sideral year                       rsiyea = ', rsiyea, ' s'
         WRITE(numout,*) '      sideral day                        rsiday = ', rsiday, ' s'
         WRITE(numout,*) '      omega                              omega  = ', omega,  ' s^-1'
         WRITE(numout,*)
         WRITE(numout,*) '      nb of months per year               raamo = ', raamo, ' months'
         WRITE(numout,*) '      nb of hours per day                 rjjhh = ', rjjhh, ' hours'
         WRITE(numout,*) '      nb of minutes per hour              rhhmm = ', rhhmm, ' mn'
         WRITE(numout,*) '      nb of seconds per minute            rmmss = ', rmmss, ' s'
         WRITE(numout,*)
         WRITE(numout,*) '      earth radius                         ra   = ', ra   , ' m'
         WRITE(numout,*) '      gravity                              grav = ', grav , ' m/s^2'
         WRITE(numout,*)
         WRITE(numout,*) '      freezing point of water              rt0  = ', rt0  , ' K'
         WRITE(numout,*)
         WRITE(numout,*) '   reference density and heat capacity now defined in eosbn2.f90'
         WRITE(numout,*)
         WRITE(numout,*) '      thermal conductivity of pure ice          = ', rcnd_i  , ' J/s/m/K'
         WRITE(numout,*) '      thermal conductivity of snow is defined in a namelist '
         WRITE(numout,*) '      fresh ice specific heat                   = ', rcpi    , ' J/kg/K'
         WRITE(numout,*) '      latent heat of fusion of fresh ice / snow = ', rLfus   , ' J/kg'
         WRITE(numout,*) '      latent heat of subl.  of fresh ice / snow = ', rLsub   , ' J/kg'
         WRITE(numout,*) '      density of sea ice                        = ', rhoi    , ' kg/m^3'
         WRITE(numout,*) '      density of snow                           = ', rhos    , ' kg/m^3'
         WRITE(numout,*) '      density of freshwater (in melt ponds)     = ', rhow    , ' kg/m^3'
         WRITE(numout,*) '      salinity of ice (for pisces)              = ', sice    , ' psu'
         WRITE(numout,*) '      salinity of sea (for pisces and isf)      = ', soce    , ' psu'
         WRITE(numout,*) '      latent heat of evaporation (water)        = ', rLevap  , ' J/m^3'
         WRITE(numout,*) '      von Karman constant                       = ', vkarmn
         WRITE(numout,*) '      Stefan-Boltzmann constant                 = ', stefan  , ' J/s/m^2/K^4'
         WRITE(numout,*)
         WRITE(numout,*) '      conversion: degre ==> radian          rad = ', rad
         WRITE(numout,*)
         WRITE(numout,*) '      smallest real computer value       rsmall = ', rsmall
      ENDIF


      !$acc update device ( rsiyea, rsiday, omega, r1_rhoi, r1_rhos, r1_rcpi )
      
   END SUBROUTINE phy_cst

   !!======================================================================
END MODULE phycst
