! $Id: ropp_fm_constants.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

MODULE ropp_fm_constants

!****m* Modules/ropp_fm_constants *
!
! NAME
!    ropp_fm_constants - Module providing meteorological and physical constants.
!
! SYNOPSIS
!    use ropp_fm_constants
! 
! DESCRIPTION
!    This module provides numerical values for meteorological and other
!    physical constants which are used throughout the ROPP package.
!
! SEE ALSO
!
!    Thermodynamical constants: R_dry, R_vap, C_p,
!                               mw_dry_air, mw_water, epsilon_water
!
!    Specific constants for RO: kappa1, kappa2, kappa3
!
!    Other physical constants:  g_wmo
!
!    Mathematical constants:    pi
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the GRAS SAF
!   Helpdesk at http://www.grassaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

  USE typesizes, ONLY: wp => EightByteReal

  IMPLICIT NONE

!  private :: wp

!-------------------------------------------------------------------------------
! 1. Thermodynamical Constants
!-------------------------------------------------------------------------------

!****t* Constants/Thermodynamics
!
! DESCRIPTION
!    Thermodynamical constants used in the ROPP Forward Model library.
!
! SEE ALSO
!    R_dry
!    R_vap
!    C_p
!    mw_dry_air
!    mw_water
!    epsilon_water
!
!****

!****p* Thermodynamics/R_dry *
!
! NAME
!    R_dry - Gas constant of dry air
!
! SOURCE
!
  REAL(wp), PARAMETER :: R_dry = 287.0597_wp        !   J K^{-1} kg^{-1}
!
!****

!****p* Thermodynamics/R_vap *
!
! NAME
!    R_vap - Gas constant of water vapor
!
! SOURCE
!
  REAL(wp), PARAMETER :: R_vap = 461.5250_wp        !   J K^{-1} kg^{-1}
!
!****

!****p* Thermodynamics/C_p *
!
! NAME
!    C_p - Specific heat capacity of dry air at constant pressure
!
! SOURCE
!
  REAL(wp), PARAMETER :: C_p = 1005.7_wp       !   J K^{-1} kg^{-1}
!
!****


!****p* Thermodynamics/mw_dry_air *
!
! NAME
!    mw_dry_air - Molecular weight of dry air
!
! SOURCE
!
  REAL(wp), PARAMETER :: mw_dry_air = 28.9648e-3_wp !   kg mol^{-1}
!
!****

!****p* Thermodynamics/mw_water *
!
! NAME
!    mw_water - Molecular weight of water
!
! SOURCE
!
  REAL(wp), PARAMETER :: mw_water = 18.01528e-3_wp  !   kg mol^{-1}
!
!****

!****p* Thermodynamics/epsilon_water *
!
! NAME
!    epsilon_water - Ratio of molecular weight of water to that of dry air.
!
! SOURCE
!
  REAL(wp), PARAMETER :: epsilon_water = mw_water / mw_dry_air
!
!****

!-------------------------------------------------------------------------------
! 2. RO specific constants
!-------------------------------------------------------------------------------

!****t* Constants/RadioOccultation
!
! DESCRIPTION
!    Physical constants specific for radio occultation measurement as used in
!    the ROPP Forward Model library.
!
! SEE ALSO
!    kappa1
!    kappa2
!    kappa3
!
!****

!****p* RadioOccultation/kappa1 *
!
! NAME
!    kappa1 - Coefficient to calculate dry refractivity contribution.
!             N.B. Rueger (2002) value, 77.6890e-2_wp 
!             kappa_comp = (77.689/Compressibilty) at 1013.25 hPa and T=273.15K  
!
! SOURCE
!
  REAL(wp), PARAMETER :: kappa1      = 77.60e-2_wp      !   K Pa^{-1}
  REAL(wp), PARAMETER :: kappa1_comp = 77.643e-2_wp     !   K Pa^{-1}
!
!****

!****p* RadioOccultation/kappa2 *
!
! NAME
!    kappa2 - Coefficient to calculate moist refractivity contribution.
!             N.B. Rueger (2002) value, 3.75463e3_wp
! 
! SOURCE
!
  REAL(wp), PARAMETER :: kappa2      = 3.73e3_wp        !   K^2 Pa^{-1}
  REAL(wp), PARAMETER :: kappa2_comp = 3.75463e3_wp     !   K^2 Pa^{-1}
!
!****

!****p* RadioOccultation/kappa3 *
!
! NAME
!    kappa3 - Coefficient to calculate 2nd moist refractivity contribution.
!             N.B. Rueger (2002) value, -6.3938e-2_wp
!
! SOURCE
!
  REAL(wp), PARAMETER :: kappa3      = 77.60e-2_wp        !   K Pa^{-1}
  REAL(wp), PARAMETER :: kappa3_comp = 71.2952e-2_wp        !   K Pa^{-1}
!
!****

!-------------------------------------------------------------------------------
! 3. Other physical constants
!-------------------------------------------------------------------------------

!****t* Constants/Earth
!
! DESCRIPTION
!    Physical constants related to the Earth as used in the ROPP Forward Model
!    library. 
!
! SEE ALSO
!    g_wmo
!
!****

!****p* Earth/g_wmo *
!
! NAME
!    g_wmo - Gravitational acceleration (WMO standard value)
!
! SOURCE
!
  REAL(wp), PARAMETER :: g_wmo = 9.80665_wp         !   m s^{-2}
!
!****

!-------------------------------------------------------------------------------
! 4. Mathematical constants
!------------------------------------------------------------------------------

!****t* Constants/Maths
!
! DESCRIPTION
!    Mathematical constants as used in the ROPP Forward Model library. 
!
! SEE ALSO
!    pi
!
!****

!****p* Maths/pi *
!
! NAME
!    pi  
!
! SOURCE
!
  REAL(wp), PARAMETER :: pi = 3.141592653589793238_wp   
!
!****

END MODULE ropp_fm_constants
