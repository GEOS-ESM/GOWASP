! $Id: ropp_r_eff.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****f* Geodesy/R_eff
!
! NAME
!    R_eff - Effective radius of Earth for geopotential height conversion
!               following Somigliana's equation
!
! SYNOPSIS
!    R_eff = R_eff(lat)
! 
! DESCRIPTION
!    This function calculates the effective radius of Earth relative to 
!    WGS-84 reference ellipsoid
!
! INPUTS
!    lat    geographical latitude (in degree north).
!
! OUTPUT
!    R_lat  effective Earth radius (in m).
!
! NOTES
!
! SEE ALSO
!    gravity_msl
!    geopotential2geometric
!    geometric2geopotential
!
! REFERENCES
!    Mahoney, M.J., A discussion of various measures of altitudes,
!       URL: http://mtp.jpl.nasa.gov/notes/altitude/altitude.html,
!       cited 29/06/2004.
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

!-------------------------------------------------------------------------------
! 1. Double, scalar argument
!-------------------------------------------------------------------------------

function R_eff_double_0d(lat) result(R_eff)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  implicit none

  real(wp), intent(in) :: lat
  real(wp)             :: R_eff

  real(wp), parameter  :: a_earth = 6378137.0_wp      ! semi-major axis (m)
  real(wp), parameter  :: flatt = 0.003352811_wp      ! flattening
  real(wp), parameter  :: gm_ratio = 0.003449787_wp   ! gravitational ratio
  real(wp)             :: sinlat2
  real(wp), parameter  :: deg2rad = 0.0174532925_wp   ! degrees to radians 

! 1.2 Useful constant
! -------------------

  sinlat2 = sin(lat*deg2rad)**2

! 1.3 Calculate effective Earth radius
! ------------------------------------

  R_eff = a_earth / (1.0_wp + flatt + gm_ratio - 2.0_wp*flatt*sinlat2)

end function R_eff_double_0d


!-------------------------------------------------------------------------------
! 2. Double, array argument
!-------------------------------------------------------------------------------

function R_eff_double_1d(lat) result(R_eff)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  implicit none

  real(wp), dimension(:), intent(in) :: lat
  real(wp), dimension(size(lat))     :: R_eff

  real(wp), parameter  :: a_earth = 6378137.0_wp      ! semi-major axis (m)
  real(wp), parameter  :: flatt = 0.003352811_wp      ! flattening
  real(wp), parameter  :: gm_ratio = 0.003449787_wp   ! gravitational ratio
  real(wp), dimension(size(lat))     :: sinlat2
  real(wp), parameter  :: deg2rad = 0.0174532925_wp   ! degrees to radians 

! 2.2 Useful constant
! -------------------

  sinlat2 = sin(lat*deg2rad)**2

! 2.3 Calculate effective Earth radius
! ------------------------------------

  R_eff = a_earth / (1.0_wp + flatt + gm_ratio - 2.0_wp*flatt*sinlat2)

end function R_eff_double_1d

