! $Id: ropp_gravity.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****f* Geodesy/gravity *
!
! NAME
!    gravity - Calculate gravity following the Somigliana equation
!
! SYNOPSIS
!    g_lat = gravity( lat, [ h ] )
! 
! DESCRIPTION
!    This function returns gravity as function of geographic latitude and 
!    (optionally) altitude above the WGS-84 reference ellipsoid as derived 
!    from the Somigliana equation 
!    [see http://mtp.jpl.nasa.gov/notes/altitude/altitude.html for details]
!
! INPUTS
!    lat     geographical latitude (in degree north).
!    h       altitude above reference ellipsod (in m, optional). If z is not
!               given, the function returns the gravity at the ellipsoid.
!
! OUTPUT
!    g       gravity (in m/s^2).
!
! NOTES
!
! SEE ALSO
!    R_eff
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

function gravity_double_00d(lat, h) result(g)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use geodesy,   only: R_eff

  implicit none

  real(wp), intent(in) :: lat
  real(wp), optional   :: h
  real(wp)             :: g
  real(wp)             :: R_lat
  real(wp), parameter  :: g_equat = 9.7803253359_wp  ! equatorial gravity
  real(wp), parameter  :: k_somig = 1.931853e-3_wp   ! Somagliana's constant
  real(wp), parameter  :: ecc = 0.081819_wp          ! eccentricity
  real(wp)             :: sinlat2
  real(wp), parameter  :: deg2rad = 0.0174532925_wp  ! degrees to radians

! 1.2 Useful constant
! --------------------
               
  sinlat2  = sin(lat*deg2rad)**2
      
! 1.2 Calculate gravity - normal gravity on surface of ellipsoid
! ---------------------

  g = g_equat*(1.0_wp + k_somig*sinlat2) / SQRT(1.0_wp - (ecc**2)*sinlat2)

! 1.3 Altitude dependence
! -----------------------

  if (present(h)) then
     R_lat = R_eff(lat)
     g = g * R_lat**2 / (R_lat + h)**2
  endif

end function gravity_double_00d


!-------------------------------------------------------------------------------
! 2. Double, array argument for latitude
!-------------------------------------------------------------------------------

function gravity_double_10d(lat, h) result(g)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use geodesy,   only: R_eff

  implicit none

  real(wp), dimension(:), intent(in) :: lat
  real(wp),               optional   :: h
  real(wp), dimension(size(lat))     :: g

  real(wp), dimension(size(lat))     :: R_lat
  real(wp), dimension(size(lat))     :: sinlat2

  real(wp), parameter  :: g_equat = 9.7803253359_wp  ! equatorial gravity
  real(wp), parameter  :: k_somig = 1.931853e-3_wp   ! Somagliana's constant
  real(wp), parameter  :: ecc = 0.081819_wp          ! eccentricity
  real(wp), parameter  :: deg2rad = 0.0174532925_wp  ! degrees to radians 

! 2.2 Useful constant
! -------------------

  sinlat2  = sin(lat*deg2rad)**2

! 2.3 Calculate gravity
! ---------------------

  g = g_equat*(1.0_wp + k_somig*sinlat2) / SQRT(1.0_wp - (ecc**2)*sinlat2)

! 2.4 Altitude dependence
! -----------------------

  if (present(h)) then
     R_lat = R_eff(lat)
     g = g * R_lat**2 / (R_lat + h)**2
  endif

end function gravity_double_10d


!-------------------------------------------------------------------------------
! 3. Double, array argument for altitude
!-------------------------------------------------------------------------------

function gravity_double_01d(lat, h) result(g)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use geodesy,   only: R_eff

  implicit none

  real(wp),               intent(in) :: lat
  real(wp), dimension(:), intent(in) :: h
  real(wp), dimension(size(h))       :: g

  real(wp)             :: R_lat
  real(wp), parameter  :: g_equat = 9.7803253359_wp  ! equatorial gravity
  real(wp), parameter  :: k_somig = 1.931853e-3_wp   ! Somagliana's constant
  real(wp), parameter  :: ecc = 0.081819_wp          ! eccentricity
  real(wp)             :: sinlat2
  real(wp), parameter  :: deg2rad = 0.0174532925_wp  ! degrees to radians 


! 2.2 Useful constant
! -------------------

  sinlat2  = sin(lat*deg2rad)**2

! 2.3 Calculate gravity
! ---------------------

  g = g_equat*(1.0_wp + k_somig*sinlat2) / SQRT(1.0_wp - (ecc**2)*sinlat2)

! 2.4 Altitude dependence
! -----------------------

  R_lat = R_eff(lat)
  g = g * R_lat**2 / (R_lat + h)**2

end function gravity_double_01d
