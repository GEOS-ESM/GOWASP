! $Id: ropp_geodesy.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

module geodesy

!****m* Geodesy/geodesy *
!
! NAME
!    geodesy - Module interfacing various geodetic routines and functions.
!
! SYNOPSIS
!    use geodesy
! 
! DESCRIPTION
!    This module provides interfaces to some geodetic funtions and
!    subroutines of the ropp_utils library.
!
! SEE ALSO
!    R_eff
!    gravity
!    geometric2geopotential
!    geopotential2geometric
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
! 1. Effective radius of Earth according to Somigliana equation
!-------------------------------------------------------------------------------

  interface R_eff
     function R_eff_double_0d(lat) result(R_eff)
       use typesizes,  only: wp => EightByteReal
       real(wp), intent(in) :: lat
       real(wp)             :: R_eff
     end function R_eff_double_0d
     function R_eff_double_1d(lat) result(R_eff)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:), intent(in) :: lat
       real(wp), dimension(size(lat))     :: R_eff
     end function R_eff_double_1d
end interface

!-------------------------------------------------------------------------------
! 2. Gravity according to Somigliana equation
!-------------------------------------------------------------------------------

  interface gravity
     function gravity_double_00d(lat, h) result(g)
       use typesizes,  only: wp => EightByteReal
       real(wp), intent(in) :: lat
       real(wp), optional   :: h
       real(wp)             :: g
     end function gravity_double_00d
     function gravity_double_10d(lat, h) result(g)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:), intent(in) :: lat
       real(wp),               optional   :: h
       real(wp), dimension(size(lat))     :: g
     end function gravity_double_10d
     function gravity_double_01d(lat, h) result(g)
       use typesizes,  only: wp => EightByteReal
       real(wp),               intent(in) :: lat
       real(wp), dimension(:), intent(in) :: h
       real(wp), dimension(size(h))       :: g
     end function gravity_double_01d
  end interface

!-------------------------------------------------------------------------------
! 3. Altitude -> geopotential height conversion
!-------------------------------------------------------------------------------

  interface geometric2geopotential
     function geometric2geopotential_00d(lat, h) result(geop)
       use typesizes,  only: wp => EightByteReal
       real(wp), intent(in) :: lat
       real(wp), intent(in) :: h
       real(wp)             :: geop
     end function geometric2geopotential_00d
     function geometric2geopotential_01d(lat, h) result(geop)
       use typesizes,  only: wp => EightByteReal
       real(wp),               intent(in) :: lat
       real(wp), dimension(:), intent(in) :: h
       real(wp), dimension(size(h))       :: geop
     end function geometric2geopotential_01d
     function geometric2geopotential_10d(lat, h) result(geop)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:), intent(in) :: lat
       real(wp),               intent(in) :: h
       real(wp), dimension(size(lat))     :: geop
     end function geometric2geopotential_10d
     function geometric2geopotential_11d(lat, h) result(geop)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:), intent(in)      :: lat
       real(wp), dimension(:), intent(in)      :: h
       real(wp), dimension(size(lat), size(h)) :: geop
     end function geometric2geopotential_11d
  end interface

!-------------------------------------------------------------------------------
! 4. Geopotential height -> altitude conversion
!-------------------------------------------------------------------------------

  interface geopotential2geometric
     function geopotential2geometric_00d(lat, geop) result(h)
       use typesizes,  only: wp => EightByteReal
       real(wp), intent(in) :: lat
       real(wp), intent(in) :: geop
       real(wp)             :: h
     end function geopotential2geometric_00d
     function geopotential2geometric_01d(lat, geop) result(h)
       use typesizes,  only: wp => EightByteReal
       real(wp),                       intent(in) :: lat
       real(wp), dimension(:),         intent(in) :: geop
       real(wp), dimension(size(geop))            :: h
     end function geopotential2geometric_01d
     function geopotential2geometric_10d(lat, geop) result(h)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:), intent(in) :: lat
       real(wp),               intent(in) :: geop
       real(wp), dimension(size(lat))     :: h
     end function geopotential2geometric_10d
     function geopotential2geometric_11d(lat, geop) result(h)
       use typesizes,  only: wp => EightByteReal
       real(wp), dimension(:),         intent(in) :: lat
       real(wp), dimension(:),         intent(in) :: geop
       real(wp), dimension(size(lat), size(geop)) :: h
     end function geopotential2geometric_11d
  end interface

!-------------------------------------------------------------------------------
! 5. Great circle distance
!-------------------------------------------------------------------------------

  interface great_circle_distance
     function great_circle_distance_sca(lon1, lat1, lon2, lat2) result(distance)
       use typesizes, only: wp => EightByteReal
       real(wp), intent(in) :: lon1
       real(wp), intent(in) :: lat1
       real(wp), intent(in) :: lon2
       real(wp), intent(in) :: lat2
       real(wp)             :: distance
     end function great_circle_distance_sca
     function great_circle_distance_arr(lon1, lat1, lon2, lat2) result(distance)
       use typesizes, only: wp => EightByteReal
       real(wp), dimension(:), intent(in) :: lon1
       real(wp), dimension(:), intent(in) :: lat1
       real(wp), dimension(:), intent(in) :: lon2
       real(wp), dimension(:), intent(in) :: lat2
       real(wp), dimension(size(lon1))    :: distance
     end function great_circle_distance_arr
  end interface

end module geodesy
