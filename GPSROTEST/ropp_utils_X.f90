! $Id: ropp_utils_X.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

module ropp_utils

!****m* Modules/ropp_utils *
!
! NAME
!    ropp_utils - The ROPP utils library
!
! SYNOPSIS
!    use ropp_utils
!
! DESCRIPTION
!    This Fortran module provides interfaces required for the use of the
!    ROPP utils library
!
! SEE ALSO
!    arrays
!    geodesy
!    coordinates
!    datetimeprogs
!    messages
!    unitconvert
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
! 1. Include all ropp_utils modules
!-------------------------------------------------------------------------------

  use typesizes, only: wp => EightByteReal

!qqqq  use arrays
!qqqq   use geodesy
!qqqq   use coordinates
!qqqq   use datetimeprogs
!qqqq   use messages
!qqqq   use earthmod
!qqqq   use unitconvert



!-------------------------------------------------------------------------------
! 2. ROPP missing value definitions
!-------------------------------------------------------------------------------

!****ip* Initialisation/ropp_mdfv
!
! NAME
!    ropp_MDFV - Internal global 'Missing Data Flag/Test Value(s)'
!
! NOTES
!    These parameters are used to indicate/test a 'missing' (invalid)
!    data value.
!    ropp_MDFV should be used to set invalid data for most
!      ROPP parameters, but a single universal value is not suitable for
!      all; some - e.g. (X,Y,Z) coordinate vectors - are set to zero.
!      For others, parameter-specific values may be more appropriate.
!    ropp_ZERO can be used to set parameters to zero.
!    ropp_MDTV can used for testing for invalid parameter values;
!      anything less than this value can be assumed to be set 'missing',
!      though again, some parameters may have specific values to test for.
!    ropp_ZDTV can be used to test for (almost) zero, e.g.
!      if ( abs(value) < ropp_ZDTV ) then ...
!         ! value can be considered to be zero
!
! SOURCE
!
  REAL(wp), PARAMETER              :: ropp_MDFV = -99999000.0_wp
  REAL(wp), PARAMETER              :: ropp_ZERO =     0.0_wp

  REAL(wp), PARAMETER              :: ropp_MDTV = -9999.0_wp
  REAL(wp), PARAMETER              :: ropp_ZDTV =   1e-10_wp
!
!****
!qqqq various interfaces removed

end module ropp_utils
