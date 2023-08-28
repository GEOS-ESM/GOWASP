! $Id: ropp_fm_2d_plane.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****BendingAngle2d/ropp_fm_2d_plane
!
! NAME
!    ropp_fm_2d_plane  - Tool for calculating points in 2D plane for 2D
!                        Operator
!                            .
!
! SYNOPSIS
!    call ropp_fm_2d_plane(plat,plon,pazim,dtheta,n_horiz,plat_2d,plon_2d,kerror)
! 
! DESCRIPTION
!    This routine calculates the points in a 2D plane, given a lat,lon and
!    azimuthal angle (all in degrees). The points are separated by "dtheta" (radians). 
!    plat and plon are assumed to be the coordinates of the tangent point
!    and are placed at the centre of the plane.   
!
! INPUTS
!    plat,plon            ! location of tangent point (degs)
!    pazim                ! azimuthal angle
!    dtheta               ! angular separation between points in plane (rad)
!    n_horiz              ! number of points in plane
!    plat_2d,plon_2d      ! lat and lon of points in 2D plane
!    
! OUTPUT
!    plat_2d,plon_2d     :: vector of co-ords in 2d plane     ! Observation vector
!    kerror              :: error flag, 0 = ok
!
! NOTES
!
! SEE ALSO
!
! AUTHOR
!   ECMWF, Reading, UK.
!   Any comments on this software should be given via the GRAS SAF
!   Helpdesk at http://www.grassaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_fm_2d_plane(plat,plon,pazim,dtheta,n_horiz,plat_2d,plon_2d,kerror)

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_constants
  
  IMPLICIT NONE
  
  REAL(wp), INTENT(IN)       :: plat               ! Lat of tangent pt [deg]
  REAL(wp), INTENT(IN)       :: plon               ! Lon of tangent pt [deg]
  REAL(wp), INTENT(IN)       :: pazim              ! Angle [deg to north] 
  REAL(wp), INTENT(IN)       :: dtheta             ! angular sep in plane
                                                   ! is eastward
  INTEGER,  INTENT(IN)       :: n_horiz            ! Number of points         
  REAL(wp), INTENT(OUT)      :: plat_2d(n_horiz)   ! Lats of points [deg] 
  REAL(wp), INTENT(OUT)      :: plon_2d(n_horiz)   ! Lons of points [deg]
  INTEGER, INTENT(OUT)            :: kerror             ! Error code
  
  ! local
  
  REAL(wp)          :: plat_rad    ! in radians
  REAL(wp)          :: plon_rad    ! in radians
  REAL(wp)          :: pazim_rad   ! in radians
  REAL(wp)          :: ztheta(n_horiz)
  REAL(wp)          :: zcos_pazim
  REAL(wp)          :: zcos_plon
  REAL(wp)          :: zcos_plat
  REAL(wp)          :: zsin_plon
  REAL(wp)          :: zsin_plat
  REAL(wp)          :: zsin_theta
  REAL(wp)          :: zcos_theta
  REAL(wp)          :: zsin_plat_2d
  REAL(wp)          :: zcos_delta
  REAL(wp)          :: zrplat
  REAL(wp)          :: zazim
  INTEGER           :: i,imid

!-------------------------------------------------------------------------------
! 2. Convert to radians lat,lon and azimuthal angle to radians for calculation
!-------------------------------------------------------------------------------

  plat_rad  = plat  * pi/180.0_wp
  plon_rad  = plon  * pi/180.0_wp
  pazim_rad = pazim * pi/180.0_wp
  
  ! central profile
  
  imid = n_horiz/2 + 1  ! the central profile
  
!-------------------------------------------------------------------------------
! 3. Calculate the theta values in the plane
!-------------------------------------------------------------------------------

  DO i = 1, n_horiz
    ztheta(i) = REAL(i-imid)*dtheta  ! dtheta in radians
  ENDDO

  kerror = 0
  IF ( ABS(plat_rad) == 0.5_wp*pi ) THEN
    kerror = 1
    RETURN
  ENDIF
  
  zazim = pazim_rad
  IF (pazim_rad > pi ) zazim =  pazim_rad - 2.0_wp*pi  
  
  zcos_pazim = COS(zazim)
  zcos_plon  = COS(plon_rad)
  zcos_plat  = COS(plat_rad)
  zsin_plon  = SIN(plon_rad)
  zsin_plat  = SIN(plat_rad)
  
  DO i = 1, n_horiz
    
    zsin_theta = SIN(ztheta(i))  ! dtheta in radians
    zcos_theta = COS(ztheta(i))
    
    zsin_plat_2d = zsin_theta * zcos_pazim * zcos_plat + zsin_plat * zcos_theta
    
    zrplat = ASIN(zsin_plat_2d)
    
    IF( ABS(zrplat) == pi/2.0_wp ) THEN
      zcos_delta = 1.0_wp
    ELSE
      zcos_delta = ( zcos_theta - zsin_plat * zsin_plat_2d ) /   &
                      ( zcos_plat * COS(zrplat) )
      IF( zcos_delta >  1.0_wp)  zcos_delta =  1.0_wp
      IF( zcos_delta < -1.0_wp)  zcos_delta = -1.0_wp
    ENDIF
    
    IF( SIGN(1.0_wp,ztheta(i)) == SIGN(1.0_wp,zazim) ) THEN
      plon_2d(i) = plon_rad + ACOS(zcos_delta)
    ELSE
      plon_2d(i) = plon_rad - ACOS(zcos_delta)
    ENDIF
    
!-------------------------------------------------------------------------------
! 4. Locations in 2D plane, converted back into degrees
!-------------------------------------------------------------------------------
    
    plat_2d(i) = zrplat*180.0_wp/pi 
    
    plon_2d(i) = (MODULO(plon_2d(i)+pi,2.0_wp*pi) - pi)*180.0_wp/pi
    
  ENDDO
  
END SUBROUTINE ropp_fm_2d_plane
