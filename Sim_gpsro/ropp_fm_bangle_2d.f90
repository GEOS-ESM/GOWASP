! $Id: ropp_fm_bangle_2d.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****s* BendingAngle2d/ropp_fm_bangle_2d *
!
! NAME
!    ropp_fm_bangle_2d - Forward model calculating a bending
!                        angle profile from planar information.
!
! SYNOPSIS
!    call ropp_fm_bangle_2d(2D_state, obs)
! 
! DESCRIPTION
!    This routine is a forward model calculating a vertical profile of 
!    bending angles from temperature, humidity and surface pressure. Bending
!    angle values are calculated for the impact parameters given in the
!    observation vector.
!
! INPUTS
!    type(State2dFM)       :: x     ! State vector
!    type(Obs1dBangle)     :: y     ! Observation vector (levels required)
!
! OUTPUT
!    type(Obs1dBangle)     :: y     ! Observation vector
!
! NOTES
!    The 2D forward model evaluates bending angles on the observed impact
!    parameter values using a ray tracer to solve the ray equations. The
!    operator neglects tangent point drift and uses the observed impact
!    parameter to determine the ray tangent height. 
! 
!
! SEE ALSO
!    ropp_fm_types
!    ropp_fm_bangle_2d_ad
!    ropp_fm_bangle_2d_tl
!
! AUTHOR
!   ECMWF, UK.
!   Any comments on this software should be given via the GRAS SAF
!   Helpdesk at http://www.grassaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_fm_bangle_2d(azim_flag, x, y)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
!qqq  USE ropp_fm,   not_this => ropp_fm_bangle_2d
  USE geodesy,   ONLY: gravity, R_eff
  USE ropp_fm_types
  USE ropp_fm_constants

  IMPLICIT NONE

  INTEGER, INTENT(in) :: azim_flag  
  TYPE(state2dFM),   INTENT(inout) :: x
  TYPE(Obs1dBangle), INTENT(inout) :: y

  INTEGER                                  :: i
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: pwvp
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: pdry
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: z
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: rad
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: impact
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)   :: refrac
  REAL(wp), DIMENSION(y%nobs,2)            :: a_path
  REAL(wp)                                 :: glat,g_ratio,rad_eff   

! non-ideal gas 

  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)     :: z_geop          ! geopotential height of model levels
  REAL(wp), DIMENSION(x%n_lev,x%n_horiz)     :: zcomp_dry_inv, zcomp_wet_inv ! compressibilities
  
  REAL(wp)                         :: kap1,kap2,kap3  ! refractivity coefficients used in routine
      
  ! Set maximum height for 2D calulation, perform 1D calculation above
  REAL(wp), PARAMETER                      :: z_xd(2) = (/1.0e0_wp, 1.0e4_wp/) 
  REAL(wp)                                 :: z_2d 
!QQQ  REAL(wp), PARAMETER                      :: z_2d = 1.0e4_wp  

  ! Set step size in ray-tracer
  INTEGER,  PARAMETER                      :: msplit = 4 



!-------------------------------------------------------------------------------
! 1. Non-ideal gas options
!-------------------------------------------------------------------------------

! set inverse of compressibilities

  zcomp_dry_inv(:,:) = 1.0_wp
  zcomp_wet_inv(:,:) = 1.0_wp
   
! initialise geopotential heights
   
  z_geop(:,:) = x%geop(:,:)
 
  IF (x%non_ideal) THEN

! if non ideal gas calculation, use adjusted coefficients

     kap1 = kappa1_comp
     kap2 = kappa2_comp
     kap3 = kappa3_comp   

!    calculate compressibilty and adjust geopotential heights in z_geop
     CALL ropp_fm_compress_2d(x,z_geop,zcomp_dry_inv,zcomp_wet_inv)  

  ELSE

     kap1 = kappa1
     kap2 = kappa2
     kap3 = kappa3
     
  ENDIF  


! 1.1 Calculate water vapor  and dry pressure
! ----------------------------------

  pwvp(:,:) = x%pres(:,:) * x%shum(:,:) / (epsilon_water + (1.0_wp - epsilon_water)*x%shum(:,:))

  pdry(:,:) = x%pres(:,:) - pwvp(:,:)

 ! 1.2 Calculate refractivity
! --------------------------

  refrac(:,:) = kap1 * pdry * zcomp_dry_inv/ x%temp +  &
           kap2 * pwvp * zcomp_wet_inv/ x%temp**2 + &
           kap3 * pwvp * zcomp_wet_inv/ x%temp
 
! 1.3 Calculate geometric height at each location in plane
! --------------------------------------------------------
  DO i = 1, x%n_horiz
  
! g value and effective radius used in the geopotential -> geometric transform
    
    glat=gravity(x%lat(i))   
    g_ratio = glat/g_wmo  
    rad_eff = R_eff(x%lat(i))

! geometric height
  
    z(:,i) =  z_geop(:,i)/ (g_ratio - z_geop(:,i)/rad_eff)

  ENDDO 
  
! 1.4 Calculate radius
! ------------------------------

  rad(:,:) = z(:,:) + y%r_curve + y%undulation
  
! 1.5 Calculate radius-refractive index product 
! ----------------------------------------------
 
  impact(:,:) = (1.0_wp + 1.e-6_wp*refrac(:,:)) * rad(:,:) 

! store the refractivity and nr product in x. Useful for diagnostics
! These variables are not active, so they don't appear in TL and AD routines.    

  x%refrac(:,:) = refrac(:,:)
  x%nr(:,:)     = impact(:,:)

! 1.6 Calculate bending angles with a 2D operator
! -----------------------------------------------

! call 2D operator based on a 4th order Runge-Kutta solution of the ray equations. 

  z_2d=z_xd(azim_flag) ! 2d calc will be done below 1m (i.e., not done) or below 10km
  CALL ropp_fm_alpha2drk(y%nobs, x%n_lev, x%n_horiz, msplit, x%dtheta, y%impact, &
                         refrac, rad, impact, y%r_curve, z_2d, a_path, y%bangle, y%rtan)

! 1.7 Store the a=nr*sin*(phi) value at ray endpoints 
! ---------------------------------------------------
  y%a_path(:,:) = a_path(:,:)

END SUBROUTINE ropp_fm_bangle_2d
