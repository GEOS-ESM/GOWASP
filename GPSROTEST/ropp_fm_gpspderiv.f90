! $Id: ropp_fm_gpspderiv.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****s* BendingAngle2d/ropp_fm_gpspderivs *
!
! NAME
!    ropp_fm_gpspderivs - calculates the refractivity derivatives used in
!                 the ray-tracer 
!          
!
! SYNOPSIS
!    call ropp_fm_gpspderivs(klev,khoriz, ...)
! 
! DESCRIPTION
!    This routine calculates the refractivity gradients used in the
!    raytracer.  
!
! INPUTS
!
!           klev   =  number of vertical levels
!           khoriz =  number of horizontal locations
!           kk      = in the kkth box in the horizontal
!           pdsep  =  angular spacing
!	    ptheta_min = minium theta 
!           ptheta_max = maximum theta 
!	    ptheta_tan = theta at tangent point & 
!	    prtan = radius of tangent point
!           pamult = 1 or -1 depending on direction
!           prefrac = refractivity
!           pradius = radius value
!           py = current location
!
! OUTPUT
!
!           pdydh = derivs used in the RK calculation
!
! NOTES
!    
!     We limit the magnitude of the radial refractivity gradient to 1.5e-7, to stop the ray
!     getting stuck 
!
!     More information on method in Healy et al, (2007), Q.J.R. Meteorol.Soc., vol 133, 1213-1217, 
!     section 3.2  
!    
!
! SEE ALSO
!    ropp_fm_types
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

SUBROUTINE ropp_fm_gpspderivs(klev,   & ! no.of observations
		      khoriz, & ! no. of horizontal layers  ODD
		      kk, &
		      pdsep, &
		      ptheta_min, & 
		      ptheta_max, & 
		      ptheta_tan, & 
		      prtan, &
		      pamult, &
		      prefrac, &
		      pradius, &
		      py, &
		      pdydh)

USE typesizes, only: wp => EightByteReal
		   
IMPLICIT NONE
		   
INTEGER, INTENT(IN)  :: klev           ! no. of refractivity levels
INTEGER, INTENT(IN)  :: khoriz         ! no. of horizontal locations
INTEGER, INTENT(IN)  :: kk 
REAL(KIND=wp),    INTENT(IN)  :: pdsep           ! angular spacing of grid
REAL(KIND=wp),    INTENT(IN)  :: ptheta_min
REAL(KIND=wp),    INTENT(IN)  :: ptheta_max
REAL(KIND=wp),    INTENT(IN)  :: ptheta_tan
REAL(KIND=wp),    INTENT(IN)  :: prtan
REAL(KIND=wp),    INTENT(IN)  :: pamult          !
REAL(KIND=wp),    INTENT(IN)  :: prefrac(klev,khoriz)   ! refractivity values on levels
REAL(KIND=wp),    INTENT(IN)  :: pradius(klev,khoriz)   ! radius values
REAL(KIND=wp),    INTENT(IN)  :: py(4)   !  current location
REAL(KIND=wp),    INTENT(OUT) :: pdydh(4)

! local

INTEGER :: ik,ikp1
REAL(KIND=wp) :: zhwt1,zhwt2
REAL(KIND=wp) :: zref_up,zref_low
REAL(KIND=wp) :: zrad_up,zrad_low
REAL(KIND=wp) :: zkval,zdndr,zdndr2
REAL(KIND=wp) :: zed



! the easy bits


pdydh(1) = MAX(1.0E-10_wp,COS(py(3)))

pdydh(2) = pamult*SIN(py(3))/(py(1)+prtan)


IF ( py(2) >= ptheta_min .AND. py(2) <= ptheta_max) THEN

   ik = INT((py(2) + ptheta_tan)/pdsep)+1
   ik = MIN(MAX(1,ik),khoriz)
   ikp1 = MIN(khoriz,ik+1)	

! horizontal weighting factor	       
	       
   zhwt1 = (REAL(ik)*pdsep - (py(2)+ptheta_tan))/pdsep   
   zhwt2 = 1.0_wp - zhwt1

ELSE IF (py(2) < ptheta_min ) THEN

   ik = 1
   ikp1 = 2
   zhwt1 = 1.0_wp
   zhwt2 = 0.0_wp
   
ELSE IF  (py(2) > ptheta_max) THEN

   ik = khoriz -1 
   ikp1 = khoriz
   zhwt1 = 0.0_wp
   zhwt2 = 1.0_wp
   
ENDIF    

! now calculate the local radial gradient of n
      
zref_up  = zhwt1*prefrac(kk+1,ik)+zhwt2*prefrac(kk+1,ikp1)
zref_low = zhwt1*prefrac(kk,ik) + zhwt2*prefrac(kk,ikp1)

zrad_up = zhwt1*pradius(kk+1,ik)+zhwt2*pradius(kk+1,ikp1)
zrad_low = zhwt1*pradius(kk,ik) + zhwt2*pradius(kk,ikp1)


zkval = LOG(zref_low/zref_up)/(zrad_up-zrad_low)

zed = MAX(0.0_wp,(py(1)+prtan-zrad_low))

zdndr = - 1.0E-6_wp*zkval*zref_low*EXP(-zkval*zed)

zdndr2 = MAX(-1.5E-7_wp,zdndr)


!!if (abs(zdndr) > 0.15E-6_wp) then

!!   write (*,*) '################################'
   
!   write (*,*) 'biggradients',kk,ik,zhwt1,zhwt2,zkval
!   write (*,*) 'big-gradients',py(1),prtan,zrad_low,(py(1)+prtan-zrad_low)
!   write (*,*) 'big grads',pradius(kk,ik),pradius(kk,ikp1)
!   write (*,*) 'big rads2',pradius(kk+1,ik),pradius(kk+1,ikp1)
!   write (*,*) 'big ref1',prefrac(kk,ik),prefrac(kk,ikp1)
!   write (*,*) 'big ref2',prefrac(kk+1,ik),prefrac(kk+1,ikp1)
!   write (*,*) 'rads',zrad_up,zrad_low
!   write (*,*) 'refs',zref_up,zref_low
!   write (*,*) 'ref-gradient',1.0E6_wp*zdndr
!   write (*,*) 'deriv3',-SIN(py(3))*(1.0_wp/(py(1)+prtan) + zdndr)
!   write (*,*) 'deriv4',- SIN(py(3))*zdndr   
   

!   write (*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 

!endif 

pdydh(3) = -SIN(py(3))*(1.0_wp/(py(1)+prtan) + zdndr2)

pdydh(4) = - SIN(py(3))*zdndr



RETURN

END SUBROUTINE ropp_fm_gpspderivs
