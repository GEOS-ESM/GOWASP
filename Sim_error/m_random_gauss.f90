     module m_random_gauss  
!
! Module for drawing random numbers from a Gaussian distribution 
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
     public random_gauss_r8
!
     contains
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
     subroutine random_gauss_r8 (mean,stdv,xmax,draw1)
!
! Retrieve an 8-byte random number for a Gaussian distribution with 
! mean and standard deviation specified as input. 
! Any random draw more than +/- xmax standard deviations is set to +/- xmax 
! standard deviations.  This algorithm uses the central limit theorem by 
! averaging over ndraw random variables drwan from a uniform distribution 
! between 0 and 1.  Although this is not the best algorithm to use a Gaussian
! distribution, it is sufficient for the purpose here.
!
     implicit none
     integer, parameter :: r8=8      ! precision of real variables
     integer, parameter :: ndraw=12  ! # of rand # from uniform dist to avg
!
     real(r8), intent(in)   :: mean  ! specified mean of dist.
     real(r8), intent(in)   :: stdv  ! specified standard deviation of dist. 
     real(r8), intent(in)   :: xmax  ! don't allow more than this many stdevs
     real(r8), intent(out)  :: draw1 ! random draw from Gaussian distribution
!
     integer      :: n
     real(r8) :: x, xsum
!
     xsum=0._r8
     do n=1,ndraw
       call random_number (x)
       xsum=xsum+x
     enddo
! The following x is normal(0,1).  The factor 12 appearing below is
! the inverse of the variance of a uniform distribution spanning the range 
! -0.5<=x<=0.5.      
     x=(xsum-0.5_r8*ndraw)*sqrt(12._r8/ndraw)
     x=min( xmax,x)     
     x=max(-xmax,x)
     draw1=x*stdv+mean
!
     end subroutine random_gauss_r8
!      
     end module m_random_gauss  
