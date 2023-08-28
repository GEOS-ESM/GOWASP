  program test_ft
!
! Test generation of random fields using random spectra for 
! desired correlation function shape in 1-D using Fourier transforms
!
  implicit none
!
  integer, parameter :: rkind=8
  integer, parameter :: ngrid=128
  integer, parameter :: nreal=40000
  integer, parameter :: nc=ngrid/2+1
  real(rkind) :: lscale=8.d0  ! length scale as inverse fraction of circle
!
  integer(8) :: iseed(1)
  integer :: nd(ngrid)
  integer ::ng,n,n2,nx,m
  real(rkind) :: pi2
  real(rkind) :: ds, deld, d, drds
  real(rkind) :: s1, xr, xi, ssum
  real(rkind) :: xg(ngrid,1), x(ngrid,nreal)
  complex(rkind) :: cg(nc,1), c(nc,nreal)
  real(rkind) :: s(nc), cc(nc)
  real(rkind) :: xm(ngrid), xv(ngrid), xc(ngrid,ngrid)
  real(rkind) :: xd(ngrid)
!
   pi2=8.d0*atan(1.d0)  ! = 2*pi
!
! Set random seed
   iseed(1)=1111
   call random_seed(put=iseed(1:1))
!
! Test Fourier transform (FT) on small problem 
   call ctest 
!
! Set desired correlation as function of separation distance
   ds=pi2/lscale          ! desired length scale
   deld=pi2/real(ngrid)   ! grid spacing
   do ng=1,ngrid/2+1 
     d=(ng-1)*deld - pi2/2.0d0
     drds=abs(d/ds)
     xg(ng,1)=exp(-0.5d0*(drds**2))                   ! Gaussian
!    xg(ng,1)=exp(-drds)                              ! Exponential
!    xg(ng,1)=(1.d0+drds+(drds**2)/3.0d0)*exp(-drds)  ! TOAR 
   enddo
   do ng=ngrid/2+2,ngrid  ! ensures symmetry about d=0
     xg(ng,1)=xg(ngrid-ng+2,1)
   enddo
!
   print *,'xg=',xg
!
! Compute FT coefs corresponding to xg
   call x2c (ngrid,nc,1,-1,xg,cg)
   print *,'cg=',cg
!
! Normalize coefficients so that the variance at any point among 
! all the realizations will be 1 
   ssum=0.d0
   do n=1,nc
     s(n)=sqrt(abs(real(cg(n,1))))
     if (n==1) then 
       ssum=ssum+s(n)**2
     else   ! account for presence of - wave number with conjugate coef.
       ssum=ssum+2.d0*s(n)**2
     endif
   enddo
   s(:)=s(:)/sqrt(ssum)
   print *,'ssum=',ssum
   print *,'s=',s
!
! Checking the average variance
   ssum=0.d0
   do n=1,nc
     if (n==1) then 
       ssum=ssum+s(n)**2
     else
       ssum=ssum+2.0d0*s(n)**2
     endif
   enddo
   print *,'ssum2=',ssum  ! should be 1
!
! Draw random values for coefficients for each realization, 
! assuming the spectral standard deviations s(n)
! Called routine uses mean=0, stdv=1. with max deviation=5 st. devs.
   do m=1,nreal
     do n=1,nc
       call random_gauss_r8 (0.d0,1.d0,5.d0,xr)   ! real part
       if (n == 1 .or. 2*(n-1) == ngrid) then     ! imag part =0
         c(n,m)=s(n)*cmplx(xr,0.d0)
       else
         call random_gauss_r8 (0.d0,1.d0,5.d0,xi) ! imaginary part
         c(n,m)=sqrt(0.5d0)*s(n)*cmplx(xr,xi)     ! same var for re and im parts
       endif
     enddo
   enddo
!
! Compute mean squared moduli of spectral coefs, assuming 0 mean.
   cc(:)=0.d0
   do n=1,nc
     do m=1,nreal
       cc(n)=cc(n)+c(n,m)*conjg(c(n,m))
     enddo
     cc(n)=cc(n)/real(nreal,rkind)
   enddo
!
   print *,'s**2, cc'
   do n=1,nc
     print ('(i4,2f14.8)'),n,s(n)**2,cc(n)
   enddo
!  
! Project realizations of random coefs into grid space 
   call x2c (ngrid,nc,nreal,1,x,c)
!
! Check mean of realizations at each point (should be 0 everywhere)
   xm(:)=0.d0
   do m=1,nreal
     do n=1,ngrid
       xm(n)=xm(n)+x(n,m)
     enddo
   enddo
   xm(:)=xm(:)/real(nreal)
   print *,'xm=',xm
!
! Check variance of realizations at each point (should be 1 everywhere)
   xv(:)=0.d0
   do m=1,nreal
     do n=1,ngrid
       xv(n)=xv(n)+(x(n,m)-xm(n))**2
     enddo
   enddo
   xv(:)=xv(:)/real(nreal-1)
   print *,'xv=',xv
!
! Compute correlations between all pairs of points   
   xc(:,:)=0.d0    
   do n=1,ngrid
     do n2=1,ngrid
       do m=1,nreal
         xc(n,n2)=xc(n,n2)+(x(n,m)-xm(n))*(x(n2,m)-xm(n2))
       enddo
       xc(n,n2)=xc(n,n2)/real(nreal-1)
       xc(n,n2)=xc(n,n2)/sqrt(xv(n)*xv(n2))
     enddo
   enddo  
!
   print *,'xc='
   do n=1,16
     print ('(a,i2,16f6.2)'),'n=',n,(xc(n,n2),n2=1,16)
   enddo  
!
! Compute average correlations for each separation distance 
! Need only consider 0 <= |d| <= pi, since field is periodic.
   xd(:)=0.d0
   nd(:)=0
   do n=1,ngrid
     do n2=n,ngrid
       nx=n2-n+1
       if (nx>ngrid/2+1) then
         nx=n+(ngrid-n2)+1
       endif
       xd(nx)=xd(nx)+xc(n,n2)
       nd(nx)=nd(nx)+1
     enddo
   enddo
!
   do nx=1,ngrid/2+1
     xd(nx)=xd(nx)/nd(nx)
   enddo
   print *,'nd=',nd(1:ngrid/2)
   print *,'xd=',xd(1:ngrid/2)
!
! Print desired and computed correlations as function of distance
   print *,' ' 
   print *,'xg,d'
   nx=ngrid/2
   do n=1,nx
     print ('(a,i3,2f10.4)'),'n ',n,xg(n+nx,1),xd(n)
   enddo
!  
   end program test_ft
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
     subroutine random_gauss_r8 (mean,stdv,xmax,draw1)
!
! Retrieve a random number for a Gaussian distribution with mean and  
! standard deviation specified as input.  Any random draw more than +/- 
! xmax standard deviations is set to +/- xmax standard deviations.  
! This algorithm uses the central limit theorem by averaging over ndraw 
! random variables drwan from a uniform distribution between 0 and 1.
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
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
!
     subroutine x2c (ng,nc,mm,ix,x,c)
!
!  Fourier transform or its inverse
!
     implicit none
!
     integer, parameter :: rkind=8
     integer :: ng,nc,mm,ix
     integer :: k,n,m
     real(rkind) :: x(ng,mm)
     complex(rkind) :: c(nc,mm)
     real(rkind) :: pi,pi2,ag,xc,xs,cc,fac,rng,xf,km1
!
     pi=4.0d0*atan(1.0d0)
     pi2=2.0d0*pi
     rng=real(ng,rkind)
     fac=pi2/rng
!
     if (ix == -1) then ! x to c
       c(:,:)=cmplx(0.d0,0.d0,rkind)
       do k=1,nc
         km1=real(k-1,rkind)
!
! If ng is even and m=k-1 is largest resolved wavenumber, then divide 
! derived coef by 2 since for this grid and m, the wavenumbers m and -m 
! are degenerate (i.e., both the same). this factor defines the coef as
! if this were not the case, so that all m>0 are treated the same; i.e.,
! as if the present wavenumbers were -M, -(M-1), ..., -1, 0, 1, ..., M-1, M 
         if (2*(k-1)==ng) then
           rng=2.*real(ng,rkind)
         else
           rng=real(ng,rkind)
         endif
         do n=1,ng
           ag=km1*(real(n-1,rkind)*fac-pi)
           xc=cos(ag)/rng
           xs=sin(ag)/rng
           do m=1,mm
             c(k,m)=c(k,m)+x(n,m)*cmplx(xc,-xs,rkind)
           enddo
         enddo
       enddo
!
     else              ! c to x
       x(:,:)=0.d0
       do k=1,nc
         km1=real(k-1,rkind)
!
! The following assumes wavenumber -M and M must be considered (see above)
         if (k==1) then
           xf=1.d0
         else
           xf=2.d0
         endif
         do n=1,ng
           ag=km1*(real(n-1,rkind)*fac-pi)
           xc=cos(ag)
           xs=sin(ag)
           do m=1,mm
             cc=real(c(k,m)*cmplx(xc,xs,rkind))
             x(n,m)=x(n,m)+cc*xf
           enddo
         enddo
       enddo
!
     endif
!     
     end subroutine x2c
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
     subroutine ctest 
!
     implicit none
!
     integer, parameter :: rkind=8
     integer, parameter :: nn=12
     integer, parameter :: nc=1+nn/2
     integer, parameter :: mm=2
     real(rkind) :: x(nn,mm)
     complex(rkind) :: c(nc,mm)
!
     c(:,:)=cmplx(0.d0,0.d0)
     c(2,1)=cmplx(1.d0,0.d0)
     c(3,2)=cmplx(1.d0,2.d0)
     c(1,2)=cmplx(3.d0,0.d0)
     c(nc,2)=cmplx(5.d0,0.d0)
     x(:,:)=0.d0
!
     call x2c (nn,nc,mm,1,x,c)
     print *,' '    
     print *,'Test x1:'
     print *,x(:,1)
     print *,'Test x2:'
     print *,x(:,2)
!
     c(:,:)=cmplx(0.d0,0.d0)
     call x2c (nn,nc,mm,-1,x,c)
     print *,' '    
     print *,'Test c1:'
     print *,c(:,1)
     print *,'Test c2:'
     print *,c(:,2)
     print *,' ' 
!
     end subroutine ctest
