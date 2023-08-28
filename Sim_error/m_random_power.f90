   module m_random_power
!
! Determine expected spherical-harmonic power spectra for desired 
! homogeneous, isotropic correlation function on the sphere. 
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1
   use m_parameters, only : earthr
!
   implicit none
!
   private
   public :: random_power_compute
!
   private :: rp_alp0,rp_emns0,rp_bssl,rp_gauss
!
   integer, parameter :: rkindc=8  ! precision of spectral code used here
   integer, parameter :: r4=4      ! precision of real vars on file_err_corr 
!
   real(rkindc), parameter :: zero=0._rkindc
   real(rkindc), parameter :: one=1._rkindc
   real(rkindc), parameter :: two=2._rkindc
   real(rkindc), parameter :: three=3._rkindc
!
   contains
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine random_power_compute (nmax,nlevels,lev1,lev2,lprint, &
                                    corr_shape,hcorr_lengths,power,ierr)
!
!  Determine expected spherical-harmonic power spectra for desired 
!  homogeneous, isotropic correlation function on the sphere.
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(in) :: nmax
   integer, intent(in) :: nlevels,lev1,lev2
   integer, intent(out) :: ierr
   real(r4), intent(in)  :: hcorr_lengths(nlevels)
   real(rkind1), intent(out) :: power(0:nmax,nlevels)
   character(len=*), intent(in) :: corr_shape
!
   integer, parameter :: iordlt=0
   integer :: nlats
   integer :: j,n,lev 
!
   real(rkindc) :: piover2 ! pi/2 
   real(rkindc) :: d,ds,corr,psum,xc,earthrkm
   real(rkindc), allocatable :: sinp(:)
   real(rkindc), allocatable :: theta(:)
   real(rkindc), allocatable :: weights(:)
   real(rkindc), allocatable :: alp(:)
   real(rkindc), allocatable :: epsi(:)
   real(rkindc), allocatable :: trans_mat(:,:)
!
   if (trim(corr_shape) == 'WHITE' .or. &
       trim(corr_shape) == 'GAUSS' .or. &
       trim(corr_shape) == 'EXP'   .or. &
       trim(corr_shape) == 'TOAR' ) then 
     ierr=0
   else
     ierr=-1
     if (lprint) then
       print *,' '
       print *,'ERROR: Specified Correlation function not acceptable'
       print *,'  User specified in file = ',trim(corr_shape)
       print *,'  Acceptible= WHITE, GAUSS, EXP, or TOAR'
       return
     endif
   endif
!
   earthrkm=0.001*earthr
   piover2=two*atan(one)
   nlats=2*(nmax/2)+4
!
   allocate (sinp(nlats),theta(nlats),weights(nlats))
   allocate (alp(0:nmax),epsi(0:nmax+1),trans_mat(0:nmax,nlats))
!
! Compute weights and lats for Legendre transform
! (use a set of Gaussian latitudes for this purpose)
   call rp_gauss (sinp,weights,nlats,iordlt)
   do j=1,nlats
     theta(j)=asin(sinp(j))  ! determine "grid" lat given its sine
   enddo
!
   call rp_emns0 (nmax,epsi)
   do j=1,nlats
     call rp_alp0 (nmax,sinp(j),epsi,alp) ! determine P_(m,n) for m=0 
     do n=0,nmax
       trans_mat(n,j)=alp(n)*weights(j)    ! matrix used for inverse transform
     enddo     
   enddo
!
! Determine required power spectrum from values of legendre coefficients for 
! P(m=0,n) by using inverse transform applied to the desired correlation 
! function shape.  Values of power here are non-normalized. 
   if (trim(corr_shape) /= 'WHITE') then
     power(:,:)=0. 
     do lev=lev1,lev2
       do j=1,nlats  
!
         d=earthrkm*(piover2-theta(j)) ! distance to N pole
         if (hcorr_lengths(lev) > 1.0) then
           ds=d/hcorr_lengths(lev)   
           if (trim(corr_shape) == 'GAUSS') then 
             corr=exp(-ds*ds/two)
           elseif (trim(corr_shape) == 'EXP') then 
             corr=exp(-ds) 
           elseif (trim(corr_shape) == 'TOAR') then 
             corr=(one+ds+ds*ds/three)*exp(-ds) 
           endif
!
         else ! small corr length so set corr as delta funct next to N. pole
           if (j==1 .and. theta(1) > 0.)  then            ! j=1 is nearest N.P.
             corr=one
           elseif (j==nlats .and. theta(nlats) > 0.) then ! j=nlats nearest N.P.
             corr=one
           else
             corr=zero
           endif
         endif  
 !               
         do n=0,nmax
           power(n,lev)=power(n,lev)+trans_mat(n,j)*corr
         enddo
!   
       enddo  ! loop over j
     enddo    ! loop over lev
!
! The following scaling by the factor xc appears to be necessary for the 
! desired shapes and the shapes computed from the resulting random fields 
! to be the same. This is unlike what is described in the literature cited at
! the beginning of this module. I suspect this concerns the normalization of
! the Legendre polynomials: whether they are normalized such that the integral 
! of their squares is a constant (e.g., 1 or 2) over the domain or 1/(2n+1).  
! Here, the polynomials are normalized such that it is 2. (It does not 
! matter what this constant is since it is normalized out in this context to 
! yield a global mean variance of 1 for the random field, just so long as it
! is m and n independent.)  The addition theorem of spherical harmonics also
! involves a factor of 2n+1, depending on the normalization used. It is 
! the consideration of these 2 aspects that yields the required net factor 
! here of 1/sqrt(2n+1).
! 
     do n=0,nmax
       xc=sqrt(real(2*n+1,rkindc))
       do lev=lev1,lev2
         power(n,lev)=abs(power(n,lev))/xc
       enddo
     enddo
!
   else       ! desired spectrum is WHITE
     power(:,:)=1.  ! will result in all normalized values being the same
   endif      ! check if not WHITE 
!
! Normalize expected variances so that global-integrated variance is 1. 
! xc=2*n+1 is the number of coefficients m,n for each n, including m<0
   do lev=lev1,lev2
     psum=zero
     do n=0,nmax
       xc=real(2*n+1,rkindc)
       psum=psum+power(n,lev)*xc
     enddo    
     power(:,lev)=power(:,lev)/psum
   enddo
!
   deallocate (sinp,theta,weights)
   deallocate (alp,epsi,trans_mat)
!
   if (lprint) then
     print ('(2a)'),' Expected power spectrum specified for corr function ', &
                    trim(corr_shape)
   endif
!
   end subroutine random_power_compute
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rp_alp0 (nlast,sinp,epsi,alp)
!
!  Compute associated Legendre polynomials for m=0 and latitude.
!  The index of alp is the value n-m. The polynomials are normalized such
!  that the integral of their squared value from pole to pole equals 1. 
!  This renormalization is consistent with a weighting that represents a 
!  fraction of the surface area on the sphere. This normalization is 
!  specified by the value of abfac.
!
!  If this is for a vector field at the poles, disregard m=0 and 
!  instead consider alp(n) = lim (P_m^(n)/cos(lat)) as lat => pole
!  (This yields non zero values only for m=1)
!
   implicit none
!
   integer, intent(in) :: nlast      
   real(rkindc), intent(in) :: epsi(0:nlast+1)
   real(rkindc), intent(in) :: sinp
   real(rkindc), intent(out) :: alp(0:nlast)
!
   integer :: m
   integer :: nm  ! = n-m
   integer :: np  ! indicator for order of alp storage 
   real(rkindc) :: fm  ! = real(m)  
!
   m=0
   fm=real(m,rkindc)
   alp(0)=one   ! specifically for for m=0
!
! Compute values for n-m=1
   if (nlast > 0) then
     alp(1)=sqrt(two*fm+three)*sinp*alp(0)
   endif
!
! Compute values for n-m>1 using recursion
   if (nlast > 1) then
     do nm=2,nlast
       alp(nm)=(sinp*alp(nm-1)-epsi(nm-1)*alp(nm-2))/epsi(nm)
     enddo
   endif
!
   end subroutine rp_alp0
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rp_emns0 (nmax,epsi)   
!
!  Compute factors required for computing associated Legendre polynomials
!
   implicit none
!
   integer :: nmax
   real(rkindc) :: epsi(0:nmax+1)
   integer :: m,n,n1,nm
   real(rkindc) :: fnum,fden
!
   do m=0,0
     if (m == 0) then
       epsi(0)=zero    
       n1=1
     else
       n1=0 
     endif
     do nm=n1,nmax+1
       n=m+nm
       fnum=real(n**2-m**2, rkindc)
       fden=real(4*n**2-1, rkindc)
       epsi(nm)=sqrt(fnum/fden)
     enddo
   enddo
!
   end subroutine rp_emns0
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine rp_bssl (bes,n)
!
!  Compute zeros of some Bessel functions used for determining Gaussian 
!  latitudes
!
      implicit none
!
      integer, intent(in)   :: n
      real(rkindc), intent(out) :: bes(n)
!
      real(rkindc)              :: bz(50)
      integer                   :: nn,j
      real(rkindc)              :: pi
      data pi/3.14159265358979/
      data bz / 2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391, &
   14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308, &
   27.4934791320, 30.6346064684, 33.7758202136, 36.9170983537, 40.0584257646, &
   43.1997917132, 46.3411883717, 49.4826098974, 52.6240518411, &
   55.7655107550, 58.9069839261, 62.0484691902, 65.1899648002, &
   68.3314693299, 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
   84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, & 
   96.6052679510, 99.7468198587, 102.8883742542, 106.0299309165, & 
   109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309, &
   121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752, &
   134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737, &
   146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/

      nn=min0(n,50)
      do 1 j=1,nn
      bes(j)=bz(j)
    1 continue
      if (n.le.50) return
      do 2 j=51,n
      bes(j)=bes(j-1)+pi
    2 continue
!
      end subroutine rp_bssl 
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine rp_gauss (a,w,k,iordlt)
!
! Compute Gaussian latitudes and weights used for Gaussian quadrature 
! on a sphere
!
      implicit none
!
      integer, intent(in)    :: k       ! number of latitudes
      integer, intent(in)    :: iordlt  ! ordering instruction
      real(rkindc), intent(out)  :: a(k)    ! sines of Gaussian latitudes
      real(rkindc), intent(out)  :: w(k)    ! Gaussian weights
!
! The weights are renormalized at the end of this routine such that their
! sum from pole to pole equals 1.  
! This is consistent with a change in normalization of Legendre Polynomials
! that are normalized such that the integral of their squared value from 
! pole to pole equals 1. This renormalization is consistent with a weighting
! that represents a fraction of the surface area on the sphere. 
!
!  pre  1981  NCAR       Initial algorithm acquired from ECMWF or BMRC?
!  ?????1981  R. Errico  Develped for normal mode software
!  04Jun2003  R. Errico  Initial algorithm at GMAO
!  04Oct2004  R. Errico  Module development
!
!-------------------------------------------------------------------------
!
!    local variables

!  iordlt=0 if order is north to south;
!  iordlt=1 if order is south to north;
!  iordlt=2 if order is alternating north, south, starting at poles.

      integer     :: kk,is,iter,n,l,n1,n2
      real(rkindc)    :: worka(k),facsin,eps,c,fk
      real(rkindc)    :: xz,pkm2,pkm1,fn,pk,pkmrk,sp,avsp
!     
      facsin=45./atan(one)
      if (rkindc == 4) then 
        eps=1.e-8
      else
        eps=1.d-14
      endif 
      c=(one-(two/3.14159265358979)**2)/(two*two)
      fk=k
      kk=k/2
      call rp_bssl (a,kk)

      do 3 is=1,kk
      xz=cos(a(is)/sqrt((fk+one/two)**2+c))
      iter=0
    1 pkm2=1.
      pkm1=xz
      iter=iter+1
      if (iter.gt.10) then
         print 100
         stop
      endif
      do 2 n=2,k
      fn=n
      pk=((two*fn-one)*xz*pkm1-(fn-one)*pkm2)/fn
      pkm2=pkm1
      pkm1=pk
    2 continue
      pkm1=pkm2
      pkmrk=(fk*(pkm1-xz*pk))/(one-xz**2)
      sp=pk/pkmrk
      xz=xz-sp
      avsp=abs(sp)
      if (avsp.gt.eps) goto 1
      a(is)=xz
!  original code
!           w(is)=(two*(one-xz**2))/(fk*pkm1)**2
!  code from M. Ehrendorfer:
      w(is)=two/((one-xz*xz)*pkmrk*pkmrk)

    3 continue

      if (k.ne.kk*2) then
         a(kk+1)=0.
         pk=2./fk**2
         do 4 n=2,k,2
         fn=n
         pk=pk*fn**2/(fn-one)**2
    4    continue
         w(kk+1)=pk
      endif

      do 5 n=1,kk
      l=k+1-n
! replaces sin lat by lat itself in degrees
!  a(n)=facsin*asin(a(n))
      a(l)=-a(n)
      w(l)=w(n)
      worka(n)=a(n)
      worka(l)=w(n)
    5 continue
      if (k.ne.kk*2) worka(kk+1)=w(kk+1)

      if (iordlt.eq.0) return
      if (iordlt.eq.1) then
         do 11 n=1,kk
         l=k+1-n
         a(n)=-worka(n)
         a(l)= worka(n)
         w(n)= worka(l)
         w(l)= worka(l)
   11    continue
         if (k.ne.kk*2) then
            a(kk+1)=0.
            w(kk+1)=worka(kk+1)
         endif
      elseif (iordlt.eq.2) then
         do 12 n=1,kk
         l=k+1-n
         n1=n*2
         n2=n1-1
         a(n2)= worka(n)
         a(n1)=-worka(n)
         w(n2)= worka(l)
         w(n1)= worka(l)
   12    continue
         if (k.ne.kk*2) then
            a(k)=0.
            w(k)=worka(kk+1)
         endif
      else
         print 101,iordlt
         stop
      endif
!
! Renormalize weights such that sum from pole to pole = 1
!
      w(:)=w(:)/two
!
  100 format(//,5x,'x x x x x  convergence failure in shgaus')
  101 format(//,5x,'x x x x x  bad option in shgaus      iordlt=',i4)
      end subroutine rp_gauss
!
!
  end module m_random_power
