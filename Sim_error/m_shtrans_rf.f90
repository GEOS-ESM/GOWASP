   module m_shtrans_rf
!
!  Module for projecting spectral coefs onto scalar or vector fields 
!  (Used to create random fields)
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use m_kinds, only: rkind1, rkind2
   use m_parameters, only : earthr
!
   implicit none
!
   private
   public :: sh_init_factors
   public :: sh_init_lats
   public :: sh_clean
   public :: sh_trans
   public :: sh_calc_power
!
   integer, parameter :: rkindc=rkind2  ! kind for computational factors
   integer, parameter, public :: rkinds=rkind1  ! kind for spectral coefs 
   integer, parameter, public :: rkindf=rkind1  ! kind for output fields
   integer, parameter :: icoft=0        ! indicates ordering of spectral coefs
   integer, public :: nmax,mmax,kmax,imax,jmax,nspects
!
   integer, allocatable, public :: mtrunc(:),ntrunc(:),jindex(:,:)
   real(rkindc) :: zero,one,two,three,half
   real(rkindc), allocatable :: epsi(:,:)   
   real(rkindc), allocatable :: cosp(:),sinp(:)  ! values on lat grid
   real(rkindc), allocatable :: abfaclats(:,:)
   complex(rkindc), allocatable :: fm(:,:,:)  ! half_transformed variables
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_clean
!
!  Deallocate arrays used to define spectral truncation and FFT constants
!
   use m_fft, only: clean_fft
!
   deallocate (mtrunc,ntrunc,jindex)
   call clean_fft
!
   end subroutine sh_clean
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_init_factors (rf_nmax)
!
!  Call routines to determine lat-lon grid and all factors used to 
!  evaluate associated Legendre polynomials
!
   use m_fft, only: set_fft
!
   implicit none
!
   integer :: rf_nmax
!
   integer :: mnk,m1n1   
!
   zero=0.0_rkind2
   one=1.0_rkind2
   two=2.0_rkind2
   three=3.0_rkind2
   half=one/two
!
! Set spectral truncation to those for triangular truncation at wavenumber
! N=rf_nmax
   nmax=rf_nmax
   mmax=nmax
   kmax=nmax
!
!  Set grid dimensions 
   jmax=2*(nmax/2)+3  ! will include equator since odd
   call sh_imax      

! Compute number of spectral coefficients for general pentagonal truncation
   mnk=mmax+nmax-kmax
   m1n1=(mmax+1)*(nmax+1)
   nspects=m1n1-mnk*(mnk+1)/2
!
! Compute arrays describing truncation and spectral ordering
   allocate (mtrunc(0:nmax))
   allocate (ntrunc(0:mmax))
   allocate (jindex(0:mmax,0:nmax))
   call sh_trunc 
!
! Compute factorials used for computing Legendre polynomials
   allocate (epsi(0:nmax+1,0:mmax))
   call sh_emns 
!
! Compute factors required by FFT routines
   call set_fft (imax)
!
! Compute factors that depend on latitude
   call sh_init_lats 
!
   end subroutine sh_init_factors
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_init_lats 
!
!  Call routines to Determine latitude dependent factors used for
!  computing values of associated Legendre polynomials
!
   implicit none
!
   allocate (cosp(jmax),sinp(jmax))
   allocate (abfaclats(jmax,0:mmax))
!
! Set sines and cosines of latitude points (centers and edges)
   call sh_trig
!     
! Compute latitude and m factors used for computing associated 
! Legendre polynomials. 
   call sh_abfaclats 
!
   end subroutine sh_init_lats
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_imax 
!
!  Find imax that is the minimum integer that has factors of powers of 
!  only 2, 3, or 5 and that is greater than the max of 2*mmax+2 or 16. 
!
   implicit none
! 
   integer :: imin,n2,n3,n5,ig2,ig23,ig235,i1,i3,i5   
   real(rkind1) :: rmax
!
   imin=max((mmax+1)*2,16)
   rmax=real(imin)
   n2=2+int(log(rmax)/log(2.))
   n3=2+int(log(rmax)/log(3.))
   n5=2+int(log(rmax)/log(5.))
   imax=imin*2
   do i1=1,n2
     ig2=2**i1
     do i3=0,n3
       ig23=ig2*(3**i3)
       do i5=0,n5
         ig235=ig23*(5**i5)
         if (ig235 > imin .and. ig235 < imax) then
           imax=ig235
         endif
       enddo
     enddo
   enddo    
!    
   end subroutine sh_imax 
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_abfaclats 
!
!  Compute some latitude-dependent factors used to compute associated 
!  Legendre polynomials 
!
   implicit none
!
   integer :: j,m
   real (rkindc) :: cos2,a,b,prod
!
   abfaclats(:,0)=one
   abfaclats(1,1:mmax)=zero      ! value for pole point
   abfaclats(jmax,1:mmax)=zero   ! value for pole point
   do j=2,jmax-1
     cos2=one-sinp(j)*sinp(j)
     prod=one
     a=one
     b=zero
     do m=1,mmax      
       a=a+two
       b=b+two
       prod=prod*cos2*a/b
       abfaclats(j,m)=sqrt(prod)  
     enddo    ! loop over m
   enddo      ! loop over j   
!
   end subroutine sh_abfaclats
!
! 
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_alp (nlast,m,j,name,alp)
!
!  Compute associated Legendre polynomials for given m and latitude.
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
   integer, intent(in) :: m
   integer, intent(in) :: j
   real(rkindc), intent(out) :: alp(0:nlast)
   character(len=*), intent(in) :: name
!
   integer :: nm  ! = n-m
   integer :: np  ! indicator for order of alp storage 
   real(rkindc) :: fm  ! = real(m)  
!
   fm=real(m,rkindc)
!
! Compute values for n-m=0 (special treatment of vector fields at pole)
   if (name /= 'S' .and. (j == 1 .or. j == jmax)) then
     if (m /= 1) then
       alp(:)=zero
     else
       alp(0)=sqrt(three/two)
     endif 
   else    
     alp(0)=abfaclats(j,m)
   endif
!
! Compute values for n-m=1
   if (nlast > 0) then
     alp(1)=sqrt(two*fm+three)*sinp(j)*alp(0)
   endif
!
! Compute values for n-m>1 using recursion
   if (nlast > 1) then
     do nm=2,nlast
       alp(nm)=(sinp(j)*alp(nm-1)-epsi(nm-1,m)*alp(nm-2))/epsi(nm,m)
     enddo
   endif
!
   end subroutine sh_alp
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_emns 
!
!  Compute some latitude independent factors used to compute associated
!  Legendre polynomials
!
   implicit none
!
   integer :: m,n,n1,nm
   real(rkindc) :: fnum,fden
!
   do m=0,mmax
     if (m == 0) then
       epsi(0,0)=zero    
       n1=1
     else
       n1=0 
     endif
     do nm=n1,nmax+1
       n=m+nm
       fnum=real(n**2-m**2, rkindc)
       fden=real(4*n**2-1, rkindc)
       epsi(nm,m)=dsqrt(fnum/fden)
     enddo
   enddo
!
   end subroutine sh_emns
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_trunc 
!
!  Determine ordering of associated Legendre polynomials for general
!  pentagonal trunction (of which triangular and rhomboidal are special 
!  examples)
!
   implicit none
!
   integer :: j, ksubm, ksubn, m, n
!
   ksubm=kmax-mmax
   ksubn=kmax-nmax
!
   do n=0,nmax
     mtrunc(n)=mmax
     if (n > ksubm) mtrunc(n)=mtrunc(n-1)-1
   enddo
!
   do m=0,mmax
     ntrunc(m)=nmax
     if (m > ksubn) ntrunc(m)=ntrunc(m-1)-1
   enddo
!
!  specify ordering of modes
!     icoft=0 is ordering along vertical lines of m,n diagram
!     icoft=1 is ordering along diagonal lines of m,n diagram
!
   if (icoft == 0) then
     j=0
     do m=0,mmax
       do n=0,ntrunc(m)
         j=j+1
         jindex(m,n)=j
       enddo
     enddo
   elseif (icoft == 1) then  
     j=0
     do n=0,nmax
       do m=1,mtrunc(n)
         j=j+1
         jindex(m,n)=j
       enddo
     enddo
   endif
!
   end subroutine sh_trunc
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_trig 
!
! Determine sines and cosines at edges of latitude bands
!
   implicit none
!
   integer :: j,jmax1
   real(rkindc) :: pi_fac,pi_lat,dp,half_pi
!
   jmax1=jmax-1
   half_pi=two*atan(one)
   pi_fac=two*half_pi
   dp=pi_fac/real(jmax1,rkind2)  ! latitudinal spacing in radians 
!
! Define sin at edges of lat bands
   sinp(1)=-one
   sinp(jmax)=one
   do j=2,jmax1
     pi_lat=-half_pi+real(j-1,rkind2)*dp
     sinp(j)=sin(pi_lat)
   enddo
!
   cosp(1)=zero
   cosp(jmax)=zero
   do j=2,jmax1
     cosp(j)=sqrt(one-sinp(j)*sinp(j))
   enddo
!
   end subroutine sh_trig
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_trans (nlevs,jfirst,jlast,scoefs,stypes,f,name)
!
!  Sequence routines to transform from spherical harmonics to global 
!  scalar or vector fields
!
   use m_fft, only: run_fft
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: jfirst,jlast
   integer, intent(in) :: stypes
   real(rkindf), intent(out) :: f(imax,jfirst:jlast,nlevs)
   complex(rkinds), intent(in) :: scoefs(nspects,nlevs,stypes)
   character(len=*), intent(in) :: name
!
   integer :: imax2p1,jlats
   integer :: i,j,k,m,j1
   real(rkindc), allocatable :: fxj(:,:)  
   complex(rkindc), allocatable :: fmj(:,:)  
!
   imax2p1=imax/2+1
   jlats=jlast-jfirst+1
!
   allocate (fm(0:mmax,jfirst:jlast,nlevs))
   fm(:,:,:)=cmplx(zero,zero,rkindc)
!
! Perform half transform from coefs to zonal waves on lats 
   call sh_s2fm (nlevs,jfirst,jlast,scoefs,stypes,name)
!
! If wind field, divide by coslat
   if (name /= 'S') then 
     call sh_uvcos (nlevs,jfirst,jlast)
   endif
!
! Do zonal half-transform on all latitudes and levels
   allocate (fmj(imax2p1,jlats))
   allocate (fxj(imax,jlats))
   do k=1,nlevs
!
     do j=1,jlats
       j1=jfirst+j-1
       do m=1,mmax+1
         fmj(m,j)=fm(m-1,j1,k)
       enddo
       if (mmax+1 < imax2p1) then 
         do m=mmax+2,imax2p1
           fmj(m,j)=cmplx(zero,zero,rkindc)
         enddo
       endif
     enddo
!   
     call run_fft (imax,imax,imax2p1,jlats,fxj,fmj)
!
     do j=1,jlats
       j1=jfirst+j-1
       do i=1,imax
         f(i,j1,k)=fxj(i,j)
       enddo
     enddo
!
   enddo  ! loop over levels
   deallocate (fxj,fmj,fm)
!
   end subroutine sh_trans 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_s2fm (nlevs,jfirst,jlast,scoefs,stypes,name)
!
! Half transform from spherical harmonic coefficients to zonal Fourier 
! coefficients
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: jfirst,jlast
   integer, intent(in) :: stypes
   complex(rkinds), intent(in) :: scoefs(nspects,nlevs,stypes)
   character(len=*), intent(in) :: name
!
   integer :: j,m
   integer :: nlast ! max of n-m given n (+1 if vector field)
   real(rkindc), allocatable :: alp(:)   ! Legendre polynomials for 1 m 
   complex(rkindc), allocatable :: cwork(:,:)
!
! Loop over zonal wave numbers m
   do m=0,mmax
     nlast=ntrunc(m) 
     if ( name /= 'S') then ! truncation required for vector field
       nlast=nlast+1
     endif
!
     allocate (alp(0:nlast))
     allocate (cwork(0:nlast,nlevs))   
!
! If scalar fields, copy spectral coefficients for one m to an array
! If vector fields, compute spectral coefficients for one m from 
! spectral coeficients of vorticity and divergence
     if (name == 'S') then
       call sh_s2c_s (nlevs,m,nlast,scoefs,cwork)
     else 
       call sh_s2c_v (nlevs,m,nlast,scoefs,cwork,name)
     endif
!
! Project polynomials onto latitudes for 1 m  
! First Compute alp for one m but all n and latitudes

     do j=jfirst,jlast
       call sh_alp (nlast,m,j,name,alp)
       call sh_c2fm (nlevs,m,nlast,j,alp,cwork)
     enddo
!
     deallocate (alp)
     deallocate (cwork)
!
   enddo ! loop over m
!
   end subroutine sh_s2fm
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_s2c_s (nlevs,m,nlast,scoefs,cwork)
!
! Copy spectral coefficients for a scalar for one value of m
!
   implicit none
!
   integer, intent(in) ::  nlevs
   integer, intent(in) ::  m
   integer, intent(in) ::  nlast
   complex(rkinds), intent(in) :: scoefs(nspects,nlevs)
   complex(rkindc), intent(out)  :: cwork(0:nlast,nlevs)
!
   integer :: n,i,nf
!
   do n=0,nlast
     i=jindex(m,n)
     do nf=1,nlevs           
       cwork(n,nf)=scoefs(i,nf)
     enddo
   enddo
!
   end subroutine sh_s2c_s
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_s2c_v (nlevs,m,nlast,scoefs,cwork,name)
!
!  Compute spherical harmonic coefficients for scalar u*cos(lat) or 
!  v*cos(lat) fields from spherical harmonic coeficients for stream 
!  function and velocity potential 
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: m
   integer, intent(in) :: nlast
   complex(rkinds), intent(in)  :: scoefs(nspects,nlevs,2)
   complex(rkindc), intent(out) :: cwork(0:nlast,nlevs)
   character(len=*), intent(in)  :: name
!
   integer :: n,i,nf
   integer :: truen  ! the true index n
   integer :: id1, id2, is
   real(rkindc) :: rfacm, rfac
   complex(rkindc) :: cfacm
!
   rfacm=-real(m,rkindc)*earthr
   if (name == 'U') then
     id1=2
     id2=1
     is=1
   else      ! name='V' assumed
     id1=1
     id2=2
     is=-1
   endif
!
   do n=0,nlast-1
     truen=n+m
!
     i=jindex(m,n)
     if (truen == 0) then
       cfacm=cmplx(zero,zero,rkindc)
     else
       cfacm=cmplx(zero,rfacm,rkindc)/((truen+1)*truen)
     endif
     do nf=1,nlevs           
       cwork(n,nf)=cfacm*scoefs(i,nf,id1)
     enddo
!      
     if ( n > 0) then
       i=jindex(m,n-1)
       rfac=-is*epsi(n,m)*earthr/truen
       do nf=1,nlevs           
         cwork(n,nf)=cwork(n,nf)+rfac*scoefs(i,nf,id2)
       enddo
     endif  
!      
     if ( n < nlast -1 ) then
       i=jindex(m,n+1)
       rfac=is*epsi(n+1,m)*earthr/(truen+1)
       do nf=1,nlevs           
         cwork(n,nf)=cwork(n,nf)+rfac*scoefs(i,nf,id2)
       enddo
     endif  
!
   enddo ! loop over n
!
   i=jindex(m,nlast-1)
   truen=nlast+m
   rfac=-is*epsi(nlast,m)*earthr/truen
   do nf=1,nlevs           
     cwork(nlast,nf)=rfac*scoefs(i,nf,id2)
   enddo
!
   end subroutine sh_s2c_v
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!    
   subroutine sh_c2fm (nlevs,m,nlast,j,alpj,cwork)
!
!  Perform half transform from spherical-harmonic coefficients for a scalar
!  field to Fourier coefficients for each grid latitude
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: m
   integer, intent(in) :: nlast
   integer, intent(in) :: j
   real(rkindc), intent(in) :: alpj(0:nlast)
   complex(rkindc), intent(in)  :: cwork(0:nlast,nlevs)
!
   integer :: n,nf
!
   do nf=1,nlevs
     fm(m,j,nf)=cmplx(zero,zero,rkindc)   
     do n=0,nlast
       fm(m,j,nf)=fm(m,j,nf)+alpj(n)*cwork(n,nf)
     enddo
   enddo
!
   end subroutine sh_c2fm 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_uvcos (nlevs,jfirst,jlast)
!
!  Divide u*cos(lat) and v*cos(lat) by cos(lat)
!
   implicit none
!
   integer, intent(in) :: nlevs,jfirst,jlast
!
   integer :: j,m,k
   real(rkindc) :: fac
!
   do j=max(jfirst,2),min(jlast,jmax-1)
     fac=one/cosp(j)
     do m=0,mmax
       do k=1,nlevs
         fm(m,j,k)=fm(m,j,k)*fac
       enddo
     enddo
   enddo
!
   end subroutine sh_uvcos
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sh_calc_power (nlevs,nfields,scoefs,power,field_type)
!
! Routine for calc of power spectra from scoefs
!
   implicit none
!
   integer,  intent(in)  :: nlevs  ! number of horiz surfaces
   integer,  intent(in)  :: nfields ! number of distinct field types 
   real(rkind1), intent(out) :: power(0:kmax,nlevs,nfields) 
   complex(rkinds), intent(in) :: scoefs(nspects,nlevs,nfields)
   character(len=*), intent(in)  :: field_type ! value of S or V 
!
   logical :: lvector
   integer  :: nfield,nlev,nm,m1,n,m,j
   real(rkind1) :: aa,afac,xmfac,cr,ci
!
!  Compute power spectra from spectral coefficients, by summing over all
!  zonal wavenumbers m for each value of n. Formally, the power spectra 
!  are determined by summing over -m,...,0,...,m  but since thhe spectral
!  coefficients for m<0 are complex conjugates of those m >= 0, the 
!  summation need only be done for non-negative values, but for 
!  contributions by all m /= 0, a factor of 2 is introduced. If the 
!  spectral coefficients are identified as determined from a vector 
!  field (by specifying input name='V'), then the power is determined as 
!  if for kinetic energy from spectral coefficients of vorticity and 
!  divergence (except for the factor 1/2 in the definition of energy, 
!  as in KE = 1/2 * (u**2 +v**2)).  If spectral coefficients of vorticity 
!  are input, but the field is indicated as a scalar one (name='S' input)  
!  instead, then enstropy spectra rather than kinetic energy spectra 
!  will be calculated.
!
!  Note that if 2*mmax=imax, then coefs for zonal wavenumber -max are
!  identical to those for mmax, and should not be counted separately
!  (i.e., the contribution by this zonal wavenumber to the power 
!  should not be multiplied by 2).  This value of mmax, however, is not
!  permitted if the vector wind is a transformed field. 
!
!  Note that the power spectra is only computed for n=0,...,nmax.
!  Thus, if thespectral truncation is not triangular, power that resides
!  in other spectral components will not be counted.
!
   aa=earthr*earthr
   lvector=(field_type == 'V') 
   power(:,:,:)=0.
!
   do nfield=1,nfields
     do nlev=1,nlevs
       do nm=0,nmax      ! n-m
!
! Compute power for vector (scoefs of vort and divg) fields 
         if (lvector) then 
           if (nm == 0) then
             m1=1              ! no power in mode n=m=0 for vectors 
           else
             m1=0
           endif      
           do m=m1,mtrunc(nm) 
             j=jindex(m,nm)
             n=nm+m         
             afac=aa/(n*(n+1)) ! change e.g., enstropy to energy
             if (m == 0) then  ! account for absence of m=-0
               xmfac=1.
             elseif (2*m == imax) then 
               xmfac=1.        ! account for redundancy of m=-mmax
             else
               xmfac=2.        ! account for m=-m conjugate waves 
             endif
             cr=real( scoefs(j,nlev,nfield))
             ci=aimag(scoefs(j,nlev,nfield))
             power(n,nlev,nfield)=power(n,nlev,nfield)+afac* &
                  xmfac*(cr*cr+ci*ci)
           enddo   
!
!  Compute power for scalar fields            
         else  
           do m=0,mtrunc(nm)
             j=jindex(m,nm)                   
             n=nm+m
             if (m == 0) then
               xmfac=1.         ! account for absence of m=-0
             elseif (2*m == imax) then 
               xmfac=1.         ! account for redundancy of m=-mmax
             else
               xmfac=2.         ! account for m=-m conjugate waves 
             endif
             cr=real( scoefs(j,nlev,nfield))
             ci=aimag(scoefs(j,nlev,nfield))
             power(n,nlev,nfield)=power(n,nlev,nfield)+ &
                  xmfac*(cr*cr+ci*ci)
           enddo               
!
         endif  ! test on lvector
       enddo    ! end loop over nm
     enddo      ! end loop over nlev
   enddo        ! end loop over nfield
!
   end subroutine sh_calc_power
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_shtrans_rf

