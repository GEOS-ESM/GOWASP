      module m_fft
!
! Calls either a slow or fast Fourier transform from Fourier coefficients to   
! grid-point values. Length of field should be at least one factor of 2.
! Further restrictions on length of field apply to fast transform.
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
      implicit none
      private
      public set_fft
      public run_fft
      public clean_fft
      public rkind_fft
      public imax2
      integer, parameter :: rkind_fft=8
      integer, parameter :: ifirst=1  ! 0 if 1st pt is x=0; 1 if it is x=1/imax
! ifirst=1 is used only to make results agree with previous scheme
      integer :: imax2
      real(rkind_fft), parameter :: zero=0._rkind_fft
      real(rkind_fft), parameter :: one=1._rkind_fft
      real(rkind_fft), parameter :: two=2._rkind_fft
      character(len=4), parameter :: which_fft='F991'  ! 'SLOW' or 'F991'
!  FFT99 and FFT991 are 1980 ECMWF F77 routines
      real(rkind_fft), allocatable :: fft_trigs(:,:,:)  
      integer :: fft_ifax(13) ! used for FFT99
!
      contains
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine set_fft (imax)
!
! Set some trigonometric factors used in the Fourier transforms
!
      integer :: imax
      integer :: ndim
      if (which_fft == 'F991') then       ! use 1980 ECMWF FFT
        ndim=3*imax/2+1
        allocate (fft_trigs(ndim,1,1)) 
        call set99 (fft_trigs(:,1,1),fft_ifax,imax)      
      elseif  (which_fft == 'SLOW') then  ! use slow transform
        imax2=imax/2                   
        if (imax2*2 /= imax) then
          imax2=imax2+1  ! coefs span k=0,...,ndim
        endif 
        allocate (fft_trigs(2,0:imax2,imax))
        call sh_ft_slow_init (imax)
      endif
      end subroutine set_fft
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine clean_fft
!
! Deallocate array of trigonometric factors used by FFT
!
      if (allocated(fft_trigs)) then
        deallocate (fft_trigs)
      endif
      end subroutine clean_fft
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine run_fft (imax,ldx,ldy,nums,x,y)
!
! Call routine to transform from Fourier coefficients to grid values
!
      integer :: imax
      integer :: nums
      integer :: ldx, ldy
      real(rkind_fft) :: x(ldx,nums) 
      complex(rkind_fft) :: y(ldy,nums) 
!
      integer :: i,i1,i2,j,k
      integer :: jump, inc, ndim
      real(rkind_fft), allocatable :: f(:)
      real(rkind_fft), allocatable :: work(:)
!
      if (which_fft == 'SLOW') then
        call sh_ft_slow_m2f (imax,nums,y,ldy,x,ldx)
!
      elseif (which_fft == 'F991') then
        inc=1
        jump=ldx+2
        ndim=nums*(ldx+2)
        allocate (f(ndim))
        ndim=nums*(imax+1)
        allocate (work(ndim))
!
        do j=1,nums  
          do k=1,ldy
            i2=k*2+jump*(j-1)
            i1=i2-1
            f(i1)=real(y(k,j))
            f(i2)=aimag(y(k,j))
          enddo
        enddo
!
        call fft991 (f,work,fft_trigs(:,1,1),fft_ifax,inc,jump, &
                     imax,nums,1)
!      
        if (ifirst == 0) then ! first point location is x=0
          do j=1,nums
            do k=1,ldx
              i1=k+jump*(j-1)
              x(k,j)=f(i1)
            enddo
          enddo
        else  ! first point location for x is shifted compared to FFT991
          do j=1,nums
            do k=1,ldx-1
              i1=k+jump*(j-1)+1        ! shift by 1
              x(k,j)=f(i1)
            enddo
            x(ldx,j)=f(1+jump*(j-1))   ! given periodicity
          enddo
        endif
!
        deallocate (f,work)
!
      else
        print *,' '
        print *,' ERROR IN M_FFT PARMATER:'
        print *,' which_fft=',which_fft,' not a valid option'
        stop
!
      endif  ! test of which_fft
!
      end subroutine run_fft 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine sh_ft_slow_init (imax)
!
! Set trig factors for slow Fourier trans.
!
      implicit none
!
      integer,  intent(in) :: imax
! 
      integer             :: m      ! integral zonal wavenumber
      integer             :: j      ! zonal grid point index
      integer             :: jshift 
      real(rkind_fft)            :: pifac  ! 2*pi/imax
      real(rkind_fft)            :: pms    ! k*x
!      
! Pre-compute sines and cosines of k*x;  k=m*2*pi ; x=j/imax
! j=1,...,imax; m=0,..., imax2=imax/2 (+1 if imax is odd)
! for use in slow Fourier transform routine.
!
      pifac=(two**3)*datan(one)/imax
      if (ifirst == 0) then ! first point location is x=0
        jshift=1
      else                  ! first point location is x=1/imax
        jshift=0
      endif
!
      do m=0,imax2
        pms=pifac*m
        do j=1,imax
          fft_trigs(1,m,j)=cos(pms*(j-jshift))
          fft_trigs(2,m,j)=sin(pms*(j-jshift))
        enddo
      enddo
!
      end subroutine sh_ft_slow_init 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
      subroutine sh_ft_slow_m2f (imax,lot,y,ldy,x,ldx)
!
! A slow Fourier trans.: zonal coefs to fields
!
      integer,  intent(in) :: imax
      integer,  intent(in) :: lot
      integer,  intent(in) :: ldx
      integer,  intent(in) :: ldy
      complex(rkind_fft),  intent(in) :: y(0:ldy-1,lot)
      real(rkind_fft),intent(out) :: x(0:ldx-1,lot)
!
      integer             :: m,j,n
      real(rkind_fft)            :: sfac
      complex(rkind_fft)         :: wjm
!
! A slow algorithm to compute zonal transforms: from complex zonal wave 
! coefficients to real fields. 
! An FFT that allows imax with factors of 2*(2**m)*(3**n)*(5**k) should
! be called in place of this slow algorithm. 
!
      x=zero
      sfac=one
!
      do m=0,ldy-1
        do j=0,ldx-1
          wjm=cmplx(fft_trigs(1,m,j+1),fft_trigs(2,m,j+1), rkind_fft)
          do n=1,lot
            x(j,n)=x(j,n)+sfac*real(wjm*y(m,n))
          enddo
        enddo
        sfac=two ! factor applied when m>0
      enddo
!
      end subroutine sh_ft_slow_m2f 
!
!
      end module m_fft
