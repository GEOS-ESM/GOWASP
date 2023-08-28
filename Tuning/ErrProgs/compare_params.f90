  Program Comapre_corr_params
!
! Visually compare estimated correlation parameters from several 
! instruments or kx
!  
  implicit none
!
  logical, parameter :: luvavg=.true.
  integer, parameter :: nfiles=11
  integer, parameter :: iunit=10
  integer, parameter :: kmax=1000   ! max nk permitted
!
  integer :: n,k
  integer :: nk
  integer :: lencf
  integer :: ic2
  integer :: kx(nfiles)
  integer :: icount(kmax,nfiles)
  real :: f2,h2,s2
  real :: f(kmax,nfiles)
  real :: hlength(kmax,nfiles)
  real :: stdv(kmax,nfiles)
  real :: plev(kmax)
  character(len=120) :: cpath
  character(len=120) :: cfiles(nfiles)
  character(len=240) :: cf
  character(len=40)  :: cdum(2)
!
  cpath='/discover/nobackup/rerrico/ODSstats/G514osse/Est_conv_July14'
  cfiles(1)='est_G514osse_conv_242_u.txt'
  cfiles(2)='est_G514osse_conv_243_u.txt'
  cfiles(3)='est_G514osse_conv_245_u.txt'
  cfiles(4)='est_G514osse_conv_246_u.txt'
  cfiles(5)='est_G514osse_conv_250_u.txt'
  cfiles(6)='est_G514osse_conv_252_u.txt'
  cfiles(7)='est_G514osse_conv_253_u.txt'
  cfiles(8)='est_G514osse_conv_254_u.txt'
  cfiles(9)='est_G514osse_conv_257_u.txt'
  cfiles(10)='est_G514osse_conv_258_u.txt'
  cfiles(11)='est_G514osse_conv_259_u.txt'
!
  do n=1,nfiles
    cf=trim(cpath)//'/'//trim(cfiles(n))
    open (iunit,file=trim(cf))
    read (iunit,'(a1)') cdum(1)
    read (iunit,'(54x,2i5)') nk,kx(n)
    print *,'n,file,nk,kx=',n,trim(cf),'  ',nk,kx(n)
    do k=1,7
      read (iunit,'(a1)') cdum(1)
    enddo
    do k=1,min(nk,kmax)
      read (iunit,'(4x,f12.2,i8,13x,e13.3,7x,f7.3,7x,f7.0)') &
          plev(k),icount(k,n),stdv(k,n),f(k,n),hlength(k,n)
    enddo
    close (iunit)
!
!  If file is for u or v and table is to have avg value for i and v, then ...
!
    lencf=len(trim(cf))     
    if (luvavg .and. (cf(lencf-5:lencf) == '_u.txt' .or. &
                      cf(lencf-5:lencf) == '_v.txt')) then
      if (cf(lencf-5:lencf) == '_u.txt') then 
        cf(lencf-5:lencf) = '_v.txt' 
      else
        cf(lencf-5:lencf) = '_u.txt' 
      endif
!
      open (iunit,file=trim(cf))
      read (iunit,'(a1)') cdum(1)
      read (iunit,'(54x,2i5)') nk,kx(n)
      print *,'n,file,nk,kx=',n,trim(cf),'  ',nk,kx(n)
      do k=1,7
        read (iunit,'(a1)') cdum(1)
      enddo
      do k=1,min(nk,kmax)
      read (iunit,'(4x,f12.2,i8,13x,e13.3,7x,f7.3,7x,f7.0)') &
        plev(k),ic2,s2,f2,h2
        icount(k,n)=0.5*(icount(k,n)+ic2)
        f(k,n)=0.5*(f(k,n)+f2)
        hlength(k,n)=0.5*(hlength(k,n)+h2)
        stdv(k,n)=0.5*(stdv(k,n)+s2)
      enddo
      close (iunit)
!
    endif
!
  enddo
!
  print *,' '
  print *,'Table of ic'
  print ('(13x,12i8)'),kx(1:nfiles)  
  do k=1,min(nk,kmax)
    print ('(i3,f10.2,12i8)'),k,plev(k),icount(k,:)
  enddo
!
  print *,' '
  print *,'Table of f'
  print ('(13x,12i8)'),kx(1:nfiles)  
  do k=1,min(nk,kmax)
    print ('(i3,f10.2,12f8.3)'),k,plev(k),f(k,:)
  enddo
!
  print *,' '
  print *,'Table of hlength'
  print ('(13x,12i8)'),kx(1:nfiles)  
  do k=1,min(nk,kmax)
    print ('(i3,f10.2,12f8.0)'),k,plev(k),hlength(k,:)
  enddo
!
  print *,' '
  print *,'Table of stdv'
  print ('(13x,12i8)'),kx(1:nfiles)  
  do k=1,min(nk,kmax)
    print ('(i3,f10.2,12f8.3)'),k,plev(k),stdv(k,:)
  enddo
!
  end program Comapre_corr_params
