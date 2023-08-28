    program time_histogram
!
    integer, parameter :: iunitin=11
    integer, parameter :: iunitout=12
    integer, parameter :: nlons=36
    integer, parameter :: nkxs=19
    integer, parameter :: ntimes=48
    integer :: ntx,kid
    integer :: nt,k
    integer :: kx(nkxs)
    integer :: satid(nkxs)
    integer :: nclon(nlons) 
    integer :: nctot(ntimes,nkxs)
    real(4) :: hour(ntimes)
    real(4), parameter :: dhour=0.5 
    character(len=17) :: char1
    character(len=7)  :: char2
    character(len=*), parameter  :: filein='satloc_osse.txt'
    character(len=*), parameter  :: fileout='time_series_osse.txt'
!
    do nt=1,ntimes
      hour(nt)=(nt-1)*dhour
    enddo
!  
    open (iunitin,file=trim(filein))
    read (iunitin,('(a)')) char2   ! skip header record
    do nt=1,ntimes
      do k=1,nkxs
        read (iunitin,('(a,4i6)')) char1,ntx,kid,kx(k),satid(k)
        read (iunitin,('(7x,12i6)')) nclon   
        nctot(nt,k)=sum(nclon)
        read (iunitin,('(7x,12i6)')) nclon   
        nctot(nt,k)=nctot(nt,k)+sum(nclon)
      enddo
    enddo
    close (iunitin)
!
    open (iunitout,file=trim(fileout))
    write (iunitout,('(11x,19i6)')) kx(:)
    write (iunitout,('(11x,19i6)')) satid(:)
    write (iunitout,('(a)')) ' '
    do nt=1,ntimes
      write (iunitout,('(i4,f7.1,19i6)')) nt,hour(nt),nctot(nt,:)
    enddo
    close (iunitout)
!
    end program time_histogram
