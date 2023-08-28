!
    subroutine read_hcovs (max_nbins,nk,n_count,file_cov,  &
                           nlevs,nbins,ierr,xbin,          &
                           ibvalues,pvalues,binvalues)
!
!  Read stats previously written to file produced by hcorr_tables program
!
    use m_kinds, only : rkind1
    implicit none
!
    integer, intent(in)  :: max_nbins
    integer, intent(in)  :: nk
    integer, intent(out) :: n_count  ! number of assimilation periods averaged
    integer, intent(out) :: nlevs
    integer, intent(out) :: nbins
    integer, intent(out) :: ierr
    integer, intent(out) :: ibvalues(max_nbins)
    real(rkind1), intent(out) :: xbin
    real(rkind1), intent(out) :: pvalues(nk)
    real(rkind1), intent(out) :: binvalues(max_nbins,nk,3)
    character(len=*), intent(in) :: file_cov
!
    integer, parameter :: iunit=10
    integer, parameter :: ikind8=8
    integer, parameter :: rkind8=8
    integer :: minlevs
    integer :: n_ksmax,n_date,n_time,kt_look,kx_look
    integer(ikind8), allocatable :: iaccept(:) 
    real(rkind8) :: x_max,x1,x_bin,latN,latS,lonW,lonE,delp
    real(rkind8), allocatable :: xmean(:), xstdv(:), plevs(:)
    real(rkind8), allocatable :: xbins(:,:)  
    real(rkind8), allocatable :: corbins(:,:)  
    real(rkind8), allocatable :: covbins(:,:)  
    character(len=120) :: c_file
!
    open (unit=iunit,file=trim(file_cov),form='unformatted')
    print *,'Hcovs file opened: ',trim(file_cov)
    read (iunit) n_count,n_ksmax,n_date,n_time,nbins,x_max, &
              x1,x_bin,latN,latS,lonW,lonE,                 &
              delp,kt_look,kx_look,nlevs,c_file
    xbin=x_bin
    print ('(a)'),'n_count,n_ksmax,n_date,n_time,nbins,kt_look,kx_look,nlevs:'
    print ('(8i10)'),n_count,n_ksmax,n_date,n_time,nbins,kt_look,kx_look,nlevs
    print ('(a)'),'x_max,x1,xbin,latN,latS,lonW,lonE,delp'
    print ('(8f10.2)'),x_max,x1,x_bin,latN,latS,lonW,lonE,delp
    allocate (iaccept(nk))
    allocate (xmean(nk))
    allocate (xstdv(nk))
    allocate (plevs(nlevs))
    allocate (xbins(nbins,nk))
    allocate (corbins(nbins,nk))
    allocate (covbins(nbins,nk))
    read (iunit) iaccept,xmean,xstdv,plevs
    read (iunit) xbins
    read (iunit) corbins
    read (iunit) covbins
    close (iunit)     
!
    if (nbins > max_nbins .or. nk /= n_ksmax) then
      ierr=100

      print *,'ERROR: either nbins>max_nbins or nk/=n_ksmax'
      print ('(a,5i4)'),'nbins,max_nbins,nk,n_ksmax,nlevs= ', &
                         nbins,max_nbins,nk,n_ksmax,nlevs
    else
      ierr=0
      minlevs=min(nk,nlevs)  ! for radiances, pvalues array may be smaller
      ibvalues(1:nk)=iaccept(:)
      pvalues(1:minlevs)=plevs(1:minlevs)
      binvalues(1:nbins,1:nk,1)=covbins(1:nbins,1:nk)
      binvalues(1:nbins,1:nk,2)=corbins(1:nbins,1:nk)
      binvalues(1:nbins,1:nk,3)=xbins(1:nbins,1:nk)*n_count
    endif
!
    deallocate (iaccept,xmean,xstdv,plevs)
    deallocate (xbins,corbins,covbins)
!
    end subroutine read_hcovs 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine read_counts (nk,n_count,file_in,icounts,plevs,stdvs)
!
!  Read stats previously written to file produced by counobs_table
!
    use m_kinds, only : rkind1 
    implicit none
!
    integer, intent(in) :: nk
    integer, intent(out) :: n_count  ! number of assimilation periods averaged
    integer, intent(out) :: icounts(nk)  ! avg num obs per assim period
    real(rkind1), intent(out) :: stdvs(nk)
    real(rkind1), intent(out) :: plevs(nk)
    character(len=*), intent(in) :: file_in
!
    integer, parameter :: iunit=10
    integer :: n
    character(len=1) :: cdum
!
! Open input file 
    open (unit=iunit,file=trim(file_in),form='formatted')
!
! Skip header except for n_count
    read (iunit,'(a)') cdum
    read (iunit,'(10x,i3)') n_count
    do n=1,8
      read (iunit,'(a)') cdum
    enddo
!
    do n=1,nk
      read (iunit,'(5x,f12.2,5x,i8,36x,f12.6)') plevs(n),icounts(n),stdvs(n) 
    enddo
!
    close (iunit)
    print *,'Counts and stdvs read from ',trim(file_in)
!
    end subroutine read_counts
!   
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine read_chcovs (nk,n_count,file_cov,ierr,ibvalues,binvalues)
!
! Read stats previously written to file produced by chcorr_tables program
!
    use m_kinds, only : rkind1
    implicit none
!
    integer, intent(in)  :: nk
    integer, intent(out) :: ierr
    integer, intent(out) :: n_count  ! number of assimilation periods averaged
    integer, intent(out) :: ibvalues(nk)
    real(rkind1), intent(out) :: binvalues(nk,nk,3)
    character(len=*), intent(in) :: file_cov
!
    integer, parameter :: iunit=10
    integer, parameter :: ikind8=8
    integer, parameter :: rkind8=8
    integer :: n_date,n_time,nbins,count_used
    integer(ikind8), allocatable :: iaccept(:) 
    real(rkind8) :: x_max,x1,x_bin,latN,latS,lonW,lonE,delp
    real(rkind8), allocatable :: xmean(:), xstdv(:), plevs(:)
    real(rkind8), allocatable :: xbins(:,:)  
    real(rkind8), allocatable :: corbins(:,:)  
    real(rkind8), allocatable :: covbins(:,:)  
    character(len=120) :: c_file
!
    open (unit=iunit,file=trim(file_cov),form='unformatted')
    print *,'Chcovs file opened: ',trim(file_cov)
    read (iunit) n_count,n_date,n_time,nbins,count_used, &
                 latN,latS,lonW,lonE
    print ('(a,5i10)'),'n_count,n_date,n_time,nbins,count_used:', &
                        n_count,n_date,n_time,nbins,count_used
    print ('(a,4f10.2)'),'latN,latS,lonW,lonE',latN,latS,lonW,lonE
    allocate (iaccept(nbins))
    allocate (xmean(nbins))
    allocate (xstdv(nbins))
    allocate (xbins(nbins,nbins))
    allocate (corbins(nbins,nbins))
    allocate (covbins(nbins,nbins))
    read (iunit) iaccept,xmean,xstdv,plevs
    read (iunit) xbins
    read (iunit) corbins
    read (iunit) covbins
    close (iunit)     
!
    if (nbins /= nk) then
      ierr=100
      print *,'ERROR: read chcorr nbins /= program nk =',nk
    else
      ierr=0
      ibvalues(1:nk)=iaccept(:)
      binvalues(1:nbins,1:nk,1)=covbins(1:nbins,1:nk)
      binvalues(1:nbins,1:nk,2)=corbins(1:nbins,1:nk)
      binvalues(1:nbins,1:nk,3)=xbins(1:nbins,1:nk)*n_count
    endif
!
    deallocate (iaccept,xmean,xstdv)
    deallocate (xbins,corbins,covbins)
!
    end subroutine read_chcovs 
