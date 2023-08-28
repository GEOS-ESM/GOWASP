      program countobs_tables
!
! Determine tables of channel counts, means, stdvs, and lat distributions
! from file of previously accumulated sums.  
!
! At the end of the program, for each level or channel, icount is the average 
! obs count per assimilation period. 
!
      implicit none
!
      integer, parameter :: rkind=8
      integer, parameter :: r4=4           ! r*4 so that grads can use these
      integer, parameter :: ikind=8
      integer, parameter :: n_c_select=15  
      integer(ikind), allocatable :: iaccept(:)
      integer(ikind), allocatable :: ibins(:,:)
      integer :: ib
      integer :: count_used  ! number of channels actually used
      integer :: i1, i2
      integer :: i
      integer :: i_list(n_c_select)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(n_c_select)
      integer :: c_select(n_c_select)
      integer :: n
      integer :: n_bins   ! number of channels on file
      integer :: n_count
      integer :: n_date
      integer :: n_time
      integer :: nplevs
      integer :: argc
      integer :: nlats
      integer :: kt,kx
      integer :: count_sum
!
      real(rkind), parameter :: zero=0._rkind
      real(rkind) :: dlat, xlat
      real(rkind) :: xfac ! change from 1/n for mean to 1/(n-1) for var
      real(rkind) :: latN, latS, lonW, lonE
      real(rkind) :: xa, xb, xm1, xm2, xs1, xs2
      real(rkind) :: z_gps(3)
      real(rkind) :: x_values(n_c_select)
      real(rkind), allocatable :: plevs(:)
      real(rkind), allocatable :: aoverb(:)
      real(rkind), allocatable :: xsum(:,:)
      real(rkind), allocatable :: xsumsq(:,:)
      real(rkind), allocatable :: xmean(:,:)
      real(rkind), allocatable :: xmsqr(:,:)
      real(rkind), allocatable :: xstdv(:,:)
!
      character(len=4) :: c_surf
      character(len=220) :: c_file
      character(len=220) :: file_in
      character(len=220) :: file_out
      character(len=224) :: file_out_1
!
! This is not set properly:
      do n=1,n_c_select
        c_select(n)=120+n
      enddo
!
      argc = iargc()
      if (argc .ne. 2) then
        print *,' usage must be: vcorr_tables.x filein fileout'
        stop
      endif
      call GetArg( 1_4, file_in)
      call GetArg( 2_4, file_out)
      file_out_1=trim(file_out)//'.txt'
!
! read file of previously accumulated values
      open (unit=10,file=trim(file_in),form='unformatted')
      read (10) n_count,n_bins,n_date,n_time,nlats, &
                latN,latS,lonW,lonE,kt,kx,c_surf,c_file
      allocate (plevs(n_bins))
      allocate (xsum(n_bins,3))
      allocate (xsumsq(n_bins,3))
      allocate (iaccept(n_bins))
      allocate (ibins(nlats,n_bins))
      read (10) plevs
      read (10) iaccept,xsum,xsumsq
      read (10) ibins 
      close (10)     
      print *,'Input file read: ',trim(file_in)
!
      allocate (xmean(n_bins,3))
      allocate (xmsqr(n_bins,3))
      allocate (xstdv(n_bins,3))
      allocate (aoverb(n_bins))
      xmean(:,:)=zero
      xstdv(:,:)=zero
      aoverb(:)=zero
!
! Reset c_select if not AIRS
      if (n_bins < 183) then  ! 183 is bounds because max(c_select)=183 as set
        do ib=1,min(n_bins,n_c_select)
          c_select(ib)=ib
        enddo
      endif
!
! Compute means and standard deviations
      do ib=1,n_bins
        if (iaccept(ib) > 0) then
          if (iaccept(ib) > 1) then
            xfac=real(iaccept(ib))/real(iaccept(ib)-1)
          else  
            xfac=0.
          endif
!
          do n=1,3
            xmean(ib,n)=xsum(ib,n)/iaccept(ib)
            xmsqr(ib,n)=xfac*(xsumsq(ib,n)/iaccept(ib) - &
                       xmean(ib,n)*xmean(ib,n))
            if (xmsqr(ib,n) > zero) then
              xstdv(ib,n)=sqrt(xmsqr(ib,n))
            endif
          enddo
          if (xstdv(ib,1) > zero) then
            aoverb(ib)=xstdv(ib,2)/xstdv(ib,1)
          endif
        endif
      enddo
!
! Replace total number of obs by average number per synoptic period
      iaccept(:)=iaccept(:)/n_count
      print ('(a,4i10)'),'Header info: ',n_count,n_bins,n_date,n_time
!
! Determine parameters for printing a subset of distributions
      i=1
      do ib=1,n_bins
        if (ib == c_select(i)) then
          i_list(i)=ib
          i_list_last=i
          if (i<n_c_select) then
            i=i+1
          endif  
        endif
      enddo
!
! Open output file 
      open (unit=10,file=trim(file_out_1),form='formatted')
      print *,' '
      print *,'Output file opened: ',trim(file_out_1)
!
! Write header information
      write (10,'(2a)') 'Channel Counts and Statistics'
      write (10,'(a,i3,a)') 'Data from ',n_count,' times'     
      write (10,'(a,2i12)') 'Last date and time: ',n_date,n_time
      write (10,'(a,a)')    'Last file of data: ',trim(c_file)
      write (10,'(a,i3)')   'Number of channels: ',n_bins
      write (10,'(a,2f8.1)') 'LatN, LatS: ',latN,latS 
      write (10,'(a,2f8.1)') 'LonW, LonE: ',lonW,lonE 
      write (10,'(2a)') 'Surface type: ',c_surf
      print *,'Header written to output file'
!
! Write table of average counts, means, and stdvs
! count used is the number of channels actually used
! this considers only channels averaging 1 or more per time slot.
      write (10,'(a)') ' '
      write (10,'(3a)') 'Table 1: in, plevs, index_used, avg counts, ', &
                        'means(omf,oma,amf), stdvs(omf,oma,amf), ', &
                        'r=stdv(oma)/stdv(omf)'
      count_sum=0   ! count # of obs used
      count_used=0  ! count # of channels used
      do ib=1,n_bins
        count_sum=count_sum+iaccept(ib)
        if (iaccept(ib) >0) then
          count_used=count_used+1  
          write (10,'(i5,f12.2,i5,i8,7f12.6)') ib,plevs(ib), & 
                count_used,iaccept(ib),  &
                (xmean(ib,n),n=1,3),(xstdv(ib,n),n=1,3),aoverb(ib)
        else
          write (10,'(i5,f12.2,i5,i8,7f12.6)') ib,plevs(ib), &
                iaccept(ib),iaccept(ib), &
                (xmean(ib,n),n=1,3),(xstdv(ib,n),n=1,3),aoverb(ib)
        endif 
      enddo
      print *,'Means, stdvs written to output file'
!
! Write sample info about numbers of observations in each bin
      write (10,'(a)')  ' '
      write (10,'(2a)'), 'Table 2: Average number of obs in lat range ', &
                         'per synop period'
      write (10,'(a)')  ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      write (10,'(a)')  ' '
      dlat=180./(nlats-1)
      do n=1,nlats
        xlat=-90.+dlat*(n-1)
        do i=1,i_list_last
          i_values(i)=ibins(n,i_list(i))/n_count
        enddo 
        write (10,'(i4,f10.3,15i7)') n,xlat,i_values(1:i_list_last)
      enddo
      print *,'Bin counts   written to output file'
!
      write (10,'(a)') ' '
      write (10,'(2a)') 'last data file ',trim(c_file)
      write (10,'(a,i8)') 'sum of counts = ',count_sum
      close (10)
!
      print *,'Program completed'
!
      end program countobs_tables
