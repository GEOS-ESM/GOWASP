      program chcorr_tables
!
! Determine tables of channel correlations and other statistics 
! from file of previously accumulated sums.  
! Print a sample of correlations and covariances
! Write full covariance and correlations to a file
!
! At the end of the program, icount is the average obs count for each channel 
! per assimilation period and xbins is the average number of independent 
! pairs of obs in each channel pair bin per assimilation period 
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
      integer :: argc
!
      real(rkind), parameter :: zero=0._rkind
      real(rkind) :: latN, latS, lonW, lonE
      real(rkind) :: xfac ! change from 1/n for mean to 1/(n-1) for var
      real(rkind) :: xa, xb, xm1, xm2, xs1, xs2
      real(rkind) :: x_values(n_c_select)
      real(rkind), allocatable :: cbins(:,:,:)  ! sums read
      real(rkind), allocatable :: corbins(:,:)  ! correlations
      real(rkind), allocatable :: covbins(:,:)  ! covariances
      real(r4),    allocatable :: cor_used(:,:) ! correlations of only ch used
      real(r4),    allocatable :: cov_used(:,:) ! covariances  of only ch used
      real(rkind), allocatable :: xbins(:,:)    ! counts in bins
      real(rkind), allocatable :: xsum(:)
      real(rkind), allocatable :: xsumsq(:)
      real(rkind), allocatable :: xmean(:)
      real(rkind), allocatable :: xmsqr(:)
      real(rkind), allocatable :: xstdv(:)
!
      character(len=3)   :: c_comp
      character(len=220) :: c_file
      character(len=220) :: file_in
      character(len=220) :: file_out
      character(len=224) :: file_out_1
      character(len=224) :: file_out_2
!
! Specify a subset of channels for which tables will be printed.  
! The following is for AIRS, otherwise this selection will be overwritten
       data c_select/96,97,98,99,100,101,102,103,104,176,177,178,181,182,183/ 
!
      argc = iargc()
      if (argc .ne. 2) then
        print *,' usage must be: vcorr_tables.x filein fileout'
        stop
      endif
      call GetArg( 1_4, file_in)
      call GetArg( 2_4, file_out)
      file_out_1=trim(file_out)//'.txt'
      file_out_2=trim(file_out)//'.bin'
!
! read file of previously accumulated values
      open (unit=10,file=trim(file_in),form='unformatted')
      read (10) n_count,n_bins,n_date,n_time, &
                latN,latS,lonW,lonE,c_comp,c_file
      allocate (xsum(n_bins))
      allocate (xsumsq(n_bins))
      allocate (iaccept(n_bins))
      allocate (ibins(n_bins,n_bins))
      allocate (cbins(n_bins,n_bins,5))
      allocate (corbins(n_bins,n_bins))
      allocate (covbins(n_bins,n_bins))
      allocate (xbins(n_bins,n_bins))
      read (10) iaccept,xsum,xsumsq
      read (10) ibins 
      read (10) cbins     
      close (10)     
      print *,'Input file read: ',trim(file_in)
!
      allocate (xmean(n_bins))
      allocate (xmsqr(n_bins))
      allocate (xstdv(n_bins))
      xmean(:)=zero
      xstdv(:)=zero
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

          xmean(ib)=xsum(ib)/iaccept(ib)
          xmsqr(ib)=xsumsq(ib)/iaccept(ib) - &
                       xmean(ib)*xmean(ib)
          if (xmsqr(ib) > zero) then
            xstdv(ib)=sqrt(xmsqr(ib))
          endif
        endif
      enddo
!
! Compute correlations and covariances
      do ib=1,n_bins
        if (iaccept(ib) > 0) then
          do n=1,n_bins
            if (ibins(n,ib) > 0) then
              xm1=cbins(ib,n,2)/ibins(ib,n)
              xm2=cbins(ib,n,4)/ibins(ib,n)
              xs1=cbins(ib,n,3)/ibins(ib,n)
              xs2=cbins(ib,n,5)/ibins(ib,n)
              xa=(xs1-xm1*xm1)*(xs2-xm2*xm2)
              if (ibins(n,ib) > 1) then
                xfac=real(ibins(n,ib))/real(ibins(n,ib)-1)
              else  
                xfac=0.
              endif
              if (xa > zero) then   
                xb=cbins(n,ib,1)/ibins(n,ib)-xm1*xm2
                corbins(n,ib)=xb/sqrt(xa)   ! correlation
                covbins(n,ib)=xfac*xb       ! covariance
              else
                corbins(n,ib)=zero
                covbins(n,ib)=zero
              endif
            else
              corbins(n,ib)=zero
              covbins(n,ib)=zero
            endif   
          enddo
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
      write (10,'(2a)') 'Channel Correlations and Covariances of ',c_comp
      write (10,'(a,i3,a)') 'Data from ',n_count,' times'     
      write (10,'(a,2i12)') 'Last date and time: ',n_date,n_time
      write (10,'(a,a)')    'Last file of data: ',trim(c_file)
      write (10,'(a,i3)')   'Number of channels: ',n_bins
      write (10,'(a,2f8.1)') 'LatN, LatS: ',latN,latS 
      write (10,'(a,2f8.1)') 'LonW, LonE: ',lonW,lonE 
      print *,'Header written to output file'
!
! Write table of average counts, means, and stdvs
! count used is the number of channels actually used
! this considers only channels averaging 1 or more per time slot.
      write (10,'(a)') ' '
      write (10,'(a)') 'Table 1: index_used, avg counts, means, stdvs'
      count_used=0
      do ib=1,n_bins
        if (iaccept(ib) >0) then
          count_used=count_used+1  
          write (10,'(2i6,i8,2f12.6)') ib,count_used,iaccept(ib), &
                                       xmean(ib),xstdv(ib)
        else
          write (10,'(i6,6x,i8,2f12.6)') ib,iaccept(ib),   &
                                       xmean(ib),xstdv(ib)
        endif 
      enddo
      print *,'Means, stdvs written to output file'
!
! Write sample info about numbers of observations in each bin
      write (10,'(a)')  ' '
      write (10,'(2a)'), 'Table 2: Average bin counts per '&
                         'synoptic period (in hundreds)'
      write (10,'(a)')  ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      write (10,'(a)')  ' '
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i))/n_count
      enddo 
      write (10,'(6x,a,15i7)') 'iaccept=',i_values(1:i_list_last)
      write (10,'(a)') ' '
      do n=1,n_bins
        do i=1,i_list_last
          i_values(i)=ibins(n,i_list(i))/100
        enddo 
        write (10,'(i4,10x,15i7)') n,i_values(1:i_list_last)
      enddo
      print *,'Bin counts   written to output file'
!
! Write a subset of correlations 
      write (10,'(a)')  ' '
      write (10,'(2a)') 'Table 3: Channel Correlations of ',c_comp
      write (10,'(a)')  ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      write (10,'(a)')  ' '
      do n=1,n_bins
        do i=1,i_list_last
          x_values(i)=corbins(n,i_list(i))
        enddo 
        write (10,'(i4,10x,15f7.3)') n,x_values(1:i_list_last)
      enddo
      print *,'Correlations written to output file'
!
! Write a subset of covariances
      write (10,'(a)')  ' '
      write (10,'(2a)') 'Table 4: Channel Covariances of ',c_comp
      write (10,'(a)')  ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      write (10,'(a)')  ' '
      do n=1,n_bins
        do i=1,i_list_last
          x_values(i)=covbins(n,i_list(i))
        enddo 
        write (10,'(i4,10x,15f7.3)') n,x_values(1:i_list_last)
      enddo
      print *,'Covariances  written to output file'
!
      write (10,'(a)') ' '
      write (10,'(2a)') 'last data file ',trim(c_file)
      close (10)
!
! Fill cor and cov matrices included only channels whouse counts average
! more than 1 per time slot 
!
      allocate (cor_used(count_used,count_used))
      allocate (cov_used(count_used,count_used))
!
! Extract array elements only for channels used
      i1=0
      do ib=1,n_bins
        if (iaccept(ib) > 0) then
          i1=i1+1
          i2=0
          do n=1,n_bins
            if (iaccept(n) > 0) then
              i2=i2+1
              cor_used(i1,i2)=corbins(ib,n)          
              cov_used(i1,i2)=covbins(ib,n)          
            endif
          enddo
        endif
      enddo
!
!
! Write binary file
      do i=1,n_bins 
        do n=1,n_bins 
          xbins(i,n)=real(ibins(i,n))/real(n_count)
        enddo
      enddo
!
      open (unit=30,file=trim(file_out_2),form='unformatted')
      write (30) n_count,n_date,n_time,n_bins,count_used, &
                 latN,latS,lonW,lonE   
      write (30) iaccept,xmean,xstdv
      write (30) xbins
      write (30) corbins
      write (30) covbins
      write (30) cor_used
      write (30) cov_used
      close (30) 
      print *,'Binary File written'
!
      deallocate (cor_used, cov_used)
      deallocate (xsum,xsumsq,iaccept,ibins,cbins,corbins,covbins)
      deallocate (xmean,xmsqr,xstdv,xbins)
!
      print *,'Program completed'
!
      end program chcorr_tables
