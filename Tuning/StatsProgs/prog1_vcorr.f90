      program compute_stats
!
      implicit none
!
      integer, parameter :: rkind=8
      integer, parameter :: ikind=8
      integer, parameter :: n_p_select=14
      integer(ikind), allocatable :: iaccept(:,:)
      integer(ikind), allocatable :: ibins(:,:,:)
      integer :: ib, id
      integer :: i
      integer :: i_list(n_p_select)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(n_p_select)
      integer :: n
      integer :: n_bins
      integer :: n_count
      integer :: n_date
      integer :: n_time
      integer :: argc
!
      real(rkind), parameter :: zero=0._rkind
      real(rkind) :: p
      real(rkind) :: xa, xb, xm1, xm2, xs1, xs2
      real(rkind) :: xlog1, xlog2, xbin
      real(rkind) :: x_values(n_p_select)
      real(rkind) :: p_list(n_p_select)
      real(rkind) :: p_select(n_p_select)
      real(rkind), allocatable :: cbins(:,:,:,:)  ! sums read
      real(rkind), allocatable :: corbins(:,:,:)  ! correlations
      real(rkind), allocatable :: covbins(:,:,:)  ! covariances
      real(rkind), allocatable :: xsum(:,:)
      real(rkind), allocatable :: xsumsq(:,:)
      real(rkind), allocatable :: xmean(:,:)
      real(rkind), allocatable :: xmsqr(:,:)
      real(rkind), allocatable :: xstdv(:,:)
!
      character(len=4), parameter :: f_names='Tquv'
      character(len=220) :: c_file
      character(len=220) :: file_in
      character(len=220) :: file_out
!
      data p_select/10., 20., 50., 100., 200., 300., 400., 500., 600., 700., &
                    800., 850., 925., 1000./ 
!
      argc = iargc()
      if (argc .ne. 2) then
        print *,' usage must be: prog1 filein fileout'
        stop
      endif
      call GetArg( 1_4, file_in)
      call GetArg( 2_4, file_out)
!
! read file of previously accumulated values
      open (unit=10,file=trim(file_in),form='unformatted')
      read (10) n_count,n_date,n_time,n_bins,xlog1,xlog2,xbin,c_file
      allocate (xsum(n_bins,4))
      allocate (xsumsq(n_bins,4))
      allocate (iaccept(n_bins,4))
      allocate (ibins(n_bins,n_bins,4))
      allocate (cbins(n_bins,n_bins,5,4))
      allocate (corbins(n_bins,n_bins,4))
      allocate (covbins(n_bins,n_bins,4))
      read (10) iaccept,xsum,xsumsq
      read (10) ibins 
      read (10) cbins     
      close (10)     
      print *,'Input file read: ',trim(file_in)
!
      allocate (xmean(n_bins,4))
      allocate (xmsqr(n_bins,4))
      allocate (xstdv(n_bins,4))
      xmean(:,:)=zero
      xstdv(:,:)=zero
!
! Loop over field types 
      do id=1,4 
!
! Compute means and standard deviations
        do ib=1,n_bins
          if (iaccept(ib,id) > 0) then
            xmean(ib,id)=xsum(ib,id)/iaccept(ib,id)
            xmsqr(ib,id)=xsumsq(ib,id)/iaccept(ib,id) - &
                         xmean(ib,id)*xmean(ib,id)
            if (xmsqr(ib,id) > zero) then
              xstdv(ib,id)=sqrt(xmsqr(ib,id))
            endif
          endif
        enddo
!
! Compute correlations and covariances
        do ib=1,n_bins
          if (iaccept(ib,id) > 0) then
            do n=1,n_bins
              if (ibins(n,ib,id) > 0) then
                xm1=cbins(ib,n,2,id)/ibins(ib,n,id)
                xm2=cbins(ib,n,4,id)/ibins(ib,n,id)
                xs1=cbins(ib,n,3,id)/ibins(ib,n,id)
                xs2=cbins(ib,n,5,id)/ibins(ib,n,id)
                xa=(xs1-xm1*xm1)*(xs2-xm2*xm2)
                if (xa > zero) then   
                  xb=cbins(n,ib,1,id)/ibins(n,ib,id)-xm1*xm2
                  corbins(n,ib,id)=xb/sqrt(xa)   ! correlation
                  covbins(n,ib,id)=xb            ! covariance
                else
                  corbins(n,ib,id)=zero
                  covbins(n,ib,id)=zero
                endif
              else
                corbins(n,ib,id)=zero
                covbins(n,ib,id)=zero
              endif   
            enddo
          endif
        enddo 
!
      enddo ! loop over field types
!
! Replace total number of obs by average number per synoptic period
      iaccept(:,:)=iaccept(:,:)/n_count
      print ('(a,4i10)'),'Header info: ',n_count,n_bins,n_date,n_time
!
! Determine parameters for printing a subset of distributions
      i=1
      do ib=1,n_bins
        p=exp(xlog1+xbin*(ib-1))
        if ( (p > p_select(i)+0.1) .and. (i < n_p_select) ) then
          i=i+1  
        endif
        i_list(i)=ib
        p_list(i)=p
      enddo
      i_list_last=i
!
! Open output file 
      open (unit=10,file=trim(file_out),form='formatted')
      print *,' '
      print *,'Output file opened: ',trim(file_out)
!
! Write header information
      write (10,'(a)') 'Vertical Correlations and Covariances of O-F'
      write (10,'(a,i3,a)') 'Data from ',n_count,' times'     
      write (10,'(a,2i12)') 'Last date and time: ',n_date,n_time
      write (10,'(a,a)')    'Last file of data: ',trim(c_file)
      write (10,'(a,i3)')   'Number of bins: ',n_bins
      write (10,'(a,f8.1)') 'Bin width in log(p): ',xbin 
      print *,'Header written to output file'
!
! Loop over field types
      do id=1,4
        print *,'Writing results for field=',f_names(id:id)
!
! Write table of average counts, means, and variances
        write (10,'(a)') ' '
        write (10,'(2a)') 'Table 1: avg counts, means, stdvs for field=', &
                           f_names(id:id)
        do ib=1,n_bins
          p=exp(xlog1+xbin*(ib-1))
          write (10,'(i8,f6.0,i8,2f12.6)') ib,p,iaccept(ib,id), &
                                           xmean(ib,id),xstdv(ib,id)
        enddo
        print *,'Means, stdvs written to output file for id=',id 
!
! Write sample info about numbers of observations in each bin
        write (10,'(a)') ' '
        write (10,'(3a)'),'Table 2: Average bin counts per synoptic ', &
                          'period for field=',f_names(id:id)
        write (10,'(a)') ' '
        write (10,'(6x,a,15i7)') 'bin=    ',i_list(1:i_list_last)
        do i=1,i_list_last
          x_values(i)=exp(xlog1+xbin*(i_list(i)-1))
        enddo 
        write (10,'(6x,a,15f7.0)') 'plev=   ',x_values(1:i_list_last)
        do i=1,i_list_last
          i_values(i)=iaccept(i_list(i),id)/n_count
        enddo 
        write (10,'(6x,a,15i7)') 'iaccept=',i_values(1:i_list_last)
        write (10,'(a)') ' '
        do n=1,n_bins
          do i=1,i_list_last
            i_values(i)=ibins(n,i_list(i),id)
          enddo 
          p=exp(xlog1+xbin*(n-1))
          write (10,'(i4,f10.0,15i7)') n,p,i_values(1:i_list_last)
        enddo
        print *,'Bin counts   written to output file for id=',id
!
! Write a subset of correlations 
        write (10,'(a)') ' '
        write (10,'(2a)') &
             'Table 3: Vertical Correlations of O-F for field=',f_names(id:id)
        write (10,'(a)') ' '
        write (10,'(6x,a,15i7)') 'bin=    ',i_list(1:i_list_last)
        do i=1,i_list_last
          x_values(i)=exp(xlog1+xbin*(i_list(i)-1))
        enddo 
        write (10,'(6x,a,15f7.0)') 'plev=   ',x_values(1:i_list_last)
        write (10,'(a)') ' '
        do n=1,n_bins
          do i=1,i_list_last
            x_values(i)=corbins(n,i_list(i),id)
          enddo 
          p=exp(xlog1+xbin*(n-1))
          write (10,'(i4,f10.0,15f7.3)') n,p,x_values(1:i_list_last)
        enddo
        print *,'Correlations written to output file for id=',id
!
! Write a subset of covariances
        write (10,'(a)') ' '
        write (10,'(2a)') &
                'Table 4: Vertical Covariances of O-F for field=', &
                f_names(id:id)
        write (10,'(a)') ' '
        write (10,'(6x,a,15i7)') 'bin=    ',i_list(1:i_list_last)
        do i=1,i_list_last
          x_values(i)=exp(xlog1+xbin*(i_list(i)-1))
        enddo 
        write (10,'(6x,a,15f7.0)') 'plev=   ',x_values(1:i_list_last)
        write (10,'(a)') ' '
        do n=1,n_bins
          do i=1,i_list_last
            x_values(i)=covbins(n,i_list(i),id)
          enddo 
          p=exp(xlog1+xbin*(n-1))
          write (10,'(i4,f10.0,15f7.3)') n,p,x_values(1:i_list_last)
        enddo
        print *,'Covariances  written to output file for id=',id
!
! End loop over field types
      enddo 
!
      write (10,'(a)') ' '
      write (10,'(2a)') 'last data file ',trim(c_file)
      close (10)
!
      print *,'Program completed'
!
      end program compute_stats
