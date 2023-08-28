      program hcorr_tables
!
! Determine tables of horizontal correlations and other statistics 
! from file of previously accumulated sums
!
! At the end of the program, for each level or channel, icount is the average 
! obs count per assimilation period and xbins is the average number of 
! independent pairs of obs in each separation distance interval bin per 
! assimilation period. 
!
      implicit none
!
      integer, parameter :: rkind=8
      integer, parameter :: ikind=8
      integer(ikind), allocatable :: iaccept(:)
      integer(ikind), allocatable :: ibins(:,:)
      integer :: ib
      integer :: iks,     iks2
      integer :: i
      integer :: i_list(15)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(15)
      integer :: kt_look, kx_look
      integer :: n
      integer :: n_bins
      integer :: n_count
      integer :: n_ksmax
      integer :: n_date
      integer :: n_time
      integer :: nlevs
      integer :: argc
!
      real(rkind), parameter :: zero=0._rkind
      real(rkind) :: x
      real(rkind) :: xm1, xm2, xs1, xs2, xa, xb
      real(rkind) :: latN, latS, lonW, lonE
      real(rkind) :: delp
      real(rkind) :: x_max
      real(rkind) :: x1  
      real(rkind) :: xfac ! change from 1/n for mean to 1/(n-1) for var
      real(rkind) :: xbin
      real(rkind) :: x_values(15)
      real(rkind), allocatable :: xbins(:,:)
      real(rkind), allocatable :: cbins(:,:,:)  ! sums of products read
      real(rkind), allocatable :: corbins(:,:)  ! correlations
      real(rkind), allocatable :: covbins(:,:)  ! covariances
      real(rkind), allocatable :: xsum(:)
      real(rkind), allocatable :: xsumsq(:)
      real(rkind), allocatable :: xmean(:)
      real(rkind), allocatable :: xmsqr(:)
      real(rkind), allocatable :: xstdv(:)
      real(rkind), allocatable :: plevs(:)
!
      character(len=3)   :: c_comp
      character(len=220) :: c_file
      character(len=220) :: file_in
      character(len=220) :: file_out
!
      argc = iargc()
      if (argc .ne. 2) then
        print *,' usage must be: hcorr_tables.x filein fileout'
        stop
      endif
      call GetArg( 1_4, file_in)
      call GetArg( 2_4, file_out)
!
! read file of previously accumulated values
      open (unit=10,file=trim(file_in),form='unformatted')
      read (10) n_count,n_ksmax,n_date,n_time,n_bins,x_max, &
                x1,xbin,latN,latS,lonW,lonE,                &
                delp,kt_look,kx_look,nlevs,c_comp,c_file
      allocate (xsum(n_ksmax))
      allocate (xsumsq(n_ksmax))
      allocate (iaccept(n_ksmax))
      allocate (plevs(nlevs))      
      allocate (ibins(n_bins,n_ksmax))
      allocate (cbins(n_bins,5,n_ksmax))
      allocate (corbins(n_bins,n_ksmax))
      allocate (covbins(n_bins,n_ksmax))
      read (10) plevs
      read (10) iaccept,xsum,xsumsq
      read (10) ibins 
      read (10) cbins     
      close (10)     
      print *,'Input file read: ',trim(file_in)
!
      allocate (xmean(n_ksmax))
      allocate (xmsqr(n_ksmax))
      allocate (xstdv(n_ksmax))
      allocate (xbins(n_bins,n_ksmax))
      xmean(:)=zero
      xstdv(:)=zero
      xbins(:,:)=zero
!
      do iks=1,n_ksmax
        if (iaccept(iks) > 0) then
          xmean(iks)=xsum(iks)/iaccept(iks)
          xmsqr(iks)=xsumsq(iks)/iaccept(iks) - xmean(iks)*xmean(iks)
          if (xmsqr(iks) > zero) then
            xstdv(iks)=sqrt(xmsqr(iks))
          endif
          do n=1,n_bins
            xbins(n,iks)=real(ibins(n,iks))/real(n_count)
            if (ibins(n,iks) > 0) then
              xm1=cbins(n,2,iks)/ibins(n,iks)
              xm2=cbins(n,4,iks)/ibins(n,iks)
              xs1=cbins(n,3,iks)/ibins(n,iks)
              xs2=cbins(n,5,iks)/ibins(n,iks)
              xa=(xs1-xm1*xm1)*(xs2-xm2*xm2) 
              if (ibins(n,iks) > 1) then
                xfac=real(ibins(n,iks))/real(ibins(n,iks)-1)
              else  
                xfac=0.
              endif
              if (xa > zero) then   
                xb=cbins(n,1,iks)/ibins(n,iks)-xm1*xm2
                corbins(n,iks)=xb/sqrt(xa)  ! correlation
                covbins(n,iks)=xfac*xb      ! covariance
              else
                corbins(n,iks)=zero
                covbins(n,iks)=zero
              endif
            else
              corbins(n,iks)=zero
              covbins(n,iks)=zero
            endif   
          enddo
        endif
      enddo 
!
! Replace total number of obs by average number per synoptic period
      iaccept(:)=iaccept(:)/n_count
!
! Print table of average counts, means, and variances
      print ('(a,4i10)'),'Header info: ',n_count,n_ksmax,n_date,n_time
      print *,' '
      print *,'Table 1: mean count per synoptic time, means, stdvs'
      do iks=1,n_ksmax
        print ('(2i8,2f12.6)'),iks,iaccept(iks),xmean(iks),xstdv(iks)
      enddo
!
! Determine parameters for printing a subset of distributions
      i_list_skip=1+(n_ksmax-1)/15
      i=0
      do iks=1,n_ksmax,i_list_skip
        if (i < 16) then
          i=i+1
          i_list(i)=iks
        endif 
      enddo
      i_list_last=i
!
! Open output file 
      open (unit=10,file=trim(file_out),form='formatted')
      print *,' '
      print *,'Output file opened: ',trim(file_out)
!
! Write header information
      write (10,'(2a)') 'Horizontal Correlations and Covariances of ',c_comp
      write (10,'(a,i3,a)') 'Data from ',n_count,' times'     
      write (10,'(a,2i12)') 'Last date and time: ',n_date,n_time
      write (10,'(a,a)')    'Last file of data: ',trim(c_file)
      write (10,'(a,i3)')   'Number of bins: ',n_bins
      write (10,'(a,i3)')   'Number of channels or p-levels: ',n_ksmax
      write (10,'(a,f8.1)') 'Bin width in kilometers: ',xbin 
      write (10,'(a,2f8.1)') 'LatN, LatS=',latN,latS 
      write (10,'(a,2f8.1)') 'LonW, LonE=',lonW,lonE 
      if ( (kx_look > 99) .and. (kx_look < 300 ) ) then ! conv obs
        write (10,'(a,f8.1,2i5)') 'delp, kt, kx=',delp,kt_look,kx_look
      endif 
      print *,'Header written to output file'
!
! Write average counts, means and standard deviations of O-F
      write (10,'(a)') ' '
      write (10,'(a)') 'Table 1: mean count per synoptic time, means, stdvs'
      do iks=1,n_ksmax
        write (10,'(2i8,2f12.6)') iks,iaccept(iks),xmean(iks),xstdv(iks)
      enddo
      print *,'Means and standard deviations written to output file'
!
! Write sample info about numbers of observations in each bin
      write (10,'(a)') ' '
      write (10,'(a)'),'Table 2: Average bin counts per synoptic period ', &
                       '(divided by iaccept)' 
      write (10,'(a)') ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i))
      enddo 
      write (10,'(6x,a,15i7)') 'iaccept=',i_values(1:i_list_last)
      write (10,'(a)') ' '
      do n=1,n_bins
        do i=1,i_list_last
          iks=i_list(i)
          x_values(i)=xbins(n,i_list(i))
        enddo 
        if (n == 1) then
          x=0.
        else  
          x=(n-1.5)*xbin
        endif
        write (10,'(i4,f10.4,15f7.0)') n,x,x_values(1:i_list_last)
      enddo
      print *,'Average bin counts written to output file'
!
! Write a subset of correlations 
      write (10,'(a)') ' '
      write (10,'(2a)') 'Table 3: Horizontal Correlations of ',c_comp
      write (10,'(a)') ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i))
      enddo 
      write (10,'(6x,a,15i7)') 'iaccept=',i_values(1:i_list_last)
      write (10,'(a)') ' '
      do n=1,n_bins
        do i=1,i_list_last
          iks=i_list(i)
          x_values(i)=corbins(n,i_list(i))
        enddo 
        if (n == 1) then
          x=0.
        else  
          x=(n-1.5)*xbin
        endif
        write (10,'(i4,f10.4,15f7.3)') n,x,x_values(1:i_list_last)
      enddo
      print *,'Correlations written to output file'
!
! Write a subset of covariances
      write (10,'(a)') ' '
      write (10,'(2a)') 'Table 4: Horizontal Covariances of ',c_comp
      write (10,'(a)') ' '
      write (10,'(6x,a,15i7)') 'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i))
      enddo 
      write (10,'(6x,a,15i7)') 'iaccept=',i_values(1:i_list_last)
      write (10,'(a)') ' '
      do n=1,n_bins
        do i=1,i_list_last
          iks=i_list(i)
          x_values(i)=covbins(n,i_list(i))
        enddo 
        if (n == 1) then
          x=0.
        else  
          x=(n-1.5)*xbin
        endif
        write (10,'(i4,f10.4,15f7.3)') n,x,x_values(1:i_list_last)
      enddo
      print *,'Covariances written to output file'
!
      write (10,'(a)') ' '
      write (10,'(2a)') 'last data file ',trim(c_file)
      close (10)
!
! Write binary file
      open (unit=30,file='corr.bin',form='unformatted')
      write (30) n_count,n_ksmax,n_date,n_time,n_bins,x_max, &
                x1,xbin,latN,latS,lonW,lonE,                 &
                delp,kt_look,kx_look,nlevs,c_file   
      write (30) iaccept,xmean,xstdv,plevs
      write (30) xbins
      write (30) corbins
      write (30) covbins
      close (30)     
!
      print *,'Program completed'
!
      end program hcorr_tables
