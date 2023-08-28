   module m_satloc
!
   use m_kinds, only : rkind1
   implicit none
!
   private
   public satloc_filter
   public satloc_setup
   public satloc_bins_get
   public satloc_bins_write
   public satloc_bins_print
   public satloc_clean
!
   integer, parameter, public :: satloc_nbins=36.
   integer, allocatable, public :: satloc_bins(:,:,:)
!
   real(rkind1), parameter, public :: satloc_geost_view=60.
   real(rkind1), allocatable, public :: satloc_dlon(:,:)  ! lon1, dlon, lon2
!
   logical :: satloc_file_binary
   integer, parameter :: iunit=33
   integer :: satloc_ntimes
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_clean
   deallocate (satloc_bins,satloc_dlon)
   end subroutine satloc_clean
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_setup (kx_num,kx_type,kx_locs)  
!
   implicit none
   integer, intent(in) :: kx_num
   real(rkind1), intent(in) :: kx_locs(kx_num)
   character(len=*), intent(in) :: kx_type(kx_num)
!
   integer :: k
   real(rkind1) :: range
   real(rkind1) :: lon_min
   real(rkind1) :: dlon
!
   allocate (satloc_dlon(3,kx_num))     
   allocate (satloc_bins(satloc_nbins,2,kx_num))
!     
   do k=1,kx_num
     if (kx_type(k)(1:1) == 'G') then ! geostationary   
       range=2.*satloc_geost_view
       lon_min=kx_locs(k)-satloc_geost_view
     else
       range=360.
       lon_min=0.
     endif
!
     dlon=range/real(satloc_nbins)
     if (lon_min < 0.)then
       lon_min= 360.+lon_min
     endif
!
     satloc_dlon(1,k)=lon_min
     satloc_dlon(2,k)=dlon
     satloc_dlon(3,k)=mod(lon_min+dlon*satloc_nbins,360.)
!
   enddo
!   
   satloc_ntimes=1  ! default (or first) value
!
   end subroutine satloc_setup    
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_skip_times (ntime,kx_num,dhours,cdtime0_prog, &
                                 cdtime0_file,lprint)
!
! Skip to a file time based on a modular=arithmic claculation:
! Essentially, the times written on the file comprise a time span that 
! is then referenced as though it is periodic throughout the calendar year. 
!
! It is assumed that the spacing between times in the file is the same as the 
! time-spacing between interpolated fields in the program, which would always 
! be the case if the same kx_table_in file is used to create and use the 
! satloc file.
!
   use m_time_compute, only : time_compute_unpack
   use m_time_compute, only : time_compute_dhours
   use m_time_compute, only : rkindh
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(in) :: kx_num
   integer, intent(inout) :: ntime 
   real(rkindh), intent(in) :: dhours   ! t-spacing of interpolated fields
   character(len=*), intent(in) :: cdtime0_prog
   character(len=*), intent(in) :: cdtime0_file
!
   integer :: n,jx
   integer :: nskip
   integer :: ierr
   real(rkind1) :: time_prog(6),time_file(6)
   real(rkindh) :: diff_hours 
   character(len=14) :: cdum
!
! Compute difference between 1st data time in program and 1st time on file
! First change from year on ods file to NR year
   cdum=cdtime0_file
   cdum(1:4)=cdtime0_prog(1:4)  
   call time_compute_unpack (cdum,time_file)
   call time_compute_unpack (cdtime0_prog,time_prog)
   call time_compute_dhours (time_file,time_prog,diff_hours,ierr)
!
   if (diff_hours >= 0._rkindh) then 
     nskip=int(1.e-8+diff_hours/dhours)
     nskip=mod(nskip,satloc_ntimes)
   else
     nskip=int(-1.e-8+diff_hours/dhours)
     nskip=mod(nskip,satloc_ntimes)
     if (nskip < 0) then
       nskip=satloc_ntimes+nskip
     endif 
   endif
   jx=kx_num*(1+2*(1+(satloc_nbins-1)/12)) ! 12=num of values 1 file record  
!
   do n=1,nskip*jx
     if (satloc_file_binary) then 
       read (iunit) 
     else
       read (iunit,'(a1)') cdum(1:1)
     endif
   enddo
   ntime=ntime+nskip
!
   if (lprint) then
     print ('(a,i8,f12.1)'),'times skipped in satloc-file: ', &
           nskip,diff_hours
   endif
!
   end subroutine satloc_skip_times  
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_bins_get (kx_num,kx_type,kx_satid,ntime,thours, &
                               dhours,cdtime0,filename,lprint,ierr)  
!
! Get counts for desired time from satloc file.
!
   use m_time_compute, only : rkindh
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(inout) :: ntime
   integer, intent(in) :: kx_num
   integer, intent(out) :: ierr
   integer, intent(in) :: kx_satid(kx_num)
   real(rkindh), intent(in) :: thours   ! time since beginning of program
   real(rkindh), intent(in) :: dhours   ! t-spacing of interpolated fields
   character(len=*), intent(in) :: kx_type(kx_num)
   character(len=*), intent(in) :: cdtime0
   character(len=*), intent(in) :: filename
!
   integer :: len_fname
   integer :: ix1,ix2,nx
   integer :: nt,kx1,kx2,kx3
   integer :: k,nlat
   integer :: j,jx,j1,j2
   character(len=17) :: cdum
   character(len=14) :: cdtime
!
   ierr=0  ! default indicates OK
!
   if (trim(filename) == 'none') then  ! simulate locations
     call satloc_simulate (kx_num,kx_type,kx_satid,thours)  
!
   else        
! read locations
! 1st, if this is the first time to be processed or if the first time on the file 
! should be used, then the satloc file must be opened or re-opened and its header
! read. When program starts, satloc_ntimes is set to 1; thereafter it is read from
! the header in the satloc file. 
     if (satloc_ntimes == 1 .or. mod(ntime,satloc_ntimes) == 1) then 
!
! Determine if formated or binary file (the latter is much smaller for same data)  
       len_fname=len_trim(filename) 
       satloc_file_binary= filename(len_fname-2:len_fname) == 'bin'
!
! Open file and read header for desired format      
       if (satloc_file_binary) then 
         open (iunit,file=trim(filename),form='unformatted')
         read (iunit) ix1,ix2,satloc_ntimes,cdtime 
       else
         open (iunit,file=trim(filename),form='formatted')
         read (iunit,'(3i4,2x,a)') ix1,ix2,satloc_ntimes,cdtime 
       endif
!
       if (lprint) then 
         print *,'Satloc file opened = ',trim(filename) 
       endif
       ierr=0 ! default value
!
! Check if some times on file need to be skipped.  If so, skip them and reset ntime
! to the cosseponding time index in the satloc file.
       if (ntime == 1) then 
         if (ix1 == satloc_nbins .and. ix2 == kx_num) then
           call satloc_skip_times (ntime,kx_num,dhours,cdtime0, &
                                   cdtime,lprint)
         else 
           ierr=1
           if (lprint) then
             print *,' '
             print *,'Incompatible satloc_bins file:'
             print *,'ix1,ix2,satloc_nbins,kx_num=',ix1,ix2, &
                   satloc_nbins,kx_num
           endif
         endif
       endif   ! test on ntime 
!
     endif     ! test on whether satloc file need to be opened.
!     
     if (ierr == 0) then
       jx=1+(satloc_nbins-1)/12  ! 12=num of values in 1 file record  
       do k=1,kx_num
         if (satloc_file_binary) then 
           read (iunit) nt,kx1,kx2,kx3,cdtime
         else
           read (iunit,'(a17,4i6,2x,a)') cdum(1:17),nt,kx1,kx2,kx3,cdtime
         endif
         do nlat=1,2   ! loop over lat 2 hemispheres S/N
           j2=0
           do j=1,jx   ! loop over records to read
             j1=j2+1
             j2=min(j1+11,satloc_nbins)
             if (j == 1) then 
               if (satloc_file_binary) then 
                 read (iunit) nx,satloc_bins(j1:j2,nlat,k)
               else
                 read (iunit,'(a5,i2,12i6)') cdum(1:5),nx, &
                                             satloc_bins(j1:j2,nlat,k)
               endif
             else  ! format differs
               if (satloc_file_binary) then 
                 read (iunit) satloc_bins(j1:j2,nlat,k)
               else
                 read (iunit,'(7x,12i6)') satloc_bins(j1:j2,nlat,k)
               endif
             endif 
           enddo
         enddo     ! loop over lat (hemispheres)
       enddo       ! loop over k 
       if (lprint) then 
         print *,' '
         print ('(a,i5,2x,a)'),'Satlocs read for ntime,cdtime=',ntime,cdtime
       endif
!
! If all times on satloc file have been read, close it 
       if (mod(ntime,satloc_ntimes) == 0) then 
         close (iunit)
       endif
       ntime=ntime+1
!  
     endif
   endif
!
   end subroutine satloc_bins_get 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_bins_write (ntime,ntimes,kx_num,cdtime, &
                                 filename,kx_list,kx_satid)  
!
   implicit none
!
   integer, intent(in) :: ntime
   integer, intent(in) :: ntimes
   integer, intent(in) :: kx_num
   integer, intent(in) :: kx_list(kx_num)
   integer, intent(in) :: kx_satid(kx_num)
   character(len=*), intent(in) :: cdtime
   character(len=*), intent(in) :: filename
!
   integer :: len_fname
   integer :: k,nlat
   integer :: j,jx,j1,j2
!
   if (ntime == 1) then
!
! Determine if formated or binary file (the latter is much smaller for same data)  
     len_fname=len_trim(filename) 
     satloc_file_binary= filename(len_fname-2:len_fname) == 'bin'
!
! Open file with desired format and write header
     if (satloc_file_binary) then 
       open (iunit,file=trim(filename),form='unformatted')
       write (iunit) satloc_nbins,kx_num,ntimes,cdtime 
     else
       open (iunit,file=trim(filename),form='formatted')
       write (iunit,'(3i4,2x,a)') satloc_nbins,kx_num,ntimes,cdtime 
     endif
   endif
!
   jx=1+(satloc_nbins-1)/12
   do k=1,kx_num
!
     if (satloc_file_binary) then 
       write (iunit) ntime,k,kx_list(k),kx_satid(k),cdtime 
     else
       write (iunit,'(a17,4i6,2x,a)') 'ntime,k,kx,satid=', &
                      ntime,k,kx_list(k),kx_satid(k),cdtime 
     endif
!
     do nlat=1,2
       j2=0
       do j=1,jx
         j1=j2+1
         j2=min(j1+11,satloc_nbins)
         if (j == 1) then 
           if (satloc_file_binary) then 
             write (iunit) nlat,satloc_bins(j1:j2,nlat,k)
           else
             write (iunit,'(a5,i2,12i6)') 'nlat=',nlat,satloc_bins(j1:j2,nlat,k)
           endif
         else
           if (satloc_file_binary) then 
             write (iunit) satloc_bins(j1:j2,nlat,k)
           else
             write (iunit,'(7x,12i6)') satloc_bins(j1:j2,nlat,k)
           endif
         endif 
       enddo
     enddo
   enddo 
!  
   if (ntime == ntimes) then
     close (iunit)
   endif
!
   print ('(a,i4,2a)'),'satloc_bins written for ntime=',ntime, &
                       '  cdtime=',trim(cdtime)
!
   end subroutine satloc_bins_write 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_bins_print (ntime,cdtime,kx_num,kx_list,kx_satid)  
!
   implicit none
!
   integer, intent(in) :: ntime
   integer, intent(in) :: kx_num
   integer, intent(in) :: kx_list(kx_num)
   integer, intent(in) :: kx_satid(kx_num)
   character(len=*), intent(in) :: cdtime
!
   integer :: k,nlat
   integer :: j,jx,j1,j2
!
   print *,' '
   print ('(a,i4,2a)'),' ntime=',ntime,'  cdtime=',trim(cdtime)
!
   jx=1+(satloc_nbins-1)/12
  
   do k=1,kx_num
     print ('(a,4i6)'),'ntime,k,kx,satid=',ntime,k,kx_list(k),kx_satid(k) 
     do nlat=1,2
       j2=0
       do j=1,jx
         j1=j2+1
         j2=min(j1+11,satloc_nbins)
         if (j == 1) then 
           print ('(a,i2,12i6)'),'nlat=',nlat,satloc_bins(j1:j2,nlat,k)
         else
           print ('(7x,12i6)'),satloc_bins(j1:j2,nlat,k)
         endif 
       enddo
     enddo
   enddo 
!  
   end subroutine satloc_bins_print 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_simulate (kx_num,kx_type,kx_satid,thours)  
!
   use m_time_compute, only : rkindh
   implicit none
!
   integer, intent(in) :: kx_num
   integer, intent(in) :: kx_satid(kx_num)
   real(rkindh), intent(in) :: thours
   character(len=*), intent(in) :: kx_type(kx_num)
!
   integer :: k
   integer :: icut
   integer :: i,ib
   integer :: ibin1,ispan,ilatshift
   real(rkind1), parameter :: himawari_cutlon=10.
   real(rkind1), parameter :: polar_rate=10.  ! degrees/hour
   real(rkind1), parameter :: polar_span1=140.
   real(rkind1), parameter :: polar_span2=80.
   real(rkind1), parameter :: polar_break=70.
   real(rkind1), parameter :: polar_latshift=60.
   real(rkind1) :: xtime
   real(rkind1) :: lon1
!
   do k=1,kx_num
     if (kx_type(k)(1:1) == 'G') then ! geostationary
       satloc_bins(:,:,k)=10000      ! default is to consider all lons
       if (kx_satid(k) == 173) then  ! for Himawari, cut edge lons
         icut=int(himawari_cutlon/satloc_dlon(2,k))
         if (icut > 0) then 
           satloc_bins(1:icut,:,k)=0
           icut=satloc_nbins-icut+1         
           satloc_bins(icut:satloc_nbins,:,k)=0
         endif         
       endif
!
     else                    ! polar
       satloc_bins(1:satloc_nbins,:,k)=0  ! default is no obs
       lon1=mod(thours*polar_rate,360.)
       ilatshift=int(polar_latshift/satloc_dlon(2,k))
!
       ibin1=1+int(lon1/satloc_dlon(2,k))
       ispan=1+int(polar_span1/satloc_dlon(2,k))
       do i=1,ispan
         ib=1+mod(ibin1+i-2,satloc_nbins)
         satloc_bins(ib,1,k)=10000       
         ib=1+mod(ib+ilatshift-1,satloc_nbins)
         satloc_bins(ib,2,k)=10000       
       enddo
!
       ibin1=ibin1+ispan+int(polar_break/satloc_dlon(2,k))
       ispan=1+int(polar_span2/satloc_dlon(2,k))
       do i=1,ispan
         ib=1+mod(ibin1+i-2,satloc_nbins)
         satloc_bins(ib,1,k)=10000       
         ib=1+mod(ib+ilatshift-1,satloc_nbins)
         satloc_bins(ib,2,k)=10000       
       enddo
!
     endif   ! check on sat type
   enddo     ! loop over k
!
   end subroutine satloc_simulate  
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine satloc_filter (kx_nfilters,kx_num,kx_filt_mcnt, & 
                             kx_type,kx_filters)
!
   implicit none
!
   integer, intent(in) :: kx_nfilters
   integer, intent(in) :: kx_num
   integer, intent(in) :: kx_filt_mcnt
   real(rkind1), intent(in) :: kx_filters(kx_nfilters,kx_num) 
   character(len=*), intent(in) :: kx_type(kx_num)
!
   integer :: k,j,n     
   integer :: count_min  ! min time-mean count for each kx*sat to be considered
   integer :: count_max  ! max time-mean count in any bin for each kx*sat
   integer :: count_sum  ! time-mean total count for each kx*sat
   integer :: ctemp(satloc_nbins) 
!
   do k=1,kx_num
     count_min=nint(kx_filters(kx_filt_mcnt,k))
     do j=1,2  ! S/N hemisphere
!
! Find max value in bins
       count_max=0
       count_sum=0
       do n=1,satloc_nbins
         count_sum=count_sum+satloc_bins(n,j,k)
         count_max=max(satloc_bins(n,j,k),count_max)
       enddo
!
! For bins with small counts, reset count to 0
       count_max=count_max/100
       do n=1,satloc_nbins
         if (count_sum < count_min .or. &
             satloc_bins(n,j,k) < count_max) then
           satloc_bins(n,j,k)=0
         endif
       enddo
!
! Fill in empty bins when adjacent bins are populated
       ctemp(:)=satloc_bins(:,j,k)
       do n=2,satloc_nbins-1
         if (satloc_bins(n-1,j,k) > 0 .and. satloc_bins(n+1,j,k) > 0 &
             .and. satloc_bins(n,j,k) == 0) then 
           ctemp(n)=(satloc_bins(n-1,j,k)+satloc_bins(n+1,j,k))/2
         endif
       enddo
!
       if (kx_type(k)(1:1) == 'P') then ! polar viewing satellite
         n=satloc_nbins
         if (satloc_bins(1,j,k) > 0 .and. satloc_bins(n-1,j,k) > 0 &
             .and. satloc_bins(n,j,k) == 0) then
           ctemp(n)=(satloc_bins(1,j,k)+satloc_bins(n-1,j,k))/2
         endif
         if (satloc_bins(2,j,k) > 0 .and. satloc_bins(n,j,k) > 0 & 
             .and. satloc_bins(1,j,k) == 0) then 
           ctemp(1)=(satloc_bins(2,j,k)+satloc_bins(n,j,k))/2
         endif
       endif
!
       satloc_bins(:,j,k)=ctemp(:)
!   
     enddo  ! loop over hemispheres
   enddo    ! loop over kx*sats
!
   end subroutine satloc_filter 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_satloc

     
