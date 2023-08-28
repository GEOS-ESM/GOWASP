  program count_amv_obs
!
! Compute counts of assimilated satwind obs (as indicated in ODS files) 
! at each NR time, p-layer, lat band, lon band, and sea/land region. 
! Also compute 6-hr avearge of such such counts for each region/layer.
! Also, optionally, save location data in a formatted file for post-processing
! to detmine distributions of closest neighbors to each obs for each obs type
! (for this application, the number of considered times should be small; e.g 
! only 1 day). 
!
! NOTE: This code does not work properly when using .ods files created from 
! an OSSE experiment, presumably because for those, the obs are read in using 
! read_prepbufr.f90 rather than read_satwind.f90. The problem is, the ks values 
! in the .ods files for the OSSE do not indicate the particular satellite, unlike
! for those generated from experiments using real data. Thus, unless the ks_subtype
! values in Table 2 of the the kx_table_in file are 0, the data will not be 
! considered. Thus, when determining results from OSSE files, data from all 
! satellites corresponding to a particular kx can only be considered together. 
!
   use m_ods_RE
!
   use m_kinds, only : rkind1
!
   use m_nr_fields_info, only : nr_fields_setup
   use m_nr_fields_info, only : field_imax,field_jmax,field_time_slots
   use m_nr_fields_info, only : field_time_delta,field_time_first
!
   use m_time_compute, only : time_compute_pack,time_compute_unpack
   use m_time_compute, only : time_compute_add
   use m_time_compute, only : rkindh
!
   use m_kx_table 
   use m_satloc
!
   implicit none
!
   logical :: lprint
   logical :: lsim_counts
   logical :: lontest1
   logical :: lontest2
   logical :: dxmin_save
!
   integer :: n_obs
   integer :: ierr
   integer :: n_missing
   integer :: argc
   integer(4) :: iargc
   integer :: ib1, ib2
   integer :: np, jbin
   integer :: nband_lat,nband_lors,nband_time
   integer :: kx,ks,kid
   integer :: k,n
   integer :: ix,jx
   integer :: idata_time
   integer :: dxmin_ncnt
   integer :: dxmin_iunit
!
   real(rkind1) :: xlat,xlon
   real(rkind1) :: xband, dlonx
   real(rkind1) :: plev,delp
   real(rkind1), allocatable :: fld_seaf(:,:)
!
! Variables concerning time keeping: 
! ana here refers to ods data file times
   integer :: nxtime1,nxtime2
   integer :: ntime_max          ! max number of ntime_counts
   integer :: ntime_count        ! number of distinct times considered
   integer :: ntime_periods      ! number of ana or BUFR dataset periods to process
   integer :: ntime_ana          ! count of analysis period
   integer :: ntime_kx           ! count of interpolation times in current ana period
   integer :: ntime_satloc       ! count of index time slot in staloc file
   integer :: ods_time_offset
   integer :: n_date             ! date input to read ods  yyyymmdd
   integer :: n_time             ! time input to read odse hhmmss
!
   real(rkindh) :: dhours_ana    ! time between analysis (or BUFR file) periods 
   real(rkindh) :: dhours_nr     ! time between nr data files used
   real(rkindh) :: dhours_kx     ! spacing between (possibly) interpolated tumes
   real(rkindh) :: thours_ana    ! ana time period hours from first time considered
   real(rkindh) :: rhours_kx     ! interpolated field time relative to ana time
   real(rkind1) :: time0_ana(6)  ! (center) time of first analysis period
   real(rkind1) :: time1_ana(6)  ! (center) time of current analysis period 
   real(rkind1) :: time0_kx(6)   ! reference t for 1st interp field in this run
   real(rkind1) :: time1_kx(6)   ! ref t for 1st interp field in current ana period
!
   character(len=14) :: cdtime0_ana  ! yyyymmddhh correspomding to time0_ana
   character(len=14) :: cdtime1_ana  ! yyyymmddhh correspomding to time1_ana
   character(len=14) :: cdtime0_kx   ! yyyymmddhh correspomding to time0_kx
   character(len=14) :: cdtime1_kx   ! yyyymmddhh correspomding to time1_kx
   character(len=4)   :: coffset     ! offset required to read some .ods files
   character(len=14)  :: cdtimeNR    ! time for NR file data containing sea fraction
!
   character(len=240) :: field_list_file
   character(len=120) :: ods_path
   character(len=120) :: ods_file_template
   character(len=240) :: ods_file       ! Must be 240 since hardwired in set_field_file_name
   character(len=30)  :: type  ! output from Get_ODS
   character(len=240) :: kx_file_in,kx_file_out
   character(len=240) :: dxmin_file_out
!
   type (ods_vect) :: ods
!
! Read and check arguments
   argc = iargc()
   if (argc .ne. 8) then
     print *,' usage must be: prog.x cdimeNR field_list_file ', &
             ' ods_path ods_file_template kx_file_in',          &
             ' kx_file_out ods_time_offset dxmin_file_out'
     stop
   endif
   call GetArg( 1_4, cdtimeNR)           ! time of NR file to get frac of sea 
   call GetArg( 2_4, field_list_file)    ! NR rc file to get grid info, templates
   call GetArg( 3_4, ods_path)           ! dir path for ODS files
   call GetArg( 4_4, ods_file_template)  ! tempate for ods file name
   call GetArg( 5_4, kx_file_in)         ! kx table
   call GetArg( 6_4, kx_file_out)        ! will contain 6-hr mean counts
   call GetArg( 7_4, coffset)            ! reference t in ods files (-180 or 0 min)
   call GetArg( 8_4, dxmin_file_out)     ! file for saving info to compute min dx
!
! ods_time_offset is the value (minutes) that denotes what the time that is 
! recorded for each datum in the ods file refers to when that time is 0 and 
! the reference is with respect to the synoptic (i.e., center) time of the 
! ods dataset. In some ods data sets, 0 as recorded for a data value means 
! the start of a 6-hr assimilation window (-180 min) whereas for other data 
! sets, it means the center of a 6-hr assimilation window. The range of 
! times recorded for data in a sample ods files should be examined to determine 
! what is meant.
!
   read (coffset,'(i4)') ods_time_offset ! reference time for data in ods files
   lprint=.true.
!
! Read kx table to get variables and descriptors of obs types
! (...,3) here means only need to read fist 3 parts of table file
   lsim_counts=.false.  
   call kx_table_read (kx_file_in,lprint,3,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to kx_table_read: ierr=',ierr 
     stop
   endif      
!
   ntime_periods=kx_times_satloc ! number of 6-hr periods to process 
   cdtime0_ana=kx_cdtime0_satloc
   call time_compute_unpack (cdtime0_ana,time0_ana)
   kx_obs_count(:,:,:,:)=0
   nband_time=0 
   delp=1000./kx_pbins  ! thickness of each p-layer band (mb)
!
! Get field containing grid dims fraction of sea at each grid point
   call nr_fields_setup ('SATWIND',field_list_file,lprint,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to nr_fields_setup: ierr=',ierr 
     stop
   endif      
   allocate (fld_seaf(field_imax,field_jmax))
   call read_nr_data_2d ('FROCEAN',cdtimeNR,fld_seaf,ierr)
!
   dhours_nr=real(field_time_delta,rkindh)
   dhours_ana=dhours_nr*real(field_time_slots,rkindh)
   dhours_kx=dhours_ana/real(kx_field_slots,rkindh)
!
   if (lprint) then 
     print ('(a,3f9.4)'),'dhours_ana,dhours_nr,dhours_kx =', &
                          dhours_ana,dhours_nr,dhours_kx
     print ('(a,i4,2x,a)'),'ods_time_offset,cdtime0_ana =', &
                            ods_time_offset,cdtime0_ana 
   endif
!
! Setup array to hold obs location array
   call satloc_setup (kx_num,kx_type,kx_locs)
!
! If dxmin file is to be output, set some values and open it
   dxmin_save = (trim(dxmin_file_out) /= 'none') 
   if (dxmin_save) then
     dxmin_iunit=22
     dxmin_ncnt=0
     open (dxmin_iunit,file=trim(dxmin_file_out),form='formatted')
     write (dxmin_iunit,'(3i4)') ntime_periods,kx_field_slots,kx_num
   endif
!
! Loop over ods data set periods
   n_missing=0 
   ntime_count=0
   ntime_max=ntime_periods*kx_field_slots
   do ntime_ana=1,ntime_periods  
     thours_ana=real(ntime_ana-1,rkindh)*dhours_ana
     call time_compute_add (thours_ana,time0_ana,time1_ana,ierr)
     call time_compute_pack (time1_ana,cdtime1_ana)
!
     call set_field_file_name ('NONE',ods_file_template, &
                               ods_path,cdtime1_ana,ods_file,ierr)
     read (cdtime1_ana(1:8), '(i8)') n_date
     read (cdtime1_ana(9:14),'(i6)') n_time
!
! Read in all obs on this ODS file for this 6-hr period
     call ods_get_RE (.true.,trim(ods_file),n_time,n_obs,ods,ierr)
     if (n_obs == 0 .or. ierr /= 0) then
       n_missing=n_missing+1
       print ('(2a)'),'ODS FILE EMPTY OR MISSING: file=',trim(ods_file)
       print *,'number of ods reading errors =',ierr
     endif
     print *,'number of obs to be considered =',n_obs
!
! Loop over separate times within each 6-hour period 
! nxtime1,2 specifies a time span, in  minutes
     do ntime_kx=1,kx_field_slots
       satloc_bins(:,:,:)=0
       ntime_count=ntime_count+1  
       nxtime1=nint((field_time_first+dhours_kx*(ntime_kx-1))*60.) 
       nxtime2=nint((field_time_first+dhours_kx*ntime_kx)*60.)
!
       rhours_kx=field_time_first+dhours_kx*(ntime_kx-1)
       call time_compute_add (rhours_kx,time1_ana,time1_kx,ierr)
       call time_compute_pack (time1_kx,cdtime1_kx)
!
       print ('(2a,4i8,f8.2,1x,a)'),'ntime_kx,ntime_count,nxtime1,',  &
                                     'nxtime2,rhours_kx,cdtime1_kx=', &
             ntime_kx,ntime_count,nxtime1,nxtime2,rhours_kx,cdtime1_kx
!
! Check each obs in the ODS file to see if it is of a type to be considered
! and, if so, what location bin it belongs in.
       do n=1,n_obs
!
! Check if this obs is an amv        
         if (ods%data%kx(n) > 239 .and. ods%data%kx(n) < 261) then 
!
! Check if obs is for u (so don't double count)
         if (ods%data%kt(n) == 4) then 
!
! Check if obs was accepted by QC
         if (ods%data%qcexcl(n) == 0) then
!
! Check if obs in time sub-interval with the 6 hr period on file 
         idata_time=ods%data%time(n)+ods_time_offset
         if (idata_time >= nxtime1 .and. idata_time < nxtime2) then
!
! Determine which kx table index this kx,ks combo corresponds to.
         kx=ods%data%kx(n) 
         ks=ods%data%ks(n)
         kid=0
         do k=1,kx_num
           if (kx == kx_list(k) .and. &
              (ks == kx_subtype(k) .or. kx_subtype(k) == 0)) then
             kid=k
           endif
         enddo
!
         if (kid > 0) then ! obs type found in list in kx_table
!
           plev=ods%data%lev(n)
           xlat=ods%data%lat(n)
           xlon=mod(ods%data%lon(n),360.)
           if (xlon < 0.)  xlon=xlon+360.
           if (xlon >= 360.) xlon=0.       ! accounts for round-off of xlon
!
! Output records for computing dx min in a separate program
           if (dxmin_save) then
             dxmin_ncnt=dxmin_ncnt+1
             write (dxmin_iunit,'(3i4,3f12.3)') &
                    kid,ntime_ana,ntime_kx,plev,xlat,xlon
           endif
!
! Determine if land or sea
           call nearest_gridpt (xlon,xlat,ix,jx)
           if (fld_seaf(ix,jx) > 0.5) then
             nband_lors=1   ! sea point
           else
             nband_lors=0   ! land point
           endif
! Determine jbin index
           xband=kx_jbins_lats*(xlat+90.)/180.
           nband_lat=1+min(int(xband),kx_jbins_lats-1) 
           jbin=nband_lat+kx_jbins_lats* &
                          (nband_lors+kx_jbins_lors*nband_time)
!
! Determine pbin index
           np=1+min(int(plev/delp),kx_pbins-1)     ! p is in this p-bin
!
! Determine satloc_bin index
           lontest1= xlon >= satloc_dlon(1,kid)
           lontest2= xlon <= satloc_dlon(3,kid)
           ib1=0  ! flag if lon not found to be in assigned observable range
!
! First check if the accepted longitudinal range for this kx does not 
! straddle the prime meridion.
           if (satloc_dlon(1,kid) < satloc_dlon(3,kid)) then 
             if (lontest1 .and. lontest2) then  
               dlonx=(xlon-satloc_dlon(1,kid))/satloc_dlon(2,kid)
               ib1=1+int(dlonx)
             endif
           else       ! range straddles the prime meridion
             if (lontest1 .or. lontest2) then  
               if (lontest1) then  ! satloc_dlon(1,kid) <= xlon < 360.
                 dlonx=(xlon-satloc_dlon(1,kid))/satloc_dlon(2,kid)
               else                ! 0. <= xlon <= satloc_dlon(1,kid)
                 dlonx=(360.+xlon-satloc_dlon(1,kid))/satloc_dlon(2,kid)
               endif
               ib1=1+int(dlonx)
             endif
           endif
!
           if (xlat >= 0.) then
             ib2=2 
           else
             ib2=1
           endif
!
! Increment counters
           if (ib1 > 0 .and. ib1 <= satloc_nbins) then
             satloc_bins(ib1,ib2,kid)=satloc_bins(ib1,ib2,kid)+1
             kx_obs_count(np,jbin,kid,1)=kx_obs_count(np,jbin,kid,1)+1
             call kx_latlon_range (kid,xlon,xlat)
           endif
!
         endif  ! check kid 
         endif  ! check time
         endif  ! check QC mark
         endif  ! check kt
         endif  ! check kx
!
       enddo    ! loop over obs in ods file
       call satloc_filter (kx_nfilters,kx_num,kx_filt_mcnt, & 
                           kx_type,kx_filters)
       call satloc_bins_write (ntime_count,ntime_max,kx_num,cdtime1_kx, &
                               kx_satloc_file,kx_list,kx_satid)
!
     enddo      ! loop over slots within an ana period window
     call ods_clean_RE (ods)
   enddo        ! loop over ana period windows
!
! Normalize counts by 6 hours and write out
   kx_obs_count(:,:,:,2)=kx_obs_count(:,:,:,1)/ntime_periods   
   call kx_table_write (kx_file_out)
!   
   call kx_table_clean
   call satloc_clean
   deallocate (fld_seaf)
!
! Write last record and close if the dxmin_file was created
   if (dxmin_save) then 
     kid=999  ! used as end of file indicator
     dxmin_ncnt=dxmin_ncnt+1
     write (dxmin_iunit,'(3i4,3f12.3)') kid,ntime_ana,ntime_kx,plev,xlat,xlon
     close (dxmin_iunit)
     print *,'dxmin_file_out closed with dxmin_ncnt=',dxmin_ncnt
   endif 
!
   print ('(a,i4)'),'Number of missing or empty ODS files =',n_missing
   print *,'Program complete'
!
   end program count_amv_obs
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
  subroutine nearest_gridpt (xlon,xlat,ix,jx)
!
  use m_kinds, only : rkind1
  use m_nr_fields_info, only : field_imax, field_jmax
  use m_nr_fields_info, only : field_lon_first
!
  implicit none
!
  integer, intent(out) :: ix,jx
  real(rkind1), intent(in) :: xlon,xlat
  real(rkind1) :: dlon,dlat,dx
!
  dlat=180./(field_jmax-1)
  dlon=360./field_imax
!
  dx=mod(xlon-field_lon_first,360.)
  ix=min(1+nint(dx/dlon),field_imax)
  jx=1+nint((xlat+90.)/dlat)
! 
  end subroutine nearest_gridpt
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine check(status,loc)
!
   use netcdf
   implicit none 
   integer, intent(in) :: status
   character(*), intent(in), optional :: loc
   if (status /= nf90_noerr) then
     if (present(loc)) print *,'Error ',status,' at ',loc
     stop
   end if
!
   end subroutine check
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_nr_data_2d (f_name,cdtime,field_in,iers)
!
!  Read 2d field on netcdf file. 
!
   use netcdf           ! for reading the NR files
   use m_kinds, only : rkind1
   use m_nr_fields_info, only : field_imax, field_jmax, field_num_2d
   use m_nr_fields_info, only : field_lon_first, field_common_path
   use m_nr_fields_info, only : field_names, field_types, field_files
!
     implicit none
!
     character(len=*), intent(in) :: cdtime
     character(len=*), intent(in) :: f_name                
     integer, intent(out) :: iers
     real(rkind1), intent(out) :: field_in(field_imax,field_jmax)
!                         
     integer :: ier
     integer :: imx, jmx
     integer :: ncid,varid,id
     integer :: id_start(3)
     integer :: id_count(3)
     character(len=120) :: c_notice
     character(len=240) :: file_name
     character(len=*), parameter :: subname='read_nr_data_2d'
!
     iers=0
     call find_name (field_num_2d,field_names(:,1,2),.false., &
                     subname,f_name,id) 
     if (id == 0) then 
       iers=1000
       print *,'FATAL ERROR: Required name not found in ',subname
     else 
       call set_field_file_name (field_names(id,2,2),field_files(id,2),  &
                               field_common_path,cdtime,file_name,ier) 
       iers=iers+ier
     endif 
!
     if (iers == 0) then
       c_notice='Opening file for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid),         &  
                  trim(c_notice),ier)
!
! Get dimension information to check
       call check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
       call check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                    'nf90_inq 2',ier)
       call check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
       call check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                    'nf90_inq 4',ier)
       if (imx /= field_imax .or. jmx /= field_jmax) then 
         print *,'Grid dimension mismatch in routine : read_shmem_data'
         print *,'file_name=',trim(file_name)
         print *,'imax, jmax on file = ',imx,jmx
         print *,'imax, jmax in program = ',field_imax,field_jmax
         iers=iers+10
       endif
!
       c_notice='Getting vari for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_inq_varid(ncid,trim(field_names(id,2,2)),varid), &
                  trim(c_notice),ier)
!
       c_notice='reading field for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_get_var(ncid,varid,field_in(:,:)),trim(c_notice),ier) 
!
       c_notice='Closing file for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_close(ncid),trim(c_notice),ier)
     endif
!
   end subroutine read_nr_data_2d 
