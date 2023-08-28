   module m_pert_loop
!
! Module that loops through messages and reports on bufr file and 
! calls appropriate routines to add errors
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use m_kinds, only : rkind1, rkind2
!
   private
   public :: pert_loop_setup 
   public :: pert_loop_do
   public :: pert_loop_clean
!
   integer, parameter :: obs_max_chan_or_levs=1000
   integer, parameter :: obs_max_fields=2
   integer, parameter :: luin_clrsky=25 
!
   logical :: lrad
   logical :: lrelhum
   logical :: leof ! EOF flag 
   logical :: clear_sky_file ! true if corresp .bin file of clear-sky BT used
   logical :: obs_psflag(obs_max_chan_or_levs)
!
   integer :: luin, luout
   integer :: nvalues
   integer :: index_itype,index,nlevs,index_lat,index_lon,index_nlevs
   integer :: obs_nfields
!
   real(rkind1) :: clrsky_bt(obs_max_chan_or_levs)
   real(rkind2) :: bmiss_data
   real(rkind2) :: obs_data(obs_max_chan_or_levs,obs_max_fields)
   real(rkind2) :: obs_levs(obs_max_chan_or_levs)
   character(len=4)  :: obs_file_type   ! 'GTXT' if obs file generic rad text
   character(len=12) :: obs_file_format ! 'formatted' or 'unformatted'
   character(len=240) :: clear_sky_name ! name for clear-sky BT .bin file
   character(len=*), parameter :: myname='m_pert_loop'
!
   contains
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_setup (bufr_in_file,bufr_out_file,bufr_tab_file, &
                               lprint,dtype,clear_sky_suffix,ierr)
!
! Calls routines to set some variables and open read/write BUFR files for
! requested obs type.  
!
   use m_count_types, only : count_types_setup
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_tab_file
   character(len=*), intent(in) :: bufr_out_file
   character(len=*), intent(in) :: clear_sky_suffix  ! add to sim_obs file name for clear BT file
!
! leof is always false, unless data type is GENRADTXT, in which case
! its value is true when the obs file EOF is reached.
   leof=.false.            ! default value
   obs_file_type='BUFR'    ! default value for obs data file type
   obs_file_format='unformatted'  ! default value for obs data format
!   
   if (trim(dtype) == 'PREPBUFR') then    
     call pert_loop_setup_prepbufr (bufr_in_file,bufr_out_file, &    
                                    lprint,ierr)
   elseif (trim(dtype) == 'GPSRO') then 
     call pert_loop_setup_gpsro (bufr_in_file,bufr_out_file, &
                                 lprint,ierr)    
   else   ! rad obs assumed
     call pert_loop_setup_rad (bufr_in_file,bufr_out_file, &
                               bufr_tab_file,lprint,dtype,ierr)
!
! Open additional .bin file of corresponding clear-sky BT if statistics of 
! errors to be added depend on differences between clear and all-sky BT
! (The all-sky values are assumed to be in the BUFR files).     
     if (trim(clear_sky_suffix) /= 'none') then
       clear_sky_file=.true.
       clear_sky_name=trim(bufr_in_file)//trim(clear_sky_suffix)
print *,'QQQ1 ',trim(bufr_in_file)
print *,'QQQ2 ',trim(clear_sky_suffix)
print *,'QQQ3 ',trim(clear_sky_name)

       open (luin_clrsky,file=trim(clear_sky_name),form='unformatted') 
       if (lprint) then
         print ('(2a)'), 'File of clear-sky BT opened: ',trim(clear_sky_name)
       endif
     else
       clear_sky_file=.false.
     endif
!
   endif
!
   call count_types_setup
!
   end subroutine pert_loop_setup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_setup_prepbufr (bufr_in_file,bufr_out_file, &
                                        lprint,ierr)
!
! Set some variables and open read/write BUFR files for conventional obs
! Includes any obs in PREPBUFR format files.
!
   use m_conv_names, only : conv_names_setup
   use m_conv_names, only : conv_nhead, conv_nfields, conv_max_levs
   use m_conv_names, only : ixob,iyob,ityp
!
   use m_bufr_conv, only  : conv_rw_setup
   use m_bufr_conv, only  : bufr_unit_in, bufr_unit_out, conv_bmiss
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_out_file
!
   integer :: ier
   character(len=*), parameter :: mysub=myname//'::pert_loop_setup_prepbufr'
!
   call conv_names_setup (ier)
   ierr=ier
   if (lprint .and. ierr /= 0) then
     print *,'Error detected in call to conv_names_setup: ierr=',ierr
   endif
!
   index_itype=ityp
   index_lat=iyob
   index_lon=ixob
   obs_nfields=2   ! either T,q or u,v (p is carried in different array)
   lrad=.false.
!
! Open bufr files 
! Also setup some indexes in the m_bufr_conv module
! 'SIMERR' signifies that call is in conjunction with simulation of error
   call conv_rw_setup (bufr_in_file,bufr_out_file,lprint,'SIMERR',ier)
   ierr=ierr+abs(ier)
   luin=bufr_unit_in
   luout=bufr_unit_out
   bmiss_data=real(conv_bmiss,rkind2)
!
! Compress messages in BUFR files written (so several records fit in 1 msg)
   call maxout(200000)            ! increase size of ouput bufr
!
   end subroutine pert_loop_setup_prepbufr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_setup_gpsro (bufr_in_file,bufr_out_file, &
                                     lprint,ierr)
!
! Set some variables and open read/write BUFR files for GPSRO obs
!
   use m_gpsro_names, only : gpsro_names_setup
   use m_gpsro_names, only : bblat, bblon, bbsaid, bbnlold
!
   use m_bufr_gpsro, only  : gpsro_rw_setup
   use m_bufr_gpsro, only  : bufr_unit_in, bufr_unit_out, obs_bmiss
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_out_file
!
   integer :: ier
   character(len=*), parameter :: mysub=myname//'::pert_loop_setup_gpsro'
!
   call gpsro_names_setup (ier)
   ierr=ier
   if (lprint .and. ierr /= 0) then
     print *,'Error detected in call to gpsro_names_setup: ierr=',ierr
    endif
!
   index_itype=bbsaid
   index_lat=bblat
   index_lon=bblon
   index_nlevs=bbnlold
   obs_nfields=1   ! only 1 data type (bending angle)
   lrad=.false.
!
! Open bufr files 
   call gpsro_rw_setup (bufr_in_file,bufr_out_file,lprint,'SIMERR',ier)
   ierr=ierr+abs(ier)
   luin=bufr_unit_in
   luout=bufr_unit_out
   bmiss_data=obs_bmiss
!
   end subroutine pert_loop_setup_gpsro
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_setup_rad (bufr_in_file,bufr_out_file, &
                                   bufr_tab_file,lprint,dtype,ierr)
!
! Set some variables and open read/write BUFR files for radiance obs
!
   use m_rad_obs_arrays, only : rad_obs_arrays_setup 
   use m_rad_obs_arrays, only : obs_info_num,obs_info_names,obs_info_hdr
   use m_bufr_rad, only : read_write_obs_rad 
   use m_bufr_rad, only : bmiss
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_tab_file
   character(len=*), intent(in) :: bufr_out_file
!
   logical, parameter :: lstop=.false.
   logical :: ldum   ! EOF flag argument not used here
   integer :: nobs
   integer :: ier
   integer :: luin_tab
   character(len=*), parameter :: mysub=myname//'::pert_loop_setup_rad'
!
   luin=20
   luout=21
   luin_tab=22  ! default value; changed if separate bufr_tab_file required
!
! Setup max sizes for arrays to hold read observation information 
   call rad_obs_arrays_setup (700,5,5,100,5)  
!
! Open both read and write obs data files
! Extract some information from read/write routines
! If File s generic text, read and write file header record
   call open_obs_files (luin,luin_tab,luout,.true.,.true.,       &
                        lprint,dtype,bufr_in_file,bufr_tab_file, &
                        bufr_out_file,'none',obs_file_format,    &
                        obs_file_type,ierr)
!  
! Get index for lat in header     
   call find_name (obs_info_hdr,obs_info_names,lstop,' ','CLAT',index_lat)  
   if (index_lat == 0) then
     call find_name (obs_info_hdr,obs_info_names,lstop,mysub,'CLATH',index_lat)
   endif
!
! Get index for lon in header     
   call find_name (obs_info_hdr,obs_info_names,lstop,' ','CLON',index_lon)  
   if (index_lon == 0) then
     call find_name (obs_info_hdr,obs_info_names,lstop,mysub,'CLONH',index_lon)
   endif
!
! Get index for sat or, if AQUA, instrument
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then 
     call find_name (obs_info_hdr,obs_info_names,lstop,mysub,'SIID',index_itype)
   else
     call find_name (obs_info_hdr,obs_info_names,lstop,mysub,'SAID',index_itype)
   endif
!
   bmiss_data=bmiss
   obs_nfields=1   ! only 1 data type (either Tb or Rad)
   lrad=.true.
   ierr=0
!
   end subroutine pert_loop_setup_rad
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_do (lprint,change_year,ltest,crtm_coef_dir,dtype)
!
!  Within a do-loop over all obs on file. sequence calls to read original 
!  obs, create and add perturbations, and write new obs. 
!   
   use m_bufr_gpsro, only : gpsro_rw
   use m_bufr_gpsro, only : gpsro_info, gpsro_values, obs_nlevs 
   use m_gpsro_names, only : obs_info_num, obs_nvalues
!
   use m_bufr_conv, only : conv_rw    
   use m_bufr_conv, only : conv_nlevs
   use m_bufr_conv, only : conv_info,conv_values
!
   use m_bufr_rad, only : read_write_obs_rad   
   use m_bufr_rad, only : read_write_gmi_1st_msg
   use m_rad2bt, only : rad2bt_setup
   use m_rad2bt, only : rad2bt_transforms
   use m_rad2bt, only : rad2bt_cleanup 
   use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
   use m_rad_obs_arrays, only : obs_n_channels
!
   use m_obs_error_table, only : error_table_find_corr_id
   use m_obs_error_table, only : et_n_err1, et_n_err2, et_n_err3
   use m_obs_error_table, only : et_err_itype, et_err_tab
   use m_obs_error_table, only : et_icolumn_tq, et_icolumn_uv, et_icolumn_ps
   use m_obs_error_table, only : et_icolumn_gpsro, et_icolumn_rad
   use m_obs_error_table, only : et_pert_fac, et_vcorr_dist
   use m_obs_error_table, only : et_itypes_corr
!
   use m_count_types, only : count_subtypes
   use m_count_types, only : count_subtypes_get_rpts
!
   use m_obs_pert, only : pert_find_itype
   use m_obs_pert, only : pert_obs
   use m_obs_pert, only : pert_ps
!
   implicit none
!
   logical, intent(in) :: ltest
   logical, intent(in) :: lprint
   integer, intent(in) :: change_year     
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: crtm_coef_dir
!
! local variables
!
   logical :: lprint_sample
   logical :: lmsg1
   logical :: lread
   logical :: lerror  
   logical :: generic_ireadsb ! generalization of ireadsb BUFR lib function
   logical :: generic_ireadmg ! generalization of ireadmg BUFR lib function
!
   integer, parameter :: nerrors=15
   integer :: err_clrsky
   integer :: err_index   ! index for data subtype
   integer :: nobs
   integer :: idate   
   integer :: ierrors(nerrors)
   integer :: ier
   integer :: nrpts
   integer :: itype
   integer :: olevs
   integer :: onums,onums1
   integer :: n_mesg
   integer :: icolumns(2)
   integer :: hcorr_id   ! id of group of kx or ks to have horiz corr perts
   integer :: prev_itype ! previous sat id (used for IASI rad2tb setup)
!    
   real(rkind2) :: vcorr_dist_flds(10)
   real(rkind2) :: lat, lon
   character(len=8) :: psubset, subset
!
! Begin loop to read input observations
   nobs=0  
   psubset =''
   ierrors(:)=0
   n_mesg=0
   itype=0  
   prev_itype=0
!
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'PREPBUFR' .or.  &
       trim(dtype) == 'GPSRO' .or. obs_file_type == 'GTXT' ) then
     lmsg1=.false.   ! do not need to create first message with 0 reports
   else
     lmsg1=.true.    
   endif
!
! Special instructions to read/write 1st message in GMI bufr data
    if (trim(dtype) == 'GMI') then
      call read_write_gmi_1st_msg (luin,.true.,idate,ier)
      if (change_year > 0) then
        idate=idate+(change_year-(idate/1000000))*1000000
      endif
      call read_write_gmi_1st_msg (luout,.false.,idate,ier)
    endif
!
! Loop over messages
!
   do while (generic_ireadmg(luin,leof,obs_file_type,subset,idate))
!
     if (change_year > 0) then
       idate=idate+(change_year-(idate/1000000))*1000000
     endif
!
! If new data source (compare new subset name with previous subset name)
     if (subset .ne. psubset ) then
       if (lprint) then   
         print *,'Processing subset ',subset,' for idate ',idate
       endif
!
! Create first message as just date/time stamp with no reports
! GSI is designed to skip the first message on file for some data types, 
! so this must be message filled if the first set of reports are to be 
! read by GSI
       if (lmsg1) then
         call openmg(luout,subset,idate)  ! idate here is yyyymmddhh
         call minimg(luout,0)             ! add minutes=00 to idate write
         call closmg(luout)
         lmsg1=.false. 
       endif
       psubset=subset ! previous subset name updated to current subset name
     endif
!
     if (obs_file_type == 'BUFR') then
       call openmb (luout,subset,idate)
     endif
     do while (generic_ireadsb(luin,leof,obs_file_type))
!
! Used to limit printing of sample output with the hope of printing a
! sample that will include all prepbufr data subtypes (kx values)
! This will be unsuccessful unless occurances of each data type occur 
! consecutively in the bufr file for at least the first few of each type.
       lprint_sample=ltest .and. lprint
       if (trim(dtype) == 'PREPBUFR') then  
         call count_subtypes_get_rpts (itype,nrpts)
         if (nrpts > 4) then 
           lprint_sample=.false.
         endif
       else
         if (nobs > 2) then  
           lprint_sample=.false.
         endif
       endif
!
! Read in next obs record
       lread=.true. ! indicates bufr data to be read
       if (trim(dtype) == 'PREPBUFR') then  
         nobs=nobs+1
         call conv_rw (lprint,lprint_sample,'EITHER',nobs, &
                       lread,nerrors,ierrors,lerror,ier)
       elseif (trim(dtype) == 'GPSRO') then
         nobs=nobs+1  
         call gpsro_rw (lprint_sample,lread,'9999999999',nobs,nerrors,ierrors, &
                        lerror,ier) 
       else  ! rad data  (nobs incremented in call to routine)
         call read_write_obs_rad (luin,lprint_sample,dtype,nobs,.true., &
                                  lread,leof,ier)
!
! Read additional .bin file of corresponding clear-sky BT if required
         if (clear_sky_file) then
           call pert_loop_clear_sky (err_clrsky)
           if (err_clrsky /= 0) then
             print *,'Sequence differs in clear-sky .bin and all-sky bufr files'
             stop
           endif
         else
           err_clrsky=0
         endif
!
       endif
!
! If EOF on input file not encountered, set olevs to number of valid 
! data records (levels or channels) in obs data report
       if (.not. leof) then 
         if (trim(dtype) == 'PREPBUFR') then  
           olevs=conv_nlevs
         elseif (trim(dtype) == 'GPSRO') then  
           olevs=obs_nlevs
         else   ! radiance  
           olevs=obs_n_channels
         endif
       else
         olevs=0  ! default value indicating report has no data va
       endif
!
! If EOF on input file not encountered and 1 or more valid data records present
! Specify parameters to locate required information
! Also copy data from BUFR array format to pert array format 
       if (.not. leof .and. olevs > 0) then 
         vcorr_dist_flds(1:obs_nfields)=et_vcorr_dist(1:obs_nfields)
         lrelhum=.false.
         if (trim(dtype) == 'PREPBUFR') then  
           itype=nint(conv_info(index_itype))
           lat=conv_info(index_lat)
           lon=conv_info(index_lon)
           onums=olevs*2
           if (itype < 200) then 
             lrelhum=.true.
             icolumns(1:2)=et_icolumn_tq(1:2)
             call pert_loop_copy_values (olevs,dtype,'MASS','B2P')
           else
             vcorr_dist_flds(1:2)=et_vcorr_dist(3:4)  ! set values for u,v
             icolumns(1:2)=et_icolumn_uv(1:2)
             call pert_loop_copy_values (olevs,dtype,'WIND','B2P')
           endif
!
         elseif (trim(dtype) == 'GPSRO') then  
           itype=nint(gpsro_info(index_itype))
           onums=olevs
           lat=gpsro_info(index_lat)
           icolumns(1)=et_icolumn_gpsro
           call pert_loop_copy_values (olevs,dtype,'ALL','B2P')
!
         else   ! radiance  
           itype=nint(obs_info(index_itype)) 
           lat=obs_info(index_lat)
           lon=obs_info(index_lon)
           onums=olevs
           icolumns(1)=et_icolumn_rad
           call pert_loop_copy_values (olevs,dtype,'ALL','B2P')
         endif
!
! Count the number of obs for each separate obs typ
         if (obs_psflag(obs_max_chan_or_levs)) then
           onums1=onums+1 ! value augmented by 1 to indicate ps present
         else
           onums1=onums
         endif
         call count_subtypes (itype,onums1)
!
! Find the error table index for this data subtype 
! NOTE: Search of the sat info table will have assigned the sat id numbers 
! in this table to the id found for the first occurance of the corresponding 
! sat name found in the sat info file, so if that id is set wrong, it will 
! be wrongly set for all later sats of the same name.
         if (trim(dtype) == 'GPSRO') then 
           err_index=1
         else
           call pert_find_itype (et_n_err3,itype,et_err_itype,err_index)
         endif
         if (err_index == 0) then
           print *,' '
           print ('(a,i4,a)'),'REQUESTED SUBTYPE itype=',itype, &
                   ' NOT FOUND AMONG SUBTYPES READ IN ERROR TABLE'
           print ('(a,20i5)'),'TYPES FOUND IN TABLE:',et_err_itype(1:et_n_err3)
           print ('(a)'),'EXECUTION HALTED'
           stop
         endif
!      
! If this is an observation type that should have horizontally corrrelated 
! errors, then determine what set of correlation parameters should be used.
         if (et_itypes_corr > 0) then
           call error_table_find_corr_id (itype,hcorr_id)
         else
           hcorr_id=0
         endif
!
! For IASI or CRIS, convert radiances to brightness temperatures here
! This requires use of appropriate CRTM coefficient files. Since these 
! differ for IASI on metop-a and metop-b, the file must be changed if 
! the sat id (itype) differs from the previously considered id.
         if (trim(dtype) == 'IASI' .or. trim(dtype) == 'CRIS' .or. &
             trim(dtype) == 'CRISFSR') then
           if (prev_itype /= itype) then 
             if (prev_itype /= 0) then   ! first remove previously set file
               call rad2bt_cleanup (ier)  
             endif 
             call rad2bt_setup (crtm_coef_dir,dtype,itype,ier)
           endif
           call rad2bt_transforms (olevs,obs_data,'R2TB')
           prev_itype=itype             ! reset previous value to current one
         endif
!
         if (obs_psflag(obs_max_chan_or_levs)) then ! data for valid ps value
           call pert_ps (et_n_err1,et_n_err2,obs_max_chan_or_levs,nobs, &
                         olevs,et_pert_fac,et_err_tab(:,:,err_index),   &
                         obs_levs,obs_psflag,lprint_sample)
         endif
!
         call pert_obs (et_n_err1,et_n_err2,obs_max_chan_or_levs,          &
                        obs_max_fields,nobs,                               &
                        olevs,obs_nfields,icolumns,bmiss_data,et_pert_fac, &
                        et_err_tab(:,:,err_index),obs_levs,obs_data,       &
                        vcorr_dist_flds,lat,lon,hcorr_id,lprint_sample,    &
                        lrelhum,dtype,lrad)
!
! For IASI or CRIS, convert brightness temperatures to radiances here
         if (trim(dtype) == 'IASI' .or. trim(dtype) == 'CRIS' .or. &
             trim(dtype) == 'CRISFSR') then
           call rad2bt_transforms (olevs,obs_data,'TB2R')
         endif
!
! Copy data from pert array format to BUFR array format
         if (trim(dtype) == 'PREPBUFR') then  
           if (itype < 200) then 
             call pert_loop_copy_values (olevs,dtype,'MASS','P2B')
           else
             call pert_loop_copy_values (olevs,dtype,'WIND','P2B')
           endif
         elseif (trim(dtype) == 'GPSRO') then  
           call pert_loop_copy_values (olevs,dtype,'ALL','P2B')
         else  ! Assumed to be radiance data
           call pert_loop_copy_values (olevs,dtype,'ALL','P2B')
         endif
!
! Write next obs record
         lread=.false.  ! indicates bufr data to be written
         if (trim(dtype) == 'PREPBUFR') then  
           call conv_rw (lprint,lprint_sample,'EITHER',nobs, &
                         lread,nerrors,ierrors,lerror,ier)
         elseif (trim(dtype) == 'GPSRO') then  
           call gpsro_rw (lprint_sample,lread,'9999999999',nobs,nerrors, &
                          ierrors,lerror,ier) 
         else  ! rad data
           call read_write_obs_rad (luout,lprint_sample,dtype,nobs,.true., &
                                    .false.,leof,ier)
         endif
!
       endif  ! check on leof and olevs (no EOF and >0 records in report)
!
! End loop over buffer read
     enddo               ! loop over do while ( ireadsb(luin) .eq. 0 )
     if (obs_file_type == 'BUFR') then
       call closmg(luout)  
     endif
     n_mesg=n_mesg+1       ! count messages read in original file
     enddo                   ! loop over (ireadmg(luin,subset,idate).eq. 0)
     if (obs_file_type == 'BUFR') then   ! BUFR file format
       call closmg(luout)
       call closbf(luout)        
       call closbf(luin)
     else                                  ! not BUFR file format
       close (luout)
       close (luin)
     endif
!
     if (clear_sky_file) then
       close (luin_clrsky)
     endif
!
     if (trim(dtype) == 'IASI' .or. trim(dtype) == 'CRIS' .or. &
         trim(dtype) == 'CRISFSR') then 
       call rad2bt_cleanup (ier)
     endif
!
! Print out number of obs records processed
   if (lprint) then
     print *,' '
     print ('(i7,2a)'),n_mesg,' observation message-groups read all types'
     print ('(i7,2a)'),ierrors(1),' observation reports found with no',  &
                       ' corresponding entry in error table'
     print ('(i7,a)' ),ierrors(9),' read/write errors detected'
     print ('(i7,2a)'),nobs,' observation reports processed for type = ' &
                       ,trim(dtype)
   endif
!
   end subroutine pert_loop_do
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_copy_values (olevs,dtype,stype,func)
!
! Copy observation values and obs headers between arrays used by 
! create_error routines and obs read/write routines. 
!
   use  m_conv_names, only : bbp, bbt, bbq, bbu, bbv, bbc, bbzq
   use m_bufr_conv, only : conv_values
!
   use m_gpsro_names, only : bbimpp, bbbang, bbcurve
   use m_bufr_gpsro, only : gpsro_values, gpsro_info  
!
   use m_rad_obs_arrays, only : obs_values, obs_channels
!
   implicit none
!
   integer, intent(in) :: olevs
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: stype
   character(len=*), intent(in) :: func
!  
   integer :: n
!
   if (trim(func) == 'B2P') then  ! copy BUFR to pert arrays
!
     obs_psflag(:)=.false. ! default set to no valid ps values in data
!
     if (trim(dtype) == 'PREPBUFR') then
       do n=1,olevs
         obs_levs(n)=conv_values(bbp,n)
         if (stype == 'WIND') then           ! WIND data
           obs_data(n,1)=conv_values(bbu,n)
           obs_data(n,2)=conv_values(bbv,n)
         else                                ! MASS data
           obs_data(n,1)=conv_values(bbt,n)
           obs_data(n,2)=conv_values(bbq,n)
!
! check if valid ps value is present and set flag
           if (nint(conv_values(bbc,n)) == 0 .and. &
                            conv_values(bbzq,n) < 3.5) then
             obs_psflag(n)=.true.                    ! value valid at this lev
             obs_psflag(obs_max_chan_or_levs)=.true. ! means a valid value found
           endif             
         endif
       enddo
!
! For gpsro errors are specified as function of height above the local radius
! of curvature specified for this location in the data header.
     elseif (trim(dtype) == 'GPSRO') then
       do n=1,olevs
         obs_data(n,1)=gpsro_values(bbbang,n)
         obs_levs(n)=gpsro_values(bbimpp,n)-gpsro_info(bbcurve) ! z=r-refsurf
       enddo
!
     else  ! assumed to be radiance     
       do n=1,olevs
         obs_data(n,1)=obs_values(n,1)
         obs_levs(n)=obs_channels(n,1)
       enddo
     endif
!
   else
!
     if (trim(dtype) == 'PREPBUFR') then
       do n=1,olevs
         conv_values(bbp,n)=obs_levs(n) 
         if (stype == 'WIND') then           ! WIND data
           conv_values(bbu,n)=obs_data(n,1)
           conv_values(bbv,n)=obs_data(n,2)
           conv_values(bbt,n)=bmiss_data
           conv_values(bbq,n)=bmiss_data
         else                                ! MASS data
           conv_values(bbt,n)=obs_data(n,1)
           conv_values(bbq,n)=obs_data(n,2)
           conv_values(bbu,n)=bmiss_data
           conv_values(bbv,n)=bmiss_data
         endif
       enddo
!
     elseif (trim(dtype) == 'GPSRO') then
       do n=1,olevs
         gpsro_values(bbbang,n)=obs_data(n,1)
       enddo
!
     else  ! assumed to be radiance     
       do n=1,olevs
         obs_values(n,1)=obs_data(n,1)
       enddo
     endif
!
   endif
!
   end subroutine pert_loop_copy_values 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_clear_sky (clrsky_err)
!
! Read corresponding set of clear-sky BT if these will be used to compare
! with all-sky BT values for determining all-sky error statistics
!
   use m_rad_obs_arrays, only : obs_n_channels, obs_info
!
   implicit none
!
   integer, intent(out) :: clrsky_err
   integer :: clrsky_n
   real(rkind2) :: clrsky_lat,clrsky_lon,clrsky_time
!  
   clrsky_err=0
   read (luin_clrsky) clrsky_n,clrsky_lat,clrsky_lon,clrsky_time, &
                      clrsky_bt(1:obs_n_channels)
   if (abs(clrsky_lat-obs_info(index_lat)) > 0.02 .or. &
       abs(clrsky_lon-obs_info(index_lon)) > 0.02 ) then
     print ('(2a,i6,4f8.2)'),'Obs locations in clear-sky .bin and all-sky .bufr ', &
             'files do not correspond: n, lats, lons=',clrsky_n,clrsky_lat,        &
             obs_info(index_lat),clrsky_lon,obs_info(index_lon)
     clrsky_err=2
   endif          

   if (mod(clrsky_n,1000) == 1) then
     print ('(a,i6,4f8.2)'),'QQQ ',clrsky_n,clrsky_lat,        &
             obs_info(index_lat),clrsky_lon,obs_info(index_lon)
   endif          

!
   end subroutine pert_loop_clear_sky
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine pert_loop_clean 
!
! Deallocate some arrays
!
   use m_rad_obs_arrays, only : rad_obs_arrays_clean
!
   if (lrad) then 
     call rad_obs_arrays_clean
   endif
!
   end subroutine pert_loop_clean
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_pert_loop
