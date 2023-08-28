   module m_read_loop
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1, rkind2
!
   private
   public :: read_loop_setup 
   public :: read_loop_do
   public :: read_loop_clean
!
   logical :: lrad
   logical :: lrelhum
   logical :: leof ! EOF flag required when reading text obs data files
   integer, parameter :: obs_max_chan_or_levs=1000
   integer, parameter :: obs_max_fields=2
   integer, public :: nobs
   integer :: i_file
   integer :: luin, luout
   integer :: nvalues
   integer :: index_itype,index,nlevs,index_lat,index_lon,index_nlevs
   integer :: index_curve  ! index for gpsro local radius of curvature
   integer :: obs_nfields
   logical :: obs_psflag(obs_max_chan_or_levs)
   real(rkind2) :: gpsro_curve   ! local radius of curvature of geoid
   real(rkind2) :: obs_ps
   real(rkind2) :: bmiss_data
   real(rkind2) :: obs_data(obs_max_chan_or_levs,obs_max_fields)
   real(rkind2) :: obs_levs(obs_max_chan_or_levs)
   character(len=4)  :: obs_file_type   ! 'GTXT' if obs file generic rad text
   character(len=12) :: obs_file_format ! 'formatted' or 'unformatted'
   character(len=*), parameter :: myname='m_read_loop'
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_setup (bufr_in_file,bufr_out_file,bufr_tab_file, &
                               crtm_coef_dir,lprint,dtype,i_file_in,ierr)
!
   use m_count_types, only : count_types_setup
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(in)  :: i_file_in
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_tab_file
   character(len=*), intent(in) :: bufr_out_file
   character(len=*), intent(in) :: crtm_coef_dir
!   
   i_file=i_file_in
   call count_types_setup
!
! leof is always false, unless data type is GENRADTXT, in which case
! its value is true when the obs file EOF is reached.
   leof=.false.  ! default value
!
   if (trim(dtype) == 'PREPBUFR') then    
     call read_loop_setup_prepbufr (bufr_in_file,bufr_out_file, &    
                                    lprint,ierr)
   elseif (trim(dtype) == 'GPSRO') then 
     call read_loop_setup_gpsro (bufr_in_file,bufr_out_file, &
                                 lprint,ierr)
   else   ! rad obs assumed
     call read_loop_setup_rad (bufr_in_file,bufr_out_file,bufr_tab_file, &
                               crtm_coef_dir,lprint,dtype,ierr)    
   endif
!
   end subroutine read_loop_setup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_setup_prepbufr (bufr_in_file,bufr_out_file, &
                                        lprint,ierr)
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
   character(len=*), parameter :: mysub=myname//'::read_loop_setup_prepbufr'
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
   call conv_rw_setup (bufr_in_file,bufr_out_file,lprint,'SIMOBS',ier)
   luin=bufr_unit_in
   luout=bufr_unit_out
   bmiss_data=real(conv_bmiss,rkind2)
!
! Compress messages in BUFR files written (so several records fit in 1 msg)
   call maxout(200000)            ! increase size of ouput bufr
!
   end subroutine read_loop_setup_prepbufr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_setup_gpsro (bufr_in_file,bufr_out_file, &
                                     lprint,ierr)
!
   use m_gpsro_names, only : gpsro_names_setup
   use m_gpsro_names, only : bblat, bblon, bbsaid, bbnlold, bbcurve
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
   character(len=*), parameter :: mysub=myname//'::read_loop_setup_gpsro'
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
   index_curve=bbcurve
   index_nlevs=bbnlold
   obs_nfields=1   ! only 1 data type (bending angle)
   lrad=.false.
!
! Open bufr files 
   call gpsro_rw_setup (bufr_in_file,bufr_out_file,lprint,'SIMOBS',ier)
   luin=bufr_unit_in
   luout=bufr_unit_out
   bmiss_data=obs_bmiss
!
   end subroutine read_loop_setup_gpsro
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_setup_rad (bufr_in_file,bufr_out_file,bufr_tab_file, &
                                   crtm_coef_dir,lprint,dtype,ierr)
!
   use m_rad_obs_arrays, only : rad_obs_arrays_setup 
   use m_rad_obs_arrays, only : obs_info_num,obs_info_names,obs_info_hdr
   use m_iasi_radiance, only : iasi_radiance_setup
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
   character(len=*), intent(in) :: crtm_coef_dir
!
   logical, parameter :: lstop=.false.
   logical :: ldum   ! EOF flag argument not used here
   integer :: ier
   integer :: luin_tab
   character(len=*), parameter :: mysub=myname//'::read_loop_setup_rad'
!
   luin=20
   luout=21
   luin_tab=22  ! default value; changed if separate bufr_tab_file required
!
! Setup max sizes for arrays to hold read observation information 
   call rad_obs_arrays_setup (700,2,2,100,5)  
!
! Open files and get some information from them or from reading routines.
! 0 and 'none' 
   call open_obs_files (luin,luin_tab,0,.true.,.true.,           &   
                        lprint,dtype,bufr_in_file,bufr_tab_file, &
                        'none','none',obs_file_format,           &
                        obs_file_type,ier)
!
! If IASI, get info regarding transforms between radiance and brightness T
   if (trim(dtype) == 'IASI') then
     call iasi_radiance_setup (crtm_coef_dir,ierr)
   endif
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
!
   end subroutine read_loop_setup_rad
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_do (lprint,change_year,ltest,dtype)
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
   use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
   use m_rad_obs_arrays, only : obs_n_channels
   use m_iasi_radiance, only: iasi_radiance_transforms
!
   use m_count_types, only : count_subtypes
   use m_count_types, only : count_subtypes_get_rpts
!
   use m_saved_data, only : sd_copy
!
   implicit none
!
   logical, intent(in) :: ltest
   logical, intent(in) :: lprint
   integer, intent(in) :: change_year     
   character(len=*), intent(in) :: dtype
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
   integer :: err_index   ! index for data subtype
   integer :: idate   
   integer :: ierrors(nerrors)
   integer :: ier
   integer :: nrpts
   integer :: itype
   integer :: olevs
   integer :: onums,onums1
   integer :: n_mesg
   integer :: hcorr_id   ! id of group of kx or ks to have horiz corr perts
   integer :: ireadsb,ireadmg ! reading functions in bufr lib 
!    
   real(rkind2) :: lat, lon
   character(len=8) :: psubset, subset
!
! Begin loop to read input observations
   nobs=0  
   psubset =''
   ierrors(:)=0
   n_mesg=0
   itype=0  
!
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'PREPBUFR' .or.  &
       trim(dtype) == 'GPSRO') then
     lmsg1=.false.   ! do not need to create first message with 0 reports
   else
     lmsg1=.true.    
   endif
!
!QQQQQ
   obs_file_type='BUFR'

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
       psubset=subset ! previous subset name updated to current subset name
     endif
!
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
       endif
!
       if (.not. leof) then 
!
! Specify parameters to locate required information
! Also copy data from BUFR array format to pert array format
       lrelhum=.false.
       if (trim(dtype) == 'PREPBUFR') then  
         itype=nint(conv_info(index_itype))
         olevs=conv_nlevs
         lat=conv_info(index_lat)
         lon=conv_info(index_lon)
         onums=olevs*2
         if (itype < 200) then 
           lrelhum=.true.
           call read_loop_copy_values (olevs,dtype,'MASS','B2P')
         else
           call read_loop_copy_values (olevs,dtype,'WIND','B2P')
         endif
       elseif (trim(dtype) == 'GPSRO') then  
         itype=nint(gpsro_info(index_itype))
         olevs=obs_nlevs
         lat=gpsro_info(index_lat)
         lon=gpsro_info(index_lon)
         gpsro_curve=gpsro_info(index_curve)
         onums=olevs
         call read_loop_copy_values (olevs,dtype,'ALL','B2P')
       else  
         itype=nint(obs_info(index_itype)) 
         olevs=obs_n_channels
         lat=obs_info(index_lat)
         lon=obs_info(index_lon)
         onums=olevs
         call read_loop_copy_values (olevs,dtype,'ALL','B2P')
       endif
!
! Count the number of obs for each separate obs typ
       if (obs_psflag(obs_max_chan_or_levs)) then
         onums1=onums
       else
         onums1=onums+1   ! value augmented by 1 to indicate ps present
       endif
       call count_subtypes (itype,onums1)
!
! For IASI, convert radiances to brightness temperatures here
       if (trim(dtype) == 'IASI') then
         call iasi_radiance_transforms (olevs,obs_data,'R2TB')
       endif
!
! Copy data into arrays for saving
       call sd_copy (obs_max_chan_or_levs,obs_max_fields,dtype,i_file, &
                     itype,olevs,nobs,lat,lon,obs_ps,obs_data,obs_levs)
!
     endif   ! Test on leof
!
! End loop over buffer read
     enddo                 ! loop over do while ( ireadsb(luin) .eq. 0 )
     n_mesg=n_mesg+1       ! count messages read in original file
   enddo                   ! loop over (ireadmg(luin,subset,idate).eq. 0)
!
   if (obs_file_type == 'BUFR') then
     call closbf(luin)     ! close bufr file   
   else
     close (luin)          ! close text file
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
   end subroutine read_loop_do
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_copy_values (olevs,dtype,stype,func)
!
   use m_conv_names, only : bbp, bbt, bbq, bbu, bbv, bbc, bbzq
   use m_bufr_conv, only : conv_values
!
   use m_gpsro_names, only : bbimpp, bbbang
   use m_bufr_gpsro, only : gpsro_values
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
     obs_ps=1.e11
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
             obs_ps=obs_levs(n)                      ! ps
             obs_psflag(obs_max_chan_or_levs)=.true. ! means a valid value found
           endif             
         endif
       enddo
!
     elseif (trim(dtype) == 'GPSRO') then
       do n=1,olevs
         obs_data(n,1)=gpsro_values(bbbang,n)
         obs_levs(n)=gpsro_values(bbimpp,n)-gpsro_curve
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
   end subroutine read_loop_copy_values 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_loop_clean (dtype)
!
   use m_rad_obs_arrays, only : rad_obs_arrays_clean
   use m_iasi_radiance, only : iasi_radiance_cleanup 
!
   character(len=*), intent(in) :: dtype
   integer :: ierr
!
   if (lrad) then 
     call rad_obs_arrays_clean
     if (trim(dtype) == 'IASI') then 
       call iasi_radiance_cleanup (ierr)
     endif
   endif

!
   end subroutine read_loop_clean
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_read_loop
