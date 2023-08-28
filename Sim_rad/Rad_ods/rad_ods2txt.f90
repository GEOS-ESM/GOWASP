   program satwind_ods2txt 
!
! Input ods file.  Output text file of satwind meta data
! (lat, lon, said) and logical flags indicating what channels 
! were used by the DAS at each location.  Only include locations 
! where at least one channel was used (i.e., after thinning and QC).
!  
   use m_ods
   use m_sat_info_table, only : sat_info_table_read
   use m_sat_info_table, only : sat_info_table_get_info
!
   implicit none
!
   logical :: ods_exists
!
   integer, parameter :: ndmeta1=3
   integer, parameter :: txt_unit_out=21
   integer, parameter :: nsats_max=10 ! max number of satellites for any instr.
   integer :: nlocs_max
   integer :: nchans_max
   integer :: kx, nchans
   integer :: n, id, ichan
   integer :: ndate, ntime
   integer :: nobs
   integer :: nlocs
   integer :: argc
   integer(4) :: iargc
   integer :: ierr
   integer :: nsat, nsats
   integer :: said(nsats_max)
!
   real(4) :: xlon,xlat
!
   character(len=14)  :: cdtime0
   character(len=100) :: txt_file_out     ! prepbufr file of satwind obs
   character(len=120) :: ods_path
   character(len=200) :: ods_file
   character(len=30)  :: ods_expid
   character(len=240) :: sat_info_file
   character(len=30)  :: cdum_type  ! output from Get_ODS 
   character(len=20)  :: dtype
   character(len=20)  :: inst_sat(nsats_max)
   character(len=1)   :: ctest
!
   logical, allocatable :: lchan(:,:)
   real(4), allocatable :: dmeta(:,:)
!
   type (ods_vect) :: ods
!
! Read arguments defined when program is executed in script
   argc = iargc()
   if (argc /= 6) then
     print *,'Number of arguments /= 6'
     stop
   endif
!
   call GetArg( 1_4, cdtime0)    ! to get date and time of data on ods file
   call GetArg( 2_4, sat_info_file)
   call GetArg( 3_4, ods_path)   ! path for ods file (before /Yxxxx/Mxx/etc)
   call GetArg( 4_4, ods_expid)  ! exp id as part of ods file name 
   call GetArg( 5_4, dtype)      ! data type
   call GetArg( 6_4, ctest)      ! T or F if sample to be printed for testing
!
   print *,'Setting real locations for dtype=',trim(dtype)
!
! set a maximum number of thinned and QC'd locations per satellite
! (used to define array dimension) 
   if (trim(dtype) == 'ATMS' .or. trim(dtype) == 'SSMIS' .or. &
       trim(dtype) == 'AMSR2' ) then 
     nlocs_max=1000000
   else 
     nlocs_max=14000
   endif
!
! set a maxinum channel id as the largest of the channels used.
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'IASI' ) then 
     nchans_max=230
   else if (trim(dtype) == 'CRIS' ) then
     nchans_max=320
   else if (trim(dtype) == 'CRISFSR' ) then
     nchans_max=450
   else
     nchans_max=999  ! will be overwritten by actual channel numbers
   endif
!
! Get list of satellite ids for obs of requested data type
   call sat_info_table_read (sat_info_file,.true.)   
   call sat_info_table_get_info (nsats_max,dtype,nsats,nchans, &
         inst_sat,said,ierr)
!
   nchans=min(nchans_max,nchans) 
   allocate (lchan(nchans,nlocs_max))
   allocate (dmeta(ndmeta1,nlocs_max))
   print *,'nsats_max, nlocs_max, nchans set to ', &
            nsats_max, nlocs_max, nchans
!
! Open output file
   txt_file_out='radlocs_'//trim(dtype)//'.txt'
   open (txt_unit_out,file=trim(txt_file_out))
   print ('(3a,i3)'),' text output file=',trim(txt_file_out),   &
                ' opened on unit=',txt_unit_out
!
   read (cdtime0(1:8),'(i8)')  ndate                                         
   read (cdtime0(9:14),'(i6)') ntime 
!
! Read ods files to construct list of locations and channels used
   do nsat=1,nsats
!   
     ods_file=trim(ods_path)//'/Y'//cdtime0(1:4)//'/M'       &
              //cdtime0(5:6)//'/'//trim(ods_expid)//'.diag_' &
              //trim(inst_sat(nsat))//'.'//cdtime0(1:8)//'_' &
              //cdtime0(9:10)//'z.ods'
     print *,'ODS_file=',trim(ods_file)
!
     call ODS_Get (trim(ods_file),ndate,ntime,cdum_type,ods,ierr)
     nlocs=0
     if (ierr /= 0) then 
       print *,'ods file does not exist: ierr=',ierr
       ods_exists=.false.
     else
       ods_exists=.true.
       nobs=ods%data%nobs  ! number of obs in this file  
       print *,'n obs in ods file =',nobs
!
! loop over all obs in ods file       
       lchan(:,:)=.false.
       do n=1,nobs  
         if (ods%data%qcexcl(n) == 0) then  ! check if obs actually used 
!
           xlon=mod(ods%data%lon(n),360.)
           if (xlon<0.) xlon=xlon+360.
           xlat=ods%data%lat(n)
!
! Determine if obs n is at a location not already found
           call search_obs (ndmeta1,nlocs_max,nlocs,xlat,xlon,dmeta,id)
!
           ichan=nint(ods%data%lev(n))
           if (ichan <= nchans) then 
             lchan(ichan,id)=.true.    ! flag channel at this locat as used   
           endif
           if (id > nlocs .and. id <= nlocs_max) then ! save as new location  
             dmeta(1,id)=xlon
             dmeta(2,id)=xlat
             dmeta(3,id)=ods%data%time(n)
             nlocs=id
           endif
!
         endif
       enddo 
!    
     endif 
!
! Write meta data and channel list used for current satellite id
     write (txt_unit_out,'(5i8)') nsat,said(nsat),nlocs,nsats,nchans
     do n=1,nlocs
       write (txt_unit_out,'(3f7.2,1x,320l1)') dmeta(:,n),lchan(1:nchans,n)
     enddo
!
     if (ctest == 'T') then 
       print ('(a,5i8)'),'nsat,said(nsat),nlocs,nsats,nchans = ', &
                          nsat,said(nsat),nlocs,nsats,nchans
       print ('(a,3f8.2)'),'lon, lat, time for first location:',dmeta(:,1)
       print ('(a)'),'channels used at this location:'
       print ('(80l1)'),lchan(1:nchans,1)
     endif
!
     print *,' '
     print *,'Output text data written for said, nlocs=',said(nsat),nlocs   
     call ods_clean (ods,ierr)
   enddo  ! loop over sats 
!
   close (txt_unit_out)   
   print *,' '
   print *,'program completed' 
!
   end program satwind_ods2txt 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine search_obs (ndmeta1,ndmeta2,nlocs,xlat,xlon,dmeta,id)
!
! 
!
   integer, intent(in) :: ndmeta1,ndmeta2,nlocs
   real(4), intent(in) :: xlat,xlon,dmeta(ndmeta1,ndmeta2)
   integer, intent(out) :: id
!
   integer :: n
!
! default values are that this is a new location
   id=nlocs+1             
!
! see if this location was previously found
   do n=1,nlocs 
     if (abs(xlon-dmeta(1,n)) <= .01 .and. abs(xlat-dmeta(2,n)) <= .01) then 
       id=n
       exit
     endif
   enddo
!       
   end subroutine search_obs 
