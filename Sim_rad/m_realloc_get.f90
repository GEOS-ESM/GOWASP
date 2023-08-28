   module m_realloc_get
!
! Get text file of thinned and QC'd obs locations obtained from 
! the ods files and feed them to the create_rad_obs_list program
!
   use m_kinds, only : rkind1
!
   implicit none
!   
   private
   public :: realloc_get_read
   public :: realloc_get_find
!
   integer, parameter :: ndmeta1=2
   integer, parameter :: nsats_max=20
   integer :: nlocs_max
   integer :: nsats
   integer :: glist_lat_id, glist_lon_id, glist_subset_id 
   integer :: isatinfo(2,nsats_max)
   real(rkind1), allocatable :: dmeta(:,:,:)
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
!
   subroutine realloc_get_read (dtype,c_test) 
!
!  Read text file made from filtered (i.e., list of used) ods file data 
!  for all satellite ids for the requested data type. Also determine some 
!  indexes for pulling lat and lon from meta data array.
!
   use m_rad_obs_arrays, only : obs_info_num,obs_info_names
!
   implicit none
!
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: c_test
!
   logical, parameter :: lstop=.false.
   integer :: n, ns, nx
   integer :: nsfile, nlocs
   integer :: said
   integer, parameter :: txt_unit_in=51
   character(len=100) :: txt_file_in
   character(len=1) :: cdum
!
   if (allocated(dmeta)) then
     deallocate (dmeta)
   endif
!    
   txt_file_in='radlocs_'//trim(dtype)//'.txt'
   open (txt_unit_in,file=trim(txt_file_in))
   print ('(3a,i3)'),' text input file=',trim(txt_file_in),   &
                ' opened on unit=',txt_unit_in
!
! Read text file 1st of 2 times to obtain a list of SAID 
! (satellite ID numbers) and corresponding counts of numbers of 
! obs locations for each 
   nsats=0
   nlocs_max=0
   do ns=1,100
     read (txt_unit_in,'(5i8)') nx,said,nlocs,nsfile     
     if (nlocs > 0) then 
       nsats=nsats+1
       do n=1,nlocs
         read (txt_unit_in,'(a1)') cdum
       enddo
       nlocs_max=max(nlocs_max,nlocs)
       isatinfo(1,nsats)=said
       isatinfo(2,nsats)=nlocs
     endif
     if (ns == nsfile) exit
   enddo
!
   if (c_test == 'T') then
     print *,'isatinfo said and nlocs:'
     print ('(10i8)'),isatinfo(1,1:nsats) 
     print ('(10i8)'),isatinfo(2,1:nsats) 
   endif
!
   close (txt_unit_in)

   allocate (dmeta(ndmeta1,nlocs_max,nsats))  
!
! Read text file 2nd of 2 times to obtain, for each separate SAID, 
! the lat and lon of each obs and the list of which channels were used 
! at that location. 
   open (txt_unit_in,file=trim(txt_file_in))
   nsats=0 
   do ns=1,nsats_max
     read (txt_unit_in,'(5i8)') nx,said,nlocs,nsfile     
     if (nlocs > 0) then 
       nsats=nsats+1
       do n=1,nlocs
         read (txt_unit_in,'(2f7.2)') dmeta(:,n,nsats)
       enddo
     endif
     if (ns == nsfile) exit
   enddo
   close (txt_unit_in)
!
   if (c_test == 'T') then
     print *,'dmeta lon lat for first occurance of each satellite'
     print ('(10f8.2)'),dmeta(1,1,1:nsats)
     print ('(10f8.2)'),dmeta(2,1,1:nsats) 
   endif
!
   call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLAT',glist_lat_id)
   if (glist_lat_id == 0) then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLATH',glist_lat_id)
   endif
!
   call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLON',glist_lon_id)
   if (glist_lon_id == 0) then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLONH',glist_lon_id)
   endif
!
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA' ) then
     glist_subset_id=-999
   else
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                     lstop,' ','SAID',glist_subset_id)
   endif
!
   if (c_test == 'T') then
     print ('(a,3i4)'),'m_realloc_get ids:', &
                        glist_lon_id,glist_lat_id,glist_subset_id
   endif

!
   end subroutine realloc_get_read
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine realloc_get_find (lremove)
!
! Check if obs in obs_info and obs_values array is in the ods text file list.
! If not, overwrite obs location to values that will flag current obs as 
! not OK.
!
   use m_rad_obs_arrays, only : obs_info
!
   implicit none
!
   logical, intent(in) :: lremove ! prevent duplicates from being found
!
   logical :: lfound
   integer :: nsat
   integer :: n
   real(rkind1), parameter :: dmin=0.02
   real(rkind1) :: xlon
   real(rkind1) :: dlon,dlat   
!   
   if (glist_subset_id < 0) then
     nsat=1
   else
     nsat=0
     do n=1,nsats
       if (isatinfo(1,n) == obs_info(glist_subset_id) ) then
         nsat=n
       endif
     enddo
   endif
!
   xlon=mod(obs_info(glist_lon_id),360.)
   if (xlon < 0.) then
     xlon=xlon+360.
   endif
!
   lfound=.false.
   if (nsat > 0) then
     do n=1,isatinfo(2,nsat)
       dlon=xlon-dmeta(1,n,nsat)
       dlat=obs_info(glist_lat_id)-dmeta(2,n,nsat)
       if (abs(dlon) < dmin .and. abs(dlat) < dmin) then 
         lfound=.true.
         if (lremove) then ! remove this location from further consideration
           dmeta(1:2,n,nsat)=990.         
         endif
         exit
       endif 
     enddo
   endif
!
! If the location is not one that was used after thinning and QC, 
! then set lat and lon to bad values that will cause the obs to get
! discarded.
   if (.not. lfound) then
     obs_info(glist_lat_id)=999.
     obs_info(glist_lon_id)=999.
   endif
!
   end subroutine realloc_get_find
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_realloc_get
