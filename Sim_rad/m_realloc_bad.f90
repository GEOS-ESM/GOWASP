   module m_realloc_bad
!
! Get text file of thinned and QC'd obs locations and list of channels 
! used at each location to set unused channels to a bad value so that
! the channel will not be used.
!
   use m_kinds, only : rkind1
!
   private
   public :: realloc_bad_read
   public :: realloc_bad_set 
   public :: realloc_bad_counter
!
   integer, parameter :: ndmeta1=3
   integer :: said_copy
   integer :: nlocs
   integer :: nchans
   integer :: no_match
   integer :: glist_lat_id, glist_lon_id
   logical, allocatable :: lchan(:,:)
   real(rkind1) :: dmin  ! precision for matching lats and lons
   real(rkind1), parameter :: bad_factor=1.15 ! factor to inflate bad values
   real(rkind1), allocatable :: dmeta(:,:)
!
! bad_factor must be large enough that any so inflated obs values will be 
! recognized as bad so that it is discarded by the das QC but not be too 
! large as to make the complete set of channels at an indvidual location 
! suspect to the point the entire location is discarded.  
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
!
   subroutine realloc_bad_read (dtype,said,ctest,lprint)
!
!  Read text file made from filtered (i.e., list of used) ods file data 
!  for reqested data type and satellite id.  Also determine some 
!  indexes for pulling lat and lon from meta data array.
!
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c
   use m_obs_list, only : obs_list_names
! 
   implicit none
! 
   logical, intent(in) :: lprint
   integer, intent(in) :: said
   character (len=*), intent(in) :: dtype
   character (len=*), intent(in) :: ctest
!
   logical, parameter :: lstop=.false.
   logical :: said_found
   integer :: nsat, nsats, n, ns
   integer :: n1, m1, m2
   integer :: said_file
   integer, parameter :: txt_unit_in=51
   integer, parameter :: nsats_max=1000
   character(len=100) :: txt_file_in
   character(len=1) :: cdum
!
   if (allocated(lchan)) then
     deallocate (lchan)
   endif
   if (allocated(dmeta)) then
     deallocate (dmeta)
   endif
!
   said_copy=said    
!
   txt_file_in='radlocs_'//trim(dtype)//'.txt'
   open (txt_unit_in,file=trim(txt_file_in))
   if (lprint) then 
     print ('(3a,i3)'),' text input file=',trim(txt_file_in),   &
                ' opened on unit=',txt_unit_in
   endif 
!
! Find set of obs for requested satellite ID. 
! For AIRS and AMSUAAQUA, SAID is not available in the profile data,
! but then nsats=1. For nsats=1, the assumption is that the proper set
! of obs has been found.
! 
   said_found=.false.
   nsats=1
   do nsat=1,nsats_max
     if (said_found .or. nsat > nsats) then
       exit
     else
       read (txt_unit_in,'(5i8)') ns,said_file,nlocs,nsats,nchans
       if (said == said_file .or. nsats == 1) then
         said_found=.true.
       else    ! skip to next said on file
         do n=1,nlocs
           read (txt_unit_in,'(a1)') cdum
         enddo
       endif
     endif
   enddo
!
   no_match=0  ! initialize counter
!
! If subtype found in list pulled from .ods file, then read meta 
! data and channel usage for all locations for this subtype       
   if (said_found) then  
     allocate (lchan(nchans,nlocs))
     allocate (dmeta(ndmeta1,nlocs))
     do n=1,nlocs
       read (txt_unit_in,'(3f7.2,1x,320l1)') dmeta(:,n),lchan(1:nchans,n)
     enddo
!
     if (ctest == 'T' .and. lprint) then
       print ('(a,i4,i8)'),'requested said its nlocs on file:',said,nlocs
       print ('(a,2f8.2)'),'dmeta lon lat for first location:',dmeta(1:2,1)
       print ('(a)'),'channels used at first location:'
       print ('(80l1)'),lchan(1:nchans,1)
     endif
!
   elseif (lprint) then 
     print *,'Requested SAID not found in m_realloc_bad: SAID=',said
   endif
   close (txt_unit_in)
!
! Set precision for matching lat and lons between BUFR and ODS files
! If the precision in the BUFR file is low, inflate this value a bit
   dmin=0.02
!
   m1=obs_list_len_i+1
   m2=obs_list_len_i+obs_list_len_r
   n1=obs_list_len_r
   call find_name (n1,obs_list_names(m1:m2),lstop,' ', &
                   'CLAT',glist_lat_id)
   if (glist_lat_id == 0) then
     call find_name (n1,obs_list_names(m1:m2),lstop,' ', &
                   'CLATH',glist_lat_id)
   else
     dmin=.05  ! since BUFR precision is low
   endif
!
   call find_name (n1,obs_list_names(m1:m2),lstop,' ', &
                   'CLON',glist_lon_id)
   if (glist_lon_id == 0) then
     call find_name (n1,obs_list_names(m1:m2),lstop,' ', &
                   'CLONH',glist_lon_id)
   else
     dmin=.05 ! since BUFR precision is low
   endif
! 
   end subroutine realloc_bad_read
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine realloc_bad_set (dtype)
!
! Pull list of cchannels used at the current bs location and distort 
! obs value if channel not used there. 
!
   use m_rad_obs_arrays, only : obs_info,obs_values
!
   implicit none
!
   character(len=*), intent(in) :: dtype
!
   logical :: lremove
   integer :: n, nfound
   real(rkind1) :: dlon,dlat,xlon   
!   
   xlon=mod(obs_info(glist_lon_id),360.)
   if (xlon < 0.) then
     xlon=xlon+360.
   endif
!
   nfound=0  ! means not found
   do n=1,nlocs
     dlon=xlon-dmeta(1,n)
     dlat=obs_info(glist_lat_id)-dmeta(2,n)
     if (abs(dlon) < dmin .and. abs(dlat) < dmin) then 
       nfound=n
       exit
     endif 
   enddo
!
! Increase obs values of unused channels by bad_factor
! This is done by a mutiplicative factor rather than addition because
! most values are brightness temperatures but for IASI they are scaled 
! radiances. 
   do n=1,nchans
     if (nfound > 0) then
       lremove= .not. lchan(n,nfound)   ! remove this channel
     else
       lremove=.true.  ! remove all channels since obs loc not in accepted list
     endif
!
     if (lremove) then  ! spoil the value for this channel
       if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then 
         obs_values(n,2)=4  ! set QC flag to non zero value
       else
         obs_values(n,1)=obs_values(n,1)*bad_factor
       endif
     endif
   enddo
!
   if (nfound == 0) then
     no_match=no_match+1  ! increment counter
   endif
!
   end subroutine realloc_bad_set
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine realloc_bad_counter (dtype)
!
! Print counter if non zero
!
   implicit none
   character(len=*), intent(in) :: dtype
!
   if (no_match > 0) then
     print *,'Number of locations not matching those QC-accepted in .ods file:'
     print ('(3a,i4,a,i6)'),'dtype=',trim(dtype),'   sat id=',said_copy, &
                    '   number not matching=',no_match
   endif
!
   end subroutine realloc_bad_counter 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_realloc_bad
