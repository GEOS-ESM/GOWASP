   module m_obs_list
!
! Module to setup and read obs list (obs header) information, including
! file header information and allowing for multiple rad types to be 
! stored simulatneously.
!
! The public variables in this module contain all the informaion read
! from or written to the obs_list files that are required by the various 
! programs. (The obs-list files are those that contain obs header 
! information and some field information used to select thinned data.) 
! They also contain data lengths and counts. Copies of these variables 
! are held on each processor.
!
! This module was originally developed to deal with a single data type,
! with its own data structures, counts, etc. It was later modified to
! treat multiple data types in the same application. This was
! accomplished by using a data type variable that contains separate
! substrings for all the variables in the original module, with the one
! exception being the list of real values in the observation headers
! that is instead stored in a proceesor-distributed array, Integer and
! character header variables are not distibuted in that way because the
! shared memory modules only distribute real variables. Copies of the
! type data exist on all processors.
!
! When data for any type need to be used, the data for that type stored
! as the type variable are copied to the corresponding original
! variables. In this way, none of the other modules needed to be changed
! to reference the particular desired type.  This is only a tiny amount of
! information and the change of requested type is infrequent.
!
! Initial Code  Ronald Errico  July 15 2014
! Obs_list_type variable added by R. Errico Aug. 21 2014
!
   use m_kinds, only : rkind1, rkind2
   implicit none
!
   private 
!
! public subroutines
   public :: obs_list_setup
   public :: obs_list_read_header
   public :: obs_list_read_recs
   public :: obs_list_print_info 
   public :: obs_list_clean
   public :: obs_list_header_allocate 
   public :: obs_list_types_allocate
   public :: obs_list_types_clean
   public :: obs_list_types_put_setup
   public :: obs_list_types_get_setup
   public :: obs_list_types_put_recs
   public :: obs_list_types_get_recs
   public :: obs_list_types_rec_ids

! public variables
   integer, public :: obs_list_len_i, obs_list_len_r, obs_list_len_c
   integer, public :: obs_list_extra_recs
   integer, public :: obs_list_extra_dim1
   integer, public :: obs_list_dim2
   integer, public :: obs_list_id_time  ! index of list_r array holding time
   integer, public :: obs_list_id_lat   ! index of list_r array holding lat
   integer, public :: obs_list_id_lon   ! index of list_r array holding lon 
   integer, public :: obs_list_format_header  ! indicator for header format
   integer, public :: obs_list_format_recs    
   integer, public :: obs_list_header_nrecs 
   integer, public :: obs_list_tslots         ! number of time slots 
   integer, public :: obs_list_n_channels     ! number of instrument channels
   integer, public :: obs_list_subtypes      
   integer, public :: obs_list_unit_in
   integer, allocatable, public :: obs_list_tslots_recs(:)   
   integer, allocatable, public :: obs_list_i(:,:)  ! integer variables in obs header   
   integer, allocatable, public :: obs_list_counter(:,:)   
   real(rkind1), public :: obs_list_time_delta  ! time period of each time slot
   real(rkind1), public :: obs_list_time_first  
   real(rkind1), allocatable, public :: obs_list_r(:,:)       ! real variables in obs header   
   real(rkind2), allocatable, public :: obs_list_extra_info(:,:) 
   character(len=16), allocatable, public  :: obs_list_c(:,:) ! char variables in obs header   
   character(len=16), allocatable, public  :: obs_list_names(:) ! names of variables in obs header
   character(len=16), allocatable, public  :: obs_list_extra_names(:) 
   character(len=240), allocatable, public :: obs_list_header_crecs(:)
!
! private variables
   character(len=*), parameter :: my_name='m_obs_list'
!
   type obs_list_type
     integer :: ol_len_i, ol_len_r, ol_len_c   
     integer :: ol_extra_recs
     integer :: ol_extra_dim1
     integer :: ol_dim2
     integer :: ol_id_time
     integer :: ol_id_lat
     integer :: ol_id_lon
     integer :: ol_format_header
     integer :: ol_format_recs
     integer :: ol_header_nrecs 
     integer :: ol_tslots
     integer :: ol_n_channels
     integer :: ol_subtypes
     integer :: ol_unit_in
     integer, allocatable :: ol_tslots_recs(:)   
     integer, allocatable :: ol_i(:,:)
     integer, allocatable :: ol_counter(:,:)   
     real(rkind1) :: ol_time_delta
     real(rkind1) :: ol_time_first
     real(rkind2), allocatable :: ol_extra_info(:,:) 
     character(len=16), allocatable :: ol_c(:,:)
     character(len=16), allocatable :: ol_names(:)
     character(len=16), allocatable :: ol_extra_names(:) 
     character(len=240), allocatable :: ol_header_crecs(:)
   end type obs_list_type
!
   type(obs_list_type), allocatable :: obs_list_types(:)
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine obs_list_setup (obs_list_file,lprint,iers)
! 
! Create variables and arrays required for storing observation header 
! information for an observation type based on info in the file header.
! (a file produced by create_rad_obs_list)
!
   logical, intent(in) :: lprint
   integer, intent(out) :: iers
   character(len=*), intent(in) :: obs_list_file
!
   integer :: ier
   integer :: nsum
   integer :: len_name
   character(len=*), parameter :: my_sub=my_name//'::obs_list_setup'
   character(len=1), parameter :: cblank1=' '
   character(len=240) :: cname1
!
! The following surpresses printing of all names in searched list 
! when lprint=.false. and the desired name is not found. this printing is
! always surpressed when 'CLAT or 'CLON' are the desired names, because 
! when not found, the names 'CLATH' and 'CLONH" are sought instead. 
   if (lprint) then   
     len_name=len(my_sub)
     cname1(1:len_name)=my_sub(1:len_name)
   else
     cname1(1:1)=cblank1
   endif
!
   call obs_list_read_header (obs_list_file,ier)
   nsum=obs_list_len_i+obs_list_len_r+obs_list_len_c
   iers=ier
!
! Allocate variables that will hold all obs info for 1 time slot
   obs_list_dim2=maxval(obs_list_tslots_recs(:))
   allocate (obs_list_i(obs_list_len_i,obs_list_dim2))   
   allocate (obs_list_r(obs_list_len_r,obs_list_dim2))   
   allocate (obs_list_c(obs_list_len_c,obs_list_dim2))   
!
! Set variables specifying ids for required header variables
   call find_name (nsum,obs_list_names,.false.,cblank1,  &
                   'CLAT',obs_list_id_lat)
   if (obs_list_id_lat == 0) then
     call find_name (nsum,obs_list_names,.false.,cname1, &
                   'CLATH',obs_list_id_lat)
   endif
!
   call find_name (nsum,obs_list_names,.false.,cblank1,  &
                   'CLON',obs_list_id_lon)
   if (obs_list_id_lon == 0) then
     call find_name (nsum,obs_list_names,.false.,cname1, &
                   'CLONH',obs_list_id_lon)
   endif
!
   call find_name (nsum,obs_list_names,.false.,cname1,   &
                   'DHOURS',obs_list_id_time)
!
   obs_list_id_lat=obs_list_id_lat-obs_list_len_i
   obs_list_id_lon=obs_list_id_lon-obs_list_len_i
   obs_list_id_time=obs_list_id_time-obs_list_len_i
   if (obs_list_id_lat <= 0 .or. obs_list_id_lon <= 0 .or. &
       obs_list_id_time <= 0) then
     iers=iers+1
     if (lprint) then
       print *,' '
       print *,'Error in routine ',my_sub
       print ('(a,3i4)'), &
              'obs_list_id_lat,obs_list_id_lon,obs_list_id_time = ', &
               obs_list_id_lat,obs_list_id_lon,obs_list_id_time 
     endif
   endif
!
   end subroutine obs_list_setup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine obs_list_read_header (obs_list_file,ier)
! 
! Read file header information from a file created by create_rad_obs_list
! (This file contains a list of possibly thinned observation meta data 
! along with some NR field values at the obs locations.)
!
   integer, intent(out) :: ier
   character(len=*), intent(in) :: obs_list_file
!
   integer :: n, nsum
   integer :: ierr
   integer :: iunit
   character(len=*), parameter :: my_sub=my_name//'::obs_list_read_header'
!  
   ier=0 
   iunit=obs_list_unit_in 
   open (iunit,file=trim(obs_list_file), form='unformatted')
!
! Read obs_list_scalars
   read (iunit) obs_list_format_header, obs_list_format_recs
   read (iunit) obs_list_n_channels, obs_list_header_nrecs
   read (iunit) obs_list_extra_dim1,obs_list_extra_recs
   read (iunit) obs_list_len_i, obs_list_len_r, obs_list_len_c, &
                obs_list_tslots, obs_list_time_delta,           & 
                obs_list_time_first, obs_list_subtypes
!
! Allocate obs_list_ arrays to be read from header 
   nsum=obs_list_len_i+obs_list_len_r+obs_list_len_c
   call obs_list_header_allocate (nsum,ierr) 
!
! Read obs_list_ arrays
   do n=1,obs_list_header_nrecs
     read (iunit) obs_list_header_crecs(n)
   enddo
!
   read (iunit) obs_list_extra_names(1:obs_list_extra_recs)
   do n=1,obs_list_extra_recs
     read (iunit) obs_list_extra_info(:,n)
   enddo
!
   read (iunit) obs_list_tslots_recs,obs_list_counter
   read (iunit) obs_list_names(1:nsum)
!
   end subroutine obs_list_read_header
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine obs_list_header_allocate (nsum,ier) 
!
! Allocate arrays need to store information in the header of a file created 
! by the program create_rad_obs_list.
!
   implicit none
   integer, intent(in)  :: nsum
   integer, intent(out) :: ier
!
   allocate (obs_list_header_crecs(obs_list_header_nrecs))
   allocate (obs_list_extra_names(obs_list_extra_recs))
   allocate (obs_list_extra_info(obs_list_extra_dim1,obs_list_extra_recs))
   allocate (obs_list_tslots_recs(obs_list_tslots-1))
   allocate (obs_list_names(nsum))
   allocate (obs_list_counter(obs_list_tslots,obs_list_subtypes))
   ier=0
!
   end subroutine obs_list_header_allocate
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine obs_list_read_recs (k_tslot)
!
!  Read records of data for all observations of a particular obs type and 
!  time slot from a file created by the program create_rad_obs_list.
! 
   integer, intent(in) :: k_tslot
!
   integer :: n
   integer :: iunit
   character(len=*), parameter :: my_name='obs_list_read_recs'
!  
   iunit=obs_list_unit_in
   do n=1,obs_list_tslots_recs(k_tslot)
     read (iunit) obs_list_i(:,n),obs_list_r(:,n),obs_list_c(:,n)
   enddo
!
   if (k_tslot == obs_list_tslots) then 
     close (iunit)
   endif
!
   end subroutine obs_list_read_recs
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine obs_list_print_info (call_loc)
! 
!  Print some information read from the header of a file created by the 
!  program create_rad_obs_list.
!
   implicit none
   character(len=*), intent(in) :: call_loc  ! calling location
!
   integer :: n,num
!
   print ('(2a)'),'Obs list info: called from ',trim(call_loc)
   do n=1,obs_list_header_nrecs
     print ('(a,i2,2a)'),'Header ',n,': ',trim(obs_list_header_crecs(n))
   enddo
   print ('(a,4i4)'),'obs_list_len_i,r,c= ',obs_list_len_i, &
        obs_list_len_r,obs_list_len_c
   print ('(a,i4,2f8.4)'),'obs_list_tslots,time_delta,time_first=', &
        obs_list_tslots,obs_list_time_delta,obs_list_time_first
   print ('(a,25i8)'),'obs_list_tslots_recs= ',obs_list_tslots_recs(:)
   num=obs_list_len_i+obs_list_len_r+obs_list_len_c
   print ('(a,30(1x,a))'),'obs_list_names=', &
        (trim(obs_list_names(n)),n=1,num)
!
   end subroutine obs_list_print_info
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_clean
!
!  Deallocate arrays used to pass information between this module and
!  other mudules for a single obs data type
!
   deallocate (obs_list_i,obs_list_r,obs_list_c,obs_list_tslots_recs)
   deallocate (obs_list_header_crecs, obs_list_counter, obs_list_names)
   deallocate (obs_list_extra_info, obs_list_extra_names)
   end subroutine obs_list_clean
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_allocate (ntypes)
!
! Allocate array to contain type-dependent obs information local to the module 
! m_obs_list  
!
   implicit none
   integer, intent(in) :: ntypes
   allocate (obs_list_types(ntypes))
   end subroutine obs_list_types_allocate
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_clean
!
! The reverse of obs_list_types_allocate
!
   deallocate (obs_list_types)
   end subroutine obs_list_types_clean
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_put_setup (ntype)
!
!  Copy obs_list_file header information values from a private type-
!  dependent array to a public array to be used by other modules.  
!  Each call to this routine concerns a particular obs type and time slot.
!  Also, allocate some arrays based on specific size requirement for 
!  this obs type and time slot.  
!
   implicit none
!   
   integer, intent(in) :: ntype
   integer :: nsum
!
   obs_list_types(ntype)%ol_len_i=obs_list_len_i
   obs_list_types(ntype)%ol_len_r=obs_list_len_r
   obs_list_types(ntype)%ol_len_c=obs_list_len_c
   obs_list_types(ntype)%ol_extra_recs=obs_list_extra_recs
   obs_list_types(ntype)%ol_extra_dim1=obs_list_extra_dim1
   obs_list_types(ntype)%ol_dim2=obs_list_dim2
   obs_list_types(ntype)%ol_id_time=obs_list_id_time
   obs_list_types(ntype)%ol_id_lat=obs_list_id_lat
   obs_list_types(ntype)%ol_id_lon=obs_list_id_lon
   obs_list_types(ntype)%ol_format_header=obs_list_format_header
   obs_list_types(ntype)%ol_format_recs=obs_list_format_recs
   obs_list_types(ntype)%ol_header_nrecs=obs_list_header_nrecs 
   obs_list_types(ntype)%ol_tslots=obs_list_tslots
   obs_list_types(ntype)%ol_n_channels=obs_list_n_channels
   obs_list_types(ntype)%ol_subtypes=obs_list_subtypes
   obs_list_types(ntype)%ol_unit_in=obs_list_unit_in
   obs_list_types(ntype)%ol_time_delta=obs_list_time_delta
   obs_list_types(ntype)%ol_time_first=obs_list_time_first
!
   allocate (obs_list_types(ntype)%ol_header_crecs(obs_list_header_nrecs))
   allocate (obs_list_types(ntype)%ol_extra_names(obs_list_extra_recs))
   allocate (obs_list_types(ntype)%ol_tslots_recs(obs_list_tslots-1))
   nsum=obs_list_len_i+obs_list_len_r+obs_list_len_c
   allocate (obs_list_types(ntype)%ol_names(nsum))
   allocate (obs_list_types(ntype)%ol_counter(obs_list_tslots, &
                                              obs_list_subtypes))
   allocate (obs_list_types(ntype)%ol_extra_info(obs_list_extra_dim1, &
                                                 obs_list_extra_recs))
   allocate (obs_list_types(ntype)%ol_i(obs_list_len_i,obs_list_dim2))   
   allocate (obs_list_types(ntype)%ol_c(obs_list_len_c,obs_list_dim2))   
!
   obs_list_types(ntype)%ol_extra_info(:,:)=obs_list_extra_info(:,:)   
   obs_list_types(ntype)%ol_header_crecs(:)=obs_list_header_crecs(:) 
   obs_list_types(ntype)%ol_extra_names(:)=obs_list_extra_names(:)     
   obs_list_types(ntype)%ol_tslots_recs(:)=obs_list_tslots_recs(:)   
   obs_list_types(ntype)%ol_names(:)=obs_list_names(:)   
   obs_list_types(ntype)%ol_counter(:,:)=obs_list_counter(:,:) 
!
   end subroutine obs_list_types_put_setup 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_put_recs (ntype)
!
!  Copy all obs_list_file record arrays from a private type-dependent 
!  array to a public array to be used by other modules for a particular 
!  obs type and time slot.
!
   implicit none
!   
   integer, intent(in) :: ntype
!
   obs_list_types(ntype)%ol_i(:,:)=obs_list_i(:,:) 
   obs_list_types(ntype)%ol_c(:,:)=obs_list_c(:,:) 
   obs_list_types(ntype)%ol_counter(:,:)=obs_list_counter(:,:)   
!
   end subroutine obs_list_types_put_recs 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_get_setup (ntype)
!
!  The reverse of obs_list_types_put_setup
!
   implicit none
!   
   integer, intent(in) :: ntype
   integer :: ierr 
   integer :: nsum
!
   obs_list_len_i=obs_list_types(ntype)%ol_len_i
   obs_list_len_r=obs_list_types(ntype)%ol_len_r
   obs_list_len_c=obs_list_types(ntype)%ol_len_c
   obs_list_extra_recs=obs_list_types(ntype)%ol_extra_recs
   obs_list_extra_dim1=obs_list_types(ntype)%ol_extra_dim1
   obs_list_dim2=obs_list_types(ntype)%ol_dim2
   obs_list_id_time=obs_list_types(ntype)%ol_id_time
   obs_list_id_lat=obs_list_types(ntype)%ol_id_lat
   obs_list_id_lon=obs_list_types(ntype)%ol_id_lon
   obs_list_format_header=obs_list_types(ntype)%ol_format_header
   obs_list_format_recs=obs_list_types(ntype)%ol_format_recs
   obs_list_header_nrecs=obs_list_types(ntype)%ol_header_nrecs
   obs_list_tslots=obs_list_types(ntype)%ol_tslots
   obs_list_n_channels=obs_list_types(ntype)%ol_n_channels
   obs_list_subtypes=obs_list_types(ntype)%ol_subtypes
   obs_list_unit_in=obs_list_types(ntype)%ol_unit_in
   obs_list_time_delta=obs_list_types(ntype)%ol_time_delta
   obs_list_time_first=obs_list_types(ntype)%ol_time_first
!
   allocate (obs_list_i(obs_list_len_i,obs_list_dim2))   
   allocate (obs_list_r(obs_list_len_r,obs_list_dim2))   
   allocate (obs_list_c(obs_list_len_c,obs_list_dim2))   
!
   nsum=obs_list_len_i+obs_list_len_r+obs_list_len_c
   call obs_list_header_allocate (nsum,ierr) 
!
   obs_list_extra_info(:,:)=obs_list_types(ntype)%ol_extra_info(:,:)
   obs_list_header_crecs(:)=obs_list_types(ntype)%ol_header_crecs(:)
   obs_list_extra_names(:)=obs_list_types(ntype)%ol_extra_names(:)
   obs_list_tslots_recs(:)=obs_list_types(ntype)%ol_tslots_recs(:)
   obs_list_names(:)=obs_list_types(ntype)%ol_names(:)
   obs_list_counter(:,:)=obs_list_types(ntype)%ol_counter(:,:)
!
   end subroutine obs_list_types_get_setup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_get_recs (ntype)
!
!  The reverse of obs_list_types_put_recs
!
   implicit none
!   
   integer, intent(in) :: ntype
!
   obs_list_i(:,:)=obs_list_types(ntype)%ol_i(:,:)
   obs_list_c(:,:)=obs_list_types(ntype)%ol_c(:,:)
   obs_list_counter(:,:)=obs_list_types(ntype)%ol_counter(:,:)
!
   end subroutine obs_list_types_get_recs 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   subroutine obs_list_types_rec_ids (ntypes,rec_ids_dim3,ntypes_max, &
                                      rec_r_dim1,rec_ids_max,rec_ids)
!
! Determine variables to be used as dimensions and indexes for the arrays 
! rec_list_*.  These are required to accomodate use of different obs types 
! having different obs meta info and their array requirements within a
! single program execution. 
!
   implicit none
!   
   integer, intent(in) :: ntypes
   integer, intent(in) :: rec_ids_dim3
   integer, intent(in) :: ntypes_max
   integer, intent(out) :: rec_r_dim1
   integer, intent(out) :: rec_ids_max
   integer, intent(out) :: rec_ids(2,ntypes_max,rec_ids_dim3)
!
   integer ntime,ntype
!
! The values of obs_list_types(n)%ol_tslots-1 are assumed the same for all n
!
   do ntime=1,obs_list_types(1)%ol_tslots-1
     rec_ids(1,1,ntime)=1
     rec_ids(2,1,ntime)=obs_list_types(1)%ol_tslots_recs(ntime)
     do ntype=2,ntypes
       rec_ids(1,ntype,ntime)=rec_ids(2,ntype-1,ntime)+1
       rec_ids(2,ntype,ntime)=rec_ids(1,ntype,ntime)-1+ &
                              obs_list_types(ntype)%ol_tslots_recs(ntime)
     enddo
   enddo
!
   rec_ids_max=maxval(rec_ids(2,ntypes,1:obs_list_types(1)%ol_tslots-1)) 
!
   rec_r_dim1=obs_list_types(1)%ol_len_r
   do ntype=2,ntypes
     rec_r_dim1=max(rec_r_dim1,obs_list_types(ntype)%ol_len_r)
   enddo
!
   end subroutine obs_list_types_rec_ids 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
! 
   end module m_obs_list
