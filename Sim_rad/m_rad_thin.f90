!
   module m_rad_thin
!
! Module used to store required obs bufr input data header information, and if 
! requested, to also thin the data.
!
! Thining is performed by selecting a single "best"
! observation for each geogaphical "box." It essentially saves 
! the complete header information required to simulate the 
! observation in an array named "box_info" and assigns it a score.  
! Scores for ther subsequently considered obs that fall into the 
! same box are compared, and only the the information for the obs with
! the best score is kept.
! 
! Initial code by Ronald Errico July 2014
! Modified to include non-thiined data Oct 2015     
!
   use m_kinds, only : rkind1,rkind2
!
   use m_rad_obs_arrays, only : obs_info_num,obs_info_names,obs_info
   use m_rad_obs_arrays, only : obs_time_slots,obs_time_delta
   use m_rad_obs_arrays, only : obs_time_first
!
   implicit none      
!
   private
   public :: rad_thin_clean       ! deallocate arrays
   public :: rad_thin_setup       ! initilize module 
   public :: rad_thin_put         ! consider new obs for saving in box
   public :: rad_thin_write       ! write info in boxes to a file
   public :: rad_thin_box_count   ! count numbers of obs in time slots for each subtype
!
   private :: rad_thin_set_boxes 
   private :: rad_thin_check_hdr
   private :: rad_thin_find_subset
   private :: rad_thin_penalty      
   private :: rad_thin_find_box_id
   private :: rad_thin_box_update
!
   public :: thin_info_dim,thin_info_names,thin_info
!
   logical :: l_thin_data   ! true if data are thinned 
!
   integer(8), parameter :: no_thin_max_bytes=450000000  ! max num of bytes to store headers
   integer :: no_thin_id    ! obs index if no thinning is performed
   integer :: thin_info_dim
   integer :: thin_box_dim(3)
   integer :: thin_nlats
   integer :: thin_lat_id, thin_lon_id, thin_year_id, thin_subset_id
   integer :: thin_nobs_id, thin_penalty_id, thin_dhours_id, thin_boxcnt_id
   integer :: thin_fovn_id, thin_qmrkh_id
   integer, allocatable :: thin_box_counter(:,:)
   integer, allocatable :: thin_nlonP(:)
   integer, allocatable :: thin_nlons(:)
   real(rkind1), parameter :: r90=90._rkind1
   real(rkind1), parameter :: r180=180._rkind1
   real(rkind1), parameter :: r360=360._rkind1
   real(rkind1) :: thin_dlat
   real(rkind1) :: thin_box_size
   real(rkind1), allocatable :: thin_info(:) 
   real(rkind1), allocatable :: thin_dlons(:) 
   real(rkind1), allocatable :: thin_box_info(:,:,:)
   character(len=16), allocatable :: thin_info_names(:) 
!
   character(len=*),parameter :: myname='m_rad_thin'
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_setup (dtype,rcfile,cdtime,lstop,lprint,ier) 
!
! Set variables and arrays used for thinning obs based on .rc files
!
   use m_rad_thin_flds, only : rad_thin_flds_setup
   use m_rad_thin_flds, only : fld_thin_box_size
   use m_rad_thin_flds, only : fld_thin_box_subtypes
   use m_rad_thin_flds, only : fld_num_tvar,fld_num_tcon,fld_names
!
   implicit none  
!
   logical, intent(in) :: lstop
   logical, intent(in) :: lprint
   character(len=*), intent(in) :: dtype   
   character(len=*), intent(in) :: rcfile
   character(len=*), intent(in) :: cdtime
   integer, intent(out) :: ier
!
   integer :: ier1, ier2
   integer :: n,m
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_setup'
!
! Read resource file that describes thinning parameters
   call rad_thin_flds_setup (dtype,rcfile,cdtime,lstop,lprint,ier)
   if (ier /= 0) then
     print *,' ' 
     print ('(a,i5,2a)'),'Error flag =',ier, & 
            ' returned from call to rad_thin_flds_setup by ',myname_sub
     return
   endif
!
! Determine if thinning is to be performed 
   if (fld_thin_box_size >= 1.) then
     l_thin_data=.true.
   else
     l_thin_data=.false.
   endif
!
! Determine dimensions of array required to hold header info for thinned obs
   thin_box_size=fld_thin_box_size 
   thin_info_dim=obs_info_num+fld_num_tcon+fld_num_tvar+4  ! see 4 added below
   thin_box_dim(1)=thin_info_dim
   thin_box_dim(3)=fld_thin_box_subtypes
   call rad_thin_set_boxes (lprint,ier2)
!
   if (ier2 == 0) then
     ier=0
   else
     ier=20
   endif
!
   allocate (thin_info(thin_info_dim))
   allocate (thin_info_names(thin_info_dim))
!
! add names of values to be added to obs_info
   m=obs_info_num
   thin_info_names(1:m)=obs_info_names(1:m) 
   do n=1,fld_num_tcon  ! add names of temporally constant fields
     m=m+1
     thin_info_names(m)=fld_names(1,n,1)
   enddo
   do n=1,fld_num_tvar  ! add names of temporally varying fields
     m=m+1
     thin_info_names(m)=fld_names(1,n,2)
   enddo
   thin_info_names(m+1)='NOBS'
   thin_info_names(m+2)='PENALTY'
   thin_info_names(m+3)='DHOURS'
   thin_info_names(m+4)='BOXCNT'
   thin_nobs_id=m+1
   thin_penalty_id=m+2
   thin_dhours_id=m+3
   thin_boxcnt_id=m+4
!
   call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLAT',thin_lat_id)
   if (thin_lat_id == 0) then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLATH',thin_lat_id)
   endif
   if (thin_lat_id == 0) then
     ier=ier+1
     print *,'Neither CLAT or CLATH found in obs_info_names' 
   endif  
!
   call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLON',thin_lon_id)
   if (thin_lon_id == 0) then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','CLONH',thin_lon_id)
   endif
   if (thin_lat_id == 0) then
     ier=ier+1
     print *,'Neither CLON or CLONH found in obs_info_names' 
   endif  
!
!
   call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','YEAR',thin_year_id)
   if (thin_year_id == 0) then
     ier=ier+1
     print *,'YEAR not found in obs_info_names' 
   endif  
!
   if (dtype == 'AIRS' .or. dtype == 'AMSUAAQUA') then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','SIID',thin_subset_id)
     if (thin_subset_id == 0) then
       ier=ier+1
       print *,'SIID not found in obs_info_names' 
     endif
   else
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','SAID',thin_subset_id)
     if (thin_subset_id == 0) then
       ier=ier+1
       print *,'SAID not found in obs_info_names' 
     endif
   endif
!
   if (dtype == 'CRIS' .or. dtype == 'CRISFSR' ) then
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','FOVN',thin_fovn_id)
     if (thin_fovn_id == 0) then
       ier=ier+1
       print *,'FOVN not found in obs_info_names' 
     endif
! 
     call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                   lstop,' ','QMRKH',thin_qmrkh_id)
     if (thin_qmrkh_id == 0) then
       ier=ier+1
       print *,'QMRKH not found in obs_info_names' 
     endif 
   endif
!
   if (ier /= 0) then 
     print ('(a,i5,2a)'),'Error: ier=',ier,' in ',myname_sub
     print ('(a)'),'List of available obs_info_names is:'
     print ('(8(1x,a12))'),obs_info_names(1:obs_info_num)
     if (lstop) then
       stop
     endif
   endif
!
   end subroutine rad_thin_setup 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_set_boxes (lprint,ier)
!
!  Define locations of thining boxes and allocate arrays to hold locations 
!  and obs info 
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
!
   integer(8) :: n1,n2
   integer :: n
   integer :: nboxes
   real(rkind1), parameter :: earth_r=6370.  ! radius of earth in km
   real(rkind1) :: pi, lat, circum
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_set_boxes'
!
   ier=0
!
   if (l_thin_data) then  ! determine arrays based on thinning specifications
     pi=4.*atan(1.)
     thin_nlats=int(pi*earth_r/thin_box_size)
     thin_dlat=r180/thin_nlats
     allocate (thin_nlonP(thin_nlats))
     allocate (thin_nlons(thin_nlats))
     allocate (thin_dlons(thin_nlats))
!      
     nboxes=0
     do n=1,thin_nlats
       lat=-r90+(n-.5)*thin_dlat
       thin_nlonP(n)=nboxes
       circum=2.*pi*earth_r*cos(pi*lat/r180)
       thin_nlons(n)=max(1,int(circum/thin_box_size))
       thin_dlons(n)=r360/thin_nlons(n)
       nboxes=nboxes+thin_nlons(n)
     enddo
!
   else  ! no thinning to be performed
     no_thin_id=0  ! initialize counter to 0
     n1=rkind1*thin_box_dim(1)*thin_box_dim(3)
     n2=no_thin_max_bytes/n1
     nboxes=n2
     allocate (thin_nlonP(1))  ! not used
     allocate (thin_nlons(1))  ! not used
     allocate (thin_dlons(1))  ! not used
   endif  
!
   thin_box_dim(2)=nboxes
!
   allocate (thin_box_info(thin_box_dim(1),thin_box_dim(2),thin_box_dim(3)))
   thin_box_info(:,:,:)=0.
!
   allocate (thin_box_counter(obs_time_slots,thin_box_dim(3)))
!
   if (lprint) then
     print *,' '
     if (l_thin_data) then 
       print ('(a,i10,a)'),' Thinning boxes defined for ', &
               thin_box_dim(2),' boxes'
       print ('(a,f8.1,a,i4,f6.2,a,i3,a,i3)'),' box_size=',thin_box_size, &
               ', nlats,dlat=',thin_nlats,thin_dlat,', ntypes=', &
               thin_box_dim(3),', info size=',thin_box_dim(1)
     else
       print ('(a)'),' No obs thinning will be performed by this program'
       print ('(2a,i7)'),' The max number of obs locations to be', &
                         ' considered is ',thin_box_dim(2)
     endif
   endif       
!
   end subroutine rad_thin_set_boxes
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_clean
!
!  Deallocate arrays defined in module rad_thin
!
   use m_rad_thin_flds, only : rad_thin_flds_clean
!
   deallocate (thin_dlons,thin_nlons,thin_nlonP)
   deallocate (thin_box_info,thin_info,thin_info_names)
   deallocate (thin_box_counter)
!
   call rad_thin_flds_clean
!
   end subroutine rad_thin_clean
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_put (cdtime_old,cdtime_new,n_errors,dtype, &
                            nobs,ierrors,lerror)
!              
! Copy required obs header into thinning box array based on computed penalty.
! Also replace obs time if requested. This change in time maintains the 
! difference between the time in the obs report header and idatetime_old.
!
   use m_time_compute, only  : time_compute_add, time_compute_unpack 
   use m_rad_thin_flds, only : rad_thin_flds_interp
   use m_rad_thin_flds, only : fld_num_all 
!
   implicit none
!
   integer,intent(in) :: n_errors
   integer,intent(in) :: nobs
   integer,intent(inout) :: ierrors(n_errors)
   logical,intent(out) :: lerror
   character(len=*),intent(in) :: cdtime_old
   character(len=*),intent(in) :: cdtime_new
   character(len=*),intent(in) :: dtype
!
! Local variables
   logical :: lerror1
   integer :: ier
   integer :: box_id
   integer :: type_id
   integer, parameter :: nfintp=100  ! should be bigger than 10 
   real(rkind1) :: fintp(nfintp)
   real(rkind1) :: dhours
   real(rkind1) :: penalty
   real(rkind1) :: obs_ymdhms(6)
   real(rkind1) :: obs_ymdhms_new(6)
   real(rkind1) :: obs_lat
   real(rkind1) :: obs_lon
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_put'
!
   thin_info(1:obs_info_num)=obs_info(1:obs_info_num)
   obs_lat=thin_info(thin_lat_id)
   obs_lon=thin_info(thin_lon_id)
   obs_ymdhms(1:6)=thin_info(thin_year_id:thin_year_id+5)
!
! Compute departure from snoptic time and test some obs header info 
   call rad_thin_check_hdr (n_errors,cdtime_old,obs_lat,obs_lon, &
                            obs_ymdhms,dtype,dhours,ierrors,lerror1)
   lerror=lerror1
!
! Determine satellite or instrument subset id if obs header is OK
   if (.not. lerror1) then 
     call rad_thin_find_subset (type_id)
     if (type_id == 0) then
       ierrors(7)=ierrors(7)+1
       lerror1=.true.
     endif
   endif
!
! Replace time information if obs header is OK
   if (.not. lerror1 .and. cdtime_old /= cdtime_new) then
     call time_compute_unpack (cdtime_new,obs_ymdhms)
     call time_compute_add (dhours,obs_ymdhms,obs_ymdhms_new,ier) 
     thin_info(thin_year_id:thin_year_id+5)=obs_ymdhms_new(1:6)
   endif
!
! If obs header Ok, then interpolate required fields
   if (.not.lerror1) then  
     call rad_thin_flds_interp (obs_lat,obs_lon,dhours,nfintp,fintp)
     call rad_thin_penalty (nfintp,fintp,penalty)
!
     thin_info(obs_info_num+1:obs_info_num+fld_num_all)= &
           fintp(1:fld_num_all)   
     thin_info(thin_nobs_id)=nobs
     thin_info(thin_penalty_id)=penalty
     thin_info(thin_dhours_id)=dhours
!
! Determin box index based either on location (if data is to be thinned)
! or simply on obs count  
     if (l_thin_data) then  ! determine box_id based on geographical location 
       call rad_thin_find_box_id (obs_lat,obs_lon,box_id)
     else                   ! determine box_id based nly on OK obs count
       no_thin_id=min(no_thin_id+1,thin_box_dim(2))
       box_id=no_thin_id
     endif
!
! Compare penalties to update saved thinning-box values with best choice
     call rad_thin_box_update (box_id,type_id)
     ierrors(n_errors)=ierrors(n_errors)+1  ! number of OK obs
!
   endif ! test on lerror1
!
   end subroutine rad_thin_put 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_check_hdr (n_errors,cdtime,obs_lat,obs_lon, &
                                  obs_ymdhms,dtype,dhours,ierrors, &
                                  lerror)
!              
! Check that some obs header info is within proper limits. 
!
   use m_time_compute, only : time_compute_unpack 
   use m_time_compute, only : time_compute_dhours 
!
   implicit none
!
   integer,intent(in) :: n_errors
   integer,intent(inout) :: ierrors(n_errors)
   real(rkind1),intent(in) :: obs_ymdhms(6)
   real(rkind1),intent(in) :: obs_lat
   real(rkind1),intent(in) :: obs_lon
   real(rkind1), intent(out) :: dhours
   logical,intent(out) :: lerror
   character(len=*),intent(in) :: cdtime
   character(len=*),intent(in) :: dtype
!
   logical :: lerrors1, lerrors2
   integer :: id
   integer :: ier
   real(rkind1) :: obs_tref(6)
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_check_hdr'
!
   lerror=.false.
!
   if (obs_lat > 1.e9) then  ! missing data
     thin_info(thin_info_dim)=1.e10
     ierrors(1)=ierrors(1)+1
     lerror=.true.
!
   else                          ! treatment of non-missing data
!
! Compute time difference between observation time in report header 
! record and time in message headers in bufr file
     call time_compute_unpack (cdtime,obs_tref)
     call time_compute_dhours (obs_tref,obs_ymdhms,dhours,ier)
     if (ier /= 0) then
       lerror=.true.
     else
       lerror=.false.  
     endif
!
! Check if time range of observation is acceptable
     if (dhours > 3.) then
       ierrors(2)=ierrors(2)+1
       lerror=.true.
     endif
     if (dhours < -3.) then
       ierrors(3)=ierrors(3)+1
       lerror=.true.
     endif
!
! Check lat and lon 

     if (obs_lon < -r180 .or. obs_lon > r360) then
       ierrors(4)=ierrors(4)+1
       lerror=.true.
     endif
!
     if (obs_lat < -r90 .or. obs_lat > r90) then
       ierrors(5)=ierrors(5)+1
       lerror=.true.
     endif
!
! If AIRS, check that sensor zenith angle <=80 (otherwise crtm crashes)
     if (trim(dtype) == 'AIRS') then 
       call find_name (obs_info_num,obs_info_names(1:obs_info_num), &
                       .false.,myname_sub,'SAZA',id)
       if (thin_info(id) > 80. ) then
         ierrors(6)=ierrors(6)+1 
         lerror=.true.
       endif
     endif        
!
! If CRIS, only include obs with FOVN = 5  (center of 9-spot view) 
! and obs with quality mark =0 in report header
     if (trim(dtype) == 'CRIS' .or. trim(dtype) == 'CRISFSR') then 
       if (nint(obs_info(thin_fovn_id)) /= 5 .or. &
           abs(obs_info(thin_qmrkh_id)) > .001 ) then 
         ierrors(6)=ierrors(6)+1
         lerror=.true.
       endif
     endif
! 
   endif ! test on whether header is readable 
!
   end subroutine rad_thin_check_hdr
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_find_subset (id) 
!
   use m_rad_thin_flds, only : fld_thin_box_subtypes, fld_thin_box_ids 
!
   implicit none
   integer, intent(out) :: id
!
   integer :: n
!
   id=0
   do n=1,fld_thin_box_subtypes
     if (nint(thin_info(thin_subset_id)) == fld_thin_box_ids(n)) then
       id=n
       exit
     endif
   enddo  
!
   end subroutine rad_thin_find_subset  
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_penalty (nfintp,fintp,penalty)
!
! Compute obs thining penalty as linear function of extracted NR values
!
   use m_rad_thin_flds, only : fld_num_tvar,fld_num_tcon,fld_weights
!
   implicit none
   integer, intent(in) :: nfintp
   real(rkind1), intent(in) :: fintp(nfintp)
   real(rkind1), intent(out) :: penalty
!
   integer :: k,n
!   
   penalty=0._rkind1
   n=0
   do k=1,fld_num_tcon
     n=n+1
     penalty=penalty+fld_weights(k,1)*fintp(n)
   enddo
!
   do k=1,fld_num_tvar
     n=n+1
     penalty=penalty+fld_weights(k,2)*fintp(n)
   enddo
!
   end subroutine rad_thin_penalty      
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_find_box_id (lat,lon,id)
!
! Find thining-box index based on obs location
! 
   implicit none
   real(rkind1), intent(in) :: lat
   real(rkind1), intent(in) :: lon
   integer, intent(out)     :: id
!
   integer :: nlat, nlon
   real(rkind1) :: lonx
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_find_box_id'
!
   nlat=1+int((lat+r90)/thin_dlat)
   nlat=min(nlat,thin_nlats)
!
   if (lon < 0._rkind1 ) then
     lonx=r360+lon
   else
     lonx=lon
   endif
   nlon=1+int(lonx/thin_dlons(nlat))
   nlon=min(nlon,thin_nlons(nlat))
!
   id=nlon+thin_nlonP(nlat)
!
   end subroutine rad_thin_find_box_id
! 
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_box_update (box_id,type_id)
!
! Place (or over-write) obs information saved in a thinning box based 
! on assigned penalty and other conditions
!
   implicit none
   integer, intent(in)  :: box_id   ! index of box number (location)
   integer, intent(in)  :: type_id  ! index of observation sub-type
!
   logical :: cond1, cond2, cond3
   real(rkind1) :: cnt
   real(rkind1) :: dp
   real(rkind1) :: dh
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_box_update'
!
! the following conditional tests whether
!  (1) this is the first obs falling in the thinning box, and therefore
!      this obs's scalar info should be saved, or
!  (2) this obs is in a clearer profile than any previous obs considered
!      for this thinning box (i.e., penalty smaller), and therefore this 
!      new obs should be saved in place of any earlier one.
!  (3) if the latest best obs and the current one considered have the same 
!      penalty value, then retain the obs closer to the central synoptic time.
!
! First increment counter of obs falling within this box
   thin_box_info(thin_boxcnt_id,box_id,type_id)= &
           thin_box_info(thin_boxcnt_id,box_id,type_id)+1.0_rkind1
   thin_info(thin_boxcnt_id)= &
           thin_box_info(thin_boxcnt_id,box_id,type_id)
!
   dp=thin_box_info(thin_penalty_id,box_id,type_id)-thin_info(thin_penalty_id)
   dh=abs(thin_box_info(thin_dhours_id,box_id,type_id))- &
      abs(thin_info(thin_dhours_id))
   cond1=nint(thin_box_info(thin_boxcnt_id,box_id,type_id)) == 1
   cond2=dp > 0._rkind1
   cond3=dp == 0._rkind1 .and. dh > 0._rkind1 
   if (cond1 .or. cond2 .or. cond3) then 
     thin_box_info(:,box_id,type_id)=thin_info(:)
   endif
!
   end subroutine rad_thin_box_update 
! 
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_box_count (lprint)
!
! Count numbers of thinning boxes that actually have observation 
! information as functions of obs subtype and time slot
!
   use m_rad_thin_flds, only : fld_thin_box_subtypes,fld_thin_box_ids
   implicit none
   logical, intent(in) :: lprint
!
   logical :: lcheck
   integer :: id,m,n
   integer :: nonempty_box
   real(rkind1) :: fracs(fld_thin_box_subtypes)  
   real(rkind1) :: fracsall
   real(rkind1) :: t0
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_box_count'
!
   thin_box_counter(:,:)=0.
!
! Separate counts for each intended obs time slot and each obs subtype
   do n=1,thin_box_dim(3)
     do m=1,thin_box_dim(2)
       if (thin_box_info(thin_boxcnt_id,m,n) > 0) then
         id=1+int((thin_box_info(thin_dhours_id,m,n)-obs_time_first)/  &
                  obs_time_delta)
         id=min(id,obs_time_slots-1)
         thin_box_counter(id,n)=thin_box_counter(id,n)+1
       endif
     enddo
   enddo 
!
! Sum counts over all obs time slots for each obs subtype
   do n=1,thin_box_dim(3)
     do id=1,obs_time_slots-1
       thin_box_counter(obs_time_slots,n)=thin_box_counter(id,n)+ &
           thin_box_counter(obs_time_slots,n)
     enddo
   enddo        
!
! Count nonempty_boxes (no obs of any subtype)
   nonempty_box=0
   do m=1,thin_box_dim(2)
     lcheck=.false.
     do n=1,thin_box_dim(3)
       if (thin_box_info(thin_boxcnt_id,m,n) > 0) then
         lcheck=.true.
       endif
     enddo
     if (lcheck) then
       nonempty_box=nonempty_box+1
     endif
   enddo
!
! Compute fractions of non-empty boxes for each subtype 
   do n=1,thin_box_dim(3)
     fracs(n)=real(thin_box_counter(obs_time_slots,n)) &
              /real(thin_box_dim(2))    
   enddo
   fracsall=real(nonempty_box)/real(thin_box_dim(2))    
!
   if (lprint) then
     print *,' '
     print *,'Table of counts of retained observations:'
     print ('(a)'),'slot_id,   t_end, subtype_id= '
     print ('(16x,10i8)'),fld_thin_box_ids(1:thin_box_dim(3))
     do id=1,obs_time_slots-1
       t0=obs_time_first+id*obs_time_delta
       print ('(i7,f9.4,10i9)'), id,t0, &
           thin_box_counter(id,1:thin_box_dim(3)) 
     enddo
     print ('(a,12x,10i9)'),'sum=',     &
           thin_box_counter(obs_time_slots,1:thin_box_dim(3)) 
!
     if (l_thin_data) then  ! print fractions of thinning boxes with obs
       print *,'fractions of non-empty boxes for each separate sub type:'
       print ('(16x,10f10.5)'),fracs(:)
       print *,'fraction of boxes with some kind of obs =',fracsall
     endif
   endif
!
   end subroutine rad_thin_box_count
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine rad_thin_write (file_name,obs_list_header_nrecs,lprint, &
                              obs_list_header_crecs,modtype,ier)
!
! Write thinned obs list to file ordered by time slot by pulling observation
! information saved in thinning box arrays.
!
   use m_rad_obs_arrays, only : obs_n_channels
   use m_rad_obs_arrays, only : obs_ndim1,obs_info_extra_recs
   use m_rad_obs_arrays, only : obs_info_extra_names,obs_info_extra
!
   use m_set_unit_nums, only : un_list_write
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(in)  :: obs_list_header_nrecs
   integer, intent(out) :: ier
   character(len=*), intent(in)   :: file_name
   character(len=*), intent(in)   :: modtype
   character(len=240), intent(in) :: &
                                 obs_list_header_crecs(obs_list_header_nrecs)
!
   integer, parameter :: obs_list_format_header=1
   integer, parameter :: obs_list_format_recs=1
   integer, parameter :: obs_list_len_i=1
   integer, parameter :: obs_list_len_c=1
   integer, parameter :: iunit=un_list_write
   integer, allocatable :: obs_list_tslot_recs(:)
   integer :: id,m,n,nc,nt
   integer :: obs_list_tslots
   integer :: obs_list_len_r
   integer :: obs_list_i(obs_list_len_i)
   character(len=16) :: obs_list_names_i(obs_list_len_i)
   character(len=16) :: obs_list_names_c(obs_list_len_c)
   character(len=16) :: obs_list_c(obs_list_len_c)
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_write'
!
   ier=0 
   open (iunit,file=trim(file_name), form='unformatted')
!
! Write file format indicators 
   write (iunit) obs_list_format_header, obs_list_format_recs
   write (iunit) obs_n_channels,obs_list_header_nrecs
   write (iunit) obs_ndim1,obs_info_extra_recs
   obs_list_len_r=thin_box_dim(1)
   write (iunit) obs_list_len_i, obs_list_len_r, obs_list_len_c, &
                 obs_time_slots, obs_time_delta, obs_time_first, &
                 thin_box_dim(3)
!
!  Write obs_list_arrays
   do n=1,obs_list_header_nrecs
     write (iunit) obs_list_header_crecs(n)
   enddo
!
   write (iunit) obs_info_extra_names(1:obs_info_extra_recs)
   do n=1,obs_info_extra_recs
     write (iunit) obs_info_extra(:,n)
   enddo
!
   obs_list_tslots=obs_time_slots-1
   allocate (obs_list_tslot_recs(obs_list_tslots))
   do id=1,obs_list_tslots
     obs_list_tslot_recs(id)=0
     do n=1,thin_box_dim(3)
       obs_list_tslot_recs(id)=obs_list_tslot_recs(id)+ &
              thin_box_counter(id,n)
     enddo
   enddo
   write (iunit) obs_list_tslot_recs,thin_box_counter
!
! Write names of time-varying information
   obs_list_names_i(1)='num_in_tslot'
   obs_list_names_c(1)='bufrmsglab'
   write (iunit) obs_list_names_i(:), &
                 thin_info_names(1:obs_list_len_r),obs_list_names_c(:)
!   
! Write time-varing information, sorted into separate groups for
! each time slot and obs subtype
   do nt=1,obs_time_slots-1 
     nc=0
     do n=1,thin_box_dim(3)
       do m=1,thin_box_dim(2)
         if (thin_box_info(thin_boxcnt_id,m,n) > 0) then
           id=1+int((thin_box_info(thin_dhours_id,m,n)-obs_time_first)/  &
                     obs_time_delta)
           id=min(id,obs_time_slots-1)
           if (id == nt) then
             nc=nc+1
             obs_list_i(1)=nc
             obs_list_c(1)=trim(modtype) 
             write (iunit) obs_list_i,thin_box_info(:,m,n),obs_list_c
           endif
         endif
       enddo
     enddo
   enddo 
!
   close (iunit)
!
   if (lprint) then
     print *,' '
     print *,'Obs file written in ',myname_sub
     print ('(a,2i5)'),'file format =',obs_list_format_header, &
                  obs_list_format_recs 
     print ('(a,i5)'),'number of obs channels = ',obs_n_channels
     print ('(a,i5)'),'number of descriptor records = ',obs_list_header_nrecs
     do n=1,obs_list_header_nrecs
       print ('(a,i1,2a)'),'descriptor ',n,': ', &
                   trim(obs_list_header_crecs(n))
     enddo  
     if (obs_info_extra_recs > 0) then
       print ('(11a)'),'t-constant info list: ', &
                       obs_info_extra_names(1:obs_info_extra_recs)
     else
       print ('(a)'),'t-constant info list: none'
     endif
     print ('(2a)'),'obs_list_len_i, obs_list_len_r, obs_list_len_c, ', &
                    'obs_list_tslots, obs_subtypes : '               
     print ('(5i5)'),obs_list_len_i,obs_list_len_r,obs_list_len_c,      &   
                     obs_list_tslots,thin_box_dim(3) 
     print ('(5a20)'),'integer variables: ',obs_list_names_i(:)
     print ('(5a20)'),'real variables: ', &
                   thin_info_names(1:obs_list_len_r)
     print ('(5a20)'),'character variables: ',obs_list_names_c(:)

     print ('(a,12i8)'),'num obs in tslot: ',obs_list_tslot_recs(:)
   endif
!
   deallocate (obs_list_tslot_recs)
!
   end subroutine rad_thin_write 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   end module m_rad_thin

