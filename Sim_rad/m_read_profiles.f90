   module m_read_profiles 
!
! Module to read file of profiles created by program create_rad_profs
! These pofiles are sets of 2d and 3d fields interpolated to observation
! locations for sets of radiance observations.  Header information about 
! each observation is included
!
! Initial Code: Ronald Errico  August 15 2014
!
   use m_kinds, only : rkind1
!   
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c  
   use m_obs_list, only : obs_list_i, obs_list_r, obs_list_c, obs_list_names
   use m_obs_list, only : obs_list_format_header, obs_list_format_recs
   use m_obs_list, only : obs_list_unit_in, obs_list_dim2
!
   implicit none
!
   private 
   public :: read_profiles_setup
   public :: read_profiles_recs
   public :: read_profiles_cleanup
!
! Public variables 
   integer, public :: prof_format_header
   integer, public :: prof_format_recs
   integer, public :: prof_num_2d, prof_num_3d, prof_kmax
   integer, public :: prof_dim1,prof_dim2
   real(rkind1), allocatable, public :: prof_akbk_int(:,:)
   real(rkind1), allocatable, public :: prof_all(:,:)
   character(len=240), public :: prof_common_path
   character(len=12), allocatable, public :: prof_names(:,:)
!
! Private variables 

   integer :: prof_unit
   character(len=*), parameter :: my_name='m_read_profiles'
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_profiles_setup (profile_file,ndim2,iers)
!
! Read header from file of profiles and allocate arrays to read profiles
!
   use m_obs_list, only: obs_list_read_header
   implicit none
! 
   integer, intent(in) :: ndim2  
   character(len=*), intent(in) :: profile_file
   integer, intent(out) :: iers
!
   integer :: n, n_max, i
   integer :: iunit
   character(len=*), parameter :: my_sub=my_name//'::write_profiles_header'
!
!  Open file of profiles and read header portion common to the obs_list_file  
   call obs_list_read_header (profile_file,iers)
   prof_unit=obs_list_unit_in
   prof_format_header=obs_list_format_header
   prof_format_recs=obs_list_format_recs
   obs_list_dim2=ndim2
   allocate (obs_list_i(obs_list_len_i,obs_list_dim2))   
   allocate (obs_list_r(obs_list_len_r,obs_list_dim2))   
   allocate (obs_list_c(obs_list_len_c,obs_list_dim2))   
!
!  Read header information concerning profile fields
   iunit=prof_unit
   read (iunit) prof_num_2d,prof_num_3d,prof_kmax
   allocate (prof_akbk_int(prof_kmax+1,2)) 
   read (iunit) prof_akbk_int   ! read ak, bk at data interfaces 
!
   n_max=max(prof_num_2d,prof_num_3d)
   allocate (prof_names(n_max,2:3))
   read (iunit) (prof_names(i,2),i=1,prof_num_2d), &
                (prof_names(i,3),i=1,prof_num_3d)
   read (iunit) prof_common_path
   prof_dim1=prof_num_2d+prof_num_3d*prof_kmax
   prof_dim2=ndim2
   allocate (prof_all(prof_dim1,prof_dim2))
!
   end subroutine read_profiles_setup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_profiles_recs (nob,lprint,ier)
!
! Read profiles for one geographic observing location
!
   implicit none
! 
   logical, intent(in)  :: lprint
   integer, intent(in)  :: nob
   integer, intent(out) :: ier
!
   integer :: iunit
   character(len=*), parameter :: my_sub=my_name//'::write_profiles_recs'
!  
   if (nob > prof_dim2) then
     print *,'Error in ',my_sub
     print *,'Attempt to read more profiles than allocated for array'
     print *,'nob,prof_dim2 = ',nob,prof_dim2
     ier=1
   else
     ier=0
     read (prof_unit) obs_list_i(:,nob),obs_list_r(:,nob),obs_list_c(:,nob) 
     read (prof_unit) prof_all(:,nob)
   endif
!

   end subroutine read_profiles_recs
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_profiles_cleanup (lclose)
!
!  Dealloacte arrays created by read_profiles_setup and obs_list_read_header
!
   use m_obs_list, only: obs_list_clean 
   implicit none
   logical, intent(in) :: lclose
   if (lclose) close (prof_unit)
   call obs_list_clean
   deallocate (prof_akbk_int,prof_names,prof_all)
!
   end subroutine read_profiles_cleanup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_read_profiles 
