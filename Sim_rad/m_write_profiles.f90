   module m_write_profiles 
!
! Module to write a file of profiles created by program create_rad_profs
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
   use m_set_unit_nums, only : un_prof_out
!
   implicit none
!
   private 
   public :: write_profiles_header
   public :: write_profiles_recs
   public :: write_profiles_sample
   public :: write_profiles_close
   public :: write_profiles_tv2t
   public :: write_profiles_o3units
!
   integer :: profile_unit
   integer :: profile_format_header
   integer :: profile_format_recs
!
   contains
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_header (profile_file,format_hdr,format_rec)
!
!  Write file header to file of profiles for one observation type.
!  This contains index values, counts, and other info common to all obs of
!  this type, including names of the fields included.
!
   use m_obs_list, only : obs_list_tslots, obs_list_tslots_recs
   use m_obs_list, only : obs_list_time_delta, obs_list_time_first
   use m_obs_list, only : obs_list_header_nrecs, obs_list_header_crecs
   use m_obs_list, only : obs_list_n_channels, obs_list_subtypes
   use m_obs_list, only : obs_list_extra_dim1, obs_list_extra_recs
   use m_obs_list, only : obs_list_extra_names, obs_list_extra_info
   use m_obs_list, only : obs_list_counter
!
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_kmax, field_akbk, field_names
   use m_nr_fields_info, only : field_common_path
!
   implicit none
! 
   integer, intent(in) :: format_hdr
   integer, intent(in) :: format_rec
   character(len=*), intent(in) :: profile_file
!
   integer, parameter :: default_format_header=1
   integer, parameter :: default_format_recs=1
   integer :: n
   integer :: iunit
   character(len=*), parameter :: my_name='write_profiles_header'
!
   profile_unit=un_prof_out
! 
! Indicate either specified format (e.g., necessary when copying header) 
! Or set to dafault value
   if (format_hdr > 0) then
     profile_format_header=format_hdr
   else
     profile_format_header=default_format_header
   endif
!
   if (format_rec > 0) then
     profile_format_recs=format_rec
   else
     profile_format_recs=default_format_recs
   endif
!
!  Open file 
   iunit=profile_unit
   open (iunit,file=trim(profile_file), form='unformatted')
!
!  Write obs_list_ scalars
   write (iunit) profile_format_header, profile_format_recs
   write (iunit) obs_list_n_channels, obs_list_header_nrecs
   write (iunit) obs_list_extra_dim1,obs_list_extra_recs
   write (iunit) obs_list_len_i,obs_list_len_r,obs_list_len_c, &
                 obs_list_tslots,obs_list_time_delta,          & 
                 obs_list_time_first,obs_list_subtypes
!
!  Write obs_list_ arrays
   do n=1,obs_list_header_nrecs
     write (iunit) obs_list_header_crecs(n)
   enddo
!
   write (iunit) obs_list_extra_names(1:obs_list_extra_recs)
   do n=1,obs_list_extra_recs
     write (iunit) obs_list_extra_info(:,n)
   enddo
! 
   write (iunit) obs_list_tslots_recs,obs_list_counter
   write (iunit) obs_list_names
!
!  Write field info
   write (iunit) field_num_2d,field_num_3d,field_kmax 
   write (iunit) field_akbk(1:field_kmax+1,1:2)
   write (iunit) field_names(1:field_num_2d,1,2),field_names(1:field_num_3d,1,3)
   write (iunit) field_common_path
!
   end subroutine write_profiles_header
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_recs (ndim,nob,profile)
!
! Write observation report header and profile fields for one geographic 
! location.  Note that this header generally contains both obs meta-data
! and values for some 2-d fields at the obs location. 
!
   implicit none
! 
   integer, intent(in) :: ndim 
   integer, intent(in) :: nob
   real(rkind1), intent(in) :: profile(ndim)
!
   integer :: iunit
   character(len=*), parameter :: my_name='write_profiles_recs'
!  
   profile_unit=un_prof_out
   iunit=profile_unit
   write (iunit) obs_list_i(:,nob),obs_list_r(:,nob),obs_list_c(:,nob) 
   write (iunit) profile
!
   end subroutine write_profiles_recs
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_sample (ndim,nob,dtype,profile)
!
!  Print profile field values for one observation location
!  (Used for printing a sample of results when requested)
!
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_kmax, field_names
!
   implicit none
! 
   integer, intent(in) :: ndim 
   integer, intent(in) :: nob
   real(rkind1), intent(in) :: profile(ndim)
   character(len=*), intent(in) :: dtype
!
   integer :: k,n,nlast
   real(rkind1), allocatable :: prof_2d(:)
   real(rkind1), allocatable :: prof_3d(:,:)
   character(len=*), parameter :: my_name='write_profiles_recs'
!  
   allocate (prof_2d(field_num_2d))
   allocate (prof_3d(field_kmax,field_num_3d))
!
   call prof_all_2d3d (field_num_2d,field_num_3d,field_kmax,ndim,profile, &
                       prof_2d,prof_3d)
!
   print *,' '
   print ('(3a,i7)'),'Sample profile for ',trim(dtype),' obs number =',nob
   print *,'values for 2-d field follow:'
   do n=1,field_num_2d
     print ('(a12,a,1p1e15.5)'),field_names(n,1,2),'=',prof_2d(n)
   enddo
   print *,'values for 3-d field follow:'
   do n=1,field_num_3d,6
     nlast=min(n+5,field_num_3d)
     print ('(a,8(1x,a12))'),'   k',field_names(n:nlast,1,3)
     do k=1,field_kmax
       print ('(i3,1p8e13.5)'),k,prof_3d(k,n:nlast)
     enddo
   enddo
!   
   deallocate (prof_2d,prof_3d)
!
   end subroutine write_profiles_sample
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_tv2t (ndim,nob,profile)
!
!  Change virtual temperature to temperature if the field values read and
!  stored as the named "T" field were actually virtual T. This change is 
!  only made if the corresponding QV field is also present.
!
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_kmax, field_names
!
   implicit none
! 
   integer, intent(in) :: ndim
   integer, intent(in) :: nob 
   real(rkind1), intent(inout) :: profile(ndim)
!
   integer :: ntpoint,nqpoint,nt,nq,n 
   real(rkind1), parameter :: eps1=0.622  ! ratio of weight of dry air to wvap
   character(len=*), parameter :: my_name='write_profiles_tv2t'
!
   ntpoint=0  
   nqpoint=0  
   do n=1,field_num_3d
     if (field_names(n,1,3) == 'T' .and. (field_names(n,2,3) == 'TV' .or. &
         field_names(n,2,3) == 'tv' .or.  field_names(n,2,3) == 'Tv')) then 
       ntpoint=n
     endif  
     if (field_names(n,1,3) == 'QV') then
       nqpoint=n
     endif  
   enddo
!
   if (ntpoint /= 0 .and. nqpoint /= 0) then
     ntpoint=field_num_2d+field_kmax*(ntpoint-1)      
     nqpoint=field_num_2d+field_kmax*(nqpoint-1)
     do n=1,field_kmax
       nt=n+ntpoint     
       nq=n+nqpoint     
       profile(nt)=profile(nt)/(1.+eps1*profile(nq))
     enddo
     if (nob == 1) then
       print ('(a,2i6)'),'TV converted to T: ntpoint,nqpoint=', ntpoint,nqpoint
     endif
   endif 
!
   end subroutine write_profiles_tv2t
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_o3units (ndim,nob,profile)
!
!  Change units of ozone from ozone_volume_mixing_ratio to 
!  ozone_mass_mixing_ratio if indicated by ozone variable name
!  The OSSE software m_crtm_interface assumes ozone units are 
!  input as kg/kg and then changed to ppmv before passing to the CRTM.
!  Here, if the variable name on file is ozone, as it is in the GMAO GSI 
!  background files, the read values are assumed to be ppmv as in those files,
!  whereas if it is O3, then the units are assumed to be kg/kg, 
!  as required here.
!
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_kmax, field_names
!
   implicit none
! 
   integer, intent(in) :: ndim
   integer, intent(in) :: nob 
   real(rkind1), intent(inout) :: profile(ndim)
!
   integer :: nopoint,no,n 
   real(rkind1), parameter :: fact=1./604229.  ! factor to change units
   character(len=*), parameter :: my_name='write_profiles_o3units'
!
   nopoint=0  
   do n=1,field_num_3d
     if (field_names(n,1,3) == 'O3' .and. field_names(n,2,3) == 'ozone') then
       nopoint=n
     endif  
   enddo
!
   if (nopoint /= 0) then
     nopoint=field_num_2d+field_kmax*(nopoint-1)      
     do n=1,field_kmax
       no=n+nopoint     
       profile(no)=fact*profile(no)
     enddo
     if (nob == 1) then
       print ('(a,i6)'),'O3 converted from ppmv to kg/kg: nopoint=', nopoint
     endif
   endif 
!
   end subroutine write_profiles_o3units
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_profiles_close
!
!  Close FORTRAN unit to which profiles have been written 
!
   profile_unit=un_prof_out
   close (profile_unit)
!
   end subroutine write_profiles_close
!
! 
   end module m_write_profiles 
