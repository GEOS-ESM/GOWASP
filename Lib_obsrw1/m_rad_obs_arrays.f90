   module m_rad_obs_arrays
!
!  Module containing arrays and variables used to hold observation 
!  information for radiance data types
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1, rkind2
!
   private
   public :: rad_obs_arrays_setup 
   public :: rad_obs_arrays_clean
!
! Specify maximum dimensions for obs data arrays
! obs_ndim1: maximum number channels for observations 
! obs_ndim2: maximum number of observation data info in buffer 
!            (2 includes both brightness temperature and quality mark)
! obs_ndim3: number of channel info
!            (2 includes both channel number and log reciprical wavelength)
! obs_ndim4: maximum number of single value report info (lat, lon etc.)
! obs_ndim5: maximum number of extra info arrays considered
!
   integer, public :: obs_ndim1
   integer, public :: obs_ndim2
   integer, public :: obs_ndim3
   integer, public :: obs_ndim4
   integer, public :: obs_ndim5
!
   integer, public :: obs_info_extra_recs
   integer, public :: obs_info_hdr
   integer, public :: obs_info_num
   integer, public :: obs_n_channels
   integer, public :: obs_n_data
   integer, public :: obs_time_slots
   integer, public :: obs_generic_int(3) 
!
   real(rkind1), public :: obs_time_first
   real(rkind1), public :: obs_time_delta
   real(rkind2), public, allocatable :: obs_info(:) ! id,lon,lat,time,type,etc.
   real(rkind2), public, allocatable :: obs_channels(:,:) ! channel id and freq
   real(rkind2), public, allocatable :: obs_values(:,:)   ! B-temp and qc flag
   real(rkind2), public, allocatable :: obs_info_extra(:,:) ! extra info
!
   character(len=16), public, allocatable :: obs_info_names(:) 
   character(len=16), public, allocatable :: obs_info_extra_names(:) 
   character(len=20), public :: obs_generic_char(5)  ! len to accommodate sis 
!
   contains 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad_obs_arrays_setup (ndim1,ndim2,ndim3,ndim4,ndim5)
!
!  Allocate some arrays to hold rad obs info
!
   implicit none
   integer, intent(in) :: ndim1,ndim2,ndim3,ndim4,ndim5
!
   obs_ndim1=ndim1
   obs_ndim2=ndim2   
   obs_ndim3=ndim3
   obs_ndim4=ndim4
   obs_ndim5=ndim5
!  
   allocate(obs_info(obs_ndim4))
   allocate(obs_channels(obs_ndim1,obs_ndim3))
   allocate(obs_values(obs_ndim1,obs_ndim2))
   allocate(obs_info_extra(obs_ndim1,obs_ndim5))
   allocate(obs_info_names(obs_ndim4))
   allocate(obs_info_extra_names(obs_ndim5))
!
   end subroutine rad_obs_arrays_setup 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad_obs_arrays_clean 
!
!  Deallocate arrays declared in rad_obs_arrays_setup
! 
   implicit none
!
   if (allocated(obs_info)) deallocate(obs_info)
   if (allocated(obs_values)) deallocate(obs_values)
   if (allocated(obs_channels)) deallocate(obs_channels)
   if (allocated(obs_info_extra)) deallocate(obs_info_extra)
   if (allocated(obs_info_names)) deallocate(obs_info_names)
   if (allocated(obs_info_extra_names)) deallocate(obs_info_extra_names)
!
   end subroutine rad_obs_arrays_clean
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_rad_obs_arrays  
