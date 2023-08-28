   module m_iasi_radiance
!
! IASI BUFR data is in terms of radiances, but pert code and CRTM requires
! brightness T.  This module uses CRTM routines to perform the required
! transforms.
!
   use crtm_module, only: crtm_destroy, crtm_init
   use crtm_planck_functions, only: crtm_planck_temperature
   use crtm_planck_functions, only: crtm_planck_radiance
   use crtm_module, only: crtm_channelinfo_type
!
   implicit none
   private
!
   public :: iasi_radiance_setup
   public :: iasi_radiance_cleanup
   public :: iasi_radiance_transforms
!
! For CRTM initialization to use planck function calls
   type(crtm_channelinfo_type),dimension(1) :: channelinfo
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine iasi_radiance_setup (crtm_coef_dir,error_status)
!
   implicit none
!
   integer, intent(out) :: error_status 
   character(len=*), intent(in) :: crtm_coef_dir
!
   character(len=256) :: sensor_id(1)
   character(len=256) :: file_prefix='iasi616_metop-a'
!                                                 !
   sensor_id(1)=ADJUSTL(file_prefix)
   error_status=crtm_init(ChannelInfo=channelinfo,            &
           Sensor_ID=sensor_id,File_Path=trim(crtm_coef_dir), &
           Load_CloudCoeff=.false.,Load_AerosolCoeff=.false.)
!   
   end subroutine iasi_radiance_setup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine iasi_radiance_cleanup (error_status)
!
   implicit none
!
   integer :: error_status 
!
   error_status=crtm_destroy(channelinfo)
!
   end subroutine iasi_radiance_cleanup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine iasi_radiance_transforms (nchans,obs_values,func)
!
   use m_kinds, only : rkind2
   implicit none
!
   integer :: nchans
   real(rkind2) :: obs_values(nchans)
   character(len=*) :: func
!
   integer :: i
   real(rkind2) :: x

!
   if (func == 'R2TB') then   ! transform radiance to BT 
     do i=1,nchans
       call crtm_planck_temperature (1,i,obs_values(i),x)
       obs_values(i)=x
     enddo
   else
     do i=1,nchans            ! transform BT to radiance 
       call crtm_planck_radiance (1,i,obs_values(i),x)
       obs_values(i)=x
     enddo
   endif
!
   end subroutine iasi_radiance_transforms
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_iasi_radiance
