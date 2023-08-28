   module m_rad2bt
!
! Module to change between radiances and brightness temperatures using 
! the CRTM.
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use crtm_module, only: crtm_destroy, crtm_init
   use crtm_planck_functions, only: crtm_planck_temperature
   use crtm_planck_functions, only: crtm_planck_radiance
   use crtm_module, only: crtm_channelinfo_type
!
   implicit none
   private
!
   public :: rad2bt_setup
   public :: rad2bt_cleanup
   public :: rad2bt_transforms
!
! For CRTM initialization to use planck function calls
   type(crtm_channelinfo_type),dimension(1) :: channelinfo
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad2bt_setup (crtm_coef_dir,dtype,said,error_status)
!
!  Set up CRTM to use proper spec and tau coefficients for transforms 
!  between radiances and TB.
!
   implicit none
!
   integer, intent(out) :: error_status 
   integer, intent(in) :: said                ! satellite id number
   character(len=*), intent(in) :: crtm_coef_dir
   character(len=*), intent(in) :: dtype
!
   character(len=256) :: sensor_id(1)
   character(len=256) :: file_prefix
!                    
   if (trim(dtype) == 'IASI' .and. said == 4) then 
       file_prefix='iasi616_metop-a'
   elseif (trim(dtype) == 'IASI' .and. said == 3) then 
       file_prefix='iasi616_metop-b'
   elseif (trim(dtype) == 'CRIS') then
     file_prefix='cris399_npp'
   elseif (trim(dtype) == 'CRISFSR') then
     file_prefix='cris-fsr_n20'
   else
     print *,' '
     print *,' ERROR IN rad2bt_setup: unaccepted request for dtype =', &
             trim(dtype)
     file_prefix='NONE'
   endif   
!                             !
   sensor_id(1)=ADJUSTL(file_prefix)
   error_status=crtm_init(ChannelInfo=channelinfo,            &
           Sensor_ID=sensor_id,File_Path=trim(crtm_coef_dir), &
           Load_CloudCoeff=.false.,Load_AerosolCoeff=.false.)
!   
   end subroutine rad2bt_setup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad2bt_cleanup (error_status)
!
!  The reverse of rad2bt_setup
!
   implicit none
!
   integer :: error_status 
!
   error_status=crtm_destroy(channelinfo)
!
   end subroutine rad2bt_cleanup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad2bt_transforms (nchans,obs_values,func)
!
! Transforms back and forth between radiances and brightness temperatures.
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
   end subroutine rad2bt_transforms
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_rad2bt
