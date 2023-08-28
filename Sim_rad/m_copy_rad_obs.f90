   module m_copy_rad_obs
!
! Module used to copy arrays defined within the m_rad_obs_arrays module 
! to corresponding arrays defined within the m_bufr_rad module
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!   
   use m_kinds, only : rkind1, rkind2
   use m_rad_obs_arrays, only : obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4
   use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
   use m_rad_obs_arrays, only : obs_info_num,obs_info_names,obs_info_hdr
   use m_rad_obs_arrays, only : obs_n_channels, obs_n_data
   use m_rad_obs_arrays, only : obs_ndim5,obs_info_extra_recs
   use m_rad_obs_arrays, only : obs_info_extra_names,obs_info_extra
!
   use m_obs_list, only : obs_list_extra_recs,obs_list_extra_info
!
   implicit none
   private
! 
   public :: copy_rad_obs
!
   contains 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine copy_rad_obs (n_chan,b_chan,n_info,dtype,quality, &
                            crtm_obs,prof_info,c_realloc)

   use m_rad_index, only : nf_chnm, nf_spec, nf_fcph
   use m_bufr_rad, only : bmiss
!
   implicit none
!
   integer :: n_chan     ! number of channels considered by CRTM
   integer :: n_info     ! size of array prof_info (aka obs_list_r)
   integer :: b_chan     ! number of channels in BUFR file
   real(rkind1), intent(in) :: quality
   real(rkind1), intent(in) :: prof_info(n_info)
   real(rkind1), intent(in) :: crtm_obs(n_chan)
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: c_realloc
!
   integer :: n
!
   obs_info(1:n_info)=prof_info(1:n_info)         ! copy header info
   obs_values(1:n_chan,1)=crtm_obs(1:n_chan)      ! copy rad or BT
   obs_n_data=1
   obs_n_channels=b_chan
!
! copy channel numbers (generally in this array)
   if (nf_chnm > 0) then                    
     obs_info_extra(1:n_chan,1)=obs_list_extra_info(1:n_chan,nf_chnm)
   else 
     do n=1,n_chan
       obs_info_extra(n,1)=real(n,rkind2)
     enddo 
   endif 
   obs_channels(1:n_chan,1)=obs_info_extra(1:n_chan,1)
!     
! copy other location-independent channel-dependent information 
  if (nf_spec > 0) then
     obs_info_extra(1:n_chan,2)=obs_list_extra_info(1:n_chan,nf_spec)
   else 
     obs_info_extra(1:n_chan,2)=0._rkind2
   endif 
   obs_channels(1:n_chan,2)=obs_info_extra(1:n_chan,2)
!
! If GMI, there are 4 channel-dependent variables stored in array 
! obs_list_extra_info.  
  if (trim(dtype) == 'GMI') then
    obs_info_extra_recs=obs_list_extra_recs
    obs_info_extra(1:n_chan,1:obs_info_extra_recs)= &
        obs_list_extra_info(1:n_chan,1:obs_info_extra_recs)
  endif 
!      
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then
     obs_n_data=2
     obs_values(1:n_chan,2)=0.  ! quality flag
   endif
!
! Fill channel-independent quality info if assigned
   if (nf_fcph > 0) then  
     if (trim(dtype) == 'IASI' .and. c_realloc == 'T') then ! real locations for IASI
       obs_info(nf_fcph)=99. ! set almost clear so this will not affect QC
     else
       obs_info(nf_fcph)=quality
     endif              
   endif
!
! Fill empty portions of channel-dependent array with missing values
   if (b_chan > n_chan) then
     obs_values(n_chan+1:b_chan,1:obs_n_data)=bmiss
     obs_channels(n_chan+1:b_chan,1:2)=bmiss
   endif
!
   end subroutine copy_rad_obs 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_copy_rad_obs 

