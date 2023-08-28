   module m_rad_index
!
!  Module for setting indexes to denote locations of required fields 
!  in NR profile arrays input to the crtm interface (m_interface_crtm).
!
!  Initial Code: Ronald Errico August 15 2014
!
   use m_read_profiles, only : prof_num_2d, prof_num_3d, prof_names
   use m_read_profiles, only : prof_kmax
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c  
   use m_obs_list, only : obs_list_names
   use m_obs_list, only : obs_list_extra_names, obs_list_extra_recs
   use m_rad_prob, only : list_cloud_maxn, list_aerosol_maxn
!
   implicit none
!
   private
   public :: rad_index_setup
!
   integer, public :: nf_fovn, nf_saza, nf_soza, nf_hmsl
   integer, public :: nf_flnd, nf_flic, nf_fsic, nf_fsno
   integer, public :: nf_cldh, nf_cldm, nf_cldl
   integer, public :: nf_ps,   nf_ts,   nf_u10m, nf_v10m
   integer, public :: nf_sndp, nf_vegt, nf_vegf 
   integer, public :: nf_z0m,  nf_tsno, nf_swet
   integer, public :: nf_t,    nf_qv,   nf_o3,   nf_rh
   integer, public :: nf_siid, nf_said
   integer, public :: nf_lat, nf_lon, nf_time
   integer, public :: nf_satazm, nf_solazm
   integer, public :: nf_chnm, nf_spec, nf_fcph, nf_bmsg, nf_pang
   integer, public :: nf_aerosol_list(list_aerosol_maxn) ! indexes for a fields 
   integer, public :: nf_cloud_list(list_cloud_maxn)     ! indexes for c fields 
!
   character(len=*), parameter :: my_name='m_rad_index'
! 
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine rad_index_setup (dtype,lprint,ier) 
!
! Set indexes nf_.... used by m_crtm_interface to point to particular 
! values held in common arrays based on names of fields or values in prof_file
!
   use m_rad_prob, only : list_cloud_nums, list_aerosol_nums
   use m_rad_prob, only : list_cloud_names, list_aerosol_names 
!
   implicit none
!
   logical, intent(in) :: lprint
   character(len=*), intent(in) :: dtype
   integer, intent(out) :: ier 
!
   logical, parameter :: lstop=.false.
   logical :: lfound
   integer :: m1,m2,n1,n2,n3,n4,n5,n,i
   character(len=*), parameter :: mysub=my_name//'::setup_names'
!
   ier=0
!
!  Set required ids for required names in obs_list (profile header info)
!  In the calls to find_name, if the 4th argument is a character blank, 
!  then no error information will be printed if name not found in list. 
   m1=obs_list_len_i+1
   m2=obs_list_len_i+obs_list_len_r
   n1=obs_list_len_r
   call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SAZA',     nf_saza)
   call find_name (n1,obs_list_names(m1:m2),lstop,' ',  'FRLAND',   nf_flnd)
   call find_name (n1,obs_list_names(m1:m2),lstop,' ',  'FRLANDICE',nf_flic)
   call find_name (n1,obs_list_names(m1:m2),lstop,' ',  'FRSEAICE', nf_fsic)
   call find_name (n1,obs_list_names(m1:m2),lstop,' ',  'SIID',     nf_siid) 
   if (dtype /= 'AIRS' .and. dtype /= 'AMSUAAQUA') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SAID',nf_said)
   endif
   if (dtype == 'CRIS' .or. dtype == 'CRISFSR' ) then ! field of view index
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'FORN',nf_fovn) 
   else
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'FOVN',nf_fovn)
   endif
   if (dtype == 'GENRADTXT') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'PANGLR',nf_pang)
   endif
!
! Set indexes required when .bin file of clear-skly BT is written
! (In place of mysub a blank string is used so that no warning is printed
! until CLATH and CLONH are also searched)
     call find_name (n1,obs_list_names(m1:m2),.false.,' ','CLAT',nf_lat)
   if (nf_lat == 0) then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'CLATH',nf_lat)
   endif
!
   call find_name (n1,obs_list_names(m1:m2),.false.,' ','CLON',nf_lon)
   if (nf_lon == 0) then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'CLONH',nf_lon)
   endif
!
   call find_name (n1,obs_list_names(m1:m2),.false.,mysub,'DHOURS',nf_time)
   if (nf_lat == 0 .or. nf_lon == 0 .or. nf_time == 0 ) then 
     ier=3
     print *,'Error in ',mysub
     print ('(2a,7i4)'),'nf_lat,nf_lon,nf_time = ',nf_lat,nf_lon,nf_time 
   endif 
!
! Set QC flag index
   if (dtype == 'IASI') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'FCPH',nf_fcph)
   elseif (dtype == 'SSMIS') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'RFLAG',nf_fcph)
   elseif (dtype == 'CRIS' .or. dtype == 'CRISFSR') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'QMRKH',nf_fcph)
   elseif (dtype == 'AVCSAM' .or. dtype == 'AVCSPM') then
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'CLAVR',nf_fcph)
   else
     nf_fcph=0 ! no value used
   endif 
!
! Check solar zenith and azimuth angles 
   if (trim(dtype) == 'GMI') then 
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SZA',   nf_soza)   
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SMA',   nf_solazm)
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SAMA',  nf_satazm)
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'HMSL',  nf_hmsl)
   else
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SOZA',  nf_soza)
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'SOLAZI',nf_solazm) 
     call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'BEARAZ',nf_satazm)
     nf_hmsl=999  ! not used
   endif
!
   m1=obs_list_len_i+obs_list_len_r+1
   m2=obs_list_len_i+obs_list_len_r+obs_list_len_c
   n1=obs_list_len_c
   call find_name (n1,obs_list_names(m1:m2),lstop,mysub,'bufrmsglab',nf_bmsg) 
!  
   if (nf_fovn == 0 .or. nf_saza == 0 .or. nf_satazm == 0 .or. & 
       nf_bmsg == 0 .or. nf_soza == 0 .or. nf_solazm == 0 .or. &
       nf_hmsl == 0) then
     ier=ier+1
     if (lprint) then 
       print *,'Error in ',mysub
       print ('(2a,7i4)'),'nf_fovn,nf_saza,nf_soza,nf_bmsg,nf_satazm,', &
                          'nf_solazm,nf_hmsl =',                        &
                           nf_fovn,nf_saza,nf_soza,nf_bmsg,nf_satazm,   &
                           nf_solazm,nf_hmsl 
     endif 
   endif
!
!  Set required ids for required names in 2d field list (profile header info)
   n2=prof_num_2d 
   call find_name (n2,prof_names(1:n2,2),lstop,mysub,'PS',     nf_ps)
   call find_name (n2,prof_names(1:n2,2),lstop,mysub,'TS',     nf_ts)
   call find_name (n2,prof_names(1:n2,2),lstop,mysub,'U10M',   nf_u10m)
   call find_name (n2,prof_names(1:n2,2),lstop,mysub,'V10M',   nf_v10m)
   call find_name (n2,prof_names(1:n2,2),lstop,' ',   'vegtype',nf_vegt)
   call find_name (n2,prof_names(1:n2,2),lstop,' ',   'vegfrac',nf_vegf)
!
! If fields not among obs_list_names then look among prof_names
   if (nf_flnd == 0) then
     call find_name (n2,prof_names(1:n2,2),lstop,mysub,'FRLAND',nf_flnd)
     nf_flnd=nf_flnd+100
   endif 
   if (nf_flic == 0) then
     call find_name (n2,prof_names(1:n2,2),lstop,mysub,'FRLANDICE',nf_flic)
     nf_flic=nf_flic+100
   endif 
   if (nf_fsic == 0) then
     call find_name (n2,prof_names(1:n2,2),lstop,mysub,'FRSEAICE',nf_fsic)
     nf_fsic=nf_fsic+100
   endif 
!
   if (nf_flnd == 0 .or. nf_flic == 0 .or. nf_fsic == 0 .or. &
       nf_ps == 0 .or. nf_ts == 0) then
     ier=ier+10
     if (lprint) then 
       print *,'Error in ',mysub
       print ('(a,5i4)'),'nf_flnd,nf_flic,nf_fsic,nf_ps,nf_ts=', &
                          nf_flnd,nf_flic,nf_fsic,nf_ps,nf_ts
     endif
   endif
!
!  Set other ids for optional names in 2d field list 
   call find_name (n2,prof_names(1:n2,2),lstop,' ','SNODP',  nf_sndp)
   call find_name (n2,prof_names(1:n2,2),lstop,' ','TPSNOW', nf_tsno)
   call find_name (n2,prof_names(1:n2,2),lstop,' ','SFMC',   nf_swet)
   call find_name (n2,prof_names(1:n2,2),lstop,' ','FRSNO',  nf_fsno)
   call find_name (n2,prof_names(1:n2,2),lstop,' ','Z0M',    nf_z0m)
   if ((nf_sndp == 0 .or. nf_tsno == 0 .or. nf_swet == 0 .or. &
       nf_fsno == 0 .or. nf_z0m == 0) .and. lprint) then
     print *,'Fields with 0 values in the following list are not present ', &
        'among the profile fields; therefore default fields values will be used'
     print ('(a,5i4)'),'nf_sndp, nf_tsno, nf_swet, nf_fsno, nf_z0m =', &
        nf_sndp,nf_tsno,nf_swet,nf_fsno,nf_z0m
     print *,'Names of 2d fields on file (excluding in the obs headers ', &
        'written by create_rad_list):'
     print ('(10(1x,a10))'),(prof_names(i,2),i=1,n2)
   endif
!
!  Set required ids for required names in 3d field list 
   n3=prof_num_3d 
   call find_name (n3,prof_names(1:n3,3),lstop,mysub,'T',nf_t)
   call find_name (n3,prof_names(1:n3,3),lstop,mysub,'QV',nf_qv)
   call find_name (n3,prof_names(1:n3,3),lstop,mysub,'O3',nf_o3)
!
   if (nf_t == 0 .or. nf_qv == 0) then
     ier=ier+100
     if (lprint) then 
       print *,'Error in ',mysub
       print *,'nf_t,nf_qv=',nf_t,nf_qv
     endif
   endif
!
!  Set profile indexes for any cloud fields to be used
   n3=prof_num_3d
   lfound=.true. 
   if (list_cloud_nums > 0) then
     do n1=1,list_cloud_nums
       call find_name (n3,prof_names(1:n3,3),lstop,mysub, &
                       trim(list_cloud_names(n1)),nf_cloud_list(n1))
       if (nf_cloud_list(n1) == 0) then
         lfound=.false.
         ier=ier+1000
       endif
     enddo
     if (.not. lfound) then
       print *,'SOME REQUESTED CLOUD FIELDS NOT FOUND'
     elseif (lprint) then
       print *,'Cloud fields to be used by CRTM are: '
       print ('(10a10)'),list_cloud_names(1:list_cloud_nums)     
     endif
   endif  
!
!  Set profile indexes for any aerosol fields to be used
!  Also check that rh field present
   n3=prof_num_3d
   lfound=.true. 
   if (list_aerosol_nums > 0) then
     do n1=1,list_aerosol_nums
       call find_name (n3,prof_names(1:n3,3),lstop,mysub, &
                       trim(list_aerosol_names(n1)),nf_aerosol_list(n1))
       if (nf_aerosol_list(n1) == 0) then
         lfound=.false.
         ier=ier+1000
       endif
     enddo
     if (.not. lfound) then
       print *,'SOME REQUESTED AEROSOL FIELDS NOT FOUND'
     else
       print *,'Aerosol fields to be used by RTM are: '
       print ('(10a10)'),list_aerosol_names(1:list_aerosol_nums)     
     endif
!
     call find_name (n3,prof_names(1:n3,3),lstop,mysub,'RH',nf_rh)
     if (nf_rh == 0) then
       ier=ier+1000
       print *,'RH REQUIRED FOR USING AEROSOLS NOT FOUND'
     endif     
   endif    
!
!  Set required ids for required variables in extra info list
   n4=obs_list_extra_recs
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA' .or. &
       trim(dtype) == 'IASI' .or. trim(dtype) == 'CRIS' .or.      &
       trim(dtype) == 'CRISFSR' ) then     
     call find_name (n4,obs_list_extra_names,lstop,mysub,'CHNM',nf_chnm)
   endif
!
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then
       call find_name (n4,obs_list_extra_names,lstop,mysub,'LOGRCW',nf_spec)
   endif
!
   if (dtype == 'IASI') then
     call find_name (n4,obs_list_extra_names,lstop,mysub,'CSCALE',nf_spec)
   endif
!     
   end subroutine rad_index_setup 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_rad_index
