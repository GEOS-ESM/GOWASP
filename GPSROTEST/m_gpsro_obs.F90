   module m_gpsro_obs
!
!  module that uses Shared memory obs info and obs value arrays
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
!
   use m_kinds, only : rkind1,rkind2
   use m_counter  !QQQQQ
!
   implicit none
   include "mpif.h"
!
   private
   public :: gpsro_obs_setup
   public :: gpsro_obs_clean
   public :: gpsro_obs_read_bufr
   public :: gpsro_obs_in_tslot  
   public :: gpsro_obs_process
   public :: gpsro_obs_determine_nobs 
   public :: gpsro_obs_output
!
! Shared memory arrays containing obs info and data
! obs_info contains obs header and other information
! obs_data_old contains original obs values 
! obs_data_new contains new obs values 
   real(rkind2), pointer :: obs_info(:,:)            => null()
   real(rkind1), pointer :: obs_data_old(:,:,:) => null()
   real(rkind1), pointer :: obs_data_new(:,:,:) => null()
!
   integer, parameter :: obs_data_max=10000   ! max # of multi-lev reports
   integer, public :: obs_count_all
   integer, public :: obs_in_tslot_num
   integer, public :: obs_in_tslot(obs_data_max)
!
   character(len=8) :: obs_msg_type(obs_data_max)
   character(len=*), parameter :: my_name='m_gpsro_obs'
!
   contains 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_obs_setup (lprint,ier) 
!              
! Allocate space for the global obs info and data fields using SHMEM
!
   use m_gpsro_names, only : obs_max_levs, obs_nvalues, obs_info_num
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
!
   integer :: dim2(2)
   integer :: dim3(3)
   integer :: ierr
   character(len=*), parameter :: mysub=my_name//'::gpsro_obs_setup'
!
   ier=0
   dim2=(/obs_info_num,obs_data_max/) 
   call MAPL_AllocNodeArray(obs_info,dim2,rc=ierr)
   if (ierr /=0) ier=ier+1
   dim3=(/obs_nvalues,obs_max_levs,obs_data_max/) 
   call MAPL_AllocNodeArray(obs_data_old,dim3,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_AllocNodeArray(obs_data_new,dim3,rc=ierr)
   if (ierr /=0) ier=ier+1
!
   if (ier ==0) then 
     obs_info(:,:)=0._rkind2
     obs_data_old(:,:,:)=0._rkind1
     obs_data_new(:,:,:)=0._rkind1
   elseif (lprint) then
     print *,'ERROR IN ',mysub,'  ier=',ier
   endif
!
   end subroutine gpsro_obs_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_obs_clean (lprint,ier)
!
! Deallocate Shmem arrays holding observations
!
   implicit none
   logical, intent(in) :: lprint
   integer, intent(out) :: ier
!
   integer :: ierr
   character(len=*), parameter :: mysub=my_name//'::gpsro_obs_clean'
!
   call MAPL_DeallocNodeArray(obs_info,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_DeallocNodeArray(obs_data_old,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_DeallocNodeArray(obs_data_new,rc=ierr)
   if (ierr /=0) ier=ier+1
!
   if (lprint .and. ier /= 0) then
     print *,'ERROR IN ',mysub,'  ier=',ier
   endif
!
   end subroutine gpsro_obs_clean
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine gpsro_obs_read_bufr (lprint,cdtime,dtype,ier)
!  
! Read all observations in BUFR file. Determine relative time of obs
!
   use m_nr_fields_info, only : field_time_slots
   use m_nr_fields_info, only : field_time_delta, field_time_first
!
   use m_gpsro_names, only : obs_nvalues, obs_info_num
   use m_gpsro_names, only : bbnum, bbnlold, bbnlnew, bbnmsg, bbdhr, bbtslot
!
   use m_bufr_gpsro, only : gpsro_rw
   use m_bufr_gpsro, only : obs_bmiss, obs_nlevs, bufr_unit_in
   use m_bufr_gpsro, only : gpsro_info, gpsro_values
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(out) :: ier
   character(len=*), intent(in) :: cdtime
   character(len=*), intent(in) :: dtype
!
   logical, parameter :: lread=.true.     
   logical :: lfound           ! true if an obs eport is found in msg
   logical :: lerror
   logical :: modtype
   logical :: check_type_gpsro ! logical function in m_buf_conv.f90 file
   logical :: lmult            ! obs with multi-levels possible
!
   integer, parameter :: nerrors=7
   integer, parameter :: mods=300 ! used for printing a sample of input reports
   integer :: ierrors(nerrors)
   integer :: nob
   integer :: obs_id     ! next index for either single- or multi-lev report
   integer :: itype
   integer :: k,n
   integer :: nmsg       ! message counter
   integer :: idate_old  ! date time in messages in bufr file read        
   integer :: nread      ! to pass urrent obs id before it is augmented
   integer :: ireadsb,ireadmg ! reading functions in bufr lib 
! 
   real(rkind1) :: bmiss99
   real(rkind1) :: pmin, p
   real(rkind1) :: tmin, obs_time
   character(len=*), parameter :: reqtype='EITHER' ! requested type  
   character(len=8) :: subset
!     
   bmiss99=obs_bmiss*0.999
   ierrors(:)=0
   nob=0
   nmsg=0
!
   do while (ireadmg(bufr_unit_in,subset,idate_old) == 0)
!
! The function check_type_gpsro is in the file m_bufr_gpsro.f90 
     modtype = check_type_gpsro (subset,dtype) 
     if (modtype) then
       lfound=.false. 
       do while (ireadsb(bufr_unit_in) == 0)
!
! Get next report
         nread=nob+1
         call gpsro_rw (lprint,lread,cdtime,nread,nerrors,ierrors,lerror,ier)
         if (ier == 0 .and. (.not. lerror)) then ! no errors detected in header
           nob=nob+1  ! count only reports of requested type with data
!
! If this is first accepted report in msg, then increment msg counter and 
! change found flag 
           if (.not. lfound) then 
             nmsg=nmsg+1
             lfound=.true.
           endif
!
! Copy header info that has been passed
           obs_info(1:obs_info_num,nob)=gpsro_info(1:obs_info_num)
           obs_data_old(1:obs_nvalues,1:obs_nlevs,nob)= & 
              gpsro_values(1:obs_nvalues,1:obs_nlevs) 
!
! Copy indexes to indicated arrays and original number of obs levels
           obs_info(bbnmsg,nob)=real(nmsg)
           obs_msg_type(nmsg)=subset
           obs_info(bbnum,nob)=real(nob)     ! id for obs_info array 
           obs_info(bbnlold,nob)=real(obs_nlevs)
!
! Determine initial time slot for obs
           n=int((obs_info(bbdhr,nob)-field_time_first)/field_time_delta)
           obs_info(bbtslot,nob)=min(n+1,field_time_slots-1)
!
         endif    ! check on nlevs>0
       enddo      ! while ireadsb
     endif        ! check on modtype
   enddo          ! while ireadmg
!
   if (lprint) then
     print *,'Number of obs locations of requested types found =',nob
     print *,'Number of messages for requested types found =',nmsg
!
! Print out error count summary
     if (sum(ierrors(1:nerrors)) > 0) then
       print *,' '
       print *,' Numbers of bad reports or other errors detected:'
       print ('(i8,a)'),ierrors(1), &
           ' observation reports with missing header data' 
       print ('(i8,a)'),ierrors(2), &
            ' observation reports found where time<tmin'
       print ('(i8,a)'),ierrors(3), &
            ' observation reports found where time>tmax' 
       print ('(i8,a)'),ierrors(4), &
            ' observation reports found where longitude out of range' 
       print ('(i8,a)'),ierrors(5), &
            ' observation reports found where latitude out of range' 
       print ('(i8,a)'),ierrors(6), &
            ' observation reports with obs on 0 levels' 
     endif  ! test on sum of errors
   endif    ! test on lprint
!      
   end subroutine gpsro_obs_read_bufr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine gpsro_obs_in_tslot (ntslot)
!   
!  Count number of obs in time slot and set indexes of obs in each time slot
!
   use m_gpsro_names, only : bbtslot
!
   implicit none
!
   integer, intent(in) :: ntslot
   integer :: n

   obs_in_tslot_num=0
   do n=1,obs_count_all
     if (nint(obs_info(bbtslot,n)) == ntslot) then
       obs_in_tslot_num=obs_in_tslot_num+1
       obs_in_tslot(obs_in_tslot_num)=n
     endif
   enddo
!
   end subroutine gpsro_obs_in_tslot  
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x      
!
   subroutine gpsro_obs_process (nob,lprint,rtime1,rtime2,ier) 
!
! Call ropp-sequence for each obs  
!
   use m_gpsro_names, only : obs_info_num, obs_nvalues, obs_max_levs
   use m_gpsro_names, only : bbnlnew, bbnlold, bbimpp, bbbang
   use m_gpsro_names, only : bblat, bblon, bbbaz
   use m_gpsro_names, only : bblatk, bblonk, bbbazk
   use m_gpsro_ropp, only  : gpsro_ropp_sequence
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(in)  :: nob           ! observation set counter
   integer, intent(out) :: ier
   real(rkind1), intent(in) :: rtime1  ! rel. time(hrs) for NR data begin tslot
   real(rkind1), intent(in) :: rtime2  ! rel. time(hrs) for NR data end tslot
!
   integer :: ierr_azim
   integer :: azim_flag
   integer :: ilevs
   integer :: i
   integer :: j1  ! always = 1 
   real(rkind2) :: obs_info1(obs_info_num)
   real(rkind1) :: obs_values1(obs_nvalues,obs_max_levs)
!
   ier=0
   obs_info(bbnlnew,nob)=obs_info(bbnlold,nob)
   obs_info1(:)=obs_info(:,nob)
   ilevs=nint(obs_info1(bbnlnew))
   obs_info1(bbnlold)=1
   obs_info1(bbnlnew)=1
!
! First determine if azimuth angle missing from most obs in a profile, 
! in which case the obs will be treated useing a 1-d approach at all levels
   ierr_azim=0
   do i=1,ilevs
     if (abs(obs_data_old(bbbazk,i,nob)) >= 360.) then 
       ierr_azim=ierr_azim+1
     endif
   enddo
   if (ierr_azim > ilevs/3) then 
     azim_flag=1  ! flag means treat profile as 1-d at all levels
   else
     azim_flag=2  ! flag means treat profile as 2-d in lower trop.
   endif 
!
   if (m_counter_test > 0) azim_flag=1 ! QQQQSS1
!
   do i=1,ilevs
     obs_info1(bblat)=obs_data_old(bblatk,i,nob)
     obs_info1(bblon)=obs_data_old(bblonk,i,nob)
     obs_info1(bbbaz)=obs_data_old(bbbazk,i,nob)
     obs_data_new(bbimpp,i,nob)=obs_data_old(bbimpp,i,nob)
     obs_data_new(bblatk,i,nob)=obs_data_old(bblatk,i,nob)
     obs_data_new(bblonk,i,nob)=obs_data_old(bblonk,i,nob)
     obs_data_new(bbbang,i,nob)=1.e11  ! init to missing value
     obs_values1(bbimpp,1)=obs_data_old(bbimpp,i,nob)
!
! If azim flag indicates that profile has viewing azimith defined at most 
! levels, then use the given azim, otherwise replace by 90. but do 1-d
! calculation (The azim must still be defined since the code will still
! attempt to produce a plane of values, but in the 1-d case will only use 
! the center of that plane).
     if (azim_flag == 2) then
       obs_data_new(bbbazk,i,nob)=obs_data_old(bbbazk,i,nob)
     else
       obs_data_new(bbbazk,i,nob)=90.
     endif
     obs_info1(bbbaz)=obs_data_new(bbbazk,i,nob)
!
     j1=1
     if (obs_data_old(bbbang,i,nob) < 9.e10 .and. &
         obs_data_old(bbimpp,i,nob) < 9.e10 .and. &
         abs(obs_info1(bblat)) <= 90.       .and. &
         abs(obs_info1(bblon)) < 360.       .and. &
         abs(obs_info1(bbbaz)) < 360.) then 
       call gpsro_ropp_sequence (obs_info1,obs_values1,lprint, &
                                 azim_flag,rtime1,rtime2,j1,ier)
       obs_data_new(bbbang,i,nob)=obs_values1(bbbang,1)
     endif
!
   enddo 
!
   end subroutine gpsro_obs_process
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine gpsro_obs_determine_nobs (lprint)
!
!  Determine the number of obs to be processed.
!  This is accomplished by checking how many obs_info have information.
!  In this way, all processors have the information.
!
   use m_gpsro_names, only : bbnum
   implicit none
!
   logical, intent(in) :: lprint
!
   integer :: k
!
   obs_count_all=0
   do k=1,obs_data_max
     if (nint(obs_info(bbnum,k)) == 0) then
       exit
     else
       obs_count_all=k
     endif
   enddo
!
   if (obs_count_all == 0 .and. lprint) then
     print *,'NO OBS FOUND TO PROCESS IN BUFR INPUT FILE'
   endif
!
   end subroutine gpsro_obs_determine_nobs
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine gpsro_obs_output (lprint,lsample,cdtime)  
!
! Output all GPSRO obs
!
   use m_gpsro_names, only : obs_info_num, obs_nvalues, obs_max_levs
   use m_gpsro_names, only : obs_info_names
   use m_gpsro_names, only : bbnmsg, bbnlold, bbnlnew
   use m_gpsro_names, only : bbimpp, bbbang
   use m_gpsro_names, only : bblatk, bblonk, bbbazk
!
   use m_bufr_gpsro, only : gpsro_rw
   use m_bufr_gpsro, only : bufr_unit_out, obs_nlevs, obs_bmiss
   use m_bufr_gpsro, only : gpsro_info, gpsro_values
!   
   implicit none
!
   logical, intent(in) :: lprint
   logical, intent(in) :: lsample
   character(len=*), intent(in) :: cdtime
!
   logical, parameter :: lread=.false.   ! .false. means write bufr data
   logical :: ldum                       ! argument not used here   
   integer, parameter :: ndum=1          ! dimension of dummy argument array
   integer, parameter :: isample=4       ! print isample of each obs itype  
   integer :: i,k,nob
   integer :: nlevs_old
   integer :: obs_id
   integer :: nmsg
   integer :: ier
   integer :: ilook
   integer :: lunit
   integer :: idate
   integer :: idum(ndum)                 ! dummy argument array
   integer :: nsum(2)                    ! numbers of old/new obs values  
   character(len=8) :: subset
!
   nmsg=0
   lunit=bufr_unit_out
   ilook=max(obs_count_all/isample,1)
   nsum(:)=0
   read(cdtime,'(i10)') idate
!
   do nob=1,obs_count_all
!
     if (nint(obs_info(bbnmsg,nob)) /= nmsg) then
       if (nmsg /= 0) then
          call closmg(lunit)
       endif 
       nmsg=nint(obs_info(bbnmsg,nob))
       subset=obs_msg_type(nmsg)
       call openmb(lunit,subset,idate)
     endif 
!
     nlevs_old=nint(obs_info(bbnlold,nob))
     obs_nlevs=nint(obs_info(bbnlnew,nob))
     nsum(1)=nsum(1)+nlevs_old
     nsum(2)=nsum(2)+obs_nlevs
     gpsro_info(1:obs_info_num)=obs_info(1:obs_info_num,nob)
!
     do k=1,min(nlevs_old,obs_nlevs)
       if (obs_data_old(bbbang,k,nob) > 0.9*obs_bmiss) then
         obs_data_new(bbbang,k,nob)=obs_bmiss
       endif  
     enddo
!
     gpsro_values(:,:)=obs_data_new(:,:,nob)
!
     call gpsro_rw (lprint,lread,cdtime,nob,ndum,idum,ldum,ier)
     obs_info(1:obs_info_num,nob)=gpsro_info(1:obs_info_num)
!
     if (lsample .and. mod(nob,ilook) == 1) then
!
       print *,' '
       print ('(2a,i7)'),'Sample of input and output obs to be', &
                            ' printed for nob=',nob
       print ('(a)'),'obs_info (after any changes) ='
       print ('(4(a10,1pe13.4))'), &
            (obs_info_names(i),obs_info(i,nob),i=1,obs_info_num)
!
       print ('(3a)'),'   k', &
                ' OLD: lat     lon    bear       impact       bangle', &   
                ' NEW: lat     lon    bear       impact       bangle'    
       do k=1,max(obs_nlevs,nlevs_old) 
         if (k <= nlevs_old) then 
           if (k <= obs_nlevs) then   
             print ('(i4,2(3x,f6.2,2f8.2,e13.6,e13.5))'),k,            &
                obs_data_old(bblatk,k,nob),obs_data_old(bblonk,k,nob), &
                obs_data_old(bbbazk,k,nob),                            &
                obs_data_old(bbimpp,k,nob),obs_data_old(bbbang,k,nob), &
                obs_data_new(bblatk,k,nob),obs_data_new(bblonk,k,nob), &
                obs_data_new(bbbazk,k,nob),                            &
                obs_data_new(bbimpp,k,nob),obs_data_new(bbbang,k,nob)  
           else
             print ('(i4,3x,f6.2,2f8.2,e13.6,e13.5)'),k,               &
                obs_data_old(bblatk,k,nob),obs_data_old(bblonk,k,nob), &
                obs_data_old(bbbazk,k,nob),                            &
                obs_data_old(bbimpp,k,nob),obs_data_old(bbbang,k,nob)
           endif
         else
           print ('(i4,53x,f6.2,2f8.2,e13.6,e13.5)'),k,              &
              obs_data_new(bblatk,k,nob),obs_data_new(bblonk,k,nob), &
              obs_data_new(bbbazk,k,nob),                            &
              obs_data_new(bbimpp,k,nob),obs_data_new(bbbang,k,nob)
         endif   
       enddo
     endif  ! check on whether to print sample
! 
   enddo    ! loop over obs
!
   call closmg (lunit)   ! close message
   call closbf (lunit)   ! close output buffer file
!
   if (lprint) then 
     print *,' '
     print *,'Number of obs reports =',obs_count_all
     print *,'Number of obs messages =',nmsg
     print *,'Numbers of old and new obs =',nsum(:)
   endif
!
   end subroutine gpsro_obs_output
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   end module m_gpsro_obs
