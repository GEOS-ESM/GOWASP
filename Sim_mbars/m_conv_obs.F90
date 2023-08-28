   module m_conv_obs
!
!  Module that uses Shared memory obs info and obs value arrays
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
!
   use m_kinds, only : rkind1,rkind2
!
   use m_conv_types, only: footprint_report
!
   implicit none
   include "mpif.h"
!
   private
   public :: conv_obs_setup
   public :: conv_obs_clean
   public :: conv_obs_read_bufr
   public :: conv_obs_in_tslot  
   public :: conv_obs_process
   public :: conv_obs_determine_nobs 
   public :: conv_obs_output
!
! Shared memory arrays containing obs info and data
! Obs_info contains obs header and other information
! obs_data_1lev contains original and new obs values for single-level reports 
! obs_data_mlev_old contains original obs values for multi-level reports 
! obs_data_mlev_new contains new obs values for multi-level reports
   real(rkind2), pointer :: obs_info(:,:)          => null()
   real(rkind1), pointer :: obs_data_1lev(:,:,:)   => null()
   real(rkind1), pointer :: obs_data_mlev_old(:,:,:) => null()
   real(rkind1), pointer :: obs_data_mlev_new(:,:,:) => null()
!
   integer, parameter :: obs_data_1lev_max=500000 ! max # of single-lev reports
   integer, parameter :: obs_data_mlev_max=3500   ! max # of multi-lev reports
   integer, parameter :: obs_data_max=obs_data_1lev_max+obs_data_mlev_max
   integer, public :: obs_count_all
   integer, public :: obs_in_tslot_num
   integer, public :: obs_in_tslot(obs_data_max)
   integer :: icounts(2,2,100:299)  ! (reports/obs,old/new.kx)   
   integer :: report_kind_count(2)  ! (1) for single-lev (2) for multi-lev    
!
   character(len=8) :: obs_msg_type(obs_data_max)
   character(len=*), parameter :: my_name='m_conv_obs'
!
   contains 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine conv_obs_setup (lprint,ier) 
!
! Allocate and initialize Shmem arrays to hold obs information              
!
   use m_conv_names, only : conv_max_levs, conv_nfields
   use m_conv_names, only : obs_info_num
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
!
   integer :: dim2(2)
   integer :: dim3(3)
   integer :: ierr
   character(len=*), parameter :: mysub=my_name//'::conv_obs_setup'
!
! Allocate space for the global obs info and data fields using SHMEM
   ier=0
   dim2=(/obs_info_num,obs_data_max/) 
   call MAPL_AllocNodeArray(obs_info,dim2,rc=ierr)
   if (ierr /=0) ier=ier+1
   dim3=(/conv_nfields,2,obs_data_1lev_max/)   ! (1)=old, (2)=new
   call MAPL_AllocNodeArray(obs_data_1lev,dim3,rc=ierr)
   if (ierr /=0) ier=ier+1
   dim3=(/conv_nfields,conv_max_levs,obs_data_mlev_max/) 
   call MAPL_AllocNodeArray(obs_data_mlev_old,dim3,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_AllocNodeArray(obs_data_mlev_new,dim3,rc=ierr)
   if (ierr /=0) ier=ier+1
!
   if (ier == 0) then 
     obs_info(:,:)=0._rkind2
     obs_data_1lev(:,:,:)=0._rkind1
     obs_data_mlev_old(:,:,:)=0._rkind1
     obs_data_mlev_new(:,:,:)=0._rkind1
   elseif (lprint) then
     print *,'ERROR IN ',mysub,'  ier=',ier
   endif
!
   end subroutine conv_obs_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine conv_obs_clean (lprint,ier)
!
!  Deallocate arrays created by conv_obs_setup 
!
   implicit none
   logical, intent(in) :: lprint
   integer, intent(out) :: ier
!
   integer :: ierr
   character(len=*), parameter :: mysub=my_name//'::conv_obs_clean'
!
   call MAPL_DeallocNodeArray(obs_info,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_DeallocNodeArray(obs_data_1lev,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_DeallocNodeArray(obs_data_mlev_old,rc=ierr)
   if (ierr /=0) ier=ier+1
   call MAPL_DeallocNodeArray(obs_data_mlev_new,rc=ierr)
   if (ierr /=0) ier=ier+1
!
   if (lprint .and. ier /= 0) then
     print *,'ERROR IN ',mysub,'  ier=',ier
   endif
!
   end subroutine conv_obs_clean
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine conv_obs_read_bufr (lprint,ltest_sample,dtype,ier)
!
!  Read and save all conventional observations 
!  
   use m_nr_fields_info, only : field_time_slots
   use m_nr_fields_info, only : field_time_delta, field_time_first
!
   use m_conv_names, only : conv_nfields, conv_nhead
   use m_conv_names, only : bbp, bbnall, bbnknd, bbnlold, bbpmin, bbtslot
   use m_conv_names, only : idhr, ityp, bbr, bbtmin, bbnmsg, bboslv, ielv
!
   use m_bufr_conv, only : conv_rw, bufr_unit_in
   use m_bufr_conv, only : conv_bmiss, conv_nlevs
   use m_bufr_conv, only : conv_info, conv_values
!
   implicit none
!
   logical, intent(in) :: lprint
   logical, intent(in) :: ltest_sample
   integer, intent(out) :: ier
   character(len=*), intent(in) :: dtype
!
   logical, parameter :: lread=.true.     
   logical :: lfound           ! true if an obs eport is found in msg
   logical :: lerror
   logical :: ok_obs
   logical :: lsample
   logical :: modtype
   logical :: check_type_conv  ! logical function in m_buf_conv.f90 file
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
   bmiss99=conv_bmiss*0.999
   ierrors(:)=0
   report_kind_count(:)=0
   icounts(:,:,:)=0
   nob=0
   nmsg=0
!
   do while (ireadmg(bufr_unit_in,subset,idate_old) == 0)
!
! The function check_type_conv is in the file m_bufr_conv.f90 
     modtype = check_type_conv (subset,dtype,'SIMOBS') 
     if (modtype) then
       lfound=.false. 
       do while (ireadsb(bufr_unit_in) == 0)
!
! Determine if data for next report read should be printed as a sample
! (-1 here turns this off, since a sample is instead printed at program end) 
         if (mod(nob,mods) == -1 .and. ltest_sample) then
           lsample=.true.
         else
           lsample=.false.
         endif
!
! Get next report
         nread=nob+1
         call conv_rw (lprint,lsample,reqtype,nread,lread, &
                       nerrors,ierrors,lerror,ier)
!
! Only consider obs if (a) no errors in detected in report, 
! (b) space exists in obs_info array, and (c) space exists in 
! appropriate obs_data array. 
! 
         ok_obs=.false. 
         if (ier == 0 .and. (.not. lerror) .and. nob < obs_data_max ) then 
           itype=nint(conv_info(ityp))
           lmult= conv_nlevs > 1 .or.                     &
                  (itype >= 120 .and. itype <= 129) .or.  &
                  (itype >= 220 .and. itype <= 229) 
           if ((lmult .and.                                      &
                  report_kind_count(2) < obs_data_mlev_max) .or. &
               ((.not. lmult) .and.                              &
                  report_kind_count(1) < obs_data_1lev_max)) then  
             ok_obs=.true.
           endif 
         endif 
!
! Only consider and save info and data for ok_obs that fit in arrays
         if (ok_obs) then 
           nob=nob+1  ! count only reports of requested types with OK data
!
! If this is first accepted report in msg, then increment msg counter and 
! change found flag 
           if (.not. lfound) then 
             nmsg=nmsg+1
             lfound=.true.
           endif
!
! Copy header info that has been passed
! Copy reported station elevation since the read value will later be replaced
           obs_info(1:conv_nhead,nob)=conv_info(1:conv_nhead)
           obs_info(bboslv,nob)=obs_info(ielv,nob)
!
! Copy field data to save all in memory
           if (.not. lmult) then ! save as single-level report 
             report_kind_count(1)=report_kind_count(1)+1
             obs_id=report_kind_count(1)
             obs_data_1lev(1:conv_nfields,1,obs_id)= &
                 conv_values(1:conv_nfields,1) 
             tmin=max(obs_info(idhr,nob),field_time_first)
             obs_info(bbtmin,nob)=min(tmin,field_time_first+ &
                         field_time_delta*(field_time_slots-1))
             obs_info(bbpmin,nob)=1. ! not required here for single-level obs
           else                   ! save as multi-level report   
             report_kind_count(2)=report_kind_count(2)+1
             obs_id=report_kind_count(2)
             obs_data_mlev_old(1:conv_nfields,1:conv_nlevs,obs_id)= &
                 conv_values(1:conv_nfields,1:conv_nlevs)
!
! Determine minumum pressure and time in report (needed for sondes)
! This will be the minimum of the time reported in the header (obs_info(idhr))
! and in the reported "drift" in sonde reports. When drift info is indicated 
! as missing values, the minimum will be the header value.
! 
             pmin=1.e6
             tmin=obs_info(idhr,nob) 
             do k=1,conv_nlevs
               p=obs_data_mlev_old(bbp,k,obs_id) 
               pmin=min(pmin,p)
               obs_time=obs_data_mlev_old(bbr,k,obs_id) 
               tmin=min(tmin,obs_time)
             enddo
             tmin=max(tmin,field_time_first)
             tmin=min(tmin,field_time_first+field_time_delta* &
                         (field_time_slots-1))
             obs_info(bbtmin,nob)=tmin
             obs_info(bbpmin,nob)=pmin-1. ! reduce so min. lev included
!
           endif       ! test on whether multi-level obs
!
! Copy indexes to indicated arrays and original number of obs levels
           obs_info(bbnmsg,nob)=real(nmsg)
           obs_msg_type(nmsg)=subset
           obs_info(bbnall,nob)=real(nob)     ! id for obs_info array 
           obs_info(bbnknd,nob)=real(obs_id)  ! id for approp. obs_data array
           obs_info(bbnlold,nob)=real(conv_nlevs)
           icounts(1,1,itype)=icounts(1,1,itype)+1
           icounts(2,1,itype)=icounts(2,1,itype)+conv_nlevs*2
!
! Determine initial time slot for obs
           n=int((obs_info(bbtmin,nob)-field_time_first)/field_time_delta)
           obs_info(bbtslot,nob)=min(n+1,field_time_slots-1)
!
         endif    ! check on if obs report and type OK and space available
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
       print ('(i8,a)'),ierrors(7), &
           ' observation reports found but sat or instrument ID not used' 
     endif  ! test on sum of errors
   endif    ! test on lprint
!      
   end subroutine conv_obs_read_bufr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine conv_obs_in_tslot (ntslot)
!
!  Count the numbers of obs in each time slot
!   
   use m_conv_names, only : bbtslot
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
   end subroutine conv_obs_in_tslot  
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x      
!
   subroutine conv_obs_process (nob,rtime1,rtime2,ctype,ier) 
!
!  Call appropriate routine to create each type of conventional obs.
!  See list of m_conv_types below.
!
   use m_conv_types, only : sonde_drift_wind
   use m_conv_types, only : sonde_drift_mass
   use m_conv_types, only : multi_level_report
   use m_conv_types, only : single_level_report
   use m_conv_types, only : surface_level_report
!
   use m_conv_names, only : conv_nfields, conv_max_levs
   use m_conv_names, only : obs_info_num
   use m_conv_names, only : ityp, iwsd, bbnknd, bbnlold, bbnlnew
!
   use m_nr_fields_info, only : field_kdim, field_num_3d
!   
   implicit none
!
   integer, intent(in)  :: nob           ! observation set counter
   integer, intent(out) :: ier
   real(rkind1), intent(in) :: rtime1  ! rel. time(hrs) for NR data begin tslot
   real(rkind1), intent(in) :: rtime2  ! rel. time(hrs) for NR data end tslot
   character(len=*), intent(in) :: ctype   ! 'PWD' or 'TQ'
!
   logical  lwind20m  ! need 20m wind
   integer :: itype
   integer :: klevs
   integer :: obs_id
   integer :: sub_index
   real(rkind2) :: r_sub_index
   real(rkind2) :: obs_info1(obs_info_num)
   real(rkind1) :: obs_data_n(conv_nfields,conv_max_levs)
   real(rkind1) :: obs_data_o(conv_nfields,conv_max_levs) 
!
   ier=0
!
   obs_info1(:)=obs_info(:,nob)
   obs_id=obs_info1(bbnknd)      ! index for obs_data arrays
   klevs=nint(obs_info1(bbnlold))
   itype=nint(obs_info1(ityp))
   if (itype == 280) then 
     r_sub_index=min(9999._rkind2,obs_info1(iwsd))
     sub_index=nint(r_sub_index)    ! ship type if < 9999
     lwind20m=(sub_index == 522 .or. sub_index == 523 .or. sub_index == 531)
   else
     lwind20m=.false.
   endif
!  
! Process raob, pibal, or dropsonde 
   if (mod(itype,100) == 20 .or. itype == 221 .or. mod(itype,100) == 32) then
     obs_data_n(:,:)=obs_data_mlev_new(:,:,obs_id) 
     if (ctype == 'PWD') then           ! fill in drift locations and p,z,u,v
       call sonde_drift_wind (obs_info_num,conv_nfields,field_kdim, &
                              rtime1,rtime2,obs_info1,obs_data_n)
     else                               ! fill in t,q
       call sonde_drift_mass (obs_info_num,conv_nfields,field_kdim, &
                              rtime1,rtime2,obs_info1,obs_data_n)  
     endif
     obs_data_mlev_new(:,:,obs_id)=obs_data_n(:,:)
!
! Process multi-level non-sonde report (RASS, VADWIND, PROFILER)
   elseif (mod(itype,100) >= 23  .and. mod(itype,100) <= 29) then 
     obs_data_o(:,:)=obs_data_mlev_old(:,:,obs_id)
     obs_data_n(:,:)=obs_data_mlev_new(:,:,obs_id)
     call multi_level_report (obs_info_num,conv_nfields,field_kdim,ctype, &
                              rtime1,rtime2,obs_info1,obs_data_o,obs_data_n)
     obs_data_mlev_new(:,:,obs_id)=obs_data_n(:,:)
!
! Process surface obs 
   elseif (mod(itype,100) >= 80 .and. mod(itype,100) <= 90 .and. &
           (.not. lwind20m)) then
     obs_data_o(:,1)=obs_data_1lev(:,1,obs_id)   
     obs_data_n(:,1)=obs_data_1lev(:,2,obs_id)   
     call surface_level_report (obs_info_num,conv_nfields,ctype,rtime1, &
                          rtime2,obs_info1,obs_data_o(:,1),obs_data_n(:,1))
     obs_data_1lev(:,2,obs_id)=obs_data_n(:,1)   
! process mbars
   elseif (itype == 191) then
     obs_data_o(:,1)=obs_data_1lev(:,1,obs_id)
     obs_data_n(:,1)=obs_data_1lev(:,2,obs_id)
     call footprint_report (obs_info_num,conv_nfields,ctype,rtime1, &
                  rtime2,obs_info1,obs_data_o(:,1),obs_data_n(:,1),50.)
     obs_data_1lev(:,2,obs_id)=obs_data_n(:,1)
! 
! Process single-level report that is not a surface or 10m wind report      
   elseif (klevs == 1) then 
     obs_data_o(:,1)=obs_data_1lev(:,1,obs_id)   
     obs_data_n(:,1)=obs_data_1lev(:,2,obs_id)   
     call single_level_report (obs_info_num,conv_nfields,field_kdim,ctype, &
                      lwind20m,rtime1,rtime2,obs_info1,obs_data_o(:,1),    &
                      obs_data_n(:,1))
     obs_data_1lev(:,2,obs_id)=obs_data_n(:,1)   
   endif
!
! Update header
   obs_info(:,nob)=obs_info1(:)
!
   end subroutine conv_obs_process
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine conv_obs_determine_nobs (lprint)
!
!  Determine the number of obs to be processed.
!  This is accomplished by checking how many obs_info have information.
!  In this way, all processors have the information.
!
   use m_conv_names, only : bbnall
   implicit none
!
   logical, intent(in) :: lprint
!
   integer :: k
!
   obs_count_all=0
   do k=1,obs_data_max
     if (nint(obs_info(bbnall,k)) == 0) then
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
   end subroutine conv_obs_determine_nobs
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine conv_obs_output (lprint,lsample,cdtime_new)  
!
! Finish processing sonde reports and output all obs
!
   use m_conv_names, only : conv_nfields, conv_max_levs
   use m_conv_names, only : obs_info_num
   use m_conv_names, only : ityp, bbnknd, bbnmsg, bbnlold
!
   use m_bufr_conv,  only : bufr_unit_out
!   
   implicit none
!
   logical, intent(in)  :: lprint
   logical, intent(in)  :: lsample
   character(len=*), intent(in) :: cdtime_new
!
   integer :: k
   integer :: itype
   integer :: nlevs
   integer :: obs_id
   integer :: nmsg
   integer :: lunit
   integer :: idate_new
   real(rkind2) :: obs_info1(obs_info_num)
   real(rkind1) :: obs_old(conv_nfields,conv_max_levs)
   real(rkind1) :: obs_new(conv_nfields,conv_max_levs)
   character(len=8) :: subset
!
   nmsg=0
   lunit=bufr_unit_out
   read (cdtime_new,'(i10)') idate_new
!
   do k=1,report_kind_count(1)+report_kind_count(2)
!
     if (nint(obs_info(bbnmsg,k)) /= nmsg) then
       if (nmsg /= 0) then
          call closmg(lunit)
       endif 
       nmsg=nint(obs_info(bbnmsg,k))
       subset=obs_msg_type(nmsg)
       call openmb(lunit,subset,idate_new)
     endif 
!
     itype=nint(obs_info(ityp,k))
     nlevs=nint(obs_info(bbnlold,k))
     obs_info1(:)=obs_info(:,k)
     obs_id=nint(obs_info(bbnknd,k))
!
! Output for sondes
     if (mod(itype,100) == 20 .or. mod(itype,100) == 21 .or. &
         mod(itype,100) == 32) then  
       obs_old(:,:)=obs_data_mlev_old(:,:,obs_id)
       obs_new(:,:)=obs_data_mlev_new(:,:,obs_id)
!
       call sondes_output (conv_nfields,conv_max_levs,obs_info_num,lunit, &
                           itype,k,lprint,lsample,icounts(:,:,itype), & 
                           obs_info1,obs_old,obs_new) 
!
! Output for obs other than sondes
     else                   
!
       if (nlevs > 1) then        
         obs_old(:,:)=obs_data_mlev_old(:,:,obs_id)
         obs_new(:,:)=obs_data_mlev_new(:,:,obs_id)
       else
         obs_old(:,1)=obs_data_1lev(:,1,obs_id)
         obs_new(:,1)=obs_data_1lev(:,2,obs_id)
       endif 
!
       call non_sonde_output (conv_nfields,conv_max_levs,obs_info_num,lunit, &
                           itype,k,lprint,lsample,icounts(:,:,itype), & 
                           obs_info1,obs_old,obs_new)
! 
     endif  ! test on whether obs is sonde
   enddo    ! loop over obs
!
   call closmg (lunit)   ! close message
   call closbf (lunit)   ! close output buffer file
!
   if (lprint) then 
     call summary_counts (icounts)
   endif
!
   end subroutine conv_obs_output
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   end module m_conv_obs
