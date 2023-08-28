   program merge_bufr
!
! Merge several GPSRO BUFR files into one BUFR file changing: 
! times (so all have same 6 hour period)
! additional levels below those in the original file
! index in meta-data indicating which numbered file was used to create the data
! Make all obs with same satid = 745 and qc=100 percent OK
!
   use m_kinds, only : rkind1,rkind2
! 
   use m_gpsro_names, only : obs_max_levs, obs_nvalues, obs_info_num
   use m_gpsro_names, only : gpsro_names_setup
   use m_gpsro_names, only : bbnum, bbnlold, bbnlnew, bbnmsg, bbdhr, bbtslot
   use m_gpsro_names, only : bbsaid, bbyear, bbcurve, bbimpp, bbbang 
   use m_gpsro_names, only : bbptid, bbpccf
!
   use m_nr_fields_info, only : nr_fields_setup  
   use m_nr_fields_info, only : field_time_slots
   use m_nr_fields_info, only : field_time_delta, field_time_first
!
   use m_bufr_gpsro, only : gpsro_rw 
   use m_bufr_gpsro, only : gpsro_info, gpsro_values
   use m_bufr_gpsro, only : bufr_unit_in, bufr_unit_out, obs_nlevs, obs_bmiss
!
   use m_time_compute, only : time_compute_pack, time_compute_unpack
   use m_time_compute, only : time_compute_add, time_compute_dhours
   use m_time_compute, only : rkindh
!
   implicit none
!
   logical, parameter :: lsample=.true. ! print out a sample of input/output 
   logical :: lprint
   logical :: modtype
   logical :: lfound
   logical :: lerror
   logical :: ldum
!
   integer, parameter :: obs_data_max=25000   ! max # of multi-lev reports  
   integer, parameter :: nerrors=7            ! leave as is
   integer, parameter :: ndum=1               ! leave as is
   integer, allocatable :: obs_counts(:)
   integer :: ierrors(nerrors)
   integer :: idum(ndum)                 ! dummy argument array
   integer :: k,n
   integer :: obslevs_new 
   integer :: nob, nobs
   integer :: nskip_bufr
   integer :: ntimes_bufr
   integer :: ier
   integer :: ntime_b
   integer :: nread
   integer :: nmsg       ! message counter
   integer :: idate      ! date time in bufr messages 
   integer :: ireadsb,ireadmg ! reading functions in bufr lib 
   integer :: argc
   integer(4)  :: iargc
!
   real(rkindh), parameter :: t6h=6._rkindh   ! analysis time period = 6hr
!OLD   real(rkind1), parameter :: zmin=1500.      ! lowest new levvel will be roc+zmin
   real(rkind1), parameter :: zmin=1.e7  ! NO new levels
   real(rkind1), parameter :: zdel=150.       ! spacing (m) btween added new levels
   real(rkind1) :: zrmin
   real(rkindh) :: dtime_bufr
   real(rkindh) :: time_bufr(6), time_bufr_old(6), tref(6), tobs(6)
   real(rkind2) :: bmiss99
   real(rkindh) :: dhours
   real(rkind1) :: roc  ! radius of curvature as read
   real(rkind1), allocatable :: obs_info(:,:)
   real(rkind1), allocatable :: obs_data(:,:,:)
   real(rkind1), allocatable :: znew(:)
!
   character(len=8) :: subset
   character(len=1)   :: c_test       ! T or F indicates whether to print more
   character(len=14)  :: cdtime_old   ! yyyymmddhhmmss of original bufr data 
   character(len=14)  :: cdtime_new   ! yyyymmddhhmmss of new bufr data (NR time)
   character(len=14)  :: cdtime_table ! cdtime of file with non-missing BUFR table
   character(len=14)  :: cdtime
   character(len=3)   :: cntimes_bufr
   character(len=4)   :: cnskip_bufr
   character(len=240) :: field_list_file ! .rc file containing field info
   character(len=240) :: bufr_in_path    
   character(len=240) :: bufr_out_file
   character(len=240) :: bufr_file_name
   character(len=*), parameter ::  my_name='main_prog'
   character(len=*), parameter :: bufr_in_template= &
         '/Y$yyyy#/M$mm#/gdas1.$yymmdd#.t$hh#z.gpsro.tm00.bufr_d'
!
! Read arguments
   argc = iargc()
   if (lprint .and. argc /= 8) then
     print *,' usage must be: prog.x cdtime_old cdtime_new', &
             ' cntimes_bufr cnskip_bufr field_list_file bufr_in_path', &
             ' bufr_out_file c_test'
     stop
   endif
   call GetArg( 1_4, cdtime_old)     ! datetime for first input BUFR file
   call GetArg( 2_4, cdtime_new)     ! datetime for merged BUFR file to create 
   call GetArg( 3_4, cntimes_bufr)   ! number of BUFR files to merge
   call GetArg( 4_4, cnskip_bufr)    ! number of 6hr periods between input BUFR files
   call GetArg( 5_4, field_list_file) ! .rc file (needed but will not impact output)
   call GetArg( 6_4, bufr_in_path)    ! dir path for BUFR files input
   call GetArg( 7_4, bufr_out_file)   ! output BUFR file name
   call GetArg( 8_4, c_test)          ! T or F (T if sample i/o to be printed)
!
   read (cntimes_bufr,'(i3)') ntimes_bufr
   read (cnskip_bufr,'(i3)') nskip_bufr
   lprint = (c_test == 'T')
   print *,' '
   print *,'ntimes_bufr,nskip_bufr= ',ntimes_bufr,nskip_bufr
!
   call gpsro_names_setup (ier)
   allocate (znew(obs_max_levs))
   allocate (obs_counts(ntimes_bufr))
   allocate (obs_info(obs_info_num,obs_data_max))
   allocate (obs_data(obs_nvalues,obs_max_levs,obs_data_max))
!
   call nr_fields_setup ('GPSRO',field_list_file,lprint,ier)
!
   bmiss99=obs_bmiss*0.999
   dtime_bufr=real(nskip_bufr,rkindh)*t6h  
   call time_compute_unpack (cdtime_old,time_bufr_old)
   call time_compute_unpack (cdtime_new,tref)
!
! Loop over multiple BUFR files 
!
   nmsg=0
   nob=0
   obs_counts(:)=0
   do ntime_b=1,ntimes_bufr
!
     dhours=real(ntime_b-1,rkindh)*dtime_bufr
     call time_compute_add (dhours,time_bufr_old,time_bufr,ier)
     call time_compute_pack (time_bufr,cdtime)
     print *,'ntime_b,cdtime,dhours= ',ntime_b,cdtime,dhours
!
! Open input BUFR file 
     call set_field_file_name ('NONE',bufr_in_template, &
                               bufr_in_path,cdtime,bufr_file_name,ier)
     open(unit=bufr_unit_in,file=trim(bufr_file_name),form='unformatted')
     print ('(3a,i3)'),' bufr input file=',trim(bufr_file_name),     &
              ' opened on unit=',bufr_unit_in
     call openbf(bufr_unit_in,'IN ',bufr_unit_in)
!   
! Read all reports from BUFR file 
! 
     ierrors(:)=0
!
     do while (ireadmg(bufr_unit_in,subset,idate) == 0)
       nmsg=nmsg+1
       do while (ireadsb(bufr_unit_in) == 0)
!
! Get next report
         nread=nob+1
         call gpsro_rw (lprint,.true.,cdtime,nread,nerrors,ierrors,lerror,ier)
         if (ier == 0 .and. (.not. lerror)) then ! no errors detected in header
           nob=nob+1  ! count only reports of requested type with data
           nob=min(nob,obs_data_max)
           obs_counts(ntime_b)=obs_counts(ntime_b)+1
!
! Copy header info that has been passed
           obs_info(1:obs_info_num,nob)=gpsro_info(1:obs_info_num)
           obs_data(1:obs_nvalues,1:obs_nlevs,nob)= &
                     gpsro_values(1:obs_nvalues,1:obs_nlevs)
!
! Copy indexes to indicated arrays and original number of obs levels
           obs_info(bbnmsg,nob)=real(nmsg)
           obs_info(bbnum,nob)=real(nob)     ! id for obs_info array 
           obs_info(bbnlold,nob)=real(obs_nlevs)
           obs_info(11,nob)=real(ntime_b)    ! replace ptid by bufr file counter
!
           if (lsample .and. mod(nob,5000) == 1) then
             call sample (nob,obs_nlevs,'INPUT')
           endif
!
! Add levels lower down (assume data level 1 is nearest the surface)
           zrmin=zmin+obs_info(bbcurve,nob)
!
! First count the number of levels needed to fill in
           if (obs_data(bbimpp,1,nob) > zrmin) then 
             do k=1,obs_max_levs-obs_nlevs
               obslevs_new=k
               znew(k)=obs_data(bbimpp,1,nob)-k*zdel
               if (znew(k) < zrmin) exit
             enddo
!
! Shift old values up in data array to fill in lower values first            
             do k=obs_nlevs,1,-1
               n=k+obslevs_new
               obs_data(1:obs_nvalues,n,nob)=obs_data(1:obs_nvalues,k,nob)
             enddo
!
! Copy additional new values into data array in reverse order
! Set all new data vales to same corresponding values in lowest level of
! BUFR data, except for impp
             do k=obslevs_new,1,-1
               n=obslevs_new-k+1
               obs_data(1:obs_nvalues,n,nob)= &
                    obs_data(1:obs_nvalues,k,nob) 
               obs_data(bbimpp,n,nob)=znew(k) 
             enddo 
!
! Repace any missing bang values by OK value 
             obslevs_new=obs_nlevs+obslevs_new
             do k=1,obslevs_new
               if (obs_data(bbbang,k,nob) > bmiss99) then 
                 obs_data(bbbang,k,nob)=9.9999e-6
               endif
             enddo 
!
           else
             obslevs_new=obs_nlevs ! new and old num of levels same
           endif  ! test on whether additional lower levs required
           obs_info(bbnlnew,nob)=real(obslevs_new)
!
! Replace time (dhours is time relative to central time of file)
           call time_compute_dhours (time_bufr,obs_info(bbyear:bbyear+5,nob), &
                                     dhours,ier)
           call time_compute_add (dhours,tref,tobs,ier)
           obs_info(bbdhr,nob)=dhours
           obs_info(bbyear:bbyear+5,nob)=tobs(1:6)       ! new obs time   
!
! Replace SAT-ID, Transmitter-ID, quality flag
! replace transmitter id with flag with index of original bufr-file used
! (Then this index can be used to exclude sets of obs in GSI read bufr)
           obs_info(bbsaid,nob)=745  ! a COSMIC FM6 obs
           obs_info(bbpccf,nob)=100. ! 100 percent confidence Q OK
           obs_info(bbptid,nob)=real(ntime_b) 
!
! Determine which time slot the obs belongs in 
           n=int((obs_info(bbdhr,nob)-field_time_first)/field_time_delta)
           obs_info(bbtslot,nob)=min(n+1,field_time_slots-1)
!
         endif    ! check on nlevs>0
       enddo      ! while ireadsb
     enddo        ! while ireadmg
     call closbf (bufr_unit_in)   ! close input buffer file
     if (obs_counts(ntime_b) >0) then
       cdtime_table=cdtime
     endif
!
   enddo ! end loop over ntime_b
!
   nobs=nob
   print *,'nmsg,nobs=',nmsg,nobs
   print *,'obs_counts=',obs_counts
!
   if (nobs == 0) then
     print *,'NO OBS TO COPY AT ANT TIME'
!
! Open input BUFR file to copy BUFR table
   else
     call set_field_file_name ('NONE',bufr_in_template, &
                               bufr_in_path,cdtime_table,bufr_file_name,ier)
     open(unit=bufr_unit_in,file=trim(bufr_file_name),form='unformatted')
     print ('(3a,i3,a)'),' bufr input file=',trim(bufr_file_name),     &
              ' opened on unit=',bufr_unit_in,'  to read BUFR table'
     call openbf(bufr_unit_in,'IN ',bufr_unit_in)
!
! Open output BUFR file 
     print *,' '
     print *,'Attempting to output new file: nobs=',nobs
     open(unit=bufr_unit_out,file=trim(bufr_out_file),form='unformatted')
     print ('(3a,i3)'),' bufr output file=',trim(bufr_out_file),     &
              ' opened on unit=',bufr_unit_out
     call openbf(bufr_unit_out,'OUT',bufr_unit_in)
!
     call maxout(200000) 
     call datelen (10)
     read (cdtime_new,'(i10)') idate
     subset='NC003010'
!
     print *,'cdtime_new,subset',cdtime_new,' ',subset,' ',idate
    
     do nob=1,nobs
!
       if (lsample .and. mod(nob,5000) == 1) then
         obslevs_new=nint(obs_info(bbnlnew,nob))
         call sample (nob,obslevs_new,'OUTPUT')
       endif
!
       call openmb(bufr_unit_out,subset,idate)
       gpsro_info(1:obs_info_num)=obs_info(1:obs_info_num,nob)
       gpsro_values(:,:)=obs_data(:,:,nob)
       call gpsro_rw (lprint,.false.,cdtime_new,nob,ndum,idum,ldum,ier)
       call closmg (bufr_unit_out)
     enddo ! loop over obd
     call closbf (bufr_unit_out)   ! close output buffer file
   endif
!
   deallocate (znew,obs_info,obs_data,obs_counts)
!
   contains
!
   subroutine sample (nx,levs,cx)
   implicit none
   integer :: nx,kx,levs
   character(len=*) :: cx
!
   print *,' '
   print *,'Sample: ',nx,cx,levs
!
   print *,'Info: ',obs_info(:,nx)
   do kx=1,levs
     print ('(i3,1p8e20.8)'),kx,obs_data(:,kx,nx)
   enddo
!
   end subroutine sample
! 
   end program merge_bufr
