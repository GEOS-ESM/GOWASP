   program create_gpsro
!
! Create simulated GPSRO bending angles observations
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
!
   use m_kinds, only : rkind1,rkind2
!
   use m_gpsro_fields, only : gpsro_fields_setup
   use m_gpsro_fields, only : gpsro_fields_read
   use m_gpsro_fields, only : gpsro_fields_tv2t
!
! Use of module that defines info for required fields
   use m_nr_fields_info, only : nr_fields_setup  
   use m_nr_fields_info, only : field_time_slots
   use m_nr_fields_info, only : field_time_delta, field_time_first
   use m_nr_fields_info, only : field_ranseed
!
   use m_time_compute, only : time_compute_new_cdtime
!
   use m_gpsro_names, only : gpsro_names_setup
   use m_gpsro_names, only : obs_nvalues, obs_nhead
   use m_gpsro_names, only : obs_max_levs, obs_info_num
!
   use m_bufr_gpsro, only : gpsro_rw_setup
!
   use m_gpsro_obs, only : gpsro_obs_setup    
   use m_gpsro_obs, only : gpsro_obs_read_bufr    
   use m_gpsro_obs, only : gpsro_obs_in_tslot    
   use m_gpsro_obs, only : gpsro_obs_process    
   use m_gpsro_obs, only : gpsro_obs_determine_nobs    
   use m_gpsro_obs, only : gpsro_obs_output    
   use m_gpsro_obs, only : obs_in_tslot_num, obs_in_tslot     
   use m_gpsro_obs, only : obs_count_all
!
   use m_gpsro_ropp, only : gpsro_ropp_setup
!
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
   include "mpif.h"
! 
   integer, parameter :: ROOT=0
   integer :: ierr          ! returned error flag
   integer :: myid          ! processor id number 0, ... ,npet
   integer :: npet          ! number of processors used
   integer :: CoresPerNode
   integer :: nobs
   integer :: n,n1,n2
   integer :: ierr_read
   integer :: ntime         ! time slot index
   integer     :: argc
   integer(4)  :: iargc
!
   logical, parameter :: lruntime=.true. ! true if read times to be printed
   logical :: lprint
   logical :: ltest_sample
!
   real(rkind1) :: rtime1  ! relative time(hrs) for NR data begin time slot
   real(rkind1) :: rtime2  ! relative time(hrs) for NR data end of time slot
!
   character(len=20)  :: dtype        ! data type name
   character(len=1)   :: c_test       ! T or F indicates whether to print more
   character(len=14)  :: cdtime_old   ! yyyymmddhhmmss of original bufr data 
   character(len=14)  :: cdtime_new   ! yyyymmddhhmmss of new bufr data (NR time)
   character(len=14)  :: cdtime1      ! yyyymmddhhmmss of begin of time slot
   character(len=14)  :: cdtime2      ! yyyymmddhhmmss of end of time slot
   character(len=240) :: field_list_file ! .rc file containing field info
   character(len=240) :: bufr_in_file    ! 
   character(len=240) :: bufr_out_file   ! 
   character(len=8)   :: run_date(2)     ! run date from system clock
   character(len=10)  :: run_time(2)     ! run time from system clock
   character(len=*), parameter :: my_name='main_program'
!                                       ---
!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'
   die_proc_id=myid 
!
!  ---------------------------------------------------------
!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   call MAPL_InitializeShmem(rc=ierr)
!
!  Initialize other variables
   if (myid == ROOT) then
     lprint=.true.
   else
     lprint=.false.
   endif   
!
! Read arguments
   argc = iargc()
   if (lprint .and. argc /= 7) then
     print *,' usage must be: prog.x dtype cdtime_old cdtime_new', & 
             ' field_list_file bufr_in_file bufr_out_file c_test'
     call mpi_die (my_name,77)
   endif
   call GetArg( 1_4, dtype)
   call GetArg( 2_4, cdtime_old)
   call GetArg( 3_4, cdtime_new)
   call GetArg( 4_4, field_list_file)
   call GetArg( 5_4, bufr_in_file)
   call GetArg( 6_4, bufr_out_file)
   call GetArg( 7_4, c_test)
!
! ltest_sample is true if a sample of output is to be printed for testing
   if (lprint .and. c_test == 'T') then
     ltest_sample=.true.
   else
     ltest_sample=.false.
   endif
!
!  Get the field requirement info (field names and file templates)
   call nr_fields_setup (dtype,field_list_file,lprint,ierr)
   if (lprint .and. ierr /= 0) then
     print *,'Error detected in call to gpsro_fields_setup: ierr=',ierr
     call mpi_die (my_name,ierr)
   endif
!
! Setup names associated with obs header and data fields
   call gpsro_names_setup (ierr)
   if (lprint .and. ierr /= 0) then
     print *,'Error detected in call to gpsro_names_setup: ierr=',ierr
     call mpi_die (my_name,ierr)
   endif
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
! Allocate space for the global NR fields using SHMEM
   call gpsro_fields_setup
   call MPI_Barrier(MPI_COMM_WORLD,ierr)             
!
! Allocate space for the global obs info and obs value arrays using SHMEM
   call gpsro_obs_setup (lprint,ierr)
   if (ierr /= 0) then
     call mpi_die (my_name,73)
   endif  
!
!  Allocate space for some ROPP data structures
   call gpsro_ropp_setup (1)
!
   if (myid == ROOT) then
!
! Open bufr files (the input one is used to get the bufr table)
     call gpsro_rw_setup (bufr_in_file,bufr_out_file,lprint,'SIMOBS',ierr)     
     if (ierr /= 0) then
       call mpi_die (my_name,ierr)
     endif  
     call gpsro_obs_read_bufr (lprint,cdtime_old,dtype,ierr)
!
   endif  ! check on ROOT   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
! Determine number of obs (=obs_count_all) to be considered 
   call gpsro_obs_determine_nobs (lprint)  
   if (obs_count_all == 0) then
     call mpi_die (my_name,72)
   endif
!
!  Read phis (cdtime_old not used here)
   call date_and_time (date=run_date(1),time=run_time(1))  
   call gpsro_fields_read (myid,cdtime_new,cdtime_old,'PHIS',ierr_read) 
   if (ierr_read /= 0) then
     call mpi_die (my_name,ierr_read)
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                        
   call date_and_time (date=run_date(2),time=run_time(2))    
   if (lprint) then
     print *,'PHIS field read'
     if (lruntime) print *,'run times: ',run_time(1),' ',run_time(2)
   endif
!
!  Loop over time intervals in the data assimilation period
   do ntime=1,field_time_slots-1      
!
! Determine indexes and count of number of obs in tslot
     call gpsro_obs_in_tslot (ntime)
     if (obs_in_tslot_num /= 0) then
!   
!  Determine time of field data for this time
!  cdtime1 is the date time in format 'yyyymmddhhmmss'
!  rtime1 is the time relative to the central time of the assimilation period
       rtime1=(ntime-1)*field_time_delta+field_time_first
       rtime2=rtime1+field_time_delta
       call time_compute_new_cdtime (cdtime_new,cdtime1,rtime1,ierr)
       call time_compute_new_cdtime (cdtime_new,cdtime2,rtime2,ierr)
!
! Read first set of fields for this time
       call date_and_time (date=run_date(1),time=run_time(1))    
       call gpsro_fields_read (myid,cdtime1,cdtime2,'PTQ',ierr_read)
       if (ierr_read /= 0) then
         call mpi_die (my_name,ierr_read)
       endif 
       call MPI_Barrier(MPI_COMM_WORLD,ierr) 
!  
! Print CPU time to read files
       call date_and_time (date=run_date(2),time=run_time(2))       
       if (lprint) then
         print *,'PS, T, and QV fields read for time1=',cdtime1
         if (lruntime) print *,'run times: ',run_time(1),' ',run_time(2)
       endif
!
! If virtual temperature read, replace by temperature
! ierr is used as a flag here, not an error indicator
       if (myid == ROOT) then                       
         call gpsro_fields_tv2t (ierr) 
         if (ierr > 0) then  
           print *,'TV field read converted to T' 
         endif
       endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr) 
!
! Process obs in present time slot
       do n1=1,obs_in_tslot_num,npet
         n2=min(n1+npet-1,obs_in_tslot_num) 
         do n=n1,n2
           if (myid == n-n1) then 
             nobs=obs_in_tslot(n)
             call gpsro_obs_process (nobs,lprint,rtime1,rtime2,ierr)
             if (ierr /= 0) then
               print *,'ERROR in call to gpsro_obs_process from main: ', &
                        'myid=',myid
               call mpi_die (my_name,ierr)
             endif
           endif
         enddo
         call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
       enddo       
!
     endif  ! check if any obs are in time slot
   enddo    ! loop over time
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
!  Write saved obs_info
!
   if (myid == ROOT ) then
     call gpsro_obs_output (lprint,ltest_sample,cdtime_new)
   endif 
   call MPI_Barrier(MPI_COMM_WORLD,ierr)           
!
! De-allocate arrays
   call shutdown (lprint)
!
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine shutdown (lprnt)
!
! Deallocate arrays and terminate MPI
!
   use m_gpsro_fields, only : gpsro_fields_clean
   use m_gpsro_obs, only    : gpsro_obs_clean
   use m_gpsro_ropp, only   : gpsro_ropp_clean
!
   implicit none
   logical, intent(in) :: lprnt
   integer :: ier
!
   call gpsro_fields_clean (ier)
   call gpsro_obs_clean (lprnt,ier)
   call gpsro_ropp_clean
!
   call MAPL_FinalizeShmem (rc=ier)
   call MPI_Finalize (ier)
!
   end subroutine shutdown
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   end program create_gpsro
