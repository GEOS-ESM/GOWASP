  program create_rad_profs
!
! This program extracts atmospheric profiles from netcdf files of NR
! fields for subsequent calculations of radiances. Mutiple types (classes) 
! of radiance observations are handeled in a single execution. It uses shared
! memory arrays created by MAPL_ShmemMod (i.e., arrays distributed
! among all processors and accessible by all using standard FORTRAN 
! without explcitly invoking MPI calls.  Files of profiles along
! with observation header information are output to a binary file.
! 
! (Shmem framework provided through Atanas Trayanov and Arlinda DaSilva)
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
   use netcdf           ! for reading the NR files
!
   use m_kinds, only : rkind1
!
! Use of module that defines info for required fields
   use m_nr_fields_info, only : nr_fields_setup
   use m_nr_fields_info, only : field_imax, field_jmax, field_kmax
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_lon_first, field_common_path
   use m_nr_fields_info, only : field_names, field_types, field_files
   use m_nr_fields_info, only : field_tavg_offset
! 
! Use of module that defines the observation data set
   use m_obs_list, only : obs_list_setup
   use m_obs_list, only : obs_list_read_recs
   use m_obs_list, only : obs_list_print_info
   use m_obs_list, only : obs_list_clean
   use m_obs_list, only : obs_list_types_put_setup
   use m_obs_list, only : obs_list_types_get_setup
   use m_obs_list, only : obs_list_types_put_recs
   use m_obs_list, only : obs_list_types_get_recs
   use m_obs_list, only : obs_list_types_rec_ids
   use m_obs_list, only : obs_list_types_allocate
   use m_obs_list, only : obs_list_types_clean
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c
   use m_obs_list, only : obs_list_i, obs_list_r, obs_list_c
   use m_obs_list, only : obs_list_tslots, obs_list_tslots_recs
   use m_obs_list, only : obs_list_time_delta, obs_list_time_first
   use m_obs_list, only : obs_list_id_time, obs_list_id_lat, obs_list_id_lon   
   use m_obs_list, only : obs_list_unit_in
!
! Use of module that writes out profiles at observation locations
   use m_write_profiles, only : write_profiles_header
   use m_write_profiles, only : write_profiles_recs
   use m_write_profiles, only : write_profiles_sample
   use m_write_profiles, only : write_profiles_tv2t
   use m_write_profiles, only : write_profiles_o3units
   use m_write_profiles, only : write_profiles_close
!
! Use module for computing time increments
   use m_time_compute, only : time_compute_new_cdtime
!
! Use modeule for assigning unit numbers
   use m_set_unit_nums, only : un_list_in, un_list_in0
   use m_set_unit_nums, only : un_prof_out, un_prof_out0
!
! Use module for printing info when aborting jobs
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
   include "mpif.h"
!
!  Global arrays to be allocated using SHMEM
!  ---------------------------------------------
   real(rkind1), pointer :: field1(:,:,:) => null()
   real(rkind1), pointer :: field2(:,:,:) => null()
   real(rkind1), pointer :: field3(:,:,:) => null()
   real(rkind1), pointer :: field4(:,:,:) => null()
   real(rkind1), pointer :: profiles(:,:) => null()
   real(rkind1), pointer :: rec_list_r(:,:) => null()
!
!  Array of profile information to output
   integer :: profiles_1dim
   integer :: profiles_2dim
   real(rkind1), allocatable :: profiles_1fld(:)
   real(rkind1), allocatable :: profiles_1loc(:)
 
!  Miscellaneous
!  -------------
   integer, parameter :: ntypes_max=20    ! max number of dtypes allowed 
   integer, parameter :: rec_ids_dim3=100 ! max number of time slots allowed
   integer :: ntypes        ! number of dtypes requested
   integer :: ntype         ! index for dtype requested
   integer :: i
   integer :: ierr          ! error flag
   integer :: ierr_read     ! error flag when reading fields
   integer :: myid          ! processor id number 0, ... ,npet
   integer :: npet          ! number of processors used
   integer :: CoresPerNode
   integer :: obs_list_set  
   integer :: obs_num       ! obs counter in current time slot
   integer :: nlist      
   integer :: ntime         ! time slot index
   integer :: nfield1
   integer :: nfield2
   integer :: nread   
   integer :: nfid1   
   integer :: nfid2
   integer :: nfid3
   integer :: nfid4
   integer :: nfid5
   integer :: profiles_dim1
   integer :: profiles_dim2
   integer :: tot_num_obs(ntypes_max)
   integer :: itv
   integer :: rec_id
   integer :: rec_r_dim1
   integer :: rec_ids_max
   integer :: rec_ids(2,ntypes_max,rec_ids_dim3)
   integer :: dim2(2)
   integer :: dim3(3)
   integer :: argc
   integer(4) :: iargc
!
   logical :: print_sample
   logical :: lprint
!
   real(rkind1) :: rtime1  ! relative time(hrs) for NR data begin time slot
   real(rkind1) :: rtime2  ! relative time(hrs) for NR data end of time slot
   real(rkind1) :: rtime3  ! offset hours needed to define tavg file times
   real(rkind1) :: obs_time
   real(rkind1) :: time_diff(2) ! time of obs compared to times bounding t slot
   real(rkind1), parameter :: bmiss=1.e10_rkind1
!
   character(len=440) :: dtype  ! data type name (reset to 'ALL' if muliple)
   character(len=20)  :: dtypes(ntypes_max)  ! individual data type names
   character(len=1)   :: c_test       ! T or F indicates whether to print more
   character(len=14)  :: cdtime_ref   ! yyyymmddhhmmss of center of assim window
   character(len=14)  :: cdtime1(2)   ! yyyymmddhhmmss of begin of time slot
   character(len=14)  :: cdtime2(2)   ! yyyymmddhhmmss of end of time slot
   character(len=240) :: obs_list_file_dir           ! dir with input obs info 
   character(len=240) :: obs_list_files(ntypes_max)  ! files of input obs info 
   character(len=240) :: field_list_file             ! rc file with field info
   character(len=240) :: profile_file_dir          ! dir of profiles to output
   character(len=240) :: profile_files(ntypes_max) ! files of profiles to output
   character(len=*), parameter :: my_name='main_program'
!                                       ---
!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   die_proc_id=myid
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'
!
!  ---------------------------------------------------------
!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   if (ierr /= 0) call mpi_die ('main:CoresPerNode',ierr)
   call MAPL_InitializeShmem(rc=ierr)
   if (ierr /= 0) call mpi_die ('main:InitializeShmem',ierr)
!
!  Only ask for printing when myid=0
   if (myid == 0) then
     lprint=.true. 
   else  
     lprint=.false.
   endif   
!
! Read arguments
   argc = iargc()
   if (argc /= 6) then
     if (lprint) then
       print *,' usage must be: create_rad_profs.x dtype cdtime_ref', & 
               ' obs_list_file_dir field_list_file profile_file_dir c_test'
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
     call mpi_die (my_name,77)
   endif
   call GetArg( 1_4, dtype)
   call GetArg( 2_4, cdtime_ref)
   call GetArg( 3_4, obs_list_file_dir)
   call GetArg( 4_4, field_list_file)
   call GetArg( 5_4, profile_file_dir)
   call GetArg( 6_4, c_test)
!
! Indicate that some extra info will be printed
   if (lprint .and. c_test == 'T') then
     print_sample=.true.
   else 
     print_sample=.false. 
   endif
!
! Determine names of i/o files for each data type requested
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   call prof_io_files (ntypes_max,dtype,obs_list_file_dir,profile_file_dir, & 
                       cdtime_ref,lprint,ntypes,obs_list_files,profile_files, &
                       dtypes,ierr)
   call obs_list_types_allocate (ntypes)
   if (ierr > 0) then
       call mpi_die (my_name,78)
   endif
!
!  Get the obs_list header info (from file of obs locations, time, etc)
   do ntype=1,ntypes
     tot_num_obs(ntype)=0 
     obs_list_unit_in=un_list_in0+ntype
     call obs_list_setup (obs_list_files(ntype),lprint,ierr)
     call obs_list_types_put_setup (ntype)
     if (lprint) then
       print *,' '
       print ('(3a,i2)'),'dtype=',trim(dtypes(ntype)),'  ntype=',ntype
       call obs_list_print_info(my_name)
     endif 
     call obs_list_clean
   enddo 
   call obs_list_types_rec_ids (ntypes,rec_ids_dim3,ntypes_max, &
                                rec_r_dim1,rec_ids_max,rec_ids)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
!  Get the field requirement info (field names and file templates)
   call nr_fields_setup (dtype,field_list_file,lprint,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to fields_read_rad_rc: ierr=',ierr
     call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
     call mpi_die (my_name,ierr)
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  Write header for output file of profiles
   if (myid == 0) then 
     print *,' '
     do ntype=1,ntypes
       un_prof_out=un_prof_out0+ntype
       call obs_list_types_get_setup (ntype)
       call write_profiles_header (profile_files(ntype),0,0)
       call obs_list_clean
       print ('(3a,i2,a,i2,a)'),'Profile header records written for dtype=', &
                                trim(dtypes(ntype)),'  ntype=',ntype,        &
                                '  unit=',un_prof_out,' and file:'
       print ('(a)'),trim(profile_files(ntype)) 
     enddo
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  Allocate local arrays to store profiles
   profiles_dim1=field_num_2d+field_num_3d*field_kmax
   profiles_dim2=rec_ids_max
   allocate (profiles_1fld(field_kmax)) 
   allocate (profiles_1loc(profiles_dim1)) 
!
!  Allocate space for the global fields using SHMEM
   dim3=(/field_imax,field_jmax,field_kmax/)
   call MAPL_AllocNodeArray(field1,dim3,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_1',ierr)
   call MAPL_AllocNodeArray(field2,dim3,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_2',ierr)
   call MAPL_AllocNodeArray(field3,dim3,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_3',ierr)
   call MAPL_AllocNodeArray(field4,dim3,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_4',ierr)
!
   dim2=(/profiles_dim1,profiles_dim2/)
   call MAPL_AllocNodeArray(profiles,dim2,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_5',ierr)
!
   dim2=(/rec_r_dim1,profiles_dim2/)
   call MAPL_AllocNodeArray(rec_list_r,dim2,rc=ierr)
   if (ierr /= 0) call mpi_die ('main:AllocNodeArray_6',ierr)
! 
!  Loop over time intervals in the data assimilation period
   call obs_list_types_get_setup (1) ! need obs_list_tslots, etc 
   do ntime=1,obs_list_tslots-1      
!   
!  Determine pair of times of field data bracketing this time slot
!  cdtime# is the date time in format 'yyyymmddhhmmss'
!  rtime# is the time relative to the central time of the assimilation period
!  #=1 is the beginning of time slot; #=2 the end
!  (1) is for inst or const files; (2) for tavg files that have a time offset
!  This formulation of rtime3 allows for use of nonconsecutive output times
     rtime1=(ntime-1)*obs_list_time_delta+obs_list_time_first
     rtime2=rtime1+obs_list_time_delta
     call time_compute_new_cdtime (cdtime_ref,cdtime1(1),rtime1,ierr)
     call time_compute_new_cdtime (cdtime_ref,cdtime2(1),rtime2,ierr)
     rtime3=field_tavg_offset
     call time_compute_new_cdtime (cdtime1(1),cdtime1(2),rtime3,ierr)
     rtime3=-field_tavg_offset 
     call time_compute_new_cdtime (cdtime2(1),cdtime2(2),rtime3,ierr)
!
     if (lprint) then
       print *,' '
       print ('(a,i2,2(2x,a))'),'LOOP over time: ',ntime, &
              cdtime1(1),cdtime2(1)
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
!  Read in all obs locations that occur within this time slot
     do ntype=1,ntypes
       obs_list_unit_in=un_list_in0+ntype
       call obs_list_clean
       call obs_list_types_get_setup (ntype)
       call obs_list_read_recs (ntime) 
       call obs_list_types_put_recs (ntype)
       if (myid == 0) then 
         call rec_list_r_put    ! save list of real values 
       endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
     enddo
!
!  Consideration of 2-d fields
!
     if (field_num_2d > 0) then
!
! Read all 2d fields at 2 successive times
       do nfid1=1,field_num_2d,npet
         nfid2=min(field_num_2d,nfid1+npet-1)
         nread=nfid2-nfid1
         do nfield1=0,nread
           ierr_read=0
           if (myid == nfield1) then
             nfield2=myid+nfid1
             itv=1+(field_types(nfield2,2)-1)/2 
             call read_shmem_data (trim(field_names(nfield2,2,2)),       &
                                   trim(field_files(nfield2,2)),         &
                                   trim(field_common_path),cdtime1(itv), &
                                   2,nfield2,field1,ierr_read) 
             call read_shmem_data (trim(field_names(nfield2,2,2)),       &
                                   trim(field_files(nfield2,2)),         &
                                   trim(field_common_path),cdtime2(itv), &
                                   2,nfield2,field2,ierr_read)
           endif
         enddo    ! loop over all fields in subset
!
         if (ierr_read /= 0) then
           print ('(2(a,i3))'),'Error in a call to read_shmem_data on proc=', &
                  myid,'  with error=',ierr_read
         endif
         call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
         if (ierr_read /= 0) then
           call mpi_die (my_name,ierr_read)
         endif
!
       enddo      ! loop over all successive subsets of fields
!
! Loop over observation list (locations) within this time slot
! This assumes that obs have already been sorted into time slots, but
! performs a check to make sure that the time_diff values are valid
! (which they may not be due to round-off errors).
!
       do ntype=1,ntypes
         call obs_list_clean
         call obs_list_types_get_setup (ntype)         
         call obs_list_types_get_recs (ntype)                 
         call rec_list_r_get 
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
         do obs_list_set=1,obs_list_tslots_recs(ntime),npet
           nread=min(npet,obs_list_tslots_recs(ntime)-obs_list_set+1)
           do nlist=1,nread
             if (myid == nlist-1) then
               obs_num=obs_list_set+nlist-1
               rec_id=obs_num+rec_ids(1,ntype,ntime)-1
               obs_time=obs_list_r(obs_list_id_time,obs_num)
               time_diff(1)=obs_time-rtime1 
               time_diff(2)=rtime2-obs_time
               if (time_diff(1) >= -.0001 .and. time_diff(2) >= -.0001) then
                 time_diff(1)=max(time_diff(1),0._rkind1) 
                 time_diff(2)=max(time_diff(2),0._rkind1) 
                 call interpolate_fields_2d (obs_num,field_num_2d,   &
                                            profiles_1fld)
                 profiles(1:field_num_2d,rec_id)= &
                    profiles_1fld(1:field_num_2d)
               else
                 print ('(a,i6,3a,i2,a,f10.6,a,4f10.6)'),'Obs ',obs_num, &
                      ' for dtype=',trim(dtypes(ntype)),' ntype=',ntype, &
                      ' at time=',obs_time,'outside range =',rtime1,rtime2
                 profiles(1:field_num_2d,rec_id)=bmiss                
               endif
             endif
           enddo
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
         enddo    ! loop over sets of npet headers
       enddo      ! loop over types
!
     endif      ! check on field_num_2d 
!
!
!  Consideration of 3-d fields
!
!  Loop over all required 3-D fields 
     nfid5=field_num_2d
     do nfield1=1,field_num_3d,2
!
       nfid1=nfid5+1
       nfid2=nfid1+field_kmax-1
       if (nfield1 < field_num_3d) then
         nfield2=nfield1+1
         nfid3=nfid2+1
         nfid4=nfid3+field_kmax-1
         nfid5=nfid4
       else   
         nfield2=0
         nfid3=nfid1
         nfid4=nfid2
         nfid5=nfid2
       endif
!
! Read either 1 or 2 fields at 2 successive times
       ierr_read=0
       if (myid == 0) then
         itv=1+(field_types(nfield1,3)-1)/2 
         call read_shmem_data (trim(field_names(nfield1,2,3)),       &
                               trim(field_files(nfield1,3)),         &
                               trim(field_common_path),cdtime1(itv), &
                               3,0,field1,ierr_read) 
       elseif (myid == 1) then
         itv=1+(field_types(nfield1,3)-1)/2 
         call read_shmem_data (trim(field_names(nfield1,2,3)),       &
                               trim(field_files(nfield1,3)),         &
                               trim(field_common_path),cdtime2(itv), &
                               3,0,field2,ierr_read)
       elseif (myid == 2 .and. nfield2 > 0) then
         itv=1+(field_types(nfield2,3)-1)/2  
         call read_shmem_data (trim(field_names(nfield2,2,3)),       &
                               trim(field_files(nfield2,3)),         &
                               trim(field_common_path),cdtime1(itv), &
                               3,0,field3,ierr_read)
       elseif (myid == 3 .and. nfield2 > 0) then
         itv=1+(field_types(nfield2,3)-1)/2   
         call read_shmem_data (trim(field_names(nfield2,2,3)),       &
                               trim(field_files(nfield2,3)),         &
                               trim(field_common_path),cdtime2(itv), &
                               3,0,field4,ierr_read)
       endif
!
       if (ierr_read /= 0) then
         print ('(2(a,i3))'),'Error in a call to read_shmem data on proc=', &
                myid,'  with error=',ierr_read
       endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
       if (ierr_read /= 0) then
         call mpi_die (my_name,ierr_read)
       endif
!
! Loop over observation list (locations) within this time slot
!
       do ntype=1,ntypes
         call obs_list_clean
         call obs_list_types_get_setup (ntype)         
         call obs_list_types_get_recs (ntype)                 
         call rec_list_r_get 
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         do obs_list_set=1,obs_list_tslots_recs(ntime),npet
           nread=min(npet,obs_list_tslots_recs(ntime)-obs_list_set+1)
           do nlist=1,nread
             if (myid == nlist-1) then
               obs_num=obs_list_set+nlist-1
               rec_id=obs_num+rec_ids(1,ntype,ntime)-1
               obs_time=obs_list_r(obs_list_id_time,obs_num)
               time_diff(1)=obs_time-rtime1 
               time_diff(2)=rtime2-obs_time 
               if (time_diff(1) >= -.0001 .and. time_diff(2) >= -.0001) then
                 time_diff(1)=max(time_diff(1),0._rkind1) 
                 time_diff(2)=max(time_diff(2),0._rkind1) 
                 call interpolate_fields_3d (obs_num,nfield1,field1,field2,   &
                                             profiles_1fld)
                 profiles(nfid1:nfid2,rec_id)=profiles_1fld(:)
                 if (nfield2 > 0) then          
                   call interpolate_fields_3d (obs_num,nfield2,field3,field4, &
                                               profiles_1fld)
                   profiles(nfid3:nfid4,rec_id)=profiles_1fld(:)
                 endif
               else
                 print ('(a,i6,3a,i2,a,f10.6,a,4f10.6)'),'Obs ',obs_num, &
                      ' for dtype=',trim(dtypes(ntype)),' ntype=',ntype, &
                      ' at time=',obs_time,'outside range =',rtime1,rtime2
                 profiles(nfid1:nfid4,rec_id)=bmiss                
               endif
             endif
           enddo
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
         enddo    ! loop over sets of npet headers
       enddo      ! loop over ntypes
!
     enddo ! loop over 3-d fields
!
!  Write proiles for all observation locations within this time slot
!  If T field is actually TV, then convert to T if QV is also present
     if (myid == 0 ) then 
       do ntype=1,ntypes
         call obs_list_clean
         call obs_list_types_get_setup (ntype)         
         call obs_list_types_get_recs (ntype)                 
         call rec_list_r_get 
         un_prof_out=un_prof_out0+ntype                  
!
         obs_num=obs_list_tslots_recs(ntime)
         if (obs_num > 0 .and. nfid5 > 0) then
           do nlist=1,obs_num
             rec_id=nlist+rec_ids(1,ntype,ntime)-1
             profiles_1loc(1:nfid5)=profiles(1:nfid5,rec_id)
             call write_profiles_tv2t (nfid5,nlist,profiles_1loc(1:nfid5))
             call write_profiles_o3units (nfid5,nlist,profiles_1loc(1:nfid5))
             call write_profiles_recs (nfid5,nlist,profiles_1loc(1:nfid5)) 
             if (lprint .and. print_sample .and. nlist == 1) then             
               call write_profiles_sample (nfid5,nlist,dtypes(ntype), &
                                           profiles_1loc(1:nfid5))
             endif
           enddo
           tot_num_obs(ntype)=tot_num_obs(ntype)+obs_num
           print *,' '
           print ('(3a,i2,a,i2,a,i8)'),'Profiles written for dtype=',         &
                                 trim(dtypes(ntype)),' ntype=',ntype,         &
                                 ' time slot=',ntime,' for number profiles=', &
                                 obs_num
         else
           print *,' '  
           print ('(3a,i2,a)'),'NO OBS FOR dtype=',trim(dtypes(ntype)), &
                               ' ntype=',ntype,' FOR THIS TIME SLOT:'
           print ('(a,2f12.4,a,i6)'),'rtime1,2 = ',rtime1,rtime2, &
                                     '  nfid5 = ',nfid5      
         endif
!
       enddo ! loop over types
     endif   ! check on myid=0
     
       call MPI_Barrier(MPI_COMM_WORLD,ierr)

!
   enddo   ! loop over nr file times
!
! Close profile output file
   if (myid == 0) then
     do ntype=1,ntypes
       un_prof_out=un_prof_out0+ntype                  
       call write_profiles_close
       print ('(3a,i2,a,i6)'),'Total number of profiles written for dtype=', & 
                              trim(dtypes(ntype)),' ntype=',ntype,           &
                              ' tot_num_obs=',tot_num_obs(ntype)
       print ('(a,i3)'),'Profile file closed for unit=',un_prof_out 
     enddo
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
! De-allocate arrays
   deallocate (profiles_1fld,profiles_1loc)
   call obs_list_clean
   call obs_list_types_clean
   call shutdown()
!
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine check(status, loc, ier)
!
! Check return status of NETCDF command function
!
     integer, intent(in) :: status
     integer, intent(inout) :: ier
     character(len=*), intent(in) :: loc
!
     if (status /= NF90_NOERR) then
       ier=ier+1
       if (loc /= ' ') then 
         print *,'Error at ', loc
         print *,NF90_STRERROR(status)
       endif
     endif
!
   end subroutine check
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine shutdown()
!
! shmem must deallocate shared memory arrays
!
     call MAPL_DeallocNodeArray(field1,rc=ierr)
     call MAPL_DeallocNodeArray(field2,rc=ierr)
     call MAPL_DeallocNodeArray(field3,rc=ierr)
     call MAPL_DeallocNodeArray(field4,rc=ierr)
     call MAPL_DeallocNodeArray(profiles,rc=ierr)
     call MAPL_DeallocNodeArray(rec_list_r,rc=ierr)
!
     call MAPL_FinalizeShmem (rc=ierr)
     call MPI_Finalize(ierr)
!
   end subroutine shutdown
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_shmem_data (f_name,f_file,f_common,cdtime, &
                              field_dim,k,field_in,iers) 
!
!  Read field on netcdf file into shared memory
!
     implicit none
!
     character(len=*), intent(in) :: cdtime
     character(len=*), intent(in) :: f_name                
     character(len=*), intent(in) :: f_file
     character(len=*), intent(in) :: f_common
     integer, intent(in) :: field_dim
     integer, intent(in) :: k
!
     integer, intent(out)  :: iers 
     real(rkind1), intent(out) :: field_in(field_imax,field_jmax,field_kmax)
!
     integer :: imx,jmx,kmx
     integer :: ier, ier1, ier2, ier3
     integer :: ncid,varid
     character(len=120) :: c_notice
     character(len=240) :: file_name
!
     call set_field_file_name (f_name,f_file,f_common,cdtime, &
                               file_name,ier) 
     iers=ier
!
     c_notice='Opening file for f='//trim(f_name)//' t='//cdtime
     call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid), &
                trim(c_notice),iers)
     if (iers /= 0) then
       print ('(4a)'),' ERROR attempting to open file of ',trim(f_name), &
                      ' field: ',trim(file_name)
     else  ! proceed since file successfully opened
!
! Get dimension information to check (lon)
       ier1=0
       ier2=0 
       call check (nf90_inq_dimid(ncid,'lon',varid),' ',ier1)
       if (ier1 /= 0) then 
         call check (nf90_inq_dimid(ncid,'longitude',varid),' ',ier2)
       endif
       if (ier1 == 0 .or. ier2 == 0) then
         call check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                  'nf90_inq 2',iers)
       else
         print *,'Error: cannot find one of lon or longitude on file'
         print *,'       file=',trim(file_name)
         print *,'       ier1, ier2, iers = ',ier1,ier2,iers 
         iers=iers+1
         imx=0
       endif
!
! Get dimension information to check (lat)
       ier1=0
       ier2=0 
       call check (nf90_inq_dimid(ncid,'lat',varid),' ',ier1)
       if (ier1 /= 0) then 
         call check (nf90_inq_dimid(ncid,'latitude',varid),' ',ier2)
       endif
       if (ier1 == 0 .or. ier2 == 0) then
         call check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                     'nf90_inq 4',iers)
       else
         print *,'Error: cannot find one of lat or latitude on file'
         print *,'       file=',trim(file_name)
         print *,'       ier1, ier2, iers = ',ier1,ier2,iers 
         iers=iers+1
         jmx=0
       endif
!
! Compare grid dimensions on file with those specified in program
       if (imx /= field_imax .or. jmx /= field_jmax) then 
         print *,'Grid dimension mismatch in routine : read_shmem_data'
         print *,'file_name=',trim(file_name)
         print *,'imax, jmax on file = ',imx,jmx
         print *,'imax, jmax in program = ',field_imax,field_jmax
         iers=iers+1
       endif
!
       if (field_dim == 3) then   ! check vertical dimension 
!
         ier1=0
         ier2=0 
         ier3=0       
         call check (nf90_inq_dimid(ncid,'lev',varid),' ',ier1)
         if (ier1 /= 0) then 
           call check (nf90_inq_dimid(ncid,'vertical level',varid),' ',ier2)
         endif
         if (ier2 /= 0) then 
           call check (nf90_inq_dimid(ncid,'index',varid),' ',ier3)
         endif
         if (ier1 == 0 .or. ier2 == 0 .or. ier3 == 0) then
           call check (nf90_inquire_dimension(ncid,varid,c_notice,kmx), &
                       'nf90_inq 6',iers)
         else
           print *,'Error: cannot find one of lev, vertical level, or index ', &
                   'on file'
           print *,'       file=',trim(file_name)
           print *,'       ier1, ier2, ier3, iers = ',ier1,ier2,ier3,iers 
           iers=iers+1
           kmx=0
         endif
         if (kmx /= field_kmax) then
           print *,'Grid dimension mismatch in routine : read_shmem_data'
           print *,'file_name=',trim(file_name)
           print *,'kmax on file=',kmx,'  kmax in program=',field_kmax
           iers=iers+1
         endif
       endif   ! test on whether 3D field
!
       if (iers == 0) then ! no errors in reading yet, so now read field
         c_notice='Getting vari for f='//trim(f_name)//' t='//cdtime
         call check (nf90_inq_varid(ncid,trim(f_name),varid),      &
                     trim(c_notice),iers)
         c_notice='reading field for f='//trim(f_name)//' t='//cdtime
         if (field_dim == 3) then ! read 3d field
           call check (nf90_get_var(ncid,varid,field_in),          &
                       trim(c_notice),iers)
         else                  ! read 2 d field but place in 3d array
           call check (nf90_get_var(ncid,varid,field_in(:,:,k)),   &
                       trim(c_notice),iers)
         endif
       endif
!
       c_notice='Closing file for f='//trim(f_name)//' t='//cdtime
       call check (nf90_close(ncid),trim(c_notice),iers)
!
     endif  ! test on whether file opened successfully
!
   end subroutine read_shmem_data 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interpolate_fields_2d (obs_num,nflds,profile1)
!
! Interpolate one or more 2d fields of data horizontally and temporally
! fld_t1 and fld_t2 are the same fields at successive times. Includes an 
! option for using area-mean rather than interpolated fields. 
!
     use m_nr_fields_info, only : field_avg 
     implicit none
!
     integer, intent(in) :: obs_num
     integer, intent(in) :: nflds
     real(rkind1), intent(out) :: profile1(field_kmax)
!
     integer, parameter :: n_area_dim=100
     integer :: ij_area(2,n_area_dim)
     integer :: h_index(2,2)
     integer :: n_area_max,n_area
     integer :: ipa,jpa,k
     integer :: nk
     integer :: time_max
     real(rkind1) :: f_int(2)              ! horz intp values 2 times
     real(rkind1) :: h_weights(2,2)        
     real(rkind1) :: hw(2,2)        
     real(rkind1) :: w_area(2,n_area_dim)
     real(rkind1) :: time_weights(2)        
     real(rkind1) :: tw(2)        
     real(rkind1) :: obs_lat               ! latitude of observation 
     real(rkind1) :: obs_lon               ! longitude of observation 
!
! Perform horizontal interpolation
!
     nk=nflds
     obs_lat=obs_list_r(obs_list_id_lat,obs_num)
     obs_lon=mod(obs_list_r(obs_list_id_lon,obs_num),360._rkind1)
     if (obs_lon < 0._rkind1) obs_lon=obs_lon+360._rkind1
!
! Determine grid box location that surrounds observation location
! and the interpolation weights for each such point
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                  obs_lat,obs_lon,h_index,h_weights)
!
! Determine weights for 2 times
     time_weights(1)=time_diff(2)/(time_diff(1)+time_diff(2))
     time_weights(2)=1._rkind1-time_weights(1)
!
! Loop over fields
     do k=1,nk
!
! Determine whether to use interpolated values (use of closest 
! point is included in this option) or to use area averaged values
       if (field_avg(k,2) < 1.0 .or. mod(field_types(k,2),2) == 0) then      
!
! Interpolate horizontally. Use weighted values at 2 times at 4 surrounding 
! points to the SW, NW, SE and NE. If indicated by field_type, use only 
! nearest point  horizontally by reassigning weight values to a single pont. 
         hw(:,:)=h_weights(:,:) 
         if (mod(field_types(k,2),2) == 0) then  
           call get_interp_horiz_nearest (hw)  
         endif
         f_int(1)=hw(1,1)*field1(h_index(1,1),h_index(1,2),k) &
                 +hw(1,2)*field1(h_index(1,1),h_index(2,2),k) &
                 +hw(2,1)*field1(h_index(2,1),h_index(1,2),k) &
                 +hw(2,2)*field1(h_index(2,1),h_index(2,2),k) 
         f_int(2)=hw(1,1)*field2(h_index(1,1),h_index(1,2),k) &
                 +hw(1,2)*field2(h_index(1,1),h_index(2,2),k) &
                 +hw(2,1)*field2(h_index(2,1),h_index(1,2),k) &
                 +hw(2,2)*field2(h_index(2,1),h_index(2,2),k)
!
       else
!
! Define values as a local sample area mean 
         call integrate_area_params (field_imax,field_jmax,n_area_dim, &
                                     obs_lat,obs_lon,field_lon_first,  &
                                     field_avg(k,2),n_area_max,ij_area,w_area)
         f_int(:)=0.
         do n_area=1,n_area_max
           ipa=ij_area(1,n_area)
           jpa=ij_area(2,n_area)
           f_int(1)=f_int(1)+w_area(1,n_area)*field1(ipa,jpa,k)
           f_int(2)=f_int(2)+w_area(1,n_area)*field2(ipa,jpa,k)
         enddo              
!
       endif ! check on whether to compute area average
!
! Interpolate in time
! If indicated by field_type, use only nearest point temporally
       tw(:)=time_weights(:) 
       if (field_types(k,2) > 1) then   
         call get_interp_time_nearest (tw)
       endif
       profile1(k)=tw(1)*f_int(1)+tw(2)*f_int(2) 
!
     enddo ! loop over 2 d field indexes
!
   end subroutine interpolate_fields_2d
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interpolate_fields_3d (obs_num,nfld,fld_t1,fld_t2,profile1)
!
! Interpolate a single 3d fields of data horizontally and temporally
! fld_t1 and fld_t2 are the same fields at successive times. Includes an 
! option for using area-mean rather than interpolated fields.
!
     use m_nr_fields_info, only : field_avg 
     implicit none
!
     integer, intent(in) :: obs_num
     integer, intent(in) :: nfld
     real(rkind1), intent(in)  :: fld_t1(field_imax,field_jmax,field_kmax) 
     real(rkind1), intent(in)  :: fld_t2(field_imax,field_jmax,field_kmax) 
     real(rkind1), intent(out) :: profile1(field_kmax)
!
     integer, parameter :: n_area_dim=100
     integer :: ij_area(2,n_area_dim)
     integer :: h_index(2,2)
     integer :: n_area_max,n_area
     integer :: ipa,jpa,k
     integer :: nk
     integer :: time_max
     real(rkind1) :: f_int_horiz(field_kmax,2) ! horz intp values 2 times
     real(rkind1) :: h_weights(2,2)        
     real(rkind1) :: w_area(2,n_area_dim)
     real(rkind1) :: time_weights(2)        
     real(rkind1) :: obs_lat               ! latitude of observation 
     real(rkind1) :: obs_lon               ! longitude of observation 
!
! Perform horizontal interpolation or integration
!
     nk=field_kmax
     obs_lat=obs_list_r(obs_list_id_lat,obs_num)
     obs_lon=mod(obs_list_r(obs_list_id_lon,obs_num),360._rkind1)
     if (obs_lon < 0._rkind1) obs_lon=obs_lon+360._rkind1
!
! Determine weights for 2 times
     time_weights(1)=time_diff(2)/(time_diff(1)+time_diff(2))
     time_weights(2)=1._rkind1-time_weights(1)
!
! If indicated by field_type, use only nearest temporally by reassigning 
! weight values to a single time
     if (field_types(nfld,3) > 1) then          
       call get_interp_time_nearest (time_weights)
     endif
!
! Determine whether to use interpolated values (use of closest 
! pont is included in this option) or to use area averaged values
     if (field_avg(nfld,3) < 1.0 .or. mod(field_types(nfld,3),2) == 0) then      
!
! Determine grid box location that surrounds observation location
! and the interpolation weights for each such point
       call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                    obs_lat,obs_lon,h_index,h_weights)
!
! If indicated by field_type, use only nearest temporally by reassigning 
! weight values to a single point
       if (mod(field_types(nfld,3),2) == 0) then  
         call get_interp_horiz_nearest (h_weights)
       endif
!
! Compute horizontally interpolated values of fields at 2 times
! This adds contributions by points to SW, NW, SE, NE         
       do k=1,nk
         f_int_horiz(k,1)=h_weights(1,1)*                     &
                          fld_t1(h_index(1,1),h_index(1,2),k) &
                         +h_weights(1,2)*                     &
                          fld_t1(h_index(1,1),h_index(2,2),k) &
                         +h_weights(2,1)*                     &
                          fld_t1(h_index(2,1),h_index(1,2),k) &
                         +h_weights(2,2)*                     & 
                          fld_t1(h_index(2,1),h_index(2,2),k) 
         f_int_horiz(k,2)=h_weights(1,1)*                     &
                          fld_t2(h_index(1,1),h_index(1,2),k) &
                         +h_weights(1,2)*                     &
                          fld_t2(h_index(1,1),h_index(2,2),k) &
                         +h_weights(2,1)*                     &
                          fld_t2(h_index(2,1),h_index(1,2),k) &
                         +h_weights(2,2)*                     &
                          fld_t2(h_index(2,1),h_index(2,2),k) 
       enddo
!
     else 
!
! define a profile by a sample area mean 
       call integrate_area_params (field_imax,field_jmax,n_area_dim, &
                                   obs_lat,obs_lon,field_lon_first,  &
                                   field_avg(nfld,3),n_area_max,ij_area,w_area)
       f_int_horiz(:,:)=0.
       do n_area=1,n_area_max
         ipa=ij_area(1,n_area)
         jpa=ij_area(2,n_area)
         do k=1,nk
           f_int_horiz(k,1)=f_int_horiz(k,1)+w_area(1,n_area)*fld_t1(ipa,jpa,k)
           f_int_horiz(k,2)=f_int_horiz(k,2)+w_area(1,n_area)*fld_t2(ipa,jpa,k)
         enddo              
       enddo
!
     endif ! check on whether to compute interpolated or area average values

!
! Perform interpolation in time
     do k=1,nk
       profile1(k)=time_weights(1)*f_int_horiz(k,1)  &
                  +time_weights(2)*f_int_horiz(k,2) 
     enddo
!
   end subroutine interpolate_fields_3d
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rec_list_r_put 
! 
!  Put all real-valued observation information in all time slots for a 
!  single observation type into a common array containg real-valued
!  information for all time slots and all observation types. 
!  
   implicit none
   integer :: nr1, nr2, nr3
!      
   do nr2=1,obs_list_tslots_recs(ntime)
     nr3=rec_ids(1,ntype,ntime)+nr2-1
     do nr1=1,obs_list_len_r
       rec_list_r(nr1,nr3)=obs_list_r(nr1,nr2)
     enddo
   enddo
!
   end subroutine rec_list_r_put 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rec_list_r_get 
!
!  The reverse of rec_list_r_put, extracting information for a single 
!  observation type.
!   
   implicit none
   integer :: nr1, nr2, nr3
!      
   do nr2=1,obs_list_tslots_recs(ntime)
     nr3=rec_ids(1,ntype,ntime)+nr2-1
     do nr1=1,obs_list_len_r
       obs_list_r(nr1,nr2)=rec_list_r(nr1,nr3)
     enddo
   enddo
!
   end subroutine rec_list_r_get 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
!      
 end program create_rad_profs
