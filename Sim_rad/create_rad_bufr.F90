! 
   program create_rad_bufr
!
!  Program driver for creating radiances from profiles (read from a 
!  file) of fields at observation locations and outputing them in 
!  a BUFR format readable by GSI. 
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
! 
      	!!--- use m_kinds           ! precision specification module 
   use m_kinds, only : rkind1

      	!!--- use m_bufr_rad        ! Read/Write BUFR data
   use m_bufr_rad, only: read_write_obs_rad
   use m_bufr_rad, only: check_rad_type
   use m_bufr_rad, only: read_write_gmi_1st_msg
!
      	!!--- use m_realloc_bad ! set which channels to use at each loc.
   use m_realloc_bad, only : realloc_bad_read
   use m_realloc_bad, only : realloc_bad_set
   use m_realloc_bad, only : realloc_bad_counter 
!
      	!!--- use m_bufr_rad        ! arrays for holding some obs info
   use m_rad_obs_arrays, only: rad_obs_arrays_setup
   use m_rad_obs_arrays, only: obs_generic_int

      	!!--- use m_read_profiles   ! observation and profile information
   use m_read_profiles, only: read_profiles_setup
   use m_read_profiles, only: read_profiles_recs
   use m_read_profiles, only: read_profiles_cleanup
   use m_read_profiles, only: prof_kmax,prof_akbk_int
   use m_read_profiles, only: prof_num_2d,prof_num_3d,prof_dim1,prof_all

      	!!--- use m_copy_rad_obs   ! copies obs into arrays used by m_bufr_rw
   use m_copy_rad_obs, only: copy_rad_obs

      	!!--- use m_time_compute   ! compute new date time for bufr messages
   use m_time_compute, only: time_compute_new_cdtime
      
      	!!--- use m_rad_prob        ! used to compute probability 
   use m_rad_prob, only: rad_prob_setup 
   use m_rad_prob, only: rad_prob_compute
   use m_rad_prob, only: list_aerosol_nums, list_cloud_nums
   use m_rad_prob, only: prob_seed, prob_names, prob_vars_flds
   use m_rad_prob, only: prob_vars, prob_vars_nums
    
      	!!--- use m_rad_index       ! setup indexs for all requied named fields
   use m_rad_index, only: rad_index_setup
   use m_rad_index, only: nf_ps, nf_siid, nf_said, nf_bmsg

      	!!--- use m_obs_list        ! observation header information
   use m_obs_list, only: obs_list_tslots, obs_list_subtypes
   use m_obs_list, only: obs_list_counter, obs_list_n_channels
   use m_obs_list, only: obs_list_i, obs_list_r, obs_list_c
   use m_obs_list, only: obs_list_len_r, obs_list_extra_recs
   use m_obs_list, only: obs_list_time_delta, obs_list_time_first
   use m_obs_list, only: obs_list_unit_in

      	!!--- use m_set_unit_nums ! set unit numbers
   use m_set_unit_nums, only: un_bufrin, un_bufrout  
   use m_set_unit_nums, only: un_prof_in1

       	!!--- use m_crtm_interface   ! interface module to call crtm
   use m_crtm_interface, only : crtm_interface_init 
   use m_crtm_interface, only : crtm_interface_comp_rad

       	!!--- use m_sat_info_table   ! read sat/inst info table
   use m_sat_info_table, only : sat_info_table_read

	!!--- use m_mympi   ! for typical message passing tasks.
   use m_mympi, only : COMM => mympi_comm_world
   use m_mympi, only : mympi_setup
   use m_mympi, only : mympi_bcast
   use m_mympi, only : mympi_gather
   use m_mympi, only : mympi_scatter
   use m_mympi, only : mympi_comm_split
   use m_mympi, only : mympi_comm_free
   use m_mympi, only : mympi_finalize

      	!!--- program will "die()" to kill all processors, instead of a single
	!!--- processor STOP.
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
!
   logical :: ltest_sample    ! true if some sample results to be printed
   logical :: lprint          ! controls if printing to be done
   logical :: subtype_present ! true if some obs for this subtype are present
   logical :: lcrtm_init      ! true if crtm has been initilized
   logical :: lcheck_dtype    ! true if dtype is in allowed list
   logical :: ldum            ! EOF flag argument not used here
   logical :: clear_sky_file  ! compute both clear-sky BT all-sky BT
!
   integer, parameter :: u_bufrin=un_bufrin    ! input unit for bufr data
   integer, parameter :: u_bufrin2=un_bufrin+1 ! input unit for separate table 
   integer, parameter :: u_bufrout=un_bufrout  ! unit for bufr output
   integer, parameter :: u_clear_sky=un_bufrout+1  ! unit for optional clear-sky BT
   integer, parameter :: ROOT=0
   integer, parameter :: n_sfctypes=4
   integer, parameter :: count_mod0=4 ! print sample of 4 for each subtype
   integer :: ntype,ntime,nfirst,nlast,nob,nobs
   integer :: n_channels    ! adjusted channel numbers
   integer :: effective_ks  ! index of level of effective radiating surface 
   integer :: subtype_id    ! id for obs subtype (either sat or instr id)
   integer :: ierr,ier
   integer :: count_all     ! used to number complete set of obs
   integer :: count_type    ! used for determining if obs to print as sample 
   integer :: count_mod1    ! used for determining if obs to print as sample 
   integer :: count_missing 
   integer :: count_missing_sum
   integer :: i,n,nid
   integer :: nsize         ! the size of an array
   integer :: subtype_count ! obs count for current subtype being read
   integer :: n_processors  ! number of processors used (for MPI)
   integer :: myid          ! individual processor id
   integer :: u_bufrtab     ! unit specifying bufr table to be copied
   integer :: argc
   integer :: sfc_type      ! flag: 1=ocean, 2=ice, 3=land, 4=snow
   integer :: sfc_types(n_sfctypes)  ! counter of surfaces for each type
   integer :: idatetime     ! integer equal to yyyymmddhh
   integer :: idatetime4    ! integer equal to yymmddhh
   integer :: newCOMM       ! sub-communicator where a profile presents.
   integer :: n_id
   integer :: prob_index
   integer :: prob_sum
   integer :: i_random_seed(2)
   integer(4) :: iargc
   integer, allocatable :: prob_counts(:) 
   integer, allocatable :: count_missing_M(:) 
   integer, allocatable :: sfc_type_M(:) 
   integer, allocatable :: prob_index_M(:) 
!
   real(rkind1), parameter :: bmiss=1.e10_rkind1  ! same as in create_rad_profs
   real(rkind1) :: effective_ps ! pressure of effective radiating surface
   real(rkind1) :: obs_quality  ! quality value required by IASI
   real(rkind1) :: dhours
   real(8)      :: x_random
   real(rkind1), allocatable :: prob_fracs(:)
   real(rkind1), allocatable :: prof_2d(:)   ! values of 2d fields at obs loc.
   real(rkind1), allocatable :: prof_3d(:,:) ! values of 3d fields at obs loc.
   real(rkind1), allocatable :: plevels(:,:) ! p at interfaces and data levels
   real(rkind1), allocatable :: crtm_obs(:)  ! brightness T or IASI radiance
   real(rkind1), allocatable :: crtm_clear(:)  ! brightness T for clear-sky 
   real(rkind1), allocatable :: obs_quality_M(:)
   real(rkind1), allocatable :: crtm_obs_M(:,:)
   real(rkind1), allocatable :: crtm_clear_M(:,:)
!
   character(len=1)   :: ctest           ! 
   character(len=1)   :: c_realloc       ! 'T' if real locat and channel used
   character(len=4)   :: obs_file_type   ! 'GTXT' or 'BUFR' file format 
   character(len=12)  :: obs_file_format ! 'formatted' or 'unformatted'
   character(len=16)  :: dtype           ! data type
   character(len=14)  :: cdtime_ref      ! synoptic date time : yyyymmddhhmmss
   character(len=14)  :: cdtime1         ! date time obs subset: yyyymmddhhmmss
   character(len=240) :: rcfile          ! file with probability params
   character(len=240) :: sat_info_file   ! file of some sat and instr params
   character(len=240) :: prof_file       ! file of profiles for obs
   character(len=240) :: bufr_in_file    ! file with BUFR table to copy
   character(len=240) :: bufr_tab_file   ! optional separate BUFR table file
   character(len=240) :: bufr_out_file   ! BUFR output with simulated obs
   character(len=240) :: crtm_coef_dir   ! directory containing crtm coef files
   character(len=240) :: clear_sky_suffix  ! add to sim_obs file name for clear BT file
   character(len=240) :: clear_sky_name    ! .bin file name of clear sky BT 
   character(len=*),parameter :: myname='main_program'
!
   call mympi_setup(COMM,n_processors,myid)
   die_proc_id=myid
   if (myid == ROOT) then
     lprint=.true.
   else
     lprint=.false.
   endif 
!
   if (lprint) then
     print *,' '
     print ('(a,i4,a)'),'BEGIN PROGRAM sim_obs_rad ON ',n_processors, &
                        ' PROCESSORS'
     print *,'COMM=',COMM
   endif
!
! Read and check arguments
   argc = iargc()
   if (argc /= 11 .and. argc /= 12) then
     if (lprint) then 
       print *,' usage must be: prog.x dtype cdtime rcfile sat_info_file',  &
               ' prof_file bufr_in_file bufr_tab_file bufr_out_file ctest', &
               ' c_realloc clear_sky_suffix'
     endif
     call MPI_Barrier (COMM,ier)                 
     call mpi_die (myname,71)
   endif
!
   call GetArg( 1_4, dtype)
   call GetArg( 2_4, cdtime_ref)
   call GetArg( 3_4, rcfile)
   call GetArg( 4_4, sat_info_file)
   call GetArg( 5_4, prof_file)
   call GetArg( 6_4, bufr_in_file)
   call GetArg( 7_4, bufr_tab_file)
   call GetArg( 8_4, bufr_out_file)
   call GetArg( 9_4, crtm_coef_dir)
   call GetArg( 10_4, ctest)
   call GetArg( 11_4, c_realloc)
   if (argc ==  12) then 
     call GetArg( 12_4, clear_sky_suffix)
   else
     clear_sky_suffix='none'
   endif
!
   if (lprint) then
     print *,' '
     call check_rad_type ('ANY',dtype,lcheck_dtype)     
     if (lcheck_dtype) then           ! This obs data type is allowed 
       print *,'Begin processing for type=',trim(dtype)  
     else
       print *,'Specification of d_type=',trim(dtype),' is not in allowed list'
       print *,'See routine check_rad_type in m_rad_bufr.f90 for allowed dtype'
       call MPI_Barrier (COMM,ier)                 
       call mpi_die (myname,72)
     endif
   endif 
!
   if (lprint .and. c_realloc == 'T') then 
     print *,'Only data thinned and QC OK locations and channels'
     print *,'actually used in an ealier experiment will be used here'
   endif
!
! Read and save table of sat/instr info
   call sat_info_table_read (sat_info_file,lprint,ierr)
   if (ierr /= 0) then 
     call MPI_Barrier (COMM,ier)                 
     call mpi_die ('main:sat_info_table_read',ierr)
   endif
!
! Specify whether obs file format is BUFR or generic text for radiances
    if (trim(dtype) == 'GENRADTXT') then
      obs_file_type='GTXT'
      obs_file_format='formatted'
    else
      obs_file_type='BUFR'
      obs_file_format='unformatted'
    endif
!
! Setup required profile info (read some info from profile file) 
   obs_list_unit_in=un_prof_in1
   if (lprint) print *,'Setup profile and obs arrays'
   call read_profiles_setup (prof_file,n_processors,ierr)
   call rad_obs_arrays_setup (obs_list_n_channels,3,2,100,5)
!
! Read rc file and set up information for computing probabilities of an 
! observation being affected by clouds or precip.
   if (lprint) print *,'Setup cloud/precip effect probabilities'
   call rad_prob_setup (dtype,rcfile,lprint,ierr)
   if (ierr /= 0) then
     call MPI_Barrier (COMM,ier)                 
     call mpi_die (myname,73)
   endif 
!
! Consruct .bin file name for clear-sky BT corresponding to sim all-sky BT  
! (only required if the obs sim file has all-sky values)
    if (trim(clear_sky_suffix) == 'none' .or. list_cloud_nums == 0) then
      clear_sky_file=.false.
      clear_sky_name='none'
    else
      clear_sky_file=.true.
      clear_sky_name=trim(bufr_out_file)//trim(clear_sky_suffix)
      open(u_clear_sky,file=trim(clear_sky_name),form='unformatted')
      if (lprint) then
        print ('(2a)'),'Clear-sky BT file opened for binary output: ', &
                        trim(clear_sky_name)
      endif
    endif
!
! setup variable indexes to be associated with required named values of fields
   if (lprint) print *,'Setup indexes to required obs header info'
   call rad_index_setup (dtype,lprint,ierr) 
   if (ierr /= 0) then
     call MPI_Barrier (COMM,ier)                 
     call mpi_die (myname,74)
   endif 
!
! ltest_sample is true if a sample of output is to be printed for testing
   if (ctest == 'T') then
     ltest_sample=.true.
   else
     ltest_sample=.false.
   endif
!
! Set seed for random number generator
   read (cdtime_ref(3:10),'(i8)') idatetime4    
   i_random_seed(1)=idatetime4
   i_random_seed(2)=prob_seed+myid
   call random_seed (put=i_random_seed(1:2))
   if (lprint) then
     print ('(a,2i10,a,i10)'),' Two Seeds for random number generator = ' &
            ,i_random_seed(1:2),'  for idatetime4=',idatetime4 
   endif
!
! Adjust channel numbers to be compatible with CRTM
   if (dtype (1:4) == 'HIRS' ) then
     n_channels=19
   else
     n_channels=obs_list_n_channels
   endif
!
! Allocates
   allocate (plevels(prof_kmax+1,2))
   allocate (prof_2d(prof_num_2d+1))
   allocate (prof_3d(prof_kmax,prof_num_3d+1))
   allocate (crtm_obs(n_channels))
   allocate (crtm_clear(n_channels))
   allocate (crtm_obs_M(n_channels,n_processors))
   allocate (crtm_clear_M(n_channels,n_processors))
   allocate (sfc_type_M(n_processors))   
   allocate (obs_quality_M(n_processors)) 
   allocate (prob_index_M(n_processors)) 
   allocate (count_missing_M(n_processors)) 
   allocate (prob_counts(prob_vars_flds+1))
   allocate (prob_fracs(prob_vars_flds+1))
   count_missing_sum=0
   sfc_types(:)=0  
   prob_counts(:)=0
   count_all=0
!
!  Open file that contains bufr table to be copied
!  (This may require a separate file with the table; e.g., for AIRS)
!  Also open Bufr file for output. If obs data file format is text, then
!  also read and write file header.
!
   if (myid == ROOT) then
!====================
! 
! Open pair of obs data files for reading and writing
     u_bufrtab=u_bufrin2 ! may change in the next call
     call open_obs_files (u_bufrin,u_bufrtab,u_bufrout,.true.,.false., &
                         lprint,dtype,bufr_in_file,bufr_tab_file,      &
                         bufr_out_file,cdtime_ref,obs_file_format,     &
                         obs_file_type,ierr)
!
   endif  ! I am ROOT
!
! For obs data files with format generic rad text, the file header must 
! be read on all processors because some values will be used to setup
! the CRTM interface on each processor. 
!
   if (obs_file_type == 'GTXT' .and. myid /= ROOT) then
     open(unit=u_bufrin,file=trim(bufr_in_file),form=trim(obs_file_format))
     n=-1
     call read_write_obs_rad (u_bufrin,.false.,dtype,n,.false.,  &
                              .true.,ldum,ierr)
     close (u_bufrin)
   endif
   call MPI_Barrier (COMM,ierr)                 
!
!
! Start loop over sets of profiles
!
   lcrtm_init=.false.   
   do ntype=1,obs_list_subtypes   ! loop over
!
! Read in first obs if this subtype if it exists, to get obs subtype
     count_type=0
     count_mod1=max(2,1+obs_list_counter(obs_list_tslots,ntype)/count_mod0)
     subtype_count=0
     if (sum(obs_list_counter(:,ntype)) > 0) then  
       subtype_present=.true.
       call read_profiles_recs (1,lprint,ierr)
     else
       subtype_present=.false.
     endif
!
!  Initalize CRTM  
     if (subtype_present) then
       if (nf_said /= 0) then
         subtype_id=obs_list_r(nf_said,1)    ! set subtype id as satellite id
       else
         subtype_id=obs_list_r(nf_siid,1)    ! set subtype id as instrument id
       endif 
!
! If requested, only use channels flagged as actually used (after QC) 
! from a previously-run assimilation Setup and read of a list of 
! such used channels occurs here.
       if (c_realloc == 'T') then 
         call realloc_bad_read (dtype,subtype_id,ctest,lprint)
       endif
!
! Initilaize the CRTM for this instrument and satellite
       call crtm_interface_init (dtype,lcrtm_init,subtype_id,lprint, &
                                 n_channels,crtm_coef_dir,ierr)
       if (ierr /= 0) then 
         call MPI_Barrier (COMM,ier)                 
         call mpi_die (myname,76)
       else
         lcrtm_init=.true.
       endif 
     endif
!
! Loop over obs times slots
     do ntime=1,obs_list_tslots-1
       dhours=obs_list_time_first+(ntime-1)*obs_list_time_delta
       call time_compute_new_cdtime (cdtime_ref,cdtime1,dhours,ierr)
       read (cdtime1(1:10),'(i10)') idatetime
!  
! Loop over groups of obs (number in each group upto n_processors)
       do nfirst=1,obs_list_counter(ntime,ntype),n_processors
         nlast=min(nfirst+n_processors-1,obs_list_counter(ntime,ntype))
         nobs=nlast-nfirst+1
!
!  Initialize to zero since some processors may have no obs to calculate
         count_missing=0
         obs_quality=0
         sfc_type=0
         prob_index=0
         crtm_obs(:)=0.
!    
! Loop over obs within current group being considered
         do nob=1,nobs
!
! Read obs (first of each subtype is read earlier)
           subtype_count=subtype_count+1
           if (subtype_count > 1) then
             call read_profiles_recs (nob,lprint,ierr) 
           endif
!
           if (nob == myid+1) then
             call prof_all_2d3d (prof_num_2d,prof_num_3d,prof_kmax, &
                                 prof_dim1,prof_all(:,nob),prof_2d,prof_3d)
! 
! The following condition checks if a read profile has been labeled as 
! outside its proper time span
             if (prof_all(1,nob) > 0.99*bmiss) then
               count_missing=1
             else
               call compute_plevs (prof_kmax,prof_akbk_int, &
                                   prof_2d(nf_ps),plevels)    
               call rad_prob_compute (nob,prof_kmax,plevels,prob_index, &
                                 dtype,effective_ks,effective_ps,obs_quality)
!
               call crtm_interface_comp_rad (obs_list_len_r,prof_num_2d, &
                      prof_num_3d,prof_kmax,n_channels,effective_ks,     &
                      effective_ps,dtype,subtype_id,lprint,              &
                      obs_list_r(:,nob),prof_2d(:),prof_3d(:,:),plevels, &
                      sfc_type,crtm_obs(:),crtm_clear(:),clear_sky_file,ierr) 
               if (ierr /=0) then
                 call mpi_die (myname,77)
               endif
             endif ! check on whether missing values are indicated        
           endif   ! consideration of processor id
!
         enddo     ! loop over obs within considered set
!
!       
         call MPI_Barrier (COMM,ierr)                 
!
         call mympi_gather (crtm_obs,crtm_obs_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(crtm_obs)',ierr)
         call mympi_gather (crtm_clear,crtm_clear_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(crtm_clear)',ierr)
         call mympi_gather (sfc_type,sfc_type_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(sfc_type)',ierr)
         call mympi_gather (prob_index,prob_index_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(prob_index)',ierr)
         call mympi_gather (obs_quality,obs_quality_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(obs_quality)',ierr)
         call mympi_gather (count_missing,count_missing_M,ROOT,COMM,stat=ierr)
         if (ierr /= 0) call mpi_die('main:mympi_gather(count_missing)',ierr)
         call MPI_Barrier (COMM,ierr)                 
!
         if (myid == ROOT) then
!
! Special instructions to write 1st message to GMI bufr data
! (obs_quality_M and crtm_obs_M are not used here)
           if (trim(dtype) == 'GMI' .and. ntype == 1 .and. & 
                       ntime == 1 .and. nfirst == 1) then
             call copy_rad_obs (n_channels,obs_list_n_channels,        &
                                obs_list_len_r,dtype,obs_quality_M(1), &
                                crtm_obs_M(:,1),obs_list_r(:,1),c_realloc)
             call read_write_gmi_1st_msg (u_bufrout,.false.,ierr,idatetime)
           endif
!
! Loop over obs considered by all processors 
           do nob=1,nobs
!
! Get index denoting sfc type (ocean,ice,land,snow)
             nid=sfc_type_M(nob)
             if (nid > 0 .and. nid <= n_sfctypes) then
               sfc_types(nid)=sfc_types(nid)+1
             endif
!
! Get index denoting whether there is a cloud or precip effect on this obs 
             nid=prob_index_M(nob)
             if (nid > 0) then
               prob_counts(nid)=prob_counts(nid)+1
             endif
!
! Write obs to bufr if it did not have missing values set in its profile
! (This occurs when a profile was labled as outside its proper time span)
             if (count_missing_M(nob) == 1) then
               count_missing_sum=count_missing_sum+1  ! profile missing flag 
             else
               count_type=count_type+1
               count_all=count_all+1
               call copy_rad_obs (n_channels,obs_list_n_channels,          &
                                  obs_list_len_r,dtype,obs_quality_M(nob), &
                                  crtm_obs_M(:,nob),obs_list_r(:,nob),c_realloc)
               if (obs_file_type == 'BUFR') then ! write BUFR message
                 call openmb (u_bufrout,trim(obs_list_c(nf_bmsg,nob)),idatetime)
               endif  
!
! If the assimilation to be run is to use the same channels at each location
! as a previously-run experiment, then check the data set of real locations 
! and channels. If a channel was not used, then distort this OSSE simulated 
! observation value so that it will later be elimated by QC when running 
! the OSSE assimilation. 
!
               if (c_realloc == 'T') then 
                 call realloc_bad_set (dtype) 
               endif
!
! Write obs to BUFR file
               call read_write_obs_rad (u_bufrout,.false.,dtype,           &
                                        count_all,.true.,.false.,ldum,ierr)
!
               if (clear_sky_file) then
                 call write_clear (count_all,u_clear_sky,crtm_clear_M(:,nob)) 
               endif
!
               if (ltest_sample .and. mod(count_type,count_mod1) == 1) then 
                 call print_sample (dtype,sfc_type_M(nob),                 &
                                    obs_quality_M(nob),obs_list_i(:,nob),  &
                                    obs_list_c(:,nob),crtm_clear_M(:,nob))
               endif
             endif    ! test on count_missing
           enddo      ! loop over obs
         endif        ! test on myid 
         call MPI_Barrier (COMM,ierr)                 
!
       enddo     ! loop over groups of n_processor obs
     enddo     ! loop over times
!
     if (lprint) then
       print *,' '
       print ('(a,i8)'),'Number of obs processed for this subtype=', &
                         subtype_count
     endif
!
! Print count of any locations not found in input .ods file if real locs used 
     if (c_realloc == 'T' .and. myid==ROOT) then
       call realloc_bad_counter (dtype)    
     endif
! 
   enddo       ! loop over subtypes
!
   if (myid==ROOT) then
     if (obs_file_type == 'BUFR') then     
       call closmg(u_bufrout)    ! final close out of messages
       call closbf(u_bufrin) 
       call closbf(u_bufrout)
     else
       close (u_bufrout)
     endif
     if (clear_sky_file) then
       close (u_clear_sky)
       print *,' '
       print ('(a,i3,2a)'),'Clear-sky file written: unit=',u_clear_sky, &
                           '  name=',trim(clear_sky_name)
     endif
   endif	
!
! For comparing random x from different runs:
   if (ltest_sample .and. lprint) then 
     call random_number (x_random)
     print *,' '  
     print *,'random number after last drawn = ',x_random
   endif
!
   if (lprint) then
!
! Print fractions of simulated observations determined for elevated surfaces
     prob_sum=sum(prob_counts(:))
     do n=1,prob_vars_flds+1
       prob_fracs(n)=real(prob_counts(n))/prob_sum
     enddo
     print *,' '  
     print ('(a,i8)'),'Number of observations produced=',prob_sum
     print *,'Fractions of simulated observation with surface set as:'
     do n=2,prob_vars_flds+1
       print ('(a,f8.4,a,f8.4,a,a)'),'frac=',prob_fracs(n),'  sigma=', &
            prob_vars(prob_vars_nums,n-1),'  determining field=',      &
            trim(prob_names(n-1))
     enddo 
     print ('(a,f8.4,a)'),'frac=',prob_fracs(1),'  sigma=  1.0000' 
!
! Print surface types specified for obs (all IR treated as ocean)
     print *,' '
     print *,' Numbers of profiles for each surface type:'
     print ('(4(i8,a))'), sfc_types(1),'=ocean',sfc_types(2),'=ice', &
                          sfc_types(3),'=land',sfc_types(4),'=snow'

     print *,' '
     print *,'Number of profiles processed (no missing values):', &
             sum(sfc_types(:))
     print *,'Number of profiles found with missing values =', &
             count_missing_sum
   endif  ! test on lprint
!
! Deallocate all arrays        
   call read_profiles_cleanup (.true.)  
   deallocate (prof_2d,prof_3d)
   deallocate (crtm_obs,crtm_obs_M)
   deallocate (sfc_type_M)
   deallocate (prob_counts,prob_index_M)
   deallocate (obs_quality_M,count_missing_M)
   call mympi_finalize()
!
   end program create_rad_bufr
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
    subroutine print_sample (dtype,sfc,quality,obs_i,obs_c,obs_clear)
!
!  Print sample of output 
!
    use m_kinds, only : rkind1
    use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
    use m_rad_obs_arrays, only : obs_n_channels, obs_n_data
    use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c  
    use m_obs_list, only : obs_list_names
!
    implicit none      
    integer, intent(in) :: sfc
    integer, intent(in) :: obs_i(obs_list_len_i)
    real(rkind1), intent(in) :: quality
    real(rkind1), intent(in) :: obs_clear(obs_n_channels)
    character(len=*), intent(in) :: dtype
    character(len=*), intent(in) :: obs_c(obs_list_len_c)
!
    integer :: n,n1,n2,n3
    real(rkind1) :: output_info(10)
!
    print *,' '
    print ('(3a)'),' Sample output for obs type = ',trim(dtype), &
                   ': header and derived info follow' 
    print ('(a,i8)'),(obs_list_names(n),obs_i(n),n=1,obs_list_len_i)
    n3=obs_list_len_i+obs_list_len_r
    print ('(a,2x,a)'),(obs_list_names(n+n3),obs_c(n),n=1,obs_list_len_c)
    n3=obs_list_len_i
    do n1=1,obs_list_len_r,4
      n2=min(n1+3,obs_list_len_r)
      print ('(4(2x,a16,f13.5))'),(obs_list_names(n+n3),obs_info(n),n=n1,n2)
    enddo
    print ('(a,i2,a,f7.2)'),' sfc_type_index=',sfc,'  quality=',quality
    print ('(a,i3,a)'),'  Information for ',obs_n_channels,' channels:'
    do n=1,obs_n_channels
      output_info(1:2)=obs_channels(n,1:2)
      output_info(3:2+obs_n_data)=obs_values(n,1:obs_n_data)
      if (obs_clear(n) > 0.) then
        print ('(i4,1p6e16.7)'),n,output_info(1:2+obs_n_data),obs_clear(n)
      else
        print ('(i4,1p5e16.7)'),n,output_info(1:2+obs_n_data)
      endif
    enddo
!
    end subroutine print_sample 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
    subroutine write_clear (ncount,iunit,crtm_clear)
!
!  Write clear-sky BT and identifying info for one observation to a .bin file
!
    use m_kinds, only : rkind1
    use m_rad_obs_arrays, only : obs_info, obs_n_channels
    use m_rad_index, only : nf_lat, nf_lon, nf_time
!
    implicit none
!
    integer, intent(in) :: ncount, iunit
    real(rkind1), intent(in) :: crtm_clear(obs_n_channels)
!
    write (iunit) ncount,obs_info(nf_lat),obs_info(nf_lon),obs_info(nf_time), &
                  crtm_clear(1:obs_n_channels)
    end subroutine write_clear         
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
