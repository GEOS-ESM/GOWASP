   program create_amv
!
! Create simulated amv obs using probabilistic model based on cloud 
! fraction and water vapor fields. Uses a file of probability model 
! parameters read. 
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
!
! Use of module that defines info for required fields
! This module also reads the .rc file
   use m_nr_fields_info, only : nr_fields_setup
   use m_nr_fields_info, only : field_imax, field_jmax, field_kdim 
   use m_nr_fields_info, only : field_lon_first, field_time_slots
   use m_nr_fields_info, only : field_time_delta, field_time_first
   use m_nr_fields_info, only : field_akbk_dlev, field_akbk
!
   use m_kinds, only : rkind1,rkind2
!
   use m_time_compute, only : time_compute_pack,time_compute_unpack
   use m_time_compute, only : time_compute_add
   use m_time_compute, only : rkindh
!
   use m_amv_fields
   use m_kx_table 
   use m_amv_view
   use m_satloc
!
! Use of the module responsible for some MPI-specific instructions
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
   include "mpif.h"
!
   logical :: lprint
   logical :: ltest
   logical :: ltest_1time
   logical :: lsim_counts
   logical :: visir_here
   logical :: wv_here
   logical :: l_interp_horiz
   logical :: l_interp_time
   logical :: l_interp_vert
   logical, parameter :: ltestgrid=.false.
!
   integer, parameter :: ncldmax=5          ! max number of cloud layers
   integer, parameter :: obs_data_max=40000 ! max obs per processor per time
   integer, parameter :: obs_data_dim1=10    ! number of obs variables to save
   integer :: nlats,nlons
   integer :: nband_lat,nband_lors,nband_time
   integer :: jbin
   integer :: i,j,k
   integer :: argc
   integer :: iseed
   integer :: i_random_seed(2)
   integer(4) :: iargc
   integer :: ncld                     ! number of cloud layers at location
   integer :: cld_klev(ncldmax,2)      ! cloud layer top and bot id
   integer :: id_ncld(10)              ! indicates if obs present for layer
   integer :: interp_nfields
!
   real(rkind1) :: lon, lat 
   real(rkind1) :: dlon, dlat          ! spacing between grid points (degrees)
   real(rkind1) :: xband
   real(rkind1) :: cld_frac(ncldmax)   ! max cloud frac in layer
   real(rkind1) :: cld_obsc(ncldmax)   ! prob of cloud is obscured from above
   real(rkind1) :: cld_plev(ncldmax,2) ! cloud layer top and bot p
   real(rkind1) :: cld_puvs(ncldmax,5) ! p,u,v,speed at obs amv level
!
! Variables concerning time keeping:
   integer :: ntime_periods      ! number of ana or BUFR dataset periods to process
   integer :: ntime_ana          ! count of analysis period
   integer :: ntime_nr           ! count of nr time spans in current analysis period
   integer :: ntime_kx           ! count of interpolation times in current ana period
   integer :: ntime_count        ! number of interp times considered
   integer :: ntime_satloc       ! count of index time slot in staloc file
   integer :: idatetime4         ! yymmddhh of initial center analysis period
   integer :: idatetime8         ! integer yyyymmddhh
!
   real(rkindh) :: dhours_ana    ! time between analysis (or BUFR file) periods 
   real(rkindh) :: dhours_nr     ! time between nr data files used
   real(rkindh) :: dhours_kx     ! spacing between (possibly) interpolated tumes
   real(rkindh) :: thours_ana    ! ana time period hours from first time considered
   real(rkindh) :: thours_kx     ! current t rel to beginning time
   real(rkindh) :: rhours_nr     ! nr data time relative to ana time
   real(rkindh) :: rhours_nr_next  ! time to check if nr time to be incremented 
   real(rkindh) :: rhours_kx     ! interpolated field time relative to ana time
   real(rkindh) :: time0_offset  ! first t in current ana period rel to ana time
   real(rkind1) :: tnorm         ! 
   real(rkind1) :: time0_ana(6)  ! (center) time of first analysis period
   real(rkind1) :: time1_ana(6)  ! (center) time of current analysis period 
   real(rkind1) :: time1_nr(6)   ! nr time at start of current interpolation period 
   real(rkind1) :: time2_nr(6)   ! nr time at end of current interpolation period 
   real(rkind1) :: time0_kx(6)   ! time for 1st interp field in this run
   real(rkind1) :: time1_kx(6)   ! time for 1st interp field in current ana period
   real(rkind1) :: time2_kx(6)   ! time for current interp field 
!
   character(len=14) :: cdtime0_ana  ! yyyymmddhh correspomding to time0_ana
   character(len=14) :: cdtime1_ana  ! yyyymmddhh correspomding to time1_ana
   character(len=14) :: cdtime1_nr   ! yyyymmddhh correspomding to time1_nr
   character(len=14) :: cdtime2_nr   ! yyyymmddhh correspomding to time2_nr
   character(len=14) :: cdtime0_kx   ! yyyymmddhh correspomding to time0_kx
   character(len=14) :: cdtime1_kx   ! yyyymmddhh correspomding to time1_kx
   character(len=14) :: cdtime2_kx   ! yyyymmddhh correspomding to time2_kx
   character(len=4)  :: cntimes_to_create ! argument = number of bufr file periods 
!
   character(len=*), parameter :: my_name='main_program'
   character(len=1)  :: ctest
   character(len=4)  :: ciseed
   character(len=240) :: field_list_file   ! rc file of NR field info
   character(len=240) :: kx_table_file_in  ! file of obs type info & prob params
   character(len=240) :: kx_table_file_out ! file of above + obs counts
   character(len=240) :: bufr_file_in      ! file with prepbufr table
   character(len=240) :: bufr_path_out     ! prepbufr file of satwind obs
   character(len=240) :: text_file_out     ! optional obs file for 1st 6hr 
!
!  Global arrays to be allocated using SHMEM
!  ---------------------------------------------
   real(rkind1), pointer :: obs_data(:,:,:) => null()
   real(rkind1), pointer :: obs_data_count(:) => null()
!
   integer :: dim1(1)
   integer :: dim3(3)
   integer :: ierr          ! returned error flag
   integer :: ierr_read     ! returned error flag from field reading routines
   integer :: myid          ! processor id number 0, ... ,npet
   integer :: myidp1        ! myid+1
   integer :: npet          ! number of processors used
   integer :: CoresPerNode
   integer :: die_proc
!
   real(rkind2), allocatable :: x_random(:,:)
   real(rkind1), allocatable :: ijtest(:,:,:)
   real(rkind1), allocatable :: prof_p(:,:)
   real(rkind1), allocatable :: prof_c(:)
   real(rkind1), allocatable :: prof_q(:)
   real(rkind1), allocatable :: prof_z(:)
   real(rkind1), allocatable :: prof_u(:), prof_v(:)
   real(rkind1), allocatable :: prof_ipwg(:,:)
   real(rkind1), allocatable :: ipw_puvs(:,:)
   real(rkind1), allocatable :: ipw_obsc(:)
!
!                                       ---
!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   if (myid == 0) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'
   die_proc_id=myid 
   myidp1=myid+1
!
!  ---------------------------------------------------------
!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   call MAPL_InitializeShmem(rc=ierr)
!
!  Initialize lprint so only info on processor 0 will be printed
   if (myid == 0) then
     lprint=.true.   
   else
     lprint=.false.
   endif   
!
! Read arguments defined when program is executed in script
   argc = iargc()
   if (lprint .and. argc /= 10) then
     print *,' usage must be: compute_params.x cdtime0, ntimes', & 
             ' field_list_file kx_table_file_in kx_table_file_out', &
             ' ciseed bufr_in bufr_path_out text_out ctest' 
     call mpi_die (my_name,77)
   endif
!
   call GetArg( 1_4, cdtime0_ana)  ! first BUFR dataset (ana) time to create
   call GetArg( 2_4, cntimes_to_create)   ! number of BUFR times to create
   call GetArg( 3_4, field_list_file)
   call GetArg( 4_4, kx_table_file_in)
   call GetArg( 5_4, kx_table_file_out)
   call GetArg( 6_4, ciseed)
   call GetArg( 7_4, bufr_file_in)
   call GetArg( 8_4, bufr_path_out)
   call GetArg( 9_4, text_file_out)
   call GetArg(10_4, ctest)
!
! Set some variables based on the arguments
   read (cntimes_to_create,'(i4)') ntime_periods
   read (ciseed,'(i4)') iseed
!
   if (ctest == 'T') then
     ltest=.true.
     ltest_1time=.true.
   else
     ltest=.false.
     ltest_1time=.false.
   endif
!
!  Get the field requirement info (field names and file templates)
   call nr_fields_setup ('SATWIND',field_list_file,lprint,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to nr_fields_setup: ierr=',ierr
     call mpi_die (my_name,ierr)
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
! 
! Read kx table and satloc distribution
! (...,6) means read all 6 parts of the table file
   call kx_table_read (kx_table_file_in,lprint,6,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to kx_table_read: ierr=',ierr
     call mpi_die (my_name,ierr)
   endif
!
   kx_obs_count(:,:,:,:)=0
   call satloc_setup (kx_num,kx_type,kx_locs)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)  
!
! Set variables concerning time based on kx_table and .rc file inputs
   dhours_nr=real(field_time_delta,rkindh)
   dhours_ana=dhours_nr*real(field_time_slots,rkindh)
   dhours_kx=dhours_ana/real(kx_field_slots,rkindh)
!
   call time_compute_unpack (cdtime0_ana,time0_ana)
   time0_offset=real(field_time_first,rkindh)
   call time_compute_add (time0_offset,time0_ana,time0_kx,ierr) 
   call time_compute_pack (time0_kx,cdtime0_kx)  
!
   if (lprint) then 
     print ('(a,3f9.4)'),'dhours_ana,dhours_nr,dhours_kx =', &
                          dhours_ana,dhours_nr,dhours_kx
     print ('(a,f8.3,2(2x,a))'),'time0_offset,cdtime0_ana,cdtime0_kx =', &
                              time0_offset,cdtime0_ana,cdtime0_kx
     print *,' '
   endif
!
   if (kx_field_imax == field_imax .and. kx_field_jmax == field_jmax) then
     l_interp_horiz=.false.
   else
     l_interp_horiz=.true.     
     if (myid == 0) then
       print ('(a,2i6,a,2i6)'),'Interpolate from nr nlons,nlats=', &
                  field_imax,field_jmax,' to',kx_field_imax,kx_field_jmax
     endif
   endif
!
   if (kx_field_slots == field_time_slots) then
     l_interp_time=.false.
     interp_nfields=7
   else
     l_interp_time=.true.
     interp_nfields=13
     if (myid == 0) then
       print ('(a,i3,a,i3)'),'Interpolate from nr time slots=', &
                         field_time_slots,' to',kx_field_slots
     endif
   endif
!
   if (kx_field_kdim == field_kdim) then
     l_interp_vert=.false.
   else
     l_interp_vert=.true.
     if (myid == 0) then
       print ('(a,i3,a,i3)'),'Interpolate from nr used nlevs =', &
                         field_kdim,' to',kx_field_kdim
       if (kx_akbk_dlev(1,1) < field_akbk_dlev(1,1) .or. &
           kx_akbk_dlev(kx_field_kdim,2) > field_akbk_dlev(field_kdim,2) ) then
         print *,' '
         print *,'Unacceptible kx_akbk:'
         print ('(2a,2f12.5)'),'field_akbk_dlev(1,1), ', &  
                      'field_akbk_dlev(field_kdim,2) =', &   
                      field_akbk_dlev(1,1),field_akbk_dlev(field_kdim,2)      
         print ('(2a,2f12.5)'),'kx_akbk_dlev(1,1), ',    &  
                      'kx_akbk_dlev(kx_field_kdim,2) =', &    
                      kx_akbk_dlev(1,1),kx_akbk_dlev(kx_field_kdim,2) 
         print *,'New data levels must lie within NR data levels used'         
         call mpi_die ('akbk def',92)
       endif
     endif
   endif
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr) 
!
! Allocate space for globally held fields 
   call amv_fields_allocate (interp_nfields,ierr)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  Allocate space for globally held obs_data
   dim3=(/obs_data_dim1,obs_data_max,npet/) 
   call MAPL_AllocNodeArray(obs_data,dim3,rc=ierr)
   dim1=(/npet/)
   call MAPL_AllocNodeArray(obs_data_count,dim1,rc=ierr)
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
! Allocate arrays of profiles at a single lat,lon
   allocate (prof_p(kx_field_kdim,3))
   allocate (prof_c(kx_field_kdim))
   allocate (prof_q(kx_field_kdim))
   allocate (prof_z(kx_field_kdim))
   allocate (prof_u(kx_field_kdim))
   allocate (prof_v(kx_field_kdim))
   allocate (prof_ipwg(kx_ipw_nlevs,2))
   allocate (ipw_puvs(kx_ipw_nlevs,5))
   allocate (ipw_obsc(kx_ipw_nlevs))
   prof_ipwg(:,:)=0.
   ipw_obsc(:)=0.  ! always 0 since there is no obscuring of ipw
!
! Initialize some variables and counters. This change allows printing
! of test grid results using a more managable size grid
   if (ltestgrid) then
     nlats=181
     nlons=360
   else                  
     nlats=kx_field_jmax
     nlons=kx_field_imax
   endif
   dlat=180./real(nlats-1)
   dlon=360./real(nlons)
   allocate (x_random(kx_num,nlats))
!
   if (ltestgrid) then
     call check_thining (0)
   endif
!
! Loop over analysis or BUFR dataset time periods
!
   ntime_count=0
   ntime_satloc=1         ! init time slot index used/modified by m_satloc
   do ntime_ana=1,ntime_periods
     thours_ana=real(ntime_ana-1,rkindh)*dhours_ana  ! ana hours since 1st
!
! Set time array and datetime to current analysis time     
     call time_compute_add (thours_ana,time0_ana,time1_ana,ierr)  
     call time_compute_pack (time1_ana,cdtime1_ana)
!
! Set reference time for 1st interpolated field in this ana period               
     call time_compute_add (time0_offset,time1_ana,time1_kx,ierr) 
     call time_compute_pack (time1_kx,cdtime1_kx)  
!
     rhours_nr_next=-1.e9    ! any large negative number
     ntime_nr=0              ! counter for nr times accessed in this ana period
! 
     if (lprint) then 
       print ('(a,i3,2x,a)'),'Processing ana period:',ntime_ana,cdtime1_ana
     endif
!
! Set seed for random number generator based on iseed and cdtime of ana period.
! (By using the ana time, results should not depend on how different job runs
! are used to create data for different sets of days.)
     read (cdtime1_ana(1:10),'(i10)') idatetime8
     read (cdtime1_ana(3:10),'(i8)')  idatetime4
     i_random_seed(1)=idatetime4
     i_random_seed(2)=iseed
     call random_seed (put=i_random_seed(1:2))
     if (lprint) then 
       print ('(a,2i10,a,i10,a)'),' Two Seeds for random number generator = ' &
              ,i_random_seed(1:2),'  for idatetime4=',idatetime4,' on proc=0' 
     endif
!
! Loop over interpolation times within each ana or BUFR dataset period
     do ntime_kx=1,kx_field_slots
!
       ntime_count=ntime_count+1
       rhours_kx=time0_offset+real(ntime_kx-1,rkindh)*dhours_kx
       thours_kx=rhours_kx+thours_ana
       call time_compute_add (rhours_kx,time1_ana,time2_kx,ierr)
       call time_compute_pack (time2_kx,cdtime2_kx)
!
       obs_data_count(:)=0.  ! initialize counter at each time
!
! Re set NR data times that sandwich the time for the interpolated field       
       if (rhours_kx > rhours_nr_next-0.001_rkindh) then ! go to next nr time 
         ntime_nr=ntime_nr+1 
         rhours_nr=time0_offset+real(ntime_nr-1,rkindh)*dhours_nr 
         rhours_nr_next=rhours_nr+dhours_nr
         call time_compute_add (rhours_nr,time1_ana,time1_nr,ierr)
         call time_compute_add (dhours_nr,time1_nr,time2_nr,ierr)
         call time_compute_pack (time1_nr,cdtime1_nr)  ! new begin slot time 
         call time_compute_pack (time2_nr,cdtime2_nr)  ! new end slot time
       endif
       nband_time=mod(ntime_kx,kx_jbins_time)
!
       if (ltest .and. lprint) then 
         print *,' '
         print ('(a,i3,2x,a,f8.2)'),'Processing for field time: ', &
               ntime_kx,cdtime2_kx,thours_kx
         print ('(a,i3,2(2x,a))'),'Using NR time slot nslot, t1, t2 ', & 
               ntime_nr,cdtime1_nr,cdtime2_nr
       endif
!
! Read NR data and compute ipw
       call amv_fields_read (myid,npet,nlats,nlons,cdtime1_nr,cdtime2_nr, &
                       interp_nfields,l_interp_time,l_interp_horiz,       &
                       ltest_1time,dlat,dlon,rhours_nr,dhours_nr,rhours_kx)
!
! Fill array of random numbers for use in determining set of thinned locs:
! first set all processors to the same seed; then preset the random numbers to
! use for setting the indexes for the first latitude and longitude considered
! such that these indexes will be the same on all processor.
       call random_fields_synch_ranf 
       do k=1,kx_num
         do j=1,nlats
           call random_number (x_random(k,j))
         enddo
       enddo
!
! Initialize stride.  This may depend on time slot. For some data types, 
! the i_stride is independent of lat, so it must also be determined here)  
       call amv_stride_j (dlat,x_random(:,1))
       call amv_stride_i (99.0_rkind1,dlon,x_random(:,nlats))
!
! Simulate longitudinal ranges of viewing by each satellite
       call satloc_bins_get (kx_num,kx_type,kx_satid,ntime_satloc, &
                             thours_kx,dhours_kx,cdtime0_kx,       &
                             kx_satloc_file,lprint,ierr_read)
       if (ierr_read /= 0) then
         call mpi_die (my_name,ierr_read)
       endif
!
! Loop over lat (exclude poles)
       do j=2,nlats-1
         if (mod(j-2,npet) == myid) then   ! distribute onto different procs
           lat=-90.+(j-1)*dlat
           xband=kx_jbins_lats*(lat+90.)/180.
           nband_lat=1+min(int(xband),kx_jbins_lats-1) 
! 
! Set longitudinal stride and first value
           call amv_stride_i (lat,dlon,x_random(:,j))
! 
! Loop over longitudes
           do i=1,nlons
             lon=field_lon_first+(i-1)*dlon
             if (lon < 0.) then
               lon=lon+360.
             endif
             if (lon >= 360.) then
               lon=lon-360.
             endif  
!
! set nband_lors dependent on fraction of sea coverage
             call amv_fields_sea_or_land (i,j,nband_lors)
!
! Determine scalar index jbins as though it was ordered as a 
! 3-d  array with dimensions jbins_lats,jbins_lors,jbins_time)
             jbin=nband_lat+kx_jbins_lats* &
                          (nband_lors+kx_jbins_lors*nband_time)
!
! Check if this is a possible location for an obs given the 
! obs thinning, satellite viewing (including geometry and temporal 
! characteristics) and solar elevation.
!
             call amv_check_loc (i,j,jbin,lat,lon,time2_kx, &  
                                 visir_here,wv_here)
!
             if (ltestgrid) then
               call check_thining (ntime_kx)
             endif
!
! If VIS/IR obs are present here for one or more obs types or sats
! extract required profiles from 3-d fields and determine cloud layers 
             if (visir_here .or. wv_here) then
               call amv_fields_extract_profiles (i,j,l_interp_vert,     &
                                      visir_here,wv_here,prof_p,prof_c, &
                                      prof_q,prof_z,prof_u,prof_v)
             endif
!      
             if (visir_here) then
               call amv_layers (kx_field_kdim,ncldmax,prof_p,prof_c,ncld, &
                                cld_klev,cld_frac,cld_plev)
               call amv_obscured (ncld,ncldmax,cld_frac,cld_obsc)
             endif
!
             do k=1,kx_num
               if (kx_present(k)) then 
!
                 if (kx_type(k)(2:2) == 'V' .or. kx_type(k)(2:2) == 'I' &
                     .or. kx_type(k)(2:2) == 'B') then
!
! Compute pressure, wind vector, and speed at assigned level of cloud.
! This may depend on k since different obs types use different algorithms
! for estimating the observed cloud level.
                   call amv_speed_cld (ncld,ncldmax,kx_field_kdim,k,cld_klev, &
                                       cld_plev,prof_u,prof_v,prof_z,cld_puvs)
                   call amv_use_params (kx_ifunc(k),k,jbin,ncldmax,ncld, &
                                        cld_frac,cld_obsc,cld_puvs,id_ncld)
                   call fill_obs_values (ncld,ncldmax,cld_puvs) 
!
                 else ! obs type is WV AMV
                   call amv_fields_ipw_grad (i,j,k,dlon,dlat,lon,lat,prof_ipwg)
                   call amv_speed_ipw (kx_field_kdim,k,prof_p,prof_q, &
                                       prof_u,prof_v,prof_z,ipw_puvs)
                   call amv_use_params (kx_ifunc(k),k,jbin,kx_ipw_nlevs, &
                                        kx_ipw_nlevs,prof_ipwg(:,1),     &
                                        ipw_obsc,ipw_puvs,id_ncld)  
                   call fill_obs_values (kx_ipw_nlevs,kx_ipw_nlevs,ipw_puvs) 
                 endif
!  
               endif   ! test on kx_present
             enddo     ! loop over k
!
           enddo   ! i
         endif     ! myid
       enddo       ! j
       call MPI_Barrier(MPI_COMM_WORLD,ierr)             
!
!  Write saved obs_info
       if (myid == 0 ) then
         call output_data
       endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)             
!
       ltest_1time=.false.
     enddo         ! sub time intervals
   enddo           ! BUFR dataset periods
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr) 
!
! For myid=0, get sum of values of kx_pt_count(:,:,1) residing on 
! all processors. Place the sum in array kx_pt_count(:,:,2)
   i=kx_jbins*kx_num  ! size of kx_pt_count(:,:,1)         
   call MPI_ALLreduce (kx_pt_count(:,:,1),kx_pt_count(:,:,2),i, &
                       MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
! Do the same for kx_obs_count
   if (lsim_counts) then  ! use table value read initially 
     kx_obs_count(:,:,:,2)=kx_obs_count(:,:,:,1)
   else                   ! use accumulated counts
     i=kx_pbins*kx_jbins*kx_num  ! size of kx_obs_count(:,:,1)         
     call MPI_ALLreduce (kx_obs_count(:,:,:,1),kx_obs_count(:,:,:,2),i, &
                         MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
   endif   
!
! Write out new counts
   if (myid == 0) then
     kx_obs_count(:,:,:,2)=kx_obs_count(:,:,:,2)/(ntime_periods)
     call kx_table_write (kx_table_file_out)
   endif
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)           
!
   if (ltestgrid) then
     call check_thining (-1)
   endif
!
! For comparing random x from different runs:
   if (ltest .and. lprint) then 
     call random_number (x_random(1,1))  
     print *,' '
     print *,'random number after last drawn = ',x_random(1,1)
   endif
!
   if (lprint) then
     print *,' '
     print *,'program completed' 
   endif
!
! De-allocate arrays
   call kx_table_clean  
   call MPI_Barrier(MPI_COMM_WORLD,ierr)           
   call shutdown()
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine check_thining (nt1)
!
   implicit none
   integer, intent(in) :: nt1
!   
   integer, parameter :: iunit=10
   integer :: mjx,mkx,ijsize
   character(len=*), parameter :: &
          filename='/discover/nobackup/rerrico/test_locs' 
!
   if (nt1 == 0) then
     allocate (ijtest(nlons,nlats,kx_num+1))
     ijtest(:,:,:)=0
   endif
!
   if (nt1 == 1) then
     do mkx=1,kx_num
       if (kx_present(mkx)) then
         ijtest(i,j,mkx)=ijtest(i,j,mkx)+1
       endif
     enddo
   endif
!
   if (nt1 == -1) then 
     ijsize=nlons*nlats
     do mkx=1,kx_num
       call MPI_Barrier(MPI_COMM_WORLD,ierr) 
       call MPI_ALLreduce (ijtest(:,:,mkx),ijtest(:,:,kx_num+1),ijsize, &
                       MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)      
       ijtest(:,:,mkx)=ijtest(:,:,kx_num+1)
     enddo
!
     if (myid == 0) then
       open (iunit,file=filename,form='unformatted')
       write (iunit) nlons,nlats,kx_num
!
       do mkx=1,kx_num
         write (iunit) mkx 
         do mjx=1,nlats
           write (iunit) mjx,ijtest(:,mjx,mkx)
         enddo
       enddo
!
       close (iunit)
       print *,' '
       print *,'Test file written: ',filename
     endif
     deallocate (ijtest)
!
   endif
!
   end subroutine check_thining 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine shutdown()
!
     deallocate (prof_p,prof_c,prof_q,prof_z,prof_u,prof_v)
     deallocate (x_random,prof_ipwg,ipw_puvs,ipw_obsc)
!
! shmem must deallocate shared memory arrays
!
     call amv_fields_shutdown()
     call MAPL_DeallocNodeArray(obs_data, rc=ierr)
     call MAPL_DeallocNodeArray(obs_data_count,rc=ierr)
!
     call MAPL_FinalizeShmem (rc=ierr)
     call MPI_Finalize(ierr)
!
   end subroutine shutdown
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine fill_obs_values (nused,nusedmax,puvs)
!
! Fill array of obs info for reports at this location for this obs sub type k
! There may be mutiple reports if there are obs at more than 1 p-level for 
! this k and location. These will be passed to the bufr write routine.
!
   implicit none 
!
   integer, intent(in) :: nused,nusedmax
   real(rkind1), intent(in) :: puvs(nusedmax,5)
!
   integer :: nuse
   integer :: ncnt_obs
   real(rkind1) :: hours99
!
! Multiplying by 0.99 here ensures that rhours_kx is within the interval
! between +/- 3 hours 
   hours99=real(rhours_kx,rkind1)*0.99 
!
   do nuse=1,nused
     if (id_ncld(nuse) /= 0) then
       obs_data_count(myidp1)=obs_data_count(myidp1)+1.
       ncnt_obs=min(nint(obs_data_count(myidp1)),obs_data_max)
       obs_data_count(myidp1)=real(ncnt_obs)
       obs_data(1,ncnt_obs,myidp1)=ncnt_obs
       obs_data(2,ncnt_obs,myidp1)=lon
       obs_data(3,ncnt_obs,myidp1)=lat
       obs_data(4,ncnt_obs,myidp1)=hours99       ! hrs from center window
       obs_data(5,ncnt_obs,myidp1)=kx_list(k)    ! kx
       obs_data(6,ncnt_obs,myidp1)=puvs(nuse,5)  ! z
       obs_data(7,ncnt_obs,myidp1)=puvs(nuse,1)  ! p
       obs_data(8,ncnt_obs,myidp1)=puvs(nuse,2)  ! u
       obs_data(9,ncnt_obs,myidp1)=puvs(nuse,3)  ! v
       obs_data(10,ncnt_obs,myidp1)=kx_satid(k)  ! WMO sat id# or 0
     endif
   enddo
!
   end subroutine fill_obs_values
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine output_data
!
!  Output amv obs and meta data (incl. call to write 1 BUFR report) 
!
   implicit none
!
   integer, parameter :: np_test=5 ! num of printed obs per proc when test 
   integer, parameter :: bufr_unit_in=20
   integer, parameter :: bufr_unit_out=21
   integer, parameter :: text_unit_out=22
   integer, parameter :: pts_in_msg=200 ! max # of obs reports per message
   integer :: n1,n2,n3
   integer :: np_skip
!
   character(len=244) :: text_file_out_augmented
   character(len=240) :: bufr_file_out
   character(len=*), parameter :: bufr_template= &
            '/Y$yyyy#/M$mm#/satwind.$yymmdd#.t$hh#z.prepbufr' 
!
! Open bufr file that has the bufr table to use
   if (ntime_ana == 1 .and. ntime_kx == 1) then 
     open(unit=bufr_unit_in,file=trim(bufr_file_in),form='unformatted')
     call openbf(bufr_unit_in,'IN ',bufr_unit_in)
     print *,' '
     print ('(3a,i3)'),' bufr input file=',trim(bufr_file_in),     &
            ' opened on unit=',bufr_unit_in
   endif
!
! Open output bufr file
   if (ntime_kx == 1) then 
     call set_field_file_name ('none',bufr_template,bufr_path_out, &
                               cdtime1_ana,bufr_file_out,ierr)
     open(unit=bufr_unit_out,file=trim(bufr_file_out),form='unformatted')
     print ('(3a,i3)'),' bufr output file=',trim(bufr_file_out),   &
            ' opened on unit=',bufr_unit_out
     call openbf(bufr_unit_out,'OUT',bufr_unit_in)  
     call maxout(200000)            ! increase size of ouput bufr
     call datelen (10) ! returns a 10 digit value of idate YYYYMMDDHH
   endif
!
! Open optional text obs file used to facilitate tuning
   if (ntime_kx == 1 .and. trim(text_file_out) /= 'none') then
     write (text_file_out_augmented,'(a,i3.3)'), &
          trim(text_file_out),ntime_ana
     open (unit=text_unit_out,file=trim(text_file_out_augmented))
     print *,' '
     print *,'text obs file opened to aid tuning on unit=',text_unit_out
     print *,'text_file_out name =', trim(text_file_out_augmented)
   endif
!
   n3=0  ! obs counter for this particular NR time 
   do n2=1,npet  ! loop over processors
     print ('(a,i2,1(2x,a,i8))'),'Obs count on proc=',n2,'obs=', &
            nint(obs_data_count(n2))
     np_skip=1+nint(obs_data_count(n2))/np_test
     do n1=1,nint(obs_data_count(n2))  ! loop over obs on processor
       n3=n3+1
       if (mod(n3,pts_in_msg) == 1) then   ! group reports into messages 
         call openmb (bufr_unit_out,'SATWND',idatetime8)
       endif          
!
! Write one obs bufr report 
       call bufr_satwind_write (bufr_unit_out,obs_data_dim1, &
                                obs_data(:,n1,n2))
!
! Print sample of output when testing
       if (ltest .and. mod(n1,np_skip) == 1) then
         call print_sample (n1,n2,n3)
       endif
!
! Write obs location to text file if requested (lon,lat,hours,kx,p(Pa),said)
       if (trim(text_file_out) /= 'none') then  
         write(unit=text_unit_out,fmt='(6f10.3)') obs_data(2:5,n1,n2), & 
              obs_data(7,n1,n2), obs_data(10,n1,n2)
       endif
!
       if (mod(n3,pts_in_msg) == 0) then  
         call closmg(bufr_unit_out)  
       endif 
!
     enddo
   enddo
!
! Close bufr after final report if not closed already
   if (ntime_kx == kx_field_slots) then
     if (mod(n3,pts_in_msg) /= 0) then   
       call closmg (bufr_unit_out)  
     endif 
     call closbf (bufr_unit_out)  
     print *,' '
     print *,'Output bufr file written'   
   endif
!
! Close optional text obs file if used
   if (ntime_kx == kx_field_slots .and. &
       trim(text_file_out) /= 'none') then
     obs_data(:,1,1)=999.
     write(unit=text_unit_out,fmt='(6f10.3)') obs_data(1:6,1,1)
     close (text_unit_out)
     print *,' '
     print *,' text_file_out closed'
   endif
!
   end subroutine output_data
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine print_sample (n1,n2,n3)
!
   implicit none
   integer, intent(in) :: n1,n2,n3   
!
   print ('(a,i8,i4,i8)'),'sample for n1,n2,n3=',n1,n2,n3
   print ('(10f14.3)'),obs_data(:,n1,n2)
!
   end subroutine print_sample 
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      

   subroutine random_fields_synch_ranf 
!
!  Synchronize random numbers so all processors create the same sequence
!  at this point
!
   implicit none
!
   integer :: ierrx
   integer :: nproc
   integer :: itag
   integer :: iseedx
   integer :: stat(MPI_STATUS_SIZE)
   integer :: i_random_seedx(2)
   real(8) :: x
   real(8) :: dum
!
   if (myid == 0) then
     call random_number (x)
     iseedx=int(1.e6*x)
     do nproc=1,npet-1
       itag=nproc
       call MPI_SEND (iseedx,1,MPI_INTEGER,nproc,itag,MPI_COMM_WORLD,ierrx)
     enddo
   endif
! 
   do nproc=1,npet-1
     if (myid == nproc) then
       itag=nproc
       call MPI_RECV (iseedx,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD,stat,ierrx)
     endif
   enddo
!
   i_random_seedx(1)=iseedx
   i_random_seedx(2)=1111    ! an arbitrary integer here
   call random_seed (put=i_random_seedx(1:2))
!
   end subroutine random_fields_synch_ranf 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end program create_amv



