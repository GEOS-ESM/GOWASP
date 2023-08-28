!
   program check_bufr
!
!  Main program for adding random observational error to simulated obs
!
   use m_kinds, only : rkind1
!
   use m_read_loop, only : read_loop_setup
   use m_read_loop, only : read_loop_do
   use m_read_loop, only : read_loop_clean
!
   use m_count_types, only : count_subtypes_print
!
   use m_saved_data, only : sd_setup
   use m_saved_data, only : sd_set_subtype
   use m_saved_data, only : sd_diff
   use m_saved_data, only : sd_clean
!
   use m_comp_stats, only : comp_stats_meanvar
   use m_comp_stats, only : comp_stats_chan_corr
   use m_comp_stats, only : comp_stats_hcorr_rad
   use m_comp_stats, only : comp_stats_hcorr_conv1
   use m_comp_stats, only : comp_stats_vcorr

!
   implicit none
!
   include "mpif.h"   ! defines some standard MPI variables
!
   logical :: check_class
   logical :: lprint
   logical :: ltest
   logical :: lrelhum
   logical :: lcorrfld ! true if random correlated fields required
!
   integer, parameter :: nerrors=11    ! number of types of counters
   integer, parameter :: ntypes=70     ! max number of subtypes counted
   integer, parameter :: change_year=0  ! 0 if year on file to be unchanged
!
   integer :: ROOT         
   integer :: i_file
   integer :: ierr          ! error flag
   integer :: iseed 
   integer :: myid          ! processor id number 0, ... ,npet
   integer :: npet          ! number of processors used
   integer :: n_mesg        ! counter of messages in file
   integer :: n_test        ! counter of obs selected for test printout
   integer :: CoresPerNode
   integer :: argc
   integer(4) :: iargc
   integer*8 :: i_random(1)
   integer*8 :: idatetime
   integer :: count_types(ntypes,3)  ! counter for obs sub-types
   integer :: ierrors(nerrors)
!
   integer :: obs_ids_dim1
   integer :: obs_ids_dim2
   integer :: obs_sing_dim2
   integer :: obs_multi_dim1
   integer :: obs_multi_dim2
   integer :: obs_multi_dim3
   real(rkind1), allocatable :: obs_ids(:,:,:)
   real(rkind1), allocatable :: obs_sing(:,:,:)
   real(rkind1), allocatable :: obs_multi(:,:,:,:)
!
   character(len=1)   :: ctest
   character(len=12)  :: dtype
   character(len=10)  :: cdatetime
   character(len=256) :: bufr_in_file(2)
   character(len=256) :: bufr_tab_file
   character(len=256) :: crtm_coef_dir
   character(len=12), parameter :: obs_classes(15) = (/ 'PREPBUFR',  &
      'AIRS', 'AMSUA', 'AMSUB', 'MSU', 'AMSUAAQUA', 'IASI', 'ATMS', &
      'HIRS2', 'HIRS3', 'HIRS4', 'MHS', 'GPSRO', 'SATWIND', 'GENRADTXT' /)
   character(len=*), parameter :: my_name='main_program'
!
!
   lprint=.true. 
!
! Read and check arguments
   argc = iargc()
   if (argc .ne. 7) then
     if (lprint) then
       print ('(3a)'),'usage must be: check_bufr.x dtype datetime ', &
          'file_bufr_1 file_bufr_2 bufr_tab_file crtm_coef_dir ctest'
     endif
   endif
   call GetArg( 1_4, dtype)      ! data type
   call GetArg( 2_4, cdatetime)  ! date time yyyymmddhh
   call GetArg( 3_4, bufr_in_file(1))  ! 1st input bufr file
   call GetArg( 4_4, bufr_in_file(2))  ! 2nd input bufr file
   call GetArg( 5_4, bufr_tab_file) ! 'none' or bufr table file, if required
   call GetArg( 6_4, crtm_coef_dir)
   call GetArg( 7_4, ctest)
!
   check_class=any(trim(dtype) == obs_classes)
   if (check_class) then
     if (lprint) then 
       print *,'Processing for d_type=',trim(dtype),'  cdatetime=', &
               trim(cdatetime)
     endif
   else
     if (lprint) then 
       print *,' Specification of d_type=',trim(dtype), &
               ' not one that is implimented'
     endif
     stop
   endif
!
   call sd_setup (dtype)
!
! Set some variables
   if (ctest == 'T') then
     ltest=.true.
   else
     ltest=.false.
   endif
!
! Open bufr files and initialize pert loop
   do i_file=1,2
     call read_loop_setup (bufr_in_file(i_file),'none',bufr_tab_file, &
                           crtm_coef_dir,lprint,dtype,i_file,ierr)
     call read_loop_do (lprint,change_year,ltest,dtype)
     call count_subtypes_print 
     call read_loop_clean (dtype)
   enddo
!
   call sd_diff (dtype) 
   call sd_set_subtype 
   call comp_stats_meanvar (dtype)
   if (trim(dtype) == 'IASI' .or. trim(dtype) == 'AIRS') then
     call comp_stats_chan_corr (dtype,2)
   endif
!
!  Do not compute horiz correl of radiances if there are too many 
!  obs locations, such as for MHS or AMSUA unless comp time increased
   if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'IASI' .or. &
       trim(dtype) == 'AMSUAAQUA') then
     call comp_stats_hcorr_rad (dtype,2)
   endif
!
   if (trim(dtype) == 'PREPBUFR') then
     call comp_stats_hcorr_conv1 (dtype,2)
   endif
!
   if (trim(dtype) == 'PREPBUFR' .or. trim(dtype) == 'GPSRO') then
     call comp_stats_vcorr (dtype,2)
   endif
!
   if (lprint) then
     print *,' '
     print *,'program completed' 
   endif
!
   call sd_clean
!
   end program check_bufr
   
