    Program create_rad_obs_list
!
! Program to read a bufr file of radiance data, thin the data set, and 
! write out a file of just observation locations and other header information
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
    use m_kinds, only : rkind1
!
    use m_realloc_get, only : realloc_get_read
    use m_realloc_get, only : realloc_get_find
!
    use m_rad_thin, only : rad_thin_clean
    use m_rad_thin, only : rad_thin_setup
    use m_rad_thin, only : rad_thin_put
    use m_rad_thin, only : rad_thin_write
    use m_rad_thin, only : rad_thin_box_count
!
    use m_bufr_rad, only : read_write_obs_rad
    use m_bufr_rad, only : check_rad_type
    use m_bufr_rad, only : read_write_gmi_1st_msg
!
    use m_rad_obs_arrays, only : rad_obs_arrays_setup 
    use m_rad_obs_arrays, only : rad_obs_arrays_clean 
    use m_rad_obs_arrays, only : obs_generic_char, obs_generic_int
!
    use m_set_unit_nums, only : un_bufrin 
!
    implicit none
!
    logical, parameter :: ltest_all=.true.
    logical, parameter :: lstop=.true.
    logical :: test_sample
    logical :: lprint
    logical :: ltest
    logical :: lerror
    logical :: lcheck_dtype    ! true if dtype is in allowed list
    logical :: msg_check
    logical :: leof            ! EOF flag only used if obs_file_type='GTXT' 
    logical :: generic_ireadsb ! generalization of BUFR ireadsb
    logical :: generic_ireadmg ! generalization of BUFR ireadmg 
!
    integer, parameter :: lunit_bufr=un_bufrin
    integer, parameter :: n_errors=11 ! number of error counters
    integer, parameter :: header_info_nums=2  ! # of header info text read
    integer :: n_test_sample          ! every n-th obs will be printed as sample
    integer :: ierrors(n_errors)      ! error counter (last is number of OK obs)
    integer :: lunit_bufr_tab         ! unit number for bufr table
    integer :: idate_msg              ! date/time YYYYMMDDHH read as message
    integer :: n_mesg(1)              ! number of bufr messages read
    integer :: nobs                   ! number of obs reports read    
    integer :: ier                   
    integer :: ireadmg
    integer :: ireadsb
    integer :: argc                    ! number of arguments
    integer(4) :: iargc
!
! Arguments 
    character(len=16)  :: dtype           ! data type
    character(len=14)  :: cdatetime_old   ! yyyymmddhh of read in bufr data
    character(len=14)  :: cdatetime_new   ! yyyymmddhh of write out in bufr data
    character(len=240) :: rcfile          ! file of thin params & needed fields
    character(len=240) :: inputfile       ! BUFR input with obs locations
    character(len=240) :: outputfile      ! BUFR output with sim. obs
    character(len=240) :: bufr_tab_file   ! optional seeparate BUFR table file
    character(len=1)   :: c_test          ! 'T' if extra info to be printed
    character(len=1)   :: c_realloc       ! 'T' if real locat and channels to use
    character(len=8)   :: msg_set_prev ! previous message set read in msg hdr
    character(len=8)   :: msg_set_curr ! current message set read in msg hdr
    character(len=8)   :: msg_set_1    ! first message set indicated in file
    character(len=4)   :: obs_file_type   ! 'GTXT' or 'BUFR'
    character(len=12)  :: obs_file_format ! 'formatted' or 'unformatted'
    character(len=240) :: header_info(header_info_nums)
    character(len=*),parameter :: myname='create_rad_obs_list'
!
    print *,' '
    print ('(2a)'),'BEGIN PROGRAM ',myname
!
! Read and check arguments
    argc = iargc()
    if (argc /= 11) then
      print *,' usage must be: sim_rad.x dtype old_datetime',      &
              ' new_datetime rcfile inputbufr outputbufr',         &
              ' bufr_table_file header_info1 header_info2 c_test', &
              ' c_realloc'
      if (lstop) then
        stop
      endif
    endif
    call GetArg( 1_4, dtype)
    call GetArg( 2_4, cdatetime_old)
    call GetArg( 3_4, cdatetime_new)
    call GetArg( 4_4, rcfile)
    call GetArg( 5_4, inputfile)
    call GetArg( 6_4, outputfile)
    call GetArg( 7_4, bufr_tab_file)
    call GetArg( 8_4, header_info(1))
    call GetArg( 9_4, header_info(2))
    call GetArg(10_4, c_test)
    call GetArg(11_4, c_realloc)
    print *,' '
    call check_rad_type ('ANY',dtype,lcheck_dtype)     
    if (lcheck_dtype) then           ! This obs data type is allowed 
      print *,'Begin processing for type=',trim(dtype)  
    else
      print *,'Specification of d_type=',trim(dtype),' is not in allowed list'
      print *,'See routine check_rad_type in m_rad_bufr.f90 for allowed dtype'
      if (lstop) then 
        stop
      endif 
   endif
!
   if (c_realloc == 'T') then 
     print *,'Only data thinned and QC OK locations and channels'
     print *,'actually used in an ealier experiment will be used here'
   endif
!
! test_sample is true if a sample of output is to be printed for testing
    if (c_test == 'T') then
      test_sample=.true.
      lprint=.true.
    else
      test_sample=.false.
      lprint=.true.
    endif
!
! if test_sample=.true., then every n_test_sample obs read will be printed
    if (dtype(1:4) == 'AMSU' .or. dtype(1:4) == 'HIRS' ) then 
      n_test_sample=25000
    elseif (dtype(1:4) == 'ATMS') then 
      n_test_sample=50000
    elseif (dtype(1:6) == 'AVCSAM' .or. dtype(1:6) == 'AVCSPM') then 
      n_test_sample=1000000
    else 
      n_test_sample=150000
    endif
!
! Setup max sizes for arrays to hold read observation information 
    if (dtype(1:9) == 'GENRADTXT') then
      call rad_obs_arrays_setup (3500,3,2,20,5)
    else
      call rad_obs_arrays_setup (700,3,2,100,5)
    endif
!
! Open obs data input files.  Also extract information from obs data 
! read/write routines.  For dtype='GENRADTXT' also read file header.
! Arguments 0 and 'none' indicate that no obs data file is written.
    lunit_bufr_tab=lunit_bufr+1 
    call open_obs_files (lunit_bufr,lunit_bufr_tab,0,.true.,.true., &
                         lprint,dtype,inputfile,bufr_tab_file,      &
                         'none',cdatetime_new,obs_file_format,      &
                         obs_file_type,ier)
!
! If requested, use locations of observations used by a previously-run
! data assimilation (after its data thinning and QC) as given by an 
! assimilation diagnostic output file. Here, read a text extract from 
! that file for reference later. 
!
    if (c_realloc == 'T') then 
      call realloc_get_read (dtype,c_test) 
    endif
!
! Declare and initialize arrays, including fields, required for data thinning

    call rad_thin_setup (dtype,rcfile,cdatetime_new,lstop,lprint,ier)
    if (ier /= 0) stop
!
! Initialize some variables and counters
    n_mesg(:)=0
    ierrors(:)=0
    ltest=ltest_all
    nobs=0  
    msg_set_prev(1:8)=' ' 
!
! Begin loop to read observations from existing file
! 
     do while (generic_ireadmg(lunit_bufr,leof,obs_file_type, &
                              msg_set_curr,idate_msg))
!
! If a generic rad obs text file, reset idate and subset names. 
      if (obs_file_type /= 'BUFR') then
        idate_msg=obs_generic_int(2) 
        msg_set_curr=obs_generic_char(5)(1:8)
      endif
!
! If new data type, set previous type to new one
      if (trim(msg_set_curr) /= trim(msg_set_prev )) then
        print ('(3a,i10)'),' Processing message subset= ',trim(msg_set_curr), &
                           ' for msg date =',idate_msg
        msg_set_prev=msg_set_curr
        if (nobs==0) then ! Save first value of message header with info
          msg_set_1=msg_set_curr
        endif
      endif
!
      call check_rad_type (msg_set_curr,dtype,msg_check)
      if (msg_check) then         
        do while (generic_ireadsb(lunit_bufr,leof,obs_file_type))
!
! Read next observation report
          call read_write_obs_rad (lunit_bufr,.false.,dtype,nobs,.true., &
                                   .true.,leof,ier)
!
          if (.not. leof) then  ! not at end of file
!
! If requested, only use obs locations for this satellite and instrument 
! that are in the list of previously recorded as used by an assimilation 
! after QC and thining. If the BUFR data is not found in that list, the 
! BUFR meta-data location is given an unacceptable value so that it will 
! be discarded by what follows. The arguement .true. effectively removes 
! duplicates (or near duplicates) from the allowed BUFR data. 
            if (c_realloc == 'T') then 
              call realloc_get_find (.true.)
            endif
!
! Save required information after applying data thinning 
            call rad_thin_put (cdatetime_old,cdatetime_new,n_errors,dtype, &
                               nobs,ierrors,lerror)
            if (.not. lerror) then
              if (test_sample .and. &
                        mod(ierrors(n_errors),n_test_sample) == 1) then
                call print_sample (dtype)
              endif
            endif   ! check on lerror
          endif     ! check on leof
!
! End loops to read in observation buffers
        enddo               ! loop over do while ( ireadsb(luin) .eq. 0 )
      endif                 ! check on modtype
      n_mesg(1)=n_mesg(1)+1 ! count messages read in original file
    enddo                   ! loop over (ireadmg(luin,subset,idate).eq. 0)
!
    if (obs_file_type == 'BUFR') then
      call closbf(lunit_bufr) 
      if (lunit_bufr /= lunit_bufr_tab) then 
        call closbf(lunit_bufr_tab) 
      endif
    else   ! close generic file
      close (lunit_bufr)
    endif
!
    print ('(3a)'),'SUMMARY TABLE FOR ',trim(dtype),':' 
    print ('(i8,a)'),nobs,' observation reports read' 
    print ('(4a)'),'bufr datetime=',trim(cdatetime_old), &
                   ';  nr datetime=',trim(cdatetime_new)
    print ('(i8,2a)'),n_mesg(1),' observation message-groups read all types'
    print ('(a,i8)'),'number of OK obs=',ierrors(n_errors)
!
! Print out error count summary
    if (sum(ierrors(1:n_errors-1)) > 0) then
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
           ' other problems or thinning using instrument header found' 
      print ('(i8,a)'),ierrors(7), &
           ' observation reports found but sat or instrument ID not used' 
    endif
!
! Compute and print numbers of obs for each time slot and type
! Dump out all thinned observations (ierrors(n_errors) is number of OK obs read)
    if (ierrors(n_errors) > 0) then 
      call rad_thin_box_count (lprint)
      call rad_thin_write (outputfile,header_info_nums,lprint, &
                           header_info,msg_set_1,ier)
    else
      print *,' '
      print *,'NO OBS TO PROCESS FOR ',trim(dtype)
    endif
!
    call rad_thin_clean
    call rad_obs_arrays_clean 
    print *,' '
    print *,'Program complete'
!
!
    end program create_rad_obs_list
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
    subroutine print_sample (dtype)
!
!  Print sample of output 
!
    use m_kinds, only : rkind1
    use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
    use m_rad_obs_arrays, only : obs_n_channels, obs_n_data
    use m_rad_thin, only : thin_info_dim,thin_info_names,thin_info
!
    implicit none      
    character(len=*), intent(in) :: dtype
!
    integer :: n,n1,n2
    real(rkind1) :: output_info(10)
!
    print *,' '
    print ('(3a)'),' Sample output for obs type = ',trim(dtype), &
                   ': header and derived info follow' 
    do n1=1,thin_info_dim,4
      n2=min(n1+3,thin_info_dim)
      print ('(4(2x,a16,f12.4))'), ((thin_info_names(n),thin_info(n)),n=n1,n2)
    enddo
    print ('(a,i3,a)'),'  Information for ',obs_n_channels,' channels:'
    do n=1,obs_n_channels
      output_info(1:2)=obs_channels(n,1:2)
      output_info(3:2+obs_n_data)=obs_values(n,1:obs_n_data)
      print ('(i4,1p5e16.7)'),n,output_info(1:2+obs_n_data)
    enddo
!
    end subroutine print_sample 



