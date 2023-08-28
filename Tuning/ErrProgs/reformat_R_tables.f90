    program refort_R_tables
!
! Copy part of GSI sat error table (slightly reformatted) to OSSE format
! The GSI table needs to be modified first by:
!   a) copy the first 3 lines from an OSSE sat error file header 
!   b) the value in ngroups should be the number of new groups to copy; 
!     a "group" is a set of channels for a distinct data type and satellite
!   c) before each group begins add a counter in lines 1-2 (right justified). 
!   d) after each counter line add
!      1) left justified in columns 2-17, the value of dtype for the OSSE code
!      2) left justified in columns 21-37, the name of the particular satellite
!      3) right justified in columns 39-42, the number of channels to follow 
!
! This is intended to be applied to only new groups to be used by the OSSE in order to
! have proper names, channel umbers, and starting values for error std. devs. to consider.
! Note that values of dtype and sat names must match a combination appearing in the first
! 2 columns of the GOWASP sat_info.rc file to be used.
!
! Set the 4 parameters below 
!
    implicit none
!
    integer, parameter :: n_chan_max=600 ! max umber of channels for any group read
    integer, parameter :: n_groups_max=8 ! >= number ofgroups to be read
    real, parameter :: err_factor=0.5    ! mult evalues by this before writing
    logical, parameter :: allsky=.false. ! .true. if all sky values to be copied
!  
    integer :: argc
    integer :: i, j, n
    integer :: i_type
    integer :: ihead   ! number of header records to skip in sat_error file
    integer :: igroups ! number of sat platforms or instruments in table
    integer :: n_channels(n_groups_max)
    integer :: n_chan(n_chan_max,n_groups_max)
    integer :: n_use(n_chan_max,n_groups_max)
    integer :: inum(n_groups_max)
!
    real(4) :: rvalues(5,n_chan_max,n_groups_max)
!
    character(len=16) :: sat_name(n_groups_max)
    character(len=16) :: d_type(n_groups_max)
!
    integer, parameter :: luin=11
    integer, parameter :: luout=12
    character(len=18) :: c_name(n_chan_max,n_groups_max)
    character(len=1)  :: adum1
    character(len=70) :: file_table_1
    character(len=70) :: file_table_2
!
    print *,' '
    print *,'BEGIN PROGRAM reformat R tables'
!
! Read arguments
    argc = iargc() 
    if (argc .ne. 2) then
      print *,'useage must be: compare_R_tables.x TabIN TaOUT'
      stop
    endif
    call GetArg( 1_4, file_table_1)
    call GetArg( 2_4, file_table_2)
!
! Read modified GEOS table 
    open(unit=luin,file=trim(file_table_1),form='formatted')
    print ('(3a,i4)'),' input file=',trim(file_table_1), &
              ' opened on unit=',luin
    read (luin,'(7x,i2,30x,i2)') ihead,igroups  
    print *,'Table 1: ihead,igroups = ',ihead,igroups
    do i=1,ihead                     ! skip these header records
      read (luin,'(a1)') adum1
    enddo
!
    do i=1,igroups   
      read (luin,'(i2)') inum(i)
      read (luin,'(1x,a16,3x,a16,1x,i4)') &
              d_type(i),sat_name(i),n_channels(i)
      print ('(a,i3,1x,2a16,i5)'),'i,dtype,sat,nchan=', &
              inum(i),d_type(i),sat_name(i),n_channels(i)
      do j=1,n_channels(i)    
        read (luin,'(1x,a18,2i5,5f8.3)') &
                c_name(j,i),n_chan(j,i),n_use(j,i),rvalues(1:5,j,i)
        if (allsky) then
          rvalues(1,j,i)=err_factor*rvalues(2,j,i)  ! scale and move value
        else
          rvalues(1,j,i)=err_factor*rvalues(1,j,i)  ! scale value
        endif
      enddo
!
    enddo     
    close (luin)
    print *,'Table 1 read'
!
! Write reformatted table
!
    open(unit=luout,file=trim(file_table_2),form='formatted')
    print ('(3a,i4)'),' output file=',trim(file_table_2), &
              ' opened on unit=',luout
    write (luout,'(a7,i2,a30,i2,a1)') 'header=',ihead, &
                  ' more lines. Number of groups=',igroups,' '
    write (luout,'(a58)') &   
                  '123456789012345678x1234xxx1212345671234567123456712345678'
    write (luout,'(a)') &  
                  ' sensor/instr/sat   chan iuse  error  ermax  var_b  var_pg'
!
    do i=1,igroups   
      write (luout,'(i2)') inum(i)
      write (luout,'(1x,a16,3x,a16,1x,i4)') &
              d_type(i),sat_name(i),n_channels(i)
      do j=1,n_channels(i)    
        write (luout,'(1x,a18,2i5,2f7.3,f7.2,f8.3)') &
                c_name(j,i),n_chan(j,i),n_use(j,i),rvalues(1,j,i),rvalues(3:5,j,i)
      enddo
    enddo     
    write (luout,'(a1)') ' ' 
    write (luout,'(a)') &  
                  ' sensor/instr/sat   chan iuse  error  ermax  var_b  var_pg'
    close (luout)
    print *,'Table written'
!
    end program reformat_R_tables
