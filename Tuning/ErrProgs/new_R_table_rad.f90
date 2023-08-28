    program make_new_table
!
    implicit none
!
    logical, parameter :: passive=.false. ! change values for passive channels
! Once the R-values for passive channels have been changed by this routine,
! they should probably not be changed again since that would just reduce them
! by an arbitrry factor further.
    integer, parameter :: unit_R=10
    integer, parameter :: luin=11
    integer, parameter :: luout=12
    integer, parameter :: max_chans=700
    integer :: n_channels
    real, parameter :: pchan_airs=0.45, pchan_other=0.7
    real    :: err_tab1(4) 
    real    :: R_n(max_chans)
    character(len=16) :: sat_name
    character(len=16) :: d_type
!
    integer :: argc
    integer :: i, j, n
    integer :: i_type
    integer :: inum
    integer :: ihead   ! number of header records to skip in sat_error file
    integer :: igroups ! number of sat platforms or instruments in table
    integer :: n_chan
    integer :: n_use
    integer :: n_channels_tab
    character(len=19) :: c_name
    character(len=16)  :: d_type_tab
    character(len=16)  :: sat_name_tab
    character(len=1)  :: adum1
    character(len=120) :: file_new_R
    character(len=120) :: file_table_old
    character(len=120) :: file_table_new
    character(len=16)  :: c_space
!
    print *,' '
    print *,'BEGIN PROGRAM make_new_table'
    c_space='                '
!
! Read arguments
    argc = iargc()
    if (argc .ne. 3) then
      print *,'useage must be: new_R_table.x file_new_R file_table_old file_table_new'
      stop
    endif
    call GetArg( 1_4, file_new_R)
    call GetArg( 2_4, file_table_old)
    call GetArg( 3_4, file_table_new)

    open(unit_R,file=trim(file_new_R),form='formatted')
    print ('(3a,i4)'),' input file=',trim(file_new_R),' opened on unit=',unit_R
!    
    read (unit_R,'(41x,a16,2x,a16)'),d_type,sat_name
    read (unit_R,'(54x,i5)'),n_channels
    print ('(5a,i5)'),'instrument_type=',trim(d_type),'  satellite_name=', &
                      trim(sat_name),'  number_channels=',n_channels
    do i=1,7  ! skip header
      read (unit_R,'(a1)') adum1
    enddo
    do n=1,n_channels
      read (unit_R,'(37x,e13.3)') R_n(n)
    enddo
    close (unit_R)
!
! Read table of obs error standard deviations for conventional obs.
    open(unit=luin,file=trim(file_table_old),form='formatted')
    open(unit=luout,file=trim(file_table_new),form='formatted')
    print *,' '
    print ('(3a,i4)'),' input file=',trim(file_table_old), &
              ' opened on unit=',luin
    print ('(3a,i4)'),' outputfile=',trim(file_table_new), &
              ' opened on unit=',luout
!
! Read header: 
! igroups=number of distinct
    read (luin,'(7x,i2,30x,i2)') ihead,igroups  
    write (luout,'(a7,i2,a30,i2,a1)') 'header=',ihead, &
                  ' more lines. Number of groups=',igroups,' '
    write (luout,'(a58)') &   
                  'x123456789012345678x1234xxx1212345671234567123456712345678'
    write (luout,'(a58)') &  
                  ' sensor/instr/sat   chan iuse  error  ermax  var_b  var_pg'
    do i=1,ihead                     ! skip these header records
      read (luin,'(a1)') adum1
    enddo
!
! Check the following records and skip or read as required
    do i=1,igroups   
      read (luin,'(i2)') inum
      read (luin,'(1x,a16,1x,a16,1x,i4)') &
              d_type_tab,sat_name_tab,n_channels_tab
      write (luout,'(i2,a)') inum,c_space
      write (luout,'(1x,a16,1x,a16,1x,i4)') &
              d_type_tab,sat_name_tab,n_channels_tab

      if (trim(d_type_tab) == trim(d_type) .and. &
          trim(sat_name_tab) == trim(sat_name)) then
        print *,'Table being changed for ',trim(d_type),' ',trim(sat_name)
!
        do j=1,n_channels_tab
          read (luin,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
              c_name,n_chan,n_use,err_tab1
          if (R_n(j) > 0. .and. err_tab1(1) < 999.) then
            err_tab1(1)=R_n(j)
          else if (passive .and. (err_tab1(1) < 999.)) then ! change passive channel 
            err_tab1(1)=err_tab1(1)*pchan_other
          endif
          if (err_tab1(1) > 999.) then 
            write (luout,'(1x,a19,i4,i5,2f7.0,f7.2,f8.3)') &
                c_name,n_chan,n_use,err_tab1
          else
            write (luout,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
              c_name,n_chan,n_use,err_tab1
          endif
        enddo
!
      else  ! copy these records        
        do j=1,n_channels_tab
          read   (luin,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
              c_name,n_chan,n_use,err_tab1
          if (err_tab1(1) > 999.) then 
            write (luout,'(1x,a19,i4,i5,2f7.0,f7.2,f8.3)') &
               c_name,n_chan,n_use,err_tab1
          else
            write (luout,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
              c_name,n_chan,n_use,err_tab1
          endif
        enddo
!
      endif  ! test on d_type  
!
    enddo    ! loop over all groups of entries in table   
!
    close (luin)
!
    write (luout,'(a1)') ' ' 
    write (luout,'(a58)') &  
                  ' sensor/instr/sat   chan iuse  error  ermax  var_b  var_pg'
    close (luout)
!
    print *,'Sat Error Table read and copied'  
!
    end program make_new_table 
