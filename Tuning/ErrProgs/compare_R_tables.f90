    program compare_R_tables
!
! Compare 2 R tables for radiance obs
! Compare second file with first as Table2/Table1
!
    implicit none
!
    integer, parameter :: n_chan_max=600
    integer, parameter :: n_groups_max=30
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
    integer :: n_chan2
!
    real(4) :: rvalues(4,n_chan_max,n_groups_max,2)
    real(4) :: ratio
!
    character(len=4) :: sat_name(n_groups_max)
    character(len=5) :: d_type(n_groups_max)
    character(len=9) :: instr_name
!
    integer, parameter :: luin=11
    integer, parameter :: luout=12
    character(len=19) :: c_name(n_chan_max,n_groups_max), c_name2
    character(len=1)  :: adum1
    character(len=70) :: file_table_1
    character(len=70) :: file_table_2
    character(len=70) :: file_table_3
!
    print *,' '
    print *,'BEGIN PROGRAM compare R tables'
!
! Read arguments
   argc = iargc() 
   if (argc .ne. 3) then
      print *,'useage must be: compare_R_tables.x Tab1 Tab2 Tab2/Tab1'
      stop
    endif
    call GetArg( 1_4, file_table_1)
    call GetArg( 2_4, file_table_2)
    call GetArg( 3_4, file_table_3)
!
! Read table 1
    open(unit=luin,file=trim(file_table_1),form='formatted')
    print ('(3a,i4)'),' input file=',trim(file_table_1), &
              ' opened on unit=',luin
    read (luin,'(7x,i2,30x,i2)') ihead,igroups  
    print *,'Table 1: ihead,igroups = ',ihead,igroups
    do i=1,ihead                     ! skip these header records
      read (luin,'(a1)') adum1
    enddo
    do i=1,igroups   
      read (luin,'(i2)') inum(i)
      read (luin,'(1x,a5,1x,a4,1x,i3)') &
              d_type(i),sat_name(i),n_channels(i)
      do j=1,n_channels(i)    
        read (luin,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
                c_name(j,i),n_chan(j,i),n_use(j,i),rvalues(1:4,j,i,1)
      enddo
    enddo     
    close (luin)
    print *,'Table 1 read'
!
! Read table 2 
    open(unit=luin,file=trim(file_table_2),form='formatted')
    print ('(3a,i4)'),' input file=',trim(file_table_2), &
              ' opened on unit=',luin
    read (luin,'(7x,i2,30x,i2)') ihead,igroups 
    print *,'Table 2: ihead,igroups = ',ihead,igroups
    do i=1,ihead                     ! skip these header records
      read (luin,'(a1)') adum1
    enddo
    do i=1,igroups   
      read (luin,'(i2)') inum(i)
      read (luin,'(1x,a5,1x,a4,1x,i3)') &
              d_type(i),sat_name(i),n_channels(i)
      do j=1,n_channels(i)    
        read (luin,'(1x,a19,i4,i5,2f7.3,f7.2,f8.3)') &
                c_name2,n_chan2,n_use(j,i),rvalues(1:4,j,i,2)
        if ((c_name(j,i) /= c_name2) .or. (n_chan(j,i) /= n_chan2)) then
          print *,'c_name1,c_name2,n_chan1,n_chan2 = ', &
                   c_name(j,i),c_name2,n_chan(j,i),n_chan2
        endif 
      enddo
    enddo     
    close (luin)
    print *,'Table 2 read'
!
! Write table 3 
    open(unit=luout,file=trim(file_table_3),form='formatted')
    print ('(3a,i4)'),' output file=',trim(file_table_3), &
              ' opened on unit=',luout
    write (luout,'(a7,i2,a30,i2,a1)') 'header=',ihead, &
                  ' more lines. Number of groups=',igroups,' '
    write (luout,'(a58)') &   
                  '123456789012345678x1234xxx1212345671234567123456712345678'
    write (luout,'(a58)') &  
                  'sensor/instr/sat   chan iuse   GsiTab   SimErr      S/G  '
    do i=1,igroups   
      write (luout,'(i2)') inum(i)
      write (luout,'(1x,a5,1x,a4,1x,i3)') &
              d_type(i),sat_name(i),n_channels(i)
      do j=1,n_channels(i)    
        ratio=rvalues(1,j,i,2)/rvalues(1,j,i,1)
        write (luout,'(1x,a19,i4,i5,3f9.3)')    &
                c_name(j,i),n_chan(j,i),n_use(j,i), &
                rvalues(1,j,i,1),rvalues(1,j,i,2),ratio
      enddo
    enddo     
    write (luout,'(a1)') ' ' 
    write (luout,'(a58)') &  
                  'sensor/instr/sat   chan iuse   GsiTab   SimErr      S/G  '
    close (luout)
    print *,'Table 3 written'
!
    end program compare_R_tables
