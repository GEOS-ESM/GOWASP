    program compare_C_tables
!
! Compare 2 R tables for conv obs
! Compare second file with first as Table2/Table1
!
    implicit none
!
    integer, parameter :: unit_R_est=10
    integer, parameter :: unit_rc=11
    integer, parameter :: luin=12
    integer, parameter :: luout=13
    integer, parameter :: n_est_values=10
    integer, parameter :: n_tab_values=33
    integer :: i,j,k,k1,n,ng
    integer :: j_above,j_below
    integer :: n_kx
    integer :: ngroups,nheader
    integer :: n_levels
    integer :: argc
    integer :: kx, kx_read
    integer :: kx_list(10)
    integer :: n_col_tab
!
    real    :: plev
    real    :: weight
    real    :: err_tab(6,n_tab_values,100:299,3) 
    real    :: R_new(n_tab_values)
    real    :: R_est(n_est_values)
    real    :: p_est_levs(n_est_values)
!
    character(len=12)  :: c_group
    character(len=1)   :: c_dum
    character(len=70) :: file_R_est
    character(len=70) :: file_table_1
    character(len=70) :: file_table_2
    character(len=70) :: file_table_3
    character(len=1)   :: c_field
!   
    data p_est_levs(1:5)  /1000., 925., 850., 700., 500./
    data p_est_levs(6:10) / 400., 300., 200., 100.,  50./
!
    print *,' '
    print *,'BEGIN PROGRAM make_new_table'
!
! Read arguments
    argc = iargc()
    if (argc .ne. 3) then
      print *,'useage must be: compare_C_tables.x Tab1 Tab2 Tab2/Tab1'
      stop
    endif
    call GetArg( 1_4, file_table_1)
    call GetArg( 2_4, file_table_2)
    call GetArg( 3_4, file_table_3)
!
! Read table 1 of obs error standard deviations for conventional obs.
    open (unit=luin,file=trim(file_table_1),form='formatted')
    print ('(3a,i4)'),'Input file=',trim(file_table_1), &
              ' opened on unit=',luin
    do i=100,299
      read (luin,'(i4)') kx_read
      do j=1,n_tab_values
        read (luin,'(1x,6e12.5)') err_tab(:,j,kx_read,1)
      enddo  
    enddo
    print *,'Tables 1 read'
    close (luin)
!
! Read table 2 of obs error standard deviations for conventional obs.
    open (unit=luin,file=trim(file_table_2),form='formatted')
    print ('(3a,i4)'),'Input file=',trim(file_table_2), &
              ' opened on unit=',luin
    do i=100,299
      read (luin,'(i4)') kx_read
      do j=1,n_tab_values
        read (luin,'(1x,6e12.5)') err_tab(:,j,kx_read,2)
      enddo  
    enddo
    print *,'Tables 2 read'
    close (luin)
! 
! Compute differences
    do i=100,299
      do j=2,n_tab_values
       err_tab(1,j,i,3)=err_tab(1,j,i,1) 
       err_tab(2:6,j,i,3)=err_tab(2:6,j,i,2)/err_tab(2:6,j,i,1) 
      enddo  
    enddo
    print *,'Table ratios computed'
!
! Write new file of R values
    open(unit=luout,file=trim(file_table_3),form='formatted')
    print ('(3a,i4)'),' outputfile=',trim(file_table_3), &
              ' opened on unit=',luout
    do i=100,299
      write (luout,'(i4,a17)') i,' OBSERVATION TYPE'
      do j=1,n_tab_values
        write (luout,'(1x,6e12.5)') err_tab(:,j,i,3)
      enddo  
    enddo
    close (luout)
!
    print *,'Table of ratios created'
!
    end program compare_C_tables
