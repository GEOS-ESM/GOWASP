    program make_new_table
!
! Makes new table of conventional R values by interpolating values from
! estimates of R produced by estimate_Rr_conv.f90. 
!
    implicit none
!
    logical :: luvavg ! true if u and v input stats on 2 files to be averaged.
!
    integer, parameter :: unit_R_est=10
    integer, parameter :: unit_rc=11
    integer, parameter :: luin=12
    integer, parameter :: luout=13
    integer, parameter :: n_tab_values=33
    integer :: i,j,k,k1,n,ng
    integer :: j_above,j_below
    integer :: n_kx
    integer :: ngroups,nheader
    integer :: n_est_values
    integer :: argc
    integer :: kx, kx_read
    integer :: kx_list(10)
    integer :: n_col_tab
    integer :: clength
!
    real    :: R_est2,p_est2,V_est
    real    :: plev
    real    :: weight
    real    :: err_tab(6,n_tab_values,100:299) 
    real    :: R_new(n_tab_values)
    real, allocatable :: R_est(:)
    real, allocatable :: p_est_levs(:)
!
    character(len=20)  :: c_names(2)
    character(len=1)   :: c_dum
    character(len=5)   :: c_tail
    character(len=80)  :: file_R_est
    character(len=120) :: file_table_old
    character(len=120) :: file_table_new
    character(len=120) :: file_rc
    character(len=3)   :: ckx
    character(len=16)  :: instr_name
    character(len=16)  :: field_name
!   
    print *,' '
    print *,'BEGIN PROGRAM make_new_table'
!
! Read arguments
    argc = iargc()
    if (argc .ne. 3) then
      print *,'useage must be: prog.x file_rc file_table_old file_table_new'
      stop
    endif
    call GetArg( 1_4, file_rc)
    call GetArg( 2_4, file_table_old)
    call GetArg( 3_4, file_table_new)
!
! Read table of obs error standard deviations for conventional obs.
    open (unit=luin,file=trim(file_table_old),form='formatted')
    print ('(3a,i4)'),'Input file=',trim(file_table_old), &
              ' opened on unit=',luin
    do i=100,299
      read (luin,'(i4)') kx_read
      do j=1,n_tab_values
        read (luin,'(1x,6e12.5)') err_tab(:,j,kx_read)
      enddo  
    enddo
    print *,'Old R table read'
    close (luin)
!
! Read resource file
    open (unit=unit_rc,file=trim(file_rc),form='formatted')   
    print ('(a)'),' Input file opened = ',trim(file_rc)
    read (unit_rc,*) c_names(1),ngroups,c_names(2),nheader
    do i=1,nheader
      read (unit_rc,'(a1)') c_dum
    enddo 
    read (unit_rc,*) c_names(1),luvavg
!
    do ng=1,ngroups
      read (unit_rc,'(a1)') c_dum
      read (unit_rc,*) c_names(1:3),n_kx
      print ('(a,i3,3(2x,a),i2)'),'Group=',ng,trim(c_names(1)), &
                                trim(c_names(2)),' Number of copies=',n_kx 
      read (unit_rc,'(a)') file_R_est
      read (unit_rc,*) kx_list(1:n_kx)         
!
! Read file of estimated values 
      open(unit_R_est,file=trim(file_R_est),form='formatted')
      print ('(3a,i4)'),' input file=',trim(file_R_est), &
                        ' opened on unit=',unit_R_est
      read (unit_R_est,'(57x,2(2x,a16))') field_name,instr_name
      read (unit_R_est,'(54x,i5)') n_est_values
      ckx=instr_name(1:3)
      read (ckx,'(i3)') kx
      print ('(a,i3,3a,i3)'),'kx=',kx,'  field=',trim(field_name), &
                             '  levels=',n_est_values
!
      allocate (R_est(n_est_values))
      allocate (p_est_levs(n_est_values))
      do i=1,7  ! skip header
        read (unit_R_est,'(a1)') c_dum
      enddo
      do n=1,n_est_values
        read (unit_R_est,'(4x,f12.2,21x,e13.3)') p_est_levs(n), R_est(n)
      enddo
      close (unit_R_est)
!
!  If file is for u or v and table is to have avg value for i and v, then ...
!
      clength=len(trim(file_R_est))     
      if (luvavg .and. (file_R_est(clength-5:clength) == '_u.txt' .or. &
                        file_R_est(clength-5:clength) == '_v.txt')) then
        if (file_R_est(clength-5:clength) == '_u.txt') then 
          file_R_est(clength-5:clength) = '_v.txt' 
        else
          file_R_est(clength-5:clength) = '_u.txt' 
        endif
!
! Read 2nd file of estimated values 
        open(unit_R_est,file=trim(file_R_est),form='formatted')
        print ('(3a,i4)'),' input file=',trim(file_R_est), &
                          ' opened on unit=',unit_R_est
        read (unit_R_est,'(57x,2(2x,a16))') field_name,instr_name
        read (unit_R_est,'(54x,i5)') n_est_values
        ckx=instr_name(1:3)
        read (ckx,'(i3)') kx
        print ('(a,i3,3a,i3)'),'kx=',kx,'  field=',trim(field_name), &
                               '  levels=',n_est_values
!
        do i=1,7  ! skip header
          read (unit_R_est,'(a1)') c_dum
        enddo
        do n=1,n_est_values
          read (unit_R_est,'(4x,f12.2,21x,e13.3)') p_est2,R_est2
          V_est=0.5*(R_est2**2 + R_est(n)**2) 
          if (V_est > 0.) then
            R_est(n)=sqrt(V_est)    ! replace by sqrt avg var of u and v 
          endif                        ! otherwise leave R_est(n) unchanged
        enddo
        close (unit_R_est)
!
      endif  ! end of consideration of 2nd (complementary) file
!
      if (trim(field_name) == 'T') then
        n_col_tab=2
      else if (trim(field_name) == 'q') then
        n_col_tab=3
      else if (trim(field_name) == 'u') then
        n_col_tab=4
      else if (trim(field_name) == 'v') then
        n_col_tab=4
      else if (trim(field_name) == 'p') then
        n_col_tab=5
      endif
!
! Compute new R values in table for selected kx
      do i=1,n_kx
        if ((kx_list(i) < 300) .and. (kx_list(i) > 99) ) then
          do j=1,n_tab_values
            plev=err_tab(1,j,kx_list(i))
            j_above=0
            j_below=0 
            do k=1,n_est_values
              if ((R_est(k) > 0.) .and. (plev < p_est_levs(k)) ) then
                j_below=k
              endif  
              k1=n_est_values+1-k
              if ((R_est(k1) > 0.) .and. (plev >= p_est_levs(k1)) ) then
                j_above=k1
              endif  
            enddo
!        
            if ((j_above == 0) .and. (j_below > 0)) then
              err_tab(n_col_tab,j,kx_list(i))=R_est(j_below)
            else if ((j_above > 0) .and. (j_below == 0)) then
              err_tab(n_col_tab,j,kx_list(i))=R_est(j_above)  
            else if ((j_above > 0) .and. (j_below > 0)) then
              if (j_above > j_below) then  ! value can be interpolated
                weight=(p_est_levs(j_below)-plev)/ &
                       (p_est_levs(j_below)-p_est_levs(j_above))
                err_tab(n_col_tab,j,kx_list(i))=   &
                       weight*R_est(j_above)+(1.-weight)*R_est(j_below)
              else                          ! value at level
                err_tab(n_col_tab,j,kx_list(i))=R_est(j_above)
              endif
            endif                           ! else, old value kept
!  
          enddo  ! loop over p-levels in table
        endif
      enddo      ! loop over selected kx in current group    
!
     deallocate (R_est,p_est_levs)
    enddo ! loop over kx groups specified in rc file
    close (unit_rc)
!
! Write new file of R values
    open(unit=luout,file=trim(file_table_new),form='formatted')
    print ('(3a,i4)'),' outputfile=',trim(file_table_new), &
              ' opened on unit=',luout
    do i=100,299
      write (luout,'(i4,a17)') i,' OBSERVATION TYPE'
      do j=1,n_tab_values
        write (luout,'(1x,6e12.5)') err_tab(:,j,i)
      enddo  
    enddo
    close (luout)
!
    print *,'Conv Error Table read and copied'  
!
    end program make_new_table 
