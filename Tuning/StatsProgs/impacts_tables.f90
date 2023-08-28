  program impact_tables
!
! Program for constructing summary tables of obs impact estimates
! from previously produced output from impacts_rad and impacts_conv 
! as instructed by a resource file. 
!
  implicit none
!  
  logical :: levels
  logical :: leof     ! eof found in file
!
  integer, parameter :: ndatamax=100   ! max number allowed obs types im rcfile
  integer, parameter :: nlevmax=700    ! max number of channels, p, or z levels
  integer(8), parameter :: izero=0_8
  integer :: n,k,id,nt,i,j
  integer :: nrcdata
  integer :: ngroups
  integer :: nlevs
  integer :: tot_count
  integer :: itypes
  integer :: ntimes
  integer :: ntimes_read(ndatamax)
  integer :: ikx,ikt
  integer :: idgroups(ndatamax,2)
  integer :: idtypes(ndatamax,4)   
  integer :: itab(nlevmax,2)
  integer(8) :: itable(nlevmax,3,ndatamax)
!
! idtypes(:,1) is to which group a type of obs is assigned
! idtypes(:,2) is this obs type (1) rad, (2) conv, or (3) gpsro
! idtypes(:,2) is 0 or conv kx value (300 is for gpsro)
! idtypes(:,4) is 0 or if conv then 1 for T or u, 2 for q or v, 3 for ps or w
! idgroups(:,1) like idtypes(:,2) but for a group of obs types
! idgroups(:,2) is number of levels or channels for a group
!
  real :: pc  
  real :: tot_impact
  real :: fracneg
  real :: plev(nlevmax) 
  real :: plevs_conv(nlevmax) , plevs_gps(nlevmax) 
  real :: rtable(nlevmax,3,ndatamax)
  real :: rtab(nlevmax,2)
!
  character (len=220) :: file_rc
  character (len=220) :: filepath
  character (len=220) :: filename
  character (len=20)  :: rcdata(ndatamax,3) 
  character (len=20)  :: cgroups(ndatamax)
  character (len=20)  :: cfield(ndatamax) 
  character (len=40)  :: cexpname
  character (len=40)  :: cdum
  character (len=20)  :: c_kx, c_kt 
  character (len=20)  :: c_date, c_time
  character (len=20)  :: ctrim1, ctrim2
  character (len=20), parameter  :: file_out='file_table_output'
!
! Get arguments
  call GetArg( 1_4, file_rc)
!
  idtypes(:,:)=0
  idgroups(:,:)=0
  cfield(:)=' '
  cgroups(:)=' '
  itable(:,:,:)=izero
  rtable(:,:,:)=0.
!
! Read resource file
  open (unit=10,file=file_rc,form='formatted')
  read (10,'(a10,a)') cdum,filepath
  read (10,*) cdum,cexpname
  read (10,*) cdum,ntimes
  read (10,*) cdum,levels
  do n=1,100
    read (10,*) rcdata(n,:)
    if (trim(rcdata(n,1)) == 'EOF') exit
    nrcdata=n
  enddo 
  close (10)
!
! Determine number of distinct groups  
  ngroups=1
  cgroups(1)=rcdata(1,3)
  idtypes(1,1)=1
  do n=2,nrcdata
    id=0
    do k=1,ngroups
      if (trim(rcdata(n,3)) == trim(cgroups(k))) then 
        id=k
      endif
    enddo
    if (id == 0) then
      id=ngroups+1 
      ngroups=id
      cgroups(id)=rcdata(n,3)  ! name of newly found distinct group in rcfile
    endif
    idtypes(n,1)=id            ! index of group containing rcfile type n
  enddo
!
  print *,'number of groups found = ',ngroups
!
! Determine wether each obs type in rcfile is rad, conv, or gps
! For conv, set field index based on field name 
  do n=1,nrcdata
    idtypes(n,2)=1    ! rad obs default
    ctrim1=trim(rcdata(n,1))
    ctrim2=trim(rcdata(n,2))    
    if (ctrim1(1:1) == '1' .or. ctrim1(1:1) == '2') then
      idtypes(n,2)=2  ! conv obs
      if (ctrim2 == 'T') then 
        cfield(n)='1'
      elseif (ctrim2 == 'q') then 
        cfield(n)='2'
      elseif (ctrim2 == 'ps') then 
        cfield(n)='3'
      elseif (ctrim2 == 'u') then 
        cfield(n)='1'
      elseif (ctrim2 == 'v') then 
        cfield(n)='2'
      elseif (ctrim2 == 'w') then 
        cfield(n)='3'
      endif
    elseif (ctrim1(1:1) == '3' ) then 
      idtypes(n,2)=3  ! gps obs
      cfield(n)='1'
    endif
!
! For conv or gps obs, interpret kx and cfield as integers
    if (idtypes(n,2) > 1) then 
      read (ctrim1,*) idtypes(n,3)
      read (cfield(n),*) idtypes(n,4)
    endif 
!
  enddo   ! loop over types of obs in rcfile
!
  ntimes_read(:)=0 ! counter for number of times read

!   
! Read previous impact tables for each requested obs type in rc file 
  do n=1,nrcdata
!
    leof=.false.
!
    if (idtypes(n,2) == 1) then ! rad obs
      filename=trim(filepath)//'/impacts_'//trim(cexpname)//'_'// &
               trim(rcdata(n,1))//'_'//trim(rcdata(n,2))//'.txt'
      open (unit=10,file=trim(filename),form='formatted')
      do nt=1,ntimes
        if (.not. leof) then 
          read (10,'(a4)') cdum
          if (cdum(1:4) == '#EOF') then
            leof=.true.
          endif
        endif
        if (.not. leof) then 
!
          read (10,'(a38,i5)') cdum,nlevs
          read (10,'(a1)') cdum
          read (10,'(a1)') cdum
!
          do k=1,nlevs+1
            read (10,'(i5,2i8,1p2e14.4)') j,itab(k,1:2),rtab(k,1:2)
          enddo
!
          id=idtypes(n,1)
          ntimes_read(id)=ntimes_read(id)+1
          idgroups(id,1)=1
          idgroups(id,2)=max(nlevs,idgroups(id,2))
          do k=1,nlevs+1
            itable(k,1:2,id)=itable(k,1:2,id)+itab(k,1:2)
            rtable(k,1,id)=rtable(k,1,id)+rtab(k,1)
            rtable(k,2,id)=rtable(k,2,id)+itab(k,1)*rtab(k,2)**2 ! sum of sqrs
          enddo
!
        endif  ! test on EOF
      enddo
      close (10)
!     
    else      
!
! Obs type is conv or gpsro
      filename=trim(filepath)//'/impacts_'//trim(cexpname)//'_conv.txt'
      open (unit=10,file=trim(filename),form='formatted')
      do nt=1,ntimes
        if (.not. leof) then 
          read (10,'(a8,i5)') cdum,itypes
          if (cdum(1:4) == '#EOF') then
            leof=.true.
          endif
        endif
        if (.not. leof) then 
!
          read (10,'(a30,i5)') cdum,nlevs 
          read (10,'(a1)') cdum
          read (10,'(a1)') cdum
          id=idtypes(n,1)
          ntimes_read(id)=ntimes_read(id)+1
          idgroups(id,2)=max(nlevs,idgroups(id,2))
          do i=1,itypes
            read (10,'(2a4,2a9)') c_kx,c_kt,c_date,c_time
            read (c_kx,*) ikx          
            read (c_kt,*) ikt
            do k=1,nlevs+1
              read (10,'(i5,f8.0,2i8,1p2e14.4)') j,plev(k), &
                                   itab(k,1:2),rtab(k,1:2)
            enddo
            if (ikx == idtypes(n,3) .and. ikt == idtypes(n,4)) then
              if (ikx == 300) then   ! gpsro
                idgroups(id,1)=3
                plevs_gps(1:nlevs+1)=plev(1:nlevs+1)
              else                   ! conv
                idgroups(id,1)=2 
                plevs_conv(1:nlevs+1)=plev(1:nlevs+1)
              endif
              do k=1,nlevs+1
                itable(k,1:2,id)=itable(k,1:2,id)+itab(k,1:2)
                rtable(k,1,id)=rtable(k,1,id)+rtab(k,1)
                rtable(k,2,id)=rtable(k,2,id)+itab(k,1)*rtab(k,2)**2 ! sum sqrs
              enddo
            endif
          enddo
        endif
!
      enddo      
      close (10)
!
    endif 
  enddo
!
! Compute time means and rms from mean square
  do id=1,ngroups
    nlevs=idgroups(id,2)
    do k=1,nlevs+1
      itable(k,1:2,id)=itable(k,1:2,id)/ntimes
      rtable(k,1:2,id)=rtable(k,1:2,id)/ntimes
      if (itable(k,1,id) > izero) then
        rtable(k,3,id)=rtable(k,1,id)/itable(k,1,id)
        if (rtable(k,2,id) > 0) then
          rtable(k,2,id)=sqrt(rtable(k,2,id)/itable(k,1,id))
        endif
      endif
    enddo
  enddo
!
! Write table of sums over channels or levels
!
  open (unit=10,file=trim(file_out),form='formatted')   
  write (10,'(2(a,i5))') 'ngroups= ',ngroups,'  ntimes= ',ntimes

  write (10,'(2a)') '  k               group    count  negvals', &   
                    '      impact   impact/ob      rmsO-F   frac-'
  tot_impact=0.
  tot_count=0
  do id=1,ngroups
    k=idgroups(id,2)+1
    if (itable(k,1,id) > 0_8) then
      fracneg=real(itable(k,2,id))/real(itable(k,1,id))
    else
      fracneg=0.
    endif
    tot_impact=tot_impact+rtable(k,1,id)
    tot_count=tot_count+itable(k,1,id)
    write (10,'(i3,a20,2i9,1p3e12.3,0pf8.3)') id,trim(cgroups(id)),      &
           itable(k,1:2,id),rtable(k,1,id),rtable(k,3,id),rtable(k,2,id), &
           fracneg
  enddo
  write (10,'(a5,18x,i9,9x,1pe12.3)') 'total',tot_count,tot_impact

! Write table of levels or channels
  if (levels) then 
    do id=1,ngroups
      write (10,'(a)') ' '
      write (10,'(i3,a20,i3)') id,trim(cgroups(id)),idgroups(id,2)
      write (10,'(2a)') '  k LevChan    count  negvals', &
                        '      impact   impact/ob      rmsO-F   frac-'
      do k=1,idgroups(id,2)
        pc=real(k)
        if (idgroups(id,1) == 2) then
          pc=plevs_conv(k)
        elseif (idgroups(id,1) == 3) then 
          pc=plevs_gps(k)
        endif
        if (itable(k,1,id) > 0_8) then
          fracneg=real(itable(k,2,id))/real(itable(k,1,id))
        else
          fracneg=0.
        endif
        write (10,'(i3,f8.0,2i9,1p3e12.3,0pf8.3)') k,pc,itable(k,1:2,id), &
               rtable(k,1,id), rtable(k,3,id),rtable(k,2,id),fracneg
      enddo 
    enddo
  endif
!
  write (10,'(a)') ' '
  close (10)
!
  end program impact_tables
