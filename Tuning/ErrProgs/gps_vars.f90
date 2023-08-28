    program gps_err_var
!
!  Construct new GPSRO error table from Vcorr stats for GPSRO
!
    implicit none
    integer, parameter :: line_gps=1536   ! line for table 1 field=G label in stats file
    integer, parameter :: line_lats=7     ! line having lat info in stats files
    integer, parameter :: nlevs_stats=90  ! number of z levels in table 1
    integer, parameter :: nlevs_tab=102   ! number + 1 of z levels in err table
    integer, parameter :: max_ifiles=12   ! maximum number of lat zones
    integer, parameter :: iunit=10
    integer, parameter :: m1_min=5         ! do not consider level < m1_min
    character(len=*), parameter :: ctype='var'
!
    integer :: i,n,m1,m2,mx
    integer :: ier
    integer :: ifiles
    integer :: argc
    integer(4) :: iargc
!
    real(4), allocatable :: stats_osse(:,:)
    real(4), allocatable :: stats_real(:,:)
    real(4), allocatable :: stats_diff(:,:)
    real(4), allocatable :: stats_levs(:)
    real(4), allocatable :: tab_old(:,:)
    real(4), allocatable :: tab_new(:,:)
    real(4), allocatable :: tab_levs(:)
    real(4), allocatable :: lats(:,:)
    real(4) :: w1, w2, wd
    real(4) :: var
    real(4), parameter :: var_min=1.e-16
!
    integer(4) :: iarg
    character(len=200) :: c_osse_omf
    character(len=200) :: c_real_omf
    character(len=8), allocatable :: c_names(:)
    character(len=20)  :: c_endname
    character(len=200) :: c_tab_old
    character(len=200) :: c_tab_new
    character(len=200) :: c_file
! 
! Read arguments
    argc = iargc()
    ifiles=argc-5
    if (ifiles > 12) then
      print *,'ifiles=',ifiles,' too large'
      stop
    endif 
    allocate (c_names(ifiles))
!
    call GetArg( 1_4, c_real_omf) ! 
    call GetArg( 2_4, c_osse_omf) ! 
    call GetArg( 3_4, c_endname)  ! 
    call GetArg( 4_4, c_tab_old)   ! table of old error factors
    call GetArg( 5_4, c_tab_new)  ! table of new error factors
    do i=1,ifiles
      iarg=i+5
      call GetArg( iarg, c_names(i))  
    enddo
!
    allocate (stats_osse(nlevs_stats,ifiles))
    allocate (stats_real(nlevs_stats,ifiles))
    allocate (stats_diff(nlevs_stats,ifiles))
    allocate (stats_levs(nlevs_stats))
    allocate (tab_old(nlevs_tab,max_ifiles))
    allocate (tab_new(nlevs_tab,max_ifiles))
    allocate (tab_levs(nlevs_tab))
    allocate (lats(3,ifiles))
    tab_new(:,:)=0.
!
    do i=1,ifiles
      c_file=trim(c_real_omf)//trim(c_names(i))//trim(c_endname)        
      call read_vcorr_stats (nlevs_stats,stats_levs,stats_real(:,i), &
          lats(:,i),c_file,line_lats,line_gps,ctype)
      c_file=trim(c_osse_omf)//trim(c_names(i))//trim(c_endname)        
      call read_vcorr_stats (nlevs_stats,stats_levs,stats_osse(:,i), &
          lats(:,i),c_file,line_lats,line_gps,ctype)
    enddo
!
    print *,' '
    print *,'Lats in stats files'
    do n=1,3
      print ('(12f7.2)'),(lats(n,i),i=1,ifiles)
    enddo
!
    ier=0
    do i=2,ifiles
      if (lats(2,i) >= lats(2,i-1)) then
        print *,'*** lat order problem: ',lats(2,i),' >= ',lats(2,i-1)  
        ier=ier+1
      endif
    enddo
    if (ier > 0) stop    
!
    stats_diff(:,:)=stats_real(:,:)-stats_osse(:,:)    
    call adjust_diff (nlevs_stats,ifiles,stats_diff)
!
    do n=1,nlevs_stats
      stats_levs(n)=1000.*stats_levs(n) ! units m
    enddo
!
    tab_levs(1)=real(ifiles) 
    do n=2,nlevs_tab
      tab_levs(n)=1000.*(n-2)
    enddo
!
    do i=1,ifiles
      tab_new(1,i)=lats(2,i)
    enddo    
!
    if (trim(c_tab_old) /= 'none') then 
      call read_table (nlevs_tab,max_ifiles,tab_old,c_tab_old)
      ier=0
      do i=1,ifiles
        if (tab_old(1,i) /= tab_new(1,i)) then
          print *,'*** old and new values of lat differ ',tab_old(1,i),tab_new(1,i)
          ier=ier+1
        endif
      enddo
      if (ier > 0) stop    
    else
      tab_old(:,:)=0.
    endif
!
    do n=2,nlevs_tab
      if (stats_levs(1) >= tab_levs(n)) then
        m1=1
        m2=m1
      elseif (stats_levs(nlevs_stats) <= tab_levs(n)) then
        m1=nlevs_stats
        m2=m1
      else  
        do mx=1,nlevs_stats-1
          if (stats_levs(mx) <= tab_levs(n) .and. stats_levs(mx+1) > tab_levs(n)) then
            m1=mx
            m2=mx+1
          endif
        enddo
      endif
!  
      if (m1 /= m2) then
        wd=stats_levs(m2)-stats_levs(m1)
        w1=(stats_levs(m2)-tab_levs(n))/wd
      else
        w1=1.
      endif
      w2=1.-w1
!
      do i=1,ifiles
        var=w1*stats_diff(m1,i)+w2*stats_diff(m2,i)+tab_old(n,i)**2
        tab_new(n,i)=sqrt(max(var,var_min))
      enddo
!
    enddo
!
    call write_table (nlevs_tab,max_ifiles,tab_levs,tab_new,c_tab_new)
!
! test interpolation
    call gpsro_lat_test (nlevs_tab,max_ifiles+1,c_tab_new) 
!
    end program gps_err_var
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine read_table (nlevs,imax,sfac,cfile)
    implicit none
    integer, intent(in) :: nlevs,imax 
    character(len=*), intent(in) :: cfile
    real(4), intent(out) :: sfac(nlevs,imax)
    integer :: i, n, nx
    integer, parameter :: iunit=10
    real(4) :: tablev
    character(len=1) :: c1
!
    open (iunit,file=trim(cfile),form='formatted')
    print *,' '       
    print *,'File opened: read table  file=',trim(cfile)
!
    read (iunit,'(a1)') c1 
!
    print *,'OLD table values'       
    read (iunit,*) nx,tablev,(sfac(1,i),i=1,imax)
    print ('(i3,f9.1,12f9.2)'),nx,tablev,(sfac(n,i),i=1,imax)
    do n=2,nlevs
      read (iunit,*) nx,tablev,(sfac(n,i),i=1,imax)
      print ('(i3,f9.1,12e9.2)'),nx,tablev,(sfac(n,i),i=1,imax)
    enddo 
    close (iunit)
!
    end subroutine read_table 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine write_table (nlevs,imax,tablevs,sfac,cfile)
    implicit none
    integer, intent(in) :: nlevs,imax 
    character(len=*), intent(in) :: cfile
    real(4), intent(in) :: tablevs(nlevs)
    real(4), intent(in) :: sfac(nlevs,imax)
    integer :: i, n, nx
    integer, parameter :: iunit=10
!
    open (iunit,file=trim(cfile),form='formatted')
    print *,' '       
    print *,'File opened: write table  file=',trim(cfile)
!
    n=1
    write (iunit,'(i4)') n
!
    print *,'NEW table values'       
    nx=1
    write (iunit,'(i3,f9.1,12f9.2)') nx,tablevs(nx),(sfac(nx,i),i=1,imax)
    print ('(i3,f9.1,12f9.2)'),nx,tablevs(nx),(sfac(nx,i),i=1,imax)

    do nx=2,nlevs
      write (iunit,'(i3,f9.1,1p12e9.2)') nx,tablevs(nx),(sfac(nx,i),i=1,imax)
      print ('(i3,f9.1,1p12e9.2)'),nx,tablevs(nx),(sfac(nx,i),i=1,imax)
    enddo 
    close (iunit)
!
    end subroutine write_table 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    subroutine read_vcorr_stats (nlevs,statlevs,stats,lats,cfile, &
                                 line_lats,line_gps,ctype)   
    implicit none
    integer, intent(in) :: nlevs
    integer, intent(in) :: line_lats,line_gps 
    real(4), intent(out) :: statlevs(nlevs),lats(3)
    real(4), intent(out) :: stats(nlevs)
    character(len=*), intent(in) :: cfile
    character(len=*), intent(in) :: ctype
!
    integer :: n, nx, ncount, nskip
    integer, parameter :: iunit=10
    real(4) :: s1, s2
    character(len=1) :: c1
!
    open (iunit,file=trim(cfile),form='formatted')
    print *,'File opened: type=',ctype,'  file=',trim(cfile)
!
    nskip=line_lats-1
    do n=1,nskip
      read (iunit,'(a1)') c1 
    enddo
    read (iunit,'(11x,2f9.1)') lats(1),lats(3)
    lats(2)=0.5*(lats(1)+lats(3))
    print *,'lats=',lats(:)
!        
    nskip=line_gps-line_lats
    do n=1,nskip
      read (iunit,'(a1)') c1 
    enddo
!
    do n=1,nlevs
      read (iunit,*) nx,statlevs(n),ncount,s1,s2
      if (n==1) print *,nx,statlevs(n),ncount,s1,s2,ctype
      stats(n)=s2**2
    enddo
!
    close (iunit)
!
    end subroutine read_vcorr_stats 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    subroutine adjust_diff (nlevs,nlats,stats)
    implicit none
    integer :: nlevs,nlats
    real(4) :: stats(nlevs,nlats)
    integer k,j
    real(4) :: smax,sfac
!
    print *,' '
    print *,'DIFF MEAN SQR BEFORE ADJUSTMENT'
    do k=1,nlevs
      print ('(i3,1p12e9.2)'),k,stats(k,:)
    enddo
!
! find maximum value in 10 highest levels
    smax=0.
    do j=1,nlats
      do k=nlevs-9,nlevs
        smax=max(smax,stats(k,j))
      enddo
    enddo        
!
! ensure that all values at uppermost level are within factor of 10 
    smax=0.1*smax
    do j=1,nlats
      stats(nlevs,j)=max(smax,stats(nlevs,j))
    enddo
!
! ensure that values do not increase downward
    sfac=1.
    do j=1,nlats
      do k=nlats-1,1,-1
        stats(k,j)=max(sfac*stats(k+1,j),stats(k,j))
      enddo
    enddo
!
! ensure that values at each level in all columns within factor of 20
    do k=1,nlevs
      smax=0.
      do j=1,nlats
        smax=max(stats(k,j),smax)
      enddo
      smax=0.05*smax
      do j=1,nlats
        stats(k,j)=max(stats(k,j),smax)
      enddo
    enddo    
!
    print *,' '
    print *,'DIFF MEAN SQR AFTER ADJUSTMENT'
    do k=1,nlevs
      print ('(i3,1p12e9.2)'),k,stats(k,:)
    enddo
!
    end subroutine adjust_diff 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_lat_test (nd1,nd2,cfile)
!
   implicit none
   integer :: nd1,nd2
   real(8) :: lat,height,stdv
   logical, parameter :: ltest=.true.
   integer :: i,j,ix,nd2p1
   character(len=*) :: cfile
   real(8) :: err_tab(nd1,nd2)
   print *,' '
   print *,'cfile=',trim(cfile)
   open (10,file=trim(cfile))     
   read (10,'(i4)') ix
   print *,'ix=',ix
   read (10,'(i3,f9.1,12f9.2)') ix,err_tab(1,1),err_tab(1,2:nd2)
   do j=2,nd1
     read (10,'(i3,f9.1,12e9.2)') ix,err_tab(j,1),err_tab(j,2:nd2)
   enddo
   close (10)
!
   lat=88.
   height=0.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=200000.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=10200.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   lat=-88.
   height=0.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=200000.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=10200.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   lat=60.
   height=0.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=200000.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=10200.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   lat=-60.
   height=0.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=200000.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
   height=200.
   call gpsro_lat_func (nd1,nd2,lat,err_tab,height,stdv,ltest)
!
   end subroutine gpsro_lat_test 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine gpsro_lat_func (nerr1,nerr2,lat,errtab,height,stdv,ltest)
!
!  Compute function of height and lat to multiply obs error table value
!  
!QQQQ   use m_kinds, only : rkind2         ! precision 
!
   integer, parameter :: rkind2=8 !QQQQQ
   integer, intent(in) :: nerr1,nerr2
   logical, intent(in) :: ltest
   real(rkind2), intent(in)  :: errtab(nerr1,nerr2)
   real(rkind2), intent(in)  :: lat
   real(rkind2), intent(in)  :: height ! height above geoid
   real(rkind2), intent(out) :: stdv
!
   integer :: ilats,ilat1,ilat2,ilev1,ilev2
   integer :: i, nlat
   real(rkind2) :: vd,v1,v2,hd,h1,h2,s1,s2
!
   ilats=nint(errtab(1,1))
!
   ilat1=2
   do i=3,ilats+1
     if (lat < errtab(1,i)) then
       ilat1=i
     endif
   enddo
   if (lat >= errtab(1,2) .or. lat < errtab(1,ilats+1)) then
     ilat2=ilat1
   else
     ilat2=ilat1+1
   endif
!
   ilev1=2
   do i=3,nerr1   
     if (height > errtab(i,1)) then
       ilev1=i
     endif
   enddo
   if (height <= errtab(2,1) .or. height > errtab(nerr1,1)) then
     ilev2=ilev1
   else
     ilev2=ilev1+1 
   endif
!
   if (ilev1 /= ilev2) then
     vd=errtab(ilev2,1)-errtab(ilev1,1)
     v1=(errtab(ilev2,1)-height)/vd
   else
     vd=0.
     v1=1.
   endif
   v2=1.-v1
!
   if (ilat1 /= ilat2) then
     hd=errtab(1,ilat1)-errtab(1,ilat2)
     h2=(errtab(1,ilat1)-lat)/hd
   else
     hd=0.
     h2=1.
   endif
   h1=1.-h2
!
   s1=h1*errtab(ilev1,ilat1)+h2*errtab(ilev1,ilat2)
   s2=h1*errtab(ilev2,ilat1)+h2*errtab(ilev2,ilat2)
   stdv=v1*s1+v2*s2
!
   if (ltest) then
     print *,' '
     print *,'TEST OF gpsro_lat_func'
     print ('(a,3i5,2f12.2)'),'nerr1,nerr2,ilats,lat,height=', &
            nerr1,nerr2,ilats,lat,height
     print ('(a,4i6)'),'ilev1,ilev2,ilat1,ilat2=',ilev1,ilev2,ilat1,ilat2
     print ('(a,3f10.2,2f7.2)'),'lev1,lev2,vd,v1,v2=', &
            errtab(ilev1,1),errtab(ilev2,1),vd,v1,v2
     print ('(a,3f10.2,2f7.2)'),'lat2,lat1,hd,h,h2=', &
            errtab(1,ilat2),errtab(1,ilat1),hd,h1,h2
     print ('(a,1p4e12.2)'),'tab(ilev1,ilat1:2),tab(ilev2,ilat1:2)=', &
            errtab(ilev1,ilat1),errtab(ilev1,ilat2),errtab(ilev2,ilat1),errtab(ilev2,ilat2)
     print ('(a,1p3e12.3)'),'s1,s2,stdv=',s1,s2,stdv
   endif
!
   end subroutine gpsro_lat_func 
