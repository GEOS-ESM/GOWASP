      program countobs_sum
!
! Augment pre-computed sums of ods data for previous times by adding 
! contributions for a new time 
!
      use m_ods_RE
      use m_sat_info_table, only : sat_info_table_read
      use m_sat_info_table, only : sat_info_table_get_1i
!
      implicit none
!
      logical :: lev_OK
      logical :: print_lats
      logical :: l_conv, l_gps
      logical :: lon_range1, lon_range2
!
      integer, parameter :: mlats=2881
      integer, parameter :: mlons=5760
      integer, parameter :: rkinds=8
      integer, parameter :: ikinds=8
      integer, parameter :: nlook=7200000
      integer, parameter :: nplevs=12
      integer(ikinds), allocatable :: iaccept(:,:)
      integer(ikinds), allocatable :: ibins(:,:,:)
      integer :: ib
      integer :: i
      integer :: i_list(15)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(15)
      integer :: kt_look, kx_look
      integer :: kt_look2, kx_look2
      integer :: iqc
      integer :: ifield, itype
      integer :: k, n, n1, n2
      integer :: nch1, nch2
      integer :: nlat, nlats, nlats2
      integer :: nsave
      integer :: n_bins,  n_bins2
      integer :: n_count, n_count2
      integer :: n_ksmax, n_ksmax2
      integer :: n_date,  n_date2
      integer :: n_time,  n_time2
      integer :: n_obs
      integer :: argc
      integer :: ierr
!
      real(rkinds), parameter :: zero=0._rkinds
      real(rkinds), parameter :: dlat=0.25  ! to determine spatial dist of obs
      real(rkinds) :: clat
      real(rkinds) :: earthr2
      real(rkinds) :: d, d1, d2, d3
      real(rkinds) :: latN, latS, lonW, lonE
      real(rkinds) :: latN2, latS2, lonW2, lonE2
      real(rkinds) :: ods_lev
      real(rkinds) :: x
      real(rkinds) :: xlev, xlon, xlat
      real(rkinds) :: xm1, xm2, xs1, xs2, xa, xb, xc
      real(rkinds) :: x1,    x12
      real(rkinds) :: xbin,  xbin2
      real(rkinds) :: x_values(15)
      real(rkinds) :: xsave(nlook,6)
      real(rkinds) :: dplevs(nplevs)
      real(rkinds) :: z_gps(3)
      real(4) :: vegtype(mlons,mlats)
      real(4) :: obs_lat,obs_lon,veg_obs
      real(rkinds), allocatable :: xsum(:,:,:)
      real(rkinds), allocatable :: plevs(:)
      real(rkinds), allocatable :: xsumsq(:,:,:)
!
      character(len=4) :: c_surf, c_surf2
      character(len=3) :: c_count
      character(len=3) :: c_nbins
      character(len=8) :: c_xmax
      character(len=6) :: c_latN, c_latS, c_lonW, c_lonE
      character(len=3) :: c_kt, c_kx
      character(len=8) :: c_date
      character(len=6) :: c_time
      character(len=20)  :: c_sat_inst
      character(len=220) :: c_sat_info_file
      character(len=220) :: c_file, c_file2
      character(len=11), parameter :: file_in='input_stats'
      character(len=12), parameter :: file_out='output_stats'
      character(len=*), parameter :: veg_file= &
              '/discover/nobackup/nprive/postproc/geosnr/land.7km.nc'
!
      type (ods_vect) :: ods
      character (len=30) :: type
!       
! Read and check arguments
      argc = iargc()
      if (argc .ne. 13) then
        print *,' usage must be: channel_sumods.x ncount date time',     &
                ' latN latS lonW lonE kt kx sfctype filename sat_instr', &
                '  sat_info_file'
        stop
      endif
      call GetArg( 1_4, c_count)
      call GetArg( 2_4, c_date)
      call GetArg( 3_4, c_time)
      call GetArg( 4_4, c_latN)
      call GetArg( 5_4, c_latS)
      call GetArg( 6_4, c_lonW)
      call GetArg( 7_4, c_lonE)
      call GetArg( 8_4, c_kt)
      call GetArg( 9_4, c_kx)
      call GetArg(10_4, c_surf)
      call GetArg(11_4, c_file)
      call GetArg(12_4, c_sat_inst)
      call GetArg(13_4, c_sat_info_file)

!
      read (c_count,'(i3)')   n_count
      read (c_date, '(i8)')   n_date
      read (c_time, '(i6)')   n_time
      read (c_kt,   '(i3)')   kt_look
      read (c_kx,   '(i3)')   kx_look
      read (c_latN, '(f6.1)') latN  
      read (c_latS, '(f6.1)') latS  
      read (c_lonW, '(f6.1)') lonW
      read (c_lonE, '(f6.1)') lonE
!
      call read_nc4_2dfield (mlons,mlats,veg_file,'vegtype',  &
                             vegtype,.true.,.true.,ierr) 
!
      if (kx_look > 99 .and. kx_look < 300) then  ! conventional data
        n_ksmax=nplevs
      else
        call sat_info_table_read (c_sat_info_file,.true.,ierr)
        if (ierr /= 0) stop
        call sat_info_table_get_1i ('instr','nchan',c_sat_inst,n_ksmax,ierr)
        if (ierr /= 0) stop
      endif
      allocate (plevs(n_ksmax)) 
!
! Set plevs for conv data  (assumes nplevs=12) 
      if (kx_look > 99 .and. kx_look < 300) then  ! conventional data
        plevs(1:12)= & 
         (/0.0001, 0.00002, 0.00003, 0.0004, 0.0005, 71.428, 214.286, 357.143, & 
           500., 642.857, 785.714, 928.571 /)
!         (/5., 10., 50., 100., 200., 300., 400., 500., 700., 850., 925., 1000./)
        dplevs(1)=0.
        do n=2,nplevs
          dplevs(n)=0.5*(plevs(n)+plevs(n-1)) 
        enddo
        l_conv=.true.
      else 
        l_conv=.false.
      endif
!
      if (kx_look == 89) then ! GPS impact heights
        z_gps(1)=6.37e6     ! min. radius of curvature of earth considered (m)
        z_gps(2)=6.0e4      ! max. height of obs. considered (m)
        z_gps(3)=z_gps(2)/n_ksmax
        do n=1,n_ksmax
          plevs(n)=(n_ksmax-n+0.5)*z_gps(3)
        enddo
        l_gps=.true.
      else
        l_gps=.false.
      endif
!
      if ((.not. l_conv) .and. (.not. l_gps)) then    
        do n=1,n_ksmax
          plevs(n)=n        ! channel index
        enddo 
      endif
!
! Determine direction of longitude range
      if (lonW < 0.) then
        lonW=lonW+360.
      endif
      if (lonE < 0.) then
        lonE=lonE+360.
      endif
      if (lonW < lonE) then
        lon_range1=.true.
      else
        lon_range1=.false.
      endif
!
      print *,' '
      print ('(2i5,2i10)'),n_count,n_ksmax,n_date,n_time
      print ('(4f8.1,2x,a3)'),latN,latS,lonW,lonE
      print ('(a)'),trim(c_file)
!
! Check that lats requested are OK
      if (latN <= latS) then
        print *,'ERROR: latN=',latN,' must be > latS=',latS
        stop
      endif
!
! Read ods file for 1 time
      call ods_get_RE (.true.,trim(c_file),n_time,n_obs,ods,ierr)
      if (n_obs == 0 .or. ierr /= 0) then
        print *,'Job stopping since either there are no obs on ods file ', &
                '(nobs=',n_obs,') or reading errors detected (ierr=',ierr,')'    
        stop
      else
        print *,'number of obs in file =',n_obs
      endif
!
! allocate all arrays
      nlats=1+int(180./dlat)
      allocate (xsum(n_ksmax,3,2))
      allocate (xsumsq(n_ksmax,3,2))
      allocate (iaccept(n_ksmax,2))
      allocate (ibins(nlats,n_ksmax,2))
!
! initialize arrays
      ibins(:,:,1)=0
      iaccept(:,1)=0
      xsum(:,:,1)=zero
      xsumsq(:,:,1)=zero
!
! save obs in regional subset
      nsave=0  ! counts # observations in desired region, channels, or levels
      do n=1,n_obs
!
! reset longitude so that it is in the range 0-360 
        xlon=mod(ods%data%lon(n),360.)
        if (xlon < zero)  xlon=xlon+360.
        if (xlon >= 360.) xlon=zero       ! accounts for round-off of xlon
!
! check for obs level set too large (nint function below otherwise fails)
! 1.e7 is just some lareg number that still converts properly to an integer
        if (l_gps) then  
          ods_lev=min(1.e7,ods%data%xm(n))   ! impact parameter (height) 
        else
          ods_lev=min(1.e7,ods%data%lev(n))  ! pressure or channel number
        endif
!
! Determine if obs longitude in desired range
        lon_range2=(lon_range1 .and. xlon >= lonW .and. xlon <= lonE) .or. & 
           ((.not. lon_range1) .and. (xlon >= lonW .or. xlon <= lonE)) 
!
        lev_OK=.false.
        if ( (nsave < nlook)                    .and. & 
             (ods%data%qcexcl(n) == 0)          .and. &
             (abs(ods%data%omf(n)) < 1.e8)      .and. &    
             (ods%data%lat(n) > latS)           .and. &
             (ods%data%lat(n) < latN)           .and. &
             lon_range2                        ) then
          lev_OK=.true.
        endif
!
! Check for proper requested data type
        if (l_conv .and. (int(ods%data%kx(n)) /= kx_look .or. &
                          int(ods%data%kt(n)) /= kt_look) ) then
            lev_OK=.false.
        elseif (l_gps .and. int(ods%data%kt(n)) /= kt_look) then 
            lev_OK=.false.
        endif   
!
        if (lev_OK) then    ! assign pressure, height, or channel index 
          if (l_conv) then      ! use pressure level
            do k=1,nplevs
              if (ods_lev > dplevs(k)) then
                xlev=k          ! pressure level bin index
              endif
            enddo
          else if (l_gps) then  ! use impact parameter (height)    
            k=min(n_ksmax-1,int(ods_lev/z_gps(3)))
            xlev=n_ksmax-k      ! height level bin index
          else 
            xlev=ods_lev        ! sat as channel index
          endif
        endif
!
        if (lev_OK) then    ! check index
          if (nint(xlev) > n_ksmax .or. nint(xlev) < 1 ) then
            lev_OK=.false.
          endif
        endif
!
        if (lev_OK) then    ! check land/sea
          obs_lon=xlon
          obs_lat=ods%data%lat(n)
          call get_nearest_field_value (mlons,mlats,obs_lat,obs_lon, &
                                       vegtype,veg_obs)
          if ( (c_surf(1:3) == 'SEA'  .and. nint(veg_obs) /= 17) .or. &
               (c_surf(1:4) == 'LAND' .and. nint(veg_obs) == 17) ) then 
             lev_OK=.false. 
          endif
        endif
!
        if (lev_OK) then
          nsave=nsave+1
          xsave(nsave,1)=ods%data%lat(n)
          xsave(nsave,2)=ods%data%lon(n)
          xsave(nsave,3)=xlev 
          xsave(nsave,4)=ods%data%time(n)
          xsave(nsave,5)=ods%data%omf(n)
          xsave(nsave,6)=ods%data%oma(n)
        endif
!
      enddo  ! loop over all obs
!
      clat=90.0+dlat*0.5
      do n1=1,nsave                   ! loop over all saved obs
        nch1=nint(xsave(n1,3))
        iaccept(nch1,1)=iaccept(nch1,1)+1
        xsum(nch1,1,1)=xsum(nch1,1,1)+xsave(n1,5)               ! omf
        xsum(nch1,2,1)=xsum(nch1,2,1)+xsave(n1,6)               ! oma
        xsum(nch1,3,1)=xsum(nch1,3,1)+xsave(n1,6)-xsave(n1,5)   ! amf
        xsumsq(nch1,1,1)=xsumsq(nch1,1,1)+xsave(n1,5)*xsave(n1,5)
        xsumsq(nch1,2,1)=xsumsq(nch1,2,1)+xsave(n1,6)*xsave(n1,6)
        xsumsq(nch1,3,1)=xsumsq(nch1,3,1)+(xsave(n1,6)-xsave(n1,5))**2
!
        nlat=1+int((xsave(n1,1)+clat)/dlat)
        ibins(nlat,nch1,1)=ibins(nlat,nch1,1)+1      
      enddo      
!
!
! read file pf previously accumulated values
      if (n_count > 1) then
        open (unit=10,file=file_in,form='unformatted')
        read (10) n_count2,n_ksmax2,n_date2,n_time2,nlats2, &
                  latN2,latS2,lonW2,lonE2,kt_look2,kx_look2,c_surf2,c_file2
        read (10) plevs
        read (10) iaccept(:,2),xsum(:,:,2),xsumsq(:,:,2)
        read (10) ibins(:,:,2) 
        close (10)     
        print *,'File read'      
!
! check that parameters on accumulation file to be updated agree with program
        if ( (n_ksmax2 /= n_ksmax)           .or. &
             (latN2 /= latN)                 .or. &
             (latS2 /= latS)                 .or. &
             (lonW2 /= lonW)                 .or. &
             (lonE2 /= lonE)                 .or. &
             (kt_look2 /= kt_look)           .or. &
             (kx_look2 /= kx_look)           .or. &
             (c_surf2 /= c_surf)             .or. &
             (nlats2 /= nlats)             ) then
          print *,' '
          print *,'WRONG FILE ACCESSED X X X X X' 
          print *,'Parameters on file are: '
          print *,'n_ksmax,nlats,kt,kx=',n_ksmax2,nlats2,kt_look2,kx_look2
          print *,'latN,latS,lonW,lonE,c_surf=',latN2,latS2,lonW2,lonE2,c_surf2
          print *,'Parameters specified in current program are: '
          print *,'n_ksmax,nlats,kt,kx=',n_ksmax,nlats,kt_look,kx_look
          print *,'latN,latS,lonW,lonE,c_surf=',latN,latS,lonW,lonE,c_surf
          stop
        endif
!
! update accumulations
        iaccept(:,1)=iaccept(:,1)+iaccept(:,2)
        xsum(:,:,1)=xsum(:,:,1)+xsum(:,:,2)
        xsumsq(:,:,1)=xsumsq(:,:,1)+xsumsq(:,:,2)
        ibins(:,:,1)=ibins(:,:,1)+ibins(:,:,2)
!
      endif
!
! write file of updated accumulated values
      open (unit=10,file=file_out,form='unformatted')
      write (10) n_count,n_ksmax,n_date,n_time,nlats, &
                 latN,latS,lonW,lonE,kt_look,kx_look,c_surf,c_file
      write (10) plevs
      write (10) iaccept(:,1),xsum(:,:,1),xsumsq(:,:,1)
      write (10) ibins(:,:,1) 
      close (10)     
!
      print *,'Program end for time=',n_count,n_date,n_time
!
      end program countobs_sum
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine read_nc4_2dfield (imax,jmax,file_name,field_name, &
                                field,lstop,lprint,ier)
!
! Read a single 2-d field in netcdf nc4 format
!
   use netcdf
   implicit none
!
   integer, parameter :: rkind1=4
   logical, intent(in) :: lprint
   logical, intent(in) :: lstop
   integer, intent(in) :: imax
   integer, intent(in) :: jmax
   integer, intent(out) :: ier
   real(rkind1), intent(out) :: field(imax,jmax)
   character(len=*), intent(in)  :: file_name
   character(len=*), intent(in)  :: field_name
! 
   integer :: ncid        ! assignd Fortran unit number
   integer :: rc          ! error return code
   integer :: im,jm       ! grid dimension information on file 
   integer :: varid       ! id of requested variable in file
   character(len=60) :: c_dum  ! returned value not used
   character(len=*), parameter :: myname='read_nc4_2dfield'
!
! Open file
   rc=nf90_open(trim(file_name), nf90_nowrite, ncid)
   if ( rc < 0 ) then
     print *,'Error in call to GFio_Open from ',myname
     print *,'return code =',rc
     print *,'file name to open=',trim(file_name)
     stop
   elseif (lprint) then
     print *,'Reading file ',trim(file_name)
   endif
!
! Get dimension information
   call nccheck(nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,im),'nf90_inq 2')
   call nccheck(nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,jm),'nf90_inq 4')
   if (im /= imax .or. jm /= jmax) then 
     print *,'Grid dimension mismatch in routine: ',myname
     print *,'file_name=',trim(file_name)
     print *,'imax, jmax on file = ',im,jm
     print *,'imax, jmax in program = ',imax,jmax
     stop
   endif
!
! Get field
   call nccheck(nf90_inq_varid(ncid,trim(field_name),varid),'nf90_inq 5')
   call nccheck(nf90_get_var(ncid,varid,field),'nf90_get_var')
!
!  Close file
   call nccheck(nf90_close(ncid))
   ier=0 
!
   end subroutine read_nc4_2dfield
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine nccheck(status,loc)
!
   use netcdf
   implicit none 
   integer, intent(in) :: status
   character(*), intent(in), optional :: loc
   if (status /= nf90_noerr) then
     if (present(loc)) print *,'Error ',status,' at ',loc
     stop
   end if
!
   end subroutine nccheck
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_nearest_field_value (imax,jmax,obs_lat,obs_lon, &
                                       field,field_value)
!
      implicit none
      integer, parameter :: rkind1=4
      real(rkind1), parameter :: lon_first= -180.  ! first lon for each lat
! input
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      real(rkind1), intent(in) :: obs_lat
      real(rkind1), intent(in) :: obs_lon 
      real(rkind1), intent(in) :: field(imax,jmax)
! output
      real(rkind1), intent(out) :: field_value
! local
      integer :: i,j
      integer :: h_index(2,2)
      real(rkind1) :: h_weights(2,2)
!
      call get_interp_horiz_index (imax,jmax,lon_first,obs_lat,obs_lon, &
                                      h_index,h_weights)
      call get_interp_horiz_nearest (h_weights)
!
      do i=1,2
        do j=1,2
          if (h_weights(i,j) > 0.9) then 
            field_value=field(h_index(i,1),h_index(j,2))
          endif
        enddo
      enddo
!
   end subroutine get_nearest_field_value 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_horiz_index (imax,jmax,lon_first,obs_lat,obs_lon, &
                                      h_index,h_weights)
!
! Compute grid box location that surrounds observation location
!
! Initial Code: Ronald Errico July 15 2014 
!
      implicit none
      integer, parameter :: rkind1=4
! input
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      real(rkind1), intent(in) :: lon_first ! first longitude on each lat circle
      real(rkind1), intent(in) :: obs_lat
      real(rkind1), intent(in) :: obs_lon 
! output
      integer, intent(out) :: h_index(2,2)
      real(rkind1), intent(out) :: h_weights(2,2)
! local
      real(rkind1) :: hwmax
      real(rkind1) :: xlat    ! latitude in units of 180/(jmax-1) rel. to SP
      real(rkind1) :: xlon    ! longitude in units of 180/imax rel. to lon_first
      real(rkind1) :: rel_lon ! longitude relative to lon_first
      real(rkind1) :: lat_weights(2)
      real(rkind1) :: lon_weights(2)
!
      xlat=real(jmax-1,rkind1)*(obs_lat+90._rkind1)/180._rkind1
      h_index(1,2)=min(int(xlat)+1,jmax-1)   ! j value of closest point S
      h_index(2,2)=h_index(1,2)+1            ! j value of closest point N
      lat_weights(1)=real(h_index(2,2)-1,rkind1)-xlat      
      lat_weights(2)=1._rkind1-lat_weights(1)
!      
      rel_lon=mod(obs_lon-lon_first,360._rkind1)
      if (rel_lon < 0._rkind1) rel_lon=rel_lon+360._rkind1
      if (rel_lon >= 360._rkind1) rel_lon=0._rkind1  ! accounts for round off
      xlon=real(imax,rkind1)*rel_lon/360._rkind1
!
      h_index(1,1)=min(int(xlon)+1,imax)      ! i value closest point W
      if (h_index(1,1) == imax) then
        h_index(2,1)=1                        ! i value closest point E
      else
        h_index(2,1)=h_index(1,1)+1           ! i value closest point W
      endif
      lon_weights(2)=xlon-real(h_index(1,1)-1,rkind1)
      lon_weights(1)=1._rkind1-lon_weights(2)
!
      h_weights(1,1)=lon_weights(1)*lat_weights(1)
      h_weights(1,2)=lon_weights(1)*lat_weights(2)
      h_weights(2,1)=lon_weights(2)*lat_weights(1)
      h_weights(2,2)=lon_weights(2)*lat_weights(2)
!
   end subroutine get_interp_horiz_index 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_horiz_nearest (h_weights)
!
! re-define interpolation weights so that the nearest point is given 
! weight 1 and others re-set to 0.
!
      implicit none
      integer, parameter :: rkind1=4
!
      real(rkind1), intent(inout) :: h_weights(2,2)
      integer :: loc_max(2)
!
      loc_max=maxloc(h_weights)
      h_weights(:,:)=0._rkind1
      h_weights(loc_max(1),loc_max(2))=1._rkind1
!
   end subroutine get_interp_horiz_nearest
