      program impacts_rad
!
! Compute estimated obs impacts for radiances from GSI adjoint results
! Impacts are summed for one instrument/satellite/time as function of channel
!
      use m_ods_RE
      use m_sat_info_table, only : sat_info_table_read
      use m_sat_info_table, only : sat_info_table_get_1i
!
      implicit none
!
      logical :: lev_OK
      logical :: print_lats
      logical :: lon_range1, lon_range2
!
      integer, parameter :: mlats=2881  ! in veg-type file
      integer, parameter :: mlons=5760  ! in veg-type file
      integer, parameter :: rkinds=8
      integer :: ichan, nchannels
      integer :: k, n
      integer :: n_obs
      integer :: argc
      integer :: ierr
      integer :: n_date, n_time
      integer, allocatable :: itable(:,:)
!
      real(rkinds), parameter :: zero=0._rkinds
      real(rkinds) :: latN, latS, lonW, lonE
      real(rkinds) :: xlon
      real(rkinds) :: ximpact
      real(4) :: vegtype(mlons,mlats)
      real(4) :: obs_lat,obs_lon,veg_obs
      real(rkinds), allocatable :: rtable(:,:)
!
      character(len=4) :: c_surf
      character(len=6) :: c_latN, c_latS, c_lonW, c_lonE
      character(len=8) :: c_date, c_time
      character(len=20)  :: c_sat_inst, c_sat_platform
      character(len=220) :: c_sat_info_file
      character(len=220) :: c_ods_file
      character(len=12), parameter :: file_out='output_stats'
      character(len=*), parameter :: veg_file= &
              '/discover/nobackup/nprive/postproc/geosnr/land.7km.nc'
!
      type (ods_vect) :: ods
      character (len=30) :: type
!       
! Read and check arguments
      argc = iargc()
      if (argc .ne. 11) then
        print *,' usage must be: impacts_rad.x date time',     &
                ' latN latS lonW lonE sfctype sat_info_file',  &
                ' instr platform ods_file' 
        stop
      endif
      call GetArg( 1_4, c_date)
      call GetArg( 2_4, c_time)
      call GetArg( 3_4, c_latN)
      call GetArg( 4_4, c_latS)
      call GetArg( 5_4, c_lonW)
      call GetArg( 6_4, c_lonE)
      call GetArg( 7_4, c_surf)
      call GetArg( 8_4, c_sat_info_file)
      call GetArg( 9_4, c_sat_inst)
      call GetArg(10_4, c_sat_platform)
      call GetArg(11_4, c_ods_file)
!
      read (c_latN, '(f6.1)') latN  
      read (c_latS, '(f6.1)') latS  
      read (c_lonW, '(f6.1)') lonW
      read (c_lonE, '(f6.1)') lonE
      read (c_date, '(i8)')   n_date
      read (c_time, '(i6)')   n_time
!
! The vegtype file is used to specify whether surface is SEA or not
      call read_nc4_2dfield (mlons,mlats,veg_file,'vegtype',  &
                             vegtype,.true.,.true.,ierr) 
!
! The sat info file is used to get the number of channels for a given rad instr.
      call sat_info_table_read (c_sat_info_file,.true.,ierr)
      if (ierr /= 0) stop
      call sat_info_table_get_1i ('instr','nchan',c_sat_inst,nchannels,ierr)
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
      print ('(4(a,1x))'),c_sat_inst,c_sat_platform,c_date,c_time
      print ('(a,i5,2x,4f8.1,2x,a)'),'channels, latN/S, lonW/E c_surf=', &
              nchannels,latN,latS,lonW,lonE,c_surf
      print ('(a)'),trim(c_ods_file)
!
! Check that lats and lons requested are OK
      if (latN <= latS) then
        print *,'ERROR: latN=',latN,' must be > latS=',latS
        stop
      endif
!
! Read ods file for 1 time
      call ods_get_RE (.true.,trim(c_ods_file),n_time,n_obs,ods,ierr)
      if (ierr > 0) then
        print *,'ERROR WHEN READING ODS FILE: ierr=',ierr
        stop
      endif
      if (n_obs == 0) then
        print *,'NO OBS ON THIS FILE'
      else
        print *,'number of obs in file =',n_obs
      endif 
!
! allocate all arrays
      allocate (itable(nchannels+1,2))
      allocate (rtable(nchannels+1,2))
!
! initialize arrays
      itable(:,:)=0
      rtable(:,:)=zero
!
! extract obs actually used in the desired region
      do n=1,n_obs
!
        lev_OK=.false.  ! default value
        ichan=nint(ods%data%lev(n))
        xlon=mod(ods%data%lon(n),360.)
        if (xlon < zero)  xlon=xlon+360.
        if (xlon >= 360.) xlon=zero       ! accounts for round-off of xlon
!
! Determine if obs longitude in desired range
        lon_range2=(lon_range1 .and. xlon >= lonW .and. xlon <= lonE) .or. & 
           ((.not. lon_range1) .and. (xlon >= lonW .or. xlon <= lonE)) 
!
        if ( (abs(ods%data%xvec(n)) > 1.d-30)   .and. &
             (ods%data%qcexcl(n) == 0)          .and. &
             (abs(ods%data%omf(n)) < 1.d8)      .and. &    
             (ods%data%lat(n) >= latS)          .and. &
             (ods%data%lat(n) <= latN)          .and. &
             lon_range2                        ) then
          lev_OK=.true.
        endif
!
        if (lev_OK) then    ! check index
          if (ichan > nchannels .or. ichan < 1 ) then
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
! accumulate counts, impacts, and sum of squared OMF
        if (lev_OK) then
          itable(ichan,1)=itable(ichan,1)+1
          ximpact=ods%data%xvec(n) 
          rtable(ichan,1)=rtable(ichan,1)+ximpact
          rtable(ichan,2)=rtable(ichan,2)+ods%data%omf(n)**2
          if (ximpact < zero) then  ! count obs that decreas J as desired
            itable(ichan,2)=itable(ichan,2)+1 
          endif
        endif
!
      enddo  ! loop over all obs
!
!  Compute sums oer all channels
      do k=1,2
        itable(nchannels+1,k)=sum(itable(1:nchannels,k))
        rtable(nchannels+1,k)=sum(rtable(1:nchannels,k))
      enddo
!
!  Compute rms OMF from sum of squares of OMF
      do k=1,nchannels+1
        if (itable(k,1) > 0) then
          rtable(k,2)=rtable(k,2)/itable(k,1)
          if (rtable(k,2) > zero) then
            rtable(k,2)=sqrt(rtable(k,2))
          endif
        endif
      enddo
!
! write file of computed values
      open (unit=10,file=file_out,form='formatted')
      write (10,'(4(a,1x))') c_sat_inst,c_sat_platform,c_date,c_time
      write (10,'(a,i5,2x,4f8.1,2x,a,i9)') 'channels, latN/S, lonW/E c_surf n_obs=', &
              nchannels,latN,latS,lonW,lonE,c_surf,n_obs
      write (10,'(a)') trim(c_ods_file)
      write (10,'(a)') '    k obsused  impneg     totimpact        rmsomf'
      do k=1,nchannels+1
        write (10,'(i5,2i8,1p2e14.4)') k,itable(k,1:2),rtable(k,1:2)
      enddo
      close (10)     
!
      call ods_clean_RE (ods)
      print *,'Program end'
!
      end program impacts_rad
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
