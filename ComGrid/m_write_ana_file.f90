   module m_write_ana
!
! For use with common_grid_ana.f90
! http://www.unidata.ucar.edu/software/netcdf/examples/programs/
!
   use netcdf
   implicit none
   private
   public :: write_nc4_setup
   public :: write_nc4_2dfld
   public :: write_nc4_3dfld
   public :: write_nc4_close
   integer, parameter :: rkind1=4
   integer :: ncid
   integer :: ps_varid,phis_varid,t_varid,q_varid,u_varid,v_varid
   integer :: delp_varid,ozone_varid,qitot_varid,qltot_varid,ts_varid
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_setup (nlons,nlats,nlevs,file_name)
!
   implicit none
!   
   integer, intent(in) :: nlons,nlats,nlevs
   character(len=*), intent(in) :: file_name
!
   integer :: n
   integer :: lat_dimid,lon_dimid,lev_dimid,time_dimid
   integer :: lat_varid,lon_varid,lev_varid,time_varid
   integer :: dim2d(3),dim3d(4)
   real(rkind1) :: lats(nlats)
   real(rkind1) :: lons(nlons)
   real(rkind1) :: levs(nlevs)
   real(rkind1) :: times(1)
   real(rkind1) :: dlon,dlat
   character(len=*), parameter :: UNITS='units'
   character(len=*), parameter :: lat_name='lat'
   character(len=*), parameter :: lon_name='lon'
   character(len=*), parameter :: lev_name='lev'
   character(len=*), parameter :: time_name='time'
   character(len=*), parameter :: lat_units='degrees_north'
   character(len=*), parameter :: lon_units='degrees_east'
   character(len=*), parameter :: lev_units='eta'
   character(len=*), parameter :: time_units='minutes'
   character(len=*), parameter :: ps_units='hPa'
   character(len=*), parameter :: phis_units='m**2/s**2'
   character(len=*), parameter :: t_units='K'
   character(len=*), parameter :: q_units='g/g'
   character(len=*), parameter :: u_units='m/s'
   character(len=*), parameter :: v_units='m/s'
   character(len=*), parameter :: delp_units='hPa'   
   character(len=*), parameter :: ozone_units='g/g' 
!
   dlat=180./(nlats-1)
   lats(1)=-90.
   lats(nlats)=90.
   do n=2,nlats-1
     lats(n)=-90.+(n-1)*dlat
   enddo
!
   dlon=360./nlons
   do n=1,nlons
     lons(n)=-180.+(n-1)*dlon
   enddo
!
   do n=1,nlevs
     levs(n)=real(n)
   enddo     
!
  times(1)=0.
!
   call check( nf90_create(file_name, nf90_clobber, ncid) )
!
   call check( nf90_def_dim(ncid, lat_name, nlats, lat_dimid) )
   call check( nf90_def_dim(ncid, lon_name, nlons, lon_dimid) )
   call check( nf90_def_dim(ncid, lev_name, nlevs, lev_dimid) )
   call check( nf90_def_dim(ncid, time_name, 1, time_dimid) )

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  call check( nf90_def_var(ncid, lat_name, NF90_REAL, lat_dimid, lat_varid) )
  call check( nf90_def_var(ncid, lon_name, NF90_REAL, lon_dimid, lon_varid) )
  call check( nf90_def_var(ncid, lev_name, NF90_REAL, lev_dimid, lev_varid) )
  call check( nf90_def_var(ncid, time_name, NF90_REAL, time_dimid, time_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, lat_varid, units, lat_units) )
  call check( nf90_put_att(ncid, lon_varid, units, lon_units) )
  call check( nf90_put_att(ncid, lev_varid, units, lev_units) )
  call check( nf90_put_att(ncid, time_varid, units, time_units) )

  !
  dim2d=(/ lon_dimid, lat_dimid, time_dimid /)
  dim3d=(/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.

  call check( nf90_def_var(ncid, 'ps',   NF90_REAL, dim2d, ps_varid) )
  call check( nf90_def_var(ncid, 'phis', NF90_REAL, dim2d, phis_varid) )
  call check( nf90_def_var(ncid, 'ts',   NF90_REAL, dim2d, ts_varid) )
  call check( nf90_def_var(ncid, 'tv',   NF90_REAL, dim3d, t_varid) )
  call check( nf90_def_var(ncid, 'sphu', NF90_REAL, dim3d, q_varid) )
  call check( nf90_def_var(ncid, 'u',    NF90_REAL, dim3d, u_varid) )
  call check( nf90_def_var(ncid, 'v',    NF90_REAL, dim3d, v_varid) )
  call check( nf90_def_var(ncid, 'delp', NF90_REAL, dim3d, delp_varid) )
  call check( nf90_def_var(ncid, 'ozone',NF90_REAL, dim3d, ozone_varid) )
  call check( nf90_def_var(ncid, 'qitot',NF90_REAL, dim3d, qitot_varid) )
  call check( nf90_def_var(ncid, 'qltot',NF90_REAL, dim3d, qltot_varid) )

  ! Assign units attributes to the pressure and temperature netCDF
  ! variables.
  call check( nf90_put_att(ncid, ps_varid,   UNITS, ps_units) )
  call check( nf90_put_att(ncid, phis_varid, UNITS, phis_units) )
  call check( nf90_put_att(ncid, ts_varid,   UNITS, T_units) )
  call check( nf90_put_att(ncid, t_varid,    UNITS, T_units) )
  call check( nf90_put_att(ncid, q_varid,    UNITS, q_units) )
  call check( nf90_put_att(ncid, u_varid,    UNITS, u_units) )
  call check( nf90_put_att(ncid, v_varid,    UNITS, v_units) )
  call check( nf90_put_att(ncid, delp_varid, UNITS, delp_units) )
  call check( nf90_put_att(ncid, ozone_varid,UNITS, ozone_units) )
  call check( nf90_put_att(ncid, qitot_varid,UNITS, q_units) )
  call check( nf90_put_att(ncid, qltot_varid,UNITS, q_units) )

  ! End define mode.
  call check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, lat_varid, lats) )
  call check( nf90_put_var(ncid, lon_varid, lons) )
  call check( nf90_put_var(ncid, lev_varid, levs) )
  call check( nf90_put_var(ncid, time_varid, times) )
!
  print *,'Grid and variable defs written to netcdf file'
!
  end subroutine write_nc4_setup
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_2dfld (nlons,nlats,f,fname)
   implicit none
   integer, intent(in) :: nlons,nlats
   real(rkind1), intent(in) :: f(nlons,nlats)
   character(len=*), intent(in) :: fname
!
   integer :: varid
!
  print *,'Attempting to write ',trim(fname),' to netcdf file'
!
   varid=0
   if (trim(fname) == 'ps' ) then
     varid=ps_varid 
   elseif (trim(fname) == 'ts' ) then
     varid=ts_varid 
   else
     varid=phis_varid
   endif 
!
   call check( nf90_put_var(ncid, varid, f) )

   end subroutine write_nc4_2dfld
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_3dfld (nlons,nlats,nlevs,f,fname)
   implicit none
   integer, intent(in) :: nlons,nlats,nlevs
   real(rkind1), intent(in) :: f(nlons,nlats,nlevs)
   character(len=*), intent(in) :: fname
!
   integer :: varid
!
  print *,'Attempting to write ',trim(fname),' to netcdf file'
!
   varid=0
   if (trim(fname) == 'tv' ) then
     varid=t_varid 
   else if (trim(fname) == 'sphu' ) then
     varid=q_varid
   else if (trim(fname) == 'u' ) then
     varid=u_varid
   else if (trim(fname) == 'v' ) then
     varid=v_varid
   else if (trim(fname) == 'delp' ) then
     varid=delp_varid
   else if (trim(fname) == 'ozone' ) then
     varid=ozone_varid
   else if (trim(fname) == 'qitot' ) then
     varid=qitot_varid
   else if (trim(fname) == 'qltot' ) then
     varid=qltot_varid
   endif 
!
   call check( nf90_put_var(ncid, varid, f) )

   end subroutine write_nc4_3dfld
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_close
   implicit none
   call check( nf90_close(ncid) )
   end subroutine write_nc4_close
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check  
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_write_ana
