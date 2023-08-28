   subroutine get_grid_info (imax,jmax,kmax,nfiles,lon_first,file_names,ier)
!
! Get some grid info and compare with that on other files
!
   use netcdf
   implicit none
!
   integer, parameter :: rkind1=4
   integer, parameter :: rkind2=8
!
   integer, intent(in)  :: nfiles
   integer, intent(out) :: imax
   integer, intent(out) :: jmax
   integer, intent(out) :: kmax
   integer, intent(out) :: ier
   real(rkind1), intent(out) :: lon_first
   character(len=*), intent(in) :: file_names(nfiles)
! 
   integer :: n
   integer :: ncid        ! assignd Fortran unit number
   integer :: rc          ! error return code
   integer :: im,jm,km    ! grid dimension information on file 
   integer :: varid       ! id of requested variable in file
!
   real(rkind2), allocatable :: lons(:)
!
   character(len=60) :: c_dum  ! returned value not used
   character(len=240) :: file_name
   character(len=*), parameter :: myname='get_grid_info'
!
   ier=0
   do n=1,nfiles
     file_name=trim(file_names(n))
!
! Open file
     rc=nf90_open(trim(file_name), nf90_nowrite, ncid)
     if ( rc < 0 ) then
       print *,'Error in call to open file in routine=',trim(myname)
       print *,'return code =',rc, ' for n=',n
       print *,'file name to open=',trim(file_name)
       ier=ier+10**n
     else
       print *,'Reading file to get grid info ',trim(file_name)
     endif
!
! Get dimension information
     if (ier < 10**n) then
       call nccheck(nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1')
       call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,im),'nf90_inq 2')
       call nccheck(nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3')
       call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,jm),'nf90_inq 4')
       call nccheck(nf90_inq_dimid(ncid,'lev',varid),'nf90_inq 5')
       call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,km),'nf90_inq 6')
!
! Get lon_first
       allocate (lons(im))
       call nccheck(nf90_inq_varid(ncid,'lon',varid),'nf90_inq 7')
       call nccheck(nf90_get_var(ncid,varid,lons),'nf90_get_var')
!
       print ('(a,3i5,f10.2)'),'im,jm,km,lon_first=',im,jm,km,lons(1)
       if (n == 1) then
         imax=im
         jmax=jm
         kmax=km
         lon_first=real(lons(1),rkind1)
       else
         if (imax /= im) ier=ier+1
         if (jmax /= jm) ier=ier+1
         if (kmax /= km) ier=ier+1
         if (lon_first /= real(lons(1),rkind1)) ier=ier+1
       endif
       deallocate (lons)
!
       call nccheck(nf90_close(ncid))
     endif
   enddo
!
   end subroutine get_grid_info
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
!
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
   subroutine read_nc4_3dfield (imax,jmax,nlevs,file_name,field_name, &
                                field,lstop,lprint,ier)
!
! Read a single 3-d field in netcdf nc4 format
!
   use netcdf
   implicit none
!
   integer, parameter :: rkind1=4
!
   logical, intent(in) :: lprint
   logical, intent(in) :: lstop
   integer, intent(in) :: imax
   integer, intent(in) :: jmax
   integer, intent(in) :: nlevs
   integer, intent(out) :: ier
   real(rkind1), intent(out) :: field(imax,jmax,nlevs)
   character(len=*), intent(in)  :: file_name
   character(len=*), intent(in)  :: field_name
! 
   integer :: ncid        ! assignd Fortran unit number
   integer :: rc          ! error return code
   integer :: im,jm,km    ! grid dimension information on file 
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
   call nccheck(nf90_inq_dimid(ncid,'lev',varid),'nf90_inq 5')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,km),'nf90_inq 6')
   if (im /= imax .or. jm /= jmax .or. km /= nlevs) then 
     print *,'Grid dimension mismatch in routine: ',myname
     print *,'file_name=',trim(file_name)
     print *,'imax, jmax, nlevs on file = ',im,jm,km
     print *,'imax, jmax in program = ',imax,jmax,nlevs
     stop
   endif
!
! Get field
   call nccheck(nf90_inq_varid(ncid,trim(field_name),varid),'nf90_inq 7')
   call nccheck(nf90_get_var(ncid,varid,field),'nf90_get_var')
!
!  Close file
   call nccheck(nf90_close(ncid))
   ier=0 
!
   end subroutine read_nc4_3dfield
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_2dfld (imax,jmax, &
                               cdatetime,field_name,file_name,f)
!
! overwrite single 2d field in existing netcdf file
!
   use netcdf
   implicit none
   integer, intent(in) :: imax,jmax
   real(4), intent(in) :: f(imax,jmax)
   character(len=*), intent(in) :: cdatetime
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: file_name
!
   integer :: itime, idate
   integer :: ncid        ! assignd Fortran unit number
   integer :: im,jm       ! grid dimension information on file 
   integer :: varid       ! id of requested variable in file
   character(len=60) :: c_dum  ! returned value not used
!
   print *,'Attempting to write ',trim(field_name),' to file=',trim(file_name)
!
   call nccheck (nf90_open(trim(file_name), nf90_write, ncid),'open file')
   call nccheck (nf90_redef(ncid))
!
   call nccheck(nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 2d vim')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,im),'nf90_inq 2d im')
   call nccheck(nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 2d vjm')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,jm),'nf90_inq 2d jm')
   if (im /= imax .or. jm /= jmax) then 
     print *,'Grid dimension mismatch in routine: write_nc4_2dfld'
     print *,'file_name=',trim(file_name)
     print *,'imax, jmax on file = ',im,jm
     print *,'imax, jmax in program = ',imax,jmax
     stop
   endif
!
   read (cdatetime(1:8),'(i8)') idate
   read (cdatetime(9:14),'(i6)') itime
   call nccheck (nf90_inq_varid(ncid,'time',varid),'nf90_inq time vid')
   call nccheck (nf90_put_att(ncid,varid,'begin_date',idate),'begin_date')
   call nccheck (nf90_put_att(ncid,varid,'begin_time',itime),'begin_time')
!
   call nccheck (nf90_inq_varid(ncid,trim(field_name),varid),'nf90_inq var 2d')
   call nccheck (nf90_put_var(ncid,varid,f))
!
!   call nccheck (nf90_enddef(ncid),'nf_enddef')
   call nccheck (nf90_close(ncid),'nf90_close')
!
   print *,'2d field written for ',trim(field_name),' to netcdf file'
!
   end subroutine write_nc4_2dfld
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine write_nc4_3dfld (imax,jmax,kmax, &
                               cdatetime,field_name,file_name,f)
!
! overwrite single 3d field in existing netcdf file
!
   use netcdf
   implicit none
   integer, intent(in) :: imax,jmax,kmax
   real(4), intent(inout) :: f(imax,jmax,kmax)
   character(len=*), intent(in) :: cdatetime
   character(len=*), intent(in) :: field_name
   character(len=*), intent(in) :: file_name
!
   integer :: itime, idate
   integer :: ncid        ! assignd Fortran unit number
   integer :: im,jm,km    ! grid dimension information on file 
   integer :: varid       ! id of requested variable in file
   character(len=60) :: c_dum  ! returned value not used
!
   print *,'Attempting to write ',trim(field_name),' to file=',trim(file_name)
!
   call nccheck (nf90_open(trim(file_name), nf90_write, ncid),'nf90_open')
   call nccheck (nf90_redef(ncid),'nf90_redef')
!
   call nccheck(nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 3d vim')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,im),'nf90_inq 3d im')
   call nccheck(nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3d vjm')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,jm),'nf90_inq 3d jm')
   call nccheck(nf90_inq_dimid(ncid,'lev',varid),'nf90_inq 3d vkm')
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,km),'nf90_inq 3d km')
   if (im /= imax .or. jm /= jmax .or. km /= kmax) then 
     print *,'Grid dimension mismatch in routine: write_nc4_3dfld'
     print *,'file_name=',trim(file_name)
     print *,'imax, jmax, nlevs on file = ',im,jm,km
     print *,'imax, jmax, nlevs in program = ',imax,jmax,kmax
     stop
   endif
!
   read (cdatetime(1:8),'(i8)') idate
   read (cdatetime(9:14),'(i6)') itime
   call nccheck (nf90_inq_varid(ncid,'time',varid),'nf90_inq time vid')
   call nccheck (nf90_put_att(ncid,varid,'begin_date',idate))
   call nccheck (nf90_put_att(ncid,varid,'begin_time',itime))
!
   call nccheck (nf90_inq_varid(ncid,trim(field_name),varid),'nf90_inq var 3d')
   call nccheck (nf90_put_var(ncid,varid,f))
!
   call nccheck (nf90_close(ncid),'nf90_close')
!
   print *,'3d field written for ',trim(field_name),' to netcdf file'
!
   end subroutine write_nc4_3dfld
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
