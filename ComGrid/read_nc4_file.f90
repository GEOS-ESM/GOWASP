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

