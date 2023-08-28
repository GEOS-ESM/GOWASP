   subroutine read_nc4_2dfield (imax,jmax,file_name,field_name, &
                                field,lstop,lprint,ier)
!
! Read a single 2-d field in netcdf nc4 format
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1
   use netcdf
   implicit none
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
   integer :: ier1, ier2, iers
   integer :: ncid        ! assignd Fortran unit number
   integer :: rc          ! error return code
   integer :: im,jm       ! grid dimension information on file 
   integer :: varid       ! id of requested variable in file
   character(len=60) :: c_dum  ! returned value not used
   character(len=*), parameter :: myname='read_nc4_2dfield'
!
! Open file
   rc=nf90_open(trim(file_name), nf90_nowrite, ncid)
   if ( rc /= 0 ) then
     print *,'ERROR in call to GFio_Open from ',myname
     print *,'return code =',rc
     print *,'file name to open=',trim(file_name)
     print *,'check ALL file names requested in rad_thin...rc file'
     stop
   elseif (lprint) then
     print *,'Reading file ',trim(file_name)
   endif
!
! Get dimension information (lon)
   iers=0
   ier1=0
   ier2=0
   call nccheck(nf90_inq_dimid(ncid,'lon',varid),' ',ier1)
   if (ier1 /= 0) then
     call nccheck(nf90_inq_dimid(ncid,'longitude',varid),' ',ier2)
   endif   
   if (ier1 == 0 .and. ier2 == 0) then
     call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,im),'nf90_inq 2',iers)
   else
     print *,'Error: cannot find one of lon or longitude on file'
     print *,'       file=',trim(file_name)
     print *,'       ier1, ier2, iers = ',ier1,ier2,iers 
     iers=iers+1
     im=0
   endif
!
! Get dimension information (lat)
   ier1=0
   ier2=0
   call nccheck(nf90_inq_dimid(ncid,'lat',varid),' ',ier1)
   if (ier1 /= 0) then
     call nccheck(nf90_inq_dimid(ncid,'latitude',varid),' ',ier2)
   endif   
   if (ier1 == 0 .and. ier2 == 0) then
     call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,jm),'nf90_inq 4',iers)
   else
     print *,'Error: cannot find one of lat or latitude on file'
     print *,'       file=',trim(file_name)
     print *,'       ier1, ier2, iers = ',ier1,ier2,iers 
     iers=iers+1
     jm=0
   endif
!
! Compare grid diensions in program and in file
   if (im /= imax .or. jm /= jmax) then 
     print *,'Grid dimension mismatch in routine: ',myname
     print *,'file_name=',trim(file_name)
     print *,'imax, jmax on file = ',im,jm
     print *,'imax, jmax in program = ',imax,jmax
     stop
   endif
!
! Get field
   call nccheck(nf90_inq_varid(ncid,trim(field_name),varid),'nf90_inq 5',iers)
   call nccheck(nf90_get_var(ncid,varid,field),'nf90_get_var',iers)
!
   if (iers /= 0) then
     print *,'Error in reading field=',trim(field_name),' with varid=',varid
     print *,'file_name=',trim(file_name)
     stop
   endif
!
!  Close file
   call nccheck(nf90_close(ncid),iers)
   ier=iers 
!
   end subroutine read_nc4_2dfield
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine nccheck(status,loc,ier)
!
!  Check the return status of a netcdf command
!
   use netcdf
   implicit none 
   integer, intent(in) :: status
   character(*), intent(in) :: loc
   integer, intent(inout) :: ier
   if (status /= nf90_noerr) then
     if (loc /= ' ') then 
       print *,'Error ',status,' at ',loc
     endif
     ier=ier+1
   endif
!
   end subroutine nccheck
!
