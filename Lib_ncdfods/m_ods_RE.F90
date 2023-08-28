   module m_ods_RE
!
! Read data in ods format
!
   use netcdf
   use m_kinds, only : rkind1
!
   implicit none
!
   private 
   public :: ods_get_RE
   public :: ods_clean_RE
   private :: ods_getR_RE
   private :: ods_getI_RE
   public  :: obs_vect
   public  :: ods_vect
!
!  ----------------------------------------------------                             !               
   type obs_vect
      integer :: nobs      ! actual number of observations
      integer :: nvct      ! vector size allocated (nobs .le. nvct)
      integer, allocatable :: kx(:)     ! data source index 
      integer, allocatable :: kt(:)     ! data type   index 
      integer, allocatable :: ks(:)     ! sounding    index 
      integer, allocatable :: time(:)   ! time (relative to the input/output)  
      integer, allocatable :: qcexcl(:) ! On-line QC exclusion flag   
      integer, allocatable :: qchist(:) ! On-line QC history flag
      real(rkind1), allocatable :: lat(:)    ! latitute of obs (degrees) 
      real(rkind1), allocatable :: lon(:)    ! longitude    of obs (degrees)
      real(rkind1), allocatable :: lev(:)    ! level        of obs (hPa)
      real(rkind1), allocatable :: xm(:)     ! atomic metadata (depends on kx) 
      real(rkind1), allocatable :: obs(:)    ! observation value (units depend on kt)
      real(rkind1), allocatable :: omf(:)    ! obs minus forecast (O-F) innovations  
      real(rkind1), allocatable :: oma(:)    ! obs minus analysis (O-A) residuals
      real(rkind1), allocatable :: xvec(:)   ! PSAS CG solution vector                 
    end type obs_vect
!
    type ods_vect
      type(obs_vect) data
    end type ods_vect
!
    contains
! 
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X   
!
    Subroutine ods_get_RE (lprint,file_name,n_time,nd,ods,iers)   
!
   implicit none
!
   logical, intent(in)          :: lprint
   integer, intent(in)          :: n_time  ! time of data set hhmmss
   integer, intent(out)         :: nd      ! number of observations
   integer, intent(out)         :: iers      
   type(ods_vect), intent(out)  ::  ods    
   character(len=*), intent(in) :: file_name
!
   integer :: n
   integer :: ntime_offset
   integer :: ncid        ! assignd Fortran unit number
   integer :: nbatches,batchlen
   integer :: ndays,nsyn,ndim
   integer :: rc
   integer :: varid
   integer, allocatable :: syn_len(:,:)
   character(len=60) :: c_dum  ! returned value not used
! 
   iers=0
   rc=nf90_open(trim(file_name),nf90_nowrite,ncid)
   if (rc == 0) then
     print *,'File opened =',trim(file_name)
   else
     print *,'UNABLE TO OPEN FILE =',trim(file_name)
     iers=999
     return
   endif
!
! Get information describing the netcdf variable syn_len
   call nccheck(nf90_inq_dimid(ncid,'ndays',varid),'netcdf_A1',rc)     
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,ndays),'netcdf_A2',rc)
   call nccheck(nf90_inq_dimid(ncid,'nsyn',varid),'netcdf_A3',rc)
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,nsyn),'netcdf_A4',rc)
!
! Get variable array syn_len to obtain obs count
   call nccheck(nf90_inq_varid(ncid,'syn_len',varid),'netcdf_B1',rc)
   allocate (syn_len(nsyn,ndays))
   call nccheck(nf90_get_var(ncid,varid,syn_len),'netcdf_B2',rc)
   if (rc /= 0) then
     print *,'UNABLE TO DETERMINE NUMBERS OF OBS in file ',trim(file_name)
     iers=888
     return
   endif
   nd=sum(syn_len(1:nsyn,1))
   deallocate (syn_len)
!   
! Get dimensions of arrays decribing batches
   call nccheck(nf90_inq_dimid(ncid,'nbatches',varid),'netcdf_C1',rc)
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,nbatches),'netcdf_C2',rc)
   call nccheck(nf90_inq_dimid(ncid,'batchlen',varid),'netcdf_C3',rc)
   call nccheck(nf90_inquire_dimension(ncid,varid,c_dum,batchlen),'netcdf_C4',rc)
   if (rc /= 0) then
     print *,'UNABLE TO DETERMINE DIMS of VARIABLES in file ',trim(file_name)
     iers=777
     return
   endif
!
! Set obs array dimensions and allocate obs data arrays  
   ndim=max(nd,1)
   ods%data%nobs=nd
   ods%data%nvct=ndim
!
! Allocate obs data arrays  
   allocate (ods%data%kx(ndim))
   allocate (ods%data%kt(ndim))
   allocate (ods%data%ks(ndim))
   allocate (ods%data%time(ndim))
   allocate (ods%data%qcexcl(ndim))
   allocate (ods%data%qchist(ndim))
   allocate (ods%data%lat(ndim))
   allocate (ods%data%lon(ndim))
   allocate (ods%data%lev(ndim))
   allocate (ods%data%obs(ndim))
   allocate (ods%data%omf(ndim))
   allocate (ods%data%oma(ndim))
   allocate (ods%data%xm(ndim))
   allocate (ods%data%xvec(ndim))
!
   if (nd > 0) then
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%kt,'kt',iers)
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%kx,'kx',iers)
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%ks,'ks',iers)
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%time,'time',iers)
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%qcexcl,'qcexcl',iers)
     call ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%qchist,'qchist',iers)
!
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%lat,'lat',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%lon,'lon',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%lev,'lev',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%obs,'obs',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%omf,'omf',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%oma,'oma',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%xm,'xm',iers)
     call ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,ods%data%xvec,'xvec',iers)
!
   endif
   call nccheck(nf90_close(ncid),'close',iers)
!
! Adjust obs time since what is on file is hours with respect to start of 
! current day rather than with respect to central time of current data set.
! The GSI version of m_ods determines this adjustment based on information in 
! the netcdf .ods header.  Here it is hard-wired. 
   if (n_time >= 30000 .and. n_time < 90000) then
     ntime_offset=-360
   elseif (n_time >= 90000 .and. n_time < 150000) then
     ntime_offset=-720
   elseif (n_time >= 150000 .and. n_time < 210000) then
     ntime_offset=-1080
   else 
     ntime_offset=0
   endif
   print ('(a,i6,a)'),'Times on file offset by ',ntime_offset,' min.'
!   
   do n=1,nd
     ods%data%time(n)=ods%data%time(n)+ntime_offset
   enddo
!
   if (iers > 0) then
     print *,'Number of errors detected when attempting read of ods file =',iers
   endif
!
   end subroutine ods_get_RE
! 
!                                                                                           
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X                
!                                  
   subroutine ods_getI_RE (batchlen,nbatches,ncid,nd,ndim,odsI,cname,ier)
   implicit none
   integer, intent(in) :: batchlen,nbatches
   integer, intent(in) :: ncid
   integer, intent(in) :: nd,ndim
   integer, intent(out) :: odsI(ndim)
   integer, intent(inout) :: ier
   character(len=*), intent(in) :: cname
!   
   integer :: varid
   integer :: offset
   integer :: ier1  ! number of unsuccessful nccheck calls
   integer :: i,j,ie,je,k
   integer(1), allocatable :: idata1(:,:)
   integer(2), allocatable :: idata2(:,:)
   real(4)   , allocatable :: rdata(:,:)
!
   ier1=0
   je=nd/batchlen
   ie=nd-je*batchlen
!
   call nccheck(nf90_inq_varid(ncid,cname,varid),'netcdf_I1',ier1)
!
   if (cname /= 'qcexcl') then 
     call nccheck(nf90_get_att(ncid,varid,'add_offset',offset),'netcdf_I2',ier1)
   else
     offset=0
   endif
!
   if (cname == 'kt' .or. cname == 'qcexcl') then 
     allocate (idata1(batchlen,nbatches)) 
     call nccheck(nf90_get_var(ncid,varid,idata1),'netcdf_I3',ier1)
     k=0
     do j=1,je
       do i=1,batchlen
         k=k+1
         odsI(k)=idata1(i,j)+offset
       enddo
     enddo
!
     do i=1,ie
       k=k+1
       odsI(k)=idata1(i,je+1)+offset
     enddo  
     deallocate (idata1)
   endif
!
   if (cname == 'kx' .or. cname == 'qchist'  .or. cname == 'time') then 
     allocate (idata2(batchlen,nbatches))
     call nccheck(nf90_get_var(ncid,varid,idata2),'netcdf_I4',ier1)
     k=0
     do j=1,je
       do i=1,batchlen
         k=k+1
         odsI(k)=idata2(i,j)+offset
       enddo
     enddo
!
     do i=1,ie
       k=k+1
       odsI(k)=idata2(i,je+1)+offset
     enddo  
     deallocate (idata2)
   endif
!            
   if (cname == 'ks') then 
     allocate (rdata(batchlen,nbatches))
     call nccheck(nf90_get_var(ncid,varid,rdata),'netcdf_I5',ier1)
     k=0
     do j=1,je
       do i=1,batchlen
         k=k+1
         odsI(k)=nint(rdata(i,j)+offset)
       enddo
     enddo
!
     do i=1,ie
       k=k+1
       odsI(k)=nint(rdata(i,je+1)+offset)
     enddo  
     deallocate (rdata)
   endif
!
   if (ier1 /= 0) then
     ier=ier+ier1
     print *,'Problem acquiring values for ods data: ',cname,' ier1=',ier1
   endif 
!
   end subroutine ods_getI_RE  
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine ods_getR_RE (batchlen,nbatches,ncid,nd,ndim,odsR,cname,ier)
   implicit none
   integer, intent(in) :: batchlen,nbatches
   integer, intent(in) :: ncid
   integer, intent(in) :: nd,ndim
   integer, intent(inout) :: ier
   real(rkind1), intent(out) :: odsR(ndim)
   character(len=*), intent(in) :: cname
!   
   integer :: varid
   integer :: ier1  ! number of unsuccessful nccheck calls
   integer :: i,j,ie,je,k
   real(4) :: sfac
   real(4), allocatable :: rdata(:,:)
!
   ier1=0
   je=nd/batchlen
   ie=nd-je*batchlen
!
   call nccheck(nf90_inq_varid(ncid,cname,varid),'netcdf_R1',ier1)
   if (cname == 'lat' .or. cname == 'lon') then
     call nccheck(nf90_get_att(ncid,varid,'scale_factor',sfac),'netcdf_R2',ier1)
   else
     sfac=1.
   endif
!
   allocate (rdata(batchlen,nbatches))
   call nccheck(nf90_get_var(ncid,varid,rdata),'netcdf_R3',ier1)
   k=0
   do j=1,je
     do i=1,batchlen
       k=k+1
       odsR(k)=rdata(i,j)*sfac
     enddo
   enddo
!
   do i=1,ie
     k=k+1
     odsR(k)=rdata(i,je+1)*sfac
   enddo  
   deallocate (rdata)
!
   if (ier1 /= 0) then
     ier=ier+ier1
     print *,'Problem acquiring values for ods data: ',cname,' ier1=',ier1
   endif 
!
   end subroutine ods_getR_RE  
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
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine ods_clean_RE (ods)
!
! Deallocate obs data arrays  
!
   type(ods_vect) ::  ods    
!
   deallocate (ods%data%kx)
   deallocate (ods%data%kt)
   deallocate (ods%data%ks)
   deallocate (ods%data%time)
   deallocate (ods%data%qcexcl)
   deallocate (ods%data%qchist)
   deallocate (ods%data%lat)
   deallocate (ods%data%lon)
   deallocate (ods%data%lev)
   deallocate (ods%data%obs)
   deallocate (ods%data%omf)
   deallocate (ods%data%oma)
   deallocate (ods%data%xm)
   deallocate (ods%data%xvec)
!
   end subroutine ods_clean_RE
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_ods_RE
