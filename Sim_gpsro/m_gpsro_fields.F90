   module m_gpsro_fields
!
! Module for reading and using NR fields required by create_gpsro
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
   use netcdf           ! for reading the NR files
!
   use m_kinds, only : rkind1
   use m_nr_fields_info, only : field_imax, field_jmax, field_kmax
   use m_nr_fields_info, only : field_kdim
!
   implicit none
   include "mpif.h"
!
   private
   public :: gpsro_fields_setup
   public :: gpsro_fields_read
   public :: gpsro_fields_clean
   public :: gpsro_fields_tv2t
   public :: horiz_interp_2dfld 
   public :: horiz_interp_3dfld 
!
!  Global arrays to be allocated using SHMEM
!  ---------------------------------------------
   real(rkind1), pointer :: fld_phis(:,:)  => null()  ! phis
   real(rkind1), pointer :: fld_pst1(:,:)  => null()  ! ps at t1
   real(rkind1), pointer :: fld_pst2(:,:)  => null()  ! ps at t2
   real(rkind1), pointer :: fld_1t1(:,:,:) => null()  ! t at t1
   real(rkind1), pointer :: fld_1t2(:,:,:) => null()  ! t at t2
   real(rkind1), pointer :: fld_2t1(:,:,:) => null()  ! q at t1
   real(rkind1), pointer :: fld_2t2(:,:,:) => null()  ! q at t2
! 
   character(len=*), parameter :: my_name='m_shmem_fields'
!                                       ---
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_fields_setup 
!
!  Allocate space for the global fields using SHMEM
!
   implicit none
!   
   integer :: dim2(2)
   integer :: dim3(3)
   integer :: ierr          ! returned error flag
!
   dim2=(/field_imax,field_jmax/)
   call MAPL_AllocNodeArray(fld_phis,dim2,rc=ierr)
   call MAPL_AllocNodeArray(fld_pst1,dim2,rc=ierr)
   call MAPL_AllocNodeArray(fld_pst2,dim2,rc=ierr)
!
   dim3=(/field_imax,field_jmax,field_kdim/)
   call MAPL_AllocNodeArray(fld_1t1,dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_1t2,dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_2t1,dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_2t2,dim3,rc=ierr)
!
   end subroutine gpsro_fields_setup    
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_fields_read (myid,cdtime1,cdtime2,cset,ierr_read)
!
! Call Shared-memory reader for requested set of NR fields
! Read 1 field per processor
!
   implicit none
!
   integer, intent(in)  :: myid              ! processor id 
   integer, intent(out) :: ierr_read         ! error flag 
   character(len=*), intent(in) :: cdtime1   ! 1st file datetime 
   character(len=*), intent(in) :: cdtime2   ! 2nd file datetime 
   character(len=*), intent(in) :: cset      ! set of fields to read
!
   character(len=*), parameter :: my_sub=my_name//'::fields_read'
!
   ierr_read=0
!
   if (cset == 'PHIS') then  ! read PHIS field 
     if (myid == 0) then 
       call read_data_2d ('PHIS',cdtime1,fld_phis,ierr_read)
     endif 
   endif
!
   if (cset == 'PTQ') then  ! Read PS, T and QV  
     if (myid == 0) then
       call read_data_2d ('PS',cdtime1,fld_pst1,ierr_read)
     elseif (myid == 1) then
       call read_data_2d ('PS',cdtime2,fld_pst2,ierr_read)
     elseif (myid == 2) then
       call read_data_3d ('T',cdtime1,fld_1t1,ierr_read)
     elseif (myid == 3) then
       call read_data_3d ('T',cdtime2,fld_1t2,ierr_read)
     elseif (myid == 4) then
       call read_data_3d ('QV',cdtime1,fld_2t1,ierr_read)
     elseif (myid == 5) then
       call read_data_3d ('QV',cdtime2,fld_2t2,ierr_read)
     endif
   endif
!
   if (ierr_read /= 0) then
     print ('(a,i2,3a,i4)'),'Error in read_data on proc=', &
              myid,' for set=',cset,' with error=',ierr_read
   endif
!
   end subroutine gpsro_fields_read
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine check (status, loc, ier)
!
   integer, intent(in) :: status
   integer, intent(inout) :: ier
   character(len=*), intent(in) :: loc
!
   if (status /= NF90_NOERR) then
     print *,'Error at ', loc
     print *,NF90_STRERROR(status)
     ier=ier+1
   endif
!
   end subroutine check
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_fields_clean (ierr)
!
! shmem must deallocate shared memory arrays
!
   implicit none
   integer, intent(out) :: ierr 
!
   call MAPL_DeallocNodeArray(fld_phis,rc=ierr)
   call MAPL_DeallocNodeArray(fld_pst1,rc=ierr)
   call MAPL_DeallocNodeArray(fld_pst2,rc=ierr)
   call MAPL_DeallocNodeArray(fld_1t1, rc=ierr)
   call MAPL_DeallocNodeArray(fld_1t2, rc=ierr)
   call MAPL_DeallocNodeArray(fld_2t1, rc=ierr)
   call MAPL_DeallocNodeArray(fld_2t2, rc=ierr)
!
   end subroutine gpsro_fields_clean
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_fields_tv2t (iflag)
!
! Replace virtual temperature by temperature if former was read field
! Returned value of iflag = -1 if field not found, =0 if no conversion 
! necessary, =1 if conversion performed
!
   use m_parameters, only : ratio4R
   use m_nr_fields_info, only : field_num_3d, field_names
!
   implicit none
   integer, intent(out) :: iflag
!
   integer :: id
   character(len=*), parameter :: my_sub=my_name//'::gpsro_fields_tv2t'
!
   iflag=-1
   call find_name (field_num_3d,field_names(:,1,3),.false.,my_sub,'T',id)
   if (id > 0) then
     iflag=0
     if (trim(field_names(id,2,3)) == 'TV' .or. &
         trim(field_names(id,2,3)) == 'tv' .or. &
         trim(field_names(id,2,3)) == 'Tv') then ! change field values to T
       iflag=1 
       fld_1t1(:,:,:)=fld_1t1(:,:,:)/(1.+ratio4R*fld_2t1(:,:,:))
       fld_1t2(:,:,:)=fld_1t2(:,:,:)/(1.+ratio4R*fld_2t2(:,:,:))
     endif
   endif
!
   end subroutine gpsro_fields_tv2t 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_data_2d (f_name,cdtime,field_in,iers)
!
!  Read part of 2d field on netcdf file. 
!
     use m_nr_fields_info, only : field_num_2d, field_names 
     use m_nr_fields_info, only : field_common_path, field_files
!
     implicit none
!
     character(len=*), intent(in) :: cdtime
     character(len=*), intent(in) :: f_name                
     integer, intent(out) :: iers
     real(rkind1), intent(out) :: field_in(field_imax,field_jmax)
!                         
     integer :: ier
     integer :: imx, jmx
     integer :: ncid,varid,id
     integer :: id_start(3)
     integer :: id_count(3)
     character(len=120) :: c_notice
     character(len=240) :: file_name
     character(len=*), parameter :: my_sub=my_name//'::read_data_2d'
!
     iers=0
     call find_name (field_num_2d,field_names(:,1,2),.false., &
                     my_sub,f_name,id)
     if (id == 0) then 
       iers=1000
       print *,'FATAL ERROR: Required name not found in ',my_sub
     else 
       call set_field_file_name (field_names(id,2,2),field_files(id,2),  &
                                 field_common_path,cdtime,file_name,ier) 
       iers=iers+ier
     endif 
!
     if (iers == 0) then
       ier=0
       c_notice='Opening file for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid),         &  
                  trim(c_notice),ier)
       if (ier /= 0) then
          iers=999
          print ('(2a)'),' ERROR attempting to open file of fields: ',trim(file_name)
       else  ! proceed since file successfully opened
!
! Get dimension information to check
         call check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
         call check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                     'nf90_inq 2',ier)
         call check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
         call check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                     'nf90_inq 4',ier)
         if (imx /= field_imax .or. jmx /= field_jmax) then 
           print *,'Grid dimension mismatch in routine = ',my_sub
           print *,'file_name=',trim(file_name)
           print *,'imax, jmax on file = ',imx,jmx
           print *,'imax, jmax in program = ',field_imax,field_jmax
           iers=iers+10
!
         else  ! read fields if dimensions OK
           c_notice='Getting vari for f='//trim(field_names(id,2,2))// &
                    ' t='//cdtime
           call check (nf90_inq_varid(ncid,trim(field_names(id,2,2)),varid), &
                       trim(c_notice),ier)
           c_notice='reading field for f='//trim(field_names(id,2,2))// &
                    ' t='//cdtime
           call check (nf90_get_var(ncid,varid,field_in(:,:)),trim(c_notice),ier) 
         endif
!
         c_notice='Closing file for f='//trim(field_names(id,2,2))//' t='//cdtime
         call check (nf90_close(ncid),trim(c_notice),ier)
       endif  ! check on whether file was successfully openend
!
     endif    ! name of file successfully constructed
!
   end subroutine read_data_2d 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_data_3d (f_name,cdtime,field_in,iers)
!
!  Read part of 3d field on netcdf file. 
!
     use m_nr_fields_info, only : field_num_3d, field_names
     use m_nr_fields_info, only : field_common_path, field_files

!
     implicit none
!
     character(len=*), intent(in) :: cdtime
     character(len=*), intent(in) :: f_name                
     integer, intent(out) :: iers
     real(rkind1), intent(out) :: field_in(field_imax,field_jmax,field_kdim)
!                         
     integer :: ier
     integer :: imx, jmx, kmx
     integer :: ncid,varid,id
     integer :: id_start(3)
     integer :: id_count(3)
     character(len=120) :: c_notice
     character(len=240) :: file_name
     character(len=*), parameter :: my_sub=my_name//'::read_data_3d'
!
     iers=0
     call find_name (field_num_3d,field_names(:,1,3),.false., &
                     my_sub,f_name,id)
!
     if (id == 0) then 
       iers=1000
       print *,'FATAL ERROR: Required name not found in ',my_sub
     else 
       call set_field_file_name (field_names(id,2,3),field_files(id,3),  &
                                 field_common_path,cdtime,file_name,ier) 
       iers=iers+ier
     endif 
!
     if (iers == 0) then
       ier=0
       c_notice='Opening file for f='//trim(field_names(id,2,3))//' t='//cdtime
       call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid), &  
                  trim(c_notice),ier)
       if (ier /= 0) then
          iers=999
          print ('(4a)'),' ERROR attempting to open file of ',trim(f_name), &
                         ' field: ',trim(file_name)
       else  ! proceed since file successfully opened
!
! Get dimension information to check
         call check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
         call check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                     'nf90_inq 2',ier)
         call check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
         call check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                     'nf90_inq 4',ier)
         call check (nf90_inq_dimid(ncid,'lev',varid),'nf90_inq 5',ier)
         call check (nf90_inquire_dimension(ncid,varid,c_notice,kmx), &
                     'nf90_inq 6',ier)
         if (imx /= field_imax .or. jmx /= field_jmax       &
                               .or. kmx /= field_kmax) then 
           print *,'Grid dimension mismatch in routine = ',my_sub
           print *,'file_name=',trim(file_name)
           print *,'imax, jmax, kmax on file = ',imx,jmx,kmx
           print *,'imax, jmax in program = ',field_imax,field_jmax,field_kmax
           iers=iers+10
!
         else ! read fields if dimensions OK
           c_notice='Getting vari for f='//trim(field_names(id,2,3))// &
                    ' t='//cdtime
           call check (nf90_inq_varid(ncid,trim(field_names(id,2,3)),varid), &
                       trim(c_notice),ier)
           c_notice='reading field for f='//trim(field_names(id,2,3))// &
                    ' t='//cdtime
           call check( nf90_get_var(ncid,varid,field_in(:,:,:)),trim(c_notice),ier)
         endif
!
         c_notice='Closing file for f='//trim(f_name)//' t='//cdtime
         call check( nf90_close(ncid),trim(c_notice),ier)
       endif  ! check on whether file was successfully openend
!
     endif    ! name of file successfully constructed
!
   end subroutine read_data_3d 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine horiz_interp_2dfld (h_index,hw,cfld,ntimes,fh)
!
! Horizontally Interpolate requested 2d fields at 1 or 2 times
!
   implicit none
!
   integer, intent(in) :: ntimes
   integer, intent(in) :: h_index(2,2)
   real(rkind1), intent(in) :: hw(2,2)        
   real(rkind1), intent(out) :: fh(1,ntimes)
   character(len=*), intent(in) :: cfld
!
! Interpolate horizontally at 2 times and then interpolate temporally
! (Except for PHIS that is at 1 time only)
   if (cfld == 'PHIS') then
     fh(1,1)=hw(1,1)*fld_phis(h_index(1,1),h_index(1,2)) &  
            +hw(1,2)*fld_phis(h_index(1,1),h_index(2,2)) &
            +hw(2,1)*fld_phis(h_index(2,1),h_index(1,2)) &
            +hw(2,2)*fld_phis(h_index(2,1),h_index(2,2)) 
   elseif (cfld == 'PS') then 
     fh(1,1)=hw(1,1)*fld_pst1(h_index(1,1),h_index(1,2)) &
            +hw(1,2)*fld_pst1(h_index(1,1),h_index(2,2)) &
            +hw(2,1)*fld_pst1(h_index(2,1),h_index(1,2)) &
            +hw(2,2)*fld_pst1(h_index(2,1),h_index(2,2)) 
   endif
!
   if (ntimes == 2) then 
     if (cfld == 'PS') then 
       fh(1,2)=hw(1,1)*fld_pst2(h_index(1,1),h_index(1,2)) &
              +hw(1,2)*fld_pst2(h_index(1,1),h_index(2,2)) &
              +hw(2,1)*fld_pst2(h_index(2,1),h_index(1,2)) &
              +hw(2,2)*fld_pst2(h_index(2,1),h_index(2,2)) 
     endif
   endif
!
   end subroutine horiz_interp_2dfld
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine horiz_interp_3dfld (nfh,h_index,hw,cfld,ntimes,k1,k2,fh)
!
! Horizontally Interpolate requested 3d fields at 1 or 2 times
!
   implicit none
!
   integer, intent(in) :: nfh
   integer, intent(in) :: ntimes
   integer, intent(in) :: k1,k2
   integer, intent(in) :: h_index(2,2)
   real(rkind1), intent(in) :: hw(2,2)        
   real(rkind1), intent(out) :: fh(nfh,ntimes)
   character(len=*), intent(in) :: cfld
!
   integer :: k,k2x
!
   k2x=min(k2,field_kdim)  
!
   if (cfld == 'T') then
     do k=k1,k2x
       fh(k,1)=hw(1,1)*fld_1t1(h_index(1,1),h_index(1,2),k) &
              +hw(1,2)*fld_1t1(h_index(1,1),h_index(2,2),k) &
              +hw(2,1)*fld_1t1(h_index(2,1),h_index(1,2),k) &
              +hw(2,2)*fld_1t1(h_index(2,1),h_index(2,2),k) 
     enddo
   elseif (cfld == 'Q') then     
     do k=k1,k2x
       fh(k,1)=hw(1,1)*fld_2t1(h_index(1,1),h_index(1,2),k) &
              +hw(1,2)*fld_2t1(h_index(1,1),h_index(2,2),k) &
              +hw(2,1)*fld_2t1(h_index(2,1),h_index(1,2),k) &
              +hw(2,2)*fld_2t1(h_index(2,1),h_index(2,2),k) 
     enddo
   endif
!
   if (ntimes == 2) then
     if (cfld == 'T') then
       do k=k1,k2x
         fh(k,2)=hw(1,1)*fld_1t2(h_index(1,1),h_index(1,2),k) &
                +hw(1,2)*fld_1t2(h_index(1,1),h_index(2,2),k) &
                +hw(2,1)*fld_1t2(h_index(2,1),h_index(1,2),k) &
                +hw(2,2)*fld_1t2(h_index(2,1),h_index(2,2),k) 
       enddo
     elseif (cfld == 'Q') then     
       do k=k1,k2x
         fh(k,2)=hw(1,1)*fld_2t2(h_index(1,1),h_index(1,2),k) &
                +hw(1,2)*fld_2t2(h_index(1,1),h_index(2,2),k) &
                +hw(2,1)*fld_2t2(h_index(2,1),h_index(1,2),k) &
                +hw(2,2)*fld_2t2(h_index(2,1),h_index(2,2),k) 
       enddo
     endif
   endif

!
   end subroutine horiz_interp_3dfld 
!
!
   end module m_gpsro_fields
