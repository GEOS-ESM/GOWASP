 program create_gp_obs
!
! Peel off atmospheric profiles from netcdf files of NR fields for subsequent 
! calculations of raobs at selected grid points using a shared memory array 
! using MAPL_ShmemMod
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure
   use netcdf           ! for reading the NR files
!
   use m_kinds, only : rkind1
!
! Use of m_nr_fields_info module that defines info for required fields
   use m_nr_fields_info, only : nr_fields_setup
   use m_nr_fields_info, only : nr_fields_print_info
   use m_nr_fields_info, only : field_imax, field_jmax, field_kmax
   use m_nr_fields_info, only : field_kdim, field_num_3d
   use m_nr_fields_info, only : field_time_slots, field_num_2d
   use m_nr_fields_info, only : field_time_delta, field_time_first
   use m_nr_fields_info, only : field_lon_first, field_lev_first 
   use m_nr_fields_info, only : field_common_path
   use m_nr_fields_info, only : field_names, field_types, field_files
   use m_nr_fields_info, only : field_stride_i, field_stride_j
   use m_nr_fields_info, only : field_akbk_dlev, field_akbk
!
   use m_time_compute, only : time_compute_new_cdtime
!
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
   include "mpif.h"
!
!  Global arrays to be allocated using SHMEM
!  ---------------------------------------------
   real(rkind1), pointer :: fld_ps(:,:)   => null()
   real(rkind1), pointer :: fld_phis(:,:) => null()
   real(rkind1), pointer :: fld_t(:,:,:) => null()
   real(rkind1), pointer :: fld_qv(:,:,:) => null()
   real(rkind1), pointer :: fld_z(:,:,:) => null()
   real(rkind1), pointer :: fld_u(:,:,:) => null()
   real(rkind1), pointer :: fld_v(:,:,:) => null()
   real(rkind1), pointer :: prof_info(:,:) => null()
   real(rkind1), pointer :: prof_data(:,:,:) => null()
!
!  Local variables
!  -------------
   logical :: lprint
   logical :: ltest
!
   integer, parameter :: prof_data_max=1040000  ! max # obs (1/4 deg) 
   integer, parameter :: prof_info_dim1=5
   integer, parameter :: prof_data_dim2=6
   integer, parameter :: obs_ndim2=2   ! # of data fields to copy
   integer, parameter :: obs_ndim3=2   ! # of level specifiers 
   integer, parameter :: obs_ndim4=7   ! # of header info 
   integer, parameter :: bufr_unit_in=20
   integer, parameter :: bufr_unit_out=21
   integer, parameter :: pts_in_msg=40 ! max # of obs reports per message
!
   integer :: nobs  ! obs counter
   integer :: obs_ndim1
   integer :: n3
   integer :: dim2(3)
   integer :: dim3(3)
   integer :: ierr          ! returned error flag
   integer :: ierr_read     ! returned error flag from field reading routines
   integer :: myid          ! processor id number 0, ... ,npet
   integer :: npet          ! number of processors used
   integer :: CoresPerNode
   integer :: obs_list_set  
   integer :: nlist      
   integer :: idatetime
   integer :: ntime         ! time slot index
   integer :: tot_num_obs
   integer     :: argc
   integer(4)  :: iargc
!
   real(rkind1) :: rtime1  ! relative time(hrs) for NR data begin time slot
   real(rkind1) :: obs_time
   real(rkind1) :: obs_info(obs_ndim4)  
   real(rkind1), allocatable :: obs_levels(:,:)  
   real(rkind1), allocatable :: obs_data(:,:)  
!
   character(len=20)  :: dtype        ! data type name
   character(len=1)   :: c_test       ! T or F indicates whether to print more
   character(len=14)  :: cdtime_ref   ! yyyymmddhhmmss of center of assim window
   character(len=14)  :: cdtime1      ! yyyymmddhhmmss of begin of time slot
   character(len=240) :: field_list_file ! .rc file containing field info
   character(len=240) :: bufr_in_file    ! 
   character(len=240) :: bufr_out_file   ! 
   character(len=*), parameter :: my_name='main_program'
!                                       ---
!  Initialize MPI
!  --------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   if (myid == 0) then 
     write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'
     lprint=.true.  ! only print from processor 0   
   else
     lprint=.false.
   endif
   die_proc_id=myid 
!
!  ---------------------------------------------------------
!  Initialize SHMEM
!  ----------------
   CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
   call MAPL_InitializeShmem(rc=ierr)
!
! Read arguments
   argc = iargc()
   if (myid == 0 .and. argc .lt. 6) then
     print *,' usage must be: create_rad_profs.x dtype cdtime_ref', & 
             ' field_list_file bufr_in_file bufr_out_file c_test'
     call mpi_die (my_name,77)
   endif
   call GetArg( 1_4, dtype)
   call GetArg( 2_4, cdtime_ref)
   call GetArg( 3_4, field_list_file)   ! rc file describing fields
   call GetArg( 4_4, bufr_in_file)      ! used to copy bufr table
   call GetArg( 5_4, bufr_out_file)     ! bufr file to be created
   call GetArg( 6_4, c_test)            ! T or F used to print samples
!
   if (myid == 0) then
     print *,'dtype=',trim(dtype)
   endif
!
   if (c_test == 'T') then
     ltest=.true.
   else
     ltest=.false.
   endif
!
!  Get the field requirement info (field names and file templates)
   call nr_fields_setup ('GRIDPT',field_list_file,lprint,ierr)
   if (ierr /= 0) then
     print *,'Error detected in call to gp_fields_setup: ierr=',ierr
     call mpi_die (my_name,ierr)
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!  Allocate space for the global fields using SHMEM
   dim2=(/field_imax,field_jmax,1/)
   call MAPL_AllocNodeArray(fld_ps,  dim2,rc=ierr)
   call MAPL_AllocNodeArray(fld_phis,dim2,rc=ierr)
!
   dim3=(/field_imax,field_jmax,field_kdim/)
   call MAPL_AllocNodeArray(fld_t, dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_qv,dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_z, dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_u, dim3,rc=ierr)
   call MAPL_AllocNodeArray(fld_v, dim3,rc=ierr)
! 
!  Allocate space for globally held prof_data
   dim2=(/prof_info_dim1,prof_data_max,1/) 
   call MAPL_AllocNodeArray(prof_info,dim2,rc=ierr)
   dim3=(/field_kdim,prof_data_dim2,prof_data_max/) 
   call MAPL_AllocNodeArray(prof_data,dim3,rc=ierr)
!
   nobs=0  ! initialize obs location counter
!
!  Loop over time intervals in the data assimilation period
   do ntime=1,field_time_slots      
!   
!  Determine time of field data for this time
!  cdtime1 is the date time in format 'yyyymmddhhmmss'
!  rtime1 is the time relative to the central time of the assimilation period
     rtime1=field_time_first+(ntime-1)*field_time_delta
     call time_compute_new_cdtime (cdtime_ref,cdtime1,rtime1,ierr)
     read (cdtime1(1:10),'(i10)') idatetime    
!
     if (myid == 0 .and. lprint) then
       print ('(a,i2,2x,a)'),'LOOP over time: ',ntime,cdtime1
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
! Read all fields for this time
     ierr_read=0 ! default for unused processors
     if (myid == 0) then
       call read_shem_data_3d ('T',cdtime1,fld_t,ierr_read)
     elseif (myid == 1) then
       call read_shem_data_3d ('QV',cdtime1,fld_qv,ierr_read)
     elseif (myid == 2) then
       call read_shem_data_3d ('Z',cdtime1,fld_z,ierr_read)
     elseif (myid == 3) then
       call read_shem_data_3d ('U',cdtime1,fld_u,ierr_read)
     elseif (myid == 4) then
       call read_shem_data_3d ('V',cdtime1,fld_v,ierr_read)
     elseif (myid == 5) then
       call read_shem_data_2d ('PS',cdtime1,fld_ps,ierr_read)
     elseif (myid == 6) then
       call read_shem_data_2d ('PHIS',cdtime1,fld_phis,ierr_read)
     endif
!
     if (myid < 7 .and. ierr_read /= 0) then
       print ('(2(a,i3))'),'Error in a call to read_shem_data on proc=', &
                myid,'  with error=',ierr_read
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
     if (ierr_read /= 0) then
       call mpi_die (my_name,ierr_read)
     endif
!
! Pull observations from selected points
     call pull_obs_from_points
     call MPI_Barrier(MPI_COMM_WORLD,ierr)           
!
   enddo  ! loop over time
!
!  Write saved obs_info
!
   if (myid == 0 ) then
     obs_ndim1=field_kdim
     allocate (obs_data(obs_ndim1,obs_ndim2))
     allocate (obs_levels(obs_ndim1,obs_ndim3))
!
! Open bufr files (the input one is used to get the bufr table)
     open(unit=bufr_unit_in,file=trim(bufr_in_file),form='unformatted')
     print *,' '
     print ('(3a,i3)'),' bufr input file=',trim(bufr_in_file),     &
            ' opened on unit=',bufr_unit_in
     open(unit=bufr_unit_out,file=trim(bufr_out_file),form='unformatted')
     print ('(3a,i3)'),' bufr output file=',trim(bufr_out_file),   &
            ' opened on unit=',bufr_unit_out
     call openbf(bufr_unit_in,'IN ',bufr_unit_in)
     call openbf(bufr_unit_out,'OUT',bufr_unit_in)  
     call maxout(200000)            ! increase size of ouput bufr
     call datelen (10) ! returns a 10 digit value of idate YYYYMMDDHH
!
     do n3=1,nobs
!
       if (mod(n3,pts_in_msg) == 1) then   ! group reports into messages 
         call openmb (bufr_unit_out,'ADPUPA',idatetime)
       endif          
!
       obs_info(1)=real(nobs)           ! serves as id number
       obs_info(2)=prof_info(1,n3)    ! lon
       obs_info(3)=prof_info(2,n3)    ! lat
       obs_info(4)=prof_info(3,n3)    ! relative time
       obs_info(6)=prof_info(5,n3)    ! station elevation
       obs_info(7)=11.                  ! 'T29' value for raobs 
       obs_levels(:,1)=prof_data(:,3,n3)   ! p levels
       obs_levels(:,2)=prof_data(:,4,n3)   ! z levels
!
       obs_info(5)=120.                   ! kx id for raob mass report 
       obs_data(:,1)=prof_data(:,1,n3)    ! T levels
       obs_data(:,2)=prof_data(:,2,n3)    ! q levels
       if (trim(dtype) /= 'GPWIND') then
         call write_obs_tq (obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4,   &
                            obs_ndim1,obs_ndim2,n3,nobs,bufr_unit_out, &
                            ltest,obs_info,obs_levels,obs_data)
       endif
!
       obs_info(5)=220.                   ! kx id for raob wind report 
       obs_data(:,1)=prof_data(:,5,n3)    ! u levels
       obs_data(:,2)=prof_data(:,6,n3)    ! v levels
       if (trim(dtype) /= 'GPMASS') then
         call write_obs_uv (obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4,   &
                            obs_ndim1,obs_ndim2,n3,nobs,bufr_unit_out, &
                            ltest,obs_info,obs_levels,obs_data)
       endif
!
       if (mod(n3,pts_in_msg) == 0) then  
         call closmg(bufr_unit_out)  
       endif 
!
     enddo  ! loop over obs index
!
     deallocate (obs_levels,obs_data)
     print ('(a,i7,a,i12)'),'number of obs locations=',nobs, &
                            ' number of obs=',nobs*4*field_kdim
!
! Close bufr after final report if not closed already
     if (mod(n3,pts_in_msg) /= 0) then   
        call closmg (bufr_unit_out)  
     endif 
     call closbf (bufr_unit_out)  
     print *,' '
     print *,'Output file written'    
   endif  ! check on processor = ROOT
   call MPI_Barrier(MPI_COMM_WORLD,ierr)           
!
! De-allocate arrays
   call shutdown()
!
!
   contains
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine check(status, loc, ier)
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
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine shutdown()
!
! shmem must deallocate shared memory arrays
!
     call MAPL_DeallocNodeArray(fld_ps,  rc=ierr)
     call MAPL_DeallocNodeArray(fld_phis,rc=ierr)
     call MAPL_DeallocNodeArray(fld_t,   rc=ierr)
     call MAPL_DeallocNodeArray(fld_qv,  rc=ierr)
     call MAPL_DeallocNodeArray(fld_z,   rc=ierr)
     call MAPL_DeallocNodeArray(fld_u,   rc=ierr)
     call MAPL_DeallocNodeArray(fld_v,   rc=ierr)
     call MAPL_DeallocNodeArray(prof_info,rc=ierr)
     call MAPL_DeallocNodeArray(prof_data,rc=ierr)
!
     call MAPL_FinalizeShmem (rc=ierr)
     call MPI_Finalize(ierr)
!
   end subroutine shutdown
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_shem_data_2d (f_name,cdtime,field_in,iers)
!
!  Read part of 2d field on netcdf file. 
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
     character(len=*), parameter :: subname='read_shem_data_2d'
!
     iers=0
     call find_name (field_num_2d,field_names(:,1,2),.false., &
                     subname,f_name,id)
     if (id == 0) then 
       iers=1000
       print *,'FATAL ERROR: Required name not found in ',subname
     else 
       call set_field_file_name (field_names(id,2,2),field_files(id,2),  &
                               field_common_path,cdtime,file_name,ier) 
       iers=iers+ier
     endif 
!
     if (iers == 0) then
       c_notice='Opening file for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid),         &  
                  trim(c_notice),ier)
       if (ier /= 0) then
         iers=999
         print ('(4a)'),' ERROR attempting to open file of ',trim(f_name), &
                      ' field: ',trim(file_name)
         return
       endif
!
! Get dimension information to check
       call check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
       call check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                    'nf90_inq 2',ier)
       call check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
       call check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                    'nf90_inq 4',ier)
       if (imx /= field_imax .or. jmx /= field_jmax) then 
         print *,'Grid dimension mismatch in routine : read_shem_data'
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
!
         c_notice='reading field for f='//trim(field_names(id,2,2))// &
                  ' t='//cdtime
         call check (nf90_get_var(ncid,varid,field_in(:,:)),trim(c_notice),ier)
       endif
!
       c_notice='Closing file for f='//trim(field_names(id,2,2))//' t='//cdtime
       call check (nf90_close(ncid),trim(c_notice),ier)
     endif
!
   end subroutine read_shem_data_2d 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine read_shem_data_3d (f_name,cdtime,field_in,iers)
!
!  Read part of 3d field on netcdf file. 
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
     character(len=*), parameter :: subname='read_shem_data_3d'
!
     iers=0
     call find_name (field_num_3d,field_names(:,1,3),.false., &
                     subname,f_name,id)
!
     if (id == 0) then 
       iers=1000
       print *,'FATAL ERROR: Required name not found in ',subname
     else 
       call set_field_file_name (field_names(id,2,3),field_files(id,3),  &
                                 field_common_path,cdtime,file_name,ier) 
       iers=iers+ier
     endif 
!
     if (iers == 0) then
       c_notice='Opening file for f='//trim(field_names(id,2,3))//' t='//cdtime
       call check (nf90_open(trim(file_name),NF90_NOWRITE,ncid), &  
                  trim(c_notice),ier)
       if (ier /= 0) then
         iers=999
         print ('(4a)'),' ERROR attempting to open file of ',trim(f_name), &
                      ' field: ',trim(file_name)
         return
       endif
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
         print *,'Grid dimension mismatch in routine : read_shem_data'
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
         id_start(:)=(/1,1,field_lev_first/)                                         
         id_count(:)=(/field_imax,field_jmax,field_kdim/)                            
         call check( nf90_get_var(ncid,varid,field_in(:,:,:),start=id_start, &       
                                  count=id_count),trim(c_notice),ier)                
       endif                                                                         
!
       c_notice='Closing file for f='//trim(f_name)//' t='//cdtime
       call check( nf90_close(ncid),trim(c_notice),ier)
!
     endif
!
   end subroutine read_shem_data_3d 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine pull_obs_from_points
!
! Pull observations from selected grid points, equally spaced in lats 
! or lons.
!
   use m_parameters, only   : grav
!
   implicit none
!
   integer :: iy,jy,ky
   real(rkind1) :: dlat,dlon,xlon
!
   dlat=180._rkind1/real(field_jmax)
   dlon=360._rkind1/real(field_imax)
!
   do jy=1,field_jmax,field_stride_j                  ! lat index  
     do iy=1,field_imax,field_stride_i                ! lon index
       if (nobs < prof_data_max) then 
         nobs=nobs+1
         if (mod(nobs-1,npet) == myid) then
           xlon=field_lon_first+(iy-1)*dlon
           if (xlon > 360.) xlon=xlon-360.
           if (xlon < 0.)   xlon=xlon+360.
           prof_info(1,nobs)=xlon                     ! lon of point
           prof_info(2,nobs)=-90._rkind1+(jy-1)*dlat  ! lat of point
           prof_info(3,nobs)=rtime1*0.99999   ! adjusted so that never = +/- 3
           prof_info(4,nobs)=fld_ps(iy,jy)            ! ps
           prof_info(5,nobs)=fld_phis(iy,jy)/grav     ! zs
!
           do ky=1,field_kdim                         ! lev index
             prof_data(ky,1,nobs)=fld_t(iy,jy,ky)
             prof_data(ky,2,nobs)=fld_qv(iy,jy,ky)        
             prof_data(ky,3,nobs)=field_akbk_dlev(ky,1)+ &  
                     field_akbk_dlev(ky,2)*fld_ps(iy,jy)  ! p at data level
             prof_data(ky,4,nobs)=fld_z(iy,jy,ky)        
             prof_data(ky,5,nobs)=fld_u(iy,jy,ky)        
             prof_data(ky,6,nobs)=fld_v(iy,jy,ky)        
           enddo
!
         endif   ! test on myid
       endif     ! test on nobs
!
     enddo       ! loop over iy
   enddo         ! loop over jy
!
   end subroutine pull_obs_from_points
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
 end program create_gp_obs
