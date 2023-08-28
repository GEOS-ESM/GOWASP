!
   module m_amv_fields
!
! Contains all routines referencing NR or interpolated fields used by
! the programs compute_amv_params and create_amv
!
   use MAPL_ShmemMod    ! The SHMEM infrastructure   
!
   use m_kinds, only : rkind1
   use m_time_compute, only : rkindh!
!
! Use of the module responsible for some MPI-specific instructions
   use m_die, only : mpi_die
   use m_die, only : die_proc_id
!
   implicit none
   private
!
   public amv_fields_allocate
   public amv_fields_read
   public amv_fields_extract_profiles
   public amv_fields_ipw_grad
   public amv_fields_compute_q_table
   public amv_fields_sea_or_land 
   public amv_fields_shutdown
   private amv_fields_interpolate_1level
   private amv_fields_read_2d
   private amv_fields_read_3d
   private amv_fields_compute_ipw
   private amv_fields_check
!
   logical :: ltest_one_time

!  Global arrays to be allocated using SHMEM
!  ---------------------------------------------
   real(rkind1), pointer :: fld_ps(:,:)     => null() ! surface p (Pa)
   real(rkind1), pointer :: fld_seaf(:,:)   => null() ! frac of sea coverage
   real(rkind1), pointer :: fld_ipw(:,:,:)  => null() ! integrat precip water
   real(rkind1), pointer :: fld_cldf(:,:,:) => null() ! cloud fraction
   real(rkind1), pointer :: fld_qv(:,:,:)   => null() ! water vapor
   real(rkind1), pointer :: fld_u(:,:,:)    => null()    
   real(rkind1), pointer :: fld_v(:,:,:)    => null()
   real(rkind1), pointer :: fld_z(:,:,:)    => null()  ! height of grid point
   real(rkind1), pointer :: flev_read(:,:,:) => null() ! array to read fields
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_allocate (interp_nfields,ier)
!
! Allocate NR and interpolated fields using SHMEM
!
   use m_nr_fields_info, only : field_imax, field_jmax, field_kdim
   use m_kx_table, only : kx_field_imax, kx_field_jmax
   use m_kx_table, only : kx_ipw_nlevs
!
   implicit none
!
   integer, intent(in) :: interp_nfields
   integer, intent(out) :: ier
!
   integer, parameter :: nerr=9
   integer :: n
   integer :: dim1(1)
   integer :: dim2(2)
   integer :: dim3(3)
   integer :: ierr(nerr)          ! returned error flag
!
!  Allocate space for the global fields using SHMEM
   dim2=(/kx_field_imax,kx_field_jmax/)
   call MAPL_AllocNodeArray(fld_ps,  dim2,rc=ierr(1))
   call MAPL_AllocNodeArray(fld_seaf,dim2,rc=ierr(2))
!
   dim3=(/kx_field_imax,kx_field_jmax,kx_ipw_nlevs/)
   call MAPL_AllocNodeArray(fld_ipw, dim3,rc=ierr(3))
!
   dim3=(/kx_field_imax,kx_field_jmax,field_kdim/)
   call MAPL_AllocNodeArray(fld_qv,  dim3,rc=ierr(4))
   call MAPL_AllocNodeArray(fld_cldf,dim3,rc=ierr(5))
   call MAPL_AllocNodeArray(fld_u,   dim3,rc=ierr(6))
   call MAPL_AllocNodeArray(fld_v,   dim3,rc=ierr(7))
   call MAPL_AllocNodeArray(fld_z,   dim3,rc=ierr(8))
!
   dim3=(/field_imax,field_jmax,interp_nfields/)
   call MAPL_AllocNodeArray(flev_read,dim3,rc=ierr(9))
!
   ier=0
   do n=1,nerr
     if (ierr(n) /= 0) then
       ier=ier+1
     endif
   enddo
!
   end subroutine amv_fields_allocate
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_compute_ipw (myid,npet)
!
!  Compute integrated precipitable water between p1 <= p <= p2
!  Units are mm water per m**2 area if p is in Pa.
!
   use m_nr_fields_info, only : field_kdim, field_akbk
   use m_kx_table, only : kx_field_imax, kx_field_jmax
   use m_kx_table, only : kx_ipw_nlevs, kx_ipw_plevs  
   use m_parameters, only : grav
!
   implicit none
!
   integer, intent(in) :: myid,npet
!
   integer :: ix,jx,j1x,kx,n  
   real(rkind1) :: dp,p_edge_above,p_edge_below
   real(rkind1) :: ipw_top,ipw_bot
!
   do n=1,kx_ipw_nlevs     ! consider different sets of levels
     ipw_top=100.*kx_ipw_plevs(1,n)
     ipw_bot=100.*kx_ipw_plevs(2,n)
!
     do jx=1,kx_field_jmax
!
       if (mod(jx,npet) == myid) then
         do ix=1,kx_field_imax
           fld_ipw(ix,jx,n)=0._rkind1
           do kx=1,field_kdim
             p_edge_above=field_akbk(kx,  1)+field_akbk(kx,2)*fld_ps(ix,jx)
             p_edge_below=field_akbk(kx+1,1)+field_akbk(kx+1,2)*fld_ps(ix,jx)
!
! Check if some of current layer overlaps range.
! The first value of dp is the full p-thickness of the layer.
! It is then modified to remove any portion of the layer that does not
! overlap the desired range.
             if (p_edge_above <= ipw_bot .and. p_edge_below >= ipw_top) then
               dp=p_edge_below-p_edge_above
               dp=dp-max(p_edge_below-ipw_bot,0.)-max(ipw_top-p_edge_above,0.)
               fld_ipw(ix,jx,n)=fld_ipw(ix,jx,n)+dp*fld_qv(ix,jx,kx)/grav
             endif
           enddo   ! loop over k
         enddo     ! loop over i
       endif       ! check on processor id
     enddo         ! loop over jx
   enddo           ! loop over n
!
   end subroutine amv_fields_compute_ipw
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_read (myid,npet,nlats,nlons,cdtime1_nr,cdtime2_nr, &
                         interp_nfields,l_interp_time,l_interp_horiz,       &
                         ltest_1time,dlat,dlon,rhours_nr,dhours_nr,rhours_kx)
!
! Read NR fields, interpolate horizontally and temporally if requested, and 
! compute integrated precipitable water.
!
   use m_nr_fields_info, only : field_kdim
!
   implicit none
   include "mpif.h"
!
   logical, intent(in) :: l_interp_time,l_interp_horiz,ltest_1time
   integer, intent(in) :: interp_nfields
   integer, intent(in) :: myid,npet
   integer, intent(in) :: nlats,nlons
   real(rkind1), intent(in) :: dlat,dlon
   real(rkindh), intent(in) :: rhours_nr,dhours_nr,rhours_kx
   character(len=*), intent(in) :: cdtime1_nr,cdtime2_nr 
!
   integer :: k
   integer :: ierr_read  ! returned error flag on each distinct proc
   integer :: ierr
   character(len=*), parameter :: my_name='amv_fields_read'
!
   ltest_one_time=ltest_1time
!
! Read NR data
! Read 1 level at a time using different processors dfor different fields
! and times. Doing this by level reduces the amount of memory required, since 
! NR fields for only 1 level need to be stored simulataneously.
!
   do k=1,field_kdim
     ierr_read=0
!
! Read first required time
     if (myid == 0) then
       call amv_fields_read_3d (     'QV',cdtime1_nr,1,k,myid,ierr_read)
     elseif (myid == 1) then
       call amv_fields_read_3d (  'CLDFR',cdtime1_nr,2,k,myid,ierr_read)
     elseif (myid == 2) then
       call amv_fields_read_3d (      'U',cdtime1_nr,3,k,myid,ierr_read)
     elseif (myid == 3) then
       call amv_fields_read_3d (      'V',cdtime1_nr,4,k,myid,ierr_read)
     elseif (myid == 4) then
       call amv_fields_read_3d (      'Z',cdtime1_nr,5,k,myid,ierr_read)
     elseif (myid == 5 .and. k == 1) then
       call amv_fields_read_2d (     'PS',cdtime1_nr,6,myid,ierr_read)
     elseif (myid == 6 .and. k == 1) then
       call amv_fields_read_2d ('FROCEAN',cdtime1_nr,7,myid,ierr_read)
     endif
!
! Read second required time if interpolation to be performed
     if (l_interp_time) then
       if (myid == 7) then
         call amv_fields_read_3d (   'QV',cdtime2_nr, 8,k,myid,ierr_read)
       elseif (myid == 8) then
         call amv_fields_read_3d ('CLDFR',cdtime2_nr, 9,k,myid,ierr_read)
       elseif (myid == 9) then
         call amv_fields_read_3d (    'U',cdtime2_nr,10,k,myid,ierr_read)
       elseif (myid == 10) then
         call amv_fields_read_3d (    'V',cdtime2_nr,11,k,myid,ierr_read)
       elseif (myid == 11) then
         call amv_fields_read_3d (    'Z',cdtime2_nr,12,k,myid,ierr_read)
       elseif (myid == 12 .and. k == 1) then
         call amv_fields_read_2d (   'PS',cdtime2_nr,13,myid,ierr_read)
       endif
     endif  ! test on l_interp_time
!
     if (ierr_read /= 0) then
       print ('(2(a,i3))'),'Error in a call to read_shmem_data on proc=', &
               myid,'  with error=',ierr_read
       call mpi_die (my_name,70)
     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
! Interpolate all fields horizontally and temporally
     call amv_fields_interpolate_1level (k,myid,npet,nlats,nlons,       &
                           interp_nfields,l_interp_time,l_interp_horiz, &
                           dlat,dlon,rhours_nr,dhours_nr,rhours_kx)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
! end loop over read k 
   enddo
!
! Compute integrated precipitable water over requested p-layers
   call amv_fields_compute_ipw (myid,npet)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                 
!
   end subroutine amv_fields_read
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_ipw_grad (i,j,k,dlon,dlat,lon,lat,prof_ipwg)
!
! Determine range of values in a geographical box centered on the current
! grid point (i,j). The size of the box considered is kx_dx(k) except near
! the poles where it may be smaller. The range is defined as the difference 
! between the max and min values within the box, scaled by the value of 
! value of kx_ipw_bparm(1,kir,k) read from the kx_table file.
!
   use m_kx_table, only : kx_field_imax, kx_field_jmax, kx_ipw_nlevs
   use m_kx_table, only : kx_ipw_bparm, kx_dx, kx_ipw_ids
   use m_parameters, only : earthr,pifac_k1
!
   implicit none
!
   integer, intent(in) :: i,j,k
   real(rkind1), intent(in) :: dlon,dlat
   real(rkind1), intent(in) :: lon,lat
   real(rkind1), intent(out) :: prof_ipwg(kx_ipw_nlevs,2)
!
   integer :: irange,jrange
   integer :: ir1,jr1,kir,ir,jr,jk,ik
   integer :: k_bparm
   real(rkind1) :: wmax,wmin,diff
   real(rkind1) :: dfac
!
! Indicate which set of ipw_bparm values should be used for this obs type k 
   k_bparm=kx_ipw_ids(k)
!
   dfac=earthr*0.001*pifac_k1   ! distance corrresponding to 1 deg lat (km)
   irange=1+kx_dx(k)/(dlon*dfac*cos(lat*pifac_k1))
   irange=min(kx_field_imax/10,irange)
   jrange=1+kx_dx(k)/(dlat*dfac)
!
! SW corner of box to consider 
   ir1=1+mod(kx_field_imax+i-(irange+1)/2,kx_field_imax)
   jr1=max(j-(jrange+1)/2,1)
!
! Loop over ranges for which different ipw are defined.
   do kir=1,kx_ipw_nlevs
     wmax=0.
     wmin=1.e20
!
! Loop over points ir, jr in box 
     do jr=1,jrange
       jk=min(jr1+jr,kx_field_jmax)
       do ir=1,irange
         ik=1+mod(ir1+ir,kx_field_imax)
         wmax=max(fld_ipw(ik,jk,kir),wmax)
         wmin=min(fld_ipw(ik,jk,kir),wmin)
       enddo
     enddo
!
! prof_ipwg(kir,1) is the scaled range for this box with abs(kx_ipw_bparm) 
! the scaling; the sign of kx_ipw_bparm is a flag used elsewhere. 
! prof_ipwg(kir,2) is the unscaled range for this box, used for diagnostic 
! purposes   
     diff=wmax-wmin
     prof_ipwg(kir,1)=min(diff/abs(kx_ipw_bparm(1,kir,k_bparm)),1.)
     prof_ipwg(kir,2)=diff
   enddo
!
   end subroutine amv_fields_ipw_grad
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_interpolate_1level (kid,myid,npet,nlats,nlons, &
                           interp_nfields,l_interp_time,l_interp_horiz, &
                           dlat,dlon,rhours_nr,dhours_nr,rhours_kx)
!
! Interpolate 1 level for all fields, horizontally and temporally as may 
! be requested; otherwise just copy the fields to the appropriate 3-D array 
!
   use m_nr_fields_info, only : field_imax, field_jmax, field_lon_first
!
   implicit none 
!
   logical, intent(in) :: l_interp_time,l_interp_horiz
   integer, intent(in) :: interp_nfields
   integer, intent(in) :: myid,npet
   integer, intent(in) :: kid
   integer, intent(in) :: nlats,nlons
   real(rkind1), intent(in) :: dlat,dlon
   real(rkindh), intent(in) :: rhours_nr,dhours_nr,rhours_kx
!
   integer :: kiz,i,j
   integer :: h_index(2,2) ! grid point indexes used for horiz interp 
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: fh(interp_nfields)
   real(rkind1) :: time1,timei,time2
   real(rkind1) :: lat1,lon1
!
   if (l_interp_time) then
     time1=rhours_nr
     time2=time1+dhours_nr
     timei=rhours_kx
     call get_interp_time_weights (time1,timei,time2,t_weights)
   else
     t_weights(:)=(/1.,0./)
   endif
!
! interpolate (or copy) horizontally
   fh(:)=0.
   do j=1,nlats
     if (mod(j,npet) == myid) then
       lat1=-90.+(j-1)*dlat
       do i=1,nlons
!
! interpolate horizontally
!
         if (l_interp_horiz) then
           lon1=field_lon_first+(i-1)*dlon
           call get_interp_horiz_index (field_imax,field_jmax, &
                              field_lon_first,lat1,lon1,h_index,h_weights)
           do kiz=1,interp_nfields
             fh(kiz)=h_weights(1,1)*flev_read(h_index(1,1),h_index(1,2),kiz) &
                    +h_weights(1,2)*flev_read(h_index(1,1),h_index(2,2),kiz) &
                    +h_weights(2,1)*flev_read(h_index(2,1),h_index(1,2),kiz) &
                    +h_weights(2,2)*flev_read(h_index(2,1),h_index(2,2),kiz)
           enddo
         else
           fh(:)=flev_read(i,j,:)
         endif
!
! copy if 1 time or temporally interpolate if 2 times
         if (l_interp_time ) then
           fld_qv(i,j,kid)=t_weights(1)*fh(1)+t_weights(2)*fh(8)
           fld_cldf(i,j,kid)=t_weights(1)*fh(2)+t_weights(2)*fh(9)
           fld_u(i,j,kid)=t_weights(1)*fh(3)+t_weights(2)*fh(10)
           fld_v(i,j,kid)=t_weights(1)*fh(4)+t_weights(2)*fh(11)
           fld_z(i,j,kid)=t_weights(1)*fh(5)+t_weights(2)*fh(12)
           if (kid == 1) then
             fld_ps(i,j)=t_weights(1)*fh(6)+t_weights(2)*fh(13)
             fld_seaf(i,j)=fh(7)
           endif
         else  ! interpolate in time                                                  
           fld_qv(i,j,kid)=fh(1)
           fld_cldf(i,j,kid)=fh(2)
           fld_u(i,j,kid)=fh(3)
           fld_v(i,j,kid)=fh(4)
           fld_z(i,j,kid)=fh(5)
           if (kid == 1) then
             fld_ps(i,j)=fh(6)
             fld_seaf(i,j)=fh(7)
           endif
         endif
!
       enddo  ! loop over i
     endif  ! test on myid
   enddo  ! loop over j
!
   end subroutine amv_fields_interpolate_1level
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_shutdown()
!
! shmem must deallocate shared memory arrays
!
   implicit none
   integer :: ierr
!
   call MAPL_DeallocNodeArray(fld_ps,   rc=ierr)
   call MAPL_DeallocNodeArray(fld_seaf, rc=ierr)
   call MAPL_DeallocNodeArray(fld_ipw,  rc=ierr)
   call MAPL_DeallocNodeArray(fld_cldf, rc=ierr)
   call MAPL_DeallocNodeArray(fld_qv,   rc=ierr)
   call MAPL_DeallocNodeArray(fld_u,    rc=ierr)
   call MAPL_DeallocNodeArray(fld_v,    rc=ierr)
   call MAPL_DeallocNodeArray(fld_z,    rc=ierr)
   call MAPL_DeallocNodeArray(flev_read,rc=ierr)
!
   end subroutine amv_fields_shutdown
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_read_2d (f_name,cdtime,idfld,myid,iers)
!
!  Read 2d field on netcdf file.
!
   use m_nr_fields_info, only : field_imax, field_jmax
   use m_nr_fields_info, only : field_names, field_types, field_files
   use m_nr_fields_info, only : field_num_2d, field_common_path
   use netcdf           ! for reading the NR files
!
   implicit none
!
   character(len=*), intent(in) :: cdtime
   character(len=*), intent(in) :: f_name                
   integer, intent(in)  :: myid
   integer, intent(in)  :: idfld
   integer, intent(out) :: iers!
   integer :: ier
   integer :: imx, jmx
   integer :: ncid,varid,id
   character(len=120) :: c_notice
   character(len=240) :: file_name
   character(len=*), parameter :: subname='read_shmem_data_2d'
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
     ier=0
     c_notice='Opening file for f='//trim(field_names(id,2,2))// &
              't='//cdtime
     call amv_fields_check (nf90_open(trim(file_name),NF90_NOWRITE,ncid),   &
                            trim(c_notice),ier)
     if (ier /= 0) then
       iers=999
       print ('(2a)'),' ERROR attempting to open file of fields: ',trim(file_name)
     else  ! proceed since file successfully opened
!
! Get dimension information to check
       call amv_fields_check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
       call amv_fields_check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                    'nf90_inq 2',ier)
       call amv_fields_check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
       call amv_fields_check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                    'nf90_inq 4',ier)
!
       if (imx /= field_imax .or. jmx /= field_jmax) then
         print *,'Grid dimension mismatch in routine : read_shmem_data'
         print *,'file_name=',trim(file_name)
         print *,'imax, jmax on file = ',imx,jmx
         print *,'imax, jmax in program = ',field_imax,field_jmax
         iers=iers+10
!
       else  ! read fields if dimensions OK 
         c_notice='Getting vari for f='//trim(field_names(id,2,2))// &
                  ' t='//cdtime
         call amv_fields_check (nf90_inq_varid(ncid,trim(field_names(id,2,2)),varid), &
                    trim(c_notice),ier)
         c_notice='reading field for f='//trim(field_names(id,2,2))// &
                  ' t='//cdtime
         call amv_fields_check (nf90_get_var(ncid,varid,flev_read(:,:,idfld)), &
                   trim(c_notice),ier)
       endif
!
       c_notice='Closing file for f='//trim(field_names(id,2,2))// &
                ' t='//cdtime
       call amv_fields_check (nf90_close(ncid),trim(c_notice),ier)
!
       if (ltest_one_time) then
         print ('(2a,2a16,3i4)'),'Read 2D field: ',         &
                        'f_name,cdtime,idfld,iers,myid =', &
                         f_name,cdtime,idfld,iers,myid 
       endif
!
     endif  ! check on whether file was successfully openend
!
   endif    ! name of file successfully constructed
!
   end subroutine amv_fields_read_2d
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_read_3d (f_name,cdtime,idfld,kid,myid,iers)
!
!  Read part of 3d field on netcdf file.
!
   use m_nr_fields_info, only : field_imax, field_jmax, field_kmax
   use m_nr_fields_info, only : field_lev_first, field_num_3d
   use m_nr_fields_info, only : field_names, field_types, field_files
   use m_nr_fields_info, only : field_num_2d, field_common_path
   use netcdf           ! for reading the NR files
!
   implicit none
!
   character(len=*), intent(in) :: cdtime
   character(len=*), intent(in) :: f_name
   integer, intent(in)  :: myid
   integer, intent(in)  :: idfld
   integer, intent(in)  :: kid                
   integer, intent(out) :: iers
!
   integer :: ier
   integer :: k_get
   integer :: imx, jmx, kmx
   integer :: ncid,varid,id
   integer :: id_start(3)
   integer :: id_count(3)
   character(len=120) :: c_notice
   character(len=240) :: file_name
   character(len=*), parameter :: subname='read_shmem_data_3d'
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
     ier=0
     c_notice='Opening file for f='//trim(field_names(id,2,3))// &
              ' t='//cdtime
     call amv_fields_check (nf90_open(trim(file_name),NF90_NOWRITE,ncid), &
                trim(c_notice),ier)
     if (ier /= 0) then
       iers=999
       print ('(2a)'),' ERROR attempting to open file of fields: ',trim(file_name)
     else  ! proceed since file successfully opened
!
! Get dimension information to check 
       call amv_fields_check (nf90_inq_dimid(ncid,'lon',varid),'nf90_inq 1',ier)
       call amv_fields_check (nf90_inquire_dimension(ncid,varid,c_notice,imx), &
                    'nf90_inq 2',ier)
       call amv_fields_check (nf90_inq_dimid(ncid,'lat',varid),'nf90_inq 3',ier)
       call amv_fields_check (nf90_inquire_dimension(ncid,varid,c_notice,jmx), &
                    'nf90_inq 4',ier)
       call amv_fields_check (nf90_inq_dimid(ncid,'lev',varid),'nf90_inq 5',ier)
       call amv_fields_check (nf90_inquire_dimension(ncid,varid,c_notice,kmx), &
                    'nf90_inq 6',ier)
       if (imx /= field_imax .or. jmx /= field_jmax       &
                             .or. kmx /= field_kmax) then
         print *,'Grid dimension mismatch in routine : read_shmem_data'
         print *,'file_name=',trim(file_name)
         print *,'imax, jmax, kmax on file = ',imx,jmx,kmx
         print *,'imax, jmax in program = ',field_imax,field_jmax,field_kmax
         iers=iers+10
!
       else ! read fields if dimensions OK 
         c_notice='Getting vari for f='//trim(field_names(id,2,3))// &
                  ' t='//cdtime
         call amv_fields_check (nf90_inq_varid(ncid,trim(field_names(id,2,3)),varid), &
                     trim(c_notice),ier)
         c_notice='reading field for f='//trim(field_names(id,2,3))// &
                  ' t='//cdtime
         k_get=field_lev_first+kid-1
         id_start(:)=(/1,1,k_get/)
         id_count(:)=(/field_imax,field_jmax,1/)
         call amv_fields_check (nf90_get_var(ncid,varid,flev_read(:,:,idfld:idfld), &
                     start=id_start,count=id_count),trim(c_notice),ier)
       endif
!
       c_notice='Closing file for f='//trim(f_name)//' t='//cdtime
       call amv_fields_check (nf90_close(ncid),trim(c_notice),ier)
!
       if (ltest_one_time) then
         print ('(2a,2a16,5i4)'),'Read 3D field: ',               &
                  'f_name,cdtime,idfld,kid,k_get,iers,myid  =', &
                   f_name,cdtime,idfld,kid,k_get,iers,myid 
       endif
!
     endif  ! check on whether file was successfully openend
!
   endif    ! name of file successfully constructed
!
   end subroutine amv_fields_read_3d
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_extract_profiles (i,j,l_interp_vert,visir_here,  &
                                           wv_here,prof_p,prof_c,prof_q, &
                                           prof_z,prof_u,prof_v)
!
! Extract required atmos profiles at this location.  Perform vertical 
! interpolation if requested.
!
   use m_nr_fields_info, only : field_akbk, field_akbk_dlev, field_kdim
   use m_kx_table, only : kx_field_kdim, kx_akbk, kx_akbk_dlev
!
   implicit none
!
   logical, intent(in) :: l_interp_vert,visir_here,wv_here
   integer, intent(in) :: i,j
   real(rkind1), intent(out) :: prof_p(kx_field_kdim,3)
   real(rkind1), intent(out) :: prof_c(kx_field_kdim)
   real(rkind1), intent(out) :: prof_q(kx_field_kdim)
   real(rkind1), intent(out) :: prof_z(kx_field_kdim)
   real(rkind1), intent(out) :: prof_u(kx_field_kdim), prof_v(kx_field_kdim)
!
   integer :: kk,kk_nr_test,kk_nr,kk_nrm1
   real(rkind1) :: pdata(field_kdim)
   real(rkind1) :: weight(2)
!
   do kk=1,field_kdim
     pdata(kk)=field_akbk_dlev(kk,1)+ &
               field_akbk_dlev(kk,2)*fld_ps(i,j) ! p at data level 
   enddo
!
   if (l_interp_vert) then
     do kk=1,kx_field_kdim
       prof_p(kk,1)=kx_akbk(kk,1)+ &
                    kx_akbk(kk,2)*fld_ps(i,j)      ! p at interface above 
       prof_p(kk,2)=kx_akbk_dlev(kk,1)+ &
                    kx_akbk_dlev(kk,2)*fld_ps(i,j) ! p at data level 
       prof_p(kk,3)=kx_akbk(kk+1,1)+ &
                    kx_akbk(kk+1,2)*fld_ps(i,j)    ! p at interface below  
!
       do kk_nr_test=2,field_kdim
         if (prof_p(kk,2) < pdata(kk_nr_test)) then
           exit
         endif
         kk_nr=kk_nr_test
       enddo
!
       kk_nrm1=kk_nr-1
       weight(1)=log(pdata(kk_nr)/prof_p(kk,2))/log(pdata(kk_nr)/pdata(kk_nrm1))
       weight(2)=1.-weight(1)
!
       prof_u(kk)=weight(1)*fld_u(i,j,kk_nrm1)+weight(2)*fld_u(i,j,kk_nr)
       prof_v(kk)=weight(1)*fld_v(i,j,kk_nrm1)+weight(2)*fld_v(i,j,kk_nr)
       prof_z(kk)=weight(1)*fld_z(i,j,kk_nrm1)+weight(2)*fld_z(i,j,kk_nr)
       if (visir_here) then
         prof_c(kk)=weight(1)*fld_cldf(i,j,kk_nrm1)+ &
                    weight(2)*fld_cldf(i,j,kk_nr)
       endif
       if (wv_here) then
         prof_q(kk)=weight(1)*fld_qv(i,j,kk_nrm1)+ &
                    weight(2)*fld_qv(i,j,kk_nr)
       endif
     enddo
!
   else  ! no interpolation 
!
     do kk=1,field_kdim
       prof_p(kk,1)=field_akbk(kk,1)+ &
                    field_akbk(kk,2)*fld_ps(i,j)      ! p at interface above
       prof_p(kk,2)=pdata(kk)
       prof_p(kk,3)=field_akbk(kk+1,1)+ &
                    field_akbk(kk+1,2)*fld_ps(i,j)    ! p at interface below 
       prof_u(kk)=fld_u(i,j,kk)
       prof_v(kk)=fld_v(i,j,kk)
       prof_z(kk)=fld_z(i,j,kk)
       if (visir_here) then
         prof_c(kk)=fld_cldf(i,j,kk)
       endif
       if (wv_here) then
         prof_q(kk)=fld_qv(i,j,kk)
       endif
     enddo
!
   endif
!
   end subroutine amv_fields_extract_profiles
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_check (status, loc, ier)
!
   use netcdf           ! for reading the NR files
!
   implicit none
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
   end subroutine amv_fields_check
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_compute_q_table (nlons,nlats)
!
! Compute and print table of max values of q for each lat band and p layer.
! Also for each layer over which an ipw value is computed.
!
   use m_nr_fields_info, only : field_akbk_dlev, field_kdim
   use m_kx_table, only : kx_pbins, kx_jbins_lats, kx_ipw_nlevs 
   use m_kx_table, only : kx_ipw_plevs, kx_ipw_nlevs 
!
   implicit none
!
   integer, intent(in) :: nlons,nlats
!
   integer :: ic,jc,kc,mc
   integer :: nband_latc,npm1c,npc
   real(rkind1) :: dlatc,latc,xbandc
   real(rkind1) :: dpc,pc,pbc  
   real(rkind1), allocatable :: qmax_p_table(:,:)
   real(rkind1), allocatable :: qmax_i_table(:,:)
   real(rkind1), allocatable :: xlatsc(:,:)
!
   allocate (qmax_p_table(kx_pbins,kx_jbins_lats))
   allocate (qmax_i_table(kx_ipw_nlevs,kx_jbins_lats))
   allocate (xlatsc(kx_jbins_lats,2))
   qmax_p_table(:,:)=0.
   qmax_i_table(:,:)=0.
!
   dlatc=180./real(nlats-1)
   npm1c=kx_pbins-1
   dpc=100000./real(kx_pbins)   ! dp in units of Pa
!
   do jc=1,nlats
     latc=-90.+(jc-1)*dlatc
     xbandc=kx_jbins_lats*(latc+90.)/180.
     nband_latc=1+min(int(xbandc),kx_jbins_lats-1)
     do ic=1,nlons
       do kc=1,field_kdim
         pc=field_akbk_dlev(kc,1)+field_akbk_dlev(kc,2)*fld_ps(ic,jc)
         npc=1+min(int(pc/dpc),npm1c)
         qmax_p_table(npc,nband_latc)= &
             max(qmax_p_table(npc,nband_latc),fld_qv(ic,jc,kc))
         do mc=1,kx_ipw_nlevs
           if (pc > 100.*kx_ipw_plevs(1,mc) .and. &
               pc < 100.*kx_ipw_plevs(2,mc)) then
             qmax_i_table(mc,nband_latc)= &
                   max(qmax_i_table(mc,nband_latc),fld_qv(ic,jc,kc))
           endif
         enddo
       enddo
     enddo
   enddo
!
   dlatc=180./real(kx_jbins_lats)
   do jc=1,kx_jbins_lats
     xlatsc(jc,1)=-90.+(jc-1)*dlatc
     xlatsc(jc,2)=xlatsc(jc,1)+dlatc
   enddo
!
  print *,' '
   print *,'qmax_p_table:'
   print ('(22x,20f11.2)'),xlatsc(:,1)
   print ('(22x,20f11.2)'),xlatsc(:,2)
   do npc=1,kx_pbins
     pc=(npc-1)*dpc
     pbc=pc+dpc
     print ('(i2,2f10.1,1p20e11.1)'),npc,pc,pbc,qmax_p_table(npc,:)
   enddo
!
   print *,' '
   print *,'qmax_i_table:'
   print ('(22x,20f11.2)'),xlatsc(:,1)
   print ('(22x,20f11.2)'),xlatsc(:,2)
   do mc=1,kx_ipw_nlevs
     pc=100.*kx_ipw_plevs(1,mc)
     pbc=100.*kx_ipw_plevs(2,mc)
     print ('(i2,2f10.1,1p20e11.1)'),mc,pc,pbc,qmax_i_table(mc,:)
   enddo
!
   deallocate (xlatsc,qmax_p_table,qmax_i_table)
!
   end subroutine amv_fields_compute_q_table
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine amv_fields_sea_or_land (i,j,nband_lors)
!
! Return a flag that indicates whether this point should be cosidered as 
! sea or not.
!
   implicit none
!
   integer, intent(in) :: i,j
   integer, intent(out) :: nband_lors
!
   if (fld_seaf(i,j) > 0.5) then
     nband_lors=1   ! sea point
   else
     nband_lors=0   ! land point
   endif
!
   end subroutine amv_fields_sea_or_land 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x


   end module m_amv_fields
