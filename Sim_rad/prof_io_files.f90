   subroutine prof_io_files (ntypes_max,dtype,dir_in,dir_out,cdtime,lprint, &
                             ntypes,input_file_names,output_file_names, &
                             dtypes,ierr)
!
! Construct list of input obs_list and output obs_prof file names from
! date/time and list of obs types. 
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_set_unit_nums, only : un_info
! 
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(in) :: ntypes_max
   integer, intent(out) :: ntypes
   integer, intent(out) :: ierr
   character(len=*), intent(inout) :: dtype
   character(len=*), intent(in) :: dir_in
   character(len=*), intent(in) :: dir_out
   character(len=*), intent(in) :: cdtime
   character(len=*), intent(out) :: input_file_names(ntypes_max)
   character(len=*), intent(out) :: output_file_names(ntypes_max)
   character(len=*), intent(out) :: dtypes(ntypes_max)
!
   integer, parameter :: iunit=un_info
   integer :: n,n1,n2,n3
   integer :: dtype_len,dtype_slashes,dtypes_len
   integer :: dtype_slash(ntypes_max*2)
   integer :: ier
   character(len=*), parameter :: prof_io_tmpl='prof_io_file.tmpl'
   character(len=240) :: input_tmpl
   character(len=240) :: output_tmpl
!
   ierr=0
   open (iunit,file=trim(prof_io_tmpl))
   read (iunit,'(a)') input_tmpl
   read (iunit,'(a)') output_tmpl
   if (lprint) print *,'prof_io_file.tmpl read'
!
   dtypes_len=len(dtypes)
   dtypes(:)(1:dtypes_len)=' '
!
   dtype_len=len_trim(dtype)
   dtype_slashes=1
   dtype_slash(dtype_slashes)=0
   do n=1,dtype_len
     if (dtype(n:n) == '/') then
       dtype_slashes=dtype_slashes+1
       dtype_slash(dtype_slashes)=n
     endif
   enddo
   dtype_slashes=dtype_slashes+1
   dtype_slash(dtype_slashes)=dtype_len+1
   ntypes=dtype_slashes-1
   if (ntypes > ntypes_max) then
     ierr=100 
     print *,'Number of requested data types exceeds maximum allowed'
     print *,'Number requested =',ntypes,'   Max allowed=',ntypes_max
     print *,'This will cause a problem of specification of i/o units for each'
     print *,'Terminating Program'
   endif 
!
   if (ntypes > 1) then
     do n=1,ntypes
       n1=dtype_slash(n)+1
       n2=dtype_slash(n+1)-1
       n3=n2-n1+1
       dtypes(n)(1:n3)=dtype(n1:n2)
     enddo
     dtype(1:dtype_len)=' '
     dtype='ALL'
   else
     dtypes(1)=trim(dtype)
   endif
   if (lprint) print ('(8(a,1x))'),'dtypes=',(trim(dtypes(n)),n=1,ntypes) 
!    
! Apply the routine set_field_file_name here to set required specific 
! input and output obs_list and obs_prof files names from therir generic forms.
   ier=0
   do n=1,min(ntypes,ntypes_max)
     call set_field_file_name (dtypes(n),input_tmpl,dir_in,cdtime, &
            input_file_names(n),ier)
     ierr=ierr+abs(ier) 
     call set_field_file_name (dtypes(n),output_tmpl,dir_out,cdtime, &
            output_file_names(n),ier)
     ierr=ierr+abs(ier)  
   enddo  
!     
   end subroutine prof_io_files
