   subroutine set_field_file_name (file_field_name,file_name_in, &
                                   file_common_path,c_datetime,  &
                                   file_name_out,ier)
!
!  Construct specific file name from template
!
!  Initial Code:  Ronald Errico  July 15 2014
!
   implicit none
!
   integer, intent(out) :: ier  ! =0 of no errors detected
   character(len=*), intent(in)  :: file_field_name
   character(len=*), intent(in)  :: file_name_in
   character(len=*), intent(in)  :: file_common_path
   character(len=14), intent(in)  :: c_datetime
   character(len=240), intent(out) :: file_name_out
!
   integer :: nlength       ! length of file name read
   integer :: nlength_path  ! length of common file path
   integer :: noccur(20,2)
   integer :: nfound(2)
   integer :: n, n1, n2, j1, j2
!
   character(len=2), parameter :: field_file_char='$#'  
   character(len=20) :: char_copy
   character(len=*), parameter :: my_name='nr_fields_file_1time'
!   
   ier=0
   file_name_out(1:240)=' '
!
! Determine length of file name   
   nlength=len_trim(file_name_in)
   nlength_path=len_trim(file_common_path)
!   
! Find locations of all occurrances of demarcation characters
   noccur(:,:)=0
   nfound(:)=0
   do n2=1,2
     do n=1,nlength
       if (file_name_in(n:n) == field_file_char(n2:n2)) then
         nfound(n2)=nfound(n2)+1
         noccur(nfound(n2),n2)=n
       endif
     enddo
   enddo
   noccur(nfound(1)+1,1)=nlength+1     
!
   if (nfound(1) /= nfound(2)) then
     print *,' '
     print *,'Error in ',my_name
     print *,'Special characters in file name not paired ',nfound(:)
     ier=ier+1
   endif
!
! Construct file name from template 
   if (nfound(1) == 0) then
!
! Vegetation file held in different path than NR files
     if (trim(file_field_name) == 'vegtype' .or. &     
         trim(file_field_name) == 'vegfrac') then
       file_name_out=trim(file_name_in)     
     else
       file_name_out=trim(file_common_path)//trim(file_name_in)
     endif
!
   else
     n1=1
     n2=nlength_path
     file_name_out(n1:n2)=file_common_path(n1:n2)
     n1=n2+1
     n2=n1+noccur(1,1)-2
     if (n2 >= n1) then
       file_name_out(n1:n2)=file_name_in(1:noccur(1,1)-1)
     endif  
!
     do n=1,nfound(1)
       n1=n2+1
       j1=noccur(n,1)+1 
       j2=noccur(n,2)-1
       char_copy(1:20)=' '
       if (file_name_in(j1:j2) == 'ffff') then
         char_copy=trim(file_field_name)
       else if (file_name_in(j1:j2) == 'yyyy') then
         char_copy(1:4)=c_datetime(1:4)
       else if (file_name_in(j1:j2) == 'mm') then
         char_copy(1:2)=c_datetime(5:6)
       else if (file_name_in(j1:j2) == 'dd') then
         char_copy(1:2)=c_datetime(7:8)
       else if (file_name_in(j1:j2) == 'hh') then
         char_copy(1:2)=c_datetime(9:10)
       else if (file_name_in(j1:j2) == 'yy') then
         char_copy(1:2)=c_datetime(3:4)
       else if (file_name_in(j1:j2) == 'hhmm') then
         char_copy(1:2)=c_datetime(9:12)
       else if (file_name_in(j1:j2) == 'yyyymmdd') then
         char_copy(1:8)=c_datetime(1:8)
       else if (file_name_in(j1:j2) == 'yyyymmdd_hh') then
         char_copy(1:11)=c_datetime(1:8)//'_'//c_datetime(9:10)
       else if (file_name_in(j1:j2) == 'yyyymmdd_hhmm') then
         char_copy(1:13)=c_datetime(1:8)//'_'//c_datetime(9:12)
       endif
!
! Replace variable part of file name by specific values
       n2=n1+len_trim(char_copy)-1
       file_name_out(n1:n2)=trim(char_copy)
!
! Copy next part of file name, up to next special character or the end
       n1=n2+1
       n2=n1+noccur(n+1,1)-noccur(n,2)-2
       if (n2 >= n1) then 
         file_name_out(n1:n2)=file_name_in(noccur(n,2)+1:noccur(n+1,1)-1)
       endif 
!       
     enddo
   endif      
!
   end subroutine set_field_file_name
