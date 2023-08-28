   module m_sat_info_table
!
!  Read sat info table and get requested values for desired dtype/sat/inst
!
   use m_kinds, only : rkind1
!
   implicit none
   private 
   public :: sat_info_table_read 
   public :: sat_info_table_get_ic
   public :: sat_info_table_get_1c
   public :: sat_info_table_get_2c
   public :: sat_info_table_get_1i
   public :: sat_info_table_get_2i
   public :: sat_info_table_get_1r
   public :: sat_info_table_get_2r
   public :: sat_info_table_get_info_4
!
   integer, parameter :: nchar=5  ! num of char variables in table
   integer, parameter :: nintg=4  ! num of intg variables in table
   integer, parameter :: nreal=4  ! num of real variables in table
   integer, parameter :: nmax=200 ! max number of sat/inst pairs in table 
   integer, parameter :: ni0=nchar       ! id prev to 1st intg var in tab
   integer, parameter :: nr0=nchar+nintg ! id prev to 1st real var in tab
   integer, parameter :: nall=nchar+nintg+nreal ! tot num vars in table
   integer :: ntable                            ! num sat/inst pairs in tab
   integer :: itable(nintg,nmax)                ! intg vars in table
   real(rkind1) :: rtable(nreal,nmax)           ! real vars in table
   character(len=20) :: ctable(nchar,nmax)         ! char vars in table
   character(len=20) :: col_names(nall)            ! names of vars in table
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_read (file_sat_info,lprint,ier)
!
! Read the sat_info table file and save values for later reference
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(out) :: ier
   character(len=*), intent(in) :: file_sat_info
!
   integer, parameter :: iunit=10
   integer :: n, ios
!
   ier=0
!   
   open (iunit,file=trim(file_sat_info),status='old',iostat=ios)
   if (ios /= 0) then 
     ier=ios
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open file_sat_info for iunit=', &
                     iunit,' iostat=',ios,' and file name=',trim(file_sat_info)
     ier=1000
     return  
   elseif (lprint) then
     print *,' '
     print ('(3a,i4)'),' Sat_file_info opened as file=',trim(file_sat_info), &
              ' on unit=',iunit
   endif

   read (iunit,'(4(a16,1x),(a20,1x),4(a5,1x),3(a10,1x),a10)') col_names(:) 
   do n=1,nmax   
     read (iunit,'(4(a16,1x),1(a20,1x),4(i5,1x),3(f10.5,1x),f10.5)') &
                 ctable(:,n),itable(:,n),rtable(:,n)
     if (trim(ctable(1,n)) == 'EOF') then
       ntable=n-1
       exit
     endif
   enddo  
   close (iunit)
!
   if (lprint) then 
     print *,' '
     print ('(a,i3)'),'number of sat/inst read: ',ntable
   endif
!
   end subroutine sat_info_table_read 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_ic (clook1,cwant,in1,cout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is character valued based on matching 1 integer value
!
   implicit none
   integer, intent(out) :: ierr           ! return flag (OK =0)
   integer, intent(in)  :: in1            ! value to match
   character(len=*), intent(in) :: clook1 ! name of column to match
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(out) :: cout  ! character string to output
!
   integer :: i,n
   integer :: idin1   ! id of column with character name to match
   integer :: idout   ! 1= first table column with character values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nintg,ni0,clook1,idin1,ierr)
!
! Determine id of column with wanted info
   call sat_info_table_find_name (nchar,0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (in1 == itable(idin1,n) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       cout=trim(ctable(idout,nout))
     endif
   endif
!
   end subroutine sat_info_table_get_ic
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_1c (clook1,cwant,cin1,cout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is character valued based on matching 1 input name
!
   implicit none
   integer, intent(out) :: ierr           ! return flag (OK =0)
   character(len=*), intent(in) :: clook1 ! name of column to match
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to match in 1st column
   character(len=*), intent(out) :: cout  ! character string to output
!
   integer :: n
   integer :: idin1   ! id of column with character name to match
   integer :: idout   ! 1= first table column with character values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of column with wanted info
   call sat_info_table_find_name (nchar,0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       cout=trim(ctable(idout,nout))
     endif
   endif   
!
   end subroutine sat_info_table_get_1c
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_2c (clook1,clook2,cwant,cin1,cin2,cout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is character valued based on matching 2 input names
!
   implicit none
   integer, intent(out) :: ierr   ! return flag (OK =0)
   character(len=*), intent(in) :: clook1 ! name of 1st column to match
   character(len=*), intent(in) :: clook2 ! name of 2nd column to match  
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to find in 1st column
   character(len=*), intent(in) :: cin2   ! name to find in 2nd column
   character(len=*), intent(out) :: cout  ! character string to output
!
   integer :: n
   integer :: idin1, idin2 ! ids of columns with character names to match
   integer :: idout   ! 1= first table column with character values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of 2nd column to search
   call sat_info_table_find_name (nchar,0,clook2,idin2,ierr)
!
! Determine id of column with wanted info
   call sat_info_table_find_name (nchar,0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) .and. &
           trim(cin2) == trim(ctable(idin2,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       cout=trim(ctable(idout,nout))
     endif
   endif   
!
   end subroutine sat_info_table_get_2c
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_1i (clook1,cwant,cin1,iout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is integer valued based on matching 1 input names
!
   implicit none
   integer, intent(out) :: ierr   ! return flag (OK =0)
   integer, intent(out) :: iout   ! integer table value to output
   character(len=*), intent(in) :: clook1 ! name of column to match
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to find in 1st column
!
   integer :: n
   integer :: idin1   ! id of column with character name to match
   integer :: idout   ! 1= first table column with integer values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of integer column with wanted info
   call sat_info_table_find_name (nintg,ni0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       iout=itable(idout,nout)
     endif
   endif   
!
   end subroutine sat_info_table_get_1i
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_2i (clook1,clook2,cwant,cin1,cin2,iout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is integer valued based on matching 2 input names
!
   implicit none
   integer, intent(out) :: ierr   ! return flag (OK =0)
   integer, intent(out) :: iout  ! integer value to output
   character(len=*), intent(in) :: clook1 ! name of 1st column to match
   character(len=*), intent(in) :: clook2 ! name of 2nd column to match  
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to find in 1st column
   character(len=*), intent(in) :: cin2   ! name to find in 2nd column
!
   integer :: n
   integer :: idin1, idin2 ! ids of columns with character names to match
   integer :: idout        ! 1= first table column with integer values
   integer :: nout         ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of 2nd column to search
   call sat_info_table_find_name (nchar,0,clook2,idin2,ierr)
!
! Determine id of integer column with wanted info
   call sat_info_table_find_name (nintg,ni0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) .and. &
           trim(cin2) == trim(ctable(idin2,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       iout=itable(idout,nout)
     endif
   endif   
!
   end subroutine sat_info_table_get_2i
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_1r (clook1,cwant,cin1,rout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is real valued based on matching 1 input names
!
   implicit none
   integer, intent(out) :: ierr           ! return flag (OK =0)
   real(rkind1),intent(out) :: rout       ! real value to putput from table
   character(len=*), intent(in) :: clook1 ! name of column to match 
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to find in 1st column
!
   integer :: n
   integer :: idin1   ! id of column with character name to match
   integer :: idout   ! 1= first table column with real values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of column with wanted info
   call sat_info_table_find_name (nreal,nr0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       rout=rtable(idout,nout)
     endif
   endif   
!
   end subroutine sat_info_table_get_1r
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_2r (clook1,clook2,cwant,cin1,cin2,rout,ierr)
!
! Get information from previously read table of instr/sat information
! when requested info is real valued based on matching 2 input names
!
   implicit none
   integer, intent(out) :: ierr           ! return flag (OK =0)
   real(rkind1),intent(out) :: rout       ! real value to putput from table
   character(len=*), intent(in) :: clook1 ! name of 1st column to match 
   character(len=*), intent(in) :: clook2 ! name of 2nd column to match  
   character(len=*), intent(in) :: cwant  ! name of column to get
   character(len=*), intent(in) :: cin1   ! name to find in 1st column
   character(len=*), intent(in) :: cin2   ! name to find in 2nd column
!
   integer :: n
   integer :: idin1, idin2 ! ids of columns with character names to match
   integer :: idout        ! 1= first table column with real values
   integer :: nout         ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,clook1,idin1,ierr)
!
! Determine id of 2nd column to search
   call sat_info_table_find_name (nchar,0,clook2,idin2,ierr)
!
! Determine id of column with wanted real-valued info
   call sat_info_table_find_name (nreal,nr0,cwant,idout,ierr)
!
! Get wanted info
   if (ierr == 0) then
     nout=0 
     do n=1,ntable   
       if (trim(cin1) == trim(ctable(idin1,n)) .and. &
           trim(cin2) == trim(ctable(idin2,n)) ) then
         nout=n
         exit
       endif
     enddo
     if (nout == 0) then
       ierr=ierr+1
     else
       rout=rtable(idout,nout)
     endif
   endif   
!
   end subroutine sat_info_table_get_2r
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_get_info_4 (nsat_max,ctype,nfound,nchans, &
                                         inst_sat,said,ierr)
!
!  Find names for all satellites that have an instrument of a given type(=dtype)
!  This includes sat names and (WMO) sat-id numbers.
!
   implicit none
   integer, intent(in) :: nsat_max
   integer, intent(out) :: ierr           ! return flag (OK =0)
   integer, intent(out) :: nfound
   integer, intent(out) :: nchans
   integer, intent(out) :: said(nsat_max)
   character(len=*), intent(in) :: ctype ! name of column to match
   character(len=*), intent(out) :: inst_sat(nsat_max)
!
   integer :: i,n
   integer :: idin1   ! id of column with character name to match
   integer :: idout(5)   ! 1= first table column with character values
   integer :: nout    ! row number in table for desired dtype/sat/instr
!
   ierr=0
!
! Determine id of 1st column to search
   call sat_info_table_find_name (nchar,0,'dtype',idin1,ierr)
!
! Determine ids of columns with wanted info
   call sat_info_table_find_name (nchar,0,'instr',idout(1),ierr)
   call sat_info_table_find_name (nchar,0,'sat',idout(2),ierr)
   call sat_info_table_find_name (nintg,ni0,'said',idout(3),ierr)
   call sat_info_table_find_name (nintg,ni0,'nchan',idout(4),ierr)
!
! Get wanted info
   if (ierr == 0) then
     nfound=0
     do n=1,ntable   
       if (trim(ctype) == ctable(idin1,n) ) then
         nfound=nfound+1
         if (nfound > nsat_max) then
           print *,'Too many satellites found: nmax,nfound=',nsat_max,nfound
           stop
         endif
         inst_sat(nfound)=trim(ctable(idout(1),n))//'_'//trim(ctable(idout(2),n))
         said(nfound)=itable(idout(3),n)
         nchans=itable(idout(4),n)  ! assumed the same for all satellites for the same instr.
       endif
     enddo
     if (nfound == 0) then
       ierr=ierr+1
     endif
   endif   
!
   end subroutine sat_info_table_get_info_4
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sat_info_table_find_name (nlook,noffset,clook,id,ierr)
!
!  Find the id of a desired column in the sat info table
!
   implicit none
   integer, intent(in) :: nlook           ! number of columns to ssearch
   integer, intent(in) :: noffset         ! first coulmn to search is noffset+1
   integer, intent(out) :: id             ! id found for wanted column
   integer, intent(inout) :: ierr           
   character(len=*), intent(in) :: clook  ! name of wanted column 
!
   integer :: n
!
   id=0
   do n=1,nlook  
     if (trim(clook) == trim(col_names(n+noffset))) then
       id=n
       exit 
     endif
   enddo
   if (id == 0) then
     print *,'Requested name not found in sat_info_table: ',trim(clook)
     ierr=ierr+1
   endif
!
   end subroutine sat_info_table_find_name 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_sat_info_table
