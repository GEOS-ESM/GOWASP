   subroutine find_name (nlist,flist,lstop,myname,fname,id)
!
!  Find first index of occurance of fname in character array flist
!  Return 0 if not found.  If not found and lstop=.true., then stop
!  exection before returning.
!
!  If myname(1)=' ', then no error info is printed if fname not found
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   implicit none
   logical, intent(in) :: lstop
   integer, intent(in) :: nlist   
   character(len=*), intent(in) :: myname
   character(len=*), intent(in) :: flist(nlist)
   character(len=*), intent(in) :: fname
   integer, intent(out) :: id
!
   integer :: n
!
   id=0
   do n=1,nlist
     if (id == 0 .and. trim(fname) == trim(flist(n))) then
       id=n
     endif
   enddo
!
   if (id == 0 .and. myname(1:1) /= ' ') then
     print *,'Error in find_name called from ',myname
     print *,'name to find = ',trim(fname)
     print *,'list to search ='     
     print ('(10(1x,a10))'),flist(:)
     if (lstop) then
       stop
     endif
   endif
!
   end subroutine find_name 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine find_name_2 (nlist,flist,lstop,myname,fname,id)
!
!  Find first index of occurance of fname in the initial character string 
!  (i.e., before a blank) within a character array flist.
!  Return 0 if not found.  If not found and lstop=.true., then stop
!  exection before returning.
!
! Initial Code: Ronald Errico  July 15 2014
!
   implicit none
   logical, intent(in) :: lstop
   integer, intent(in) :: nlist   
   character(len=*), intent(in) :: myname
   character(len=*), intent(in) :: flist(nlist)
   character(len=*), intent(in) :: fname
   integer, intent(out) :: id
!
   integer :: n
   character(len=30) :: cdum
!
   id=0
   do n=1,nlist
     read (flist(n),*) cdum
     if (id ==0 .and. trim(cdum) == trim(fname)) then 
       id=n
     endif
   enddo
!
   if (id == 0 .and. myname(1:1) /= ' ') then
     print *,'Error in find_name_2 called from ',myname
     print *,'name to find = ',trim(fname)
     print *,'list to search ='     
     do n=1,nlist
       read (flist(n),*) cdum
       print ('(i3,2x,a)'),n,trim(cdum)
     enddo
     if (lstop) then
       stop
     endif
   endif
!
   end subroutine find_name_2 
