   module m_count_types
!
!  Count and then print the numbers of obs processed for each data subtype
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   implicit none
   private
   public :: count_types_setup
   public :: count_subtypes
   public :: count_subtypes_get_rpts
   public :: count_subtypes_print 
!
   integer, parameter :: ntypes=50   ! maximum number of types
   integer :: count_types(ntypes,3)
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
     subroutine count_types_setup 
     count_types(:,:)=0
     end subroutine count_types_setup 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
     subroutine count_subtypes (itype,obs_num)
!
!  Fill table of subtype counters
!
     implicit none   
     integer, intent(in) :: itype, obs_num
!     
     integer :: n, n1
!
     n1=0                                 ! flag (0 means type not found yet)
     do n=1,ntypes                        ! search through list
       if (n1 == 0) then                  ! correct type not found in list yet
         if (count_types(n,1) == 0) then  ! fill type since slot is empty 
           count_types(n,1)=itype
         endif 
         if (itype == count_types(n,1)) then ! check if type found in list 
           n1=n                           ! type found in list
         endif
       endif
     enddo
!
     count_types(n1,2)=count_types(n1,2)+1
     count_types(n1,3)=count_types(n1,3)+obs_num
!
     end subroutine count_subtypes 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
     subroutine count_subtypes_get_rpts (itype,num_rpts)
!
!  Get number of reports for requested type
!
     implicit none   
     integer, intent(in)  :: itype
     integer, intent(out) :: num_rpts
!     
     integer :: n, n1
!
     n1=0                                 ! flag (0 means type not found yet)
     do n=1,ntypes                        ! search through list
       if (n1 == 0) then                  ! correct type not found in list yet
         if (itype == count_types(n,1)) then ! check if type found in list 
           n1=n                           ! type found in list
         endif
       endif
     enddo
!
     if (n1 > 0) then 
       num_rpts=count_types(n1,2)
     else
       num_rpts=-1   ! type not found in list provided in count_types(:,1)
     endif
!
     end subroutine count_subtypes_get_rpts
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
     subroutine count_subtypes_print 
!
!  Print table of subtype counters
!
     implicit none   
!     
     integer :: n, n1, n_next_min, n_filled
     integer :: itype, itype_max
     integer :: count_types_reord(ntypes,3)
!
! Count number of filled table entries
     do n=1,ntypes
       if (count_types(n,1) > 0) then
         n_filled=n
       endif
     enddo
!    
! Reorder filled entries 
     itype_max=100000 
     do n=1,n_filled     ! index for re-ordered table entries
       itype=itype_max 
       do n1=1,n_filled  ! search filled values in original table
         if (count_types(n1,1) < itype) then
           n_next_min=n1
           itype=count_types(n1,1) 
         endif
       enddo
       count_types_reord(n,:)=count_types(n_next_min,:)
       count_types(n_next_min,1)=itype_max
     enddo
!
! Print table
     print *,' '
     print *,' Table of subtype report and observation value counters:'
     do n=1,n_filled                    
       print ('(3(a,i9))'), '   Subtype=',count_types_reord(n,1),   &   
                            '   Rep Count=',count_types_reord(n,2), &
                            '   Obs Count=',count_types_reord(n,3)
     enddo
!
     end subroutine count_subtypes_print 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
     end module m_count_types
