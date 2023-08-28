   subroutine prof_all_2d3d (num_2d,num_3d,kmax,num_all,prof_all, &
                             prof_2d,prof_3d)
!  
!  Create separate arrays of 2d fields and 3d fields from a single array 
!  of profile data.
!
!  Initial Code:  Ronald Errico  September 15 2014
!
   use m_kinds, only : rkind1
   implicit none
   integer, intent(in) :: num_2d,num_3d,kmax,num_all
   real(rkind1), intent(in) :: prof_all(num_all)
   real(rkind1), intent(out) :: prof_2d(num_2d), prof_3d(kmax,num_3d)
!
   integer :: k,k1,nf
! 
! Copy 2d field values
   prof_2d(1:num_2d)=prof_all(1:num_2d)
!
! Copy 3d field values
   k1=num_2d
   do nf=1,num_3d
     do k=1,kmax
       k1=k1+1
       prof_3d(k,nf)=prof_all(k1)
     enddo
   enddo
!
   end subroutine prof_all_2d3d
