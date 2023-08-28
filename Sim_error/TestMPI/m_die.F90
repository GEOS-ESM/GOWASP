   module m_die
!
   implicit none
   private
   public :: mpi_die
   public :: die_proc_id
!
   include "mpif.h"
   integer :: die_proc_id
!
   contains
!      
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine mpi_die (d_name,d_ier)
!
   implicit none
   integer, intent(in) :: d_ier
   character(len=*), intent(in) :: d_name
   integer :: ier
!
   print *,' '
   print *,' Fatal error detected: Program aborting for all processors'
   print ('(a,i2,a,i3,2a)'),' proc_id=',die_proc_id,'  error=',d_ier, &
                            '  name=',d_name
   call mpi_abort(MPI_comm_world,d_ier,ier)
   if (ier /= 0) then
     print *,'XXXXXXXXXXXXXXXX'
     print *,'Attempt to abort all processors failed: ier=',ier
   endif
!
   end subroutine mpi_die 
!
   end module m_die
