  module m_set_unit_nums
!
! Set some unit numbers for Fortran reading/writing 
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
  implicit none
  public
!
! In the following reserve +1 for separate bufr table file if needed
  integer, parameter :: un_bufrin=12           ! input bufr file 
  integer, parameter :: un_bufrin2=14          ! input 2nd bufr file
  integer, parameter :: un_bufrout=16          ! output bufr file  
!
  integer, parameter :: un_info=10             ! input for rc files
!
  integer, parameter :: un_list_write=11
  integer, parameter :: un_list_in0=20
  integer, parameter :: un_prof_out0=40
  integer, parameter :: un_prof_in1=18
  integer, parameter :: un_prof_in2=19
  integer :: un_list_in
  integer :: un_prof_out
!
  end module m_set_unit_nums
