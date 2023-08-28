   program sum_table
!
! Sum counts in kx_table (sum over all p for each j, and also 
! compute the total; separate sums for each kx,ks) 
!
   use m_kx_table 
   integer :: ier
   character(len=240) :: cfilein
   character(len=240) :: cfileout
!
   call GetArg( 1_4, cfilein)
   call GetArg( 2_4, cfileout)
   call kx_table_read (cfilein,.true.,5,ier)
   call kx_table_sum (cfilein,cfileout)
! 
   end program sum_table
