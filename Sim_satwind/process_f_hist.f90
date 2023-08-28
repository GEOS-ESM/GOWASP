   program process_f_hist
!
! Print the computed histograms of cf or sacled delta ipw for 
! diagnostic purposes.  This NEEDS SOME CHANGES
!
   implicit none
!
   integer, parameter :: fileformat=1
   integer, parameter :: iunitin=10
   integer, parameter :: iunitout=11
   integer :: i,j,k,n
   integer :: ip,jp,kp,np
   integer :: kx_list,kx_satid,hsum3
   integer :: kx_cbins,kx_pbins,kx_jbins,kx_num
   integer, allocatable :: kx_hist(:)
   integer :: hsum1,hsum2
   real(4), allocatable :: percnt(:)

   character(len=30) :: cdum
   character(len=*), parameter :: filein='histogram_test' 
   character(len=*), parameter :: fileout='histogram_pcnt.txt' 
!
   open (iunitin,file=trim(filein))
   open (iunitout,file=trim(fileout))
!   read (iunitin,'(*)') cdum,fileformat
!   read (iunitin,'(*)') cdum,kx_cbins,kx_pbins,kx_jbins,kx_num
   kx_cbins=21
   kx_pbins=7
   kx_jbins=12
   kx_num=19
   allocate (kx_hist(kx_cbins))
   allocate (percnt(kx_cbins))
   read (iunitin,'(a1)') cdum(1:1)  ! skip labels
   do k=1,kx_num
     do j=1,kx_jbins 
       do n=1,kx_pbins 
         read (iunitin,'(i2,2i3,2i5,i9)') kp,jp,np,kx_list,kx_satid,hsum3
print *,kp,jp,np,kx_list,kx_satid,hsum3
         read (iunitin,'(12x,21i6)') kx_hist(2:) ! do not read first 
print *,kx_hist(1),kx_hist(kx_cbins)
         read (iunitin,'(a1)') cdum(1:1)
         read (iunitin,'(a1)') cdum(1:1)
print *,cdum(1:1)
         hsum2=sum(kx_hist(2:))
         hsum1=hsum3
         do i=2,kx_cbins
           percnt(i)=kx_hist(i)/real(max(hsum2,1))
         enddo 
         percnt(1)=(hsum1-hsum2)/real(max(hsum1,1))
         write (iunitout,'(i2,2i3,2i5,i9,f8.4)') k,j,n,kx_list,kx_satid,hsum3,percnt(1)
         write (iunitout,'(20f6.3)') percnt(2:)
       enddo
     enddo
   enddo
   close (iunitin)
   close (iunitout)
   print *,'Histogram printed to file :',trim(fileout)
!   
   end program process_f_hist
