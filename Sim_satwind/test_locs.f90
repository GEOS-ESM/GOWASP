   program test_locs
   implicit none
   integer :: nlons,nlats,kx_num
   integer :: i,j,k,jx,kx
   integer, parameter :: iunit=10
   integer, allocatable :: ijlocs(:,:,:)
   character (len=1), allocatable :: cx(:,:,:)
   character(len=1), parameter :: c1='*'
   character(len=1), parameter :: c0=' '
   character(len=*), parameter :: &
               filename='/discover/nobackup/rerrico/test_locs'
!
   open (iunit,file=filename,form='unformatted')
   read (iunit) nlons,nlats,kx_num
   print *,'nlons,nlats,kx_num=',nlons,nlats,kx_num
   allocate (ijlocs(nlons,nlats,kx_num))
   allocate (cx(nlons,nlats,kx_num))
!   
   do k=1,kx_num
     read (iunit) kx
     do j=1,nlats
       read (iunit) jx,ijlocs(:,j,k)
     enddo
     print *,'kx read=',kx
   enddo  
   close (iunit)
!
   do k=1,kx_num
     do j=1,nlats
       do i=1,nlons
         if (ijlocs(i,j,k) > 0) then
           cx(i,j,k)=c1
         else
           cx(i,j,k)=c0
         endif 
       enddo
     enddo
   enddo
!
   open (iunit,file='/discover/nobackup/rerrico/test_locs.txt')
   do k=1,kx_num
     write (iunit,'(a1)') ' '
     write (iunit,'(a,i2)') 'KXS=',k
     do i=1,nlons
       write (iunit,'(a2,i4,361a1)') 'i=',i,cx(i,:,k)
     enddo
   enddo
   close (iunit)
!
   end program test_locs
