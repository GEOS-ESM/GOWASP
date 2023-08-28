   program test
   use m_sft, only : sft_setup
   use m_sft, only : sft_apply
   implicit none
   integer, parameter :: rkind1=4
   integer, parameter :: imax1=32
   integer, parameter :: imax2=32
   integer, parameter :: jmax=2
   integer, parameter :: nlevs=4
!   
   integer :: i
   real(rkind1) :: pi, fac
   real(rkind1) :: f1(imax1,jmax,nlevs)
   real(rkind1) :: f2(imax2,jmax,nlevs)
!
   pi=4._rkind1*atan(1._rkind1)
   do i=1,imax1
     fac=(i-1)*2_rkind1*pi/real(imax1) 
     f1(i,1,1)=cos(fac)
     f1(i,2,1)=sin(fac)
     f1(i,1,2)=cos(fac*4.)
     f1(i,2,2)=sin(fac*4.)
     f1(i,1,3)=sin(fac*16.)
     f1(i,2,3)=cos(fac*16.)
     f1(i,1,4)=cos(fac*0.)
     f1(i,2,4)=sin(fac*0.)
   enddo
!
   print *,' '
   print *,'F1'
   do i=1,imax1
     print ('(i2,8f12.8)'),i,f1(i,:,:)
   enddo
 !
   call sft_setup (imax1,imax2)
   call sft_apply (jmax,nlevs,f1,f2)
 !
   print *,' '
   print *,'F2'
   do i=1,imax2
     print ('(i2,8f12.8)'),i,f2(i,:,:)
   enddo
!
   end program test

