   module m_prnt_test
!
   implicit none
   private
   public :: ptest_open
   public :: ptest_fill
   public :: ptest_close
!
   integer, parameter :: iunit=41
   integer :: icount
   real(4) :: fs(5,100)
   character(len=12) :: fnames(100)    
   character(len=6) :: cnames(100)    
!
   contains   
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!   
   subroutine ptest_open (ltvq,caorb,scaleT,scaleQ,scaleW,scaleP)
!
   implicit none
!                          
   logical :: ltvq
   real(4) :: scaleT,scaleQ,scaleW,scaleP
   character(len=*) :: caorb
!
   open (iunit,file='enorm_prnt_test')
   write (iunit,'(a,l8,2a)') 'ltvq= ',ltvq,'  caorb=',caorb 
   write (iunit,'(a,1p5e15.4)') 'scalePQTWW= ', &
          scaleP,scaleQ,scaleT,scaleW,scaleW
! scales written in order of corresponding field names
   cnames(:)=' '    
   fnames(:)=' '    
   fs(:,:)=0.    
   icount=0
!
   end subroutine ptest_open 
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!   
   subroutine ptest_fill (iprnt,imax,jmax,kmax,cname,fname,f)
!                          
   implicit none 
!
   integer :: iprnt
   integer :: imax,jmax,kmax
   real(4) :: f(imax,jmax,kmax)
   character(len=*) :: fname
   character(len=*) :: cname
!
   integer :: kp
!
   icount=max(icount,iprnt)
   kp=min(60,kmax)
   cnames(iprnt)=trim(cname)
   fnames(iprnt)=trim(fname)
   fs(1,iprnt)=f(1,100,kp)
   fs(2,iprnt)=f(1,150,kp)
   fs(3,iprnt)=f(1,200,kp)
   fs(4,iprnt)=f(1,250,kp)
   fs(5,iprnt)=f(1,300,kp)
!
   end subroutine ptest_fill
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!   
   subroutine ptest_close 
!
   implicit none 
   integer :: i
!                          
   do i=1,icount
     write (iunit,'(2a,1p5e14.5)') cnames(i),fnames(i),fs(:,i)
   enddo
   close (iunit)
!
   end subroutine ptest_close 
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!   
   end module m_prnt_test
