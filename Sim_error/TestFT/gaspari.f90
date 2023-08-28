   program gaspari
   integer, parameter :: nmax=7
   integer, parameter :: rkind=8
   real(rkind) :: f(0:nmax,2,0:2)  ! term order, region, order integrated  
   real(rkind) :: x(2,0:2)         ! log term factor
   real(rkind) :: s(4,0:2)         ! val of fuc at z/c=0, 1(reg 1), 1(reg 2), 2
   real(rkind) :: zc,a,d
   real(rkind) :: c12,c1,c2,sx
!
   f(:,:,:)=0.d0
   x(:,:)=0.d0
!
   m=1
   f(5,m,0)=-0.25d0
   f(4,m,0)= 0.50d0
   f(3,m,0)= 5.0d0/8.0d0
   f(2,m,0)=-5.0d0/3.0d0 
   f(1,m,0)= 0.0d0
   f(0,m,0)= 1.0d0
!
   m=2
   f(5,m,0)= 1.0d0/12.0d0
   f(4,m,0)=-0.50d0
   f(3,m,0)= 5.0d0/8.0d0
   f(2,m,0)= 5.0d0/3.0d0 
   f(1,m,0)=-5.0d0
   f(0,m,0)= 4.0d0
   x(m,0)=  -2.0d0/3.0d0
!
! perform 1st integral
   do n=0,5
     a=real(n+1,rkind)
     do m=1,2
       f(n+1,m,1)=f(n,m,0)/a
     enddo
   enddo
   f(0,1,1)=1.d0   ! So that function = 1 at z/c = 0
   x(2,1)=x(2,0)
!
! perform 2nd integral
   do n=0,6
     a=real(n+1,rkind)
     do m=1,2
       f(n+1,m,2)=f(n,m,1)/a
     enddo
   enddo
   f(0,1,2)=1.d0   ! So that function = 1 at z/c = 0
   x(2,2)=x(2,0)
!   
   print *,'Region 0<= z/c <= 1'
   do n=7,0,-1
     print ('(a,i1,7x,3f15.8)'),' (z/c)**',n,f(n,1,:)    
   enddo
!
   print *,' '
   print *,'Region 1<= z/c <= 2'
   do n=7,0,-1
     print ('(a,i1,10x,3f15.8)'),' (z/c)**',n,f(n,2,:)    
   enddo
   print ('(a,7x,f15.8)'),' (z/c)**(-1)',x(2,0)    
   print ('(a,22x,f15.8)'),' log(z/c)   ',x(2,1)    
   print ('(a,31x,f15.8)'),' (z/c)*(log(z/c)-1)',x(2,2)    
!
! Evaluate each integral function (0,1,2) at z/c=0 
   s(1,:)=f(0,1,:)   ! = term of order 0
! Evaluate each integral function (0,1,2) at z/c=1 for region 1
   s(2,:)=0.d0       ! initialize accumulator to 0 
! Evaluate each integral function (0,1,2) at z/c=1 for region 2
   s(3,0)=x(2,0)     ! contribution by (z/c)**-1 term
   s(3,1)=0.d0       ! contribution by log(z/c) term 
   s(3,2)=-x(2,2)    ! contribution by (z/c)(log(z/c)-1) term 
! Evaluate each integral function (0,1,2) at z/c=2
   s(4,0)=x(2,0)/2.d0       ! contribution by (z/c)**-1 term 
   s(4,1)=x(2,1)*log(2.d0)  ! contribution by log(z/c) term 
   s(4,2)=x(2,2)*2.d0*(log(2.d0)-1.d0)  ! contribution by (z/c)(log(z/c)-1)
   do n=0,7
     zc=real(2**n,rkind)
     do i=0,2
       s(2,i)=s(2,i)+f(n,1,i)
       s(3,i)=s(3,i)+f(n,2,i)
       s(4,i)=s(4,i)+f(n,2,i)*zc
     enddo
   enddo
   print *,' '
   print *,' '
   print ('(a,4f14.8)'),'Orig function edges: ',s(1:4,0)
   print ('(a,4f14.8)'),'1st integral  edges: ',s(1:4,1)
   print ('(a,4f14.8)'),'2nd integral  edges: ',s(1:4,2)
!
   c12=s(2,2)-s(3,2)
   c1=-(c12+s(4,2))
   c2=c12-c1
   f(0,2,1)=c1
   f(0,2,2)=c2
   f(1,2,2)=c1
! 
   print *,' '  
   print *,' '  
   print *,'Region 0<= z/c <= 1'
   do n=7,0,-1
     print ('(a,i1,7x,3f15.8)'),' (z/c)**',n,f(n,1,:)    
   enddo
!
   print *,' '
   print *,'Region 1<= z/c <= 2'
   do n=7,0,-1
     print ('(a,i1,10x,3f15.8)'),' (z/c)**',n,f(n,2,:)    
   enddo
   print ('(a,7x,f15.8)'),' (z/c)**(-1)',x(2,0)    
   print ('(a,22x,f15.8)'),' log(z/c)   ',x(2,1)    
   print ('(a,31x,f15.8)'),' (z/c)*(log(z/c)-1)',x(2,2)    
!
! Evaluate each integral function (0,1,2) at z/c=0 
   s(1,:)=f(0,1,:)   ! = term of order 0
! Evaluate each integral function (0,1,2) at z/c=1 for region 1
   s(2,:)=0.d0       ! initialize accumulator to 0 
! Evaluate each integral function (0,1,2) at z/c=1 for region 2
   s(3,0)=x(2,0)     ! contribution by (z/c)**-1 term
   s(3,1)=0.d0       ! contribution by log(z/c) term 
   s(3,2)=-x(2,2)    ! contribution by (z/c)(log(z/c)-1) term 
! Evaluate each integral function (0,1,2) at z/c=2
   s(4,0)=x(2,0)/2.d0       ! contribution by (z/c)**-1 term 
   s(4,1)=x(2,1)*log(2.d0)  ! contribution by log(z/c) term 
   s(4,2)=x(2,2)*2.d0*(log(2.d0)-1.d0)  ! contribution by (z/c)(log(z/c)-1)
   do n=0,7
     zc=real(2**n,rkind)
     do i=0,2
       s(2,i)=s(2,i)+f(n,1,i)
       s(3,i)=s(3,i)+f(n,2,i)
       s(4,i)=s(4,i)+f(n,2,i)*zc
     enddo
   enddo
   print *,' '
   print *,' '
   print ('(a,4f14.8)'),'Orig function edges: ',s(1:4,0)
   print ('(a,4f14.8)'),'1st integral  edges: ',s(1:4,1)
   print ('(a,4f14.8)'),'2nd integral  edges: ',s(1:4,2)
!
! print final 2nd integral function
   print *,' '
   print *,' '
   print *,' final 2nd integral function='
   nn=20
   d=1./real(nn,rkind)
   do i=0,nn 
     zc=i*d
     sx=0.d0
     do n=0,7
       sx=sx+f(n,1,2)*zc**n
     enddo
     print ('(i4,2f15.8)'),i,zc,sx
   enddo
   do i=nn,2*nn 
     zc=i*d
     sx=x(2,2)*zc*(log(zc)-1.d0)
     do n=0,7
       sx=sx+f(n,2,2)*zc**n
     enddo
     print ('(i4,2f15.8)'),i,zc,sx
   enddo
!
   end program gaspari



