program test
integer, parameter :: nmax=180*4 
real(8) x(0:nmax),a,L,z,s,y,w

L=700.
a=6371.
z=-0.5*(L/a)**2
s=0.
do n=0,nmax
y=real(n*n+n,8)
w=real(2*n+1,8)
x(n)=sqrt(w)*exp(z*y)
s=s+x(n)
s=1.
enddo
x=x/s
open (10,file='testx.txt')
do n=0,nmax
write (10,'(i4,1p1e10.2)'),n,x(n)
enddo
close (10)
end program test 
