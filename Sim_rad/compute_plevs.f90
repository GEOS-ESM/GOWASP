   subroutine compute_plevs (kmax,akbk,ps,p)
!
! Compute p at layer interfaces (p(:,1) and at data
! levels (p(:,2)) using ps and the eta-coordinate 
! parameters a(k), b(k) for the vertical grid
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1
   implicit none
!
   real(rkind1), parameter :: Rdryair=287.04
   real(rkind1), parameter :: Cpdryair=1000.5
   real(rkind1), parameter :: Xkappa=Rdryair/Cpdryair
   real(rkind1), parameter :: XkappaI=1./Xkappa
!
   integer, intent(in) :: kmax
   real(rkind1), intent(in)  :: akbk(kmax+1,2)
   real(rkind1), intent(in)  :: ps
   real(rkind1), intent(out) :: p(kmax+1,2)
!
   integer :: k
   real(rkind1) :: logp(kmax+1), pk(kmax+1), dlog, dpk
!
! Determine p at interfaces of vertical grid layers
! using definition of eta coordinate
   do k=1,kmax+1
     p(k,1)=akbk(k,1)+akbk(k,2)*ps
   enddo
!
! Set p on kmax data levels 
! An additional value of a data level kmax+1 is defined as ps
   p(kmax+1,2)=ps
!
   do k=1,kmax+1
     logp(k)=log(p(k,1))
     pk(k)=p(k,1)**Xkappa
   enddo
!
! This formula is based on constraining the integrals of R*T*logp
! to be identical if T=theta*(p/po)**(R/Cp) with theta constant through the 
! layer or if T is constant through the layer, with the relationship between 
! the constants theta and T defined at the particular data level p 
! determined here.
!
   do k=1,kmax
     dlog=logp(k+1)-logp(k)    
     dpk=pk(k+1)-pk(k)
     p(k,2)=(XkappaI*dpk/dlog)**XkappaI
   enddo  
!
   end subroutine compute_plevs 
