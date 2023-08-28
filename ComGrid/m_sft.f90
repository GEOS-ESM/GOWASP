!
   module m_sft

   implicit none
   private
   public :: sft_setup
   public :: sft_apply
   integer, parameter :: rkind1=4
   integer, parameter :: rkind2=8
   integer :: imaxNR,kmaxNR
   integer :: imaxCG,kmaxCG
   real(rkind2), allocatable :: csNR(:,:,:)
   real(rkind2), allocatable :: csCG(:,:,:)
   real(rkind2) :: rfac
!
   contains
!
   subroutine sft_setup (nlonNR,nlonCG)
!
   implicit none
!
   integer, intent(in) :: nlonNR,nlonCG
   integer :: i,k
   real(rkind2) :: pi, fac1, fac2, pfac
!   
   pi=4._rkind2*atan(1._rkind2)
   imaxNR=nlonNR
   imaxCG=nlonCG
   kmaxNR=imaxNR/2
   kmaxCG=imaxCG/2
   rfac=1._rkind2/real(imaxNR)
   allocate (csNR(imaxNR,0:kmaxCG,2))
   allocate (csCG(imaxCG,0:kmaxCG,2))
!
   fac1=2._rkind2*pi/real(imaxNR)
   do i=1,imaxNR
     fac2=fac1*real(i-1)
     do k=0,kmaxCG  
       pfac=fac2*real(k)
       csNR(i,k,1)=cos(pfac)
       csNR(i,k,2)=sin(pfac)
     enddo
   enddo
!
   fac1=2._rkind2*pi/real(imaxCG)
   do i=1,imaxCG
     fac2=fac1*real(i-1)
     do k=0,kmaxCG  
       pfac=fac2*real(k)
       csCG(i,k,1)=cos(pfac)
       csCG(i,k,2)=sin(pfac)
     enddo
   enddo
!
   end subroutine sft_setup 
!
!
   subroutine sft_apply (jmax,nlevs,fNR,fCG)
!
   implicit none
   integer, intent(in) :: jmax,nlevs
   real(rkind1) :: fNR(imaxNR,jmax,nlevs)
   real(rkind1) :: fCG(imaxCG,jmax,nlevs)
!
   integer :: i,j,k,n
   real(rkind2) :: coefs(2), sfac
!
   fCG(:,:,:)=0._rkind1
   do n=1,nlevs
     do j=1,jmax
       do k=0,kmaxCG
!
         if (k==0 .or. k==kmaxNR) then
           sfac=1._rkind2
         else
           sfac=2._rkind2
         endif
!
         coefs(:)=0._rkind2
         do i=1,imaxNR
           coefs(1)=coefs(1)+fNR(i,j,n)*csNR(i,k,1)
           coefs(2)=coefs(2)+fNR(i,j,n)*csNR(i,k,2)
         enddo
!
         coefs(:)=rfac*coefs(:)
!
         do i=1,imaxCG
           fCG(i,j,n)=fCG(i,j,n)+sfac* &
               (coefs(1)*csCG(i,k,1)+coefs(2)*csCG(i,k,2))
         enddo
!        
       enddo
     enddo
   enddo
!
   end subroutine sft_apply 
!
!
   end module m_sft
