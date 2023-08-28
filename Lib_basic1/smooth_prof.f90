     subroutine smooth_profile (f,klevs,ksmooth)
!
!  Smooth a part of a vertical profile of field f
!  Points at levels k=ksmooth(1),ksmooth(2) are smoothed 
!  Smooth is simple mean of original values j=k-ksmooth(3),k+ksmooth(3)
!  Only works for ksmooth(3)>1
!  Smoothing is only applied in the range of k for which j refers to a valid level
!  Smoothing is only applied in ranges were values are all monotonically 
!  increasing or decreasing with k.
!
     use m_kinds, only: rkind1, zero_k1
     implicit none
!
     integer, intent(in) :: klevs 
     integer, intent(in) :: ksmooth(3)
     real(rkind1) :: f(klevs)
!
     logical :: modf
     integer :: kmin,kmax
     integer :: kspan, kspan2
     integer :: k, k1, k2
     real(rkind1) :: rspan
     real(rkind1) :: fnew(klevs)
     real(rkind1) :: fdiff(klevs)
!
     kspan=ksmooth(3)
     kspan2=kspan*2+1
     rspan=1./real(kspan2)     
     kmin=max(ksmooth(1),kspan+1)
     kmax=min(ksmooth(2),klevs-kspan)
!
     do k=kmin-kspan,kmax+kspan-1
       fdiff(k)=f(k+1)-f(k)
     enddo
!
     do k=kmin,kmax
       k1=k-kspan
       k2=k+kspan-1
       modf=all(fdiff(k1:k2) >= zero_k1) .or. all(fdiff(k1:k2) <= zero_k1)
       if (modf) then
         fnew(k)=rspan*sum(f(k1:k2+1))
       else
         fnew(k)=f(k) 
       endif
     enddo
!
     f(kmin:kmax)=fnew(kmin:kmax)
!
     end subroutine smooth_profile
