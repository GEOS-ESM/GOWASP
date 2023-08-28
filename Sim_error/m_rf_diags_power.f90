   module m_rf_diags_power
!
! Module used to compute diagnostics of power spectra of correlated random 
! fields, for examination or code testing purposes
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use m_kinds, only : rkind1
   use m_shtrans_rf, only : sh_calc_power
   use m_shtrans_rf, only : rkinds
   use m_shtrans_rf, only : nmax,kmax,nspects
!
   implicit none
!
   private
   public :: rf_diags_power_calc
!
   contains
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_power_calc (nlevs,nfields,kpid,scoefs,field_type,info)
!
! Call routines to calculate and print power spectrum from coefficients of 
! spherical harmonics. 
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: nfields
   integer, intent(in) :: kpid
   complex(rkinds), intent(in) :: scoefs(nspects,nlevs,nfields)
   character(len=*), intent(in) :: field_type
   character(len=*), intent(in) :: info
!
   real(rkind1) :: power(0:kmax,nlevs,nfields)
!
   call sh_calc_power (nlevs,nfields,scoefs,power,field_type)
   call rf_diags_power_print (nlevs,nfields,kpid,power,info)
!
   end subroutine rf_diags_power_calc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
  subroutine rf_diags_power_print (nlevs,nfields,kpid,power,info)
!
! Print power spectra and log of power spectra for up to 7 fields         
!                
   implicit none  
!
   integer :: nlevs,nfields,kpid 
   real(rkind1) :: power(0:kmax,nlevs,nfields)
   character(len=*) :: info
!                      
   integer, parameter :: kpx=9
   integer :: k,klev,km,kmd
   integer :: nf,n
   integer :: klev_p(kpx)
   real(rkind1) :: pk(kpx),sump(kpx)
   real(rkind1) :: pk_avg,sump_avg
!
! Select sample of levels to print
   kmd=max((nlevs-kpx)/kpx,1)
   do k=1,kpx
     klev=1+(k-1)*kmd
     if (klev <= nlevs) then
       km=k
       klev_p(k)=klev
     endif
   enddo
!
   print *,' '
   print ('(a,9i4)'),' Power for levels',klev_p(1:kpx)
   print ('(2a,i4)'),info,':  kpid=',kpid
   print ('(2a,i4)'),' Additional last column is for avg over all', &
                     ' levels: nlevs=',nlevs
   do nf=1,nfields
     print ('(a,i2)'),' Field',nf
     sump_avg=0.
     sump(:)=0.
     do n=0,nmax
       pk_avg=0.
       do klev=1,nlevs
         pk_avg=pk_avg+power(n,klev,nf)
       enddo
       pk_avg=pk_avg/nlevs
       sump_avg=sump_avg+pk_avg
       do k=1,km
         klev=klev_p(k)
         pk(k)=power(n,klev,nf)
         sump(k)=sump(k)+pk(k)
       enddo
       print ('(i4,1p10e11.2)'),n,pk(1:km),pk_avg
     enddo
     print ('(a,1p10e11.2)'),'sum=',sump(1:km),sump_avg
   enddo
!
   end subroutine rf_diags_power_print
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_rf_diags_power
