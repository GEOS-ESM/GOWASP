   module m_time_compute
!
!  Module to interpret and modify variables that indicate times 
!  associated with observation or field data 
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1, rkind2
   implicit none      

   private
   public :: time_compute_pack
   public :: time_compute_unpack
   public :: time_compute_dhours
   public :: time_compute_add
   public :: time_compute_new_cdtime
   public :: time_compute_secs_since_yd 
!
   public :: rkindh
!
   integer, parameter :: rkindh=rkind1
   integer :: mdays_normal(12)
   character(len=*), parameter :: myname='time_comp'
   data mdays_normal /31,28,31,30,31,30,31,31,30,31,30,31/
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_new_cdtime (cdtime_in,cdtime_out,dhours,ier)
!
!  Add hours to previous date-time variable to create new date-time variable
!
   implicit none
!
   real(rkind1), intent(in) :: dhours
   character(len=*), intent(in)  :: cdtime_in
   character(len=*), intent(out) :: cdtime_out
   integer, intent(out) :: ier
!
   real(rkind1) :: t_in(6)
   real(rkind1) :: t_out(6)
   character(len=*), parameter :: my_name_sub= &
        myname//'::time_compute_new_cdtime'
!
   call time_compute_unpack (cdtime_in,t_in)     
   call time_compute_add (dhours,t_in,t_out,ier)
   call time_compute_pack (t_out,cdtime_out)
!        
   end subroutine time_compute_new_cdtime 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_unpack (cdtime,t)
!
! Extract year, month, day, hour, minutes, seconds from date--time 
! character variable.
!
   implicit none
   real(rkind1), intent(out) :: t(6)
   character(len=*), intent(in) :: cdtime
!
   integer :: lenc
   integer :: iyear,imonth,iday,ihour,iminu,isec
!
   lenc=len(trim(cdtime))
!
   read (cdtime( 1: 4),'(i4)') iyear
   read (cdtime( 5: 6),'(i2)') imonth
   read (cdtime( 7: 8),'(i2)') iday
   read (cdtime( 9:10),'(i2)') ihour
   t(1)=iyear
   t(2)=imonth
   t(3)=iday
   t(4)=ihour
   if (lenc >= 12) then
     read (cdtime(11:12),'(i2)') iminu
     t(5)=iminu
   else
     t(5)=0._rkind1  
   endif
   if (lenc == 14) then
     read (cdtime(13:14),'(i2)') isec
     t(6)=isec
   else
     t(6)=0._rkind1
   endif
!
   end subroutine time_compute_unpack 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_pack (t,cdtime)
!
! Construct date--time character variable from separate variables for 
! year, month, day, hour, minutes, seconds.
!
   implicit none
   real(rkind1), intent(in) :: t(6)
   character(len=*), intent(inout) :: cdtime
!
   integer :: lenc
   integer :: iyear,imonth,iday,ihour,iminu,isec
   integer :: n
!
   lenc=len(cdtime)
!
   iyear =nint(t(1))
   imonth=nint(t(2))
   iday  =nint(t(3))
   ihour =nint(t(4))
   iminu =nint(t(5))
   isec  =nint(t(6))
!
   if (lenc == 10) then 
     write (cdtime,'(i4,3i2)') iyear,imonth,iday,ihour
   elseif (lenc == 12) then 
     write (cdtime,'(i4,4i2)') iyear,imonth,iday,ihour,iminu
   elseif (lenc == 14) then 
     write (cdtime,'(i4,5i2)') iyear,imonth,iday,ihour,iminu,isec
   endif
!
   do n=1,lenc
     if (cdtime(n:n) == ' ') then
       cdtime(n:n)='0'
     endif
   enddo 
!
   end subroutine time_compute_pack 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_dhours (t1,t2,dhours,ier)
!
!  Compute difference t2-t1 in hours
!
   implicit none
!
   real(rkind1), intent(in)  :: t1(6)  ! y,m,d,h,m,s
   real(rkind1), intent(in)  :: t2(6)
   real(rkindh), intent(out) :: dhours
   integer, intent(out) :: ier
!
   real(rkind2) :: dsecs1
   real(rkind2) :: dsecs2
   integer :: ier1,ier2
!
   call time_compute_secs_since_yd (t1,dsecs1,ier1)
   call time_compute_secs_since_yd (t2,dsecs2,ier2)
   ier=ier1+ier2
   dhours=(dsecs2-dsecs1)/3600.
!
   end subroutine time_compute_dhours 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_secs_since_yd (t,dsecs,ier)
!
! Compute seconds between date-time variables and a reference time
! (1 Jan 1901, 0UTC).
!
   implicit none
!
   real(rkind1), intent(in)  :: t(6)
   real(rkind2), intent(out) :: dsecs
   integer, intent(out) :: ier
!
   integer, parameter :: yref=1901 !  1 Jan 1901  0Z is ref. time here
   integer :: i
   integer :: iyear,imonth,iday,ihour,iminu
   integer :: leapyears
   integer :: mdays(12)
   real(rkind2), parameter :: r24=24.0_rkind2
   real(rkind2), parameter :: r60=60.0_rkind2
   real(rkind2), parameter :: r365=365.0_rkind2
   real(rkind2) :: dhours
!
   ier=0
   dsecs=0._rkind2
   iyear=nint(t(1))
!
! Determine if input year is a leap year: works for 1901 <= jyear <= 2099
   mdays=mdays_normal
   if (mod(iyear,4) == 0) then   ! leap year
     mdays(2)=29
   endif
!   
! Determine secs between input year and yref
   iyear=iyear-yref
   leapyears=iyear/4 
   if (iyear < 0 ) then
     ier=1
   else 
     dhours=r24*(real(iyear,rkind2)*r365+real(leapyears,rkind2))
   endif
!
   imonth=nint(t(2))
   if (imonth > 1) then
     do i=1,imonth-1
       dhours=dhours+r24*real(mdays(i),rkind2)          
     enddo
   endif
!
   iday=nint(t(3))
   ihour=nint(t(4))
   dhours=dhours+r24*real(iday-1,rkind2)+real(ihour,rkind2)
!
   iminu=nint(t(5))
   dsecs=t(6)+r60*(real(iminu,rkind2)+r60*dhours)
!
   end subroutine time_compute_secs_since_yd 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine time_compute_add (dhours,tin,tout,ier)
!
! Add a number of hours to previous date-time values to create new 
! date-time values 
!
   implicit none
!
   real(rkindh), intent(in) :: dhours
   real(rkind1), intent(in) :: tin(6)
   real(rkind1), intent(out) :: tout(6)
   integer, intent(out) :: ier
!
   integer :: n
   integer :: iyear,imonth,iday
   integer :: ndays,nhours,nmins
   integer :: mdays(12)
   real(rkindh), parameter :: zero_rh=0._rkindh
   real(rkindh), parameter :: one_rh=1._rkindh
   real(rkindh), parameter :: r24=24._rkind1
   real(rkindh), parameter :: r60=60._rkind1
   real(rkindh), parameter :: r3600=3600._rkind1
   real(rkindh) :: rhours
   real(rkindh) :: rmins
   real(rkindh) :: rsecs
!
   ier=0
!
! Determine if reference is leap year: works for 1901 <= jyear <= 2099
   iyear=nint(tin(1))
   imonth=nint(tin(2))
   mdays=mdays_normal
   if (mod(iyear,4) == 0) then   ! leap year
     mdays(2)=29
   endif
!
   ndays=int(dhours/r24)
   rhours=dhours-real(ndays,rkindh)*r24
   if (rhours < 0._rkindh) then 
     ndays=ndays-1 
     rhours=r24+rhours
   endif
!
   nhours=int(rhours)
   rmins=(rhours-real(nhours,rkindh))*r60
   if (rmins < 0._rkindh) then 
     nhours=nhours-1 
     rmins=r60+rmins
   endif
!
   nmins=int(rmins)
   rsecs=(rmins-real(nmins,rkindh))*r60
   if (rsecs < 0._rkindh) then 
     nmins=nmins-1 
     rsecs=r60+rsecs
   endif
!
   tout(6)=tin(6)+rsecs   
   tout(5)=tin(5)+real(nmins)   
   tout(4)=tin(4)+real(nhours)   
!
! adjust sec
   if (tout(6) < zero_rh) then
     tout(5)=tout(5)-one_rh
     tout(6)=tout(6)+r60
   elseif (tout(6) >= r60) then
     tout(5)=tout(5)+one_rh
     tout(6)=tout(6)-r60
   endif
!
! adjust min
   if (tout(5) < zero_rh) then
     tout(4)=tout(4)-one_rh
     tout(5)=tout(5)+r60
   elseif (tout(5) >= r60) then
     tout(4)=tout(4)+one_rh
     tout(5)=tout(5)-r60
   endif
!
! adjust hour
   if (tout(4) < zero_rh) then
     ndays=ndays-1
     tout(4)=tout(4)+r24
   elseif (tout(4) >= r24) then
     ndays=ndays+1
     tout(4)=tout(4)-r24
   endif
!
! adjust day, month, year
   iday=nint(tin(3))
   if (ndays == 0) then
     tout(3)=tin(3)
     tout(2)=tin(2)
     tout(1)=tin(1)
   elseif (ndays > 0) then
     do n=1,ndays   
       iday=iday+1
       if (iday > mdays(imonth)) then
         iday=1
         imonth=imonth+1
         if (imonth > 12) then
           imonth=1
           iyear=iyear+1
           if (mod(iyear,4) == 0) then   ! leap year
              mdays(2)=29
           else
              mdays(2)=28
           endif
         endif
       endif
     enddo
   elseif (ndays < 0) then
     do n=1,-ndays   
       iday=iday-1
       if (iday == 0) then
         if (imonth > 1) then 
           imonth=imonth-1
           iday=mdays(imonth)
         else 
           imonth=12
           iday=mdays(12)
           iyear=iyear-1
           if (mod(iyear,4) == 0) then   ! leap year
              mdays(2)=29
           else
              mdays(2)=28
           endif
         endif
       endif
     enddo
   endif
!
   tout(3)=real(iday,rkind1)
   tout(2)=real(imonth,rkind1)
   tout(1)=real(iyear,rkind1)
!
   end subroutine time_compute_add 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   end module m_time_compute
