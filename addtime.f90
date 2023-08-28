
      Program Add_time
!
!  increment date/time in the form yyyymmddhh
!  e.g., datetime 2006010100 +036 will tield 2006 01 02 12 
!
      implicit none
!
      character(len=10) :: c_datetime
      character(len=13) :: c_datetime_out
      character(len=4 ) :: c_addhours
      logical :: month_error
      integer :: nd
      integer :: jyear
      integer :: jmonth
      integer :: jday
      integer :: jhour
      integer :: addhours
      integer :: adddays
      integer :: ierror
!
      integer :: ichar
      integer :: iyear
      integer :: mdays(12)
      data mdays /31,28,31,30,31,30,31,31,30,31,30,31/
!
      ierror=0
!
      call GetArg( 1_4,c_datetime )
      call GetArg( 2_4, c_addhours)
      read (c_datetime,'(i4,3i2)') jyear,jmonth,jday,jhour
      read (c_addhours,'(i4)')  addhours
!
! Determine if leap year: works for 1901 <= jyear <= 2099
      if (mod(jyear,4) == 0) then   ! leap year
        mdays(2)=29
      endif
!
! Check that input date time is ok 
      if ( (jmonth < 1) .or. (jmonth > 12) ) then
        ierror=ierror+2
        month_error=.true.
      else
        month_error=.false.
      endif
!
      if (.not. month_error) then  ! only check day if month is ok
        if ( (jday < 1) .or. (jday > mdays(jmonth)) ) then
          ierror=ierror+4
        endif
      endif    
!
      if ( (jhour < 0) .or. (jhour > 24) )then
        ierror=ierror+8
      endif
!
! To decipher the error value, convert to binary: Then
!    unit diget = 1 if year  is unacceptable
!    2's  diget = 1 if month is unacceptable
!    4's  diget = 1 if day   is unacceptable
!    8's  diget = 1 if hour  is unacceptable
!
      if (ierror /= 0) then
        write (c_datetime,'(a6,i4)') 'error=',ierror
        write (*,*) c_datetime
        stop
      endif
      adddays=addhours/24
      addhours=mod(addhours,24)
!
      if (adddays>0) then
        do nd=1,adddays
          jday=jday+1
          if (jday>mdays(jmonth)) then
            jmonth=jmonth+1
            jday=1
            if (jmonth>12) then
              jmonth=jmonth-12
              jyear=jyear+1
            endif
          endif
        enddo
      endif
!
      if (adddays<0) then
        do nd=adddays,-1
          jday=jday-1
          if (jday<1) then
            jmonth=jmonth-1
            if (jmonth<1) then
              jmonth=12
              jyear=jyear-1
            endif
            jday=mdays(jmonth)
          endif  
        enddo
      endif       
!      
      jhour=jhour+addhours
      if (jhour > 23) then
        jhour=jhour-24 
        jday=jday+1
        if (jday>mdays(jmonth)) then
          jmonth=jmonth+1
          jday=1
          if (jmonth>12) then
            jmonth=jmonth-12
            jyear=jyear+1
          endif
        endif
      endif
      if (jhour < 0) then
        jhour=jhour+24 
        jday=jday-1
        if (jday<1) then
          jmonth=jmonth-1
          if (jmonth<1) then
            jmonth=12
            jyear=jyear-1
          endif
          jday=mdays(jmonth)
        endif  
      endif       
!         
      write (c_datetime_out,'(i4,3i3)') jyear,jmonth,jday,jhour
      do ichar=6,12,3
        if (c_datetime_out(ichar:ichar)==' ') then
          c_datetime_out(ichar:ichar)='0'
        endif
      enddo
      write (*,*) c_datetime_out
!
      end program Add_time


