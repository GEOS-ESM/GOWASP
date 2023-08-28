   program satwind_ods2txt 
!
! Input ods file.  Output text file of satwind meta data.
!
   use m_ods_RE
!
   implicit none
!
   integer :: n
   integer :: n_date, n_time
   integer :: n_obs
   integer :: n_cnt
   integer :: argc
   integer(4) :: iargc
   integer :: ierr
   integer, parameter :: txt_unit_out=21
   integer :: iout(5)
!
   real(4) :: x
!
   character(len=14)  :: cdtime0
   character(len=240) :: txt_file_out     ! prepbufr file of satwind obs
   character(len=240) :: ods_file_in
   character(len=30)  :: cdum_type  ! output from Get_ODS 
!
   type (ods_vect) :: ods
!
! Read arguments defined when program is executed in script
   argc = iargc()
   if (argc /= 3) then
     print *,'Number of arguments /= 3'
     stop
   endif
!
   call GetArg( 1_4, cdtime0)
   call GetArg( 2_4, ods_file_in)
   call GetArg( 3_4, txt_file_out)
!
   read (cdtime0(1:8),'(i8)')  n_date                                         
   read (cdtime0(9:14),'(i6)') n_time 
!
   call ods_get_RE (.true.,trim(ods_file_in),n_time,n_obs,ods,ierr)
   if (ierr /= 0) then
     print *,'Problems detected when attempting read of ods file: ierr=',ierr
     stop
   endif
   print *,'number of obs in file =',n_obs
!
   open(unit=txt_unit_out,file=trim(txt_file_out),form='formatted')
   print ('(3a,i3)'),' bin output file=',trim(txt_file_out),   &
            ' opened on unit=',txt_unit_out
!
   n_cnt=0
   do n=1,n_obs
     if ( (ods%data%kx(n) > 239) .and. (ods%data%kx(n) < 261) .and.    &
          (ods%data%kt(n) == 4)  .and. (ods%data%qcexcl(n) == 0) ) then
!
       n_cnt=n_cnt+1
       x=mod(ods%data%lon(n),360.)
       if (x<0.) x=x+360.
       iout(1)=ods%data%kx(n)-200
       iout(2)=nint(x*100.)
       iout(3)=nint(ods%data%lat(n)*100.)
       iout(4)=ods%data%time(n)
       iout(5)=nint(ods%data%lev(n))
       write (txt_unit_out,'(i2,2i5,2i4)') iout(:)
!    
     endif ! check if acceptible obs in ods
   enddo   ! loop over n_obs
!
   iout(1)=-9
   iout(2:3)=-9999
   iout(4:5)=-999
   write (txt_unit_out,'(i2,2i5,2i4)') iout(:)
   close (txt_unit_out)   
   print *,' '
   print *,'Output bufr file written: n_cnt=',n_cnt   
!
   call ods_clean_RE (ods)
!
   print *,' '
   print *,'program completed' 
!
   end program satwind_ods2txt 
