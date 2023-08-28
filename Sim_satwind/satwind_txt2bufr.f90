   program satwind_txt2bufr 
!
   implicit none
!
   integer, parameter :: rkind1=4
   integer :: n
   integer :: idtime
   integer :: ncnt
   integer :: argc
   integer(4) :: iargc
   integer :: ierr
   integer :: iout(5)
   integer, parameter :: obs_data_dim1=9
   integer, parameter :: txt_unit_in=10
   integer, parameter :: bufr_unit_in=20
   integer, parameter :: bufr_unit_out=21
   integer, parameter :: pts_in_msg=200 ! max # of obs reports per message
!
   real(rkind1) :: obs_data(obs_data_dim1)
!
   character(len=14) :: cdtime0
   character(len=240) :: bufr_file_in      ! file with prepbufr table
   character(len=240) :: bufr_file_out     ! prepbufr file of satwind obs
   character(len=240) :: txt_file_in
!
! Read arguments defined when program is executed in script
   argc = iargc()
   if (argc /= 4) then
     print *,'Wrong number of arguments'
     stop
   endif
!
   call GetArg( 1_4, cdtime0)
   call GetArg( 2_4, txt_file_in)
   call GetArg( 3_4, bufr_file_in)
   call GetArg( 4_4, bufr_file_out)
!
   open(unit=bufr_unit_in,file=trim(bufr_file_in),form='unformatted')
   call openbf(bufr_unit_in,'IN ',bufr_unit_in)
   print *,' '
   print ('(3a,i3)'),' bufr input file=',trim(bufr_file_in),     &
            ' opened on unit=',bufr_unit_in
!
! Open output bufr file
   open(unit=bufr_unit_out,file=trim(bufr_file_out),form='unformatted')
   print ('(3a,i3)'),' bufr output file=',trim(bufr_file_out),   &
            ' opened on unit=',bufr_unit_out
   call openbf(bufr_unit_out,'OUT',bufr_unit_in)  
!
   open(unit=txt_unit_in,file=trim(txt_file_in),form='formatted')
!
   ncnt=0
   do n=1,10000000
     read (txt_unit_in,'(i2,2i5,2i4)') iout(:)   
     if (iout(1) > 0) then
       ncnt=ncnt+1
       obs_data(1)=ncnt
       obs_data(2)=real(iout(2))/100.   ! lon 
       obs_data(3)=real(iout(3))/100.   ! lat
       obs_data(4)=real(iout(4))/60.    ! time (hrs) relative to center
       obs_data(5)=iout(1)+200          ! kx
       obs_data(6)=0.                   ! z
       obs_data(7)=real(iout(5))*100.   ! p (Pa)
       obs_data(8)=0.                   ! u
       obs_data(9)=0.                   ! v
!
       if (mod(ncnt,pts_in_msg) == 1) then   ! group reports into messages 
         call openmb (bufr_unit_out,'SATWND',idtime)
       endif          
!
! Write one obs bufr report 
       call bufr_satwind_write (bufr_unit_out,obs_data_dim1,obs_data)
       if (mod(ncnt,pts_in_msg) == 0) then  
         call closmg(bufr_unit_out)
       endif 
!    
     else
       exit
     endif
   enddo
!
   if (mod(ncnt,pts_in_msg) /= 0) then   
     call closmg (bufr_unit_out)  
   endif 
!
   call closbf (bufr_unit_out)  
   call closbf (bufr_unit_in)  
   print *,' '
   print *,'Output bufr file written for ncnt=',ncnt   
!
   print *,' '
   print *,'program completed' 
!
   end program satwind_txt2bufr 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine bufr_satwind_write (luout,obs_ndim,obs_info) 
!
!  Write satwind observations in .prepbufr (BUFR) format
!
!  Initial code Ronald Errico July 2014
!
   implicit none
!
   integer, parameter :: rkind1=4
   integer,intent(in) :: luout    ! unit number for bufr ouput
   integer,intent(in) :: obs_ndim ! dim of obs_info array
   real(rkind1), intent(in) :: obs_info(obs_ndim)
!
! local variables
   integer, parameter :: nhead=8      !  size of header
   integer :: llev
   real(rkind1), parameter :: pscale=1.e2
   real(8) :: hdr(nhead)    !  location and metadata
   real(8) :: uvout(5)      !  output wind array
   real(8) :: boutp(5)      !  output for pevstr 
   real(8) :: boutz(4)      !  output for zevstr 
!
   character(len=8)  :: station_id
   character(len=44) :: hdstr,uevstr,pevstr,zevstr
      data hdstr     /'SID XOB YOB DHR TYP ELV T29 SAID            ' /
      data pevstr    /'POB PQM PPC PRC CAT                         ' /
      data uevstr    /'UOB VOB WQM WPC WRC                         ' /
      data zevstr    /'ZOB ZQM ZPC ZRC                             ' /
!
     write (station_id,'(i8)') nint(obs_info(1))  ! expess id as character
     hdr(1)=transfer(station_id,hdr(1))           ! reformat as type hdr
     hdr(2:6)=obs_info(2:6)
     hdr(7)=63.                                   ! indicates satwind obs
     hdr(8)=1.d10                                 ! missing value 
!
! Write out p values
     boutp(1)=obs_info(7)/pscale
     boutp(2)=2._8
     boutp(3)=1._8               
     boutp(4)=0._8               
     boutp(5)=6._8               
     call ufbint (luout,boutp,5,1,llev,pevstr)
!
! Write out u and v values
     uvout(1)=obs_info(8)
     uvout(2)=obs_info(9)
     uvout(3)=2._8
     uvout(4)=4._8
     uvout(5)=0._8                 ! reason code
     call ufbint(luout,uvout,5,1,llev,uevstr)
!
! write out z values
     boutz(1)=obs_info(6)
     boutz(2)=2._8
     boutz(3)=1._8
     boutz(4)=0._8
     call ufbint(luout,boutz,4,1,llev,zevstr)
!
!  Finally, write out the header information.
     call ufbint(luout,hdr,nhead,1,llev,hdstr)
!
!  Finished updates; so write output buffer to output file for this report
!     
     call writsb(luout)
!
     end subroutine bufr_satwind_write

