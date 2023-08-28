   subroutine bufr_satwind_write (luout,obs_ndim,obs_info) 
!
!  Write satwind observations in .prepbufr (BUFR) format
!
!  Initial code Ronald Errico July 2014
!
   use m_kinds, only : rkind1
   implicit none
!
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
     if (obs_info(10) > 0.1) then
       hdr(8)=obs_info(10)             ! sat id#
     else
       hdr(8)=1.d10                    ! missing value (polar sats) 
     endif               
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
