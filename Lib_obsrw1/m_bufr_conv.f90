!
   Module m_bufr_conv
!
! Module to read and write conventional observation data in BUFR format.
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1,rkind2
!
   use m_conv_names, only : nheadc, nheadm, nheadw
   use m_conv_names, only : conv_nhead, conv_nfields, conv_max_levs
   use m_conv_names, only : bbx, bby, bbr, bbp, bbz, bbu, bbv, bbt, bbq
   use m_conv_names, only : bbc, bbf, bbpq, bbzq, bbwq, bbtq, bbqq
   use m_conv_names, only : ityp, idhr, ixob, iyob, isid, ielv, ilzf, iwsd
   use m_conv_names, only : hdstrc, hdstrm, hdstrw, hdstrm_extra, hdstrw_extra
   use m_conv_names, only : hdstr_merge
   use m_conv_names, only : obs_info_names, conv_value_names
!
   implicit none
!
   private
   public :: conv_rw_setup
   public :: conv_rw
   public :: conv_check_obs
   public :: conv_print_sample
!
   logical :: l_z_same    ! indicates that station and obs elev same
   logical :: l_z_present ! indicates that heights are present 
   logical :: l_z_zero    ! true if sation elev =0 when nlevs=1
   logical :: l_sim_obs   ! true if setup for simulating obs rather than errors
!
   integer, parameter, public :: bufr_unit_in=20
   integer, parameter, public :: bufr_unit_out=21
   integer, parameter :: rkind8=8  ! must be 8 since bufr interface is r8
   integer, public :: conv_nlevs   ! number of levels of data in report
   integer :: int_vtcd             ! vtcd flag converted to integer
!
   real(rkind1), parameter, public :: conv_bmiss=1.0e11
   real(rkind8), parameter :: bmiss8=1.0d11
   real(rkind8), parameter :: pscale=1.e2  ! to convert p units from hPa to Pa
   real(rkind8), parameter :: zero8=0.d0
   real(rkind8), parameter :: one8=1.d0
   real(rkind8), parameter :: ten8=10.d0
   real(rkind8), parameter :: bmiss99=bmiss8*0.99        ! allows for round-off
   real(rkind8), parameter :: code_qual_good=one8        ! default good quality
   real(rkind8), parameter :: code_qual_bad=9._rkind8    ! default bad quality
   real(rkind8), parameter :: code_program=one8          ! program code default
   real(rkind8), parameter :: code_reason_good=zero8     ! reason for good qual
   real(rkind8), parameter :: code_reason_bad=13._rkind8 ! reason for bad qual
   real(rkind8), public :: conv_info(conv_nhead)
   real(rkind8), public :: conv_values(conv_nfields,conv_max_levs)
   real(rkind8) :: hdrc(nheadc)
   real(4) :: vtcd   ! value of flag indicating Tv
!
! _OB indicates observation value, _QM indicates quality mark     
! _RC indicates "reason code"; i.e., the reason for a particular quality mark 
! _PC indicates "program code"; i.e., how obs values were computed 
! TYP is the ncep prepbufr report type #. T29 is the data dump report type #.
! The order here for driftstr,pevstr,zevstr,uevstr,tevstr,qevstr cannot
! be changed without also changeing some indexes inside the conv_rw_* routines 
   character(len=*), parameter :: driftstr='XDR YDR HRDR'
   character(len=*), parameter :: pevstr='POB PQM PPC PRC CAT'
   character(len=*), parameter :: zevstr='ZOB ZQM ZPC ZRC'
   character(len=*), parameter :: uevstr='UOB VOB WQM WPC WRC' 
   character(len=*), parameter :: tevstr='TOB TQM TPC TRC'
   character(len=*), parameter :: qevstr='QOB QQM QPC QRC' 
   character(len=*), parameter :: myname='m_bufr_conv'
!
   contains   
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_rw_setup (bufr_in_file,bufr_out_file,lprint,sim_prog, &
                             ier)
   implicit none
!
!  Open BUFR files for read and write and set some variables 
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_out_file
   character(len=*), intent(in) :: sim_prog      ! 'SIMERR' or 'SIMOBS' 
!
   integer :: n, ios
   character(len=8) :: hdrw_names(nheadw)
!
   ier=0
!
! Some functions of these subroutines are disabled if obs errors rather than 
! obs are being simulated. Their disabling is indicated by l_sim_obs=.false.
   if (sim_prog == 'SIMOBS') then ! setup called from error program
     l_sim_obs=.true.
   else
     l_sim_obs=.false.
   endif
!
! Open bufr files 
   if (trim(bufr_in_file) /= 'none') then 
     open(unit=bufr_unit_in,file=trim(bufr_in_file),form='unformatted', &
          status='old',iostat=ios)
     if (ios /= 0) then 
       ier=10
       print *,' '
       print ('(2a,i3,a,i4,2a)'),' ERROR attempting to open existing BUFR ', &
                   'file for unit=',bufr_unit_in,' iostat=',ios, &
                   ' and file name=',trim(bufr_in_file)
       return
     elseif (lprint) then 
       print *,' '
       print ('(3a,i3)'),' bufr input file=',trim(bufr_in_file),     &
              ' opened on unit=',bufr_unit_in
       call openbf(bufr_unit_in,'IN ',bufr_unit_in)
     endif
   endif
!
   if (trim(bufr_out_file) /= 'none') then
     open(unit=bufr_unit_out,file=trim(bufr_out_file),form='unformatted', &
          status='replace',iostat=ios)
     if (ios /= 0) then
       ier=20 
       print *,' '
       print ('(2a,i3,a,i4,2a)'),' ERROR attempting to open new BUFR ', &
                   'file for unit=',bufr_unit_out,' iostat=',ios, &
                   ' and file name=',trim(bufr_out_file)
       return
     elseif (lprint) then 
       print *,' '
       print ('(3a,i3)'),' bufr output file=',trim(bufr_out_file),     &
            ' opened on unit=',bufr_unit_out
     endif
     call openbf(bufr_unit_out,'OUT',bufr_unit_in)
   endif 
!
   if (ier == 0 .and. trim(bufr_in_file) /= 'none') then 
     call datelen (10)    ! formats returned idate as YYYYMMDDHH
     call ufbqcd (bufr_unit_in,'VIRTMP',vtcd) ! retrives flag value for tv
     int_vtcd=nint(vtcd) 
   endif  ! check on initialization
!
   end subroutine conv_rw_setup 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_check_obs (nerrors,ierrors,lerror)
!
!  Check if parts of obs header appears OK and if obs sybtype is in 
!  list to be considered. 
!
   use m_nr_fields_info, only : field_time_delta, field_time_first
   use m_nr_fields_info, only : field_time_slots
   use m_nr_fields_info, only : field_obs_types_num, field_obs_types
!
   implicit none
!
   logical, intent(out)   :: lerror
   integer, intent(in)    :: nerrors
   integer, intent(inout) :: ierrors(nerrors)
!
   logical :: lfound
   integer :: n 
   real(rkind8), parameter :: r180=180._rkind8
   real(rkind8), parameter :: r360=360._rkind8
   real(rkind8), parameter :: r90=90._rkind8
   real(rkind1) :: time_last
!
   lerror=.false. 
!
! Check if time range of observation is acceptable
! (only relevant if this concerns an observation simulation)
   if (l_sim_obs) then  
     if (conv_info(idhr) < field_time_first) then
       ierrors(2)=ierrors(2)+1
       lerror=.true.
     endif
     time_last=field_time_first+field_time_delta*(field_time_slots-1)
     if (conv_info(idhr) > time_last) then
       ierrors(3)=ierrors(3)+1
       lerror=.true.
     endif
   endif 
!
! Check lon 
   if (conv_info(ixob) < -r180 .or. conv_info(ixob) > r360) then
     ierrors(4)=ierrors(4)+1
     lerror=.true.
   endif
!
! Check lat
   if (abs(conv_info(iyob)) > r90) then
     ierrors(5)=ierrors(5)+1
     lerror=.true.
   endif
!
! Check nlevs
   if (conv_nlevs < 1)  then
     ierrors(6)=ierrors(6)+1
     lerror=.true.
   endif
!
! Check if obs type is in the list requested
   if (l_sim_obs) then  ! (Check only if this is an obs error simulation)  
     lfound=.false.
     do n=1,field_obs_types_num   
       if (conv_info(ityp) == field_obs_types(n)) then
         lfound=.true.
       endif 
     enddo
!
     if (.not. lfound) then 
       ierrors(7)=ierrors(7)+1
       lerror=.true.
     endif
!
   endif ! second check on l_sim_obs   
!
   end subroutine conv_check_obs
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_rw (lprint,lprint_sample,request_type,nob, &
                       lread,nerrors,ierrors,lerror,ier)
!
! Read or write one conventional observation report (may include mutiple 
! pressure levels).
!  
   logical, intent(in) :: lprint
   logical, intent(in) :: lread
   logical, intent(in) :: lprint_sample
   logical, intent(out) :: lerror
!
   integer, intent(in) :: nob
   integer, intent(in) :: nerrors
   integer, intent(out) :: ier
   integer, intent(inout) :: ierrors(nerrors)
   character(len=*), intent(in) :: request_type 
!
   integer :: lunit
   integer :: ier1
   integer :: kread
   integer :: nlevs
   integer :: ipack
   integer :: itype
   character(len=*), parameter :: mysub=myname//'::conv_rw'
!
   ier=0   ! initialize error flag to 0 (meaning no error)
!
   if (lread) then
     lunit=bufr_unit_in
     call ufbint (lunit,hdrc,nheadc,1,kread,hdstrc)
     if (kread /= 1) then 
       if (lprint) then  
         print *,' '
         print *,'ERROR DETECTED IN ',mysub
         print *,'Unable to read bufr header: kread=',kread
       endif
       ier=100
     endif
!
     if (ier == 0) then    ! can continue since header has been read
       itype=nint(hdrc(ityp))
       conv_info(:)=zero8
       nlevs=0 
!
       ier=10  ! will be replace by =0 if read is OK
       if ((request_type == 'WIND' .or. request_type == 'EITHER') .and. &
           itype > 199) then
         call conv_rw_uv (lunit,nlevs,lprint,lread,nob,ier)
       elseif ((request_type == 'MASS' .or. request_type == 'EITHER') .and. &
           itype < 200) then
         call conv_rw_tq (lunit,nlevs,lprint,lread,nob,ier)
       endif
       conv_nlevs=nlevs
!
       if (ier == 0) then
         call conv_check_obs (nerrors,ierrors,lerror)
       else
         lerror=.true.
       endif
!   
       if (.not. lerror) then ! save some flag information
         conv_info(ilzf)=0.
         if (l_z_same) conv_info(ilzf)=conv_info(ilzf)+1.
         if (l_z_zero) conv_info(ilzf)=conv_info(ilzf)+2. 
         if (l_z_present) conv_info(ilzf)=conv_info(ilzf)+4.
       else       
         ier=9
       endif
!
     endif  ! check on whether header read
!
   elseif ((.not. lread) .and. conv_nlevs > 0)  then ! attempt to write results
     lunit=bufr_unit_out
     nlevs=conv_nlevs
!     
! first, unpack some flags
     ipack=nint(conv_info(ilzf))
     if (mod(ipack,2) < 1) then
       l_z_same=.false.
     else
       l_z_same=.true.
     endif
     if (mod(ipack,4) < 2) then
       l_z_zero=.false.
     else
       l_z_zero=.true.
     endif
     if (mod(ipack,8) < 4) then
       l_z_present=.false.
     else
       l_z_present=.true.
     endif
!
! Write report to BUFR file
     itype=nint(conv_info(ityp))  
     if (itype > 199) then
       call conv_rw_uv (lunit,nlevs,lprint,lread,nob,ier1)
       ier=ier+ier1
     elseif (itype < 200) then
       call conv_rw_tq (lunit,nlevs,lprint,lread,nob,ier1)
       ier=ier+1000*ier1
     endif
!
   endif  ! test on whether read or write
!
   if (ier == 0 .and. lprint_sample) then
     call conv_print_sample (nlevs,nob,lread,' sub=conv_rw')
   endif
!
   end subroutine conv_rw 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_print_sample (nlevs,nob,lread,cinfo)
!
!  Print header and values for one observation.
!
   implicit none
!
   logical, intent(in) :: lread
   integer, intent(in) :: nlevs,nob
   character(len=*), intent(in) :: cinfo
!
   integer :: i,n
   real(rkind1) :: vprint(conv_nfields)
   character(len=8) :: station_id
!
   print *,' '
   if (lread) then
     print *,' Read observation buffer for ',cinfo
   else
     print *,' Write observation buffer ',cinfo
   endif
!
   station_id=transfer(conv_info(isid),'abcdefgh')
   print ('(a,i6,2x,a)'),'Obs number and station id:',nob,station_id
   print *,'Obs header:'
   print ('(4(a10,f12.5))'), &
           (obs_info_names(i),conv_info(i),i=1,conv_nhead)
   if (nlevs > 0) then
     print ('(a3,4x,3a8,3x,5a8,a4,7a5)'),'lev',conv_value_names(bbx),       &
       conv_value_names(bby),conv_value_names(bbr),conv_value_names(bbp),   &
       conv_value_names(bbz),conv_value_names(bbu),conv_value_names(bbv),   &
       conv_value_names(bbt),conv_value_names(bbq),conv_value_names(bbc),   &
       conv_value_names(bbf),conv_value_names(bbpq),conv_value_names(bbzq), &
       conv_value_names(bbwq),conv_value_names(bbtq),conv_value_names(bbqq)
     do i=1,nlevs
       vprint(1)=min(conv_values(bbx,i),999.999)
       vprint(2)=min(conv_values(bby,i),999.999)
       vprint(3)=min(conv_values(bbr,i),99.9999)
       vprint(4)=min(conv_values(bbp,i),999999.9)
       vprint(5)=min(conv_values(bbz,i),99999.9)
       vprint(6)=min(conv_values(bbu,i),99999.9)
       vprint(7)=min(conv_values(bbv,i),99999.9)
       vprint(8)=min(conv_values(bbt,i),99999.9)
       vprint(9)=min(conv_values(bbq,i),9.99999)
       vprint(10)=min(conv_values(bbc,i),999.)
       vprint(11)=min(conv_values(bbf,i),999.)
       vprint(12)=min(conv_values(bbpq,i),999.)
       vprint(13)=min(conv_values(bbzq,i),999.)
       vprint(14)=min(conv_values(bbwq,i),999.)
       vprint(15)=min(conv_values(bbtq,i),999.)
       vprint(16)=min(conv_values(bbqq,i),999.)
       print ('(i3,f9.3,f8.3,f8.4,f9.1,4f8.1,f8.5,7f5.0)'),i,vprint(1:16)
     enddo
   endif
!
   end subroutine conv_print_sample 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_rw_uv (lunit,nlevs,lprint,lread,nob,ier)
!
!  Read/write wind observations in NCEP .prepbufr (BUFR) format
!
   implicit none
!
   logical, intent(in)    :: lprint
   logical, intent(in)    :: lread
   integer, intent(in)    :: lunit
   integer, intent(inout) :: nlevs
   integer, intent(in)    :: nob
   integer, intent(out)   :: ier
!
! local variables
   logical :: ltype(4)
   logical :: ldrift
!
   integer, parameter :: nhead_extra=nheadw-nheadc
   integer :: said   ! obs subytpe index
   integer :: itype  ! kx (TYP) of obs
   integer :: i      ! loop counters
   integer :: nw     ! # of wind levels in a report
   integer :: nz     ! # of height values in a report
   integer :: levw, klev, levd, llev, levdch
   integer :: indw(conv_max_levs)            !  index of wind obs 
   real(rkind8) :: hdr(nheadw)               !  header info
   real(rkind8) :: obsdat(11,conv_max_levs)  !  obs profile from bufr file
   real(rkind8) :: drift(3,conv_max_levs)    !  balloon drift
   real(rkind8) :: uvout(5,conv_max_levs)    !  output wind array
   real(rkind8) :: boutp(5,conv_max_levs)    !  output for pevstr 
   real(rkind8) :: boutz(4,conv_max_levs)    !  output for zevstr 
   character(len=*), parameter :: wndstr= &
             'POB PQM UOB VOB WQM CAT ZOB ZQM PPC WPC ZPC'
!
   ier=0
!
! This first section is if a bufr read has been requested
   if (lread) then  
     nw=0            ! initialize wind level counter
     nz=0            ! initialize height level counter
     hdr(1:nheadc)=hdrc(1:nheadc)
     itype=nint(hdr(ityp))
! 
! Read extra wind information required from header
     if (nhead_extra > 0) then 
       call ufbint (lunit,hdr(nheadc+1:nheadw),nhead_extra,1,klev,&
                    hdstrw_extra)
     endif          
!
! Read obs values for all levels
     call ufbint (lunit,obsdat,11,conv_max_levs,levw,wndstr)
!
! Read drift information if RAOB, PIBAL or DROPSONDE 
     ldrift = (itype == 220) .or. (itype == 221) .or. (itype == 232) 
     if (ldrift) then
       call ufbint (lunit,drift,3,conv_max_levs,levd,driftstr)
       if (levd /= levw) then 
         if (lprint) then 
           print *,'incompatible level number returned in read of drift', &
                   nob,levd,levw
         endif
         ier=ier+1
       endif
     endif
!
! Check for levels with wind observations (and non-missing pressure)
! Also check if some observed heights are present
     do i=1,levw
       if (obsdat(1,i) < bmiss99 .and. obsdat(3,i) < bmiss99 .and.  &
           obsdat(4,i) < bmiss99) then   ! p,u,v values all present
         nw = nw + 1
         indw(nw) = i
         if (obsdat(7,i) < bmiss99) then ! z values present
           nz = nz + 1
         endif
       endif
     enddo
!
     if (nz == 0) then 
       l_z_present=.false.   ! no z values found
     else      
       l_z_present=.true.    ! some z data are present 
     endif
!
! Check if station elev set to obs. z if this is a single-level obs.
     if (l_z_present .and. nw == 1 .and. &
         abs(hdr(ielv)-obsdat(7,indw(1))) < one8) then  ! two values within 1m
       l_z_same=.true.
     else
       l_z_same=.false.
     endif 
!
! Check if station elev set to 0 if this is a surface obs.
     if (nw == 1 .and. abs(hdr(ielv)) < 0.5 .and.  &
                (itype >= 280 .or. itype <= 290) ) then
       l_z_zero=.true.
     else
       l_z_zero=.false.
     endif
!
! copy into main program arrays
     conv_values(:,:)=conv_bmiss
     if (nw > 0) then        ! found wind info, continue processing
       do i=1,nw            
         conv_values(bbu,i)=obsdat(3,indw(i))         ! u wind
         conv_values(bbv,i)=obsdat(4,indw(i))         ! v wind
         conv_values(bbp,i)=obsdat(1,indw(i))*pscale  ! p in Pa
         if (l_z_present) then
           conv_values(bbz,i)=obsdat(7,indw(i))       ! z    
         endif
! 
! Copy drift information
         if (ldrift) then
           conv_values(bbx,i)=drift(1,indw(i))  ! drifted lon
           conv_values(bby,i)=drift(2,indw(i))  ! drifted lat
           conv_values(bbr,i)=drift(3,indw(i))  ! relative time (hours)
         endif
!
! Copy CAT code
         conv_values(bbc,i)=obsdat(6,indw(i))   ! copy existing CAT code
!
! Copy quality marks
         conv_values(bbpq,i)=obsdat(2,indw(i))  ! p quality mark
         conv_values(bbwq,i)=obsdat(5,indw(i))  ! wind quality mark
         conv_values(bbzq,i)=obsdat(8,indw(i))  ! z quality mark
!
       enddo  ! loop over number of present obs   
     endif    ! test on whether any obs are present
!
! nw=0 means no data of type wndstr found
     nlevs=nw
     conv_info(1:nheadw)=hdr(1:nheadw)
!
! End of read of observation buffer
!
   elseif ((.not. lread) .and. nlevs > 0) then  
!
! This is a call to write observation buffer
!
! Write the changes to wind obs into the output buffer
! Only write report if there are valid observations in it. 
!
     hdr(1:nheadw)=conv_info(1:nheadw)
     itype=nint(hdr(ityp))
!
!  First change station elevation as necessary to agree with specifications 
!  in the NCEP .prepbufr files
     if (nlevs == 1) then
       if (l_z_zero .and. (.not. (itype == 281 .or. itype == 287))) then
         hdr(ielv)=zero8         ! preserve 0 if original sfc value was 0 
         conv_info(ielv)=zero8
         if (l_z_same) then      ! set obs lev 1 z to station elev
           conv_values(bbz,1)=zero8
         endif
       endif
       if (l_z_same) then        ! e.g., for AIRCAR or SATWND
         hdr(ielv)=conv_values(bbz,1)     
         conv_info(ielv)=conv_values(bbz,1)     
       endif
       if (itype == 285 .or. itype >= 288) then ! SCAT wind
         hdr(ielv)=ten8
         conv_info(ielv)=ten8
       endif
     endif   ! check is single-level obs
!
! Write out p event 
! (note special setting for SSMI and QKSCAT Winds in NCEP .prebufr files)
     do i=1,nlevs
       if (conv_values(bbp,i) > bmiss99) then
         boutp(1:2,i)=bmiss8
         boutp(4,i)=code_reason_bad   
       elseif (itype == 283 .or. itype == 285 .or. &
               itype >= 288 .or. itype == 282 ) then
         boutp(1,i)=1013._8        ! case for SCAT and buoy winds
         boutp(2,i)=conv_values(bbpq,i) 
         boutp(4,i)=code_reason_good 
         conv_values(bbp,i)=101300._8
       else 
         boutp(1,i)=conv_values(bbp,i)/pscale
         boutp(2,i)=conv_values(bbpq,i)      
         boutp(4,i)=code_reason_good       
       endif
       if (boutp(2,i) > bmiss99) then   ! if close, reset value to bmiss
         boutp(2,i)=bmiss8
       endif
       boutp(3,i)=code_program     
       boutp(5,i)=conv_values(bbc,i)  ! CAT value that indicates if sfc value
     enddo
     call ufbint (lunit,boutp,5,nlevs,llev,pevstr)
     if (llev /= nlevs) then 
       if (lprint) then 
         print *,'incompatible level number returned in write of pevstr ', &
                  llev,nlevs
       endif 
       ier=ier+1
     endif
!
! Write out u and v event
     do i=1,nlevs
       if (conv_values(bbu,i) > bmiss99 .or. conv_values(bbv,i) > bmiss99) then 
         uvout(1:3,i)=bmiss8
         uvout(5,i)=code_reason_bad   
       else
         uvout(1,i)=conv_values(bbu,i)   ! u wind
         uvout(2,i)=conv_values(bbv,i)   ! v wind
         uvout(3,i)=conv_values(bbwq,i)  ! wind quality mark
         uvout(5,i)=code_reason_good
       endif
       if (uvout(3,i) > bmiss99) then   ! if close, reset value to bmiss
         uvout(3,i)=bmiss8
       endif
       uvout(4,i)=code_program
     enddo
     call ufbint (lunit,uvout,5,nlevs,llev,uevstr)
     if (llev /= nlevs) then 
       if (lprint) then
         print *,'incompatible level number returned in write of uvout', &
                  llev,nlevs
       endif
       ier=ier+1
     endif
!
! write out z event 
     if (l_z_present) then  
       do i=1,nlevs
         if (conv_values(bbz,i) > bmiss99) then
           boutz(1:2,i)=bmiss8
           boutz(4,i)=code_reason_bad
         else
           boutz(1,i)=conv_values(bbz,i)
           boutz(2,i)=conv_values(bbzq,i)
           boutz(4,i)=code_reason_good
         endif  
         if (boutz(2,i) > bmiss99) then   ! if close, reset value to bmiss
           boutz(2,i)=bmiss8
         endif
         boutz(3,i)=code_program
       enddo
       call ufbint (lunit,boutz,4,nlevs,llev,zevstr)
       if (llev /= nlevs) then
         if (lprint) then 
           print *,'incompatible level number returned in write of zevstr', &
                    llev,nlevs
         endif 
         ier=ier+1
       endif
     endif   ! check on l_z_present
!
! reset drift information for radiosonde, pibal, or dropsonde
     if (itype == 220 .or. itype == 221 .or. itype == 232) then          
       do i=1,nlevs
         drift(1,i)=conv_values(bbx,i)
         drift(2,i)=conv_values(bby,i)
         drift(3,i)=conv_values(bbr,i)
       enddo
       call ufbint (lunit,drift,3,nlevs,levdch,driftstr)
       if (levdch /= nlevs) then
         if (lprint) then 
           print *,'incompatible level number returned in write of drift', &
                    levdch,nlevs
         endif
         ier=ier+1
       endif
     endif   ! test on ldrift
!
!  Finally, write out the header information.
     call ufbint(lunit,hdr,nheadw,1,klev,hdstrw)
!
!  Finished updates; so write output buffer to output file for this report
!     
     call writsb(lunit)
!
   endif  ! test on whether to read or write observation buffer
!
   end subroutine conv_rw_uv
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine conv_rw_tq (lunit,nlevs,lprint,lread,nob,ier)
!
!  Read/write t, q, ps observations in NCEP .prepbufr (BUFR) format
!
   use m_parameters, only : ratio4R 
   implicit none
!
   logical, intent(in)    :: lprint
   logical, intent(in)    :: lread
   integer, intent(in)    :: lunit
   integer, intent(inout) :: nlevs
   integer, intent(in)    :: nob
   integer, intent(out)   :: ier
!
! local variables
   logical :: ldrift
!
   integer, parameter :: nhead_extra=nheadm-nheadc
   integer, parameter :: nhead_all=nheadm+nheadw-nheadc
   integer :: itype  ! kx (TYP) of obs
   integer :: i      ! loop counters
   integer :: ntq    ! # of valid tq levels in a report
   integer :: nz     ! # of non-missing z values
   integer :: levtq, klev, levd, llev, levdch
   integer :: indm(conv_max_levs)            ! index of obs in vertical
   real(rkind8) :: obsdat(13,conv_max_levs)  ! obs profile from bufr file
   real(rkind8) :: drift(3,conv_max_levs)    ! balloon drift
   real(rkind8) :: hdr(nhead_all)            ! location and metadata
   real(rkind8) :: hdr_out(nheadm)           ! location and metadata
   real(rkind8) :: bout(4,conv_max_levs)     ! output event
   real(rkind8) :: boutp(5,conv_max_levs)    ! output event for pevstr
   real(rkind2), parameter :: c2kelvin=273.15 
   real(rkind2), parameter :: qscale=1.d-6
   character(len=*), parameter :: tqdstr= &
           'POB PQM ZOB ZQM TOB TQM QOB QQM CAT PPC ZPC TPC QPC '
!
   ier=0
!
! Check if this is a read or write
   if (lread) then  
     ntq=0            ! initialize t,q level counter
     nz=0             ! initialize height level counter
     hdr(1:nheadc)=hdrc(1:nheadc)
     itype=nint(hdr(ityp))
! 
! Read mass wind information required from header
     if (nhead_extra > 0) then 
       call ufbint (lunit,hdr(nheadw+1:nheadw+nhead_extra),nhead_extra, &
                    1,klev,hdstrm_extra)
     endif          
!
! Read obs values for all levels
     call ufbint(lunit,obsdat,13,conv_max_levs,levtq,tqdstr)
!
! Read drift information if RAOB or DROPSONDE 
     ldrift = (itype == 120) .or. (itype == 132)
     if (ldrift) then
       call ufbint(lunit,drift,3,conv_max_levs,levd,driftstr)
       if (levd /= levtq) then 
         if (lprint) then
           print *,'incompatible level number returned in read of drift', &
                    nob,levd,levtq
         endif 
         ier=ier+1
       endif
     endif
!
! Check for levels with observations (non-missing pressure and 
! one of T, q, and/or surface level).
! Also check if some observed heights are present
! If T or q is missing, make sure that quality mark is > 3.
     do i=1,levtq
       if (obsdat(1,i) < bmiss99 .and. (obsdat(5,i) < bmiss99 .or. &
               obsdat(7,i) < bmiss99 .or. abs(obsdat(9,i)) < 0.5)) then
         ntq=ntq+1
         indm(ntq)=i
         if (obsdat(3,i) < bmiss99) then
           nz=nz+1
         endif
         if (obsdat(5,i) > bmiss99) then
           obsdat(6,i)=14.
         endif
         if (obsdat(7,i) > bmiss99) then
           obsdat(8,i)=14.
         endif
       endif
     enddo
!
     if (nz == 0) then 
       l_z_present=.false. ! no valid z data found
     else      
       l_z_present=.true.  ! some z data are present 
     endif
!
! Check if station elev set to obs. z if this is a single-level obs.
     if (l_z_present .and. ntq == 1 .and. &
            abs(hdr(ielv)-obsdat(3,indm(1))) < one8) then
       l_z_same=.true.
     else
       l_z_same=.false.
     endif
!
! Check if station elev set to 0 if this is a surface obs.
     if (ntq == 1 .and. abs(hdr(ielv)) < 0.5 .and. & 
         (itype >= 180 .or. itype <= 189) ) then
       l_z_zero=.true.
     else
       l_z_zero=.false.
     endif
!
! copy into main program arrays
     if (.not. l_sim_obs) then 
       conv_values(:,:)=conv_bmiss 
     endif
!
     if (ntq > 0) then        ! found ptq info, continue processing
       do i=1,ntq 
!
         if (obsdat(5,indm(i)) < bmiss99) then           ! T not missing
           conv_values(bbt,i)=obsdat(5,indm(i))+c2kelvin ! change to Kelvin
         else
           conv_values(bbt,i)=conv_bmiss
         endif
!
         if (obsdat(7,indm(i)) < bmiss99) then          ! q not missing
           conv_values(bbq,i)=obsdat(7,indm(i))*qscale  ! change units of q 
         else
           conv_values(bbq,i)=conv_bmiss
         endif 
!
! If obs is Tv then change to T if appropriate and q value available
         if (abs(obsdat(12,indm(i))-vtcd) < 0.5 .and. & ! code = Tv
             obsdat(1,indm(i)) > 300.           .and. & ! below 300 mb
             obsdat(5,indm(i)) < bmiss99        .and. & ! T exists
             obsdat(7,indm(i)) < bmiss99        .and. & ! q exists
             obsdat(8,indm(i)) < 3.5        )   then    ! q quality mark OK
           conv_values(bbt,i)=conv_values(bbt,i)/ &
                  (one8+ratio4R*conv_values(bbq,i)) 
         endif
!
         conv_values(bbp,i)=obsdat(1,indm(i))*pscale ! change p from hPa to Pa
!
         if (l_z_present) then
           conv_values(bbz,i)=obsdat(3,indm(i))      ! z obs value
         else
           conv_values(bbz,i)=conv_bmiss
         endif
! 
! Copy drift information
         if (ldrift) then
           conv_values(bbx,i)=drift(1,indm(i))  ! drifted lon
           conv_values(bby,i)=drift(2,indm(i))  ! drifted lat
           conv_values(bbr,i)=drift(3,indm(i))  ! relative time (hours)
         else 
           conv_values(bbx,i)=conv_bmiss
           conv_values(bby,i)=conv_bmiss
           conv_values(bbr,i)=conv_bmiss
         endif
!
! Copy CAT and TV codes
         conv_values(bbc,i)=obsdat(9,indm(i))     ! CAT
         if (.not. l_sim_obs) then 
           conv_values(bbf,i)=obsdat(12,indm(i))  ! TV flag
         endif
!
! Copy quality marks
         conv_values(bbpq,i)=obsdat(2,indm(i))   ! p quality mark
         conv_values(bbzq,i)=obsdat(4,indm(i))   ! z quality mark
         conv_values(bbtq,i)=obsdat(6,indm(i))   ! T quality mark
         conv_values(bbqq,i)=obsdat(8,indm(i))   ! q quality mark
!
       enddo   ! loop over number of present obs    
     endif     ! test on whether any obs are present
!
! ntq=0 means no T,Q,Z data found
     nlevs=ntq
     if (nheadw > nheadc) then
       hdr(nheadc+1:nheadw)=bmiss8
     endif 
     conv_info(1:nhead_all)=hdr(1:nhead_all)   
!
! End of read of observation buffer
!
   elseif ((.not. lread) .and. nlevs > 0) then
!
! This is a call to write observation buffer
! Only write report if there are valid observations in it. 
!
     hdr(1:nhead_all)=conv_info(1:nhead_all)
     itype=nint(hdr(ityp))
!
! Change surface elevation to NR elevation in line of report
     do i=1,nlevs
       if (nint(conv_values(bbc,i)) == 0) then  ! CAT code for obs
         conv_values(bbz,i)=conv_info(ielv)     ! set lev z to surface height
       endif                
     enddo
!
! For single-level reports: set to station elev if surface report, unless 
! original value of station elev was 0, in which case the 0 value is preserved.
! Also, if station elev and obs z were originally the same, set the former to
! the latter computed from the nature run, unless both were originally 0.
!
     if (nlevs == 1) then
       if (itype >= 181 .or. itype <= 189 ) then
         if (l_z_zero) then
           hdr(ielv)=zero8           ! preserve 0 if original value was 0 
           conv_info(ielv)=zero8
         else 
           hdr(ielv)=conv_info(ielv) ! reset to sfc z in nature run
         endif
       endif
       if (l_z_same) then  
         if (l_z_zero) then
           conv_values(bbz,1)=zero8          ! hdr(ilev) already set to 0 above.
         else
           hdr(ielv)=conv_values(bbz,1)      ! e.g., for AIRCAR
           conv_info(ielv)=conv_values(bbz,1)
         endif
       endif
     endif   ! test on whether obs is single-level
!
! Write out p event 
     do i=1,nlevs
       if (conv_values(bbp,i) > bmiss99) then
         boutp(1:2,i)=bmiss8
         boutp(4,i)=code_reason_bad   
       else 
         boutp(1,i)=conv_values(bbp,i)/pscale
         boutp(2,i)=conv_values(bbpq,i)      
         boutp(4,i)=code_reason_good       
       endif
       if (boutp(2,i) > bmiss99) then   ! if close, reset value to bmiss
         boutp(2,i)=bmiss8
       endif
       boutp(3,i)=code_program     
       boutp(5,i)=conv_values(bbc,i)   ! CAT value that indicates if sfc value
     enddo
     call ufbint (lunit,boutp,5,nlevs,llev,pevstr)
     if (llev /= nlevs) then 
       if (lprint) then
         print *,'incompatible level number returned in write of pevstr ', &
                  llev,nlevs
       endif 
       ier=ier+1
     endif
!
! write out T event (either T or Tv as requested)
     do i=1,nlevs
       if (conv_values(bbt,i) > bmiss99) then 
         bout(1:3,i)=bmiss8
         bout(4,i)=code_reason_bad
       elseif (nint(conv_values(bbf,i)) /= int_vtcd .or. & ! code /= Tv
                    conv_values(bbp,i) < 30000.     .or. & ! obs above 300 mb
                    conv_values(bbq,i) > bmiss99    .or. & ! new q missing
                    conv_values(bbqq,i) > 3.5 )     then   ! old q value bad 
         bout(1,i)=conv_values(bbt,i)-c2kelvin    ! write as T
         bout(2,i)=conv_values(bbtq,i)            ! copy quality mark
         bout(4,i)=zero8                          ! reason code
         if (nint(conv_values(bbf,i)) == int_vtcd) then
           bout(3,i)=one8                         ! replace code /= Tv
         else 
           bout(3,i)=conv_values(bbf,i)           ! copy program code
         endif
       else                                       ! write as Tv
         bout(1,i)=conv_values(bbt,i)* &
                    (one8+ratio4R*conv_values(bbq,i))-c2kelvin
         bout(2,i)=conv_values(bbtq,i)            ! copy quality mark
         bout(3,i)=real(int_vtcd,8)               ! set code = Tv 
         bout(4,i)=zero8                          ! reason code
       endif
       if (bout(2,i) > bmiss99) then   ! if close, reset value to bmiss
         bout(2,i)=bmiss8
       endif
     enddo
     call ufbint(lunit,bout,4,nlevs,llev,tevstr)
     if (llev /= nlevs) then 
       if (lprint) then
         print *,'incompatible level number returned in write of tevstr ', &
                  llev,nlevs
       endif
       ier=ier+1
     endif
!
! write out q event
     do i=1,nlevs
       if (conv_values(bbq,i) > bmiss99) then  
         bout(1:3,i)=bmiss8
         bout(4,i)=code_reason_bad
       else
         bout(1,i)=conv_values(bbq,i)/qscale  ! change scale of q
         bout(2,i)=conv_values(bbqq,i)        ! copy quality mark
         bout(3,i)=code_program               ! program code
         bout(4,i)=code_reason_good           ! reason code
       endif
       if (bout(2,i) > bmiss99) then   ! if close, reset value to bmiss
         bout(2,i)=bmiss8
       endif
     enddo
     call ufbint (lunit,bout,4,nlevs,llev,qevstr)
     if (llev /= nlevs) then
       if (lprint) then 
         print *,'incompatible level number returned in write of qevstr ', &
                  llev,nlevs
       endif 
       ier=ier+1
     endif
!
! write out z event if either surface value or z values present 
     do i=1,nlevs
       if (abs(conv_values(bbc,i)) < 0.5)  then  ! this is a surface value
         bout(1,i)=conv_values(bbz,i)            ! new (NR) surface value
         bout(2,i)=conv_values(bbzq,i)           ! copy original quality marker
         bout(4,i)=code_reason_good              ! reason code
       elseif (l_z_present .and. conv_values(bbz,i) < bmiss99) then
         bout(1,i)=conv_values(bbz,i)    ! computed z
         bout(2,i)=conv_values(bbzq,i)   ! copy original quality marker
         bout(4,i)=code_reason_good      ! reason code
       else                              ! observation bad
         bout(1,i)=bmiss8        
         bout(2,i)=bmiss8        
         bout(4,i)=code_reason_bad       ! reason code
       endif
       if (bout(2,i) > bmiss99) then     ! if close, reset value to bmiss
         bout(2,i)=bmiss8
       endif
       bout(3,i)=code_program            ! copy program code
     enddo
     call ufbint (lunit,bout,4,nlevs,llev,zevstr)
     if (llev .ne. nlevs) then
       if (lprint) then 
         print *,'incompatible level number returned in write of zevstr ',&
                  llev,nlevs
       endif
       ier=ier+1
     endif
!
! reset drift information for radiosonde or dropsonde
     if (itype == 120 .or. itype == 132) then          
       do i=1,nlevs
         drift(1,i)=conv_values(bbx,i)
         drift(2,i)=conv_values(bby,i)
         drift(3,i)=conv_values(bbr,i)
       enddo
       call ufbint (lunit,drift,3,nlevs,levdch,driftstr)
       if (levdch /= nlevs) then
         if (lprint) then 
           print *,'incompatible level number returned in write of drift', &
                    levdch,nlevs
         endif
         ier=ier+1
       endif
     endif   ! test on ldrift

!  Finally, write out the header information
     hdr_out(1:nheadc)=hdr(1:nheadc)
     if (nhead_extra > 0) then
       hdr_out(nheadc+1:nheadc+nhead_extra)=hdr(nheadw+1:nheadw+nhead_extra)
     endif
     call ufbint(lunit,hdr_out,nheadm,1,klev,hdstrm)
!
!  finished updates; so write output buffer to output file for this report
!     
     call writsb(lunit)
!
   endif  ! test on write and nlevs>0
!
   end subroutine conv_rw_tq
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_bufr_conv
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   logical function check_type_conv (subset,type,sim_prog)
   character(len=*) subset
   character(len=*) type
   character(len=*), intent(in) :: sim_prog      ! 'SIMERR' or 'SIMOBS' 

! Code for making data selection of subtypes given the data type requested:
! subset is the data subset type name read from the BUFR files
! type is the data type specified in the main program argument list.
!
! NOTE: if type is indicated as ALL, then both mass and wind observations are
! included, except for SATWND if obs rather than obs errors are to be simulated.
!
   logical :: check_type
   character(len=8), parameter :: temp_types(7) = (/ 'ADPUPA', &
       'RASSDA', 'SFCSHP', 'ADPSFC', 'MSONET', 'AIRCAR', 'AIRCFT' /)
!
   character(len=8), parameter :: wind_types(15) = (/ 'SYNDAT', &
      'ADPUPA','PROFLR','VADWND','AIRCFT','AIRCAR','SATWND',    &
      'SFCSHP','SPSSMI','ADPSFC','QKSWND','ERS1DA','MSONET',    &
      'ASCATW','WDSATR' /)
!
   select case ( trim(type) )
   case ('MASS')
     check_type = any( subset == temp_types )
   case ('WIND')
     check_type = any( subset == wind_types )
   case ('ALL')
     check_type = .true.
     if (sim_prog == 'SIMOBS' .and. subset == 'SATWND') then 
       check_type=.false.
     endif 
   end select
   check_type_conv=check_type
   return
   end function check_type_conv 
!

