!
   Module m_bufr_gpsro
!
! Module to read and write GPSRO observation data in bufr format.
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1,rkind2
!
   use m_gpsro_names, only : obs_info_num, obs_nvalues, obs_max_levs
   use m_gpsro_names, only : bbnum, bbnlold, bbnlnew, bbnmsg
   use m_gpsro_names, only : bbsiid, bbbaz, bbyear, bblat, bblon
   use m_gpsro_names, only : bbimpp, bbbang, bbdhr
   use m_gpsro_names, only : bblatk, bblonk, bbbazk
   use m_gpsro_names, only : obs_nhead, obs_hdstr1, obs_hdstr2
!
   implicit none
!
   private
   public :: gpsro_rw_setup
   public :: gpsro_rw
   public :: gpsro_check_obs
!
   logical :: l_sim_obs   ! true if setup for simulating obs rather than errors
!
   integer, parameter, public :: bufr_unit_in=20
   integer, parameter, public :: bufr_unit_out=21
   integer, parameter :: rkind8=8  ! must be 8 since bufr interface is r8
   integer, public :: obs_nlevs   ! number of levels of data in report
!
   real(rkind1), parameter, public :: obs_bmiss=10.e10
   real(rkind8), parameter :: bmiss8=10.d10
   real(rkind8), parameter :: pscale=1.e2  ! to convert p units from hPa to Pa
   real(rkind8), parameter :: zero8=0.d0
   real(rkind8), parameter :: one8=1.d0
   real(rkind8), parameter :: bmiss99=bmiss8*0.99    ! allows for round-off
   real(rkind8), public :: gpsro_info(obs_info_num)
   real(rkind8), public :: gpsro_values(obs_nvalues,obs_max_levs)
   character(len=*), parameter :: myname='m_bufr_gpsro'
!
   contains   
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine gpsro_rw_setup (bufr_in_file,bufr_out_file,lprint, &
                              sim_prog,ier)
!
!  Open BUFR file for GPSRO obs and set some variables
!
   implicit none
!
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_out_file
   character(len=*), intent(in) :: sim_prog      ! 'SIMERR' or 'SIMOBS' 
!
   integer :: ios
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
     call maxout(200000)            ! increase size of ouput bufr
   endif 
!
   if (ier == 0 .and. trim(bufr_in_file) /= 'none') then 
     call datelen (10)    ! formats returned idate as YYYYMMDDHH
   endif  ! check on initialization
!
   end subroutine gpsro_rw_setup 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine gpsro_check_obs (nerrors,ierrors,lerror)
!
!  Check that some obs header parameters are in proper ranges
!
   use m_nr_fields_info, only : field_time_delta, field_time_first
   use m_nr_fields_info, only : field_time_slots
!
   implicit none
!
   logical, intent(out)   :: lerror
   integer, intent(in)    :: nerrors
   integer, intent(inout) :: ierrors(nerrors)
!
   real(rkind8), parameter :: r180=180._rkind8
   real(rkind8), parameter :: r360=360._rkind8
   real(rkind8), parameter :: r90=90._rkind8
   real(rkind1) :: time_last
!
   lerror=.false. 
!
! Check if time range of observation is acceptable
! (only relevant if this concerns an observation simulation)
   if (l_sim_obs) then  ! 
     if (gpsro_info(bbdhr) < field_time_first) then  
       ierrors(2)=ierrors(2)+1
       lerror=.true.
     endif
     time_last=field_time_first+field_time_delta*(field_time_slots-1)
     if (gpsro_info(bbdhr) > time_last) then
       ierrors(3)=ierrors(3)+1
       lerror=.true.
     endif
   endif
!
! Check lon 
   if (gpsro_info(bblon) < -r180 .or. gpsro_info(bblon) > r360) then
     ierrors(4)=ierrors(4)+1
     lerror=.true.
   endif
!
! Check lat
   if (abs(gpsro_info(bblat)) > r90) then
     ierrors(5)=ierrors(5)+1
     lerror=.true.
   endif
!
! Check nlevs
   if (obs_nlevs < 1)  then
     ierrors(6)=ierrors(6)+1
     lerror=.true.
   endif
!
   end subroutine gpsro_check_obs
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine gpsro_rw (lprint,lread,cdtime,nob,nerrors,ierrors,lerror,ier)
! 
! Call routine to read/write and check obs for one GPSRO for one location 
! 
   use m_time_compute, only : time_compute_unpack
   use m_time_compute, only : time_compute_dhours
   use m_time_compute, only : time_compute_add
   use m_time_compute, only : rkindh
!
   implicit none
!
   logical, intent(in) :: lprint
   logical, intent(in) :: lread
   logical, intent(out) :: lerror
   integer, intent(in) :: nob
   integer, intent(in) :: nerrors
   integer, intent(out) :: ier
   integer, intent(inout) :: ierrors(nerrors)
   character(len=*), intent(in) :: cdtime
!
   real(rkindh) :: dhours
   real(rkind1) :: tref(6)
   real(rkind1) :: tobs(6)
   character(len=*), parameter :: mysub=myname//'::gpsro_rw'
!
   ier=0   ! initialize error flag to 0 (meaning no error)
!
   if (lread) then
     call gpsro_rw_bang (lprint,lread,nob,ier)
     if (ier == 0) then    ! can continue since header has been read
!
! Compute obs time - ref time
       if (l_sim_obs) then  
         tobs(1:6)=gpsro_info(bbyear:bbyear+5)
         call time_compute_unpack (cdtime,tref)   ! old reference time
         call time_compute_dhours (tref,tobs,dhours,ier)
         gpsro_info(bbdhr)=dhours
       endif  
       call gpsro_check_obs (nerrors,ierrors,lerror)
     endif    ! check on whether obs read properly
!
   elseif ((.not. lread) .and. obs_nlevs > 0)  then 
!
! Write report to BUFR file
     if (l_sim_obs) then  
       dhours=gpsro_info(bbdhr)     
       call time_compute_unpack (cdtime,tref)      ! new reference time
       call time_compute_add (dhours,tref,tobs,ier)  
       gpsro_info(bbyear:bbyear+5)=tobs(1:6)       ! new obs time
     endif
     call gpsro_rw_bang (lprint,lread,nob,ier)
!
   endif  ! test on whether read or write
!
   end subroutine gpsro_rw 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine gpsro_rw_bang (lprint,lread,nob,ier)
!
!  Use calls to BUFR library to read file of GPSRO bending angle obs
!
      implicit none
!
      logical, intent(in)    :: lprint
      logical, intent(in)    :: lread
      integer, intent(in)    :: nob
      integer, intent(out)   :: ier
!
! local variables
      integer, parameter :: ndval=25     ! max number of data values per level
      integer, parameter :: nval_out=10  ! number of obs values output
      integer :: lunit      
      integer :: khead1,khead2,klevs1,klevs2
      integer :: i, j, j6
      integer :: roseq2(obs_max_levs)
      real(rkind8) :: hdr(obs_nhead)               ! instrument header information
      real(rkind8) :: obsdat(ndval,obs_max_levs)   ! obs profile from bufr file
      real(rkind8) :: obsdat_out(nval_out,obs_max_levs)  ! obs profile from bufr file
      real(rkind8) :: bout(ndval,obs_max_levs)     ! output event 
      real(rkind8) :: nfreqs(obs_max_levs)         ! number of frequencies for imapct lev
      character(len=*),parameter :: mysub=myname//"::gpsro_rw"
!
      ier=0
!
! Check if this is a read or write
      if (lread) then
        lunit=bufr_unit_in
!
! Read header information for 1 report
        j=obs_nhead-8
        call ufbint(lunit,hdr(1:8),8,1,khead1,obs_hdstr1)
        call ufbint(lunit,hdr(9:obs_nhead),j,1,khead2,obs_hdstr2)
        if (khead1 /= 1 .or. khead2 /= 1) then 
          if (lprint) then  
            print *,' '
            print *,'ERROR DETECTED IN ',mysub
            print *,'Unable to read bufr header: khead1,2=',khead1,khead2
          endif
          ier=1
        endif
!
! Read data for all impact levels
        call ufbint(lunit,nfreqs,1,obs_max_levs,klevs1,'{ROSEQ2}')
        call ufbseq(lunit,obsdat,ndval,obs_max_levs,klevs2,'ROSEQ1') 
        if (klevs1 /= klevs2) then
          if (lprint) then  
            print *,' '
            print *,'ERROR DETECTED IN ',mysub
            print *,'klevs1 /= klevs2 :',klevs1,klevs2
          endif
          ier=ier+10
          klevs1=0
        endif 
!
        obs_nlevs=0
        do i=1,klevs1
!
! The following are alt, lon, bearing that vary with impact parameter
          gpsro_values(bblatk,i)=obsdat(1,i)
          gpsro_values(bblonk,i)=obsdat(2,i)
          gpsro_values(bbbazk,i)=obsdat(3,i)
          do j=1,int(nfreqs(i))      ! loop over frequencies present 
            j6=6*j-2
            if (nint(obsdat(j6,i))==0) then          ! only consider freq 0 
              gpsro_values(bbimpp,i)=obsdat(j6+1,i)  ! impact parameter (IMPP)
              gpsro_values(bbbang,i)=obsdat(j6+2,i)  ! bending angle (BNDA)
              obs_nlevs=obs_nlevs+1
            endif 
          enddo
        enddo
        gpsro_info(1:obs_nhead)=hdr(1:obs_nhead)    
!
! End of read of observation buffer
!
      elseif (obs_nlevs > 0) then
        lunit=bufr_unit_out
!
! This is a call to write observation buffer, and observation values exist
!
! Write out the header information
        hdr(1:obs_nhead)=gpsro_info(1:obs_nhead)
        j=obs_nhead-8
        call ufbint(lunit,hdr(1:8),8,1,khead1,obs_hdstr1)
        call ufbint(lunit,hdr(9:obs_nhead),j,1,khead2,obs_hdstr2)
!
        call drfini(lunit,obs_nlevs,1,'(ROSEQ1)') 
        roseq2(:)=1
        call drfini(lunit,roseq2,1,'{ROSEQ2}') ! only one frequency simulated
!
! Write out data sequence
        do i=1,obs_nlevs
          obsdat_out(1,i)=gpsro_values(bblatk,i) ! lat (vary with imp param) 
          obsdat_out(2,i)=gpsro_values(bblonk,i) ! lon (vary with imp param) 
          obsdat_out(3,i)=gpsro_values(bbbazk,i) ! azimuth (vary with imp param)
          obsdat_out(4,i)=0.0_8                  ! set frequency 0
          obsdat_out(5,i)=gpsro_values(bbimpp,i) ! impact parameter (IMPP)
          obsdat_out(6,i)=gpsro_values(bbbang,i) ! bending angle (BNDA)
          obsdat_out(7,i)=13.0_8                 ! FOST 
          obsdat_out(8,i)=gpsro_values(bbbang,i) ! bending angle (BNDA)
          obsdat_out(9,i)=13.0_8                 ! FOST 
          obsdat_out(10,i)=hdr(15)               ! PCCF
          do j=1,nval_out
            if (obsdat_out(j,i) > bmiss99) then
              obsdat_out(j,i)=bmiss8 
            endif
          enddo
        enddo
        call ufbseq(lunit,obsdat_out,nval_out,obs_nlevs,klevs2,'ROSEQ1') 
!
! Write output buffer to output file for this report
!     
        call writsb(lunit)
!
      endif   ! test on whether to read or write observation buffer
!
      end subroutine gpsro_rw_bang
!
!
      end module m_bufr_gpsro
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      logical function check_type_gpsro (subset,type)
      character(len=*) subset
      character(len=*) type

! Code for making data selection of subtypes given the data type requested:
! subset is the data subset type name read from the BUFR files
! type is the data type specified in the main program argument list
!
      logical :: check_type
      character(len=8), parameter :: gpsro_types(1) = 'NC003010'
!
      select case ( trim(type) )
      case ('GPSRO')
        check_type = any( subset == gpsro_types )
      end select
      check_type_gpsro=check_type
      return
      end function check_type_gpsro 
!

 
