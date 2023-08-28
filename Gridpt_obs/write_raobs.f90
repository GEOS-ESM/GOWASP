!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine write_obs_tq (obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4, &
                               obs_nlevs,obs_nfields,nob,nobs,luout,    &
                               ltest,obs_info,obs_levels,obs_values)
!
!  Write temperature, moisture, or pressure obs 
!  in NCEP .prepbufr (BUFR) format
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
      use m_kinds, only : rkind1
      use m_parameters
      implicit none
!
      logical, intent(in) :: ltest
!
      integer, intent(in) :: obs_ndim1 
      integer, intent(in) :: obs_ndim2 
      integer, intent(in) :: obs_ndim3 
      integer, intent(in) :: obs_ndim4 
      integer, intent(in) :: obs_nlevs
      integer, intent(in) :: obs_nfields
      integer, intent(in) :: nob,nobs
      integer, intent(in) :: luout
!
      real(rkind1), intent(in) :: obs_info(obs_ndim4)
      real(rkind1), intent(in) :: obs_levels(obs_ndim1,obs_ndim3)
      real(rkind1), intent(in) :: obs_values(obs_ndim1,obs_ndim2)
!
! local variables
      integer, parameter :: nhead=7       !  size of header
      integer :: i                        !  loop counters
      integer :: levtq, klev, levd, llev, levdch
      integer :: nobsx
!
      real(8), parameter :: bmiss8=10.d10 ! indicates missing value
      real(8), parameter :: c2kelvin=273.15 
      real(8), parameter :: pscale=1.e2
      real(8), parameter :: qscale=1.e-6
      real(4) :: vtcd ! value of flag that indicated Tv rather than T
      real(8) :: drift(3,255)       !  balloon drift
      real(8) :: hdr(nhead)         !  location and metadata
      real(8) :: bout(4,255)        !  output event
      real(8) :: boutp(5,255)       !  output event for pevstr
!
      character(len=8)   :: station_id
      character(len=52)  :: hdstr, driftstr
      character(len=32)  :: pevstr, tevstr, qevstr, zevstr
      data hdstr     /'SID XOB YOB DHR TYP ELV T29' /
      data driftstr  /'XDR YDR HRDR               ' /
      data pevstr    /'POB  PQM  PPC PRC CAT ' /
      data tevstr    /'TOB  TQM  TPC TRC     ' /
      data qevstr    /'QOB  QQM  QPC QRC     ' /
      data zevstr    /'ZOB  ZQM  ZPC ZRC     ' /
!
! _OB indicates observation value, _QM indicates quality mark     
! _RC indicates "reason code"; i.e., the reason for a particular quality mark 
! _PC indicates "program code"; i.e., how obs values were computed 
! Here, the _QM are copied from the original since (1) the quality-accepted obs
! count for the simulated data will then be similar to the original data
! and (2) any "missing data" created by the simulation program (e.g., 
! below ground values in a raob report) will still be detectable as a 
! missing value flag.
! The _PC values are copied from the original in case the DAS expects some
! particular values (in the GSI used at time of this development, _PC values
! are not used, but this can change...)
!
      call ufbqcd(luout,'VIRTMP',vtcd)
!
      write (station_id,'(i8)') nint(obs_info(1))  ! expess id as character
      hdr(1)=transfer(station_id,hdr(1))           ! reformat as type hdr
      hdr(2)=obs_info(2)       ! lon of obs
      hdr(3)=obs_info(3)       ! lat of obs
      hdr(4)=obs_info(4)       ! time of obs
      hdr(5)=obs_info(5)       ! type of obs (ncep kx #)
      hdr(6)=obs_info(6)       ! height of surface 
      hdr(7)=obs_info(7)       ! type of obs 
!
! write out p event
      do i = 1,obs_nlevs
        boutp(1,i) = obs_levels(i,1)/pscale
        boutp(2,i) = 1.                  ! quality mark
        boutp(3,i) = 1.                  ! program code
        boutp(4,i) = 0.                  ! reason code
        boutp(5,i) = 1.                  ! CAT code (not surface)
      enddo 
      boutp(5,1) = 0.                    ! CAT code (surface)
      call ufbint(luout,boutp,5,obs_nlevs,llev,pevstr)
!
! write out t event
      do i = 1,obs_nlevs
        if ( (obs_levels(i,1) < 30000.)     .or. & ! above 300 m
             (obs_values(i,2)==bmiss8) ) then ! new q missing
          bout(1,i) = obs_values(i,1)-c2kelvin     ! write as T
          bout(3,i) = 1.                           ! must be not vtcd
        else                                       ! write as Tv
          bout(1,i) = obs_values(i,1) * &
                      (1. + ratio4R*obs_values(i,2)) - c2kelvin
          bout(3,i) = vtcd                  ! program code indicating Tv
        endif
        bout(2,i) = 1.                      ! quality mark
        bout(4,i) = 0.                      ! reason code
      enddo
      call ufbint(luout,bout,4,obs_nlevs,llev,tevstr)
!
! write out q event
      do i = 1,obs_nlevs
        bout(1,i) = obs_values(i,2)/qscale  ! change scale of q
        bout(2,i) = 1.                      ! quality mark
        bout(3,i) = 0.                      ! program code
        bout(4,i) = 0.                      ! reason code
      enddo
      call ufbint(luout,bout,4,obs_nlevs,llev,qevstr)
!
! write out z event 
      do i = 1,obs_nlevs
        bout(1,i) = obs_levels(i,2)      ! computed z
        bout(2,i) = 1.                   ! quality marker
        bout(3,i) = 1.                   ! program code
        bout(4,i) = 0.                   ! reason code
      enddo
      call ufbint(luout,bout,4,obs_nlevs,llev,zevstr)
!
! Set drift information for radiosonde and pibal obs (120, 220, 221)     
! as though there is no drift
      if (hdr(5)==120..or.hdr(5)==220..or.hdr(5)==221.) then          
        do i = 1,obs_nlevs
          drift(1,i)=hdr(2)
          drift(2,i)=hdr(3)
          drift(3,i)=hdr(4)
        enddo
        call ufbint(luout,drift,3,obs_nlevs,levdch,driftstr)
      endif   ! test on obs_type hdr(5)

!  Finally, write out the header information
      call ufbint(luout,hdr,nhead,1,klev,hdstr)
!
!  finished updates; so write output buffer to output file for this report
      call writsb(luout)
!
! Optionally print some info if testing mode
      nobsx=nobs/5 ! will yield 5 samples for output
      if (ltest .and. mod(nob,nobsx) == 1) then 
        print *,' '
        print *,' Write observation buffer'
        print ('(a,3i7,a8,6f10.2)'),' OBS:',nob,obs_nlevs,obs_nfields, &
                                    station_id,hdr(2:nhead)
        if (obs_nlevs > 0) then 
          do i=1,obs_nlevs
            print ('(i4,2f10.2,15f12.5)'),i,obs_levels(i,1:2), &
                                obs_values(i,1:obs_nfields)
          enddo
        endif
      endif  ! test on ltest
!
      end subroutine write_obs_tq
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine write_obs_uv (obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4, &
                               obs_nlevs,obs_nfields,nob,nobs,luout,    &
                               ltest,obs_info,obs_levels,obs_values)
!
!  Read/write wind observations in NCEP .prepbufr (BUFR) format
!
      use m_kinds, only : rkind1
      implicit none
!
      logical, intent(in) :: ltest
!
      integer, intent(in) :: obs_ndim1 
      integer, intent(in) :: obs_ndim2 
      integer, intent(in) :: obs_ndim3 
      integer, intent(in) :: obs_ndim4 
      integer, intent(in) :: obs_nlevs
      integer, intent(in) :: obs_nfields
      integer, intent(in) :: nob,nobs
      integer, intent(in) :: luout
!
      real(rkind1), intent(in) :: obs_info(obs_ndim4)
      real(rkind1), intent(in) :: obs_levels(obs_ndim1,obs_ndim3)
      real(rkind1), intent(in) :: obs_values(obs_ndim1,obs_ndim2)
!
! local variables
      logical :: l_z_same ! true if station elev=z of 1-level observation
      logical :: l_z_zero ! true if sation elev =0 for 1-level obs.
!
      integer, parameter :: nhead=7      !  size of header
      integer :: i                ! loop counters
      integer :: nw               ! # of wind levels in a report
      integer :: nz               ! # of height values in a report
      integer :: levw, klev, levd, llev, levdch
      integer :: nobsx
!
      real(8), parameter :: bmiss8=10.d10 ! indicates missing value
      real(8), parameter :: pscale=1.e2
      real(8) :: drift(3,255)      !  balloon drift
      real(8) :: hdr(nhead)        !  location and metadata
      real(8) :: uvout(5,255)      !  output wind array
      real(8) :: boutp(5,255)      !  output for pevstr 
      real(8) :: boutz(4,255)      !  output for zevstr 
!
      character(len=8)  :: station_id
      character(len=44) :: hdstr, driftstr, uevstr, pevstr, zevstr
      data hdstr     /'SID XOB YOB DHR TYP ELV T29                 ' /
      data driftstr  /'XDR YDR HRDR                                ' /
      data pevstr    /'POB PQM PPC PRC CAT                         ' /
      data uevstr    /'UOB VOB WQM WPC WRC                         ' /
      data zevstr    /'ZOB ZQM ZPC ZRC                             ' /
!
! _OB indicates observation value, _QM indicates quality mark     
! _RC indicates "reason code"; i.e., the reason for a particular quality mark 
! _PC indicates "program code"; i.e., how obs values were computed 
! Here, the _QM are copied from the original since (1) the quality-accepted obs
! count for the simulated data will then be similar to the original data
! and (2) any "missing data" created by the simulation program (e.g., 
! below ground values in a raob report) will still be detectable as a 
! missing value flag.
! The _PC values are copied from the original in case the DAS expects some
! particular values (in the GSI used at time of this development, _PC values
! are not used, but this can change...)
!
      write (station_id,'(i8)') nint(obs_info(1))  ! expess id as character
      hdr(1)=transfer(station_id,hdr(1))           ! reformat as type hdr
      hdr(2)=obs_info(2)       ! lon of obs
      hdr(3)=obs_info(3)       ! lat of obs
      hdr(4)=obs_info(4)       ! time of obs
      hdr(5)=obs_info(5)       ! type of obs (ncep kx #)
      hdr(6)=obs_info(6)       ! height of surface 
      hdr(7)=obs_info(7)       ! type of obs
!
! write out p event
      do i = 1,obs_nlevs
        boutp(1,i) = obs_levels(i,1)/pscale
        boutp(2,i) = 1.                  ! copy quality mark
        boutp(3,i) = 1.                  ! copy program code
        boutp(4,i) = 0.                  ! reason code
        boutp(5,i) = 1.                  ! CAT code (not surface)
      enddo 
      boutp(5,1) = 0.                    ! CAT code (surface)
      call ufbint(luout,boutp,5,obs_nlevs,llev,pevstr)
!
! Write out u and v event
      do i = 1,obs_nlevs
        uvout(1,i) = obs_values(i,1)
        uvout(2,i) = obs_values(i,2)
        uvout(3,i) = 1.                 ! copy quality mark
        uvout(4,i) = 0.                 ! program code
        uvout(5,i) = 0.                 ! reason code
      enddo
      call ufbint(luout,uvout,5,obs_nlevs,llev,uevstr)
!
! write out z event 
      do i = 1,obs_nlevs
        boutz(1,i) = obs_levels(i,2)      ! computed z
        boutz(2,i) = 1.                   ! quality marker
        boutz(3,i) = 1.                   ! program code
        boutz(4,i) = 0.                   ! reason code
      enddo
      call ufbint(luout,boutz,4,obs_nlevs,llev,zevstr)
!
! Set drift information for radiosonde and pibal obs (120, 220, 221)     
! as though there is no drift
      if (hdr(5)==120..or.hdr(5)==220..or.hdr(5)==221.) then          
        do i = 1,obs_nlevs
          drift(1,i)=hdr(2)
          drift(2,i)=hdr(3)
          drift(3,i)=hdr(4)
        enddo
        call ufbint(luout,drift,3,obs_nlevs,levdch,driftstr)
      endif   ! test on obs_type hdr(5)

!  Finally, write out the header information
      call ufbint(luout,hdr,nhead,1,klev,hdstr)
!
!  finished updates; so write output buffer to output file for this report
!     
      call writsb(luout)
!
! Optionally print some info if testing mode
      nobsx=nobs/5 ! will yield 5 samples for output
      if (ltest .and. mod(nob,nobsx) == 1) then 
        print *,' '
        print *,' Write observation buffer'
        print ('(a,3i7,a8,6f10.2)'),' OBS:',nob,obs_nlevs,obs_nfields, &
                                    station_id,hdr(2:nhead)
        if (obs_nlevs > 0) then 
          do i=1,obs_nlevs
            print ('(i4,2f10.2,15f12.5)'),i,obs_levels(i,1:2), &
                                obs_values(i,1:obs_nfields)
          enddo
        endif
      endif  ! test on ltest
!
      end subroutine write_obs_uv
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      logical function check_type_raob(subset,type)
      character(len=*) subset
      character(len=*) type

! Code for making data selection of subtypes given the data type requested:
! subset is the data subset type name read from the BUFR files
! type is the data type specified in the main program argument list
!
      character(len=8), parameter :: types(1) = (/ 'ADPUPA' /)
!
      check_type_raob = (types(1) == subset)
!
      end
!
!
