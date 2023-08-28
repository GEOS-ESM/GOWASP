!
  Module m_bufr_rad
!
! Module to read and write radiance observation data in bufr format.
! Also includes routine to read generic radiance text files.
!
! This is explicitly written to be compatible with datasets used at the GMAO
! for input into GSI (July 2014 versions)
!
! Initial Code provided by Meta Sienkiewicz (Jan. 2008)
! Initial Code thereafter extensively modified by Ronald Errico  
! This version initiated by Ronald Errico (July 2014)
! Initial GOWASP-3 code: Ronald Errico  Oct 1 2016  
!
  use m_kinds, only : rkind1, rkind2
  use m_rad_obs_arrays, only : obs_ndim1,obs_ndim2,obs_ndim3,obs_ndim4
  use m_rad_obs_arrays, only : obs_info,obs_channels,obs_values
  use m_rad_obs_arrays, only : obs_info_num,obs_info_names,obs_info_hdr
  use m_rad_obs_arrays, only : obs_n_channels, obs_n_data
  use m_rad_obs_arrays, only : obs_ndim5,obs_info_extra_recs
  use m_rad_obs_arrays, only : obs_info_extra_names,obs_info_extra
  use m_rad_obs_arrays, only : obs_generic_char, obs_generic_int
!
  implicit none
!
  private
  public :: read_write_obs_rad
  public :: read_write_gmi_1st_msg
  public :: check_rad_type
  public :: bmiss 
!
  real(rkind2), parameter :: bmiss=10.0d10
  character(len=*),parameter :: myname="m_bufr_rw"
!
  contains   
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_rad (lunit,lprint,dtype,nobs,lcopy, &
                                   lread,leof,ier)
!
!  Call appropriate routine to read/write BUFR or generic text file 
!  of radiance observations
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    logical, intent(out)   :: leof
    character(len=*), intent(in) :: dtype   
!   
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_rad'
!
    ier=0
    leof=.false. ! only set true if eof encountered in read genradtxt
!
    if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then
      call read_write_obs_aqua (lunit,ier,lprint,dtype,nobs,lcopy,lread)
!
    elseif (trim(dtype) == 'IASI') then
      call read_write_obs_iasi (lunit,ier,lprint,dtype,nobs,lcopy,lread)
!
    elseif (trim(dtype) == 'GMI') then
     call read_write_obs_gmi (lunit,ier,lprint,dtype,nobs,lcopy,lread)
!
    elseif (trim(dtype) == 'CRIS' .or. trim(dtype) == 'CRISFSR' ) then
     call read_write_obs_cris (lunit,ier,lprint,dtype,nobs,lcopy,lread)
!
    elseif (trim(dtype) == 'SSMIS') then
     call read_write_obs_ssmis (lunit,ier,lprint,dtype,nobs,lcopy,lread)
!
    elseif (dtype(1:4) == 'HIRS' .or. dtype(1:3) == 'MHS' .or.  &
            dtype(1:4) == 'AMSU' .or. dtype(1:3) == 'MSU' .or.  &   
            dtype(1:4) == 'ATMS') then
      call read_write_obs_tovs (lunit,ier,lprint,dtype,nobs,lcopy,lread)
! 
    elseif (trim(dtype) == 'GENRADTXT') then
      call read_write_obs_genradtxt (lunit,ier,lprint,dtype,nobs,lcopy, &
                                     leof,lread)
!
    elseif (trim(dtype) == 'AVCSAM' .or. trim(dtype) == 'AVCSPM' ) then
      call read_write_obs_avhrr (lunit,ier,lprint,dtype,nobs,lcopy,lread)      
!
    elseif (trim(dtype) == 'AMSR2') then
      call read_write_obs_amsr2 (lunit,ier,lprint,dtype,nobs,lcopy,lread)      
!
    else
      print *,' '
      print *,'Invalid option ',trim(dtype),' in ',myname_sub
      ier=66
    endif
!
    end subroutine read_write_obs_rad 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_avhrr (lunit,ier,lprint,dtype,nobs,lcopy, &
                                     lread)
!
!  Read/write bufr file of AVHRR
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=3      ! number of data values per channel
    integer, parameter :: nhead1=13    ! number of header values for hdstr1
    integer, parameter :: nhead2=2     ! number of header values for hdstr2
    integer, parameter :: nhead=nhead1+nhead2
    integer, parameter :: max_ch=5     ! maximum number of channels
    integer :: i,ix,j
    integer :: obs_n_offset
    integer :: khead, kchan
    real(rkind2) :: obsdat(ndval,max_ch)   !  obs profile from bufr file
    real(rkind2) :: hdr(nhead)             !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)     !  output event 
!
    character(len=67) :: hdstr1            ! header info (values actually read)
    character(len=13) :: hdstr2            ! header values not read (set to bmiss)
    character(len=14)  :: str_chan         ! channel data
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_avhrr'
!
    data hdstr1 &
      /'SAID CLATH CLONH YEAR MNTH DAYS HOUR MINU SECO FOVN SAZA SOZA CLAVR'/
    data hdstr2 /'SOLAZI BEARAZ'/
    data str_chan /'INCN ALBD TMBR'/
!
    ier=0 
!
! Set offset since GDAS only uses channels 3,4,5, with n_channels=3
    obs_n_offset=2
    obs_n_channels=3
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1   ! set as flag here; reset later
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*) obs_info_names(1:nhead1)
      read (hdstr2,*) obs_info_names(nhead1+1:nhead)
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
! The values for hdr(nhead1+1:nhead) are looked for by other routines 
      call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
      call ufbrep(lunit,obsdat,ndval,max_ch,kchan,str_chan)
      hdr(nhead1+1:nhead)=bmiss  
!
! copy into main program arrays
      do i=1,obs_n_channels
        ix=i+obs_n_offset
        obs_channels(i,1)=obsdat(1,ix)      ! channel number
        obs_values(i,1)=obsdat(3,ix)        ! brightness temp
        obs_values(i,2)=obsdat(2,ix)        ! albedo
      enddo
      obs_n_data=2  
      obs_info(1:nhead)=hdr(1:nhead)   
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
          obs_info(2) < 100.) then  
        obs_info_extra_recs=1
        obs_info_extra_names(1)='CHNM'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
      endif
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
! 
        hdr(1:nhead)=obs_info(1:nhead)
!
! Set channel-dependent data values 
! NOTE: albedo here are original read values, independent of NR fields
        do i=1,obs_n_channels
          ix=i+obs_n_offset
          bout(1,ix)=obs_channels(i,1)      ! channel numbers 
          bout(3,ix)=obs_values(i,1)        ! brightness temps
          bout(2,ix)=obs_values(i,2)        ! albedo
        enddo 
!
! Fill values not used (channels 1,2) to missing
        do i=1,obs_n_offset
          bout(1,i)=real(i,rkind2)
          bout(3,i)=bmiss  
          bout(2,i)=bmiss
        enddo 
!
! Write out the header information
        call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)

! Write out channel-dependent location-dependent values, incl. TMBR
        call ufbrep(lunit,bout,ndval,max_ch,kchan,str_chan)
        if (max_ch /= kchan) then
          print *,'Error 3 in ',myname_sub
          ier=1
        endif
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for AVHRR'
      else
        print *,' Write observation buffer for AVHRR'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,7f10.3,f12.3)'),' OBS HDR ( 9:13):',hdr(9:13)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,f15.5, 2e15.5)'), &
                 i,obs_channels(i,1),(obs_values(i,j),j=1,2)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_avhrr
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_amsr2 (lunit,ier,lprint,dtype,nobs,lcopy, &
                                     lread)
!
!  Read/write bufr file of AMSR2
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=5      ! number of data values per channel
    integer, parameter :: nhead1=10    ! number of header values in string1
    integer, parameter :: nhead2=8     ! number of header values in string2
    integer, parameter :: nhead3=4     ! number of header values in string3
    integer, parameter :: nhead1p2=nhead1+nhead2
    integer, parameter :: nhead=nhead1p2+nhead3 
    integer, parameter :: max_ch=14    ! maximum number of channels
    integer :: i,j
    integer :: khead, kchan
    real(rkind2) :: salfr
    real(rkind2) :: obsdat(ndval,max_ch)   !  obs profile from bufr file
    real(rkind2) :: hdr(nhead)             !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)     !  output event 
!
! Regarding hdstr3: This is included because although SAZA, BEARAZ and SOZA  
! are not in the BUFR table, subsequent GOWASP_3 routines look for these names
! and reference their values for input to the CRTM. So here, they are included 
! in the complete header string (hdstr1 + hdstr2 + hdstr3) but not in any 
! ufbint call. The values for hdstr3 variables are set following their 
! assignment in GEOS-5 read_amsr2.f90: SAZA is set to that for IANG (but earlier 
! in that routine, before it is used, it is set to 55., but then overwritten by 
! this value; BEARAZ is set to that for AANG; SOZA is set to that for 90.-SOEL.
! SALFR contains a land fraction scalar derived as mean ALFR over channels.
!
    character(len=51) :: hdstr1            ! header variable names string1 
    character(len=41) :: hdstr2            ! header variable names string2
    character(len=22) :: hdstr3            ! header variable names string3
    character(len=26) :: str_chan          ! channel data
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_amsr2'
!
    data hdstr1 /'SIID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH AANG'/
    data hdstr2 /'SAID ORBN SOLAZI SOEL IANG FOVN SLNM ACQF'/
    data hdstr3 /'SAZA BEARAZ SOZA SALFR'/  ! Not in AMSR2 BUFR table 
    data str_chan /'SCCF ALFR VIIRSQ ANPO TMBR'/ 
!
    ier=0 
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1   ! set as flag here; reset later
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*) obs_info_names(1:nhead1)
      read (hdstr2,*) obs_info_names(nhead1+1:nhead1p2)
      read (hdstr3,*) obs_info_names(nhead1p2+1:nhead)
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
      call ufbint(lunit,hdr(nhead1+1:nhead1p2),nhead2,1,khead,hdstr2)
      call ufbrep(lunit,obsdat,ndval,max_ch,obs_n_channels,str_chan)
!
! copy into main program arrays 
      do i=1,obs_n_channels
        obs_channels(i,1)=real(i,rkind2) ! channel report index
        obs_channels(i,2)=obsdat(1,i)    ! frequency (Hz) 
        obs_values(i,1)=obsdat(5,i)      ! brightness temp
        obs_values(i,2)=obsdat(2,i)      ! land fraction
        obs_values(i,3)=obsdat(4,i)      ! antenna polarization
      enddo
      obs_n_data=1                       
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
          obs_info(2) < 100.) then  
        obs_info_extra_recs=2
        obs_info_extra_names(1)='CHNM'
        obs_info_extra_names(2)='SCCF'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
        obs_info_extra(1:obs_n_channels,2)=obs_channels(1:obs_n_channels,2)
      endif
!
! Save info regarding land fraction
      if (obs_n_channels > 0) then 
        salfr=sum(obs_values(1:obs_n_channels,2))/real(obs_n_channels,rkind2)
      else
        salfr=0.0_rkind2
      endif
!
      hdr(19)=hdr(15)             ! SAZA=IANG
      hdr(20)=hdr(10)             ! BEARAZ=AANG
      hdr(21)=90.0_rkind2-hdr(14) ! SOZA=90-SOEL
      hdr(22)=salfr
      obs_info(1:nhead)=hdr(1:nhead)   
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
! 
        hdr(1:nhead)=obs_info(1:nhead)
        salfr=hdr(22)
        obs_channels(1:obs_n_channels,1)=obs_info_extra(1:obs_n_channels,1)
        obs_channels(1:obs_n_channels,2)=obs_info_extra(1:obs_n_channels,2)
!
! Set channel-dependent data values 
        do i=1,obs_n_channels
          bout(1,i)=obs_info_extra(i,2)   ! channel central frequency (Hz) 
          bout(2,i)=salfr                 ! land fraction (as avg over channels) 
          bout(3,i)=real(0.,rkind2)       ! quality flag (set to 0 here)
          bout(4,i)=real(mod(i,2),rkind2) ! antenna polarization flag (assumed)
          bout(5,i)=obs_values(i,1)       ! brightness temp
        enddo 
!
! Write out the header information
        call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
        call ufbint(lunit,hdr(nhead1+1:nhead1p2),nhead2,1,khead,hdstr2)
        call drfini(lunit,obs_n_channels,1,'{AMSRCH}')  

! Write out channel-dependent location-dependent values, incl. TMBR
        call ufbseq(lunit,bout,ndval,obs_n_channels,kchan,'AMSRCH')
        if (obs_n_channels /= kchan) then
          print *,'Error 3 in ',myname_sub
          ier=1
        endif
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for AMSR2'
      else
        print *,' Write observation buffer for AMSR2'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,7f10.3,f12.3)'),' OBS HDR ( 9:16):',hdr(9:16)
      print ('(a,3f12.3)'),' OBS HDR (17:22):',hdr(17:22)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,f15.5,e15.5,f15.5)'), &
                  i,(obs_channels(i,j),j=1,2),obs_values(i,1)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_amsr2
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_aqua (lunit,ier,lprint,dtype,nobs,lcopy, &
                                    lread)
!
!  Read/write bufr file of AIRS obs
!  Each record is assumed to contain possible information for
!  AIRS, AMSUA, and HSB data, although the latter is missing. 
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: n_airs_channels=281 
    integer, parameter :: n_amsu_channels=15
    integer, parameter :: n_max_channels=320
    integer, parameter :: ndval=4  ! number of data values per channel
    integer, parameter :: nhead1=12 ! number of header values in 
    integer, parameter :: nhead2=2 ! number of additional header values
    integer, parameter :: nhead=nhead1+nhead2
    integer :: i, i1, i2
    integer :: instr_type
    real(rkind2) :: hdr1(nhead1,3)               ! instrument header 'SPOT'
    real(rkind2) :: hdr2(nhead2)                 ! instrument header set 2
    real(rkind2) :: obsdat(ndval,n_max_channels) ! channel data
    real(rkind2) :: bout(ndval,n_max_channels)   ! output event 
    integer :: khead1, kchan, khead2
    character(len=8)  :: name_instr_chan, name_instr_spot
    character(len=70) :: hdstr1       ! header string
    data hdstr1 /'SIID CLATH CLONH YEAR MNTH DAYS HOUR MINU SECO SAZA FOVN BEARAZ'/
    character(len=11) :: hdstr2       ! header string 2
    data hdstr2 /'SOZA SOLAZI'/
    character(len=*), parameter :: name_amsu='AMSUAAQUA' 
    character(len=*), parameter :: name_airs='AIRS' 
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_aqua'
!
    ier=0
!
    if (trim(dtype) ==  name_amsu) then
      instr_type=2
      name_instr_chan='AMSUCHAN'
      name_instr_spot='AMSUSPOT'
    else 
      instr_type=1
      name_instr_chan='SCBTSEQN'
      name_instr_spot='AIRSSPOT'
    endif
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1    ! set as flag here; reset later 
      obs_info_hdr=nhead
      obs_info_num=obs_info_hdr
      read (hdstr1,*) obs_info_names(1:nhead1)       ! copy hdstr1 to ...names
      read (hdstr2,*) obs_info_names(nhead1+1:nhead) ! copy hdstr2 to ...names 
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1

! Read header information for 1 report
      call ufbrep(lunit,hdr1,nhead1,3,khead1,hdstr1)
      call ufbint(lunit,hdr2,nhead2,1,khead2,hdstr2)
      if (khead1 /= 3 .or. khead2 /= 1 ) then 
        print *,'Error 1 in in ',myname_sub
        print *,' Error in reading bfr report header in aqua format:'
        print *,' nhead1=',nhead1,' khead1=',khead1,' (should be 3)'
        print *,' nhead2=',nhead2,' khead2=',khead2,' (should be 1)'
	ier=1
      endif
!
! Copy header information
      obs_n_data=1  ! only one type of data: brightness T.       
      obs_info(1:nhead1)=hdr1(1:nhead1,instr_type)    
      obs_info(nhead1+1:nhead)=hdr2(1:nhead2)
!
! Read brightness temperatures
      call ufbseq(lunit,obsdat,ndval,n_max_channels,obs_n_channels, &
                 name_instr_chan)
!
! Copy channel info
      do i=1,obs_n_channels
        obs_values(i,1)=obsdat(4,i)     ! Tb
        obs_values(i,2)=obsdat(3,i)     ! quality flag
        obs_channels(i,1)=obsdat(1,i)   ! Channel number
        obs_channels(i,2)=obsdat(2,i)   ! log-10 of radiance central wavenumber
      enddo
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
          obs_info(2) < 100.) then
        obs_info_extra_recs=2
        obs_info_extra_names(1)='CHNM'
        obs_info_extra_names(2)='LOGRCW'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
        obs_info_extra(1:obs_n_channels,2)=obs_channels(1:obs_n_channels,2)
      endif
!
! End of read of observation buffer
!
    else
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
!
! Write header: order values in hdr1 here so that order is as in 
! 'AIRSPOT' and 'AMSUSPOT'
        hdr1(:,:)=bmiss
        hdr2(:)=bmiss
        hdr1(1,1)=obs_info(1)
        do i2=2,7
          hdr1(i2,1)=obs_info(2+i2)
        enddo
        hdr1( 8,1)=obs_info( 2)
        hdr1( 9,1)=obs_info( 3)
        hdr1(10,1)=obs_info(10)
        hdr1(11,1)=obs_info(12)
        hdr1(12,1)=obs_info(11)
        hdr2(1:nhead2)=obs_info(1+nhead1:nhead)
        call ufbseq(lunit,hdr1(1:nhead1,1),nhead1,1,khead1,name_instr_spot)
        call ufbint(lunit,hdr2(1:nhead2),nhead2,1,khead2,hdstr2)
!
! Extract some information stored elsewhere
        if (lcopy) then
          obs_channels(1:obs_n_channels,1)=obs_info_extra(1:obs_n_channels,1)
          obs_channels(1:obs_n_channels,2)=obs_info_extra(1:obs_n_channels,2)
        endif
!
! Write out channel information (order values as CHNM LOGRCW ACQF TMBR)
        do i=1,obs_n_channels
          if (obs_values (i,1) < bmiss) then
            bout(4,i)=obs_values(i,1)
          else
            bout(4,i)=bmiss
          endif
          bout(3,i)=obs_values(i,2)
          bout(1,i)=obs_channels(i,1)
          bout(2,i)=obs_channels(i,2)
        enddo
!
        if (instr_type == 1) then 
          call drfini(lunit,n_airs_channels,1,'(SCBTSEQN)')  
        else
          call drfini(lunit,obs_n_channels,1,'(SCBTSEQN)')  
        endif
        call ufbseq(lunit,bout,ndval,obs_n_channels,kchan,name_instr_chan)
!
! Check that all required data has been written
        if (kchan /= obs_n_channels) then
          print *,'Error 2 in in ',myname_sub
          print *,' Error in write of airs observation report:' 
          print *,' obs_n_channels=',obs_n_channels,' and kchan=',kchan
          print *,' nhead1=',nhead1,' khead1=',khead1
          print *,' nhead2=',nhead2,' khead2=',khead2
          print *,' name_instr_chan=',name_instr_chan
          print *,' name_instr_spot=',name_instr_spot
	  ier=2
        endif
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer:'
      else
        print *,' Write observation buffer:'
      endif
      print ('(2a,2i8)'),' OBS type, number, channels: ',trim(dtype), &
                          nobs,obs_n_channels
      print ('(a,10f10.3)'),' OBS INFO ( 1:10):',obs_info( 1:10)
      print ('(a,10f10.3)'),' OBS INFO (11:14):',obs_info(11:nhead)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,4f15.5)'),i,obs_channels(i,1:2),obs_values(i,1:2)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_aqua
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_tovs (lunit,ier,lprint,dtype,nobs,lcopy, &
                                    lread)
!
!  Read/write bufr file of HIRS2/3/4, AMSUA/B, MSU, MHS, or ATMS obs
!  (Unlike in the GSI code, there is no data selection here, so the 
!  read/write of ATMS is very similar to that for TOVS and a single code
!  can be used with some special accounting for ATMS.)
!
!  Modifification for ATMS Oct 15 2015 by R. Errico
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=2      ! number of data values per channel
    integer, parameter :: nhead_max=16 ! max num header values for TOVS or ATMS
    integer, parameter :: max_ch=30    ! maximum number of channels
    integer :: i                
    integer :: nhead1,nhead2,nhead
    real(rkind2) :: obsdat(ndval,max_ch)   !  obs profile from bufr file
    real(rkind2) :: hdr(nhead_max)         !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)     !  output event 
    integer :: khead, kchan
    character(len=80) ::  hdstr1,hdstr2    ! header strings
    character(len=10) ::  chstr            ! channel string
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_tovs'
!
    ier=0 
!
    if (trim(dtype) /= 'ATMS') then 
      nhead1=9
      nhead2=7
      hdstr1='SIID CLAT CLON YEAR MNTH DAYS HOUR MINU SECO '
      hdstr2='SAZA FOVN SOZA SAID HOLS BEARAZ SOLAZI' 
      chstr='CHNM TMBR '  ! order matter here
    else  ! for ATMS
      nhead1=10
      nhead2=5
      hdstr1='SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH '
      hdstr2='SAZA SOZA BEARAZ SOLAZI SIID'
      chstr='CHNM TMANT'  ! order matter here
    endif
    nhead=nhead1+nhead2
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1  ! set as flag here; reset later 
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*) obs_info_names(1:nhead1)
      read (hdstr2,*) obs_info_names(nhead1+1:nhead)
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
      call ufbint(lunit,hdr(1+nhead1:nhead),nhead2,1,khead,hdstr2)
      call ufbrep(lunit,obsdat,ndval,max_ch,obs_n_channels,chstr)
!
! copy into main program arrays
      do i=1,obs_n_channels
        obs_channels(i,1)=obsdat(1,i)
        obs_values(i,1)=obsdat(2,i)
      enddo
      obs_n_data=1 ! only one type of data: brightness T.               
      obs_info(1:nhead)=hdr(1:nhead)   
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
                      obs_info(2) < 100.) then
        obs_info_extra_recs=1
        obs_info_extra_names(1)='CHNM'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
      endif
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
!
! Write out channel observations 
        bout=bmiss
        do i=1,obs_n_channels
          bout(2,i)=obs_values(i,1)
          if (lcopy) then 
            bout(1,i)=obs_info_extra(i,1)
          else 
            bout(1,i)=obs_channels(i,1)
          endif
        enddo
!
        if (trim(dtype) == 'ATMS' ) then
          call drfini(lunit,obs_n_channels,1,'(ATMSCH)')  
        endif
        call ufbrep(lunit,bout,ndval,obs_n_channels,kchan,chstr)
        if (obs_n_channels /= kchan) then
          print *,'Error 1 in in ',myname_sub
          ier=1
        endif
!
! Write out the header information
        hdr(1:nhead)=obs_info(1:nhead)
        call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
        call ufbint(lunit,hdr(1+nhead1:nhead),nhead2,1,khead,hdstr2)
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for dtype=',trim(dtype)
      else
        print *,' Write observation buffer for dtype=',trim(dtype)
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,7f10.3,f12.3)'),' OBS HDR ( 9:16):',hdr(9:nhead)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,2f15.5)'),i,obs_channels(i,1),obs_values(i,1)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_tovs
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_iasi (lunit,ier,lprint,dtype,nobs,lcopy, &
                                    lread)
!
!  Read/write bufr file of IASI obs
!  BUFR data are in radiances
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=2     ! number of data values per channel
    integer, parameter :: nhead=19    ! number of header values
    integer, parameter :: max_ch=616  ! max number of channels to be read
    integer :: i, n, n1, n2 
    integer :: iret                   !  number of bufr values returned
    integer :: khead, kchan
    integer :: iscale
    real(rkind2) :: radiance
    real(rkind2) :: obsdat(ndval,max_ch)  !  obs profile from bufr file
    real(rkind2) :: hdr(nhead)            !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)    !  output event
    real(rkind2) :: cldfrac(6)            !  estimated cloud fraction 
    real(rkind2) :: cscale(3,10)          !  coefs for scaling radiances
    real(rkind2) :: scale_fac(10)         !  scaling factor for radiances    
    real(rkind2) :: xdum(10)              !  copied vales for printing
    real(rkind1) :: rad_max(10)           !  max radiance in each scaling group
    real(rkind1) :: alog10_c, alog10_r, alog10_f
    character(len=46) ::  hdstr1  ! header string (cannot be too long)
    character(len=53) ::  hdstr2  ! header string
    character(len=10) ::  chstr   ! channel radiance string
    character(len= 4) ::  cfstr   ! fractional cloud cover string
    character(len=15) ::  csstr   ! radiance scaling coefficients
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_iasi'
    data hdstr1 /'SIID CLATH CLONH YEAR MNTH DAYS HOUR MINU SECO'/
    data hdstr2 /'SAZA FOVN SOZA SOLAZI SAID BEARAZ SLNM QGFQ MJFC SELV'/
    data chstr /'CHNM SCRA '/          ! channel rads; order matters here
    data cfstr /'FCPH'/                ! estimated fractions of cloud cover
    data csstr /'STCH ENCH CHSF '/     ! radiance scaling coefficients
!
    ier=0
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1   ! set as flag here; reset later
      obs_info_hdr=nhead
      read (hdstr1,*) obs_info_names(1:9)
      read (hdstr2,*) obs_info_names(10:nhead)
      do n=1,6
        read (cfstr,*) obs_info_names(nhead+n)
      enddo
      obs_info_num=nhead+6
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      call ufbint(lunit,hdr(1:9),9,1,khead,hdstr1)
      call ufbint(lunit,hdr(10:19),10,1,khead,hdstr2)
      call ufbint(lunit,obsdat,ndval,max_ch,obs_n_channels,chstr)
      call ufbrep(lunit,cscale,3,10,iret,csstr)
      call ufbrep(lunit,cldfrac,1,6,iret,cfstr)
!
! Compute radiance scaling factor from coefficients on file
! The extra exponent of 5 is used to convert W/m2 to mW/m2 as indicated in 
! GSI read_iasi routine. Sometimes very large values of cscale appear and the 
! check is intended to fill those radiances with bad values so they are not used. 
      do n=1,10            
        scale_fac(n)=10.**(5-nint(min(cscale(3,n),100.)))
      enddo
!
! Convert from scaled radiance to unscaled radiance 
      n1=1
      do i=1,obs_n_channels
        obs_channels(i,1)=obsdat(1,i)
        do n=n1,10
          if ( (obsdat(1,i) >= cscale(1,n))   .and. &
               (obsdat(1,i) <= cscale(2,n)) ) then
            radiance=obsdat(2,i)*scale_fac(n)
            n1=n
            exit
          endif
        enddo
        obs_values(i,1)=radiance
      enddo   ! loop over channels
!
      obs_n_data=1  ! only one type of data: brightness T.               
      obs_info(1:nhead)=hdr(1:nhead)    
!
! Copy cldfrac into obs_info array
      n1=nhead+1
      n2=n1+5
      obs_info(n1:n2)=cldfrac(1:6)
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
          obs_info(2) < 100.) then
        obs_info_extra_recs=2
        obs_info_extra_names(1)='CHNM'
        obs_info_extra_names(2)='CSCALE'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
        n2=0
        do n=1,10
          n1=n2+1
          n2=n1+2
          obs_info_extra(n1:n2,2)=cscale(1:3,n)
        enddo
      endif
!
! If array obs_channels has space, copy cscale so if a 
! sample of obs are printed, these values will be printed also
     if (obs_ndim3 > 1) then 
       obs_channels(:,2)=0.
       n2=0
       do n=1,10
         n1=n2+1
         n2=n1+2
         obs_channels(n1:n2,2)=cscale(1:3,n)
       enddo
     endif 
!
! End of read of observation buffer
!
    else
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
!
! Extract some common report information stored elsewhere
! The values cscale(3,:) will be overwritten below
        if (lcopy) then
          obs_channels(1:obs_n_channels,1)=obs_info_extra(1:obs_n_channels,1)
          n1=1
          do n=1,10
            n2=n1+2
            cscale(1:3,n)=obs_info_extra(n1:n2,2)
            n1=n1+3
          enddo
        endif
!
! Extract cldfrac from obs_info array
        n1=nhead+1
        n2=n1+5
        cldfrac(1:6)=obs_info(n1:n2)
!
! Compute radiance scaling factor as in IASI write bufr software
! The extra exponent of 5 in iscale is used to convert W/m2 to mW/m2 
! as indicated in the GSI read_iasi routine
        rad_max(:)=0._rkind1
        alog10_f=1.0_rkind1/log(10._rkind1)   ! conversion factor loge -> log10
        alog10_c=alog10_f*log(2._rkind1**16)  ! log_10 (2**16)
        do n=1,10
          rad_max(n)=0._rkind1
          do i=1,obs_n_channels          
            if (nint(obs_channels(i,1)) >= nint(cscale(1,n)) .and. & 
                nint(obs_channels(i,1)) <= nint(cscale(2,n)) ) then
              rad_max(n)=max(rad_max(n),obs_values(i,1))
            endif
          enddo              
          iscale=5+int(alog10_c-alog10_f*log(rad_max(n)))
          cscale(3,n)=real(iscale)
        enddo
!
        do n=1,10            
          scale_fac(n)=10.**(5-nint(cscale(3,n)))
        end do
!
! Scale radiances and write out 
        bout=bmiss
        n1=1
        do i=1,obs_n_channels
          bout(1,i)=obs_channels(i,1)
          bout(2,i)=obs_values(i,1)
          do n=n1,10
            if ( (bout(1,i) >= cscale(1,n))   .and. &
                 (bout(1,i) <= cscale(2,n)) ) then
              bout(2,i)=bout(2,i)/scale_fac(n)
              n1=n
              exit
            endif
          enddo   ! loop over scaling subgroups
        enddo     ! loop over channels
!
        call drfini(lunit,obs_n_channels,1,'(IASICHN)')  
        call ufbint(lunit,bout,ndval,obs_n_channels,kchan,chstr)
        if (obs_n_channels /= kchan) then
          print *,'Error 1 in in ',myname_sub
          ier=1
        endif
!
! Write out the header information
        hdr(1:nhead)=obs_info(1:nhead)
        call ufbint(lunit,hdr(1:9),9,1,khead,hdstr1)
        call ufbint(lunit,hdr(10:19),10,1,khead,hdstr2)
!
! Write out cloud fraction and radinace scaling parameters
        call ufbrep(lunit,cscale,3,10,iret,csstr)
        call ufbrep(lunit,cldfrac,1,6,iret,cfstr)
!
! If array obs_channels has space, copy cscale so if a 
! sample of obs are printed, these values will be printed also
        if (obs_ndim3 > 1) then 
          obs_channels(:,2)=0.
          n1=1
          do n=1,10
            n2=n1+2
            obs_channels(n1:n2,2)=cscale(1:3,n)
            n1=n1+3
          enddo
        endif 
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer'
      else
        print *,' Write observation buffer'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1:10):',hdr(1:10)
      print ('(a,8f10.3,f12.3)'),' OBS HDR (11:19):',hdr(11:19)
      print ('(a,6f10.3)'),' OBS CLDFRAC:',cldfrac(1:6)
      do i=1,3
        xdum(:)=cscale(i,:)
        print ('(a,i2,10f10.2)'),' OBS CSACALE i=',i,xdum(:)
      enddo
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,1p2e15.5)'),i,obs_channels(i,1),obs_values(i,1)
        enddo
      endif
    endif  ! test on lprint
!
   end subroutine read_write_obs_iasi
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_cris (lunit,ier,lprint,dtype,nobs,lcopy, &
                                   lread)
!
!  Read/write bufr file of CRIS or CRISFSR
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=2      ! number of data values per channel
    integer, parameter :: nhead1=9     ! number of header values for hdstr1
    integer, parameter :: nhead2=11    ! number of header values for hdstr2
    integer, parameter :: nhead=nhead1+nhead2
    integer, parameter :: max_ch=600   ! must be larger than number of channels 
    integer :: i
    integer :: khead, kchan
! rfac changes from W m**2 SR**(-1) cm  to mW m**2 SR**(-1) cm 
    real(rkind2), parameter :: rfac=1000._rkind2 
    real(rkind2) :: obsdat(ndval,max_ch)   !  obs profile from bufr file
    real(rkind2) :: hdr(nhead)             !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)     !  output event 
!
    character(len=46) :: hdstr1    ! header info 
    character(len=59) :: hdstr2    ! header info 
    character(len=9)  :: str_chan  ! channel data
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_cris'
    character(len=6)  :: ncpnum    ! remove this
!
    data hdstr1 /'SAID CLATH CLONH YEAR MNTH DAYS HOUR MINU SECO'/
    data hdstr2 /'FOVN SLNM QMRKH MJFC HMSL FORN SAZA BEARAZ SOZA SOLAZI MTYP'/
    data str_chan /'CHNM SRAD'/
!
    ier=0 
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1   ! set as flag here; reset later
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*) obs_info_names(1:nhead1)
      read (hdstr2,*) obs_info_names(nhead1+1:nhead)
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      call ufbint(lunit,hdr(1:nhead1),      nhead1,1,khead,hdstr1)
      call ufbint(lunit,hdr(nhead1+1:nhead),nhead2,1,khead,hdstr2)
      call ufbrep(lunit,obsdat,ndval,max_ch,kchan,str_chan)
!
! set number of channels to only those used by GDAS 
      if (trim(dtype) == 'CRIS') then
        obs_n_channels=399
      else
        obs_n_channels=431
      endif
!
! copy into main program arrays
      do i=1,obs_n_channels
        obs_channels(i,1)=obsdat(1,i)      ! channel number
        obs_values(i,1)=obsdat(2,i)*rfac   ! radiances (change units)
      enddo
      obs_n_data=1  ! only one type of data: brightness T.         
      obs_info(1:nhead)=hdr(1:nhead)   
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
          obs_info(2) < 100.) then  
        obs_info_extra_recs=1
        obs_info_extra_names(1)='CHNM'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
      endif
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
! 
        hdr(1:nhead)=obs_info(1:nhead)
!
! Special treatment for CRIS FSR   added 2/22
!
        if (obs_n_channels == 431) then      ! assume MTYP='FSR'
           hdr(nhead)=transfer('FSR',1.0_rkind2)    ! assume nhead is index for hdr value MTYP
           if (nobs< 10) print *,'MTYP=',transfer(hdr(nhead),'abc') ! TESTING ONLY
        endif
!
! Set channel-dependent data values 
        do i=1,obs_n_channels
          bout(1,i)=obs_channels(i,1)      ! channel numbers 
          bout(2,i)=obs_values(i,1)/rfac   ! radiances (change units) 
        enddo 
!
! Write out the header information
        call ufbint(lunit,hdr(1:nhead1),      nhead1,1,khead,hdstr1)
        call ufbint(lunit,hdr(nhead1+1:nhead),nhead2,1,khead,hdstr2)

! Write out channel-dependent location-dependent values, incl. TMBR
        if (trim(dtype) == 'CRIS') then 
          call drfini(lunit,obs_n_channels,1,'(CRCHN)')  
        else
          call drfini(lunit,obs_n_channels,1,'(CRCHNM)')  
        endif
        call ufbrep(lunit,bout,ndval,obs_n_channels,kchan,str_chan)
        if (obs_n_channels /= kchan) then
          print *,'Error 3 in in ',myname_sub
          ier=1
        endif
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for ',trim(dtype)
      else
        print *,' Write observation buffer for ',trim(dtype)
      endif
      print ('(a,2i8)'),' OBS number, obs_n_channels :', &
                        nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,7f10.3,f12.3)'),' OBS HDR ( 9:16):',hdr(9:16)
      print ('(a,7f10.3,f12.3)'),' OBS HDR (17:19):',hdr(17:nhead)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,3f15.5)'),i,obs_channels(i,1),obs_values(i,1:2)
        enddo
      endif
     endif  ! lprint
!
    end subroutine read_write_obs_cris
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_gmi (lunit,ier,lprint,dtype,nobs,lcopy, &
                                   lread)
!
!  Read/write bufr file of GMI
!
!  This must be used and be in agreement with routine read_write_gmi_1st_msg.
!  Note that this is a conical scanning instrument with 2 sets of 
!  obs channels, each with their own scan angles
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=4      ! number of data values per channel
    integer, parameter :: nhead1=5     ! number of header values for hdstr1
    integer, parameter :: nhead2=12    ! number of header values for hdstr2
    integer, parameter :: nhead_loc=3  ! number of header values for pix loc
    integer, parameter :: nhead_ang=5  ! number of header values for angles
    integer, parameter :: nhead_qc=3   ! number of derived qc values
!
    integer, parameter :: nhead=nhead1+nhead2+nhead_loc+2*nhead_ang+nhead_qc
    integer, parameter :: max_ch=13     ! maximum number of channels
    integer :: iqcflag(nhead_qc)
    integer :: i,i2,i3,n1,n2                
    integer :: khead, kchan
    integer :: hindex(0:6)  ! indexes for header broken into parts
    real(rkind2), parameter :: default_scan_ang1=52.74_rkind2
    real(rkind2), parameter :: default_scan_ang2=49.2_rkind2
    real(rkind2) :: bmiss9
    real(rkind2) :: default_scan_angs(2)
    real(rkind2) :: obsdat(ndval,max_ch)   !  obs profile from bufr file
    real(rkind2) :: hdr(nhead)             !  location and metadata
    real(rkind2) :: bout(ndval,max_ch)     !  output event 
    real(rkind2) :: val_angles(nhead_ang,2)  ! 2 sets of angle values
!
! Note that the expectation in subroutine crtm_interface_set_gmi is
! that nhead_angles=5 (it uses 5 as an offset for the corresponding 
! indexes between the 2 sets of angles)  
    
    character(len=24) :: hdstr1    ! header info found in NC021200 
    character(len=62) :: hdstr2    ! header info found in NC021204
    character(len=16) :: str_loc   ! pixel location info
    character(len=21) :: str_ang1  ! view geometry for chan 1-9
    character(len=23) :: str_chan  ! channel and loc dependent data
    character(len=26) :: str_ang2  ! view geometry for chan 10-13
    character(len=17) :: str_qc    ! QC flags derived from GMICHQ and GMIRFI
!
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_gmi'
!
    data hdstr1 /'SAID SIID OGCE GSES SACV'/
    data hdstr2 & 
      /'ORBN HMSL SCLAT SCLON GMISQ YEAR MNTH DAYS HOUR MINU SECO SLNM'/
    data str_loc /'CLATH CLONH FOVN'/
    data str_ang1 /'SAZA SAMA SZA SMA SGA'/
    data str_chan /'CHNM TMBR GMICHQ GMIRFI'/
    data str_ang2 /'SAZA2 SAMA2 SZA2 SMA2 SGA2'/
    data str_qc /'CHQC1 CHQC2 CHQC3'/ ! These are not in BUFR Table but pass info in GOWASP
!
    ier=0 
    bmiss9=0.9*bmiss
!
! Set indexs for header info that has been broken into pieces
! hindex(n) is the last index of hdr for the nth broken part
    hindex(0)=0  
    hindex(1)=hindex(0)+nhead1  
    hindex(2)=hindex(1)+nhead2
    hindex(3)=hindex(2)+nhead_loc
    hindex(4)=hindex(3)+nhead_ang
    hindex(5)=hindex(4)+nhead_ang
    hindex(6)=hindex(5)+nhead_qc
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*)   obs_info_names(1:hindex(1))
      read (hdstr2,*)   obs_info_names(hindex(1)+1:hindex(2))
      read (str_loc,*)  obs_info_names(hindex(2)+1:hindex(3))
      read (str_ang1,*) obs_info_names(hindex(3)+1:hindex(4))
      read (str_ang2,*) obs_info_names(hindex(4)+1:hindex(5))
      read (str_qc,*)   obs_info_names(hindex(5)+1:hindex(6))
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      hdr(1:hindex(1))=obs_info(1:hindex(1))  ! already read in ...gmi_1st_msg 
      call ufbint(lunit,hdr(hindex(1)+1:hindex(2)),nhead2,1,khead,hdstr2)
      call ufbint(lunit,hdr(hindex(2)+1:hindex(3)),nhead_loc,1,khead,str_loc)
      call ufbrep(lunit,val_angles,nhead_ang,2,khead,str_ang1) 
      call ufbrep(lunit,obsdat,ndval,max_ch,obs_n_channels,str_chan)
!
! copy into main program arrays
      do i=1,obs_n_channels
        obs_channels(i,1)=obsdat(1,i)   ! channel number
        obs_channels(i,2)=obsdat(3,i)   ! quality mark (BUFR values for GMICHQ) 
        obs_values(i,1)=obsdat(2,i)     ! brightness T
        obs_values(i,2)=obsdat(4,i)     ! 2nd quality mark (BUFR values for GMIRFI)
      enddo
!
! Simply copy scan angle information for the two channel sets (1-9, 10-13)
! Any required modifications will be performed in crtm_interface
      hdr(hindex(3)+1:hindex(4))=val_angles(:,1)
      hdr(hindex(4)+1:hindex(5))=val_angles(:,2)
!        
! Set the channel-dependent quality control marker as a single scalar  
      iqcflag(:)=0
      do i=1,obs_n_channels
        i3=2**(i-1)
        if (obs_channels(i,2) > 1.5 .or. obs_channels(i,2) < -0.5 ) then 
          iqcflag(1)=iqcflag(1)+i3
        endif
        if (obs_values(i,2) > 0.1) then 
          iqcflag(2)=iqcflag(2)+i3
        endif
        if (obs_values(i,1) > bmiss9) then 
          iqcflag(3)=iqcflag(3)+i3
        endif
      enddo
!
      do i=1,nhead_qc
        hdr(hindex(5)+i)=real(iqcflag(i))
      enddo
!
      obs_n_data=1  ! only one type of data: brightness T.         
      obs_info(hindex(1)+1:nhead)=hdr(hindex(1)+1:nhead)   
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
! 
        hdr(1:nhead)=obs_info(1:nhead)
        bout=bmiss  ! initialize default as value missing
!
! Set channel-dependent data and QC values based on scalar quality flag 
! Non-acceptable QC values are set to -5 indicating multiple missing values 
        do i=1,nhead_qc
          iqcflag(i)=nint(hdr(hindex(5)+i))
        enddo
! 
        do i=obs_n_channels,1,-1
          i3=2**(i-1)
          if (iqcflag(1)/i3 > 0) then 
            obs_channels(i,2)=real(-2)  
            iqcflag(1)=iqcflag(1)-i3
          else
            obs_channels(i,2)=0.        ! QC=OK
          endif
!
          if (iqcflag(2)/i3 > 0) then 
            obs_values(i,2)=real(-2)  
            iqcflag(2)=iqcflag(2)-i3
          else
            obs_values(i,2)=0.          ! QC=OK
          endif
!
          if (iqcflag(3)/i3 > 0) then 
            obs_values(i,1)=bmiss
            iqcflag(3)=iqcflag(3)-i3
          endif
        enddo
!
        do i=1,obs_n_channels
          bout(1,i)=obs_channels(i,1)   
          bout(3,i)=obs_channels(i,2)   ! quality flag
          bout(4,i)=obs_values(i,2)     ! quality flag
          bout(2,i)=obs_values(i,1)     ! TB; otherwise, default is bmiss
        enddo 
!
! Write out the header information
        call ufbint(lunit,hdr(hindex(1)+1:hindex(2)),nhead2,1,khead,hdstr2)
        call ufbint(lunit,hdr(hindex(2)+1:hindex(3)),nhead_loc,1,khead,str_loc)
!
! Write out viewing and sun angle variables for 2 sets of channels 
        val_angles(:,1)=hdr(hindex(3)+1:hindex(4))
        val_angles(:,2)=hdr(hindex(4)+1:hindex(5))
        call ufbrep(lunit,val_angles,nhead_ang,2,khead,str_ang1)    
!
! Write out channel-dependent location-dependent values, incl. TMBR
        call ufbrep(lunit,bout,ndval,obs_n_channels,kchan,str_chan)
        if (obs_n_channels /= kchan) then
          print *,'Error 3 in in ',myname_sub
          ier=1
        endif
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for GMI'
      else
        print *,' Write observation buffer for GMI'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,8f10.3)'),' OBS HDR ( 9:16):',hdr(9:16)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,4f15.5)'),i,obs_channels(i,1:2),obs_values(i,1:2)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_gmi
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    subroutine read_write_gmi_1st_msg (lunit,lread,idate,ier)
!
! Special instructions required for handling GMI bufr files:
! Unlike other bufr files for radiances, the GMI files contain 2 distinct
! message types. Subtype NC021200 includes only header and channel information 
! that is identical for all obs locations. That subtype is read or written here.
! Subtype NC021204 includes all location dependent information. 
!
    implicit none
!
    logical, intent(in)  :: lread
    integer, intent(in)  :: lunit
    integer, intent(in)  :: idate
    integer, intent(out) :: ier
!
    integer, parameter :: ncval=4      ! number of channel descriptor info
    integer, parameter :: nhead1=5     ! number of header values for hdstr1
    integer, parameter :: max_ch=13     ! maximum number of channels
    integer :: i
    integer :: idate1
    integer :: ier1
    integer :: ireadmg, ireadsb
    real(rkind2) :: hdr(nhead1)             
    real(rkind2) :: chnm(ncval,max_ch)             
    character(len=8) :: msg,msg1
    character(len=24) :: hdstr1    ! header info found in NC021200 
    character(len=19) :: chnmstr   ! channel info that is loc independent
!
    data hdstr1  /'SAID SIID OGCE GSES SACV'/
    data msg1    /'NC021200'/  ! label that shoud be for 1st GMI msg
    data chnmstr /'CHNM SCCF SCBW ANPO'/
!        
    if (lread) then 
      ier1=ireadmg(lunit,msg,idate1) 
      if (msg /= msg1 .or. ier1 /= 0) then
        ier=ier1
        print *,' '
        print *,'ERROR: Unexpected subtype in 1st message in GMI file: ', &
                'ier=',ier1,'  msg=',msg,' not ',msg1
      else  
        do while (ireadsb(lunit) == 0)
          call ufbint(lunit,hdr,nhead1,1,ier1,hdstr1) 
          call ufbrep(lunit,chnm,ncval,max_ch,obs_n_channels,chnmstr)
!
          read (chnmstr,*) obs_info_extra_names(1:ncval)  ! save variable names
          obs_info(1:nhead1)=hdr(1:nhead1)   
          obs_info_extra_recs=ncval
          do i=1,obs_n_channels
            obs_info_extra(i,1:ncval)=chnm(1:ncval,i)
          enddo
          ier=0 
        enddo
      endif
!
      print *,' '
      print *,'First msg read from GMI bufr file'
      print ('(a,8f13.3)'),'hdr=',hdr(1:nhead1)
!
    else ! written 1st message containing report for loc-independent header info
!
      hdr(1:nhead1)=obs_info(1:nhead1)
      do i=1,obs_n_channels
        chnm(1:ncval,i)=obs_info_extra(i,1:ncval)
      enddo
!
      call openmg(lunit,msg1,idate)   
      call ufbint(lunit,hdr,nhead1,1,ier1,hdstr1) 
      call ufbrep(lunit,chnm,ncval,obs_n_channels,ier1,chnmstr)
      call writsb(lunit)              
      call closmg(lunit)
!
      print *,'  '
      print *,'First msg written to GMI bufr file'
      print ('(a,8f13.3)'),'hdr=',hdr(1:nhead1)
!
    endif 
!
    end subroutine read_write_gmi_1st_msg 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_ssmis (lunit,ier,lprint,dtype,nobs,lcopy, &
                                    lread)
!
!  Read/write bufr file of SSMIS.
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    character(len=*), intent(in) :: dtype   
!
! local variables
!
    integer, parameter :: ndval=2      ! number of data values per channel
    integer, parameter :: nhead1=9
    integer, parameter :: nhead2=8
    integer, parameter :: nhead3=4
    integer, parameter :: nhead12=nhead1+nhead2
    integer, parameter :: nhead=nhead12+nhead3
    integer, parameter :: nhead1p1=nhead1+1
    integer, parameter :: max_ch=24    ! number of channels
    integer :: i                
    integer :: khead, kchan
!
    real(rkind2) :: data_ch(2,max_ch)  ! chan info in     
    real(rkind2) :: data_ymd(3,5)      ! year month day in
    real(rkind2) :: data_hm(2,2)       ! hour min in
    real(rkind2) :: data_xy(2,29)      ! lat lon in
    real(rkind2) :: hdr(nhead)         
!
    character(len=80) :: hdstr1,hdstr2,hdstr3  ! header strings
    character(len=80) :: chstr                 ! channel string
    character(len=*), parameter :: myname_sub=myname//'::read_write_obs_ssmis'
!
    ier=0 
!
    hdstr1='SAID FOVN SFLG RFLAG SLNM ORBN TSIG SSID DIMS'
    hdstr2='CLAT CLON YEAR MNTH DAYS HOUR MINU SECO'  ! must be in this order
    hdstr3='SAZA BEARAZ SOZA SOLAZI'                  ! must be in this order
    chstr='CHNM TMBR'                                 ! must be in this order
!
! nobs<0 signifies that only header and other dtype-general information 
! are to be retrieved from this subroutine.
    if (nobs < 0) then
      obs_info_extra_recs=-1  ! set as flag here; reset later 
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdstr1,*) obs_info_names(1:nhead1)
      read (hdstr2,*) obs_info_names(nhead1+1:nhead12)
      read (hdstr3,*) obs_info_names(nhead12+1:nhead)
      obs_info_extra(:,:)=0._rkind2
!
! Check if this is a read or write
    elseif (lread) then
      nobs=nobs+1
!
! read information from the input file needed for interpolation
      call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
      call ufbint(lunit,hdr(nhead12),1,1,khead,'SECO')
      call ufbrep(lunit,data_ymd,3,5,khead,'YEAR MNTH DAYS')
      call ufbrep(lunit,data_hm,2,2,khead,'HOUR MINU')
      call ufbrep(lunit,data_xy,2,29,khead,'CLAT CLON')
      call ufbrep(lunit,data_ch,2,max_ch,obs_n_channels,'CHNM TMBR')
!
! copy into main program arrays
      obs_n_data=1  ! only one type of data: brightness T.               
      do i=1,obs_n_channels
        obs_channels(i,1)=data_ch(1,i)
        obs_values(i,1)=data_ch(2,i)
      enddo
!
! Copy header info 
      obs_info(1:nhead1)=hdr(1:nhead1)   
      obs_info(nhead1+1:nhead1+2)=data_xy(1:2,1)  ! CLAT CLON
      obs_info(nhead1+3:nhead1+5)=data_ymd(1:3,1) ! YEAR MNTH DAYS
      obs_info(nhead1+6:nhead1+7)=data_hm(1:2,1)  ! HOUR MINU
      obs_info(nhead1+8)=hdr(nhead1+8)            ! SECO      
!
! Set look angles and sun angles
      obs_info(nhead12+1)=53.   ! SAZA
      obs_info(nhead12+2)=0.    ! BEARAZ
      call angles_ssmis (obs_info(nhead1+3:nhead1+8),   &
               obs_info(nhead1+1),obs_info(nhead1+2),   &
               obs_info(nhead12+3),obs_info(nhead12+4)) ! SOZA and SOLAZI
!
! Copy some information for use by other routines
! (This is done only once, using first report with non-missing information)
      if (lcopy .and. obs_info_extra_recs < 0 .and. &
                      obs_info(2) < 100.) then
        obs_info_extra_recs=1
        obs_info_extra_names(1)='CHNM'
        obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
      endif
!
      if (lprint) then  
        hdr(1:nhead)=obs_info(1:nhead)  ! copied here for printing below 
      endif
!
! End of read of observation buffer
!
    else 
!
! This is a call to write observation buffer
!
      if (obs_n_channels > 0) then
!
! Write header
        hdr(1:nhead)=obs_info(1:nhead)
        call ufbint(lunit,hdr(1:nhead1),nhead1,1,khead,hdstr1)
        call ufbint(lunit,hdr(nhead1+nhead2),1,1,khead,'SECO')
!
        data_xy(1,:)=obs_info(nhead1+1)     ! CLAT
        data_xy(2,:)=obs_info(nhead1+2)     ! CLON
        call ufbrep(lunit,data_xy,2,29,khead,'CLAT CLON')
!
        data_ymd(1,:)=obs_info(nhead1+3)    ! YEAR
        data_ymd(2,:)=obs_info(nhead1+4)    ! MNTH
        data_ymd(3,:)=obs_info(nhead1+5)    ! DAYS
        call ufbrep(lunit,data_ymd,3,5,khead,'YEAR MNTH DAYS')
!
        data_hm(1,:)=obs_info(nhead1+6)     ! HOUR
        data_hm(2,:)=obs_info(nhead1+7)     ! MINU
        call ufbrep(lunit,data_hm,2,2,khead,'HOUR MINU')
!
        do i=1,obs_n_channels
          if (lcopy) then 
            data_ch(1,i)=obs_info_extra(i,1)
          else 
            data_ch(1,i)=obs_channels(i,1)
          endif
          data_ch(2,i)=obs_values(i,1)
        enddo
        call ufbrep(lunit,data_ch,2,max_ch,obs_n_channels,'CHNM TMBR')
!
!  finished updates; so write output buffer to output file for this report
!
        call writsb(lunit)
!
      endif  ! test on obs_n_channels
!
    endif  ! test on whether to read or write observation buffer
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for SSMIS'
      else
        print *,' Write observation buffer for SSMIS'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,10f10.3)'),' OBS HDRSTR1:',hdr(1:nhead1)
      print ('(a,10f10.3)'),' OBS HDRSTR2:',hdr(nhead1+1:nhead12)
      print ('(a,10f10.3)'),' OBS HDRSTR3:',hdr(nhead12+1:nhead)
      if (obs_n_channels > 0) then
        do i=1,obs_n_channels
          print ('(i4,2f15.5)'),i,obs_channels(i,1),obs_values(i,1)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_ssmis
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
    Subroutine read_write_obs_genradtxt (lunit,ier,lprint,dtype,nobs, &
                                         lcopy,leof,lread)
!
!  Read/write generic text file for radiance observations
!
    implicit none
!
    integer, intent(inout) :: nobs
    integer, intent(out)   :: ier
    integer, intent(in)    :: lunit
    logical, intent(in)    :: lprint
    logical, intent(in)    :: lcopy
    logical, intent(in)    :: lread
    logical, intent(out)    :: leof
    character(len=*), intent(in) :: dtype   
!
! For generic radiances, the sat name and intrument will be provided 
! by names provided in the file, not using WMO id numbers
    integer, parameter :: siid=999.  ! flag for generic
    integer, parameter :: nhead_max=20
    integer :: nhead1,nhead2,nhead
    integer :: ifov
    integer :: n
    integer :: itime(6)
    integer :: kidsat,idate,nchanl  ! integer variables in file header
    real(rkind2) :: hdr(nhead_max)         !  location and metadata
    character(len=80)  ::  hdrstr1,hdrstr2   
    character(len=100) :: instfmt,reclocfmtr,reclocfmtw,recgeofmt,obfmt
    character(len=20)  :: sis
    character(len=10)  :: obstype,jsatid
!
! set text formatting strings
    instfmt   = '(A20,1x,A10,3x,A10,3x,I5,3x,i10,3x,i5)'
    reclocfmtr = '(F9.4,1x,F9.4,10x,I4,"-",I2,"-",I2,1x,I2,":",I2,":",I2)'
    reclocfmtw = '(F9.4,1x,F9.4,10x,I4.4,"-",I2.2,"-",I2.2,1x,I2.2,":",I2.2,":",I2.2)'
    recgeofmt = '(I4,3x,5(F8.3,1x))'
    nhead1=9
    nhead2=7
    nhead=nhead1+nhead2
    leof=.false.                   ! default value
    ier=0
!
    if (nobs < 0) then             ! initialization block
!
! hdrstr1 variables are in obs report header record 1 (except SAID)
! hdrstr1 variables are in obs report header record 2 (except SIID)
      hdrstr1   = 'SAID CLONH CLATH YEAR MNTH DAYS HOUR MINU SECO'
      hdrstr2   = 'FOVN PANGLR SAZA BEARAZ SOZA SOLAZI SIID'
      obs_info_hdr=nhead
      obs_info_num=nhead
      read (hdrstr1,*) obs_info_names(1:nhead1)
      read (hdrstr2,*) obs_info_names(nhead1+1:nhead)
      if (lread) then                   ! read file header
        read (lunit,fmt=instfmt) sis,obstype,jsatid,kidsat,idate,nchanl
        obs_generic_char(1)=sis
        obs_generic_char(2)='NOTUSED'   ! no longer used
        obs_generic_char(3)=obstype
        obs_generic_char(4)=jsatid
        obs_generic_char(5)='GENRAD'
        obs_generic_int(1)=kidsat
        obs_generic_int(2)=idate
        obs_generic_int(3)=nchanl
     else                               ! write file header
        sis=obs_generic_char(1)
        obstype=obs_generic_char(3)
        jsatid=obs_generic_char(4)
        kidsat=obs_generic_int(1)
        idate=obs_generic_int(2)
        nchanl=obs_generic_int(3)
        write (lunit,fmt=instfmt) sis,obstype,jsatid,kidsat,idate,nchanl
      endif
    endif   
!
    obs_n_channels=obs_generic_int(3) 
    obs_n_data=1  ! only one type of data: brightness T.          
    kidsat=obs_generic_int(1)     
    do n=1,obs_n_channels
      obs_channels(n,1)=real(n)  ! simple setting of channel numbers
    enddo
    write (obfmt,fmt='(a,i,a)') '(8x,',obs_n_channels,'(f10.4))'
!
    if (lcopy .and. lread .and. nobs < 0) then  ! initialization block
      obs_info_extra_recs=1
      obs_info_extra_names(1)='CHNM'
      obs_info_extra(1:obs_n_channels,1)=obs_channels(1:obs_n_channels,1)
    endif
! 
    if (nobs > -1 .and. lread) then    ! read data records
      leof=.true.                      ! reset below if eof not encountered
      read(lunit,fmt=reclocfmtr,end=500) hdr(2:3),itime(1:6)
      read(lunit,fmt=recgeofmt,end=500) ifov,hdr(nhead1+2:nhead-1)
      read(lunit,fmt=obfmt) obs_values(1:obs_n_channels,1)
      nobs=nobs+1                      ! increment counter
      hdr(1)=real(kidsat)
      hdr(nhead1+1)=real(ifov)
      hdr(nhead)=siid
      do n=4,9
        hdr(n)=real(itime(n-3))
      enddo
      obs_info(1:nhead)=hdr(1:nhead)   ! copy obs header info
      leof=.false.                     ! if this reached, eof not encountered
500   continue                         ! if this reached, eof found
    endif
!
    if (nobs > -1 .and. (.not. lread)) then  ! write data records
      hdr(1:nhead)=obs_info(1:nhead)   ! copy obs header info
      ifov=nint(hdr(nhead1+1))
      do n=1,6
        itime(n)=nint(hdr(n+3))
      enddo
      write (lunit,fmt=reclocfmtw) hdr(2:3),itime(1:6)
      write (lunit,fmt=recgeofmt) ifov,hdr(nhead1+2:nhead-1)
      write (lunit,fmt=obfmt) obs_values(1:obs_n_channels,1)
    endif
!
! Optionally print some info if testing mode
    if (lprint .and. nobs < 3  .and. nobs > 0 .and. (.not. leof) ) then
      print *,' '
      if (lread) then
        print *,' Read observation buffer for GENRADTXT'
      else
        print *,' Write observation buffer for GENRADTXT'
      endif
      print ('(a,2i8)'),' OBS number, channels:',nobs,obs_n_channels
      print ('(a,8f10.3)'),' OBS HDR ( 1: 8):',hdr(1:8)
      print ('(a,7f10.3,f12.3)'),' OBS HDR ( 9:16):',hdr(9:16)

      if (obs_n_channels > 0) then
        do n=1,obs_n_channels
          print ('(i4,2f15.5)'),n,obs_channels(n,1),obs_values(n,1)
        enddo
      endif
    endif  ! test on lprint
!
    end subroutine read_write_obs_genradtxt
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine check_rad_type (subset,dtype,lcheck)
!
! Code for making data selection of subtypes given the data type requested:
! subset is the data subset type name read from the BUFR files
! dtype is the data type specified in the main program argument list
! The subset name 'ANY' is included so that this routine can be used to 
! simply test if the data type requested is any one of these here. 
!
   implicit none   
   character(len=*), intent(in) :: subset
   character(len=*), intent(in) :: dtype
   logical, intent(out) :: lcheck
!
! For GMI, NC021200 contains loc-independent channel and header info
!
   character(len=8), parameter ::   gmi_types(3) = &
                                        (/'NC021204','NC021200','ANY'/)
   character(len=8), parameter ::   atms_types(2) = (/'NC021203','ANY'/)
   character(len=8), parameter ::  hirs2_types(2) = (/'NC021021','ANY'/)
   character(len=8), parameter ::  hirs3_types(2) = (/'NC021025','ANY'/)
   character(len=8), parameter ::  hirs4_types(2) = (/'NC021028','ANY'/)
   character(len=8), parameter ::  amsua_types(2) = (/'NC021023','ANY'/)
   character(len=8), parameter ::  amsub_types(2) = (/'NC021024','ANY'/)
   character(len=8), parameter ::    mhs_types(2) = (/'NC021027','ANY'/)
   character(len=8), parameter ::    msu_types(2) = (/'NC021022','ANY'/)
   character(len=8), parameter ::   airs_types(2) = (/'NC021250','ANY'/)
   character(len=8), parameter ::   iasi_types(2) = (/'NC021241','ANY'/)
   character(len=8), parameter ::   cris_types(2) = (/'NC021202','ANY'/)
   character(len=8), parameter :: crisf4_types(2) = (/'NC021206','ANY'/)
   character(len=8), parameter ::  ssmis_types(2) = (/'NC021201','ANY'/)
   character(len=8), parameter ::  avhrr_types(3) = (/'NC021051','NC021053','ANY'/)
   character(len=8), parameter ::  amsr2_types(2) = (/'NC021248','ANY'/)
   character(len=8), parameter :: amsua_aqua_types(2) = (/'NC021249','ANY'/)  
   character(len=8), parameter :: genradtxt_types(2) = (/'GENRAD','ANY'/)  
!
   lcheck=.false.
   select case (trim(dtype))
     case ('GMI')
       lcheck = any( subset == gmi_types )
     case ('ATMS')
       lcheck = any( subset == atms_types )
     case ('HIRS2')
       lcheck = any( subset == hirs2_types )
     case ('HIRS3')
       lcheck = any( subset == hirs3_types )
     case ('HIRS4')
       lcheck = any( subset == hirs4_types )
     case ('AIRS')
       lcheck = any( subset == airs_types )
     case ('IASI')
       lcheck = any( subset == iasi_types )
     case ('AMSUA')
       lcheck = any( subset == amsua_types )
     case ('AMSUB')
       lcheck = any( subset == amsub_types )
     case ('MHS')
       lcheck = any( subset == mhs_types )
     case ('MSU')
       lcheck = any( subset == msu_types )
     case ('CRIS')
       lcheck = any( subset == cris_types )
     case ('CRISFSR')
       lcheck = any( subset == crisf4_types )
     case ('AVCSAM')
       lcheck = any( subset == avhrr_types )
     case ('AVCSPM')
       lcheck = any( subset == avhrr_types )
     case ('AMSR2')
       lcheck = any( subset == amsr2_types )
     case ('SSMIS')
       lcheck = any( subset == ssmis_types )
     case ('AMSUAAQUA')
       lcheck = any( subset == amsua_aqua_types )
     case ('GENRADTXT')
       lcheck = any( subset == genradtxt_types )
   end select
!
   end subroutine check_rad_type 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_bufr_rad
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   logical function generic_ireadsb (iunit,leof,ctype)
!
!  A generalization of the BUFR library function ireadsb that is suitable 
!  for application to either BUFR or generic text files.
!
   implicit none
   logical, intent(in) :: leof
   integer, intent(in) :: iunit
   character(len=*), intent(in) :: ctype
!
   integer ::   ireadsb  ! a data read function in BUFR lib 
   if (trim(ctype) == 'BUFR') then
     generic_ireadsb=(ireadsb(iunit) == 0)
   else
     generic_ireadsb=(.not. leof)
   endif
   return
   end function generic_ireadsb
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   logical function generic_ireadmg (iunit,leof,ctype,subset,idate)
!
!  A generalization of the BUFR library function ireadmg that is suitable 
!  for application to either BUFR or generic text files.
!
   use m_rad_obs_arrays, only: obs_generic_char, obs_generic_int
   implicit none
   logical, intent(in) :: leof
   integer, intent(in) :: iunit
   integer, intent(out) :: idate
   character(len=*), intent(in) :: ctype
   character(len=8), intent(out) :: subset
!
   integer ::   ireadmg  ! a message read function in BUFR lib 
   if (trim(ctype) == 'BUFR') then
     generic_ireadmg=(ireadmg(iunit,subset,idate) == 0)
   else
     subset=obs_generic_char(5)(1:8)
     idate=obs_generic_int(2)
     generic_ireadmg=(.not. leof)
   endif
!
   return
   end function generic_ireadmg
  





