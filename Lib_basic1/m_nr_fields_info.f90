    module m_nr_fields_info
!
! Module to read .rc file and setup field requirements for 
! NR fields to be read elsewhere
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
    use m_kinds, only : rkind1 
    implicit none
!
    private
!
    public :: nr_fields_setup
    public :: nr_fields_print_info
!
    private :: nr_fields_read_rc
!
    public :: field_imax, field_jmax, field_kmax 
    public :: field_kdim, field_lev_first, field_top 
    public :: field_lon_first, field_tavg_offset, field_avg
    public :: field_max_akbk, field_akbk, field_akbk_dlev
    public :: field_time_slots, field_time_delta, field_time_first
    public :: field_max_num, field_num_2d, field_num_3d
    public :: field_common_path, field_names, field_files, field_types
    public :: field_obs_types_num, field_obs_types 
    public :: field_ranseed
    public :: field_pole_lat, field_pole_lons
    public :: field_ipw_top, field_ipw_bot
    public :: field_stride_i, field_stride_j, field_satwind_types
    public :: field_gpsro_ksmooth, field_gpsro_algorithm
    public :: field_conv_ksmooth
!
! Information describing fields needed to construct profiles
! Last index on _nums, _names, _files indicates (2) 2d or (3) 3d
! First index on _names is (1) name used by osse sim software, (2) name on file 
    integer, parameter :: field_max_num=60      ! max number of 2d or 3d fields
    integer, parameter :: field_max_akbk=201    ! max number of eta levels
!
    integer(8) :: field_ranseed  ! random seed to use if required 
    integer :: field_imax        ! num lons on grid 
    integer :: field_jmax        ! num lats on grid
    integer :: field_kmax        ! num of levels for full 3D fields 
    integer :: field_kdim        ! number of retained levels of 3D fields
    integer :: field_lev_first   ! first level index of retained 3D levels 
    integer :: field_format      ! index for format of rc file
    integer :: field_time_slots  ! number of time slots in assim. window
    integer :: field_num_2d      ! number of 2d fields to be processed
    integer :: field_num_3d      ! number of 3d fields to be processed
    integer :: field_obs_types_num    ! number of obs subtypes to consider
    integer :: field_conv_ksmooth(3) ! kmin,kmax,ksmooth (avg 2*ksmooth+1 levs) 
    integer :: field_gpsro_ksmooth(3) ! kmin,kmax,ksmooth (avg 2*ksmooth+1 levs) 
    integer :: field_gpsro_algorithm  ! 0=2-d ROPP, 1=1-d Ropp, 2=GSI, 3=mod GSI
    integer :: field_types(field_max_num,2:3) ! =nk, with - sign if no interp
    integer :: field_obs_types(200)           ! to hold list of requested types
    real(rkind1) :: field_top    ! top-most p (units Pa) for sounding pressure
    real(rkind1) :: field_tavg_offset  ! offset (hr) for time-averaged files
    real(rkind1) :: field_time_delta   ! hours between files  
    real(rkind1) :: field_time_first   ! hours for 1st file relat to assim time 
    real(rkind1) :: field_lon_first    ! longitude for index i=1
    real(rkind1) :: field_avg(field_max_num,2:3)  ! size (km) for area avg of f 
    real(rkind1) :: field_akbk(field_max_akbk,2)       ! 1 for ak, 2 for bk
    real(rkind1) :: field_akbk_dlev(field_max_akbk,2)  ! akbk for data levels
    character(len=240) :: field_common_path            ! path for NR dat files
    character(len=12)  :: field_names(field_max_num,2,2:3)
    character(len=240) :: field_files(field_max_num,2:3)
!
! The following are only used when making conventional obs
    real(rkind1) :: field_pole_lons(12) ! selected longitudes for pole advect
    real(rkind1) :: field_pole_lat     ! lat adjacent to N. pole
!
! The following are only used when making satwind obs
    real(rkind1), parameter :: field_ipw_top=10000. ! top p (Pa) for ipw calc
    real(rkind1), parameter :: field_ipw_bot=60000. ! bottom p (Pa) for ipw calc
    integer :: field_stride_i    ! consider every nth lon for satwind loc
    integer :: field_stride_j    ! consider every nth lat for satwind loc
    character(len=4)   :: field_satwind_types  ! selects satwind types to create
!
    character(len=*),parameter :: myname='m_nr_fields_info'
!
    contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine nr_fields_setup (dtype,rcfile,lprint,ier)
!
!  Gather and set information required to read and use NR field information 
!
    implicit none
!
! arguments
    logical, intent(in)  :: lprint  
    integer, intent(out) :: ier
    character(len=*), intent(in) :: dtype
    character(len=*), intent(in) :: rcfile
!
! local variables
    integer :: k,k1
    integer :: i12
    real(rkind1), parameter :: ps_rep=1.0e5 ! a "representative" ps
    real(rkind1) :: p,p1,p2
    real(rkind1) :: pole_lon
    character(len=*), parameter :: mysub=myname//'::nr_fields_setup'
!
    call nr_fields_read_rc (dtype,rcfile,lprint,ier)
!
    if (ier == 0) then 
!
! Find lev idexes for top of the required portion of all 3d fields  
      field_lev_first=1
      do k=1,field_kmax-2  ! willl be at least 3 data levels
        p=field_akbk(k,1)+field_akbk(k,2)*ps_rep
        if (p < field_top) then
          field_lev_first=k
        endif
      enddo
      field_kdim=field_kmax-field_lev_first+1
!
! Adjust indexes of levels to be smoothed for raobs since not all levs used for sim conv
! (This concerns correcting for noisy stratosphere T's for the GEOS-5 NR)
      if (field_conv_ksmooth(3) > 0) then
        field_conv_ksmooth(1)=max(field_conv_ksmooth(1)-field_lev_first+1,field_conv_ksmooth(3)+1)
        field_conv_ksmooth(2)=field_conv_ksmooth(2)-field_lev_first+1
        if (field_conv_ksmooth(2) < field_conv_ksmooth(1)) then
          field_conv_ksmooth(:)=0
        endif
      endif
!
! Shift ak, bk levs to start with field_lev_first
      k=0
      do k1=field_lev_first,field_kmax+1 
        k=k+1
        field_akbk(k,1:2)=field_akbk(k1,1:2)
      enddo
!
! Compute ak and bk at data levels
      do k=1,field_kdim
        field_akbk_dlev(k,:)=0.5*(field_akbk(k,:)+field_akbk(k+1,:))
      enddo
      k=field_kdim+1
      field_akbk_dlev(k,:)=field_akbk(k,:)
!
! The following is required for conventional obs only.
! Specifically, these are required for computing advection near the poles
! Specify value of grid latitude adjacent to north pole
! Also specify 12 selected longitudes used if advection near pole
      field_pole_lat=90.-180./(field_jmax-1)
      i12=max(1,field_imax/12)    
      do k=1,12
        pole_lon=field_lon_first+(k-1)*i12*360./real(field_imax)
        field_pole_lons(k)=mod(pole_lon,360.)
      enddo
!
! Print summary of information read from .rc file
      if (lprint) then
        call nr_fields_print_info (dtype,mysub)
      endif
!
    endif  ! test on ier=0
!
    end subroutine nr_fields_setup 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine nr_fields_read_rc (dtype,rcfile,lprint,ier)
!
!  Read .rc file describing properties and selection of required NR fields 
!
    implicit none
!
! arguments
    logical, intent(in) :: lprint
    integer, intent(out) :: ier
    character(len=*), intent(in) :: dtype
    character(len=*), intent(in) :: rcfile
!
! local variables
    logical :: found                    ! true if dtype found in rc file
    logical :: not_rad                  ! true if dtype not a radiance dtype
    logical, parameter :: lstop=.false. ! must be false when mpi is used
!
    integer, parameter :: read_max=100  ! >= # of rc header records
    integer :: field_format_header
    integer :: field_format_recs
    integer, parameter :: iunit=10
    integer :: ier1, ios
    integer :: n, n1, n2, nf, nrecs
    integer :: obs_num_mass  ! number of mass (kx values) subtypes to process
    integer :: obs_num_wind  ! number of wind (kx values) subtypes to process
!
    real(rkind1) :: area_average_diameter ! avg fields within circle of diameter (km)
    character(len=30)  :: read_name
    character(len=30)  :: cdum
    character(len=260) :: read_info(read_max)
    character(len=260) :: akbk_filename
    character(len=1)   :: mysub0
!
    ier=0  ! error counter
!
! Set msub0 to eliminate error handeling in sub=find_name_2; Instead errors 
! are handeled here so that defaults are set without error indicated 
    mysub0(1:1)=' '
!
! Set some default values that may, however, not be used for all data types
    field_ranseed=1111
    field_satwind_types='TTTT'
    field_stride_i=1
    field_stride_j=1
    field_obs_types_num=1
    field_obs_types(:)=0
    not_rad= trim(dtype) == 'SATWIND' .or.  trim(dtype) == 'CONV' .or. &
             trim(dtype) == 'GRIDPT' .or.  trim(dtype) == 'GPSRO' 
    field_conv_ksmooth(:)=0  ! default values
!
! Open .rc file (exit routine if problem opening file)
    open (iunit,file=trim(rcfile),form='formatted',status='old',iostat=ios)
    if (ios /= 0) then 
      ier=99
      print *,' '
      print ('(a,i3,a,i4,2a)'),' ERROR attempting to open Rc file for iunit=', &
                     iunit,' iostat=',ios,' and file name=',trim(rcfile)
      return  
    elseif (lprint) then
      print *,' '
      print ('(3a,i4)'),' Rc file=',trim(rcfile), &
              ' opened on unit=',iunit
    endif
!  
! read header format indicator
    read (iunit,*) cdum,field_format_header,field_format_recs
!
! skip file descriptor records
    do n1=1,20  
      read (iunit,'(a)') read_name
      if (read_name(1:1) == '#') exit
    enddo
!
! Read header info
    do n1=1,read_max
      read (iunit,'(a)') read_info(n1)
      if (read_info(n1)(1:3) == '#') exit
      nrecs=n1
    enddo     
!
! Look for value of field_common_path for files of NR fields
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_common_path',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (lprint) then
        print *,' '
        print *,'ERROR reading RC file: ', &
                'field_common_path not found in common portion of header'
      endif
      ier=ier+1
    else 
      read (read_info(n),'(a18,a)') cdum,field_common_path
    endif
!
! Look for value of field_dimensions
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_dimensions',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (lprint) then
        print *,' '
        print *,'ERROR reading RC file: ', &
               'field_dimensions not found in common portion of header'
      endif
      ier=ier+1
    else 
      read (read_info(n),*) cdum,field_imax,field_jmax,field_kmax
    endif
!
! Look for value of lon_first
    call find_name_2 (nrecs,read_info,lstop,mysub0,'lon_first',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (lprint) then
        print *,' '
        print *,'ERROR reading RC file: ', &
                'lon_first not found in common portion of header'
      endif
      ier=ier+1
    else 
      read (read_info(n),*) cdum,field_lon_first
    endif
!
! Look for value of field_time_slots
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_time_slots',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (not_rad) then     ! this is a required variable
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'field_time_slots not found in common portion of header'
        endif
        ier=ier+1
      else              ! variable not required; set to default non-use flag
        field_time_slots=-99
      endif
    else 
      read (read_info(n),*) cdum,field_time_slots
    endif
!
! Look for value of field_time_delta
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_time_delta',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (not_rad) then     ! this is a required variable
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'field_time_delta not found in common portion of header'
        endif
        ier=ier+1
      else             ! variable not required; set to default non-use flag
        field_time_delta=-99.00
      endif
    else 
      read (read_info(n),*) cdum,field_time_delta
    endif
!
! Look for value of field_time_first
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_time_first',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (not_rad) then     ! this is a required variable
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'field_time_first not found in common portion of header'
        endif
        ier=ier+1
      else             ! variable not required; set to default non-use flag
        field_time_first=-99.00
      endif
    else 
      read (read_info(n),*) cdum,field_time_first
    endif
!
! Look for value of field_top
    call find_name_2 (nrecs,read_info,lstop,mysub0,'field_top(Pa)',n)
    if (n == 0) then        ! named variable not found in list in rc file
      field_top=0.          ! default value in Pa
      if (lprint) then
        print *,' '
        print *,'field_top not found in common portion of header'
        print ('(a,f7.1)'),'Default value set as field_top(Pa)=',field_top
      endif
    else 
      read (read_info(n),*) cdum,field_top
    endif
!
! Look for value of tavg_offset
    call find_name_2 (nrecs,read_info,lstop,mysub0,'tavg_offset',n)
    if (n == 0) then        ! named variable not found in list in rc file
      field_tavg_offset=0.  ! default value = no offset
      if (lprint) then
        print *,' '
        print *,'tavg_offset not found in common portion of header'
        print ('(a,f7.1)'),'Default value set as field_tavg_offset=', &
              field_tavg_offset
      endif
    else 
      read (read_info(n),*) cdum,field_tavg_offset
    endif
!
! Look for value of akbk_filename
    call find_name_2 (nrecs,read_info,lstop,mysub0,'akbk_filename',n)
    if (n == 0) then        ! named variable not found in list in rc file
      if (lprint) then
        print *,' '
        print *,'ERROR reading RC file: ', &
                'akbk_filename not found in common portion of header'
      endif
      ier=ier+1
    else 
      read (read_info(n),'(a14,a240)') cdum,akbk_filename
    endif
!
! Special requirements for SATWIND
!
    if (trim(dtype) == 'SATWIND') then ! Read info required only for Satwinds
!
! Look for value of random_seed 
      call find_name_2 (nrecs,read_info,lstop,mysub0,'random_seed',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'WARNING reading RC file: ', & 
                  'random_seed not found in common portion of header'
          print ('(a,i8)'),'Default vale set as field_ranseed=',field_ranseed
        endif
      else 
        read (read_info(n),*) cdum,field_ranseed
      endif
!
! Look for value of field_satwind_types
      call find_name_2 (nrecs,read_info,lstop,mysub0,'obs_dtypes',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'WARNING reading RC file: ', &
                  'obs_dtypes not found in common portion of header'
          print ('(a,i8)'),'Default value set as field_satwind_types=', &
                            field_satwind_types
        endif
      else 
        read (read_info(n),*) cdum,field_satwind_types
      endif
!
    endif  ! test on dtype=SATWIND
!
! Look for value of field_stride_ij 
    if (trim(dtype) == 'SATWIND' .or. trim(dtype) == 'GRIDPT' ) then
      call find_name_2 (nrecs,read_info,lstop,mysub0,'field_stride_ij',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'field_stride_ij not found in common portion of header'
        endif
        ier=ier+1
      else 
        read (read_info(n),*) cdum,field_stride_i,field_stride_j
      endif
    endif    ! test on dtype
!
! Special requirements for CONV
!
    if (trim(dtype) == 'CONV') then 
!
! Read in field profile smoothing parameters if present
      call find_name_2 (nrecs,read_info,lstop,mysub0,'conv_ksmooth',n)
      if (n /= 0) then            ! named variable found in list in rc file
        read (read_info(n),*) cdum,field_conv_ksmooth(:)
      endif
!
! Look for list of subtypes (kx values) to process
      call find_name_2 (nrecs,read_info,lstop,mysub0,'obs_types_num2',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'obs_types_num2 not found in common portion of header'
        endif
        ier=ier+1
      else 
        read (read_info(n),*) cdum,obs_num_mass,obs_num_wind
      endif
!
      field_obs_types(:)=0
      call find_name_2 (nrecs,read_info,lstop,mysub0,'obs_types_mass',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'obs_types_mass not found in common portion of header'
        endif
        ier=ier+1
      else 
        read (read_info(n),*) cdum,field_obs_types(1:obs_num_mass)
      endif
!
      n1=obs_num_mass+1
      n2=obs_num_mass+obs_num_wind
      field_obs_types_num=n2
      call find_name_2 (nrecs,read_info,lstop,mysub0,'obs_types_wind',n)
      if (n == 0) then        ! named variable not found in list in rc file
        if (lprint) then
          print *,' '
          print *,'ERROR reading RC file: ', &
                  'obs_types_wind not found in common portion of header'
        endif
        ier=ier+1
      else 
        read (read_info(n),*) cdum,field_obs_types(n1:n2)
      endif
!
    endif  ! test on dtype=CONV
!
! Special requirements for GPSRO
!
    if (trim(dtype) == 'GPSRO') then ! Read info required only for GPSRO
!
! Read in field profile smoothing parameters if present
      call find_name_2 (nrecs,read_info,lstop,mysub0,'gpsro_ksmooth',n)
      if (n == 0) then            ! named variable not found in list in rc file
        field_gpsro_ksmooth(:)=0  ! default is no vertical smoothing
      else 
        read (read_info(n),*) cdum,field_gpsro_ksmooth(:)
      endif
!
! Read in indicator for selection of GPSRO algorithm
      call find_name_2 (nrecs,read_info,lstop,mysub0,'gpsro_algorithm',n)
      if (n == 0) then           ! named variable not found in list in rc file
        field_gpsro_algorithm=0  ! default is 2-D ROPP
      else 
        read (read_info(n),*) cdum,field_gpsro_algorithm
      endif
!
    endif  ! test on dtype=GPSRO
!
!
! continue reading rc file
!
! find and read information for requested radiance type in rc file
    read (iunit,*) cdum  ! skip spacer record
    read (iunit,*) cdum  ! skip spacer record
    found=.false.
    do n1=1,10000
      read (iunit,'(a)') read_name
      if (trim(read_name) == 'EOF') exit 
      if (trim(read_name) == trim(dtype)) then
        found=.true.
        nrecs=0
        do n2=1,100 ! read info for found data type
          read (iunit,'(a)') read_info(n2)
          if (read_info(n2)(1:3) == '---') exit
          nrecs=n2
        enddo
        exit
      endif
    enddo
! 
    if (.not. found) then
      ier=ier+1
      if (lprint) then 
        print *,'field file info for ',trim(dtype),' not found on .rc file'
      endif
    endif
!
    close (iunit)
!
! read file of a(k) and b(k) values if akbk_filename provided.
    if (akbk_filename /= 'none') then
      call read_akbk (field_max_akbk,field_kmax,iunit,akbk_filename, &
                      lprint,field_akbk,ier1)
      ier=ier+ier1
    endif
!      
! Check if any necessary values are missing 
    if (ier /= 0) then
      if (lprint) then 
        print *,' '
        print ('(i3,a)'),ier,' required variables on rc file not found'
        print ('(a)'),' Attempt to read rc file further stopped'
      endif
      return
    endif
!
! Set parameter (if present) for how fields are to be area averaged
    n=1  ! record index
    if (read_info(n)(1:12) == 'area_average') then
      read (read_info(n),*) cdum,area_average_diameter
      field_avg(:,:)=area_average_diameter  ! km 
      n=n+1 
    else
      field_avg(:,:)=0.0  ! Do not use area averages for field profiles
    endif
!
! Interpret info for 2-d fields required for profiles
    read (read_info(n),*) cdum,field_num_2d 
    do nf=1,field_num_2d
      n=n+1
      read (read_info(n),'(2a12,i3,2x,a)') field_names(nf,1,2),      &
              field_names(nf,2,2),field_types(nf,2),field_files(nf,2)
    enddo
!
! Interpret info for 3-d fields required for profiles
    n=n+1 
    read (read_info(n),*) cdum   ! skip spacer
    n=n+1 
    read (read_info(n),*) cdum,field_num_3d
    do nf=1,field_num_3d
      n=n+1
      read (read_info(n),'(2a12,i3,2x,a)') field_names(nf,1,3),      &
              field_names(nf,2,3),field_types(nf,3),field_files(nf,3)
    enddo
!
    end subroutine nr_fields_read_rc 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine nr_fields_print_info (dtype,call_loc)
!
!  Print some information extracted from the field list .rc file. 
!
    implicit none
    character(len=*), intent(in) :: dtype
    character(len=*), intent(in) :: call_loc  ! name of calling routine
!    
    integer :: n
!
    print *,' '
    print ('(2a)'),'Information regarding required fields: ',trim(call_loc)
    print ('(a)'),'(-99 indicates a non-used variable for this data type)' 
    print ('(a)'),'(non-specified required variables may be set to default)' 
    print ('(a,3i6,a,f12.4)'),'Grid dimensions:',field_imax,field_jmax, &
            field_kmax,'  first longitude on grid =',field_lon_first
    print ('(a,2i4)'),'field_stride_i,j =',field_stride_i,field_stride_j
    print ('(a,i3,a,f8.0)'),'field_lev_first=',field_lev_first, &
                            '  field_top=',field_top
    print ('(a,i3,2f8.2)'),'Time slots: number, delta, and time-first =', &
            field_time_slots, field_time_delta, field_time_first
    print ('(a,f8.2,2x,a)'),'field_tavg_offset=',field_tavg_offset, &
                            '(0.=default)'
    print ('(2a)'),'field_common_path=',trim(field_common_path)
!
    if (trim(dtype) == 'SATWIND') then 
      print ('(a,2f8.0,2x,a,i10)'),'field_ipw_top,_bot,_types,_ranseed=', &
            field_ipw_top,field_ipw_bot,field_satwind_types,field_ranseed
    elseif (trim(dtype) == 'CONV') then 
      print ('(a,60i4)'),'field_obs_types=', &
           field_obs_types(1:field_obs_types_num)
      print ('(a,3i4)'),'field_conv_ksmooth=',field_conv_ksmooth(:) 
    elseif (trim(dtype) == 'GPSRO') then
      print ('(a,3i4)'),'field_gpsro_ksmooth=',field_gpsro_ksmooth(:) 
      print ('(a,i2)'),'field_gpsro_algorithm=',field_gpsro_algorithm
    endif
!
    print *,' '
    print ('(a,2i3)'),'Info for 2d fields in profiles:'
    do n=1,field_num_2d

      print ('(i2,2(2x,a12),i4,f5.1,2x,a)'),n,trim(field_names(n,1,2)), &
             trim(field_names(n,2,2)),field_types(n,2),field_avg(n,2),  &
             trim(field_files(n,2)) 
    enddo
    print ('(a,2i3)'),'Info for 3d fields in profiles:'
    do n=1,field_num_3d
      print ('(i2,2(2x,a12),i4,f5.1,2x,a)'),n,trim(field_names(n,1,3)), &
             trim(field_names(n,2,3)),field_types(n,3),field_avg(n,3),  &
             trim(field_files(n,3)) 
    enddo
!
    end subroutine nr_fields_print_info
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    end module m_nr_fields_info
