    module m_rad_thin_flds
!
!  Module for reading fields used for observation thinning.
!
!  Initial code Ronald Errico July 2014 
!
    use m_kinds, only : rkind1
!
    implicit none
!
    private
!
! public routines
    public :: rad_thin_flds_setup
    public :: rad_thin_flds_clean
    public :: rad_thin_flds_interp
!
! private routines
    private :: rad_thin_flds_rc
    private :: rad_thin_flds_read
    private :: rad_thin_flds_print_info
!
! public variables
    public :: fld_imax, fld_jmax, fld_lon_first 
    public :: fld_thin_box_size, fld_thin_box_subtypes, fld_thin_box_ids
    public :: fld_time_slots, fld_time_delta, fld_time_first
    public :: fld_max_num, fld_num_tvar, fld_num_tcon, fld_num_all
    public :: fld_common_path, fld_names, fld_files, fld_types
    public :: fld_weights, fld_tavg_offset, fld_avg
!     
! Information describing fields needed to construct profiles
! Last index on _nums, _names, _files indicates (2) 2d or (3) 3d
! First index on _names is (1) name used by osse sim software, (2) name on file 
    integer, parameter :: fld_max_num=5
    integer :: fld_thin_box_subtypes
    integer :: fld_thin_box_ids(100)
    integer :: fld_imax
    integer :: fld_jmax
    integer :: fld_time_slots  ! number of time slots in assim. window
    integer :: fld_num_tcon    ! number of temporally constant 2d fields to use
    integer :: fld_num_tvar    ! number of time-varying 2d fields to be used
    integer :: fld_num_all     ! sum of fld_num_tvar + fld_num_tcon
    integer :: fld_types(fld_max_num,2) ! see below
    real(rkind1) :: fld_thin_box_size   ! thinning box width (km)
    real(rkind1) :: fld_time_delta   ! hours between files  
    real(rkind1) :: fld_time_first   ! hours for 1st file relat to assim time 
    real(rkind1) :: fld_lon_first    ! longitude for index i=1
    real(rkind1) :: fld_tavg_offset  ! hours offset for tavg files requested
    real(rkind1) :: fld_weights(fld_max_num,2) ! weights for thinning
    real(rkind1) :: fld_avg(fld_max_num,2)     ! avg flds over box this size (km) 
    real, allocatable :: fld_tcon(:,:,:)       ! time-independent NR fields
    real, allocatable :: fld_tvar(:,:,:,:)     ! time-dependent NR fields
    character(len=240) :: fld_common_path
    character(len=16)  :: fld_names(2,fld_max_num,2)
    character(len=240) :: fld_files(fld_max_num,2)
!
! fld_types = 1 if field is inst, with horiz and temporal interp requested
! fld_types = 2 if field is inst, with no horiz or temporal interp requested
! fld_types = 3 if field is tavg, with horiz interp requested
! fld_types = 4 if field is tavg, with no horiz interp requested
! fld_types = 1 or 2 if field const; no temporal interp done
!
    character(len=*),parameter :: myname='m_rad_thin_flds'
!
    contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    subroutine rad_thin_flds_setup (dtype,rcfile,cdtime,lstop,lprint,ier)
!
! Read and print info about fields to use for computing thinning penalty.
! Also allocate arrays to hold required NR fields. 
!
    implicit none   
!    
    logical, intent(in) :: lstop
    logical, intent(in) :: lprint
    character(len=*), intent(in) :: dtype
    character(len=*), intent(in) :: rcfile 
    character(len=*), intent(in) :: cdtime
    integer, intent(out) :: ier
!   
    integer :: ier1, ier2, ier3
    character(len=*),parameter :: myname_sub=myname//'::rad_thin_flds_setup'
!
    call rad_thin_flds_rc (dtype,rcfile,lprint,ier)
    if (ier /= 0) return
!
    fld_num_all=fld_num_tcon+fld_num_tvar
    if (lprint) then 
      call rad_thin_flds_print_info (myname_sub)
    endif
!
    if (fld_num_tcon > 0) then
      allocate (fld_tcon(fld_imax,fld_jmax,fld_num_tcon))
    endif    
    if (fld_num_tvar > 0) then
      allocate (fld_tvar(fld_imax,fld_jmax,fld_time_slots,fld_num_tvar))
    endif
    if (fld_num_all > 0) then
      call rad_thin_flds_read (cdtime,lstop,lprint,ier)
    endif
!
    end subroutine rad_thin_flds_setup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine rad_thin_flds_clean
!
! Deallocate arrays created by rad_thin_flds_setup
!
    if (fld_num_tcon > 0) then 
      deallocate (fld_tcon)
    endif
    if (fld_num_tvar > 0) then 
      deallocate (fld_tvar)
    endif
!
    end subroutine rad_thin_flds_clean
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine rad_thin_flds_rc (dtype,rcfile,lprint,ier)
!
! Routine to read .rc file for radiance field requirements for creating
! (possibly) thinned list of obs
!
    use m_rad_obs_arrays, only : obs_time_slots, obs_time_delta
    use m_rad_obs_arrays, only : obs_time_first
    use m_set_unit_nums, only : un_info
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
    logical :: found
    integer :: fld_format_header
    integer :: fld_format_recs
    integer, parameter :: iunit=un_info
    integer :: n1, n2, nrecs, nf, n, ios 
!
    real(rkind1) :: area_average_diameter ! avg fields within circle of diameter (km)
    character(len=30)  :: read_name
    character(len=260) :: read_info(100)
    character(len=30)  :: cdum
!
    ier=0  ! error counter
!
    open (iunit,file=trim(rcfile),form='formatted',status='old',iostat=ios)
    if (ios /= 0) then 
      ier=99
      print *,' '
      print ('(2a,i3,a,i4,2a)'),' ERROR attempting to open rad_thin_flds Rc file', &
                   ' for iunit=',iunit,' iostat=',ios,' and file name=',trim(rcfile)
      return
    elseif (lprint) then
      print *,' '
      print ('(3a,i4)'),' Rc file=',trim(rcfile),' opened on unit=',iunit
    endif
!  
! skip header records 
    read (iunit,*) cdum,fld_format_header,fld_format_recs
!
    do n1=1,20  
      read (iunit,'(a)') read_name
      if (read_name(1:1) == '#') exit
    enddo
!
! read directory and grid information for fields
!
    read (iunit,'(a18,a240)') cdum,fld_common_path
    if (trim(cdum) /= 'field_common_path') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not= field_common_path'
    endif
    read (iunit,*) cdum,fld_imax,fld_jmax
    if (trim(cdum) /= 'field_dimensions') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not= field_dimensions'
    endif
    read (iunit,*) cdum,fld_lon_first
    if (trim(cdum) /= 'lon_first') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = lon_first'
    endif
    read (iunit,*) cdum,obs_time_slots
    if (trim(cdum) /= 'obs_time_slots') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = obs_time_slots'
    endif
    read (iunit,*) cdum,obs_time_delta
    if (trim(cdum) /= 'obs_time_delta') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = obs_time_delta'
    endif
    read (iunit,*) cdum,obs_time_first
    if (trim(cdum) /= 'obs_time_first') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = obs_time_first'
    endif
    read (iunit,*) cdum,fld_time_slots
    if (trim(cdum) /= 'field_time_slots') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = field_time_slots'
    endif
    read (iunit,*) cdum,fld_time_delta
    if (trim(cdum) /= 'field_time_delta') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = field_time_delta'
    endif
    read (iunit,*) cdum,fld_time_first
    if (trim(cdum) /= 'field_time_first') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = field_time_first'
    endif
    read (iunit,*) cdum,fld_tavg_offset
    if (trim(cdum) /= 'tavg_offset') then
      ier=ier+1
      print *,'cdum=',trim(cdum),' not = tavg_offset'
    endif
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
    close (iunit)
!
    if (.not. found) then
      ier=ier+1
      print *,'field info for ',trim(dtype),' not found'
!      
    else
!
! Interpret info for thinning box size
! fld_thin_box_size < 1. indicates no thinning will be performed
      n=1
      read (read_info(n),*) cdum,fld_thin_box_size
      n=n+1
      read (read_info(n),*) cdum,fld_thin_box_subtypes
      n=n+1
      read (read_info(n),*) cdum,fld_thin_box_ids(1:fld_thin_box_subtypes)
!
! Interpret info for 2-d time-invariant fields required 
      n=n+1    ! record index
      nf=0   ! field index
!
!
! Set parameter (if present) for how fields are to be area averaged
      if (read_info(n)(1:12) == 'area_average') then
        read (read_info(n),*) cdum,area_average_diameter
        fld_avg(:,:)=area_average_diameter  ! km 
        n=n+1 
      else
        fld_avg(:,:)=0.0  ! Do not use area averages for field profiles
      endif
!
      read (read_info(n),*) cdum,fld_num_tcon
      do n2=1,fld_num_tcon
        nf=nf+1 
        n=n+1
        read (read_info(n),'(2a16,i3,f7.2,2x,a240)') fld_names(:,nf,1), &
                 fld_types(nf,1),fld_weights(nf,1),fld_files(nf,1)
      enddo
!
! Interpret info for 2-d time-varying fields required 
      n=n+1 
      nf=0
      read (read_info(n),*) cdum,fld_num_tvar
      do n2=1,fld_num_tvar
        nf=nf+1 
        n=n+1
        read (read_info(n),'(2a16,i3,f7.2,2x,a240)') fld_names(:,nf,2), &
                 fld_types(nf,2),fld_weights(nf,2),fld_files(nf,2)
      enddo
!
    endif ! check on found
!
    end subroutine rad_thin_flds_rc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine rad_thin_flds_print_info (call_loc)
!
! Print some info read from the resource files 
!
    use m_rad_obs_arrays, only : obs_time_slots, obs_time_delta
    use m_rad_obs_arrays, only : obs_time_first
    implicit none
    character(len=*) :: call_loc
!    
    integer :: n
!
    print *,' '
    print ('(2a)'),'Information regarding required fields: ',trim(call_loc)
    print ('(a,2i6,a,f12.4)'),'Grid dimensions:',fld_imax,fld_jmax, &
           '  first longitude on grid =',fld_lon_first
    print ('(2a,i3,3f8.2)'),'Time slots for fields: ',  &
            'number, delta, first, tavg_offset =',      &
            fld_time_slots,fld_time_delta,fld_time_first,fld_tavg_offset
    print ('(2a,i3,2f8.2)'),'Time slots for sorted obs: ',  &
            'number, delta, and first =',                   &
            obs_time_slots,obs_time_delta,obs_time_first
    print ('(a,10i8)'),'Obs subset IDs: ', &
            fld_thin_box_ids(1:fld_thin_box_subtypes)
    print *,'Field_common_path=',trim(fld_common_path)
    print ('(a,2i3)'),'Info for 2d time constant fields required:'
    do n=1,fld_num_tcon
      print ('(i2,2(2x,a10),i4,2f7.2,2x,a)'),n,trim(fld_names(1,n,1)), &
             trim(fld_names(2,n,1)),fld_types(n,1),fld_weights(n,1),   &
             fld_avg(n,1),trim(fld_files(n,1)) 
    enddo
    print ('(a,2i3)'),'Info for 2d time varying fields required:'
    do n=1,fld_num_tvar
      print ('(i2,2(2x,a10),i4,2f7.2,2x,a)'),n,trim(fld_names(1,n,2)), &
             trim(fld_names(2,n,2)),fld_types(n,2),fld_weights(n,2),   &
             fld_avg(n,2),trim(fld_files(n,2)) 
    enddo
!
    end subroutine rad_thin_flds_print_info
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
    subroutine rad_thin_flds_read (cdtime_ref,lstop,lprint,ier)
!
! Sequence commands to read NR fields required to compute thinning penalty 
!
    use m_time_compute, only : time_compute_new_cdtime
    implicit none 
!
    logical, intent(in) :: lstop
    logical, intent(in) :: lprint
    character(len=14), intent(in) :: cdtime_ref
    integer, intent(out) :: ier
!
    integer :: nf,nt
    integer :: ier1, ier2
    real(rkind1) :: dhours
    character(len=14) :: cdtime1
    character(len=240) :: file_name
!
    ier=0
!
    do nf=1,fld_num_tcon
      call set_field_file_name (fld_names(2,nf,1),fld_files(nf,1),           &
                                fld_common_path,cdtime_ref,file_name,ier) 
      call read_nc4_2dfield (fld_imax,fld_jmax,file_name,fld_names(2,nf,1),  &
                             fld_tcon(:,:,nf),lstop,lprint,ier) 
    enddo
!
    do nt=1,fld_time_slots
      do nf=1,fld_num_tvar
        dhours=fld_time_first+(nt-1)*fld_time_delta
        if (fld_types(nf,2) > 2) then
          dhours=dhours+fld_tavg_offset
        endif 
        call time_compute_new_cdtime (cdtime_ref,cdtime1,dhours,ier)
        call set_field_file_name (fld_names(2,nf,2),fld_files(nf,2),          &
                                  fld_common_path,cdtime1,file_name,ier) 
        call read_nc4_2dfield (fld_imax,fld_jmax,file_name,fld_names(2,nf,2), &
                               fld_tvar(:,:,nt,nf),lstop,lprint,ier) 
      enddo
    enddo
!
    end subroutine rad_thin_flds_read
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    subroutine rad_thin_flds_interp (obs_lat,obs_lon,obs_time,nflds,f_int)
!
! Interpolate all time varying 2d fields of data horizontally and temporally
! Interpolate all time constant 2d fields horizontally
!
    implicit none
!
    integer, intent(in)       :: nflds     ! total number of field values to get
    real(rkind1), intent(in)  :: obs_lat   ! latitude of observation 
    real(rkind1), intent(in)  :: obs_lon   ! longitude of observation 
    real(rkind1), intent(in)  :: obs_time  ! longitude of observation 
    real(rkind1), intent(out) :: f_int(nflds)
!
    integer, parameter :: n_area_dim=100
    integer :: ij_area(2,n_area_dim)
    integer :: h_index(2,2)
    integer :: t_index
    integer :: n_area_max,n_area
    integer :: ipa,jpa,k
    integer :: ka
    integer :: nt, nt1, n
    real(rkind1) :: fh(2)
    real(rkind1) :: td1, td2
    real(rkind1) :: h_weights(2,2)        
    real(rkind1) :: hw(2,2)        
    real(rkind1) :: w_area(2,n_area_dim)
    real(rkind1) :: time_weights(2)        
    real(rkind1) :: tw(2)        
!
   character(len=*),parameter :: myname_sub=myname//'::rad_thin_flds_interp'
!
! Determine grid box location that surrounds observation location
! and the interpolation weights for each such point
    call get_interp_horiz_index (fld_imax,fld_jmax,fld_lon_first, &
                                 obs_lat,obs_lon,h_index,h_weights)
!
! First consider temporally constant fields 
    ka=0 ! field counter: time-constant fields followed by time-varying fields
    do k=1,fld_num_tcon
      ka=ka+1
!
! Determine whether to use interpolated values (use of closest 
! point is included in this option) or to use area averaged values
      if (fld_avg(k,1) < 1.01 .or. mod(fld_types(k,1),2) == 0) then      
!
! If indicated by field_type, use only nearest point horizonally 
        hw(:,:)=h_weights(:,:) 
        if (mod(fld_types(k,1),2) == 0) then  
          call get_interp_horiz_nearest (hw)
        endif
!
        f_int(ka)=hw(1,1)*fld_tcon(h_index(1,1),h_index(1,2),k) &
                 +hw(1,2)*fld_tcon(h_index(1,1),h_index(2,2),k) &
                 +hw(2,1)*fld_tcon(h_index(2,1),h_index(1,2),k) &
                 +hw(2,2)*fld_tcon(h_index(2,1),h_index(2,2),k) 
!
      else
!
! Define values as a local sample area mean 
        call integrate_area_params (fld_imax,fld_jmax,n_area_dim, &
                                    obs_lat,obs_lon,fld_lon_first,  &
                                    fld_avg(k,1),n_area_max,ij_area,w_area)
        f_int(ka)=0.
        do n_area=1,n_area_max
          ipa=ij_area(1,n_area)
          jpa=ij_area(2,n_area)
          f_int(ka)=f_int(ka)+w_area(1,n_area)*fld_tcon(ipa,jpa,k)
        enddo              
!
      endif ! check on whether to compute area average
!
    enddo
!
! Next consider temporally varying fields 
    if (fld_num_tvar > 0) then
!
! Determine weights for 2 times
      td1=obs_time-fld_time_first    ! time between obs time and first file time
      nt1=1+int(td1/fld_time_delta)  ! file time slot pair in which obs occurs
      nt1=min(nt1,fld_time_slots-1)
!
      td2=nt1*fld_time_delta-td1    ! time between obs time and next file time
      time_weights(1)=td2/fld_time_delta
      time_weights(2)=1._rkind1-time_weights(1)
!
      do k=1,fld_num_tvar
        ka=ka+1
!
! If indicated by field_type, use only nearest point temporally
        tw(:)=time_weights(:) 
        if (fld_types(k,2) == 2) then   
          call get_interp_time_nearest (tw)
        elseif (fld_types(k,2) > 2) then  ! use value from tavg file
          tw(1)=1.
          tw(2)=0.
        endif   
!
! Determine whether to use interpolated values (use of closest 
! point is included in this option) or to use area averaged values
        if (fld_avg(k,2) < 1.0 .or. mod(fld_types(k,2),2) == 0) then      
!
! Use interpolation. If indicated by field_type, use only nearest 
! point horizonally and assign weight to only a single point. For 
! interpolation consider 4 surrounding points: SW, NW, SE, NE
          hw(:,:)=h_weights(:,:) 
          if (mod(fld_types(k,2),2) == 0) then  
            call get_interp_horiz_nearest (hw)  
          endif
!
! interpolate horizontally at 2 times 
          do n=1,2
            nt=nt1+(n-1)
            fh(n)=hw(1,1)*fld_tvar(h_index(1,1),h_index(1,2),nt,k) &
                 +hw(1,2)*fld_tvar(h_index(1,1),h_index(2,2),nt,k) &
                 +hw(2,1)*fld_tvar(h_index(2,1),h_index(1,2),nt,k) &
                 +hw(2,2)*fld_tvar(h_index(2,1),h_index(2,2),nt,k) 
          enddo
!
        else
!
! Define values as a local sample area mean 
          call integrate_area_params (fld_imax,fld_jmax,n_area_dim,   &
                                      obs_lat,obs_lon,fld_lon_first,  &
                                      fld_avg(k,2),n_area_max,ij_area,w_area)
          fh(:)=0.
          do n_area=1,n_area_max
            ipa=ij_area(1,n_area)
            jpa=ij_area(2,n_area)
            do n=1,2
              nt=nt1+(n-1)
              fh(n)=fh(n)+w_area(1,n_area)*fld_tvar(ipa,jpa,nt,k)
            enddo
          enddo              
!
        endif ! check on whether to compute area average
!
! Interpolate in time
        f_int(ka)=tw(1)*fh(1)+tw(2)*fh(2)
!
      enddo ! loop over time-varying field indexes
    endif   ! consideration of time-varying fields
!
    if (ka /= fld_num_all) then
      print *,'Error in ',myname_sub
      print *,' ka /= fld_num_all ',ka,fld_num_all
      stop
    endif 
!
    end subroutine rad_thin_flds_interp
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
    end module m_rad_thin_flds
