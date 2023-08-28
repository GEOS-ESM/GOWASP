    program eparams_rad_hcorr
!
! Program to estimate obs error table parameters for radiance obs
! Includes parameters for stdv, fraction of var of correlated error, 
! and corrleation lengths.
!
    use m_kinds, only : rkind1
!
    use m_sat_info_table, only : sat_info_table_read
    use m_sat_info_table, only : sat_info_table_get_2c
    use m_sat_info_table, only : sat_info_table_get_1i
!
    use m_obs_error_table, only : error_table_setup
    use m_obs_error_table, only : error_table_clean
    use m_obs_error_table, only : error_table_find_stdv_id
    use m_obs_error_table, only : error_table_find_corr_id
    use m_obs_error_table, only : et_pert_fac
    use m_obs_error_table, only : et_err_tab, et_err_itype
    use m_obs_error_table, only : et_l_vc_corr, et_nlevels    
    use m_obs_error_table, only : et_cov_matrix
    use m_obs_error_table, only : et_hcorr_lengths, et_frac_corr
!
    implicit none
!
!  In the following,
!   _t refers to a target value (generally values from a DAS of real obs)
!   _o refers to an osse value used by or obtained from a prev DAS osse result
!   _d refers to previously used or here-estimated added error statistic
!   -n refers to a new estimated statistic for errors to add
!
    logical :: lcorr_prev  ! previous added errors were correlated
    logical :: lcorr_both  ! correl for both traget and osse avialable 
    logical :: lerrtable   ! indicates that prev table of added stdv available
    logical, parameter :: lwrite_corr=.true.  ! T if corrs written to 2nd file 
!
    integer, parameter :: rkind=8   ! precision (8 because of variable in file)
    integer, parameter :: min_count_cov=1 ! min count to consider cov realiable   
    integer, parameter :: iunit=11          ! unit number for output
    integer :: k,n,j
    integer :: islash_t,islash_o ! character position of last \ in file names
    integer :: argc
    integer :: nk         ! number of channels or p-levels
    integer :: ierr       ! return code from subroutines
    integer :: id_type    ! id for desired obs subtype in extracted error table
    integer :: hc_type    ! id for desired obs subtype in file of corr params
    integer :: itype      ! either sat or inst id number to match 
    integer :: nbins      ! min (nbins_t, nbins_o)
    integer :: nbins_t    ! number of bins of horiz. distances for target stats
    integer :: nbins_o    ! number of bins of horiz. distance for OSSE stats
    integer :: ntimes_o   ! number of assimilation periods used to comp stats
    integer :: ntimes_t   ! number of assimilation periods used to comp stats
    integer :: ncnt
    integer :: idum2(2) 
    integer :: count_corr_fix(2)  ! number of correlations adjusted
    integer(4) :: iargc
    integer, allocatable :: count_t(:)  ! mean num of real obs per assim time
    integer, allocatable :: count_o(:)  ! mean num of osse obs per assim time
    integer, allocatable :: count_d(:)  ! min(count_t, count_o)
!
    real(rkind1), parameter :: corr_max=0.95 ! max allowed inter channel corr
    real(rkind1) :: pert_fac2  ! square of err pert_fac used in prev osse
    real(rkind1) :: dvar       
    real(rkind1) :: xfac,xfac1,xfac2
    real(rkind1) :: xmin_count ! real value of min_count_var
    real(rkind1) :: xdum(4)
    real(rkind1), allocatable :: frac_corr(:) ! frac of added var that is correl
    real(rkind1), allocatable :: plevs_t(:) ! used as dummy argument
    real(rkind1), allocatable :: plevs_o(:) ! used as dummy argument
    real(rkind1), allocatable :: Var_t(:)   ! var of target innovations     
    real(rkind1), allocatable :: Var_o(:)   ! var of prev OSSE innovations
    real(rkind1), allocatable :: Rsqrt_o(:) ! stdv from prev OSSE err table 
    real(rkind1), allocatable :: Rsqrt_n(:) ! stdv for new OSSE error table
    real(rkind1), allocatable :: Rvar_o(:)  ! var of prev OSSE added errors 
    real(rkind1), allocatable :: Rvar_n(:)  ! updated estimate for Rvar
    real(rkind1), allocatable :: Bvar_o(:)  ! var of osse obs-space bkg error
    real(rkind1), allocatable :: bfrac_o(:) ! ratio Bvar_e/Var_o
    real(rkind1), allocatable :: C_t(:,:,:) ! cov, corr, bin_counts for target
    real(rkind1), allocatable :: C_o(:,:,:) ! cov, corr, bin_counts prev OSSE
    real(rkind1), allocatable :: C_d(:,:,:) ! cov, corr, bin_counts for estimate
!
    character(len=4), allocatable :: qualflag(:)
    character(len=16)  :: cname(2)
    character(len=16)  :: dtype
    character(len=16)  :: sat_name
    character(len=16)  :: instr_name
    character(len=16)  :: c_sat_name
    character(len=240) :: file_target_stats
    character(len=240) :: file_osse_stats
    character(len=240) :: file_error_rc
    character(len=240) :: file_sat_info
    character(len=240) :: file_eparams_new
    character(len=240) :: file_cov_new
!
    print *,' '
    print *,'BEGIN NEW PROGRAM estimate_eparams_rad_chcorr'
!
! Read arguments
    argc = iargc()
    if (argc .ne. 8) then
      print *,'useage must be: estimate.x '                             
      stop
    endif
    call GetArg( 1_4, instr_name)
    call GetArg( 2_4, c_sat_name)
    call GetArg( 3_4, file_error_rc)
    call GetArg( 4_4, file_target_stats)
    call GetArg( 5_4, file_osse_stats)
    call GetArg( 6_4, file_sat_info)
    call GetArg( 7_4, file_eparams_new)
    call GetArg( 8_4, file_cov_new)
!
    print *,'instrument=',trim(instr_name),'  satellite=',trim(c_sat_name)
!
! Read and save table of sat/instr info
    call sat_info_table_read (file_sat_info,.true.,ierr)
!
! Find dtype, SATNAME, and n_channels for c_sat_name and instr_name  
    call  sat_info_table_get_2c ('instr','sat','dtype',instr_name, &
                                 c_sat_name,dtype,ierr)
    call  sat_info_table_get_2c ('instr','sat','platform',instr_name, &
                                 c_sat_name,sat_name,ierr)
    call  sat_info_table_get_1i ('instr','nchan',instr_name,nk,ierr)
!
    xmin_count=real(min_count_cov)
!
! Allocate arrays whose dimennsions depend only on n_channels
    allocate (qualflag(nk))
    allocate (frac_corr(nk))  
    allocate (plevs_t(nk))  
    allocate (plevs_o(nk))  
    allocate (count_t(nk))
    allocate (count_o(nk))
    allocate (count_d(nk))
    allocate (Var_t(nk))
    allocate (Var_o(nk))
    allocate (Rsqrt_o(nk))
    allocate (Rsqrt_n(nk))
    allocate (Rvar_o(nk))
    allocate (Rvar_n(nk))
    allocate (Bvar_o(nk))
    allocate (bfrac_o(nk))
    allocate (C_t(nk,nk,3))        ! 1-3: cov, corr, bin count measure
    allocate (C_o(nk,nk,3))        ! 1-3: cov, corr, bin count measure 
    allocate (C_d(nk,nk,3))        ! 1-3: cov, corr, bin count measure 
!
print *,'before error.rc'
      !
! Read error.rc file and the files possibly specified within it, including file 
! of prev error stdv table and file of added error correlation parameters.
    if (trim(file_error_rc) == 'none') then 
      lcorr_prev=.false.
      lerrtable=.false.      
    else
      call error_table_setup (.true.,dtype,file_error_rc,lerrtable, &
                              lcorr_prev,idum2(1),ierr)
!
      if (lcorr_prev .and. nk /= et_nlevels) then
        print *,'Error: mismatch between nk in sat table and et_n_levels' 
        print *,'nk=',nk,'  et_n_levels=',et_nlevels
        print *,'Program stopping'
        stop
      endif
! 
   print *,'after error.rc'
! Find id of desired sat/instr in extracted portions of error table 
! First get either sat id or, if sat=AQUA, inst id number.   
      if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then
        call sat_info_table_get_1i ('dtype','siid',dtype,itype,ierr)
      else
        call sat_info_table_get_1i ('platform','said',sat_name,itype,ierr)
      endif
!
      call error_table_find_stdv_id (itype,id_type)
      if (id_type == 0) then
        print *,'ERROR in finding sat/inst in extracted part of errror table:'
        print ('(a,i5,a)'),'itype wanted =',itype,'   itypes in list ='  
        print ('(10i5)'),et_err_itype(:)
        stop
      endif
!
    endif 
!
  print *,'dtype=',dtype
! Determine last character in file path part of file name 
    islash_t=0
    islash_o=0
    do n=1,len(file_target_stats)    
      if (file_target_stats(n:n) == '/' ) islash_t=n
      if (file_osse_stats(n:n) == '/' ) islash_o=n
    enddo
!
! Read target statistics
!
    if (file_target_stats(islash_t+1:islash_t+5) == 'count') then
!
! Get counts and standard deviations of innovations from count_tables
! (No horizontal covariances read)
      call read_counts (nk,ntimes_t,file_target_stats,count_t,plevs_t,Var_t)
      Var_t(:)=Var_t(:)**2  ! replace read stdv by var
      nbins_t=0
      do k=1,nk
        C_t(k,k,3)=count_t(k)*ntimes_t
      enddo
!
    else   
!
! Read binned horizontal counts and covariances of innovations 
! from hcorr_tables
      call read_chcovs (nk,ntimes_t,file_target_stats,ierr,count_t,C_t)
      nbins_t=nk
      if (ierr == 0) then
        do k=1,nk
          Var_t(k)=C_t(k,k,1)
        enddo    
      else
        print *,'ERROR in 1st call to read_chcovs: ierr=',ierr
        stop
      endif 
    endif
!
! Read osse statistics 
!
    if (file_osse_stats(islash_o+1:islash_o+5) == 'count') then
!
! Get counts and standard deviations of innovations from count_tables
! (No horizontal covariances read)
      call read_counts (nk,ntimes_o,file_osse_stats,count_o,plevs_o,Var_o)
      Var_o(:)=Var_o(:)**2  ! replace read stdv by var
      nbins_o=0
      do k=1,nk  
        C_o(k,k,3)=count_o(k)*ntimes_o
      enddo
!
    else   
!
! Read binned horizontal counts and covariances of innovations 
! from hcorr_tables
      call read_chcovs (nk,ntimes_o,file_osse_stats,ierr,count_o,C_o)
      nbins_o=nk
      if (ierr == 0) then
        do k=1,nk
          Var_o(k)=C_o(k,k,1)
        enddo    
      else
        print *,'ERROR in 2nd call to read_chcovs: ierr=',ierr
        stop
      endif 
    endif
!
! Determine if it is possible to compare target and osse covariances
    nbins=min(nbins_t,nbins_o)
    lcorr_both=nbins > 0 
    If (.not. lcorr_both) then 
      print *,' '
      print *,'Correlations cannot be compared because covariance '
      print *,'have not been read for both target and osse:' 
      print ('(a,2i4,2f8.2)'),'nbins_t,nbins_o = ',nbins_t,nbins_o
      lcorr_prev=.false.  ! need not consider previous cov stats
    endif
!
! Assign the min of target and osse mean numbers of obs per DAS cycle as a 
! measure of statistical reliability (larger measure means more reliable).
    if (lcorr_both) then
      do k=1,nk
        do n=1,nbins
          C_d(n,k,3)=min(C_t(n,k,3),C_o(n,k,3))
        enddo
      enddo
    endif
!  
! If previously added errors were correlated, then addjust the difference in 
! covariances so that it reflects what the covariances would be without any 
! and added error.
! Retrieve previous correlation parameters for this data type
    if (lcorr_prev .and. et_l_vc_corr) then 
      call error_table_find_corr_id (itype,hc_type)  
      frac_corr(:)=et_frac_corr(:,hc_type)
    else
      hc_type=id_type  
      frac_corr(:)=0.     ! no previous channel correlation to consider
    endif
!
! Extract stdv from error table if present
    if (trim(file_error_rc) == 'none') then 
      Rvar_o(:)=0.
      et_pert_fac=0.
    else
      Rvar_o(:)=(1.-frac_corr(:))*et_err_tab(:,2,hc_type)**2
      if (et_l_vc_corr) then 
        do k=1,nk
          Rvar_o(k)=Rvar_o(k)+frac_corr(k)*et_cov_matrix(k,k,hc_type)
        enddo
      endif
    endif 
!
    do k=1,nk
      if (Rvar_o(k) > 0.) then
        Rsqrt_o(k)=sqrt(Rvar_o(k))
      else
        Rsqrt_o(k)=0.
      endif
    enddo
    pert_fac2=et_pert_fac**2
!
! Compute new estimate of stdvd of obs error to add
! Only compute if obs counts >= some min number
    do k=1,nk
      if (Var_o(k) > 1.e-20 .and. C_o(k,k,3) >= xmin_count) then 
        Bvar_o(k)=max(Var_o(k)-Rvar_o(k)*pert_fac2, 0.)
        bfrac_o(k)=Bvar_o(k)/Var_o(k)
      else
        Bvar_o(k)=0.
        bfrac_o(k)=0.
      endif
!
! Only compute if Bvar is minimially reasonable.
! Do not allow new Rvar to be less than 0.25 times previous value
      if (Var_t(k) > 1.e-20 .and. C_t(k,k,3) >= xmin_count &
          .and. Bvar_o(k) > 0.) then  
        Rvar_n(k)=max(Var_t(k)-Bvar_o(k), 0.25*Rvar_o(k))
      else
        Rvar_n(k)=Rvar_o(k)
      endif
!
      if (Rvar_n(k) > 0.) then 
        Rsqrt_n(k)=sqrt(Rvar_n(k))
      else
        Rsqrt_n(k)=0.
      endif
!
! Set quality flag for estimate
      if (C_t(k,k,3) < xmin_count .or. C_o(k,k,3) < xmin_count) then 
        qualflag(k)='QVCT'  ! counts too small
      elseif (Var_t(k) <= 1.e-20 .or. Var_o(k) <= 1.e-20 ) then 
        qualflag(k)='QV00'  ! innovations exceptionally small
      elseif (Var_o(k)-Rvar_o(k)*pert_fac2 < 0.) then
        qualflag(k)='QVBN'  ! estimated B is negative
      elseif (Var_t(k)-Bvar_o(k) < 0.25*Rvar_o(k)) then
        qualflag(k)='QVBS'  ! estimated R is much smaller than previously
      else
        qualflag(k)='    '  ! no obvious problems detected
      endif
!
    enddo
!
! Compute difference in binned covariances if available. 
! Any previously added error variance is added to the difference, so it 
! represents what the difference would be if no added were added.
! Only consider covariances that appear to be sampled well enough
    if (lcorr_both) then
      do k=1,nk
        if (C_d(k,k,3) >= xmin_count) then
          C_d(k,k,1)=max(C_t(k,k,1)-C_o(k,k,1)+pert_fac2*Rvar_o(k),0.)
        else
          C_d(k,k,1)=0.
        endif
        do n=1,nk
          if (n /= k) then
            if (C_d(n,k,3) >= xmin_count) then       
              C_d(n,k,1)=C_t(n,k,1)-C_o(n,k,1)
            else
              C_d(n,k,1)=0.
            endif
          endif
        enddo 
      enddo
!
    endif
!  
! If previously added errors were correlated, then addjust the difference in 
! covariances so that it reflects the total added covariances required. This 
! total is not just what is additionaly required but includes what has been  
! previously added. This total is required if the simulated error is added to
! the original error-free observations rather than to the previous iterate of 
! the observations with added errors.  C_d(:,:,2) is the correlation for 
! each bin. The factor pert_fac2 accounts for its possible adjustment in 
! the previous experiment.
    if (lcorr_prev) then
      do k=1,nk
        do n=1,nk
          if (k /= n) then   ! value for n=k (i.e., variances) comp earlier
            xfac1=sqrt(max((1.-frac_corr(k))*(1.-frac_corr(n)),1.e-20)) 
            xfac2=sqrt(max(frac_corr(k)*frac_corr(n),1.e-20)) 
            xfac=pert_fac2*(xfac1+xfac2)
            C_d(n,k,1)=C_d(n,k,1)+xfac*et_cov_matrix(n,k,hc_type)
          endif
        enddo
      enddo
    endif    
!
! Compute new estimate of correlations of obs errors to add 
! The estimated Cov matrix may require adjutment.
    count_corr_fix(:)=0
    if (lcorr_both) then
      do k=1,nk
        do n=1,nk
          dvar=C_d(k,k,1)*C_d(n,n,1)
          if (dvar > 0.) then
            dvar=1./sqrt(dvar)
          else
            dvar=0.
          endif
          C_d(n,k,2)=dvar*C_d(n,k,1)
!
! Do not allow abs value of correlations greater than corr_max.
! Count the number of non-trivial correls examined and the number 
! of adjustments made.  Also adjust the corresp. covariances. 
          if (n /= k .and. dvar > 0.) then 
            count_corr_fix(1)=count_corr_fix(1)+1 
            if (abs(C_d(n,k,2)) > corr_max) then
              count_corr_fix(2)=count_corr_fix(2)+1 
              if (C_d(n,k,2) > 0.) then   ! adjust corr 
                C_d(n,k,2)=corr_max
              else
                C_d(n,k,2)=-corr_max
              endif
              C_d(n,k,1)=C_d(n,k,2)/dvar  ! adjust cov
            endif
          endif
        enddo
      enddo
      print *,' '
      print ('(a,i8)'), &
            'Number of off-diagonal elements of COV matrix = ', &
            count_corr_fix(1)
      print ('(a,i8,a,f8.5)'), &
            'Number of off-diagonal elements of COV matrix altered = ', &
             count_corr_fix(2),' to render all abs corr <= ',corr_max
    endif    
!
! Optional write corr matrices to a separate file
    if (lwrite_corr) then 
      open (iunit,file='corr_matrix_3',form='unformatted')
      write (iunit) C_t(:,:,2)
      write (iunit) C_o(:,:,2)
      write (iunit) C_d(:,:,2)
      close (iunit)
      print *,' '
      print *,'File corr_matrix_3 written for diagnostic use'
    endif
!
! Write out table of fiting parameters for new corr function
!
    xdum(:)=0.
    write (cname(1),'(l8)') lcorr_prev
    write (cname(2),'(l8)') et_l_vc_corr
!
    open (iunit,file=trim(file_eparams_new))
!
! Write header information
    write (iunit,'(a,4(2x,a16))') & 
                            'dtype,sat_name,instr_name,c_sat_name = ', &
                             dtype,sat_name,instr_name,c_sat_name
    write (iunit,'(a,5i5,f6.0,2x,a)') &
             'nk,itype,nbins,ntimes_t,ntimes_o,xbin,corr_function = ', &
              nk,itype,nk,ntimes_t,ntimes_o,xdum(1),trim(cname(2))
    write (iunit,'(a,3i5,f6.0,2x,a)') &
                            'nk,itype,nbins,xbin,corr_function = ', &
                             nk,itype,nk,xdum(1),trim(cname(2)) 
    write (iunit,'(2a)') 'file_target_stats = ',trim(file_target_stats)
    write (iunit,'(2a)') 'file_osse_stats   = ',trim(file_osse_stats)
    write (iunit,'(2a)') 'file_error_rc     = ',trim(file_error_rc)
!
! Write table of old and new error parameters
!
    write (iunit,'(a)') ' '
    write (iunit,'(a)') 'Table 1: Old and new error parameters: ' 
    write (iunit,'(a,2(2x,a))') 'Previous available, and et_l_vc_corr:  ', &  
                                trim(cname(1)),trim(cname(2)) 
    write (iunit,'(11a)') '   k','           0','   count',           &
                          '   AddSDV_old','   AddSDV_new','   none',  &
                          '   none','   none','   none','  qflag'     
    do k=1,nk
      ncnt=nint(C_d(k,k,3))   
      write (iunit,'(i4,f12.2,i8,2e13.3,2f7.3,2f7.0,2x,a4)')        &
                    k,xdum(1),ncnt,Rsqrt_o(k),Rsqrt_n(k),xdum(1:4), &
                    qualflag(k)
    enddo
!   
! Write table of other read and derived statistics
!
    write (iunit,'(a)') ' '
    write (iunit,'(a)') 'Table 2: Other read and derived statistics: '
    write (iunit,'(10a)') '   k','           0',' count_t',' count_o',   & 
                       '   VarInnov_t','   VarInnov_o','    BkgVarEst',  &
                       '  AddVar_prev','   AddVar_new','  BVar/InnVar'
    do k=1,nk
      write (iunit,'(i4,f12.2,2i8,1p6e13.3)')                       &
                    k,xdum(1),count_t(k),count_o(k),Var_t(k),    &
                    Var_o(k),Bvar_o(k),Rvar_o(k),Rvar_n(k),bfrac_o(k)
    enddo
    close (iunit)
!
    print *,' '
    print *,'Eparams txt file written: ',trim(file_eparams_new)
!
! Write out new covariances to binary file
!
    if (lcorr_both) then
      open (iunit,file=trim(file_cov_new),form='unformatted')
      write (iunit) dtype,sat_name,instr_name,c_sat_name
      write (iunit) nk,itype,lcorr_prev,et_l_vc_corr
      write (iunit) file_target_stats
      write (iunit) file_osse_stats
      write (iunit) file_error_rc
      write (iunit) C_d(:,:,3)
      write (iunit) C_d(:,:,1)
      write (iunit) C_d(:,:,2)
      close (iunit)
      print *,' '
      print *,'Binary file of new covs written: ',trim(file_cov_new)
    endif
!
    call error_table_clean
!
    print *,'Program complete'
!
    end program eparams_rad_hcorr

