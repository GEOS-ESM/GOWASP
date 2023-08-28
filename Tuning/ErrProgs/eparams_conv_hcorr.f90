    program eparams_conv_hcorr
!
! Program to estimate obs error table parameters for radiance obs
! Includes parameters for stdv, fraction of var of correlated error, 
! and corrleation lengths.
!
    use m_kinds, only : rkind1
!
    use m_obs_error_table, only : error_table_setup
    use m_obs_error_table, only : error_table_clean
    use m_obs_error_table, only : error_table_find_stdv_id
    use m_obs_error_table, only : error_table_find_corr_id
    use m_obs_error_table, only : et_pert_fac, et_ks_list
    use m_obs_error_table, only : et_err_tab, et_err_itype, et_n_err1
    use m_obs_error_table, only : et_corr_shapes, et_itypes_corr
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
    logical :: lcorr_both  ! correl for both target and osse avialable 
    logical :: lerrtable   ! indicates that prev table of added stdv available
    logical :: luse_cntfls ! indicates to use count rather than hcorr files
    logical :: luse_t,luse_o 
!
    integer, parameter :: rkind=8   ! precision (8 because of variable in file)
    integer, parameter :: max_corr_funcs=3  ! max(2,num corr shapes considered)
    integer, parameter :: max_nbins=100     ! max num of horiz-distance bins
    integer, parameter :: min_count_cov=100 ! min count to consider cov realiable   
    integer, parameter :: iunit=11          ! unit number for output
!
    integer :: k,n
    integer :: islash_t,islash_o ! character position of last \ in file names
    integer :: ke,ka,kb
    integer :: argc
    integer :: nk         ! number of channels or p-levels
    integer :: ierr       ! return code from subroutines
    integer :: id_type    ! id for desired error column in error table
    integer :: hc_type    ! id for desired obs subtype in file of corr params
    integer :: itype      ! either sat or inst id number to match 
    integer :: nbins      ! min (nbins_t, nbins_o)
    integer :: nbins_t    ! number of bins of horiz. distances for target stats
    integer :: nbins_o    ! number of bins of horiz. distance for OSSE stats
    integer :: nlevs_t    ! number of plevs for target stats
    integer :: nlevs_o    ! number of plevs for osse stats
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
    real(rkind1), parameter :: corr_max=0.95 ! max allowed horiz corr
    real(rkind1) :: xka,xkb,xkc
    real(rkind1) :: pert_fac2  ! square of err pert_fac used in prev osse
    real(rkind1) :: dvar       ! difference in variance 
    real(rkind1) :: xbin_t     ! horizontal width of bins for target sats
    real(rkind1) :: xbin_o     ! horizontal width of bins for osse stats
    real(rkind1) :: xmin_count ! real value of min_count_var
    real(rkind1), allocatable :: plevs_t(:) ! plevs or chan nums for target
    real(rkind1), allocatable :: plevs_o(:) ! plevs or chan nums for osse
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
    real(rkind1), allocatable :: C_d(:,:,:) ! est cov, corr, counts added error 
    real(rkind1), allocatable :: fits(:,:,:)  ! parama for fit of added correl 
!
    character(len=4), allocatable :: qualflag(:)
    character(len=3)   :: ckx
    character(len=16)  :: dtype
    character(len=16)  :: corr_func_n
    character(len=16)  :: corr_func_o
    character(len=16)  :: sat_name
    character(len=16)  :: instr_name      ! PREPBUFR kx number
    character(len=16)  :: field_name
    character(len=240) :: file_target_stats
    character(len=240) :: file_osse_stats
    character(len=240) :: file_error_rc
    character(len=240) :: file_eparams_new
!
    print *,' '
    print *,'BEGIN PROGRAM estimate_eparams_conv_hcorr'
!
! Read arguments
    argc = iargc()
    if (argc .ne. 7) then
      print *,'useage must be: estimate.x '                             
      stop
    endif
    call GetArg( 1_4, instr_name)
    call GetArg( 2_4, field_name)
    call GetArg( 3_4, file_error_rc)
    call GetArg( 4_4, file_target_stats)
    call GetArg( 5_4, file_osse_stats)
    call GetArg( 6_4, file_eparams_new)
    call GetArg( 7_4, corr_func_n)
!
    print *,'kx number=',trim(instr_name),' field=',trim(field_name)
!
    dtype='PREPBUFR'
    ckx=instr_name(1:3)
    read (ckx,'(i3)') itype  ! kx value for obs in table  
    if (trim(field_name) == 'p' ) then
      id_type=5                          ! column in PREPBUFR error table
    else if (trim(field_name) == 'T' ) then
      id_type=2
    else if (trim(field_name) == 'u' ) then
      id_type=4
    else if (trim(field_name) == 'v' ) then
      id_type=4
    else
      print *,'Option ',field_name,' not allowed'
      stop
    endif
!
! Determine last character in file path part of file name 
    islash_t=0
    islash_o=0
    do n=1,len(file_target_stats)    
      if (file_target_stats(n:n) == '/' ) islash_t=n
      if (file_osse_stats(n:n) == '/' ) islash_o=n
    enddo
!
! Determine if pair of stat file types is acceptable
    luse_t=file_target_stats(islash_t+1:islash_t+5) == 'count' 
    luse_o=file_osse_stats(islash_o+1:islash_o+5) == 'count'
    if (luse_t .and. luse_o) then 
      luse_cntfls=.true.
    elseif  ((.not. luse_t) .and. (.not. luse_o)) then
      luse_cntfls=.false.
    else
      print *,' '
      print *,'Mismatched input stat file types : Both must be' 
      print *,'hcorr*bin or count*txt files due to different nk values'
      print *,'use_t, use_o =',luse_t,luse_o
      stop
    endif
!
! Number of p-levels of data determined here by stat file type
    if (luse_cntfls) then 
      nk=12   
    else 
      nk=10
    endif
    xmin_count=real(min_count_cov)
!
! Allocate arrays whose dimennsions depend only on n_channels
    allocate (qualflag(nk))
    allocate (count_t(nk))
    allocate (count_o(nk))
    allocate (count_d(nk))
    allocate (plevs_t(nk))
    allocate (plevs_o(nk))
    allocate (Var_t(nk))
    allocate (Var_o(nk))
    allocate (Rsqrt_o(nk))
    allocate (Rsqrt_n(nk))
    allocate (Rvar_o(nk))
    allocate (Rvar_n(nk))
    allocate (Bvar_o(nk))
    allocate (bfrac_o(nk))
    allocate (fits(nk,3,max_corr_funcs))  ! 1-3: fvarcorr, corrlength, scorefit
    allocate (C_t(max_nbins,nk,3))        ! 1-3: cov, corr, bin count measure
    allocate (C_o(max_nbins,nk,3))        ! 1-3: cov, corr, bin count measure 
    allocate (C_d(max_nbins,nk,4))        ! 1-3: as for C_o; 4: corr prev eparms
!
! Read error.rc file and the files possibly specified within it, including file 
! of prev error stdv table and file of added error correlation parameters.
    if (trim(file_error_rc) == 'none') then 
      lcorr_prev=.false.
      lerrtable=.false.      
    else
      call error_table_setup (.true.,dtype,file_error_rc,lerrtable, &
                              lcorr_prev,idum2(1),ierr)
      if (lcorr_prev) then
        lcorr_prev=.false.
        do n=1,et_itypes_corr
          do k=1,10
            if (itype == et_ks_list(k,n)) lcorr_prev=.true.
          enddo
        enddo
      endif 
    endif 
!
! Read target statistics from 1 of 2 types of stat files 
!
    if (luse_t) then
!
! Get counts and standard deviations of innovations from count_tables
! (No horizontal covariances read)
      call read_counts (nk,ntimes_t,file_target_stats,count_t,plevs_t,Var_t)
      Var_t(:)=Var_t(:)**2  ! replace read stdv by var
      nbins_t=0
      xbin_t=0.
      do k=1,nk  
        C_t(1,k,3)=count_t(k)*ntimes_t
      enddo
!
    else   
!
! Read binned horizontal counts and covariances of innovations 
! from hcorr_tables
      call read_hcovs (max_nbins,nk,ntimes_t,file_target_stats, &
                       nlevs_t,nbins_t,ierr,xbin_t,             &
                       count_t,plevs_t,C_t)
      if (ierr == 0) then
        Var_t(:)=C_t(1,:,1)    
      else
        print *,'ERROR in 1st call to read_hcovs: ierr=',ierr
        stop
      endif 
    endif
!
! Read osse statistics from 1 of 2 types of stat files 
!
    if (luse_o) then
!
! Get counts and standard deviations of innovations from count_tables
! (No horizontal covariances read)
      call read_counts (nk,ntimes_o,file_osse_stats,count_o,plevs_o,Var_o)
      Var_o(:)=Var_o(:)**2  ! replace read stdv by var
      nbins_o=0
      xbin_o=0.
      do k=1,nk  
        C_o(1,k,3)=count_o(k)*ntimes_o
      enddo
!
    else   
!
! Read binned horizontal counts and covariances of innovations 
! from hcorr_tables
      call read_hcovs (max_nbins,nk,ntimes_o,file_osse_stats, &
                       nlevs_o,nbins_o,ierr,xbin_o,           &
                       count_o,plevs_o,C_o)
      if (ierr == 0) then
        Var_o(:)=C_o(1,:,1)    
      else
        print *,'ERROR in 2nd call to read_hcovs: ierr=',ierr
        stop
      endif 
    endif
!
! Determine if it is possible to compare target and osse covariances
    nbins=min(nbins_t,nbins_o)
    lcorr_both=nbins > 0 .and. xbin_t == xbin_o
    If (.not. lcorr_both) then 
      print *,' '
      print *,'Correlations cannot be compared because either ', &
              'covariance have not'
      print *,'been read for both target and osse or their bin widths differ:' 
      print ('(a,2i4,2f8.2)'),'nbins_t,nbins_o,xbin_t,xbin_o = ', &
                               nbins_t,nbins_o,xbin_t,xbin_o
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
! Extract stdv from error table if present
    if (trim(file_error_rc) == 'none') then 
      Rsqrt_o(:)=0.  
      et_pert_fac=0.
    else
      do k=1,nk
        kb=1
        do ke=1,et_n_err1-1
          if (plevs_o(k) < 0.01*et_err_tab(ke,1,itype-99)) then
            kb=ke
          endif
        enddo
        ka=kb+1
!
! Either set Rsqrt to 2nd value from top in table or linearly 
! interpolate in p using table values
        xka=plevs_o(k)-0.01*et_err_tab(ka,1,itype-99)
        xkb=0.01*et_err_tab(kb,1,itype-99)-plevs_o(k)
        xkc=1.0/(xka+xkb)
        Rsqrt_o(k)=xkc*(xkb*et_err_tab(ka,id_type,itype-99)+ &
                        xka*et_err_tab(kb,id_type,itype-99))
      enddo
    endif 
    Rvar_o(:)=Rsqrt_o(:)**2
    pert_fac2=et_pert_fac**2
!
! Compute new estimate of stdvd of obs error to add
! Only compute if obs counts >= some min number
!
    do k=1,nk
!
! Compute estimate of Bkg err var in prev OSSE as a residual
! (Note that this estimate includes implicit rep. error.)
      if (Var_o(k) > 1.e-20 .and. C_o(1,k,3) >= xmin_count) then 
        Bvar_o(k)=max(Var_o(k)-Rvar_o(k)*pert_fac2, 0.)
        bfrac_o(k)=Bvar_o(k)/Var_o(k)
      else
        Bvar_o(k)=0.
        bfrac_o(k)=0.
      endif
!
! Only compute if Bvar is minimially reasonable.
! Do not allow new Rvar to be less than 0.25 times previous value
      if (Var_t(k) > 1.e-20 .and. C_t(1,k,3) >= xmin_count &
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
      if (C_t(1,k,3) < xmin_count .or. C_o(1,k,3) < xmin_count) then 
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
    if (lcorr_both) then
      do k=1,nk
        C_d(1,k,1)=max(C_t(1,k,1)-C_o(1,k,1)+pert_fac2*Rvar_o(k),0.)
        do n=2,nbins
          C_d(n,k,1)=C_t(n,k,1)-C_o(n,k,1)
        enddo
      enddo
    endif
!  
! If previously added errors were correlated, then adjust the difference in 
! covariances so that it reflects what the new covariances would be without 
! any added error.
    corr_func_o='WHITE'      ! value if lcorr_prev=.false.
    fits(:,:,:)=0.           ! default value
    if (lcorr_prev) then
!
! Retrieve previous correlation parameters for this data type
! If hc_type=0 then obs type not found in list of correlated obs, and
! therefore it is assumed WHITE
      call error_table_find_corr_id (itype,hc_type)  
      if (hc_type == 0) then 
        print *,'NO CORRELATION PARAMETERS PREVIOUSLY SPECIFIED'
      else  
        corr_func_o=et_corr_shapes(hc_type)  ! shape of corr prev used
        fits(:,1,1)=et_frac_corr(:,hc_type)
        fits(:,2,1)=et_hcorr_lengths(:,hc_type)
        print *,'Previous correlation shape = ',trim(corr_func_o)
!
! Construct the prev. added error cov from the previous correl function 
! parameters.  C_d(:,:,4) is the correlation for each bin. The factor 
! Rvar changes the correl to covariances and the factor pert_fac2 accounts
! for its possible adjustment in the previous experiment.
        call construct_corr (max_nbins,nbins,nk,corr_func_o, &
                             xbin_o,C_d(:,:,4),fits(:,:,1),ierr)
        do k=1,nk
          do n=2,nbins  ! value for n=1 (i.e., variances) comp earlier
            C_d(n,k,1)=C_d(n,k,1)+pert_fac2*Rvar_o(k)*C_d(n,k,4)
          enddo
        enddo
!  
      endif   ! hc_type
    endif     ! check on lcorr_prev 
!
    if (lcorr_both) then
!
!  Compute new estimate of correlations of obs errors to add 
      do k=1,nk
        if (C_d(1,k,1) > 0.) then 
          dvar=C_d(1,k,1)
        else
          dvar=0.
        endif
        C_d(1,k,2)=1.
        if (dvar == 0.) then
          C_d(2:nbins,k,1:2)=0. 
        else
!
! Do not allow abs value of correlations greater than corr_max.
! Count the number of non-trivial correls examined and the number 
! of adjustments made.  Also adjust the corresp. covariances. 
          do n=2,nbins
            count_corr_fix(1)=count_corr_fix(1)+1 
            C_d(n,k,2)=C_d(n,k,1)/dvar
            if (abs(C_d(n,k,2)) > corr_max) then 
              count_corr_fix(2)=count_corr_fix(2)+1 
              if (C_d(n,k,2) > 0.) then 
                C_d(n,k,2)=corr_max
              else
                C_d(n,k,2)=-corr_max
              endif
              C_d(n,k,1)=C_d(n,k,2)*dvar
            endif
          enddo
!
        endif   ! check on dvar
      enddo     ! loop over k
!
      print *,' '
      print ('(a,i8)'), &
            'Number of non-trivial correl bins = ',count_corr_fix(1)
      print ('(a,i8,a,f8.5)'), &
            'Number of correl bins altered = ',count_corr_fix(2), &
            ' to render all abs corr <= ',corr_max
!
! Fit correlations to prescribed funtion shapes
      call fit_corr (max_nbins,nbins,nk,min_count_cov,corr_func_n, &
                     xbin_t,C_d(:,:,1:3),fits(:,:,2),ierr)
    endif    
!
! Write out table of fiting parameters for new corr function
!
    open (iunit,file=trim(file_eparams_new))
!
! Write header information
    write (iunit,'(a,4(2x,a16))') & 
                            'dtype,field_name,instr_name,field_nm = ', &
                             dtype,field_name,instr_name,field_name
    write (iunit,'(a,5i5,f6.0,2x,a)') &
             'nk,itype,nbins,ntimes_t,ntimes_o,xbin,corr_function = ', &
              nk,itype,nbins,ntimes_t,ntimes_o,xbin_t,trim(corr_func_n)
    write (iunit,'(2a)') 'file_target_stats = ',trim(file_target_stats)
    write (iunit,'(2a)') 'file_osse_stats   = ',trim(file_osse_stats)
    write (iunit,'(2a)') 'file_error_rc     = ',trim(file_error_rc)
!
! Write table of old and new error parameters
!
    write (iunit,'(a)') ' '
    write (iunit,'(a)') 'Table 1: Old and new error parameters: ' 
    write (iunit,'(a,2(2x,a))') 'Previous and new corr function shapes: ', &
                                trim(corr_func_o),trim(corr_func_n) 
    write (iunit,'(11a)') '   k','     p_or_ch','   count',           &
                          '   AddSDV_old','   AddSDV_new','  f_old',  &
                          '  f_new','  L_old','  L_new','  qflag'         
    do k=1,nk
      ncnt=min(count_t(k),count_o(k))   
      write (iunit,'(i4,f12.2,i8,2e13.3,2f7.3,2f7.0,2x,a4)')         &
                    k,plevs_t(k),ncnt,Rsqrt_o(k),Rsqrt_n(k),         &
                    fits(k,1,1),fits(k,1,2),fits(k,2,1),fits(k,2,2), &
                    qualflag(k)
    enddo
!   
! Write table of other read and derived statistics
!
    write (iunit,'(a)') ' '
    write (iunit,'(a)') 'Table 2: Other read and derived statistics: '
    write (iunit,'(10a)') '   k','     p_or_ch',' count_t',' count_o',   & 
                       '   VarInnov_t','   VarInnov_o','    BkgVarEst',  &
                       '  AddVar_prev','   AddVar_new','  BVar/InnVar'
    do k=1,nk
      write (iunit,'(i4,f12.2,2i8,1p6e13.3)')                       &
                    k,plevs_t(k),count_t(k),count_o(k),Var_t(k),    &
                    Var_o(k),Bvar_o(k),Rvar_o(k),Rvar_n(k),bfrac_o(k)
    enddo
!
! Write table of fits of Corr to be added, considering all function shapes
    if (lcorr_both) then 
      call fit_corr (max_nbins,nbins,nk,min_count_cov,'GAUSS', &
                     xbin_t,C_d,fits(:,:,1),ierr)
      call fit_corr (max_nbins,nbins,nk,min_count_cov,'EXP',   &
                     xbin_t,C_d,fits(:,:,2),ierr)
      call fit_corr (max_nbins,nbins,nk,min_count_cov,'TOAR',  &
                     xbin_t,C_d,fits(:,:,3),ierr)
!        
      write (iunit,'(a)') ' '
      write (iunit,'(a)') 'Table 3: Fits to est obs error cov:'
      write (iunit,'(24x,3(11x,a,10x))') 'GAUSS',' EXP ',' TOAR' 
      write (iunit,'(12a)') '   k','     p_or_ch','   count', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore'
      do k=1,nk
        ncnt=min(count_t(k),count_o(k))   
        write (iunit,'(i4,f12.2,i8,3(2x,f7.3,f7.0,e10.2))')  &
                      k,plevs_t(k),ncnt,fits(k,:,:)
      enddo
    endif
!
! Write table of fits of Corr for target innov, considering all function shapes
    if (nbins_t) then 
      call fit_corr (max_nbins,nbins_t,nk,min_count_cov,'GAUSS', &
                     xbin_t,C_t,fits(:,:,1),ierr)
      call fit_corr (max_nbins,nbins_t,nk,min_count_cov,'EXP',   &
                     xbin_t,C_t,fits(:,:,2),ierr)
      call fit_corr (max_nbins,nbins_t,nk,min_count_cov,'TOAR',  &
                     xbin_t,C_t,fits(:,:,3),ierr)
      write (iunit,'(a)') ' '
      write (iunit,'(a)') 'Table 4: Fits to target covs:'
      write (iunit,'(24x,3(11x,a,10x))') 'GAUSS',' EXP ',' TOAR' 
      write (iunit,'(12a)') '   k','     p_or_ch','   count', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore'
      do k=1,nk
        write (iunit,'(i4,f12.2,i8,3(2x,f7.3,f7.0,e10.2))')  &
                      k,plevs_t(k),count_t(k),fits(k,:,:)
      enddo
    endif
!
! Write table of fits of Corr for OSSE innov, considering all function shapes
    if (nbins_o) then 
      call fit_corr (max_nbins,nbins_o,nk,min_count_cov,'GAUSS', &
                     xbin_o,C_o,fits(:,:,1),ierr)
      call fit_corr (max_nbins,nbins_o,nk,min_count_cov,'EXP',   &
                     xbin_o,C_o,fits(:,:,2),ierr)
      call fit_corr (max_nbins,nbins_o,nk,min_count_cov,'TOAR',  &
                     xbin_o,C_o,fits(:,:,3),ierr)
      write (iunit,'(a)') ' '
      write (iunit,'(a)') 'Table 5: Fits to OSSE covs:'
      write (iunit,'(24x,3(11x,a,10x))') 'GAUSS',' EXP ',' TOAR' 
      write (iunit,'(12a)') '   k','     p_or_ch','   count', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore', &
                          '    fcorr','  Lcorr','  fitscore'
      do k=1,nk
        write (iunit,'(i4,f12.2,i8,3(2x,f7.3,f7.0,e10.2))')  &
                      k,plevs_o(k),count_o(k),fits(k,:,:)
      enddo
    endif
!
    close (iunit)
!
    print *,' '
    print *,'Eparams txt file written: ',trim(file_eparams_new)
    call error_table_clean
    print *,'Program complete'
!
    end program eparams_conv_hcorr

