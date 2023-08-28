      program new_corr_err_params
!
! Write file containing parameters and statistics required for specifying 
! the random fields.  Separate files are required for each data type.
! (All onventional observations are handeled by a single file with up to
! 5 different sets of parameter specifications.)
!
      use m_sat_info_table, only : sat_info_table_read
!     
      implicit none
      logical, parameter :: lprint_mat=.false. ! true if large matrics printed
      integer, parameter :: r4=4
      integer, parameter :: r8=8
      integer, parameter :: iunit=10
      integer, parameter :: nchars=140
!
      logical :: l_vc_corr
      integer :: ier
      integer :: id,k
      integer :: itypes_corr
      integer :: file_format_in
      integer :: file_format_out
      integer :: nlevels
      integer :: nlevels_read
      integer, allocatable :: ks_list(:,:)
      integer :: argc
      integer :: nmax
!
      real(r4) :: r_version_num         ! version number of file to be created
      real(r4) :: pmax                  ! if conv obs, max p-level 
      real(r4) :: pmin                  ! if conv obs, min p-level
      real(r4) :: h_min                 ! min value of horiz corr allowed
      real(r4), allocatable :: fix_params(:,:)
      real(r4), allocatable :: v_lengths(:)
      real(r4), allocatable :: e_vects(:,:,:) !  eofs of vert(chan) cov m
      real(r4), allocatable :: h_lengths(:,:) ! horiz corr lengths (km)
      real(r4), allocatable :: frac_corr(:,:) ! fraction of tot var correl
      real(r4), allocatable :: e_values(:,:)  ! eigenvalues of cov matrix
      real(r4), allocatable :: cov(:,:,:)
      real(r4), allocatable :: sigma_read(:,:)
      real(r4), allocatable :: plevs_read(:)
      real(r4), allocatable :: plevels(:)
      real(r4), allocatable :: sigma(:,:)
      real(r8), allocatable :: e_vects_r8(:,:) 
      real(r8), allocatable :: e_values_r8(:) 
      real(r8), allocatable :: cov_r8(:,:) 
!
      character(len=20), allocatable :: corr_shapes(:)
      character(len=13) :: c_input(20)
      character(len=16) :: d_type
      character(len=6)  :: c_version_num   ! version number (compatible with f6.0)
      character(len=nchars) :: c_info      ! description of file to create
      character(len=nchars) :: orig_path   ! path of original file
      character(len=nchars) :: sigma_file  ! name of file of sigma values to read
      character(len=nchars) :: param_file_input     ! name of file of input parameters
      character(len=nchars) :: file_sat_info        ! sat info file for rad obs
      character(len=nchars) :: hcorr_file_output    ! name of file to be written
      character(len=nchars) :: hcorr_file_params(2) ! name of horiz corr params
      character(len=nchars) :: file_names(2)        ! file path + names
      character(len=nchars) :: adjust_name
      character(len=nchars), allocatable :: vcov_files(:)  ! vert cov files to be read
!
      data file_format_out /2/ 
      data h_min /10./ 
!
      argc = iargc()
      if (argc .ne. 3) then
        print *,' usage must be: '
        print *,' prog.x param_file_input hcorr_file_output sat_info_file'
        stop
      endif
      call GetArg( 1_4, param_file_input)
      call GetArg( 2_4, hcorr_file_output)
      call GetArg( 3_4, file_sat_info)
!     
      open (unit=32,file=trim(param_file_input),form='formatted')
      read (32,*) c_input(1), d_type
      read (32,*) c_input(2), file_format_in
      read (32,*) c_input(3), itypes_corr  ! num of param sets (<6)
      read (32,*) c_input(4), nlevels      ! num of chans or levs 
      read (32,*) c_input(5), pmin, pmax    ! 
      read (32,*) c_input(6), nmax         ! spect trunc ran fields
      read (32,*) c_input(7), c_version_num ! version # for output
      read (32,*) c_input(8), c_info       ! info to write to file
      read (32,'(2a)') c_input(9), sigma_file   ! file of sigma values
      read (32,'(2a)') c_input(10),orig_path    ! path of est files used
!
      read (c_version_num,'(f6.0)') r_version_num
!
      print ('(a,1x,a)'),  c_input(1),  trim(d_type)
      print ('(a,i3)'),    c_input(2),  file_format_in
      print ('(a,i3)'),    c_input(3),  itypes_corr    
      print ('(a,i4)'),    c_input(4),  nlevels        
      print ('(a,2f9.2)'), c_input(5),  pmin,pmax
      print ('(a,i5)'),    c_input(6),  nmax
      print ('(a,f6.2)'),  c_input(7),  r_version_num  
      print ('(a,1x,a)'),  c_input(8),  trim(c_info)    
      print ('(a,1x,a)'),  c_input(9),  trim(sigma_file) 
      print ('(a,1x,a)'),  c_input(10), trim(orig_path) 
!
! Set number of levels to read from file of error standard deviations
! Set pmax and pmin (units Pascals). For radiances, simply set pmax<pmin 
! so that these values are not used (since there are only channel, not 
! vertical correlations).  For CONV obs, set pmin<pmax so that vert range 
! of obs is covered. Set v_lengths to the vertical correlation lengths in 
! units of meters, one value for each data type subset; for radiance set 
! to 0 since not relevant.  
      if ((trim(d_type) == 'WIND') .or. (trim(d_type) == 'MASS')) then
        nlevels_read=33
        if (pmin > pmax) then
          print *,'ERROR: pmin > pmax for conventional data not allowed'
          stop
        endif
        allocate (plevels(nlevels))
        call define_plevels (nlevels,pmax,pmin,plevels)
      else
        call sat_info_table_read (file_sat_info,.true.,ier)
        nlevels_read=nlevels
        pmin=10.
        pmax= 1. 
        print *,' pmax, pmin changed so that no interp is done'
      endif
!
      allocate (ks_list(10,itypes_corr))
      allocate (fix_params(4,itypes_corr))
      allocate (v_lengths(itypes_corr))
      allocate (corr_shapes(itypes_corr))
      allocate (e_vects_r8(nlevels,nlevels))
      allocate (cov_r8(nlevels,nlevels))
      allocate (e_values_r8(nlevels))
      allocate (e_vects(nlevels,nlevels,itypes_corr))
      allocate (cov(nlevels,nlevels,itypes_corr))
      allocate (e_values(nlevels,itypes_corr))
      allocate (h_lengths(nlevels,itypes_corr))
      allocate (frac_corr(nlevels,itypes_corr))
      allocate (sigma(nlevels,itypes_corr))
      allocate (sigma_read(nlevels_read,itypes_corr))
      allocate (plevs_read(nlevels_read))
      allocate (vcov_files(itypes_corr))
      e_values=0.
!      
! Read some params and file names containing hcorr parameters or chan covs
      fix_params(:,:)=0.
      ks_list(:,:)=0
      do id=1,itypes_corr
        read (32,*) c_input(11)    ! skip
        read (32,*) c_input(12), ks_list(:,id)
        read (32,*) c_input(13), v_lengths(id)  ! vert corr len (m)
        read (32,*) c_input(14), fix_params(:,id)
        read (32,*) c_input(15), vcov_files(id)
        read (32,*) c_input(16), hcorr_file_params(1)
        read (32,*) c_input(17), hcorr_file_params(2)
!
        print *,' '
        print ('(a,i2)'),'Index for set of read correlation parameters=',id
        print ('(a,10i5)'),  c_input(12),ks_list(:,id)
        print ('(a,f6.0)'),  c_input(13),v_lengths(id)
        print ('(a,4f8.2)'), c_input(14),fix_params(:,id)
        print ('(a,a)'),     c_input(15),trim(vcov_files(id))
        print ('(a,a)'),     c_input(16),trim(hcorr_file_params(1))
        print ('(a,a)'),     c_input(17),trim(hcorr_file_params(2))
!
! Read file(s) of horiz correl parameters if requested
! (can be 2 files if u and v statistics are to be averaged)
! Also adjust parameters as requested
        if (trim(hcorr_file_params(1)) /= 'none') then
          call read_hcorr_params (d_type,hcorr_file_params,orig_path, &
                                  nlevels,pmin,pmax,plevels,          &
                                  corr_shapes(id),frac_corr(:,id),    &
                                  h_lengths(:,id))
        else  ! will be set by adjust_hcorr_params
          frac_corr(:,id)=0.   
          h_lengths(:,id)=0.
        endif    
        call adjust_hcorr_params (nlevels,fix_params(:,id), &
                                  frac_corr(:,id),h_lengths(:,id))
      enddo
      close (32)
!
! Set l_vc_corr true if either vertical or channel correlations considered 
      l_vc_corr=.false. 
      do id=1,itypes_corr
        if (trim(vcov_files(id)) /= 'none') l_vc_corr=.true.     
        if (v_lengths(id) > 0.) l_vc_corr=.true.     
      enddo
!
! Create vertical or channel covariance matrix if requested
! Adjust this matrix so that it pertains only to the corrrelated part 
! of the error. No adjustment if conv obs since cov matrix computed assuming
! it describes only uncorrelated part of error. In contrast, the cov matrix 
! read for radiance obs types includes an estimate of combined correlated
! and cuncorrelated parts. 
      if (l_vc_corr) then
        call read_sigmas (nlevels_read,itypes_corr,d_type,sigma_file,&
                          file_sat_info,plevs_read,sigma_read,ks_list) 
        do id=1,itypes_corr
          if ((trim(d_type) == 'WIND') .or. (trim(d_type) == 'MASS')) then
            call fill_cov_matrix_conv (nlevels_read,nlevels,pmax,pmin,  &
                           v_lengths(id),sigma_read(:,id),plevs_read,   &
                           plevels,sigma(:,id),cov(:,:,id))
          else 
            if (trim(vcov_files(id)) /= 'none') then
              call fill_cov_matrix_rad (nlevels,vcov_files(id),orig_path, &
                                        cov(:,:,id))
              call adjust_matrix (id,itypes_corr,d_type,nlevels,    &
                                  frac_corr(:,id),sigma_read(:,id), &
                                  cov(:,:,id))
            else
              cov(:,:,id)=0.
            endif
          endif
        enddo
!
        call print_p_and_v (nlevels_read,itypes_corr,plevs_read, &
                            sigma_read,'plevs_read','sigma_read_all_types')
        if (allocated(plevels)) then   
          call print_p_and_v (nlevels,itypes_corr,plevels, &
                              sigma,'plevels','sigma_all_types')
        endif
        if (lprint_mat) then 
          do id=1,itypes_corr
            call print_m (nlevels,nlevels,cov(:,:,id),'cov matrix',id)
          enddo
        endif
!
! If required, compute eigen values and normalized eigen vectors of cov matrix
! Also, do not allow relatively small or negative eigenvalues
        do id=1,itypes_corr
          cov_r8(:,:)=cov(:,:,id)
          call pert_eigen (nlevels,cov_r8,e_vects_r8,e_values_r8)
          call pert_normalize_evects (nlevels,e_vects_r8)
          call test_evalues (nlevels,cov(:,:,id),e_vects_r8,e_values_r8,id)
          call pert_fix_evalues (nlevels,e_values_r8)
          e_values(:,id)=e_values_r8(:)
          e_vects(:,:,id)=e_vects_r8(:,:)
        enddo
!       
        call print_m (nlevels,itypes_corr,e_values,'e_values',999)
        if (lprint_mat) then 
          do id=1,itypes_corr
            call print_m (nlevels,nlevels,e_vects(1,1,id),'EOF matrix',id)
          enddo
        endif  
!
! Adjust h_lengths, since on input these are assumed to be
! indexed to levels or channels but in this vc correlated case must be 
! indexed to the eigenvectors (EOFs or PCs). 
        do id=1,itypes_corr
          call adjust_h_vc (nlevels,h_min,e_vects(1,1,id),h_lengths(:,id))
        enddo
!
      else
        print *,'No vertical or channel correlations considered'
      endif ! check on l_vc_corr
!
! Write ouput file
      open (iunit,file=trim(hcorr_file_output),form='unformatted')
      write (iunit) itypes_corr, d_type, file_format_out
      write (iunit) r_version_num,c_info,param_file_input
      write (iunit) nlevels,pmax,pmin,l_vc_corr
      write (iunit) v_lengths
      write (iunit) ks_list
      write (iunit) nmax
      write (iunit) corr_shapes
!
      print ('(2a,i2)'),' Number of types of obs to have horizontally ', &
              'correlated errors for this class =',itypes_corr
      print ('(a,f6.2)'),' Version number =',r_version_num
      print *,' ',trim(c_info)
      print ('(a,i4)'),' Number of levels or channels =', nlevels
      print ('(a,l8)'),' Vertical or channel correlations also considered =', &
              l_vc_corr
      print ('(a,2f10.2,i5)'),' pmax, pmin, nmax =',pmax,pmin,nmax
      if (pmax > pmin) then 
        do id=1,itypes_corr
          print ('(a,i2,a,10i6)'),'  ks_list(:,',id,')=',ks_list(:,id)
        enddo 
      endif
!
      do id=1,itypes_corr
        print ('(a,i2,a,4f10.3)'),'  fix_params(:,',id,')=',fix_params(:,id)
      enddo 
!
      write (iunit) frac_corr
      write (iunit) h_lengths
!
      if (l_vc_corr) then 
        write (iunit) cov
        write (iunit) e_values
        write (iunit) e_vects
      endif
      close (iunit)
!
      print *,' Table of h_lengths, frac_corr, and e_values follows'
      do id=1,itypes_corr
        print ('(a,i2,2a)'),' itype=',id,'   corr_shape=',trim(corr_shapes(id))
        do k=1, nlevels
          print ('(a,i4,f10.1,f10.4,1p1e14.3)'),'k=',k, &
               h_lengths(k,id),frac_corr(k,id),e_values(k,id)
        enddo
      enddo        
!
      print *,' All information on file=',trim(hcorr_file_output),' written'
!
      end program new_corr_err_params
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine fill_cov_matrix_conv (nlevels_read,nlevels,pmax,pmin,  &
                                       v_lengths,sigma_read,plevs_read, & 
                                       plevels,sigma,cov)
!
      implicit none
      integer, parameter :: r4=4
!
      integer :: nlevels
      integer :: nlevels_read
      real(r4) :: pmax, pmin
      real(r4) :: v_lengths
      real(r4) :: sigma_read(nlevels_read)
      real(r4) :: plevs_read(nlevels_read)
      real(r4) :: plevels(nlevels)
      real(r4) :: sigma(nlevels)
      real(r4) :: cov(nlevels,nlevels)
!
      integer  :: j, k, k1, k1s
      real(r4), parameter :: T_ref=270.
      real(r4), parameter :: R_gas=287.
      real(r4), parameter :: Grav=9.8
      real(r4) :: dlogp, pfac
      real(r4) :: w1, w2
      real(r4) :: x
      real(r4) :: xfac
!
! Define pressure levels for random fields

      pfac=(pmin/pmax)**(1./real(nlevels-1))
      plevels(1)=pmax
      do k=2,nlevels-1
        plevels(k)=plevels(k-1)*pfac
      enddo 
      plevels(nlevels)=pmin
!
! Define sigma on levels of random fields by interpolation from table values.
! The interpolation is linear in log(p)
      k1s=1
      do k=1,nlevels
!
! Find the table level just above the desired level
        do k1=k1s,nlevels_read
          if (plevs_read(k1) < plevels(k)) then
            exit
          endif
        enddo
        k1s=k1
!       
        dlogp=log(plevs_read(k1-1)/plevs_read(k1))
        w1=log(plevs_read(k1-1)/plevels(k))/dlogp
        w2=1.-w1
        sigma(k)=w1*sigma_read(k1)+w2*sigma_read(k1-1)
      enddo  
!
! Fill symmetric covariance matrix
      xfac=-0.5*(R_gas*T_ref/(Grav*v_lengths))**2
      do k=1,nlevels
        do j=k,nlevels
          x=log(plevels(k)/plevels(j))
          cov(k,j)=sigma(k)*sigma(j)*exp(xfac*x*x)
          cov(j,k)=cov(k,j)
        enddo
      enddo
!
      end subroutine fill_cov_matrix_conv 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine fill_cov_matrix_rad (nlevels,vcov_file,file_path,cov)
!
      implicit none
      integer, parameter :: r4=4
!
      integer :: nlevels
      real(r4) :: cov(nlevels,nlevels)
      character(len=*) :: vcov_file
      character(len=*) :: file_path
!
      integer, parameter  :: iunit=23
      integer  :: nk_read
      character(len=16)  :: cread(4)
      character(len=240) :: file_name
      character(len=240) :: file_target_stats 
      character(len=240) :: file_osse_stats 
      character(len=240) :: file_error_rc 
!
      file_name=trim(file_path)//'/'//trim(vcov_file)
      open (iunit,file=trim(file_name),form='unformatted')      
      print *,'File opened to read channel covariance matrix:'
      print *,'file_name= ',trim(file_name)
! 
      read (iunit) cread(:)
      read (iunit) nk_read
      read (iunit) file_target_stats
      read (iunit) file_osse_stats
      read (iunit) file_error_rc
      read (iunit) ! skip xbins counts
!
      print *,'cread= ',cread(:)
      print *,'file_target_stats= ',trim(file_target_stats)
      print *,'file_osse_stats= ',trim(file_osse_stats)
      print *,'file_error_rcs= ',trim(file_error_rc)
!
      if (nk_read /= nlevels) then
        print *,'ERROR: nk_read /= nlevels in fill_cov_matrix_rad:'
        print *,'nk_read, nlevels = ',nk_read,nlevels
        stop
      endif 
!
      read (iunit) cov
      close (iunit)
!
      end subroutine fill_cov_matrix_rad 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine read_sigmas (nlevels,itypes,d_type,sigma_file, &
                              file_sat_info,plevs,sigma,ks_list) 
!
      use m_obs_error_table, only : error_table_read_stdv
      use m_obs_error_table, only : error_table_clean
      use m_obs_error_table, only : error_table_find_stdv_id
      use m_obs_error_table, only : et_file_err_stdv
      use m_obs_error_table, only : et_err_tab
!
      implicit none
      integer, parameter :: r4=4
!
      integer :: nlevels
      integer :: itypes
      integer :: ks_list(10,10)
      real(r4) :: plevs(nlevels)
      real(r4) :: sigma(nlevels,itypes)
      character(len=*) :: d_type
      character(len=*) :: sigma_file
      character(len=*) :: file_sat_info
!
      integer :: k,id
      integer :: ier1
      integer :: id_field, id_type
      character(len=16) :: dtype_table        
!
      print *,' '
      print *,' Subroutine read_sigmas entered'
!
      if (trim(d_type) == 'WIND' .or. trim(d_type) == 'SATWIND' .or. & 
          trim(d_type) == 'MASS' .or. trim(d_type) == 'PREPBUFR') then
        dtype_table='PREPBUFR'
      else
        dtype_table=d_type
      endif
!
!
      et_file_err_stdv=sigma_file
      call error_table_read_stdv (.true.,dtype_table,ier1)
!
!
      if (dtype_table == 'PREPBUFR') then
        do id=1,itypes
          if (ks_list(1,id) < 200) then 
            id_field=2    ! values for T
          else
            id_field=4    ! values for wind
          endif
!
          id_type=ks_list(1,id)-99 
          do k=1,nlevels
            plevs(k)=0.01*et_err_tab(k,1,id_type)
            sigma(k,id)=et_err_tab(k,id_field,id_type)
          enddo
        enddo
!
      else  ! a radiance type
!
        do id=1,itypes
          call error_table_find_stdv_id (ks_list(1,id),id_type)
          do k=1,nlevels
            sigma(k,id)=et_err_tab(k,2,id_type)
          enddo
        enddo
!
      endif
      call error_table_clean
      print ('(3a,i2,a,i3,a)'),' File=',trim(sigma_file),' read for ', &
              itypes,' types and ',nlevels,' p-levels or channels'  
!
      end subroutine read_sigmas
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine print_p_and_v (nn1,nn2,p,v,pname,vname)
!
      integer, parameter :: r4=4
      integer :: nn1, nn2
      real(r4) :: p(nn1) 
      real(r4) :: v(nn1,nn2)
      character(len=*) :: pname 
      character(len=*) :: vname 
!
      integer :: n
!
      print *,' '
      print *,'1-d arrays: pname= ',pname,'  vname=',vname
      do n=1,nn1
        print ('(i3,1p1e12.3,2x,5e12.3)'),n,p(n),v(n,:)
      enddo       
!
      end subroutine print_p_and_v 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine print_m (nn1,nn2,v,vname,id)
!
      implicit none
      integer, parameter :: r4=4
      integer :: nn1, nn2, id
      real(r4) :: v(nn1,nn2)
      character(len=*) :: vname 
!
      integer :: n, n2, n2e
      integer :: ns(10)      
!
      print *,' '
      do n2=1,nn2,10
        n2e=min(nn2,n2+9)
        if (id <=5 ) then 
          print ('(3a,i2)'),'array=',vname,'  for id=',id
        else
          print ('(2a)'),'array=',vname
        endif
        do n=1,n2e-n2+1
          ns(n)=n+n2-1
        enddo
        print ('(a,10i12)'),'  n',ns(1:n2e-n2+1)
        do n=1,nn1
          print ('(i3,1p10e12.3)'),n,v(n,n2:n2e)
        enddo
      enddo       
!
      end subroutine print_m
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_eigen (nk,matrix,evects,evalues)
!
! Call library routine to compute eigenvalues and vectors of a symmetric 
! matrix (stored as an array of 8-byte values).
!
      implicit none
      integer, parameter :: rkind3=8
!
      integer,  intent(in)   :: nk
      real(rkind3), intent(in)   :: matrix(nk,nk)
      real(rkind3), intent(out)  :: evects(nk,nk)
      real(rkind3), intent(out)  :: evalues(nk)
!
      integer(4)     :: nwork1  
      integer(4)     :: info 
      integer(4)     :: nk4
      real(rkind3), allocatable :: work1(:) ! must be size >= 3*nk-1
      external dsyev
!
      nk4=nk
      nwork1=3*nk-1
      allocate (work1(nwork1))
      evects=matrix
      call dsyev ( 'V','U',nk4,evects,nk4,evalues,work1,nwork1,info)
      if (info /= 0) then
        print *,'   '
        print *,' * * * * * * * * * * * * * * * * * * * * * * *  '
        print *,' info from dsyev =  ',info
        print *,'   '
      endif
!
      deallocate (work1)  
!
      end subroutine pert_eigen
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_normalize_evects (nk,evects)
!
!  Normalize each eigenvector such that the sum of its squared 
!  components is 1
!
      implicit none
      integer, parameter :: rkind3=8
!
      integer,  intent(in)    :: nk
      real(rkind3), intent(inout) :: evects(nk,nk)
!
      integer  :: k, i
      real(rkind3) :: sum  
     
      do k=1,nk
        sum=0._rkind3
        do i=1,nk
          sum=sum+evects(i,k)*evects(i,k)
        enddo
        sum=1./sqrt(sum)
        evects(:,k)=evects(:,k)*sum
      enddo
!
      end subroutine pert_normalize_evects
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_fix_evalues (nk,evalues)
!
!  Set small (relative to largest values) or negative values to 
!  a small value positive value
!
      implicit none
      integer, parameter :: rkind3=8
!
      integer,  intent(in)    :: nk
      real(rkind3), intent(inout) :: evalues(nk)
!
      integer  :: k
      real(rkind3) :: emax, emin
!
      emax=0.
      do k=1,nk
        if (abs(evalues(k)) > emax) emax=abs(evalues(k))
      enddo
!
      if (rkind3==8) then
        emin=emax*1.e-12
      else
        emin=emax*1.e-6
      endif
!
      do k=1,nk
        if (evalues(k) < emin) evalues(k)=emin
      enddo
!
      end subroutine pert_fix_evalues
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine test_evalues (nlevels,cov,e_vects,e_values,id)
!
      implicit none
      integer, parameter :: r4=4
      integer, parameter :: r8=8
      integer :: nlevels, id
      real(r4) :: cov(nlevels,nlevels)
      real(r8) :: e_vects(nlevels,nlevels)
      real(r8) :: e_values(nlevels)
!
      integer :: i,j,k
      real(r8) :: xm(nlevels,nlevels)   
      real(r8) :: ym(nlevels,nlevels)   
      real(r4) :: zm(nlevels,nlevels)   
!    
      xm=0._r8
      do i=1,nlevels
        do j=1,nlevels
          xm(i,j)=xm(i,j)+e_values(i)*e_vects(j,i)
        enddo
      enddo
!
      ym=0._r8
      do i=1,nlevels
        do j=1,nlevels
          do k=1,nlevels
            ym(i,j)=ym(i,j)+e_vects(i,k)*xm(k,j)
          enddo
        enddo
      enddo
!
      do i=1,nlevels
        do j=1,nlevels
          zm(i,j)=ym(i,j)-cov(i,j)
        enddo
      enddo
!
      if (nlevels <= 10) then 
        call print_m (nlevels,nlevels,zm,'test eigen',id)      
      endif
!     
      end subroutine test_evalues 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine adjust_h_vc (nlevels,h_min,e_vects,h_old)
!
! Set h for vertically or channel correlated errors to a weighted
! mean of the level or channel values provided.  The weights depend on
! the structures of the eigenvectors since the h values are supposed
! to in this correlated case apply to the eigenvectors.  This is not the 
! correct way, but an easy interim way.  Better is to estimate the 
! h-lengths for each eigenvector rather than each channel or levele.
!
      implicit none
      integer, parameter :: r4=4
!
      integer :: nlevels
      real(r4) :: h_min    
      real(r4) :: e_vects(nlevels,nlevels)
      real(r4) :: h_old(nlevels)
!
      integer :: i, j
      real(r4) :: w, w_min, w_sum
      real(r4) :: h_new(nlevels)
!
      w_min=(0.1/nlevels)**2
!
      do j=1,nlevels
        h_new(j)=0.
        w_sum=0.
        do i=1,nlevels
          w=e_vects(i,j)*e_vects(i,j)
          if (w>w_min) then  
            h_new(j)=h_new(j)+w*h_old(i)
            w_sum=w_sum+w
          endif
        enddo
        if (w_sum>0.) then
          h_new(j)=h_new(j)/w_sum
        endif
!
      enddo    
!
! Replace old by new values
      do j=1,nlevels
        h_old(j)=max(h_new(j),h_min)
      enddo
!
      end subroutine adjust_h_vc
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine define_plevels (nlevels,pmax,pmin,plevels)
!
! Define pressure levels for random fields
!
      implicit none
      integer, parameter :: r4=4
      integer, intent(in) :: nlevels
      real(r4), intent(in) :: pmax,pmin
      real(r4), intent(out) :: plevels(nlevels)
!
      integer  :: k
      real(r4) :: pfac
!
      pfac=(pmin/pmax)**(1./real(nlevels-1))
      plevels(1)=pmax
      do k=2,nlevels-1
        plevels(k)=plevels(k-1)*pfac
      enddo 
      plevels(nlevels)=pmin
!
      end subroutine define_plevels
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine read_hcorr_params (dtype,file_eparams,file_path, &
                                    nlevels,pmin,pmax,plevels,    &
                                    corr_func,frac_corr,h_lengths)
!
      implicit none
      integer, parameter :: r4=4
      integer, parameter :: iunit=40
      integer, intent(in) :: nlevels
      real(r4), intent(in) :: pmin, pmax
      real(r4), intent(in) :: plevels(nlevels)
      real(r4), intent(out) :: frac_corr(nlevels)
      real(r4), intent(out) :: h_lengths(nlevels)
      character(len=*), intent(in)  :: dtype
      character(len=*), intent(in)  :: file_eparams(2)   
      character(len=*), intent(in)  :: file_path
      character(len=*), intent(out) :: corr_func
!
      integer :: k,kr,n,nf
      integer :: k_above,k_below
      integer :: nlevf
      integer :: nfiles
      integer, allocatable :: cntf(:,:)
      real(r4), parameter :: frac_min=0.01
      real(r4) :: w1,w2 
      real(r4), allocatable :: plevf(:)
      real(r4), allocatable :: fracf(:,:)
      real(r4), allocatable :: hlenf(:,:)
      character(len=16)  :: cread(4)
      character(len=1)   :: cdum
      character(len=240) :: cfile
      character(len=240) :: file_name
!
!  Read 1 or 2 files of error parameters.
      nfiles=0
      do nf=1,2
        if (trim(file_eparams(nf)) /= 'none') then
          nfiles=nfiles+1
          file_name=trim(file_path)//'/'//trim(file_eparams(nf))
          open (iunit,file=trim(file_name))
!
          read (iunit,'(39x,4(2x,a16))') cread(:)
          read (iunit,'(54x,i5,28x,a)') nlevf,corr_func
          print *,' ' 
          print ('(a,i2,2a)'),'Eparams file',nf,' opened: ',file_name
          print ('(5a)'),' Type on file=',cread(:)
          print ('(a,i3,2a)'),'nlevf=',nlevf,'   corr_func=',corr_func
!
          do n=1,3
            read (iunit,'(a)') cfile
            print ('(a,i2,2x,a)'),'Stats used file',n,trim(cfile)
          enddo
          do n=1,4 
            read (iunit,'(a)') cdum
          enddo
!
          if (nfiles == 1) then 
            allocate (cntf(nlevf,2))
            allocate (plevf(nlevf))
            allocate (fracf(nlevf,2))
            allocate (hlenf(nlevf,2))
          endif
          do n=1,nlevf
            read (iunit,'(4x,f12.2,i8,33x,f7.3,7x,f7.0)') &
                               plevf(n),cntf(n,nfiles), &
                               fracf(n,nfiles),hlenf(n,nfiles)
            print ('(i3,f12.2,i8,2f12.2)'),n,plevf(n),cntf(n,nfiles), &
                               fracf(n,nfiles),hlenf(n,nfiles)
          enddo
          close (iunit)
        endif
      enddo  
!      
! average read params
      if (nfiles == 2) then
        do n=1,nlevf
          cntf(n,1)=(cntf(n,1)+cntf(n,2))/2
          fracf(n,1)=(fracf(n,1)+fracf(n,2))/2.
          hlenf(n,1)=(hlenf(n,1)+hlenf(n,2))/2.
        enddo
      endif
!
! If pmin>pmax, then simply copy values to output argument, as appropriate
! for radiances
      if (pmin > pmax) then 
        frac_corr(1:nlevf)=fracf(1:nlevf,1)
        h_lengths(1:nlevf)=hlenf(1:nlevf,1)
!        
! Otherwise, interpolate to new set of plevels. Only interpolate from 
! levels with relevant data (fracf(k,1) > frac_min, thereby excluding 
! levels with undetermined or statistically unreliable estimates.  
! Interpolation is linear in p. Extrapolation below or above top and bottom
! assumes zero gradient. 
      else
!
        do n=1,nlevels
          k_above=0
          k_below=0
          do k=1,nlevf
!
! Determine levels read with valid data above and below desired level n
            if (fracf(k,1) > frac_min .and. plevels(n) < plevf(k)) then
              k_below=k
            endif
            kr=nlevf+1-k  ! reverse order of search
            if (fracf(kr,1) > frac_min .and. plevels(n) >= plevf(kr)) then
              k_above=kr
            endif
          enddo
!
          if (k_above == 0 .and. k_below > 0 ) then   
            frac_corr(n)=fracf(k_below,1)
            h_lengths(n)=hlenf(k_below,1)
          else if (k_above > 0 .and. k_below == 0 ) then  
            frac_corr(n)=fracf(k_above,1)
            h_lengths(n)=hlenf(k_above,1)
          else if (k_above > 0 .and. k_below > 0 ) then               
            if (k_above > k_below) then ! value can be interpolated
              w1=(plevels(n)-plevf(k_above))/(plevf(k_below)-plevf(k_above))
              w2=1.-w1
              frac_corr(n)=w1*fracf(k_below,1)+w2*fracf(k_above,1)
              h_lengths(n)=w1*hlenf(k_below,1)+w2*hlenf(k_above,1)
            else
              frac_corr(n)=fracf(k_above,1)
              h_lengths(n)=hlenf(k_above,1)
            endif
          endif
!
        enddo
      endif
!
      if (nfiles > 0) then      
        deallocate (cntf,plevf,fracf,hlenf)
      endif
!        
      end subroutine read_hcorr_params 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine adjust_hcorr_params (nlevels,fix_params,frac_corr,x_lengths)
!
      implicit none
      integer, parameter :: r4=4
      integer, intent(in) :: nlevels
      real(r4), intent(in) :: fix_params(4)
      real(r4), intent(inout) :: frac_corr(nlevels)
      real(r4), intent(inout) :: x_lengths(nlevels)
!
      integer :: n
!
      do n=1,nlevels
        frac_corr(n)=frac_corr(n)*fix_params(1)+fix_params(2)
        x_lengths(n)=x_lengths(n)*fix_params(3)+fix_params(4)
      enddo
!
      end subroutine adjust_hcorr_params 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine adjust_matrix (id,itypes_corr,d_type,nlevels, &
                                frac_corr,sigma,cov)
!    
!  Adjust cov matrix so that it accounts for frac_corr factors and 
!  is a true cov, such that implied correlations are bounded. 
!  Here, those bounds are set to a little less than 1. 
!
!  This can be considerded as equivalent to constructing a new channel or
!  vert correl matrix that applies only to the part of the error modeled
!  as being correlated. 
!
      implicit none
      integer, parameter :: r4=4
!
      integer, intent(in) :: id,itypes_corr
      integer, intent(in) :: nlevels
      real(r4), intent(in) ::frac_corr(nlevels) 
      real(r4), intent(in) :: sigma(nlevels) 
      real(r4), intent(inout) :: cov(nlevels,nlevels)
      character(len=*), intent(in) :: d_type
!
      logical, parameter :: lwrite_corr=.true.  ! T if corrs written to 2nd file
      integer :: n,k
      integer :: count_fix(2)
      real(r4), parameter :: corr_max_conv=0.95
      real(r4), parameter :: corr_max_rad=0.95
      real(r4) :: corr_max
      real(r4) :: xfac,xfac1,xfac2
      real(r4) :: cnew(nlevels,nlevels)
      real(r4) :: corr(nlevels,nlevels)
!
! Compute new trial COV accounting for frac_corr factors
      do k=1,nlevels
        do n=1,nlevels
          xfac1=sqrt(max((1.-frac_corr(k))*(1.-frac_corr(n)),1.e-20)) 
          xfac2=sqrt(max(frac_corr(k)*frac_corr(n),1.e-20)) 
          xfac=xfac1+xfac2
          cnew(k,n)=cov(k,n)/xfac
        enddo
      enddo
!
! Compute trial correlations from trial covariances
      do k=1,nlevels
        do n=1,nlevels
          xfac1=cnew(n,n)*cnew(k,k)
          if (xfac1 > 0.) then
            xfac2=1./sqrt(xfac1)
          else
            xfac2=0.
          endif
          corr(k,n)=xfac2*cnew(k,n)
        enddo
      enddo
!
      if (trim(d_type) == 'WIND' .or. trim(d_type) == 'MASS' &
          .or. trim(d_type) == 'PREPBUFR') then
        corr_max=corr_max_conv
      else 
        corr_max=corr_max_rad
      endif
!
! Adjust correlations so that cnew is a true cov matrix
! (i.e., correlations are bounded).
      count_fix(:)=0. 
      do k=1,nlevels
        do n=1,nlevels
          xfac1=sqrt(max(cnew(n,n)*cnew(k,k),1.e-20))
          if (k /= n) then
            count_fix(1)=count_fix(1)+1
            if (abs(corr(k,n)) > corr_max) then
              count_fix(2)=count_fix(2)+1
              if (corr(k,n) > 0.) then
                corr(k,n)=corr_max
              else
                corr(k,n)=-corr_max
              endif
              cnew(k,n)=xfac1*corr(k,n)
            endif
            cov(k,n)=cnew(k,n)
          endif
        enddo
      enddo
!
! Optional write corr matrices to a separate file for diagnostic use
      if (lwrite_corr) then 
        if (id == 1) then
          open (92,file='corr_matrix_1',form='unformatted')
        endif
        write (92) corr(:,:)
        if (id == itypes_corr) then
          close (92)
          print *,' '
          print *,'File corr_matrix_1 written for diagnostic use'
          print *,'number of matrices written =',id
        endif
      endif
!
      print *,' '
      print *,'Matrix adjusted'
      print ('(a,i8)'),' number of non-diagonal matrix elements = ',count_fix(1)
      print ('(a,i8)'),' number of adjusted matrix elements = ',count_fix(2)
      print ('(a,f8.5)'),' max absolute correl value allowed = ',corr_max
!
      end subroutine adjust_matrix 

                

