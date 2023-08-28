   module m_obs_error_table
!
!  Module for reading and and computing error parameters from 3 files: 
!  resource file, table of standard deviations, and correlation parameters
!  Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use m_kinds, only : rkind1, rkind2
!
   implicit none
   private
   public :: error_table_setup
   public :: error_table_clean
   public :: error_table_find_corr_id
   public :: error_table_find_stdv_id 
   public :: error_table_read_stdv 
   private :: error_table_read_rc
   private :: error_table_read_corr_params 
!
   integer, parameter :: r4=4 ! precision of real variables on et_file_err_corr
!
   logical, public :: et_l_vc_corr   ! true if also channel or vert correl
!
   integer, public :: et_n_err1, et_n_err2, et_n_err3
   integer, public :: et_vcorr_dist_num
   integer, public :: et_icolumn_tq(2)
   integer, public :: et_icolumn_uv(2)
   integer, public :: et_icolumn_gpsro
   integer, public :: et_icolumn_rad
   integer, public :: et_icolumn_ps
   integer, public :: et_nlevels ! number of p-levels, channels, or principal   
                                 ! components to create for random fields 
   integer, public :: et_itypes_corr ! number of distinct sets of random fields 
   integer, public :: et_ran_fields_kinds ! =2 so u,v differ; =1 otherwise
   integer, public :: et_rf_nmax     ! spectral truncation of horiz correl func
   integer, public :: et_hcorr_ks(2) ! index and number of valid ks ids
   integer, public, allocatable :: et_err_itype(:)
   integer, public, allocatable :: et_ks_list(:,:) ! list of ks ids per itype 
!
   real(rkind2), public :: et_pert_fac
   real(rkind2), public :: et_vcorr_dist(10) ! params for vert corr dist
   real(rkind2), public, allocatable :: et_err_tab(:,:,:) ! stdvs of added errs
   real(r4), public :: et_pmax               ! if conv obs, max p-level 
   real(r4), public :: et_pmin               ! if conv obs, min p-level
   real(r4), public, allocatable :: et_e_vects(:,:,:) ! eofs of vert or chan cov
   real(r4), public, allocatable :: et_e_v_sqrt(:,:)  ! sqrt of e-values of cov
   real(r4), public, allocatable :: et_hcorr_lengths(:,:) ! horiz corr Len (km)
   real(r4), public, allocatable :: et_frac_corr(:,:)     ! frac of tot var correl
   real(r4), public, allocatable :: et_cov_matrix(:,:,:)  ! channel covariances
!
   character(len=240), public :: et_file_err_stdv
   character(len=240), public :: et_file_err_corr
   character(len=20), public, allocatable :: et_corr_shapes(:)
   character(len=*), parameter :: myname='m_obs_err_table'
!
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_setup (lprint,dtype,file_rc,lerrtable, &
                                 lcorrfld,iseed,ierr) 
!
!  Read error (1) resource file, (2) stdv table, (3) correl params
!
   implicit none
!
   logical, intent(in) :: lprint
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: file_rc
!
   logical, intent(out) :: lerrtable
   logical, intent(out) :: lcorrfld
   integer, intent(out) :: iseed 
   integer, intent(out) :: ierr
!
   integer :: ier1
!
   ierr=0
!
! Read resource file
   call error_table_read_rc (file_rc,iseed,lprint,dtype,lerrtable, &
                             lcorrfld,ier1)
   if (ier1 /= 0) then
     print *,' ERROR in attempt to read resorce file: ',ier1, &
             trim(file_rc)
     ierr=ier1+ierr
   endif 
!
! Read error stdv table file
   if (ierr == 0 .and. lerrtable) then 
     call error_table_read_stdv (lprint,dtype,ier1)
     if (ier1 /= 0) then
       print *,' ERROR in attempt to read error table stdv file: ',ier1
       ierr=ier1+ierr
     endif
   endif 
!
! Read error correlation parameters file
   if (ierr == 0 .and. lcorrfld) then
     call error_table_read_corr_params (lprint,dtype,ier1)
     if (ier1 /= 0) then
       print *,' ERROR in attempt to read error corr params file: ',ier1
       ierr=ier1+ierr
     endif
   endif 
!
   end subroutine error_table_setup
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_clean
!
!  Deallocate error table arrays allocated in error_table_read_stdv
!  or in error_table_read_corr_params
!
   if (allocated(et_err_tab))       deallocate (et_err_tab)
   if (allocated(et_err_itype))     deallocate (et_err_itype) 
   if (allocated(et_hcorr_lengths)) deallocate (et_hcorr_lengths)
   if (allocated(et_frac_corr))     deallocate (et_frac_corr)
   if (allocated(et_e_vects))       deallocate (et_e_vects)
   if (allocated(et_e_v_sqrt))      deallocate (et_e_v_sqrt)
   if (allocated(et_corr_shapes))   deallocate (et_corr_shapes)
   if (allocated(et_ks_list))       deallocate (et_ks_list)
   if (allocated(et_cov_matrix))    deallocate (et_cov_matrix)
!
   end subroutine error_table_clean 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_fix 
!
! Change some units and large values in table of error standard deviations.
! Values of 0.1e10 at unobserved levels will produce problems when used to
! create vertical correlation matrix.
!
   implicit none
!
! The following are set to signify "large" errors
   real(rkind1), parameter :: T_err=4.
   real(rkind1), parameter :: q_err=0.05   ! % rel. humidity
   real(rkind1), parameter :: w_err=5.     ! each wind component
   real(rkind1), parameter :: p_err=200.   ! Pa
   integer :: i,j
!
! Change units of p from hPa to Pa
   do i=1,et_n_err3
     do j=1,et_n_err1
       et_err_tab(j,1,i)=100.*et_err_tab(j,1,i)    ! p levels     
       if (et_err_tab(j,5,i) > 1.e8) then       ! replace large p errors
         et_err_tab(j,5,i)= p_err
       else
         et_err_tab(j,5,i)=100.*et_err_tab(j,5,i)  ! stdv of p error  
       endif 
     enddo
   enddo
!
! replace values of 0.1e10 read in table by reasonably large values 
   do i=1,et_n_err3
     do j=1,et_n_err1
       if (et_err_tab(j,2,i) > 1.e8) then   ! replace large T errors
         et_err_tab(j,2,i)= T_err
       endif
       if (et_err_tab(j,4,i) > 1.e8) then   ! replace large wind errors
         et_err_tab(j,4,i)= w_err
       endif
     enddo
   enddo
!
! Change q error to fraction of rh
! Also, don't allow large errors at very low p where rh computed poorly
   do i=1,et_n_err3
     do j=1,et_n_err1
       if (et_err_tab(j,3,i) > 1.e8) then      ! replace large q errors
         et_err_tab(j,3,i)= q_err
       else
         et_err_tab(j,3,i)=0.1*et_err_tab(j,3,i)  ! rescale error to fraction rh  
       endif
       if (et_err_tab(j,1,i) < 3.e4) then   ! reduce values at low pressure
         et_err_tab(j,3,i)=q_err*(et_err_tab(j,1,i)/3.e4)**3
       endif            
     enddo
   enddo
! 
   end subroutine error_table_fix 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_read_rc (file_rc,iseed,lprint,dtype, &
                                   lerrtable,lcorrfld,ierr)
!
! Read resource file for obtaining error parameters and error file names
!
   implicit none      
   logical, intent(in)  :: lprint
   logical, intent(out) :: lerrtable
   logical, intent(out) :: lcorrfld 
   integer, intent(out) :: ierr
   integer, intent(out) :: iseed
   character(len=*), intent(in) :: file_rc
   character(len=*), intent(in) :: dtype
!
   logical, parameter :: lstop=.false. ! must be false when mpi is used
   logical :: found
!
   integer, parameter  :: iunit=40
   integer :: i, n, n1, n2
   integer :: nrecs, ios
   integer :: format_header
   integer :: iseed_exp, iseed_dtype
   integer(1) :: n_lines
   character(len=30)  :: read_name
   character(len=30)  :: cdum
   character(len=260) :: read_info(100)
   character(len=*), parameter :: mysub=myname//'::read_error_rc'
   character(len=1) :: mysub0
!
   ierr=0
!
! Set msub0 to eliminate error handeling in sub=find_name_2; Instead errors 
! are handeled here so that defaults are set without error indicated 
   mysub0(1:1)=' '
!
! Open rc file
   open(unit=iunit,file=trim(file_rc),form='formatted',status='old',iostat=ios)
   if (ios /= 0) then 
     ierr=99
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open error.rc file for iunit=', &
                     iunit,' iostat=',ios,' and file name=',trim(file_rc)
     return  
   elseif (lprint) then
     print *,' '
     print ('(3a,i4)'),' Rc file=',trim(file_rc),' opened on unit=',iunit
   endif
!  
! Read format indicators
   read (iunit,*) cdum,format_header
!
! Read common header  (i.e., portion of header valid for all obs types)
   do n1=1,100 
     read (iunit,'(a)') read_info(n1)
     if (read_info(n1)(1:3) == '---') exit
     nrecs=n1
   enddo     
!
! Look for value of common random seed for exp
   call find_name_2 (nrecs,read_info,lstop,mysub0,'seed_for_exp=',n)
   if (n == 0) then
     if (lprint) then
       print *,' '
       print *,'ERROR reading RC file:'
       print *,'seed_for_exp not found in common portion of header'
     endif
     ierr=ierr+200
     return
   else 
     read (read_info(n),*) cdum,iseed_exp
   endif
!
! Find parameters for desired dtype
   found=.false.
   do n1=1,1000
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
     ierr=ierr+300
     if (lprint) then 
       print *,' '
       print *,'ERROR reading RC file:'
       print *,'info for dtype=',trim(dtype),' not found'
     endif
     return
   endif
!
!  Look for value of pert_fac
   call find_name_2 (nrecs,read_info,lstop,mysub0,'pert_fac=',n)
   if (n == 0) then
     et_pert_fac=1.
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR pert_fac'
     endif
   else 
     read (read_info(n),*) cdum,et_pert_fac
   endif
!
!  Look for value of iseedtype
   call find_name_2 (nrecs,read_info,lstop,mysub0,'seed_data_type=',n)
   if (n == 0) then
     call error_table_default_seed (iseed_dtype,dtype)
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR seed_data_type'
     endif
   else 
     read (read_info(n),*) cdum,iseed_dtype
   endif
!
!  Look for value of vertical correlation lengths
   call find_name_2 (nrecs,read_info,lstop,mysub0,'vcorr_dist=',n)
   if (n == 0) then
     if (trim(dtype) == 'PREPBUFR') then
       et_vcorr_dist_num=4
       et_vcorr_dist(:)=1.
     elseif (trim(dtype) == 'SATWIND') then 
       et_vcorr_dist_num=2
       et_vcorr_dist(:)=1.
     elseif (trim(dtype) == 'GPSRO') then 
       et_vcorr_dist_num=1
       et_vcorr_dist(:)=1.
     else
       et_vcorr_dist_num=0
       et_vcorr_dist(:)=0.
     endif
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR vcorr_dist'
     endif
   else 
     read (read_info(n),*) cdum,et_vcorr_dist_num
     read (read_info(n),*) cdum,n2,et_vcorr_dist(1:et_vcorr_dist_num)
   endif
!
!  Look for value of et_file_err_stdv
   call find_name_2 (nrecs,read_info,lstop,mysub0,'file_err_var=',n)
   if (n == 0) then
     et_file_err_stdv='NONE'
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR file_err_stdv name'
     endif
   else 
     read (read_info(n),*) cdum,et_file_err_stdv
   endif
!
!  Look for value of et_file_err_corr
   call find_name_2 (nrecs,read_info,lstop,mysub0,'file_err_corr=',n)
   if (n == 0) then
     et_file_err_corr='NONE'
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR file_err_correlations name'
     endif
   else 
     read (read_info(n),*) cdum,et_file_err_corr
   endif
!
   if (trim(et_file_err_stdv) == 'none' .or. &
       trim(et_file_err_stdv) == 'NONE') then
     lerrtable=.false.
   else
     lerrtable=.true. 
   endif
!
   if (trim(et_file_err_corr) == 'none' .or. &
       trim(et_file_err_corr) == 'NONE') then
     lcorrfld=.false.
   else
     lcorrfld=.true. 
   endif
! 
! Set random seed based on seeds for experiment and obs type 
! (This will later be augmented by the datetime) 
   iseed=iseed_exp+iseed_dtype
!
! Print some information read from the rc file
   if (lprint) then
     print ('(a)'),' Values read from error.rc file:'
     print ('(a,f7.3,14x,a)'),'   pert_fac=',  &
                 et_pert_fac,' (factor to multiply table stdv by)'
     print ('(3(a,i7))'),'   random_seed=',       &
                 iseed_exp,' (exp) + ',iseed_dtype,' (dtype) = ',iseed
     if (et_vcorr_dist_num > 0) then 
       print ('(a,10f7.3)'),'   vert_corr_dist=', &
                            et_vcorr_dist(1:et_vcorr_dist_num)
     else
       print ('(a)'),'   vert_corr_dist= none'
     endif
     print ('(2a)'),'   file_err_stdv = ',trim(et_file_err_stdv)
     print ('(2a)'),'   file_err_corr= ',trim(et_file_err_corr)
   endif
! 
   end subroutine error_table_read_rc 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_default_seed (iseed,dtype)
!
!  Create a default value of iseed_dtype based on the characters in 
!  the obs type name dtype 
!
   implicit none
   integer :: iseed
   character(len=*) :: dtype
!
   integer :: len1
   integer :: m,n
   character(len=26) :: alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
   len1=len(trim(dtype))
   iseed=len1
   do m=1,len1
     do n=1,26
       if (dtype(m:m) == alphabet(n:n)) then
         iseed=iseed+n
       endif 
     enddo
   enddo
!
   end subroutine error_table_default_seed   
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_read_stdv (lprint,dtype,ierr)
!
! Read table of error standard deviations for desired data type
!
   use m_sat_info_table, only :  sat_info_table_get_1i
   use m_sat_info_table, only :  sat_info_table_get_1r
!
   implicit none
   logical :: lprint
   integer :: ierr
   character(len=*) :: dtype
!
   integer, parameter :: iunit=10
   integer :: i, j, j1
   integer :: ier, ios
   integer :: ic, inum, i1, i2
   integer :: ihead   ! number of header records to skip in sat_error file
   integer :: igroups ! number of sat platforms or instruments in table
   integer :: n_channels
   real(rkind1), allocatable :: err_tab_2(:) 
   character(len=18) :: adum
   character(len=16) :: i_type
   character(len=16) :: sat_name
   character(len=1)  :: adum1
!
   ierr=0
!
! Set index to denote in which column of the table desired stdv value  is 
   et_icolumn_tq=(/2,3/)
   et_icolumn_uv=(/4,4/)
   et_icolumn_gpsro=2
   et_icolumn_rad=2
   et_icolumn_ps=1 
!
! Allocate array for table of errors
   if (trim(dtype) == 'PREPBUFR' .or. trim(dtype) == 'SATWIND') then
     et_n_err1=33      ! max number of p-levels in table
     et_n_err2=6       ! number of fields in table
     et_n_err3=200     ! max number of sub-types of observations in table
   elseif (trim(dtype) == 'GPSRO') then
     et_n_err1=102
     et_n_err2=13        
     et_n_err3=1        
   else   ! radiance obs type
     call sat_info_table_get_1i ('dtype','nchan',dtype,et_n_err1,ier)
     et_n_err2=2     
     et_n_err3=8       ! max number of sub-types of observations in table
   endif
!
   allocate (et_err_itype(et_n_err3))
   allocate (et_err_tab(et_n_err1,et_n_err2,et_n_err3))
   allocate (err_tab_2(et_n_err2))
!
! Read table of obs error standard deviations for conventional obs.
! The columns are arranged as: (1) p-level, (2) T, (3) r.h., (4) wind, 
! (5) ps; (6) ?
   open(unit=iunit,file=trim(et_file_err_stdv),form='formatted',status='old',iostat=ios)
   if (ios /= 0) then 
     ierr=98
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open error_stdv file for iunit=', &
                    iunit,' iostat=',ios,' and file name=',trim(et_file_err_stdv)
     return  
   elseif (lprint) then
     print *,' '
     print ('(3a,i4)'),' Input file=',trim(et_file_err_stdv),' opened on unit=',iunit
   endif
!
   if (trim(dtype)=='PREPBUFR') then
     do i=1,et_n_err3  ! loop over observation types
       read (iunit,'(i4)') et_err_itype(i)
       do j=1,et_n_err1
         read (iunit,'(1x,6e12.5)') err_tab_2(:)
         et_err_tab(j,:,i)=err_tab_2(:)
       enddo
     enddo
!
   elseif (trim(dtype)=='GPSRO') then
     do i=1,et_n_err3  ! loop over observation types (maybe only 1 here) 
       read (iunit,'(i4)') et_err_itype(i)
       read (iunit,'(i3,f9.1,12f9.2)') ic,err_tab_2(1),err_tab_2(2:et_n_err2)
       et_err_tab(1,:,i)=err_tab_2(:)
       do j=2,et_n_err1
         read (iunit,'(i3,f9.1,12e9.2)') ic,err_tab_2(1),err_tab_2(2:et_n_err2)
         et_err_tab(j,:,i)=err_tab_2(:)
       enddo
     enddo
!
   else  ! prepare table for satellites
!        
     ic=0
     et_err_itype(:)=0
!
! Read header: 
! igroups=number of distinct
     read (iunit,'(7x,i2,30x,i2)') ihead,igroups  
     do i=1,ihead                     ! skip these header records
       read (iunit,'(a1)') adum1
     enddo
!
! Check the following records and skip or read as required
     do i=1,igroups   
       read (iunit,'(i2)') inum
       read (iunit,'(1x,a16,1x,a16,1x,i4)') i_type,sat_name,n_channels
       if (trim(i_type) == trim(dtype)) then  ! interpret required records
         ic=ic+1
         if (ic > et_n_err3  .or. n_channels > et_n_err1)  then
           if (lprint) then
             print *,' '
             print *,'Problem reading sat_error file'
             print ('(i3,a,i3)'),ic,'=ic must not be > n_err3=',et_n_err3
             print ('(i3,2a,i3)'),n_channels,'=n_channels must not be', &
                     ' > n_err1=', et_n_err1
           endif
           ierr=1
         endif   
         do j=1,n_channels
           read (iunit,'(1x,a18,2i5,f7.3)') adum,i1,i2,et_err_tab(j,2,ic)
           et_err_tab(j,1,ic)=j
         enddo
!
! Set err_itype to either sat ID or to instrument ID numbers
! (Get values from table)
         if (trim(dtype) == 'AIRS' .or. trim(dtype) == 'AMSUAAQUA') then
           call sat_info_table_get_1i ('dtype','siid',dtype, &
                                       et_err_itype(ic),ier) 
         else
           call sat_info_table_get_1i ('platform','said',sat_name, &
                                       et_err_itype(ic),ier) 
         endif
       else  ! skip these records
         do j=1,n_channels
           read (iunit,'(a1)') adum1  
         enddo
       endif  ! test on dtype  
     enddo    ! loop over all groups of entries in sat error table      
!
   endif      ! test on dtype 
   close (iunit)
!
! Change some values in the error table
   if (trim(dtype) == 'PREPBUFR') then
     call error_table_fix 
   endif
!
   deallocate (err_tab_2)
!
   if (lprint) then
     print *,' Table of error standard deviations read'
   endif
!
   end subroutine error_table_read_stdv 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine error_table_read_corr_params (lprint,dtype,iret)
!
! Read file containing parameters and statistics required for specifying 
! the random fields.  Separate files are required for each data type.
! (All onventional observations are handeled by a single file with up to
! 5 different sets of parameter specifications.)
!
   implicit none
!
   logical :: lprint
   integer, intent(out) :: iret              ! return code 
   character(len=*), intent(in) :: dtype     ! name of type of obs data
!
   integer, parameter :: iunit=33
   integer :: file_format          ! index for future changes to file_format
   integer :: i, k
   integer :: k1, k_max, ios
   real(r4) :: version_number
   real(r4), allocatable :: v_lengths_i(:)    ! vert correl lengths each field
   real(r4) :: e_max  
   real(r4), allocatable :: e_values(:,:)     ! eigenvalues of cov matrix
   character(len=16)  :: dtype_file
   character(len=140) :: c_info       ! description of what is on file
   character(len=140) :: c_file_orig  ! file used to create param file
!
   iret=0
!
! Open optional file of correlation parameters
! (This file name is supplied in the error.rc file)
   open (iunit,file=trim(et_file_err_corr),form='unformatted',status='old',iostat=ios)
   if (ios /= 0) then 
     iret=97
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open error_corr file for iunit=', &
                    iunit,' iostat=',ios,' and file name=',trim(et_file_err_corr)
     return  
   elseif (lprint) then 
     print *,' File=',trim(et_file_err_corr),' opened to read correlation parameters'
   endif
!
! Read number of subtype groups, data-type, and format version number
   read (iunit) et_itypes_corr,dtype_file,file_format
!
   if (trim(dtype_file) == 'MASS' .or. trim(dtype_file) == 'WIND') then
     if (lprint) then
       print ('(3a)'),'dtype on file is ',trim(dtype_file), &
                      ' is considered as dtype PREPBUFR'
     endif
     dtype_file='PREPBUFR'  ! treat as this data type
   endif
!  
! Check header    
   if (et_itypes_corr <= 0) then
     iret=-1 
     if (lprint) print *,' No horiz correlation information on file'
     return
   elseif (trim(dtype) /= trim(dtype_file)) then
     iret=-2
     if (lprint) print *,' Requested data type = ',trim(dtype), &
                         ' differs from data type = ', &
                         trim(dtype_file),' on hcorr parameter file'
     return
   endif
!
! Read info on file (num=record number):
! (2) info specifing this version of error parameters
! (3) number of levels (or channels), pressure range required to cover all 
!     data if levels, and flag whether vert or chan correl should also be 
!     considered 
! (4) vertical correlation lengths for each subtype group
! (5) ks or sat id list of each subgroup (up to 10 in ech group)
   allocate (v_lengths_i(et_itypes_corr)) 
   allocate (et_ks_list(10,et_itypes_corr))  ! only upto 10 subtypes in a group 
   read (iunit) version_number,c_info,c_file_orig   
   read (iunit) et_nlevels,et_pmax,et_pmin,et_l_vc_corr
   read (iunit) v_lengths_i
   read (iunit) et_ks_list
   et_pmax=et_pmax*100.  ! change mb to Pa
   et_pmin=et_pmin*100.  ! change mb to Pa
!
! Print some info read
   if (lprint) then
     print ('(a,i2)'),'  Number of correlation parameter sets = ',et_itypes_corr
     print ('(a,f6.2)'),'  Version number for file format = ',version_number
     print ('(2a)'),'  Info=',trim(c_info)
     print ('(2a)'),'  Orig_file=',trim(c_file_orig)
     print ('(a,i3)'),'  Number of levels or channels = ', et_nlevels
     print ('(a,l8)'),'  Vertical or channel correlations also considered = ', &
                        et_l_vc_corr
     print ('(a,2f12.3)'),'  pmax, pmin =',et_pmax,et_pmin 
     print ('(a,5f7.1)'),'  v_lengths(id) (meters) =',v_lengths_i(1:et_itypes_corr)
     if (et_pmax > et_pmin) then 
       do i=1,et_itypes_corr
         print ('(a,i2,a,10i6)'),'  ks_list(*,',i,'):',et_ks_list(:,i)
       enddo 
     endif
   endif
   deallocate (v_lengths_i) ! This is not used but is for information only
!
! (6) spectral truncation for horizontal random fields to be produced (should
!     be enough to very well resolve the correlation scale)
! (7) correlation shape requested for each subtype group
   allocate (et_corr_shapes(et_itypes_corr))
   if (file_format > 1) then
     read (iunit) et_rf_nmax
     read (iunit) et_corr_shapes
   else
     et_corr_shapes(:)='WHITE'
     et_rf_nmax=180*3
   endif
   if (lprint) then
     print ('(a,i4)'),' Spectral truncation for random fields = ',et_rf_nmax 
     print ('(a,10(1x,a))'),                                    &
         ' Correlation function shapes for random field sets:', &
         (trim(et_corr_shapes(i)),i=1,et_itypes_corr)
     if (file_format < 2) then
       print ('(a)'),' Shape and truncation above were set as default values'
     endif 
   endif
!
! Read parameters that define correlations as functions of level or channel
! (8) fractions of variance that is correlated for each level or channel and 
!     for subtype group
! (9) like (8) but for horiz. correlation lengths
   allocate (et_frac_corr(et_nlevels,et_itypes_corr)) 
   allocate (et_hcorr_lengths(et_nlevels,et_itypes_corr)) 
   read (iunit) et_frac_corr
   read (iunit) et_hcorr_lengths
!
! Read information about vertical level or channel correlations if on file
! (10) covariance matrix from which EOFs determined
! (11) eigenvalues of covariance matrix (variance explained by each EOF)
! (12) matrix whose columns are EOFS 
   if (et_l_vc_corr) then 
     allocate (et_cov_matrix(et_nlevels,et_nlevels,et_itypes_corr)) 
     allocate (e_values(et_nlevels,et_itypes_corr)) 
     allocate (et_e_v_sqrt(et_nlevels,et_itypes_corr)) 
     allocate (et_e_vects(et_nlevels,et_nlevels,et_itypes_corr)) 
     read (iunit) et_cov_matrix
     read (iunit) e_values
     read (iunit) et_e_vects
     do i=1,et_itypes_corr
       do k=1,et_nlevels
         et_e_v_sqrt(k,i)=sqrt(e_values(k,i))
       enddo
     enddo
     deallocate (e_values)
   endif
!
   close (iunit)
!
! Print table of some level or channel dependent inputs
! In the following 'peak-level' is the level or channel at which the 
! corresponding eigenvector has a maximum amplitude
   if (lprint) then 
     if (et_l_vc_corr) then 
       do i=1,et_itypes_corr
         print ('(2a,i2)'),' Table of h_lengths, frac_corr,', &
                           ' and sqrt(e_values) and peak-level for itype=',i
         do k=1,et_nlevels
           e_max=0.
           k_max=0
           do k1=1,et_nlevels
             if (e_max < abs(et_e_vects(k1,k,i)) ) then
               k_max=k1
               e_max=abs(et_e_vects(k1,k,i))
             endif
           enddo
           print ('(a,i4,f10.1,f10.4,1p1e14.3,i5)'),'k=',k, &
                  et_hcorr_lengths(k,i),et_frac_corr(k,i),  &
                  et_e_v_sqrt(k,i),k_max
         enddo
       enddo
     else  ! no arrays e_values or e_vects have been allocated
       do i=1,et_itypes_corr
         print ('(2a,i2)'),' Table of h_lengths and frac_corr', &
                            ' for itype=',i
         do k=1,et_nlevels
           print ('(a,i4,f10.1,f10.4)'),'k=',k, &
                 et_hcorr_lengths(k,i),et_frac_corr(k,i)
         enddo
       enddo
     endif
!
     print *,' All information on file hcorr_stats.bin read'
   endif  ! check on lprint
!
   end subroutine error_table_read_corr_params
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine error_table_find_corr_id (itype,hcorr_id)
!
! Determine index for the set of horizontal correlation parameters that 
! should be used, if any. Also determine number of subtypes included in
! the same group of ks (The latter is required for rotating the random 
! field to reduce the correlation among different subtypes in each group). 
!
   implicit none
!
   integer, intent(in)  :: itype       ! subtype of data (e.g., kx value)
   integer, intent(out) :: hcorr_id    ! group in which subtype found 
!
   integer :: i, k, n 
!
! Find first occurance of value itype in ks_list
   hcorr_id=0  ! default indicates no horizontal correl for this type
   if (et_itypes_corr > 0) then
     do k=1,et_itypes_corr
       if (hcorr_id == 0) then ! ks not found in list yet 
         do i=1,10 
           if (itype == et_ks_list(i,k) .or. et_ks_list(i,k) == 9999 ) then
             hcorr_id=k
             et_hcorr_ks(1)=i ! index of found ks in list for type index k
           endif 
         enddo
       endif  
     enddo
   endif
!
   et_hcorr_ks(2)=0  ! default indicates no ks values in list
   if (hcorr_id > 0) then
     do i=1,10
       if (et_ks_list(i,hcorr_id) > 0) then
         et_hcorr_ks(2)=i   ! num of valid ks values in list for type index k
       endif
     enddo
   endif
!
   end subroutine error_table_find_corr_id
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine error_table_find_stdv_id (itype,err_index)
!
! Find index of desired subtype in list of subtypes in extracted error 
! tables
!
      implicit none
!
      integer, intent(in) :: itype
      integer, intent(out) :: err_index
      integer :: k  
!     
      err_index=0
      do k=1,et_n_err3
        if (itype == et_err_itype(k)) then
          err_index=k
        endif
      enddo
!
   end subroutine error_table_find_stdv_id 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   end module m_obs_error_table
