   module m_kx_table
!
! Read/write table of observation and satellite types, data counts,
! and probability function params                   
!
   use m_kinds, only : rkind1
   implicit none
!                                                     
   private
   public kx_table_read
   public kx_table_write
   public kx_table_clean
   public kx_table_sum
   public kx_histogram_print
   public kx_latlon_range 
!
   logical, allocatable, public :: kx_present(:)
! 
   integer, public :: kx_cbins   ! num of bins for distribution histograms
   integer, public :: kx_pbins   ! num of pressure layers for indpend. histogram
   integer, public :: kx_jbins_lats ! num of lat ranges for indpend. histograms
   integer, public :: kx_jbins_lors ! separations into indep bins for land, sea
   integer, public :: kx_jbins_time ! num of separate time slots for ind. histograms
   integer, public :: kx_jbins      ! kx_jbins_lats+lors+time  
   integer, public :: kx_nparams    ! num of prob params to determine 
   integer, public :: kx_num        ! num of independ pairs of kxs and sats
   integer, public :: kx_ipw_nlevs  ! num of layers of ipw to consider           
   integer, public :: kx_ipw_types  ! num of sets of scale and qmin values
   integer, public :: kx_field_imax ! num of lons in grid interp from nr grid
   integer, public :: kx_field_jmax ! num of lats in grid interp from nr grid
   integer, public :: kx_field_kdim ! num of levs in grid interp from nr grid
   integer, public :: kx_field_slots  ! num of interp grid per ana period
   integer, public :: kx_times_satloc ! num of ana periods used to comp satloc
   integer, public :: kx_times_params ! num of ana periods used to comp params
!
! filter values:
!   flat is a fixed lat to use for determining i_stride
!   slev refers to the min solar elevation permitted for this obs type
!   torb is a pressure which indicates that for cloud tops below this level,
!     the cloud bottom should be used for the observation level; otherwise the
!     cloud top level should be used to define the obs level.
!   mcnt is the min number of obs in a satloc region for that region to be
!     considered viewed by that kx at an NR time.
   integer, public :: kx_nfilters   ! num of filtering parameters
   integer, public :: kx_filt_flat  ! index for lat used to set i-stride
   integer, public :: kx_filt_slev  ! index for min solar elevation
   integer, public :: kx_filt_torb  ! index for torb pressure cutoff
   integer, public :: kx_filt_mcnt  ! index for min counts in satloc bins
!
   integer, allocatable, public :: kx_list(:)       ! kx values
   integer, allocatable, public :: kx_satid(:)      ! WMO sat id number
   integer, allocatable, public :: kx_ncnum(:)      ! NC number in message header
   integer, allocatable, public :: kx_ijgrid(:)     ! type of grid pattern     
   integer, allocatable, public :: kx_nwcsm(:)      ! value of NCWSM BUFR variable
   integer, allocatable, public :: kx_ifunc(:)      ! id for form of prob function
   integer, allocatable, public :: kx_subtype(:)    ! data source index (ks)    
   integer, allocatable, public :: kx_ipw_ids(:)    ! see below          
   integer, allocatable, public :: kx_i_stride(:,:) ! i spacing and 1st index
   integer, allocatable, public :: kx_j_stride(:,:) ! j spacing and 1st index
   integer, allocatable, public :: kx_histogram(:,:,:,:,:)
   integer, allocatable, public :: kx_pt_count(:,:,:)
   integer, allocatable, public :: kx_obs_count(:,:,:,:)
!
! kx_idw_ids is derived from 3 character in kx_type. A value n means consider
! the integrated precipitable water (ipw) for the single nth range of
! integration; a value 1+kx_ipw_nlevs indicates use all ranges.                 
!
   real(rkind1), public :: kx_speed_min      ! min speed allowed
   real(rkind1), public :: kx_cloud_obscure  ! min prob obscured above allowed
!
! kx_locs: either lon of geost sat or min viewing lat of polar
! kx_dx: desired separation (km) between considered pts                        
! kx_params: probability function parameters to be determined
! kx_ipw_plevs: p-pairs designating ipw layers. These are min and max p in (mb)
!            for each successive ipw layer defined
! kx_ipw_bparam: pairs of values designating layer-dependent bins for
!            histograms and minimum q value to use for determining obs p-level
!            (See explanation of kx_ipw_... below)        
! kx_qcfac: factor to increase counts to compensate for QC rejections
   real(rkind1), allocatable, public :: kx_filters(:,:)
   real(rkind1), allocatable, public :: kx_locs(:) 
   real(rkind1), allocatable, public :: kx_dx(:)   
   real(rkind1), allocatable, public :: kx_qcfac(:) 
   real(rkind1), allocatable, public :: kx_params(:,:,:,:)
   real(rkind1), allocatable, public :: kx_ipw_plevs(:,:) 
   real(rkind1), allocatable, public :: kx_ipw_bparm(:,:,:) 
   real(rkind1), allocatable, public :: kx_akbk(:,:) 
   real(rkind1), allocatable, public :: kx_akbk_dlev(:,:) 
   real(rkind1), allocatable, public :: kx_range(:,:,:)  ! lat/lon range of obs
!
! kx_type:  2 vars for each k: (either G,P) (either V,I,B,W)
   character(len=10), allocatable, public :: kx_satname(:)
   character(len=2), allocatable, public :: kx_type(:) 
   character(len=25)  :: kx_varnames(30)
   character(len=120) :: kx_path_name   ! common path for kx_akbk and kx_satloc 
   character(len=120) :: kx_akbk_file_name   ! file for a,b(k) for vert interp
   character(len=120) :: kx_satloc_file_name ! file for counts by lon/hemis/time 
   character(len=14), public  :: kx_cdtime0_satloc  ! 1st ana t used to comp satloc
   character(len=14), public  :: kx_cdtime0_params  ! 1st ana t used to comp params
   character(len=240), public :: kx_akbk_file   ! path +  kx_akbk file name
   character(len=240), public :: kx_satloc_file ! path +  kx_satloc file name
!
! kx_ipw_... explanation:
!   These all refer to how layers of integated precipitable water (ipw) are
!       defined and used for various kx types and subtypes.
!   Each AMV type and subtype designed to simulate feature tracking based on
!       images of water vapor will consider one or more ipw layers. Each layer is
!       defined by a pair of real numbers indicating its top and bottom p in mb.     
!       The number of separate, perhaps overlapping, layers considered by any and
!       and all of the types and subtypes is set as kx_ipw_nlevs. Separate ipw
!       values will be computed for each defined layer.
!   Different AMV types and subtypes may use different sets of ipw values (i.e.,
!       different layers or combinations of layers). Also, they may use different
!       scalings of ipw values or minimum values of q for determining the vertical
!       location of observations of its considered layers. The number of distinct
!       layer combinations, scalings, and minimum q values is set as kx_ipw_itypes.
!       Each specified combination is refered to by an index 0,1,...,kx_ipw_itypes,
!       with 0 meaning none considered (i.e., for not a WV-based observation).
!       These latter indexes are specified in Table 2 for each kx,ks.
!   The specific combinations are defined in Table 3 and stored in the array
!       kx_ipw_bparms(1:2,1:kx_ipw_nlevs,1:kx_ipw_types). Its first index refers
!       to the two values q-scaling and q-min. These must be defined for each of
!       the kx_ipw_nlevs layers as ordered when kx_ipw_plevs was specified (they
!       must both be set to 0. for any layers not to be included in this
!       combination. There are kx_ipw_types of combinations that must be defined.
!   For example, for kx_ipw_nlevs=2, kx_ipw_plevs=  100. 700. 600. 900. means
!       that ipw will be determined for each of the 2 layers 100<p<700 and
!       600<p<900. Then for kx_ipw_types=4, the combinations
!       "1  3.00 .00010 0.  0.", "2 0. 0. 2.00 .0005", "3  3.00 .00020 1. .001",
!       means that for any kx,ks WV type refering to combination 1, only
!       ipw for layer 1 will be considered, with scaling of delta q by a factor 
!       of 3. and qmin set to .0001, to combination 2, only
!       ipw for layer 2 will be considered, with scaling of delta q by a factor
!       of 2. and qmin set to .0005 (gnerally the respecive lower and higher
!       values compared to combination 1 because the q values are larger
!       lower in the atmosphere), and combination 3 meaning that this kx,ks
!       type considers 2 different ipw layers, each with its own scaling and
!       q min). Note that although different kx types may refer to the same
!       layers the scaling and q min may need to differ, e.g., if one concerns
!       observations near the poles and the other near the equator, different
!       scaling and qmin m,ay be required because the q field has such different
!       magnitudes in the 2 regions.
!  As coded this appears the most general way of specifying which and how
!       layers are to be used for distinct kx,ks.
!
   contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
!
   subroutine kx_table_clean
!
!  De-allocate all kx arrays
!
   deallocate (kx_list,kx_satid,kx_ncnum,kx_type,kx_subtype)
   deallocate (kx_present,kx_i_stride,kx_j_stride,kx_ifunc)
   deallocate (kx_filters,kx_locs,kx_dx,kx_histogram,kx_nwcsm)
   deallocate (kx_params,kx_satname,kx_obs_count,kx_pt_count)
   deallocate (kx_ipw_plevs,kx_ipw_bparm,kx_ipw_ids,kx_ijgrid)
   deallocate (kx_qcfac,kx_akbk,kx_akbk_dlev,kx_range)
!
   end subroutine kx_table_clean
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine kx_table_read (filename,lprint,iread,ier)
!
! Read file with kx variables and tables:
!  table of obs types 
!  table of obs counts for all pbins given k, jbin
!  table of probability function parameters
!
   implicit none
   logical, intent(in) :: lprint
   integer, intent(in) :: iread
   integer, intent(out) :: ier
   character(len=*), intent(in) :: filename
!
   integer, parameter :: iunit=10
   integer :: j,k,n,jr,kr,nr
   integer :: nakbk,ios
   character(len=100) :: cdum
!
   ier=0
!
   open (iunit,file=trim(filename),status='old',iostat=ios)
   if (ios /= 0) then 
     ier=99
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open kx_table file for iunit=', &
                    iunit,' iostat=',ios,' and file name=',trim(filename)
     return
   elseif (lprint) then
     print *,' '
     print ('(2a)'),'kx_table_file opened: name=',trim(filename)
   endif
!
! Skip 2 lines and then read 1st part of table
   read (iunit,'(a1)') cdum
   read (iunit,'(a1)') cdum
!
! Read 1st part of table
   read (iunit,*) kx_varnames( 1),kx_num
   read (iunit,*) kx_varnames( 2),kx_cbins
   read (iunit,*) kx_varnames( 3),kx_pbins
   read (iunit,*) kx_varnames( 4),kx_jbins_lats
   read (iunit,*) kx_varnames( 5),kx_jbins_lors
   read (iunit,*) kx_varnames( 6),kx_jbins_time
   read (iunit,*) kx_varnames( 7),kx_nparams
   read (iunit,*) kx_varnames( 8),kx_nfilters
   read (iunit,*) kx_varnames( 9),kx_filt_flat
   read (iunit,*) kx_varnames(10),kx_filt_slev
   read (iunit,*) kx_varnames(11),kx_filt_torb
   read (iunit,*) kx_varnames(12),kx_filt_mcnt
   read (iunit,*) kx_varnames(13),kx_ipw_nlevs
   read (iunit,*) kx_varnames(14),kx_ipw_types
   read (iunit,*) kx_varnames(15),kx_times_satloc
   read (iunit,*) kx_varnames(16),kx_times_params
   read (iunit,*) kx_varnames(17),kx_field_kdim
   read (iunit,*) kx_varnames(18),kx_field_slots
   read (iunit,*) kx_varnames(19),kx_field_imax
   read (iunit,*) kx_varnames(20),kx_field_jmax
   read (iunit,*) kx_varnames(21),kx_speed_min
   read (iunit,*) kx_varnames(22),kx_cloud_obscure
   read (iunit,'(a18,1x,a)') kx_varnames(23),kx_cdtime0_satloc
   read (iunit,'(a18,1x,a)') kx_varnames(24),kx_cdtime0_params
   read (iunit,'(a13,1x,a)') kx_varnames(25),kx_path_name
   read (iunit,'(a20,1x,a)') kx_varnames(26),kx_satloc_file_name
   read (iunit,'(a18,1x,a)') kx_varnames(27),kx_akbk_file_name
!
   kx_akbk_file=trim(kx_path_name)//trim(kx_akbk_file_name) 
   kx_satloc_file=trim(kx_path_name)//trim(kx_satloc_file_name) 
   kx_jbins=kx_jbins_lats*kx_jbins_lors*kx_jbins_time
   nakbk=kx_field_kdim+1
!
   allocate (kx_present(kx_num))
   allocate (kx_list(kx_num))
   allocate (kx_satid(kx_num))
   allocate (kx_ncnum(kx_num))
   allocate (kx_type(kx_num))
   allocate (kx_subtype(kx_num))
   allocate (kx_i_stride(2,kx_num))
   allocate (kx_j_stride(2,kx_num))
   allocate (kx_filters(kx_nfilters,kx_num))
   allocate (kx_locs(kx_num))
   allocate (kx_dx(kx_num))
   allocate (kx_histogram(kx_cbins,kx_pbins,kx_jbins,kx_num,2))
   allocate (kx_params(kx_nparams,kx_pbins,kx_jbins,kx_num))
   allocate (kx_obs_count(kx_pbins,kx_jbins,kx_num,2))
   allocate (kx_pt_count(kx_jbins,kx_num,2))
   allocate (kx_satname(kx_num))
   allocate (kx_nwcsm(kx_num))
   allocate (kx_ifunc(kx_num))
   allocate (kx_ijgrid(kx_num))
   allocate (kx_qcfac(kx_num))
   allocate (kx_akbk(nakbk,2))
   allocate (kx_akbk_dlev(nakbk,2))
   allocate (kx_ipw_plevs(2,kx_ipw_nlevs))
   allocate (kx_ipw_bparm(2,kx_ipw_nlevs,kx_ipw_types))
   allocate (kx_ipw_ids(kx_num))
   allocate (kx_range(2,2,kx_num))
!
! Skip 3 lines and then read 2nd part of table
   read (iunit,'(a1)') cdum
   read (iunit,'(a1)') cdum
   read (iunit,'(a1)') cdum
!
! n_ kx_ said NC#___ satname__ type location dx_thin filters
   do k=1,kx_num
     read (iunit,'(i3,2i5,i7,2i3,1x,a10,a2,i3,2f6.1,2i3,f6.2,4f5.0)') &
           kr,kx_list(k),kx_satid(k),kx_ncnum(k),kx_subtype(k),       &
           kx_nwcsm(k),kx_satname(k),kx_type(k),kx_ipw_ids(k),        &
           kx_locs(k),kx_dx(k),kx_ifunc(k),kx_ijgrid(k),kx_qcfac(k),  &
           kx_filters(:,k)
   enddo
!
   if (iread > 2) then  ! Read the third portion of the table
!
! Skip 3 lines and then read 2nd part of table
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
!
     do k=1,kx_ipw_types
       read (iunit,*) kr,kx_ipw_bparm(:,:,k)
     enddo
     read (iunit,*) kx_varnames(28),kx_ipw_plevs(:,:)
!
   endif
!
!
   if (iread > 3) then  ! Read the fourth portion of the table
!
! Skip 3 lines
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
!
!  Read lat lon range table
     do k=1,kx_num
       read (iunit,'(i3,4f9.2)') kr,kx_range(:,:,k)
     enddo
!
   else
     do k=1,kx_num
       kx_range(1,1,k)= 999.
       kx_range(2,1,k)=-999.
       if (kx_type(k)(1:1) == 'P') then ! poler orbiter        
         kx_range(1,2,k)=-999.
         kx_range(2,2,k)= 999.
       else                             ! geostationary
         kx_range(1,2,k)= 999.
         kx_range(2,2,k)=-999.
       endif
     enddo
   endif
!
!
   if (iread > 4) then  ! Read the fifth portion of the table
!
! Skip 3 lines
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
!
!  Read obs count table
     do k=1,kx_num
       do j=1,kx_jbins
         read (iunit,'(2i3,10i8)') jr,kr,kx_obs_count(:,j,k,1)
       enddo
     enddo
!
   else
     kx_obs_count(:,:,:,1)=0
   endif
!
!
   if (iread > 5) then  ! Read the sixth portion of the table
!
! Skip 3 lines
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
     read (iunit,'(a1)') cdum
!
! Read parameters for probability functions                                                                      
     do k=1,kx_num
       do j=1,kx_jbins
         do n=1,kx_pbins
           read (iunit,'(3i3,10f10.0)') nr,jr,kr,kx_params(:,n,j,k)
         enddo
       enddo
     enddo
!
   else
     kx_params(:,:,:,:)=0.
   endif
!
   close (iunit)
!                                                                                     
! Initialize some arrays to 0
   kx_pt_count(:,:,:)=0
   kx_histogram(:,:,:,:,:)=0
!
   if (lprint) then
     print ('(a)'),'kx_table_file read'
   endif
!
! Read kx_akbk file if requested and compute values at data levels
   if (trim(kx_akbk_file_name) /= 'none') then
     call read_akbk (nakbk,kx_field_kdim,iunit,kx_akbk_file,lprint, &
                     kx_akbk,ios)
     if (ios /= 0) then
       ier=ier+abs(ios)
     else
       do k=1,kx_field_kdim                                                     
         kx_akbk_dlev(k,:)=0.5*(kx_akbk(k,:)+kx_akbk(k+1,:))               
       enddo                                                                 
       kx_akbk_dlev(nakbk,:)=kx_akbk(nakbk,:)
     endif 
   endif
!
   end subroutine kx_table_read
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine kx_table_write (filename)
!
! Write file of kx varaianles and tables
!
   implicit none
   character(len=*), intent(in) :: filename
!
   integer, parameter :: iunit=10
   integer :: j,k,n
   character(len=10) :: cdum
   character(len=120) :: ckid 
!
   open (iunit,file=trim(filename))
!
! Write 2 line spacers for table 1
   write (iunit,'(a)') '#'
   write (iunit,'(a)') '# Table 1: scalar values and ipw_plevs'
!
! Write table 1 values
   write (iunit,'(a,i4)') trim(kx_varnames( 1)),kx_num
   write (iunit,'(a,i4)') trim(kx_varnames( 2)),kx_cbins
   write (iunit,'(a,i4)') trim(kx_varnames( 3)),kx_pbins
   write (iunit,'(a,i4)') trim(kx_varnames( 4)),kx_jbins_lats
   write (iunit,'(a,i4)') trim(kx_varnames( 5)),kx_jbins_lors
   write (iunit,'(a,i4)') trim(kx_varnames( 6)),kx_jbins_time
   write (iunit,'(a,i4)') trim(kx_varnames( 7)),kx_nparams
   write (iunit,'(a,i4)') trim(kx_varnames( 8)),kx_nfilters
   write (iunit,'(a,i4)') trim(kx_varnames( 9)),kx_filt_flat
   write (iunit,'(a,i4)') trim(kx_varnames(10)),kx_filt_slev
   write (iunit,'(a,i4)') trim(kx_varnames(11)),kx_filt_torb
   write (iunit,'(a,i4)') trim(kx_varnames(12)),kx_filt_mcnt
   write (iunit,'(a,i4)') trim(kx_varnames(13)),kx_ipw_nlevs
   write (iunit,'(a,i4)') trim(kx_varnames(14)),kx_ipw_types
   write (iunit,'(a,i4)') trim(kx_varnames(15)),kx_times_satloc
   write (iunit,'(a,i4)') trim(kx_varnames(16)),kx_times_params
   write (iunit,'(a,i4)') trim(kx_varnames(17)),kx_field_kdim
   write (iunit,'(a,i4)') trim(kx_varnames(18)),kx_field_slots
   write (iunit,'(a,i5)') trim(kx_varnames(19)),kx_field_imax
   write (iunit,'(a,i5)') trim(kx_varnames(20)),kx_field_jmax
   write (iunit,'(a,f6.2)') trim(kx_varnames(21)),kx_speed_min
   write (iunit,'(a,f6.3)') trim(kx_varnames(22)),kx_cloud_obscure
   write (iunit,'(a18,1x,a)') trim(kx_varnames(23)),trim(kx_cdtime0_satloc)
   write (iunit,'(a18,1x,a)') trim(kx_varnames(24)),trim(kx_cdtime0_params)
   write (iunit,'(a13,1x,a)') trim(kx_varnames(25)),trim(kx_path_name)
   write (iunit,'(a20,1x,a)') trim(kx_varnames(26)),trim(kx_satloc_file_name)
   write (iunit,'(a18,1x,a)') trim(kx_varnames(27)),trim(kx_akbk_file_name)
!
! Write 3 line spacers for table 2
   write (iunit,'(a)') '#'
   write (iunit,'(a)') '# Table 2: define all AMV observation types'
   ckid='# n   kx said    NC# ks nw  satname typ it locat    dx fu'// &
        ' gd qsfac flat slev torb mcnt'
   write (iunit,'(a)') trim(ckid)
!
   do k=1,kx_num
     write (iunit,'(i3,2i5,i7,2i3,1x,a10,a2,i3,2f6.1,2i3,f6.2,4f5.0)') &
           k,kx_list(k),kx_satid(k),kx_ncnum(k),kx_subtype(k),         &
           kx_nwcsm(k),kx_satname(k),kx_type(k),kx_ipw_ids(k),         &
           kx_locs(k),kx_dx(k),kx_ifunc(k),kx_ijgrid(k),kx_qcfac(k),   &
           kx_filters(:,k)
   enddo
!
! Write 3 line spacers for table 3
   write (iunit,'(a)') '#'
   write (iunit,'(a)') '# Table 3: define ipw_bparams'
   ckid='# q scale1    qmin1 scale2    qmin2 ...'
   write (iunit,'(a)') trim(ckid)
!
! Write ipw_bparm table of q scales and qmin values
   do k=1,kx_ipw_types
     write (iunit,'(i3,5(f7.2,f9.5))') k,kx_ipw_bparm(:,:,k)
   enddo
   write (iunit,'(a,10f8.0)') trim(kx_varnames(28)),kx_ipw_plevs(:,:)
!
! Write 3 line spacer for table 4 (lon/lat ranges)
   write (iunit,'(a)') '#'
   write (iunit,'(a)') '# Table 4 lat/lon ranges'
   write (iunit,'(a)') '# k     lonW     lonE     lat1     lat2'
!
   do k=1,kx_num
     write (iunit,'(i3,4f9.2)') k,kx_range(:,:,k)
   enddo
!
! Write 3 line spacers for table 5
   write (iunit,'(a)') '#'
   write (iunit,'(2a)') '# Table 5: define target or simulation counts', &
                        ' (inflated if param table present and not all 0.)'
   write (iunit,'(a)') '# obs count for indicated jbin, k for all plevs'
!
! Write obs count
   do k=1,kx_num
     do j=1,kx_jbins
       write (iunit,'(2i3,10i8)') j,k,kx_obs_count(:,j,k,2)
     enddo
   enddo
!
! Write 3 line spacers for table 6
   write (iunit,'(a)') '#'
   write (iunit,'(a)') '# Table 6: table of probability function parameters'
   write (iunit,'(a)') '# params  for indicated  plevs, jbin, k indexes'
!
! Write parameters
   do k=1,kx_num
     do j=1,kx_jbins
       do n=1,kx_pbins
         write (iunit,'(3i3,f10.5,9f10.3)') n,j,k,kx_params(:,n,j,k)
       enddo
     enddo
   enddo
!
   close (iunit)
   print *,' '
   print *,'kx_table_file written: name=',trim(filename)
!                                                                                     
   end subroutine kx_table_write
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine kx_table_sum (cfilein,cfileout)
!
! Sum counts in kx_table (sum over all p for each j, and also 
! compute the total; separate sums for each kx,ks)
!
   implicit none
!
   character(len=*), intent(in) :: cfilein
   character(len=*), intent(in) :: cfileout
!
   integer, parameter :: iunit=10
   integer :: k,j,n
   integer :: npsum(kx_pbins+1,kx_num)   
!
   npsum(:,:)=0
   do k=1,kx_num
     do n=1,kx_pbins
       do j=1,kx_jbins
         npsum(n,k)=npsum(n,k)+kx_obs_count(n,j,k,1)
       enddo
       npsum(kx_pbins+1,k)=npsum(kx_pbins+1,k)+npsum(n,k)
     enddo
   enddo
!
   open (iunit,file=trim(cfileout))
   write (iunit,'(a)') 'Summed average counts for each p-level band'
   write (iunit,'(2a)') 'file in = ',trim(cfilein)
   write (iunit,'(2x,20i6)') kx_list(:)
   write (iunit,'(2x,20i6)') kx_satid(:)
   write (iunit,'(a)') ' '
!
   do n=1,kx_pbins
     write (iunit,'(i2,20i6)') n,npsum(n,:)
   enddo
!
   write (iunit,'(a)') ' '
   write (iunit,'(2x,20i6)') npsum(kx_pbins+1,:)
   close (iunit)
!
   end subroutine kx_table_sum
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine kx_histogram_print (cfileout)
!
! Print the computed histograms of cf or scaled delta ipw for
! diagnostic purposes
!
   implicit none
!
   character(len=*), intent(in) :: cfileout
!
   integer, parameter :: fileformat=1  ! indicator for file format
   integer, parameter :: iunit=10
   integer :: i,j,k,n
   integer :: hsum
   real(rkind1) :: percnt(kx_cbins)
   real(rkind1) :: xmlog10(kx_cbins)
!
   open (iunit,file=trim(cfileout))
   write (iunit,'(a,i5)') 'fileformat=',fileformat
   write (iunit,'(a,4i5)') 'ncbins,npbins,njbins,kxnum=', &
                       kx_cbins,kx_pbins,kx_jbins,kx_num
   write (iunit,'(a)') ' k  j  n   kx  sid  histsum'
   do k=1,kx_num
     do j=1,kx_jbins
       do n=1,kx_pbins
         hsum=sum(kx_histogram(:,n,j,k,2))
         do i=1,kx_cbins
           percnt(i)=kx_histogram(i,n,j,k,2)/real(max(hsum,1))
           xmlog10(i)=-log10(max(percnt(i),1.e-9))
         enddo
         write (iunit,'(i2,2i3,2i5,i9)') k,j,n,kx_list(k),kx_satid(k),hsum
         write (iunit,'(3x,i9,20i6)') kx_histogram(:,n,j,k,2)
         write (iunit,'(5x,f7.2,20f6.3)') percnt(:)
         write (iunit,'(5x,f7.2,20f6.3)') xmlog10(:)
       enddo
     enddo
   enddo
   close (iunit)
   print *,'Histogram printed to file :',trim(cfileout)
!
   end subroutine kx_histogram_print
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine kx_latlon_range (k,lon,lat)
!
!  Determine lat and lon ranges for each kx, ks type
!  (See m_amv_view:amv_latlon_range for use)
!
   implicit none
!
   integer, intent(in) :: k
   real(rkind1), intent(in) :: lon,lat
!
   real(rkind1) :: x
!
   if (kx_type(k)(1:1) == 'P') then       ! polar orbiter
     if (lat < 0.) then 
       kx_range(1,2,k)=max(kx_range(1,2,k),lat) ! most northern lat in SH
     else
       kx_range(2,2,k)=min(kx_range(2,2,k),lat) ! most southern lat in NH
     endif
!
   else                                   ! geostationary
     kx_range(1,2,k)=min(kx_range(1,2,k),lat)  !  min latitude    
     kx_range(2,2,k)=max(kx_range(2,2,k),lat)  !  max latitude   
!
! Determine max longitudinal separations from nadir
! (first adjust for case when separation crosses the prime meridian)
     x=lon-kx_locs(k)
     if (x < -180.) then
       x=x+360.
     elseif (x > 180.) then
       x=x-360.
     endif
!
     if (x < 0.) then  ! max separation westward from nadir  
       kx_range(1,1,k)=min(kx_range(1,1,k),x)
     else              ! max separation eastward from nadir  
       kx_range(2,1,k)=max(kx_range(2,1,k),x)
     endif
!
   endif ! test on plor or geostationary
!
   end subroutine kx_latlon_range
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_kx_table

