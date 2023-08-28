   module m_rad_prob
!
!  Determine effective radiative surface for radiance observations by using 
!  a probability function based on some NR field value. Also set field names 
!  used to compute effects of hydrometeorsor aerosols.
!
   use m_kinds, only : rkind1
   use m_set_unit_nums, only  : un_info
!
   implicit none
!
   private
!
! public routines
   public :: rad_prob_setup
   public :: rad_prob_compute
!
! private routines
   private :: rad_prob_read_rc
!
! public variables
   integer, parameter, public :: prob_flds_max=10
   integer, parameter, public :: list_cloud_maxn=10
   integer, parameter, public :: list_aerosol_maxn=30
   integer, public :: prob_seed         ! to be added to random seed 
   integer, public :: prob_vars_nums
   integer, public :: prob_vars_flds
   integer, public :: list_cloud_nums   ! number of cloud types for RTM to use
   integer, public :: list_aerosol_nums ! number of aerosol types for RTM to use
   real(rkind1), public :: prob_vars(10,prob_flds_max)
   character(len=16), public :: prob_names(prob_flds_max) 
   character(len=8),  public :: list_cloud_names(list_cloud_maxn)   
   character(len=8),  public :: list_aerosol_names(list_aerosol_maxn) 
! prob_names lists fields used to determine prob of rad sfc being elevated
! list_cloud_names are NR fields used to input clouds into RTM
! list_aerosol_names are NR fields used to input aerosols into RTM
!
   integer :: prob_seed0      ! part of random seed common to all instruments
   integer :: prob_seed1      ! part of random seed for particular instrument
   integer :: prob_k_iter
   integer :: prob_nps
   integer :: prob_ids(prob_flds_max)
!
   character(len=*), parameter :: myname='m_rad_prob'
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rad_prob_setup (dtype,rcfile,lprint,ier)
!
!  Set parameters that will be used for determination of cloud, precip, 
!  and aerosol aeffects on radiances 
!
   use m_read_profiles, only : prof_num_2d, prof_names, prof_kmax
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c  
   use m_obs_list, only : obs_list_names
!
   implicit none   
!    
   logical, intent(in) :: lprint
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: rcfile 
   integer, intent(out) :: ier
!   
   logical, parameter :: lstop=.false.
   integer :: ier1
   integer :: fid,m1,m2,n,n1,n2
   character(len=*),parameter :: mysub=myname//'::rad_thin_flds_setup'
!
   call rad_prob_read_rc (dtype,lprint,rcfile,ier1)
   ier=ier1
   if (ier /= 0) return
!   
   m1=obs_list_len_i+1
   m2=obs_list_len_i+obs_list_len_r
   n1=obs_list_len_r
   n2=prof_num_2d
   call find_name (n2,prof_names(1:n2,2),lstop,mysub,'PS',prob_nps)
!
   do n=1,prob_vars_flds      ! first check in prof header info
     call find_name (n1,obs_list_names(m1:m2),lstop,' ', &
                     trim(prob_names(n)),fid)
     if (fid == 0) then       ! check in 2d profile fields
       call find_name (n2,prof_names(1:n2,2),lstop,mysub, &
                     trim(prob_names(n)),fid)
       fid=-fid   ! flag with - sign to indicate look in prof 2d fields
     endif
     if (fid == 0) then
       ier=ier+1 
       if (lprint) then
         print *,'Error in ',mysub
         print *,'Requested data type ',trim(prob_names(n)),' not found in lists'
       endif
     else
       prob_ids(n)=fid
     endif
   enddo    
!
   prob_k_iter=int(log(real(prof_kmax+1))/log(2.)) + 1
   prob_seed=prob_seed0+prob_seed1
!
   if (lprint) then
     call rad_prob_print 
   endif
!
   end subroutine rad_prob_setup
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad_prob_read_rc (dtype,lprint,rcfile,ier)
!
! Read resource file for instructions on computing radiances that depend on
! instrument type (particularly regarding cloud and precip effects)
!
   implicit none
! arguments
   logical, intent(in)  :: lprint
   integer, intent(out) :: ier
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: rcfile
!
! local variables
   logical, parameter :: lstop=.false. ! must be false when mpi is used
   logical :: found
   integer :: ios
   integer :: prob_format_header
   integer :: prob_format_recs
   integer, parameter :: iunit=un_info
   integer :: m,n,n1,n2,nrecs,nf
   character(len=30)  :: read_name
   character(len=260) :: read_info(100)
   character(len=30)  :: cdum
   character(len=*), parameter :: mysub=myname//'::rad_prob_read_rc'
   character(len=120) :: mysub0
!
   ier=0  ! error counter
!
   open (iunit,file=trim(rcfile),form='formatted',status='old',iostat=ios)
   if (ios /= 0) then 
     ier=99
     print *,' '
     print ('(a,i3,a,i4,2a)'),' ERROR attempting to open rad_prob rc file for iunit=', &
                    iunit,' iostat=',ios,' and file name=',trim(rcfile)
     return  
   elseif (lprint) then
     print *,' '
     print ('(3a,i4)'),' Rc file=',trim(rcfile),' opened on unit=',iunit
   endif
!  
! Read format indicators
   read (iunit,*) cdum,prob_format_header,prob_format_recs
!
! Read common random seed
   read (iunit,*) cdum,prob_seed0
!
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
     ier=1
     if (lprint) then 
       print *,'info for dtype=',trim(dtype),' not found'
     endif
   else
     ier=0
!
! Set msub0 to instruct error handeling in sub=find_name_2 
   if (lprint) then  
     mysub0=mysub
   else
     mysub0(1:1)=' '
   endif
!
! Look for values of list_cloud_nums and list_cloud_names
! (These specify cloud types used by the RTM)
   call find_name_2 (nrecs,read_info,lstop,mysub0,'list_cloud_nums',n)
   if (n == 0) then
     list_cloud_nums=0
     if (lprint) then
       print ('(a)'),'DEFAULT VALUE USED FOR list_cloud_nums'
     endif
   else 
     read (read_info(n),*) cdum,list_cloud_nums
     if (list_cloud_nums > 0) then
       call find_name_2 (nrecs,read_info,lstop,mysub0,'list_cloud_names',n)
       if (n == 0 .or. list_cloud_nums > list_cloud_maxn) then
         ier=ier+1
         if (list_cloud_nums > list_cloud_maxn) then
           print ('(2a,i2)'),'REQUSTED NUMBER OF OF CLOUD FIELDS= ',    &
                  list_cloud_nums,' GREATER THAN MAX NUMBER ALLOWED= ', &
                  list_cloud_maxn
         endif
       else 
         read (read_info(n),*) cdum,list_cloud_names(1:list_cloud_nums)
       endif
     endif
   endif
!
! Look for values of list_aerosol_nums and list_aerosol_names
! (These specify aerosol types used by the RTM)
   call find_name_2 (nrecs,read_info,lstop,mysub0,'list_aerosol_nums',n)
   if (n == 0) then
     list_aerosol_nums=0
     if (lprint) then
       print ('(a)'),'DEFAULT VALUE USED FOR list_aerosol_nums'
     endif
   else 
     read (read_info(n),*) cdum,list_aerosol_nums
     if (list_aerosol_nums > 0) then 
       call find_name_2 (nrecs,read_info,lstop,mysub0,'list_aerosol_names',n)
       if (n == 0 .or. list_aerosol_nums > list_aerosol_maxn) then
         ier=ier+1
         if (list_aerosol_nums > list_aerosol_maxn) then
           print ('(2a,i2)'),'REQUSTED NUMBER OF OF AEROSOL FIELDS= ',    &
                  list_aerosol_nums,' GREATER THAN MAX NUMBER ALLOWED= ', &
                  list_aerosol_maxn
         endif
       else 
         do n1=1,list_aerosol_nums,8
           n2=min(n1+7,list_aerosol_nums)
           read (read_info(n),*) cdum,list_aerosol_names(n1:n2)
           n=n+1
         enddo
       endif
     endif
   endif
!
! Look for value of random seed
   call find_name_2 (nrecs,read_info,lstop,mysub0,'random_seed1',n)
   if (n == 0) then
     ier=ier+1
   else 
     read (read_info(n),*) cdum,prob_seed1
   endif
!
! Look for value of prob_vars_nums and prob_vars_flds
   call find_name_2 (nrecs,read_info,lstop,mysub0,'prob_num_vars&fields',n)
   if (n == 0) then
     prob_vars_nums=3
     prob_vars_flds=0
     if (lprint) then
       print *,'DEFAULT VALUE USED FOR prob_vars_nums and prob_vars_flds'
     endif
   else 
     read (read_info(n),*) cdum,prob_vars_nums,prob_vars_flds
     if (prob_vars_flds > 0) then
       do m=1,prob_vars_flds
         n=n+1
         read (read_info(n),'(a16,10f6.3)') prob_names(m), &
                                          prob_vars(1:prob_vars_nums,m)
       enddo
     endif
   endif    
!
   endif ! check on found
!
   end subroutine rad_prob_read_rc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rad_prob_print 
!
!  Print some information read from the rad_prob resource file
!
   implicit none
!    
   integer :: n
!
   print *,' '
   print ('(a)'),'Setup of m_rad_prob:'
   print ('(a,i2)'),'Number of cloud fields to consider = ',list_cloud_nums
   if (list_cloud_nums > 0) then
     print ('(a,10a10)'),'cloud fields: ',list_cloud_names(1:list_cloud_nums)
   endif
   print ('(a,i2)'),'Number of aerosol fields to consider = ',list_aerosol_nums
   if (list_aerosol_nums > 0) then
     print ('(a,10a10)'),'aerosol fields: ', &
           list_aerosol_names(1:list_aerosol_nums)
   endif
   print ('(a,2i8)'),'random seeds=',prob_seed0,prob_seed1
   print ('(2a,i2)'),'Number of fields used to determine elevated ', &
           'radiative surface = ',prob_vars_flds
   if (prob_vars_flds > 0) then 
     print ('(a)'),' n  variable_name     id       probability specifiers'  
     do n=1,prob_vars_flds
       print ('(i2,2x,a16,i4,10f7.4)'),n,prob_names(n),prob_ids(n), &
                                       prob_vars(1:prob_vars_nums,n)
     enddo
   endif
!
   end subroutine rad_prob_print
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rad_prob_compute (nob,nlevs,plevs,prob_index, &
                                dtype,effective_ks,effective_ps,obs_quality)
!
!  Compute effective radiative surface and an obs quality indicator 
!  based on probability of radiances being obstructed by clouds or other fields
!
   use m_obs_list, only : obs_list_r
   use m_read_profiles, only: prof_all      

   implicit none
!    
   integer, intent(in) :: nob
   integer, intent(in) :: nlevs
   real(rkind1),intent(in) :: plevs(nlevs)
   character(len=*), intent(in) :: dtype
!
   integer, intent(out) :: effective_ks
   integer, intent(out) :: prob_index
   real(rkind1),intent(out) :: effective_ps
   real(rkind1),intent(out) :: obs_quality  ! quality flag, esp. for IASI QC
!
   integer :: i,j,k_below
   real(rkind1) :: sigma
   real(rkind1) :: val
   real(rkind1) :: prob
   real(rkind1) :: x
   character(len=*),parameter :: mysub=myname//"::rad_prob_compute"
!
   sigma=1.
!
   do i=1,prob_vars_flds
!
!  determine value of id specified field from which to determine probabilities 
     if (prob_ids(i) > 0) then 
       val=obs_list_r(prob_ids(i),nob)
     else
       val=prof_all(-prob_ids(i),nob)
     endif 
!
!  determine probability of radiance being affected given value=val
     if (val >= prob_vars(2,i)) then
       prob=1.
     elseif (val <= prob_vars(1,i)) then
       prob=0.
     else
       prob=(val-prob_vars(1,i))/(prob_vars(2,i)-prob_vars(1,i))
     endif
!
     call random_number(x)  
     if (x < prob) then
       sigma=prob_vars(prob_vars_nums,i)
     endif
     if (sigma /= 1.) then 
       prob_index=i+1
       exit
     endif
   enddo
! 
   if (sigma == 1.) then
     prob_index=1
   endif     
   effective_ps=sigma*prof_all(prob_nps,nob)
!
! Find the grid point index for the pressure level just above pcldtop by
! successively dividing the range it is in by factors of 2.
   effective_ks=1
   k_below=nlevs
   if (effective_ps >= plevs(nlevs)) then
     effective_ks=nlevs
   else
     do i=1,prob_k_iter
       j=(effective_ks+k_below)/2 
       if (effective_ps >= plevs(j)) then
         effective_ks=j
       else
         k_below=j
       endif
     enddo
   endif     ! test if effective top between surface and lowest level
!
! Set variable that designates quality of observations based on the 
! criteria determined above. The default formulation is to set obs_quality  
! to the percent of the scene that is clear, set here to the value of 
! sigma determined randomly above. For SSMI, the value is set for the 
! bufr header value RFLAG (1=rain, 0=norain).  
!
   if (trim(dtype) == 'SSMIS' .or. trim(dtype) == 'AVCSAM' .or. &
       trim(dtype) == 'AVCSPM') then
     if (sigma < 0.9) then
       obs_quality=1.
     else    
       obs_quality=0.
     endif
   elseif (trim(dtype) == 'CRIS' .or. trim(dtype) == 'CRISFSR' ) then 
     obs_quality=0.
   else
     obs_quality=min(max(100.*sigma,0.),100.)
   endif
!
   end subroutine rad_prob_compute
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_rad_prob
