   module m_saved_data
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
  use m_kinds, only : rkind1, rkind2
!
   private
   public :: sd_setup
   public :: sd_set_subtype 
   public :: sd_copy
   public :: sd_diff
   public :: sd_clean
!
   integer, public, parameter :: obs_types_max=200
   integer, public :: obs_ids_dim1
   integer, public :: obs_ids_dim2
   integer, public :: obs_sing_dim2
   integer, public :: obs_multi_dim1
   integer, public :: obs_multi_dim2
   integer, public :: obs_multi_dim3
   integer, public :: nobs_sing_2(2),nobs_multi_2(2),nobs_all_2(2) 
   integer, public :: count_types_num
   integer, public :: count_types(obs_types_max,2)
   real(rkind1), allocatable, public :: obs_ids(:,:,:)
   real(rkind1), allocatable, public :: obs_sing(:,:,:)
   real(rkind1), allocatable, public :: obs_multi(:,:,:,:)
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sd_clean
   deallocate (obs_ids,obs_sing,obs_multi)
   end subroutine sd_clean 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sd_setup (dtype)
!
   implicit none
!
   character(len=*), intent(in) :: dtype
!
   obs_ids_dim1=8
   if (trim(dtype) == 'PREPBUFR') then
     obs_multi_dim1=255
     obs_multi_dim2=3
     obs_multi_dim3=3000
     obs_sing_dim2=1000000
     obs_ids_dim2=obs_sing_dim2+obs_multi_dim3
   else
     if (trim(dtype) == 'GPSRO') then
       obs_multi_dim1=350     
     elseif (trim(dtype) == 'IASI') then
       obs_multi_dim1=620     
     elseif (trim(dtype) == 'AIRS') then
       obs_multi_dim1=300     
     else
       obs_multi_dim1=20
     endif
     obs_sing_dim2=1    ! array obs_sing not used except for PREPBUFR
     obs_multi_dim2=2   
     obs_ids_dim2=int(8.0e6/real(obs_multi_dim1))
     obs_multi_dim3=obs_ids_dim2
   endif
! 
   allocate (obs_ids(obs_ids_dim1,obs_ids_dim2,2))
   allocate (obs_sing(obs_multi_dim2,obs_sing_dim2,2))
   allocate (obs_multi(obs_multi_dim1,obs_multi_dim2,obs_multi_dim3,2))
!
   nobs_all_2(:)=0
   nobs_sing_2(:)=0
   nobs_multi_2(:)=0
!
   end subroutine sd_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sd_copy (obs_max_chan_or_levs,obs_max_fields,dtype,i_file, &
                       itype,olevs,nobs,lat,lon,obs_ps,obs_data,obs_levs)
!
! Copy data into arrays for saving
! 
   implicit none
!
   integer :: i_file
   integer :: itype
   integer :: olevs, nobs
   integer :: obs_max_chan_or_levs,obs_max_fields
   real(rkind2) :: lat,lon
   real(rkind2) :: obs_ps
   real(rkind2) :: obs_levs(obs_max_chan_or_levs)
   real(rkind2) :: obs_data(obs_max_chan_or_levs,obs_max_fields)
   character(len=*), intent(in) :: dtype
!
   logical :: obs_ok
   integer :: nobs_sing,nobs_multi,nobs_all 
   integer :: k
!
   nobs_all=nobs_all_2(i_file)
   nobs_sing=nobs_sing_2(i_file)
   nobs_multi=nobs_multi_2(i_file)
!
   obs_ok=.true.
   if (nobs > obs_ids_dim2 .or. olevs < 1) then
     obs_ok=.false.
   endif
   if (trim(dtype) == 'PREPBUFR') then
     if ( (olevs == 1 .and. nobs_sing  == obs_sing_dim2) .or. &
          (olevs > 1  .and. nobs_multi == obs_multi_dim3) ) then
       obs_ok=.false.
     endif
   else
     if (nobs_multi == obs_multi_dim3) then
       obs_ok=.false.
     endif
   endif
!
   if (obs_ok) then
     nobs_all=nobs_all+1
     obs_ids(1,nobs_all,i_file)=itype
     obs_ids(2,nobs_all,i_file)=lat
     obs_ids(3,nobs_all,i_file)=lon
     obs_ids(4,nobs_all,i_file)=0.   ! no check on time will be performed
     obs_ids(5,nobs_all,i_file)=olevs
     obs_ids(6,nobs_all,i_file)=0.      ! default value 
     obs_ids(obs_ids_dim1,nobs_all,i_file)=obs_ps    
!
     if (trim(dtype) == 'PREPBUFR' .and. olevs ==1) then
       nobs_sing=min(nobs_sing+1,obs_sing_dim2)
       obs_ids(6,nobs_all,i_file)=nobs_sing
       obs_sing(1,nobs_sing,i_file)=obs_levs(1) 
       obs_sing(2:obs_multi_dim2,nobs_sing,i_file)= &
                           obs_data(1,1:obs_multi_dim2-1)
     else
       nobs_multi=min(nobs_multi+1,obs_multi_dim3)
       obs_ids(6,nobs_all,i_file)=nobs_multi
       do k=1,olevs
         obs_multi(k,1,nobs_multi,i_file)=obs_levs(k) 
         obs_multi(k,2:obs_multi_dim2,nobs_multi,i_file)= &
                           obs_data(k,1:obs_multi_dim2-1)
       enddo
     endif
!
   endif
!
   nobs_all_2(i_file)=nobs_all
   nobs_sing_2(i_file)=nobs_sing
   nobs_multi_2(i_file)=nobs_multi
!
   end subroutine sd_copy 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sd_diff (dtype)
!
! Difference data
! 
   implicit none
!
   character(len=*), intent(in) :: dtype
!
   logical:: lsame
   integer :: i,n,olevs,nid
   integer :: nobs_sing,nobs_multi,nobs_all 
!
   if (nobs_all_2(1) /= nobs_all_2(2)) then
     print *,'NOBS_ALL NOT SAME ON 2 FILES:', nobs_all_2
     stop
   endif
!
   if (trim(dtype) == 'PREPBUFR' .and. nobs_sing_2(1) /= nobs_sing_2(2)) then
     print *,'NOBS_SING NOT SAME ON 2 FILES:', nobs_sing_2
     stop
   endif
!
   if (nobs_multi_2(1) /= nobs_multi_2(2)) then
     print *,'NOBS_MULTI NOT SAME ON 2 FILES:', nobs_multi_2
     stop
   endif
!
   nobs_all=nobs_all_2(1)
   nobs_sing=nobs_sing_2(1)
   nobs_multi=nobs_multi_2(1)
!
   do n=1,nobs_all
!
     lsame=.true.
     do i=1,obs_ids_dim1-1  ! -1 because ps need not be same
       if (obs_ids(i,n,1) /= obs_ids(i,n,2)) then
         lsame=.false.
       endif
     enddo
!
     if (.not. lsame) then
       print *,'HEADERS DIFFER for n=',n
       do i=1,obs_ids_dim1
         print *,i,obs_ids(i,n,1),obs_ids(i,n,2)
       enddo
       exit
     endif
!
     olevs=nint(obs_ids(5,n,1))   
     nid=nint(obs_ids(6,n,1))
     obs_ids(8,n,2)=obs_ids(8,n,2)-obs_ids(8,n,1)
     if (trim(dtype) == 'PREPBUFR' .and. olevs == 1) then
       obs_sing(2:obs_multi_dim2,nid,2)=obs_sing(2:obs_multi_dim2,nid,2)- &
                                        obs_sing(2:obs_multi_dim2,nid,1)
     else
       obs_multi(1:olevs,2:obs_multi_dim2,nid,2)=                         &
                               obs_multi(1:olevs,2:obs_multi_dim2,nid,2)- &
                               obs_multi(1:olevs,2:obs_multi_dim2,nid,1)
     endif
!
   enddo
!
   end subroutine sd_diff 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine sd_set_subtype 
!
   implicit none
!
   integer :: n,n1,k,itype
!
   count_types(:,:)=0
!
   do n=1,nobs_all_2(1)
     itype=int(obs_ids(1,n,1))
     n1=0                                 ! flag (0 means type not found yet)
     do k=1,obs_types_max                 ! search through list
       if (n1 == 0) then                  ! requested type not found in list yet
         if (count_types(k,1) == 0) then  ! fill type since slot is empty 
           count_types(k,1)=itype
           count_types_num=k              ! number of found types
         endif 
         if (itype == count_types(k,1)) then ! check if type found in list 
           n1=k                              ! type found in list
         endif
       endif
     enddo
     obs_ids(7,n,1)=n1
     count_types(n1,2)=count_types(n1,2)+1
   enddo
!
   end subroutine sd_set_subtype 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_saved_data
