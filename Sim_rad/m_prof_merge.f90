   module m_prof_merge 
!
! Module containing routines and arrays used to merge 2 observation
! location profile files containing different sets of fields
!
!  Initial code: Ronald Errico  October 15 2014 
!
   use m_kinds, only : rkind1
!
   implicit none
!
   private
   public :: merge_save_header
   public :: merge_compare_headers
   public :: merge_prof_name_list 
   public :: merge_allocate_prof_arrays
   public :: merge_profile_recs 
!
   integer, public :: merge_num_2d
   integer, public :: merge_num_3d
   character(len=12), public :: file1_list_names(200)
   character(len=12), public :: file1_prof_names(100,2:3)
   character(len=12), public :: merge_prof_names(100,2:3)
!
   logical :: merge_copy(100,2:3)
!
   integer :: file1_num_2d, file1_num_3d
   integer :: file2_num_2d, file2_num_3d
   integer :: file1_list_nsum
   integer :: file1_list_len_i, file1_list_len_r, file1_list_len_c
   integer :: file1_kmax, file1_list_format_recs, file_list_nsum
   integer :: file1_list_tslots, file1_list_subtypes
   integer :: file1_list_counter(100,50)
   integer, allocatable :: merge_list_i(:,:)
!
   real(rkind1), allocatable :: merge_list_r(:,:)
   real(rkind1), allocatable :: merge_file1_2d(:)
   real(rkind1), allocatable :: merge_file2_2d(:)
   real(rkind1), allocatable :: merge_file3_2d(:)
   real(rkind1), allocatable :: merge_file1_3d(:,:)
   real(rkind1), allocatable :: merge_file2_3d(:,:)
   real(rkind1), allocatable :: merge_file3_3d(:,:)
!
   character(len=16), allocatable :: merge_list_c(:,:)
!
   contains
!
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine merge_profile_recs (unit1,unit2,unit3,n1,n2,n3,ctest)
!
!  Merge 2 sets of profiles 
!
   implicit none
   integer, intent(in) :: unit1,unit2,unit3
   integer, intent(in) :: n1,n2,n3
   character(len=*), intent(in) :: ctest
!
   logical :: agree
   integer :: n,m
!
   read (unit1) merge_list_i(:,1),merge_list_r(:,1),merge_list_c(:,1)
   read (unit1) merge_file1_2d(1:file1_num_2d), &
                merge_file1_3d(:,1:file1_num_3d)
!
   read (unit2) merge_list_i(:,2),merge_list_r(:,2),merge_list_c(:,2)
   read (unit2) merge_file2_2d(1:file2_num_2d), &
                merge_file2_3d(:,1:file2_num_3d)
!
! Check obs info list
   agree=.true.
   do n=1,file1_list_len_i
     if (merge_list_i(n,1) /= merge_list_i(n,2)) agree=.false.
   enddo
   do n=1,file1_list_len_r
     if (merge_list_r(n,1) /= merge_list_r(n,2)) agree=.false.
   enddo
   do n=1,file1_list_len_c
     if (merge_list_c(n,1) /= merge_list_c(n,2)) agree=.false.
   enddo
!
   if (.not. agree) then
     print *,'Obs info do not agree for n1,n2,n3=',n1,n2,n3
     do n=1,file1_list_len_i
       print ('(i3,2i20)'),n,merge_list_i(n,1:2)
     enddo
     do n=1,file1_list_len_r
       print ('(i3,1p2e20.8)'),n,merge_list_r(n,1:2)
     enddo
     do n=1,file1_list_len_c
       print ('(i3,2a20)'),n,merge_list_c(n,1:2)
     enddo
     stop
   endif
!
! Merge profiles
   m=file1_num_2d
   merge_file3_2d(1:m)=merge_file1_2d(1:m)
   if (merge_num_2d > m) then
     do n=1,file2_num_2d 
       if (merge_copy(n,2)) then
         m=m+1 
         merge_file3_2d(m)=merge_file2_2d(n)
       endif
     enddo
   endif
!
   m=file1_num_3d
   merge_file3_3d(:,1:m)=merge_file1_3d(:,1:m)
   if (merge_num_3d > m) then
     do n=1,file2_num_3d 
       if (merge_copy(n,3)) then
         m=m+1 
         merge_file3_3d(:,m)=merge_file2_3d(:,n)
       endif
     enddo
   endif
!
! Write merged records
   write (unit3) merge_list_i(:,1),merge_list_r(:,1),merge_list_c(:,1)
   write (unit3) merge_file3_2d(1:merge_num_2d), &
                 merge_file3_3d(:,1:merge_num_3d)
!
! Printing for testing (only first obs in each tslot)
   if (ctest == 'T' .and. n3 == 1) then
     print *,'Print to test merge for n1,n2=',n1,n2 
     print *,'2d fields in file 1 ',file1_num_2d,merge_file1_2d(1:file1_num_2d)
     print *,'2d fields in file 2 ',file2_num_2d,merge_file2_2d(1:file2_num_2d)
     print *,'2d merged fields ',merge_num_2d,merge_file3_2d(1:merge_num_2d)
     do n=1,file1_num_3d
       print *,'3d fields in file 1 ',n,merge_file1_3d(:,n)
     enddo
     do n=1,file2_num_3d
       print *,'3d fields in file 2 ',n,merge_file2_3d(:,n)
     enddo
     do n=1,merge_num_3d
       print *,'3d merged fields ',n,merge_file3_3d(:,n)
     enddo
   endif
!
   end subroutine merge_profile_recs 
!
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine merge_allocate_prof_arrays 
!
!  Allocate arrays required to hold profiles for 2 input files and the 
!  merged set.
!
   use m_read_profiles, only : prof_num_2d, prof_num_3d, prof_kmax
!
   implicit none
!
   allocate (merge_list_i(file1_list_len_i,2))
   allocate (merge_list_r(file1_list_len_r,2))
   allocate (merge_list_c(file1_list_len_c,2))
   allocate (merge_file1_2d(file1_num_2d+1))
   allocate (merge_file2_2d(file2_num_2d+1))
   allocate (merge_file3_2d(merge_num_2d+1))
   allocate (merge_file1_3d(prof_kmax,file1_num_3d+1))
   allocate (merge_file2_3d(prof_kmax,file2_num_3d+1))
   allocate (merge_file3_3d(prof_kmax,merge_num_3d+1))
!
   end subroutine merge_allocate_prof_arrays 
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine merge_compare_headers (ier)
!
! Compare observation headers (meta-data) to ensure that the sets of 
! profiles to merge are for the same observation
!
   use m_read_profiles, only : prof_num_2d, prof_num_3d, prof_kmax
   use m_read_profiles, only : prof_names
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c
   use m_obs_list, only : obs_list_tslots, obs_list_subtypes
   use m_obs_list, only : obs_list_format_recs
   use m_obs_list, only : obs_list_names, obs_list_counter
!   
   implicit none
!
   integer, intent(out) :: ier   
!
   logical :: lok
   integer :: n,m
   integer :: nsum
!
   ier=0
!
   if (file1_list_len_i /= obs_list_len_i) then
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_len_i=', &
            file1_list_len_i,obs_list_len_i        
   endif
!  
   if (file1_list_len_r /= obs_list_len_r) then
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_len_r=', &
            file1_list_len_r,obs_list_len_r        
   endif
!  
   if (file1_list_len_c /= obs_list_len_c) then 
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_len_c=', &
            file1_list_len_c,obs_list_len_c        
   endif
!  
   if (file1_kmax /= prof_kmax) then 
     ier=ier+1
     print *,'ERROR: Headers differ: prof_kmax=',file1_kmax,prof_kmax
   endif
!  
   if (file1_list_format_recs /= obs_list_format_recs) then
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_format_recs=', &
            file1_list_format_recs,obs_list_format_recs
   endif
!  
   if (file1_list_tslots /= obs_list_tslots) then
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_tslots=', &
            file1_list_tslots,obs_list_tslots
   endif
!  
   if (file1_list_subtypes /= obs_list_subtypes) then
     ier=ier+1
     print *,'ERROR: Headers differ: obs_list_subtypes=', &
            file1_list_subtypes,obs_list_subtypes
   endif
!
   nsum=obs_list_len_i+obs_list_len_r+obs_list_len_c
   if (file1_list_nsum == nsum) then
     lok=.true.
     do n=1,file1_list_nsum
       if (file1_list_names(n) /= obs_list_names(n)) then
         lok=.false.
         print *,'ERROR: Headers differ: obs_list_names(n)=', &
                n,' ',file1_list_names(n),' ',obs_list_names(n)
       endif
     enddo
     if (.not. lok) then
       ier=ier+1
     endif
   endif
! 
   if (file1_list_tslots == obs_list_tslots .and. &
       file1_list_subtypes == obs_list_subtypes) then 
     lok=.true.
     do n=1,file1_list_subtypes
       do m=1,obs_list_tslots
         if (file1_list_counter(m,n) /= obs_list_counter(m,n)) then
           lok=.false.
           print *,'ERROR: Headers differ: obs_list_counter(n,m))=', &
                n,m,file1_list_counter(n,m),obs_list_counter(n,m)
         endif
       enddo
     enddo
     if (.not. lok) then
       ier=ier+1
     endif
   endif
!
   end subroutine merge_compare_headers
!
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine merge_save_header
!
!  Merge 2 headers of NR field profile files 
!
   use m_read_profiles, only : prof_num_2d, prof_num_3d, prof_kmax
   use m_read_profiles, only : prof_names
   use m_obs_list, only : obs_list_len_i, obs_list_len_r, obs_list_len_c
   use m_obs_list, only : obs_list_tslots, obs_list_subtypes
   use m_obs_list, only : obs_list_format_recs
   use m_obs_list, only : obs_list_names, obs_list_counter
!   
   implicit none
!
   file1_list_len_i=obs_list_len_i
   file1_list_len_r=obs_list_len_r
   file1_list_len_c=obs_list_len_c
   file1_kmax=prof_kmax
   file1_list_format_recs=obs_list_format_recs
   file1_num_2d=prof_num_2d
   file1_num_3d=prof_num_3d
   file1_list_tslots=obs_list_tslots
   file1_list_subtypes=obs_list_subtypes
! 
   file1_list_nsum=file1_list_len_i+file1_list_len_r+file1_list_len_c
   file1_list_names(1:file1_list_nsum)=obs_list_names(1:file1_list_nsum)
   file1_list_counter(1:obs_list_tslots,1:obs_list_subtypes)= &
          obs_list_counter(1:obs_list_tslots,1:obs_list_subtypes)
   file1_prof_names(1:file1_num_2d,2)=prof_names(1:file1_num_2d,2)
   file1_prof_names(1:file1_num_3d,3)=prof_names(1:file1_num_3d,3)
!
   end subroutine merge_save_header
!
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine merge_prof_name_list 
!
!  Merge field name lists from 2 file headers, removing duplicates
!
   use m_read_profiles, only : prof_num_2d, prof_num_3d, prof_kmax
   use m_read_profiles, only : prof_names
!
   implicit none
!
   logical :: duplicate
   integer :: n,m
!
   file2_num_2d=prof_num_2d
   file2_num_3d=prof_num_3d
!
   merge_num_2d=file1_num_2d
   merge_prof_names(1:file1_num_2d,2)=file1_prof_names(1:file1_num_2d,2)
   do n=1,prof_num_2d
     duplicate=.false.
     do m=1,file1_num_2d
       if (prof_names(n,2) == file1_prof_names(m,2)) then
         duplicate=.true.
       endif
     enddo
!
     if (.not. duplicate) then
       merge_num_2d=merge_num_2d+1
       merge_prof_names(merge_num_2d,2)=prof_names(n,2)
       merge_copy(n,2)=.true.
     else
       merge_copy(n,2)=.false.
     endif
!
   enddo
!
   merge_num_3d=file1_num_3d
   merge_prof_names(1:file1_num_3d,3)=file1_prof_names(1:file1_num_3d,3)
   do n=1,prof_num_3d
     duplicate=.false.
     do m=1,file1_num_3d
       if (prof_names(n,3) == file1_prof_names(m,3)) then
         duplicate=.true.
       endif
     enddo
!
     if (.not. duplicate) then
       merge_num_3d=merge_num_3d+1
       merge_prof_names(merge_num_3d,3)=prof_names(n,3)
       merge_copy(n,3)=.true.
     else
       merge_copy(n,3)=.false.
     endif
!
   enddo
!
   print ('(a,3i4)'),'Number of 2d fields (file1,file2,merged)=', &
            file1_num_2d,prof_num_2d,merge_num_2d
   print ('(a,3i4)'),'Number of 3d fields (file1,file2,merged)=', &
            file1_num_3d,prof_num_3d,merge_num_3d
   print ('(8a)'),' 2D field names in file 1: ', &
                  file1_prof_names(1:file1_num_2d,2)
   print ('(8a)'),' 2D field names in file 2: ', &
                  prof_names(1:prof_num_2d,2)
   print ('(8a)'),' 2D field names in merged file: ', &
                  merge_prof_names(1:merge_num_2d,2)
   print ('(8a)'),' 3D field names in file 1: ', &
                  file1_prof_names(1:file1_num_3d,3)
   print ('(8a)'),' 3D field names in file 2: ', &
                  prof_names(1:prof_num_3d,3)
   print ('(8a)'),' 3D field names in merged file: ', &
                  merge_prof_names(1:merge_num_3d,3)
!
   end subroutine merge_prof_name_list 
!
!   
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_prof_merge










 
