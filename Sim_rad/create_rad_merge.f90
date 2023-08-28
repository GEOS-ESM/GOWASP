   program rad_merge
!
! Merge two profile data sets that have been produced for the same obs set
! For duplicate field names in the 2 files, only the profiles of those fields 
! in the first file are copied. The obs_info corresponding to each profile pair 
! to be merged must be identical. 
!
! Initial code: Ronald Errico October 8 2014 
! 
   use m_read_profiles, only: read_profiles_setup 
   use m_read_profiles, only: read_profiles_recs
   use m_read_profiles, only: read_profiles_cleanup
!
   use m_read_profiles, only: prof_num_2d, prof_num_3d, prof_kmax
   use m_read_profiles, only: prof_akbk_int, prof_all
   use m_read_profiles, only: prof_names, prof_common_path, prof_dim1
   use m_read_profiles, only: prof_format_header, prof_format_recs
!
   use m_write_profiles, only: write_profiles_header
   use m_write_profiles, only: write_profiles_recs
   use m_write_profiles, only: write_profiles_close
!
   use m_nr_fields_info, only : field_num_2d, field_num_3d
   use m_nr_fields_info, only : field_kmax, field_akbk, field_names
   use m_nr_fields_info, only : field_common_path
!
   use m_obs_list, only : obs_list_tslots, obs_list_subtypes
   use m_obs_list, only : obs_list_counter, obs_list_unit_in
!
   use m_set_unit_nums, only : un_prof_in1, un_prof_in2 
   use m_set_unit_nums, only : un_prof_out0, un_prof_out
!
   use m_prof_merge, only : merge_save_header 
   use m_prof_merge, only : merge_compare_headers
   use m_prof_merge, only : merge_allocate_prof_arrays 
   use m_prof_merge, only : merge_prof_name_list 
   use m_prof_merge, only : merge_profile_recs
   use m_prof_merge, only : merge_num_2d, merge_num_3d
   use m_prof_merge, only : merge_prof_names
!
   implicit none
!
   integer, parameter :: unit_in_1=un_prof_in1
   integer, parameter :: unit_in_2=un_prof_in2
   integer, parameter :: unit_out=un_prof_out0
   integer     :: n1,n2,n3,n4
   integer     :: iers
   integer     :: argc
   integer(4)  :: iargc
   integer     :: ntypes
   character(len=240) :: filein1,filein2,fileout
   character(len=1)   :: ctest
!  
   print *,' '
   argc = iargc()
   if (argc /= 4) then
     print *,' Usage must be: prog.x filein1 filein2 fileout ctest'
     stop
   endif
   call GetArg( 1_4, filein1)
   call GetArg( 2_4, filein2)
   call GetArg( 3_4, fileout)
   call GetArg( 4_4, ctest)
   print *,'filein=',trim(filein1)
   print *,'filein=',trim(filein2)
   print *,'fileout=',trim(fileout)
!
   if (filein1 == fileout .or. filein2 == fileout) then
     print *,'Output file name same as an input file name'
     stop
   endif
!
! Set unit for output profile file (the following identity must be held)
   un_prof_out=unit_out 
!
! Read header from first file
   obs_list_unit_in=unit_in_1
   call read_profiles_setup (filein1,1,iers)
   if (iers /= 0) then
     print *,'Error in reading header on file 1: iers=',iers
     stop
   endif
!
! Save some information from header of first file
   call merge_save_header
   call read_profiles_cleanup (.false.) ! don't close unit
!
! Read header from second file
   obs_list_unit_in=unit_in_2
   call read_profiles_setup (filein2,1,iers)
   if (iers /= 0) then
     print *,'Error in reading header on file 2: iers=',iers
     stop
   endif
!   
! Check compatability of files
   call merge_compare_headers (iers)
   if (iers /= 0) then
     print *,'Headers on 2 files incompatable: iers=',iers
     stop
   endif
!
! Merge names of profiles
   call merge_prof_name_list
!
! Write new header
   ntypes=obs_list_subtypes
   field_num_2d=merge_num_2d
   field_num_3d=merge_num_3d
   field_kmax=prof_kmax
   field_akbk(1:field_kmax+1,1:2)=prof_akbk_int(1:field_kmax+1,1:2)
   field_names(1:field_num_2d,1,2)=merge_prof_names(1:field_num_2d,2)
   field_names(1:field_num_3d,1,3)=merge_prof_names(1:field_num_3d,3)
   field_common_path=prof_common_path
   call write_profiles_header (fileout,prof_format_header, &
                               prof_format_recs)
!    
! Allocate arrays to hold profiles
   call merge_allocate_prof_arrays 
!
   do n1=1,ntypes
     print ('(a,i3,a,i6)'),'Loop 1: n1=',n1,' of ',ntypes
     do n2=1,obs_list_tslots-1
       print ('(a,i3,a,i6)'),'Loop 2: n2=',n2,' of ',obs_list_tslots-1
       do n3=1,obs_list_counter(n2,n1)
         call merge_profile_recs (unit_in_1,unit_in_2,unit_out,n1,n2,n3,ctest)
       enddo 
     enddo
   enddo
!
   print *,' '
   print *,'Files merged'
!
   call read_profiles_cleanup (.true.)
   call write_profiles_close
!
   end program rad_merge
