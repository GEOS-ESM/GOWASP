   program rad_reord
!
! Reorder profiles from an order based on time slot to an
! order based on obs subtype (e.g., satellite platform)
!
!  Initial Code: Ronald Errico August 30 2014
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
   use m_set_unit_nums, only : un_prof_in1
   use m_set_unit_nums, only : un_prof_out0, un_prof_out
!
   implicit none
!
   integer     :: n1,n2,n3,n4
   integer     :: iers
   integer     :: argc
   integer(4)  :: iargc
   integer ntypes
   character(len=240) filein,fileout
 !  
   print *,' '
   argc = iargc()
   if (argc /= 2) then
     print *,' Usage must be: prog.x filein fileout'
     stop
   endif
   call GetArg( 1_4, filein)
   call GetArg( 2_4, fileout)
   print *,'filein=',trim(filein)
   print *,'fileout=',trim(fileout)
!
! Set unit number for output
   un_prof_out=un_prof_out0
!
! copy header
   obs_list_unit_in=un_prof_in1
   call read_profiles_setup (filein,1,iers)
   if (iers /= 0) then
     print *,'Error in reading header: iers=',iers
     stop
   endif
!
   ntypes=obs_list_subtypes
   field_num_2d=prof_num_2d
   field_num_3d=prof_num_3d
   field_kmax=prof_kmax
   field_akbk(1:field_kmax+1,1:2)=prof_akbk_int(1:field_kmax+1,1:2)
   field_names(1:field_num_2d,1,2)=prof_names(1:field_num_2d,2)
   field_names(1:field_num_3d,1,3)=prof_names(1:field_num_3d,3)
   field_common_path=prof_common_path
   call write_profiles_header (fileout,prof_format_header, &
                               prof_format_recs)
   call read_profiles_cleanup (.true.)  ! .true. causes unit to be closed
!    
   do n1=1,ntypes
     call read_profiles_setup (filein,1,iers)
     print ('(a,i3,a,i8)'),'Loop 1: n1=',n1,' of ',ntypes
     do n2=1,obs_list_tslots-1
       print ('(a,i3,a,i8)'),'Loop 2: n2=',n2,' of ',obs_list_tslots-1
       do n3=1,ntypes
         print ('(a,i3,a,2i8)'),'Loop 3: n3=',n3,' of ',ntypes, &
                                obs_list_counter(n2,n3)
         do n4=1,obs_list_counter(n2,n3)
           call read_profiles_recs (1,.true.,iers)
           if (n3 == n1) then
             call write_profiles_recs (prof_dim1,1,prof_all)
           endif
         enddo
       enddo 
     enddo
     call read_profiles_cleanup (.true.) ! .true. causes unit to be closed
   enddo
!
   print *,' '
   print *,'File reordered'
!
   call write_profiles_close
!
   end program rad_reord
