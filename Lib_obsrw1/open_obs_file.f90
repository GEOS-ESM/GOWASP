   subroutine open_obs_files (luin,luin_tab,luout,lrad,lgetp,          &
                              lprint,dtype,bufr_in_file,bufr_tab_file, &
                              bufr_out_file,cdtime,obs_file_format,    &
                              obs_file_type,ierr)
!
!  Open observation data files (particularly for rad data types) 
!  for reading and/or writing. 
!  These are BUFR format files, unless the data type is GENRADTXT, in which 
!  case it is a text file. If a separate BUFR table is required, that is 
!  also opened here. Some parameter values that depend on the data type 
!  are obtained from their respective obs data read/write routines.  For 
!  GENRADTXT, extract information from the file header and, if required, 
!  copy it to a file to be written (with the desired value of idate).
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_rad_obs_arrays, only : obs_generic_int
   use m_bufr_rad, only : read_write_obs_rad
   use m_bufr_rad, only : read_write_gmi_1st_msg
!
   implicit none
!
   logical, intent(in)  :: lrad   ! true if obs are radiances
   logical, intent(in)  :: lprint
   logical, intent(in)  :: lgetp  ! true means extract params from r/w routines
   integer, intent(in)  :: luin, luout
   integer, intent(inout) :: luin_tab
   integer, intent(out) :: ierr
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: bufr_in_file
   character(len=*), intent(in) :: bufr_tab_file
   character(len=*), intent(in) :: bufr_out_file
   character(len=*), intent(in) :: cdtime
   character(len=12), intent(out) :: obs_file_format
   character(len=4), intent(out)  :: obs_file_type
!
   logical :: ldum
   integer :: nobs
   integer :: idum
!
   ierr=0
!
   if (trim(dtype) == 'GENRADTXT') then ! Generic radiance text file
     obs_file_type='GTXT'    ! obs data is formatted txt file
     obs_file_format='formatted'  
   else
     obs_file_type='BUFR'    ! obs data in bufr format
     obs_file_format='unformatted'  
   endif  
!
! Open separate bufr table file if required
   if (trim(bufr_tab_file) /= 'none') then 
     open(unit=luin_tab,file =trim(bufr_tab_file), form='formatted')
     if (lprint) then 
       print *,' '
       print ('(3a,i4)'),' Input file=',trim(bufr_tab_file), &
                         ' opened on unit=',luin_tab       
     endif
   else 
     if (luin /= 0) then 
       luin_tab=luin ! same unit as bufr file to be read in
     endif
   endif
! 
! Setup files for reading observations if requested
   if (luin /= 0 .and. trim(bufr_in_file) /= 'none') then 
     open(unit=luin,file=trim(bufr_in_file),form=trim(obs_file_format))
     if (lprint) then
       print *,' '
       print ('(3a,i4)'),' Input file=',trim(bufr_in_file), &
                         ' opened on unit=',luin
     endif
!
     if (obs_file_type == 'BUFR') then ! obs file in BUFR format 
       call datelen (10) ! sets idate to 10 digit (YYYYMMDDHH)
       call openbf(luin, 'IN ',luin_tab)
     endif
!
! Special instructions to read 1st message in GMI bufr data
      if (trim(dtype) == 'GMI') then
        call read_write_gmi_1st_msg (luin,.true.,idum,ierr)
      endif
!
! Get some obs class dependent parameters (as requested by nobs=-1).
! No reading of a file is performed, unless the format is generic rad text.
! In the latter case, only the file header is read. 
     if (lrad .and. (obs_file_type == 'GTXT' .or. &
                       (obs_file_type == 'BUFR' .and. lgetp) ) ) then  
       nobs=-1  
       call read_write_obs_rad (luin,lprint,dtype,nobs,.true., &
                               .true.,ldum,ierr)
     endif
   endif
!
! Setup files for writing observations if requested
   if (luout /= 0 .and. trim(bufr_out_file) /= 'none') then 
     open(unit=luout,file=trim(bufr_out_file),form=trim(obs_file_format))
     print ('(3a,i3)'),' output file=',trim(bufr_out_file), &
             ' opened on unit=',luout
!
     if (obs_file_type == 'BUFR') then ! obs file in BUFR format 
       call openbf(luout,'OUT',luin_tab)
       call cmpmsg('Y')     ! compress output messages
       call maxout(200000)  ! increase size of ouput bufr
     endif
!
     if (obs_file_type == 'GTXT') then  ! write gen rad text file header
       nobs=-1  ! flag that indicates call to get header names only
       if (trim(cdtime) /= 'none') then ! change idate on new obs
         read (cdtime,'(i10)') obs_generic_int(2) 
       endif
       call read_write_obs_rad (luout,lprint,dtype,nobs,.true., &
                                .false.,ldum,ierr)
     endif 
!
   endif   
!
   end subroutine open_obs_files 
