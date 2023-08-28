!
      subroutine sondes_output (ndim1,ndim2,ndim3,unit_out,       &
                                itype,nob,lprint,lsample,icounts, & 
                                obs_info,snd_old,snd_raw) 
!
! Create observation reports from raw soundings for RAOB, PIBAL, and DROPSONDE 
! Output reports to buffer file
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
       use m_kinds, only      : rkind1, rkind2
       use m_conv_names, only : bbnlold, bbnlnew
       use m_conv_names, only : obs_info_names, conv_nhead
       use m_bufr_conv, only  : conv_print_sample, conv_rw
       use m_bufr_conv, only  : conv_nlevs, conv_bmiss
       use m_bufr_conv, only  : conv_values, conv_info
!
       implicit none
!
       logical, intent(in) :: lprint
       logical, intent(in) :: lsample
       integer, intent(in) :: ndim1,ndim2,ndim3
       integer, intent(in) :: unit_out
       integer, intent(in) :: itype   ! obs type (ncep kx )
       integer, intent(in) :: nob
       integer, intent(inout) :: icounts(2,2)
       real(rkind1), intent(in) :: snd_old(ndim1,ndim2)     ! read sounding
       real(rkind2), intent(inout) :: obs_info(ndim3)
       real(rkind1), intent(inout) :: snd_raw(ndim1,ndim2)  ! raw sounding
!
       logical, parameter :: lread=.false.   ! .false. means write bufr data
       logical, parameter :: lpbufr=.false.  ! also print output in bufr sub
       logical :: ldum                       ! argument not used here
!
       integer, parameter :: ndum=1    ! dimension of dummy argument array
       integer, parameter :: isample=4 ! print isample of each obs itype  
       integer :: i
       integer :: nlold
       integer :: nlnew
       integer :: nraw
       integer :: ier
       integer :: ilook
       integer :: idum(ndum)         ! dummy argument array
       real(rkind1) :: snd_rpt(ndim1,ndim2)
       character(len=4) :: dtype
!
       nlold=nint(obs_info(bbnlold))
       nraw=nint(obs_info(bbnlnew))
!
       if (itype > 199) then
         call sonde_report_wind (ndim1,ndim2,itype,nraw,nlnew, &
                                 snd_raw,snd_rpt)
         dtype='WIND'
       else
         call sonde_report_mass (ndim1,ndim2,itype,nraw,nlnew, &
                                 snd_raw,snd_rpt)
         dtype='MASS'
       endif
       call siglevs_qc (ndim1,ndim2,nlold,nlnew,conv_bmiss,dtype, &
                        snd_old,snd_rpt)
!
       icounts(1,2)=icounts(1,2)+1
       icounts(2,2)=icounts(2,2)+nlnew*2
!
       conv_nlevs=nlnew
       conv_info(1:conv_nhead)=obs_info(1:conv_nhead)
       conv_values(:,:)=snd_rpt(:,:)
!
       call conv_rw (lprint,lpbufr,dtype,nob,lread,ndum,idum,ldum,ier)
!
       ilook=max(icounts(1,1)/isample,1)
       if (lsample .and. mod(icounts(1,2),ilook) == 1) then
!
          print *,' '
          print ('(2a,i7)'),'Sample of input and output obs to be', &
                            ' printed for nob=',nob
          print ('(a)'),'obs_info_extra='
          print ('(4(a10,f12.5))'), &
            (obs_info_names(i),obs_info(i),i=conv_nhead+1,ndim3)
!
          call conv_print_sample (nlnew,nob,lread, &
                                  'simulated observation report')
          conv_nlevs=nint(obs_info(bbnlnew))
          conv_values(:,:)=snd_raw(:,:)
          call conv_print_sample (nraw,nob,lread, &
                                  'simulated raw sounding')
          conv_nlevs=nint(obs_info(bbnlold))
          conv_values(:,:)=snd_old(:,:)
          call conv_print_sample (nlold,nob,lread, &
                                  'original read (real) sounding')
       endif
! 
       end subroutine sondes_output
! 
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine non_sonde_output (ndim1,ndim2,ndim3,unit_out,       &
                                   itype,nob,lprint,lsample,icounts, & 
                                   obs_info,obs_old,obs_new) 
!
! Output all requested conventional reports to a BUFR file that are not 
! sondes of some type  
!
! Code History:
!
!  Ronald Errico       23 Sept 2014    Initial Code
!
       use m_kinds, only      : rkind1, rkind2
       use m_conv_names, only : bbnlold, bbnlnew
       use m_conv_names, only : obs_info_names, conv_nhead
       use m_bufr_conv, only  : conv_print_sample, conv_rw
       use m_bufr_conv, only  : conv_nlevs, conv_bmiss
       use m_bufr_conv, only  : conv_values, conv_info
!
       implicit none
!
       logical, intent(in) :: lprint
       logical, intent(in) :: lsample
       integer, intent(in) :: ndim1,ndim2,ndim3
       integer, intent(in) :: unit_out
       integer, intent(in) :: itype   ! obs type (ncep kx )
       integer, intent(in) :: nob
       integer, intent(inout) :: icounts(2,2)
       real(rkind1), intent(in) :: obs_old(ndim1,ndim2)     ! read sounding
       real(rkind2), intent(inout) :: obs_info(ndim3)
       real(rkind1), intent(inout) :: obs_new(ndim1,ndim2)  ! raw sounding
!
       logical, parameter :: lread=.false.   ! .false. means write bufr data
       logical, parameter :: lpbufr=.false.  ! also print output in bufr sub
       logical :: ldum 
!
       integer, parameter :: ndum=1    ! dimension of dummy array
       integer, parameter :: isample=4 ! print isample of each obs itype  
       integer :: i
       integer :: nlold
       integer :: nlnew
       integer :: ier
       integer :: ilook
       integer :: idum(ndum)         ! dummy array
       real(rkind1) :: snd_rpt(ndim1,ndim2)
       character(len=4) :: dtype
!
       nlold=nint(obs_info(bbnlold))
       nlnew=nint(obs_info(bbnlnew))
!
       if (itype > 199) then
         dtype='WIND'
       else
         dtype='MASS'
       endif
       call siglevs_qc (ndim1,ndim2,nlold,nlnew,conv_bmiss,'NOSONDE', &
                        obs_old,obs_new)
!
       icounts(1,2)=icounts(1,2)+1
       icounts(2,2)=icounts(2,2)+nlnew*2
!
       conv_nlevs=nlnew
       conv_info(1:conv_nhead)=obs_info(1:conv_nhead)
       conv_values(:,:)=obs_new(:,:)
!
       call conv_rw (lprint,lpbufr,dtype,nob,lread,ndum,idum,ldum,ier)
!
       ilook=max(icounts(1,1)/isample,1)
       if (lsample .and. mod(icounts(1,2),ilook) == 1) then
!
          print *,' '
          print ('(2a,i7)'),'Sample of input and output obs to be', &
                            ' printed for nob=',nob
          print ('(a)'),'obs_info_extra='
          print ('(4(a10,f12.5))'), &
            (obs_info_names(i),obs_info(i),i=conv_nhead+1,ndim3)
!
          call conv_print_sample (nlnew,nob,lread, &
                                  'simulated observation report')
          conv_nlevs=nint(obs_info(bbnlold))
          conv_values(:,:)=obs_old(:,:)
          call conv_print_sample (nlold,nob,lread, &
                                  'original read (real) sounding')
       endif
! 
       end subroutine non_sonde_output
! 
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
       subroutine sonde_report_wind (ndim1,ndim2,itype,nraw,nrpt, &
                                     snd_raw,snd_rpt)
!
! Process raw report of sonde wind sounding to create obs on reported levels
! 
! Code History:
!
!  Ronald Errico       23 Sept 2014    Initial Code
!
       use m_kinds, only : rkind1
       use m_parameters, only : pifac_k1,rpifac_k1  
       use m_conv_names, only : bbp,bbu,bbv,bbc
!
       implicit none
!
       integer, intent(in)  :: ndim1,ndim2
       integer, intent(in)  :: nraw
       integer, intent(in)  :: itype   ! obs type (ncep kx )
       integer, intent(out) :: nrpt
       real(rkind1), intent(out) :: snd_rpt(ndim1,ndim2)   ! processed sounding
       real(rkind1), intent(inout) :: snd_raw(ndim1,ndim2) ! raw sounding
!
       logical :: lascent
       integer :: n
       real(rkind1) :: wspd2,wspd,wdir,utmp,vtmp
       real(rkind1) :: dum  ! used as dummy argument
! 
       snd_rpt(:,:)=0.
       nrpt=0                 ! initialize counter
!
! Sort sounding levels so descending in p values. like raob
       call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'des')
!
! convert wind vectors to wind speed and direction; also change Pa to mb
       do n=1,nraw
         wspd2=snd_raw(bbu,n)*snd_raw(bbu,n)+snd_raw(bbv,n)*snd_raw(bbv,n)
         wspd=sqrt(max(wspd2,1.e-10))
         wdir=atan2(snd_raw(bbu,n),snd_raw(bbv,n))*rpifac_k1
         snd_raw(bbu,n)=wspd
         snd_raw(bbv,n)=mod(wdir,360.)
         snd_raw(bbp,n)=snd_raw(bbp,n)*0.01
       enddo
!
! Determine standard levels within sounding
       if (nraw > 1) then
         call siglevs_mandatory (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt)
       endif
!
! Include first and last levels, mostly used for interpolation checking
! may want to remove these from the wind sounding at the end
       nrpt=nrpt+1
       snd_rpt(:,nrpt)=snd_raw(:,1)
       snd_rpt(bbc,nrpt)=0.     ! flag for surface
!
       if (nraw > 1) then
         nrpt=nrpt+1
         snd_rpt(:,nrpt) = snd_raw(:,nraw)
         snd_rpt(bbc,nrpt)=1.   ! top of sounding designated as a mandatory level       
       endif
!
! Include winds that are relative maxima with respect to adjacent levels.
! Only assign as max wind level if speed is 15 m/s or greater
       do n=2,nraw-1
         if (nrpt == ndim2) exit  ! no array space for additional obs
         if (snd_raw(bbu,n) > snd_raw(bbu,n-1) .and. &
             snd_raw(bbu,n) > snd_raw(bbu,n+1) .and. &
             snd_raw(bbu,n) >= 15. ) then
           nrpt=nrpt+1
           snd_rpt(:,nrpt)=snd_raw(:,n)
           snd_rpt(bbc,nrpt)=4   ! this is improperly assigned as 4
         endif
       enddo
!
! Include sounding wind speeds that cannot be linearly reconstruncted by obs
       call siglevs_linear (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt, & 
                            dum,'spd')
!
! Include sounding wind directions that cannot be linearly reconstruncted by obs
       call siglevs_linear (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt, & 
                            dum,'dir')
!
! Perform final sort, with order appropriate for sond type
       if (itype == 232) then  ! dropsonde
         call siglevs_sort (ndim1,ndim2,nrpt,bbp,snd_rpt,'asc')
         call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'asc')
       else                    ! raob or pibal
         call siglevs_sort (ndim1,ndim2,nrpt,bbp,snd_rpt,'des')
         call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'des')
       endif 
!
! Convert speed and direction back to vector wind and p back to Pa from mb
       do n=1,nrpt
         utmp=snd_rpt(bbu,n)*sin(snd_rpt(bbv,n)*pifac_k1)
         vtmp=snd_rpt(bbu,n)*cos(snd_rpt(bbv,n)*pifac_k1)
         snd_rpt(bbu,n)=utmp
         snd_rpt(bbv,n)=vtmp
         snd_rpt(bbp,n)=snd_rpt(bbp,n)*100.
       enddo
!
! Convert same as above for raw sounding
       do n=1,nraw
         utmp=snd_raw(bbu,n)*sin(snd_raw(bbv,n)*pifac_k1)
         vtmp=snd_raw(bbu,n)*cos(snd_raw(bbv,n)*pifac_k1)
         snd_raw(bbu,n)=utmp
         snd_raw(bbv,n)=vtmp
         snd_raw(bbp,n)=snd_raw(bbp,n)*100.
       enddo
!
       end subroutine sonde_report_wind
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
       subroutine sonde_report_mass (ndim1,ndim2,itype,nraw,nrpt, &
                                     snd_raw,snd_rpt)
!
!  Process raw report of sonde mass sounding to create obs on reported levels
! 
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
       use m_kinds, only : rkind1
       use m_conv_names, only : bbp,bbt,bbq,bbc
!
       implicit none
!
       integer, intent(in)  :: ndim1,ndim2
       integer, intent(in)  :: nraw
       integer, intent(in)  :: itype   ! obs type (ncep kx )
       integer, intent(out) :: nrpt
       real(rkind1), intent(out) :: snd_rpt(ndim1,ndim2)   ! processed sounding
       real(rkind1), intent(inout) :: snd_raw(ndim1,ndim2) ! raw sounding
!
       integer :: n_trop
       integer :: n
       real(rkind1) :: rh
       real(rkind1) :: p_trop  ! p at tropopause
!
       logical :: found_1
       character(len=*), parameter :: myname='siglevs_mass'
! 
       snd_rpt(:,:)=0.
       nrpt=0                 ! initialize counter
!
! Sort soundind levels so descending in p values. like raob
       call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'des')
!
! convert q to rh and change Pa to mb
       do n=1,nraw
         call transform_q_rh (snd_raw(bbp,n),snd_raw(bbt,n), &
                              snd_raw(bbq,n),rh,'q2rh')
         snd_raw(bbq,n)=rh
         snd_raw(bbp,n)=snd_raw(bbp,n)*0.01
       enddo
!
! Determine standard levels within sounding
       if (nraw > 1) then
         call siglevs_mandatory (ndim1,ndim2,nraw,nrpt, &
                                 snd_raw,snd_rpt)
       endif
!
! Include first and last levels, mostly used for interpolation checking
! may want to remove these from the wind sounding at the end
       nrpt=nrpt+1
       snd_rpt(:,nrpt)=snd_raw(:,1)
       snd_rpt(bbc,nrpt)=0.     ! flag for surface
!
       if (nraw > 1) then
         nrpt=nrpt+1
         snd_rpt(:,nrpt) = snd_raw(:,nraw)
         snd_rpt(bbc,nrpt)=1.   ! top of sounding designated as a mandatory level       
       endif
!
! Add obs at tropopause level if level found
       p_trop=snd_raw(bbp,1)   ! default value
       if (nraw > 3) then
         call siglevs_trop (ndim1,ndim2,nraw,snd_raw,n_trop)
         if (n_trop > 0) then
           nrpt=nrpt+1
           snd_rpt(:,nrpt)=snd_raw(:,n_trop)
           snd_rpt(bbc,nrpt)=5.
           p_trop=snd_raw(bbp,n_trop)
         endif
       endif
!
! Add a level between 110 and 100 mb if there is any in raw sounding
! Use the raw sounding level closest to 105 mb
       found_1=.false.
       do n=1,nraw
         if (snd_raw(bbp,n) > 100. .and. snd_raw(bbp,n) <= 110.) then
           if (.not. found_1) then
             nrpt=nrpt+1
             snd_rpt(:,nrpt)=snd_raw(:,n)
             snd_rpt(10,nrpt)=1
             found_1=.true.
           elseif (abs(snd_rpt(bbp,nrpt)-105.) > &
                   abs(snd_raw(bbp,n)-105.)) then
             snd_rpt(:,nrpt)=snd_raw(:,n)
             snd_rpt(bbc,nrpt)=1
           endif
         endif
       enddo
!
! Include sounding T that cannot be linearly reconstruncted by obs
       call siglevs_linear (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt, & 
                            p_trop,'t')
!
! Include sounding q that cannot be linearly reconstruncted by obs
       call siglevs_linear (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt, & 
                            p_trop,'q')
!
! Perform final sort, with order appropriate for sond type
       if (itype == 132) then  ! dropsonde
         call siglevs_sort (ndim1,ndim2,nrpt,bbp,snd_rpt,'asc')
         call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'asc')
       else                    ! raob 
         call siglevs_sort (ndim1,ndim2,nrpt,bbp,snd_rpt,'des')
         call siglevs_sort (ndim1,ndim2,nraw,bbp,snd_raw,'des')
       endif 
!
! Convert p back to Pa from mb and rh to q
       do n=1,nrpt
         snd_rpt(bbp,n)=snd_rpt(bbp,n)*100.
         rh=snd_rpt(bbq,n)
         call transform_q_rh (snd_rpt(bbp,n),snd_rpt(bbt,n), &
                              snd_rpt(bbq,n),rh,'rh2q')
       enddo 
!
! Convert same as above for raw sounding
       do n=1,nraw
         snd_raw(bbp,n)=snd_raw(bbp,n)*100.
         rh=snd_raw(bbq,n)
         call transform_q_rh (snd_raw(bbp,n),snd_raw(bbt,n), &
                              snd_raw(bbq,n),rh,'rh2q')
       enddo
!
       end subroutine sonde_report_mass
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!

       subroutine siglevs_mandatory (ndim1,ndim2,nraw,nrpt, &
                                     snd_raw,snd_rpt)
!
! Determine sonde report at mandatory pressure levels from raw sounding
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
       use m_kinds, only : rkind1
       use m_parameters, only : pifac_k1,rpifac_k1
       use m_conv_names, only : bbp, bbz, bbt, bbq, bbu, bbv
       use m_conv_names, only : bbx, bby, bbr, bbc
!
       implicit none
!
       integer, intent(in)  :: ndim1,ndim2
       integer, intent(in)  :: nraw
       integer, intent(inout) :: nrpt
       real(rkind1), intent(in)  :: snd_raw(ndim1,ndim2)  ! raw sounding
       real(rkind1), intent(out) :: snd_rpt(ndim1,ndim2)  ! processed sounding
!
       integer, parameter :: stnd_dim=16
       integer :: m,n,n1
       real(rkind1) :: fac,fac1
       real(rkind1), dimension(stnd_dim),parameter :: stnd=(/1000.,925.,850., &
             700.,500.,400.,300.,250.,200.,150.,100.,70.,50.,30.,20.,10./)
! stnd = mandatory levels in mb in descending order
!
! Find first mandatory level in sounding
       m=1
       do n=1,stnd_dim-1
         if (snd_raw(bbp,1) < stnd(n)) then
           m=n+1
         endif
       enddo 
!
! Ascend through sounding looking sequentially for mandatory levels
       do n=2,nraw
         n1=n-1
!
         if (snd_raw(bbp,n) <= stnd(m) ) then
           fac=log(stnd(m)/snd_raw(bbp,n))/log(snd_raw(bbp,n1)/snd_raw(bbp,n))
           fac1=1.-fac
           nrpt=nrpt+1
           snd_rpt(bbx,nrpt)=fac1*snd_raw(bbx,n)+fac*snd_raw(bbx,n1)
           snd_rpt(bby,nrpt)=fac1*snd_raw(bby,n)+fac*snd_raw(bby,n1)
           snd_rpt(bbr,nrpt)=fac1*snd_raw(bbr,n)+fac*snd_raw(bbr,n1)
           snd_rpt(bbp,nrpt)=stnd(m)
           snd_rpt(bbz,nrpt)=fac1*snd_raw(bbz,n)+fac*snd_raw(bbz,n1)
           snd_rpt(bbt,nrpt)=fac1*snd_raw(bbt,n)+fac*snd_raw(bbt,n1)
           snd_rpt(bbq,nrpt)=fac1*snd_raw(bbq,n)+fac*snd_raw(bbq,n1)
           snd_rpt(bbu,nrpt)=fac1*snd_raw(bbu,n)+fac*snd_raw(bbu,n1)
           snd_rpt(bbv,nrpt)=fac1*snd_raw(bbv,n)+fac*snd_raw(bbv,n1)
           snd_rpt(bbc,nrpt)=1.      ! flag for standard levels
           m=m+1
           if (m > stnd_dim) exit
         endif
       enddo
!
       end subroutine siglevs_mandatory  
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_linear (ndim1,ndim2,nraw,nrpt,snd_raw,snd_rpt, & 
                                p_trop,cfield)
!
! Add observations for sounding levels that cannot be reconstructed 
! well by linear interpolation between other observation levels  
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1
     use m_conv_names, only : bbp, bbz, bbu, bbv, bbt, bbq, bbc
     implicit none
!
     integer, intent(in) :: ndim1,ndim2
     integer, intent(in) :: nraw
     integer, intent(inout) :: nrpt
     real(rkind1), intent(in)  :: snd_raw(ndim1,ndim2)
     real(rkind1), intent(in) :: p_trop
     real(rkind1), intent(out) :: snd_rpt(ndim1,ndim2)
     character(len=*), intent(in) :: cfield  ! field to reconstruct
!
     logical :: lcond(4)
     integer :: klast
     integer :: nf
     integer :: k, m, n, mlev, diff_id
     integer :: cat_code
     real*4  :: diff, diff_max, fac, int_val
     real(rkind1), parameter :: crit_ws=1.5   ! wind spd diff max (speed)
     real(rkind1), parameter :: crit_wd=2.0   ! wind dir diff max (degrees)
     real(rkind1), parameter :: crit_tb=0.08  ! temp diff max below trop
     real(rkind1), parameter :: crit_ta=0.4   ! temp diff max above trop
     real(rkind1), parameter :: crit_q=0.015  ! rh diff max 
!
     if (cfield == 'spd') then      ! wind speed carried in u field array
       nf=bbu    
       cat_code=3
     elseif (cfield == 'dir') then  ! wind dir carried in v field array
       nf=bbv    
       cat_code=3
     elseif (cfield == 't') then 
       nf=bbt    
       cat_code=2
     elseif (cfield == 'q') then
       nf=bbq    
       cat_code=2
     endif
!
! Loop until all soundings are within criteria
! (ntsnd2-ct)/2 is the number of possible additional significant levels
! considering that this routine is called twice
     klast=(ndim2-nrpt)/2  
     do k=1,klast 
       if (nrpt == ndim2) exit  ! no more array space for additional obs
!
! First sort previously defined obs by ascending heights
       call siglevs_sort (ndim1,ndim2,nrpt,bbz,snd_rpt,'asc')
!
! Find sounding level that exhibis the worst approximation when  
! reconstructed by linear interpolation in z from existing obs.
       diff_id=1
       diff_max=0.
       mlev=1             
       do n=1,nrpt-1  ! index of raw soundings to consider
! 
! Look for reported obs just below sounding level considered
! (this will be the last value of m determined)
         do m=1,nrpt-1
           if (snd_raw(bbz,n) >= snd_rpt(bbz,m) ) then
             mlev=m
           endif
         enddo
!
! Create and compare interpolated value with sounding value for field nf
         fac=(snd_rpt(bbz,mlev+1)-snd_raw(bbz,n))/(snd_rpt(bbz,mlev+1)-snd_rpt(bbz,mlev))
         int_val=(1.-fac)*snd_rpt(nf,mlev+1)+fac*snd_rpt(nf,mlev)
!
         if (nf == bbv) then   ! for wind direction check 
           if (snd_raw(nf,n) >= int_val) then
             diff=min(snd_raw(nf,n)-int_val,360.+int_val-snd_raw(nf,n))
           else
             diff=min(int_val-snd_raw(nf,n),360.+snd_raw(nf,n)-int_val)
           endif
         else 
           diff=abs(snd_raw(nf,n)-int_val)
         endif
!
         if (diff > diff_max) then
           diff_max=diff
           diff_id=n
         endif
! 
       enddo    ! loop over sounding levels
!
! Add obs at new level if criteria met
       lcond(1)= nf == bbu .and. diff_max > crit_ws .and. snd_raw(bbp,diff_id) > 1.
       lcond(2)= nf == bbv .and. diff_max > crit_wd .and. snd_raw(bbp,diff_id) > 1.
       lcond(3)= nf == bbq .and. diff_max > crit_q 
       lcond(4)= nf == bbt .and. (                                               &
                 ( diff_max > crit_tb .and. snd_raw(bbp,diff_id) >= p_trop) .or. & 
                 ( diff_max > crit_ta .and. snd_raw(bbp,diff_id) < p_trop) ) 
       if (lcond(1) .or. lcond(2) .or. lcond(3) .or. lcond(4)) then
         nrpt=nrpt+1
         snd_rpt(:,nrpt)=snd_raw(:,diff_id)
         snd_rpt(bbc,nrpt)=cat_code
       else 
         exit  ! finished, since no new obs found
       endif
!
     enddo     ! loop over levels in raw sounding
!
     end subroutine siglevs_linear
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_sort (ndim1,ndim2,ndat,nf,dat,direction)
!
!  Sort array according to values in field index nf
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1
     implicit none
!
     character(len=*), intent(in) :: direction  ! 'des' for descending values
     integer, intent(in) :: ndim1,ndim2
     integer, intent(in) :: ndat                ! number of levels to consider
     integer, intent(in) :: nf                  ! index of field to base order
     real(rkind1), intent(inout) :: dat(ndim1,ndim2)    ! array to re-order
!
     integer :: i,j,iswap
     real(rkind1) :: datm,d1
     real(rkind1) :: temp(ndim1)
     intrinsic maxloc,minloc
!
     do i=1,ndat-1
       if (direction == 'des') then  ! look for next max value
         datm=-1.e10     
         do j=i,ndat
           d1=dat(nf,j)
           if (d1 > datm) then
             datm=d1
             iswap=j
           endif
         enddo
       else                          ! look for next min value
         datm=1.e10  
         do j=i,ndat
           d1=dat(nf,j)
           if (d1 < datm) then
             datm=d1
             iswap=j
           endif
         enddo
       endif
!   
       if (iswap /= i) then
         temp(:)=dat(:,i)
         dat(:,i)=dat(:,iswap)
         dat(:,iswap)=temp(:)
       endif
     enddo
!
     end subroutine siglevs_sort
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_trop (ndim1,ndim2,nsnd,snd,n_trop)
!
! Find level of tropopause
!
! Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1
     use m_conv_names, only : bbp,bbz,bbt
!
     implicit none 
!
     integer, intent(in)  :: ndim1,ndim2,nsnd     
     real(rkind1), intent(in)  :: snd(ndim1,ndim2)
     integer, intent(out)  :: n_trop
!
     integer  :: n,k
     real(rkind1) :: lapse(nsnd)
     real(rkind1) :: lapse_avg
     real(rkind1) :: sig_slope
     real(rkind1) :: xtemp5(5), ytemp5(5)
!
     n_trop = 0
!
! Calculate lapse rates (units K/km) for levels above 500 mb 
     lapse(:)=-100.
     do n=3,(nsnd-2)
       if (snd(bbp,n) <= 500.) then
         do k=1,5
           xtemp5(k)=snd(bbz,n+k-3)
           ytemp5(k)=snd(bbt,n+k-3)
         enddo
         call siglevs_slope (5,xtemp5,ytemp5,sig_slope)
         lapse(n)=1000.*sig_slope
       endif
     enddo
     if (snd(bbp,1) <= 500.) then 
       lapse(1)=1000.*(snd(bbt,2)-snd(bbt,1))/ &
                      (snd(bbz,2)-snd(bbz,1))
     endif
     if (snd(bbp,2) <= 500.) then 
       do k=1,3
         xtemp5(k)=snd(bbz,k)
         ytemp5(k)=snd(bbt,k)
       enddo
       call siglevs_slope (3,xtemp5,ytemp5,sig_slope)
       lapse(2)=1000.*sig_slope
     endif
     if (snd(bbp,nsnd-1) <= 500.) then
       do k=1,3
         xtemp5(k)=snd(bbz,nsnd+k-3)
         ytemp5(k)=snd(bbt,nsnd+k-3)
       enddo
       call siglevs_slope (3,xtemp5,ytemp5,sig_slope)
       lapse(nsnd-1)=1000.*sig_slope
     endif
     if (snd(bbp,nsnd) <= 500.) then
       lapse(nsnd)=1000.*(snd(bbt,nsnd)-snd(bbt,nsnd-1))/ &
                         (snd(bbz,nsnd)-snd(bbz,nsnd-1))
     endif
!
! If no levels meet above criteria, then find largest p such 
! that dT/dz >= -2K/km
     n_trop=nsnd
     do n=1,nsnd-1,1
       if (lapse(n) >= -2.) then
         n_trop=n
         exit
       endif
     enddo
!
     end subroutine siglevs_trop
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_slope (ndim,x,y,sig_slope)
!
! least square fit to find slope
!
! Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1   
     implicit none
!
     integer, intent(in) :: ndim
     real(rkind1), intent(in)  :: x(ndim)
     real(rkind1), intent(in)  :: y(ndim)
     real(rkind1), intent(out) :: sig_slope
!
     integer :: i
     real(rkind1) :: a,b,c,d
!
     a=0.
     b=0.
     c=0.
     d=0.
!
     do i=1,ndim
       a = a+x(i)
       b = b+y(i)
       c = c+x(i)*x(i)
       d = d+x(i)*y(i)
     enddo
!
     sig_slope=(ndim*d-a*b)/(ndim*c-a**2)
!
     end subroutine siglevs_slope
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_qc (ndim1,ndim2,nraw,nrpt,bmiss,ckind, &
                            snd_raw,snd_rpt)
!
!  Copy new obs to required arrays and set qc and program codes. 
!  Use QC flags in corresponding real obs to set QC in simulated obs.  
!  This will help ensure that the number and distribution of qc accepted 
!  obs is similar in the simulated and real obs.
!
!  Note, if the parameter set2miss is set to .true., then data values 
!  associated with 'bad' QC marks are set to missing values. If so, however, 
!  although reports will be written to the bufr file, when they are 
!  subsequently read, the reports will be flagged as having a missing data 
!  record. For single-level reports, there will then be 0 data levels noted. 
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1
     use m_conv_names, only : bbp,bbz,bbt,bbq,bbu,bbv
     use m_conv_names, only : bbpq,bbzq,bbtq,bbqq,bbwq
     use m_bufr_conv, only : conv_bmiss
!
     implicit none
!
     integer, intent(in) :: ndim1,ndim2
     integer, intent(in) :: nraw
     integer, intent(inout) :: nrpt
     real(rkind1), intent(in) :: bmiss
     real(rkind1), intent(in)    :: snd_raw(ndim1,ndim2)
     real(rkind1), intent(inout) :: snd_rpt(ndim1,ndim2)
     character(len=*), intent(in) :: ckind
!
     logical, parameter :: set2miss=.false. 
     integer :: n, m, m1
     real(rkind1), parameter :: wind_max=200.
     real(rkind1), parameter :: wind_min=-wind_max
     real(rkind1), parameter :: temp_max=340.
     real(rkind1), parameter :: sphu_max=0.10
     real(rkind1) :: diff,diffmin
     real(rkind1) :: uq  ! u quality
!
     do n=1,nrpt
!
! Find index of old obs whose p is closest to that of real obs
       diffmin=1.e10
       do m1=1,nraw
         diff=abs(snd_raw(bbp,m1)-snd_rpt(bbp,n))
         if (diff < diffmin) then
           diffmin=diff
           m=m1 
         endif
       enddo        
!
! Check quality of obs based on qualities of both new and old obs
! For u and v, set  qualities identically to largest quality flag value
       call siglevs_range_check (snd_raw(bbu,m),snd_raw(bbwq,n), &
                bmiss,wind_min,wind_max,snd_rpt(bbu,m),uq,set2miss)        
       call siglevs_range_check (snd_raw(bbv,m),snd_raw(bbwq,n), &
                bmiss,wind_min,wind_max,snd_rpt(bbv,m),snd_rpt(bbwq,n), &
                set2miss)
       snd_rpt(bbwq,n)=max(snd_rpt(bbwq,n),uq)        
       call siglevs_range_check (snd_raw(bbt,m),snd_raw(bbtq,n), &
                bmiss,0.,temp_max,snd_rpt(bbt,m),snd_rpt(bbtq,n),set2miss)      
       call siglevs_range_check (snd_raw(bbq,m),snd_raw(bbqq,n), &
                bmiss,0.,sphu_max,snd_rpt(bbq,m),snd_rpt(bbqq,n),set2miss)
!
! reset q quality values above 300 mb
       if (snd_rpt(bbp,n) <= 30000.) then
         snd_rpt(bbqq,n)=9.
       endif
!
! check z quality
       if (snd_raw(bbz,m) > 1.e5 .or. snd_raw(bbzq,m) > 3.5) then
         if (set2miss) then 
           snd_rpt(bbz,n)=bmiss
         endif
         snd_rpt(bbzq,n)=11.
       else
         snd_rpt(bbzq,n)=1.    ! z value OK
       endif
!
! check p quality
       if (snd_raw(bbp,m) > 1.e6 .or. snd_raw(bbpq,m) > 3.5) then
         if (set2miss) then 
           snd_rpt(bbp,n)=bmiss
         endif
         snd_rpt(bbpq,n)=11.
       else
         snd_rpt(bbpq,n)=1.    ! p value OK
       endif
!
     enddo  ! loop over set of new obs
!
! For sonde winds, remove first and last values if nrpt > 1
     if (ckind == 'WIND') then
       if (nrpt < 3) then
         m1=1
       else
         m1=2
       endif
       n=0
       do m=m1,max(1,nrpt-1)
         n=n+1
         snd_rpt(:,n)=snd_rpt(:,m)
       enddo
       nrpt=n
     endif
!
     end subroutine siglevs_qc         
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine siglevs_range_check (f_raw,q_raw,bmiss,fmin,fmax, &
                                     f_rpt,q_rpt,set2miss)
!
!  Assign quality based on whether value is in range and on original qc mark
!  A quality mark set to 10 indicates problem with sounding in NR.
!  A quality mark set to 11 indicates problem with corresponding real sounding
!  that has been copied over to the simulated sounding. 
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     use m_kinds, only : rkind1
     implicit none
!
     logical, intent(in) :: set2miss
     real(rkind1), intent(in) :: f_raw
     real(rkind1), intent(in) :: q_raw
     real(rkind1), intent(in) :: bmiss
     real(rkind1), intent(in) :: fmin
     real(rkind1), intent(in) :: fmax
     real(rkind1), intent(inout) :: f_rpt
     real(rkind1), intent(out) :: q_rpt
!
     if (f_rpt < fmin .or. f_rpt > fmax) then  
       q_rpt=10.
     elseif (f_raw < fmin .or. f_raw > fmax .or. q_raw > 3.5) then  
       q_rpt=11.
     else
       q_rpt=1.
     endif
!
!  If requested, set data values with bad qc mark to a missing value
     if (q_rpt /= 1. .and. set2miss) then 
       f_rpt=bmiss
     endif 
!
     end subroutine siglevs_range_check 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
! 
     subroutine summary_counts (icounts)
!     
!  For each type (NCEP kx) of observation, print a table of 
!  kx index
!  number of reports of this type read and accepted in input bufr file
!  number of obs values of this type read and accepted in input file file
!  number of simulated reports of this type written to output bufr file
!  number of simulated obs values of this type written to output bufr file
!
!  Note: the number of obs is simply determined as 2x number of obs levels,
!  and therefore does not properly account for ps values or T obs without q
!
!  Code History:
!
!    Ronald Errico       23 Sept 2014    Initial Code
!
     implicit none
!
     integer, intent(in) :: icounts(2,2,100:299)
!
     integer :: kx
!
     print *,' ' 
     print *,'Counts of obs by type:'
     print *,' kx     rpts_read      obs_read  rpts_written   obs_written'
     do kx=100,299
       if (icounts(1,1,kx) > 0) then 
         print ('(i3,4i14)'),kx,icounts(:,:,kx)
       endif
     enddo
!
     end subroutine summary_counts 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sonde_drift_location (itype,tmax,wlat,wlon,nru,nrv,t1,z1, &
                                    lat1,lon1,u1,v1,t2,z2,lat2,lon2)
!   
!  Determine new location (lat,lon,height) of sonde after short drift interval
!  Use simple first-order forward advection scheme with short time step
!
!  Code History:
!
!    Nikki Prive       23 Sept 2014    Initial Code
!
   use m_kinds, only : rkind1,rkind2
   use m_parameters, only : earthr, pifac_k2, rpifac_k2
!
   implicit none
!
   integer, intent(in) :: itype
   real(rkind1), intent(in) :: tmax
   real(rkind1), intent(in) :: t1
   real(rkind1), intent(in) :: z1
   real(rkind1), intent(in) :: lat1
   real(rkind1), intent(in) :: lon1
   real(rkind1), intent(in) :: u1
   real(rkind1), intent(in) :: v1
   real(rkind1), intent(in) :: wlat    ! grid latitude adjacent to pole
   real(rkind1),dimension(12), intent(in) :: wlon ! selected longitudes
   real(rkind1),dimension(12), intent(in) :: nru  ! u at select pts adj to pole
   real(rkind1),dimension(12), intent(in) :: nrv  ! v at select pts adj to pole
   real(rkind1), intent(out) :: t2
   real(rkind1), intent(out) :: z2
   real(rkind1), intent(out) :: lat2
   real(rkind1), intent(out) :: lon2
!
   real(rkind1), parameter :: zero=0._rkind1
   real(rkind1), parameter :: delt=30.    ! time between interp (seconds)
   real(rkind1), parameter :: ascent=5.   ! ascent rate for raobs/pibals in m/s
   real(rkind1), parameter :: descent=-5. ! descent rate dropsondes
   real(rkind1), parameter :: sec2hr=1./3600.  ! fac to change sec to hr.
   real(rkind1), parameter :: r90=90._rkind1
   real(rkind1), parameter :: r180=180._rkind1
   real(rkind1), parameter :: r360=360._rkind1
   real(rkind1), parameter :: pi=3.1415926
   real(rkind1) :: dt,dz,dy   ! time and vertical increment
   real(rkind2) :: nrlat,nrlon
   real(rkind2) :: latrad,lonrad   ! lat and lon in radians
   real(rkind1) :: th1             ! same as theat except rkind1
   real(rkind1) :: cx, cy, cux, cvy
   real(rkind2) :: ang, xf, yf, nxf, nyf, r
   real(rkind1), dimension(12) :: pu, pv, px, py
   integer :: np, npt, pt, i
   real, external :: interpolloc
!
! Determine new time of sounding   
   t2=min(t1+delt*sec2hr,tmax)
   if (tmax-t2 < 1.*sec2hr) then
     t2=tmax
   endif
   dt=(t2-t1)/sec2hr
!
! Determine new height of sounding
   if (mod(itype,100) == 32) then ! dropsonde
     dz=dt*descent
   else
     dz=dt*ascent
   endif
   z2=z1+dz
!
   if (abs(lat1) <= abs(wlat)) then  ! not within 1 grid pt of pole
!
! northward advection   
     lat2=lat1+rpifac_k2*dt*v1/earthr
     lat2=min(lat2,r90)
     lat2=max(lat2,-r90)
     latrad=lat2*pifac_k2
!
! eastward advection   
     lon2=lon1+rpifac_k2*dt*u1/(earthr*cos(pifac_k2*lat1))
     if (lon2 < 0.) then 
       lon2=lon2+r360 
     elseif (lon2 >= r360) then
       lon2=lon2-r360
     endif
!
   else  ! initial point near pole
!
     npt=0
     do pt=1,12 
       npt=npt+1
       nrlat=abs(wlat)
       r=earthr*((90-nrlat)*pifac_k2) 
       px(npt)=r*cos(pifac_k2*wlon(pt))
       py(npt)=r*sin(pifac_k2*wlon(pt))
       pu(npt)=nru(pt)*cos(pifac_k2*(wlon(pt)+90.))- &
               nrv(pt)*cos(pifac_k2*wlon(pt))*sign(1.,lat1)
       pv(npt)=nru(pt)*sin(pifac_k2*(wlon(pt)+90))- &
               nrv(pt)*sin(pifac_k2*wlon(pt))*sign(1.,lat1)
     end do
!
     cx=(earthr*((90-abs(lat1))*pifac_k2))*cos(pifac_k2*lon1)
     cy=(earthr*((90-abs(lat1))*pifac_k2))*sin(pifac_k2*lon1)
     cux=interpolloc(px,py,pu,cx,cy);
     cvy=interpolloc(px,py,pv,cx,cy);
     nxf=cx+(dt*cux)
     nyf=cy+(dt*cvy)
     lon2=rpifac_k2*atan2(nyf,nxf)
     if (lon2 < 0) then
       lon2 = lon2+360.
     elseif (lon2 >= 360.) then
       lon2 = lon2-360.
     endif
     if (lat1 > 0) then
       lat2=90.-(180.*(sqrt(nxf*nxf+nyf*nyf)/(earthr*pi)))
     else
       lat2=-(90.-(180.*(sqrt(nxf*nxf+nyf*nyf)/(earthr*pi))))
     end if
!
   endif
!
   end subroutine sonde_drift_location
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   real function interpolloc (x1,y1,dat,x2,y2)
!
!  Code History:
!
!    Nikki Prive       23 Sept 2014    Initial Code
!
!   
   use m_kinds, only : rkind1,rkind2
   use m_parameters, only : earthr, pifac_k2, rpifac_k2
!
   implicit none
!
   real(rkind1),dimension(12), intent(in) :: x1
   real(rkind1),dimension(12), intent(in) :: y1
   real(rkind1),dimension(12), intent(in) :: dat
   real(rkind1), intent(in) :: x2
   real(rkind1), intent(in) :: y2
!
   real(rkind1) :: npts    ! number of points used around pole
   real(rkind1) :: twt     ! weighting total
   real(rkind1) :: summ    ! sum
   integer :: skp, i       ! skip
   real(rkind1) :: d       ! distance
!
   npts = size(x1)
!
   twt = 0.
   summ = 0.
!
   do i=1,npts
     d = sqrt((x2-x1(i))**2 + (y2-y1(i))**2)
     twt = twt+ (1./d)
     summ = summ+(dat(i)/d) 
   enddo
!
   interpolloc = summ/twt
! 
   end function interpolloc
!


