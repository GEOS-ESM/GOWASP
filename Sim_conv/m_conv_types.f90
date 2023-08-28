   module m_conv_types
!
!  module for sequencing operations to create obs for each category of 
!  conventional observations, as indicated by the public routines herein
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1, rkind2
!
   use m_nr_fields_info, only : field_imax, field_jmax
   use m_nr_fields_info, only : field_lon_first
   use m_nr_fields_info, only : field_max_akbk, field_akbk
!
   use m_conv_names, only : conv_max_levs
   use m_conv_names, only : bbx,bby,bbr,bbp,bbz,bbt,bbq,bbu,bbv,bbc
   use m_conv_names, only : ityp,ixob,iyob,idhr,ielv,bbtmin,bbnlold
   use m_conv_names, only : bbnall,bbnlnew,bbpmin,bbtslot,bblints
!
   use m_shmem_fields, only : horiz_interp_2dfld
   use m_shmem_fields, only : horiz_interp_3dfld
!
   use m_bufr_conv, only : conv_bmiss
!
!
   implicit none
!
   private
   public :: sonde_drift_wind
   public :: sonde_drift_mass
   public :: multi_level_report
   public :: single_level_report
   public :: surface_level_report
!
   contains 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine surface_level_report (dim_info,dim_fields,ctype,time1, &
                              time2,obsinfo,obs_data_o,obs_data_n)
!
! Create Surface Obs: ps, 2m T,q,, or 10m wind
!
   implicit none
!
! arguments
   integer,intent(in) :: dim_info     ! dim of obsinfo array (incl. header)
   integer,intent(in) :: dim_fields   ! number of values in level data 
   real(rkind1), intent(in)     :: time1,time2 ! braketing time for tslot
   real(rkind2), intent(inout)  :: obsinfo(dim_info)
   real(rkind1), intent(in)     :: obs_data_o(dim_fields)
   real(rkind1), intent(out)    :: obs_data_n(dim_fields)
   character(len=*), intent(in) :: ctype
!
   integer :: itype        ! obs type (NCEP kx)
   integer :: h_index(2,2) ! grid point indexes used for horiz interp
   real(rkind1) :: lat1    ! obs latitude 
   real(rkind1) :: lon1    ! obs longitude 
   real(rkind1) :: rt1     ! obs time
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: zs2t(1,2)      ! (1) is horiz interpolated zs; (2) not used
   real(rkind1) :: zs             ! scalar version of zs2t(1,1)
   real(rkind1) :: ps2t(1,2)      ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: ps             ! scalar version of ps1t
   real(rkind1) :: fs1t(1)        ! time and horiz interp sfc wind component
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: us,vs,ts,qs    ! 10 m u,v, 2m t, qv
!
! Get con_info (including obs_header) 
   lat1=obsinfo(iyob)        ! obs latitude
   lon1=obsinfo(ixob)        ! obs longitude  
   rt1=obsinfo(idhr)         ! obs time relative to (central) symoptic time
   itype=nint(obsinfo(ityp)) ! type index (kx)
!
! Get indexes and weights for horizontal and temporal interpolation 
   call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                lat1,lon1,h_index,h_weights)
   call get_interp_time_weights (time1,rt1,time2,t_weights)
!
! Get zs and ps 
   call horiz_interp_2dfld (h_index,h_weights,'zs',2,1,zs2t)
   zs=zs2t(1,1)
   call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
   call interp_time (1,1,1,t_weights,ps2t,ps1t)   
   ps=ps1t(1)
!
   obs_data_n(bbp)=ps
   obs_data_n(bbz)=zs
   obs_data_n(bbc)=obs_data_o(bbc)  ! copy "CAT" code in report
   obs_data_n(bbx)=conv_bmiss  ! no "drift" info in report
   obs_data_n(bby)=conv_bmiss  ! no "drift" info in report
   obs_data_n(bbr)=conv_bmiss  ! no "drift" info in report
   obsinfo(ielv)=zs
   obsinfo(bbnlnew)=1._rkind2
!
   if (ctype == 'PWD') then 
!
! Get surface (10m) winds at obs location
     call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     us=fs1t(1)     ! surface value
     call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     vs=fs1t(1)     ! surface value
!
     obs_data_n(bbu)=us
     obs_data_n(bbv)=vs
!
   else
!  
! Get surface (2m) T, q obs location
     call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     ts=fs1t(1)     ! surface value
     call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     qs=fs1t(1)     ! surface value
!
     obs_data_n(bbt)=ts
     obs_data_n(bbq)=qs
!
   endif
!
   end subroutine surface_level_report
! 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine single_level_report (dim_info,dim_fields,kdim,ctype,lwind20m, &
                                   time1,time2,obsinfo,obs_data_o,obs_data_n)
!
! Create single-level obs (excluding some surface types).  Examples here are  
! AIRCRAFT or SHIP that report 20m wind.  If read obs pressure or height is 
! below the nature run surface, then assign it to the surface value 
! (10m wind, 2m T,q)
!
   implicit none
!
! arguments
   logical, intent(in) :: lwind20m     ! compute 20m wind
   integer, intent(in) :: dim_info     ! dim of obsinfo array (incl. header)
   integer, intent(in) :: dim_fields   ! number of values in level data 
   integer, intent(in) :: kdim         ! number of field levels considered
   real(rkind1), intent(in)    :: time1,time2 ! braketing time for tslot
   real(rkind2), intent(inout) :: obsinfo(dim_info)
   real(rkind1), intent(in)    :: obs_data_o(dim_fields)
   real(rkind1), intent(inout) :: obs_data_n(dim_fields)
   character(len=*), intent(in) :: ctype
!
   integer :: itype        ! obs type (NCEP kx)
   integer :: h_index(2,2) ! grid point indexes used for horiz interp
   integer :: ndl          ! number of field data levels
   integer :: ndi          ! number of data lev interfaces (=ndl+1)
   integer :: k1,k2        ! levels for interp between data levels
!
   real(rkind1) :: lat1    ! obs latitude 
   real(rkind1) :: lon1    ! obs longitude 
   real(rkind1) :: rt1     ! obs time
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: v_weights(2)   ! weights used for vertical interp
   real(rkind1) :: zs2t(1,2)      ! (1) is horiz interpolated zs; (2) not used
   real(rkind1) :: zs             ! scalar version of zs2t(1,1)
   real(rkind1) :: ps2t(1,2)      ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: ps             ! scalar version of ps1t
   real(rkind1) :: f2t(kdim+1,2)  ! horiz interpolated field at 2 times
                                  ! (only values for levs k1,...,k2 computed
   real(rkind1) :: f1t(kdim+1)    ! temporal interp of f2t and sfc value
   real(rkind1) :: fs1t(1)        ! time and horiz interp sfc wind component
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: pk(kdim+1,2)   ! p at (1) interface and (2) data levels 
   real(rkind1) :: zk(kdim+1)     ! z at interface levels
   real(rkind1) :: p1,z1,u1,v1,t1,q1  ! fields at obs levels
!
! Get con_info (including obs_header) 
   ndl=kdim
   ndi=ndl+1
   lat1=obsinfo(iyob)        ! obs latitude
   lon1=obsinfo(ixob)        ! obs longitude  
   rt1=obsinfo(idhr)         ! obs time relative to (central) symoptic time
   itype=nint(obsinfo(ityp)) ! type index (kx)
!
! Get indexes and weights for horizontal and temporal interpolation 
   call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                lat1,lon1,h_index,h_weights)
   call get_interp_time_weights (time1,rt1,time2,t_weights)
!
! First get p and z of report if this is 'PWD'
!
   if (ctype == 'PWD') then 
!
! Get zs and ps 
     call horiz_interp_2dfld (h_index,h_weights,'zs',2,1,zs2t)
     zs=zs2t(1,1)
     call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
     obsinfo(ielv)=zs
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
!
! Get z (first need w2t=rho at all levels) 
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f3',2,2,1,ndl,f2t)
     call interp_time (ndi,1,ndl,t_weights,f2t,f1t)
     call hydros_z (ndl,zs,pk(:,1),f1t,zk) 
!
! Determine corresponding p and z of obs, consistant with NR fields
     if (lwind20m) then  ! Get p and z corresponding to 20m wind
       z1=20.+zs         ! level of obs
       call get_interp_vert_index (ndi,'z',zk,z1,k1,k2,v_weights)
       call interp_hydros_p (ndi,k1,k2,pk(:,1),zk,z1,p1)
!
     else
!
! Get z corresponding to obs value of p if observed value present;
! otherwise get p corresponding to observed z if present
       p1=obs_data_o(bbp)
       z1=obs_data_o(bbz)
       if (p1 < 0.99*conv_bmiss) then 
         call get_interp_vert_index (ndi,'p',pk(:,1),p1,k1,k2,v_weights)
         call interp_hydros_z (ndi,k1,k2,pk(:,1),zk,p1,z1)
       elseif (z1 <0.99*conv_bmiss) then
         call get_interp_vert_index (ndi,'z',zk,z1,k1,k2,v_weights)
         call interp_hydros_p (ndi,k1,k2,pk(:,1),zk,z1,p1)
       endif
!
     endif                       ! test on whether a 20m obs wind
!
     obsinfo(bbpmin)=ps          ! save ps for use when making T,q values
     obs_data_n(bbp)=p1
     obs_data_n(bbz)=z1
     obs_data_n(bbc)=obs_data_o(bbc)  ! copy "CAT" code in report
     obs_data_n(bbx)=conv_bmiss  ! no "drift" info in report
     obs_data_n(bby)=conv_bmiss  ! no "drift" info in report
     obs_data_n(bbr)=conv_bmiss  ! no "drift" info in report
     obsinfo(bbnlnew)=1._rkind2  ! number of obs levels in new report
!
! Next, get u,v if this is a wind type
! Note that for interpolations within the bottom-half of the lowest data 
! layer, it is assumed that the the 10m wind provided by the NR is actually 
! at the surface.
!
     if (itype > 199) then ! need winds since this is a wind report
       call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
!
! Get u and v
       call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
       if (k2 == ndi) then  ! get surface winds since obs near surface
         call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
         f2t(ndi,1:2)=fs2t(1,1:2)
       else
         f2t(ndi,1:2)=0.    ! value actually not used  
       endif 
       call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
       call interp_vert (ndi,k1,k2,v_weights,f1t,u1)
!
       call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t) 
       if (k2 == ndi) then  ! get surface winds since obs near surface
         call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
         f2t(ndi,1:2)=fs2t(1,1:2)
       else
         f2t(ndi,1:2)=0.    ! value actually not used  
       endif
       call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
       call interp_vert (ndi,k1,k2,v_weights,f1t,v1)
!
       obs_data_n(bbu)=u1
       obs_data_n(bbv)=v1
!
     else   ! this is not a wind report
       obs_data_n(bbu)=conv_bmiss
       obs_data_n(bbv)=conv_bmiss
!
     endif  ! test on whether this is a wind report
!
! Consider ctype = 'TQ'
!
   elseif (itype < 200) then  ! get T, q values if mass report 
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     ps=obsinfo(bbpmin)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
!
! Get indexes and weights for vertical interpolation between data levels
     p1=obs_data_n(bbp)
     call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
!
! Get T
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
     call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,t1)
!
! Get q
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,q1)
!
     obs_data_n(bbt)=t1
     obs_data_n(bbq)=q1
!
   else ! ctype='TQ' but itype not mass report
     obs_data_n(bbt)=conv_bmiss 
     obs_data_n(bbq)=conv_bmiss
! 
   endif   ! check on ctype
!
   end subroutine single_level_report
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine multi_level_report (dim_info,dim_fields,kdim,ctype,time1, &
                                  time2,obsinfo,obs_data_o,obs_data_n)
!
! Create Multi-level obs for non-sonde instruments
! (e.g., RASS, PROFILER, VADWIND)
!
   implicit none
!
! arguments
   integer,intent(in) :: dim_info     ! dim of obsinfo array (incl. header)
   integer,intent(in) :: dim_fields   ! number of values in level data 
   integer,intent(in) :: kdim         ! number of field levels considered
   real(rkind1), intent(in)    :: time1,time2 ! braketing time for tslot
   real(rkind2), intent(inout) :: obsinfo(dim_info)
   real(rkind1), intent(in)    :: obs_data_o(dim_fields,conv_max_levs)
   real(rkind1), intent(inout) :: obs_data_n(dim_fields,conv_max_levs)
   character(len=*), intent(in) :: ctype
!
   integer :: itype        ! obs type (NCEP kx)
   integer :: h_index(2,2) ! grid point indexes used for horiz interp
   integer :: kold         ! obs level index in old report
   integer :: klevs        ! number of obs levels in report
   integer :: ndl          ! number of field data levels
   integer :: ndi          ! number of data lev interfaces (=ndl+1)
   integer :: k1,k2        ! levels for interp between data levels
   integer :: knew         ! level index in new report
!
   real(rkind1) :: lat1    ! obs latitude 
   real(rkind1) :: lon1    ! obs longitude 
   real(rkind1) :: rt1     ! obs time
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: v_weights(2)   ! weights used for vertical interp
   real(rkind1) :: zs2t(1,2)      ! (1) is horiz interpolated zs; (2) not used
   real(rkind1) :: zs             ! scalar version of zs2t(1,1)
   real(rkind1) :: ps2t(1,2)      ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: ps             ! scalar version of ps1t
   real(rkind1) :: f2t(kdim+1,2)  ! horiz interpolated field at 2 times
                                  ! (only values for levs k1,...,k2 computed
   real(rkind1) :: f1t(kdim+1)    ! temporal interp of f2t and sfc value
   real(rkind1) :: fs1t(1)        ! time and horiz interp sfc wind component
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: pk(kdim+1,2)   ! p at (1) interface and (2) data levels 
   real(rkind1) :: zk(kdim+1)     ! z at interface levels
   real(rkind1) :: p1,z1,u1,v1,t1,q1  ! fields at obs levels
!
! Get con_info (including obs_header) 
   ndl=kdim
   ndi=ndl+1
   lat1=obsinfo(iyob)           ! obs latitude
   lon1=obsinfo(ixob)           ! obs longitude  
   rt1=obsinfo(idhr)            ! obs time relative to (central) symoptic time
   itype=nint(obsinfo(ityp))    ! type index (kx)
   klevs=nint(obsinfo(bbnlold)) ! number of obs levels in original report
!
! Get indexes and weights for horizontal and temporal interpolation 
   call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                lat1,lon1,h_index,h_weights)
   call get_interp_time_weights (time1,rt1,time2,t_weights)
!
   if (ctype == 'PWD') then 
!
! Get zs and ps 
     call horiz_interp_2dfld (h_index,h_weights,'zs',2,1,zs2t)
     zs=zs2t(1,1)
     call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
     obsinfo(ielv)=zs
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
!
! Get z (first need w2t=rho at all levels) 
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f3',2,2,1,ndl,f2t)
     call interp_time (ndi,1,ndl,t_weights,f2t,f1t)
     call hydros_z (ndl,zs,pk(:,1),f1t,zk) 
!
     knew=0
     do kold=1,klevs                  

! Get z corresponding to obs value of p if observed value present;
! otherwise get p corresponding to observed z if present
       p1=obs_data_o(bbp,kold)
       z1=obs_data_o(bbz,kold)
       if (p1 < 0.99*conv_bmiss) then 
         call get_interp_vert_index (ndi,'p',pk(:,1),p1,k1,k2,v_weights)
         call interp_hydros_z (ndi,k1,k2,pk(:,1),zk,p1,z1)
       elseif (z1 <0.99*conv_bmiss) then
         call get_interp_vert_index (ndi,'z',zk,z1,k1,k2,v_weights)
         call interp_hydros_p (ndi,k1,k2,pk(:,1),zk,z1,p1)
       endif
!         
       if (k1 < ndi) then ! observation not below surface
         knew=knew+1
!
         obs_data_n(bbp,knew)=p1
         obs_data_n(bbz,knew)=z1
         obs_data_n(bbc,knew)=obs_data_o(bbc,kold) ! copy "CAT" code in report
         obs_data_n(bbx,knew)=conv_bmiss  ! no "drift" info in report
         obs_data_n(bby,knew)=conv_bmiss  ! no "drift" info in report
         obs_data_n(bbr,knew)=conv_bmiss  ! no "drift" info in report
!
! Get u and v if required 
         if (itype > 199) then  
           call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
!
           call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
           call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
           if (k2 == ndi) then  ! get surface winds since obs near surface
             call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
             f2t(ndi,1:2)=fs2t(1,1:2)
           else
             f2t(ndi,1:2)=0.    ! value actually not used  
           endif 
           call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
           call interp_vert (ndi,k1,k2,v_weights,f1t,u1)
!
           call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t) 
           if (k2 == ndi) then  ! get surface winds since obs near surface
             call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
             f2t(ndi,1:2)=fs2t(1,1:2)
           else
             f2t(ndi,1:2)=0.    ! value actually not used  
           endif
           call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
           call interp_vert (ndi,k1,k2,v_weights,f1t,v1)
!
           obs_data_n(bbu,knew)=u1
           obs_data_n(bbv,knew)=v1
!
         else  ! this is not a wind report
           obs_data_n(bbu,knew)=conv_bmiss
           obs_data_n(bbv,knew)=conv_bmiss
!
         endif ! check on obs type
       endif   ! check if obs above ground
     enddo     ! loop over obs levels
     obsinfo(bbnlnew)=real(knew) ! number of obs levels in new report
     obsinfo(bbpmin)=ps          ! save ps for use when making T,q values
!
!  Case for ctype = 'TQ'
!
   elseif (itype < 200) then     ! add T, q values
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     ps=obsinfo(bbpmin)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
!
     klevs=nint(obsinfo(bbnlnew))
     do knew=1,klevs
! 
! Get indexes and weights for vertical interpolation between data levels
       p1=obs_data_n(bbp,knew)
       call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
!
! Get T
       call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
       if (k2 == ndi) then  ! get surface winds since obs near surface
         call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
         f2t(ndi,1:2)=fs2t(1,1:2)
       else
         f2t(ndi,1:2)=0.    ! value actually not used  
       endif 
       call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
       call interp_vert (ndi,k1,k2,v_weights,f1t,t1)
!
! Get q
       call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t)
       if (k2 == ndi) then  ! get surface winds since obs near surface
         call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
         f2t(ndi,1:2)=fs2t(1,1:2)
       else
         f2t(ndi,1:2)=0.    ! value actually not used  
       endif 
       call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
       call interp_vert (ndi,k1,k2,v_weights,f1t,q1)
!
       obs_data_n(bbt,knew)=t1
       obs_data_n(bbq,knew)=q1
!
     enddo ! loop over obs levels for new T,q 
!
   else    ! this is not a T or q report
     klevs=nint(obsinfo(bbnlnew))
     do knew=1,klevs
       obs_data_n(bbt,knew)=conv_bmiss
       obs_data_n(bbq,knew)=conv_bmiss
     enddo
!
   endif   ! check on ctype
!
   end subroutine multi_level_report
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine sonde_drift_wind (dim_info,dim_fields,kdim,time1,time2, &
                                obsinfo,obs_data)
!
! Create Raob, Pibal, or Dropsonde observations of u,v, z, p 
! by ascent or descent through the NR atmosphere and considering 
! its wind-blown motion
!
   use m_nr_fields_info, only : field_pole_lat, field_pole_lons
!
   implicit none
!
! arguments
   integer,intent(in) :: dim_info     ! dim of obsinfo array (incl. header)
   integer,intent(in) :: dim_fields   ! number of values in level data 
   integer,intent(in) :: kdim         ! number of field levels considered
   real(rkind1), intent(in)    :: time1,time2 ! braketing time for tslot
   real(rkind2), intent(inout) :: obsinfo(dim_info)
   real(rkind1), intent(inout) :: obs_data(dim_fields,conv_max_levs)
!
   integer :: j           ! loop index = updated klast
   integer :: ndl         ! number of field data levels
   integer :: ndi         ! number of data lev interfaces (=ndl+1)
   integer :: k1,k2       ! levels for interp between data levels
   integer :: itype       ! obs type (NCEP kx)
   integer :: klast       ! number of obs levels previously considered; 
                          ! e.g., in previous time slots
   integer :: h_index(2,2)        ! grid point indexes used for horiz interp
!
   real(rkind1) :: lat1   ! 1st or previous latitude reported
   real(rkind1) :: lon1   ! 1st or previous longitude reported
   real(rkind1) :: lat2   ! new latiitude after sonde drift
   real(rkind1) :: lon2   ! new longitude after sonde drift
   real(rkind1) :: rt1    ! 1st or previous relative-time (hrs) reported
   real(rkind1) :: rt2    ! new obs time during sonde ascent or descent
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: v_weights(2)   ! weights used for vertical interp
   real(rkind1) :: pmin           ! p when balloon breaks or sonde dropped
   real(rkind1) :: zs2t(1,2)      ! (1) is horiz interpolated zs; (2) not used
   real(rkind1) :: zs             ! scalar version of zs2t(1,1)
   real(rkind1) :: ps2t(1,2)      ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: ps             ! scalar version of ps1t
   real(rkind1) :: f2t(kdim+1,2)  ! horiz interpolated field at 2 times
                                  ! (only values for levs k1,...,k2 computed
   real(rkind1) :: f1t(kdim+1)    ! temporal interp of f2t and sfc value
   real(rkind1) :: fs1t(1)        ! time and horiz interp sfc wind component
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: u1,v1,z1,p1    ! 1st or previous u,v,z,p reported
   real(rkind1) :: z2,p2          ! new sonde height after ascend or descend
   real(rkind1) :: pk(kdim+1,2)   ! p at (1) interface and (2) data levels 
   real(rkind1) :: zk(kdim+1)     ! z at interface levels
   real(rkind1) :: ptop           ! p
   real(rkind1) :: pole_lat       ! lat adjacent to pole
   real(rkind1) :: pole_u(12)     ! u at selected pts adjacent to pole
   real(rkind1) :: pole_v(12)     ! u at selected pts adjacent to pole
!
! Set the lowest pressure of ascent for balloons or the initial release 
! of dropsondes to lowest pressure in corresponding real obs, but not higher 
! than the p at the top of the subset of NR levels used.
   ptop=field_akbk(1,1)+field_akbk(1,2)*5.e4
   pmin=max(obsinfo(bbpmin),ptop)
!
! Get con_info (including obs_header) 
   ndl=kdim
   ndi=ndl+1
   lat1=obsinfo(iyob)       ! release latitude
   lon1=obsinfo(ixob)       ! release longitude  
   rt1 =obsinfo(bbtmin)     ! release time
   itype=nint(obsinfo(ityp))    ! type index (kx)
   klast=nint(obsinfo(bbnlnew)) ! previous obs level considered
   obsinfo(bblints)=0.d0     ! counter for number of obs in tslot
   pole_lat=sign(field_pole_lat,lat1) ! lat closest to pole insame hemisphere
!
! If initial time for balloon assumed launched at surface, then :
   if (klast == 0 .and. mod(itype,100) /= 32) then 
!
! Get indexes and weights for horizontal and temporal interpolation 
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                  lat1,lon1,h_index,h_weights)
     call get_interp_time_weights (time1,rt1,time2,t_weights)
!
! Get zs and ps (required to compute z and p at all levels)
     call horiz_interp_2dfld (h_index,h_weights,'zs',2,1,zs2t)
     zs=zs2t(1,1)
     call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
!
! Get surface (10m) winds at new location
     call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     u1=fs1t(1)     ! surface value
     call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
     call interp_time (1,1,1,t_weights,fs2t,fs1t)  
     v1=fs1t(1)     ! surface value
!
     klast=1
     obs_data(bbp,klast)=ps
     obs_data(bbz,klast)=zs
     obs_data(bbu,klast)=u1
     obs_data(bbv,klast)=v1
     obs_data(bbx,klast)=lon1
     obs_data(bby,klast)=lat1
     obs_data(bbr,klast)=rt1
     obsinfo(bblints)=1.d0
     obsinfo(ielv)=zs
!
   endif    ! test on whether first observed level for balloon ascent
!
!  Compute sounding for subsequent levels (including initial dropsonde)
!
   do j=klast+1,conv_max_levs
!     
     if (klast == 0) then         ! set initial dropsonde location
       lat2=lat1
       lon2=lon1
       rt2=rt1
       p2=pmin
     else                         !  Determine new sonde location after drift
       p1=obs_data(bbp,klast)
       z1=obs_data(bbz,klast)
       u1=obs_data(bbu,klast)
       v1=obs_data(bbv,klast)
       lon1=obs_data(bbx,klast)
       lat1=obs_data(bby,klast)
       rt1=obs_data(bbr,klast)
!
! Get more field values if advection near the pole is required
       if (abs(lat1) >= field_pole_lat) then ! need field desc. near pole
         if (j == klast+1) then   ! compute v_weights, t_weights, and k1,k2
           call get_interp_horiz_index (field_imax,field_jmax, &
                                  field_lon_first,             &
                                  lat1,lon1,h_index,h_weights)
           call get_interp_time_weights (time1,rt1,time2,t_weights)
           call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
           call interp_time (1,1,1,t_weights,ps2t,ps1t)   
           call compute_pk (field_max_akbk,ndi,field_akbk,ps1t(1),pk)
           call get_interp_vert_index (ndi,'p',pk(:,1),p1,k1,k2,v_weights)
         endif
!
         call polar_winds (ndi,k1,k2,v_weights,t_weights,pole_lat, &
                           pole_u,pole_v)
       endif
!
       call sonde_drift_location (itype,time2,pole_lat,field_pole_lons, &
                                  pole_u,pole_v,rt1,z1, &
                                  lat1,lon1,u1,v1,rt2,z2,lat2,lon2)
     endif
!
! Get indexes and weights for horizontal and temporal interpolation 
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                  lat2,lon2,h_index,h_weights)
     call get_interp_time_weights (time1,rt2,time2,t_weights)
!
! Get zs and ps at new location (required to compute z and p at all levels)
     call horiz_interp_2dfld (h_index,h_weights,'zs',2,1,zs2t)
     zs=zs2t(1,1)
     call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
!
! Get z (first need w2t=rho at all levels) 
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f3',2,2,1,ndl,f2t)
     call interp_time (ndi,1,ndl,t_weights,f2t,f1t)
     call hydros_z (ndl,zs,pk(:,1),f1t,zk)                   
!
     if (klast == 0) then ! find initial z2 corresponding to pmin for dropsonde 
       call get_interp_vert_index (ndi,'p',pk(:,1),p2,k1,k2,v_weights)
       call interp_hydros_z (ndi,k1,k2,pk(:,1),zk,p2,z2)
     else                 ! find p2 corresponding to hydrostatic z2 
       call get_interp_vert_index (ndi,'z',zk,z2,k1,k2,v_weights)
       call interp_hydros_p (ndi,k1,k2,pk(:,1),zk,z2,p2)
     endif
!
! Check if balloon is at top of its ascent (given by corres. real ascent)
     if (mod(itype,100) /= 32 .and. p2 <= pmin) then
       obsinfo(bbtslot)=999.d0
       exit     ! balloon breaks
     endif
! 
! Get indexes and weights for vertical interpolation between data levels
     call get_interp_vert_index (ndi,'p',pk(:,2),p2,k1,k2,v_weights)
!
! Get u and v
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,u1)
!
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t) 
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,v1)
!
     klast=klast+1

     obs_data(bbp,klast)=p2
     obs_data(bbz,klast)=z2
     obs_data(bbu,klast)=u1
     obs_data(bbv,klast)=v1
     obs_data(bbx,klast)=lon2
     obs_data(bby,klast)=lat2
     obs_data(bbr,klast)=rt2
     obsinfo(bblints)=obsinfo(bblints)+1.d0
!
! Check if dropsonde has hit the surface
     if (mod(itype,100) == 32 .and. p2 >= ps) then
       obsinfo(bbtslot)=999.d0
       exit   ! dropsonde at surface
     endif
!
! Check if sonde time has moved to next time slot
     if (rt2 >= time2) then
       obsinfo(bbtslot)=obsinfo(bbtslot)+1.d0
       exit
     endif    
!
   enddo  ! loop over levels
!
   obsinfo(bbnlnew)=real(klast)
   if (klast == conv_max_levs) then  ! no more levels to consider
     obsinfo(bbtslot)=999.d0
   endif
!
   end subroutine sonde_drift_wind 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine sonde_drift_mass (dim_info,dim_fields,kdim,time1,time2, & 
                                obsinfo,obs_data)
!
! Create Raob, Pibal, or Dropsonde observations of T or q at prescribed 
! locations, such as determined by sonde_drift_wind.
!
   use m_nr_fields_info, only : field_lev_first, field_conv_ksmooth

   implicit none
!
! arguments
   integer,intent(in) :: dim_info     ! dim of obsinfo array (incl. header)
   integer,intent(in) :: dim_fields   ! number of values in level data 
   integer,intent(in) :: kdim         ! number of field levels considered
   real(rkind1), intent(in)    :: time1,time2 ! braketing time for tslot
   real(rkind2), intent(inout) :: obsinfo(dim_info)
   real(rkind1), intent(inout) :: obs_data(dim_fields,conv_max_levs)
!
   integer :: ndl         ! number of field data levels
   integer :: ndi         ! number of data lev interfaces (=ndl+1)
   integer :: k1,k2       ! levels for interp between data levels
   integer :: k1s,k2s,k3s ! either =k1,k2,k2 or =1,kdim,kdim+1
   integer :: kinslot     ! number of obs levels in current tslot 
   integer :: klast       ! total number of obs levels computed thus far 
   integer :: kfirst      ! first obs level in current tslot
   integer :: klev
   integer :: h_index(2,2)        ! grid point indexes used for horiz interp
!
   real(rkind1) :: ps2t(1,2)        ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: ps             ! scalar version of ps1t
   real(rkind1) :: lat1   ! 1st or previous latitude reported
   real(rkind1) :: lon1   ! 1st or previous longitude reported
   real(rkind1) :: rt1    ! 1st or previous relative-time (hrs) reported
   real(rkind1) :: h_weights(2,2) ! weights used for horiz interp
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: v_weights(2)   ! weights used for vertical interp
   real(rkind1) :: f2t(kdim+1,2)  ! horiz interpolated field at 2 times
                                  ! (only for levs k1,...,k2)
   real(rkind1) :: f1t(kdim+1)    ! temporal and horiz interp fields
   real(rkind1) :: fs1t(1)        ! time and horiz interp sfc wind component
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: t1,q1,p1       ! t,q,p at current level
   real(rkind1) :: pk(kdim+1,2)   ! p at (1) interface and (2) data levels 
!
   ndl=kdim
   ndi=ndl+1
   klast=nint(obsinfo(bbnlnew))   ! last obs level considered for u,v data
   kinslot=nint(obsinfo(bblints)) ! number of obs layers in tslot
   kfirst=klast-kinslot+1
!
   do klev=kfirst,klast 
     p1=obs_data(bbp,klev)
     lon1=obs_data(bbx,klev)
     lat1=obs_data(bby,klev)
     rt1=obs_data(bbr,klev)
!
! Get indexes and weights for horizontal and temporal interpolation 
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                  lat1,lon1,h_index,h_weights)
     call get_interp_time_weights (time1,rt1,time2,t_weights)
!
! Get ps at new location (required to compute p at all levels)
     call horiz_interp_2dfld (h_index,h_weights,'ps',2,2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
!
! Get p at all interface and data levels (pk(ndi,1:2)=ps)
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
! 
! Get indexes and weights for vertical interpolation between data levels
     call get_interp_vert_index (ndi,'p',pk(:,2),p1,k1,k2,v_weights)
!
! Get T
! Smooth T if NR T profile has "noise" that should be removed as requested
     if (field_conv_ksmooth(3) > 0) then  ! whole column required
       k1s=1
       k2s=ndl
       k3s=ndi
     else                                 ! only 2 levels required
       k1s=k1
       k2s=k2
       k3s=k2
     endif 
!
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1s,k2s,f2t)
!
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
!
     call interp_time (ndi,k1s,k3s,t_weights,f2t,f1t)
!
     if (field_conv_ksmooth(3) > 0) then  ! smooth profile
       call smooth_profile (f1t,ndi,field_conv_ksmooth)
     endif
     call interp_vert (ndi,k1,k2,v_weights,f1t,t1)
!
! Get q
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,q1)
!
     obs_data(bbt,klev)=t1
     obs_data(bbq,klev)=q1
!
   enddo
!
   end subroutine sonde_drift_mass 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      
   subroutine polar_winds (ndi,k1,k2,v_weights,t_weights, &
                           pole_lat,pole_u,pole_v)
!
! Extract "sample" of u and v winds along grid latitude adjacent to pole
! Note that the same levels k1 ,k2 and v_weights are used for all points
!
   use m_nr_fields_info, only : field_pole_lons
!
   implicit none
!
   integer, intent(in) :: ndi
   integer, intent(in) :: k1,k2
   real(rkind1), intent(in)  :: v_weights(2)
   real(rkind1), intent(in)  :: t_weights(2)
   real(rkind1), intent(in)  :: pole_lat
   real(rkind1), intent(out) :: pole_u(12)
   real(rkind1), intent(out) :: pole_v(12)
!
   integer :: j
   integer :: h_index(2,2)        ! grid point indexes used for horiz interp
   real(rkind1) :: h_weights(2,2) 
   real(rkind1) :: fs2t(1,2)      ! horiz interp. sfc wind component at 2 times
   real(rkind1) :: f2t(ndi,2)     ! horiz interpolated field at 2 times
   real(rkind1) :: f1t(ndi)       ! temporal interp of f2t and sfc value
!
   do j=1,12
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                          pole_lat,field_pole_lons(j),h_index,h_weights)
!
! get u next to pole
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f1',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf1',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,pole_u(j))
!
! get v next to pole
     call horiz_interp_3dfld (ndi,h_index,h_weights,'f2',2,2,k1,k2,f2t)
     if (k2 == ndi) then  ! get surface winds since obs near surface
       call horiz_interp_2dfld (h_index,h_weights,'surf2',2,2,fs2t)
       f2t(ndi,1:2)=fs2t(1,1:2)
     else
       f2t(ndi,1:2)=0.    ! value actually not used  
     endif 
     call interp_time (ndi,k1,k2,t_weights,f2t,f1t)
     call interp_vert (ndi,k1,k2,v_weights,f1t,pole_v(j))
!
   enddo
!
   end subroutine polar_winds 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_conv_types
