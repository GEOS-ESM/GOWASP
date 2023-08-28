   module m_gpsro_ropp
!
!  Module to interface with ropp routines
!
   use m_kinds, only : rkind1, rkind2
   use typesizes, only: wp => eightbytereal
   use ropp_fm_types
   use ropp_fm_constants
!
   implicit none
!
   private
   public :: gpsro_ropp_setup
   public :: gpsro_ropp_clean
   public :: gpsro_ropp_sequence 
!
   integer :: ropp_nrun
   integer :: ropp_nlevs
   integer :: ropp_nx
   integer, public, parameter :: ropp_nobs1=1   
!
   real(rkind1) , parameter :: dx_plane=7.0e3      ! grid spacing (m) for plane
   real(rkind1) , parameter :: span_plane=500.0e3 ! span (m) of plane 
   real(wp) :: ropp_dx         ! grid spacing for plane in m (=dx_plane)
   real(wp) :: ropp_dtheta     ! grid spacing for plane in radians
!
   type(state2dFM), allocatable :: nr_plane(:)
   type(Obs1dBangle), allocatable :: obs_set1(:)
   character(len=*), parameter :: myname='m_gpsro_ropp'
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_setup (nrun)
!
! Set some variables and allocate memory used by ROPP
!
   use m_nr_fields_info, only : field_kdim
   use m_parameters, only  : earthr 
   implicit none
   integer, intent(in) :: nrun
!
   ropp_nrun=nrun
   ropp_nlevs=field_kdim
   ropp_dx=dx_plane
   ropp_nx=1+2*int(0.5*span_plane/dx_plane)  ! ensures number is odd
   ropp_dtheta=dx_plane/earthr
!
   allocate (nr_plane(ropp_nrun)) 
   allocate (obs_set1(ropp_nrun)) 
   call gpsro_ropp_alloc_strucs 
!
   end subroutine gpsro_ropp_setup
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_sequence (obs_info1,obs_values1,lprint, & 
                                   azim_flag,time1,time2,j1,ier)
!
!   Call sequence of routines to create propagation plane and bending angles
!
   use m_nr_fields_info, only : field_gpsro_algorithm, field_kmax
   use m_gpsro_names, only : obs_info_num, obs_max_levs, obs_nvalues
   use m_gpsro_names, only : bbnlold, bbdhr 
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(in) :: j1
   integer, intent(in) :: azim_flag ! flag to indicate 1 or 2 d calc
   integer, intent(out) :: ier
   real(rkind1), intent(in) :: time1,time2 ! braketing times for tslot
   real(rkind2), intent(in) :: obs_info1(obs_info_num)
   real(rkind1), intent(inout) :: obs_values1(obs_nvalues,obs_max_levs)
!   
   integer :: obs_nlevs
   real(rkind1) :: obs_time
!
! The following are used by the GSI GPSRO algorithm 
   real(rkind1) :: tges(field_kmax), qges(field_kmax), hges(field_kmax)
   real(rkind1) :: prsltmp(field_kmax+1), zsges  
   real(rkind1) :: tpdpres(1), rocprof, unprof, xlat, xlon, obshght
   real(rkind1) :: gsibend, hob
!
   ier=0
   obs_nlevs=obs_info1(bbnlold)
   obs_time=obs_info1(bbdhr)
   call gpsro_ropp_plane_grid (lprint,obs_info1,j1)
   call gpsro_ropp_plane_values (time1,obs_time,time2,j1)
   call gpsro_ropp_set_impp (obs_info1,obs_values1,obs_set1,j1)
   if (field_gpsro_algorithm < 2) then
     call ropp_fm_bangle_2d (azim_flag,nr_plane(j1),obs_set1(j1)) 
   else                                 ! use GSI algorithm for GPSRO
     call gpsro_ropp_gsiprof (tges,qges,hges,prsltmp,zsges, &
                              field_gpsro_algorithm,        &
                              xlat,xlon,tpdpres,rocprof,unprof,obshght,j1)
     call setupbend  (field_kmax,xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                      field_gpsro_algorithm,                             &
                      tges,qges,hges,prsltmp,gsibend,hob)
     obs_set1(j1)%bangle(1)=gsibend
   endif
   call gpsro_ropp_get_bbang (obs_values1,obs_set1,j1)
!
   end subroutine gpsro_ropp_sequence 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_plane_values (time1,obs_time,time2,j1)
!
!  Interpolate NR values to a GPS propagation plane.
!
   use m_nr_fields_info, only : field_imax, field_jmax, field_kdim
   use m_nr_fields_info, only : field_lon_first
   use m_nr_fields_info, only : field_max_akbk, field_akbk
   use m_nr_fields_info, only : field_gpsro_ksmooth
!
   use m_gpsro_fields, only : horiz_interp_2dfld
   use m_gpsro_fields, only : horiz_interp_3dfld
!
   use m_parameters, only : grav
  
!
   implicit none
!
   integer, intent(in) :: j1
   real(rkind1), intent(in) :: obs_time  
   real(rkind1), intent(in) :: time1,time2 ! braketing time for tslot
!
   integer :: ndd,ndi
   integer :: i,k,k1
   integer :: h_index(2,2)
   real(rkind1) :: plane_lat      ! lat of selected point in propagation plane
   real(rkind1) :: plane_lon      ! lon of selected point in propagation plane
   real(rkind1) :: h_weights(2,2)
   real(rkind1) :: t_weights(2)   ! weights used for temporal interp
   real(rkind1) :: phis2t(1,1)    ! (1) is horiz interpolated phis; (2) not used
   real(rkind1) :: ps2t(1,2)      ! horiz interpolated ps at 2 times
   real(rkind1) :: ps1t(1)        ! time and horiz interpolated ps
   real(rkind1) :: f2t(field_kdim,2)  ! horiz interpolated field at 2 times
                                  ! (only values for levs k1,...,k2 computed
   real(rkind1) :: ps,phis
   real(rkind1) :: t(field_kdim),q(field_kdim)  ! columns of t, q fields
   real(rkind1) :: pk(field_kdim+1,2)   ! p at (1) interface, (2) data levels 
   real(rkind1) :: phik(field_kdim+1,2) ! phi at (1) interface, (2) data levels 
!
   ndd=field_kdim    ! number of data levels
   ndi=field_kdim+1  ! number of data interfaces
   call get_interp_time_weights (time1,obs_time,time2,t_weights) 
!
   do i=1,ropp_nx ! loop over horizontal locations of 2-D plane
     plane_lat=nr_plane(j1)%lat(i)
     plane_lon=nr_plane(j1)%lon(i)
     call get_interp_horiz_index (field_imax,field_jmax,field_lon_first, &
                                  plane_lat,plane_lon,h_index,h_weights)  
!
! Get phis and ps
     call horiz_interp_2dfld (h_index,h_weights,'PHIS',1,phis2t)
     phis=phis2t(1,1)
     call horiz_interp_2dfld (h_index,h_weights,'PS',2,ps2t)
     call interp_time (1,1,1,t_weights,ps2t,ps1t)   
     ps=ps1t(1)
!
! Get T          
     call horiz_interp_3dfld (ndd,h_index,h_weights,'T',2,1,ndd,f2t)
     call interp_time (ndd,1,ndd,t_weights,f2t,t)
!
! Smooth T if NR T profile has "noise" that should be removed
     if (field_gpsro_ksmooth(3) > 0) then 
       call smooth_profile (t,ndd,field_gpsro_ksmooth)
     endif
!
! Get q
     call horiz_interp_3dfld (ndd,h_index,h_weights,'Q',2,1,ndd,f2t)
     call interp_time (ndd,1,ndd,t_weights,f2t,q)
!
! Compute pk and phis at both interfaces and data levels
     call compute_pk (field_max_akbk,ndi,field_akbk,ps,pk)
     call hydros_phis_tq (ndd,phis,pk,t,q,phik)  
!
     do k=1,ndd 
       k1=ndd+1-k                       ! reverses order to bottom to top
       nr_plane(j1)%temp(k,i)=t(k1)
       nr_plane(j1)%shum(k,i)=q(k1)
       nr_plane(j1)%pres(k,i)=pk(k1,2)        ! p at data levels
       nr_plane(j1)%geop(k,i)=phik(k1,2)/grav ! geop height at data levels
     enddo
     nr_plane(j1)%pres_sfc(i)=ps
     nr_plane(j1)%geop_sfc(i)=phis/grav       ! topographic height  
!    
    enddo
!       
   end subroutine gpsro_ropp_plane_values 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_plane_grid (lprint,obs_info1,j1)
!  
! Determine grid locations in a GPS propgation plane
!
   use m_gpsro_names, only : obs_info_num, bblat, bblon, bbbaz 
!
   implicit none
!
   logical, intent(in) :: lprint
   integer, intent(in) :: j1
   real(rkind2), intent(in) :: obs_info1(obs_info_num)
!
   integer :: ierror
   real(wp) :: obs_lat,obs_lon,obs_azim     
   real(wp) :: dtheta
   character(len=*), parameter :: mysub=myname//'::gpsro_ropp_plane_grid'
!
   obs_lat=obs_info1(bblat)
   obs_lon=obs_info1(bblon)
   nr_plane(j1)%dtheta=ropp_dtheta 
   obs_azim=obs_info1(bbbaz)
   if (obs_azim > 180.0_wp) then 
     obs_azim=obs_azim-360.0_wp
   endif
!
   call ropp_fm_2d_plane (obs_lat,obs_lon,obs_azim,nr_plane(j1)%dtheta, &
         nr_plane(j1)%n_horiz,nr_plane(j1)%lat,nr_plane(j1)%lon,ierror) 
   if (ierror /= 0) then  
     if (lprint) then 
       print *,'Errror in sub= ',mysub
       print *,'ierror=',ierror
     endif
   endif
!
   end subroutine gpsro_ropp_plane_grid 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_set_impp (obs_info1,obs_values1,obs_set,j1)
!
! Copy impact parameters to use real obs locations in vertical.   
!
   use m_gpsro_names, only : obs_info_num, obs_max_levs, obs_nvalues
   use m_gpsro_names, only : bblat, bblon, bbbaz, bbdhr, bbgeodu, bbcurve
   use m_gpsro_names, only : bbimpp
   use m_parameters, only  : earthr 
!
   implicit none
!
   integer, intent(in) :: j1
   real(rkind2), intent(in) :: obs_info1(obs_info_num)
   real(rkind1), intent(in) :: obs_values1(obs_nvalues,obs_max_levs)
   type(obs1dbangle), intent(inout) :: obs_set(ropp_nrun)
!
   integer :: i
!
   obs_set(j1)%lat=obs_info1(bblat)
   obs_set(j1)%lon=obs_info1(bblon)
   obs_set(j1)%time=obs_info1(bbdhr)
   obs_set(j1)%azimuth=obs_info1(bbbaz)
   obs_set(j1)%undulation=obs_info1(bbgeodu)
   obs_set(j1)%r_curve=obs_info1(bbcurve)
   obs_set(j1)%r_earth=earthr
!
   do i=1,obs_set(j1)%nobs
     obs_set(j1)%impact(i)=obs_values1(bbimpp,i)
   enddo
!
   end subroutine gpsro_ropp_set_impp
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_get_bbang (obs_values1,obs_set,j1)
!
! Copy bending angle from ROPP array to obs_values array
!
   use m_gpsro_names, only : obs_max_levs, obs_nvalues, bbbang
!
   implicit none
!
   integer, intent(in) :: j1
   real(rkind1), intent(out) :: obs_values1(obs_nvalues,obs_max_levs)
   type(obs1dbangle), intent(in) :: obs_set(ropp_nrun)
!
   integer :: i
!
   do i=1,obs_set(j1)%nobs
     obs_values1(bbbang,i)=obs_set(j1)%bangle(i)
   enddo
!
   end subroutine gpsro_ropp_get_bbang
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_gsiprof (tges,qges,hges,prsltmp,zsges, &
                           modify_alg,                          &
                           xlat,xlon,tpdpres,rocprof,unprof,obshght,j1)
!
!  Compute field profile for input to GSI GPSRO algorithm from ROPP 
!  profile arrays
!
   use m_parameters, only : ratio4R, grav
   use m_nr_fields_info, only : field_kmax, field_max_akbk, field_akbk
!
   implicit none      
!
   integer, intent(in) :: j1
   integer, intent(in) :: modify_alg
   real(rkind1) :: tges(field_kmax), qges(field_kmax), hges(field_kmax)
   real(rkind1) :: prsltmp(field_kmax+1), zsges  
   real(rkind1) :: tpdpres(1), rocprof, unprof, xlat, xlon, obshght
   real(rkind1) :: ps, phis
   real(rkind1) :: pk(field_kmax+1,2)   ! p at (1) interface, (2) data levels 
   real(rkind1) :: phik(field_kmax+1,2) ! phi at (1) interface, (2) data levels
   real(rkind1) :: rev_tges(field_kmax), rev_qges(field_kmax)
   integer :: k, kc, km
   integer :: re_nsig, re_nsig_p1
!
   re_nsig=field_kmax
   re_nsig_p1=re_nsig+1
   kc=(ropp_nx+1)/2
   zsges=nr_plane(j1)%geop_sfc(kc)
   ps=nr_plane(1)%pres_sfc(kc)
   phis=zsges*grav
   unprof=obs_set1(j1)%undulation
   rocprof=obs_set1(j1)%r_curve
   tpdpres(1)=obs_set1(j1)%impact(1)
   obshght=tpdpres(1)-rocprof  ! As defined for .ods files
   do k=1,re_nsig
     qges(k)=nr_plane(j1)%shum(k,kc)
     tges(k)=nr_plane(j1)%temp(k,kc)*(1.+ratio4R*qges(k)) ! change to tv
     xlat=nr_plane(j1)%lat(kc)
     xlon=nr_plane(j1)%lon(kc)
   enddo
!
! hdros needs fields in reverse order
   do k=1,re_nsig
     km=re_nsig-k+1
     rev_qges(k)=qges(km)
     rev_tges(k)=tges(km)/(1.+ratio4R*qges(km)) ! hydros needs t 
   enddo  
!
   phis=0.  ! gsi code assumes the heights are above sfc. 
   call compute_pk (field_max_akbk,re_nsig+1,field_akbk,ps,pk)
   call hydros_phis_tq (re_nsig,phis,pk,rev_tges,rev_qges,phik)      
   do k=1,re_nsig
     km=re_nsig-k+1
     prsltmp(k)=log(pk(km+1,1)/1000.)   ! change from Pa to cb
     hges(k)=phik(km+1,1)/grav
     if (modify_alg > 3) then       ! change to values on data levels instead of interfaces
       prsltmp(k)=log(pk(km,2)/1000.)   ! change from Pa to cb
       hges(k)=phik(km,2)/grav
     endif
   enddo
   prsltmp(re_nsig+1)=log(pk(1,1)/1000.)
!
   end subroutine gpsro_ropp_gsiprof
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine gpsro_ropp_alloc_strucs 
!
! Allocate some arrays input to ROPP
!
   implicit none
!
   integer :: i
!
   do i=1,ropp_nrun
!
! allocation of propagation plane
     nr_plane(i)%n_horiz=ropp_nx
     nr_plane(i)%n_lev=ropp_nlevs
     allocate(nr_plane(i)%lat(ropp_nx))
     allocate(nr_plane(i)%lon(ropp_nx))
     allocate(nr_plane(i)%geop_sfc(ropp_nx))
     allocate(nr_plane(i)%pres_sfc(ropp_nx))
     allocate(nr_plane(i)%geop(ropp_nlevs,ropp_nx))
     allocate(nr_plane(i)%pres(ropp_nlevs,ropp_nx))
     allocate(nr_plane(i)%shum(ropp_nlevs,ropp_nx))
     allocate(nr_plane(i)%temp(ropp_nlevs,ropp_nx))
     allocate(nr_plane(i)%nr(ropp_nlevs,ropp_nx))
     allocate(nr_plane(i)%refrac(ropp_nlevs,ropp_nx))
!
! bending angle and impact parameter values for set 1
     obs_set1(i)%nobs=ropp_nobs1
     allocate(obs_set1(i)%bangle(ropp_nobs1))
     allocate(obs_set1(i)%weights(ropp_nobs1))
     allocate(obs_set1(i)%impact(ropp_nobs1))
     allocate(obs_set1(i)%rtan(ropp_nobs1 ))
     allocate(obs_set1(i)%a_path(ropp_nobs1,2))
!
   enddo  
!
   end subroutine gpsro_ropp_alloc_strucs
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine gpsro_ropp_clean 
!
! Deallocate arrays allocated in gpsro_ropp_alloc_strucs
!
   implicit none
!
   integer :: i
!
   do i = 1,ropp_nrun
!
     deallocate (nr_plane(i)%lat,nr_plane(i)%lon,nr_plane(i)%geop_sfc)
     deallocate (nr_plane(i)%pres_sfc,nr_plane(i)%temp,nr_plane(i)%nr)
     deallocate (nr_plane(i)%geop,nr_plane(i)%pres,nr_plane(i)%shum)
     deallocate (nr_plane(i)%refrac)
     deallocate (nr_plane)
!
     deallocate (obs_set1(i)%bangle,obs_set1(i)%weights,obs_set1(i)%impact)
     deallocate (obs_set1(i)%rtan,obs_set1(i)%a_path)
     deallocate (obs_set1)
!
   enddo
!
   end subroutine gpsro_ropp_clean 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_gpsro_ropp

