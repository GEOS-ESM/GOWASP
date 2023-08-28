!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
!
   subroutine get_interp_horiz_index (imax,jmax,lon_first,obs_lat,obs_lon, &
                                      h_index,h_weights)
!
! Compute grid point indexes that surround a specific geographic location and 
! weights for bilinear interpolation using those points.
!
      use m_kinds, only : rkind1
      implicit none
! input
      integer, intent(in) :: imax
      integer, intent(in) :: jmax
      real(rkind1), intent(in) :: lon_first ! first longitude on each lat circle
      real(rkind1), intent(in) :: obs_lat
      real(rkind1), intent(in) :: obs_lon 
! output
      integer, intent(out) :: h_index(2,2)
      real(rkind1), intent(out) :: h_weights(2,2)
! local
      real(rkind1) :: hwmax
      real(rkind1) :: xlat    ! latitude in units of 180/(jmax-1) rel. to SP
      real(rkind1) :: xlon    ! longitude in units of 180/imax rel. to lon_first
      real(rkind1) :: rel_lon ! longitude relative to lon_first
      real(rkind1) :: lat_weights(2)
      real(rkind1) :: lon_weights(2)
!
      xlat=real(jmax-1,rkind1)*(obs_lat+90._rkind1)/180._rkind1
      h_index(1,2)=min(int(xlat)+1,jmax-1)   ! j value of closest point S
      h_index(2,2)=h_index(1,2)+1            ! j value of closest point N
      lat_weights(1)=real(h_index(2,2)-1,rkind1)-xlat      
      lat_weights(2)=1._rkind1-lat_weights(1)
!      
      rel_lon=mod(obs_lon-lon_first,360._rkind1)
      if (rel_lon < 0._rkind1) rel_lon=rel_lon+360._rkind1
      if (rel_lon >= 360._rkind1) rel_lon=0._rkind1  ! accounts for round off
      xlon=real(imax,rkind1)*rel_lon/360._rkind1
!
      h_index(1,1)=min(int(xlon)+1,imax)      ! i value closest point W
      if (h_index(1,1) == imax) then
        h_index(2,1)=1                        ! i value closest point E
      else
        h_index(2,1)=h_index(1,1)+1           ! i value closest point W
      endif
      lon_weights(2)=xlon-real(h_index(1,1)-1,rkind1)
      lon_weights(1)=1._rkind1-lon_weights(2)
!
      h_weights(1,1)=lon_weights(1)*lat_weights(1)
      h_weights(1,2)=lon_weights(1)*lat_weights(2)
      h_weights(2,1)=lon_weights(2)*lat_weights(1)
      h_weights(2,2)=lon_weights(2)*lat_weights(2)
!
   end subroutine get_interp_horiz_index 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_horiz_nearest (h_weights)
!
! re-define interpolation weights so that the nearest point is given 
! weight 1 and others re-set to 0.
!
      use m_kinds, only : rkind1
      implicit none
!
      real(rkind1), intent(inout) :: h_weights(2,2)
      integer :: loc_max(2)
!
      loc_max=maxloc(h_weights)
      h_weights(:,:)=0._rkind1
      h_weights(loc_max(1),loc_max(2))=1._rkind1
!
   end subroutine get_interp_horiz_nearest
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_time_weights (t1,t,t2,time_weights)
!
! determin wieights for linear interpolation in time
! t1 <= t <= t2
!
      use m_kinds, only : rkind1
      implicit none
!
      real(rkind1), intent(in) :: t1,t,t2  
      real(rkind1), intent(out) :: time_weights(2)
!
      time_weights(1)=(t2-t)/(t2-t1)
      time_weights(2)=1.-time_weights(1)
!
   end subroutine get_interp_time_weights
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_time_nearest (time_weights)
!
! re-define interpolation weights so that the nearest point is given 
! weight 1 and others re-set to 0.
!
      use m_kinds, only : rkind1
      implicit none
!
      real(rkind1), intent(inout) :: time_weights(2)
!
      if (time_weights(1) >= time_weights(2)) then
        time_weights(1)=1._rkind1  
        time_weights(2)=0._rkind1
      else  
        time_weights(2)=1._rkind1  
        time_weights(1)=0._rkind1
      endif 
!
   end subroutine get_interp_time_nearest
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interp_time (nk,k1,k2,time_weights,f2,f1)
!
! Perform interpolation in time
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk,k1,k2
      real(rkind1), intent(in) :: time_weights(2)
      real(rkind1), intent(in) :: f2(nk,2)
      real(rkind1), intent(out) :: f1(nk)
!
      integer :: k
! 
      do k=k1,k2
        f1(k)=time_weights(1)*f2(k,1)+time_weights(2)*f2(k,2)
      enddo
!
   end subroutine interp_time
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_interp_vert_index (nk,cf,fk,f,k1,k2,vweights)
!

! if cf='z', find pair of levels k1,k2 such that f(k1) > f > f(k2)
!    and then determine weights linear in z
! if cf='p', find pair of levels k1,k2 such that f(k1) < f < f(k2)
!    and then determine weights linear in log z
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk
      integer, intent(out) :: k1,k2
      real(rkind1), intent(in) :: fk(nk)
      real(rkind1), intent(in) :: f
      real(rkind1), intent(out) :: vweights(2)
      character(len=*), intent(in) :: cf
!
      integer :: k
! 
      if (cf == 'z') then
        k2=0
        do k=1,nk
          if (fk(k) < f) then
            k2=k
            exit
          endif
        enddo
        if (k2 == 1) then      ! value of f above top level
          k1=1
          vweights(1)=1.
        elseif (k2 == 0) then  ! value of f below bottom level
          k2=nk
          k1=nk
          vweights(1)=1.
        else                   ! value of f between levels
          k1=k2-1        
          vweights(1)=(f-fk(k2))/(fk(k1)-fk(k2)) 
        endif
      endif
!
      if (cf == 'p') then
        k2=0
        do k=1,nk
          if (fk(k) > f) then
            k2=k
            exit
          endif
        enddo
        if (k2 == 1) then      ! value of f above top level
          k1=1
          vweights(1)=1.
        elseif (k2 == 0) then  ! value of f below bottom level
          k2=nk
          k1=nk
          vweights(1)=1.
        else                   ! value of f between levels
          k1=k2-1         
          vweights(1)=log(fk(k2)/f)/log(fk(k2)/fk(k1))
        endif
      endif
!
      vweights(2)=1.-vweights(1)      

!
   end subroutine get_interp_vert_index
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interp_vert (nk,k1,k2,vweights,fk,f)
!
! Interpolate vertically between fk(k1) and fk(k2)
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk,k1,k2
      real(rkind1), intent(in)  :: vweights(2)
      real(rkind1), intent(in)  :: fk(nk)
      real(rkind1), intent(out) :: f
!
      f=vweights(1)*fk(k1)+vweights(2)*fk(k2)
!
   end subroutine interp_vert
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine hydros_phis_tq (nk,phis,p,t,q,phi)
!
! Compute hydrostatic phi at interfaces and data levels from p,t,q
!
      use m_kinds, only : rkind1
      use m_parameters, only : rdryair, ratio4R
      implicit none
!
      integer, intent(in) :: nk
      real(rkind1), intent(in)  :: phis
      real(rkind1), intent(in)  :: p(nk+1,2)
      real(rkind1), intent(in)  :: t(nk)
      real(rkind1), intent(in)  :: q(nk)
      real(rkind1), intent(out) :: phi(nk+1,2)
!
      integer :: k
      real(rkind1) :: rtv   ! R*TV
      real(rkind1) :: pfac  ! delta(log(p))
!
! First compute phi at interfaces
      phi(nk+1,1)=phis
      do k=nk,1,-1
        rtv=rdryair*t(k)*(1+ratio4R*q(k))
        pfac=log(p(k+1,1)/p(k,1))
        phi(k,1)=phi(k+1,1)+rtv*pfac
      enddo
!
! Next compute phi at data levels 
      phi(nk+1,2)=phis   ! this is not strictly a data level
      do k=1,nk
        rtv=rdryair*t(k)*(1+ratio4R*q(k))
        pfac=log(p(k+1,1)/p(k,2))
        phi(k,2)=phi(k+1,1)+rtv*pfac
      enddo     
!        
   end subroutine hydros_phis_tq
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine hydros_z (nk,zs,p_if,rho,z_if)
!
! Compute hydrostatic z at interfaces based using GEOS-5 layer-mean 
! density (defined for this purpose). 
!
      use m_kinds, only : rkind1
      use m_parameters, only : grav
      implicit none
!
      integer, intent(in) :: nk
      real(rkind1), intent(in)  :: zs
      real(rkind1), intent(in)  :: p_if(nk+1)
      real(rkind1), intent(in)  :: rho(nk)
      real(rkind1), intent(out) :: z_if(nk+1)
!
      integer :: k
!     
      z_if(nk+1)=zs
      do k=nk,1,-1
        z_if(k)=z_if(k+1)+(p_if(k+1)-p_if(k))/(grav*rho(k))
      enddo
!
   end subroutine hydros_z
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interp_hydros_z (nk,k1,k2,p_if,z_if,p,z)
!
! Interpolate z between data interfaces assuming constant virtual T 
! within a layer.
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk,k1,k2
      real(rkind1), intent(in)  :: p
      real(rkind1), intent(in)  :: p_if(nk)
      real(rkind1), intent(in)  :: z_if(nk)
      real(rkind1), intent(out) :: z
!
      real(rkind1) :: fac
!     
      if (k2 == 1) then
        fac=(z_if(1)-z_if(2))/log(p_if(2)/p_if(1))
      elseif (k1 == nk) then
        fac=0.
      else
        fac=(z_if(k1)-z_if(k2))/log(p_if(k2)/p_if(k1))
      endif
      z=z_if(k2)+fac*log(p_if(k2)/p)
!
   end subroutine interp_hydros_z
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interp_hydros_p (nk,k1,k2,p_if,z_if,z,p)
!
! Interpolate p between data interfaces assuming constant Tv
! within layer.
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk,k1,k2
      real(rkind1), intent(in)  :: z
      real(rkind1), intent(in)  :: p_if(nk)
      real(rkind1), intent(in)  :: z_if(nk)
      real(rkind1), intent(out) :: p
!
      real(rkind1) :: fac   ! = density/grav
      real(rkind1) :: ez    ! = exp(g*dz/density)
!     
      if (k2 == 1) then
        fac=(z_if(1)-z_if(2))/log(p_if(2)/p_if(1))
        ez=exp((z-z_if(2))/fac)  ! will be >= 1.
        p=p_if(2)/ez
      elseif (k1 == nk) then
        p=p_if(k1)  ! surface value
      else
        fac=(z_if(k1)-z_if(k2))/log(p_if(k2)/p_if(k1))
        ez=exp((z-z_if(k2))/fac)  ! will be >= 1.
        p=p_if(k2)/ez
      endif
!
   end subroutine interp_hydros_p
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine compute_pk (ndim,nk,akbk,ps,pk)
!
! Compute p at vertical grid interfaces (pk(:,1)) and grid data levels
! (pk(:,2)) using ps and the eta-coordinate parameters a(k), b(k) for 
! the vertical grid
!
! The formula for defining p at data levels is  based on constraining 
! the integrals of R*T*logp through a layer to be identical if 
! T=theta*(p/po)**(R/Cp) with theta constant through the layer or if, 
! instead, T is constant through the layer, with the relationship between 
! the constants theta and T defined using the data level values of p 
! determined here. 
!
      use m_kinds, only : rkind1
      use m_parameters, only : Xkappa, XkappaI
!
      implicit none
!
      integer, intent(in) :: ndim
      integer, intent(in) :: nk  
      real(rkind1), intent(in)  :: akbk(ndim,2) ! ak and bk values
      real(rkind1), intent(in)  :: ps
      real(rkind1), intent(out) :: pk(nk,2)
!
      integer :: k
      real(rkind1) :: logp(nk), pkap(nk), dlog, dpkap
!
! First, determine p at interfaces of vertical grid layers using definition 
! of eta coordinate
      do k=1,nk
        pk(k,1)=akbk(k,1)+akbk(k,2)*ps
      enddo
!
! Determine p on nk-1 data levels.  An additional value defined with a level 
! index nk is set as ps.
      pk(nk,2)=ps     
!
      do k=1,nk
        logp(k)=log(pk(k,1))
        pkap(k)=pk(k,1)**Xkappa
      enddo
!
      do k=1,nk-1
        dlog=logp(k+1)-logp(k)    
        dpkap=pkap(k+1)-pkap(k)
        pk(k,2)=(XkappaI*dpkap/dlog)**XkappaI
      enddo  
!
   end subroutine compute_pk
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine interp_rh (nk,k1,k2,vweights,pk,tk,qk,p1,t1,q1)
!
! Interpolate q vertically assumining linear in rh
!
      use m_kinds, only : rkind1
      implicit none
!
      integer, intent(in) :: nk,k1,k2
      real(rkind1), intent(in)  :: vweights(2)
      real(rkind1), intent(in)  :: pk(nk), tk(nk), qk(nk)
      real(rkind1), intent(in)  :: p1, t1 ! p and t at interpolated level
      real(rkind1), intent(out) :: q1     ! q at interpolated level
!
      real(rkind1) :: rh1, rh2, rhi
      call transform_q_rh (pk(k1),tk(k1),qk(k1),rh1,'q2rh')
      call transform_q_rh (pk(k2),tk(k2),qk(k2),rh2,'q2rh')
      rhi=vweights(1)*rh1+vweights(2)*rh2
      call transform_q_rh (p1,t1,q1,rhi,'rh2q')
!
   end subroutine interp_rh
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine transform_q_rh (p,t,q,rh,cf)
!
! Transform between specific and relative humidity
!
      use m_kinds, only : rkind1
      use m_parameters, only : ratio4R 
      implicit none
!
      real(rkind1), intent(in)      :: p,t
      real(rkind1), intent(inout)   :: q,rh
      character(len=*), intent(in) :: cf
!
      real(rkind1), parameter :: svpt0=273.15, svp1=611.2
      real(rkind1), parameter :: svp2=17.67, svp3=29.65
      real(rkind1) :: a,svp,qs
!
      a=svp2*(t-svpt0)/(t-svp3)
      svp=svp1*exp(a)             ! saturation vapor pressure
      qs=svp*ratio4R/(p-svp)      ! saturation specific humidity
      if (cf == 'q2rh') then 
        rh=q/qs
      else
        q=rh*qs
      endif
!
   end subroutine transform_q_rh 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine integrate_area_params (imax,jmax,ncdim,lat1,lon1,lon_first, &
                                     dx,nc,ij_index,w)
!
! Determine parameters for selecting a sample of Cartesian grid points that 
! fall within some distance (dx/2) of a specified location (lat1,lon1). For 
! longitudes near the poles, where meridional convergence has caused 
! longitudinal separations to be very small, consider only a sample of points 
! that are separated by approximately half the latitudinal separation. 
! Grid point values will be weighted by the area associated with each grid 
! point included, without considering its distance from the location.    
! If <4 points fall within that distance, use the 4 points that surround 
! the location and use as weights those for bi-linear interpolation. 
!
   use m_kinds, only : rkind1
   use m_parameters, only : earth_c, pifac_k1, earthr, pifac_k2
!
   implicit none
!
   integer, intent(in) :: imax,jmax    ! dimensions of grid (including pole points)
   integer, intent(in) :: ncdim        ! max permitted value of nc
   integer, intent(out) :: ij_index(2,ncdim) ! i,j indeces of points to use in avg
   integer, intent(out) :: nc   ! number of points to average over
   real(rkind1), intent(in) :: lat1,lon1 ! central point
   real(rkind1), intent(in) :: lon_first ! lon corresponding to index i=1
   real(rkind1), intent(in) :: dx        ! width of area to average over
   real(rkind1), intent(out) :: w(2,ncdim) ! weights on points in sample
!
   integer, parameter :: r8=8            ! need high precision to calc d
   integer :: i_in_box     
   integer :: j_in_box
   integer :: imax2
   integer :: j,jx,j1,j2
   integer :: i,ix,i1,i2,iz,iskip
   integer :: h_index(2,2),loc_max(2)
   real(r8) :: half,dxhalf8,earthr2
   real(r8) :: d1,d2,d3,d4,d
   real(rkind1) :: rlon1,rlat1
   real(rkind1) :: dlat_deg,dlat_km,dlat_half
   real(rkind1) :: dlon,dlon_eq,dlon_deg,dxhalf
   real(rkind1) :: lat2,lon2,rlat2,rlon2,lonx
   real(rkind1) :: lat2mh,lat2ph,wlat,wsum
   real(rkind1) :: h_weights(2,2)        !
!
   imax2=imax/2
   dxhalf=0.5*dx
   dlat_deg=180./real(jmax-1)       ! separation between grid lats (degrees)
   dlat_km=0.5*earth_c/real(jmax-1) ! separation between grid lats (km)
   dlat_half=0.5*dlat_km
   dlon_deg=360./real(imax)         ! separation between grid lons (degrees)
   dlon_eq=earth_c/real(imax)       ! sep.  between grid lons at equator (km)
   rlat1=pifac_k2*lat1              ! central point lat (radians)
   rlon1=pifac_k2*lon1              ! central point lon (radians)
   half=0.5_r8
   dxhalf8=real(dxhalf,r8)
   earthr2=real(earthr*.002,r8)     ! 2* earth radius (changed to km)
!
! Determine the 4 grid points surrounding the provided location
   call get_interp_horiz_index (imax,jmax,lon_first,lat1,lon1, &
                                h_index,h_weights)
!
! Determine max number of lat points to consider (min=2)   
   j_in_box=2+2*int(dxhalf/dlat_km)   ! number of points to consider
   j1=h_index(1,2)+1-j_in_box/2       ! S-most index to consider
   j2=j1+j_in_box-1                   ! N-most index to consider
!
! Loop over possible lats  
   nc=0  ! set counter to 0
   do jx=j1,j2
!
! Adjust for overlapping of poles
     if (jx < 1) then
       j=2-jx
     elseif (jx > jmax) then
       j=2*jmax-jx
     else
       j=jx
     endif
!
     lat2=dlat_deg*(j-1)-90.   ! lat of grid point (deg)
     rlat2=pifac_k2*lat2       ! lat of grid point (rad)
     lat2mh=max(dlat_deg*(j-1.5)-90.,-90.) ! 1/2 grid distance to S (deg)
     lat2ph=min(dlat_deg*(j-0.5)-90., 90.) ! 1/2 grid distance to N (deg)
     wlat=sin(pifac_k1*lat2ph)-sin(pifac_k1*lat2mh)  ! weight of this lat
!
     if (j > 1 .and. j < jmax) then
       dlon=dlon_eq*cos(rlat2)       ! distance between lon points (km) 
       iskip=1+int(dlat_half/dlon)   ! consider lon points spaced dlat/2 apart      
       i_in_box=2+2*int(dxhalf/dlon) ! number of points to consider if no skip
     else                            ! at pole point, consider 1 value only
       dlon=0.
       iskip=2
       i_in_box=1
     endif
!
     i1=h_index(1,1)+1-i_in_box/2     
     i2=i1+i_in_box-1   
     do ix=i1,i2,iskip
!
! Specify grid index i as that for the point across the pole when lat set
! crosses the pole
       if (jx > jmax .or. jx < 1) then
         if (ix > imax2) then
           iz=ix-imax2
         else
           iz=imax2+ix
         endif
       else
         iz=ix
       endif
!
! re-assign if wrap around
       if (iz < 1) then
         i=imax+iz
       elseif (iz > imax) then
         i=iz-imax 
       else
         i=iz
       endif
!
! Determine lon of grid point
       lonx=(i-1)*dlon_deg+lon_first
       lon2=mod(lonx,360.)
       if (lon2 < 0.) then
         lon2=360.-lon2
       endif
       rlon2=pifac_k2*lon2
! 
! Determine great-circle distance between grid and central points
       d1=sin(half*(rlat1-rlat2))
       d2=sin(half*(rlon1-rlon2))
       d3=cos(rlat1)*cos(rlat2)
       d4=d1*d1+d2*d2*d3
       if (d4 > 0.) then 
         d=earthr2*asin(sqrt(d4))
       else
         d=0._r8
       endif
!
! If distance with 1/2 dx of central point and the number of such points
! has not been exceeded, save the point indices and weight.
       if (d < dxhalf8) then
         nc=min(nc+1,ncdim)
         ij_index(1,nc)=i
         ij_index(2,nc)=j
         w(1,nc)=wlat
         w(2,nc)=d
         wsum=wsum+wlat
       endif
!
     enddo
   enddo 
!
! If <4 points found within the distance 1/2 dx, then use 4 
! surrounding points with weights those of bi-linear interpolation
   if (nc < 4) then 
     nc=4
     ij_index(1,1)=h_index(1,1)
     ij_index(1,2)=h_index(1,1)
     ij_index(1,3)=h_index(2,1)
     ij_index(1,4)=h_index(2,1)
!
     ij_index(2,1)=h_index(1,2)
     ij_index(2,2)=h_index(2,2)
     ij_index(2,3)=h_index(1,2)
     ij_index(2,4)=h_index(2,2)
!
     w(1,1)=h_weights(1,1)
     w(1,2)=h_weights(1,2)
     w(1,3)=h_weights(2,1)
     w(1,4)=h_weights(2,2)
     wsum=1. ! since already normalized h_weights
     w(2,1:4)=0.  ! not actual distances here
!
   endif
!
! Normalize weights such that sum=1
   w(1,1:nc)=w(1,1:nc)/wsum
!
   end subroutine integrate_area_params
