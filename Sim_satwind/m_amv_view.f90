!
   module m_amv_view
!
!  Determine if an obsevation is present based on lat, lon, 
!  plev, cloud, wv, wind speed, and probability function
!
   use m_kinds, only : rkind1,rkind2,zero_k1,one_k1
!
   use m_kx_table, only : kx_num,kx_jbins,kx_cbins,kx_pbins
   use m_kx_table, only : kx_nfilters,kx_nparams
   use m_kx_table, only : kx_present,kx_pt_count,kx_obs_count
   use m_kx_table, only : kx_histogram,kx_params
   use m_kx_table, only : kx_i_stride,kx_j_stride,kx_dx
   use m_kx_table, only : kx_type,kx_locs,kx_filters
   use m_kx_table, only : kx_filt_flat, kx_filt_slev, kx_filt_torb
   use m_kx_table, only : kx_speed_min, kx_cloud_obscure
   use m_kx_table, only : kx_ijgrid, kx_qcfac, kx_range
!
   implicit none
   private
!
   public amv_check_loc
   public amv_stride_i
   public amv_stride_j
   public amv_normalize_counts
   public amv_layers
   public amv_tune_params
   public amv_use_params
   public amv_accum_bins
   public amv_speed_cld
   public amv_speed_ipw
   public amv_obscured
!
   contains            
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_latlon_range (k,lat,lon,lrange)
!
! Check lat/lon range
! For polar orbiter, check that lats poleward of kx_range values
!
   implicit none
!
   integer, intent(in)  :: k
   real(rkind1), intent(in) :: lat,lon
   logical, intent(out) :: lrange
!
   real(rkind1) :: x
!
   lrange=.false.  ! default value indicats loc is out of viewing range
!
   if (kx_type(k)(1:1) == 'P') then       ! polar orbiter
     if (lat .le. kx_range(1,2,k) .or. lat .ge. kx_range(2,2,k)) then
       lrange=.true.
     endif
!
   else                                   ! geostationary    
     if (lat .ge. kx_range(1,2,k) .and. lat .le. kx_range(2,2,k)) then 
!
! Lat within range so check longitude relative to nadir
! (first adjust for case when separation crosses the prime meridian)
       x=lon-kx_locs(k)
       if (x < -180.) then
         x=x+360.
       elseif (x > 180.) then
         x=x-360.
       endif
!
       if (x < 0.) then  ! separation is westward from nadir  
         if (x .ge. kx_range(1,1,k)) then
           lrange=.true.
         endif
       else              ! separation is eastward from nadir
         if (x .le. kx_range(2,1,k)) then  
           lrange=.true.
         endif
       endif
     endif   ! check on lat for geost. view
   endif     ! check on whether polar or geost. view
!
   end subroutine amv_latlon_range 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_great_circle (lat1,lon1,lat2,lon2,dist)
!
! Calculate great circle distance between 2 points on globe
!
   use m_parameters, only : pifac_k1,earthr
   implicit none
!
   real(rkind1), intent(in)  ::  lat1,lon1,lat2,lon2 ! lat/lon pair
   real(rkind1), intent(out) ::  dist                ! great circle distance
!
   real(rkind1) ::  d1,d2,d3                 ! temp holders
   real(rkind1) ::  rlat1,rlon1,rlat2,rlon2  ! radian lat/lon
!
! convert to radians      
   rlat1 = pifac_k1*lat1
   rlon1 = pifac_k1*lon1
   rlat2 = pifac_k1*lat2
   rlon2 = pifac_k1*lon2
!
! estimate distance in km
   d1=sin(0.5*(rlat1-rlat2))
   d2=sin(0.5*(rlon1-rlon2))
   d3=cos(rlat1)*cos(rlat2)
   dist=2.*.001*earthr*asin(sqrt(d1*d1+d2*d2*d3))  
!     
   end subroutine amv_great_circle
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_check_loc (i,j,jbin,lat,lon,time_now, &
                             visir_here,wv_here)
!
! Consider this location for an AMV wind observation
! 1) Check if this point is in the subset considered for thinning
! 2) Check if this point can be viewed by the sat at this time
! 3) Check if solar elevation too low (or night) for this kx
! 4) Check if atmos profile will be required for either visir or wv obs
!
   use m_satloc, only : satloc_nbins,satloc_dlon,satloc_bins
   use m_parameters, only : earth_c
!   
   implicit none
!
   logical, intent(out) :: visir_here,wv_here
   integer, intent(in)  :: i,j,jbin
   real(rkind1), intent(in) :: lat,lon
   real(rkind1), intent(in) :: time_now(6)
!            
   logical :: lrange           
   integer :: k, ibin
   integer :: nlat
   integer :: nsum
   real(rkind1) :: xlon
   real(rkind1) :: geost_range
   real(rkind1) :: dist2nadir
   real(rkind1) :: solar_elev,solar_azim
!
! First check if this is a location that should be considered for
! this sat or data type, given the desired data thinning or density
   do k=1,kx_num
     if (mod(i-kx_i_stride(2,k),kx_i_stride(1,k)) == 0 .and. &
         mod(j-kx_j_stride(2,k),kx_j_stride(1,k)) == 0) then
       kx_present(k)=.true.  
     else
       kx_present(k)=.false.  
     endif
   enddo
!   
! Determine which hemsiphere point is in
   if (lat >= 0.) then
     nlat=2
   else
     nlat=1
   endif
!
   geost_range=60.*earth_c/360.  ! max distance from nadir (km) at equator
!
! Check whether this is a viewable location
   do k=1,kx_num
     if (kx_present(k)) then    ! This is a thinned point to consider
       kx_present(k)=.false.    ! Replace with default value of no view
!
       xlon=lon-satloc_dlon(1,k)
       if (xlon < 0.) then
         xlon=xlon+360. 
       endif
!
! Check if lat and lon within same rnage as real obs of this type
       call amv_latlon_range (k,lat,lon,lrange)
       if (lrange) then 
!
         ibin=1+int(xlon/satloc_dlon(2,k))
         if (ibin <= satloc_nbins) then
           if (satloc_bins(ibin,nlat,k) > 0) then   ! longitude OK
             if (kx_type(k)(1:1) == 'G') then       ! geostationary
               call amv_great_circle (lat,lon,0.,kx_locs(k),dist2nadir) 
               if (dist2nadir < geost_range) then   ! check on dist from nadir
                 kx_present(k)=.true.    
               endif
             else                                   ! polar
               kx_present(k)=.true.    
             endif  ! check on obs type
           endif    ! check on longitude
         endif      ! check on ibin allowed
!
       endif        ! check on lat/lon range
     endif          ! check on thinned grid
   enddo            ! loop over k
!
!
! Exclude points where sun elevation is too low for some obs types 
   call get_solar_elevation (time_now,lat,lon,solar_elev,solar_azim)
   do k=1,kx_num
     if (solar_elev < kx_filters(kx_filt_slev,k)) then  
       kx_present(k)=.false. 
     endif
   enddo
!
! Increment counter of OK obs locations
! Aslo determine if any of type VIS/IR or WV gradient are present here
   visir_here=.false.
   wv_here=.false.
   do k=1,kx_num
     if (kx_present(k)) then
       kx_pt_count(jbin,k,1)=kx_pt_count(jbin,k,1)+1
!
       if (kx_type(k)(2:2) == 'V' .or. kx_type(k)(2:2) == 'I' .or. &
           kx_type(k)(2:2) == 'B' ) then  ! VIS, IR, or both type
         visir_here=.true.
       else                              ! WV type
         wv_here=.true.
       endif
! 
     endif
   enddo     
!
   end subroutine amv_check_loc 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_stride_j (dlat,x_random)
!
! Determine latitudinal strides and starting indexes that define the 
! sets of points to consider as possible obs locations. These set vary 
! with the values of dx assigned to each kx*sat (index k). 
! If kx_igrid(k) > 1, the starting index is determined as a random  
! value within the length of the stride.  Otherwise it is set to 2 so that
! the starting value is the same for all times.
!
   use m_parameters, only : earth_c, pi_k1, pifac_k1
   implicit none
!
   real(rkind1), intent(in) :: dlat        ! spacing of lats in degrees
   real(rkind2), intent(in) :: x_random(kx_num) ! set of random numbers
!
   integer :: k,jx
   real(rkind1) :: xlat
!
   xlat=dlat*earth_c/360.   ! spacing of lats in km
   do k=1,kx_num            ! loop over all kx*sats
     kx_j_stride(1,k)=max(nint(kx_dx(k)/xlat),1)  ! approx dx in terms of xlat
     if (kx_ijgrid(k) == 1) then  ! same value for all times
       kx_j_stride(2,k)=2         
     else                         ! different (random) value for each time
       jx=int(real(kx_j_stride(1,k))*x_random(k))  
       kx_j_stride(2,k)=2+min(jx,kx_j_stride(1,k)-1) 
     endif
   enddo
!   
   end subroutine amv_stride_j 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_stride_i (lat,dlon,x_random)  
!
! Determine longitudinal strides and starting indexes that define the 
! sets of points to consider as possible obs locations. These sets 
! vary with the values of dx assigned to each kx*sat (index k).  They 
! can vary with lat also because the grid spacing varies in distance as 
! the meridians converge towards the pole. 
!
! An input value lat> 90. is used to initialize these values, particularly
! as required when starting i-indexes are to be independent of lat j. 
! If kx_igrid(k) = 3, the starting index for i is determined as a random  
! value for each time and lat j within the length of the stride, for which the 
! length in terms of numbers of grid points depends on the cosine of the lat. 
! If kx_igrid(k) = 1, the starting index is set to 1 so that it is the 
! same for all times and lats. 
! If kx_igrid(k) = 2, the starting index will vary with time but not lat.
!
   use m_parameters, only : earth_c, pifac_k1
   implicit none
!
   real(rkind1), intent(in) :: lat
   real(rkind1), intent(in) :: dlon
   real(rkind2), intent(in) :: x_random(kx_num) ! set of random numbers
!
   integer :: k,ix
   real(rkind1) :: slat           ! lat to use for determining stride
   real(rkind1) :: xlon           ! distance corresponding to dlon
!
   do k=1,kx_num
!
! If this either preceeds consideration of the first lat (so that the 
! starting index i is the same for each lat) or the starting index should 
! vary with ech lat, then ...
     if (lat > 90. .or. kx_ijgrid(k) == 3) then
       if (lat > 90.) then   ! set initial (or fixed) value of lat
         slat=abs(kx_filters(kx_filt_flat,k))
       else                  ! use current value of lat
         slat=lat
       endif
       xlon=dlon*earth_c*cos(slat*pifac_k1)/360.
       kx_i_stride(1,k)=max(nint(kx_dx(k)/xlon),1)
       if (kx_ijgrid(k) == 1) then    ! same for all time
         kx_i_stride(2,k)=1
       else                           ! varies randomly with time
         ix=int(real(kx_i_stride(1,k))*x_random(k))
         kx_i_stride(2,k)=1+min(ix,kx_i_stride(1,k)-1)
       endif
     endif
!
   enddo
!
   end subroutine amv_stride_i 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_normalize_counts (hnorm6) 
!
! Normalize counts so they are averages over 6 hour assimilation periods
! 
   use m_kinds, only : rkind1
   implicit none
!
   real(rkind1), intent(in) :: hnorm6
!
   kx_histogram(:,:,:,:,2)=kx_histogram(:,:,:,:,2)/hnorm6
   kx_pt_count(:,:,2)=kx_pt_count(:,:,2)/hnorm6
!
   end subroutine amv_normalize_counts
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_use_params (ifunc,k,jbin,ncldmax,ncld,cld_frac, &
                              cld_obsc,cld_puvs,id_ncld)
!
! Use a parameterized probabilistic model to determine if an observation
! is present:
! For ifunc=1, the probability is defined as p1*(1-abs(cf-0.5))**p2 where
! p1 and p2 are kx, p-layer, and lat/land-bin parameters. If p2<1, then
! the probability is one for p3<cf<p4.  Note also the condition on minimal 
! wind speed and obscuring cloud cover.
!
   implicit none
!
   integer, intent(in) :: k,jbin
   integer, intent(in) :: ifunc
   integer, intent(in) :: ncldmax
   integer, intent(in) :: ncld
   integer, intent(out) :: id_ncld(10)
   real(rkind1), intent(in) :: cld_frac(ncldmax)
   real(rkind1), intent(in) :: cld_obsc(ncldmax)
   real(rkind1), intent(in) :: cld_puvs(ncldmax,5)
!
   integer :: n
   integer :: np
   integer :: npm1
   real(rkind1), parameter :: epsilon=1.e-3 
   real(rkind1) :: dp
   real(rkind1) :: func
   real(rkind1) :: prob
   real(rkind1) :: x4
   real(rkind2) :: x
!
   npm1=kx_pbins-1
   dp=100000./real(kx_pbins)   ! dp in units of Pa
   id_ncld(:)=0
!
   if (ifunc == 1 .or. ifunc == 2) then ! only 1,2 implimented
!
     do n=1,ncld
       np=1+min(int(cld_puvs(n,1)/dp),npm1)     ! p-layer index
       if (kx_params(3,np,jbin,k)-epsilon < cld_frac(n) .and. & 
           kx_params(4,np,jbin,k)+epsilon > cld_frac(n) .and. &
           cld_obsc(n) < kx_cloud_obscure       .and. & 
           cld_puvs(n,4) > kx_speed_min) then    ! Cloud is observed 
!
! Set probability based on variable cf and parameters 
! cf here is either cloud fraction or scaled ipw difference in region   
         if (ifunc == 1) then 
           func=1.-abs(cld_frac(n)-0.5)
         elseif (ifunc == 2) then 
           func=cld_frac(n)
         endif
!
         if (kx_params(2,np,jbin,k) > .999) then  
           prob=kx_params(1,np,jbin,k)*(func**kx_params(2,np,jbin,k))
         else 
           prob=kx_params(1,np,jbin,k)
         endif
!
         call random_number (x)
         x4=real(x,rkind1)
         if (x4 < prob) then
           id_ncld(n)=1  ! flag
           kx_obs_count(np,jbin,k,1)=kx_obs_count(np,jbin,k,1)+1
         endif
       endif
     enddo
!
   endif  
!
   end subroutine amv_use_params
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_tune_params (ifunc,k)
!
!  Compute parameters for determining probabilities of an observation 
!  being present based on the distributions and values of a field-derived
!  value (such as cloud fraction).  
! 
   implicit none
!
   integer, intent(in) :: k
   integer, intent(in) :: ifunc
!
   logical :: lfound ! true if acceptable parameters found
   integer :: i,j,n,mp
   integer :: nc,nhalf,nc_min,nc_max
   integer :: histsum
   integer :: mpmax
   real(rkind1), parameter :: alpha_max=999.999 ! limted due to print format
   real(rkind1) :: cf
   real(rkind1) :: half_bin
   real(rkind1) :: frac
   real(rkind1) :: dpsum
   real(rkind1) :: prob
   real(rkind1) :: alpha
   real(rkind1) :: xcount
   real(rkind1) :: prob_max
   real(rkind1) :: func(kx_cbins)
   real(rkind1) :: cbins(kx_cbins,2)
!
! Inflate desired counts to adjust for GSI QC rejections
   do j=1,kx_jbins     ! loop over lat bands, sfc types, time slots
     do i=1,kx_pbins   ! loop over pressure bands
       xcount=real(kx_obs_count(i,j,k,2))
       kx_obs_count(i,j,k,2)=int(xcount*kx_qcfac(k))
     enddo
   enddo
!
   if (ifunc > 0) then  
!
! Define a value of cf for each c-bin
     half_bin=0.5/real(kx_cbins-1)
     do n=1,kx_cbins 
       cf=real(n-1)/real(kx_cbins-1)
       cbins(n,1)=cf          ! value of cf at center of distribution bin 
       cbins(n,2)=cf+half_bin ! value of cf at high-end edge of bin
     enddo
   endif
!
! Define functional forms being tested (exclusive of the exponent parameter m) 
   if (ifunc == 1) then  ! consider functions P=alpha*(1.-abs(cf-0.5))**m
     mpmax=5
     func(1)=0.          ! bin with cf=0 is unobservable
     do n=2,kx_cbins
       func(n)=1.-abs(cbins(n,1)-0.5)
     enddo
!
   elseif (ifunc == 2) then  ! consider functions P=alpha*cf**m
     mpmax=3                 
     func(1)=0.              ! bin with cf=0 is unobservable
     do n=2,kx_cbins
       func(n)=cbins(n,1)
     enddo
   endif  
!
   do j=1,kx_jbins     ! loop over lat bands, sfc types, time slots
     do i=1,kx_pbins   ! loop over pressure bands
!
       histsum=sum(kx_histogram(:,i,j,k,2))
!
       lfound=.false.       ! initial value
       if (ifunc <= 2) then
         do mp=mpmax,1,-1   ! consider possible values of param 2
!
! Determine sum of products of conditional probability function being 
! tested (P of obs existing given cloud fraction at this location, 
! excluding the proportionality factor to be determined) with number of  
! points having this cloud fraction.
           dpsum=0.
           prob_max=0.
           do n=1,kx_cbins
             if (func(n) > 0.) then 
               prob=func(n)**mp  ! function being tested
             else
               prob=0.
             endif
             prob_max=max(prob,prob_max) 
             dpsum=dpsum+prob*kx_histogram(n,i,j,k,2) 
           enddo
!
! Determine the required factor to obtain the desired counts
! Check whether this factor is suitable    
           if (dpsum > 0.) then  
             alpha=kx_obs_count(i,j,k,2)/dpsum
             if (alpha <= alpha_max .and. &
                 alpha*prob_max < 1.0001) then ! required for a probability
               lfound=.true.
               kx_params(1,i,j,k)=alpha
               kx_params(2,i,j,k)=real(mp)
               kx_params(3,i,j,k)=cbins(1,2)         ! allow all bins except 
               kx_params(4,i,j,k)=cbins(kx_cbins,1)  ! the first that contains c=0
             endif
           endif
           if (lfound) exit  ! no need for futher consideration
         enddo               ! end loop over mp 
       endif                 ! test of function form id 
!
! If suitable form of considered function not yet found, then 
! consider a conditional P that is a constant within a range 
! centered on the middle of possible values.
       if (.not. lfound) then 
         dpsum=0.
         nhalf=1+kx_cbins/2  ! index of center of range of values
         nc_max=nhalf        ! lowest index
!
         do n=1,kx_cbins-1 
!
! For ifunc=1,
! alternately consider bins >= nhalf and < nhalf, starting at nc=nhalf
! Do not consider bin 1 that generally contains value 0 (e.g., clear sky) 
           if (ifunc == 1) then
             if (mod(n,2) == 1) then
               nc=nhalf-n/2
               nc_min=nc
             else
               nc=nhalf+n/2
               nc_max=nc
             endif
!
! For ifunc=2,
! alternately consider bins in reverse order (i.e., for largest cf first)
! Do not consider bin 1 that generally contains value 0 (e.g., clear sky) 
           else
             nc=kx_cbins+1-n 
             nc_min=nc
             nc_max=kx_cbins 
           endif 
!
! Sum distribution with range of bins considered            
           if (nc > 1 .and. nc <= kx_cbins) then  
             dpsum=dpsum+kx_histogram(nc,i,j,k,2)
!
! If considered sunbet of bins include a total fraction of 
! occurance greater than desired frac of points with observations, 
! then restrict range to this subset and compute alpha based
! on ratio of these fractions.  
             if (dpsum > kx_obs_count(i,j,k,2)) then
               lfound=.true.  
               kx_params(1,i,j,k)=kx_obs_count(i,j,k,2)/dpsum
               kx_params(2,i,j,k)=0.
               kx_params(3,i,j,k)=cbins(nc_min-1,2)
               kx_params(4,i,j,k)=cbins(nc_max,2)
             endif
           endif   ! consideration of bin index nc
!
           if (lfound) exit
         enddo     ! loop over numbers of bins
!
! If no suitable function or range found, then set P=1 for all bin values
! exept for values within the first bin.           
         if (.not. lfound) then
           lfound=.true.  
           if (kx_obs_count(i,j,k,2) == 0) then 
             kx_params(1,i,j,k)=0.
           else
             kx_params(1,i,j,k)=1. 
           endif
           kx_params(2,i,j,k)=0.
           kx_params(3,i,j,k)=cbins(1,2)
           kx_params(4,i,j,k)=cbins(kx_cbins,1)
         endif  ! 2nd test on lfound
       endif    ! 1st test on lfound  
!  
! Ensure that alpha is not too large to write values to text file
!
       kx_params(1,i,j,k)=min(kx_params(1,i,j,k), alpha_max)
!
     enddo      ! loop over pressure bands
   enddo        ! loop over other (lat, land/se, time?) bands
!
   end subroutine amv_tune_params
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_accum_bins (ncldmax,ncld,k,jbin,cld_obsc,cld_puvs, &
                              cld_frac)
!
! Accumulate counts in cloud distribution histogram.
! Exclude obs obscured from above or cloud speeds les than a min value.
! If there are no clouds within a particular pressure bin at this location
! that have at least the min speed, then count it as a cloud free point. 
! 
   implicit none
!
   integer, intent(in) :: ncldmax,ncld
   integer, intent(in) :: k,jbin
   real(rkind1), intent(in) :: cld_puvs(ncldmax,5) ! p,u,v,speed of cloud
   real(rkind1), intent(in) :: cld_frac(ncldmax)   ! max cloud frac in layer
   real(rkind1), intent(in) :: cld_obsc(ncldmax)   ! opacity above cloud layer
!
   integer :: npm1,ncm1,n,np,nc
   integer :: icount(kx_cbins,kx_pbins)
   integer :: isum
   real(rkind1) :: dp,db,db2
!
   npm1=kx_pbins-1
   dp=100000./real(kx_pbins)   ! dp in units of Pa
   ncm1=kx_cbins-1
   db=1./real(ncm1)
   db2=db/2.
!
   icount(:,:)=0    ! default value: no amv for any bin for any p-level bin
   do n=1,ncld      ! consider all cloud layers at this loc
     if (cld_obsc(n) < kx_cloud_obscure .and. & 
         cld_puvs(n,4) > kx_speed_min) then     ! Cloud is observed 
       np=1+min(int(cld_puvs(n,1)/dp),npm1)     ! p is in this p-bin
       nc=1+min(int((cld_frac(n)+db2)/db),ncm1) ! cf is in this c-bin
       icount(nc,np)=icount(nc,np)+1
     endif
   enddo
!
! Note that for kx_histogram, the total over all bins can vary with the 
! pbin index (np) because at any considered grid point, the numbers of 
! cloud levels viewed can vary between the pressure ranges.

   do np=1,kx_pbins          ! loop over p-level bins
     isum=sum(icount(:,np))  ! see if any bins populated for this p-level bin
     if (isum == 0) then     ! no amv withim this p-bin (c-bin 1 for clear sky
       kx_histogram(1,np,jbin,k,1)=kx_histogram(1,np,jbin,k,1)+1
     else
       do nc=1,kx_cbins      ! increment counts for all c-bins for this p-bin
         kx_histogram(nc,np,jbin,k,1)=kx_histogram(nc,np,jbin,k,1)+icount(nc,np)
       enddo
     endif
   enddo
!
   end subroutine amv_accum_bins
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_speed_cld (ncld,ncldmax,klevs,k,cld_klev,cld_plev, &
                             u,v,z,cld_puvs)
!
! Compute pressure, wind vector, and speed at assigned level of cloud.
! The assigned level may be the cloud top or bottom, depending on kx.
! This dependence exists because of the way the cloud height is 
! estimated from real images.
!
   use m_kinds, only : rkind1
   implicit none
!
   integer, intent(in) :: ncldmax, klevs
   integer, intent(in) :: ncld
   integer, intent(in) :: k
   integer, intent(in) :: cld_klev(ncldmax,2) ! top and bot of each layer
   real(rkind1), intent(in) :: cld_plev(ncldmax,2) ! top and bot of each layer
   real(rkind1), intent(in) :: u(klevs),v(klevs),z(klevs)   ! u,v,z profiles
   real(rkind1), intent(out) :: cld_puvs(ncldmax,5)! p,u,v,speed,z of cloud
!
   integer :: n, id
   real(rkind1) :: speed2
!
! Below, kx_filters(kx_filt_torb,k) is assumed in hPa, but cld_plev in Pa
   do n=1,ncld
     if (cld_plev(n,1) > 100.*kx_filters(kx_filt_torb,k)) then
       id=2  ! bottom of cloud seen by images
     else        
       id=1  ! top of cloud seen by images
     endif
     cld_puvs(n,1)=cld_plev(n,id)
     cld_puvs(n,5)=z(cld_klev(n,id))
     cld_puvs(n,2)=u(cld_klev(n,id))
     cld_puvs(n,3)=v(cld_klev(n,id))
     speed2=cld_puvs(n,2)*cld_puvs(n,2)+cld_puvs(n,3)*cld_puvs(n,3)
     if (speed2 > 1.e-20) then 
       cld_puvs(n,4)=sqrt(speed2)
     else
       cld_puvs(n,4)=0.
     endif
   enddo
!
   end subroutine amv_speed_cld
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_speed_ipw (klevs,k,p,q,u,v,z,ipw_puvs)
!
! Compute pressure, wind vector, and speed at assigned level of cloud.
! The assigned level may be the cloud top or bottom, depending on kx.
! This dependence exists because of the way the cloud height is 
! estimated from real images.
!
   use m_kinds, only : rkind1
   use m_kx_table, only : kx_ipw_nlevs,kx_ipw_plevs,kx_ipw_bparm,kx_ipw_ids
!
  implicit none
!
   integer, intent(in) :: klevs  ! number of levels in profiles
   integer, intent(in) :: k      ! id number for kx, sat id pair in table
   real(rkind1), intent(in) :: p(klevs,3) ! p at data and 2 interface levels
   real(rkind1), intent(in) :: u(klevs),v(klevs),q(klevs),z(klevs) ! profile
   real(rkind1), intent(out) :: ipw_puvs(kx_ipw_nlevs,5) ! p,u,v,speed of cloud
!
   integer :: n,id,nplev
   integer :: k_bparm
   integer :: kid(kx_ipw_nlevs) 
   real(rkind1) :: pa,pb,speed2
!
! Loop over the number of ipw layers defined in the kx_table file
! If either this obs type (with kx_table id k) uses specifically this layer 
! (as indicated by kx_ipw_bparm(1,n,kx_ipw_ids(k))>0.) then find the 
! profile level where the first value of q is found for which 
! q > kx_ipw_bparm(2,n,kx_ipw_ids(k)), searching starting from  
! the top of the layer. This defines the level assigned to the obs at this 
! location. Note that if more than one ipw layer is considered, the first 
! layer for which a desired q is found determines which layer is considered
! for the observation.
!
! Indicate which set of ipw_bparm values should be used for this obs type k
   k_bparm=kx_ipw_ids(k)
!
   do n=1,kx_ipw_nlevs
     nplev=0
     if (kx_ipw_bparm(1,n,k_bparm) > 0.) then
       pa=100.*kx_ipw_plevs(1,n)  ! change to units of Pa
       pb=100.*kx_ipw_plevs(2,n)  
       do id=1,klevs
         if (nplev == 0 .and. p(id,2) >= pa .and. p(id,2) <= pb) then   
           if ( q(id) > kx_ipw_bparm(2,n,k_bparm)) then
             nplev=id
           endif
         endif 
       enddo
     endif 
     kid(n)=nplev
   enddo
!
   do n=1,kx_ipw_nlevs
     id=kid(n)      
     if (id > 0) then  
       ipw_puvs(n,1)=p(id,2)
       ipw_puvs(n,5)=z(id)
       ipw_puvs(n,2)=u(id)
       ipw_puvs(n,3)=v(id)
       speed2=ipw_puvs(n,2)*ipw_puvs(n,2)+ipw_puvs(n,3)*ipw_puvs(n,3)
       if (speed2 > 1.e-20) then 
         ipw_puvs(n,4)=sqrt(speed2)
       else
         ipw_puvs(n,4)=0.
       endif
     else
       ipw_puvs(n,1)=p(klevs,2)
       ipw_puvs(n,5)=0.
       ipw_puvs(n,2)=0.
       ipw_puvs(n,3)=0.
       ipw_puvs(n,4)=0.
     endif
   enddo
!
   end subroutine amv_speed_ipw
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_obscured (ncld,ncldmax,cld_frac,cld_obsc)
!
! Compute probability that cloud layer is obscured by clouds above by 
! assuming random overlap between the different cloud layers. The combination 
! of assuming maximumum overlap within a cloud layer but random overlap 
! between layers is called the maximum random overlap assumption. 
!
   implicit none
!
   integer, intent(in) :: ncld, ncldmax
   real(rkind1), intent(in) :: cld_frac(ncldmax)   ! max cloud frac in layer
   real(rkind1), intent(out) :: cld_obsc(ncldmax)  ! opacity above cloud layer
!
   integer :: k
   real(rkind1) :: cld_clear(ncldmax)  ! work space   
!
   cld_clear(:)=one_k1
   do k=2,ncld 
     cld_clear(k)=cld_clear(k-1)*(one_k1-cld_frac(k-1))
   enddo
!
   do k=1,ncld
     cld_obsc(k)=min(max(one_k1-cld_clear(k),zero_k1),one_k1)
   enddo
!
   end subroutine amv_obscured 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   subroutine amv_layers (klevs,ncldmax,plevs,cf,ncld, &
                          cld_klev,cld_frac,cld_plev)
!
   use m_kinds, only : rkind1
   implicit none
!
! Find cloud layers (clouds at consecutive levels)
! Within each layer, maximum overlap of cloud fraction is assumed, so that
! the cloud fraction assigned to the layer is the maximum cloud fraction 
! of all the levels within that layer. 
!
   integer, intent(in) :: klevs, ncldmax
   real(rkind1), intent(in) :: cf(klevs) 
   real(rkind1), intent(in) :: plevs(klevs,3) 
!
   integer, intent(out) :: ncld
   integer, intent(out) :: cld_klev(ncldmax,2)      ! cloud layer top and bot id
   real(rkind1), intent(out) :: cld_frac(ncldmax)   ! max cloud frac in layer
   real(rkind1), intent(out) :: cld_plev(ncldmax,2) ! cloud layer top and bot p
!
   logical :: cf_above
   integer :: k
   real(rkind1), parameter :: cf_min=0.02  ! min value to consider for cf
! 
   ncld=0            ! initialize counter of cloud layers
   cf_above=.false.  ! indicates no clouds in level above
!
! The top of each cloud layer is defined as that of the highest model 
! data level in which the cloud resides. The bottom is first defined 
! defined by same level, but then is overwritten by that for successive 
! model layers below that are still within the contiguous cloud. 
! plevs(k,2) refers to p at data level k, determined as the average of the
! p at the adcacent interface levels above and below. 
!
   do k=1,klevs                
     if (cf(k) > cf_min ) then       ! cloud present at this k
       if (cf_above) then            ! within present cloud layer
         cld_plev(ncld,2)=plevs(k,2) ! will be last level in layer  
         cld_frac(ncld)=max(cf(k),cld_frac(ncld)) ! max cf in layer  
         cld_klev(ncld,2)=k          ! will be last level in layer   
       elseif (ncld < ncldmax) then  ! top of new cloud layer
         cf_above=.true.             ! now within cloud layer
         ncld=ncld+1                 ! new cloud layer
         cld_plev(ncld,1)=plevs(k,2) ! cloud top p 
         cld_frac(ncld)=cf(k)        ! cloud frac at layer top
         cld_klev(ncld,1:2)=k        ! index of cloud layer top
         cld_plev(ncld,2)=plevs(k,2) ! p of first data level in layer
       endif
     else
       cf_above=.false.               ! now within cloud free layer
     endif             
   enddo
!
   end subroutine amv_layers
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
   end module m_amv_view
