  program overlap 
!
  implicit none
!
  logical :: ltest
  integer :: i,j
  inTeger :: nlon1,nlat1,nlon2,nlat2
  integer :: i1,j1           ! grid 1 indices
  integer :: i2a,i2b,j2a,j2b ! grid 2 indeces
  real :: deg2rad            ! (pi/4) radians / 45 degrees 
  real :: rnlon1              ! 1/nlon1
  real :: sumw               ! sumed weights for testing (=2)
  real :: wja,wjb,wia,wib    ! 1-d weighting factors
  real :: xa,xb,ya,yb        ! lat and lons 
  real :: dlon1,dlat1,dlon2,dlat2
  real, allocatable :: lon1(:),lat1(:),sinlat1(:),area1(:)
  real, allocatable :: lon2(:),lat2(:),sinlat2(:),area2(:)
  real, allocatable :: weights_id(:,:,:),weights(:,:,:)
  real, allocatable :: test(:,:)
!
  ltest=.true.
  nlon1=36
  nlat1=7
  nlat2=7
  nlon2=12
  nlon1=5760
  nlat1=2881
  nlat2=1152
  nlon2=721
!
  allocate (lon1(0:nlon1))
  allocate (lat1(0:nlat1))
  allocate (sinlat1(0:nlat1))
  allocate (area1(nlat1))
!
  allocate (lon2(0:nlon2))
  allocate (lat2(0:nlat2))
  allocate (sinlat2(0:nlat2))
  allocate (area2(nlat2))
!
  allocate (weights_id(4,nlon1,nlat1))
  allocate (weights(4,nlon1,nlat1))
  weights(:,:,:)=0.
!
  dlon1=360./nlon1
  dlat1=180./(nlat1-1)
  dlon2=360./nlon2
  dlat2=180./(nlat2-1)
!
  lat1(0)=-90.-0.5*dlat1
  do j=1,nlat1-1
   lat1(j)=lat1(j-1)+dlat1
  enddo
  lat1(0)=-90.
  lat1(nlat1)=90.
!
  lat2(0)=-90.-0.5*dlat2
  do j=1,nlat2-1
   lat2(j)=lat2(j-1)+dlat2
  enddo
  lat2(0)=-90.
  lat2(nlat2)=90.
!
  deg2rad=atan(1.)/45.
  do j=1,nlat1-1
    sinlat1(j)=sin(deg2rad*lat1(j))
  enddo
  sinlat1(0)=-1.  
  sinlat1(nlat1)=1.  
!
  do j=1,nlat2-1
    sinlat2(j)=sin(deg2rad*lat2(j))
  enddo
  sinlat2(0)=-1.  
  sinlat2(nlat2)=1.  
!
  do j=1,nlat1
   area1(j)=(sinlat1(j)-sinlat1(j-1))/nlon1
  enddo
!
  do j=1,nlat2
   area2(j)=(sinlat2(j)-sinlat2(j-1))/nlon2
  enddo
!
  do j=1,nlat1
!    print ('(a,i3,4e15.5)'),'GRID1',j,lat1(j),sinlat1(j),area1(j)
  enddo
  sumw=nlon1*sum(area1(:))
  print *,'SUM1=',sumw
  do j=1,nlat2
!    print ('(a,i3,4e15.5)'),'GRID2',j,lat2(j),sinlat2(j),area2(j)
  enddo
  sumw=nlon2*sum(area2(:))
  print *,'SUM2=',sumw
!
  lon1(0)=-0.5*dlon1
  do i=1,nlon1
    lon1(i)=lon1(i-1)+dlon1
  enddo
  rnlon1=1./real(nlon1)
!
  lon2(0)=-0.5*dlon2
  do i=1,nlon2
    lon2(i)=lon2(i-1)+dlon2
  enddo
!
  sumw=0.
!
  do j1=1,nlat1     ! loop over grid 1 lats
    ya=lat1(j1-1)   ! S. edge of box on grid 1 
    yb=lat1(j1)     ! N. edge of box on grid 1 
    j2a=1+int(0.5+(ya+90.)/dlat2) ! lat index for grid 2 box containing ya
    j2b=1+int(0.5+(yb+90.)/dlat2) ! lat index for grid 2 box containing yb
    if (j2a==j2b) then               ! lats of grid 1 box fully in grid 2 box 
      wja=sinlat1(j1)-sinlat1(j1-1)  ! no need to consider 2 boxes in N-S 
      wjb=0.                         ! so second box contribution can be 0
    else                             ! compute N-S distance of part of 
      wja=sinlat2(j2a)-sinlat1(j1-1) ! grid 1 box j1 in grid 2 box j2a
      wjb=sinlat1(j1)-sinlat2(j2a)   ! grid 1 box j1 in grid 2 box j2b
    endif
!
    do i1=1,nlon1     ! loop over grid 1 lons 
      xa=lon1(i1-1)   ! W. edge of box on grid 1      
      xb=lon1(i1)     ! E. edge of box on grid 1      
      i2a=1+int(0.5+xa/dlon2)  ! lon index for grid 2 box containing xa
      i2b=1+int(0.5+xb/dlon2)  ! lon index for grid 2 box containing xb
      if (i2a>nlon2) i2a=1
      if (i2b>nlon2) i2b=1
      if (i2a==i2b) then                ! grid 1 box fully in grid 2 box 
        wia=rnlon1                      ! no need to consider 2 boxes in E-W 
        wib=0.                          ! so second box contribution can be 0
      else                              ! compute E-W distance of part of 
        wia=(lon2(i2a)-lon1(i1-1))/360. ! grid 1 box i1 within grid 2 box i2a
        wib=(lon1(i1)-lon2(i2a))/360.   ! grid 1 box i1 within grid 2 box i2b
      endif
!
! weights_id are the 4 i,j indeces of the 4 possible adjacent grid 2 boxes 
! overlayed by the grid 1 box i1,j1 
!
      weights_id(1,i1,j1)=real(i2a)
      weights_id(2,i1,j1)=real(i2b)
      weights_id(3,i1,j1)=real(j2a)
      weights_id(4,i1,j1)=real(j2b)
!
! weights are the fractions of area of grid 1 box i1,j1 that contribute to the 
! areas of each of the 4 possible grid 2 boxes 1 box (i2a,j2a), (i2b,j2a),  
! (i2a,j2b), (i2b,j2b)
!
      weights(1,i1,j1)=wja*wia/area2(j2a)  
      weights(2,i1,j1)=wja*wib/area2(j2a)
      weights(3,i1,j1)=wjb*wia/area2(j2b)
      weights(4,i1,j1)=wjb*wib/area2(j2b)
      sumw=sumw+(wja+wjb)*(wia+wib)
    enddo
  enddo 
!
  print ('(a,e15.5)'),'SUM=',sumw
!
  if (ltest) then
    allocate (test(nlon2,nlat2))
    test(:,:)=0.
    do j1=1,nlat1
      do i1=1,nlon1
        i2a=nint(weights_id(1,i1,j1))
        i2b=nint(weights_id(2,i1,j1))
        j2a=nint(weights_id(3,i1,j1))
        j2b=nint(weights_id(4,i1,j1))
        test(i2a,j2a)=test(i2a,j2a)+weights(1,i1,j1)
        test(i2b,j2a)=test(i2b,j2a)+weights(2,i1,j1)
        test(i2a,j2b)=test(i2a,j2b)+weights(3,i1,j1)
        test(i2b,j2b)=test(i2b,j2b)+weights(4,i1,j1)
      enddo
    enddo
!
    print *,'Test sample of weights:  All should equal 1'
    print ('(a)'),'Lat 1'  
    print ('(10f10.6)'),test(:,1)
    print ('(a,i5)'),'Lat',nlat2
    print ('(10f10.6)'),test(:,nlat2)  
    print ('(a)'),'Lon 1'  
    print ('(10f10.6)'),test(1,:)
    print ('(a,i5)'),'Lon',nlon2
    print ('(10f10.6)'),test(nlon2,:)  
    print ('(a,i5)'),'Lat',nlat2/3  
    print ('(10f10.6)'),test(:,nlat2/3)
    print ('(a,i5)'),'Lat',3*nlat2/4
    print ('(10f10.6)'),test(:,3*nlat2/4)  
    print ('(a,i5)'),'Lon',nlon2/3  
    print ('(10f10.6)'),test(nlon2/3,:)
    print ('(a,i5)'),'Lon 2'
    print ('(10f10.6)'),test(2,:)  
!
    deallocate (test)
  endif ! test on ltest
!
  deallocate (lon1,lat1,area1,sinlat1)
  deallocate (lon2,lat2,area2,sinlat2)
  deallocate (weights_id,weights)

  end program overlap   

