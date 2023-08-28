   program test_gsi
!
   implicit none
!   
  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
  integer, parameter :: nsig=72
!
  real(r_kind) :: xlon, xlat
  real(r_kind) :: tges(nsig)  ! tv                (bottom to top)
  real(r_kind) :: qges(nsig)  ! specific humidity (bottom to top)
  real(r_kind) :: hges(nsig)  ! height (m) at interface below
  real(r_kind) :: prsltmp(nsig+1)  ! log ps(cb) at interfaces
  real(r_kind) :: zsges       ! sfc h (m)
  real(r_kind) :: rocprof     ! local radius of curvature (m)
  real(r_kind) :: unprof      ! geoid undulation (m)
  real(r_kind) :: tpdpres(1)  ! impact parameter (m)
  real(r_kind) :: bend_angle  
!
  integer :: i
!
  call get_data (nsig,xlon,xlat,tges,qges,hges,prsltmp,zsges, &
                 rocprof,unprof,tpdpres)
  tpdpres(1)=6.410037E+06
  call setupbend (xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                  tges,qges,hges,prsltmp,bend_angle)
  print *,'QQQ ZZZ',tpdpres(1),bend_angle

  tpdpres(1)=6.420062E+06
  call setupbend (xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                  tges,qges,hges,prsltmp,bend_angle)
  print *,'QQQ ZZZ',tpdpres(1),bend_angle
  tpdpres(1)=6.430184E+06
  call setupbend (xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                  tges,qges,hges,prsltmp,bend_angle)
  print *,'QQQ ZZZ',tpdpres(1),bend_angle
  tpdpres(1)=6.440018E+06
  call setupbend (xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                  tges,qges,hges,prsltmp,bend_angle)
  print *,'QQQ ZZZ',tpdpres(1),bend_angle


!
  end program test_gsi
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
  subroutine get_data (nsig,xlon,xlat,tges,qges,hges,prsltmp,zsges, &
                       rocprof,unprof,tpdpres)
  implicit none
!   
  integer, parameter :: r_kind=8
  integer, parameter :: i_kind=4 
  integer, parameter :: iunit=10
!
  integer :: nsig
  integer :: i,ix
!
  real(r_kind) :: xlon, xlat
  real(r_kind) :: tges(nsig)  ! tv                (bottom to top)
  real(r_kind) :: qges(nsig)  ! specific humidity (bottom to top)
  real(r_kind) :: hges(nsig)  ! height (m) at interface below
  real(r_kind) :: prsltmp(nsig+1)  ! log ps(cb) at interfaces
  real(r_kind) :: zsges       ! sfc h (m)
  real(r_kind) :: rocprof     ! local radius of curvature (m)
  real(r_kind) :: unprof      ! geoid undulation (m)
  real(r_kind) :: tpdpres(1)  ! impact parameter (m)
  real(r_kind) :: pges(nsig)  ! pressure at interface below (Pa)
!
  character(len=4) :: cdum1
!
  xlon=180.312
  xlat=51.
  rocprof=6.39242E+06
  unprof=-44.49_r_kind
!
  open (iunit,file='profdata.txt')
  do i=1,nsig
    read (iunit,*) cdum1,ix,tges(i),qges(i),hges(i),prsltmp(i)
    pges(i)=10._r_kind*exp(prsltmp(i)) ! p in Pa
  enddo
  close (iunit)
  prsltmp(nsig+1)=log(1./10.)
!
  print ('(a)'),'Field Profile read'
  i=2
  print ('(5e15.5)'),tges(i),qges(i),hges(i),prsltmp(i)
!
  end subroutine get_data 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
  subroutine setupbend (xlon,xlat,rocprof,unprof,tpdpres,zsges, &
                        tges,qges,hges,prsltmp,bend_angle)
!
  implicit none
!
  logical, parameter :: use_compress=.false. 
!
  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
  integer, parameter :: nsig=72
  integer, parameter :: nsig_ext=13
  integer, parameter :: grids_dim=80
  integer, parameter :: nobs=1
!
! parameters replacing gsi use variables
  real(r_kind), parameter :: zero=0._r_kind
  real(r_kind), parameter :: one=1._r_kind
  real(r_kind), parameter :: two=2._r_kind
  real(r_kind), parameter :: three=3._r_kind
  real(r_kind), parameter :: four=4._r_kind
  real(r_kind), parameter :: five=5._r_kind
  real(r_kind), parameter :: half=0.5_r_kind
  real(r_kind), parameter :: eccentricity=8.18192E-02
  real(r_kind), parameter :: semi_major_axis=6.37814E+06
  real(r_kind), parameter :: semi_minor_axis=6356.7523142e3_r_kind 
  real(r_kind), parameter :: grav_equator=9.78033E+00
  real(r_kind), parameter :: grav_polar= 9.8321849378_r_kind 
  real(r_kind), parameter :: grav_ratio=3.44979E-03
  real(r_kind), parameter :: grav=9.80665e+0_r_kind
  real(r_kind), parameter :: rd=2.8706e+2_r_kind  ! 2.8705e+2_r_kind
  real(r_kind), parameter :: rv=4.6150e+2_r_kind
  real(r_kind), parameter :: eps=0.622              ! rd/rv
  real(r_kind), parameter :: fv= (one/eps -one )      ! rv/rd-one
  real(r_kind), parameter :: pifac=3.14159265358979_r_kind
  real(r_kind), parameter :: deg2rad=pifac/180._r_kind
!
! other parameters
  real(r_kind),parameter::  r240 = 240.0_r_kind
  real(r_kind),parameter:: six = 6.0_r_kind
  real(r_kind),parameter:: ten = 10.0_r_kind
  real(r_kind),parameter:: eight = 8.0_r_kind
  real(r_kind),parameter:: nine = 9.0_r_kind
  real(r_kind),parameter:: eleven = 11.0_r_kind
  real(r_kind),parameter:: ds=10000.0_r_kind
  real(r_kind),parameter:: r12=12.0_r_kind
  real(r_kind),parameter:: r18=18.0_r_kind
  real(r_kind),parameter:: r20=20.0_r_kind
  real(r_kind),parameter:: r40=40.0_r_kind
  real(r_kind),parameter:: r1em3 = 1.0e-3_r_kind
  real(r_kind),parameter:: r1em6 = 1.0e-6_r_kind
  character(len=*),parameter :: myname='setupbend'
  real(r_kind),parameter:: crit_grad = 157.0_r_kind
!
  real(r_kind) :: n_a,n_b,n_c
  real(r_kind) :: somigliana
  real(r_kind) :: flattening


! Declare passed variables

  real(r_kind) :: xlon, xlat
  real(r_kind) :: tges(nsig)  ! tv                (bottom to top)
  real(r_kind) :: qges(nsig)  ! specific humidity (bottom to top)
  real(r_kind) :: hges(nsig)  ! height (m) at interface below
  real(r_kind) :: prsltmp(nsig+1)  ! log ps(cb) at interfaces
  real(r_kind) :: zsges       ! sfc h (m)
  real(r_kind) :: rocprof     ! local radius of curvature (m)
  real(r_kind) :: unprof      ! geoid undulation (m)
  real(r_kind) :: tpdpres(1)  ! impact parameter (m)
  real(r_kind) :: bend_angle  ! bend angle output

 !QQQQQQQQ
  logical :: lprintobs                                                               
  real(r_kind),parameter:: ron_lat=51.   ! -80.  51.                                 
  real(r_kind),parameter:: ron_lon=180.31     !   90. 180.31                         
  real(r_kind) :: ron_x1, ron_x2, ron_x3                                             
 


! Declare local parameters
! Declare local variables

  real(r_kind) sin2
  real(r_kind),dimension(grids_dim):: ddnj,grid_s,ref_rad_s

  real(r_kind) rsig,rsig_up,ddbend,tmean,qmean
  real(r_kind) termg,termr,termrg,hob,dbend,grad_mod
  real(r_kind) fact,pw,nrefges1,nrefges2,nrefges3,k4,delz
  real(r_kind) ratio,residual,obserror,obserrlm

  real(r_kind),dimension(nsig):: dbenddn,dbenddxi
  real(r_kind) pressure,hob_s,d_ref_rad,hob_s_top
  real(r_kind),dimension(4) :: w4,dw4
  real(r_kind),dimension(nsig,1) :: rges
  real(r_kind) :: nrefges(200,1)
  real(r_kind) :: xj(200,1)
  real(r_kind) :: dbend_loc(200,1)
  
  integer(i_kind) i,j,k,kk,mreal,nreal,jj,ikxx,ibin
  integer(i_kind) mm1,nsig_up,ihob,istatus
  integer(i_kind) kprof,istat,k1,k2,nobs_out,top_layer_SR,bot_layer_SR,count_SR
  integer(i_kind),dimension(4) :: gps_ij

  real(r_kind),dimension(3,nsig+nsig_ext) :: q_w
  real(r_kind),dimension(nsig) :: irefges,zges
  real(r_kind),dimension(0:nsig+nsig_ext+1) :: ref_rad
  real(r_kind),dimension(nsig+nsig_ext+20) :: ref_rad_out
  real(r_kind) :: dlat,dlon,dtime,dpressure,trefges,qrefges
  real(r_kind) :: dbetan,dbetaxi,rdog  

  somigliana = &
       (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
  flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis


!   Define refractive index coefficients here    
    if (use_compress) then
       ! Constants for gpsro data (Rueger 2002) 
       n_a = 77.6890_r_kind   ! K/mb             
       n_b = 3.75463e+5_r_kind  ! K^2/mb      
       n_c = 71.2952_r_kind   ! K/mb   
    else
       ! Constants for gpsro data (Bevis et al 1994)    
       n_a = 77.60_r_kind     ! K/mb  
       n_b = 3.739e+5_r_kind  ! K^2/mb 
       n_c = 70.4_r_kind      ! K/mb   
    endif


! Intialize variables
  nsig_up=nsig+nsig_ext ! extend nsig_ext levels above interface level nsig
  rsig=float(nsig)
  rdog=rd/grav
  rsig_up=float(nsig_up)
  nobs_out=0
  hob_s_top=one

! define new equally spaced grid s
  do j=0,grids_dim-1
     grid_s(j+1)=j*ds
  enddo

! A loop over all obs.
  loopoverobs1: &
  do i=1,1 ! loop over obs 

     sin2  = sin(xlat*deg2rad)**2

! Compute refractivity index-radius product at interface
!
!    Convert geopotential height at layer midpoints to geometric height using
!    equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
!    measures of altitude" (2001).  Available on the web at
!    http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!
!    termg  = equation 17
!    termr  = equation 21
!    termrg = first term in the denominator of equation 23
!    zges   = equation 23
 
     termg = grav_equator * &
             ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
     termr = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sin2)
     termrg = (termg/grav)*termr


 print ('(a,i4,1p6e15.5)'), 'QQQ PLANE',nsig,zsges,xlat,xlon
 do k=1,nsig                                                                  
   print ('(a,i4,1p4e15.5)'),'QQQA',k,tges(k),qges(k),hges(k),prsltmp(k)    
 enddo                                                                        
 print ('(a,1p5e15.5)'), 'QQQ GRAV A', termg,termr,termrg,grav_equator,eccentricity
 print ('(a,1p3e15.5)'), 'QQQ GRAV B', semi_major_axis,flattening,grav_ratio
   


     do k=nsig,1,-1 
        zges(k) = (termr*hges(k)) / (termrg-hges(k)) ! eq (23) at
        rges(k,i) = zges(k) + zsges + unprof + rocprof ! radius r_i

print ('(a,i4,1p6e15.5)'),'QQQ RAD', k,zges(k),rges(k,i),unprof,rocprof,zsges


        if(k>1) then
           qmean=(qges(k)+qges(k-1))/two
           tmean=(tges(k)+tges(k-1))/two
        else
           qmean=qges(1)
           tmean=tges(1)
        endif
        fact=(one+fv*qmean)
        pw=eps+qmean*(one-eps)
        k4=n_c-n_a
        pressure=ten*exp(prsltmp(k)) ! pressure of interface level in mb
        nrefges1=n_a*(pressure/tmean)*fact
        nrefges2=n_b*qmean*pressure*fact**2/(tmean**2*pw)
        nrefges3=k4*fact*qmean*pressure/(tmean*pw)
        nrefges(k,i)=nrefges1+nrefges2+nrefges3 ! refractivity N_i
        irefges(k)= one+(r1em6*nrefges(k,i))    ! index of refractivity n_i
        ref_rad(k)=irefges(k)*rges(k,i)         ! refractivity index-radius product x_i
 
print ('(a,i4,1p6e15.5)'),'QQQ NNN', k, nrefges(k,i), irefges(k),  ref_rad(k)

    end do 

!        Extending atmosphere above interface level nsig
         d_ref_rad=ref_rad(nsig)-ref_rad(nsig-1)
         do k=1,nsig_ext
            ref_rad(nsig+k)=ref_rad(nsig)+ k*d_ref_rad ! extended x_i
            nrefges(nsig+k,i)=nrefges(nsig+k-1,i)**2/nrefges(nsig+k-2,i) ! exended N_i
         end do

!        Lagrange coefficients
         ref_rad(0)=ref_rad(3)
         ref_rad(nsig_up+1)=ref_rad(nsig_up-2)
         do k=1,nsig_up
            call setq(q_w(:,k),ref_rad(k-1:k+1),3)
         enddo

!        Get refractivity index-radius and [d(ln(n))/dx] in new grid.
         intloop: do j=1,grids_dim
           ref_rad_s(j)=sqrt(grid_s(j)*grid_s(j)+tpdpres(i)*tpdpres(i)) !x_j
           xj(j,i)=ref_rad_s(j)
           hob_s=ref_rad_s(j)
           call grdcrd1(hob_s,ref_rad(1),nsig_up,1)
           dbend_loc(j,i)=hob_s  !location of x_j with respect to extended x_i
 

           if (hob_s < rsig_up) then  !obs inside the new grid
              ihob=hob_s

!             Compute refractivity and derivative at target points 
!             using Lagrange interpolators
              call slagdw(ref_rad(ihob-1:ihob+2),ref_rad_s(j),&
                   q_w(:,ihob),q_w(:,ihob+1),&
                   w4,dw4,4)
              if(ihob==1) then
                 w4(4)=w4(4)+w4(1); w4(1:3)=w4(2:4);w4(4)=zero
                 dw4(4)=dw4(4)+dw4(1);dw4(1:3)=dw4(2:4);dw4(4)=zero
                 ihob=ihob+1
              endif
              if(ihob==nsig_up-1) then
                 w4(1)=w4(1)+w4(4); w4(2:4)=w4(1:3);w4(1)=zero
                 dw4(1)=dw4(1)+dw4(4); dw4(2:4)=dw4(1:3);dw4(1)=zero
                 ihob=ihob-1
              endif
              ddnj(j)=dot_product(dw4,nrefges(ihob-1:ihob+2,i))!derivative (dN/dx)_j
              ddnj(j)=max(zero,abs(ddnj(j)))
           else
              print *,'obs outside new grid i,j=',i,j
           endif !obs in new grid
         end do intloop

!        bending angle (radians)
         dbend=ds*ddnj(1)/ref_rad_s(1)

print ('(a,2i4,1p6e16.6)'), 'QQQ DDD0', i,grids_dim,ds,tpdpres(i)
print ('(a,2i4,1p6e16.6)'),'QQQ DDD1', i,1,grid_s(1),ddnj(1),ref_rad_s(1),dbend



         do j=2,grids_dim
            ddbend=ds*ddnj(j)/ref_rad_s(j)
            dbend=dbend+two*ddbend

!QQQQQ start                                                                              
    if (abs(xlat-ron_lat) < 0.1 .and. abs(xlon-ron_lon) < 0.1 &    
        .and. ( abs(tpdpres(i)-6410037.) < 10. .or. &                                
                abs(tpdpres(i)-6420062.) < 10. .or. &                                
                abs(tpdpres(i)-6430184.) < 10. .or. &                                
                abs(tpdpres(i)-6440018.) < 10. ) ) then                              
                                                                                     
         ron_x1=sqrt(max(ref_rad_s(j)-tpdpres(i),1.))                                
         if (j<grids_dim) then                                                       
            ron_x2=sqrt(max(ref_rad_s(j+1)-tpdpres(i),1.))                           
         else                                                                        
            ron_x2=ron_x1                                                            
         endif                                                                       
         ron_x1=ron_x2-ron_x1                                                        
         ron_x3=sqrt(2.*tpdpres(i))*1.e-6*ddnj(j)*ron_x1  ! actually dbend on 1 side
                                                                                     
                                                                                     
         ron_x1=r1em6*tpdpres(i)*ddbend                                              
         ron_x2=r1em6*tpdpres(i)*dbend                                               
print ('(a,2i4,1p7e15.6)'),'QQQ DDD2', i,j,grid_s(j),ddnj(j),ref_rad_s(j), &
ddbend,ron_x3,ron_x1,ron_x2
    endif
!QQQQQQQQQ end

         end do
         dbend=r1em6*tpdpres(i)*dbend  

print ('(a,1p6e15.5)'), 'QQQ BBB', &
  dbend,nrefges(nsig+nsig_ext,i),ref_rad(nsig+nsig_ext)



         bend_angle=dbend !innovation vector

3000 continue
  end do loopoverobs1 ! end of loop over observations
  return
end subroutine setupbend


subroutine setq(q,x,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    setq
!
!   prgrmmr:
!
! abstract:      Precompute the N constant denominator factors of the N-point 
!                Lagrange polynomial interpolation formula.
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     X -    The N abscissae.
!     N -    The number of points involved.
!
!   output argument list:
!     Q -    The N denominator constants.
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

IMPLICIT NONE

  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
  real(r_kind), parameter :: zero=0._r_kind
  real(r_kind), parameter :: one=1._r_kind

INTEGER(i_kind)          ,INTENT(in   ) :: n
REAL(r_kind),DIMENSION(n),INTENT(  out) :: q
REAL(r_kind),DIMENSION(n),INTENT(in   ) :: x
!-----------------------------------------------------------------------------
INTEGER(i_kind)                      :: i,j
!=============================================================================
DO i=1,n
   q(i)=one
   DO j=1,n
      IF(j /= i)q(i)=q(i)/(x(i)-x(j))
   ENDDO
ENDDO
end subroutine setq



!============================================================================
subroutine lagdw(x,xt,q,w,dw,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    lagdw
!
!   prgrmmr:
!
! abstract:      Construct the Lagrange weights and their derivatives when 
!                target abscissa is known and denominators Q have already 
!                been precomputed
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     X   - Grid abscissae
!     XT  - Target abscissa
!     Q   - Q factors (denominators of the Lagrange weight formula)
!     N   - Number of grid points involved in the interpolation
!
!   output argument list:
!     W   - Lagrange weights
!     DW  - Derivatives, dW/dX, of Lagrange weights W
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
IMPLICIT NONE
  integer, parameter :: r_kind=8
  integer, parameter :: i_kind=4 
  real(r_kind), parameter :: zero=0._r_kind
  real(r_kind), parameter :: one=1._r_kind

INTEGER(i_kind)          ,INTENT(IN   ) :: n
REAL(r_kind)             ,INTENT(IN   ) :: xt
REAL(r_kind),DIMENSION(n),INTENT(IN   ) :: x,q
REAL(r_kind),DIMENSION(n),INTENT(  OUT) :: w,dw
!-----------------------------------------------------------------------------
REAL(r_kind),DIMENSION(n)            :: d,pa,pb,dpa,dpb
INTEGER(i_kind)                      :: j
!============================================================================
pa(1)=one
dpa(1)=zero
do j=1,n-1
   d(j)=xt-x(j)
   pa (j+1)=pa (j)*d(j)
   dpa(j+1)=dpa(j)*d(j)+pa(j)
enddo
d(n)=xt-x(n)

pb(n)=one
dpb(n)=zero
do j=n,2,-1
   pb (j-1)=pb (j)*d(j)
   dpb(j-1)=dpb(j)*d(j)+pb(j)
enddo
do j=1,n
   w (j)= pa(j)*pb (j)*q(j)
   dw(j)=(pa(j)*dpb(j)+dpa(j)*pb(j))*q(j)
enddo
end subroutine lagdw


!============================================================================
subroutine slagdw(x,xt,aq,bq,w,dw,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    slagdw
!
!   prgrmmr:
!
! abstract:      Construct weights and their derivatives for interpolation 
!                to a given target based on a linearly weighted mixture of 
!                the pair of Lagrange interpolators which omit the respective 
!                end points of the source nodes provided. The number of source 
!                points provided must be even and the nominal target interval 
!                is the unique central one. The linear weighting pair is 
!                determined by the relative location of the target within 
!                this central interval, or else the extreme values, 0 and 1, 
!                when target lies outside this interval.  The objective is to 
!                provide an interpolator whose derivative is continuous.
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     X     - Grid abscissae (N points)
!     XT    - Target abscissa
!     AQ,BQ - Q factors (denominators of the Lagrange weight formula for
!             the two respective consecutive subsets of N-1 grid points)
!     N     - Even number of grid points involved in the interpolation
!
!   output argument list:
!     W     - Final weights (N values)
!     DW    - Final derivatives, dW/dX, of weights W (N values)
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
IMPLICIT NONE
  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
  real(r_kind), parameter :: zero=0._r_kind
  real(r_kind), parameter :: one=1._r_kind

INTEGER(i_kind)            ,INTENT(IN   ) :: n
REAL(r_kind)               ,INTENT(IN   ) :: xt
REAL(r_kind),DIMENSION(n)  ,INTENT(IN   ) :: x
REAL(r_kind),DIMENSION(n-1),INTENT(IN   ) :: aq,bq
REAL(r_kind),DIMENSION(n)  ,INTENT(  OUT) :: w,dw
!-----------------------------------------------------------------------------
REAL(r_kind),DIMENSION(n)               :: aw,bw,daw,dbw
REAL(r_kind)                            :: xa,xb,dwb,wb
INTEGER(i_kind)                         :: na
!============================================================================
CALL lagdw(x(1:n-1),xt,aq,aw(1:n-1),daw(1:n-1),n-1)
CALL lagdw(x(2:n  ),xt,bq,bw(2:n  ),dbw(2:n  ),n-1)
aw(n)=zero
daw(n)=zero
bw(1)=zero
dbw(1)=zero
na=n/2
IF(na*2 /= n)STOP 'In slagdw; n must be even'
xa =x(na     )
xb =x(na+1)
dwb=one/(xb-xa)
wb =(xt-xa)*dwb
IF    (wb>one )THEN
   wb =one
   dwb=zero
ELSEIF(wb<zero)THEN
   wb =zero
   dwb=zero
ENDIF
bw =bw -aw
dbw=dbw-daw

w =aw +wb*bw
dw=daw+wb*dbw+dwb*bw
end subroutine slagdw

subroutine grdcrd1(d,x,nx,flg)

  implicit none

  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
  real(r_kind), parameter :: zero=0._r_kind
  real(r_kind), parameter :: one=1._r_kind

  integer(i_kind)           ,intent(in   ) :: nx
  integer(i_kind)           ,intent(in   ) :: flg
  real(r_kind)              ,intent(inout) :: d
  real(r_kind),dimension(nx),intent(in   ) :: x

  integer(i_kind) ix,isrchf

! Treat "normal" case in which nx>1
  if(nx>1) then
     if (flg==1) then
!       Case in which x is in increasing order
        if(d<=x(1)) then
           ix=1
        else
           ix=isrchf(nx-1,x,d,flg)-1
        end if
        if(ix==nx) ix=ix-1

     else if (flg==(-1)) then

!       Case in which x is in decreasing order     
        if(d>=x(1)) then
           ix=1
        else
           ix=isrchf(nx-1,x,d,flg)-1
        end if
     end if
     d=float(ix)+(d-x(ix))/(x(ix+1)-x(ix))

! Treat special case of nx=1                       
  elseif (nx==1) then
     d = one
  endif

  return
end subroutine grdcrd1


function isrchf(nx1,x,y,flg)

  integer, parameter :: r_kind=8 
  integer, parameter :: i_kind=4 
                                                                                    
  integer(i_kind)            ,intent(in   ) :: nx1       
  integer(i_kind)            ,intent(in   ) :: flg       
  real(r_kind)               ,intent(in   ) :: y         
  real(r_kind),dimension(nx1),intent(in   ) :: x         
  integer(i_kind) k                                   
  integer(i_kind):: isrchf                            

  if(flg==1) then                                     
     do k=1,nx1                  
        if(y<=x(k)) then         
           isrchf=k                                              
           go to 100                                             
        end if               
     end do                  
  else                                                   

     do k=1,nx1                  
        if(y>=x(k)) then         
           isrchf=k                    
           go to 100                   
        end if                
     end do                   
  end if                             
  isrchf=nx1+1                    
  if(nx1<=0) isrchf=0             
100 continue                                        
  return                                            
end function isrchf
