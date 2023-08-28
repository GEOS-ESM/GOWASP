   subroutine angles_ssmis (datetime,lat,lon,soza,solazi)
!
!  Determine solar zenith and azimuth angles from location and datetime
!  (Required for ssmis since this information is not in BUFR file)
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only: rkind1, rkind2
!
   implicit none
!
   real(rkind2), intent(in)  :: datetime(6)  ! y,m,d,h,m,s
   real(rkind2), intent(in)  :: lat,lon
   real(rkind2), intent(out) :: soza,solazi
!
   integer :: m
   integer :: mlen(12)
   integer :: jday           ! Julian day
   real(rkind1), parameter :: r60=1./60.
   real(rkind1) :: utc_hour  ! UTC hour
   real(rkind1) :: clon,clat
   real(rkind1) :: solar_elev,solar_azim
! 
   data mlen /31,28,31,30,31,30,31,31,30,31,30,31/
!
! Determine Julian day
   jday=nint(datetime(3))
   do m=1,nint(datetime(2))-1
     jday=jday+mlen(m)
   enddo
   if ((mod(nint(datetime(1)),4)==0).and.(nint(datetime(2)) > 2))  then
     jday=jday+1
   end if
!
! Determine time of day (in fractions of hours)
   utc_hour=datetime(4)+r60*(datetime(5)+r60*datetime(6))
!
!  Calculate solar zenith/azimuth angle                                         
   clat=lat
   clon=lon
   if (clon > 180.0_rkind1) then
     clon=clon-360.0_rkind1
   endif
!
   call zensun (jday,utc_hour,clat,clon,solar_elev,solar_azim)        
   soza= 90.-solar_elev  ! change from solar elev to solar zenith angle                              
   solazi=solar_azim
!
   end subroutine angles_ssmis 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine get_solar_elevation (datetime,lat,lon,solar_elev,solar_azim) 
!
!  Determine solar elevation and azimuth angles from location and datetime
!  as required for sat winds, similar to routine angles-ssmis, 
!
   use m_kinds, only: rkind1
!
   implicit none
!
   real(rkind1), intent(in)  :: datetime(6)  ! y,m,d,h,m,s
   real(rkind1), intent(in)  :: lat,lon
   real(rkind1), intent(out) :: solar_elev,solar_azim
!
   integer :: m
   integer :: mlen(12)
   integer :: jday           ! Julian day
   real(rkind1), parameter :: r60=1./60.
   real(rkind1) :: utc_hour  ! UTC hour
   real(rkind1) :: clon,clat
! 
   data mlen /31,28,31,30,31,30,31,31,30,31,30,31/
!
! Determine Julian day
   jday=nint(datetime(3))
   do m=1,nint(datetime(2))-1
     jday=jday+mlen(m)
   enddo
   if ((mod(nint(datetime(1)),4)==0).and.(nint(datetime(2)) > 2))  then
     jday=jday+1
   endif
!
! Determine time of day (in fractions of hours)
   utc_hour=datetime(4)+r60*(datetime(5)+r60*datetime(6))
!
!  Calculate solar zenith/azimuth angle                                         
   clat=lat
   clon=lon
   if (clon > 180.0_rkind1) then
     clon=clon-360.0_rkind1
   endif
!
   call zensun (jday,utc_hour,clat,clon,solar_elev,solar_azim)        
!
   end subroutine get_solar_elevation
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
subroutine zensun(jday,time,lat,lon,sun_elev,sun_azimuth)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  zensun         make sun elevation and sun azimuth angle
!
!   prgmmr: Paul Ricchiazzi org: Earth Space Research Group,UCSB  date: 1992-10-23
!
! abstract: 
!       Compute solar position information as a function of
!      geographic coordinates, date and time.
!
! program history log:
!   2005-10-21  kazumori - reformatted for GSI
!   2013-02-13  eliu    - bug fix for solar zenith calculation 
!   2014-06-04  sienkiewicz - take out doy**3 and replace with simple
!                               linear interpolation of tables
!
!   input argument list:
!     day -     Julian day (positive scalar or vector)
!               (spring equinox =  80)
!               (summer solstice= 171)
!               (fall equinox   = 266)
!               (winter solstice= 356)
!     time -    Universal Time in hours (scalar or vector)
!     lat  -    geographic latitude of point on earth's surface (degrees)
!     lon  -    geographic longitude of point on earth's surface (degrees)
!
!   output argument list:
!     sun_zenith  - solar zenith angle
!     sun_azimuth - solar azimuth angle
!
!   comments:
!
!
!     PROCEDURE:
!
!  1. Calculate the subsolar point latitude and longitude, based on
!     DAY and TIME. Since each year is 365.25 days long the exact
!     value of the declination angle changes from year to year.  For
!     precise values consult THE AMERICAN EPHEMERIS AND NAUTICAL
!     ALMANAC published yearly by the U.S. govt. printing office.  The
!     subsolar coordinates used in this code were provided by a
!     program written by Jeff Dozier.
!
!  2. Given the subsolar latitude and longitude, spherical geometry is
!     used to find the solar zenith, azimuth and flux multiplier.
!
!  eqt = equation of time (minutes)  ! solar longitude correction = -15*eqt
!  dec = declination angle (degrees) = solar latitude
!
! LOWTRAN v7 data (25 points)
!     The LOWTRAN solar position data is characterized by only 25 points.
!     This should predict the subsolar angles within one degree.  For
!     increased accuracy add more data points.
!
!nday=[   1.,    9.,   21.,   32.,   44.,   60.,  91.,  121.,  141.,  152.,$
!       160.,  172.,  182.,  190.,  202.,  213., 244.,  274.,  305.,  309.,$
!       325.,  335.,  343.,  355.,  366.]
!
!eqt=[ -3.23, -6.83,-11.17,-13.57,-14.33,-12.63, -4.2,  2.83,  3.57,  2.45,$
!       1.10, -1.42, -3.52, -4.93, -6.25, -6.28,-0.25, 10.02, 16.35, 16.38,$
!       14.3, 11.27,  8.02,  2.32, -3.23]
!
!dec=[-23.07,-22.22,-20.08,-17.32,-13.62, -7.88, 4.23, 14.83, 20.03, 21.95,$
!      22.87, 23.45, 23.17, 22.47, 20.63, 18.23, 8.58, -2.88,-14.18,-15.45,$
!     -19.75,-21.68,-22.75,-23.43,-23.07]
!
! Analemma information from Jeff Dozier
!     This data is characterized by 74 points
!
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!
  use m_kinds, only: rkind1  

  implicit none

  integer, parameter :: i_kind=4
  integer, parameter :: r_single=4
  integer, parameter :: r_kind=rkind1
  real(r_kind), parameter :: one=1.
  real(r_kind), parameter :: zero=0.
  integer, intent(in) :: jday
  real(r_kind) :: deg2rad,rad2deg,r60inv
  real(r_kind), intent(in   ) :: time,lat,lon
  real(r_kind), intent(  out) :: sun_elev,sun_azimuth
  integer(i_kind)   day
  integer(i_kind)   di
  real(r_kind)      ut,noon
  real(r_kind)      y(5),y2(5),x(2,5),Tx(5,2),xTx(2,2),aTx(5,2),det
  real(r_kind)      tt,eqtime,decang,latsun,lonsun,frac
  real(r_kind)      nday(74),eqt(74),dec(74)
  real(r_kind)      beta(2), beta2(2), a(2,2)
  real(r_kind)      t0,t1,p0,p1,zz,xx,yy

  data   nday/one,   6.0_r_kind,  11.0_r_kind,  16.0_r_kind,  21.0_r_kind,  26.0_r_kind,  &
                                    31.0_r_kind,  36.0_r_kind,  41.0_r_kind,  46.0_r_kind,&
             51.0_r_kind,  56.0_r_kind,  61.0_r_kind,  66.0_r_kind,  71.0_r_kind,  &
                76.0_r_kind,  81.0_r_kind,  86.0_r_kind,  91.0_r_kind,  96.0_r_kind,&
             101.0_r_kind, 106.0_r_kind, 111.0_r_kind, 116.0_r_kind, 121.0_r_kind, &
              126.0_r_kind, 131.0_r_kind, 136.0_r_kind, 141.0_r_kind, 146.0_r_kind,&
             151.0_r_kind, 156.0_r_kind, 161.0_r_kind, 166.0_r_kind, 171.0_r_kind, &
              176.0_r_kind, 181.0_r_kind, 186.0_r_kind, 191.0_r_kind, 196.0_r_kind,&
             201.0_r_kind, 206.0_r_kind, 211.0_r_kind, 216.0_r_kind, 221.0_r_kind, &
             226.0_r_kind, 231.0_r_kind, 236.0_r_kind, 241.0_r_kind, 246.0_r_kind,&
             251.0_r_kind, 256.0_r_kind, 261.0_r_kind, 266.0_r_kind, 271.0_r_kind, &
             276.0_r_kind, 281.0_r_kind, 286.0_r_kind, 291.0_r_kind, 296.0_r_kind,&
             301.0_r_kind, 306.0_r_kind, 311.0_r_kind, 316.0_r_kind, 321.0_r_kind, &
             326.0_r_kind, 331.0_r_kind, 336.0_r_kind, 341.0_r_kind, 346.0_r_kind,&
             351.0_r_kind, 356.0_r_kind, 361.0_r_kind, 366.0_r_kind/

  data  eqt/ -3.23_r_kind, -5.49_r_kind, -7.60_r_kind, -9.48_r_kind,-11.09_r_kind,&
             -12.39_r_kind,-13.34_r_kind,-13.95_r_kind,-14.23_r_kind,-14.19_r_kind,&
            -13.85_r_kind,-13.22_r_kind,-12.35_r_kind,-11.26_r_kind,-10.01_r_kind, &
            -8.64_r_kind, -7.18_r_kind, -5.67_r_kind, -4.16_r_kind, -2.69_r_kind,&
             -1.29_r_kind, -0.02_r_kind,  1.10_r_kind,  2.05_r_kind,  2.80_r_kind,  &
              3.33_r_kind,  3.63_r_kind,  3.68_r_kind,  3.49_r_kind,  3.09_r_kind,&
              2.48_r_kind,  1.71_r_kind,  0.79_r_kind, -0.24_r_kind, -1.33_r_kind, &
              -2.41_r_kind, -3.45_r_kind, -4.39_r_kind, -5.20_r_kind, -5.84_r_kind,&
             -6.28_r_kind, -6.49_r_kind, -6.44_r_kind, -6.15_r_kind, -5.60_r_kind, &
              -4.82_r_kind, -3.81_r_kind, -2.60_r_kind, -1.19_r_kind,  0.36_r_kind,&
              2.03_r_kind,  3.76_r_kind,  5.54_r_kind,  7.31_r_kind,  9.04_r_kind, &
             10.69_r_kind, 12.20_r_kind, 13.53_r_kind, 14.65_r_kind, 15.52_r_kind,&
             16.12_r_kind, 16.41_r_kind, 16.36_r_kind, 15.95_r_kind, 15.19_r_kind, &
             14.09_r_kind, 12.67_r_kind, 10.93_r_kind,  8.93_r_kind,  6.70_r_kind,&
              4.32_r_kind,  1.86_r_kind, -0.62_r_kind, -3.23_r_kind/

  data dec/ -23.06_r_kind,-22.57_r_kind,-21.91_r_kind,-21.06_r_kind,-20.05_r_kind,&
            -18.88_r_kind,-17.57_r_kind,-16.13_r_kind,-14.57_r_kind,-12.91_r_kind,&
            -11.16_r_kind, -9.34_r_kind, -7.46_r_kind, -5.54_r_kind, -3.59_r_kind,&
             -1.62_r_kind,  0.36_r_kind,  2.33_r_kind,  4.28_r_kind,  6.19_r_kind,&
              8.06_r_kind,  9.88_r_kind, 11.62_r_kind, 13.29_r_kind, 14.87_r_kind,&
             16.34_r_kind, 17.70_r_kind, 18.94_r_kind, 20.04_r_kind, 21.00_r_kind,&
             21.81_r_kind, 22.47_r_kind, 22.95_r_kind, 23.28_r_kind, 23.43_r_kind,&
             23.40_r_kind, 23.21_r_kind, 22.85_r_kind, 22.32_r_kind, 21.63_r_kind,&
             20.79_r_kind, 19.80_r_kind, 18.67_r_kind, 17.42_r_kind, 16.05_r_kind,&
             14.57_r_kind, 13.00_r_kind, 11.33_r_kind,  9.60_r_kind,  7.80_r_kind,&
              5.95_r_kind,  4.06_r_kind,  2.13_r_kind,  0.19_r_kind, -1.75_r_kind,&
             -3.69_r_kind, -5.62_r_kind, -7.51_r_kind, -9.36_r_kind,-11.16_r_kind,&
            -12.88_r_kind,-14.53_r_kind,-16.07_r_kind,-17.50_r_kind,-18.81_r_kind,&
            -19.98_r_kind,-20.99_r_kind,-21.85_r_kind,-22.52_r_kind,-23.02_r_kind,&
            -23.33_r_kind,-23.44_r_kind,-23.35_r_kind,-23.06_r_kind/

!
! compute the subsolar coordinates
! 
  rad2deg=45./atan(one)
  deg2rad=one/rad2deg
  r60inv=one/60. 
  day=jday  ! change of kind
!
  tt= mod(real((int(day)+time/24._r_kind-one)),365.25_r_single) + one  ! fractional day number
                                                                       ! with 12am 1jan = 1.
  do di = 1, 73
     if ((tt >= nday(di)) .and. (tt < nday(di+1))) exit
  end do


! using simple linear interpolation
  frac = (tt-nday(di))/(nday(di+1)-nday(di))
  eqtime = ((one-frac)*eqt(di)+frac*eqt(di+1))*r60inv
  decang = (one-frac)*dec(di) + frac*dec(di+1)

  latsun=decang

  ut=time
  noon=12._r_kind-lon/15._r_kind                      ! universal time of noon

  lonsun=-15._r_kind*(ut-12._r_kind+eqtime)

  t0=(90._r_kind-lat)*deg2rad
  t1=(90._r_kind-latsun)*deg2rad

  p0=lon*deg2rad
  p1=lonsun*deg2rad

  zz=cos(t0)*cos(t1)+sin(t0)*sin(t1)*cos(p1-p0)
  xx=sin(t1)*sin(p1-p0)
  yy=sin(t0)*cos(t1)-cos(t0)*sin(t1)*cos(p1-p0)

  sun_elev=90_r_kind-acos(zz)*rad2deg
  sun_azimuth=atan2(xx,yy)*rad2deg
  if (sun_azimuth < zero) sun_azimuth = sun_azimuth + 360.0_r_kind

  return
end subroutine zensun
