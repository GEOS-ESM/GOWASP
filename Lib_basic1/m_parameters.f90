   module m_parameters
!
! Set some physical parameters and methematical constants
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1,rkind2
!
   real(rkind1), parameter :: grav=9.80665     ! accel. gravity m/s**2
   real(rkind1), parameter :: earthr=6.371e6   ! earth's radius in m
   real(rkind1), parameter :: ratio4R=0.622    ! weight water vap / dry air
   real(rkind1), parameter :: Rdryair=287.06
   real(rkind1), parameter :: Cpdryair=1005.7
   real(rkind1), parameter :: Xkappa=Rdryair/Cpdryair
   real(rkind1), parameter :: XkappaI=1./Xkappa
   real(rkind1), parameter :: pi_k1=3.14159265
   real(rkind2), parameter :: pi_k2=3.14159265358979
   real(rkind1), parameter :: pifac_k1= pi_k1/180._rkind1  ! degrees to radians
   real(rkind2), parameter :: pifac_k2= pi_k2/180._rkind2  ! degrees to radians
   real(rkind1), parameter :: rpifac_k1=1._rkind1/pifac_k1
   real(rkind2), parameter :: rpifac_k2=1._rkind2/pifac_k2
   real(rkind1), parameter :: earth_c=2.0_rkind1*pi_k1*earthr*0.001 ! C (km)
!
   end module m_parameters
