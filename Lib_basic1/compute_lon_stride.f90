  subroutine compute_lon_stride (imax,jmax,j,stride_eq,stride)
!
! Compute the number of gridpoints at the requested latitude  
! corresponding to the requested stride at the equator accouting 
! for convergence of the meridians.
! Do not allow greater than imax/2 since otherwise points are closer than
! that given the periodic nature of the points, except at the poles where 
! the stride is set to imax.
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
  use m_kinds, only : rkind1
!
  implicit none
! 
  integer, intent(in)  :: imax,jmax  ! numbers of lon, lat points
  integer, intent(in)  :: j          ! index of lat
  integer, intent(in)  :: stride_eq  ! requested stride at equator
  integer, intent(out) :: stride     
!
  real(rkind1) :: dlon    
  real(rkind1) :: dlat   ! spacing of lats in radians  
  real(rkind1) :: rlat   ! lat corresponding to j (radians)
  real(rkind1) :: pifac  ! pi
! 
! Compute distance dlon of requested stride in units of number of 
! longitude points at the requested latitude
  if (j ==1 .or. j==jmax) then  ! consider only 1 lon at poles 
    stride=imax
  else
    pifac=4*atan(1.)
    dlat=pifac/real(jmax-1)
    rlat=real(j-1)*dlat-0.5*pifac
    dlon=real(stride_eq)/cos(rlat)
    stride=min(imax/2,max(stride_eq,int(dlon)))
  endif  
!
  end subroutine compute_lon_stride 
