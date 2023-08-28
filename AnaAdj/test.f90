  program test
!
! Compute test results for enorm.f90 by computing exact and estimated norms
! using output from that program in test mode
!
  logical :: ltv 
  integer, parameter :: nfields=5
  integer, parameter :: nvalues=5
  real(4) :: f(nvalues,3,nfields,2)
  real(4) :: x(nvalues,nfields,2)
  real(4) :: t(nvalues,2,2)
  real(4) :: scales(nfields)
  character(len=12) :: fnames(3,nfields+1,2)    
  character(len=6)  :: cnames(3,nfields+1,2)
  character(len=3)  :: caorb
  character(len=20) :: cdum(3)
!    
! read test output
  print *,'Test program'
  open (10,file='file_prnt_test.txt')
  print *,'file_opened'
!
  do nt=1,2  ! loop over ana then bkg
    read (10,'(a6,l8,a8,a3)') cdum(1),ltv,cdum(2),caorb 
    print *,cdum(1),ltv,cdum(2),caorb
    read (10,'(a12,5e15.4)') cdum(3),scales(:)
    print *,cdum(3),scales(:)
!
    do nf=1,nfields  ! loop over fields p,q,t,u,v 
      do n=1,3   ! loop over ana, nr, diff
        read (10,'(a6,a12,5e14.5)') &
               cnames(n,nf,nt),fnames(n,nf,nt),f(:,n,nf,nt)
        print ('(3i4,2(1x,a))'),n,nf,nt,cnames(n,nf,nt),fnames(n,nf,nt)
      enddo 
    enddo 
    if (ltv) then ! read T derived from TV
      nf=nfields+1
      do n=1,2
        read (10,'(a6,a12,5e14.5)') cnames(n,nf,nt),fnames(n,nf,nt),t(:,n,nt)
        print ('(3i4,2(1x,a))'),n,nf,nt,cnames(n,nf,nt),fnames(n,nf,nt)
      enddo 
    endif    
!
  enddo 
!
! compute norms for each separate data field
  print *,' '
  do nf=1,nfields
    do n=1,nvalues
! Compute exact difference in norms
      amt=f(n,1,nf,1)-f(n,2,nf,1)
      bmt=f(n,1,nf,2)-f(n,2,nf,2)
      x(n,nf,1)=scales(nf)*(amt*amt-bmt*bmt)
!
! If Tv control variable, then recompute ET from T rather than Tv 
      if (ltv .and. nf==3) then  
        amt=t(n,1,1)-t(n,2,1)
        bmt=t(n,1,2)-t(n,2,2)
        x(n,nf,1)=scales(nf)*(amt*amt-bmt*bmt)  
      endif
!
! Comput estimate of difference in norms using adjoint fields
! Note use of factor 0.5
      amb=f(n,1,nf,1)-f(n,1,nf,2)
      x(n,nf,2)=0.5*amb*(f(n,3,nf,1)+f(n,3,nf,2)) 
    enddo
!
! Print reults for true and adj estimates
    print *,' '
    print ('(i2,1x,a,5e14.5)'),nf,fnames(1,nf,1),x(:,nf,1) 
    print ('(i2,1x,a,5e14.5)'),nf,fnames(3,nf,1),x(:,nf,2) 
  enddo
!
! print out results for EQ and ET summed.  Only this sum should agree
! when Tv is used as control variable but ET=cT*T
  nf=nfields+1
  do n=1,nvalues
    if (ltv) then ! adjust norm to what should be eact comparison
      z=(1.+.622*f(n,1,2,1))*(1.+.622*f(n,1,2,2))
      fnames(1,nf,1)='sum EQ+Z*ET '
    else
      z=1.
      fnames(1,nf,1)='sum EQ+ET   '
    endif 
    fnames(3,nf,1)='adjest EQ+ET'
    x(n,3,1)=x(n,2,1)+x(n,3,1)*z
    x(n,3,2)=x(n,2,2)+x(n,3,2)
  enddo
  print *,' '
  print ('(i2,1x,a,5e14.5)'),nf,fnames(1,nf,1),x(:,3,1) 
  print ('(i2,1x,a,5e14.5)'),nf,fnames(3,nf,1),x(:,3,2) 
!
  end program test
