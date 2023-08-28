    subroutine read_akbk (ndim,nlevs,iunit,fname,lprint,akbk,ier) 
!
! Read file of a(k) and b(k) values for defining eta coordinates
! (pressures at interfaces of NR data-layers)
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
    implicit none
!
    logical, intent(in)  :: lprint
    integer, parameter :: rkind1=4
    integer, intent(in)  :: ndim   
    integer, intent(in)  :: nlevs  ! number of data levels for NR 3D fields
    integer, intent(in)  :: iunit
    integer, intent(out) :: ier
    real(rkind1), intent(out) :: akbk(ndim,2)
    character(len=*), intent(in) :: fname
!
    integer :: n, n1
    integer :: nlevsp1
    integer :: iformat, nvalues
    character(len=100) :: cdum
!
! read ak, bk for vertical eta grid structure
    nlevsp1=nlevs+1
    open (iunit,file=trim(fname),form='formatted')
    read (iunit,*) cdum,iformat,nvalues
!
    if (nvalues == nlevsp1) then
      read (iunit,'(a1)') cdum   ! skip descriptor record  
      do n=1,nvalues
       read (iunit,*) n1,akbk(n,1:2)
      enddo
!
      close (iunit)
      if (lprint) then
        print ('(a,i3,2a)'),'File of akbk read:  values=',nvalues, &
                            '  file name=',trim(fname)
      endif
      ier=0
    else
      n=nlevs+1
      if (lprint) then
        print *,'ERROR: PROBLEM WITH CHOICE OF eta_akbk FILE:'     
        print ('(a,i3)'),'number of values on chosenfile =',nvalues
        print ('(a,i3)'),'nlevs + 1 of values required for chosen NR =',nlevsp1
      endif
      ier=100
    endif
! 
    end subroutine read_akbk
!
