   program bkg_fields
!
! read selected bkg or ana fields
!
   implicit none
!
   integer, parameter :: nlons=576
   integer, parameter :: nlats=361
   integer, parameter :: nlevs=72
   integer, parameter :: nfields2d=5
   integer, parameter :: nfields3d=0
   integer, parameter :: nfields=nfields2d+nfields3d
   integer :: ier
   integer :: nf,nx,i,j
   integer :: argc
   integer(4) :: iargc
!  
   real :: xi,xj,s1,s2,di,dj
   real(4), allocatable :: fields2d(:,:,:)
   real(4), allocatable :: fields3d(:,:,:,:)
!
   character(len=240) :: file_name
   character(len=10)   :: field_name(nfields)
!
! set names of all 2d fields, followed by all 3d 
   field_name(1)='frlake'
   field_name(2)='frland'
   field_name(3)='frlandice'
   field_name(4)='frocean'
   field_name(5)='frseaice'
!
   if (nfields2d > 0) then
     allocate (fields2d(nlons,nlats,nfields2d))
   endif
   if (nfields3d > 0) then
     allocate (fields3d(nlons,nlats,nlevs,nfields3d))
   endif
!
! Get program arguments
   argc=iargc()
   if (argc /= 1) then
     print *,' usage must be: prog.x filename'
     stop
   endif
   call GetArg(1_4, file_name)
!
! Loop over input 2d fields
  if (nfields2d > 0) then 
    do nf=1,nfields2d
      call read_nc4_2dfield (nlons,nlats,file_name,field_name(nf), &
                    fields2d(:,:,nf),.true.,.true.,ier)
    enddo
  endif 
!
! Loop over input 3d fields
  if (nfields3d > 0) then  
    do nf=1,nfields3d
      nx=nf+nfields2d
      call read_nc4_2dfield (nlons,nlats,nlevs,file_name,field_name(nx), &
                    fields3d(:,:,:,nf),.true.,.true.,ier)
    enddo
  endif
!
! Output some data
  xi=180.-85.
  di=360./nlons
  dj=180./(nlats-1)
  i=1+xi/di
  do j=1,nlats
    xj=-90.+(j-1)*dj
    s1=fields2d(i,j,1)+fields2d(i,j,2)+fields2d(i,j,3)+fields2d(i,j,4)
    s2=s1+fields2d(i,j,5)
    print ('(2i4,2f10.5,8f8.5)'),i,j,xi,xj,fields2d(i,j,:),s1,s2
  enddo
!
  print *,'Program Done'
!
  end program bkg_fields
