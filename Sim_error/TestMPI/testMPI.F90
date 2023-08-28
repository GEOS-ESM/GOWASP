      program testMPI
!         use MAPL_ShmemMod    ! The SHMEM infrastructure

   implicit none
   include "mpif.h"
!
   integer, parameter :: r4=4
   integer :: myid          ! processor id number 0, ... ,npet                             
   integer :: npet          ! number of processors used        
   integer :: ierr
   integer :: i,j,k,m,n
   integer :: ks,itag
   integer(8) :: i_random(1)
   real(r4), allocatable :: f(:,:,:,:)
   real(r4) :: fld(5,2,2)
   real(8) :: x
   integer :: iseed
   integer :: stat(MPI_STATUS_SIZE)

    !  Initialize MPI                                                                                                          
!  --------------                                                                                                          
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,npet,ierr)
   print *,'Starting MPI on ',npet, ' processors: myid=',myid
!
   
   if (myid==0) then
     allocate (f(5,2,3,1:2))
     m=0
     do k=1,5
       do n=1,2
         do i=1,3 
           do j=1,2
             m=m+1
             f(k,n,i,j)=m
           enddo
         enddo
       enddo
     enddo
   endif 
   if (myid==1) then
     allocate (f(5,2,3,3:4))
     m=0
     do k=1,5
       do n=1,2
         do i=1,3 
           do j=3,4
             m=m-1
             f(k,n,i,j)=m
           enddo
         enddo
       enddo
     enddo
     print *,'SET ',myid,f(2:3,1,2,3)
   endif 
!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)                                                                                                                          
! Set seed for random number generator to datetime                                                                         
   i_random(1)=111
 call random_seed (put=i_random(1:1))
 call random_number(x)
 print *,'R1 ',myid,x
if (myid==1) then
  call random_number(x)
 print *,'R2 ',myid,x
endif
  call random_number(x)
 print *,'R3 ',myid,x
 itag=1
 if (myid==0) then
  iseed=int(1111+1000.*x)
  call MPI_SEND (iseed,1,MPI_INTEGER,1,itag,MPI_COMM_WORLD,ierr) 
 endif  
 if (myid==1) then
     call MPI_RECV  (iseed,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD,stat,ierr)
 endif
 print *,'IS ',myid,iseed
   i_random(1)=iseed
  call random_seed (put=i_random(1:1))
  call random_number(x)
 print *,'R4 ',myid,x
 
!   
   ks=2 
   itag=2


   if (myid==1) then 
     call MPI_SEND (f(2,1,2,3),ks,MPI_REAL4,0,itag,MPI_COMM_WORLD,ierr)
     print *,'C1 ',ierr
   endif
    
    if (myid==0) then
     call MPI_RECV  (fld(1,1,1),ks,MPI_REAL4,1,itag,MPI_COMM_WORLD,stat,ierr)
     print *,'C2 ',ierr
   endif
!
   print *,'FLD=',myid,itag,fld(1:ks,1,1) 
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
   print *,'END ',myid
!
   end program testMPI
