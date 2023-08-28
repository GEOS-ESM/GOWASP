module m_mympi
  use goss_mpeu,only : die,perr
  use goss_mpeu,only : IK => i_kind
  use goss_mpeu,only : LK => l_kind
  use goss_mpeu,only : SP => SP_kind
  use goss_mpeu,only : DP => DP_kind

  implicit none
  private
  public :: mympi_comm_world

  public :: mympi_success
  public :: mympi_setup
  public :: mympi_finalize

  public :: mympi_bcast
  public :: mympi_gather
  public :: mympi_scatter

  public :: mympi_barrier

  public :: mympi_comm_split
  public :: mympi_comm_free
  public :: die

!  integer,parameter :: IK=kind(1)            ! kind of default INTEGER
!  integer,parameter :: LK=kind(.true.)       ! kind of default LOGICAL
!  integer,parameter :: SP=4                  ! kind of default REAL
!  integer,parameter :: DP=8                  ! kind of DOUBLE PRECISION real


#ifndef NO_MPI
#define USE_MPI
#include "mpif.h"
#endif

  interface mympi_setup; module procedure &
    setup_		; end interface

  interface mympi_bcast; module procedure &
    bcast_i0,		&
    bcast_i1,		&
    bcast_i2,		&
    bcast_r0,		&
    bcast_r1,		&
    bcast_r2,		&
    bcast_c1		; end interface

  interface mympi_gather; module procedure &
    gather_i0,		&
    gather_d0,		&
    gather_d1,		&
    gather_d2,		&
    gather_r0,		&
    gather_r1,		&
    gather_r2		; end interface

  interface mympi_scatter; module procedure &
    scatter_i0,		&
    scatter_d0,		&
    scatter_d1,		&
    scatter_d2,		&
    scatter_r0,		&
    scatter_r1,		&
    scatter_r2		; end interface

  interface mympi_finalize; module procedure &
    finalize_		; end interface

  interface mympi_barrier; module procedure &
    barrier_		; end interface
  interface mympi_comm_split; module procedure &
	! call mympi_comm_split(comm,	! communicator
	!			color,	! ==.true.,key,newcomm,stat)
    comm_split_l	; end interface
  interface mympi_comm_free; module procedure &
    comm_free_		; end interface

  integer,parameter:: mympi_comm_world	= MPI_comm_world
  integer,parameter:: mympi_comm_null	= MPI_comm_null
  integer,parameter:: mympi_success	= MPI_success

  character(len=*), parameter :: myname='m_mympi'

#ifdef USE_MPI
  interface MP_type; module procedure &
    typeI0_, &
    typeL0_, &
    typeC0_, &
    typeR0_, &
    typeD0_, &
    typeI1_, &
    typeL1_, &
    typeC1_, &
    typeR1_, &
    typeD1_, &
    typeI2_, &
    typeL2_, &
    typeC2_, &
    typeR2_, &
    typeD2_, &
    typeI3_, &
    typeL3_, &
    typeC3_, &
    typeR3_, &
    typeD3_
  end interface
  interface MP_perr; module procedure &
    MP_perr_
  end interface
#endif
contains

subroutine setup_(comm,size,mype,stat)
  implicit none

!-- setup_() initializes MPI and gets the size of the
!-- "world" and the ID of a given process.

  integer,intent(in) :: comm	! the communicator (the Space)
  integer,intent(out) :: size	! the size of the comm
  integer,intent(out) :: mype	! who I am
  integer,optional,intent(out) :: stat	! check the status?

!-- local parameters
  character(len=*),parameter :: myname_=myname//"::setup_"
  integer :: ier

  if(present(stat)) stat=0

#ifdef USE_MPI
  call MPI_init(ier)
  	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_init()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call MPI_comm_rank(comm,mype,ier)
  	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  call MPI_comm_size(comm,size,ier)
  	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  myPE=0
  size=1
#endif
end subroutine setup_

subroutine finalize_(stat)
  implicit none
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::finalize_"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_finalize(ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MPI_finalize()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine finalize_

subroutine barrier_(comm,stat)
  implicit none
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat
!--
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::barrier_"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_barrier(comm,ier)
  	if(ier/=0) then
	  call MP_perr(myname_,'MPI_barrier()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine barrier_

subroutine bcast_c1(argv,root,comm,stat)
!-- args_bcast_() broadcasts a message of an array of CHARACTERs.
  implicit none

  character(len=*),dimension(:),intent(inout) :: argv	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_c1"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(argv,len(argv)*size(argv),MP_type(argv),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_c1

subroutine bcast_i0(ibufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of a scalar of INTERGER.
  implicit none

  integer,intent(inout) :: ibufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_i0"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(ibufr,1,MP_type(ibufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_i0

subroutine bcast_i1(ibufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of an array of INTEGERs.
  implicit none

  integer,dimension(:),intent(inout) :: ibufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_i1"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_i1

subroutine bcast_i2(ibufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of an array of INTEGERs.
  implicit none

  integer,dimension(:,:),intent(inout) :: ibufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_i2"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(ibufr,size(ibufr),MP_type(ibufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_i2

subroutine bcast_r0(rbufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of a scalar of REAL.
  implicit none

  real(SP),intent(inout) :: rbufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_r0"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(rbufr,1,MP_type(rbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_r0

subroutine bcast_r1(rbufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of an array of REALs.
  implicit none

  real(SP),dimension(:),intent(inout) :: rbufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_r1"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(rbufr,size(rbufr),MP_type(rbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_r1

subroutine bcast_r2(rbufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of a rank 2 array of REALs.
  implicit none

  real(SP),dimension(:,:),intent(inout) :: rbufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_r2"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(rbufr,size(rbufr),MP_type(rbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_r2

subroutine bcast_d0(rbufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of a scalar of DOUBLE PRECISION.
  implicit none

  real(DP),intent(inout) :: rbufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_d0"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(rbufr,1,MP_type(rbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_d0

subroutine bcast_d1(rbufr,root,comm,stat)
!-- args_bcast_() broadcasts a message of an array of DOUBLE PRECISIONs.
  implicit none

  real(DP),dimension(:),intent(inout) :: rbufr	! the message
  integer,intent(in) :: root		! the message source (PE).
  integer,intent(in) :: comm		! where the message is sent to.
  integer,optional,intent(out) :: stat	! status of the operation.

!-- local variable
  integer :: ier
  character(len=*),parameter :: myname_=myname//"::bcast_d1"
!-- handle the messages
  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_bcast(rbufr,size(rbufr),MP_type(rbufr),root,comm,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MPI_bcast()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#endif
end subroutine bcast_d1

subroutine gather_i0(sbufr,rbufr,root,comm,stat)
  implicit none
  integer,intent(in) :: sbufr
  integer,dimension(:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::gather_i0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_gather(sbufr,1,MP_type(sbufr),	&
  		  rbufr,1,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(1)=sbufr
#endif
end subroutine gather_i0

subroutine gather_r0(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),intent(in) :: sbufr
  real(SP),dimension(:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::gather_r0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_gather(sbufr,1,MP_type(sbufr),	&
  		  rbufr,1,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(1)=sbufr
#endif
end subroutine gather_r0

subroutine gather_r1(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),dimension(:),intent(in) :: sbufr
  real(SP),dimension(:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::gather_r1"

  if(present(stat)) stat=0
  lsize=size(sbufr)
#ifdef USE_MPI
  call MPI_gather(sbufr,lsize,MP_type(sbufr),	&
  		  rbufr,lsize,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,1)=sbufr(:)
#endif
end subroutine gather_r1

subroutine gather_r2(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),dimension(:,:),intent(in) :: sbufr
  real(SP),dimension(:,:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::gather_r2"

  if(present(stat)) stat=0
  lsize=size(sbufr)
#ifdef USE_MPI
  call MPI_gather(sbufr,lsize,MP_type(sbufr),	&
  		  rbufr,lsize,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,:,1)=sbufr(:,:)
#endif
end subroutine gather_r2

subroutine scatter_i0(sbufr,rbufr,root,comm,stat)
  implicit none
  integer,dimension(:),intent(in) :: sbufr
  integer,intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::scatter_i0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_scatter(sbufr,1,MP_type(sbufr),	&
  		   rbufr,1,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr=sbufr(1)
#endif
end subroutine scatter_i0

subroutine scatter_r0(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),dimension(:),intent(in) :: sbufr
  real(SP),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::scatter_r0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_scatter(sbufr,1,MP_type(sbufr),	&
  		   rbufr,1,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr=sbufr(1)
#endif
end subroutine scatter_r0

subroutine scatter_r1(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),dimension(:,:),intent(in) :: sbufr
  real(SP),dimension(:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::scatter_r1"

  if(present(stat)) stat=0
  lsize=size(rbufr)
#ifdef USE_MPI
  call MPI_scatter(sbufr,lsize,MP_type(sbufr),	&
  		   rbufr,lsize,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:)=sbufr(:,1)
#endif
end subroutine scatter_r1

subroutine scatter_r2(sbufr,rbufr,root,comm,stat)
  implicit none
  real(SP),dimension(:,:,:),intent(in) :: sbufr
  real(SP),dimension(:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::scatter_r2"

  if(present(stat)) stat=0
  lsize=size(rbufr)
#ifdef USE_MPI
  call MPI_scatter(sbufr,lsize,MP_type(sbufr),	&
  		   rbufr,lsize,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,:)=sbufr(:,:,1)
#endif
end subroutine scatter_r2

subroutine gather_d0(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),intent(in) :: sbufr
  real(DP),dimension(:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::gather_d0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_gather(sbufr,1,MP_type(sbufr),	&
  		  rbufr,1,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(1)=sbufr
#endif
end subroutine gather_d0

subroutine gather_d1(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),dimension(:),intent(in) :: sbufr
  real(DP),dimension(:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::gather_d1"

  if(present(stat)) stat=0
  lsize=size(sbufr)
#ifdef USE_MPI
  call MPI_gather(sbufr,lsize,MP_type(sbufr),	&
  		  rbufr,lsize,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,1)=sbufr(:)
#endif
end subroutine gather_d1

subroutine gather_d2(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),dimension(:,:),intent(in) :: sbufr
  real(DP),dimension(:,:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::gather_d2"

  if(present(stat)) stat=0
  lsize=size(sbufr)
#ifdef USE_MPI
  call MPI_gather(sbufr,lsize,MP_type(sbufr),	&
  		  rbufr,lsize,MP_type(rbufr),	&
		  root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_gather()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,:,1)=sbufr(:,:)
#endif
end subroutine gather_d2

subroutine scatter_d0(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),dimension(:),intent(in) :: sbufr
  real(DP),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::scatter_d0"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_scatter(sbufr,1,MP_type(sbufr),	&
  		   rbufr,1,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr=sbufr(1)
#endif
end subroutine scatter_d0

subroutine scatter_d1(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),dimension(:,:),intent(in) :: sbufr
  real(DP),dimension(:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::scatter_d1"

  if(present(stat)) stat=0
  lsize=size(rbufr)
#ifdef USE_MPI
  call MPI_scatter(sbufr,lsize,MP_type(sbufr),	&
  		   rbufr,lsize,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:)=sbufr(:,1)
#endif
end subroutine scatter_d1

subroutine scatter_d2(sbufr,rbufr,root,comm,stat)
  implicit none
  real(DP),dimension(:,:,:),intent(in) :: sbufr
  real(DP),dimension(:,:),intent(out) :: rbufr
  integer,intent(in) :: root
  integer,intent(in) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier,lsize
  character(len=*),parameter :: myname_=myname//"::scatter_d2"

  if(present(stat)) stat=0
  lsize=size(rbufr)
#ifdef USE_MPI
  call MPI_scatter(sbufr,lsize,MP_type(sbufr),	&
  		   rbufr,lsize,MP_type(rbufr),	&
		   root,comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  rbufr(:,:)=sbufr(:,:,1)
#endif
end subroutine scatter_d2

subroutine comm_split_l(comm,color,newcomm,rank,stat)
  implicit none
  integer,intent(in) :: comm	! a given communicator
  logical,intent(in) :: color	! split according to color is .true. or .false.
  integer,intent(out) :: newcomm
				! newcomm is the new communicator after this
				! spliting.  In fact, there are two newcomm
				! created.  One contains all PEs with .true.
				! (color).  Another contains all PEs with
				! .false.

  integer,optional,intent(in ) :: rank
  	! If rank is not present, myPE is used for rank assignment control.
  integer,optional,intent(out) :: stat
  	! If stat is not present, die() is invoked upon an error.

  integer :: ier,icolor,irank
  character(len=*),parameter :: myname_=myname//"::comm_split_l"

  if(present(stat)) stat=0

  if(present(rank)) then
    irank=rank
  else
    call MPI_comm_rank(comm,irank,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

#ifdef USE_MPI
  icolor=0
  if(color) icolor=1
  call MPI_comm_split(comm,icolor,irank,newcomm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_comm_split()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  newcomm=comm
#endif
end subroutine comm_split_l

subroutine comm_free_(comm,stat)
  implicit none
  integer,intent(inout) :: comm
  integer,optional,intent(out) :: stat

  integer :: ier
  character(len=*),parameter :: myname_=myname//"::comm_free_"

  if(present(stat)) stat=0
#ifdef USE_MPI
  call MPI_comm_free(comm,ier)
	if(ier/=MYMPI_SUCCESS) then
	  call MP_perr(myname_,'MPI_comm_free()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
#else
  comm=mympi_comm_null
#endif
end subroutine comm_free_

#ifdef USE_MPI
function typeI0_(m) result(mytype)
  implicit none
  integer(IK),intent(in):: m
  integer:: mytype
  mytype=MPI_integer
  end function typeI0_
function typeL0_(m) result(mytype)
  implicit none
  logical(LK),intent(in):: m
  integer:: mytype
  mytype=MPI_logical
  end function typeL0_
function typeC0_(m) result(mytype)
  implicit none
  character(len=*),intent(in):: m
  integer:: mytype
  mytype=MPI_character
  end function typeC0_
function typeR0_(m) result(mytype)
  implicit none
  real(SP),intent(in):: m
  integer:: mytype
  mytype=MPI_real
  end function typeR0_
function typeD0_(m) result(mytype)
  implicit none
  real(DP),intent(in):: m
  integer:: mytype
  mytype=MPI_double_precision
  end function typeD0_

function typeI1_(m) result(mytype)
  implicit none
  integer(IK),dimension(:),intent(in):: m
  integer:: mytype
  mytype=MPI_integer
  end function typeI1_
function typeL1_(m) result(mytype)
  implicit none
  logical(LK),dimension(:),intent(in):: m
  integer:: mytype
  mytype=MPI_logical
  end function typeL1_
function typeC1_(m) result(mytype)
  implicit none
  character(len=*),dimension(:),intent(in):: m
  integer:: mytype
  mytype=MPI_character
  end function typeC1_
function typeR1_(m) result(mytype)
  implicit none
  real(SP),dimension(:),intent(in):: m
  integer:: mytype
  mytype=MPI_real
  end function typeR1_
function typeD1_(m) result(mytype)
  implicit none
  real(DP),dimension(:),intent(in):: m
  integer:: mytype
  mytype=MPI_double_precision
  end function typeD1_

function typeI2_(m) result(mytype)
  implicit none
  integer(IK),dimension(:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_integer
  end function typeI2_
function typeL2_(m) result(mytype)
  implicit none
  logical(LK),dimension(:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_logical
  end function typeL2_
function typeC2_(m) result(mytype)
  implicit none
  character(len=*),dimension(:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_character
  end function typeC2_
function typeR2_(m) result(mytype)
  implicit none
  real(SP),dimension(:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_real
  end function typeR2_
function typeD2_(m) result(mytype)
  implicit none
  real(DP),dimension(:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_double_precision
  end function typeD2_

function typeI3_(m) result(mytype)
  implicit none
  integer(IK),dimension(:,:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_integer
  end function typeI3_
function typeL3_(m) result(mytype)
  implicit none
  logical(LK),dimension(:,:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_logical
  end function typeL3_
function typeC3_(m) result(mytype)
  implicit none
  character(len=*),dimension(:,:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_character
  end function typeC3_
function typeR3_(m) result(mytype)
  implicit none
  real(SP),dimension(:,:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_real
  end function typeR3_
function typeD3_(m) result(mytype)
  implicit none
  real(DP),dimension(:,:,:),intent(in):: m
  integer:: mytype
  mytype=MPI_double_precision
  end function typeD3_

subroutine MP_perr_(proc,MP_proc,ierror)
  implicit none
  character(len=*),intent(in) :: proc
  character(len=*),intent(in) :: MP_proc
  integer,intent(in) :: ierror

  character(len=*),parameter :: myname_=myname//'::MP_perr_'
  character(len=MPI_MAX_ERROR_STRING) :: estr
  integer :: ln,ier

  call perr(proc,MP_proc//', ierror =',ierror)
  call MPI_error_string(ierror,estr,ln,ier)
  if(ier == 0 .and. ln>0) &
    call perr(proc,MP_proc//', "'//estr(1:ln)//'"')
  end subroutine MP_perr_
#endif
end module m_mympi
