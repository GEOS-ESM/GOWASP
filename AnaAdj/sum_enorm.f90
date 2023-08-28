   program sum_enorm
!   
! Sum norm values produced by enorms.f90
!
   implicit none
!
   integer, parameter :: ifield_max=20
   integer, parameter :: iunit=10
   integer :: ntimes
   integer :: ifields
   integer :: nt, nf, nk
!
   real(4) :: enorm(ifield_max)
   real(4) :: enorm_sum(ifield_max,2)
!
   character(len=6)   :: cdum
   character(len=44)  :: cdum2
   character(len=3)   :: caorb
   character(len=8)   :: cnorm
   character(len=4)   :: ctimes
   character(len=4)   :: field_name(ifield_max)   
   character(len=220)   :: file_in
!
! read arguments
   call GetArg( 1_4, file_in)
   call GetArg( 2_4, ctimes)
   read (ctimes,'(i4)') ntimes
!
   enorm_sum(:,:)=0.
!
! write norms to separate file
   open (iunit,file=trim(file_in))
   do nt=1,ntimes*2
     nk=mod(nt+1,2)+1
     read (iunit,'(a44,i2)') cdum2,ifields
     read (iunit,'(a6,a4,e13.4)') cdum,field_name(1),enorm(1)      
     do nf=2,ifields+1
       read (iunit,'(a6,a4,e13.4)') cdum,field_name(nf),enorm(nf)      
     enddo
     enorm_sum(:,nk)=enorm_sum(:,nk)+enorm(:)
   enddo
   close (iunit)
!
   enorm_sum(:,:)=enorm_sum(:,:)/ntimes
!
   print ('(a,i3,a)'),'mean for ana over ',ntimes,' times'
   do nf=1,ifields+1
     print ('(a,1p1e13.4)'),field_name(nf),enorm_sum(nf,1)
   enddo
   print ('(a,i3,a)'),'mean for bkg over ',ntimes,' times'
   do nf=1,ifields+1
     print ('(a,1p1e13.4)'),field_name(nf),enorm_sum(nf,2)
     enorm_sum(nf,1)=enorm_sum(nf,1)-enorm_sum(nf,2)
   enddo
   print ('(a,i3,a)'),'mean for ana-bkg over ',ntimes,' times'
   do nf=1,ifields+1
     print ('(a,1p1e13.4)'),field_name(nf),enorm_sum(nf,1)
   enddo
!
   end program sum_enorm
