   program compute_dx
!   
! Detmine distributions of closest neighbors to each obs for each obs type. 
!
! For each kx,ks type and subtype, determine a histogram of closest separation distances
! between observations of the same type,subtype and same NR time interval. 
! Do not distinguish between different p-layers. 
! For each observation within each such subset, consider all other obs within the same subset
! to find the distance to its nearest neighbor. Bin that distance (here, 5 km bin widths).
! Then count the number of values that fall within each bin.  Sum over many NR times.
! Then determine the fraction of values that fall within each bin compared with the total number
! of obs considered for this subtype. Print that table.
! 

   implicit none
!
   integer, parameter :: rkind1=4
   integer, parameter :: nlevs=1
   integer, parameter :: ncmax=20000
   integer, parameter :: nend=999999
   integer, parameter :: nbins=40
   integer, parameter :: iunit=21
!
   integer :: nkmax
   integer :: ntmax
   integer :: kid
   integer :: n6h
   integer :: n6h_sub
   integer :: i,n,n1,n2
   integer :: nt,nk
   integer :: ibsum
   integer :: max_6hr,max_n6h_sub
   integer, allocatable :: nc(:,:)
   integer, allocatable :: ibins(:,:,:)
   real(rkind1) :: dmin
   real(rkind1) :: earthr2,d,d1,d2,d3
   real(rkind1) :: rbsum
   real(rkind1) :: plev,xlat,xlon
   real(rkind1) :: read3(3)
   real(rkind1) :: rbins(nbins)
   real(rkind1), allocatable :: data3(:,:,:,:)
!
   real(rkind1), parameter :: dbin=5.
   real(rkind1), parameter :: earthr=6.371e3
   real(rkind1), parameter :: pifac=3.14159265/180.
!
   character(len=240) :: file_in  ! file output by count_amv.f90
   character(len=240) :: file_out ! distibution of dx values 
!
   earthr2=2.*earthr
!
   call GetArg(1_4, file_in)  ! input file of obs locations from count_amv 
   call GetArg(2_4, file_out) ! output file of dx min distributions
!
   open (iunit,file=trim(file_in),form='formatted')
   print *,trim(file_in),' opened'
   read (iunit,'(3i4)') max_6hr,max_n6h_sub,nkmax
   ntmax=max_6hr*max_n6h_sub
   allocate (data3(3,ncmax,ntmax,nkmax))
   allocate (nc(ntmax,nkmax))
   allocate (ibins(nbins,ntmax,nkmax))

   print *,'data3 allocated with ncmax,ntmax,nkmax=',ncmax,ntmax,nkmax
!
! Read all data into an array
   nc(:,:)=0
   do n=1,nend
     read (iunit,'(3i4,3f12.3)') kid,n6h,n6h_sub,read3
     if (kid == 999) then  ! indicates eof
       exit
     else
       nt=n6h_sub+max_n6h_sub*(n6h-1)   ! time counter
       nk=kid                           ! index for kx type and subtype pair
       nc(nt,nk)=min(nc(nt,nk)+1,ncmax) ! counter of obs for each nt, kid
       data3(1,nc(nt,nk),nt,nk)=read3(1)         ! plev of obs
       data3(2,nc(nt,nk),nt,nk)=read3(2)*pifac   ! lat of obs (radians)
       data3(3,nc(nt,nk),nt,nk)=read3(3)*pifac   ! lon of obs (radians)
     endif
   enddo
   close (iunit)
!
   ibins(:,:,:)=0
!
   do nk=1,nkmax
   do nt=1,ntmax
     print *,'considering obs for nk,nt,ncnt=',nk,nt,nc(nt,nk)
     if (nc(nt,nk) > 1) then
       do n1=1,nc(nt,nk)      ! consider each obs of this type and time
         dmin=1.e10
         do n2=1,nc(nt,nk)    ! consider all pairs with obs n1 of same nk,nt
           if (n1 /= n2) then ! do not pair with self
             d1=sin(0.5*(data3(2,n1,nt,nk)-data3(2,n2,nt,nk))) ! sin 1/2 dlat
             d2=sin(0.5*(data3(3,n1,nt,nk)-data3(3,n2,nt,nk))) ! sin 1/2 dlon
             d3=cos(data3(2,n1,nt,nk))*cos(data3(2,n2,nt,nk))
             d=earthr2*asin(sqrt(d1*d1+d2*d2*d3))  ! d between pair (km)
             dmin=min(d,dmin)  ! d of closest obs to this obs 
           endif
         enddo
         i=min(1+int(dmin/dbin),nbins)
         ibins(i,nt,nk)=ibins(i,nt,nk)+1  
       enddo
     endif
   enddo
   enddo
!
! accumulate over all sub intervals nt (so result like an average in time)
   do nk=1,nkmax
   do nt=2,ntmax
     ibins(:,1,nk)=ibins(:,1,nk)+ibins(:,nt,nk)
   enddo
   enddo
!
! Write result to file
!
   open (iunit,file=trim(file_out),form='formatted')
   do nk=1,nkmax
     ibsum=0
     nt=1          ! this contains sum over all times
     do n=1,nbins
       ibsum=ibsum+ibins(n,nt,nk)
     enddo
     if (ibsum > 0) then
       rbsum=real(ibsum)
       do i=1,nbins 
         rbins(i)=real(ibins(i,nt,nk))/rbsum
       enddo
       write (iunit,'(2i3,i10)') nk,nt,ibsum
       write (iunit,'(10f10.5)') rbins  
     endif
   enddo
   close (iunit)
   print *,trim(file_out),' written'
!
   end program compute_dx
