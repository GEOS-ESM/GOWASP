   program diff_stats
   implicit none
   integer, parameter :: nbins=90
   integer, parameter :: nrmax=1000000
   integer :: nread,n,ib,nbad,ndatamax,ndataib,nc
   real(4) :: data(6,nrmax)
   real(4) :: omf(2,nbins)
   real(4) :: xdata1
   real(8) :: xbins(nbins,6), z, dz, xr, xx
   integer :: ibins(nbins)
   character(len=8) :: cnread
!
   xdata1=0.
   call GetArg( 1_4, cnread)
   read (cnread,'(i8)') nread
   print *,'nread=',nread,nrmax
!
!
   open (10,file='/discover/home/rerrico/GOWASP_3/GPSROTEST/omf.txt')
   do n=1,90 ! nbins
     read (10,*) ib,z,nc,omf(:,n)
   enddo
   close (10)
!
!   open (93,file='/discover/home/rerrico/GOWASP_3/GPSROTEST/xxgg45.txt',form='unformatted')
   open (93,file='/discover/nobackup/rerrico/TEST_GPS6/xxgg1pt.txt',form='unformatted')
   do n=1,min(nread,nrmax)
     read (93) data(:,n)

   if (n> 300    .and. n< 700   ) print ('(a,i7,1p6e15.5)'),'NNN',n,data(:,n) 
   if (n> 152162 .and. n< 152169) print ('(a,i7,1p6e15.5)'),'NNN',n,data(:,n) 
   enddo
!
   ibins(:)=0
   xbins(:,:)=0.
   dz=6.0e4/real(nbins)
   ndataib=int(1+16.*1000./dz)
   nbad=0
   do n=1,min(nread,nrmax)
     if (data(5,n) > .999 .or. data(5,n) < 0. .or. &
         data(4,n) > .999 .or. data(4,n) < 0. .or. data(3,n) < 0.) then 
       nbad=nbad+1
     else
       ib=1+int(data(3,n)/dz)
       ib=min(ib,nbins)    
       ibins(ib)=ibins(ib)+1
       xbins(ib,1)=xbins(ib,1)+data(4,n)
       xbins(ib,2)=xbins(ib,2)+data(5,n)
       xbins(ib,3)=xbins(ib,3)+data(6,n)
       xbins(ib,4)=xbins(ib,4)+data(4,n)*data(4,n)
       xbins(ib,5)=xbins(ib,5)+data(5,n)*data(5,n)
       xbins(ib,6)=xbins(ib,6)+data(6,n)*data(6,n)
       if (ib == ndataib) then
         if (abs(data(6,n)) > xdata1) then
           xdata1=abs(data(6,n))
           ndatamax=n
         endif
       endif
     endif 
   enddo
!
   do ib=1,nbins
     if (ibins(ib) > 0) then
       xbins(ib,:)=xbins(ib,:)/ibins(ib)
     endif
   enddo
!
   do ib=1,nbins
     if (ibins(ib) > 0) then
       do n=4,6
         xbins(ib,n)=xbins(ib,n)-xbins(ib,n-3)*xbins(ib,n-3)
         if (xbins(ib,n) > 0.) xbins(ib,n)=sqrt(xbins(ib,n))
       enddo
     endif
   enddo
!
   print *,' '
   do ib=nbins,1,-1
     z=(ib-1)*dz/1000.
     if (xbins(ib,1) > 0.) then
       xr=xbins(ib,3)/xbins(ib,1)
       xx=xbins(ib,6)/xbins(ib,1)
     else
       xr=0.
       xx=0.
     endif
     print ('(i3,f7.2,i6,2f8.4,1p6e12.3)'),ib,z,ibins(ib),xr,xx,xbins(ib,:)
   enddo
   print *,' '
   print *,'nbad=',nbad
   print *,'ndatamax',ndatamax,ndataib,xdata1
!
!
   print *,' '
   do ib=nbins,1,-1
     z=(ib-1)*dz/1000.
     if (omf(2,ib) > 0.) then
       xx=xbins(ib,6)/omf(2,ib)
     else
       xx=0.
     endif
     print ('(i3,f7.2,f8.4,1p2e12.3)'),ib,z,xx,omf(2,ib),xbins(ib,6)
   enddo


   end program diff_stats  


