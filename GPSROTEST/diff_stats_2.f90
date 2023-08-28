   program diff_stats
   implicit none
   integer, parameter :: nbins=90
   integer, parameter :: nrmax=1000000
   integer :: nread,n,ib,nbad,ndatamax,ndataib,nc,nb,nf
   integer :: idata(nrmax), idataBKG(nrmax)
   real(4) :: data(6,nrmax), dataBKG(6,nrmax), datadiff(6)
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
   open (93,file='/discover/nobackup/rerrico/TEST_GPS_517/xxgg45shNR0z.txt',form='unformatted')
   open (94,file='/discover/nobackup/rerrico/TEST_GPS_517/xxgg45shBKG0z.txt',form='unformatted')
   nbad=0
   do n=1,min(nread,nrmax)
     read (93) idata(n),data(:,n)
     read (94) idataBKG(n),dataBKG(:,n)

    if (abs(data(1,n)+50.2)<0.3 .and. abs(data(2,n)-159.4)<0.3) then
      print ('(a,i7,i8,1p6e14.5)'),'NNNR',n,idata(n),data(:,n) 
      print ('(a,i7,i8,1p6e14.5)'),'NNNB',n,idataBKG(n),dataBKG(:,n) 
    endif
!    if (n> 300    .and. n< 700   ) print ('(a,i7,i8,1p6e14.5)'),'NNNR',n,idata(n),data(:,n) 
!    if (n> 300    .and. n< 700   ) print ('(a,i7,i8,1p6e14.5)'),'NNNB',n,idataBKG(n),dataBKG(:,n) 
   
!if (idata(n) /= idataBKG(n) .or. data(1,n) /= dataBKG(1,n) .or. &
!    data(2,n) /= dataBKG(2,n) .or. data(3,n) /= dataBKG(3,n) ) then 
!  nbad=nbad+1
!   if (abs(idataBKG(n) - 4751) <5 .or. abs(idata(n) - 4751) <5) then
!  if (nbad < 100) then
!   print ('(a,i4,3i8,4f10.3,2f10.1)'),'HHH',nbad,n,idata(n),idataBKG(n),data(1,n),dataBKG(1,n), &
!         data(2,n),dataBKG(2,n),data(3,n),dataBKG(3,n)
!  endif

   enddo
!
   ibins(:)=0
   xbins(:,:)=0.
   dz=6.0e4/real(nbins)
   ndataib=int(1+16.*1000./dz)
   nbad=0
   do n=1,min(nread,nrmax)
 
   nf=0
   do  nb=1,min(nread,nrmax)
     if (data(1,n) == dataBKG(1,nb) .and. data(2,n) == dataBKG(2,nb) .and. data(3,n) == dataBKG(3,nb)) then
       nf=nb
     endif
   enddo

    if (nf /= 0) then 
      if (data(5,n) > .999 .or. data(5,n) < 0. .or. &
          data(4,n) > .999 .or. data(4,n) < 0. .or. data(3,n) < 0. .or. &
          dataBKG(5,nf) > .999 .or. dataBKG(5,nf) < 0. .or. &
          dataBKG(4,nf) > .999 .or. dataBKG(4,nf) < 0. .or. dataBKG(3,nf) < 0.) then 
        nf=0
      endif
    endif
    
    if (nf == 0) then    
       nbad=nbad+1
     else
!       datadiff(4)=data(4,n)-data(5,n)    ! ROPP_NR-GSI_NR
!       datadiff(5)=dataBKG(4,n)-data(5,n) ! ROPP_BKG-GSI_NR
!       datadiff(6)=data(4,n)-dataBKG(5,n) ! O-F  ROPP_NR-GSI_BKG
       datadiff(4)=data(4,n)              ! O  ROPP NR
       datadiff(5)=dataBKG(5,nf)           ! F  GSI BKG 
       datadiff(6)=data(4,n)-dataBKG(5,nf) ! O-F  ROPP_NR-GSI_BKG
       ib=1+int(data(3,n)/dz)
       ib=min(ib,nbins)    
       ibins(ib)=ibins(ib)+1
       xbins(ib,1)=xbins(ib,1)+datadiff(4)
       xbins(ib,2)=xbins(ib,2)+datadiff(5)
       xbins(ib,3)=xbins(ib,3)+datadiff(6)
       xbins(ib,4)=xbins(ib,4)+datadiff(4)*datadiff(4)
       xbins(ib,5)=xbins(ib,5)+datadiff(5)*datadiff(5)
       xbins(ib,6)=xbins(ib,6)+datadiff(6)*datadiff(6)

if (n<1000) print ('(a,i6,2f8.2,f8.0,1p4e14.3)'),'XXX',n,data(1:4,n),dataBKG(5,n),datadiff(4),datadiff(6)
!if (ib==14) print ('(a,i6,2f9.2,f9.0,1p4e14.3)'),'XXX',n,data(1:4,n),dataBKG(5,n),datadiff(4),datadiff(6)

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
   print *,'STATS'
   do ib=nbins,1,-1
     z=(ib-1)*dz/1000.
     if (xbins(ib,1) > 0.) then
       xr=xbins(ib,3)/xbins(ib,1)
       xx=xbins(ib,6)/xbins(ib,1)
     else
       xr=0.
       xx=0.
     endif
     print ('(i3,f7.2,i6,1p6e12.3)'),ib,z,ibins(ib),xbins(ib,:)
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


