
   subroutine gpsro_lat_func (nerr1,nerr2,lat,errtab,height,stdv,ltest)
!
!  Linearly interpolate in lat and lev from table values
!  For lats beyond edges of the table, use values at edges
!  
   use m_kinds, only : rkind2         ! precision 
!
   integer, intent(in) :: nerr1,nerr2
   logical, intent(in) :: ltest
   real(rkind2), intent(in)  :: errtab(nerr1,nerr2)
   real(rkind2), intent(in)  :: lat
   real(rkind2), intent(in)  :: height ! height above geoid
   real(rkind2), intent(out) :: stdv
!
   integer :: ilats,ilat1,ilat2,ilev1,ilev2
   integer :: i, nlat
   real(rkind2) :: vd,v1,v2,hd,h1,h2,s1,s2
!
   ilats=nint(errtab(1,1)) ! numbers of lats actually specified in table
!
! Determine which lats to interpolate between (lats ordered N to S)
! lats start in second column of table (excluding line index)
   ilat1=2 ! column for northern-most lat encompassing desired lat
   do i=3,ilats+1
     if (lat < errtab(1,i)) then
       ilat1=i
     endif
   enddo
   if (lat >= errtab(1,2) .or. lat < errtab(1,ilats+1)) then
     ilat2=ilat1
   else
     ilat2=ilat1+1  ! column for southern-most lat encompassing desired lat
   endif
!
! Determine which levs to interpolate between (lats ordered lowest to highest)
! levs start in second row of table (first row contains values of specified lats)
   ilev1=2 ! row for table lev below encompassing desired lev
   do i=3,nerr1   
     if (height > errtab(i,1)) then
       ilev1=i
     endif
   enddo
   if (height <= errtab(2,1) .or. height > errtab(nerr1,1)) then
     ilev2=ilev1
   else
     ilev2=ilev1+1 ! row for table lev above encompassing desired lev
   endif
!
! Determine interpolation weights for table levels below (v1) and above (v2)
   if (ilev1 /= ilev2) then
     vd=errtab(ilev2,1)-errtab(ilev1,1)
     v1=(errtab(ilev2,1)-height)/vd
   else  ! above highest or below lowest table levs
     vd=0.
     v1=1.
   endif
   v2=1.-v1
!
! Determine interpolation weights for lats to the north (h1) and south (h2)
   if (ilat1 /= ilat2) then
     hd=errtab(1,ilat1)-errtab(1,ilat2)
     h2=(errtab(1,ilat1)-lat)/hd
   else ! north of northern-most or south of southern-most table lats
     hd=0.
     h2=1.
   endif
   h1=1.-h2
!
! Linearly Interpolate
   s1=h1*errtab(ilev1,ilat1)+h2*errtab(ilev1,ilat2)
   s2=h1*errtab(ilev2,ilat1)+h2*errtab(ilev2,ilat2)
   stdv=v1*s1+v2*s2
!
! Print some values if testing this routine
   if (ltest) then
     print *,' '
     print *,'TEST OF gpsro_lat_func'
     print ('(a,3i5,2f12.2)'),'nerr1,nerr2,ilats,lat,height=', &
            nerr1,nerr2,ilats,lat,height
     print ('(a,4i6)'),'ilev1,ilev2,ilat1,ilat2=',ilev1,ilev2,ilat1,ilat2
     print ('(a,3f10.2,2f7.2)'),'lev1,lev2,vd,v1,v2=', &
            errtab(ilev1,1),errtab(ilev2,1),vd,v1,v2
     print ('(a,3f10.2,2f7.2)'),'lat2,lat1,hd,h,h2=', &
            errtab(1,ilat2),errtab(1,ilat1),hd,h1,h2
     print ('(a,1p4e12.2)'),'tab(ilev1,ilat1:2),tab(ilev2,ilat1:2)=', &
            errtab(ilev1,ilat1:ilat2),errtab(ilev2,ilat1:ilat2)
     print ('(a,1p3e12.3)'),'s1,s2,stdv=',s1,s2,stdv
   endif
!
   end subroutine gpsro_lat_func 
