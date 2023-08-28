! $Id: ropp_fm_alpha2drk.f90,v 1.1 2013-09-03 21:45:32 rerrico Exp $

!****s* BendingAngle2d/ropp_fm_alpha2drk *
!
! NAME
!    ropp_fm_alpha2drk - Forward model calculating a bending
!                        angle profile from planar information.
!
! SYNOPSIS
!    call ropp_fm_alpha2drk(kobs, klev, ...)
! 
! DESCRIPTION
!    This routine is a forward model calculating the bending angle profile
!    from planar refractivity information.  
!
! INPUTS
!
!           kobs   =  number of observed bending angles 
!           klev   =  number of vertical levels
!           khoriz =  number of horizontal locations
!           ksplit =  splitting of model levels
!           pdsep  =  angular spacing
!           pa     =  impact parameters
!           proc   =  radius of curvature 
!           pz_2d  =  2D impact height (do a 1D calculation above pz_2d)
!           prefrac=  refractivity values
!           pradius=  radius values
!           pnr    =  nr product 
!
! OUTPUT
! 
!           pa_path = impact parameter at end points of ray path
!           palpha  = bending angle values
!           prtan   = radius of tangent point 
!
! NOTES
!    The forward model calculate the bending angle as a function of
!    impact parameter. Below "pz_2d" a Runge-Kutta solver is used
!    calculate the ray-path and bending angle. Above "pz_2d" we use
!    the 1D method, based on the error function solution of the 
!    bending angle integral.  
!
! SEE ALSO
!    ropp_fm_types
!
! AUTHOR
!   ECMWF, UK.
!   Any comments on this software should be given via the GRAS SAF
!   Helpdesk at http://www.grassaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****


SUBROUTINE ropp_fm_alpha2drk(kobs,    & ! no.of observations
                             klev,    & ! no. of vertical levels
		             khoriz,  & ! no. of horizontal layers  ODD
		             ksplit,  &  
		             pdsep,   & ! the angular spacing 
                             pa,      & ! impact parameter values
   		             prefrac, & ! refractivity
		             pradius, & ! radius values
		             pnr,     &
		             proc,    &
		             pz_2d,   &
		             pa_path, &
		             palpha,  &
		             prtan)  

  USE ropp_utils, ONLY: ropp_MDFV
  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_fm_constants, ONLY : pi
  
  IMPLICIT NONE
  
!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)  :: kobs           ! size of ob. vector
  INTEGER,  INTENT(IN)  :: klev           ! no. of refractivity levels
  INTEGER,  INTENT(IN)  :: khoriz         ! no. of horizontal locations
  INTEGER,  INTENT(IN)  :: ksplit         ! 
  REAL(wp), INTENT(IN)  :: pdsep          ! angular spacing of grid
  REAL(wp), INTENT(IN)  :: pa(kobs)       ! impact parameter 
  REAL(wp), INTENT(IN)  :: prefrac(klev,khoriz)  ! refractivity values on levels
  REAL(wp), INTENT(IN)  :: pradius(klev,khoriz)   ! radius values
  REAL(wp), INTENT(IN)  :: pnr(klev,khoriz)
  REAL(wp), INTENT(IN)  :: proc                   ! radius of curvature
  REAL(wp), INTENT(IN)  :: pz_2d
  REAL(wp), INTENT(OUT) :: pa_path(kobs,2)        
  REAL(wp), INTENT(OUT) :: palpha(kobs)   ! bending angle
  REAL(wp), INTENT(OUT) :: prtan(kobs)    ! radius of tangent point
  
!
! local variables
!

  INTEGER :: i,j,in,ibot,jbot,ikbot,iside,ik,ikp1,in_2d
  INTEGER :: ikcen
  REAL(wp), PARAMETER :: zhmax = 3.0E4_wp
  REAL(wp), PARAMETER :: zhmin = 1.0E2_wp
  REAL(wp) :: zrad,zdndr
  REAL(wp) :: zhwt1,zhwt2
  REAL(wp) :: zamult
  REAL(wp) :: zh,zh2,zhuse,zh_up
  REAL(wp) :: zy(4),zyt(4)
  REAL(wp) :: zdydh(4),zdydht(4,4)
  REAL(wp) :: ztheta_tan,ztheta_min,ztheta_max
  REAL(wp) :: zdr_max,zdr_dtheta,zrtan,zdr
  REAL(wp) :: zalpha_half(2)
  REAL(wp) :: zkval(klev-1,khoriz)
  REAL(wp) :: ztlow,ztup,zdalpha,zRoot_halfPI
  REAL(wp) :: zerf_up,zerf_low,zt,zdiff_erf,znr_low,zref_low,zaval
  LOGICAL :: lfirst_1d,LEAVING,lone_d_calc
 
!-------------------------------------------------------------------------------
! 2. Set up the central profile kcen 
!-------------------------------------------------------------------------------

  ikcen = khoriz/2 + 1
  ztheta_tan = REAL(ikcen-1)*pdsep 
  ztheta_min = -ztheta_tan
  ztheta_max =  ztheta_tan
  
!-------------------------------------------------------------------------------
! 3. Set the kvals used in the 1d calculation
!-------------------------------------------------------------------------------

  DO i = 1,klev-1
    DO j = 1, khoriz
      
      zkval(i,j) = LOG(prefrac(i,j)/prefrac(i+1,j)) /  &
                      MAX((pnr(i+1,j) - pnr(i,j)),1.0_wp)
      zkval(i,j) = MAX(1.0E-6_wp,zkval(i,j))
      
    ENDDO
  ENDDO
  
!-------------------------------------------------------------------------------
! 4. Set n_2d level. For levels below n_2d we do 2D ray bending calculation
! above n_2d we do the 1D calculation
!-------------------------------------------------------------------------------

  in_2d = 0 

  DO WHILE ((pnr(in_2d+1,ikcen)-proc < pz_2d) .AND. (in_2d < klev - 1)) 
    
    in_2d = in_2d + 1
    
  ENDDO
      
  jbot = 1
  
  DO
    
    IF (prefrac(jbot,ikcen) > 0.0_wp .AND. pnr(jbot,ikcen) > 0.0_wp) EXIT
    jbot = jbot + 1
    
  ENDDO

  ikbot = klev
  
  DO i=klev,jbot+1,-1

    IF (pnr(ikbot,ikcen) < pnr(ikbot-1,ikcen)) EXIT 
    
    ikbot = ikbot - 1
    
  ENDDO
  
  jbot = MAX(jbot,ikbot)
  
!-------------------------------------------------------------------------------
! 5. Set the outputs to missing
!-------------------------------------------------------------------------------

  palpha(:)=ropp_MDFV
  pa_path(:,:)=ropp_MDFV

!-------------------------------------------------------------------------------
! 6. 2D bending angle calculation
!-------------------------------------------------------------------------------

  zRoot_halfPI = SQRT(0.5_wp*pi)
  
  OBLOOP: DO in=1,kobs
    
    IF (pa(in) < pnr(jbot,ikcen) .OR. pa(in) >= pnr(klev,ikcen)) CYCLE  
    
    
    ibot = jbot
    
    DO 
      
      IF (pnr(ibot+1,ikcen) - pa(in) > 1.0_wp) EXIT 
      
      ibot=ibot+1
      
    ENDDO
    
! 6.1 Calculate the radius at tangent point   
! -----------------------------------------
  
    zrad = 0.5_wp*(pradius(ibot,ikcen)+pradius(ibot+1,ikcen))
    zdndr = 1.0E-6_wp*(prefrac(ibot+1,ikcen)-prefrac(ibot,ikcen))/ &
                 (pradius(ibot+1,ikcen)-pradius(ibot,ikcen)) 
  
    IF ( zrad*zdndr > -1.0_wp) THEN
      
      zrtan = pradius(ibot,ikcen) + &
                (pa(in)-pnr(ibot,ikcen))/(1.0_wp + zrad*zdndr)
      
    ELSE
      
      zrtan = zrad   ! probably in a super-refracting layer
       
    ENDIF
   
! 6.2 If zrtan is within a 1 m of upper level set to upper level
! --------------------------------------------------------------

   IF ((zrtan - pradius(ibot+1,ikcen)) > -1.0_wp) THEN

       ibot = ibot + 1
          
       zrtan = pradius(ibot,ikcen)
       
   ENDIF     

! 6.3 Save the radius of the tangent point. prtan is not active. 
! Just storing for diagnostics. 
! --------------------------------------------------------------

   prtan(in) = zrtan
          	  
! 6.4 Set bending angle value  
! ---------------------------  

   zalpha_half(:) = 0.0_wp
   

   DO iside = 1,2 
      
      pa_path(in,iside) = pa(in)
      lfirst_1d = .TRUE.
      lone_d_calc = .FALSE.
      zamult = 1.0_wp
      IF (iside == 2) zamult  = -1.0_wp
      
! 6.5 Initialise vector
! ---------------------

      zy(1) = 0.0_wp           ! height above tangent point           
      zy(2) = 0.0_wp           ! theta
      zy(3) = ASIN(1.0_wp)     ! thi
      zy(4) = 0.0_wp           ! bending angle 
        
      DO i = ibot,klev-1
      
	
	 IF ( i < MIN(in_2d,klev-1) .AND. khoriz > 1) THEN  
	
            IF (i == ibot) THEN
	 
	       ik = ikcen	 
               zdr_max = (pradius(i+1,ik)- zrtan)/REAL(ksplit) 
               zh = SQRT(2.0_wp*pradius(ibot,ik)*zdr_max)
	    
            ELSE 
 
               ik = INT((zy(2) + ztheta_tan)/pdsep)+1
	       ik = MIN(MAX(1,ik),khoriz)
               zdr_max = (pradius(i+1,ik)- pradius(i,ik))/REAL(ksplit)	       
	       zh_up = SQRT(2.0_wp*pradius(i,ik)*zdr_max)
	       zh = zdr_max/MAX(COS(zy(3)),1.0E-10_wp)	
	       zh = MIN(zh_up,zh) 
	
            ENDIF

! 6.6 Estimate the step-length        
! -----------------------------
	         
! 6.6.1 Limit to horizontal distance between grid points

 	    zh = MAX(MIN(zh,zhmax),zhmin) 	 
	 
	    zh2 = 0.5_wp*zh

! 6.7 Now calculate the path with a RUNGE-KUTTA
! ---------------------------------------------

            LEAVING = .FALSE. 
	 
	    DO  j = 1,ksplit 
	  	  	    	    
   	       zyt(:) = zy(:)  
	    
! 6.7.1 First calculation of derivs

	       CALL ropp_fm_gpspderivs(klev,khoriz,i,pdsep,ztheta_min,  &
                                       ztheta_max,ztheta_tan,zrtan,     &
                                       zamult,prefrac,pradius,zyt,zdydht(:,1)) 
               
	       zyt(:) = zy(:) + zdydht(:,1)*zh2

! 6.7.2 Second call at new yt
	    
	       CALL ropp_fm_gpspderivs(klev,khoriz,i,pdsep,ztheta_min,  &
                                       ztheta_max,ztheta_tan,zrtan,     &
                                       zamult,prefrac,pradius,zyt,zdydht(:,2)) 
               
               zyt(:) = zy(:) + zdydht(:,2)*zh2
	    
! 6.7.3 Third call at new yt
               
	       CALL ropp_fm_gpspderivs(klev,khoriz,i,pdsep,ztheta_min,  &
                                       ztheta_max,ztheta_tan,zrtan,     &
                                       zamult,prefrac,pradius,zyt,zdydht(:,3)) 
               
               zyt(:) = zy(:) + zdydht(:,3)*zh

! 6.7.4 Fourth last call

	       CALL ropp_fm_gpspderivs(klev,khoriz,i,pdsep,ztheta_min,  &
                                       ztheta_max,ztheta_tan,zrtan,     &
                                       zamult,prefrac,pradius,zyt,zdydht(:,4)) 
               
               zdydh(:) = (zdydht(:,1) + zdydht(:,4) + 2.0_wp*(zdydht(:,2) +  &
                             zdydht(:,3)))/6.0_wp
	    
	       zyt(:) = zy(:) + zdydh(:)*zh
	    
! 6.8 Check the radius - have we exited the level
! -----------------------------------------------	    
	    
	       ik = INT((zyt(2) + ztheta_tan)/pdsep)+1
	       ik = MIN(MAX(1,ik),khoriz-1)
	       ikp1 = ik+1
	 
! 6.9 Horizontal weighting factor	      
! ------------------------------- 
	    
	       IF ( zyt(2) < ztheta_max .AND. zyt(2) > ztheta_min) THEN
	    	       
                  zhwt1 = (REAL(ik)*pdsep - (zyt(2)+ztheta_tan))/pdsep    
	          zhwt2 = 1.0_wp - zhwt1
	    
	       ELSE IF (zyt(2) < ztheta_min) THEN	    
	    
	          zhwt1 = 1.0_wp
	          zhwt2 = 0.0_wp
	       
	       ELSE IF (zyt(2) > ztheta_max) THEN
	    
	          zhwt1 = 0.0_wp
	          zhwt2 = 1.0_wp
	       
	       ENDIF    
	          	       	
! 6.10 Radius of pressure level	     
! -----------------------------

  	       zrad = zhwt1*pradius(i+1,ik)+zhwt2*pradius(i+1,ikp1)   

! 6.11 If gone over the boundary scale h
! --------------------------------------	
    
	       IF ( j == ksplit .OR. (zy(1)+zrtan - zrad) > 0.0_wp ) THEN
	    
	          LEAVING = .TRUE.
	    
	          zdr_dtheta = 0.0_wp
	    
	          IF (zyt(2) < ztheta_max .AND. zyt(2) > ztheta_min) &	    
	          zdr_dtheta = (pradius(i+1,ikp1)-pradius(i+1,ik))/pdsep
	
  	          zhuse = zh - (zyt(1)+zrtan-zrad)/(zdydh(1) - zdr_dtheta*zdydh(2))

                  zhuse = MAX(MIN(zhuse,zhmax),zhmin)
	       	  		    
	       ELSE 
	    
	          zhuse = zh  
	       	    
	       ENDIF 
	
! 6.12 Update the position vector
! -------------------------------
	    	 
	       zy(:) = zy(:) + zdydh(:)*zhuse 
	       	    
               IF (LEAVING) EXIT  

! 6.13 Try to maintain roughly the same radial increment by adjusting h
! ---------------------------------------------------------------------

               IF (j < ksplit) THEN
	 
	       zdr = (zrad-zy(1)-zrtan)/REAL(ksplit-j)
	 
	       zh = MIN(zh,zdr/MAX(COS(zy(3)),1.0E-10_wp))
	    
	       zh = MAX(MIN(zh,zhmax),zhmin)
	    
	       zh2 = 0.5_wp*zh
	        		
               ENDIF 

  	    ENDDO  ! complete path thru ith layer

      ELSE 

! ------------------------------------------------------------------------------
! 7. Do 1D calculation
! ------------------------------------------------------------------------------

         IF (lfirst_1d) THEN
	 
            ik = NINT((zy(2) + ztheta_tan)/pdsep)+1
	    ik = MIN(MAX(1,ik),khoriz-1)
	    pa_path(in,iside) = &
                   (1.0_wp+1.0E-6_wp*prefrac(i,ik))*((zy(1)+zrtan)*SIN(zy(3)))
 	    zaval = pa_path(in,iside)
	    zalpha_half(iside) = zy(4)
	    lfirst_1d = .FALSE.
	    
 	 ENDIF     

! 7.1 Continue with 1D bending angle calculation
! ----------------------------------------------

         IF ( i == ibot) THEN 

! we are doing a 1d calc for entire ray path
      
            lone_d_calc = .TRUE. 
            zref_low = prefrac(ibot,ik)*EXP(-zkval(ibot,ik)*(pa(in)-pnr(ibot,ik)))
	    zaval = pa(in)
	    pa_path(in,iside) = pa(in)  
	    znr_low = pa(in)
	
         ELSE 
      
            zref_low = prefrac(i,ik) 
	    znr_low  = pnr(i,ik) 
	 
	 ENDIF

         ztlow = 0.0_wp
	 IF (i > ibot) ztlow = SQRT(MAX(zkval(i,ik)*(pnr(i,ik) - zaval),1.0E-10_wp))

	 ztup = SQRT(MAX(zkval(i,ik)*(pnr(i+1,ik)-zaval),1.0E-10_wp))

! 7.2 Calculate the error functions within this routine rather than an external function call.
! -----------------------------------------------------------------------------

         IF (i == ibot) THEN
	 	 		 
	    zerf_low = 0.0_wp	 
		 
	    zt = 1.0_wp/(1.0_wp+0.47047_wp*ztup)
	    
	    zerf_up= &
           &1.0_wp-(0.3480242_wp-(0.0958798_wp-0.7478556_wp*zt)*zt)*zt*EXP(-(ztup*ztup))	    
	 
	 ELSE IF (i > ibot .AND. i < klev-1) THEN
	 	   
! lower
	    zt = 1.0_wp/(1.0_wp+0.47047_wp*ztlow) 
	    
	    zerf_low = -(0.3480242_wp-(0.0958798_wp-0.7478556_wp*zt)*zt)*zt*EXP(-(ztlow*ztlow)) 

! upper

	    zt = 1.0_wp/(1.0_wp+0.47047_wp*ztup)
	    
	    zerf_up= -(0.3480242_wp-(0.0958798_wp-0.7478556_wp*zt)*zt)*zt*EXP(-(ztup*ztup))	    
	    	    
         ELSE
	 
	    zerf_up = 0.0_wp 
	 
	    zt = 1.0_wp/(1.0_wp+0.47047_wp*ztlow) 
	    
	    zerf_low = -(0.3480242_wp-(0.0958798_wp-0.7478556_wp*zt)*zt)*zt*EXP(-(ztlow*ztlow)) 
	    
	 ENDIF     
          
	  
	 zdiff_erf = zerf_up - zerf_low 

! bending angle	  
	 
         zdalpha    =  &
        & 1.0E-6_wp * zRoot_halfPI* SQRT(zaval*zkval(i,ik)) & 
        & * zref_low*EXP(zkval(i,ik)*(znr_low-zaval))*zdiff_erf 
 
         zalpha_half(iside) = zalpha_half(iside) + zdalpha 

          
	  

      ENDIF       
          
      ENDDO  ! i the layers
      
!
! if we performed a 1d calculation don't evaluate iside = 2
!
     
      IF (lone_d_calc) THEN

         zalpha_half(2) = zalpha_half(1)
	 pa_path(in,2)  = pa_path(in,1)
	  
	 EXIT  ! exiting the iside loop
     
     ENDIF      

   ENDDO ! iside

! ------------------------------------------------------------------------------
! 8. The total bending angle adding both sides
! ------------------------------------------------------------------------------
       
   palpha(in) = zalpha_half(1) + zalpha_half(2) 

ENDDO OBLOOP

RETURN

          
END SUBROUTINE ropp_fm_alpha2drk
