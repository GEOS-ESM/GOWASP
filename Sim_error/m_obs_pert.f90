!!
      module m_obs_pert
!
! Module to perturb observations with random perturbations
! For use with observation simulation software for the NASA version of the 
! ECMWF/NCEP/NASA/etc. OSSE system.
!
! Initial code provided by Ronald Errico (Jan. 2008) 
!
      use m_kinds, only : rkind1, rkind2, zero_k2, one_k2
      use m_random_gauss, only : random_gauss_r8
!
      implicit none
!
      private
      public pert_find_itype
      public pert_obs
      public pert_ps
!
      real(rkind2), parameter :: stdv_max=5.0_rkind2  ! max st dev for perts
!
      contains
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_eigen (nk,matrix,evects,evalues)
!
! Call library routine to compute eigenvalues and vectors of a symmetric 
! matrix (stored as an array of 8-byte values).
!
      implicit none
!
      integer,  intent(in)   :: nk
      real(rkind2), intent(in)   :: matrix(nk,nk)
      real(rkind2), intent(out)  :: evects(nk,nk)
      real(rkind2), intent(out)  :: evalues(nk)
!
      integer(4)     :: nwork1  
      integer(4)     :: info 
      integer(4)     :: nk4
      real(rkind2), allocatable :: work1(:) ! must be size >= 3*nk-1
      external dsyev
!
      nk4=nk
      nwork1=3*nk-1
      allocate (work1(nwork1))
      evects=matrix
      call dsyev ( 'V','U',nk4,evects,nk4,evalues,work1,nwork1,info)
      if (info /= 0) then
        print *,'   '
        print *,' * * * * * * * * * * * * * * * * * * * * * * *  '
        print *,' info from dsyev =  ',info
        print *,'   '
      endif
!
      deallocate (work1)  
!
      end subroutine pert_eigen
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_find_itype (n_err3,itype,err_itype,err_index)
!
! Find index of desired subtype in list of subtypes in error tables
!
      implicit none
      integer :: n_err3
      integer :: itype
      integer :: err_itype(n_err3)
      integer :: err_index
      integer :: k  
!     
      err_index=0
      do k=1,n_err3
        if (itype == err_itype(k)) then
          err_index=k
        endif
      enddo
!
      end subroutine pert_find_itype 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_fix_evalues (nk,evalues)
!
!  Set small (relative to largest values) or negative eigenvalues to 
!  a small positive value
!
      integer,  intent(in)    :: nk
      real(rkind2), intent(inout) :: evalues(nk)
!
      integer  :: k
      real(rkind2) :: emax, emin
!
      emax=zero_k2
      do k=1,nk
        if (abs(evalues(k)) > emax) emax=abs(evalues(k))
      enddo
!
      if (rkind2==8) then
        emin=emax*1.e-16
      else
        emin=emax*1.e-8
      endif
!
      do k=1,nk
        if (abs(evalues(k)) < emin) evalues(k)=emin
      enddo
!
      end subroutine pert_fix_evalues
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_interp_stdv (nptab,p_tab,stdv_tab,p,stdv) 
!  
! Interpolate stdv from table defined on p-levels defined by p_tab.
!
      implicit none
!
      integer :: nptab
      real(rkind2) :: p_tab(nptab)
      real(rkind2) :: stdv_tab(nptab)
      real(rkind2) :: p
      real(rkind2) :: stdv
      real(rkind2) :: r1,r2
      integer :: k, klev
!
! first find klev as p-level just below p (i.e., next largers p_lev)

      if (p > p_tab(1) ) then 
        stdv=stdv_tab(1)
      elseif (p < p_tab(nptab) ) then
        stdv=stdv_tab(nptab)
      else
        klev=1
        do k=2,nptab-1
          if (p > p_tab(k)) exit
          klev=k
        enddo 
!
        r2=(p_tab(klev)-p)/(p_tab(klev)-p_tab(klev+1))
        r1=one_k2-r2
        stdv=r1*stdv_tab(klev)+r2*stdv_tab(klev+1)
!       
      endif
!
      end subroutine pert_interp_stdv 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine pert_normalize_evects (nk,evects)
!
!  Normalize each eigenvector such that the sum of its squared 
!  components is 1
!
      implicit none
!
      integer,  intent(in)    :: nk
      real(rkind2), intent(inout) :: evects(nk,nk)
!
      integer  :: k, i
      real(rkind2) :: sum  
     
      do k=1,nk
        sum=zero_k2
        do i=1,nk
          sum=sum+evects(i,k)*evects(i,k)
        enddo
        sum=1./sqrt(sum)
        evects(:,k)=evects(:,k)*sum
      enddo
!
      end subroutine pert_normalize_evects
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
       Subroutine pert_obs (n_err1,n_err2,obs_ndim1,obs_ndim2,nobs,      &
                      obs_nlevs,obs_nfields,icolumns,bmiss,pert_fac,     &
                      err_tab,obs_levels,obs_values,corr_dist,xlat,xlon, &
                      hcorr_id,ltest,lrelhum,d_type,lrad)
!
!  Compute and add random perturbations to observations.
!  Perturbations are drawn from a Gaussian distribution
!  For observations at multiple levels (e.g., for RAOBS), the perts are 
!  drawn independentaly for each eigenvector of the covariance matrix, 
!  assuming a vertical correlation function that has shape that is 
!  approximately Gaussian in vertical separation.  No correlation is assumed 
!  between different fields of the same observation (i.e., between t,q or
!  u,v.
!
!  Specifically, the verical correlation function is c=e(a*z**2) where 
!  z= log(p1/p2) for obs at p-levels p1 and p2, 
!  and a=log(xfac_r)/log(xfac_p)**2;
!  such that when either p1/p2 or p2/p1 = xfac_p, then c=xfac_r.
!  This formulation does not depend on height since for raobs, e.g.,
!  height values for the obs may be unavailable although p values are. 
!  In contrast for GPSRO where the height of the obs are available, 
!  the expression is in terms of height rather than log(p).  Note, however, 
!  that in either case, the user-specified correlation lengths are those that
!  yield correlations of 0.1, not the usual 1/sqrt(e).
!
!  Standard deviations of errors of conventional observations are determined
!  by interpolation from tabular values to the p-level of the observation.
!  For radiances, no interpolation between channels is required to 
!  determine the standard deviations of the corresponding observation 
!  values.  Therefore the standard deviations are simply copied from the
!  table provided.
!
! Note the special treatments of q and some conventional single-level obs
!
      use m_obs_error_table, only : et_l_vc_corr
      use m_random_fields, only: random_fields_get_values
!
      implicit none
!
      integer :: n_err1
      integer :: n_err2
      integer :: obs_ndim1 
      integer :: obs_ndim2 
      integer :: obs_nlevs
      integer :: obs_nfields
      integer :: hcorr_id
      integer :: icolumns(obs_nfields)
      integer :: nobs
      real(rkind2) :: pert_fac
      real(rkind2) :: bmiss
      real(rkind2) :: err_tab(n_err1,n_err2)
      real(rkind2) :: obs_levels(obs_ndim1)
      real(rkind2) :: obs_values(obs_ndim1,obs_ndim2)
      real(rkind2) :: corr_dist(obs_nfields)
      real(rkind2) :: horiz_corr1
      real(rkind2) :: xlon, xlat
      character(len=*) :: d_type
!
      logical :: ltest
      logical :: lrelhum
      logical :: l_conv_obs ! true if conventional obs err. to be horiz. corr.
      logical :: lrad       ! true if this is for radiance obs
!
      integer      :: i, j, nf, icol
      real(rkind2) :: corr_d
      real(rkind2) :: cov(obs_nlevs,obs_nlevs)
      real(rkind2) :: evalues(obs_nlevs)
      real(rkind2) :: evects(obs_nlevs,obs_nlevs)
      real(rkind2) :: obs_old(obs_nlevs)
      real(rkind2) :: pert
      real(rkind2) :: xpert(obs_nlevs)
      real(rkind2) :: sqrtev
      real(rkind2) :: stdv(obs_nlevs)
      real(rkind2) :: x
      real(rkind2) :: xfac
      real(rkind2) :: horiz_uncor
      real(rkind2) :: h_corr(obs_nlevs)
      real(rkind2) :: weights_cu(obs_nlevs,2)  ! weights for cor and uncor part
      real(rkind2), parameter :: xfac_r=0.1_rkind2
!
! Loop over obs fields
      do nf=1,obs_nfields
!
! save original values for later comparison and to save q for bmiss
        obs_old(1:obs_nlevs)=obs_values(1:obs_nlevs,nf)
!    
! If obs. field #2 is for q, perturb as relative humidity
! Therefore convert specific humidity obs to relative humidity obs
        if (lrelhum .and. (nf==2)) then
          call pert_q2rh (obs_nlevs,obs_values(1:obs_nlevs,1),     &
                obs_levels(1:obs_nlevs),obs_values(1:obs_nlevs,2),'Q2RH')
        endif 
!
! Fill arrays of standard deviations for each obs error 
        icol=icolumns(nf)
        do i=1,obs_nlevs
          if (lrad) then   ! copy values for corresponding channels
            stdv(i)=err_tab(i,icol) 
          elseif (trim(d_type)=='GPSRO') then  ! interpolate both lat and lev 
            call gpsro_lat_func (n_err1,n_err2,xlat,err_tab, &
                   obs_levels(i),stdv(i),.false.)
          else  ! interpolate vertically between levels in table
            call pert_interp_stdv (n_err1,err_tab(:,1), &
                   err_tab(:,icol),obs_levels(i),stdv(i))
          endif
        enddo
!
! Determine horizontally correlated parts of errors if requested.
! If channel correlated for radiances, compute random components for 
! horizontally uncorrelated but channel correlated part.
        if (hcorr_id > 0)  then
          call random_fields_get_values (obs_ndim1,obs_nlevs,obs_nfields, &
                               nf,nobs,obs_levels,xlat,xlon,stdv,h_corr,  & 
                               weights_cu,hcorr_id)
          if (lrad .and. et_l_vc_corr) then
            call pert_rad_chcorr (obs_nlevs,hcorr_id,xpert)
          endif
        else
          h_corr(:)=0.       ! default value is no horizontal correlation
          weights_cu(:,1)=0.
          weights_cu(:,2)=1. ! only use uncorrelatd error
        endif
!
! Determine if obs:
! (1) are single-level conventional soundings that never need vertical 
!      correlations;
! (2) have vertical correlation lengths 0; 
! (3) are radiances that never have vertical spatial correlation;
! (4) have vertical correlations handeled in conjunction with horizontal
!     correlations.
        if ( (obs_nlevs <= 1) .or. (corr_dist(nf) <= 1.) .or. lrad .or. &
             (hcorr_id > 0 ) ) then 
!          
          do i=1,obs_nlevs 
            if (obs_values(i,nf) >= bmiss) then  ! original obs is missing value
              obs_values(i,nf)=bmiss             ! retain missing value flag
            else ! add horizontally correlated and uncorrelated parts of perts
              if (hcorr_id > 0 .and. lrad .and. et_l_vc_corr) then
                pert=xpert(i) 
              else
                call random_gauss_r8 (zero_k2,stdv(i),stdv_max,pert) 
              endif
              pert=pert_fac*(weights_cu(i,1)*h_corr(i)+weights_cu(i,2)*pert)
              obs_values(i,nf)=obs_values(i,nf)+pert 
            endif
          enddo        
!
        else  ! multi-level conventional or gpsro obs errors are correlated in 
!               the vertical but not horizontal 
!
! Fill symmetric covariance matrix
          if (trim(d_type) == 'GPSRO') then  
            xfac=log(xfac_r)/corr_dist(nf)**2              
          else
            xfac=log(xfac_r)/(log(corr_dist(nf)))**2 
          endif
!
          do i=1,obs_nlevs
            do j=i,obs_nlevs
              if (trim(d_type) == 'GPSRO') then  
                x=obs_levels(i)-obs_levels(j)          ! difference in height
              else 
                x=log(obs_levels(i)/obs_levels(j))     ! difference in log p
              endif
              cov(i,j)=stdv(i)*stdv(j)*exp(xfac*x*x)
              cov(j,i)=cov(i,j)
            enddo
          enddo
!
! Compute eigen values and normalized eigen vectors of cov matrix
! Also, do not allow relatively small or negative eigenvalues
          call pert_eigen (obs_nlevs,cov,evects,evalues)
          call pert_normalize_evects (obs_nlevs,evects)
          call pert_fix_evalues (obs_nlevs,evalues)
!
! Perturb each eigenvector independently and sum result
          do i=1,obs_nlevs
            sqrtev=sqrt(evalues(i))
            call random_gauss_r8 (zero_k2,sqrtev,stdv_max,pert)
            pert=pert*pert_fac
            do j=1,obs_nlevs
              obs_values(j,nf)=obs_values(j,nf)+pert*evects(j,i)
            enddo
          enddo
!
        endif      ! test on whether vertically correlated    
!
! If obs. field #2 is for q, perturb as relative humidity
! Therefore convert back to specific humdity here
        if (lrelhum .and. (nf==2)) then
          call pert_q2rh (obs_nlevs,obs_values(1:obs_nlevs,1), &
                          obs_levels(1:obs_nlevs),obs_values(1:obs_nlevs,2), & 
                          'RH2Q')
        endif 
!
! Restore missing values in original data
        do i=1,obs_nlevs
          if (obs_old(i) >= bmiss) then
            obs_values(i,nf)=bmiss
          endif
        enddo  
! 
! Print some arrays if testing
        if (ltest) then
          print *,' '
          print ('(a,i6,a,i1,a,f6.2)'),' obs number=',nobs,' nf=',nf,      &
                                       ' vcorr_dist=',corr_dist(nf)
          print ('(2a)'),'  k,      press,       old obs,       new obs,', &
                         '       diff,       stdv'
          do i=1,obs_nlevs
            x=obs_values(i,nf)-obs_old(i)
            print ('(i3,f12.2,1p2e15.4,2e12.2)'),i,obs_levels(i),          &
                                obs_old(i),obs_values(i,nf),x,stdv(i)
          enddo
!          
          if (hcorr_id > 0) then
            print ('(a,2i2)'),' Test hcorr output for id, nf = ',hcorr_id,nf
            print ('(a)'),'  k,    hc_pert,       w_hc,     w_nohc,   pert-fac' 
            do i=1,obs_nlevs 
              print ('(i3,1p4e12.2)'),i,h_corr(i),weights_cu(i,:),pert_fac
            enddo
          endif  ! test on hcorr_id       
!
        endif    ! test on ltest
 
      enddo      ! loop over nfields
!
      end subroutine pert_obs 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
       Subroutine pert_ps (n_err1,n_err2,obs_ndim1,nobs,obs_nlevs, &
                           pert_fac,err_tab,obs_levels,psflag,ltest)
!
!  Compute and add random perturbations to surface pressure 
!
      implicit none
!
      integer :: n_err1
      integer :: n_err2
      integer :: obs_ndim1 
      integer :: obs_nlevs
      integer :: nobs
      logical :: psflag(obs_ndim1)
      real(rkind2) :: pert_fac
      real(rkind2) :: err_tab(n_err1,n_err2)
      real(rkind2) :: obs_levels(obs_ndim1)
      logical :: ltest
!
      integer      :: i
      integer      :: ncount        ! counter for number of surface values 
      real(rkind2) :: obs_old(obs_nlevs)
      real(rkind2) :: pert
      real(rkind2) :: stdv(obs_nlevs)
      real(rkind2) :: x
      integer, parameter :: icol=5  ! column 5 in error table is for ps
!
      ncount=0
      do i=1,obs_nlevs
        if (psflag(i)) then
          ncount=ncount+1
          obs_old(i)=obs_levels(i)
          call pert_interp_stdv (n_err1,err_tab(:,1), &
               err_tab(:,icol),obs_levels(i),stdv(i))
          call random_gauss_r8 (zero_k2,stdv(i),stdv_max,pert)
          obs_levels(i)=obs_levels(i)+pert*pert_fac
        endif
      enddo
! 
! Print some arrays if testing and ps data in observation set
      if (ltest .and. (ncount>0)) then
        print *,' '
        print *,' obs number=',nobs
        print *,' Perturbation of ps observations: old value,', &
                 ' new value, diff, stdv'
        do i=1,obs_nlevs
          if (psflag(i)) then
            x=obs_levels(i)-obs_old(i)
            print ('(i2,1p4e15.4)'),i,obs_old(i),obs_levels(i),x,stdv(i)
          endif 
        enddo
      endif    ! test on ltest
!
      end subroutine pert_ps 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_print_matrix (n1,n2,np,matrix,cname)
!
! Print matrix
!
      implicit none
!
      integer,  intent(in)   :: n1, n2, np
      real(rkind2), intent(in)    :: matrix(n1,n2)
      character(len=*)       :: cname
!
      integer     :: n  
!
      print *,cname,' np,n2=',np,n2
      do n=1,n2
        print ('(1p10e12.2)'),matrix(1:np,n)
      enddo
!
      end subroutine pert_print_matrix
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_print_vector (n1,vector,cname)
!
! Print vector
!
      implicit none
!
      integer,  intent(in)   :: n1
      real(rkind2), intent(in)    :: vector(n1)
      character(len=*)       :: cname
!
      print *,cname,' n1=',n1
      print ('(1p10e12.2)'),vector(:)
!
      end subroutine pert_print_vector
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_q2rh (nlevs,t,p,q,func)
!
! Convert q to rh or rh to q
!
      implicit none
!
      integer,      intent(in)    :: nlevs
      real(rkind2), intent(in)    :: t(nlevs), p(nlevs)
      real(rkind2), intent(inout) :: q(nlevs)
      character(len=*) :: func
!
      integer      :: k, ier
      real(rkind1) :: t2, p2, q2
      real(rkind1), parameter :: ratio4R=0.622 ! mol weight water vap / dry air
      real(rkind1), parameter :: svpt0=273.15, svp1=611.2
      real(rkind1), parameter :: svp2=17.67, svp3=29.65
      real(rkind1) :: a,svp,qs
!    
      do k=1,nlevs  
        ier=0
        t2=t(k)      
        p2=p(k)      
        q2=q(k)
        a=svp2*(t2-svpt0)/(t2-svp3)
        svp=svp1*exp(a)             ! saturation vapor pressure
        qs=svp*ratio4R/(p2-svp)      ! saturation specific humidity
        if (func == 'Q2RH') then
          if (p2 > svp*1.1 .and. t2 > 180.) then 
            q(k)=q(k)/qs
          else
            q(k)=1.e-4
          endif  
        else
          if (p2 > svp*1.1 .and. t2 > 180.) then 
            q(k)=q(k)*qs
          else
            q(k)=1.e-8
          endif  
        endif
      enddo 
!
      end subroutine pert_q2rh
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      subroutine pert_rad_chcorr (obs_nlevs,id,xpert)
!
! Create channel correlated errors that are horizontally uncorrelated.
!
      use m_obs_error_table, only : et_e_vects, et_e_v_sqrt
      use m_kinds, only : zero_k2
!
      implicit none
!     
      integer, intent(in) :: obs_nlevs,id
      real(rkind2), intent(out) :: xpert(obs_nlevs)
!
      integer :: i,j
      real(rkind2) :: stdv_pc
      real(rkind2) :: xp(obs_nlevs)
!
      do i=1,obs_nlevs
        stdv_pc=et_e_v_sqrt(i,id)
        call random_gauss_r8 (zero_k2,stdv_pc,stdv_max,xp(i)) 
      enddo
!
      do i=1,obs_nlevs
        xpert(i)=0.
        do j=1,obs_nlevs
          xpert(i)=xpert(i)+et_e_vects(i,j,id)*xp(j)
        enddo
      enddo
!
      end subroutine pert_rad_chcorr 
!
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      end module m_obs_pert



