      subroutine fit_corr (max_nbins,n_bins,nk,min_count,c_func, &
                           xbin,covs,fits,ierr)
!
! Compute correlation distance and fraction of variance correlated   
! Output the distances and fractions as functions of channel.
!
      use m_kinds, only : rkind1 
!
      implicit none
!
      integer, intent(in) :: max_nbins
      integer, intent(in) :: n_bins
      integer, intent(in) :: nk
      integer, intent(in) :: min_count
      integer, intent(out) :: ierr
      real(rkind1), intent(in) :: xbin 
      real(rkind1), intent(in) :: covs(max_nbins,nk,3)  ! cov, corr, bincount
      real(rkind1), intent(out) :: fits(nk,3)  ! frac_var, length, score
      character(len=*), intent(in) :: c_func

      integer, parameter :: n_afacs=100
      integer, parameter :: n_xlengths=60  ! max length =10*this km
      integer, parameter :: n_funcs=4
      integer :: k1, k2, n, nok, k, nf
      integer :: iscore(2)
!
      real(rkind1), parameter :: xmin_corr=0.05 ! min correl to consider
      real(rkind1) :: afacs(n_afacs), xlengths(n_xlengths)
      real(rkind1) :: sum_weights, score, score1
      real(rkind1) :: xd
      real(rkind1) :: xmin_count
      real(rkind1) :: xc(n_funcs)
      real(rkind1) :: weights(n_bins)
      real(rkind1) :: func_table(n_bins,n_xlengths)
      real(rkind1) :: d(n_bins)
      real(rkind1) :: dmin(nk,n_funcs)
      real(rkind1) :: vmin(nk,n_funcs)
      real(rkind1) :: smin(nk,n_funcs)
!
      character(len=12) :: c_funcs(n_funcs)
!
      c_funcs=(/'GAUSS','EXP','TOAR','WHITE'/)
      xmin_count=real(min_count)
!
      ierr=1
      do nf=1,n_funcs
        if (ierr == 1 .and. trim(c_func) == trim(c_funcs(nf))) then
          ierr=0
        endif
      enddo
      if (ierr /= 0) then
        print ('(a)'),'ERROR: unacceptable function requested in fit_corr:'
        print ('(2a)'),'requested= ',trim(c_func)
        print ('(5(a,1x))'),'acceptable=',c_funcs(:)
        return
      endif 
!
! Set separation distances in km corresponding to mid-range for each bin
      d(1)=0.     
      do n=2,n_bins
        d(n)=(n-1.5)*xbin
      enddo
!
! Set values of afacs and xlengths
      do n=1,n_afacs
        afacs(n)=real(n)/real(n_afacs)
      enddo
      do n=1,n_xlengths
        xlengths(n)=10.*real(n)
      enddo
!
! Create table of functional values
      do n=1,n_bins
        do k1=1,n_xlengths
          xd=d(n)/xlengths(k1)
          if (trim(c_func) == 'GAUSS') then 
            func_table(n,k1)=exp(-0.5*xd*xd) 
          elseif (trim(c_func) == 'EXP') then 
            func_table(n,k1)=exp(-xd) 
          elseif (trim(c_func) == 'TOAR') then 
            func_table(n,k1)=(1+xd+xd*xd/3.)*exp(-xd)   
          elseif (trim(c_func) == 'WHITE') then 
            if (n == 1) then
              func_table(n,k1)=1.
            else
              func_table(n,k1)=0.
            endif
          endif
        enddo
      enddo 
!
! First set weights that define the metric used to measure the 
! separation between the actual correlations and the seelected function
! shaped correlation. The weight here is the same for all distances 
! except that for distance 0. 
      weights(:)=1.
      weights(1)=0.
!
! Computes score that measures the fit between actual correl and a
! selected functioan-shaped correlation for ranges of afac and xlength
! parameters. Save values of parameters that define the best fit.
! Only correlations > 0.05 are considered significant enough to fit.
! The weights used to define each score are normalized here so that the 
! their sum over all bins (separation distances) considered is 1.
      do k=1,nk       
        score=1.e30   
        do k1=1,n_afacs
          do k2=1,n_xlengths
            score1=0.
            sum_weights=0.               
            do n=2,n_bins
!
! only consider bins with non negligible correl and obs-pair counts
              if (covs(n,k,2) >= xmin_corr .and. covs(n,k,3) >= xmin_count ) then 
                sum_weights=sum_weights+weights(n)
                score1=score1+weights(n)* &
                     (afacs(k1)*func_table(n,k2)-covs(n,k,2))**2
              endif
            enddo
            if (sum_weights > 0.) then
              score1=score1/sum_weights
            endif
            if (score1 < score) then
              score=score1
              iscore(1)=k1
              iscore(2)=k2
            endif
          enddo 
        enddo
        fits(k,1)=afacs(iscore(1))
        fits(k,2)=xlengths(iscore(2))
        fits(k,3)=sqrt(score+1.e-20)
      enddo 
!
      end subroutine fit_corr
!
! 
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine construct_corr (max_nbins,n_bins,nk,c_func, &
                                 xbin,corr,fits,ierr)
!
! Compute correlation for prescribed function and parameters
!
      use m_kinds, only : rkind1 
!
      implicit none
!
      integer, intent(in) :: max_nbins
      integer, intent(in) :: n_bins
      integer, intent(in) :: nk
      integer, intent(out) :: ierr
      real(rkind1), intent(in) :: xbin 
      real(rkind1), intent(in) :: fits(nk,3)  ! frac_var, length, score
      real(rkind1), intent(out) :: corr(max_nbins,nk)  ! corr
      character(len=*), intent(in) :: c_func
!
      integer, parameter :: n_funcs=4
      integer :: k,n,nf
      real(rkind1) :: d, xd
      character(len=12) :: c_funcs(n_funcs)
!
      c_funcs=(/'GAUSS','EXP','TOAR','WHITE'/)
!
      ierr=1
      do nf=1,n_funcs
        if (ierr == 1 .and. trim(c_func) == c_funcs(nf)) then
          ierr=0
        endif
      enddo
      if (ierr /= 0) then
        print ('(2a)'),'ERROR: unacceptable function requested in ', &
                      'construct_corr:'
        print ('(2a)'),'requested= ',trim(c_func)
        print ('(5(a,1x))'),'acceptable=',c_funcs(:)
        return
      endif 
!
      do k=1,nk
        corr(1,k)=1.
      enddo
!
      do n=2,n_bins      
        d=(n-1.5)*xbin
        do k=1,nk 
          xd=d/fits(k,2)
          if (trim(c_func) == 'GAUSS') then 
            corr(n,k)=fits(k,1)*exp(-0.5*xd*xd) 
          elseif (trim(c_func) == 'EXP') then 
            corr(n,k)=fits(k,1)*exp(-xd) 
          elseif (trim(c_func) == 'TOAR') then 
            corr(n,k)=fits(k,1)*(1+xd+xd*xd/3.)*exp(-xd)   
          elseif (trim(c_func) == 'WHITE') then 
            corr(n,k)=0.
          endif
        enddo
      enddo 
!
      end subroutine construct_corr
