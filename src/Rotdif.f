!     Rotdif.f: Originally randvec.f +++++++++++++++++++++++++++++++++++++++++++
!     --------------------------------------------------------------------------
      double precision function random(seed)
      implicit none
      real*4 mbig,mseed,mz,fac
      parameter(mbig=4000000.,mseed=1618033.,mz=0,fac=1./mbig)
      real*4 ma(55),mj,mk
      integer iff,i,ii,k,inext,seed,inextp
      save iff,inext,inextp,ma
      data iff /0/
      if(seed.lt.0.or.iff.eq.0) then
        iff=1
        mj=mseed-iabs(seed)
        mj=amod(mj,mbig)
        ma(55)=mj
        mk=1
        do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
        enddo
        do k=1,4
          do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz) ma(i)=ma(i)+mbig
          enddo
        enddo
        inext=0
        inextp=31
        seed=iabs(seed)
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      random=mj*fac
      return
      end function random

! ==============================================================================
!     Rotdif.f: Originally tensorfit.f +++++++++++++++++++++++++++++++++++++++++
!     --------------------------------------------------------------------------
      subroutine tensorfit(random_vectors, nvecs, deff_in, ndeff, lflag, delqfrac,infoflag)
      implicit none

!     Fits the rotational diffusion tensor to ensemble of local diffusion
!     constants.  Uses SVD fit as starting point for minimization of 
!     chi-squared using any of a number of options (simplex, Powell's method)
!     The SVD fit is borrowed from the code rotdif_v5.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Comments from rotdif_v5 (tensor fitting using SVD for small
!     anistropy limit).  

!     use SVD to solve A*Q=Deff for Q
!     SVD subroutines lifted from sec. 2.6 Numerical Recipes in Fortran
!     (Press, Vetterling, Teukolsky, Flannery)
!     m will be the number of local diffusion constants - generated from
!     fit to ln[C(tau)] v tau; n = 6

!     CHANGES IN VERSION 2:
!     (1)  histograms absolute deviations in SVD back calculated Deff 
!          (A*Q=Deff) from MD derived Deff; histograms tau(l=1)/tau(l=2)
!          when l=1 and l=2 data were were used to construct Deff from MD
!     (2)  uses fitted diffusion tensors to compute Deff from Woessner
!          type expressions for <P(l)>, i.e. tau(l)=integral <P(l)> and
!          histograms deviations from MD derived Deff; histograms ratio
!          tau(l=1)/tau(l=2) derived from Woessner like model 
!     (3)  generates A matrix internally from input vectors (subroutine
!          matgen)
!     (4)  no longer back calculates local isotropic time constants for
!          each vector (~ inverse of local isotropic difffusion constant)

!     VERSION 3:
!     (1) computes Diso, anisotropy, and rhombicity
!     (2) changed subroutine woessnertimes to asymtop (identical 
!         routines)


!     VERSION 4:
!     (1) writes out ratios of 2tau(1)/6tau(2) for each vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8, dimension(*), intent(in) :: random_vectors, deff_in
      integer, intent(in) :: nvecs, ndeff, lflag,infoflag
      real*8, intent(in) :: delqfrac

      integer vecidx
      !character(len=80) vecs,deffs,arg
      integer nvec,n,mp,np,maxbins,iarg,idummy,ios,last_arg
      parameter (n=6,mp=2000,np=6,maxbins=200)

      ! Variables for LAPACK dsyev
      integer LWMAX
      parameter (LWMAX = 1000)
      integer info, lwork
      real*8 work( LWMAX )
      external dsyev

      integer svd_chk,itermax,sim_chk,back_cal
      integer rdatflag
      real*8 cut_ratio

      !integer label(mp)
      real*8 deff(mp),sig(mp),ycalc(mp),sumc2(mp)

      real*8 mag,x(mp,3),sgn

      real*8 a(mp,np),copy_a(mp,np),u(mp,np),w(np),v(np,np),q(np)

      real*8 d(3,3),copy1_d(3,3),copy2_d(3,3)
      real*8 dia(3),off(3)

      integer i,j

      real*8 dshape(3)

      integer seed
      real*8 ftol,delqfrac_save
      integer nsimp

      external matgen,svdsolv,svdchk,qtod,similar_trans

      external convertd,simpmin

      common/ chsq/ deff,sig,ycalc,sumc2,nvec

      common/ difftheo/ x,work,lwork
     
!     control variables for calculation:
!     lflag - lflag=1 => only l=1 diffusion constants used to
!             generate input Deff
!             lflag=2 => only l=2 constants are used  (default)
!             lflag=3 => both l=1 and l=2 used
!     cut_ratio - threshold ratio for removing small singular
!                 values in SVD
!     itermax - max number of iteration in tqli diagonalization
!               of diffusion tensor
!     svd_chk - svd_chk=1 => verifies A = U*W*V^T
!     sim_chk - sim_chk=1 => verifies eigenvectors of diffusion
!               tensor by similarity transformation of
!               undiagonalized tensor
!     infoflag - writes out miscillaneous data concerning SVD,
!                undiagonalized diffusion tensor
!     back_cal - writes out all back-calculated local effective
!                diffusion constants from SVD and Woessner type
!                diffusion model
!     deffs  - filename for local diffusion constant vector Deff from
!              MD (which may contain 1/(2*tau(l=1)), 1/(6*tau(l=2)), or both)
!     vecs   - filename for Cartesian coordinates of vectors used to
!              generate Deff from MD (in pdb frame)
!     nsimp  - number of times simplex minimization is done
!     seed   - random number seed for getting initial simplexes
!     delqfrac - how to scale simplexes
!     ftol   - fractional tolerance to end simplex minimization
!             
!     Input defaults:

      !lflag = 2
      cut_ratio = 1.d-6
      itermax = 789
      svd_chk = 0
      sim_chk = 0
      !infoflag = 0
      back_cal = 1
      !deffs = 'deffs'
      !vecs = 'vecs'
      nsimp = 1
      ! NOTE: Setting the random seed to a negative value implicitly 
      !       re-initializes the random number generator. This is
      !       required in order for results to match up with original
      !       tensorfit program.
      seed = -3001796
      !delqfrac = 0.5d0
      ftol = 0.0000001

!     deffs:  local diffusion constant vector Deff from MD (which may
!             contain 1/(2*tau(l=1)), 1/(6*tau(l=2)), or both)

      if (ndeff > mp) then
        write(6,*) 'tensorfit(): too many input deffs'
        return
      endif
      do i=1,ndeff
         deff(i) = deff_in(i)
         sig(i) =1.d0
      end do
!   95 nvec=i-1

      delqfrac_save = delqfrac

!     Read in vectors and generate  A in subroutine matgen.
!     In the event that Deff contains both l=1 and l=2 data, only
!     m/2 vectors exist, but Deff is of length m; the second half
!     of x is a copy of the first.

      vecidx=1 ! index into the random_vectors array
      do i=1,nvecs
        mag=0d0
        do j=1,3
          x(i,j) = random_vectors(vecidx)
          mag = mag + (random_vectors(vecidx) * random_vectors(vecidx))
          vecidx = vecidx + 1
        enddo
        mag = dsqrt(mag);
        do j=1,3
          x(i,j) = x(i,j) / mag
        enddo
      enddo
      if (lflag==3) then
        vecidx=1 ! index into first half of x
        do i=nvecs+1,ndeff
          do j=1,3
            x(i,j) = x(vecidx,j)
          enddo
          vecidx = vecidx + 1
        enddo
      endif
      nvec = nvecs
        
!     generate matrix A(=e*e^T, where e are unit vectors whose 
!     orientational correlation functions were used to extract
!     Deff from MD); to be used to solve A*Q=Deff
      call matgen(nvec,x,a,mp)

!     in svdsolv, a is preserved; a=u*w*v^T
      call svdsolv(a,deff,nvec,n,mp,np,cut_ratio,u,w,v,q)

!     write out results of svd solution
      if(infoflag.gt.0)then
        do i=1,nvec
           write(4,1) (a(i,j),j=1,n)
        end do
1       format(6(e15.8,2x))
        write(4,1)
        do i=1,nvec
           write(4,1) (u(i,j),j=1,n)
        end do
        write(4,1)
        write(4,1) (w(i),i=1,n)
        write(4,1)
        do i=1,n
           write(4,1) (v(i,j),j=1,n)
        end do
      end if

!     array q can be used for the jackknife analysis as an initial guess 
      write(6,*) 'svd q values:'
      write(6,1) (q(i),i=1,n)
      !write(6,'(6(f12.4))') (q(i),i=1,n)
         
!     verify svd a=u*w*v^T
!     WARNING: elements of a are set to 0 inside svdchk prior to final
!     multiplication

      do i=1,nvec
         do j=1,n
            copy_a(i,j)=a(i,j)
         end do
      end do
      if(svd_chk==1)then
        call svdchk(u,w,v,nvec,n,mp,np,copy_a)
        do i=1,nvec
           write(8,1) (copy_a(i,j),j=1,n)
        end do
      end if

!     Tr(Q)=Tr(D) so compute D usin D=3*Diso*I-2*Q; Diso=Tr(D)/3=Tr(Q)/3
!     J. Biomol. NMR 9, 287 (1997) L. K. Lee, M. Rance, W. J. Chazin,
!     A. G. Palmer
!     it has been assumed that the Q(l=1) has the same relationship to
!     D(l=1) as in the l=2 case; we have not yet proved this for the
!     non-symmetric case, but it is true for the axially symmetric case
!     also, Deff(l=1) = 1/(2*tau(l=1)) = e^T * Q(l=1) * e yields the 
!     correct values for tau(l=1) (known from Woessner model type 
!     correlation functions) along principal axes

      call qtod(q,d)
      if(infoflag==0)then
        do i=1,3
           write(9,2) (d(i,j),j=1,3)
        end do
2       format(3(e15.8,2x))
      end if

!     diagonalize D to find principal values and axes
!     tred2/tqli lifted from sec. 11.2 11.3 Numerical Recipes in Fortran
!     d is destroyed by tred2/tqli (columns contain eigenvectors); 
!     store d if desired

      do i=1,3
         do j=1,3
            copy1_d(i,j)=d(i,j)
            copy2_d(i,j)=d(i,j)
         end do
      end do
      lwork=-1
      call dsyev('Vectors','Upper',3,d,3,dia,work,lwork,info)
      lwork = min( LWMAX, int( work( 1 ) ) )
      call dsyev('Vectors','Upper',3,d,3,dia,work,lwork,info)

      call convertd(dia,dshape)

!     principal values of d (eigenvalues) stored in dia
!     principal axes (eigenvectors) stored in columns of d

      write(6,*) 'Results of small anisotropy (SVD) analysis:'
      write(6,'(a,3f10.5)') 'Dav, aniostropy, rhombicity:',(dshape(i),i=1,3)
      write(6,'(a,3f10.5)') 'D tensor eigenvalues       :',(dia(i),i=1,3)
      do i=1,3
         write(6,'(a,i5,a,3f10.5)')'D tensor eigenvector  ',i,':',(d(i,j),j=1,3)
      end do

!     verify eigenvectors by similarity transformation of matrix d

      if(sim_chk==1)then
        call similar_trans(3,3,copy1_d,d)
        do i=1,3
           write(11,2) (copy1_d(i,j),j=1,3)
        end do
      end if

!     back calculate local diffusion constants using SVD results
!     (via A*Q=Deff):

      call locdiff(a,copy2_d,nvec,n,mp,np,ycalc)
      write(6,'(a)') '     taueff(obs) taueff(calc)'
      sgn=0.d0
      do i=1,nvec
         ! for the following chisq fits, convert deff to taueff:
         deff(i) = 1.d0/(6.d0*deff(i))
         ycalc(i) = 1.d0/(6.d0*ycalc(i))
         write(6,'(i5,2f10.5)') i,deff(i),ycalc(i)
         sgn = sgn + ((ycalc(i)-deff(i))/sig(i))**2
      end do
      write(6,'(a,f15.5)')  '  chisq for above is ',sgn
      write(6,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     The chi-squared minimization procedure using
!     diffusion theory (not in the small anisotropy 
!     limit) begins here.

!     Search will be conducted in the space Qxx, Qyy, Qzz, Qxy, etc.
!     The initial guess will be the values found by SVD. 

!#if 0
!!     dac: have initial values just be Diso = Qiso
!
!      q(1) = (3.d0*dshape(1) - dia(1))/2.d0
!      q(2) = (3.d0*dshape(1) - dia(2))/2.d0
!      q(3) = (3.d0*dshape(1) - dia(3))/2.d0
!      q(4) = 0.d0
!      q(5) = 0.d0
!      q(6) = 0.d0
!#endif

      call simpmin(q,n,seed,delqfrac,nsimp,np,ftol,itermax)
      !call gridsrch(q,n,delqfrac_save,np,itermax) !DEBUG
      
      !stop
      end subroutine tensorfit

!     --------------------------------------------------------------------------
      subroutine gridsrch(q,nn,delqfrac,np,tqliitmax)
      implicit none

      integer nn,np,tqliitmax
      real*8 q(np),delqfrac

      integer i,j,k,l,m,n,ii
      real*8 xsmplx(np+1,np),sgn,random

      real*8 xsearch(np),ysearch(np+1),chisq,dia(3),dshape(3),d(3,3)
      integer tqliitmax2
      integer SIZE
      parameter(SIZE=5)

      integer iter
     
      external random,chisq,convertd,amoeba

      integer nvmx,nvec
      parameter(nvmx=2000)
      real*8 ycalc(nvmx),yobs(nvmx),sig(nvmx),sumc2(nvmx),sgn0
      common/ chsq/ yobs,sig,ycalc,sumc2,nvec
      common/ diag/ dia,d,tqliitmax2

      write(6,*)  'grid search:'
      sgn0 = chisq(q)
      write(6,'(a,f15.5)')  '  starting chisq is ',sgn0
!#ifdef TENSORFIT_DEBUG
!      write(6,'(a)') '     taueff(obs) taueff(calc)'
!      do i=1,nvec
!         write(6,'(i5,3f10.5)') i,yobs(i),ycalc(i),sumc2(i)
!      end do
!#endif

!#define SIZE 5

      do i=1,nn
        xsearch(i) = q(i)
      end do

      do i=-SIZE,SIZE
        write(6,*) 'i = ', i
        call flush(6)

        xsearch(1) = q(1) + i*delqfrac/100.

         do j=-SIZE,SIZE
           xsearch(2) = q(2) + j*delqfrac/100.

         do k=-SIZE,SIZE
           xsearch(3) = q(3) + k*delqfrac/100.

         do l=-SIZE,SIZE
           xsearch(4) = q(4) + l*delqfrac/100.

         do m=-SIZE,SIZE
           xsearch(5) = q(5) + m*delqfrac/100.

         do n=-SIZE,SIZE
           xsearch(6) = q(6) + n*delqfrac/100.

         sgn=chisq(xsearch)         
         if( sgn < sgn0 ) write(6,'(f12.5,6f9.5)') sgn, (xsearch(ii),ii=1,6)

!#ifdef TENSORFIT_DEBUG 
!         call convertd(dia,dshape)
!         write(6,*) 'input to amoeba - average at cycle',i
!         write(6,'(a,f15.5)') '   initial chisq = ', sgn
!         write(6,'(a,3f10.5)') 'Dav, aniostropy, rhombicity:',(dshape(k),k=1,3)
!         write(6,'(a,3f10.5)') 'D tensor eigenvalues       :',(dia(k),k=1,3)
!         do k=1,3
!            write(6,'(a,i5,a,3f10.5)')'D tensor eigenvector  ',k,':',(d(k,j),j=1,3)
!         end do
!#endif

        end do
      end do  
        end do
      end do  
        end do
      end do  

      return
      end subroutine gridsrch
        
!     --------------------------------------------------------------------------
      subroutine simpmin(q,n,seed,delqfrac,nsearch,np,ftol,tqliitmax)
      implicit none

!     Driver routine that calls the simplex method optimizer, amoeba.

      integer n,seed,nsearch,np,tqliitmax
      real*8 q(np),delqfrac,ftol

      integer i,j
      real*8 xsmplx(np+1,np),sgn,random

      integer k,l
      real*8 xsearch(np),ysearch(np+1),chisq,dia(3),dshape(3),d(3,3)
      integer tqliitmax2

      integer iter
     
      external random,chisq,convertd,amoeba

      integer nvmx,nvec
      parameter(nvmx=2000)
      real*8 ycalc(nvmx),yobs(nvmx),sig(nvmx),sumc2(nvmx)
      common/ chsq/ yobs,sig,ycalc,sumc2,nvec
      common/ diag/ dia,d,tqliitmax2

      tqliitmax2=tqliitmax

!     first, back-calculate with the SVD tensor, but with the full anisotropy:

      write(6,*)  'same diffusion tensor, but full anisotropy:'
      sgn = chisq(q)
      write(6,'(a,f15.5)')  '  chisq for SVD tensor is ',sgn
      write(6,'(a)') '     taueff(obs) taueff(calc)'
      do i=1,nvec
         write(6,'(i5,3f10.5)') i,yobs(i),ycalc(i),sumc2(i)
      end do
      !return ! DEBUG
!     In the simplex method, N+1 initial points (where N is the 
!     dimension of the search space) must be chosen; the SVD
!     solution provides one of these. 
!     Initial points are stored in rows of xsmplx.
!     Components of vector delq should be of the order of the 
!     characteristic "lengthscales" over which the Q 
!     tensor varies.  delqfrac determines the size of variation
!     for each of the components of Q; the sign of the variation
!     is randomly chosen.

!     Now execute the simplex search method with initial vertices,
!     xsimplx, and chi-squared values, ysearch.
!     We restart the minimization a pre-set number of times to avoid
!     an anomalous result.

      do i=1,n
        xsearch(i) = q(i)
      end do

      do i=1,nsearch

         do j=1,n
            xsmplx(1,j)=xsearch(j)
         end do
         do j=1,n
            do k=1,n
               if(j==k)then
                 sgn=dsign(1d0,random(seed)-0.5d0)
                 xsmplx(j+1,k)=xsmplx(1,k)*(1d0+sgn*delqfrac)
               else
                 xsmplx(j+1,k)=xsmplx(1,k)
               end if
            end do
         end do

         ! As  to amoeba, chi-squared must be evaluated for all
         ! vertices in the initial simplex.

         do j=1,n+1
            do k=1,n
               xsearch(k)=xsmplx(j,k)
            end do
            ysearch(j)=chisq(xsearch)
         end do

         !  Average the vertices and compute details of the average.

         do j=1,n
            xsearch(j)=0d0
         end do
         do j=1,n
            do k=1,n+1 
               xsearch(j)=xsearch(j)+xsmplx(k,j)
            end do
            xsearch(j)=xsearch(j)/dble(n+1)
         end do
         sgn=chisq(xsearch)         
         call convertd(dia,dshape)

         write(6,*) 'input to amoeba - average at cycle',i
         write(6,'(a,f15.5)') '   initial chisq = ', sgn
         write(6,'(a,3f10.5)') 'Dav, aniostropy, rhombicity:',(dshape(k),k=1,3)
         write(6,'(a,3f10.5)') 'D tensor eigenvalues       :',(dia(k),k=1,3)
         do k=1,3
            write(6,'(a,i5,a,3f10.5)')'D tensor eigenvector  ',k,':',(d(k,j),j=1,3)
         end do

         call amoeba(xsmplx,ysearch,np+1,np,n,ftol,chisq,iter)

         do j=1,n+1
            do k=1,n
               xsearch(k)=xsmplx(j,k)
            end do
            ysearch(j)=chisq(xsearch)
         end do

         ! Average the vertices and compute details of the average.

         do j=1,n
            xsearch(j)=0d0
         end do
         do j=1,n
            do k=1,n+1
               xsearch(j)=xsearch(j)+xsmplx(k,j)
            end do
            xsearch(j)=xsearch(j)/dble(n+1)
         end do
         sgn=chisq(xsearch)          
         call convertd(dia,dshape)

         write(6,*) 'output from amoeba - average at cycle',i
         write(6,'(a,f15.5)') '   final   chisq = ', sgn
         write(6,'(a,3f10.5)') 'Dav, aniostropy, rhombicity:',(dshape(k),k=1,3)
         write(6,'(a,3f10.5)') 'D tensor eigenvalues       :',(dia(k),k=1,3)
         do k=1,3
            write(6,'(a,i5,a,3f10.5)')'D tensor eigenvector  ',k,':',(d(k,j),j=1,3)
         end do
         write(6,'(a)') '     taueff(obs) taueff(calc)'
         do k=1,nvec
            write(6,'(i5,3f10.5)')  k,yobs(k),ycalc(k),sumc2(k)
         end do
           
         !  cycle over main loop, but first reduce the size of delqfrac:

         delqfrac = 0.750 * delqfrac
         write(6,'(a,f15.7)') 'setting delqfrac to ',delqfrac

      end do  

      !  set q vector to the final average result from simpmin:

      q(1:n) = xsearch(1:n)

      return
      end subroutine simpmin

!     --------------------------------------------------------------------------
      subroutine convertd(dxyz,dshape)
      implicit none

!     compute Dav, anisotropy, rhombicity
!     definitions found in J. Biomol. NMR, 27, 261 (2003), Hall J B and
!     Fushman D

      real*8 dxyz(3),dshape(3)     

!     Diso = (Dx + Dy + Dz)/3
      dshape(1)=(dxyz(1)+dxyz(2)+dxyz(3))/3d0

!     anistropy = (2*Dz)/(Dx + Dy)
!     rhombicity = 1.5*(Dy - Dx)/[Dz - 0.5*(Dx + Dy)]

      if((dxyz(1)>dxyz(2)).and.(dxyz(1)>dxyz(3)))then
        dshape(2)=(2d0*dxyz(1))/(dxyz(2)+dxyz(3))
        dshape(3)=(1.5d0*(dabs(dxyz(2)-dxyz(3))))/  &
                  (dxyz(1)-0.5d0*(dxyz(2)+dxyz(3)))
      elseif((dxyz(2)>dxyz(1)).and.(dxyz(2)>dxyz(3)))then
        dshape(2)=(2d0*dxyz(2))/(dxyz(1)+dxyz(3))
        dshape(3)=(1.5d0*(dabs(dxyz(1)-dxyz(3))))/  &
                  (dxyz(2)-0.5d0*(dxyz(1)+dxyz(3)))
      elseif((dxyz(3)>dxyz(2)).and.(dxyz(3)>dxyz(1)))then
        dshape(2)=(2d0*dxyz(3))/(dxyz(1)+dxyz(2))
        dshape(3)=(1.5d0*(dabs(dxyz(1)-dxyz(2))))/  &
                  (dxyz(3)-0.5d0*(dxyz(1)+dxyz(2)))
      end if

      return
      end subroutine convertd

!     --------------------------------------------------------------------------
      real*8 function chisq(x)
      implicit none

      real*8 x(6)

      real*8 d(3,3),dia(3),off(3)
      integer itermax
      
      integer nvmx,nvec
      parameter(nvmx=2000)
      real*8 vec(nvmx,3),tau1(nvmx),tau2(nvmx),sumc2(nvmx)

      ! Parameters for dsyev
      integer LWMAX
      parameter(LWMAX=1000)
      real*8 work(LWMAX)
      integer info, lwork
      external dsyev

      integer i
      real*8 ycalc(nvmx),yobs(nvmx),sig(nvmx)
      common/ chsq/ yobs,sig,ycalc,sumc2,nvec

      common/ diag/ dia,d,itermax
      common/ difftheo/ vec,work,lwork

      call qtod(x,d)
      call dsyev('Vectors','Upper',3,d,3,dia,work,lwork,info)
      call asymtop(dia,d,vec,nvec,nvmx,tau1,tau2,sumc2)
      do i=1,nvec
!     for the "tau-eff fit", don't convert to deff
!        ycalc(i)=1.d0/(6.d0*tau2(i))
         ycalc(i) = tau2(i)
      end do
      chisq=0d0
      do i=1,nvec
         chisq=chisq+((yobs(i)-ycalc(i))/sig(i))**2
      end do

      return
      end function chisq

!     --------------------------------------------------------------------------
      subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      implicit none

      integer iter,mp,ndim,np,nmax,itmax
      real*8 ftol,p(mp,np),y(mp),funk
!     parameter (itmax=100000)
      parameter (itmax=10000)
!     parameter (nmax=20,itmax=5000)
      external funk
!     uses amotry,funk

      integer i,ihi,ilo,inhi,j,m,n
      real*8 rtol,sum,swap,ysave,ytry,psum(np),amotry
!     real*8 rtol,sum,swap,ysave,ytry,psum(nmax),amotry

      iter=0
!1     write (6,'(a,i6)') 'hit loop one ',iter 
1     do n=1,ndim
         sum=0d0
         do m=1,ndim+1
            !write (6,'(a,i6,i6,f10.5)') 'xsmplx ',m,n,p(m,n)
            !write (6,'(a,i6,f10.5)') 'ysearch ',m,y(m)
            sum=sum+p(m,n)
         end do
         psum(n)=sum
      end do
!2     write (6,'(a,i6)') 'hit loop two ',iter
      !do n=1,ndim
      !   write (6,'(8x,a,i6,f10.5)') 'psum ',n,psum(n)
      !enddo 
2     ilo=1
      if(y(1)>y(2))then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      end if
      do i=1,ndim+1
         if(y(i)<=y(ilo)) ilo=i
         if(y(i)>y(ihi))then
           inhi=ihi
           ihi=i
         else if(y(i)>y(inhi))then
           if(i/=ihi) inhi=i
         end if
      end do
      !write(6,'(a,2(f10.5))') 'yihi yilo = ',y(ihi),y(ilo)
      !write(6,'(8x,a,2(i6))') 'yihi yilo = ',ihi,ilo
      rtol=2d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo)))
      if(rtol<ftol)then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do n=1,ndim
           swap=p(1,n)
           p(1,n)=p(ilo,n)
           p(ilo,n)=swap
        end do
        return
      end if
      write(6,'(8x,a,i7,e15.6)') 'in amoeba, iter, rtol = ', iter, rtol
      if(iter>=itmax)then
        write(6,*) 'itmax exceeded in amoeba:', itmax
        stop
      end if
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1d0)
      !write (6,'(8x,a,i6,f10.5)') 'ytry ',iter,ytry
      if(ytry<=y(ilo))then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2d0)
        !write(6,'(8x,a,f10.5)') 'case 1 ',ytry
      else if(ytry>=y(inhi))then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)
        !write(6,'(8x,a,f10.5)') 'case 2 ',ytry
        if(ytry>=ysave)then
          do i=1,ndim+1
             if(i/=ilo)then
               do j=1,ndim
                  psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                  p(i,j)=psum(j)
               end do
               y(i)=funk(psum)
             end if
          end do
          iter=iter+ndim
          go to 1
        end if
      else
        !write(6,'(8x,a)') 'case 3 '
        iter=iter-1
      end if
      go to 2
      
      end subroutine amoeba

!     --------------------------------------------------------------------------
      real*8 function amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      implicit none

      integer ihi,mp,ndim,np,nmax
      real*8 fac,p(mp,np),psum(np),y(mp),funk
!     parameter (nmax=20)
      external funk
!     uses funk

      integer j
      real*8 fac1,fac2,ytry,ptry(np)
!     real*8 fac1,fac2,ytry,ptry(nmax)
      
      fac1=(1d0-fac)/ndim
      fac2=fac1-fac
      do j=1,ndim
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
         !write(6,'(16x,a,i6,f10.5)') 'amotry: ',j,ptry(j)
      end do
      ytry=funk(ptry)
      if(ytry<y(ihi))then
        y(ihi)=ytry
        do j=1,ndim
           psum(j)=psum(j)-p(ihi,j)+ptry(j)
           p(ihi,j)=ptry(j)
           !write(6,'(16x,a,i6,i6,f10.5)') 'amotryx: ',ihi,j,p(ihi,j)
        end do
      end if
      amotry=ytry
      
      return
      end function amotry
           
!     --------------------------------------------------------------------------
      subroutine svdsolv(a,b,m,n,mp,np,cut_ratio,u,w,v,x)
      implicit none

      integer m,n,mp,np
      real*8 a(mp,np),w(np),v(np,np),b(mp),x(np)

      real*8 u(mp,np)
      
      integer i,j
      real*8 wmax,wmin,cut_ratio

!     solve for the tensor Q as an intermediate step to computing 
!     the rotational diffusion tensor; see J. Mag. Res. 149, 204 (2001)
!     R. Ghose, D. Fushman, D. Cowburn and related material in
!     J. Biomol. NMR 9, 287 (1997) L. K. Lee, M. Rance, W. J. Chazin,
!     A. G. Palmer
!     Science 268, 886 (1995) R. Bruschweiler, X. Liao, P. E. Wright

!     uses svbksb and svdcmp to implement the svd/backsubstituion method
!     for solving linear system as detailed in sec. 2.6 of
!     Numerical Recipes in Fortran 
!     uses driver routine from Numerical Recipes to set the 
!     threshold for keeping singular values

!     problem cast into form A*Q=Deff (or a*x=b)

!     since a is destroyed by svdcmp, save here
      do i=1,m
         !write(6,'(a,i6,a,6(f12.6))') 'A[',i,']: ',a(i,1),a(i,2),a(i,3),a(i,4),a(i,5),a(i,6)
         do j=1,n
            !write(6,'(2(a,i6),a,f12.4)') 'A[',i,',',j,']=',a(i,j)
            u(i,j)=a(i,j)
         end do
      end do
      call svdcmp(u,m,n,mp,np,w,v)
      !do i=1,n
      !  write(6,'(a,i6,f12.6)') 'svdcmp_w1 ', i, w(i)
      !enddo
      wmax=0d0
      do i=1,n
         if(w(i)>wmax) wmax=w(i)
      end do
      wmin=wmax*cut_ratio
      do i=1,n
         if(w(i)<wmin) w(i)=0d0
      end do
      do i=1,n
        write(6,'(a,i6,f12.6)') 'svdcmp_w2 ', i, w(i)
      enddo
      !do i=1,m
      !  do j=1,m
      !    write(6,'(2(a,i6),a,f12.4)') 'U[',i,',',j,']=',u(i,j)
      !  enddo
      !enddo
      !do i=1,n
      !  do j=1,n
      !    write(6,'(2(a,i6),a,f12.4)') 'Vt[',i,',',j,']=',v(j,i)
      !  enddo
      !enddo
      call svbksb(u,w,v,m,n,mp,np,b,x)

      return
      end

!     --------------------------------------------------------------------------
      subroutine svdcmp(a,m,n,mp,np,w,v)

!     subroutine svdcmp: source code from sec. 2.6 on singular value
!                        decomposition from Numerical Recipes in Fortran, 
!                        2nd edition, Press/Vetterling/Flannery/Teukolsky
!     given matrix a(m,n) (physical dimensions mp,np), computes svd
!     a=u*w*v^T; u replaces a on output; diagonal matrix of singular 
!     values w is output as vector w(n); matrix v (not v^T) is output
!     as v(n,n)
!     calls pythag

      implicit none
      integer m,mp,n,np,nmax
      real*8 a(mp,np),v(np,np),w(np)
      parameter (nmax=500)

      integer i,its,j,jj,k,l,nm
      real*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(nmax),pythag

      g=0d0
      scale=0d0
      anorm=0d0
      do i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0d0
         s=0d0
         scale=0d0
         if(i<=m)then
           do k=i,m
              scale=scale+dabs(a(k,i))
           end do
           if(scale/=0d0)then
             do k=i,m
                a(k,i)=a(k,i)/scale
                s=s+a(k,i)*a(k,i)
             end do
             f=a(i,i)
             g=-dsign(dsqrt(s),f)
             h=f*g-s
             a(i,i)=f-g
             do j=l,n
                s=0d0
                do k=i,m
                   s=s+a(k,i)*a(k,j)
                end do
                f=s/h
                do k=i,m
                   a(k,j)=a(k,j)+f*a(k,i)
                end do
             end do
             do k=i,m
                a(k,i)=scale*a(k,i)
             end do
           end if
         end if
         w(i)=scale*g
         g=0d0
         s=0d0
         scale=0d0
         if((i<=m).and.(i/=n))then
           do k=l,n
              scale=scale+dabs(a(i,k))
           end do
           if(scale/=0d0)then
             do k=l,n
                a(i,k)=a(i,k)/scale
                s=s+a(i,k)*a(i,k)
             end do
             f=a(i,l)
             g=-dsign(dsqrt(s),f)
             h=f*g-s
             a(i,l)=f-g
             do k=l,n
                rv1(k)=a(i,k)/h
             end do
             do j=l,m
                s=0d0
                do k=l,n
                   s=s+a(j,k)*a(i,k)
                end do
                do k=l,n
                   a(j,k)=a(j,k)+s*rv1(k)
                end do
             end do
             do k=l,n
                a(i,k)=scale*a(i,k)
             end do
           end if
         end if
         anorm=dmax1(anorm,(dabs(w(i))+dabs(rv1(i))))
      end do
      do i=n,1,-1
         if(i<n)then
           if(g/=0d0)then
             do j=l,n
                v(j,i)=(a(i,j)/a(i,l))/g
             end do
             do j=l,n
                s=0d0
                do k=l,n
                   s=s+a(i,k)*v(k,j)
                end do
                do k=l,n
                   v(k,j)=v(k,j)+s*v(k,i)
                end do
             end do
           end if
           do j=l,n
              v(i,j)=0d0
              v(j,i)=0d0
           end do
         end if
         v(i,i)=1d0
         g=rv1(i)
         l=i
      end do
      do i=min0(m,n),1,-1
         l=i+1
         g=w(i)
         do j=l,n
            a(i,j)=0d0
         end do
         if(g/=0d0)then
           g=1d0/g
           do j=l,n
              s=0d0
              do k=l,m
                 s=s+a(k,i)*a(k,j)
              end do
              f=(s/a(i,i))*g
              do k=i,m
                 a(k,j)=a(k,j)+f*a(k,i)
              end do
           end do
           do j=i,m
              a(j,i)=a(j,i)*g
           end do
         else
           do j=i,m
              a(j,i)=0d0
           end do
         end if
         a(i,i)=a(i,i)+1d0
      end do
      do k=n,1,-1
         do its=1,30
            do l=k,1,-1
               nm=l-1
               if((dabs(rv1(l))+anorm)==anorm) goto 2
               if((dabs(w(nm))+anorm)==anorm) goto 1
            end do
1           c=0d0
            s=1d0
            do i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if((dabs(f)+anorm)==anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=1d0/h
               c=(g*h)
               s=-(f*h)
               do j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
               end do
            end do
2           z=w(k)
            if(l==k)then
              if(z<0d0)then
                w(k)=-z
                do j=1,n
                   v(j,k)=-v(j,k)
                end do        
              end if
              go to 3
            end if
            if(its==30)then
              write(6,*) 'no convergence in svdcmp'
              stop
            end if
!           if(its==30) pause 'no convergence in svdcmp'
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2d0*h*y)
            g=pythag(f,1d0)
            f=((x-z)*(x+z)+h*((y/(f+dsign(g,f)))-h))/x
            c=1d0
            s=1d0
            do j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f=(x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)=(x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
               end do
               z=pythag(f,h)
               w(j)=z
               if(z/=0d0)then
                 z=1d0/z
                 c=f*z
                 s=h*z
               end if
               f=(c*g)+(s*y)
               x=-(s*g)+(c*y)
               do jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)=(y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
               end do
            end do
            rv1(l)=0d0
            rv1(k)=f
            w(k)=x
         end do
3        continue
      end do
     
      return
      end subroutine svdcmp

!     --------------------------------------------------------------------------
      real*8 function pythag(a,b)
      implicit none

      real*8 a,b

      real*8 absa,absb

      absa=dabs(a)
      absb=dabs(b)
      if(absa>absb)then
        pythag=absa*dsqrt(1d0+(absb/absa)**2)
      else
        if(absb==0d0)then
          pythag=0d0
        else
          pythag=absb*dsqrt(1d0+(absa/absb)**2)
        end if
      end if

      return
      end 

!     --------------------------------------------------------------------------
      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
      
!     subroutine svdcmp: source code from sec. 2.6 on singular value
!                        decomposition from Numerical Recipes in Fortran,
!                        2nd edition, Press/Vetterling/Flannery/Teukolsky
!     solves a*x=b for vector x; a is specified by arrays u,w,v returned
!     from svdcmp; m,n are dimensions of a (mp,np physical dimensions); 
!     b is the RHS of a*x=b; x is the output solution vector; no inputs
!     are destroyed

      implicit none
      
      integer m,mp,n,np,nmax
      real*8 b(mp),u(mp,np),v(np,np),w(np),x(np)
      parameter (nmax=500)

      integer i,j,jj
      real*8 s,tmp(nmax)

      do j=1,n
         s=0d0
         if(w(j)/=0d0)then
           do i=1,m
              s=s+u(i,j)*b(i)
           end do
           s=s/w(j)
         end if
         tmp(j)=s
      end do
      do j=1,n
         s=0d0
         do jj=1,n
            s=s+v(j,jj)*tmp(jj)
         end do
         x(j)=s
      end do

      return
      end subroutine svbksb

!     --------------------------------------------------------------------------
      subroutine svdchk(u,w,v,m,n,mp,np,a)
      implicit none

      integer m,n,mp,np
      real*8 u(mp,np),w(np),v(np,np)

      integer i,j,k
      real*8 scratch(np,np),a(mp,np)

!     verify svd a=u*w*v^T in two steps

!     compute w*v^T; since w is diagonal, only one term per element
!     of (w*v^T); row of (w*v^T) determined by w and column determined
!     by column (row) of v^T (v)
      do i=1,n
         do j=1,n
            scratch(j,i)=w(j)*v(i,j)
         end do
      end do

!     compute u*(w*v^T)
      do i=1,m
         do j=1,n
            a(i,j)=0d0
            do k=1,n
               a(i,j)=a(i,j)+u(i,k)*scratch(k,j)
            end do
         end do
      end do

      return
      end

!     --------------------------------------------------------------------------
      subroutine qtod(q,d)
      implicit none

!     given tensor Q from svd solution to A*Q=Deff, compute rotational 
!     diffusion tensor D from D=3*Diso*I-2*Q, Diso=Tr(D)/3=Tr(Q)/3
!     J. Biomol. NMR 9, 287 (1997) L. K. Lee, M. Rance, W. J. Chazin,
!     A. G. Palmer
!     J. Mag. Res. 149, 204 (2001) R. Ghose, D. Fushman, D. Cowburn

      real*8 q(6),d(3,3)
 
      real*8 trq
      integer i,j

!     recall that here Q is written as a column vector ordered as
!     (Qxx Qyy Qzz Qxy Qxz Qyz)^T

!     compute Tr(Q)
      trq=0d0
      do i=1,3
         trq=trq+q(i)
      end do
   
!     tred2/tqli will be used to diagonalize D, so it must be stored
!     in conventional format (i.e. not as a column vector)
      d(1,1) = trq - (2d0 * q(1))
      d(1,2) = -2d0 * q(4)
      d(1,3) = -2d0 * q(6)
      d(2,1) = d(1,2)
      d(2,2) = trq - (2d0 * q(2))
      d(2,3) = -2d0 * q(5)
      d(3,1) = d(1,3)
      d(3,2) = d(2,3)
      d(3,3) = trq - (2d0 * q(3))
      !do i=1,3
      !   d(i,i)=trq-2d0*q(i)
      !   do j=i+1,3
      !      d(i,j)=-2d0*q(i+j+1)
      !      d(j,i)=d(i,j)
      !   end do
      !end do

      return
      end

!     --------------------------------------------------------------------------
      subroutine similar_trans(ndiag,nmax,mat,z)
      implicit none

      integer i,j,k,ndiag,nmax
      real*8 scratch_mat(nmax,nmax),mat(nmax,nmax),z(nmax,nmax)

      do i=1,ndiag
         do j=1,ndiag
            scratch_mat(i,j)=0d0
         end do
      end do

!     evaluate product mat*z
!     the k loop multiplies the ith row of mat by the jth vector
      do i=1,ndiag
         do j=1,ndiag
            do k=1,ndiag
               scratch_mat(i,j)=scratch_mat(i,j)+mat(i,k)*z(k,j)
            end do
         end do
      end do

      do i=1,ndiag
         do j=1,ndiag
            mat(i,j)=0d0
         end do
      end do

!     for an orthonormal basis, the inverse is the transpose
!     evaluate product z^T*(mat*z)
!     the k loop multiplies the ith row of z^T (ith column of z) by the
!     jth column of (mat*z)
      do i=1,ndiag
         do j=1,ndiag
            do k=1,ndiag
               mat(i,j)=mat(i,j)+z(k,i)*scratch_mat(k,j)
            end do
         end do
      end do

      return
      end subroutine similar_trans

!     --------------------------------------------------------------------------
      subroutine locdiff(a,d,m,n,mp,np,deff)
      implicit none

!     given diffusion tensor D, computes Q from
!     Q=(3*Diso*I-D)/2, then computes local diffusion constants Deff
!     using A*Q=Deff
!     see [21] and [29] in J. Mag. Res. 149, 204 (2001) R. Ghose,
!     D. Fushman, A. G. Palmer

      integer m,n,mp,np
      real*8 a(mp,np),d(3,3),deff(mp)

      real*8 q(np)

      integer i,j

!     compute Q
      call dtoq(d,q)
      do i=1,6
        write(6,'(a,i6,f10.5)') 'dtoq ',i,q(i)
      enddo

!     multiply A*Q (=Deff)
      do i=1,m
         deff(i)=0d0
         do j=1,n
            deff(i)=deff(i)+a(i,j)*q(j)
         end do
      end do

      return
      end

!     --------------------------------------------------------------------------
      subroutine dtoq(d,q)
      implicit none

!     given rotational diffusion tensor D, compute Q from
!     Q=(3*Diso*I-D)/2; see [21] in J. Mag. Res. 149, 204 (2001)
!     R. Ghose, D. Fushman, A. G. Palmer

      real*8 q(6),d(3,3)

      real*8 trd
      integer i,j

!     Diso=Tr(D)/3; compute trace(D)
      trd=0d0
      do i=1,3
         trd=trd+d(i,i)
      end do

!     compute Q
      q(1) = 0.5d0 * (trd - d(1,1))
      q(2) = 0.5d0 * (trd - d(2,2))
      q(3) = 0.5d0 * (trd - d(3,3))
      q(4) = -0.5d0 * d(1,2)
      q(5) = -0.5d0 * d(2,3)
      q(6) = -0.5d0 * d(1,3)
      !do i=1,3
      !   q(i)=0.5d0*(trd-d(i,i))
      !   do j=i+1,3
      !      q(i+j+1)=-0.5d0*d(i,j)
      !   end do
      !end do

      return
      end

!     --------------------------------------------------------------------------
      subroutine matgen(nvec,x,mat,nmax)
      implicit none

      integer nvec,i,j,k,nmax
      real*8 x(nmax,3),mat(nmax,6)

      real*8 mag

!     given x,y,z coordinates of a vector, generate elements of 
!     matrix A (see [29] in J. Mag. Res., 149, 204, (2001) R. Ghose,
!     D. Fushman, D. Cowburn) from N-H or global vectors used to estimate
!     local rotational diffusion times used to compute Deff
!     to be SVD'd in eventual computation rotational diffusion tensor
!     related material found in J. Biomol. NMR, 9, 287, (1997) L. K. Lee,
!     M. Rance, W. J. Chazin, A. G. Palmer 
!     Science, 268, 886 (1995) R. Bruschweiler, X. Liao, P. E. Wright

!     normalize vectors
      do i=1,nvec
         mag=0d0
         do j=1,3
            mag=mag+x(i,j)*x(i,j)
         end do
         mag=dsqrt(mag)
         do j=1,3
            x(i,j)=x(i,j)/mag
         end do
      end do

      do i=1,nvec
         mat(i,1) = x(i,1) * x(i,1) ! x^2
         mat(i,2) = x(i,2) * x(i,2) ! y^2
         mat(i,3) = x(i,3) * x(i,3) ! z^2
         mat(i,4) = 2 * x(i,1) * x(i,2) ! 2xy
         mat(i,5) = 2 * x(i,2) * x(i,3) ! 2yz
         mat(i,6) = 2 * x(i,1) * x(i,3) ! 2xz
          
!         do j=1,3
!            mat(i,j)=x(i,j)*x(i,j)
!            write(6,'(4(a,i6),a)') 'A[',i,',',j,']=x(',i,',',j,')^2'
!            do k=j+1,3
!               mat(i,j+k+1)=2d0*x(i,j)*x(i,k)
!               write(6,'(6(a,i6),a)') 'A[',i,',',j+k+1,']=x(',i,',',j, &
!                     ')*x(',i,',',k,')'
!            end do
!         end do
      end do

      return
      end 

!     --------------------------------------------------------------------------
      subroutine asymtop(d,pa,r,nvec,nmax,tau1,tau2,sumc2)
      implicit none

!     computes tau(l=1) and tau(l=2) given principal components
!     of rotational diffusion tensor and angular coordinates of
!     a vector (in the PA frame, n*ez = cos(theta),
!     n*ey = sin(theta)*sin(phi), n*ex = sin(theta)*cos(phi)

      integer nmax,nvec
      real*8 d(3),pa(3,3),r(nmax,3)

      real*8 dx,dy,dz,theta,phi,pi,da,ea,epsx,epsy,epsz
      parameter (pi=3.141592653589793d0)

      real*8 sth,cth,sphi,cphi
      real*8 sth2,cth2,sphi2,cphi2,tau1(*),tau2(*),sumc2(*)

      real*8 sth4,s2phi2
      real*8 amp(8),lam(8)

      real*8 s2th2

      real*8 dav,dpr,u,delta,w,n

      real*8 thrcth2m1,thrcth2m12,c2phi,c2phi2

      real*8 x0(3),y0(3),z0(3),x,y,z,mag,dot1,dot2,dot3
      real*8 copy_d(3),copy_pa(3,3)

      integer i,j

!     since it is assumed that Dx <= Dy <= Dz, 
!     diffusion tensor arrays must be sorted
!     leave originals alone to avoid altering downstream results
      do i=1,3
         copy_d(i)=d(i)
         !write(6,'(a,i6,f10.5)') 'D',i,d(i)
         !write(6,'(a,i6,3(f10.5))') 'Vec ',i,pa(i,1),pa(i,2),pa(i,3)
         do j=1,3
            copy_pa(i,j)=pa(i,j)
         end do
      end do
      call sort(3,3,3,copy_d,copy_pa)

      dx=copy_d(1)
      dy=copy_d(2)
      dz=copy_d(3)
      do i=1,3
         x0(i)=copy_pa(1,i)
         y0(i)=copy_pa(2,i)
         z0(i)=copy_pa(3,i)
      end do

      do i=1,3
         mag=x0(i)*x0(i)+y0(i)*y0(i)+z0(i)*z0(i)
      if( mag < 0.d0 ) then
         write(6,*) 'sqrt 3:', mag
         stop
      end if
         mag=dsqrt(mag)
         x0(i)=x0(i)/mag
         y0(i)=y0(i)/mag
         z0(i)=z0(i)/mag
      end do

!     tau = sum(m){amp(l,m)/lambda(l,m)}
!     see Korzhnev DM, Billeter M, Arseniev AS, Orekhov VY;
!     Prog. Nuc. Mag. Res. Spec., 38, 197 (2001) for details
!     only weights need to be computed for each vector, 
!     decay constants can be computed once
!     first three decay constants correspond to l=1, m=-1,0,+1
!     next five correspond to l=2, m=-2,-1,0,+1,+2

!     l=1, m=-1 term:
!     lambda(1,-1) = Dy + Dz
      lam(1)=dy+dz

!     l=1, m=0 term:
!     lambda(1,0) = Dx + Dy
      lam(2)=dx+dy

!     l=1, m=+1 term:
!     lambda(1,+1) = Dx + Dz
      lam(3)=dx+dz

!     l=2, m=-2 term:
!     lambda(2,-2) = Dx + Dy + 4*Dz
      lam(4)=dx+dy+4d0*dz

!     l=2, m=-1 term:
!     lambda(2,-1) = Dx + 4*Dy + Dz
      lam(5)=dx+4d0*dy+dz

!     l=2, m=0 term:
!     lambda(2,0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
!     Dav = (Dx + Dy + Dz)/3, Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
!     u = sqrt(3)*(Dx - Dy)
!     delta = 3*sqrt(Dav*Dav - Dpr*Dpr)
!     w = 2*Dz - Dx - Dy + 2*delta
!     N = 2*sqrt(delta*w)

      dav=(dx+dy+dz)/3d0
      if( dx*dy+dy*dz+dx*dz < 0.d0 ) then
         !  write(6,*) 'sqrt 4:', dx,dy,dz,dx*dy+dy*dz+dx*dz
         dpr = 0.0
      else
         dpr=dsqrt((dx*dy+dy*dz+dx*dz)/3.d0)
      endif
      u=(dx-dy)*dsqrt(3d0)
      if( dav*dav-dpr*dpr < 0.d0 ) then
         write(6,*) 'sqrt 1:', dav, dpr
         stop
      end if
      delta=3d0*dsqrt(dav*dav-dpr*dpr)
      w=2.d0*dz-dx-dy+2.d0*delta
      if( w*delta < 0.d0 ) then
         write(6,*) 'sqrt 2:', w,delta
         stop
      end if
      n=2d0*dsqrt(delta*w)

      lam(6)=dav-dsqrt(dav*dav-dpr*dpr)
      lam(6)=6d0*lam(6)

!     l=2, m=+1 term:
!     lambda(2,+1) = 4*Dx + Dy + Dz
      lam(7)=4d0*dx+dy+dz

!     l=2, m=+2 term:
!     lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav = Dpr*Dpr)]
      lam(8)=dav+dsqrt(dav*dav-dpr*dpr)
      lam(8)=6d0*lam(8)

      do i=1,nvec
         x=r(i,1)
         y=r(i,2)
         z=r(i,3)
   
         mag=x*x+y*y+z*z
         if( mag < 0.d0 ) then
            write(6,*) 'sqrt 5:', mag
            stop
         end if
         mag=dsqrt(mag)
         x=x/mag
         y=y/mag
         z=z/mag

         dot1=x*x0(1)+y*y0(1)+z*z0(1)
         dot2=x*x0(2)+y*y0(2)+z*z0(2)
         dot3=x*x0(3)+y*y0(3)+z*z0(3)
      !write(6,'(a,3(f10.5))') 'DBG: dot1-3= ',dot1,dot2,dot3

!     assuming e(3)*n = cos(theta), e(1)*n = sin(theta)*cos(phi),
!     e(2)*n = sin(theta)*sin(phi), theta >= 0;
!     sin(theta) = sqrt(1 - cos(theta)^2) and theta = tan^-1[sin(theta)/cos(theta)];
!     phi = tan^-1[sin(phi)/cos(phi)] = tan^-1[(e(2)*n)/(e(1)*n)]

      if( 1d0-dot3*dot3 < 0.d0 ) then
         write(6,*) 'sqrt 6:', dot3
         stop
      end if
         theta=datan2(dsqrt(1d0-dot3*dot3),dot3)
         phi=datan2(dot2,dot1)

!     theta=theta*(pi/180d0)
!     phi=phi*(pi/180d0)
         sth=dsin(theta)
         cth=dcos(theta)
         sphi=dsin(phi)
         cphi=dcos(phi)

!     tau(l=1) = [sin(theta)^2*cos(phi)^2]/(Dy+Dz) +
!                + [sin(theta)^2*sin(phi)^2]/(Dx+Dz) +
!                + [cos(theta)^2]/(Dx+Dy)
!     compute correlation time for l=1
         sth2=sth*sth
         cth2=cth*cth
         sphi2=sphi*sphi
         cphi2=cphi*cphi
         amp(1)=sth2*cphi2
         amp(2)=cth2
         amp(3)=sth2*sphi2

         tau1(i)=0d0
         do j=1,3
            tau1(i)=tau1(i)+amp(j)/lam(j)
         end do

!     m=-2 term:
!     lambda(2,-2) = Dx + Dy + 4*Dz
!     amp(2,-2) = 0.75*sin(theta)^4*sin(2*phi)^2 = 3l^2m^2
!        sth4=sth2*sth2
!        s2phi2=dsin(2d0*phi)
!        s2phi2=s2phi2*s2phi2
!        amp(4)=0.75d0*sth4*s2phi2
         amp(4) = 3.d0*dot1*dot1*dot2*dot2

!     m=-1 term:
!     lambda(2,-1) = Dx + 4*Dy + Dz
!     amp(2,-1) = 0.75*sin(2*theta)^2*cos(phi)^2 = 3l^2n^2
!        s2th2=dsin(2d0*theta)
!        s2th2=s2th2*s2th2
!        amp(5)=0.75d0*s2th2*cphi2
         amp(5) = 3.d0*dot1*dot1*dot3*dot3

!     m=0 term:
!     lambda(2,0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
!     Dav = (Dx + Dy + Dz)/3, Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
!     amp(2,0) = (w/N)^2*0.25*[3*cos(theta)^2-1]^2 +
!     (u/N)^2*0.75*sin(theta)^4*cos(2*phi)^2 -
!     (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
!     u = sqrt(3)*(Dx - Dy)
!     delta = 3*sqrt(Dav*Dav - Dpr*Dpr)
!     w = 2*Dz - Dx - Dy + 2*delta
!     N = 2*sqrt(delta*w)
!        thrcth2m1=3d0*cth2-1d0
!        thrcth2m12=thrcth2m1*thrcth2m1
!        c2phi=dcos(2d0*phi)
!        c2phi2=c2phi*c2phi
!        amp(6)=(w/n)*(w/n)*0.25d0*thrcth2m12
!        amp(6)=amp(6)+(u/n)*(u/n)*0.75d0*sth4*c2phi2
!        amp(6)=amp(6)-(u/delta)*0.125d0*dsqrt(3d0)*thrcth2m1*sth2*c2phi
         da = 0.25*(3.d0*(dot1**4 + dot2**4 + dot3**4) -1.d0)
         if( delta > 1.d-8) then
            epsx = 3.d0*(Dx-Dav)/delta
            epsy = 3.d0*(Dy-Dav)/delta
            epsz = 3.d0*(Dz-Dav)/delta
            ea = epsx*(3.d0*dot1**4 + 6.d0*(dot2*dot3)**2 -1.d0) &
               + epsy*(3.d0*dot2**4 + 6.d0*(dot1*dot3)**2 -1.d0) &
               + epsz*(3.d0*dot3**4 + 6.d0*(dot1*dot2)**2 -1.d0) 
            ea = ea/12.d0
         else
            ea = 0.d0
         endif
         amp(6) = da + ea

!     m=+1 term:
!     lambda(2,+1) = 4*Dx + Dy + Dz
!     amp(2,+1) = 0.75*sin(2*theta)^2*sin(phi)^2
!        amp(7)=0.75d0*s2th2*sphi2
         amp(7)=3.d0*dot2*dot2*dot3*dot3

!     m=+2 term:
!     lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav = Dpr*Dpr)]
!     amp(2,+2) = (u/n)^2*0.25*[3*cos(theta)^2-1]^2 +
!     (w/n)^2*0.75*sin(theta)^4*cos(2*phi)^2 +
!     (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
!        amp(8)=(u/n)*(u/n)*0.25d0*thrcth2m12
!        amp(8)=amp(8)+(w/n)*(w/n)*0.75d0*sth4*c2phi2
!        amp(8)=amp(8)+(u/delta)*0.125d0*dsqrt(3d0)*thrcth2m1*sth2*c2phi
         amp(8)=da-ea

         tau2(i)=0d0
         sumc2(i)=0.d0
         do j=4,8
            tau2(i)=tau2(i)+amp(j)/lam(j)
            sumc2(i)=sumc2(i)+amp(j)
         end do
      end do

      return
      end
              
!     --------------------------------------------------------------------------
      subroutine sort(n,xdim,ydim,x,y)
      implicit none
      integer i,j,k,l,m,n,o,ln2n,fl,xdim,ydim
      real*8 x(xdim),temp
      real*8 y(ydim,ydim)

      integer p

      ln2n=int(dlog(dble(n))/0.693147181+1d-5)
      m=n
      do o=1,ln2n
        m=m/2
        k=n-m
        do j=1,k
          i=j
          fl=1
          do while(i.ge.1.and.fl.eq.1)
            l=i+m
            fl=0
            if(x(l).lt.x(i)) then
              temp=x(i)
              x(i)=x(l)
              x(l)=temp

!             additional 2-D array to sort in same order as array x
              do p=1,n
                 temp=y(p,i)
                 y(p,i)=y(p,l)
                 y(p,l)=temp
              end do

              i=i-m
              fl=1
            endif
          enddo
        enddo
      enddo
      return
      end subroutine sort

