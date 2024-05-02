      program driver_mg
c      
c     Test program for the different EMT routines
c
c     Warning: This is a research model/code. Please use with caution.
c              There is no guarantee that the results are correct.
c
c     License GNU v 3.0
c
      implicit none
c
c --- define variables
      integer max_nwave,max_ninc
      parameter(max_nwave=1000,max_ninc=10)
c
      integer nwave,ninc,iter,niter,error,i
      double precision L1(1:max_ninc),L2(1:max_ninc),
     *                 vol_frac(1:max_ninc),tol,chisq
      complex e_mat(1:max_nwave),e_inc(1:max_ninc,1:max_nwave)
      complex e_eff(1:max_nwave)
c     test data
      nwave = 1
      ninc  = 2
      e_mat(1)   = (1.0,0.2)
      e_inc(1,1) = (1.1,0.4)
      e_inc(2,1) = (1.1,0.1)
      vol_frac(1)= 0.02d0
      vol_frac(2)= 0.05d0
      L1(1)      = 0.3
      L2(1)      = 0.3
      L1(2)      = 0.3
      L2(2)      = 0.3
      tol        = 1.0d-3
      niter      = 100
c      
c     log of a complex number = principal value
c
      write(*,*) 'matrix =',e_mat(1)
      do i=1,ninc
	write(*,*) 'inclusion ',i,' =',e_inc(i,1)
      enddo
      write(*,*) 'vol_frac(1) =',vol_frac(1)
      write(*,*) 'vol_frac(2) =',vol_frac(2)
	write(*,*) 'L1(1) =',L1(1)
      write(*,*) 'L2(1) =',L2(1) 
	write(*,*) 'L1(2) =',L1(2)
      write(*,*) 'L2(2) =',L2(2) 

      call maxwell_garnett(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc
     *  ,vol_frac,L1,L2,tol,niter,e_eff,iter,error,chisq) 
      print*,'Generalized Maxwell-Garnett:',
     *        error,iter,chisq,e_eff(1)
c
      call bruggeman(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,
     *                     L1,L2,tol,niter,e_eff,iter,error,chisq) 

      print*,'Bruggeman ellipsoids:',error,iter,chisq,e_eff(1)
c
      call hunderi(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,
     *                     L1,L2,tol,niter,e_eff,iter,error,chisq) 

      print*,'Hunderi ellipsoids:',error,iter,chisq,e_eff(1)
c
      call mg_cde(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,
     *            tol,niter,e_eff,iter,error,chisq)
      print*,'MG CDE:',error,iter,chisq,e_eff(1)
c
      call br_cde(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,
     *            tol,niter,e_eff,iter,error,chisq)
      print*,'Br CDE:',error,iter,chisq,e_eff(1)
c
      call hu_cde(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,
     *            tol,niter,e_eff,iter,error,chisq)
      print*,'Hunderi CDE:',error,iter,chisq,e_eff(1)
c
      end
