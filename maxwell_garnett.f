c      subroutine maxwell_garnett(max_nwave,max_ninc,nwave,ninc,e_mat,e_inc,vol_frac,L1,L2,
c     *                           tol,niter,e_mg,iter,error,chisq) 
c
c     Purpose: Compute the effective permeability using the generalized
c              Maxwell-Garnett formulation for abitrary shape and number
c              of inclusions types
c
c     Warning: This is a research model/code. Please use with caution.
c
c     Reference: Granqvist C. G. & Hunderi O. Phys. Rev. B, 16, 3513 (1977)
c                Bohren & Huffman
c
c     Version: 0.1 May 18 2004 Wing-Fai Thi
c              0.2 May 19 2004 Wing-Fai Thi use max_nwave and max_ninc to define size of arrays
c
c     License GNU v 3.0
c
c     input
c     max_nwave             max size for the wave array
c     max_ninc              max size for the inc array
c     nwave                 total number (integer) of wavelength bins
c     ninc                  total number (integer) of inclusions types
c     e_mat(1:max_nwave)        complex dielectric function of the matrix 
c     e_inc(1:max_ninc,1:max_nwave) complex dielectric function of the inclusions
c     vol_frac(1:max_ninc)      volumn vol_fraction (double) of each inclusion type
c     L1(1:max_ninc)            polarization factor in the first direction; sphere L1=L2=L3=1.0D0/3.0D0
c     L2(1:max_ninc)            polarization factor in the second direction
c
c     tol                   tolerance for convergence (~1d-2)
c     niter                 maximum number of iterations (integer~100)       
c
c     output
c     e_mg(1:max_nwave)   complex effective permeability of the effective medium
c                     (dielectric function)
c     error           error handeling variables error = 0 sucess
c                                               error = 1 not converging
c                                               error = 2 error in the input data
c
c     local variables
c     alpha1(1:max_nwave) complex polarizability of inclusion inc polarization 1
c     alpha2(1:max_nwave) complex polarizability of inclusion inc polarization 2
c     alpha3(1:max_nwave) complex polarizability of inclusion inc polarization 3
c     alpha(1:max_nwave)  total complex polarizability of inclusion inc
c     L3(1:max_ninc)             polarization factor in the third direction L3=1.0d0-L1-L2
c     e_mg(1:max_nwave)          local complex effective permeability use for the iteration
c     inc                    inclusion loop index
c     wave                   wavelength loop index
c     iter                   iteration loop index
c     chisq                  double chi square
c     chisq2                 temporary value of chi square
c     mchisq                 relative difference in chisq between two iterations
c     vol_tot                total volumn fration of the inclusions
c 
      subroutine maxwell_garnett(max_nwave,max_ninc,
     *                           nwave,ninc,e_mat,e_inc,vol_frac,L1,L2,
     *                           tol,niter,e_mg,iter,error,chisq) 
      implicit none
c --- parameters
      integer max_nwave,max_ninc
c
c --- define variables
c     input
      integer nwave,ninc,niter,error
      double precision L1(1:max_ninc),L2(1:max_ninc),
     *                 vol_frac(1:max_ninc),tol

      complex e_mat(1:max_nwave),e_inc(1:max_ninc,1:max_nwave)
c     output
      integer iter
      complex e_mg(1:max_nwave)
c     local variables
      integer inc,wave
      double precision chisq,chisq2,mchisq,L3(1:max_ninc),min_chisq,
     *   vol_tot
      complex alpha1(1:max_nwave),alpha2(1:max_nwave),
     *        alpha3(1:max_nwave),alpha(1:max_nwave),
     *        sum_inc,e_mg2(1:max_nwave)
      parameter(min_chisq=1d-30)

c --- declare variables
      error = 0
      chisq = 100.0d0
      vol_tot = 0.0d0
c --- compute total volumn occupied by the inclusions
      do inc=1,ninc
         vol_tot = vol_tot + vol_frac(inc)
      enddo   
      if (vol_tot.gt.1.0d0) then
         error = 2
         goto 2
      endif   
c --- compute the polarizabilities with as first guess (e_mg=e_mat)

      do wave = 1,nwave
         sum_inc = (0.0,0.0)
         do inc = 1,ninc
c           L1 + L2 + L3 = 1.0
            L3(inc) = 1.0d0 - L2(inc) - L1(inc)
            if (L3(inc).lt.0.0d0) then
               error = 2
               goto 2
            endif   
c
            alpha1(wave) = (e_inc(inc,wave)-e_mat(wave))/
     *   (e_mat(wave)+L1(inc)*(e_inc(inc,wave)-e_mat(wave)))
            alpha2(wave) = (e_inc(inc,wave)-e_mat(wave))/
     *   (e_mat(wave)+L2(inc)*(e_inc(inc,wave)-e_mat(wave)))
            alpha3(wave) = (e_inc(inc,wave)-e_mat(wave))/
     *   (e_mat(wave)+L3(inc)*(e_inc(inc,wave)-e_mat(wave)))
            alpha(wave) = (alpha1(wave)+alpha2(wave)+alpha3(wave))/3.0d0
            sum_inc = sum_inc + vol_frac(inc)*alpha(wave)
         enddo
         e_mg(wave) = e_mat(wave)*((1.0d0+(2.0d0/3.0d0)*sum_inc)/
     *                             (1.0d0-(1.0d0/3.0d0)*sum_inc))

      enddo   
c --- the actual iteration starts here
      do iter=1,niter
         chisq2= 0.0d0
         do wave = 1,nwave
           sum_inc = (0.0,0.0)
           do inc = 1,ninc
            alpha1(wave) = (e_inc(inc,wave)-e_mg(wave))/
     *   (e_mat(wave)+L1(inc)*(e_inc(inc,wave)-e_mg(wave)))
            alpha2(wave) = (e_inc(inc,wave)-e_mat(wave))/
     *   (e_mat(wave)+L2(inc)*(e_inc(inc,wave)-e_mg(wave)))
            alpha3(wave) = (e_inc(inc,wave)-e_mat(wave))/
     *   (e_mat(wave)+L3(inc)*(e_inc(inc,wave)-e_mg(wave)))
            alpha(wave) = (alpha1(wave)+alpha2(wave)+alpha3(wave))/3.0d0
            sum_inc = sum_inc + vol_frac(inc)*alpha(wave)
           enddo
           e_mg2(wave) = e_mat(wave)*((1.0d0+(2.0d0/3.0d0)*sum_inc)/
     *                               (1.0d0-(1.0d0/3.0d0)*sum_inc))
           chisq2 = chisq2+(abs(e_mg(wave)-e_mg2(wave))**2.0d0)
           e_mg(wave) = e_mg2(wave)
        enddo
           chisq2 = chisq2 / nwave
           mchisq = abs(chisq-chisq2)/(0.5d0*(chisq+chisq2))
           chisq  = chisq2
           if (mchisq.le.tol.or.chisq.lt.min_chisq) goto 1
      enddo
c     convergence not reached
      error = 1
c
 1    continue
      do wave = 1,nwave
         e_mg(wave) = e_mg2(wave)
      enddo   
c
 2    continue
      return
      end
c
