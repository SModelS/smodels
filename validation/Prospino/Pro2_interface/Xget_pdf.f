cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c GET_PDF(inlo,x,mu,pdf)                                               c
c        interface calling the parton densities                        c
c                                                                      c
c   input  inlo=0,1 for LO, NLO set, dummy parameter for CTEQ calls    c
c          x        momentum fraction (cut for safety)                 c
c          mu       factorization scale                                c
c                                                                      c
c   output pdf(-6:6) just like CteQ: g,u,d,s,c,b,(t)                   c
c                                                                      c
c n.b. lambda_qcd has to be set accordingly in GET_LAMBDA_QCD          c
c      cteq arrays have to initialized using INIT_PDF                  c
c                                                                      c
c If for some reason you do not like Cteq5, you can switch to Cteq6:   c
c   (1) change the actual pdf call below                               c
c   (2) change the lambda_QCD values below                             c
c   (3) change the initialization call below                           c
c   note for (1)-(3) the Cteq6 values are commented out using `ctq6'   c 
c   (4) link the proper routine in the Makefile                        c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GET_PDF(inlo,x,mu,pdf)

      implicit none
      
      integer inlo,i1
      real*8  x,mu,pdf(-6:6),Ctq6Pdf,Ctq5Pdf

      logical ldebug
      parameter( ldebug=.false. )

      if (ldebug) print*, " GET_PDF: input parameters ",x,mu

c              initialize the array 
      do i1=-6,6,1
         pdf(i1) = 0.D0
      end do

c              check if range makes sense, otherwise return zero
      if ( (x.ge.1.D0) .or. (x.lt.1.D-5) ) then
         print*, " GET_PDF: problem with x ",x
         return
      end if

c              if scale too low, return zero
      if (mu.lt.5.D0) then 
         print*, " GET_PDF: problem with too small muF ",mu
         return
      end if

c              if scale too high, set to maximum
      if (mu.gt.1.5D3) then
         print*, " GET_PDF: problem with too large muF ",mu
         mu = 1.5D3
      end if
         
     
      do i1=-5,5,1
         pdf(i1) = Ctq6Pdf(i1,x,mu)
ctq5         pdf(i1) = Ctq5Pdf(i1,x,mu)
         if (ldebug) print*, " GET_PDF: pdf call ",i1,pdf(i1)
      end do
   
      return
      end

c -----------------------------------------------
c   this initialization is needed for a consistent alpha_s
      subroutine GET_LAMBDA_QCD(inlo,lambda_qcd)
      
      implicit none

      integer inlo
      real*8  lambda_qcd

      if (inlo.eq.0) then 
         lambda_qcd = 0.165D0
ctq5         lambda_qcd = 0.146D0
      else if (inlo.eq.1) then 
         lambda_qcd = 0.226D0
ctq5         lambda_qcd = 0.226D0
      else
         print*, " GET_LAMBDA_QCD: inlo set wrongly ",inlo
      end if
         
      return
      end 

c -----------------------------------------------
c   this initialization function is only needed for Cteq grids 
      subroutine INIT_PDF(inlo)

      implicit none

      integer inlo

      if (inlo.eq.0) then 
         call SetCtq6(4)
ctq5         call SetCtq5(3)
      else if (inlo.eq.1) then 
         call SetCtq6(400)
ctq6         call SetCtq6(1)
ctq5         call SetCtq5(8)
      else
         print*, " INIT_PDF: inlo set wrongly ",inlo
      end if
         
      return
      end 


