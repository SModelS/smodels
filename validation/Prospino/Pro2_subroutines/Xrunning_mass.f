CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C        SPIRIX' SUBROUTINE ALPHA_S AND THE RUNNING QUARK MASSES       C
C        tp: make it stand alone                                       C
C        tp: chose the order of running masses correctly               C
C                                                                      C
C           ALSINI_RUNM(ACC,LAMBDA,MC,MB,MT,NF)                        C
C           RUNM_EXT(Q,NF)                                             C
C           ALPHAS_RUNM(Q,N)                                           C
C           XITER_ORIG(*) ONLY INTERNALLY USED, NOT INCLUDED           C
C                                                                      C
C        FIRST CALL ALSINI_RUNM FOR ALPHAS_RUNM AND FOR RUNM_EXT !!!   C
C                                                                      C
C        * NO EXTERNAL COMMON BLOCKS                                   C
C        * INTERNAL COMMON BLOCKS: ALS,ALSLAM [SET BY ALSINI_RUNM]     C
C        * STRANGE QUARK MASS SET INSIDE THE PROGRAM                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ALSINI_RUNM(acc,xlambda_in,amc_in,amb_in,amt_in,n0_in)

      implicit none

      integer n0_in,i
      real*8  acc,xlambda_in,amc_in,amb_in,amt_in,xlb(6),XITER_ORIG

      real*8             xlb1(6),xlb2(6)
      COMMON/ALSLAM_RUNM/xlb1   ,xlb2

      integer n0
      real*8          xlambda,amc,amb,amt
      COMMON/ALS_RUNM/xlambda,amc,amb,amt,n0
ctp [begin]
c               fill the common input als by arguments
      xlambda = xlambda_in
      amc = amc_in
      amb = amb_in
      amt = amt_in
      n0  = n0_in
ctp [end]      
c               compute the different values of lambda_qcd 
      xlb1(1)=0d0
      xlb1(2)=0d0
      xlb2(1)=0d0
      xlb2(2)=0d0

      if(n0.eq.3)then
         xlb(3)=xlambda
         xlb(4)=xlb(3)*(xlb(3)/amc)**(2.d0/25.d0)
         xlb(5)=xlb(4)*(xlb(4)/amb)**(2.d0/23.d0)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
      elseif(n0.eq.4)then
         xlb(4)=xlambda
         xlb(5)=xlb(4)*(xlb(4)/amb)**(2.d0/23.d0)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
      elseif(n0.eq.5)then
         xlb(5)=xlambda
         xlb(4)=xlb(5)*(xlb(5)/amb)**(-2.d0/25.d0)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
      elseif(n0.eq.6)then
         xlb(6)=xlambda
         xlb(5)=xlb(6)*(xlb(6)/amt)**(-2.d0/23.d0)
         xlb(4)=xlb(5)*(xlb(5)/amb)**(-2.d0/25.d0)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
      endif
      
      do i=1,6
         xlb1(i)=xlb(i)
      end do
      
      if(n0.eq.3)then
         xlb(3)=xlambda
         xlb(4)=xlb(3)*(xlb(3)/amc)**(2.d0/25.d0)
     .                *(2.d0*log(amc/xlb(3)))**(-107.d0/1875.d0)
         xlb(4)=xiter_orig(amc,xlb(3),3,xlb(4),4,acc)
         xlb(5)=xlb(4)*(xlb(4)/amb)**(2.d0/23.d0)
     .                *(2.d0*log(amb/xlb(4)))**(-963.d0/13225.d0)
         xlb(5)=xiter_orig(amb,xlb(4),4,xlb(5),5,acc)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
     .                *(2.d0*log(amt/xlb(5)))**(-321.d0/3381.d0)
         xlb(6)=xiter_orig(amt,xlb(5),5,xlb(6),6,acc)
      elseif(n0.eq.4)then
         xlb(4)=xlambda
         xlb(5)=xlb(4)*(xlb(4)/amb)**(2.d0/23.d0)
     .                *(2.d0*log(amb/xlb(4)))**(-963.d0/13225.d0)
         xlb(5)=xiter_orig(amb,xlb(4),4,xlb(5),5,acc)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
     .                *(2.d0*log(amc/xlb(4)))**(107.d0/2025.d0)
         xlb(3)=xiter_orig(amc,xlb(4),4,xlb(3),3,acc)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
     .                *(2.d0*log(amt/xlb(5)))**(-321.d0/3381.d0)
         xlb(6)=xiter_orig(amt,xlb(5),5,xlb(6),6,acc)
      elseif(n0.eq.5)then
         xlb(5)=xlambda
         xlb(4)=xlb(5)*(xlb(5)/amb)**(-2.d0/25.d0)
     .                *(2.d0*log(amb/xlb(5)))**(963.d0/14375.d0)
         xlb(4)=xiter_orig(amb,xlb(5),5,xlb(4),4,acc)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
     .                *(2.d0*log(amc/xlb(4)))**(107.d0/2025.d0)
         xlb(3)=xiter_orig(amc,xlb(4),4,xlb(3),3,acc)
         xlb(6)=xlb(5)*(xlb(5)/amt)**(2.d0/21.d0)
     .                *(2.d0*log(amt/xlb(5)))**(-321.d0/3381.d0)
         xlb(6)=xiter_orig(amt,xlb(5),5,xlb(6),6,acc)
      elseif(n0.eq.6)then
         xlb(6)=xlambda
         xlb(5)=xlb(6)*(xlb(6)/amt)**(-2.d0/23.d0)
     .                *(2.d0*log(amt/xlb(6)))**(321.d0/3703.d0)
         xlb(5)=xiter_orig(amt,xlb(6),6,xlb(5),5,acc)
         xlb(4)=xlb(5)*(xlb(5)/amb)**(-2.d0/25.d0)
     .                *(2.d0*log(amb/xlb(5)))**(963.d0/14375.d0)
         xlb(4)=xiter_orig(amb,xlb(5),5,xlb(4),4,acc)
         xlb(3)=xlb(4)*(xlb(4)/amc)**(-2.d0/27.d0)
     .                *(2.d0*log(amc/xlb(4)))**(107.d0/2025.d0)
         xlb(3)=xiter_orig(amc,xlb(4),4,xlb(3),3,acc)
      endif

      do i=1,6
         xlb2(i)=xlb(i)
      end do
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function RUNM_EXT(q,nf,n)
      
      implicit none 
      
      integer nf,n,nn,n0,i
      real*8  q,zeta3,amsb,amc,amb,amt
     &       ,x,xk,xkfac,xmsb,pi,q0,nnlo
     &       ,B0,B1,B2,G0,G1,G2,C1,C2,TRAN,cq
     &       ,ALPHAS_RUNM

      parameter (nn=6)
      parameter (zeta3 = 1.202056903159594d0)
      real*8    am(nn),ymsb(nn)

      integer    n0a
      real*8          xlambda,amca,amba,amta
      COMMON/ALS_RUNM/xlambda,amca,amba,amta,n0a

      B0(nf)=(33.d0-2.d0*nf)/12d0
      B1(nf) = (102d0-38d0/3d0*nf)/16d0
      B2(nf) = (2857d0/2d0-5033d0/18d0*nf+325d0/54d0*nf**2)/64d0
      G0(nf) = 1d0
      G1(nf) = (202d0/3d0-20d0/9d0*nf)/16d0
      G2(nf) = (1249d0-(2216d0/27d0+160d0/3d0*zeta3)*nf
     .       - 140d0/81d0*nf**2)/64d0
      C1(nf,n) = dmin1( dble(n-1), 1.D0) 
     .          *( G1(nf)/B0(nf) - B1(nf)*G0(nf)/B0(nf)**2 )
      C2(nf,n) = dmax1( dble(n-2), 0.D0)
     .          *( ((G1(nf)/B0(nf) - B1(nf)*G0(nf)/B0(nf)**2)**2
     .       + G2(nf)/B0(nf) + B1(nf)**2*G0(nf)/B0(nf)**3
     .       - B1(nf)*G1(nf)/B0(nf)**2 - B2(nf)*G0(nf)/B0(nf)**2)/2d0 )
      TRAN(x,xk,n)=1d0+4d0/3d0*ALPHAS_RUNM(x,n)/pi
     .              +xk*(ALPHAS_RUNM(x,n)/pi)**2
      CQ(x,nf,n)=(2d0*B0(nf)*x)**(G0(nf)/B0(nf))
     .            *(1d0+C1(nf,n)*x+C2(nf,n)*x**2)

      pi=4d0*atan(1d0)
ctp [begin]
      amsb = .3d0
      amc  = amca
      amb  = amba 
      amt  = amta 
ctp [end]
      if (n.lt.3) then 
         nnlo = 0
      else if (n.eq.3) then 
         nnlo = 1
      end if 

      am(1) = 0
      am(2) = 0
      am(3) = amsb
      am(4) = amc
      am(5) = amb
      am(6) = amt
      
      xk = 16.11d0
      do i=1,nf-1
         xk = xk - 1.04d0*(1.d0-am(i)/am(nf))
      end do
      
      if(nf.ge.4)then
         xmsb = am(nf)/TRAN(am(nf),0d0,n)
      else
         xmsb = 0
      endif
      ymsb(3) = amsb
      if(nf.eq.3)then
         ymsb(4) = ymsb(3)*CQ(ALPHAS_RUNM(am(4),n)/pi,3,n)/
     .                     CQ(ALPHAS_RUNM(1.d0 ,n)/pi,3,n)
         ymsb(5) = ymsb(4)*CQ(ALPHAS_RUNM(am(5),n)/pi,4,n)/
     .                     CQ(ALPHAS_RUNM(am(4),n)/pi,4,n)
         ymsb(6) = ymsb(5)*CQ(ALPHAS_RUNM(am(6),n)/pi,5,n)/
     .                     CQ(ALPHAS_RUNM(am(5),n)/pi,5,n)
      elseif(nf.eq.4)then
         ymsb(4) = xmsb
         ymsb(5) = ymsb(4)*CQ(ALPHAS_RUNM(am(5),n)/pi,4,n)/
     .                     CQ(ALPHAS_RUNM(am(4),n)/pi,4,n)
         ymsb(6) = ymsb(5)*CQ(ALPHAS_RUNM(am(6),n)/pi,5,n)/
     .                     CQ(ALPHAS_RUNM(am(5),n)/pi,5,n)
      elseif(nf.eq.5)then
         ymsb(5) = xmsb
         ymsb(4) = ymsb(5)*CQ(ALPHAS_RUNM(am(4),n)/pi,4,n)/
     .                     CQ(ALPHAS_RUNM(am(5),n)/pi,4,n)
         ymsb(6) = ymsb(5)*CQ(ALPHAS_RUNM(am(6),n)/pi,5,n)/
     .                     CQ(ALPHAS_RUNM(am(5),n)/pi,5,n)
      elseif(nf.eq.6)then
         ymsb(6) = xmsb
         ymsb(5) = ymsb(6)*CQ(ALPHAS_RUNM(am(5),n)/pi,5,n)/
     .                     CQ(ALPHAS_RUNM(am(6),n)/pi,5,n)
         ymsb(4) = ymsb(5)*CQ(ALPHAS_RUNM(am(4),n)/pi,4,n)/
     .                     CQ(ALPHAS_RUNM(am(5),n)/pi,4,n)
      endif
      if(q.lt.amc)then
         n0=3
         q0 = 1.d0
      elseif(q.le.amb)then
         n0=4
         q0 = amc
      elseif(q.le.amt)then
         n0=5
         q0 = amb
      else
         n0=6
         q0 = amt
      endif
      if(nnlo.eq.1.and.nf.gt.3)then
         xkfac = TRAN(am(nf),0d0,n)/TRAN(am(nf),xk,n)
      else
         xkfac = 1d0
      endif
      RUNM_EXT = ymsb(n0)*CQ(ALPHAS_RUNM(q,n)/pi,n0,n)/
     .                    CQ(ALPHAS_RUNM(q0,n)/pi,n0,n)
     .       * xkfac
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ALPHAS_RUNM(q,n)

      implicit none

      integer n,i,nf
      real*8  q,B0,B1,ALS1,ALS2,pi,x,xlb(6)
      
      real*8             xlb1(6),xlb2(6)
      COMMON/ALSLAM_RUNM/xlb1   ,xlb2

      integer n0
      real*8          xlambda,amc,amb,amt
      COMMON/ALS_RUNM/xlambda,amc,amb,amt,n0

      B0(nf)=33.d0-2.d0*nf
      B1(nf)=6.d0*(153.d0-19.d0*nf)/B0(nf)**2
      ALS1(nf,x)=12.d0*pi/(B0(nf)*log(x**2/xlb(nf)**2))
      ALS2(nf,x)=12.d0*pi/(B0(nf)*log(x**2/xlb(nf)**2))
     .          *(1.d0-B1(nf)*log(log(x**2/xlb(nf)**2))
     .           /log(x**2/xlb(nf)**2))
      pi=4.d0*atan(1.d0)
      if (n.eq.1) then
         do i=1,6
            xlb(i)=xlb1(i)
         end do
      else if (n.gt.1) then 
         do i=1,6
            xlb(i)=xlb2(i)
         end do
      else 
         print*, " ALPHAS_RUNM: order not set correctly ",n
         stop
      end if

      if(q.lt.amc)then
       nf=3
      elseif(q.le.amb)then
       nf=4
      elseif(q.le.amt)then
       nf=5
      else
       nf=6
      endif

      if (n.eq.1) then
         ALPHAS_RUNM = ALS1(nf,q)
      else if (n.gt.1) then 
         ALPHAS_RUNM = ALS2(nf,q)
      else 
         print*, " ALPHAS_RUNM: order not set correctly ",n
         stop
      endif

      return
      end


