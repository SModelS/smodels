CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C        SPIRIX' SUBROUTINE ALPHA_S AND THE RUNNING QUARK MASSES       C
C                                                                      C
C           ALSINI(ACC,LAMBDA,MC,MB,MT,NF)                             C
C           RUNM_ORIG(Q,NF)                                            C
C           ALPHAS(Q,N)                                                C
C           XITER_ORIG(*) ONLY INTERNALLY USED                         C
C                                                                      C
C        FIRST CALL ALSINI FOR ALPHAS AND FOR RUNM_ORIG !!!            C
C                                                                      C
C        * NO EXTERNAL COMMON BLOCKS                                   C
C        * INTERNAL COMMON BLOCKS: ALS_ORIG [SET BY ALSINI]            C
C                                  ALSLAM_ORIG [SET BY ALSINI]         C
C        * STRANGE QUARK MASS SET INSIDE THE PROGRAM                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ALSINI(acc,xlambda_in,amc_in,amb_in,amt_in,n0_in)
ctp      SUBROUTINE ALSINI(ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM_ORIG/XLB1(6),XLB2(6)
      COMMON/ALS_ORIG/XLAMBDA,AMC,AMB,AMT,N0
ctp [begin]
c               fill the common input ALS by arguments
      xlambda = xlambda_in
      amc = amc_in
      amb = amb_in
      amt = amt_in
      n0  = n0_in
ctp [end]      
ctp      PI=4.D0*ATAN(1.D0)
      XLB1(1)=0D0
      XLB1(2)=0D0
      XLB2(1)=0D0
      XLB2(2)=0D0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
      ENDIF
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*LOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER_ORIG(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*LOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER_ORIG(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*LOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_ORIG(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*LOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER_ORIG(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*LOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_ORIG(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*LOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_ORIG(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*LOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER_ORIG(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*LOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_ORIG(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*LOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_ORIG(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*LOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER_ORIG(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*LOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER_ORIG(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*LOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_ORIG(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION RUNM_ORIG(Q,NF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/ALS_ORIG/XLAMBDA,AMCA,AMBA,AMTA,N0A
ctp      COMMON/QMASSES/AMS,AMC,AMB,AMT
      SAVE ISTRANGE
      B0(NF)=(33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X,XK)=1D0+4D0/3D0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      PI=4D0*ATAN(1D0)
ctp [begin]
      ams  = .3D0
      amc  = amca
      amb  = amba 
      amt  = amta 
ctp [end]
      AMSB = AMS
      NNLO = 0
      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
C--------------------------------------------
      IMSBAR = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(ALPHAS(AMS,2)/PI,3)
     .            *CQ(ALPHAS(1.D0,2)/PI,3)/TRAN(AMS,0D0)
        DD = (XMSB-AMSB)/AMSB
        IF(ABS(DD).GE.ACC)THEN
         IF(DD.LE.0.D0)THEN
          AMSD = AM(3)
         ELSE
          AMSU = AM(3)
         ENDIF
         GOTO 123
        ENDIF
        ISTRANGE=1
       ENDIF
       AM(3) = AMSB
      ELSE
       AMS=AMSB
       AM(3) = AMS
      ENDIF
C--------------------------------------------
      AM(3) = AMSB
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT
      XK = 16.11D0
      DO 1 I=1,NF-1
       XK = XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF),0D0)
ctp       XMHAT = XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0
ctp       XMHAT = 0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .                   CQ(ALPHAS(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .                   CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      IF(NNLO.EQ.1.AND.NF.GT.3)THEN
       XKFAC = TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC = 1D0
      ENDIF
      RUNM_ORIG = YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .               CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM_ORIG/XLB1(6),XLB2(6)
      COMMON/ALS_ORIG/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS1(NF,X)=12.D0*PI/(B0(NF)*LOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*LOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*LOG(LOG(X**2/XLB(NF)**2))
     .           /LOG(X**2/XLB(NF)**2))
      PI=4.D0*ATAN(1.D0)
      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1      CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2      CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMB)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSE
        ALPHAS=ALS2(NF,Q)
      ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION XITER_ORIG(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*LOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*LOG(LOG(X**2/XLB**2))
     .              /LOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+SQRT(1D0-4D0*B*LOG(X)))
      PI=4.D0*ATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=LOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2.D0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=ABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER_ORIG=XLB2
      RETURN
      END

