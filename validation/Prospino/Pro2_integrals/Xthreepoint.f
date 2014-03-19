c$$$************************************************************************
c$$$        FUNCTION C03(P1,P2,P3,M1,M2,M3)
c$$$************************************************************************
c$$$*  SCALAR 3-POINT FUNCTION                                             *
c$$$*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
c$$$*----------------------------------------------------------------------*
c$$$*  5.12.96  M. SPIRA    					       *
c$$$* Es muss die Dittmaier-Routine D04 hinzugelinkt werden fuer die eta's *
c$$$*                   und die Spence-Funktionen                          *
c$$$************************************************************************
c$$$      IMPLICIT REAL*8 (A-H,O-Z)
c$$$      REAL*8 M1,M2,M3
c$$$      COMPLEX*16 C03,Cc
c$$$
c$$$      C03 = Cc(p1,p2,p3,m1**2,m2**2,m3**2)
c$$$
c$$$      return
c$$$      end
c$$$
C called as   C03(p02,p12,p22;m2,m0,m1)
C             Cc (p02,p12,p22;m0,m1,m2)


************************************************************************
        FUNCTION C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
* Es muss die Dittmaier-Routine D04 hinzugelinkt werden fuer die eta's *
*                   und die Spence-Funktionen                          *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
c ----------------------------------------------
      complex*16 cx,cy
c ----------------------------------------------
      COMPLEX*16 C03,CSPEN,ETA,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM
C     REAL*8 KAPPA
      COMPLEX*16 KAPPA
C     KAPPA(A,B,C) = SQRT(A**2+B**2+C**2-2*(A*B+A*C+B*C))
C     KAPPA(A,B,C) = SQRT(ABS(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      KAPPA(A,B,C) = SQRT(DCMPLX(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      EPS = 1.D-8
      IM = DCMPLX(0.D0,1.D0)
      IEPS = DCMPLX(0.D0,1.D-17)
      PI = 4*ATAN(1.D0)
      XX = 0.D0
C     IF(P1.LT.0.D0.OR.P2.LT.0.D0.OR.P3.LT.0.D0) XX=1.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
c ----------------------------------------------
      alpha = abs( kappa(q10,q21,q20) )
c      ALPHA = KAPPA(Q10,Q21,Q20)
c ----------------------------------------------
      ALP(0) = KAPPA(Q21,SM1,SM2)*(1+IEPS*Q21)
      ALP(1) = KAPPA(Q20,SM2,SM0)*(1+IEPS*Q20)
      ALP(2) = KAPPA(Q10,SM0,SM1)*(1+IEPS*Q10)
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
c ----------------------------------------------
      cdum = dcmplx(0.D0)
c      CDUM=0.D0
c ----------------------------------------------
      DO I=0,2
       DO J=1,2
        CDUM = CDUM + CSPEN((Y0(I)-1)/Y(I,J)) - CSPEN(Y0(I)/Y(I,J))
        CX = ETA(1-X(I,J),1/Y(I,J))
c ----------------------------------------------
        if (abs(cx).ne.0.D0) then 
c        IF(CX.NE.0.D0)THEN
c ----------------------------------------------
         CDUM = CDUM + CX*LOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = ETA(-X(I,J),1/Y(I,J))
c ----------------------------------------------
        if (abs(cy).ne.0.D0) then 
c        IF(CY.NE.0.D0)THEN 
c ----------------------------------------------
         CDUM = CDUM - CY*LOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = ETA(-X(I,1),-X(I,2))
c ----------------------------------------------
       if (abs(cx).ne.0.D0) then 
c       IF(CX.NE.0.D0)THEN
c ----------------------------------------------
        CDUM = CDUM - CX*LOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       CY = ETA(Y(I,1),Y(I,2))
c ----------------------------------------------
       if (abs(cy).ne.0.D0) then 
c       IF(CY.NE.0.D0)THEN
c ----------------------------------------------
        CDUM = CDUM + CY*LOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -AIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2*PI*IM*LOG((1-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      C03 = CDUM/ALPHA
      RETURN
      END

