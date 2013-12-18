C**********************************************************************
C***
C***  THIS ARE THE FUNCTIONS FOR THE SQUARK-GLUINO PROCESS
C***
C**********************************************************************

      REAL*8 FUNCTION DSGQG1(ALPHAS,S,TG,MS,MG,SCA)
C***  GIVES THE SCALE DEPENDENCE OF VIRTUAL
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,U1,TG,UG,MS,MG,SCA,PI,BETAL,PGG,PQQ,DSGQGB,NS

      PI = 4.D0*ATAN(1.D0)
      NS = 6.D0
      U1 = -S -TG
      UG = U1 +MS**2 -MG**2
      T1 = TG +MG**2 -MS**2
      BETAL = -11.D0 + 2.D0/3.D0*(NS -1.D0)
      PGG = 3.D0*( -2*LOG(-U1/MS**2) +11.D0/6.D0) -1.D0/3.D0*(NS-1.D0)
      PQQ = 4.D0/3.D0 * ( -2*LOG(-T1/MS**2) + 3.D0/2.D0)
      DSGQG1 = DSGQGB(ALPHAS,S,TG,MS,MG)*
     +     ALPHAS/2.D0/PI*LOG(1.D0/SCA)*( BETAL +PGG +PQQ)
      RETURN
      END


      REAL*8 FUNCTION DSGQG2(ALPHAS,S,TG,S4,MS,MG,DEL,S4MAX,SCA)
C***  GIVES THE SCALE DEPENDENCE OF LOG(DEL)
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,TG,S4,S4MAX,MS,MG,SCA,PI,DSGQGB,DLDEL1,DEL
      PI = 4.D0*ATAN(1.D0)
      DLDEL1 = LOG(S4MAX/MS**2) -(S4MAX -DEL)/S4
      DSGQG2 = 
     +     (2.D0*3.D0 +8.D0/3.D0)*DLDEL1*DSGQGB(ALPHAS,S,TG,MS,MG)
     +     *ALPHAS/2.D0/PI*LOG(1.D0/SCA)
      RETURN
      END


      REAL*8 FUNCTION DSGQG3(ALPHAS,S,TG,S4,MS,MG,SCA)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,TG,T1,U1,UG,S4,MS,MG,PI
      REAL*8 DSGQGB,PQQ1,PGG2,SCA,X1,X2
      PI = 4.D0*ATAN(1.D0)
      U1 = S4 -S -TG
      T1 = TG +MG**2 -MS**2
      UG = U1 -MG**2 +MS**2
      X1 = -T1/(S+UG)
      X2 = -U1/(S+TG)
      PQQ1 = 4.D0/3.D0*(1.D0 +X1**2)/(1.D0 -X1)
      PGG2 = 3.D0*(2.D0/(1.D0-X2) + 2.D0/X2 -4 +2*X2 -2*X2**2)
      DSGQG3 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) *
     +     (-1/T1*DSGQGB(ALPHAS,X1*S,TG,MS,MG)*PQQ1*X1**2
     +      -1/U1*DSGQGB(ALPHAS,X2*S,X2*TG,MS,MG)*PGG2*X2**2 )
      RETURN
      END



      REAL*8 FUNCTION DSGQGR(ALPHAS,S,TG,MS,MG)
C***  GIVES THE CHANGE FROM MSBAR TO DRBAR FOR THE YUKAWA COUPLINGS
      IMPLICIT NONE
      REAL*8 ALPHAS,S,TG,MS,MG,PI,DRMS2,DSGQGB

      PI = 4.D0*ATAN(1.D0)
      DRMS2 = ALPHAS *2.D0/3.D0/PI

      DSGQGR = DRMS2 * DSGQGB(ALPHAS,S,TG,MS,MG)
      RETURN
      END


      REAL*8 FUNCTION DSGQGB(ALPHAS,S,TG,MS,MG)
C***  BORN CROSS SECTIONS FOR Q + G -> SQ + GL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,U1,TG,MS,MG,MS2,MG2,M2
      REAL*8 M2QGB, NS,CONV,PI,N,CO,CK,AVG
      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0)/N
      U1 = -S -TG
      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2

      M2QGB = 0.D0
      M2QGB = M2QGB + CO * (  - 2 + 2*S**(-1)*TG*U1**(-2)*M2*MS2 - 
     +    S**(-1)*TG + 2*S**(-1)*U1**(-1)*M2*MG2 - 2*S**(-1)*M2 - 2*S*
     +    TG**(-1) - 4*TG**(-2)*M2*MG2 - 4*TG**(-1)*U1**(-1)*M2*MS2 - 4
     +    *TG**(-1)*M2 - 2*U1**(-1)*M2 )
     +
      M2QGB = M2QGB + CK * (  - 2*S**(-1)*TG*U1**(-2)*M2*MS2 + S**(-1)*
     +    TG - 2*S**(-1)*U1**(-1)*M2*MG2 + 2*S**(-1)*M2 + 2*U1**(-1)*M2
     +     )

      M2QGB = 2.D0* M2QGB
      DSGQGB = ALPHAS**2*AVG*M2QGB *PI/S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSGQGV(ALPHAS,S,TG,MS,MG,MT)
C***  CROSS SECTIONS FOR Q + G -> SQ + GL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T,T1,TG,U,U1,UG,SB,S1,MS,MG,MT,MS2,MG2,MT2,M2     
      REAL*8 NS,CONV,PI,ZETA2,N,CO,CK,CQED,AVG,BETA,XS,SPENCE,THREE
      REAL*8 XLAM,XSG,Q2MS2,FOUR
      REAL*8 SK1B0A(1:3), SK1B0B(1:5), SK1B0C(1:3), SK1B0D(1:4,1:2)
      REAL*8 SK1B0E(1:3), SK1BP(1:7), SK1C0A(1:8), SK1C0B(1:6)
      REAL*8 SK1C0C(1:8,1:2), SK1D0(1:12,1:2), SK2C0B(1:4)
      REAL*8 SK2C0C(1:6,1:2), SK2C0D(1:8,1:2), SK2D0(1:10,1:2)
      REAL*8 SL1B0A, SL1B0B, SL1B0C, SL1B0D, SL1B0E
      REAL*8 SL1BP, SL1C0A, SL1C0B, SL1C0C, SL1D0
      REAL*8 SL2C0B, SL2C0C, SL2C0D, SL2D0
      REAL*8 SOF2(1:8)
      REAL*8 M2QGV,DSGQGB, DSGQGR, DSGQG1
      COMPLEX*16 KBETAG,KXG,KBETAT,KXT            
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      CONV = 389379660.D0

      NS = 6
      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0) /N
      THREE = 3.D0
      FOUR = 4.D0

      MS2 = MS**2
      MG2 = MG**2
      MT2 = MT**2
      M2 = MG2 - MS2

      XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
      BETA = SQRT(1.D0 -4*MG*MS/(S-(MG-MS)**2))
      XSG = (1.D0 -BETA)/(1.D0 +BETA)

      U1 = -S -TG
      T = TG +MG2
      U = U1 +MS2
      T1 = TG +M2
      UG = U1 -M2
      S1 = S -MS2 -MG2
      
C     BETA = SQRT(1 -4*MS2/S)
      BETA = REAL(SQRT( DCMPLX(1.D0 - 4.D0*MS2/S)))
      SB = S * BETA
      XS = (1.D0 -BETA)/(1.D0 +BETA)
      KBETAG = SQRT( DCMPLX(1.D0 - 4.D0*MG2/S)) 
      KXG = (1.D0 -KBETAG)/(1.D0+KBETAG)
      KBETAT = SQRT( DCMPLX(1.D0 - 4.D0*MT2/S)) 
      KXT = (1.D0 -KBETAT)/(1.D0+KBETAT)


      IF (MS.EQ.MG) THEN

C$$$      SK1B0A(1) = SL1B0A(MS2,MG2,MT2,1)
      SK1B0A(2) = SL1B0A(MS2,MG2,MT2,2)
C$$$      SK1B0A(3) = SL1B0A(MS2,MG2,MT2,3)

      SK1B0B(1) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,1)
      SK1B0B(2) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,2)
C$$$  SK1B0B(3) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,3)
C$$$  SK1B0B(4) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,4)
C$$$  SK1B0B(5) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,5)

      SK1B0C(1) = SL1B0C(MS2,MG2,MT2,1)
C$$$      SK1B0C(2) = SL1B0C(MS2,MG2,MT2,2)
C$$$      SK1B0C(3) = SL1B0C(MS2,MG2,MT2,3)

      SK1B0D(1,1) = SL1B0D(MS2,MG2,MT2,T,1)
C$$$      SK1B0D(2,1) = SL1B0D(MS2,MG2,MT2,T,2)
      SK1B0D(3,1) = SL1B0D(MS2,MG2,MT2,T,3)
C$$$      SK1B0D(4,1) = SL1B0D(MS2,MG2,MT2,T,4)
      SK1B0D(1,2) = SL1B0D(MS2,MG2,MT2,U,1)
C$$$      SK1B0D(2,2) = SL1B0D(MS2,MG2,MT2,U,2)
C$$$      SK1B0D(3,2) = SL1B0D(MS2,MG2,MT2,U,3)
C$$$      SK1B0D(4,2) = SL1B0D(MS2,MG2,MT2,U,4)

C$$$      SK1B0E(1) = SL1B0E(MS2,MG2,MT2,1)
C$$$      SK1B0E(2) = SL1B0E(MS2,MG2,MT2,2)
      SK1B0E(3) = SL1B0E(MS2,MG2,MT2,3)

      SK1BP(1) = SL1BP(MS2,MG2,MT2,1)
C$$$      SK1BP(2) = SL1BP(MS2,MG2,MT2,2)
C$$$      SK1BP(3) = SL1BP(MS2,MG2,MT2,3)
C$$$      SK1BP(4) = SL1BP(MS2,MG2,MT2,4)
      SK1BP(5) = SL1BP(MG2,MS2,MT2,5)
C$$$      SK1BP(6) = SL1BP(MS2,MG2,MT2,6)
C$$$      SK1BP(7) = SL1BP(MS2,MG2,MT2,7)


      SK1C0A(1) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,1)
      SK1C0A(2) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,2)
C$$$      SK1C0A(3) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,3)
C$$$      SK1C0A(4) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,4)
C$$$      SK1C0A(5) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,5)
C$$$      SK1C0A(6) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,6)
C$$$      SK1C0A(7) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,7)
C$$$      SK1C0A(8) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,8)

      SK1C0B(1) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,1)
C$$$      SK1C0B(2) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,2)
      SK1C0B(3) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,3)
C$$$      SK1C0B(4) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,4)
C$$$      SK1C0B(5) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,5)
C$$$      SK1C0B(6) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,6)

      SK1C0C(1,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$      SK1C0C(2,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,2)
      SK1C0C(3,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,3)
C$$$      SK1C0C(4,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,4)
C$$$      SK1C0C(5,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,5)
C$$$      SK1C0C(6,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,6)
C$$$      SK1C0C(7,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,7)
C$$$      SK1C0C(8,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,8)

      SK1C0C(1,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,1)
C$$$      SK1C0C(2,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,2)
      SK1C0C(3,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,3)
C$$$      SK1C0C(4,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,4)
C$$$      SK1C0C(5,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,5)
C$$$      SK1C0C(6,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,6)
C$$$      SK1C0C(7,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,7)
C$$$      SK1C0C(8,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,8)


      SK1D0(1,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,1)
      SK1D0(2,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,2)
      SK1D0(3,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,3)
C$$$      SK1D0(4,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,4)
C$$$      SK1D0(5,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,5)
C$$$      SK1D0(6,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,6)
C$$$      SK1D0(7,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,7)
C$$$      SK1D0(8,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,8)
C$$$      SK1D0(9,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,9)
C$$$      SK1D0(10,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,10)
C$$$      SK1D0(11,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,11)
C$$$      SK1D0(12,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,12)

      SK1D0(1,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,1)
      SK1D0(2,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,2)
      SK1D0(3,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,3)
C$$$      SK1D0(4,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,4)
C$$$      SK1D0(5,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,5)
C$$$      SK1D0(6,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,6)
C$$$      SK1D0(7,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,7)
C$$$      SK1D0(8,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,8)
C$$$      SK1D0(9,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,9)
C$$$      SK1D0(10,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,10)
C$$$      SK1D0(11,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,11)
C$$$      SK1D0(12,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,12)


C$$$      SK2C0B(1) = SL2C0B(MS2,MG2,MT2,S,ZETA2,1)
C$$$      SK2C0B(2) = SL2C0B(MS2,MG2,MT2,S,ZETA2,2)
C$$$      SK2C0B(3) = SL2C0B(MS2,MG2,MT2,S,ZETA2,3)
C$$$      SK2C0B(4) = SL2C0B(MS2,MG2,MT2,S,ZETA2,4)

C$$$      SK2C0C(1,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$      SK2C0C(2,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,2)
C$$$      SK2C0C(3,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,3)
C$$$      SK2C0C(4,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,4)
C$$$      SK2C0C(5,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,5)
C$$$      SK2C0C(6,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,6)

C$$$      SK2C0C(1,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,1)
C$$$      SK2C0C(2,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,2)
C$$$      SK2C0C(3,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,3)
C$$$      SK2C0C(4,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,4)
C$$$      SK2C0C(5,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,5)
C$$$      SK2C0C(6,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,6)


C$$$      SK2C0D(1,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,1)
C$$$      SK2C0D(2,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,2)
C$$$      SK2C0D(3,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,3)
C$$$      SK2C0D(4,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,4)
C$$$      SK2C0D(5,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,5)
C$$$      SK2C0D(6,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,6)
      SK2C0D(7,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,7)
      SK2C0D(8,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,8)

C$$$      SK2C0D(1,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,1)
C$$$      SK2C0D(2,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,2)
C$$$      SK2C0D(3,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,3)
C$$$      SK2C0D(4,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,4)
C$$$      SK2C0D(5,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,5)
C$$$      SK2C0D(6,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,6)
C$$$      SK2C0D(7,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,7)
C$$$      SK2C0D(8,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,8)

C$$$      SK2D0(1,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,1)
C$$$      SK2D0(2,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,2)
C$$$      SK2D0(3,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,3)
C$$$      SK2D0(4,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,4)
C$$$      SK2D0(5,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,5)
C$$$      SK2D0(6,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,6)
C$$$      SK2D0(7,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,7)
C$$$      SK2D0(8,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,8)
C$$$      SK2D0(9,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,9)
C$$$      SK2D0(10,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,10)

C$$$      SK2D0(1,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,1)
C$$$      SK2D0(2,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,2)
C$$$      SK2D0(3,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,3)
C$$$      SK2D0(4,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,4)
C$$$      SK2D0(5,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,5)
C$$$      SK2D0(6,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,6)
C$$$      SK2D0(7,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,7)
C$$$      SK2D0(8,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,8)
C$$$      SK2D0(9,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,9)
C$$$      SK2D0(10,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,10)

      ELSE

      SK1B0A(1) = SL1B0A(MS2,MG2,MT2,1)
      SK1B0A(2) = SL1B0A(MS2,MG2,MT2,2)
      SK1B0A(3) = SL1B0A(MS2,MG2,MT2,3)

      SK1B0B(1) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,1)
C$$$  SK1B0B(2) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,2)
C$$$  SK1B0B(3) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,3)
C$$$  SK1B0B(4) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,4)
      SK1B0B(5) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,5)

      SK1B0C(1) = SL1B0C(MS2,MG2,MT2,1)
      SK1B0C(2) = SL1B0C(MS2,MG2,MT2,2)
C$$$      SK1B0C(3) = SL1B0C(MS2,MG2,MT2,3)

      SK1B0D(1,1) = SL1B0D(MS2,MG2,MT2,T,1)
      SK1B0D(2,1) = SL1B0D(MS2,MG2,MT2,T,2)
      SK1B0D(3,1) = SL1B0D(MS2,MG2,MT2,T,3)
C$$$      SK1B0D(4,1) = SL1B0D(MS2,MG2,MT2,T,4)
      SK1B0D(1,2) = SL1B0D(MS2,MG2,MT2,U,1)
      SK1B0D(2,2) = SL1B0D(MS2,MG2,MT2,U,2)
C$$$      SK1B0D(3,2) = SL1B0D(MS2,MG2,MT2,U,3)
C$$$      SK1B0D(4,2) = SL1B0D(MS2,MG2,MT2,U,4)

      SK1B0E(1) = SL1B0E(MS2,MG2,MT2,1)
      SK1B0E(2) = SL1B0E(MS2,MG2,MT2,2)
      SK1B0E(3) = SL1B0E(MS2,MG2,MT2,3)

      SK1BP(1) = SL1BP(MS2,MG2,MT2,1)
      SK1BP(2) = SL1BP(MS2,MG2,MT2,2)
      SK1BP(3) = SL1BP(MS2,MG2,MT2,3)
C$$$      SK1BP(4) = SL1BP(MS2,MG2,MT2,4)
      SK1BP(5) = SL1BP(MG2,MS2,MT2,5)
      SK1BP(6) = SL1BP(MS2,MG2,MT2,6)
      SK1BP(7) = SL1BP(MS2,MG2,MT2,7)


      SK1C0A(1) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,1)
C$$$      SK1C0A(2) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,2)
C$$$      SK1C0A(3) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,3)
C$$$      SK1C0A(4) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,4)
C$$$      SK1C0A(5) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,5)
C$$$      SK1C0A(6) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,6)
      SK1C0A(7) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,7)
      SK1C0A(8) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,8)

C$$$      SK1C0B(1) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,1)
C$$$      SK1C0B(2) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,2)
C$$$      SK1C0B(3) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,3)
C$$$      SK1C0B(4) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,4)
C$$$      SK1C0B(5) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,5)
C$$$      SK1C0B(6) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,6)

C$$$      SK1C0C(1,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$      SK1C0C(2,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,2)
C$$$      SK1C0C(3,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,3)
C$$$      SK1C0C(4,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,4)
C$$$      SK1C0C(5,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,5)
C$$$      SK1C0C(6,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,6)
C$$$      SK1C0C(7,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,7)
C$$$      SK1C0C(8,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,8)

C$$$      SK1C0C(1,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,1)
C$$$      SK1C0C(2,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,2)
C$$$      SK1C0C(3,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,3)
C$$$      SK1C0C(4,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,4)
C$$$      SK1C0C(5,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,5)
C$$$      SK1C0C(6,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,6)
C$$$      SK1C0C(7,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,7)
C$$$      SK1C0C(8,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,8)


C$$$      SK1D0(1,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,1)
C$$$      SK1D0(2,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,2)
C$$$      SK1D0(3,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,3)
C$$$      SK1D0(4,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,4)
C$$$      SK1D0(5,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,5)
C$$$      SK1D0(6,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,6)
C$$$      SK1D0(7,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,7)
C$$$      SK1D0(8,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,8)
C$$$      SK1D0(9,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,9)
C$$$      SK1D0(10,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,10)
C$$$      SK1D0(11,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,11)
C$$$      SK1D0(12,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,12)

C$$$      SK1D0(1,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,1)
C$$$      SK1D0(2,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,2)
C$$$      SK1D0(3,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,3)
C$$$      SK1D0(4,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,4)
C$$$      SK1D0(5,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,5)
C$$$      SK1D0(6,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,6)
C$$$      SK1D0(7,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,7)
C$$$      SK1D0(8,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,8)
C$$$      SK1D0(9,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,9)
C$$$      SK1D0(10,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,10)
C$$$      SK1D0(11,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,11)
C$$$      SK1D0(12,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,12)


      SK2C0B(1) = SL2C0B(MS2,MG2,MT2,S,ZETA2,1)
      SK2C0B(2) = SL2C0B(MS2,MG2,MT2,S,ZETA2,2)
      SK2C0B(3) = SL2C0B(MS2,MG2,MT2,S,ZETA2,3)
      SK2C0B(4) = SL2C0B(MS2,MG2,MT2,S,ZETA2,4)

C$$$      SK2C0C(1,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$      SK2C0C(2,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,2)
      SK2C0C(3,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,3)
      SK2C0C(4,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,4)
      SK2C0C(5,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,5)
      SK2C0C(6,1) = SL2C0C(MS2,MG2,MT2,T,ZETA2,6)

      SK2C0C(1,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,1)
      SK2C0C(2,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,2)
      SK2C0C(3,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,3)
      SK2C0C(4,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,4)
C$$$      SK2C0C(5,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,5)
C$$$      SK2C0C(6,2) = SL2C0C(MS2,MG2,MT2,U,ZETA2,6)


      SK2C0D(1,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,1)
      SK2C0D(2,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,2)
      SK2C0D(3,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,3)
      SK2C0D(4,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,4)
C$$$      SK2C0D(5,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,5)
C$$$      SK2C0D(6,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,6)
      SK2C0D(7,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,7)
      SK2C0D(8,1) = SL2C0D(MS2,MG2,MT2,T,ZETA2,8)

C$$$      SK2C0D(1,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,1)
C$$$      SK2C0D(2,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,2)
      SK2C0D(3,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,3)
      SK2C0D(4,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,4)
      SK2C0D(5,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,5)
      SK2C0D(6,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,6)
C$$$      SK2C0D(7,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,7)
C$$$      SK2C0D(8,2) = SL2C0D(MS2,MG2,MT2,U,ZETA2,8)

      SK2D0(1,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,1)
      SK2D0(2,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,2)
C$$$      SK2D0(3,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,3)
      SK2D0(4,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,4)
C$$$      SK2D0(5,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,5)
      SK2D0(6,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,6)
      SK2D0(7,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,7)
      SK2D0(8,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,8)
      SK2D0(9,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,9)
      SK2D0(10,1) = SL2D0(MS2,MG2,MT2,S,T,U,ZETA2,10)

      SK2D0(1,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,1)
      SK2D0(2,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,2)
      SK2D0(3,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,3)
C$$$      SK2D0(4,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,4)
      SK2D0(5,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,5)
C$$$      SK2D0(6,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,6)
C$$$      SK2D0(7,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,7)
C$$$      SK2D0(8,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,8)
C$$$      SK2D0(9,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,9)
C$$$      SK2D0(10,2) = SL2D0(MS2,MG2,MT2,S,U,T,ZETA2,10)

      END IF

      SOF2(1) = 0.125D0/XLAM * ( +4*SPENCE(XSG) +4*SPENCE(-XSG) 
     +     +4*LOG(XSG)*LOG(1-XSG**2)-2*LOG(XSG)**2 -2*ZETA2 )
      SOF2(2) = 0.125D0/MG2* (+2*S1/XLAM*LOG(XSG))
      SOF2(3) = 0.125D0/MS2*(-2)
      SOF2(4) = 0.125D0/U1*(-1.5D0*ZETA2)
      SOF2(5) = 0.125D0/T1*(-1.5D0*ZETA2)
      SOF2(6) = 0.125D0/TG*(-2*LOG(XSG)**2 + LOG(XSG*MG*U1/MS/TG)**2
     +     +2*SPENCE(1 -MS*TG/MG/U1/XSG) -2*SPENCE(1-MG*U1/MS/TG/XSG)
     +     -1.5D0*ZETA2 )
      SOF2(7) = 0.125D0/UG*(-2*LOG(XSG)**2 + LOG(XSG*MG*T1/MS/UG)**2
     +     +2*SPENCE(1 -MS*UG/MG/T1/XSG) -2*SPENCE(1-MG*T1/MS/UG/XSG)
     +     -1.5D0*ZETA2 )
      SOF2(8) = 0.125D0/S*( -2*SPENCE(1 -T1*U1/MS2/S) -3*ZETA2 )


      IF (MS.EQ.MG) THEN

      M2QGV = 0.D0
      M2QGV = M2QGV + N*CO * ( 4 - 2*S**(-1)*T1*T**(-1)*MS2 + 2*S**(-1)
     +    *T1 - 2*S**(-1)*T**(-1)*MS2**2 + 2*S**(-1)*MS2 - 8*S*T1**(-3)
     +    *T**(-1)*MS2**3 + 8*S*T1**(-3)*MS2**2 - 24*S*T1**(-2)*T**(-1)
     +    *MS2**2 + 20*S*T1**(-2)*MS2 - 16*S*T1**(-1)*T**(-1)*MS2 + 4*S
     +    *T1**(-1) - 4*T1**(-2)*T**(-1)*MS2**3 + 4*T1**(-2)*MS2**2 - 
     +    20*T1**(-1)*T**(-1)*MS2**2 + 14*T1**(-1)*MS2 - 2*U1**(-1)*
     +    T**(-1)*MS2**2 - 12*T**(-1)*MS2 )
     +
      M2QGV = M2QGV + N*CK * (  - 2*S**(-1)*T1 )
     +
      M2QGV = M2QGV + CO * ( 2*S**(-1)*T1*NS*T**(-1)*MS2 - S**(-1)*T1*
     +    NS - 2*S**(-1)*T1*T**(-1)*MT2 + S**(-1)*T1*MS2**(-1)*MT2 + 2*
     +    S**(-1)*NS*T**(-1)*MS2**2 - 2*S**(-1)*NS*MS2 - 2*S**(-1)*
     +    T**(-1)*MS2*MT2 + 2*S**(-1)*MT2 + 8*S*T1**(-3)*NS*T**(-1)*
     +    MS2**3 - 8*S*T1**(-3)*NS*MS2**2 - 8*S*T1**(-3)*T**(-1)*MS2**2
     +    *MT2 + 8*S*T1**(-3)*MS2*MT2 + 24*S*T1**(-2)*NS*T**(-1)*MS2**2
     +     - 20*S*T1**(-2)*NS*MS2 - 24*S*T1**(-2)*T**(-1)*MS2*MT2 + 20*
     +    S*T1**(-2)*MT2 + 16*S*T1**(-1)*NS*T**(-1)*MS2 - 2*S*T1**(-1)*
     +    NS - 16*S*T1**(-1)*T**(-1)*MT2 + 2*S*T1**(-1)*MS2**(-1)*MT2
     +     + 4*T1**(-2)*NS*T**(-1)*MS2**3 - 4*T1**(-2)*NS*MS2**2 - 4*
     +    T1**(-2)*T**(-1)*MS2**2*MT2 + 4*T1**(-2)*MS2*MT2 + 20*
     +    T1**(-1)*NS*T**(-1)*MS2**2 - 14*T1**(-1)*NS*MS2 - 20*T1**(-1)
     +    *T**(-1)*MS2*MT2 + 14*T1**(-1)*MT2 + 2*U1**(-1)*NS*T**(-1)*
     +    MS2**2 - 2*U1**(-1)*T**(-1)*MS2*MT2 + 12*NS*T**(-1)*MS2 - 2*
     +    NS - 12*T**(-1)*MT2 + 2*MS2**(-1)*MT2 )
     +
      M2QGV = M2QGV + CK * ( S**(-1)*T1*NS - S**(-1)*T1*MS2**(-1)*MT2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CO*THREE**(-1) * ( 4 - 6*S**(-1)*T1*
     +    T**(-1)*MT2 + 3*S**(-1)*T1*MS2**(-1)*MT2 + 2*S**(-1)*T1 - 6*
     +    S**(-1)*T**(-1)*MS2*MT2 + 6*S**(-1)*MT2 - 24*S*T1**(-3)*
     +    T**(-1)*MS2**2*MT2 + 24*S*T1**(-3)*MS2*MT2 - 72*S*T1**(-2)*
     +    T**(-1)*MS2*MT2 + 60*S*T1**(-2)*MT2 - 48*S*T1**(-1)*T**(-1)*
     +    MT2 + 6*S*T1**(-1)*MS2**(-1)*MT2 + 4*S*T1**(-1) - 12*T1**(-2)
     +    *T**(-1)*MS2**2*MT2 + 12*T1**(-2)*MS2*MT2 - 60*T1**(-1)*
     +    T**(-1)*MS2*MT2 + 42*T1**(-1)*MT2 - 6*U1**(-1)*T**(-1)*MS2*
     +    MT2 - 36*T**(-1)*MT2 + 6*MS2**(-1)*MT2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CK*THREE**(-1) * (  - 3*S**(-1)*T1*
     +    MS2**(-1)*MT2 - 2*S**(-1)*T1 )
     +
      M2QGV = M2QGV + SK1B0B(1)*N*CO * ( S**(-1)*T1 - 4*T1**(-1)*MS2**2
     +    *(S-4*MS2)**(-1) - T1*(S-4*MS2)**(-1) - 4*MS2*(S-4*MS2)**(-1)
     +     )
     +
      M2QGV = M2QGV + SK1B0B(1)*N*CK * (  - 2 - 3*S**(-1)*T1 - 2*
     +    T1**(-1)*MS2 - 4*T1**(-1)*MS2**2*(S-4*MS2)**(-1) - T1*
     +    (S-4*MS2)**(-1) + 2*U1**(-1)*MS2 - 4*MS2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1B0B(1)*CQED * ( S**(-1)*T1 + T1*
     +    (S-4*MS2)**(-1) - 2*U1**(-1)*MS2 - 4*U1**(-1)*MS2**2*
     +    (S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1B0B(2)*N*CO * ( 2 + S**(-1)*T1 + 2*T1**(-1)*
     +    MS2 + 4*T1**(-1)*MS2**2*(S-4*MS2)**(-1) + T1*(S-4*MS2)**(-1)
     +     + 4*MS2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1B0B(2)*N*CK * (  - 3*S**(-1)*T1 + 4*T1**(-1)*
     +    MS2**2*(S-4*MS2)**(-1) + T1*(S-4*MS2)**(-1) + 2*U1**(-1)*MS2
     +     + 4*MS2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1B0B(2)*CQED * ( S**(-1)*T1 - T1*
     +    (S-4*MS2)**(-1) + 4*U1**(-1)*MS2**2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CO * (  - 10 - 5*S**(-1)*T1 + 4*S*
     +    T1**(-3)*MS2**2 - 8*S*T1**(-2)*MS2 - 10*S*T1**(-1) + 6*
     +    T1**(-2)*MS2**2 - 2*T1**(-1)*U1**(-1)*MS2**2 - 6*T1**(-1)*MS2
     +     - 4*U1**(-2)*MS2**2 - 2*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CK * ( 2 + 7*S**(-1)*T1 + 4*S*
     +    T1**(-2)*MS2 + 2*S*T1**(-1) - 8*T1**(-2)*MS2**2 + 4*T1**(-1)*
     +    MS2 + 12*U1**(-2)*MS2**2 + 6*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*CO * (  - 4 + 2*S**(-1)*T1*NS - 2*
     +    S**(-1)*T1 - 4*S*T1**(-3)*NS*MS2**2 + 4*S*T1**(-3)*MS2**2 + 4
     +    *S*T1**(-2)*NS*MS2 - 4*S*T1**(-2)*MS2 + 4*S*T1**(-1)*NS - 4*S
     +    *T1**(-1) + 2*T1**(-2)*NS*MS2**2 - 2*T1**(-2)*MS2**2 + 2*
     +    T1**(-1)*U1**(-1)*NS*MS2**2 - 2*T1**(-1)*U1**(-1)*MS2**2 + 2*
     +    T1**(-1)*NS*MS2 - 2*T1**(-1)*MS2 + 4*NS )
     +
      M2QGV = M2QGV + SK1B0C(1)*CK * (  - 2*S**(-1)*T1*NS + 2*S**(-1)*
     +    T1 )
     +
      M2QGV = M2QGV + SK1B0C(1)*CQED * (  - S**(-1)*T1 - 4*U1**(-2)*
     +    MS2**2 - 2*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*N*CO * ( 2 + 2*S**(-1)*T1*T**(-1)*MS2
     +     + 2*S**(-1)*T**(-1)*MS2**2 - 2*S**(-1)*MS2 + 8*S*T1**(-3)*
     +    T**(-1)*MS2**3 - 12*S*T1**(-3)*MS2**2 + 24*S*T1**(-2)*T**(-1)
     +    *MS2**2 - 12*S*T1**(-2)*MS2 + 16*S*T1**(-1)*T**(-1)*MS2 + 4*S
     +    *T1**(-1) + 4*T1**(-2)*T**(-1)*MS2**3 - 10*T1**(-2)*MS2**2 - 
     +    2*T1**(-1)*U1**(-1)*MS2**2 + 20*T1**(-1)*T**(-1)*MS2**2 - 14*
     +    T1**(-1)*MS2 + 2*U1**(-1)*T**(-1)*MS2**2 - 2*U1**(-1)*MS2 + 
     +    12*T**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*N*CK * (  - 2 - 4*S*T1**(-2)*MS2 - 4*
     +    S*T1**(-1) + 8*T1**(-2)*MS2**2 + 4*T1**(-1)*U1**(-1)*MS2**2
     +     + 2*T1**(-1)*MS2 + 2*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO * (  - 2*S**(-1)*T1*NS*T**(-1)*MS2
     +     + 2*S**(-1)*T1*T**(-1)*MS2 - 2*S**(-1)*NS*T**(-1)*MS2**2 + 2
     +    *S**(-1)*NS*MS2 + 2*S**(-1)*T**(-1)*MS2**2 - 2*S**(-1)*MS2 - 
     +    8*S*T1**(-3)*NS*T**(-1)*MS2**3 + 12*S*T1**(-3)*NS*MS2**2 + 8*
     +    S*T1**(-3)*T**(-1)*MS2**3 - 12*S*T1**(-3)*MS2**2 - 24*S*
     +    T1**(-2)*NS*T**(-1)*MS2**2 + 16*S*T1**(-2)*NS*MS2 + 24*S*
     +    T1**(-2)*T**(-1)*MS2**2 - 16*S*T1**(-2)*MS2 - 16*S*T1**(-1)*
     +    NS*T**(-1)*MS2 + 16*S*T1**(-1)*T**(-1)*MS2 - 4*T1**(-2)*NS*
     +    T**(-1)*MS2**3 + 2*T1**(-2)*NS*MS2**2 + 4*T1**(-2)*T**(-1)*
     +    MS2**3 - 2*T1**(-2)*MS2**2 - 2*T1**(-1)*U1**(-1)*NS*MS2**2 + 
     +    2*T1**(-1)*U1**(-1)*MS2**2 - 20*T1**(-1)*NS*T**(-1)*MS2**2 + 
     +    12*T1**(-1)*NS*MS2 + 20*T1**(-1)*T**(-1)*MS2**2 - 12*T1**(-1)
     +    *MS2 - 2*U1**(-1)*NS*T**(-1)*MS2**2 + 2*U1**(-1)*T**(-1)*
     +    MS2**2 - 12*NS*T**(-1)*MS2 + 12*T**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CO * ( 4*T1**(-1)*U1**(-1)*MS2**2
     +     + 4*T1**(-1)*MS2 + 4*U1**(-2)*MS2**2 + 4*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CK * (  - 4*T1**(-1)*U1**(-1)*
     +    MS2**2 - 4*T1**(-1)*MS2 - 12*U1**(-2)*MS2**2 - 12*U1**(-1)*
     +    MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*CQED * ( 4*U1**(-2)*MS2**2 + 4*
     +    U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO * (  - 2*S**(-1)*T1*T**(-1)*MS2 + 
     +    2*S**(-1)*T1*T**(-1)*MT2 + 2*S**(-1)*T**(-1)*MS2*MT2 - 2*
     +    S**(-1)*T**(-1)*MS2**2 + 2*S**(-1)*MS2 - 2*S**(-1)*MT2 + 8*S*
     +    T1**(-3)*T**(-1)*MS2**2*MT2 - 8*S*T1**(-3)*T**(-1)*MS2**3 - 
     +    12*S*T1**(-3)*MS2*MT2 + 12*S*T1**(-3)*MS2**2 + 24*S*T1**(-2)*
     +    T**(-1)*MS2*MT2 - 24*S*T1**(-2)*T**(-1)*MS2**2 + 16*S*
     +    T1**(-2)*MS2 - 16*S*T1**(-2)*MT2 - 16*S*T1**(-1)*T**(-1)*MS2
     +     + 16*S*T1**(-1)*T**(-1)*MT2 + 4*T1**(-2)*T**(-1)*MS2**2*MT2
     +     - 4*T1**(-2)*T**(-1)*MS2**3 - 2*T1**(-2)*MS2*MT2 + 2*
     +    T1**(-2)*MS2**2 + 2*T1**(-1)*U1**(-1)*MS2*MT2 - 2*T1**(-1)*
     +    U1**(-1)*MS2**2 + 20*T1**(-1)*T**(-1)*MS2*MT2 - 20*T1**(-1)*
     +    T**(-1)*MS2**2 + 12*T1**(-1)*MS2 - 12*T1**(-1)*MT2 + 2*
     +    U1**(-1)*T**(-1)*MS2*MT2 - 2*U1**(-1)*T**(-1)*MS2**2 - 12*
     +    T**(-1)*MS2 + 12*T**(-1)*MT2 )
     +
      M2QGV = M2QGV + SK1B0E(3)*CO * ( 4 - S**(-1)*T1*MS2**(-1)*MT2 + 2
     +    *S**(-1)*T1 + 4*S*T1**(-3)*MS2*MT2 - 4*S*T1**(-3)*MS2**2 + 4*
     +    S*T1**(-2)*MS2 - 4*S*T1**(-2)*MT2 - 2*S*T1**(-1)*MS2**(-1)*
     +    MT2 + 4*S*T1**(-1) - 2*T1**(-2)*MS2*MT2 + 2*T1**(-2)*MS2**2
     +     - 2*T1**(-1)*U1**(-1)*MS2*MT2 + 2*T1**(-1)*U1**(-1)*MS2**2
     +     + 2*T1**(-1)*MS2 - 2*T1**(-1)*MT2 - 2*MS2**(-1)*MT2 )
     +
      M2QGV = M2QGV + SK1B0E(3)*CK * ( S**(-1)*T1*MS2**(-1)*MT2 - 2*
     +    S**(-1)*T1 )
     +
      M2QGV = M2QGV + SK1BP(1)*N*CO * (  - 6*S**(-1)*T1*MS2 - 12*S*
     +    T1**(-1)*MS2 - 12*MS2 )
     +
      M2QGV = M2QGV + SK1BP(1)*N*CK * ( 10*S**(-1)*T1*MS2 + 4*S*
     +    T1**(-1)*MS2 + 4*MS2 )
     +
      M2QGV = M2QGV + SK1BP(1)*CQED * (  - 2*S**(-1)*T1*MS2 )
     +
      M2QGV = M2QGV + SK1BP(5)*CO * ( 2*S**(-1)*T1*MT2 + 4*S*T1**(-1)*
     +    MT2 + 4*MT2 )
     +
      M2QGV = M2QGV + SK1BP(5)*CK * (  - 2*S**(-1)*T1*MT2 )
     +
      M2QGV = M2QGV + SK1C0A(1)*N*CO * ( 2*S + S**2*T1**(-1) + T1 )
     +
      M2QGV = M2QGV + SK1C0A(1)*N*CK * (  - S*T1**(-1)*MS2 + S*U1**(-1)
     +    *MS2 + S + T1 )
     +
      M2QGV = M2QGV + SK1C0A(1)*CQED * (  - S*U1**(-1)*MS2 + S + S**2*
     +    U1**(-1) - T1 )
     +
      M2QGV = M2QGV + SK1C0A(2)*N*CO * ( 2*S*T1**(-1)*MS2 - S**2*
     +    T1**(-1) - 2*T1**(-1)*MS2**2 + MS2 )
     +
      M2QGV = M2QGV + SK1C0A(2)*N*CK * (  - S*T1**(-1)*MS2 - S*U1**(-1)
     +    *MS2 - S - 2*T1**(-1)*MS2**2 + MS2 )
     +
      M2QGV = M2QGV + SK1C0A(2)*CQED * (  - 2*S*U1**(-1)*MS2 - 2*
     +    U1**(-1)*MS2**2 - MS2 )
     +
      M2QGV = M2QGV + SK2C0D(7,1)*CO * (  - 4*S*T1**(-2)*MS2*MT2 + 2*
     +    T1**(-1)*MS2*MT2 + 2*U1**(-1)*MS2*MT2 )
     +
      M2QGV = M2QGV + SK2C0D(8,1)*CO * ( 4*S*T1**(-2)*MS2**2 - 2*
     +    T1**(-1)*MS2**2 - 2*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SK1C0B(1)*N*CO * (  - 4*S**(-1)*T1*MS2 - 2*S*
     +    T1**(-1)*MS2 + 4*S + S**2*T1**(-1) + 2*T1 - 9*MS2 )
     +
      M2QGV = M2QGV + SK1C0B(1)*N*CK * ( 4*S**(-1)*T1*MS2 - S*T1**(-1)*
     +    MS2 + S*U1**(-1)*MS2 - S + 4*T1**(-1)*MS2**2 - 2*T1 - 4*
     +    U1**(-1)*MS2**2 + 3*MS2 )
     +
      M2QGV = M2QGV + SK1C0B(1)*CQED * (  - 2*S*U1**(-1)*MS2 + 4*
     +    U1**(-1)*MS2**2 - MS2 )
     +
      M2QGV = M2QGV + SK1C0B(3)*N*CO * (  - 2*S - S**2*T1**(-1) - 2*
     +    T1**(-1)*MS2**2 - 8*T1**(-1)*MS2**3*(S-4*MS2)**(-1) - 2*T1*
     +    MS2*(S-4*MS2)**(-1) - T1 - 2*MS2 - 8*MS2**2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1C0B(3)*N*CK * (  - S*T1**(-1)*MS2 - S*U1**(-1)
     +    *MS2 + S - 2*T1**(-1)*MS2**2 - 8*T1**(-1)*MS2**3*
     +    (S-4*MS2)**(-1) - 2*T1*MS2*(S-4*MS2)**(-1) + 3*T1 - 2*MS2 - 8
     +    *MS2**2*(S-4*MS2)**(-1) )
     +
      M2QGV = M2QGV + SK1C0B(3)*CQED * (  - 3*S*U1**(-1)*MS2 + S + S**2
     +    *U1**(-1) + 2*T1*MS2*(S-4*MS2)**(-1) - T1 - 2*U1**(-1)*MS2**2
     +     - 8*U1**(-1)*MS2**3*(S-4*MS2)**(-1) - 2*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,1)*N*CO * (  - 4*S*T1**(-2)*MS2**2 - S*
     +    U1**(-1)*MS2 - 4*S - 2*T1**(-1)*MS2**2 - T1 + 3*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,1)*N*CK * ( S*U1**(-1)*MS2 - 4*T1**(-1)*
     +    MS2**2 - T1 - 2*U1**(-1)*MS2**2 - 7*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,1)*CO * ( 4*S*T1**(-2)*NS*MS2**2 - 4*S*
     +    T1**(-2)*MS2**2 - 2*T1**(-1)*NS*MS2**2 + 2*T1**(-1)*MS2**2 - 
     +    2*U1**(-1)*NS*MS2**2 + 2*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SK1C0C(1,1)*CQED * (  - 3*S**(-1)*T1*MS2 - 
     +    S**(-1)*T1**2 + 3*S*U1**(-1)*MS2 + 3*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,2)*N*CO * (  - 6*S*T1**(-1)*MS2 + 3*S + 
     +    2*S**2*T1**(-1) + 2*T1**(-1)*MS2**2 + T1 + 2*U1**(-1)*MS2**2
     +     - MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,2)*N*CK * ( 3*S + 2*T1**(-1)*MS2**2 + 3*
     +    T1 + 2*U1**(-1)*MS2**2 + MS2 )
     +
      M2QGV = M2QGV + SK1C0C(1,2)*CQED * ( 3*S**(-1)*T1*MS2 + S**(-1)*
     +    T1**2 + T1 - 2*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SK1C0C(3,1)*N*CO * ( S*U1**(-1)*MS2 - T1 + MS2 )
     +
      M2QGV = M2QGV + SK1C0C(3,1)*N*CK * (  - S*U1**(-1)*MS2 - 4*S - 5*
     +    T1 - 5*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(3,1)*CQED * (  - 3*S**(-1)*T1*MS2 - 
     +    S**(-1)*T1**2 + 3*S*U1**(-1)*MS2 - 2*S - 2*S**2*U1**(-1) + 2*
     +    T1 + 3*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(3,2)*N*CO * (  - S - 2*S**2*T1**(-1) + T1
     +     - MS2 )
     +
      M2QGV = M2QGV + SK1C0C(3,2)*N*CK * ( 2*S*T1**(-1)*MS2 - S - T1 + 
     +    5*MS2 )
     +
      M2QGV = M2QGV + SK1C0C(3,2)*CQED * ( 3*S**(-1)*T1*MS2 + S**(-1)*
     +    T1**2 + T1 + 2*MS2 )
     +
      M2QGV = M2QGV + SK1D0(1,1)*N*CO * ( 2*S*T1 + S**2 + T1**2 )
     +
      M2QGV = M2QGV + SK1D0(1,1)*N*CK * ( S*MS2 - T1**2 )
     +
      M2QGV = M2QGV + SK1D0(1,1)*CQED * (  - S*T1 - S*MS2 - S**2*
     +    U1**(-1)*MS2 + S**2 + S**3*U1**(-1) + T1**2 )
     +
      M2QGV = M2QGV + SK1D0(1,2)*N*CK * ( 3*S*T1 - S*MS2 + S**2 + 2*
     +    T1**2 )
     +
      M2QGV = M2QGV + SK1D0(2,1)*N*CO * (  - 2*S*MS2 + S**2 + T1*MS2 )
     +
      M2QGV = M2QGV + SK1D0(2,1)*N*CK * ( S*MS2 - T1*MS2 - 4*MS2**2 )
     +
      M2QGV = M2QGV + SK1D0(2,1)*CQED * ( 4*S*U1**(-1)*MS2**2 - 2*S*MS2
     +     - 2*S**2*U1**(-1)*MS2 + T1*MS2 + 4*MS2**2 )
     +
      M2QGV = M2QGV + SK1D0(2,2)*N*CK * (  - S*T1 + S*MS2 - S**2 + 2*T1
     +    *MS2 + 4*MS2**2 )
     +
      M2QGV = M2QGV + SK1D0(3,1)*N*CO * (  - 3*S*T1 - 2*S**2 + T1*MS2
     +     - T1**2 )
     +
      M2QGV = M2QGV + SK1D0(3,1)*N*CK * (  - S*T1 - 2*S*MS2 - 3*T1*MS2
     +     - T1**2 )
     +
      M2QGV = M2QGV + SK1D0(3,1)*CQED * (  - 3*S**(-1)*T1**2*MS2 - 
     +    S**(-1)*T1**3 - 2*T1*MS2 - T1**2 )
     +
      M2QGV = M2QGV + SOF2(1)*N*CO * ( 16*S**(-1)*T1*MS2 + 16*S*
     +    T1**(-1)*MS2 - 16*S - 8*S**2*T1**(-1) - 8*T1 + 32*MS2 )
     +
      M2QGV = M2QGV + SOF2(1)*N*CK * (  - 16*S**(-1)*T1*MS2 + 8*T1 )
     +
      M2QGV = M2QGV + SOF2(2)*N*CO * ( 16*S**(-1)*T1*MS2 + 32*S*
     +    T1**(-1)*MS2 + 32*MS2 )
     +
      M2QGV = M2QGV + SOF2(2)*N*CK * (  - 16*S**(-1)*T1*MS2 )
     +
      M2QGV = M2QGV + SOF2(3)*N*CO * ( 8*S**(-1)*T1*MS2 + 16*S*T1**(-1)
     +    *MS2 + 16*MS2 )
     +
      M2QGV = M2QGV + SOF2(3)*N*CK * (  - 24*S**(-1)*T1*MS2 - 16*S*
     +    T1**(-1)*MS2 - 16*MS2 )
     +
      M2QGV = M2QGV + SOF2(3)*CQED * ( 8*S**(-1)*T1*MS2 )
     +
      M2QGV = M2QGV + SOF2(4)*N*CO * ( 8*S + 8*S**2*T1**(-1) )
     +
      M2QGV = M2QGV + SOF2(4)*N*CK * (  - 8*S**(-1)*T1**2 - 8*T1 )
     +
      M2QGV = M2QGV + SOF2(5)*N*CK * ( 8*S**(-1)*T1**2 + 16*S + 16*T1 )
     +
      M2QGV = M2QGV + SOF2(5)*CQED * (  - 8*S**(-1)*T1**2 )
     +
      M2QGV = M2QGV + SOF2(6)*N*CO * (  - 8*S**(-1)*T1**2 - 16*S - 16*
     +    T1 )
     +
      M2QGV = M2QGV + SOF2(7)*N*CO * ( 8*S + 8*S**2*T1**(-1) )
     +
      M2QGV = M2QGV + SOF2(7)*N*CK * (  - 8*S**(-1)*T1**2 - 8*T1 )
     +
      M2QGV = M2QGV + SOF2(8)*N*CO * (  - 16*S - 8*S**2*T1**(-1) - 8*T1
     +     )
     +
      M2QGV = M2QGV + SOF2(8)*N*CK * ( 8*T1 )

      ELSE

      M2QGV = 0.D0
      M2QGV = M2QGV + N*CO*THREE**(-1)*FOUR**(-1) * ( 24 + 36*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 - 36*S**(-2)*TG**(-1)*MS2**2*MG2 + 12*
     +    S**(-2)*TG**(-1)*MS2**3 - 12*S**(-2)*TG**(-1)*MG2**3 - 12*
     +    S**(-2)*UG*MS2 + 12*S**(-2)*UG*MG2 - 12*S**(-2)*UG**2*
     +    U1**(-1)*MS2 + 12*S**(-2)*UG**2*U1**(-1)*MG2 - 24*S**(-2)*MS2
     +    *MG2 + 12*S**(-2)*MS2**2 + 12*S**(-2)*MG2**2 - 96*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2**3 + 48*S**(-1)*TG**(-1)*
     +    U1**(-1)*T**(-1)*MS2**2*MG2**2 + 48*S**(-1)*TG**(-1)*U1**(-1)
     +    *T**(-1)*MG2**4 + 60*S**(-1)*TG**(-1)*U1**(-1)*MS2*MG2**2 - 
     +    12*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2 - 12*S**(-1)*TG**(-1)
     +    *U1**(-1)*MS2**3 - 36*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 + 24*
     +    S**(-1)*TG**(-1)*T**(-1)*MS2*MG2**2 - 24*S**(-1)*TG**(-1)*
     +    T**(-1)*MS2**2*MG2 - 180*S**(-1)*TG**(-1)*MS2*MG2 + 96*
     +    S**(-1)*TG**(-1)*MS2**2 + 84*S**(-1)*TG**(-1)*MG2**2 + 24*
     +    S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MG2 - 24*S**(-1)*UG*U1**(-1)*
     +    T**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + N*CO*THREE**(-1)*FOUR**(-1) * ( 48*S**(-1)*UG*
     +    U1**(-1)*MS2 - 48*S**(-1)*UG*U1**(-1)*MG2 + 24*S**(-1)*UG*
     +    T**(-1)*MG2 - 24*S**(-1)*UG - 96*S**(-1)*U1**(-1)*T**(-1)*MS2
     +    *MG2**2 + 96*S**(-1)*U1**(-1)*T**(-1)*MG2**3 + 72*S**(-1)*
     +    U1**(-1)*MS2*MG2 + 12*S**(-1)*U1**(-1)*MS2**2 - 84*S**(-1)*
     +    U1**(-1)*MG2**2 + 24*S**(-1)*T**(-1)*MS2*MG2 - 48*S**(-1)*
     +    T**(-1)*MG2**2 - 84*S**(-1)*MS2 + 108*S**(-1)*MG2 - 96*S*
     +    TG**(-3)*T**(-1)*MG2**3 + 96*S*TG**(-3)*MG2**2 - 96*S*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2**2 + 96*S*TG**(-2)*U1**(-1)
     +    *T**(-1)*MG2**3 + 96*S*TG**(-2)*U1**(-1)*MS2*MG2 - 96*S*
     +    TG**(-2)*U1**(-1)*MG2**2 - 288*S*TG**(-2)*T**(-1)*MG2**2 + 
     +    240*S*TG**(-2)*MG2 - 72*S*TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2
     +     + 72*S*TG**(-1)*U1**(-1)*T**(-1)*MG2**2 - 192*S*TG**(-1)*
     +    T**(-1)*MG2 + 48*S*TG**(-1) + 192*TG**(-3)*T**(-1)*MS2*MG2**3
     +     - 192*TG**(-3)*T**(-1)*MG2**4 - 192*TG**(-3)*MS2*MG2**2 + 
     +    192*TG**(-3)*MG2**3 )
     +
      M2QGV = M2QGV + N*CO*THREE**(-1)*FOUR**(-1) * (  - 96*TG**(-2)*
     +    U1**(-1)*T**(-1)*MS2*MG2**3 + 96*TG**(-2)*U1**(-1)*T**(-1)*
     +    MS2**2*MG2**2 + 96*TG**(-2)*U1**(-1)*MS2*MG2**2 - 96*TG**(-2)
     +    *U1**(-1)*MS2**2*MG2 + 336*TG**(-2)*T**(-1)*MS2*MG2**2 - 384*
     +    TG**(-2)*T**(-1)*MG2**3 - 240*TG**(-2)*MS2*MG2 + 288*TG**(-2)
     +    *MG2**2 - 192*TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2**2 + 72*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MG2 + 120*TG**(-1)*U1**(-1)*
     +    T**(-1)*MG2**3 + 24*TG**(-1)*U1**(-1)*MS2*MG2 + 12*TG**(-1)*
     +    U1**(-1)*MS2**2 - 36*TG**(-1)*U1**(-1)*MG2**2 + 120*TG**(-1)*
     +    T**(-1)*MS2*MG2 - 360*TG**(-1)*T**(-1)*MG2**2 - 108*TG**(-1)*
     +    MS2 + 276*TG**(-1)*MG2 + 48*U1**(-2)*MS2*MG2 - 48*U1**(-2)*
     +    MS2**2 - 48*U1**(-1)*T**(-1)*MS2*MG2 + 24*U1**(-1)*T**(-1)*
     +    MG2**2 - 48*U1**(-1)*MS2 + 48*U1**(-1)*MG2 - 120*T**(-1)*MG2
     +     )
     +
      M2QGV = M2QGV + N*CK*THREE**(-1)*FOUR**(-1) * ( 24 - 108*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 + 108*S**(-2)*TG**(-1)*MS2**2*MG2 - 36*
     +    S**(-2)*TG**(-1)*MS2**3 + 36*S**(-2)*TG**(-1)*MG2**3 + 36*
     +    S**(-2)*UG*MS2 - 36*S**(-2)*UG*MG2 + 36*S**(-2)*UG**2*
     +    U1**(-1)*MS2 - 36*S**(-2)*UG**2*U1**(-1)*MG2 + 72*S**(-2)*MS2
     +    *MG2 - 36*S**(-2)*MS2**2 - 36*S**(-2)*MG2**2 - 36*S**(-1)*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 + 36*S**(-1)*TG**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**3 + 12*S**(-1)
     +    *TG**(-1)*U1**(-1)*MG2**3 + 60*S**(-1)*TG**(-1)*MS2*MG2 - 24*
     +    S**(-1)*TG**(-1)*MS2**2 - 36*S**(-1)*TG**(-1)*MG2**2 + 24*
     +    S**(-1)*UG*U1**(-1)*MS2 - 24*S**(-1)*UG*U1**(-1)*MG2 + 24*
     +    S**(-1)*UG - 48*S**(-1)*U1**(-1)*MS2*MG2 + 12*S**(-1)*
     +    U1**(-1)*MS2**2 + 36*S**(-1)*U1**(-1)*MG2**2 + 84*S**(-1)*MS2
     +     - 84*S**(-1)*MG2 + 48*S*TG**(-1)*U1**(-1)*MS2 - 48*S*
     +    TG**(-1)*U1**(-1)*MG2 + 72*TG**(-1)*U1**(-1)*MS2*MG2 - 36*
     +    TG**(-1)*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + N*CK*THREE**(-1)*FOUR**(-1) * (  - 36*TG**(-1)*
     +    U1**(-1)*MG2**2 + 60*TG**(-1)*MS2 - 60*TG**(-1)*MG2 - 48*
     +    U1**(-2)*MS2*MG2 + 48*U1**(-2)*MS2**2 + 96*U1**(-1)*MS2 - 96*
     +    U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * (  - 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*NS*T**(-1)*MS2*MG2**3 + 96*S**(-1)*TG**(-1)
     +    *U1**(-1)*NS*T**(-1)*MS2**2*MG2**2 - 48*S**(-1)*TG**(-1)*
     +    U1**(-1)*NS*T**(-1)*MS2**3*MG2 + 48*S**(-1)*TG**(-1)*U1**(-1)
     +    *NS*MS2*MG2**2 - 96*S**(-1)*TG**(-1)*U1**(-1)*NS*MS2**2*MG2
     +     + 48*S**(-1)*TG**(-1)*U1**(-1)*NS*MS2**3 - 96*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2*MG2**2 + 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MT2*MG2 + 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2**3 + 96*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MT2*MG2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*
     +    MT2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MT2*MG2**2 - 192*S**(-1)*
     +    TG**(-1)*NS*T**(-1)*MS2*MG2**2 + 72*S**(-1)*TG**(-1)*NS*
     +    T**(-1)*MS2**2*MG2 + 24*S**(-1)*TG**(-1)*NS*T**(-1)*MS2**3 + 
     +    96*S**(-1)*TG**(-1)*NS*T**(-1)*MG2**3 + 216*S**(-1)*TG**(-1)*
     +    NS*MS2*MG2 - 120*S**(-1)*TG**(-1)*NS*MS2**2 - 96*S**(-1)*
     +    TG**(-1)*NS*MG2**2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * ( 24*S**(-1)*TG**(-1)
     +    *T**(-1)*MS2*MT2*MG2 - 24*S**(-1)*TG**(-1)*T**(-1)*MS2**2*MT2
     +     + 24*S**(-1)*TG**(-1)*MS2*MT2 - 24*S**(-1)*TG**(-1)*MT2*MG2
     +     + 120*S**(-1)*UG*U1**(-1)*NS*T**(-1)*MS2*MG2 - 120*S**(-1)*
     +    UG*U1**(-1)*NS*T**(-1)*MG2**2 - 120*S**(-1)*UG*U1**(-1)*NS*
     +    MS2 + 120*S**(-1)*UG*U1**(-1)*NS*MG2 + 24*S**(-1)*UG*U1**(-1)
     +    *T**(-1)*MS2*MT2 - 24*S**(-1)*UG*U1**(-1)*T**(-1)*MT2*MG2 - 
     +    24*S**(-1)*UG*U1**(-1)*MS2*MT2*MG2**(-1) + 24*S**(-1)*UG*
     +    U1**(-1)*MT2 - 24*S**(-1)*UG*NS*T**(-1)*MG2 - 12*S**(-1)*UG*
     +    NS*MS2*MG2**(-1) + 24*S**(-1)*UG*NS + 24*S**(-1)*UG*T**(-1)*
     +    MT2 - 12*S**(-1)*UG*MT2*MG2**(-1) - 24*S**(-1)*UG**2*U1**(-1)
     +    *NS*T**(-1)*MS2 + 24*S**(-1)*UG**2*U1**(-1)*NS*T**(-1)*MG2 + 
     +    24*S**(-1)*UG**2*U1**(-1)*NS*MS2*MG2**(-1) - 24*S**(-1)*UG**2
     +    *U1**(-1)*NS + 48*S**(-1)*U1**(-1)*NS*T**(-1)*MS2**2*MG2 - 48
     +    *S**(-1)*U1**(-1)*NS*T**(-1)*MG2**3 - 48*S**(-1)*U1**(-1)*NS*
     +    MS2**2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * ( 48*S**(-1)*U1**(-1)
     +    *NS*MG2**2 - 96*S**(-1)*U1**(-1)*T**(-1)*MS2*MT2*MG2 + 96*
     +    S**(-1)*U1**(-1)*T**(-1)*MT2*MG2**2 + 96*S**(-1)*U1**(-1)*MS2
     +    *MT2 - 96*S**(-1)*U1**(-1)*MT2*MG2 - 96*S**(-1)*NS*T**(-1)*
     +    MS2*MG2 - 24*S**(-1)*NS*T**(-1)*MS2**2 + 144*S**(-1)*NS*
     +    T**(-1)*MG2**2 + 132*S**(-1)*NS*MS2 - 12*S**(-1)*NS*MS2**2*
     +    MG2**(-1) - 144*S**(-1)*NS*MG2 + 24*S**(-1)*T**(-1)*MS2*MT2
     +     - 48*S**(-1)*T**(-1)*MT2*MG2 + 12*S**(-1)*MS2*MT2*MG2**(-1)
     +     + 12*S**(-1)*MT2 + 96*S*TG**(-3)*NS*T**(-1)*MS2*MG2**2 - 96*
     +    S*TG**(-3)*NS*MS2*MG2 - 96*S*TG**(-3)*T**(-1)*MT2*MG2**2 + 96
     +    *S*TG**(-3)*MT2*MG2 - 192*S*TG**(-2)*U1**(-1)*NS*T**(-1)*MS2*
     +    MG2**2 + 192*S*TG**(-2)*U1**(-1)*NS*T**(-1)*MS2**2*MG2 + 192*
     +    S*TG**(-2)*U1**(-1)*NS*MS2*MG2 - 192*S*TG**(-2)*U1**(-1)*NS*
     +    MS2**2 - 96*S*TG**(-2)*U1**(-1)*T**(-1)*MS2*MT2*MG2 + 96*S*
     +    TG**(-2)*U1**(-1)*T**(-1)*MT2*MG2**2 + 96*S*TG**(-2)*U1**(-1)
     +    *MS2*MT2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * (  - 96*S*TG**(-2)*
     +    U1**(-1)*MT2*MG2 + 192*S*TG**(-2)*NS*T**(-1)*MS2*MG2 + 96*S*
     +    TG**(-2)*NS*T**(-1)*MG2**2 - 144*S*TG**(-2)*NS*MS2 - 96*S*
     +    TG**(-2)*NS*MG2 - 288*S*TG**(-2)*T**(-1)*MT2*MG2 + 240*S*
     +    TG**(-2)*MT2 - 408*S*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2*MG2 + 
     +    144*S*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**2 + 264*S*TG**(-1)*
     +    U1**(-1)*NS*T**(-1)*MG2**2 + 384*S*TG**(-1)*U1**(-1)*NS*MS2
     +     - 96*S*TG**(-1)*U1**(-1)*NS*MS2**2*MG2**(-1) - 288*S*
     +    TG**(-1)*U1**(-1)*NS*MG2 - 72*S*TG**(-1)*U1**(-1)*T**(-1)*MS2
     +    *MT2 + 72*S*TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2 + 48*S*TG**(-1)
     +    *U1**(-1)*MS2*MT2*MG2**(-1) - 48*S*TG**(-1)*U1**(-1)*MT2 + 
     +    120*S*TG**(-1)*NS*T**(-1)*MS2 + 72*S*TG**(-1)*NS*T**(-1)*MG2
     +     + 24*S*TG**(-1)*NS*MS2*MG2**(-1) - 48*S*TG**(-1)*NS - 192*S*
     +    TG**(-1)*T**(-1)*MT2 + 24*S*TG**(-1)*MT2*MG2**(-1) - 72*S*
     +    U1**(-1)*NS*T**(-1)*MS2 + 72*S*U1**(-1)*NS*T**(-1)*MG2 + 48*S
     +    *U1**(-1)*NS*MS2*MG2**(-1) )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * (  - 48*S*U1**(-1)*NS
     +     - 96*S**2*TG**(-2)*U1**(-1)*NS*T**(-1)*MS2*MG2 + 96*S**2*
     +    TG**(-2)*U1**(-1)*NS*T**(-1)*MG2**2 + 96*S**2*TG**(-2)*
     +    U1**(-1)*NS*MS2 - 96*S**2*TG**(-2)*U1**(-1)*NS*MG2 - 72*S**2*
     +    TG**(-1)*U1**(-1)*NS*T**(-1)*MS2 + 72*S**2*TG**(-1)*U1**(-1)*
     +    NS*T**(-1)*MG2 + 48*S**2*TG**(-1)*U1**(-1)*NS*MS2*MG2**(-1)
     +     - 48*S**2*TG**(-1)*U1**(-1)*NS + 192*TG**(-3)*NS*T**(-1)*MS2
     +    *MG2**3 - 192*TG**(-3)*NS*T**(-1)*MS2**2*MG2**2 - 192*
     +    TG**(-3)*NS*MS2*MG2**2 + 192*TG**(-3)*NS*MS2**2*MG2 + 192*
     +    TG**(-3)*T**(-1)*MS2*MT2*MG2**2 - 192*TG**(-3)*T**(-1)*MT2*
     +    MG2**3 - 192*TG**(-3)*MS2*MT2*MG2 + 192*TG**(-3)*MT2*MG2**2
     +     + 96*TG**(-2)*U1**(-1)*NS*T**(-1)*MS2**2*MG2**2 - 96*
     +    TG**(-2)*U1**(-1)*NS*T**(-1)*MS2**3*MG2 - 96*TG**(-2)*
     +    U1**(-1)*NS*MS2**2*MG2 + 96*TG**(-2)*U1**(-1)*NS*MS2**3 - 96*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MT2*MG2**2 + 96*TG**(-2)*
     +    U1**(-1)*T**(-1)*MS2**2*MT2*MG2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * ( 96*TG**(-2)*
     +    U1**(-1)*MS2*MT2*MG2 - 96*TG**(-2)*U1**(-1)*MS2**2*MT2 + 288*
     +    TG**(-2)*NS*T**(-1)*MS2*MG2**2 - 240*TG**(-2)*NS*T**(-1)*
     +    MS2**2*MG2 - 144*TG**(-2)*NS*MS2*MG2 + 96*TG**(-2)*NS*MS2**2
     +     + 336*TG**(-2)*T**(-1)*MS2*MT2*MG2 - 384*TG**(-2)*T**(-1)*
     +    MT2*MG2**2 - 192*TG**(-2)*MS2*MT2 + 240*TG**(-2)*MT2*MG2 - 
     +    408*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2*MG2**2 + 384*TG**(-1)*
     +    U1**(-1)*NS*T**(-1)*MS2**2*MG2 - 72*TG**(-1)*U1**(-1)*NS*
     +    T**(-1)*MS2**3 + 96*TG**(-1)*U1**(-1)*NS*T**(-1)*MG2**3 + 432
     +    *TG**(-1)*U1**(-1)*NS*MS2*MG2 - 384*TG**(-1)*U1**(-1)*NS*
     +    MS2**2 + 48*TG**(-1)*U1**(-1)*NS*MS2**3*MG2**(-1) - 96*
     +    TG**(-1)*U1**(-1)*NS*MG2**2 - 192*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2*MT2*MG2 + 72*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MT2 + 120*
     +    TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2**2 + 192*TG**(-1)*U1**(-1)*
     +    MS2*MT2 - 48*TG**(-1)*U1**(-1)*MS2**2*MT2*MG2**(-1) - 144*
     +    TG**(-1)*U1**(-1)*MT2*MG2 )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * ( 120*TG**(-1)*NS*
     +    T**(-1)*MS2*MG2 - 48*TG**(-1)*NS*T**(-1)*MS2**2 + 168*
     +    TG**(-1)*NS*T**(-1)*MG2**2 + 72*TG**(-1)*NS*MS2 - 48*TG**(-1)
     +    *NS*MS2**2*MG2**(-1) - 192*TG**(-1)*NS*MG2 + 120*TG**(-1)*
     +    T**(-1)*MS2*MT2 - 360*TG**(-1)*T**(-1)*MT2*MG2 + 168*TG**(-1)
     +    *MT2 - 48*UG*U1**(-2)*NS*MS2 + 24*UG*U1**(-2)*NS*MS2**2*
     +    MG2**(-1) + 24*UG*U1**(-2)*NS*MG2 + 48*UG*U1**(-1)*NS*T**(-1)
     +    *MS2 - 48*UG*U1**(-1)*NS*T**(-1)*MG2 - 48*UG*U1**(-1)*NS*MS2*
     +    MG2**(-1) + 48*UG*U1**(-1)*NS - 72*U1**(-2)*NS*MS2*MG2 + 48*
     +    U1**(-2)*NS*MS2**2 + 24*U1**(-2)*NS*MG2**2 + 24*U1**(-2)*MS2*
     +    MT2 - 24*U1**(-2)*MS2**2*MT2*MG2**(-1) - 168*U1**(-1)*NS*
     +    T**(-1)*MS2*MG2 + 72*U1**(-1)*NS*T**(-1)*MS2**2 + 120*
     +    U1**(-1)*NS*T**(-1)*MG2**2 + 216*U1**(-1)*NS*MS2 - 48*
     +    U1**(-1)*NS*MS2**2*MG2**(-1) - 168*U1**(-1)*NS*MG2 - 48*
     +    U1**(-1)*T**(-1)*MS2*MT2 + 24*U1**(-1)*T**(-1)*MT2*MG2 + 24*
     +    U1**(-1)*MS2*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + CO*THREE**(-1)*FOUR**(-1) * (  - 24*U1**(-1)*MT2
     +     + 72*NS*T**(-1)*MS2 + 48*NS*T**(-1)*MG2 + 36*NS*MS2*
     +    MG2**(-1) - 48*NS - 120*T**(-1)*MT2 + 12*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + CK*THREE**(-1)*FOUR**(-1) * (  - 24*S**(-1)*UG*
     +    U1**(-1)*NS*MS2 + 24*S**(-1)*UG*U1**(-1)*NS*MG2 + 24*S**(-1)*
     +    UG*U1**(-1)*MS2*MT2*MG2**(-1) - 24*S**(-1)*UG*U1**(-1)*MT2 + 
     +    12*S**(-1)*UG*NS*MS2*MG2**(-1) - 24*S**(-1)*UG*NS + 12*
     +    S**(-1)*UG*MT2*MG2**(-1) - 24*S**(-1)*UG**2*U1**(-1)*NS*MS2*
     +    MG2**(-1) + 24*S**(-1)*UG**2*U1**(-1)*NS - 12*S**(-1)*NS*MS2
     +     + 12*S**(-1)*NS*MS2**2*MG2**(-1) - 12*S**(-1)*MS2*MT2*
     +    MG2**(-1) + 12*S**(-1)*MT2 + 48*UG*U1**(-2)*NS*MS2 - 24*UG*
     +    U1**(-2)*NS*MS2**2*MG2**(-1) - 24*UG*U1**(-2)*NS*MG2 + 72*
     +    U1**(-2)*NS*MS2*MG2 - 48*U1**(-2)*NS*MS2**2 - 24*U1**(-2)*NS*
     +    MG2**2 - 24*U1**(-2)*MS2*MT2 + 24*U1**(-2)*MS2**2*MT2*
     +    MG2**(-1) - 24*U1**(-1)*NS*MS2 + 24*U1**(-1)*NS*MG2 + 24*
     +    U1**(-1)*MS2*MT2*MG2**(-1) - 24*U1**(-1)*MT2 - 12*NS*MS2*
     +    MG2**(-1) + 12*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + CQED*THREE**(-1)*FOUR**(-1) * ( 36*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 - 36*S**(-2)*TG**(-1)*MS2**2*MG2 + 12*
     +    S**(-2)*TG**(-1)*MS2**3 - 12*S**(-2)*TG**(-1)*MG2**3 - 12*
     +    S**(-2)*UG*MS2 + 12*S**(-2)*UG*MG2 - 12*S**(-2)*UG**2*
     +    U1**(-1)*MS2 + 12*S**(-2)*UG**2*U1**(-1)*MG2 - 24*S**(-2)*MS2
     +    *MG2 + 12*S**(-2)*MS2**2 + 12*S**(-2)*MG2**2 + 36*S**(-1)*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 - 36*S**(-1)*TG**(-1)*U1**(-1)*
     +    MS2**2*MG2 + 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**3 - 12*S**(-1)
     +    *TG**(-1)*U1**(-1)*MG2**3 + 72*S**(-1)*TG**(-1)*MS2*MG2 - 36*
     +    S**(-1)*TG**(-1)*MS2**2 - 36*S**(-1)*TG**(-1)*MG2**2 - 36*
     +    S**(-1)*UG*U1**(-1)*MS2 + 36*S**(-1)*UG*U1**(-1)*MG2 + 36*
     +    S**(-1)*U1**(-1)*MS2*MG2 - 12*S**(-1)*U1**(-1)*MS2**2 - 24*
     +    S**(-1)*U1**(-1)*MG2**2 + 24*S*TG**(-1)*U1**(-1)*MS2 - 24*S*
     +    TG**(-1)*U1**(-1)*MG2 + 72*TG**(-1)*U1**(-1)*MS2*MG2 - 36*
     +    TG**(-1)*U1**(-1)*MS2**2 - 36*TG**(-1)*U1**(-1)*MG2**2 + 24*
     +    TG**(-1)*MS2 )
     +
      M2QGV = M2QGV + CQED*THREE**(-1)*FOUR**(-1) * (  - 24*TG**(-1)*
     +    MG2 + 24*U1**(-1)*MS2 - 24*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(1)*N*CO*THREE**(-1)*FOUR**(-1) * ( 20 + 24
     +    *S**(-2)*TG**(-1)*MS2*MG2**2 - 12*S**(-2)*TG**(-1)*MS2**2*MG2
     +     - 12*S**(-2)*TG**(-1)*MG2**3 + 12*S**(-2)*UG*U1**(-1)*MS2*
     +    MG2 - 12*S**(-2)*UG*U1**(-1)*MG2**2 + 24*S**(-2)*UG*MG2 - 12*
     +    S**(-2)*MS2*MG2 + 12*S**(-2)*MG2**2 - 96*S**(-1)*TG**(-1)*
     +    U1**(-1)*T**(-1)*MS2*MG2**3 + 48*S**(-1)*TG**(-1)*U1**(-1)*
     +    T**(-1)*MS2**2*MG2**2 + 48*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*
     +    MG2**4 + 72*S**(-1)*TG**(-1)*U1**(-1)*MS2*MG2**2 - 36*S**(-1)
     +    *TG**(-1)*U1**(-1)*MS2**2*MG2 - 36*S**(-1)*TG**(-1)*U1**(-1)*
     +    MG2**3 + 24*S**(-1)*TG**(-1)*T**(-1)*MS2*MG2**2 - 24*S**(-1)*
     +    TG**(-1)*T**(-1)*MS2**2*MG2 + 12*S**(-1)*TG**(-1)*MS2*MG2 + 
     +    24*S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MG2 - 24*S**(-1)*UG*
     +    U1**(-1)*T**(-1)*MG2**2 - 40*S**(-1)*UG*U1**(-1)*MS2 + 40*
     +    S**(-1)*UG*U1**(-1)*MG2 + 24*S**(-1)*UG*T**(-1)*MG2 - 20*
     +    S**(-1)*UG - 96*S**(-1)*U1**(-1)*T**(-1)*MS2*MG2**2 + 96*
     +    S**(-1)*U1**(-1)*T**(-1)*MG2**3 )
     +
      M2QGV = M2QGV + SK1B0A(1)*N*CO*THREE**(-1)*FOUR**(-1) * ( 72*
     +    S**(-1)*U1**(-1)*MS2*MG2 - 72*S**(-1)*U1**(-1)*MG2**2 + 24*
     +    S**(-1)*T**(-1)*MS2*MG2 - 48*S**(-1)*T**(-1)*MG2**2 + 20*
     +    S**(-1)*MS2 + 16*S**(-1)*MG2 - 96*S*TG**(-3)*T**(-1)*MG2**3
     +     + 96*S*TG**(-3)*MG2**2 - 96*S*TG**(-2)*U1**(-1)*T**(-1)*MS2*
     +    MG2**2 + 96*S*TG**(-2)*U1**(-1)*T**(-1)*MG2**3 + 96*S*
     +    TG**(-2)*U1**(-1)*MS2*MG2 - 96*S*TG**(-2)*U1**(-1)*MG2**2 - 
     +    288*S*TG**(-2)*T**(-1)*MG2**2 + 240*S*TG**(-2)*MG2 - 72*S*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2 + 72*S*TG**(-1)*U1**(-1)*
     +    T**(-1)*MG2**2 + 80*S*TG**(-1)*U1**(-1)*MS2 - 80*S*TG**(-1)*
     +    U1**(-1)*MG2 - 192*S*TG**(-1)*T**(-1)*MG2 + 40*S*TG**(-1) + 
     +    192*TG**(-3)*T**(-1)*MS2*MG2**3 - 192*TG**(-3)*T**(-1)*MG2**4
     +     - 192*TG**(-3)*MS2*MG2**2 + 192*TG**(-3)*MG2**3 - 96*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2**3 + 96*TG**(-2)*U1**(-1)*
     +    T**(-1)*MS2**2*MG2**2 + 96*TG**(-2)*U1**(-1)*MS2*MG2**2 - 96*
     +    TG**(-2)*U1**(-1)*MS2**2*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(1)*N*CO*THREE**(-1)*FOUR**(-1) * ( 336*
     +    TG**(-2)*T**(-1)*MS2*MG2**2 - 384*TG**(-2)*T**(-1)*MG2**3 - 
     +    224*TG**(-2)*MS2*MG2 + 272*TG**(-2)*MG2**2 - 192*TG**(-1)*
     +    U1**(-1)*T**(-1)*MS2*MG2**2 + 72*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2**2*MG2 + 120*TG**(-1)*U1**(-1)*T**(-1)*MG2**3 + 200*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 80*TG**(-1)*U1**(-1)*MS2**2 - 120
     +    *TG**(-1)*U1**(-1)*MG2**2 + 120*TG**(-1)*T**(-1)*MS2*MG2 - 
     +    360*TG**(-1)*T**(-1)*MG2**2 + 180*TG**(-1)*MG2 + 40*U1**(-2)*
     +    MS2*MG2 - 40*U1**(-2)*MS2**2 - 48*U1**(-1)*T**(-1)*MS2*MG2 + 
     +    24*U1**(-1)*T**(-1)*MG2**2 + 40*U1**(-1)*MS2 - 40*U1**(-1)*
     +    MG2 - 120*T**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(1)*N*CK*THREE**(-1)*FOUR**(-1) * ( 20 - 72
     +    *S**(-2)*TG**(-1)*MS2*MG2**2 + 36*S**(-2)*TG**(-1)*MS2**2*MG2
     +     + 36*S**(-2)*TG**(-1)*MG2**3 - 36*S**(-2)*UG*U1**(-1)*MS2*
     +    MG2 + 36*S**(-2)*UG*U1**(-1)*MG2**2 - 72*S**(-2)*UG*MG2 + 36*
     +    S**(-2)*MS2*MG2 - 36*S**(-2)*MG2**2 - 24*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 + 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     + 12*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 - 36*S**(-1)*TG**(-1)*
     +    MS2*MG2 + 24*S**(-1)*TG**(-1)*MG2**2 + 40*S**(-1)*UG*U1**(-1)
     +    *MS2 - 40*S**(-1)*UG*U1**(-1)*MG2 + 20*S**(-1)*UG - 48*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 72*S**(-1)*U1**(-1)*MG2**2 - 20*
     +    S**(-1)*MS2 - 64*S**(-1)*MG2 - 24*TG**(-1)*U1**(-1)*MS2*MG2
     +     + 24*TG**(-1)*U1**(-1)*MG2**2 - 12*TG**(-1)*MG2 - 40*
     +    U1**(-2)*MS2*MG2 + 40*U1**(-2)*MS2**2 + 40*U1**(-1)*MS2 - 40*
     +    U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(1)*CQED*THREE**(-1)*FOUR**(-1) * ( 24*
     +    S**(-2)*TG**(-1)*MS2*MG2**2 - 12*S**(-2)*TG**(-1)*MS2**2*MG2
     +     - 12*S**(-2)*TG**(-1)*MG2**3 + 12*S**(-2)*UG*U1**(-1)*MS2*
     +    MG2 - 12*S**(-2)*UG*U1**(-1)*MG2**2 + 24*S**(-2)*UG*MG2 - 12*
     +    S**(-2)*MS2*MG2 + 12*S**(-2)*MG2**2 + 24*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 12*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 + 24*S**(-1)*TG**(-1)*
     +    MS2*MG2 - 24*S**(-1)*TG**(-1)*MG2**2 + 36*S**(-1)*U1**(-1)*
     +    MS2*MG2 - 48*S**(-1)*U1**(-1)*MG2**2 + 36*S**(-1)*MG2 + 24*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 24*TG**(-1)*U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CO*THREE**(-1)*FOUR**(-1) * ( 8 - 96*
     +    S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2*MG2**2 + 48*S**(-1)
     +    *TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MT2*MG2 + 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2**3 + 96*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MT2*MG2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*
     +    MT2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MT2*MG2**2 + 24*S**(-1)*
     +    TG**(-1)*T**(-1)*MS2*MT2*MG2 - 24*S**(-1)*TG**(-1)*T**(-1)*
     +    MS2**2*MT2 + 24*S**(-1)*TG**(-1)*MS2*MT2 - 24*S**(-1)*
     +    TG**(-1)*MT2*MG2 + 24*S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MT2 - 
     +    24*S**(-1)*UG*U1**(-1)*T**(-1)*MT2*MG2 - 24*S**(-1)*UG*
     +    U1**(-1)*MS2*MT2*MG2**(-1) - 16*S**(-1)*UG*U1**(-1)*MS2 + 24*
     +    S**(-1)*UG*U1**(-1)*MT2 + 16*S**(-1)*UG*U1**(-1)*MG2 + 24*
     +    S**(-1)*UG*T**(-1)*MT2 - 12*S**(-1)*UG*MT2*MG2**(-1) - 8*
     +    S**(-1)*UG - 96*S**(-1)*U1**(-1)*T**(-1)*MS2*MT2*MG2 + 96*
     +    S**(-1)*U1**(-1)*T**(-1)*MT2*MG2**2 + 96*S**(-1)*U1**(-1)*MS2
     +    *MT2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CO*THREE**(-1)*FOUR**(-1) * (  - 96*
     +    S**(-1)*U1**(-1)*MT2*MG2 + 24*S**(-1)*T**(-1)*MS2*MT2 - 48*
     +    S**(-1)*T**(-1)*MT2*MG2 + 12*S**(-1)*MS2*MT2*MG2**(-1) + 8*
     +    S**(-1)*MS2 + 12*S**(-1)*MT2 - 8*S**(-1)*MG2 - 96*S*TG**(-3)*
     +    T**(-1)*MT2*MG2**2 + 96*S*TG**(-3)*MT2*MG2 - 96*S*TG**(-2)*
     +    U1**(-1)*T**(-1)*MS2*MT2*MG2 + 96*S*TG**(-2)*U1**(-1)*T**(-1)
     +    *MT2*MG2**2 + 96*S*TG**(-2)*U1**(-1)*MS2*MT2 - 96*S*TG**(-2)*
     +    U1**(-1)*MT2*MG2 - 288*S*TG**(-2)*T**(-1)*MT2*MG2 + 240*S*
     +    TG**(-2)*MT2 - 72*S*TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2 + 72*S*
     +    TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2 + 48*S*TG**(-1)*U1**(-1)*
     +    MS2*MT2*MG2**(-1) + 32*S*TG**(-1)*U1**(-1)*MS2 - 48*S*
     +    TG**(-1)*U1**(-1)*MT2 - 32*S*TG**(-1)*U1**(-1)*MG2 - 192*S*
     +    TG**(-1)*T**(-1)*MT2 + 24*S*TG**(-1)*MT2*MG2**(-1) + 16*S*
     +    TG**(-1) + 192*TG**(-3)*T**(-1)*MS2*MT2*MG2**2 - 192*TG**(-3)
     +    *T**(-1)*MT2*MG2**3 - 192*TG**(-3)*MS2*MT2*MG2 + 192*TG**(-3)
     +    *MT2*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CO*THREE**(-1)*FOUR**(-1) * (  - 96*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MT2*MG2**2 + 96*TG**(-2)*
     +    U1**(-1)*T**(-1)*MS2**2*MT2*MG2 + 96*TG**(-2)*U1**(-1)*MS2*
     +    MT2*MG2 - 96*TG**(-2)*U1**(-1)*MS2**2*MT2 + 336*TG**(-2)*
     +    T**(-1)*MS2*MT2*MG2 - 384*TG**(-2)*T**(-1)*MT2*MG2**2 - 192*
     +    TG**(-2)*MS2*MT2 - 32*TG**(-2)*MS2*MG2 + 240*TG**(-2)*MT2*MG2
     +     + 32*TG**(-2)*MG2**2 - 192*TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2
     +    *MG2 + 72*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MT2 + 120*TG**(-1)
     +    *U1**(-1)*T**(-1)*MT2*MG2**2 + 192*TG**(-1)*U1**(-1)*MS2*MT2
     +     + 32*TG**(-1)*U1**(-1)*MS2*MG2 - 48*TG**(-1)*U1**(-1)*MS2**2
     +    *MT2*MG2**(-1) - 32*TG**(-1)*U1**(-1)*MS2**2 - 144*TG**(-1)*
     +    U1**(-1)*MT2*MG2 + 120*TG**(-1)*T**(-1)*MS2*MT2 - 360*
     +    TG**(-1)*T**(-1)*MT2*MG2 + 168*TG**(-1)*MT2 + 24*U1**(-2)*MS2
     +    *MT2 + 16*U1**(-2)*MS2*MG2 - 24*U1**(-2)*MS2**2*MT2*MG2**(-1)
     +     - 16*U1**(-2)*MS2**2 - 48*U1**(-1)*T**(-1)*MS2*MT2 + 24*
     +    U1**(-1)*T**(-1)*MT2*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(2)*CO*THREE**(-1)*FOUR**(-1) * ( 24*
     +    U1**(-1)*MS2*MT2*MG2**(-1) + 16*U1**(-1)*MS2 - 24*U1**(-1)*
     +    MT2 - 16*U1**(-1)*MG2 - 120*T**(-1)*MT2 + 12*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0A(2)*CK*THREE**(-1)*FOUR**(-1) * ( 8 + 24*
     +    S**(-1)*UG*U1**(-1)*MS2*MT2*MG2**(-1) + 16*S**(-1)*UG*
     +    U1**(-1)*MS2 - 24*S**(-1)*UG*U1**(-1)*MT2 - 16*S**(-1)*UG*
     +    U1**(-1)*MG2 + 12*S**(-1)*UG*MT2*MG2**(-1) + 8*S**(-1)*UG - 
     +    12*S**(-1)*MS2*MT2*MG2**(-1) - 8*S**(-1)*MS2 + 12*S**(-1)*MT2
     +     + 8*S**(-1)*MG2 - 24*U1**(-2)*MS2*MT2 - 16*U1**(-2)*MS2*MG2
     +     + 24*U1**(-2)*MS2**2*MT2*MG2**(-1) + 16*U1**(-2)*MS2**2 + 24
     +    *U1**(-1)*MS2*MT2*MG2**(-1) + 16*U1**(-1)*MS2 - 24*U1**(-1)*
     +    MT2 - 16*U1**(-1)*MG2 + 12*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0A(3)*N*CO*FOUR**(-1) * ( 2 - 24*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-2)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-2)*TG**(-1)*MS2**3 + 8*S**(-2)*TG**(-1)*MG2**3 + 8*
     +    S**(-2)*UG*MS2 - 8*S**(-2)*UG*MG2 - 16*S**(-2)*MS2*MG2 + 8*
     +    S**(-2)*MS2**2 + 8*S**(-2)*MG2**2 - 12*S**(-1)*TG**(-1)*MS2*
     +    MG2 + 4*S**(-1)*TG**(-1)*MS2**2 + 8*S**(-1)*TG**(-1)*MG2**2
     +     - 4*S**(-1)*UG*U1**(-1)*MS2 + 4*S**(-1)*UG*U1**(-1)*MG2 - 2*
     +    S**(-1)*UG + 6*S**(-1)*MS2 - 6*S**(-1)*MG2 + 8*S*TG**(-1)*
     +    U1**(-1)*MS2 - 8*S*TG**(-1)*U1**(-1)*MG2 + 4*S*TG**(-1) - 8*
     +    TG**(-2)*MS2*MG2 + 8*TG**(-2)*MG2**2 + 8*TG**(-1)*U1**(-1)*
     +    MS2*MG2 - 8*TG**(-1)*U1**(-1)*MS2**2 + 4*TG**(-1)*MS2 - 4*
     +    TG**(-1)*MG2 + 4*U1**(-2)*MS2*MG2 - 4*U1**(-2)*MS2**2 + 4*
     +    U1**(-1)*MS2 - 4*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(3)*N*CK*FOUR**(-1) * ( 2 + 72*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 - 72*S**(-2)*TG**(-1)*MS2**2*MG2 + 24*
     +    S**(-2)*TG**(-1)*MS2**3 - 24*S**(-2)*TG**(-1)*MG2**3 - 24*
     +    S**(-2)*UG*MS2 + 24*S**(-2)*UG*MG2 + 48*S**(-2)*MS2*MG2 - 24*
     +    S**(-2)*MS2**2 - 24*S**(-2)*MG2**2 + 48*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     + 16*S**(-1)*TG**(-1)*U1**(-1)*MS2**3 - 16*S**(-1)*TG**(-1)*
     +    U1**(-1)*MG2**3 + 140*S**(-1)*TG**(-1)*MS2*MG2 - 64*S**(-1)*
     +    TG**(-1)*MS2**2 - 76*S**(-1)*TG**(-1)*MG2**2 - 56*S**(-1)*UG*
     +    U1**(-1)*MS2 + 56*S**(-1)*UG*U1**(-1)*MG2 + 6*S**(-1)*UG + 32
     +    *S**(-1)*U1**(-1)*MS2*MG2 - 16*S**(-1)*U1**(-1)*MS2**2 - 16*
     +    S**(-1)*U1**(-1)*MG2**2 + 34*S**(-1)*MS2 - 34*S**(-1)*MG2 + 
     +    36*S*TG**(-1)*U1**(-1)*MS2 - 36*S*TG**(-1)*U1**(-1)*MG2 - 4*S
     +    *TG**(-1) + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2 + 120*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 52*TG**(-1)*U1**(-1)*MS2**2 - 68*
     +    TG**(-1)*U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0A(3)*N*CK*FOUR**(-1) * ( 40*TG**(-1)*MS2 - 
     +    40*TG**(-1)*MG2 - 12*U1**(-2)*MS2*MG2 + 12*U1**(-2)*MS2**2 + 
     +    48*U1**(-1)*MS2 - 48*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(3)*CQED*FOUR**(-1) * (  - 2 - 24*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-2)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-2)*TG**(-1)*MS2**3 + 8*S**(-2)*TG**(-1)*MG2**3 + 8*
     +    S**(-2)*UG*MS2 - 8*S**(-2)*UG*MG2 - 16*S**(-2)*MS2*MG2 + 8*
     +    S**(-2)*MS2**2 + 8*S**(-2)*MG2**2 - 24*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 8*S**(-1)*TG**(-1)*U1**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*
     +    U1**(-1)*MG2**3 - 64*S**(-1)*TG**(-1)*MS2*MG2 + 30*S**(-1)*
     +    TG**(-1)*MS2**2 + 34*S**(-1)*TG**(-1)*MG2**2 + 30*S**(-1)*UG*
     +    U1**(-1)*MS2 - 30*S**(-1)*UG*U1**(-1)*MG2 - 2*S**(-1)*UG - 16
     +    *S**(-1)*U1**(-1)*MS2*MG2 + 8*S**(-1)*U1**(-1)*MS2**2 + 8*
     +    S**(-1)*U1**(-1)*MG2**2 - 20*S**(-1)*MS2 + 20*S**(-1)*MG2 - 
     +    22*S*TG**(-1)*U1**(-1)*MS2 + 22*S*TG**(-1)*U1**(-1)*MG2 - 64*
     +    TG**(-1)*U1**(-1)*MS2*MG2 + 30*TG**(-1)*U1**(-1)*MS2**2 + 34*
     +    TG**(-1)*U1**(-1)*MG2**2 - 22*TG**(-1)*MS2 + 22*TG**(-1)*MG2
     +     + 4*U1**(-2)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0A(3)*CQED*FOUR**(-1) * (  - 4*U1**(-2)*
     +    MS2**2 - 26*U1**(-1)*MS2 + 26*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0B(1)*N*CO*THREE**(-1)*FOUR**(-1) * ( 24*
     +    S**(-1)*TG**(-1)*MS2*MG2 - 12*S**(-1)*TG**(-1)*MS2**2 - 12*
     +    S**(-1)*TG**(-1)*MG2**2 - 12*S**(-1)*UG*U1**(-1)*MS2 + 12*
     +    S**(-1)*UG*U1**(-1)*MG2 + 12*S**(-1)*MS2 - 12*S**(-1)*MG2 + 
     +    12*S*TG**(-1)*U1**(-1)*MS2 - 12*S*TG**(-1)*U1**(-1)*MG2 + 48*
     +    S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 24*S**2*TG**(-1)*XLAM**(-2)*
     +    MG2 + 24*TG**(-1)*U1**(-1)*MS2*MG2 - 12*TG**(-1)*U1**(-1)*
     +    MS2**2 - 12*TG**(-1)*U1**(-1)*MG2**2 + 48*TG**(-1)*XLAM**(-2)
     +    *MS2*MG2**2 - 24*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 24*TG**(-1)
     +    *XLAM**(-2)*MG2**3 + 12*TG**(-1)*MS2 + 12*TG**(-1)*MG2 + 48*
     +    UG*XLAM**(-2)*MG2 + 12*U1**(-1)*MS2 - 12*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0B(1)*N*CK*THREE**(-1)*FOUR**(-1) * (  - 12
     +     - 72*S**(-1)*TG**(-1)*MS2*MG2 + 36*S**(-1)*TG**(-1)*MS2**2
     +     + 36*S**(-1)*TG**(-1)*MG2**2 + 36*S**(-1)*UG*U1**(-1)*MS2 - 
     +    36*S**(-1)*UG*U1**(-1)*MG2 - 36*S**(-1)*MS2 + 36*S**(-1)*MG2
     +     - 36*S*TG**(-1)*U1**(-1)*MS2 + 36*S*TG**(-1)*U1**(-1)*MG2 - 
     +    12*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 18*S*TG**(-1)*XLAM**(-2)*
     +    MS2**2 - 18*S*TG**(-1)*XLAM**(-2)*MG2**2 + 6*S*TG**(-1) + 48*
     +    S*UG*XLAM**(-2) + 96*S*U1**(-1)*XLAM**(-2)*MG2**2 - 12*S*
     +    U1**(-1) - 48*S*XLAM**(-2)*MS2 - 96*S*XLAM**(-2)*MG2 + 18*
     +    S**2*TG**(-1)*XLAM**(-2)*MS2 + 6*S**2*TG**(-1)*XLAM**(-2)*MG2
     +     - 6*S**2*U1**(-1)*XLAM**(-2)*MS2 - 66*S**2*U1**(-1)*
     +    XLAM**(-2)*MG2 + 36*S**2*XLAM**(-2) - 6*S**3*TG**(-1)*
     +    XLAM**(-2) + 12*S**3*U1**(-1)*XLAM**(-2) - 72*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 36*TG**(-1)*U1**(-1)*MS2**2 + 36*TG**(-1)*
     +    U1**(-1)*MG2**2 - 30*TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 6*
     +    TG**(-1)*XLAM**(-2)*MS2**2*MG2 )
     +
      M2QGV = M2QGV + SK1B0B(1)*N*CK*THREE**(-1)*FOUR**(-1) * ( 6*
     +    TG**(-1)*XLAM**(-2)*MS2**3 + 18*TG**(-1)*XLAM**(-2)*MG2**3 - 
     +    42*TG**(-1)*MS2 + 18*TG**(-1)*MG2 - 48*UG*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2 + 48*UG*U1**(-1)*XLAM**(-2)*MG2**2 - 42*UG*
     +    XLAM**(-2)*MS2 - 102*UG*XLAM**(-2)*MG2 - 6*UG**2*U1**(-1)*
     +    XLAM**(-2)*MS2 + 6*UG**2*U1**(-1)*XLAM**(-2)*MG2 - 30*
     +    U1**(-1)*MS2 + 78*U1**(-1)*MG2 + 12*XLAM**(-2)*MS2*MG2 + 18*
     +    XLAM**(-2)*MS2**2 - 30*XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0B(1)*CQED*THREE**(-1)*FOUR**(-1) * ( 6 + 24*
     +    S**(-1)*TG**(-1)*MS2*MG2 - 12*S**(-1)*TG**(-1)*MS2**2 - 12*
     +    S**(-1)*TG**(-1)*MG2**2 - 12*S**(-1)*UG*U1**(-1)*MS2 + 12*
     +    S**(-1)*UG*U1**(-1)*MG2 + 12*S**(-1)*MS2 - 12*S**(-1)*MG2 + 
     +    12*S*TG**(-1)*U1**(-1)*MS2 - 12*S*TG**(-1)*U1**(-1)*MG2 - 24*
     +    S*UG*XLAM**(-2) - 48*S*U1**(-1)*XLAM**(-2)*MG2**2 + 12*S*
     +    U1**(-1) + 36*S*XLAM**(-2)*MS2 + 36*S*XLAM**(-2)*MG2 + 6*S**2
     +    *U1**(-1)*XLAM**(-2)*MS2 + 42*S**2*U1**(-1)*XLAM**(-2)*MG2 - 
     +    30*S**2*XLAM**(-2) - 12*S**3*U1**(-1)*XLAM**(-2) + 24*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 12*TG**(-1)*U1**(-1)*MS2**2 - 12*
     +    TG**(-1)*U1**(-1)*MG2**2 + 12*TG**(-1)*MS2 - 12*TG**(-1)*MG2
     +     + 24*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 - 24*UG*U1**(-1)*
     +    XLAM**(-2)*MG2**2 + 18*UG*XLAM**(-2)*MS2 + 30*UG*XLAM**(-2)*
     +    MG2 + 6*UG**2*U1**(-1)*XLAM**(-2)*MS2 - 6*UG**2*U1**(-1)*
     +    XLAM**(-2)*MG2 + 6*U1**(-1)*MS2 - 30*U1**(-1)*MG2 - 12*
     +    XLAM**(-2)*MS2**2 )
     +
      M2QGV = M2QGV + SK1B0B(1)*CQED*THREE**(-1)*FOUR**(-1) * ( 12*
     +    XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0B(5)*N*CO*FOUR**(-1) * ( 12*S**(-2)*TG**(-1)
     +    *MS2*MG2**2 - 12*S**(-2)*TG**(-1)*MS2**2*MG2 + 4*S**(-2)*
     +    TG**(-1)*MS2**3 - 4*S**(-2)*TG**(-1)*MG2**3 - 4*S**(-2)*UG*
     +    MS2 + 4*S**(-2)*UG*MG2 + 4*S**(-2)*UG**2*U1**(-1)*MS2 - 4*
     +    S**(-2)*UG**2*U1**(-1)*MG2 + 24*S**(-2)*MS2*MG2 - 12*S**(-2)*
     +    MS2**2 - 12*S**(-2)*MG2**2 + 12*S**(-1)*TG**(-1)*U1**(-1)*MS2
     +    *MG2**2 - 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2 + 4*S**(-1)
     +    *TG**(-1)*U1**(-1)*MS2**3 - 4*S**(-1)*TG**(-1)*U1**(-1)*
     +    MG2**3 + 48*S**(-1)*TG**(-1)*MS2*MG2 - 20*S**(-1)*TG**(-1)*
     +    MS2**2 - 28*S**(-1)*TG**(-1)*MG2**2 - 16*S**(-1)*UG*U1**(-1)*
     +    MS2 + 16*S**(-1)*UG*U1**(-1)*MG2 + 8*S**(-1)*U1**(-1)*MS2*MG2
     +     - 4*S**(-1)*U1**(-1)*MS2**2 - 4*S**(-1)*U1**(-1)*MG2**2 + 16
     +    *S**(-1)*MS2 - 16*S**(-1)*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2 - 
     +    16*S*TG**(-1)*U1**(-1)*MG2 + 12*S*TG**(-1)*XLAM**(-2)*MS2**2
     +     + 4*S*TG**(-1)*XLAM**(-2)*MG2**2 - 4*S*TG**(-1) - 8*S*UG*
     +    XLAM**(-2) )
     +
      M2QGV = M2QGV + SK1B0B(5)*N*CO*FOUR**(-1) * (  - 12*S**2*TG**(-1)
     +    *XLAM**(-2)*MS2 - 4*S**2*TG**(-1)*XLAM**(-2)*MG2 + 4*S**3*
     +    TG**(-1)*XLAM**(-2) + 40*TG**(-1)*U1**(-1)*MS2*MG2 - 20*
     +    TG**(-1)*U1**(-1)*MS2**2 - 20*TG**(-1)*U1**(-1)*MG2**2 + 4*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 4*TG**(-1)*XLAM**(-2)*MS2**2
     +    *MG2 - 4*TG**(-1)*XLAM**(-2)*MS2**3 - 4*TG**(-1)*XLAM**(-2)*
     +    MG2**3 + 20*TG**(-1)*MS2 - 12*TG**(-1)*MG2 + 8*UG*XLAM**(-2)*
     +    MS2 + 8*UG*XLAM**(-2)*MG2 + 16*U1**(-1)*MS2 - 16*U1**(-1)*MG2
     +     )
     +
      M2QGV = M2QGV + SK1B0B(5)*N*CK*FOUR**(-1) * (  - 36*S**(-2)*
     +    TG**(-1)*MS2*MG2**2 + 36*S**(-2)*TG**(-1)*MS2**2*MG2 - 12*
     +    S**(-2)*TG**(-1)*MS2**3 + 12*S**(-2)*TG**(-1)*MG2**3 + 12*
     +    S**(-2)*UG*MS2 - 12*S**(-2)*UG*MG2 - 12*S**(-2)*UG**2*
     +    U1**(-1)*MS2 + 12*S**(-2)*UG**2*U1**(-1)*MG2 - 72*S**(-2)*MS2
     +    *MG2 + 36*S**(-2)*MS2**2 + 36*S**(-2)*MG2**2 - 36*S**(-1)*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 + 36*S**(-1)*TG**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**3 + 12*S**(-1)
     +    *TG**(-1)*U1**(-1)*MG2**3 - 128*S**(-1)*TG**(-1)*MS2*MG2 + 56
     +    *S**(-1)*TG**(-1)*MS2**2 + 72*S**(-1)*TG**(-1)*MG2**2 + 60*
     +    S**(-1)*UG*U1**(-1)*MS2 - 60*S**(-1)*UG*U1**(-1)*MG2 - 16*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 12*S**(-1)*U1**(-1)*MS2**2 + 4*
     +    S**(-1)*U1**(-1)*MG2**2 - 60*S**(-1)*MS2 + 60*S**(-1)*MG2 - 
     +    44*S*TG**(-1)*U1**(-1)*MS2 + 44*S*TG**(-1)*U1**(-1)*MG2 - 16*
     +    S*TG**(-1)*XLAM**(-2)*MS2*MG2 + 4*S*UG*U1**(-1)*XLAM**(-2)*
     +    MS2 )
     +
      M2QGV = M2QGV + SK1B0B(5)*N*CK*FOUR**(-1) * (  - 4*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2 + 8*S*UG*XLAM**(-2) + 8*S*U1**(-1)*XLAM**(-2)*
     +    MS2*MG2 + 24*S*U1**(-1)*XLAM**(-2)*MG2**2 - 4*S*U1**(-1) - 12
     +    *S*XLAM**(-2)*MS2 - 20*S*XLAM**(-2)*MG2 + 8*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2 - 4*S**2*U1**(-1)*XLAM**(-2)*MS2 - 20*S**2*
     +    U1**(-1)*XLAM**(-2)*MG2 + 8*S**2*XLAM**(-2) + 4*S**3*U1**(-1)
     +    *XLAM**(-2) - 120*TG**(-1)*U1**(-1)*MS2*MG2 + 56*TG**(-1)*
     +    U1**(-1)*MS2**2 + 64*TG**(-1)*U1**(-1)*MG2**2 - 16*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 8*TG**(-1)*XLAM**(-2)*MS2**2*MG2 + 8*
     +    TG**(-1)*XLAM**(-2)*MG2**3 - 44*TG**(-1)*MS2 + 36*TG**(-1)*
     +    MG2 - 16*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 16*UG*U1**(-1)*
     +    XLAM**(-2)*MG2**2 - 4*UG*XLAM**(-2)*MS2 - 44*UG*XLAM**(-2)*
     +    MG2 - 4*UG**2*U1**(-1)*XLAM**(-2)*MS2 + 4*UG**2*U1**(-1)*
     +    XLAM**(-2)*MG2 - 40*U1**(-1)*MS2 + 56*U1**(-1)*MG2 + 8*
     +    XLAM**(-2)*MS2*MG2 + 4*XLAM**(-2)*MS2**2 - 12*XLAM**(-2)*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK1B0B(5)*CQED*FOUR**(-1) * ( 12*S**(-2)*TG**(-1)
     +    *MS2*MG2**2 - 12*S**(-2)*TG**(-1)*MS2**2*MG2 + 4*S**(-2)*
     +    TG**(-1)*MS2**3 - 4*S**(-2)*TG**(-1)*MG2**3 - 4*S**(-2)*UG*
     +    MS2 + 4*S**(-2)*UG*MG2 + 4*S**(-2)*UG**2*U1**(-1)*MS2 - 4*
     +    S**(-2)*UG**2*U1**(-1)*MG2 + 24*S**(-2)*MS2*MG2 - 12*S**(-2)*
     +    MS2**2 - 12*S**(-2)*MG2**2 + 12*S**(-1)*TG**(-1)*U1**(-1)*MS2
     +    *MG2**2 - 12*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2 + 4*S**(-1)
     +    *TG**(-1)*U1**(-1)*MS2**3 - 4*S**(-1)*TG**(-1)*U1**(-1)*
     +    MG2**3 + 40*S**(-1)*TG**(-1)*MS2*MG2 - 18*S**(-1)*TG**(-1)*
     +    MS2**2 - 22*S**(-1)*TG**(-1)*MG2**2 - 22*S**(-1)*UG*U1**(-1)*
     +    MS2 + 22*S**(-1)*UG*U1**(-1)*MG2 + 4*S**(-1)*U1**(-1)*MS2*MG2
     +     - 4*S**(-1)*U1**(-1)*MS2**2 + 22*S**(-1)*MS2 - 22*S**(-1)*
     +    MG2 + 14*S*TG**(-1)*U1**(-1)*MS2 - 14*S*TG**(-1)*U1**(-1)*MG2
     +     - 16*S*U1**(-1)*XLAM**(-2)*MG2**2 + 16*S*XLAM**(-2)*MG2 + 8*
     +    S**2*U1**(-1)*XLAM**(-2)*MG2 + 40*TG**(-1)*U1**(-1)*MS2*MG2
     +     - 18*TG**(-1)*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SK1B0B(5)*CQED*FOUR**(-1) * (  - 22*TG**(-1)*
     +    U1**(-1)*MG2**2 + 14*TG**(-1)*MS2 - 14*TG**(-1)*MG2 + 8*UG*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*UG*U1**(-1)*XLAM**(-2)*MG2**2
     +     + 16*UG*XLAM**(-2)*MG2 + 14*U1**(-1)*MS2 - 22*U1**(-1)*MG2
     +     - 8*XLAM**(-2)*MS2*MG2 + 8*XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CO*FOUR**(-1) * (  - 8 + 8*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 8*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 + 4*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4
     +     - 4*S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 + 4*S**(-1)*UG*
     +    U1**(-1)*MS2 - 4*S**(-1)*UG*U1**(-1)*MG2 - 4*S**(-1)*UG*
     +    XLAM**(-2)*MS2**2 + 4*S**(-1)*UG*XLAM**(-2)*MG2**2 + 4*
     +    S**(-1)*UG - 8*S**(-1)*U1**(-1)*MS2*MG2 + 8*S**(-1)*U1**(-1)*
     +    MG2**2 - 4*S**(-1)*MS2 - 4*S**(-1)*MG2 - 16*S*TG**(-1)*
     +    T1**(-1)*MG2 - 4*S*TG**(-1)*U1**(-1)*MS2 + 4*S*TG**(-1)*
     +    U1**(-1)*MG2 - 8*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 16*S*
     +    TG**(-1)*XLAM**(-2)*MG2**2 + 4*S*TG**(-1) + 4*S*UG*XLAM**(-2)
     +     - 16*S*T1**(-1) + 8*S**2*TG**(-1)*XLAM**(-2)*MS2 + 12*S**2*
     +    TG**(-1)*XLAM**(-2)*MG2 - 4*S**3*TG**(-1)*XLAM**(-2) - 32*
     +    TG**(-2)*MG2**2 + 16*TG**(-1)*T1**(-1)*MG2**2 - 32*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 4*TG**(-1)*U1**(-1)*MS2**2 - 4*TG**(-1)*
     +    U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CO*FOUR**(-1) * (  - 8*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 4*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 8*
     +    TG**(-1)*XLAM**(-2)*MS2**3 + 12*TG**(-1)*XLAM**(-2)*MG2**3 - 
     +    4*TG**(-1)*MS2 - 48*TG**(-1)*MG2 - 8*UG*XLAM**(-2)*MG2 + 8*
     +    T1**(-1)*U1**(-1)*MS2**2 + 8*T1**(-1)*MS2 + 16*T1**(-1)*MG2
     +     - 16*U1**(-2)*MS2*MG2 + 4*U1**(-1)*MS2 - 12*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CK*FOUR**(-1) * ( 6 - 88*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 96*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 40*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 + 4*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 + 28*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 + 32*S**(-1)*TG**(-1)*MS2*
     +    MG2 - 16*S**(-1)*TG**(-1)*MS2**2 - 16*S**(-1)*TG**(-1)*MG2**2
     +     - 12*S**(-1)*UG*U1**(-1)*MS2 + 12*S**(-1)*UG*U1**(-1)*MG2 + 
     +    22*S**(-1)*UG*XLAM**(-2)*MS2*MG2 + S**(-1)*UG*XLAM**(-2)*
     +    MS2**2 - 23*S**(-1)*UG*XLAM**(-2)*MG2**2 - S**(-1)*UG + 24*
     +    S**(-1)*U1**(-1)*MS2*MG2 - 24*S**(-1)*U1**(-1)*MG2**2 + 12*
     +    S**(-1)*MS2 + 12*S**(-1)*MG2 + 16*S*TG**(-1)*T1**(-1)*MG2 - 
     +    64*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 32*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 32*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 + 12*S*TG**(-1)*U1**(-1)*MS2 - 12*S*
     +    TG**(-1)*U1**(-1)*MG2 - 4*S*TG**(-1)*XLAM**(-2)*MS2*MG2 + 26*
     +    S*TG**(-1)*XLAM**(-2)*MS2**2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CK*FOUR**(-1) * ( 34*S*TG**(-1)*
     +    XLAM**(-2)*MG2**2 - 2*S*TG**(-1) - 12*S*UG*U1**(-1)*
     +    XLAM**(-2)*MS2 + 12*S*UG*U1**(-1)*XLAM**(-2)*MG2 - 23*S*UG*
     +    XLAM**(-2) - 8*S*T1**(-1)*U1**(-1)*MS2 + 8*S*T1**(-1)*
     +    U1**(-1)*MG2 + 16*S*T1**(-1) - 40*S*U1**(-1)*XLAM**(-2)*MS2*
     +    MG2 - 40*S*U1**(-1)*XLAM**(-2)*MG2**2 + 8*S*U1**(-1) + 24*S*
     +    XLAM**(-2)*MS2 + 16*S*XLAM**(-2)*MG2 - 16*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 + 16*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**2 - 14*S**2*TG**(-1)*XLAM**(-2)*MS2 - 6*S**2*
     +    TG**(-1)*XLAM**(-2)*MG2 + 10*S**2*U1**(-1)*XLAM**(-2)*MS2 + 
     +    38*S**2*U1**(-1)*XLAM**(-2)*MG2 - 14*S**2*XLAM**(-2) + 2*S**3
     +    *TG**(-1)*XLAM**(-2) - 8*S**3*U1**(-1)*XLAM**(-2) + 32*
     +    TG**(-2)*MG2**2 - 16*TG**(-1)*T1**(-1)*MG2**2 - 48*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**3 + 48*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 16*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*N*CK*FOUR**(-1) * ( 16*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MG2**4 + 48*TG**(-1)*U1**(-1)*MS2*MG2 - 
     +    12*TG**(-1)*U1**(-1)*MS2**2 - 4*TG**(-1)*U1**(-1)*MG2**2 - 54
     +    *TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 66*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 - 18*TG**(-1)*XLAM**(-2)*MS2**3 + 6*TG**(-1)*
     +    XLAM**(-2)*MG2**3 + 22*TG**(-1)*MS2 + 38*TG**(-1)*MG2 + 32*UG
     +    *U1**(-1)*XLAM**(-2)*MS2*MG2 - 32*UG*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 12*UG*XLAM**(-2)*MS2 + 56*UG*XLAM**(-2)*MG2 + 10*
     +    UG**2*U1**(-1)*XLAM**(-2)*MS2 - 10*UG**2*U1**(-1)*XLAM**(-2)*
     +    MG2 + 16*T1**(-1)*U1**(-1)*MS2*MG2 - 16*T1**(-1)*U1**(-1)*
     +    MS2**2 - 8*T1**(-1)*U1**(-1)*MG2**2 - 16*T1**(-1)*MS2 - 8*
     +    T1**(-1)*MG2 + 48*U1**(-2)*MS2*MG2 - 80*U1**(-1)*XLAM**(-2)*
     +    MS2*MG2**2 + 16*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 64*U1**(-1)*
     +    XLAM**(-2)*MG2**3 - 14*U1**(-1)*MS2 + 22*U1**(-1)*MG2 - 8*
     +    XLAM**(-2)*MS2**2 + 8*XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*CQED*FOUR**(-1) * (  - 4 + 40*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 48*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 24*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 4*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 - 12*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 - 16*S**(-1)*TG**(-1)*MS2*
     +    MG2 + 12*S**(-1)*TG**(-1)*MS2**2 + 4*S**(-1)*TG**(-1)*MG2**2
     +     + 8*S**(-1)*UG*U1**(-1)*MS2 - 8*S**(-1)*UG*U1**(-1)*MG2 - 16
     +    *S**(-1)*UG*XLAM**(-2)*MS2*MG2 + 4*S**(-1)*UG*XLAM**(-2)*
     +    MS2**2 + 12*S**(-1)*UG*XLAM**(-2)*MG2**2 - 4*S**(-1)*UG - 8*
     +    S**(-1)*MS2 + 68*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 
     +    52*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 12*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**3 - 28*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 - 12*S*TG**(-1)*U1**(-1)*MS2 + 12*S*
     +    TG**(-1)*U1**(-1)*MG2 + 32*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 12
     +    *S*TG**(-1)*XLAM**(-2)*MS2**2 - 20*S*TG**(-1)*XLAM**(-2)*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(1)*CQED*FOUR**(-1) * ( 6*S*UG*U1**(-1)*
     +    XLAM**(-2)*MS2 - 6*S*UG*U1**(-1)*XLAM**(-2)*MG2 + 12*S*UG*
     +    XLAM**(-2) + 48*S*U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*S*U1**(-1)*
     +    XLAM**(-2)*MS2**2 + 16*S*U1**(-1)*XLAM**(-2)*MG2**2 - 2*S*
     +    U1**(-1) - 22*S*XLAM**(-2)*MS2 - 18*S*XLAM**(-2)*MG2 + 32*
     +    S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2 - 12*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2 - 20*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**2 + 4*S**2*TG**(-1)*XLAM**(-2)*MS2 - 4*S**2*
     +    TG**(-1)*XLAM**(-2)*MG2 + 2*S**2*U1**(-1)*XLAM**(-2)*MS2 - 22
     +    *S**2*U1**(-1)*XLAM**(-2)*MG2 + 12*S**2*XLAM**(-2) + 4*S**3*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2 - 4*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 + 2*S**3*U1**(-1)*XLAM**(-2) + 40*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**3 - 48*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 24*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 4*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**4 - 12*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**4 )
     +
      M2QGV = M2QGV + SK1B0C(1)*CQED*FOUR**(-1) * (  - 16*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 12*TG**(-1)*U1**(-1)*MS2**2 + 4*TG**(-1)*
     +    U1**(-1)*MG2**2 + 68*TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 52*
     +    TG**(-1)*XLAM**(-2)*MS2**2*MG2 + 12*TG**(-1)*XLAM**(-2)*
     +    MS2**3 - 28*TG**(-1)*XLAM**(-2)*MG2**3 - 12*TG**(-1)*MS2 + 12
     +    *TG**(-1)*MG2 - 44*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 4*UG*
     +    U1**(-1)*XLAM**(-2)*MS2**2 + 40*UG*U1**(-1)*XLAM**(-2)*MG2**2
     +     - 6*UG*XLAM**(-2)*MS2 - 34*UG*XLAM**(-2)*MG2 - 10*UG**2*
     +    U1**(-1)*XLAM**(-2)*MS2 + 10*UG**2*U1**(-1)*XLAM**(-2)*MG2 - 
     +    16*U1**(-2)*MS2*MG2 + 52*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 20*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 4*U1**(-1)*XLAM**(-2)*MS2**3
     +     - 36*U1**(-1)*XLAM**(-2)*MG2**3 - 6*U1**(-1)*MS2 + 10*
     +    U1**(-1)*MG2 + 8*XLAM**(-2)*MS2*MG2 + 10*XLAM**(-2)*MS2**2 - 
     +    18*XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(2)*N*CO*FOUR**(-1) * ( 4 + 56*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 72*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 40*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 - 16*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 - 16*S**(-1)*TG**(-1)*MS2*
     +    MG2 + 8*S**(-1)*TG**(-1)*MS2**2 + 8*S**(-1)*TG**(-1)*MG2**2
     +     - 8*S**(-1)*UG*U1**(-1)*MS2 + 8*S**(-1)*UG*U1**(-1)*MG2 - 36
     +    *S**(-1)*UG*XLAM**(-2)*MS2*MG2 + 14*S**(-1)*UG*XLAM**(-2)*
     +    MS2**2 + 22*S**(-1)*UG*XLAM**(-2)*MG2**2 - 18*S**(-1)*UG - 8*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 8*S**(-1)*U1**(-1)*MG2**2 + 4*
     +    S**(-1)*MS2 - 12*S**(-1)*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2 - 
     +    16*S*TG**(-1)*U1**(-1)*MG2 - 8*S*TG**(-1)*XLAM**(-2)*MS2*MG2
     +     - 36*S*TG**(-1)*XLAM**(-2)*MS2**2 - 12*S*TG**(-1)*XLAM**(-2)
     +    *MG2**2 + 8*S*TG**(-1) + 14*S*UG*XLAM**(-2) + 14*S**2*
     +    TG**(-1)*XLAM**(-2)*MS2 + 2*S**2*TG**(-1)*XLAM**(-2)*MG2 - 16
     +    *TG**(-2)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0C(2)*N*CO*FOUR**(-1) * (  - 16*TG**(-2)*
     +    MG2**2 + 16*TG**(-1)*T1**(-1)*MG2**2 - 8*TG**(-1)*U1**(-1)*
     +    MS2*MG2 - 16*TG**(-1)*U1**(-1)*MS2**2 - 8*TG**(-1)*U1**(-1)*
     +    MG2**2 + 2*TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 58*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2 + 30*TG**(-1)*XLAM**(-2)*MS2**3 + 26*
     +    TG**(-1)*XLAM**(-2)*MG2**3 - 14*TG**(-1)*MS2 - 42*TG**(-1)*
     +    MG2 - 28*UG*XLAM**(-2)*MS2 - 36*UG*XLAM**(-2)*MG2 + 8*
     +    T1**(-1)*U1**(-1)*MS2*MG2 + 16*T1**(-1)*MG2 - 8*U1**(-2)*MS2*
     +    MG2 - 8*U1**(-2)*MS2**2 + 8*U1**(-1)*MS2 - 24*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0C(2)*N*CK*FOUR**(-1) * ( 4 - 104*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 120*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 56*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 + 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 + 32*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 + 16*S**(-1)*TG**(-1)*MS2*
     +    MG2 - 8*S**(-1)*TG**(-1)*MS2**2 - 8*S**(-1)*TG**(-1)*MG2**2
     +     + 24*S**(-1)*UG*U1**(-1)*MS2 - 24*S**(-1)*UG*U1**(-1)*MG2 + 
     +    32*S**(-1)*UG*XLAM**(-2)*MS2*MG2 - 4*S**(-1)*UG*XLAM**(-2)*
     +    MS2**2 - 28*S**(-1)*UG*XLAM**(-2)*MG2**2 + 16*S**(-1)*UG + 24
     +    *S**(-1)*U1**(-1)*MS2*MG2 - 24*S**(-1)*U1**(-1)*MG2**2 - 12*
     +    S**(-1)*MS2 + 36*S**(-1)*MG2 - 136*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 104*S*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 - 24*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3 + 56*S*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**3 - 8*S*TG**(-1)*U1**(-1)*
     +    MS2 + 8*S*TG**(-1)*U1**(-1)*MG2 - 56*S*TG**(-1)*XLAM**(-2)*
     +    MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0C(2)*N*CK*FOUR**(-1) * ( 24*S*TG**(-1)*
     +    XLAM**(-2)*MS2**2 + 56*S*TG**(-1)*XLAM**(-2)*MG2**2 - 8*S*
     +    TG**(-1) - 16*S*UG*U1**(-1)*XLAM**(-2)*MS2 + 16*S*UG*U1**(-1)
     +    *XLAM**(-2)*MG2 - 4*S*UG*XLAM**(-2) - 96*S*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 + 16*S*U1**(-1)*XLAM**(-2)*MS2**2 + 16*S*
     +    XLAM**(-2)*MS2 - 64*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2
     +     + 24*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2 + 40*S**2*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**2 - 8*S**2*TG**(-1)*
     +    XLAM**(-2)*MS2 + 24*S**2*U1**(-1)*XLAM**(-2)*MG2 - 8*S**3*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2 + 8*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 + 16*TG**(-2)*MS2*MG2 + 16*TG**(-2)*MG2**2 - 
     +    16*TG**(-1)*T1**(-1)*MG2**2 - 80*TG**(-1)*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2**3 + 96*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2**2
     +     - 48*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3*MG2 + 8*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**4 + 24*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**4 )
     +
      M2QGV = M2QGV + SK1B0C(2)*N*CK*FOUR**(-1) * ( 24*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 8*TG**(-1)*U1**(-1)*MS2**2 - 128*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 112*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 
     +    24*TG**(-1)*XLAM**(-2)*MS2**3 + 40*TG**(-1)*XLAM**(-2)*MG2**3
     +     + 8*TG**(-1)*MS2 + 40*TG**(-1)*MG2 + 72*UG*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 8*UG*U1**(-1)*XLAM**(-2)*MS2**2 - 64*UG*
     +    U1**(-1)*XLAM**(-2)*MG2**2 - 8*UG*XLAM**(-2)*MS2 + 48*UG*
     +    XLAM**(-2)*MG2 + 16*UG**2*U1**(-1)*XLAM**(-2)*MS2 - 16*UG**2*
     +    U1**(-1)*XLAM**(-2)*MG2 - 8*T1**(-1)*U1**(-1)*MS2*MG2 - 16*
     +    T1**(-1)*MG2 + 24*U1**(-2)*MS2*MG2 + 24*U1**(-2)*MS2**2 - 104
     +    *U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 40*U1**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 - 8*U1**(-1)*XLAM**(-2)*MS2**3 + 72*U1**(-1)*
     +    XLAM**(-2)*MG2**3 + 8*U1**(-1)*MS2 + 16*U1**(-1)*MG2 - 24*
     +    XLAM**(-2)*MS2*MG2 - 8*XLAM**(-2)*MS2**2 + 32*XLAM**(-2)*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK1B0C(2)*CQED*FOUR**(-1) * (  - 4 + 24*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 24*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 8*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 8*S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 + 8*S**(-1)
     +    *TG**(-1)*MS2*MG2 - 8*S**(-1)*TG**(-1)*MG2**2 - 8*S**(-1)*UG*
     +    U1**(-1)*MS2 + 8*S**(-1)*UG*U1**(-1)*MG2 - 8*S**(-1)*UG*
     +    XLAM**(-2)*MS2*MG2 + 8*S**(-1)*UG*XLAM**(-2)*MG2**2 - 4*
     +    S**(-1)*UG + 4*S**(-1)*MS2 - 12*S**(-1)*MG2 + 32*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 16*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2 - 16*S*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**3 + 8*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 8*S*TG**(-1)*
     +    XLAM**(-2)*MG2**2 + 8*S*U1**(-1)*XLAM**(-2)*MS2*MG2 + 16*S*
     +    U1**(-1)*XLAM**(-2)*MG2**2 - 8*S*XLAM**(-2)*MG2 + 8*S**2*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MG2**2 - 8*S**2*U1**(-1)*XLAM**(-2)*MG2
     +     + 24*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3 )
     +
      M2QGV = M2QGV + SK1B0C(2)*CQED*FOUR**(-1) * (  - 24*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 8*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 - 8*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**4
     +     + 8*TG**(-1)*U1**(-1)*MS2*MG2 - 8*TG**(-1)*U1**(-1)*MG2**2
     +     + 32*TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 16*TG**(-1)*XLAM**(-2)
     +    *MS2**2*MG2 - 16*TG**(-1)*XLAM**(-2)*MG2**3 - 8*UG*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 + 8*UG*U1**(-1)*XLAM**(-2)*MG2**2 - 8*UG*
     +    XLAM**(-2)*MG2 - 8*U1**(-2)*MS2*MG2 - 8*U1**(-2)*MS2**2 + 40*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 8*U1**(-1)*XLAM**(-2)*MS2**2
     +    *MG2 - 32*U1**(-1)*XLAM**(-2)*MG2**3 - 8*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*N*CK*FOUR**(-1) * (  - 8 - 16*S**(-1)
     +    *TG**(-1)*MS2*MG2 + 16*S**(-1)*TG**(-1)*MG2**2 - 16*S**(-1)*
     +    U1**(-1)*MS2*MG2 + 16*S**(-1)*U1**(-1)*MG2**2 - 16*S*TG**(-1)
     +    *T1**(-1)*MG2 + 8*S*T1**(-1)*U1**(-1)*MS2 - 8*S*T1**(-1)*
     +    U1**(-1)*MG2 - 16*S*T1**(-1) + 32*TG**(-1)*T1**(-1)*MG2**2 - 
     +    16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*U1**(-1)*MG2**2 - 
     +    32*TG**(-1)*MG2 - 8*T1**(-1)*U1**(-1)*MS2*MG2 + 16*T1**(-1)*
     +    U1**(-1)*MS2**2 + 8*T1**(-1)*U1**(-1)*MG2**2 + 16*T1**(-1)*
     +    MS2 + 24*T1**(-1)*MG2 + 16*U1**(-1)*MS2 - 8*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * ( 48*S**(-1)*TG**(-2)
     +    *NS*MS2*MG2**2 - 48*S**(-1)*TG**(-2)*NS*MS2**2*MG2 + 16*
     +    S**(-1)*TG**(-2)*NS*MS2**3 - 16*S**(-1)*TG**(-2)*NS*MG2**3 - 
     +    48*S**(-1)*TG**(-2)*MS2*MG2**2 + 48*S**(-1)*TG**(-2)*MS2**2*
     +    MG2 - 16*S**(-1)*TG**(-2)*MS2**3 + 16*S**(-1)*TG**(-2)*MG2**3
     +     + 16*S**(-1)*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2*MG2**3 - 32*
     +    S**(-1)*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**2*MG2**2 + 16*
     +    S**(-1)*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**3*MG2 + 32*S**(-1)*
     +    TG**(-1)*U1**(-1)*NS*MS2*MG2**2 - 16*S**(-1)*TG**(-1)*
     +    U1**(-1)*NS*MS2**2*MG2 - 16*S**(-1)*TG**(-1)*U1**(-1)*NS*
     +    MG2**3 - 16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2**3 + 32
     +    *S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MG2**2 - 16*S**(-1)
     +    *TG**(-1)*U1**(-1)*T**(-1)*MS2**3*MG2 - 32*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 + 16*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     + 16*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 + 64*S**(-1)*TG**(-1)*
     +    NS*T**(-1)*MS2*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * (  - 24*S**(-1)*
     +    TG**(-1)*NS*T**(-1)*MS2**2*MG2 - 8*S**(-1)*TG**(-1)*NS*
     +    T**(-1)*MS2**3 - 32*S**(-1)*TG**(-1)*NS*T**(-1)*MG2**3 + 72*
     +    S**(-1)*TG**(-1)*NS*MS2*MG2 - 32*S**(-1)*TG**(-1)*NS*MS2**2
     +     - 40*S**(-1)*TG**(-1)*NS*MG2**2 - 64*S**(-1)*TG**(-1)*
     +    T**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*T**(-1)*MS2**2*MG2
     +     + 8*S**(-1)*TG**(-1)*T**(-1)*MS2**3 + 32*S**(-1)*TG**(-1)*
     +    T**(-1)*MG2**3 - 72*S**(-1)*TG**(-1)*MS2*MG2 + 32*S**(-1)*
     +    TG**(-1)*MS2**2 + 40*S**(-1)*TG**(-1)*MG2**2 - 40*S**(-1)*UG*
     +    U1**(-1)*NS*T**(-1)*MS2*MG2 + 40*S**(-1)*UG*U1**(-1)*NS*
     +    T**(-1)*MG2**2 - 8*S**(-1)*UG*U1**(-1)*NS*MS2 + 8*S**(-1)*UG*
     +    U1**(-1)*NS*MG2 + 40*S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MG2 - 40
     +    *S**(-1)*UG*U1**(-1)*T**(-1)*MG2**2 + 8*S**(-1)*UG*U1**(-1)*
     +    MS2 - 8*S**(-1)*UG*U1**(-1)*MG2 + 8*S**(-1)*UG*NS*T**(-1)*MG2
     +     - 8*S**(-1)*UG*T**(-1)*MG2 + 8*S**(-1)*UG**2*U1**(-1)*NS*
     +    T**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * (  - 8*S**(-1)*UG**2*
     +    U1**(-1)*NS*T**(-1)*MG2 - 8*S**(-1)*UG**2*U1**(-1)*T**(-1)*
     +    MS2 + 8*S**(-1)*UG**2*U1**(-1)*T**(-1)*MG2 - 16*S**(-1)*
     +    U1**(-1)*NS*T**(-1)*MS2**2*MG2 + 16*S**(-1)*U1**(-1)*NS*
     +    T**(-1)*MG2**3 + 32*S**(-1)*U1**(-1)*NS*MS2*MG2 - 32*S**(-1)*
     +    U1**(-1)*NS*MG2**2 + 16*S**(-1)*U1**(-1)*T**(-1)*MS2**2*MG2
     +     - 16*S**(-1)*U1**(-1)*T**(-1)*MG2**3 - 32*S**(-1)*U1**(-1)*
     +    MS2*MG2 + 32*S**(-1)*U1**(-1)*MG2**2 + 32*S**(-1)*NS*T**(-1)*
     +    MS2*MG2 + 8*S**(-1)*NS*T**(-1)*MS2**2 - 48*S**(-1)*NS*T**(-1)
     +    *MG2**2 + 16*S**(-1)*NS*MS2 - 8*S**(-1)*NS*MG2 - 32*S**(-1)*
     +    T**(-1)*MS2*MG2 - 8*S**(-1)*T**(-1)*MS2**2 + 48*S**(-1)*
     +    T**(-1)*MG2**2 - 16*S**(-1)*MS2 + 8*S**(-1)*MG2 - 32*S*
     +    TG**(-3)*NS*T**(-1)*MS2*MG2**2 + 48*S*TG**(-3)*NS*MS2*MG2 + 
     +    32*S*TG**(-3)*T**(-1)*MS2*MG2**2 - 48*S*TG**(-3)*MS2*MG2 + 64
     +    *S*TG**(-2)*U1**(-1)*NS*T**(-1)*MS2*MG2**2 - 64*S*TG**(-2)*
     +    U1**(-1)*NS*T**(-1)*MS2**2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * ( 32*S*TG**(-2)*
     +    U1**(-1)*NS*MS2**2 - 32*S*TG**(-2)*U1**(-1)*NS*MG2**2 - 64*S*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2**2 + 64*S*TG**(-2)*U1**(-1)
     +    *T**(-1)*MS2**2*MG2 - 32*S*TG**(-2)*U1**(-1)*MS2**2 + 32*S*
     +    TG**(-2)*U1**(-1)*MG2**2 - 64*S*TG**(-2)*NS*T**(-1)*MS2*MG2
     +     - 32*S*TG**(-2)*NS*T**(-1)*MG2**2 + 48*S*TG**(-2)*NS*MS2 + 
     +    16*S*TG**(-2)*NS*MG2 + 64*S*TG**(-2)*T**(-1)*MS2*MG2 + 32*S*
     +    TG**(-2)*T**(-1)*MG2**2 - 48*S*TG**(-2)*MS2 - 16*S*TG**(-2)*
     +    MG2 + 136*S*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2*MG2 - 48*S*
     +    TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**2 - 88*S*TG**(-1)*U1**(-1)*
     +    NS*T**(-1)*MG2**2 - 136*S*TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2
     +     + 48*S*TG**(-1)*U1**(-1)*T**(-1)*MS2**2 + 88*S*TG**(-1)*
     +    U1**(-1)*T**(-1)*MG2**2 - 40*S*TG**(-1)*NS*T**(-1)*MS2 - 24*S
     +    *TG**(-1)*NS*T**(-1)*MG2 + 40*S*TG**(-1)*T**(-1)*MS2 + 24*S*
     +    TG**(-1)*T**(-1)*MG2 + 24*S*U1**(-1)*NS*T**(-1)*MS2 - 24*S*
     +    U1**(-1)*NS*T**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * (  - 24*S*U1**(-1)*
     +    T**(-1)*MS2 + 24*S*U1**(-1)*T**(-1)*MG2 + 32*S**2*TG**(-2)*
     +    U1**(-1)*NS*T**(-1)*MS2*MG2 - 32*S**2*TG**(-2)*U1**(-1)*NS*
     +    T**(-1)*MG2**2 - 16*S**2*TG**(-2)*U1**(-1)*NS*MS2 + 16*S**2*
     +    TG**(-2)*U1**(-1)*NS*MG2 - 32*S**2*TG**(-2)*U1**(-1)*T**(-1)*
     +    MS2*MG2 + 32*S**2*TG**(-2)*U1**(-1)*T**(-1)*MG2**2 + 16*S**2*
     +    TG**(-2)*U1**(-1)*MS2 - 16*S**2*TG**(-2)*U1**(-1)*MG2 + 24*
     +    S**2*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2 - 24*S**2*TG**(-1)*
     +    U1**(-1)*NS*T**(-1)*MG2 - 24*S**2*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2 + 24*S**2*TG**(-1)*U1**(-1)*T**(-1)*MG2 - 64*TG**(-3)*NS*
     +    T**(-1)*MS2*MG2**3 + 64*TG**(-3)*NS*T**(-1)*MS2**2*MG2**2 + 
     +    64*TG**(-3)*NS*MS2*MG2**2 - 64*TG**(-3)*NS*MS2**2*MG2 + 64*
     +    TG**(-3)*T**(-1)*MS2*MG2**3 - 64*TG**(-3)*T**(-1)*MS2**2*
     +    MG2**2 - 64*TG**(-3)*MS2*MG2**2 + 64*TG**(-3)*MS2**2*MG2 - 32
     +    *TG**(-2)*U1**(-1)*NS*T**(-1)*MS2**2*MG2**2 + 32*TG**(-2)*
     +    U1**(-1)*NS*T**(-1)*MS2**3*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * ( 48*TG**(-2)*
     +    U1**(-1)*NS*MS2*MG2**2 - 16*TG**(-2)*U1**(-1)*NS*MS2**2*MG2
     +     - 16*TG**(-2)*U1**(-1)*NS*MS2**3 - 16*TG**(-2)*U1**(-1)*NS*
     +    MG2**3 + 32*TG**(-2)*U1**(-1)*T**(-1)*MS2**2*MG2**2 - 32*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2**3*MG2 - 48*TG**(-2)*U1**(-1)*
     +    MS2*MG2**2 + 16*TG**(-2)*U1**(-1)*MS2**2*MG2 + 16*TG**(-2)*
     +    U1**(-1)*MS2**3 + 16*TG**(-2)*U1**(-1)*MG2**3 - 96*TG**(-2)*
     +    NS*T**(-1)*MS2*MG2**2 + 80*TG**(-2)*NS*T**(-1)*MS2**2*MG2 + 
     +    104*TG**(-2)*NS*MS2*MG2 - 64*TG**(-2)*NS*MS2**2 - 32*TG**(-2)
     +    *NS*MG2**2 + 96*TG**(-2)*T**(-1)*MS2*MG2**2 - 80*TG**(-2)*
     +    T**(-1)*MS2**2*MG2 - 104*TG**(-2)*MS2*MG2 + 64*TG**(-2)*
     +    MS2**2 + 32*TG**(-2)*MG2**2 + 136*TG**(-1)*U1**(-1)*NS*
     +    T**(-1)*MS2*MG2**2 - 128*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**2*
     +    MG2 + 24*TG**(-1)*U1**(-1)*NS*T**(-1)*MS2**3 - 32*TG**(-1)*
     +    U1**(-1)*NS*T**(-1)*MG2**3 + 64*TG**(-1)*U1**(-1)*NS*MS2*MG2
     +     - 72*TG**(-1)*U1**(-1)*NS*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * (  - 136*TG**(-1)*
     +    U1**(-1)*T**(-1)*MS2*MG2**2 + 128*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2**2*MG2 - 24*TG**(-1)*U1**(-1)*T**(-1)*MS2**3 + 32*
     +    TG**(-1)*U1**(-1)*T**(-1)*MG2**3 - 64*TG**(-1)*U1**(-1)*MS2*
     +    MG2 + 72*TG**(-1)*U1**(-1)*MG2**2 - 40*TG**(-1)*NS*T**(-1)*
     +    MS2*MG2 + 16*TG**(-1)*NS*T**(-1)*MS2**2 - 56*TG**(-1)*NS*
     +    T**(-1)*MG2**2 + 64*TG**(-1)*NS*MS2 - 16*TG**(-1)*NS*MG2 + 40
     +    *TG**(-1)*T**(-1)*MS2*MG2 - 16*TG**(-1)*T**(-1)*MS2**2 + 56*
     +    TG**(-1)*T**(-1)*MG2**2 - 64*TG**(-1)*MS2 + 16*TG**(-1)*MG2
     +     - 16*UG*U1**(-1)*NS*T**(-1)*MS2 + 16*UG*U1**(-1)*NS*T**(-1)*
     +    MG2 + 16*UG*U1**(-1)*T**(-1)*MS2 - 16*UG*U1**(-1)*T**(-1)*MG2
     +     + 56*U1**(-1)*NS*T**(-1)*MS2*MG2 - 24*U1**(-1)*NS*T**(-1)*
     +    MS2**2 - 40*U1**(-1)*NS*T**(-1)*MG2**2 + 16*U1**(-1)*NS*MS2
     +     - 16*U1**(-1)*NS*MG2 - 56*U1**(-1)*T**(-1)*MS2*MG2 + 24*
     +    U1**(-1)*T**(-1)*MS2**2 + 40*U1**(-1)*T**(-1)*MG2**2 - 16*
     +    U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0D(1,1)*CO*FOUR**(-1) * ( 16*U1**(-1)*MG2 - 
     +    24*NS*T**(-1)*MS2 - 16*NS*T**(-1)*MG2 + 24*T**(-1)*MS2 + 16*
     +    T**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CO*FOUR**(-1) * (  - 2*S*TG**(-1)*
     +    U1**(-2)*U**(-1)*MS2*MG2**2 - 2*S*TG**(-1)*U1**(-2)*U**(-1)*
     +    MS2**2*MG2 + 2*S*TG**(-1)*U1**(-2)*U**(-1)*MS2**3 + 2*S*
     +    TG**(-1)*U1**(-2)*U**(-1)*MG2**3 + 2*S*TG**(-1)*U1**(-2)*MS2*
     +    MG2 - S*TG**(-1)*U1**(-2)*MS2**2 - S*TG**(-1)*U1**(-2)*MG2**2
     +     + 2*S*TG**(-1)*U1**(-1)*U**(-1)*MS2*MG2 - 2*S*TG**(-1)*
     +    U1**(-1)*U**(-1)*MG2**2 + 2*S*U1**(-2)*U**(-1)*MS2*MG2 - S*
     +    U1**(-2)*U**(-1)*MS2**2 - S*U1**(-2)*U**(-1)*MG2**2 + 2*S**2*
     +    TG**(-1)*U1**(-2)*U**(-1)*MS2*MG2 - S**2*TG**(-1)*U1**(-2)*
     +    U**(-1)*MS2**2 - S**2*TG**(-1)*U1**(-2)*U**(-1)*MG2**2 + 8*
     +    TG**(-1)*UG**(-1)*MG2**2 + TG**(-1)*U1**(-2)*U**(-1)*MS2**2*
     +    MG2**2 - TG**(-1)*U1**(-2)*U**(-1)*MS2**4 - TG**(-1)*U1**(-2)
     +    *MS2*MG2**2 + TG**(-1)*U1**(-2)*MS2**3 - TG**(-1)*U1**(-1)*
     +    U**(-1)*MS2*MG2**2 - 2*TG**(-1)*U1**(-1)*U**(-1)*MS2**2*MG2
     +     + 2*TG**(-1)*U1**(-1)*U**(-1)*MG2**3 + 2*TG**(-1)*U1**(-1)*
     +    MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CO*FOUR**(-1) * (  - TG**(-1)*
     +    U1**(-1)*MG2**2 - TG**(-1)*U**(-1)*MG2**2 + 8*TG**(-1)*MG2 + 
     +    8*UG**(-1)*U1**(-1)*MG2**2 - 6*UG*U1**(-3)*U**(-1)*MS2*MG2**2
     +     + 6*UG*U1**(-3)*U**(-1)*MS2**2*MG2 - 2*UG*U1**(-3)*U**(-1)*
     +    MS2**3 + 2*UG*U1**(-3)*U**(-1)*MG2**3 + 2*UG*U1**(-2)*U**(-1)
     +    *MS2*MG2 - UG*U1**(-2)*U**(-1)*MS2**2 - UG*U1**(-2)*U**(-1)*
     +    MG2**2 - 8*U1**(-3)*U**(-1)*MS2*MG2**3 + 12*U1**(-3)*U**(-1)*
     +    MS2**2*MG2**2 - 6*U1**(-3)*U**(-1)*MS2**3*MG2 + 2*U1**(-3)*
     +    U**(-1)*MG2**4 - 2*U1**(-3)*MS2**2*MG2 + 2*U1**(-3)*MS2**3 + 
     +    7*U1**(-2)*U**(-1)*MS2*MG2**2 - 8*U1**(-2)*U**(-1)*MS2**2*MG2
     +     + U1**(-2)*U**(-1)*MS2**3 - U1**(-2)*U**(-1)*MG2**3 + 
     +    U1**(-2)*MS2*MG2 + U1**(-2)*MS2**2 - U1**(-2)*MG2**2 - 
     +    U1**(-1)*U**(-1)*MS2*MG2 + 8*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CK*FOUR**(-1) * (  - 8*S**(-1)*
     +    TG**(-1)*MS2**2 + 8*S**(-1)*TG**(-1)*MG2**2 - 12*S**(-1)*UG*
     +    U1**(-2)*U**(-1)*MS2*MG2**2 + 6*S**(-1)*UG*U1**(-2)*U**(-1)*
     +    MS2**2*MG2 + 6*S**(-1)*UG*U1**(-2)*U**(-1)*MG2**3 + 6*S**(-1)
     +    *UG*U1**(-2)*MS2*MG2 - 3*S**(-1)*UG*U1**(-2)*MS2**2 - 3*
     +    S**(-1)*UG*U1**(-2)*MG2**2 + 3*S**(-1)*UG*U1**(-1)*U**(-1)*
     +    MS2*MG2 - 3*S**(-1)*UG*U1**(-1)*U**(-1)*MG2**2 - 8*S**(-1)*UG
     +    *U1**(-1)*MS2 + 8*S**(-1)*UG*U1**(-1)*MG2 - 6*S**(-1)*UG**2*
     +    U1**(-2)*U**(-1)*MS2*MG2 + 3*S**(-1)*UG**2*U1**(-2)*U**(-1)*
     +    MS2**2 + 3*S**(-1)*UG**2*U1**(-2)*U**(-1)*MG2**2 - 6*S**(-1)*
     +    U1**(-2)*U**(-1)*MS2*MG2**3 + 3*S**(-1)*U1**(-2)*U**(-1)*
     +    MS2**2*MG2**2 + 3*S**(-1)*U1**(-2)*U**(-1)*MG2**4 + 6*S**(-1)
     +    *U1**(-2)*MS2*MG2**2 - 3*S**(-1)*U1**(-2)*MS2**2*MG2 - 3*
     +    S**(-1)*U1**(-2)*MG2**3 + 3*S**(-1)*U1**(-1)*U**(-1)*MS2*
     +    MG2**2 - 3*S**(-1)*U1**(-1)*U**(-1)*MG2**3 - 19*S**(-1)*
     +    U1**(-1)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CK*FOUR**(-1) * ( 19*S**(-1)*
     +    U1**(-1)*MG2**2 + 8*S**(-1)*MS2 - 8*S**(-1)*MG2 + 4*S*
     +    TG**(-1)*U1**(-2)*U**(-1)*MS2*MG2**2 + 4*S*TG**(-1)*U1**(-2)*
     +    U**(-1)*MS2**2*MG2 - 4*S*TG**(-1)*U1**(-2)*U**(-1)*MS2**3 - 4
     +    *S*TG**(-1)*U1**(-2)*U**(-1)*MG2**3 - 4*S*TG**(-1)*U1**(-2)*
     +    MS2*MG2 + 2*S*TG**(-1)*U1**(-2)*MS2**2 + 2*S*TG**(-1)*
     +    U1**(-2)*MG2**2 - 4*S*TG**(-1)*U1**(-1)*U**(-1)*MS2*MG2 + 4*S
     +    *TG**(-1)*U1**(-1)*U**(-1)*MG2**2 + 8*S*TG**(-1)*U1**(-1)*MS2
     +     - 8*S*TG**(-1)*U1**(-1)*MG2 - 4*S*U1**(-2)*U**(-1)*MS2*MG2
     +     + 2*S*U1**(-2)*U**(-1)*MS2**2 + 2*S*U1**(-2)*U**(-1)*MG2**2
     +     - 4*S**2*TG**(-1)*U1**(-2)*U**(-1)*MS2*MG2 + 2*S**2*TG**(-1)
     +    *U1**(-2)*U**(-1)*MS2**2 + 2*S**2*TG**(-1)*U1**(-2)*U**(-1)*
     +    MG2**2 - 8*TG**(-1)*UG**(-1)*MG2**2 - 2*TG**(-1)*U1**(-2)*
     +    U**(-1)*MS2**2*MG2**2 + 2*TG**(-1)*U1**(-2)*U**(-1)*MS2**4 + 
     +    2*TG**(-1)*U1**(-2)*MS2*MG2**2 - 2*TG**(-1)*U1**(-2)*MS2**3
     +     + 2*TG**(-1)*U1**(-1)*U**(-1)*MS2*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CK*FOUR**(-1) * ( 4*TG**(-1)*
     +    U1**(-1)*U**(-1)*MS2**2*MG2 - 4*TG**(-1)*U1**(-1)*U**(-1)*
     +    MG2**3 - 4*TG**(-1)*U1**(-1)*MS2*MG2 - 8*TG**(-1)*U1**(-1)*
     +    MS2**2 + 10*TG**(-1)*U1**(-1)*MG2**2 + 2*TG**(-1)*U**(-1)*
     +    MG2**2 + 8*TG**(-1)*MS2 - 16*TG**(-1)*MG2 - 24*UG**(-1)*
     +    U1**(-1)*MG2**2 + 30*UG*U1**(-3)*U**(-1)*MS2*MG2**2 - 30*UG*
     +    U1**(-3)*U**(-1)*MS2**2*MG2 + 10*UG*U1**(-3)*U**(-1)*MS2**3
     +     - 10*UG*U1**(-3)*U**(-1)*MG2**3 - 16*UG*U1**(-2)*U**(-1)*MS2
     +    *MG2 + 8*UG*U1**(-2)*U**(-1)*MS2**2 + 8*UG*U1**(-2)*U**(-1)*
     +    MG2**2 + 40*U1**(-3)*U**(-1)*MS2*MG2**3 - 60*U1**(-3)*U**(-1)
     +    *MS2**2*MG2**2 + 30*U1**(-3)*U**(-1)*MS2**3*MG2 - 10*U1**(-3)
     +    *U**(-1)*MG2**4 + 10*U1**(-3)*MS2**2*MG2 - 10*U1**(-3)*MS2**3
     +     - 50*U1**(-2)*U**(-1)*MS2*MG2**2 + 43*U1**(-2)*U**(-1)*
     +    MS2**2*MG2 - 2*U1**(-2)*U**(-1)*MS2**3 + 14*U1**(-2)*U**(-1)*
     +    MG2**3 + U1**(-2)*MS2*MG2 - 8*U1**(-2)*MS2**2 + 2*U1**(-2)*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*N*CK*FOUR**(-1) * ( 11*U1**(-1)*
     +    U**(-1)*MS2*MG2 - 6*U1**(-1)*U**(-1)*MG2**2 + 8*U1**(-1)*MS2
     +     - 32*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*CQED*FOUR**(-1) * ( 8*S**(-1)*UG*
     +    U1**(-2)*U**(-1)*MS2*MG2**2 - 4*S**(-1)*UG*U1**(-2)*U**(-1)*
     +    MS2**2*MG2 - 4*S**(-1)*UG*U1**(-2)*U**(-1)*MG2**3 - 4*S**(-1)
     +    *UG*U1**(-2)*MS2*MG2 + 2*S**(-1)*UG*U1**(-2)*MS2**2 + 2*
     +    S**(-1)*UG*U1**(-2)*MG2**2 - 2*S**(-1)*UG*U1**(-1)*U**(-1)*
     +    MS2*MG2 + 2*S**(-1)*UG*U1**(-1)*U**(-1)*MG2**2 + 4*S**(-1)*
     +    UG**2*U1**(-2)*U**(-1)*MS2*MG2 - 2*S**(-1)*UG**2*U1**(-2)*
     +    U**(-1)*MS2**2 - 2*S**(-1)*UG**2*U1**(-2)*U**(-1)*MG2**2 + 4*
     +    S**(-1)*U1**(-2)*U**(-1)*MS2*MG2**3 - 2*S**(-1)*U1**(-2)*
     +    U**(-1)*MS2**2*MG2**2 - 2*S**(-1)*U1**(-2)*U**(-1)*MG2**4 - 4
     +    *S**(-1)*U1**(-2)*MS2*MG2**2 + 2*S**(-1)*U1**(-2)*MS2**2*MG2
     +     + 2*S**(-1)*U1**(-2)*MG2**3 - 2*S**(-1)*U1**(-1)*U**(-1)*MS2
     +    *MG2**2 + 2*S**(-1)*U1**(-1)*U**(-1)*MG2**3 + 2*S**(-1)*
     +    U1**(-1)*MS2*MG2 - 2*S**(-1)*U1**(-1)*MG2**2 + 4*S*UG**(-1)*
     +    T1**(-1)*MS2 - 4*S*UG**(-1)*T1**(-1)*MG2 + 8*UG**(-1)*
     +    U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(1,2)*CQED*FOUR**(-1) * ( 4*UG**(-1)*MS2 - 
     +    4*UG**(-1)*MG2 - 12*UG*U1**(-3)*U**(-1)*MS2*MG2**2 + 12*UG*
     +    U1**(-3)*U**(-1)*MS2**2*MG2 - 4*UG*U1**(-3)*U**(-1)*MS2**3 + 
     +    4*UG*U1**(-3)*U**(-1)*MG2**3 + 8*UG*U1**(-2)*U**(-1)*MS2*MG2
     +     - 4*UG*U1**(-2)*U**(-1)*MS2**2 - 4*UG*U1**(-2)*U**(-1)*
     +    MG2**2 + 4*T1**(-1)*MS2 - 4*T1**(-1)*MG2 - 16*U1**(-3)*
     +    U**(-1)*MS2*MG2**3 + 24*U1**(-3)*U**(-1)*MS2**2*MG2**2 - 12*
     +    U1**(-3)*U**(-1)*MS2**3*MG2 + 4*U1**(-3)*U**(-1)*MG2**4 - 4*
     +    U1**(-3)*MS2**2*MG2 + 4*U1**(-3)*MS2**3 + 24*U1**(-2)*U**(-1)
     +    *MS2*MG2**2 - 18*U1**(-2)*U**(-1)*MS2**2*MG2 - 8*U1**(-2)*
     +    U**(-1)*MG2**3 - 2*U1**(-2)*MS2*MG2 + 4*U1**(-2)*MS2**2 - 6*
     +    U1**(-1)*U**(-1)*MS2*MG2 + 4*U1**(-1)*U**(-1)*MG2**2 + 8*
     +    U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(2,1)*N*CO*FOUR**(-1) * ( 8 - 64*S**(-1)*
     +    TG**(-2)*MS2*MG2**2 + 32*S**(-1)*TG**(-2)*MS2**2*MG2 + 32*
     +    S**(-1)*TG**(-2)*MG2**3 + 32*S**(-1)*TG**(-1)*U1**(-1)*
     +    T**(-1)*MS2*MG2**3 - 16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2**2*MG2**2 - 16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MG2**4
     +     - 96*S**(-1)*TG**(-1)*U1**(-1)*MS2*MG2**2 + 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 + 48*S**(-1)*TG**(-1)*U1**(-1)*
     +    MG2**3 - 8*S**(-1)*TG**(-1)*T**(-1)*MS2*MG2**2 + 8*S**(-1)*
     +    TG**(-1)*T**(-1)*MS2**2*MG2 - 72*S**(-1)*TG**(-1)*MS2*MG2 + 8
     +    *S**(-1)*TG**(-1)*MS2**2 + 64*S**(-1)*TG**(-1)*MG2**2 - 8*
     +    S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MG2 + 8*S**(-1)*UG*U1**(-1)*
     +    T**(-1)*MG2**2 + 8*S**(-1)*UG*U1**(-1)*MS2 - 8*S**(-1)*UG*
     +    U1**(-1)*MG2 - 8*S**(-1)*UG*T**(-1)*MG2 + 32*S**(-1)*U1**(-1)
     +    *T**(-1)*MS2*MG2**2 - 32*S**(-1)*U1**(-1)*T**(-1)*MG2**3 - 80
     +    *S**(-1)*U1**(-1)*MS2*MG2 + 80*S**(-1)*U1**(-1)*MG2**2 - 8*
     +    S**(-1)*T**(-1)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(2,1)*N*CO*FOUR**(-1) * ( 16*S**(-1)*
     +    T**(-1)*MG2**2 - 8*S**(-1)*MS2 + 32*S*TG**(-3)*T**(-1)*MG2**3
     +     - 48*S*TG**(-3)*MG2**2 + 32*S*TG**(-2)*U1**(-1)*T**(-1)*MS2*
     +    MG2**2 - 32*S*TG**(-2)*U1**(-1)*T**(-1)*MG2**3 - 64*S*
     +    TG**(-2)*U1**(-1)*MS2*MG2 + 64*S*TG**(-2)*U1**(-1)*MG2**2 + 
     +    96*S*TG**(-2)*T**(-1)*MG2**2 - 64*S*TG**(-2)*MG2 + 16*S*
     +    TG**(-1)*T1**(-1)*MG2 + 24*S*TG**(-1)*U1**(-1)*T**(-1)*MS2*
     +    MG2 - 24*S*TG**(-1)*U1**(-1)*T**(-1)*MG2**2 - 8*S*TG**(-1)*
     +    U1**(-1)*MS2 + 8*S*TG**(-1)*U1**(-1)*MG2 + 64*S*TG**(-1)*
     +    T**(-1)*MG2 + 16*S*T1**(-1) - 64*TG**(-3)*T**(-1)*MS2*MG2**3
     +     + 64*TG**(-3)*T**(-1)*MG2**4 + 64*TG**(-3)*MS2*MG2**2 - 64*
     +    TG**(-3)*MG2**3 + 32*TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2**3 - 
     +    32*TG**(-2)*U1**(-1)*T**(-1)*MS2**2*MG2**2 - 96*TG**(-2)*
     +    U1**(-1)*MS2*MG2**2 + 64*TG**(-2)*U1**(-1)*MS2**2*MG2 + 32*
     +    TG**(-2)*U1**(-1)*MG2**3 - 112*TG**(-2)*T**(-1)*MS2*MG2**2 + 
     +    128*TG**(-2)*T**(-1)*MG2**3 )
     +
      M2QGV = M2QGV + SK1B0D(2,1)*N*CO*FOUR**(-1) * ( 32*TG**(-2)*MS2*
     +    MG2 - 40*TG**(-2)*MG2**2 - 32*TG**(-1)*T1**(-1)*MG2**2 + 64*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2**2 - 24*TG**(-1)*U1**(-1)*
     +    T**(-1)*MS2**2*MG2 - 40*TG**(-1)*U1**(-1)*T**(-1)*MG2**3 - 
     +    136*TG**(-1)*U1**(-1)*MS2*MG2 + 8*TG**(-1)*U1**(-1)*MS2**2 + 
     +    136*TG**(-1)*U1**(-1)*MG2**2 - 40*TG**(-1)*T**(-1)*MS2*MG2 + 
     +    120*TG**(-1)*T**(-1)*MG2**2 - 8*TG**(-1)*MS2 - 8*TG**(-1)*MG2
     +     - 8*T1**(-1)*U1**(-1)*MS2*MG2 - 8*T1**(-1)*U1**(-1)*MS2**2
     +     - 8*T1**(-1)*MS2 - 32*T1**(-1)*MG2 + 16*U1**(-1)*T**(-1)*MS2
     +    *MG2 - 8*U1**(-1)*T**(-1)*MG2**2 - 16*U1**(-1)*MS2 + 8*
     +    U1**(-1)*MG2 + 40*T**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(2,2)*N*CO*FOUR**(-1) * ( 8*TG**(-1)*
     +    UG**(-1)*MG2**2 + 8*TG**(-1)*MG2 + 8*UG**(-1)*U1**(-1)*MG2**2
     +     + 8*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(2,2)*N*CK*FOUR**(-1) * (  - 16*S**(-1)*
     +    TG**(-1)*MS2*MG2 + 16*S**(-1)*TG**(-1)*MG2**2 - 16*S**(-1)*
     +    U1**(-1)*MS2*MG2 + 16*S**(-1)*U1**(-1)*MG2**2 - 8*TG**(-1)*
     +    UG**(-1)*MG2**2 - 16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*
     +    U1**(-1)*MG2**2 - 8*TG**(-1)*MG2 - 24*UG**(-1)*U1**(-1)*
     +    MG2**2 - 24*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(2,2)*CQED*FOUR**(-1) * ( 8*UG**(-1)*
     +    U1**(-1)*MG2**2 + 8*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * ( 32*S**(-1)*TG**(-2)
     +    *MS2*MT2*MG2 + 48*S**(-1)*TG**(-2)*MS2*MG2**2 - 16*S**(-1)*
     +    TG**(-2)*MS2**2*MT2 - 48*S**(-1)*TG**(-2)*MS2**2*MG2 + 16*
     +    S**(-1)*TG**(-2)*MS2**3 - 16*S**(-1)*TG**(-2)*MT2*MG2**2 - 16
     +    *S**(-1)*TG**(-2)*MG2**3 + 32*S**(-1)*TG**(-1)*U1**(-1)*
     +    T**(-1)*MS2*MT2*MG2**2 + 16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)
     +    *MS2*MG2**3 - 16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MT2
     +    *MG2 - 32*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MG2**2 + 
     +    16*S**(-1)*TG**(-1)*U1**(-1)*T**(-1)*MS2**3*MG2 - 16*S**(-1)*
     +    TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2**3 + 32*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 16*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 16*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 - 8*S**(-1)*TG**(-1)*
     +    T**(-1)*MS2*MT2*MG2 + 64*S**(-1)*TG**(-1)*T**(-1)*MS2*MG2**2
     +     + 8*S**(-1)*TG**(-1)*T**(-1)*MS2**2*MT2 - 24*S**(-1)*
     +    TG**(-1)*T**(-1)*MS2**2*MG2 - 8*S**(-1)*TG**(-1)*T**(-1)*
     +    MS2**3 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * (  - 32*S**(-1)*
     +    TG**(-1)*T**(-1)*MG2**3 + 24*S**(-1)*TG**(-1)*MS2*MT2 + 72*
     +    S**(-1)*TG**(-1)*MS2*MG2 - 32*S**(-1)*TG**(-1)*MS2**2 - 24*
     +    S**(-1)*TG**(-1)*MT2*MG2 - 40*S**(-1)*TG**(-1)*MG2**2 - 8*
     +    S**(-1)*UG*U1**(-1)*T**(-1)*MS2*MT2 - 40*S**(-1)*UG*U1**(-1)*
     +    T**(-1)*MS2*MG2 + 8*S**(-1)*UG*U1**(-1)*T**(-1)*MT2*MG2 + 40*
     +    S**(-1)*UG*U1**(-1)*T**(-1)*MG2**2 - 8*S**(-1)*UG*U1**(-1)*
     +    MS2 + 8*S**(-1)*UG*U1**(-1)*MG2 - 8*S**(-1)*UG*T**(-1)*MT2 + 
     +    8*S**(-1)*UG*T**(-1)*MG2 + 8*S**(-1)*UG**2*U1**(-1)*T**(-1)*
     +    MS2 - 8*S**(-1)*UG**2*U1**(-1)*T**(-1)*MG2 + 32*S**(-1)*
     +    U1**(-1)*T**(-1)*MS2*MT2*MG2 - 16*S**(-1)*U1**(-1)*T**(-1)*
     +    MS2**2*MG2 - 32*S**(-1)*U1**(-1)*T**(-1)*MT2*MG2**2 + 16*
     +    S**(-1)*U1**(-1)*T**(-1)*MG2**3 + 32*S**(-1)*U1**(-1)*MS2*MG2
     +     - 32*S**(-1)*U1**(-1)*MG2**2 - 8*S**(-1)*T**(-1)*MS2*MT2 + 
     +    32*S**(-1)*T**(-1)*MS2*MG2 + 8*S**(-1)*T**(-1)*MS2**2 + 16*
     +    S**(-1)*T**(-1)*MT2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * (  - 48*S**(-1)*
     +    T**(-1)*MG2**2 + 16*S**(-1)*MS2 - 8*S**(-1)*MT2 - 8*S**(-1)*
     +    MG2 - 32*S*TG**(-3)*T**(-1)*MS2*MG2**2 + 32*S*TG**(-3)*
     +    T**(-1)*MT2*MG2**2 + 48*S*TG**(-3)*MS2*MG2 - 48*S*TG**(-3)*
     +    MT2*MG2 + 32*S*TG**(-2)*U1**(-1)*T**(-1)*MS2*MT2*MG2 + 64*S*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2**2 - 64*S*TG**(-2)*U1**(-1)
     +    *T**(-1)*MS2**2*MG2 - 32*S*TG**(-2)*U1**(-1)*T**(-1)*MT2*
     +    MG2**2 - 16*S*TG**(-2)*U1**(-1)*MS2*MT2 + 32*S*TG**(-2)*
     +    U1**(-1)*MS2**2 + 16*S*TG**(-2)*U1**(-1)*MT2*MG2 - 32*S*
     +    TG**(-2)*U1**(-1)*MG2**2 - 64*S*TG**(-2)*T**(-1)*MS2*MG2 + 96
     +    *S*TG**(-2)*T**(-1)*MT2*MG2 - 32*S*TG**(-2)*T**(-1)*MG2**2 + 
     +    48*S*TG**(-2)*MS2 - 64*S*TG**(-2)*MT2 + 16*S*TG**(-2)*MG2 + 
     +    24*S*TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2 + 136*S*TG**(-1)*
     +    U1**(-1)*T**(-1)*MS2*MG2 - 48*S*TG**(-1)*U1**(-1)*T**(-1)*
     +    MS2**2 - 24*S*TG**(-1)*U1**(-1)*T**(-1)*MT2*MG2 - 88*S*
     +    TG**(-1)*U1**(-1)*T**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * (  - 40*S*TG**(-1)*
     +    T**(-1)*MS2 + 64*S*TG**(-1)*T**(-1)*MT2 - 24*S*TG**(-1)*
     +    T**(-1)*MG2 + 24*S*U1**(-1)*T**(-1)*MS2 - 24*S*U1**(-1)*
     +    T**(-1)*MG2 + 32*S**2*TG**(-2)*U1**(-1)*T**(-1)*MS2*MG2 - 32*
     +    S**2*TG**(-2)*U1**(-1)*T**(-1)*MG2**2 - 16*S**2*TG**(-2)*
     +    U1**(-1)*MS2 + 16*S**2*TG**(-2)*U1**(-1)*MG2 + 24*S**2*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2 - 24*S**2*TG**(-1)*U1**(-1)*
     +    T**(-1)*MG2 - 64*TG**(-3)*T**(-1)*MS2*MT2*MG2**2 - 64*
     +    TG**(-3)*T**(-1)*MS2*MG2**3 + 64*TG**(-3)*T**(-1)*MS2**2*
     +    MG2**2 + 64*TG**(-3)*T**(-1)*MT2*MG2**3 + 64*TG**(-3)*MS2*MT2
     +    *MG2 + 64*TG**(-3)*MS2*MG2**2 - 64*TG**(-3)*MS2**2*MG2 - 64*
     +    TG**(-3)*MT2*MG2**2 + 32*TG**(-2)*U1**(-1)*T**(-1)*MS2*MT2*
     +    MG2**2 - 32*TG**(-2)*U1**(-1)*T**(-1)*MS2**2*MT2*MG2 - 32*
     +    TG**(-2)*U1**(-1)*T**(-1)*MS2**2*MG2**2 + 32*TG**(-2)*
     +    U1**(-1)*T**(-1)*MS2**3*MG2 + 48*TG**(-2)*U1**(-1)*MS2*MG2**2
     +     + 16*TG**(-2)*U1**(-1)*MS2**2*MT2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * (  - 16*TG**(-2)*
     +    U1**(-1)*MS2**2*MG2 - 16*TG**(-2)*U1**(-1)*MS2**3 - 16*
     +    TG**(-2)*U1**(-1)*MT2*MG2**2 - 16*TG**(-2)*U1**(-1)*MG2**3 - 
     +    112*TG**(-2)*T**(-1)*MS2*MT2*MG2 - 96*TG**(-2)*T**(-1)*MS2*
     +    MG2**2 + 80*TG**(-2)*T**(-1)*MS2**2*MG2 + 128*TG**(-2)*
     +    T**(-1)*MT2*MG2**2 + 80*TG**(-2)*MS2*MT2 + 104*TG**(-2)*MS2*
     +    MG2 - 64*TG**(-2)*MS2**2 - 88*TG**(-2)*MT2*MG2 - 32*TG**(-2)*
     +    MG2**2 + 64*TG**(-1)*U1**(-1)*T**(-1)*MS2*MT2*MG2 + 136*
     +    TG**(-1)*U1**(-1)*T**(-1)*MS2*MG2**2 - 24*TG**(-1)*U1**(-1)*
     +    T**(-1)*MS2**2*MT2 - 128*TG**(-1)*U1**(-1)*T**(-1)*MS2**2*MG2
     +     + 24*TG**(-1)*U1**(-1)*T**(-1)*MS2**3 - 40*TG**(-1)*U1**(-1)
     +    *T**(-1)*MT2*MG2**2 - 32*TG**(-1)*U1**(-1)*T**(-1)*MG2**3 + 8
     +    *TG**(-1)*U1**(-1)*MS2*MT2 + 64*TG**(-1)*U1**(-1)*MS2*MG2 - 
     +    72*TG**(-1)*U1**(-1)*MG2**2 - 40*TG**(-1)*T**(-1)*MS2*MT2 - 
     +    40*TG**(-1)*T**(-1)*MS2*MG2 + 16*TG**(-1)*T**(-1)*MS2**2 + 
     +    120*TG**(-1)*T**(-1)*MT2*MG2 )
     +
      M2QGV = M2QGV + SK1B0D(3,1)*CO*FOUR**(-1) * (  - 56*TG**(-1)*
     +    T**(-1)*MG2**2 + 64*TG**(-1)*MS2 - 48*TG**(-1)*MT2 - 16*
     +    TG**(-1)*MG2 - 16*UG*U1**(-1)*T**(-1)*MS2 + 16*UG*U1**(-1)*
     +    T**(-1)*MG2 + 16*U1**(-1)*T**(-1)*MS2*MT2 + 56*U1**(-1)*
     +    T**(-1)*MS2*MG2 - 24*U1**(-1)*T**(-1)*MS2**2 - 8*U1**(-1)*
     +    T**(-1)*MT2*MG2 - 40*U1**(-1)*T**(-1)*MG2**2 + 16*U1**(-1)*
     +    MS2 - 16*U1**(-1)*MG2 - 24*T**(-1)*MS2 + 40*T**(-1)*MT2 - 16*
     +    T**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0E(1)*N*CO*FOUR**(-1) * (  - 16 + 64*S**(-1)*
     +    TG**(-2)*MS2*MG2**2 - 32*S**(-1)*TG**(-2)*MS2**2*MG2 - 32*
     +    S**(-1)*TG**(-2)*MG2**3 + 64*S**(-1)*TG**(-1)*U1**(-1)*MS2*
     +    MG2**2 - 32*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MG2 - 32*S**(-1)
     +    *TG**(-1)*U1**(-1)*MG2**3 - 64*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2*MG2**3 + 72*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**2*MG2**2 - 
     +    32*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**3*MG2 + 4*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2**4 + 20*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MG2**4 + 120*S**(-1)*TG**(-1)*MS2*MG2 - 36*S**(-1)*TG**(-1)*
     +    MS2**2 - 84*S**(-1)*TG**(-1)*MG2**2 + 36*S**(-1)*UG*
     +    XLAM**(-2)*MS2*MG2 - 10*S**(-1)*UG*XLAM**(-2)*MS2**2 - 26*
     +    S**(-1)*UG*XLAM**(-2)*MG2**2 + 26*S**(-1)*UG + 64*S**(-1)*
     +    U1**(-1)*MS2*MG2 - 64*S**(-1)*U1**(-1)*MG2**2 + 16*S**(-1)*
     +    MS2 + 16*S*TG**(-3)*MG2**2 + 32*S*TG**(-2)*U1**(-1)*MS2*MG2
     +     - 32*S*TG**(-2)*U1**(-1)*MG2**2 - 16*S*TG**(-2)*MG2 - 32*S*
     +    TG**(-1)*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0E(1)*N*CO*FOUR**(-1) * ( 32*S*TG**(-1)*
     +    U1**(-1)*MG2 + 24*S*TG**(-1)*XLAM**(-2)*MS2**2 + 24*S*
     +    TG**(-1)*XLAM**(-2)*MG2**2 - 32*S*TG**(-1) - 10*S*UG*
     +    XLAM**(-2) - 10*S**2*TG**(-1)*XLAM**(-2)*MS2 - 2*S**2*
     +    TG**(-1)*XLAM**(-2)*MG2 + 64*TG**(-2)*U1**(-1)*MS2*MG2**2 - 
     +    32*TG**(-2)*U1**(-1)*MS2**2*MG2 - 32*TG**(-2)*U1**(-1)*MG2**3
     +     + 80*TG**(-2)*MS2*MG2 - 24*TG**(-2)*MG2**2 - 16*TG**(-1)*
     +    UG**(-1)*MG2**2 + 120*TG**(-1)*U1**(-1)*MS2*MG2 + 32*TG**(-1)
     +    *U1**(-1)*MS2**2 - 96*TG**(-1)*U1**(-1)*MG2**2 - 14*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 58*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 
     +    18*TG**(-1)*XLAM**(-2)*MS2**3 - 26*TG**(-1)*XLAM**(-2)*MG2**3
     +     + 42*TG**(-1)*MS2 - 6*TG**(-1)*MG2 - 16*UG**(-1)*U1**(-1)*
     +    MG2**2 + 20*UG*XLAM**(-2)*MS2 + 20*UG*XLAM**(-2)*MG2 + 32*
     +    U1**(-2)*MS2**2 + 16*U1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK1B0E(1)*N*CK*FOUR**(-1) * (  - 12 + 64*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 72*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 32*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 4*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 - 20*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 - 56*S**(-1)*TG**(-1)*MS2*
     +    MG2 + 36*S**(-1)*TG**(-1)*MS2**2 + 20*S**(-1)*TG**(-1)*MG2**2
     +     - 12*S**(-1)*UG*XLAM**(-2)*MS2*MG2 - 2*S**(-1)*UG*XLAM**(-2)
     +    *MS2**2 + 14*S**(-1)*UG*XLAM**(-2)*MG2**2 - 14*S**(-1)*UG - 
     +    16*S**(-1)*MS2 + 100*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*
     +    MG2**2 - 68*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 12*S*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3 - 44*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 - 36*S*TG**(-1)*U1**(-1)*MS2 + 36*S*
     +    TG**(-1)*U1**(-1)*MG2 + 40*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 12
     +    *S*TG**(-1)*XLAM**(-2)*MS2**2 - 28*S*TG**(-1)*XLAM**(-2)*
     +    MG2**2 + 16*S*UG*U1**(-1)*XLAM**(-2)*MS2 - 16*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2 )
     +
      M2QGV = M2QGV + SK1B0E(1)*N*CK*FOUR**(-1) * (  - 2*S*UG*
     +    XLAM**(-2) + 72*S*U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*S*U1**(-1)*
     +    XLAM**(-2)*MS2**2 - 16*S*U1**(-1)*XLAM**(-2)*MG2**2 - 8*S*
     +    XLAM**(-2)*MS2 + 24*S*XLAM**(-2)*MG2 + 40*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 - 12*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2 - 28*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 4*S**2*TG**(-1)*XLAM**(-2)*MS2 - 4*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2 - 4*S**2*U1**(-1)*XLAM**(-2)*MS2 - 8*S**2*
     +    U1**(-1)*XLAM**(-2)*MG2 - 4*S**2*XLAM**(-2) + 4*S**3*TG**(-1)
     +    *U1**(-1)*XLAM**(-2)*MS2 - 4*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 + 64*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3
     +     - 72*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 32*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3*MG2 - 4*TG**(-1)*U1**(-1)
     +    *XLAM**(-2)*MS2**4 - 20*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**4
     +     - 56*TG**(-1)*U1**(-1)*MS2*MG2 + 36*TG**(-1)*U1**(-1)*MS2**2
     +     + 20*TG**(-1)*U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(1)*N*CK*FOUR**(-1) * ( 100*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 - 68*TG**(-1)*XLAM**(-2)*MS2**2*MG2 + 
     +    12*TG**(-1)*XLAM**(-2)*MS2**3 - 44*TG**(-1)*XLAM**(-2)*MG2**3
     +     - 36*TG**(-1)*MS2 + 36*TG**(-1)*MG2 + 16*UG**(-1)*U1**(-1)*
     +    MG2**2 - 36*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 4*UG*U1**(-1)*
     +    XLAM**(-2)*MS2**2 + 32*UG*U1**(-1)*XLAM**(-2)*MG2**2 + 16*UG*
     +    XLAM**(-2)*MS2 - 8*UG*XLAM**(-2)*MG2 - 12*UG**2*U1**(-1)*
     +    XLAM**(-2)*MS2 + 12*UG**2*U1**(-1)*XLAM**(-2)*MG2 - 32*
     +    U1**(-2)*MS2**2 + 92*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 28*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 4*U1**(-1)*XLAM**(-2)*MS2**3
     +     - 68*U1**(-1)*XLAM**(-2)*MG2**3 - 60*U1**(-1)*MS2 + 56*
     +    U1**(-1)*MG2 + 4*XLAM**(-2)*MS2*MG2 + 4*XLAM**(-2)*MS2**2 - 8
     +    *XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*N*CK*FOUR**(-1) * ( 2 + 128*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 - 144*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 64*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 - 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 - 40*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 + 16*S**(-1)*TG**(-1)*MS2*
     +    MG2 + 8*S**(-1)*TG**(-1)*MS2**2 - 24*S**(-1)*TG**(-1)*MG2**2
     +     - 42*S**(-1)*UG*XLAM**(-2)*MS2*MG2 + 5*S**(-1)*UG*XLAM**(-2)
     +    *MS2**2 + 37*S**(-1)*UG*XLAM**(-2)*MG2**2 - 5*S**(-1)*UG - 32
     +    *S**(-1)*MG2 + 100*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2
     +     - 68*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 12*S*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3 - 44*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 - 4*S*TG**(-1)*U1**(-1)*MS2 + 4*S*TG**(-1)*
     +    U1**(-1)*MG2 + 40*S*TG**(-1)*XLAM**(-2)*MS2*MG2 - 32*S*
     +    TG**(-1)*XLAM**(-2)*MS2**2 - 56*S*TG**(-1)*XLAM**(-2)*MG2**2
     +     + 8*S*UG*U1**(-1)*XLAM**(-2)*MS2 - 8*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*N*CK*FOUR**(-1) * ( 5*S*UG*XLAM**(-2)
     +     + 56*S*U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*S*U1**(-1)*XLAM**(-2)
     +    *MS2**2 - 4*S*XLAM**(-2)*MS2 + 12*S*XLAM**(-2)*MG2 + 40*S**2*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2 - 12*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2 - 28*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**2 + 12*S**2*TG**(-1)*XLAM**(-2)*MS2 - 12*S**2
     +    *U1**(-1)*XLAM**(-2)*MG2 - 2*S**2*XLAM**(-2) + 4*S**3*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2 - 4*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 - 64*TG**(-2)*MG2**2 + 16*TG**(-1)*UG**(-1)*
     +    MG2**2 + 64*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3 - 72*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 32*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**3*MG2 - 4*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**4 - 20*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**4 - 
     +    24*TG**(-1)*U1**(-1)*MS2*MG2 + 4*TG**(-1)*U1**(-1)*MS2**2 - 
     +    44*TG**(-1)*U1**(-1)*MG2**2 + 108*TG**(-1)*XLAM**(-2)*MS2*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*N*CK*FOUR**(-1) * (  - 120*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2 + 28*TG**(-1)*XLAM**(-2)*MS2**3 - 16*
     +    TG**(-1)*XLAM**(-2)*MG2**3 - 12*TG**(-1)*MS2 - 32*TG**(-1)*
     +    MG2 + 32*UG**(-1)*U1**(-1)*MG2**2 - 36*UG*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2 + 4*UG*U1**(-1)*XLAM**(-2)*MS2**2 + 32*UG*U1**(-1)*
     +    XLAM**(-2)*MG2**2 - 2*UG*XLAM**(-2)*MS2 - 18*UG*XLAM**(-2)*
     +    MG2 - 8*UG**2*U1**(-1)*XLAM**(-2)*MS2 + 8*UG**2*U1**(-1)*
     +    XLAM**(-2)*MG2 - 64*U1**(-2)*MS2*MG2 + 92*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2**2 - 28*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 4*U1**(-1)*
     +    XLAM**(-2)*MS2**3 - 68*U1**(-1)*XLAM**(-2)*MG2**3 - 20*
     +    U1**(-1)*MG2 + 8*XLAM**(-2)*MS2*MG2 + 2*XLAM**(-2)*MS2**2 - 
     +    10*XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CO*FOUR**(-1) * (  - 20 - 48*S**(-1)*
     +    TG**(-2)*NS*MS2*MG2**2 + 48*S**(-1)*TG**(-2)*NS*MS2**2*MG2 - 
     +    16*S**(-1)*TG**(-2)*NS*MS2**3 + 16*S**(-1)*TG**(-2)*NS*MG2**3
     +     + 48*S**(-1)*TG**(-2)*MS2*MG2**2 - 48*S**(-1)*TG**(-2)*
     +    MS2**2*MG2 + 16*S**(-1)*TG**(-2)*MS2**3 - 16*S**(-1)*TG**(-2)
     +    *MG2**3 - 48*S**(-1)*TG**(-1)*U1**(-1)*NS*MS2*MG2**2 + 48*
     +    S**(-1)*TG**(-1)*U1**(-1)*NS*MS2**2*MG2 - 16*S**(-1)*TG**(-1)
     +    *U1**(-1)*NS*MS2**3 + 16*S**(-1)*TG**(-1)*U1**(-1)*NS*MG2**3
     +     + 48*S**(-1)*TG**(-1)*U1**(-1)*MS2*MG2**2 - 48*S**(-1)*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 + 16*S**(-1)*TG**(-1)*U1**(-1)*
     +    MS2**3 - 16*S**(-1)*TG**(-1)*U1**(-1)*MG2**3 - 128*S**(-1)*
     +    TG**(-1)*NS*MS2*MG2 + 64*S**(-1)*TG**(-1)*NS*MS2**2 + 64*
     +    S**(-1)*TG**(-1)*NS*MG2**2 + 128*S**(-1)*TG**(-1)*MS2*MG2 - 
     +    64*S**(-1)*TG**(-1)*MS2**2 - 64*S**(-1)*TG**(-1)*MG2**2 + 32*
     +    S**(-1)*UG*U1**(-1)*NS*MS2 - 32*S**(-1)*UG*U1**(-1)*NS*MG2 - 
     +    32*S**(-1)*UG*U1**(-1)*MS2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CO*FOUR**(-1) * ( 32*S**(-1)*UG*
     +    U1**(-1)*MG2 + 4*S**(-1)*UG*NS*MS2*MG2**(-1) - 12*S**(-1)*UG*
     +    NS - 4*S**(-1)*UG*MS2*MG2**(-1) + 12*S**(-1)*UG - 8*S**(-1)*
     +    UG**2*U1**(-1)*NS*MS2*MG2**(-1) + 8*S**(-1)*UG**2*U1**(-1)*NS
     +     + 8*S**(-1)*UG**2*U1**(-1)*MS2*MG2**(-1) - 8*S**(-1)*UG**2*
     +    U1**(-1) - 32*S**(-1)*U1**(-1)*NS*MS2*MG2 + 16*S**(-1)*
     +    U1**(-1)*NS*MS2**2 + 16*S**(-1)*U1**(-1)*NS*MG2**2 + 32*
     +    S**(-1)*U1**(-1)*MS2*MG2 - 16*S**(-1)*U1**(-1)*MS2**2 - 16*
     +    S**(-1)*U1**(-1)*MG2**2 - 48*S**(-1)*NS*MS2 + 4*S**(-1)*NS*
     +    MS2**2*MG2**(-1) + 44*S**(-1)*NS*MG2 + 48*S**(-1)*MS2 - 4*
     +    S**(-1)*MS2**2*MG2**(-1) - 44*S**(-1)*MG2 - 16*S*TG**(-3)*NS*
     +    MS2*MG2 + 16*S*TG**(-3)*MS2*MG2 - 64*S*TG**(-2)*U1**(-1)*NS*
     +    MS2*MG2 + 32*S*TG**(-2)*U1**(-1)*NS*MS2**2 + 32*S*TG**(-2)*
     +    U1**(-1)*NS*MG2**2 + 64*S*TG**(-2)*U1**(-1)*MS2*MG2 - 32*S*
     +    TG**(-2)*U1**(-1)*MS2**2 - 32*S*TG**(-2)*U1**(-1)*MG2**2 + 16
     +    *S*TG**(-2)*NS*MG2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CO*FOUR**(-1) * (  - 16*S*TG**(-2)*MG2
     +     - 104*S*TG**(-1)*U1**(-1)*NS*MS2 + 32*S*TG**(-1)*U1**(-1)*NS
     +    *MS2**2*MG2**(-1) + 72*S*TG**(-1)*U1**(-1)*NS*MG2 + 104*S*
     +    TG**(-1)*U1**(-1)*MS2 - 32*S*TG**(-1)*U1**(-1)*MS2**2*
     +    MG2**(-1) - 72*S*TG**(-1)*U1**(-1)*MG2 - 8*S*TG**(-1)*NS*MS2*
     +    MG2**(-1) + 24*S*TG**(-1)*NS + 8*S*TG**(-1)*MS2*MG2**(-1) - 
     +    24*S*TG**(-1) - 16*S*U1**(-1)*NS*MS2*MG2**(-1) + 16*S*
     +    U1**(-1)*NS + 16*S*U1**(-1)*MS2*MG2**(-1) - 16*S*U1**(-1) - 
     +    16*S**2*TG**(-2)*U1**(-1)*NS*MS2 + 16*S**2*TG**(-2)*U1**(-1)*
     +    NS*MG2 + 16*S**2*TG**(-2)*U1**(-1)*MS2 - 16*S**2*TG**(-2)*
     +    U1**(-1)*MG2 - 16*S**2*TG**(-1)*U1**(-1)*NS*MS2*MG2**(-1) + 
     +    16*S**2*TG**(-1)*U1**(-1)*NS + 16*S**2*TG**(-1)*U1**(-1)*MS2*
     +    MG2**(-1) - 16*S**2*TG**(-1)*U1**(-1) - 48*TG**(-2)*U1**(-1)*
     +    NS*MS2*MG2**2 + 48*TG**(-2)*U1**(-1)*NS*MS2**2*MG2 - 16*
     +    TG**(-2)*U1**(-1)*NS*MS2**3 + 16*TG**(-2)*U1**(-1)*NS*MG2**3
     +     + 48*TG**(-2)*U1**(-1)*MS2*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CO*FOUR**(-1) * (  - 48*TG**(-2)*
     +    U1**(-1)*MS2**2*MG2 + 16*TG**(-2)*U1**(-1)*MS2**3 - 16*
     +    TG**(-2)*U1**(-1)*MG2**3 - 72*TG**(-2)*NS*MS2*MG2 + 32*
     +    TG**(-2)*NS*MS2**2 + 48*TG**(-2)*NS*MG2**2 + 72*TG**(-2)*MS2*
     +    MG2 - 32*TG**(-2)*MS2**2 - 48*TG**(-2)*MG2**2 - 176*TG**(-1)*
     +    U1**(-1)*NS*MS2*MG2 + 104*TG**(-1)*U1**(-1)*NS*MS2**2 - 16*
     +    TG**(-1)*U1**(-1)*NS*MS2**3*MG2**(-1) + 96*TG**(-1)*U1**(-1)*
     +    NS*MG2**2 + 176*TG**(-1)*U1**(-1)*MS2*MG2 - 104*TG**(-1)*
     +    U1**(-1)*MS2**2 + 16*TG**(-1)*U1**(-1)*MS2**3*MG2**(-1) - 96*
     +    TG**(-1)*U1**(-1)*MG2**2 - 80*TG**(-1)*NS*MS2 + 16*TG**(-1)*
     +    NS*MS2**2*MG2**(-1) + 72*TG**(-1)*NS*MG2 + 80*TG**(-1)*MS2 - 
     +    16*TG**(-1)*MS2**2*MG2**(-1) - 72*TG**(-1)*MG2 + 16*UG*
     +    U1**(-2)*NS*MS2 - 8*UG*U1**(-2)*NS*MS2**2*MG2**(-1) - 8*UG*
     +    U1**(-2)*NS*MG2 - 16*UG*U1**(-2)*MS2 + 8*UG*U1**(-2)*MS2**2*
     +    MG2**(-1) + 8*UG*U1**(-2)*MG2 + 16*UG*U1**(-1)*NS*MS2*
     +    MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(2)*CO*FOUR**(-1) * (  - 16*UG*U1**(-1)*NS
     +     - 16*UG*U1**(-1)*MS2*MG2**(-1) + 16*UG*U1**(-1) + 32*
     +    U1**(-2)*NS*MS2*MG2 - 24*U1**(-2)*NS*MS2**2 - 8*U1**(-2)*NS*
     +    MG2**2 - 32*U1**(-2)*MS2*MG2 + 24*U1**(-2)*MS2**2 + 8*
     +    U1**(-2)*MG2**2 - 72*U1**(-1)*NS*MS2 + 16*U1**(-1)*NS*MS2**2*
     +    MG2**(-1) + 56*U1**(-1)*NS*MG2 + 72*U1**(-1)*MS2 - 16*
     +    U1**(-1)*MS2**2*MG2**(-1) - 56*U1**(-1)*MG2 - 12*NS*MS2*
     +    MG2**(-1) + 20*NS + 12*MS2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(2)*CK*FOUR**(-1) * (  - 4 + 16*S**(-1)*UG*
     +    U1**(-1)*NS*MS2 - 16*S**(-1)*UG*U1**(-1)*NS*MG2 - 16*S**(-1)*
     +    UG*U1**(-1)*MS2 + 16*S**(-1)*UG*U1**(-1)*MG2 - 4*S**(-1)*UG*
     +    NS*MS2*MG2**(-1) + 12*S**(-1)*UG*NS + 4*S**(-1)*UG*MS2*
     +    MG2**(-1) - 12*S**(-1)*UG + 8*S**(-1)*UG**2*U1**(-1)*NS*MS2*
     +    MG2**(-1) - 8*S**(-1)*UG**2*U1**(-1)*NS - 8*S**(-1)*UG**2*
     +    U1**(-1)*MS2*MG2**(-1) + 8*S**(-1)*UG**2*U1**(-1) - 4*S**(-1)
     +    *NS*MS2**2*MG2**(-1) + 4*S**(-1)*NS*MG2 + 4*S**(-1)*MS2**2*
     +    MG2**(-1) - 4*S**(-1)*MG2 - 16*UG*U1**(-2)*NS*MS2 + 8*UG*
     +    U1**(-2)*NS*MS2**2*MG2**(-1) + 8*UG*U1**(-2)*NS*MG2 + 16*UG*
     +    U1**(-2)*MS2 - 8*UG*U1**(-2)*MS2**2*MG2**(-1) - 8*UG*U1**(-2)
     +    *MG2 - 32*U1**(-2)*NS*MS2*MG2 + 24*U1**(-2)*NS*MS2**2 + 8*
     +    U1**(-2)*NS*MG2**2 + 32*U1**(-2)*MS2*MG2 - 24*U1**(-2)*MS2**2
     +     - 8*U1**(-2)*MG2**2 + 16*U1**(-1)*NS*MS2 - 16*U1**(-1)*NS*
     +    MG2 - 16*U1**(-1)*MS2 + 16*U1**(-1)*MG2 + 4*NS*MS2*MG2**(-1)
     +     + 4*NS )
     +
      M2QGV = M2QGV + SK1B0E(2)*CK*FOUR**(-1) * (  - 4*MS2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(2)*CQED*FOUR**(-1) * ( 2 - 64*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 72*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 32*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 + 4*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4 + 20*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**4 - 8*S**(-1)*TG**(-1)*MS2*
     +    MG2 - 4*S**(-1)*TG**(-1)*MS2**2 + 12*S**(-1)*TG**(-1)*MG2**2
     +     + 24*S**(-1)*UG*XLAM**(-2)*MS2*MG2 - 4*S**(-1)*UG*XLAM**(-2)
     +    *MS2**2 - 20*S**(-1)*UG*XLAM**(-2)*MG2**2 + 4*S**(-1)*UG + 16
     +    *S**(-1)*MG2 - 100*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2
     +     + 68*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2 - 12*S*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3 + 44*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 + 4*S*TG**(-1)*U1**(-1)*MS2 - 4*S*TG**(-1)*
     +    U1**(-1)*MG2 - 40*S*TG**(-1)*XLAM**(-2)*MS2*MG2 + 12*S*
     +    TG**(-1)*XLAM**(-2)*MS2**2 + 28*S*TG**(-1)*XLAM**(-2)*MG2**2
     +     - 8*S*UG**(-1)*T1**(-1)*MS2 + 8*S*UG**(-1)*T1**(-1)*MG2 - 6*
     +    S*UG*U1**(-1)*XLAM**(-2)*MS2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CQED*FOUR**(-1) * ( 6*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2 - 4*S*UG*XLAM**(-2) - 56*S*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2 + 8*S*U1**(-1)*XLAM**(-2)*MS2**2 - 2*S*U1**(-1) + 10
     +    *S*XLAM**(-2)*MS2 - 2*S*XLAM**(-2)*MG2 - 40*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 + 12*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2 + 28*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**2 - 4*S**2*TG**(-1)*XLAM**(-2)*MS2 + 4*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2 - 4*S**2*U1**(-1)*XLAM**(-2)*MS2 + 8*S**2*
     +    U1**(-1)*XLAM**(-2)*MG2 - 2*S**2*XLAM**(-2) - 4*S**3*TG**(-1)
     +    *U1**(-1)*XLAM**(-2)*MS2 + 4*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 + 2*S**3*U1**(-1)*XLAM**(-2) - 64*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**3 + 72*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 32*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MS2**3*MG2 + 4*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**4 + 20*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**4 - 8*TG**(-1)*U1**(-1)*MS2
     +    *MG2 )
     +
      M2QGV = M2QGV + SK1B0E(2)*CQED*FOUR**(-1) * (  - 4*TG**(-1)*
     +    U1**(-1)*MS2**2 + 12*TG**(-1)*U1**(-1)*MG2**2 - 100*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 + 68*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 
     +    12*TG**(-1)*XLAM**(-2)*MS2**3 + 44*TG**(-1)*XLAM**(-2)*MG2**3
     +     + 4*TG**(-1)*MS2 - 4*TG**(-1)*MG2 - 16*UG**(-1)*U1**(-1)*
     +    MG2**2 - 8*UG**(-1)*MS2 + 8*UG**(-1)*MG2 + 36*UG*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 4*UG*U1**(-1)*XLAM**(-2)*MS2**2 - 32*UG*
     +    U1**(-1)*XLAM**(-2)*MG2**2 + 16*UG*XLAM**(-2)*MG2 + 8*UG**2*
     +    U1**(-1)*XLAM**(-2)*MS2 - 8*UG**2*U1**(-1)*XLAM**(-2)*MG2 - 8
     +    *T1**(-1)*MS2 + 8*T1**(-1)*MG2 + 32*U1**(-2)*MS2*MG2 - 92*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 28*U1**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 - 4*U1**(-1)*XLAM**(-2)*MS2**3 + 68*U1**(-1)*
     +    XLAM**(-2)*MG2**3 + 4*U1**(-1)*MG2 - 6*XLAM**(-2)*MS2**2 + 6*
     +    XLAM**(-2)*MG2**2 )
     +
      M2QGV = M2QGV + SK1B0E(3)*CO*FOUR**(-1) * ( 20 - 32*S**(-1)*
     +    TG**(-2)*MS2*MT2*MG2 - 48*S**(-1)*TG**(-2)*MS2*MG2**2 + 16*
     +    S**(-1)*TG**(-2)*MS2**2*MT2 + 48*S**(-1)*TG**(-2)*MS2**2*MG2
     +     - 16*S**(-1)*TG**(-2)*MS2**3 + 16*S**(-1)*TG**(-2)*MT2*
     +    MG2**2 + 16*S**(-1)*TG**(-2)*MG2**3 - 32*S**(-1)*TG**(-1)*
     +    U1**(-1)*MS2*MT2*MG2 - 48*S**(-1)*TG**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 16*S**(-1)*TG**(-1)*U1**(-1)*MS2**2*MT2 + 48*S**(-1)
     +    *TG**(-1)*U1**(-1)*MS2**2*MG2 - 16*S**(-1)*TG**(-1)*U1**(-1)*
     +    MS2**3 + 16*S**(-1)*TG**(-1)*U1**(-1)*MT2*MG2**2 + 16*S**(-1)
     +    *TG**(-1)*U1**(-1)*MG2**3 - 32*S**(-1)*TG**(-1)*MS2*MT2 - 128
     +    *S**(-1)*TG**(-1)*MS2*MG2 + 64*S**(-1)*TG**(-1)*MS2**2 + 32*
     +    S**(-1)*TG**(-1)*MT2*MG2 + 64*S**(-1)*TG**(-1)*MG2**2 + 8*
     +    S**(-1)*UG*U1**(-1)*MS2*MT2*MG2**(-1) + 32*S**(-1)*UG*
     +    U1**(-1)*MS2 - 8*S**(-1)*UG*U1**(-1)*MT2 - 32*S**(-1)*UG*
     +    U1**(-1)*MG2 + 4*S**(-1)*UG*MS2*MG2**(-1) + 4*S**(-1)*UG*MT2*
     +    MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(3)*CO*FOUR**(-1) * (  - 12*S**(-1)*UG - 8*
     +    S**(-1)*UG**2*U1**(-1)*MS2*MG2**(-1) + 8*S**(-1)*UG**2*
     +    U1**(-1) - 32*S**(-1)*U1**(-1)*MS2*MT2 - 32*S**(-1)*U1**(-1)*
     +    MS2*MG2 + 16*S**(-1)*U1**(-1)*MS2**2 + 32*S**(-1)*U1**(-1)*
     +    MT2*MG2 + 16*S**(-1)*U1**(-1)*MG2**2 - 4*S**(-1)*MS2*MT2*
     +    MG2**(-1) - 48*S**(-1)*MS2 + 4*S**(-1)*MS2**2*MG2**(-1) + 4*
     +    S**(-1)*MT2 + 44*S**(-1)*MG2 - 16*S*TG**(-3)*MS2*MG2 + 16*S*
     +    TG**(-3)*MT2*MG2 - 16*S*TG**(-2)*U1**(-1)*MS2*MT2 - 64*S*
     +    TG**(-2)*U1**(-1)*MS2*MG2 + 32*S*TG**(-2)*U1**(-1)*MS2**2 + 
     +    16*S*TG**(-2)*U1**(-1)*MT2*MG2 + 32*S*TG**(-2)*U1**(-1)*
     +    MG2**2 - 16*S*TG**(-2)*MT2 + 16*S*TG**(-2)*MG2 - 16*S*
     +    TG**(-1)*U1**(-1)*MS2*MT2*MG2**(-1) - 104*S*TG**(-1)*U1**(-1)
     +    *MS2 + 32*S*TG**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 16*S*
     +    TG**(-1)*U1**(-1)*MT2 + 72*S*TG**(-1)*U1**(-1)*MG2 - 8*S*
     +    TG**(-1)*MS2*MG2**(-1) - 8*S*TG**(-1)*MT2*MG2**(-1) + 24*S*
     +    TG**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(3)*CO*FOUR**(-1) * (  - 16*S*U1**(-1)*MS2*
     +    MG2**(-1) + 16*S*U1**(-1) - 16*S**2*TG**(-2)*U1**(-1)*MS2 + 
     +    16*S**2*TG**(-2)*U1**(-1)*MG2 - 16*S**2*TG**(-1)*U1**(-1)*MS2
     +    *MG2**(-1) + 16*S**2*TG**(-1)*U1**(-1) - 32*TG**(-2)*U1**(-1)
     +    *MS2*MT2*MG2 - 48*TG**(-2)*U1**(-1)*MS2*MG2**2 + 16*TG**(-2)*
     +    U1**(-1)*MS2**2*MT2 + 48*TG**(-2)*U1**(-1)*MS2**2*MG2 - 16*
     +    TG**(-2)*U1**(-1)*MS2**3 + 16*TG**(-2)*U1**(-1)*MT2*MG2**2 + 
     +    16*TG**(-2)*U1**(-1)*MG2**3 - 16*TG**(-2)*MS2*MT2 - 72*
     +    TG**(-2)*MS2*MG2 + 32*TG**(-2)*MS2**2 + 8*TG**(-2)*MT2*MG2 + 
     +    48*TG**(-2)*MG2**2 - 72*TG**(-1)*U1**(-1)*MS2*MT2 - 176*
     +    TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*U1**(-1)*MS2**2*MT2*
     +    MG2**(-1) + 104*TG**(-1)*U1**(-1)*MS2**2 - 16*TG**(-1)*
     +    U1**(-1)*MS2**3*MG2**(-1) + 48*TG**(-1)*U1**(-1)*MT2*MG2 + 96
     +    *TG**(-1)*U1**(-1)*MG2**2 - 80*TG**(-1)*MS2 + 16*TG**(-1)*
     +    MS2**2*MG2**(-1) - 8*TG**(-1)*MT2 + 72*TG**(-1)*MG2 + 16*UG*
     +    U1**(-2)*MS2 )
     +
      M2QGV = M2QGV + SK1B0E(3)*CO*FOUR**(-1) * (  - 8*UG*U1**(-2)*
     +    MS2**2*MG2**(-1) - 8*UG*U1**(-2)*MG2 + 16*UG*U1**(-1)*MS2*
     +    MG2**(-1) - 16*UG*U1**(-1) - 8*U1**(-2)*MS2*MT2 + 32*U1**(-2)
     +    *MS2*MG2 + 8*U1**(-2)*MS2**2*MT2*MG2**(-1) - 24*U1**(-2)*
     +    MS2**2 - 8*U1**(-2)*MG2**2 - 8*U1**(-1)*MS2*MT2*MG2**(-1) - 
     +    72*U1**(-1)*MS2 + 16*U1**(-1)*MS2**2*MG2**(-1) + 8*U1**(-1)*
     +    MT2 + 56*U1**(-1)*MG2 - 12*MS2*MG2**(-1) - 4*MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1B0E(3)*CK*FOUR**(-1) * ( 4 - 8*S**(-1)*UG*
     +    U1**(-1)*MS2*MT2*MG2**(-1) + 16*S**(-1)*UG*U1**(-1)*MS2 + 8*
     +    S**(-1)*UG*U1**(-1)*MT2 - 16*S**(-1)*UG*U1**(-1)*MG2 - 4*
     +    S**(-1)*UG*MS2*MG2**(-1) - 4*S**(-1)*UG*MT2*MG2**(-1) + 12*
     +    S**(-1)*UG + 8*S**(-1)*UG**2*U1**(-1)*MS2*MG2**(-1) - 8*
     +    S**(-1)*UG**2*U1**(-1) + 4*S**(-1)*MS2*MT2*MG2**(-1) - 4*
     +    S**(-1)*MS2**2*MG2**(-1) - 4*S**(-1)*MT2 + 4*S**(-1)*MG2 - 16
     +    *UG*U1**(-2)*MS2 + 8*UG*U1**(-2)*MS2**2*MG2**(-1) + 8*UG*
     +    U1**(-2)*MG2 + 8*U1**(-2)*MS2*MT2 - 32*U1**(-2)*MS2*MG2 - 8*
     +    U1**(-2)*MS2**2*MT2*MG2**(-1) + 24*U1**(-2)*MS2**2 + 8*
     +    U1**(-2)*MG2**2 - 8*U1**(-1)*MS2*MT2*MG2**(-1) + 16*U1**(-1)*
     +    MS2 + 8*U1**(-1)*MT2 - 16*U1**(-1)*MG2 + 4*MS2*MG2**(-1) - 4*
     +    MT2*MG2**(-1) )
     +
      M2QGV = M2QGV + SK1BP(1)*N*CO*FOUR**(-1) * ( 16*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 16*S**(-1)*UG*U1**(-1)*MG2**2 - 8*S**(-1)*
     +    UG*MS2 + 16*S**(-1)*UG*MG2 + 16*S**(-1)*UG**2*U1**(-1)*MS2 - 
     +    16*S**(-1)*UG**2*U1**(-1)*MG2 + 8*S**(-1)*MS2*MG2 - 8*S**(-1)
     +    *MS2**2 + 64*S*TG**(-1)*U1**(-1)*MS2*MG2 - 64*S*TG**(-1)*
     +    U1**(-1)*MS2**2 + 16*S*TG**(-1)*MS2 - 32*S*TG**(-1)*MG2 + 32*
     +    S*U1**(-1)*MS2 - 32*S*U1**(-1)*MG2 + 32*S**2*TG**(-1)*
     +    U1**(-1)*MS2 - 32*S**2*TG**(-1)*U1**(-1)*MG2 - 32*TG**(-2)*
     +    MS2*MG2**2 + 32*TG**(-2)*MS2**2*MG2 - 32*TG**(-1)*U1**(-1)*
     +    MS2**2*MG2 + 32*TG**(-1)*U1**(-1)*MS2**3 + 32*TG**(-1)*MS2*
     +    MG2 - 32*TG**(-1)*MS2**2 - 32*UG*U1**(-2)*MS2*MG2 + 16*UG*
     +    U1**(-2)*MS2**2 + 16*UG*U1**(-2)*MG2**2 - 32*UG*U1**(-1)*MS2
     +     + 32*UG*U1**(-1)*MG2 - 48*U1**(-2)*MS2*MG2**2 + 32*U1**(-2)*
     +    MS2**2*MG2 + 16*U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*MG2 - 32*
     +    U1**(-1)*MS2**2 + 16*U1**(-1)*MG2**2 + 24*MS2 - 32*MG2 )
     +
      M2QGV = M2QGV + SK1BP(1)*N*CK*FOUR**(-1) * (  - 48*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 + 48*S**(-1)*UG*U1**(-1)*MG2**2 + 24*S**(-1)
     +    *UG*MS2 - 48*S**(-1)*UG*MG2 - 48*S**(-1)*UG**2*U1**(-1)*MS2
     +     + 48*S**(-1)*UG**2*U1**(-1)*MG2 - 24*S**(-1)*MS2*MG2 + 24*
     +    S**(-1)*MS2**2 - 64*S*TG**(-1)*U1**(-1)*MS2*MG2 + 64*S*
     +    TG**(-1)*U1**(-1)*MS2**2 - 16*S*TG**(-1)*MS2 + 32*S*TG**(-1)*
     +    MG2 - 32*S*U1**(-1)*MS2 + 32*S*U1**(-1)*MG2 - 32*S**2*
     +    TG**(-1)*U1**(-1)*MS2 + 32*S**2*TG**(-1)*U1**(-1)*MG2 + 32*
     +    TG**(-2)*MS2*MG2**2 - 32*TG**(-2)*MS2**2*MG2 + 32*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 32*TG**(-1)*U1**(-1)*MS2**3 - 32*
     +    TG**(-1)*MS2*MG2 + 32*TG**(-1)*MS2**2 + 96*UG*U1**(-2)*MS2*
     +    MG2 - 48*UG*U1**(-2)*MS2**2 - 48*UG*U1**(-2)*MG2**2 + 32*UG*
     +    U1**(-1)*MS2 - 32*UG*U1**(-1)*MG2 + 144*U1**(-2)*MS2*MG2**2
     +     - 96*U1**(-2)*MS2**2*MG2 - 48*U1**(-2)*MG2**3 - 48*U1**(-1)*
     +    MS2*MG2 + 32*U1**(-1)*MS2**2 + 16*U1**(-1)*MG2**2 - 40*MS2 + 
     +    32*MG2 )
     +
      M2QGV = M2QGV + SK1BP(1)*CQED*FOUR**(-1) * ( 16*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 16*S**(-1)*UG*U1**(-1)*MG2**2 - 8*S**(-1)*
     +    UG*MS2 + 16*S**(-1)*UG*MG2 + 16*S**(-1)*UG**2*U1**(-1)*MS2 - 
     +    16*S**(-1)*UG**2*U1**(-1)*MG2 + 8*S**(-1)*MS2*MG2 - 8*S**(-1)
     +    *MS2**2 - 32*UG*U1**(-2)*MS2*MG2 + 16*UG*U1**(-2)*MS2**2 + 16
     +    *UG*U1**(-2)*MG2**2 - 48*U1**(-2)*MS2*MG2**2 + 32*U1**(-2)*
     +    MS2**2*MG2 + 16*U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*MG2 - 16*
     +    U1**(-1)*MG2**2 + 8*MS2 )
     +
      M2QGV = M2QGV + SK1BP(2)*N*CO*FOUR**(-1) * ( 4*S**(-1)*UG*MS2 - 4
     +    *S**(-1)*UG*MG2 - 8*S**(-1)*UG**2*U1**(-1)*MS2 + 8*S**(-1)*
     +    UG**2*U1**(-1)*MG2 - 8*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2 + 4
     +    *S**(-1)*MG2**2 - 48*S*TG**(-1)*U1**(-1)*MS2*MG2 + 32*S*
     +    TG**(-1)*U1**(-1)*MS2**2 + 16*S*TG**(-1)*U1**(-1)*MG2**2 - 8*
     +    S*TG**(-1)*MS2 + 8*S*TG**(-1)*MG2 - 16*S*U1**(-1)*MS2 + 16*S*
     +    U1**(-1)*MG2 - 16*S**2*TG**(-1)*U1**(-1)*MS2 + 16*S**2*
     +    TG**(-1)*U1**(-1)*MG2 + 32*TG**(-2)*MS2*MG2**2 - 16*TG**(-2)*
     +    MS2**2*MG2 - 16*TG**(-2)*MG2**3 - 16*TG**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 32*TG**(-1)*U1**(-1)*MS2**2*MG2 - 16*TG**(-1)*
     +    U1**(-1)*MS2**3 - 16*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MS2**2 + 
     +    16*UG*U1**(-2)*MS2*MG2 - 8*UG*U1**(-2)*MS2**2 - 8*UG*U1**(-2)
     +    *MG2**2 + 16*UG*U1**(-1)*MS2 - 16*UG*U1**(-1)*MG2 + 16*
     +    U1**(-2)*MS2*MG2**2 - 8*U1**(-2)*MS2**2*MG2 - 8*U1**(-2)*
     +    MG2**3 - 16*U1**(-1)*MS2*MG2 + 16*U1**(-1)*MS2**2 - 12*MS2 + 
     +    12*MG2 )
     +
      M2QGV = M2QGV + SK1BP(2)*N*CK*FOUR**(-1) * (  - 12*S**(-1)*UG*MS2
     +     + 12*S**(-1)*UG*MG2 + 24*S**(-1)*UG**2*U1**(-1)*MS2 - 24*
     +    S**(-1)*UG**2*U1**(-1)*MG2 + 24*S**(-1)*MS2*MG2 - 12*S**(-1)*
     +    MS2**2 - 12*S**(-1)*MG2**2 + 48*S*TG**(-1)*U1**(-1)*MS2*MG2
     +     - 32*S*TG**(-1)*U1**(-1)*MS2**2 - 16*S*TG**(-1)*U1**(-1)*
     +    MG2**2 + 8*S*TG**(-1)*MS2 - 8*S*TG**(-1)*MG2 + 16*S*U1**(-1)*
     +    MS2 - 16*S*U1**(-1)*MG2 + 16*S**2*TG**(-1)*U1**(-1)*MS2 - 16*
     +    S**2*TG**(-1)*U1**(-1)*MG2 - 32*TG**(-2)*MS2*MG2**2 + 16*
     +    TG**(-2)*MS2**2*MG2 + 16*TG**(-2)*MG2**3 + 16*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 32*TG**(-1)*U1**(-1)*MS2**2*MG2 + 16*
     +    TG**(-1)*U1**(-1)*MS2**3 + 16*TG**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    MS2**2 - 48*UG*U1**(-2)*MS2*MG2 + 24*UG*U1**(-2)*MS2**2 + 24*
     +    UG*U1**(-2)*MG2**2 - 16*UG*U1**(-1)*MS2 + 16*UG*U1**(-1)*MG2
     +     - 48*U1**(-2)*MS2*MG2**2 + 24*U1**(-2)*MS2**2*MG2 + 24*
     +    U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*MG2 - 16*U1**(-1)*MS2**2 + 
     +    20*MS2 )
     +
      M2QGV = M2QGV + SK1BP(2)*N*CK*FOUR**(-1) * (  - 20*MG2 )
     +
      M2QGV = M2QGV + SK1BP(2)*CQED*FOUR**(-1) * ( 4*S**(-1)*UG*MS2 - 4
     +    *S**(-1)*UG*MG2 - 8*S**(-1)*UG**2*U1**(-1)*MS2 + 8*S**(-1)*
     +    UG**2*U1**(-1)*MG2 - 8*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2 + 4
     +    *S**(-1)*MG2**2 + 16*UG*U1**(-2)*MS2*MG2 - 8*UG*U1**(-2)*
     +    MS2**2 - 8*UG*U1**(-2)*MG2**2 + 16*U1**(-2)*MS2*MG2**2 - 8*
     +    U1**(-2)*MS2**2*MG2 - 8*U1**(-2)*MG2**3 - 4*MS2 + 4*MG2 )
     +
      M2QGV = M2QGV + SK1BP(3)*N*CO*FOUR**(-1) * (  - 2*S**(-1)*UG*MS2
     +     + 2*S**(-1)*UG*MG2 + 4*S**(-1)*UG**2*U1**(-1)*MS2 - 4*
     +    S**(-1)*UG**2*U1**(-1)*MG2 + 4*S**(-1)*MS2*MG2 - 2*S**(-1)*
     +    MS2**2 - 2*S**(-1)*MG2**2 + 24*S*TG**(-1)*U1**(-1)*MS2*MG2 - 
     +    16*S*TG**(-1)*U1**(-1)*MS2**2 - 8*S*TG**(-1)*U1**(-1)*MG2**2
     +     + 4*S*TG**(-1)*MS2 - 4*S*TG**(-1)*MG2 + 8*S*U1**(-1)*MS2 - 8
     +    *S*U1**(-1)*MG2 + 8*S**2*TG**(-1)*U1**(-1)*MS2 - 8*S**2*
     +    TG**(-1)*U1**(-1)*MG2 - 16*TG**(-2)*MS2*MG2**2 + 8*TG**(-2)*
     +    MS2**2*MG2 + 8*TG**(-2)*MG2**3 + 8*TG**(-1)*U1**(-1)*MS2*
     +    MG2**2 - 16*TG**(-1)*U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*
     +    U1**(-1)*MS2**3 + 8*TG**(-1)*MS2*MG2 - 8*TG**(-1)*MS2**2 - 8*
     +    UG*U1**(-2)*MS2*MG2 + 4*UG*U1**(-2)*MS2**2 + 4*UG*U1**(-2)*
     +    MG2**2 - 8*UG*U1**(-1)*MS2 + 8*UG*U1**(-1)*MG2 - 8*U1**(-2)*
     +    MS2*MG2**2 + 4*U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 + 8*
     +    U1**(-1)*MS2*MG2 - 8*U1**(-1)*MS2**2 + 6*MS2 - 6*MG2 )
     +
      M2QGV = M2QGV + SK1BP(3)*N*CK*FOUR**(-1) * ( 6*S**(-1)*UG*MS2 - 6
     +    *S**(-1)*UG*MG2 - 12*S**(-1)*UG**2*U1**(-1)*MS2 + 12*S**(-1)*
     +    UG**2*U1**(-1)*MG2 - 12*S**(-1)*MS2*MG2 + 6*S**(-1)*MS2**2 + 
     +    6*S**(-1)*MG2**2 - 24*S*TG**(-1)*U1**(-1)*MS2*MG2 + 16*S*
     +    TG**(-1)*U1**(-1)*MS2**2 + 8*S*TG**(-1)*U1**(-1)*MG2**2 - 4*S
     +    *TG**(-1)*MS2 + 4*S*TG**(-1)*MG2 - 8*S*U1**(-1)*MS2 + 8*S*
     +    U1**(-1)*MG2 - 8*S**2*TG**(-1)*U1**(-1)*MS2 + 8*S**2*TG**(-1)
     +    *U1**(-1)*MG2 + 16*TG**(-2)*MS2*MG2**2 - 8*TG**(-2)*MS2**2*
     +    MG2 - 8*TG**(-2)*MG2**3 - 8*TG**(-1)*U1**(-1)*MS2*MG2**2 + 16
     +    *TG**(-1)*U1**(-1)*MS2**2*MG2 - 8*TG**(-1)*U1**(-1)*MS2**3 - 
     +    8*TG**(-1)*MS2*MG2 + 8*TG**(-1)*MS2**2 + 24*UG*U1**(-2)*MS2*
     +    MG2 - 12*UG*U1**(-2)*MS2**2 - 12*UG*U1**(-2)*MG2**2 + 8*UG*
     +    U1**(-1)*MS2 - 8*UG*U1**(-1)*MG2 + 24*U1**(-2)*MS2*MG2**2 - 
     +    12*U1**(-2)*MS2**2*MG2 - 12*U1**(-2)*MG2**3 - 8*U1**(-1)*MS2*
     +    MG2 + 8*U1**(-1)*MS2**2 - 10*MS2 + 10*MG2 )
     +
      M2QGV = M2QGV + SK1BP(3)*CQED*FOUR**(-1) * (  - 2*S**(-1)*UG*MS2
     +     + 2*S**(-1)*UG*MG2 + 4*S**(-1)*UG**2*U1**(-1)*MS2 - 4*
     +    S**(-1)*UG**2*U1**(-1)*MG2 + 4*S**(-1)*MS2*MG2 - 2*S**(-1)*
     +    MS2**2 - 2*S**(-1)*MG2**2 - 8*UG*U1**(-2)*MS2*MG2 + 4*UG*
     +    U1**(-2)*MS2**2 + 4*UG*U1**(-2)*MG2**2 - 8*U1**(-2)*MS2*
     +    MG2**2 + 4*U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 + 2*MS2 - 
     +    2*MG2 )
     +
      M2QGV = M2QGV + SK1BP(5)*CO*FOUR**(-1) * (  - 16*S**(-1)*UG*
     +    U1**(-1)*MS2*MT2 + 16*S**(-1)*UG*U1**(-1)*MT2*MG2 - 8*S**(-1)
     +    *UG*MS2 - 8*S**(-1)*UG*MT2 + 8*S**(-1)*UG*MG2 + 16*S**(-1)*
     +    UG**2*U1**(-1)*MS2 - 16*S**(-1)*UG**2*U1**(-1)*MG2 + 8*
     +    S**(-1)*MS2*MT2 + 16*S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2 - 8*
     +    S**(-1)*MT2*MG2 - 8*S**(-1)*MG2**2 + 32*S*TG**(-1)*U1**(-1)*
     +    MS2*MT2 + 96*S*TG**(-1)*U1**(-1)*MS2*MG2 - 64*S*TG**(-1)*
     +    U1**(-1)*MS2**2 - 32*S*TG**(-1)*U1**(-1)*MT2*MG2 - 32*S*
     +    TG**(-1)*U1**(-1)*MG2**2 + 16*S*TG**(-1)*MS2 + 16*S*TG**(-1)*
     +    MT2 - 16*S*TG**(-1)*MG2 + 32*S*U1**(-1)*MS2 - 32*S*U1**(-1)*
     +    MG2 + 32*S**2*TG**(-1)*U1**(-1)*MS2 - 32*S**2*TG**(-1)*
     +    U1**(-1)*MG2 - 32*TG**(-2)*MS2*MT2*MG2 - 64*TG**(-2)*MS2*
     +    MG2**2 + 32*TG**(-2)*MS2**2*MG2 + 32*TG**(-2)*MT2*MG2**2 + 32
     +    *TG**(-2)*MG2**3 + 32*TG**(-1)*U1**(-1)*MS2*MT2*MG2 + 32*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 - 32*TG**(-1)*U1**(-1)*MS2**2*
     +    MT2 )
     +
      M2QGV = M2QGV + SK1BP(5)*CO*FOUR**(-1) * (  - 64*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 32*TG**(-1)*U1**(-1)*MS2**3 + 32*
     +    TG**(-1)*MS2*MG2 - 32*TG**(-1)*MS2**2 - 32*UG*U1**(-2)*MS2*
     +    MG2 + 16*UG*U1**(-2)*MS2**2 + 16*UG*U1**(-2)*MG2**2 - 32*UG*
     +    U1**(-1)*MS2 + 32*UG*U1**(-1)*MG2 + 16*U1**(-2)*MS2*MT2*MG2
     +     - 32*U1**(-2)*MS2*MG2**2 - 16*U1**(-2)*MS2**2*MT2 + 16*
     +    U1**(-2)*MS2**2*MG2 + 16*U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*
     +    MT2 + 32*U1**(-1)*MS2*MG2 - 32*U1**(-1)*MS2**2 - 16*U1**(-1)*
     +    MT2*MG2 + 24*MS2 + 8*MT2 - 24*MG2 )
     +
      M2QGV = M2QGV + SK1BP(5)*CK*FOUR**(-1) * ( 16*S**(-1)*UG*U1**(-1)
     +    *MS2*MT2 - 16*S**(-1)*UG*U1**(-1)*MT2*MG2 + 8*S**(-1)*UG*MS2
     +     + 8*S**(-1)*UG*MT2 - 8*S**(-1)*UG*MG2 - 16*S**(-1)*UG**2*
     +    U1**(-1)*MS2 + 16*S**(-1)*UG**2*U1**(-1)*MG2 - 8*S**(-1)*MS2*
     +    MT2 - 16*S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 + 8*S**(-1)*MT2*
     +    MG2 + 8*S**(-1)*MG2**2 + 32*UG*U1**(-2)*MS2*MG2 - 16*UG*
     +    U1**(-2)*MS2**2 - 16*UG*U1**(-2)*MG2**2 - 16*U1**(-2)*MS2*MT2
     +    *MG2 + 32*U1**(-2)*MS2*MG2**2 + 16*U1**(-2)*MS2**2*MT2 - 16*
     +    U1**(-2)*MS2**2*MG2 - 16*U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*
     +    MT2 - 16*U1**(-1)*MT2*MG2 - 8*MS2 + 8*MT2 + 8*MG2 )
     +
      M2QGV = M2QGV + SK1BP(6)*N*CO*FOUR**(-1) * ( 32*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 32*S**(-1)*UG*U1**(-1)*MG2**2 + 16*S**(-1)
     +    *UG*MG2 - 16*S**(-1)*MS2*MG2 + 16*S**(-1)*MG2**2 - 64*S*
     +    TG**(-1)*U1**(-1)*MS2*MG2 + 64*S*TG**(-1)*U1**(-1)*MG2**2 - 
     +    32*S*TG**(-1)*MG2 + 64*TG**(-2)*MS2*MG2**2 - 64*TG**(-2)*
     +    MG2**3 - 64*TG**(-1)*U1**(-1)*MS2*MG2**2 + 64*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 32*U1**(-2)*MS2*MG2**2 + 32*U1**(-2)*
     +    MS2**2*MG2 - 32*U1**(-1)*MS2*MG2 + 32*U1**(-1)*MG2**2 - 16*
     +    MG2 )
     +
      M2QGV = M2QGV + SK1BP(6)*N*CK*FOUR**(-1) * (  - 32*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 + 32*S**(-1)*UG*U1**(-1)*MG2**2 - 16*S**(-1)
     +    *UG*MG2 + 16*S**(-1)*MS2*MG2 - 16*S**(-1)*MG2**2 + 32*
     +    U1**(-2)*MS2*MG2**2 - 32*U1**(-2)*MS2**2*MG2 - 32*U1**(-1)*
     +    MS2*MG2 + 32*U1**(-1)*MG2**2 - 16*MG2 )
     +
      M2QGV = M2QGV + SK1BP(7)*CO*FOUR**(-1) * (  - 8*S**(-1)*UG*NS*MS2
     +     + 8*S**(-1)*UG*NS*MG2 + 8*S**(-1)*UG*MS2 - 8*S**(-1)*UG*MG2
     +     + 16*S**(-1)*UG**2*U1**(-1)*NS*MS2 - 16*S**(-1)*UG**2*
     +    U1**(-1)*NS*MG2 - 16*S**(-1)*UG**2*U1**(-1)*MS2 + 16*S**(-1)*
     +    UG**2*U1**(-1)*MG2 + 16*S**(-1)*NS*MS2*MG2 - 8*S**(-1)*NS*
     +    MS2**2 - 8*S**(-1)*NS*MG2**2 - 16*S**(-1)*MS2*MG2 + 8*S**(-1)
     +    *MS2**2 + 8*S**(-1)*MG2**2 + 96*S*TG**(-1)*U1**(-1)*NS*MS2*
     +    MG2 - 64*S*TG**(-1)*U1**(-1)*NS*MS2**2 - 32*S*TG**(-1)*
     +    U1**(-1)*NS*MG2**2 - 96*S*TG**(-1)*U1**(-1)*MS2*MG2 + 64*S*
     +    TG**(-1)*U1**(-1)*MS2**2 + 32*S*TG**(-1)*U1**(-1)*MG2**2 + 16
     +    *S*TG**(-1)*NS*MS2 - 16*S*TG**(-1)*NS*MG2 - 16*S*TG**(-1)*MS2
     +     + 16*S*TG**(-1)*MG2 + 32*S*U1**(-1)*NS*MS2 - 32*S*U1**(-1)*
     +    NS*MG2 - 32*S*U1**(-1)*MS2 + 32*S*U1**(-1)*MG2 + 32*S**2*
     +    TG**(-1)*U1**(-1)*NS*MS2 - 32*S**2*TG**(-1)*U1**(-1)*NS*MG2
     +     - 32*S**2*TG**(-1)*U1**(-1)*MS2 + 32*S**2*TG**(-1)*U1**(-1)*
     +    MG2 )
     +
      M2QGV = M2QGV + SK1BP(7)*CO*FOUR**(-1) * (  - 64*TG**(-2)*NS*MS2*
     +    MG2**2 + 32*TG**(-2)*NS*MS2**2*MG2 + 32*TG**(-2)*NS*MG2**3 + 
     +    64*TG**(-2)*MS2*MG2**2 - 32*TG**(-2)*MS2**2*MG2 - 32*TG**(-2)
     +    *MG2**3 + 32*TG**(-1)*U1**(-1)*NS*MS2*MG2**2 - 64*TG**(-1)*
     +    U1**(-1)*NS*MS2**2*MG2 + 32*TG**(-1)*U1**(-1)*NS*MS2**3 - 32*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 + 64*TG**(-1)*U1**(-1)*MS2**2*
     +    MG2 - 32*TG**(-1)*U1**(-1)*MS2**3 + 32*TG**(-1)*NS*MS2*MG2 - 
     +    32*TG**(-1)*NS*MS2**2 - 32*TG**(-1)*MS2*MG2 + 32*TG**(-1)*
     +    MS2**2 - 32*UG*U1**(-2)*NS*MS2*MG2 + 16*UG*U1**(-2)*NS*MS2**2
     +     + 16*UG*U1**(-2)*NS*MG2**2 + 32*UG*U1**(-2)*MS2*MG2 - 16*UG*
     +    U1**(-2)*MS2**2 - 16*UG*U1**(-2)*MG2**2 - 32*UG*U1**(-1)*NS*
     +    MS2 + 32*UG*U1**(-1)*NS*MG2 + 32*UG*U1**(-1)*MS2 - 32*UG*
     +    U1**(-1)*MG2 - 32*U1**(-2)*NS*MS2*MG2**2 + 16*U1**(-2)*NS*
     +    MS2**2*MG2 + 16*U1**(-2)*NS*MG2**3 + 32*U1**(-2)*MS2*MG2**2
     +     - 16*U1**(-2)*MS2**2*MG2 - 16*U1**(-2)*MG2**3 + 32*U1**(-1)*
     +    NS*MS2*MG2 )
     +
      M2QGV = M2QGV + SK1BP(7)*CO*FOUR**(-1) * (  - 32*U1**(-1)*NS*
     +    MS2**2 - 32*U1**(-1)*MS2*MG2 + 32*U1**(-1)*MS2**2 + 24*NS*MS2
     +     - 24*NS*MG2 - 24*MS2 + 24*MG2 )
     +
      M2QGV = M2QGV + SK1BP(7)*CK*FOUR**(-1) * ( 8*S**(-1)*UG*NS*MS2 - 
     +    8*S**(-1)*UG*NS*MG2 - 8*S**(-1)*UG*MS2 + 8*S**(-1)*UG*MG2 - 
     +    16*S**(-1)*UG**2*U1**(-1)*NS*MS2 + 16*S**(-1)*UG**2*U1**(-1)*
     +    NS*MG2 + 16*S**(-1)*UG**2*U1**(-1)*MS2 - 16*S**(-1)*UG**2*
     +    U1**(-1)*MG2 - 16*S**(-1)*NS*MS2*MG2 + 8*S**(-1)*NS*MS2**2 + 
     +    8*S**(-1)*NS*MG2**2 + 16*S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2
     +     - 8*S**(-1)*MG2**2 + 32*UG*U1**(-2)*NS*MS2*MG2 - 16*UG*
     +    U1**(-2)*NS*MS2**2 - 16*UG*U1**(-2)*NS*MG2**2 - 32*UG*
     +    U1**(-2)*MS2*MG2 + 16*UG*U1**(-2)*MS2**2 + 16*UG*U1**(-2)*
     +    MG2**2 + 32*U1**(-2)*NS*MS2*MG2**2 - 16*U1**(-2)*NS*MS2**2*
     +    MG2 - 16*U1**(-2)*NS*MG2**3 - 32*U1**(-2)*MS2*MG2**2 + 16*
     +    U1**(-2)*MS2**2*MG2 + 16*U1**(-2)*MG2**3 - 8*NS*MS2 + 8*NS*
     +    MG2 + 8*MS2 - 8*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(1)*N*CO*FOUR**(-1) * (  - 8*S*TG**(-1)*MS2
     +     + 8*S*TG**(-1)*MG2 + 4*S + 4*S**2*TG**(-1) - 16*TG**(-1)*MS2
     +    *MG2 + 8*TG**(-1)*MS2**2 + 8*TG**(-1)*MG2**2 - 4*UG - 4*MS2
     +     + 4*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(1)*N*CK*FOUR**(-1) * (  - 4*S*TG**(-1)*MG2
     +     + 8*S*U1**(-1)*MS2 - 4*S*U1**(-1)*MG2 - 16*TG**(-1)*MS2*MG2
     +     + 8*TG**(-1)*MS2**2 + 8*TG**(-1)*MG2**2 - 4*UG - 4*MS2 + 4*
     +    MG2 )
     +
      M2QGV = M2QGV + SK1C0A(1)*CQED*FOUR**(-1) * ( 8*S*U1**(-1)*MS2 - 
     +    12*S*U1**(-1)*MG2 + 8*S + 4*S**2*U1**(-1) + 8*UG*U1**(-1)*MS2
     +     - 8*UG*U1**(-1)*MG2 + 4*UG - 4*MS2 + 4*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(7)*N*CO*FOUR**(-1) * ( 8*S**(-1)*MS2*MG2
     +     - 8*S**(-1)*MG2**2 + 8*S*TG**(-1)*MS2 - 4*S**2*TG**(-1) - 8*
     +    TG**(-1)*MS2**2 + 4*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(7)*N*CK*FOUR**(-1) * (  - 8*S**(-1)*MS2*
     +    MG2 + 8*S**(-1)*MG2**2 - 4*S*U1**(-1)*MG2 - 24*U1**(-1)*MS2*
     +    MG2 + 16*U1**(-1)*MG2**2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(8)*N*CK*FOUR**(-1) * (  - 16*S**(-1)*MS2*
     +    MG2 + 16*S**(-1)*MS2**2 - 4*S*TG**(-1)*MG2 - 8*S*U1**(-1)*MS2
     +     + 8*S*U1**(-1)*MG2 - 4*S - 24*TG**(-1)*MS2*MG2 + 16*TG**(-1)
     +    *MG2**2 + 16*UG*U1**(-1)*MS2 - 16*UG*U1**(-1)*MG2 + 32*
     +    U1**(-1)*MS2*MG2 - 24*U1**(-1)*MG2**2 - 16*MS2 + 24*MG2 )
     +
      M2QGV = M2QGV + SK1C0A(8)*CQED*FOUR**(-1) * ( 8*S**(-1)*MS2*MG2
     +     - 8*S**(-1)*MS2**2 - 8*S*U1**(-1)*MG2 - 8*UG*U1**(-1)*MS2 + 
     +    8*UG*U1**(-1)*MG2 - 32*U1**(-1)*MS2*MG2 + 24*U1**(-1)*MG2**2
     +     + 8*MS2 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CO*FOUR**(-1) * ( 12*S**(-1)*TG**(-1)
     +    *XLAM**(-2)*MS2*MG2**4 - 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**2
     +    *MG2**3 - 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**3*MG2**2 + 12*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4*MG2 - 4*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**5 - 4*S**(-1)*TG**(-1)*XLAM**(-2)*MG2**5 + 4*
     +    S**(-1)*TG**(-1)*MS2*MG2**2 + 4*S**(-1)*TG**(-1)*MS2**2*MG2
     +     - 4*S**(-1)*TG**(-1)*MS2**3 - 4*S**(-1)*TG**(-1)*MG2**3 - 4*
     +    S**(-1)*UG*XLAM**(-2)*MS2*MG2**2 - 4*S**(-1)*UG*XLAM**(-2)*
     +    MS2**2*MG2 + 4*S**(-1)*UG*XLAM**(-2)*MS2**3 + 4*S**(-1)*UG*
     +    XLAM**(-2)*MG2**3 + 4*S**(-1)*UG*MS2 + 4*S**(-1)*UG*MG2 - 4*
     +    S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 - 4*S**(-1)*MG2**2 + 8*S*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 4*S*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 - 16*S*TG**(-1)*XLAM**(-2)*MS2**3 - 28*S*TG**(-1)*
     +    XLAM**(-2)*MG2**3 - 8*S*TG**(-1)*MS2 + 12*S*TG**(-1)*MG2 + 12
     +    *S*UG*XLAM**(-2)*MS2 + 12*S*UG*XLAM**(-2)*MG2 + 8*S + 12*S**2
     +    *TG**(-1)*XLAM**(-2)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CO*FOUR**(-1) * ( 16*S**2*TG**(-1)*
     +    XLAM**(-2)*MS2**2 + 28*S**2*TG**(-1)*XLAM**(-2)*MG2**2 - 4*
     +    S**2*UG*XLAM**(-2) - 12*S**3*TG**(-1)*XLAM**(-2)*MS2 - 16*
     +    S**3*TG**(-1)*XLAM**(-2)*MG2 + 4*S**4*TG**(-1)*XLAM**(-2) - 
     +    20*TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 4*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2**2 - 12*TG**(-1)*XLAM**(-2)*MS2**3*MG2 + 12*
     +    TG**(-1)*XLAM**(-2)*MS2**4 + 16*TG**(-1)*XLAM**(-2)*MG2**4 - 
     +    12*TG**(-1)*MS2*MG2 + 12*TG**(-1)*MS2**2 - 8*TG**(-1)*MG2**2
     +     - 8*UG*XLAM**(-2)*MS2*MG2 - 12*UG*XLAM**(-2)*MS2**2 - 12*UG*
     +    XLAM**(-2)*MG2**2 - 4*UG - 16*MS2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CK*FOUR**(-1) * (  - 12*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**4 + 8*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**3 + 8*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2**2 - 12*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4*MG2 + 4
     +    *S**(-1)*TG**(-1)*XLAM**(-2)*MS2**5 + 4*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MG2**5 - 4*S**(-1)*TG**(-1)*MS2*MG2**2 - 4*S**(-1)
     +    *TG**(-1)*MS2**2*MG2 + 4*S**(-1)*TG**(-1)*MS2**3 + 4*S**(-1)*
     +    TG**(-1)*MG2**3 - 2*S**(-1)*UG*XLAM**(-2)*MS2*MG2**2 + 10*
     +    S**(-1)*UG*XLAM**(-2)*MS2**2*MG2 - 6*S**(-1)*UG*XLAM**(-2)*
     +    MS2**3 - 2*S**(-1)*UG*XLAM**(-2)*MG2**3 - 2*S**(-1)*UG*MS2 - 
     +    6*S**(-1)*UG*MG2 + 4*S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2 + 4*
     +    S**(-1)*MG2**2 - 8*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3
     +     - 24*S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 40*S*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**3*MG2 - 16*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**4 + 8*S*TG**(-1)*U1**(-1)*XLAM**(-2)
     +    *MG2**4 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CK*FOUR**(-1) * ( 8*S*TG**(-1)*
     +    U1**(-1)*MS2*MG2 - 8*S*TG**(-1)*U1**(-1)*MS2**2 + 24*S*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 48*S*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 + 24*S*TG**(-1)*XLAM**(-2)*MS2**3 + 4*S*TG**(-1)*
     +    MS2 - 4*S*TG**(-1)*MG2 - 20*S*UG*U1**(-1)*XLAM**(-2)*MS2*MG2
     +     + 8*S*UG*U1**(-1)*XLAM**(-2)*MS2**2 + 12*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2**2 - 10*S*UG*XLAM**(-2)*MS2 - 6*S*UG*
     +    XLAM**(-2)*MG2 + 4*S*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 20*S*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 12*S*U1**(-1)*XLAM**(-2)*
     +    MS2**3 - 28*S*U1**(-1)*XLAM**(-2)*MG2**3 + 12*S*U1**(-1)*MS2
     +     + 28*S*XLAM**(-2)*MS2*MG2 + 4*S*XLAM**(-2)*MS2**2 + 24*S**2*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 48*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 24*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**3 + 4*S**2*TG**(-1)*U1**(-1)*MS2 - 4*S**2*
     +    TG**(-1)*U1**(-1)*MG2 + 24*S**2*TG**(-1)*XLAM**(-2)*MS2*MG2
     +     - 16*S**2*TG**(-1)*XLAM**(-2)*MS2**2 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CK*FOUR**(-1) * (  - 8*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2**2 + 2*S**2*UG*XLAM**(-2) + 36*S**2*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 12*S**2*U1**(-1)*XLAM**(-2)*MS2**2 + 32*
     +    S**2*U1**(-1)*XLAM**(-2)*MG2**2 - 4*S**2*U1**(-1) - 12*S**2*
     +    XLAM**(-2)*MS2 - 12*S**2*XLAM**(-2)*MG2 + 24*S**3*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 - 16*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2 - 8*S**3*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 4*S**3*TG**(-1)*XLAM**(-2)*MS2 - 4*S**3*TG**(-1)*
     +    XLAM**(-2)*MG2 - 4*S**3*U1**(-1)*XLAM**(-2)*MS2 - 24*S**3*
     +    U1**(-1)*XLAM**(-2)*MG2 + 4*S**3*XLAM**(-2) + 4*S**4*TG**(-1)
     +    *U1**(-1)*XLAM**(-2)*MS2 - 4*S**4*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2 + 4*S**4*U1**(-1)*XLAM**(-2) - 12*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**4 + 8*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**3 + 8*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MS2**3*MG2**2 - 12*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**4*MG2 + 
     +    4*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**5 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CK*FOUR**(-1) * ( 4*TG**(-1)*U1**(-1)
     +    *XLAM**(-2)*MG2**5 - 4*TG**(-1)*U1**(-1)*MS2*MG2**2 - 4*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 + 4*TG**(-1)*U1**(-1)*MS2**3 + 4
     +    *TG**(-1)*U1**(-1)*MG2**3 - 8*TG**(-1)*XLAM**(-2)*MS2*MG2**3
     +     - 24*TG**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 40*TG**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 - 16*TG**(-1)*XLAM**(-2)*MS2**4 + 8*
     +    TG**(-1)*XLAM**(-2)*MG2**4 + 8*TG**(-1)*MS2*MG2 - 8*TG**(-1)*
     +    MS2**2 + 44*UG*U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 4*UG*U1**(-1)
     +    *XLAM**(-2)*MS2**2*MG2 - 4*UG*U1**(-1)*XLAM**(-2)*MS2**3 - 44
     +    *UG*U1**(-1)*XLAM**(-2)*MG2**3 - 4*UG*U1**(-1)*MS2 + 4*UG*
     +    U1**(-1)*MG2 - 8*UG*XLAM**(-2)*MS2*MG2 + 14*UG*XLAM**(-2)*
     +    MS2**2 + 26*UG*XLAM**(-2)*MG2**2 + 6*UG + 28*UG**2*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 4*UG**2*U1**(-1)*XLAM**(-2)*MS2**2 - 24*
     +    UG**2*U1**(-1)*XLAM**(-2)*MG2**2 - 4*UG**2*XLAM**(-2)*MS2 + 4
     +    *UG**2*XLAM**(-2)*MG2 + 4*UG**3*U1**(-1)*XLAM**(-2)*MS2 - 4*
     +    UG**3*U1**(-1)*XLAM**(-2)*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(1)*N*CK*FOUR**(-1) * (  - 8*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2**3 + 8*U1**(-1)*XLAM**(-2)*MS2**3*MG2 - 4*
     +    U1**(-1)*XLAM**(-2)*MS2**4 + 4*U1**(-1)*XLAM**(-2)*MG2**4 - 
     +    12*U1**(-1)*MS2*MG2 - 4*U1**(-1)*MS2**2 + 8*U1**(-1)*MG2**2
     +     - 20*XLAM**(-2)*MS2**2*MG2 + 20*XLAM**(-2)*MG2**3 + 12*MS2
     +     - 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(2)*N*CK*FOUR**(-1) * ( 64*S**(-1)*TG**(-1)
     +    *XLAM**(-2)*MS2*MG2**4 - 96*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2**3 + 64*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**3*MG2**2
     +     - 16*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4*MG2 - 16*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MG2**5 + 32*S**(-1)*TG**(-1)*MS2*MG2**2
     +     - 16*S**(-1)*TG**(-1)*MS2**2*MG2 - 16*S**(-1)*TG**(-1)*
     +    MG2**3 - 32*S**(-1)*UG*XLAM**(-2)*MS2*MG2**2 + 16*S**(-1)*UG*
     +    XLAM**(-2)*MS2**2*MG2 + 16*S**(-1)*UG*XLAM**(-2)*MG2**3 + 16*
     +    S**(-1)*UG*MG2 + 24*S**(-1)*MS2*MG2 - 24*S**(-1)*MG2**2 + 72*
     +    S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3 - 72*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 + 24*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 - 24*S*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**4 + 8*S*TG**(-1)*U1**(-1)*MS2*MG2 - 8*S*TG**(-1)*
     +    U1**(-1)*MG2**2 + 64*S*TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 8*S*
     +    TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 24*S*TG**(-1)*XLAM**(-2)*
     +    MG2**3 )
     +
      M2QGV = M2QGV + SK2C0B(2)*N*CK*FOUR**(-1) * (  - 12*S*TG**(-1)*
     +    MG2 - 8*S*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 8*S*UG*U1**(-1)*
     +    XLAM**(-2)*MG2**2 - 16*S*UG*XLAM**(-2)*MG2 + 16*S*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2**2 - 16*S*U1**(-1)*XLAM**(-2)*MS2**2*MG2
     +     - 32*S*U1**(-1)*XLAM**(-2)*MG2**3 + 12*S*U1**(-1)*MG2 + 24*S
     +    *XLAM**(-2)*MS2*MG2 + 8*S*XLAM**(-2)*MG2**2 + 48*S**2*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 24*S**2*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2 - 24*S**2*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**3 - 16*S**2*TG**(-1)*XLAM**(-2)*MS2*MG2 - 24*
     +    S**2*TG**(-1)*XLAM**(-2)*MG2**2 + 24*S**2*U1**(-1)*XLAM**(-2)
     +    *MS2*MG2 + 16*S**2*U1**(-1)*XLAM**(-2)*MG2**2 - 8*S**2*
     +    XLAM**(-2)*MG2 + 8*S**3*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2
     +     - 8*S**3*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**2 + 8*S**3*
     +    TG**(-1)*XLAM**(-2)*MG2 - 8*S**3*U1**(-1)*XLAM**(-2)*MG2 + 32
     +    *TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**4 - 48*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2**3 )
     +
      M2QGV = M2QGV + SK2C0B(2)*N*CK*FOUR**(-1) * ( 32*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**3*MG2**2 - 8*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**4*MG2 - 8*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**5
     +     + 16*TG**(-1)*U1**(-1)*MS2*MG2**2 - 8*TG**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 8*TG**(-1)*U1**(-1)*MG2**3 + 80*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2**3 - 88*TG**(-1)*XLAM**(-2)*MS2**2*MG2**2
     +     + 32*TG**(-1)*XLAM**(-2)*MS2**3*MG2 - 24*TG**(-1)*XLAM**(-2)
     +    *MG2**4 + 20*TG**(-1)*MS2*MG2 + 4*TG**(-1)*MG2**2 - 24*UG*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 8*UG*U1**(-1)*XLAM**(-2)*
     +    MS2**2*MG2 + 16*UG*U1**(-1)*XLAM**(-2)*MG2**3 + 24*UG*
     +    XLAM**(-2)*MS2*MG2 + 40*UG*XLAM**(-2)*MG2**2 + 8*UG**2*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*UG**2*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 24*U1**(-1)*XLAM**(-2)*MS2*MG2**3 - 24*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 + 8*U1**(-1)*XLAM**(-2)*MS2**3*MG2
     +     - 8*U1**(-1)*XLAM**(-2)*MG2**4 + 12*U1**(-1)*MS2*MG2 - 36*
     +    U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0B(2)*N*CK*FOUR**(-1) * ( 24*XLAM**(-2)*MS2*
     +    MG2**2 - 8*XLAM**(-2)*MS2**2*MG2 - 16*XLAM**(-2)*MG2**3 + 8*
     +    MG2 )
     +
      M2QGV = M2QGV + SK2C0B(2)*CQED*FOUR**(-1) * (  - 32*S**(-1)*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**4 + 48*S**(-1)*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**3 - 32*S**(-1)*TG**(-1)*XLAM**(-2)*
     +    MS2**3*MG2**2 + 8*S**(-1)*TG**(-1)*XLAM**(-2)*MS2**4*MG2 + 8*
     +    S**(-1)*TG**(-1)*XLAM**(-2)*MG2**5 - 16*S**(-1)*TG**(-1)*MS2*
     +    MG2**2 + 8*S**(-1)*TG**(-1)*MS2**2*MG2 + 8*S**(-1)*TG**(-1)*
     +    MG2**3 + 16*S**(-1)*UG*XLAM**(-2)*MS2*MG2**2 - 8*S**(-1)*UG*
     +    XLAM**(-2)*MS2**2*MG2 - 8*S**(-1)*UG*XLAM**(-2)*MG2**3 - 8*
     +    S**(-1)*UG*MG2 - 12*S**(-1)*MS2*MG2 + 12*S**(-1)*MG2**2 - 72*
     +    S*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**3 + 72*S*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2**2 - 24*S*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 + 24*S*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**4 - 8*S*TG**(-1)*U1**(-1)*MS2*MG2 + 8*S*TG**(-1)*
     +    U1**(-1)*MG2**2 - 48*S*TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 24*S*
     +    TG**(-1)*XLAM**(-2)*MS2**2*MG2 + 24*S*TG**(-1)*XLAM**(-2)*
     +    MG2**3 )
     +
      M2QGV = M2QGV + SK2C0B(2)*CQED*FOUR**(-1) * ( 8*S*UG*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 8*S*UG*U1**(-1)*XLAM**(-2)*MG2**2 + 8*S*
     +    UG*XLAM**(-2)*MG2 - 16*S*U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 16*
     +    S*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 32*S*U1**(-1)*XLAM**(-2)*
     +    MG2**3 - 16*S*U1**(-1)*MG2 - 24*S*XLAM**(-2)*MS2*MG2 - 8*S*
     +    XLAM**(-2)*MG2**2 - 48*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*
     +    MG2**2 + 24*S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MS2**2*MG2 + 24
     +    *S**2*TG**(-1)*U1**(-1)*XLAM**(-2)*MG2**3 - 8*S**2*TG**(-1)*
     +    XLAM**(-2)*MS2*MG2 + 8*S**2*TG**(-1)*XLAM**(-2)*MG2**2 - 24*
     +    S**2*U1**(-1)*XLAM**(-2)*MS2*MG2 - 16*S**2*U1**(-1)*
     +    XLAM**(-2)*MG2**2 + 8*S**2*XLAM**(-2)*MG2 - 8*S**3*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 + 8*S**3*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MG2**2 + 8*S**3*U1**(-1)*XLAM**(-2)*MG2 - 32*
     +    TG**(-1)*U1**(-1)*XLAM**(-2)*MS2*MG2**4 + 48*TG**(-1)*
     +    U1**(-1)*XLAM**(-2)*MS2**2*MG2**3 - 32*TG**(-1)*U1**(-1)*
     +    XLAM**(-2)*MS2**3*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0B(2)*CQED*FOUR**(-1) * ( 8*TG**(-1)*U1**(-1)
     +    *XLAM**(-2)*MS2**4*MG2 + 8*TG**(-1)*U1**(-1)*XLAM**(-2)*
     +    MG2**5 - 16*TG**(-1)*U1**(-1)*MS2*MG2**2 + 8*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*U1**(-1)*MG2**3 - 72*
     +    TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 72*TG**(-1)*XLAM**(-2)*
     +    MS2**2*MG2**2 - 24*TG**(-1)*XLAM**(-2)*MS2**3*MG2 + 24*
     +    TG**(-1)*XLAM**(-2)*MG2**4 - 8*TG**(-1)*MS2*MG2 + 8*TG**(-1)*
     +    MG2**2 + 24*UG*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 8*UG*U1**(-1)
     +    *XLAM**(-2)*MS2**2*MG2 - 16*UG*U1**(-1)*XLAM**(-2)*MG2**3 - 8
     +    *UG*XLAM**(-2)*MS2*MG2 - 24*UG*XLAM**(-2)*MG2**2 - 8*UG**2*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 + 8*UG**2*U1**(-1)*XLAM**(-2)*
     +    MG2**2 - 24*U1**(-1)*XLAM**(-2)*MS2*MG2**3 + 24*U1**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 8*U1**(-1)*XLAM**(-2)*MS2**3*MG2
     +     + 8*U1**(-1)*XLAM**(-2)*MG2**4 - 16*U1**(-1)*MS2*MG2 + 40*
     +    U1**(-1)*MG2**2 - 24*XLAM**(-2)*MS2*MG2**2 + 8*XLAM**(-2)*
     +    MS2**2*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(2)*CQED*FOUR**(-1) * ( 16*XLAM**(-2)*
     +    MG2**3 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(3)*N*CK*FOUR**(-1) * ( 48*S**(-1)*TG**(-1)
     +    *MS2*MG2**2 - 48*S**(-1)*TG**(-1)*MS2**2*MG2 + 16*S**(-1)*
     +    TG**(-1)*MS2**3 - 16*S**(-1)*TG**(-1)*MG2**3 - 8*S**(-1)*UG*
     +    MS2 + 8*S**(-1)*UG*MG2 + 16*S**(-1)*MS2*MG2 - 8*S**(-1)*
     +    MS2**2 - 8*S**(-1)*MG2**2 + 32*S*TG**(-1)*U1**(-1)*MS2*MG2 - 
     +    16*S*TG**(-1)*U1**(-1)*MS2**2 - 16*S*TG**(-1)*U1**(-1)*MG2**2
     +     - 26*S*TG**(-1)*XLAM**(-2)*MS2*MG2**2 - 14*S*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2 + 10*S*TG**(-1)*XLAM**(-2)*MS2**3 + 30*
     +    S*TG**(-1)*XLAM**(-2)*MG2**3 + 6*S*TG**(-1)*MS2 - 18*S*
     +    TG**(-1)*MG2 - 28*S*UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 28*S*UG*
     +    U1**(-1)*XLAM**(-2)*MG2**2 - 4*S*UG*XLAM**(-2)*MS2 - 48*S*UG*
     +    XLAM**(-2)*MG2 - 14*S*UG**2*U1**(-1)*XLAM**(-2)*MS2 + 14*S*
     +    UG**2*U1**(-1)*XLAM**(-2)*MG2 + 24*S*U1**(-1)*XLAM**(-2)*MS2*
     +    MG2**2 - 24*S*U1**(-1)*XLAM**(-2)*MG2**3 - 2*S*U1**(-1)*MS2
     +     + 14*S*U1**(-1)*MG2 - 16*S*XLAM**(-2)*MS2*MG2 + 14*S*
     +    XLAM**(-2)*MS2**2 )
     +
      M2QGV = M2QGV + SK2C0B(3)*N*CK*FOUR**(-1) * ( 2*S*XLAM**(-2)*
     +    MG2**2 - 20*S + 8*S**2*TG**(-1)*U1**(-1)*MS2 - 8*S**2*
     +    TG**(-1)*U1**(-1)*MG2 - 4*S**2*TG**(-1)*XLAM**(-2)*MS2*MG2 - 
     +    10*S**2*TG**(-1)*XLAM**(-2)*MS2**2 - 26*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2**2 + 2*S**2*TG**(-1) + 12*S**2*UG*U1**(-1)*
     +    XLAM**(-2)*MS2 - 12*S**2*UG*U1**(-1)*XLAM**(-2)*MG2 + 17*S**2
     +    *UG*XLAM**(-2) + 24*S**2*U1**(-1)*XLAM**(-2)*MS2*MG2 + 32*
     +    S**2*U1**(-1)*XLAM**(-2)*MG2**2 - 4*S**2*U1**(-1) - 26*S**2*
     +    XLAM**(-2)*MS2 - 14*S**2*XLAM**(-2)*MG2 + 6*S**3*TG**(-1)*
     +    XLAM**(-2)*MS2 + 10*S**3*TG**(-1)*XLAM**(-2)*MG2 - 6*S**3*
     +    U1**(-1)*XLAM**(-2)*MS2 - 22*S**3*U1**(-1)*XLAM**(-2)*MG2 + 
     +    12*S**3*XLAM**(-2) - 2*S**4*TG**(-1)*XLAM**(-2) + 4*S**4*
     +    U1**(-1)*XLAM**(-2) + 24*TG**(-1)*U1**(-1)*MS2*MG2**2 - 24*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*U1**(-1)*MS2**3 - 8
     +    *TG**(-1)*U1**(-1)*MG2**3 + 40*TG**(-1)*XLAM**(-2)*MS2*MG2**3
     +     - 48*TG**(-1)*XLAM**(-2)*MS2**2*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0B(3)*N*CK*FOUR**(-1) * ( 24*TG**(-1)*
     +    XLAM**(-2)*MS2**3*MG2 - 4*TG**(-1)*XLAM**(-2)*MS2**4 - 12*
     +    TG**(-1)*XLAM**(-2)*MG2**4 + 12*TG**(-1)*MS2*MG2 - 20*
     +    TG**(-1)*MS2**2 + 8*TG**(-1)*MG2**2 - 12*UG*U1**(-1)*MS2 + 12
     +    *UG*U1**(-1)*MG2 - 2*UG*XLAM**(-2)*MS2*MG2 - 3*UG*XLAM**(-2)*
     +    MS2**2 + 5*UG*XLAM**(-2)*MG2**2 - 25*UG - 8*UG**2*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 + 8*UG**2*U1**(-1)*XLAM**(-2)*MG2**2 - 4*
     +    UG**2*XLAM**(-2)*MS2 + 4*UG**2*XLAM**(-2)*MG2 + 4*UG**3*
     +    U1**(-1)*XLAM**(-2)*MS2 - 4*UG**3*U1**(-1)*XLAM**(-2)*MG2 + 
     +    40*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MS2**2 - 32*U1**(-1)*MG2**2
     +     - 22*XLAM**(-2)*MS2*MG2**2 + 14*XLAM**(-2)*MS2**2*MG2 - 2*
     +    XLAM**(-2)*MS2**3 + 10*XLAM**(-2)*MG2**3 + 2*MS2 + 6*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(3)*CQED*FOUR**(-1) * (  - 24*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*TG**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*MG2**3 + 4*
     +    S**(-1)*UG*MS2 - 4*S**(-1)*UG*MG2 - 8*S**(-1)*MS2*MG2 + 4*
     +    S**(-1)*MS2**2 + 4*S**(-1)*MG2**2 - 32*S*TG**(-1)*U1**(-1)*
     +    MS2*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2**2 + 16*S*TG**(-1)*
     +    U1**(-1)*MG2**2 - 8*S*TG**(-1)*MS2 + 8*S*TG**(-1)*MG2 + 24*S*
     +    UG*U1**(-1)*XLAM**(-2)*MS2*MG2 - 24*S*UG*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 24*S*UG*XLAM**(-2)*MG2 + 8*S*UG**2*U1**(-1)*
     +    XLAM**(-2)*MS2 - 8*S*UG**2*U1**(-1)*XLAM**(-2)*MG2 - 24*S*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2**2 + 24*S*U1**(-1)*XLAM**(-2)*
     +    MG2**3 + 4*S*U1**(-1)*MS2 - 24*S*U1**(-1)*MG2 + 4*S*
     +    XLAM**(-2)*MS2*MG2 - 6*S*XLAM**(-2)*MS2**2 + 2*S*XLAM**(-2)*
     +    MG2**2 + 18*S - 8*S**2*TG**(-1)*U1**(-1)*MS2 + 8*S**2*
     +    TG**(-1)*U1**(-1)*MG2 - 8*S**2*UG*XLAM**(-2) + 8*S**2*
     +    U1**(-1)*XLAM**(-2)*MS2*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(3)*CQED*FOUR**(-1) * (  - 48*S**2*U1**(-1)
     +    *XLAM**(-2)*MG2**2 + 8*S**2*U1**(-1) + 10*S**2*XLAM**(-2)*MS2
     +     + 22*S**2*XLAM**(-2)*MG2 + 24*S**3*U1**(-1)*XLAM**(-2)*MG2
     +     - 10*S**3*XLAM**(-2) - 4*S**4*U1**(-1)*XLAM**(-2) - 24*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 + 24*TG**(-1)*U1**(-1)*MS2**2*
     +    MG2 - 8*TG**(-1)*U1**(-1)*MS2**3 + 8*TG**(-1)*U1**(-1)*MG2**3
     +     - 32*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MS2**2 + 16*TG**(-1)*
     +    MG2**2 + 12*UG*U1**(-1)*MS2 - 12*UG*U1**(-1)*MG2 - 8*UG*
     +    XLAM**(-2)*MS2*MG2 + 4*UG*XLAM**(-2)*MS2**2 + 4*UG*XLAM**(-2)
     +    *MG2**2 + 12*UG + 8*UG**2*U1**(-1)*XLAM**(-2)*MS2*MG2 - 8*
     +    UG**2*U1**(-1)*XLAM**(-2)*MG2**2 + 4*UG**2*XLAM**(-2)*MS2 - 4
     +    *UG**2*XLAM**(-2)*MG2 - 4*UG**3*U1**(-1)*XLAM**(-2)*MS2 + 4*
     +    UG**3*U1**(-1)*XLAM**(-2)*MG2 - 44*U1**(-1)*MS2*MG2 + 8*
     +    U1**(-1)*MS2**2 + 36*U1**(-1)*MG2**2 + 22*XLAM**(-2)*MS2*
     +    MG2**2 - 14*XLAM**(-2)*MS2**2*MG2 + 2*XLAM**(-2)*MS2**3 - 10*
     +    XLAM**(-2)*MG2**3 )
     +
      M2QGV = M2QGV + SK2C0B(3)*CQED*FOUR**(-1) * (  - 6*MS2 - 2*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(4)*N*CO*FOUR**(-1) * ( 24*S**(-1)*TG**(-1)
     +    *MS2*MG2**2 - 24*S**(-1)*TG**(-1)*MS2**2*MG2 + 8*S**(-1)*
     +    TG**(-1)*MS2**3 - 8*S**(-1)*TG**(-1)*MG2**3 - 4*S**(-1)*UG*
     +    MS2 + 4*S**(-1)*UG*MG2 + 8*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2
     +     - 4*S**(-1)*MG2**2 + 50*S*TG**(-1)*XLAM**(-2)*MS2*MG2**2 + 6
     +    *S*TG**(-1)*XLAM**(-2)*MS2**2*MG2 - 18*S*TG**(-1)*XLAM**(-2)*
     +    MS2**3 - 38*S*TG**(-1)*XLAM**(-2)*MG2**3 + 14*S*TG**(-1)*MS2
     +     + 2*S*TG**(-1)*MG2 + 4*S*UG*XLAM**(-2)*MS2 + 12*S*UG*
     +    XLAM**(-2)*MG2 - 4*S + 16*S**2*TG**(-1)*XLAM**(-2)*MS2*MG2 + 
     +    12*S**2*TG**(-1)*XLAM**(-2)*MS2**2 + 28*S**2*TG**(-1)*
     +    XLAM**(-2)*MG2**2 - 4*S**2*TG**(-1) - 2*S**2*UG*XLAM**(-2) - 
     +    2*S**3*TG**(-1)*XLAM**(-2)*MS2 - 14*S**3*TG**(-1)*XLAM**(-2)*
     +    MG2 - 32*TG**(-1)*XLAM**(-2)*MS2*MG2**3 + 48*TG**(-1)*
     +    XLAM**(-2)*MS2**2*MG2**2 - 32*TG**(-1)*XLAM**(-2)*MS2**3*MG2
     +     + 8*TG**(-1)*XLAM**(-2)*MS2**4 + 8*TG**(-1)*XLAM**(-2)*
     +    MG2**4 )
     +
      M2QGV = M2QGV + SK2C0B(4)*N*CO*FOUR**(-1) * ( 32*TG**(-1)*MS2*MG2
     +     - 24*TG**(-1)*MS2**2 - 8*TG**(-1)*MG2**2 - 4*UG*XLAM**(-2)*
     +    MS2*MG2 - 2*UG*XLAM**(-2)*MS2**2 + 6*UG*XLAM**(-2)*MG2**2 + 6
     +    *UG + 8*MS2 - 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0B(4)*N*CK*FOUR**(-1) * (  - 24*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*TG**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*MG2**3 + 4*
     +    S**(-1)*UG*MS2 - 4*S**(-1)*UG*MG2 - 8*S**(-1)*MS2*MG2 + 4*
     +    S**(-1)*MS2**2 + 4*S**(-1)*MG2**2 - 32*S*TG**(-1)*U1**(-1)*
     +    MS2*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2**2 + 16*S*TG**(-1)*
     +    U1**(-1)*MG2**2 - 8*S*TG**(-1)*MS2 + 8*S*TG**(-1)*MG2 - 40*S*
     +    UG*U1**(-1)*XLAM**(-2)*MS2*MG2 + 40*S*UG*U1**(-1)*XLAM**(-2)*
     +    MG2**2 + 16*S*UG*XLAM**(-2)*MS2 - 40*S*UG*XLAM**(-2)*MG2 - 24
     +    *S*UG**2*U1**(-1)*XLAM**(-2)*MS2 + 24*S*UG**2*U1**(-1)*
     +    XLAM**(-2)*MG2 + 8*S*U1**(-1)*XLAM**(-2)*MS2*MG2**2 - 8*S*
     +    U1**(-1)*XLAM**(-2)*MG2**3 + 12*S*U1**(-1)*MG2 + 24*S*
     +    XLAM**(-2)*MS2**2 - 24*S*XLAM**(-2)*MG2**2 - 8*S**2*TG**(-1)*
     +    U1**(-1)*MS2 + 8*S**2*TG**(-1)*U1**(-1)*MG2 + 24*S**2*UG*
     +    U1**(-1)*XLAM**(-2)*MS2 - 24*S**2*UG*U1**(-1)*XLAM**(-2)*MG2
     +     + 4*S**2*UG*XLAM**(-2) )
     +
      M2QGV = M2QGV + SK2C0B(4)*N*CK*FOUR**(-1) * ( 56*S**2*U1**(-1)*
     +    XLAM**(-2)*MS2*MG2 - 24*S**2*XLAM**(-2)*MS2 + 16*S**2*
     +    XLAM**(-2)*MG2 - 8*S**3*U1**(-1)*XLAM**(-2)*MS2 - 8*S**3*
     +    U1**(-1)*XLAM**(-2)*MG2 - 24*TG**(-1)*U1**(-1)*MS2*MG2**2 + 
     +    24*TG**(-1)*U1**(-1)*MS2**2*MG2 - 8*TG**(-1)*U1**(-1)*MS2**3
     +     + 8*TG**(-1)*U1**(-1)*MG2**3 - 32*TG**(-1)*MS2*MG2 + 16*
     +    TG**(-1)*MS2**2 + 16*TG**(-1)*MG2**2 - 8*UG*U1**(-1)*MS2 + 8*
     +    UG*U1**(-1)*MG2 + 16*UG*XLAM**(-2)*MS2*MG2 - 4*UG*XLAM**(-2)*
     +    MS2**2 - 12*UG*XLAM**(-2)*MG2**2 - 8*UG - 8*UG**2*XLAM**(-2)*
     +    MS2 + 8*UG**2*XLAM**(-2)*MG2 + 8*UG**3*U1**(-1)*XLAM**(-2)*
     +    MS2 - 8*UG**3*U1**(-1)*XLAM**(-2)*MG2 - 36*U1**(-1)*MS2*MG2
     +     + 8*U1**(-1)*MS2**2 + 28*U1**(-1)*MG2**2 - 24*XLAM**(-2)*MS2
     +    *MG2**2 + 24*XLAM**(-2)*MS2**2*MG2 - 8*XLAM**(-2)*MS2**3 + 8*
     +    XLAM**(-2)*MG2**3 + 4*MS2 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(1,2)*N*CK*FOUR**(-1) * ( 4*S**(-1)*UG*MG2
     +     - 4*S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 - 4*S**(-1)*MG2**2 + 
     +    4*S*TG**(-1)*MG2 - 16*TG**(-1)*MS2*MG2 + 8*TG**(-1)*MS2**2 + 
     +    8*TG**(-1)*MG2**2 - 8*UG - 8*MS2 + 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(1,2)*CQED*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MG2 + 4*S**(-1)*UG**2 - 4*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2
     +     + 4*UG + 4*MS2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(2,2)*N*CO*FOUR**(-1) * (  - 8*S*TG**(-1)*
     +    MS2 + 4*S + 4*S**2*TG**(-1) - 16*TG**(-1)*MS2*MG2 + 8*
     +    TG**(-1)*MS2**2 + 8*TG**(-1)*MG2**2 - 4*UG - 4*MS2 )
     +
      M2QGV = M2QGV + SK2C0C(2,2)*N*CK*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MG2 - 12*S**(-1)*MS2*MG2 + 12*S**(-1)*MG2**2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(3,1)*N*CK*FOUR**(-1) * ( 4*S**(-1)*UG*MS2
     +     - 4*S**(-1)*UG*MG2 - 8*S**(-1)*UG**2*U1**(-1)*MS2 + 8*
     +    S**(-1)*UG**2*U1**(-1)*MG2 - 4*S**(-1)*UG**2 - 48*S*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 32*S*TG**(-1)*U1**(-1)*MS2**2 + 16*S*
     +    TG**(-1)*U1**(-1)*MG2**2 - 24*S*U1**(-1)*MS2 + 24*S*U1**(-1)*
     +    MG2 - 16*S**2*TG**(-1)*U1**(-1)*MS2 + 16*S**2*TG**(-1)*
     +    U1**(-1)*MG2 + 64*TG**(-2)*MS2*MG2**2 - 32*TG**(-2)*MS2**2*
     +    MG2 - 32*TG**(-2)*MG2**3 - 16*TG**(-1)*U1**(-1)*MS2*MG2**2 + 
     +    32*TG**(-1)*U1**(-1)*MS2**2*MG2 - 16*TG**(-1)*U1**(-1)*MS2**3
     +     + 56*TG**(-1)*MS2*MG2 - 56*TG**(-1)*MG2**2 + 16*UG*U1**(-1)*
     +    MS2 - 16*UG*U1**(-1)*MG2 + 12*UG + 16*U1**(-1)*MS2**2 - 16*
     +    U1**(-1)*MG2**2 + 4*MS2 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(3,1)*CQED*FOUR**(-1) * ( 4*S**(-1)*UG*MG2
     +     - 4*S*U1**(-1)*MS2 + 8*S*U1**(-1)*MG2 - 8*S - 4*S**2*
     +    U1**(-1) - 4*UG + 4*U1**(-1)*MS2*MG2 - 4*U1**(-1)*MG2**2 + 8*
     +    MG2 )
     +
      M2QGV = M2QGV + SK2C0C(3,2)*N*CO*FOUR**(-1) * ( 8*S*TG**(-1)*MS2
     +     - 8*S*TG**(-1)*MG2 - 4*S - 4*S**2*TG**(-1) - 8*TG**(-1)*
     +    MS2**2 + 8*TG**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0C(3,2)*N*CK*FOUR**(-1) * ( 4*S**(-1)*UG**2
     +     - 8*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2 + 4*S**(-1)*MG2**2 + 
     +    4*UG + 4*MS2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(4,1)*N*CO*FOUR**(-1) * ( 4*S**(-1)*UG*MS2
     +     - 4*S**(-1)*UG*MG2 - 8*S**(-1)*UG**2*U1**(-1)*MS2 + 8*
     +    S**(-1)*UG**2*U1**(-1)*MG2 - 4*S**(-1)*UG**2 - 48*S*TG**(-1)*
     +    U1**(-1)*MS2*MG2 + 32*S*TG**(-1)*U1**(-1)*MS2**2 + 16*S*
     +    TG**(-1)*U1**(-1)*MG2**2 - 8*S*TG**(-1)*MS2 + 8*S*TG**(-1)*
     +    MG2 - 16*S*U1**(-1)*MS2 + 20*S*U1**(-1)*MG2 - 16*S**2*
     +    TG**(-1)*U1**(-1)*MS2 + 16*S**2*TG**(-1)*U1**(-1)*MG2 + 64*
     +    TG**(-2)*MS2*MG2**2 - 32*TG**(-2)*MS2**2*MG2 - 32*TG**(-2)*
     +    MG2**3 - 16*TG**(-1)*U1**(-1)*MS2*MG2**2 + 32*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 16*TG**(-1)*U1**(-1)*MS2**3 + 16*
     +    TG**(-1)*MS2*MG2 + 16*TG**(-1)*MS2**2 - 32*TG**(-1)*MG2**2 + 
     +    24*UG*U1**(-1)*MS2 - 24*UG*U1**(-1)*MG2 + 4*UG + 4*U1**(-1)*
     +    MS2*MG2 + 16*U1**(-1)*MS2**2 - 20*U1**(-1)*MG2**2 - 12*MS2 + 
     +    16*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(4,2)*N*CK*FOUR**(-1) * (  - 4*S**(-1)*
     +    UG**2 + 8*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 4*S**(-1)*
     +    MG2**2 + 4*S*TG**(-1)*MG2 - 16*TG**(-1)*MS2*MG2 + 16*TG**(-1)
     +    *MG2**2 + 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(4,2)*CQED*FOUR**(-1) * (  - 8*S**(-1)*UG*
     +    MG2 - 8*S**(-1)*MS2*MG2 + 8*S**(-1)*MG2**2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(5,1)*N*CO*FOUR**(-1) * ( 16*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 16*S**(-1)*UG*U1**(-1)*MG2**2 - 8*S**(-1)*
     +    UG*MS2 + 12*S**(-1)*UG*MG2 + 8*S**(-1)*UG**2*U1**(-1)*MS2 - 8
     +    *S**(-1)*UG**2*U1**(-1)*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2*MG2
     +     - 32*S*TG**(-1)*U1**(-1)*MS2**2 + 16*S*TG**(-1)*U1**(-1)*
     +    MG2**2 + 8*S*TG**(-1)*MS2 - 24*S*TG**(-1)*MG2 - 16*S*T1**(-1)
     +    *U1**(-1)*MS2*MG2 + 4*S*T1**(-1)*U1**(-1)*MS2**2 + 12*S*
     +    T1**(-1)*U1**(-1)*MG2**2 + 12*S*T1**(-1)*MS2 + 4*S*T1**(-1)*
     +    MG2 + 4*S*U1**(-1)*MS2 - 4*S*U1**(-1)*MG2 - 8*S + 16*S**2*
     +    TG**(-1)*U1**(-1)*MS2 - 16*S**2*TG**(-1)*U1**(-1)*MG2 - 4*
     +    S**2*T1**(-1)*U1**(-1)*MS2 + 4*S**2*T1**(-1)*U1**(-1)*MG2 + 
     +    32*TG**(-2)*MS2**2*MG2 - 32*TG**(-2)*MG2**3 - 16*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 + 16*TG**(-1)*U1**(-1)*MS2**3 + 48*
     +    TG**(-1)*MS2*MG2 - 16*TG**(-1)*MS2**2 - 32*TG**(-1)*MG2**2 - 
     +    12*UG*U1**(-1)*MS2 + 12*UG*U1**(-1)*MG2 + 40*T1**(-1)*
     +    U1**(-1)*MS2*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0C(5,1)*N*CO*FOUR**(-1) * (  - 24*T1**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 16*T1**(-1)*U1**(-1)*MG2**3 - 32*
     +    T1**(-1)*MS2*MG2 + 16*T1**(-1)*MG2**2 - 16*U1**(-1)*MS2*MG2
     +     - 20*U1**(-1)*MS2**2 + 28*U1**(-1)*MG2**2 + 20*MS2 - 16*MG2
     +     )
     +
      M2QGV = M2QGV + SK2C0C(6,1)*N*CK*FOUR**(-1) * ( 16*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 16*S**(-1)*UG*U1**(-1)*MG2**2 + 4*S**(-1)*
     +    UG*MG2 - 32*S*TG**(-1)*U1**(-1)*MS2*MG2 + 32*S*TG**(-1)*
     +    U1**(-1)*MG2**2 - 8*S*T1**(-1)*U1**(-1)*MS2*MG2 + 8*S*
     +    T1**(-1)*U1**(-1)*MG2**2 + 4*S*U1**(-1)*MG2 + 64*TG**(-2)*MS2
     +    *MG2**2 - 64*TG**(-2)*MG2**3 - 32*TG**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 32*TG**(-1)*U1**(-1)*MS2**2*MG2 + 56*TG**(-1)*MS2*
     +    MG2 - 56*TG**(-1)*MG2**2 + 24*T1**(-1)*U1**(-1)*MS2*MG2**2 - 
     +    16*T1**(-1)*U1**(-1)*MS2**2*MG2 - 8*T1**(-1)*U1**(-1)*MG2**3
     +     - 24*T1**(-1)*MS2*MG2 + 8*T1**(-1)*MG2**2 - 28*U1**(-1)*MS2*
     +    MG2 + 20*U1**(-1)*MG2**2 - 16*MG2 )
     +
      M2QGV = M2QGV + SK2C0C(6,1)*CQED*FOUR**(-1) * ( 4*S**(-1)*UG*MG2
     +     + 4*S*U1**(-1)*MG2 + 4*U1**(-1)*MS2*MG2 - 4*U1**(-1)*MG2**2
     +     + 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(1,1)*N*CO*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MG2 + 4*S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2 + 4*S**(-1)*MG2**2
     +     - 16*S*TG**(-2)*MG2**2 - 4*S*U1**(-1)*MG2 - 4*S + 16*
     +    TG**(-1)*MS2*MG2 - 8*TG**(-1)*MG2**2 + 8*UG*U1**(-1)*MS2 - 8*
     +    UG*U1**(-1)*MG2 + 4*UG + 8*U1**(-1)*MS2*MG2 + 4*MS2 + 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(2,1)*N*CK*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MG2 - 12*S**(-1)*MS2*MG2 + 12*S**(-1)*MG2**2 + 8*S*U1**(-1)*
     +    MS2 - 8*S*U1**(-1)*MG2 + 4*S + 8*UG*U1**(-1)*MS2 - 8*UG*
     +    U1**(-1)*MG2 + 4*UG - 4*MS2 - 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(2,1)*CO*FOUR**(-1) * ( 16*S*TG**(-2)*NS*
     +    MS2*MG2 - 16*S*TG**(-2)*MS2*MG2 + 8*TG**(-1)*NS*MS2*MG2 - 16*
     +    TG**(-1)*NS*MS2**2 - 8*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MS2**2
     +     - 8*UG*U1**(-1)*NS*MS2 + 8*UG*U1**(-1)*NS*MG2 + 8*UG*
     +    U1**(-1)*MS2 - 8*UG*U1**(-1)*MG2 - 16*U1**(-1)*NS*MS2*MG2 + 8
     +    *U1**(-1)*NS*MG2**2 + 16*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MG2**2
     +     + 8*NS*MS2 - 8*NS*MG2 - 8*MS2 + 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(2,1)*CQED*FOUR**(-1) * ( 8*S**(-1)*UG*MG2
     +     - 4*S**(-1)*UG**2 + 16*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 
     +    12*S**(-1)*MG2**2 - 8*S*U1**(-1)*MS2 + 16*S*U1**(-1)*MG2 - 4*
     +    S - 8*UG*U1**(-1)*MS2 + 8*UG*U1**(-1)*MG2 - 8*UG + 16*
     +    U1**(-1)*MS2*MG2 - 16*U1**(-1)*MG2**2 + 16*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(3,1)*N*CO*FOUR**(-1) * ( 4*S**(-1)*UG**2
     +     - 8*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2 + 4*S**(-1)*MG2**2 + 
     +    8*S*U1**(-1)*MS2 - 8*S*U1**(-1)*MG2 + 4*S - 8*UG*U1**(-1)*MS2
     +     + 8*UG*U1**(-1)*MG2 - 16*U1**(-1)*MS2*MG2 + 16*U1**(-1)*
     +    MG2**2 + 8*MS2 - 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(3,2)*N*CO*FOUR**(-1) * ( 8*S*TG**(-1)*
     +    UG**(-1)*MS2*MG2 - 4*S*TG**(-1)*UG**(-1)*MS2**2 - 4*S*
     +    TG**(-1)*UG**(-1)*MG2**2 - 16*S*TG**(-1)*U1**(-1)*MS2*MG2 + 
     +    16*S*TG**(-1)*U1**(-1)*MS2**2 + 8*S*TG**(-1)*MS2 - 8*S*
     +    TG**(-1)*MG2 + 4*S*UG**(-1)*MS2 - 4*S*UG**(-1)*MG2 - 8*S*
     +    U1**(-1)*MS2 + 8*S*U1**(-1)*MG2 - 4*S + 4*S**2*TG**(-1)*
     +    UG**(-1)*MS2 - 4*S**2*TG**(-1)*UG**(-1)*MG2 - 8*S**2*TG**(-1)
     +    *U1**(-1)*MS2 + 8*S**2*TG**(-1)*U1**(-1)*MG2 - 4*S**2*
     +    TG**(-1) + 8*TG**(-1)*U1**(-1)*MS2*MG2**2 + 8*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 8*TG**(-1)*U1**(-1)*MS2**3 - 8*TG**(-1)
     +    *U1**(-1)*MG2**3 + 32*UG*U1**(-2)*MS2*MG2 - 16*UG*U1**(-2)*
     +    MS2**2 - 16*UG*U1**(-2)*MG2**2 - 8*UG*U1**(-1)*MS2 + 8*UG*
     +    U1**(-1)*MG2 - 4*UG + 32*U1**(-2)*MS2*MG2**2 - 16*U1**(-2)*
     +    MS2**2*MG2 - 16*U1**(-2)*MG2**3 - 20*U1**(-1)*MS2*MG2 + 8*
     +    U1**(-1)*MS2**2 + 12*U1**(-1)*MG2**2 - 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(3,2)*N*CK*FOUR**(-1) * (  - 24*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*TG**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*MG2**3 + 4*
     +    S**(-1)*UG*MS2 - 4*S**(-1)*UG*MG2 + 4*S**(-1)*UG**2 - 16*
     +    S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 + 8*S**(-1)*MG2**2 - 32*S*
     +    TG**(-1)*U1**(-1)*MS2*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2**2 + 
     +    16*S*TG**(-1)*U1**(-1)*MG2**2 - 8*S*TG**(-1)*MS2 + 8*S*
     +    TG**(-1)*MG2 - 8*S*U1**(-1)*MS2 + 8*S*U1**(-1)*MG2 - 8*S**2*
     +    TG**(-1)*U1**(-1)*MS2 + 8*S**2*TG**(-1)*U1**(-1)*MG2 - 24*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 + 24*TG**(-1)*U1**(-1)*MS2**2*
     +    MG2 - 8*TG**(-1)*U1**(-1)*MS2**3 + 8*TG**(-1)*U1**(-1)*MG2**3
     +     - 32*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MS2**2 + 16*TG**(-1)*
     +    MG2**2 - 32*UG*U1**(-2)*MS2*MG2 + 16*UG*U1**(-2)*MS2**2 + 16*
     +    UG*U1**(-2)*MG2**2 + 24*UG*U1**(-1)*MS2 - 24*UG*U1**(-1)*MG2
     +     + 8*UG - 32*U1**(-2)*MS2*MG2**2 + 16*U1**(-2)*MS2**2*MG2 + 
     +    16*U1**(-2)*MG2**3 )
     +
      M2QGV = M2QGV + SK2C0D(3,2)*N*CK*FOUR**(-1) * ( 4*U1**(-1)*MS2*
     +    MG2 + 8*U1**(-1)*MS2**2 - 12*U1**(-1)*MG2**2 - 8*MS2 + 12*MG2
     +     )
     +
      M2QGV = M2QGV + SK2C0D(4,1)*N*CK*FOUR**(-1) * ( 4*S**(-1)*UG**2
     +     - 8*S**(-1)*MS2*MG2 + 4*S**(-1)*MS2**2 + 4*S**(-1)*MG2**2 - 
     +    4*S*U1**(-1)*MG2 + 4*S + 8*UG - 16*U1**(-1)*MS2*MG2 + 16*
     +    U1**(-1)*MG2**2 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(4,1)*CQED*FOUR**(-1) * ( 8*S**(-1)*UG*MG2
     +     - 4*S**(-1)*UG**2 + 16*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 
     +    12*S**(-1)*MG2**2 - 8*S*U1**(-1)*MS2 + 16*S*U1**(-1)*MG2 - 12
     +    *S - 4*S**2*U1**(-1) - 8*UG*U1**(-1)*MS2 + 8*UG*U1**(-1)*MG2
     +     - 12*UG + 16*U1**(-1)*MS2*MG2 - 16*U1**(-1)*MG2**2 + 4*MS2
     +     + 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(4,2)*N*CK*FOUR**(-1) * ( 24*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 - 24*S**(-1)*TG**(-1)*MS2**2*MG2 + 8*
     +    S**(-1)*TG**(-1)*MS2**3 - 8*S**(-1)*TG**(-1)*MG2**3 - 4*
     +    S**(-1)*UG*MS2 + 4*S**(-1)*UG*MG2 - 4*S**(-1)*UG**2 + 16*
     +    S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2 - 8*S**(-1)*MG2**2 + 30*S*
     +    TG**(-1)*UG**(-2)*MS2*MG2**2 - 30*S*TG**(-1)*UG**(-2)*MS2**2*
     +    MG2 + 10*S*TG**(-1)*UG**(-2)*MS2**3 - 10*S*TG**(-1)*UG**(-2)*
     +    MG2**3 - 40*S*TG**(-1)*UG**(-1)*MS2*MG2 + 20*S*TG**(-1)*
     +    UG**(-1)*MS2**2 + 20*S*TG**(-1)*UG**(-1)*MG2**2 + 16*S*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 16*S*TG**(-1)*U1**(-1)*MG2**2 - 
     +    30*S*TG**(-1)*MS2 + 34*S*TG**(-1)*MG2 + 20*S*UG**(-2)*MS2*MG2
     +     - 10*S*UG**(-2)*MS2**2 - 10*S*UG**(-2)*MG2**2 - 30*S*
     +    UG**(-1)*MS2 + 30*S*UG**(-1)*MG2 + 20*S**2*TG**(-1)*UG**(-2)*
     +    MS2*MG2 - 10*S**2*TG**(-1)*UG**(-2)*MS2**2 - 10*S**2*TG**(-1)
     +    *UG**(-2)*MG2**2 - 30*S**2*TG**(-1)*UG**(-1)*MS2 + 30*S**2*
     +    TG**(-1)*UG**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(4,2)*N*CK*FOUR**(-1) * ( 32*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 16*TG**(-1)*U1**(-1)*MS2**2*MG2 - 16*
     +    TG**(-1)*U1**(-1)*MG2**3 + 12*TG**(-1)*MS2*MG2 - 8*TG**(-1)*
     +    MS2**2 - 4*TG**(-1)*MG2**2 + 64*UG*U1**(-2)*MS2*MG2 - 32*UG*
     +    U1**(-2)*MS2**2 - 32*UG*U1**(-2)*MG2**2 - 16*UG*U1**(-1)*MS2
     +     + 16*UG*U1**(-1)*MG2 - 8*UG + 64*U1**(-2)*MS2*MG2**2 - 32*
     +    U1**(-2)*MS2**2*MG2 - 32*U1**(-2)*MG2**3 - 16*U1**(-1)*MS2*
     +    MG2 + 16*U1**(-1)*MG2**2 + 4*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(4,2)*CQED*FOUR**(-1) * (  - 24*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 24*S**(-1)*TG**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*TG**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*MG2**3 + 4*
     +    S**(-1)*UG*MS2 - 8*S**(-1)*UG*MG2 + 4*S**(-1)*UG**2 - 16*
     +    S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 + 8*S**(-1)*MG2**2 + 18*S*
     +    TG**(-1)*UG**(-2)*MS2*MG2**2 - 18*S*TG**(-1)*UG**(-2)*MS2**2*
     +    MG2 + 6*S*TG**(-1)*UG**(-2)*MS2**3 - 6*S*TG**(-1)*UG**(-2)*
     +    MG2**3 - 16*S*TG**(-1)*UG**(-1)*MS2*MG2 + 8*S*TG**(-1)*
     +    UG**(-1)*MS2**2 + 8*S*TG**(-1)*UG**(-1)*MG2**2 - 32*S*
     +    TG**(-1)*U1**(-1)*MS2*MG2 + 16*S*TG**(-1)*U1**(-1)*MS2**2 + 
     +    16*S*TG**(-1)*U1**(-1)*MG2**2 - 22*S*TG**(-1)*MS2 + 22*S*
     +    TG**(-1)*MG2 + 12*S*UG**(-2)*MS2*MG2 - 6*S*UG**(-2)*MS2**2 - 
     +    6*S*UG**(-2)*MG2**2 - 16*S*UG**(-1)*T1**(-1)*MS2*MG2 + 12*S*
     +    UG**(-1)*T1**(-1)*MS2**2 + 4*S*UG**(-1)*T1**(-1)*MG2**2 - 8*S
     +    *UG**(-1)*MS2 + 8*S*UG**(-1)*MG2 + 6*S*T1**(-1)*MS2 - 6*S*
     +    T1**(-1)*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(4,2)*CQED*FOUR**(-1) * (  - 8*S*U1**(-1)*
     +    MS2 + 8*S*U1**(-1)*MG2 + 12*S**2*TG**(-1)*UG**(-2)*MS2*MG2 - 
     +    6*S**2*TG**(-1)*UG**(-2)*MS2**2 - 6*S**2*TG**(-1)*UG**(-2)*
     +    MG2**2 - 14*S**2*TG**(-1)*UG**(-1)*MS2 + 14*S**2*TG**(-1)*
     +    UG**(-1)*MG2 - 8*S**2*TG**(-1)*U1**(-1)*MS2 + 8*S**2*TG**(-1)
     +    *U1**(-1)*MG2 + 6*S**2*UG**(-1)*T1**(-1)*MS2 - 6*S**2*
     +    UG**(-1)*T1**(-1)*MG2 - 24*TG**(-1)*U1**(-1)*MS2*MG2**2 + 24*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 - 8*TG**(-1)*U1**(-1)*MS2**3 + 8
     +    *TG**(-1)*U1**(-1)*MG2**3 - 32*TG**(-1)*MS2*MG2 + 16*TG**(-1)
     +    *MS2**2 + 16*TG**(-1)*MG2**2 - 16*UG**(-1)*MS2*MG2 + 12*
     +    UG**(-1)*MS2**2 + 4*UG**(-1)*MG2**2 - 32*UG*U1**(-2)*MS2*MG2
     +     + 16*UG*U1**(-2)*MS2**2 + 16*UG*U1**(-2)*MG2**2 + 16*UG*
     +    U1**(-1)*MS2 - 16*UG*U1**(-1)*MG2 + 4*UG - 16*T1**(-1)*MS2*
     +    MG2 + 12*T1**(-1)*MS2**2 + 4*T1**(-1)*MG2**2 - 32*U1**(-2)*
     +    MS2*MG2**2 + 16*U1**(-2)*MS2**2*MG2 + 16*U1**(-2)*MG2**3 + 8*
     +    U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SK2C0D(4,2)*CQED*FOUR**(-1) * (  - 8*U1**(-1)*
     +    MG2**2 - 8*MS2 + 8*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(5,2)*N*CO*FOUR**(-1) * ( 9*S*TG**(-1)*
     +    UG**(-2)*MS2*MG2**2 - 9*S*TG**(-1)*UG**(-2)*MS2**2*MG2 + 3*S*
     +    TG**(-1)*UG**(-2)*MS2**3 - 3*S*TG**(-1)*UG**(-2)*MG2**3 - 12*
     +    S*TG**(-1)*UG**(-1)*MS2*MG2 + 6*S*TG**(-1)*UG**(-1)*MS2**2 + 
     +    6*S*TG**(-1)*UG**(-1)*MG2**2 - 16*S*TG**(-1)*U1**(-1)*MS2**2
     +     + 16*S*TG**(-1)*U1**(-1)*MG2**2 - 13*S*TG**(-1)*MS2 - 3*S*
     +    TG**(-1)*MG2 + 6*S*UG**(-2)*MS2*MG2 - 3*S*UG**(-2)*MS2**2 - 3
     +    *S*UG**(-2)*MG2**2 - 9*S*UG**(-1)*MS2 + 9*S*UG**(-1)*MG2 + 8*
     +    S*U1**(-1)*MS2 - 8*S*U1**(-1)*MG2 + 4*S + 6*S**2*TG**(-1)*
     +    UG**(-2)*MS2*MG2 - 3*S**2*TG**(-1)*UG**(-2)*MS2**2 - 3*S**2*
     +    TG**(-1)*UG**(-2)*MG2**2 - 9*S**2*TG**(-1)*UG**(-1)*MS2 + 9*
     +    S**2*TG**(-1)*UG**(-1)*MG2 + 8*S**2*TG**(-1)*U1**(-1)*MS2 - 8
     +    *S**2*TG**(-1)*U1**(-1)*MG2 + 4*S**2*TG**(-1) - 8*TG**(-1)*
     +    UG**(-1)*MS2*MG2**2 + 8*TG**(-1)*UG**(-1)*MG2**3 - 8*TG**(-1)
     +    *U1**(-1)*MS2*MG2**2 + 8*TG**(-1)*U1**(-1)*MS2**2*MG2 + 8*
     +    TG**(-1)*U1**(-1)*MS2**3 )
     +
      M2QGV = M2QGV + SK2C0D(5,2)*N*CO*FOUR**(-1) * (  - 8*TG**(-1)*
     +    U1**(-1)*MG2**3 - 16*TG**(-1)*MS2*MG2 + 24*TG**(-1)*MG2**2 + 
     +    8*UG**(-1)*MG2**2 - 32*UG*U1**(-2)*MS2*MG2 + 16*UG*U1**(-2)*
     +    MS2**2 + 16*UG*U1**(-2)*MG2**2 - 64*U1**(-2)*MS2*MG2**2 + 48*
     +    U1**(-2)*MS2**2*MG2 + 16*U1**(-2)*MG2**3 + 16*U1**(-1)*MS2*
     +    MG2 - 8*U1**(-1)*MS2**2 - 8*U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SK2C0D(5,2)*N*CK*FOUR**(-1) * (  - 8*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 - 8*S**(-1)*TG**(-1)*MS2**2*MG2 + 8*
     +    S**(-1)*TG**(-1)*MS2**3 + 8*S**(-1)*TG**(-1)*MG2**3 - 4*
     +    S**(-1)*UG*MG2 - 8*S**(-1)*MS2**2 + 8*S**(-1)*MG2**2 + 16*S*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 16*S*TG**(-1)*U1**(-1)*MS2**2 + 8
     +    *S*TG**(-1)*MS2 - 8*S*TG**(-1)*MG2 + 8*S*U1**(-1)*MS2 - 8*S*
     +    U1**(-1)*MG2 + 8*S**2*TG**(-1)*U1**(-1)*MS2 - 8*S**2*TG**(-1)
     +    *U1**(-1)*MG2 - 8*TG**(-1)*U1**(-1)*MS2*MG2**2 - 8*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*U1**(-1)*MS2**3 + 8*TG**(-1)
     +    *U1**(-1)*MG2**3 + 16*TG**(-1)*MS2*MG2 - 16*TG**(-1)*MS2**2
     +     - 8*UG**(-1)*MG2**2 + 32*UG*U1**(-2)*MS2*MG2 - 16*UG*
     +    U1**(-2)*MS2**2 - 16*UG*U1**(-2)*MG2**2 - 16*UG*U1**(-1)*MS2
     +     + 16*UG*U1**(-1)*MG2 - 4*UG + 64*U1**(-2)*MS2*MG2**2 - 48*
     +    U1**(-2)*MS2**2*MG2 - 16*U1**(-2)*MG2**3 - 32*U1**(-1)*MS2*
     +    MG2 - 8*U1**(-1)*MS2**2 + 40*U1**(-1)*MG2**2 + 8*MS2 - 24*MG2
     +     )
     +
      M2QGV = M2QGV + SK2C0D(6,2)*N*CK*FOUR**(-1) * ( 32*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 - 16*S**(-1)*TG**(-1)*MS2**2*MG2 - 16*
     +    S**(-1)*TG**(-1)*MG2**3 + 4*S**(-1)*UG*MG2 + 16*S**(-1)*MS2*
     +    MG2 - 16*S**(-1)*MG2**2 - 4*S*TG**(-1)*MG2 - 8*TG**(-1)*
     +    UG**(-1)*MS2*MG2**2 + 8*TG**(-1)*UG**(-1)*MG2**3 + 32*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 - 32*TG**(-1)*U1**(-1)*MG2**3 + 
     +    4*TG**(-1)*MS2*MG2 + 4*TG**(-1)*MG2**2 + 16*UG**(-1)*MG2**2
     +     - 64*U1**(-2)*MS2*MG2**2 + 64*U1**(-2)*MS2**2*MG2 + 56*
     +    U1**(-1)*MS2*MG2 - 56*U1**(-1)*MG2**2 + 20*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(6,2)*CQED*FOUR**(-1) * (  - 32*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 16*S**(-1)*TG**(-1)*MS2**2*MG2 + 16*
     +    S**(-1)*TG**(-1)*MG2**3 - 8*S**(-1)*UG*MG2 - 16*S**(-1)*MS2*
     +    MG2 + 16*S**(-1)*MG2**2 - 16*S*TG**(-1)*U1**(-1)*MS2*MG2 + 16
     +    *S*TG**(-1)*U1**(-1)*MG2**2 - 32*TG**(-1)*U1**(-1)*MS2*MG2**2
     +     + 16*TG**(-1)*U1**(-1)*MS2**2*MG2 + 16*TG**(-1)*U1**(-1)*
     +    MG2**3 - 16*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MG2**2 - 8*
     +    UG**(-1)*MG2**2 + 32*U1**(-2)*MS2*MG2**2 - 32*U1**(-2)*MS2**2
     +    *MG2 - 44*U1**(-1)*MS2*MG2 + 44*U1**(-1)*MG2**2 - 12*MG2 )
     +
      M2QGV = M2QGV + SK2C0D(7,1)*CO*FOUR**(-1) * (  - 16*S*TG**(-2)*
     +    MT2*MG2 + 16*TG**(-1)*MS2*MT2 - 8*TG**(-1)*MT2*MG2 + 8*
     +    U1**(-1)*MS2*MT2 )
     +
      M2QGV = M2QGV + SK2C0D(8,1)*CO*FOUR**(-1) * ( 16*S*TG**(-2)*MS2*
     +    MG2 + 8*TG**(-1)*MS2*MG2 - 16*TG**(-1)*MS2**2 - 8*UG*U1**(-1)
     +    *MS2 + 8*UG*U1**(-1)*MG2 - 16*U1**(-1)*MS2*MG2 + 8*U1**(-1)*
     +    MG2**2 + 8*MS2 - 8*MG2 )
     +
      M2QGV = M2QGV + SK2D0(1,1)*N*CK*FOUR**(-1) * (  - 4*S*TG**(-1)*
     +    MS2*MG2 + 4*S*TG**(-1)*MG2**2 - 8*S*UG - 4*S*MS2 + 8*S*MG2 - 
     +    4*S**2 + 24*TG**(-1)*MS2*MG2**2 - 24*TG**(-1)*MS2**2*MG2 + 8*
     +    TG**(-1)*MS2**3 - 8*TG**(-1)*MG2**3 - 4*UG*MS2 + 4*UG*MG2 - 4
     +    *UG**2 + 16*MS2*MG2 - 8*MS2**2 - 8*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(1,1)*CQED*FOUR**(-1) * ( 16*S*UG*U1**(-1)*
     +    MS2 - 16*S*UG*U1**(-1)*MG2 + 12*S*UG - 4*S*U1**(-1)*MS2*MG2
     +     + 4*S*U1**(-1)*MG2**2 - 4*S*MS2 + 12*S**2*U1**(-1)*MS2 - 16*
     +    S**2*U1**(-1)*MG2 + 12*S**2 + 4*S**3*U1**(-1) - 4*UG*MS2 + 4*
     +    UG*MG2 + 8*UG**2*U1**(-1)*MS2 - 8*UG**2*U1**(-1)*MG2 + 4*
     +    UG**2 )
     +
      M2QGV = M2QGV + SK2D0(1,2)*N*CK*FOUR**(-1) * ( 16*S*UG*U1**(-1)*
     +    MS2 - 16*S*UG*U1**(-1)*MG2 + 4*S*UG + 16*S*U1**(-1)*MS2*MG2
     +     - 16*S*U1**(-1)*MG2**2 - 12*S*MS2 + 12*S*MG2 + 4*UG**2 - 8*
     +    MS2*MG2 + 4*MS2**2 + 4*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(2,1)*N*CO*FOUR**(-1) * (  - 16*S*TG**(-1)*
     +    MS2*MG2 + 16*S*TG**(-1)*MG2**2 + 4*UG**2 - 8*MS2*MG2 + 4*
     +    MS2**2 + 4*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(2,2)*N*CK*FOUR**(-1) * (  - 4*S*U1**(-1)*
     +    MS2*MG2 + 4*S*U1**(-1)*MG2**2 - 4*S*MG2 - 4*UG*MS2 + 4*UG*MG2
     +     + 8*UG**2*U1**(-1)*MS2 - 8*UG**2*U1**(-1)*MG2 + 4*UG**2 )
     +
      M2QGV = M2QGV + SK2D0(3,2)*N*CK*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MS2*MG2 + 4*S**(-1)*UG*MG2**2 + 24*S**(-1)*MS2*MG2**2 - 12*
     +    S**(-1)*MS2**2*MG2 - 12*S**(-1)*MG2**3 - 4*S*U1**(-1)*MS2*MG2
     +     + 4*S*U1**(-1)*MG2**2 - 4*S*MG2 + 16*UG*U1**(-1)*MS2*MG2 - 
     +    16*UG*U1**(-1)*MG2**2 + 4*UG*MG2 + 32*U1**(-1)*MS2*MG2**2 - 
     +    32*U1**(-1)*MG2**3 - 12*MS2*MG2 + 28*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(4,1)*N*CK*FOUR**(-1) * ( 4*S**(-1)*UG*MS2*
     +    MG2 - 4*S**(-1)*UG*MG2**2 - 24*S**(-1)*MS2*MG2**2 + 12*
     +    S**(-1)*MS2**2*MG2 + 12*S**(-1)*MG2**3 - 4*S*TG**(-1)*MS2*MG2
     +     + 4*S*TG**(-1)*MG2**2 + 8*S*MG2 + 16*TG**(-1)*MS2**2*MG2 - 
     +    16*TG**(-1)*MG2**3 + 4*UG*MG2 + 8*MS2*MG2 - 24*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(4,1)*CQED*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MS2*MG2 + 4*S**(-1)*UG*MG2**2 + 24*S**(-1)*MS2*MG2**2 - 12*
     +    S**(-1)*MS2**2*MG2 - 12*S**(-1)*MG2**3 - 16*S*U1**(-1)*MS2*
     +    MG2 + 32*S*U1**(-1)*MG2**2 - 12*S*MG2 - 8*S**2*U1**(-1)*MG2
     +     - 4*UG*MG2 + 32*U1**(-1)*MS2*MG2**2 - 32*U1**(-1)*MG2**3 - 
     +    12*MS2*MG2 + 28*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(5,2)*N*CK*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MS2*MG2 + 4*S**(-1)*UG*MG2**2 + 12*S**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*MS2**3 - 4*S**(-1)*MG2**3 + 16*S*UG*U1**(-1)*MS2 - 16
     +    *S*UG*U1**(-1)*MG2 + 4*S*UG + 16*S*U1**(-1)*MS2*MG2 - 16*S*
     +    U1**(-1)*MG2**2 - 12*S*MS2 + 12*S*MG2 - 48*UG*U1**(-1)*MS2*
     +    MG2 + 48*UG*U1**(-1)*MG2**2 + 12*UG*MS2 - 24*UG*MG2 - 16*
     +    UG**2*U1**(-1)*MS2 + 16*UG**2*U1**(-1)*MG2 - 32*U1**(-1)*MS2*
     +    MG2**2 + 32*U1**(-1)*MG2**3 + 4*MS2*MG2 + 20*MS2**2 - 24*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(6,1)*N*CO*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MS2*MG2 + 4*S**(-1)*UG*MG2**2 + 12*S**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*MS2**3 - 4*S**(-1)*MG2**3 - 16*S*TG**(-1)*MS2*MG2 + 
     +    16*S*TG**(-1)*MG2**2 - 12*S*MS2 + 4*S**2 + 16*TG**(-1)*MS2**2
     +    *MG2 - 16*TG**(-1)*MG2**3 - 4*UG*MG2 - 16*MS2*MG2 + 16*MS2**2
     +     )
     +
      M2QGV = M2QGV + SK2D0(7,1)*N*CK*FOUR**(-1) * ( 4*S*TG**(-1)*MS2*
     +    MG2 - 4*S*TG**(-1)*MG2**2 + 4*S*U1**(-1)*MS2*MG2 - 4*S*
     +    U1**(-1)*MG2**2 + 4*S*MG2 + 32*TG**(-1)*MS2*MG2**2 - 16*
     +    TG**(-1)*MS2**2*MG2 - 16*TG**(-1)*MG2**3 + 16*UG*U1**(-1)*MS2
     +    *MG2 - 16*UG*U1**(-1)*MG2**2 + 8*UG*MG2 + 8*MS2*MG2 - 8*
     +    MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(7,1)*CQED*FOUR**(-1) * (  - 16*S**(-1)*UG*
     +    MS2*MG2 + 16*S**(-1)*UG*MG2**2 - 8*S**(-1)*UG**2*MG2 + 16*
     +    S**(-1)*MS2*MG2**2 - 8*S**(-1)*MS2**2*MG2 - 8*S**(-1)*MG2**3
     +     - 4*S*U1**(-1)*MS2*MG2 + 4*S*U1**(-1)*MG2**2 - 4*S*MG2 - 16*
     +    UG*U1**(-1)*MS2*MG2 + 16*UG*U1**(-1)*MG2**2 - 12*UG*MG2 - 4*
     +    MS2*MG2 + 4*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(8,1)*N*CO*FOUR**(-1) * ( 16*S*TG**(-1)*MS2*
     +    MG2 - 16*S*TG**(-1)*MG2**2 + 16*S*UG*U1**(-1)*MS2 - 16*S*UG*
     +    U1**(-1)*MG2 + 4*S*UG + 16*S*U1**(-1)*MS2*MG2 - 16*S*U1**(-1)
     +    *MG2**2 - 12*S*MS2 + 12*S*MG2 + 8*MS2**2 - 8*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(9,1)*N*CK*FOUR**(-1) * ( 4*S*TG**(-1)*MS2*
     +    MG2 - 4*S*TG**(-1)*MG2**2 - 16*S*UG*U1**(-1)*MS2 + 16*S*UG*
     +    U1**(-1)*MG2 - 4*S*UG - 16*S*U1**(-1)*MS2*MG2 + 16*S*U1**(-1)
     +    *MG2**2 + 12*S*MS2 - 12*S*MG2 + 8*TG**(-1)*MS2*MG2**2 + 8*
     +    TG**(-1)*MS2**2*MG2 - 8*TG**(-1)*MS2**3 - 8*TG**(-1)*MG2**3
     +     - 16*UG*U1**(-1)*MS2*MG2 + 16*UG*U1**(-1)*MG2**2 + 12*UG*MS2
     +     - 8*UG*MG2 - 16*UG**2*U1**(-1)*MS2 + 16*UG**2*U1**(-1)*MG2
     +     - 4*UG**2 + 8*MS2**2 - 8*MG2**2 )
     +
      M2QGV = M2QGV + SK2D0(9,1)*CQED*FOUR**(-1) * (  - 4*S**(-1)*UG*
     +    MS2*MG2 + 4*S**(-1)*UG*MS2**2 - 4*S**(-1)*UG**2*MG2 + 4*
     +    S**(-1)*UG**3 + 16*S*UG*U1**(-1)*MS2 - 16*S*UG*U1**(-1)*MG2
     +     + 4*S*UG + 16*S*U1**(-1)*MS2*MG2 - 16*S*U1**(-1)*MG2**2 - 12
     +    *S*MS2 + 12*S*MG2 + 16*UG*U1**(-1)*MS2*MG2 - 16*UG*U1**(-1)*
     +    MG2**2 - 12*UG*MS2 + 8*UG*MG2 + 16*UG**2*U1**(-1)*MS2 - 16*
     +    UG**2*U1**(-1)*MG2 + 8*UG**2 - 4*MS2*MG2 + 4*MS2**2 )
     +
      M2QGV = M2QGV + SK2D0(10,1)*N*CO*FOUR**(-1) * ( 16*S*TG**(-1)*MS2
     +    *MG2 - 16*S*TG**(-1)*MG2**2 - 4*S*U1**(-1)*MS2*MG2 + 4*S*
     +    U1**(-1)*MG2**2 - 4*S*MG2 + 32*TG**(-1)*MS2*MG2**2 - 16*
     +    TG**(-1)*MS2**2*MG2 - 16*TG**(-1)*MG2**3 - 16*UG*U1**(-1)*MS2
     +    *MG2 + 16*UG*U1**(-1)*MG2**2 + 4*UG*MS2 - 8*UG*MG2 - 8*UG**2*
     +    U1**(-1)*MS2 + 8*UG**2*U1**(-1)*MG2 - 4*UG**2 + 16*MS2*MG2 - 
     +    16*MG2**2 )
     +
      M2QGV = M2QGV + SOF2(1)*N*CO*FOUR**(-1) * (  - 64*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 - 64*S**(-1)*TG**(-1)*MS2**2*MG2 + 64*
     +    S**(-1)*TG**(-1)*MS2**3 + 64*S**(-1)*TG**(-1)*MG2**3 - 32*
     +    S**(-1)*UG*MS2 - 32*S**(-1)*UG*MG2 - 32*S**(-1)*MS2**2 + 32*
     +    S**(-1)*MG2**2 + 64*S*TG**(-2)*MS2*MG2 - 64*S*TG**(-2)*MG2**2
     +     + 96*S*TG**(-1)*MS2 - 32*S*TG**(-1)*MG2 - 32*S - 32*S**2*
     +    TG**(-1) - 64*TG**(-2)*MS2**2*MG2 + 64*TG**(-2)*MG2**3 + 128*
     +    TG**(-1)*MS2*MG2 - 128*TG**(-1)*MS2**2 + 32*UG + 64*MS2 )
     +
      M2QGV = M2QGV + SOF2(1)*N*CK*FOUR**(-1) * ( 128*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 128*S**(-1)*UG*U1**(-1)*MG2**2 - 32*
     +    S**(-1)*UG*MS2 + 96*S**(-1)*UG*MG2 + 64*S**(-1)*UG**2*
     +    U1**(-1)*MS2 - 64*S**(-1)*UG**2*U1**(-1)*MG2 - 32*S**(-1)*
     +    MS2**2 + 32*S**(-1)*MG2**2 + 64*S*U1**(-2)*MS2*MG2 - 64*S*
     +    U1**(-2)*MS2**2 - 64*S*U1**(-1)*MS2 + 64*S*U1**(-1)*MG2 - 32*
     +    S - 128*UG*U1**(-2)*MS2*MG2 + 64*UG*U1**(-2)*MS2**2 + 64*UG*
     +    U1**(-2)*MG2**2 - 32*UG - 256*U1**(-2)*MS2*MG2**2 + 192*
     +    U1**(-2)*MS2**2*MG2 + 64*U1**(-2)*MG2**3 + 256*U1**(-1)*MS2*
     +    MG2 - 64*U1**(-1)*MS2**2 - 192*U1**(-1)*MG2**2 + 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(2)*N*CO*FOUR**(-1) * (  - 256*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 128*S**(-1)*TG**(-1)*MS2**2*MG2 + 128*
     +    S**(-1)*TG**(-1)*MG2**3 - 64*S**(-1)*UG*MG2 - 64*S**(-1)*MS2*
     +    MG2 + 64*S**(-1)*MG2**2 + 128*S*TG**(-1)*U1**(-1)*MS2*MG2 - 
     +    128*S*TG**(-1)*U1**(-1)*MG2**2 + 128*S*TG**(-1)*MG2 - 256*
     +    TG**(-2)*MS2*MG2**2 + 256*TG**(-2)*MG2**3 - 128*TG**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 128*TG**(-1)*U1**(-1)*MG2**3 - 128*
     +    TG**(-1)*MS2*MG2 + 128*TG**(-1)*MG2**2 + 128*U1**(-2)*MS2*
     +    MG2**2 - 128*U1**(-2)*MS2**2*MG2 + 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(2)*N*CK*FOUR**(-1) * ( 128*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 - 128*S**(-1)*UG*U1**(-1)*MG2**2 + 64*
     +    S**(-1)*UG*MG2 - 64*S**(-1)*MS2*MG2 + 64*S**(-1)*MG2**2 - 128
     +    *U1**(-2)*MS2*MG2**2 + 128*U1**(-2)*MS2**2*MG2 + 128*U1**(-1)
     +    *MS2*MG2 - 128*U1**(-1)*MG2**2 + 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(3)*N*CO*FOUR**(-1) * ( 64*S**(-1)*TG**(-1)*
     +    MS2*MG2**2 - 128*S**(-1)*TG**(-1)*MS2**2*MG2 + 64*S**(-1)*
     +    TG**(-1)*MS2**3 - 32*S**(-1)*UG*MS2 + 32*S**(-1)*MS2*MG2 - 32
     +    *S**(-1)*MS2**2 - 64*S*TG**(-1)*U1**(-1)*MS2*MG2 + 128*S*
     +    TG**(-1)*U1**(-1)*MS2**2 - 64*S*TG**(-1)*U1**(-1)*MG2**2 + 64
     +    *S*TG**(-1)*MG2 - 64*S*U1**(-1)*MS2 + 64*S*U1**(-1)*MG2 - 64*
     +    S**2*TG**(-1)*U1**(-1)*MS2 + 64*S**2*TG**(-1)*U1**(-1)*MG2 + 
     +    128*TG**(-2)*MS2*MG2**2 - 128*TG**(-2)*MS2**2*MG2 + 64*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 - 64*TG**(-1)*U1**(-1)*MS2**3 + 
     +    64*TG**(-1)*MS2*MG2 - 64*TG**(-1)*MG2**2 + 128*UG*U1**(-2)*
     +    MS2*MG2 - 64*UG*U1**(-2)*MS2**2 - 64*UG*U1**(-2)*MG2**2 + 64*
     +    UG*U1**(-1)*MS2 - 64*UG*U1**(-1)*MG2 + 192*U1**(-2)*MS2*
     +    MG2**2 - 128*U1**(-2)*MS2**2*MG2 - 64*U1**(-2)*MG2**3 + 64*
     +    U1**(-1)*MS2**2 - 64*U1**(-1)*MG2**2 - 32*MS2 + 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(3)*N*CK*FOUR**(-1) * (  - 64*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 128*S**(-1)*TG**(-1)*MS2**2*MG2 - 64*
     +    S**(-1)*TG**(-1)*MS2**3 + 128*S**(-1)*UG*U1**(-1)*MS2*MG2 - 
     +    128*S**(-1)*UG*U1**(-1)*MG2**2 - 32*S**(-1)*UG*MS2 + 128*
     +    S**(-1)*UG*MG2 + 128*S**(-1)*UG**2*U1**(-1)*MS2 - 128*S**(-1)
     +    *UG**2*U1**(-1)*MG2 + 32*S**(-1)*MS2*MG2 - 32*S**(-1)*MS2**2
     +     + 64*S*TG**(-1)*U1**(-1)*MS2*MG2 - 128*S*TG**(-1)*U1**(-1)*
     +    MS2**2 + 64*S*TG**(-1)*U1**(-1)*MG2**2 - 64*S*TG**(-1)*MG2 + 
     +    64*S*U1**(-1)*MS2 - 64*S*U1**(-1)*MG2 + 64*S**2*TG**(-1)*
     +    U1**(-1)*MS2 - 64*S**2*TG**(-1)*U1**(-1)*MG2 - 128*TG**(-2)*
     +    MS2*MG2**2 + 128*TG**(-2)*MS2**2*MG2 - 64*TG**(-1)*U1**(-1)*
     +    MS2*MG2**2 + 64*TG**(-1)*U1**(-1)*MS2**3 - 64*TG**(-1)*MS2*
     +    MG2 + 64*TG**(-1)*MG2**2 - 384*UG*U1**(-2)*MS2*MG2 + 192*UG*
     +    U1**(-2)*MS2**2 + 192*UG*U1**(-2)*MG2**2 - 64*UG*U1**(-1)*MS2
     +     + 64*UG*U1**(-1)*MG2 - 576*U1**(-2)*MS2*MG2**2 + 384*
     +    U1**(-2)*MS2**2*MG2 )
     +
      M2QGV = M2QGV + SOF2(3)*N*CK*FOUR**(-1) * ( 192*U1**(-2)*MG2**3
     +     + 128*U1**(-1)*MS2*MG2 - 64*U1**(-1)*MS2**2 - 64*U1**(-1)*
     +    MG2**2 + 96*MS2 - 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(3)*CQED*FOUR**(-1) * (  - 64*S**(-1)*UG*
     +    U1**(-1)*MS2*MG2 + 64*S**(-1)*UG*U1**(-1)*MG2**2 + 32*S**(-1)
     +    *UG*MS2 - 64*S**(-1)*UG*MG2 - 64*S**(-1)*UG**2*U1**(-1)*MS2
     +     + 64*S**(-1)*UG**2*U1**(-1)*MG2 - 32*S**(-1)*MS2*MG2 + 32*
     +    S**(-1)*MS2**2 + 128*UG*U1**(-2)*MS2*MG2 - 64*UG*U1**(-2)*
     +    MS2**2 - 64*UG*U1**(-2)*MG2**2 + 192*U1**(-2)*MS2*MG2**2 - 
     +    128*U1**(-2)*MS2**2*MG2 - 64*U1**(-2)*MG2**3 - 64*U1**(-1)*
     +    MS2*MG2 + 64*U1**(-1)*MG2**2 - 32*MS2 )
     +
      M2QGV = M2QGV + SOF2(4)*N*CO*FOUR**(-1) * (  - 64*S*TG**(-2)*MS2*
     +    MG2 + 64*S*TG**(-2)*MG2**2 - 64*S*TG**(-1)*MS2 + 64*S*
     +    TG**(-1)*MG2 + 32*S + 32*S**2*TG**(-1) - 64*TG**(-1)*MS2*MG2
     +     + 64*TG**(-1)*MS2**2 + 64*UG*U1**(-1)*MS2 - 64*UG*U1**(-1)*
     +    MG2 + 64*U1**(-1)*MS2*MG2 - 64*U1**(-1)*MG2**2 - 64*MS2 + 64*
     +    MG2 )
     +
      M2QGV = M2QGV + SOF2(4)*N*CK*FOUR**(-1) * (  - 32*S**(-1)*UG**2
     +     + 64*S**(-1)*MS2*MG2 - 32*S**(-1)*MS2**2 - 32*S**(-1)*MG2**2
     +     - 64*UG*U1**(-1)*MS2 + 64*UG*U1**(-1)*MG2 - 32*UG - 64*
     +    U1**(-1)*MS2*MG2 + 64*U1**(-1)*MG2**2 + 32*MS2 - 32*MG2 )
     +
      M2QGV = M2QGV + SOF2(5)*N*CK*FOUR**(-1) * (  - 192*S**(-1)*
     +    TG**(-1)*MS2*MG2**2 + 192*S**(-1)*TG**(-1)*MS2**2*MG2 - 64*
     +    S**(-1)*TG**(-1)*MS2**3 + 64*S**(-1)*TG**(-1)*MG2**3 + 32*
     +    S**(-1)*UG*MS2 - 32*S**(-1)*UG*MG2 + 32*S**(-1)*UG**2 - 128*
     +    S**(-1)*MS2*MG2 + 64*S**(-1)*MS2**2 + 64*S**(-1)*MG2**2 + 128
     +    *S*TG**(-1)*U1**(-1)*MS2*MG2 - 128*S*TG**(-1)*U1**(-1)*MS2**2
     +     - 64*S*U1**(-2)*MS2*MG2 + 64*S*U1**(-2)*MS2**2 + 128*S*
     +    U1**(-1)*MS2 - 128*S*U1**(-1)*MG2 + 32*S + 64*S**2*TG**(-1)*
     +    U1**(-1)*MS2 - 64*S**2*TG**(-1)*U1**(-1)*MG2 - 256*TG**(-2)*
     +    MS2*MG2**2 + 128*TG**(-2)*MS2**2*MG2 + 128*TG**(-2)*MG2**3 - 
     +    64*TG**(-1)*U1**(-1)*MS2*MG2**2 - 64*TG**(-1)*U1**(-1)*MS2**2
     +    *MG2 + 64*TG**(-1)*U1**(-1)*MS2**3 + 64*TG**(-1)*U1**(-1)*
     +    MG2**3 - 256*TG**(-1)*MS2*MG2 + 256*TG**(-1)*MG2**2 - 128*UG*
     +    U1**(-2)*MS2*MG2 + 64*UG*U1**(-2)*MS2**2 + 64*UG*U1**(-2)*
     +    MG2**2 - 128*UG*U1**(-1)*MS2 + 128*UG*U1**(-1)*MG2 - 128*
     +    U1**(-2)*MS2*MG2**2 )
     +
      M2QGV = M2QGV + SOF2(5)*N*CK*FOUR**(-1) * ( 64*U1**(-2)*MS2**2*
     +    MG2 + 64*U1**(-2)*MG2**3 - 192*U1**(-1)*MS2*MG2 + 192*
     +    U1**(-1)*MG2**2 + 96*MS2 - 96*MG2 )
     +
      M2QGV = M2QGV + SOF2(5)*CQED*FOUR**(-1) * ( 32*S**(-1)*UG*MS2 - 
     +    32*S**(-1)*UG*MG2 - 64*S**(-1)*UG**2*U1**(-1)*MS2 + 64*
     +    S**(-1)*UG**2*U1**(-1)*MG2 - 32*S**(-1)*UG**2 + 64*S*U1**(-2)
     +    *MS2*MG2 - 64*S*U1**(-2)*MS2**2 - 64*S*U1**(-1)*MS2 + 64*S*
     +    U1**(-1)*MG2 - 32*S + 128*UG*U1**(-2)*MS2*MG2 - 64*UG*
     +    U1**(-2)*MS2**2 - 64*UG*U1**(-2)*MG2**2 - 64*UG*U1**(-1)*MS2
     +     + 64*UG*U1**(-1)*MG2 - 64*UG + 128*U1**(-2)*MS2*MG2**2 - 64*
     +    U1**(-2)*MS2**2*MG2 - 64*U1**(-2)*MG2**3 + 64*U1**(-1)*MS2*
     +    MG2 - 64*U1**(-1)*MS2**2 - 32*MS2 + 32*MG2 )
     +
      M2QGV = M2QGV + SOF2(6)*N*CO*FOUR**(-1) * (  - 32*S**(-1)*UG**2
     +     + 64*S**(-1)*MS2*MG2 - 32*S**(-1)*MS2**2 - 32*S**(-1)*MG2**2
     +     + 64*S*U1**(-2)*MS2*MG2 - 64*S*U1**(-2)*MS2**2 - 64*S*
     +    U1**(-1)*MS2 + 64*S*U1**(-1)*MG2 - 32*S + 128*TG**(-1)*MS2*
     +    MG2 - 128*TG**(-1)*MG2**2 + 64*UG*U1**(-1)*MS2 - 64*UG*
     +    U1**(-1)*MG2 + 192*U1**(-1)*MS2*MG2 - 64*U1**(-1)*MS2**2 - 
     +    128*U1**(-1)*MG2**2 - 64*MS2 + 64*MG2 )
     +
      M2QGV = M2QGV + SOF2(7)*N*CO*FOUR**(-1) * (  - 64*S*TG**(-2)*MS2*
     +    MG2 + 64*S*TG**(-2)*MG2**2 + 128*S*TG**(-1)*U1**(-1)*MS2*MG2
     +     - 128*S*TG**(-1)*U1**(-1)*MS2**2 - 32*S*TG**(-1)*MS2 + 32*S*
     +    TG**(-1)*MG2 + 64*S*U1**(-1)*MS2 - 64*S*U1**(-1)*MG2 + 32*S
     +     + 64*S**2*TG**(-1)*U1**(-1)*MS2 - 64*S**2*TG**(-1)*U1**(-1)*
     +    MG2 + 32*S**2*TG**(-1) - 128*TG**(-2)*MS2*MG2**2 + 64*
     +    TG**(-2)*MS2**2*MG2 + 64*TG**(-2)*MG2**3 - 64*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 64*TG**(-1)*U1**(-1)*MS2**2*MG2 + 64*
     +    TG**(-1)*U1**(-1)*MS2**3 + 64*TG**(-1)*U1**(-1)*MG2**3 - 64*
     +    TG**(-1)*MS2*MG2 + 64*TG**(-1)*MG2**2 - 128*UG*U1**(-2)*MS2*
     +    MG2 + 64*UG*U1**(-2)*MS2**2 + 64*UG*U1**(-2)*MG2**2 - 128*
     +    U1**(-2)*MS2*MG2**2 + 64*U1**(-2)*MS2**2*MG2 + 64*U1**(-2)*
     +    MG2**3 + 64*U1**(-1)*MS2*MG2 - 64*U1**(-1)*MS2**2 )
     +
      M2QGV = M2QGV + SOF2(7)*N*CK*FOUR**(-1) * ( 32*S**(-1)*UG*MS2 - 
     +    32*S**(-1)*UG*MG2 - 64*S**(-1)*UG**2*U1**(-1)*MS2 + 64*
     +    S**(-1)*UG**2*U1**(-1)*MG2 - 32*S**(-1)*UG**2 + 128*UG*
     +    U1**(-2)*MS2*MG2 - 64*UG*U1**(-2)*MS2**2 - 64*UG*U1**(-2)*
     +    MG2**2 - 64*UG*U1**(-1)*MS2 + 64*UG*U1**(-1)*MG2 - 32*UG + 
     +    128*U1**(-2)*MS2*MG2**2 - 64*U1**(-2)*MS2**2*MG2 - 64*
     +    U1**(-2)*MG2**3 - 64*U1**(-1)*MS2*MG2 + 64*U1**(-1)*MG2**2 )
     +
      M2QGV = M2QGV + SOF2(8)*N*CO*FOUR**(-1) * ( 64*S*TG**(-2)*MS2*MG2
     +     - 64*S*TG**(-2)*MG2**2 + 64*S*TG**(-1)*MS2 - 64*S*TG**(-1)*
     +    MG2 - 32*S - 32*S**2*TG**(-1) + 128*TG**(-1)*MS2*MG2 - 64*
     +    TG**(-1)*MS2**2 - 64*TG**(-1)*MG2**2 + 32*UG + 32*MS2 - 32*
     +    MG2 )
     +
      M2QGV = M2QGV + SOF2(8)*N*CK*FOUR**(-1) * ( 64*S*U1**(-2)*MS2*MG2
     +     - 64*S*U1**(-2)*MS2**2 - 64*S*U1**(-1)*MS2 + 64*S*U1**(-1)*
     +    MG2 - 32*S - 32*UG + 128*U1**(-1)*MS2*MG2 - 64*U1**(-1)*
     +    MS2**2 - 64*U1**(-1)*MG2**2 - 32*MS2 + 32*MG2 )


      END IF
      
      M2QGV = 2.D0 * M2QGV

      Q2MS2 = (MS + MG)**2/4.D0/MS**2
            
      DSGQGV = 
     +     ALPHAS**3 * AVG * M2QGV /4.D0 /S**2 *CONV
     +     + ALPHAS/PI*DSGQGB(ALPHAS,S,TG,MS,MG) 
     +     * (LOG(MG2/MS2)*(-1.D0) + LOG(MT2/MS2)*(-1.D0/3.D0))
     +     +  DSGQGR(ALPHAS,S,TG,MS,MG)
     +     +  DSGQG1(ALPHAS,S,TG,MS,MG,Q2MS2)
C***  CHANGES TO SUBTRACTED MSBAR, TO DRBAR FOR THE YUKAWA-COUPLING 
C***  AND TO THE SCALE Q**2 = (MS+MG)**2/4
      RETURN
      END



      REAL*8 FUNCTION DSGQGH(ALPHAS,S,TG,S4,MS,MG)
C***  CROSS SECTIONS FOR Q + G -> SQ + GL +G
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,S4G,T1,TG,U1,UG,MS,MG,MS2,MG2,M2
      REAL*8 NS,CONV,N,CO,CK,CQED,AVG,TWO
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 AHP1P1, ABP1P1, ABP1M1, ABP1P2, A4P2P2, A4P1P2, A4P1P1
      REAL*8 A4M1P2, A4M1P1, A4P0P2, A4P1P0, ABP2P2
      REAL*8 ABP2P1, ABP2P0
      REAL*8 M2QGH,Q2MS2,DSGQG3
      REAL*8 ANG2(1:44), COLO2(1:9)
           
      CONV = 389379660.D0
      NS = 6.D0

      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0) /N
      TWO = 2.D0

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6) + M2
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9) 
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) + MS2 +MG2
      ANB(1,2) = -ANB(1,1)
      ANC(1,2) = -ANC(1,1)
      ANA(1,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(1,3) = -ANA(1,3)
      ANC(1,3) =  0.D0
      ANA(1,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(1,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)-2.D0*ANGDEF(3)*ANGDEF(4)
      ANC(1,4) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(1,5) = -ANB(1,3)
      ANC(1,5) = -ANC(1,3)
      ANA(1,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(1,6) = -ANB(1,4)
      ANC(1,6) = -ANC(1,4)
      ANA(1,7) = +ANA(1,1) -M2
      ANB(1,7) = +ANB(1,1)
      ANC(1,7) = +ANC(1,1)
      ANA(1,8) = +ANA(1,5) -M2
      ANB(1,8) = +ANB(1,5)
      ANC(1,8) = +ANC(1,5)
      ANA(1,9) = +ANA(1,6) -M2
      ANB(1,9) = +ANB(1,6)
      ANC(1,9) = +ANC(1,6)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(3,2) = -ANB(3,1)
      ANC(3,2) = -ANC(3,1)
      ANA(3,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(3,3) = 
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10) -2.D0*ANGDEF(2)*ANGDEF(4)
      ANC(3,3) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(3,4) = -ANA(3,4)
      ANC(3,4) =  0.D0
      ANA(3,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(3,5) = -ANB(3,3)
      ANC(3,5) = -ANC(3,3)
      ANA(3,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(3,6) = -ANB(3,4)
      ANC(3,6) = -ANC(3,4)
      ANA(3,7) = +ANA(3,1) -M2
      ANB(3,7) = +ANB(3,1)
      ANC(3,7) = +ANC(3,1)
      ANA(3,8) = +ANA(3,5) -M2
      ANB(3,8) = +ANB(3,5)
      ANC(3,8) = +ANC(3,5)
      ANA(3,9) = +ANA(3,6) -M2
      ANB(3,9) = +ANB(3,6)
      ANC(3,9) = +ANC(3,6)

      ANG2(1) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(2) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(3) = A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG2(4) = A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(5) = A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG2(6) = A4M1P1(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG2(7) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(8) = A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG2(9) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(10)= A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))

      ANG2(11)= A4P0P2(ANA(1,7),ANB(1,7),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG2(12)= A4P1P0(ANA(2,7),ANB(2,7),ANA(2,7),ANB(2,7),ANC(2,7))
      ANG2(13)= A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(14)= A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(15)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(16)= A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(17)= ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(18)= A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(19)= A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(20)= AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))

      ANG2(21)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG2(22)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG2(23)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(24)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(25)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(26)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(27)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(28)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(29)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(30)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))

      ANG2(31)= ABP1M1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(32)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(33)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(34)= ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(35)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(36)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(37)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(38)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(39)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(40)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))

      ANG2(41)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(42)= A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(43)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(44)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))

      COLO2(9) = LOG(S4**2/MS2/(S4+MS2))

      M2QGH = 0.D0
      M2QGH = M2QGH + N*CO*(TG*UG-M2*S)**(-2)*(S4+MS2)*TWO**(-1) * ( 8*
     +    S*TG*MS2*S4**(-1) + 8*S*M2*MS2*S4**(-1) - 8*S*MS2 + 8*TG*M2*
     +    MS2*S4**(-1) - 16*TG*MS2 + 8*TG**2*MS2*S4**(-1) - 8*M2*MS2 + 
     +    8*MS2*S4 )
     +
      M2QGH = M2QGH + N*CO*(TG*UG-M2*S)**(-1)*(S4+MS2)*TWO**(-1) * ( 
     +     - 4 + 4*S*TG**(-1)*M2*S4**(-1) + 8*S*TG**(-1)*MS2*S4**(-1)
     +     + 4*S*S4**(-1) - 8*TG**(-1)*M2 - 16*TG**(-1)*MS2 + 4*TG*
     +    S4**(-1) + 4*M2*S4**(-1) + 16*MS2*S4**(-1) )
     +
      M2QGH = M2QGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +    *T1**(-1)*M2 - 20*S**(-1)*TG**(-1) - 8*S**(-1)*TG*T1**(-2)*
     +    M2**2*S4**(-2) + 8*S**(-1)*TG*T1**(-1)*M2*S4**(-2) - 4*
     +    S**(-1)*TG*T1**(-1)*S4**(-1) + 16*S**(-1)*TG*(S+TG)**(-2) + 
     +    12*S**(-1)*TG*(S+TG)**(-1)*S4**(-1) + 2*S**(-1)*TG*
     +    (S+U1)**(-1)*S4**(-1) - 4*S**(-1)*TG*(T1+UG)**(-1)*S4**(-1)
     +     - 4*S**(-1)*TG*S4**(-2) - 4*S**(-1)*TG**2*T1**(-2)*M2*
     +    S4**(-2) + 4*S**(-1)*TG**2*T1**(-1)*S4**(-2) - 16*S**(-1)*
     +    TG**2*(S+TG)**(-2)*S4**(-1) + 4*S**(-1)*TG**2*(S+TG)**(-1)*
     +    S4**(-2) + 4*S**(-1)*T1**(-2)*M2 - 4*S**(-1)*T1**(-2)*M2**3*
     +    S4**(-2) - 4*S**(-1)*T1**(-1)*M2*S4**(-1) + 4*S**(-1)*
     +    T1**(-1)*M2**2*S4**(-2) - 32*S**(-1)*U1**(-4)*M2**2*S4 + 32*
     +    S**(-1)*U1**(-3)*M2*S4 + 32*S**(-1)*U1**(-3)*M2**2 - 32*
     +    S**(-1)*U1**(-2)*M2 - 4*S**(-1)*U1**(-1) + 2*S**(-1)*
     +    (S+U1)**(-1) + 2*S**(-1)*(T1+UG)**(-1) + 8*S**(-1)*S4**(-1)
     +     + 4*S*TG**(-2)*(S+U1)**(-1) )
     +
      M2QGH = M2QGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S*TG**(-2)*
     +    S4**(-1) + 16*S*TG**(-1)*M2*(S+U1)**(-1)*(S+UG)**(-2) - 8*S*
     +    TG**(-1)*M2*(S+U1)**(-1)*(S+UG)**(-1)*S4**(-1) - 16*S*
     +    TG**(-1)*M2*(S+UG)**(-2)*S4**(-1) - 16*S*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*(S+UG)**(-2)*S4**(-1) + 8*S*TG**(-1)*M2**2*
     +    (S+UG)**(-2)*S4**(-2) + 32*S*TG**(-1)*(S+TG)**(-2) + 40*S*
     +    TG**(-1)*(S+TG)**(-1)*S4**(-1) - 12*S*TG**(-1)*(S+U1)**(-1)*
     +    (S+UG)**(-1) - 16*S*TG**(-1)*(S+U1)**(-1)*S4**(-1) + 4*S*
     +    TG**(-1)*(S+UG)**(-1)*S4**(-1) + 2*S*TG**(-1)*(T1+UG)**(-1)*
     +    S4**(-1) + 16*S*TG**(-1)*S4**(-2) - 16*S*TG*(S+U1)**(-1)*
     +    (S+UG)**(-2)*S4**(-1) + 8*S*TG*(S+UG)**(-2)*S4**(-2) - 14*S*
     +    U1**(-1)*(T1+UG)**(-1)*S4**(-1) - 32*S*M2*(S+U1)**(-1)*
     +    (S+UG)**(-2)*S4**(-1) + 16*S*M2*(S+UG)**(-2)*S4**(-2) - 72*S*
     +    (S+TG)**(-2)*S4**(-1) + 4*S*(S+TG)**(-1)*S4**(-2) + 16*S*
     +    (S+U1)**(-1)*(S+UG)**(-2) - 4*S*(S+U1)**(-1)*(S+UG)**(-1)*
     +    S4**(-1) )
     +
      M2QGH = M2QGH + N*CO*(S4+MS2)*TWO**(-1) * ( 20*S*(S+U1)**(-1)*
     +    S4**(-2) - 16*S*(S+UG)**(-2)*S4**(-1) + 18*S*(T1+UG)**(-1)*
     +    S4**(-2) - 32*S**2*TG**(-1)*(S+TG)**(-2)*S4**(-1) - 4*S**2*
     +    TG**(-1)*(S+TG)**(-1)*S4**(-2) - 64*TG**(-2)*U1**(-4)*M2*MS2*
     +    S4**2 - 64*TG**(-2)*U1**(-4)*M2**2*S4**2 + 128*TG**(-2)*
     +    U1**(-3)*M2*MS2*S4 + 128*TG**(-2)*U1**(-3)*M2**2*S4 - 64*
     +    TG**(-2)*U1**(-2)*M2*MS2 - 64*TG**(-2)*U1**(-2)*M2**2 + 2*
     +    TG**(-2)*U1**(-1)*(S+U1)**(-1)*S4**2 - 2*TG**(-2)*U1**(-1)*S4
     +     + 16*TG**(-2)*M2*(S+UG)**(-1) + 16*TG**(-2)*MS2*(S+UG)**(-1)
     +     - 2*TG**(-2)*(S+U1)**(-1)*S4 + 2*TG**(-2) - 12*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) - 8*TG**(-1)*T1**(-1)*MS2*S4**(-1) + 64*
     +    TG**(-1)*U1**(-4)*M2*MS2*S4 - 64*TG**(-1)*U1**(-3)*M2*MS2 + 
     +    64*TG**(-1)*U1**(-3)*M2*S4 - 64*TG**(-1)*U1**(-2)*M2 - 8*
     +    TG**(-1)*U1**(-1)*(S+TG)**(-1)*S4 + 10*TG**(-1)*U1**(-1)*
     +    (S+U1)**(-1)*S4 - 8*TG**(-1)*U1**(-1) + 6*TG**(-1)*M2*
     +    (S+U1)**(-1)*(S+UG)**(-1) )
     +
      M2QGH = M2QGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 6*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4**(-1) - 16*TG**(-1)*MS2*(S+UG)**(-1)*S4**(-1)
     +     - 8*TG**(-1)*(S+TG)**(-1) - 6*TG**(-1)*(S+U1)**(-1)*
     +    (S+UG)**(-1)*S4 + 2*TG**(-1)*(S+U1)**(-1) + 6*TG**(-1)*
     +    (S+UG)**(-1) + 2*TG**(-1)*(T1+UG)**(-1) + 10*TG**(-1)*
     +    S4**(-1) - 8*TG*T1**(-1)*(S+U1)**(-1)*S4**(-1) + 8*TG*
     +    T1**(-1)*S4**(-2) - 32*TG*M2*(S+U1)**(-1)*(S+UG)**(-2)*
     +    S4**(-1) - 56*TG*(S+TG)**(-2)*S4**(-1) + 12*TG*(S+TG)**(-1)*
     +    S4**(-2) + 32*TG*(S+U1)**(-1)*(S+UG)**(-2) - 14*TG*
     +    (S+U1)**(-1)*(S+UG)**(-1)*S4**(-1) + 16*TG*(S+U1)**(-1)*
     +    S4**(-2) - 16*TG*(S+UG)**(-2)*S4**(-1) - 8*TG*(S+UG)**(-1)*
     +    S4**(-2) + 4*TG*(T1+UG)**(-1)*S4**(-2) - 16*TG**2*
     +    (S+U1)**(-1)*(S+UG)**(-2)*S4**(-1) + 8*T1**(-1)*M2*S4**(-2)
     +     + 8*T1**(-1)*(S+U1)**(-1) - 12*T1**(-1)*S4**(-1) - 32*
     +    U1**(-4)*M2*MS2 - 32*U1**(-3)*M2 - 16*U1**(-2)*M2*
     +    (S+TG)**(-1) )
     +
      M2QGH = M2QGH + N*CO*(S4+MS2)*TWO**(-1) * ( 16*U1**(-2)*M2*
     +    S4**(-1) + 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 12*U1**(-1)
     +    *(S+TG)**(-1) - 8*U1**(-1)*(S+U1)**(-1) + 18*U1**(-1)*
     +    (T1+UG)**(-1) - 26*U1**(-1)*S4**(-1) + 32*M2*(S+U1)**(-1)*
     +    (S+UG)**(-2) - 14*M2*(S+U1)**(-1)*(S+UG)**(-1)*S4**(-1) - 16*
     +    M2*(S+UG)**(-2)*S4**(-1) - 8*M2*(S+UG)**(-1)*S4**(-2) - 16*
     +    M2**2*(S+U1)**(-1)*(S+UG)**(-2)*S4**(-1) + 48*(S+TG)**(-2) + 
     +    32*(S+TG)**(-1)*S4**(-1) - 16*(S+U1)**(-1)*(S+UG)**(-2)*S4 + 
     +    20*(S+U1)**(-1)*(S+UG)**(-1) - 8*(S+U1)**(-1)*S4**(-1) + 16*
     +    (S+UG)**(-2) - 6*(S+UG)**(-1)*S4**(-1) - 20*(T1+UG)**(-1)*
     +    S4**(-1) + 10*S4**(-2) )
     +
      M2QGH = M2QGH + N*CO*TWO**(-1) * (  - 20*S**(-1)*TG**(-1)*M2*MS2*
     +    S4**(-1) - 20*S**(-1)*TG**(-1)*M2 + 16*S**(-1)*TG**(-1)*M2**2
     +    *MS2*S4**(-2) + 12*S**(-1)*TG**(-1)*M2**2*S4**(-1) + 4*
     +    S**(-1)*TG**(-1)*MS2 + 8*S**(-1)*TG**(-1)*S4 + 8*S**(-1)*TG*
     +    MS2*S4**(-2) + 4*S**(-1)*TG*S4**(-1) + 16*S**(-1)*M2*MS2*
     +    S4**(-2) + 12*S**(-1)*M2*S4**(-1) - 8*S**(-1)*MS2*S4**(-1) - 
     +    8*S**(-1) + 16*S*TG**(-1)*MS2*S4**(-2) + 8*S*TG**(-1)*
     +    S4**(-1) + 32*TG**(-2)*M2*MS2*S4**(-1) + 32*TG**(-2)*M2*
     +    MS2**2*S4**(-2) - 16*TG**(-2)*M2 + 32*TG**(-2)*M2**2*MS2*
     +    S4**(-2) + 32*TG**(-2)*M2**2*S4**(-1) - 16*TG**(-2)*MS2 + 32*
     +    TG**(-1)*U1**(-1)*M2*MS2*S4**(-1) + 32*TG**(-1)*U1**(-1)*M2*
     +    MS2**2*S4**(-2) + 16*TG**(-1)*U1**(-1)*M2**2*MS2*S4**(-2) + 8
     +    *TG**(-1)*U1**(-1)*M2**2*S4**(-1) - 12*TG**(-1)*U1**(-1)*MS2
     +     - 4*TG**(-1)*U1**(-1)*S4 + 32*TG**(-1)*M2*MS2*S4**(-2) + 44*
     +    TG**(-1)*M2*S4**(-1) + 16*TG**(-1)*MS2*S4**(-1) - 12*TG**(-1)
     +     + 16*U1**(-2)*M2*MS2*S4**(-1) )
     +
      M2QGH = M2QGH + N*CO*TWO**(-1) * ( 16*U1**(-2)*M2*MS2**2*S4**(-2)
     +     + 8*U1**(-2)*M2 + 16*U1**(-1)*M2*MS2*S4**(-2) + 8*U1**(-1)*
     +    M2*S4**(-1) + 8*U1**(-1)*MS2*S4**(-1) + 16*MS2*S4**(-2) + 16*
     +    S4**(-1) )
     +
      M2QGH = M2QGH + N*CK*(TG*UG-M2*S)**(-2)*(S4+MS2)*TWO**(-1) * ( 
     +     - 24*S*TG*MS2*S4**(-1) - 24*S*M2*MS2*S4**(-1) + 24*S*MS2 - 
     +    24*TG*M2*MS2*S4**(-1) + 48*TG*MS2 - 24*TG**2*MS2*S4**(-1) + 
     +    24*M2*MS2 - 24*MS2*S4 )
     +
      M2QGH = M2QGH + N*CK*(TG*UG-M2*S)**(-1)*(S4+MS2)*TWO**(-1) * ( 
     +     - 4 + 4*S**(-1)*TG*M2*S4**(-1) - 8*S**(-1)*TG + 4*S**(-1)*
     +    TG**2*S4**(-1) - 4*S**(-1)*M2 + 4*S**(-1)*S4 - 4*S*TG**(-1)*
     +    M2*S4**(-1) + 8*S*TG**(-1)*MS2*S4**(-1) + 8*S*TG*(S+TG)**(-1)
     +    *S4**(-1) + 8*S*M2*(S+TG)**(-1)*S4**(-1) - 4*S*S4**(-1) + 16*
     +    TG**(-1)*M2 + 16*TG**(-1)*MS2 + 8*TG*M2*(S+TG)**(-1)*S4**(-1)
     +     - 8*TG*(S+TG)**(-1) + 8*TG**2*(S+TG)**(-1)*S4**(-1) - 16*MS2
     +    *S4**(-1) + 4*(S+TG)**(-1)*S4 )
     +
      M2QGH = M2QGH + N*CK*(S4+MS2)*TWO**(-1) * ( 12*S**(-1)*TG**(-1)*
     +    T1**(-1)*M2 - 4*S**(-1)*TG**(-1)*M2*(S+UG)**(-1) + 4*S**(-1)*
     +    TG**(-1)*(S+UG)**(-1)*S4 + 24*S**(-1)*TG*T1**(-2)*M2**2*
     +    S4**(-2) - 4*S**(-1)*TG*T1**(-1)*M2*(T1+UG)**(-1)*S4**(-1) - 
     +    32*S**(-1)*TG*T1**(-1)*M2*S4**(-2) + 12*S**(-1)*TG*T1**(-1)*
     +    S4**(-1) - 8*S**(-1)*TG*M2*(S+UG)**(-2)*S4**(-1) + 4*S**(-1)*
     +    TG*M2**2*(S+UG)**(-2)*S4**(-2) - 8*S**(-1)*TG*(S+TG)**(-2) - 
     +    16*S**(-1)*TG*(S+TG)**(-1)*S4**(-1) + 4*S**(-1)*TG*
     +    (S+UG)**(-2) + 4*S**(-1)*TG*(S+UG)**(-1)*S4**(-1) + 4*S**(-1)
     +    *TG*(T1+UG)**(-1)*S4**(-1) + 12*S**(-1)*TG**2*T1**(-2)*M2*
     +    S4**(-2) - 12*S**(-1)*TG**2*T1**(-1)*S4**(-2) + 8*S**(-1)*
     +    TG**2*M2*(S+UG)**(-2)*S4**(-2) + 8*S**(-1)*TG**2*(S+TG)**(-2)
     +    *S4**(-1) + 4*S**(-1)*TG**2*(S+TG)**(-1)*S4**(-2) - 8*S**(-1)
     +    *TG**2*(S+UG)**(-2)*S4**(-1) + 4*S**(-1)*TG**3*(S+UG)**(-2)*
     +    S4**(-2) - 12*S**(-1)*T1**(-2)*M2 + 12*S**(-1)*T1**(-2)*M2**3
     +    *S4**(-2) )
     +
      M2QGH = M2QGH + N*CK*(S4+MS2)*TWO**(-1) * ( 4*S**(-1)*T1**(-1)*M2
     +    *(T1+UG)**(-1) + 8*S**(-1)*T1**(-1)*M2*S4**(-1) - 4*S**(-1)*
     +    T1**(-1)*M2**2*(T1+UG)**(-1)*S4**(-1) - 20*S**(-1)*T1**(-1)*
     +    M2**2*S4**(-2) + 32*S**(-1)*U1**(-4)*M2**2*S4 - 32*S**(-1)*
     +    U1**(-3)*M2*S4 - 32*S**(-1)*U1**(-3)*M2**2 + 32*S**(-1)*
     +    U1**(-2)*M2 + 4*S**(-1)*U1**(-1)*(S+TG)**(-1)*S4 - 4*S**(-1)*
     +    U1**(-1) + 4*S**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 4*S**(-1)*M2*
     +    (T1+UG)**(-1)*S4**(-1) + 8*S**(-1)*M2*S4**(-2) + 8*S**(-1)*
     +    (S+TG)**(-1) - 8*S**(-1)*(S+UG)**(-1) - 4*S**(-1)*
     +    (T1+UG)**(-1) + 8*S**(-1)*S4**(-1) + 8*S*TG**(-1)*M2*
     +    (S+UG)**(-2)*S4**(-1) - 8*S*TG**(-1)*M2*(S+UG)**(-1)*S4**(-2)
     +     + 16*S*TG**(-1)*(S+UG)**(-1)*S4**(-1) - 8*S*TG**(-1)*
     +    S4**(-2) - 4*S*(S+TG)**(-1)*S4**(-2) + 8*S*(S+UG)**(-2)*
     +    S4**(-1) - 8*S*(S+UG)**(-1)*S4**(-2) - 16*TG**(-2)*M2*
     +    (S+UG)**(-1) - 16*TG**(-2)*MS2*(S+UG)**(-1) + 12*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) )
     +
      M2QGH = M2QGH + N*CK*(S4+MS2)*TWO**(-1) * ( 24*TG**(-1)*T1**(-1)*
     +    MS2*S4**(-1) - 4*TG**(-1)*U1**(-1)*M2*(S+UG)**(-2)*S4 + 4*
     +    TG**(-1)*U1**(-1)*M2*(S+UG)**(-1) + 4*TG**(-1)*U1**(-1)*M2**2
     +    *(S+UG)**(-2) + 4*TG**(-1)*M2*(S+UG)**(-2) + 12*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4**(-1) - 4*TG**(-1)*M2**2*(S+UG)**(-2)*
     +    S4**(-1) + 16*TG**(-1)*MS2*(S+UG)**(-1)*S4**(-1) - 16*
     +    TG**(-1)*(S+UG)**(-1) + 24*TG*M2*(S+UG)**(-2)*S4**(-2) + 8*TG
     +    *(S+TG)**(-2)*S4**(-1) - 8*TG*(S+UG)**(-2)*S4**(-1) + 8*TG*
     +    (S+UG)**(-1)*S4**(-2) - 4*TG*(T1+UG)**(-1)*S4**(-2) + 12*
     +    TG**2*(S+UG)**(-2)*S4**(-2) - 4*T1**(-1)*M2*(T1+UG)**(-1)*
     +    S4**(-1) + 4*T1**(-1)*S4**(-1) + 32*U1**(-4)*M2*MS2 + 32*
     +    U1**(-3)*M2 + 16*U1**(-2)*M2*(S+TG)**(-1) - 16*U1**(-2)*M2*
     +    S4**(-1) - 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 4*U1**(-1)*
     +    M2*(S+UG)**(-2) - 8*U1**(-1)*(S+TG)**(-1) + 8*U1**(-1)*
     +    S4**(-1) - 12*M2*(S+UG)**(-2)*S4**(-1) + 8*M2*(S+UG)**(-1)*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + N*CK*(S4+MS2)*TWO**(-1) * ( 12*M2**2*(S+UG)**(-2)
     +    *S4**(-2) - 8*(S+TG)**(-2) + 8*(S+TG)**(-1)*S4**(-1) - 4*
     +    (S+UG)**(-2) + 16*(S+UG)**(-1)*S4**(-1) + 4*(T1+UG)**(-1)*
     +    S4**(-1) )
     +
      M2QGH = M2QGH + N*CK*TWO**(-1) * ( 20*S**(-1)*TG**(-1)*M2*MS2*
     +    S4**(-1) + 20*S**(-1)*TG**(-1)*M2 - 16*S**(-1)*TG**(-1)*M2**2
     +    *MS2*S4**(-2) - 12*S**(-1)*TG**(-1)*M2**2*S4**(-1) - 4*
     +    S**(-1)*TG**(-1)*MS2 - 8*S**(-1)*TG**(-1)*S4 - 24*S**(-1)*TG*
     +    MS2*S4**(-2) - 4*S**(-1)*TG*S4**(-1) - 8*S**(-1)*U1**(-1)*M2*
     +    MS2*S4**(-1) - 8*S**(-1)*U1**(-1)*M2 + 32*S**(-1)*U1**(-1)*
     +    M2**2*MS2*S4**(-2) + 8*S**(-1)*U1**(-1)*M2**2*S4**(-1) - 48*
     +    S**(-1)*M2*MS2*S4**(-2) - 4*S**(-1)*M2*S4**(-1) + 32*S**(-1)*
     +    MS2*S4**(-1) + 8*S**(-1) - 16*S*TG**(-1)*MS2*S4**(-2) - 8*S*
     +    TG**(-1)*S4**(-1) - 8*S*U1**(-1)*S4**(-1) - 32*TG**(-2)*M2*
     +    MS2*S4**(-1) - 32*TG**(-2)*M2*MS2**2*S4**(-2) + 16*TG**(-2)*
     +    M2 - 32*TG**(-2)*M2**2*MS2*S4**(-2) - 32*TG**(-2)*M2**2*
     +    S4**(-1) + 16*TG**(-2)*MS2 - 32*TG**(-1)*U1**(-1)*M2*MS2*
     +    S4**(-1) - 32*TG**(-1)*U1**(-1)*M2*MS2**2*S4**(-2) - 16*
     +    TG**(-1)*U1**(-1)*M2**2*MS2*S4**(-2) - 8*TG**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) )
     +
      M2QGH = M2QGH + N*CK*TWO**(-1) * ( 12*TG**(-1)*U1**(-1)*MS2 + 4*
     +    TG**(-1)*U1**(-1)*S4 - 32*TG**(-1)*M2*MS2*S4**(-2) - 44*
     +    TG**(-1)*M2*S4**(-1) - 16*TG**(-1)*MS2*S4**(-1) + 12*TG**(-1)
     +     - 48*U1**(-2)*M2*MS2*S4**(-1) - 48*U1**(-2)*M2*MS2**2*
     +    S4**(-2) - 24*U1**(-2)*M2 - 48*U1**(-1)*M2*MS2*S4**(-2) - 8*
     +    U1**(-1)*M2*S4**(-1) + 8*U1**(-1) - 16*MS2*S4**(-2) - 24*
     +    S4**(-1) )
     +
      M2QGH = M2QGH + CQED*(TG*UG-M2*S)**(-2)*(S4+MS2)*TWO**(-1) * ( 8*
     +    S*TG*MS2*S4**(-1) + 8*S*M2*MS2*S4**(-1) - 8*S*MS2 + 8*TG*M2*
     +    MS2*S4**(-1) - 16*TG*MS2 + 8*TG**2*MS2*S4**(-1) - 8*M2*MS2 + 
     +    8*MS2*S4 )
     +
      M2QGH = M2QGH + CQED*(TG*UG-M2*S)**(-1)*(S4+MS2)*TWO**(-1) * ( 4
     +     - 4*S**(-1)*TG*M2*S4**(-1) + 8*S**(-1)*TG - 4*S**(-1)*TG**2*
     +    S4**(-1) + 4*S**(-1)*M2 - 4*S**(-1)*S4 - 8*S*TG**(-1)*MS2*
     +    S4**(-1) - 4*TG**(-1)*M2 - 4*TG*S4**(-1) - 4*M2*S4**(-1) )
     +
      M2QGH = M2QGH + CQED*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +    *T1**(-1)*M2 - 8*S**(-1)*TG*T1**(-2)*M2**2*S4**(-2) + 16*
     +    S**(-1)*TG*T1**(-1)*M2*S4**(-2) - 4*S**(-1)*TG*T1**(-1)*
     +    S4**(-1) + 4*S**(-1)*TG*M2*(S+UG)**(-1)*S4**(-2) - 4*S**(-1)*
     +    TG*(S+UG)**(-1)*S4**(-1) - 4*S**(-1)*TG**2*T1**(-2)*M2*
     +    S4**(-2) + 4*S**(-1)*TG**2*T1**(-1)*S4**(-2) + 4*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4**(-2) + 4*S**(-1)*T1**(-2)*M2 - 4*
     +    S**(-1)*T1**(-2)*M2**3*S4**(-2) - 4*S**(-1)*T1**(-1)*M2*
     +    S4**(-1) + 12*S**(-1)*T1**(-1)*M2**2*S4**(-2) - 8*S**(-1)*M2*
     +    S4**(-2) + 4*S*U1**(-1)*(S+UG)**(-1)*S4**(-1) - 8*TG**(-1)*
     +    T1**(-1)*MS2*S4**(-1) - 4*U1**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 4*U1**(-1)*S4**(-1) + 4*(S+UG)**(-1)*S4**(-1) )
     +
      M2QGH = M2QGH + CQED*TWO**(-1) * ( 8*S**(-1)*TG*MS2*S4**(-2) + 4*
     +    S**(-1)*U1**(-1)*M2*MS2*S4**(-1) + 4*S**(-1)*U1**(-1)*M2 - 16
     +    *S**(-1)*U1**(-1)*M2**2*MS2*S4**(-2) - 4*S**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) + 16*S**(-1)*M2*MS2*S4**(-2) - 4*S**(-1)*M2*
     +    S4**(-1) - 12*S**(-1)*MS2*S4**(-1) + 4*S*U1**(-1)*S4**(-1) + 
     +    16*U1**(-2)*M2*MS2*S4**(-1) + 16*U1**(-2)*M2*MS2**2*S4**(-2)
     +     + 8*U1**(-2)*M2 + 16*U1**(-1)*M2*MS2*S4**(-2) - 4*U1**(-1)*
     +    MS2*S4**(-1) - 4*U1**(-1) + 4*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(1)*N*CO*TWO**(-1) * ( 32*M2*MS2**2 + 64*
     +    M2**2*MS2 + 32*M2**3 )
     +
      M2QGH = M2QGH + ANG2(2)*N*CO*TWO**(-1) * ( 32*TG**(-1)*M2*MS2**2
     +     + 64*TG**(-1)*M2**2*MS2 + 32*TG**(-1)*M2**3 - 8*TG*M2 + 32*
     +    M2*MS2 + 32*M2**2 )
     +
      M2QGH = M2QGH + ANG2(3)*N*CO*TWO**(-1) * (  - 32*S**(-1)*M2*
     +    MS2**2 - 32*S**(-1)*M2**2*MS2 + 16*S*M2 + 16*S*MS2 + 32*M2*
     +    MS2 + 32*M2**2 )
     +
      M2QGH = M2QGH + ANG2(4)*N*CO*TWO**(-1) * (  - 32*S**(-2)*M2*
     +    MS2**2 - 48*S**(-2)*M2**2*MS2 - 16*S**(-2)*M2**3 - 32*S**(-1)
     +    *TG**(-1)*M2*MS2*S4 - 16*S**(-1)*TG**(-1)*M2*MS2**2 - 16*
     +    S**(-1)*TG**(-1)*M2**2*MS2 - 32*S**(-1)*TG**(-1)*M2**2*S4 + 
     +    18*S**(-1)*TG*M2*MS2*(S+U1)**(-1) + 2*S**(-1)*TG*M2*MS2*
     +    (T1+UG)**(-1) + 2*S**(-1)*TG*M2*MS2*S4**(-1) + 6*S**(-1)*TG*
     +    M2*(S+U1)**(-1)*S4 - 2*S**(-1)*TG*M2 + 6*S**(-1)*TG*M2**2*
     +    (S+U1)**(-1) + 2*S**(-1)*TG*M2**2*(T1+UG)**(-1) + 2*S**(-1)*
     +    TG*M2**2*S4**(-1) - 18*S**(-1)*TG*MS2*(S+U1)**(-1)*S4 + 2*
     +    S**(-1)*TG*MS2 - 12*S**(-1)*TG*(S+U1)**(-1)*S4**2 + 8*S**(-1)
     +    *TG*S4 - 2*S**(-1)*TG**2*M2*(S+U1)**(-1) + 10*S**(-1)*TG**2*
     +    MS2*(S+U1)**(-1) + 12*S**(-1)*TG**2*(S+U1)**(-1)*S4 - 4*
     +    S**(-1)*TG**2 - 4*S**(-1)*TG**3*(S+U1)**(-1) + 2*S**(-1)*
     +    U1**(-1)*M2*MS2*S4 - 16*S**(-1)*U1**(-1)*M2*MS2**2 - 16*
     +    S**(-1)*U1**(-1)*M2**2*MS2 + 2*S**(-1)*U1**(-1)*M2**2*S4 - 4*
     +    S**(-1)*U1**(-1)*M2**3 )
     +
      M2QGH = M2QGH + ANG2(4)*N*CO*TWO**(-1) * (  - 16*S**(-1)*M2*MS2*
     +    (S+U1)**(-1)*S4 - 4*S**(-1)*M2*MS2*(T1+UG)**(-1)*S4 + 18*
     +    S**(-1)*M2*MS2 - 4*S**(-1)*M2*(S+U1)**(-1)*S4**2 + 4*S**(-1)*
     +    M2*S4 + 8*S**(-1)*M2**2*MS2*(S+U1)**(-1) + 8*S**(-1)*M2**2*
     +    MS2*S4**(-1) - 4*S**(-1)*M2**2*(S+U1)**(-1)*S4 + 18*S**(-1)*
     +    M2**2 + 4*S**(-1)*M2**3*(S+U1)**(-1) - 4*S**(-1)*M2**3*
     +    (T1+UG)**(-1) + 4*S**(-1)*M2**3*S4**(-1) + 8*S**(-1)*MS2*
     +    (S+U1)**(-1)*S4**2 + 4*S**(-1)*(S+U1)**(-1)*S4**3 - 4*S**(-1)
     +    *S4**2 - 8*S*TG**(-1)*M2*MS2*(S+U1)**(-1) - 4*S*TG**(-1)*M2*
     +    MS2*(T1+UG)**(-1) - 8*S*TG**(-1)*M2*MS2*S4**(-1) + 8*S*
     +    TG**(-1)*M2*(S+U1)**(-1)*S4 - 2*S*TG**(-1)*M2 - 8*S*TG**(-1)*
     +    M2**2*(S+U1)**(-1) - 4*S*TG**(-1)*M2**2*(T1+UG)**(-1) - 8*S*
     +    TG**(-1)*M2**2*S4**(-1) + 8*S*TG**(-1)*MS2*(S+U1)**(-1)*S4 - 
     +    2*S*TG**(-1)*MS2 - 6*S*TG*(S+U1)**(-1) + 6*S*U1**(-1)*MS2 - 
     +    14*S*M2*(S+U1)**(-1) - 5*S*M2*(T1+UG)**(-1) - 2*S*M2*S4**(-1)
     +     - 8*S*MS2*(S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(4)*N*CO*TWO**(-1) * ( 6*S*MS2*(T1+UG)**(-1)
     +     + 4*S*MS2*S4**(-1) + 6*S*(S+U1)**(-1)*S4 + S*(T1+UG)**(-1)*
     +    S4 - 7*S - 2*S**2*TG**(-1)*M2*(T1+UG)**(-1) - 2*S**2*TG**(-1)
     +    *MS2*(T1+UG)**(-1) - 2*S**2*U1**(-1) - 2*S**2*(S+U1)**(-1) - 
     +    3*S**2*(T1+UG)**(-1) - 2*S**2*S4**(-1) - 8*TG**(-1)*U1**(-1)*
     +    M2*MS2*S4 - 8*TG**(-1)*U1**(-1)*M2**2*S4 + 16*TG**(-1)*M2*MS2
     +    *(S+U1)**(-1)*S4 - 4*TG**(-1)*M2*MS2*(T1+UG)**(-1)*S4 + 28*
     +    TG**(-1)*M2*MS2 + 16*TG**(-1)*M2*MS2**2*S4**(-1) - 8*TG**(-1)
     +    *M2*(S+U1)**(-1)*S4**2 + 2*TG**(-1)*M2*(T1+UG)**(-1)*S4**2 - 
     +    2*TG**(-1)*M2*S4 - 8*TG**(-1)*M2**2*MS2*(S+U1)**(-1) + 24*
     +    TG**(-1)*M2**2*MS2*S4**(-1) + 16*TG**(-1)*M2**2*(S+U1)**(-1)*
     +    S4 - 4*TG**(-1)*M2**2*(T1+UG)**(-1)*S4 + 28*TG**(-1)*M2**2 - 
     +    8*TG**(-1)*M2**3*(S+U1)**(-1) + 8*TG**(-1)*M2**3*S4**(-1) - 8
     +    *TG**(-1)*MS2*(S+U1)**(-1)*S4**2 + 2*TG**(-1)*MS2*
     +    (T1+UG)**(-1)*S4**2 - 2*TG**(-1)*MS2*S4 - 14*TG*M2*
     +    (S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(4)*N*CO*TWO**(-1) * ( 2*TG*M2*(T1+UG)**(-1)
     +     - 2*TG*M2*S4**(-1) + 2*TG*MS2*(S+U1)**(-1) + 2*TG*MS2*
     +    (T1+UG)**(-1) - 2*TG*MS2*S4**(-1) + 16*TG*(S+U1)**(-1)*S4 - 8
     +    *TG - 8*TG**2*(S+U1)**(-1) + 14*U1**(-1)*M2*MS2 + 2*U1**(-1)*
     +    M2*S4 + 6*U1**(-1)*M2**2 + 2*U1**(-1)*MS2*S4 - 8*M2*MS2*
     +    (S+U1)**(-1) + 12*M2*MS2*(T1+UG)**(-1) + 8*M2*MS2*S4**(-1) + 
     +    24*M2*(S+U1)**(-1)*S4 - M2*(T1+UG)**(-1)*S4 + 15*M2 - 16*
     +    M2**2*(S+U1)**(-1) + 8*MS2*(S+U1)**(-1)*S4 - 6*MS2*
     +    (T1+UG)**(-1)*S4 + 28*MS2 - 8*(S+U1)**(-1)*S4**2 + 8*S4 )
     +
      M2QGH = M2QGH + ANG2(5)*N*CO*TWO**(-1) * ( 8*S**(-1)*M2 + 8*
     +    S**(-1)*MS2 )
     +
      M2QGH = M2QGH + ANG2(5)*N*CK*TWO**(-1) * (  - 8*S**(-1)*M2 - 8*
     +    S**(-1)*MS2 )
     +
      M2QGH = M2QGH + ANG2(6)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG**(-1)*M2*
     +    MS2*S4**(-1) - 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*TG**(-1)*MS2 - 6*S**(-1)*M2*
     +    (T1+UG)**(-1) + 2*S**(-1)*M2*S4**(-1) - 4*S**(-1)*MS2*
     +    (T1+UG)**(-1) + 4*S**(-1)*MS2*S4**(-1) + 2*S**(-1)*
     +    (T1+UG)**(-1)*S4 - 2*S**(-1) - 2*TG**(-1)*M2*S4**(-1) - 2*
     +    TG**(-1)*MS2*S4**(-1) - 2*(T1+UG)**(-1) - 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(6)*N*CK*TWO**(-1) * ( 2*S**(-1)*U1**(-1)*M2*
     +    MS2*S4**(-1) - 2*S**(-1)*U1**(-1)*M2 + 2*S**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*U1**(-1)*MS2 + 6*S**(-1)*M2*
     +    (T1+UG)**(-1) - 2*S**(-1)*M2*S4**(-1) + 4*S**(-1)*MS2*
     +    (T1+UG)**(-1) - 4*S**(-1)*MS2*S4**(-1) - 2*S**(-1)*
     +    (T1+UG)**(-1)*S4 + 2*S**(-1) + 2*S*U1**(-1)*S4**(-1) - 4*
     +    U1**(-1)*M2*S4**(-1) - 6*U1**(-1)*MS2*S4**(-1) - 2*U1**(-1)
     +     + 2*(T1+UG)**(-1) + 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(7)*N*CO*TWO**(-1) * ( 16*M2*MS2**2 + 16*
     +    M2**2*MS2 )
     +
      M2QGH = M2QGH + ANG2(7)*N*CK*TWO**(-1) * (  - 16*M2*MS2**2 - 16*
     +    M2**2*MS2 )
     +
      M2QGH = M2QGH + ANG2(8)*N*CO*TWO**(-1) * (  - 32*S**(-1)*M2*
     +    MS2**2 - 48*S**(-1)*M2**2*MS2 - 16*S**(-1)*M2**3 + 16*M2*MS2
     +     + 16*M2**2 )
     +
      M2QGH = M2QGH + ANG2(8)*N*CK*TWO**(-1) * ( 16*S**(-1)*M2**2*MS2
     +     + 16*S**(-1)*M2**3 - 16*M2*MS2 - 16*M2**2 )
     +
      M2QGH = M2QGH + ANG2(9)*N*CO*TWO**(-1) * ( 8*S*MS2 + 8*TG*MS2 + 
     +    16*U1**(-1)*M2*MS2**2 + 16*M2*MS2 - 8*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(9)*N*CK*TWO**(-1) * (  - 8*S*MS2 - 8*TG*MS2
     +     - 16*U1**(-1)*M2*MS2**2 - 16*M2*MS2 + 8*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(10)*N*CO*TWO**(-1) * (  - 32*S**(-2)*M2*
     +    MS2**2 - 48*S**(-2)*M2**2*MS2 - 16*S**(-2)*M2**3 - 8*S**(-1)*
     +    TG**(-1)*M2*MS2*S4 - 16*S**(-1)*TG**(-1)*M2*MS2**2 - 24*
     +    S**(-1)*TG**(-1)*M2**2*MS2 - 8*S**(-1)*TG**(-1)*M2**2*S4 - 8*
     +    S**(-1)*TG**(-1)*M2**3 + 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) + 6
     +    *S**(-1)*TG*M2**2*(S+TG)**(-1) + 2*S**(-1)*TG*MS2*
     +    (S+TG)**(-1)*S4 - 6*S**(-1)*TG*MS2 + 2*S**(-1)*TG*
     +    (S+TG)**(-1)*S4**2 + 2*S**(-1)*TG*S4 - 2*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) - 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 + 4*S**(-1)*
     +    U1**(-1)*M2*MS2*S4 - 16*S**(-1)*U1**(-1)*M2*MS2**2 - 8*
     +    S**(-1)*U1**(-1)*M2**2*MS2 + 4*S**(-1)*U1**(-1)*M2**2*S4 - 8*
     +    S**(-1)*M2*MS2*(S+TG)**(-1)*S4 - 24*S**(-1)*M2*MS2 + 4*
     +    S**(-1)*M2*S4 - 8*S**(-1)*M2**2*(S+TG)**(-1)*S4 - 18*S**(-1)*
     +    M2**2 + 6*S**(-1)*MS2*S4 - 2*S**(-1)*S4**2 + 6*S*TG**(-1)*MS2
     +     + 4*S*TG**(-1)*S4 + 2*S*TG*(S+TG)**(-1) - 8*S*U1**(-1)*M2*
     +    MS2*(S+TG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(10)*N*CO*TWO**(-1) * (  - 4*S*M2*
     +    (S+TG)**(-1) - 4*S - 2*S**2*TG**(-1) + 2*TG**(-1)*U1**(-1)*M2
     +    *MS2*S4 - 8*TG**(-1)*U1**(-1)*M2**2*MS2 + 2*TG**(-1)*U1**(-1)
     +    *M2**2*S4 - 4*TG**(-1)*U1**(-1)*M2**3 + 14*TG**(-1)*M2*MS2 + 
     +    2*TG**(-1)*M2*S4 + 6*TG**(-1)*M2**2 - 4*TG**(-1)*MS2*S4 - 2*
     +    TG**(-1)*S4**2 - 2*TG*MS2*(S+TG)**(-1) - 4*TG*(S+TG)**(-1)*S4
     +     - 2*TG + 2*TG**2*(S+TG)**(-1) - 8*U1**(-2)*M2**2*MS2 + 4*
     +    U1**(-1)*M2*MS2*(S+TG)**(-1)*S4 + 20*U1**(-1)*M2*MS2 - 2*
     +    U1**(-1)*M2**2*(S+TG)**(-1)*S4 - 2*U1**(-1)*M2**2 + 8*M2*MS2*
     +    (S+TG)**(-1) + 4*M2*(S+TG)**(-1)*S4 + 16*M2 + 12*M2**2*
     +    (S+TG)**(-1) + 2*MS2*(S+TG)**(-1)*S4 + 12*MS2 + 6*S4 )
     +
      M2QGH = M2QGH + ANG2(10)*N*CK*TWO**(-1) * ( 2*S**(-1)*TG*M2*MS2*
     +    (T1+UG)**(-1) + 2*S**(-1)*TG*M2*MS2*S4**(-1) + 4*S**(-1)*TG*
     +    M2 + 2*S**(-1)*TG*M2**2*(T1+UG)**(-1) + 2*S**(-1)*TG*M2**2*
     +    S4**(-1) - 4*S**(-1)*U1**(-1)*M2*MS2*S4 + 8*S**(-1)*U1**(-1)*
     +    M2**2*MS2 - 4*S**(-1)*U1**(-1)*M2**2*S4 - 4*S**(-1)*M2*MS2*
     +    (T1+UG)**(-1)*S4 + 24*S**(-1)*M2*MS2 - 8*S**(-1)*M2*S4 + 8*
     +    S**(-1)*M2**2*MS2*S4**(-1) + 36*S**(-1)*M2**2 - 4*S**(-1)*
     +    M2**3*(T1+UG)**(-1) + 4*S**(-1)*M2**3*S4**(-1) + 2*S*TG*
     +    S4**(-1) - 4*S*U1**(-1)*M2*MS2*(T1+UG)**(-1) + 8*S*U1**(-1)*
     +    M2*MS2*S4**(-1) - 4*S*M2*(T1+UG)**(-1) - 6*S*MS2*S4**(-1) - 6
     +    *S + 2*S**2*S4**(-1) - 2*TG*M2*(T1+UG)**(-1) - 4*TG*M2*
     +    S4**(-1) - 2*TG*MS2*(T1+UG)**(-1) - 6*TG*MS2*S4**(-1) - 4*TG
     +     + 8*U1**(-2)*M2**2*MS2 + 2*U1**(-1)*M2*MS2*(T1+UG)**(-1)*S4
     +     - 28*U1**(-1)*M2*MS2 - 16*U1**(-1)*M2*MS2**2*S4**(-1) - 8*
     +    U1**(-1)*M2**2*MS2*S4**(-1) - 2*U1**(-1)*M2**2*(T1+UG)**(-1)*
     +    S4 )
     +
      M2QGH = M2QGH + ANG2(10)*N*CK*TWO**(-1) * ( 4*U1**(-1)*M2**2 + 4*
     +    U1**(-1)*M2**3*(T1+UG)**(-1) + 16*M2*MS2*(T1+UG)**(-1) - 6*M2
     +    *MS2*S4**(-1) + 2*M2*(T1+UG)**(-1)*S4 - 14*M2 + 8*M2**2*
     +    (T1+UG)**(-1) - 6*M2**2*S4**(-1) - 2*MS2*(T1+UG)**(-1)*S4 - 6
     +    *MS2 + 4*S4 )
     +
      M2QGH = M2QGH + ANG2(11)*N*CO*TWO**(-1) * ( 16*S**(-1)*M2*MS2 + 
     +    16*S**(-1)*M2**2 + 16*M2 + 16*MS2 )
     +
      M2QGH = M2QGH + ANG2(11)*N*CK*TWO**(-1) * (  - 16*S**(-1)*M2*MS2
     +     - 16*S**(-1)*M2**2 )
     +
      M2QGH = M2QGH + ANG2(12)*N*CO*TWO**(-1) * ( 3 - 8*S**(-1)*
     +    TG**(-1)*M2*MS2 - 6*S**(-1)*TG**(-1)*M2*S4 + 8*S**(-1)*
     +    TG**(-1)*M2**2*MS2*S4**(-1) + 4*S**(-1)*TG**(-1)*M2**3*
     +    S4**(-1) + 2*S**(-1)*TG**(-1)*S4**2 + 4*S**(-1)*TG*M2*
     +    (T1+UG)**(-1) + 2*S**(-1)*TG*M2*S4**(-1) + 4*S**(-1)*TG*MS2*
     +    (T1+UG)**(-1) + 4*S**(-1)*TG*MS2*S4**(-1) - S**(-1)*TG*
     +    (T1+UG)**(-1)*S4 + 3*S**(-1)*TG - 2*S**(-1)*M2*MS2*
     +    (T1+UG)**(-1) + 12*S**(-1)*M2*MS2*S4**(-1) + 3*S**(-1)*M2*
     +    (T1+UG)**(-1)*S4 + S**(-1)*M2 - 6*S**(-1)*M2**2*(T1+UG)**(-1)
     +     + 4*S**(-1)*M2**2*S4**(-1) - 2*S**(-1)*MS2*(T1+UG)**(-1)*S4
     +     - 2*S**(-1)*MS2 - S**(-1)*(T1+UG)**(-1)*S4**2 - 3*S**(-1)*S4
     +     - 2*S*TG**(-1)*M2*S4**(-1) + 4*S*TG**(-1)*MS2*S4**(-1) + 6*S
     +    *TG**(-1) - 5*S*(T1+UG)**(-1) - 8*S*S4**(-1) - 2*S**2*
     +    TG**(-1)*S4**(-1) + 8*TG**(-1)*M2*MS2*S4**(-1) + 4*TG**(-1)*
     +    M2 - 8*TG**(-1)*MS2 - 6*TG**(-1)*S4 + 3*TG*(T1+UG)**(-1) - 2*
     +    TG*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(12)*N*CO*TWO**(-1) * (  - 11*M2*
     +    (T1+UG)**(-1) + 4*M2*S4**(-1) - 2*MS2*(T1+UG)**(-1) + 20*MS2*
     +    S4**(-1) + 6*(T1+UG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(12)*N*CK*TWO**(-1) * (  - 4*S**(-1)*TG*M2*
     +    (T1+UG)**(-1) - 2*S**(-1)*TG*M2*S4**(-1) - 4*S**(-1)*TG*MS2*
     +    (T1+UG)**(-1) - 4*S**(-1)*TG*MS2*S4**(-1) - 4*S**(-1)*TG*
     +    (T1+UG)**(-1)*S4 + 2*S**(-1)*TG - 12*S**(-1)*U1**(-1)*M2*MS2
     +     + 2*S**(-1)*U1**(-1)*M2*S4 + 8*S**(-1)*U1**(-1)*M2**2*MS2*
     +    S4**(-1) - 8*S**(-1)*U1**(-1)*M2**2 + 4*S**(-1)*U1**(-1)*
     +    M2**3*S4**(-1) + 2*S**(-1)*U1**(-1)*MS2*S4 + 2*S**(-1)*M2*MS2
     +    *(T1+UG)**(-1) - 12*S**(-1)*M2*MS2*S4**(-1) - 6*S**(-1)*M2*
     +    (T1+UG)**(-1)*S4 + 4*S**(-1)*M2 + 6*S**(-1)*M2**2*
     +    (T1+UG)**(-1) - 4*S**(-1)*M2**2*S4**(-1) + 2*S**(-1)*MS2*
     +    (T1+UG)**(-1)*S4 + 2*S**(-1)*MS2 + 4*S**(-1)*(T1+UG)**(-1)*
     +    S4**2 - 2*S**(-1)*S4 - 6*S*U1**(-1)*MS2*S4**(-1) - 2*S*
     +    U1**(-1) + 4*S*(T1+UG)**(-1) + 4*S*S4**(-1) + 2*S**2*U1**(-1)
     +    *S4**(-1) + 2*TG*(T1+UG)**(-1) + 2*TG*S4**(-1) - 6*U1**(-1)*
     +    M2*MS2*S4**(-1) - 2*U1**(-1)*M2 - 6*U1**(-1)*M2**2*S4**(-1)
     +     + 4*U1**(-1)*MS2 )
     +
      M2QGH = M2QGH + ANG2(12)*N*CK*TWO**(-1) * ( 2*M2*(T1+UG)**(-1) - 
     +    6*MS2*(T1+UG)**(-1) - 12*MS2*S4**(-1) - 8*(T1+UG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(13)*N*CO*TWO**(-1) * (  - 2*S**(-1)*TG*M2 + 
     +    2*S**(-1)*M2*S4 + 4*M2 )
     +
      M2QGH = M2QGH + ANG2(13)*N*CK*TWO**(-1) * ( 6*S**(-1)*TG*M2 - 6*
     +    S**(-1)*M2*S4 - 4*M2 )
     +
      M2QGH = M2QGH + ANG2(13)*CQED*TWO**(-1) * (  - 2*S**(-1)*TG*M2 + 
     +    2*S**(-1)*M2*S4 )
     +
      M2QGH = M2QGH + ANG2(14)*N*CO*TWO**(-1) * (  - 2 + 4*S**(-1)*TG*
     +    M2*(T1+UG)**(-1) + 4*S**(-1)*TG*MS2*(T1+UG)**(-1) - 2*S**(-1)
     +    *TG*(T1+UG)**(-1)*S4 + 2*S**(-1)*TG - 2*S**(-1)*M2*MS2*
     +    (T1+UG)**(-1) - 4*S**(-1)*M2*(T1+UG)**(-1)*S4 + 2*S**(-1)*M2
     +     - 6*S**(-1)*M2**2*(T1+UG)**(-1) - 6*S**(-1)*MS2*
     +    (T1+UG)**(-1)*S4 + 4*S**(-1)*MS2 + 2*S**(-1)*(T1+UG)**(-1)*
     +    S4**2 - 2*S**(-1)*S4 - 2*S*(T1+UG)**(-1) + 4*TG*(T1+UG)**(-1)
     +     - 4*M2*(T1+UG)**(-1) + 2*MS2*(T1+UG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(14)*N*CK*TWO**(-1) * (  - 8 + 2*S**(-1)*
     +    TG**(-1)*M2*MS2 + 4*S**(-1)*TG**(-1)*M2*S4 - 6*S**(-1)*
     +    TG**(-1)*M2**2 + 4*S**(-1)*TG**(-1)*M2**3*S4**(-1) - 2*
     +    S**(-1)*TG**(-1)*MS2*S4 - 2*S**(-1)*TG**(-1)*S4**2 - 4*
     +    S**(-1)*TG*M2*(T1+UG)**(-1) - 4*S**(-1)*TG*MS2*(T1+UG)**(-1)
     +     - 8*S**(-1)*TG*MS2*S4**(-1) - 6*S**(-1)*U1**(-1)*M2*S4 + 8*
     +    S**(-1)*U1**(-1)*M2**2 - 4*S**(-1)*U1**(-1)*M2**3*S4**(-1) + 
     +    2*S**(-1)*U1**(-1)*S4**2 + 2*S**(-1)*M2*MS2*(T1+UG)**(-1) - 4
     +    *S**(-1)*M2*MS2*S4**(-1) + 4*S**(-1)*M2*(T1+UG)**(-1)*S4 + 2*
     +    S**(-1)*M2 + 6*S**(-1)*M2**2*(T1+UG)**(-1) + 4*S**(-1)*M2**2*
     +    S4**(-1) + 6*S**(-1)*MS2*(T1+UG)**(-1)*S4 + 8*S**(-1)*MS2 - 2
     +    *S**(-1)*(T1+UG)**(-1)*S4**2 + 2*S**(-1)*S4 + 2*S*TG**(-1)*M2
     +    *S4**(-1) + 2*S*TG**(-1)*MS2*S4**(-1) - 4*S*U1**(-1)*M2*
     +    S4**(-1) + 4*S*U1**(-1)*MS2*S4**(-1) - 4*S*U1**(-1) - 2*S*
     +    (T1+UG)**(-1) - 12*TG**(-1)*M2*MS2*S4**(-1) - 8*TG**(-1)*
     +    M2**2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(14)*N*CK*TWO**(-1) * ( 4*TG**(-1)*MS2 + 2*
     +    TG**(-1)*S4 - 2*TG*(T1+UG)**(-1) - 4*TG*S4**(-1) + 8*U1**(-1)
     +    *M2*MS2*S4**(-1) + 8*U1**(-1)*M2 + 8*U1**(-1)*M2**2*S4**(-1)
     +     + 4*U1**(-1)*MS2 + 2*U1**(-1)*S4 - 8*M2*(T1+UG)**(-1) - 10*
     +    M2*S4**(-1) - 10*MS2*(T1+UG)**(-1) + 4*(T1+UG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(14)*CQED*TWO**(-1) * ( 2 + 4*S**(-1)*TG*MS2*
     +    S4**(-1) - 2*S**(-1)*U1**(-1)*M2*MS2 + 8*S**(-1)*U1**(-1)*M2*
     +    S4 - 10*S**(-1)*U1**(-1)*M2**2 + 4*S**(-1)*U1**(-1)*M2**3*
     +    S4**(-1) + 2*S**(-1)*U1**(-1)*MS2*S4 - 2*S**(-1)*U1**(-1)*
     +    S4**2 + 2*S**(-1)*M2*MS2*S4**(-1) - 2*S**(-1)*M2 - 2*S**(-1)*
     +    M2**2*S4**(-1) - 6*S**(-1)*MS2 + 8*S*U1**(-1)*M2*S4**(-1) + 2
     +    *S*U1**(-1)*MS2*S4**(-1) + 4*S*U1**(-1) - 2*S*S4**(-1) - 2*
     +    S**2*U1**(-1)*S4**(-1) + 2*TG*S4**(-1) - 10*U1**(-1)*M2*MS2*
     +    S4**(-1) - 2*U1**(-1)*M2 - 10*U1**(-1)*M2**2*S4**(-1) + 4*
     +    U1**(-1)*MS2 + 6*M2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(15)*N*CO*TWO**(-1) * (  - 8*S*TG**(-1)*M2 + 
     +    32*TG**(-2)*M2*MS2**2 + 64*TG**(-2)*M2**2*MS2 + 32*TG**(-2)*
     +    M2**3 + 32*TG**(-1)*M2*MS2 + 32*TG**(-1)*M2**2 + 8*M2 )
     +
      M2QGH = M2QGH + ANG2(16)*N*CO*TWO**(-1) * (  - 3 - 8*S**(-1)*
     +    TG**(-1)*M2*MS2 - 6*S**(-1)*TG**(-1)*M2*S4 + 8*S**(-1)*
     +    TG**(-1)*M2**2*MS2*S4**(-1) + 4*S**(-1)*TG**(-1)*M2**3*
     +    S4**(-1) + 2*S**(-1)*TG**(-1)*S4**2 + 2*S**(-1)*TG*M2*
     +    S4**(-1) + 2*S**(-1)*TG*MS2*S4**(-1) + 12*S**(-1)*M2*MS2*
     +    S4**(-1) + 2*S**(-1)*M2 + 8*S**(-1)*M2**2*S4**(-1) - 2*
     +    S**(-1)*MS2 - 2*S**(-1)*S4 + 8*S*TG**(-1)*MS2*S4**(-1) + 6*S*
     +    TG**(-1)*(S+U1)**(-1)*S4 - 6*S*TG**(-1) - 4*S*U1**(-1)*M2*
     +    S4**(-1) + 4*S*U1**(-1)*MS2*S4**(-1) - 2*S*U1**(-1) - 2*S*
     +    (S+U1)**(-1) + 4*S*S4**(-1) + 32*TG**(-2)*M2*MS2*(S+U1)**(-1)
     +    *S4 + 32*TG**(-2)*M2*MS2**2*S4**(-1) - 16*TG**(-2)*M2*
     +    (S+U1)**(-1)*S4**2 - 16*TG**(-2)*M2**2*MS2*(S+U1)**(-1) + 48*
     +    TG**(-2)*M2**2*MS2*S4**(-1) + 32*TG**(-2)*M2**2*(S+U1)**(-1)*
     +    S4 - 16*TG**(-2)*M2**3*(S+U1)**(-1) + 16*TG**(-2)*M2**3*
     +    S4**(-1) - 16*TG**(-2)*MS2*(S+U1)**(-1)*S4**2 + 6*TG**(-1)*
     +    U1**(-1)*M2*MS2 )
     +
      M2QGH = M2QGH + ANG2(16)*N*CO*TWO**(-1) * ( 16*TG**(-1)*U1**(-1)*
     +    M2*MS2**2*S4**(-1) + 16*TG**(-1)*U1**(-1)*M2**2*MS2*S4**(-1)
     +     - 2*TG**(-1)*U1**(-1)*M2**2 + 4*TG**(-1)*U1**(-1)*M2**3*
     +    S4**(-1) - 6*TG**(-1)*U1**(-1)*MS2*S4 - 2*TG**(-1)*U1**(-1)*
     +    S4**2 - 36*TG**(-1)*M2*MS2*(S+U1)**(-1) + 36*TG**(-1)*M2*MS2*
     +    S4**(-1) + 47*TG**(-1)*M2*(S+U1)**(-1)*S4 - 3*TG**(-1)*M2 - 
     +    44*TG**(-1)*M2**2*(S+U1)**(-1) + 28*TG**(-1)*M2**2*S4**(-1)
     +     + 36*TG**(-1)*MS2*(S+U1)**(-1)*S4 + 8*TG**(-1)*MS2 - 3*
     +    TG**(-1)*(S+U1)**(-1)*S4**2 + 3*TG**(-1)*S4 - 3*TG*
     +    (S+U1)**(-1) + 8*U1**(-1)*M2*MS2*S4**(-1) + 8*U1**(-1)*M2 + 8
     +    *U1**(-1)*M2**2*S4**(-1) + 6*U1**(-1)*MS2 + 2*U1**(-1)*S4 - 
     +    31*M2*(S+U1)**(-1) + 2*M2*S4**(-1) - 20*MS2*(S+U1)**(-1) + 6*
     +    (S+U1)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(17)*N*CO*TWO**(-1) * ( 8*U1**(-2)*M2*S4**2
     +     - 16*U1**(-1)*M2*S4 + 8*M2 )
     +
      M2QGH = M2QGH + ANG2(17)*N*CK*TWO**(-1) * (  - 8*U1**(-2)*M2*
     +    S4**2 + 16*U1**(-1)*M2*S4 - 8*M2 )
     +
      M2QGH = M2QGH + ANG2(18)*N*CO*TWO**(-1) * ( 8*U1**(-2)*M2*MS2**2
     +     + 8*U1**(-1)*M2*MS2 - 4*MS2 )
     +
      M2QGH = M2QGH + ANG2(18)*N*CK*TWO**(-1) * (  - 24*U1**(-2)*M2*
     +    MS2**2 - 24*U1**(-1)*M2*MS2 + 12*MS2 )
     +
      M2QGH = M2QGH + ANG2(18)*CQED*TWO**(-1) * ( 8*U1**(-2)*M2*MS2**2
     +     + 8*U1**(-1)*M2*MS2 - 4*MS2 )
     +
      M2QGH = M2QGH + ANG2(19)*N*CO*TWO**(-1) * ( 7 - 2*S*U1**(-1)*
     +    (S+TG)**(-1)*S4 + 2*S*U1**(-1) - 3*S*(S+TG)**(-1) - 5*TG*
     +    (S+TG)**(-1) + 10*U1**(-2)*M2*MS2*(S+TG)**(-1)*S4 - 2*
     +    U1**(-2)*M2*MS2 - U1**(-2)*M2*(S+TG)**(-1)*S4**2 + U1**(-2)*
     +    M2*S4 - 2*U1**(-1)*M2*MS2*(S+TG)**(-1) + 12*U1**(-1)*M2*
     +    (S+TG)**(-1)*S4 - 3*U1**(-1)*M2 - 4*U1**(-1)*MS2*(S+TG)**(-1)
     +    *S4 + 4*U1**(-1)*MS2 + U1**(-1)*(S+TG)**(-1)*S4**2 - 3*
     +    U1**(-1)*S4 - 3*M2*(S+TG)**(-1) + 4*MS2*(S+TG)**(-1) - 5*
     +    (S+TG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(19)*N*CK*TWO**(-1) * ( 4 + 4*S**(-1)*TG*M2*
     +    S4**(-1) + 2*S**(-1)*TG*MS2*S4**(-1) + 2*S**(-1)*TG + 8*
     +    S**(-1)*U1**(-1)*M2*MS2 - 4*S**(-1)*U1**(-1)*M2*S4 + 8*
     +    S**(-1)*U1**(-1)*M2**2*MS2*S4**(-1) + 8*S**(-1)*U1**(-1)*
     +    M2**2 - 2*S**(-1)*U1**(-1)*MS2*S4 + 8*S**(-1)*M2*MS2*S4**(-1)
     +     + 8*S**(-1)*M2 + 8*S**(-1)*M2**2*S4**(-1) + 2*S**(-1)*MS2 - 
     +    4*S**(-1)*S4 - 4*S*TG**(-1)*M2*S4**(-1) - 6*S*TG**(-1)*MS2*
     +    S4**(-1) - 4*S*S4**(-1) - 16*TG**(-1)*U1**(-1)*M2*MS2 - 16*
     +    TG**(-1)*U1**(-1)*M2*MS2**2*S4**(-1) - 8*TG**(-1)*U1**(-1)*
     +    M2**2*MS2*S4**(-1) - 8*TG**(-1)*U1**(-1)*M2**2 + 6*TG**(-1)*
     +    U1**(-1)*MS2*S4 + 2*TG**(-1)*U1**(-1)*S4**2 - 8*TG**(-1)*M2*
     +    MS2*S4**(-1) + 4*TG**(-1)*M2 - 8*TG**(-1)*M2**2*S4**(-1) + 2*
     +    TG**(-1)*MS2 - 2*TG**(-1)*S4 - 4*TG*S4**(-1) - 10*U1**(-2)*M2
     +    *MS2*(S+TG)**(-1)*S4 - 14*U1**(-2)*M2*MS2 - 32*U1**(-2)*M2*
     +    MS2**2*S4**(-1) + U1**(-2)*M2*(S+TG)**(-1)*S4**2 - U1**(-2)*
     +    M2*S4 )
     +
      M2QGH = M2QGH + ANG2(19)*N*CK*TWO**(-1) * ( 2*U1**(-1)*M2*MS2*
     +    (S+TG)**(-1) - 16*U1**(-1)*M2*MS2*S4**(-1) - 12*U1**(-1)*M2*
     +    (S+TG)**(-1)*S4 + 3*U1**(-1)*M2 + U1**(-1)*MS2*(S+TG)**(-1)*
     +    S4 - 9*U1**(-1)*MS2 + U1**(-1)*(S+TG)**(-1)*S4**2 - 3*
     +    U1**(-1)*S4 + 3*M2*(S+TG)**(-1) - 8*M2*S4**(-1) - MS2*
     +    (S+TG)**(-1) + 4*MS2*S4**(-1) + 3*(S+TG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(19)*CQED*TWO**(-1) * ( 2 - 8*S**(-1)*TG*M2*
     +    S4**(-1) - 2*S**(-1)*TG*MS2*S4**(-1) + 4*S**(-1)*TG - 2*
     +    S**(-1)*TG**2*S4**(-1) + 8*S**(-1)*U1**(-1)*M2*S4 - 8*S**(-1)
     +    *U1**(-1)*M2**2*MS2*S4**(-1) - 8*S**(-1)*U1**(-1)*M2**2 + 2*
     +    S**(-1)*U1**(-1)*MS2*S4 - 2*S**(-1)*U1**(-1)*S4**2 + 4*
     +    S**(-1)*M2 - 8*S**(-1)*M2**2*S4**(-1) - 2*S**(-1)*MS2 - 2*
     +    S**(-1)*S4 - 2*TG*S4**(-1) + 8*U1**(-2)*M2*MS2 + 16*U1**(-2)*
     +    M2*MS2**2*S4**(-1) + 8*U1**(-1)*M2*MS2*S4**(-1) + 4*U1**(-1)*
     +    MS2 + 2*U1**(-1)*S4 - 4*M2*S4**(-1) - 6*MS2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(20)*N*CO*TWO**(-1) * ( 8*S*TG**(-1)*M2*MS2*
     +    (S+U1)**(-1) - 6*S*TG**(-1)*M2*(S+U1)**(-1)*S4 + 2*S*TG**(-1)
     +    *M2 + 8*S*TG**(-1)*M2**2*(S+U1)**(-1) - 4*S*TG**(-1)*MS2*
     +    (S+U1)**(-1)*S4 + 4*S*TG**(-1)*MS2 + 3*S*TG**(-1)*
     +    (S+U1)**(-1)*S4**2 + S*TG**(-1)*S4 + 8*S*TG*(S+U1)**(-1) + 10
     +    *S*M2*(S+U1)**(-1) + 4*S*MS2*(S+U1)**(-1) - 7*S*(S+U1)**(-1)*
     +    S4 - 2*S**2*TG**(-1)*(S+U1)**(-1)*S4 + 4*S**2*(S+U1)**(-1) + 
     +    4*TG**(-1)*M2*(S+U1)**(-1)*S4**2 - 2*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 2*TG**(-1)*M2**2 - 2*TG**(-1)*(S+U1)**(-1)*
     +    S4**3 + 10*TG*M2*(S+U1)**(-1) - 5*TG*(S+U1)**(-1)*S4 + 4*
     +    TG**2*(S+U1)**(-1) - 10*M2*(S+U1)**(-1)*S4 + 2*M2 + 6*M2**2*
     +    (S+U1)**(-1) + 5*(S+U1)**(-1)*S4**2 + S4 )
     +
      M2QGH = M2QGH + ANG2(20)*N*CK*TWO**(-1) * ( 2*S*TG*(S+TG)**(-1)
     +     - 8*S*U1**(-1)*M2*MS2*(S+TG)**(-1) + 4*S*U1**(-1)*M2*
     +    (S+TG)**(-1)*S4 - 4*S*U1**(-1)*MS2*(S+TG)**(-1)*S4 + 4*S*
     +    U1**(-1)*MS2 + 2*S*U1**(-1)*(S+TG)**(-1)*S4**2 - 2*S*M2*
     +    (S+TG)**(-1) + 4*S*MS2*(S+TG)**(-1) - 3*S*(S+TG)**(-1)*S4 + S
     +     - S**2*(S+TG)**(-1) + 6*TG*M2*(S+TG)**(-1) - 3*TG*
     +    (S+TG)**(-1)*S4 + TG + 3*TG**2*(S+TG)**(-1) - 4*U1**(-1)*M2*
     +    S4 + U1**(-1)*M2**2*(S+TG)**(-1)*S4 - 5*U1**(-1)*M2**2 - 2*
     +    U1**(-1)*S4**2 - 4*M2*(S+TG)**(-1)*S4 + 2*M2 + 3*M2**2*
     +    (S+TG)**(-1) + 2*(S+TG)**(-1)*S4**2 + S4 )
     +
      M2QGH = M2QGH + ANG2(21)*N*CO*TWO**(-1) * ( 16*S**(-1)*TG**(-1)*
     +    M2*MS2*S4 + 16*S**(-1)*TG**(-1)*M2**2*S4 + 4*S**(-1)*TG*M2*
     +    MS2*(S+TG)**(-1) - 4*S**(-1)*TG*M2 + 6*S**(-1)*TG*M2**2*
     +    (S+TG)**(-1) + 2*S**(-1)*TG*MS2*(S+TG)**(-1)*S4 + 2*S**(-1)*
     +    TG*MS2 + 2*S**(-1)*TG*(S+TG)**(-1)*S4**2 + 6*S**(-1)*TG*S4 - 
     +    2*S**(-1)*TG**2*MS2*(S+TG)**(-1) - 2*S**(-1)*TG**2*
     +    (S+TG)**(-1)*S4 - 2*S**(-1)*TG**2 - 8*S**(-1)*M2*MS2*
     +    (S+TG)**(-1)*S4 - 4*S**(-1)*M2*MS2 + 4*S**(-1)*M2*S4 - 8*
     +    S**(-1)*M2**2*(S+TG)**(-1)*S4 - 10*S**(-1)*M2**2 - 2*S**(-1)*
     +    MS2*S4 - 4*S**(-1)*S4**2 + 4*S*TG**(-1)*M2*MS2*(T1+UG)**(-1)
     +     + 2*S*TG**(-1)*M2 + 4*S*TG**(-1)*M2**2*(T1+UG)**(-1) + 2*S*
     +    TG**(-1)*MS2 + 2*S*TG*(S+TG)**(-1) + 7*S*TG*(T1+UG)**(-1) - 8
     +    *S*U1**(-1)*M2*MS2*(S+TG)**(-1) + 4*S*U1**(-1)*M2 + 2*S*
     +    U1**(-1)*S4 - 4*S*M2*(S+TG)**(-1) + 2*S*M2*(T1+UG)**(-1) - 2*
     +    S*MS2*(T1+UG)**(-1) + 3*S*(T1+UG)**(-1)*S4 - 8*S + 2*S**2*
     +    TG**(-1)*M2*(T1+UG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(21)*N*CO*TWO**(-1) * ( 2*S**2*TG**(-1)*MS2*
     +    (T1+UG)**(-1) + 8*TG**(-1)*U1**(-1)*M2*MS2*S4 + 8*TG**(-1)*
     +    U1**(-1)*M2**2*S4 + 4*TG**(-1)*M2*MS2*(T1+UG)**(-1)*S4 - 12*
     +    TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*(T1+UG)**(-1)*S4**2 + 2*
     +    TG**(-1)*M2*S4 + 4*TG**(-1)*M2**2*(T1+UG)**(-1)*S4 - 12*
     +    TG**(-1)*M2**2 - 2*TG**(-1)*MS2*(T1+UG)**(-1)*S4**2 + 2*
     +    TG**(-1)*MS2*S4 + 4*TG*M2*(T1+UG)**(-1) - 2*TG*MS2*
     +    (S+TG)**(-1) - 4*TG*(S+TG)**(-1)*S4 - 3*TG*(T1+UG)**(-1)*S4
     +     - 3*TG + 2*TG**2*(S+TG)**(-1) + 2*TG**2*(T1+UG)**(-1) + 4*
     +    U1**(-1)*M2*MS2*(S+TG)**(-1)*S4 - 4*U1**(-1)*M2*MS2 - 2*
     +    U1**(-1)*M2**2*(S+TG)**(-1)*S4 - 2*U1**(-1)*M2**2 + 8*M2*MS2*
     +    (S+TG)**(-1) + 4*M2*(S+TG)**(-1)*S4 + 2*M2*(T1+UG)**(-1)*S4
     +     - 6*M2 + 12*M2**2*(S+TG)**(-1) + 4*M2**2*(T1+UG)**(-1) + 2*
     +    MS2*(S+TG)**(-1)*S4 + 2*MS2*(T1+UG)**(-1)*S4 - 2*MS2 - 
     +    (T1+UG)**(-1)*S4**2 + 9*S4 )
     +
      M2QGH = M2QGH + ANG2(22)*N*CO*TWO**(-1) * ( 18*S**(-1)*TG*M2*MS2*
     +    (S+U1)**(-1) + 6*S**(-1)*TG*M2*(S+U1)**(-1)*S4 - 2*S**(-1)*TG
     +    *M2 + 6*S**(-1)*TG*M2**2*(S+U1)**(-1) - 18*S**(-1)*TG*MS2*
     +    (S+U1)**(-1)*S4 + 10*S**(-1)*TG*MS2 - 12*S**(-1)*TG*
     +    (S+U1)**(-1)*S4**2 + 8*S**(-1)*TG*S4 - 2*S**(-1)*TG**2*M2*
     +    (S+U1)**(-1) + 10*S**(-1)*TG**2*MS2*(S+U1)**(-1) + 12*S**(-1)
     +    *TG**2*(S+U1)**(-1)*S4 - 4*S**(-1)*TG**2 - 4*S**(-1)*TG**3*
     +    (S+U1)**(-1) - 2*S**(-1)*U1**(-1)*M2*MS2*S4 - 8*S**(-1)*
     +    U1**(-1)*M2**2*MS2 - 2*S**(-1)*U1**(-1)*M2**2*S4 - 4*S**(-1)*
     +    U1**(-1)*M2**3 - 16*S**(-1)*M2*MS2*(S+U1)**(-1)*S4 + 18*
     +    S**(-1)*M2*MS2 - 4*S**(-1)*M2*(S+U1)**(-1)*S4**2 + 4*S**(-1)*
     +    M2*S4 + 8*S**(-1)*M2**2*MS2*(S+U1)**(-1) - 4*S**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 + 6*S**(-1)*M2**2 + 4*S**(-1)*M2**3*
     +    (S+U1)**(-1) + 8*S**(-1)*MS2*(S+U1)**(-1)*S4**2 - 8*S**(-1)*
     +    MS2*S4 + 4*S**(-1)*(S+U1)**(-1)*S4**3 - 4*S**(-1)*S4**2 - 8*S
     +    *TG**(-1)*M2*MS2*(S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(22)*N*CO*TWO**(-1) * ( 8*S*TG**(-1)*M2*
     +    (S+U1)**(-1)*S4 - 2*S*TG**(-1)*M2 - 8*S*TG**(-1)*M2**2*
     +    (S+U1)**(-1) + 8*S*TG**(-1)*MS2*(S+U1)**(-1)*S4 - 8*S*
     +    TG**(-1)*MS2 - 4*S*TG**(-1)*S4 - 6*S*TG*(S+U1)**(-1) - 14*S*
     +    M2*(S+U1)**(-1) - 8*S*MS2*(S+U1)**(-1) + 6*S*(S+U1)**(-1)*S4
     +     - 2*S + 2*S**2*TG**(-1) - 2*S**2*(S+U1)**(-1) - 2*TG**(-1)*
     +    U1**(-1)*M2*MS2*S4 + 8*TG**(-1)*U1**(-1)*M2**2*MS2 - 2*
     +    TG**(-1)*U1**(-1)*M2**2*S4 + 4*TG**(-1)*U1**(-1)*M2**3 + 16*
     +    TG**(-1)*M2*MS2*(S+U1)**(-1)*S4 - 16*TG**(-1)*M2*MS2 - 8*
     +    TG**(-1)*M2*(S+U1)**(-1)*S4**2 + 4*TG**(-1)*M2*S4 - 8*
     +    TG**(-1)*M2**2*MS2*(S+U1)**(-1) + 16*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 8*TG**(-1)*M2**2 - 8*TG**(-1)*M2**3*
     +    (S+U1)**(-1) - 8*TG**(-1)*MS2*(S+U1)**(-1)*S4**2 + 10*
     +    TG**(-1)*MS2*S4 + 2*TG**(-1)*S4**2 - 14*TG*M2*(S+U1)**(-1) + 
     +    2*TG*MS2*(S+U1)**(-1) + 16*TG*(S+U1)**(-1)*S4 - 8*TG - 8*
     +    TG**2*(S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(22)*N*CO*TWO**(-1) * ( 8*U1**(-2)*M2**2*MS2
     +     - 8*U1**(-1)*M2*MS2 + 8*U1**(-1)*M2**2 - 8*M2*MS2*
     +    (S+U1)**(-1) + 24*M2*(S+U1)**(-1)*S4 - 14*M2 - 16*M2**2*
     +    (S+U1)**(-1) + 8*MS2*(S+U1)**(-1)*S4 + 2*MS2 - 8*(S+U1)**(-1)
     +    *S4**2 + 6*S4 )
     +
      M2QGH = M2QGH + ANG2(22)*N*CK*TWO**(-1) * (  - 6*S**(-1)*TG*M2 + 
     +    4*S**(-1)*TG*S4 - 2*S**(-1)*TG**2 + 2*S**(-1)*U1**(-1)*M2*MS2
     +    *S4 + 2*S**(-1)*U1**(-1)*M2**2*S4 + 4*S**(-1)*U1**(-1)*M2**3
     +     + 4*S**(-1)*M2*S4 - 8*S**(-1)*M2**2 - 2*S**(-1)*MS2*S4 - 2*
     +    S**(-1)*S4**2 + 6*S*TG*(T1+UG)**(-1) + 4*S*U1**(-1)*M2*MS2*
     +    (T1+UG)**(-1) + 8*S*M2*(T1+UG)**(-1) - 4*S*MS2*(T1+UG)**(-1)
     +     - 8*S*(T1+UG)**(-1)*S4 + 4*S + 4*S**2*(T1+UG)**(-1) + 6*TG*
     +    M2*(T1+UG)**(-1) - 6*TG*(T1+UG)**(-1)*S4 + 2*TG + 2*TG**2*
     +    (T1+UG)**(-1) - 8*U1**(-2)*M2**2*MS2 - 2*U1**(-1)*M2*MS2*
     +    (T1+UG)**(-1)*S4 + 12*U1**(-1)*M2*MS2 + 2*U1**(-1)*M2**2*
     +    (T1+UG)**(-1)*S4 - 8*U1**(-1)*M2**2 - 4*U1**(-1)*M2**3*
     +    (T1+UG)**(-1) - 10*M2*(T1+UG)**(-1)*S4 + 8*M2 + 8*M2**2*
     +    (T1+UG)**(-1) + 2*MS2*(T1+UG)**(-1)*S4 - 4*MS2 + 4*
     +    (T1+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2QGH = M2QGH + ANG2(23)*N*CO*TWO**(-1) * ( 32*S**(-2)*M2*MS2**2
     +     + 48*S**(-2)*M2**2*MS2 + 16*S**(-2)*M2**3 + 8*S**(-1)*
     +    TG**(-1)*M2*MS2*S4 + 16*S**(-1)*TG**(-1)*M2*MS2**2 + 24*
     +    S**(-1)*TG**(-1)*M2**2*MS2 + 8*S**(-1)*TG**(-1)*M2**2*S4 + 8*
     +    S**(-1)*TG**(-1)*M2**3 - 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) - 
     +    18*S**(-1)*TG*M2*MS2*(S+U1)**(-1) - 6*S**(-1)*TG*M2*
     +    (S+U1)**(-1)*S4 + 2*S**(-1)*TG*M2 - 6*S**(-1)*TG*M2**2*
     +    (S+TG)**(-1) - 6*S**(-1)*TG*M2**2*(S+U1)**(-1) - 2*S**(-1)*TG
     +    *MS2*(S+TG)**(-1)*S4 + 18*S**(-1)*TG*MS2*(S+U1)**(-1)*S4 - 4*
     +    S**(-1)*TG*MS2 - 2*S**(-1)*TG*(S+TG)**(-1)*S4**2 + 12*S**(-1)
     +    *TG*(S+U1)**(-1)*S4**2 - 10*S**(-1)*TG*S4 + 2*S**(-1)*TG**2*
     +    M2*(S+U1)**(-1) + 2*S**(-1)*TG**2*MS2*(S+TG)**(-1) - 10*
     +    S**(-1)*TG**2*MS2*(S+U1)**(-1) + 2*S**(-1)*TG**2*(S+TG)**(-1)
     +    *S4 - 12*S**(-1)*TG**2*(S+U1)**(-1)*S4 + 4*S**(-1)*TG**2 + 4*
     +    S**(-1)*TG**3*(S+U1)**(-1) - 2*S**(-1)*U1**(-1)*M2*MS2*S4 + 
     +    16*S**(-1)*U1**(-1)*M2*MS2**2 )
     +
      M2QGH = M2QGH + ANG2(23)*N*CO*TWO**(-1) * ( 16*S**(-1)*U1**(-1)*
     +    M2**2*MS2 - 2*S**(-1)*U1**(-1)*M2**2*S4 + 4*S**(-1)*U1**(-1)*
     +    M2**3 + 8*S**(-1)*M2*MS2*(S+TG)**(-1)*S4 + 16*S**(-1)*M2*MS2*
     +    (S+U1)**(-1)*S4 + 6*S**(-1)*M2*MS2 + 4*S**(-1)*M2*
     +    (S+U1)**(-1)*S4**2 - 8*S**(-1)*M2*S4 - 8*S**(-1)*M2**2*MS2*
     +    (S+U1)**(-1) + 8*S**(-1)*M2**2*(S+TG)**(-1)*S4 + 4*S**(-1)*
     +    M2**2*(S+U1)**(-1)*S4 + 12*S**(-1)*M2**2 - 4*S**(-1)*M2**3*
     +    (S+U1)**(-1) - 8*S**(-1)*MS2*(S+U1)**(-1)*S4**2 + 2*S**(-1)*
     +    MS2*S4 - 4*S**(-1)*(S+U1)**(-1)*S4**3 + 6*S**(-1)*S4**2 - 2*S
     +    *TG**(-1)*M2*MS2*(S+U1)**(-1) + 2*S*TG**(-1)*M2*(S+U1)**(-1)*
     +    S4 - 4*S*TG**(-1)*M2 - 2*S*TG**(-1)*M2**2*(S+U1)**(-1) + 4*S*
     +    TG**(-1)*MS2*(S+U1)**(-1)*S4 - 6*S*TG**(-1)*MS2 - 4*S*TG*
     +    (S+TG)**(-1) + 2*S*TG*(S+U1)**(-1) + 8*S*U1**(-1)*M2*MS2*
     +    (S+TG)**(-1) + 4*S*U1**(-1)*M2 - 4*S*U1**(-1)*MS2 - 2*S*
     +    U1**(-1)*S4 - 4*S*M2*(S+TG)**(-1) - 4*S*M2*(S+U1)**(-1) - 4*S
     +    *MS2*(S+TG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(23)*N*CO*TWO**(-1) * (  - 6*S*MS2*
     +    (S+U1)**(-1) + 4*S*(S+TG)**(-1)*S4 - 2*S*(S+U1)**(-1)*S4 + 2*
     +    S + 16*TG**(-1)*M2*MS2*(S+U1)**(-1)*S4 - 8*TG**(-1)*M2*MS2 - 
     +    4*TG**(-1)*M2*(S+U1)**(-1)*S4**2 + 6*TG**(-1)*M2*S4 - 8*
     +    TG**(-1)*M2**2*MS2*(S+U1)**(-1) + 8*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 8*TG**(-1)*M2**2 - 4*TG**(-1)*M2**3*
     +    (S+U1)**(-1) - 8*TG**(-1)*MS2*(S+U1)**(-1)*S4**2 + 4*TG**(-1)
     +    *MS2*S4 - 2*TG**(-1)*S4**2 - 8*TG*M2*(S+TG)**(-1) - 8*TG*M2*
     +    (S+U1)**(-1) - 2*TG*MS2*(S+TG)**(-1) - 22*TG*MS2*(S+U1)**(-1)
     +     + 8*TG*(S+TG)**(-1)*S4 - 12*TG*(S+U1)**(-1)*S4 + 10*TG - 4*
     +    TG**2*(S+TG)**(-1) + 6*TG**2*(S+U1)**(-1) + 8*U1**(-1)*M2*MS2
     +    *(S+TG)**(-1)*S4 - 4*U1**(-1)*M2*S4 + 8*U1**(-1)*M2**2*
     +    (S+TG)**(-1)*S4 - 2*U1**(-1)*MS2*S4 - 4*M2*MS2*(S+TG)**(-1)
     +     - 28*M2*MS2*(S+U1)**(-1) + 4*M2*(S+TG)**(-1)*S4 + 10*M2*
     +    (S+U1)**(-1)*S4 - 10*M2**2*(S+TG)**(-1) - 16*M2**2*
     +    (S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(23)*N*CO*TWO**(-1) * (  - 2*MS2*(S+TG)**(-1)
     +    *S4 + 28*MS2*(S+U1)**(-1)*S4 - 18*MS2 - 4*(S+TG)**(-1)*S4**2
     +     + 6*(S+U1)**(-1)*S4**2 - 12*S4 )
     +
      M2QGH = M2QGH + ANG2(24)*N*CO*TWO**(-1) * (  - 8*TG**(-1)*M2 )
     +
      M2QGH = M2QGH + ANG2(25)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG**(-1)*M2
     +    *MS2*S4**(-1) - 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*TG**(-1)*MS2 + 2*S**(-1)*M2*
     +    S4**(-1) + 2*S**(-1)*MS2*S4**(-1) + 2*S*U1**(-1)*S4**(-1) - 2
     +    *TG**(-1)*U1**(-1)*M2*MS2*S4**(-1) + 2*TG**(-1)*U1**(-1)*M2
     +     - 2*TG**(-1)*U1**(-1)*M2**2*S4**(-1) + 2*TG**(-1)*U1**(-1)*
     +    MS2 + 4*TG**(-1)*M2*(S+U1)**(-1) - 4*TG**(-1)*M2*S4**(-1) - 4
     +    *TG**(-1)*(S+U1)**(-1)*S4 + 4*TG**(-1) - 4*U1**(-1)*M2*
     +    S4**(-1) - 2*U1**(-1)*MS2*S4**(-1) + 4*(S+U1)**(-1) - 2*
     +    S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(26)*N*CK*TWO**(-1) * ( 2*S**(-1)*TG*S4**(-1)
     +     + 8*S**(-1)*U1**(-1)*M2*MS2*S4**(-1) + 4*S**(-1)*U1**(-1)*M2
     +     - 2*S**(-1)*U1**(-1)*S4 + 4*S**(-1)*M2*S4**(-1) - 2*S**(-1)
     +     - 2*S*TG**(-1)*S4**(-1) - 8*TG**(-1)*U1**(-1)*M2*MS2*
     +    S4**(-1) - 4*TG**(-1)*U1**(-1)*M2 + 2*TG**(-1)*U1**(-1)*S4 - 
     +    4*TG**(-1)*M2*S4**(-1) + 2*TG**(-1) )
     +
      M2QGH = M2QGH + ANG2(26)*CQED*TWO**(-1) * (  - 2*S**(-1)*TG*
     +    S4**(-1) - 8*S**(-1)*U1**(-1)*M2*MS2*S4**(-1) - 4*S**(-1)*
     +    U1**(-1)*M2 + 2*S**(-1)*U1**(-1)*S4 - 4*S**(-1)*M2*S4**(-1)
     +     + 2*S**(-1) - 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(27)*N*CO*TWO**(-1) * ( 16*TG**2*M2 )
     +
      M2QGH = M2QGH + ANG2(28)*N*CO*TWO**(-1) * ( 16*TG*M2 + 32*M2*MS2
     +     + 32*M2**2 )
     +
      M2QGH = M2QGH + ANG2(29)*N*CO*TWO**(-1) * (  - 16*S*M2 - 16*S**2*
     +    U1**(-1)*M2 - 16*TG*M2 + 16*U1**(-1)*M2*S4**2 - 16*M2*S4 )
     +
      M2QGH = M2QGH + ANG2(30)*N*CO*TWO**(-1) * (  - 16*S**(-1)*
     +    TG**(-1)*M2*MS2*S4 - 16*S**(-1)*TG**(-1)*M2**2*S4 - 4*S**(-1)
     +    *TG*M2*MS2*(S+TG)**(-1) - 6*S**(-1)*TG*M2**2*(S+TG)**(-1) - 2
     +    *S**(-1)*TG*MS2*(S+TG)**(-1)*S4 - 2*S**(-1)*TG*MS2 - 2*
     +    S**(-1)*TG*(S+TG)**(-1)*S4**2 - 2*S**(-1)*TG*S4 + 2*S**(-1)*
     +    TG**2*MS2*(S+TG)**(-1) + 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 + 8*
     +    S**(-1)*M2*MS2*(S+TG)**(-1)*S4 + 12*S**(-1)*M2*MS2 + 4*
     +    S**(-1)*M2*S4 + 8*S**(-1)*M2**2*(S+TG)**(-1)*S4 + 10*S**(-1)*
     +    M2**2 + 2*S**(-1)*MS2*S4 + 8*S*TG**(-1)*M2*MS2*(S+U1)**(-1)
     +     + 8*S*TG**(-1)*M2*MS2*S4**(-1) - 6*S*TG**(-1)*M2*
     +    (S+U1)**(-1)*S4 + 6*S*TG**(-1)*M2 + 8*S*TG**(-1)*M2**2*
     +    (S+U1)**(-1) + 8*S*TG**(-1)*M2**2*S4**(-1) - 4*S*TG**(-1)*MS2
     +    *(S+U1)**(-1)*S4 + 4*S*TG**(-1)*MS2 + 3*S*TG**(-1)*
     +    (S+U1)**(-1)*S4**2 - 3*S*TG**(-1)*S4 - 4*S*TG*(S+TG)**(-1) + 
     +    8*S*TG*(S+U1)**(-1) + 2*S*TG*S4**(-1) + 16*S*U1**(-2)*M2*S4
     +     + 8*S*U1**(-1)*M2*MS2*(S+TG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(30)*N*CO*TWO**(-1) * ( 8*S*U1**(-1)*M2*MS2*
     +    S4**(-1) - 24*S*U1**(-1)*M2 - 8*S*U1**(-1)*S4 - 4*S*M2*
     +    (S+TG)**(-1) + 10*S*M2*(S+U1)**(-1) + 4*S*M2*S4**(-1) - 4*S*
     +    MS2*(S+TG)**(-1) + 4*S*MS2*(S+U1)**(-1) + 4*S*(S+TG)**(-1)*S4
     +     - 7*S*(S+U1)**(-1)*S4 + 12*S - 2*S**2*TG**(-1)*(S+U1)**(-1)*
     +    S4 + 2*S**2*TG**(-1) - 16*S**2*U1**(-2)*M2 + 4*S**2*
     +    (S+U1)**(-1) + 2*S**2*S4**(-1) - 32*TG**(-1)*U1**(-1)*M2*MS2*
     +    S4 - 32*TG**(-1)*U1**(-1)*M2**2*S4 + 16*TG**(-1)*M2*MS2 + 4*
     +    TG**(-1)*M2*(S+U1)**(-1)*S4**2 - 4*TG**(-1)*M2*S4 - 2*
     +    TG**(-1)*M2**2*(S+U1)**(-1)*S4 + 18*TG**(-1)*M2**2 - 2*
     +    TG**(-1)*(S+U1)**(-1)*S4**3 + 2*TG**(-1)*S4**2 - 8*TG*M2*
     +    (S+TG)**(-1) + 10*TG*M2*(S+U1)**(-1) - 2*TG*MS2*(S+TG)**(-1)
     +     + 8*TG*(S+TG)**(-1)*S4 - 5*TG*(S+U1)**(-1)*S4 + 10*TG - 4*
     +    TG**2*(S+TG)**(-1) + 4*TG**2*(S+U1)**(-1) + 8*U1**(-1)*M2*MS2
     +    *(S+TG)**(-1)*S4 + 56*U1**(-1)*M2*MS2 + 24*U1**(-1)*M2*S4 + 8
     +    *U1**(-1)*M2**2*(S+TG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(30)*N*CO*TWO**(-1) * ( 8*U1**(-1)*M2**2 - 8*
     +    U1**(-1)*S4**2 - 4*M2*MS2*(S+TG)**(-1) + 8*M2*MS2*S4**(-1) + 
     +    4*M2*(S+TG)**(-1)*S4 - 10*M2*(S+U1)**(-1)*S4 + 18*M2 - 10*
     +    M2**2*(S+TG)**(-1) + 6*M2**2*(S+U1)**(-1) + 4*M2**2*S4**(-1)
     +     - 2*MS2*(S+TG)**(-1)*S4 + 4*MS2 - 4*(S+TG)**(-1)*S4**2 + 5*
     +    (S+U1)**(-1)*S4**2 - 7*S4 )
     +
      M2QGH = M2QGH + ANG2(31)*N*CO*TWO**(-1) * (  - 2*S*TG**(-1)*
     +    S4**(-1) - 8*TG**(-1)*U1**(-1)*M2*MS2*S4**(-1) - 4*TG**(-1)*
     +    U1**(-1)*M2 + 2*TG**(-1)*U1**(-1)*S4 - 4*TG**(-1)*M2*S4**(-1)
     +     + 2*TG**(-1) - 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(31)*N*CK*TWO**(-1) * (  - 2*S**(-1)*TG*
     +    S4**(-1) - 8*S**(-1)*U1**(-1)*M2*MS2*S4**(-1) - 4*S**(-1)*
     +    U1**(-1)*M2 + 2*S**(-1)*U1**(-1)*S4 - 4*S**(-1)*M2*S4**(-1)
     +     + 2*S**(-1) - 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(32)*N*CO*TWO**(-1) * ( 4*S*MS2 + 4*TG*MS2 - 
     +    8*U1**(-1)*M2**2*MS2 + 8*M2*MS2 - 4*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(32)*N*CK*TWO**(-1) * (  - 12*S*MS2 - 12*TG*
     +    MS2 + 24*U1**(-1)*M2**2*MS2 - 24*M2*MS2 + 12*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(32)*CQED*TWO**(-1) * ( 4*S*MS2 + 4*TG*MS2 - 
     +    8*U1**(-1)*M2**2*MS2 + 8*M2*MS2 - 4*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(33)*N*CO*TWO**(-1) * (  - 18*S**(-1)*TG*M2*
     +    MS2*(S+U1)**(-1) - 6*S**(-1)*TG*M2*(S+U1)**(-1)*S4 + 2*
     +    S**(-1)*TG*M2 - 6*S**(-1)*TG*M2**2*(S+U1)**(-1) + 18*S**(-1)*
     +    TG*MS2*(S+U1)**(-1)*S4 - 10*S**(-1)*TG*MS2 + 12*S**(-1)*TG*
     +    (S+U1)**(-1)*S4**2 - 8*S**(-1)*TG*S4 + 2*S**(-1)*TG**2*M2*
     +    (S+U1)**(-1) - 10*S**(-1)*TG**2*MS2*(S+U1)**(-1) - 12*S**(-1)
     +    *TG**2*(S+U1)**(-1)*S4 + 4*S**(-1)*TG**2 + 4*S**(-1)*TG**3*
     +    (S+U1)**(-1) + 2*S**(-1)*U1**(-1)*M2*MS2*S4 + 8*S**(-1)*
     +    U1**(-1)*M2**2*MS2 + 2*S**(-1)*U1**(-1)*M2**2*S4 + 4*S**(-1)*
     +    U1**(-1)*M2**3 + 16*S**(-1)*M2*MS2*(S+U1)**(-1)*S4 - 18*
     +    S**(-1)*M2*MS2 + 4*S**(-1)*M2*(S+U1)**(-1)*S4**2 - 4*S**(-1)*
     +    M2*S4 - 8*S**(-1)*M2**2*MS2*(S+U1)**(-1) + 4*S**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 6*S**(-1)*M2**2 - 4*S**(-1)*M2**3*
     +    (S+U1)**(-1) - 8*S**(-1)*MS2*(S+U1)**(-1)*S4**2 + 8*S**(-1)*
     +    MS2*S4 - 4*S**(-1)*(S+U1)**(-1)*S4**3 + 4*S**(-1)*S4**2 - 2*S
     +    *TG**(-1)*M2*MS2*(S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(33)*N*CO*TWO**(-1) * ( 2*S*TG**(-1)*M2*
     +    (S+U1)**(-1)*S4 - 2*S*TG**(-1)*M2**2*(S+U1)**(-1) + 4*S*
     +    TG**(-1)*MS2*(S+U1)**(-1)*S4 + 2*S*TG*(S+U1)**(-1) - 4*S*M2*
     +    (S+U1)**(-1) - 6*S*MS2*(S+U1)**(-1) - 2*S*(S+U1)**(-1)*S4 + 4
     +    *S + 2*TG**(-1)*U1**(-1)*M2*MS2*S4 - 8*TG**(-1)*U1**(-1)*
     +    M2**2*MS2 + 2*TG**(-1)*U1**(-1)*M2**2*S4 - 4*TG**(-1)*
     +    U1**(-1)*M2**3 + 16*TG**(-1)*M2*MS2*(S+U1)**(-1)*S4 - 2*
     +    TG**(-1)*M2*MS2 - 4*TG**(-1)*M2*(S+U1)**(-1)*S4**2 - 8*
     +    TG**(-1)*M2**2*MS2*(S+U1)**(-1) + 8*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 2*TG**(-1)*M2**2 - 4*TG**(-1)*M2**3*
     +    (S+U1)**(-1) - 8*TG**(-1)*MS2*(S+U1)**(-1)*S4**2 - 8*TG*M2*
     +    (S+U1)**(-1) - 22*TG*MS2*(S+U1)**(-1) - 12*TG*(S+U1)**(-1)*S4
     +     + 8*TG + 6*TG**2*(S+U1)**(-1) - 8*U1**(-2)*M2**2*MS2 + 8*
     +    U1**(-1)*M2*MS2 - 8*U1**(-1)*M2**2 - 28*M2*MS2*(S+U1)**(-1)
     +     + 10*M2*(S+U1)**(-1)*S4 + 4*M2 - 16*M2**2*(S+U1)**(-1) + 28*
     +    MS2*(S+U1)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(33)*N*CO*TWO**(-1) * (  - 12*MS2 + 6*
     +    (S+U1)**(-1)*S4**2 - 8*S4 )
     +
      M2QGH = M2QGH + ANG2(33)*N*CK*TWO**(-1) * ( 8*S**(-1)*TG*M2 + 2*
     +    S**(-1)*TG*MS2 - 4*S**(-1)*TG*S4 + 2*S**(-1)*TG**2 - 2*
     +    S**(-1)*U1**(-1)*M2*MS2*S4 - 2*S**(-1)*U1**(-1)*M2**2*S4 - 4*
     +    S**(-1)*U1**(-1)*M2**3 + 2*S**(-1)*M2*MS2 - 6*S**(-1)*M2*S4
     +     + 10*S**(-1)*M2**2 + 2*S**(-1)*S4**2 + 2*S*TG**(-1)*M2*MS2*
     +    S4**(-1) - 2*S*TG**(-1)*M2 + 2*S*TG**(-1)*M2**2*S4**(-1) + 4*
     +    S*TG**(-1)*MS2 + 2*S*TG*(S+TG)**(-1) + 2*S*TG*S4**(-1) - 8*S*
     +    U1**(-1)*M2*MS2*(S+TG)**(-1) - 8*S*U1**(-1)*M2*MS2*S4**(-1)
     +     + 4*S*U1**(-1)*M2*(S+TG)**(-1)*S4 - 4*S*U1**(-1)*M2 - 4*S*
     +    U1**(-1)*MS2*(S+TG)**(-1)*S4 + 4*S*U1**(-1)*MS2 + 2*S*
     +    U1**(-1)*(S+TG)**(-1)*S4**2 - 2*S*U1**(-1)*S4 - 2*S*M2*
     +    (S+TG)**(-1) + 4*S*M2*S4**(-1) + 4*S*MS2*(S+TG)**(-1) + 2*S*
     +    MS2*S4**(-1) - 3*S*(S+TG)**(-1)*S4 - S - S**2*(S+TG)**(-1) - 
     +    6*TG**(-1)*U1**(-1)*M2*MS2*S4 + 8*TG**(-1)*U1**(-1)*M2**2*MS2
     +     - 6*TG**(-1)*U1**(-1)*M2**2*S4 + 12*TG**(-1)*U1**(-1)*M2**3
     +     - 10*TG**(-1)*M2*MS2 )
     +
      M2QGH = M2QGH + ANG2(33)*N*CK*TWO**(-1) * ( 8*TG**(-1)*M2*S4 + 8*
     +    TG**(-1)*M2**2*MS2*S4**(-1) - 10*TG**(-1)*M2**2 + 4*TG**(-1)*
     +    M2**3*S4**(-1) + 8*TG**(-1)*MS2*S4 + 6*TG*M2*(S+TG)**(-1) + 8
     +    *TG*M2*S4**(-1) + 2*TG*MS2*S4**(-1) - 3*TG*(S+TG)**(-1)*S4 - 
     +    3*TG + 3*TG**2*(S+TG)**(-1) + 2*TG**2*S4**(-1) + 8*U1**(-2)*
     +    M2**2*MS2 - 8*U1**(-1)*M2*MS2 + 8*U1**(-1)*M2**2*MS2*S4**(-1)
     +     + U1**(-1)*M2**2*(S+TG)**(-1)*S4 + 15*U1**(-1)*M2**2 + 2*M2*
     +    MS2*S4**(-1) - 4*M2*(S+TG)**(-1)*S4 - 12*M2 + 3*M2**2*
     +    (S+TG)**(-1) + 10*M2**2*S4**(-1) - 4*MS2 + 2*(S+TG)**(-1)*
     +    S4**2 - S4 )
     +
      M2QGH = M2QGH + ANG2(33)*CQED*TWO**(-1) * ( 2*S**(-1)*TG*M2*MS2*
     +    S4**(-1) - 14*S**(-1)*TG*M2 + 10*S**(-1)*TG*M2**2*S4**(-1) - 
     +    2*S**(-1)*TG*MS2 + 6*S**(-1)*TG*S4 + 8*S**(-1)*TG**2*M2*
     +    S4**(-1) + 2*S**(-1)*TG**2*MS2*S4**(-1) - 6*S**(-1)*TG**2 + 2
     +    *S**(-1)*TG**3*S4**(-1) + 6*S**(-1)*M2*S4 - 8*S**(-1)*M2**2
     +     + 4*S**(-1)*M2**3*S4**(-1) - 2*S**(-1)*S4**2 - 4*S*TG**(-1)*
     +    MS2 + 8*S*U1**(-1)*M2*MS2*S4**(-1) + 2*TG**(-1)*U1**(-1)*M2*
     +    MS2*S4 + 2*TG**(-1)*U1**(-1)*M2**2*S4 - 4*TG**(-1)*U1**(-1)*
     +    M2**3 - 2*TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*S4 + 2*TG**(-1)*
     +    M2**2 + 4*TG*M2*S4**(-1) + 2*TG*MS2*S4**(-1) - 4*TG + 2*TG**2
     +    *S4**(-1) - 8*U1**(-1)*M2**2*MS2*S4**(-1) - 4*U1**(-1)*M2**2
     +     + 10*M2*MS2*S4**(-1) - 2*M2 + 2*M2**2*S4**(-1) - 2*MS2 + 2*
     +    S4 )
     +
      M2QGH = M2QGH + ANG2(34)*N*CO*TWO**(-1) * ( 4*TG**(-1)*M2*
     +    (S+U1)**(-1) - 4*TG**(-1)*(S+U1)**(-1)*S4 + 4*TG**(-1) + 4*
     +    (S+U1)**(-1) )
     +
      M2QGH = M2QGH + ANG2(34)*N*CK*TWO**(-1) * ( 2*S**(-1)*TG**(-1)*M2
     +    *MS2*S4**(-1) - 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*TG**(-1)*MS2 + 2*S**(-1)*M2*
     +    S4**(-1) + 2*S**(-1)*MS2*S4**(-1) + 2*S*U1**(-1)*S4**(-1) - 2
     +    *TG**(-1)*U1**(-1)*M2*MS2*S4**(-1) + 2*TG**(-1)*U1**(-1)*M2
     +     - 2*TG**(-1)*U1**(-1)*M2**2*S4**(-1) + 2*TG**(-1)*U1**(-1)*
     +    MS2 - 4*TG**(-1)*M2*S4**(-1) - 4*U1**(-1)*M2*S4**(-1) - 2*
     +    U1**(-1)*MS2*S4**(-1) - 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(35)*N*CO*TWO**(-1) * ( 8*S**2*M2 )
     +
      M2QGH = M2QGH + ANG2(35)*N*CK*TWO**(-1) * (  - 8*S**2*M2 )
     +
      M2QGH = M2QGH + ANG2(36)*N*CO*TWO**(-1) * ( 8*S*M2 )
     +
      M2QGH = M2QGH + ANG2(36)*N*CK*TWO**(-1) * (  - 8*S*M2 )
     +
      M2QGH = M2QGH + ANG2(37)*N*CO*TWO**(-1) * ( 16*S*U1**(-1)*M2*S4
     +     - 16*S*M2 - 16*S**2*U1**(-1)*M2 )
     +
      M2QGH = M2QGH + ANG2(37)*N*CK*TWO**(-1) * (  - 16*S*U1**(-1)*M2*
     +    S4 + 16*S*M2 )
     +
      M2QGH = M2QGH + ANG2(38)*N*CO*TWO**(-1) * ( 4*S*TG**(-1)*M2*MS2*
     +    (T1+UG)**(-1) + 2*S*TG**(-1)*M2 + 4*S*TG**(-1)*M2**2*
     +    (T1+UG)**(-1) + 2*S*TG**(-1)*MS2 + 7*S*TG*(T1+UG)**(-1) + 16*
     +    S*U1**(-2)*M2*S4 - 16*S*U1**(-1)*M2 + 2*S*M2*(T1+UG)**(-1) - 
     +    2*S*MS2*(T1+UG)**(-1) + 3*S*(T1+UG)**(-1)*S4 + 2*S + 2*S**2*
     +    TG**(-1)*M2*(T1+UG)**(-1) + 2*S**2*TG**(-1)*MS2*(T1+UG)**(-1)
     +     - 16*S**2*U1**(-2)*M2 - 16*TG**(-1)*U1**(-1)*M2*MS2*S4 - 16*
     +    TG**(-1)*U1**(-1)*M2**2*S4 + 4*TG**(-1)*M2*MS2*(T1+UG)**(-1)*
     +    S4 + 12*TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*(T1+UG)**(-1)*S4**2
     +     + 6*TG**(-1)*M2*S4 + 4*TG**(-1)*M2**2*(T1+UG)**(-1)*S4 + 8*
     +    TG**(-1)*M2**2 - 2*TG**(-1)*MS2*(T1+UG)**(-1)*S4**2 + 2*
     +    TG**(-1)*MS2*S4 - 2*TG**(-1)*S4**2 + 4*TG*M2*(T1+UG)**(-1) - 
     +    3*TG*(T1+UG)**(-1)*S4 + 7*TG + 2*TG**2*(T1+UG)**(-1) + 16*
     +    U1**(-1)*M2*MS2 + 8*U1**(-1)*M2*S4 - 16*U1**(-1)*M2**2 - 8*
     +    U1**(-1)*S4**2 + 2*M2*(T1+UG)**(-1)*S4 + 6*M2 + 4*M2**2*
     +    (T1+UG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(38)*N*CO*TWO**(-1) * ( 2*MS2*(T1+UG)**(-1)*
     +    S4 - 2*MS2 - (T1+UG)**(-1)*S4**2 + 5*S4 )
     +
      M2QGH = M2QGH + ANG2(38)*N*CK*TWO**(-1) * (  - 6*S*TG*
     +    (S+TG)**(-1) - 2*S*TG*S4**(-1) + 8*S*U1**(-1)*M2*MS2*
     +    (S+TG)**(-1) + 8*S*U1**(-1)*M2*MS2*S4**(-1) - 8*S*U1**(-1)*M2
     +     - 8*S*U1**(-1)*S4 - 4*S*M2*(S+TG)**(-1) + 2*S - 2*S**2*
     +    (S+TG)**(-1) - 8*TG*M2*(S+TG)**(-1) - 4*TG*M2*S4**(-1) + 2*TG
     +    *(S+TG)**(-1)*S4 - 4*TG**2*(S+TG)**(-1) - 2*TG**2*S4**(-1) - 
     +    8*U1**(-1)*M2*S4 + 8*U1**(-1)*M2**2*(S+TG)**(-1)*S4 + 24*
     +    U1**(-1)*M2**2 + 8*U1**(-1)*S4**2 - 4*M2*(S+TG)**(-1)*S4 - 4*
     +    M2**2*(S+TG)**(-1) - 4*M2**2*S4**(-1) + 2*(S+TG)**(-1)*S4**2
     +     - 6*S4 )
     +
      M2QGH = M2QGH + ANG2(39)*N*CO*TWO**(-1) * (  - 2*S*M2 )
     +
      M2QGH = M2QGH + ANG2(39)*N*CK*TWO**(-1) * ( 6*S*M2 )
     +
      M2QGH = M2QGH + ANG2(39)*CQED*TWO**(-1) * (  - 2*S*M2 )
     +
      M2QGH = M2QGH + ANG2(40)*N*CO*TWO**(-1) * ( 2*S*TG**(-1)*M2*MS2*
     +    (S+U1)**(-1) - 2*S*TG**(-1)*M2*(S+U1)**(-1)*S4 - 2*S*TG**(-1)
     +    *M2 + 2*S*TG**(-1)*M2**2*(S+U1)**(-1) - 4*S*TG**(-1)*MS2*
     +    (S+U1)**(-1)*S4 - 4*S*TG*(S+U1)**(-1) + 4*S*M2*(S+U1)**(-1)
     +     + 6*S*MS2*(S+U1)**(-1) + 4*S*(S+U1)**(-1)*S4 - 4*S - 4*
     +    TG**(-1)*M2*(S+U1)**(-1)*S4**2 + 2*TG**(-1)*M2*S4 + 8*
     +    TG**(-1)*M2**2*(S+U1)**(-1)*S4 - 4*TG**(-1)*M2**2 - 4*
     +    TG**(-1)*M2**3*(S+U1)**(-1) - 19*TG*M2*(S+U1)**(-1) + 7*TG*
     +    (S+U1)**(-1)*S4 - 5*TG - 7*TG**2*(S+U1)**(-1) + 17*M2*
     +    (S+U1)**(-1)*S4 - 11*M2 - 16*M2**2*(S+U1)**(-1) - 2*
     +    (S+U1)**(-1)*S4**2 + 2*S4 )
     +
      M2QGH = M2QGH + ANG2(40)*N*CK*TWO**(-1) * (  - 2*S*TG**(-1)*M2*
     +    MS2*S4**(-1) + 4*S*TG**(-1)*M2 - 2*S*TG**(-1)*M2**2*S4**(-1)
     +     + 12*S*TG**(-1)*MS2 + 6*S*TG*(T1+UG)**(-1) + 4*S*U1**(-1)*M2
     +    *MS2*(T1+UG)**(-1) - 8*S*U1**(-1)*M2 - 2*S*U1**(-1)*MS2 - 4*S
     +    *U1**(-1)*S4 + 8*S*M2*(T1+UG)**(-1) - 2*S*M2*S4**(-1) - 4*S*
     +    MS2*(T1+UG)**(-1) - 2*S*MS2*S4**(-1) - 8*S*(T1+UG)**(-1)*S4
     +     + 6*S + 2*S**2*U1**(-1) + 4*S**2*(T1+UG)**(-1) - 4*TG**(-1)*
     +    U1**(-1)*M2*MS2*S4 - 4*TG**(-1)*U1**(-1)*M2**2*S4 + 8*
     +    TG**(-1)*U1**(-1)*M2**3 + 4*TG**(-1)*M2*MS2 + 6*TG**(-1)*M2*
     +    S4 - 8*TG**(-1)*M2**2 + 4*TG**(-1)*M2**3*S4**(-1) + 6*TG*M2*
     +    (T1+UG)**(-1) + 6*TG*M2*S4**(-1) - 6*TG*(T1+UG)**(-1)*S4 + 4*
     +    TG + 2*TG**2*(T1+UG)**(-1) + 2*TG**2*S4**(-1) - 2*U1**(-1)*M2
     +    *MS2*(T1+UG)**(-1)*S4 + 6*U1**(-1)*M2*MS2 + 6*U1**(-1)*M2*S4
     +     + 2*U1**(-1)*M2**2*(T1+UG)**(-1)*S4 + 10*U1**(-1)*M2**2 - 4*
     +    U1**(-1)*M2**3*(T1+UG)**(-1) + 2*U1**(-1)*S4**2 - 10*M2*
     +    (T1+UG)**(-1)*S4 )
     +
      M2QGH = M2QGH + ANG2(40)*N*CK*TWO**(-1) * (  - 2*M2 + 8*M2**2*
     +    (T1+UG)**(-1) + 8*M2**2*S4**(-1) + 2*MS2*(T1+UG)**(-1)*S4 - 2
     +    *MS2 + 4*(T1+UG)**(-1)*S4**2 - 6*S4 )
     +
      M2QGH = M2QGH + ANG2(40)*CQED*TWO**(-1) * (  - 4*S*TG**(-1)*MS2
     +     + 2*S*TG*S4**(-1) - 2*S*U1**(-1)*M2*MS2*S4**(-1) - 6*S*
     +    U1**(-1)*M2 - 10*S*U1**(-1)*M2**2*S4**(-1) - 2*S*U1**(-1)*S4
     +     + 8*S*M2*S4**(-1) + 2*S*MS2*S4**(-1) + 2*S + 8*S**2*U1**(-1)
     +    *M2*S4**(-1) + 2*S**2*U1**(-1)*MS2*S4**(-1) + 4*S**2*U1**(-1)
     +     - 2*S**2*S4**(-1) - 2*S**3*U1**(-1)*S4**(-1) + 2*TG**(-1)*
     +    U1**(-1)*M2*MS2*S4 + 2*TG**(-1)*U1**(-1)*M2**2*S4 - 4*
     +    TG**(-1)*U1**(-1)*M2**3 - 2*TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*
     +    S4 + 2*TG**(-1)*M2**2 - 6*TG*M2*S4**(-1) - 2*TG**2*S4**(-1)
     +     - 2*U1**(-1)*M2*MS2 - 2*U1**(-1)*M2**2 + 4*U1**(-1)*M2**3*
     +    S4**(-1) + 2*M2 - 8*M2**2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(41)*N*CO*TWO**(-1) * (  - 6*S**(-1)*M2*
     +    (T1+UG)**(-1) - 4*S**(-1)*MS2*(T1+UG)**(-1) + 2*S**(-1)*
     +    (T1+UG)**(-1)*S4 - 2*S**(-1) - 2*(T1+UG)**(-1) )
     +
      M2QGH = M2QGH + ANG2(41)*N*CK*TWO**(-1) * ( 2*S**(-1)*TG**(-1)*M2
     +    *MS2*S4**(-1) - 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*TG**(-1)*MS2 - 2*S**(-1)*U1**(-1)*
     +    M2*MS2*S4**(-1) + 2*S**(-1)*U1**(-1)*M2 - 2*S**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) + 2*S**(-1)*U1**(-1)*MS2 + 6*S**(-1)*M2*
     +    (T1+UG)**(-1) + 4*S**(-1)*M2*S4**(-1) + 4*S**(-1)*MS2*
     +    (T1+UG)**(-1) + 8*S**(-1)*MS2*S4**(-1) - 2*S**(-1)*
     +    (T1+UG)**(-1)*S4 + 2*S**(-1) - 2*S*U1**(-1)*S4**(-1) - 2*
     +    TG**(-1)*M2*S4**(-1) - 2*TG**(-1)*MS2*S4**(-1) + 4*U1**(-1)*
     +    M2*S4**(-1) + 6*U1**(-1)*MS2*S4**(-1) + 2*U1**(-1) + 2*
     +    (T1+UG)**(-1) - 4*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(41)*CQED*TWO**(-1) * ( 2*S**(-1)*U1**(-1)*M2
     +    *MS2*S4**(-1) - 2*S**(-1)*U1**(-1)*M2 + 2*S**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) - 2*S**(-1)*U1**(-1)*MS2 - 2*S**(-1)*M2*
     +    S4**(-1) - 4*S**(-1)*MS2*S4**(-1) + 2*S*U1**(-1)*S4**(-1) - 4
     +    *U1**(-1)*M2*S4**(-1) - 6*U1**(-1)*MS2*S4**(-1) - 2*U1**(-1)
     +     + 2*S4**(-1) )
     +
      M2QGH = M2QGH + ANG2(42)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG*M2*MS2*
     +    (T1+UG)**(-1) - 2*S**(-1)*TG*M2 + 2*S**(-1)*TG*M2**2*
     +    (T1+UG)**(-1) - 2*S**(-1)*TG*MS2 - 4*S**(-1)*M2*MS2*
     +    (T1+UG)**(-1)*S4 + 2*S**(-1)*M2*MS2 - 4*S**(-1)*M2*S4 + 6*
     +    S**(-1)*M2**2 - 4*S**(-1)*M2**3*(T1+UG)**(-1) + 2*S**(-1)*MS2
     +    *S4 + 2*S**(-1)*S4**2 + 2*S*TG**(-1)*M2*MS2*(S+U1)**(-1) - 4*
     +    S*TG**(-1)*M2*MS2*(T1+UG)**(-1) - 2*S*TG**(-1)*M2*
     +    (S+U1)**(-1)*S4 - 2*S*TG**(-1)*M2 + 2*S*TG**(-1)*M2**2*
     +    (S+U1)**(-1) - 4*S*TG**(-1)*M2**2*(T1+UG)**(-1) - 4*S*
     +    TG**(-1)*MS2*(S+U1)**(-1)*S4 - 4*S*TG*(S+U1)**(-1) - 16*S*
     +    U1**(-2)*M2*S4 + 16*S*U1**(-1)*M2 + 4*S*M2*(S+U1)**(-1) - 5*S
     +    *M2*(T1+UG)**(-1) + 6*S*MS2*(S+U1)**(-1) + 6*S*MS2*
     +    (T1+UG)**(-1) + 4*S*(S+U1)**(-1)*S4 + S*(T1+UG)**(-1)*S4 - 7*
     +    S - 2*S**2*TG**(-1)*M2*(T1+UG)**(-1) - 2*S**2*TG**(-1)*MS2*
     +    (T1+UG)**(-1) + 16*S**2*U1**(-2)*M2 - 3*S**2*(T1+UG)**(-1) + 
     +    16*TG**(-1)*U1**(-1)*M2*MS2*S4 )
     +
      M2QGH = M2QGH + ANG2(42)*N*CO*TWO**(-1) * ( 16*TG**(-1)*U1**(-1)*
     +    M2**2*S4 - 4*TG**(-1)*M2*MS2*(T1+UG)**(-1)*S4 - 14*TG**(-1)*
     +    M2*MS2 - 4*TG**(-1)*M2*(S+U1)**(-1)*S4**2 + 2*TG**(-1)*M2*
     +    (T1+UG)**(-1)*S4**2 - 2*TG**(-1)*M2*S4 + 8*TG**(-1)*M2**2*
     +    (S+U1)**(-1)*S4 - 4*TG**(-1)*M2**2*(T1+UG)**(-1)*S4 - 14*
     +    TG**(-1)*M2**2 - 4*TG**(-1)*M2**3*(S+U1)**(-1) + 2*TG**(-1)*
     +    MS2*(T1+UG)**(-1)*S4**2 + 2*TG**(-1)*S4**2 - 19*TG*M2*
     +    (S+U1)**(-1) + 2*TG*M2*(T1+UG)**(-1) + 2*TG*MS2*(T1+UG)**(-1)
     +     + 7*TG*(S+U1)**(-1)*S4 - 7*TG - 7*TG**2*(S+U1)**(-1) - 16*
     +    U1**(-1)*M2*MS2 - 8*U1**(-1)*M2*S4 + 16*U1**(-1)*M2**2 + 8*
     +    U1**(-1)*S4**2 + 12*M2*MS2*(T1+UG)**(-1) + 17*M2*(S+U1)**(-1)
     +    *S4 - M2*(T1+UG)**(-1)*S4 - 16*M2 - 16*M2**2*(S+U1)**(-1) - 6
     +    *MS2*(T1+UG)**(-1)*S4 + 6*MS2 - 2*(S+U1)**(-1)*S4**2 - 2*S4 )
     +
      M2QGH = M2QGH + ANG2(43)*N*CO*TWO**(-1) * ( 2*S**(-1)*M2 )
     +
      M2QGH = M2QGH + ANG2(43)*N*CK*TWO**(-1) * (  - 6*S**(-1)*M2 )
     +
      M2QGH = M2QGH + ANG2(43)*CQED*TWO**(-1) * ( 2*S**(-1)*M2 )
     +
      M2QGH = M2QGH + ANG2(44)*N*CK*TWO**(-1) * ( 2*S**(-1)*TG*M2*MS2*
     +    (T1+UG)**(-1) - 2*S**(-1)*TG*M2 + 2*S**(-1)*TG*M2**2*
     +    (T1+UG)**(-1) - 2*S**(-1)*TG*MS2 - 4*S**(-1)*M2*MS2*
     +    (T1+UG)**(-1)*S4 + 4*S**(-1)*M2*MS2 - 6*S**(-1)*M2*S4 + 8*
     +    S**(-1)*M2**2 - 4*S**(-1)*M2**3*(T1+UG)**(-1) + 2*S**(-1)*
     +    S4**2 - 8*S*TG**(-1)*MS2 - 6*S*TG*(S+TG)**(-1) + 8*S*U1**(-1)
     +    *M2*MS2*(S+TG)**(-1) - 4*S*U1**(-1)*M2*MS2*(T1+UG)**(-1) + 4*
     +    S*U1**(-1)*M2 - 4*S*U1**(-1)*MS2 - 2*S*U1**(-1)*S4 - 4*S*M2*
     +    (S+TG)**(-1) - 4*S*M2*(T1+UG)**(-1) + 2*S - 2*S**2*
     +    (S+TG)**(-1) + 4*TG**(-1)*U1**(-1)*M2*MS2*S4 + 4*TG**(-1)*
     +    U1**(-1)*M2**2*S4 - 8*TG**(-1)*U1**(-1)*M2**3 - 4*TG**(-1)*M2
     +    *MS2 - 4*TG**(-1)*M2*S4 + 4*TG**(-1)*M2**2 - 8*TG*M2*
     +    (S+TG)**(-1) - 2*TG*M2*(T1+UG)**(-1) - 2*TG*MS2*(T1+UG)**(-1)
     +     + 2*TG*(S+TG)**(-1)*S4 + 4*TG - 4*TG**2*(S+TG)**(-1) + 2*
     +    U1**(-1)*M2*MS2*(T1+UG)**(-1)*S4 - 4*U1**(-1)*M2*MS2 - 8*
     +    U1**(-1)*M2*S4 )
     +
      M2QGH = M2QGH + ANG2(44)*N*CK*TWO**(-1) * ( 8*U1**(-1)*M2**2*
     +    (S+TG)**(-1)*S4 - 2*U1**(-1)*M2**2*(T1+UG)**(-1)*S4 + 4*
     +    U1**(-1)*M2**3*(T1+UG)**(-1) - 2*U1**(-1)*MS2*S4 + 2*U1**(-1)
     +    *S4**2 + 16*M2*MS2*(T1+UG)**(-1) - 4*M2*(S+TG)**(-1)*S4 + 2*
     +    M2*(T1+UG)**(-1)*S4 + 4*M2 - 4*M2**2*(S+TG)**(-1) + 8*M2**2*
     +    (T1+UG)**(-1) - 2*MS2*(T1+UG)**(-1)*S4 - 4*MS2 + 2*
     +    (S+TG)**(-1)*S4**2 - 2*S4 )
     +
      M2QGH = M2QGH + ANG2(44)*CQED*TWO**(-1) * (  - 2*S**(-1)*U1**(-1)
     +    *M2*MS2*S4 + 8*S**(-1)*U1**(-1)*M2*S4**2 - 10*S**(-1)*
     +    U1**(-1)*M2**2*S4 + 4*S**(-1)*U1**(-1)*M2**3 + 2*S**(-1)*
     +    U1**(-1)*MS2*S4**2 - 2*S**(-1)*U1**(-1)*S4**3 + 2*S**(-1)*M2*
     +    MS2 - 2*S**(-1)*M2*S4 + 2*S**(-1)*M2**2 - 2*S**(-1)*MS2*S4 + 
     +    4*S*TG**(-1)*MS2 - 2*TG**(-1)*U1**(-1)*M2*MS2*S4 - 2*TG**(-1)
     +    *U1**(-1)*M2**2*S4 + 4*TG**(-1)*U1**(-1)*M2**3 + 2*TG**(-1)*
     +    M2*MS2 + 2*TG**(-1)*M2*S4 - 2*TG**(-1)*M2**2 - 8*U1**(-1)*M2*
     +    MS2 - 4*U1**(-1)*M2*S4 + 4*U1**(-1)*MS2*S4 + 2*U1**(-1)*S4**2
     +     )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(TG*UG-M2*S)**(-2)*(S4+MS2)*
     + TWO**(-1) * (  - 16*S*TG*T1**(-2)*M2**2*MS2*S4**(-1) + 16*S*TG*
     +    T1**(-2)*M2**3*MS2*S4**(-2) + 8*S*TG*M2*MS2*S4**(-2) + 8*S*
     +    TG**2*T1**(-2)*M2**2*MS2*S4**(-2) + 8*S*T1**(-2)*M2**2*MS2 - 
     +    16*S*T1**(-2)*M2**3*MS2*S4**(-1) + 8*S*T1**(-2)*M2**4*MS2*
     +    S4**(-2) + 8*S*M2**2*MS2*S4**(-2) - 8*S**2*TG*T1**(-2)*M2**2*
     +    MS2*S4**(-2) + 8*S**2*T1**(-2)*M2**2*MS2*S4**(-1) - 8*S**2*
     +    T1**(-2)*M2**3*MS2*S4**(-2) + 8*S**3*T1**(-1)*M2*MS2*S4**(-2)
     +     - 24*TG*T1**(-2)*M2**2*MS2 + 48*TG*T1**(-2)*M2**3*MS2*
     +    S4**(-1) - 24*TG*T1**(-2)*M2**4*MS2*S4**(-2) - 8*TG*M2**2*MS2
     +    *S4**(-2) + 24*TG**2*T1**(-2)*M2**2*MS2*S4**(-1) - 24*TG**2*
     +    T1**(-2)*M2**3*MS2*S4**(-2) - 8*TG**3*T1**(-2)*M2**2*MS2*
     +    S4**(-2) + 8*T1**(-2)*M2**2*MS2*S4 - 24*T1**(-2)*M2**3*MS2 + 
     +    24*T1**(-2)*M2**4*MS2*S4**(-1) - 8*T1**(-2)*M2**5*MS2*
     +    S4**(-2) + 8*M2**2*MS2*S4**(-1) - 8*M2**3*MS2*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * ( 24*S**(-1)*TG*T1**(-2)*M2**2 - 48*S**(-1)*TG*
     +    T1**(-2)*M2**3*S4**(-1) + 24*S**(-1)*TG*T1**(-2)*M2**4*
     +    S4**(-2) + 8*S**(-1)*TG*M2**2*S4**(-2) - 24*S**(-1)*TG**2*
     +    T1**(-2)*M2**2*S4**(-1) + 24*S**(-1)*TG**2*T1**(-2)*M2**3*
     +    S4**(-2) + 8*S**(-1)*TG**3*T1**(-2)*M2**2*S4**(-2) - 8*
     +    S**(-1)*T1**(-2)*M2**2*S4 + 24*S**(-1)*T1**(-2)*M2**3 - 24*
     +    S**(-1)*T1**(-2)*M2**4*S4**(-1) + 8*S**(-1)*T1**(-2)*M2**5*
     +    S4**(-2) - 8*S**(-1)*M2**2*S4**(-1) + 8*S**(-1)*M2**3*
     +    S4**(-2) - 8*S*TG*T1**(-2)*M2*MS2*S4**(-2) + 8*S*T1**(-2)*M2*
     +    MS2*S4**(-1) - 8*S*T1**(-2)*M2**2*MS2*S4**(-2) + 8*S**2*
     +    T1**(-2)*M2*MS2*S4**(-2) - 16*TG**(-1)*T1**(-1)*M2*MS2 + 32*
     +    TG**(-1)*T1**(-1)*M2**2*MS2*S4**(-1) - 16*TG**(-1)*T1**(-1)*
     +    M2**3*MS2*S4**(-2) - 16*TG**(-1)*M2**2*MS2*S4**(-2) - 16*TG*
     +    T1**(-2)*M2*MS2*S4**(-1) + 16*TG*T1**(-2)*M2**2*MS2*S4**(-2)
     +     - 16*TG*T1**(-1)*M2*MS2*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * ( 16*TG*T1**(-1)*M2*S4**(-1) - 16*TG*T1**(-1)*M2**2*
     +    S4**(-2) - 8*TG*M2*S4**(-2) + 8*TG**2*T1**(-2)*M2*MS2*
     +    S4**(-2) - 8*TG**2*T1**(-1)*M2*S4**(-2) + 8*T1**(-2)*M2*MS2
     +     - 16*T1**(-2)*M2**2*MS2*S4**(-1) + 8*T1**(-2)*M2**3*MS2*
     +    S4**(-2) + 32*T1**(-1)*M2*MS2*S4**(-1) - 8*T1**(-1)*M2 - 32*
     +    T1**(-1)*M2**2*MS2*S4**(-2) + 16*T1**(-1)*M2**2*S4**(-1) - 8*
     +    T1**(-1)*M2**3*S4**(-2) - 8*M2*MS2*S4**(-2) - 8*M2**2*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * ( 8*S**(-1)*TG
     +    *T1**(-2)*M2*S4**(-1) - 8*S**(-1)*TG*T1**(-2)*M2**2*S4**(-2)
     +     - 8*S**(-1)*TG*T1**(-1)*M2*S4**(-2) + 8*S**(-1)*TG*T1**(-1)*
     +    S4**(-1) - 16*S**(-1)*TG*(S+TG)**(-2) - 16*S**(-1)*TG*
     +    (S+TG)**(-1)*S4**(-1) - 20*S**(-1)*TG*S4**(-2) - 4*S**(-1)*
     +    TG**2*T1**(-2)*M2*S4**(-2) - 4*S**(-1)*TG**2*T1**(-1)*
     +    S4**(-2) + 16*S**(-1)*TG**2*(S+TG)**(-2)*S4**(-1) - 4*S**(-1)
     +    *T1**(-2)*M2 + 8*S**(-1)*T1**(-2)*M2**2*S4**(-1) - 4*S**(-1)*
     +    T1**(-2)*M2**3*S4**(-2) + 8*S**(-1)*T1**(-1)*M2*S4**(-1) - 4*
     +    S**(-1)*T1**(-1)*M2**2*S4**(-2) - 4*S**(-1)*T1**(-1) - 32*
     +    S**(-1)*U1**(-4)*M2**2*S4 + 32*S**(-1)*U1**(-3)*M2*S4 + 32*
     +    S**(-1)*U1**(-3)*M2**2 + 32*S**(-1)*U1**(-2)*M2*MS2*
     +    (S+TG)**(-1) - 32*S**(-1)*U1**(-2)*M2*MS2*S4**(-1) - 32*
     +    S**(-1)*U1**(-2)*M2 - 64*S**(-1)*U1**(-2)*M2**2*S4**(-1) - 16
     +    *S**(-1)*U1**(-2)*S4 + 32*S**(-1)*U1**(-1)*M2*MS2*
     +    (S+TG)**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 64*
     +    S**(-1)*U1**(-1)*M2*MS2*(S+TG)**(-1)*S4**(-1) + 64*S**(-1)*
     +    U1**(-1)*M2*S4**(-1) - 32*S**(-1)*U1**(-1)*M2**2*(S+TG)**(-1)
     +    *S4**(-1) + 32*S**(-1)*U1**(-1)*M2**2*S4**(-2) + 16*S**(-1)*
     +    U1**(-1) - 32*S**(-1)*M2*MS2*(S+TG)**(-2)*S4**(-1) + 32*
     +    S**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 40*S**(-1)*M2*S4**(-2) - 8
     +    *S*TG**(-1)*M2**2*(S+UG)**(-2)*S4**(-2) - 32*S*TG**(-1)*
     +    (S+TG)**(-2) - 32*S*TG**(-1)*(S+TG)**(-1)*S4**(-1) - 40*S*
     +    TG**(-1)*S4**(-2) - 8*S*TG*(S+UG)**(-2)*S4**(-2) - 16*S*M2*
     +    (S+UG)**(-2)*S4**(-2) + 64*S*(S+TG)**(-2)*S4**(-1) + 32*S**2*
     +    TG**(-1)*(S+TG)**(-2)*S4**(-1) + 16*TG**(-2)*T1**(-1)*M2*MS2*
     +    S4**(-1) - 16*TG**(-2)*T1**(-1)*M2**2*MS2*S4**(-2) + 16*
     +    TG**(-2)*T1**(-1)*M2**2*S4**(-1) - 16*TG**(-2)*T1**(-1)*M2**3
     +    *S4**(-2) - 64*TG**(-2)*U1**(-4)*M2*MS2*S4**2 - 64*TG**(-2)*
     +    U1**(-4)*M2**2*S4**2 + 128*TG**(-2)*U1**(-3)*M2*MS2*S4 + 128*
     +    TG**(-2)*U1**(-3)*M2**2*S4 )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 192*
     +    TG**(-2)*U1**(-2)*M2*MS2 - 192*TG**(-2)*U1**(-2)*M2**2 + 128*
     +    TG**(-2)*U1**(-1)*M2*MS2*S4**(-1) + 128*TG**(-2)*U1**(-1)*
     +    M2**2*S4**(-1) - 64*TG**(-2)*M2*MS2*S4**(-2) + 16*TG**(-2)*
     +    M2**2*MS2*(S+UG)**(-1)*S4**(-2) - 64*TG**(-2)*M2**2*S4**(-2)
     +     + 16*TG**(-2)*M2**3*(S+UG)**(-1)*S4**(-2) - 16*TG**(-1)*
     +    T1**(-1)*M2*MS2*S4**(-2) + 16*TG**(-1)*T1**(-1)*M2*S4**(-1)
     +     - 32*TG**(-1)*T1**(-1)*M2**2*S4**(-2) + 64*TG**(-1)*U1**(-4)
     +    *M2*MS2*S4 - 64*TG**(-1)*U1**(-3)*M2*MS2 + 64*TG**(-1)*
     +    U1**(-3)*M2*S4 + 128*TG**(-1)*U1**(-2)*M2*MS2*S4**(-1) - 64*
     +    TG**(-1)*U1**(-2)*M2 - 32*TG**(-1)*U1**(-2)*S4 + 64*TG**(-1)*
     +    U1**(-1)*M2*MS2*(S+TG)**(-1)*S4**(-1) - 64*TG**(-1)*U1**(-1)*
     +    M2*MS2*S4**(-2) + 128*TG**(-1)*U1**(-1)*M2*S4**(-1) + 32*
     +    TG**(-1)*U1**(-1) + 16*TG**(-1)*M2*MS2*(S+UG)**(-1)*S4**(-2)
     +     + 64*TG**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 64*TG**(-1)*M2*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * ( 32*TG**(-1)*
     +    M2**2*(S+UG)**(-1)*S4**(-2) - 8*TG*T1**(-1)*S4**(-2) + 48*TG*
     +    (S+TG)**(-2)*S4**(-1) + 8*TG*(S+UG)**(-1)*S4**(-2) - 24*
     +    T1**(-1)*M2*S4**(-2) + 8*T1**(-1)*S4**(-1) - 32*U1**(-4)*M2*
     +    MS2 - 32*U1**(-3)*M2 - 32*U1**(-2)*M2*MS2*(S+TG)**(-1)*
     +    S4**(-1) - 32*U1**(-2)*M2*MS2*S4**(-2) + 16*U1**(-2) - 32*
     +    U1**(-1)*M2*MS2*(S+TG)**(-2)*S4**(-1) - 32*U1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 32*U1**(-1)*M2*S4**(-2) - 32*M2*
     +    (S+TG)**(-2)*S4**(-1) + 24*M2*(S+UG)**(-1)*S4**(-2) - 32*
     +    (S+TG)**(-2) - 32*(S+TG)**(-1)*S4**(-1) - 32*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(TG*UG-M2*S)**(-2)*(S4+MS2)*
     + TWO**(-1) * ( 48*S*TG*T1**(-2)*M2**2*MS2*S4**(-1) - 48*S*TG*
     +    T1**(-2)*M2**3*MS2*S4**(-2) - 24*S*TG*M2*MS2*S4**(-2) - 24*S*
     +    TG**2*T1**(-2)*M2**2*MS2*S4**(-2) - 24*S*T1**(-2)*M2**2*MS2
     +     + 48*S*T1**(-2)*M2**3*MS2*S4**(-1) - 24*S*T1**(-2)*M2**4*MS2
     +    *S4**(-2) - 24*S*M2**2*MS2*S4**(-2) + 24*S**2*TG*T1**(-2)*
     +    M2**2*MS2*S4**(-2) - 24*S**2*T1**(-2)*M2**2*MS2*S4**(-1) + 24
     +    *S**2*T1**(-2)*M2**3*MS2*S4**(-2) - 24*S**3*T1**(-1)*M2*MS2*
     +    S4**(-2) + 72*TG*T1**(-2)*M2**2*MS2 - 144*TG*T1**(-2)*M2**3*
     +    MS2*S4**(-1) + 72*TG*T1**(-2)*M2**4*MS2*S4**(-2) + 24*TG*
     +    M2**2*MS2*S4**(-2) - 72*TG**2*T1**(-2)*M2**2*MS2*S4**(-1) + 
     +    72*TG**2*T1**(-2)*M2**3*MS2*S4**(-2) + 24*TG**3*T1**(-2)*
     +    M2**2*MS2*S4**(-2) - 24*T1**(-2)*M2**2*MS2*S4 + 72*T1**(-2)*
     +    M2**3*MS2 - 72*T1**(-2)*M2**4*MS2*S4**(-1) + 24*T1**(-2)*
     +    M2**5*MS2*S4**(-2) - 24*M2**2*MS2*S4**(-1) + 24*M2**3*MS2*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * (  - 72*S**(-1)*TG*T1**(-2)*M2**2 + 144*S**(-1)*TG*
     +    T1**(-2)*M2**3*S4**(-1) - 72*S**(-1)*TG*T1**(-2)*M2**4*
     +    S4**(-2) - 24*S**(-1)*TG*M2**2*S4**(-2) + 72*S**(-1)*TG**2*
     +    T1**(-2)*M2**2*S4**(-1) - 72*S**(-1)*TG**2*T1**(-2)*M2**3*
     +    S4**(-2) - 24*S**(-1)*TG**3*T1**(-2)*M2**2*S4**(-2) + 24*
     +    S**(-1)*T1**(-2)*M2**2*S4 - 72*S**(-1)*T1**(-2)*M2**3 + 72*
     +    S**(-1)*T1**(-2)*M2**4*S4**(-1) - 24*S**(-1)*T1**(-2)*M2**5*
     +    S4**(-2) + 24*S**(-1)*M2**2*S4**(-1) - 24*S**(-1)*M2**3*
     +    S4**(-2) + 24*S*TG*T1**(-2)*M2*MS2*S4**(-2) - 24*S*T1**(-2)*
     +    M2*MS2*S4**(-1) + 24*S*T1**(-2)*M2**2*MS2*S4**(-2) - 24*S**2*
     +    T1**(-2)*M2*MS2*S4**(-2) + 16*TG**(-1)*T1**(-1)*M2*MS2 - 32*
     +    TG**(-1)*T1**(-1)*M2**2*MS2*S4**(-1) + 16*TG**(-1)*T1**(-1)*
     +    M2**3*MS2*S4**(-2) + 16*TG**(-1)*M2**2*MS2*S4**(-2) + 48*TG*
     +    T1**(-2)*M2*MS2*S4**(-1) - 48*TG*T1**(-2)*M2**2*MS2*S4**(-2)
     +     + 16*TG*T1**(-1)*M2*MS2*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * (  - 48*TG*T1**(-1)*M2*S4**(-1) + 48*TG*T1**(-1)*
     +    M2**2*S4**(-2) + 24*TG*M2*S4**(-2) - 24*TG**2*T1**(-2)*M2*MS2
     +    *S4**(-2) + 24*TG**2*T1**(-1)*M2*S4**(-2) - 24*T1**(-2)*M2*
     +    MS2 + 48*T1**(-2)*M2**2*MS2*S4**(-1) - 24*T1**(-2)*M2**3*MS2*
     +    S4**(-2) - 32*T1**(-1)*M2*MS2*S4**(-1) + 24*T1**(-1)*M2 + 32*
     +    T1**(-1)*M2**2*MS2*S4**(-2) - 48*T1**(-1)*M2**2*S4**(-1) + 24
     +    *T1**(-1)*M2**3*S4**(-2) - 8*M2*MS2*S4**(-2) + 24*M2**2*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * (  - 24*
     +    S**(-1)*TG*T1**(-2)*M2*S4**(-1) + 24*S**(-1)*TG*T1**(-2)*
     +    M2**2*S4**(-2) + 24*S**(-1)*TG*T1**(-1)*M2*S4**(-2) - 24*
     +    S**(-1)*TG*T1**(-1)*S4**(-1) + 16*S**(-1)*TG*(S+TG)**(-2) + 
     +    16*S**(-1)*TG*(S+TG)**(-1)*S4**(-1) + 28*S**(-1)*TG*S4**(-2)
     +     + 12*S**(-1)*TG**2*T1**(-2)*M2*S4**(-2) + 12*S**(-1)*TG**2*
     +    T1**(-1)*S4**(-2) - 16*S**(-1)*TG**2*(S+TG)**(-2)*S4**(-1) + 
     +    12*S**(-1)*T1**(-2)*M2 - 24*S**(-1)*T1**(-2)*M2**2*S4**(-1)
     +     + 12*S**(-1)*T1**(-2)*M2**3*S4**(-2) - 24*S**(-1)*T1**(-1)*
     +    M2*S4**(-1) + 12*S**(-1)*T1**(-1)*M2**2*S4**(-2) + 12*S**(-1)
     +    *T1**(-1) + 32*S**(-1)*U1**(-4)*M2**2*S4 - 32*S**(-1)*
     +    U1**(-3)*M2*S4 - 32*S**(-1)*U1**(-3)*M2**2 - 32*S**(-1)*
     +    U1**(-2)*M2*MS2*(S+TG)**(-1) + 32*S**(-1)*U1**(-2)*M2*MS2*
     +    S4**(-1) + 32*S**(-1)*U1**(-2)*M2 + 64*S**(-1)*U1**(-2)*M2**2
     +    *S4**(-1) + 16*S**(-1)*U1**(-2)*S4 - 32*S**(-1)*U1**(-1)*M2*
     +    MS2*(S+TG)**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * ( 64*S**(-1)*
     +    U1**(-1)*M2*MS2*(S+TG)**(-1)*S4**(-1) - 64*S**(-1)*U1**(-1)*
     +    M2*S4**(-1) + 32*S**(-1)*U1**(-1)*M2**2*(S+TG)**(-1)*S4**(-1)
     +     - 32*S**(-1)*U1**(-1)*M2**2*S4**(-2) - 16*S**(-1)*U1**(-1)
     +     + 32*S**(-1)*M2*MS2*(S+TG)**(-2)*S4**(-1) - 32*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) + 56*S**(-1)*M2*S4**(-2) + 8*S*TG**(-1)
     +    *M2**2*(S+UG)**(-2)*S4**(-2) + 8*S*TG**(-1)*S4**(-2) + 8*S*TG
     +    *(S+UG)**(-2)*S4**(-2) + 16*S*M2*(S+UG)**(-2)*S4**(-2) - 16*
     +    TG**(-2)*T1**(-1)*M2*MS2*S4**(-1) + 16*TG**(-2)*T1**(-1)*
     +    M2**2*MS2*S4**(-2) - 16*TG**(-2)*T1**(-1)*M2**2*S4**(-1) + 16
     +    *TG**(-2)*T1**(-1)*M2**3*S4**(-2) - 16*TG**(-2)*M2**2*MS2*
     +    (S+UG)**(-1)*S4**(-2) - 16*TG**(-2)*M2**3*(S+UG)**(-1)*
     +    S4**(-2) + 16*TG**(-1)*T1**(-1)*M2*MS2*S4**(-2) - 16*TG**(-1)
     +    *T1**(-1)*M2*S4**(-1) + 32*TG**(-1)*T1**(-1)*M2**2*S4**(-2)
     +     - 16*TG**(-1)*M2*MS2*(S+UG)**(-1)*S4**(-2) - 32*TG**(-1)*
     +    M2**2*(S+UG)**(-1)*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * ( 8*TG*
     +    T1**(-1)*S4**(-2) - 16*TG*(S+TG)**(-2)*S4**(-1) - 8*TG*
     +    (S+UG)**(-1)*S4**(-2) + 24*T1**(-1)*M2*S4**(-2) - 8*T1**(-1)*
     +    S4**(-1) + 32*U1**(-4)*M2*MS2 + 32*U1**(-3)*M2 + 32*U1**(-2)*
     +    M2*MS2*(S+TG)**(-1)*S4**(-1) + 32*U1**(-2)*M2*MS2*S4**(-2) - 
     +    16*U1**(-2) + 32*U1**(-1)*M2*MS2*(S+TG)**(-2)*S4**(-1) + 32*
     +    U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 32*U1**(-1)*M2*S4**(-2)
     +     + 32*M2*(S+TG)**(-2)*S4**(-1) - 24*M2*(S+UG)**(-1)*S4**(-2)
     +     )
     +
      M2QGH = M2QGH + COLO2(9)*CQED*(TG*UG-M2*S)**(-2)*(S4+MS2)*
     + TWO**(-1) * (  - 16*S*TG*T1**(-2)*M2**2*MS2*S4**(-1) + 16*S*TG*
     +    T1**(-2)*M2**3*MS2*S4**(-2) + 8*S*TG*M2*MS2*S4**(-2) + 8*S*
     +    TG**2*T1**(-2)*M2**2*MS2*S4**(-2) + 8*S*T1**(-2)*M2**2*MS2 - 
     +    16*S*T1**(-2)*M2**3*MS2*S4**(-1) + 8*S*T1**(-2)*M2**4*MS2*
     +    S4**(-2) + 8*S*M2**2*MS2*S4**(-2) - 8*S**2*TG*T1**(-2)*M2**2*
     +    MS2*S4**(-2) + 8*S**2*T1**(-2)*M2**2*MS2*S4**(-1) - 8*S**2*
     +    T1**(-2)*M2**3*MS2*S4**(-2) + 8*S**3*T1**(-1)*M2*MS2*S4**(-2)
     +     - 24*TG*T1**(-2)*M2**2*MS2 + 48*TG*T1**(-2)*M2**3*MS2*
     +    S4**(-1) - 24*TG*T1**(-2)*M2**4*MS2*S4**(-2) - 8*TG*M2**2*MS2
     +    *S4**(-2) + 24*TG**2*T1**(-2)*M2**2*MS2*S4**(-1) - 24*TG**2*
     +    T1**(-2)*M2**3*MS2*S4**(-2) - 8*TG**3*T1**(-2)*M2**2*MS2*
     +    S4**(-2) + 8*T1**(-2)*M2**2*MS2*S4 - 24*T1**(-2)*M2**3*MS2 + 
     +    24*T1**(-2)*M2**4*MS2*S4**(-1) - 8*T1**(-2)*M2**5*MS2*
     +    S4**(-2) + 8*M2**2*MS2*S4**(-1) - 8*M2**3*MS2*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*CQED*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * ( 24*S**(-1)*TG*T1**(-2)*M2**2 - 48*S**(-1)*TG*
     +    T1**(-2)*M2**3*S4**(-1) + 24*S**(-1)*TG*T1**(-2)*M2**4*
     +    S4**(-2) + 8*S**(-1)*TG*M2**2*S4**(-2) - 24*S**(-1)*TG**2*
     +    T1**(-2)*M2**2*S4**(-1) + 24*S**(-1)*TG**2*T1**(-2)*M2**3*
     +    S4**(-2) + 8*S**(-1)*TG**3*T1**(-2)*M2**2*S4**(-2) - 8*
     +    S**(-1)*T1**(-2)*M2**2*S4 + 24*S**(-1)*T1**(-2)*M2**3 - 24*
     +    S**(-1)*T1**(-2)*M2**4*S4**(-1) + 8*S**(-1)*T1**(-2)*M2**5*
     +    S4**(-2) - 8*S**(-1)*M2**2*S4**(-1) + 8*S**(-1)*M2**3*
     +    S4**(-2) - 8*S*TG*T1**(-2)*M2*MS2*S4**(-2) + 8*S*T1**(-2)*M2*
     +    MS2*S4**(-1) - 8*S*T1**(-2)*M2**2*MS2*S4**(-2) + 8*S**2*
     +    T1**(-2)*M2*MS2*S4**(-2) - 16*TG*T1**(-2)*M2*MS2*S4**(-1) + 
     +    16*TG*T1**(-2)*M2**2*MS2*S4**(-2) + 16*TG*T1**(-1)*M2*
     +    S4**(-1) - 16*TG*T1**(-1)*M2**2*S4**(-2) - 8*TG*M2*S4**(-2)
     +     + 8*TG**2*T1**(-2)*M2*MS2*S4**(-2) - 8*TG**2*T1**(-1)*M2*
     +    S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*CQED*(TG*UG-M2*S)**(-1)*(S4+MS2)*
     + TWO**(-1) * ( 8*T1**(-2)*M2*MS2 - 16*T1**(-2)*M2**2*MS2*S4**(-1)
     +     + 8*T1**(-2)*M2**3*MS2*S4**(-2) - 8*T1**(-1)*M2 + 16*
     +    T1**(-1)*M2**2*S4**(-1) - 8*T1**(-1)*M2**3*S4**(-2) + 8*M2*
     +    MS2*S4**(-2) - 8*M2**2*S4**(-2) )
     +
      M2QGH = M2QGH + COLO2(9)*CQED*(S4+MS2)*TWO**(-1) * ( 8*S**(-1)*TG
     +    *T1**(-2)*M2*S4**(-1) - 8*S**(-1)*TG*T1**(-2)*M2**2*S4**(-2)
     +     - 8*S**(-1)*TG*T1**(-1)*M2*S4**(-2) + 8*S**(-1)*TG*T1**(-1)*
     +    S4**(-1) - 4*S**(-1)*TG*S4**(-2) - 4*S**(-1)*TG**2*T1**(-2)*
     +    M2*S4**(-2) - 4*S**(-1)*TG**2*T1**(-1)*S4**(-2) - 4*S**(-1)*
     +    T1**(-2)*M2 + 8*S**(-1)*T1**(-2)*M2**2*S4**(-1) - 4*S**(-1)*
     +    T1**(-2)*M2**3*S4**(-2) + 8*S**(-1)*T1**(-1)*M2*S4**(-1) - 4*
     +    S**(-1)*T1**(-1)*M2**2*S4**(-2) - 4*S**(-1)*T1**(-1) - 8*
     +    S**(-1)*M2*S4**(-2) )


      M2QGH = 2.D0 * M2QGH

      Q2MS2 = (MS + MG)**2/4.D0/MS**2

      DSGQGH = S4/(S4+MS2)/2.D0 *
     +     ALPHAS**3 * AVG * M2QGH /4.D0 /S**2 *CONV
     +     +  DSGQG3(ALPHAS,S,TG,S4,MS,MG,Q2MS2) 

C***  CHANGES TO THE SCALE Q**2 = (MS+MG)**2/4

      RETURN
      END


      REAL*8 FUNCTION DSGQGD(ALPHAS,S,TG,S4,MS,MG,DEL,S4MAX)
C***  LOG(DEL) PART OF CROSS SECTIONS FOR Q +G -> SQ + GL +G 
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,MS,MG,MS2,MG2,M2,CO,CK,CQED
      REAL*8 M2QGD,NS,CONV,N,AVG,XLAM,BETA,XSG,DEL,S4MAX,DLDEL1,DLDEL2
      REAL*8 COLO2(1:9),Q2MS2,DSGQG2

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CO = (N**2-1.D0)*N
      CK = (N**2-1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /N /(N**2 -1.D0)

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      T1 = TG +M2

      XLAM = SQRT((S -MS2 -MG2)**2 -4*MS2*MG2)
      BETA = SQRT(1.D0 -4*MG*MS/(S-(MG-MS)**2))
      XSG = (1.D0 -BETA)/(1.D0 +BETA)

      DLDEL1 = LOG(S4MAX/MS2) -(S4MAX-DEL)/S4
      DLDEL2 = LOG(S4MAX/MS2)**2 -2.D0*(S4MAX -DEL)/S4*LOG(S4/MS2)

      COLO2(1) = LOG(XSG)
      COLO2(2) = LOG(MG2/MS2)
      COLO2(3) = LOG(S/MS2)
      COLO2(4) = LOG(-TG/MS2)
      COLO2(5) = LOG(-T1/MS2)
      COLO2(6) = LOG((S+T1)/MS2)
      COLO2(7) = LOG((S+TG)/MS2)


      M2QGD = 0.D0
      M2QGD = M2QGD + N*CO*DLDEL1 * ( 12 - 24*S**(-1)*TG**(-1)*MS2*MG2
     +     + 12*S**(-1)*TG**(-1)*MS2**2 + 12*S**(-1)*TG**(-1)*MG2**2 + 
     +    6*S**(-1)*T1 - 6*S**(-1)*MS2 + 6*S**(-1)*MG2 - 12*S*TG**(-1)*
     +    MS2*(S+TG)**(-1) + 12*S*TG**(-1)*MG2*(S+TG)**(-1) + 12*S*
     +    TG**(-1) - 24*TG**(-2)*MS2*MG2 + 24*TG**(-2)*MG2**2 - 12*
     +    TG**(-1)*MS2 + 12*TG**(-1)*MS2**2*(S+TG)**(-1) + 12*TG**(-1)*
     +    MG2 - 12*TG**(-1)*MG2**2*(S+TG)**(-1) + 12*MS2*MG2*
     +    (S+TG)**(-2) - 12*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + N*CO*DLDEL2 * (  - 12 + 24*S**(-1)*TG**(-1)*MS2*
     +    MG2 - 12*S**(-1)*TG**(-1)*MS2**2 - 12*S**(-1)*TG**(-1)*MG2**2
     +     - 6*S**(-1)*T1 + 6*S**(-1)*MS2 - 6*S**(-1)*MG2 + 12*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) - 12*S*TG**(-1)*MG2*(S+TG)**(-1) - 
     +    12*S*TG**(-1) + 24*TG**(-2)*MS2*MG2 - 24*TG**(-2)*MG2**2 + 12
     +    *TG**(-1)*MS2 - 12*TG**(-1)*MS2**2*(S+TG)**(-1) - 12*TG**(-1)
     +    *MG2 + 12*TG**(-1)*MG2**2*(S+TG)**(-1) - 12*MS2*MG2*
     +    (S+TG)**(-2) + 12*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + N*CK*DLDEL1 * (  - 4 + 8*S**(-1)*TG**(-1)*MS2*MG2
     +     - 4*S**(-1)*TG**(-1)*MS2**2 - 4*S**(-1)*TG**(-1)*MG2**2 + 16
     +    *S**(-1)*T1*MS2*(S+TG)**(-1) - 16*S**(-1)*T1*MG2*(S+TG)**(-1)
     +     - 10*S**(-1)*T1 - 6*S**(-1)*MS2 + 6*S**(-1)*MG2 + 4*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) - 4*S*TG**(-1)*MG2*(S+TG)**(-1) - 4
     +    *S*TG**(-1) + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2 + 4*
     +    TG**(-1)*MS2 - 4*TG**(-1)*MS2**2*(S+TG)**(-1) - 4*TG**(-1)*
     +    MG2 + 4*TG**(-1)*MG2**2*(S+TG)**(-1) - 20*MS2*MG2*
     +    (S+TG)**(-2) + 20*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + N*CK*DLDEL2 * ( 4 - 8*S**(-1)*TG**(-1)*MS2*MG2 + 
     +    4*S**(-1)*TG**(-1)*MS2**2 + 4*S**(-1)*TG**(-1)*MG2**2 - 16*
     +    S**(-1)*T1*MS2*(S+TG)**(-1) + 16*S**(-1)*T1*MG2*(S+TG)**(-1)
     +     + 10*S**(-1)*T1 + 6*S**(-1)*MS2 - 6*S**(-1)*MG2 - 4*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) + 4*S*TG**(-1)*MG2*(S+TG)**(-1) + 4
     +    *S*TG**(-1) - 4*S*T1**(-1)*MS2*(S+TG)**(-1) + 4*S*T1**(-1)*
     +    MG2*(S+TG)**(-1) - 8*TG**(-2)*MS2*MG2 + 8*TG**(-2)*MG2**2 - 4
     +    *TG**(-1)*MS2 + 4*TG**(-1)*MS2**2*(S+TG)**(-1) + 4*TG**(-1)*
     +    MG2 - 4*TG**(-1)*MG2**2*(S+TG)**(-1) + 8*T1**(-1)*MS2*MG2*
     +    (S+TG)**(-1) + 4*T1**(-1)*MS2 - 4*T1**(-1)*MS2**2*
     +    (S+TG)**(-1) - 4*T1**(-1)*MG2 - 4*T1**(-1)*MG2**2*
     +    (S+TG)**(-1) + 20*MS2*MG2*(S+TG)**(-2) - 4*MS2*(S+TG)**(-1)
     +     - 20*MS2**2*(S+TG)**(-2) + 4*MG2*(S+TG)**(-1) )
     +
      M2QGD = M2QGD + CQED*DLDEL1 * (  - 4*S**(-1)*T1*MS2*(S+TG)**(-1)
     +     + 4*S**(-1)*T1*MG2*(S+TG)**(-1) + 2*S**(-1)*T1 + 2*S**(-1)*
     +    MS2 - 2*S**(-1)*MG2 + 4*MS2*MG2*(S+TG)**(-2) - 4*MS2**2*
     +    (S+TG)**(-2) )
     +
      M2QGD = M2QGD + CQED*DLDEL2 * ( 4*S**(-1)*T1*MS2*(S+TG)**(-1) - 4
     +    *S**(-1)*T1*MG2*(S+TG)**(-1) - 2*S**(-1)*T1 - 2*S**(-1)*MS2
     +     + 2*S**(-1)*MG2 - 4*MS2*MG2*(S+TG)**(-2) + 4*MS2**2*
     +    (S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(1)*N*CO*DLDEL1 * ( 24*S**(-1)*TG**(-1)*
     +    XLAM**(-3)*MS2*MG2**4 - 16*S**(-1)*TG**(-1)*XLAM**(-3)*MS2**2
     +    *MG2**3 - 16*S**(-1)*TG**(-1)*XLAM**(-3)*MS2**3*MG2**2 + 24*
     +    S**(-1)*TG**(-1)*XLAM**(-3)*MS2**4*MG2 - 8*S**(-1)*TG**(-1)*
     +    XLAM**(-3)*MS2**5 - 8*S**(-1)*TG**(-1)*XLAM**(-3)*MG2**5 + 4*
     +    S**(-1)*T1*XLAM**(-3)*MS2*MG2**2 + 4*S**(-1)*T1*XLAM**(-3)*
     +    MS2**2*MG2 - 4*S**(-1)*T1*XLAM**(-3)*MS2**3 - 4*S**(-1)*T1*
     +    XLAM**(-3)*MG2**3 + 8*S**(-1)*XLAM**(-3)*MS2*MG2**3 - 8*
     +    S**(-1)*XLAM**(-3)*MS2**3*MG2 + 4*S**(-1)*XLAM**(-3)*MS2**4
     +     - 4*S**(-1)*XLAM**(-3)*MG2**4 - 8*S*TG**(-2)*XLAM**(-3)*MS2*
     +    MG2**3 + 8*S*TG**(-2)*XLAM**(-3)*MS2**2*MG2**2 - 24*S*
     +    TG**(-2)*XLAM**(-3)*MS2**3*MG2 + 24*S*TG**(-2)*XLAM**(-3)*
     +    MG2**4 + 20*S*TG**(-1)*XLAM**(-3)*MS2*MG2**2 + 36*S*TG**(-1)*
     +    XLAM**(-3)*MS2**2*MG2 - 52*S*TG**(-1)*XLAM**(-3)*MS2**3 - 4*S
     +    *TG**(-1)*XLAM**(-3)*MG2**3 - 12*S*T1*XLAM**(-3)*MS2 - 12*S*
     +    T1*XLAM**(-3)*MG2 )
     +
      M2QGD = M2QGD + COLO2(1)*N*CO*DLDEL1 * ( 16*S*XLAM**(-3)*MS2*MG2
     +     + 36*S*XLAM**(-3)*MS2**2 + 12*S*XLAM**(-3)*MG2**2 + 24*S**2*
     +    TG**(-2)*XLAM**(-3)*MS2**2*MG2 - 24*S**2*TG**(-2)*XLAM**(-3)*
     +    MG2**3 - 8*S**2*TG**(-1)*XLAM**(-3)*MS2*MG2 + 44*S**2*
     +    TG**(-1)*XLAM**(-3)*MS2**2 - 4*S**2*TG**(-1)*XLAM**(-3)*
     +    MG2**2 + 4*S**2*T1*XLAM**(-3) - 28*S**2*XLAM**(-3)*MS2 - 20*
     +    S**2*XLAM**(-3)*MG2 - 8*S**3*TG**(-2)*XLAM**(-3)*MS2*MG2 + 8*
     +    S**3*TG**(-2)*XLAM**(-3)*MG2**2 - 20*S**3*TG**(-1)*XLAM**(-3)
     +    *MS2 - 4*S**3*TG**(-1)*XLAM**(-3)*MG2 + 8*S**3*XLAM**(-3) + 4
     +    *S**4*TG**(-1)*XLAM**(-3) + 16*TG**(-2)*XLAM**(-3)*MS2*MG2**4
     +     - 16*TG**(-2)*XLAM**(-3)*MS2**3*MG2**2 + 8*TG**(-2)*
     +    XLAM**(-3)*MS2**4*MG2 - 8*TG**(-2)*XLAM**(-3)*MG2**5 - 16*
     +    TG**(-1)*XLAM**(-3)*MS2*MG2**3 + 16*TG**(-1)*XLAM**(-3)*
     +    MS2**2*MG2**2 - 48*TG**(-1)*XLAM**(-3)*MS2**3*MG2 + 32*
     +    TG**(-1)*XLAM**(-3)*MS2**4 + 16*TG**(-1)*XLAM**(-3)*MG2**4 + 
     +    8*T1*XLAM**(-3)*MS2*MG2 )
     +
      M2QGD = M2QGD + COLO2(1)*N*CO*DLDEL1 * ( 12*T1*XLAM**(-3)*MS2**2
     +     + 12*T1*XLAM**(-3)*MG2**2 + 4*XLAM**(-3)*MS2*MG2**2 + 12*
     +    XLAM**(-3)*MS2**2*MG2 - 20*XLAM**(-3)*MS2**3 + 4*XLAM**(-3)*
     +    MG2**3 )
     +
      M2QGD = M2QGD + COLO2(1)*N*CK*DLDEL1 * ( 4*S**(-1)*T1*XLAM**(-3)*
     +    MS2*MG2**2 + 4*S**(-1)*T1*XLAM**(-3)*MS2**2*MG2 - 4*S**(-1)*
     +    T1*XLAM**(-3)*MS2**3 - 4*S**(-1)*T1*XLAM**(-3)*MG2**3 + 8*
     +    S**(-1)*T1**2*XLAM**(-3)*MS2**2 - 8*S**(-1)*T1**2*XLAM**(-3)*
     +    MG2**2 - 16*S**(-1)*T1**3*XLAM**(-3)*MS2*MG2*(S+TG)**(-1) - 8
     +    *S**(-1)*T1**3*XLAM**(-3)*MS2 + 8*S**(-1)*T1**3*XLAM**(-3)*
     +    MG2 + 16*S**(-1)*T1**3*XLAM**(-3)*MG2**2*(S+TG)**(-1) + 8*
     +    S**(-1)*T1**4*XLAM**(-3)*MS2*(S+TG)**(-1) - 8*S**(-1)*T1**4*
     +    XLAM**(-3)*MG2*(S+TG)**(-1) + 8*S**(-1)*XLAM**(-3)*MS2*MG2**3
     +     - 8*S**(-1)*XLAM**(-3)*MS2**3*MG2 + 4*S**(-1)*XLAM**(-3)*
     +    MS2**4 - 4*S**(-1)*XLAM**(-3)*MG2**4 - 40*S*T1*XLAM**(-3)*MS2
     +    *MG2*(S+TG)**(-1) + 272*S*T1*XLAM**(-3)*MS2*MG2**2*
     +    (S+TG)**(-2) - 92*S*T1*XLAM**(-3)*MS2 - 136*S*T1*XLAM**(-3)*
     +    MS2**2*MG2*(S+TG)**(-2) - 48*S*T1*XLAM**(-3)*MS2**2*
     +    (S+TG)**(-1) + 116*S*T1*XLAM**(-3)*MG2 + 40*S*T1*XLAM**(-3)*
     +    MG2**2*(S+TG)**(-1) )
     +
      M2QGD = M2QGD + COLO2(1)*N*CK*DLDEL1 * (  - 136*S*T1*XLAM**(-3)*
     +    MG2**3*(S+TG)**(-2) - 96*S*T1**2*XLAM**(-3)*MS2*MG2*
     +    (S+TG)**(-2) + 144*S*T1**2*XLAM**(-3)*MS2*(S+TG)**(-1) + 48*S
     +    *T1**2*XLAM**(-3)*MS2**2*(S+TG)**(-2) - 144*S*T1**2*
     +    XLAM**(-3)*MG2*(S+TG)**(-1) + 48*S*T1**2*XLAM**(-3)*MG2**2*
     +    (S+TG)**(-2) - 32*S*XLAM**(-3)*MS2*MG2 - 336*S*XLAM**(-3)*MS2
     +    *MG2**2*(S+TG)**(-1) - 224*S*XLAM**(-3)*MS2*MG2**3*
     +    (S+TG)**(-2) + 88*S*XLAM**(-3)*MS2**2*MG2*(S+TG)**(-1) + 144*
     +    S*XLAM**(-3)*MS2**2*MG2**2*(S+TG)**(-2) + 52*S*XLAM**(-3)*
     +    MS2**2 + 28*S*XLAM**(-3)*MG2**2 + 248*S*XLAM**(-3)*MG2**3*
     +    (S+TG)**(-1) + 80*S*XLAM**(-3)*MG2**4*(S+TG)**(-2) - 192*S**2
     +    *T1*XLAM**(-3)*MS2*MG2*(S+TG)**(-2) + 184*S**2*T1*XLAM**(-3)*
     +    MS2*(S+TG)**(-1) + 96*S**2*T1*XLAM**(-3)*MS2**2*(S+TG)**(-2)
     +     - 176*S**2*T1*XLAM**(-3)*MG2*(S+TG)**(-1) + 96*S**2*T1*
     +    XLAM**(-3)*MG2**2*(S+TG)**(-2) - 4*S**2*T1*XLAM**(-3) + 8*
     +    S**2*XLAM**(-3)*MS2*MG2*(S+TG)**(-1) )
     +
      M2QGD = M2QGD + COLO2(1)*N*CK*DLDEL1 * ( 320*S**2*XLAM**(-3)*MS2*
     +    MG2**2*(S+TG)**(-2) - 76*S**2*XLAM**(-3)*MS2 - 184*S**2*
     +    XLAM**(-3)*MS2**2*MG2*(S+TG)**(-2) - 24*S**2*XLAM**(-3)*
     +    MS2**2*(S+TG)**(-1) + 68*S**2*XLAM**(-3)*MG2 - 32*S**2*
     +    XLAM**(-3)*MG2**2*(S+TG)**(-1) - 136*S**2*XLAM**(-3)*MG2**3*
     +    (S+TG)**(-2) - 120*S**3*XLAM**(-3)*MS2*MG2*(S+TG)**(-2) + 72*
     +    S**3*XLAM**(-3)*MS2*(S+TG)**(-1) + 64*S**3*XLAM**(-3)*MS2**2*
     +    (S+TG)**(-2) - 64*S**3*XLAM**(-3)*MG2*(S+TG)**(-1) + 56*S**3*
     +    XLAM**(-3)*MG2**2*(S+TG)**(-2) - 64*T1*XLAM**(-3)*MS2*MG2 - 
     +    224*T1*XLAM**(-3)*MS2*MG2**2*(S+TG)**(-1) - 32*T1*XLAM**(-3)*
     +    MS2*MG2**3*(S+TG)**(-2) + 88*T1*XLAM**(-3)*MS2**2*MG2*
     +    (S+TG)**(-1) + 16*T1*XLAM**(-3)*MS2**2*MG2**2*(S+TG)**(-2) + 
     +    36*T1*XLAM**(-3)*MS2**2 - 4*T1*XLAM**(-3)*MG2**2 + 136*T1*
     +    XLAM**(-3)*MG2**3*(S+TG)**(-1) + 16*T1*XLAM**(-3)*MG2**4*
     +    (S+TG)**(-2) - 8*T1**2*XLAM**(-3)*MS2*MG2*(S+TG)**(-1) + 48*
     +    T1**2*XLAM**(-3)*MS2*MG2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(1)*N*CK*DLDEL1 * (  - 40*T1**2*XLAM**(-3)*
     +    MS2 - 24*T1**2*XLAM**(-3)*MS2**2*MG2*(S+TG)**(-2) - 24*T1**2*
     +    XLAM**(-3)*MS2**2*(S+TG)**(-1) + 40*T1**2*XLAM**(-3)*MG2 + 32
     +    *T1**2*XLAM**(-3)*MG2**2*(S+TG)**(-1) - 24*T1**2*XLAM**(-3)*
     +    MG2**3*(S+TG)**(-2) - 16*T1**3*XLAM**(-3)*MS2*MG2*
     +    (S+TG)**(-2) + 48*T1**3*XLAM**(-3)*MS2*(S+TG)**(-1) + 8*T1**3
     +    *XLAM**(-3)*MS2**2*(S+TG)**(-2) - 48*T1**3*XLAM**(-3)*MG2*
     +    (S+TG)**(-1) + 8*T1**3*XLAM**(-3)*MG2**2*(S+TG)**(-2) + 132*
     +    XLAM**(-3)*MS2*MG2**2 + 256*XLAM**(-3)*MS2*MG2**3*
     +    (S+TG)**(-1) + 4*XLAM**(-3)*MS2**2*MG2 - 128*XLAM**(-3)*
     +    MS2**2*MG2**2*(S+TG)**(-1) - 20*XLAM**(-3)*MS2**3 - 116*
     +    XLAM**(-3)*MG2**3 - 128*XLAM**(-3)*MG2**4*(S+TG)**(-1) )
     +
      M2QGD = M2QGD + COLO2(2)*N*CO*DLDEL1 * ( 4 - 8*S**(-1)*TG**(-1)*
     +    MS2*MG2 + 4*S**(-1)*TG**(-1)*MS2**2 + 4*S**(-1)*TG**(-1)*
     +    MG2**2 + 2*S**(-1)*T1 - 2*S**(-1)*MS2 + 2*S**(-1)*MG2 - 8*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) + 8*S*TG**(-1)*MG2*(S+TG)**(-1) + 6
     +    *S*TG**(-1) - 12*TG**(-2)*MS2*MG2 + 12*TG**(-2)*MG2**2 - 4*
     +    TG**(-1)*MS2 + 8*TG**(-1)*MS2**2*(S+TG)**(-1) + 4*TG**(-1)*
     +    MG2 - 8*TG**(-1)*MG2**2*(S+TG)**(-1) + 8*MS2*MG2*(S+TG)**(-2)
     +     - 8*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(2)*N*CK*DLDEL1 * ( 4*S**(-1)*T1*MS2*
     +    (S+TG)**(-1) - 4*S**(-1)*T1*MG2*(S+TG)**(-1) - 2*S**(-1)*T1
     +     - 2*S**(-1)*MS2 + 2*S**(-1)*MG2 - 4*MS2*MG2*(S+TG)**(-2) + 4
     +    *MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(3)*N*CO*DLDEL1 * (  - 8 + 16*S**(-1)*
     +    TG**(-1)*MS2*MG2 - 8*S**(-1)*TG**(-1)*MS2**2 - 8*S**(-1)*
     +    TG**(-1)*MG2**2 - 4*S**(-1)*T1 + 4*S**(-1)*MS2 - 4*S**(-1)*
     +    MG2 - 4*S*TG**(-1) + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2
     +     + 8*TG**(-1)*MS2 - 8*TG**(-1)*MG2 )
     +
      M2QGD = M2QGD + COLO2(3)*N*CK*DLDEL1 * (  - 8*S**(-1)*T1*MS2*
     +    (S+TG)**(-1) + 8*S**(-1)*T1*MG2*(S+TG)**(-1) + 4*S**(-1)*T1
     +     + 4*S**(-1)*MS2 - 4*S**(-1)*MG2 + 8*MS2*MG2*(S+TG)**(-2) - 8
     +    *MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(4)*N*CO*DLDEL1 * (  - 8 + 16*S**(-1)*
     +    TG**(-1)*MS2*MG2 - 8*S**(-1)*TG**(-1)*MS2**2 - 8*S**(-1)*
     +    TG**(-1)*MG2**2 - 4*S**(-1)*T1 + 4*S**(-1)*MS2 - 4*S**(-1)*
     +    MG2 + 8*S*TG**(-1)*MS2*(S+TG)**(-1) - 8*S*TG**(-1)*MG2*
     +    (S+TG)**(-1) - 8*S*TG**(-1) + 16*TG**(-2)*MS2*MG2 - 16*
     +    TG**(-2)*MG2**2 + 8*TG**(-1)*MS2 - 8*TG**(-1)*MS2**2*
     +    (S+TG)**(-1) - 8*TG**(-1)*MG2 + 8*TG**(-1)*MG2**2*
     +    (S+TG)**(-1) - 8*MS2*MG2*(S+TG)**(-2) + 8*MS2**2*(S+TG)**(-2)
     +     )
     +
      M2QGD = M2QGD + COLO2(5)*N*CO*DLDEL1 * ( 8 - 16*S**(-1)*TG**(-1)*
     +    MS2*MG2 + 8*S**(-1)*TG**(-1)*MS2**2 + 8*S**(-1)*TG**(-1)*
     +    MG2**2 + 4*S**(-1)*T1 - 4*S**(-1)*MS2 + 4*S**(-1)*MG2 - 8*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) + 8*S*TG**(-1)*MG2*(S+TG)**(-1) + 8
     +    *S*TG**(-1) - 16*TG**(-2)*MS2*MG2 + 16*TG**(-2)*MG2**2 - 8*
     +    TG**(-1)*MS2 + 8*TG**(-1)*MS2**2*(S+TG)**(-1) + 8*TG**(-1)*
     +    MG2 - 8*TG**(-1)*MG2**2*(S+TG)**(-1) + 8*MS2*MG2*(S+TG)**(-2)
     +     - 8*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(5)*N*CK*DLDEL1 * ( 16*S**(-1)*T1*MS2*
     +    (S+TG)**(-1) - 16*S**(-1)*T1*MG2*(S+TG)**(-1) - 8*S**(-1)*T1
     +     - 8*S**(-1)*MS2 + 8*S**(-1)*MG2 - 16*MS2*MG2*(S+TG)**(-2) + 
     +    16*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(6)*N*CO*DLDEL1 * ( 8*S*TG**(-1)*MS2*
     +    (S+TG)**(-1) - 8*S*TG**(-1)*MG2*(S+TG)**(-1) - 4*S*TG**(-1)
     +     + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2 - 8*TG**(-1)*MS2**2
     +    *(S+TG)**(-1) + 8*TG**(-1)*MG2**2*(S+TG)**(-1) - 8*MS2*MG2*
     +    (S+TG)**(-2) + 8*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(6)*N*CK*DLDEL1 * (  - 8*S**(-1)*T1*MS2*
     +    (S+TG)**(-1) + 8*S**(-1)*T1*MG2*(S+TG)**(-1) + 4*S**(-1)*T1
     +     + 4*S**(-1)*MS2 - 4*S**(-1)*MG2 + 8*MS2*MG2*(S+TG)**(-2) - 8
     +    *MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(7)*N*CO*DLDEL1 * ( 16 - 32*S**(-1)*TG**(-1)
     +    *MS2*MG2 + 16*S**(-1)*TG**(-1)*MS2**2 + 16*S**(-1)*TG**(-1)*
     +    MG2**2 + 8*S**(-1)*T1 - 8*S**(-1)*MS2 + 8*S**(-1)*MG2 - 8*S*
     +    TG**(-1)*MS2*(S+TG)**(-1) + 8*S*TG**(-1)*MG2*(S+TG)**(-1) + 
     +    12*S*TG**(-1) - 24*TG**(-2)*MS2*MG2 + 24*TG**(-2)*MG2**2 - 16
     +    *TG**(-1)*MS2 + 8*TG**(-1)*MS2**2*(S+TG)**(-1) + 16*TG**(-1)*
     +    MG2 - 8*TG**(-1)*MG2**2*(S+TG)**(-1) + 8*MS2*MG2*(S+TG)**(-2)
     +     - 8*MS2**2*(S+TG)**(-2) )
     +
      M2QGD = M2QGD + COLO2(7)*N*CK*DLDEL1 * ( 8*S**(-1)*T1*MS2*
     +    (S+TG)**(-1) - 8*S**(-1)*T1*MG2*(S+TG)**(-1) - 4*S**(-1)*T1
     +     - 4*S**(-1)*MS2 + 4*S**(-1)*MG2 - 8*MS2*MG2*(S+TG)**(-2) + 8
     +    *MS2**2*(S+TG)**(-2) )


      M2QGD = 2.D0 * M2QGD

      Q2MS2 = (MS + MG)**2/4.D0/MS**2
           
      DSGQGD = ALPHAS**3 * AVG * M2QGD /4.D0 /S**2 *CONV
     +     + DSGQG2(ALPHAS,S,TG,S4,MS,MG,DEL,S4MAX,Q2MS2)

C***  CHANGES TO THE SCALE Q**2 = (MS+MG)**2/4

      RETURN
      END



      REAL*8 FUNCTION DSGGG3(ALPHAS,S,TG,S4,MS,MG,SCA)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,TG,U1,UG,S4,MS,MG,SCA,X1,X2
      REAL*8 PI,N,TF,NS, DSGQGB, PQG1, PQG2

      N = 3.D0
      TF = 0.5D0
      NS = 6.D0
      PI = 4.D0*ATAN(1.D0)
      U1 = S4 -S -TG
      T1 = TG +MG**2 -MS**2
      UG = U1 -MG**2 +MS**2
      X1 = -T1/(S+UG)
      X2 = -U1/(S+TG)
      PQG1 = TF*(X1**2 + (1.D0 -X1)**2)
      PQG2 = TF*(X2**2 + (1.D0 -X2)**2)
      DSGGG3 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) * (NS - 1.D0) *
     +     (-1.D0/T1*DSGQGB(ALPHAS,X1*S,TG,MS,MG)*PQG1*X1**2
     +      -1.D0/U1*DSGQGB(ALPHAS,X2*S,UG,MS,MG)*PQG2*X2**2 )

      RETURN
      END


      REAL*8 FUNCTION DSGGGH(ALPHAS,S,TG,S4,MS,MG,EPSS,EPSG)
C***  CROSS SECTIONS FOR G + G -> SQ + GL +QB
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,UG,MS,MG,MS2,MG2,M2,EPSS,EPSG
      REAL*8 NS,CONV,N,CO,CK,CQED,AVG,S4G, S4G2,TWO,DEL,SYMBT,SYMBU
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 AHP1P1, ABP1P1, ABP1M1, ABP1P2, A4P1P2, A4P1P1
      REAL*8 A4M1P2, A4M1P1, A4P0P2, A4P1P0, A4M2P1, A4P2M2
      REAL*8 ABP1M2,A4M1P0
      REAL*8 C4P1P2, C4P1P1, C4M1P1, C4P1P0, CBP1P1, C4M2P1
      REAL*8 M2GGH, ANG2(1:72), COLO2(1:9),DSGGG3,Q2MS2
           
      CONV = 389379660.D0
      NS = 6.D0

      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0)**2

      TWO = 2.D0

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2
      S4G2 = S4G**2 + EPSG*MG**4
      SYMBT = (M2*S-TG*S4G)/((M2*S-TG*S4G)**2 +(S+TG)*EPSS*MS**4 )
      SYMBU = (M2*S-UG*S4G)/((M2*S-UG*S4G)**2 +(S+UG)*EPSS*MS**4 )


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6) + M2
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9) 
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) + MS2 +MG2
      ANB(1,2) = -ANB(1,1)
      ANC(1,2) = -ANC(1,1)
      ANA(1,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(1,3) = -ANA(1,3)
      ANC(1,3) =  0.D0
      ANA(1,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(1,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)-2.D0*ANGDEF(3)*ANGDEF(4)
      ANC(1,4) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(1,5) = -ANB(1,3)
      ANC(1,5) = -ANC(1,3)
      ANA(1,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(1,6) = -ANB(1,4)
      ANC(1,6) = -ANC(1,4)
      ANA(1,7) = +ANA(1,1) -M2
      ANB(1,7) = +ANB(1,1)
      ANC(1,7) = +ANC(1,1)
      ANA(1,8) = +ANA(1,5) -M2
      ANB(1,8) = +ANB(1,5)
      ANC(1,8) = +ANC(1,5)
      ANA(1,9) = +ANA(1,6) -M2
      ANB(1,9) = +ANB(1,6)
      ANC(1,9) = +ANC(1,6)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(3,2) = -ANB(3,1)
      ANC(3,2) = -ANC(3,1)
      ANA(3,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(3,3) = 
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10) -2.D0*ANGDEF(2)*ANGDEF(4)
      ANC(3,3) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(3,4) = -ANA(3,4)
      ANC(3,4) =  0.D0
      ANA(3,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(3,5) = -ANB(3,3)
      ANC(3,5) = -ANC(3,3)
      ANA(3,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(3,6) = -ANB(3,4)
      ANC(3,6) = -ANC(3,4)
      ANA(3,7) = +ANA(3,1) -M2
      ANB(3,7) = +ANB(3,1)
      ANC(3,7) = +ANC(3,1)
      ANA(3,8) = +ANA(3,5) -M2
      ANB(3,8) = +ANB(3,5)
      ANC(3,8) = +ANC(3,5)
      ANA(3,9) = +ANA(3,6) -M2
      ANB(3,9) = +ANB(3,6)
      ANC(3,9) = +ANC(3,6)

C$$$      ANG2(1) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(2) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(3) = A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(4) = A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(5) = A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(6) = A4M1P1(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(7) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(8) = A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(9) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(10)= A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))

C$$$      ANG2(11)= A4P0P2(ANA(1,7),ANB(1,7),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(12)= A4P1P0(ANA(2,7),ANB(2,7),ANA(2,7),ANB(2,7),ANC(2,7))
      ANG2(13)= A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(14)= A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(15)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(16)= A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(17)= ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(18)= A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(19)= A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(20)= AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))

C$$$      ANG2(21)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(22)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(23)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(24)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(25)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(26)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(27)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(28)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(29)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(30)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))

      ANG2(31)= ABP1M1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(32)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(33)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(34)= ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(35)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(36)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG2(37)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(38)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(39)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(40)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))

      ANG2(41)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(42)= A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(43)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(44)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(45)= A4P2M2(ANA(2,2),ANB(2,2),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(46)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(47)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(48)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(49)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
      ANG2(50)= A4M1P0(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))

      ANG2(51)= A4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(52)= A4P1P1(ANA(2,2),ANB(2,2),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(53)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(54)= A4P1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(55)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(56)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(57)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(59)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1))
      ANG2(60)= ABP1M2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))

      ANG2(61)= ABP1M2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(62)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(63)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(64)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(65)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(66)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
      ANG2(67)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(68)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(69)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(70)= A4M1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))

      ANG2(71)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(72)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))

      DEL = EPSS * MS**4
      IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
         ANG2(47)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
         ANG2(48)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
         ANG2(49)= C4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
         ANG2(51)= C4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
         ANG2(56)= C4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
         ANG2(57)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= CBP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1),DEL)
      ANG2(59)= CBP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1),DEL)
         ANG2(66)= C4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
         ANG2(69)= C4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
         ANG2(72)= C4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      END IF


      COLO2(9) = LOG(S4**2/MS2/(S4+MS2))

      M2GGH = 0D0
      M2GGH = M2GGH + N*CO*S4G2**(-1)*(S4+MS2)*TWO**(-1)*S4G * (  - 6*
     +    S**(-1)*TG*T1**(-1)*M2*S4**(-1) + 2*S**(-1)*TG*T1**(-1) + 8*
     +    S**(-1)*TG*M2*(S+TG)**(-2) + 6*S**(-1)*TG*M2*(S+TG)**(-1)*
     +    S4**(-1) - 24*S**(-1)*TG*M2*(S+UG)**(-2) + 8*S**(-1)*TG*M2*
     +    (S+UG)**(-1)*S4**(-1) - 4*S**(-1)*TG*M2**2*(S+UG)**(-3) + 24*
     +    S**(-1)*TG*M2**2*(S+UG)**(-2)*S4**(-1) + 4*S**(-1)*TG*M2**3*
     +    (S+UG)**(-3)*S4**(-1) - 10*S**(-1)*TG*(S+TG)**(-1) - 8*
     +    S**(-1)*TG*(S+UG)**(-1) + 4*S**(-1)*TG*S4**(-1) - 4*S**(-1)*
     +    TG**2*T1**(-1)*S4**(-1) - 4*S**(-1)*TG**2*M2*(S+TG)**(-3) - 8
     +    *S**(-1)*TG**2*M2*(S+TG)**(-2)*S4**(-1) - 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-3) + 22*S**(-1)*TG**2*M2*(S+UG)**(-2)*S4**(-1) + 4*
     +    S**(-1)*TG**2*M2**2*(S+UG)**(-3)*S4**(-1) + 4*S**(-1)*TG**2*
     +    (S+TG)**(-2) + 4*S**(-1)*TG**2*(S+TG)**(-1)*S4**(-1) - 22*
     +    S**(-1)*TG**2*(S+UG)**(-2) + 4*S**(-1)*TG**2*(S+UG)**(-1)*
     +    S4**(-1) - 4*S**(-1)*TG**3*M2*(S+UG)**(-3)*S4**(-1) + 8*
     +    S**(-1)*TG**3*(S+TG)**(-3) )
     +
      M2GGH = M2GGH + N*CO*S4G2**(-1)*(S4+MS2)*TWO**(-1)*S4G * (  - 4*
     +    S**(-1)*TG**3*(S+TG)**(-2)*S4**(-1) + 4*S**(-1)*TG**3*
     +    (S+UG)**(-3) + 4*S**(-1)*TG**3*(S+UG)**(-2)*S4**(-1) - 4*
     +    S**(-1)*TG**4*(S+TG)**(-3)*S4**(-1) - 4*S**(-1)*TG**4*
     +    (S+UG)**(-3)*S4**(-1) + 2*S**(-1)*T1**(-1)*M2 - 2*S**(-1)*
     +    T1**(-1)*M2**2*S4**(-1) - 4*S**(-1)*M2*(S+TG)**(-1) - 4*
     +    S**(-1)*M2*(S+UG)**(-1) + 2*S**(-1)*M2*S4**(-1) - 6*S**(-1)*
     +    M2**2*(S+UG)**(-2) + 4*S**(-1)*M2**2*(S+UG)**(-1)*S4**(-1) + 
     +    6*S**(-1)*M2**3*(S+UG)**(-2)*S4**(-1) - 2*S**(-1) - 4*S*
     +    TG**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 40*S*TG**(-1)*M2*
     +    (S+UG)**(-2) + 24*S*TG**(-1)*M2*(S+UG)**(-1)*S4**(-1) - 4*S*
     +    TG**(-1)*M2**2*(S+UG)**(-3) + 32*S*TG**(-1)*M2**2*
     +    (S+UG)**(-2)*S4**(-1) + 4*S*TG**(-1)*M2**3*(S+UG)**(-3)*
     +    S4**(-1) + 4*S*TG**(-1)*(S+TG)**(-1) - 16*S*TG**(-1)*
     +    (S+UG)**(-1) - 8*S*TG*M2*(S+UG)**(-3)*S4**(-1) - 8*S*TG*
     +    (S+TG)**(-3) )
     +
      M2GGH = M2GGH + N*CO*S4G2**(-1)*(S4+MS2)*TWO**(-1)*S4G * (  - 4*S
     +    *TG*(S+TG)**(-2)*S4**(-1) - 8*S*TG*(S+UG)**(-3) - 8*S*TG*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 4*S*TG**2*(S+TG)**(-3)*S4**(-1)
     +     - 8*S*TG**2*(S+UG)**(-3)*S4**(-1) - 8*S*UG**(-1)*M2*
     +    (S+TG)**(-2) + 8*S*UG**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 4*S*
     +    UG**(-1)*M2**2*(S+UG)**(-2)*S4**(-1) + 4*S*M2*(S+TG)**(-2)*
     +    S4**(-1) - 16*S*M2*(S+UG)**(-3) + 32*S*M2*(S+UG)**(-2)*
     +    S4**(-1) - 8*S*M2*(TG*UG-M2*S)**(-1)*S4**(-1) + 12*S*M2**2*
     +    (S+UG)**(-3)*S4**(-1) - 8*S*(S+TG)**(-1)*S4**(-1) - 40*S*
     +    (S+UG)**(-2) + 8*S*(S+UG)**(-1)*S4**(-1) + 4*S**2*TG*
     +    (S+TG)**(-3)*S4**(-1) + 4*S**2*UG**(-1)*M2**2*(S+UG)**(-3)*
     +    S4**(-1) - 6*TG**(-1)*M2*(S+UG)**(-1) - 6*TG**(-1)*M2**2*
     +    (S+UG)**(-2) + 6*TG**(-1)*M2**2*(S+UG)**(-1)*S4**(-1) + 6*
     +    TG**(-1)*M2**3*(S+UG)**(-2)*S4**(-1) + 4*TG*M2*(S+TG)**(-3)
     +     + 4*TG*M2*(S+TG)**(-2)*S4**(-1) - 12*TG*M2*(S+UG)**(-3) + 66
     +    *TG*M2*(S+UG)**(-2)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*S4G2**(-1)*(S4+MS2)*TWO**(-1)*S4G * (  - 8*
     +    TG*M2*(TG*UG-M2*S)**(-1)*S4**(-1) + 12*TG*M2**2*(S+UG)**(-3)*
     +    S4**(-1) + 4*TG*(S+TG)**(-2) - 12*TG*(S+TG)**(-1)*S4**(-1) - 
     +    66*TG*(S+UG)**(-2) + 20*TG*(S+UG)**(-1)*S4**(-1) + 8*TG*
     +    (TG*UG-M2*S)**(-1) - 8*TG**2*M2*(S+UG)**(-3)*S4**(-1) - 8*
     +    TG**2*(S+TG)**(-2)*S4**(-1) + 4*TG**2*(S+UG)**(-2)*S4**(-1)
     +     - 8*TG**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 4*TG**3*(S+TG)**(-3)
     +    *S4**(-1) - 12*TG**3*(S+UG)**(-3)*S4**(-1) - 8*UG**(-1)*M2*
     +    (S+TG)**(-1) - 8*T1**(-1)*M2*S4**(-1) + 8*T1**(-1) - 8*
     +    U1**(-1)*M2*S4**(-1) + 8*U1**(-1) - 12*M2*(S+TG)**(-2) - 16*
     +    M2*(S+TG)**(-1)*S4**(-1) - 52*M2*(S+UG)**(-2) + 42*M2*
     +    (S+UG)**(-1)*S4**(-1) - 8*M2**2*(S+UG)**(-3) + 52*M2**2*
     +    (S+UG)**(-2)*S4**(-1) + 8*M2**3*(S+UG)**(-3)*S4**(-1) + 20*
     +    (S+TG)**(-1) - 34*(S+UG)**(-1) + 4*S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4 + 40*S**(-2)*
     +    TG*M2 + 12*S**(-1)*TG**(-1)*M2*MS2 + 12*S**(-1)*TG**(-1)*
     +    M2**2 + 20*S**(-1)*UG**(-1)*M2*MS2 + 20*S**(-1)*UG**(-1)*
     +    M2**2 - 56*S**(-1)*M2 - 48*S**(-1)*MS2 - 32*TG**(-2)*M2*MS2
     +     - 32*TG**(-2)*M2**2 - 48*TG**(-1)*UG**(-1)*M2*MS2 - 48*
     +    TG**(-1)*UG**(-1)*M2**2 - 76*TG**(-1)*M2 - 84*TG**(-1)*MS2 - 
     +    64*UG**(-2)*M2*MS2 - 64*UG**(-2)*M2**2 - 90*UG**(-1)*M2 - 84*
     +    UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +    *T1**(-1)*M2 + 4*S**(-1)*TG**(-1)*T1**(-1)*M2**2*S4**(-1) + 8
     +    *S**(-1)*TG**(-1)*M2*(S+UG)**(-1) - 4*S**(-1)*TG**(-1)*M2**2*
     +    (S+UG)**(-1)*S4**(-1) - 4*S**(-1)*TG**(-1)*(S+UG)**(-1)*S4 - 
     +    8*S**(-1)*TG*T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 8*
     +    S**(-1)*TG*T1**(-2)*M2*S4**(-1) + 8*S**(-1)*TG*T1**(-2) + 8*
     +    S**(-1)*TG*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 + 10*S**(-1)
     +    *TG*T1**(-1)*S4**(-1) - 16*S**(-1)*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2) + 20*S**(-1)*TG*M2*(S+TG)**(-2)*S4**(-1)
     +     + 8*S**(-1)*TG*M2*(S+UG)**(-3) - 4*S**(-1)*TG*M2*
     +    (S+UG)**(-2)*S4**(-1) + 16*S**(-1)*TG*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 4*S**(-1)*TG*M2**2*(S+TG)**(-3)
     +    *S4**(-1) - 2*S**(-1)*TG*(S+TG)**(-1)*S4**(-1) - 4*S**(-1)*TG
     +    *(S+UG)**(-3)*S4 + 28*S**(-1)*TG*(S+UG)**(-2) + 12*S**(-1)*TG
     +    *(S+UG)**(-1)*S4**(-1) + 16*S**(-1)*TG**2*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)*TG**2*M2
     +    *(S+TG)**(-3)*S4**(-1) - 8*S**(-1)*TG**2*(S+TG)**(-3) - 12*
     +    S**(-1)*TG**2*(S+TG)**(-2)*S4**(-1) + 4*S**(-1)*TG**2*
     +    (S+UG)**(-3) + 2*S**(-1)*TG**2*(S+UG)**(-2)*S4**(-1) + 4*
     +    S**(-1)*TG**3*(S+TG)**(-3)*S4**(-1) - 4*S**(-1)*UG**(-1)*
     +    U1**(-1)*M2 + 4*S**(-1)*UG**(-1)*U1**(-1)*M2**2*S4**(-1) - 12
     +    *S**(-1)*UG**(-1)*M2*(S+TG)**(-2)*S4 + 20*S**(-1)*UG**(-1)*M2
     +    *(S+TG)**(-1) - 4*S**(-1)*UG**(-1)*M2**2*(S+TG)**(-3)*S4 + 32
     +    *S**(-1)*UG**(-1)*M2**2*(S+TG)**(-2) - 20*S**(-1)*UG**(-1)*
     +    M2**2*(S+TG)**(-1)*S4**(-1) + 8*S**(-1)*UG**(-1)*M2**3*
     +    (S+TG)**(-3) - 20*S**(-1)*UG**(-1)*M2**3*(S+TG)**(-2)*
     +    S4**(-1) - 4*S**(-1)*UG**(-1)*M2**4*(S+TG)**(-3)*S4**(-1) - 4
     +    *S**(-1)*UG**(-1)*(S+TG)**(-1)*S4 + 8*S**(-1)*T1**(-2)*M2*MS2
     +    *(TG*UG-M2*S)**(-1)*S4 + 16*S**(-1)*T1**(-2)*M2 + 8*S**(-1)*
     +    T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*T1**(-2)*
     +    M2**2*S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 22*S**(-1)*T1**(-1)*
     +    M2*S4**(-1) - 4*S**(-1)*T1**(-1) + 8*S**(-1)*U1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 + 8*S**(-1)*U1**(-2)*M2 + 8*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 + 8*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*U1**(-1)*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**2 + 4*S**(-1)*U1**(-1)*M2*S4**(-1)
     +     - 8*S**(-1)*U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 + 4*
     +    S**(-1)*U1**(-1) + 8*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 + 
     +    16*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) + 12*S**(-1)*M2
     +    *(S+TG)**(-2) - 12*S**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 24*
     +    S**(-1)*M2*(S+UG)**(-2) + 8*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     + 4*S**(-1)*M2**2*(S+TG)**(-3) - 20*S**(-1)*M2**2*
     +    (S+TG)**(-2)*S4**(-1) - 6*S**(-1)*M2**2*(S+UG)**(-2)*S4**(-1)
     +     + 16*S**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 4*S**(-1)*
     +    M2**3*(S+TG)**(-3)*S4**(-1) + 20*S**(-1)*(S+TG)**(-1) - 12*
     +    S**(-1)*(S+UG)**(-2)*S4 )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 12*S**(-1)*
     +    (S+UG)**(-1) + 10*S**(-1)*S4**(-1) + 48*S*TG**(-1)*M2*
     +    (S+UG)**(-2)*S4**(-1) + 4*S*TG**(-1)*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) + 28*S*TG**(-1)*M2**2*(S+UG)**(-3)*S4**(-1) + 8*S*
     +    TG**(-1)*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 4*S*TG**(-1)*
     +    (S+TG)**(-1)*S4**(-1) + 16*S*TG**(-1)*(S+UG)**(-1)*S4**(-1)
     +     - 64*S*TG*UG**(-1)*(S+TG)**(-3) - 32*S*TG*UG**(-1)*
     +    (S+TG)**(-2)*S4**(-1) + 16*S*TG*MS2*(TG*UG-M2*S)**(-2)*
     +    S4**(-1) - 4*S*TG*(S+TG)**(-3)*S4**(-1) + 32*S*TG*
     +    (S+UG)**(-3)*S4**(-1) + 32*S*TG**2*UG**(-1)*(S+TG)**(-3)*
     +    S4**(-1) + 32*S*UG**(-2)*M2*MS2*(S+TG)**(-2)*S4**(-1) + 32*S*
     +    UG**(-2)*M2**2*(S+TG)**(-2)*S4**(-1) - 32*S*UG**(-1)*M2*MS2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 16*S*UG**(-1)*M2*
     +    (S+TG)**(-2)*S4**(-1) + 8*S*UG**(-1)*M2*(S+UG)**(-2)*S4**(-1)
     +     + 4*S*UG**(-1)*M2*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*S*UG**(-1)
     +    *M2**2*(S+TG)**(-3)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S*UG**(-1)*M2**2
     +    *(S+UG)**(-3)*S4**(-1) + 8*S*UG**(-1)*MS2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) + 32*S*UG**(-1)*(S+TG)**(-3)*S4 + 32*S*UG**(-1)*
     +    (S+TG)**(-2) + 16*S*M2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 16*S
     +    *M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 64*S*M2*
     +    (S+UG)**(-3)*S4**(-1) - 8*S*MS2*(TG*UG-M2*S)**(-2) + 40*S*
     +    (S+TG)**(-2)*S4**(-1) + 48*S*(S+UG)**(-2)*S4**(-1) + 4*S*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 64*S**2*TG*UG**(-1)*
     +    (S+TG)**(-3)*S4**(-1) - 64*S**2*UG**(-1)*(S+TG)**(-3) - 32*
     +    S**2*UG**(-1)*(S+TG)**(-2)*S4**(-1) + 32*S**3*UG**(-1)*
     +    (S+TG)**(-3)*S4**(-1) - 16*TG**(-2)*T1**(-1)*M2*MS2*S4**(-1)
     +     - 16*TG**(-2)*T1**(-1)*M2**2*S4**(-1) - 48*TG**(-2)*M2*MS2*
     +    (S+UG)**(-1)*S4**(-1) + 16*TG**(-2)*M2*(S+UG)**(-1) - 32*
     +    TG**(-2)*M2**2*MS2*(S+UG)**(-2)*S4**(-1) - 48*TG**(-2)*M2**2*
     +    (S+UG)**(-1)*S4**(-1) - 32*TG**(-2)*M2**3*(S+UG)**(-2)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 16*TG**(-2)*MS2*
     +    (S+UG)**(-1) + 4*TG**(-1)*UG**(-1)*M2**2*(S+UG)**(-3)*S4 - 4*
     +    TG**(-1)*UG**(-1)*M2**2*(S+UG)**(-2) - 8*TG**(-1)*UG**(-1)*
     +    M2**3*(S+UG)**(-3) + 4*TG**(-1)*UG**(-1)*M2**3*(S+UG)**(-2)*
     +    S4**(-1) + 4*TG**(-1)*UG**(-1)*M2**4*(S+UG)**(-3)*S4**(-1) + 
     +    16*TG**(-1)*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) - 28*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) - 8*TG**(-1)*T1**(-1)*MS2*S4**(-1) - 32*
     +    TG**(-1)*M2*MS2*(S+UG)**(-2)*S4**(-1) + 32*TG**(-1)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 12*TG**(-1)*M2*(S+UG)**(-2) - 
     +    70*TG**(-1)*M2*(S+UG)**(-1)*S4**(-1) - 8*TG**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) + 32*TG**(-1)*M2**2*MS2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 4*TG**(-1)*M2**2*(S+UG)**(-3)
     +     - 66*TG**(-1)*M2**2*(S+UG)**(-2)*S4**(-1) + 8*TG**(-1)*M2**2
     +    *(TG*UG-M2*S)**(-1)*S4**(-1) + 4*TG**(-1)*M2**3*(S+UG)**(-3)*
     +    S4**(-1) - 16*TG**(-1)*MS2*(S+UG)**(-1)*S4**(-1) - 16*
     +    TG**(-1)*MS2*(TG*UG-M2*S)**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 16*TG**(-1)*
     +    (S+UG)**(-1) + 32*TG*UG**(-2)*M2*MS2*(S+TG)**(-2)*S4**(-1) + 
     +    32*TG*UG**(-2)*M2**2*(S+TG)**(-2)*S4**(-1) - 32*TG*UG**(-1)*
     +    M2*MS2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 32*TG*
     +    UG**(-1)*M2*(S+TG)**(-2)*S4**(-1) + 32*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 4*TG*M2*(S+TG)**(-3)*S4**(-1)
     +     - 16*TG*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 16*TG*
     +    M2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*MS2*
     +    (TG*UG-M2*S)**(-2) + 8*TG*(S+TG)**(-3) + 52*TG*(S+TG)**(-2)*
     +    S4**(-1) + 12*TG*(S+UG)**(-3) - 14*TG*(S+UG)**(-2)*S4**(-1)
     +     + 16*TG**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 16*UG**(-2)*
     +    U1**(-1)*M2*MS2*S4**(-1) - 16*UG**(-2)*U1**(-1)*M2**2*
     +    S4**(-1) - 32*UG**(-2)*M2*MS2*(S+TG)**(-2) - 48*UG**(-2)*M2*
     +    MS2*(S+TG)**(-1)*S4**(-1) + 16*UG**(-2)*M2*(S+TG)**(-1) - 32*
     +    UG**(-2)*M2**2*(S+TG)**(-2) - 48*UG**(-2)*M2**2*(S+TG)**(-1)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 16*UG**(-2)*MS2*
     +    (S+TG)**(-1) + 16*UG**(-1)*U1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)
     +     - 28*UG**(-1)*U1**(-1)*M2*S4**(-1) - 8*UG**(-1)*U1**(-1)*MS2
     +    *S4**(-1) + 32*UG**(-1)*M2*MS2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1) + 32*UG**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) + 4*UG**(-1)*M2*(S+TG)**(-2) - 92*UG**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 8*UG**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 8*UG**(-1)*M2*(TG*UG-M2*S)**(-1) + 12*UG**(-1)*M2**2*
     +    (S+TG)**(-3) - 40*UG**(-1)*M2**2*(S+TG)**(-2)*S4**(-1) - 4*
     +    UG**(-1)*M2**2*(S+UG)**(-3) + 4*UG**(-1)*M2**2*(S+UG)**(-2)*
     +    S4**(-1) + 8*UG**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 12*
     +    UG**(-1)*M2**3*(S+TG)**(-3)*S4**(-1) + 4*UG**(-1)*M2**3*
     +    (S+UG)**(-3)*S4**(-1) - 16*UG**(-1)*MS2*(S+TG)**(-1)*S4**(-1)
     +     - 16*UG**(-1)*MS2*(TG*UG-M2*S)**(-1) + 16*UG**(-1)*
     +    (S+TG)**(-1) + 8*T1**(-1)*M2*(TG*UG-M2*S)**(-1) - 20*T1**(-1)
     +    *S4**(-1) )
     +
      M2GGH = M2GGH + N*CO*(S4+MS2)*TWO**(-1) * ( 8*U1**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) - 20*U1**(-1)*S4**(-1) + 32*M2*MS2*
     +    (S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 24*M2*MS2*
     +    (TG*UG-M2*S)**(-2) - 20*M2*(S+TG)**(-2)*S4**(-1) + 16*M2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 8*M2*(S+UG)**(-3) - 52*M2*
     +    (S+UG)**(-2)*S4**(-1) + 36*M2*(TG*UG-M2*S)**(-1)*S4**(-1) + 
     +    16*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 8*M2**2*
     +    (S+TG)**(-3)*S4**(-1) + 16*M2**2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 8*MS2*(TG*UG-M2*S)**(-2)*S4 + 
     +    32*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 40*(S+TG)**(-2) - 84*
     +    (S+TG)**(-1)*S4**(-1) + 36*(S+UG)**(-2) - 38*(S+UG)**(-1)*
     +    S4**(-1) - 4*(TG*UG-M2*S)**(-1) )
     +
      M2GGH = M2GGH + N*CO*TWO**(-1) * (  - 16*S**(-2)*M2 + 8*S**(-1)*
     +    TG**(-1)*M2 + 12*S**(-1)*TG**(-1)*MS2 + 4*S**(-1)*TG**(-1)*S4
     +     + 10*S**(-1)*UG**(-1)*M2 + 12*S**(-1)*UG**(-1)*MS2 + 8*
     +    S**(-1)*UG**(-1)*S4 - 36*S**(-1) + 48*TG**(-1)*UG**(-1)*MS2
     +     + 24*TG**(-1)*UG**(-1)*S4 - 34*TG**(-1) - 38*UG**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * ( 12*S**(-1)*TG**(-1)*
     +    T1**(-1)*M2 - 12*S**(-1)*TG**(-1)*T1**(-1)*M2**2*S4**(-1) - 
     +    12*S**(-1)*TG**(-1)*M2*S4**(-1) + 12*S**(-1)*TG**(-1) + 24*
     +    S**(-1)*TG*T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 + 24*
     +    S**(-1)*TG*T1**(-2)*M2*S4**(-1) - 24*S**(-1)*TG*T1**(-2) - 24
     +    *S**(-1)*TG*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 - 4*S**(-1)
     +    *TG*T1**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 6*S**(-1)*TG*T1**(-1)
     +    *M2*(TG+UG)**(-1)*S4**(-1) + 8*S**(-1)*TG*T1**(-1)*
     +    (S+UG)**(-1) - 2*S**(-1)*TG*T1**(-1)*(TG+UG)**(-1) - 38*
     +    S**(-1)*TG*T1**(-1)*S4**(-1) + 48*S**(-1)*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2) + 8*S**(-1)*TG*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 8*S**(-1)*TG*M2*SYMBU*S4**(-1) - 48*S**(-1)*TG*
     +    M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 40*S**(-1)*TG*
     +    (S+TG)**(-1)*S4**(-1) - 36*S**(-1)*TG*(S+UG)**(-1)*S4**(-1)
     +     - 4*S**(-1)*TG*(TG+UG)**(-1)*S4**(-1) - 8*S**(-1)*TG*
     +    (TG*UG-M2*S)**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * ( 8*S**(-1)*TG*SYMBU - 
     +    4*S**(-1)*TG**2*T1**(-1)*(S+UG)**(-1)*S4**(-1) + 4*S**(-1)*
     +    TG**2*T1**(-1)*(TG+UG)**(-1)*S4**(-1) - 48*S**(-1)*TG**2*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 8*S**(-1)*TG**2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 4*S**(-1)*TG**2*SYMBT*S4**(-1)
     +     - 4*S**(-1)*TG**2*SYMBU*S4**(-1) - 4*S**(-1)*UG**(-1)*
     +    T1**(-1)*M2 + 4*S**(-1)*UG**(-1)*U1**(-1)*M2 - 4*S**(-1)*
     +    UG**(-1)*U1**(-1)*M2**2*S4**(-1) - 16*S**(-1)*UG**(-1)*M2*
     +    S4**(-1) + 20*S**(-1)*UG**(-1) - 24*S**(-1)*T1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 - 48*S**(-1)*T1**(-2)*M2 - 24*S**(-1)*
     +    T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 + 24*S**(-1)*T1**(-2)*
     +    M2**2*S4**(-1) - 2*S**(-1)*T1**(-1)*M2*(TG+UG)**(-1) - 62*
     +    S**(-1)*T1**(-1)*M2*S4**(-1) + 2*S**(-1)*T1**(-1)*M2**2*
     +    (TG+UG)**(-1)*S4**(-1) + 12*S**(-1)*T1**(-1) - 24*S**(-1)*
     +    U1**(-2)*M2*MS2*(TG*UG-M2*S)**(-1)*S4 - 24*S**(-1)*U1**(-2)*
     +    M2 )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 24*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 24*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 + 24*S**(-1)*U1**(-1)*M2
     +    *MS2*(TG*UG-M2*S)**(-2)*S4**2 - 8*S**(-1)*U1**(-1)*M2*
     +    (S+TG)**(-1) + 4*S**(-1)*U1**(-1)*M2*S4**(-1) + 24*S**(-1)*
     +    U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 - 20*S**(-1)*
     +    U1**(-1) - 24*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 - 48*
     +    S**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 12*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 52*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 2*S**(-1)*M2*(TG+UG)**(-1)*S4**(-1) - 8*S**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) + 8*S**(-1)*M2*SYMBU - 44*S**(-1)*M2**2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 4*S**(-1)*M2**2*SYMBU*S4**(-1)
     +     - 28*S**(-1)*(S+TG)**(-1) + 4*S**(-1)*(S+UG)**(-1) + 2*
     +    S**(-1)*(TG+UG)**(-1) + 4*S**(-1)*(TG*UG-M2*S)**(-1)*S4 - 4*
     +    S**(-1)*SYMBU*S4 - 38*S**(-1)*S4**(-1) + 24*S*TG**(-1)*M2*
     +    (S+UG)**(-3) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 56*S*TG**(-1)*M2*
     +    (S+UG)**(-2)*S4**(-1) - 4*S*TG**(-1)*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 40*S*TG**(-1)*M2**2*(S+UG)**(-3)*S4**(-1) + 8*S*
     +    TG**(-1)*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) + 4*S*TG**(-1)*
     +    (S+TG)**(-1)*S4**(-1) + 16*S*TG**(-1)*(S+UG)**(-2) - 16*S*
     +    TG**(-1)*(S+UG)**(-1)*S4**(-1) + 4*S*TG**(-1)*SYMBU + 64*S*TG
     +    *UG**(-1)*(S+TG)**(-3) + 32*S*TG*UG**(-1)*(S+TG)**(-2)*
     +    S4**(-1) - 48*S*TG*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 8*S*TG*
     +    (S+TG)**(-3)*S4**(-1) - 40*S*TG*(S+UG)**(-3)*S4**(-1) - 10*S*
     +    TG*(TG+UG)**(-1)*SYMBT*S4**(-1) + 2*S*TG*(TG+UG)**(-1)*SYMBU*
     +    S4**(-1) - 32*S*TG**2*UG**(-1)*(S+TG)**(-3)*S4**(-1) - 32*S*
     +    UG**(-2)*M2*MS2*(S+TG)**(-2)*S4**(-1) - 32*S*UG**(-2)*M2**2*
     +    (S+TG)**(-2)*S4**(-1) + 32*S*UG**(-1)*M2*MS2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 24*S*UG**(-1)*M2*(S+TG)**(-3)
     +     - 56*S*UG**(-1)*M2*(S+TG)**(-2)*S4**(-1) - 4*S*UG**(-1)*M2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 8*S*UG**(-1)*M2**2
     +    *(S+TG)**(-3)*S4**(-1) + 8*S*UG**(-1)*MS2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 32*S*UG**(-1)*(S+TG)**(-3)*S4 - 16*S*UG**(-1)*
     +    (S+TG)**(-2) - 16*S*UG**(-1)*(S+TG)**(-1)*S4**(-1) + 4*S*
     +    UG**(-1)*(S+UG)**(-1)*S4**(-1) + 4*S*UG**(-1)*SYMBT - 48*S*M2
     +    *MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 48*S*M2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 80*S*M2*(S+UG)**(-3)*S4**(-1)
     +     - 2*S*M2*(TG+UG)**(-1)*SYMBT*S4**(-1) + 2*S*M2*(TG+UG)**(-1)
     +    *SYMBU*S4**(-1) + 24*S*MS2*(TG*UG-M2*S)**(-2) - 40*S*
     +    (S+TG)**(-2)*S4**(-1) + 24*S*(S+UG)**(-3) - 56*S*(S+UG)**(-2)
     +    *S4**(-1) + 4*S*(S+UG)**(-1)*(TG*UG-M2*S)**(-1) + 8*S*
     +    (TG+UG)**(-1)*SYMBT - 8*S*(TG*UG-M2*S)**(-1)*S4**(-1) - 6*S*
     +    SYMBT*S4**(-1) - 64*S**2*TG*UG**(-1)*(S+TG)**(-3)*S4**(-1) + 
     +    64*S**2*UG**(-1)*(S+TG)**(-3) + 32*S**2*UG**(-1)*(S+TG)**(-2)
     +    *S4**(-1) - 6*S**2*(TG+UG)**(-1)*SYMBT*S4**(-1) - 32*S**3*
     +    UG**(-1)*(S+TG)**(-3)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * ( 16*TG**(-2)*T1**(-1)*
     +    M2*MS2*S4**(-1) + 16*TG**(-2)*T1**(-1)*M2**2*S4**(-1) + 48*
     +    TG**(-2)*M2*MS2*(S+UG)**(-1)*S4**(-1) - 16*TG**(-2)*M2*
     +    (S+UG)**(-1) + 32*TG**(-2)*M2**2*MS2*(S+UG)**(-2)*S4**(-1) + 
     +    48*TG**(-2)*M2**2*(S+UG)**(-1)*S4**(-1) + 32*TG**(-2)*M2**3*
     +    (S+UG)**(-2)*S4**(-1) - 16*TG**(-2)*MS2*(S+UG)**(-1) - 16*
     +    TG**(-1)*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) + 28*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) + 24*TG**(-1)*T1**(-1)*MS2*S4**(-1) + 32
     +    *TG**(-1)*M2*MS2*(S+UG)**(-2)*S4**(-1) - 32*TG**(-1)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 8*TG**(-1)*M2*(S+TG)**(-1)*
     +    S4**(-1) + 68*TG**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 16*TG**(-1)
     +    *M2*(TG*UG-M2*S)**(-1) + 8*TG**(-1)*M2*SYMBU - 32*TG**(-1)*
     +    M2**2*MS2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 64*
     +    TG**(-1)*M2**2*(S+UG)**(-2)*S4**(-1) - 16*TG**(-1)*M2**2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 4*TG**(-1)*M2**2*SYMBU*S4**(-1)
     +     + 16*TG**(-1)*MS2*(S+UG)**(-1)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * ( 16*TG**(-1)*MS2*
     +    (TG*UG-M2*S)**(-1) - 8*TG**(-1)*(S+TG)**(-1) - 20*TG**(-1)*
     +    (S+UG)**(-1) - 4*TG**(-1)*SYMBU*S4 - 4*TG**(-1)*S4**(-1) - 32
     +    *TG*UG**(-2)*M2*MS2*(S+TG)**(-2)*S4**(-1) - 32*TG*UG**(-2)*
     +    M2**2*(S+TG)**(-2)*S4**(-1) + 32*TG*UG**(-1)*M2*MS2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 32*TG*UG**(-1)*M2*
     +    (S+TG)**(-2)*S4**(-1) + 4*TG*T1**(-1)*(S+UG)**(-1)*S4**(-1)
     +     - 2*TG*T1**(-1)*(TG+UG)**(-1)*S4**(-1) - 96*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 8*TG*M2*(S+TG)**(-3)*S4**(-1)
     +     + 48*TG*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*
     +    M2*(S+UG)**(-3)*S4**(-1) - 48*TG*M2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 2*TG*M2*(TG+UG)**(-1)*SYMBT*
     +    S4**(-1) - 6*TG*M2*(TG+UG)**(-1)*SYMBU*S4**(-1) + 48*TG*MS2*
     +    (TG*UG-M2*S)**(-2) - 16*TG*(S+TG)**(-3) - 16*TG*(S+TG)**(-2)*
     +    S4**(-1) - 4*TG*(S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 32*TG*
     +    (S+UG)**(-3) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * ( 16*TG*(S+UG)**(-2)*
     +    S4**(-1) + 4*TG*(S+UG)**(-1)*(TG*UG-M2*S)**(-1) + 6*TG*
     +    (TG+UG)**(-1)*SYMBT + 2*TG*(TG+UG)**(-1)*SYMBU + 8*TG*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 10*TG*SYMBT*S4**(-1) + 2*TG*
     +    SYMBU*S4**(-1) - 48*TG**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 8
     +    *TG**2*(S+TG)**(-3)*S4**(-1) - 8*TG**2*(S+UG)**(-3)*S4**(-1)
     +     - 4*TG**2*(TG+UG)**(-1)*SYMBT*S4**(-1) - 4*TG**2*
     +    (TG+UG)**(-1)*SYMBU*S4**(-1) + 16*UG**(-2)*U1**(-1)*M2*MS2*
     +    S4**(-1) + 16*UG**(-2)*U1**(-1)*M2**2*S4**(-1) + 32*UG**(-2)*
     +    M2*MS2*(S+TG)**(-2) + 48*UG**(-2)*M2*MS2*(S+TG)**(-1)*
     +    S4**(-1) - 16*UG**(-2)*M2*(S+TG)**(-1) + 32*UG**(-2)*M2**2*
     +    (S+TG)**(-2) + 48*UG**(-2)*M2**2*(S+TG)**(-1)*S4**(-1) - 16*
     +    UG**(-2)*MS2*(S+TG)**(-1) + 4*UG**(-1)*T1**(-1)*M2*S4**(-1)
     +     - 16*UG**(-1)*U1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) + 44*
     +    UG**(-1)*U1**(-1)*M2*S4**(-1) + 24*UG**(-1)*U1**(-1)*MS2*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 32*UG**(-1)*M2*MS2
     +    *(S+TG)**(-1)*(TG*UG-M2*S)**(-1) - 32*UG**(-1)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 24*UG**(-1)*M2*(S+TG)**(-3)*S4
     +     + 72*UG**(-1)*M2*(S+TG)**(-2) + 52*UG**(-1)*M2*(S+TG)**(-1)*
     +    S4**(-1) + 8*UG**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 16*UG**(-1)*
     +    M2*(TG*UG-M2*S)**(-1) + 8*UG**(-1)*M2*SYMBT + 32*UG**(-1)*
     +    M2**2*(S+TG)**(-3) - 24*UG**(-1)*M2**2*(S+TG)**(-2)*S4**(-1)
     +     - 16*UG**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 4*UG**(-1)
     +    *M2**2*SYMBT*S4**(-1) - 8*UG**(-1)*M2**3*(S+TG)**(-3)*
     +    S4**(-1) + 16*UG**(-1)*MS2*(S+TG)**(-1)*S4**(-1) + 16*
     +    UG**(-1)*MS2*(TG*UG-M2*S)**(-1) - 20*UG**(-1)*(S+TG)**(-1) - 
     +    8*UG**(-1)*(S+UG)**(-1) - 4*UG**(-1)*SYMBT*S4 - 20*UG**(-1)*
     +    S4**(-1) + 4*T1**(-1)*M2*(S+UG)**(-1)*S4**(-1) - 2*T1**(-1)*
     +    M2*(TG+UG)**(-1)*S4**(-1) - 24*T1**(-1)*M2*(TG*UG-M2*S)**(-1)
     +     - 8*T1**(-1)*(S+UG)**(-1) + 20*T1**(-1)*S4**(-1) - 24*
     +    U1**(-1)*M2*(TG*UG-M2*S)**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 8*U1**(-1)*
     +    (S+TG)**(-1) + 36*U1**(-1)*S4**(-1) - 32*M2*MS2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 72*M2*MS2*(TG*UG-M2*S)**(-2) + 
     +    24*M2*(S+TG)**(-3) - 24*M2*(S+TG)**(-2)*S4**(-1) - 52*M2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 32*M2*(S+UG)**(-3) + 48*M2*
     +    (S+UG)**(-2)*S4**(-1) + 2*M2*(TG+UG)**(-1)*SYMBT + 2*M2*
     +    (TG+UG)**(-1)*SYMBU - 76*M2*(TG*UG-M2*S)**(-1)*S4**(-1) - 4*
     +    M2*SYMBT*S4**(-1) - 2*M2*SYMBU*S4**(-1) - 48*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) - 8*M2**2*(S+TG)**(-3)*S4**(-1)
     +     - 8*M2**2*(S+UG)**(-3)*S4**(-1) - 48*M2**2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 2*M2**2*(TG+UG)**(-1)*SYMBU*
     +    S4**(-1) - 24*MS2*(TG*UG-M2*S)**(-2)*S4 - 32*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 32*(S+TG)**(-2) + 4*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4 + 100*(S+TG)**(-1)*
     +    S4**(-1) - 24*(S+UG)**(-3)*S4 + 16*(S+UG)**(-2) + 56*
     +    (S+UG)**(-1)*S4**(-1) )
     +
      M2GGH = M2GGH + N*CK*(S4+MS2)*TWO**(-1) * (  - 2*(TG+UG)**(-1)*
     +    SYMBT*S4 - 6*(TG+UG)**(-1)*S4**(-1) - 4*(TG*UG-M2*S)**(-1) + 
     +    6*SYMBT + 4*SYMBU )
     +
      M2GGH = M2GGH + N*CK*TWO**(-1) * ( 8*S**(-2)*TG*(TG+UG)**(-1)*S4
     +     - 8*S**(-2)*TG - 16*S**(-2)*M2*(TG+UG)**(-1)*S4 + 8*S**(-2)*
     +    M2 - 8*S**(-2)*M2**2*(TG+UG)**(-1) - 8*S**(-2)*(TG+UG)**(-1)*
     +    S4**2 + 8*S**(-2)*S4 + 6*S**(-1)*M2*(TG+UG)**(-1) - 16*
     +    S**(-1)*MS2*(TG+UG)**(-1) - 6*S**(-1)*(TG+UG)**(-1)*S4 + 4*
     +    S**(-1) + 14*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + CQED*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +    *T1**(-1)*M2 + 4*S**(-1)*TG**(-1)*T1**(-1)*M2**2*S4**(-1) - 8
     +    *S**(-1)*TG*T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 8*
     +    S**(-1)*TG*T1**(-2)*M2*S4**(-1) + 8*S**(-1)*TG*T1**(-2) + 8*
     +    S**(-1)*TG*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 + 4*S**(-1)*
     +    TG*T1**(-1)*M2*(S+UG)**(-1)*S4**(-1) - 8*S**(-1)*TG*T1**(-1)*
     +    (S+UG)**(-1) + 16*S**(-1)*TG*T1**(-1)*S4**(-1) - 16*S**(-1)*
     +    TG*M2*MS2*(TG*UG-M2*S)**(-2) + 16*S**(-1)*TG*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) - 12*S**(-1)*TG*(S+TG)**(-1)*
     +    S4**(-1) + 8*S**(-1)*TG*(S+UG)**(-1)*S4**(-1) + 4*S**(-1)*
     +    TG**2*T1**(-1)*(S+UG)**(-1)*S4**(-1) + 16*S**(-1)*TG**2*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 4*S**(-1)*UG**(-1)*M2*
     +    S4**(-1) - 4*S**(-1)*UG**(-1) + 8*S**(-1)*T1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 + 16*S**(-1)*T1**(-2)*M2 + 8*S**(-1)*
     +    T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*T1**(-2)*
     +    M2**2*S4**(-1) )
     +
      M2GGH = M2GGH + CQED*(S4+MS2)*TWO**(-1) * ( 20*S**(-1)*T1**(-1)*
     +    M2*S4**(-1) - 4*S**(-1)*T1**(-1) + 8*S**(-1)*U1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 + 8*S**(-1)*U1**(-2)*M2 + 8*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 + 8*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*U1**(-1)*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**2 + 8*S**(-1)*U1**(-1)*M2*
     +    (S+TG)**(-1) - 8*S**(-1)*U1**(-1)*M2*S4**(-1) - 8*S**(-1)*
     +    U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 + 8*S**(-1)*U1**(-1)
     +     + 8*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 + 16*S**(-1)*M2*MS2
     +    *(TG*UG-M2*S)**(-1)*S4**(-1) + 4*S**(-1)*M2*(S+TG)**(-1)*
     +    S4**(-1) + 16*S**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 16*S**(-1)*
     +    M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) + 4*S**(-1)*(S+TG)**(-1) + 
     +    8*S**(-1)*S4**(-1) - 8*S*TG**(-1)*MS2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) + 16*S*TG*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 4*S*TG*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 4*S*TG*(TG+UG)**(-1)*SYMBU*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + CQED*(S4+MS2)*TWO**(-1) * (  - 8*S*UG**(-1)*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 16*S*M2*MS2*(TG*UG-M2*S)**(-2)*
     +    S4**(-1) - 16*S*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1)
     +     - 4*S*M2*(TG+UG)**(-1)*SYMBU*S4**(-1) - 8*S*MS2*
     +    (TG*UG-M2*S)**(-2) + 4*S*(S+TG)**(-1)*(TG*UG-M2*S)**(-1) - 4*
     +    S*(TG+UG)**(-1)*SYMBT + 4*S*SYMBT*S4**(-1) + 4*S**2*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 8*TG**(-1)*T1**(-1)*MS2*
     +    S4**(-1) - 4*TG**(-1)*M2*(TG*UG-M2*S)**(-1) + 4*TG**(-1)*
     +    M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) + 4*TG*T1**(-1)*
     +    (TG+UG)**(-1)*S4**(-1) + 32*TG*M2*MS2*(TG*UG-M2*S)**(-2)*
     +    S4**(-1) - 16*TG*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1)
     +     + 16*TG*M2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*
     +    MS2*(TG*UG-M2*S)**(-2) + 4*TG*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)
     +     - 4*TG*(S+UG)**(-1)*(TG*UG-M2*S)**(-1) + 4*TG*SYMBT*S4**(-1)
     +     - 4*TG*SYMBU*S4**(-1) + 16*TG**2*MS2*(TG*UG-M2*S)**(-2)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + CQED*(S4+MS2)*TWO**(-1) * (  - 8*UG**(-1)*
     +    U1**(-1)*M2*S4**(-1) - 8*UG**(-1)*U1**(-1)*MS2*S4**(-1) - 4*
     +    UG**(-1)*M2*(TG*UG-M2*S)**(-1) + 4*UG**(-1)*M2**2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 8*UG**(-1)*S4**(-1) + 4*
     +    T1**(-1)*M2*(TG+UG)**(-1)*S4**(-1) + 8*T1**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) + 8*U1**(-1)*M2*(TG*UG-M2*S)**(-1) - 8*
     +    U1**(-1)*S4**(-1) - 24*M2*MS2*(TG*UG-M2*S)**(-2) + 20*M2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 24*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 4*M2*SYMBU*S4**(-1) + 16*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 16*M2**2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 8*MS2*(TG*UG-M2*S)**(-2)*S4 - 4
     +    *(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4 - 12*(S+TG)**(-1)*
     +    S4**(-1) + 4*(TG+UG)**(-1)*S4**(-1) + 4*SYMBU )
     +
      M2GGH = M2GGH + ANG2(13)*N*CO*TWO**(-1) * ( 16*S**(-2)*TG*M2*S4
     +     - 8*S**(-2)*TG**2*M2 - 8*S**(-2)*M2*S4**2 - 8*S**(-1)*TG*M2
     +     + 8*S**(-1)*M2*S4 - 8*M2 )
     +
      M2GGH = M2GGH + ANG2(13)*N*CK*TWO**(-1) * (  - 16*S**(-2)*TG*M2*
     +    S4 + 8*S**(-2)*TG**2*M2 + 8*S**(-2)*M2*S4**2 + 8*S**(-1)*TG*
     +    M2 - 8*S**(-1)*M2*S4 + 16*M2 )
     +
      M2GGH = M2GGH + ANG2(13)*CQED*TWO**(-1) * (  - 4*M2 )
     +
      M2GGH = M2GGH + ANG2(14)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 16*
     +    S**(-2)*TG*M2**2 - 16*S**(-2)*TG**2*M2 - 8*S**(-1)*TG**(-1)*
     +    M2**2*MS2 - 8*S**(-1)*TG**(-1)*M2**3 - 16*S**(-1)*TG*M2 + 8*
     +    S**(-1)*UG**(-1)*M2**2*MS2 + 8*S**(-1)*UG**(-1)*M2**3 + 16*
     +    S**(-1)*M2*MS2 - 16*S**(-1)*M2**2 - 2*S*TG**(-1)*M2 - 2*S*
     +    TG**(-1)*MS2 - 2*S*UG**(-1)*M2 - 2*S*UG**(-1)*MS2 - 2*S + 8*
     +    TG**(-1)*M2*MS2 + 4*TG**(-1)*M2**2 + 8*UG**(-1)*M2*MS2 - 8*
     +    UG**(-1)*M2**2 - 4*M2 - 8*MS2 )
     +
      M2GGH = M2GGH + ANG2(14)*N*CO*TWO**(-1) * ( 4 + 32*S**(-2)*TG*M2
     +     - 16*S**(-2)*M2*S4 - 8*S**(-1)*TG**(-1)*M2*MS2 - 6*S**(-1)*
     +    TG**(-1)*M2**2 + 4*S**(-1)*UG**(-1)*M2*S4 + 4*S**(-1)*
     +    UG**(-1)*M2**2 + 2*S**(-1)*UG**(-1)*S4**2 - 4*S**(-1)*M2 - 4*
     +    S**(-1)*S4 - 4*TG**(-1)*M2 - 2*TG**(-1)*MS2 - 12*UG**(-1)*M2
     +     - 2*UG**(-1)*MS2 - 2*UG**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(14)*N*CK*TWO**(-1) * ( 22 + 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 - 12*S**(-2)*TG*M2 + 4*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) + 4*S**(-2)*TG*(TG+UG)**(-1)*S4**2 - 4*S**(-2)*
     +    TG*S4 - 8*S**(-2)*M2*(TG+UG)**(-1)*S4**2 + 12*S**(-2)*M2*S4
     +     - 4*S**(-2)*M2**2*(TG+UG)**(-1)*S4 - 4*S**(-2)*(TG+UG)**(-1)
     +    *S4**3 + 4*S**(-2)*S4**2 + 8*S**(-1)*M2*MS2*(TG+UG)**(-1) + 
     +    24*S**(-1)*M2*(TG+UG)**(-1)*S4 - 22*S**(-1)*M2 + 8*S**(-1)*
     +    M2**2*(TG+UG)**(-1) + 8*S**(-1)*MS2*(TG+UG)**(-1)*S4 - 8*
     +    S**(-1)*MS2 + 16*S**(-1)*(TG+UG)**(-1)*S4**2 - 10*S**(-1)*S4
     +     + 2*S*TG**(-1) + 2*S*UG**(-1) + 16*S*(TG+UG)**(-1) + 2*
     +    TG**(-1)*M2 - 2*TG**(-1)*S4 + 2*UG**(-1)*M2 - 2*UG**(-1)*S4
     +     - 4*M2*(TG+UG)**(-1) + 8*MS2*(TG+UG)**(-1) - 28*
     +    (TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(14)*CQED*TWO**(-1) * (  - 4 - 4*S*
     +    (TG+UG)**(-1) + 4*M2*(TG+UG)**(-1) + 4*(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(18)*N*CO*TWO**(-1) * ( 8*S*UG**(-1)*MS2 - 16
     +    *UG**(-2)*M2*MS2**2 - 16*UG**(-2)*M2**2*MS2 - 8*UG**(-1)*M2*
     +    MS2 )
     +
      M2GGH = M2GGH + ANG2(18)*N*CK*TWO**(-1) * (  - 8*S*UG**(-1)*MS2
     +     + 16*UG**(-2)*M2*MS2**2 + 16*UG**(-2)*M2**2*MS2 + 8*UG**(-1)
     +    *M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(19)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 2*
     +    S**(-1)*TG*M2 + 8*S**(-1)*UG**(-1)*M2**2*MS2 + 8*S**(-1)*
     +    UG**(-1)*M2**3 + 8*S**(-1)*M2*MS2 + 6*S**(-1)*M2**2 - 2*S*
     +    TG**(-1)*M2 - 2*S*TG**(-1)*MS2 + 8*S*UG**(-1)*MS2 - 2*S - 16*
     +    TG**(-1)*UG**(-1)*M2*MS2**2 - 32*TG**(-1)*UG**(-1)*M2**2*MS2
     +     - 16*TG**(-1)*UG**(-1)*M2**3 - 4*TG**(-1)*M2*MS2 - 2*TG - 32
     +    *UG**(-2)*M2*MS2**2 - 48*UG**(-2)*M2**2*MS2 - 16*UG**(-2)*
     +    M2**3 - 8*UG**(-1)*M2*MS2 - 4*UG**(-1)*M2**2 + 2*M2 + 10*MS2
     +     )
     +
      M2GGH = M2GGH + ANG2(19)*N*CO*TWO**(-1) * ( 8*S**(-1)*UG**(-1)*M2
     +    *MS2 + 2*S**(-1)*UG**(-1)*M2*S4 + 4*S**(-1)*UG**(-1)*M2**2 + 
     +    2*S**(-1)*M2 - 2*S*TG**(-1) - 18*TG**(-1)*UG**(-1)*M2*MS2 - 
     +    10*TG**(-1)*UG**(-1)*M2**2 + 6*TG**(-1)*UG**(-1)*MS2*S4 + 2*
     +    TG**(-1)*UG**(-1)*S4**2 - 6*TG**(-1)*M2 - 14*TG**(-1)*MS2 - 2
     +    *TG**(-1)*S4 - 16*UG**(-2)*M2*MS2 - 16*UG**(-2)*M2**2 - 16*
     +    UG**(-1)*M2 - 12*UG**(-1)*MS2 - 4*UG**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(19)*N*CK*TWO**(-1) * ( 12 + 2*S**(-1)*TG*M2*
     +    (TG+UG)**(-1) - 2*S**(-1)*TG*(TG+UG)**(-1)*S4 + 4*S**(-1)*M2*
     +    MS2*(TG+UG)**(-1) - 10*S**(-1)*M2*(TG+UG)**(-1)*S4 + 4*
     +    S**(-1)*M2**2*(TG+UG)**(-1) - 4*S**(-1)*MS2*(TG+UG)**(-1)*S4
     +     + 4*S**(-1)*MS2 + 6*S**(-1)*(TG+UG)**(-1)*S4**2 - 2*S**(-1)*
     +    S4 + 2*S*TG**(-1) + 4*S*UG**(-1)*M2*(S+TG)**(-1) - 8*S*
     +    UG**(-1)*(S+TG)**(-1)*S4 + 4*S*UG**(-1) + 6*S*(TG+UG)**(-1)
     +     + 4*TG**(-1)*M2 - 4*TG**(-1)*S4 - 4*TG*(S+TG)**(-1) + 2*TG*
     +    (TG+UG)**(-1) + 16*UG**(-2)*M2*(S+TG)**(-1)*S4**2 - 16*
     +    UG**(-2)*M2*S4 + 16*UG**(-2)*MS2*(S+TG)**(-1)*S4**2 - 16*
     +    UG**(-2)*MS2*S4 + 8*UG**(-1)*M2*MS2*(S+TG)**(-1) - 28*
     +    UG**(-1)*M2*(S+TG)**(-1)*S4 + 14*UG**(-1)*M2 + 8*UG**(-1)*
     +    M2**2*(S+TG)**(-1) - 24*UG**(-1)*MS2*(S+TG)**(-1)*S4 + 8*
     +    UG**(-1)*MS2 + 12*UG**(-1)*(S+TG)**(-1)*S4**2 - 14*UG**(-1)*
     +    S4 + 8*M2*(S+TG)**(-1) + 10*M2*(TG+UG)**(-1) + 8*MS2*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(19)*N*CK*TWO**(-1) * ( 4*MS2*(TG+UG)**(-1)
     +     - 12*(S+TG)**(-1)*S4 - 12*(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(20)*N*CK*TWO**(-1) * (  - 2*S*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4 - 2*S*TG**(-1)*MS2*(S+UG)**(-1)*S4 - 2*S*
     +    UG**(-1)*M2*(S+TG)**(-1)*S4 - 2*S*UG**(-1)*MS2*(S+TG)**(-1)*
     +    S4 - 2*S*(S+TG)**(-1)*S4 - 2*TG**(-1)*M2*S4 + 4*TG**(-1)*
     +    (S+UG)**(-1)*S4**3 - 2*TG**(-1)*S4**2 - 2*TG*(S+TG)**(-1)*S4
     +     + 2*TG*(S+UG)**(-1)*S4 - 2*UG**(-1)*M2*S4 + 4*UG**(-1)*
     +    (S+TG)**(-1)*S4**3 - 2*UG**(-1)*S4**2 - 2*M2*(S+TG)**(-1)*S4
     +     - 2*(S+TG)**(-1)*S4**2 - 4*(S+UG)**(-1)*S4**2 + 4*S4 )
     +
      M2GGH = M2GGH + ANG2(20)*CQED*TWO**(-1) * ( 2*S*M2*(S+TG)**(-1)*
     +    (S+UG)**(-1)*S4 + 2*S*MS2*(S+TG)**(-1)*(S+UG)**(-1)*S4 - 4*S*
     +    (S+TG)**(-1)*(S+UG)**(-1)*S4**2 + 4*S*(S+TG)**(-1)*S4 + 2*S*
     +    (S+UG)**(-1)*S4 - 2*S**2*(S+TG)**(-1)*(S+UG)**(-1)*S4 + 2*TG*
     +    (S+TG)**(-1)*S4 - 2*TG*(S+UG)**(-1)*S4 + 2*M2*(S+TG)**(-1)*S4
     +     - 4*(S+TG)**(-1)*(S+UG)**(-1)*S4**3 + 2*(S+TG)**(-1)*S4**2
     +     + 4*(S+UG)**(-1)*S4**2 - 4*S4 )
     +
      M2GGH = M2GGH + ANG2(26)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 2 + 2*
     +    S**(-1)*TG + 8*S**(-1)*UG**(-1)*M2*MS2 + 8*S**(-1)*UG**(-1)*
     +    M2**2 - 8*TG**(-1)*UG**(-1)*M2*MS2 - 8*TG**(-1)*UG**(-1)*
     +    M2**2 + 2*TG**(-1)*M2 + 4*TG**(-1)*MS2 + 4*UG**(-1)*M2 + 8*
     +    UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(26)*N*CO*TWO**(-1) * ( 2*S**(-1)*UG**(-1)*M2
     +     - 2*S**(-1)*UG**(-1)*S4 - 2*S**(-1) - 2*TG**(-1)*UG**(-1)*M2
     +     + 4*TG**(-1)*UG**(-1)*MS2 + 4*TG**(-1)*UG**(-1)*S4 - 4*
     +    TG**(-1) - 4*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(26)*N*CK*TWO**(-1) * ( 4*UG**(-1)*M2*
     +    (S+TG)**(-1) - 12*UG**(-1)*(S+TG)**(-1)*S4 + 8*UG**(-1) + 4*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(31)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4 - 8*
     +    S**(-1)*UG**(-1)*M2*MS2 - 8*S**(-1)*UG**(-1)*M2**2 + 8*
     +    TG**(-1)*UG**(-1)*M2*MS2 + 8*TG**(-1)*UG**(-1)*M2**2 + 2*
     +    TG**(-1)*M2 + 4*UG**(-1)*M2 )
     +
      M2GGH = M2GGH + ANG2(31)*N*CO*TWO**(-1) * ( 4*S**(-1)*UG**(-1)*M2
     +     - 4*S**(-1)*UG**(-1)*S4 + 4*S**(-1) + 2*TG**(-1)*UG**(-1)*M2
     +     - 4*TG**(-1) + 8*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(31)*N*CK*TWO**(-1) * ( 4*UG**(-1)*M2*
     +    (S+TG)**(-1) - 12*UG**(-1)*(S+TG)**(-1)*S4 + 8*UG**(-1) + 4*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(32)*N*CO*TWO**(-1) * ( 4*S*MS2 + 4*TG*MS2 - 
     +    8*UG**(-1)*M2**2*MS2 - 4*M2*MS2 - 4*MS2*S4 )
     +
      M2GGH = M2GGH + ANG2(32)*N*CK*TWO**(-1) * (  - 12*S*MS2 - 12*TG*
     +    MS2 + 24*UG**(-1)*M2**2*MS2 + 12*M2*MS2 + 12*MS2*S4 )
     +
      M2GGH = M2GGH + ANG2(32)*CQED*TWO**(-1) * ( 4*S*MS2 + 4*TG*MS2 - 
     +    8*UG**(-1)*M2**2*MS2 - 4*M2*MS2 - 4*MS2*S4 )
     +
      M2GGH = M2GGH + ANG2(33)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 2*S*
     +    TG**(-1)*M2*MS2 + 2*S*TG**(-1)*M2**2 + 2*S*UG**(-1)*M2*MS2 + 
     +    2*S*UG**(-1)*M2**2 + 2*S*M2 + 2*S*MS2 - 8*TG**(-1)*M2**2*MS2
     +     - 4*TG**(-1)*M2**3 - 8*UG**(-1)*M2**2*MS2 - 4*UG**(-1)*M2**3
     +     - 8*M2*MS2 - 4*M2**2 )
     +
      M2GGH = M2GGH + ANG2(33)*N*CO*TWO**(-1) * ( 2*S*TG**(-1)*M2 + 4*S
     +    *TG**(-1)*MS2 + 2*S - 2*TG**(-1)*UG**(-1)*M2*MS2*S4 - 6*
     +    TG**(-1)*UG**(-1)*M2**2*MS2 - 2*TG**(-1)*UG**(-1)*M2**2*S4 - 
     +    2*TG**(-1)*UG**(-1)*M2**3 - 6*TG**(-1)*M2*MS2 - 4*TG**(-1)*M2
     +    *S4 - 2*TG**(-1)*M2**2 - 8*TG**(-1)*MS2*S4 + 2*TG + 10*
     +    UG**(-1)*M2*MS2 - 6*UG**(-1)*M2**2 - 2*M2 + 8*MS2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(33)*N*CK*TWO**(-1) * ( 4*S**(-1)*TG*M2*MS2*
     +    (S+UG)**(-1) - 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 + 4*S**(-1)*TG
     +    *M2 + 2*S**(-1)*TG*M2**2*(S+UG)**(-1) - 12*S**(-1)*TG*MS2*
     +    (S+UG)**(-1)*S4 + 4*S**(-1)*TG*MS2 + 10*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**2 - 8*S**(-1)*TG*S4 + 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-1) + 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) - 8*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4 + 2*S**(-1)*TG**2 + 2*S**(-1)*TG**3*
     +    (S+UG)**(-1) - 16*S**(-1)*UG**(-1)*M2*MS2*S4 + 8*S**(-1)*
     +    UG**(-1)*M2**2*MS2 + 8*S**(-1)*M2*MS2 + 4*S**(-1)*M2*
     +    (S+UG)**(-1)*S4**2 + 2*S**(-1)*M2**2 + 8*S**(-1)*MS2*
     +    (S+UG)**(-1)*S4**2 - 8*S**(-1)*MS2*S4 - 4*S**(-1)*
     +    (S+UG)**(-1)*S4**3 + 6*S**(-1)*S4**2 + 2*S*TG**(-1)*M2*MS2*
     +    (S+UG)**(-1) - 2*S*TG**(-1)*M2*(S+UG)**(-1)*S4 - 4*S*TG**(-1)
     +    *MS2*(S+UG)**(-1)*S4 + 8*S*TG**(-1)*MS2 - 2*S*UG**(-1)*M2*
     +    (S+TG)**(-1)*S4 + 2*S*UG**(-1)*M2 - 2*S*UG**(-1)*MS2*
     +    (S+TG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(33)*N*CK*TWO**(-1) * ( 2*S*UG**(-1)*MS2 + 2*
     +    S*MS2*(S+UG)**(-1) - 2*S*(S+TG)**(-1)*S4 - 2*S*(S+UG)**(-1)*
     +    S4 + 6*TG**(-1)*UG**(-1)*M2*MS2*S4 + 2*TG**(-1)*UG**(-1)*
     +    M2**2*MS2 + 6*TG**(-1)*UG**(-1)*M2**2*S4 + 6*TG**(-1)*
     +    UG**(-1)*M2**3 - 6*TG**(-1)*M2*MS2 + 4*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4**2 + 4*TG**(-1)*M2*S4 - 2*TG**(-1)*M2**2 + 8*
     +    TG**(-1)*MS2*(S+UG)**(-1)*S4**2 + 4*TG*M2*(S+UG)**(-1) + 6*TG
     +    *MS2*(S+UG)**(-1) - 2*TG*(S+TG)**(-1)*S4 - 8*TG*(S+UG)**(-1)*
     +    S4 + 2*TG + 2*TG**2*(S+UG)**(-1) + 16*UG**(-2)*M2**2*MS2 - 14
     +    *UG**(-1)*M2*MS2 - 4*UG**(-1)*M2*S4 + 6*UG**(-1)*M2**2 + 4*
     +    UG**(-1)*(S+TG)**(-1)*S4**3 - 4*UG**(-1)*S4**2 + 6*M2*MS2*
     +    (S+UG)**(-1) - 2*M2*(S+TG)**(-1)*S4 - 8*M2*(S+UG)**(-1)*S4 + 
     +    2*M2 + 2*M2**2*(S+UG)**(-1) - 14*MS2*(S+UG)**(-1)*S4 - 4*MS2
     +     - 2*(S+TG)**(-1)*S4**2 + 6*(S+UG)**(-1)*S4**2 - 6*S4 )
     +
      M2GGH = M2GGH + ANG2(33)*CQED*TWO**(-1) * (  - 4*S**(-1)*TG*M2*
     +    MS2*(S+UG)**(-1) + 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 - 2*
     +    S**(-1)*TG*M2 - 2*S**(-1)*TG*M2**2*(S+UG)**(-1) + 12*S**(-1)*
     +    TG*MS2*(S+UG)**(-1)*S4 - 4*S**(-1)*TG*MS2 - 10*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**2 + 6*S**(-1)*TG*S4 - 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-1) - 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) + 8*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4 - 2*S**(-1)*TG**2 - 2*S**(-1)*TG**3*
     +    (S+UG)**(-1) + 8*S**(-1)*UG**(-1)*M2*MS2*S4 - 4*S**(-1)*M2*
     +    (S+UG)**(-1)*S4**2 - 8*S**(-1)*MS2*(S+UG)**(-1)*S4**2 + 8*
     +    S**(-1)*MS2*S4 + 4*S**(-1)*(S+UG)**(-1)*S4**3 - 4*S**(-1)*
     +    S4**2 - 4*S*TG**(-1)*MS2 + 2*S*M2*(S+TG)**(-1)*(S+UG)**(-1)*
     +    S4 + 2*S*MS2*(S+TG)**(-1)*(S+UG)**(-1)*S4 - 4*S*(S+TG)**(-1)*
     +    (S+UG)**(-1)*S4**2 + 4*S*(S+TG)**(-1)*S4 + 2*S*(S+UG)**(-1)*
     +    S4 - 2*S**2*(S+TG)**(-1)*(S+UG)**(-1)*S4 - 2*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 + 2*TG**(-1)*UG**(-1)*M2**2*MS2 - 2*
     +    TG**(-1)*UG**(-1)*M2**2*S4 )
     +
      M2GGH = M2GGH + ANG2(33)*CQED*TWO**(-1) * (  - 2*TG**(-1)*
     +    UG**(-1)*M2**3 + 2*TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*S4 - 2*TG*
     +    MS2*(S+UG)**(-1) + 2*TG*(S+TG)**(-1)*S4 - 8*UG**(-2)*M2**2*
     +    MS2 + 2*UG**(-1)*M2*MS2 - 2*UG**(-1)*M2**2 - 2*M2*MS2*
     +    (S+UG)**(-1) + 2*M2*(S+TG)**(-1)*S4 + 2*MS2*(S+UG)**(-1)*S4
     +     - 2*MS2 - 4*(S+TG)**(-1)*(S+UG)**(-1)*S4**3 + 2*(S+TG)**(-1)
     +    *S4**2 + 4*(S+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(34)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4 - 8*
     +    S**(-1)*TG**(-1)*M2*MS2 - 8*S**(-1)*TG**(-1)*M2**2 - 4*
     +    S**(-1)*M2 + 8*TG**(-1)*UG**(-1)*M2*MS2 + 8*TG**(-1)*UG**(-1)
     +    *M2**2 + 4*TG**(-1)*M2 + 6*UG**(-1)*M2 )
     +
      M2GGH = M2GGH + ANG2(34)*N*CO*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +    *S4 + 4*S**(-1) + 6*TG**(-1)*UG**(-1)*M2 + 8*TG**(-1) - 4*
     +    UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(34)*N*CK*TWO**(-1) * ( 4*TG**(-1)*M2*
     +    (S+UG)**(-1) - 12*TG**(-1)*(S+UG)**(-1)*S4 + 8*TG**(-1) + 4*
     +    (S+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(36)*N*CO*TWO**(-1) * (  - 2*S*M2 )
     +
      M2GGH = M2GGH + ANG2(36)*N*CK*TWO**(-1) * ( 6*S*M2 )
     +
      M2GGH = M2GGH + ANG2(36)*CQED*TWO**(-1) * (  - 2*S*M2 )
     +
      M2GGH = M2GGH + ANG2(38)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 4*S
     +    *TG - 2*S*UG**(-1)*M2*MS2 - 2*S*UG**(-1)*M2**2 + 4*S*M2 - 2*S
     +    *MS2 - 2*S**2 + 6*TG*M2 - 2*TG**2 - 4*UG**(-1)*M2**3 - 8*
     +    M2**2 )
     +
      M2GGH = M2GGH + ANG2(38)*N*CO*TWO**(-1) * (  - 4*S*UG**(-1)*M2 - 
     +    4*S*UG**(-1)*MS2 + 4*S + 4*TG - 2*UG**(-1)*M2*S4 - 2*UG**(-1)
     +    *M2**2 - 2*M2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(38)*N*CK*TWO**(-1) * (  - 4*S*TG**(-1)*M2*
     +    MS2*(TG+UG)**(-1) + 14*S*TG**(-1)*M2*(TG+UG)**(-1)*S4 - 6*S*
     +    TG**(-1)*M2 - 4*S*TG**(-1)*M2**2*(TG+UG)**(-1) + 6*S*TG**(-1)
     +    *MS2*(TG+UG)**(-1)*S4 - 4*S*TG**(-1)*S4 - 6*S*TG*(S+TG)**(-1)
     +     + 9*S*TG*(TG+UG)**(-1) - 2*S*UG**(-1)*M2*MS2*(S+TG)**(-1) + 
     +    2*S*UG**(-1)*M2*(S+TG)**(-1)*S4 + 2*S*UG**(-1)*M2 + 4*S*
     +    UG**(-1)*MS2*(S+TG)**(-1)*S4 + 8*S*UG**(-1)*MS2 - 2*S*M2*
     +    (S+TG)**(-1) + 3*S*M2*(TG+UG)**(-1) - 2*S*MS2*(S+TG)**(-1) + 
     +    4*S*MS2*(TG+UG)**(-1) - 14*S*(TG+UG)**(-1)*S4 + 15*S - 6*S**2
     +    *TG**(-1)*M2*(TG+UG)**(-1) - 2*S**2*TG**(-1)*MS2*
     +    (TG+UG)**(-1) + 2*S**2*TG**(-1) - 4*S**2*(S+TG)**(-1) + 7*
     +    S**2*(TG+UG)**(-1) + 4*TG**(-1)*UG**(-1)*M2*MS2*S4 - 4*
     +    TG**(-1)*UG**(-1)*M2**2*MS2 + 4*TG**(-1)*UG**(-1)*M2**2*S4 + 
     +    4*TG**(-1)*UG**(-1)*M2**3 + 4*TG**(-1)*M2*MS2*(TG+UG)**(-1)*
     +    S4 - 4*TG**(-1)*M2*MS2 - 8*TG**(-1)*M2*(TG+UG)**(-1)*S4**2 + 
     +    10*TG**(-1)*M2*S4 )
     +
      M2GGH = M2GGH + ANG2(38)*N*CK*TWO**(-1) * ( 4*TG**(-1)*M2**2*
     +    (TG+UG)**(-1)*S4 + 4*TG**(-1)*M2**2 - 4*TG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4**2 + 4*TG**(-1)*MS2*S4 + 2*TG**(-1)*S4**2 - 
     +    2*TG*M2*(S+TG)**(-1) + 3*TG*M2*(TG+UG)**(-1) - 9*TG*
     +    (TG+UG)**(-1)*S4 + 11*TG - 2*TG**2*(S+TG)**(-1) + 2*TG**2*
     +    (TG+UG)**(-1) - 4*UG**(-1)*M2*MS2 + 4*UG**(-1)*M2*
     +    (S+TG)**(-1)*S4**2 + 2*UG**(-1)*M2*S4 - 2*UG**(-1)*M2**2 + 2*
     +    M2*MS2*(TG+UG)**(-1) - 2*M2*(S+TG)**(-1)*S4 - M2*
     +    (TG+UG)**(-1)*S4 - 10*M2 + 2*M2**2*(TG+UG)**(-1) - 2*MS2*
     +    (TG+UG)**(-1)*S4 + 2*MS2 + 2*(S+TG)**(-1)*S4**2 + 7*
     +    (TG+UG)**(-1)*S4**2 - 13*S4 )
     +
      M2GGH = M2GGH + ANG2(38)*CQED*TWO**(-1) * ( 4*S*TG**(-1)*M2*MS2*
     +    (TG+UG)**(-1) - 14*S*TG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*
     +    TG**(-1)*M2 + 4*S*TG**(-1)*M2**2*(TG+UG)**(-1) - 6*S*TG**(-1)
     +    *MS2*(TG+UG)**(-1)*S4 + 2*S*TG**(-1)*MS2 - 2*S*TG*M2*
     +    (S+TG)**(-1)*(TG+UG)**(-1) - 2*S*TG*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) + 2*S*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*TG
     +    *(S+TG)**(-1) - 2*S*TG*(TG+UG)**(-1) - 2*S*TG**2*(S+TG)**(-1)
     +    *(TG+UG)**(-1) - 4*S*UG**(-1)*MS2 - 2*S*M2*(S+TG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S*M2*(S+TG)**(-1) - 2*S*MS2*(S+TG)**(-1)
     +    *(TG+UG)**(-1)*S4 + 8*S*(S+TG)**(-1)*(TG+UG)**(-1)*S4**2 - 2*
     +    S*(S+TG)**(-1)*S4 + 4*S*(TG+UG)**(-1)*S4 - 4*S + 6*S**2*
     +    TG**(-1)*M2*(TG+UG)**(-1) + 2*S**2*TG**(-1)*MS2*(TG+UG)**(-1)
     +     - 6*S**2*TG*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S**2*M2*
     +    (S+TG)**(-1)*(TG+UG)**(-1) - 2*S**2*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) - 4*S**3*(S+TG)**(-1)*(TG+UG)**(-1) - 2*
     +    TG**(-1)*UG**(-1)*M2*MS2*S4 )
     +
      M2GGH = M2GGH + ANG2(38)*CQED*TWO**(-1) * ( 2*TG**(-1)*UG**(-1)*
     +    M2**2*MS2 - 2*TG**(-1)*UG**(-1)*M2**2*S4 - 2*TG**(-1)*
     +    UG**(-1)*M2**3 - 4*TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 2*
     +    TG**(-1)*M2*MS2 + 8*TG**(-1)*M2*(TG+UG)**(-1)*S4**2 - 8*
     +    TG**(-1)*M2*S4 - 4*TG**(-1)*M2**2*(TG+UG)**(-1)*S4 - 2*
     +    TG**(-1)*M2**2 + 4*TG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 4*
     +    TG**(-1)*MS2*S4 + 2*TG*M2*(S+TG)**(-1) - 2*TG*M2*
     +    (TG+UG)**(-1) + 4*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4**2 + 2*TG*
     +    (TG+UG)**(-1)*S4 - 6*TG + 2*TG**2*(S+TG)**(-1) + 2*UG**(-1)*
     +    M2*MS2 - 2*UG**(-1)*M2*S4 - 2*M2*MS2*(TG+UG)**(-1) + 2*M2*
     +    (S+TG)**(-1)*S4 + 2*M2*(TG+UG)**(-1)*S4 - 2*M2**2*
     +    (TG+UG)**(-1) + 2*MS2*(TG+UG)**(-1)*S4 - 2*MS2 - 4*
     +    (S+TG)**(-1)*(TG+UG)**(-1)*S4**3 + 2*(S+TG)**(-1)*S4**2 - 4*
     +    (TG+UG)**(-1)*S4**2 + 4*S4 )
     +
      M2GGH = M2GGH + ANG2(39)*N*CO*TWO**(-1) * (  - 2*S*M2 )
     +
      M2GGH = M2GGH + ANG2(39)*N*CK*TWO**(-1) * ( 6*S*M2 )
     +
      M2GGH = M2GGH + ANG2(39)*CQED*TWO**(-1) * (  - 2*S*M2 )
     +
      M2GGH = M2GGH + ANG2(40)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 2*S
     +    *TG**(-1)*M2*MS2 - 2*S*TG**(-1)*M2**2 - 2*S*M2 - 2*S*MS2 - 4*
     +    TG**(-1)*M2**3 - 6*TG*M2 - 2*TG**2 - 8*M2**2 )
     +
      M2GGH = M2GGH + ANG2(40)*N*CO*TWO**(-1) * (  - 4*S*TG**(-1)*M2 - 
     +    4*S*TG**(-1)*MS2 - 2*TG**(-1)*M2*S4 - 2*TG**(-1)*M2**2 + 2*M2
     +     )
     +
      M2GGH = M2GGH + ANG2(40)*N*CK*TWO**(-1) * (  - 2*S*TG**(-1)*M2*
     +    MS2*(S+UG)**(-1) + 2*S*TG**(-1)*M2*(S+UG)**(-1)*S4 + 2*S*
     +    TG**(-1)*M2 + 4*S*TG**(-1)*MS2*(S+UG)**(-1)*S4 + 8*S*TG**(-1)
     +    *MS2 + 2*S*TG*(S+UG)**(-1) - 5*S*TG*(TG+UG)**(-1) - 4*S*
     +    UG**(-1)*M2*MS2*(TG+UG)**(-1) + 14*S*UG**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 6*S*UG**(-1)*M2 - 4*S*UG**(-1)*M2**2*
     +    (TG+UG)**(-1) + 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 - 4*S*
     +    UG**(-1)*S4 + 2*S*M2*(S+UG)**(-1) - 5*S*M2*(TG+UG)**(-1) - 2*
     +    S*MS2*(S+UG)**(-1) + 4*S*MS2*(TG+UG)**(-1) - 2*S*(S+UG)**(-1)
     +    *S4 + 4*S - 6*S**2*UG**(-1)*M2*(TG+UG)**(-1) - 2*S**2*
     +    UG**(-1)*MS2*(TG+UG)**(-1) + 2*S**2*UG**(-1) + 4*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 - 4*TG**(-1)*UG**(-1)*M2**2*MS2 + 4*
     +    TG**(-1)*UG**(-1)*M2**2*S4 + 4*TG**(-1)*UG**(-1)*M2**3 - 4*
     +    TG**(-1)*M2*MS2 + 4*TG**(-1)*M2*(S+UG)**(-1)*S4**2 + 2*
     +    TG**(-1)*M2*S4 - 2*TG**(-1)*M2**2 - 2*TG*M2*(S+UG)**(-1) + TG
     +    *M2*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(40)*N*CK*TWO**(-1) * ( 4*TG*(S+UG)**(-1)*S4
     +     + 5*TG*(TG+UG)**(-1)*S4 - 11*TG - 2*TG**2*(S+UG)**(-1) + 2*
     +    TG**2*(TG+UG)**(-1) + 4*UG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 - 4*
     +    UG**(-1)*M2*MS2 - 8*UG**(-1)*M2*(TG+UG)**(-1)*S4**2 + 10*
     +    UG**(-1)*M2*S4 + 4*UG**(-1)*M2**2*(TG+UG)**(-1)*S4 + 4*
     +    UG**(-1)*M2**2 - 4*UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 + 4*
     +    UG**(-1)*MS2*S4 + 2*UG**(-1)*S4**2 + 2*M2*MS2*(TG+UG)**(-1)
     +     + 7*M2*(TG+UG)**(-1)*S4 - 21*M2 + M2**2*(TG+UG)**(-1) - 2*
     +    MS2*(TG+UG)**(-1)*S4 + 2*MS2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(40)*CQED*TWO**(-1) * (  - 4*S*TG**(-1)*MS2
     +     - 2*S*TG*M2*(S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*MS2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S*TG*(S+UG)**(-1) + 2*S*TG*(TG+UG)**(-1)
     +     - 2*S*TG**2*(S+UG)**(-1)*(TG+UG)**(-1) + 4*S*UG**(-1)*M2*MS2
     +    *(TG+UG)**(-1) - 14*S*UG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*
     +    UG**(-1)*M2 + 4*S*UG**(-1)*M2**2*(TG+UG)**(-1) - 6*S*UG**(-1)
     +    *MS2*(TG+UG)**(-1)*S4 + 2*S*UG**(-1)*MS2 + 2*S*M2*MS2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) - 2*S*M2*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S*M2*(S+UG)**(-1) + 4*S*M2*(TG+UG)**(-1)
     +     - 4*S*MS2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 4*S*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 - 4*S*(S+UG)**(-1)*S4 + 2*S + 2*S**2*TG*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 6*S**2*UG**(-1)*M2*(TG+UG)**(-1)
     +     + 2*S**2*UG**(-1)*MS2*(TG+UG)**(-1) + 2*S**2*M2*(S+UG)**(-1)
     +    *(TG+UG)**(-1) - 4*S**2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*
     +    S**2*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(40)*CQED*TWO**(-1) * (  - 2*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 + 2*TG**(-1)*UG**(-1)*M2**2*MS2 - 2*
     +    TG**(-1)*UG**(-1)*M2**2*S4 - 2*TG**(-1)*UG**(-1)*M2**3 + 2*
     +    TG**(-1)*M2*MS2 - 2*TG**(-1)*M2*S4 + 2*TG*M2*(S+UG)**(-1) + 2
     +    *TG*M2*(TG+UG)**(-1) - 4*TG*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2
     +     - 4*TG*(S+UG)**(-1)*S4 - 2*TG*(TG+UG)**(-1)*S4 + 6*TG + 2*
     +    TG**2*(S+UG)**(-1) - 4*UG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 2*
     +    UG**(-1)*M2*MS2 + 8*UG**(-1)*M2*(TG+UG)**(-1)*S4**2 - 8*
     +    UG**(-1)*M2*S4 - 4*UG**(-1)*M2**2*(TG+UG)**(-1)*S4 - 2*
     +    UG**(-1)*M2**2 + 4*UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 4*
     +    UG**(-1)*MS2*S4 - 2*M2*MS2*(TG+UG)**(-1) - 4*M2*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 - 2*M2*(TG+UG)**(-1)*S4 + 6*M2 + 2*MS2*
     +    (TG+UG)**(-1)*S4 - 2*MS2 + 4*(S+UG)**(-1)*S4**2 - 2*
     +    (TG+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(41)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 16*
     +    S**(-2)*TG*M2 - 8*S**(-1)*TG**(-1)*M2*MS2 - 8*S**(-1)*
     +    TG**(-1)*M2**2 + 8*S**(-1)*UG**(-1)*M2*MS2 + 8*S**(-1)*
     +    UG**(-1)*M2**2 - 8*S**(-1)*M2 + 2*TG**(-1)*M2 - 10*UG**(-1)*
     +    M2 )
     +
      M2GGH = M2GGH + ANG2(41)*N*CO*TWO**(-1) * (  - 16*S**(-2)*M2 + 4*
     +    S**(-1)*UG**(-1)*M2 + 4*S**(-1)*UG**(-1)*S4 - 8*S**(-1) - 4*
     +    UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(41)*N*CK*TWO**(-1) * ( 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) + 8*S**(-2)*TG*(TG+UG)**(-1)*S4 - 8*S**(-2)*TG
     +     - 16*S**(-2)*M2*(TG+UG)**(-1)*S4 + 12*S**(-2)*M2 - 4*S**(-2)
     +    *M2**2*(TG+UG)**(-1) - 12*S**(-2)*(TG+UG)**(-1)*S4**2 + 12*
     +    S**(-2)*S4 + 12*S**(-1)*M2*(TG+UG)**(-1) + 12*S**(-1)*
     +    (TG+UG)**(-1)*S4 - 4*S**(-1) )
     +
      M2GGH = M2GGH + ANG2(43)*N*CO*TWO**(-1) * ( 16*S**(-2)*TG*M2 - 16
     +    *S**(-2)*M2*S4 + 8*S**(-1)*M2 )
     +
      M2GGH = M2GGH + ANG2(43)*N*CK*TWO**(-1) * (  - 16*S**(-2)*TG*M2
     +     + 16*S**(-2)*M2*S4 - 8*S**(-1)*M2 )
     +
      M2GGH = M2GGH + ANG2(44)*N*CK*TWO**(-1) * ( 4*S**(-1)*TG*M2*MS2*
     +    (TG+UG)**(-1) + 2*S**(-1)*TG*M2**2*(TG+UG)**(-1) + 4*S**(-1)*
     +    TG*MS2*(TG+UG)**(-1)*S4 - 4*S**(-1)*TG*MS2 - 2*S**(-1)*TG*
     +    (TG+UG)**(-1)*S4**2 - 4*S**(-1)*M2*(TG+UG)**(-1)*S4**2 + 2*
     +    S**(-1)*M2**2 + 4*S**(-1)*(TG+UG)**(-1)*S4**3 - 2*S**(-1)*
     +    S4**2 - 8*S*TG**(-1)*MS2 - 6*S*TG*(S+TG)**(-1) + 2*S*TG*
     +    (TG+UG)**(-1) - 2*S*UG**(-1)*M2*MS2*(S+TG)**(-1) + 4*S*
     +    UG**(-1)*M2*MS2*(TG+UG)**(-1) + 2*S*UG**(-1)*M2*(S+TG)**(-1)*
     +    S4 - 14*S*UG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*UG**(-1)*M2 + 4*
     +    S*UG**(-1)*M2**2*(TG+UG)**(-1) + 4*S*UG**(-1)*MS2*
     +    (S+TG)**(-1)*S4 - 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 - 2*S*M2*
     +    (S+TG)**(-1) + 2*S*M2*(TG+UG)**(-1) - 2*S*MS2*(S+TG)**(-1) + 
     +    2*S*MS2*(TG+UG)**(-1) - 4*S*(TG+UG)**(-1)*S4 + 8*S + 6*S**2*
     +    UG**(-1)*M2*(TG+UG)**(-1) + 2*S**2*UG**(-1)*MS2*(TG+UG)**(-1)
     +     - 4*S**2*(S+TG)**(-1) + 4*S**2*(TG+UG)**(-1) - 4*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 )
     +
      M2GGH = M2GGH + ANG2(44)*N*CK*TWO**(-1) * ( 4*TG**(-1)*UG**(-1)*
     +    M2**2*MS2 - 4*TG**(-1)*UG**(-1)*M2**2*S4 - 4*TG**(-1)*
     +    UG**(-1)*M2**3 + 4*TG**(-1)*M2*MS2 - 4*TG**(-1)*M2*S4 - 2*TG*
     +    M2*(S+TG)**(-1) + 2*TG*M2*(TG+UG)**(-1) - 2*TG*(TG+UG)**(-1)*
     +    S4 + 4*TG - 2*TG**2*(S+TG)**(-1) - 4*UG**(-1)*M2*MS2*
     +    (TG+UG)**(-1)*S4 + 12*UG**(-1)*M2*MS2 + 4*UG**(-1)*M2*
     +    (S+TG)**(-1)*S4**2 + 8*UG**(-1)*M2*(TG+UG)**(-1)*S4**2 - 12*
     +    UG**(-1)*M2*S4 - 4*UG**(-1)*M2**2*(TG+UG)**(-1)*S4 + 4*
     +    UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 4*UG**(-1)*MS2*S4 + 4*M2*
     +    MS2*(TG+UG)**(-1) - 2*M2*(S+TG)**(-1)*S4 + 6*M2*(TG+UG)**(-1)
     +    *S4 + 2*M2**2*(TG+UG)**(-1) + 4*MS2*(TG+UG)**(-1)*S4 - 6*MS2
     +     + 2*(S+TG)**(-1)*S4**2 - 4*(TG+UG)**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(44)*CQED*TWO**(-1) * ( 4*S*TG**(-1)*MS2 - 2*
     +    S*TG*M2*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S*TG*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) + 2*S*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*TG
     +    *(S+TG)**(-1) - 2*S*TG**2*(S+TG)**(-1)*(TG+UG)**(-1) - 4*S*
     +    UG**(-1)*M2*MS2*(TG+UG)**(-1) + 14*S*UG**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 6*S*UG**(-1)*M2 - 4*S*UG**(-1)*M2**2*
     +    (TG+UG)**(-1) + 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 - 2*S*
     +    UG**(-1)*MS2 - 2*S*M2*(S+TG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*M2*
     +    (S+TG)**(-1) - 2*S*M2*(TG+UG)**(-1) - 2*S*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 4*S*MS2*(TG+UG)**(-1) + 8*S*(S+TG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 - 2*S*(S+TG)**(-1)*S4 + 4*S*(TG+UG)**(-1)
     +    *S4 - 4*S - 6*S**2*TG*(S+TG)**(-1)*(TG+UG)**(-1) - 6*S**2*
     +    UG**(-1)*M2*(TG+UG)**(-1) - 2*S**2*UG**(-1)*MS2*(TG+UG)**(-1)
     +     - 2*S**2*M2*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S**2*MS2*
     +    (S+TG)**(-1)*(TG+UG)**(-1) - 4*S**3*(S+TG)**(-1)*
     +    (TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(44)*CQED*TWO**(-1) * ( 2*TG**(-1)*UG**(-1)*
     +    M2*MS2*S4 - 2*TG**(-1)*UG**(-1)*M2**2*MS2 + 2*TG**(-1)*
     +    UG**(-1)*M2**2*S4 + 2*TG**(-1)*UG**(-1)*M2**3 - 2*TG**(-1)*M2
     +    *MS2 + 2*TG**(-1)*M2*S4 + 2*TG*M2*(S+TG)**(-1) - 2*TG*M2*
     +    (TG+UG)**(-1) + 4*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4**2 + 2*TG*
     +    (TG+UG)**(-1)*S4 - 4*TG + 2*TG**2*(S+TG)**(-1) + 4*UG**(-1)*
     +    M2*MS2*(TG+UG)**(-1)*S4 - 2*UG**(-1)*M2*MS2 - 8*UG**(-1)*M2*
     +    (TG+UG)**(-1)*S4**2 + 8*UG**(-1)*M2*S4 + 4*UG**(-1)*M2**2*
     +    (TG+UG)**(-1)*S4 + 2*UG**(-1)*M2**2 - 4*UG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4**2 + 4*UG**(-1)*MS2*S4 - 2*M2*MS2*
     +    (TG+UG)**(-1) + 2*M2*(S+TG)**(-1)*S4 + 2*M2*(TG+UG)**(-1)*S4
     +     - 2*M2 - 2*M2**2*(TG+UG)**(-1) - 6*MS2*(TG+UG)**(-1)*S4 + 6*
     +    MS2 - 4*(S+TG)**(-1)*(TG+UG)**(-1)*S4**3 + 2*(S+TG)**(-1)*
     +    S4**2 - 4*(TG+UG)**(-1)*S4**2 + 2*S4 )
     +
      M2GGH = M2GGH + ANG2(45)*N*CO*TWO**(-1) * (  - 8*S**(-2)*M2 )
     +
      M2GGH = M2GGH + ANG2(45)*N*CK*TWO**(-1) * ( 8*S**(-2)*M2 )
     +
      M2GGH = M2GGH + ANG2(46)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 8*
     +    S**(-1) - 2*TG**(-1) - 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(46)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG**(-1) + 
     +    2*S**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(46)*N*CK*TWO**(-1) * (  - 8*S**(-2)*M2*
     +    (TG+UG)**(-1) - 8*S**(-2)*(TG+UG)**(-1)*S4 + 8*S**(-2) + 8*
     +    S**(-1)*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(47)*N*CO*TWO**(-1) * (  - 16*TG**(-1)*M2*
     +    MS2**2 - 8*TG**(-1)*M2**2*MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(47)*N*CK*TWO**(-1) * ( 16*TG**(-1)*M2*MS2**2
     +     + 24*TG**(-1)*M2**2*MS2 + 24*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(47)*CQED*TWO**(-1) * (  - 8*TG**(-1)*M2**2*
     +    MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(48)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 8*
     +    S**(-1)*TG*M2*MS2 - 2*S*TG**(-1)*M2*MS2 - 2*S*TG**(-1)*M2**2
     +     - 2*S*TG - 4*S*M2 - 2*S*MS2 - 16*TG**(-1)*M2*MS2**2 - 16*
     +    TG**(-1)*M2**2*MS2 - 4*TG**(-1)*M2**3 - 4*TG*M2 + 4*TG*MS2 - 
     +    8*M2*MS2 - 8*M2**2 )
     +
      M2GGH = M2GGH + ANG2(48)*N*CO*TWO**(-1) * ( 8*S**(-1)*TG**(-1)*M2
     +    *MS2*S4 - 8*S**(-1)*TG**(-1)*M2**2*MS2 - 8*S**(-1)*M2*MS2 - 
     +    16*TG**(-1)*M2*MS2 - 4*TG**(-1)*M2**2 - 4*M2 - 4*MS2 )
     +
      M2GGH = M2GGH + ANG2(48)*N*CK*TWO**(-1) * (  - 32*S**(-2)*M2*
     +    MS2**2 - 24*S**(-1)*TG**(-1)*M2*MS2*S4 - 16*S**(-1)*TG**(-1)*
     +    M2*MS2**2 + 8*S**(-1)*TG**(-1)*M2**2*MS2 - 4*S**(-1)*TG*M2*
     +    MS2*(S+UG)**(-1) + 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 - 8*
     +    S**(-1)*TG*M2*(TG+UG)**(-1)*S4 - 2*S**(-1)*TG*M2 - 2*S**(-1)*
     +    TG*M2**2*(S+UG)**(-1) + 12*S**(-1)*TG*MS2*(S+UG)**(-1)*S4 - 8
     +    *S**(-1)*TG*MS2*(TG+UG)**(-1)*S4 + 4*S**(-1)*TG*MS2 - 10*
     +    S**(-1)*TG*(S+UG)**(-1)*S4**2 + 8*S**(-1)*TG*(TG+UG)**(-1)*
     +    S4**2 - 2*S**(-1)*TG*S4 - 4*S**(-1)*TG**2*M2*(S+UG)**(-1) - 4
     +    *S**(-1)*TG**2*MS2*(S+UG)**(-1) + 8*S**(-1)*TG**2*
     +    (S+UG)**(-1)*S4 - 2*S**(-1)*TG**2 - 2*S**(-1)*TG**3*
     +    (S+UG)**(-1) - 8*S**(-1)*UG**(-1)*M2*MS2*S4 - 16*S**(-1)*
     +    UG**(-1)*M2*MS2**2 - 4*S**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 8*
     +    S**(-1)*M2*MS2 - 4*S**(-1)*M2*(S+UG)**(-1)*S4**2 + 8*S**(-1)*
     +    M2*(TG+UG)**(-1)*S4**2 - 4*S**(-1)*M2*S4 - 4*S**(-1)*M2**2*
     +    (TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(48)*N*CK*TWO**(-1) * (  - 8*S**(-1)*MS2*
     +    (S+UG)**(-1)*S4**2 + 4*S**(-1)*MS2*(TG+UG)**(-1)*S4**2 + 4*
     +    S**(-1)*MS2*S4 + 4*S**(-1)*(S+UG)**(-1)*S4**3 - 4*S**(-1)*
     +    (TG+UG)**(-1)*S4**3 + 4*S*TG**(-1)*M2*MS2*(TG+UG)**(-1) + 2*S
     +    *TG**(-1)*M2*(S+UG)**(-1)*S4 - 14*S*TG**(-1)*M2*(TG+UG)**(-1)
     +    *S4 + 4*S*TG**(-1)*M2 + 4*S*TG**(-1)*M2**2*(TG+UG)**(-1) + 2*
     +    S*TG**(-1)*MS2*(S+UG)**(-1)*S4 - 6*S*TG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4 + 2*S*TG*(TG+UG)**(-1) - 2*S*UG**(-1)*M2 + 6
     +    *S*UG**(-1)*MS2 + 2*S*UG**(-1)*S4 - 2*S*M2*(TG+UG)**(-1) - 2*
     +    S*MS2*(TG+UG)**(-1) + 2*S*(S+UG)**(-1)*S4 - 6*S*(TG+UG)**(-1)
     +    *S4 - 2*S + 6*S**2*TG**(-1)*M2*(TG+UG)**(-1) + 2*S**2*
     +    TG**(-1)*MS2*(TG+UG)**(-1) - 2*S**2*UG**(-1) + 2*S**2*
     +    (TG+UG)**(-1) + 16*TG**(-2)*M2**2*MS2 + 2*TG**(-1)*UG**(-1)*
     +    M2*MS2*S4 + 6*TG**(-1)*UG**(-1)*M2**2*MS2 + 2*TG**(-1)*
     +    UG**(-1)*M2**2*S4 + 2*TG**(-1)*UG**(-1)*M2**3 - 4*TG**(-1)*M2
     +    *MS2*(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(48)*N*CK*TWO**(-1) * ( 6*TG**(-1)*M2*MS2 + 8
     +    *TG**(-1)*M2*(S+UG)**(-1)*S4**2 + 8*TG**(-1)*M2*(TG+UG)**(-1)
     +    *S4**2 - 12*TG**(-1)*M2*S4 - 4*TG**(-1)*M2**2*(TG+UG)**(-1)*
     +    S4 + 2*TG**(-1)*M2**2 + 8*TG**(-1)*MS2*(S+UG)**(-1)*S4**2 + 4
     +    *TG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 12*TG**(-1)*MS2*S4 - 4*
     +    TG**(-1)*(S+UG)**(-1)*S4**3 + 4*TG**(-1)*S4**2 + 4*TG*M2*
     +    (S+UG)**(-1) + 4*TG*M2*(TG+UG)**(-1) + 2*TG*MS2*(S+UG)**(-1)
     +     + 4*TG*MS2*(TG+UG)**(-1) - 8*TG*(S+UG)**(-1)*S4 - 8*TG*
     +    (TG+UG)**(-1)*S4 + 2*TG + 2*TG**2*(S+UG)**(-1) + 14*UG**(-1)*
     +    M2*MS2 + 4*UG**(-1)*M2*S4 + 2*UG**(-1)*M2**2 + 2*M2*MS2*
     +    (S+UG)**(-1) - 8*M2*(S+UG)**(-1)*S4 - 2*M2 + 2*M2**2*
     +    (S+UG)**(-1) - 10*MS2*(S+UG)**(-1)*S4 + 4*MS2*(TG+UG)**(-1)*
     +    S4 + 6*MS2 + 10*(S+UG)**(-1)*S4**2 + 8*(TG+UG)**(-1)*S4**2 - 
     +    10*S4 )
     +
      M2GGH = M2GGH + ANG2(48)*CQED*TWO**(-1) * ( 16*S**(-2)*M2*MS2**2
     +     + 8*S**(-1)*TG**(-1)*M2*MS2*S4 + 4*S**(-1)*TG*M2*MS2*
     +    (S+UG)**(-1) - 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 + 2*S**(-1)*TG
     +    *M2 + 2*S**(-1)*TG*M2**2*(S+UG)**(-1) - 12*S**(-1)*TG*MS2*
     +    (S+UG)**(-1)*S4 + 4*S**(-1)*TG*MS2 + 10*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**2 - 6*S**(-1)*TG*S4 + 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-1) + 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) - 8*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4 + 2*S**(-1)*TG**2 + 2*S**(-1)*TG**3*
     +    (S+UG)**(-1) + 4*S**(-1)*M2*(S+UG)**(-1)*S4**2 + 8*S**(-1)*
     +    MS2*(S+UG)**(-1)*S4**2 - 8*S**(-1)*MS2*S4 - 4*S**(-1)*
     +    (S+UG)**(-1)*S4**3 + 4*S**(-1)*S4**2 - 4*S*TG**(-1)*M2*MS2*
     +    (TG+UG)**(-1) + 14*S*TG**(-1)*M2*(TG+UG)**(-1)*S4 - 6*S*
     +    TG**(-1)*M2 - 4*S*TG**(-1)*M2**2*(TG+UG)**(-1) + 6*S*TG**(-1)
     +    *MS2*(TG+UG)**(-1)*S4 - 2*S*TG**(-1)*MS2 - 2*S*TG*M2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*MS2*(S+UG)**(-1)*
     +    (TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(48)*CQED*TWO**(-1) * ( 2*S*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S*TG*(S+UG)**(-1) - 2*S*TG**2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 2*S*M2*MS2*(S+UG)**(-1)*
     +    (TG+UG)**(-1) - 2*S*M2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*M2
     +    *(S+UG)**(-1) - 4*S*MS2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 4*S*
     +    MS2*(TG+UG)**(-1) + 4*S*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*
     +    S*(S+UG)**(-1)*S4 + 2*S*(TG+UG)**(-1)*S4 - 6*S**2*TG**(-1)*M2
     +    *(TG+UG)**(-1) - 2*S**2*TG**(-1)*MS2*(TG+UG)**(-1) + 2*S**2*
     +    TG*(S+UG)**(-1)*(TG+UG)**(-1) + 2*S**2*M2*(S+UG)**(-1)*
     +    (TG+UG)**(-1) - 4*S**2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 - 8*
     +    TG**(-2)*M2**2*MS2 + 4*TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 - 8*
     +    TG**(-1)*M2*(TG+UG)**(-1)*S4**2 + 8*TG**(-1)*M2*S4 + 4*
     +    TG**(-1)*M2**2*(TG+UG)**(-1)*S4 - 4*TG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4**2 + 4*TG**(-1)*MS2*S4 + 2*TG*M2*
     +    (TG+UG)**(-1) + 2*TG*MS2*(S+UG)**(-1) - 4*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(48)*CQED*TWO**(-1) * (  - 2*TG*(TG+UG)**(-1)
     +    *S4 + 2*TG + 2*M2*MS2*(S+UG)**(-1) - 2*M2*MS2*(TG+UG)**(-1)
     +     - 4*M2*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 2*M2*(TG+UG)**(-1)
     +    *S4 - 2*MS2*(S+UG)**(-1)*S4 - 6*MS2*(TG+UG)**(-1)*S4 + 4*MS2
     +     + 4*(S+UG)**(-1)*S4**2 - 2*(TG+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(49)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 12*
     +    S**(-1)*M2 + 24*S**(-1)*MS2 + S*TG**(-1) + 3*TG**(-1)*M2 + 8*
     +    TG**(-1)*MS2 - UG**(-1)*M2 + 4*UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(49)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG**(-1)*M2
     +     - 6*S**(-1)*TG**(-1)*MS2 - 4*S**(-1)*TG**(-1)*S4 - 2*S**(-1)
     +    *UG**(-1)*M2 - 6*S**(-1)*UG**(-1)*MS2 - 2*S**(-1)*UG**(-1)*S4
     +     + 16*S**(-1) + 3*TG**(-1) )
     +
      M2GGH = M2GGH + ANG2(49)*N*CK*TWO**(-1) * ( 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) + 8*S**(-2)*TG*(TG+UG)**(-1)*S4 - 8*S**(-2)*TG
     +     + 9*S**(-1)*M2*(TG+UG)**(-1) + 8*S**(-1)*MS2*(TG+UG)**(-1)
     +     + 7*S**(-1)*(TG+UG)**(-1)*S4 - 6*S**(-1) - 7*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(50)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 4*
     +    S**(-2)*TG + 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*MS2
     +     - 2*S**(-1)*UG**(-1)*M2 - 2*S**(-1)*UG**(-1)*MS2 - 8*S**(-1)
     +     + 16*TG**(-2)*M2 + 16*TG**(-2)*MS2 + 8*TG**(-1)*UG**(-1)*M2
     +     + 8*TG**(-1)*UG**(-1)*MS2 - TG**(-1) - 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(50)*N*CO*TWO**(-1) * (  - 4*S**(-1)*TG**(-1)
     +     - 5*S**(-1)*UG**(-1) + 12*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(50)*N*CK*TWO**(-1) * (  - 4*S**(-2)*TG*
     +    (TG+UG)**(-1) - 8*S**(-2)*M2*(TG+UG)**(-1) - 8*S**(-2)*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2) + 8*S**(-1)*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(51)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 4*
     +    S**(-2)*TG*M2**2 + 2*S**(-1)*TG**(-1)*M2**2*MS2 + 2*S**(-1)*
     +    TG**(-1)*M2**3 - 2*S**(-1)*UG**(-1)*M2**2*MS2 - 2*S**(-1)*
     +    UG**(-1)*M2**3 - 40*S**(-1)*M2*MS2 - 28*S**(-1)*M2**2 + 4*S*
     +    TG**(-1)*MS2 - S*UG**(-1)*M2 + 4*S*UG**(-1)*MS2 + 2*S - 4*
     +    TG**(-1)*M2*MS2 - 4*TG**(-1)*M2**2 - 8*UG**(-1)*M2*MS2 - 5*
     +    UG**(-1)*M2**2 - 4*M2 )
     +
      M2GGH = M2GGH + ANG2(51)*N*CO*TWO**(-1) * ( 8 + 2*S**(-1)*
     +    TG**(-1)*M2*MS2 + 4*S**(-1)*TG**(-1)*M2**2 + 2*S**(-1)*
     +    UG**(-1)*M2*MS2 + 2*S**(-1)*UG**(-1)*M2*S4 - S**(-1)*UG**(-1)
     +    *M2**2 - 8*S**(-1)*M2 + S*TG**(-1) + S*UG**(-1) + 3*TG**(-1)*
     +    M2 - 2*TG**(-1)*S4 - 2*UG**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(51)*N*CK*TWO**(-1) * ( 5 + 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 - 8*S**(-2)*TG*M2 + 4*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) + 37*S**(-1)*M2*(TG+UG)**(-1)*S4 - 27*S**(-1)*
     +    M2 + 2*S**(-1)*M2**2*(TG+UG)**(-1) + 24*S**(-1)*MS2*
     +    (TG+UG)**(-1)*S4 - 24*S**(-1)*MS2 + 9*S**(-1)*(TG+UG)**(-1)*
     +    S4**2 - 11*S**(-1)*S4 - 2*S*TG**(-1) - 2*S*UG**(-1) + 6*S*
     +    (TG+UG)**(-1) - 2*TG**(-1)*M2 + 2*TG**(-1)*S4 - 2*UG**(-1)*M2
     +     + 2*UG**(-1)*S4 - 20*M2*(TG+UG)**(-1) - 15*(TG+UG)**(-1)*S4
     +     )
     +
      M2GGH = M2GGH + ANG2(51)*CQED*TWO**(-1) * (  - 4 - 4*S*
     +    (TG+UG)**(-1) + 4*M2*(TG+UG)**(-1) + 4*(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(52)*N*CK*TWO**(-1) * (  - 4*S**(-1)*TG*M2*
     +    MS2*(TG+UG)**(-1) - 2*S**(-1)*TG*M2**2*(TG+UG)**(-1) - 4*
     +    S**(-1)*TG*MS2*(TG+UG)**(-1)*S4 + 4*S**(-1)*TG*MS2 + 2*
     +    S**(-1)*TG*(TG+UG)**(-1)*S4**2 + 4*S**(-1)*M2*MS2 - 2*S**(-1)
     +    *M2*(TG+UG)**(-1)*S4**2 - 4*S**(-1)*M2**2*MS2*(TG+UG)**(-1)
     +     + 2*S**(-1)*M2**2*(TG+UG)**(-1)*S4 + 2*S**(-1)*M2**2 - 2*
     +    S**(-1)*M2**3*(TG+UG)**(-1) + 4*S**(-1)*MS2*(TG+UG)**(-1)*
     +    S4**2 - 4*S**(-1)*MS2*S4 + 2*S**(-1)*(TG+UG)**(-1)*S4**3 - 2*
     +    S**(-1)*S4**2 - 2*S*TG**(-1)*M2*MS2*(S+UG)**(-1) + 4*S*
     +    TG**(-1)*M2*MS2*(TG+UG)**(-1) + 2*S*TG**(-1)*M2*(S+UG)**(-1)*
     +    S4 - 14*S*TG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*TG**(-1)*M2 + 4*
     +    S*TG**(-1)*M2**2*(TG+UG)**(-1) + 4*S*TG**(-1)*MS2*
     +    (S+UG)**(-1)*S4 - 6*S*TG**(-1)*MS2*(TG+UG)**(-1)*S4 + 2*S*TG*
     +    (S+UG)**(-1) - 2*S*TG*(TG+UG)**(-1) - 8*S*UG**(-1)*MS2 + 2*S*
     +    M2*(S+UG)**(-1) - 2*S*M2*(TG+UG)**(-1) - 2*S*MS2*(S+UG)**(-1)
     +     + 2*S*MS2*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(52)*N*CK*TWO**(-1) * (  - 2*S*(S+UG)**(-1)*
     +    S4 + 4*S + 6*S**2*TG**(-1)*M2*(TG+UG)**(-1) + 2*S**2*TG**(-1)
     +    *MS2*(TG+UG)**(-1) + 2*S**2*(TG+UG)**(-1) - 4*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 + 4*TG**(-1)*UG**(-1)*M2**2*MS2 - 4*
     +    TG**(-1)*UG**(-1)*M2**2*S4 - 4*TG**(-1)*UG**(-1)*M2**3 - 4*
     +    TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 12*TG**(-1)*M2*MS2 + 4*
     +    TG**(-1)*M2*(S+UG)**(-1)*S4**2 + 8*TG**(-1)*M2*(TG+UG)**(-1)*
     +    S4**2 - 12*TG**(-1)*M2*S4 - 4*TG**(-1)*M2**2*(TG+UG)**(-1)*S4
     +     + 4*TG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 4*TG**(-1)*MS2*S4 - 2
     +    *TG*M2*(S+UG)**(-1) - 2*TG*M2*(TG+UG)**(-1) + 4*TG*
     +    (S+UG)**(-1)*S4 + 2*TG*(TG+UG)**(-1)*S4 - 4*TG - 2*TG**2*
     +    (S+UG)**(-1) + 4*UG**(-1)*M2*MS2 - 4*UG**(-1)*M2*S4 + 10*M2*
     +    (TG+UG)**(-1)*S4 - 4*M2 - 2*M2**2*(TG+UG)**(-1) - 2*MS2 - 4*
     +    (TG+UG)**(-1)*S4**2 + 4*S4 )
     +
      M2GGH = M2GGH + ANG2(52)*CQED*TWO**(-1) * (  - 4*S*TG**(-1)*M2*
     +    MS2*(TG+UG)**(-1) + 14*S*TG**(-1)*M2*(TG+UG)**(-1)*S4 - 6*S*
     +    TG**(-1)*M2 - 4*S*TG**(-1)*M2**2*(TG+UG)**(-1) + 6*S*TG**(-1)
     +    *MS2*(TG+UG)**(-1)*S4 - 2*S*TG**(-1)*MS2 - 2*S*TG*M2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*MS2*(S+UG)**(-1)*
     +    (TG+UG)**(-1) + 2*S*TG*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*TG
     +    *(S+UG)**(-1) - 2*S*TG**2*(S+UG)**(-1)*(TG+UG)**(-1) + 4*S*
     +    UG**(-1)*MS2 + 2*S*M2*MS2*(S+UG)**(-1)*(TG+UG)**(-1) - 2*S*M2
     +    *(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*M2*(S+UG)**(-1) - 4*S*
     +    MS2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 4*S*MS2*(TG+UG)**(-1) + 4
     +    *S*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*S*(S+UG)**(-1)*S4 + 2
     +    *S*(TG+UG)**(-1)*S4 - 6*S**2*TG**(-1)*M2*(TG+UG)**(-1) - 2*
     +    S**2*TG**(-1)*MS2*(TG+UG)**(-1) + 2*S**2*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1) + 2*S**2*M2*(S+UG)**(-1)*(TG+UG)**(-1) - 4*S**2
     +    *(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*TG**(-1)*UG**(-1)*M2*MS2*
     +    S4 )
     +
      M2GGH = M2GGH + ANG2(52)*CQED*TWO**(-1) * (  - 2*TG**(-1)*
     +    UG**(-1)*M2**2*MS2 + 2*TG**(-1)*UG**(-1)*M2**2*S4 + 2*
     +    TG**(-1)*UG**(-1)*M2**3 + 4*TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4
     +     - 2*TG**(-1)*M2*MS2 - 8*TG**(-1)*M2*(TG+UG)**(-1)*S4**2 + 8*
     +    TG**(-1)*M2*S4 + 4*TG**(-1)*M2**2*(TG+UG)**(-1)*S4 + 2*
     +    TG**(-1)*M2**2 - 4*TG**(-1)*MS2*(TG+UG)**(-1)*S4**2 + 4*
     +    TG**(-1)*MS2*S4 + 2*TG*M2*(S+UG)**(-1) + 2*TG*M2*
     +    (TG+UG)**(-1) - 4*TG*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*TG*
     +    (S+UG)**(-1)*S4 - 2*TG*(TG+UG)**(-1)*S4 + 4*TG + 2*TG**2*
     +    (S+UG)**(-1) - 2*UG**(-1)*M2*MS2 + 2*UG**(-1)*M2*S4 - 2*M2*
     +    MS2*(TG+UG)**(-1) - 4*M2*(S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 2
     +    *M2*(TG+UG)**(-1)*S4 + 2*M2 - 6*MS2*(TG+UG)**(-1)*S4 + 6*MS2
     +     + 4*(S+UG)**(-1)*S4**2 - 2*(TG+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(53)*N*CO*TWO**(-1) * ( 8*S*TG**(-1)*MS2 - 16
     +    *TG**(-2)*M2*MS2**2 - 16*TG**(-2)*M2**2*MS2 - 16*TG**(-1)*M2*
     +    MS2 )
     +
      M2GGH = M2GGH + ANG2(53)*N*CK*TWO**(-1) * (  - 8*S*TG**(-1)*MS2
     +     + 16*TG**(-2)*M2*MS2**2 + 16*TG**(-2)*M2**2*MS2 + 16*
     +    TG**(-1)*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(54)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4*
     +    S**(-1)*TG*M2 + 8*S**(-1)*M2*MS2 + 4*S**(-1)*M2**2 + 8*S*
     +    TG**(-1)*MS2 - 2*S*UG**(-1)*M2 - 2*S*UG**(-1)*MS2 - 32*
     +    TG**(-2)*M2*MS2**2 - 48*TG**(-2)*M2**2*MS2 - 16*TG**(-2)*
     +    M2**3 - 16*TG**(-1)*UG**(-1)*M2*MS2**2 - 24*TG**(-1)*UG**(-1)
     +    *M2**2*MS2 - 8*TG**(-1)*UG**(-1)*M2**3 - 16*TG**(-1)*M2*MS2
     +     - 8*TG**(-1)*M2**2 + 2*TG - 8*UG**(-1)*M2*MS2 + 4*M2 + 10*
     +    MS2 )
     +
      M2GGH = M2GGH + ANG2(54)*N*CO*TWO**(-1) * (  - 2 + 8*S**(-1)*
     +    TG**(-1)*M2*MS2 + 4*S**(-1)*TG**(-1)*M2*S4 - 2*S*UG**(-1) - 
     +    16*TG**(-2)*M2*MS2 - 16*TG**(-2)*M2**2 - 22*TG**(-1)*UG**(-1)
     +    *M2*MS2 - 4*TG**(-1)*UG**(-1)*M2*S4 - 6*TG**(-1)*UG**(-1)*
     +    M2**2 + 6*TG**(-1)*UG**(-1)*MS2*S4 + 2*TG**(-1)*UG**(-1)*
     +    S4**2 - 12*TG**(-1)*M2 - 12*TG**(-1)*MS2 - 4*TG**(-1)*S4 - 2*
     +    UG**(-1)*M2 - 14*UG**(-1)*MS2 - 2*UG**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(54)*N*CK*TWO**(-1) * ( 12 - 2*S**(-1)*TG*M2*
     +    (TG+UG)**(-1) + 2*S**(-1)*TG*(TG+UG)**(-1)*S4 + 4*S**(-1)*M2*
     +    MS2*(TG+UG)**(-1) - 6*S**(-1)*M2*(TG+UG)**(-1)*S4 + 2*S**(-1)
     +    *M2**2*(TG+UG)**(-1) - 4*S**(-1)*MS2*(TG+UG)**(-1)*S4 + 4*
     +    S**(-1)*MS2 + 4*S**(-1)*(TG+UG)**(-1)*S4**2 - 2*S**(-1)*S4 + 
     +    4*S*TG**(-1)*M2*(S+UG)**(-1) - 8*S*TG**(-1)*(S+UG)**(-1)*S4
     +     + 4*S*TG**(-1) + 2*S*UG**(-1) + 4*S*(S+UG)**(-1) + 4*S*
     +    (TG+UG)**(-1) + 16*TG**(-2)*M2*(S+UG)**(-1)*S4**2 - 16*
     +    TG**(-2)*M2*S4 + 16*TG**(-2)*MS2*(S+UG)**(-1)*S4**2 - 16*
     +    TG**(-2)*MS2*S4 + 8*TG**(-1)*M2*MS2*(S+UG)**(-1) - 16*
     +    TG**(-1)*M2*(S+UG)**(-1)*S4 + 6*TG**(-1)*M2 + 4*TG**(-1)*
     +    M2**2*(S+UG)**(-1) - 24*TG**(-1)*MS2*(S+UG)**(-1)*S4 + 8*
     +    TG**(-1)*MS2 + 12*TG**(-1)*(S+UG)**(-1)*S4**2 - 14*TG**(-1)*
     +    S4 + 4*TG*(S+UG)**(-1) - 2*TG*(TG+UG)**(-1) + 4*UG**(-1)*M2
     +     - 4*UG**(-1)*S4 + 8*M2*(S+UG)**(-1) + 6*M2*(TG+UG)**(-1) + 8
     +    *MS2*(S+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(54)*N*CK*TWO**(-1) * ( 4*MS2*(TG+UG)**(-1)
     +     - 16*(S+UG)**(-1)*S4 - 8*(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(55)*N*CK*TWO**(-1) * ( 32*S**(-2)*M2*MS2**2
     +     + 8*S**(-1)*TG**(-1)*M2*MS2*S4 + 16*S**(-1)*TG**(-1)*M2*
     +    MS2**2 + 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) + 4*S**(-1)*TG*M2*
     +    MS2*(S+UG)**(-1) - 4*S**(-1)*TG*M2*(S+TG)**(-1)*S4 - 8*
     +    S**(-1)*TG*M2*(S+UG)**(-1)*S4 + 4*S**(-1)*TG*M2 + 2*S**(-1)*
     +    TG*M2**2*(S+UG)**(-1) + 4*S**(-1)*TG*MS2*(S+TG)**(-1)*S4 - 12
     +    *S**(-1)*TG*MS2*(S+UG)**(-1)*S4 + 10*S**(-1)*TG*(S+UG)**(-1)*
     +    S4**2 - 4*S**(-1)*TG*S4 - 2*S**(-1)*TG**2*M2*(S+TG)**(-1) + 4
     +    *S**(-1)*TG**2*M2*(S+UG)**(-1) + 4*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) + 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) - 2*S**(-1)*
     +    TG**2*(S+TG)**(-1)*S4 - 8*S**(-1)*TG**2*(S+UG)**(-1)*S4 + 4*
     +    S**(-1)*TG**2 - 2*S**(-1)*TG**3*(S+TG)**(-1) + 2*S**(-1)*
     +    TG**3*(S+UG)**(-1) + 8*S**(-1)*UG**(-1)*M2*MS2*S4 + 16*
     +    S**(-1)*UG**(-1)*M2*MS2**2 + 8*S**(-1)*M2*MS2*(S+TG)**(-1)*S4
     +     - 4*S**(-1)*M2*MS2 + 4*S**(-1)*M2*(S+UG)**(-1)*S4**2 + 4*
     +    S**(-1)*M2*S4 )
     +
      M2GGH = M2GGH + ANG2(55)*N*CK*TWO**(-1) * ( 8*S**(-1)*MS2*
     +    (S+UG)**(-1)*S4**2 - 12*S**(-1)*MS2*S4 - 4*S**(-1)*
     +    (S+UG)**(-1)*S4**3 + 4*S**(-1)*S4**2 + 2*S*TG**(-1)*M2*MS2*
     +    (S+UG)**(-1) - 2*S*TG**(-1)*M2*(S+UG)**(-1)*S4 - 4*S*TG**(-1)
     +    *MS2*(S+UG)**(-1)*S4 - 2*S*TG**(-1)*MS2 - 2*S*TG**(-1)*S4 - 2
     +    *S*TG*(S+TG)**(-1) + 2*S*UG**(-1)*M2*MS2*(S+TG)**(-1) - 2*S*
     +    UG**(-1)*M2*(S+TG)**(-1)*S4 - 4*S*UG**(-1)*MS2*(S+TG)**(-1)*
     +    S4 - 2*S*UG**(-1)*MS2 - 2*S*UG**(-1)*S4 - 2*S*M2*(S+TG)**(-1)
     +     + 2*S*MS2*(S+UG)**(-1) - 2*S*(S+UG)**(-1)*S4 + 4*TG**(-1)*M2
     +    *(S+UG)**(-1)*S4**2 + 8*TG**(-1)*MS2*(S+UG)**(-1)*S4**2 - 8*
     +    TG**(-1)*MS2*S4 - 4*TG*M2*(S+TG)**(-1) + 4*TG*M2*(S+UG)**(-1)
     +     + 2*TG*MS2*(S+TG)**(-1) + 6*TG*MS2*(S+UG)**(-1) - 8*TG*
     +    (S+UG)**(-1)*S4 + 4*TG - 4*TG**2*(S+TG)**(-1) + 2*TG**2*
     +    (S+UG)**(-1) + 4*UG**(-1)*M2*(S+TG)**(-1)*S4**2 + 8*UG**(-1)*
     +    MS2*(S+TG)**(-1)*S4**2 - 8*UG**(-1)*MS2*S4 + 4*M2*MS2*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(55)*N*CK*TWO**(-1) * ( 6*M2*MS2*(S+UG)**(-1)
     +     - 4*M2*(S+TG)**(-1)*S4 - 8*M2*(S+UG)**(-1)*S4 + 4*M2 + 2*
     +    M2**2*(S+UG)**(-1) - 4*MS2*(S+TG)**(-1)*S4 - 14*MS2*
     +    (S+UG)**(-1)*S4 + 6*(S+UG)**(-1)*S4**2 - 8*S4 )
     +
      M2GGH = M2GGH + ANG2(55)*CQED*TWO**(-1) * (  - 16*S**(-2)*M2*
     +    MS2**2 - 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) - 4*S**(-1)*TG*M2*
     +    MS2*(S+UG)**(-1) + 4*S**(-1)*TG*M2*(S+TG)**(-1)*S4 + 8*
     +    S**(-1)*TG*M2*(S+UG)**(-1)*S4 - 4*S**(-1)*TG*M2 - 2*S**(-1)*
     +    TG*M2**2*(S+UG)**(-1) - 4*S**(-1)*TG*MS2*(S+TG)**(-1)*S4 + 12
     +    *S**(-1)*TG*MS2*(S+UG)**(-1)*S4 - 10*S**(-1)*TG*(S+UG)**(-1)*
     +    S4**2 + 4*S**(-1)*TG*S4 + 2*S**(-1)*TG**2*M2*(S+TG)**(-1) - 4
     +    *S**(-1)*TG**2*M2*(S+UG)**(-1) - 4*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) - 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) + 2*S**(-1)*
     +    TG**2*(S+TG)**(-1)*S4 + 8*S**(-1)*TG**2*(S+UG)**(-1)*S4 - 4*
     +    S**(-1)*TG**2 + 2*S**(-1)*TG**3*(S+TG)**(-1) - 2*S**(-1)*
     +    TG**3*(S+UG)**(-1) - 8*S**(-1)*M2*MS2*(S+TG)**(-1)*S4 + 4*
     +    S**(-1)*M2*MS2 - 4*S**(-1)*M2*(S+UG)**(-1)*S4**2 - 4*S**(-1)*
     +    M2*S4 - 8*S**(-1)*MS2*(S+UG)**(-1)*S4**2 + 12*S**(-1)*MS2*S4
     +     + 4*S**(-1)*(S+UG)**(-1)*S4**3 - 4*S**(-1)*S4**2 + 6*S*TG*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(55)*CQED*TWO**(-1) * ( 2*S*M2*(S+TG)**(-1)*
     +    (S+UG)**(-1)*S4 + 2*S*M2*(S+TG)**(-1) + 2*S*MS2*(S+TG)**(-1)*
     +    (S+UG)**(-1)*S4 - 2*S*MS2*(S+TG)**(-1) - 4*S*(S+TG)**(-1)*
     +    (S+UG)**(-1)*S4**2 + 4*S*(S+TG)**(-1)*S4 + 2*S*(S+UG)**(-1)*
     +    S4 - 2*S - 2*S**2*(S+TG)**(-1)*(S+UG)**(-1)*S4 + 2*S**2*
     +    (S+TG)**(-1) + 4*TG*M2*(S+TG)**(-1) - 6*TG*MS2*(S+TG)**(-1)
     +     - 2*TG*MS2*(S+UG)**(-1) + 4*TG*(S+TG)**(-1)*S4 - 4*TG + 6*
     +    TG**2*(S+TG)**(-1) - 4*M2*MS2*(S+TG)**(-1) - 2*M2*MS2*
     +    (S+UG)**(-1) + 4*M2*(S+TG)**(-1)*S4 - 2*M2 - 4*MS2*
     +    (S+TG)**(-1)*S4 + 2*MS2*(S+UG)**(-1)*S4 + 4*MS2 - 4*
     +    (S+TG)**(-1)*(S+UG)**(-1)*S4**3 + 4*(S+TG)**(-1)*S4**2 + 4*
     +    (S+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(56)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 8*
     +    S**(-1)*TG*M2*MS2 + 2*S*TG - 2*S*UG**(-1)*M2*MS2 - 2*S*
     +    UG**(-1)*M2**2 - 6*S*MS2 + 2*S**2 + 4*TG*M2 - 4*TG*MS2 - 16*
     +    UG**(-1)*M2*MS2**2 - 16*UG**(-1)*M2**2*MS2 - 4*UG**(-1)*M2**3
     +     - 16*M2*MS2 - 8*M2**2 )
     +
      M2GGH = M2GGH + ANG2(56)*N*CO*TWO**(-1) * ( 8*S**(-1)*UG**(-1)*M2
     +    *MS2*S4 - 8*S**(-1)*UG**(-1)*M2**2*MS2 - 2*S - 16*UG**(-1)*M2
     +    *MS2 - 4*UG**(-1)*M2**2 - 8*M2 )
     +
      M2GGH = M2GGH + ANG2(56)*N*CK*TWO**(-1) * (  - 32*S**(-2)*M2*
     +    MS2**2 - 8*S**(-1)*TG**(-1)*M2*MS2*S4 - 16*S**(-1)*TG**(-1)*
     +    M2*MS2**2 - 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) + 4*S**(-1)*TG*
     +    M2*(S+TG)**(-1)*S4 + 8*S**(-1)*TG*M2*(TG+UG)**(-1)*S4 - 2*
     +    S**(-1)*TG*M2 - 4*S**(-1)*TG*MS2*(S+TG)**(-1)*S4 + 8*S**(-1)*
     +    TG*MS2*(TG+UG)**(-1)*S4 - 4*S**(-1)*TG*MS2 - 8*S**(-1)*TG*
     +    (TG+UG)**(-1)*S4**2 + 6*S**(-1)*TG*S4 + 2*S**(-1)*TG**2*M2*
     +    (S+TG)**(-1) - 4*S**(-1)*TG**2*MS2*(S+TG)**(-1) + 2*S**(-1)*
     +    TG**2*(S+TG)**(-1)*S4 - 2*S**(-1)*TG**2 + 2*S**(-1)*TG**3*
     +    (S+TG)**(-1) - 24*S**(-1)*UG**(-1)*M2*MS2*S4 - 16*S**(-1)*
     +    UG**(-1)*M2*MS2**2 + 8*S**(-1)*UG**(-1)*M2**2*MS2 - 8*S**(-1)
     +    *M2*MS2*(S+TG)**(-1)*S4 + 4*S**(-1)*M2*MS2*(TG+UG)**(-1)*S4
     +     + 4*S**(-1)*M2*MS2 - 8*S**(-1)*M2*(TG+UG)**(-1)*S4**2 + 4*
     +    S**(-1)*M2**2*(TG+UG)**(-1)*S4 - 4*S**(-1)*MS2*(TG+UG)**(-1)*
     +    S4**2 + 8*S**(-1)*MS2*S4 + 4*S**(-1)*(TG+UG)**(-1)*S4**3 - 4*
     +    S**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(56)*N*CK*TWO**(-1) * (  - 2*S*TG**(-1)*M2 + 
     +    6*S*TG**(-1)*MS2 + 2*S*TG**(-1)*S4 + 10*S*TG*(S+TG)**(-1) - 2
     +    *S*TG*(TG+UG)**(-1) + 4*S*UG**(-1)*M2*MS2*(TG+UG)**(-1) + 2*S
     +    *UG**(-1)*M2*(S+TG)**(-1)*S4 - 14*S*UG**(-1)*M2*(TG+UG)**(-1)
     +    *S4 + 4*S*UG**(-1)*M2 + 4*S*UG**(-1)*M2**2*(TG+UG)**(-1) + 2*
     +    S*UG**(-1)*MS2*(S+TG)**(-1)*S4 - 6*S*UG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4 + 2*S*M2*(S+TG)**(-1) - 8*S*M2*(TG+UG)**(-1)
     +     - 6*S*MS2*(S+TG)**(-1) - 6*S*MS2*(TG+UG)**(-1) + 8*S*
     +    (S+TG)**(-1)*S4 + 4*S*(TG+UG)**(-1)*S4 - 6*S - 2*S**2*
     +    TG**(-1) + 6*S**2*UG**(-1)*M2*(TG+UG)**(-1) + 2*S**2*UG**(-1)
     +    *MS2*(TG+UG)**(-1) + 4*S**2*(S+TG)**(-1) + 2*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 + 6*TG**(-1)*UG**(-1)*M2**2*MS2 + 2*
     +    TG**(-1)*UG**(-1)*M2**2*S4 + 2*TG**(-1)*UG**(-1)*M2**3 + 14*
     +    TG**(-1)*M2*MS2 + 4*TG**(-1)*M2*S4 + 2*TG**(-1)*M2**2 + 4*TG*
     +    M2*(S+TG)**(-1) - 4*TG*M2*(TG+UG)**(-1) - 10*TG*MS2*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(56)*N*CK*TWO**(-1) * (  - 4*TG*MS2*
     +    (TG+UG)**(-1) + 8*TG*(S+TG)**(-1)*S4 + 8*TG*(TG+UG)**(-1)*S4
     +     - 6*TG + 8*TG**2*(S+TG)**(-1) + 16*UG**(-2)*M2**2*MS2 - 4*
     +    UG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 6*UG**(-1)*M2*MS2 + 8*
     +    UG**(-1)*M2*(S+TG)**(-1)*S4**2 + 8*UG**(-1)*M2*(TG+UG)**(-1)*
     +    S4**2 - 12*UG**(-1)*M2*S4 - 4*UG**(-1)*M2**2*(TG+UG)**(-1)*S4
     +     + 2*UG**(-1)*M2**2 + 8*UG**(-1)*MS2*(S+TG)**(-1)*S4**2 + 4*
     +    UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 - 12*UG**(-1)*MS2*S4 - 4*
     +    UG**(-1)*(S+TG)**(-1)*S4**3 + 4*UG**(-1)*S4**2 - 4*M2*MS2*
     +    (S+TG)**(-1) - 4*M2*MS2*(TG+UG)**(-1) + 4*M2*(S+TG)**(-1)*S4
     +     + 20*M2*(TG+UG)**(-1)*S4 - 6*M2 - 4*M2**2*(TG+UG)**(-1) - 12
     +    *MS2*(S+TG)**(-1)*S4 + 16*MS2*(TG+UG)**(-1)*S4 + 2*MS2 + 4*
     +    (S+TG)**(-1)*S4**2 - 8*(TG+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(56)*CQED*TWO**(-1) * ( 16*S**(-2)*M2*MS2**2
     +     + 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) - 4*S**(-1)*TG*M2*
     +    (S+TG)**(-1)*S4 + 2*S**(-1)*TG*M2 + 4*S**(-1)*TG*MS2*
     +    (S+TG)**(-1)*S4 - 4*S**(-1)*TG*MS2 + 2*S**(-1)*TG*S4 - 2*
     +    S**(-1)*TG**2*M2*(S+TG)**(-1) + 4*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) - 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 + 2*S**(-1)*
     +    TG**2 - 2*S**(-1)*TG**3*(S+TG)**(-1) + 8*S**(-1)*UG**(-1)*M2*
     +    MS2*S4 + 8*S**(-1)*M2*MS2*(S+TG)**(-1)*S4 - 4*S**(-1)*M2*MS2
     +     + 4*S**(-1)*M2*S4 - 4*S**(-1)*MS2*S4 - 2*S*TG*M2*
     +    (S+TG)**(-1)*(TG+UG)**(-1) - 2*S*TG*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) + 2*S*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4 - 8*S*TG
     +    *(S+TG)**(-1) - 2*S*TG**2*(S+TG)**(-1)*(TG+UG)**(-1) - 4*S*
     +    UG**(-1)*M2*MS2*(TG+UG)**(-1) + 14*S*UG**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 6*S*UG**(-1)*M2 - 4*S*UG**(-1)*M2**2*
     +    (TG+UG)**(-1) + 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 - 2*S*
     +    UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(56)*CQED*TWO**(-1) * (  - 2*S*M2*
     +    (S+TG)**(-1)*(TG+UG)**(-1)*S4 - 2*S*M2*(S+TG)**(-1) - 2*S*M2*
     +    (TG+UG)**(-1) - 2*S*MS2*(S+TG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*
     +    MS2*(S+TG)**(-1) + 4*S*MS2*(TG+UG)**(-1) + 8*S*(S+TG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 - 4*S*(S+TG)**(-1)*S4 + 4*S*(TG+UG)**(-1)
     +    *S4 - 6*S**2*TG*(S+TG)**(-1)*(TG+UG)**(-1) - 6*S**2*UG**(-1)*
     +    M2*(TG+UG)**(-1) - 2*S**2*UG**(-1)*MS2*(TG+UG)**(-1) - 2*S**2
     +    *M2*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S**2*MS2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) - 4*S**2*(S+TG)**(-1) - 4*S**3*(S+TG)**(-1)*
     +    (TG+UG)**(-1) - 4*TG*M2*(S+TG)**(-1) - 2*TG*M2*(TG+UG)**(-1)
     +     + 6*TG*MS2*(S+TG)**(-1) + 4*TG*(S+TG)**(-1)*(TG+UG)**(-1)*
     +    S4**2 - 4*TG*(S+TG)**(-1)*S4 + 2*TG*(TG+UG)**(-1)*S4 + 2*TG
     +     - 6*TG**2*(S+TG)**(-1) - 8*UG**(-2)*M2**2*MS2 + 4*UG**(-1)*
     +    M2*MS2*(TG+UG)**(-1)*S4 - 8*UG**(-1)*M2*(TG+UG)**(-1)*S4**2
     +     + 8*UG**(-1)*M2*S4 + 4*UG**(-1)*M2**2*(TG+UG)**(-1)*S4 - 4*
     +    UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(56)*CQED*TWO**(-1) * ( 4*UG**(-1)*MS2*S4 + 4
     +    *M2*MS2*(S+TG)**(-1) - 2*M2*MS2*(TG+UG)**(-1) - 4*M2*
     +    (S+TG)**(-1)*S4 + 2*M2*(TG+UG)**(-1)*S4 - 2*M2**2*
     +    (TG+UG)**(-1) + 4*MS2*(S+TG)**(-1)*S4 - 6*MS2*(TG+UG)**(-1)*
     +    S4 - 4*(S+TG)**(-1)*(TG+UG)**(-1)*S4**3 + 4*(S+TG)**(-1)*
     +    S4**2 - 4*(TG+UG)**(-1)*S4**2 + 2*S4 )
     +
      M2GGH = M2GGH + ANG2(57)*N*CO*TWO**(-1) * (  - 16*UG**(-1)*M2*
     +    MS2**2 - 8*UG**(-1)*M2**2*MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(57)*N*CK*TWO**(-1) * ( 16*UG**(-1)*M2*MS2**2
     +     + 24*UG**(-1)*M2**2*MS2 + 24*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(57)*CQED*TWO**(-1) * (  - 8*UG**(-1)*M2**2*
     +    MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(58)*N*CK*TWO**(-1) * ( 16*S**(-1)*TG**(-1)*
     +    M2*MS2*S4 - 8*S**(-1)*TG**(-1)*M2**2*MS2 - 4*S**(-1)*TG*M2*
     +    MS2*(S+TG)**(-1) + 4*S**(-1)*TG*M2*(S+TG)**(-1)*S4 - 4*
     +    S**(-1)*TG*M2 - 4*S**(-1)*TG*MS2*(S+TG)**(-1)*S4 + 4*S**(-1)*
     +    TG*MS2 - 4*S**(-1)*TG*S4 + 2*S**(-1)*TG**2*M2*(S+TG)**(-1) - 
     +    4*S**(-1)*TG**2*MS2*(S+TG)**(-1) + 2*S**(-1)*TG**2*
     +    (S+TG)**(-1)*S4 + 2*S**(-1)*TG**3*(S+TG)**(-1) - 8*S**(-1)*M2
     +    *MS2*(S+TG)**(-1)*S4 + 4*S**(-1)*M2*MS2 + 4*S**(-1)*MS2*S4 - 
     +    4*S*TG**(-1)*M2*MS2*(TG+UG)**(-1) + 14*S*TG**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 6*S*TG**(-1)*M2 - 4*S*TG**(-1)*M2**2*
     +    (TG+UG)**(-1) + 6*S*TG**(-1)*MS2*(TG+UG)**(-1)*S4 - 2*S*
     +    TG**(-1)*MS2 + 2*S*TG*(S+TG)**(-1) + 9*S*TG*(TG+UG)**(-1) + 2
     +    *S*UG**(-1)*M2*(S+TG)**(-1)*S4 + 2*S*UG**(-1)*MS2*
     +    (S+TG)**(-1)*S4 + 4*S*M2*(S+TG)**(-1) + 3*S*M2*(TG+UG)**(-1)
     +     - 4*S*MS2*(S+TG)**(-1) + 4*S*MS2*(TG+UG)**(-1) - 14*S*
     +    (TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(58)*N*CK*TWO**(-1) * ( 7*S - 6*S**2*TG**(-1)
     +    *M2*(TG+UG)**(-1) - 2*S**2*TG**(-1)*MS2*(TG+UG)**(-1) + 7*
     +    S**2*(TG+UG)**(-1) - 16*TG**(-2)*M2**2*MS2 - 2*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 - 6*TG**(-1)*UG**(-1)*M2**2*MS2 - 2*
     +    TG**(-1)*UG**(-1)*M2**2*S4 - 2*TG**(-1)*UG**(-1)*M2**3 + 4*
     +    TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 2*TG**(-1)*M2*MS2 - 8*
     +    TG**(-1)*M2*(TG+UG)**(-1)*S4**2 + 8*TG**(-1)*M2*S4 + 4*
     +    TG**(-1)*M2**2*(TG+UG)**(-1)*S4 - 2*TG**(-1)*M2**2 - 4*
     +    TG**(-1)*MS2*(TG+UG)**(-1)*S4**2 + 4*TG**(-1)*MS2*S4 + 6*TG*
     +    M2*(S+TG)**(-1) + 3*TG*M2*(TG+UG)**(-1) - 8*TG*MS2*
     +    (S+TG)**(-1) - 9*TG*(TG+UG)**(-1)*S4 + 7*TG + 4*TG**2*
     +    (S+TG)**(-1) + 2*TG**2*(TG+UG)**(-1) + 8*UG**(-1)*M2*MS2*
     +    (S+TG)**(-1)*S4 - 2*UG**(-1)*M2*(S+TG)**(-1)*S4**2 + 2*
     +    UG**(-1)*M2*S4 - 2*UG**(-1)*M2**2*MS2*(S+TG)**(-1) + 8*
     +    UG**(-1)*M2**2*(S+TG)**(-1)*S4 - 2*UG**(-1)*M2**3*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(58)*N*CK*TWO**(-1) * ( 2*UG**(-1)*MS2*
     +    (S+TG)**(-1)*S4**2 - 2*UG**(-1)*MS2*S4 - 6*M2*MS2*
     +    (S+TG)**(-1) + 2*M2*MS2*(TG+UG)**(-1) + 10*M2*(S+TG)**(-1)*S4
     +     - M2*(TG+UG)**(-1)*S4 - 4*M2 - 2*M2**2*(S+TG)**(-1) + 2*
     +    M2**2*(TG+UG)**(-1) - 6*MS2*(S+TG)**(-1)*S4 - 2*MS2*
     +    (TG+UG)**(-1)*S4 + 8*MS2 + 7*(TG+UG)**(-1)*S4**2 - 7*S4 )
     +
      M2GGH = M2GGH + ANG2(58)*CQED*TWO**(-1) * (  - 8*S**(-1)*TG**(-1)
     +    *M2*MS2*S4 + 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) - 4*S**(-1)*TG*
     +    M2*(S+TG)**(-1)*S4 + 2*S**(-1)*TG*M2 + 4*S**(-1)*TG*MS2*
     +    (S+TG)**(-1)*S4 - 4*S**(-1)*TG*MS2 + 2*S**(-1)*TG*S4 - 2*
     +    S**(-1)*TG**2*M2*(S+TG)**(-1) + 4*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) - 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 + 2*S**(-1)*
     +    TG**2 - 2*S**(-1)*TG**3*(S+TG)**(-1) + 8*S**(-1)*M2*MS2*
     +    (S+TG)**(-1)*S4 - 4*S**(-1)*M2*MS2 + 4*S**(-1)*M2*S4 - 4*
     +    S**(-1)*MS2*S4 + 4*S*TG**(-1)*M2*MS2*(TG+UG)**(-1) - 14*S*
     +    TG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*TG**(-1)*M2 + 4*S*TG**(-1)
     +    *M2**2*(TG+UG)**(-1) - 6*S*TG**(-1)*MS2*(TG+UG)**(-1)*S4 + 2*
     +    S*TG**(-1)*MS2 - 2*S*TG*M2*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S*
     +    TG*MS2*(S+TG)**(-1)*(TG+UG)**(-1) + 2*S*TG*(S+TG)**(-1)*
     +    (TG+UG)**(-1)*S4 - 8*S*TG*(S+TG)**(-1) - 2*S*TG*(TG+UG)**(-1)
     +     - 2*S*TG**2*(S+TG)**(-1)*(TG+UG)**(-1) - 2*S*M2*(S+TG)**(-1)
     +    *(TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(58)*CQED*TWO**(-1) * (  - 2*S*M2*
     +    (S+TG)**(-1) - 2*S*MS2*(S+TG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*
     +    MS2*(S+TG)**(-1) + 8*S*(S+TG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*S
     +    *(S+TG)**(-1)*S4 + 4*S*(TG+UG)**(-1)*S4 + 6*S**2*TG**(-1)*M2*
     +    (TG+UG)**(-1) + 2*S**2*TG**(-1)*MS2*(TG+UG)**(-1) - 6*S**2*TG
     +    *(S+TG)**(-1)*(TG+UG)**(-1) - 2*S**2*M2*(S+TG)**(-1)*
     +    (TG+UG)**(-1) - 2*S**2*MS2*(S+TG)**(-1)*(TG+UG)**(-1) - 4*
     +    S**2*(S+TG)**(-1) - 4*S**3*(S+TG)**(-1)*(TG+UG)**(-1) + 8*
     +    TG**(-2)*M2**2*MS2 - 4*TG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 8*
     +    TG**(-1)*M2*(TG+UG)**(-1)*S4**2 - 8*TG**(-1)*M2*S4 - 4*
     +    TG**(-1)*M2**2*(TG+UG)**(-1)*S4 + 4*TG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4**2 - 4*TG**(-1)*MS2*S4 - 4*TG*M2*
     +    (S+TG)**(-1) - 2*TG*M2*(TG+UG)**(-1) + 6*TG*MS2*(S+TG)**(-1)
     +     + 4*TG*(S+TG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*TG*(S+TG)**(-1)*
     +    S4 + 2*TG*(TG+UG)**(-1)*S4 - 6*TG**2*(S+TG)**(-1) + 4*M2*MS2*
     +    (S+TG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(58)*CQED*TWO**(-1) * (  - 2*M2*MS2*
     +    (TG+UG)**(-1) - 4*M2*(S+TG)**(-1)*S4 + 2*M2*(TG+UG)**(-1)*S4
     +     + 2*M2 - 2*M2**2*(TG+UG)**(-1) + 4*MS2*(S+TG)**(-1)*S4 + 2*
     +    MS2*(TG+UG)**(-1)*S4 - 4*MS2 - 4*(S+TG)**(-1)*(TG+UG)**(-1)*
     +    S4**3 + 4*(S+TG)**(-1)*S4**2 - 4*(TG+UG)**(-1)*S4**2 + 4*S4 )
     +
      M2GGH = M2GGH + ANG2(59)*N*CK*TWO**(-1) * (  - 4*S**(-1)*TG*M2*
     +    MS2*(S+UG)**(-1) + 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 + 4*
     +    S**(-1)*TG*M2 - 2*S**(-1)*TG*M2**2*(S+UG)**(-1) + 12*S**(-1)*
     +    TG*MS2*(S+UG)**(-1)*S4 - 4*S**(-1)*TG*MS2 - 10*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**2 + 4*S**(-1)*TG*S4 - 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-1) - 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) + 8*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4 - 2*S**(-1)*TG**3*(S+UG)**(-1) + 16*
     +    S**(-1)*UG**(-1)*M2*MS2*S4 - 8*S**(-1)*UG**(-1)*M2**2*MS2 - 4
     +    *S**(-1)*M2*(S+UG)**(-1)*S4**2 + 4*S**(-1)*M2**2 - 8*S**(-1)*
     +    MS2*(S+UG)**(-1)*S4**2 + 8*S**(-1)*MS2*S4 + 4*S**(-1)*
     +    (S+UG)**(-1)*S4**3 - 4*S**(-1)*S4**2 + 2*S*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4 + 2*S*TG**(-1)*MS2*(S+UG)**(-1)*S4 - 5*S*TG*
     +    (TG+UG)**(-1) - 4*S*UG**(-1)*M2*MS2*(TG+UG)**(-1) + 14*S*
     +    UG**(-1)*M2*(TG+UG)**(-1)*S4 - 6*S*UG**(-1)*M2 - 4*S*UG**(-1)
     +    *M2**2*(TG+UG)**(-1) + 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 - 2*
     +    S*UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(59)*N*CK*TWO**(-1) * (  - 5*S*M2*
     +    (TG+UG)**(-1) + 4*S*MS2*(TG+UG)**(-1) + 2*S*(S+UG)**(-1)*S4
     +     - 6*S**2*UG**(-1)*M2*(TG+UG)**(-1) - 2*S**2*UG**(-1)*MS2*
     +    (TG+UG)**(-1) - 2*TG**(-1)*UG**(-1)*M2*MS2*S4 - 6*TG**(-1)*
     +    UG**(-1)*M2**2*MS2 - 2*TG**(-1)*UG**(-1)*M2**2*S4 - 2*
     +    TG**(-1)*UG**(-1)*M2**3 + 8*TG**(-1)*M2*MS2*(S+UG)**(-1)*S4
     +     - 2*TG**(-1)*M2*(S+UG)**(-1)*S4**2 + 2*TG**(-1)*M2*S4 - 2*
     +    TG**(-1)*M2**2*MS2*(S+UG)**(-1) + 8*TG**(-1)*M2**2*
     +    (S+UG)**(-1)*S4 - 2*TG**(-1)*M2**3*(S+UG)**(-1) + 2*TG**(-1)*
     +    MS2*(S+UG)**(-1)*S4**2 - 2*TG**(-1)*MS2*S4 - 6*TG*M2*
     +    (S+UG)**(-1) + TG*M2*(TG+UG)**(-1) + 8*TG*(S+UG)**(-1)*S4 + 5
     +    *TG*(TG+UG)**(-1)*S4 - 7*TG - 2*TG**2*(S+UG)**(-1) + 2*TG**2*
     +    (TG+UG)**(-1) - 16*UG**(-2)*M2**2*MS2 + 4*UG**(-1)*M2*MS2*
     +    (TG+UG)**(-1)*S4 + 2*UG**(-1)*M2*MS2 - 8*UG**(-1)*M2*
     +    (TG+UG)**(-1)*S4**2 + 8*UG**(-1)*M2*S4 + 4*UG**(-1)*M2**2*
     +    (TG+UG)**(-1)*S4 )
     +
      M2GGH = M2GGH + ANG2(59)*N*CK*TWO**(-1) * (  - 2*UG**(-1)*M2**2
     +     - 4*UG**(-1)*MS2*(TG+UG)**(-1)*S4**2 + 4*UG**(-1)*MS2*S4 - 2
     +    *M2*MS2*(S+UG)**(-1) + 2*M2*MS2*(TG+UG)**(-1) + 16*M2*
     +    (S+UG)**(-1)*S4 + 7*M2*(TG+UG)**(-1)*S4 - 7*M2 - 6*M2**2*
     +    (S+UG)**(-1) + M2**2*(TG+UG)**(-1) - 2*MS2*(S+UG)**(-1)*S4 - 
     +    2*MS2*(TG+UG)**(-1)*S4 + 4*MS2 - 6*(S+UG)**(-1)*S4**2 + 4*S4
     +     )
     +
      M2GGH = M2GGH + ANG2(59)*CQED*TWO**(-1) * ( 4*S**(-1)*TG*M2*MS2*
     +    (S+UG)**(-1) - 8*S**(-1)*TG*M2*(S+UG)**(-1)*S4 + 2*S**(-1)*TG
     +    *M2 + 2*S**(-1)*TG*M2**2*(S+UG)**(-1) - 12*S**(-1)*TG*MS2*
     +    (S+UG)**(-1)*S4 + 4*S**(-1)*TG*MS2 + 10*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**2 - 6*S**(-1)*TG*S4 + 4*S**(-1)*TG**2*M2*
     +    (S+UG)**(-1) + 4*S**(-1)*TG**2*MS2*(S+UG)**(-1) - 8*S**(-1)*
     +    TG**2*(S+UG)**(-1)*S4 + 2*S**(-1)*TG**2 + 2*S**(-1)*TG**3*
     +    (S+UG)**(-1) - 8*S**(-1)*UG**(-1)*M2*MS2*S4 + 4*S**(-1)*M2*
     +    (S+UG)**(-1)*S4**2 + 8*S**(-1)*MS2*(S+UG)**(-1)*S4**2 - 8*
     +    S**(-1)*MS2*S4 - 4*S**(-1)*(S+UG)**(-1)*S4**3 + 4*S**(-1)*
     +    S4**2 - 2*S*TG*M2*(S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*MS2*
     +    (S+UG)**(-1)*(TG+UG)**(-1) + 2*S*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S*TG*(S+UG)**(-1) + 2*S*TG*(TG+UG)**(-1)
     +     - 2*S*TG**2*(S+UG)**(-1)*(TG+UG)**(-1) + 4*S*UG**(-1)*M2*MS2
     +    *(TG+UG)**(-1) - 14*S*UG**(-1)*M2*(TG+UG)**(-1)*S4 + 6*S*
     +    UG**(-1)*M2 )
     +
      M2GGH = M2GGH + ANG2(59)*CQED*TWO**(-1) * ( 4*S*UG**(-1)*M2**2*
     +    (TG+UG)**(-1) - 6*S*UG**(-1)*MS2*(TG+UG)**(-1)*S4 + 2*S*
     +    UG**(-1)*MS2 + 2*S*M2*MS2*(S+UG)**(-1)*(TG+UG)**(-1) - 2*S*M2
     +    *(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 2*S*M2*(S+UG)**(-1) + 4*S*M2
     +    *(TG+UG)**(-1) - 4*S*MS2*(S+UG)**(-1)*(TG+UG)**(-1)*S4 + 4*S*
     +    (S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 4*S*(S+UG)**(-1)*S4 + 2*S
     +     + 2*S**2*TG*(S+UG)**(-1)*(TG+UG)**(-1) + 6*S**2*UG**(-1)*M2*
     +    (TG+UG)**(-1) + 2*S**2*UG**(-1)*MS2*(TG+UG)**(-1) + 2*S**2*M2
     +    *(S+UG)**(-1)*(TG+UG)**(-1) - 4*S**2*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4 + 2*S**2*(TG+UG)**(-1) + 2*TG*M2*
     +    (TG+UG)**(-1) + 2*TG*MS2*(S+UG)**(-1) - 4*TG*(S+UG)**(-1)*
     +    (TG+UG)**(-1)*S4**2 - 2*TG*(TG+UG)**(-1)*S4 + 4*TG + 8*
     +    UG**(-2)*M2**2*MS2 - 4*UG**(-1)*M2*MS2*(TG+UG)**(-1)*S4 + 8*
     +    UG**(-1)*M2*(TG+UG)**(-1)*S4**2 - 8*UG**(-1)*M2*S4 - 4*
     +    UG**(-1)*M2**2*(TG+UG)**(-1)*S4 + 4*UG**(-1)*MS2*
     +    (TG+UG)**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(59)*CQED*TWO**(-1) * (  - 4*UG**(-1)*MS2*S4
     +     + 2*M2*MS2*(S+UG)**(-1) - 2*M2*MS2*(TG+UG)**(-1) - 4*M2*
     +    (S+UG)**(-1)*(TG+UG)**(-1)*S4**2 - 2*M2*(TG+UG)**(-1)*S4 + 4*
     +    M2 - 2*MS2*(S+UG)**(-1)*S4 + 2*MS2*(TG+UG)**(-1)*S4 + 4*
     +    (S+UG)**(-1)*S4**2 - 2*(TG+UG)**(-1)*S4**2 - 2*S4 )
     +
      M2GGH = M2GGH + ANG2(60)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 2*
     +    S**(-1) - 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(60)*N*CO*TWO**(-1) * ( 2*S**(-1)*TG**(-1) - 
     +    2*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(61)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 2*
     +    S**(-1) - 2*TG**(-1) )
     +
      M2GGH = M2GGH + ANG2(61)*N*CO*TWO**(-1) * ( 2*S**(-1)*UG**(-1) - 
     +    2*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(62)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 2*S*
     +    TG**(-1)*M2*MS2 + 2*S*TG**(-1)*M2**2 + 2*S*UG**(-1)*M2*MS2 + 
     +    2*S*UG**(-1)*M2**2 + 2*S*M2 + 2*S*MS2 - 8*TG**(-1)*M2**2*MS2
     +     - 4*TG**(-1)*M2**3 - 8*UG**(-1)*M2**2*MS2 - 4*UG**(-1)*M2**3
     +     - 8*M2*MS2 - 4*M2**2 )
     +
      M2GGH = M2GGH + ANG2(62)*N*CO*TWO**(-1) * ( 2*S*UG**(-1)*M2 + 4*S
     +    *UG**(-1)*MS2 - 2*TG**(-1)*UG**(-1)*M2*MS2*S4 - 6*TG**(-1)*
     +    UG**(-1)*M2**2*MS2 - 2*TG**(-1)*UG**(-1)*M2**2*S4 - 2*
     +    TG**(-1)*UG**(-1)*M2**3 + 10*TG**(-1)*M2*MS2 - 6*TG**(-1)*
     +    M2**2 - 2*TG - 6*UG**(-1)*M2*MS2 - 4*UG**(-1)*M2*S4 - 2*
     +    UG**(-1)*M2**2 - 8*UG**(-1)*MS2*S4 - 4*M2 + 8*MS2 )
     +
      M2GGH = M2GGH + ANG2(62)*N*CK*TWO**(-1) * (  - 16*S**(-1)*
     +    TG**(-1)*M2*MS2*S4 + 8*S**(-1)*TG**(-1)*M2**2*MS2 + 4*S**(-1)
     +    *TG*M2*MS2*(S+TG)**(-1) - 4*S**(-1)*TG*M2*(S+TG)**(-1)*S4 + 4
     +    *S**(-1)*TG*MS2*(S+TG)**(-1)*S4 - 4*S**(-1)*TG*MS2 + 4*
     +    S**(-1)*TG*S4 - 2*S**(-1)*TG**2*M2*(S+TG)**(-1) + 4*S**(-1)*
     +    TG**2*MS2*(S+TG)**(-1) - 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 + 2*
     +    S**(-1)*TG**2 - 2*S**(-1)*TG**3*(S+TG)**(-1) + 8*S**(-1)*M2*
     +    MS2*(S+TG)**(-1)*S4 + 4*S**(-1)*M2*MS2 + 8*S**(-1)*M2*S4 - 4*
     +    S**(-1)*MS2*S4 - 2*S*TG**(-1)*M2*(S+UG)**(-1)*S4 + 2*S*
     +    TG**(-1)*M2 - 2*S*TG**(-1)*MS2*(S+UG)**(-1)*S4 + 2*S*TG**(-1)
     +    *MS2 - 2*S*TG*(S+TG)**(-1) + 2*S*UG**(-1)*M2*MS2*(S+TG)**(-1)
     +     - 2*S*UG**(-1)*M2*(S+TG)**(-1)*S4 - 4*S*UG**(-1)*MS2*
     +    (S+TG)**(-1)*S4 + 8*S*UG**(-1)*MS2 - 2*S*M2*(S+TG)**(-1) + 16
     +    *TG**(-2)*M2**2*MS2 + 6*TG**(-1)*UG**(-1)*M2*MS2*S4 + 2*
     +    TG**(-1)*UG**(-1)*M2**2*MS2 + 6*TG**(-1)*UG**(-1)*M2**2*S4 + 
     +    6*TG**(-1)*UG**(-1)*M2**3 )
     +
      M2GGH = M2GGH + ANG2(62)*N*CK*TWO**(-1) * (  - 14*TG**(-1)*M2*MS2
     +     - 4*TG**(-1)*M2*S4 + 6*TG**(-1)*M2**2 + 4*TG**(-1)*
     +    (S+UG)**(-1)*S4**3 - 4*TG**(-1)*S4**2 - 4*TG*M2*(S+TG)**(-1)
     +     + 2*TG*MS2*(S+TG)**(-1) + 2*TG*(S+UG)**(-1)*S4 + 2*TG - 4*
     +    TG**2*(S+TG)**(-1) - 6*UG**(-1)*M2*MS2 + 4*UG**(-1)*M2*
     +    (S+TG)**(-1)*S4**2 + 4*UG**(-1)*M2*S4 - 2*UG**(-1)*M2**2 + 8*
     +    UG**(-1)*MS2*(S+TG)**(-1)*S4**2 + 4*M2*MS2*(S+TG)**(-1) - 4*
     +    M2*(S+TG)**(-1)*S4 - 4*MS2*(S+TG)**(-1)*S4 - 8*MS2 - 4*
     +    (S+UG)**(-1)*S4**2 )
     +
      M2GGH = M2GGH + ANG2(62)*CQED*TWO**(-1) * ( 8*S**(-1)*TG**(-1)*M2
     +    *MS2*S4 - 4*S**(-1)*TG*M2*MS2*(S+TG)**(-1) + 4*S**(-1)*TG*M2*
     +    (S+TG)**(-1)*S4 - 2*S**(-1)*TG*M2 - 4*S**(-1)*TG*MS2*
     +    (S+TG)**(-1)*S4 + 4*S**(-1)*TG*MS2 - 2*S**(-1)*TG*S4 + 2*
     +    S**(-1)*TG**2*M2*(S+TG)**(-1) - 4*S**(-1)*TG**2*MS2*
     +    (S+TG)**(-1) + 2*S**(-1)*TG**2*(S+TG)**(-1)*S4 - 2*S**(-1)*
     +    TG**2 + 2*S**(-1)*TG**3*(S+TG)**(-1) - 8*S**(-1)*M2*MS2*
     +    (S+TG)**(-1)*S4 + 4*S**(-1)*M2*MS2 - 4*S**(-1)*M2*S4 + 4*
     +    S**(-1)*MS2*S4 + 6*S*TG*(S+TG)**(-1) - 4*S*UG**(-1)*MS2 + 2*S
     +    *M2*(S+TG)**(-1)*(S+UG)**(-1)*S4 + 2*S*M2*(S+TG)**(-1) + 2*S*
     +    MS2*(S+TG)**(-1)*(S+UG)**(-1)*S4 - 2*S*MS2*(S+TG)**(-1) - 4*S
     +    *(S+TG)**(-1)*(S+UG)**(-1)*S4**2 + 4*S*(S+TG)**(-1)*S4 + 2*S*
     +    (S+UG)**(-1)*S4 - 2*S - 2*S**2*(S+TG)**(-1)*(S+UG)**(-1)*S4
     +     + 2*S**2*(S+TG)**(-1) - 8*TG**(-2)*M2**2*MS2 - 2*TG**(-1)*
     +    UG**(-1)*M2*MS2*S4 + 2*TG**(-1)*UG**(-1)*M2**2*MS2 - 2*
     +    TG**(-1)*UG**(-1)*M2**2*S4 )
     +
      M2GGH = M2GGH + ANG2(62)*CQED*TWO**(-1) * (  - 2*TG**(-1)*
     +    UG**(-1)*M2**3 + 2*TG**(-1)*M2*MS2 - 2*TG**(-1)*M2**2 + 4*TG*
     +    M2*(S+TG)**(-1) - 6*TG*MS2*(S+TG)**(-1) + 4*TG*(S+TG)**(-1)*
     +    S4 - 2*TG*(S+UG)**(-1)*S4 - 4*TG + 6*TG**2*(S+TG)**(-1) + 2*
     +    UG**(-1)*M2*MS2 - 2*UG**(-1)*M2*S4 - 4*M2*MS2*(S+TG)**(-1) + 
     +    4*M2*(S+TG)**(-1)*S4 - 2*M2 - 4*MS2*(S+TG)**(-1)*S4 + 2*MS2
     +     - 4*(S+TG)**(-1)*(S+UG)**(-1)*S4**3 + 4*(S+TG)**(-1)*S4**2
     +     + 4*(S+UG)**(-1)*S4**2 - 4*S4 )
     +
      M2GGH = M2GGH + ANG2(63)*N*CO*TWO**(-1) * (  - 8*TG**(-1)*M2**2*
     +    MS2 - 4*TG*MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(63)*N*CK*TWO**(-1) * ( 24*TG**(-1)*M2**2*MS2
     +     + 12*TG*MS2 + 24*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(63)*CQED*TWO**(-1) * (  - 8*TG**(-1)*M2**2*
     +    MS2 - 4*TG*MS2 - 8*M2*MS2 )
     +
      M2GGH = M2GGH + ANG2(64)*N*CO*TWO**(-1) * ( 8*UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(64)*N*CK*TWO**(-1) * (  - 8*UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(65)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 2*
     +    S**(-1) + 2*TG**(-1) )
     +
      M2GGH = M2GGH + ANG2(65)*N*CO*TWO**(-1) * (  - 2*S**(-1)*UG**(-1)
     +     + 2*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(66)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4*
     +    S**(-2)*TG - 2*S**(-1)*TG**(-1)*M2 - 2*S**(-1)*TG**(-1)*MS2
     +     + 2*S**(-1)*UG**(-1)*M2 + 2*S**(-1)*UG**(-1)*MS2 + 8*S**(-1)
     +     + 3*TG**(-1) )
     +
      M2GGH = M2GGH + ANG2(66)*N*CO*TWO**(-1) * (  - 2*S**(-1)*TG**(-1)
     +     - S**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(66)*N*CK*TWO**(-1) * ( 4*S**(-2)*TG*
     +    (TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(67)*N*CO*TWO**(-1) * ( 8*TG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(67)*N*CK*TWO**(-1) * (  - 8*TG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(68)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 8*
     +    S**(-1)*TG**(-1)*M2*MS2 + 8*S**(-1)*TG**(-1)*M2**2 - 2*
     +    S**(-1)*TG + 4*S**(-1)*M2 - 8*TG**(-1)*UG**(-1)*M2*MS2 - 8*
     +    TG**(-1)*UG**(-1)*M2**2 + 4*TG**(-1)*M2 + 8*TG**(-1)*MS2 - 2*
     +    UG**(-1)*M2 + 4*UG**(-1)*MS2 )
     +
      M2GGH = M2GGH + ANG2(68)*N*CO*TWO**(-1) * ( 6*S**(-1)*TG**(-1)*M2
     +     - 2*S**(-1)*TG**(-1)*S4 - 6*TG**(-1)*UG**(-1)*M2 + 4*
     +    TG**(-1)*UG**(-1)*MS2 + 4*TG**(-1)*UG**(-1)*S4 - 4*TG**(-1)
     +     - 4*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(68)*N*CK*TWO**(-1) * ( 4*TG**(-1)*M2*
     +    (S+UG)**(-1) - 12*TG**(-1)*(S+UG)**(-1)*S4 + 8*TG**(-1) + 4*
     +    (S+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(69)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 4
     +     + 8*S**(-2)*TG*M2 - 4*S**(-1)*TG**(-1)*M2*MS2 - 4*S**(-1)*
     +    TG**(-1)*M2**2 + 4*S**(-1)*UG**(-1)*M2*MS2 + 4*S**(-1)*
     +    UG**(-1)*M2**2 + 4*S**(-1)*M2 + 24*S**(-1)*MS2 - S*TG**(-1)
     +     - TG**(-1)*M2 + 4*TG**(-1)*MS2 - 3*UG**(-1)*M2 + 8*UG**(-1)*
     +    MS2 )
     +
      M2GGH = M2GGH + ANG2(69)*N*CO*TWO**(-1) * (  - 6*S**(-1)*TG**(-1)
     +    *MS2 - 2*S**(-1)*TG**(-1)*S4 + 6*S**(-1)*UG**(-1)*M2 - 6*
     +    S**(-1)*UG**(-1)*MS2 - 4*S**(-1)*UG**(-1)*S4 + 12*S**(-1) + 
     +    TG**(-1) + 4*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(69)*N*CK*TWO**(-1) * (  - 8*S**(-2)*TG*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2)*TG - 8*S**(-2)*M2*(TG+UG)**(-1)*
     +    S4 + 8*S**(-2)*(TG+UG)**(-1)*S4**2 - 8*S**(-2)*S4 + 5*S**(-1)
     +    *M2*(TG+UG)**(-1) + 8*S**(-1)*MS2*(TG+UG)**(-1) - 13*S**(-1)*
     +    (TG+UG)**(-1)*S4 + 10*S**(-1) + 5*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(70)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * ( 4*
     +    S**(-2)*TG - 2*S**(-1)*TG**(-1)*M2 - 2*S**(-1)*TG**(-1)*MS2
     +     + 2*S**(-1)*UG**(-1)*M2 + 2*S**(-1)*UG**(-1)*MS2 + 8*S**(-1)
     +     + 8*TG**(-1)*UG**(-1)*M2 + 8*TG**(-1)*UG**(-1)*MS2 + 
     +    TG**(-1) + 16*UG**(-2)*M2 + 16*UG**(-2)*MS2 + 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(70)*N*CO*TWO**(-1) * (  - 8*S**(-1)*TG**(-1)
     +     - 7*S**(-1)*UG**(-1) + 12*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(70)*N*CK*TWO**(-1) * ( 4*S**(-2)*TG*
     +    (TG+UG)**(-1) + 8*S**(-2)*M2*(TG+UG)**(-1) + 8*S**(-2)*
     +    (TG+UG)**(-1)*S4 - 8*S**(-2) - 8*S**(-1)*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + ANG2(71)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 2*
     +    S**(-1) + 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(71)*N*CO*TWO**(-1) * (  - 2*S**(-1)*TG**(-1)
     +     + 2*TG**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(72)*N*CO*S4G2**(-1)*TWO**(-1)*S4G * (  - 4*
     +    S**(-2)*TG + 2*S**(-1)*TG**(-1)*M2 + 2*S**(-1)*TG**(-1)*MS2
     +     - 2*S**(-1)*UG**(-1)*M2 - 2*S**(-1)*UG**(-1)*MS2 - TG**(-1)
     +     + 2*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(72)*N*CO*TWO**(-1) * (  - S**(-1)*UG**(-1) )
     +
      M2GGH = M2GGH + ANG2(72)*N*CK*TWO**(-1) * (  - 4*S**(-2)*TG*
     +    (TG+UG)**(-1) - 8*S**(-2)*M2*(TG+UG)**(-1) - 8*S**(-2)*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2) + 8*S**(-1)*(TG+UG)**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * ( 8*S**(-1)*TG
     +    *T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 4*S**(-1)*TG*
     +    T1**(-2) - 8*S**(-1)*TG*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4
     +     - 4*S**(-1)*TG*T1**(-1)*S4**(-1) + 16*S**(-1)*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2) - 16*S**(-1)*TG*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 8*S**(-1)*TG*(S+TG)**(-1)*
     +    S4**(-1) - 8*S**(-1)*TG*(S+UG)**(-1)*S4**(-1) - 16*S**(-1)*
     +    TG**2*M2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 8*S**(-1)*T1**(-2)
     +    *M2*MS2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*T1**(-2)*M2 - 8*
     +    S**(-1)*T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*
     +    T1**(-1)*M2*S4**(-1) - 8*S**(-1)*U1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 - 4*S**(-1)*U1**(-2)*M2 - 8*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 8*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 + 8*S**(-1)*U1**(-1)*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**2 - 4*S**(-1)*U1**(-1)*M2*S4**(-1)
     +     + 8*S**(-1)*U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 )
     +
      M2GGH = M2GGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)
     +    *U1**(-1) - 8*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 - 16*
     +    S**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 16*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 16*S**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*S**(-1)*
     +    (S+TG)**(-1) - 4*S**(-1)*S4**(-1) - 16*S*TG**(-1)*M2*
     +    (S+UG)**(-2)*S4**(-1) - 16*S*TG**(-1)*M2**2*(S+UG)**(-3)*
     +    S4**(-1) - 8*S*TG**(-1)*(S+UG)**(-1)*S4**(-1) + 32*S*TG*
     +    UG**(-1)*(S+TG)**(-3) + 16*S*TG*UG**(-1)*(S+TG)**(-2)*
     +    S4**(-1) - 16*S*TG*(S+UG)**(-3)*S4**(-1) - 16*S*TG**2*
     +    UG**(-1)*(S+TG)**(-3)*S4**(-1) - 32*S*UG**(-2)*M2*MS2*
     +    (S+TG)**(-2)*S4**(-1) - 32*S*UG**(-2)*M2**2*(S+TG)**(-2)*
     +    S4**(-1) + 32*S*UG**(-1)*M2*MS2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 32*S*UG**(-1)*M2*(S+TG)**(-2)*
     +    S4**(-1) - 16*S*UG**(-1)*(S+TG)**(-3)*S4 - 16*S*UG**(-1)*
     +    (S+TG)**(-2) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 8*S*
     +    UG**(-1)*(S+TG)**(-1)*S4**(-1) - 8*S*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 16*S*M2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 32*S*M2*(S+UG)**(-3)*S4**(-1)
     +     - 16*S*(S+TG)**(-2)*S4**(-1) - 16*S*(S+UG)**(-2)*S4**(-1) - 
     +    32*S**2*TG*UG**(-1)*(S+TG)**(-3)*S4**(-1) + 32*S**2*UG**(-1)*
     +    (S+TG)**(-3) + 16*S**2*UG**(-1)*(S+TG)**(-2)*S4**(-1) - 16*
     +    S**3*UG**(-1)*(S+TG)**(-3)*S4**(-1) + 16*TG**(-2)*T1**(-1)*M2
     +    *MS2*S4**(-1) + 16*TG**(-2)*T1**(-1)*M2**2*S4**(-1) + 32*
     +    TG**(-2)*M2*MS2*(S+UG)**(-1)*S4**(-1) + 32*TG**(-2)*M2**2*MS2
     +    *(S+UG)**(-2)*S4**(-1) + 32*TG**(-2)*M2**2*(S+UG)**(-1)*
     +    S4**(-1) + 32*TG**(-2)*M2**3*(S+UG)**(-2)*S4**(-1) - 16*
     +    TG**(-1)*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) + 16*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) + 32*TG**(-1)*M2*MS2*(S+UG)**(-2)*
     +    S4**(-1) - 16*TG**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) + 
     +    32*TG**(-1)*M2*(S+UG)**(-1)*S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 32*
     +    TG**(-1)*M2**2*MS2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1)
     +     + 64*TG**(-1)*M2**2*(S+UG)**(-2)*S4**(-1) - 32*TG*UG**(-2)*
     +    M2*MS2*(S+TG)**(-2)*S4**(-1) - 32*TG*UG**(-2)*M2**2*
     +    (S+TG)**(-2)*S4**(-1) + 32*TG*UG**(-1)*M2*MS2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 32*TG*UG**(-1)*M2*(S+TG)**(-2)*
     +    S4**(-1) - 16*TG*M2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 16*TG*
     +    M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*M2*
     +    (S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*(S+TG)**(-2)
     +    *S4**(-1) + 16*TG*(S+UG)**(-2)*S4**(-1) + 16*UG**(-2)*
     +    U1**(-1)*M2*MS2*S4**(-1) + 16*UG**(-2)*U1**(-1)*M2**2*
     +    S4**(-1) + 32*UG**(-2)*M2*MS2*(S+TG)**(-2) + 32*UG**(-2)*M2*
     +    MS2*(S+TG)**(-1)*S4**(-1) + 32*UG**(-2)*M2**2*(S+TG)**(-2) + 
     +    32*UG**(-2)*M2**2*(S+TG)**(-1)*S4**(-1) - 16*UG**(-1)*
     +    U1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) + 16*UG**(-1)*U1**(-1)*M2*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CO*(S4+MS2)*TWO**(-1) * (  - 32*
     +    UG**(-1)*M2*MS2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1) - 16*UG**(-1)
     +    *M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) + 32*UG**(-1)*M2*
     +    (S+TG)**(-2) + 32*UG**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 8*
     +    T1**(-1)*M2*(TG*UG-M2*S)**(-1) + 8*T1**(-1)*S4**(-1) - 8*
     +    U1**(-1)*M2*(TG*UG-M2*S)**(-1) + 8*U1**(-1)*S4**(-1) - 32*M2*
     +    MS2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 8*M2*MS2*
     +    (TG*UG-M2*S)**(-2) - 16*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 
     +    48*M2*(S+UG)**(-2)*S4**(-1) - 16*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 8*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 16*M2**2
     +    *(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 16*(S+TG)**(-2)
     +     + 24*(S+TG)**(-1)*S4**(-1) + 16*(S+UG)**(-1)*S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * (  - 24*
     +    S**(-1)*TG*T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 + 12*
     +    S**(-1)*TG*T1**(-2) + 24*S**(-1)*TG*T1**(-1)*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4 + 12*S**(-1)*TG*T1**(-1)*S4**(-1) - 48*
     +    S**(-1)*TG*M2*MS2*(TG*UG-M2*S)**(-2) + 48*S**(-1)*TG*M2**2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 24*S**(-1)*TG*(S+TG)**(-1)*
     +    S4**(-1) + 24*S**(-1)*TG*(S+UG)**(-1)*S4**(-1) + 48*S**(-1)*
     +    TG**2*M2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 24*S**(-1)*
     +    T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-1)*S4 + 24*S**(-1)*T1**(-2)*
     +    M2 + 24*S**(-1)*T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 + 24*
     +    S**(-1)*T1**(-1)*M2*S4**(-1) + 24*S**(-1)*U1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 + 12*S**(-1)*U1**(-2)*M2 + 24*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 + 24*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 24*S**(-1)*U1**(-1)*M2
     +    *MS2*(TG*UG-M2*S)**(-2)*S4**2 + 12*S**(-1)*U1**(-1)*M2*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * (  - 24*
     +    S**(-1)*U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 + 12*S**(-1)
     +    *U1**(-1) + 24*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 + 48*
     +    S**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) + 24*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) + 48*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     + 48*S**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) + 24*S**(-1)*
     +    (S+TG)**(-1) + 12*S**(-1)*S4**(-1) + 16*S*TG**(-1)*M2*
     +    (S+UG)**(-2)*S4**(-1) + 16*S*TG**(-1)*M2**2*(S+UG)**(-3)*
     +    S4**(-1) + 8*S*TG**(-1)*(S+UG)**(-1)*S4**(-1) - 32*S*TG*
     +    UG**(-1)*(S+TG)**(-3) - 16*S*TG*UG**(-1)*(S+TG)**(-2)*
     +    S4**(-1) + 16*S*TG*(S+UG)**(-3)*S4**(-1) + 16*S*TG**2*
     +    UG**(-1)*(S+TG)**(-3)*S4**(-1) + 32*S*UG**(-2)*M2*MS2*
     +    (S+TG)**(-2)*S4**(-1) + 32*S*UG**(-2)*M2**2*(S+TG)**(-2)*
     +    S4**(-1) - 32*S*UG**(-1)*M2*MS2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 32*S*UG**(-1)*M2*(S+TG)**(-2)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * ( 16*S*
     +    UG**(-1)*(S+TG)**(-3)*S4 + 16*S*UG**(-1)*(S+TG)**(-2) + 8*S*
     +    UG**(-1)*(S+TG)**(-1)*S4**(-1) + 24*S*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) - 48*S*M2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) + 32*S*M2*(S+UG)**(-3)*S4**(-1)
     +     + 16*S*(S+TG)**(-2)*S4**(-1) + 16*S*(S+UG)**(-2)*S4**(-1) + 
     +    32*S**2*TG*UG**(-1)*(S+TG)**(-3)*S4**(-1) - 32*S**2*UG**(-1)*
     +    (S+TG)**(-3) - 16*S**2*UG**(-1)*(S+TG)**(-2)*S4**(-1) + 16*
     +    S**3*UG**(-1)*(S+TG)**(-3)*S4**(-1) - 16*TG**(-2)*T1**(-1)*M2
     +    *MS2*S4**(-1) - 16*TG**(-2)*T1**(-1)*M2**2*S4**(-1) - 32*
     +    TG**(-2)*M2*MS2*(S+UG)**(-1)*S4**(-1) - 32*TG**(-2)*M2**2*MS2
     +    *(S+UG)**(-2)*S4**(-1) - 32*TG**(-2)*M2**2*(S+UG)**(-1)*
     +    S4**(-1) - 32*TG**(-2)*M2**3*(S+UG)**(-2)*S4**(-1) + 16*
     +    TG**(-1)*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) - 16*TG**(-1)*
     +    T1**(-1)*M2*S4**(-1) - 32*TG**(-1)*M2*MS2*(S+UG)**(-2)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * ( 16*TG**(-1)*
     +    M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 32*TG**(-1)*M2*
     +    (S+UG)**(-1)*S4**(-1) + 32*TG**(-1)*M2**2*MS2*(S+UG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 64*TG**(-1)*M2**2*(S+UG)**(-2)*
     +    S4**(-1) + 32*TG*UG**(-2)*M2*MS2*(S+TG)**(-2)*S4**(-1) + 32*
     +    TG*UG**(-2)*M2**2*(S+TG)**(-2)*S4**(-1) - 32*TG*UG**(-1)*M2*
     +    MS2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 32*TG*UG**(-1)
     +    *M2*(S+TG)**(-2)*S4**(-1) + 48*TG*M2*MS2*(TG*UG-M2*S)**(-2)*
     +    S4**(-1) - 48*TG*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1)
     +     + 48*TG*M2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 16*TG*
     +    (S+TG)**(-2)*S4**(-1) - 16*TG*(S+UG)**(-2)*S4**(-1) - 16*
     +    UG**(-2)*U1**(-1)*M2*MS2*S4**(-1) - 16*UG**(-2)*U1**(-1)*
     +    M2**2*S4**(-1) - 32*UG**(-2)*M2*MS2*(S+TG)**(-2) - 32*
     +    UG**(-2)*M2*MS2*(S+TG)**(-1)*S4**(-1) - 32*UG**(-2)*M2**2*
     +    (S+TG)**(-2) - 32*UG**(-2)*M2**2*(S+TG)**(-1)*S4**(-1) + 16*
     +    UG**(-1)*U1**(-1)*M2*MS2*(TG*UG-M2*S)**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*N*CK*(S4+MS2)*TWO**(-1) * (  - 16*
     +    UG**(-1)*U1**(-1)*M2*S4**(-1) + 32*UG**(-1)*M2*MS2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1) + 16*UG**(-1)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 32*UG**(-1)*M2*(S+TG)**(-2) - 
     +    32*UG**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 24*T1**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) - 8*T1**(-1)*S4**(-1) + 24*U1**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) - 8*U1**(-1)*S4**(-1) + 32*M2*MS2*
     +    (S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 24*M2*MS2*
     +    (TG*UG-M2*S)**(-2) + 48*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1) - 
     +    48*M2*(S+UG)**(-2)*S4**(-1) + 48*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) + 24*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) + 48*
     +    M2**2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 16*
     +    (S+TG)**(-2) - 40*(S+TG)**(-1)*S4**(-1) - 16*(S+UG)**(-1)*
     +    S4**(-1) )
     +
      M2GGH = M2GGH + COLO2(9)*CQED*(S4+MS2)*TWO**(-1) * ( 8*S**(-1)*TG
     +    *T1**(-2)*M2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 4*S**(-1)*TG*
     +    T1**(-2) - 8*S**(-1)*TG*T1**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4
     +     - 4*S**(-1)*TG*T1**(-1)*S4**(-1) + 16*S**(-1)*TG*M2*MS2*
     +    (TG*UG-M2*S)**(-2) - 16*S**(-1)*TG*M2**2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 8*S**(-1)*TG*(S+TG)**(-1)*
     +    S4**(-1) - 8*S**(-1)*TG*(S+UG)**(-1)*S4**(-1) - 16*S**(-1)*
     +    TG**2*M2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 8*S**(-1)*T1**(-2)
     +    *M2*MS2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*T1**(-2)*M2 - 8*
     +    S**(-1)*T1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 - 8*S**(-1)*
     +    T1**(-1)*M2*S4**(-1) - 8*S**(-1)*U1**(-2)*M2*MS2*
     +    (TG*UG-M2*S)**(-1)*S4 - 4*S**(-1)*U1**(-2)*M2 - 8*S**(-1)*
     +    U1**(-2)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**2 - 8*S**(-1)*
     +    U1**(-2)*M2**2*(TG*UG-M2*S)**(-1)*S4 + 8*S**(-1)*U1**(-1)*M2*
     +    MS2*(TG*UG-M2*S)**(-2)*S4**2 - 4*S**(-1)*U1**(-1)*M2*S4**(-1)
     +     + 8*S**(-1)*U1**(-1)*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4 )
     +
      M2GGH = M2GGH + COLO2(9)*CQED*(S4+MS2)*TWO**(-1) * (  - 4*S**(-1)
     +    *U1**(-1) - 8*S**(-1)*M2*MS2*(TG*UG-M2*S)**(-2)*S4 - 16*
     +    S**(-1)*M2*MS2*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 16*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 16*S**(-1)*M2**2*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*S**(-1)*
     +    (S+TG)**(-1) - 4*S**(-1)*S4**(-1) - 8*S*M2*MS2*
     +    (TG*UG-M2*S)**(-2)*S4**(-1) + 16*S*M2*(S+TG)**(-1)*
     +    (TG*UG-M2*S)**(-1)*S4**(-1) - 16*TG*M2*MS2*(TG*UG-M2*S)**(-2)
     +    *S4**(-1) + 16*TG*M2*(S+TG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1)
     +     - 16*TG*M2*(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) - 8*
     +    T1**(-1)*M2*(TG*UG-M2*S)**(-1) - 8*U1**(-1)*M2*
     +    (TG*UG-M2*S)**(-1) + 8*M2*MS2*(TG*UG-M2*S)**(-2) - 16*M2*
     +    (S+TG)**(-1)*(TG*UG-M2*S)**(-1) - 16*M2*(TG*UG-M2*S)**(-1)*
     +    S4**(-1) - 8*M2**2*MS2*(TG*UG-M2*S)**(-2)*S4**(-1) - 16*M2**2
     +    *(S+UG)**(-1)*(TG*UG-M2*S)**(-1)*S4**(-1) + 8*(S+TG)**(-1)*
     +    S4**(-1) )


      Q2MS2 = (MS + MG)**2/4.D0/MS**2

      M2GGH = 2.D0*(NS -1.D0) * M2GGH

      DSGGGH = S4/(S4+MS2)/2.D0 *
     +     ALPHAS**3 * AVG * M2GGH /4.D0 /S**2 *CONV
     +     + DSGGG3(ALPHAS,S,TG,S4,MS,MG,Q2MS2)

      RETURN
      END


      REAL*8 FUNCTION DSGGGS(ALPHAS,S,TG,S4,MS,MG)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR G +G -> SQ + GL +QB
C***  THE 1/S4G**2 PART
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,TG,U1,UG,MS,MG,MS2,MG2,M2,CO,CK,CQED
      REAL*8 M2GGS, NS,CONV,N,AVG

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2

      U1 = S4 -S -TG
      UG = U1 -M2
      

      M2GGS = 0.D0
      M2GGS = M2GGS + N*CO * (  - 16*S**(-2)*TG**2*M2 - 16*S**(-1)*TG*
     +    M2 - 32*S**(-1)*M2*MS2 - 32*S**(-1)*M2**2 - 8*S*TG**(-1)*M2
     +     - 8*S*UG**(-1)*M2 - 32*TG**(-2)*M2*MS2**2 - 64*TG**(-2)*
     +    M2**2*MS2 - 32*TG**(-2)*M2**3 - 32*TG**(-1)*UG**(-1)*M2*
     +    MS2**2 - 64*TG**(-1)*UG**(-1)*M2**2*MS2 - 32*TG**(-1)*
     +    UG**(-1)*M2**3 - 32*TG**(-1)*M2*MS2 - 32*TG**(-1)*M2**2 - 32*
     +    UG**(-2)*M2*MS2**2 - 64*UG**(-2)*M2**2*MS2 - 32*UG**(-2)*
     +    M2**3 - 32*UG**(-1)*M2*MS2 - 32*UG**(-1)*M2**2 - 24*M2 )


      M2GGS = 2.D0* (NS -1.D0)* M2GGS

      DSGGGS = ALPHAS**3 *AVG *M2GGS *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END





      REAL*8 FUNCTION DSGGGT(ALPHAS,S,TG,S4,S3,MS,MG)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR G +G -> SQ + GL +QB
C***  THE 1/S3**2 PART

      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,S4G,S3,T1,TG,U1,UG,MS,MG,MS2,MG2,M2
      REAL*8 M2GGT, NS,CONV,N,AVG,CO,CK,CQED
      REAL*8 ANGDEF(1:11), ANA(2:2,1:9), ANB(2:2,1:9), ANC(2:2,1:9)
      REAL*8 ANGS3(1:6),XX(5:9),YY2(5:9),XPHI

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CO = (N**2 -1.D0)*N
      CK = (N**2 -1.D0)/N
      CQED = (N**4 -1.D0)/N**2
      AVG = (1.D0/2.D0)**2 /(N**2 -1.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2

      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)

      XPHI = (S3 -ANA(2,1))/ANB(2,1)

      XX(5) = (ANA(2,5) + ANB(2,5)*XPHI)
      YY2(5)= ANC(2,5)**2 * (1.D0-XPHI**2)
      XX(6) = (ANA(2,6) + ANB(2,6)*XPHI)
      YY2(6)= ANC(2,6)**2 * (1.D0-XPHI**2)
      XX(8) = (ANA(2,8) + ANB(2,8)*XPHI)
      YY2(8)= ANC(2,8)**2 * (1.D0-XPHI**2)
      XX(9) = (ANA(2,9) + ANB(2,9)*XPHI)
      YY2(9)= ANC(2,9)**2 * (1.D0-XPHI**2)

      ANGS3(1) = -XX(6)/(XX(6)**2 - YY2(6))**(1.5D0)
      ANGS3(2) = -1.D0/SQRT(XX(6)**2 - YY2(6))
      ANGS3(3) = -XX(5)/(XX(5)**2 - YY2(5))**(1.5D0)
      ANGS3(4) = -1.D0/SQRT(XX(5)**2 - YY2(5))
      ANGS3(5) = XX(5)
      ANGS3(6) = XX(5)**2 +0.5D0*YY2(5)

      M2GGT = 0.D0
      M2GGT = M2GGT + N*CO * (  - 8*S**(-1)*M2*MS2 - 2*M2 )
     +
      M2GGT = M2GGT + N*CK * ( 8*S**(-1)*M2*MS2 + 6*M2 )
     +
      M2GGT = M2GGT + CQED * (  - 2*M2 )
     +
      M2GGT = M2GGT + ANGS3(1)*N*CO * (  - 4*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(1)*N*CK * ( 12*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(1)*CQED * (  - 4*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(2)*N*CO * (  - 4*M2*MS2 )
     +
      M2GGT = M2GGT + ANGS3(2)*N*CK * (  - 16*S**(-1)*M2*MS2**2 + 12*M2
     +    *MS2 )
     +
      M2GGT = M2GGT + ANGS3(2)*CQED * ( 8*S**(-1)*M2*MS2**2 - 4*M2*MS2
     +     )
     +
      M2GGT = M2GGT + ANGS3(3)*N*CO * (  - 4*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(3)*N*CK * ( 12*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(3)*CQED * (  - 4*M2*MS2**2 )
     +
      M2GGT = M2GGT + ANGS3(4)*N*CO * (  - 4*M2*MS2 )
     +
      M2GGT = M2GGT + ANGS3(4)*N*CK * (  - 16*S**(-1)*M2*MS2**2 + 12*M2
     +    *MS2 )
     +
      M2GGT = M2GGT + ANGS3(4)*CQED * ( 8*S**(-1)*M2*MS2**2 - 4*M2*MS2
     +     )
     +
      M2GGT = M2GGT + ANGS3(5)*N*CO * (  - 4*S**(-1)*M2 )
     +
      M2GGT = M2GGT + ANGS3(5)*N*CK * ( 4*S**(-1)*M2 )
     +
      M2GGT = M2GGT + ANGS3(6)*N*CO * (  - 4*S**(-2)*M2 )
     +
      M2GGT = M2GGT + ANGS3(6)*N*CK * ( 4*S**(-2)*M2 )


      M2GGT = 2.D0* (NS -1.D0)* M2GGT

      DSGGGT = ALPHAS**3 *AVG *M2GGT *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSGQB3(ALPHAS,S,TG,S4,MS,MG,SCA,IFL)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,TG,U1,UG,S4,MS,MG,SCA,X2
      REAL*8 PI,N,CF,NS, DSGQGB, PGQ2
      INTEGER IFL

      N = 3.D0
      CF = (N**2 -1.D0)/2.D0/N
      NS = 6.D0
      PI = 4.D0*ATAN(1.D0)
      U1 = S4 -S -TG
      T1 = TG +MG**2 -MS**2
      UG = U1 -MG**2 +MS**2
      X2 = -U1/(S+TG)
      PGQ2 = CF*(1.D0 + (1.D0 -X2)**2)/X2

      DSGQB3 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) * 
     +     (-1.D0/U1*DSGQGB(ALPHAS,X2*S,X2*TG,MS,MG)*PGQ2*X2**2 )

      RETURN
      END


      REAL*8 FUNCTION DSGQBH(ALPHAS,S,TG,S4,MS,MG,IFL,EPSS,EPSG)
C***  CROSS SECTIONS FOR Q + QB -> SQ + GL +QB
C***  M2QBH1: Q QBP --> GL SQ QBP
C***  M2QBH2: Q QB  --> GL SQP QBP
C***  M2QBH3: INTERFERENCE 

      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,UG,MS,MG,MS2,MG2,M2,EPSS,EPSG
      REAL*8 NS,CONV,N,CF,AVG,S4G, S4G2,DEL,SYMBT,SYMBU,SYMBY
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 ABP1P1, ABP1M1, ABP1P2, A4P1P2, A4P1P1
      REAL*8 A4M1P2, A4M1P1, A4P0P2, A4P1P0, A4P2M2, A4M2P1,ABP2P2
      REAL*8 ABP2P1, ABP2P0, A4M1P0
      REAL*8 C4P1P2, C4P1P1, C4M1P1, C4P1P0, CBP1P1, C4M2P1
      REAL*8 M2QBH,M2QBH1,M2QBH2,M2QBH3
      REAL*8 ANG2(1:74), COLO2(1:9),DSGQB3,Q2MS2
      INTEGER IFL
           
      CONV = 389379660.D0
      NS = 6.D0

      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0

      AVG = (1.D0/2.D0)**2 /N**2
      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2
      S4G2 = S4G**2 + EPSG*MG**4
      SYMBT = (M2*S-TG*S4G)/((M2*S-TG*S4G)**2 +(S+TG)*EPSS*MS**4 )
      SYMBU = (M2*S-UG*S4G)/((M2*S-UG*S4G)**2 +(S+UG)*EPSS*MS**4 )
      SYMBY = 1.D0/((TG-M2)**2 +EPSS*MS**4)

      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6) + M2
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9) 
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) + MS2 +MG2
      ANB(1,2) = -ANB(1,1)
      ANC(1,2) = -ANC(1,1)
      ANA(1,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(1,3) = -ANA(1,3)
      ANC(1,3) =  0.D0
      ANA(1,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(1,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)-2.D0*ANGDEF(3)*ANGDEF(4)
      ANC(1,4) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(1,5) = -ANB(1,3)
      ANC(1,5) = -ANC(1,3)
      ANA(1,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(1,6) = -ANB(1,4)
      ANC(1,6) = -ANC(1,4)
      ANA(1,7) = +ANA(1,1) -M2
      ANB(1,7) = +ANB(1,1)
      ANC(1,7) = +ANC(1,1)
      ANA(1,8) = +ANA(1,5) -M2
      ANB(1,8) = +ANB(1,5)
      ANC(1,8) = +ANC(1,5)
      ANA(1,9) = +ANA(1,6) -M2
      ANB(1,9) = +ANB(1,6)
      ANC(1,9) = +ANC(1,6)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(3,2) = -ANB(3,1)
      ANC(3,2) = -ANC(3,1)
      ANA(3,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(3,3) = 
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10) -2.D0*ANGDEF(2)*ANGDEF(4)
      ANC(3,3) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(3,4) = -ANA(3,4)
      ANC(3,4) =  0.D0
      ANA(3,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(3,5) = -ANB(3,3)
      ANC(3,5) = -ANC(3,3)
      ANA(3,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(3,6) = -ANB(3,4)
      ANC(3,6) = -ANC(3,4)
      ANA(3,7) = +ANA(3,1) -M2
      ANB(3,7) = +ANB(3,1)
      ANC(3,7) = +ANC(3,1)
      ANA(3,8) = +ANA(3,5) -M2
      ANB(3,8) = +ANB(3,5)
      ANC(3,8) = +ANC(3,5)
      ANA(3,9) = +ANA(3,6) -M2
      ANB(3,9) = +ANB(3,6)
      ANC(3,9) = +ANC(3,6)

C$$$      ANG2(1) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(2) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(3) = A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(4) = A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(5) = A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(6) = A4M1P1(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(7) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(8) = A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(9) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(10)= A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))

C$$$      ANG2(11)= A4P0P2(ANA(1,7),ANB(1,7),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(12)= A4P1P0(ANA(2,7),ANB(2,7),ANA(2,7),ANB(2,7),ANC(2,7))
      ANG2(13)= A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(14)= A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(15)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(16)= A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(17)= ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(18)= A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(19)= A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(20)= AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))

C$$$      ANG2(21)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(22)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(23)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(24)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(25)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(26)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(27)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(28)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(29)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(30)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))

      ANG2(31)= ABP1M1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(32)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(33)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(34)= ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(35)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(36)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(37)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(38)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG2(39)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(40)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))

      ANG2(41)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(42)= A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(43)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(44)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(45)= A4P2M2(ANA(2,2),ANB(2,2),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(46)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(47)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(48)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG2(49)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
      ANG2(50)= A4M1P0(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))

      ANG2(51)= A4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(52)= A4P1P1(ANA(2,2),ANB(2,2),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(53)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(54)= A4P1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(55)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(56)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(57)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(59)= ABP1P1(ANA(1,4),ANB(1,4),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(60)= ABP1M2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))

C$$$      ANG2(61)= ABP1M2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(62)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(63)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(64)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(65)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(66)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG2(67)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(68)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(69)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(70)= A4M1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG2(71)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG2(72)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(73)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(74)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))

      DEL = EPSS * MS**4
      IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
C$$$  ANG2(47)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$  ANG2(48)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
         ANG2(49)= C4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
         ANG2(51)= C4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$  ANG2(56)= C4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$  ANG2(57)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= CBP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1),DEL)
C$$$  ANG2(59)= CBP1P1(ANA(1,4),ANB(1,4),ANA(1,1),ANB(1,1),ANC(1,1),DEL)
         ANG2(66)= C4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
         ANG2(69)= C4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
         ANG2(72)= C4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
         ANG2(73)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
         ANG2(74)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
      END IF


      COLO2(9) = LOG(S4**2/MS2/(S4+MS2))

      M2QBH1 = 0D0
      M2QBH1 = M2QBH1 + N*CF**2*(S4+MS2) * (  - 32*S**(-1)*U1**(-4)*
     +    M2**2*S4 + 32*S**(-1)*U1**(-3)*M2*S4 + 32*S**(-1)*U1**(-3)*
     +    M2**2 + 16*S**(-1)*U1**(-2)*M2*MS2*(S+TG)**(-1) - 16*S**(-1)*
     +    U1**(-2)*M2*MS2*S4**(-1) - 32*S**(-1)*U1**(-2)*M2 - 16*
     +    S**(-1)*U1**(-2)*M2**2*S4**(-1) - 16*S**(-1)*U1**(-1)*M2*MS2*
     +    (S+TG)**(-1)*S4**(-1) + 16*S**(-1)*U1**(-1)*M2*S4**(-1) - 8*
     +    S**(-1)*S4**(-1) + 16*S*TG**(-1)*M2*(TG+UG)**(-1)*SYMBY - 16*
     +    S*TG**(-1)*M2*SYMBY*S4**(-1) - 16*S*TG**(-1)*(TG+UG)**(-1)*
     +    S4**(-1) - 16*S*TG*(TG+UG)**(-1)*SYMBT*SYMBY*S4 - 16*S*TG*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) + 16*S*TG*SYMBT*SYMBY - 16*S*
     +    TG**2*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 16*S*TG**3*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*S*T1**(-1)*
     +    (S+TG)**(-1)*S4**(-1) + 16*S*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4
     +     - 16*S*M2*SYMBT*SYMBY - 16*S*(TG+UG)**(-1)*SYMBY - 16*S*
     +    SYMBT*S4**(-1) + 16*S*SYMBY*S4**(-1) - 16*S**2*TG**(-1)*M2*
     +    (TG+UG)**(-1)*SYMBY*S4**(-1) )
     +
      M2QBH1 = M2QBH1 + N*CF**2*(S4+MS2) * ( 16*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 32*S**2*TG*(TG+UG)**(-1)
     +    *SYMBT*SYMBY - 16*S**2*TG*SYMBT*SYMBY*S4**(-1) - 16*S**2*
     +    TG**2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 32*S**2*M2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY + 16*S**2*M2*SYMBT*SYMBY*S4**(-1)
     +     + 16*S**2*(TG+UG)**(-1)*SYMBY*S4**(-1) - 16*S**3*TG*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 16*S**3*M2*(TG+UG)**(-1)
     +    *SYMBT*SYMBY*S4**(-1) - 16*TG*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4
     +     + 16*TG*M2*(TG+UG)**(-1)*SYMBY*S4**(-1) + 16*TG*M2*SYMBT*
     +    SYMBY + 16*TG*(TG+UG)**(-1)*SYMBT + 16*TG*(TG+UG)**(-1)*SYMBY
     +     - 16*TG*SYMBT*S4**(-1) - 16*TG*SYMBY*S4**(-1) + 32*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY - 16*TG**2*M2*SYMBT*SYMBY*S4**(-1)
     +     + 16*TG**2*(TG+UG)**(-1)*SYMBT*SYMBY*S4 - 16*TG**2*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 16*TG**2*(TG+UG)**(-1)*SYMBY*
     +    S4**(-1) - 16*TG**2*SYMBT*SYMBY - 16*TG**3*M2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4**(-1) )
     +
      M2QBH1 = M2QBH1 + N*CF**2*(S4+MS2) * (  - 32*TG**3*(TG+UG)**(-1)*
     +    SYMBT*SYMBY + 16*TG**3*SYMBT*SYMBY*S4**(-1) + 16*TG**4*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 16*T1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) + 16*T1**(-1)*S4**(-1) - 32*U1**(-4)*M2
     +    *MS2 - 32*U1**(-3)*M2 - 16*U1**(-2)*M2*MS2*(S+TG)**(-1)*
     +    S4**(-1) - 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) - 16*M2*
     +    (TG+UG)**(-1)*SYMBY + 16*M2*SYMBY*S4**(-1) - 8*(S+TG)**(-1)*
     +    S4**(-1) )
     +
      M2QBH1 = M2QBH1 + N**2*CF*(S4+MS2) * (  - 4*S*TG**(-1)*M2*
     +    (TG+UG)**(-1)*SYMBY + 4*S*TG**(-1)*M2*SYMBY*S4**(-1) + 4*S*
     +    TG**(-1)*(TG+UG)**(-1)*S4**(-1) + 4*S*TG*M2*SYMBT*SYMBY*
     +    S4**(-1) + 4*S*TG*(TG+UG)**(-1)*SYMBT*SYMBY*S4 + 4*S*TG*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 4*S*TG*SYMBT*SYMBY + 4*S*TG**2
     +    *M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 4*S*TG**2*SYMBT*
     +    SYMBY*S4**(-1) - 4*S*TG**3*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1)
     +     + 4*S*T1**(-1)*(S+TG)**(-1)*S4**(-1) - 4*S*M2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4 + 4*S*M2*SYMBT*SYMBY + 4*S*(TG+UG)**(-1)*SYMBY
     +     + 8*S*SYMBT*S4**(-1) - 4*S*SYMBY*S4**(-1) + 4*S**2*TG**(-1)*
     +    M2*(TG+UG)**(-1)*SYMBY*S4**(-1) - 4*S**2*TG*M2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4**(-1) - 8*S**2*TG*(TG+UG)**(-1)*SYMBT*SYMBY + 
     +    4*S**2*TG*SYMBT*SYMBY*S4**(-1) + 4*S**2*TG**2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4**(-1) + 8*S**2*M2*(TG+UG)**(-1)*SYMBT*SYMBY - 
     +    4*S**2*M2*SYMBT*SYMBY*S4**(-1) - 4*S**2*(TG+UG)**(-1)*SYMBY*
     +    S4**(-1) )
     +
      M2QBH1 = M2QBH1 + N**2*CF*(S4+MS2) * ( 4*S**3*TG*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4**(-1) - 4*S**3*M2*(TG+UG)**(-1)*SYMBT*SYMBY*
     +    S4**(-1) - 32*TG**(-2)*U1**(-4)*M2*MS2*S4**2 - 32*TG**(-2)*
     +    U1**(-4)*M2**2*S4**2 + 64*TG**(-2)*U1**(-3)*M2*MS2*S4 + 64*
     +    TG**(-2)*U1**(-3)*M2**2*S4 - 48*TG**(-2)*U1**(-2)*M2*MS2 - 48
     +    *TG**(-2)*U1**(-2)*M2**2 + 16*TG**(-2)*U1**(-1)*M2*MS2*
     +    S4**(-1) + 16*TG**(-2)*U1**(-1)*M2**2*S4**(-1) + 32*TG**(-1)*
     +    U1**(-4)*M2*MS2*S4 - 32*TG**(-1)*U1**(-3)*M2*MS2 + 32*
     +    TG**(-1)*U1**(-3)*M2*S4 + 16*TG**(-1)*U1**(-2)*M2*MS2*
     +    S4**(-1) - 32*TG**(-1)*U1**(-2)*M2 + 16*TG**(-1)*U1**(-1)*M2*
     +    S4**(-1) - 8*TG**(-1)*S4**(-1) + 4*TG*M2*(TG+UG)**(-1)*SYMBT*
     +    SYMBY*S4 - 4*TG*M2*(TG+UG)**(-1)*SYMBY*S4**(-1) - 8*TG*M2*
     +    SYMBT*SYMBY - 4*TG*(TG+UG)**(-1)*SYMBT - 4*TG*(TG+UG)**(-1)*
     +    SYMBY + 8*TG*SYMBT*S4**(-1) + 8*TG*SYMBY*S4**(-1) - 8*TG**2*
     +    M2*(TG+UG)**(-1)*SYMBT*SYMBY + 8*TG**2*M2*SYMBT*SYMBY*
     +    S4**(-1) )
     +
      M2QBH1 = M2QBH1 + N**2*CF*(S4+MS2) * (  - 4*TG**2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4 + 4*TG**2*(TG+UG)**(-1)*SYMBT*S4**(-1) + 4*
     +    TG**2*(TG+UG)**(-1)*SYMBY*S4**(-1) + 8*TG**2*SYMBT*SYMBY + 4*
     +    TG**3*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 8*TG**3*
     +    (TG+UG)**(-1)*SYMBT*SYMBY - 8*TG**3*SYMBT*SYMBY*S4**(-1) - 4*
     +    TG**4*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 4*T1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 4*T1**(-1)*S4**(-1) + 4*M2*
     +    (TG+UG)**(-1)*SYMBY - 8*M2*SYMBY*S4**(-1) + 4*(S+TG)**(-1)*
     +    S4**(-1) )
     +
      M2QBH1 = M2QBH1 + ANG2(13)*N*CF**2 * ( 4*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(14)*N*CF**2 * ( 4 + 4*S*U1**(-1) + 8*S*
     +    (TG+UG)**(-1) - 16*T1**(-1)*M2 - 16*T1**(-1)*MS2 - 4*U1**(-1)
     +    *M2 - 8*U1**(-1)*MS2 - 4*U1**(-1)*S4 - 8*M2*(TG+UG)**(-1) - 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(14)*N**2*CF * (  - 2*S*U1**(-1) - 2*S*
     +    (TG+UG)**(-1) + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 2*U1**(-1)*
     +    M2 + 4*U1**(-1)*MS2 + 2*U1**(-1)*S4 + 2*M2*(TG+UG)**(-1) + 2*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(15)*N*CF**2 * ( 4*T1**(-2)*M2**3 + 8*
     +    T1**(-1)*M2*MS2 + 8*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(16)*N*CF**2 * (  - 4 - 8*S*U1**(-1) - 4*
     +    T1**(-2)*M2*S4 + 4*T1**(-2)*M2**2 - 8*T1**(-1)*U1**(-1)*M2*S4
     +     - 16*T1**(-1)*U1**(-1)*MS2*S4 - 8*T1**(-1)*U1**(-1)*S4**2 + 
     +    4*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 4*T1**(-1)*S4 + 16*U1**(-1)*
     +    M2 + 16*U1**(-1)*MS2 + 8*U1**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(16)*N**2*CF * ( 2 + 2*S*U1**(-1) + 2*
     +    T1**(-1)*U1**(-1)*M2*S4 + 4*T1**(-1)*U1**(-1)*MS2*S4 + 2*
     +    T1**(-1)*U1**(-1)*S4**2 - 4*U1**(-1)*M2 - 4*U1**(-1)*MS2 - 2*
     +    U1**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(17)*N*CF**2 * ( 8*U1**(-2)*M2*S4**2 - 16*
     +    U1**(-1)*M2*S4 + 8*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(27)*N**2*CF * ( 8*TG**2*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(28)*N**2*CF * ( 12*TG*M2 + 4*TG**2*M2**2*
     +    SYMBY - 4*TG**3*M2*SYMBY + 4*T1**(-1)*M2**3 + 8*M2*MS2 + 4*
     +    M2**2 )
     +
      M2QBH1 = M2QBH1 + ANG2(29)*N**2*CF * (  - 8*S*M2 - 8*S**2*
     +    U1**(-1)*M2 - 8*TG*M2 + 8*U1**(-1)*M2*S4**2 - 8*M2*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(30)*N*CF**2 * ( 24*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 8*S*TG*M2*SYMBY + 8*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 16*S*TG*(TG+UG)**(-1) - 8*S*TG*SYMBY*S4 - 16*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 32*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 16*
     +    S*TG**2*SYMBY + 16*S*TG**3*(TG+UG)**(-1)*SYMBY - 16*S*
     +    U1**(-1)*M2*(TG+UG)**(-1)*S4 - 32*S*U1**(-1)*M2*SYMBY*S4**2
     +     - 16*S*U1**(-1)*M2 + 16*S*U1**(-1)*M2**2*SYMBY*S4 - 8*S*M2*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 8*S*M2*(TG+UG)**(-1) + 24*S*M2*
     +    SYMBY*S4 + 8*S*M2**2*(TG+UG)**(-1)*SYMBY*S4 - 8*S*M2**2*SYMBY
     +     + 8*S*(TG+UG)**(-1)*S4 - 8*S**2*TG*M2*(TG+UG)**(-1)*SYMBY - 
     +    8*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 + 16*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBY + 16*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 16*
     +    S**2*U1**(-1)*M2*SYMBY*S4 + 8*S**2*M2*(TG+UG)**(-1)*SYMBY*S4
     +     - 8*S**2*M2**2*(TG+UG)**(-1)*SYMBY - 8*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4**2 - 8*TG*M2*(TG+UG)**(-1) - 8*TG*M2*SYMBY*S4 + 8*TG
     +    *M2**2*(TG+UG)**(-1)*SYMBY*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(30)*N*CF**2 * ( 8*TG*M2**2*SYMBY + 8*TG*
     +    (TG+UG)**(-1)*S4 - 8*TG**2*M2*SYMBY - 8*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*
     +    TG**2*SYMBY*S4 + 8*TG**3*M2*(TG+UG)**(-1)*SYMBY - 8*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 + 16*T1**(-1)*U1**(-1)*M2**3 + 16*
     +    U1**(-1)*M2*SYMBY*S4**3 + 16*U1**(-1)*M2*S4 - 16*U1**(-1)*
     +    M2**2*SYMBY*S4**2 - 16*U1**(-1)*M2**2 - 16*M2*SYMBY*S4**2 - 
     +    16*M2 + 16*M2**2*SYMBY*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(30)*N**2*CF * (  - 6*S*TG*M2*(TG+UG)**(-1)
     +    *SYMBY*S4 + 4*S*TG*M2*SYMBY - 2*S*TG*(TG+UG)**(-1)*SYMBY*
     +    S4**2 + 4*S*TG*(TG+UG)**(-1) + 2*S*TG*SYMBY*S4 + 4*S*TG**2*M2
     +    *(TG+UG)**(-1)*SYMBY + 8*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 - 6*S
     +    *TG**2*SYMBY - 4*S*TG**3*(TG+UG)**(-1)*SYMBY + 8*S*U1**(-2)*
     +    M2*S4 + 4*S*U1**(-1)*M2*(TG+UG)**(-1)*S4 + 16*S*U1**(-1)*M2*
     +    SYMBY*S4**2 - 4*S*U1**(-1)*M2 - 8*S*U1**(-1)*M2**2*SYMBY*S4
     +     - 4*S*U1**(-1)*S4 + 2*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 2*S*
     +    M2*(TG+UG)**(-1) - 10*S*M2*SYMBY*S4 - 2*S*M2**2*(TG+UG)**(-1)
     +    *SYMBY*S4 + 2*S*M2**2*SYMBY - 2*S*(TG+UG)**(-1)*S4 + 4*S + 2*
     +    S**2*TG*M2*(TG+UG)**(-1)*SYMBY + 2*S**2*TG*(TG+UG)**(-1)*
     +    SYMBY*S4 - 4*S**2*TG**2*(TG+UG)**(-1)*SYMBY - 8*S**2*U1**(-2)
     +    *M2 - 4*S**2*U1**(-1)*M2*(TG+UG)**(-1) - 8*S**2*U1**(-1)*M2*
     +    SYMBY*S4 - 2*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 + 2*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 2*TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 + 2*
     +    TG*M2*(TG+UG)**(-1) )
     +
      M2QBH1 = M2QBH1 + ANG2(30)*N**2*CF * ( 6*TG*M2*SYMBY*S4 - 2*TG*
     +    M2**2*(TG+UG)**(-1)*SYMBY*S4 - 8*TG*M2**2*SYMBY - 2*TG*
     +    (TG+UG)**(-1)*S4 + 2*TG + 6*TG**2*M2*SYMBY + 2*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 2*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 2*
     +    TG**2*SYMBY*S4 - 2*TG**3*M2*(TG+UG)**(-1)*SYMBY + 2*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 - 2*TG**3*SYMBY - 4*T1**(-1)*U1**(-1)*
     +    M2**2*S4 - 8*T1**(-1)*U1**(-1)*M2**3 - 2*T1**(-1)*M2*S4 + 6*
     +    T1**(-1)*M2**2 + 8*U1**(-1)*M2*MS2 - 8*U1**(-1)*M2*SYMBY*
     +    S4**3 + 8*U1**(-1)*M2*S4 + 8*U1**(-1)*M2**2*SYMBY*S4**2 + 8*
     +    U1**(-1)*M2**2 - 4*U1**(-1)*S4**2 + 8*M2*SYMBY*S4**2 - 2*M2
     +     - 8*M2**2*SYMBY*S4 + 2*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(35)*N*CF**2 * ( 8*S**2*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(36)*N*CF**2 * ( 8*S*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(37)*N*CF**2 * ( 16*S*U1**(-1)*M2*S4 - 16*S
     +    *M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(37)*N**2*CF * (  - 8*S**2*U1**(-1)*M2 )
     +
      M2QBH1 = M2QBH1 + ANG2(38)*N*CF**2 * ( 8*S*TG*(TG+UG)**(-1) - 8*S
     +    *T1**(-1)*M2 + 8*S*T1**(-1)*S4 - 16*S*U1**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 8*S*U1**(-1)*M2 + 8*S*U1**(-1)*S4 + 8*S*M2
     +    *(TG+UG)**(-1) + 16*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 8*S**2*
     +    (TG+UG)**(-1) + 4*TG + 16*T1**(-1)*U1**(-1)*M2**2*S4 + 16*
     +    T1**(-1)*U1**(-1)*M2**3 - 16*T1**(-1)*M2**2 + 8*U1**(-1)*M2*
     +    S4 - 24*U1**(-1)*M2**2 - 8*U1**(-1)*S4**2 + 4*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(38)*N**2*CF * (  - 2*S*TG*(TG+UG)**(-1) + 
     +    4*S*T1**(-1)*M2 - 4*S*T1**(-1)*S4 + 8*S*U1**(-2)*M2*S4 + 4*S*
     +    U1**(-1)*M2*(TG+UG)**(-1)*S4 - 4*S*U1**(-1)*M2 - 4*S*U1**(-1)
     +    *S4 - 2*S*M2*(TG+UG)**(-1) + 2*S - 8*S**2*U1**(-2)*M2 - 4*
     +    S**2*U1**(-1)*M2*(TG+UG)**(-1) - 2*S**2*(TG+UG)**(-1) - 8*
     +    T1**(-1)*U1**(-1)*M2**2*S4 - 8*T1**(-1)*U1**(-1)*M2**3 + 8*
     +    T1**(-1)*M2**2 + 8*U1**(-1)*M2**2 )
     +
      M2QBH1 = M2QBH1 + ANG2(42)*N*CF**2 * ( 16*S*U1**(-1)*M2*
     +    (TG+UG)**(-1)*S4 + 16*S*U1**(-1)*M2 + 16*S*MS2*(TG+UG)**(-1)
     +     - 16*S**2*U1**(-1)*M2*(TG+UG)**(-1) - 16*T1**(-1)*U1**(-1)*
     +    M2**2*S4 - 16*T1**(-1)*U1**(-1)*M2**3 - 16*T1**(-1)*M2*MS2 - 
     +    8*T1**(-1)*M2**2 - 16*T1**(-1)*MS2*S4 - 8*T1**(-1)*S4**2 + 16
     +    *U1**(-1)*M2**2 + 16*M2 - 8*M2**2*(TG+UG)**(-1) + 16*MS2 - 8*
     +    (TG+UG)**(-1)*S4**2 )
     +
      M2QBH1 = M2QBH1 + ANG2(42)*N**2*CF * (  - 8*S*U1**(-2)*M2*S4 - 4*
     +    S*U1**(-1)*M2*(TG+UG)**(-1)*S4 - 4*S*MS2*(TG+UG)**(-1) + 8*
     +    S**2*U1**(-2)*M2 + 4*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 8*
     +    T1**(-1)*U1**(-1)*M2**2*S4 + 8*T1**(-1)*U1**(-1)*M2**3 + 8*
     +    T1**(-1)*M2*MS2 + 4*T1**(-1)*M2**2 + 8*T1**(-1)*MS2*S4 + 4*
     +    T1**(-1)*S4**2 - 4*U1**(-1)*M2*S4 - 4*U1**(-1)*M2**2 + 4*
     +    U1**(-1)*S4**2 - 4*M2 + 2*M2**2*(TG+UG)**(-1) - 4*MS2 + 2*
     +    (TG+UG)**(-1)*S4**2 )
     +
      M2QBH1 = M2QBH1 + ANG2(51)*N*CF**2 * ( 12 + 8*S*U1**(-1) + 8*S*
     +    (TG+UG)**(-1) - 8*T1**(-1)*M2 - 8*T1**(-1)*MS2 - 8*U1**(-1)*
     +    M2 - 16*U1**(-1)*MS2 - 8*U1**(-1)*S4 - 8*M2*(TG+UG)**(-1) - 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(51)*N**2*CF * (  - 4 - 4*S*U1**(-1) - 2*S*
     +    (TG+UG)**(-1) + 4*T1**(-1)*M2 + 4*T1**(-1)*MS2 + 4*U1**(-1)*
     +    M2 + 8*U1**(-1)*MS2 + 4*U1**(-1)*S4 + 2*M2*(TG+UG)**(-1) + 2*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(58)*N*CF**2 * (  - 24*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 8*S*TG*M2*SYMBY - 8*S*TG*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 24*S*TG*(TG+UG)**(-1) + 8*S*TG*
     +    SYMBY*S4 + 16*S*TG**2*M2*(TG+UG)**(-1)*SYMBY + 32*S*TG**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 16*S*TG**2*SYMBY - 16*S*TG**3*
     +    (TG+UG)**(-1)*SYMBY + 32*S*U1**(-1)*M2*SYMBY*S4**2 - 16*S*
     +    U1**(-1)*M2**2*SYMBY*S4 + 8*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 
     +    24*S*M2*SYMBY*S4 - 8*S*M2**2*(TG+UG)**(-1)*SYMBY*S4 + 8*S*
     +    M2**2*SYMBY - 8*S*(TG+UG)**(-1)*S4 + 8*S + 8*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBY + 8*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 - 16*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBY - 16*S**2*U1**(-1)*M2*SYMBY*S4
     +     - 8*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 + 8*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 8*S**2*(TG+UG)**(-1) + 8*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 8*TG*M2*(TG+UG)**(-1) + 8*TG*M2*
     +    SYMBY*S4 - 8*TG*M2**2*(TG+UG)**(-1)*SYMBY*S4 - 8*TG*M2**2*
     +    SYMBY )
     +
      M2QBH1 = M2QBH1 + ANG2(58)*N*CF**2 * (  - 8*TG*(TG+UG)**(-1)*S4
     +     + 8*TG + 8*TG**2*M2*SYMBY + 8*TG**2*M2**2*(TG+UG)**(-1)*
     +    SYMBY - 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 8*TG**2*SYMBY*S4
     +     - 8*TG**3*M2*(TG+UG)**(-1)*SYMBY + 8*TG**3*(TG+UG)**(-1)*
     +    SYMBY*S4 - 16*U1**(-1)*M2*SYMBY*S4**3 + 8*U1**(-1)*M2*S4 + 16
     +    *U1**(-1)*M2**2*SYMBY*S4**2 - 8*U1**(-1)*S4**2 + 16*M2*SYMBY*
     +    S4**2 - 8*M2 - 16*M2**2*SYMBY*S4 + 8*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(58)*N**2*CF * ( 6*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 4*S*TG*M2*SYMBY + 2*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 6*S*TG*(TG+UG)**(-1) - 2*S*TG*SYMBY*S4 - 4*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 8*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 6*S*
     +    TG**2*SYMBY + 4*S*TG**3*(TG+UG)**(-1)*SYMBY - 16*S*U1**(-1)*
     +    M2*SYMBY*S4**2 + 8*S*U1**(-1)*M2**2*SYMBY*S4 - 2*S*M2*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 10*S*M2*SYMBY*S4 + 2*S*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 2*S*M2**2*SYMBY + 2*S*(TG+UG)**(-1)*
     +    S4 - 4*S - 2*S**2*TG*M2*(TG+UG)**(-1)*SYMBY - 2*S**2*TG*
     +    (TG+UG)**(-1)*SYMBY*S4 + 4*S**2*TG**2*(TG+UG)**(-1)*SYMBY + 8
     +    *S**2*U1**(-1)*M2*SYMBY*S4 + 2*S**2*M2*(TG+UG)**(-1)*SYMBY*S4
     +     - 2*S**2*M2**2*(TG+UG)**(-1)*SYMBY - 2*S**2*(TG+UG)**(-1) - 
     +    2*TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 2*TG*M2*(TG+UG)**(-1) - 6
     +    *TG*M2*SYMBY*S4 + 2*TG*M2**2*(TG+UG)**(-1)*SYMBY*S4 + 8*TG*
     +    M2**2*SYMBY + 2*TG*(TG+UG)**(-1)*S4 - 4*TG - 6*TG**2*M2*SYMBY
     +     - 2*TG**2*M2**2*(TG+UG)**(-1)*SYMBY )
     +
      M2QBH1 = M2QBH1 + ANG2(58)*N**2*CF * ( 2*TG**2*(TG+UG)**(-1)*
     +    SYMBY*S4**2 - 2*TG**2*SYMBY*S4 + 2*TG**3*M2*(TG+UG)**(-1)*
     +    SYMBY - 2*TG**3*(TG+UG)**(-1)*SYMBY*S4 + 2*TG**3*SYMBY + 8*
     +    U1**(-1)*M2*SYMBY*S4**3 - 4*U1**(-1)*M2*S4 - 8*U1**(-1)*M2**2
     +    *SYMBY*S4**2 + 4*U1**(-1)*S4**2 - 8*M2*SYMBY*S4**2 + 4*M2 + 8
     +    *M2**2*SYMBY*S4 - 4*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(73)*N*CF**2 * ( 24*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 8*S*TG*M2*SYMBY + 8*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 16*S*TG*(TG+UG)**(-1) - 8*S*TG*SYMBY*S4 - 16*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 32*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 16*
     +    S*TG**2*SYMBY + 16*S*TG**3*(TG+UG)**(-1)*SYMBY - 32*S*
     +    U1**(-1)*M2*SYMBY*S4**2 + 8*S*U1**(-1)*M2 + 16*S*U1**(-1)*
     +    M2**2*SYMBY*S4 - 16*S*U1**(-1)*MS2 - 8*S*M2*(TG+UG)**(-1)*
     +    SYMBY*S4**2 + 8*S*M2*(TG+UG)**(-1) + 24*S*M2*SYMBY*S4 + 8*S*
     +    M2**2*(TG+UG)**(-1)*SYMBY*S4 - 8*S*M2**2*SYMBY + 16*S*MS2*
     +    (TG+UG)**(-1) + 8*S*(TG+UG)**(-1)*S4 + 4*S - 8*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBY - 8*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 + 16*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBY + 16*S**2*U1**(-1)*M2*SYMBY*S4
     +     + 8*S**2*U1**(-1) + 8*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 - 8*
     +    S**2*M2**2*(TG+UG)**(-1)*SYMBY - 8*TG*M2*(TG+UG)**(-1)*SYMBY*
     +    S4**2 - 8*TG*M2*(TG+UG)**(-1) - 8*TG*M2*SYMBY*S4 + 8*TG*M2**2
     +    *(TG+UG)**(-1)*SYMBY*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(73)*N*CF**2 * ( 8*TG*M2**2*SYMBY + 8*TG*
     +    (TG+UG)**(-1)*S4 - 8*TG**2*M2*SYMBY - 8*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*
     +    TG**2*SYMBY*S4 + 8*TG**3*M2*(TG+UG)**(-1)*SYMBY - 8*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 - 24*T1**(-1)*M2*MS2 - 24*T1**(-1)*
     +    M2**2 - 16*U1**(-1)*M2*MS2 + 16*U1**(-1)*M2*SYMBY*S4**3 - 16*
     +    U1**(-1)*M2*S4 - 16*U1**(-1)*M2**2*SYMBY*S4**2 - 16*U1**(-1)*
     +    M2**2 - 16*M2*SYMBY*S4**2 + 20*M2 - 8*M2**2*(TG+UG)**(-1) + 
     +    16*M2**2*SYMBY*S4 + 8*MS2 - 8*(TG+UG)**(-1)*S4**2 + 8*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(73)*N**2*CF * (  - 6*S*TG*M2*(TG+UG)**(-1)
     +    *SYMBY*S4 + 4*S*TG*M2*SYMBY - 2*S*TG*(TG+UG)**(-1)*SYMBY*
     +    S4**2 + 4*S*TG*(TG+UG)**(-1) + 2*S*TG*SYMBY*S4 + 4*S*TG**2*M2
     +    *(TG+UG)**(-1)*SYMBY + 8*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 - 6*S
     +    *TG**2*SYMBY - 4*S*TG**3*(TG+UG)**(-1)*SYMBY + 16*S*U1**(-1)*
     +    M2*SYMBY*S4**2 - 4*S*U1**(-1)*M2 - 8*S*U1**(-1)*M2**2*SYMBY*
     +    S4 + 8*S*U1**(-1)*MS2 + 2*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 2*
     +    S*M2*(TG+UG)**(-1) - 10*S*M2*SYMBY*S4 - 2*S*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 2*S*M2**2*SYMBY - 4*S*MS2*
     +    (TG+UG)**(-1) - 2*S*(TG+UG)**(-1)*S4 + 2*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBY + 2*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 - 4*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBY - 8*S**2*U1**(-1)*M2*SYMBY*S4
     +     - 4*S**2*U1**(-1) - 2*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 + 2*
     +    S**2*M2**2*(TG+UG)**(-1)*SYMBY + 2*TG*M2*(TG+UG)**(-1)*SYMBY*
     +    S4**2 + 2*TG*M2*(TG+UG)**(-1) + 6*TG*M2*SYMBY*S4 - 2*TG*M2**2
     +    *(TG+UG)**(-1)*SYMBY*S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(73)*N**2*CF * (  - 8*TG*M2**2*SYMBY - 2*TG
     +    *(TG+UG)**(-1)*S4 + 2*TG + 6*TG**2*M2*SYMBY + 2*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 2*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 2*
     +    TG**2*SYMBY*S4 - 2*TG**3*M2*(TG+UG)**(-1)*SYMBY + 2*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 - 2*TG**3*SYMBY + 12*T1**(-1)*M2*MS2
     +     + 12*T1**(-1)*M2**2 + 8*U1**(-1)*M2*MS2 - 8*U1**(-1)*M2*
     +    SYMBY*S4**3 + 8*U1**(-1)*M2*S4 + 8*U1**(-1)*M2**2*SYMBY*S4**2
     +     + 8*U1**(-1)*M2**2 + 8*M2*SYMBY*S4**2 - 6*M2 + 2*M2**2*
     +    (TG+UG)**(-1) - 8*M2**2*SYMBY*S4 + 2*(TG+UG)**(-1)*S4**2 - 2*
     +    S4 )
     +
      M2QBH1 = M2QBH1 + ANG2(74)*N*CF**2 * (  - 16*T1**(-1)*M2**2*MS2
     +     - 16*T1**(-1)*M2**3 + 8*M2*MS2 )
     +
      M2QBH1 = M2QBH1 + ANG2(74)*N**2*CF * ( 4*TG*M2 + 4*TG**2*M2**2*
     +    SYMBY - 4*TG**3*M2*SYMBY + 8*T1**(-1)*M2**2*MS2 + 8*T1**(-1)*
     +    M2**3 + 4*M2**2 )
     +
      M2QBH1 = M2QBH1 + COLO2(9)*N*CF**2*(S4+MS2) * (  - 8*S**(-1)*TG*
     +    (S+TG)**(-1)*S4**(-1) - 32*S**(-1)*U1**(-4)*M2**2*S4 + 32*
     +    S**(-1)*U1**(-3)*M2*S4 + 32*S**(-1)*U1**(-3)*M2**2 + 16*
     +    S**(-1)*U1**(-2)*M2*MS2*(S+TG)**(-1) - 16*S**(-1)*U1**(-2)*M2
     +    *MS2*S4**(-1) - 32*S**(-1)*U1**(-2)*M2 - 16*S**(-1)*U1**(-2)*
     +    M2**2*S4**(-1) - 16*S**(-1)*U1**(-2)*S4 - 16*S**(-1)*U1**(-1)
     +    *M2*MS2*(S+TG)**(-1)*S4**(-1) + 16*S**(-1)*U1**(-1)*M2*
     +    S4**(-1) + 16*S**(-1)*U1**(-1) - 32*U1**(-4)*M2*MS2 - 32*
     +    U1**(-3)*M2 - 16*U1**(-2)*M2*MS2*(S+TG)**(-1)*S4**(-1) + 16*
     +    U1**(-2) - 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) )
     +
      M2QBH1 = M2QBH1 + COLO2(9)*N**2*CF*(S4+MS2) * (  - 8*S*TG**(-1)*
     +    (S+TG)**(-1)*S4**(-1) - 32*TG**(-2)*U1**(-4)*M2*MS2*S4**2 - 
     +    32*TG**(-2)*U1**(-4)*M2**2*S4**2 + 64*TG**(-2)*U1**(-3)*M2*
     +    MS2*S4 + 64*TG**(-2)*U1**(-3)*M2**2*S4 - 48*TG**(-2)*U1**(-2)
     +    *M2*MS2 - 48*TG**(-2)*U1**(-2)*M2**2 + 16*TG**(-2)*U1**(-1)*
     +    M2*MS2*S4**(-1) + 16*TG**(-2)*U1**(-1)*M2**2*S4**(-1) + 32*
     +    TG**(-1)*U1**(-4)*M2*MS2*S4 - 32*TG**(-1)*U1**(-3)*M2*MS2 + 
     +    32*TG**(-1)*U1**(-3)*M2*S4 + 16*TG**(-1)*U1**(-2)*M2*MS2*
     +    S4**(-1) - 32*TG**(-1)*U1**(-2)*M2 - 16*TG**(-1)*U1**(-2)*S4
     +     + 16*TG**(-1)*U1**(-1)*M2*S4**(-1) + 16*TG**(-1)*U1**(-1) - 
     +    8*(S+TG)**(-1)*S4**(-1) )

      IF (IFL.EQ.1) THEN

      M2QBH2 = 0D0
      M2QBH2 = M2QBH2 + N*CF**2*S4G2**(-1)*S4G * (  - 32*S**(-1)*
     +    T1**(-1)*M2**2 - 16*S**(-1)*U1**(-1)*M2**2 + 48*S**(-1)*M2 - 
     +    48*T1**(-1)*U1**(-1)*M2*MS2 - 48*T1**(-1)*U1**(-1)*M2**2 + 16
     +    *T1**(-1)*M2 + 16*T1**(-1)*MS2 + 8*U1**(-2)*M2**2 + 24*
     +    U1**(-1)*M2 - 16*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + N*CF**2 * ( 8*S**(-2)*TG*(TG+UG)**(-1)*S4 - 8*
     +    S**(-2)*TG - 16*S**(-2)*M2*(TG+UG)**(-1)*S4 + 8*S**(-2)*M2 - 
     +    8*S**(-2)*M2**2*(TG+UG)**(-1) - 8*S**(-2)*(TG+UG)**(-1)*S4**2
     +     + 8*S**(-2)*S4 - 16*S**(-1)*MS2*(TG+UG)**(-1) - 8*S**(-1)*
     +    (TG+UG)**(-1)*S4 - 16*S**(-1) - 16*T1**(-1)*U1**(-1)*M2 - 16*
     +    T1**(-1)*U1**(-1)*MS2 + 16*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + N**2*CF*S4G2**(-1)*S4G * (  - 20*S**(-2)*TG*M2
     +     + 16*S**(-1)*T1**(-1)*M2**2 + 12*S**(-1)*U1**(-1)*M2**2 - 12
     +    *S**(-1)*M2 + 8*S**(-1)*MS2 + 24*T1**(-1)*U1**(-1)*M2*MS2 + 
     +    24*T1**(-1)*U1**(-1)*M2**2 - 4*U1**(-1)*M2 + 8*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + N**2*CF * (  - 4*S**(-2)*TG*(TG+UG)**(-1)*S4 + 
     +    4*S**(-2)*TG + 8*S**(-2)*M2*(TG+UG)**(-1)*S4 + 4*S**(-2)*M2
     +     + 4*S**(-2)*M2**2*(TG+UG)**(-1) + 4*S**(-2)*(TG+UG)**(-1)*
     +    S4**2 - 4*S**(-2)*S4 + 8*S**(-1)*MS2*(TG+UG)**(-1) + 4*
     +    S**(-1)*(TG+UG)**(-1)*S4 + 8*S**(-1) + 8*T1**(-1)*U1**(-1)*M2
     +     + 8*T1**(-1)*U1**(-1)*MS2 - 8*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(13)*N*CF**2 * (  - 16*S**(-2)*TG*M2*S4 + 8
     +    *S**(-2)*TG**2*M2 + 8*S**(-2)*M2*S4**2 + 8*S**(-1)*TG*M2 - 8*
     +    S**(-1)*M2*S4 + 4*M2 )
     +
      M2QBH2 = M2QBH2 + ANG2(14)*N*CF**2*S4G2**(-1)*S4G * (  - 32*
     +    S**(-1)*TG*M2 - 32*S**(-1)*T1**(-1)*M2**3 + 32*S**(-1)*M2**2
     +     - 32*T1**(-1)*M2*MS2 - 16*T1**(-1)*M2**2 + 16*U1**(-1)*M2**2
     +     - 16*M2 )
     +
      M2QBH2 = M2QBH2 + ANG2(14)*N*CF**2 * ( 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 - 12*S**(-2)*TG*M2 + 4*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) + 4*S**(-2)*TG*(TG+UG)**(-1)*S4**2 - 4*S**(-2)*
     +    TG*S4 - 8*S**(-2)*M2*(TG+UG)**(-1)*S4**2 + 12*S**(-2)*M2*S4
     +     - 4*S**(-2)*M2**2*(TG+UG)**(-1)*S4 - 4*S**(-2)*(TG+UG)**(-1)
     +    *S4**3 + 4*S**(-2)*S4**2 - 16*S**(-1)*T1**(-1)*M2**2 + 16*
     +    S**(-1)*M2*(TG+UG)**(-1)*S4 + 16*S**(-1)*M2 + 4*S**(-1)*M2**2
     +    *(TG+UG)**(-1) + 4*S**(-1)*(TG+UG)**(-1)*S4**2 - 8*S**(-1)*S4
     +     - 16*T1**(-1)*M2 - 16*T1**(-1)*MS2 + 16*U1**(-1)*M2 + 8*MS2*
     +    (TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(14)*N**2*CF*S4G2**(-1)*S4G * (  - 8*
     +    S**(-2)*TG*M2**2 + 8*S**(-2)*TG**2*M2 + 20*S**(-1)*TG*M2 + 16
     +    *S**(-1)*T1**(-1)*M2**3 - 12*S**(-1)*M2**2 + 16*T1**(-1)*M2*
     +    MS2 + 8*T1**(-1)*M2**2 - 4*U1**(-1)*M2**2 + 8*M2 )
     +
      M2QBH2 = M2QBH2 + ANG2(14)*N**2*CF * (  - 4*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 - 10*S**(-2)*TG*M2 - 2*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) - 2*S**(-2)*TG*(TG+UG)**(-1)*S4**2 + 2*S**(-2)*
     +    TG*S4 + 4*S**(-2)*M2*(TG+UG)**(-1)*S4**2 + 2*S**(-2)*M2*S4 + 
     +    2*S**(-2)*M2**2*(TG+UG)**(-1)*S4 + 2*S**(-2)*(TG+UG)**(-1)*
     +    S4**3 - 2*S**(-2)*S4**2 + 8*S**(-1)*T1**(-1)*M2**2 - 8*
     +    S**(-1)*M2*(TG+UG)**(-1)*S4 - 8*S**(-1)*M2 - 2*S**(-1)*M2**2*
     +    (TG+UG)**(-1) - 2*S**(-1)*(TG+UG)**(-1)*S4**2 + 4*S**(-1)*S4
     +     + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 - 4*U1**(-1)*M2 - 4*MS2*
     +    (TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(41)*N*CF**2*S4G2**(-1)*S4G * (  - 16*
     +    S**(-1)*T1**(-1)*M2**2 - 16*S**(-1)*U1**(-1)*M2**2 + 32*
     +    S**(-1)*M2 - 16*T1**(-1)*M2 - 16*T1**(-1)*MS2 + 16*U1**(-1)*
     +    M2 - 16*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(41)*N*CF**2 * ( 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) + 8*S**(-2)*TG*(TG+UG)**(-1)*S4 - 8*S**(-2)*TG
     +     - 16*S**(-2)*M2*(TG+UG)**(-1)*S4 + 12*S**(-2)*M2 - 4*S**(-2)
     +    *M2**2*(TG+UG)**(-1) - 12*S**(-2)*(TG+UG)**(-1)*S4**2 + 12*
     +    S**(-2)*S4 + 8*S**(-1)*T1**(-1)*M2 + 8*S**(-1)*T1**(-1)*S4 - 
     +    8*S**(-1)*U1**(-1)*M2 - 8*S**(-1)*U1**(-1)*S4 + 12*S**(-1)*M2
     +    *(TG+UG)**(-1) + 12*S**(-1)*(TG+UG)**(-1)*S4 - 20*S**(-1) - 8
     +    *T1**(-1) + 8*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(41)*N**2*CF*S4G2**(-1)*S4G * (  - 8*
     +    S**(-2)*TG*M2 + 8*S**(-1)*T1**(-1)*M2**2 + 4*S**(-1)*U1**(-1)
     +    *M2**2 - 8*S**(-1)*M2 + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 - 4*
     +    U1**(-1)*M2 + 4*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(41)*N**2*CF * (  - 4*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) - 4*S**(-2)*TG*(TG+UG)**(-1)*S4 + 4*S**(-2)*TG
     +     + 8*S**(-2)*M2*(TG+UG)**(-1)*S4 + 2*S**(-2)*M2 + 2*S**(-2)*
     +    M2**2*(TG+UG)**(-1) + 6*S**(-2)*(TG+UG)**(-1)*S4**2 - 6*
     +    S**(-2)*S4 - 4*S**(-1)*T1**(-1)*M2 - 4*S**(-1)*T1**(-1)*S4 + 
     +    2*S**(-1)*U1**(-1)*M2 + 2*S**(-1)*U1**(-1)*S4 - 6*S**(-1)*M2*
     +    (TG+UG)**(-1) - 6*S**(-1)*(TG+UG)**(-1)*S4 + 12*S**(-1) + 4*
     +    T1**(-1) - 2*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(43)*N*CF**2 * (  - 16*S**(-2)*TG*M2 + 16*
     +    S**(-2)*M2*S4 - 8*S**(-1)*M2 )
     +
      M2QBH2 = M2QBH2 + ANG2(45)*N*CF**2 * ( 8*S**(-2)*M2 )
     +
      M2QBH2 = M2QBH2 + ANG2(46)*N*CF**2*S4G2**(-1)*S4G * (  - 8*
     +    T1**(-1) + 8*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(46)*N*CF**2 * (  - 8*S**(-2)*M2*
     +    (TG+UG)**(-1) - 8*S**(-2)*(TG+UG)**(-1)*S4 + 8*S**(-2) + 8*
     +    S**(-1)*T1**(-1) - 8*S**(-1)*U1**(-1) + 8*S**(-1)*
     +    (TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(46)*N**2*CF*S4G2**(-1)*S4G * ( 4*S**(-1)
     +     + 4*T1**(-1) - 2*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(46)*N**2*CF * ( 4*S**(-2)*M2*(TG+UG)**(-1)
     +     + 4*S**(-2)*(TG+UG)**(-1)*S4 - 4*S**(-2) - 4*S**(-1)*
     +    T1**(-1) + 2*S**(-1)*U1**(-1) - 4*S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(49)*N*CF**2*S4G2**(-1)*S4G * ( 8*S**(-1)*
     +    T1**(-1)*M2**2 - 8*S**(-1)*U1**(-1)*M2**2 - 16*T1**(-1)*M2 - 
     +    16*T1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(49)*N*CF**2 * ( 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) + 8*S**(-2)*TG*(TG+UG)**(-1)*S4 - 8*S**(-2)*TG
     +     - 8*S**(-1)*U1**(-1)*M2 + 12*S**(-1)*M2*(TG+UG)**(-1) + 8*
     +    S**(-1)*MS2*(TG+UG)**(-1) + 8*S**(-1)*(TG+UG)**(-1)*S4 - 4*
     +    S**(-1) - 8*T1**(-1) - 8*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(49)*N**2*CF*S4G2**(-1)*S4G * (  - 2*
     +    S**(-1)*T1**(-1)*M2**2 + 4*S**(-1)*U1**(-1)*M2**2 - 4*S**(-1)
     +    *M2 - 4*S**(-1)*MS2 + 4*T1**(-1)*M2 + 4*T1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(49)*N**2*CF * (  - 4*S**(-2)*TG*M2*
     +    (TG+UG)**(-1) - 4*S**(-2)*TG*(TG+UG)**(-1)*S4 + 4*S**(-2)*TG
     +     + 4*S**(-1)*U1**(-1)*M2 - 6*S**(-1)*M2*(TG+UG)**(-1) - 4*
     +    S**(-1)*MS2*(TG+UG)**(-1) - 4*S**(-1)*(TG+UG)**(-1)*S4 + 2*
     +    T1**(-1) + 4*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(50)*N*CF**2*S4G2**(-1)*S4G * (  - 8*
     +    S**(-1)*T1**(-1)*M2 + 8*S**(-1) - 4*T1**(-2)*M2 + 4*T1**(-1)
     +     + 8*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(50)*N*CF**2 * (  - 4*S**(-2)*TG*
     +    (TG+UG)**(-1) - 8*S**(-2)*M2*(TG+UG)**(-1) - 8*S**(-2)*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2) + 8*S**(-1)*T1**(-1) - 8*S**(-1)
     +    *U1**(-1) + 8*S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(50)*N**2*CF*S4G2**(-1)*S4G * ( 2*S**(-2)*
     +    TG + 2*S**(-1) - 2*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(50)*N**2*CF * ( 2*S**(-2)*TG*(TG+UG)**(-1)
     +     + 4*S**(-2)*M2*(TG+UG)**(-1) + 4*S**(-2)*(TG+UG)**(-1)*S4 - 
     +    4*S**(-2) - 4*S**(-1)*T1**(-1) + 2*S**(-1)*U1**(-1) - 4*
     +    S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(51)*N*CF**2*S4G2**(-1)*S4G * (  - 8*
     +    S**(-1)*U1**(-1)*M2**3 + 8*S**(-1)*M2**2 - 16*T1**(-1)*M2*MS2
     +     - 8*U1**(-1)*M2**2 )
     +
      M2QBH2 = M2QBH2 + ANG2(51)*N*CF**2 * ( 8 + 8*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 - 8*S**(-2)*TG*M2 + 4*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) - 8*S**(-1)*U1**(-1)*M2**2 + 16*S**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 12*S**(-1)*M2 + 4*S**(-1)*M2**2*
     +    (TG+UG)**(-1) + 8*S**(-1)*MS2*(TG+UG)**(-1)*S4 - 8*S**(-1)*
     +    MS2 + 8*S**(-1)*(TG+UG)**(-1)*S4**2 - 8*S**(-1)*S4 - 8*
     +    U1**(-1)*M2 + 8*MS2*(TG+UG)**(-1) - 8*(TG+UG)**(-1)*S4 )
     +
      M2QBH2 = M2QBH2 + ANG2(51)*N**2*CF*S4G2**(-1)*S4G * ( 2*S**(-2)*
     +    TG*M2**2 + 4*S**(-1)*U1**(-1)*M2**3 + 4*S**(-1)*M2*MS2 + 4*
     +    T1**(-1)*M2*MS2 + 4*U1**(-1)*M2**2 )
     +
      M2QBH2 = M2QBH2 + ANG2(51)*N**2*CF * (  - 4 - 4*S**(-2)*TG*M2*
     +    (TG+UG)**(-1)*S4 + 4*S**(-2)*TG*M2 - 2*S**(-2)*TG*M2**2*
     +    (TG+UG)**(-1) + 4*S**(-1)*U1**(-1)*M2**2 - 8*S**(-1)*M2*
     +    (TG+UG)**(-1)*S4 + 8*S**(-1)*M2 - 2*S**(-1)*M2**2*
     +    (TG+UG)**(-1) - 4*S**(-1)*MS2*(TG+UG)**(-1)*S4 + 4*S**(-1)*
     +    MS2 - 4*S**(-1)*(TG+UG)**(-1)*S4**2 + 4*S**(-1)*S4 + 4*
     +    U1**(-1)*M2 - 4*MS2*(TG+UG)**(-1) + 4*(TG+UG)**(-1)*S4 )
     +
      M2QBH2 = M2QBH2 + ANG2(66)*N*CF**2*S4G2**(-1)*S4G * ( 8*S**(-1)*
     +    T1**(-1)*M2 - 8*S**(-1) - 8*T1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(66)*N*CF**2 * ( 4*S**(-2)*TG*(TG+UG)**(-1)
     +     )
     +
      M2QBH2 = M2QBH2 + ANG2(66)*N**2*CF*S4G2**(-1)*S4G * (  - 2*
     +    S**(-2)*TG - 2*S**(-1)*T1**(-1)*M2 + 2*T1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(66)*N**2*CF * (  - 2*S**(-2)*TG*
     +    (TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(69)*N*CF**2*S4G2**(-1)*S4G * ( 8*S**(-1)*
     +    T1**(-1)*M2**2 + 8*S**(-1)*U1**(-1)*M2**2 - 16*S**(-1)*M2 - 8
     +    *T1**(-1)*M2 + 8*U1**(-1)*M2 + 16*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(69)*N*CF**2 * (  - 8*S**(-2)*TG*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2)*TG - 8*S**(-2)*M2*(TG+UG)**(-1)*
     +    S4 + 8*S**(-2)*(TG+UG)**(-1)*S4**2 - 8*S**(-2)*S4 + 8*S**(-1)
     +    *T1**(-1)*M2 + 8*S**(-1)*M2*(TG+UG)**(-1) + 8*S**(-1)*MS2*
     +    (TG+UG)**(-1) - 12*S**(-1)*(TG+UG)**(-1)*S4 + 12*S**(-1) + 8*
     +    T1**(-1) + 4*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(69)*N**2*CF*S4G2**(-1)*S4G * (  - 4*
     +    S**(-2)*TG*M2 - 2*S**(-1)*T1**(-1)*M2**2 - 4*S**(-1)*U1**(-1)
     +    *M2**2 + 4*S**(-1)*M2 - 4*S**(-1)*MS2 + 2*T1**(-1)*M2 - 4*
     +    U1**(-1)*M2 - 8*U1**(-1)*MS2 )
     +
      M2QBH2 = M2QBH2 + ANG2(69)*N**2*CF * ( 4*S**(-2)*TG*(TG+UG)**(-1)
     +    *S4 - 4*S**(-2)*TG + 4*S**(-2)*M2*(TG+UG)**(-1)*S4 - 4*
     +    S**(-2)*(TG+UG)**(-1)*S4**2 + 4*S**(-2)*S4 - 2*S**(-1)*
     +    T1**(-1)*M2 - 4*S**(-1)*M2*(TG+UG)**(-1) - 4*S**(-1)*MS2*
     +    (TG+UG)**(-1) + 6*S**(-1)*(TG+UG)**(-1)*S4 - 6*S**(-1) - 2*
     +    T1**(-1) - 2*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(70)*N*CF**2*S4G2**(-1)*S4G * ( 8*S**(-1)*
     +    T1**(-1)*M2 - 8*S**(-1) - 4*U1**(-2)*M2 - 4*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(70)*N*CF**2 * ( 4*S**(-2)*TG*(TG+UG)**(-1)
     +     + 8*S**(-2)*M2*(TG+UG)**(-1) + 8*S**(-2)*(TG+UG)**(-1)*S4 - 
     +    8*S**(-2) - 8*S**(-1)*T1**(-1) + 8*S**(-1)*U1**(-1) - 8*
     +    S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(70)*N**2*CF*S4G2**(-1)*S4G * (  - 2*
     +    S**(-2)*TG - 2*S**(-1)*T1**(-1)*M2 - 2*S**(-1)*U1**(-1)*M2 + 
     +    2*S**(-1) + 2*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(70)*N**2*CF * (  - 2*S**(-2)*TG*
     +    (TG+UG)**(-1) - 4*S**(-2)*M2*(TG+UG)**(-1) - 4*S**(-2)*
     +    (TG+UG)**(-1)*S4 + 4*S**(-2) + 2*S**(-1)*T1**(-1) - 4*S**(-1)
     +    *U1**(-1) + 4*S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(72)*N*CF**2*S4G2**(-1)*S4G * (  - 8*
     +    S**(-1)*T1**(-1)*M2 + 8*S**(-1) + 8*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(72)*N*CF**2 * (  - 4*S**(-2)*TG*
     +    (TG+UG)**(-1) - 8*S**(-2)*M2*(TG+UG)**(-1) - 8*S**(-2)*
     +    (TG+UG)**(-1)*S4 + 8*S**(-2) + 8*S**(-1)*T1**(-1) - 8*S**(-1)
     +    *U1**(-1) + 8*S**(-1)*(TG+UG)**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(72)*N**2*CF*S4G2**(-1)*S4G * ( 2*S**(-2)*
     +    TG + 2*S**(-1)*T1**(-1)*M2 - 4*S**(-1) - 4*U1**(-1) )
     +
      M2QBH2 = M2QBH2 + ANG2(72)*N**2*CF * ( 2*S**(-2)*TG*(TG+UG)**(-1)
     +     + 4*S**(-2)*M2*(TG+UG)**(-1) + 4*S**(-2)*(TG+UG)**(-1)*S4 - 
     +    4*S**(-2) - 2*S**(-1)*T1**(-1) + 4*S**(-1)*U1**(-1) - 4*
     +    S**(-1)*(TG+UG)**(-1) )

      M2QBH3 = 0D0
      M2QBH3 = M2QBH3 + N*CF*S4G2**(-1)*(S4+MS2)*S4G * ( 4*S**(-1)*TG*
     +    M2*(S+TG)**(-1)*S4**(-1) + 4*S**(-1)*TG*S4**(-1) - 4*S**(-1)*
     +    TG**2*(S+TG)**(-1)*S4**(-1) - 4*S**(-1)*M2*S4**(-1) - 4*S*
     +    (S+TG)**(-1)*S4**(-1) - 8*TG*(S+TG)**(-1)*S4**(-1) - 8*
     +    U1**(-1)*M2*S4**(-1) + 8*U1**(-1)*M2**2*(S+TG)**(-1)*S4**(-1)
     +     - 4*M2*(S+TG)**(-1)*S4**(-1) + 4*S4**(-1) )
     +
      M2QBH3 = M2QBH3 + N*CF*S4G2**(-1)*S4G * ( 4*S**(-1)*TG - 4*
     +    S**(-1)*T1**(-1)*M2*MS2 - 4*S**(-1)*T1**(-1)*M2**2 - 4*
     +    S**(-1)*U1**(-1)*M2*MS2 - 4*S**(-1)*U1**(-1)*M2**2 + 8*
     +    S**(-1)*M2 + 8*S**(-1)*MS2 + 4*S*U1**(-1) - 8*T1**(-1)*
     +    U1**(-1)*M2*MS2 - 8*T1**(-1)*U1**(-1)*M2**2 - 8*T1**(-1)*M2
     +     - 8*T1**(-1)*MS2 - 8*U1**(-1)*M2 - 8*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + N*CF*(S4+MS2) * ( 8*S**(-1)*TG*M2*SYMBT*
     +    S4**(-1) - 4*S**(-1)*TG*(S+TG)**(-1)*S4**(-1) - 4*S**(-1)*TG*
     +    (TG+UG)**(-1)*S4**(-1) - 4*S**(-1)*TG*SYMBT + 4*S**(-1)*TG**2
     +    *(S+TG)**(-2)*S4**(-1) - 4*S**(-1)*TG**2*SYMBT*S4**(-1) + 4*
     +    S**(-1)*U1**(-1)*M2*(S+TG)**(-1) - 4*S**(-1)*U1**(-1)*M2*
     +    S4**(-1) - 4*S**(-1)*U1**(-1)*(S+TG)**(-1)*S4 + 4*S**(-1)*
     +    U1**(-1) - 4*S**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 4*S**(-1)*
     +    (S+TG)**(-1) + 4*S**(-1)*(TG+UG)**(-1) - 8*S**(-1)*S4**(-1)
     +     - 8*S*TG*(TG+UG)**(-1)*SYMBT*S4**(-1) + 8*S*(TG+UG)**(-1)*
     +    SYMBT - 4*S*SYMBT*S4**(-1) - 4*S**2*(TG+UG)**(-1)*SYMBT*
     +    S4**(-1) + 4*TG*(S+TG)**(-2)*S4**(-1) + 8*TG*(TG+UG)**(-1)*
     +    SYMBT - 8*TG*SYMBT*S4**(-1) - 4*TG**2*(TG+UG)**(-1)*SYMBT*
     +    S4**(-1) + 8*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 8*M2*SYMBT*
     +    S4**(-1) - 4*(TG+UG)**(-1)*SYMBT*S4 - 4*(TG+UG)**(-1)*
     +    S4**(-1) + 4*SYMBT )
     +
      M2QBH3 = M2QBH3 + N*CF * (  - 4*S**(-1)*TG*(TG+UG)**(-1) - 16*
     +    S**(-1)*T1**(-1)*M2 - 16*S**(-1)*T1**(-1)*MS2 - 12*S**(-1)*
     +    U1**(-1)*M2 - 4*S**(-1)*M2*(TG+UG)**(-1) - 4*S**(-1)*
     +    (TG+UG)**(-1)*S4 - 8*T1**(-1)*U1**(-1)*M2 - 8*T1**(-1)*
     +    U1**(-1)*MS2 + 4*(TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + CF**2*S4G2**(-1)*S4G * (  - 8 - 8*S*U1**(-1) - 
     +    16*T1**(-2)*M2*MS2 - 16*T1**(-2)*M2**2 + 16*T1**(-1)*U1**(-1)
     +    *M2*MS2 + 16*T1**(-1)*U1**(-1)*M2**2 + 8*T1**(-1)*M2 + 8*
     +    T1**(-1)*MS2 - 16*U1**(-2)*M2*MS2 - 16*U1**(-2)*M2**2 + 16*
     +    U1**(-1)*M2 + 8*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + CF**2*(S4+MS2) * (  - 8*S**(-1)*TG*M2*SYMBT*
     +    S4**(-1) + 8*S**(-1)*TG*(TG+UG)**(-1)*S4**(-1) + 8*S**(-1)*
     +    TG**2*SYMBT*S4**(-1) - 8*S**(-1)*U1**(-1)*M2*(S+TG)**(-1) + 8
     +    *S**(-1)*U1**(-1)*M2*S4**(-1) + 8*S**(-1)*U1**(-1)*
     +    (S+TG)**(-1)*S4 - 8*S**(-1)*U1**(-1) + 8*S**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 8*S**(-1)*(S+TG)**(-1) - 8*S**(-1)*
     +    (TG+UG)**(-1) + 8*S**(-1)*S4**(-1) + 16*S*TG*(TG+UG)**(-1)*
     +    SYMBT*S4**(-1) - 16*S*(TG+UG)**(-1)*SYMBT + 8*S*SYMBT*
     +    S4**(-1) + 8*S**2*(TG+UG)**(-1)*SYMBT*S4**(-1) - 16*TG*
     +    (TG+UG)**(-1)*SYMBT + 16*TG*SYMBT*S4**(-1) + 8*TG**2*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 8*M2*SYMBT*S4**(-1) + 8*
     +    (TG+UG)**(-1)*SYMBT*S4 + 8*(TG+UG)**(-1)*S4**(-1) - 8*SYMBT )
     +
      M2QBH3 = M2QBH3 + CF**2 * ( 8*S**(-1)*TG*(TG+UG)**(-1) + 16*
     +    S**(-1)*T1**(-1)*M2 + 16*S**(-1)*T1**(-1)*MS2 + 8*S**(-1)*
     +    U1**(-1)*M2 + 8*S**(-1)*M2*(TG+UG)**(-1) + 8*S**(-1)*
     +    (TG+UG)**(-1)*S4 - 8*S**(-1) + 16*T1**(-1)*U1**(-1)*M2 + 16*
     +    T1**(-1)*U1**(-1)*MS2 - 16*U1**(-2)*M2 + 8*U1**(-1) - 8*
     +    (TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(13)*CF**2 * ( 8*S**(-1)*TG*M2 - 8*S**(-1)*
     +    M2*S4 - 8*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(14)*N*CF*S4G2**(-1)*S4G * ( 4*S**(-1)*TG*
     +    M2 - 8*S*T1**(-1)*M2 - 8*S*T1**(-1)*MS2 + 4*S*U1**(-1)*M2 - 4
     +    *S*U1**(-1)*MS2 - 4*S - 2*TG - 16*T1**(-1)*M2*MS2 - 16*
     +    T1**(-1)*M2**2 - 4*U1**(-1)*M2**2 + 2*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(14)*N*CF * ( 2*S**(-1)*TG*M2*(TG+UG)**(-1)
     +     + 2*S**(-1)*TG*(TG+UG)**(-1)*S4 - 16*S**(-1)*T1**(-1)*M2*MS2
     +     - 16*S**(-1)*T1**(-1)*M2*S4 - 16*S**(-1)*T1**(-1)*M2**2 - 16
     +    *S**(-1)*T1**(-1)*MS2*S4 - 2*S**(-1)*U1**(-1)*M2**2 - 2*
     +    S**(-1)*U1**(-1)*S4**2 - 2*S**(-1)*M2*(TG+UG)**(-1)*S4 + 8*
     +    S**(-1)*M2 - 4*S**(-1)*M2**2*(TG+UG)**(-1) + 8*S**(-1)*MS2 - 
     +    6*S**(-1)*(TG+UG)**(-1)*S4**2 + 2*S*U1**(-1) + 2*S*
     +    (TG+UG)**(-1) + 2*TG*(TG+UG)**(-1) - 8*T1**(-1)*M2 - 8*
     +    T1**(-1)*MS2 - 4*U1**(-1)*MS2 + 6*M2*(TG+UG)**(-1) + 8*MS2*
     +    (TG+UG)**(-1) + 4*(TG+UG)**(-1)*S4 )
     +
      M2QBH3 = M2QBH3 + ANG2(14)*CF**2*S4G2**(-1)*S4G * ( 8*S*T1**(-1)*
     +    M2 + 8*S*T1**(-1)*MS2 - 8*S*U1**(-1)*M2 + 8*S*U1**(-1)*MS2 + 
     +    4*S + 16*T1**(-1)*M2*MS2 + 16*T1**(-1)*M2**2 + 8*U1**(-1)*
     +    M2**2 - 8*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(14)*CF**2 * (  - 4 - 4*S**(-1)*TG*M2*
     +    (TG+UG)**(-1) - 4*S**(-1)*TG*(TG+UG)**(-1)*S4 + 16*S**(-1)*
     +    T1**(-1)*M2*MS2 + 16*S**(-1)*T1**(-1)*M2*S4 + 16*S**(-1)*
     +    T1**(-1)*M2**2 + 16*S**(-1)*T1**(-1)*MS2*S4 + 4*S**(-1)*
     +    U1**(-1)*M2**2 + 4*S**(-1)*U1**(-1)*S4**2 + 4*S**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 16*S**(-1)*M2 + 8*S**(-1)*M2**2*
     +    (TG+UG)**(-1) - 8*S**(-1)*MS2 + 12*S**(-1)*(TG+UG)**(-1)*
     +    S4**2 - 8*S**(-1)*S4 - 4*S*U1**(-1) - 4*S*(TG+UG)**(-1) - 4*
     +    TG*(TG+UG)**(-1) + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 8*
     +    U1**(-1)*MS2 - 12*M2*(TG+UG)**(-1) - 16*MS2*(TG+UG)**(-1) - 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QBH3 = M2QBH3 + ANG2(16)*N*CF*S4G2**(-1)*S4G * (  - 2*S**(-1)*
     +    TG*M2 - 4*S**(-1)*TG*MS2 - 6*S**(-1)*T1**(-1)*M2**2*MS2 - 6*
     +    S**(-1)*T1**(-1)*M2**3 + 6*S**(-1)*M2*MS2 + 6*S**(-1)*M2**2
     +     + 6*S*U1**(-1)*M2 + 8*S*U1**(-1)*MS2 - 4*S - 4*S**2*U1**(-1)
     +     + 4*TG - 6*T1**(-1)*U1**(-1)*M2**2*MS2 - 6*T1**(-1)*U1**(-1)
     +    *M2**3 + 2*U1**(-1)*M2*MS2 + 4*U1**(-1)*M2**2 + 2*M2 + 8*MS2
     +     )
     +
      M2QBH3 = M2QBH3 + ANG2(16)*N*CF * ( 4*S**(-1)*TG - 8*S**(-1)*
     +    T1**(-1)*M2*MS2 - 6*S**(-1)*T1**(-1)*M2*S4 - 8*S**(-1)*
     +    T1**(-1)*M2**2 - 6*S**(-1)*T1**(-1)*MS2*S4 + 10*S**(-1)*M2 + 
     +    10*S**(-1)*MS2 - 4*S**(-1)*S4 + 4*S*U1**(-1) - 6*T1**(-1)*
     +    U1**(-1)*M2*MS2 - 2*T1**(-1)*U1**(-1)*M2*S4 - 6*T1**(-1)*
     +    U1**(-1)*M2**2 - 2*T1**(-1)*U1**(-1)*MS2*S4 + 2*U1**(-1)*M2
     +     + 2*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(16)*CF**2*S4G2**(-1)*S4G * (  - 4*S*
     +    U1**(-1)*M2 - 12*S*U1**(-1)*MS2 + 4*S + 4*S**2*U1**(-1) - 16*
     +    T1**(-2)*M2**2*MS2 - 16*T1**(-2)*M2**3 + 12*T1**(-1)*U1**(-1)
     +    *M2**2*MS2 + 12*T1**(-1)*U1**(-1)*M2**3 + 4*T1**(-1)*M2*MS2
     +     + 4*T1**(-1)*M2**2 - 8*U1**(-1)*M2*MS2 - 12*U1**(-1)*M2**2
     +     - 4*M2 - 8*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(16)*CF**2 * (  - 4*S**(-1)*TG + 4*S**(-1)*
     +    T1**(-1)*M2*S4 + 4*S**(-1)*T1**(-1)*MS2*S4 - 8*S**(-1)*M2 - 4
     +    *S**(-1)*MS2 + 4*S**(-1)*S4 - 4*S*U1**(-1) - 8*T1**(-2)*M2*
     +    MS2 - 8*T1**(-2)*M2**2 + 12*T1**(-1)*U1**(-1)*M2*MS2 + 4*
     +    T1**(-1)*U1**(-1)*M2*S4 + 12*T1**(-1)*U1**(-1)*M2**2 + 4*
     +    T1**(-1)*U1**(-1)*MS2*S4 - 8*U1**(-1)*M2 - 4*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(25)*N*CF*S4G2**(-1)*S4G * ( 2 - 2*S**(-1)*
     +    TG + 2*S**(-1)*T1**(-1)*M2*MS2 + 2*S**(-1)*T1**(-1)*M2**2 - 2
     +    *S**(-1)*M2 - 2*S**(-1)*MS2 + 2*S*U1**(-1) + 2*T1**(-1)*
     +    U1**(-1)*M2*MS2 + 2*T1**(-1)*U1**(-1)*M2**2 - 4*U1**(-1)*M2
     +     - 2*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(25)*N*CF * ( 2*S**(-1)*T1**(-1)*M2 + 2*
     +    S**(-1)*T1**(-1)*MS2 + 2*S**(-1) + 2*T1**(-1)*U1**(-1)*M2 + 2
     +    *T1**(-1)*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(25)*CF**2*S4G2**(-1)*S4G * (  - 4 - 4*S*
     +    U1**(-1) - 4*T1**(-1)*U1**(-1)*M2*MS2 - 4*T1**(-1)*U1**(-1)*
     +    M2**2 + 8*U1**(-1)*M2 + 4*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(25)*CF**2 * (  - 4*T1**(-1)*U1**(-1)*M2 - 
     +    4*T1**(-1)*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(30)*N*CF*S4G2**(-1)*S4G * ( 2*S*TG - 2*S*
     +    T1**(-1)*M2*MS2 - 2*S*T1**(-1)*M2**2 + 4*S*U1**(-1)*M2**2 + 2
     +    *S*M2 + 2*S*MS2 - 2*TG*MS2 + 2*TG**2 - 4*T1**(-1)*M2**2*MS2
     +     - 4*T1**(-1)*M2**3 - 4*U1**(-1)*M2**3 + 4*M2*MS2 + 8*M2**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(30)*N*CF * (  - 2*S**(-1)*TG*M2*MS2*SYMBY*
     +    S4 - 2*S**(-1)*TG*M2 - 2*S**(-1)*TG*M2**2*SYMBY*S4 + 2*
     +    S**(-1)*TG**2*M2*MS2*SYMBY + 2*S**(-1)*TG**2*M2*SYMBY*S4 - 2*
     +    S**(-1)*TG**2*M2**2*SYMBY + 2*S**(-1)*TG**2*MS2*SYMBY*S4 + 2*
     +    S**(-1)*TG**2 + 2*S**(-1)*TG**3*M2*SYMBY - 2*S**(-1)*TG**3*
     +    MS2*SYMBY + 8*S*U1**(-1)*M2 - 2*TG*M2*MS2*SYMBY - 2*TG*M2**2*
     +    SYMBY + 2*TG + 2*TG**2*M2*SYMBY + 2*TG**2*MS2*SYMBY - 2*
     +    T1**(-1)*M2*MS2 - 2*T1**(-1)*M2**2 - 8*U1**(-1)*M2*S4 - 4*
     +    U1**(-1)*M2**2 + 10*M2 + 2*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(31)*N*CF*S4G2**(-1)*S4G * (  - 2*S**(-1)*
     +    U1**(-1)*M2*MS2 - 2*S**(-1)*U1**(-1)*M2**2 + 2*S**(-1)*M2 + 2
     +    *S**(-1)*MS2 - 2*T1**(-1)*U1**(-1)*M2*MS2 - 2*T1**(-1)*
     +    U1**(-1)*M2**2 + 2*T1**(-1)*M2 + 2*T1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(31)*N*CF * (  - 2*S**(-1)*U1**(-1)*M2 - 2*
     +    S**(-1)*U1**(-1)*MS2 - 2*T1**(-1)*U1**(-1)*M2 - 2*T1**(-1)*
     +    U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(31)*CF**2*S4G2**(-1)*S4G * ( 4*T1**(-1)*
     +    U1**(-1)*M2*MS2 + 4*T1**(-1)*U1**(-1)*M2**2 - 4*T1**(-1)*M2
     +     - 4*T1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(31)*CF**2 * ( 4*T1**(-1)*U1**(-1)*M2 + 4*
     +    T1**(-1)*U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(36)*CF**2 * (  - 8*S*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(38)*N*CF*S4G2**(-1)*S4G * (  - 6*S*TG + 4*
     +    S*U1**(-1)*M2**2 - 6*S*M2 - 4*S*MS2 - 8*S**2*T1**(-1)*M2 - 8*
     +    S**2*T1**(-1)*MS2 - 4*S**2 - 2*TG**2 - 2*M2**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(38)*N*CF * ( 4*S*TG*(TG+UG)**(-1) + 8*S*
     +    U1**(-1)*M2 + 2*S*M2*(TG+UG)**(-1) + 4*S*MS2*(TG+UG)**(-1) - 
     +    2*S*(TG+UG)**(-1)*S4 + 2*S + 2*S**2*(TG+UG)**(-1) + 2*TG*M2*
     +    (TG+UG)**(-1) - 2*TG*(TG+UG)**(-1)*S4 + 2*TG + 2*TG**2*
     +    (TG+UG)**(-1) + 2*M2*(TG+UG)**(-1)*S4 - 2*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(38)*CF**2*S4G2**(-1)*S4G * ( 4*S*TG - 8*S*
     +    U1**(-1)*M2**2 + 4*S*M2 + 8*S**2*T1**(-1)*M2 + 8*S**2*
     +    T1**(-1)*MS2 + 4*S**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(38)*CF**2 * (  - 8*S*TG*(TG+UG)**(-1) - 8*
     +    S*U1**(-1)*M2 - 4*S*M2*(TG+UG)**(-1) - 8*S*MS2*(TG+UG)**(-1)
     +     + 4*S*(TG+UG)**(-1)*S4 - 4*S - 4*S**2*(TG+UG)**(-1) - 4*TG*
     +    M2*(TG+UG)**(-1) + 4*TG*(TG+UG)**(-1)*S4 - 4*TG - 4*TG**2*
     +    (TG+UG)**(-1) - 8*U1**(-1)*M2*S4 - 4*M2*(TG+UG)**(-1)*S4 + 12
     +    *M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(41)*N*CF*S4G2**(-1)*S4G * (  - 2 + 2*S*
     +    U1**(-1) - 8*T1**(-1)*M2 - 8*T1**(-1)*MS2 - 4*U1**(-1)*M2 - 4
     +    *U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(41)*N*CF * (  - 8*S**(-1)*T1**(-1)*M2 - 8*
     +    S**(-1)*T1**(-1)*MS2 - 2*S**(-1)*U1**(-1)*M2 - 4*S**(-1)*
     +    U1**(-1)*MS2 - 2*S**(-1)*U1**(-1)*S4 - 4*S**(-1)*M2*
     +    (TG+UG)**(-1) - 4*S**(-1)*(TG+UG)**(-1)*S4 + 2*S**(-1) + 4*
     +    (TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(41)*CF**2*S4G2**(-1)*S4G * (  - 4*S*
     +    U1**(-1) + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 8*U1**(-1)*M2 + 8
     +    *U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(41)*CF**2 * ( 8*S**(-1)*T1**(-1)*M2 + 8*
     +    S**(-1)*T1**(-1)*MS2 + 4*S**(-1)*U1**(-1)*M2 + 8*S**(-1)*
     +    U1**(-1)*MS2 + 4*S**(-1)*U1**(-1)*S4 + 8*S**(-1)*M2*
     +    (TG+UG)**(-1) + 8*S**(-1)*(TG+UG)**(-1)*S4 - 8*S**(-1) - 8*
     +    (TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(42)*N*CF * (  - 4*S**(-1)*TG*M2 + 2*
     +    S**(-1)*TG*M2**2*(TG+UG)**(-1) - 4*S**(-1)*TG*MS2 + 2*S**(-1)
     +    *TG*(TG+UG)**(-1)*S4**2 - 16*S**(-1)*T1**(-1)*M2*MS2*S4 - 8*
     +    S**(-1)*T1**(-1)*M2*S4**2 - 8*S**(-1)*T1**(-1)*M2**2*MS2 - 16
     +    *S**(-1)*T1**(-1)*M2**2*S4 - 8*S**(-1)*T1**(-1)*M2**3 - 8*
     +    S**(-1)*T1**(-1)*MS2*S4**2 + 8*S**(-1)*M2*MS2 + 12*S**(-1)*M2
     +    *S4 - 2*S**(-1)*M2**2*(TG+UG)**(-1)*S4 + 6*S**(-1)*M2**2 + 12
     +    *S**(-1)*MS2*S4 - 2*S**(-1)*(TG+UG)**(-1)*S4**3 - 2*S**(-1)*
     +    S4**2 - 4*S*U1**(-1)*M2 - 4*TG*MS2*(TG+UG)**(-1) + 4*U1**(-1)
     +    *M2*S4 + 4*M2*(TG+UG)**(-1)*S4 - 4*M2 + 2*M2**2*(TG+UG)**(-1)
     +     + 4*MS2*(TG+UG)**(-1)*S4 + 2*(TG+UG)**(-1)*S4**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(42)*CF**2 * ( 8*S**(-1)*TG*M2 - 4*S**(-1)*
     +    TG*M2**2*(TG+UG)**(-1) + 8*S**(-1)*TG*MS2 - 4*S**(-1)*TG*
     +    (TG+UG)**(-1)*S4**2 + 16*S**(-1)*T1**(-1)*M2*MS2*S4 + 8*
     +    S**(-1)*T1**(-1)*M2*S4**2 + 8*S**(-1)*T1**(-1)*M2**2*MS2 + 16
     +    *S**(-1)*T1**(-1)*M2**2*S4 + 8*S**(-1)*T1**(-1)*M2**3 + 8*
     +    S**(-1)*T1**(-1)*MS2*S4**2 - 8*S**(-1)*M2*MS2 - 16*S**(-1)*M2
     +    *S4 + 4*S**(-1)*M2**2*(TG+UG)**(-1)*S4 - 8*S**(-1)*M2**2 - 16
     +    *S**(-1)*MS2*S4 + 4*S**(-1)*(TG+UG)**(-1)*S4**3 + 8*TG*MS2*
     +    (TG+UG)**(-1) - 8*M2*(TG+UG)**(-1)*S4 - 4*M2**2*(TG+UG)**(-1)
     +     - 8*MS2*(TG+UG)**(-1)*S4 - 4*(TG+UG)**(-1)*S4**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(43)*CF**2 * (  - 8*S**(-1)*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(49)*N*CF * ( 2*S**(-1)*TG*(TG+UG)**(-1) - 
     +    8*S**(-1)*U1**(-1)*MS2 - 4*S**(-1)*U1**(-1)*S4 - 2*S**(-1)*M2
     +    *(TG+UG)**(-1) - 2*S**(-1)*(TG+UG)**(-1)*S4 + 4*S**(-1) + 2*
     +    (TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(49)*CF**2 * (  - 4*S**(-1)*TG*
     +    (TG+UG)**(-1) + 8*S**(-1)*U1**(-1)*MS2 + 4*S**(-1)*U1**(-1)*
     +    S4 + 4*S**(-1)*M2*(TG+UG)**(-1) + 4*S**(-1)*(TG+UG)**(-1)*S4
     +     - 4*S**(-1) - 4*(TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(51)*N*CF*S4G2**(-1)*S4G * (  - 4*S**(-1)*
     +    TG*M2 + 2*S*T1**(-1)*M2 + 2*S*T1**(-1)*MS2 - 4*S*U1**(-1)*M2
     +     - 12*S*U1**(-1)*MS2 + 4*S + 4*S**2*U1**(-1) - 2*TG - 8*
     +    T1**(-1)*M2*MS2 - 8*T1**(-1)*M2**2 - 8*U1**(-1)*M2**2 + 2*M2
     +     - 8*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(51)*N*CF * ( 4 + 2*S**(-1)*TG*M2*
     +    (TG+UG)**(-1) + 4*S**(-1)*TG*(TG+UG)**(-1)*S4 - 2*S**(-1)*TG
     +     - 6*S**(-1)*T1**(-1)*M2*MS2 - 2*S**(-1)*T1**(-1)*M2*S4 - 6*
     +    S**(-1)*T1**(-1)*M2**2 - 2*S**(-1)*T1**(-1)*MS2*S4 - 4*
     +    S**(-1)*U1**(-1)*M2*MS2 - 8*S**(-1)*U1**(-1)*M2*S4 - 4*
     +    S**(-1)*U1**(-1)*MS2*S4 - 4*S**(-1)*U1**(-1)*S4**2 + 4*
     +    S**(-1)*M2 - 4*S**(-1)*M2**2*(TG+UG)**(-1) - 4*S**(-1)*
     +    (TG+UG)**(-1)*S4**2 + 8*S**(-1)*S4 - 4*S*U1**(-1) + 2*S*
     +    (TG+UG)**(-1) + 2*TG*(TG+UG)**(-1) - 4*U1**(-1)*M2 - 8*
     +    U1**(-1)*MS2 + 6*M2*(TG+UG)**(-1) + 8*MS2*(TG+UG)**(-1) + 2*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QBH3 = M2QBH3 + ANG2(51)*CF**2*S4G2**(-1)*S4G * (  - 4*S*
     +    T1**(-1)*M2 - 4*S*T1**(-1)*MS2 + 4*S*U1**(-1)*M2 + 12*S*
     +    U1**(-1)*MS2 - 4*S - 4*S**2*U1**(-1) + 16*T1**(-1)*M2*MS2 + 
     +    16*T1**(-1)*M2**2 + 8*U1**(-1)*M2**2 - 8*M2 )
     +
      M2QBH3 = M2QBH3 + ANG2(51)*CF**2 * (  - 4 - 4*S**(-1)*TG*M2*
     +    (TG+UG)**(-1) - 8*S**(-1)*TG*(TG+UG)**(-1)*S4 + 4*S**(-1)*TG
     +     + 12*S**(-1)*T1**(-1)*M2*MS2 + 4*S**(-1)*T1**(-1)*M2*S4 + 12
     +    *S**(-1)*T1**(-1)*M2**2 + 4*S**(-1)*T1**(-1)*MS2*S4 + 4*
     +    S**(-1)*U1**(-1)*M2*MS2 + 8*S**(-1)*U1**(-1)*M2*S4 + 4*
     +    S**(-1)*U1**(-1)*MS2*S4 + 4*S**(-1)*U1**(-1)*S4**2 - 8*
     +    S**(-1)*M2 + 8*S**(-1)*M2**2*(TG+UG)**(-1) - 8*S**(-1)*MS2 + 
     +    8*S**(-1)*(TG+UG)**(-1)*S4**2 - 12*S**(-1)*S4 + 4*S*U1**(-1)
     +     - 4*S*(TG+UG)**(-1) - 4*TG*(TG+UG)**(-1) + 4*U1**(-1)*M2 + 8
     +    *U1**(-1)*MS2 - 12*M2*(TG+UG)**(-1) - 16*MS2*(TG+UG)**(-1) - 
     +    4*(TG+UG)**(-1)*S4 )
     +
      M2QBH3 = M2QBH3 + ANG2(58)*N*CF * ( 2*S**(-1)*TG*M2*MS2*SYMBY*S4
     +     + 2*S**(-1)*TG*M2**2*SYMBY*S4 - 2*S**(-1)*TG*MS2 + 4*S**(-1)
     +    *TG*S4 - 2*S**(-1)*TG**2*M2*MS2*SYMBY - 2*S**(-1)*TG**2*M2*
     +    SYMBY*S4 + 2*S**(-1)*TG**2*M2**2*SYMBY - 2*S**(-1)*TG**2*MS2*
     +    SYMBY*S4 + 2*S**(-1)*TG**2 - 2*S**(-1)*TG**3*M2*SYMBY + 2*
     +    S**(-1)*TG**3*MS2*SYMBY + 4*S**(-1)*U1**(-1)*M2*MS2*S4 - 4*
     +    S**(-1)*U1**(-1)*M2*S4**2 - 4*S**(-1)*U1**(-1)*MS2*S4**2 - 4*
     +    S**(-1)*U1**(-1)*S4**3 - 4*S**(-1)*M2*MS2 + 4*S**(-1)*M2*S4
     +     + 4*S**(-1)*MS2*S4 + 4*S**(-1)*S4**2 + 4*S*TG*(TG+UG)**(-1)
     +     + 2*S*M2*(TG+UG)**(-1) + 4*S*MS2*(TG+UG)**(-1) - 2*S*
     +    (TG+UG)**(-1)*S4 + 2*S + 2*S**2*(TG+UG)**(-1) + 2*TG*M2*MS2*
     +    SYMBY + 2*TG*M2*(TG+UG)**(-1) + 2*TG*M2**2*SYMBY - 2*TG*
     +    (TG+UG)**(-1)*S4 + 4*TG - 2*TG**2*M2*SYMBY - 2*TG**2*MS2*
     +    SYMBY + 2*TG**2*(TG+UG)**(-1) - 4*U1**(-1)*MS2*S4 + 4*
     +    U1**(-1)*S4**2 + 2*M2*(TG+UG)**(-1)*S4 + 6*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(58)*CF**2 * (  - 4*S**(-1)*TG*S4 - 4*
     +    S**(-1)*TG**2 - 4*S**(-1)*U1**(-1)*M2*MS2*S4 + 4*S**(-1)*
     +    U1**(-1)*M2*S4**2 + 4*S**(-1)*U1**(-1)*MS2*S4**2 + 4*S**(-1)*
     +    U1**(-1)*S4**3 + 4*S**(-1)*M2*MS2 - 4*S**(-1)*M2*S4 - 4*
     +    S**(-1)*MS2*S4 - 4*S**(-1)*S4**2 - 8*S*TG*(TG+UG)**(-1) - 4*S
     +    *M2*(TG+UG)**(-1) - 8*S*MS2*(TG+UG)**(-1) + 4*S*(TG+UG)**(-1)
     +    *S4 - 4*S - 4*S**2*(TG+UG)**(-1) - 4*TG*M2*(TG+UG)**(-1) + 4*
     +    TG*(TG+UG)**(-1)*S4 - 8*TG - 4*TG**2*(TG+UG)**(-1) + 4*
     +    U1**(-1)*MS2*S4 - 4*U1**(-1)*S4**2 - 4*M2*(TG+UG)**(-1)*S4 - 
     +    8*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(69)*N*CF*S4G2**(-1)*S4G * (  - 2 - 4*S*
     +    U1**(-1) + 4*T1**(-1)*M2 + 4*T1**(-1)*MS2 + 8*U1**(-1)*M2 + 8
     +    *U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(69)*N*CF * ( 2*S**(-1)*TG*(TG+UG)**(-1) + 
     +    4*S**(-1)*T1**(-1)*M2 + 4*S**(-1)*T1**(-1)*MS2 + 4*S**(-1)*
     +    U1**(-1)*M2 + 2*S**(-1)*M2*(TG+UG)**(-1) + 2*S**(-1)*
     +    (TG+UG)**(-1)*S4 - 2*S**(-1) - 2*(TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(69)*CF**2*S4G2**(-1)*S4G * ( 4*S*U1**(-1)
     +     - 8*T1**(-1)*M2 - 8*T1**(-1)*MS2 - 8*U1**(-1)*M2 - 8*
     +    U1**(-1)*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(69)*CF**2 * (  - 4*S**(-1)*TG*
     +    (TG+UG)**(-1) - 8*S**(-1)*T1**(-1)*M2 - 8*S**(-1)*T1**(-1)*
     +    MS2 - 4*S**(-1)*U1**(-1)*M2 - 4*S**(-1)*M2*(TG+UG)**(-1) - 4*
     +    S**(-1)*(TG+UG)**(-1)*S4 + 4*S**(-1) + 4*(TG+UG)**(-1) )
     +
      M2QBH3 = M2QBH3 + ANG2(73)*N*CF*S4G2**(-1)*S4G * (  - 4*S**(-1)*
     +    TG*M2**2 - 2*S*TG - 2*S*T1**(-1)*M2*MS2 - 2*S*T1**(-1)*M2**2
     +     - 12*S*U1**(-1)*M2*MS2 + 12*S*U1**(-1)*M2*S4 - 24*S*U1**(-1)
     +    *M2**2 + 6*S*M2 - 10*S*MS2 + 4*S**2*U1**(-1)*M2 - 12*S**2*
     +    U1**(-1)*MS2 - 4*S**2*U1**(-1)*S4 + 4*S**2 + 4*S**3*U1**(-1)
     +     - 4*TG*M2 + 4*TG*MS2 - 8*T1**(-1)*M2**2*MS2 - 8*T1**(-1)*
     +    M2**3 - 8*U1**(-1)*M2**2*S4 - 12*M2*MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(73)*N*CF * (  - 2*S**(-1)*TG*M2*MS2*SYMBY*
     +    S4 + 4*S**(-1)*TG*M2*(TG+UG)**(-1)*S4 - 6*S**(-1)*TG*M2 - 2*
     +    S**(-1)*TG*M2**2*SYMBY*S4 + 2*S**(-1)*TG*MS2 + 2*S**(-1)*
     +    TG**2*M2*MS2*SYMBY + 2*S**(-1)*TG**2*M2*SYMBY*S4 - 2*S**(-1)*
     +    TG**2*M2**2*SYMBY + 2*S**(-1)*TG**2*MS2*SYMBY*S4 + 2*S**(-1)*
     +    TG**3*M2*SYMBY - 2*S**(-1)*TG**3*MS2*SYMBY - 2*S**(-1)*
     +    T1**(-1)*M2*MS2*S4 - 6*S**(-1)*T1**(-1)*M2**2*MS2 - 2*S**(-1)
     +    *T1**(-1)*M2**2*S4 - 6*S**(-1)*T1**(-1)*M2**3 + 2*S**(-1)*M2*
     +    MS2 - 2*S**(-1)*M2*(TG+UG)**(-1)*S4**2 - 2*S**(-1)*M2**3*
     +    (TG+UG)**(-1) - 2*S**(-1)*MS2*S4 - 4*S*TG*(TG+UG)**(-1) - 20*
     +    S*U1**(-1)*M2 + 4*S*U1**(-1)*MS2 - 2*S*M2*(TG+UG)**(-1) + 2*S
     +    *(TG+UG)**(-1)*S4 - 2*TG*M2*MS2*SYMBY - 6*TG*M2*(TG+UG)**(-1)
     +     - 2*TG*M2**2*SYMBY - 4*TG*MS2*(TG+UG)**(-1) + 6*TG*
     +    (TG+UG)**(-1)*S4 - 6*TG + 2*TG**2*M2*SYMBY + 2*TG**2*MS2*
     +    SYMBY - 4*T1**(-1)*M2*MS2 - 4*T1**(-1)*M2**2 + 4*U1**(-1)*M2*
     +    MS2 )
     +
      M2QBH3 = M2QBH3 + ANG2(73)*N*CF * ( 4*U1**(-1)*M2**2 + 8*M2*
     +    (TG+UG)**(-1)*S4 - 6*M2 - 2*M2**2*(TG+UG)**(-1) + 4*MS2*
     +    (TG+UG)**(-1)*S4 - 6*MS2 - 2*(TG+UG)**(-1)*S4**2 + 2*S4 )
     +
      M2QBH3 = M2QBH3 + ANG2(73)*CF**2*S4G2**(-1)*S4G * ( 4*S*T1**(-1)*
     +    M2*MS2 + 4*S*T1**(-1)*M2**2 + 12*S*U1**(-1)*M2*MS2 - 12*S*
     +    U1**(-1)*M2*S4 + 24*S*U1**(-1)*M2**2 - 8*S*M2 + 8*S*MS2 - 4*
     +    S**2*U1**(-1)*M2 + 12*S**2*U1**(-1)*MS2 + 4*S**2*U1**(-1)*S4
     +     - 4*S**2 - 4*S**3*U1**(-1) + 16*T1**(-1)*M2**2*MS2 + 16*
     +    T1**(-1)*M2**3 + 8*U1**(-1)*M2**2*S4 - 8*M2**2 )
     +
      M2QBH3 = M2QBH3 + ANG2(73)*CF**2 * (  - 8*S**(-1)*TG*M2*
     +    (TG+UG)**(-1)*S4 + 8*S**(-1)*TG*M2 + 4*S**(-1)*T1**(-1)*M2*
     +    MS2*S4 + 12*S**(-1)*T1**(-1)*M2**2*MS2 + 4*S**(-1)*T1**(-1)*
     +    M2**2*S4 + 12*S**(-1)*T1**(-1)*M2**3 - 8*S**(-1)*M2*MS2 + 4*
     +    S**(-1)*M2*(TG+UG)**(-1)*S4**2 - 4*S**(-1)*M2*S4 + 4*S**(-1)*
     +    M2**2 + 4*S**(-1)*M2**3*(TG+UG)**(-1) + 8*S*TG*(TG+UG)**(-1)
     +     + 20*S*U1**(-1)*M2 - 4*S*U1**(-1)*MS2 + 4*S*M2*(TG+UG)**(-1)
     +     - 4*S*(TG+UG)**(-1)*S4 + 12*TG*M2*(TG+UG)**(-1) + 8*TG*MS2*
     +    (TG+UG)**(-1) - 12*TG*(TG+UG)**(-1)*S4 + 12*TG + 8*T1**(-1)*
     +    M2*MS2 + 8*T1**(-1)*M2**2 - 4*U1**(-1)*M2*MS2 - 4*U1**(-1)*
     +    M2**2 - 16*M2*(TG+UG)**(-1)*S4 + 8*M2 + 4*M2**2*(TG+UG)**(-1)
     +     - 8*MS2*(TG+UG)**(-1)*S4 + 4*MS2 + 4*(TG+UG)**(-1)*S4**2 - 4
     +    *S4 )

      END IF

      Q2MS2 = (MS + MG)**2/4.D0/MS**2

      IF (IFL.EQ.1) M2QBH = 
     +     2.D0*((NS -1.D0)*M2QBH2 + M2QBH1 + M2QBH3)
      IF (IFL.EQ.0) M2QBH = 2.D0*M2QBH1 

      DSGQBH = S4/(S4+MS2)/2.D0 *
     +     ALPHAS**3 * AVG * M2QBH /4.D0 /S**2 *CONV
     +     + DSGQB3(ALPHAS,S,TG,S4,MS,MG,Q2MS2,IFL)

      RETURN
      END


      REAL*8 FUNCTION DSGQBS(ALPHAS,S,TG,S4,MS,MG,IFL)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR Q +QB -> SQ + GL +QB
C***  THE 1/S4G**2 PART
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,UG,MS,MG,MS2,MG2,M2,CF
      REAL*8 M2QBS, M2QBS2, NS,CONV,N,AVG
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      
      U1 = S4 -S -TG
      T1 = TG +M2
      UG = U1 -M2
      
      M2QBS2 = 0.D0
      M2QBS2 = M2QBS2 + N*CF**2 * ( 8*T1**(-2)*M2**3 - 32*T1**(-1)*
     +    U1**(-1)*M2**2*MS2 - 32*T1**(-1)*U1**(-1)*M2**3 + 16*T1**(-1)
     +    *M2*MS2 + 8*U1**(-2)*M2**3 + 16*U1**(-1)*M2*MS2 + 16*M2 )
     +
      M2QBS2 = M2QBS2 + N**2*CF * ( 16*S**(-2)*TG**2*M2 + 16*S**(-1)*TG
     +    *M2 + 8*S**(-1)*T1**(-1)*M2**3 + 8*S**(-1)*U1**(-1)*M2**3 + 
     +    16*S**(-1)*M2*MS2 + 16*T1**(-1)*U1**(-1)*M2**2*MS2 + 16*
     +    T1**(-1)*U1**(-1)*M2**3 )

      IF (IFL.EQ.1) M2QBS = 2.D0*(NS -1.D0)*M2QBS2 
      IF (IFL.EQ.0) M2QBS = 0.D0

      DSGQBS = ALPHAS**3 *AVG *M2QBS *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END



      REAL*8 FUNCTION DSGQBT(ALPHAS,S,TG,S4,S3,MS,MG,IFL)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR Q +QB -> SQ + GL +QB
C***  THE 1/S3**2 PART

      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,S4G,S3,T1,TG,U1,UG,MS,MG,MS2,MG2,M2
      REAL*8 NS,CONV,N,AVG,CF, M2QBT, M2QBT1, M2QBT2, M2QBT3
      REAL*8 ANGDEF(1:11), ANA(2:2,1:9), ANB(2:2,1:9), ANC(2:2,1:9)
      REAL*8 ANGS3(1:8),XX(5:9),YY2(5:9),XPHI
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2

      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)

      XPHI = (S3 -ANA(2,1))/ANB(2,1)

      XX(5) = (ANA(2,5) + ANB(2,5)*XPHI)
      YY2(5)= ANC(2,5)**2 * (1.D0-XPHI**2)
      XX(6) = (ANA(2,6) + ANB(2,6)*XPHI)
      YY2(6)= ANC(2,6)**2 * (1.D0-XPHI**2)
      XX(8) = (ANA(2,8) + ANB(2,8)*XPHI)
      YY2(8)= ANC(2,8)**2 * (1.D0-XPHI**2)
      XX(9) = (ANA(2,9) + ANB(2,9)*XPHI)
      YY2(9)= ANC(2,9)**2 * (1.D0-XPHI**2)

      ANGS3(1) = -XX(6)/(XX(6)**2 - YY2(6))**(1.5D0)
      ANGS3(2) = -1.D0/SQRT(XX(6)**2 - YY2(6))
      ANGS3(3) = -XX(5)/(XX(5)**2 - YY2(5))**(1.5D0)
      ANGS3(4) = -1.D0/SQRT(XX(5)**2 - YY2(5))
      ANGS3(5) = XX(5)
      ANGS3(6) = XX(5)**2 +0.5D0*YY2(5)
      ANGS3(7) = -XX(9)/(XX(9)**2 - YY2(9))**(1.5D0)
      ANGS3(8) = -1.D0/SQRT(XX(9)**2 - YY2(9))

      M2QBT1 = 0.D0
      M2QBT1 = M2QBT1 + N*CF**2 * ( 4*M2 )
     +
      M2QBT1 = M2QBT1 + ANGS3(7)*N*CF**2 * ( 4*M2**3 )
     +
      M2QBT1 = M2QBT1 + ANGS3(8)*N*CF**2 * ( 4*S*M2 + 8*M2**2 )


      M2QBT2 = 0.D0
      M2QBT2 = M2QBT2 + N*CF**2 * ( 8*S**(-1)*M2*MS2 )
     +
      M2QBT2 = M2QBT2 + ANGS3(5)*N*CF**2 * ( 8*S**(-1)*M2 )
     +
      M2QBT2 = M2QBT2 + ANGS3(6)*N*CF**2 * ( 8*S**(-2)*M2 )

      M2QBT3 = 0.D0
      M2QBT3 = M2QBT3 + CF**2 * (  - 8*S**(-1)*M2**2 )
     +
      M2QBT3 = M2QBT3 + ANGS3(5)*CF**2 * ( 8*S**(-1)*M2 )
     +
      M2QBT3 = M2QBT3 + ANGS3(8)*CF**2 * (  - 8*S**(-1)*M2**3 - 8*M2*
     +    MS2 - 8*M2**2 )


      IF (IFL.EQ.1) M2QBT = 
     +     2.D0*((NS -1.D0)*M2QBT2 + M2QBT1 + M2QBT3)
      IF (IFL.EQ.0) M2QBT = 2.D0*M2QBT1 

      DSGQBT = ALPHAS**3 *AVG *M2QBT *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSGQQ3(ALPHAS,S,TG,S4,MS,MG,SCA,IFL)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,TG,U1,UG,S4,MS,MG,SCA,X1,X2
      REAL*8 PI,N,CF,NS, DSGQGB, PGQ1, PGQ2
      INTEGER IFL

      N = 3.D0
      CF = (N**2 -1.D0)/2.D0/N
      NS = 6.D0
      PI = 4.D0*ATAN(1.D0)
      U1 = S4 -S -TG
      T1 = TG +MG**2 -MS**2
      UG = U1 -MG**2 +MS**2
      X1 = -T1/(S+UG)
      X2 = -U1/(S+TG)
      PGQ1 = CF*(1.D0 + (1.D0 -X1)**2)/X1
      PGQ2 = CF*(1.D0 + (1.D0 -X2)**2)/X2

      DSGQQ3 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) * 
     +     (-1.D0/T1*DSGQGB(ALPHAS,X1*S,X1*UG,MS,MG)*PGQ1*X1**2 
     +      -1.D0/U1*DSGQGB(ALPHAS,X2*S,X2*TG,MS,MG)*PGQ2*X2**2 )

      RETURN
      END


      REAL*8 FUNCTION DSGQQH(ALPHAS,S,TG,S4,MS,MG,IFL,EPSS,EPSG)
C***  CROSS SECTIONS FOR Q + Q -> SQ + GL +Q
C***  M2QQH1: Q QP --> GL SQ QP
C***  M2QQH2: Q QP --> GL SQP Q
C***  M2QQH3: INTERFERENCE 

      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,UG,MS,MG,MS2,MG2,M2,EPSS,EPSG
      REAL*8 NS,CONV,N,CF,AVG,S4G, S4G2,DEL,SYMBT,SYMBU,SYMBX,SYMBY
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 AHP1P1, ABP1P1, ABP1P2, A4P1P2, A4P1P1
      REAL*8 A4P0P2, A4P1P0, ABP2P2
      REAL*8 ABP2P1, ABP2P0
      REAL*8 C4P1P2, C4P1P1, C4P1P0, CBP1P1
      REAL*8 M2QQH,M2QQH1,M2QQH2,M2QQH3
      REAL*8 ANG2(1:87), COLO2(1:9),DSGQQ3,Q2MS2
      INTEGER IFL
           
      CONV = 389379660.D0
      NS = 6.D0

      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0

      AVG = (1.D0/2.D0)**2 /N**2
      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2
      S4G2 = S4G**2 + EPSG*MG**4
      SYMBT = (M2*S-TG*S4G)/((M2*S-TG*S4G)**2 +(S+TG)*EPSS*MS**4 )
      SYMBU = (M2*S-UG*S4G)/((M2*S-UG*S4G)**2 +(S+UG)*EPSS*MS**4 )
      SYMBX = 1.D0/((UG-M2)**2 +EPSS*MS**4)
      SYMBY = 1.D0/((TG-M2)**2 +EPSS*MS**4)


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6) + M2
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9) 
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) + MS2 +MG2
      ANB(1,2) = -ANB(1,1)
      ANC(1,2) = -ANC(1,1)
      ANA(1,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(1,3) = -ANA(1,3)
      ANC(1,3) =  0.D0
      ANA(1,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(1,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)-2.D0*ANGDEF(3)*ANGDEF(4)
      ANC(1,4) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(1,5) = -ANB(1,3)
      ANC(1,5) = -ANC(1,3)
      ANA(1,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(1,6) = -ANB(1,4)
      ANC(1,6) = -ANC(1,4)
      ANA(1,7) = +ANA(1,1) -M2
      ANB(1,7) = +ANB(1,1)
      ANC(1,7) = +ANC(1,1)
      ANA(1,8) = +ANA(1,5) -M2
      ANB(1,8) = +ANB(1,5)
      ANC(1,8) = +ANC(1,5)
      ANA(1,9) = +ANA(1,6) -M2
      ANB(1,9) = +ANB(1,6)
      ANC(1,9) = +ANC(1,6)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(3,2) = -ANB(3,1)
      ANC(3,2) = -ANC(3,1)
      ANA(3,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(3,3) = 
     +     2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10) -2.D0*ANGDEF(2)*ANGDEF(4)
      ANC(3,3) = +2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(3,4) = -ANA(3,4)
      ANC(3,4) =  0.D0
      ANA(3,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(3,5) = -ANB(3,3)
      ANC(3,5) = -ANC(3,3)
      ANA(3,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(3,6) = -ANB(3,4)
      ANC(3,6) = -ANC(3,4)
      ANA(3,7) = +ANA(3,1) -M2
      ANB(3,7) = +ANB(3,1)
      ANC(3,7) = +ANC(3,1)
      ANA(3,8) = +ANA(3,5) -M2
      ANB(3,8) = +ANB(3,5)
      ANC(3,8) = +ANC(3,5)
      ANA(3,9) = +ANA(3,6) -M2
      ANB(3,9) = +ANB(3,6)
      ANC(3,9) = +ANC(3,6)

C$$$      ANG2(1) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(2) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(3) = A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(4) = A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(5) = A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(6) = A4M1P1(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(7) = A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(8) = A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(9) = A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(10)= A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))

C$$$      ANG2(11)= A4P0P2(ANA(1,7),ANB(1,7),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(12)= A4P1P0(ANA(2,7),ANB(2,7),ANA(2,7),ANB(2,7),ANC(2,7))
      ANG2(13)= A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(14)= A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(15)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(16)= A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(17)= ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(18)= A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG2(19)= A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(20)= AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))

C$$$      ANG2(21)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG2(22)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG2(23)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(24)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(25)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(26)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(27)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(28)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(29)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(30)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))

C$$$      ANG2(31)= ABP1M1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(32)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(33)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(34)= ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(35)= ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(36)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(37)= ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(38)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(39)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(40)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))

C$$$      ANG2(41)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(42)= A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(43)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(44)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG2(45)= A4P2M2(ANA(2,2),ANB(2,2),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG2(46)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG2(47)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(48)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(49)= A4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG2(50)= A4M1P0(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))

      ANG2(51)= A4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(52)= A4P1P1(ANA(2,2),ANB(2,2),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(53)= A4P0P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG2(54)= A4P1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(55)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(56)= A4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(57)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(59)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG2(60)= ABP1M2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))

C$$$      ANG2(61)= ABP1M2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG2(62)= ABP1P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(63)= ABP1P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(64)= A4M1P2(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(65)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG2(66)= A4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG2(67)= A4M1P2(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(68)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(69)= A4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG2(70)= A4M1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG2(71)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG2(72)= A4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG2(73)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(74)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
      ANG2(75)= A4P1P1(ANA(1,8),ANB(1,8),ANA(1,2),ANB(1,2),ANC(1,2))
      ANG2(76)= A4P1P1(ANA(1,8),ANB(1,8),ANA(1,9),ANB(1,9),ANC(1,9))
      ANG2(77)= A4P1P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
      ANG2(78)= A4P1P1(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
      ANG2(79)= A4P0P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
      ANG2(80)= A4P1P0(ANA(1,8),ANB(1,8),ANA(1,2),ANB(1,2),ANC(1,2))

      ANG2(81)= ABP2P0(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG2(82)= ABP2P2(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8))
      ANG2(83)= ABP1P2(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8))
      ANG2(84)= ABP2P1(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8))
      ANG2(85)= ABP1P1(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8))
      ANG2(86)= ABP2P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
      ANG2(87)= ABP2P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))

      DEL = EPSS * MS**4
      IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
C$$$  ANG2(47)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$  ANG2(48)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$  ANG2(49)= C4M1P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
         ANG2(51)= C4P1P0(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$  ANG2(56)= C4P1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$  ANG2(57)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
      ANG2(58)= CBP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1),DEL)
      ANG2(59)= CBP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1),DEL)
C$$$  ANG2(66)= C4M2P1(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$  ANG2(69)= C4M1P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$  ANG2(72)= C4M2P1(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
         ANG2(73)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
         ANG2(74)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
         ANG2(77)= C4P1P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
         ANG2(78)= C4P1P1(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
      END IF

      COLO2(9) = LOG(S4**2/MS2/(S4+MS2))

      M2QQH1 = 0D0
      M2QQH1 = M2QQH1 + N*CF**2*(S4+MS2) * (  - 32*S**(-1)*U1**(-4)*
     +    M2**2*S4 + 32*S**(-1)*U1**(-3)*M2*S4 + 32*S**(-1)*U1**(-3)*
     +    M2**2 + 16*S**(-1)*U1**(-2)*M2*MG2*(S+TG)**(-1) - 16*S**(-1)*
     +    U1**(-2)*M2*MG2*S4**(-1) - 32*S**(-1)*U1**(-2)*M2 - 16*
     +    S**(-1)*U1**(-2)*M2**2*(S+TG)**(-1) - 16*S**(-1)*U1**(-1)*M2*
     +    MG2*(S+TG)**(-1)*S4**(-1) + 16*S**(-1)*U1**(-1)*M2*S4**(-1)
     +     + 16*S**(-1)*U1**(-1)*M2**2*(S+TG)**(-1)*S4**(-1) - 8*
     +    S**(-1)*S4**(-1) - 16*S*TG**(-1)*M2*(TG+UG)**(-1)*SYMBY + 16*
     +    S*TG**(-1)*M2*SYMBY*S4**(-1) + 16*S*TG**(-1)*(TG+UG)**(-1)*
     +    S4**(-1) + 16*S*TG*(TG+UG)**(-1)*SYMBT*SYMBY*S4 + 16*S*TG*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) - 16*S*TG*SYMBT*SYMBY + 16*S*
     +    TG**2*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*S*TG**3*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 16*S*T1**(-1)*
     +    (S+TG)**(-1)*S4**(-1) - 16*S*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4
     +     + 16*S*M2*SYMBT*SYMBY + 16*S*(TG+UG)**(-1)*SYMBY + 16*S*
     +    SYMBT*S4**(-1) )
     +
      M2QQH1 = M2QQH1 + N*CF**2*(S4+MS2) * (  - 16*S*SYMBY*S4**(-1) + 
     +    16*S**2*TG**(-1)*M2*(TG+UG)**(-1)*SYMBY*S4**(-1) - 16*S**2*TG
     +    *M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 32*S**2*TG*
     +    (TG+UG)**(-1)*SYMBT*SYMBY + 16*S**2*TG*SYMBT*SYMBY*S4**(-1)
     +     + 16*S**2*TG**2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 32*S**2
     +    *M2*(TG+UG)**(-1)*SYMBT*SYMBY - 16*S**2*M2*SYMBT*SYMBY*
     +    S4**(-1) - 16*S**2*(TG+UG)**(-1)*SYMBY*S4**(-1) + 16*S**3*TG*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*S**3*M2*(TG+UG)**(-1)
     +    *SYMBT*SYMBY*S4**(-1) + 16*TG*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4
     +     - 16*TG*M2*(TG+UG)**(-1)*SYMBY*S4**(-1) - 16*TG*M2*SYMBT*
     +    SYMBY - 16*TG*(TG+UG)**(-1)*SYMBT - 16*TG*(TG+UG)**(-1)*SYMBY
     +     + 16*TG*SYMBT*S4**(-1) + 16*TG*SYMBY*S4**(-1) - 32*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY + 16*TG**2*M2*SYMBT*SYMBY*S4**(-1)
     +     - 16*TG**2*(TG+UG)**(-1)*SYMBT*SYMBY*S4 + 16*TG**2*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) + 16*TG**2*(TG+UG)**(-1)*SYMBY*
     +    S4**(-1) )
     +
      M2QQH1 = M2QQH1 + N*CF**2*(S4+MS2) * ( 16*TG**2*SYMBT*SYMBY + 16*
     +    TG**3*M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 32*TG**3*
     +    (TG+UG)**(-1)*SYMBT*SYMBY - 16*TG**3*SYMBT*SYMBY*S4**(-1) - 
     +    16*TG**4*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*T1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) - 16*T1**(-1)*S4**(-1) - 32*U1**(-4)*M2
     +    *MG2 + 32*U1**(-4)*M2**2 - 32*U1**(-3)*M2 - 16*U1**(-2)*M2*
     +    MG2*(S+TG)**(-1)*S4**(-1) + 16*U1**(-2)*M2**2*(S+TG)**(-1)*
     +    S4**(-1) - 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) + 16*M2*
     +    (TG+UG)**(-1)*SYMBY - 16*M2*SYMBY*S4**(-1) + 24*(S+TG)**(-1)*
     +    S4**(-1) )
     +
      M2QQH1 = M2QQH1 + N**2*CF*(S4+MS2) * ( 8*S*TG**(-1)*M2*
     +    (TG+UG)**(-1)*SYMBY - 8*S*TG**(-1)*M2*SYMBY*S4**(-1) - 8*S*
     +    TG**(-1)*(TG+UG)**(-1)*S4**(-1) + 4*S*TG*M2*SYMBT*SYMBY*
     +    S4**(-1) - 8*S*TG*(TG+UG)**(-1)*SYMBT*SYMBY*S4 - 8*S*TG*
     +    (TG+UG)**(-1)*SYMBT*S4**(-1) + 8*S*TG*SYMBT*SYMBY - 8*S*TG**2
     +    *M2*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 4*S*TG**2*SYMBT*
     +    SYMBY*S4**(-1) + 8*S*TG**3*(TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1)
     +     - 8*S*T1**(-1)*(S+TG)**(-1)*S4**(-1) + 8*S*M2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4 - 8*S*M2*SYMBT*SYMBY - 8*S*(TG+UG)**(-1)*SYMBY
     +     - 4*S*SYMBT*S4**(-1) + 8*S*SYMBY*S4**(-1) - 8*S**2*TG**(-1)*
     +    M2*(TG+UG)**(-1)*SYMBY*S4**(-1) + 8*S**2*TG*M2*(TG+UG)**(-1)*
     +    SYMBT*SYMBY*S4**(-1) + 16*S**2*TG*(TG+UG)**(-1)*SYMBT*SYMBY
     +     - 8*S**2*TG*SYMBT*SYMBY*S4**(-1) - 8*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*S**2*M2*(TG+UG)**(-1)
     +    *SYMBT*SYMBY + 8*S**2*M2*SYMBT*SYMBY*S4**(-1) + 8*S**2*
     +    (TG+UG)**(-1)*SYMBY*S4**(-1) )
     +
      M2QQH1 = M2QQH1 + N**2*CF*(S4+MS2) * (  - 8*S**3*TG*(TG+UG)**(-1)
     +    *SYMBT*SYMBY*S4**(-1) + 8*S**3*M2*(TG+UG)**(-1)*SYMBT*SYMBY*
     +    S4**(-1) - 32*TG**(-2)*U1**(-4)*M2*MG2*S4**2 + 64*TG**(-2)*
     +    U1**(-3)*M2*MG2*S4 - 48*TG**(-2)*U1**(-2)*M2*MG2 + 16*
     +    TG**(-2)*U1**(-1)*M2*MG2*S4**(-1) + 32*TG**(-1)*U1**(-4)*M2*
     +    MG2*S4 - 32*TG**(-1)*U1**(-4)*M2**2*S4 - 32*TG**(-1)*U1**(-3)
     +    *M2*MG2 + 32*TG**(-1)*U1**(-3)*M2*S4 + 32*TG**(-1)*U1**(-3)*
     +    M2**2 + 16*TG**(-1)*U1**(-2)*M2*MG2*S4**(-1) - 32*TG**(-1)*
     +    U1**(-2)*M2 - 16*TG**(-1)*U1**(-2)*M2**2*S4**(-1) + 16*
     +    TG**(-1)*U1**(-1)*M2*S4**(-1) - 8*TG**(-1)*S4**(-1) - 8*TG*M2
     +    *(TG+UG)**(-1)*SYMBT*SYMBY*S4 + 8*TG*M2*(TG+UG)**(-1)*SYMBY*
     +    S4**(-1) + 4*TG*M2*SYMBT*SYMBY + 8*TG*(TG+UG)**(-1)*SYMBT + 8
     +    *TG*(TG+UG)**(-1)*SYMBY - 4*TG*SYMBT*S4**(-1) - 4*TG*SYMBY*
     +    S4**(-1) + 16*TG**2*M2*(TG+UG)**(-1)*SYMBT*SYMBY - 4*TG**2*M2
     +    *SYMBT*SYMBY*S4**(-1) + 8*TG**2*(TG+UG)**(-1)*SYMBT*SYMBY*S4
     +     - 8*TG**2*(TG+UG)**(-1)*SYMBT*S4**(-1) )
     +
      M2QQH1 = M2QQH1 + N**2*CF*(S4+MS2) * (  - 8*TG**2*(TG+UG)**(-1)*
     +    SYMBY*S4**(-1) - 4*TG**2*SYMBT*SYMBY - 8*TG**3*M2*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) - 16*TG**3*(TG+UG)**(-1)*
     +    SYMBT*SYMBY + 4*TG**3*SYMBT*SYMBY*S4**(-1) + 8*TG**4*
     +    (TG+UG)**(-1)*SYMBT*SYMBY*S4**(-1) + 8*T1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) + 8*T1**(-1)*S4**(-1) - 8*M2*
     +    (TG+UG)**(-1)*SYMBY + 4*M2*SYMBY*S4**(-1) - 8*(S+TG)**(-1)*
     +    S4**(-1) )
     +
      M2QQH1 = M2QQH1 + ANG2(13)*N*CF**2 * ( 4*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(14)*N*CF**2 * ( 4 + 4*S*U1**(-1) - 8*S*
     +    (TG+UG)**(-1) + 16*T1**(-1)*M2 + 16*T1**(-1)*MS2 - 4*U1**(-1)
     +    *M2 - 8*U1**(-1)*MS2 - 4*U1**(-1)*S4 + 8*M2*(TG+UG)**(-1) + 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(14)*N**2*CF * (  - 2*S*U1**(-1) + 4*S*
     +    (TG+UG)**(-1) - 4*T1**(-1)*M2 - 4*T1**(-1)*MS2 + 2*U1**(-1)*
     +    M2 + 4*U1**(-1)*MS2 + 2*U1**(-1)*S4 - 4*M2*(TG+UG)**(-1) - 4*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(15)*N*CF**2 * ( 4*T1**(-2)*M2**3 + 8*
     +    T1**(-1)*M2*MG2 - 8*T1**(-1)*M2**2 + 8*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(16)*N*CF**2 * ( 12 + 8*S*U1**(-1) - 4*
     +    T1**(-2)*M2*S4 + 4*T1**(-2)*M2**2 + 8*T1**(-1)*U1**(-1)*M2*S4
     +     + 16*T1**(-1)*U1**(-1)*MS2*S4 + 8*T1**(-1)*U1**(-1)*S4**2 + 
     +    4*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 4*T1**(-1)*S4 - 16*U1**(-1)*
     +    M2 - 16*U1**(-1)*MS2 - 8*U1**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(16)*N**2*CF * (  - 4 - 4*S*U1**(-1) - 4*
     +    T1**(-1)*U1**(-1)*M2*S4 - 8*T1**(-1)*U1**(-1)*MS2*S4 - 4*
     +    T1**(-1)*U1**(-1)*S4**2 + 8*U1**(-1)*M2 + 8*U1**(-1)*MS2 + 4*
     +    U1**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(17)*N*CF**2 * ( 8*U1**(-2)*M2*S4**2 - 16*
     +    U1**(-1)*M2*S4 + 8*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(27)*N**2*CF * ( 8*TG**2*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(28)*N**2*CF * ( 12*TG*M2 + 4*TG**2*M2**2*
     +    SYMBY - 4*TG**3*M2*SYMBY + 4*T1**(-1)*M2**3 + 8*M2*MG2 - 4*
     +    M2**2 )
     +
      M2QQH1 = M2QQH1 + ANG2(29)*N**2*CF * (  - 8*S*M2 - 8*S**2*
     +    U1**(-1)*M2 - 8*TG*M2 + 8*U1**(-1)*M2*S4**2 - 8*M2*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(30)*N*CF**2 * (  - 24*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 8*S*TG*M2*SYMBY - 8*S*TG*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 16*S*TG*(TG+UG)**(-1) + 8*S*TG*
     +    SYMBY*S4 + 16*S*TG**2*M2*(TG+UG)**(-1)*SYMBY + 32*S*TG**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 16*S*TG**2*SYMBY - 16*S*TG**3*
     +    (TG+UG)**(-1)*SYMBY + 16*S*U1**(-1)*M2*(TG+UG)**(-1)*S4 + 32*
     +    S*U1**(-1)*M2*SYMBY*S4**2 + 16*S*U1**(-1)*M2 - 16*S*U1**(-1)*
     +    M2**2*SYMBY*S4 + 8*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*S*M2*
     +    (TG+UG)**(-1) - 24*S*M2*SYMBY*S4 - 8*S*M2**2*(TG+UG)**(-1)*
     +    SYMBY*S4 + 8*S*M2**2*SYMBY - 8*S*(TG+UG)**(-1)*S4 + 8*S**2*TG
     +    *M2*(TG+UG)**(-1)*SYMBY + 8*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 - 
     +    16*S**2*TG**2*(TG+UG)**(-1)*SYMBY - 16*S**2*U1**(-1)*M2*
     +    (TG+UG)**(-1) - 16*S**2*U1**(-1)*M2*SYMBY*S4 - 8*S**2*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 8*S**2*M2**2*(TG+UG)**(-1)*SYMBY + 8
     +    *TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 + 8*TG*M2*(TG+UG)**(-1) + 8*
     +    TG*M2*SYMBY*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(30)*N*CF**2 * (  - 8*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 8*TG*M2**2*SYMBY - 8*TG*
     +    (TG+UG)**(-1)*S4 + 8*TG**2*M2*SYMBY + 8*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 8*
     +    TG**2*SYMBY*S4 - 8*TG**3*M2*(TG+UG)**(-1)*SYMBY + 8*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 - 16*T1**(-1)*U1**(-1)*M2**3 - 16*
     +    U1**(-1)*M2*SYMBY*S4**3 - 16*U1**(-1)*M2*S4 + 16*U1**(-1)*
     +    M2**2*SYMBY*S4**2 + 16*U1**(-1)*M2**2 + 16*M2*SYMBY*S4**2 + 
     +    16*M2 - 16*M2**2*SYMBY*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(30)*N**2*CF * ( 12*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 2*S*TG*M2*SYMBY + 4*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 8*S*TG*(TG+UG)**(-1) - 4*S*TG*SYMBY*S4 - 8*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 16*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 6*S
     +    *TG**2*SYMBY + 8*S*TG**3*(TG+UG)**(-1)*SYMBY + 8*S*U1**(-2)*
     +    M2*S4 - 8*S*U1**(-1)*M2*(TG+UG)**(-1)*S4 - 8*S*U1**(-1)*M2*
     +    SYMBY*S4**2 - 16*S*U1**(-1)*M2 + 4*S*U1**(-1)*M2**2*SYMBY*S4
     +     - 4*S*U1**(-1)*S4 - 4*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 + 4*S*
     +    M2*(TG+UG)**(-1) + 8*S*M2*SYMBY*S4 + 4*S*M2**2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 4*S*M2**2*SYMBY + 4*S*(TG+UG)**(-1)*S4 + 4*S - 4*
     +    S**2*TG*M2*(TG+UG)**(-1)*SYMBY - 4*S**2*TG*(TG+UG)**(-1)*
     +    SYMBY*S4 + 8*S**2*TG**2*(TG+UG)**(-1)*SYMBY - 8*S**2*U1**(-2)
     +    *M2 + 8*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 4*S**2*U1**(-1)*M2*
     +    SYMBY*S4 + 4*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 - 4*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 4*TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 4*
     +    TG*M2*(TG+UG)**(-1) )
     +
      M2QQH1 = M2QQH1 + ANG2(30)*N**2*CF * ( 4*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 2*TG*M2**2*SYMBY + 4*TG*(TG+UG)**(-1)*S4 + 2*TG - 
     +    4*TG**2*M2**2*(TG+UG)**(-1)*SYMBY + 4*TG**2*(TG+UG)**(-1)*
     +    SYMBY*S4**2 - 4*TG**2*SYMBY*S4 + 4*TG**3*M2*(TG+UG)**(-1)*
     +    SYMBY - 4*TG**3*(TG+UG)**(-1)*SYMBY*S4 - 2*TG**3*SYMBY - 4*
     +    T1**(-1)*U1**(-1)*M2**2*S4 + 4*T1**(-1)*U1**(-1)*M2**3 - 2*
     +    T1**(-1)*M2*S4 + 6*T1**(-1)*M2**2 + 8*U1**(-1)*M2*MG2 + 4*
     +    U1**(-1)*M2*SYMBY*S4**3 + 20*U1**(-1)*M2*S4 - 4*U1**(-1)*
     +    M2**2*SYMBY*S4**2 - 12*U1**(-1)*M2**2 - 4*U1**(-1)*S4**2 - 4*
     +    M2*SYMBY*S4**2 - 14*M2 + 4*M2**2*SYMBY*S4 + 2*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(35)*N*CF**2 * ( 8*S**2*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(36)*N*CF**2 * ( 8*S*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(37)*N*CF**2 * ( 16*S*U1**(-1)*M2*S4 - 16*S
     +    *M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(37)*N**2*CF * (  - 8*S**2*U1**(-1)*M2 )
     +
      M2QQH1 = M2QQH1 + ANG2(38)*N*CF**2 * (  - 8*S*TG*(TG+UG)**(-1) + 
     +    8*S*T1**(-1)*M2 - 8*S*T1**(-1)*S4 + 16*S*U1**(-1)*M2*
     +    (TG+UG)**(-1)*S4 + 24*S*U1**(-1)*M2 + 8*S*U1**(-1)*S4 - 8*S*
     +    M2*(TG+UG)**(-1) - 16*S**2*U1**(-1)*M2*(TG+UG)**(-1) - 8*S**2
     +    *(TG+UG)**(-1) + 4*TG - 16*T1**(-1)*U1**(-1)*M2**2*S4 - 16*
     +    T1**(-1)*U1**(-1)*M2**3 + 16*T1**(-1)*M2**2 + 8*U1**(-1)*M2*
     +    S4 + 8*U1**(-1)*M2**2 - 8*U1**(-1)*S4**2 + 4*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(38)*N**2*CF * ( 4*S*TG*(TG+UG)**(-1) - 2*S
     +    *T1**(-1)*M2 + 2*S*T1**(-1)*S4 + 8*S*U1**(-2)*M2*S4 - 8*S*
     +    U1**(-1)*M2*(TG+UG)**(-1)*S4 - 16*S*U1**(-1)*M2 - 4*S*
     +    U1**(-1)*S4 + 4*S*M2*(TG+UG)**(-1) + 2*S - 8*S**2*U1**(-2)*M2
     +     + 8*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 4*S**2*(TG+UG)**(-1) + 
     +    4*T1**(-1)*U1**(-1)*M2**2*S4 + 4*T1**(-1)*U1**(-1)*M2**3 - 4*
     +    T1**(-1)*M2**2 - 4*U1**(-1)*M2**2 )
     +
      M2QQH1 = M2QQH1 + ANG2(42)*N*CF**2 * (  - 16*S*U1**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 16*S*U1**(-1)*M2 - 16*S*MS2*(TG+UG)**(-1)
     +     + 16*S**2*U1**(-1)*M2*(TG+UG)**(-1) + 16*T1**(-1)*U1**(-1)*
     +    M2**2*S4 + 16*T1**(-1)*U1**(-1)*M2**3 + 16*T1**(-1)*M2*MG2 - 
     +    8*T1**(-1)*M2**2 + 16*T1**(-1)*MS2*S4 + 8*T1**(-1)*S4**2 - 16
     +    *U1**(-1)*M2**2 - 16*M2 + 8*M2**2*(TG+UG)**(-1) - 16*MS2 + 8*
     +    (TG+UG)**(-1)*S4**2 )
     +
      M2QQH1 = M2QQH1 + ANG2(42)*N**2*CF * (  - 8*S*U1**(-2)*M2*S4 + 8*
     +    S*U1**(-1)*M2*(TG+UG)**(-1)*S4 + 12*S*U1**(-1)*M2 + 8*S*MS2*
     +    (TG+UG)**(-1) + 8*S**2*U1**(-2)*M2 - 8*S**2*U1**(-1)*M2*
     +    (TG+UG)**(-1) - 4*T1**(-1)*U1**(-1)*M2**2*S4 - 4*T1**(-1)*
     +    U1**(-1)*M2**3 - 4*T1**(-1)*M2*MG2 + 2*T1**(-1)*M2**2 - 4*
     +    T1**(-1)*MS2*S4 - 2*T1**(-1)*S4**2 - 4*U1**(-1)*M2*S4 + 8*
     +    U1**(-1)*M2**2 + 4*U1**(-1)*S4**2 + 8*M2 - 4*M2**2*
     +    (TG+UG)**(-1) + 8*MS2 - 4*(TG+UG)**(-1)*S4**2 )
     +
      M2QQH1 = M2QQH1 + ANG2(51)*N*CF**2 * (  - 20 - 8*S*U1**(-1) - 8*S
     +    *(TG+UG)**(-1) - 8*T1**(-1)*M2 - 8*T1**(-1)*MS2 + 8*U1**(-1)*
     +    M2 + 16*U1**(-1)*MS2 + 8*U1**(-1)*S4 + 8*M2*(TG+UG)**(-1) + 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(51)*N**2*CF * ( 8 + 2*S*U1**(-1) + 4*S*
     +    (TG+UG)**(-1) + 4*T1**(-1)*M2 + 4*T1**(-1)*MS2 - 2*U1**(-1)*
     +    M2 - 4*U1**(-1)*MS2 - 2*U1**(-1)*S4 - 4*M2*(TG+UG)**(-1) - 4*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(58)*N*CF**2 * ( 24*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 8*S*TG*M2*SYMBY + 8*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 24*S*TG*(TG+UG)**(-1) - 8*S*TG*SYMBY*S4 - 16*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 32*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 16*
     +    S*TG**2*SYMBY + 16*S*TG**3*(TG+UG)**(-1)*SYMBY - 32*S*
     +    U1**(-1)*M2*SYMBY*S4**2 + 16*S*U1**(-1)*M2**2*SYMBY*S4 - 8*S*
     +    M2*(TG+UG)**(-1)*SYMBY*S4**2 + 24*S*M2*SYMBY*S4 + 8*S*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 8*S*M2**2*SYMBY + 8*S*(TG+UG)**(-1)*
     +    S4 - 8*S - 8*S**2*TG*M2*(TG+UG)**(-1)*SYMBY - 8*S**2*TG*
     +    (TG+UG)**(-1)*SYMBY*S4 + 16*S**2*TG**2*(TG+UG)**(-1)*SYMBY + 
     +    16*S**2*U1**(-1)*M2*SYMBY*S4 + 8*S**2*M2*(TG+UG)**(-1)*SYMBY*
     +    S4 - 8*S**2*M2**2*(TG+UG)**(-1)*SYMBY - 8*S**2*(TG+UG)**(-1)
     +     - 8*TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*TG*M2*(TG+UG)**(-1)
     +     - 8*TG*M2*SYMBY*S4 + 8*TG*M2**2*(TG+UG)**(-1)*SYMBY*S4 + 8*
     +    TG*M2**2*SYMBY + 8*TG*(TG+UG)**(-1)*S4 - 8*TG - 8*TG**2*M2*
     +    SYMBY )
     +
      M2QQH1 = M2QQH1 + ANG2(58)*N*CF**2 * (  - 8*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*
     +    TG**2*SYMBY*S4 + 8*TG**3*M2*(TG+UG)**(-1)*SYMBY - 8*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 + 16*U1**(-1)*M2*SYMBY*S4**3 - 8*
     +    U1**(-1)*M2*S4 - 16*U1**(-1)*M2**2*SYMBY*S4**2 + 8*U1**(-1)*
     +    S4**2 - 16*M2*SYMBY*S4**2 + 8*M2 + 16*M2**2*SYMBY*S4 - 8*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(58)*N**2*CF * (  - 12*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 2*S*TG*M2*SYMBY - 4*S*TG*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 12*S*TG*(TG+UG)**(-1) + 4*S*TG*
     +    SYMBY*S4 + 8*S*TG**2*M2*(TG+UG)**(-1)*SYMBY + 16*S*TG**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 6*S*TG**2*SYMBY - 8*S*TG**3*
     +    (TG+UG)**(-1)*SYMBY + 8*S*U1**(-1)*M2*SYMBY*S4**2 - 4*S*
     +    U1**(-1)*M2**2*SYMBY*S4 + 4*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 
     +    8*S*M2*SYMBY*S4 - 4*S*M2**2*(TG+UG)**(-1)*SYMBY*S4 + 4*S*
     +    M2**2*SYMBY - 4*S*(TG+UG)**(-1)*S4 + 2*S + 4*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBY + 4*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 - 8*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBY - 4*S**2*U1**(-1)*M2*SYMBY*S4
     +     - 4*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 + 4*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY + 4*S**2*(TG+UG)**(-1) + 4*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 4*TG*M2*(TG+UG)**(-1) - 4*TG*
     +    M2**2*(TG+UG)**(-1)*SYMBY*S4 + 2*TG*M2**2*SYMBY - 4*TG*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(58)*N**2*CF * ( 2*TG + 4*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 4*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 4*
     +    TG**2*SYMBY*S4 - 4*TG**3*M2*(TG+UG)**(-1)*SYMBY + 4*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 + 2*TG**3*SYMBY - 4*U1**(-1)*M2*SYMBY*
     +    S4**3 + 2*U1**(-1)*M2*S4 + 4*U1**(-1)*M2**2*SYMBY*S4**2 - 2*
     +    U1**(-1)*S4**2 + 4*M2*SYMBY*S4**2 - 2*M2 - 4*M2**2*SYMBY*S4
     +     + 2*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(73)*N*CF**2 * (  - 24*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 8*S*TG*M2*SYMBY - 8*S*TG*
     +    (TG+UG)**(-1)*SYMBY*S4**2 + 16*S*TG*(TG+UG)**(-1) + 8*S*TG*
     +    SYMBY*S4 + 16*S*TG**2*M2*(TG+UG)**(-1)*SYMBY + 32*S*TG**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 16*S*TG**2*SYMBY - 16*S*TG**3*
     +    (TG+UG)**(-1)*SYMBY + 32*S*U1**(-1)*M2*SYMBY*S4**2 - 8*S*
     +    U1**(-1)*M2 - 16*S*U1**(-1)*M2**2*SYMBY*S4 + 16*S*U1**(-1)*
     +    MS2 + 8*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 - 8*S*M2*(TG+UG)**(-1)
     +     - 24*S*M2*SYMBY*S4 - 8*S*M2**2*(TG+UG)**(-1)*SYMBY*S4 + 8*S*
     +    M2**2*SYMBY - 16*S*MS2*(TG+UG)**(-1) - 8*S*(TG+UG)**(-1)*S4
     +     - 12*S + 8*S**2*TG*M2*(TG+UG)**(-1)*SYMBY + 8*S**2*TG*
     +    (TG+UG)**(-1)*SYMBY*S4 - 16*S**2*TG**2*(TG+UG)**(-1)*SYMBY - 
     +    16*S**2*U1**(-1)*M2*SYMBY*S4 - 8*S**2*U1**(-1) - 8*S**2*M2*
     +    (TG+UG)**(-1)*SYMBY*S4 + 8*S**2*M2**2*(TG+UG)**(-1)*SYMBY + 8
     +    *TG*M2*(TG+UG)**(-1)*SYMBY*S4**2 + 8*TG*M2*(TG+UG)**(-1) + 8*
     +    TG*M2*SYMBY*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(73)*N*CF**2 * (  - 8*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 8*TG*M2**2*SYMBY - 8*TG*
     +    (TG+UG)**(-1)*S4 + 8*TG**2*M2*SYMBY + 8*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBY - 8*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 + 8*
     +    TG**2*SYMBY*S4 - 8*TG**3*M2*(TG+UG)**(-1)*SYMBY + 8*TG**3*
     +    (TG+UG)**(-1)*SYMBY*S4 - 24*T1**(-1)*M2*MG2 + 16*U1**(-1)*M2*
     +    MG2 - 16*U1**(-1)*M2*SYMBY*S4**3 + 16*U1**(-1)*M2*S4 + 16*
     +    U1**(-1)*M2**2*SYMBY*S4**2 + 16*M2*SYMBY*S4**2 - 28*M2 + 8*
     +    M2**2*(TG+UG)**(-1) - 16*M2**2*SYMBY*S4 + 8*MS2 + 8*
     +    (TG+UG)**(-1)*S4**2 - 8*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(73)*N**2*CF * ( 12*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBY*S4 - 2*S*TG*M2*SYMBY + 4*S*TG*(TG+UG)**(-1)*SYMBY*S4**2
     +     - 8*S*TG*(TG+UG)**(-1) - 4*S*TG*SYMBY*S4 - 8*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBY - 16*S*TG**2*(TG+UG)**(-1)*SYMBY*S4 + 6*S
     +    *TG**2*SYMBY + 8*S*TG**3*(TG+UG)**(-1)*SYMBY - 8*S*U1**(-1)*
     +    M2*SYMBY*S4**2 + 2*S*U1**(-1)*M2 + 4*S*U1**(-1)*M2**2*SYMBY*
     +    S4 - 4*S*U1**(-1)*MS2 - 4*S*M2*(TG+UG)**(-1)*SYMBY*S4**2 + 4*
     +    S*M2*(TG+UG)**(-1) + 8*S*M2*SYMBY*S4 + 4*S*M2**2*
     +    (TG+UG)**(-1)*SYMBY*S4 - 4*S*M2**2*SYMBY + 8*S*MS2*
     +    (TG+UG)**(-1) + 4*S*(TG+UG)**(-1)*S4 + 6*S - 4*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBY - 4*S**2*TG*(TG+UG)**(-1)*SYMBY*S4 + 8*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBY + 4*S**2*U1**(-1)*M2*SYMBY*S4
     +     + 2*S**2*U1**(-1) + 4*S**2*M2*(TG+UG)**(-1)*SYMBY*S4 - 4*
     +    S**2*M2**2*(TG+UG)**(-1)*SYMBY - 4*TG*M2*(TG+UG)**(-1)*SYMBY*
     +    S4**2 - 4*TG*M2*(TG+UG)**(-1) + 4*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBY*S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(73)*N**2*CF * (  - 2*TG*M2**2*SYMBY + 4*TG
     +    *(TG+UG)**(-1)*S4 + 2*TG - 4*TG**2*M2**2*(TG+UG)**(-1)*SYMBY
     +     + 4*TG**2*(TG+UG)**(-1)*SYMBY*S4**2 - 4*TG**2*SYMBY*S4 + 4*
     +    TG**3*M2*(TG+UG)**(-1)*SYMBY - 4*TG**3*(TG+UG)**(-1)*SYMBY*S4
     +     - 2*TG**3*SYMBY + 12*T1**(-1)*M2*MG2 - 4*U1**(-1)*M2*MG2 + 4
     +    *U1**(-1)*M2*SYMBY*S4**3 - 4*U1**(-1)*M2*S4 - 4*U1**(-1)*
     +    M2**2*SYMBY*S4**2 - 4*M2*SYMBY*S4**2 + 12*M2 - 4*M2**2*
     +    (TG+UG)**(-1) + 4*M2**2*SYMBY*S4 - 4*(TG+UG)**(-1)*S4**2 + 4*
     +    S4 )
     +
      M2QQH1 = M2QQH1 + ANG2(74)*N*CF**2 * (  - 16*T1**(-1)*M2**2*MG2
     +     + 8*M2*MG2 - 8*M2**2 )
     +
      M2QQH1 = M2QQH1 + ANG2(74)*N**2*CF * ( 4*TG*M2 + 4*TG**2*M2**2*
     +    SYMBY - 4*TG**3*M2*SYMBY + 8*T1**(-1)*M2**2*MG2 + 4*M2**2 )
     +
      M2QQH1 = M2QQH1 + COLO2(9)*N*CF**2*(S4+MS2) * (  - 8*S**(-1)*TG*
     +    (S+TG)**(-1)*S4**(-1) - 32*S**(-1)*U1**(-4)*M2**2*S4 + 32*
     +    S**(-1)*U1**(-3)*M2*S4 + 32*S**(-1)*U1**(-3)*M2**2 + 16*
     +    S**(-1)*U1**(-2)*M2*MG2*(S+TG)**(-1) - 16*S**(-1)*U1**(-2)*M2
     +    *MG2*S4**(-1) - 32*S**(-1)*U1**(-2)*M2 - 16*S**(-1)*U1**(-2)*
     +    M2**2*(S+TG)**(-1) - 16*S**(-1)*U1**(-2)*S4 - 16*S**(-1)*
     +    U1**(-1)*M2*MG2*(S+TG)**(-1)*S4**(-1) + 16*S**(-1)*U1**(-1)*
     +    M2*S4**(-1) + 16*S**(-1)*U1**(-1)*M2**2*(S+TG)**(-1)*S4**(-1)
     +     + 16*S**(-1)*U1**(-1) - 32*U1**(-4)*M2*MG2 + 32*U1**(-4)*
     +    M2**2 - 32*U1**(-3)*M2 - 16*U1**(-2)*M2*MG2*(S+TG)**(-1)*
     +    S4**(-1) + 16*U1**(-2)*M2**2*(S+TG)**(-1)*S4**(-1) + 16*
     +    U1**(-2) - 16*U1**(-1)*M2*(S+TG)**(-1)*S4**(-1) )
     +
      M2QQH1 = M2QQH1 + COLO2(9)*N**2*CF*(S4+MS2) * (  - 8*S*TG**(-1)*
     +    (S+TG)**(-1)*S4**(-1) - 32*TG**(-2)*U1**(-4)*M2*MG2*S4**2 + 
     +    64*TG**(-2)*U1**(-3)*M2*MG2*S4 - 48*TG**(-2)*U1**(-2)*M2*MG2
     +     + 16*TG**(-2)*U1**(-1)*M2*MG2*S4**(-1) + 32*TG**(-1)*
     +    U1**(-4)*M2*MG2*S4 - 32*TG**(-1)*U1**(-4)*M2**2*S4 - 32*
     +    TG**(-1)*U1**(-3)*M2*MG2 + 32*TG**(-1)*U1**(-3)*M2*S4 + 32*
     +    TG**(-1)*U1**(-3)*M2**2 + 16*TG**(-1)*U1**(-2)*M2*MG2*
     +    S4**(-1) - 32*TG**(-1)*U1**(-2)*M2 - 16*TG**(-1)*U1**(-2)*
     +    M2**2*S4**(-1) - 16*TG**(-1)*U1**(-2)*S4 + 16*TG**(-1)*
     +    U1**(-1)*M2*S4**(-1) + 16*TG**(-1)*U1**(-1) - 8*(S+TG)**(-1)*
     +    S4**(-1) )


      M2QQH2 = 0D0
      M2QQH2 = M2QQH2 + N*CF**2*(S4+MS2) * (  - 32*S**(-1)*T1**(-4)*
     +    M2**2*S4 + 32*S**(-1)*T1**(-3)*M2*S4 + 32*S**(-1)*T1**(-3)*
     +    M2**2 + 16*S**(-1)*T1**(-2)*M2*MG2*(S+UG)**(-1) - 16*S**(-1)*
     +    T1**(-2)*M2*MG2*S4**(-1) - 32*S**(-1)*T1**(-2)*M2 - 16*
     +    S**(-1)*T1**(-2)*M2**2*(S+UG)**(-1) - 16*S**(-1)*T1**(-1)*M2*
     +    MG2*(S+UG)**(-1)*S4**(-1) + 16*S**(-1)*T1**(-1)*M2*S4**(-1)
     +     + 16*S**(-1)*T1**(-1)*M2**2*(S+UG)**(-1)*S4**(-1) - 8*
     +    S**(-1)*S4**(-1) + 96*S*TG*M2*(TG+UG)**(-1)*SYMBU*SYMBX + 128
     +    *S*TG*M2*SYMBU*SYMBX*S4**(-1) - 208*S*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBU*SYMBX*S4**(-1) + 16*S*TG*(TG+UG)**(-1)*SYMBU*S4**(-1)
     +     + 32*S*TG*(TG+UG)**(-1)*SYMBX*S4**(-1) - 48*S*TG*SYMBU*SYMBX
     +     - 176*S*TG**2*M2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) + 48*S*
     +    TG**2*(TG+UG)**(-1)*SYMBU*SYMBX + 48*S*TG**2*SYMBU*SYMBX*
     +    S4**(-1) - 48*S*TG**3*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 16
     +    *S*UG**(-1)*M2*(TG+UG)**(-1)*SYMBX + 16*S*UG**(-1)*M2*SYMBX*
     +    S4**(-1) )
     +
      M2QQH2 = M2QQH2 + N*CF**2*(S4+MS2) * ( 16*S*UG**(-1)*
     +    (TG+UG)**(-1)*S4**(-1) + 16*S*U1**(-1)*(S+UG)**(-1)*S4**(-1)
     +     + 16*S*M2*(TG+UG)**(-1)*SYMBU*S4**(-1) + 48*S*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**(-1) - 48*S*M2*SYMBU*SYMBX + 48*S*
     +    M2**2*(TG+UG)**(-1)*SYMBU*SYMBX + 80*S*M2**2*SYMBU*SYMBX*
     +    S4**(-1) - 80*S*M2**3*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 32
     +    *S*SYMBX*S4**(-1) - 64*S**2*TG*M2*(TG+UG)**(-1)*SYMBU*SYMBX*
     +    S4**(-1) + 32*S**2*TG*SYMBU*SYMBX*S4**(-1) - 32*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) + 16*S**2*UG**(-1)*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**(-1) + 32*S**2*M2*SYMBU*SYMBX*
     +    S4**(-1) - 32*S**2*M2**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1)
     +     - 32*TG*M2*(TG+UG)**(-1)*SYMBU*SYMBX*S4 + 32*TG*M2*
     +    (TG+UG)**(-1)*SYMBU*S4**(-1) + 48*TG*M2*(TG+UG)**(-1)*SYMBX*
     +    S4**(-1) - 80*TG*M2*SYMBU*SYMBX + 128*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBU*SYMBX + 80*TG*M2**2*SYMBU*SYMBX*S4**(-1) - 112*TG*M2**3
     +    *(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) )
     +
      M2QQH2 = M2QQH2 + N*CF**2*(S4+MS2) * (  - 16*TG*(TG+UG)**(-1)*
     +    SYMBU - 16*TG*(TG+UG)**(-1)*SYMBX + 16*TG*SYMBU*SYMBX*S4 - 16
     +    *TG*SYMBU*S4**(-1) - 16*TG*SYMBX*S4**(-1) + 112*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBU*SYMBX + 64*TG**2*M2*SYMBU*SYMBX*S4**(-1)
     +     - 144*TG**2*M2**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 16*
     +    TG**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4 + 16*TG**2*(TG+UG)**(-1)*
     +    SYMBU*S4**(-1) + 16*TG**2*(TG+UG)**(-1)*SYMBX*S4**(-1) - 32*
     +    TG**2*SYMBU*SYMBX - 80*TG**3*M2*(TG+UG)**(-1)*SYMBU*SYMBX*
     +    S4**(-1) + 32*TG**3*(TG+UG)**(-1)*SYMBU*SYMBX + 16*TG**3*
     +    SYMBU*SYMBX*S4**(-1) - 16*TG**4*(TG+UG)**(-1)*SYMBU*SYMBX*
     +    S4**(-1) - 32*T1**(-4)*M2*MG2 + 32*T1**(-4)*M2**2 - 32*
     +    T1**(-3)*M2 - 16*T1**(-2)*M2*MG2*(S+UG)**(-1)*S4**(-1) + 16*
     +    T1**(-2)*M2**2*(S+UG)**(-1)*S4**(-1) - 16*T1**(-1)*M2*
     +    (S+UG)**(-1)*S4**(-1) - 16*U1**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 16*U1**(-1)*S4**(-1) - 16*M2*(TG+UG)**(-1)*SYMBU - 16*M2*
     +    (TG+UG)**(-1)*SYMBX )
     +
      M2QQH2 = M2QQH2 + N*CF**2*(S4+MS2) * ( 16*M2*SYMBU*SYMBX*S4 - 16*
     +    M2*SYMBU*S4**(-1) - 32*M2*SYMBX*S4**(-1) - 16*M2**2*
     +    (TG+UG)**(-1)*SYMBU*SYMBX*S4 + 16*M2**2*(TG+UG)**(-1)*SYMBU*
     +    S4**(-1) + 32*M2**2*(TG+UG)**(-1)*SYMBX*S4**(-1) - 48*M2**2*
     +    SYMBU*SYMBX + 48*M2**3*(TG+UG)**(-1)*SYMBU*SYMBX + 32*M2**3*
     +    SYMBU*SYMBX*S4**(-1) - 32*M2**4*(TG+UG)**(-1)*SYMBU*SYMBX*
     +    S4**(-1) + 24*(S+UG)**(-1)*S4**(-1) + 16*SYMBU + 16*SYMBX )
     +
      M2QQH2 = M2QQH2 + N**2*CF*(S4+MS2) * (  - 48*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBU*SYMBX - 44*S*TG*M2*SYMBU*SYMBX*S4**(-1)
     +     + 104*S*TG*M2**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 8*S*TG
     +    *(TG+UG)**(-1)*SYMBU*S4**(-1) - 16*S*TG*(TG+UG)**(-1)*SYMBX*
     +    S4**(-1) + 16*S*TG*SYMBU*SYMBX + 88*S*TG**2*M2*(TG+UG)**(-1)*
     +    SYMBU*SYMBX*S4**(-1) - 24*S*TG**2*(TG+UG)**(-1)*SYMBU*SYMBX
     +     - 16*S*TG**2*SYMBU*SYMBX*S4**(-1) + 24*S*TG**3*(TG+UG)**(-1)
     +    *SYMBU*SYMBX*S4**(-1) + 8*S*UG**(-1)*M2*(TG+UG)**(-1)*SYMBX
     +     - 8*S*UG**(-1)*M2*SYMBX*S4**(-1) - 8*S*UG**(-1)*
     +    (TG+UG)**(-1)*S4**(-1) - 8*S*U1**(-1)*(S+UG)**(-1)*S4**(-1)
     +     - 8*S*M2*(TG+UG)**(-1)*SYMBU*S4**(-1) - 24*S*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**(-1) + 16*S*M2*SYMBU*SYMBX - 24*S*
     +    M2**2*(TG+UG)**(-1)*SYMBU*SYMBX - 28*S*M2**2*SYMBU*SYMBX*
     +    S4**(-1) + 40*S*M2**3*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) + 12
     +    *S*SYMBX*S4**(-1) + 32*S**2*TG*M2*(TG+UG)**(-1)*SYMBU*SYMBX*
     +    S4**(-1) )
     +
      M2QQH2 = M2QQH2 + N**2*CF*(S4+MS2) * (  - 12*S**2*TG*SYMBU*SYMBX*
     +    S4**(-1) + 16*S**2*TG**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1)
     +     - 8*S**2*UG**(-1)*M2*(TG+UG)**(-1)*SYMBX*S4**(-1) - 12*S**2*
     +    M2*SYMBU*SYMBX*S4**(-1) + 16*S**2*M2**2*(TG+UG)**(-1)*SYMBU*
     +    SYMBX*S4**(-1) + 16*TG*M2*(TG+UG)**(-1)*SYMBU*SYMBX*S4 - 16*
     +    TG*M2*(TG+UG)**(-1)*SYMBU*S4**(-1) - 24*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**(-1) + 20*TG*M2*SYMBU*SYMBX - 64*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBU*SYMBX - 20*TG*M2**2*SYMBU*SYMBX*S4**(-1)
     +     + 56*TG*M2**3*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) + 8*TG*
     +    (TG+UG)**(-1)*SYMBU + 8*TG*(TG+UG)**(-1)*SYMBX - 4*TG*SYMBU*
     +    SYMBX*S4 + 4*TG*SYMBU*S4**(-1) + 4*TG*SYMBX*S4**(-1) - 56*
     +    TG**2*M2*(TG+UG)**(-1)*SYMBU*SYMBX - 16*TG**2*M2*SYMBU*SYMBX*
     +    S4**(-1) + 72*TG**2*M2**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1)
     +     + 8*TG**2*(TG+UG)**(-1)*SYMBU*SYMBX*S4 - 8*TG**2*
     +    (TG+UG)**(-1)*SYMBU*S4**(-1) - 8*TG**2*(TG+UG)**(-1)*SYMBX*
     +    S4**(-1) )
     +
      M2QQH2 = M2QQH2 + N**2*CF*(S4+MS2) * ( 8*TG**2*SYMBU*SYMBX + 40*
     +    TG**3*M2*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 16*TG**3*
     +    (TG+UG)**(-1)*SYMBU*SYMBX - 4*TG**3*SYMBU*SYMBX*S4**(-1) + 8*
     +    TG**4*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 32*UG**(-2)*
     +    T1**(-4)*M2*MG2*S4**2 + 64*UG**(-2)*T1**(-3)*M2*MG2*S4 - 48*
     +    UG**(-2)*T1**(-2)*M2*MG2 + 16*UG**(-2)*T1**(-1)*M2*MG2*
     +    S4**(-1) + 32*UG**(-1)*T1**(-4)*M2*MG2*S4 - 32*UG**(-1)*
     +    T1**(-4)*M2**2*S4 - 32*UG**(-1)*T1**(-3)*M2*MG2 + 32*UG**(-1)
     +    *T1**(-3)*M2*S4 + 32*UG**(-1)*T1**(-3)*M2**2 + 16*UG**(-1)*
     +    T1**(-2)*M2*MG2*S4**(-1) - 32*UG**(-1)*T1**(-2)*M2 - 16*
     +    UG**(-1)*T1**(-2)*M2**2*S4**(-1) + 16*UG**(-1)*T1**(-1)*M2*
     +    S4**(-1) - 4*UG**(-1)*U1**(-1)*M2*S4**(-1) - 4*UG**(-1)*
     +    S4**(-1) + 8*U1**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 4*U1**(-1)*
     +    S4**(-1) + 8*M2*(TG+UG)**(-1)*SYMBU + 8*M2*(TG+UG)**(-1)*
     +    SYMBX - 4*M2*SYMBU*SYMBX*S4 + 4*M2*SYMBU*S4**(-1) + 8*M2*
     +    SYMBX*S4**(-1) )
     +
      M2QQH2 = M2QQH2 + N**2*CF*(S4+MS2) * ( 8*M2**2*(TG+UG)**(-1)*
     +    SYMBU*SYMBX*S4 - 8*M2**2*(TG+UG)**(-1)*SYMBU*S4**(-1) - 16*
     +    M2**2*(TG+UG)**(-1)*SYMBX*S4**(-1) + 12*M2**2*SYMBU*SYMBX - 
     +    24*M2**3*(TG+UG)**(-1)*SYMBU*SYMBX - 8*M2**3*SYMBU*SYMBX*
     +    S4**(-1) + 16*M2**4*(TG+UG)**(-1)*SYMBU*SYMBX*S4**(-1) - 8*
     +    (S+UG)**(-1)*S4**(-1) - 4*SYMBU - 4*SYMBX )
     +
      M2QQH2 = M2QQH2 + ANG2(13)*N*CF**2 * ( 4*M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(14)*N*CF**2 * ( 4 + 4*S*T1**(-1) - 8*S*
     +    (TG+UG)**(-1) - 4*T1**(-1)*M2 - 8*T1**(-1)*MS2 - 4*T1**(-1)*
     +    S4 + 16*U1**(-1)*M2 + 16*U1**(-1)*MS2 + 8*M2*(TG+UG)**(-1) + 
     +    8*(TG+UG)**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(14)*N**2*CF * (  - 2*S*T1**(-1) + 4*S*
     +    (TG+UG)**(-1) + 2*T1**(-1)*M2 + 4*T1**(-1)*MS2 + 2*T1**(-1)*
     +    S4 - 4*U1**(-1)*M2 - 4*U1**(-1)*MS2 - 4*M2*(TG+UG)**(-1) - 4*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(39)*N*CF**2 * ( 8*S*M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(40)*N*CF**2 * ( 8*S*TG*(TG+UG)**(-1) + 16*
     +    S*T1**(-1)*M2*(TG+UG)**(-1)*S4 + 24*S*T1**(-1)*M2 + 8*S*
     +    T1**(-1)*S4 + 8*S*U1**(-1)*M2 - 8*S*U1**(-1)*S4 - 8*S*
     +    (TG+UG)**(-1)*S4 - 4*S - 16*S**2*T1**(-1)*M2*(TG+UG)**(-1) - 
     +    4*TG - 16*T1**(-1)*U1**(-1)*M2**2*S4 - 16*T1**(-1)*U1**(-1)*
     +    M2**3 + 8*T1**(-1)*M2*S4 + 8*T1**(-1)*M2**2 - 8*T1**(-1)*
     +    S4**2 + 16*U1**(-1)*M2**2 - 4*M2 + 8*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(40)*N**2*CF * (  - 4*S*TG*(TG+UG)**(-1) + 
     +    8*S*T1**(-2)*M2*S4 - 8*S*T1**(-1)*M2*(TG+UG)**(-1)*S4 - 16*S*
     +    T1**(-1)*M2 - 4*S*T1**(-1)*S4 - 2*S*U1**(-1)*M2 + 2*S*
     +    U1**(-1)*S4 + 4*S*(TG+UG)**(-1)*S4 + 2*S - 8*S**2*T1**(-2)*M2
     +     + 8*S**2*T1**(-1)*M2*(TG+UG)**(-1) + 4*T1**(-1)*U1**(-1)*
     +    M2**2*S4 + 4*T1**(-1)*U1**(-1)*M2**3 - 4*T1**(-1)*M2**2 - 4*
     +    U1**(-1)*M2**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(51)*N*CF**2 * (  - 20 - 8*S*T1**(-1) - 8*S
     +    *(TG+UG)**(-1) + 8*T1**(-1)*M2 + 16*T1**(-1)*MS2 + 8*T1**(-1)
     +    *S4 - 8*U1**(-1)*M2 - 8*U1**(-1)*MS2 + 8*M2*(TG+UG)**(-1) + 8
     +    *(TG+UG)**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(51)*N**2*CF * ( 8 + 2*S*T1**(-1) + 4*S*
     +    (TG+UG)**(-1) - 2*T1**(-1)*M2 - 4*T1**(-1)*MS2 - 2*T1**(-1)*
     +    S4 + 4*U1**(-1)*M2 + 4*U1**(-1)*MS2 - 4*M2*(TG+UG)**(-1) - 4*
     +    (TG+UG)**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N*CF**2 * ( 136*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4 + 24*S*TG*M2*SYMBX - 144*S*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBX - 24*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 + 24*S*TG*
     +    (TG+UG)**(-1) - 40*S*TG*SYMBX*S4 - 88*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX + 40*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 + 16*
     +    S*TG**2*SYMBX - 16*S*TG**3*(TG+UG)**(-1)*SYMBX - 32*S*
     +    T1**(-1)*M2*SYMBX*S4**2 + 16*S*T1**(-1)*M2**2*SYMBX*S4 - 40*S
     +    *M2*(TG+UG)**(-1)*SYMBX*S4**2 + 32*S*M2*(TG+UG)**(-1) + 112*S
     +    *M2**2*(TG+UG)**(-1)*SYMBX*S4 - 8*S*M2**2*SYMBX - 72*S*M2**3*
     +    (TG+UG)**(-1)*SYMBX - 24*S*(TG+UG)**(-1)*S4 + 24*S*SYMBX*
     +    S4**2 - 112*S**2*TG*M2*(TG+UG)**(-1)*SYMBX + 32*S**2*TG*
     +    (TG+UG)**(-1)*SYMBX*S4 + 32*S**2*TG*SYMBX - 32*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBX + 16*S**2*T1**(-1)*M2*SYMBX*S4 + 64*S**2*
     +    M2*(TG+UG)**(-1)*SYMBX*S4 + 32*S**2*M2*SYMBX - 96*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 16*S**2*(TG+UG)**(-1) - 32*S**2*SYMBX*
     +    S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N*CF**2 * (  - 16*S**3*TG*
     +    (TG+UG)**(-1)*SYMBX - 32*S**3*M2*(TG+UG)**(-1)*SYMBX + 16*
     +    S**3*SYMBX - 48*TG*M2*(TG+UG)**(-1)*SYMBX*S4**2 + 8*TG*M2*
     +    (TG+UG)**(-1) + 8*TG*M2*SYMBX*S4 + 80*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBX*S4 - 24*TG*M2**2*SYMBX - 40*TG*M2**3*(TG+UG)**(-1)*
     +    SYMBX + 8*TG*(TG+UG)**(-1)*SYMBX*S4**3 - 8*TG*(TG+UG)**(-1)*
     +    S4 + 16*TG*SYMBX*S4**2 + 8*TG + 48*TG**2*M2*(TG+UG)**(-1)*
     +    SYMBX*S4 - 8*TG**2*M2*SYMBX - 32*TG**2*M2**2*(TG+UG)**(-1)*
     +    SYMBX - 16*TG**2*(TG+UG)**(-1)*SYMBX*S4**2 - 8*TG**2*SYMBX*S4
     +     - 8*TG**3*M2*(TG+UG)**(-1)*SYMBX + 8*TG**3*(TG+UG)**(-1)*
     +    SYMBX*S4 + 16*T1**(-1)*M2*SYMBX*S4**3 - 8*T1**(-1)*M2*S4 - 16
     +    *T1**(-1)*M2**2*SYMBX*S4**2 + 8*T1**(-1)*S4**2 + 8*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**3 - 16*M2*(TG+UG)**(-1)*S4 - 16*M2*
     +    SYMBX*S4**2 + 16*M2 - 32*M2**2*(TG+UG)**(-1)*SYMBX*S4**2 + 8*
     +    M2**2*(TG+UG)**(-1) + 40*M2**2*SYMBX*S4 + 40*M2**3*
     +    (TG+UG)**(-1)*SYMBX*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N*CF**2 * (  - 16*M2**3*SYMBX - 16*
     +    M2**4*(TG+UG)**(-1)*SYMBX + 8*(TG+UG)**(-1)*S4**2 - 8*SYMBX*
     +    S4**3 - 16*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N**2*CF * (  - 68*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 - 26*S*TG*M2*SYMBX + 72*S*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 12*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 - 12*
     +    S*TG*(TG+UG)**(-1) + 28*S*TG*SYMBX*S4 + 44*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX - 20*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 - 12*
     +    S*TG**2*SYMBX + 8*S*TG**3*(TG+UG)**(-1)*SYMBX + 8*S*T1**(-1)*
     +    M2*SYMBX*S4**2 - 4*S*T1**(-1)*M2**2*SYMBX*S4 + 20*S*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**2 - 16*S*M2*(TG+UG)**(-1) + 22*S*M2*
     +    SYMBX*S4 - 56*S*M2**2*(TG+UG)**(-1)*SYMBX*S4 - 12*S*M2**2*
     +    SYMBX + 36*S*M2**3*(TG+UG)**(-1)*SYMBX + 12*S*(TG+UG)**(-1)*
     +    S4 - 16*S*SYMBX*S4**2 + 56*S**2*TG*M2*(TG+UG)**(-1)*SYMBX - 
     +    16*S**2*TG*(TG+UG)**(-1)*SYMBX*S4 - 18*S**2*TG*SYMBX + 16*
     +    S**2*TG**2*(TG+UG)**(-1)*SYMBX - 4*S**2*T1**(-1)*M2*SYMBX*S4
     +     - 32*S**2*M2*(TG+UG)**(-1)*SYMBX*S4 - 20*S**2*M2*SYMBX + 48*
     +    S**2*M2**2*(TG+UG)**(-1)*SYMBX - 8*S**2*(TG+UG)**(-1) + 18*
     +    S**2*SYMBX*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N**2*CF * ( 8*S**3*TG*(TG+UG)**(-1)*
     +    SYMBX + 16*S**3*M2*(TG+UG)**(-1)*SYMBX - 8*S**3*SYMBX + 24*TG
     +    *M2*(TG+UG)**(-1)*SYMBX*S4**2 - 4*TG*M2*(TG+UG)**(-1) + 20*TG
     +    *M2*SYMBX*S4 - 40*TG*M2**2*(TG+UG)**(-1)*SYMBX*S4 - 8*TG*
     +    M2**2*SYMBX + 20*TG*M2**3*(TG+UG)**(-1)*SYMBX - 4*TG*
     +    (TG+UG)**(-1)*SYMBX*S4**3 + 4*TG*(TG+UG)**(-1)*S4 - 14*TG*
     +    SYMBX*S4**2 - 2*TG - 24*TG**2*M2*(TG+UG)**(-1)*SYMBX*S4 - 6*
     +    TG**2*M2*SYMBX + 16*TG**2*M2**2*(TG+UG)**(-1)*SYMBX + 8*TG**2
     +    *(TG+UG)**(-1)*SYMBX*S4**2 + 10*TG**2*SYMBX*S4 + 4*TG**3*M2*
     +    (TG+UG)**(-1)*SYMBX - 4*TG**3*(TG+UG)**(-1)*SYMBX*S4 - 2*
     +    TG**3*SYMBX - 4*T1**(-1)*M2*SYMBX*S4**3 + 2*T1**(-1)*M2*S4 + 
     +    4*T1**(-1)*M2**2*SYMBX*S4**2 - 2*T1**(-1)*S4**2 - 4*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**3 + 8*M2*(TG+UG)**(-1)*S4 - 10*M2*
     +    SYMBX*S4**2 - 4*M2 + 16*M2**2*(TG+UG)**(-1)*SYMBX*S4**2 - 4*
     +    M2**2*(TG+UG)**(-1) + 8*M2**2*SYMBX*S4 - 20*M2**3*
     +    (TG+UG)**(-1)*SYMBX*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(59)*N**2*CF * (  - 4*M2**3*SYMBX + 8*M2**4
     +    *(TG+UG)**(-1)*SYMBX - 4*(TG+UG)**(-1)*S4**2 + 6*SYMBX*S4**3
     +     + 4*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(75)*N*CF**2 * (  - 16*S*T1**(-1)*M2*
     +    (TG+UG)**(-1)*S4 - 16*S*T1**(-1)*M2 - 16*S*MS2*(TG+UG)**(-1)
     +     + 16*S**2*T1**(-1)*M2*(TG+UG)**(-1) + 16*T1**(-1)*U1**(-1)*
     +    M2**2*S4 + 16*T1**(-1)*U1**(-1)*M2**3 - 16*T1**(-1)*M2**2 + 
     +    16*U1**(-1)*M2*MG2 - 8*U1**(-1)*M2**2 + 16*U1**(-1)*MS2*S4 + 
     +    8*U1**(-1)*S4**2 - 16*M2 + 8*M2**2*(TG+UG)**(-1) - 16*MS2 + 8
     +    *(TG+UG)**(-1)*S4**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(75)*N**2*CF * (  - 8*S*T1**(-2)*M2*S4 + 8*
     +    S*T1**(-1)*M2*(TG+UG)**(-1)*S4 + 12*S*T1**(-1)*M2 + 8*S*MS2*
     +    (TG+UG)**(-1) + 8*S**2*T1**(-2)*M2 - 8*S**2*T1**(-1)*M2*
     +    (TG+UG)**(-1) - 4*T1**(-1)*U1**(-1)*M2**2*S4 - 4*T1**(-1)*
     +    U1**(-1)*M2**3 - 4*T1**(-1)*M2*S4 + 8*T1**(-1)*M2**2 + 4*
     +    T1**(-1)*S4**2 - 4*U1**(-1)*M2*MG2 + 2*U1**(-1)*M2**2 - 4*
     +    U1**(-1)*MS2*S4 - 2*U1**(-1)*S4**2 + 8*M2 - 4*M2**2*
     +    (TG+UG)**(-1) + 8*MS2 - 4*(TG+UG)**(-1)*S4**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(77)*N*CF**2 * (  - 16*U1**(-1)*M2**2*MG2
     +     + 8*M2*MG2 - 8*M2**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(77)*N**2*CF * (  - 24*S*TG*M2*SYMBX*S4 + 
     +    32*S*TG*M2**2*SYMBX + 12*S*TG**2*M2*SYMBX + 12*S*M2*SYMBX*
     +    S4**2 - 4*S*M2 - 32*S*M2**2*SYMBX*S4 + 20*S*M2**3*SYMBX + 12*
     +    S**2*TG*M2*SYMBX - 12*S**2*M2*SYMBX*S4 + 16*S**2*M2**2*SYMBX
     +     + 4*S**3*M2*SYMBX + 12*TG*M2*SYMBX*S4**2 - 4*TG*M2 - 32*TG*
     +    M2**2*SYMBX*S4 + 20*TG*M2**3*SYMBX - 12*TG**2*M2*SYMBX*S4 + 
     +    16*TG**2*M2**2*SYMBX + 4*TG**3*M2*SYMBX + 8*U1**(-1)*M2**2*
     +    MG2 - 4*M2*SYMBX*S4**3 + 4*M2*S4 + 16*M2**2*SYMBX*S4**2 - 20*
     +    M2**3*SYMBX*S4 + 8*M2**4*SYMBX )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N*CF**2 * (  - 136*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 - 24*S*TG*M2*SYMBX + 144*S*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 24*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 - 16*
     +    S*TG*(TG+UG)**(-1) + 40*S*TG*SYMBX*S4 + 88*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX - 40*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 - 16*
     +    S*TG**2*SYMBX + 16*S*TG**3*(TG+UG)**(-1)*SYMBX + 32*S*
     +    T1**(-1)*M2*SYMBX*S4**2 - 8*S*T1**(-1)*M2 - 16*S*T1**(-1)*
     +    M2**2*SYMBX*S4 + 16*S*T1**(-1)*MS2 + 40*S*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 - 32*S*M2*(TG+UG)**(-1) - 112*S*M2**2*
     +    (TG+UG)**(-1)*SYMBX*S4 + 8*S*M2**2*SYMBX + 72*S*M2**3*
     +    (TG+UG)**(-1)*SYMBX - 16*S*MS2*(TG+UG)**(-1) + 16*S*
     +    (TG+UG)**(-1)*S4 - 24*S*SYMBX*S4**2 - 12*S + 112*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBX - 32*S**2*TG*(TG+UG)**(-1)*SYMBX*S4 - 32*
     +    S**2*TG*SYMBX + 32*S**2*TG**2*(TG+UG)**(-1)*SYMBX - 16*S**2*
     +    T1**(-1)*M2*SYMBX*S4 - 8*S**2*T1**(-1) - 64*S**2*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N*CF**2 * (  - 32*S**2*M2*SYMBX + 96*
     +    S**2*M2**2*(TG+UG)**(-1)*SYMBX - 16*S**2*(TG+UG)**(-1) + 32*
     +    S**2*SYMBX*S4 + 16*S**3*TG*(TG+UG)**(-1)*SYMBX + 32*S**3*M2*
     +    (TG+UG)**(-1)*SYMBX - 16*S**3*SYMBX + 48*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 - 8*TG*M2*(TG+UG)**(-1) - 8*TG*M2*SYMBX*S4 - 80*
     +    TG*M2**2*(TG+UG)**(-1)*SYMBX*S4 + 24*TG*M2**2*SYMBX + 40*TG*
     +    M2**3*(TG+UG)**(-1)*SYMBX - 8*TG*(TG+UG)**(-1)*SYMBX*S4**3 + 
     +    8*TG*(TG+UG)**(-1)*S4 - 16*TG*SYMBX*S4**2 - 48*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 + 8*TG**2*M2*SYMBX + 32*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 16*TG**2*(TG+UG)**(-1)*SYMBX*S4**2 + 8*
     +    TG**2*SYMBX*S4 + 8*TG**3*M2*(TG+UG)**(-1)*SYMBX - 8*TG**3*
     +    (TG+UG)**(-1)*SYMBX*S4 + 16*T1**(-1)*M2*MG2 - 16*T1**(-1)*M2*
     +    SYMBX*S4**3 + 16*T1**(-1)*M2*S4 + 16*T1**(-1)*M2**2*SYMBX*
     +    S4**2 - 24*U1**(-1)*M2*MG2 - 8*M2*(TG+UG)**(-1)*SYMBX*S4**3
     +     + 16*M2*(TG+UG)**(-1)*S4 + 16*M2*SYMBX*S4**2 - 28*M2 + 32*
     +    M2**2*(TG+UG)**(-1)*SYMBX*S4**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N*CF**2 * (  - 40*M2**2*SYMBX*S4 - 40*
     +    M2**3*(TG+UG)**(-1)*SYMBX*S4 + 16*M2**3*SYMBX + 16*M2**4*
     +    (TG+UG)**(-1)*SYMBX + 8*MS2 + 8*SYMBX*S4**3 - 8*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N**2*CF * ( 68*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4 + 26*S*TG*M2*SYMBX - 72*S*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBX - 12*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 + 8*S*TG*
     +    (TG+UG)**(-1) - 28*S*TG*SYMBX*S4 - 44*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX + 20*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 + 12*
     +    S*TG**2*SYMBX - 8*S*TG**3*(TG+UG)**(-1)*SYMBX - 8*S*T1**(-1)*
     +    M2*SYMBX*S4**2 + 2*S*T1**(-1)*M2 + 4*S*T1**(-1)*M2**2*SYMBX*
     +    S4 - 4*S*T1**(-1)*MS2 - 20*S*M2*(TG+UG)**(-1)*SYMBX*S4**2 + 
     +    16*S*M2*(TG+UG)**(-1) - 22*S*M2*SYMBX*S4 + 56*S*M2**2*
     +    (TG+UG)**(-1)*SYMBX*S4 + 12*S*M2**2*SYMBX - 36*S*M2**3*
     +    (TG+UG)**(-1)*SYMBX + 8*S*MS2*(TG+UG)**(-1) - 8*S*
     +    (TG+UG)**(-1)*S4 + 16*S*SYMBX*S4**2 + 4*S - 56*S**2*TG*M2*
     +    (TG+UG)**(-1)*SYMBX + 16*S**2*TG*(TG+UG)**(-1)*SYMBX*S4 + 18*
     +    S**2*TG*SYMBX - 16*S**2*TG**2*(TG+UG)**(-1)*SYMBX + 4*S**2*
     +    T1**(-1)*M2*SYMBX*S4 + 2*S**2*T1**(-1) + 32*S**2*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N**2*CF * ( 20*S**2*M2*SYMBX - 48*S**2
     +    *M2**2*(TG+UG)**(-1)*SYMBX + 8*S**2*(TG+UG)**(-1) - 18*S**2*
     +    SYMBX*S4 - 8*S**3*TG*(TG+UG)**(-1)*SYMBX - 16*S**3*M2*
     +    (TG+UG)**(-1)*SYMBX + 8*S**3*SYMBX - 24*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 + 4*TG*M2*(TG+UG)**(-1) - 20*TG*M2*SYMBX*S4 + 40*
     +    TG*M2**2*(TG+UG)**(-1)*SYMBX*S4 + 8*TG*M2**2*SYMBX - 20*TG*
     +    M2**3*(TG+UG)**(-1)*SYMBX + 4*TG*(TG+UG)**(-1)*SYMBX*S4**3 - 
     +    4*TG*(TG+UG)**(-1)*S4 + 14*TG*SYMBX*S4**2 - 2*TG + 24*TG**2*
     +    M2*(TG+UG)**(-1)*SYMBX*S4 + 6*TG**2*M2*SYMBX - 16*TG**2*M2**2
     +    *(TG+UG)**(-1)*SYMBX - 8*TG**2*(TG+UG)**(-1)*SYMBX*S4**2 - 10
     +    *TG**2*SYMBX*S4 - 4*TG**3*M2*(TG+UG)**(-1)*SYMBX + 4*TG**3*
     +    (TG+UG)**(-1)*SYMBX*S4 + 2*TG**3*SYMBX - 4*T1**(-1)*M2*MG2 + 
     +    4*T1**(-1)*M2*SYMBX*S4**3 - 4*T1**(-1)*M2*S4 - 4*T1**(-1)*
     +    M2**2*SYMBX*S4**2 + 12*U1**(-1)*M2*MG2 + 4*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**3 - 8*M2*(TG+UG)**(-1)*S4 + 10*M2*SYMBX*S4**2 + 10*
     +    M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(78)*N**2*CF * (  - 16*M2**2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 - 8*M2**2*SYMBX*S4 + 20*M2**3*(TG+UG)**(-1)*SYMBX
     +    *S4 + 4*M2**3*SYMBX - 8*M2**4*(TG+UG)**(-1)*SYMBX - 6*SYMBX*
     +    S4**3 + 6*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(79)*N*CF**2 * ( 4*U1**(-2)*M2**3 + 8*
     +    U1**(-1)*M2*MG2 - 8*U1**(-1)*M2**2 + 8*M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(80)*N*CF**2 * ( 12 + 8*S*T1**(-1) + 8*
     +    T1**(-1)*U1**(-1)*M2*S4 + 16*T1**(-1)*U1**(-1)*MS2*S4 + 8*
     +    T1**(-1)*U1**(-1)*S4**2 - 16*T1**(-1)*M2 - 16*T1**(-1)*MS2 - 
     +    8*T1**(-1)*S4 - 4*U1**(-2)*M2*S4 + 4*U1**(-2)*M2**2 + 4*
     +    U1**(-1)*M2 + 8*U1**(-1)*MS2 + 4*U1**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(80)*N**2*CF * (  - 4 - 4*S*T1**(-1) - 4*
     +    T1**(-1)*U1**(-1)*M2*S4 - 8*T1**(-1)*U1**(-1)*MS2*S4 - 4*
     +    T1**(-1)*U1**(-1)*S4**2 + 8*T1**(-1)*M2 + 8*T1**(-1)*MS2 + 4*
     +    T1**(-1)*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(81)*N*CF**2 * ( 8*T1**(-2)*M2*S4**2 - 16*
     +    T1**(-1)*M2*S4 + 8*M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(82)*N**2*CF * ( 16*S*TG*M2 - 16*S*M2*S4 + 
     +    16*S*M2**2 + 8*S**2*M2 - 16*TG*M2*S4 + 16*TG*M2**2 + 8*TG**2*
     +    M2 + 8*M2*S4**2 - 16*M2**2*S4 + 8*M2**3 )
     +
      M2QQH2 = M2QQH2 + ANG2(83)*N**2*CF * (  - 24*S*TG*M2*SYMBX*S4 + 
     +    32*S*TG*M2**2*SYMBX + 12*S*TG**2*M2*SYMBX + 12*S*M2*SYMBX*
     +    S4**2 - 12*S*M2 - 32*S*M2**2*SYMBX*S4 + 20*S*M2**3*SYMBX + 12
     +    *S**2*TG*M2*SYMBX - 12*S**2*M2*SYMBX*S4 + 16*S**2*M2**2*SYMBX
     +     + 4*S**3*M2*SYMBX + 12*TG*M2*SYMBX*S4**2 - 12*TG*M2 - 32*TG*
     +    M2**2*SYMBX*S4 + 20*TG*M2**3*SYMBX - 12*TG**2*M2*SYMBX*S4 + 
     +    16*TG**2*M2**2*SYMBX + 4*TG**3*M2*SYMBX + 4*U1**(-1)*M2**3 + 
     +    8*M2*MG2 - 4*M2*SYMBX*S4**3 + 12*M2*S4 + 16*M2**2*SYMBX*S4**2
     +     - 16*M2**2 - 20*M2**3*SYMBX*S4 + 8*M2**4*SYMBX )
     +
      M2QQH2 = M2QQH2 + ANG2(84)*N**2*CF * (  - 8*S**2*T1**(-1)*M2 + 8*
     +    TG*M2 + 8*T1**(-1)*M2*S4**2 - 16*M2*S4 + 8*M2**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N*CF**2 * (  - 136*S*TG*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 - 24*S*TG*M2*SYMBX + 144*S*TG*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 24*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 - 16*
     +    S*TG*(TG+UG)**(-1) + 40*S*TG*SYMBX*S4 + 88*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX - 40*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 - 16*
     +    S*TG**2*SYMBX + 16*S*TG**3*(TG+UG)**(-1)*SYMBX + 16*S*
     +    T1**(-1)*M2*(TG+UG)**(-1)*S4 + 32*S*T1**(-1)*M2*SYMBX*S4**2
     +     + 16*S*T1**(-1)*M2 - 16*S*T1**(-1)*M2**2*SYMBX*S4 + 40*S*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**2 - 32*S*M2*(TG+UG)**(-1) - 112*S*
     +    M2**2*(TG+UG)**(-1)*SYMBX*S4 + 8*S*M2**2*SYMBX + 72*S*M2**3*
     +    (TG+UG)**(-1)*SYMBX + 16*S*(TG+UG)**(-1)*S4 - 24*S*SYMBX*
     +    S4**2 + 112*S**2*TG*M2*(TG+UG)**(-1)*SYMBX - 32*S**2*TG*
     +    (TG+UG)**(-1)*SYMBX*S4 - 32*S**2*TG*SYMBX + 32*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBX - 16*S**2*T1**(-1)*M2*(TG+UG)**(-1) - 16*
     +    S**2*T1**(-1)*M2*SYMBX*S4 - 64*S**2*M2*(TG+UG)**(-1)*SYMBX*S4
     +     - 32*S**2*M2*SYMBX )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N*CF**2 * ( 96*S**2*M2**2*
     +    (TG+UG)**(-1)*SYMBX - 16*S**2*(TG+UG)**(-1) + 32*S**2*SYMBX*
     +    S4 + 16*S**3*TG*(TG+UG)**(-1)*SYMBX + 32*S**3*M2*
     +    (TG+UG)**(-1)*SYMBX - 16*S**3*SYMBX + 48*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 - 8*TG*M2*(TG+UG)**(-1) - 8*TG*M2*SYMBX*S4 - 80*
     +    TG*M2**2*(TG+UG)**(-1)*SYMBX*S4 + 24*TG*M2**2*SYMBX + 40*TG*
     +    M2**3*(TG+UG)**(-1)*SYMBX - 8*TG*(TG+UG)**(-1)*SYMBX*S4**3 + 
     +    8*TG*(TG+UG)**(-1)*S4 - 16*TG*SYMBX*S4**2 - 48*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX*S4 + 8*TG**2*M2*SYMBX + 32*TG**2*M2**2*
     +    (TG+UG)**(-1)*SYMBX + 16*TG**2*(TG+UG)**(-1)*SYMBX*S4**2 + 8*
     +    TG**2*SYMBX*S4 + 8*TG**3*M2*(TG+UG)**(-1)*SYMBX - 8*TG**3*
     +    (TG+UG)**(-1)*SYMBX*S4 - 16*T1**(-1)*U1**(-1)*M2**3 - 16*
     +    T1**(-1)*M2*SYMBX*S4**3 - 16*T1**(-1)*M2*S4 + 16*T1**(-1)*
     +    M2**2*SYMBX*S4**2 + 16*T1**(-1)*M2**2 - 8*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**3 + 16*M2*(TG+UG)**(-1)*S4 + 16*M2*SYMBX*S4**2 + 16
     +    *M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N*CF**2 * ( 32*M2**2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 - 8*M2**2*(TG+UG)**(-1) - 40*M2**2*SYMBX*S4 - 40*
     +    M2**3*(TG+UG)**(-1)*SYMBX*S4 + 16*M2**3*SYMBX + 16*M2**4*
     +    (TG+UG)**(-1)*SYMBX - 8*(TG+UG)**(-1)*S4**2 + 8*SYMBX*S4**3 )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N**2*CF * ( 68*S*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4 + 26*S*TG*M2*SYMBX - 72*S*TG*M2**2*(TG+UG)**(-1)*
     +    SYMBX - 12*S*TG*(TG+UG)**(-1)*SYMBX*S4**2 + 8*S*TG*
     +    (TG+UG)**(-1) - 28*S*TG*SYMBX*S4 - 44*S*TG**2*M2*
     +    (TG+UG)**(-1)*SYMBX + 20*S*TG**2*(TG+UG)**(-1)*SYMBX*S4 + 12*
     +    S*TG**2*SYMBX - 8*S*TG**3*(TG+UG)**(-1)*SYMBX + 8*S*T1**(-2)*
     +    M2*S4 - 8*S*T1**(-1)*M2*(TG+UG)**(-1)*S4 - 8*S*T1**(-1)*M2*
     +    SYMBX*S4**2 - 16*S*T1**(-1)*M2 + 4*S*T1**(-1)*M2**2*SYMBX*S4
     +     - 4*S*T1**(-1)*S4 - 20*S*M2*(TG+UG)**(-1)*SYMBX*S4**2 + 16*S
     +    *M2*(TG+UG)**(-1) - 22*S*M2*SYMBX*S4 + 56*S*M2**2*
     +    (TG+UG)**(-1)*SYMBX*S4 + 12*S*M2**2*SYMBX - 36*S*M2**3*
     +    (TG+UG)**(-1)*SYMBX - 8*S*(TG+UG)**(-1)*S4 + 16*S*SYMBX*S4**2
     +     + 2*S - 56*S**2*TG*M2*(TG+UG)**(-1)*SYMBX + 16*S**2*TG*
     +    (TG+UG)**(-1)*SYMBX*S4 + 18*S**2*TG*SYMBX - 16*S**2*TG**2*
     +    (TG+UG)**(-1)*SYMBX - 8*S**2*T1**(-2)*M2 + 8*S**2*T1**(-1)*M2
     +    *(TG+UG)**(-1) )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N**2*CF * ( 4*S**2*T1**(-1)*M2*SYMBX*
     +    S4 + 32*S**2*M2*(TG+UG)**(-1)*SYMBX*S4 + 20*S**2*M2*SYMBX - 
     +    48*S**2*M2**2*(TG+UG)**(-1)*SYMBX + 8*S**2*(TG+UG)**(-1) - 18
     +    *S**2*SYMBX*S4 - 8*S**3*TG*(TG+UG)**(-1)*SYMBX - 16*S**3*M2*
     +    (TG+UG)**(-1)*SYMBX + 8*S**3*SYMBX - 24*TG*M2*(TG+UG)**(-1)*
     +    SYMBX*S4**2 + 4*TG*M2*(TG+UG)**(-1) - 20*TG*M2*SYMBX*S4 + 40*
     +    TG*M2**2*(TG+UG)**(-1)*SYMBX*S4 + 8*TG*M2**2*SYMBX - 20*TG*
     +    M2**3*(TG+UG)**(-1)*SYMBX + 4*TG*(TG+UG)**(-1)*SYMBX*S4**3 - 
     +    4*TG*(TG+UG)**(-1)*S4 + 14*TG*SYMBX*S4**2 - 2*TG + 24*TG**2*
     +    M2*(TG+UG)**(-1)*SYMBX*S4 + 6*TG**2*M2*SYMBX - 16*TG**2*M2**2
     +    *(TG+UG)**(-1)*SYMBX - 8*TG**2*(TG+UG)**(-1)*SYMBX*S4**2 - 10
     +    *TG**2*SYMBX*S4 - 4*TG**3*M2*(TG+UG)**(-1)*SYMBX + 4*TG**3*
     +    (TG+UG)**(-1)*SYMBX*S4 + 2*TG**3*SYMBX - 4*T1**(-1)*U1**(-1)*
     +    M2**2*S4 + 4*T1**(-1)*U1**(-1)*M2**3 + 8*T1**(-1)*M2*MG2 + 4*
     +    T1**(-1)*M2*SYMBX*S4**3 + 20*T1**(-1)*M2*S4 - 4*T1**(-1)*
     +    M2**2*SYMBX*S4**2 )
     +
      M2QQH2 = M2QQH2 + ANG2(85)*N**2*CF * (  - 12*T1**(-1)*M2**2 - 4*
     +    T1**(-1)*S4**2 - 2*U1**(-1)*M2*S4 + 6*U1**(-1)*M2**2 + 4*M2*
     +    (TG+UG)**(-1)*SYMBX*S4**3 - 8*M2*(TG+UG)**(-1)*S4 + 10*M2*
     +    SYMBX*S4**2 - 16*M2 - 16*M2**2*(TG+UG)**(-1)*SYMBX*S4**2 + 4*
     +    M2**2*(TG+UG)**(-1) - 8*M2**2*SYMBX*S4 + 20*M2**3*
     +    (TG+UG)**(-1)*SYMBX*S4 + 4*M2**3*SYMBX - 8*M2**4*
     +    (TG+UG)**(-1)*SYMBX + 4*(TG+UG)**(-1)*S4**2 - 6*SYMBX*S4**3
     +     + 4*S4 )
     +
      M2QQH2 = M2QQH2 + ANG2(86)*N*CF**2 * ( 8*S**2*M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(87)*N*CF**2 * ( 16*S*T1**(-1)*M2*S4 - 16*S
     +    *M2 )
     +
      M2QQH2 = M2QQH2 + ANG2(87)*N**2*CF * (  - 8*S**2*T1**(-1)*M2 )
     +
      M2QQH2 = M2QQH2 + COLO2(9)*N*CF**2*(S4+MS2) * ( 8*S**(-1)*TG*
     +    (S+UG)**(-1)*S4**(-1) - 32*S**(-1)*T1**(-4)*M2**2*S4 + 32*
     +    S**(-1)*T1**(-3)*M2*S4 + 32*S**(-1)*T1**(-3)*M2**2 + 16*
     +    S**(-1)*T1**(-2)*M2*MG2*(S+UG)**(-1) - 16*S**(-1)*T1**(-2)*M2
     +    *MG2*S4**(-1) - 32*S**(-1)*T1**(-2)*M2 - 16*S**(-1)*T1**(-2)*
     +    M2**2*(S+UG)**(-1) - 16*S**(-1)*T1**(-2)*S4 - 16*S**(-1)*
     +    T1**(-1)*M2*MG2*(S+UG)**(-1)*S4**(-1) + 16*S**(-1)*T1**(-1)*
     +    M2*S4**(-1) + 16*S**(-1)*T1**(-1)*M2**2*(S+UG)**(-1)*S4**(-1)
     +     + 16*S**(-1)*T1**(-1) + 8*S**(-1)*M2*(S+UG)**(-1)*S4**(-1)
     +     - 8*S**(-1)*(S+UG)**(-1) - 32*T1**(-4)*M2*MG2 + 32*T1**(-4)*
     +    M2**2 - 32*T1**(-3)*M2 - 16*T1**(-2)*M2*MG2*(S+UG)**(-1)*
     +    S4**(-1) + 16*T1**(-2)*M2**2*(S+UG)**(-1)*S4**(-1) + 16*
     +    T1**(-2) - 16*T1**(-1)*M2*(S+UG)**(-1)*S4**(-1) + 8*
     +    (S+UG)**(-1)*S4**(-1) )
     +
      M2QQH2 = M2QQH2 + COLO2(9)*N**2*CF*(S4+MS2) * (  - 8*S*UG**(-1)*
     +    (S+UG)**(-1)*S4**(-1) - 32*UG**(-2)*T1**(-4)*M2*MG2*S4**2 + 
     +    64*UG**(-2)*T1**(-3)*M2*MG2*S4 - 48*UG**(-2)*T1**(-2)*M2*MG2
     +     + 16*UG**(-2)*T1**(-1)*M2*MG2*S4**(-1) + 32*UG**(-1)*
     +    T1**(-4)*M2*MG2*S4 - 32*UG**(-1)*T1**(-4)*M2**2*S4 - 32*
     +    UG**(-1)*T1**(-3)*M2*MG2 + 32*UG**(-1)*T1**(-3)*M2*S4 + 32*
     +    UG**(-1)*T1**(-3)*M2**2 + 16*UG**(-1)*T1**(-2)*M2*MG2*
     +    S4**(-1) - 32*UG**(-1)*T1**(-2)*M2 - 16*UG**(-1)*T1**(-2)*
     +    M2**2*S4**(-1) - 16*UG**(-1)*T1**(-2)*S4 + 16*UG**(-1)*
     +    T1**(-1)*M2*S4**(-1) + 16*UG**(-1)*T1**(-1) - 8*(S+UG)**(-1)*
     +    S4**(-1) )

      IF (IFL.EQ.1) THEN

      M2QQH3 = 0D0
      M2QQH3 = M2QQH3 + N*CF*(S4+MS2) * ( 4*S**(-1)*TG*(S+T1)**(-1)*
     +    S4**(-1) - 4*S**(-1)*TG*(S+U1)**(-1)*S4**(-1) + 4*S**(-1)*
     +    T1**(-1)*M2*S4**(-1) + 4*S**(-1)*T1**(-1) + 4*S**(-1)*
     +    U1**(-1)*M2*S4**(-1) + 4*S**(-1)*U1**(-1) - 4*S**(-1)*M2*
     +    (S+U1)**(-1)*S4**(-1) - 4*S**(-1)*(S+T1)**(-1) - 8*S**(-1)*
     +    S4**(-1) - 8*S*T1**(-1)*(S+TG)**(-1)*S4**(-1) - 8*S*U1**(-1)*
     +    (S+UG)**(-1)*S4**(-1) + 4*S*(S+T1)**(-1)*(S+TG)**(-1)*
     +    S4**(-1) + 12*TG**(-1)*T1**(-1)*U1**(-1)*M2*(S+U1)**(-1)*S4
     +     - 8*TG**(-1)*T1**(-1)*U1**(-1)*M2 + 12*TG**(-1)*T1**(-1)*
     +    U1**(-1)*M2**2*(S+U1)**(-1) - 4*TG**(-1)*T1**(-1)*U1**(-1)*
     +    M2**2*S4**(-1) + 4*TG**(-1)*T1**(-1)*U1**(-1)*M2**3*
     +    (S+U1)**(-1)*S4**(-1) + 4*TG**(-1)*T1**(-1)*U1**(-1)*
     +    (S+U1)**(-1)*S4**2 - 4*TG**(-1)*T1**(-1)*U1**(-1)*S4 - 8*
     +    TG**(-1)*T1**(-1)*M2*(S+U1)**(-1) + 4*TG**(-1)*T1**(-1)*M2*
     +    S4**(-1) - 4*TG**(-1)*T1**(-1)*M2**2*(S+U1)**(-1)*S4**(-1) - 
     +    4*TG**(-1)*T1**(-1)*(S+U1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + N*CF*(S4+MS2) * ( 4*TG**(-1)*T1**(-1) - 16*
     +    TG**(-1)*U1**(-1)*M2*(S+U1)**(-1) + 8*TG**(-1)*U1**(-1)*M2*
     +    S4**(-1) - 4*TG**(-1)*U1**(-1)*M2**2*(S+U1)**(-1)*S4**(-1) - 
     +    12*TG**(-1)*U1**(-1)*(S+U1)**(-1)*S4 + 8*TG**(-1)*U1**(-1) + 
     +    4*TG**(-1)*M2*(S+U1)**(-1)*S4**(-1) + 8*TG**(-1)*(S+U1)**(-1)
     +     - 4*TG**(-1)*S4**(-1) + 4*TG*(S+T1)**(-1)*(S+TG)**(-1)*
     +    S4**(-1) - 4*TG*(S+U1)**(-1)*(S+UG)**(-1)*S4**(-1) + 8*
     +    UG**(-1)*T1**(-1)*U1**(-1)*M2*(S+T1)**(-1)*S4 - 4*UG**(-1)*
     +    T1**(-1)*U1**(-1)*M2 + 4*UG**(-1)*T1**(-1)*U1**(-1)*M2**2*
     +    (S+T1)**(-1) + 4*UG**(-1)*T1**(-1)*U1**(-1)*(S+T1)**(-1)*
     +    S4**2 - 4*UG**(-1)*T1**(-1)*U1**(-1)*S4 - 8*UG**(-1)*T1**(-1)
     +    *M2*(S+T1)**(-1) + 4*UG**(-1)*T1**(-1)*M2*S4**(-1) - 8*
     +    UG**(-1)*T1**(-1)*(S+T1)**(-1)*S4 + 4*UG**(-1)*T1**(-1) - 4*
     +    UG**(-1)*U1**(-1)*M2*(S+T1)**(-1) - 4*UG**(-1)*U1**(-1)*
     +    (S+T1)**(-1)*S4 + 4*UG**(-1)*U1**(-1) + 4*UG**(-1)*
     +    (S+T1)**(-1) )
     +
      M2QQH3 = M2QQH3 + N*CF*(S4+MS2) * ( 8*T1**(-1)*U1**(-1)*M2*
     +    (S+U1)**(-1) - 4*T1**(-1)*U1**(-1)*M2*S4**(-1) + 4*T1**(-1)*
     +    U1**(-1)*M2**2*(S+U1)**(-1)*S4**(-1) + 4*T1**(-1)*U1**(-1)*
     +    (S+U1)**(-1)*S4 - 4*T1**(-1)*U1**(-1) + 8*T1**(-1)*M2*
     +    (S+TG)**(-1)*S4**(-1) + 4*T1**(-1)*M2*(S+U1)**(-1)*
     +    (S+UG)**(-1) - 4*T1**(-1)*M2*(S+U1)**(-1)*S4**(-1) - 4*
     +    T1**(-1)*(S+U1)**(-1)*(S+UG)**(-1)*S4 + 4*T1**(-1)*
     +    (S+U1)**(-1) + 8*T1**(-1)*S4**(-1) + 4*U1**(-1)*M2*
     +    (S+T1)**(-1)*(S+TG)**(-1) + 8*U1**(-1)*M2*(S+UG)**(-1)*
     +    S4**(-1) - 4*U1**(-1)*(S+T1)**(-1)*(S+TG)**(-1)*S4 + 8*
     +    U1**(-1)*(S+T1)**(-1) + 4*U1**(-1)*S4**(-1) - 4*M2*
     +    (S+T1)**(-1)*(S+TG)**(-1)*S4**(-1) - 8*M2*(S+U1)**(-1)*
     +    (S+UG)**(-1)*S4**(-1) + 4*(S+T1)**(-1)*(S+TG)**(-1) - 4*
     +    (S+T1)**(-1)*S4**(-1) - 8*(S+TG)**(-1)*S4**(-1) + 8*
     +    (S+U1)**(-1)*(S+UG)**(-1) - 8*(S+U1)**(-1)*S4**(-1) - 8*
     +    (S+UG)**(-1)*S4**(-1) )
     +
      M2QQH3 = M2QQH3 + CF**2*(S4+MS2) * ( 8*S*T1**(-1)*(S+TG)**(-1)*
     +    S4**(-1) + 8*S*U1**(-1)*(S+UG)**(-1)*S4**(-1) - 8*T1**(-1)*M2
     +    *(S+TG)**(-1)*S4**(-1) - 8*T1**(-1)*S4**(-1) - 8*U1**(-1)*M2*
     +    (S+UG)**(-1)*S4**(-1) - 8*U1**(-1)*S4**(-1) + 8*(S+TG)**(-1)*
     +    S4**(-1) + 8*(S+UG)**(-1)*S4**(-1) )
     +
      M2QQH3 = M2QQH3 + ANG2(16)*N*CF * (  - 4*T1**(-1)*U1**(-1)*MS2*S4
     +     - 4*T1**(-1)*U1**(-1)*S4**2 + 4*U1**(-1)*MS2 + 4*U1**(-1)*S4
     +     )
     +
      M2QQH3 = M2QQH3 + ANG2(16)*CF**2 * (  - 8*T1**(-2)*M2*MG2 + 8*
     +    T1**(-2)*M2*(S+U1)**(-1)*S4**2 - 8*T1**(-2)*M2*S4 + 8*
     +    T1**(-2)*M2**2 + 4*T1**(-1)*U1**(-1)*MS2*S4 + 4*T1**(-1)*
     +    U1**(-1)*S4**2 - 16*T1**(-1)*M2*(S+U1)**(-1)*S4 - 4*T1**(-1)*
     +    MS2*(S+U1)**(-1)*S4 + 8*T1**(-1)*MS2 - 4*T1**(-1)*
     +    (S+U1)**(-1)*S4**2 + 8*T1**(-1)*S4 - 4*U1**(-1)*MS2 - 4*
     +    U1**(-1)*S4 + 8*M2*(S+U1)**(-1) + 4*MS2*(S+U1)**(-1) + 4*
     +    (S+U1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(20)*N*CF * ( 2*S*T1**(-1)*M2*(S+U1)**(-1)*
     +    S4 + 2*S*T1**(-1)*MS2*(S+U1)**(-1)*S4 + 2*S*U1**(-1)*M2*
     +    (S+T1)**(-1)*S4 + 2*S*U1**(-1)*MS2*(S+T1)**(-1)*S4 - 4*
     +    T1**(-1)*U1**(-1)*M2*MG2*S4 - 4*T1**(-1)*U1**(-1)*M2*S4**2 + 
     +    4*T1**(-1)*U1**(-1)*M2**2*S4 - 4*T1**(-1)*U1**(-1)*MS2*S4**2
     +     + 4*T1**(-1)*U1**(-1)*S4**3 + 4*T1**(-1)*M2*(S+U1)**(-1)*
     +    S4**2 + 2*T1**(-1)*M2*S4 + 4*T1**(-1)*MS2*S4 - 2*T1**(-1)*
     +    S4**2 + 4*U1**(-1)*M2*(S+T1)**(-1)*S4**2 + 2*U1**(-1)*M2*S4
     +     + 4*U1**(-1)*MS2*S4 - 2*U1**(-1)*S4**2 - 2*M2*(S+T1)**(-1)*
     +    S4 - 2*M2*(S+U1)**(-1)*S4 - 2*(S+T1)**(-1)*S4**2 - 2*
     +    (S+U1)**(-1)*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(20)*CF**2 * ( 4*T1**(-1)*U1**(-1)*M2*MG2*
     +    S4 + 4*T1**(-1)*U1**(-1)*M2*S4**2 - 4*T1**(-1)*U1**(-1)*M2**2
     +    *S4 + 4*T1**(-1)*U1**(-1)*MS2*S4**2 - 4*T1**(-1)*U1**(-1)*
     +    S4**3 - 4*T1**(-1)*M2*S4 - 4*T1**(-1)*MS2*S4 - 4*U1**(-1)*M2*
     +    S4 - 4*U1**(-1)*MS2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(30)*N*CF * ( 2*S*T1**(-1)*M2*(S+U1)**(-1)*
     +    S4 - 2*S*T1**(-1)*M2 + 2*S*T1**(-1)*MS2*(S+U1)**(-1)*S4 - 2*S
     +    *T1**(-1)*MS2 - 4*S*U1**(-1)*M2*(S+T1)**(-1)*S4 - 4*S*
     +    U1**(-1)*M2 + 4*TG*M2*MG2*(S+T1)**(-1)*(S+2*M2)**(-1) - 2*TG*
     +    M2*MG2*(S+2*M2)**(-1)*SYMBY*S4 - 2*TG*M2*MG2*SYMBY - 2*TG*M2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1)*S4 - 4*TG*M2*(S+T1)**(-1) + 4*TG*
     +    M2**2*MG2*(S+2*M2)**(-1)*SYMBY - 2*TG*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1)*S4 - 2*TG*MS2*(S+T1)**(-1) - 2*TG*(S+T1)**(-1)
     +    *S4 - 2*TG**2*M2*MG2*(S+2*M2)**(-1)*SYMBY + 2*TG**2*M2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1) + 2*TG**2*M2*(S+2*M2)**(-1)*SYMBY
     +    *S4 + 2*TG**2*M2*SYMBY + 2*TG**2*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1) + 2*TG**2*MS2*(S+2*M2)**(-1)*SYMBY*S4 + 2*
     +    TG**2*MS2*SYMBY - 2*TG**3*M2*(S+2*M2)**(-1)*SYMBY - 2*TG**3*
     +    MS2*(S+2*M2)**(-1)*SYMBY - 2*T1**(-1)*M2*MG2 + 4*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4**2 - 4*T1**(-1)*M2*S4 + 4*T1**(-1)*M2**2 + 4*
     +    U1**(-1)*M2*(S+T1)**(-1)*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(30)*N*CF * ( 4*U1**(-1)*M2*S4 - 4*M2*
     +    (S+T1)**(-1)*S4 - 2*M2*(S+U1)**(-1)*S4 - 4*M2 + 2*MS2 - 2*
     +    (S+U1)**(-1)*S4**2 + 2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(36)*CF**2 * ( 8*S*M2 )
     +
      M2QQH3 = M2QQH3 + ANG2(38)*N*CF * ( 2*S*T1**(-1)*M2 + 4*S*
     +    T1**(-1)*MS2 + 2*S*T1**(-1)*S4 - 4*S*U1**(-1)*M2*(S+T1)**(-1)
     +    *S4 - 4*S*U1**(-1)*M2 - 4*S*M2*(S+T1)**(-1) - 4*S*MS2*
     +    (S+T1)**(-1) - 8*S**2*M2*(S+T1)**(-1)*(TG+UG)**(-1) - 8*S**2*
     +    MS2*(S+T1)**(-1)*(TG+UG)**(-1) - 2*TG*M2*(S+T1)**(-1) - 2*TG*
     +    (S+T1)**(-1)*S4 - 2*T1**(-1)*M2**2 - 2*T1**(-1)*S4**2 - 2*M2*
     +    (S+T1)**(-1)*S4 + 2*M2 + 2*(S+T1)**(-1)*S4**2 + 2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(38)*CF**2 * (  - 4*S*T1**(-1)*M2 - 8*S*
     +    T1**(-1)*MS2 - 4*S*T1**(-1)*S4 + 8*S*U1**(-1)*M2*(S+T1)**(-1)
     +    *S4 - 4*S*M2*(S+T1)**(-1) - 4*S*(S+T1)**(-1)*S4 + 8*S**2*M2*
     +    (S+T1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+T1)**(-1)*
     +    (TG+UG)**(-1) + 4*T1**(-1)*M2**2 + 4*T1**(-1)*S4**2 + 8*
     +    U1**(-1)*M2*S4 - 8*M2 )
     +
      M2QQH3 = M2QQH3 + ANG2(39)*CF**2 * ( 8*S*M2 )
     +
      M2QQH3 = M2QQH3 + ANG2(40)*N*CF * (  - 4*S*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4 - 4*S*T1**(-1)*M2 + 2*S*U1**(-1)*M2 + 4*S*
     +    U1**(-1)*MS2 + 2*S*U1**(-1)*S4 - 2*S*M2*(S+U1)**(-1) - 4*S*
     +    MS2*(S+U1)**(-1) + 2*S*(S+U1)**(-1)*S4 - 8*S**2*M2*
     +    (S+U1)**(-1)*(TG+UG)**(-1) - 8*S**2*MS2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) + 2*TG*M2*(S+U1)**(-1) + 2*TG*(S+U1)**(-1)*S4
     +     - 2*U1**(-1)*M2**2 - 2*U1**(-1)*S4**2 - 2*M2*(S+U1)**(-1)*S4
     +     + 2*M2 + 2*M2**2*(S+U1)**(-1) + 2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(40)*CF**2 * ( 8*S*T1**(-1)*M2*(S+U1)**(-1)
     +    *S4 - 4*S*U1**(-1)*M2 - 8*S*U1**(-1)*MS2 - 4*S*U1**(-1)*S4 - 
     +    4*S*M2*(S+U1)**(-1) - 4*S*(S+U1)**(-1)*S4 + 8*S**2*M2*
     +    (S+U1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) + 8*T1**(-1)*M2*S4 + 4*U1**(-1)*M2**2 + 4*
     +    U1**(-1)*S4**2 - 8*M2 )
     +
      M2QQH3 = M2QQH3 + ANG2(42)*N*CF * (  - 4*S*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4 + 2*S*T1**(-1)*M2 + 4*S*T1**(-1)*MS2 + 2*S*
     +    T1**(-1)*S4 + 4*S*U1**(-1)*M2 - 2*S*M2*(S+U1)**(-1) - 4*S*MS2
     +    *(S+U1)**(-1) + 2*S*(S+U1)**(-1)*S4 - 8*S**2*M2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) - 8*S**2*MS2*(S+U1)**(-1)*(TG+UG)**(-1) + 2*TG*
     +    M2*(S+U1)**(-1) + 2*TG*(S+U1)**(-1)*S4 + 4*T1**(-1)*M2*S4 - 2
     +    *T1**(-1)*M2**2 - 2*T1**(-1)*S4**2 - 4*U1**(-1)*M2*S4 - 2*M2*
     +    (S+U1)**(-1)*S4 + 2*M2 + 2*M2**2*(S+U1)**(-1) + 2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(42)*CF**2 * ( 8*S*T1**(-1)*M2*(S+U1)**(-1)
     +    *S4 - 4*S*T1**(-1)*M2 - 8*S*T1**(-1)*MS2 - 4*S*T1**(-1)*S4 - 
     +    4*S*M2*(S+U1)**(-1) - 4*S*(S+U1)**(-1)*S4 + 8*S**2*M2*
     +    (S+U1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) - 8*T1**(-1)*M2*S4 + 4*T1**(-1)*M2**2 + 4*
     +    T1**(-1)*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(58)*N*CF * ( 2*S*U1**(-1)*M2*(S+T1)**(-1)*
     +    S4 + 2*S*U1**(-1)*MS2*(S+T1)**(-1)*S4 - 6*S*M2*(S+T1)**(-1)
     +     - 6*S*MS2*(S+T1)**(-1) - 8*S**2*M2*(S+T1)**(-1)*
     +    (TG+UG)**(-1) - 8*S**2*MS2*(S+T1)**(-1)*(TG+UG)**(-1) - 4*TG*
     +    M2*MG2*(S+T1)**(-1)*(S+2*M2)**(-1) + 2*TG*M2*MG2*
     +    (S+2*M2)**(-1)*SYMBY*S4 + 2*TG*M2*MG2*SYMBY + 2*TG*M2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1)*S4 - 4*TG*M2**2*MG2*
     +    (S+2*M2)**(-1)*SYMBY + 2*TG*MS2*(S+T1)**(-1)*(S+2*M2)**(-1)*
     +    S4 + 2*TG**2*M2*MG2*(S+2*M2)**(-1)*SYMBY - 2*TG**2*M2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1) - 2*TG**2*M2*(S+2*M2)**(-1)*SYMBY
     +    *S4 - 2*TG**2*M2*SYMBY - 2*TG**2*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1) - 2*TG**2*MS2*(S+2*M2)**(-1)*SYMBY*S4 - 2*
     +    TG**2*MS2*SYMBY + 2*TG**3*M2*(S+2*M2)**(-1)*SYMBY + 2*TG**3*
     +    MS2*(S+2*M2)**(-1)*SYMBY - 2*U1**(-1)*M2*MG2*(S+T1)**(-1)*S4
     +     + 2*U1**(-1)*M2*(S+T1)**(-1)*S4**2 + 2*U1**(-1)*MS2*
     +    (S+T1)**(-1)*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(58)*N*CF * ( 2*M2*MG2*(S+T1)**(-1) - 2*M2*
     +    (S+T1)**(-1)*S4 - 2*MS2*(S+T1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(58)*CF**2 * (  - 4*S*U1**(-1)*M2*
     +    (S+T1)**(-1)*S4 - 4*S*U1**(-1)*MS2*(S+T1)**(-1)*S4 + 8*S*M2*
     +    (S+T1)**(-1) + 8*S*MS2*(S+T1)**(-1) + 8*S**2*M2*(S+T1)**(-1)*
     +    (TG+UG)**(-1) + 8*S**2*MS2*(S+T1)**(-1)*(TG+UG)**(-1) + 4*
     +    U1**(-1)*M2*MG2*(S+T1)**(-1)*S4 - 4*U1**(-1)*M2*(S+T1)**(-1)*
     +    S4**2 - 4*U1**(-1)*MS2*(S+T1)**(-1)*S4**2 - 4*M2*MG2*
     +    (S+T1)**(-1) + 4*M2*(S+T1)**(-1)*S4 + 4*MS2*(S+T1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(59)*N*CF * (  - 10*S*TG*M2*SYMBX - 10*S*TG
     +    *MS2*SYMBX + 2*S*T1**(-1)*M2*(S+U1)**(-1)*S4 + 2*S*T1**(-1)*
     +    MS2*(S+U1)**(-1)*S4 - 6*S*M2*MG2*SYMBX - 8*S*M2*(S+U1)**(-1)
     +     + 8*S*M2*SYMBX*S4 - 8*S*MS2*(S+U1)**(-1) + 8*S*MS2*SYMBX*S4
     +     - 8*S**2*M2*(S+U1)**(-1)*(TG+UG)**(-1) - 4*S**2*M2*SYMBX - 8
     +    *S**2*MS2*(S+U1)**(-1)*(TG+UG)**(-1) - 4*S**2*MS2*SYMBX + 8*
     +    TG*M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1) - 14*TG*M2*MG2*
     +    (S+2*M2)**(-1)*SYMBX*S4 - 2*TG*M2*MG2*SYMBX + 2*TG*M2*
     +    (S+U1)**(-1)*(S+2*M2)**(-1)*S4 - 4*TG*M2*(S+U1)**(-1) - 2*TG*
     +    M2*(S+2*M2)**(-1)*SYMBX*S4**2 + 12*TG*M2*SYMBX*S4 - 6*TG*
     +    M2**2*MG2*(S+2*M2)**(-1)*SYMBX + 2*TG*MS2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1)*S4 - 4*TG*MS2*(S+U1)**(-1) - 2*TG*MS2*
     +    (S+2*M2)**(-1)*SYMBX*S4**2 + 12*TG*MS2*SYMBX*S4 + 8*TG**2*M2*
     +    MG2*(S+2*M2)**(-1)*SYMBX - 2*TG**2*M2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) + 4*TG**2*M2*(S+2*M2)**(-1)*SYMBX*S4 - 8*TG**2
     +    *M2*SYMBX )
     +
      M2QQH3 = M2QQH3 + ANG2(59)*N*CF * (  - 2*TG**2*MS2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) + 4*TG**2*MS2*(S+2*M2)**(-1)*SYMBX*S4 - 8*
     +    TG**2*MS2*SYMBX - 2*TG**3*M2*(S+2*M2)**(-1)*SYMBX - 2*TG**3*
     +    MS2*(S+2*M2)**(-1)*SYMBX - 2*T1**(-1)*M2*MG2*(S+U1)**(-1)*S4
     +     + 2*T1**(-1)*M2*(S+U1)**(-1)*S4**2 + 2*T1**(-1)*MS2*
     +    (S+U1)**(-1)*S4**2 - 6*M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4
     +     + 6*M2*MG2*(S+U1)**(-1) + 6*M2*MG2*(S+2*M2)**(-1)*SYMBX*
     +    S4**2 - 4*M2*SYMBX*S4**2 - 6*M2**2*MG2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) + 6*M2**2*MG2*(S+2*M2)**(-1)*SYMBX*S4 - 2*
     +    M2**2*MG2*SYMBX - 4*MS2*SYMBX*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(59)*CF**2 * (  - 4*S*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4 - 4*S*T1**(-1)*MS2*(S+U1)**(-1)*S4 + 8*S*M2*
     +    (S+U1)**(-1) + 8*S*MS2*(S+U1)**(-1) + 8*S**2*M2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) + 8*S**2*MS2*(S+U1)**(-1)*(TG+UG)**(-1) + 4*
     +    T1**(-1)*M2*MG2*(S+U1)**(-1)*S4 - 4*T1**(-1)*M2*(S+U1)**(-1)*
     +    S4**2 - 4*T1**(-1)*MS2*(S+U1)**(-1)*S4**2 - 4*M2*MG2*
     +    (S+U1)**(-1) + 4*M2*(S+U1)**(-1)*S4 + 4*MS2*(S+U1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(73)*N*CF * ( 2*S*T1**(-1)*M2*(S+U1)**(-1)*
     +    S4 - 2*S*T1**(-1)*M2 + 2*S*T1**(-1)*MS2*(S+U1)**(-1)*S4 - 2*S
     +    *T1**(-1)*MS2 + 2*S*U1**(-1)*M2 + 2*S*U1**(-1)*MS2 - 8*S*M2*
     +    (S+U1)**(-1) - 8*S*MS2*(S+U1)**(-1) - 8*S**2*M2*(S+U1)**(-1)*
     +    (TG+UG)**(-1) - 8*S**2*MS2*(S+U1)**(-1)*(TG+UG)**(-1) + 8*TG*
     +    M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1) - 2*TG*M2*MG2*
     +    (S+2*M2)**(-1)*SYMBY*S4 - 2*TG*M2*MG2*SYMBY + 2*TG*M2*
     +    (S+U1)**(-1)*(S+2*M2)**(-1)*S4 - 6*TG*M2*(S+U1)**(-1) + 4*TG*
     +    M2**2*MG2*(S+2*M2)**(-1)*SYMBY + 2*TG*MS2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1)*S4 - 6*TG*MS2*(S+U1)**(-1) - 2*TG**2*M2*MG2*
     +    (S+2*M2)**(-1)*SYMBY - 2*TG**2*M2*(S+U1)**(-1)*(S+2*M2)**(-1)
     +     + 2*TG**2*M2*(S+2*M2)**(-1)*SYMBY*S4 + 2*TG**2*M2*SYMBY - 2*
     +    TG**2*MS2*(S+U1)**(-1)*(S+2*M2)**(-1) + 2*TG**2*MS2*
     +    (S+2*M2)**(-1)*SYMBY*S4 + 2*TG**2*MS2*SYMBY - 2*TG**3*M2*
     +    (S+2*M2)**(-1)*SYMBY - 2*TG**3*MS2*(S+2*M2)**(-1)*SYMBY - 4*
     +    T1**(-1)*M2*MG2*(S+U1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(73)*N*CF * (  - 2*T1**(-1)*M2*MG2*
     +    (S+2*M2)**(-1)*S4 + 6*T1**(-1)*M2**2*MG2*(S+2*M2)**(-1) - 2*
     +    U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 - 4*U1**(-1)*M2*MG2 + 2*
     +    U1**(-1)*M2*S4 + 6*U1**(-1)*M2**2*MG2*(S+2*M2)**(-1) + 2*
     +    U1**(-1)*MS2*S4 - 6*M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4 + 6
     +    *M2*MG2*(S+U1)**(-1) + 14*M2*MG2*(S+2*M2)**(-1) + 4*M2*
     +    (S+U1)**(-1)*S4 - 2*M2*(S+2*M2)**(-1)*S4 - 8*M2 - 6*M2**2*MG2
     +    *(S+U1)**(-1)*(S+2*M2)**(-1) + 4*MS2*(S+U1)**(-1)*S4 - 2*MS2*
     +    (S+2*M2)**(-1)*S4 - 8*MS2 )
     +
      M2QQH3 = M2QQH3 + ANG2(73)*CF**2 * (  - 4*S*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4 + 4*S*T1**(-1)*M2 - 4*S*T1**(-1)*MS2*
     +    (S+U1)**(-1)*S4 + 4*S*T1**(-1)*MS2 - 4*S*U1**(-1)*M2 - 4*S*
     +    U1**(-1)*MS2 + 8*S*M2*(S+U1)**(-1) + 8*S*MS2*(S+U1)**(-1) + 8
     +    *S**2*M2*(S+U1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+U1)**(-1)
     +    *(TG+UG)**(-1) + 4*TG*M2*(S+U1)**(-1) + 4*TG*MS2*(S+U1)**(-1)
     +     + 8*T1**(-1)*M2*MG2*(S+U1)**(-1)*S4 + 4*T1**(-1)*M2*MG2*
     +    (S+2*M2)**(-1)*S4 - 12*T1**(-1)*M2**2*MG2*(S+2*M2)**(-1) + 4*
     +    U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 + 8*U1**(-1)*M2*MG2 - 4*
     +    U1**(-1)*M2*S4 - 12*U1**(-1)*M2**2*MG2*(S+2*M2)**(-1) - 4*
     +    U1**(-1)*MS2*S4 - 4*M2*MG2*(S+U1)**(-1) - 8*M2*MG2*
     +    (S+2*M2)**(-1) - 4*M2*(S+U1)**(-1)*S4 + 4*M2 + 16*M2**2*MG2*
     +    (S+2*M2)**(-2) - 4*MS2*(S+U1)**(-1)*S4 + 4*MS2 )
     +
      M2QQH3 = M2QQH3 + ANG2(75)*N*CF * ( 4*S*T1**(-1)*M2 - 4*S*
     +    U1**(-1)*M2*(S+T1)**(-1)*S4 + 2*S*U1**(-1)*M2 + 4*S*U1**(-1)*
     +    MS2 + 2*S*U1**(-1)*S4 - 4*S*M2*(S+T1)**(-1) - 4*S*MS2*
     +    (S+T1)**(-1) - 8*S**2*M2*(S+T1)**(-1)*(TG+UG)**(-1) - 8*S**2*
     +    MS2*(S+T1)**(-1)*(TG+UG)**(-1) - 2*TG*M2*(S+T1)**(-1) - 2*TG*
     +    (S+T1)**(-1)*S4 - 4*T1**(-1)*M2*S4 + 4*U1**(-1)*M2*S4 - 2*
     +    U1**(-1)*M2**2 - 2*U1**(-1)*S4**2 - 2*M2*(S+T1)**(-1)*S4 + 2*
     +    M2 + 2*(S+T1)**(-1)*S4**2 + 2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(75)*CF**2 * ( 8*S*U1**(-1)*M2*(S+T1)**(-1)
     +    *S4 - 4*S*U1**(-1)*M2 - 8*S*U1**(-1)*MS2 - 4*S*U1**(-1)*S4 - 
     +    4*S*M2*(S+T1)**(-1) - 4*S*(S+T1)**(-1)*S4 + 8*S**2*M2*
     +    (S+T1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+T1)**(-1)*
     +    (TG+UG)**(-1) - 8*U1**(-1)*M2*S4 + 4*U1**(-1)*M2**2 + 4*
     +    U1**(-1)*S4**2 )
     +
      M2QQH3 = M2QQH3 + ANG2(76)*N*CF * (  - 4*S*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4 + 2*S*T1**(-1)*M2 + 4*S*T1**(-1)*MS2 + 2*S*
     +    T1**(-1)*S4 - 4*S*U1**(-1)*M2*(S+T1)**(-1)*S4 + 2*S*U1**(-1)*
     +    M2 + 4*S*U1**(-1)*MS2 + 2*S*U1**(-1)*S4 + 6*S*M2*(S+U1)**(-1)
     +     + 4*S*MS2*(S+U1)**(-1) + 2*S*(S+U1)**(-1)*S4 + 4*TG*M2*MG2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1) - 8*TG*M2*MG2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) - 2*TG*M2*(S+T1)**(-1)*(S+2*M2)**(-1)*S4 - 4*
     +    TG*M2*(S+T1)**(-1) - 2*TG*M2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4
     +     + 8*TG*M2*(S+U1)**(-1) - 2*TG*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1)*S4 - 2*TG*MS2*(S+T1)**(-1) - 2*TG*MS2*
     +    (S+U1)**(-1)*(S+2*M2)**(-1)*S4 + 6*TG*MS2*(S+U1)**(-1) - 2*TG
     +    *(S+T1)**(-1)*S4 + 2*TG*(S+U1)**(-1)*S4 + 2*TG**2*M2*
     +    (S+T1)**(-1)*(S+2*M2)**(-1) + 2*TG**2*M2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) + 2*TG**2*MS2*(S+T1)**(-1)*(S+2*M2)**(-1) + 2*
     +    TG**2*MS2*(S+U1)**(-1)*(S+2*M2)**(-1) - 12*T1**(-1)*U1**(-1)*
     +    M2*MG2*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(76)*N*CF * (  - 12*T1**(-1)*U1**(-1)*M2*
     +    S4**2 + 12*T1**(-1)*U1**(-1)*M2**2*S4 - 12*T1**(-1)*U1**(-1)*
     +    MS2*S4**2 - 4*T1**(-1)*U1**(-1)*S4**3 + 2*T1**(-1)*M2*MG2*
     +    (S+2*M2)**(-1)*S4 + 2*T1**(-1)*M2*MG2 + 4*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4**2 - 6*T1**(-1)*M2**2*MG2*(S+2*M2)**(-1) - 2*
     +    T1**(-1)*M2**2 + 10*T1**(-1)*MS2*S4 + 4*T1**(-1)*S4**2 + 2*
     +    U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 + 2*U1**(-1)*M2*MG2 + 4*
     +    U1**(-1)*M2*(S+T1)**(-1)*S4**2 - 6*U1**(-1)*M2**2*MG2*
     +    (S+2*M2)**(-1) - 2*U1**(-1)*M2**2 + 10*U1**(-1)*MS2*S4 + 4*
     +    U1**(-1)*S4**2 + 6*M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4 - 2*
     +    M2*MG2*(S+U1)**(-1) - 14*M2*MG2*(S+2*M2)**(-1) - 4*M2*
     +    (S+T1)**(-1)*S4 - 8*M2*(S+U1)**(-1)*S4 + 2*M2*(S+2*M2)**(-1)*
     +    S4 + 14*M2 + 6*M2**2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1) + 2*
     +    M2**2*(S+U1)**(-1) - 4*MS2*(S+U1)**(-1)*S4 + 2*MS2*
     +    (S+2*M2)**(-1)*S4 + 10*MS2 - 2*(S+U1)**(-1)*S4**2 + 4*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(76)*CF**2 * ( 12*T1**(-1)*U1**(-1)*M2*MG2*
     +    S4 + 12*T1**(-1)*U1**(-1)*M2*S4**2 - 12*T1**(-1)*U1**(-1)*
     +    M2**2*S4 + 12*T1**(-1)*U1**(-1)*MS2*S4**2 + 4*T1**(-1)*
     +    U1**(-1)*S4**3 - 4*T1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 - 8*
     +    T1**(-1)*M2*MG2 + 12*T1**(-1)*M2**2*MG2*(S+2*M2)**(-1) + 4*
     +    T1**(-1)*M2**2 - 8*T1**(-1)*MS2*S4 - 4*T1**(-1)*S4**2 - 4*
     +    U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 - 8*U1**(-1)*M2*MG2 + 12*
     +    U1**(-1)*M2**2*MG2*(S+2*M2)**(-1) + 4*U1**(-1)*M2**2 - 8*
     +    U1**(-1)*MS2*S4 - 4*U1**(-1)*S4**2 + 8*M2*MG2*(S+2*M2)**(-1)
     +     - 8*M2 - 16*M2**2*MG2*(S+2*M2)**(-2) )
     +
      M2QQH3 = M2QQH3 + ANG2(78)*N*CF * ( 10*S*TG*M2*SYMBX + 10*S*TG*
     +    MS2*SYMBX + 2*S*T1**(-1)*M2 + 2*S*T1**(-1)*MS2 + 2*S*U1**(-1)
     +    *M2*(S+T1)**(-1)*S4 - 2*S*U1**(-1)*M2 + 2*S*U1**(-1)*MS2*
     +    (S+T1)**(-1)*S4 - 2*S*U1**(-1)*MS2 + 6*S*M2*MG2*SYMBX - 4*S*
     +    M2*(S+T1)**(-1) - 8*S*M2*SYMBX*S4 - 4*S*MS2*(S+T1)**(-1) - 8*
     +    S*MS2*SYMBX*S4 - 8*S**2*M2*(S+T1)**(-1)*(TG+UG)**(-1) + 4*
     +    S**2*M2*SYMBX - 8*S**2*MS2*(S+T1)**(-1)*(TG+UG)**(-1) + 4*
     +    S**2*MS2*SYMBX - 4*TG*M2*MG2*(S+T1)**(-1)*(S+2*M2)**(-1) + 14
     +    *TG*M2*MG2*(S+2*M2)**(-1)*SYMBX*S4 + 2*TG*M2*MG2*SYMBX + 2*TG
     +    *M2*(S+T1)**(-1)*(S+2*M2)**(-1)*S4 + 2*TG*M2*(S+T1)**(-1) + 2
     +    *TG*M2*(S+2*M2)**(-1)*SYMBX*S4**2 - 12*TG*M2*SYMBX*S4 + 6*TG*
     +    M2**2*MG2*(S+2*M2)**(-1)*SYMBX + 2*TG*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1)*S4 + 2*TG*MS2*(S+T1)**(-1) + 2*TG*MS2*
     +    (S+2*M2)**(-1)*SYMBX*S4**2 - 12*TG*MS2*SYMBX*S4 - 8*TG**2*M2*
     +    MG2*(S+2*M2)**(-1)*SYMBX - 2*TG**2*M2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1) )
     +
      M2QQH3 = M2QQH3 + ANG2(78)*N*CF * (  - 4*TG**2*M2*(S+2*M2)**(-1)*
     +    SYMBX*S4 + 8*TG**2*M2*SYMBX - 2*TG**2*MS2*(S+T1)**(-1)*
     +    (S+2*M2)**(-1) - 4*TG**2*MS2*(S+2*M2)**(-1)*SYMBX*S4 + 8*
     +    TG**2*MS2*SYMBX + 2*TG**3*M2*(S+2*M2)**(-1)*SYMBX + 2*TG**3*
     +    MS2*(S+2*M2)**(-1)*SYMBX - 2*T1**(-1)*M2*MG2*(S+2*M2)**(-1)*
     +    S4 - 4*T1**(-1)*M2*MG2 + 2*T1**(-1)*M2*S4 + 6*T1**(-1)*M2**2*
     +    MG2*(S+2*M2)**(-1) + 2*T1**(-1)*MS2*S4 - 4*U1**(-1)*M2*MG2*
     +    (S+T1)**(-1)*S4 - 2*U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 + 6*
     +    U1**(-1)*M2**2*MG2*(S+2*M2)**(-1) + 4*M2*MG2*(S+T1)**(-1) - 6
     +    *M2*MG2*(S+2*M2)**(-1)*SYMBX*S4**2 + 14*M2*MG2*(S+2*M2)**(-1)
     +     - 2*M2*(S+2*M2)**(-1)*S4 + 4*M2*SYMBX*S4**2 - 8*M2 - 6*M2**2
     +    *MG2*(S+2*M2)**(-1)*SYMBX*S4 + 2*M2**2*MG2*SYMBX - 2*MS2*
     +    (S+2*M2)**(-1)*S4 + 4*MS2*SYMBX*S4**2 - 8*MS2 )
     +
      M2QQH3 = M2QQH3 + ANG2(78)*CF**2 * (  - 4*S*T1**(-1)*M2 - 4*S*
     +    T1**(-1)*MS2 - 4*S*U1**(-1)*M2*(S+T1)**(-1)*S4 + 4*S*U1**(-1)
     +    *M2 - 4*S*U1**(-1)*MS2*(S+T1)**(-1)*S4 + 4*S*U1**(-1)*MS2 + 4
     +    *S*M2*(S+T1)**(-1) + 4*S*MS2*(S+T1)**(-1) + 8*S**2*M2*
     +    (S+T1)**(-1)*(TG+UG)**(-1) + 8*S**2*MS2*(S+T1)**(-1)*
     +    (TG+UG)**(-1) - 4*TG*M2*(S+T1)**(-1) - 4*TG*MS2*(S+T1)**(-1)
     +     + 4*T1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 + 8*T1**(-1)*M2*MG2 - 
     +    4*T1**(-1)*M2*S4 - 12*T1**(-1)*M2**2*MG2*(S+2*M2)**(-1) - 4*
     +    T1**(-1)*MS2*S4 + 8*U1**(-1)*M2*MG2*(S+T1)**(-1)*S4 + 4*
     +    U1**(-1)*M2*MG2*(S+2*M2)**(-1)*S4 - 12*U1**(-1)*M2**2*MG2*
     +    (S+2*M2)**(-1) - 8*M2*MG2*(S+T1)**(-1) - 8*M2*MG2*
     +    (S+2*M2)**(-1) + 4*M2 + 16*M2**2*MG2*(S+2*M2)**(-2) + 4*MS2 )
     +
      M2QQH3 = M2QQH3 + ANG2(80)*N*CF * (  - 4*T1**(-1)*U1**(-1)*MS2*S4
     +     - 4*T1**(-1)*U1**(-1)*S4**2 + 4*T1**(-1)*MS2 + 4*T1**(-1)*S4
     +     )
     +
      M2QQH3 = M2QQH3 + ANG2(80)*CF**2 * ( 4*T1**(-1)*U1**(-1)*MS2*S4
     +     + 4*T1**(-1)*U1**(-1)*S4**2 - 4*T1**(-1)*MS2 - 4*T1**(-1)*S4
     +     - 8*U1**(-2)*M2*MG2 + 8*U1**(-2)*M2*(S+T1)**(-1)*S4**2 - 8*
     +    U1**(-2)*M2*S4 + 8*U1**(-2)*M2**2 - 16*U1**(-1)*M2*
     +    (S+T1)**(-1)*S4 - 4*U1**(-1)*MS2*(S+T1)**(-1)*S4 + 8*U1**(-1)
     +    *MS2 - 4*U1**(-1)*(S+T1)**(-1)*S4**2 + 8*U1**(-1)*S4 + 8*M2*
     +    (S+T1)**(-1) + 4*MS2*(S+T1)**(-1) + 4*(S+T1)**(-1)*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(85)*N*CF * ( 10*S*TG*M2*SYMBX + 10*S*TG*
     +    MS2*SYMBX - 4*S*T1**(-1)*M2*(S+U1)**(-1)*S4 - 4*S*T1**(-1)*M2
     +     + 2*S*U1**(-1)*M2*(S+T1)**(-1)*S4 - 2*S*U1**(-1)*M2 + 2*S*
     +    U1**(-1)*MS2*(S+T1)**(-1)*S4 - 2*S*U1**(-1)*MS2 + 6*S*M2*MG2*
     +    SYMBX + 6*S*M2*(S+U1)**(-1) - 8*S*M2*SYMBX*S4 + 4*S*MS2*
     +    (S+U1)**(-1) - 8*S*MS2*SYMBX*S4 + 2*S*(S+U1)**(-1)*S4 + 4*
     +    S**2*M2*SYMBX + 4*S**2*MS2*SYMBX - 8*TG*M2*MG2*(S+U1)**(-1)*
     +    (S+2*M2)**(-1) + 14*TG*M2*MG2*(S+2*M2)**(-1)*SYMBX*S4 + 2*TG*
     +    M2*MG2*SYMBX - 2*TG*M2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4 + 8*TG*
     +    M2*(S+U1)**(-1) + 2*TG*M2*(S+2*M2)**(-1)*SYMBX*S4**2 - 12*TG*
     +    M2*SYMBX*S4 + 6*TG*M2**2*MG2*(S+2*M2)**(-1)*SYMBX - 2*TG*MS2*
     +    (S+U1)**(-1)*(S+2*M2)**(-1)*S4 + 6*TG*MS2*(S+U1)**(-1) + 2*TG
     +    *MS2*(S+2*M2)**(-1)*SYMBX*S4**2 - 12*TG*MS2*SYMBX*S4 + 2*TG*
     +    (S+U1)**(-1)*S4 - 8*TG**2*M2*MG2*(S+2*M2)**(-1)*SYMBX + 2*
     +    TG**2*M2*(S+U1)**(-1)*(S+2*M2)**(-1) - 4*TG**2*M2*
     +    (S+2*M2)**(-1)*SYMBX*S4 )
     +
      M2QQH3 = M2QQH3 + ANG2(85)*N*CF * ( 8*TG**2*M2*SYMBX + 2*TG**2*
     +    MS2*(S+U1)**(-1)*(S+2*M2)**(-1) - 4*TG**2*MS2*(S+2*M2)**(-1)*
     +    SYMBX*S4 + 8*TG**2*MS2*SYMBX + 2*TG**3*M2*(S+2*M2)**(-1)*
     +    SYMBX + 2*TG**3*MS2*(S+2*M2)**(-1)*SYMBX + 4*T1**(-1)*M2*
     +    (S+U1)**(-1)*S4**2 + 4*T1**(-1)*M2*S4 - 2*U1**(-1)*M2*MG2 + 4
     +    *U1**(-1)*M2*(S+T1)**(-1)*S4**2 - 4*U1**(-1)*M2*S4 + 4*
     +    U1**(-1)*M2**2 + 6*M2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1)*S4 - 2*
     +    M2*MG2*(S+U1)**(-1) - 6*M2*MG2*(S+2*M2)**(-1)*SYMBX*S4**2 - 2
     +    *M2*(S+T1)**(-1)*S4 - 8*M2*(S+U1)**(-1)*S4 + 4*M2*SYMBX*S4**2
     +     - 4*M2 + 6*M2**2*MG2*(S+U1)**(-1)*(S+2*M2)**(-1) - 6*M2**2*
     +    MG2*(S+2*M2)**(-1)*SYMBX*S4 + 2*M2**2*MG2*SYMBX + 2*M2**2*
     +    (S+U1)**(-1) - 4*MS2*(S+U1)**(-1)*S4 + 4*MS2*SYMBX*S4**2 + 2*
     +    MS2 - 2*(S+T1)**(-1)*S4**2 - 2*(S+U1)**(-1)*S4**2 + 2*S4 )


      END IF

      Q2MS2 = (MS + MG)**2/4.D0/MS**2

      IF (IFL.EQ.1) M2QQH = 2.D0*( M2QQH1 + M2QQH2 + M2QQH3)
      IF (IFL.EQ.0) M2QQH = 2.D0*( M2QQH1 + M2QQH2)

      DSGQQH = S4/(S4+MS2)/2.D0 *
     +     ALPHAS**3 * AVG * M2QQH /4.D0 /S**2 *CONV
     +     + DSGQQ3(ALPHAS,S,TG,S4,MS,MG,Q2MS2,IFL)

      RETURN
      END



      REAL*8 FUNCTION DSGQQT(ALPHAS,S,TG,S4,S3,MS,MG,IFL)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR Q +Q -> SQ + GL +Q
C***  THE 1/S3**2 PART

      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,S4G,S3,T1,TG,U1,UG,MS,MG,MS2,MG2,M2
      REAL*8 NS,CONV,N,AVG,CF, M2QQT, M2QQT1, M2QQT2, M2QQT3
      REAL*8 ANGDEF(1:11), ANA(2:2,1:9), ANB(2:2,1:9), ANC(2:2,1:9)
      REAL*8 ANGS3(1:10),XX(5:9),YY2(5:9),XPHI
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -TG
      S4G= S4 -M2
      T1 = TG +M2
      UG = U1 -M2

      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +UG)/ANGDEF(1)
      ANGDEF(3) = (S +TG)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(TG +UG +2.D0*MG2)/ANGDEF(1)
      ANGDEF(7) = SQRT((TG +UG)**2 -4.D0*MG2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (TG*S4G -S*(UG+2.D0*MG2))/(S+TG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (UG*S4G -S*(TG+2.D0*MG2))/(S+UG)/SQRT((TG+UG)**2-4.D0*MG2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6) +M2
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +MS2 +MG2
      ANB(2,2) = -ANB(2,1)
      ANC(2,2) = -ANC(2,1)
      ANA(2,3) = -2.D0*ANGDEF(3)*ANGDEF(4)
      ANB(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,3) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(2,4) =
     +     2.D0*ANGDEF(4)*ANGDEF(7) -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,4) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -ANB(2,3)
      ANC(2,5) = -ANC(2,3)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(2,6) = -ANB(2,4)
      ANC(2,6) = -ANC(2,4)
      ANA(2,7) = +ANA(2,1) -M2
      ANB(2,7) = +ANB(2,1)
      ANC(2,7) = +ANC(2,1)
      ANA(2,8) = +ANA(2,5) -M2
      ANB(2,8) = +ANB(2,5)
      ANC(2,8) = +ANC(2,5)
      ANA(2,9) = +ANA(2,6) -M2
      ANB(2,9) = +ANB(2,6)
      ANC(2,9) = +ANC(2,6)

      XPHI = (S3 -ANA(2,1))/ANB(2,1)

      XX(5) = (ANA(2,5) + ANB(2,5)*XPHI)
      YY2(5)= ANC(2,5)**2 * (1.D0-XPHI**2)
      XX(6) = (ANA(2,6) + ANB(2,6)*XPHI)
      YY2(6)= ANC(2,6)**2 * (1.D0-XPHI**2)
      XX(8) = (ANA(2,8) + ANB(2,8)*XPHI)
      YY2(8)= ANC(2,8)**2 * (1.D0-XPHI**2)
      XX(9) = (ANA(2,9) + ANB(2,9)*XPHI)
      YY2(9)= ANC(2,9)**2 * (1.D0-XPHI**2)

      ANGS3(1) = -XX(6)/(XX(6)**2 - YY2(6))**(1.5D0)
      ANGS3(2) = -1.D0/SQRT(XX(6)**2 - YY2(6))
      ANGS3(3) = -XX(5)/(XX(5)**2 - YY2(5))**(1.5D0)
      ANGS3(4) = -1.D0/SQRT(XX(5)**2 - YY2(5))
      ANGS3(5) = XX(5)
      ANGS3(6) = XX(5)**2 +0.5D0*YY2(5)
      ANGS3(7) = -XX(9)/(XX(9)**2 - YY2(9))**(1.5D0)
      ANGS3(8) = -1.D0/SQRT(XX(9)**2 - YY2(9))
      ANGS3(9) = -XX(8)/(XX(8)**2 - YY2(8))**(1.5D0)
      ANGS3(10)= -1.D0/SQRT(XX(8)**2 - YY2(8))

      M2QQT1 = 0.D0
      M2QQT1 = M2QQT1 + N*CF**2 * ( 4*M2 )
     +
      M2QQT1 = M2QQT1 + ANGS3(7)*N*CF**2 * ( 4*M2**3 )
     +
      M2QQT1 = M2QQT1 + ANGS3(8)*N*CF**2 * ( 4*S*M2 + 8*M2**2 )

      M2QQT2 = 0.D0
      M2QQT2 = M2QQT2 + N*CF**2 * ( 4*M2 )
     +
      M2QQT2 = M2QQT2 + ANGS3(9)*N*CF**2 * ( 4*M2**3 )
     +
      M2QQT2 = M2QQT2 + ANGS3(10)*N*CF**2 * ( 4*S*M2 + 8*M2**2 )

      M2QQT3 = 0.D0
      M2QQT3 = M2QQT3 + ANGS3(8)*CF**2 * (  - 8*M2*MG2 + 16*M2**2*MG2*
     +    (S+2*M2)**(-1) )
     +
      M2QQT3 = M2QQT3 + ANGS3(10)*CF**2 * (  - 8*M2*MG2 + 16*M2**2*MG2*
     +    (S+2*M2)**(-1) )

      IF (IFL.EQ.1) M2QQT = 2.D0*( M2QQT1 + M2QQT2 + M2QQT3)
      IF (IFL.EQ.0) M2QQT = 2.D0*( M2QQT1 + M2QQT2)

      DSGQQT = ALPHAS**3 *AVG *M2QQT *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END
















