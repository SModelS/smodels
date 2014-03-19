C**********************************************************************
C***
C***  THIS ARE THE FUNCTIONS FOR THE SQUARK-SQUARK PROCESS
C***
C**********************************************************************
      REAL*8 FUNCTION DSSQQ1(ALPHAS,S,T1,MS,MG,SCA,IFL)
C***  GIVES THE SCALE DEPENDENCE OF VIRTUAL
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,U1,MS,MG,PI
      REAL*8 DSSQQB,BETAL,PQQ,NS,SCA
      INTEGER IFL
      PI = 4.D0*ATAN(1.D0)
      NS = 6.D0
      U1 = -S -T1
      BETAL = -11.D0 + 2.D0/3.D0*(NS -1.D0)
      PQQ = 4.D0/3.D0*(-2.D0*LOG(-T1/MS**2) -2*LOG(-U1/MS**2) +3.D0)
      DSSQQ1 = DSSQQB(ALPHAS,S,T1,MS,MG,IFL)*
     +     ALPHAS/2.D0/PI*LOG(1.D0/SCA)*( BETAL +PQQ)
      RETURN
      END


      REAL*8 FUNCTION DSSQQ2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA,IFL)
C***  GIVES THE SCALE DEPENDENCE OF LOG(DEL)
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,S4,S4MAX,MS,MG,PI
      REAL*8 DSSQQB,SCA,DEL,DLDEL1
      INTEGER IFL
      PI = 4.D0*ATAN(1.D0)
      DLDEL1 = LOG(S4MAX/MS**2) -(S4MAX -DEL)/S4
      DSSQQ2 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) *
     +      16.D0/3.D0*DLDEL1*DSSQQB(ALPHAS,S,T1,MS,MG,IFL) 
      RETURN
      END


      REAL*8 FUNCTION DSSQQ3(ALPHAS,S,T1,S4,MS,MG,SCA,IFL)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,U1,S4,MS,MG,PI
      REAL*8 DSSQQB,PQQ1,PQQ2,SCA,X1,X2
      INTEGER IFL
      PI = 4.D0*ATAN(1.D0)
      U1 = S4 -S -T1
      X1 = -U1/(S+T1)
      X2 = -T1/(S+U1)
      PQQ1 = 4.D0/3.D0*(1.D0 +X1**2)/(1.D0 -X1)
      PQQ2 = 4.D0/3.D0*(1.D0 +X2**2)/(1.D0 -X2)
      DSSQQ3 = ALPHAS/2.D0/PI*LOG(1.D0/SCA) *
     +     (-1/U1*DSSQQB(ALPHAS,X1*S,X1*T1,MS,MG,IFL)*PQQ1*X1**2
     +      -1/T1*DSSQQB(ALPHAS,X2*S,T1,MS,MG,IFL)   *PQQ2*X2**2 )
      RETURN
      END


      REAL*8 FUNCTION DSSQQR(ALPHAS,S,T1,MS,MG,IFL)
C***  GIVES THE CHANGE FROM MSBAR TO DRBAR FOR THE YUKAWA COUPLINGS
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,MS,MG,PI
      REAL*8 DSSQQB
      INTEGER IFL
      PI = 4.D0*ATAN(1.D0)
      DSSQQR = ALPHAS * 4.D0/3.D0/PI *DSSQQB(ALPHAS,S,T1,MS,MG,IFL)
      RETURN
      END


      REAL*8 FUNCTION DSSQQB(ALPHAS,S,T1,MS,MG,IFL)
C***  BORN CROSS SECTIONS FOR Q + Q -> SQ + SQ
C***  IFL = 0 <--> Q NOT Q
C***  IFL = 1 <--> Q  =  Q
C***  EQUAL FLAVOR AND EQUAL HELICITY 
C***  HAS A FACTOR 1/2 FOR IDENTICAL PARTICLES 
C***  SUMMED OVER LL +RR +LR +RL
C***  MQPLL?: UNEQUAL FLAVORS, SAME HELICITY ( Q QP --> SQL SQPL )
C***  MQPLR?: UNEQUAL FLAVORS, DIFF HELICITY ( Q QP --> SQL SQPR )
C***  MQQLL?: EQUAL FLAVORS,   SAME HELICITY ( Q Q  --> SQL SQL  )
C***  MQQLR?: EQUAL FLAVORS,   DIFF HELICITY ( Q Q  --> SQL SQR  )

      IMPLICIT NONE
      REAL*8 ALPHAS,S,T,T1,TG,U,U1,UG,SB,S1,MS,MG,MS2,MG2,M2
      REAL*8 MQPLLB, MQPLRB, MQQLLB, MQQLRB, MQQB
      REAL*8 NS,CONV,PI,ZETA2,N,CF,AVG
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      N = 3.D0
      CF = (N**2 -1)/2.D0/N
      AVG = (1.D0/2.D0)**2 * (1.D0/3.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = -S -T1
      TG = T1 -M2
      UG = U1 -M2
      SB = S*SQRT(1 -4.D0*MS2/S)
      S1 = 4.D0*MS2 - S
      T = T1 +MS2
      U = U1 +MS2

      MQPLLB = + N*CF * ( 2*MG**2*S*TG**(-2) )

      MQPLRB = + N*CF * (  - 2*MS**2*S*TG**(-2) + 2*TG**(-2)*T1*U1 )

      MQQLLB = + N*CF * ( 2*MG**2*S*TG**(-2) + 2*MG**2*S*UG**(-2) )
     +           + CF * (  - 4*MG**2*S*TG**(-1)*UG**(-1) )

      MQQLRB = + N*CF * (  - 2*MS**2*S*TG**(-2) + 2*TG**(-2)*T1*U1 ) 
     +         + N*CF * (  - 2*MS**2*S*UG**(-2) + 2*UG**(-2)*T1*U1 )

      IF (IFL.EQ.0) MQQB = 2*MQPLLB + 2*MQPLRB 
      IF (IFL.EQ.1) MQQB = 2*0.5D0*MQQLLB + 1*MQQLRB 

      DSSQQB = ALPHAS**2 * AVG * MQQB * PI/S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSSQQB_NG(ALPHAS,S,T1,MS1,MS2,MG,ISQ1,ISQ2)

      IMPLICIT NONE
      REAL*8 ALPHAS,S,T,T1,TG,U,U2,UG,MS1,MS2,MG
      REAL*8 MQPLLB, MQPLRB, MQQLLB, MQQLRB, MQQB
      REAL*8 NS,CONV,PI,ZETA2,N,CF,AVG
      REAL*8 MMtRR, MMuRR, MMtuRR
      REAL*8 MMtLR, MMuLR, MMtuLR
      REAL*8 MMt,MMu,MMtu
      INTEGER ISQ1,ISQ2

ctp protect the isq/abs(isq) below
      IF ( (ISQ1 .EQ. 0).OR.(ISQ2 .EQ. 0) ) THEN
         PRINT*, " DSSQQB_NG: should not be here ",ISQ1,ISQ2
         CALL HARD_STOP
      END IF 

      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      N = 3.D0
      CF = (N**2 -1)/2.D0/N
      AVG = (1.D0/2.D0)**2 * (1.D0/3.D0)**2

      U2 = -S -T1
      TG = T1 + MS1**2 - MG**2
      UG = U2 + MS2**2 - MG**2 
      T = T1 +MS1**2
      U = U2 +MS2**2

      MMtRR = N*Cf*tg**(-2) * ( 2*mg**2*s )
      MMtLR = N*Cf*tg**(-2) * ( 2*u*t - 2*ms1**2*ms2**2 )

ctp note that MMu and MMt should be identical after integrated over PS
      MMuRR = N*Cf*ug**(-2) * ( 2*mg**2*s )
      MMuLR = N*Cf*ug**(-2) * ( 2*u*t - 2*ms1**2*ms2**2 )

      MMtuRR = Cf*tg**(-1)*ug**(-1) * (  - 4*mg**2*s )
      MMtuRR = MMtuRR + N*Cf*tg**(-2) * ( 2*mg**2*s )
      MMtuRR = MMtuRR + N*Cf*ug**(-2) * ( 2*mg**2*s )

      MMtuLR =          N*Cf*tg**(-2) * ( 2*u*t - 2*ms1**2*ms2**2 )
      MMtuLR = MMtuLR + N*Cf*ug**(-2) * ( 2*u*t - 2*ms1**2*ms2**2 )

ctp same helicity in final state? -> means different gluino mass insertion
ctp averaging factor for identical particles in the final state
      IF ( (ISQ1/ABS(ISQ1)) .EQ. ISQ2/ABS(ISQ2) ) THEN
         MMt  = MMtRR
         MMtu = MMtuRR/2.D0
      ELSE
         MMt  = MMtLR
         MMtu = MMtuLR
      END IF
      
ctp same flavor in final state? -> means t and u channel exist
      IF ( ABS(ISQ1) .EQ. ABS(ISQ2) ) THEN
         MQQB = MMtu
      ELSE 
         MQQB = MMt
      END IF 

      DSSQQB_NG = ALPHAS**2 * AVG * MQQB * PI/S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSSQQV(ALPHAS,S,T1,MS,MG,MT,IFL)
C***  VIRTUAL+SOFT CROSS SECTIONS FOR Q + Q -> SQ + SQ
C***  IFL = 0 <--> Q NOT Q
C***  IFL = 1 <--> Q  =  Q
C***  EQUAL FLAVOR AND EQUAL HELICITY 
C***  HAS A FACTOR 1/2 FOR IDENTICAL PARTICLES 
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T,T1,TG,U,U1,UG,SB,S1,MS,MG,MT,MS2,MG2,MT2,M2
      REAL*8 MQPLLV, MQPLRV, MQQLLV, MQQLRV, MQQV
      REAL*8 NS,CONV,PI,ZETA2,N,CF,AVG,BETA,XS, SPENCE
      REAL*8 SK1B0A(1:3), SK1B0B(1:4), SK1B0C(1:2), SK1B0D(1:3,1:2)
      REAL*8 SK1BP(1:4), SK1C0A(1:6), SK1C0B(1:4)
      REAL*8 SK1C0C(1:6,1:2), SK1D0(1:9,1:2), SK1B0E(1:3)
      REAL*8 SL1B0A, SL1B0B, SL1B0C, SL1B0D, SL1B0E
      REAL*8 SL1BP, SL1C0A, SL1C0B
      REAL*8 SL1C0C, SL1D0
      REAL*8 SOF1(1:8)
      REAL*8 DSSQQB, DSSQQR
      COMPLEX*16 KBETAG,KXG,KBETAT,KXT            
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      N = 3.D0
      CF = (N**2 -1)/2.D0/N
      AVG = (1.D0/2.D0)**2 * (1.D0/3.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      MT2 = MT**2
      M2 = MG2 - MS2
      U1 = -S -T1
      TG = T1 -M2
      UG = U1 -M2
      SB = S*SQRT(1 -4.D0*MS2/S)
      S1 = 4.D0*MS2 - S
      T = T1 +MS2
      U = U1 +MS2
      BETA = SQRT(1 -4*MS2/S)
      XS = (1.D0 -BETA)/(1.D0 +BETA)

      KBETAG = SQRT( DCMPLX(1.D0 - 4.D0*MG2/S)) 
      KXG = (DCMPLX(1.D0) -KBETAG)*(DCMPLX(1.D0)+KBETAG)**(-1)
      KBETAT = SQRT( DCMPLX(1.D0 - 4.D0*MT2/S)) 
      KXT = (DCMPLX(1.D0) -KBETAT)*(DCMPLX(1.D0)+KBETAT)**(-1)
      

      
      SK1B0A(1) = SL1B0A(MS2,MG2,MT2,1)
      SK1B0A(2) = SL1B0A(MS2,MG2,MT2,2)
      SK1B0A(3) = SL1B0A(MS2,MG2,MT2,3)

      SK1B0B(1) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,1)
      SK1B0B(2) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,2)
      SK1B0B(3) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,3)
      SK1B0B(4) = SL1B0B(MS2,MG2,MT2,S,BETA,XS,KBETAG,KXG,KBETAT,KXT,4)

      SK1B0C(1) = SL1B0C(MS2,MG2,MT2,1)
      SK1B0C(2) = SL1B0C(MS2,MG2,MT2,2)

      SK1B0D(1,1) = SL1B0D(MS2,MG2,MT2,T,1)
      SK1B0D(2,1) = SL1B0D(MS2,MG2,MT2,T,2)
      SK1B0D(3,1) = SL1B0D(MS2,MG2,MT2,T,3)
      SK1B0D(1,2) = SL1B0D(MS2,MG2,MT2,U,1)
      SK1B0D(2,2) = SL1B0D(MS2,MG2,MT2,U,2)
      SK1B0D(3,2) = SL1B0D(MS2,MG2,MT2,U,3)

      SK1B0E(1) = SL1B0E(MS2,MG2,MT2,1)
      SK1B0E(2) = SL1B0E(MS2,MG2,MT2,2)
      SK1B0E(3) = SL1B0E(MS2,MG2,MT2,3)

      SK1BP(1) = SL1BP(MS2,MG2,MT2,1)
      SK1BP(2) = SL1BP(MS2,MG2,MT2,2)
      SK1BP(3) = SL1BP(MS2,MG2,MT2,3)
      SK1BP(4) = SL1BP(MS2,MG2,MT2,4)


      IF (MS.EQ.MG) THEN

         SK1C0A(1) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,1)
         SK1C0A(2) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,2)
C$$$         SK1C0A(3) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,3)
C$$$         SK1C0A(4) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,4)
C$$$         SK1C0A(5) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,5)
C$$$         SK1C0A(6) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,6)
         
         SK1C0B(1) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,1)
C$$$         SK1C0B(2) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,2)
         SK1C0B(3) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,3)
C$$$         SK1C0B(4) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,4)
         
         SK1C0C(1,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$         SK1C0C(2,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,2)
         SK1C0C(3,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,3)
C$$$         SK1C0C(4,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,4)
C$$$         SK1C0C(5,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,5)
C$$$         SK1C0C(6,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,6)
         
         SK1C0C(1,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,1)
C$$$         SK1C0C(2,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,2)
         SK1C0C(3,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,3)
C$$$         SK1C0C(4,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,4)
C$$$         SK1C0C(5,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,5)
C$$$         SK1C0C(6,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,6)
         
         SK1D0(1,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,1)
         SK1D0(2,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,2)
         SK1D0(3,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,3)
C$$$         SK1D0(4,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,4)
C$$$         SK1D0(5,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,5)
C$$$         SK1D0(6,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,6)
C$$$         SK1D0(7,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,7)
C$$$         SK1D0(8,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,8)
C$$$         SK1D0(9,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,9)
         
         SK1D0(1,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,1)
         SK1D0(2,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,2)
         SK1D0(3,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,3)
C$$$         SK1D0(4,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,4)
C$$$         SK1D0(5,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,5)
C$$$         SK1D0(6,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,6)
C$$$         SK1D0(7,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,7)
C$$$         SK1D0(8,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,8)
C$$$         SK1D0(9,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,9)

      ELSE

         SK1C0A(1) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,1)
C$$$         SK1C0A(2) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,2)
C$$$         SK1C0A(3) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,3)
C$$$         SK1C0A(4) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,4)
         SK1C0A(5) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,5)
         SK1C0A(6) = SL1C0A(MS2,MG2,MT2,S,XS,ZETA2,KXG,KXT,6)
         
         SK1C0B(1) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,1)
C$$$         SK1C0B(2) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,2)
C$$$         SK1C0B(3) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,3)
         SK1C0B(4) = SL1C0B(MS2,MG2,MT2,S,XS,ZETA2,SB,4)
         
C$$$         SK1C0C(1,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,1)
C$$$         SK1C0C(2,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,2)
         SK1C0C(3,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,3)
         SK1C0C(4,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,4)
         SK1C0C(5,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,5)
         SK1C0C(6,1) = SL1C0C(MS2,MG2,MT2,T,ZETA2,6)
         
C$$$         SK1C0C(1,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,1)
C$$$         SK1C0C(2,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,2)
         SK1C0C(3,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,3)
         SK1C0C(4,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,4)
         SK1C0C(5,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,5)
         SK1C0C(6,2) = SL1C0C(MS2,MG2,MT2,U,ZETA2,6)
         
C$$$         SK1D0(1,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,1)
C$$$         SK1D0(2,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,2)
C$$$         SK1D0(3,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,3)
         SK1D0(4,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,4)
C$$$         SK1D0(5,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,5)
C$$$         SK1D0(6,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,6)
         SK1D0(7,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,7)
         SK1D0(8,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,8)
C$$$         SK1D0(9,1) = SL1D0(MS2,MG2,MT2,S,T,U,SB,XS,ZETA2,9)
         
C$$$         SK1D0(1,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,1)
C$$$         SK1D0(2,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,2)
C$$$         SK1D0(3,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,3)
         SK1D0(4,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,4)
C$$$         SK1D0(5,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,5)
C$$$         SK1D0(6,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,6)
         SK1D0(7,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,7)
         SK1D0(8,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,8)
C$$$         SK1D0(9,2) = SL1D0(MS2,MG2,MT2,S,U,T,SB,XS,ZETA2,9)


      END IF

      SOF1(1) = 0.125D0/SB * ( +4*SPENCE(XS) +4*SPENCE(-XS) 
     +     +4*LOG(XS)*LOG(1-XS**2)-2*LOG(XS)**2 -2*ZETA2 )
      SOF1(2) = 0.125D0/MS2* (+2*(S -2*MS2)/SB*LOG(XS))
      SOF1(3) = 0.125D0/MS2*(-2)
      SOF1(4) = 0.125D0/U1*(-1.5D0*ZETA2)
      SOF1(5) = 0.125D0/T1*(-1.5D0*ZETA2)
      SOF1(6) = 0.125D0/T1*(-2*LOG(XS)**2 +LOG(XS*U1/T1)**2 
     +     +2*SPENCE(1 -T1/XS/U1) -2*SPENCE(1 -U1/XS/T1) -1.5D0*ZETA2 )
      SOF1(7) = 0.125D0/U1*(-2*LOG(XS)**2 +LOG(XS*T1/U1)**2 
     +     +2*SPENCE(1 -U1/XS/T1) -2*SPENCE(1 -T1/XS/U1) -1.5D0*ZETA2)
      SOF1(8) = 0.125D0/S*(LOG(T1*U1/MS2/S)**2
     +     +2*SPENCE(1 -MS2*S/T1/U1) -3*ZETA2 )


      IF (MS.EQ.MG) THEN

         IF (IFL.EQ.0) THEN
      MQPLLV = 0.D0
      MQPLLV = MQPLLV + N*CF * (  - 8*S*T**(-1)*T1**(-3)*NS*MS2**3 + 8*
     +    S*T**(-1)*T1**(-3)*MS2**2*MT2 - 8*S*T**(-1)*T1**(-2)*NS*
     +    MS2**2 + 8*S*T**(-1)*T1**(-2)*MS2*MT2 + 8*S*T1**(-3)*NS*
     +    MS2**2 - 8*S*T1**(-3)*MS2*MT2 + 4*S*T1**(-2)*NS*MS2 - 4*S*
     +    T1**(-2)*MT2 )
     +
      MQPLLV = MQPLLV + N**2*CF * ( 8*S*T**(-1)*T1**(-3)*MS2**3 + 8*S*
     +    T**(-1)*T1**(-2)*MS2**2 - 8*S*T1**(-3)*MS2**2 - 8*S*T1**(-2)*
     +    MS2 )
     +
      MQPLLV = MQPLLV + SK1B0A(2)*N*CF * ( 8*S*T**(-1)*T1**(-3)*MS2**2*
     +    MT2 + 8*S*T**(-1)*T1**(-2)*MS2*MT2 - 8*S*T1**(-3)*MS2*MT2 - 4
     +    *S*T1**(-2)*MT2 )
     +
      MQPLLV = MQPLLV + SK1B0C(1)*N*CF**2 * ( 32*S*T1**(-3)*MS2**2 + 24
     +    *S*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0C(1)*N**2*CF * (  - 16*S*T1**(-3)*MS2**2
     +     - 8*S*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N*CF * ( 8*S*T**(-1)*T1**(-3)*NS*
     +    MS2**3 - 8*S*T**(-1)*T1**(-3)*MS2**3 + 8*S*T**(-1)*T1**(-2)*
     +    NS*MS2**2 - 8*S*T**(-1)*T1**(-2)*MS2**2 - 8*S*T1**(-3)*NS*
     +    MS2**2 + 8*S*T1**(-3)*MS2**2 - 8*S*T1**(-2)*NS*MS2 + 8*S*
     +    T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N*CF**2 * (  - 32*S*T1**(-3)*MS2**2
     +     - 32*S*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N**2*CF * (  - 8*S*T**(-1)*T1**(-3)
     +    *MS2**3 - 8*S*T**(-1)*T1**(-2)*MS2**2 + 24*S*T1**(-3)*MS2**2
     +     + 24*S*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0D(3,1)*N*CF * (  - 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 + 8*S*T**(-1)*T1**(-3)*MS2**3 - 8*S*T**(-1)*
     +    T1**(-2)*MS2*MT2 + 8*S*T**(-1)*T1**(-2)*MS2**2 - 8*S*T1**(-3)
     +    *MS2**2 - 8*S*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SK1B0E(3)*N*CF * ( 8*S*T1**(-3)*MS2*MT2 + 4*S*
     +    T1**(-2)*MT2 )
     +
      MQPLLV = MQPLLV + SK1BP(1)*N*CF**2 * ( 16*S*T1**(-2)*MS2**2 )
     +
      MQPLLV = MQPLLV + SK1C0C(1,1)*N*CF**2 * (  - 16*S*T1**(-2)*MS2**2
     +     - 16*S*T1**(-1)*MS2 + 16*MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(1,1)*N**2*CF * ( 16*S*T1**(-2)*MS2**2 + 
     +    16*S*T1**(-1)*MS2 - 4*MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(1,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 
     +    16*MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(1,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 + 4*
     +    MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,1)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 + 
     +    16*MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,1)*N**2*CF * ( 16*S*T1**(-1)*MS2 - 4*
     +    MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 
     +    16*MS2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 + 4*
     +    MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(1,1)*N*CF**2 * ( 16*S**2*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(1,1)*N**2*CF * (  - 8*S**2*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(2,1)*N*CF**2 * (  - 32*S*T1**(-1)*MS2**2
     +     + 16*S**2*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(2,1)*N**2*CF * ( 16*S*T1**(-1)*MS2**2 - 8
     +    *S**2*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(3,1)*N*CF**2 * ( 48*S*MS2 + 32*S**2*
     +    T1**(-1)*MS2 + 16*T1*MS2 )
     +
      MQPLLV = MQPLLV + SK1D0(3,1)*N**2*CF * (  - 12*S*MS2 - 8*S**2*
     +    T1**(-1)*MS2 - 4*T1*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(1)*N*CF**2 * ( 128*S*T1**(-2)*MS2**2 - 64*
     +    S**2*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(1)*N**2*CF * (  - 64*S*T1**(-2)*MS2**2 + 
     +    32*S**2*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(2)*N*CF**2 * (  - 32*S*T1**(-2)*MS2**2 )
     +
      MQPLLV = MQPLLV + SOF1(3)*N*CF**2 * (  - 32*S*T1**(-2)*MS2**2 )
     +
      MQPLLV = MQPLLV + SOF1(4)*N*CF**2 * (  - 64*S*T1**(-1)*MS2 - 64*
     +    S**2*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(4)*N**2*CF * ( 16*S*T1**(-1)*MS2 + 16*S**2
     +    *T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(5)*N*CF**2 * ( 32*S*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(5)*N**2*CF * (  - 16*S*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(6)*N*CF**2 * ( 32*S*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(6)*N**2*CF * (  - 16*S*T1**(-1)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(7)*N*CF**2 * (  - 64*S*T1**(-1)*MS2 - 64*
     +    S**2*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(7)*N**2*CF * ( 16*S*T1**(-1)*MS2 + 16*S**2
     +    *T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(8)*N*CF**2 * (  - 64*S**2*T1**(-2)*MS2 )
     +
      MQPLLV = MQPLLV + SOF1(8)*N**2*CF * ( 32*S**2*T1**(-2)*MS2 )

      MQPLRV = 0.D0
      MQPLRV = MQPLRV + N*CF * ( 8*S*T**(-1)*T1**(-3)*NS*MS2**3 - 8*S*
     +    T**(-1)*T1**(-3)*MS2**2*MT2 + 12*S*T**(-1)*T1**(-2)*NS*MS2**2
     +     - 12*S*T**(-1)*T1**(-2)*MS2*MT2 + 4*S*T**(-1)*T1**(-1)*NS*
     +    MS2 - 4*S*T**(-1)*T1**(-1)*MT2 - 8*S*T1**(-3)*NS*MS2**2 + 8*S
     +    *T1**(-3)*MS2*MT2 - 8*S*T1**(-2)*NS*MS2 + 8*S*T1**(-2)*MT2 + 
     +    8*T**(-1)*T1**(-1)*NS*MS2**2 - 8*T**(-1)*T1**(-1)*MS2*MT2 + 4
     +    *T**(-1)*NS*MS2 - 4*T**(-1)*MT2 - 8*T1**(-1)*NS*MS2 + 8*
     +    T1**(-1)*MT2 )
     +
      MQPLRV = MQPLRV + N**2*CF * ( 4 - 8*S*T**(-1)*T1**(-3)*MS2**3 - 
     +    12*S*T**(-1)*T1**(-2)*MS2**2 - 4*S*T**(-1)*T1**(-1)*MS2 + 8*S
     +    *T1**(-3)*MS2**2 + 12*S*T1**(-2)*MS2 + 4*S*T1**(-1) - 8*
     +    T**(-1)*T1**(-1)*MS2**2 - 4*T**(-1)*MS2 + 8*T1**(-1)*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0A(2)*N*CF * (  - 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 - 12*S*T**(-1)*T1**(-2)*MS2*MT2 - 4*S*T**(-1)*
     +    T1**(-1)*MT2 + 8*S*T1**(-3)*MS2*MT2 + 8*S*T1**(-2)*MT2 - 8*
     +    T**(-1)*T1**(-1)*MS2*MT2 - 4*T**(-1)*MT2 + 8*T1**(-1)*MT2 )
     +
      MQPLRV = MQPLRV + SK1B0C(1)*N*CF**2 * (  - 8 - 32*S*T1**(-3)*
     +    MS2**2 - 40*S*T1**(-2)*MS2 - 8*S*T1**(-1) - 32*T1**(-1)*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0C(1)*N**2*CF * ( 16*S*T1**(-3)*MS2**2 + 16
     +    *S*T1**(-2)*MS2 + 16*T1**(-1)*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N*CF * (  - 4 - 8*S*T**(-1)*
     +    T1**(-3)*NS*MS2**3 + 8*S*T**(-1)*T1**(-3)*MS2**3 - 12*S*
     +    T**(-1)*T1**(-2)*NS*MS2**2 + 12*S*T**(-1)*T1**(-2)*MS2**2 - 4
     +    *S*T**(-1)*T1**(-1)*NS*MS2 + 4*S*T**(-1)*T1**(-1)*MS2 + 8*S*
     +    T1**(-3)*NS*MS2**2 - 8*S*T1**(-3)*MS2**2 + 12*S*T1**(-2)*NS*
     +    MS2 - 12*S*T1**(-2)*MS2 + 4*S*T1**(-1)*NS - 4*S*T1**(-1) - 8*
     +    T**(-1)*T1**(-1)*NS*MS2**2 + 8*T**(-1)*T1**(-1)*MS2**2 - 4*
     +    T**(-1)*NS*MS2 + 4*T**(-1)*MS2 + 8*T1**(-1)*NS*MS2 - 8*
     +    T1**(-1)*MS2 + 4*NS )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N*CF**2 * ( 16 + 32*S*T1**(-3)*
     +    MS2**2 + 48*S*T1**(-2)*MS2 + 16*S*T1**(-1) + 32*T1**(-1)*MS2
     +     )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N**2*CF * (  - 12 + 8*S*T**(-1)*
     +    T1**(-3)*MS2**3 + 12*S*T**(-1)*T1**(-2)*MS2**2 + 4*S*T**(-1)*
     +    T1**(-1)*MS2 - 24*S*T1**(-3)*MS2**2 - 36*S*T1**(-2)*MS2 - 12*
     +    S*T1**(-1) + 8*T**(-1)*T1**(-1)*MS2**2 + 4*T**(-1)*MS2 - 24*
     +    T1**(-1)*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0D(3,1)*N*CF * ( 4 + 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 - 8*S*T**(-1)*T1**(-3)*MS2**3 + 12*S*T**(-1)*
     +    T1**(-2)*MS2*MT2 - 12*S*T**(-1)*T1**(-2)*MS2**2 - 4*S*T**(-1)
     +    *T1**(-1)*MS2 + 4*S*T**(-1)*T1**(-1)*MT2 + 8*S*T1**(-3)*
     +    MS2**2 + 12*S*T1**(-2)*MS2 + 4*S*T1**(-1) + 8*T**(-1)*
     +    T1**(-1)*MS2*MT2 - 8*T**(-1)*T1**(-1)*MS2**2 - 4*T**(-1)*MS2
     +     + 4*T**(-1)*MT2 + 8*T1**(-1)*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0E(3)*N*CF * (  - 8*S*T1**(-3)*MS2*MT2 - 8*
     +    S*T1**(-2)*MT2 - 8*T1**(-1)*MT2 )
     +
      MQPLRV = MQPLRV + SK1BP(1)*N*CF**2 * (  - 16*S*T1**(-2)*MS2**2 - 
     +    16*S*T1**(-1)*MS2 - 16*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0A(1)*N*CF**2 * (  - 8*S - 8*S**2*T1**(-1)
     +     )
     +
      MQPLRV = MQPLRV + SK1C0A(1)*N**2*CF * ( 4*S + 4*S**2*T1**(-1) )
     +
      MQPLRV = MQPLRV + SK1C0A(2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 + 8*
     +    S**2*T1**(-1) )
     +
      MQPLRV = MQPLRV + SK1C0A(2)*N**2*CF * ( 8*S*T1**(-1)*MS2 - 4*S**2
     +    *T1**(-1) )
     +
      MQPLRV = MQPLRV + SK1C0B(1)*N*CF**2 * ( 16*S*T1**(-1)*MS2 - 16*S
     +     - 8*S**2*T1**(-1) + 32*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0B(1)*N**2*CF * (  - 8*S*T1**(-1)*MS2 + 8*S
     +     + 4*S**2*T1**(-1) - 16*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0B(3)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 32
     +    *S*S1**(-1)*MS2 + 16*S + 8*S**2*T1**(-1) + 8*S**2*S1**(-1) )
     +
      MQPLRV = MQPLRV + SK1C0B(3)*N**2*CF * ( 8*S*T1**(-1)*MS2 + 16*S*
     +    S1**(-1)*MS2 - 8*S - 4*S**2*T1**(-1) - 4*S**2*S1**(-1) )
     +
      MQPLRV = MQPLRV + SK1C0C(1,1)*N*CF**2 * ( 16*S*T1**(-2)*MS2**2 + 
     +    16*S*T1**(-1)*MS2 - 16*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(1,1)*N**2*CF * (  - 16*S*T1**(-2)*MS2**2
     +     - 16*S*T1**(-1)*MS2 - 4*S + 4*T1 - 4*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(1,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 + 
     +    32*S + 16*S**2*T1**(-1) + 16*T1 - 16*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(1,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 - 8*S
     +     - 4*S**2*T1**(-1) - 4*T1 + 4*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,1)*N*CF**2 * ( 16*S*T1**(-1)*MS2 + 16*
     +    S + 16*T1 - 16*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,1)*N**2*CF * (  - 16*S*T1**(-1)*MS2 - 
     +    12*S - 12*T1 + 4*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,2)*N*CF**2 * ( 16*S*T1**(-1)*MS2 - 32*
     +    S - 16*S**2*T1**(-1) - 16*T1 + 16*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,2)*N**2*CF * (  - 4*S*T1**(-1)*MS2 + 8
     +    *S + 4*S**2*T1**(-1) + 4*T1 - 4*MS2 )
     +
      MQPLRV = MQPLRV + SK1D0(1,1)*N*CF**2 * (  - 8*S*T1 - 16*S**2*
     +    T1**(-1)*MS2 - 8*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(1,1)*N**2*CF * ( 4*S*T1 + 8*S**2*T1**(-1)
     +    *MS2 + 4*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(2,1)*N*CF**2 * ( 32*S*T1**(-1)*MS2**2 + 
     +    16*S*MS2 - 16*S**2*T1**(-1)*MS2 - 8*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(2,1)*N**2*CF * (  - 16*S*T1**(-1)*MS2**2
     +     - 8*S*MS2 + 8*S**2*T1**(-1)*MS2 + 4*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(3,1)*N*CF**2 * (  - 32*S*T1 - 48*S*MS2 - 
     +    32*S**2*T1**(-1)*MS2 - 16*S**2 - 16*T1*MS2 - 16*T1**2 )
     +
      MQPLRV = MQPLRV + SK1D0(3,1)*N**2*CF * ( 8*S*T1 + 12*S*MS2 + 8*
     +    S**2*T1**(-1)*MS2 + 4*S**2 + 4*T1*MS2 + 4*T1**2 )
     +
      MQPLRV = MQPLRV + SOF1(1)*N*CF**2 * (  - 128*S*T1**(-2)*MS2**2 - 
     +    128*S*T1**(-1)*MS2 + 64*S + 64*S**2*T1**(-2)*MS2 + 64*S**2*
     +    T1**(-1) - 128*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(1)*N**2*CF * ( 64*S*T1**(-2)*MS2**2 + 64*S
     +    *T1**(-1)*MS2 - 32*S - 32*S**2*T1**(-2)*MS2 - 32*S**2*
     +    T1**(-1) + 64*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(2)*N*CF**2 * ( 32*S*T1**(-2)*MS2**2 + 32*S
     +    *T1**(-1)*MS2 + 32*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(3)*N*CF**2 * ( 32*S*T1**(-2)*MS2**2 + 32*S
     +    *T1**(-1)*MS2 + 32*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(4)*N*CF**2 * ( 64*S*T1**(-1)*MS2 + 128*S
     +     + 64*S**2*T1**(-2)*MS2 + 64*S**2*T1**(-1) + 64*T1 )
     +
      MQPLRV = MQPLRV + SOF1(4)*N**2*CF * (  - 16*S*T1**(-1)*MS2 - 32*S
     +     - 16*S**2*T1**(-2)*MS2 - 16*S**2*T1**(-1) - 16*T1 )
     +
      MQPLRV = MQPLRV + SOF1(5)*N*CF**2 * (  - 32*S*T1**(-1)*MS2 - 32*S
     +     - 32*T1 )
     +
      MQPLRV = MQPLRV + SOF1(5)*N**2*CF * ( 16*S*T1**(-1)*MS2 + 16*S + 
     +    16*T1 )
     +
      MQPLRV = MQPLRV + SOF1(6)*N*CF**2 * (  - 32*S*T1**(-1)*MS2 - 32*S
     +     - 32*T1 )
     +
      MQPLRV = MQPLRV + SOF1(6)*N**2*CF * ( 16*S*T1**(-1)*MS2 + 16*S + 
     +    16*T1 )
     +
      MQPLRV = MQPLRV + SOF1(7)*N*CF**2 * ( 64*S*T1**(-1)*MS2 + 128*S
     +     + 64*S**2*T1**(-2)*MS2 + 64*S**2*T1**(-1) + 64*T1 )
     +
      MQPLRV = MQPLRV + SOF1(7)*N**2*CF * (  - 16*S*T1**(-1)*MS2 - 32*S
     +     - 16*S**2*T1**(-2)*MS2 - 16*S**2*T1**(-1) - 16*T1 )
     +
      MQPLRV = MQPLRV + SOF1(8)*N*CF**2 * ( 64*S + 64*S**2*T1**(-2)*MS2
     +     + 64*S**2*T1**(-1) )
     +
      MQPLRV = MQPLRV + SOF1(8)*N**2*CF * (  - 32*S - 32*S**2*T1**(-2)*
     +    MS2 - 32*S**2*T1**(-1) )

         END IF

         IF (IFL.EQ.1) THEN
      MQQLLV = 0.D0
      MQQLLV = MQQLLV + N*CF * (  - 8*S*T**(-1)*T1**(-3)*NS*MS2**3 + 8*
     +    S*T**(-1)*T1**(-3)*MS2**2*MT2 - 8*S*T**(-1)*T1**(-2)*NS*
     +    MS2**2 + 8*S*T**(-1)*T1**(-2)*MS2*MT2 - 8*S*U**(-1)*U1**(-3)*
     +    NS*MS2**3 + 8*S*U**(-1)*U1**(-3)*MS2**2*MT2 - 8*S*U**(-1)*
     +    U1**(-2)*NS*MS2**2 + 8*S*U**(-1)*U1**(-2)*MS2*MT2 + 8*S*
     +    T1**(-3)*NS*MS2**2 - 8*S*T1**(-3)*MS2*MT2 + 4*S*T1**(-2)*NS*
     +    MS2 - 4*S*T1**(-2)*MT2 + 8*S*U1**(-3)*NS*MS2**2 - 8*S*
     +    U1**(-3)*MS2*MT2 + 4*S*U1**(-2)*NS*MS2 - 4*S*U1**(-2)*MT2 + 8
     +    *T**(-1)*T1**(-2)*MS2**3 + 8*T**(-1)*T1**(-1)*U1**(-1)*MS2**3
     +     + 8*T**(-1)*T1**(-1)*MS2**2 + 8*T**(-1)*U1**(-1)*MS2**2 + 8*
     +    U**(-1)*T1**(-1)*U1**(-1)*MS2**3 + 8*U**(-1)*T1**(-1)*MS2**2
     +     + 8*U**(-1)*U1**(-2)*MS2**3 + 8*U**(-1)*U1**(-1)*MS2**2 - 8*
     +    T1**(-2)*MS2**2 - 16*T1**(-1)*U1**(-1)*MS2**2 - 16*T1**(-1)*
     +    MS2 - 8*U1**(-2)*MS2**2 - 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + N*CF**2 * ( 16*T**(-1)*T1**(-2)*NS*MS2**3 - 16*
     +    T**(-1)*T1**(-2)*MS2**2*MT2 + 16*T**(-1)*T1**(-1)*U1**(-1)*NS
     +    *MS2**3 - 16*T**(-1)*T1**(-1)*U1**(-1)*MS2**2*MT2 + 16*
     +    T**(-1)*T1**(-1)*NS*MS2**2 - 16*T**(-1)*T1**(-1)*MS2*MT2 + 16
     +    *T**(-1)*U1**(-1)*NS*MS2**2 - 16*T**(-1)*U1**(-1)*MS2*MT2 + 
     +    16*U**(-1)*T1**(-1)*U1**(-1)*NS*MS2**3 - 16*U**(-1)*T1**(-1)*
     +    U1**(-1)*MS2**2*MT2 + 16*U**(-1)*T1**(-1)*NS*MS2**2 - 16*
     +    U**(-1)*T1**(-1)*MS2*MT2 + 16*U**(-1)*U1**(-2)*NS*MS2**3 - 16
     +    *U**(-1)*U1**(-2)*MS2**2*MT2 + 16*U**(-1)*U1**(-1)*NS*MS2**2
     +     - 16*U**(-1)*U1**(-1)*MS2*MT2 - 16*T1**(-2)*NS*MS2**2 + 16*
     +    T1**(-2)*MS2*MT2 - 32*T1**(-1)*U1**(-1)*NS*MS2**2 + 32*
     +    T1**(-1)*U1**(-1)*MS2*MT2 - 16*T1**(-1)*NS*MS2 + 16*T1**(-1)*
     +    MT2 - 16*U1**(-2)*NS*MS2**2 + 16*U1**(-2)*MS2*MT2 - 16*
     +    U1**(-1)*NS*MS2 + 16*U1**(-1)*MT2 )
     +
      MQQLLV = MQQLLV + N**2*CF * ( 8*S*T**(-1)*T1**(-3)*MS2**3 + 8*S*
     +    T**(-1)*T1**(-2)*MS2**2 + 8*S*U**(-1)*U1**(-3)*MS2**3 + 8*S*
     +    U**(-1)*U1**(-2)*MS2**2 - 8*S*T1**(-3)*MS2**2 - 8*S*T1**(-2)*
     +    MS2 - 8*S*U1**(-3)*MS2**2 - 8*S*U1**(-2)*MS2 - 8*T**(-1)*
     +    T1**(-2)*NS*MS2**3 + 8*T**(-1)*T1**(-2)*MS2**2*MT2 - 8*
     +    T**(-1)*T1**(-1)*U1**(-1)*NS*MS2**3 + 8*T**(-1)*T1**(-1)*
     +    U1**(-1)*MS2**2*MT2 - 8*T**(-1)*T1**(-1)*NS*MS2**2 + 8*
     +    T**(-1)*T1**(-1)*MS2*MT2 - 8*T**(-1)*U1**(-1)*NS*MS2**2 + 8*
     +    T**(-1)*U1**(-1)*MS2*MT2 - 8*U**(-1)*T1**(-1)*U1**(-1)*NS*
     +    MS2**3 + 8*U**(-1)*T1**(-1)*U1**(-1)*MS2**2*MT2 - 8*U**(-1)*
     +    T1**(-1)*NS*MS2**2 + 8*U**(-1)*T1**(-1)*MS2*MT2 - 8*U**(-1)*
     +    U1**(-2)*NS*MS2**3 + 8*U**(-1)*U1**(-2)*MS2**2*MT2 - 8*
     +    U**(-1)*U1**(-1)*NS*MS2**2 + 8*U**(-1)*U1**(-1)*MS2*MT2 + 8*
     +    T1**(-2)*NS*MS2**2 - 8*T1**(-2)*MS2*MT2 + 16*T1**(-1)*
     +    U1**(-1)*NS*MS2**2 - 16*T1**(-1)*U1**(-1)*MS2*MT2 + 8*
     +    T1**(-1)*NS*MS2 )
     +
      MQQLLV = MQQLLV + N**2*CF * (  - 8*T1**(-1)*MT2 + 8*U1**(-2)*NS*
     +    MS2**2 - 8*U1**(-2)*MS2*MT2 + 8*U1**(-1)*NS*MS2 - 8*U1**(-1)*
     +    MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N*CF * ( 8*S*T**(-1)*T1**(-3)*MS2**2*
     +    MT2 + 8*S*T**(-1)*T1**(-2)*MS2*MT2 + 8*S*U**(-1)*U1**(-3)*
     +    MS2**2*MT2 + 8*S*U**(-1)*U1**(-2)*MS2*MT2 - 8*S*T1**(-3)*MS2*
     +    MT2 - 4*S*T1**(-2)*MT2 - 8*S*U1**(-3)*MS2*MT2 - 4*S*U1**(-2)*
     +    MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N*CF**2 * (  - 16*T**(-1)*T1**(-2)*
     +    MS2**2*MT2 - 16*T**(-1)*T1**(-1)*U1**(-1)*MS2**2*MT2 - 16*
     +    T**(-1)*T1**(-1)*MS2*MT2 - 16*T**(-1)*U1**(-1)*MS2*MT2 - 16*
     +    U**(-1)*T1**(-1)*U1**(-1)*MS2**2*MT2 - 16*U**(-1)*T1**(-1)*
     +    MS2*MT2 - 16*U**(-1)*U1**(-2)*MS2**2*MT2 - 16*U**(-1)*
     +    U1**(-1)*MS2*MT2 + 16*T1**(-2)*MS2*MT2 + 32*T1**(-1)*U1**(-1)
     +    *MS2*MT2 + 16*T1**(-1)*MT2 + 16*U1**(-2)*MS2*MT2 + 16*
     +    U1**(-1)*MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N**2*CF * ( 8*T**(-1)*T1**(-2)*MS2**2
     +    *MT2 + 8*T**(-1)*T1**(-1)*U1**(-1)*MS2**2*MT2 + 8*T**(-1)*
     +    T1**(-1)*MS2*MT2 + 8*T**(-1)*U1**(-1)*MS2*MT2 + 8*U**(-1)*
     +    T1**(-1)*U1**(-1)*MS2**2*MT2 + 8*U**(-1)*T1**(-1)*MS2*MT2 + 8
     +    *U**(-1)*U1**(-2)*MS2**2*MT2 + 8*U**(-1)*U1**(-1)*MS2*MT2 - 8
     +    *T1**(-2)*MS2*MT2 - 16*T1**(-1)*U1**(-1)*MS2*MT2 - 8*T1**(-1)
     +    *MT2 - 8*U1**(-2)*MS2*MT2 - 8*U1**(-1)*MT2 )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*N*CF * (  - 16*T1**(-2)*MS2**2 - 32*
     +    T1**(-1)*U1**(-1)*MS2**2 - 16*T1**(-1)*MS2 - 16*U1**(-2)*
     +    MS2**2 - 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*N*CF**2 * ( 32*S*T1**(-3)*MS2**2 + 24
     +    *S*T1**(-2)*MS2 + 32*S*U1**(-3)*MS2**2 + 24*S*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*N**2*CF * (  - 16*S*T1**(-3)*MS2**2
     +     - 8*S*T1**(-2)*MS2 - 16*S*U1**(-3)*MS2**2 - 8*S*U1**(-2)*MS2
     +     )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*CF**2 * ( 32*T1**(-2)*MS2**2 + 64*
     +    T1**(-1)*U1**(-1)*MS2**2 + 48*T1**(-1)*MS2 + 32*U1**(-2)*
     +    MS2**2 + 48*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF * ( 8*S*T**(-1)*T1**(-3)*NS*
     +    MS2**3 - 8*S*T**(-1)*T1**(-3)*MS2**3 + 8*S*T**(-1)*T1**(-2)*
     +    NS*MS2**2 - 8*S*T**(-1)*T1**(-2)*MS2**2 - 8*S*T1**(-3)*NS*
     +    MS2**2 + 8*S*T1**(-3)*MS2**2 - 8*S*T1**(-2)*NS*MS2 + 8*S*
     +    T1**(-2)*MS2 - 8*T**(-1)*T1**(-2)*MS2**3 - 8*T**(-1)*T1**(-1)
     +    *U1**(-1)*MS2**3 - 8*T**(-1)*T1**(-1)*MS2**2 - 8*T**(-1)*
     +    U1**(-1)*MS2**2 + 24*T1**(-2)*MS2**2 + 24*T1**(-1)*U1**(-1)*
     +    MS2**2 + 24*T1**(-1)*MS2 + 24*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF**2 * (  - 32*S*T1**(-3)*MS2**2
     +     - 32*S*T1**(-2)*MS2 - 16*T**(-1)*T1**(-2)*NS*MS2**3 + 16*
     +    T**(-1)*T1**(-2)*MS2**3 - 16*T**(-1)*T1**(-1)*U1**(-1)*NS*
     +    MS2**3 + 16*T**(-1)*T1**(-1)*U1**(-1)*MS2**3 - 16*T**(-1)*
     +    T1**(-1)*NS*MS2**2 + 16*T**(-1)*T1**(-1)*MS2**2 - 16*T**(-1)*
     +    U1**(-1)*NS*MS2**2 + 16*T**(-1)*U1**(-1)*MS2**2 + 16*T1**(-2)
     +    *NS*MS2**2 - 16*T1**(-2)*MS2**2 + 16*T1**(-1)*U1**(-1)*NS*
     +    MS2**2 - 16*T1**(-1)*U1**(-1)*MS2**2 + 16*T1**(-1)*NS*MS2 - 
     +    16*T1**(-1)*MS2 + 16*U1**(-1)*NS*MS2 - 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N**2*CF * (  - 8*S*T**(-1)*T1**(-3)
     +    *MS2**3 - 8*S*T**(-1)*T1**(-2)*MS2**2 + 24*S*T1**(-3)*MS2**2
     +     + 24*S*T1**(-2)*MS2 + 8*T**(-1)*T1**(-2)*NS*MS2**3 - 8*
     +    T**(-1)*T1**(-2)*MS2**3 + 8*T**(-1)*T1**(-1)*U1**(-1)*NS*
     +    MS2**3 - 8*T**(-1)*T1**(-1)*U1**(-1)*MS2**3 + 8*T**(-1)*
     +    T1**(-1)*NS*MS2**2 - 8*T**(-1)*T1**(-1)*MS2**2 + 8*T**(-1)*
     +    U1**(-1)*NS*MS2**2 - 8*T**(-1)*U1**(-1)*MS2**2 - 8*T1**(-2)*
     +    NS*MS2**2 + 8*T1**(-2)*MS2**2 - 8*T1**(-1)*U1**(-1)*NS*MS2**2
     +     + 8*T1**(-1)*U1**(-1)*MS2**2 - 8*T1**(-1)*NS*MS2 + 8*
     +    T1**(-1)*MS2 - 8*U1**(-1)*NS*MS2 + 8*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*CF**2 * (  - 32*T1**(-2)*MS2**2 - 
     +    32*T1**(-1)*U1**(-1)*MS2**2 - 32*T1**(-1)*MS2 - 32*U1**(-1)*
     +    MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N*CF * ( 8*S*U**(-1)*U1**(-3)*NS*
     +    MS2**3 - 8*S*U**(-1)*U1**(-3)*MS2**3 + 8*S*U**(-1)*U1**(-2)*
     +    NS*MS2**2 - 8*S*U**(-1)*U1**(-2)*MS2**2 - 8*S*U1**(-3)*NS*
     +    MS2**2 + 8*S*U1**(-3)*MS2**2 - 8*S*U1**(-2)*NS*MS2 + 8*S*
     +    U1**(-2)*MS2 - 8*U**(-1)*T1**(-1)*U1**(-1)*MS2**3 - 8*U**(-1)
     +    *T1**(-1)*MS2**2 - 8*U**(-1)*U1**(-2)*MS2**3 - 8*U**(-1)*
     +    U1**(-1)*MS2**2 + 24*T1**(-1)*U1**(-1)*MS2**2 + 24*T1**(-1)*
     +    MS2 + 24*U1**(-2)*MS2**2 + 24*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N*CF**2 * (  - 32*S*U1**(-3)*MS2**2
     +     - 32*S*U1**(-2)*MS2 - 16*U**(-1)*T1**(-1)*U1**(-1)*NS*MS2**3
     +     + 16*U**(-1)*T1**(-1)*U1**(-1)*MS2**3 - 16*U**(-1)*T1**(-1)*
     +    NS*MS2**2 + 16*U**(-1)*T1**(-1)*MS2**2 - 16*U**(-1)*U1**(-2)*
     +    NS*MS2**3 + 16*U**(-1)*U1**(-2)*MS2**3 - 16*U**(-1)*U1**(-1)*
     +    NS*MS2**2 + 16*U**(-1)*U1**(-1)*MS2**2 + 16*T1**(-1)*U1**(-1)
     +    *NS*MS2**2 - 16*T1**(-1)*U1**(-1)*MS2**2 + 16*T1**(-1)*NS*MS2
     +     - 16*T1**(-1)*MS2 + 16*U1**(-2)*NS*MS2**2 - 16*U1**(-2)*
     +    MS2**2 + 16*U1**(-1)*NS*MS2 - 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N**2*CF * (  - 8*S*U**(-1)*U1**(-3)
     +    *MS2**3 - 8*S*U**(-1)*U1**(-2)*MS2**2 + 24*S*U1**(-3)*MS2**2
     +     + 24*S*U1**(-2)*MS2 + 8*U**(-1)*T1**(-1)*U1**(-1)*NS*MS2**3
     +     - 8*U**(-1)*T1**(-1)*U1**(-1)*MS2**3 + 8*U**(-1)*T1**(-1)*NS
     +    *MS2**2 - 8*U**(-1)*T1**(-1)*MS2**2 + 8*U**(-1)*U1**(-2)*NS*
     +    MS2**3 - 8*U**(-1)*U1**(-2)*MS2**3 + 8*U**(-1)*U1**(-1)*NS*
     +    MS2**2 - 8*U**(-1)*U1**(-1)*MS2**2 - 8*T1**(-1)*U1**(-1)*NS*
     +    MS2**2 + 8*T1**(-1)*U1**(-1)*MS2**2 - 8*T1**(-1)*NS*MS2 + 8*
     +    T1**(-1)*MS2 - 8*U1**(-2)*NS*MS2**2 + 8*U1**(-2)*MS2**2 - 8*
     +    U1**(-1)*NS*MS2 + 8*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*CF**2 * (  - 32*T1**(-1)*U1**(-1)*
     +    MS2**2 - 32*T1**(-1)*MS2 - 32*U1**(-2)*MS2**2 - 32*U1**(-1)*
     +    MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N*CF * (  - 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 + 8*S*T**(-1)*T1**(-3)*MS2**3 - 8*S*T**(-1)*
     +    T1**(-2)*MS2*MT2 + 8*S*T**(-1)*T1**(-2)*MS2**2 - 8*S*T1**(-3)
     +    *MS2**2 - 8*S*T1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N*CF**2 * ( 16*T**(-1)*T1**(-2)*
     +    MS2**2*MT2 - 16*T**(-1)*T1**(-2)*MS2**3 + 16*T**(-1)*T1**(-1)
     +    *U1**(-1)*MS2**2*MT2 - 16*T**(-1)*T1**(-1)*U1**(-1)*MS2**3 + 
     +    16*T**(-1)*T1**(-1)*MS2*MT2 - 16*T**(-1)*T1**(-1)*MS2**2 + 16
     +    *T**(-1)*U1**(-1)*MS2*MT2 - 16*T**(-1)*U1**(-1)*MS2**2 + 16*
     +    T1**(-2)*MS2**2 + 16*T1**(-1)*U1**(-1)*MS2**2 + 16*T1**(-1)*
     +    MS2 + 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N**2*CF * (  - 8*T**(-1)*T1**(-2)*
     +    MS2**2*MT2 + 8*T**(-1)*T1**(-2)*MS2**3 - 8*T**(-1)*T1**(-1)*
     +    U1**(-1)*MS2**2*MT2 + 8*T**(-1)*T1**(-1)*U1**(-1)*MS2**3 - 8*
     +    T**(-1)*T1**(-1)*MS2*MT2 + 8*T**(-1)*T1**(-1)*MS2**2 - 8*
     +    T**(-1)*U1**(-1)*MS2*MT2 + 8*T**(-1)*U1**(-1)*MS2**2 - 8*
     +    T1**(-2)*MS2**2 - 8*T1**(-1)*U1**(-1)*MS2**2 - 8*T1**(-1)*MS2
     +     - 8*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N*CF * (  - 8*S*U**(-1)*U1**(-3)*
     +    MS2**2*MT2 + 8*S*U**(-1)*U1**(-3)*MS2**3 - 8*S*U**(-1)*
     +    U1**(-2)*MS2*MT2 + 8*S*U**(-1)*U1**(-2)*MS2**2 - 8*S*U1**(-3)
     +    *MS2**2 - 8*S*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N*CF**2 * ( 16*U**(-1)*T1**(-1)*
     +    U1**(-1)*MS2**2*MT2 - 16*U**(-1)*T1**(-1)*U1**(-1)*MS2**3 + 
     +    16*U**(-1)*T1**(-1)*MS2*MT2 - 16*U**(-1)*T1**(-1)*MS2**2 + 16
     +    *U**(-1)*U1**(-2)*MS2**2*MT2 - 16*U**(-1)*U1**(-2)*MS2**3 + 
     +    16*U**(-1)*U1**(-1)*MS2*MT2 - 16*U**(-1)*U1**(-1)*MS2**2 + 16
     +    *T1**(-1)*U1**(-1)*MS2**2 + 16*T1**(-1)*MS2 + 16*U1**(-2)*
     +    MS2**2 + 16*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N**2*CF * (  - 8*U**(-1)*T1**(-1)*
     +    U1**(-1)*MS2**2*MT2 + 8*U**(-1)*T1**(-1)*U1**(-1)*MS2**3 - 8*
     +    U**(-1)*T1**(-1)*MS2*MT2 + 8*U**(-1)*T1**(-1)*MS2**2 - 8*
     +    U**(-1)*U1**(-2)*MS2**2*MT2 + 8*U**(-1)*U1**(-2)*MS2**3 - 8*
     +    U**(-1)*U1**(-1)*MS2*MT2 + 8*U**(-1)*U1**(-1)*MS2**2 - 8*
     +    T1**(-1)*U1**(-1)*MS2**2 - 8*T1**(-1)*MS2 - 8*U1**(-2)*MS2**2
     +     - 8*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N*CF * ( 8*S*T1**(-3)*MS2*MT2 + 4*S*
     +    T1**(-2)*MT2 + 8*S*U1**(-3)*MS2*MT2 + 4*S*U1**(-2)*MT2 )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N*CF**2 * (  - 16*T1**(-2)*MS2*MT2 - 
     +    32*T1**(-1)*U1**(-1)*MS2*MT2 - 16*T1**(-1)*MT2 - 16*U1**(-2)*
     +    MS2*MT2 - 16*U1**(-1)*MT2 )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N**2*CF * ( 8*T1**(-2)*MS2*MT2 + 16*
     +    T1**(-1)*U1**(-1)*MS2*MT2 + 8*T1**(-1)*MT2 + 8*U1**(-2)*MS2*
     +    MT2 + 8*U1**(-1)*MT2 )
     +
      MQQLLV = MQQLLV + SK1BP(1)*N*CF**2 * ( 16*S*T1**(-2)*MS2**2 + 16*
     +    S*U1**(-2)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1BP(1)*CF**2 * ( 32*T1**(-1)*MS2**2 + 32*
     +    U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1C0A(2)*N*CF * ( 8*S*T1**(-1)*MS2 + 8*S*
     +    U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0A(2)*N*CF**2 * ( 16*S*T1**(-1)*MS2 + 16*S*
     +    U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0A(2)*N**2*CF * (  - 8*S*T1**(-1)*MS2 - 8*S
     +    *U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0A(2)*CF**2 * (  - 12*S*T1**(-1)*MS2 - 12*S
     +    *U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,1)*N*CF * (  - 20*S*U1**(-1)*MS2 + 16*
     +    T1**(-1)*MS2**2 + 16*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,1)*N*CF**2 * (  - 16*S*T1**(-2)*MS2**2
     +     - 16*S*T1**(-1)*MS2 - 16*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,1)*N**2*CF * ( 16*S*T1**(-2)*MS2**2 + 
     +    16*S*T1**(-1)*MS2 + 4*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,1)*CF**2 * ( 24*S*U1**(-1)*MS2 - 16*
     +    T1**(-1)*MS2**2 - 16*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,2)*N*CF * (  - 20*S*T1**(-1)*MS2 + 16*
     +    T1**(-1)*MS2**2 + 16*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 
     +    16*S*U1**(-2)*MS2**2 - 16*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-2)*MS2**2 + 16*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(1,2)*CF**2 * ( 24*S*T1**(-1)*MS2 - 16*
     +    T1**(-1)*MS2**2 - 16*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N*CF * (  - 20*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 
     +    16*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N**2*CF * ( 16*S*T1**(-1)*MS2 + 4*S
     +    *U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*CF**2 * ( 24*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N*CF * (  - 20*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 
     +    16*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*CF**2 * ( 24*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,1)*N*CF * ( 8*S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,1)*N*CF**2 * ( 16*S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,1)*N**2*CF * (  - 8*S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,1)*CF**2 * (  - 8*S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,2)*N*CF * ( 8*S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,2)*N*CF**2 * ( 16*S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,2)*N**2*CF * (  - 8*S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(1,2)*CF**2 * (  - 8*S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,1)*N*CF * (  - 16*S*U1**(-1)*MS2**2 + 8
     +    *S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,1)*N*CF**2 * (  - 32*S*T1**(-1)*MS2**2
     +     + 16*S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,1)*N**2*CF * ( 16*S*T1**(-1)*MS2**2 - 8
     +    *S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,1)*CF**2 * ( 16*S*U1**(-1)*MS2**2 - 8*
     +    S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,2)*N*CF * (  - 16*S*T1**(-1)*MS2**2 + 8
     +    *S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,2)*N*CF**2 * (  - 32*S*U1**(-1)*MS2**2
     +     + 16*S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,2)*N**2*CF * ( 16*S*U1**(-1)*MS2**2 - 8
     +    *S**2*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(2,2)*CF**2 * ( 16*S*T1**(-1)*MS2**2 - 8*
     +    S**2*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,1)*N*CF * (  - 8*S*MS2 - 4*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,1)*N*CF**2 * ( 48*S*MS2 + 32*S**2*
     +    T1**(-1)*MS2 + 16*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,1)*N**2*CF * (  - 12*S*MS2 - 8*S**2*
     +    T1**(-1)*MS2 - 4*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,1)*CF**2 * ( 16*S*MS2 + 8*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,2)*N*CF * (  - 4*S*MS2 + 4*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,2)*N*CF**2 * ( 32*S*MS2 + 32*S**2*
     +    U1**(-1)*MS2 - 16*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,2)*N**2*CF * (  - 8*S*MS2 - 8*S**2*
     +    U1**(-1)*MS2 + 4*T1*MS2 )
     +
      MQQLLV = MQQLLV + SK1D0(3,2)*CF**2 * ( 8*S*MS2 - 8*T1*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N*CF * ( 64*S*T1**(-1)*MS2 + 64*S*
     +    U1**(-1)*MS2 - 128*T1**(-1)*MS2**2 - 128*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N*CF**2 * ( 128*S*T1**(-2)*MS2**2 + 128
     +    *S*U1**(-2)*MS2**2 - 64*S**2*T1**(-2)*MS2 - 64*S**2*U1**(-2)*
     +    MS2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N**2*CF * (  - 64*S*T1**(-2)*MS2**2 - 
     +    64*S*U1**(-2)*MS2**2 + 32*S**2*T1**(-2)*MS2 + 32*S**2*
     +    U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*CF**2 * (  - 64*S*T1**(-1)*MS2 - 64*S*
     +    U1**(-1)*MS2 + 128*T1**(-1)*MS2**2 + 128*U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(2)*N*CF**2 * (  - 32*S*T1**(-2)*MS2**2 - 
     +    32*S*U1**(-2)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(2)*CF**2 * (  - 64*T1**(-1)*MS2**2 - 64*
     +    U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(3)*N*CF**2 * (  - 32*S*T1**(-2)*MS2**2 - 
     +    32*S*U1**(-2)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(3)*CF**2 * (  - 64*T1**(-1)*MS2**2 - 64*
     +    U1**(-1)*MS2**2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N*CF * ( 32*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N*CF**2 * (  - 64*S*T1**(-1)*MS2 + 32*S
     +    *U1**(-1)*MS2 - 64*S**2*T1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N**2*CF * ( 16*S*T1**(-1)*MS2 - 16*S*
     +    U1**(-1)*MS2 + 16*S**2*T1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*CF**2 * (  - 64*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*N*CF * ( 32*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*N*CF**2 * ( 32*S*T1**(-1)*MS2 - 64*S*
     +    U1**(-1)*MS2 - 64*S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*N**2*CF * (  - 16*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-1)*MS2 + 16*S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*CF**2 * (  - 64*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*N*CF * ( 32*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*N*CF**2 * ( 32*S*T1**(-1)*MS2 - 64*S*
     +    U1**(-1)*MS2 - 64*S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*N**2*CF * (  - 16*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-1)*MS2 + 16*S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*CF**2 * (  - 64*S*U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N*CF * ( 32*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N*CF**2 * (  - 64*S*T1**(-1)*MS2 + 32*S
     +    *U1**(-1)*MS2 - 64*S**2*T1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N**2*CF * ( 16*S*T1**(-1)*MS2 - 16*S*
     +    U1**(-1)*MS2 + 16*S**2*T1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*CF**2 * (  - 64*S*T1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N*CF * ( 64*S*T1**(-1)*MS2 + 64*S*
     +    U1**(-1)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N*CF**2 * (  - 64*S**2*T1**(-2)*MS2 - 
     +    64*S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N**2*CF * ( 32*S**2*T1**(-2)*MS2 + 32*
     +    S**2*U1**(-2)*MS2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*CF**2 * (  - 64*S*T1**(-1)*MS2 - 64*S*
     +    U1**(-1)*MS2 )


      MQQLRV = 0.D0
      MQQLRV = MQQLRV + N*CF * ( 8*S*T**(-1)*T1**(-3)*NS*MS2**3 - 8*S*
     +    T**(-1)*T1**(-3)*MS2**2*MT2 + 12*S*T**(-1)*T1**(-2)*NS*MS2**2
     +     - 12*S*T**(-1)*T1**(-2)*MS2*MT2 + 4*S*T**(-1)*T1**(-1)*NS*
     +    MS2 - 4*S*T**(-1)*T1**(-1)*MT2 + 8*S*U**(-1)*U1**(-3)*NS*
     +    MS2**3 - 8*S*U**(-1)*U1**(-3)*MS2**2*MT2 + 12*S*U**(-1)*
     +    U1**(-2)*NS*MS2**2 - 12*S*U**(-1)*U1**(-2)*MS2*MT2 + 4*S*
     +    U**(-1)*U1**(-1)*NS*MS2 - 4*S*U**(-1)*U1**(-1)*MT2 - 8*S*
     +    T1**(-3)*NS*MS2**2 + 8*S*T1**(-3)*MS2*MT2 - 8*S*T1**(-2)*NS*
     +    MS2 + 8*S*T1**(-2)*MT2 - 8*S*U1**(-3)*NS*MS2**2 + 8*S*
     +    U1**(-3)*MS2*MT2 - 8*S*U1**(-2)*NS*MS2 + 8*S*U1**(-2)*MT2 + 8
     +    *T**(-1)*T1**(-1)*NS*MS2**2 - 8*T**(-1)*T1**(-1)*MS2*MT2 + 4*
     +    T**(-1)*NS*MS2 - 4*T**(-1)*MT2 + 8*U**(-1)*U1**(-1)*NS*MS2**2
     +     - 8*U**(-1)*U1**(-1)*MS2*MT2 + 4*U**(-1)*NS*MS2 - 4*U**(-1)*
     +    MT2 - 8*T1**(-1)*NS*MS2 + 8*T1**(-1)*MT2 - 8*U1**(-1)*NS*MS2
     +     + 8*U1**(-1)*MT2 )
     +
      MQQLRV = MQQLRV + N**2*CF * ( 8 - 8*S*T**(-1)*T1**(-3)*MS2**3 - 
     +    12*S*T**(-1)*T1**(-2)*MS2**2 - 4*S*T**(-1)*T1**(-1)*MS2 - 8*S
     +    *U**(-1)*U1**(-3)*MS2**3 - 12*S*U**(-1)*U1**(-2)*MS2**2 - 4*S
     +    *U**(-1)*U1**(-1)*MS2 + 8*S*T1**(-3)*MS2**2 + 12*S*T1**(-2)*
     +    MS2 + 4*S*T1**(-1) + 8*S*U1**(-3)*MS2**2 + 12*S*U1**(-2)*MS2
     +     + 4*S*U1**(-1) - 8*T**(-1)*T1**(-1)*MS2**2 - 4*T**(-1)*MS2
     +     - 8*U**(-1)*U1**(-1)*MS2**2 - 4*U**(-1)*MS2 + 8*T1**(-1)*MS2
     +     + 8*U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0A(2)*N*CF * (  - 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 - 12*S*T**(-1)*T1**(-2)*MS2*MT2 - 4*S*T**(-1)*
     +    T1**(-1)*MT2 - 8*S*U**(-1)*U1**(-3)*MS2**2*MT2 - 12*S*U**(-1)
     +    *U1**(-2)*MS2*MT2 - 4*S*U**(-1)*U1**(-1)*MT2 + 8*S*T1**(-3)*
     +    MS2*MT2 + 8*S*T1**(-2)*MT2 + 8*S*U1**(-3)*MS2*MT2 + 8*S*
     +    U1**(-2)*MT2 - 8*T**(-1)*T1**(-1)*MS2*MT2 - 4*T**(-1)*MT2 - 8
     +    *U**(-1)*U1**(-1)*MS2*MT2 - 4*U**(-1)*MT2 + 8*T1**(-1)*MT2 + 
     +    8*U1**(-1)*MT2 )
     +
      MQQLRV = MQQLRV + SK1B0C(1)*N*CF**2 * (  - 16 - 32*S*T1**(-3)*
     +    MS2**2 - 40*S*T1**(-2)*MS2 - 8*S*T1**(-1) - 32*S*U1**(-3)*
     +    MS2**2 - 40*S*U1**(-2)*MS2 - 8*S*U1**(-1) - 32*T1**(-1)*MS2
     +     - 32*U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0C(1)*N**2*CF * ( 16*S*T1**(-3)*MS2**2 + 16
     +    *S*T1**(-2)*MS2 + 16*S*U1**(-3)*MS2**2 + 16*S*U1**(-2)*MS2 + 
     +    16*T1**(-1)*MS2 + 16*U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N*CF * (  - 4 - 8*S*T**(-1)*
     +    T1**(-3)*NS*MS2**3 + 8*S*T**(-1)*T1**(-3)*MS2**3 - 12*S*
     +    T**(-1)*T1**(-2)*NS*MS2**2 + 12*S*T**(-1)*T1**(-2)*MS2**2 - 4
     +    *S*T**(-1)*T1**(-1)*NS*MS2 + 4*S*T**(-1)*T1**(-1)*MS2 + 8*S*
     +    T1**(-3)*NS*MS2**2 - 8*S*T1**(-3)*MS2**2 + 12*S*T1**(-2)*NS*
     +    MS2 - 12*S*T1**(-2)*MS2 + 4*S*T1**(-1)*NS - 4*S*T1**(-1) - 8*
     +    T**(-1)*T1**(-1)*NS*MS2**2 + 8*T**(-1)*T1**(-1)*MS2**2 - 4*
     +    T**(-1)*NS*MS2 + 4*T**(-1)*MS2 + 8*T1**(-1)*NS*MS2 - 8*
     +    T1**(-1)*MS2 + 4*NS )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N*CF**2 * ( 16 + 32*S*T1**(-3)*
     +    MS2**2 + 48*S*T1**(-2)*MS2 + 16*S*T1**(-1) + 32*T1**(-1)*MS2
     +     )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N**2*CF * (  - 12 + 8*S*T**(-1)*
     +    T1**(-3)*MS2**3 + 12*S*T**(-1)*T1**(-2)*MS2**2 + 4*S*T**(-1)*
     +    T1**(-1)*MS2 - 24*S*T1**(-3)*MS2**2 - 36*S*T1**(-2)*MS2 - 12*
     +    S*T1**(-1) + 8*T**(-1)*T1**(-1)*MS2**2 + 4*T**(-1)*MS2 - 24*
     +    T1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N*CF * (  - 4 - 8*S*U**(-1)*
     +    U1**(-3)*NS*MS2**3 + 8*S*U**(-1)*U1**(-3)*MS2**3 - 12*S*
     +    U**(-1)*U1**(-2)*NS*MS2**2 + 12*S*U**(-1)*U1**(-2)*MS2**2 - 4
     +    *S*U**(-1)*U1**(-1)*NS*MS2 + 4*S*U**(-1)*U1**(-1)*MS2 + 8*S*
     +    U1**(-3)*NS*MS2**2 - 8*S*U1**(-3)*MS2**2 + 12*S*U1**(-2)*NS*
     +    MS2 - 12*S*U1**(-2)*MS2 + 4*S*U1**(-1)*NS - 4*S*U1**(-1) - 8*
     +    U**(-1)*U1**(-1)*NS*MS2**2 + 8*U**(-1)*U1**(-1)*MS2**2 - 4*
     +    U**(-1)*NS*MS2 + 4*U**(-1)*MS2 + 8*U1**(-1)*NS*MS2 - 8*
     +    U1**(-1)*MS2 + 4*NS )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N*CF**2 * ( 16 + 32*S*U1**(-3)*
     +    MS2**2 + 48*S*U1**(-2)*MS2 + 16*S*U1**(-1) + 32*U1**(-1)*MS2
     +     )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N**2*CF * (  - 12 + 8*S*U**(-1)*
     +    U1**(-3)*MS2**3 + 12*S*U**(-1)*U1**(-2)*MS2**2 + 4*S*U**(-1)*
     +    U1**(-1)*MS2 - 24*S*U1**(-3)*MS2**2 - 36*S*U1**(-2)*MS2 - 12*
     +    S*U1**(-1) + 8*U**(-1)*U1**(-1)*MS2**2 + 4*U**(-1)*MS2 - 24*
     +    U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0D(3,1)*N*CF * ( 4 + 8*S*T**(-1)*T1**(-3)*
     +    MS2**2*MT2 - 8*S*T**(-1)*T1**(-3)*MS2**3 + 12*S*T**(-1)*
     +    T1**(-2)*MS2*MT2 - 12*S*T**(-1)*T1**(-2)*MS2**2 - 4*S*T**(-1)
     +    *T1**(-1)*MS2 + 4*S*T**(-1)*T1**(-1)*MT2 + 8*S*T1**(-3)*
     +    MS2**2 + 12*S*T1**(-2)*MS2 + 4*S*T1**(-1) + 8*T**(-1)*
     +    T1**(-1)*MS2*MT2 - 8*T**(-1)*T1**(-1)*MS2**2 - 4*T**(-1)*MS2
     +     + 4*T**(-1)*MT2 + 8*T1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0D(3,2)*N*CF * ( 4 + 8*S*U**(-1)*U1**(-3)*
     +    MS2**2*MT2 - 8*S*U**(-1)*U1**(-3)*MS2**3 + 12*S*U**(-1)*
     +    U1**(-2)*MS2*MT2 - 12*S*U**(-1)*U1**(-2)*MS2**2 - 4*S*U**(-1)
     +    *U1**(-1)*MS2 + 4*S*U**(-1)*U1**(-1)*MT2 + 8*S*U1**(-3)*
     +    MS2**2 + 12*S*U1**(-2)*MS2 + 4*S*U1**(-1) + 8*U**(-1)*
     +    U1**(-1)*MS2*MT2 - 8*U**(-1)*U1**(-1)*MS2**2 - 4*U**(-1)*MS2
     +     + 4*U**(-1)*MT2 + 8*U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0E(3)*N*CF * (  - 8*S*T1**(-3)*MS2*MT2 - 8*
     +    S*T1**(-2)*MT2 - 8*S*U1**(-3)*MS2*MT2 - 8*S*U1**(-2)*MT2 - 8*
     +    T1**(-1)*MT2 - 8*U1**(-1)*MT2 )
     +
      MQQLRV = MQQLRV + SK1BP(1)*N*CF**2 * (  - 16*S*T1**(-2)*MS2**2 - 
     +    16*S*T1**(-1)*MS2 - 16*S*U1**(-2)*MS2**2 - 16*S*U1**(-1)*MS2
     +     - 32*MS2 )
     +
      MQQLRV = MQQLRV + SK1C0A(1)*N*CF**2 * (  - 16*S - 8*S**2*T1**(-1)
     +     - 8*S**2*U1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0A(1)*N**2*CF * ( 8*S + 4*S**2*T1**(-1) + 4
     +    *S**2*U1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0A(2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 16
     +    *S*U1**(-1)*MS2 + 8*S**2*T1**(-1) + 8*S**2*U1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0A(2)*N**2*CF * ( 8*S*T1**(-1)*MS2 + 8*S*
     +    U1**(-1)*MS2 - 4*S**2*T1**(-1) - 4*S**2*U1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0B(1)*N*CF**2 * ( 16*S*T1**(-1)*MS2 + 16*S*
     +    U1**(-1)*MS2 - 32*S - 8*S**2*T1**(-1) - 8*S**2*U1**(-1) + 64*
     +    MS2 )
     +
      MQQLRV = MQQLRV + SK1C0B(1)*N**2*CF * (  - 8*S*T1**(-1)*MS2 - 8*S
     +    *U1**(-1)*MS2 + 16*S + 4*S**2*T1**(-1) + 4*S**2*U1**(-1) - 32
     +    *MS2 )
     +
      MQQLRV = MQQLRV + SK1C0B(3)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 - 16
     +    *S*U1**(-1)*MS2 - 64*S*S1**(-1)*MS2 + 32*S + 8*S**2*T1**(-1)
     +     + 8*S**2*U1**(-1) + 16*S**2*S1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0B(3)*N**2*CF * ( 8*S*T1**(-1)*MS2 + 8*S*
     +    U1**(-1)*MS2 + 32*S*S1**(-1)*MS2 - 16*S - 4*S**2*T1**(-1) - 4
     +    *S**2*U1**(-1) - 8*S**2*S1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0C(1,1)*N*CF**2 * ( 16*S*T1**(-2)*MS2**2 + 
     +    16*S*T1**(-1)*MS2 - 16*S*U1**(-1)*MS2 + 16*S + 16*S**2*
     +    U1**(-1) - 32*T1 - 16*MS2 )
     +
      MQQLRV = MQQLRV + SK1C0C(1,1)*N**2*CF * (  - 16*S*T1**(-2)*MS2**2
     +     - 16*S*T1**(-1)*MS2 + 4*S*U1**(-1)*MS2 - 8*S - 4*S**2*
     +    U1**(-1) + 8*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(1,2)*N*CF**2 * (  - 16*S*T1**(-1)*MS2 + 
     +    16*S*U1**(-2)*MS2**2 + 16*S*U1**(-1)*MS2 + 48*S + 16*S**2*
     +    T1**(-1) + 32*T1 - 16*MS2 )
     +
      MQQLRV = MQQLRV + SK1C0C(1,2)*N**2*CF * ( 4*S*T1**(-1)*MS2 - 16*S
     +    *U1**(-2)*MS2**2 - 16*S*U1**(-1)*MS2 - 16*S - 4*S**2*T1**(-1)
     +     - 8*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,1)*N*CF**2 * ( 16*S*T1**(-1)*MS2 + 16*
     +    S*U1**(-1)*MS2 - 16*S**2*U1**(-1) + 32*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,1)*N**2*CF * (  - 16*S*T1**(-1)*MS2 - 
     +    4*S*U1**(-1)*MS2 - 8*S + 4*S**2*U1**(-1) - 16*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,2)*N*CF**2 * ( 16*S*T1**(-1)*MS2 + 16*
     +    S*U1**(-1)*MS2 - 32*S - 16*S**2*T1**(-1) - 32*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,2)*N**2*CF * (  - 4*S*T1**(-1)*MS2 - 
     +    16*S*U1**(-1)*MS2 + 8*S + 4*S**2*T1**(-1) + 16*T1 )
     +
      MQQLRV = MQQLRV + SK1D0(1,1)*N*CF**2 * (  - 8*S*T1 - 16*S**2*
     +    T1**(-1)*MS2 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(1,1)*N**2*CF * ( 4*S*T1 + 8*S**2*T1**(-1)
     +    *MS2 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(1,2)*N*CF**2 * ( 8*S*T1 - 16*S**2*
     +    U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1D0(1,2)*N**2*CF * (  - 4*S*T1 + 8*S**2*
     +    U1**(-1)*MS2 )
     +
      MQQLRV = MQQLRV + SK1D0(2,1)*N*CF**2 * ( 32*S*T1**(-1)*MS2**2 + 
     +    16*S*MS2 - 16*S**2*T1**(-1)*MS2 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(2,1)*N**2*CF * (  - 16*S*T1**(-1)*MS2**2
     +     - 8*S*MS2 + 8*S**2*T1**(-1)*MS2 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(2,2)*N*CF**2 * ( 32*S*U1**(-1)*MS2**2 + 
     +    16*S*MS2 - 16*S**2*U1**(-1)*MS2 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(2,2)*N**2*CF * (  - 16*S*U1**(-1)*MS2**2
     +     - 8*S*MS2 + 8*S**2*U1**(-1)*MS2 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(3,1)*N*CF**2 * (  - 32*S*T1 - 48*S*MS2 - 
     +    32*S**2*T1**(-1)*MS2 - 16*S**2 - 16*T1*MS2 - 16*T1**2 )
     +
      MQQLRV = MQQLRV + SK1D0(3,1)*N**2*CF * ( 8*S*T1 + 12*S*MS2 + 8*
     +    S**2*T1**(-1)*MS2 + 4*S**2 + 4*T1*MS2 + 4*T1**2 )
     +
      MQQLRV = MQQLRV + SK1D0(3,2)*N*CF**2 * (  - 32*S*MS2 - 32*S**2*
     +    U1**(-1)*MS2 + 16*T1*MS2 - 16*T1**2 )
     +
      MQQLRV = MQQLRV + SK1D0(3,2)*N**2*CF * ( 8*S*MS2 + 8*S**2*
     +    U1**(-1)*MS2 - 4*T1*MS2 + 4*T1**2 )
     +
      MQQLRV = MQQLRV + SOF1(1)*N*CF**2 * (  - 128*S*T1**(-2)*MS2**2 - 
     +    128*S*T1**(-1)*MS2 - 128*S*U1**(-2)*MS2**2 - 128*S*U1**(-1)*
     +    MS2 + 128*S + 64*S**2*T1**(-2)*MS2 + 64*S**2*T1**(-1) + 64*
     +    S**2*U1**(-2)*MS2 + 64*S**2*U1**(-1) - 256*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(1)*N**2*CF * ( 64*S*T1**(-2)*MS2**2 + 64*S
     +    *T1**(-1)*MS2 + 64*S*U1**(-2)*MS2**2 + 64*S*U1**(-1)*MS2 - 64
     +    *S - 32*S**2*T1**(-2)*MS2 - 32*S**2*T1**(-1) - 32*S**2*
     +    U1**(-2)*MS2 - 32*S**2*U1**(-1) + 128*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(2)*N*CF**2 * ( 32*S*T1**(-2)*MS2**2 + 32*S
     +    *T1**(-1)*MS2 + 32*S*U1**(-2)*MS2**2 + 32*S*U1**(-1)*MS2 + 64
     +    *MS2 )
     +
      MQQLRV = MQQLRV + SOF1(3)*N*CF**2 * ( 32*S*T1**(-2)*MS2**2 + 32*S
     +    *T1**(-1)*MS2 + 32*S*U1**(-2)*MS2**2 + 32*S*U1**(-1)*MS2 + 64
     +    *MS2 )
     +
      MQQLRV = MQQLRV + SOF1(4)*N*CF**2 * ( 64*S*T1**(-1)*MS2 - 32*S*
     +    U1**(-1)*MS2 + 128*S + 64*S**2*T1**(-2)*MS2 + 64*S**2*
     +    T1**(-1) + 96*T1 )
     +
      MQQLRV = MQQLRV + SOF1(4)*N**2*CF * (  - 16*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-1)*MS2 - 32*S - 16*S**2*T1**(-2)*MS2 - 16*S**2*
     +    T1**(-1) - 32*T1 )
     +
      MQQLRV = MQQLRV + SOF1(5)*N*CF**2 * (  - 32*S*T1**(-1)*MS2 + 64*S
     +    *U1**(-1)*MS2 + 32*S + 64*S**2*U1**(-2)*MS2 + 64*S**2*
     +    U1**(-1) - 96*T1 )
     +
      MQQLRV = MQQLRV + SOF1(5)*N**2*CF * ( 16*S*T1**(-1)*MS2 - 16*S*
     +    U1**(-1)*MS2 - 16*S**2*U1**(-2)*MS2 - 16*S**2*U1**(-1) + 32*
     +    T1 )
     +
      MQQLRV = MQQLRV + SOF1(6)*N*CF**2 * (  - 32*S*T1**(-1)*MS2 + 64*S
     +    *U1**(-1)*MS2 + 32*S + 64*S**2*U1**(-2)*MS2 + 64*S**2*
     +    U1**(-1) - 96*T1 )
     +
      MQQLRV = MQQLRV + SOF1(6)*N**2*CF * ( 16*S*T1**(-1)*MS2 - 16*S*
     +    U1**(-1)*MS2 - 16*S**2*U1**(-2)*MS2 - 16*S**2*U1**(-1) + 32*
     +    T1 )
     +
      MQQLRV = MQQLRV + SOF1(7)*N*CF**2 * ( 64*S*T1**(-1)*MS2 - 32*S*
     +    U1**(-1)*MS2 + 128*S + 64*S**2*T1**(-2)*MS2 + 64*S**2*
     +    T1**(-1) + 96*T1 )
     +
      MQQLRV = MQQLRV + SOF1(7)*N**2*CF * (  - 16*S*T1**(-1)*MS2 + 16*S
     +    *U1**(-1)*MS2 - 32*S - 16*S**2*T1**(-2)*MS2 - 16*S**2*
     +    T1**(-1) - 32*T1 )
     +
      MQQLRV = MQQLRV + SOF1(8)*N*CF**2 * ( 128*S + 64*S**2*T1**(-2)*
     +    MS2 + 64*S**2*T1**(-1) + 64*S**2*U1**(-2)*MS2 + 64*S**2*
     +    U1**(-1) )
     +
      MQQLRV = MQQLRV + SOF1(8)*N**2*CF * (  - 64*S - 32*S**2*T1**(-2)*
     +    MS2 - 32*S**2*T1**(-1) - 32*S**2*U1**(-2)*MS2 - 32*S**2*
     +    U1**(-1) )

         END IF

      ELSE

         IF (IFL.EQ.0) THEN
      MQPLLV = 0.D0
      MQPLLV = MQPLLV + N*CF * (  - 8*TG**(-3)*S*T**(-1)*NS*MS2*MG2**2
     +     + 8*TG**(-3)*S*T**(-1)*MG2**2*MT2 + 8*TG**(-3)*S*NS*MS2*MG2
     +     - 8*TG**(-3)*S*MG2*MT2 - 8*TG**(-2)*S*T**(-1)*NS*MS2*MG2 + 8
     +    *TG**(-2)*S*T**(-1)*MG2*MT2 + 4*TG**(-2)*S*NS*MS2 - 4*
     +    TG**(-2)*S*MT2 )
     +
      MQPLLV = MQPLLV + N**2*CF * ( 8*TG**(-3)*S*T**(-1)*MG2**3 - 8*
     +    TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*MG2**2 - 8*TG**(-2)*
     +    S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0A(1)*N**2*CF * ( 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 - 8*TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*MG2**2 - 
     +    4*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0A(2)*N*CF * ( 8*TG**(-3)*S*T**(-1)*MG2**2*
     +    MT2 - 8*TG**(-3)*S*MG2*MT2 + 8*TG**(-2)*S*T**(-1)*MG2*MT2 - 4
     +    *TG**(-2)*S*MT2 )
     +
      MQPLLV = MQPLLV + SK1B0A(3)*N*CF**2 * (  - 4*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0C(1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 16*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0C(2)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 8*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N*CF * ( 8*TG**(-3)*S*T**(-1)*NS*
     +    MS2*MG2**2 - 8*TG**(-3)*S*T**(-1)*MS2*MG2**2 - 8*TG**(-3)*S*
     +    NS*MG2**2 + 8*TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*NS*MS2
     +    *MG2 - 8*TG**(-2)*S*T**(-1)*MS2*MG2 - 8*TG**(-2)*S*NS*MG2 + 8
     +    *TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N*CF**2 * (  - 32*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2 - 32*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0D(1,1)*N**2*CF * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 16*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0D(2,1)*N**2*CF * (  - 8*TG**(-3)*S*T**(-1)
     +    *MG2**3 + 24*TG**(-3)*S*MG2**2 - 8*TG**(-2)*S*T**(-1)*MG2**2
     +     - 16*TG**(-2)*S*T1**(-1)*MS2*MG2 + 8*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0D(3,1)*N*CF * ( 8*TG**(-3)*S*T**(-1)*MS2*
     +    MG2**2 - 8*TG**(-3)*S*T**(-1)*MG2**2*MT2 - 8*TG**(-3)*S*
     +    MG2**2 + 8*TG**(-2)*S*T**(-1)*MS2*MG2 - 8*TG**(-2)*S*T**(-1)*
     +    MG2*MT2 - 8*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0E(1)*N**2*CF * (  - 16*TG**(-3)*S*MG2**2
     +     - 8*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0E(2)*N*CF * (  - 8*TG**(-3)*S*NS*MS2*MG2
     +     + 8*TG**(-3)*S*NS*MG2**2 + 8*TG**(-3)*S*MS2*MG2 - 8*TG**(-3)
     +    *S*MG2**2 - 4*TG**(-2)*S*NS*MS2 + 4*TG**(-2)*S*NS*MG2 + 4*
     +    TG**(-2)*S*MS2 - 4*TG**(-2)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1B0E(3)*N*CF * (  - 8*TG**(-3)*S*MS2*MG2 + 8*
     +    TG**(-3)*S*MG2*MT2 + 8*TG**(-3)*S*MG2**2 - 4*TG**(-2)*S*MS2
     +     + 4*TG**(-2)*S*MG2 + 4*TG**(-2)*S*MT2 )
     +
      MQPLLV = MQPLLV + SK1BP(1)*N*CF**2 * ( 16*TG**(-2)*S*MS2*MG2 )
     +
      MQPLLV = MQPLLV + SK1BP(2)*N*CF**2 * (  - 8*TG**(-2)*S*MS2*MG2 + 
     +    8*TG**(-2)*S*MG2**2 )
     +
      MQPLLV = MQPLLV + SK1BP(3)*N*CF**2 * ( 4*TG**(-2)*S*MS2*MG2 - 4*
     +    TG**(-2)*S*MG2**2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,1)*N*CF**2 * ( 16*TG**(-2)*S*MS2*MG2
     +     - 16*TG**(-2)*S*MG2**2 - 16*TG**(-1)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,1)*N**2*CF * (  - 8*TG**(-2)*S*MS2*MG2
     +     + 8*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,2)*N*CF**2 * (  - 16*TG**(-1)*S*MG2 + 
     +    16*TG**(-1)*MS2*MG2 - 16*TG**(-1)*MG2**2 - 16*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(3,2)*N**2*CF * ( 4*TG**(-1)*S*MG2 - 4*
     +    TG**(-1)*MS2*MG2 + 4*TG**(-1)*MG2**2 + 4*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(4,1)*N*CF**2 * (  - 16*TG**(-1)*MS2*MG2
     +     + 16*TG**(-1)*MG2**2 + 16*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(4,1)*N**2*CF * (  - 8*TG**(-2)*S*MS2*MG2
     +     + 8*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 + 4*TG**(-1)*MS2*
     +    MG2 - 4*TG**(-1)*MG2**2 - 4*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(5,1)*N*CF**2 * (  - 16*TG**(-1)*MS2*MG2
     +     + 16*TG**(-1)*MG2**2 + 16*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(5,1)*N**2*CF * ( 8*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 8*TG**(-2)*S*
     +    MG2**2 + 8*TG**(-1)*S*MG2 + 4*TG**(-1)*MS2*MG2 - 4*TG**(-1)*
     +    MG2**2 - 4*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(6,1)*N*CF**2 * (  - 16*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2**2 + 16*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 16*
     +    TG**(-2)*S*MS2*MG2 - 32*TG**(-2)*S*MG2**2 - 16*TG**(-1)*S*MG2
     +     )
     +
      MQPLLV = MQPLLV + SK1C0C(6,1)*N**2*CF * ( 8*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 - 8*TG**(-2)*S*
     +    MS2*MG2 + 16*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(6,2)*N*CF**2 * (  - 16*TG**(-1)*S*MG2 + 
     +    16*TG**(-1)*MS2*MG2 - 16*TG**(-1)*MG2**2 - 16*MG2 )
     +
      MQPLLV = MQPLLV + SK1C0C(6,2)*N**2*CF * ( 4*TG**(-1)*S*MG2 - 4*
     +    TG**(-1)*MS2*MG2 + 4*TG**(-1)*MG2**2 + 4*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(4,1)*N*CF**2 * ( 16*TG**(-1)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(4,1)*N**2*CF * (  - 8*TG**(-1)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(7,1)*N*CF**2 * (  - 32*TG**(-1)*S*MS2*MG2
     +     + 16*TG**(-1)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(7,1)*N**2*CF * ( 16*TG**(-1)*S*MS2*MG2 - 
     +    8*TG**(-1)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(8,1)*N*CF**2 * (  - 32*TG**(-1)*S*MS2*MG2
     +     + 32*TG**(-1)*S*MG2**2 + 32*TG**(-1)*S**2*MG2 + 48*S*MG2 + 
     +    16*T1*MG2 )
     +
      MQPLLV = MQPLLV + SK1D0(8,1)*N**2*CF * ( 8*TG**(-1)*S*MS2*MG2 - 8
     +    *TG**(-1)*S*MG2**2 - 8*TG**(-1)*S**2*MG2 - 12*S*MG2 - 4*T1*
     +    MG2 )
     +
      MQPLLV = MQPLLV + SOF1(1)*N*CF**2 * ( 128*TG**(-2)*S*MS2*MG2 - 64
     +    *TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(1)*N**2*CF * (  - 64*TG**(-2)*S*MS2*MG2 + 
     +    32*TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(2)*N*CF**2 * (  - 32*TG**(-2)*S*MS2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(3)*N*CF**2 * (  - 32*TG**(-2)*S*MS2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(4)*N*CF**2 * (  - 64*TG**(-2)*S*T1*MG2 - 
     +    64*TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(4)*N**2*CF * ( 16*TG**(-2)*S*T1*MG2 + 16*
     +    TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(5)*N*CF**2 * ( 32*TG**(-2)*S*T1*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(5)*N**2*CF * (  - 16*TG**(-2)*S*T1*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(6)*N*CF**2 * ( 32*TG**(-2)*S*T1*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(6)*N**2*CF * (  - 16*TG**(-2)*S*T1*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(7)*N*CF**2 * (  - 64*TG**(-2)*S*T1*MG2 - 
     +    64*TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(7)*N**2*CF * ( 16*TG**(-2)*S*T1*MG2 + 16*
     +    TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(8)*N*CF**2 * (  - 64*TG**(-2)*S**2*MG2 )
     +
      MQPLLV = MQPLLV + SOF1(8)*N**2*CF * ( 32*TG**(-2)*S**2*MG2 )

      MQPLRV = 0.D0
      MQPLRV = MQPLRV + N*CF * ( 8*TG**(-3)*S*T**(-1)*NS*MS2*MG2**2 - 8
     +    *TG**(-3)*S*T**(-1)*MG2**2*MT2 - 8*TG**(-3)*S*NS*MS2*MG2 + 8*
     +    TG**(-3)*S*MG2*MT2 + 8*TG**(-3)*T**(-1)*NS*MS2*MG2**3 - 16*
     +    TG**(-3)*T**(-1)*NS*MS2**2*MG2**2 + 8*TG**(-3)*T**(-1)*NS*
     +    MS2**3*MG2 + 16*TG**(-3)*T**(-1)*MS2*MG2**2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MS2**2*MG2*MT2 - 8*TG**(-3)*T**(-1)*MG2**3*MT2 - 8*
     +    TG**(-3)*NS*MS2*MG2**2 + 16*TG**(-3)*NS*MS2**2*MG2 - 8*
     +    TG**(-3)*NS*MS2**3 - 16*TG**(-3)*MS2*MG2*MT2 + 8*TG**(-3)*
     +    MS2**2*MT2 + 8*TG**(-3)*MG2**2*MT2 + 12*TG**(-2)*S*T**(-1)*NS
     +    *MS2*MG2 - 12*TG**(-2)*S*T**(-1)*MG2*MT2 - 8*TG**(-2)*S*NS*
     +    MS2 + 8*TG**(-2)*S*MT2 + 12*TG**(-2)*T**(-1)*T1*NS*MS2*MG2 - 
     +    4*TG**(-2)*T**(-1)*T1*NS*MS2**2 + 4*TG**(-2)*T**(-1)*T1*MS2*
     +    MT2 - 12*TG**(-2)*T**(-1)*T1*MG2*MT2 + 8*TG**(-2)*T**(-1)*NS*
     +    MS2*MG2**2 - 8*TG**(-2)*T**(-1)*NS*MS2**2*MG2 + 8*TG**(-2)*
     +    T**(-1)*MS2*MG2*MT2 - 8*TG**(-2)*T**(-1)*MG2**2*MT2 - 8*
     +    TG**(-2)*T1*NS*MS2 )
     +
      MQPLRV = MQPLRV + N*CF * ( 8*TG**(-2)*T1*MT2 - 8*TG**(-2)*NS*MS2*
     +    MG2 + 8*TG**(-2)*NS*MS2**2 - 8*TG**(-2)*MS2*MT2 + 8*TG**(-2)*
     +    MG2*MT2 + 4*TG**(-1)*S*T**(-1)*NS*MS2 - 4*TG**(-1)*S*T**(-1)*
     +    MT2 + 4*TG**(-1)*T**(-1)*T1*NS*MS2 - 4*TG**(-1)*T**(-1)*T1*
     +    MT2 )
     +
      MQPLRV = MQPLRV + N**2*CF * (  - 8*TG**(-3)*S*T**(-1)*MG2**3 + 8*
     +    TG**(-3)*S*MG2**2 + 16*TG**(-3)*T**(-1)*MS2*MG2**3 - 8*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*MG2**4 - 
     +    16*TG**(-3)*MS2*MG2**2 + 8*TG**(-3)*MS2**2*MG2 + 8*TG**(-3)*
     +    MG2**3 - 12*TG**(-2)*S*T**(-1)*MG2**2 + 12*TG**(-2)*S*MG2 + 4
     +    *TG**(-2)*T**(-1)*T1*MS2*MG2 - 12*TG**(-2)*T**(-1)*T1*MG2**2
     +     + 8*TG**(-2)*T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*MG2**3
     +     - 4*TG**(-2)*T1*MS2 + 12*TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*
     +    MG2 + 8*TG**(-2)*MG2**2 - 4*TG**(-1)*S*T**(-1)*MG2 + 4*
     +    TG**(-1)*S - 4*TG**(-1)*T**(-1)*T1*MG2 + 4*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0A(1)*N**2*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 + 8*TG**(-3)*S*MG2**2 + 16*TG**(-3)*T**(-1)*MS2*MG2**3
     +     - 8*TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*
     +    MG2**4 - 16*TG**(-3)*MS2*MG2**2 + 8*TG**(-3)*MS2**2*MG2 + 8*
     +    TG**(-3)*MG2**3 - 12*TG**(-2)*S*T**(-1)*MG2**2 + 8*TG**(-2)*S
     +    *MG2 + 4*TG**(-2)*T**(-1)*T1*MS2*MG2 - 12*TG**(-2)*T**(-1)*T1
     +    *MG2**2 + 8*TG**(-2)*T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*
     +    MG2**3 + 8*TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*MG2 + 8*TG**(-2)*
     +    MG2**2 - 4*TG**(-1)*S*T**(-1)*MG2 - 4*TG**(-1)*T**(-1)*T1*MG2
     +     )
     +
      MQPLRV = MQPLRV + SK1B0A(2)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MG2**2*MT2 + 8*TG**(-3)*S*MG2*MT2 + 16*TG**(-3)*T**(-1)*MS2*
     +    MG2**2*MT2 - 8*TG**(-3)*T**(-1)*MS2**2*MG2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MG2**3*MT2 - 16*TG**(-3)*MS2*MG2*MT2 + 8*TG**(-3)*
     +    MS2**2*MT2 + 8*TG**(-3)*MG2**2*MT2 - 12*TG**(-2)*S*T**(-1)*
     +    MG2*MT2 + 8*TG**(-2)*S*MT2 + 4*TG**(-2)*T**(-1)*T1*MS2*MT2 - 
     +    12*TG**(-2)*T**(-1)*T1*MG2*MT2 + 8*TG**(-2)*T**(-1)*MS2*MG2*
     +    MT2 - 8*TG**(-2)*T**(-1)*MG2**2*MT2 + 8*TG**(-2)*T1*MT2 - 8*
     +    TG**(-2)*MS2*MT2 + 8*TG**(-2)*MG2*MT2 - 4*TG**(-1)*S*T**(-1)*
     +    MT2 - 4*TG**(-1)*T**(-1)*T1*MT2 )
     +
      MQPLRV = MQPLRV + SK1B0A(3)*N*CF**2 * ( 4*TG**(-2)*S*MG2 - 4*
     +    TG**(-2)*T1*MS2 + 4*TG**(-2)*T1*MG2 + 4*TG**(-1)*S + 4*
     +    TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0C(1)*N*CF**2 * (  - 16*TG**(-2)*S*T1**(-1)
     +    *MS2**2 - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*S*MG2 - 16*TG**(-2)
     +    *T1*MG2 - 16*TG**(-1)*S - 16*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0C(2)*N*CF**2 * (  - 16*TG**(-2)*S*T1**(-1)
     +    *MS2*MG2 - 8*TG**(-2)*S*MG2 - 8*TG**(-2)*T1*MS2 - 8*TG**(-2)*
     +    T1*MG2 + 8*TG**(-1)*S + 8*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*NS
     +    *MS2*MG2**2 + 8*TG**(-3)*S*T**(-1)*MS2*MG2**2 + 8*TG**(-3)*S*
     +    NS*MG2**2 - 8*TG**(-3)*S*MG2**2 - 8*TG**(-3)*T**(-1)*NS*MS2*
     +    MG2**3 + 16*TG**(-3)*T**(-1)*NS*MS2**2*MG2**2 - 8*TG**(-3)*
     +    T**(-1)*NS*MS2**3*MG2 + 8*TG**(-3)*T**(-1)*MS2*MG2**3 - 16*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 + 8*TG**(-3)*T**(-1)*MS2**3*
     +    MG2 - 16*TG**(-3)*NS*MS2*MG2**2 + 8*TG**(-3)*NS*MS2**2*MG2 + 
     +    8*TG**(-3)*NS*MG2**3 + 16*TG**(-3)*MS2*MG2**2 - 8*TG**(-3)*
     +    MS2**2*MG2 - 8*TG**(-3)*MG2**3 - 12*TG**(-2)*S*T**(-1)*NS*MS2
     +    *MG2 + 12*TG**(-2)*S*T**(-1)*MS2*MG2 + 12*TG**(-2)*S*NS*MG2
     +     - 12*TG**(-2)*S*MG2 - 12*TG**(-2)*T**(-1)*T1*NS*MS2*MG2 + 4*
     +    TG**(-2)*T**(-1)*T1*NS*MS2**2 + 12*TG**(-2)*T**(-1)*T1*MS2*
     +    MG2 - 4*TG**(-2)*T**(-1)*T1*MS2**2 - 8*TG**(-2)*T**(-1)*NS*
     +    MS2*MG2**2 + 8*TG**(-2)*T**(-1)*NS*MS2**2*MG2 + 8*TG**(-2)*
     +    T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*MS2**2*MG2 - 4*
     +    TG**(-2)*T1*NS*MS2 )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N*CF * ( 12*TG**(-2)*T1*NS*MG2 + 4*
     +    TG**(-2)*T1*MS2 - 12*TG**(-2)*T1*MG2 - 8*TG**(-2)*NS*MS2*MG2
     +     + 8*TG**(-2)*NS*MG2**2 + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*
     +    MG2**2 - 4*TG**(-1)*S*T**(-1)*NS*MS2 + 4*TG**(-1)*S*T**(-1)*
     +    MS2 + 4*TG**(-1)*S*NS - 4*TG**(-1)*S - 4*TG**(-1)*T**(-1)*T1*
     +    NS*MS2 + 4*TG**(-1)*T**(-1)*T1*MS2 + 4*TG**(-1)*T1*NS - 4*
     +    TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 16*TG**(-2)*S*T1**(-1)*MS2**2 + 16*TG**(-2)*S*MS2
     +     + 32*TG**(-2)*S*MG2 + 32*TG**(-2)*T1*MG2 + 16*TG**(-1)*S + 
     +    16*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0D(1,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2 - 8*TG**(-2)*S*T1**(-1)*MS2**2 - 8*TG**(-2)*
     +    S*MS2 - 16*TG**(-2)*S*MG2 - 16*TG**(-2)*T1*MG2 - 8*TG**(-1)*S
     +     - 8*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0D(2,1)*N**2*CF * ( 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 - 24*TG**(-3)*S*MG2**2 - 16*TG**(-3)*T**(-1)*MS2*
     +    MG2**3 + 8*TG**(-3)*T**(-1)*MS2**2*MG2**2 + 8*TG**(-3)*
     +    T**(-1)*MG2**4 + 48*TG**(-3)*MS2*MG2**2 - 24*TG**(-3)*MS2**2*
     +    MG2 - 24*TG**(-3)*MG2**3 + 12*TG**(-2)*S*T**(-1)*MG2**2 + 8*
     +    TG**(-2)*S*T1**(-1)*MS2*MG2 + 8*TG**(-2)*S*T1**(-1)*MS2**2 + 
     +    8*TG**(-2)*S*MS2 - 20*TG**(-2)*S*MG2 - 4*TG**(-2)*T**(-1)*T1*
     +    MS2*MG2 + 12*TG**(-2)*T**(-1)*T1*MG2**2 - 8*TG**(-2)*T**(-1)*
     +    MS2*MG2**2 + 8*TG**(-2)*T**(-1)*MG2**3 + 12*TG**(-2)*T1*MS2
     +     - 20*TG**(-2)*T1*MG2 + 24*TG**(-2)*MS2*MG2 - 24*TG**(-2)*
     +    MG2**2 + 4*TG**(-1)*S*T**(-1)*MG2 - 4*TG**(-1)*S + 4*TG**(-1)
     +    *T**(-1)*T1*MG2 - 4*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0D(3,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MS2*MG2**2 + 8*TG**(-3)*S*T**(-1)*MG2**2*MT2 + 8*TG**(-3)*S*
     +    MG2**2 - 16*TG**(-3)*T**(-1)*MS2*MG2**2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MS2*MG2**3 + 8*TG**(-3)*T**(-1)*MS2**2*MG2*MT2 + 16*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*MS2**3*
     +    MG2 + 8*TG**(-3)*T**(-1)*MG2**3*MT2 - 16*TG**(-3)*MS2*MG2**2
     +     + 8*TG**(-3)*MS2**2*MG2 + 8*TG**(-3)*MG2**3 - 12*TG**(-2)*S*
     +    T**(-1)*MS2*MG2 + 12*TG**(-2)*S*T**(-1)*MG2*MT2 + 12*TG**(-2)
     +    *S*MG2 - 12*TG**(-2)*T**(-1)*T1*MS2*MG2 - 4*TG**(-2)*T**(-1)*
     +    T1*MS2*MT2 + 4*TG**(-2)*T**(-1)*T1*MS2**2 + 12*TG**(-2)*
     +    T**(-1)*T1*MG2*MT2 - 8*TG**(-2)*T**(-1)*MS2*MG2*MT2 - 8*
     +    TG**(-2)*T**(-1)*MS2*MG2**2 + 8*TG**(-2)*T**(-1)*MS2**2*MG2
     +     + 8*TG**(-2)*T**(-1)*MG2**2*MT2 - 4*TG**(-2)*T1*MS2 + 12*
     +    TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*MG2 + 8*TG**(-2)*MG2**2 - 4*
     +    TG**(-1)*S*T**(-1)*MS2 + 4*TG**(-1)*S*T**(-1)*MT2 + 4*
     +    TG**(-1)*S )
     +
      MQPLRV = MQPLRV + SK1B0D(3,1)*N*CF * (  - 4*TG**(-1)*T**(-1)*T1*
     +    MS2 + 4*TG**(-1)*T**(-1)*T1*MT2 + 4*TG**(-1)*T1 )
     +
      MQPLRV = MQPLRV + SK1B0E(1)*N**2*CF * ( 16*TG**(-3)*S*MG2**2 - 32
     +    *TG**(-3)*MS2*MG2**2 + 16*TG**(-3)*MS2**2*MG2 + 16*TG**(-3)*
     +    MG2**3 + 16*TG**(-2)*S*MG2 + 16*TG**(-2)*T1*MG2 - 16*TG**(-2)
     +    *MS2*MG2 + 16*TG**(-2)*MG2**2 )
     +
      MQPLRV = MQPLRV + SK1B0E(2)*N*CF * ( 8*TG**(-3)*S*NS*MS2*MG2 - 8*
     +    TG**(-3)*S*NS*MG2**2 - 8*TG**(-3)*S*MS2*MG2 + 8*TG**(-3)*S*
     +    MG2**2 + 24*TG**(-3)*NS*MS2*MG2**2 - 24*TG**(-3)*NS*MS2**2*
     +    MG2 + 8*TG**(-3)*NS*MS2**3 - 8*TG**(-3)*NS*MG2**3 - 24*
     +    TG**(-3)*MS2*MG2**2 + 24*TG**(-3)*MS2**2*MG2 - 8*TG**(-3)*
     +    MS2**3 + 8*TG**(-3)*MG2**3 + 8*TG**(-2)*S*NS*MS2 - 8*TG**(-2)
     +    *S*NS*MG2 - 8*TG**(-2)*S*MS2 + 8*TG**(-2)*S*MG2 + 8*TG**(-2)*
     +    T1*NS*MS2 - 8*TG**(-2)*T1*NS*MG2 - 8*TG**(-2)*T1*MS2 + 8*
     +    TG**(-2)*T1*MG2 + 16*TG**(-2)*NS*MS2*MG2 - 8*TG**(-2)*NS*
     +    MS2**2 - 8*TG**(-2)*NS*MG2**2 - 16*TG**(-2)*MS2*MG2 + 8*
     +    TG**(-2)*MS2**2 + 8*TG**(-2)*MG2**2 )
     +
      MQPLRV = MQPLRV + SK1B0E(3)*N*CF * ( 8*TG**(-3)*S*MS2*MG2 - 8*
     +    TG**(-3)*S*MG2*MT2 - 8*TG**(-3)*S*MG2**2 + 16*TG**(-3)*MS2*
     +    MG2*MT2 + 24*TG**(-3)*MS2*MG2**2 - 24*TG**(-3)*MS2**2*MG2 - 8
     +    *TG**(-3)*MS2**2*MT2 + 8*TG**(-3)*MS2**3 - 8*TG**(-3)*MG2**2*
     +    MT2 - 8*TG**(-3)*MG2**3 + 8*TG**(-2)*S*MS2 - 8*TG**(-2)*S*MG2
     +     - 8*TG**(-2)*S*MT2 + 8*TG**(-2)*T1*MS2 - 8*TG**(-2)*T1*MG2
     +     - 8*TG**(-2)*T1*MT2 + 16*TG**(-2)*MS2*MG2 + 8*TG**(-2)*MS2*
     +    MT2 - 8*TG**(-2)*MS2**2 - 8*TG**(-2)*MG2*MT2 - 8*TG**(-2)*
     +    MG2**2 )
     +
      MQPLRV = MQPLRV + SK1BP(1)*N*CF**2 * (  - 16*TG**(-2)*S*MS2*MG2
     +     - 16*TG**(-2)*T1*MS2*MG2 + 16*TG**(-2)*T1*MS2**2 - 16*
     +    TG**(-1)*S*MS2 - 16*TG**(-1)*T1*MS2 )
     +
      MQPLRV = MQPLRV + SK1BP(2)*N*CF**2 * ( 8*TG**(-2)*S*MS2*MG2 - 8*
     +    TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*MS2*MG2 - 8*TG**(-2)*T1*
     +    MS2**2 - 8*TG**(-2)*T1*MG2**2 + 8*TG**(-1)*S*MS2 - 8*TG**(-1)
     +    *S*MG2 + 8*TG**(-1)*T1*MS2 - 8*TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1BP(3)*N*CF**2 * (  - 4*TG**(-2)*S*MS2*MG2 + 
     +    4*TG**(-2)*S*MG2**2 - 8*TG**(-2)*T1*MS2*MG2 + 4*TG**(-2)*T1*
     +    MS2**2 + 4*TG**(-2)*T1*MG2**2 - 4*TG**(-1)*S*MS2 + 4*TG**(-1)
     +    *S*MG2 - 4*TG**(-1)*T1*MS2 + 4*TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0A(1)*N*CF**2 * ( 8*TG**(-2)*S*T1*MS2 - 8*
     +    TG**(-2)*S*T1*MG2 - 16*TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*S*
     +    MS2**2 + 8*TG**(-2)*S*MG2**2 - 8*TG**(-1)*S*T1 - 8*TG**(-1)*
     +    S**2 )
     +
      MQPLRV = MQPLRV + SK1C0A(1)*N**2*CF * (  - 4*TG**(-2)*S*T1*MS2 + 
     +    4*TG**(-2)*S*T1*MG2 + 8*TG**(-2)*S*MS2*MG2 - 4*TG**(-2)*S*
     +    MS2**2 - 4*TG**(-2)*S*MG2**2 + 4*TG**(-1)*S*T1 + 4*TG**(-1)*
     +    S**2 )
     +
      MQPLRV = MQPLRV + SK1C0A(5)*N*CF**2 * (  - 16*TG**(-2)*S*T1*MS2
     +     + 16*TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MS2**2 + 8*TG**(-2)*
     +    S**2*T1 + 8*TG**(-2)*S**2*MS2 - 8*TG**(-2)*S**2*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0A(5)*N**2*CF * ( 8*TG**(-2)*S*T1*MS2 - 8*
     +    TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*S*MS2**2 - 4*TG**(-2)*S**2*T1
     +     - 4*TG**(-2)*S**2*MS2 + 4*TG**(-2)*S**2*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0B(1)*N*CF**2 * (  - 16*TG**(-2)*S*T1*MS2
     +     + 16*TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MS2**2 + 8*TG**(-2)*
     +    S**2*T1 + 8*TG**(-2)*S**2*MS2 - 8*TG**(-2)*S**2*MG2 - 16*
     +    TG**(-1)*S*T1 + 32*TG**(-1)*S*MS2 - 16*TG**(-1)*S**2 + 32*
     +    TG**(-1)*T1*MS2 )
     +
      MQPLRV = MQPLRV + SK1C0B(1)*N**2*CF * ( 8*TG**(-2)*S*T1*MS2 - 8*
     +    TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*S*MS2**2 - 4*TG**(-2)*S**2*T1
     +     - 4*TG**(-2)*S**2*MS2 + 4*TG**(-2)*S**2*MG2 + 8*TG**(-1)*S*
     +    T1 - 16*TG**(-1)*S*MS2 + 8*TG**(-1)*S**2 - 16*TG**(-1)*T1*MS2
     +     )
     +
      MQPLRV = MQPLRV + SK1C0B(4)*N*CF**2 * (  - 8*TG**(-2)*S*T1*MS2 - 
     +    8*TG**(-2)*S*T1*MG2 - 8*TG**(-2)*S*MS2**2 + 8*TG**(-2)*S*
     +    MG2**2 + 24*TG**(-2)*S**2*T1*S1**(-1)*MS2 + 8*TG**(-2)*S**2*
     +    T1*S1**(-1)*MG2 - 8*TG**(-2)*S**2*T1 - 16*TG**(-2)*S**2*
     +    S1**(-1)*MS2*MG2 + 24*TG**(-2)*S**2*S1**(-1)*MS2**2 - 8*
     +    TG**(-2)*S**2*S1**(-1)*MG2**2 - 8*TG**(-2)*S**2*MS2 + 8*
     +    TG**(-2)*S**2*MG2 - 32*TG**(-1)*S*T1*S1**(-1)*MS2 + 16*
     +    TG**(-1)*S*T1 - 16*TG**(-1)*S*MS2 + 16*TG**(-1)*S*MG2 + 8*
     +    TG**(-1)*S**2*T1*S1**(-1) - 24*TG**(-1)*S**2*S1**(-1)*MS2 - 8
     +    *TG**(-1)*S**2*S1**(-1)*MG2 + 16*TG**(-1)*S**2 - 16*TG**(-1)*
     +    T1*MS2 + 16*TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0B(4)*N**2*CF * ( 4*TG**(-2)*S*T1*MS2 + 4*
     +    TG**(-2)*S*T1*MG2 + 4*TG**(-2)*S*MS2**2 - 4*TG**(-2)*S*MG2**2
     +     - 12*TG**(-2)*S**2*T1*S1**(-1)*MS2 - 4*TG**(-2)*S**2*T1*
     +    S1**(-1)*MG2 + 4*TG**(-2)*S**2*T1 + 8*TG**(-2)*S**2*S1**(-1)*
     +    MS2*MG2 - 12*TG**(-2)*S**2*S1**(-1)*MS2**2 + 4*TG**(-2)*S**2*
     +    S1**(-1)*MG2**2 + 4*TG**(-2)*S**2*MS2 - 4*TG**(-2)*S**2*MG2
     +     + 16*TG**(-1)*S*T1*S1**(-1)*MS2 - 8*TG**(-1)*S*T1 + 8*
     +    TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2*T1*
     +    S1**(-1) + 12*TG**(-1)*S**2*S1**(-1)*MS2 + 4*TG**(-1)*S**2*
     +    S1**(-1)*MG2 - 8*TG**(-1)*S**2 + 8*TG**(-1)*T1*MS2 - 8*
     +    TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,1)*N*CF**2 * ( 16*TG**(-2)*S*T1*MS2 - 
     +    32*TG**(-2)*S*MS2*MG2 + 16*TG**(-2)*S*MS2**2 + 16*TG**(-2)*S*
     +    MG2**2 - 16*TG**(-2)*S**2*T1 - 16*TG**(-2)*S**2*MS2 + 16*
     +    TG**(-2)*S**2*MG2 - 32*TG**(-2)*T1*MS2*MG2 + 16*TG**(-2)*T1*
     +    MS2**2 + 16*TG**(-2)*T1*MG2**2 - 32*TG**(-1)*S*MS2 + 32*
     +    TG**(-1)*S*MG2 + 16*TG**(-1)*S**2 - 32*TG**(-1)*T1*MS2 + 32*
     +    TG**(-1)*T1*MG2 + 16*S + 16*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,1)*N**2*CF * (  - 8*TG**(-2)*S*T1*MS2
     +     + 16*TG**(-2)*S*MS2*MG2 - 8*TG**(-2)*S*MS2**2 - 8*TG**(-2)*S
     +    *MG2**2 + 8*TG**(-2)*S**2*T1 + 8*TG**(-2)*S**2*MS2 - 8*
     +    TG**(-2)*S**2*MG2 + 16*TG**(-2)*T1*MS2*MG2 - 8*TG**(-2)*T1*
     +    MS2**2 - 8*TG**(-2)*T1*MG2**2 + 16*TG**(-1)*S*MS2 - 16*
     +    TG**(-1)*S*MG2 - 8*TG**(-1)*S**2 + 16*TG**(-1)*T1*MS2 - 16*
     +    TG**(-1)*T1*MG2 - 8*S - 8*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,2)*N*CF**2 * (  - 16*TG**(-2)*S*T1*MS2
     +     + 16*TG**(-2)*S*T1*MG2 + 32*TG**(-2)*S*MS2*MG2 - 16*TG**(-2)
     +    *S*MS2**2 - 16*TG**(-2)*S*MG2**2 + 16*TG**(-2)*S**2*T1 + 16*
     +    TG**(-2)*S**2*MS2 - 16*TG**(-2)*S**2*MG2 - 16*TG**(-1)*S*T1
     +     + 48*TG**(-1)*S*MS2 - 32*TG**(-1)*S*MG2 - 32*TG**(-1)*S**2
     +     + 32*TG**(-1)*T1*MS2 - 16*TG**(-1)*T1*MG2 - 16*S - 16*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(3,2)*N**2*CF * ( 4*TG**(-2)*S*T1*MS2 - 4
     +    *TG**(-2)*S*T1*MG2 - 8*TG**(-2)*S*MS2*MG2 + 4*TG**(-2)*S*
     +    MS2**2 + 4*TG**(-2)*S*MG2**2 - 4*TG**(-2)*S**2*T1 - 4*
     +    TG**(-2)*S**2*MS2 + 4*TG**(-2)*S**2*MG2 + 4*TG**(-1)*S*T1 - 
     +    12*TG**(-1)*S*MS2 + 8*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2 - 8*
     +    TG**(-1)*T1*MS2 + 4*TG**(-1)*T1*MG2 + 4*S + 4*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(4,1)*N*CF**2 * ( 16*TG**(-2)*S*T1*MG2 + 
     +    16*TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MG2**2 - 16*TG**(-1)*S*
     +    MG2 - 16*TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0C(4,1)*N**2*CF * (  - 8*TG**(-2)*S*T1*MS2
     +     - 4*TG**(-2)*S*T1*MG2 + 12*TG**(-2)*S*MS2*MG2 - 8*TG**(-2)*S
     +    *MS2**2 - 4*TG**(-2)*S*MG2**2 + 4*TG**(-2)*S**2*T1 + 4*
     +    TG**(-2)*S**2*MS2 - 4*TG**(-2)*S**2*MG2 + 16*TG**(-2)*T1*MS2*
     +    MG2 - 8*TG**(-2)*T1*MS2**2 - 8*TG**(-2)*T1*MG2**2 + 12*
     +    TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2 + 8*
     +    TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 - 4*S - 4*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(5,1)*N*CF**2 * ( 16*TG**(-2)*S*T1*MG2 + 
     +    16*TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MG2**2 + 16*TG**(-1)*S*
     +    T1 + 16*TG**(-1)*S*MS2 - 32*TG**(-1)*S*MG2 - 16*TG**(-1)*T1*
     +    MG2 - 16*S - 16*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(5,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2**2*MG2 + 8*TG**(-2)*S*T1**(-1)*MS2**3 - 8*
     +    TG**(-2)*S*T1*MS2 - 4*TG**(-2)*S*T1*MG2 - 4*TG**(-2)*S*MS2*
     +    MG2 - 4*TG**(-2)*S*MG2**2 + 4*TG**(-2)*S**2*T1 + 4*TG**(-2)*
     +    S**2*MS2 - 4*TG**(-2)*S**2*MG2 + 8*TG**(-2)*T1*MS2**2 - 8*
     +    TG**(-2)*T1*MG2**2 - 8*TG**(-1)*S*T1 + 4*TG**(-1)*S*MS2 - 4*
     +    TG**(-1)*S**2 + 8*TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 + 4*S
     +     + 4*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(6,1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 16*TG**(-2)*S*T1**(-1)*MS2**2*MG2 - 16*TG**(-2)*
     +    S*MS2*MG2 + 32*TG**(-2)*S*MG2**2 - 32*TG**(-2)*T1*MS2*MG2 + 
     +    32*TG**(-2)*T1*MG2**2 + 16*TG**(-1)*S*MG2 + 16*TG**(-1)*T1*
     +    MG2 )
     +
      MQPLRV = MQPLRV + SK1C0C(6,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2**2 + 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 8*
     +    TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*
     +    MS2*MG2 - 16*TG**(-2)*T1*MG2**2 - 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*T1*MG2 )
     +
      MQPLRV = MQPLRV + SK1C0C(6,2)*N*CF**2 * ( 16*TG**(-2)*S*T1*MS2 - 
     +    16*TG**(-2)*S*T1*MG2 - 32*TG**(-2)*S*MS2*MG2 + 16*TG**(-2)*S*
     +    MS2**2 + 16*TG**(-2)*S*MG2**2 - 16*TG**(-2)*S**2*T1 - 16*
     +    TG**(-2)*S**2*MS2 + 16*TG**(-2)*S**2*MG2 + 16*TG**(-1)*S*T1
     +     - 48*TG**(-1)*S*MS2 + 32*TG**(-1)*S*MG2 + 32*TG**(-1)*S**2
     +     - 32*TG**(-1)*T1*MS2 + 16*TG**(-1)*T1*MG2 + 16*S + 16*T1 )
     +
      MQPLRV = MQPLRV + SK1C0C(6,2)*N**2*CF * (  - 4*TG**(-2)*S*T1*MS2
     +     + 4*TG**(-2)*S*T1*MG2 + 8*TG**(-2)*S*MS2*MG2 - 4*TG**(-2)*S*
     +    MS2**2 - 4*TG**(-2)*S*MG2**2 + 4*TG**(-2)*S**2*T1 + 4*
     +    TG**(-2)*S**2*MS2 - 4*TG**(-2)*S**2*MG2 - 4*TG**(-1)*S*T1 + 
     +    12*TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 - 8*TG**(-1)*S**2 + 8*
     +    TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 - 4*S - 4*T1 )
     +
      MQPLRV = MQPLRV + SK1D0(4,1)*N*CF**2 * ( 16*TG**(-2)*S*T1*MS2*MG2
     +     - 8*TG**(-2)*S*T1*MS2**2 - 8*TG**(-2)*S*T1*MG2**2 - 24*
     +    TG**(-2)*S*MS2*MG2**2 + 24*TG**(-2)*S*MS2**2*MG2 - 8*TG**(-2)
     +    *S*MS2**3 + 8*TG**(-2)*S*MG2**3 - 8*TG**(-2)*S**2*T1*MS2 - 8*
     +    TG**(-2)*S**2*T1*MG2 - 8*TG**(-2)*S**2*MS2**2 + 8*TG**(-2)*
     +    S**2*MG2**2 + 8*TG**(-1)*S*T1*MS2 - 8*TG**(-1)*S*T1*MG2 + 8*
     +    TG**(-1)*S**2*MS2 - 8*TG**(-1)*S**2*MG2 - 8*S*T1 - 8*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(4,1)*N**2*CF * (  - 8*TG**(-2)*S*T1*MS2*
     +    MG2 + 4*TG**(-2)*S*T1*MS2**2 + 4*TG**(-2)*S*T1*MG2**2 + 12*
     +    TG**(-2)*S*MS2*MG2**2 - 12*TG**(-2)*S*MS2**2*MG2 + 4*TG**(-2)
     +    *S*MS2**3 - 4*TG**(-2)*S*MG2**3 + 4*TG**(-2)*S**2*T1*MS2 + 4*
     +    TG**(-2)*S**2*T1*MG2 + 4*TG**(-2)*S**2*MS2**2 - 4*TG**(-2)*
     +    S**2*MG2**2 - 4*TG**(-1)*S*T1*MS2 + 4*TG**(-1)*S*T1*MG2 - 4*
     +    TG**(-1)*S**2*MS2 + 4*TG**(-1)*S**2*MG2 + 4*S*T1 + 4*S**2 )
     +
      MQPLRV = MQPLRV + SK1D0(7,1)*N*CF**2 * (  - 16*TG**(-2)*S*T1*MS2*
     +    MG2 + 48*TG**(-2)*S*T1*MS2**2 + 16*TG**(-2)*S*MS2*MG2**2 - 64
     +    *TG**(-2)*S*MS2**2*MG2 + 48*TG**(-2)*S*MS2**3 - 32*TG**(-2)*
     +    S**2*T1*MS2 + 32*TG**(-2)*S**2*MS2*MG2 - 32*TG**(-2)*S**2*
     +    MS2**2 + 32*TG**(-1)*S*T1*MS2 - 16*TG**(-1)*S*T1*MG2 + 32*
     +    TG**(-1)*S*MS2*MG2 - 32*TG**(-1)*S*MS2**2 - 8*TG**(-1)*S**2*
     +    T1 + 24*TG**(-1)*S**2*MS2 - 8*TG**(-1)*S**2*MG2 + 32*TG**(-1)
     +    *T1*MS2*MG2 - 32*TG**(-1)*T1*MS2**2 )
     +
      MQPLRV = MQPLRV + SK1D0(7,1)*N**2*CF * ( 8*TG**(-2)*S*T1*MS2*MG2
     +     - 24*TG**(-2)*S*T1*MS2**2 - 8*TG**(-2)*S*MS2*MG2**2 + 32*
     +    TG**(-2)*S*MS2**2*MG2 - 24*TG**(-2)*S*MS2**3 + 16*TG**(-2)*
     +    S**2*T1*MS2 - 16*TG**(-2)*S**2*MS2*MG2 + 16*TG**(-2)*S**2*
     +    MS2**2 - 16*TG**(-1)*S*T1*MS2 + 8*TG**(-1)*S*T1*MG2 - 16*
     +    TG**(-1)*S*MS2*MG2 + 16*TG**(-1)*S*MS2**2 + 4*TG**(-1)*S**2*
     +    T1 - 12*TG**(-1)*S**2*MS2 + 4*TG**(-1)*S**2*MG2 - 16*TG**(-1)
     +    *T1*MS2*MG2 + 16*TG**(-1)*T1*MS2**2 )
     +
      MQPLRV = MQPLRV + SK1D0(8,1)*N*CF**2 * (  - 48*TG**(-2)*S*T1*MS2*
     +    MG2 + 16*TG**(-2)*S*T1*MS2**2 + 32*TG**(-2)*S*T1*MG2**2 + 80*
     +    TG**(-2)*S*MS2*MG2**2 - 64*TG**(-2)*S*MS2**2*MG2 + 16*
     +    TG**(-2)*S*MS2**3 - 32*TG**(-2)*S*MG2**3 - 16*TG**(-2)*S**2*
     +    T1*MS2 + 32*TG**(-2)*S**2*T1*MG2 + 48*TG**(-2)*S**2*MS2*MG2
     +     - 16*TG**(-2)*S**2*MS2**2 - 32*TG**(-2)*S**2*MG2**2 - 16*
     +    TG**(-1)*S*T1*MS2 - 16*TG**(-1)*S*T1*MG2 + 144*TG**(-1)*S*MS2
     +    *MG2 - 64*TG**(-1)*S*MS2**2 - 80*TG**(-1)*S*MG2**2 + 16*
     +    TG**(-1)*S**2*T1 + 32*TG**(-1)*S**2*MS2 - 80*TG**(-1)*S**2*
     +    MG2 + 64*TG**(-1)*T1*MS2*MG2 - 32*TG**(-1)*T1*MS2**2 - 32*
     +    TG**(-1)*T1*MG2**2 - 32*S*T1 + 32*S*MS2 - 48*S*MG2 - 32*S**2
     +     + 16*T1*MS2 - 32*T1*MG2 - 16*T1**2 )
     +
      MQPLRV = MQPLRV + SK1D0(8,1)*N**2*CF * ( 12*TG**(-2)*S*T1*MS2*MG2
     +     - 4*TG**(-2)*S*T1*MS2**2 - 8*TG**(-2)*S*T1*MG2**2 - 20*
     +    TG**(-2)*S*MS2*MG2**2 + 16*TG**(-2)*S*MS2**2*MG2 - 4*TG**(-2)
     +    *S*MS2**3 + 8*TG**(-2)*S*MG2**3 + 4*TG**(-2)*S**2*T1*MS2 - 8*
     +    TG**(-2)*S**2*T1*MG2 - 12*TG**(-2)*S**2*MS2*MG2 + 4*TG**(-2)*
     +    S**2*MS2**2 + 8*TG**(-2)*S**2*MG2**2 + 4*TG**(-1)*S*T1*MS2 + 
     +    4*TG**(-1)*S*T1*MG2 - 36*TG**(-1)*S*MS2*MG2 + 16*TG**(-1)*S*
     +    MS2**2 + 20*TG**(-1)*S*MG2**2 - 4*TG**(-1)*S**2*T1 - 8*
     +    TG**(-1)*S**2*MS2 + 20*TG**(-1)*S**2*MG2 - 16*TG**(-1)*T1*MS2
     +    *MG2 + 8*TG**(-1)*T1*MS2**2 + 8*TG**(-1)*T1*MG2**2 + 8*S*T1
     +     - 8*S*MS2 + 12*S*MG2 + 8*S**2 - 4*T1*MS2 + 8*T1*MG2 + 4*
     +    T1**2 )
     +
      MQPLRV = MQPLRV + SOF1(1)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 - 
     +    128*TG**(-2)*S*MS2**2 + 64*TG**(-2)*S**2*MS2 + 128*TG**(-2)*
     +    T1*U1*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(1)*N**2*CF * ( 32*TG**(-2)*S*T1*U1 + 64*
     +    TG**(-2)*S*MS2**2 - 32*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1*U1*
     +    MS2 )
     +
      MQPLRV = MQPLRV + SOF1(2)*N*CF**2 * ( 32*TG**(-2)*S*MS2**2 - 32*
     +    TG**(-2)*T1*U1*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(3)*N*CF**2 * ( 32*TG**(-2)*S*MS2**2 - 32*
     +    TG**(-2)*T1*U1*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(4)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S*T1*MS2 + 64*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1**2
     +    *U1 )
     +
      MQPLRV = MQPLRV + SOF1(4)*N**2*CF * ( 16*TG**(-2)*S*T1*U1 - 16*
     +    TG**(-2)*S*T1*MS2 - 16*TG**(-2)*S**2*MS2 + 16*TG**(-2)*T1**2*
     +    U1 )
     +
      MQPLRV = MQPLRV + SOF1(5)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MS2 + 
     +    32*TG**(-2)*T1**2*U1 )
     +
      MQPLRV = MQPLRV + SOF1(5)*N**2*CF * ( 16*TG**(-2)*S*T1*MS2 - 16*
     +    TG**(-2)*T1**2*U1 )
     +
      MQPLRV = MQPLRV + SOF1(6)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MS2 + 
     +    32*TG**(-2)*T1**2*U1 )
     +
      MQPLRV = MQPLRV + SOF1(6)*N**2*CF * ( 16*TG**(-2)*S*T1*MS2 - 16*
     +    TG**(-2)*T1**2*U1 )
     +
      MQPLRV = MQPLRV + SOF1(7)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S*T1*MS2 + 64*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1**2
     +    *U1 )
     +
      MQPLRV = MQPLRV + SOF1(7)*N**2*CF * ( 16*TG**(-2)*S*T1*U1 - 16*
     +    TG**(-2)*S*T1*MS2 - 16*TG**(-2)*S**2*MS2 + 16*TG**(-2)*T1**2*
     +    U1 )
     +
      MQPLRV = MQPLRV + SOF1(8)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S**2*MS2 )
     +
      MQPLRV = MQPLRV + SOF1(8)*N**2*CF * ( 32*TG**(-2)*S*T1*U1 - 32*
     +    TG**(-2)*S**2*MS2 )

         END IF

         IF (IFL.EQ.1) THEN
      MQQLLV = 0.D0
      MQQLLV = MQQLLV + N**(-1)*CF * ( 32*TG**(-1)*UG**(-2)*S*T**(-1)*
     +    T1*NS*MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*MG2**2
     +     - 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MG2**3 - 32*TG**(-1)*
     +    UG**(-2)*S*T1*NS*MG2 - 32*TG**(-1)*UG**(-2)*S*NS*MS2*MG2 + 32
     +    *TG**(-1)*UG**(-2)*S*NS*MG2**2 - 32*UG**(-2)*S*T**(-1)*NS*
     +    MG2**2 + 32*UG**(-2)*S*NS*MG2 )
     +
      MQQLLV = MQQLLV + N*CF * ( 8*TG**(-3)*S*T**(-1)*T1*NS*MG2**2 - 8*
     +    TG**(-3)*S*T**(-1)*NS*MG2**3 + 8*TG**(-3)*S*T**(-1)*MG2**2*
     +    MT2 - 8*TG**(-3)*S*T1*NS*MG2 + 8*TG**(-3)*S*NS*MG2**2 - 8*
     +    TG**(-3)*S*MG2*MT2 - 8*TG**(-2)*UG**(-1)*S*T**(-1)*MG2**3 + 8
     +    *TG**(-2)*UG**(-1)*S*MG2**2 - 8*TG**(-2)*S*T**(-1)*NS*MS2*MG2
     +     - 8*TG**(-2)*S*T**(-1)*NS*MG2**2 + 8*TG**(-2)*S*T**(-1)*MG2*
     +    MT2 - 4*TG**(-2)*S*T1*NS + 12*TG**(-2)*S*NS*MG2 - 4*TG**(-2)*
     +    S*MT2 - 8*TG**(-1)*UG**(-2)*S*T**(-1)*T1*NS*MG2**2 - 24*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*MG2**2 + 24*TG**(-1)*
     +    UG**(-2)*S*T**(-1)*NS*MG2**3 - 8*TG**(-1)*UG**(-2)*S*U**(-1)*
     +    MG2**3 + 8*TG**(-1)*UG**(-2)*S*T1*NS*MG2 + 24*TG**(-1)*
     +    UG**(-2)*S*NS*MS2*MG2 - 24*TG**(-1)*UG**(-2)*S*NS*MG2**2 + 8*
     +    TG**(-1)*UG**(-2)*S*MG2**2 + 8*TG**(-1)*UG**(-2)*S**2*T**(-1)
     +    *NS*MG2**2 - 8*TG**(-1)*UG**(-2)*S**2*NS*MG2 + 8*TG**(-1)*
     +    UG**(-1)*S*T**(-1)*NS*MG2**2 - 8*TG**(-1)*UG**(-1)*S*T**(-1)*
     +    MG2**2 )
     +
      MQQLLV = MQQLLV + N*CF * (  - 8*TG**(-1)*UG**(-1)*S*U**(-1)*
     +    MG2**2 - 8*TG**(-1)*UG**(-1)*S*NS*MG2 + 16*TG**(-1)*UG**(-1)*
     +    S*MG2 + 4*TG**(-1)*S*NS - 8*UG**(-3)*S*U**(-1)*T1*NS*MG2**2
     +     - 8*UG**(-3)*S*U**(-1)*NS*MG2**3 + 8*UG**(-3)*S*U**(-1)*
     +    MG2**2*MT2 + 8*UG**(-3)*S*T1*NS*MG2 + 8*UG**(-3)*S*NS*MG2**2
     +     - 8*UG**(-3)*S*MG2*MT2 - 8*UG**(-3)*S**2*U**(-1)*NS*MG2**2
     +     + 8*UG**(-3)*S**2*NS*MG2 + 8*UG**(-2)*S*T**(-1)*T1*NS*MG2 - 
     +    8*UG**(-2)*S*T**(-1)*NS*MS2*MG2 + 24*UG**(-2)*S*T**(-1)*NS*
     +    MG2**2 - 8*UG**(-2)*S*U**(-1)*T1*NS*MG2 - 16*UG**(-2)*S*
     +    U**(-1)*NS*MG2**2 + 8*UG**(-2)*S*U**(-1)*MG2*MT2 + 4*UG**(-2)
     +    *S*T1*NS - 4*UG**(-2)*S*NS*MG2 - 4*UG**(-2)*S*MT2 + 8*
     +    UG**(-2)*S**2*T**(-1)*NS*MG2 - 8*UG**(-2)*S**2*U**(-1)*NS*MG2
     +     + 4*UG**(-2)*S**2*NS + 8*UG**(-1)*S*T**(-1)*NS*MG2 - 8*
     +    UG**(-1)*S*U**(-1)*NS*MG2 + 4*UG**(-1)*S*NS )
     +
      MQQLLV = MQQLLV + N*CF**2 * ( 48*TG**(-2)*UG**(-1)*S*T**(-1)*T1*
     +    NS*MG2**2 - 16*TG**(-2)*UG**(-1)*S*T**(-1)*NS*MG2**3 + 16*
     +    TG**(-2)*UG**(-1)*S*T**(-1)*MG2**2*MT2 - 48*TG**(-2)*UG**(-1)
     +    *S*T1*NS*MG2 + 16*TG**(-2)*UG**(-1)*S*NS*MG2**2 - 16*TG**(-2)
     +    *UG**(-1)*S*MG2*MT2 + 16*TG**(-2)*UG**(-1)*S**2*T**(-1)*NS*
     +    MG2**2 - 16*TG**(-2)*UG**(-1)*S**2*NS*MG2 + 16*TG**(-2)*S*
     +    T**(-1)*NS*MG2**2 - 16*TG**(-2)*S*NS*MG2 + 16*TG**(-1)*
     +    UG**(-2)*S*U**(-1)*NS*MS2*MG2**2 - 32*TG**(-1)*UG**(-2)*S*
     +    U**(-1)*NS*MG2**3 + 16*TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2*MT2
     +     - 16*TG**(-1)*UG**(-2)*S*NS*MS2*MG2 + 32*TG**(-1)*UG**(-2)*S
     +    *NS*MG2**2 - 16*TG**(-1)*UG**(-2)*S*MG2*MT2 - 16*TG**(-1)*
     +    UG**(-2)*S**2*U**(-1)*NS*MG2**2 + 16*TG**(-1)*UG**(-2)*S**2*
     +    NS*MG2 - 48*TG**(-1)*UG**(-1)*S*T**(-1)*NS*MS2*MG2 + 16*
     +    TG**(-1)*UG**(-1)*S*T**(-1)*MG2*MT2 + 16*TG**(-1)*UG**(-1)*S*
     +    U**(-1)*NS*MS2*MG2 - 48*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MG2**2
     +     + 16*TG**(-1)*UG**(-1)*S*U**(-1)*MG2*MT2 )
     +
      MQQLLV = MQQLLV + N*CF**2 * (  - 16*TG**(-1)*UG**(-1)*S*NS*MS2 + 
     +    80*TG**(-1)*UG**(-1)*S*NS*MG2 - 16*TG**(-1)*UG**(-1)*S*MT2 + 
     +    16*TG**(-1)*UG**(-1)*S**2*T**(-1)*NS*MG2 - 16*TG**(-1)*
     +    UG**(-1)*S**2*U**(-1)*NS*MG2 + 16*TG**(-1)*UG**(-1)*S**2*NS
     +     + 16*TG**(-1)*S*T**(-1)*NS*MG2 - 16*TG**(-1)*S*U**(-1)*NS*
     +    MG2 + 16*TG**(-1)*S*NS - 16*UG**(-2)*S*U**(-1)*NS*MG2**2 + 16
     +    *UG**(-2)*S*NS*MG2 + 16*UG**(-1)*S*T**(-1)*NS*MG2 - 16*
     +    UG**(-1)*S*U**(-1)*NS*MG2 + 16*UG**(-1)*S*NS )
     +
      MQQLLV = MQQLLV + N**2*CF * ( 8*TG**(-3)*S*T**(-1)*MG2**3 - 8*
     +    TG**(-3)*S*MG2**2 - 24*TG**(-2)*UG**(-1)*S*T**(-1)*T1*NS*
     +    MG2**2 + 8*TG**(-2)*UG**(-1)*S*T**(-1)*NS*MG2**3 - 8*TG**(-2)
     +    *UG**(-1)*S*T**(-1)*MG2**2*MT2 + 24*TG**(-2)*UG**(-1)*S*T1*NS
     +    *MG2 - 8*TG**(-2)*UG**(-1)*S*NS*MG2**2 + 8*TG**(-2)*UG**(-1)*
     +    S*MG2*MT2 - 8*TG**(-2)*UG**(-1)*S**2*T**(-1)*NS*MG2**2 + 8*
     +    TG**(-2)*UG**(-1)*S**2*NS*MG2 - 8*TG**(-2)*S*T**(-1)*NS*
     +    MG2**2 + 8*TG**(-2)*S*T**(-1)*MG2**2 + 8*TG**(-2)*S*NS*MG2 - 
     +    8*TG**(-2)*S*MG2 - 8*TG**(-1)*UG**(-2)*S*U**(-1)*NS*MS2*
     +    MG2**2 + 16*TG**(-1)*UG**(-2)*S*U**(-1)*NS*MG2**3 - 8*
     +    TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2*MT2 + 8*TG**(-1)*UG**(-2)*
     +    S*NS*MS2*MG2 - 16*TG**(-1)*UG**(-2)*S*NS*MG2**2 + 8*TG**(-1)*
     +    UG**(-2)*S*MG2*MT2 + 8*TG**(-1)*UG**(-2)*S**2*U**(-1)*NS*
     +    MG2**2 - 8*TG**(-1)*UG**(-2)*S**2*NS*MG2 + 24*TG**(-1)*
     +    UG**(-1)*S*T**(-1)*NS*MS2*MG2 - 8*TG**(-1)*UG**(-1)*S*T**(-1)
     +    *MG2*MT2 )
     +
      MQQLLV = MQQLLV + N**2*CF * (  - 8*TG**(-1)*UG**(-1)*S*U**(-1)*NS
     +    *MS2*MG2 + 24*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MG2**2 - 8*
     +    TG**(-1)*UG**(-1)*S*U**(-1)*MG2*MT2 + 8*TG**(-1)*UG**(-1)*S*
     +    NS*MS2 - 40*TG**(-1)*UG**(-1)*S*NS*MG2 + 8*TG**(-1)*UG**(-1)*
     +    S*MT2 - 8*TG**(-1)*UG**(-1)*S**2*T**(-1)*NS*MG2 + 8*TG**(-1)*
     +    UG**(-1)*S**2*U**(-1)*NS*MG2 - 8*TG**(-1)*UG**(-1)*S**2*NS - 
     +    8*TG**(-1)*S*T**(-1)*NS*MG2 + 8*TG**(-1)*S*U**(-1)*NS*MG2 - 8
     +    *TG**(-1)*S*NS + 8*UG**(-3)*S*U**(-1)*MG2**3 - 8*UG**(-3)*S*
     +    MG2**2 + 8*UG**(-2)*S*U**(-1)*NS*MG2**2 + 8*UG**(-2)*S*
     +    U**(-1)*MG2**2 - 8*UG**(-2)*S*NS*MG2 - 8*UG**(-2)*S*MG2 - 8*
     +    UG**(-1)*S*T**(-1)*NS*MG2 + 8*UG**(-1)*S*U**(-1)*NS*MG2 - 8*
     +    UG**(-1)*S*NS )
     +
      MQQLLV = MQQLLV + CF**2 * ( 32*TG**(-1)*UG**(-2)*S*T**(-1)*T1*NS*
     +    MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*MG2**2 - 32*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*NS*MG2**3 - 32*TG**(-1)*UG**(-2)*
     +    S*T1*NS*MG2 - 32*TG**(-1)*UG**(-2)*S*NS*MS2*MG2 + 32*TG**(-1)
     +    *UG**(-2)*S*NS*MG2**2 - 32*UG**(-2)*S*T**(-1)*NS*MG2**2 + 32*
     +    UG**(-2)*S*NS*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0A(1)*N*CF * (  - 8*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*MG2**3 + 8*TG**(-2)*UG**(-1)*S*MG2**2 - 8*TG**(-1)*
     +    UG**(-2)*S*U**(-1)*MG2**3 + 8*TG**(-1)*UG**(-2)*S*MG2**2 - 8*
     +    TG**(-1)*UG**(-1)*S*T**(-1)*MG2**2 - 8*TG**(-1)*UG**(-1)*S*
     +    U**(-1)*MG2**2 + 8*TG**(-1)*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0A(1)*N**2*CF * ( 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 - 8*TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*MG2**2 - 
     +    4*TG**(-2)*S*MG2 + 8*UG**(-3)*S*U**(-1)*MG2**3 - 8*UG**(-3)*S
     +    *MG2**2 + 8*UG**(-2)*S*U**(-1)*MG2**2 - 4*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N*CF * ( 8*TG**(-3)*S*T**(-1)*MG2**2*
     +    MT2 - 8*TG**(-3)*S*MG2*MT2 + 8*TG**(-2)*S*T**(-1)*MG2*MT2 - 4
     +    *TG**(-2)*S*MT2 + 8*UG**(-3)*S*U**(-1)*MG2**2*MT2 - 8*
     +    UG**(-3)*S*MG2*MT2 + 8*UG**(-2)*S*U**(-1)*MG2*MT2 - 4*
     +    UG**(-2)*S*MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N*CF**2 * ( 16*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*MG2**2*MT2 - 16*TG**(-2)*UG**(-1)*S*MG2*MT2 + 16*
     +    TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2*MT2 - 16*TG**(-1)*UG**(-2)
     +    *S*MG2*MT2 + 16*TG**(-1)*UG**(-1)*S*T**(-1)*MG2*MT2 + 16*
     +    TG**(-1)*UG**(-1)*S*U**(-1)*MG2*MT2 - 16*TG**(-1)*UG**(-1)*S*
     +    MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(2)*N**2*CF * (  - 8*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*MG2**2*MT2 + 8*TG**(-2)*UG**(-1)*S*MG2*MT2 - 8*
     +    TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2*MT2 + 8*TG**(-1)*UG**(-2)*
     +    S*MG2*MT2 - 8*TG**(-1)*UG**(-1)*S*T**(-1)*MG2*MT2 - 8*
     +    TG**(-1)*UG**(-1)*S*U**(-1)*MG2*MT2 + 8*TG**(-1)*UG**(-1)*S*
     +    MT2 )
     +
      MQQLLV = MQQLLV + SK1B0A(3)*N*CF**2 * (  - 4*TG**(-2)*S*MG2 - 4*
     +    UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0A(3)*CF**2 * ( 8*TG**(-1)*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 16*TG**(-2)*S*MG2 + 16*TG**(-1)*UG**(-2)*S*U1**(-1)
     +    *MS2*MG2**2 - 16*TG**(-1)*UG**(-2)*S*U1**(-1)*MG2**3 - 16*
     +    TG**(-1)*UG**(-2)*S*MG2**2 - 16*TG**(-1)*UG**(-2)*S**2*
     +    U1**(-1)*MG2**2 + 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MS2*MG2 - 
     +    16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)
     +    *S*MG2 - 16*TG**(-1)*UG**(-1)*S**2*U1**(-1)*MG2 + 16*UG**(-2)
     +    *S*T1**(-1)*MS2*MG2 - 16*UG**(-2)*S*T1**(-1)*MG2**2 - 16*
     +    UG**(-2)*S*MG2 - 16*UG**(-2)*S**2*T1**(-1)*MG2 - 16*UG**(-1)*
     +    S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0C(1)*CF**2 * (  - 16*TG**(-1)*UG**(-1)*S*
     +    T1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 + 32
     +    *TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*S**2*T1**(-1)
     +    *MG2 + 16*TG**(-1)*S*T1**(-1)*MG2 - 16*TG**(-1)*S*U1**(-1)*
     +    MG2 - 32*UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0C(2)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 8*TG**(-2)*S*MG2 + 16*TG**(-1)*UG**(-2)*S*U1**(-1)*
     +    MS2*MG2**2 - 16*TG**(-1)*UG**(-2)*S*U1**(-1)*MG2**3 - 16*
     +    TG**(-1)*UG**(-2)*S*MG2**2 - 16*TG**(-1)*UG**(-2)*S**2*
     +    U1**(-1)*MG2**2 + 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MS2*MG2 - 
     +    16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)
     +    *S*MG2 - 16*TG**(-1)*UG**(-1)*S**2*U1**(-1)*MG2 + 16*UG**(-2)
     +    *S*T1**(-1)*MS2*MG2 - 16*UG**(-2)*S*T1**(-1)*MG2**2 - 24*
     +    UG**(-2)*S*MG2 - 16*UG**(-2)*S**2*T1**(-1)*MG2 - 16*UG**(-1)*
     +    S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0C(2)*CF**2 * (  - 16*TG**(-1)*UG**(-1)*S*
     +    T1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 + 48
     +    *TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*S**2*T1**(-1)
     +    *MG2 + 16*TG**(-1)*S*T1**(-1)*MG2 - 16*TG**(-1)*S*U1**(-1)*
     +    MG2 - 32*UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N**(-1)*CF * (  - 32*TG**(-1)*
     +    UG**(-2)*S*T**(-1)*T1*NS*MG2**2 + 32*TG**(-1)*UG**(-2)*S*
     +    T**(-1)*T1*MG2**2 - 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*
     +    MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MG2**3 + 32*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*MS2*MG2**2 - 32*TG**(-1)*UG**(-2)
     +    *S*T**(-1)*MG2**3 + 32*UG**(-2)*S*T**(-1)*NS*MG2**2 - 32*
     +    UG**(-2)*S*T**(-1)*MG2**2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*T1
     +    *NS*MG2**2 + 8*TG**(-3)*S*T**(-1)*T1*MG2**2 + 8*TG**(-3)*S*
     +    T**(-1)*NS*MG2**3 - 8*TG**(-3)*S*T**(-1)*MG2**3 - 8*TG**(-3)*
     +    S*NS*MG2**2 + 8*TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*NS*
     +    MS2*MG2 + 8*TG**(-2)*S*T**(-1)*NS*MG2**2 - 8*TG**(-2)*S*
     +    T**(-1)*MS2*MG2 - 8*TG**(-2)*S*T**(-1)*MG2**2 - 8*TG**(-2)*S*
     +    NS*MG2 + 8*TG**(-2)*S*MG2 + 8*TG**(-1)*UG**(-2)*S*T**(-1)*T1*
     +    NS*MG2**2 - 8*TG**(-1)*UG**(-2)*S*T**(-1)*T1*MG2**2 + 24*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*MG2**2 - 24*TG**(-1)*
     +    UG**(-2)*S*T**(-1)*NS*MG2**3 - 24*TG**(-1)*UG**(-2)*S*T**(-1)
     +    *MS2*MG2**2 + 24*TG**(-1)*UG**(-2)*S*T**(-1)*MG2**3 - 8*
     +    TG**(-1)*UG**(-2)*S**2*T**(-1)*NS*MG2**2 + 8*TG**(-1)*
     +    UG**(-2)*S**2*T**(-1)*MG2**2 - 8*TG**(-1)*UG**(-1)*S*T**(-1)*
     +    NS*MG2**2 + 8*TG**(-1)*UG**(-1)*S*T**(-1)*MG2**2 - 16*
     +    TG**(-1)*UG**(-1)*S*T1**(-1)*MG2**2 + 32*TG**(-1)*UG**(-1)*S*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF * ( 16*TG**(-1)*UG**(-1)*S**2*
     +    T1**(-1)*MG2 + 16*TG**(-1)*S*T1**(-1)*MG2 - 8*UG**(-2)*S*
     +    T**(-1)*T1*NS*MG2 + 8*UG**(-2)*S*T**(-1)*T1*MG2 + 8*UG**(-2)*
     +    S*T**(-1)*NS*MS2*MG2 - 24*UG**(-2)*S*T**(-1)*NS*MG2**2 - 8*
     +    UG**(-2)*S*T**(-1)*MS2*MG2 + 24*UG**(-2)*S*T**(-1)*MG2**2 - 8
     +    *UG**(-2)*S**2*T**(-1)*NS*MG2 + 8*UG**(-2)*S**2*T**(-1)*MG2
     +     - 8*UG**(-1)*S*T**(-1)*NS*MG2 + 8*UG**(-1)*S*T**(-1)*MG2 - 
     +    32*UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF**2 * (  - 48*TG**(-2)*UG**(-1)
     +    *S*T**(-1)*T1*NS*MG2**2 + 48*TG**(-2)*UG**(-1)*S*T**(-1)*T1*
     +    MG2**2 + 16*TG**(-2)*UG**(-1)*S*T**(-1)*NS*MG2**3 - 16*
     +    TG**(-2)*UG**(-1)*S*T**(-1)*MG2**3 - 16*TG**(-2)*UG**(-1)*S*
     +    NS*MG2**2 + 16*TG**(-2)*UG**(-1)*S*MG2**2 - 16*TG**(-2)*
     +    UG**(-1)*S**2*T**(-1)*NS*MG2**2 + 16*TG**(-2)*UG**(-1)*S**2*
     +    T**(-1)*MG2**2 - 16*TG**(-2)*S*T**(-1)*NS*MG2**2 + 16*
     +    TG**(-2)*S*T**(-1)*MG2**2 - 32*TG**(-2)*S*T1**(-1)*MS2*MG2 - 
     +    32*TG**(-2)*S*MG2 + 48*TG**(-1)*UG**(-1)*S*T**(-1)*NS*MS2*MG2
     +     - 48*TG**(-1)*UG**(-1)*S*T**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    UG**(-1)*S*NS*MG2 + 16*TG**(-1)*UG**(-1)*S*MG2 - 16*TG**(-1)*
     +    UG**(-1)*S**2*T**(-1)*NS*MG2 + 16*TG**(-1)*UG**(-1)*S**2*
     +    T**(-1)*MG2 - 16*TG**(-1)*S*T**(-1)*NS*MG2 + 16*TG**(-1)*S*
     +    T**(-1)*MG2 - 32*UG**(-2)*S*T1**(-1)*MS2*MG2 + 32*UG**(-2)*S*
     +    T1**(-1)*MG2**2 + 32*UG**(-2)*S*MG2 + 32*UG**(-2)*S**2*
     +    T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N*CF**2 * (  - 16*UG**(-1)*S*
     +    T**(-1)*NS*MG2 + 16*UG**(-1)*S*T**(-1)*MG2 + 32*UG**(-1)*S*
     +    T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N**2*CF * ( 24*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*T1*NS*MG2**2 - 24*TG**(-2)*UG**(-1)*S*T**(-1)*T1*
     +    MG2**2 - 8*TG**(-2)*UG**(-1)*S*T**(-1)*NS*MG2**3 + 8*TG**(-2)
     +    *UG**(-1)*S*T**(-1)*MG2**3 + 8*TG**(-2)*UG**(-1)*S*NS*MG2**2
     +     - 8*TG**(-2)*UG**(-1)*S*MG2**2 + 8*TG**(-2)*UG**(-1)*S**2*
     +    T**(-1)*NS*MG2**2 - 8*TG**(-2)*UG**(-1)*S**2*T**(-1)*MG2**2
     +     + 8*TG**(-2)*S*T**(-1)*NS*MG2**2 - 8*TG**(-2)*S*T**(-1)*
     +    MG2**2 + 16*TG**(-2)*S*T1**(-1)*MS2*MG2 + 16*TG**(-2)*S*MG2
     +     - 24*TG**(-1)*UG**(-1)*S*T**(-1)*NS*MS2*MG2 + 24*TG**(-1)*
     +    UG**(-1)*S*T**(-1)*MS2*MG2 + 8*TG**(-1)*UG**(-1)*S*NS*MG2 - 8
     +    *TG**(-1)*UG**(-1)*S*MG2 + 8*TG**(-1)*UG**(-1)*S**2*T**(-1)*
     +    NS*MG2 - 8*TG**(-1)*UG**(-1)*S**2*T**(-1)*MG2 + 8*TG**(-1)*S*
     +    T**(-1)*NS*MG2 - 8*TG**(-1)*S*T**(-1)*MG2 + 16*UG**(-2)*S*
     +    T1**(-1)*MS2*MG2 - 16*UG**(-2)*S*T1**(-1)*MG2**2 - 16*
     +    UG**(-2)*S*MG2 - 16*UG**(-2)*S**2*T1**(-1)*MG2 + 8*UG**(-1)*S
     +    *T**(-1)*NS*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*N**2*CF * (  - 8*UG**(-1)*S*T**(-1)
     +    *MG2 - 16*UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,1)*CF**2 * (  - 32*TG**(-1)*UG**(-2)*S
     +    *T**(-1)*T1*NS*MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*T1*
     +    MG2**2 - 32*TG**(-1)*UG**(-2)*S*T**(-1)*NS*MS2*MG2**2 + 32*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*NS*MG2**3 + 32*TG**(-1)*UG**(-2)*
     +    S*T**(-1)*MS2*MG2**2 - 32*TG**(-1)*UG**(-2)*S*T**(-1)*MG2**3
     +     + 32*TG**(-1)*UG**(-1)*S*T1**(-1)*MG2**2 - 64*TG**(-1)*
     +    UG**(-1)*S*MG2 - 32*TG**(-1)*UG**(-1)*S**2*T1**(-1)*MG2 - 32*
     +    TG**(-1)*S*T1**(-1)*MG2 + 32*UG**(-2)*S*T**(-1)*NS*MG2**2 - 
     +    32*UG**(-2)*S*T**(-1)*MG2**2 + 64*UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N*CF * (  - 16*TG**(-1)*UG**(-1)*S*
     +    U1**(-1)*MG2**2 - 16*TG**(-1)*S*U1**(-1)*MG2 + 8*UG**(-3)*S*
     +    U**(-1)*T1*NS*MG2**2 - 8*UG**(-3)*S*U**(-1)*T1*MG2**2 + 8*
     +    UG**(-3)*S*U**(-1)*NS*MG2**3 - 8*UG**(-3)*S*U**(-1)*MG2**3 - 
     +    8*UG**(-3)*S*NS*MG2**2 + 8*UG**(-3)*S*MG2**2 + 8*UG**(-3)*
     +    S**2*U**(-1)*NS*MG2**2 - 8*UG**(-3)*S**2*U**(-1)*MG2**2 + 8*
     +    UG**(-2)*S*U**(-1)*T1*NS*MG2 - 8*UG**(-2)*S*U**(-1)*T1*MG2 + 
     +    16*UG**(-2)*S*U**(-1)*NS*MG2**2 - 16*UG**(-2)*S*U**(-1)*
     +    MG2**2 - 8*UG**(-2)*S*NS*MG2 + 8*UG**(-2)*S*MG2 + 8*UG**(-2)*
     +    S**2*U**(-1)*NS*MG2 - 8*UG**(-2)*S**2*U**(-1)*MG2 + 8*
     +    UG**(-1)*S*U**(-1)*NS*MG2 - 8*UG**(-1)*S*U**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N*CF**2 * (  - 16*TG**(-1)*UG**(-2)
     +    *S*U**(-1)*NS*MS2*MG2**2 + 32*TG**(-1)*UG**(-2)*S*U**(-1)*NS*
     +    MG2**3 + 16*TG**(-1)*UG**(-2)*S*U**(-1)*MS2*MG2**2 - 32*
     +    TG**(-1)*UG**(-2)*S*U**(-1)*MG2**3 - 32*TG**(-1)*UG**(-2)*S*
     +    U1**(-1)*MS2*MG2**2 + 32*TG**(-1)*UG**(-2)*S*U1**(-1)*MG2**3
     +     - 16*TG**(-1)*UG**(-2)*S*NS*MG2**2 + 48*TG**(-1)*UG**(-2)*S*
     +    MG2**2 + 16*TG**(-1)*UG**(-2)*S**2*U**(-1)*NS*MG2**2 - 16*
     +    TG**(-1)*UG**(-2)*S**2*U**(-1)*MG2**2 + 32*TG**(-1)*UG**(-2)*
     +    S**2*U1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MS2*
     +    MG2 + 48*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MG2**2 + 16*TG**(-1)*
     +    UG**(-1)*S*U**(-1)*MS2*MG2 - 48*TG**(-1)*UG**(-1)*S*U**(-1)*
     +    MG2**2 - 32*TG**(-1)*UG**(-1)*S*U1**(-1)*MS2*MG2 + 32*
     +    TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)*S*
     +    NS*MG2 + 48*TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*
     +    S**2*U**(-1)*NS*MG2 - 16*TG**(-1)*UG**(-1)*S**2*U**(-1)*MG2
     +     + 32*TG**(-1)*UG**(-1)*S**2*U1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N*CF**2 * ( 16*TG**(-1)*S*U**(-1)*
     +    NS*MG2 - 16*TG**(-1)*S*U**(-1)*MG2 + 16*UG**(-2)*S*U**(-1)*NS
     +    *MG2**2 - 16*UG**(-2)*S*U**(-1)*MG2**2 + 16*UG**(-1)*S*
     +    U**(-1)*NS*MG2 - 16*UG**(-1)*S*U**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N**2*CF * ( 8*TG**(-1)*UG**(-2)*S*
     +    U**(-1)*NS*MS2*MG2**2 - 16*TG**(-1)*UG**(-2)*S*U**(-1)*NS*
     +    MG2**3 - 8*TG**(-1)*UG**(-2)*S*U**(-1)*MS2*MG2**2 + 16*
     +    TG**(-1)*UG**(-2)*S*U**(-1)*MG2**3 + 16*TG**(-1)*UG**(-2)*S*
     +    U1**(-1)*MS2*MG2**2 - 16*TG**(-1)*UG**(-2)*S*U1**(-1)*MG2**3
     +     + 8*TG**(-1)*UG**(-2)*S*NS*MG2**2 - 24*TG**(-1)*UG**(-2)*S*
     +    MG2**2 - 8*TG**(-1)*UG**(-2)*S**2*U**(-1)*NS*MG2**2 + 8*
     +    TG**(-1)*UG**(-2)*S**2*U**(-1)*MG2**2 - 16*TG**(-1)*UG**(-2)*
     +    S**2*U1**(-1)*MG2**2 + 8*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MS2*
     +    MG2 - 24*TG**(-1)*UG**(-1)*S*U**(-1)*NS*MG2**2 - 8*TG**(-1)*
     +    UG**(-1)*S*U**(-1)*MS2*MG2 + 24*TG**(-1)*UG**(-1)*S*U**(-1)*
     +    MG2**2 + 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MS2*MG2 - 16*
     +    TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 + 8*TG**(-1)*UG**(-1)*S*
     +    NS*MG2 - 24*TG**(-1)*UG**(-1)*S*MG2 - 8*TG**(-1)*UG**(-1)*
     +    S**2*U**(-1)*NS*MG2 + 8*TG**(-1)*UG**(-1)*S**2*U**(-1)*MG2 - 
     +    16*TG**(-1)*UG**(-1)*S**2*U1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*N**2*CF * (  - 8*TG**(-1)*S*U**(-1)
     +    *NS*MG2 + 8*TG**(-1)*S*U**(-1)*MG2 - 8*UG**(-2)*S*U**(-1)*NS*
     +    MG2**2 + 8*UG**(-2)*S*U**(-1)*MG2**2 - 8*UG**(-1)*S*U**(-1)*
     +    NS*MG2 + 8*UG**(-1)*S*U**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(1,2)*CF**2 * ( 32*TG**(-1)*UG**(-1)*S*
     +    U1**(-1)*MG2**2 + 32*TG**(-1)*S*U1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(2,1)*N*CF * ( 8*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*MG2**3 - 24*TG**(-2)*UG**(-1)*S*MG2**2 + 8*TG**(-1)*
     +    UG**(-1)*S*T**(-1)*MG2**2 + 16*TG**(-1)*UG**(-1)*S*T1**(-1)*
     +    MG2**2 - 56*TG**(-1)*UG**(-1)*S*MG2 - 16*TG**(-1)*UG**(-1)*
     +    S**2*T1**(-1)*MG2 - 16*TG**(-1)*S*T1**(-1)*MG2 + 32*UG**(-1)*
     +    S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(2,1)*N**2*CF * (  - 8*TG**(-3)*S*T**(-1)
     +    *MG2**3 + 24*TG**(-3)*S*MG2**2 - 8*TG**(-2)*S*T**(-1)*MG2**2
     +     - 16*TG**(-2)*S*T1**(-1)*MS2*MG2 + 8*TG**(-2)*S*MG2 - 16*
     +    UG**(-2)*S*T1**(-1)*MS2*MG2 + 16*UG**(-2)*S*T1**(-1)*MG2**2
     +     + 16*UG**(-2)*S*MG2 + 16*UG**(-2)*S**2*T1**(-1)*MG2 + 16*
     +    UG**(-1)*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(2,2)*N*CF * ( 8*TG**(-1)*UG**(-2)*S*
     +    U**(-1)*MG2**3 - 24*TG**(-1)*UG**(-2)*S*MG2**2 + 8*TG**(-1)*
     +    UG**(-1)*S*U**(-1)*MG2**2 + 16*TG**(-1)*UG**(-1)*S*U1**(-1)*
     +    MG2**2 - 24*TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*S*U1**(-1)*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(2,2)*N**2*CF * (  - 16*TG**(-1)*UG**(-2)
     +    *S*U1**(-1)*MS2*MG2**2 + 16*TG**(-1)*UG**(-2)*S*U1**(-1)*
     +    MG2**3 + 16*TG**(-1)*UG**(-2)*S*MG2**2 + 16*TG**(-1)*UG**(-2)
     +    *S**2*U1**(-1)*MG2**2 - 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MS2*
     +    MG2 + 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**2 + 16*TG**(-1)*
     +    UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*S**2*U1**(-1)*MG2 - 8*
     +    UG**(-3)*S*U**(-1)*MG2**3 + 24*UG**(-3)*S*MG2**2 - 8*UG**(-2)
     +    *S*U**(-1)*MG2**2 + 24*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N**(-1)*CF * (  - 32*TG**(-1)*
     +    UG**(-2)*S*T**(-1)*T1*MG2**2 - 32*TG**(-1)*UG**(-2)*S*T**(-1)
     +    *MS2*MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*MG2**3 + 32*
     +    UG**(-2)*S*T**(-1)*MG2**2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*T1
     +    *MG2**2 - 8*TG**(-3)*S*T**(-1)*MG2**2*MT2 + 8*TG**(-3)*S*
     +    T**(-1)*MG2**3 - 8*TG**(-3)*S*MG2**2 + 8*TG**(-2)*S*T**(-1)*
     +    MS2*MG2 - 8*TG**(-2)*S*T**(-1)*MG2*MT2 + 8*TG**(-2)*S*T**(-1)
     +    *MG2**2 - 8*TG**(-2)*S*MG2 + 8*TG**(-1)*UG**(-2)*S*T**(-1)*T1
     +    *MG2**2 + 24*TG**(-1)*UG**(-2)*S*T**(-1)*MS2*MG2**2 - 24*
     +    TG**(-1)*UG**(-2)*S*T**(-1)*MG2**3 - 8*TG**(-1)*UG**(-2)*S**2
     +    *T**(-1)*MG2**2 - 8*TG**(-1)*UG**(-1)*S*T**(-1)*MG2**2 - 8*
     +    UG**(-2)*S*T**(-1)*T1*MG2 + 8*UG**(-2)*S*T**(-1)*MS2*MG2 - 24
     +    *UG**(-2)*S*T**(-1)*MG2**2 - 8*UG**(-2)*S**2*T**(-1)*MG2 - 8*
     +    UG**(-1)*S*T**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N*CF**2 * (  - 48*TG**(-2)*UG**(-1)
     +    *S*T**(-1)*T1*MG2**2 - 16*TG**(-2)*UG**(-1)*S*T**(-1)*MG2**2*
     +    MT2 + 16*TG**(-2)*UG**(-1)*S*T**(-1)*MG2**3 - 16*TG**(-2)*
     +    UG**(-1)*S*MG2**2 - 16*TG**(-2)*UG**(-1)*S**2*T**(-1)*MG2**2
     +     - 16*TG**(-2)*S*T**(-1)*MG2**2 + 48*TG**(-1)*UG**(-1)*S*
     +    T**(-1)*MS2*MG2 - 16*TG**(-1)*UG**(-1)*S*T**(-1)*MG2*MT2 - 16
     +    *TG**(-1)*UG**(-1)*S*MG2 - 16*TG**(-1)*UG**(-1)*S**2*T**(-1)*
     +    MG2 - 16*TG**(-1)*S*T**(-1)*MG2 - 16*UG**(-1)*S*T**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*N**2*CF * ( 24*TG**(-2)*UG**(-1)*S*
     +    T**(-1)*T1*MG2**2 + 8*TG**(-2)*UG**(-1)*S*T**(-1)*MG2**2*MT2
     +     - 8*TG**(-2)*UG**(-1)*S*T**(-1)*MG2**3 + 8*TG**(-2)*UG**(-1)
     +    *S*MG2**2 + 8*TG**(-2)*UG**(-1)*S**2*T**(-1)*MG2**2 + 8*
     +    TG**(-2)*S*T**(-1)*MG2**2 - 24*TG**(-1)*UG**(-1)*S*T**(-1)*
     +    MS2*MG2 + 8*TG**(-1)*UG**(-1)*S*T**(-1)*MG2*MT2 + 8*TG**(-1)*
     +    UG**(-1)*S*MG2 + 8*TG**(-1)*UG**(-1)*S**2*T**(-1)*MG2 + 8*
     +    TG**(-1)*S*T**(-1)*MG2 + 8*UG**(-1)*S*T**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,1)*CF**2 * (  - 32*TG**(-1)*UG**(-2)*S
     +    *T**(-1)*T1*MG2**2 - 32*TG**(-1)*UG**(-2)*S*T**(-1)*MS2*
     +    MG2**2 + 32*TG**(-1)*UG**(-2)*S*T**(-1)*MG2**3 + 32*UG**(-2)*
     +    S*T**(-1)*MG2**2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N*CF * ( 8*UG**(-3)*S*U**(-1)*T1*
     +    MG2**2 - 8*UG**(-3)*S*U**(-1)*MG2**2*MT2 + 8*UG**(-3)*S*
     +    U**(-1)*MG2**3 - 8*UG**(-3)*S*MG2**2 + 8*UG**(-3)*S**2*
     +    U**(-1)*MG2**2 + 8*UG**(-2)*S*U**(-1)*T1*MG2 - 8*UG**(-2)*S*
     +    U**(-1)*MG2*MT2 + 16*UG**(-2)*S*U**(-1)*MG2**2 - 8*UG**(-2)*S
     +    *MG2 + 8*UG**(-2)*S**2*U**(-1)*MG2 + 8*UG**(-1)*S*U**(-1)*MG2
     +     )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N*CF**2 * (  - 16*TG**(-1)*UG**(-2)
     +    *S*U**(-1)*MS2*MG2**2 - 16*TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2
     +    *MT2 + 32*TG**(-1)*UG**(-2)*S*U**(-1)*MG2**3 - 16*TG**(-1)*
     +    UG**(-2)*S*MG2**2 + 16*TG**(-1)*UG**(-2)*S**2*U**(-1)*MG2**2
     +     - 16*TG**(-1)*UG**(-1)*S*U**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    UG**(-1)*S*U**(-1)*MG2*MT2 + 48*TG**(-1)*UG**(-1)*S*U**(-1)*
     +    MG2**2 - 16*TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*
     +    S**2*U**(-1)*MG2 + 16*TG**(-1)*S*U**(-1)*MG2 + 16*UG**(-2)*S*
     +    U**(-1)*MG2**2 + 16*UG**(-1)*S*U**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0D(3,2)*N**2*CF * ( 8*TG**(-1)*UG**(-2)*S*
     +    U**(-1)*MS2*MG2**2 + 8*TG**(-1)*UG**(-2)*S*U**(-1)*MG2**2*MT2
     +     - 16*TG**(-1)*UG**(-2)*S*U**(-1)*MG2**3 + 8*TG**(-1)*
     +    UG**(-2)*S*MG2**2 - 8*TG**(-1)*UG**(-2)*S**2*U**(-1)*MG2**2
     +     + 8*TG**(-1)*UG**(-1)*S*U**(-1)*MS2*MG2 + 8*TG**(-1)*
     +    UG**(-1)*S*U**(-1)*MG2*MT2 - 24*TG**(-1)*UG**(-1)*S*U**(-1)*
     +    MG2**2 + 8*TG**(-1)*UG**(-1)*S*MG2 - 8*TG**(-1)*UG**(-1)*S**2
     +    *U**(-1)*MG2 - 8*TG**(-1)*S*U**(-1)*MG2 - 8*UG**(-2)*S*
     +    U**(-1)*MG2**2 - 8*UG**(-1)*S*U**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0E(1)*N*CF * ( 16*TG**(-2)*UG**(-1)*S*
     +    MG2**2 + 16*TG**(-1)*UG**(-2)*S*MG2**2 + 16*TG**(-1)*UG**(-1)
     +    *S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0E(1)*N**2*CF * (  - 16*TG**(-3)*S*MG2**2
     +     - 8*TG**(-2)*S*MG2 - 16*UG**(-3)*S*MG2**2 - 8*UG**(-2)*S*MG2
     +     )
     +
      MQQLLV = MQQLLV + SK1B0E(2)*N**(-1)*CF * ( 32*TG**(-1)*UG**(-2)*S
     +    *T1*NS*MG2 - 32*TG**(-1)*UG**(-2)*S*T1*MG2 + 32*TG**(-1)*
     +    UG**(-2)*S*NS*MS2*MG2 - 32*TG**(-1)*UG**(-2)*S*NS*MG2**2 - 32
     +    *TG**(-1)*UG**(-2)*S*MS2*MG2 + 32*TG**(-1)*UG**(-2)*S*MG2**2
     +     - 32*UG**(-2)*S*NS*MG2 + 32*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0E(2)*N*CF * ( 8*TG**(-3)*S*T1*NS*MG2 - 8*
     +    TG**(-3)*S*T1*MG2 + 4*TG**(-2)*S*T1*NS - 4*TG**(-2)*S*T1 - 8*
     +    TG**(-2)*S*NS*MG2 + 8*TG**(-2)*S*MG2 - 8*TG**(-1)*UG**(-2)*S*
     +    T1*NS*MG2 + 8*TG**(-1)*UG**(-2)*S*T1*MG2 - 24*TG**(-1)*
     +    UG**(-2)*S*NS*MS2*MG2 + 24*TG**(-1)*UG**(-2)*S*NS*MG2**2 + 24
     +    *TG**(-1)*UG**(-2)*S*MS2*MG2 - 24*TG**(-1)*UG**(-2)*S*MG2**2
     +     + 8*TG**(-1)*UG**(-2)*S**2*NS*MG2 - 8*TG**(-1)*UG**(-2)*S**2
     +    *MG2 + 8*TG**(-1)*UG**(-1)*S*NS*MG2 - 8*TG**(-1)*UG**(-1)*S*
     +    MG2 - 4*TG**(-1)*S*NS + 4*TG**(-1)*S - 8*UG**(-3)*S*T1*NS*MG2
     +     + 8*UG**(-3)*S*T1*MG2 - 8*UG**(-3)*S**2*NS*MG2 + 8*UG**(-3)*
     +    S**2*MG2 - 4*UG**(-2)*S*T1*NS + 4*UG**(-2)*S*T1 + 8*UG**(-2)*
     +    S*NS*MG2 - 8*UG**(-2)*S*MG2 - 4*UG**(-2)*S**2*NS + 4*UG**(-2)
     +    *S**2 - 4*UG**(-1)*S*NS + 4*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(2)*N*CF**2 * ( 48*TG**(-2)*UG**(-1)*S*T1
     +    *NS*MG2 - 48*TG**(-2)*UG**(-1)*S*T1*MG2 + 16*TG**(-2)*
     +    UG**(-1)*S**2*NS*MG2 - 16*TG**(-2)*UG**(-1)*S**2*MG2 + 16*
     +    TG**(-2)*S*NS*MG2 - 16*TG**(-2)*S*MG2 + 16*TG**(-1)*UG**(-2)*
     +    S*NS*MS2*MG2 - 16*TG**(-1)*UG**(-2)*S*NS*MG2**2 - 16*TG**(-1)
     +    *UG**(-2)*S*MS2*MG2 + 16*TG**(-1)*UG**(-2)*S*MG2**2 - 16*
     +    TG**(-1)*UG**(-2)*S**2*NS*MG2 + 16*TG**(-1)*UG**(-2)*S**2*MG2
     +     + 16*TG**(-1)*UG**(-1)*S*NS*MS2 - 64*TG**(-1)*UG**(-1)*S*NS*
     +    MG2 - 16*TG**(-1)*UG**(-1)*S*MS2 + 64*TG**(-1)*UG**(-1)*S*MG2
     +     - 16*TG**(-1)*UG**(-1)*S**2*NS + 16*TG**(-1)*UG**(-1)*S**2
     +     - 16*TG**(-1)*S*NS + 16*TG**(-1)*S - 16*UG**(-2)*S*NS*MG2 + 
     +    16*UG**(-2)*S*MG2 - 16*UG**(-1)*S*NS + 16*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(2)*N**2*CF * (  - 24*TG**(-2)*UG**(-1)*S
     +    *T1*NS*MG2 + 24*TG**(-2)*UG**(-1)*S*T1*MG2 - 8*TG**(-2)*
     +    UG**(-1)*S**2*NS*MG2 + 8*TG**(-2)*UG**(-1)*S**2*MG2 - 8*
     +    TG**(-2)*S*NS*MG2 + 8*TG**(-2)*S*MG2 - 8*TG**(-1)*UG**(-2)*S*
     +    NS*MS2*MG2 + 8*TG**(-1)*UG**(-2)*S*NS*MG2**2 + 8*TG**(-1)*
     +    UG**(-2)*S*MS2*MG2 - 8*TG**(-1)*UG**(-2)*S*MG2**2 + 8*
     +    TG**(-1)*UG**(-2)*S**2*NS*MG2 - 8*TG**(-1)*UG**(-2)*S**2*MG2
     +     - 8*TG**(-1)*UG**(-1)*S*NS*MS2 + 32*TG**(-1)*UG**(-1)*S*NS*
     +    MG2 + 8*TG**(-1)*UG**(-1)*S*MS2 - 32*TG**(-1)*UG**(-1)*S*MG2
     +     + 8*TG**(-1)*UG**(-1)*S**2*NS - 8*TG**(-1)*UG**(-1)*S**2 + 8
     +    *TG**(-1)*S*NS - 8*TG**(-1)*S + 8*UG**(-2)*S*NS*MG2 - 8*
     +    UG**(-2)*S*MG2 + 8*UG**(-1)*S*NS - 8*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(2)*CF**2 * ( 32*TG**(-1)*UG**(-2)*S*T1*
     +    NS*MG2 - 32*TG**(-1)*UG**(-2)*S*T1*MG2 + 32*TG**(-1)*UG**(-2)
     +    *S*NS*MS2*MG2 - 32*TG**(-1)*UG**(-2)*S*NS*MG2**2 - 32*
     +    TG**(-1)*UG**(-2)*S*MS2*MG2 + 32*TG**(-1)*UG**(-2)*S*MG2**2
     +     - 32*UG**(-2)*S*NS*MG2 + 32*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N**(-1)*CF * ( 32*TG**(-1)*UG**(-2)*S
     +    *T1*MG2 + 32*TG**(-1)*UG**(-2)*S*MS2*MG2 - 32*TG**(-1)*
     +    UG**(-2)*S*MG2**2 - 32*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N*CF * ( 8*TG**(-3)*S*T1*MG2 + 8*
     +    TG**(-3)*S*MG2*MT2 + 4*TG**(-2)*S*T1 - 8*TG**(-2)*S*MG2 + 4*
     +    TG**(-2)*S*MT2 - 8*TG**(-1)*UG**(-2)*S*T1*MG2 - 24*TG**(-1)*
     +    UG**(-2)*S*MS2*MG2 + 24*TG**(-1)*UG**(-2)*S*MG2**2 + 8*
     +    TG**(-1)*UG**(-2)*S**2*MG2 + 8*TG**(-1)*UG**(-1)*S*MG2 - 4*
     +    TG**(-1)*S - 8*UG**(-3)*S*T1*MG2 + 8*UG**(-3)*S*MG2*MT2 - 8*
     +    UG**(-3)*S**2*MG2 - 4*UG**(-2)*S*T1 + 8*UG**(-2)*S*MG2 + 4*
     +    UG**(-2)*S*MT2 - 4*UG**(-2)*S**2 - 4*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N*CF**2 * ( 48*TG**(-2)*UG**(-1)*S*T1
     +    *MG2 + 16*TG**(-2)*UG**(-1)*S*MG2*MT2 + 16*TG**(-2)*UG**(-1)*
     +    S**2*MG2 + 16*TG**(-2)*S*MG2 + 16*TG**(-1)*UG**(-2)*S*MS2*MG2
     +     + 16*TG**(-1)*UG**(-2)*S*MG2*MT2 - 16*TG**(-1)*UG**(-2)*S*
     +    MG2**2 - 16*TG**(-1)*UG**(-2)*S**2*MG2 + 16*TG**(-1)*UG**(-1)
     +    *S*MS2 - 64*TG**(-1)*UG**(-1)*S*MG2 + 16*TG**(-1)*UG**(-1)*S*
     +    MT2 - 16*TG**(-1)*UG**(-1)*S**2 - 16*TG**(-1)*S - 16*UG**(-2)
     +    *S*MG2 - 16*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*N**2*CF * (  - 24*TG**(-2)*UG**(-1)*S
     +    *T1*MG2 - 8*TG**(-2)*UG**(-1)*S*MG2*MT2 - 8*TG**(-2)*UG**(-1)
     +    *S**2*MG2 - 8*TG**(-2)*S*MG2 - 8*TG**(-1)*UG**(-2)*S*MS2*MG2
     +     - 8*TG**(-1)*UG**(-2)*S*MG2*MT2 + 8*TG**(-1)*UG**(-2)*S*
     +    MG2**2 + 8*TG**(-1)*UG**(-2)*S**2*MG2 - 8*TG**(-1)*UG**(-1)*S
     +    *MS2 + 32*TG**(-1)*UG**(-1)*S*MG2 - 8*TG**(-1)*UG**(-1)*S*MT2
     +     + 8*TG**(-1)*UG**(-1)*S**2 + 8*TG**(-1)*S + 8*UG**(-2)*S*MG2
     +     + 8*UG**(-1)*S )
     +
      MQQLLV = MQQLLV + SK1B0E(3)*CF**2 * ( 32*TG**(-1)*UG**(-2)*S*T1*
     +    MG2 + 32*TG**(-1)*UG**(-2)*S*MS2*MG2 - 32*TG**(-1)*UG**(-2)*S
     +    *MG2**2 - 32*UG**(-2)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1BP(1)*N*CF**2 * ( 16*TG**(-2)*S*MS2*MG2 + 16
     +    *UG**(-2)*S*MS2*MG2 )
     +
      MQQLLV = MQQLLV + SK1BP(1)*CF**2 * (  - 32*TG**(-1)*UG**(-1)*S*
     +    MS2*MG2 )
     +
      MQQLLV = MQQLLV + SK1BP(2)*N*CF**2 * (  - 8*TG**(-2)*S*MS2*MG2 + 
     +    8*TG**(-2)*S*MG2**2 - 8*UG**(-2)*S*MS2*MG2 + 8*UG**(-2)*S*
     +    MG2**2 )
     +
      MQQLLV = MQQLLV + SK1BP(2)*CF**2 * ( 16*TG**(-1)*UG**(-1)*S*MS2*
     +    MG2 - 16*TG**(-1)*UG**(-1)*S*MG2**2 )
     +
      MQQLLV = MQQLLV + SK1BP(3)*N*CF**2 * ( 4*TG**(-2)*S*MS2*MG2 - 4*
     +    TG**(-2)*S*MG2**2 + 4*UG**(-2)*S*MS2*MG2 - 4*UG**(-2)*S*
     +    MG2**2 )
     +
      MQQLLV = MQQLLV + SK1BP(3)*CF**2 * (  - 8*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 + 8*TG**(-1)*UG**(-1)*S*MG2**2 )
     +
      MQQLLV = MQQLLV + SK1C0A(5)*N*CF * ( 8*TG**(-1)*S*MG2 + 8*
     +    UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0A(5)*N*CF**2 * ( 16*TG**(-1)*S*MG2 + 16*
     +    UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0A(5)*N**2*CF * (  - 8*TG**(-1)*S*MG2 - 8*
     +    UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0A(5)*CF**2 * (  - 12*TG**(-1)*S*MG2 - 12*
     +    UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N*CF * ( 24*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 - 24*TG**(-1)*UG**(-1)*S*MG2**2 - 8*TG**(-1)*UG**(-1)*
     +    S**2*MG2 - 8*TG**(-1)*S*MG2 + 4*TG**(-1)*T1*MG2 - 16*UG**(-1)
     +    *S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N*CF**2 * ( 16*TG**(-2)*S*MS2*MG2
     +     - 16*TG**(-2)*S*MG2**2 - 16*TG**(-1)*S*MG2 - 16*UG**(-2)*S*
     +    T1*MG2 + 16*UG**(-2)*S*MS2*MG2 - 16*UG**(-2)*S*MG2**2 - 16*
     +    UG**(-2)*S**2*MG2 - 16*UG**(-1)*S*MG2 + 16*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*N**2*CF * (  - 8*TG**(-2)*S*MS2*MG2
     +     + 8*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 + 8*UG**(-2)*S*T1*
     +    MG2 - 8*UG**(-2)*S*MS2*MG2 + 8*UG**(-2)*S*MG2**2 + 8*UG**(-2)
     +    *S**2*MG2 + 8*UG**(-1)*S*MG2 - 4*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,1)*CF**2 * (  - 48*TG**(-1)*UG**(-1)*S
     +    *MS2*MG2 + 48*TG**(-1)*UG**(-1)*S*MG2**2 + 16*TG**(-1)*
     +    UG**(-1)*S**2*MG2 + 16*TG**(-1)*S*MG2 - 8*TG**(-1)*T1*MG2 + 
     +    32*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N*CF * ( 8*TG**(-1)*UG**(-1)*S*MS2*
     +    MG2 - 8*TG**(-1)*UG**(-1)*S*MG2**2 - 8*TG**(-1)*S*MG2 - 4*
     +    UG**(-1)*S*MG2 - 4*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N*CF**2 * (  - 16*TG**(-1)*S*MG2 - 
     +    16*TG**(-1)*T1*MG2 + 16*UG**(-2)*S*T1*MG2 + 16*UG**(-2)*S**2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*N**2*CF * ( 4*TG**(-1)*S*MG2 + 4*
     +    TG**(-1)*T1*MG2 - 8*UG**(-2)*S*T1*MG2 - 8*UG**(-2)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(3,2)*CF**2 * (  - 16*TG**(-1)*UG**(-1)*S
     +    *MS2*MG2 + 16*TG**(-1)*UG**(-1)*S*MG2**2 + 16*TG**(-1)*S*MG2
     +     + 8*UG**(-1)*S*MG2 + 8*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,1)*N*CF * ( 24*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 - 24*TG**(-1)*UG**(-1)*S*MG2**2 - 8*TG**(-1)*UG**(-1)*
     +    S**2*MG2 - 8*TG**(-1)*S*MG2 - 16*UG**(-1)*S*MG2 + 4*UG**(-1)*
     +    T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,1)*N*CF**2 * ( 16*TG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,1)*N**2*CF * (  - 8*TG**(-2)*S*MS2*MG2
     +     + 8*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 - 4*TG**(-1)*T1*MG2
     +     + 8*UG**(-2)*S*T1*MG2 - 8*UG**(-2)*S*MS2*MG2 + 8*UG**(-2)*S*
     +    MG2**2 + 8*UG**(-2)*S**2*MG2 + 8*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,1)*CF**2 * (  - 8*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,2)*N*CF * (  - 8*TG**(-1)*UG**(-1)*S*
     +    MS2*MG2 + 8*TG**(-1)*UG**(-1)*S*MG2**2 + 8*TG**(-1)*UG**(-1)*
     +    S**2*MG2 - 4*TG**(-1)*S*MG2 - 4*TG**(-1)*T1*MG2 + 8*UG**(-1)*
     +    S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,2)*N*CF**2 * (  - 16*UG**(-1)*S*MG2 - 
     +    16*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(4,2)*N**2*CF * (  - 8*UG**(-2)*S*T1*MG2
     +     - 8*UG**(-2)*S**2*MG2 + 4*UG**(-1)*S*MG2 + 4*UG**(-1)*T1*MG2
     +     )
     +
      MQQLLV = MQQLLV + SK1C0C(4,2)*CF**2 * ( 8*TG**(-1)*S*MG2 + 8*
     +    TG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,1)*N*CF * (  - 8*TG**(-1)*UG**(-1)*S*
     +    MS2*MG2 - 24*TG**(-1)*UG**(-1)*S*MG2**2 - 32*TG**(-1)*
     +    UG**(-1)*S**2*T1**(-1)*MS2*MG2 + 24*TG**(-1)*UG**(-1)*S**2*
     +    T1**(-1)*MG2**2 - 48*TG**(-1)*UG**(-1)*S**2*MG2 - 8*TG**(-1)*
     +    UG**(-1)*S**3*T1**(-1)*MG2 - 24*TG**(-1)*S*T1**(-1)*MS2*MG2
     +     + 16*TG**(-1)*S*T1**(-1)*MG2**2 - 24*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*S**2*T1**(-1)*MG2 + 16*UG**(-1)*S*T1**(-1)*MS2*MG2
     +     - 16*UG**(-1)*S*MG2 + 32*UG**(-1)*S**2*T1**(-1)*MG2 + 4*
     +    UG**(-1)*T1*MG2 + 16*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,1)*N*CF**2 * ( 16*TG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,1)*N**2*CF * ( 8*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 8*TG**(-2)*S*
     +    MG2**2 + 8*TG**(-1)*S*MG2 - 4*TG**(-1)*T1*MG2 + 8*UG**(-2)*S*
     +    T1**(-1)*MS2*MG2**2 - 8*UG**(-2)*S*T1**(-1)*MS2**2*MG2 + 8*
     +    UG**(-2)*S*T1*MG2 + 8*UG**(-2)*S*MG2**2 + 8*UG**(-2)*S**2*
     +    T1**(-1)*MG2**2 + 16*UG**(-2)*S**2*MG2 + 8*UG**(-2)*S**3*
     +    T1**(-1)*MG2 + 8*UG**(-1)*S*T1**(-1)*MG2**2 + 16*UG**(-1)*S*
     +    MG2 + 16*UG**(-1)*S**2*T1**(-1)*MG2 + 8*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,1)*CF**2 * (  - 8*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,2)*N*CF * ( 16*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 - 32*TG**(-1)*UG**(-1)*S*MG2**2 - 12*TG**(-1)*UG**(-1)*
     +    S**2*MG2 + 8*TG**(-1)*S*U1**(-1)*MS2*MG2 - 24*TG**(-1)*S*MG2
     +     - 4*TG**(-1)*T1*MG2 - 12*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,2)*N*CF**2 * (  - 16*UG**(-1)*S*MG2 - 
     +    16*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,2)*N**2*CF * (  - 8*TG**(-1)*UG**(-1)*
     +    S*U1**(-1)*MS2*MG2**2 + 8*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**3
     +     + 8*TG**(-1)*UG**(-1)*S*MG2**2 + 8*TG**(-1)*UG**(-1)*S**2*
     +    U1**(-1)*MG2**2 - 8*TG**(-1)*S*U1**(-1)*MS2*MG2 + 8*TG**(-1)*
     +    S*U1**(-1)*MG2**2 + 8*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2*
     +    U1**(-1)*MG2 + 8*UG**(-2)*S*T1*MG2 + 16*UG**(-2)*S*MG2**2 + 8
     +    *UG**(-2)*S**2*MG2 + 28*UG**(-1)*S*MG2 + 4*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(5,2)*CF**2 * ( 8*TG**(-1)*S*MG2 + 8*
     +    TG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,1)*N*CF * ( 16*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 - 48*TG**(-1)*UG**(-1)*S*MG2**2 - 32*TG**(-1)*UG**(-1)*
     +    S**2*T1**(-1)*MS2*MG2 + 24*TG**(-1)*UG**(-1)*S**2*T1**(-1)*
     +    MG2**2 - 56*TG**(-1)*UG**(-1)*S**2*MG2 - 8*TG**(-1)*UG**(-1)*
     +    S**3*T1**(-1)*MG2 - 24*TG**(-1)*S*T1**(-1)*MS2*MG2 + 16*
     +    TG**(-1)*S*T1**(-1)*MG2**2 - 32*TG**(-1)*S*MG2 - 8*TG**(-1)*
     +    S**2*T1**(-1)*MG2 + 4*TG**(-1)*T1*MG2 + 16*UG**(-1)*S*
     +    T1**(-1)*MS2*MG2 - 24*UG**(-1)*S*MG2 + 32*UG**(-1)*S**2*
     +    T1**(-1)*MG2 + 16*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,1)*N*CF**2 * (  - 16*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2**2 + 16*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 16*
     +    TG**(-2)*S*MS2*MG2 - 32*TG**(-2)*S*MG2**2 - 16*TG**(-1)*S*MG2
     +     - 16*UG**(-2)*S*T1**(-1)*MS2*MG2**2 + 16*UG**(-2)*S*T1**(-1)
     +    *MS2**2*MG2 - 32*UG**(-2)*S*T1*MG2 + 16*UG**(-2)*S*MS2*MG2 - 
     +    32*UG**(-2)*S*MG2**2 - 16*UG**(-2)*S**2*T1**(-1)*MG2**2 - 48*
     +    UG**(-2)*S**2*MG2 - 16*UG**(-2)*S**3*T1**(-1)*MG2 - 16*
     +    UG**(-1)*S*T1**(-1)*MG2**2 - 48*UG**(-1)*S*MG2 - 32*UG**(-1)*
     +    S**2*T1**(-1)*MG2 + 16*UG**(-1)*T1*MG2 - 16*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,1)*N**2*CF * ( 8*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 - 8*TG**(-2)*S*
     +    MS2*MG2 + 16*TG**(-2)*S*MG2**2 + 8*TG**(-1)*S*MG2 + 8*
     +    UG**(-2)*S*T1**(-1)*MS2*MG2**2 - 8*UG**(-2)*S*T1**(-1)*MS2**2
     +    *MG2 + 16*UG**(-2)*S*T1*MG2 - 8*UG**(-2)*S*MS2*MG2 + 16*
     +    UG**(-2)*S*MG2**2 + 8*UG**(-2)*S**2*T1**(-1)*MG2**2 + 24*
     +    UG**(-2)*S**2*MG2 + 8*UG**(-2)*S**3*T1**(-1)*MG2 + 8*UG**(-1)
     +    *S*T1**(-1)*MG2**2 + 24*UG**(-1)*S*MG2 + 16*UG**(-1)*S**2*
     +    T1**(-1)*MG2 - 4*UG**(-1)*T1*MG2 + 8*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,1)*CF**2 * (  - 32*TG**(-1)*UG**(-1)*S
     +    *MS2*MG2 + 96*TG**(-1)*UG**(-1)*S*MG2**2 + 64*TG**(-1)*
     +    UG**(-1)*S**2*T1**(-1)*MS2*MG2 - 48*TG**(-1)*UG**(-1)*S**2*
     +    T1**(-1)*MG2**2 + 112*TG**(-1)*UG**(-1)*S**2*MG2 + 16*
     +    TG**(-1)*UG**(-1)*S**3*T1**(-1)*MG2 + 48*TG**(-1)*S*T1**(-1)*
     +    MS2*MG2 - 32*TG**(-1)*S*T1**(-1)*MG2**2 + 64*TG**(-1)*S*MG2
     +     + 16*TG**(-1)*S**2*T1**(-1)*MG2 - 8*TG**(-1)*T1*MG2 - 32*
     +    UG**(-1)*S*T1**(-1)*MS2*MG2 + 48*UG**(-1)*S*MG2 - 64*UG**(-1)
     +    *S**2*T1**(-1)*MG2 - 32*S*T1**(-1)*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,2)*N*CF * ( 16*TG**(-1)*UG**(-1)*S*MS2
     +    *MG2 - 32*TG**(-1)*UG**(-1)*S*MG2**2 - 8*TG**(-1)*UG**(-1)*
     +    S**2*MG2 + 8*TG**(-1)*S*U1**(-1)*MS2*MG2 - 16*TG**(-1)*S*MG2
     +     - 12*UG**(-1)*S*MG2 - 4*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,2)*N*CF**2 * ( 16*TG**(-1)*UG**(-1)*S*
     +    U1**(-1)*MS2*MG2**2 - 16*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**3
     +     - 16*TG**(-1)*UG**(-1)*S*MG2**2 - 16*TG**(-1)*UG**(-1)*S**2*
     +    U1**(-1)*MG2**2 + 16*TG**(-1)*S*U1**(-1)*MS2*MG2 - 16*
     +    TG**(-1)*S*U1**(-1)*MG2**2 - 32*TG**(-1)*S*MG2 - 16*TG**(-1)*
     +    S**2*U1**(-1)*MG2 - 16*TG**(-1)*T1*MG2 - 32*UG**(-2)*S*MG2**2
     +     - 32*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,2)*N**2*CF * (  - 8*TG**(-1)*UG**(-1)*
     +    S*U1**(-1)*MS2*MG2**2 + 8*TG**(-1)*UG**(-1)*S*U1**(-1)*MG2**3
     +     + 8*TG**(-1)*UG**(-1)*S*MG2**2 + 8*TG**(-1)*UG**(-1)*S**2*
     +    U1**(-1)*MG2**2 - 8*TG**(-1)*S*U1**(-1)*MS2*MG2 + 8*TG**(-1)*
     +    S*U1**(-1)*MG2**2 + 12*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2*
     +    U1**(-1)*MG2 + 4*TG**(-1)*T1*MG2 + 16*UG**(-2)*S*MG2**2 + 16*
     +    UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1C0C(6,2)*CF**2 * (  - 32*TG**(-1)*UG**(-1)*S
     +    *MS2*MG2 + 64*TG**(-1)*UG**(-1)*S*MG2**2 + 16*TG**(-1)*
     +    UG**(-1)*S**2*MG2 - 16*TG**(-1)*S*U1**(-1)*MS2*MG2 + 32*
     +    TG**(-1)*S*MG2 + 24*UG**(-1)*S*MG2 + 8*UG**(-1)*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,1)*N*CF * ( 8*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,1)*N*CF**2 * ( 16*TG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,1)*N**2*CF * (  - 8*TG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,1)*CF**2 * (  - 8*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,2)*N*CF * ( 8*TG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,2)*N*CF**2 * ( 16*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,2)*N**2*CF * (  - 8*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(4,2)*CF**2 * (  - 8*TG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,1)*N*CF * (  - 16*UG**(-1)*S*T1*MG2 - 
     +    16*UG**(-1)*S*MG2**2 - 8*UG**(-1)*S**2*MG2 - 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,1)*N*CF**2 * ( 32*TG**(-1)*S*T1*MG2 - 
     +    32*TG**(-1)*S*MG2**2 + 16*TG**(-1)*S**2*MG2 - 32*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,1)*N**2*CF * (  - 16*TG**(-1)*S*T1*MG2
     +     + 16*TG**(-1)*S*MG2**2 - 8*TG**(-1)*S**2*MG2 + 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,1)*CF**2 * ( 16*UG**(-1)*S*T1*MG2 + 16*
     +    UG**(-1)*S*MG2**2 + 8*UG**(-1)*S**2*MG2 + 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,2)*N*CF * ( 16*TG**(-1)*S*T1*MG2 - 16*
     +    TG**(-1)*S*MG2**2 + 8*TG**(-1)*S**2*MG2 - 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,2)*N*CF**2 * (  - 32*UG**(-1)*S*T1*MG2
     +     - 32*UG**(-1)*S*MG2**2 - 16*UG**(-1)*S**2*MG2 - 32*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,2)*N**2*CF * ( 16*UG**(-1)*S*T1*MG2 + 
     +    16*UG**(-1)*S*MG2**2 + 8*UG**(-1)*S**2*MG2 + 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(7,2)*CF**2 * (  - 16*TG**(-1)*S*T1*MG2 + 
     +    16*TG**(-1)*S*MG2**2 - 8*TG**(-1)*S**2*MG2 + 16*S*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,1)*N*CF * ( 20*UG**(-1)*S*T1*MG2 + 12*
     +    UG**(-1)*S**2*MG2 + 8*UG**(-1)*T1**2*MG2 + 4*S*MG2 + 4*T1*MG2
     +     )
     +
      MQQLLV = MQQLLV + SK1D0(8,1)*N*CF**2 * ( 32*TG**(-1)*S*T1*MG2 + 
     +    32*TG**(-1)*S**2*MG2 + 16*S*MG2 + 16*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,1)*N**2*CF * (  - 8*TG**(-1)*S*T1*MG2
     +     - 8*TG**(-1)*S**2*MG2 - 4*S*MG2 - 4*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,1)*CF**2 * (  - 40*UG**(-1)*S*T1*MG2 - 
     +    24*UG**(-1)*S**2*MG2 - 16*UG**(-1)*T1**2*MG2 - 8*S*MG2 - 8*T1
     +    *MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,2)*N*CF * (  - 4*TG**(-1)*S*T1*MG2 + 8*
     +    TG**(-1)*T1**2*MG2 - 4*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,2)*N*CF**2 * (  - 32*UG**(-1)*S*T1*MG2
     +     - 16*T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,2)*N**2*CF * ( 8*UG**(-1)*S*T1*MG2 + 4*
     +    T1*MG2 )
     +
      MQQLLV = MQQLLV + SK1D0(8,2)*CF**2 * ( 8*TG**(-1)*S*T1*MG2 - 16*
     +    TG**(-1)*T1**2*MG2 + 8*T1*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N*CF * ( 128*TG**(-1)*UG**(-1)*S*MS2*
     +    MG2 - 64*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MG2 + 
     +    96*TG**(-2)*S*MS2*MG2 + 32*TG**(-2)*S*MG2**2 - 64*TG**(-2)*
     +    S**2*MG2 + 32*TG**(-1)*S*MG2 + 32*UG**(-2)*S*T1*MG2 + 96*
     +    UG**(-2)*S*MS2*MG2 + 32*UG**(-2)*S*MG2**2 - 32*UG**(-2)*S**2*
     +    MG2 + 32*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*N**2*CF * (  - 64*TG**(-2)*S*MS2*MG2 + 
     +    32*TG**(-2)*S**2*MG2 - 64*UG**(-2)*S*MS2*MG2 + 32*UG**(-2)*
     +    S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(1)*CF**2 * (  - 128*TG**(-1)*UG**(-1)*S*
     +    MS2*MG2 + 64*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(2)*N*CF**2 * (  - 32*TG**(-2)*S*MS2*MG2 - 
     +    32*UG**(-2)*S*MS2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(2)*CF**2 * ( 64*TG**(-1)*UG**(-1)*S*MS2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(3)*N*CF**2 * (  - 32*TG**(-2)*S*MS2*MG2 - 
     +    32*UG**(-2)*S*MS2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(3)*CF**2 * ( 64*TG**(-1)*UG**(-1)*S*MS2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N*CF * (  - 32*TG**(-1)*UG**(-1)*S*T1*
     +    MG2 - 32*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N*CF**2 * (  - 64*TG**(-2)*S*T1*MG2 - 
     +    64*TG**(-2)*S**2*MG2 - 32*UG**(-2)*S*T1*MG2 - 32*UG**(-2)*
     +    S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*N**2*CF * ( 16*TG**(-2)*S*T1*MG2 + 16*
     +    TG**(-2)*S**2*MG2 + 16*UG**(-2)*S*T1*MG2 + 16*UG**(-2)*S**2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(4)*CF**2 * ( 64*TG**(-1)*UG**(-1)*S*T1*MG2
     +     + 64*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*N*CF * ( 32*TG**(-1)*UG**(-1)*S*T1*MG2
     +     )
     +
      MQQLLV = MQQLLV + SOF1(5)*N*CF**2 * ( 32*TG**(-2)*S*T1*MG2 + 64*
     +    UG**(-2)*S*T1*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*N**2*CF * (  - 16*TG**(-2)*S*T1*MG2 - 
     +    16*UG**(-2)*S*T1*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(5)*CF**2 * (  - 64*TG**(-1)*UG**(-1)*S*T1*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*N*CF * ( 32*TG**(-1)*UG**(-1)*S*T1*MG2
     +     )
     +
      MQQLLV = MQQLLV + SOF1(6)*N*CF**2 * ( 32*TG**(-2)*S*T1*MG2 + 64*
     +    UG**(-2)*S*T1*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*N**2*CF * (  - 16*TG**(-2)*S*T1*MG2 - 
     +    16*UG**(-2)*S*T1*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(6)*CF**2 * (  - 64*TG**(-1)*UG**(-1)*S*T1*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N*CF * (  - 32*TG**(-1)*UG**(-1)*S*T1*
     +    MG2 - 32*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N*CF**2 * (  - 64*TG**(-2)*S*T1*MG2 - 
     +    64*TG**(-2)*S**2*MG2 - 32*UG**(-2)*S*T1*MG2 - 32*UG**(-2)*
     +    S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*N**2*CF * ( 16*TG**(-2)*S*T1*MG2 + 16*
     +    TG**(-2)*S**2*MG2 + 16*UG**(-2)*S*T1*MG2 + 16*UG**(-2)*S**2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(7)*CF**2 * ( 64*TG**(-1)*UG**(-1)*S*T1*MG2
     +     + 64*TG**(-1)*UG**(-1)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N*CF * (  - 64*TG**(-1)*UG**(-1)*S**2*
     +    MG2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MG2 - 
     +    32*TG**(-2)*S*MS2*MG2 + 32*TG**(-2)*S*MG2**2 - 64*TG**(-2)*
     +    S**2*MG2 + 32*TG**(-1)*S*MG2 + 32*UG**(-2)*S*T1*MG2 - 32*
     +    UG**(-2)*S*MS2*MG2 + 32*UG**(-2)*S*MG2**2 - 32*UG**(-2)*S**2*
     +    MG2 + 32*UG**(-1)*S*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*N**2*CF * ( 32*TG**(-2)*S**2*MG2 + 32*
     +    UG**(-2)*S**2*MG2 )
     +
      MQQLLV = MQQLLV + SOF1(8)*CF**2 * ( 64*TG**(-1)*UG**(-1)*S**2*MG2
     +     )

      MQQLRV = 0.D0
      MQQLRV = MQQLRV + N*CF * ( 8*TG**(-3)*S*T**(-1)*NS*MS2*MG2**2 - 8
     +    *TG**(-3)*S*T**(-1)*MG2**2*MT2 - 8*TG**(-3)*S*NS*MS2*MG2 + 8*
     +    TG**(-3)*S*MG2*MT2 + 8*TG**(-3)*T**(-1)*NS*MS2*MG2**3 - 16*
     +    TG**(-3)*T**(-1)*NS*MS2**2*MG2**2 + 8*TG**(-3)*T**(-1)*NS*
     +    MS2**3*MG2 + 16*TG**(-3)*T**(-1)*MS2*MG2**2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MS2**2*MG2*MT2 - 8*TG**(-3)*T**(-1)*MG2**3*MT2 - 8*
     +    TG**(-3)*NS*MS2*MG2**2 + 16*TG**(-3)*NS*MS2**2*MG2 - 8*
     +    TG**(-3)*NS*MS2**3 - 16*TG**(-3)*MS2*MG2*MT2 + 8*TG**(-3)*
     +    MS2**2*MT2 + 8*TG**(-3)*MG2**2*MT2 + 12*TG**(-2)*S*T**(-1)*NS
     +    *MS2*MG2 - 12*TG**(-2)*S*T**(-1)*MG2*MT2 - 8*TG**(-2)*S*NS*
     +    MS2 + 8*TG**(-2)*S*MT2 + 12*TG**(-2)*T**(-1)*T1*NS*MS2*MG2 - 
     +    4*TG**(-2)*T**(-1)*T1*NS*MS2**2 + 4*TG**(-2)*T**(-1)*T1*MS2*
     +    MT2 - 12*TG**(-2)*T**(-1)*T1*MG2*MT2 + 8*TG**(-2)*T**(-1)*NS*
     +    MS2*MG2**2 - 8*TG**(-2)*T**(-1)*NS*MS2**2*MG2 + 8*TG**(-2)*
     +    T**(-1)*MS2*MG2*MT2 - 8*TG**(-2)*T**(-1)*MG2**2*MT2 - 8*
     +    TG**(-2)*T1*NS*MS2 )
     +
      MQQLRV = MQQLRV + N*CF * ( 8*TG**(-2)*T1*MT2 - 8*TG**(-2)*NS*MS2*
     +    MG2 + 8*TG**(-2)*NS*MS2**2 - 8*TG**(-2)*MS2*MT2 + 8*TG**(-2)*
     +    MG2*MT2 + 4*TG**(-1)*S*T**(-1)*NS*MS2 - 4*TG**(-1)*S*T**(-1)*
     +    MT2 + 4*TG**(-1)*T**(-1)*T1*NS*MS2 - 4*TG**(-1)*T**(-1)*T1*
     +    MT2 + 8*UG**(-3)*S*U**(-1)*T1*NS*MS2*MG2 - 8*UG**(-3)*S*
     +    U**(-1)*T1*MG2*MT2 + 8*UG**(-3)*S*U**(-1)*NS*MS2**2*MG2 - 8*
     +    UG**(-3)*S*U**(-1)*MS2*MG2*MT2 - 8*UG**(-3)*S*T1*NS*MS2 + 8*
     +    UG**(-3)*S*T1*MT2 - 8*UG**(-3)*S*NS*MS2**2 + 8*UG**(-3)*S*MS2
     +    *MT2 + 8*UG**(-3)*U**(-1)*T1**2*NS*MS2*MG2 - 8*UG**(-3)*
     +    U**(-1)*T1**2*MG2*MT2 - 8*UG**(-3)*T1**2*NS*MS2 + 8*UG**(-3)*
     +    T1**2*MT2 + 4*UG**(-2)*S*U**(-1)*T1*NS*MS2 - 4*UG**(-2)*S*
     +    U**(-1)*T1*MT2 + 4*UG**(-2)*S*U**(-1)*NS*MS2**2 - 4*UG**(-2)*
     +    S*U**(-1)*MS2*MT2 + 4*UG**(-2)*U**(-1)*T1**2*NS*MS2 - 4*
     +    UG**(-2)*U**(-1)*T1**2*MT2 )
     +
      MQQLRV = MQQLRV + N**2*CF * (  - 8*TG**(-3)*S*T**(-1)*MG2**3 + 8*
     +    TG**(-3)*S*MG2**2 + 16*TG**(-3)*T**(-1)*MS2*MG2**3 - 8*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*MG2**4 - 
     +    16*TG**(-3)*MS2*MG2**2 + 8*TG**(-3)*MS2**2*MG2 + 8*TG**(-3)*
     +    MG2**3 - 12*TG**(-2)*S*T**(-1)*MG2**2 + 12*TG**(-2)*S*MG2 + 4
     +    *TG**(-2)*T**(-1)*T1*MS2*MG2 - 12*TG**(-2)*T**(-1)*T1*MG2**2
     +     + 8*TG**(-2)*T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*MG2**3
     +     - 4*TG**(-2)*T1*MS2 + 12*TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*
     +    MG2 + 8*TG**(-2)*MG2**2 - 4*TG**(-1)*S*T**(-1)*MG2 + 4*
     +    TG**(-1)*S - 4*TG**(-1)*T**(-1)*T1*MG2 + 4*TG**(-1)*T1 - 8*
     +    UG**(-3)*S*U**(-1)*T1*MG2**2 - 8*UG**(-3)*S*U**(-1)*MS2*
     +    MG2**2 + 8*UG**(-3)*S*T1*MG2 + 8*UG**(-3)*S*MS2*MG2 - 8*
     +    UG**(-3)*U**(-1)*T1**2*MG2**2 + 8*UG**(-3)*T1**2*MG2 - 4*
     +    UG**(-2)*S*U**(-1)*T1*MG2 - 4*UG**(-2)*S*U**(-1)*MS2*MG2 + 4*
     +    UG**(-2)*S*T1 + 4*UG**(-2)*S*MS2 - 4*UG**(-2)*U**(-1)*T1**2*
     +    MG2 )
     +
      MQQLRV = MQQLRV + N**2*CF * ( 4*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0A(1)*N**2*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 + 8*TG**(-3)*S*MG2**2 + 16*TG**(-3)*T**(-1)*MS2*MG2**3
     +     - 8*TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*
     +    MG2**4 - 16*TG**(-3)*MS2*MG2**2 + 8*TG**(-3)*MS2**2*MG2 + 8*
     +    TG**(-3)*MG2**3 - 12*TG**(-2)*S*T**(-1)*MG2**2 + 8*TG**(-2)*S
     +    *MG2 + 4*TG**(-2)*T**(-1)*T1*MS2*MG2 - 12*TG**(-2)*T**(-1)*T1
     +    *MG2**2 + 8*TG**(-2)*T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*
     +    MG2**3 + 8*TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*MG2 + 8*TG**(-2)*
     +    MG2**2 - 4*TG**(-1)*S*T**(-1)*MG2 - 4*TG**(-1)*T**(-1)*T1*MG2
     +     - 8*UG**(-3)*S*U**(-1)*T1*MG2**2 - 8*UG**(-3)*S*U**(-1)*MS2*
     +    MG2**2 + 8*UG**(-3)*S*T1*MG2 + 8*UG**(-3)*S*MS2*MG2 - 8*
     +    UG**(-3)*U**(-1)*T1**2*MG2**2 + 8*UG**(-3)*T1**2*MG2 - 4*
     +    UG**(-2)*S*U**(-1)*T1*MG2 - 4*UG**(-2)*S*U**(-1)*MS2*MG2 - 4*
     +    UG**(-2)*U**(-1)*T1**2*MG2 )
     +
      MQQLRV = MQQLRV + SK1B0A(2)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MG2**2*MT2 + 8*TG**(-3)*S*MG2*MT2 + 16*TG**(-3)*T**(-1)*MS2*
     +    MG2**2*MT2 - 8*TG**(-3)*T**(-1)*MS2**2*MG2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MG2**3*MT2 - 16*TG**(-3)*MS2*MG2*MT2 + 8*TG**(-3)*
     +    MS2**2*MT2 + 8*TG**(-3)*MG2**2*MT2 - 12*TG**(-2)*S*T**(-1)*
     +    MG2*MT2 + 8*TG**(-2)*S*MT2 + 4*TG**(-2)*T**(-1)*T1*MS2*MT2 - 
     +    12*TG**(-2)*T**(-1)*T1*MG2*MT2 + 8*TG**(-2)*T**(-1)*MS2*MG2*
     +    MT2 - 8*TG**(-2)*T**(-1)*MG2**2*MT2 + 8*TG**(-2)*T1*MT2 - 8*
     +    TG**(-2)*MS2*MT2 + 8*TG**(-2)*MG2*MT2 - 4*TG**(-1)*S*T**(-1)*
     +    MT2 - 4*TG**(-1)*T**(-1)*T1*MT2 - 8*UG**(-3)*S*U**(-1)*T1*MG2
     +    *MT2 - 8*UG**(-3)*S*U**(-1)*MS2*MG2*MT2 + 8*UG**(-3)*S*T1*MT2
     +     + 8*UG**(-3)*S*MS2*MT2 - 8*UG**(-3)*U**(-1)*T1**2*MG2*MT2 + 
     +    8*UG**(-3)*T1**2*MT2 - 4*UG**(-2)*S*U**(-1)*T1*MT2 - 4*
     +    UG**(-2)*S*U**(-1)*MS2*MT2 - 4*UG**(-2)*U**(-1)*T1**2*MT2 )
     +
      MQQLRV = MQQLRV + SK1B0A(3)*N*CF**2 * ( 4*TG**(-2)*S*MG2 - 4*
     +    TG**(-2)*T1*MS2 + 4*TG**(-2)*T1*MG2 + 4*TG**(-1)*S + 4*
     +    TG**(-1)*T1 + 4*UG**(-2)*S*T1 + 4*UG**(-2)*S*MS2 + 4*UG**(-2)
     +    *T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0C(1)*N*CF**2 * (  - 16*TG**(-2)*S*T1**(-1)
     +    *MS2**2 - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*S*MG2 - 16*TG**(-2)
     +    *T1*MG2 - 16*TG**(-1)*S - 16*TG**(-1)*T1 - 16*UG**(-2)*S*T1
     +     - 16*UG**(-2)*S*U1**(-1)*MS2**2 - 16*UG**(-2)*S*MS2 + 16*
     +    UG**(-2)*T1*MS2 - 16*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0C(2)*N*CF**2 * (  - 16*TG**(-2)*S*T1**(-1)
     +    *MS2*MG2 - 8*TG**(-2)*S*MG2 - 8*TG**(-2)*T1*MS2 - 8*TG**(-2)*
     +    T1*MG2 + 8*TG**(-1)*S + 8*TG**(-1)*T1 + 8*UG**(-2)*S*T1 - 16*
     +    UG**(-2)*S*U1**(-1)*MS2*MG2 + 8*UG**(-2)*S*MS2 + 16*UG**(-2)*
     +    T1*MG2 + 8*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*NS
     +    *MS2*MG2**2 + 8*TG**(-3)*S*T**(-1)*MS2*MG2**2 + 8*TG**(-3)*S*
     +    NS*MG2**2 - 8*TG**(-3)*S*MG2**2 - 8*TG**(-3)*T**(-1)*NS*MS2*
     +    MG2**3 + 16*TG**(-3)*T**(-1)*NS*MS2**2*MG2**2 - 8*TG**(-3)*
     +    T**(-1)*NS*MS2**3*MG2 + 8*TG**(-3)*T**(-1)*MS2*MG2**3 - 16*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 + 8*TG**(-3)*T**(-1)*MS2**3*
     +    MG2 - 16*TG**(-3)*NS*MS2*MG2**2 + 8*TG**(-3)*NS*MS2**2*MG2 + 
     +    8*TG**(-3)*NS*MG2**3 + 16*TG**(-3)*MS2*MG2**2 - 8*TG**(-3)*
     +    MS2**2*MG2 - 8*TG**(-3)*MG2**3 - 12*TG**(-2)*S*T**(-1)*NS*MS2
     +    *MG2 + 12*TG**(-2)*S*T**(-1)*MS2*MG2 + 12*TG**(-2)*S*NS*MG2
     +     - 12*TG**(-2)*S*MG2 - 12*TG**(-2)*T**(-1)*T1*NS*MS2*MG2 + 4*
     +    TG**(-2)*T**(-1)*T1*NS*MS2**2 + 12*TG**(-2)*T**(-1)*T1*MS2*
     +    MG2 - 4*TG**(-2)*T**(-1)*T1*MS2**2 - 8*TG**(-2)*T**(-1)*NS*
     +    MS2*MG2**2 + 8*TG**(-2)*T**(-1)*NS*MS2**2*MG2 + 8*TG**(-2)*
     +    T**(-1)*MS2*MG2**2 - 8*TG**(-2)*T**(-1)*MS2**2*MG2 - 4*
     +    TG**(-2)*T1*NS*MS2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N*CF * ( 12*TG**(-2)*T1*NS*MG2 + 4*
     +    TG**(-2)*T1*MS2 - 12*TG**(-2)*T1*MG2 - 8*TG**(-2)*NS*MS2*MG2
     +     + 8*TG**(-2)*NS*MG2**2 + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*
     +    MG2**2 - 4*TG**(-1)*S*T**(-1)*NS*MS2 + 4*TG**(-1)*S*T**(-1)*
     +    MS2 + 4*TG**(-1)*S*NS - 4*TG**(-1)*S - 4*TG**(-1)*T**(-1)*T1*
     +    NS*MS2 + 4*TG**(-1)*T**(-1)*T1*MS2 + 4*TG**(-1)*T1*NS - 4*
     +    TG**(-1)*T1 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2 + 16*TG**(-2)*S*T1**(-1)*MS2**2 + 16*TG**(-2)*S*MS2
     +     + 32*TG**(-2)*S*MG2 + 32*TG**(-2)*T1*MG2 + 16*TG**(-1)*S + 
     +    16*TG**(-1)*T1 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2 - 8*TG**(-2)*S*T1**(-1)*MS2**2 - 8*TG**(-2)*
     +    S*MS2 - 16*TG**(-2)*S*MG2 - 16*TG**(-2)*T1*MG2 - 8*TG**(-1)*S
     +     - 8*TG**(-1)*T1 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N*CF * (  - 8*UG**(-3)*S*U**(-1)*T1
     +    *NS*MS2*MG2 + 8*UG**(-3)*S*U**(-1)*T1*MS2*MG2 - 8*UG**(-3)*S*
     +    U**(-1)*NS*MS2**2*MG2 + 8*UG**(-3)*S*U**(-1)*MS2**2*MG2 + 8*
     +    UG**(-3)*S*T1*NS*MG2 - 8*UG**(-3)*S*T1*MG2 + 8*UG**(-3)*S*NS*
     +    MS2*MG2 - 8*UG**(-3)*S*MS2*MG2 - 8*UG**(-3)*U**(-1)*T1**2*NS*
     +    MS2*MG2 + 8*UG**(-3)*U**(-1)*T1**2*MS2*MG2 + 8*UG**(-3)*T1**2
     +    *NS*MG2 - 8*UG**(-3)*T1**2*MG2 - 4*UG**(-2)*S*U**(-1)*T1*NS*
     +    MS2 + 4*UG**(-2)*S*U**(-1)*T1*MS2 - 4*UG**(-2)*S*U**(-1)*NS*
     +    MS2**2 + 4*UG**(-2)*S*U**(-1)*MS2**2 + 4*UG**(-2)*S*T1*NS - 4
     +    *UG**(-2)*S*T1 + 4*UG**(-2)*S*NS*MS2 - 4*UG**(-2)*S*MS2 - 4*
     +    UG**(-2)*U**(-1)*T1**2*NS*MS2 + 4*UG**(-2)*U**(-1)*T1**2*MS2
     +     + 4*UG**(-2)*T1**2*NS - 4*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N*CF**2 * ( 16*UG**(-2)*S*T1 + 16*
     +    UG**(-2)*S*U1**(-1)*MS2*MG2 + 16*UG**(-2)*S*U1**(-1)*MS2**2
     +     + 16*UG**(-2)*S*MS2 - 16*UG**(-2)*T1*MS2 - 16*UG**(-2)*T1*
     +    MG2 + 16*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0D(1,2)*N**2*CF * (  - 8*UG**(-2)*S*T1 - 8*
     +    UG**(-2)*S*U1**(-1)*MS2*MG2 - 8*UG**(-2)*S*U1**(-1)*MS2**2 - 
     +    8*UG**(-2)*S*MS2 + 8*UG**(-2)*T1*MS2 + 8*UG**(-2)*T1*MG2 - 8*
     +    UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0D(2,1)*N**2*CF * ( 8*TG**(-3)*S*T**(-1)*
     +    MG2**3 - 24*TG**(-3)*S*MG2**2 - 16*TG**(-3)*T**(-1)*MS2*
     +    MG2**3 + 8*TG**(-3)*T**(-1)*MS2**2*MG2**2 + 8*TG**(-3)*
     +    T**(-1)*MG2**4 + 48*TG**(-3)*MS2*MG2**2 - 24*TG**(-3)*MS2**2*
     +    MG2 - 24*TG**(-3)*MG2**3 + 12*TG**(-2)*S*T**(-1)*MG2**2 + 8*
     +    TG**(-2)*S*T1**(-1)*MS2*MG2 + 8*TG**(-2)*S*T1**(-1)*MS2**2 + 
     +    8*TG**(-2)*S*MS2 - 20*TG**(-2)*S*MG2 - 4*TG**(-2)*T**(-1)*T1*
     +    MS2*MG2 + 12*TG**(-2)*T**(-1)*T1*MG2**2 - 8*TG**(-2)*T**(-1)*
     +    MS2*MG2**2 + 8*TG**(-2)*T**(-1)*MG2**3 + 12*TG**(-2)*T1*MS2
     +     - 20*TG**(-2)*T1*MG2 + 24*TG**(-2)*MS2*MG2 - 24*TG**(-2)*
     +    MG2**2 + 4*TG**(-1)*S*T**(-1)*MG2 - 4*TG**(-1)*S + 4*TG**(-1)
     +    *T**(-1)*T1*MG2 - 4*TG**(-1)*T1 )
     +
      MQQLRV = MQQLRV + SK1B0D(2,2)*N**2*CF * ( 8*UG**(-3)*S*U**(-1)*T1
     +    *MG2**2 + 8*UG**(-3)*S*U**(-1)*MS2*MG2**2 - 24*UG**(-3)*S*T1*
     +    MG2 - 24*UG**(-3)*S*MS2*MG2 + 8*UG**(-3)*U**(-1)*T1**2*MG2**2
     +     - 24*UG**(-3)*T1**2*MG2 + 4*UG**(-2)*S*U**(-1)*T1*MG2 + 4*
     +    UG**(-2)*S*U**(-1)*MS2*MG2 - 4*UG**(-2)*S*T1 + 8*UG**(-2)*S*
     +    U1**(-1)*MS2*MG2 + 8*UG**(-2)*S*U1**(-1)*MS2**2 - 4*UG**(-2)*
     +    S*MS2 + 4*UG**(-2)*U**(-1)*T1**2*MG2 - 8*UG**(-2)*T1*MS2 - 8*
     +    UG**(-2)*T1*MG2 - 4*UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0D(3,1)*N*CF * (  - 8*TG**(-3)*S*T**(-1)*
     +    MS2*MG2**2 + 8*TG**(-3)*S*T**(-1)*MG2**2*MT2 + 8*TG**(-3)*S*
     +    MG2**2 - 16*TG**(-3)*T**(-1)*MS2*MG2**2*MT2 - 8*TG**(-3)*
     +    T**(-1)*MS2*MG2**3 + 8*TG**(-3)*T**(-1)*MS2**2*MG2*MT2 + 16*
     +    TG**(-3)*T**(-1)*MS2**2*MG2**2 - 8*TG**(-3)*T**(-1)*MS2**3*
     +    MG2 + 8*TG**(-3)*T**(-1)*MG2**3*MT2 - 16*TG**(-3)*MS2*MG2**2
     +     + 8*TG**(-3)*MS2**2*MG2 + 8*TG**(-3)*MG2**3 - 12*TG**(-2)*S*
     +    T**(-1)*MS2*MG2 + 12*TG**(-2)*S*T**(-1)*MG2*MT2 + 12*TG**(-2)
     +    *S*MG2 - 12*TG**(-2)*T**(-1)*T1*MS2*MG2 - 4*TG**(-2)*T**(-1)*
     +    T1*MS2*MT2 + 4*TG**(-2)*T**(-1)*T1*MS2**2 + 12*TG**(-2)*
     +    T**(-1)*T1*MG2*MT2 - 8*TG**(-2)*T**(-1)*MS2*MG2*MT2 - 8*
     +    TG**(-2)*T**(-1)*MS2*MG2**2 + 8*TG**(-2)*T**(-1)*MS2**2*MG2
     +     + 8*TG**(-2)*T**(-1)*MG2**2*MT2 - 4*TG**(-2)*T1*MS2 + 12*
     +    TG**(-2)*T1*MG2 - 8*TG**(-2)*MS2*MG2 + 8*TG**(-2)*MG2**2 - 4*
     +    TG**(-1)*S*T**(-1)*MS2 + 4*TG**(-1)*S*T**(-1)*MT2 + 4*
     +    TG**(-1)*S )
     +
      MQQLRV = MQQLRV + SK1B0D(3,1)*N*CF * (  - 4*TG**(-1)*T**(-1)*T1*
     +    MS2 + 4*TG**(-1)*T**(-1)*T1*MT2 + 4*TG**(-1)*T1 )
     +
      MQQLRV = MQQLRV + SK1B0D(3,2)*N*CF * (  - 8*UG**(-3)*S*U**(-1)*T1
     +    *MS2*MG2 + 8*UG**(-3)*S*U**(-1)*T1*MG2*MT2 + 8*UG**(-3)*S*
     +    U**(-1)*MS2*MG2*MT2 - 8*UG**(-3)*S*U**(-1)*MS2**2*MG2 + 8*
     +    UG**(-3)*S*T1*MG2 + 8*UG**(-3)*S*MS2*MG2 - 8*UG**(-3)*U**(-1)
     +    *T1**2*MS2*MG2 + 8*UG**(-3)*U**(-1)*T1**2*MG2*MT2 + 8*
     +    UG**(-3)*T1**2*MG2 - 4*UG**(-2)*S*U**(-1)*T1*MS2 + 4*UG**(-2)
     +    *S*U**(-1)*T1*MT2 + 4*UG**(-2)*S*U**(-1)*MS2*MT2 - 4*UG**(-2)
     +    *S*U**(-1)*MS2**2 + 4*UG**(-2)*S*T1 + 4*UG**(-2)*S*MS2 - 4*
     +    UG**(-2)*U**(-1)*T1**2*MS2 + 4*UG**(-2)*U**(-1)*T1**2*MT2 + 4
     +    *UG**(-2)*T1**2 )
     +
      MQQLRV = MQQLRV + SK1B0E(1)*N**2*CF * ( 16*TG**(-3)*S*MG2**2 - 32
     +    *TG**(-3)*MS2*MG2**2 + 16*TG**(-3)*MS2**2*MG2 + 16*TG**(-3)*
     +    MG2**3 + 16*TG**(-2)*S*MG2 + 16*TG**(-2)*T1*MG2 - 16*TG**(-2)
     +    *MS2*MG2 + 16*TG**(-2)*MG2**2 + 16*UG**(-3)*S*T1*MG2 + 16*
     +    UG**(-3)*S*MS2*MG2 + 16*UG**(-3)*T1**2*MG2 )
     +
      MQQLRV = MQQLRV + SK1B0E(2)*N*CF * ( 8*TG**(-3)*S*NS*MS2*MG2 - 8*
     +    TG**(-3)*S*NS*MG2**2 - 8*TG**(-3)*S*MS2*MG2 + 8*TG**(-3)*S*
     +    MG2**2 + 24*TG**(-3)*NS*MS2*MG2**2 - 24*TG**(-3)*NS*MS2**2*
     +    MG2 + 8*TG**(-3)*NS*MS2**3 - 8*TG**(-3)*NS*MG2**3 - 24*
     +    TG**(-3)*MS2*MG2**2 + 24*TG**(-3)*MS2**2*MG2 - 8*TG**(-3)*
     +    MS2**3 + 8*TG**(-3)*MG2**3 + 8*TG**(-2)*S*NS*MS2 - 8*TG**(-2)
     +    *S*NS*MG2 - 8*TG**(-2)*S*MS2 + 8*TG**(-2)*S*MG2 + 8*TG**(-2)*
     +    T1*NS*MS2 - 8*TG**(-2)*T1*NS*MG2 - 8*TG**(-2)*T1*MS2 + 8*
     +    TG**(-2)*T1*MG2 + 16*TG**(-2)*NS*MS2*MG2 - 8*TG**(-2)*NS*
     +    MS2**2 - 8*TG**(-2)*NS*MG2**2 - 16*TG**(-2)*MS2*MG2 + 8*
     +    TG**(-2)*MS2**2 + 8*TG**(-2)*MG2**2 + 8*UG**(-3)*S*T1*NS*MS2
     +     - 8*UG**(-3)*S*T1*NS*MG2 - 8*UG**(-3)*S*T1*MS2 + 8*UG**(-3)*
     +    S*T1*MG2 - 8*UG**(-3)*S*NS*MS2*MG2 + 8*UG**(-3)*S*NS*MS2**2
     +     + 8*UG**(-3)*S*MS2*MG2 - 8*UG**(-3)*S*MS2**2 + 8*UG**(-3)*
     +    T1**2*NS*MS2 - 8*UG**(-3)*T1**2*NS*MG2 - 8*UG**(-3)*T1**2*MS2
     +     + 8*UG**(-3)*T1**2*MG2 )
     +
      MQQLRV = MQQLRV + SK1B0E(3)*N*CF * ( 8*TG**(-3)*S*MS2*MG2 - 8*
     +    TG**(-3)*S*MG2*MT2 - 8*TG**(-3)*S*MG2**2 + 16*TG**(-3)*MS2*
     +    MG2*MT2 + 24*TG**(-3)*MS2*MG2**2 - 24*TG**(-3)*MS2**2*MG2 - 8
     +    *TG**(-3)*MS2**2*MT2 + 8*TG**(-3)*MS2**3 - 8*TG**(-3)*MG2**2*
     +    MT2 - 8*TG**(-3)*MG2**3 + 8*TG**(-2)*S*MS2 - 8*TG**(-2)*S*MG2
     +     - 8*TG**(-2)*S*MT2 + 8*TG**(-2)*T1*MS2 - 8*TG**(-2)*T1*MG2
     +     - 8*TG**(-2)*T1*MT2 + 16*TG**(-2)*MS2*MG2 + 8*TG**(-2)*MS2*
     +    MT2 - 8*TG**(-2)*MS2**2 - 8*TG**(-2)*MG2*MT2 - 8*TG**(-2)*
     +    MG2**2 + 8*UG**(-3)*S*T1*MS2 - 8*UG**(-3)*S*T1*MG2 - 8*
     +    UG**(-3)*S*T1*MT2 - 8*UG**(-3)*S*MS2*MG2 - 8*UG**(-3)*S*MS2*
     +    MT2 + 8*UG**(-3)*S*MS2**2 + 8*UG**(-3)*T1**2*MS2 - 8*UG**(-3)
     +    *T1**2*MG2 - 8*UG**(-3)*T1**2*MT2 )
     +
      MQQLRV = MQQLRV + SK1BP(1)*N*CF**2 * (  - 16*TG**(-2)*S*MS2*MG2
     +     - 16*TG**(-2)*T1*MS2*MG2 + 16*TG**(-2)*T1*MS2**2 - 16*
     +    TG**(-1)*S*MS2 - 16*TG**(-1)*T1*MS2 - 16*UG**(-2)*S*T1*MS2 - 
     +    16*UG**(-2)*S*MS2**2 - 16*UG**(-2)*T1**2*MS2 )
     +
      MQQLRV = MQQLRV + SK1BP(2)*N*CF**2 * ( 8*TG**(-2)*S*MS2*MG2 - 8*
     +    TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*MS2*MG2 - 8*TG**(-2)*T1*
     +    MS2**2 - 8*TG**(-2)*T1*MG2**2 + 8*TG**(-1)*S*MS2 - 8*TG**(-1)
     +    *S*MG2 + 8*TG**(-1)*T1*MS2 - 8*TG**(-1)*T1*MG2 + 8*UG**(-2)*S
     +    *T1*MS2 - 8*UG**(-2)*S*T1*MG2 - 8*UG**(-2)*S*MS2*MG2 + 8*
     +    UG**(-2)*S*MS2**2 + 8*UG**(-2)*T1**2*MS2 - 8*UG**(-2)*T1**2*
     +    MG2 )
     +
      MQQLRV = MQQLRV + SK1BP(3)*N*CF**2 * (  - 4*TG**(-2)*S*MS2*MG2 + 
     +    4*TG**(-2)*S*MG2**2 - 8*TG**(-2)*T1*MS2*MG2 + 4*TG**(-2)*T1*
     +    MS2**2 + 4*TG**(-2)*T1*MG2**2 - 4*TG**(-1)*S*MS2 + 4*TG**(-1)
     +    *S*MG2 - 4*TG**(-1)*T1*MS2 + 4*TG**(-1)*T1*MG2 - 4*UG**(-2)*S
     +    *T1*MS2 + 4*UG**(-2)*S*T1*MG2 + 4*UG**(-2)*S*MS2*MG2 - 4*
     +    UG**(-2)*S*MS2**2 - 4*UG**(-2)*T1**2*MS2 + 4*UG**(-2)*T1**2*
     +    MG2 )
     +
      MQQLRV = MQQLRV + SK1C0A(1)*N*CF**2 * (  - 8*TG**(-1)*S*T1 + 8*
     +    TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 - 8*TG**(-1)*S**2 + 8*
     +    UG**(-1)*S*T1 + 8*UG**(-1)*S*MS2 - 8*UG**(-1)*S*MG2 )
     +
      MQQLRV = MQQLRV + SK1C0A(1)*N**2*CF * ( 4*TG**(-1)*S*T1 - 4*
     +    TG**(-1)*S*MS2 + 4*TG**(-1)*S*MG2 + 4*TG**(-1)*S**2 - 4*
     +    UG**(-1)*S*T1 - 4*UG**(-1)*S*MS2 + 4*UG**(-1)*S*MG2 )
     +
      MQQLRV = MQQLRV + SK1C0A(5)*N*CF**2 * ( 8*TG**(-1)*S*T1 - 8*
     +    TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2 - 8*
     +    UG**(-2)*S*T1*MS2 + 8*UG**(-2)*S*T1*MG2 - 16*UG**(-2)*S*MS2*
     +    MG2 + 8*UG**(-2)*S*MS2**2 + 8*UG**(-2)*S*MG2**2 - 8*UG**(-2)*
     +    S**2*MS2 + 8*UG**(-2)*S**2*MG2 - 24*UG**(-1)*S*MS2 + 8*
     +    UG**(-1)*S*MG2 + 8*UG**(-1)*S**2 - 8*S )
     +
      MQQLRV = MQQLRV + SK1C0A(5)*N**2*CF * (  - 4*TG**(-1)*S*T1 + 4*
     +    TG**(-1)*S*MS2 + 4*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2 + 4*
     +    UG**(-2)*S*T1*MS2 - 4*UG**(-2)*S*T1*MG2 + 8*UG**(-2)*S*MS2*
     +    MG2 - 4*UG**(-2)*S*MS2**2 - 4*UG**(-2)*S*MG2**2 + 4*UG**(-2)*
     +    S**2*MS2 - 4*UG**(-2)*S**2*MG2 + 12*UG**(-1)*S*MS2 - 4*
     +    UG**(-1)*S*MG2 - 4*UG**(-1)*S**2 + 4*S )
     +
      MQQLRV = MQQLRV + SK1C0B(1)*N*CF**2 * (  - 8*TG**(-1)*S*T1 + 24*
     +    TG**(-1)*S*MS2 - 8*TG**(-1)*S*MG2 - 8*TG**(-1)*S**2 + 32*
     +    TG**(-1)*T1*MS2 + 24*UG**(-2)*S*T1*MS2 + 8*UG**(-2)*S*T1*MG2
     +     - 16*UG**(-2)*S*MS2*MG2 + 8*UG**(-2)*S*MS2**2 + 8*UG**(-2)*S
     +    *MG2**2 - 8*UG**(-2)*S**2*MS2 + 8*UG**(-2)*S**2*MG2 + 32*
     +    UG**(-2)*T1*MS2*MG2 - 32*UG**(-2)*T1*MS2**2 + 32*UG**(-2)*
     +    T1**2*MS2 + 16*UG**(-1)*S*T1 - 24*UG**(-1)*S*MS2 + 8*UG**(-1)
     +    *S*MG2 + 8*UG**(-1)*S**2 - 8*S )
     +
      MQQLRV = MQQLRV + SK1C0B(1)*N**2*CF * ( 4*TG**(-1)*S*T1 - 12*
     +    TG**(-1)*S*MS2 + 4*TG**(-1)*S*MG2 + 4*TG**(-1)*S**2 - 16*
     +    TG**(-1)*T1*MS2 - 12*UG**(-2)*S*T1*MS2 - 4*UG**(-2)*S*T1*MG2
     +     + 8*UG**(-2)*S*MS2*MG2 - 4*UG**(-2)*S*MS2**2 - 4*UG**(-2)*S*
     +    MG2**2 + 4*UG**(-2)*S**2*MS2 - 4*UG**(-2)*S**2*MG2 - 16*
     +    UG**(-2)*T1*MS2*MG2 + 16*UG**(-2)*T1*MS2**2 - 16*UG**(-2)*
     +    T1**2*MS2 - 8*UG**(-1)*S*T1 + 12*UG**(-1)*S*MS2 - 4*UG**(-1)*
     +    S*MG2 - 4*UG**(-1)*S**2 + 4*S )
     +
      MQQLRV = MQQLRV + SK1C0B(4)*N*CF**2 * (  - 8*TG**(-1)*S*T1*
     +    S1**(-1)*MS2 + 8*TG**(-1)*S*T1*S1**(-1)*MG2 + 8*TG**(-1)*S*T1
     +     - 16*TG**(-1)*S*S1**(-1)*MS2*MG2 + 24*TG**(-1)*S*S1**(-1)*
     +    MS2**2 - 8*TG**(-1)*S*S1**(-1)*MG2**2 - 32*TG**(-1)*S*MS2 + 
     +    16*TG**(-1)*S*MG2 - 8*TG**(-1)*S**2*S1**(-1)*MS2 + 8*TG**(-1)
     +    *S**2*S1**(-1)*MG2 + 8*TG**(-1)*S**2 - 16*TG**(-1)*T1*MS2 + 
     +    16*TG**(-1)*T1*MG2 + 16*UG**(-2)*S*T1*S1**(-1)*MS2*MG2 - 24*
     +    UG**(-2)*S*T1*S1**(-1)*MS2**2 + 8*UG**(-2)*S*T1*S1**(-1)*
     +    MG2**2 - 8*UG**(-2)*S*T1*MS2 + 8*UG**(-2)*S*T1*MG2 - 8*
     +    UG**(-2)*S*T1**2*S1**(-1)*MS2 + 8*UG**(-2)*S*T1**2*S1**(-1)*
     +    MG2 + 32*UG**(-2)*S*S1**(-1)*MS2*MG2**2 - 64*UG**(-2)*S*
     +    S1**(-1)*MS2**2*MG2 + 32*UG**(-2)*S*S1**(-1)*MS2**3 + 16*
     +    UG**(-2)*S*MS2*MG2 - 8*UG**(-2)*S*MS2**2 - 8*UG**(-2)*S*
     +    MG2**2 - 8*UG**(-2)*S**2*T1*S1**(-1)*MS2 + 8*UG**(-2)*S**2*T1
     +    *S1**(-1)*MG2 + 32*UG**(-2)*S**2*S1**(-1)*MS2*MG2 - 32*
     +    UG**(-2)*S**2*S1**(-1)*MS2**2 )
     +
      MQQLRV = MQQLRV + SK1C0B(4)*N*CF**2 * ( 8*UG**(-2)*S**2*MS2 - 8*
     +    UG**(-2)*S**2*MG2 - 32*UG**(-2)*T1*MS2*MG2 + 16*UG**(-2)*T1*
     +    MS2**2 + 16*UG**(-2)*T1*MG2**2 - 16*UG**(-2)*T1**2*MS2 + 16*
     +    UG**(-2)*T1**2*MG2 + 16*UG**(-1)*S*T1*S1**(-1)*MS2 + 16*
     +    UG**(-1)*S*T1*S1**(-1)*MG2 - 16*UG**(-1)*S*T1 + 8*UG**(-1)*S*
     +    T1**2*S1**(-1) + 32*UG**(-1)*S*S1**(-1)*MS2*MG2 - 32*UG**(-1)
     +    *S*S1**(-1)*MS2**2 - 16*UG**(-1)*S*MG2 + 8*UG**(-1)*S**2*T1*
     +    S1**(-1) + 24*UG**(-1)*S**2*S1**(-1)*MS2 + 8*UG**(-1)*S**2*
     +    S1**(-1)*MG2 - 8*UG**(-1)*S**2 + 8*S*T1*S1**(-1) - 24*S*
     +    S1**(-1)*MS2 - 8*S*S1**(-1)*MG2 + 8*S + 16*S**2*S1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0B(4)*N**2*CF * ( 4*TG**(-1)*S*T1*S1**(-1)*
     +    MS2 - 4*TG**(-1)*S*T1*S1**(-1)*MG2 - 4*TG**(-1)*S*T1 + 8*
     +    TG**(-1)*S*S1**(-1)*MS2*MG2 - 12*TG**(-1)*S*S1**(-1)*MS2**2
     +     + 4*TG**(-1)*S*S1**(-1)*MG2**2 + 16*TG**(-1)*S*MS2 - 8*
     +    TG**(-1)*S*MG2 + 4*TG**(-1)*S**2*S1**(-1)*MS2 - 4*TG**(-1)*
     +    S**2*S1**(-1)*MG2 - 4*TG**(-1)*S**2 + 8*TG**(-1)*T1*MS2 - 8*
     +    TG**(-1)*T1*MG2 - 8*UG**(-2)*S*T1*S1**(-1)*MS2*MG2 + 12*
     +    UG**(-2)*S*T1*S1**(-1)*MS2**2 - 4*UG**(-2)*S*T1*S1**(-1)*
     +    MG2**2 + 4*UG**(-2)*S*T1*MS2 - 4*UG**(-2)*S*T1*MG2 + 4*
     +    UG**(-2)*S*T1**2*S1**(-1)*MS2 - 4*UG**(-2)*S*T1**2*S1**(-1)*
     +    MG2 - 16*UG**(-2)*S*S1**(-1)*MS2*MG2**2 + 32*UG**(-2)*S*
     +    S1**(-1)*MS2**2*MG2 - 16*UG**(-2)*S*S1**(-1)*MS2**3 - 8*
     +    UG**(-2)*S*MS2*MG2 + 4*UG**(-2)*S*MS2**2 + 4*UG**(-2)*S*
     +    MG2**2 + 4*UG**(-2)*S**2*T1*S1**(-1)*MS2 - 4*UG**(-2)*S**2*T1
     +    *S1**(-1)*MG2 - 16*UG**(-2)*S**2*S1**(-1)*MS2*MG2 + 16*
     +    UG**(-2)*S**2*S1**(-1)*MS2**2 )
     +
      MQQLRV = MQQLRV + SK1C0B(4)*N**2*CF * (  - 4*UG**(-2)*S**2*MS2 + 
     +    4*UG**(-2)*S**2*MG2 + 16*UG**(-2)*T1*MS2*MG2 - 8*UG**(-2)*T1*
     +    MS2**2 - 8*UG**(-2)*T1*MG2**2 + 8*UG**(-2)*T1**2*MS2 - 8*
     +    UG**(-2)*T1**2*MG2 - 8*UG**(-1)*S*T1*S1**(-1)*MS2 - 8*
     +    UG**(-1)*S*T1*S1**(-1)*MG2 + 8*UG**(-1)*S*T1 - 4*UG**(-1)*S*
     +    T1**2*S1**(-1) - 16*UG**(-1)*S*S1**(-1)*MS2*MG2 + 16*UG**(-1)
     +    *S*S1**(-1)*MS2**2 + 8*UG**(-1)*S*MG2 - 4*UG**(-1)*S**2*T1*
     +    S1**(-1) - 12*UG**(-1)*S**2*S1**(-1)*MS2 - 4*UG**(-1)*S**2*
     +    S1**(-1)*MG2 + 4*UG**(-1)*S**2 - 4*S*T1*S1**(-1) + 12*S*
     +    S1**(-1)*MS2 + 4*S*S1**(-1)*MG2 - 4*S - 8*S**2*S1**(-1) )
     +
      MQQLRV = MQQLRV + SK1C0C(3,1)*N*CF**2 * (  - 16*TG**(-2)*S*MS2*
     +    MG2 + 16*TG**(-2)*S*MG2**2 - 32*TG**(-2)*T1*MS2*MG2 + 16*
     +    TG**(-2)*T1*MS2**2 + 16*TG**(-2)*T1*MG2**2 - 16*TG**(-1)*S*
     +    MS2 + 32*TG**(-1)*S*MG2 - 32*TG**(-1)*T1*MS2 + 32*TG**(-1)*T1
     +    *MG2 + 16*UG**(-2)*S*T1*MS2 + 16*UG**(-2)*S*T1**2 + 16*
     +    UG**(-2)*T1*MS2*MG2 - 16*UG**(-2)*T1*MS2**2 + 16*UG**(-2)*
     +    T1**2*MG2 + 16*UG**(-2)*T1**3 + 16*S + 16*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,1)*N**2*CF * ( 8*TG**(-2)*S*MS2*MG2 - 
     +    8*TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*MS2*MG2 - 8*TG**(-2)*T1*
     +    MS2**2 - 8*TG**(-2)*T1*MG2**2 + 8*TG**(-1)*S*MS2 - 16*
     +    TG**(-1)*S*MG2 + 16*TG**(-1)*T1*MS2 - 16*TG**(-1)*T1*MG2 - 4*
     +    UG**(-2)*S*T1*MS2 - 4*UG**(-2)*S*T1**2 - 4*UG**(-2)*T1*MS2*
     +    MG2 + 4*UG**(-2)*T1*MS2**2 - 4*UG**(-2)*T1**2*MG2 - 4*
     +    UG**(-2)*T1**3 - 8*S - 8*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,2)*N*CF**2 * (  - 16*TG**(-1)*S*T1 + 
     +    32*TG**(-1)*S*MS2 - 16*TG**(-1)*S*MG2 - 16*TG**(-1)*S**2 + 32
     +    *TG**(-1)*T1*MS2 - 16*TG**(-1)*T1*MG2 - 16*UG**(-2)*S*T1*MS2
     +     + 16*UG**(-2)*S*T1*MG2 - 16*UG**(-2)*S*T1**2 + 16*UG**(-2)*S
     +    *MS2*MG2 - 16*UG**(-2)*S*MS2**2 - 16*UG**(-2)*T1**3 + 16*
     +    UG**(-1)*S*T1 + 16*UG**(-1)*S*MS2 - 16*S - 16*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(3,2)*N**2*CF * ( 4*TG**(-1)*S*T1 - 8*
     +    TG**(-1)*S*MS2 + 4*TG**(-1)*S*MG2 + 4*TG**(-1)*S**2 - 8*
     +    TG**(-1)*T1*MS2 + 4*TG**(-1)*T1*MG2 + 8*UG**(-2)*S*T1*MS2 - 8
     +    *UG**(-2)*S*T1*MG2 + 8*UG**(-2)*S*T1**2 - 8*UG**(-2)*S*MS2*
     +    MG2 + 8*UG**(-2)*S*MS2**2 + 8*UG**(-2)*T1**3 - 8*UG**(-1)*S*
     +    T1 - 8*UG**(-1)*S*MS2 + 4*S + 4*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(4,1)*N*CF**2 * (  - 16*TG**(-1)*T1*MG2 )
     +
      MQQLRV = MQQLRV + SK1C0C(4,1)*N**2*CF * ( 8*TG**(-2)*S*MS2*MG2 - 
     +    8*TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*MS2*MG2 - 8*TG**(-2)*T1*
     +    MS2**2 - 8*TG**(-2)*T1*MG2**2 + 4*TG**(-1)*S*MS2 - 12*
     +    TG**(-1)*S*MG2 + 8*TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 - 4*S
     +     - 4*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(4,2)*N*CF**2 * (  - 16*UG**(-2)*S*T1*MG2
     +     + 16*UG**(-2)*T1*MS2*MG2 - 16*UG**(-2)*T1*MG2**2 - 16*
     +    UG**(-2)*T1**2*MG2 + 16*UG**(-1)*S*MG2 )
     +
      MQQLRV = MQQLRV + SK1C0C(4,2)*N**2*CF * ( 4*UG**(-2)*S*T1*MS2 + 4
     +    *UG**(-2)*S*T1**2 - 8*UG**(-2)*S*MS2*MG2 + 8*UG**(-2)*S*
     +    MS2**2 - 12*UG**(-2)*T1*MS2*MG2 + 4*UG**(-2)*T1*MS2**2 + 8*
     +    UG**(-2)*T1*MG2**2 + 4*UG**(-2)*T1**2*MG2 + 4*UG**(-2)*T1**3
     +     - 4*UG**(-1)*S*T1 - 4*UG**(-1)*S*MS2 - 8*UG**(-1)*S*MG2 )
     +
      MQQLRV = MQQLRV + SK1C0C(5,1)*N*CF**2 * (  - 16*TG**(-1)*T1*MG2
     +     - 16*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(5,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2**2*MG2 + 8*TG**(-2)*S*T1**(-1)*MS2**3 - 8*
     +    TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*S*MS2**2 - 8*TG**(-2)*S*
     +    MG2**2 + 8*TG**(-2)*T1*MS2**2 - 8*TG**(-2)*T1*MG2**2 + 4*
     +    TG**(-1)*S*MS2 - 12*TG**(-1)*S*MG2 + 8*TG**(-1)*T1*MS2 - 4*
     +    TG**(-1)*T1*MG2 - 4*S + 4*T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(5,2)*N*CF**2 * (  - 16*UG**(-2)*S*T1*MG2
     +     + 16*UG**(-2)*S*T1**2 + 32*UG**(-2)*S*MS2*MG2 - 16*UG**(-2)*
     +    S*MS2**2 - 16*UG**(-2)*S*MG2**2 + 16*UG**(-2)*S**2*MS2 - 16*
     +    UG**(-2)*S**2*MG2 - 16*UG**(-2)*T1*MS2*MG2 + 16*UG**(-2)*T1*
     +    MS2**2 - 32*UG**(-2)*T1**2*MS2 + 16*UG**(-2)*T1**2*MG2 + 16*
     +    UG**(-2)*T1**3 - 32*UG**(-1)*S*T1 + 32*UG**(-1)*S*MS2 - 16*
     +    UG**(-1)*S*MG2 - 16*UG**(-1)*S**2 )
     +
      MQQLRV = MQQLRV + SK1C0C(5,2)*N**2*CF * ( 4*UG**(-2)*S*T1*MS2 - 4
     +    *UG**(-2)*S*T1**2 - 8*UG**(-2)*S*U1**(-1)*MS2**2*MG2 + 8*
     +    UG**(-2)*S*U1**(-1)*MS2**3 - 24*UG**(-2)*S*MS2*MG2 + 8*
     +    UG**(-2)*S*MS2**2 + 8*UG**(-2)*S*MG2**2 - 8*UG**(-2)*S**2*MS2
     +     + 8*UG**(-2)*S**2*MG2 + 20*UG**(-2)*T1*MS2*MG2 - 20*UG**(-2)
     +    *T1*MS2**2 + 16*UG**(-2)*T1**2*MS2 - 12*UG**(-2)*T1**2*MG2 - 
     +    4*UG**(-2)*T1**3 + 12*UG**(-1)*S*T1 - 20*UG**(-1)*S*MS2 + 8*
     +    UG**(-1)*S*MG2 + 8*UG**(-1)*S**2 )
     +
      MQQLRV = MQQLRV + SK1C0C(6,1)*N*CF**2 * ( 16*TG**(-2)*S*T1**(-1)*
     +    MS2*MG2**2 - 16*TG**(-2)*S*T1**(-1)*MS2**2*MG2 - 16*TG**(-2)*
     +    S*MS2*MG2 + 32*TG**(-2)*S*MG2**2 - 32*TG**(-2)*T1*MS2*MG2 + 
     +    32*TG**(-2)*T1*MG2**2 + 16*TG**(-1)*S*MG2 + 16*TG**(-1)*T1*
     +    MG2 - 16*UG**(-2)*S*T1*MS2 - 16*UG**(-2)*S*T1**2 - 16*
     +    UG**(-2)*T1*MS2*MG2 + 16*UG**(-2)*T1*MS2**2 - 16*UG**(-2)*
     +    T1**2*MG2 - 16*UG**(-2)*T1**3 )
     +
      MQQLRV = MQQLRV + SK1C0C(6,1)*N**2*CF * (  - 8*TG**(-2)*S*
     +    T1**(-1)*MS2*MG2**2 + 8*TG**(-2)*S*T1**(-1)*MS2**2*MG2 + 8*
     +    TG**(-2)*S*MS2*MG2 - 16*TG**(-2)*S*MG2**2 + 16*TG**(-2)*T1*
     +    MS2*MG2 - 16*TG**(-2)*T1*MG2**2 - 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*T1*MG2 + 4*UG**(-2)*S*T1*MS2 + 4*UG**(-2)*S*T1**2 + 
     +    4*UG**(-2)*T1*MS2*MG2 - 4*UG**(-2)*T1*MS2**2 + 4*UG**(-2)*
     +    T1**2*MG2 + 4*UG**(-2)*T1**3 )
     +
      MQQLRV = MQQLRV + SK1C0C(6,2)*N*CF**2 * ( 16*TG**(-1)*S*T1 - 32*
     +    TG**(-1)*S*MS2 + 16*TG**(-1)*S*MG2 + 16*TG**(-1)*S**2 - 32*
     +    TG**(-1)*T1*MS2 + 16*TG**(-1)*T1*MG2 + 16*UG**(-2)*S*T1*MG2
     +     + 16*UG**(-2)*S*U1**(-1)*MS2*MG2**2 - 16*UG**(-2)*S*U1**(-1)
     +    *MS2**2*MG2 + 16*UG**(-2)*S*MS2*MG2 + 16*UG**(-2)*T1*MS2*MG2
     +     - 16*UG**(-2)*T1*MG2**2 + 16*UG**(-2)*T1**2*MG2 + 16*S + 16*
     +    T1 )
     +
      MQQLRV = MQQLRV + SK1C0C(6,2)*N**2*CF * (  - 4*TG**(-1)*S*T1 + 8*
     +    TG**(-1)*S*MS2 - 4*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2 + 8*
     +    TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 - 8*UG**(-2)*S*T1*MG2 - 8
     +    *UG**(-2)*S*U1**(-1)*MS2*MG2**2 + 8*UG**(-2)*S*U1**(-1)*
     +    MS2**2*MG2 - 8*UG**(-2)*S*MS2*MG2 - 8*UG**(-2)*T1*MS2*MG2 + 8
     +    *UG**(-2)*T1*MG2**2 - 8*UG**(-2)*T1**2*MG2 - 4*S - 4*T1 )
     +
      MQQLRV = MQQLRV + SK1D0(4,1)*N*CF**2 * (  - 16*TG**(-1)*S*T1*MG2
     +     + 16*TG**(-1)*S*MS2*MG2 - 16*TG**(-1)*S*MS2**2 - 16*TG**(-1)
     +    *S**2*MG2 - 8*S*T1 + 8*S*MS2 + 8*S*MG2 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(4,1)*N**2*CF * ( 8*TG**(-1)*S*T1*MG2 - 8*
     +    TG**(-1)*S*MS2*MG2 + 8*TG**(-1)*S*MS2**2 + 8*TG**(-1)*S**2*
     +    MG2 + 4*S*T1 - 4*S*MS2 - 4*S*MG2 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(4,2)*N*CF**2 * ( 16*UG**(-2)*S*T1*MS2*MG2
     +     - 16*UG**(-2)*S*T1*MG2**2 + 8*UG**(-2)*S*T1**2*MS2 - 8*
     +    UG**(-2)*S*T1**2*MG2 + 8*UG**(-2)*S*MS2*MG2**2 + 8*UG**(-2)*S
     +    *MS2**2*MG2 - 8*UG**(-2)*S*MS2**3 - 8*UG**(-2)*S*MG2**3 + 8*
     +    UG**(-2)*S**2*T1*MS2 - 8*UG**(-2)*S**2*T1*MG2 + 8*UG**(-2)*
     +    S**2*MS2**2 - 8*UG**(-2)*S**2*MG2**2 + 8*UG**(-1)*S*T1*MS2 - 
     +    8*UG**(-1)*S*T1*MG2 - 8*UG**(-1)*S*T1**2 + 16*UG**(-1)*S*MS2*
     +    MG2 - 16*UG**(-1)*S*MG2**2 - 8*UG**(-1)*S**2*T1 - 8*UG**(-1)*
     +    S**2*MS2 - 8*UG**(-1)*S**2*MG2 )
     +
      MQQLRV = MQQLRV + SK1D0(4,2)*N**2*CF * (  - 8*UG**(-2)*S*T1*MS2*
     +    MG2 + 8*UG**(-2)*S*T1*MG2**2 - 4*UG**(-2)*S*T1**2*MS2 + 4*
     +    UG**(-2)*S*T1**2*MG2 - 4*UG**(-2)*S*MS2*MG2**2 - 4*UG**(-2)*S
     +    *MS2**2*MG2 + 4*UG**(-2)*S*MS2**3 + 4*UG**(-2)*S*MG2**3 - 4*
     +    UG**(-2)*S**2*T1*MS2 + 4*UG**(-2)*S**2*T1*MG2 - 4*UG**(-2)*
     +    S**2*MS2**2 + 4*UG**(-2)*S**2*MG2**2 - 4*UG**(-1)*S*T1*MS2 + 
     +    4*UG**(-1)*S*T1*MG2 + 4*UG**(-1)*S*T1**2 - 8*UG**(-1)*S*MS2*
     +    MG2 + 8*UG**(-1)*S*MG2**2 + 4*UG**(-1)*S**2*T1 + 4*UG**(-1)*
     +    S**2*MS2 + 4*UG**(-1)*S**2*MG2 )
     +
      MQQLRV = MQQLRV + SK1D0(7,1)*N*CF**2 * (  - 16*TG**(-1)*S*T1*MG2
     +     + 48*TG**(-1)*S*MS2*MG2 - 16*TG**(-1)*S*MS2**2 - 16*TG**(-1)
     +    *S**2*MG2 + 32*TG**(-1)*T1*MS2*MG2 - 32*TG**(-1)*T1*MS2**2 + 
     +    32*S*MS2 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(7,1)*N**2*CF * ( 8*TG**(-1)*S*T1*MG2 - 24
     +    *TG**(-1)*S*MS2*MG2 + 8*TG**(-1)*S*MS2**2 + 8*TG**(-1)*S**2*
     +    MG2 - 16*TG**(-1)*T1*MS2*MG2 + 16*TG**(-1)*T1*MS2**2 - 16*S*
     +    MS2 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(7,2)*N*CF**2 * ( 8*UG**(-2)*S*T1**2*MS2
     +     - 8*UG**(-2)*S*T1**2*MG2 - 56*UG**(-2)*S*MS2*MG2**2 + 88*
     +    UG**(-2)*S*MS2**2*MG2 - 40*UG**(-2)*S*MS2**3 + 8*UG**(-2)*S*
     +    MG2**3 + 8*UG**(-2)*S**2*T1*MS2 - 8*UG**(-2)*S**2*T1*MG2 - 48
     +    *UG**(-2)*S**2*MS2*MG2 + 40*UG**(-2)*S**2*MS2**2 + 8*UG**(-2)
     +    *S**2*MG2**2 + 32*UG**(-2)*T1*MS2*MG2**2 - 64*UG**(-2)*T1*
     +    MS2**2*MG2 + 32*UG**(-2)*T1*MS2**3 + 32*UG**(-2)*T1**2*MS2*
     +    MG2 - 32*UG**(-2)*T1**2*MS2**2 - 16*UG**(-1)*S*T1*MS2 - 8*
     +    UG**(-1)*S*T1**2 - 64*UG**(-1)*S*MS2*MG2 + 88*UG**(-1)*S*
     +    MS2**2 + 8*UG**(-1)*S*MG2**2 - 8*UG**(-1)*S**2*T1 - 32*
     +    UG**(-1)*S**2*MS2 - 8*S*T1 - 8*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(7,2)*N**2*CF * (  - 4*UG**(-2)*S*T1**2*
     +    MS2 + 4*UG**(-2)*S*T1**2*MG2 + 28*UG**(-2)*S*MS2*MG2**2 - 44*
     +    UG**(-2)*S*MS2**2*MG2 + 20*UG**(-2)*S*MS2**3 - 4*UG**(-2)*S*
     +    MG2**3 - 4*UG**(-2)*S**2*T1*MS2 + 4*UG**(-2)*S**2*T1*MG2 + 24
     +    *UG**(-2)*S**2*MS2*MG2 - 20*UG**(-2)*S**2*MS2**2 - 4*UG**(-2)
     +    *S**2*MG2**2 - 16*UG**(-2)*T1*MS2*MG2**2 + 32*UG**(-2)*T1*
     +    MS2**2*MG2 - 16*UG**(-2)*T1*MS2**3 - 16*UG**(-2)*T1**2*MS2*
     +    MG2 + 16*UG**(-2)*T1**2*MS2**2 + 8*UG**(-1)*S*T1*MS2 + 4*
     +    UG**(-1)*S*T1**2 + 32*UG**(-1)*S*MS2*MG2 - 44*UG**(-1)*S*
     +    MS2**2 - 4*UG**(-1)*S*MG2**2 + 4*UG**(-1)*S**2*T1 + 16*
     +    UG**(-1)*S**2*MS2 + 4*S*T1 + 4*S**2 )
     +
      MQQLRV = MQQLRV + SK1D0(8,1)*N*CF**2 * (  - 32*TG**(-1)*S*T1*MG2
     +     + 64*TG**(-1)*S*MS2*MG2 - 32*TG**(-1)*S*MS2**2 - 32*TG**(-1)
     +    *S*MG2**2 - 32*TG**(-1)*S**2*MG2 + 64*TG**(-1)*T1*MS2*MG2 - 
     +    32*TG**(-1)*T1*MS2**2 - 32*TG**(-1)*T1*MG2**2 - 32*S*T1 + 16*
     +    S*MS2 - 32*S*MG2 - 16*S**2 + 16*T1*MS2 - 32*T1*MG2 - 16*T1**2
     +     )
     +
      MQQLRV = MQQLRV + SK1D0(8,1)*N**2*CF * ( 8*TG**(-1)*S*T1*MG2 - 16
     +    *TG**(-1)*S*MS2*MG2 + 8*TG**(-1)*S*MS2**2 + 8*TG**(-1)*S*
     +    MG2**2 + 8*TG**(-1)*S**2*MG2 - 16*TG**(-1)*T1*MS2*MG2 + 8*
     +    TG**(-1)*T1*MS2**2 + 8*TG**(-1)*T1*MG2**2 + 8*S*T1 - 4*S*MS2
     +     + 8*S*MG2 + 4*S**2 - 4*T1*MS2 + 8*T1*MG2 + 4*T1**2 )
     +
      MQQLRV = MQQLRV + SK1D0(8,2)*N*CF**2 * ( 16*UG**(-2)*S*T1*MS2*MG2
     +     - 16*UG**(-2)*S*T1*MS2**2 + 16*UG**(-2)*S*T1**2*MG2 - 16*
     +    UG**(-2)*S*T1**3 + 16*UG**(-2)*T1*MS2*MG2**2 - 32*UG**(-2)*T1
     +    *MS2**2*MG2 + 16*UG**(-2)*T1*MS2**3 - 16*UG**(-2)*T1**2*
     +    MS2**2 + 16*UG**(-2)*T1**2*MG2**2 + 16*UG**(-2)*T1**3*MS2 - 
     +    16*UG**(-2)*T1**4 + 16*UG**(-1)*S*T1*MS2 + 16*UG**(-1)*S*
     +    T1**2 )
     +
      MQQLRV = MQQLRV + SK1D0(8,2)*N**2*CF * (  - 4*UG**(-2)*S*T1*MS2*
     +    MG2 + 4*UG**(-2)*S*T1*MS2**2 - 4*UG**(-2)*S*T1**2*MG2 + 4*
     +    UG**(-2)*S*T1**3 - 4*UG**(-2)*T1*MS2*MG2**2 + 8*UG**(-2)*T1*
     +    MS2**2*MG2 - 4*UG**(-2)*T1*MS2**3 + 4*UG**(-2)*T1**2*MS2**2
     +     - 4*UG**(-2)*T1**2*MG2**2 - 4*UG**(-2)*T1**3*MS2 + 4*
     +    UG**(-2)*T1**4 - 4*UG**(-1)*S*T1*MS2 - 4*UG**(-1)*S*T1**2 )
     +
      MQQLRV = MQQLRV + SOF1(1)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 - 
     +    128*TG**(-2)*S*MS2**2 + 64*TG**(-2)*S**2*MS2 + 128*TG**(-2)*
     +    T1*U1*MS2 - 64*UG**(-2)*S*T1*U1 - 128*UG**(-2)*S*MS2**2 + 64*
     +    UG**(-2)*S**2*MS2 + 128*UG**(-2)*T1*U1*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(1)*N**2*CF * ( 32*TG**(-2)*S*T1*U1 + 64*
     +    TG**(-2)*S*MS2**2 - 32*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1*U1*
     +    MS2 + 32*UG**(-2)*S*T1*U1 + 64*UG**(-2)*S*MS2**2 - 32*
     +    UG**(-2)*S**2*MS2 - 64*UG**(-2)*T1*U1*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(2)*N*CF**2 * ( 32*TG**(-2)*S*MS2**2 - 32*
     +    TG**(-2)*T1*U1*MS2 + 32*UG**(-2)*S*MS2**2 - 32*UG**(-2)*T1*U1
     +    *MS2 )
     +
      MQQLRV = MQQLRV + SOF1(3)*N*CF**2 * ( 32*TG**(-2)*S*MS2**2 - 32*
     +    TG**(-2)*T1*U1*MS2 + 32*UG**(-2)*S*MS2**2 - 32*UG**(-2)*T1*U1
     +    *MS2 )
     +
      MQQLRV = MQQLRV + SOF1(4)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S*T1*MS2 + 64*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1**2
     +    *U1 - 64*UG**(-2)*S*T1*U1 + 32*UG**(-2)*S*T1*MS2 + 32*
     +    UG**(-2)*S*T1*MG2 - 32*UG**(-2)*S*T1**2 + 32*UG**(-2)*S*MS2*
     +    MG2 - 32*UG**(-2)*S*MS2**2 + 64*UG**(-2)*S**2*MS2 - 64*
     +    UG**(-2)*T1**2*U1 - 32*UG**(-2)*T1**3 + 32*UG**(-1)*S*T1 + 32
     +    *UG**(-1)*S*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(4)*N**2*CF * ( 16*TG**(-2)*S*T1*U1 - 16*
     +    TG**(-2)*S*T1*MS2 - 16*TG**(-2)*S**2*MS2 + 16*TG**(-2)*T1**2*
     +    U1 + 16*UG**(-2)*S*T1*U1 - 16*UG**(-2)*S*T1*MS2 - 16*UG**(-2)
     +    *S**2*MS2 + 16*UG**(-2)*T1**2*U1 )
     +
      MQQLRV = MQQLRV + SOF1(5)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MS2 + 
     +    32*TG**(-2)*T1**2*U1 - 64*UG**(-2)*S*T1*MS2 - 32*UG**(-2)*S*
     +    T1**2 + 32*UG**(-2)*T1**2*U1 - 32*UG**(-2)*T1**3 )
     +
      MQQLRV = MQQLRV + SOF1(5)*N**2*CF * ( 16*TG**(-2)*S*T1*MS2 - 16*
     +    TG**(-2)*T1**2*U1 + 16*UG**(-2)*S*T1*MS2 - 16*UG**(-2)*T1**2*
     +    U1 )
     +
      MQQLRV = MQQLRV + SOF1(6)*N*CF**2 * (  - 32*TG**(-2)*S*T1*MS2 + 
     +    32*TG**(-2)*T1**2*U1 - 64*UG**(-2)*S*T1*MS2 - 32*UG**(-2)*S*
     +    T1**2 + 32*UG**(-2)*T1**2*U1 - 32*UG**(-2)*T1**3 )
     +
      MQQLRV = MQQLRV + SOF1(6)*N**2*CF * ( 16*TG**(-2)*S*T1*MS2 - 16*
     +    TG**(-2)*T1**2*U1 + 16*UG**(-2)*S*T1*MS2 - 16*UG**(-2)*T1**2*
     +    U1 )
     +
      MQQLRV = MQQLRV + SOF1(7)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S*T1*MS2 + 64*TG**(-2)*S**2*MS2 - 64*TG**(-2)*T1**2
     +    *U1 - 64*UG**(-2)*S*T1*U1 + 32*UG**(-2)*S*T1*MS2 + 32*
     +    UG**(-2)*S*T1*MG2 - 32*UG**(-2)*S*T1**2 + 32*UG**(-2)*S*MS2*
     +    MG2 - 32*UG**(-2)*S*MS2**2 + 64*UG**(-2)*S**2*MS2 - 64*
     +    UG**(-2)*T1**2*U1 - 32*UG**(-2)*T1**3 + 32*UG**(-1)*S*T1 + 32
     +    *UG**(-1)*S*MS2 )
     +
      MQQLRV = MQQLRV + SOF1(7)*N**2*CF * ( 16*TG**(-2)*S*T1*U1 - 16*
     +    TG**(-2)*S*T1*MS2 - 16*TG**(-2)*S**2*MS2 + 16*TG**(-2)*T1**2*
     +    U1 + 16*UG**(-2)*S*T1*U1 - 16*UG**(-2)*S*T1*MS2 - 16*UG**(-2)
     +    *S**2*MS2 + 16*UG**(-2)*T1**2*U1 )
     +
      MQQLRV = MQQLRV + SOF1(8)*N*CF**2 * (  - 64*TG**(-2)*S*T1*U1 + 64
     +    *TG**(-2)*S**2*MS2 - 64*UG**(-2)*S*T1*U1 + 64*UG**(-2)*S**2*
     +    MS2 )
     +
      MQQLRV = MQQLRV + SOF1(8)*N**2*CF * ( 32*TG**(-2)*S*T1*U1 - 32*
     +    TG**(-2)*S**2*MS2 + 32*UG**(-2)*S*T1*U1 - 32*UG**(-2)*S**2*
     +    MS2 )

         END IF

      END IF


      IF (IFL.EQ.0) MQQV = 2*MQPLLV + 2*MQPLRV 
      IF (IFL.EQ.1) MQQV = 2*0.5D0*MQQLLV + 1*MQQLRV

      DSSQQV = ALPHAS**3 * AVG * MQQV /4.D0 /S**2 *CONV
     +    + ALPHAS/PI*DSSQQB(ALPHAS,S,T1,MS,MG,IFL) 
     +    * (LOG(MG2/MS2)*(-1.D0) + LOG(MT2/MS2)*(-1.D0/3.D0))
     +    +  DSSQQR(ALPHAS,S,T1,MS,MG,IFL)
C***   CHANGES TO SUBTRACTED MSBAR AND TO DRBAR FOR THE YUKAWA-COUPLING 
      RETURN
      END


      REAL*8 FUNCTION DSSQQH(ALPHAS,S,T1,S4,MS,MG,IFL)
C***  HARD CROSS SECTIONS FOR Q + Q -> SQ + SQ + G
C***  IFL = 0 <--> Q NOT Q
C***  IFL = 1 <--> Q  =  Q
C***  EQUAL FLAVOR AND EQUAL HELICITY 
C***  HAS A FACTOR 1/2 FOR IDENTICAL PARTICLES 
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T,T1,TG,U,U1,UG,SB,S1,MS,MG,MS2,MG2,M2
      REAL*8 MQPLLH, MQPLRH, MQQLLH, MQQLRH, MQQH
      REAL*8 NS,CONV,PI,ZETA2,N,CF,AVG,BETA,XS
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 AHP1P1, ABP1P1, ABP1M1, ABP1P2, A4P2P2, A4P1P2, A4P1P1
      REAL*8 A4M1P2, A4M1P1, A4P0P2, A4P1P0
      REAL*8 ANG4(1:104), COLO1(1:9)
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      N = 3.D0
      CF = (N**2 -1)/2.D0/N
      AVG = (1.D0/2.D0)**2 * (1.D0/3.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -T1
      TG = T1 -M2
      UG = U1 -M2
      SB = S*SQRT(1 -4.D0*MS2/S)
      S1 = 4.D0*MS2 - S
      T = T1 +MS2
      U = U1 +MS2
      BETA = SQRT(1 -4*MS2/S)
      XS = (1.D0 -BETA)/(1.D0 +BETA)


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +U1)/ANGDEF(1)
      ANGDEF(3) = (S +T1)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(T1 +U1 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(7) = SQRT((T1 +U1)**2 -4.D0*MS2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (T1*S4 -S*(U1+2.D0*MS2))/(S+T1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (U1*S4 -S*(T1+2.D0*MS2))/(S+U1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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

C$$$      ANG4(1) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(2) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(3) =A4P1P2(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(4) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(5) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(6) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(7) =A4P1P2(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(8) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(9) =A4P2M2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(10) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))

C$$$      ANG4(11) =A4M2P1(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(12) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(13) =A4P2M2(ANA(2,2),ANB(2,2),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(14) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(15) =A4M2P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(16) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(17) =A4P1P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(18) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(19) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(20) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,5),ANB(3,5),ANC(3,5))

C$$$      ANG4(21) =A4M1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(22) =A4M1P2(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(23) =A4M2P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(24) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(25) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(26) =A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG4(27) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(28) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(29) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(30) =A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))

C$$$      ANG4(31) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(32) =A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(33) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(34) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(35) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
       ANG4(36) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(37) =A4P1P0(ANA(2,1),ANB(2,1),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(38) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(39) =A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(40) =A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))

C$$$      ANG4(41) =A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG4(43) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(44) =A4P1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
 
C$$$      ANG4(47) =A4P1P0(ANA(2,7),ANB(2,7),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(48) =AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))
       ANG4(49) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(50) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1))

C$$$      ANG4(51) =ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(52) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(53) =ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(54) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(55) =ABP2P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(56) =ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(57) =ABP2P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(58) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(59) =ABP2P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(60) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))

C$$$      ANG4(61) =ABP2P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(62) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
       ANG4(63) =ABP1M1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(64) =ABP2P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(65) =ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(66) =ABP2P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(67) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
       ANG4(68) =ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(69) =ABP1M2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(70) =ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG4(71) =ABP2P0(ANA(3,4),ANB(3,4),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(72) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
       ANG4(73) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
       ANG4(74) =A4P1P2(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
       ANG4(75) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(76) =A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
       ANG4(77) =A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
       ANG4(78) =A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
       ANG4(79) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,9),ANB(1,9),ANC(1,9))
       ANG4(80) =A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))

       ANG4(81) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9)) 
       ANG4(82) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9)) 
C$$$      ANG4(83) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG4(84) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
       ANG4(85) = ABP1P2(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8)) 
       ANG4(86) = ABP1P1(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8)) 
       ANG4(87) = A4P2P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
       ANG4(88) = A4P1P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
       ANG4(89) = A4P1P2(ANA(1,8),ANB(1,8),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(90) = A4P1P1(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))

       ANG4(91) = A4P0P2(ANA(1,1),ANB(1,1),ANA(2,8),ANB(2,8),ANC(2,8))
       ANG4(92) = A4P1P0(ANA(1,8),ANB(1,8),ANA(1,1),ANB(1,1),ANC(1,1))
       ANG4(93) = A4M1P1(ANA(3,6),ANB(3,6),ANA(3,8),ANB(3,8),ANC(3,8))
       ANG4(94) = A4P1P1(ANA(1,8),ANB(1,8),ANA(1,9),ANB(1,9),ANC(1,9))
       ANG4(95) = A4M1P2(ANA(3,6),ANB(3,6),ANA(3,8),ANB(3,8),ANC(3,8))
C$$$      ANG4(96) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
C$$$      ANG4(97) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
C$$$      ANG4(98) =A4P1P2(ANA(1,8),ANB(1,8),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG4(99) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))

C$$$      ANG4(100)=A4P2P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(101)=A4P1P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(102)=A4P1P2(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(103)=A4P1P1(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(104)=A4P1P1(ANA(1,8),ANB(1,8),ANA(1,6),ANB(1,6),ANC(1,6))


      COLO1(9) = LOG(S4**2/MS2/(S4+MS2))

      IF (IFL.EQ.0) THEN

      MQPLLH = 0.D0
      MQPLLH = MQPLLH + N*CF**2*(S4+MS2) * ( 8*M2*TG**(-2)*S4**(-1)*
     +    (S+U1)**(-1)*MG2 + 8*M2*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*
     +    MG2 + 8*M2*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*
     +    TG**(-2)*S*S4**(-1)*(S+U1)**(-1)*MG2 - 8*TG**(-2)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2)*MG2 + 8*TG**(-1)*S4**(-1)*(S+U1)**(-1)*
     +    MG2 + 8*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 - 8*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + N*CF**2 * ( 16*M2*TG**(-2)*S4**(-2)*MS2*MG2 + 
     +    16*M2*TG**(-2)*S4**(-1)*MG2 + 16*TG**(-2)*U1*S4**(-2)*MS2*MG2
     +     + 8*TG**(-2)*U1*S4**(-1)*MG2 - 16*TG**(-2)*S4**(-1)*MS2*MG2
     +     - 8*TG**(-2)*MG2 + 16*TG**(-1)*S4**(-2)*MS2*MG2 + 16*
     +    TG**(-1)*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(48)*N*CF**2 * ( 16*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(48)*N**2*CF * (  - 8*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(49)*N*CF**2 * ( 4*M2*TG**(-2)*S*MG2 - 4*M2
     +    *TG**(-2)*U1*MG2 - 8*M2*TG**(-1)*MG2 - 4*M2**2*TG**(-2)*MG2
     +     + 4*TG**(-1)*S*MG2 - 4*TG**(-1)*U1*MG2 - 4*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(49)*N**2*CF * (  - 2*M2*TG**(-2)*S*MG2 + 2
     +    *M2*TG**(-2)*U1*MG2 + 4*M2*TG**(-1)*MG2 + 2*M2**2*TG**(-2)*
     +    MG2 - 2*TG**(-1)*S*MG2 + 2*TG**(-1)*U1*MG2 + 2*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(50)*N*CF**2 * ( 16*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 - 8*M2*(S+U1+M2)**(-1)*MG2 - 16*TG**(-1)*
     +    S*MG2 + 16*TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 8*TG**(-1)*U1*
     +    MG2 - 8*S*(S+U1+M2)**(-1)*MG2 + 8*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(50)*N**2*CF * (  - 4*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*(S+U1+M2)**(-1)*MG2 + 4*TG**(-1)*S
     +    *MG2 - 4*TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 - 2*TG**(-1)*U1*
     +    MG2 + 2*S*(S+U1+M2)**(-1)*MG2 - 2*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(68)*N*CF**2 * (  - 4*M2*TG**(-2)*S4**(-1)*
     +    MG2 - 4*TG**(-1)*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(68)*N**2*CF * ( 2*M2*TG**(-2)*S4**(-1)*MG2
     +     - 2*M2*TG**(-2)*(S+U1+M2)**(-1)*MG2 + 2*TG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQPLLH = MQPLLH + ANG4(72)*N*CF**2 * (  - 8*S*MS2*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(73)*N*CF**2 * (  - 4*M2*TG**(-1)*S*MG2 + 4
     +    *M2*TG**(-1)*U1*MG2 + 4*M2*MG2 + 4*M2**2*TG**(-1)*MG2 - 4*S*
     +    MG2 )
     +
      MQPLLH = MQPLLH + ANG4(73)*N**2*CF * (  - 8*TG**(-1)*S*MS2*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(75)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MG2 + 4
     +    *M2*TG**(-2)*U1*MG2 + 16*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 
     +    32*M2*TG**(-1)*U1*S4**(-1)*MG2 - 32*M2*TG**(-1)*S4**(-1)*MS2*
     +    MG2 + 16*M2*TG**(-1)*MG2 - 16*M2*S4**(-1)*MG2 - 8*M2*
     +    (S+U1+M2)**(-1)*MG2 + 4*M2**2*TG**(-2)*MG2 - 16*M2**2*
     +    TG**(-1)*S4**(-1)*MG2 - 20*TG**(-1)*S*MG2 + 16*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 - 32*TG**(-1)*U1*S4**(-1)*MS2*MG2 + 20*
     +    TG**(-1)*U1*MG2 - 16*TG**(-1)*U1**2*S4**(-1)*MG2 + 32*
     +    TG**(-1)*MS2*MG2 - 8*S*(S+U1+M2)**(-1)*MG2 - 16*T1*S4**(-1)*
     +    MG2 - 32*U1*S4**(-1)*MG2 - 32*S4**(-1)*MS2*MG2 + 20*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(75)*N**2*CF * ( 2*M2*TG**(-2)*S*MG2 - 2*M2
     +    *TG**(-2)*U1*MG2 - 4*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 16*
     +    M2*TG**(-1)*U1*S4**(-1)*MG2 + 16*M2*TG**(-1)*S4**(-1)*MS2*MG2
     +     - 8*M2*TG**(-1)*MG2 + 8*M2*S4**(-1)*MG2 + 2*M2*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2**2*TG**(-2)*MG2 + 8*M2**2*TG**(-1)
     +    *S4**(-1)*MG2 + 6*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 + 16*TG**(-1)*U1*S4**(-1)*MS2*MG2 - 8*
     +    TG**(-1)*U1*MG2 + 8*TG**(-1)*U1**2*S4**(-1)*MG2 - 16*TG**(-1)
     +    *MS2*MG2 + 2*S*(S+U1+M2)**(-1)*MG2 + 8*T1*S4**(-1)*MG2 + 16*
     +    U1*S4**(-1)*MG2 + 16*S4**(-1)*MS2*MG2 - 8*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(77)*N**2*CF * ( 4*M2*TG**(-2)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(78)*N**2*CF * ( 2*M2*TG**(-2)*S4**(-1)*MG2
     +     - 2*M2*TG**(-2)*(S+U1+M2)**(-1)*MG2 + 2*TG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQPLLH = MQPLLH + ANG4(79)*N**2*CF * ( 4*M2*TG**(-2)*U1*MG2 + 4*
     +    M2*TG**(-1)*MG2 + 8*M2**2*TG**(-2)*MG2 - 8*TG**(-2)*S*MS2*MG2
     +     - 4*TG**(-1)*S*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(80)*N*CF**2 * (  - 8*M2*TG**(-1)*S4**(-1)*
     +    MG2 - 8*M2*TG**(-1)*(S+U1+M2)**(-1)*MG2 - 8*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(80)*N**2*CF * ( 2*M2*TG**(-2)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*TG**(-2)*U1*S4**(-1)*MG2 + 8*M2*
     +    TG**(-2)*S4**(-1)*MS2*MG2 + 2*M2*TG**(-2)*MG2 + 10*M2*
     +    TG**(-1)*S4**(-1)*MG2 + 4*M2**2*TG**(-2)*S4**(-1)*MG2 - 2*
     +    M2**2*TG**(-2)*(S+U1+M2)**(-1)*MG2 - 4*TG**(-2)*S*MG2 + 8*
     +    TG**(-2)*U1*S4**(-1)*MS2*MG2 - 8*TG**(-2)*MS2*MG2 - 2*
     +    TG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 4*TG**(-1)*U1*S4**(-1)*MG2
     +     + 8*TG**(-1)*S4**(-1)*MS2*MG2 + 6*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(81)*N*CF**2 * (  - 4*M2*TG**(-1)*S*MG2 + 4
     +    *M2*TG**(-1)*U1*MG2 + 4*M2**2*TG**(-1)*MG2 - 4*S*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(82)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MG2 + 4
     +    *M2*TG**(-2)*U1*MG2 - 16*M2*TG**(-1)*U1*S4**(-1)*MG2 + 16*M2*
     +    TG**(-1)*MG2 + 4*M2**2*TG**(-2)*MG2 - 20*TG**(-1)*S*MG2 + 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 20*TG**(-1)*U1*MG2 - 16*
     +    TG**(-1)*U1**2*S4**(-1)*MG2 - 8*U1*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + ANG4(82)*N**2*CF * ( 2*M2*TG**(-2)*S*MG2 - 2*M2
     +    *TG**(-2)*U1*MG2 + 4*M2*TG**(-1)*U1*S4**(-1)*MG2 - 4*M2*
     +    TG**(-1)*MG2 - 2*M2**2*TG**(-2)*MG2 + 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 - 6*TG**(-1)*U1*MG2 + 4*
     +    TG**(-1)*U1**2*S4**(-1)*MG2 + 2*U1*S4**(-1)*MG2 )
     +
      MQPLLH = MQPLLH + COLO1(9)*N*CF**2*(S4+MS2) * ( 8*TG**(-2)*S*
     +    T1**2*S4**(-2)*(S+U1)**(-2)*MG2 + 8*TG**(-2)*S*S4**(-2)*MG2
     +     + 8*S*T1**2*U1**2*S4**(-2)*(S+T1)**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*T1**2*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*S**2*T1*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*S**2*T1*
     +    S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S**3*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S**3*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 )

      MQPLRH = 0.D0
      MQPLRH = MQPLRH + N*CF**2*(S4+MS2) * ( 8*M2*TG**(-2)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2) + 8*M2*TG**(-2)*U1*S4**(-1)*
     +    (S+U1)**(-1) - 8*M2*TG**(-2)*S4**(-1)*(S+U1)**(-1)*MS2 - 32*
     +    M2*TG**(-1)*S4**(-2) + 16*M2*TG**(-1)*S4**(-1)*(S+T1)**(-1)
     +     - 16*M2*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) - 8*M2*S*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 16*M2*S*S4**(-1)*
     +    (S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) - 8*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*M2*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*M2*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*M2*T1**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 32*M2*S4**(-2)*(S+U1+M2)**(-1) - 16
     +    *M2*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*MS2 + 8*TG**(-2)*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*MS2 - 32*TG**(-1)*S*S4**(-2) + 16*TG**(-1)*S*
     +    S4**(-1)*(S+T1)**(-1) )
     +
      MQPLRH = MQPLRH + N*CF**2*(S4+MS2) * ( 16*TG**(-1)*S*S4**(-1)*
     +    (S+U1)**(-1) + 8*TG**(-1)*T1*U1*S4**(-1)*(S+U1)**(-2) - 32*
     +    TG**(-1)*U1*S4**(-2) + 24*TG**(-1)*U1*S4**(-1)*(S+U1)**(-1)
     +     - 8*TG**(-1)*S4**(-1)*(S+U1)**(-1)*MS2 - 16*S*T1*S4**(-1)*
     +    (S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) - 8*S*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 + 8*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) - 8*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 32*T1*S4**(-2)*(S+U1+M2)**(-1)
     +     + 24*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 8*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) - 64*S4**(-2) + 16*S4**(-1)*
     +    (S+T1)**(-1) + 16*S4**(-1)*(S+U1)**(-1) + 8*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + N*CF**2 * (  - 16*M2*TG**(-2)*U1*S4**(-2)*MS2
     +     - 8*M2*TG**(-2)*U1*S4**(-1) - 16*M2*TG**(-2)*S4**(-2)*MS2**2
     +     - 16*TG**(-2)*U1*S4**(-2)*MS2**2 - 8*TG**(-2)*U1*S4**(-1)*
     +    MS2 + 16*TG**(-2)*S4**(-1)*MS2**2 + 8*TG**(-2)*MS2 - 16*
     +    TG**(-1)*U1*S4**(-2)*MS2 - 8*TG**(-1)*U1*S4**(-1) - 16*
     +    TG**(-1)*S4**(-2)*MS2**2 )
     +
      MQPLRH = MQPLRH + N**2*CF*(S4+MS2) * (  - 4*M2*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 4*M2*TG**(-2)*T1*
     +    U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1)*MS2 - 4*M2*TG**(-2)*
     +    T1*U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1)*MG2 - 4*M2*
     +    TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-2) + 4*M2*TG**(-2)*S4**(-1)
     +    *(S+U1)**(-1)*MS2 - 4*M2*TG**(-2)*S4**(-1)*(S+U1)**(-1)*MG2
     +     + 4*M2*TG**(-1)*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 
     +    4*M2*TG**(-1)*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*
     +    M2*TG**(-1)*S*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*
     +    TG**(-1)*S*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 - 8*M2*
     +    TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*M2*
     +    TG**(-1)*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*M2*TG**(-1)*
     +    S**2*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*TG**(-1)*
     +    S**2*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2*TG**(-1)*T1*
     +    U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1) )
     +
      MQPLRH = MQPLRH + N**2*CF*(S4+MS2) * ( 4*M2*TG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*TG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 - 4*M2*TG**(-1)*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 16*M2*TG**(-1)*S4**(-2) - 
     +    4*M2*TG**(-1)*S4**(-1)*(S+T1)**(-1) - 4*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 8*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1) - 4*M2*TG**(-1)*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)
     +    *MS2 + 4*M2*TG**(-1)*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MG2 + 4
     +    *M2*S*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) - 8*M2*S*S4**(-1)*
     +    (S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 )
     +
      MQPLRH = MQPLRH + N**2*CF*(S4+MS2) * (  - 16*M2*S4**(-2)*
     +    (S+U1+M2)**(-1) + 4*M2*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)
     +     - 8*M2*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) - 4*M2**2*TG**(-2)*S
     +    *S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) + 4*M2**2*TG**(-2)*T1*
     +    U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1) - 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1) + 4*M2**2*TG**(-1)*S*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2**2*TG**(-1)*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2**2*TG**(-1)*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2**2*TG**(-1)*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 )
     +
      MQPLRH = MQPLRH + N**2*CF*(S4+MS2) * (  - 8*M2**2*TG**(-1)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*M2**2*TG**(-1)*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*M2**2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 4*M2**3*TG**(-2)*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1) + 4*M2**3*TG**(-1)*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**3*TG**(-1)*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 4*TG**(-2)*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*MS2 + 4*TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-2)*MG2
     +     + 16*TG**(-1)*S*S4**(-2) - 4*TG**(-1)*S*S4**(-1)*
     +    (S+T1)**(-1) - 4*TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1)*MS2 + 4*TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 8*TG**(-1)*S*S4**(-1)*(S+U1)**(-1) - 4*
     +    TG**(-1)*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MS2 + 4*TG**(-1)*
     +    S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MG2 - 4*TG**(-1)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2) )
     +
      MQPLRH = MQPLRH + N**2*CF*(S4+MS2) * ( 16*TG**(-1)*U1*S4**(-2) - 
     +    8*TG**(-1)*U1*S4**(-1)*(S+U1)**(-1) + 4*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*MS2 - 4*TG**(-1)*S4**(-1)*(S+U1)**(-1)*MG2 + 8*S
     +    *T1*S4**(-1)*(S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*S*U1
     +    *S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 - 4*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1) + 4*T1*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*
     +    MS2 - 4*T1*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*T1*
     +    S4**(-2)*(S+U1+M2)**(-1) - 8*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) - 4*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) + 32*S4**(-2) - 4*S4**(-1)*
     +    (S+T1)**(-1) - 4*S4**(-1)*(S+U1)**(-1) - 4*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MS2 + 4*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 )
     +
      MQPLRH = MQPLRH + ANG4(36)*N*CF**2 * ( 8*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(37)*N*CF**2 * (  - 16*M2*TG**(-1)*S4**(-1)
     +    *MS2 - 8*M2*S4**(-1) - 8*M2**2*TG**(-1)*S4**(-1) + 16*
     +    TG**(-1)*U1*S4**(-1)*MS2 + 8*TG**(-1)*U1**2*S4**(-1) - 16*
     +    TG**(-1)*MS2 - 8*T1*S4**(-1) - 16*S4**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(37)*N**2*CF * ( 8*M2*TG**(-1)*S4**(-1)*MS2
     +     + 4*M2*S4**(-1) + 4*M2**2*TG**(-1)*S4**(-1) - 8*TG**(-1)*U1*
     +    S4**(-1)*MS2 - 4*TG**(-1)*U1**2*S4**(-1) + 8*TG**(-1)*MS2 + 4
     +    *T1*S4**(-1) + 8*S4**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(48)*N*CF**2 * ( 8*M2*TG**(-1)*S - 16*M2*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1) - 8*M2*S*(S+U1+M2)**(-1) - 16*
     +    M2**2*TG**(-1)*S*(S+U1+M2)**(-1) - 8*TG**(-1)*S*U1 - 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 8*S*T1*(S+U1+M2)**(-1) - 
     +    8*S**2*(S+U1+M2)**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(48)*N**2*CF * (  - 4*M2*TG**(-1)*S + 8*M2*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1) + 4*M2*S*(S+U1+M2)**(-1) + 8*
     +    M2**2*TG**(-1)*S*(S+U1+M2)**(-1) + 4*TG**(-1)*S*U1 + 8*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 4*S*T1*(S+U1+M2)**(-1) + 
     +    4*S**2*(S+U1+M2)**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(49)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MS2 + 4
     +    *M2*TG**(-2)*U1*MS2 - 8*M2*TG**(-1)*S + 8*M2*TG**(-1)*U1 + 8*
     +    M2*TG**(-1)*MS2 - 8*M2 - 4*M2**2*TG**(-2)*S + 4*M2**2*
     +    TG**(-2)*U1 + 4*M2**2*TG**(-2)*MS2 - 12*M2**2*TG**(-1) - 4*
     +    M2**3*TG**(-2) - 4*TG**(-1)*S*MS2 + 4*TG**(-1)*U1*MS2 - 4*S
     +     - 4*T1 + 4*U1 + 4*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(49)*N**2*CF * ( 2*M2*TG**(-2)*S*MS2 - 2*M2
     +    *TG**(-2)*U1*MS2 + 4*M2*TG**(-1)*S - 4*M2*TG**(-1)*U1 - 4*M2*
     +    TG**(-1)*MS2 + 4*M2 + 2*M2**2*TG**(-2)*S - 2*M2**2*TG**(-2)*
     +    U1 - 2*M2**2*TG**(-2)*MS2 + 6*M2**2*TG**(-1) + 2*M2**3*
     +    TG**(-2) + 2*TG**(-1)*S*MS2 - 2*TG**(-1)*U1*MS2 + 2*S + 2*T1
     +     - 2*U1 - 2*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(50)*N*CF**2 * (  - 16*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 + 16*M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) - 16*M2*TG**(-1)*U1 - 16*M2*S*(S+U1+M2)**(-1)
     +     + 8*M2*(S+U1+M2)**(-1)*MS2 + 8*M2 - 32*M2**2*TG**(-1)*S*
     +    (S+U1+M2)**(-1) + 16*M2**2*TG**(-1) - 8*M2**2*(S+U1+M2)**(-1)
     +     - 16*M2**3*TG**(-1)*(S+U1+M2)**(-1) + 16*TG**(-1)*S*MS2 - 16
     +    *TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 8*TG**(-1)*U1*MS2 + 8*
     +    TG**(-1)*U1**2 + 8*S*(S+U1+M2)**(-1)*MS2 + 8*S - 8*S**2*
     +    (S+U1+M2)**(-1) - 8*U1 - 8*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(50)*N**2*CF * ( 4*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 - 4*M2*TG**(-1)*S + 4*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) + 4*M2*TG**(-1)*U1 + 4*M2*S*(S+U1+M2)**(-1)
     +     - 2*M2*(S+U1+M2)**(-1)*MS2 - 2*M2 + 8*M2**2*TG**(-1)*S*
     +    (S+U1+M2)**(-1) - 4*M2**2*TG**(-1) + 2*M2**2*(S+U1+M2)**(-1)
     +     + 4*M2**3*TG**(-1)*(S+U1+M2)**(-1) - 4*TG**(-1)*S*MS2 + 4*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 2*TG**(-1)*U1*MS2 - 2*
     +    TG**(-1)*U1**2 - 2*S*(S+U1+M2)**(-1)*MS2 - 2*S + 2*S**2*
     +    (S+U1+M2)**(-1) + 2*U1 + 2*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(68)*N*CF**2 * ( 4*M2*TG**(-2)*S4**(-1)*MS2
     +     + 8*M2*TG**(-1)*S4**(-1) + 4*M2**2*TG**(-2)*S4**(-1) + 4*
     +    TG**(-1)*S4**(-1)*MS2 + 4*S4**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(68)*N**2*CF * (  - 2*M2*TG**(-2)*S4**(-1)*
     +    MS2 + 4*M2*TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*TG**(-2) - 4*M2*TG**(-1)*S4**(-1)
     +     + 2*M2*TG**(-1)*(S+U1+M2)**(-1) - 2*M2**2*TG**(-2)*S4**(-1)
     +     + 4*M2**2*TG**(-2)*(S+U1+M2)**(-1) - 2*TG**(-2)*MS2 + 2*
     +    TG**(-2)*MG2 - 2*TG**(-1)*S4**(-1)*MS2 - 2*S4**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(72)*N*CF**2 * ( 8*M2*S*MS2 + 8*M2**2*MS2
     +     + 8*S*MS2**2 )
     +
      MQPLRH = MQPLRH + ANG4(73)*N*CF**2 * ( 4*M2*TG**(-1)*S*MS2 - 4*M2
     +    *TG**(-1)*U1*MS2 + 4*M2*S - 4*M2*MS2 + 4*M2**2*TG**(-1)*S - 4
     +    *M2**2*TG**(-1)*U1 - 4*M2**2*TG**(-1)*MS2 + 4*M2**2 + 4*M2**3
     +    *TG**(-1) + 4*S*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(73)*N**2*CF * ( 2*M2*TG**(-1)*S*MS2 + 2*M2
     +    *TG**(-1)*S*MG2 + 2*M2*TG**(-1)*U1*MS2 - 2*M2*TG**(-1)*U1*MG2
     +     - 2*M2*S - 2*M2*MS2 + 2*M2*MG2 - 2*M2**2*TG**(-1)*S + 2*
     +    M2**2*TG**(-1)*U1 + 6*M2**2*TG**(-1)*MS2 + 2*M2**2*TG**(-1)*
     +    MG2 - 2*M2**2 - 2*M2**3*TG**(-1) + 4*TG**(-1)*S*MS2*MG2 + 4*
     +    TG**(-1)*S*MS2**2 - 2*S*MS2 + 2*S*MG2 )
     +
      MQPLRH = MQPLRH + ANG4(74)*N*CF**2 * ( 16*M2*MS2 + 8*S*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(75)*N*CF**2 * ( 4*M2*TG**(-2)*S*MS2 - 4*M2
     +    *TG**(-2)*U1*MS2 - 16*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MS2 + 24*
     +    M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*(S+U1+M2)**(-1) + 64*M2*
     +    TG**(-1)*U1*S4**(-1)*MS2 - 24*M2*TG**(-1)*U1 + 16*M2*TG**(-1)
     +    *U1**2*S4**(-1) + 32*M2*TG**(-1)*S4**(-1)*MS2**2 - 48*M2*
     +    TG**(-1)*MS2 - 16*M2*S*(S+U1+M2)**(-1) + 16*M2*U1*S4**(-1) + 
     +    16*M2*S4**(-1)*MS2 + 8*M2*(S+U1+M2)**(-1)*MS2 + 8*M2 + 4*
     +    M2**2*TG**(-2)*S - 4*M2**2*TG**(-2)*U1 - 4*M2**2*TG**(-2)*MS2
     +     - 32*M2**2*TG**(-1)*S*(S+U1+M2)**(-1) + 16*M2**2*TG**(-1)*U1
     +    *S4**(-1) + 16*M2**2*TG**(-1)*S4**(-1)*MS2 + 20*M2**2*
     +    TG**(-1) - 8*M2**2*(S+U1+M2)**(-1) + 4*M2**3*TG**(-2) - 16*
     +    M2**3*TG**(-1)*(S+U1+M2)**(-1) + 20*TG**(-1)*S*MS2 - 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 32*TG**(-1)*U1*S4**(-1)*
     +    MS2**2 - 20*TG**(-1)*U1*MS2 + 16*TG**(-1)*U1**2*S4**(-1)*MS2
     +     - 32*TG**(-1)*MS2**2 + 8*S*(S+U1+M2)**(-1)*MS2 + 12*S - 8*
     +    S**2*(S+U1+M2)**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(75)*N*CF**2 * ( 16*T1*U1*S4**(-1) + 32*T1*
     +    S4**(-1)*MS2 - 8*T1 + 8*T1**2*S4**(-1) + 48*U1*S4**(-1)*MS2
     +     - 8*U1 + 8*U1**2*S4**(-1) + 32*S4**(-1)*MS2**2 - 36*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(75)*N**2*CF * (  - 2*M2*TG**(-2)*S*MS2 + 2
     +    *M2*TG**(-2)*U1*MS2 + 4*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MS2 - 8
     +    *M2*TG**(-1)*S + 4*M2*TG**(-1)*S**2*(S+U1+M2)**(-1) - 32*M2*
     +    TG**(-1)*U1*S4**(-1)*MS2 + 10*M2*TG**(-1)*U1 - 8*M2*TG**(-1)*
     +    U1**2*S4**(-1) - 16*M2*TG**(-1)*S4**(-1)*MS2**2 + 30*M2*
     +    TG**(-1)*MS2 + 2*M2*TG**(-1)*MG2 + 4*M2*S*(S+U1+M2)**(-1) - 8
     +    *M2*U1*S4**(-1) - 8*M2*S4**(-1)*MS2 - 2*M2*(S+U1+M2)**(-1)*
     +    MS2 - 2*M2 - 2*M2**2*TG**(-2)*S + 2*M2**2*TG**(-2)*U1 + 2*
     +    M2**2*TG**(-2)*MS2 + 8*M2**2*TG**(-1)*S*(S+U1+M2)**(-1) - 8*
     +    M2**2*TG**(-1)*U1*S4**(-1) - 8*M2**2*TG**(-1)*S4**(-1)*MS2 - 
     +    6*M2**2*TG**(-1) + 2*M2**2*(S+U1+M2)**(-1) - 2*M2**3*TG**(-2)
     +     + 4*M2**3*TG**(-1)*(S+U1+M2)**(-1) - 2*TG**(-1)*S*MS2 + 4*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 16*TG**(-1)*U1*S4**(-1)*
     +    MS2**2 + 10*TG**(-1)*U1*MS2 - 2*TG**(-1)*U1*MG2 - 8*TG**(-1)*
     +    U1**2*S4**(-1)*MS2 + 16*TG**(-1)*MS2**2 - 2*S*(S+U1+M2)**(-1)
     +    *MS2 )
     +
      MQPLRH = MQPLRH + ANG4(75)*N**2*CF * (  - 4*S + 2*S**2*
     +    (S+U1+M2)**(-1) - 8*T1*U1*S4**(-1) - 16*T1*S4**(-1)*MS2 + 4*
     +    T1 - 4*T1**2*S4**(-1) - 24*U1*S4**(-1)*MS2 + 4*U1 - 4*U1**2*
     +    S4**(-1) - 16*S4**(-1)*MS2**2 + 14*MS2 + 2*MG2 )
     +
      MQPLRH = MQPLRH + ANG4(77)*N**2*CF * ( 4*M2*TG**(-2)*S + 4*M2*
     +    TG**(-2)*U1 - 4*M2*TG**(-2)*MS2 - 4*M2*TG**(-1) - 4*M2**2*
     +    TG**(-2) + 4*TG**(-2)*S*MS2 - 4*TG**(-2)*S*MG2 + 4*TG**(-2)*
     +    U1*MS2 - 4*TG**(-2)*U1*MG2 )
     +
      MQPLRH = MQPLRH + ANG4(78)*N**2*CF * (  - 2*M2*TG**(-2)*S4**(-1)*
     +    MS2 + 4*M2*TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*TG**(-2) - 4*M2*TG**(-1)*S4**(-1)
     +     + 2*M2*TG**(-1)*(S+U1+M2)**(-1) - 2*M2**2*TG**(-2)*S4**(-1)
     +     + 4*M2**2*TG**(-2)*(S+U1+M2)**(-1) - 2*TG**(-2)*MS2 + 2*
     +    TG**(-2)*MG2 - 2*TG**(-1)*S4**(-1)*MS2 - 2*S4**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(79)*N**2*CF * ( 4*M2*TG**(-2)*S*U1 + 8*M2*
     +    TG**(-2)*S*MS2 - 8*M2*TG**(-2)*S*MG2 + 4*M2*TG**(-2)*S**2 + 4
     +    *M2*TG**(-2)*U1*MS2 - 8*M2*TG**(-2)*U1*MG2 - 2*M2*TG**(-1)*S
     +     + 2*M2*TG**(-1)*U1 - 8*M2*TG**(-1)*MS2 + 4*M2*TG**(-1)*MG2
     +     + 4*M2**2*TG**(-2)*S + 4*M2**2*TG**(-2)*U1 - 4*M2**2*
     +    TG**(-2)*MS2 + 4*M2**2*TG**(-2)*MG2 - 4*M2**2*TG**(-1) - 4*
     +    M2**3*TG**(-2) + 4*TG**(-2)*S*U1*MS2 - 4*TG**(-2)*S*U1*MG2 + 
     +    4*TG**(-2)*S*MS2**2 + 4*TG**(-2)*S*MG2**2 + 4*TG**(-2)*S**2*
     +    MS2 - 4*TG**(-2)*S**2*MG2 + 2*TG**(-1)*S*MS2 + 2*TG**(-1)*S*
     +    MG2 + 2*TG**(-1)*U1*MS2 - 2*TG**(-1)*U1*MG2 )
     +
      MQPLRH = MQPLRH + ANG4(80)*N*CF**2 * ( 8*M2*TG**(-1)*S4**(-1)*MS2
     +     + 8*M2*TG**(-1)*(S+U1+M2)**(-1)*MS2 + 8*M2*S4**(-1) + 8*M2*
     +    (S+U1+M2)**(-1) + 8*M2**2*TG**(-1)*S4**(-1) + 8*M2**2*
     +    TG**(-1)*(S+U1+M2)**(-1) + 8*T1*S4**(-1) + 8*S4**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(80)*N**2*CF * (  - 2*M2*TG**(-2)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*TG**(-2)*S - 2*M2*TG**(-2)*U1*
     +    S4**(-1)*MS2 - 4*M2*TG**(-2)*U1*S4**(-1)*MG2 - 4*M2*TG**(-2)*
     +    U1 - 4*M2*TG**(-2)*S4**(-1)*MS2*MG2 - 4*M2*TG**(-2)*S4**(-1)*
     +    MS2**2 - 2*M2*TG**(-2)*MS2 + 4*M2*TG**(-2)*MG2 + 4*M2*
     +    TG**(-1)*S*(S+U1+M2)**(-1) - 10*M2*TG**(-1)*S4**(-1)*MS2 + 2*
     +    M2*TG**(-1)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 4*M2*TG**(-1) - 4*M2*S4**(-1) - 2*M2*
     +    (S+U1+M2)**(-1) + 2*M2**2*TG**(-2)*U1*S4**(-1) + 4*M2**2*
     +    TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2**2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2**2*TG**(-2) - 4*M2**2*TG**(-1)*
     +    S4**(-1) + 6*TG**(-2)*S*MS2 - 2*TG**(-2)*S*MG2 - 4*TG**(-2)*
     +    U1*S4**(-1)*MS2*MG2 - 4*TG**(-2)*U1*S4**(-1)*MS2**2 + 4*
     +    TG**(-2)*MS2*MG2 + 4*TG**(-2)*MS2**2 + 4*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 - 2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 4*
     +    TG**(-1)*U1*S4**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(80)*N**2*CF * (  - 4*TG**(-1)*U1*S4**(-1)*
     +    MG2 - 4*TG**(-1)*S4**(-1)*MS2*MG2 - 4*TG**(-1)*S4**(-1)*
     +    MS2**2 + 4*TG**(-1)*MG2 + 2*S*(S+U1+M2)**(-1) - 4*T1*S4**(-1)
     +     - 2*U1*S4**(-1) - 10*S4**(-1)*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(81)*N*CF**2 * ( 4*M2*TG**(-1)*S*MS2 - 4*M2
     +    *TG**(-1)*U1*MS2 + 4*M2*S + 4*M2*T1 + 4*M2**2*TG**(-1)*S - 4*
     +    M2**2*TG**(-1)*U1 - 4*M2**2*TG**(-1)*MS2 + 4*M2**2 + 4*M2**3*
     +    TG**(-1) + 4*S*MS2 )
     +
      MQPLRH = MQPLRH + ANG4(82)*N*CF**2 * ( 4*M2*TG**(-2)*S*MS2 - 4*M2
     +    *TG**(-2)*U1*MS2 + 24*M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) + 16*M2*TG**(-1)*U1*S4**(-1)*MS2 - 24*M2*
     +    TG**(-1)*U1 + 16*M2*TG**(-1)*U1**2*S4**(-1) - 16*M2*TG**(-1)*
     +    MS2 - 8*M2*S*(S+U1+M2)**(-1) + 16*M2 + 4*M2**2*TG**(-2)*S - 4
     +    *M2**2*TG**(-2)*U1 - 4*M2**2*TG**(-2)*MS2 - 16*M2**2*TG**(-1)
     +    *S*(S+U1+M2)**(-1) + 20*M2**2*TG**(-1) + 4*M2**3*TG**(-2) + 
     +    20*TG**(-1)*S*MS2 - 16*TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 20
     +    *TG**(-1)*U1*MS2 + 16*TG**(-1)*U1**2*S4**(-1)*MS2 - 8*S*T1*
     +    (S+U1+M2)**(-1) + 12*S - 8*S**2*(S+U1+M2)**(-1) + 12*T1 + 8*
     +    U1*S4**(-1)*MS2 - 8*U1 + 8*U1**2*S4**(-1) )
     +
      MQPLRH = MQPLRH + ANG4(82)*N**2*CF * (  - 2*M2*TG**(-2)*S*MS2 + 2
     +    *M2*TG**(-2)*U1*MS2 - 10*M2*TG**(-1)*S + 8*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) - 4*M2*TG**(-1)*U1*S4**(-1)*MS2 + 6*M2*
     +    TG**(-1)*U1 - 4*M2*TG**(-1)*U1**2*S4**(-1) + 4*M2*TG**(-1)*
     +    MS2 + 4*M2*S*(S+U1+M2)**(-1) - 6*M2 - 2*M2**2*TG**(-2)*S + 2*
     +    M2**2*TG**(-2)*U1 + 2*M2**2*TG**(-2)*MS2 + 8*M2**2*TG**(-1)*S
     +    *(S+U1+M2)**(-1) - 8*M2**2*TG**(-1) - 2*M2**3*TG**(-2) - 8*
     +    TG**(-1)*S*MS2 + 8*TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 6*
     +    TG**(-1)*U1*MS2 - 4*TG**(-1)*U1**2*S4**(-1)*MS2 + 4*S*T1*
     +    (S+U1+M2)**(-1) - 4*S + 4*S**2*(S+U1+M2)**(-1) - 4*T1 - 2*U1*
     +    S4**(-1)*MS2 + 2*U1 - 2*U1**2*S4**(-1) )
     +
      MQPLRH = MQPLRH + COLO1(9)*N*CF**2*(S4+MS2) * (  - 8*TG**(-2)*S*
     +    T1**2*S4**(-2)*(S+U1)**(-2)*MS2 - 8*TG**(-2)*S*S4**(-2)*MS2
     +     + 8*TG**(-2)*T1*U1*S4**(-2) + 8*TG**(-2)*T1**3*U1*S4**(-2)*
     +    (S+U1)**(-2) + 16*S*T1**2*U1*S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)
     +     - 8*S*T1**2*U1**2*S4**(-2)*(S+T1)**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 + 16*S*T1**2*U1**3*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2) - 8*S*T1**2*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 + 8*S**2*T1*U1*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 16*S**2*T1*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 + 8*S**2*T1*U1**3*
     +    S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2) - 16*S**2*T1*
     +    S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S**3*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S**3*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 + 8*T1**3*U1*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 8*T1**3*U1**3*S4**(-2)*(S+T1)**(-2)
     +    *(M2*(S+T1)+T1*U1)**(-2) )

      END IF

      IF (IFL.EQ.1) THEN

      MQQLLH = 0.D0
      MQQLLH = MQQLLH + N*CF**2*(S4+MS2) * ( 8*M2*TG**(-2)*S4**(-1)*
     +    (S+U1)**(-1)*MG2 + 8*M2*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*
     +    MG2 + 8*M2*S*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*M2*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*M2*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*TG**(-2)*S*S4**(-1)*
     +    (S+U1)**(-1)*MG2 - 8*TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-2)*MG2
     +     + 8*TG**(-1)*S4**(-1)*(S+U1)**(-1)*MG2 + 8*UG**(-2)*S*
     +    S4**(-1)*(S+T1)**(-1)*MG2 - 8*UG**(-2)*T1*U1*S4**(-1)*
     +    (S+T1)**(-2)*MG2 + 8*UG**(-2)*U1*S4**(-1)*(S+T1)**(-1)*MG2 + 
     +    8*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*T1*S4**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + N*CF**2*(S4+MS2) * ( 8*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 - 8*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 - 8*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + N*CF**2 * ( 16*M2*TG**(-2)*S4**(-2)*MS2*MG2 + 
     +    16*M2*TG**(-2)*S4**(-1)*MG2 + 16*TG**(-2)*U1*S4**(-2)*MS2*MG2
     +     + 8*TG**(-2)*U1*S4**(-1)*MG2 - 16*TG**(-2)*S4**(-1)*MS2*MG2
     +     - 8*TG**(-2)*MG2 + 16*TG**(-1)*S4**(-2)*MS2*MG2 + 16*
     +    TG**(-1)*S4**(-1)*MG2 + 16*UG**(-2)*T1*S4**(-2)*MS2*MG2 + 8*
     +    UG**(-2)*T1*S4**(-1)*MG2 + 16*UG**(-2)*U1*S4**(-2)*MS2*MG2 + 
     +    16*UG**(-2)*U1*S4**(-1)*MG2 - 16*UG**(-2)*S4**(-1)*MS2*MG2 - 
     +    8*UG**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + CF**2*(S4+MS2) * ( 32*M2*TG**(-1)*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*M2*UG**(-1)*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 - 16*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*MG2 + 16*UG**(-1)*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 + 16*UG**(-1)*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 - 16*UG**(-1)*S4**(-1)*
     +    (S+T1)**(-1)*MG2 + 16*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + CF**2 * (  - 32*M2*TG**(-1)*UG**(-1)*S4**(-2)*
     +    MS2*MG2 - 24*M2*TG**(-1)*UG**(-1)*S4**(-1)*MG2 - 32*TG**(-1)*
     +    UG**(-1)*U1*S4**(-2)*MS2*MG2 - 24*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 + 32*TG**(-1)*UG**(-1)*S4**(-1)*MS2*MG2 + 32*
     +    TG**(-1)*UG**(-1)*MG2 - 32*UG**(-1)*S4**(-2)*MS2*MG2 - 24*
     +    UG**(-1)*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(48)*N*CF * (  - 8*TG**(-1)*UG**(-1)*S**2*
     +    MG2 - 8*S**2*(S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(48)*N*CF**2 * ( 16*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 + 16*UG**(-1)*S**2*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(48)*N**2*CF * (  - 8*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 - 8*UG**(-1)*S**2*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(48)*CF**2 * ( 8*TG**(-1)*UG**(-1)*S**2*MG2
     +     + 8*S**2*(S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(49)*N*CF * ( 2*M2*TG**(-1)*UG**(-1)*S*MG2
     +     - 2*M2*TG**(-1)*UG**(-1)*U1*MG2 - 2*M2*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 + 2*M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 - 
     +    2*M2*UG**(-1)*MG2 + 2*M2*(S+T1+M2)**(-1)*MG2 - 2*M2**2*
     +    TG**(-1)*UG**(-1)*MG2 + 2*M2**2*TG**(-1)*(S+T1+M2)**(-1)*MG2
     +     + 2*UG**(-1)*S*MG2 - 2*UG**(-1)*T1*MG2 - 2*UG**(-1)*U1*MG2
     +     - 2*S*(S+T1+M2)**(-1)*MG2 + 2*T1*(S+T1+M2)**(-1)*MG2 + 2*U1*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(49)*N*CF**2 * ( 4*M2*TG**(-2)*S*MG2 - 4*M2
     +    *TG**(-2)*U1*MG2 - 8*M2*TG**(-1)*MG2 - 4*M2**2*TG**(-2)*MG2
     +     + 4*TG**(-1)*S*MG2 - 4*TG**(-1)*U1*MG2 - 8*UG**(-1)*S*T1*
     +    (S+T1+M2)**(-1)*MG2 + 8*UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MG2 + 
     +    8*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MG2 - 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(49)*N**2*CF * (  - 2*M2*TG**(-2)*S*MG2 + 2
     +    *M2*TG**(-2)*U1*MG2 + 4*M2*TG**(-1)*MG2 + 2*M2**2*TG**(-2)*
     +    MG2 - 2*TG**(-1)*S*MG2 + 2*TG**(-1)*U1*MG2 + 2*UG**(-1)*S*T1*
     +    (S+T1+M2)**(-1)*MG2 - 2*UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MG2 - 
     +    2*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MG2 + 2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(49)*CF**2 * (  - 4*M2*TG**(-1)*UG**(-1)*S*
     +    MG2 + 4*M2*TG**(-1)*UG**(-1)*U1*MG2 + 4*M2*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 - 4*M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 + 
     +    4*M2*UG**(-1)*MG2 - 4*M2*(S+T1+M2)**(-1)*MG2 + 4*M2**2*
     +    TG**(-1)*UG**(-1)*MG2 - 4*M2**2*TG**(-1)*(S+T1+M2)**(-1)*MG2
     +     - 4*UG**(-1)*S*MG2 + 4*UG**(-1)*T1*MG2 + 4*UG**(-1)*U1*MG2
     +     + 4*S*(S+T1+M2)**(-1)*MG2 - 4*T1*(S+T1+M2)**(-1)*MG2 - 4*U1*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(50)*N*CF * (  - 2*M2*TG**(-1)*UG**(-1)*U1*
     +    MG2 + 6*M2*UG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 2*M2*UG**(-1)*T1*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*UG**(-1)*MG2 + 2*M2**2*UG**(-1)*
     +    (S+U1+M2)**(-1)*MG2 + 2*TG**(-1)*UG**(-1)*S*U1*MG2 - 2*
     +    TG**(-1)*UG**(-1)*U1**2*MG2 - 2*UG**(-1)*S*T1*(S+U1+M2)**(-1)
     +    *MG2 - 4*UG**(-1)*S*MG2 + 4*UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2
     +     + 2*UG**(-1)*T1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(50)*N*CF**2 * ( 16*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 - 8*M2*(S+U1+M2)**(-1)*MG2 - 16*TG**(-1)*
     +    S*MG2 + 16*TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 8*TG**(-1)*U1*
     +    MG2 + 4*UG**(-2)*S*U1*MG2 - 4*UG**(-2)*T1*U1*MG2 - 4*UG**(-2)
     +    *U1**2*MG2 - 8*S*(S+U1+M2)**(-1)*MG2 + 8*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(50)*N**2*CF * (  - 4*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*(S+U1+M2)**(-1)*MG2 + 4*TG**(-1)*S
     +    *MG2 - 4*TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 - 2*TG**(-1)*U1*
     +    MG2 - 2*UG**(-2)*S*U1*MG2 + 2*UG**(-2)*T1*U1*MG2 + 2*UG**(-2)
     +    *U1**2*MG2 + 2*S*(S+U1+M2)**(-1)*MG2 - 2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(50)*CF**2 * ( 4*M2*TG**(-1)*UG**(-1)*U1*
     +    MG2 - 12*M2*UG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 4*M2*UG**(-1)*T1
     +    *(S+U1+M2)**(-1)*MG2 + 4*M2*UG**(-1)*MG2 - 4*M2**2*UG**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 4*TG**(-1)*UG**(-1)*S*U1*MG2 + 4*
     +    TG**(-1)*UG**(-1)*U1**2*MG2 + 4*UG**(-1)*S*T1*(S+U1+M2)**(-1)
     +    *MG2 + 8*UG**(-1)*S*MG2 - 8*UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2
     +     - 4*UG**(-1)*T1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(63)*N*CF * (  - 2*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 + 2*TG**(-1)*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(63)*N*CF**2 * (  - 4*UG**(-2)*U1*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(63)*N**2*CF * ( 2*UG**(-2)*S*
     +    (S+T1+M2)**(-1)*MG2 + 2*UG**(-2)*T1*(S+T1+M2)**(-1)*MG2 + 2*
     +    UG**(-2)*U1*S4**(-1)*MG2 - 2*UG**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(63)*CF**2 * ( 4*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 - 4*TG**(-1)*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(68)*N*CF * (  - 2*M2*TG**(-1)*UG**(-1)*
     +    S4**(-1)*MG2 + 2*TG**(-1)*UG**(-1)*MG2 - 2*UG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(68)*N*CF**2 * (  - 4*M2*TG**(-2)*S4**(-1)*
     +    MG2 - 4*TG**(-1)*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(68)*N**2*CF * ( 2*M2*TG**(-2)*S4**(-1)*MG2
     +     - 2*M2*TG**(-2)*(S+U1+M2)**(-1)*MG2 + 2*TG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(68)*CF**2 * ( 4*M2*TG**(-1)*UG**(-1)*
     +    S4**(-1)*MG2 - 4*TG**(-1)*UG**(-1)*MG2 + 4*UG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(72)*N*CF**2 * (  - 8*S*MS2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(73)*N*CF**2 * (  - 4*M2*TG**(-1)*S*MG2 + 4
     +    *M2*TG**(-1)*U1*MG2 + 4*M2*MG2 + 4*M2**2*TG**(-1)*MG2 - 4*S*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(73)*N**2*CF * (  - 8*TG**(-1)*S*MS2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(74)*CF**2 * (  - 16*S*(S+2*M2)**(-1)*MS2*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(75)*N*CF * (  - 2*M2*TG**(-1)*UG**(-1)*S*
     +    MG2 + 2*M2*TG**(-1)*UG**(-1)*U1*MG2 - 2*M2*UG**(-1)*T1*
     +    (S+2*M2)**(-1)*MG2 - 2*M2*UG**(-1)*U1*(S+2*M2)**(-1)*MG2 + 2*
     +    M2*UG**(-1)*MG2 - 6*M2*S*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2
     +     + 2*M2*T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 2*M2*
     +    (S+2*M2)**(-1)*MG2 + 2*M2**2*TG**(-1)*UG**(-1)*MG2 - 2*M2**2*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 8*TG**(-1)*S*
     +    (S+2*M2)**(-1)*MS2*MG2 - 2*UG**(-1)*S*T1*(S+2*M2)**(-1)*MG2
     +     - 4*UG**(-1)*S*U1*(S+2*M2)**(-1)*MG2 - 8*UG**(-1)*S*
     +    (S+2*M2)**(-1)*MS2*MG2 - 2*UG**(-1)*S*MG2 - 16*UG**(-1)*T1*U1
     +    *S4**(-1)*MG2 - 16*UG**(-1)*T1*S4**(-1)*MS2*MG2 + 8*UG**(-1)*
     +    T1*MG2 - 8*UG**(-1)*T1**2*S4**(-1)*MG2 - 16*UG**(-1)*U1*
     +    S4**(-1)*MS2*MG2 + 8*UG**(-1)*U1*MG2 - 8*UG**(-1)*U1**2*
     +    S4**(-1)*MG2 + 16*UG**(-1)*MS2*MG2 + 2*S*T1*(S+U1+M2)**(-1)*
     +    (S+2*M2)**(-1)*MG2 + 4*S*(S+2*M2)**(-1)*MG2 - 4*S**2*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(75)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MG2 + 4
     +    *M2*TG**(-2)*U1*MG2 + 16*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 
     +    32*M2*TG**(-1)*U1*S4**(-1)*MG2 - 32*M2*TG**(-1)*S4**(-1)*MS2*
     +    MG2 + 16*M2*TG**(-1)*MG2 - 16*M2*S4**(-1)*MG2 - 8*M2*
     +    (S+U1+M2)**(-1)*MG2 + 4*M2**2*TG**(-2)*MG2 - 16*M2**2*
     +    TG**(-1)*S4**(-1)*MG2 - 20*TG**(-1)*S*MG2 + 16*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 - 32*TG**(-1)*U1*S4**(-1)*MS2*MG2 + 20*
     +    TG**(-1)*U1*MG2 - 16*TG**(-1)*U1**2*S4**(-1)*MG2 + 32*
     +    TG**(-1)*MS2*MG2 - 8*S*(S+U1+M2)**(-1)*MG2 - 16*T1*S4**(-1)*
     +    MG2 - 32*U1*S4**(-1)*MG2 - 32*S4**(-1)*MS2*MG2 + 20*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(75)*N**2*CF * ( 2*M2*TG**(-2)*S*MG2 - 2*M2
     +    *TG**(-2)*U1*MG2 - 4*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 16*
     +    M2*TG**(-1)*U1*S4**(-1)*MG2 + 16*M2*TG**(-1)*S4**(-1)*MS2*MG2
     +     - 8*M2*TG**(-1)*MG2 + 8*M2*S4**(-1)*MG2 + 2*M2*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2**2*TG**(-2)*MG2 + 8*M2**2*TG**(-1)
     +    *S4**(-1)*MG2 + 6*TG**(-1)*S*MG2 - 4*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1)*MG2 + 16*TG**(-1)*U1*S4**(-1)*MS2*MG2 - 8*
     +    TG**(-1)*U1*MG2 + 8*TG**(-1)*U1**2*S4**(-1)*MG2 - 16*TG**(-1)
     +    *MS2*MG2 + 2*S*(S+U1+M2)**(-1)*MG2 + 8*T1*S4**(-1)*MG2 + 16*
     +    U1*S4**(-1)*MG2 + 16*S4**(-1)*MS2*MG2 - 8*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(75)*CF**2 * ( 4*M2*TG**(-1)*UG**(-1)*S*MG2
     +     - 4*M2*TG**(-1)*UG**(-1)*U1*MG2 - 4*M2*TG**(-1)*S*
     +    (S+2*M2)**(-1)*MG2 + 4*M2*TG**(-1)*U1*(S+2*M2)**(-1)*MG2 - 4*
     +    M2*UG**(-1)*MG2 + 12*M2*S*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2
     +     - 4*M2*T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 4*M2**2*
     +    TG**(-1)*UG**(-1)*MG2 + 4*M2**2*TG**(-1)*(S+2*M2)**(-1)*MG2
     +     + 4*M2**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 4*UG**(-1)*S*
     +    MG2 + 16*UG**(-1)*T1*U1*S4**(-1)*MG2 + 16*UG**(-1)*T1*
     +    S4**(-1)*MS2*MG2 - 8*UG**(-1)*T1*MG2 + 8*UG**(-1)*T1**2*
     +    S4**(-1)*MG2 + 16*UG**(-1)*U1*S4**(-1)*MS2*MG2 - 8*UG**(-1)*
     +    U1*MG2 + 8*UG**(-1)*U1**2*S4**(-1)*MG2 - 16*UG**(-1)*MS2*MG2
     +     - 4*S*T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 16*S*
     +    (S+2*M2)**(-2)*MS2*MG2 - 12*S*(S+2*M2)**(-1)*MG2 + 8*S**2*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(77)*N**2*CF * ( 4*M2*TG**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(78)*N*CF * (  - 2*M2*TG**(-1)*UG**(-1)*
     +    S4**(-1)*MG2 + 2*TG**(-1)*UG**(-1)*MG2 - 2*UG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(78)*N**2*CF * ( 2*M2*TG**(-2)*S4**(-1)*MG2
     +     - 2*M2*TG**(-2)*(S+U1+M2)**(-1)*MG2 + 2*TG**(-1)*S4**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(79)*N**2*CF * ( 4*M2*TG**(-2)*U1*MG2 + 4*
     +    M2*TG**(-1)*MG2 + 8*M2**2*TG**(-2)*MG2 - 8*TG**(-2)*S*MS2*MG2
     +     - 4*TG**(-1)*S*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(80)*N*CF * (  - 2*M2*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 - 8*M2*TG**(-1)*UG**(-1)*S4**(-1)*MS2*MG2 + 6*M2
     +    *TG**(-1)*UG**(-1)*MG2 + 2*M2*TG**(-1)*(S+U1+M2)**(-1)*MG2 - 
     +    4*M2*UG**(-1)*S4**(-1)*MG2 - 4*M2**2*TG**(-1)*UG**(-1)*
     +    S4**(-1)*MG2 + 4*TG**(-1)*UG**(-1)*S*MG2 - 8*TG**(-1)*
     +    UG**(-1)*U1*S4**(-1)*MS2*MG2 + 2*TG**(-1)*UG**(-1)*U1*MG2 + 8
     +    *TG**(-1)*UG**(-1)*MS2*MG2 + 2*TG**(-1)*MG2 - 6*UG**(-1)*T1*
     +    S4**(-1)*MG2 - 4*UG**(-1)*U1*S4**(-1)*MG2 - 8*UG**(-1)*
     +    S4**(-1)*MS2*MG2 + 6*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(80)*N*CF**2 * (  - 8*M2*TG**(-1)*S4**(-1)*
     +    MG2 - 8*M2*TG**(-1)*(S+U1+M2)**(-1)*MG2 - 8*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(80)*N**2*CF * ( 2*M2*TG**(-2)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*TG**(-2)*U1*S4**(-1)*MG2 + 8*M2*
     +    TG**(-2)*S4**(-1)*MS2*MG2 + 2*M2*TG**(-2)*MG2 + 10*M2*
     +    TG**(-1)*S4**(-1)*MG2 + 4*M2**2*TG**(-2)*S4**(-1)*MG2 - 2*
     +    M2**2*TG**(-2)*(S+U1+M2)**(-1)*MG2 - 4*TG**(-2)*S*MG2 + 8*
     +    TG**(-2)*U1*S4**(-1)*MS2*MG2 - 8*TG**(-2)*MS2*MG2 - 2*
     +    TG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 4*TG**(-1)*U1*S4**(-1)*MG2
     +     + 8*TG**(-1)*S4**(-1)*MS2*MG2 + 6*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(80)*CF**2 * ( 4*UG**(-1)*T1*S4**(-1)*MG2
     +     - 4*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(81)*N*CF**2 * (  - 4*M2*TG**(-1)*S*MG2 + 4
     +    *M2*TG**(-1)*U1*MG2 + 4*M2**2*TG**(-1)*MG2 - 4*S*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(82)*N*CF * ( 2*M2*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 - 2*M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 - 
     +    2*M2*UG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 2*M2*UG**(-1)*T1*
     +    (S+T1+M2)**(-1)*MG2 - 2*M2*UG**(-1)*U1*S4**(-1)*MG2 + 2*M2*
     +    UG**(-1)*MG2 - 2*M2*(S+T1+M2)**(-1)*MG2 - 2*M2**2*TG**(-1)*
     +    (S+T1+M2)**(-1)*MG2 + 4*UG**(-1)*S*T1*(S+T1+M2)**(-1)*MG2 - 4
     +    *UG**(-1)*S*MG2 + 4*UG**(-1)*S**2*(S+T1+M2)**(-1)*MG2 - 2*
     +    UG**(-1)*T1*U1*S4**(-1)*MG2 - 2*UG**(-1)*T1*U1*
     +    (S+T1+M2)**(-1)*MG2 + 4*UG**(-1)*U1*MG2 - 4*UG**(-1)*U1**2*
     +    S4**(-1)*MG2 + 4*S*(S+T1+M2)**(-1)*MG2 - 8*S**2*
     +    (S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2 - 2*U1*(S+T1+M2)**(-1)*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(82)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MG2 + 4
     +    *M2*TG**(-2)*U1*MG2 - 16*M2*TG**(-1)*U1*S4**(-1)*MG2 + 16*M2*
     +    TG**(-1)*MG2 + 4*M2**2*TG**(-2)*MG2 - 20*TG**(-1)*S*MG2 + 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 20*TG**(-1)*U1*MG2 - 16*
     +    TG**(-1)*U1**2*S4**(-1)*MG2 - 8*U1*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(82)*N**2*CF * ( 2*M2*TG**(-2)*S*MG2 - 2*M2
     +    *TG**(-2)*U1*MG2 + 4*M2*TG**(-1)*U1*S4**(-1)*MG2 - 4*M2*
     +    TG**(-1)*MG2 - 2*M2**2*TG**(-2)*MG2 + 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 - 6*TG**(-1)*U1*MG2 + 4*
     +    TG**(-1)*U1**2*S4**(-1)*MG2 + 2*U1*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(82)*CF**2 * ( 4*M2*TG**(-1)*UG**(-1)*S*MG2
     +     - 4*M2*TG**(-1)*UG**(-1)*U1*MG2 - 4*M2*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 + 4*M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 + 
     +    4*M2*UG**(-1)*U1*S4**(-1)*MG2 + 4*M2*(S+T1+M2)**(-1)*MG2 - 4*
     +    M2**2*TG**(-1)*UG**(-1)*MG2 + 4*M2**2*TG**(-1)*
     +    (S+T1+M2)**(-1)*MG2 + 8*UG**(-1)*S*MG2 + 4*UG**(-1)*T1*U1*
     +    S4**(-1)*MG2 - 4*UG**(-1)*U1*MG2 + 8*UG**(-1)*U1**2*S4**(-1)*
     +    MG2 - 8*S*(S+T1+M2)**(-1)*MG2 + 8*S**2*(S+U1+M2)**(-1)*
     +    (S+T1+M2)**(-1)*MG2 + 4*U1*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(85)*N*CF**2 * (  - 4*M2*MG2 - 4*UG**(-1)*S
     +    *U1*MG2 + 4*UG**(-1)*T1*U1*MG2 + 4*UG**(-1)*U1**2*MG2 - 4*T1*
     +    MG2 - 4*U1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(85)*N**2*CF * (  - 2*M2*UG**(-1)*S*MG2 + 2
     +    *M2*UG**(-1)*U1*MG2 + 2*UG**(-1)*S*U1*MG2 - 2*UG**(-1)*U1**2*
     +    MG2 - 2*S*MG2 + 2*U1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(86)*N*CF * ( 2*M2*TG**(-1)*UG**(-1)*S*MG2
     +     - 2*M2*TG**(-1)*UG**(-1)*U1*MG2 - 2*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*TG**(-1)*U1*S4**(-1)*MG2 + 2*M2*
     +    TG**(-1)*MG2 - 6*M2*UG**(-1)*S*(S+U1+M2)**(-1)*MG2 + 2*M2*
     +    UG**(-1)*T1*(S+U1+M2)**(-1)*MG2 + 2*M2*UG**(-1)*MG2 - 6*M2*
     +    S4**(-1)*MG2 - 6*M2**2*TG**(-1)*S4**(-1)*MG2 + 4*M2**2*
     +    TG**(-1)*(S+U1+M2)**(-1)*MG2 - 2*M2**2*UG**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 2*TG**(-1)*UG**(-1)*S*U1*MG2 + 2*
     +    TG**(-1)*UG**(-1)*U1**2*MG2 + 2*TG**(-1)*S*MG2 - 2*TG**(-1)*
     +    U1*MG2 + 2*UG**(-1)*S*T1*(S+U1+M2)**(-1)*MG2 + 4*UG**(-1)*S*
     +    MG2 - 4*UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 - 2*UG**(-1)*T1*MG2
     +     - 2*UG**(-1)*U1*MG2 + 2*S*(S+U1+M2)**(-1)*MG2 - 8*S**2*
     +    (S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2 - 4*T1*S4**(-1)*MG2 - 2*
     +    U1*S4**(-1)*MG2 + 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(86)*N*CF**2 * (  - 8*M2*UG**(-1)*T1*
     +    S4**(-1)*MG2 + 8*M2*UG**(-1)*MG2 - 4*UG**(-2)*S*U1*MG2 + 4*
     +    UG**(-2)*T1*U1*MG2 + 4*UG**(-2)*U1**2*MG2 - 16*UG**(-1)*S*MG2
     +     + 16*UG**(-1)*S**2*(S+T1+M2)**(-1)*MG2 - 8*UG**(-1)*T1*U1*
     +    S4**(-1)*MG2 + 16*UG**(-1)*T1*MG2 - 16*UG**(-1)*T1**2*
     +    S4**(-1)*MG2 - 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(86)*N**2*CF * ( 2*M2*UG**(-1)*T1*S4**(-1)*
     +    MG2 - 2*M2*UG**(-1)*MG2 + 2*UG**(-2)*S*U1*MG2 - 2*UG**(-2)*T1
     +    *U1*MG2 - 2*UG**(-2)*U1**2*MG2 + 6*UG**(-1)*S*MG2 - 8*
     +    UG**(-1)*S**2*(S+T1+M2)**(-1)*MG2 + 2*UG**(-1)*T1*U1*S4**(-1)
     +    *MG2 - 4*UG**(-1)*T1*MG2 + 4*UG**(-1)*T1**2*S4**(-1)*MG2 + 2*
     +    UG**(-1)*U1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(86)*CF**2 * (  - 4*M2*TG**(-1)*UG**(-1)*U1
     +    *MG2 + 4*M2*TG**(-1)*U1*S4**(-1)*MG2 + 4*M2*TG**(-1)*MG2 + 12
     +    *M2*UG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 4*M2*UG**(-1)*T1*
     +    (S+U1+M2)**(-1)*MG2 - 4*M2*UG**(-1)*MG2 + 12*M2*S4**(-1)*MG2
     +     + 4*M2*(S+U1+M2)**(-1)*MG2 + 12*M2**2*TG**(-1)*S4**(-1)*MG2
     +     + 4*M2**2*UG**(-1)*(S+U1+M2)**(-1)*MG2 + 4*TG**(-1)*UG**(-1)
     +    *S*U1*MG2 - 4*TG**(-1)*UG**(-1)*U1**2*MG2 + 4*TG**(-1)*S*MG2
     +     + 4*TG**(-1)*U1*MG2 - 4*UG**(-1)*S*T1*(S+U1+M2)**(-1)*MG2 - 
     +    8*UG**(-1)*S*MG2 + 8*UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 4*
     +    UG**(-1)*T1*MG2 + 8*S**2*(S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2
     +     + 8*T1*S4**(-1)*MG2 + 4*U1*S4**(-1)*MG2 - 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(87)*N*CF**2 * (  - 8*S*MS2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(88)*N*CF**2 * (  - 4*UG**(-1)*S*U1*MG2 + 4
     +    *UG**(-1)*T1*U1*MG2 + 4*UG**(-1)*U1**2*MG2 - 4*T1*MG2 - 4*U1*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(88)*N**2*CF * ( 2*M2*UG**(-1)*T1*MG2 + 2*
     +    M2*UG**(-1)*U1*MG2 - 8*UG**(-1)*S*MS2*MG2 - 2*UG**(-1)*T1*U1*
     +    MG2 - 2*UG**(-1)*U1**2*MG2 + 2*T1*MG2 + 2*U1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(89)*CF**2 * (  - 16*S*(S+2*M2)**(-1)*MS2*
     +    MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*N*CF * ( 2*M2*TG**(-1)*UG**(-1)*U1*MG2
     +     - 2*M2*TG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 16*M2*TG**(-1)*U1*
     +    S4**(-1)*MG2 + 2*M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 - 16*M2*
     +    TG**(-1)*S4**(-1)*MS2*MG2 + 4*M2*TG**(-1)*MG2 - 6*M2*UG**(-1)
     +    *S*(S+U1+M2)**(-1)*MG2 + 2*M2*UG**(-1)*T1*(S+U1+M2)**(-1)*MG2
     +     - 2*M2*UG**(-1)*T1*(S+2*M2)**(-1)*MG2 - 2*M2*UG**(-1)*U1*
     +    (S+2*M2)**(-1)*MG2 + 2*M2*UG**(-1)*MG2 - 6*M2*S*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 2*M2*T1*(S+U1+M2)**(-1)*
     +    (S+2*M2)**(-1)*MG2 - 8*M2*S4**(-1)*MG2 + 2*M2*(S+T1+M2)**(-1)
     +    *MG2 + 2*M2*(S+2*M2)**(-1)*MG2 - 8*M2**2*TG**(-1)*S4**(-1)*
     +    MG2 + 2*M2**2*TG**(-1)*(S+T1+M2)**(-1)*MG2 - 2*M2**2*UG**(-1)
     +    *(S+U1+M2)**(-1)*MG2 - 2*M2**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)
     +    *MG2 - 2*TG**(-1)*UG**(-1)*S*U1*MG2 + 2*TG**(-1)*UG**(-1)*
     +    U1**2*MG2 - 8*TG**(-1)*S*(S+2*M2)**(-1)*MS2*MG2 - 16*TG**(-1)
     +    *U1*S4**(-1)*MS2*MG2 + 4*TG**(-1)*U1*MG2 - 8*TG**(-1)*U1**2*
     +    S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*N*CF * ( 16*TG**(-1)*MS2*MG2 + 2*
     +    UG**(-1)*S*T1*(S+U1+M2)**(-1)*MG2 - 2*UG**(-1)*S*T1*
     +    (S+2*M2)**(-1)*MG2 - 4*UG**(-1)*S*U1*(S+2*M2)**(-1)*MG2 - 8*
     +    UG**(-1)*S*(S+2*M2)**(-1)*MS2*MG2 + 4*UG**(-1)*S*MG2 - 4*
     +    UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 2*UG**(-1)*U1*MG2 + 2*S*
     +    T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 2*S*(S+T1+M2)**(-1)*
     +    MG2 + 4*S*(S+2*M2)**(-1)*MG2 - 4*S**2*(S+U1+M2)**(-1)*
     +    (S+2*M2)**(-1)*MG2 - 8*T1*S4**(-1)*MG2 + 2*T1*(S+T1+M2)**(-1)
     +    *MG2 - 16*U1*S4**(-1)*MG2 + 2*U1*(S+T1+M2)**(-1)*MG2 - 16*
     +    S4**(-1)*MS2*MG2 + 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*N*CF**2 * (  - 4*UG**(-2)*S*U1*MG2 + 4
     +    *UG**(-2)*T1*U1*MG2 + 4*UG**(-2)*U1**2*MG2 - 8*UG**(-1)*S*T1*
     +    (S+T1+M2)**(-1)*MG2 - 32*UG**(-1)*T1*U1*S4**(-1)*MG2 + 8*
     +    UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MG2 - 32*UG**(-1)*T1*S4**(-1)*
     +    MS2*MG2 + 8*UG**(-1)*T1*MG2 - 16*UG**(-1)*T1**2*S4**(-1)*MG2
     +     + 8*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MG2 - 32*UG**(-1)*U1*
     +    S4**(-1)*MS2*MG2 + 8*UG**(-1)*U1*MG2 - 16*UG**(-1)*U1**2*
     +    S4**(-1)*MG2 + 32*UG**(-1)*MS2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*N**2*CF * ( 2*UG**(-2)*S*U1*MG2 - 2*
     +    UG**(-2)*T1*U1*MG2 - 2*UG**(-2)*U1**2*MG2 + 2*UG**(-1)*S*T1*
     +    (S+T1+M2)**(-1)*MG2 + 16*UG**(-1)*T1*U1*S4**(-1)*MG2 - 2*
     +    UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MG2 + 16*UG**(-1)*T1*S4**(-1)*
     +    MS2*MG2 - 4*UG**(-1)*T1*MG2 + 8*UG**(-1)*T1**2*S4**(-1)*MG2
     +     - 2*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MG2 + 16*UG**(-1)*U1*
     +    S4**(-1)*MS2*MG2 - 4*UG**(-1)*U1*MG2 + 8*UG**(-1)*U1**2*
     +    S4**(-1)*MG2 - 16*UG**(-1)*MS2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*CF**2 * (  - 4*M2*TG**(-1)*UG**(-1)*U1
     +    *MG2 + 4*M2*TG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 4*M2*TG**(-1)*S*
     +    (S+2*M2)**(-1)*MG2 + 16*M2*TG**(-1)*U1*S4**(-1)*MG2 - 4*M2*
     +    TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 + 4*M2*TG**(-1)*U1*
     +    (S+2*M2)**(-1)*MG2 + 16*M2*TG**(-1)*S4**(-1)*MS2*MG2 - 4*M2*
     +    TG**(-1)*MG2 + 12*M2*UG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 4*M2*
     +    UG**(-1)*T1*(S+U1+M2)**(-1)*MG2 - 4*M2*UG**(-1)*MG2 + 12*M2*S
     +    *(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 4*M2*T1*(S+U1+M2)**(-1)
     +    *(S+2*M2)**(-1)*MG2 + 8*M2*S4**(-1)*MG2 - 4*M2*
     +    (S+T1+M2)**(-1)*MG2 + 8*M2**2*TG**(-1)*S4**(-1)*MG2 - 4*M2**2
     +    *TG**(-1)*(S+T1+M2)**(-1)*MG2 + 4*M2**2*TG**(-1)*
     +    (S+2*M2)**(-1)*MG2 + 4*M2**2*UG**(-1)*(S+U1+M2)**(-1)*MG2 + 4
     +    *M2**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 4*TG**(-1)*
     +    UG**(-1)*S*U1*MG2 - 4*TG**(-1)*UG**(-1)*U1**2*MG2 + 16*
     +    TG**(-1)*U1*S4**(-1)*MS2*MG2 - 4*TG**(-1)*U1*MG2 + 8*TG**(-1)
     +    *U1**2*S4**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(90)*CF**2 * (  - 16*TG**(-1)*MS2*MG2 - 4*
     +    UG**(-1)*S*T1*(S+U1+M2)**(-1)*MG2 - 8*UG**(-1)*S*MG2 + 8*
     +    UG**(-1)*S**2*(S+U1+M2)**(-1)*MG2 + 4*UG**(-1)*T1*MG2 - 4*S*
     +    T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 4*S*(S+T1+M2)**(-1)*
     +    MG2 - 16*S*(S+2*M2)**(-2)*MS2*MG2 - 12*S*(S+2*M2)**(-1)*MG2
     +     + 8*S**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 8*T1*S4**(-1)*
     +    MG2 - 4*T1*(S+T1+M2)**(-1)*MG2 + 16*U1*S4**(-1)*MG2 - 4*U1*
     +    (S+T1+M2)**(-1)*MG2 + 16*S4**(-1)*MS2*MG2 - 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(91)*N**2*CF * ( 8*M2*UG**(-2)*U1*MG2 - 4*
     +    M2*UG**(-1)*MG2 - 8*UG**(-2)*S*MS2*MG2 + 4*UG**(-2)*T1*U1*MG2
     +     - 4*UG**(-1)*S*MG2 - 4*UG**(-1)*T1*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(92)*N*CF * (  - 4*M2*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 - 8*M2*TG**(-1)*UG**(-1)*S4**(-1)*MS2*MG2 + 4*M2
     +    *TG**(-1)*UG**(-1)*MG2 + 2*M2*UG**(-1)*S4**(-1)*MG2 + 2*M2**2
     +    *TG**(-1)*UG**(-1)*S4**(-1)*MG2 + 4*TG**(-1)*UG**(-1)*S*MG2
     +     - 8*TG**(-1)*UG**(-1)*U1*S4**(-1)*MS2*MG2 + 4*TG**(-1)*
     +    UG**(-1)*U1*MG2 - 4*TG**(-1)*UG**(-1)*U1**2*S4**(-1)*MG2 + 8*
     +    TG**(-1)*UG**(-1)*MS2*MG2 - 2*TG**(-1)*U1*S4**(-1)*MG2 + 2*
     +    TG**(-1)*MG2 - 2*UG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 2*UG**(-1)*
     +    T1*(S+T1+M2)**(-1)*MG2 - 4*UG**(-1)*U1*S4**(-1)*MG2 - 8*
     +    UG**(-1)*S4**(-1)*MS2*MG2 + 6*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(92)*N*CF**2 * ( 8*UG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 + 8*UG**(-1)*T1*(S+T1+M2)**(-1)*MG2 - 8*
     +    UG**(-1)*U1*S4**(-1)*MG2 - 8*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(92)*N**2*CF * (  - 2*M2*UG**(-2)*T1*
     +    S4**(-1)*MG2 + 2*M2*UG**(-2)*MG2 - 4*UG**(-2)*S*T1*
     +    (S+T1+M2)**(-1)*MG2 - 4*UG**(-2)*S**2*(S+T1+M2)**(-1)*MG2 + 4
     +    *UG**(-2)*T1*U1*S4**(-1)*MG2 + 2*UG**(-2)*T1*U1*
     +    (S+T1+M2)**(-1)*MG2 + 8*UG**(-2)*T1*S4**(-1)*MS2*MG2 + 8*
     +    UG**(-2)*U1*S4**(-1)*MS2*MG2 - 2*UG**(-2)*U1*MG2 + 4*UG**(-2)
     +    *U1**2*S4**(-1)*MG2 - 8*UG**(-2)*MS2*MG2 - 2*UG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 - 2*UG**(-1)*T1*(S+T1+M2)**(-1)*MG2 + 2*
     +    UG**(-1)*U1*S4**(-1)*MG2 + 2*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(92)*CF**2 * ( 4*TG**(-1)*U1*S4**(-1)*MG2
     +     - 4*TG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(93)*N*CF * (  - 2*TG**(-1)*UG**(-1)*U1*
     +    S4**(-1)*MG2 + 2*TG**(-1)*UG**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(93)*N**2*CF * ( 2*UG**(-2)*S*
     +    (S+T1+M2)**(-1)*MG2 + 2*UG**(-2)*T1*(S+T1+M2)**(-1)*MG2 + 2*
     +    UG**(-2)*U1*S4**(-1)*MG2 - 2*UG**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(94)*N*CF * (  - 2*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*TG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 2
     +    *M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 + 6*M2*TG**(-1)*MG2 - 2*
     +    M2*UG**(-1)*S*(S+T1+M2)**(-1)*MG2 - 2*M2*UG**(-1)*T1*
     +    (S+T1+M2)**(-1)*MG2 + 2*M2*UG**(-1)*T1*(S+2*M2)**(-1)*MG2 + 2
     +    *M2*UG**(-1)*U1*(S+2*M2)**(-1)*MG2 + 4*M2*UG**(-1)*MG2 + 6*M2
     +    *S*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 2*M2*T1*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 2*M2*(S+T1+M2)**(-1)*MG2
     +     - 2*M2*(S+2*M2)**(-1)*MG2 + 4*M2**2*TG**(-1)*(S+U1+M2)**(-1)
     +    *MG2 - 2*M2**2*TG**(-1)*(S+T1+M2)**(-1)*MG2 + 2*M2**2*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 8*TG**(-1)*S*
     +    (S+2*M2)**(-1)*MS2*MG2 + 4*TG**(-1)*S*MG2 + 2*TG**(-1)*U1*MG2
     +     + 4*UG**(-1)*S*T1*(S+T1+M2)**(-1)*MG2 + 2*UG**(-1)*S*T1*
     +    (S+2*M2)**(-1)*MG2 + 4*UG**(-1)*S*U1*(S+2*M2)**(-1)*MG2 + 8*
     +    UG**(-1)*S*(S+2*M2)**(-1)*MS2*MG2 + 4*UG**(-1)*S**2*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(94)*N*CF * (  - 2*UG**(-1)*T1*U1*
     +    (S+T1+M2)**(-1)*MG2 + 4*UG**(-1)*U1*MG2 - 2*S*T1*
     +    (S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 2*S*(S+U1+M2)**(-1)*MG2
     +     + 4*S*(S+T1+M2)**(-1)*MG2 - 4*S*(S+2*M2)**(-1)*MG2 - 8*S**2*
     +    (S+U1+M2)**(-1)*(S+T1+M2)**(-1)*MG2 + 4*S**2*(S+U1+M2)**(-1)*
     +    (S+2*M2)**(-1)*MG2 - 2*U1*(S+T1+M2)**(-1)*MG2 - 2*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(94)*CF**2 * (  - 4*M2*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 + 4*M2*TG**(-1)*S*(S+2*M2)**(-1)*MG2 + 4*
     +    M2*TG**(-1)*U1*(S+T1+M2)**(-1)*MG2 - 4*M2*TG**(-1)*U1*
     +    (S+2*M2)**(-1)*MG2 - 12*M2*S*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*
     +    MG2 + 4*M2*T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 4*M2*
     +    (S+U1+M2)**(-1)*MG2 + 4*M2*(S+T1+M2)**(-1)*MG2 + 4*M2**2*
     +    TG**(-1)*(S+T1+M2)**(-1)*MG2 - 4*M2**2*TG**(-1)*
     +    (S+2*M2)**(-1)*MG2 - 4*M2**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*
     +    MG2 + 4*S*T1*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 - 8*S*
     +    (S+T1+M2)**(-1)*MG2 + 16*S*(S+2*M2)**(-2)*MS2*MG2 + 12*S*
     +    (S+2*M2)**(-1)*MG2 + 8*S**2*(S+U1+M2)**(-1)*(S+T1+M2)**(-1)*
     +    MG2 - 8*S**2*(S+U1+M2)**(-1)*(S+2*M2)**(-1)*MG2 + 4*U1*
     +    (S+T1+M2)**(-1)*MG2 + 4*MG2 )
     +
      MQQLLH = MQQLLH + ANG4(95)*N**2*CF * ( 4*M2*UG**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + COLO1(9)*N*CF**2*(S4+MS2) * ( 8*TG**(-2)*S*
     +    T1**2*S4**(-2)*(S+U1)**(-2)*MG2 + 8*TG**(-2)*S*S4**(-2)*MG2
     +     + 8*UG**(-2)*S*U1**2*S4**(-2)*(S+T1)**(-2)*MG2 + 8*UG**(-2)*
     +    S*S4**(-2)*MG2 + 8*S*T1**2*U1**2*S4**(-2)*(S+T1)**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*T1**2*U1**2*S4**(-2)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*T1**2*S4**(-2)
     +    *(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*S4**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*T1*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*S**2*T1*
     +    S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*S**2*T1**2*U1*
     +    S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*
     +    U1*S4**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**3*T1**2*
     +    S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**3*
     +    U1**2*S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*
     +    S**3*S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 8*S**3*S4**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 )
     +
      MQQLLH = MQQLLH + COLO1(9)*CF**2*(S4+MS2) * ( 16*TG**(-1)*S*T1**2
     +    *U1*S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*
     +    TG**(-1)*S*U1*S4**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*
     +    TG**(-1)*S**2*T1**2*S4**(-2)*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S**2*S4**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*UG**(-1)*S*T1*U1**2*S4**(-2)
     +    *(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-1)*MG2 + 16*UG**(-1)*S*T1*
     +    S4**(-2)*(M2*(S+T1)+T1*U1)**(-1)*MG2 + 16*UG**(-1)*S**2*U1**2
     +    *S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-1)*MG2 + 16*
     +    UG**(-1)*S**2*S4**(-2)*(M2*(S+T1)+T1*U1)**(-1)*MG2 )


      MQQLRH = 0.D0
      MQQLRH = MQQLRH + N*CF**2*(S4+MS2) * ( 8*M2*TG**(-2)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2) + 8*M2*TG**(-2)*U1*S4**(-1)*
     +    (S+U1)**(-1) - 8*M2*TG**(-2)*S4**(-1)*(S+U1)**(-1)*MS2 - 32*
     +    M2*TG**(-1)*S4**(-2) + 16*M2*TG**(-1)*S4**(-1)*(S+T1)**(-1)
     +     - 16*M2*UG**(-1)*S*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) - 16*M2*
     +    UG**(-1)*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) - 16*M2*S*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) - 8*M2*S*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 8*M2*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 16*M2*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 16*M2*S*S4**(-1)*(S+U1+M2)**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) - 8*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*M2*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*M2*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*M2*S**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 8*M2*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) )
     +
      MQQLRH = MQQLRH + N*CF**2*(S4+MS2) * (  - 8*M2*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 8*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*M2*T1**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*M2*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*M2*U1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 32*M2*S4**(-2)*(S+U1+M2)**(-1) - 16
     +    *M2*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*MS2 + 8*TG**(-2)*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*MS2 - 32*TG**(-1)*S*S4**(-2) + 16*TG**(-1)*S*
     +    S4**(-1)*(S+T1)**(-1) + 16*TG**(-1)*S*S4**(-1)*(S+U1)**(-1)
     +     + 8*TG**(-1)*T1*U1*S4**(-1)*(S+U1)**(-2) - 32*TG**(-1)*U1*
     +    S4**(-2) + 24*TG**(-1)*U1*S4**(-1)*(S+U1)**(-1) - 8*TG**(-1)*
     +    S4**(-1)*(S+U1)**(-1)*MS2 - 8*UG**(-2)*S*S4**(-1)*
     +    (S+T1)**(-1)*MS2 + 8*UG**(-2)*T1*U1*S4**(-1)*(S+T1)**(-2)*MS2
     +     + 8*UG**(-2)*T1*U1*S4**(-1)*(S+T1)**(-1) + 8*UG**(-2)*T1*
     +    U1**2*S4**(-1)*(S+T1)**(-2) )
     +
      MQQLRH = MQQLRH + N*CF**2*(S4+MS2) * (  - 8*UG**(-2)*U1*S4**(-1)*
     +    (S+T1)**(-1)*MS2 - 32*UG**(-1)*S*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 64*UG**(-1)*S*T1*
     +    S4**(-2)*(S+T1+M2)**(-1) + 32*UG**(-1)*S*T1*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) + 16*UG**(-1)*S*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 16*UG**(-1)*S*T1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 64*UG**(-1)*S*U1*
     +    S4**(-2)*(S+T1+M2)**(-1) + 16*UG**(-1)*S*U1*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) + 32*UG**(-1)*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 16*UG**(-1)*S*U1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) + 16*UG**(-1)*S*
     +    S4**(-1)*(S+U1)**(-1) - 32*UG**(-1)*S**2*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 32*UG**(-1)*S**2*U1
     +    *S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 32*
     +    UG**(-1)*S**2*S4**(-2)*(S+T1+M2)**(-1) + 16*UG**(-1)*S**2*
     +    S4**(-1)*(S+T1)**(-1)*(S+T1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + N*CF**2*(S4+MS2) * ( 16*UG**(-1)*S**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 16*UG**(-1)*S**3*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 64*UG**(-1)*T1*U1*
     +    S4**(-2)*(S+T1+M2)**(-1) + 16*UG**(-1)*T1*U1*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) - 32*UG**(-1)*T1**2*S4**(-2)*
     +    (S+T1+M2)**(-1) + 16*UG**(-1)*T1**2*S4**(-1)*(S+T1)**(-1)*
     +    (S+T1+M2)**(-1) + 16*UG**(-1)*U1*S4**(-1)*(S+U1)**(-1) - 32*
     +    UG**(-1)*U1**2*S4**(-2)*(S+T1+M2)**(-1) + 16*UG**(-1)*U1**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) - 16*S*T1*S4**(-1)*
     +    (S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) - 8*S*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) + 8*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) )
     +
      MQQLRH = MQQLRH + N*CF**2*(S4+MS2) * (  - 8*S**2*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 32*T1*S4**(-2)*(S+U1+M2)**(-1)
     +     + 24*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 8*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) + 8*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) + 8*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 64*S4**(-2) + 16*S4**(-1)*
     +    (S+T1)**(-1) + 16*S4**(-1)*(S+U1)**(-1) + 8*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MS2 + 8*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + N*CF**2 * (  - 16*M2*TG**(-2)*U1*S4**(-2)*MS2
     +     - 8*M2*TG**(-2)*U1*S4**(-1) - 16*M2*TG**(-2)*S4**(-2)*MS2**2
     +     - 16*TG**(-2)*U1*S4**(-2)*MS2**2 - 8*TG**(-2)*U1*S4**(-1)*
     +    MS2 + 16*TG**(-2)*S4**(-1)*MS2**2 + 8*TG**(-2)*MS2 - 16*
     +    TG**(-1)*U1*S4**(-2)*MS2 - 8*TG**(-1)*U1*S4**(-1) - 16*
     +    TG**(-1)*S4**(-2)*MS2**2 - 16*UG**(-2)*T1*U1*S4**(-2)*MS2 - 8
     +    *UG**(-2)*T1*U1*S4**(-1) - 16*UG**(-2)*T1*S4**(-2)*MS2**2 - 8
     +    *UG**(-2)*T1*S4**(-1)*MS2 - 16*UG**(-2)*U1*S4**(-2)*MS2**2 + 
     +    16*UG**(-2)*S4**(-1)*MS2**2 + 8*UG**(-2)*MS2 )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 4*M2*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2*TG**(-2)*S*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 4*M2*TG**(-2)*T1*
     +    U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1)*MS2 - 4*M2*TG**(-2)*
     +    T1*U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1)*MG2 - 4*M2*
     +    TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-2) + 4*M2*TG**(-2)*S4**(-1)
     +    *(S+U1)**(-1)*MS2 - 4*M2*TG**(-2)*S4**(-1)*(S+U1)**(-1)*MG2
     +     + 4*M2*TG**(-1)*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 
     +    4*M2*TG**(-1)*S*T1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*
     +    M2*TG**(-1)*S*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*
     +    TG**(-1)*S*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 - 8*M2*
     +    TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*M2*
     +    TG**(-1)*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*M2*TG**(-1)*
     +    S**2*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*TG**(-1)*
     +    S**2*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2*TG**(-1)*T1*
     +    U1*S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * ( 4*M2*TG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*TG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 - 4*M2*TG**(-1)*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 16*M2*TG**(-1)*S4**(-2) - 
     +    4*M2*TG**(-1)*S4**(-1)*(S+T1)**(-1) - 4*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 8*M2*TG**(-1)*S4**(-1)*
     +    (S+U1)**(-1) - 4*M2*TG**(-1)*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)
     +    *MS2 + 4*M2*TG**(-1)*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MG2 + 4
     +    *M2*UG**(-1)*S*T1*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 8*M2*
     +    UG**(-1)*S*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*M2*
     +    UG**(-1)*S*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 4*M2*
     +    UG**(-1)*S**2*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 4*M2*
     +    UG**(-1)*S**2*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 4*M2*
     +    UG**(-1)*T1*U1**2*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 4*M2*
     +    UG**(-1)*U1**2*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2 )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 4*M2*UG**(-1)*U1**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 4*M2*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*M2*S*S4**(-1)*(S+U1+M2)**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) + 4*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2*T1*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 - 16*M2*S4**(-2)*(S+U1+M2)**(-1)
     +     + 4*M2*S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*M2*S4**(-1)
     +    *(M2*(S+T1)+T1*U1)**(-1) - 4*M2**2*TG**(-2)*S*S4**(-1)*
     +    (S+U1)**(-1)*(S+U1+M2)**(-1) + 4*M2**2*TG**(-2)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2)*(S+U1+M2)**(-1) - 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MS2 + 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1)*MG2 + 4*M2**2*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * ( 4*M2**2*TG**(-1)*S*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*S*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*S*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2**2*TG**(-1)*S*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 4*M2**2*TG**(-1)*S**2*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*TG**(-1)*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*M2**2*TG**(-1)*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 - 8*M2**2*TG**(-1)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) - 8*M2**2*TG**(-1)*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 8*M2**2*UG**(-1)*S*U1*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 4*M2**2*UG**(-1)*S**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 4*M2**2*UG**(-1)*U1**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 4*M2**2*S*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 4*M2**2*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 4*M2**3*TG**(-2)*
     +    S4**(-1)*(S+U1)**(-1)*(S+U1+M2)**(-1) + 4*M2**3*TG**(-1)*S*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) + 4*M2**3*TG**(-1)*T1*
     +    S4**(-1)*(M2*(S+T1)+T1*U1)**(-2) - 4*TG**(-2)*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*MS2 + 4*TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-2)*MG2
     +     + 16*TG**(-1)*S*S4**(-2) - 4*TG**(-1)*S*S4**(-1)*
     +    (S+T1)**(-1) - 4*TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1)*MS2 + 4*TG**(-1)*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 8*TG**(-1)*S*S4**(-1)*(S+U1)**(-1) - 4*
     +    TG**(-1)*S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MS2 + 4*TG**(-1)*
     +    S*S4**(-1)*(M2*(S+T1)+T1*U1)**(-1)*MG2 - 4*TG**(-1)*T1*U1*
     +    S4**(-1)*(S+U1)**(-2) + 16*TG**(-1)*U1*S4**(-2) - 8*TG**(-1)*
     +    U1*S4**(-1)*(S+U1)**(-1) + 4*TG**(-1)*S4**(-1)*(S+U1)**(-1)*
     +    MS2 - 4*TG**(-1)*S4**(-1)*(S+U1)**(-1)*MG2 - 4*UG**(-2)*S*T1*
     +    U1*S4**(-1)*(S+T1)**(-2)*(S+T1+M2)**(-1)*MS2 + 4*UG**(-2)*S*
     +    T1*U1*S4**(-1)*(S+T1)**(-2)*(S+T1+M2)**(-1)*MG2 )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 4*UG**(-2)*S*T1*U1**2*
     +    S4**(-1)*(S+T1)**(-2)*(S+T1+M2)**(-1) + 4*UG**(-2)*T1*U1*
     +    S4**(-1)*(S+T1)**(-1)*(S+T1+M2)**(-1)*MS2 - 4*UG**(-2)*T1*U1*
     +    S4**(-1)*(S+T1)**(-1)*(S+T1+M2)**(-1)*MG2 + 4*UG**(-2)*T1*
     +    U1**2*S4**(-1)*(S+T1)**(-1)*(S+T1+M2)**(-1) - 4*UG**(-2)*
     +    T1**2*U1*S4**(-1)*(S+T1)**(-2)*(S+T1+M2)**(-1)*MS2 + 4*
     +    UG**(-2)*T1**2*U1*S4**(-1)*(S+T1)**(-2)*(S+T1+M2)**(-1)*MG2
     +     - 4*UG**(-2)*T1**2*U1**2*S4**(-1)*(S+T1)**(-2)*
     +    (S+T1+M2)**(-1) + 4*UG**(-1)*S*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 4*UG**(-1)*S*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*UG**(-1)*S*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) + 32*UG**(-1)*S*T1*
     +    S4**(-2)*(S+T1+M2)**(-1) - 16*UG**(-1)*S*T1*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) - 8*UG**(-1)*S*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) + 8*UG**(-1)*S*T1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * ( 32*UG**(-1)*S*U1*S4**(-2)*
     +    (S+T1+M2)**(-1) - 8*UG**(-1)*S*U1*S4**(-1)*(S+T1)**(-1)*
     +    (S+T1+M2)**(-1) - 16*UG**(-1)*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) + 8*UG**(-1)*S*U1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) - 4*UG**(-1)*S*
     +    S4**(-1)*(S+U1)**(-1) - 4*UG**(-1)*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 + 4*UG**(-1)*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*UG**(-1)*S**2*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) + 16*UG**(-1)*S**2*U1
     +    *S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) + 16*
     +    UG**(-1)*S**2*S4**(-2)*(S+T1+M2)**(-1) - 8*UG**(-1)*S**2*
     +    S4**(-1)*(S+T1)**(-1)*(S+T1+M2)**(-1) - 8*UG**(-1)*S**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) + 8*UG**(-1)*S**3*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*(S+T1+M2)**(-1) + 32*UG**(-1)*T1*U1*
     +    S4**(-2)*(S+T1+M2)**(-1) - 8*UG**(-1)*T1*U1*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 4*UG**(-1)*T1*U1*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) + 4*UG**(-1)*T1*U1**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 4*UG**(-1)*T1*U1**2*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*UG**(-1)*T1**2*
     +    S4**(-2)*(S+T1+M2)**(-1) - 8*UG**(-1)*T1**2*S4**(-1)*
     +    (S+T1)**(-1)*(S+T1+M2)**(-1) - 4*UG**(-1)*U1*S4**(-1)*
     +    (S+U1)**(-1) - 4*UG**(-1)*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)
     +    *MS2 + 4*UG**(-1)*U1*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 
     +    16*UG**(-1)*U1**2*S4**(-2)*(S+T1+M2)**(-1) - 8*UG**(-1)*U1**2
     +    *S4**(-1)*(M2*(S+U1)+T1*U1)**(-1) + 8*S*T1*S4**(-1)*
     +    (S+U1+M2)**(-1)*(M2*(S+T1)+T1*U1)**(-1) + 4*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MS2 - 4*S*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-2)*MG2 - 4*S*S4**(-1)*(S+U1)**(-1)*
     +    (S+U1+M2)**(-1) + 4*T1*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*
     +    MS2 - 4*T1*U1*S4**(-1)*(M2*(S+T1)+T1*U1)**(-2)*MG2 + 16*T1*
     +    S4**(-2)*(S+U1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + N**2*CF*(S4+MS2) * (  - 8*T1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) - 4*U1*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1) + 32*S4**(-2) - 4*S4**(-1)*
     +    (S+T1)**(-1) - 4*S4**(-1)*(S+U1)**(-1) - 4*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MS2 + 4*S4**(-1)*
     +    (M2*(S+T1)+T1*U1)**(-1)*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(36)*N*CF**2 * ( 16*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(37)*N*CF**2 * (  - 16*M2*TG**(-1)*S4**(-1)
     +    *MS2 - 8*M2*S4**(-1) - 8*M2**2*TG**(-1)*S4**(-1) + 16*
     +    TG**(-1)*U1*S4**(-1)*MS2 + 8*TG**(-1)*U1**2*S4**(-1) - 16*
     +    TG**(-1)*MS2 + 16*UG**(-1)*T1*S4**(-1)*MS2 + 8*UG**(-1)*T1**2
     +    *S4**(-1) - 16*UG**(-1)*U1*S4**(-1)*MS2 - 8*UG**(-1)*U1**2*
     +    S4**(-1) - 16*UG**(-1)*MS2 - 8*T1*S4**(-1) - 16*S4**(-1)*MS2
     +     )
     +
      MQQLRH = MQQLRH + ANG4(37)*N**2*CF * ( 8*M2*TG**(-1)*S4**(-1)*MS2
     +     + 4*M2*S4**(-1) + 4*M2**2*TG**(-1)*S4**(-1) - 8*TG**(-1)*U1*
     +    S4**(-1)*MS2 - 4*TG**(-1)*U1**2*S4**(-1) + 8*TG**(-1)*MS2 - 8
     +    *UG**(-1)*T1*S4**(-1)*MS2 - 4*UG**(-1)*T1**2*S4**(-1) + 8*
     +    UG**(-1)*U1*S4**(-1)*MS2 + 4*UG**(-1)*U1**2*S4**(-1) + 8*
     +    UG**(-1)*MS2 + 4*T1*S4**(-1) + 8*S4**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(48)*N*CF**2 * ( 8*M2*TG**(-1)*S - 16*M2*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1) - 8*M2*S*(S+U1+M2)**(-1) - 16*
     +    M2**2*TG**(-1)*S*(S+U1+M2)**(-1) - 8*TG**(-1)*S*U1 - 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 8*UG**(-1)*S*T1**2*
     +    (S+T1+M2)**(-1) - 8*UG**(-1)*S*U1**2*(S+T1+M2)**(-1) - 8*
     +    UG**(-1)*S**2*T1*(S+T1+M2)**(-1) - 8*UG**(-1)*S**2*U1*
     +    (S+T1+M2)**(-1) - 16*UG**(-1)*S**2*(S+T1+M2)**(-1)*MS2 - 8*S*
     +    T1*(S+U1+M2)**(-1) - 8*S**2*(S+U1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(48)*N**2*CF * (  - 4*M2*TG**(-1)*S + 8*M2*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1) + 4*M2*S*(S+U1+M2)**(-1) + 8*
     +    M2**2*TG**(-1)*S*(S+U1+M2)**(-1) + 4*TG**(-1)*S*U1 + 8*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 4*UG**(-1)*S*T1**2*
     +    (S+T1+M2)**(-1) + 4*UG**(-1)*S*U1**2*(S+T1+M2)**(-1) + 4*
     +    UG**(-1)*S**2*T1*(S+T1+M2)**(-1) + 4*UG**(-1)*S**2*U1*
     +    (S+T1+M2)**(-1) + 8*UG**(-1)*S**2*(S+T1+M2)**(-1)*MS2 + 4*S*
     +    T1*(S+U1+M2)**(-1) + 4*S**2*(S+U1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(49)*N*CF**2 * (  - 4*M2*TG**(-2)*S*MS2 + 4
     +    *M2*TG**(-2)*U1*MS2 - 8*M2*TG**(-1)*S + 8*M2*TG**(-1)*U1 + 8*
     +    M2*TG**(-1)*MS2 - 8*M2 - 4*M2**2*TG**(-2)*S + 4*M2**2*
     +    TG**(-2)*U1 + 4*M2**2*TG**(-2)*MS2 - 12*M2**2*TG**(-1) - 4*
     +    M2**3*TG**(-2) - 4*TG**(-1)*S*MS2 + 4*TG**(-1)*U1*MS2 + 8*
     +    UG**(-1)*S*T1*(S+T1+M2)**(-1)*MS2 + 8*UG**(-1)*S*T1**2*
     +    (S+T1+M2)**(-1) - 8*UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MS2 - 8*
     +    UG**(-1)*T1**2*U1*(S+T1+M2)**(-1) - 8*UG**(-1)*T1**2*
     +    (S+T1+M2)**(-1)*MS2 + 8*UG**(-1)*T1**3*(S+T1+M2)**(-1) - 4*S
     +     - 4*T1 + 4*U1 + 4*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(49)*N**2*CF * ( 2*M2*TG**(-2)*S*MS2 - 2*M2
     +    *TG**(-2)*U1*MS2 + 4*M2*TG**(-1)*S - 4*M2*TG**(-1)*U1 - 4*M2*
     +    TG**(-1)*MS2 + 4*M2 + 2*M2**2*TG**(-2)*S - 2*M2**2*TG**(-2)*
     +    U1 - 2*M2**2*TG**(-2)*MS2 + 6*M2**2*TG**(-1) + 2*M2**3*
     +    TG**(-2) + 2*TG**(-1)*S*MS2 - 2*TG**(-1)*U1*MS2 - 2*UG**(-1)*
     +    S*T1*(S+T1+M2)**(-1)*MS2 - 2*UG**(-1)*S*T1**2*(S+T1+M2)**(-1)
     +     + 2*UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MS2 + 2*UG**(-1)*T1**2*U1
     +    *(S+T1+M2)**(-1) + 2*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MS2 - 2*
     +    UG**(-1)*T1**3*(S+T1+M2)**(-1) + 2*S + 2*T1 - 2*U1 - 2*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(50)*N*CF**2 * (  - 16*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 + 16*M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) - 16*M2*TG**(-1)*U1 - 16*M2*S*(S+U1+M2)**(-1)
     +     + 8*M2*(S+U1+M2)**(-1)*MS2 + 8*M2 - 32*M2**2*TG**(-1)*S*
     +    (S+U1+M2)**(-1) + 16*M2**2*TG**(-1) - 8*M2**2*(S+U1+M2)**(-1)
     +     - 16*M2**3*TG**(-1)*(S+U1+M2)**(-1) + 16*TG**(-1)*S*MS2 - 16
     +    *TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 8*TG**(-1)*U1*MS2 + 8*
     +    TG**(-1)*U1**2 - 4*UG**(-2)*S*U1*MS2 - 4*UG**(-2)*S*U1**2 + 4
     +    *UG**(-2)*T1*U1*MS2 + 4*UG**(-2)*T1*U1**2 + 4*UG**(-2)*U1**2*
     +    MS2 - 4*UG**(-2)*U1**3 + 8*S*(S+U1+M2)**(-1)*MS2 + 8*S - 8*
     +    S**2*(S+U1+M2)**(-1) - 8*U1 - 8*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(50)*N**2*CF * ( 4*M2*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 - 4*M2*TG**(-1)*S + 4*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) + 4*M2*TG**(-1)*U1 + 4*M2*S*(S+U1+M2)**(-1)
     +     - 2*M2*(S+U1+M2)**(-1)*MS2 - 2*M2 + 8*M2**2*TG**(-1)*S*
     +    (S+U1+M2)**(-1) - 4*M2**2*TG**(-1) + 2*M2**2*(S+U1+M2)**(-1)
     +     + 4*M2**3*TG**(-1)*(S+U1+M2)**(-1) - 4*TG**(-1)*S*MS2 + 4*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 2*TG**(-1)*U1*MS2 - 2*
     +    TG**(-1)*U1**2 + 2*UG**(-2)*S*U1*MS2 + 2*UG**(-2)*S*U1**2 - 2
     +    *UG**(-2)*T1*U1*MS2 - 2*UG**(-2)*T1*U1**2 - 2*UG**(-2)*U1**2*
     +    MS2 + 2*UG**(-2)*U1**3 - 2*S*(S+U1+M2)**(-1)*MS2 - 2*S + 2*
     +    S**2*(S+U1+M2)**(-1) + 2*U1 + 2*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(63)*N*CF**2 * ( 4*UG**(-2)*U1*S4**(-1)*MS2
     +     + 4*UG**(-2)*U1**2*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(63)*N**2*CF * (  - 2*UG**(-2)*S*U1*
     +    (S+T1+M2)**(-1) - 2*UG**(-2)*S*(S+T1+M2)**(-1)*MS2 - 2*
     +    UG**(-2)*T1*U1*(S+T1+M2)**(-1) - 2*UG**(-2)*T1*
     +    (S+T1+M2)**(-1)*MS2 - 2*UG**(-2)*U1*S4**(-1)*MS2 + 2*UG**(-2)
     +    *U1 - 2*UG**(-2)*U1**2*S4**(-1) + 2*UG**(-2)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(68)*N*CF**2 * ( 4*M2*TG**(-2)*S4**(-1)*MS2
     +     + 8*M2*TG**(-1)*S4**(-1) + 4*M2**2*TG**(-2)*S4**(-1) + 4*
     +    TG**(-1)*S4**(-1)*MS2 + 4*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(68)*N**2*CF * (  - 2*M2*TG**(-2)*S4**(-1)*
     +    MS2 + 4*M2*TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*TG**(-2) - 4*M2*TG**(-1)*S4**(-1)
     +     + 2*M2*TG**(-1)*(S+U1+M2)**(-1) - 2*M2**2*TG**(-2)*S4**(-1)
     +     + 4*M2**2*TG**(-2)*(S+U1+M2)**(-1) - 2*TG**(-2)*MS2 + 2*
     +    TG**(-2)*MG2 - 2*TG**(-1)*S4**(-1)*MS2 - 2*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(72)*N*CF**2 * ( 8*M2*S*MS2 + 8*M2**2*MS2
     +     + 8*S*MS2**2 )
     +
      MQQLRH = MQQLRH + ANG4(73)*N*CF**2 * ( 4*M2*TG**(-1)*S*MS2 - 4*M2
     +    *TG**(-1)*U1*MS2 + 4*M2*S - 4*M2*MS2 + 4*M2**2*TG**(-1)*S - 4
     +    *M2**2*TG**(-1)*U1 - 4*M2**2*TG**(-1)*MS2 + 4*M2**2 + 4*M2**3
     +    *TG**(-1) + 4*S*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(73)*N**2*CF * ( 2*M2*TG**(-1)*S*MS2 + 2*M2
     +    *TG**(-1)*S*MG2 + 2*M2*TG**(-1)*U1*MS2 - 2*M2*TG**(-1)*U1*MG2
     +     - 2*M2*S - 2*M2*MS2 + 2*M2*MG2 - 2*M2**2*TG**(-1)*S + 2*
     +    M2**2*TG**(-1)*U1 + 6*M2**2*TG**(-1)*MS2 + 2*M2**2*TG**(-1)*
     +    MG2 - 2*M2**2 - 2*M2**3*TG**(-1) + 4*TG**(-1)*S*MS2*MG2 + 4*
     +    TG**(-1)*S*MS2**2 - 2*S*MS2 + 2*S*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(74)*N*CF**2 * ( 16*M2*MS2 + 8*S*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(75)*N*CF**2 * ( 4*M2*TG**(-2)*S*MS2 - 4*M2
     +    *TG**(-2)*U1*MS2 - 16*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MS2 + 24*
     +    M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*(S+U1+M2)**(-1) + 64*M2*
     +    TG**(-1)*U1*S4**(-1)*MS2 - 24*M2*TG**(-1)*U1 + 16*M2*TG**(-1)
     +    *U1**2*S4**(-1) + 32*M2*TG**(-1)*S4**(-1)*MS2**2 - 48*M2*
     +    TG**(-1)*MS2 - 16*M2*S*(S+U1+M2)**(-1) + 16*M2*U1*S4**(-1) + 
     +    16*M2*S4**(-1)*MS2 + 8*M2*(S+U1+M2)**(-1)*MS2 + 8*M2 + 4*
     +    M2**2*TG**(-2)*S - 4*M2**2*TG**(-2)*U1 - 4*M2**2*TG**(-2)*MS2
     +     - 32*M2**2*TG**(-1)*S*(S+U1+M2)**(-1) + 16*M2**2*TG**(-1)*U1
     +    *S4**(-1) + 16*M2**2*TG**(-1)*S4**(-1)*MS2 + 20*M2**2*
     +    TG**(-1) - 8*M2**2*(S+U1+M2)**(-1) + 4*M2**3*TG**(-2) - 16*
     +    M2**3*TG**(-1)*(S+U1+M2)**(-1) + 20*TG**(-1)*S*MS2 - 16*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 32*TG**(-1)*U1*S4**(-1)*
     +    MS2**2 - 20*TG**(-1)*U1*MS2 + 16*TG**(-1)*U1**2*S4**(-1)*MS2
     +     - 32*TG**(-1)*MS2**2 + 8*S*(S+U1+M2)**(-1)*MS2 + 12*S - 8*
     +    S**2*(S+U1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(75)*N*CF**2 * ( 16*T1*U1*S4**(-1) + 32*T1*
     +    S4**(-1)*MS2 - 8*T1 + 8*T1**2*S4**(-1) + 48*U1*S4**(-1)*MS2
     +     - 8*U1 + 8*U1**2*S4**(-1) + 32*S4**(-1)*MS2**2 - 36*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(75)*N**2*CF * (  - 2*M2*TG**(-2)*S*MS2 + 2
     +    *M2*TG**(-2)*U1*MS2 + 4*M2*TG**(-1)*S*(S+U1+M2)**(-1)*MS2 - 8
     +    *M2*TG**(-1)*S + 4*M2*TG**(-1)*S**2*(S+U1+M2)**(-1) - 32*M2*
     +    TG**(-1)*U1*S4**(-1)*MS2 + 10*M2*TG**(-1)*U1 - 8*M2*TG**(-1)*
     +    U1**2*S4**(-1) - 16*M2*TG**(-1)*S4**(-1)*MS2**2 + 30*M2*
     +    TG**(-1)*MS2 + 2*M2*TG**(-1)*MG2 + 4*M2*S*(S+U1+M2)**(-1) - 8
     +    *M2*U1*S4**(-1) - 8*M2*S4**(-1)*MS2 - 2*M2*(S+U1+M2)**(-1)*
     +    MS2 - 2*M2 - 2*M2**2*TG**(-2)*S + 2*M2**2*TG**(-2)*U1 + 2*
     +    M2**2*TG**(-2)*MS2 + 8*M2**2*TG**(-1)*S*(S+U1+M2)**(-1) - 8*
     +    M2**2*TG**(-1)*U1*S4**(-1) - 8*M2**2*TG**(-1)*S4**(-1)*MS2 - 
     +    6*M2**2*TG**(-1) + 2*M2**2*(S+U1+M2)**(-1) - 2*M2**3*TG**(-2)
     +     + 4*M2**3*TG**(-1)*(S+U1+M2)**(-1) - 2*TG**(-1)*S*MS2 + 4*
     +    TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 16*TG**(-1)*U1*S4**(-1)*
     +    MS2**2 + 10*TG**(-1)*U1*MS2 - 2*TG**(-1)*U1*MG2 - 8*TG**(-1)*
     +    U1**2*S4**(-1)*MS2 + 16*TG**(-1)*MS2**2 - 2*S*(S+U1+M2)**(-1)
     +    *MS2 )
     +
      MQQLRH = MQQLRH + ANG4(75)*N**2*CF * (  - 4*S + 2*S**2*
     +    (S+U1+M2)**(-1) - 8*T1*U1*S4**(-1) - 16*T1*S4**(-1)*MS2 + 4*
     +    T1 - 4*T1**2*S4**(-1) - 24*U1*S4**(-1)*MS2 + 4*U1 - 4*U1**2*
     +    S4**(-1) - 16*S4**(-1)*MS2**2 + 14*MS2 + 2*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(77)*N**2*CF * ( 4*M2*TG**(-2)*S + 4*M2*
     +    TG**(-2)*U1 - 4*M2*TG**(-2)*MS2 - 4*M2*TG**(-1) - 4*M2**2*
     +    TG**(-2) + 4*TG**(-2)*S*MS2 - 4*TG**(-2)*S*MG2 + 4*TG**(-2)*
     +    U1*MS2 - 4*TG**(-2)*U1*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(78)*N**2*CF * (  - 2*M2*TG**(-2)*S4**(-1)*
     +    MS2 + 4*M2*TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2*TG**(-2) - 4*M2*TG**(-1)*S4**(-1)
     +     + 2*M2*TG**(-1)*(S+U1+M2)**(-1) - 2*M2**2*TG**(-2)*S4**(-1)
     +     + 4*M2**2*TG**(-2)*(S+U1+M2)**(-1) - 2*TG**(-2)*MS2 + 2*
     +    TG**(-2)*MG2 - 2*TG**(-1)*S4**(-1)*MS2 - 2*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(79)*N**2*CF * ( 4*M2*TG**(-2)*S*U1 + 8*M2*
     +    TG**(-2)*S*MS2 - 8*M2*TG**(-2)*S*MG2 + 4*M2*TG**(-2)*S**2 + 4
     +    *M2*TG**(-2)*U1*MS2 - 8*M2*TG**(-2)*U1*MG2 - 2*M2*TG**(-1)*S
     +     + 2*M2*TG**(-1)*U1 - 8*M2*TG**(-1)*MS2 + 4*M2*TG**(-1)*MG2
     +     + 4*M2**2*TG**(-2)*S + 4*M2**2*TG**(-2)*U1 - 4*M2**2*
     +    TG**(-2)*MS2 + 4*M2**2*TG**(-2)*MG2 - 4*M2**2*TG**(-1) - 4*
     +    M2**3*TG**(-2) + 4*TG**(-2)*S*U1*MS2 - 4*TG**(-2)*S*U1*MG2 + 
     +    4*TG**(-2)*S*MS2**2 + 4*TG**(-2)*S*MG2**2 + 4*TG**(-2)*S**2*
     +    MS2 - 4*TG**(-2)*S**2*MG2 + 2*TG**(-1)*S*MS2 + 2*TG**(-1)*S*
     +    MG2 + 2*TG**(-1)*U1*MS2 - 2*TG**(-1)*U1*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(80)*N*CF**2 * ( 8*M2*TG**(-1)*S4**(-1)*MS2
     +     + 8*M2*TG**(-1)*(S+U1+M2)**(-1)*MS2 + 8*M2*S4**(-1) + 8*M2*
     +    (S+U1+M2)**(-1) + 8*M2**2*TG**(-1)*S4**(-1) + 8*M2**2*
     +    TG**(-1)*(S+U1+M2)**(-1) + 8*T1*S4**(-1) + 8*S4**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(80)*N**2*CF * (  - 2*M2*TG**(-2)*S*
     +    (S+U1+M2)**(-1)*MG2 + 2*M2*TG**(-2)*S - 2*M2*TG**(-2)*U1*
     +    S4**(-1)*MS2 - 4*M2*TG**(-2)*U1*S4**(-1)*MG2 - 4*M2*TG**(-2)*
     +    U1 - 4*M2*TG**(-2)*S4**(-1)*MS2*MG2 - 4*M2*TG**(-2)*S4**(-1)*
     +    MS2**2 - 2*M2*TG**(-2)*MS2 + 4*M2*TG**(-2)*MG2 + 4*M2*
     +    TG**(-1)*S*(S+U1+M2)**(-1) - 10*M2*TG**(-1)*S4**(-1)*MS2 + 2*
     +    M2*TG**(-1)*(S+U1+M2)**(-1)*MS2 - 2*M2*TG**(-1)*
     +    (S+U1+M2)**(-1)*MG2 - 4*M2*TG**(-1) - 4*M2*S4**(-1) - 2*M2*
     +    (S+U1+M2)**(-1) + 2*M2**2*TG**(-2)*U1*S4**(-1) + 4*M2**2*
     +    TG**(-2)*(S+U1+M2)**(-1)*MS2 - 2*M2**2*TG**(-2)*
     +    (S+U1+M2)**(-1)*MG2 - 2*M2**2*TG**(-2) - 4*M2**2*TG**(-1)*
     +    S4**(-1) + 6*TG**(-2)*S*MS2 - 2*TG**(-2)*S*MG2 - 4*TG**(-2)*
     +    U1*S4**(-1)*MS2*MG2 - 4*TG**(-2)*U1*S4**(-1)*MS2**2 + 4*
     +    TG**(-2)*MS2*MG2 + 4*TG**(-2)*MS2**2 + 4*TG**(-1)*S*
     +    (S+U1+M2)**(-1)*MS2 - 2*TG**(-1)*S*(S+U1+M2)**(-1)*MG2 - 4*
     +    TG**(-1)*U1*S4**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(80)*N**2*CF * (  - 4*TG**(-1)*U1*S4**(-1)*
     +    MG2 - 4*TG**(-1)*S4**(-1)*MS2*MG2 - 4*TG**(-1)*S4**(-1)*
     +    MS2**2 + 4*TG**(-1)*MG2 + 2*S*(S+U1+M2)**(-1) - 4*T1*S4**(-1)
     +     - 2*U1*S4**(-1) - 10*S4**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(81)*N*CF**2 * ( 4*M2*TG**(-1)*S*MS2 - 4*M2
     +    *TG**(-1)*U1*MS2 + 4*M2*S + 4*M2*T1 + 4*M2**2*TG**(-1)*S - 4*
     +    M2**2*TG**(-1)*U1 - 4*M2**2*TG**(-1)*MS2 + 4*M2**2 + 4*M2**3*
     +    TG**(-1) + 4*S*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(82)*N*CF**2 * ( 4*M2*TG**(-2)*S*MS2 - 4*M2
     +    *TG**(-2)*U1*MS2 + 24*M2*TG**(-1)*S - 16*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) + 16*M2*TG**(-1)*U1*S4**(-1)*MS2 - 24*M2*
     +    TG**(-1)*U1 + 16*M2*TG**(-1)*U1**2*S4**(-1) - 16*M2*TG**(-1)*
     +    MS2 - 8*M2*S*(S+U1+M2)**(-1) + 16*M2 + 4*M2**2*TG**(-2)*S - 4
     +    *M2**2*TG**(-2)*U1 - 4*M2**2*TG**(-2)*MS2 - 16*M2**2*TG**(-1)
     +    *S*(S+U1+M2)**(-1) + 20*M2**2*TG**(-1) + 4*M2**3*TG**(-2) + 
     +    20*TG**(-1)*S*MS2 - 16*TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 - 20
     +    *TG**(-1)*U1*MS2 + 16*TG**(-1)*U1**2*S4**(-1)*MS2 - 8*S*T1*
     +    (S+U1+M2)**(-1) + 12*S - 8*S**2*(S+U1+M2)**(-1) + 12*T1 + 8*
     +    U1*S4**(-1)*MS2 - 8*U1 + 8*U1**2*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(82)*N**2*CF * (  - 2*M2*TG**(-2)*S*MS2 + 2
     +    *M2*TG**(-2)*U1*MS2 - 10*M2*TG**(-1)*S + 8*M2*TG**(-1)*S**2*
     +    (S+U1+M2)**(-1) - 4*M2*TG**(-1)*U1*S4**(-1)*MS2 + 6*M2*
     +    TG**(-1)*U1 - 4*M2*TG**(-1)*U1**2*S4**(-1) + 4*M2*TG**(-1)*
     +    MS2 + 4*M2*S*(S+U1+M2)**(-1) - 6*M2 - 2*M2**2*TG**(-2)*S + 2*
     +    M2**2*TG**(-2)*U1 + 2*M2**2*TG**(-2)*MS2 + 8*M2**2*TG**(-1)*S
     +    *(S+U1+M2)**(-1) - 8*M2**2*TG**(-1) - 2*M2**3*TG**(-2) - 8*
     +    TG**(-1)*S*MS2 + 8*TG**(-1)*S**2*(S+U1+M2)**(-1)*MS2 + 6*
     +    TG**(-1)*U1*MS2 - 4*TG**(-1)*U1**2*S4**(-1)*MS2 + 4*S*T1*
     +    (S+U1+M2)**(-1) - 4*S + 4*S**2*(S+U1+M2)**(-1) - 4*T1 - 2*U1*
     +    S4**(-1)*MS2 + 2*U1 - 2*U1**2*S4**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(85)*N*CF**2 * ( 4*M2*T1 + 4*M2*MS2 + 4*
     +    UG**(-1)*S*U1*MS2 + 4*UG**(-1)*S*U1**2 - 4*UG**(-1)*T1*U1*MS2
     +     - 4*UG**(-1)*T1*U1**2 - 4*UG**(-1)*U1**2*MS2 + 4*UG**(-1)*
     +    U1**3 - 4*S*U1 + 4*T1*U1 + 4*T1*MS2 + 4*U1*MS2 - 4*U1**2 )
     +
      MQQLRH = MQQLRH + ANG4(85)*N**2*CF * ( 2*M2*UG**(-1)*S*MG2 - 4*M2
     +    *UG**(-1)*T1*U1 - 4*M2*UG**(-1)*U1*MS2 + 2*M2*UG**(-1)*U1*MG2
     +     + 2*M2*S + 2*M2*U1 - 2*UG**(-1)*S*U1*MS2 - 2*UG**(-1)*S*
     +    U1**2 - 2*UG**(-1)*T1*U1*MS2 + 2*UG**(-1)*T1*U1*MG2 + 2*
     +    UG**(-1)*T1*U1**2 + 2*UG**(-1)*U1**2*MS2 - 2*UG**(-1)*U1**3
     +     + 2*S*U1 + 2*S*MS2 - 2*T1*U1 - 2*U1*MS2 + 2*U1**2 )
     +
      MQQLRH = MQQLRH + ANG4(86)*N*CF**2 * ( 8*M2*UG**(-1)*T1*S4**(-1)*
     +    MS2 - 8*M2*UG**(-1)*T1 + 8*M2*UG**(-1)*T1**2*S4**(-1) - 8*M2*
     +    UG**(-1)*MS2 + 4*UG**(-2)*S*U1*MS2 + 4*UG**(-2)*S*U1**2 - 4*
     +    UG**(-2)*T1*U1*MS2 - 4*UG**(-2)*T1*U1**2 - 4*UG**(-2)*U1**2*
     +    MS2 + 4*UG**(-2)*U1**3 + 8*UG**(-1)*S*T1 - 8*UG**(-1)*S*T1**2
     +    *(S+T1+M2)**(-1) + 8*UG**(-1)*S*U1 - 8*UG**(-1)*S*U1**2*
     +    (S+T1+M2)**(-1) + 16*UG**(-1)*S*MS2 - 8*UG**(-1)*S**2*T1*
     +    (S+T1+M2)**(-1) - 8*UG**(-1)*S**2*U1*(S+T1+M2)**(-1) - 16*
     +    UG**(-1)*S**2*(S+T1+M2)**(-1)*MS2 + 8*UG**(-1)*T1*U1*S4**(-1)
     +    *MS2 - 8*UG**(-1)*T1*U1 - 16*UG**(-1)*T1*MS2 + 8*UG**(-1)*
     +    T1**2*U1*S4**(-1) + 16*UG**(-1)*T1**2*S4**(-1)*MS2 + 8*
     +    UG**(-1)*U1**2 + 4*T1 + 4*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(86)*N**2*CF * (  - 2*M2*UG**(-1)*T1*
     +    S4**(-1)*MS2 + 2*M2*UG**(-1)*T1 - 2*M2*UG**(-1)*T1**2*
     +    S4**(-1) + 2*M2*UG**(-1)*MS2 - 2*UG**(-2)*S*U1*MS2 - 2*
     +    UG**(-2)*S*U1**2 + 2*UG**(-2)*T1*U1*MS2 + 2*UG**(-2)*T1*U1**2
     +     + 2*UG**(-2)*U1**2*MS2 - 2*UG**(-2)*U1**3 - 4*UG**(-1)*S*T1
     +     + 4*UG**(-1)*S*T1**2*(S+T1+M2)**(-1) - 4*UG**(-1)*S*U1 + 4*
     +    UG**(-1)*S*U1**2*(S+T1+M2)**(-1) - 8*UG**(-1)*S*MS2 + 2*
     +    UG**(-1)*S*MG2 + 4*UG**(-1)*S**2*T1*(S+T1+M2)**(-1) + 4*
     +    UG**(-1)*S**2*U1*(S+T1+M2)**(-1) + 8*UG**(-1)*S**2*
     +    (S+T1+M2)**(-1)*MS2 - 2*UG**(-1)*T1*U1*S4**(-1)*MS2 + 4*
     +    UG**(-1)*T1*MS2 - 2*UG**(-1)*T1**2*U1*S4**(-1) - 4*UG**(-1)*
     +    T1**2*S4**(-1)*MS2 - 4*UG**(-1)*U1*MS2 + 2*UG**(-1)*U1*MG2 - 
     +    4*UG**(-1)*U1**2 + 2*S + 2*U1 )
     +
      MQQLRH = MQQLRH + ANG4(87)*N*CF**2 * ( 8*M2*S*MS2 + 8*M2**2*MS2
     +     + 8*S*MS2**2 )
     +
      MQQLRH = MQQLRH + ANG4(88)*N*CF**2 * ( 4*M2*T1 - 4*M2*U1 + 4*
     +    UG**(-1)*S*U1*MS2 + 4*UG**(-1)*S*U1**2 - 4*UG**(-1)*T1*U1*MS2
     +     - 4*UG**(-1)*T1*U1**2 - 4*UG**(-1)*U1**2*MS2 + 4*UG**(-1)*
     +    U1**3 - 4*S*U1 + 4*T1*U1 + 4*T1*MS2 + 4*U1*MS2 - 4*U1**2 )
     +
      MQQLRH = MQQLRH + ANG4(88)*N**2*CF * ( 4*M2*UG**(-1)*S*MS2 - 2*M2
     +    *UG**(-1)*T1*MG2 + 2*M2*UG**(-1)*U1*MG2 - 2*M2*T1 + 2*M2*U1
     +     - 4*M2*MS2 + 4*M2**2*UG**(-1)*MS2 - 2*UG**(-1)*S*U1*MS2 + 2*
     +    UG**(-1)*S*U1*MG2 - 2*UG**(-1)*S*U1**2 + 4*UG**(-1)*S*MS2*MG2
     +     + 4*UG**(-1)*S*MS2**2 + 2*UG**(-1)*T1*U1*MS2 + 2*UG**(-1)*T1
     +    *U1**2 + 2*UG**(-1)*U1**2*MS2 - 2*UG**(-1)*U1**3 + 2*S*U1 - 2
     +    *T1*U1 - 2*T1*MS2 - 2*U1*MS2 + 2*U1**2 )
     +
      MQQLRH = MQQLRH + ANG4(89)*N*CF**2 * ( 16*M2*MS2 + 8*S*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(90)*N*CF**2 * ( 16*M2*UG**(-1)*T1*S4**(-1)
     +    *MS2 + 8*M2*UG**(-1)*T1**2*S4**(-1) - 16*M2*UG**(-1)*U1*
     +    S4**(-1)*MS2 - 8*M2*UG**(-1)*U1**2*S4**(-1) - 16*M2*UG**(-1)*
     +    MS2 + 4*UG**(-2)*S*U1*MS2 + 4*UG**(-2)*S*U1**2 - 4*UG**(-2)*
     +    T1*U1*MS2 - 4*UG**(-2)*T1*U1**2 - 4*UG**(-2)*U1**2*MS2 + 4*
     +    UG**(-2)*U1**3 + 16*UG**(-1)*S*T1*S4**(-1)*MS2 + 8*UG**(-1)*S
     +    *T1*(S+T1+M2)**(-1)*MS2 + 8*UG**(-1)*S*T1**2*S4**(-1) + 8*
     +    UG**(-1)*S*T1**2*(S+T1+M2)**(-1) - 16*UG**(-1)*S*U1*S4**(-1)*
     +    MS2 - 8*UG**(-1)*S*U1**2*S4**(-1) + 48*UG**(-1)*T1*U1*
     +    S4**(-1)*MS2 - 8*UG**(-1)*T1*U1*(S+T1+M2)**(-1)*MS2 + 8*
     +    UG**(-1)*T1*U1**2*S4**(-1) + 32*UG**(-1)*T1*S4**(-1)*MS2**2
     +     - 24*UG**(-1)*T1*MS2 + 16*UG**(-1)*T1**2*U1*S4**(-1) - 8*
     +    UG**(-1)*T1**2*U1*(S+T1+M2)**(-1) + 32*UG**(-1)*T1**2*
     +    S4**(-1)*MS2 - 8*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MS2 - 16*
     +    UG**(-1)*T1**2 + 8*UG**(-1)*T1**3*S4**(-1) + 8*UG**(-1)*T1**3
     +    *(S+T1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(90)*N*CF**2 * ( 32*UG**(-1)*U1*S4**(-1)*
     +    MS2**2 - 8*UG**(-1)*U1*MS2 + 16*UG**(-1)*U1**2*S4**(-1)*MS2
     +     - 32*UG**(-1)*MS2**2 + 4*T1 - 4*U1 )
     +
      MQQLRH = MQQLRH + ANG4(90)*N**2*CF * (  - 8*M2*UG**(-1)*T1*
     +    S4**(-1)*MS2 - 4*M2*UG**(-1)*T1**2*S4**(-1) + 8*M2*UG**(-1)*
     +    U1*S4**(-1)*MS2 + 4*M2*UG**(-1)*U1**2*S4**(-1) + 12*M2*
     +    UG**(-1)*MS2 - 2*UG**(-2)*S*U1*MS2 - 2*UG**(-2)*S*U1**2 + 2*
     +    UG**(-2)*T1*U1*MS2 + 2*UG**(-2)*T1*U1**2 + 2*UG**(-2)*U1**2*
     +    MS2 - 2*UG**(-2)*U1**3 - 8*UG**(-1)*S*T1*S4**(-1)*MS2 - 2*
     +    UG**(-1)*S*T1*(S+T1+M2)**(-1)*MS2 - 4*UG**(-1)*S*T1**2*
     +    S4**(-1) - 2*UG**(-1)*S*T1**2*(S+T1+M2)**(-1) + 8*UG**(-1)*S*
     +    U1*S4**(-1)*MS2 + 4*UG**(-1)*S*U1**2*S4**(-1) + 4*UG**(-1)*S*
     +    MS2 - 24*UG**(-1)*T1*U1*S4**(-1)*MS2 + 2*UG**(-1)*T1*U1*
     +    (S+T1+M2)**(-1)*MS2 + 2*UG**(-1)*T1*U1 - 4*UG**(-1)*T1*U1**2*
     +    S4**(-1) - 16*UG**(-1)*T1*S4**(-1)*MS2**2 + 14*UG**(-1)*T1*
     +    MS2 - 2*UG**(-1)*T1*MG2 - 8*UG**(-1)*T1**2*U1*S4**(-1) + 2*
     +    UG**(-1)*T1**2*U1*(S+T1+M2)**(-1) - 16*UG**(-1)*T1**2*
     +    S4**(-1)*MS2 + 2*UG**(-1)*T1**2*(S+T1+M2)**(-1)*MS2 + 6*
     +    UG**(-1)*T1**2 )
     +
      MQQLRH = MQQLRH + ANG4(90)*N**2*CF * (  - 4*UG**(-1)*T1**3*
     +    S4**(-1) - 2*UG**(-1)*T1**3*(S+T1+M2)**(-1) - 16*UG**(-1)*U1*
     +    S4**(-1)*MS2**2 + 6*UG**(-1)*U1*MS2 + 2*UG**(-1)*U1*MG2 - 8*
     +    UG**(-1)*U1**2*S4**(-1)*MS2 + 16*UG**(-1)*MS2**2 - 2*T1 + 2*
     +    U1 - 4*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(91)*N**2*CF * (  - 4*M2*UG**(-2)*S*U1 - 4*
     +    M2*UG**(-2)*T1*U1 - 8*M2*UG**(-2)*U1*MS2 + 8*M2*UG**(-2)*U1*
     +    MG2 - 4*M2*UG**(-2)*U1**2 + 4*M2*UG**(-1)*U1 - 4*M2*UG**(-1)*
     +    MS2 - 4*M2**2*UG**(-2)*U1 + 4*UG**(-2)*S*MS2**2 + 4*UG**(-2)*
     +    S*MG2**2 - 4*UG**(-2)*T1*U1*MS2 + 4*UG**(-1)*S*MG2 + 4*
     +    UG**(-1)*T1*MG2 )
     +
      MQQLRH = MQQLRH + ANG4(92)*N*CF**2 * (  - 8*UG**(-1)*S*U1*
     +    (S+T1+M2)**(-1) - 8*UG**(-1)*S*(S+T1+M2)**(-1)*MS2 - 8*
     +    UG**(-1)*T1*U1*(S+T1+M2)**(-1) - 8*UG**(-1)*T1*
     +    (S+T1+M2)**(-1)*MS2 + 8*UG**(-1)*U1*S4**(-1)*MS2 + 8*UG**(-1)
     +    *U1 + 8*UG**(-1)*U1**2*S4**(-1) + 8*UG**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(92)*N**2*CF * ( 2*M2*UG**(-2)*T1*U1*
     +    S4**(-1) + 2*M2*UG**(-2)*T1*S4**(-1)*MS2 + 4*M2*UG**(-2)*U1*
     +    S4**(-1)*MS2 - 2*M2*UG**(-2)*U1 - 2*M2*UG**(-2)*MS2 - 4*
     +    UG**(-2)*S*T1*U1*(S+T1+M2)**(-1) + 4*UG**(-2)*S*T1*
     +    (S+T1+M2)**(-1)*MS2 + 4*UG**(-2)*S**2*(S+T1+M2)**(-1)*MS2 - 4
     +    *UG**(-2)*T1*U1*S4**(-1)*MS2 - 4*UG**(-2)*T1*U1*S4**(-1)*MG2
     +     - 2*UG**(-2)*T1*U1*(S+T1+M2)**(-1)*MS2 - 2*UG**(-2)*T1*U1**2
     +    *(S+T1+M2)**(-1) - 4*UG**(-2)*T1*S4**(-1)*MS2*MG2 - 4*
     +    UG**(-2)*T1*S4**(-1)*MS2**2 - 4*UG**(-2)*T1**2*U1*
     +    (S+T1+M2)**(-1) - 4*UG**(-2)*U1*S4**(-1)*MS2*MG2 - 4*UG**(-2)
     +    *U1*S4**(-1)*MS2**2 + 2*UG**(-2)*U1*MS2 + 4*UG**(-2)*U1*MG2
     +     - 4*UG**(-2)*U1**2*S4**(-1)*MS2 - 2*UG**(-2)*U1**2 + 4*
     +    UG**(-2)*MS2*MG2 + 4*UG**(-2)*MS2**2 + 2*UG**(-1)*S*T1*
     +    (S+T1+M2)**(-1) + 4*UG**(-1)*S*U1*(S+T1+M2)**(-1) + 2*
     +    UG**(-1)*S*(S+T1+M2)**(-1)*MS2 - 2*UG**(-1)*T1*U1*S4**(-1) + 
     +    4*UG**(-1)*T1*U1*(S+T1+M2)**(-1) )
     +
      MQQLRH = MQQLRH + ANG4(92)*N**2*CF * ( 2*UG**(-1)*T1*
     +    (S+T1+M2)**(-1)*MS2 + 2*UG**(-1)*T1 + 2*UG**(-1)*T1**2*
     +    (S+T1+M2)**(-1) - 6*UG**(-1)*U1*S4**(-1)*MS2 - 4*UG**(-1)*
     +    U1**2*S4**(-1) - 2*UG**(-1)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(93)*N**2*CF * (  - 2*UG**(-2)*S*U1*
     +    (S+T1+M2)**(-1) - 2*UG**(-2)*S*(S+T1+M2)**(-1)*MS2 - 2*
     +    UG**(-2)*T1*U1*(S+T1+M2)**(-1) - 2*UG**(-2)*T1*
     +    (S+T1+M2)**(-1)*MS2 - 2*UG**(-2)*U1*S4**(-1)*MS2 + 2*UG**(-2)
     +    *U1 - 2*UG**(-2)*U1**2*S4**(-1) + 2*UG**(-2)*MS2 )
     +
      MQQLRH = MQQLRH + ANG4(95)*N**2*CF * (  - 4*M2*UG**(-2)*U1 - 4*M2
     +    *UG**(-2)*MS2 )
     +
      MQQLRH = MQQLRH + COLO1(9)*N*CF**2*(S4+MS2) * (  - 8*TG**(-2)*S*
     +    T1**2*S4**(-2)*(S+U1)**(-2)*MS2 - 8*TG**(-2)*S*S4**(-2)*MS2
     +     + 8*TG**(-2)*T1*U1*S4**(-2) + 8*TG**(-2)*T1**3*U1*S4**(-2)*
     +    (S+U1)**(-2) - 8*UG**(-2)*S*U1**2*S4**(-2)*(S+T1)**(-2)*MS2
     +     - 8*UG**(-2)*S*S4**(-2)*MS2 + 8*UG**(-2)*T1*U1*S4**(-2) + 8*
     +    UG**(-2)*T1*U1**3*S4**(-2)*(S+T1)**(-2) + 16*S*T1*U1**2*
     +    S4**(-2)*(M2*(S+U1)+T1*U1)**(-2) + 16*S*T1**2*U1*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2) - 8*S*T1**2*U1**2*S4**(-2)*
     +    (S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S*T1**2*U1**2*
     +    S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 16*S*
     +    T1**2*U1**3*S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2) - 8
     +    *S*T1**2*S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 + 16*S*T1**3*
     +    U1**2*S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) - 8*S*
     +    U1**2*S4**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*S**2*T1*U1*
     +    S4**(-2)*(M2*(S+T1)+T1*U1)**(-2) + 8*S**2*T1*U1*S4**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2) )
     +
      MQQLRH = MQQLRH + COLO1(9)*N*CF**2*(S4+MS2) * (  - 16*S**2*T1*
     +    U1**2*S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 + 8*
     +    S**2*T1*U1**3*S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)
     +     - 16*S**2*T1*S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 16*S**2*
     +    T1**2*U1*S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 
     +    8*S**2*T1**3*U1*S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)
     +     - 16*S**2*U1*S4**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*S**3*
     +    T1**2*S4**(-2)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*
     +    S**3*U1**2*S4**(-2)*(S+T1)**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2
     +     - 8*S**3*S4**(-2)*(M2*(S+T1)+T1*U1)**(-2)*MS2 - 8*S**3*
     +    S4**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*T1*U1**3*S4**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 8*T1**3*U1*S4**(-2)*
     +    (M2*(S+T1)+T1*U1)**(-2) + 8*T1**3*U1**3*S4**(-2)*(S+T1)**(-2)
     +    *(M2*(S+T1)+T1*U1)**(-2) + 8*T1**3*U1**3*S4**(-2)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) )

      END IF

      IF (IFL.EQ.0) MQQH = 2*MQPLLH + 2*MQPLRH 
      IF (IFL.EQ.1) MQQH = 2*0.5D0*MQQLLH + 1*MQQLRH

      DSSQQH = ALPHAS**3 * AVG * MQQH *S4/(S4 +MS2)/8.D0 /S**2 *CONV
      RETURN
      END


      REAL*8 FUNCTION DSSQQD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,IFL)
C***  LOG(DEL) PART OF CROSS SECTIONS FOR Q + Q -> SQ + SQ 
C***  IFL = 0 <--> Q NOT Q
C***  IFL = 1 <--> Q  =  Q
C***  EQUAL FLAVOR AND EQUAL HELICITY 
C***  HAS A FACTOR 1/2 FOR IDENTICAL PARTICLES 
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T,T1,TG,SB,S1,MS,MG,MS2,MG2,M2
      REAL*8 MQPLLD, MQPLRD, MQQLLD, MQQLRD, MQQD
      REAL*8 NS,CONV,PI,ZETA2,N,CF,AVG,BETA,XS,DEL,S4MAX,DLDEL1,DLDEL2
      REAL*8 COLO1(1:9)
      INTEGER IFL

      NS = 6.D0
      CONV = 389379660.D0
      PI = 4.D0*ATAN(1.D0)
      ZETA2 = PI**2/6.D0
      N = 3.D0
      CF = (N**2 -1)/2.D0/N
      AVG = (1.D0/2.D0)**2 * (1.D0/3.D0)**2

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      TG = T1 -M2
      SB = S*SQRT(1 -4.D0*MS2/S)
      S1 = 4.D0*MS2 - S
      T = T1 +MS2
      BETA = SQRT(1 -4*MS2/S)
      XS = (1.D0 -BETA)/(1.D0 +BETA)

      DLDEL1 = LOG(S4MAX/MS2) -(S4MAX-DEL)/S4
      DLDEL2 = LOG(S4MAX/MS2)**2 -2.D0*(S4MAX -DEL)/S4*LOG(S4/MS2)

      COLO1(1) = LOG(XS)
      COLO1(3) = LOG(S/MS2)
      COLO1(4) = LOG(-T1/MS2)
      COLO1(8) = LOG((S+T1)/MS2)

      IF (IFL.EQ.0) THEN

      MQPLLD = 0.D0
      MQPLLD = MQPLLD + N*CF**2*DLDEL1 * (  - 16*TG**(-2)*S*MG2 )
     +
      MQPLLD = MQPLLD + N*CF**2*DLDEL2 * ( 16*TG**(-2)*S*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(1)*N*CF**2*DLDEL1 * (  - 64*TG**(-2)*S*
     +    SB**(-1)*MS2*MG2 + 32*TG**(-2)*S**2*SB**(-1)*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(1)*N**2*CF*DLDEL1 * ( 32*TG**(-2)*S*
     +    SB**(-1)*MS2*MG2 - 16*TG**(-2)*S**2*SB**(-1)*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(3)*N*CF**2*DLDEL1 * (  - 32*TG**(-2)*S*
     +    MG2 )
     +
      MQPLLD = MQPLLD + COLO1(3)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(4)*N*CF**2*DLDEL1 * ( 16*TG**(-2)*S*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(4)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*
     +    MG2 )
     +
      MQPLLD = MQPLLD + COLO1(8)*N*CF**2*DLDEL1 * ( 48*TG**(-2)*S*MG2 )
     +
      MQPLLD = MQPLLD + COLO1(8)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*
     +    MG2 )

      MQPLRD = 0.D0
      MQPLRD = MQPLRD + N*CF**2*DLDEL1 * ( 16*TG**(-2)*S*T1 + 16*
     +    TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + N*CF**2*DLDEL2 * (  - 16*TG**(-2)*S*T1 - 16*
     +    TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(1)*N*CF**2*DLDEL1 * ( 64*TG**(-2)*S*T1*
     +    SB**(-1)*MS2 - 32*TG**(-2)*S*T1**2*SB**(-1) + 64*TG**(-2)*S*
     +    SB**(-1)*MS2**2 - 32*TG**(-2)*S**2*T1*SB**(-1) - 32*TG**(-2)*
     +    S**2*SB**(-1)*MS2 + 64*TG**(-2)*T1**2*SB**(-1)*MS2 )
     +
      MQPLRD = MQPLRD + COLO1(1)*N**2*CF*DLDEL1 * (  - 32*TG**(-2)*S*T1
     +    *SB**(-1)*MS2 + 16*TG**(-2)*S*T1**2*SB**(-1) - 32*TG**(-2)*S*
     +    SB**(-1)*MS2**2 + 16*TG**(-2)*S**2*T1*SB**(-1) + 16*TG**(-2)*
     +    S**2*SB**(-1)*MS2 - 32*TG**(-2)*T1**2*SB**(-1)*MS2 )
     +
      MQPLRD = MQPLRD + COLO1(3)*N*CF**2*DLDEL1 * ( 32*TG**(-2)*S*T1 + 
     +    32*TG**(-2)*S*MS2 + 32*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(3)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*T1
     +     - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(4)*N*CF**2*DLDEL1 * (  - 16*TG**(-2)*S*T1
     +     - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(4)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*T1 + 
     +    16*TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(8)*N*CF**2*DLDEL1 * (  - 48*TG**(-2)*S*T1
     +     - 48*TG**(-2)*S*MS2 - 48*TG**(-2)*T1**2 )
     +
      MQPLRD = MQPLRD + COLO1(8)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*T1 + 
     +    16*TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 )

      END IF

      IF (IFL.EQ.1) THEN
      MQQLLD = 0.D0
      MQQLLD = MQQLLD + N*CF**2*DLDEL1 * (  - 16*TG**(-2)*S*MG2 - 16*S*
     +    (S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + N*CF**2*DLDEL2 * ( 16*TG**(-2)*S*MG2 + 16*S*
     +    (S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + CF**2*DLDEL1 * (  - 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + CF**2*DLDEL2 * ( 32*TG**(-1)*S*(S+T1+M2)**(-1)*
     +    MG2 )
     +
      MQQLLD = MQQLLD + COLO1(1)*N*CF*DLDEL1 * ( 64*TG**(-1)*S*SB**(-1)
     +    *(S+T1+M2)**(-1)*MS2*MG2 - 32*TG**(-1)*S**2*SB**(-1)*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(1)*N*CF**2*DLDEL1 * ( 16*TG**(-2)*S*T1*
     +    SB**(-1)*MG2 - 48*TG**(-2)*S*SB**(-1)*MS2*MG2 - 16*TG**(-2)*S
     +    *SB**(-1)*MG2**2 + 32*TG**(-2)*S**2*SB**(-1)*MG2 - 16*
     +    TG**(-1)*S*SB**(-1)*MG2 - 16*S*T1*SB**(-1)*(S+T1+M2)**(-2)*
     +    MG2 - 48*S*SB**(-1)*(S+T1+M2)**(-2)*MS2*MG2 - 16*S*SB**(-1)*
     +    (S+T1+M2)**(-2)*MG2**2 + 16*S*SB**(-1)*(S+T1+M2)**(-1)*MG2 + 
     +    16*S**2*SB**(-1)*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(1)*N**2*CF*DLDEL1 * ( 32*TG**(-2)*S*
     +    SB**(-1)*MS2*MG2 - 16*TG**(-2)*S**2*SB**(-1)*MG2 + 32*S*
     +    SB**(-1)*(S+T1+M2)**(-2)*MS2*MG2 - 16*S**2*SB**(-1)*
     +    (S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(1)*CF**2*DLDEL1 * (  - 64*TG**(-1)*S*
     +    SB**(-1)*(S+T1+M2)**(-1)*MS2*MG2 + 32*TG**(-1)*S**2*SB**(-1)*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(3)*N*CF*DLDEL1 * ( 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(3)*N*CF**2*DLDEL1 * (  - 32*TG**(-2)*S*
     +    MG2 - 32*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(3)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*MG2
     +     + 16*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(3)*CF**2*DLDEL1 * (  - 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(4)*N*CF*DLDEL1 * (  - 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(4)*N*CF**2*DLDEL1 * ( 16*TG**(-2)*S*MG2
     +     + 48*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(4)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*
     +    MG2 - 16*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(4)*CF**2*DLDEL1 * ( 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(8)*N*CF*DLDEL1 * (  - 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(8)*N*CF**2*DLDEL1 * ( 48*TG**(-2)*S*MG2
     +     + 16*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(8)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*
     +    MG2 - 16*S*(S+T1+M2)**(-2)*MG2 )
     +
      MQQLLD = MQQLLD + COLO1(8)*CF**2*DLDEL1 * ( 32*TG**(-1)*S*
     +    (S+T1+M2)**(-1)*MG2 )

      MQQLRD = 0.D0
      MQQLRD = MQQLRD + N*CF**2*DLDEL1 * ( 16*TG**(-2)*S*T1 + 16*
     +    TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 + 16*S*T1*(S+T1+M2)**(-2)
     +     + 16*S*(S+T1+M2)**(-2)*MS2 + 16*T1**2*(S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + N*CF**2*DLDEL2 * (  - 16*TG**(-2)*S*T1 - 16*
     +    TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 - 16*S*T1*(S+T1+M2)**(-2)
     +     - 16*S*(S+T1+M2)**(-2)*MS2 - 16*T1**2*(S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(1)*N*CF**2*DLDEL1 * ( 64*TG**(-2)*S*T1*
     +    SB**(-1)*MS2 - 32*TG**(-2)*S*T1**2*SB**(-1) + 64*TG**(-2)*S*
     +    SB**(-1)*MS2**2 - 32*TG**(-2)*S**2*T1*SB**(-1) - 32*TG**(-2)*
     +    S**2*SB**(-1)*MS2 + 64*TG**(-2)*T1**2*SB**(-1)*MS2 + 64*S*T1*
     +    SB**(-1)*(S+T1+M2)**(-2)*MS2 - 32*S*T1**2*SB**(-1)*
     +    (S+T1+M2)**(-2) + 64*S*SB**(-1)*(S+T1+M2)**(-2)*MS2**2 - 32*
     +    S**2*T1*SB**(-1)*(S+T1+M2)**(-2) - 32*S**2*SB**(-1)*
     +    (S+T1+M2)**(-2)*MS2 + 64*T1**2*SB**(-1)*(S+T1+M2)**(-2)*MS2 )
     +
      MQQLRD = MQQLRD + COLO1(1)*N**2*CF*DLDEL1 * (  - 32*TG**(-2)*S*T1
     +    *SB**(-1)*MS2 + 16*TG**(-2)*S*T1**2*SB**(-1) - 32*TG**(-2)*S*
     +    SB**(-1)*MS2**2 + 16*TG**(-2)*S**2*T1*SB**(-1) + 16*TG**(-2)*
     +    S**2*SB**(-1)*MS2 - 32*TG**(-2)*T1**2*SB**(-1)*MS2 - 32*S*T1*
     +    SB**(-1)*(S+T1+M2)**(-2)*MS2 + 16*S*T1**2*SB**(-1)*
     +    (S+T1+M2)**(-2) - 32*S*SB**(-1)*(S+T1+M2)**(-2)*MS2**2 + 16*
     +    S**2*T1*SB**(-1)*(S+T1+M2)**(-2) + 16*S**2*SB**(-1)*
     +    (S+T1+M2)**(-2)*MS2 - 32*T1**2*SB**(-1)*(S+T1+M2)**(-2)*MS2 )
     +
      MQQLRD = MQQLRD + COLO1(3)*N*CF**2*DLDEL1 * ( 32*TG**(-2)*S*T1 + 
     +    32*TG**(-2)*S*MS2 + 32*TG**(-2)*T1**2 + 32*S*T1*
     +    (S+T1+M2)**(-2) + 32*S*(S+T1+M2)**(-2)*MS2 + 32*T1**2*
     +    (S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(3)*N**2*CF*DLDEL1 * (  - 16*TG**(-2)*S*T1
     +     - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 - 16*S*T1*
     +    (S+T1+M2)**(-2) - 16*S*(S+T1+M2)**(-2)*MS2 - 16*T1**2*
     +    (S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(4)*N*CF**2*DLDEL1 * (  - 16*TG**(-2)*S*T1
     +     - 16*TG**(-2)*S*MS2 - 16*TG**(-2)*T1**2 - 48*S*T1*
     +    (S+T1+M2)**(-2) - 48*S*(S+T1+M2)**(-2)*MS2 - 48*T1**2*
     +    (S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(4)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*T1 + 
     +    16*TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 + 16*S*T1*
     +    (S+T1+M2)**(-2) + 16*S*(S+T1+M2)**(-2)*MS2 + 16*T1**2*
     +    (S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(8)*N*CF**2*DLDEL1 * (  - 48*TG**(-2)*S*T1
     +     - 48*TG**(-2)*S*MS2 - 48*TG**(-2)*T1**2 - 16*S*T1*
     +    (S+T1+M2)**(-2) - 16*S*(S+T1+M2)**(-2)*MS2 - 16*T1**2*
     +    (S+T1+M2)**(-2) )
     +
      MQQLRD = MQQLRD + COLO1(8)*N**2*CF*DLDEL1 * ( 16*TG**(-2)*S*T1 + 
     +    16*TG**(-2)*S*MS2 + 16*TG**(-2)*T1**2 + 16*S*T1*
     +    (S+T1+M2)**(-2) + 16*S*(S+T1+M2)**(-2)*MS2 + 16*T1**2*
     +    (S+T1+M2)**(-2) )


      END IF

      IF (IFL.EQ.0) MQQD = 2*MQPLLD + 2*MQPLRD 
      IF (IFL.EQ.1) MQQD = 2*0.5D0*MQQLLD + 1*MQQLRD

      DSSQQD = ALPHAS**3 * AVG * MQQD /4.D0 /S**2 *CONV
      RETURN
      END

 
      REAL*8 FUNCTION DSSGQ3(ALPHAS,S,T1,S4,MS,MG,SCA)
C***  GIVES THE SCALE DEPENDENCE OF HARD
C***  SCA = Q**2/MS**2    
      IMPLICIT NONE
      REAL*8 ALPHAS,S,T1,TG,U1,S4,MS,MG,MS2,MG2,M2,SCA,N,CF,NS,CONV,AVG
      REAL*8 M2GQ3, MGPLL3, MGPLR3, MGQLL3, MGQLR3

      N = 3.D0
      CF = (N**2 -1.D0)/2.D0/N
      NS = 6.D0
      AVG = (1.D0/2.D0)**2 /N/(N**2 -1.D0)
      CONV = 389379660.D0

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      TG = T1 -M2
      U1 = S4 -S -T1

      MGQLL3 = 0D0
      MGQLL3 = MGQLL3 + N*CF**2 * ( 16*TG**(-2)*S*T1*(S+U1)**(-2)*MG2
     +     + 16*TG**(-2)*S*T1**2*(S+U1)**(-3)*MG2 + 8*TG**(-2)*S*
     +    (S+U1)**(-1)*MG2 + 16*S*T1*U1**2*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S*T1**2*U1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1*U1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1**2*U1*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*U1*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**3*T1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**3*T1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**3*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 )
     +
      MGQLL3 = MGQLL3 + CF**2 * ( 32*TG**(-1)*S*T1*U1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*TG**(-1)*S*T1**2*U1*
     +    (S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S*U1*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*TG**(-1)*S**2*
     +    T1*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*TG**(-1)*
     +    S**2*T1**2*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*
     +    TG**(-1)*S**2*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 )

      MGQLR3 = 0D0
      MGQLR3 = MGQLR3 + N*CF**2 * (  - 16*TG**(-2)*S*T1*(S+U1)**(-2)*
     +    MS2 - 16*TG**(-2)*S*T1**2*(S+U1)**(-3)*MS2 - 8*TG**(-2)*S*
     +    (S+U1)**(-1)*MS2 + 8*TG**(-2)*T1*U1*(S+U1)**(-1) + 16*
     +    TG**(-2)*T1**2*U1*(S+U1)**(-2) + 16*TG**(-2)*T1**3*U1*
     +    (S+U1)**(-3) - 16*S*T1*U1**2*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 16*S*T1*U1**2*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 16*S*T1**2*U1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 32*S*T1**2*U1**2*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 32*S*T1**3*U1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 8*S*U1**2*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 32*S**2*T1*U1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*S**2*T1*U1*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 32*S**2*T1**2*U1*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 16*S**2*T1**2*U1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 16*S**2*T1**3*U1*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2) )
     +
      MGQLR3 = MGQLR3 + N*CF**2 * (  - 16*S**2*U1*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*S**3*T1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*S**3*T1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*S**3*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*T1*U1**3*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 16*T1**2*U1**3*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 16*T1**3*U1**3*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2) )


      MGPLL3 = 0D0
      MGPLL3 = MGPLL3 + N*CF**2 * ( 16*TG**(-2)*S*T1*(S+U1)**(-2)*MG2
     +     + 16*TG**(-2)*S*T1**2*(S+U1)**(-3)*MG2 + 8*TG**(-2)*S*
     +    (S+U1)**(-1)*MG2 + 16*S*T1*U1**2*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S*T1**2*U1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1*U1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1**2*U1*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*U1*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**3*T1*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**3*T1**2*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**3*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 )
     +

      MGPLR3 = MGQLR3

      M2GQ3 = 2.D0*0.5D0*(NS -2.D0)*(MGPLL3 +MGPLR3) 
     +     +2*0.5D0*MGQLL3 + 1*MGQLR3
      DSSGQ3 = ALPHAS**3 * AVG *M2GQ3 *LOG(1.D0/SCA)/8.D0/S**2  *CONV
      RETURN
      END


      REAL*8 FUNCTION DSSGQH(ALPHAS,S,T1,S4,MS,MG,EPS)
C***  CROSS SECTIONS FOR G + Q -> SQ + SQ +QB
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,MS,MG,MS2,MG2,M2,EPS,DEL
      REAL*8 NS,CONV,N,CF,AVG,S4G, S4G2,SYMBU
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 ABP1P1, ABP1M1, ABP1P2, A4P1P2, A4P1P1
      REAL*8 A4M1P2, A4M1P1, A4P0P2, A4P1P0, A4M2P1
      REAL*8 ABP1M2
      REAL*8 C4P1P2, C4P1P1, C4M1P1, C4P1P0, CBP1P1
      REAL*8 M2GQH, MGPLLH, MGPLRH, MGQLLH, MGQLRH
      REAL*8 ANG4(1:104), COLO1(1:9)
           
      CONV = 389379660.D0
      NS = 6.D0

      N = 3.D0
      CF = (N**2 -1.D0)/N/2.D0

      AVG = (1.D0/2.D0)**2 /N /(N**2 -1.D0)

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 - MS2
      U1 = S4 -S -T1
      TG = T1 -M2
      S4G = S4 -M2
      S4G2 = S4G**2 + EPS*MS**4
      SYMBU = (M2*(S+U1)+U1*S4)/((M2*(S+U1)+U1*S4)**2 +(S+U1)*EPS*MG**4)


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +U1)/ANGDEF(1)
      ANGDEF(3) = (S +T1)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(T1 +U1 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(7) = SQRT((T1 +U1)**2 -4.D0*MS2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (T1*S4 -S*(U1+2.D0*MS2))/(S+T1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (U1*S4 -S*(T1+2.D0*MS2))/(S+U1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)

      ANA(1,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(8)
      ANC(1,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(9)
      ANA(1,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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


      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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

C$$$      ANG4(1) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(2) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(3) =A4P1P2(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(4) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(5) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(6) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(7) =A4P1P2(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(8) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(9) =A4P2M2(ANA(2,1),ANB(2,1),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(10) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))

C$$$      ANG4(11) =A4M2P1(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(12) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(13) =A4P2M2(ANA(2,2),ANB(2,2),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(14) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(15) =A4M2P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(16) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(17) =A4P1P1(ANA(3,6),ANB(3,6),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(18) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(19) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(20) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,5),ANB(3,5),ANC(3,5))

      ANG4(21) =A4M1P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG4(22) =A4M1P2(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG4(23) =A4M2P1(ANA(1,5),ANB(1,5),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(24) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(25) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(26) =A4P1P2(ANA(1,5),ANB(1,5),ANA(1,7),ANB(1,7),ANC(1,7))
C$$$      ANG4(27) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,5),ANB(2,5),ANC(2,5))
C$$$      ANG4(28) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(29) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(30) =A4P1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))

C$$$      ANG4(31) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(32) =A4M1P2(ANA(3,9),ANB(3,9),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(33) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG4(34) =A4M1P1(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
C$$$      ANG4(35) =A4P1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(36) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(37) =A4P1P0(ANA(2,1),ANB(2,1),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(38) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(39) =A4P1P0(ANA(2,2),ANB(2,2),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(40) =A4P0P2(ANA(1,1),ANB(1,1),ANA(2,5),ANB(2,5),ANC(2,5))

C$$$      ANG4(41) =A4P1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(42) =A4M1P0(ANA(1,5),ANB(1,5),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG4(43) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,6),ANB(1,6),ANC(1,6))
      ANG4(44) =A4P1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(45)= A4M1P0(ANA(3,6),ANB(3,6),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(46) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG4(47) =A4P1P0(ANA(2,7),ANB(2,7),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(48) =AHP1P1(ANA(1,3),ANB(1,3),ANA(1,4),ANB(1,4),ANC(1,4))
C$$$      ANG4(49) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(50) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,1),ANB(3,1),ANC(3,1))

C$$$      ANG4(51) =ABP2P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(52) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(53) =ABP2P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(54) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,2),ANB(1,2),ANC(1,2))
C$$$      ANG4(55) =ABP2P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(56) =ABP1P2(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(57) =ABP2P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(58) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(59) =ABP2P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(60) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))

C$$$      ANG4(61) =ABP2P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(62) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(63) =ABP1M1(ANA(1,3),ANB(1,3),ANA(1,6),ANB(1,6),ANC(1,6))
C$$$      ANG4(64) =ABP2P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(65) =ABP1P2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(66) =ABP2P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(67) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG4(68) =ABP1M1(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
      ANG4(69) =ABP1M2(ANA(3,4),ANB(3,4),ANA(3,5),ANB(3,5),ANC(3,5))
C$$$      ANG4(70) =ABP2P0(ANA(1,3),ANB(1,3),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG4(71) =ABP2P0(ANA(3,4),ANB(3,4),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(72) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(73) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(74) =A4P1P2(ANA(3,9),ANB(3,9),ANA(3,1),ANB(3,1),ANC(3,1))
C$$$      ANG4(75) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,9),ANB(2,9),ANC(2,9))
C$$$      ANG4(76) =A4P1P1(ANA(3,9),ANB(3,9),ANA(3,2),ANB(3,2),ANC(3,2))
C$$$      ANG4(77) =A4M1P2(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(78) =A4M1P1(ANA(1,5),ANB(1,5),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(79) =A4P0P2(ANA(1,1),ANB(1,1),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(80) =A4P1P0(ANA(3,9),ANB(3,9),ANA(1,1),ANB(1,1),ANC(1,1))

C$$$      ANG4(81) =ABP1P2(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(82) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(83) =ABP1P1(ANA(1,3),ANB(1,3),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG4(84) =ABP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG4(85) = ABP1P2(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8)) 
      ANG4(86) = ABP1P1(ANA(3,4),ANB(3,4),ANA(3,8),ANB(3,8),ANC(3,8)) 
C$$$      ANG4(87) =A4P2P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
C$$$      ANG4(88) =A4P1P2(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))
C$$$      ANG4(89) =A4P1P2(ANA(1,8),ANB(1,8),ANA(1,1),ANB(1,1),ANC(1,1))
C$$$      ANG4(90) =A4P1P1(ANA(2,1),ANB(2,1),ANA(2,8),ANB(2,8),ANC(2,8))

      ANG4(91) = A4P0P2(ANA(1,1),ANB(1,1),ANA(2,8),ANB(2,8),ANC(2,8))
      ANG4(92) = A4P1P0(ANA(1,8),ANB(1,8),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG4(93) = A4M1P1(ANA(3,6),ANB(3,6),ANA(3,8),ANB(3,8),ANC(3,8))
C$$$      ANG4(94) =A4P1P1(ANA(1,8),ANB(1,8),ANA(1,9),ANB(1,9),ANC(1,9))
C$$$      ANG4(95) =A4M1P2(ANA(3,6),ANB(3,6),ANA(3,8),ANB(3,8),ANC(3,8))
C$$$      ANG4(96) =A4P2P2(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
      ANG4(97) =A4P1P2(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
C$$$      ANG4(98) =A4P1P2(ANA(1,8),ANB(1,8),ANA(1,7),ANB(1,7),ANC(1,7))
      ANG4(99) =A4P1P1(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))

C$$$      ANG4(100)=A4P2P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG4(101)=A4P1P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
C$$$      ANG4(102)=A4P1P2(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
      ANG4(103)=A4P1P1(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
      ANG4(104)=A4P1P1(ANA(1,8),ANB(1,8),ANA(1,6),ANB(1,6),ANC(1,6))

      DEL = EPS * MG**4
      IF (MG.GT.MS) THEN
         ANG4(34) =C4M1P1(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
         ANG4(47) =C4P1P0(ANA(2,7),ANB(2,7),ANA(1,1),ANB(1,1),ANC(1,1))
      ANG4(84) =CBP1P1(ANA(3,4),ANB(3,4),ANA(3,7),ANB(3,7),ANC(3,7),DEL)
         ANG4(97) =C4P1P2(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
         ANG4(99) =C4P1P1(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
         ANG4(101)=C4P1P2(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
         ANG4(103)=C4P1P1(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))
       ENDIF

      COLO1(9) = LOG(S4**2/MS2/(S4+MS2))

      MGQLLH = 0D0
      MGQLLH = MGQLLH + N*CF*S4G*S4G2**(-1) * (  - 4*TG**(-1)*S**(-1)*
     +    U1*MG2 - 8*TG**(-1)*S**(-1)*MS2*MG2 + 8*TG**(-1)*S**(-1)*
     +    MG2**2 + 4*TG**(-1)*S*U1**(-1)*MG2 + 8*TG**(-1)*U1**(-1)*MS2*
     +    MG2 + 8*TG**(-1)*U1**(-1)*MG2**2 + 16*S**(-1)*T1*U1**(-1)*MG2
     +     + 16*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + N*CF * ( 8*TG**(-1)*S**(-1)*MG2 + 8*TG**(-1)*
     +    U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + N*CF**2*(S4+MS2)**(-1) * ( 16*M2*S**(-1)*
     +    U1**(-1)*S4*MS2*MG2**(-1) + 16*M2*S**(-1)*U1**(-1)*MS2 - 16*
     +    M2*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) - 4*M2*S**(-1)*MS2*
     +    MG2**(-1) - 8*M2*U1**(-2)*MS2**2*MG2**(-1) - 8*M2**2*S**(-1)*
     +    U1**(-1)*MS2*MG2**(-1) + 4*S**(-1)*T1*MS2*MG2**(-1) - 16*
     +    S**(-1)*U1**(-1)*S4*MS2 + 16*S**(-1)*U1**(-1)*S4*MS2**2*
     +    MG2**(-1) - 8*S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-1) - 8*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 16*S**(-1)*U1**(-1)*MS2**2 - 8*
     +    S**(-1)*U1**(-1)*MS2**3*MG2**(-1) + 4*S**(-1)*S4*MS2*
     +    MG2**(-1) + 4*S**(-1)*MS2 - 4*S**(-1)*MS2**2*MG2**(-1) + 8*
     +    U1**(-2)*S4*MS2**2*MG2**(-1) + 8*U1**(-2)*MS2**2 - 8*U1**(-2)
     +    *MS2**3*MG2**(-1) - 8*U1**(-1)*MS2**2*MG2**(-1) )
     +
      MGQLLH = MGQLLH + N*CF**2*(S4+MS2) * (  - 16*TG**(-2)*S*T1*
     +    S4**(-1)*(S+U1)**(-2)*MG2 - 16*TG**(-2)*S*T1**2*S4**(-1)*
     +    (S+U1)**(-3)*MG2 - 16*S*T1*U1**2*S4**(-1)*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S*T1**2*U1**2*S4**(-1)*
     +    (S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 - 8*S*U1**2*S4**(-1)*(S+U1)**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MG2 - 32*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 32*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**2*
     +    U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**2
     +    *S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**3*T1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**3*T1**2*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 8*S**3*
     +    S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 )
     +
      MGQLLH = MGQLLH + N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)*S**(-1)*
     +    MS2*MG2 - 16*TG**(-1)*S**(-1)*MG2**2 + 16*TG**(-1)*S*U1**(-1)
     +    *MG2 - 16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*U1**(-1)*
     +    MG2**2 + 16*TG**(-1)*MG2 - 4*S**(-1)*T1*MS2*MG2**(-1) - 24*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 8*S**(-1)*U1**(-1)*MS2**3*
     +    MG2**(-1) + 16*S**(-1)*U1**(-1)*MG2**2 + 4*S**(-1)*MS2**2*
     +    MG2**(-1) - 20*S**(-1)*MG2 - 8*U1**(-2)*MS2*MG2 + 8*U1**(-2)*
     +    MS2**3*MG2**(-1) + 8*U1**(-1)*MS2**2*MG2**(-1) + 8*U1**(-1)*
     +    MG2 )
     +
      MGQLLH = MGQLLH + N*CF**2 * (  - 8*M2*S**(-1)*U1**(-1)*MS2*
     +    MG2**(-1) + 8*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) + 8*S**(-1)*
     +    U1**(-1)*MS2 - 16*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 8*
     +    S**(-1)*U1**(-1)*MG2 - 4*S**(-1)*MS2*MG2**(-1) - 8*U1**(-2)*
     +    MS2**2*MG2**(-1) )
     +
      MGQLLH = MGQLLH + N**2*CF*(S4+MS2)**(-1)*S4G*S4G2**(-1) * ( 8*
     +    TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*T1*MS2*MG2 + 4*TG**(-1)*S*
     +    U1**(-1)*MS2*MG2 + 4*TG**(-1)*T1*U1**(-1)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + N**2*CF*(S4+MS2)**(-1) * (  - 4*M2*TG**(-2)*MS2
     +     + 4*M2*TG**(-2)*MS2**2*MG2**(-1) - 2*M2*TG**(-1)*S**(-2)*MS2
     +    *MG2 - 2*M2*TG**(-1)*S**(-2)*MS2*MG2**2*(S4G-S)**(-1) + 2*M2*
     +    TG**(-1)*S**(-2)*MS2**2*MG2*(S4G-S)**(-1) + 2*M2*TG**(-1)*
     +    S**(-2)*MS2**2 - 2*M2*TG**(-1)*S**(-1)*MS2*MG2*(S4G-S)**(-1)
     +     - 2*M2*TG**(-1)*S**(-1)*MS2 - 4*M2*S**(-2)*U1**(-1)*MS2*MG2
     +     - 4*M2*S**(-2)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) - 4*M2*
     +    S**(-2)*U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*
     +    U1**(-1)*MS2**2 - 8*M2*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) - 2*
     +    M2*S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 6*M2*S**(-1)*
     +    U1**(-1)*MS2 + 12*M2*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 2*M2
     +    *S**(-1)*MS2*MG2**(-1) - 2*M2*U1**(-1)*MS2*MG2**(-1) + 4*
     +    M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-1) + 12*TG**(-2)*S*MS2 + 8*
     +    TG**(-2)*T1*MS2 + 4*TG**(-2)*T1*MS2**2*MG2**(-1) + 4*TG**(-2)
     +    *S4*MS2 - 4*TG**(-2)*S4*MS2**2*MG2**(-1) + 4*TG**(-2)*MS2*MG2
     +     - 8*TG**(-2)*MS2**2 )
     +
      MGQLLH = MGQLLH + N**2*CF*(S4+MS2)**(-1) * ( 4*TG**(-2)*MS2**3*
     +    MG2**(-1) + 2*TG**(-1)*S**(-2)*T1*MS2*MG2 + 2*TG**(-1)*
     +    S**(-2)*T1*MS2*MG2**2*(S4G-S)**(-1) + 6*TG**(-1)*S**(-2)*T1*
     +    MS2**2*MG2*(S4G-S)**(-1) + 6*TG**(-1)*S**(-2)*T1*MS2**2 + 2*
     +    TG**(-1)*S**(-2)*S4*MS2*MG2 + 2*TG**(-1)*S**(-2)*S4*MS2*
     +    MG2**2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*S4*MS2**2*MG2*
     +    (S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*S4*MS2**2 + 2*TG**(-1)*
     +    S**(-2)*MS2*MG2**2 + 2*TG**(-1)*S**(-2)*MS2*MG2**3*
     +    (S4G-S)**(-1) - 4*TG**(-1)*S**(-2)*MS2**2*MG2 - 4*TG**(-1)*
     +    S**(-2)*MS2**2*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*
     +    MS2**3*MG2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*MS2**3 + 2*
     +    TG**(-1)*S**(-1)*T1*MS2 - 8*TG**(-1)*S**(-1)*T1*MS2**2*
     +    MG2**(-1) + 2*TG**(-1)*S**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 2*
     +    TG**(-1)*S**(-1)*S4*MS2 + 6*TG**(-1)*S**(-1)*MS2*MG2 + 6*
     +    TG**(-1)*S**(-1)*MS2*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*
     +    S**(-1)*MS2**2*MG2*(S4G-S)**(-1) )
     +
      MGQLLH = MGQLLH + N**2*CF*(S4+MS2)**(-1) * ( 2*TG**(-1)*S**(-1)*
     +    MS2**2 + 4*TG**(-1)*S*MS2*MG2**(-1) + 4*TG**(-1)*T1*MS2*
     +    MG2**(-1) + 2*TG**(-1)*MS2*MG2*(S4G-S)**(-1) + 4*TG**(-1)*MS2
     +     - 8*TG**(-1)*MS2**2*MG2**(-1) + 4*S**(-2)*U1**(-1)*S4*MS2*
     +    MG2 + 4*S**(-2)*U1**(-1)*S4*MS2*MG2**2*(S4G-S)**(-1) + 4*
     +    S**(-2)*U1**(-1)*S4*MS2**2*MG2*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*S4*MS2**2 + 4*S**(-2)*U1**(-1)*MS2*MG2**2 + 4*
     +    S**(-2)*U1**(-1)*MS2*MG2**3*(S4G-S)**(-1) - 4*S**(-2)*
     +    U1**(-1)*MS2**3*MG2*(S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*MS2**3
     +     - 2*S**(-2)*MS2*MG2 - 2*S**(-2)*MS2*MG2**2*(S4G-S)**(-1) - 6
     +    *S**(-2)*MS2**2*MG2*(S4G-S)**(-1) - 6*S**(-2)*MS2**2 + 2*
     +    S**(-1)*U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 6*S**(-1)*
     +    U1**(-1)*S4*MS2 - 12*S**(-1)*U1**(-1)*S4*MS2**2*MG2**(-1) + 4
     +    *S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-1) + 4*S**(-1)*U1**(-1)*
     +    MS2*MG2 + 4*S**(-1)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) - 4*
     +    S**(-1)*U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) )
     +
      MGQLLH = MGQLLH + N**2*CF*(S4+MS2)**(-1) * (  - 12*S**(-1)*
     +    U1**(-1)*MS2**2 + 8*S**(-1)*U1**(-1)*MS2**3*MG2**(-1) - 2*
     +    S**(-1)*S4*MS2*MG2**(-1) + 8*S**(-1)*MS2**2*MG2**(-1) + 2*
     +    U1**(-1)*S4*MS2*MG2**(-1) + 2*U1**(-1)*MS2*MG2*(S4G-S)**(-1)
     +     + 2*U1**(-1)*MS2 )
     +
      MGQLLH = MGQLLH + N**2*CF*S4G*S4G2**(-1) * (  - 12*TG**(-2)*S*MS2
     +     - 4*TG**(-2)*S*MG2 - 8*TG**(-2)*T1*MS2 - 4*TG**(-2)*T1*
     +    MS2**2*MG2**(-1) + 4*TG**(-2)*T1*MG2 + 12*TG**(-2)*MS2*MG2 + 
     +    4*TG**(-2)*MS2**2 - 4*TG**(-2)*MS2**3*MG2**(-1) - 12*TG**(-2)
     +    *MG2**2 - 8*TG**(-1)*S**(-2)*T1*MS2*MG2 - 6*TG**(-1)*S**(-2)*
     +    T1*MS2**2 - 2*TG**(-1)*S**(-2)*T1*MG2**2 + 2*TG**(-1)*S**(-2)
     +    *MS2*MG2**2 - 6*TG**(-1)*S**(-2)*MS2**2*MG2 - 2*TG**(-1)*
     +    S**(-2)*MS2**3 + 6*TG**(-1)*S**(-2)*MG2**3 - 2*TG**(-1)*
     +    S**(-1)*T1*MS2 + 8*TG**(-1)*S**(-1)*T1*MS2**2*MG2**(-1) - 2*
     +    TG**(-1)*S**(-1)*T1*MG2 - 8*TG**(-1)*S**(-1)*MS2*MG2 - 2*
     +    TG**(-1)*S**(-1)*MS2**2 + 10*TG**(-1)*S**(-1)*MG2**2 - 8*
     +    TG**(-1)*S*U1**(-1)*MG2 - 4*TG**(-1)*S*MS2*MG2**(-1) + 4*
     +    TG**(-1)*T1*U1**(-1)*MG2 - 4*TG**(-1)*T1*MS2*MG2**(-1) + 12*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 12*TG**(-1)*U1**(-1)*MG2**2 - 4*
     +    TG**(-1)*MS2 + 8*TG**(-1)*MS2**2*MG2**(-1) - 16*TG**(-1)*MG2
     +     - 4*S**(-2)*U1**(-1)*MS2*MG2**2 )
     +
      MGQLLH = MGQLLH + N**2*CF*S4G*S4G2**(-1) * (  - 4*S**(-2)*
     +    U1**(-1)*MS2**2*MG2 + 4*S**(-2)*U1**(-1)*MS2**3 + 4*S**(-2)*
     +    U1**(-1)*MG2**3 + 8*S**(-2)*MS2*MG2 + 6*S**(-2)*MS2**2 + 2*
     +    S**(-2)*MG2**2 + 12*S**(-1)*U1**(-1)*MS2*MG2 - 8*S**(-1)*
     +    U1**(-1)*MS2**3*MG2**(-1) - 4*S**(-1)*U1**(-1)*MG2**2 - 2*
     +    S**(-1)*MS2 - 8*S**(-1)*MS2**2*MG2**(-1) + 10*S**(-1)*MG2 - 
     +    12*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + N**2*CF * ( 2*M2*TG**(-1)*S**(-2)*MS2*MG2*
     +    (S4G-S)**(-1) - 2*M2*TG**(-1)*S**(-2)*MG2**2*(S4G-S)**(-1) - 
     +    2*M2*TG**(-1)*S**(-1)*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*
     +    U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*U1**(-1)*MG2**2
     +    *(S4G-S)**(-1) + 4*M2*S**(-1)*U1**(-1)*MS2*MG2**(-1) - 2*M2*
     +    S**(-1)*U1**(-1)*MG2*(S4G-S)**(-1) - 4*TG**(-2)*MS2 + 4*
     +    TG**(-2)*MS2**2*MG2**(-1) + 6*TG**(-1)*S**(-2)*T1*MS2*MG2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*T1*MG2**2*(S4G-S)**(-1) - 
     +    2*TG**(-1)*S**(-2)*S4*MS2*MG2*(S4G-S)**(-1) + 2*TG**(-1)*
     +    S**(-2)*S4*MG2**2*(S4G-S)**(-1) - 4*TG**(-1)*S**(-2)*MS2*
     +    MG2**2*(S4G-S)**(-1) + 10*TG**(-1)*S**(-2)*MS2**2*MG2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*MS2**2 - 2*TG**(-1)*
     +    S**(-2)*MG2**2 - 6*TG**(-1)*S**(-2)*MG2**3*(S4G-S)**(-1) + 2*
     +    TG**(-1)*S**(-1)*S4*MG2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-1)*
     +    MS2*MG2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-1)*MS2 - 2*TG**(-1)*
     +    S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + N**2*CF * (  - 6*TG**(-1)*S**(-1)*MG2**2*
     +    (S4G-S)**(-1) - 2*TG**(-1)*MG2*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 4*S**(-2)*U1**(-1)*S4*
     +    MG2**2*(S4G-S)**(-1) - 8*S**(-2)*U1**(-1)*MS2*MG2 + 4*S**(-2)
     +    *U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*
     +    MS2**2 - 4*S**(-2)*U1**(-1)*MG2**2 - 4*S**(-2)*U1**(-1)*
     +    MG2**3*(S4G-S)**(-1) - 6*S**(-2)*MS2*MG2*(S4G-S)**(-1) - 2*
     +    S**(-2)*MG2**2*(S4G-S)**(-1) - 4*S**(-1)*U1**(-1)*S4*MS2*
     +    MG2**(-1) + 2*S**(-1)*U1**(-1)*S4*MG2*(S4G-S)**(-1) - 8*
     +    S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 2*S**(-1)*U1**(-1)*
     +    MS2 + 12*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) - 6*S**(-1)*
     +    U1**(-1)*MG2 - 8*S**(-1)*U1**(-1)*MG2**2*(S4G-S)**(-1) + 2*
     +    S**(-1)*MS2*MG2**(-1) - 2*U1**(-1)*MS2*MG2**(-1) - 2*U1**(-1)
     +    *MG2*(S4G-S)**(-1) )
     +
      MGQLLH = MGQLLH + CF**2*(S4+MS2) * (  - 32*TG**(-1)*S*T1*U1*
     +    S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 - 32*
     +    TG**(-1)*S*T1**2*U1*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 - 16*TG**(-1)*S*U1*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S*
     +    S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 - 32*TG**(-1)*S**2*T1*
     +    S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 - 32*
     +    TG**(-1)*S**2*T1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 - 16*TG**(-1)*S**2*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + CF**2*S4G*S4G2**(-1) * (  - 32*S**(-1)*T1*
     +    U1**(-1)*MG2 + 32*U1**(-2)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N*CF*S4G*S4G2**(-1) * ( 4*S**(-1)*MG2
     +     + 8*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N*CF * ( 4*TG**(-1)*S**(-1)*MG2 + 8*
     +    TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N*CF**2*S4G*S4G2**(-1) * ( 8*TG**(-1)*
     +    S**(-1)*U1*MG2 - 16*TG**(-1)*S**(-1)*MS2*MG2 + 16*TG**(-1)*
     +    S**(-1)*MG2**2 + 32*TG**(-1)*S*U1**(-1)*MG2 + 32*TG**(-1)*
     +    U1**(-1)*MG2**2 + 40*TG**(-1)*MG2 + 16*S**(-1)*MG2 + 32*
     +    U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N*CF**2 * ( 8*TG**(-2)*(S+U1)**(-1)*
     +    MS2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*MG2 - 8
     +    *TG**(-1)*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-2)*
     +    S*MG2 + 2*TG**(-2)*U1*MG2 - 4*TG**(-2)*MS2*MG2 - 4*TG**(-2)*
     +    MG2**2 - 2*TG**(-1)*S**(-1)*U1*MG2 + 4*TG**(-1)*S**(-1)*MS2*
     +    MG2 - 4*TG**(-1)*S**(-1)*MG2**2 - 16*TG**(-1)*S*U1**(-1)*MG2
     +     - 16*TG**(-1)*U1**(-1)*MG2**2 - 14*TG**(-1)*MG2 - 4*S**(-1)*
     +    MG2 - 16*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*N**2*CF * (  - 4*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*MG2**2 - 2*TG**(-2)*MG2 + 
     +    4*TG**(-1)*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*CF**2*S4G*S4G2**(-1) * (  - 8*S**(-1)*
     +    MG2 - 8*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(21)*CF**2 * (  - 8*TG**(-1)*S**(-1)*MG2 - 
     +    8*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(22)*N*CF**2 * (  - 8*TG**(-2)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(23)*N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)
     +    *S**(-1)*MG2 + 16*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(23)*N**2*CF*S4G*S4G2**(-1) * (  - 4*
     +    TG**(-1)*S**(-1)*MG2 - 8*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(34)*N*CF*S4G*S4G2**(-1) * ( 2*TG**(-1)*
     +    S**(-1)*U1*MG2 + 2*TG**(-1)*S**(-1)*MS2*MG2 - 2*TG**(-1)*
     +    S**(-1)*MG2**2 - 2*S**(-1)*T1*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(34)*N*CF * (  - 4*TG**(-1)*S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(34)*CF**2*S4G*S4G2**(-1) * ( 4*S**(-1)*T1*
     +    U1**(-1)*MG2 - 4*S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(43)*N*CF**2 * (  - 8*TG**(-2)*S*MS2*MG2 - 
     +    8*TG**(-2)*MS2*MG2**2 + 8*TG**(-2)*MS2**2*MG2 - 8*TG**(-1)*
     +    MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(43)*CF**2 * ( 16*TG**(-1)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N*CF*S4G*S4G2**(-1) * (  - 2*TG**(-1)*
     +    S*MG2 - 2*TG**(-1)*U1*MG2 + 4*TG**(-1)*MS2*MG2 + 4*TG**(-1)*
     +    MG2**2 + 4*S**(-1)*T1*MG2 + 2*S**(-1)*U1*MG2 - 4*S**(-1)*MS2*
     +    MG2 + 4*S**(-1)*MG2**2 + 8*S*U1**(-1)*MG2 + 16*T1*U1**(-1)*
     +    MG2 + 8*U1**(-1)*MS2*MG2 + 8*U1**(-1)*MG2**2 + 10*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N*CF * ( 2*TG**(-1)*S**(-1)*U1*MG2 - 8
     +    *TG**(-1)*S**(-1)*MS2*MG2 + 8*TG**(-1)*S**(-1)*MG2**2 + 16*
     +    TG**(-1)*S*U1**(-1)*MG2 - 8*TG**(-1)*U1**(-1)*MS2*MG2 + 24*
     +    TG**(-1)*U1**(-1)*MG2**2 + 8*TG**(-1)*(S+U1)**(-1)*MS2*MG2 - 
     +    8*TG**(-1)*(S+U1)**(-1)*MG2**2 + 14*TG**(-1)*MG2 + 4*S**(-1)*
     +    MG2 + 16*U1**(-1)*MG2 - 8*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N*CF**2*S4G*S4G2**(-1) * ( 32*TG**(-1)
     +    *S*U1**(-1)*MG2**2 + 24*TG**(-1)*S*MG2 + 16*TG**(-1)*S**2*
     +    U1**(-1)*MG2 - 16*TG**(-1)*U1**(-1)*MS2**2*MG2 + 16*TG**(-1)*
     +    U1**(-1)*MG2**3 + 8*TG**(-1)*U1*MG2 - 24*TG**(-1)*MS2*MG2 + 
     +    24*TG**(-1)*MG2**2 + 32*S*U1**(-1)*MG2 + 16*T1*U1**(-1)*MG2
     +     + 16*U1**(-1)*MS2*MG2 + 16*U1**(-1)*MG2**2 + 24*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N*CF**2 * ( 4*TG**(-2)*S*(S+U1)**(-1)*
     +    MS2*MG2 - 4*TG**(-2)*S*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*S*MG2
     +     + 16*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 - 8*TG**(-2)*
     +    (S+U1)**(-1)*MS2**2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**3 + 8*
     +    TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MS2*MG2 - 16*
     +    TG**(-1)*(S+U1)**(-1)*MG2**2 - 8*TG**(-1)*MG2 + 8*S**(-1)*MG2
     +     + 16*U1**(-1)*MG2 - 8*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-2)*
     +    S*U1*MG2 - 6*TG**(-2)*S*MS2*MG2 - 2*TG**(-2)*S*MG2**2 + 2*
     +    TG**(-2)*S**2*MG2 + 4*TG**(-2)*MS2**2*MG2 - 4*TG**(-2)*MG2**3
     +     - 16*TG**(-1)*S*U1**(-1)*MG2**2 - 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*S**2*U1**(-1)*MG2 + 8*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 8*TG**(-1)*U1**(-1)*MG2**3 - 2*TG**(-1)*U1*MG2 + 6*
     +    TG**(-1)*MS2*MG2 - 14*TG**(-1)*MG2**2 - 16*S*U1**(-1)*MG2 - 8
     +    *T1*U1**(-1)*MG2 - 8*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MG2**2 - 
     +    10*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*N**2*CF * (  - 2*TG**(-2)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 2*TG**(-2)*S*(S+U1)**(-1)*MG2**2 - 2*
     +    TG**(-2)*S*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 + 4*
     +    TG**(-2)*(S+U1)**(-1)*MS2**2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*
     +    MG2**3 + 2*TG**(-1)*S*(S+U1)**(-1)*MG2 - 8*TG**(-1)*
     +    (S+U1)**(-1)*MS2*MG2 + 8*TG**(-1)*(S+U1)**(-1)*MG2**2 - 4*
     +    S**(-1)*MG2 - 8*U1**(-1)*MG2 + 4*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*CF**2*S4G*S4G2**(-1) * (  - 8*S**(-1)*
     +    T1*MG2 - 4*S**(-1)*U1*MG2 + 8*S**(-1)*MS2*MG2 - 8*S**(-1)*
     +    MG2**2 - 8*S*U1**(-1)*MG2 - 16*T1*U1**(-1)*MG2 - 8*U1**(-1)*
     +    MS2*MG2 - 8*U1**(-1)*MG2**2 - 12*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(44)*CF**2 * (  - 4*TG**(-1)*S**(-1)*U1*MG2
     +     + 16*TG**(-1)*S**(-1)*MS2*MG2 - 16*TG**(-1)*S**(-1)*MG2**2
     +     - 16*TG**(-1)*S*U1**(-1)*MG2 + 8*TG**(-1)*U1**(-1)*MS2*MG2
     +     - 24*TG**(-1)*U1**(-1)*MG2**2 - 16*TG**(-1)*(S+U1)**(-1)*MS2
     +    *MG2 + 16*TG**(-1)*(S+U1)**(-1)*MG2**2 - 4*TG**(-1)*MG2 - 8*
     +    S**(-1)*MG2 - 16*U1**(-1)*MG2 + 16*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(47)*N*CF*S4G*S4G2**(-1) * ( 2*TG**(-1)*
     +    S**(-1)*U1*MS2*MG2 - 2*TG**(-1)*S**(-1)*U1*MG2**2 - 16*
     +    TG**(-1)*S**(-1)*MS2*MG2**2 + 8*TG**(-1)*S**(-1)*MS2**2*MG2
     +     + 8*TG**(-1)*S**(-1)*MG2**3 - 2*TG**(-1)*U1*MG2 - 2*TG**(-1)
     +    *MS2*MG2 - 6*TG**(-1)*MG2**2 - 4*S**(-1)*T1*U1**(-1)*MS2*MG2
     +     + 4*S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*MG2 + 4*
     +    S**(-1)*T1**2*U1**(-1)*MG2 - 6*S**(-1)*MS2*MG2 + 6*S**(-1)*
     +    MG2**2 - 4*T1*U1**(-1)*MG2 - 6*U1**(-1)*MS2*MG2 - 2*U1**(-1)*
     +    MG2**2 - 2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(47)*N*CF * ( 6*TG**(-1)*S**(-1)*U1*MG2 - 
     +    18*TG**(-1)*S**(-1)*MS2*MG2 + 18*TG**(-1)*S**(-1)*MG2**2 + 4*
     +    TG**(-1)*MG2 + 10*S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(47)*N*CF**2 * ( 8*S**(-1)*T1*U1**(-1)*MG2
     +     )
     +
      MGQLLH = MGQLLH + ANG4(47)*N**2*CF * (  - 2*S**(-1)*T1*U1**(-1)*
     +    MG2 - 2*S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(47)*CF**2*S4G*S4G2**(-1) * ( 8*S**(-1)*T1*
     +    U1**(-1)*MS2*MG2 - 8*S**(-1)*T1*U1**(-1)*MG2**2 - 4*S**(-1)*
     +    T1*MG2 - 8*S**(-1)*T1**2*U1**(-1)*MG2 - 4*S**(-1)*U1*MG2 - 4*
     +    S**(-1)*MS2*MG2 + 4*S**(-1)*MG2**2 + 4*T1*U1**(-1)*MG2 + 12*
     +    U1**(-1)*MS2*MG2 + 4*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(47)*CF**2 * (  - 4*TG**(-1)*S**(-1)*U1*MG2
     +     + 4*TG**(-1)*S**(-1)*MS2*MG2 - 4*TG**(-1)*S**(-1)*MG2**2 - 8
     +    *TG**(-1)*MG2 - 4*S**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N*CF*S4G*S4G2**(-1) * ( 8*S**(-1)*MG2
     +     + 4*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N*CF * ( 8*TG**(-1)*S**(-1)*MG2 + 4*
     +    TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N*CF**2*S4G*S4G2**(-1) * (  - 32*
     +    TG**(-1)*S**(-1)*MS2*MG2 + 32*TG**(-1)*S**(-1)*MG2**2 + 8*
     +    TG**(-1)*S*U1**(-1)*MG2 - 16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*
     +    TG**(-1)*U1**(-1)*MG2**2 + 8*TG**(-1)*MG2 + 32*S**(-1)*MG2 + 
     +    16*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N*CF**2 * ( 8*TG**(-2)*(S+U1)**(-1)*
     +    MS2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*MG2 - 8
     +    *TG**(-1)*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N**2*CF*S4G*S4G2**(-1) * (  - 2*
     +    TG**(-2)*S*MG2 - 2*TG**(-2)*U1*MG2 + 4*TG**(-2)*MS2*MG2 - 4*
     +    TG**(-2)*MG2**2 + 16*TG**(-1)*S**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    S**(-1)*MG2**2 - 2*TG**(-1)*S*U1**(-1)*MG2 + 4*TG**(-1)*
     +    U1**(-1)*MS2*MG2 - 4*TG**(-1)*U1**(-1)*MG2**2 - 10*TG**(-1)*
     +    MG2 - 16*S**(-1)*MG2 - 4*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*N**2*CF * (  - 4*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*MG2**2 + 2*TG**(-2)*MG2 + 
     +    4*TG**(-1)*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*CF**2*S4G*S4G2**(-1) * (  - 8*S**(-1)*
     +    MG2 - 8*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(68)*CF**2 * (  - 8*TG**(-1)*S**(-1)*MG2 - 
     +    8*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(69)*N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)
     +    *S**(-1)*MG2 + 16*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(69)*N**2*CF*S4G*S4G2**(-1) * (  - 8*
     +    TG**(-1)*S**(-1)*MG2 - 4*TG**(-1)*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(84)*N*CF * (  - 32*TG**(-1)*S**(-1)*U1*MS2
     +    *MG2 + 32*TG**(-1)*S**(-1)*U1*MG2**2 + 8*TG**(-1)*S**(-1)*
     +    U1**2*MG2 - 64*TG**(-1)*S**(-1)*MS2*MG2**2 + 32*TG**(-1)*
     +    S**(-1)*MS2**2*MG2 + 32*TG**(-1)*S**(-1)*MG2**3 - 2*TG**(-1)*
     +    S*U1**(-1)*MS2*MG2 + 2*TG**(-1)*S*U1**(-1)*MG2**2 - 2*
     +    TG**(-1)*S*(S+U1)**(-1)*MS2*MG2 + 2*TG**(-1)*S*(S+U1)**(-1)*
     +    MG2**2 + 2*TG**(-1)*S*MG2 - 16*TG**(-1)*U1**(-1)*MS2*MG2**2
     +     + 8*TG**(-1)*U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*U1**(-1)*
     +    MG2**3 + 8*TG**(-1)*U1*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MS2*
     +    MG2**2 - 8*TG**(-1)*(S+U1)**(-1)*MS2**2*MG2 - 8*TG**(-1)*
     +    (S+U1)**(-1)*MG2**3 - 16*TG**(-1)*MS2*MG2 + 16*TG**(-1)*
     +    MG2**2 + 8*S**(-1)*T1*MG2 + 16*S**(-1)*U1*MG2 - 24*S**(-1)*
     +    MS2*MG2 + 24*S**(-1)*MG2**2 + 2*S*(S+U1)**(-1)*MG2 - 4*T1*
     +    (S+U1)**(-1)*MG2 - 4*U1**(-1)*MS2*MG2 + 4*U1**(-1)*MG2**2 + 8
     +    *(S+U1)**(-1)*MS2*MG2 - 8*(S+U1)**(-1)*MG2**2 + 4*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(84)*N*CF**2 * (  - 32*S**(-1)*T1*U1**(-1)*
     +    MS2*MG2 + 32*S**(-1)*T1*U1**(-1)*MG2**2 + 32*S**(-1)*T1*MG2
     +     + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 32*S**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 16*S**(-1)*U1**(-1)*MS2**2*MG2 + 16*S**(-1)*U1**(-1)
     +    *MG2**3 + 16*S**(-1)*U1*MG2 - 32*S**(-1)*MS2*MG2 + 32*S**(-1)
     +    *MG2**2 + 8*S*T1*U1**(-1)*(S+U1)**(-1)*MG2 + 16*T1*U1**(-1)*
     +    (S+U1)**(-1)*MS2*MG2 - 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 - 8
     +    *T1*U1**(-1)*MG2 - 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 8*
     +    U1**(-1)*MS2*MG2 + 8*U1**(-1)*MG2**2 + 8*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(84)*N**2*CF * ( 16*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 16*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1*MG2 - 8*
     +    S**(-1)*T1**2*U1**(-1)*MG2 + 16*S**(-1)*U1**(-1)*MS2*MG2**2
     +     - 8*S**(-1)*U1**(-1)*MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MG2**3
     +     - 8*S**(-1)*U1*MG2 + 16*S**(-1)*MS2*MG2 - 16*S**(-1)*MG2**2
     +     - 2*S*T1*U1**(-1)*(S+U1)**(-1)*MG2 + 2*S*U1**(-2)*MS2*MG2 - 
     +    2*S*U1**(-2)*MG2**2 - 2*S*U1**(-1)*MG2 + 4*T1*U1**(-2)*MS2*
     +    MG2 - 4*T1*U1**(-2)*MG2**2 - 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*
     +    MG2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 - 4*T1*U1**(-1)*MG2
     +     + 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 + 8*U1**(-2)*MS2*MG2**2
     +     - 4*U1**(-2)*MS2**2*MG2 - 4*U1**(-2)*MG2**3 + 12*U1**(-1)*
     +    MS2*MG2 - 12*U1**(-1)*MG2**2 - 8*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(84)*CF**2 * ( 32*TG**(-1)*S**(-1)*U1*MS2*
     +    MG2 - 32*TG**(-1)*S**(-1)*U1*MG2**2 - 8*TG**(-1)*S**(-1)*
     +    U1**2*MG2 + 64*TG**(-1)*S**(-1)*MS2*MG2**2 - 32*TG**(-1)*
     +    S**(-1)*MS2**2*MG2 - 32*TG**(-1)*S**(-1)*MG2**3 + 4*TG**(-1)*
     +    S*(S+U1)**(-1)*MS2*MG2 - 4*TG**(-1)*S*(S+U1)**(-1)*MG2**2 - 4
     +    *TG**(-1)*U1*MG2 - 32*TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 + 16*
     +    TG**(-1)*(S+U1)**(-1)*MS2**2*MG2 + 16*TG**(-1)*(S+U1)**(-1)*
     +    MG2**3 - 8*S**(-1)*T1*MG2 - 16*S**(-1)*U1*MG2 + 24*S**(-1)*
     +    MS2*MG2 - 24*S**(-1)*MG2**2 - 4*S*(S+U1)**(-1)*MG2 + 8*T1*
     +    (S+U1)**(-1)*MG2 - 16*(S+U1)**(-1)*MS2*MG2 + 16*(S+U1)**(-1)*
     +    MG2**2 + 4*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(85)*N*CF**2 * (  - 4*S*U1**(-1)*MS2*MG2 + 
     +    4*S*U1**(-1)*MG2**2 - 4*S*MG2 - 8*T1*U1**(-1)*MS2*MG2 + 8*T1*
     +    U1**(-1)*MG2**2 - 4*T1*MG2 - 16*U1**(-1)*MS2*MG2**2 + 8*
     +    U1**(-1)*MS2**2*MG2 + 8*U1**(-1)*MG2**3 + 4*MS2*MG2 - 4*
     +    MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(86)*N*CF*S4G*S4G2**(-1) * (  - 6*TG**(-1)*
     +    S*MS2*MG2 + 6*TG**(-1)*S*MG2**2 - 4*TG**(-1)*U1*MS2*MG2 + 4*
     +    TG**(-1)*U1*MG2**2 - 16*TG**(-1)*MS2*MG2**2 + 8*TG**(-1)*
     +    MS2**2*MG2 + 8*TG**(-1)*MG2**3 - 16*S**(-1)*T1*MS2*MG2 + 16*
     +    S**(-1)*T1*MG2**2 + 8*S**(-1)*T1**2*MG2 - 16*S**(-1)*MS2*
     +    MG2**2 + 8*S**(-1)*MS2**2*MG2 + 8*S**(-1)*MG2**3 - 2*S*
     +    U1**(-1)*MS2*MG2 + 2*S*U1**(-1)*MG2**2 + 4*S*MG2 - 4*T1*
     +    U1**(-1)*MS2*MG2 + 4*T1*U1**(-1)*MG2**2 + 10*T1*MG2 - 8*
     +    U1**(-1)*MS2*MG2**2 + 4*U1**(-1)*MS2**2*MG2 + 4*U1**(-1)*
     +    MG2**3 + 2*U1*MG2 - 14*MS2*MG2 + 14*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(86)*N*CF * ( 2*TG**(-1)*S*(S+U1)**(-1)*MS2
     +    *MG2 - 2*TG**(-1)*S*(S+U1)**(-1)*MG2**2 + 16*TG**(-1)*
     +    (S+U1)**(-1)*MS2*MG2**2 - 8*TG**(-1)*(S+U1)**(-1)*MS2**2*MG2
     +     - 8*TG**(-1)*(S+U1)**(-1)*MG2**3 + 4*TG**(-1)*MS2*MG2 - 4*
     +    TG**(-1)*MG2**2 - 2*S*(S+U1)**(-1)*MG2 - 4*T1*(S+U1)**(-1)*
     +    MG2 + 8*(S+U1)**(-1)*MS2*MG2 - 8*(S+U1)**(-1)*MG2**2 - 2*MG2
     +     )
     +
      MGQLLH = MGQLLH + ANG4(86)*N*CF**2 * ( 32*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 32*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1**2*
     +    U1**(-1)*MG2 + 32*S**(-1)*U1**(-1)*MS2*MG2**2 - 16*S**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 16*S**(-1)*U1**(-1)*MG2**3 - 8*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 + 4*S*U1**(-1)*MG2 - 16*T1*U1**(-1)
     +    *(S+U1)**(-1)*MS2*MG2 + 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 + 
     +    16*T1*U1**(-1)*MG2 - 16*T1*(S+U1)**(-1)*MG2 + 16*T1**2*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-1)*MS2*MG2 + 8*U1**(-1)*
     +    MG2**2 - 4*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(86)*N**2*CF * (  - 16*S**(-1)*T1*U1**(-1)*
     +    MS2*MG2 + 16*S**(-1)*T1*U1**(-1)*MG2**2 + 8*S**(-1)*T1**2*
     +    U1**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2*MG2**2 + 8*S**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MG2**3 + 2*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 2*S*U1**(-2)*MS2*MG2 + 2*S*
     +    U1**(-2)*MG2**2 - 4*T1*U1**(-2)*MS2*MG2 + 4*T1*U1**(-2)*
     +    MG2**2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2**2 + 4*T1*(S+U1)**(-1)*MG2 - 4*T1**2*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-2)*MS2*MG2**2 + 4*
     +    U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 - 4*U1**(-1)*MS2*MG2
     +     + 4*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(86)*CF**2*S4G*S4G2**(-1) * ( 16*S**(-1)*T1
     +    *MS2*MG2 - 16*S**(-1)*T1*MG2**2 - 8*S**(-1)*T1**2*MG2 + 16*
     +    S**(-1)*MS2*MG2**2 - 8*S**(-1)*MS2**2*MG2 - 8*S**(-1)*MG2**3
     +     + 4*S*U1**(-1)*MS2*MG2 - 4*S*U1**(-1)*MG2**2 + 8*T1*U1**(-1)
     +    *MS2*MG2 - 8*T1*U1**(-1)*MG2**2 - 4*T1*MG2 + 16*U1**(-1)*MS2*
     +    MG2**2 - 8*U1**(-1)*MS2**2*MG2 - 8*U1**(-1)*MG2**3 + 4*MS2*
     +    MG2 - 4*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(86)*CF**2 * ( 4*TG**(-1)*S*U1**(-1)*MS2*
     +    MG2 - 4*TG**(-1)*S*U1**(-1)*MG2**2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 4*TG**(-1)*S*(S+U1)**(-1)*MG2**2 + 8*
     +    TG**(-1)*S*MG2 + 32*TG**(-1)*U1**(-1)*MS2*MG2**2 - 16*
     +    TG**(-1)*U1**(-1)*MS2**2*MG2 - 16*TG**(-1)*U1**(-1)*MG2**3 - 
     +    32*TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 + 16*TG**(-1)*
     +    (S+U1)**(-1)*MS2**2*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MG2**3 - 
     +    16*TG**(-1)*MS2*MG2 + 16*TG**(-1)*MG2**2 + 4*S*(S+U1)**(-1)*
     +    MG2 + 8*T1*(S+U1)**(-1)*MG2 + 8*U1**(-1)*MS2*MG2 - 8*U1**(-1)
     +    *MG2**2 - 16*(S+U1)**(-1)*MS2*MG2 + 16*(S+U1)**(-1)*MG2**2 + 
     +    8*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(91)*N*CF**2 * (  - 8*S*U1**(-2)*MS2*MG2 - 
     +    4*S*U1**(-1)*MG2 - 8*T1*U1**(-2)*MS2*MG2 - 4*T1*U1**(-1)*MG2
     +     - 8*U1**(-2)*MS2*MG2**2 + 8*U1**(-2)*MS2**2*MG2 + 8*U1**(-1)
     +    *MS2*MG2 - 8*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*N*CF*S4G*S4G2**(-1) * ( 8*TG**(-1)*S*
     +    U1**(-1)*MG2**2 + 4*TG**(-1)*S*MG2 + 2*TG**(-1)*S**2*U1**(-1)
     +    *MG2 - 8*TG**(-1)*U1**(-1)*MS2**2*MG2 + 8*TG**(-1)*U1**(-1)*
     +    MG2**3 - 4*TG**(-1)*MS2*MG2 + 4*TG**(-1)*MG2**2 - 4*S**(-1)*
     +    T1*U1**(-1)*MS2*MG2 + 4*S**(-1)*T1*U1**(-1)*MG2**2 - 2*
     +    S**(-1)*T1*MG2 + 4*S**(-1)*T1**2*U1**(-1)*MG2 + 2*S**(-1)*MS2
     +    *MG2 - 2*S**(-1)*MG2**2 + 8*S*U1**(-1)*MG2 + 10*T1*U1**(-1)*
     +    MG2 - 2*U1**(-1)*MS2*MG2 + 10*U1**(-1)*MG2**2 + 2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*N*CF * ( 6*TG**(-1)*S*U1**(-1)*MG2 - 
     +    10*TG**(-1)*U1**(-1)*MS2*MG2 + 10*TG**(-1)*U1**(-1)*MG2**2 + 
     +    6*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*N*CF**2 * ( 8*S**(-1)*T1*U1**(-1)*MG2
     +     - 8*U1**(-2)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*N**2*CF * (  - 2*S**(-1)*T1*U1**(-1)*
     +    MG2 - 2*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*CF**2*S4G*S4G2**(-1) * ( 8*S**(-1)*T1*
     +    U1**(-1)*MS2*MG2 - 8*S**(-1)*T1*U1**(-1)*MG2**2 - 8*S**(-1)*
     +    T1**2*U1**(-1)*MG2 - 4*S**(-1)*MS2*MG2 + 4*S**(-1)*MG2**2 + 
     +    16*S*U1**(-2)*MS2*MG2 + 16*T1*U1**(-2)*MS2*MG2 - 8*T1*
     +    U1**(-1)*MG2 + 16*U1**(-2)*MS2*MG2**2 - 16*U1**(-2)*MS2**2*
     +    MG2 - 4*U1**(-1)*MS2*MG2 + 4*U1**(-1)*MG2**2 - 4*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(92)*CF**2 * (  - 4*TG**(-1)*S*U1**(-1)*MG2
     +     + 4*TG**(-1)*U1**(-1)*MS2*MG2 - 4*TG**(-1)*U1**(-1)*MG2**2
     +     - 4*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(93)*N*CF*S4G*S4G2**(-1) * ( 2*TG**(-1)*S*
     +    U1**(-1)*MG2 - 2*TG**(-1)*U1**(-1)*MS2*MG2 + 2*TG**(-1)*
     +    U1**(-1)*MG2**2 + 2*S**(-1)*T1*U1**(-1)*MG2 + 4*U1**(-1)*MG2
     +     )
     +
      MGQLLH = MGQLLH + ANG4(93)*CF**2*S4G*S4G2**(-1) * (  - 4*S**(-1)*
     +    T1*U1**(-1)*MG2 - 4*U1**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(97)*N**2*CF * (  - 8*S*U1**(-1)*MS2*MG2 - 
     +    4*S*MG2 - 8*T1*U1**(-1)*MS2*MG2 - 4*T1*MG2 - 8*U1**(-1)*MS2*
     +    MG2**2 + 8*U1**(-1)*MS2**2*MG2 + 8*MS2*MG2 - 8*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(99)*N*CF*S4G*S4G2**(-1) * (  - 4*S**(-1)*
     +    T1*MS2*MG2 + 4*S**(-1)*T1*MG2**2 - 2*S**(-1)*U1*MS2*MG2 + 2*
     +    S**(-1)*U1*MG2**2 - 8*S**(-1)*MS2*MG2**2 + 4*S**(-1)*MS2**2*
     +    MG2 + 4*S**(-1)*MG2**3 + 6*S*U1**(-1)*MS2*MG2 + 2*S*U1**(-1)*
     +    MG2**2 + 2*S*MG2 + 4*T1*U1**(-1)*MS2*MG2 + 4*T1*U1**(-1)*
     +    MG2**2 - 4*U1**(-1)*MS2**2*MG2 + 4*U1**(-1)*MG2**3 - 2*U1*MG2
     +     - 4*MS2*MG2 + 4*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(99)*N*CF * ( 2*TG**(-1)*S**(-1)*U1*MS2*MG2
     +     - 2*TG**(-1)*S**(-1)*U1*MG2**2 + 8*TG**(-1)*S**(-1)*MS2**2*
     +    MG2 - 8*TG**(-1)*S**(-1)*MG2**3 + 2*TG**(-1)*S*U1**(-1)*MS2*
     +    MG2 - 2*TG**(-1)*S*U1**(-1)*MG2**2 + 2*TG**(-1)*S*MG2 + 16*
     +    TG**(-1)*U1**(-1)*MS2*MG2**2 - 8*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 8*TG**(-1)*U1**(-1)*MG2**3 - 2*TG**(-1)*U1*MG2 - 4*
     +    TG**(-1)*MS2*MG2 - 4*TG**(-1)*MG2**2 - 4*S**(-1)*MS2*MG2 - 4*
     +    S**(-1)*MG2**2 + 4*U1**(-1)*MS2*MG2 - 4*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(99)*N*CF**2 * ( 16*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 48*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MS2*MG2 + 16*S**(-1)*T1*(S+U1)**(-1)*MG2**2 - 16
     +    *S**(-1)*T1**2*U1**(-1)*MG2 + 16*S**(-1)*T1**2*(S+U1)**(-1)*
     +    MG2 + 32*S**(-1)*U1**(-1)*MS2*MG2**2 - 32*S**(-1)*U1**(-1)*
     +    MG2**3 + 16*S**(-1)*MS2*MG2 - 16*S**(-1)*MG2**2 - 8*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 16*T1*U1**(-1)*(S+U1)**(-1)*MS2*
     +    MG2 + 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 + 8*T1*U1**(-1)*MG2
     +     - 8*T1*(S+U1)**(-1)*MG2 + 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2
     +     - 16*U1**(-1)*MS2*MG2 - 16*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(99)*N**2*CF * ( 4*S**(-2)*T1*MS2*MG2 + 4*
     +    S**(-2)*T1*MG2**2 - 2*S**(-2)*U1*MS2*MG2 + 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2**2*MG2 + 4*S**(-2)*MG2**3 - 12*S**(-1)
     +    *T1*U1**(-1)*MS2*MG2 + 28*S**(-1)*T1*U1**(-1)*MG2**2 + 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 4*S**(-1)*T1*(S+U1)**(-1)*
     +    MG2**2 + 8*S**(-1)*T1*MG2 + 12*S**(-1)*T1**2*U1**(-1)*MG2 - 4
     +    *S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 16*S**(-1)*U1**(-1)*MG2**3 + 2*S**(-1)*U1*MG2 - 6*
     +    S**(-1)*MS2*MG2 + 14*S**(-1)*MG2**2 + 2*S*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2 - 2*S*U1**(-2)*MS2*MG2 + 2*S*U1**(-2)*MG2**2
     +     + 2*S*U1**(-1)*MG2 - 4*T1*U1**(-2)*MS2*MG2 + 4*T1*U1**(-2)*
     +    MG2**2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2**2 + 6*T1*U1**(-1)*MG2 + 2*T1*(S+U1)**(-1)*
     +    MG2 - 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-2)*MS2*
     +    MG2**2 + 4*U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 - 6*
     +    U1**(-1)*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(99)*N**2*CF * ( 14*U1**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(101)*N*CF**2 * (  - 8*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(101)*CF**2 * ( 32*TG**(-1)*MS2*MG2**2 - 32
     +    *TG**(-1)*MS2**2*MG2 + 16*MS2*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*N*CF*S4G*S4G2**(-1) * ( 4*TG**(-1)*S*
     +    MS2*MG2 - 4*TG**(-1)*S*MG2**2 + 2*TG**(-1)*U1*MS2*MG2 - 2*
     +    TG**(-1)*U1*MG2**2 - 8*TG**(-1)*MS2**2*MG2 + 8*TG**(-1)*
     +    MG2**3 - 4*S**(-1)*T1*MS2*MG2 + 4*S**(-1)*T1*MG2**2 - 2*
     +    S**(-1)*U1*MS2*MG2 + 2*S**(-1)*U1*MG2**2 - 8*S**(-1)*MS2*
     +    MG2**2 + 4*S**(-1)*MS2**2*MG2 + 4*S**(-1)*MG2**3 - 2*S*MG2 + 
     +    16*T1*U1**(-1)*MG2**2 + 6*T1*MG2 + 8*T1**2*U1**(-1)*MG2 - 8*
     +    U1**(-1)*MS2**2*MG2 + 8*U1**(-1)*MG2**3 - 2*MS2*MG2 + 10*
     +    MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*N*CF * (  - 16*TG**(-1)*S**(-1)*MS2*
     +    MG2**2 + 16*TG**(-1)*S**(-1)*MS2**2*MG2 - 2*TG**(-1)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 2*TG**(-1)*S*(S+U1)**(-1)*MG2**2 + 16*
     +    TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 - 8*TG**(-1)*(S+U1)**(-1)*
     +    MS2**2*MG2 - 8*TG**(-1)*(S+U1)**(-1)*MG2**3 - 2*TG**(-1)*MS2*
     +    MG2 + 2*TG**(-1)*MG2**2 - 8*S**(-1)*MS2*MG2 + 2*S*
     +    (S+U1)**(-1)*MG2 - 4*T1*(S+U1)**(-1)*MG2 + 8*(S+U1)**(-1)*MS2
     +    *MG2 - 8*(S+U1)**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*N*CF**2 * (  - 32*S**(-1)*T1*U1**(-1)
     +    *MG2**2 - 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MG2**2 + 8*S**(-1)*T1*MG2 - 16*S**(-1)*T1**2*
     +    U1**(-1)*MG2 + 16*S**(-1)*T1**2*(S+U1)**(-1)*MG2 + 16*S**(-1)
     +    *U1**(-1)*MS2**2*MG2 - 16*S**(-1)*U1**(-1)*MG2**3 + 4*S**(-1)
     +    *U1*MG2 - 8*S**(-1)*MS2*MG2 + 8*S**(-1)*MG2**2 - 8*T1*
     +    (S+U1)**(-1)*MG2 - 4*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*N**2*CF * ( 4*S**(-2)*T1*MS2*MG2 + 4*
     +    S**(-2)*T1*MG2**2 - 2*S**(-2)*U1*MS2*MG2 + 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2**2*MG2 + 4*S**(-2)*MG2**3 + 16*S**(-1)
     +    *T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MG2**2 + 2*S**(-1)*T1*MG2 + 8*S**(-1)
     +    *T1**2*U1**(-1)*MG2 - 4*S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 8*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MG2**3 + 4*
     +    S**(-1)*MS2*MG2 + 4*S**(-1)*MG2**2 + 2*T1*(S+U1)**(-1)*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*CF**2*S4G*S4G2**(-1) * ( 8*S**(-1)*T1
     +    *MS2*MG2 - 8*S**(-1)*T1*MG2**2 + 4*S**(-1)*U1*MS2*MG2 - 4*
     +    S**(-1)*U1*MG2**2 + 16*S**(-1)*MS2*MG2**2 - 8*S**(-1)*MS2**2*
     +    MG2 - 8*S**(-1)*MG2**3 - 16*T1*U1**(-1)*MG2**2 - 4*T1*MG2 - 8
     +    *T1**2*U1**(-1)*MG2 + 8*U1**(-1)*MS2**2*MG2 - 8*U1**(-1)*
     +    MG2**3 + 4*MS2*MG2 - 4*MG2**2 )
     +
      MGQLLH = MGQLLH + ANG4(103)*CF**2 * ( 4*TG**(-1)*S**(-1)*U1*MS2*
     +    MG2 - 4*TG**(-1)*S**(-1)*U1*MG2**2 + 32*TG**(-1)*S**(-1)*MS2*
     +    MG2**2 - 16*TG**(-1)*S**(-1)*MS2**2*MG2 - 16*TG**(-1)*S**(-1)
     +    *MG2**3 + 4*TG**(-1)*S*(S+U1)**(-1)*MS2*MG2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MG2**2 - 32*TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 + 
     +    16*TG**(-1)*(S+U1)**(-1)*MS2**2*MG2 + 16*TG**(-1)*
     +    (S+U1)**(-1)*MG2**3 - 36*TG**(-1)*MS2*MG2 + 20*TG**(-1)*
     +    MG2**2 + 8*S**(-1)*MS2*MG2 - 8*S**(-1)*MG2**2 - 4*S*
     +    (S+U1)**(-1)*MG2 + 8*T1*(S+U1)**(-1)*MG2 - 16*(S+U1)**(-1)*
     +    MS2*MG2 + 16*(S+U1)**(-1)*MG2**2 + 12*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(104)*N*CF * (  - 2*TG**(-1)*S**(-1)*U1*MS2
     +    *MG2 + 2*TG**(-1)*S**(-1)*U1*MG2**2 - 8*TG**(-1)*S**(-1)*
     +    MS2**2*MG2 + 8*TG**(-1)*S**(-1)*MG2**3 - 16*TG**(-1)*S*
     +    U1**(-1)*MS2*MG2 + 32*TG**(-1)*S*U1**(-1)*MG2**2 + 2*TG**(-1)
     +    *S*(S+U1)**(-1)*MS2*MG2 - 2*TG**(-1)*S*(S+U1)**(-1)*MG2**2 + 
     +    8*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2*U1**(-1)*MG2 - 32*TG**(-1)
     +    *U1**(-1)*MS2*MG2**2 + 32*TG**(-1)*U1**(-1)*MG2**3 + 2*
     +    TG**(-1)*U1*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 - 8*
     +    TG**(-1)*(S+U1)**(-1)*MS2**2*MG2 - 8*TG**(-1)*(S+U1)**(-1)*
     +    MG2**3 - 10*TG**(-1)*MS2*MG2 + 18*TG**(-1)*MG2**2 + 4*S**(-1)
     +    *MS2*MG2 + 4*S**(-1)*MG2**2 + 16*S*U1**(-1)*MG2 - 2*S*
     +    (S+U1)**(-1)*MG2 + 8*T1*U1**(-1)*MG2 - 4*T1*(S+U1)**(-1)*MG2
     +     - 8*U1**(-1)*MS2*MG2 + 24*U1**(-1)*MG2**2 + 8*(S+U1)**(-1)*
     +    MS2*MG2 - 8*(S+U1)**(-1)*MG2**2 + 6*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(104)*N*CF**2 * ( 32*S**(-1)*T1*U1**(-1)*
     +    MG2**2 + 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MG2**2 + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 16*
     +    S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2**2*
     +    MG2 + 16*S**(-1)*U1**(-1)*MG2**3 - 8*S**(-1)*MS2*MG2 + 8*
     +    S**(-1)*MG2**2 + 16*S*U1**(-1)*MG2 + 32*T1*U1**(-1)*MG2 - 8*
     +    T1*(S+U1)**(-1)*MG2 + 32*U1**(-1)*MG2**2 + 8*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(104)*N**2*CF * (  - 4*S**(-2)*T1*MS2*MG2
     +     - 4*S**(-2)*T1*MG2**2 + 2*S**(-2)*U1*MS2*MG2 - 2*S**(-2)*U1*
     +    MG2**2 + 4*S**(-2)*MS2**2*MG2 - 4*S**(-2)*MG2**3 - 16*S**(-1)
     +    *T1*U1**(-1)*MG2**2 - 4*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MG2**2 - 6*S**(-1)*T1*MG2 - 8*S**(-1)
     +    *T1**2*U1**(-1)*MG2 + 4*S**(-1)*T1**2*(S+U1)**(-1)*MG2 + 8*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MG2**3 - 2*
     +    S**(-1)*U1*MG2 + 4*S**(-1)*MS2*MG2 - 12*S**(-1)*MG2**2 - 8*S*
     +    U1**(-1)*MG2 - 16*T1*U1**(-1)*MG2 + 2*T1*(S+U1)**(-1)*MG2 - 
     +    16*U1**(-1)*MG2**2 - 8*MG2 )
     +
      MGQLLH = MGQLLH + ANG4(104)*CF**2 * ( 16*TG**(-1)*S*U1**(-1)*MS2*
     +    MG2 - 32*TG**(-1)*S*U1**(-1)*MG2**2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 4*TG**(-1)*S*(S+U1)**(-1)*MG2**2 - 4*
     +    TG**(-1)*S*MG2 - 8*TG**(-1)*S**2*U1**(-1)*MG2 + 32*TG**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 32*TG**(-1)*U1**(-1)*MG2**3 - 32*
     +    TG**(-1)*(S+U1)**(-1)*MS2*MG2**2 + 16*TG**(-1)*(S+U1)**(-1)*
     +    MS2**2*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MG2**3 + 4*TG**(-1)*MS2
     +    *MG2 - 4*TG**(-1)*MG2**2 - 16*S*U1**(-1)*MG2 + 4*S*
     +    (S+U1)**(-1)*MG2 - 8*T1*U1**(-1)*MG2 + 8*T1*(S+U1)**(-1)*MG2
     +     + 8*U1**(-1)*MS2*MG2 - 24*U1**(-1)*MG2**2 - 16*(S+U1)**(-1)*
     +    MS2*MG2 + 16*(S+U1)**(-1)*MG2**2 )
     +
      MGQLLH = MGQLLH + COLO1(9)*N*CF**2*(S4+MS2) * ( 16*TG**(-2)*S*T1*
     +    S4**(-1)*(S+U1)**(-2)*MG2 + 16*TG**(-2)*S*T1**2*S4**(-1)*
     +    (S+U1)**(-3)*MG2 + 8*TG**(-2)*S*S4**(-1)*(S+U1)**(-1)*MG2 + 
     +    16*S*T1*U1**2*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 + 16*S*T1**2*U1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*S4**(-1)*(S+U1)**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*
     +    U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*
     +    S**3*T1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 
     +    16*S**3*T1**2*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 + 8*S**3*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 )
     +
      MGQLLH = MGQLLH + COLO1(9)*CF**2*(S4+MS2) * ( 32*TG**(-1)*S*T1*U1
     +    *S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*
     +    TG**(-1)*S*T1**2*U1*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S*U1*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*TG**(-1)*S**2*
     +    T1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-1)*MG2 + 32*
     +    TG**(-1)*S**2*T1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*TG**(-1)*S**2*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 )

      MGQLRH = 0D0
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2)**(-1) * ( 4*M2*S**(-1)*T1*MS2*
     +    MG2**(-2) - 48*M2*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) + 32*M2*
     +    S**(-1)*U1**(-1)*S4*MS2**2*MG2**(-2) - 24*M2*S**(-1)*U1**(-1)
     +    *S4**2*MS2*MG2**(-2) - 24*M2*S**(-1)*U1**(-1)*MS2 + 32*M2*
     +    S**(-1)*U1**(-1)*MS2**2*MG2**(-1) - 8*M2*S**(-1)*U1**(-1)*
     +    MS2**3*MG2**(-2) + 8*M2*S**(-1)*S4*MS2*MG2**(-2) + 8*M2*
     +    S**(-1)*MS2*MG2**(-1) - 4*M2*S**(-1)*MS2**2*MG2**(-2) + 16*M2
     +    *U1**(-2)*S4*MS2**2*MG2**(-2) + 16*M2*U1**(-2)*MS2**2*
     +    MG2**(-1) - 8*M2*U1**(-2)*MS2**3*MG2**(-2) - 8*M2*U1**(-1)*
     +    MS2**2*MG2**(-2) + 24*M2**2*S**(-1)*U1**(-1)*S4*MS2*MG2**(-2)
     +     + 24*M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-1) - 16*M2**2*S**(-1)
     +    *U1**(-1)*MS2**2*MG2**(-2) - 4*M2**2*S**(-1)*MS2*MG2**(-2) - 
     +    8*M2**2*U1**(-2)*MS2**2*MG2**(-2) - 8*M2**3*S**(-1)*U1**(-1)*
     +    MS2*MG2**(-2) - 4*S**(-1)*T1*S4*MS2*MG2**(-2) - 4*S**(-1)*T1*
     +    MS2*MG2**(-1) + 24*S**(-1)*U1**(-1)*S4*MS2 - 32*S**(-1)*
     +    U1**(-1)*S4*MS2**2*MG2**(-1) )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2)**(-1) * ( 8*S**(-1)*U1**(-1)*
     +    S4*MS2**3*MG2**(-2) + 24*S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-1)
     +     - 16*S**(-1)*U1**(-1)*S4**2*MS2**2*MG2**(-2) + 8*S**(-1)*
     +    U1**(-1)*S4**3*MS2*MG2**(-2) + 8*S**(-1)*U1**(-1)*MS2*MG2 - 
     +    16*S**(-1)*U1**(-1)*MS2**2 + 8*S**(-1)*U1**(-1)*MS2**3*
     +    MG2**(-1) - 8*S**(-1)*S4*MS2*MG2**(-1) + 4*S**(-1)*S4*MS2**2*
     +    MG2**(-2) - 4*S**(-1)*S4**2*MS2*MG2**(-2) - 4*S**(-1)*MS2 + 4
     +    *S**(-1)*MS2**2*MG2**(-1) - 16*U1**(-2)*S4*MS2**2*MG2**(-1)
     +     + 8*U1**(-2)*S4*MS2**3*MG2**(-2) - 8*U1**(-2)*S4**2*MS2**2*
     +    MG2**(-2) - 8*U1**(-2)*MS2**2 + 8*U1**(-2)*MS2**3*MG2**(-1)
     +     + 8*U1**(-1)*S4*MS2**2*MG2**(-2) + 8*U1**(-1)*MS2**2*
     +    MG2**(-1) )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2)*S4G*S4G2**(-1) * ( 16*TG**(-1)
     +    *S**(-1)*T1*U1**3*S4**(-1)*(S+U1)**(-2) + 16*TG**(-1)*S**(-1)
     +    *U1**2*S4**(-1)*(S+U1)**(-1)*MS2 - 16*TG**(-1)*S**(-1)*U1**2*
     +    S4**(-1)*(S+U1)**(-1)*MG2 + 16*TG**(-1)*S*T1*U1*S4**(-1)*
     +    (S+U1)**(-2) + 32*TG**(-1)*T1*U1**2*S4**(-1)*(S+U1)**(-2) + 
     +    16*TG**(-1)*U1*S4**(-1)*(S+U1)**(-1)*MS2 - 16*TG**(-1)*U1*
     +    S4**(-1)*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**2*S4**(-1)*
     +    (S+U1)**(-1) - 16*U1*S4**(-1)*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2)*SYMBU * ( 16*S**(-1)*T1*U1*
     +    S4**(-1) - 16*S**(-1)*U1*S4**(-1)*MS2 + 16*S**(-1)*U1*
     +    S4**(-1)*MG2 + 16*S**(-1)*U1**2*S4**(-1) + 16*U1*S4**(-1) - 
     +    16*S4**(-1)*MS2 + 16*S4**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2) * ( 16*TG**(-2)*S*T1*S4**(-1)*
     +    (S+U1)**(-2)*MS2 + 16*TG**(-2)*S*T1**2*S4**(-1)*(S+U1)**(-3)*
     +    MS2 - 8*TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-1) - 16*TG**(-2)*
     +    T1**2*U1*S4**(-1)*(S+U1)**(-2) - 16*TG**(-2)*T1**3*U1*
     +    S4**(-1)*(S+U1)**(-3) - 8*TG**(-2)*U1*S4**(-1)*(S+U1)**(-1)*
     +    MS2 + 8*TG**(-2)*U1*S4**(-1)*(S+U1)**(-1)*MG2 + 8*TG**(-1)*U1
     +    *S4**(-1)*(S+U1)**(-1) - 16*S**(-1)*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) + 16*S**(-1)*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 - 16*S**(-1)*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 16*S*T1*U1**2*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*S*T1*U1**2*
     +    S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 16*S*T1**2*
     +    U1**2*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 32*
     +    S*T1**2*U1**2*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)
     +     - 32*S*T1**3*U1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2) )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2) * ( 8*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*S4**(-1)*(S+U1)**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*S*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) + 32*S**2*T1*U1*S4**(-1)*(S+U1)**(-2)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2) + 32*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*S**2*
     +    T1**2*U1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) - 16*
     +    S**2*T1**3*U1*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)
     +     + 16*S**2*U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MS2 - 8*S**2*S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**3*
     +    T1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 16*
     +    S**3*T1**2*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MS2
     +     + 8*S**3*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2
     +     - 8*T1*U1**3*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2) )
     +
      MGQLRH = MGQLRH + N*CF**2*(S4+MS2) * (  - 16*T1**2*U1**3*S4**(-1)
     +    *(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) - 16*T1**3*U1**3*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2) + 8*U1*S4**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-1) + 8*U1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 8*U1**2*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 - 16*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + N*CF**2*S4G*S4G2**(-1) * (  - 16 - 16*TG**(-1)*
     +    S**(-1)*U1*MS2 + 16*TG**(-1)*S**(-1)*U1*MG2 - 16*TG**(-1)*
     +    S**(-1)*MS2*MG2 + 16*TG**(-1)*S**(-1)*MG2**2 + 16*TG**(-1)*S*
     +    U1**(-1)*MS2 - 32*TG**(-1)*S*U1**(-1)*MG2 - 32*TG**(-1)*S - 
     +    16*TG**(-1)*S**2*U1**(-1) + 16*TG**(-1)*U1**(-1)*MS2*MG2 - 16
     +    *TG**(-1)*U1**(-1)*MG2**2 - 16*TG**(-1)*U1 - 16*TG**(-1)*MG2
     +     + 12*S**(-1)*T1 - 48*S**(-1)*U1**(-1)*MS2*MG2 + 24*S**(-1)*
     +    U1**(-1)*MS2**2 + 24*S**(-1)*U1**(-1)*MG2**2 + 16*S**(-1)*U1
     +     + 8*S**(-1)*MS2 + 8*S**(-1)*MG2 - 32*S*U1**(-1) - 16*T1*
     +    U1**(-1) - 16*U1**(-2)*MS2*MG2 + 16*U1**(-2)*MS2**2 + 16*
     +    U1**(-1)*MS2 - 32*U1**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + N*CF**2 * ( 16*M2*S**(-1)*U1**(-1)*S4*MS2*
     +    MG2**(-2) + 16*M2*S**(-1)*U1**(-1)*MS2*MG2**(-1) - 16*M2*
     +    S**(-1)*U1**(-1)*MS2**2*MG2**(-2) - 8*M2*S**(-1)*U1**(-1) - 4
     +    *M2*S**(-1)*MS2*MG2**(-2) - 8*M2*U1**(-2)*MS2**2*MG2**(-2) - 
     +    8*M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-2) + 4*S**(-1)*T1*MS2*
     +    MG2**(-2) - 16*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) + 16*S**(-1)
     +    *U1**(-1)*S4*MS2**2*MG2**(-2) + 8*S**(-1)*U1**(-1)*S4 - 8*
     +    S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-2) - 32*S**(-1)*U1**(-1)*
     +    MS2 + 16*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) - 8*S**(-1)*
     +    U1**(-1)*MS2**3*MG2**(-2) + 24*S**(-1)*U1**(-1)*MG2 + 4*
     +    S**(-1)*S4*MS2*MG2**(-2) + 4*S**(-1)*MS2*MG2**(-1) - 4*
     +    S**(-1)*MS2**2*MG2**(-2) - 4*S**(-1) + 8*U1**(-2)*S4*MS2**2*
     +    MG2**(-2) - 8*U1**(-2)*MS2 + 8*U1**(-2)*MS2**2*MG2**(-1) - 8*
     +    U1**(-2)*MS2**3*MG2**(-2) - 8*U1**(-1)*MS2**2*MG2**(-2) - 8*
     +    U1**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1)*S4G*S4G2**(-1) * ( 8*
     +    TG**(-2)*S*T1*MS2 + 4*TG**(-2)*S*U1*MS2 - 4*TG**(-2)*S*S4*MS2
     +     + 4*TG**(-2)*S**2*MS2 + 4*TG**(-2)*T1*U1*MS2 + 4*TG**(-2)*T1
     +    *MS2*MG2 - 4*TG**(-2)*T1*MS2**2 - 4*TG**(-2)*S4*MS2*MG2 + 4*
     +    TG**(-2)*S4*MS2**2 + 2*TG**(-1)*S**(-1)*T1*U1*MS2 + 4*
     +    TG**(-1)*S**(-1)*T1*MS2*MG2 - 4*TG**(-1)*S**(-1)*T1*MS2**2 - 
     +    2*TG**(-1)*S**(-1)*S4*MS2*MG2 + 2*TG**(-1)*S**(-1)*S4*MS2**2
     +     - 2*TG**(-1)*S*T1*U1**(-1)*MS2 - 2*TG**(-1)*S*U1**(-1)*MS2*
     +    MG2 + 2*TG**(-1)*S*U1**(-1)*MS2**2 + 4*TG**(-1)*S*MS2 - 2*
     +    TG**(-1)*S**2*U1**(-1)*MS2 - 2*TG**(-1)*T1*U1**(-1)*MS2*MG2
     +     + 2*TG**(-1)*T1*U1**(-1)*MS2**2 + 10*TG**(-1)*T1*MS2 + 2*
     +    TG**(-1)*U1*MS2 - 6*TG**(-1)*S4*MS2 + 2*TG**(-1)*MS2*MG2 - 2*
     +    TG**(-1)*MS2**2 + 4*S**(-1)*T1*MS2 - 2*S**(-1)*S4*MS2 - 2*S*
     +    U1**(-1)*MS2 - 2*T1*U1**(-1)*MS2 + 2*MS2 )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * ( 4*M2*TG**(-2)*S*MS2*
     +    MG2**(-1) + 4*M2*TG**(-2)*MS2 - 4*M2*TG**(-2)*MS2**2*
     +    MG2**(-1) + 2*M2*TG**(-1)*S**(-2)*MS2*MG2 + 2*M2*TG**(-1)*
     +    S**(-2)*MS2*MG2**2*(S4G-S)**(-1) - 2*M2*TG**(-1)*S**(-2)*
     +    MS2**2*MG2*(S4G-S)**(-1) - 2*M2*TG**(-1)*S**(-2)*MS2**2 + 4*
     +    M2*TG**(-1)*S**(-1)*MS2*MG2*(S4G-S)**(-1) + 4*M2*TG**(-1)*
     +    S**(-1)*MS2 - 2*M2*TG**(-1)*S**(-1)*MS2**2*MG2**(-1) - 2*M2*
     +    TG**(-1)*S**(-1)*MS2**2*(S4G-S)**(-1) + 6*M2*TG**(-1)*MS2*
     +    MG2**(-1) + 2*M2*TG**(-1)*MS2*(S4G-S)**(-1) + 4*M2*S**(-2)*
     +    U1**(-1)*MS2*MG2 + 4*M2*S**(-2)*U1**(-1)*MS2*MG2**2*
     +    (S4G-S)**(-1) + 4*M2*S**(-2)*U1**(-1)*MS2**2*MG2*
     +    (S4G-S)**(-1) + 4*M2*S**(-2)*U1**(-1)*MS2**2 + 24*M2*S**(-1)*
     +    U1**(-1)*S4*MS2*MG2**(-1) - 16*M2*S**(-1)*U1**(-1)*S4*MS2**2*
     +    MG2**(-2) + 12*M2*S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-2) + 6*M2
     +    *S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) + 14*M2*S**(-1)*
     +    U1**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * (  - 16*M2*S**(-1)*
     +    U1**(-1)*MS2**2*MG2**(-1) + 4*M2*S**(-1)*U1**(-1)*MS2**2*
     +    (S4G-S)**(-1) + 4*M2*S**(-1)*U1**(-1)*MS2**3*MG2**(-2) - 8*M2
     +    *S**(-1)*S4*MS2*MG2**(-2) - 6*M2*S**(-1)*MS2*MG2**(-1) + 4*M2
     +    *S**(-1)*MS2**2*MG2**(-2) + 2*M2*S*U1**(-1)*MS2*MG2**(-2) + 8
     +    *M2*U1**(-1)*S4*MS2*MG2**(-2) + 8*M2*U1**(-1)*MS2*MG2**(-1)
     +     + 2*M2*U1**(-1)*MS2*(S4G-S)**(-1) - 4*M2*U1**(-1)*MS2**2*
     +    MG2**(-2) - 2*M2*MS2*MG2**(-2) - 12*M2**2*S**(-1)*U1**(-1)*S4
     +    *MS2*MG2**(-2) - 12*M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-1) + 8*
     +    M2**2*S**(-1)*U1**(-1)*MS2**2*MG2**(-2) + 4*M2**2*S**(-1)*MS2
     +    *MG2**(-2) - 4*M2**2*U1**(-1)*MS2*MG2**(-2) + 4*M2**3*S**(-1)
     +    *U1**(-1)*MS2*MG2**(-2) + 4*TG**(-2)*S*T1*MS2*MG2**(-1) - 4*
     +    TG**(-2)*S*S4*MS2*MG2**(-1) - 12*TG**(-2)*S*MS2 + 4*TG**(-2)*
     +    S*MS2**2*MG2**(-1) - 4*TG**(-2)*T1*MS2 - 4*TG**(-2)*T1*MS2**2
     +    *MG2**(-1) - 4*TG**(-2)*S4*MS2 + 4*TG**(-2)*S4*MS2**2*
     +    MG2**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * (  - 4*TG**(-2)*MS2*
     +    MG2 + 8*TG**(-2)*MS2**2 - 4*TG**(-2)*MS2**3*MG2**(-1) - 2*
     +    TG**(-1)*S**(-2)*T1*MS2*MG2 - 2*TG**(-1)*S**(-2)*T1*MS2*
     +    MG2**2*(S4G-S)**(-1) - 6*TG**(-1)*S**(-2)*T1*MS2**2*MG2*
     +    (S4G-S)**(-1) - 6*TG**(-1)*S**(-2)*T1*MS2**2 - 2*TG**(-1)*
     +    S**(-2)*S4*MS2*MG2 - 2*TG**(-1)*S**(-2)*S4*MS2*MG2**2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*S4*MS2**2*MG2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*S4*MS2**2 - 2*TG**(-1)*
     +    S**(-2)*MS2*MG2**2 - 2*TG**(-1)*S**(-2)*MS2*MG2**3*
     +    (S4G-S)**(-1) + 4*TG**(-1)*S**(-2)*MS2**2*MG2 + 4*TG**(-1)*
     +    S**(-2)*MS2**2*MG2**2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*
     +    MS2**3*MG2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*MS2**3 - 2*
     +    TG**(-1)*S**(-1)*T1*MS2*MG2*(S4G-S)**(-1) - 2*TG**(-1)*
     +    S**(-1)*T1*MS2**2*MG2**(-1) - 6*TG**(-1)*S**(-1)*T1*MS2**2*
     +    (S4G-S)**(-1) - 4*TG**(-1)*S**(-1)*S4*MS2*MG2*(S4G-S)**(-1)
     +     - 4*TG**(-1)*S**(-1)*S4*MS2 )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * ( 2*TG**(-1)*S**(-1)*
     +    S4*MS2**2*MG2**(-1) + 2*TG**(-1)*S**(-1)*S4*MS2**2*
     +    (S4G-S)**(-1) - 8*TG**(-1)*S**(-1)*MS2*MG2 - 8*TG**(-1)*
     +    S**(-1)*MS2*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-1)*MS2**2*
     +    MG2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-1)*MS2**2 - 2*TG**(-1)*
     +    S**(-1)*MS2**3*MG2**(-1) - 2*TG**(-1)*S**(-1)*MS2**3*
     +    (S4G-S)**(-1) - 4*TG**(-1)*S*MS2*MG2**(-1) - 2*TG**(-1)*S*MS2
     +    *(S4G-S)**(-1) + 2*TG**(-1)*T1*MS2*MG2**(-1) - 6*TG**(-1)*S4*
     +    MS2*MG2**(-1) - 2*TG**(-1)*S4*MS2*(S4G-S)**(-1) - 8*TG**(-1)*
     +    MS2*MG2*(S4G-S)**(-1) - 10*TG**(-1)*MS2 + 6*TG**(-1)*MS2**2*
     +    MG2**(-1) - 2*TG**(-1)*MS2**2*(S4G-S)**(-1) - 4*S**(-2)*
     +    U1**(-1)*S4*MS2*MG2 - 4*S**(-2)*U1**(-1)*S4*MS2*MG2**2*
     +    (S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*S4*MS2**2*MG2*
     +    (S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*S4*MS2**2 - 4*S**(-2)*
     +    U1**(-1)*MS2*MG2**2 - 4*S**(-2)*U1**(-1)*MS2*MG2**3*
     +    (S4G-S)**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * ( 4*S**(-2)*U1**(-1)*
     +    MS2**3*MG2*(S4G-S)**(-1) + 4*S**(-2)*U1**(-1)*MS2**3 + 2*
     +    S**(-2)*MS2*MG2 + 2*S**(-2)*MS2*MG2**2*(S4G-S)**(-1) + 6*
     +    S**(-2)*MS2**2*MG2*(S4G-S)**(-1) + 6*S**(-2)*MS2**2 - 6*
     +    S**(-1)*U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) - 14*S**(-1)*
     +    U1**(-1)*S4*MS2 + 16*S**(-1)*U1**(-1)*S4*MS2**2*MG2**(-1) - 4
     +    *S**(-1)*U1**(-1)*S4*MS2**2*(S4G-S)**(-1) - 4*S**(-1)*
     +    U1**(-1)*S4*MS2**3*MG2**(-2) - 12*S**(-1)*U1**(-1)*S4**2*MS2*
     +    MG2**(-1) + 8*S**(-1)*U1**(-1)*S4**2*MS2**2*MG2**(-2) - 4*
     +    S**(-1)*U1**(-1)*S4**3*MS2*MG2**(-2) - 8*S**(-1)*U1**(-1)*MS2
     +    *MG2 - 8*S**(-1)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) + 4*
     +    S**(-1)*U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) + 12*S**(-1)*
     +    U1**(-1)*MS2**2 - 4*S**(-1)*U1**(-1)*MS2**3*MG2**(-1) + 4*
     +    S**(-1)*U1**(-1)*MS2**3*(S4G-S)**(-1) + 6*S**(-1)*S4*MS2*
     +    MG2**(-1) - 4*S**(-1)*S4*MS2**2*MG2**(-2) + 4*S**(-1)*S4**2*
     +    MS2*MG2**(-2) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)**(-1) * ( 2*S**(-1)*MS2*MG2*
     +    (S4G-S)**(-1) + 2*S**(-1)*MS2 - 2*S**(-1)*MS2**2*MG2**(-1) + 
     +    6*S**(-1)*MS2**2*(S4G-S)**(-1) - 2*S*U1**(-1)*S4*MS2*
     +    MG2**(-2) - 2*S*U1**(-1)*MS2*MG2**(-1) - 2*S*U1**(-1)*MS2*
     +    (S4G-S)**(-1) - 8*U1**(-1)*S4*MS2*MG2**(-1) - 2*U1**(-1)*S4*
     +    MS2*(S4G-S)**(-1) + 4*U1**(-1)*S4*MS2**2*MG2**(-2) - 4*
     +    U1**(-1)*S4**2*MS2*MG2**(-2) - 6*U1**(-1)*MS2*MG2*
     +    (S4G-S)**(-1) - 6*U1**(-1)*MS2 + 4*U1**(-1)*MS2**2*MG2**(-1)
     +     + 4*U1**(-1)*MS2**2*(S4G-S)**(-1) + 2*S4*MS2*MG2**(-2) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)*S4G*S4G2**(-1) * (  - 8*
     +    TG**(-2)*S*T1*U1**2*S4**(-1)*(S+U1)**(-2) - 4*TG**(-2)*S*U1*
     +    S4**(-1)*(S+U1)**(-1)*MS2 + 4*TG**(-2)*S*U1*S4**(-1)*
     +    (S+U1)**(-1)*MG2 - 4*TG**(-2)*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2) - 4*TG**(-2)*T1*U1**3*S4**(-1)*(S+U1)**(-2) - 4*
     +    TG**(-2)*U1**2*S4**(-1)*(S+U1)**(-1)*MS2 + 4*TG**(-2)*U1**2*
     +    S4**(-1)*(S+U1)**(-1)*MG2 - 8*TG**(-1)*S**(-1)*T1*U1**3*
     +    S4**(-1)*(S+U1)**(-2) - 8*TG**(-1)*S**(-1)*U1**2*S4**(-1)*
     +    (S+U1)**(-1)*MS2 + 8*TG**(-1)*S**(-1)*U1**2*S4**(-1)*
     +    (S+U1)**(-1)*MG2 - 8*TG**(-1)*S*T1*U1*S4**(-1)*(S+U1)**(-2)
     +     + 4*TG**(-1)*S*U1*S4**(-1)*(S+U1)**(-1) - 16*TG**(-1)*T1*
     +    U1**2*S4**(-1)*(S+U1)**(-2) - 8*TG**(-1)*U1*S4**(-1)*
     +    (S+U1)**(-1)*MS2 + 8*TG**(-1)*U1*S4**(-1)*(S+U1)**(-1)*MG2 + 
     +    4*TG**(-1)*U1**2*S4**(-1)*(S+U1)**(-1) + 8*S**(-1)*U1**2*
     +    S4**(-1)*(S+U1)**(-1) + 8*U1*S4**(-1)*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2)*SYMBU * (  - 8*S**(-1)*T1*U1*
     +    S4**(-1) + 8*S**(-1)*U1*S4**(-1)*MS2 - 8*S**(-1)*U1*S4**(-1)*
     +    MG2 - 8*S**(-1)*U1**2*S4**(-1) + 4*S*U1**(-1)*S4**(-1)*MS2 - 
     +    4*S*U1**(-1)*S4**(-1)*MG2 - 4*S*S4**(-1) - 4*T1*S4**(-1) - 12
     +    *U1*S4**(-1) + 12*S4**(-1)*MS2 - 12*S4**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + N**2*CF*(S4+MS2) * ( 4*TG**(-2)*S*T1*U1*
     +    S4**(-1)*(S+U1)**(-2) + 4*TG**(-2)*T1*U1**2*S4**(-1)*
     +    (S+U1)**(-2) + 4*TG**(-2)*U1*S4**(-1)*(S+U1)**(-1)*MS2 - 4*
     +    TG**(-2)*U1*S4**(-1)*(S+U1)**(-1)*MG2 - 4*TG**(-1)*U1*
     +    S4**(-1)*(S+U1)**(-1) + 8*S**(-1)*T1*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 8*S**(-1)*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 + 8*S**(-1)*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 - 4*S*U1**(-1)*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MS2 + 4*S*U1**(-1)*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1)*MG2 + 4*T1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-1) - 12*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)
     +    *MS2 + 12*S4**(-1)*(M2*(S+U1)+T1*U1)**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + N**2*CF*S4G*S4G2**(-1) * ( 4 - 4*TG**(-2)*S*T1*
     +    MS2*MG2**(-1) + 4*TG**(-2)*S*T1 - 4*TG**(-2)*S*S4 + 8*
     +    TG**(-2)*S*MS2 - 4*TG**(-2)*S*MS2**2*MG2**(-1) + 4*TG**(-2)*S
     +    *MG2 + 4*TG**(-2)*T1*U1 - 4*TG**(-2)*T1*MS2 + 4*TG**(-2)*T1*
     +    MS2**2*MG2**(-1) + 4*TG**(-2)*U1*MS2 - 4*TG**(-2)*U1*MG2 + 4*
     +    TG**(-2)*S4*MS2 - 4*TG**(-2)*S4*MG2 + 20*TG**(-2)*MS2*MG2 - 
     +    20*TG**(-2)*MS2**2 + 4*TG**(-2)*MS2**3*MG2**(-1) - 4*TG**(-2)
     +    *MG2**2 + 8*TG**(-1)*S**(-2)*T1*MS2*MG2 + 6*TG**(-1)*S**(-2)*
     +    T1*MS2**2 + 2*TG**(-1)*S**(-2)*T1*MG2**2 - 10*TG**(-1)*
     +    S**(-2)*MS2*MG2**2 - 2*TG**(-1)*S**(-2)*MS2**2*MG2 + 10*
     +    TG**(-1)*S**(-2)*MS2**3 + 2*TG**(-1)*S**(-2)*MG2**3 + 2*
     +    TG**(-1)*S**(-1)*T1*U1 + 6*TG**(-1)*S**(-1)*T1*MS2 + 2*
     +    TG**(-1)*S**(-1)*T1*MS2**2*MG2**(-1) + 4*TG**(-1)*S**(-1)*T1*
     +    MG2 + 8*TG**(-1)*S**(-1)*U1*MS2 - 8*TG**(-1)*S**(-1)*U1*MG2
     +     + 2*TG**(-1)*S**(-1)*S4*MS2 - 2*TG**(-1)*S**(-1)*S4*MG2 + 6*
     +    TG**(-1)*S**(-1)*MS2*MG2 )
     +
      MGQLRH = MGQLRH + N**2*CF*S4G*S4G2**(-1) * (  - 4*TG**(-1)*
     +    S**(-1)*MS2**2 + 2*TG**(-1)*S**(-1)*MS2**3*MG2**(-1) - 4*
     +    TG**(-1)*S**(-1)*MG2**2 - 2*TG**(-1)*S*T1*U1**(-1) - 2*
     +    TG**(-1)*S*U1**(-1)*MS2 + 14*TG**(-1)*S*U1**(-1)*MG2 + 4*
     +    TG**(-1)*S*MS2*MG2**(-1) + 8*TG**(-1)*S + 6*TG**(-1)*S**2*
     +    U1**(-1) + 2*TG**(-1)*T1*U1**(-1)*MS2 - 2*TG**(-1)*T1*
     +    U1**(-1)*MG2 - 2*TG**(-1)*T1*MS2*MG2**(-1) + 4*TG**(-1)*T1 - 
     +    4*TG**(-1)*U1**(-1)*MS2*MG2 - 4*TG**(-1)*U1**(-1)*MS2**2 + 8*
     +    TG**(-1)*U1**(-1)*MG2**2 + 2*TG**(-1)*U1 - 6*TG**(-1)*S4 + 18
     +    *TG**(-1)*MS2 - 6*TG**(-1)*MS2**2*MG2**(-1) + 8*TG**(-1)*MG2
     +     - 4*S**(-2)*U1**(-1)*MS2*MG2**2 - 4*S**(-2)*U1**(-1)*MS2**2*
     +    MG2 + 4*S**(-2)*U1**(-1)*MS2**3 + 4*S**(-2)*U1**(-1)*MG2**3
     +     - 8*S**(-2)*MS2*MG2 - 6*S**(-2)*MS2**2 - 2*S**(-2)*MG2**2 - 
     +    4*S**(-1)*T1 + 20*S**(-1)*U1**(-1)*MS2*MG2 - 16*S**(-1)*
     +    U1**(-1)*MS2**2 - 4*S**(-1)*U1**(-1)*MG2**2 - 8*S**(-1)*U1 - 
     +    2*S**(-1)*S4 )
     +
      MGQLRH = MGQLRH + N**2*CF*S4G*S4G2**(-1) * (  - 10*S**(-1)*MS2 - 
     +    2*S**(-1)*MS2**2*MG2**(-1) - 4*S**(-1)*MG2 + 14*S*U1**(-1) + 
     +    6*T1*U1**(-1) + 4*U1**(-1)*MS2 + 8*U1**(-1)*MG2 + 2*MS2*
     +    MG2**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF * (  - 2*M2*TG**(-1)*S**(-2)*MS2*MG2*
     +    (S4G-S)**(-1) + 2*M2*TG**(-1)*S**(-2)*MG2**2*(S4G-S)**(-1) - 
     +    2*M2*TG**(-1)*S**(-1)*MS2*(S4G-S)**(-1) + 4*M2*TG**(-1)*
     +    S**(-1)*MG2*(S4G-S)**(-1) + 2*M2*TG**(-1)*(S4G-S)**(-1) + 4*
     +    M2*S**(-2)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) + 4*M2*S**(-2)*
     +    U1**(-1)*MG2**2*(S4G-S)**(-1) - 8*M2*S**(-1)*U1**(-1)*S4*MS2*
     +    MG2**(-2) - 8*M2*S**(-1)*U1**(-1)*MS2*MG2**(-1) + 4*M2*
     +    S**(-1)*U1**(-1)*MS2*(S4G-S)**(-1) + 8*M2*S**(-1)*U1**(-1)*
     +    MS2**2*MG2**(-2) + 6*M2*S**(-1)*U1**(-1)*MG2*(S4G-S)**(-1) + 
     +    4*M2*S**(-1)*U1**(-1) + 4*M2*S**(-1)*MS2*MG2**(-2) - 4*M2*
     +    U1**(-1)*MS2*MG2**(-2) + 2*M2*U1**(-1)*(S4G-S)**(-1) + 4*
     +    M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-2) + 4*TG**(-2)*S*MS2*
     +    MG2**(-1) + 4*TG**(-2)*S - 4*TG**(-2)*T1 - 4*TG**(-2)*MS2**2*
     +    MG2**(-1) + 4*TG**(-2)*MG2 - 6*TG**(-1)*S**(-2)*T1*MS2*MG2*
     +    (S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*T1*MG2**2*(S4G-S)**(-1) + 
     +    2*TG**(-1)*S**(-2)*S4*MS2*MG2*(S4G-S)**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF * (  - 2*TG**(-1)*S**(-2)*S4*MG2**2*
     +    (S4G-S)**(-1) + 12*TG**(-1)*S**(-2)*MS2*MG2**2*(S4G-S)**(-1)
     +     - 2*TG**(-1)*S**(-2)*MS2**2*MG2*(S4G-S)**(-1) - 2*TG**(-1)*
     +    S**(-2)*MS2**2 - 8*TG**(-1)*S**(-2)*MS2**3*(S4G-S)**(-1) + 2*
     +    TG**(-1)*S**(-2)*MG2**2 - 2*TG**(-1)*S**(-2)*MG2**3*
     +    (S4G-S)**(-1) - 6*TG**(-1)*S**(-1)*T1*MS2*(S4G-S)**(-1) - 2*
     +    TG**(-1)*S**(-1)*T1*MG2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-1)*S4
     +    *MS2*(S4G-S)**(-1) - 4*TG**(-1)*S**(-1)*S4*MG2*(S4G-S)**(-1)
     +     + 14*TG**(-1)*S**(-1)*MS2*MG2*(S4G-S)**(-1) + 2*TG**(-1)*
     +    S**(-1)*MS2 - 2*TG**(-1)*S**(-1)*MS2**2*MG2**(-1) - 2*
     +    TG**(-1)*S**(-1)*MS2**2*(S4G-S)**(-1) + 4*TG**(-1)*S**(-1)*
     +    MG2 - 4*TG**(-1)*S**(-1)*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*S*
     +    (S4G-S)**(-1) - 2*TG**(-1)*S4*(S4G-S)**(-1) + 6*TG**(-1)*MS2*
     +    MG2**(-1) + 2*TG**(-1)*MS2*(S4G-S)**(-1) + 10*TG**(-1) - 4*
     +    S**(-2)*U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) - 4*S**(-2)*
     +    U1**(-1)*S4*MG2**2*(S4G-S)**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF * ( 8*S**(-2)*U1**(-1)*MS2*MG2 + 8*
     +    S**(-2)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) + 4*S**(-2)*U1**(-1)*MS2**2
     +     - 8*S**(-2)*U1**(-1)*MS2**3*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*MG2**2 - 4*S**(-2)*U1**(-1)*MG2**3*(S4G-S)**(-1) + 6
     +    *S**(-2)*MS2*MG2*(S4G-S)**(-1) + 2*S**(-2)*MG2**2*
     +    (S4G-S)**(-1) + 8*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) - 4*
     +    S**(-1)*U1**(-1)*S4*MS2*(S4G-S)**(-1) - 8*S**(-1)*U1**(-1)*S4
     +    *MS2**2*MG2**(-2) - 6*S**(-1)*U1**(-1)*S4*MG2*(S4G-S)**(-1)
     +     - 4*S**(-1)*U1**(-1)*S4 + 4*S**(-1)*U1**(-1)*S4**2*MS2*
     +    MG2**(-2) + 16*S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) + 22*
     +    S**(-1)*U1**(-1)*MS2 - 8*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 
     +    4*S**(-1)*U1**(-1)*MS2**2*(S4G-S)**(-1) + 4*S**(-1)*U1**(-1)*
     +    MS2**3*MG2**(-2) - 6*S**(-1)*U1**(-1)*MG2 - 4*S**(-1)*
     +    U1**(-1)*MG2**2*(S4G-S)**(-1) - 4*S**(-1)*S4*MS2*MG2**(-2) - 
     +    2*S**(-1)*MS2*MG2**(-1) )
     +
      MGQLRH = MGQLRH + N**2*CF * ( 6*S**(-1)*MS2*(S4G-S)**(-1) + 4*
     +    S**(-1)*MS2**2*MG2**(-2) + 2*S**(-1)*MG2*(S4G-S)**(-1) + 4*
     +    S**(-1) + 2*S*U1**(-1)*MS2*MG2**(-2) + 2*S*U1**(-1)*
     +    (S4G-S)**(-1) + 4*U1**(-1)*S4*MS2*MG2**(-2) - 2*U1**(-1)*S4*
     +    (S4G-S)**(-1) + 4*U1**(-1)*MS2*MG2**(-1) + 8*U1**(-1)*MS2*
     +    (S4G-S)**(-1) - 4*U1**(-1)*MS2**2*MG2**(-2) + 2*U1**(-1)*MG2*
     +    (S4G-S)**(-1) + 2*U1**(-1) - 2*MS2*MG2**(-2) )
     +
      MGQLRH = MGQLRH + ANG4(21)*N*CF**2*S4G*S4G2**(-1) * (  - 32 - 8*
     +    TG**(-1)*S**(-1)*U1*MS2 - 16*TG**(-1)*S**(-1)*MS2*MG2 + 16*
     +    TG**(-1)*S**(-1)*MS2**2 - 16*TG**(-1)*S*U1**(-1)*MS2 - 32*
     +    TG**(-1)*S*U1**(-1)*MG2 - 32*TG**(-1)*S - 16*TG**(-1)*S**2*
     +    U1**(-1) - 32*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*
     +    U1**(-1)*MS2**2 - 16*TG**(-1)*U1**(-1)*MG2**2 - 16*TG**(-1)*
     +    U1 - 24*TG**(-1)*MS2 - 32*TG**(-1)*MG2 - 16*S**(-1)*MS2 - 32*
     +    S*U1**(-1) - 16*T1*U1**(-1) - 48*U1**(-1)*MS2 - 16*U1**(-1)*
     +    MG2 )
     +
      MGQLRH = MGQLRH + ANG4(21)*N*CF**2 * (  - 8*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 8*TG**(-2)*(S+U1)**(-1)*MG2**2 + 4*TG**(-2)*MG2 - 
     +    8*TG**(-1)*(S+U1)**(-1)*MS2 + 16*TG**(-1)*(S+U1)**(-1)*MG2 + 
     +    4*TG**(-1) + 8*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(21)*N**2*CF*S4G*S4G2**(-1) * ( 14 + 2*
     +    TG**(-2)*S*MG2 + 2*TG**(-2)*U1*MG2 + 4*TG**(-2)*MS2*MG2 + 4*
     +    TG**(-2)*MG2**2 + 2*TG**(-1)*S**(-1)*U1*MS2 + 4*TG**(-1)*
     +    S**(-1)*MS2*MG2 - 4*TG**(-1)*S**(-1)*MS2**2 + 8*TG**(-1)*S*
     +    U1**(-1)*MS2 + 16*TG**(-1)*S*U1**(-1)*MG2 + 14*TG**(-1)*S + 8
     +    *TG**(-1)*S**2*U1**(-1) + 16*TG**(-1)*U1**(-1)*MS2*MG2 - 8*
     +    TG**(-1)*U1**(-1)*MS2**2 + 8*TG**(-1)*U1**(-1)*MG2**2 + 6*
     +    TG**(-1)*U1 + 12*TG**(-1)*MS2 + 18*TG**(-1)*MG2 + 4*S**(-1)*
     +    MS2 + 16*S*U1**(-1) + 8*T1*U1**(-1) + 24*U1**(-1)*MS2 + 8*
     +    U1**(-1)*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(21)*N**2*CF * ( 4*TG**(-2)*(S+U1)**(-1)*
     +    MS2*MG2 - 4*TG**(-2)*(S+U1)**(-1)*MG2**2 - 2*TG**(-2)*MG2 + 4
     +    *TG**(-1)*(S+U1)**(-1)*MS2 - 8*TG**(-1)*(S+U1)**(-1)*MG2 - 2*
     +    TG**(-1) - 4*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(22)*N*CF**2 * ( 8*TG**(-2)*MS2*MG2 + 8*
     +    TG**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(23)*N*CF**2*S4G*S4G2**(-1) * (  - 8*
     +    TG**(-1)*S**(-1)*U1 - 16*TG**(-1)*S**(-1)*MG2 - 8*TG**(-1)*S*
     +    U1**(-1) - 16*TG**(-1)*U1**(-1)*MG2 - 16*TG**(-1) - 16*
     +    S**(-1) - 16*U1**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(23)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-1)*
     +    S**(-1)*U1 + 4*TG**(-1)*S**(-1)*MG2 + 4*TG**(-1)*S*U1**(-1)
     +     + 8*TG**(-1)*U1**(-1)*MG2 + 6*TG**(-1) + 4*S**(-1) + 8*
     +    U1**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(43)*N*CF**2 * ( 8*TG**(-2)*S*MS2*MG2 + 8*
     +    TG**(-2)*MS2**2*MG2 - 8*TG**(-2)*MS2**3 + 8*TG**(-1)*S*MS2 + 
     +    8*TG**(-1)*MS2**2 )
     +
      MGQLRH = MGQLRH + ANG4(44)*N*CF**2*S4G*S4G2**(-1) * (  - 32*
     +    TG**(-1)*S*U1**(-1)*MS2*MG2 + 8*TG**(-1)*S*U1**(-1)*MS2**2 - 
     +    8*TG**(-1)*S*U1**(-1)*MG2**2 - 8*TG**(-1)*S*U1 - 24*TG**(-1)*
     +    S*MS2 - 16*TG**(-1)*S*MG2 - 16*TG**(-1)*S**2*U1**(-1)*MS2 - 
     +    16*TG**(-1)*S**2*U1**(-1)*MG2 - 16*TG**(-1)*S**2 - 8*TG**(-1)
     +    *S**3*U1**(-1) - 16*TG**(-1)*U1**(-1)*MS2*MG2**2 + 16*
     +    TG**(-1)*U1**(-1)*MS2**3 - 8*TG**(-1)*U1*MS2 - 8*TG**(-1)*MS2
     +    *MG2 + 8*TG**(-1)*MS2**2 - 8*S*T1*U1**(-1) - 40*S*U1**(-1)*
     +    MS2 - 8*S*U1**(-1)*MG2 - 16*S - 16*S**2*U1**(-1) - 16*T1*
     +    U1**(-1)*MS2 - 16*U1**(-1)*MS2*MG2 - 16*U1**(-1)*MS2**2 - 8*
     +    MS2 )
     +
      MGQLRH = MGQLRH + ANG4(44)*N*CF**2 * (  - 12 - 4*TG**(-2)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 4*TG**(-2)*S*(S+U1)**(-1)*MG2**2 + 4*
     +    TG**(-2)*S*MG2 + 8*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 - 16*
     +    TG**(-2)*(S+U1)**(-1)*MS2**2*MG2 + 8*TG**(-2)*(S+U1)**(-1)*
     +    MS2**3 + 8*TG**(-2)*MS2*MG2 - 8*TG**(-2)*MS2**2 - 4*TG**(-1)*
     +    S*(S+U1)**(-1)*MS2 + 8*TG**(-1)*S*(S+U1)**(-1)*MG2 + 4*
     +    TG**(-1)*S + 16*TG**(-1)*(S+U1)**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    (S+U1)**(-1)*MS2**2 + 8*TG**(-1)*MS2 - 8*S**(-1)*T1 - 4*
     +    S**(-1)*U1 - 8*S**(-1)*MS2 - 8*S*U1**(-1) + 4*S*(S+U1)**(-1)
     +     - 16*T1*U1**(-1) - 16*U1**(-1)*MS2 + 8*(S+U1)**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(44)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-2)*
     +    S*U1*MG2 + 6*TG**(-2)*S*MS2*MG2 + 2*TG**(-2)*S*MG2**2 + 2*
     +    TG**(-2)*S**2*MG2 + 4*TG**(-2)*MS2*MG2**2 - 4*TG**(-2)*MS2**3
     +     + 16*TG**(-1)*S*U1**(-1)*MS2*MG2 - 4*TG**(-1)*S*U1**(-1)*
     +    MS2**2 + 4*TG**(-1)*S*U1**(-1)*MG2**2 + 4*TG**(-1)*S*U1 + 14*
     +    TG**(-1)*S*MS2 + 10*TG**(-1)*S*MG2 + 8*TG**(-1)*S**2*U1**(-1)
     +    *MS2 + 8*TG**(-1)*S**2*U1**(-1)*MG2 + 8*TG**(-1)*S**2 + 4*
     +    TG**(-1)*S**3*U1**(-1) + 8*TG**(-1)*U1**(-1)*MS2*MG2**2 - 8*
     +    TG**(-1)*U1**(-1)*MS2**3 + 2*TG**(-1)*U1*MS2 + 10*TG**(-1)*
     +    MS2*MG2 - 2*TG**(-1)*MS2**2 + 4*S*T1*U1**(-1) + 20*S*U1**(-1)
     +    *MS2 + 4*S*U1**(-1)*MG2 + 8*S + 8*S**2*U1**(-1) + 8*T1*
     +    U1**(-1)*MS2 + 8*U1**(-1)*MS2*MG2 + 8*U1**(-1)*MS2**2 + 6*MS2
     +     )
     +
      MGQLRH = MGQLRH + ANG4(44)*N**2*CF * ( 6 + 2*TG**(-2)*S*
     +    (S+U1)**(-1)*MS2*MG2 - 2*TG**(-2)*S*(S+U1)**(-1)*MG2**2 - 2*
     +    TG**(-2)*S*MG2 - 4*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 + 8*
     +    TG**(-2)*(S+U1)**(-1)*MS2**2*MG2 - 4*TG**(-2)*(S+U1)**(-1)*
     +    MS2**3 + 2*TG**(-1)*S*(S+U1)**(-1)*MS2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MG2 - 2*TG**(-1)*S - 8*TG**(-1)*(S+U1)**(-1)*MS2
     +    *MG2 + 8*TG**(-1)*(S+U1)**(-1)*MS2**2 + 4*S**(-1)*T1 + 2*
     +    S**(-1)*U1 + 4*S**(-1)*MS2 + 4*S*U1**(-1) - 2*S*(S+U1)**(-1)
     +     + 8*T1*U1**(-1) + 8*U1**(-1)*MS2 - 4*(S+U1)**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(47)*N*CF**2 * (  - 8*S**(-1)*T1*U1**(-1)*
     +    MS2 - 4*S**(-1)*T1 - 8*S**(-1)*T1**2*U1**(-1) + 16*U1**(-1)*
     +    MS2 )
     +
      MGQLRH = MGQLRH + ANG4(47)*N**2*CF * ( 2*S**(-1)*T1*U1**(-1)*MS2
     +     + 2*S**(-1)*T1 + 2*S**(-1)*T1**2*U1**(-1) + 2*S**(-1)*MS2 - 
     +    4*U1**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(68)*N*CF**2*S4G*S4G2**(-1) * ( 8*TG**(-1)*
     +    S**(-1)*U1**2 + 16*TG**(-1)*S**(-1)*MS2**2 - 16*TG**(-1)*
     +    S**(-1)*MG2**2 - 8*TG**(-1)*S*U1**(-1)*MS2 + 8*TG**(-1)*S - 
     +    16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*U1**(-1)*MS2**2 + 
     +    16*TG**(-1)*U1 - 8*TG**(-1)*MS2 - 16*S**(-1)*T1 - 16*S**(-1)*
     +    MS2 - 16*S**(-1)*MG2 - 16*U1**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(68)*N*CF**2 * (  - 8*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 8*TG**(-2)*(S+U1)**(-1)*MG2**2 + 4*TG**(-2)*MG2 - 
     +    8*TG**(-1)*(S+U1)**(-1)*MS2 + 16*TG**(-1)*(S+U1)**(-1)*MG2 + 
     +    4*TG**(-1) + 8*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(68)*N**2*CF*S4G*S4G2**(-1) * ( 6 - 4*
     +    TG**(-2)*S*U1 + 2*TG**(-2)*S*MG2 - 2*TG**(-2)*S**2 + 2*
     +    TG**(-2)*U1*MG2 - 2*TG**(-2)*U1**2 - 4*TG**(-2)*MS2*MG2 + 4*
     +    TG**(-2)*MG2**2 - 4*TG**(-1)*S**(-1)*U1**2 - 8*TG**(-1)*
     +    S**(-1)*MS2**2 + 8*TG**(-1)*S**(-1)*MG2**2 + 2*TG**(-1)*S*
     +    U1**(-1)*MS2 - 2*TG**(-1)*S + 4*TG**(-1)*U1**(-1)*MS2*MG2 - 4
     +    *TG**(-1)*U1**(-1)*MS2**2 - 6*TG**(-1)*U1 + 10*TG**(-1)*MG2
     +     + 8*S**(-1)*T1 + 8*S**(-1)*MS2 + 8*S**(-1)*MG2 + 4*U1**(-1)*
     +    MS2 )
     +
      MGQLRH = MGQLRH + ANG4(68)*N**2*CF * ( 2*TG**(-2)*S + 2*TG**(-2)*
     +    U1 + 4*TG**(-2)*(S+U1)**(-1)*MS2*MG2 - 4*TG**(-2)*
     +    (S+U1)**(-1)*MG2**2 - 2*TG**(-2)*MG2 + 4*TG**(-1)*
     +    (S+U1)**(-1)*MS2 - 8*TG**(-1)*(S+U1)**(-1)*MG2 - 2*TG**(-1)
     +     - 4*(S+U1)**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(69)*N*CF**2*S4G*S4G2**(-1) * (  - 8*
     +    TG**(-1)*S**(-1)*U1 - 16*TG**(-1)*S**(-1)*MG2 - 8*TG**(-1)*S*
     +    U1**(-1) - 16*TG**(-1)*U1**(-1)*MG2 - 16*TG**(-1) - 16*
     +    S**(-1) - 16*U1**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(69)*N**2*CF*S4G*S4G2**(-1) * ( 4*TG**(-1)*
     +    S**(-1)*U1 + 8*TG**(-1)*S**(-1)*MG2 + 2*TG**(-1)*S*U1**(-1)
     +     + 4*TG**(-1)*U1**(-1)*MG2 + 6*TG**(-1) + 8*S**(-1) + 4*
     +    U1**(-1) )
     +
      MGQLRH = MGQLRH + ANG4(84)*N*CF**2 * ( 16*S**(-1)*T1*U1**(-1)*
     +    MS2**2 - 16*S**(-1)*T1*U1**(-1)*MG2**2 - 8*S**(-1)*T1*U1 - 32
     +    *S**(-1)*T1*MG2 - 16*S**(-1)*T1**2*U1**(-1)*MG2 - 8*S**(-1)*
     +    T1**2 - 16*S**(-1)*U1**(-1)*MS2*MG2**2 + 32*S**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 16*S**(-1)*U1**(-1)*MS2**3 - 8*S**(-1)*U1*MS2 - 
     +    8*S**(-1)*U1*MG2 - 16*S**(-1)*MS2*MG2 + 24*S**(-1)*MS2**2 - 8
     +    *S**(-1)*MG2**2 - 8*S*T1*U1**(-1)*(S+U1)**(-1)*MS2 - 8*S*
     +    T1**2*U1**(-1)*(S+U1)**(-1) + 16*T1*U1**(-1)*(S+U1)**(-1)*MS2
     +    *MG2 - 16*T1*U1**(-1)*(S+U1)**(-1)*MS2**2 + 8*T1*U1**(-1)*MS2
     +     + 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 + 8*T1**2*U1**(-1) + 8*
     +    U1**(-1)*MS2*MG2 - 8*U1**(-1)*MG2**2 - 8*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(84)*N**2*CF * (  - 8*S**(-1)*T1*U1**(-1)*
     +    MS2**2 + 8*S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*U1 + 16*
     +    S**(-1)*T1*MG2 + 8*S**(-1)*T1**2*U1**(-1)*MG2 + 4*S**(-1)*
     +    T1**2 + 8*S**(-1)*U1**(-1)*MS2*MG2**2 - 16*S**(-1)*U1**(-1)*
     +    MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MS2**3 + 4*S**(-1)*U1*MS2 + 4
     +    *S**(-1)*U1*MG2 + 8*S**(-1)*MS2*MG2 - 12*S**(-1)*MS2**2 + 4*
     +    S**(-1)*MG2**2 + 2*S*T1*U1**(-1)*(S+U1)**(-1)*MS2 + 2*S*T1**2
     +    *U1**(-1)*(S+U1)**(-1) - 2*S*U1**(-2)*MS2*MG2 + 2*S*U1**(-2)*
     +    MG2**2 + 2*S*U1**(-1)*MG2 - 4*T1*U1**(-2)*MS2*MG2 + 4*T1*
     +    U1**(-2)*MG2**2 - 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 + 4*T1*
     +    U1**(-1)*(S+U1)**(-1)*MS2**2 - 4*T1*U1**(-1)*MS2 + 8*T1*
     +    U1**(-1)*MG2 + 2*T1 - 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 2*
     +    T1**2*U1**(-1) + 4*U1**(-2)*MS2*MG2**2 - 8*U1**(-2)*MS2**2*
     +    MG2 + 4*U1**(-2)*MS2**3 - 6*U1**(-1)*MS2**2 + 6*U1**(-1)*
     +    MG2**2 + 2*MS2 + 6*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(85)*N*CF**2 * ( 4*S*U1**(-1)*MS2*MG2 - 4*S
     +    *U1**(-1)*MG2**2 + 4*S*MG2 + 8*T1*U1**(-1)*MS2*MG2 - 8*T1*
     +    U1**(-1)*MG2**2 + 4*T1*MG2 - 8*U1**(-1)*MS2*MG2**2 + 16*
     +    U1**(-1)*MS2**2*MG2 - 8*U1**(-1)*MS2**3 - 4*U1*MS2 + 4*U1*MG2
     +     + 12*MS2*MG2 - 8*MS2**2 - 4*MG2**2 )
     +
      MGQLRH = MGQLRH + ANG4(86)*N*CF**2 * (  - 16*S**(-1)*T1*U1**(-1)*
     +    MS2**2 + 16*S**(-1)*T1*U1**(-1)*MG2**2 - 8*S**(-1)*T1*U1 + 16
     +    *S**(-1)*T1**2*U1**(-1)*MG2 - 8*S**(-1)*T1**2 + 16*S**(-1)*
     +    U1**(-1)*MS2*MG2**2 - 32*S**(-1)*U1**(-1)*MS2**2*MG2 + 16*
     +    S**(-1)*U1**(-1)*MS2**3 + 8*S**(-1)*U1*MS2 - 8*S**(-1)*U1*MG2
     +     - 16*S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 + 8*S**(-1)*MG2**2
     +     + 8*S*T1*U1**(-1)*(S+U1)**(-1)*MS2 + 8*S*T1**2*U1**(-1)*
     +    (S+U1)**(-1) + 4*S*U1**(-1)*MS2 - 8*S*U1**(-1)*MG2 + 4*S - 16
     +    *T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 + 16*T1*U1**(-1)*
     +    (S+U1)**(-1)*MS2**2 - 16*T1*U1**(-1)*MG2 + 16*T1*(S+U1)**(-1)
     +    *MS2 + 4*T1 - 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 8*T1**2*
     +    U1**(-1) + 16*T1**2*(S+U1)**(-1) - 24*U1**(-1)*MS2*MG2 + 16*
     +    U1**(-1)*MS2**2 + 8*U1**(-1)*MG2**2 + 4*U1 + 20*MS2 - 16*MG2
     +     )
     +
      MGQLRH = MGQLRH + ANG4(86)*N**2*CF * ( 8*S**(-1)*T1*U1**(-1)*
     +    MS2**2 - 8*S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*U1 - 8*
     +    S**(-1)*T1**2*U1**(-1)*MG2 + 4*S**(-1)*T1**2 - 8*S**(-1)*
     +    U1**(-1)*MS2*MG2**2 + 16*S**(-1)*U1**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*U1**(-1)*MS2**3 - 4*S**(-1)*U1*MS2 + 4*S**(-1)*U1*MG2
     +     + 8*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 4*S**(-1)*MG2**2 - 
     +    2*S*T1*U1**(-1)*(S+U1)**(-1)*MS2 - 2*S*T1**2*U1**(-1)*
     +    (S+U1)**(-1) + 2*S*U1**(-2)*MS2*MG2 - 2*S*U1**(-2)*MG2**2 - 2
     +    *S*U1**(-1)*MS2 + 2*S*U1**(-1)*MG2 + 4*T1*U1**(-2)*MS2*MG2 - 
     +    4*T1*U1**(-2)*MG2**2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4
     +    *T1*U1**(-1)*(S+U1)**(-1)*MS2**2 - 4*T1*(S+U1)**(-1)*MS2 + 2*
     +    T1 + 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 + 2*T1**2*U1**(-1) - 4
     +    *T1**2*(S+U1)**(-1) - 4*U1**(-2)*MS2*MG2**2 + 8*U1**(-2)*
     +    MS2**2*MG2 - 4*U1**(-2)*MS2**3 + 8*U1**(-1)*MS2*MG2 - 2*
     +    U1**(-1)*MS2**2 - 6*U1**(-1)*MG2**2 - 6*MS2 + 6*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(91)*N*CF**2 * ( 8*S*U1**(-2)*MS2*MG2 + 4*S
     +    *U1**(-1)*MG2 + 8*T1*U1**(-2)*MS2*MG2 + 4*T1*U1**(-1)*MG2 + 8
     +    *U1**(-2)*MS2**2*MG2 - 8*U1**(-2)*MS2**3 + 8*U1**(-1)*MS2*MG2
     +     - 8*U1**(-1)*MS2**2 - 4*MS2 + 4*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(92)*N*CF**2 * (  - 4 - 8*S**(-1)*T1*
     +    U1**(-1)*MS2 - 8*S**(-1)*T1 - 8*S**(-1)*T1**2*U1**(-1) + 8*S*
     +    U1**(-2)*MS2 - 4*S*U1**(-1) + 8*T1*U1**(-2)*MS2 - 12*T1*
     +    U1**(-1) + 8*U1**(-2)*MS2**2 )
     +
      MGQLRH = MGQLRH + ANG4(92)*N**2*CF * ( 2 + 2*S**(-1)*T1*U1**(-1)*
     +    MS2 + 2*S**(-1)*T1 + 2*S**(-1)*T1**2*U1**(-1) + 2*S*U1**(-1)
     +     + 4*T1*U1**(-1) + 2*U1**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(97)*N**2*CF * ( 8*S*U1**(-1)*MS2*MG2 + 4*S
     +    *MG2 + 8*T1*U1**(-1)*MS2*MG2 + 4*T1*MG2 + 8*U1**(-1)*MS2**2*
     +    MG2 - 8*U1**(-1)*MS2**3 - 4*U1*MS2 + 4*U1*MG2 + 8*MS2*MG2 - 8
     +    *MS2**2 )
     +
      MGQLRH = MGQLRH + ANG4(99)*N*CF**2 * ( 16*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 16*S**(-1)*T1*U1**(-1)*MS2**2 + 32*S**(-1)*T1*U1**(-1)
     +    *MG2**2 - 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MS2**2 + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 16*
     +    S**(-1)*T1**2*(S+U1)**(-1)*MG2 + 32*S**(-1)*U1**(-1)*MS2*
     +    MG2**2 - 32*S**(-1)*U1**(-1)*MS2**2*MG2 - 16*S**(-1)*MS2*MG2
     +     + 16*S**(-1)*MG2**2 + 8*S*T1*U1**(-1)*(S+U1)**(-1)*MS2 + 8*S
     +    *T1**2*U1**(-1)*(S+U1)**(-1) - 16*T1*U1**(-1)*(S+U1)**(-1)*
     +    MS2*MG2 + 16*T1*U1**(-1)*(S+U1)**(-1)*MS2**2 - 8*T1*U1**(-1)*
     +    MS2 + 8*T1*(S+U1)**(-1)*MS2 - 16*T1**2*U1**(-1)*(S+U1)**(-1)*
     +    MG2 - 8*T1**2*U1**(-1) + 8*T1**2*(S+U1)**(-1) + 16*U1**(-1)*
     +    MS2*MG2 + 16*U1**(-1)*MG2**2 )
     +
      MGQLRH = MGQLRH + ANG4(99)*N**2*CF * (  - 4*S**(-2)*T1*MS2*MG2 - 
     +    4*S**(-2)*T1*MG2**2 + 2*S**(-2)*U1*MS2*MG2 - 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2*MG2**2 + 4*S**(-2)*MS2**3 - 12*S**(-1)
     +    *T1*U1**(-1)*MS2*MG2 + 12*S**(-1)*T1*U1**(-1)*MS2**2 - 16*
     +    S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*(S+U1)**(-1)*MS2*
     +    MG2 - 4*S**(-1)*T1*(S+U1)**(-1)*MS2**2 - 4*S**(-1)*T1*MS2 - 8
     +    *S**(-1)*T1*MG2 - 12*S**(-1)*T1**2*U1**(-1)*MG2 + 4*S**(-1)*
     +    T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2*MG2**2 + 16*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 - 2*S**(-1)*U1*MG2 + 2*S**(-1)*
     +    MS2*MG2 - 10*S**(-1)*MG2**2 - 2*S*T1*U1**(-1)*(S+U1)**(-1)*
     +    MS2 - 2*S*T1**2*U1**(-1)*(S+U1)**(-1) + 2*S*U1**(-2)*MS2*MG2
     +     - 2*S*U1**(-2)*MG2**2 + 4*S*U1**(-1)*MS2 - 2*S*U1**(-1)*MG2
     +     + 4*T1*U1**(-2)*MS2*MG2 - 4*T1*U1**(-2)*MG2**2 + 4*T1*
     +    U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4*T1*U1**(-1)*(S+U1)**(-1)*
     +    MS2**2 + 6*T1*U1**(-1)*MS2 - 8*T1*U1**(-1)*MG2 - 2*T1*
     +    (S+U1)**(-1)*MS2 )
     +
      MGQLRH = MGQLRH + ANG4(99)*N**2*CF * ( 4*T1**2*U1**(-1)*
     +    (S+U1)**(-1)*MG2 + 2*T1**2*U1**(-1) - 2*T1**2*(S+U1)**(-1) - 
     +    4*U1**(-2)*MS2*MG2**2 + 8*U1**(-2)*MS2**2*MG2 - 4*U1**(-2)*
     +    MS2**3 - 6*U1**(-1)*MS2*MG2 + 8*U1**(-1)*MS2**2 - 10*U1**(-1)
     +    *MG2**2 - 4*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(101)*N*CF**2 * ( 8*T1*MS2 + 8*MS2**2 )
     +
      MGQLRH = MGQLRH + ANG4(103)*N*CF**2 * ( 32*S**(-1)*T1*U1**(-1)*
     +    MS2*MG2 - 16*S**(-1)*T1*U1**(-1)*MS2**2 + 16*S**(-1)*T1*
     +    U1**(-1)*MG2**2 - 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 16*
     +    S**(-1)*T1*(S+U1)**(-1)*MS2**2 + 8*S**(-1)*T1*MS2 - 16*
     +    S**(-1)*T1*MG2 + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 16*S**(-1)*
     +    T1**2*(S+U1)**(-1)*MG2 + 16*S**(-1)*U1**(-1)*MS2*MG2**2 - 16*
     +    S**(-1)*U1**(-1)*MS2**3 + 4*S**(-1)*U1*MS2 - 8*S**(-1)*U1*MG2
     +     - 24*S**(-1)*MS2*MG2 + 16*S**(-1)*MS2**2 + 8*S**(-1)*MG2**2
     +     - 16*T1*U1**(-1)*MS2 + 8*T1*(S+U1)**(-1)*MS2 + 4*T1 - 8*
     +    T1**2*U1**(-1) + 8*T1**2*(S+U1)**(-1) - 8*U1**(-1)*MS2**2 + 8
     +    *U1**(-1)*MG2**2 + 4*MS2 - 8*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(103)*N**2*CF * (  - 4*S**(-2)*T1*MS2*MG2
     +     - 4*S**(-2)*T1*MG2**2 + 2*S**(-2)*U1*MS2*MG2 - 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2*MG2**2 + 4*S**(-2)*MS2**3 - 16*S**(-1)
     +    *T1*U1**(-1)*MS2*MG2 + 8*S**(-1)*T1*U1**(-1)*MS2**2 - 8*
     +    S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*(S+U1)**(-1)*MS2*
     +    MG2 - 4*S**(-1)*T1*(S+U1)**(-1)*MS2**2 - 6*S**(-1)*T1*MS2 - 8
     +    *S**(-1)*T1**2*U1**(-1)*MG2 + 4*S**(-1)*T1**2*(S+U1)**(-1)*
     +    MG2 - 8*S**(-1)*U1**(-1)*MS2*MG2**2 + 8*S**(-1)*U1**(-1)*
     +    MS2**3 + 4*S**(-1)*MS2*MG2 - 6*S**(-1)*MS2**2 - 6*S**(-1)*
     +    MG2**2 + 8*T1*U1**(-1)*MS2 - 2*T1*(S+U1)**(-1)*MS2 + 4*T1**2*
     +    U1**(-1) - 2*T1**2*(S+U1)**(-1) + 4*U1**(-1)*MS2**2 - 4*
     +    U1**(-1)*MG2**2 )
     +
      MGQLRH = MGQLRH + ANG4(104)*N*CF**2 * (  - 32*S**(-1)*T1*U1**(-1)
     +    *MS2*MG2 + 16*S**(-1)*T1*U1**(-1)*MS2**2 - 16*S**(-1)*T1*
     +    U1**(-1)*MG2**2 + 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 16*
     +    S**(-1)*T1*(S+U1)**(-1)*MS2**2 - 16*S**(-1)*T1**2*U1**(-1)*
     +    MG2 + 16*S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)
     +    *MS2*MG2**2 + 16*S**(-1)*U1**(-1)*MS2**3 + 8*S**(-1)*MS2*MG2
     +     - 8*S**(-1)*MG2**2 - 16*S*T1*U1**(-1) - 16*S*U1**(-1)*MS2 - 
     +    16*S*U1**(-1)*MG2 - 8*S - 8*S**2*U1**(-1) - 16*T1*U1**(-1)*
     +    MS2 - 32*T1*U1**(-1)*MG2 + 8*T1*(S+U1)**(-1)*MS2 - 8*T1**2*
     +    U1**(-1) + 8*T1**2*(S+U1)**(-1) - 32*U1**(-1)*MS2*MG2 + 8*
     +    U1**(-1)*MS2**2 - 8*U1**(-1)*MG2**2 + 8*MS2 - 16*MG2 )
     +
      MGQLRH = MGQLRH + ANG4(104)*N**2*CF * ( 4*S**(-2)*T1*MS2*MG2 + 4*
     +    S**(-2)*T1*MG2**2 - 2*S**(-2)*U1*MS2*MG2 + 2*S**(-2)*U1*
     +    MG2**2 + 4*S**(-2)*MS2*MG2**2 - 4*S**(-2)*MS2**3 + 16*S**(-1)
     +    *T1*U1**(-1)*MS2*MG2 - 8*S**(-1)*T1*U1**(-1)*MS2**2 + 8*
     +    S**(-1)*T1*U1**(-1)*MG2**2 - 4*S**(-1)*T1*(S+U1)**(-1)*MS2*
     +    MG2 + 4*S**(-1)*T1*(S+U1)**(-1)*MS2**2 + 2*S**(-1)*T1*MS2 + 8
     +    *S**(-1)*T1*MG2 + 8*S**(-1)*T1**2*U1**(-1)*MG2 - 4*S**(-1)*
     +    T1**2*(S+U1)**(-1)*MG2 + 8*S**(-1)*U1**(-1)*MS2*MG2**2 - 8*
     +    S**(-1)*U1**(-1)*MS2**3 - 2*S**(-1)*U1*MS2 + 4*S**(-1)*U1*MG2
     +     + 4*S**(-1)*MS2*MG2 - 2*S**(-1)*MS2**2 + 6*S**(-1)*MG2**2 + 
     +    8*S*T1*U1**(-1) + 8*S*U1**(-1)*MS2 + 8*S*U1**(-1)*MG2 + 6*S
     +     + 4*S**2*U1**(-1) + 8*T1*U1**(-1)*MS2 + 16*T1*U1**(-1)*MG2
     +     - 2*T1*(S+U1)**(-1)*MS2 + 4*T1 + 4*T1**2*U1**(-1) - 2*T1**2*
     +    (S+U1)**(-1) + 16*U1**(-1)*MS2*MG2 - 4*U1**(-1)*MS2**2 + 4*
     +    U1**(-1)*MG2**2 + 2*U1 + 12*MG2 )
     +
      MGQLRH = MGQLRH + COLO1(9)*N*CF**2*(S4+MS2) * (  - 16*TG**(-2)*S*
     +    T1*S4**(-1)*(S+U1)**(-2)*MS2 - 16*TG**(-2)*S*T1**2*S4**(-1)*
     +    (S+U1)**(-3)*MS2 - 8*TG**(-2)*S*S4**(-1)*(S+U1)**(-1)*MS2 + 8
     +    *TG**(-2)*T1*U1*S4**(-1)*(S+U1)**(-1) + 16*TG**(-2)*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-2) + 16*TG**(-2)*T1**3*U1*S4**(-1)*
     +    (S+U1)**(-3) - 16*S*T1*U1**2*S4**(-1)*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 + 16*S*T1*U1**2*S4**(-1)*
     +    (S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2) - 16*S*T1**2*U1**2*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 32*S*
     +    T1**2*U1**2*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) + 
     +    32*S*T1**3*U1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2) - 8*S*U1**2*S4**(-1)*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MS2 - 32*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 + 8*S**2*T1*U1*
     +    S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2) - 32*S**2*T1**2
     +    *U1*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MS2 )
     +
      MGQLRH = MGQLRH + COLO1(9)*N*CF**2*(S4+MS2) * ( 16*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) + 16*S**2*T1**3
     +    *U1*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2) - 16*S**2*
     +    U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 16*
     +    S**3*T1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MS2 - 
     +    16*S**3*T1**2*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MS2 - 8*S**3*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MS2 + 8*T1*U1**3*S4**(-1)*(S+U1)**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2) + 16*T1**2*U1**3*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2) + 16*T1**3*U1**3*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2) )

      MGPLLH = 0D0
     +  + N*CF**2*(S4+MS2)**(-1) * ( 16*M2*S**(-1)*U1**(-1)*S4*MS2*
     +    MG2**(-1) + 16*M2*S**(-1)*U1**(-1)*MS2 - 16*M2*S**(-1)*
     +    U1**(-1)*MS2**2*MG2**(-1) - 4*M2*S**(-1)*MS2*MG2**(-1) - 8*M2
     +    *U1**(-2)*MS2**2*MG2**(-1) - 8*M2**2*S**(-1)*U1**(-1)*MS2*
     +    MG2**(-1) + 4*S**(-1)*T1*MS2*MG2**(-1) - 16*S**(-1)*U1**(-1)*
     +    S4*MS2 + 16*S**(-1)*U1**(-1)*S4*MS2**2*MG2**(-1) - 8*S**(-1)*
     +    U1**(-1)*S4**2*MS2*MG2**(-1) - 8*S**(-1)*U1**(-1)*MS2*MG2 + 
     +    16*S**(-1)*U1**(-1)*MS2**2 - 8*S**(-1)*U1**(-1)*MS2**3*
     +    MG2**(-1) + 4*S**(-1)*S4*MS2*MG2**(-1) + 4*S**(-1)*MS2 - 4*
     +    S**(-1)*MS2**2*MG2**(-1) + 8*U1**(-2)*S4*MS2**2*MG2**(-1) + 8
     +    *U1**(-2)*MS2**2 - 8*U1**(-2)*MS2**3*MG2**(-1) - 8*U1**(-1)*
     +    MS2**2*MG2**(-1) )
      MGPLLH = MGPLLH + N*CF**2*(S4+MS2) * (  - 16*TG**(-2)*S*T1*
     +    S4**(-1)*(S+U1)**(-2)*MG2 - 16*TG**(-2)*S*T1**2*S4**(-1)*
     +    (S+U1)**(-3)*MG2 - 16*S*T1*U1**2*S4**(-1)*(S+U1)**(-2)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S*T1**2*U1**2*S4**(-1)*
     +    (S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1*S4**(-1)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 - 8*S*U1**2*S4**(-1)*(S+U1)**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MG2 - 32*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 32*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**2*
     +    U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S**2
     +    *S4**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**3*T1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 16*S**3*T1**2*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 - 8*S**3*
     +    S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 )
      MGPLLH = MGPLLH + N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)*S**(-1)*
     +    MS2*MG2 - 16*TG**(-1)*S**(-1)*MG2**2 + 16*TG**(-1)*S*U1**(-1)
     +    *MG2 - 16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*TG**(-1)*U1**(-1)*
     +    MG2**2 + 16*TG**(-1)*MG2 - 4*S**(-1)*T1*MS2*MG2**(-1) - 24*
     +    S**(-1)*U1**(-1)*MS2*MG2 + 8*S**(-1)*U1**(-1)*MS2**3*
     +    MG2**(-1) + 16*S**(-1)*U1**(-1)*MG2**2 + 4*S**(-1)*MS2**2*
     +    MG2**(-1) - 20*S**(-1)*MG2 - 8*U1**(-2)*MS2*MG2 + 8*U1**(-2)*
     +    MS2**3*MG2**(-1) + 8*U1**(-1)*MS2**2*MG2**(-1) + 8*U1**(-1)*
     +    MG2 )
      MGPLLH = MGPLLH + N*CF**2 * (  - 8*M2*S**(-1)*U1**(-1)*MS2*
     +    MG2**(-1) + 8*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) + 8*S**(-1)*
     +    U1**(-1)*MS2 - 16*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 8*
     +    S**(-1)*U1**(-1)*MG2 - 4*S**(-1)*MS2*MG2**(-1) - 8*U1**(-2)*
     +    MS2**2*MG2**(-1) )
      MGPLLH = MGPLLH + N**2*CF*(S4+MS2)**(-1)*S4G*S4G2**(-1) * ( 8*
     +    TG**(-2)*S*MS2*MG2 + 8*TG**(-2)*T1*MS2*MG2 + 4*TG**(-1)*S*
     +    U1**(-1)*MS2*MG2 + 4*TG**(-1)*T1*U1**(-1)*MS2*MG2 )
      MGPLLH = MGPLLH + N**2*CF*(S4+MS2)**(-1) * (  - 4*M2*TG**(-2)*MS2
     +     + 4*M2*TG**(-2)*MS2**2*MG2**(-1) - 2*M2*TG**(-1)*S**(-2)*MS2
     +    *MG2 - 2*M2*TG**(-1)*S**(-2)*MS2*MG2**2*(S4G-S)**(-1) + 2*M2*
     +    TG**(-1)*S**(-2)*MS2**2*MG2*(S4G-S)**(-1) + 2*M2*TG**(-1)*
     +    S**(-2)*MS2**2 - 2*M2*TG**(-1)*S**(-1)*MS2*MG2*(S4G-S)**(-1)
     +     - 2*M2*TG**(-1)*S**(-1)*MS2 - 4*M2*S**(-2)*U1**(-1)*MS2*MG2
     +     - 4*M2*S**(-2)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) - 4*M2*
     +    S**(-2)*U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*
     +    U1**(-1)*MS2**2 - 8*M2*S**(-1)*U1**(-1)*S4*MS2*MG2**(-1) - 2*
     +    M2*S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 6*M2*S**(-1)*
     +    U1**(-1)*MS2 + 12*M2*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) + 2*M2
     +    *S**(-1)*MS2*MG2**(-1) - 2*M2*U1**(-1)*MS2*MG2**(-1) + 4*
     +    M2**2*S**(-1)*U1**(-1)*MS2*MG2**(-1) + 12*TG**(-2)*S*MS2 + 8*
     +    TG**(-2)*T1*MS2 + 4*TG**(-2)*T1*MS2**2*MG2**(-1) + 4*TG**(-2)
     +    *S4*MS2 - 4*TG**(-2)*S4*MS2**2*MG2**(-1) + 4*TG**(-2)*MS2*MG2
     +     - 8*TG**(-2)*MS2**2 )
      MGPLLH = MGPLLH + N**2*CF*(S4+MS2)**(-1) * ( 4*TG**(-2)*MS2**3*
     +    MG2**(-1) + 2*TG**(-1)*S**(-2)*T1*MS2*MG2 + 2*TG**(-1)*
     +    S**(-2)*T1*MS2*MG2**2*(S4G-S)**(-1) + 6*TG**(-1)*S**(-2)*T1*
     +    MS2**2*MG2*(S4G-S)**(-1) + 6*TG**(-1)*S**(-2)*T1*MS2**2 + 2*
     +    TG**(-1)*S**(-2)*S4*MS2*MG2 + 2*TG**(-1)*S**(-2)*S4*MS2*
     +    MG2**2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*S4*MS2**2*MG2*
     +    (S4G-S)**(-1) - 2*TG**(-1)*S**(-2)*S4*MS2**2 + 2*TG**(-1)*
     +    S**(-2)*MS2*MG2**2 + 2*TG**(-1)*S**(-2)*MS2*MG2**3*
     +    (S4G-S)**(-1) - 4*TG**(-1)*S**(-2)*MS2**2*MG2 - 4*TG**(-1)*
     +    S**(-2)*MS2**2*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*
     +    MS2**3*MG2*(S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*MS2**3 + 2*
     +    TG**(-1)*S**(-1)*T1*MS2 - 8*TG**(-1)*S**(-1)*T1*MS2**2*
     +    MG2**(-1) + 2*TG**(-1)*S**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 2*
     +    TG**(-1)*S**(-1)*S4*MS2 + 6*TG**(-1)*S**(-1)*MS2*MG2 + 6*
     +    TG**(-1)*S**(-1)*MS2*MG2**2*(S4G-S)**(-1) + 2*TG**(-1)*
     +    S**(-1)*MS2**2*MG2*(S4G-S)**(-1) )
      MGPLLH = MGPLLH + N**2*CF*(S4+MS2)**(-1) * ( 2*TG**(-1)*S**(-1)*
     +    MS2**2 + 4*TG**(-1)*S*MS2*MG2**(-1) + 4*TG**(-1)*T1*MS2*
     +    MG2**(-1) + 2*TG**(-1)*MS2*MG2*(S4G-S)**(-1) + 4*TG**(-1)*MS2
     +     - 8*TG**(-1)*MS2**2*MG2**(-1) + 4*S**(-2)*U1**(-1)*S4*MS2*
     +    MG2 + 4*S**(-2)*U1**(-1)*S4*MS2*MG2**2*(S4G-S)**(-1) + 4*
     +    S**(-2)*U1**(-1)*S4*MS2**2*MG2*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*S4*MS2**2 + 4*S**(-2)*U1**(-1)*MS2*MG2**2 + 4*
     +    S**(-2)*U1**(-1)*MS2*MG2**3*(S4G-S)**(-1) - 4*S**(-2)*
     +    U1**(-1)*MS2**3*MG2*(S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*MS2**3
     +     - 2*S**(-2)*MS2*MG2 - 2*S**(-2)*MS2*MG2**2*(S4G-S)**(-1) - 6
     +    *S**(-2)*MS2**2*MG2*(S4G-S)**(-1) - 6*S**(-2)*MS2**2 + 2*
     +    S**(-1)*U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 6*S**(-1)*
     +    U1**(-1)*S4*MS2 - 12*S**(-1)*U1**(-1)*S4*MS2**2*MG2**(-1) + 4
     +    *S**(-1)*U1**(-1)*S4**2*MS2*MG2**(-1) + 4*S**(-1)*U1**(-1)*
     +    MS2*MG2 + 4*S**(-1)*U1**(-1)*MS2*MG2**2*(S4G-S)**(-1) - 4*
     +    S**(-1)*U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) )
      MGPLLH = MGPLLH + N**2*CF*(S4+MS2)**(-1) * (  - 12*S**(-1)*
     +    U1**(-1)*MS2**2 + 8*S**(-1)*U1**(-1)*MS2**3*MG2**(-1) - 2*
     +    S**(-1)*S4*MS2*MG2**(-1) + 8*S**(-1)*MS2**2*MG2**(-1) + 2*
     +    U1**(-1)*S4*MS2*MG2**(-1) + 2*U1**(-1)*MS2*MG2*(S4G-S)**(-1)
     +     + 2*U1**(-1)*MS2 )
      MGPLLH = MGPLLH + N**2*CF*S4G*S4G2**(-1) * (  - 12*TG**(-2)*S*MS2
     +     - 4*TG**(-2)*S*MG2 - 8*TG**(-2)*T1*MS2 - 4*TG**(-2)*T1*
     +    MS2**2*MG2**(-1) + 4*TG**(-2)*T1*MG2 + 12*TG**(-2)*MS2*MG2 + 
     +    4*TG**(-2)*MS2**2 - 4*TG**(-2)*MS2**3*MG2**(-1) - 12*TG**(-2)
     +    *MG2**2 - 8*TG**(-1)*S**(-2)*T1*MS2*MG2 - 6*TG**(-1)*S**(-2)*
     +    T1*MS2**2 - 2*TG**(-1)*S**(-2)*T1*MG2**2 + 2*TG**(-1)*S**(-2)
     +    *MS2*MG2**2 - 6*TG**(-1)*S**(-2)*MS2**2*MG2 - 2*TG**(-1)*
     +    S**(-2)*MS2**3 + 6*TG**(-1)*S**(-2)*MG2**3 - 2*TG**(-1)*
     +    S**(-1)*T1*MS2 + 8*TG**(-1)*S**(-1)*T1*MS2**2*MG2**(-1) - 2*
     +    TG**(-1)*S**(-1)*T1*MG2 - 8*TG**(-1)*S**(-1)*MS2*MG2 - 2*
     +    TG**(-1)*S**(-1)*MS2**2 + 10*TG**(-1)*S**(-1)*MG2**2 - 8*
     +    TG**(-1)*S*U1**(-1)*MG2 - 4*TG**(-1)*S*MS2*MG2**(-1) + 4*
     +    TG**(-1)*T1*U1**(-1)*MG2 - 4*TG**(-1)*T1*MS2*MG2**(-1) + 12*
     +    TG**(-1)*U1**(-1)*MS2*MG2 - 12*TG**(-1)*U1**(-1)*MG2**2 - 4*
     +    TG**(-1)*MS2 + 8*TG**(-1)*MS2**2*MG2**(-1) - 16*TG**(-1)*MG2
     +     - 4*S**(-2)*U1**(-1)*MS2*MG2**2 )
      MGPLLH = MGPLLH + N**2*CF*S4G*S4G2**(-1) * (  - 4*S**(-2)*
     +    U1**(-1)*MS2**2*MG2 + 4*S**(-2)*U1**(-1)*MS2**3 + 4*S**(-2)*
     +    U1**(-1)*MG2**3 + 8*S**(-2)*MS2*MG2 + 6*S**(-2)*MS2**2 + 2*
     +    S**(-2)*MG2**2 + 12*S**(-1)*U1**(-1)*MS2*MG2 - 8*S**(-1)*
     +    U1**(-1)*MS2**3*MG2**(-1) - 4*S**(-1)*U1**(-1)*MG2**2 - 2*
     +    S**(-1)*MS2 - 8*S**(-1)*MS2**2*MG2**(-1) + 10*S**(-1)*MG2 - 
     +    12*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + N**2*CF * ( 2*M2*TG**(-1)*S**(-2)*MS2*MG2*
     +    (S4G-S)**(-1) - 2*M2*TG**(-1)*S**(-2)*MG2**2*(S4G-S)**(-1) - 
     +    2*M2*TG**(-1)*S**(-1)*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*
     +    U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 4*M2*S**(-2)*U1**(-1)*MG2**2
     +    *(S4G-S)**(-1) + 4*M2*S**(-1)*U1**(-1)*MS2*MG2**(-1) - 2*M2*
     +    S**(-1)*U1**(-1)*MG2*(S4G-S)**(-1) - 4*TG**(-2)*MS2 + 4*
     +    TG**(-2)*MS2**2*MG2**(-1) + 6*TG**(-1)*S**(-2)*T1*MS2*MG2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*T1*MG2**2*(S4G-S)**(-1) - 
     +    2*TG**(-1)*S**(-2)*S4*MS2*MG2*(S4G-S)**(-1) + 2*TG**(-1)*
     +    S**(-2)*S4*MG2**2*(S4G-S)**(-1) - 4*TG**(-1)*S**(-2)*MS2*
     +    MG2**2*(S4G-S)**(-1) + 10*TG**(-1)*S**(-2)*MS2**2*MG2*
     +    (S4G-S)**(-1) + 2*TG**(-1)*S**(-2)*MS2**2 - 2*TG**(-1)*
     +    S**(-2)*MG2**2 - 6*TG**(-1)*S**(-2)*MG2**3*(S4G-S)**(-1) + 2*
     +    TG**(-1)*S**(-1)*S4*MG2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-1)*
     +    MS2*MG2*(S4G-S)**(-1) - 2*TG**(-1)*S**(-1)*MS2 - 2*TG**(-1)*
     +    S**(-1)*MG2 )
      MGPLLH = MGPLLH + N**2*CF * (  - 6*TG**(-1)*S**(-1)*MG2**2*
     +    (S4G-S)**(-1) - 2*TG**(-1)*MG2*(S4G-S)**(-1) + 4*S**(-2)*
     +    U1**(-1)*S4*MS2*MG2*(S4G-S)**(-1) + 4*S**(-2)*U1**(-1)*S4*
     +    MG2**2*(S4G-S)**(-1) - 8*S**(-2)*U1**(-1)*MS2*MG2 + 4*S**(-2)
     +    *U1**(-1)*MS2**2*MG2*(S4G-S)**(-1) - 4*S**(-2)*U1**(-1)*
     +    MS2**2 - 4*S**(-2)*U1**(-1)*MG2**2 - 4*S**(-2)*U1**(-1)*
     +    MG2**3*(S4G-S)**(-1) - 6*S**(-2)*MS2*MG2*(S4G-S)**(-1) - 2*
     +    S**(-2)*MG2**2*(S4G-S)**(-1) - 4*S**(-1)*U1**(-1)*S4*MS2*
     +    MG2**(-1) + 2*S**(-1)*U1**(-1)*S4*MG2*(S4G-S)**(-1) - 8*
     +    S**(-1)*U1**(-1)*MS2*MG2*(S4G-S)**(-1) - 2*S**(-1)*U1**(-1)*
     +    MS2 + 12*S**(-1)*U1**(-1)*MS2**2*MG2**(-1) - 6*S**(-1)*
     +    U1**(-1)*MG2 - 8*S**(-1)*U1**(-1)*MG2**2*(S4G-S)**(-1) + 2*
     +    S**(-1)*MS2*MG2**(-1) - 2*U1**(-1)*MS2*MG2**(-1) - 2*U1**(-1)
     +    *MG2*(S4G-S)**(-1) )
      MGPLLH = MGPLLH + ANG4(21)*N*CF**2*S4G*S4G2**(-1) * ( 8*TG**(-1)*
     +    S**(-1)*U1*MG2 - 16*TG**(-1)*S**(-1)*MS2*MG2 + 16*TG**(-1)*
     +    S**(-1)*MG2**2 + 32*TG**(-1)*S*U1**(-1)*MG2 + 32*TG**(-1)*
     +    U1**(-1)*MG2**2 + 40*TG**(-1)*MG2 + 16*S**(-1)*MG2 + 32*
     +    U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(21)*N*CF**2 * ( 8*TG**(-2)*(S+U1)**(-1)*
     +    MS2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*MG2 - 8
     +    *TG**(-1)*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(21)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-2)*
     +    S*MG2 + 2*TG**(-2)*U1*MG2 - 4*TG**(-2)*MS2*MG2 - 4*TG**(-2)*
     +    MG2**2 - 2*TG**(-1)*S**(-1)*U1*MG2 + 4*TG**(-1)*S**(-1)*MS2*
     +    MG2 - 4*TG**(-1)*S**(-1)*MG2**2 - 16*TG**(-1)*S*U1**(-1)*MG2
     +     - 16*TG**(-1)*U1**(-1)*MG2**2 - 14*TG**(-1)*MG2 - 4*S**(-1)*
     +    MG2 - 16*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(21)*N**2*CF * (  - 4*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*MG2**2 - 2*TG**(-2)*MG2 + 
     +    4*TG**(-1)*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(22)*N*CF**2 * (  - 8*TG**(-2)*MS2*MG2 )
      MGPLLH = MGPLLH + ANG4(23)*N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)
     +    *S**(-1)*MG2 + 16*TG**(-1)*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(23)*N**2*CF*S4G*S4G2**(-1) * (  - 4*
     +    TG**(-1)*S**(-1)*MG2 - 8*TG**(-1)*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(43)*N*CF**2 * (  - 8*TG**(-2)*S*MS2*MG2 - 
     +    8*TG**(-2)*MS2*MG2**2 + 8*TG**(-2)*MS2**2*MG2 - 8*TG**(-1)*
     +    MS2*MG2 )
      MGPLLH = MGPLLH + ANG4(44)*N*CF**2*S4G*S4G2**(-1) * ( 32*TG**(-1)
     +    *S*U1**(-1)*MG2**2 + 24*TG**(-1)*S*MG2 + 16*TG**(-1)*S**2*
     +    U1**(-1)*MG2 - 16*TG**(-1)*U1**(-1)*MS2**2*MG2 + 16*TG**(-1)*
     +    U1**(-1)*MG2**3 + 8*TG**(-1)*U1*MG2 - 24*TG**(-1)*MS2*MG2 + 
     +    24*TG**(-1)*MG2**2 + 32*S*U1**(-1)*MG2 + 16*T1*U1**(-1)*MG2
     +     + 16*U1**(-1)*MS2*MG2 + 16*U1**(-1)*MG2**2 + 24*MG2 )
      MGPLLH = MGPLLH + ANG4(44)*N*CF**2 * ( 4*TG**(-2)*S*(S+U1)**(-1)*
     +    MS2*MG2 - 4*TG**(-2)*S*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*S*MG2
     +     + 16*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 - 8*TG**(-2)*
     +    (S+U1)**(-1)*MS2**2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**3 + 8*
     +    TG**(-2)*MS2*MG2 - 8*TG**(-2)*MG2**2 - 4*TG**(-1)*S*
     +    (S+U1)**(-1)*MG2 + 16*TG**(-1)*(S+U1)**(-1)*MS2*MG2 - 16*
     +    TG**(-1)*(S+U1)**(-1)*MG2**2 - 8*TG**(-1)*MG2 + 8*S**(-1)*MG2
     +     + 16*U1**(-1)*MG2 - 8*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(44)*N**2*CF*S4G*S4G2**(-1) * ( 2*TG**(-2)*
     +    S*U1*MG2 - 6*TG**(-2)*S*MS2*MG2 - 2*TG**(-2)*S*MG2**2 + 2*
     +    TG**(-2)*S**2*MG2 + 4*TG**(-2)*MS2**2*MG2 - 4*TG**(-2)*MG2**3
     +     - 16*TG**(-1)*S*U1**(-1)*MG2**2 - 8*TG**(-1)*S*MG2 - 8*
     +    TG**(-1)*S**2*U1**(-1)*MG2 + 8*TG**(-1)*U1**(-1)*MS2**2*MG2
     +     - 8*TG**(-1)*U1**(-1)*MG2**3 - 2*TG**(-1)*U1*MG2 + 6*
     +    TG**(-1)*MS2*MG2 - 14*TG**(-1)*MG2**2 - 16*S*U1**(-1)*MG2 - 8
     +    *T1*U1**(-1)*MG2 - 8*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MG2**2 - 
     +    10*MG2 )
      MGPLLH = MGPLLH + ANG4(44)*N**2*CF * (  - 2*TG**(-2)*S*
     +    (S+U1)**(-1)*MS2*MG2 + 2*TG**(-2)*S*(S+U1)**(-1)*MG2**2 - 2*
     +    TG**(-2)*S*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MS2*MG2**2 + 4*
     +    TG**(-2)*(S+U1)**(-1)*MS2**2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*
     +    MG2**3 + 2*TG**(-1)*S*(S+U1)**(-1)*MG2 - 8*TG**(-1)*
     +    (S+U1)**(-1)*MS2*MG2 + 8*TG**(-1)*(S+U1)**(-1)*MG2**2 - 4*
     +    S**(-1)*MG2 - 8*U1**(-1)*MG2 + 4*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(47)*N*CF**2 * ( 8*S**(-1)*T1*U1**(-1)*MG2
     +     )
      MGPLLH = MGPLLH + ANG4(47)*N**2*CF * (  - 2*S**(-1)*T1*U1**(-1)*
     +    MG2 - 2*S**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(68)*N*CF**2*S4G*S4G2**(-1) * (  - 32*
     +    TG**(-1)*S**(-1)*MS2*MG2 + 32*TG**(-1)*S**(-1)*MG2**2 + 8*
     +    TG**(-1)*S*U1**(-1)*MG2 - 16*TG**(-1)*U1**(-1)*MS2*MG2 + 16*
     +    TG**(-1)*U1**(-1)*MG2**2 + 8*TG**(-1)*MG2 + 32*S**(-1)*MG2 + 
     +    16*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(68)*N*CF**2 * ( 8*TG**(-2)*(S+U1)**(-1)*
     +    MS2*MG2 - 8*TG**(-2)*(S+U1)**(-1)*MG2**2 - 4*TG**(-2)*MG2 - 8
     +    *TG**(-1)*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(68)*N**2*CF*S4G*S4G2**(-1) * (  - 2*
     +    TG**(-2)*S*MG2 - 2*TG**(-2)*U1*MG2 + 4*TG**(-2)*MS2*MG2 - 4*
     +    TG**(-2)*MG2**2 + 16*TG**(-1)*S**(-1)*MS2*MG2 - 16*TG**(-1)*
     +    S**(-1)*MG2**2 - 2*TG**(-1)*S*U1**(-1)*MG2 + 4*TG**(-1)*
     +    U1**(-1)*MS2*MG2 - 4*TG**(-1)*U1**(-1)*MG2**2 - 10*TG**(-1)*
     +    MG2 - 16*S**(-1)*MG2 - 4*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(68)*N**2*CF * (  - 4*TG**(-2)*(S+U1)**(-1)
     +    *MS2*MG2 + 4*TG**(-2)*(S+U1)**(-1)*MG2**2 + 2*TG**(-2)*MG2 + 
     +    4*TG**(-1)*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(69)*N*CF**2*S4G*S4G2**(-1) * ( 16*TG**(-1)
     +    *S**(-1)*MG2 + 16*TG**(-1)*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(69)*N**2*CF*S4G*S4G2**(-1) * (  - 8*
     +    TG**(-1)*S**(-1)*MG2 - 4*TG**(-1)*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(84)*N*CF**2 * (  - 32*S**(-1)*T1*U1**(-1)*
     +    MS2*MG2 + 32*S**(-1)*T1*U1**(-1)*MG2**2 + 32*S**(-1)*T1*MG2
     +     + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 32*S**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 16*S**(-1)*U1**(-1)*MS2**2*MG2 + 16*S**(-1)*U1**(-1)
     +    *MG2**3 + 16*S**(-1)*U1*MG2 - 32*S**(-1)*MS2*MG2 + 32*S**(-1)
     +    *MG2**2 + 8*S*T1*U1**(-1)*(S+U1)**(-1)*MG2 + 16*T1*U1**(-1)*
     +    (S+U1)**(-1)*MS2*MG2 - 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 - 8
     +    *T1*U1**(-1)*MG2 - 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 8*
     +    U1**(-1)*MS2*MG2 + 8*U1**(-1)*MG2**2 + 8*MG2 )
      MGPLLH = MGPLLH + ANG4(84)*N**2*CF * ( 16*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 16*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1*MG2 - 8*
     +    S**(-1)*T1**2*U1**(-1)*MG2 + 16*S**(-1)*U1**(-1)*MS2*MG2**2
     +     - 8*S**(-1)*U1**(-1)*MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MG2**3
     +     - 8*S**(-1)*U1*MG2 + 16*S**(-1)*MS2*MG2 - 16*S**(-1)*MG2**2
     +     - 2*S*T1*U1**(-1)*(S+U1)**(-1)*MG2 + 2*S*U1**(-2)*MS2*MG2 - 
     +    2*S*U1**(-2)*MG2**2 - 2*S*U1**(-1)*MG2 + 4*T1*U1**(-2)*MS2*
     +    MG2 - 4*T1*U1**(-2)*MG2**2 - 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*
     +    MG2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 - 4*T1*U1**(-1)*MG2
     +     + 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 + 8*U1**(-2)*MS2*MG2**2
     +     - 4*U1**(-2)*MS2**2*MG2 - 4*U1**(-2)*MG2**3 + 12*U1**(-1)*
     +    MS2*MG2 - 12*U1**(-1)*MG2**2 - 8*MG2 )
      MGPLLH = MGPLLH + ANG4(85)*N*CF**2 * (  - 4*S*U1**(-1)*MS2*MG2 + 
     +    4*S*U1**(-1)*MG2**2 - 4*S*MG2 - 8*T1*U1**(-1)*MS2*MG2 + 8*T1*
     +    U1**(-1)*MG2**2 - 4*T1*MG2 - 16*U1**(-1)*MS2*MG2**2 + 8*
     +    U1**(-1)*MS2**2*MG2 + 8*U1**(-1)*MG2**3 + 4*MS2*MG2 - 4*
     +    MG2**2 )
      MGPLLH = MGPLLH + ANG4(86)*N*CF**2 * ( 32*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 32*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1**2*
     +    U1**(-1)*MG2 + 32*S**(-1)*U1**(-1)*MS2*MG2**2 - 16*S**(-1)*
     +    U1**(-1)*MS2**2*MG2 - 16*S**(-1)*U1**(-1)*MG2**3 - 8*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 + 4*S*U1**(-1)*MG2 - 16*T1*U1**(-1)
     +    *(S+U1)**(-1)*MS2*MG2 + 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 + 
     +    16*T1*U1**(-1)*MG2 - 16*T1*(S+U1)**(-1)*MG2 + 16*T1**2*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-1)*MS2*MG2 + 8*U1**(-1)*
     +    MG2**2 - 4*MG2 )
      MGPLLH = MGPLLH + ANG4(86)*N**2*CF * (  - 16*S**(-1)*T1*U1**(-1)*
     +    MS2*MG2 + 16*S**(-1)*T1*U1**(-1)*MG2**2 + 8*S**(-1)*T1**2*
     +    U1**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2*MG2**2 + 8*S**(-1)*
     +    U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MG2**3 + 2*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 2*S*U1**(-2)*MS2*MG2 + 2*S*
     +    U1**(-2)*MG2**2 - 4*T1*U1**(-2)*MS2*MG2 + 4*T1*U1**(-2)*
     +    MG2**2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2**2 + 4*T1*(S+U1)**(-1)*MG2 - 4*T1**2*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-2)*MS2*MG2**2 + 4*
     +    U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 - 4*U1**(-1)*MS2*MG2
     +     + 4*U1**(-1)*MG2**2 )
      MGPLLH = MGPLLH + ANG4(91)*N*CF**2 * (  - 8*S*U1**(-2)*MS2*MG2 - 
     +    4*S*U1**(-1)*MG2 - 8*T1*U1**(-2)*MS2*MG2 - 4*T1*U1**(-1)*MG2
     +     - 8*U1**(-2)*MS2*MG2**2 + 8*U1**(-2)*MS2**2*MG2 + 8*U1**(-1)
     +    *MS2*MG2 - 8*U1**(-1)*MG2**2 )
      MGPLLH = MGPLLH + ANG4(92)*N*CF**2 * ( 8*S**(-1)*T1*U1**(-1)*MG2
     +     - 8*U1**(-2)*MS2*MG2 )
      MGPLLH = MGPLLH + ANG4(92)*N**2*CF * (  - 2*S**(-1)*T1*U1**(-1)*
     +    MG2 - 2*U1**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(97)*N**2*CF * (  - 8*S*U1**(-1)*MS2*MG2 - 
     +    4*S*MG2 - 8*T1*U1**(-1)*MS2*MG2 - 4*T1*MG2 - 8*U1**(-1)*MS2*
     +    MG2**2 + 8*U1**(-1)*MS2**2*MG2 + 8*MS2*MG2 - 8*MG2**2 )
      MGPLLH = MGPLLH + ANG4(99)*N*CF**2 * ( 16*S**(-1)*T1*U1**(-1)*MS2
     +    *MG2 - 48*S**(-1)*T1*U1**(-1)*MG2**2 - 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MS2*MG2 + 16*S**(-1)*T1*(S+U1)**(-1)*MG2**2 - 16
     +    *S**(-1)*T1**2*U1**(-1)*MG2 + 16*S**(-1)*T1**2*(S+U1)**(-1)*
     +    MG2 + 32*S**(-1)*U1**(-1)*MS2*MG2**2 - 32*S**(-1)*U1**(-1)*
     +    MG2**3 + 16*S**(-1)*MS2*MG2 - 16*S**(-1)*MG2**2 - 8*S*T1*
     +    U1**(-1)*(S+U1)**(-1)*MG2 - 16*T1*U1**(-1)*(S+U1)**(-1)*MS2*
     +    MG2 + 16*T1*U1**(-1)*(S+U1)**(-1)*MG2**2 + 8*T1*U1**(-1)*MG2
     +     - 8*T1*(S+U1)**(-1)*MG2 + 16*T1**2*U1**(-1)*(S+U1)**(-1)*MG2
     +     - 16*U1**(-1)*MS2*MG2 - 16*U1**(-1)*MG2**2 )
      MGPLLH = MGPLLH + ANG4(99)*N**2*CF * ( 4*S**(-2)*T1*MS2*MG2 + 4*
     +    S**(-2)*T1*MG2**2 - 2*S**(-2)*U1*MS2*MG2 + 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2**2*MG2 + 4*S**(-2)*MG2**3 - 12*S**(-1)
     +    *T1*U1**(-1)*MS2*MG2 + 28*S**(-1)*T1*U1**(-1)*MG2**2 + 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 4*S**(-1)*T1*(S+U1)**(-1)*
     +    MG2**2 + 8*S**(-1)*T1*MG2 + 12*S**(-1)*T1**2*U1**(-1)*MG2 - 4
     +    *S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2*
     +    MG2**2 + 16*S**(-1)*U1**(-1)*MG2**3 + 2*S**(-1)*U1*MG2 - 6*
     +    S**(-1)*MS2*MG2 + 14*S**(-1)*MG2**2 + 2*S*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2 - 2*S*U1**(-2)*MS2*MG2 + 2*S*U1**(-2)*MG2**2
     +     + 2*S*U1**(-1)*MG2 - 4*T1*U1**(-2)*MS2*MG2 + 4*T1*U1**(-2)*
     +    MG2**2 + 4*T1*U1**(-1)*(S+U1)**(-1)*MS2*MG2 - 4*T1*U1**(-1)*
     +    (S+U1)**(-1)*MG2**2 + 6*T1*U1**(-1)*MG2 + 2*T1*(S+U1)**(-1)*
     +    MG2 - 4*T1**2*U1**(-1)*(S+U1)**(-1)*MG2 - 8*U1**(-2)*MS2*
     +    MG2**2 + 4*U1**(-2)*MS2**2*MG2 + 4*U1**(-2)*MG2**3 - 6*
     +    U1**(-1)*MS2*MG2 )
      MGPLLH = MGPLLH + ANG4(99)*N**2*CF * ( 14*U1**(-1)*MG2**2 )
      MGPLLH = MGPLLH + ANG4(101)*N*CF**2 * (  - 8*MS2*MG2 )
      MGPLLH = MGPLLH + ANG4(103)*N*CF**2 * (  - 32*S**(-1)*T1*U1**(-1)
     +    *MG2**2 - 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MG2**2 + 8*S**(-1)*T1*MG2 - 16*S**(-1)*T1**2*
     +    U1**(-1)*MG2 + 16*S**(-1)*T1**2*(S+U1)**(-1)*MG2 + 16*S**(-1)
     +    *U1**(-1)*MS2**2*MG2 - 16*S**(-1)*U1**(-1)*MG2**3 + 4*S**(-1)
     +    *U1*MG2 - 8*S**(-1)*MS2*MG2 + 8*S**(-1)*MG2**2 - 8*T1*
     +    (S+U1)**(-1)*MG2 - 4*MG2 )
      MGPLLH = MGPLLH + ANG4(103)*N**2*CF * ( 4*S**(-2)*T1*MS2*MG2 + 4*
     +    S**(-2)*T1*MG2**2 - 2*S**(-2)*U1*MS2*MG2 + 2*S**(-2)*U1*
     +    MG2**2 - 4*S**(-2)*MS2**2*MG2 + 4*S**(-2)*MG2**3 + 16*S**(-1)
     +    *T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MG2**2 + 2*S**(-1)*T1*MG2 + 8*S**(-1)
     +    *T1**2*U1**(-1)*MG2 - 4*S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 8*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MG2**3 + 4*
     +    S**(-1)*MS2*MG2 + 4*S**(-1)*MG2**2 + 2*T1*(S+U1)**(-1)*MG2 )
      MGPLLH = MGPLLH + ANG4(104)*N*CF**2 * ( 32*S**(-1)*T1*U1**(-1)*
     +    MG2**2 + 16*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 - 16*S**(-1)*T1*
     +    (S+U1)**(-1)*MG2**2 + 16*S**(-1)*T1**2*U1**(-1)*MG2 - 16*
     +    S**(-1)*T1**2*(S+U1)**(-1)*MG2 - 16*S**(-1)*U1**(-1)*MS2**2*
     +    MG2 + 16*S**(-1)*U1**(-1)*MG2**3 - 8*S**(-1)*MS2*MG2 + 8*
     +    S**(-1)*MG2**2 + 16*S*U1**(-1)*MG2 + 32*T1*U1**(-1)*MG2 - 8*
     +    T1*(S+U1)**(-1)*MG2 + 32*U1**(-1)*MG2**2 + 8*MG2 )
      MGPLLH = MGPLLH + ANG4(104)*N**2*CF * (  - 4*S**(-2)*T1*MS2*MG2
     +     - 4*S**(-2)*T1*MG2**2 + 2*S**(-2)*U1*MS2*MG2 - 2*S**(-2)*U1*
     +    MG2**2 + 4*S**(-2)*MS2**2*MG2 - 4*S**(-2)*MG2**3 - 16*S**(-1)
     +    *T1*U1**(-1)*MG2**2 - 4*S**(-1)*T1*(S+U1)**(-1)*MS2*MG2 + 4*
     +    S**(-1)*T1*(S+U1)**(-1)*MG2**2 - 6*S**(-1)*T1*MG2 - 8*S**(-1)
     +    *T1**2*U1**(-1)*MG2 + 4*S**(-1)*T1**2*(S+U1)**(-1)*MG2 + 8*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MG2**3 - 2*
     +    S**(-1)*U1*MG2 + 4*S**(-1)*MS2*MG2 - 12*S**(-1)*MG2**2 - 8*S*
     +    U1**(-1)*MG2 - 16*T1*U1**(-1)*MG2 + 2*T1*(S+U1)**(-1)*MG2 - 
     +    16*U1**(-1)*MG2**2 - 8*MG2 )
      MGPLLH = MGPLLH + COLO1(9)*N*CF**2*(S4+MS2) * ( 16*TG**(-2)*S*T1*
     +    S4**(-1)*(S+U1)**(-2)*MG2 + 16*TG**(-2)*S*T1**2*S4**(-1)*
     +    (S+U1)**(-3)*MG2 + 8*TG**(-2)*S*S4**(-1)*(S+U1)**(-1)*MG2 + 
     +    16*S*T1*U1**2*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 + 16*S*T1**2*U1**2*S4**(-1)*(S+U1)**(-3)*
     +    (M2*(S+U1)+T1*U1)**(-2)*MG2 + 8*S*U1**2*S4**(-1)*(S+U1)**(-1)
     +    *(M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1*U1*S4**(-1)*
     +    (S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 32*S**2*T1**2*U1*
     +    S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*S**2*
     +    U1*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 16*
     +    S**3*T1*S4**(-1)*(S+U1)**(-2)*(M2*(S+U1)+T1*U1)**(-2)*MG2 + 
     +    16*S**3*T1**2*S4**(-1)*(S+U1)**(-3)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 + 8*S**3*S4**(-1)*(S+U1)**(-1)*(M2*(S+U1)+T1*U1)**(-2)*
     +    MG2 )


      MGPLRH = MGQLRH 

      M2GQH = 2.D0*0.5D0*(NS -2.D0)*(MGPLLH +MGPLRH) 
     +     +2*0.5D0*MGQLLH +1*MGQLRH

      DSSGQH = S4/(S4+MS2)/2.D0 *
     +     ALPHAS**3 * AVG * M2GQH /4.D0 /S**2 *CONV

      RETURN
      END


      REAL*8 FUNCTION DSSGQS(ALPHAS,S,T1,S4,MS,MG,EPS)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR G +Q -> SQ + SQ +QB
C***  THE 1/S4G**2 PART
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,MS,MG,MS2,MG2,M2,CF,EPS
      REAL*8 M2GQS, MGPLLS, MGPLRS, MGQLLS, MGQLRS
      REAL*8 NS,CONV,N,AVG

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2-1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N /(N**2 -1.D0)

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      TG = T1 -M2
      U1 = S4 -S -T1

      MGQLLS = 0.D0
      MGQLLS = MGQLLS + N*CF**2 * ( 4*S**(-1)*T1*MS2 - 4*S**(-1)*T1*MG2
     +     - 24*S**(-1)*U1**(-1)*MS2*MG2**2 + 24*S**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MS2**3 + 8*S**(-1)*U1**(-1)*
     +    MG2**3 + 8*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 4*S**(-1)*
     +    MG2**2 - 8*U1**(-2)*MS2*MG2**2 + 16*U1**(-2)*MS2**2*MG2 - 8*
     +    U1**(-2)*MS2**3 + 16*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MS2**2 - 8
     +    *U1**(-1)*MG2**2 )
     +
      MGQLLS = MGQLLS + N**2*CF * ( 8*TG**(-2)*T1*MS2*MG2 + 4*TG**(-2)*
     +    T1*MS2**2 + 4*TG**(-2)*T1*MG2**2 + 12*TG**(-2)*MS2*MG2**2 - 4
     +    *TG**(-2)*MS2**2*MG2 + 4*TG**(-2)*MS2**3 - 12*TG**(-2)*MG2**3
     +     - 8*TG**(-1)*S**(-1)*T1*MS2*MG2 - 8*TG**(-1)*S**(-1)*T1*
     +    MS2**2 + 16*TG**(-1)*S**(-1)*MS2*MG2**2 - 16*TG**(-1)*S**(-1)
     +    *MS2**2*MG2 + 4*TG**(-1)*S*MS2 - 4*TG**(-1)*S*MG2 + 4*
     +    TG**(-1)*T1*MS2 + 4*TG**(-1)*T1*MG2 + 8*TG**(-1)*MS2*MG2 - 8*
     +    TG**(-1)*MS2**2 - 16*TG**(-1)*MG2**2 + 8*S**(-1)*U1**(-1)*MS2
     +    *MG2**2 - 16*S**(-1)*U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)
     +    *MS2**3 + 8*S**(-1)*MS2*MG2 + 8*S**(-1)*MS2**2 - 8*MG2 )

      MGQLRS = 0.D0
      MGQLRS = MGQLRS + N*CF**2 * ( 4*S**(-1)*T1*MS2 - 4*S**(-1)*T1*MG2
     +     - 24*S**(-1)*U1**(-1)*MS2*MG2**2 + 24*S**(-1)*U1**(-1)*
     +    MS2**2*MG2 - 8*S**(-1)*U1**(-1)*MS2**3 + 8*S**(-1)*U1**(-1)*
     +    MG2**3 + 8*S**(-1)*MS2*MG2 - 4*S**(-1)*MS2**2 - 4*S**(-1)*
     +    MG2**2 - 8*U1**(-2)*MS2*MG2**2 + 16*U1**(-2)*MS2**2*MG2 - 8*
     +    U1**(-2)*MS2**3 + 16*U1**(-1)*MS2*MG2 - 8*U1**(-1)*MS2**2 - 8
     +    *U1**(-1)*MG2**2 )
     +
      MGQLRS = MGQLRS + N**2*CF * (  - 8*TG**(-2)*T1*MS2*MG2 - 4*
     +    TG**(-2)*T1*MS2**2 - 4*TG**(-2)*T1*MG2**2 + 20*TG**(-2)*MS2*
     +    MG2**2 - 12*TG**(-2)*MS2**2*MG2 - 4*TG**(-2)*MS2**3 - 4*
     +    TG**(-2)*MG2**3 + 8*TG**(-1)*S**(-1)*T1*MS2*MG2 + 8*TG**(-1)*
     +    S**(-1)*T1*MS2**2 - 16*TG**(-1)*S**(-1)*MS2**2*MG2 + 16*
     +    TG**(-1)*S**(-1)*MS2**3 + 4*TG**(-1)*S*MS2 - 4*TG**(-1)*S*MG2
     +     - 4*TG**(-1)*T1*MS2 - 4*TG**(-1)*T1*MG2 + 24*TG**(-1)*MS2*
     +    MG2 - 8*TG**(-1)*MS2**2 + 8*S**(-1)*U1**(-1)*MS2*MG2**2 - 16*
     +    S**(-1)*U1**(-1)*MS2**2*MG2 + 8*S**(-1)*U1**(-1)*MS2**3 - 8*
     +    S**(-1)*MS2*MG2 - 8*S**(-1)*MS2**2 + 8*MS2 )

      MGPLLS = MGQLLS

      MGPLRS = MGQLRS

      M2GQS = 2.D0*0.5D0*(NS -2.D0)*(MGPLLS +MGPLRS) 
     +     +2*0.5D0*MGQLLS +1*MGQLRS

      DSSGQS = ALPHAS**3 *AVG *M2GQS *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END




      REAL*8 FUNCTION DSSGQT(ALPHAS,S,T1,S4,S3G,MS,MG,EPS)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR G +Q -> SQ + SQ +QB
C***  THE 1/S3G**2 PART
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,S3G,T1,TG,U1,MS,MG,MS2,MG2,M2,CF,EPS
      REAL*8 M2GQT, MGQLLT, MGQLRT, MGPLLT, MGPLRT
      REAL*8 NS,CONV,N,AVG
      REAL*8 ANGDEF(1:11), ANA(1:3,1:9), ANB(1:3,1:9), ANC(1:3,1:9)
      REAL*8 ANGS3(1:10),XX(5:9),YY2(5:9),XPHI

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2-1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N /(N**2 -1.D0)

      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      TG = T1 -M2
      U1 = S4 -S -T1


      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +U1)/ANGDEF(1)
      ANGDEF(3) = (S +T1)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(S4 -S +2.D0*MS2)/ANGDEF(1)
      ANGDEF(7) = SQRT((S -S4)**2 -4.D0*MS2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (T1*S4 -S*(U1+2.D0*MS2))/(S+T1)/SQRT((S -S4)**2-4.D0*MS2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (U1*S4 -S*(T1+2.D0*MS2))/(S+U1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)


      ANA(2,5) = -2.D0*ANGDEF(5)*ANGDEF(3)
      ANB(2,5) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,5) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,6) = -2.D0*ANGDEF(5)*ANGDEF(2) -M2
      ANB(2,6) = 
     +     -2.D0*ANGDEF(4)*ANGDEF(7) +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,6) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,7) = +2.D0*ANGDEF(4)*ANGDEF(6) -M2
      ANB(2,7) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,7) =  0.D0
      ANA(2,8) = -2.D0*ANGDEF(5)*ANGDEF(3) -M2
      ANB(2,8) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,8) = -2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)
      ANA(2,9) = -2.D0*ANGDEF(5)*ANGDEF(2) -M2
      ANB(2,9) = 
     +     -2.D0*ANGDEF(4)*ANGDEF(7) +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(8)
      ANC(2,9) = +2.D0*ANGDEF(3)*ANGDEF(4)*ANGDEF(9)

      XPHI = (S3G -ANA(2,7))/ANB(2,7)

      XX(5) = (ANA(2,5) + ANB(2,5)*XPHI)
      YY2(5)= ANC(2,5)**2 * (1.D0-XPHI**2)
      XX(6) = (ANA(2,6) + ANB(2,6)*XPHI)
      YY2(6)= ANC(2,6)**2 * (1.D0-XPHI**2)
      XX(8) = (ANA(2,8) + ANB(2,8)*XPHI)
      YY2(8)= ANC(2,8)**2 * (1.D0-XPHI**2)
      XX(9) = (ANA(2,9) + ANB(2,9)*XPHI)
      YY2(9)= ANC(2,9)**2 * (1.D0-XPHI**2)

      ANGS3(1) = -XX(5)/(XX(5)**2 - YY2(5))**(1.5D0)
      ANGS3(2) = -1.D0/SQRT(XX(5)**2 - YY2(5))
      ANGS3(3) = -XX(9)/(XX(9)**2 - YY2(9))**(1.5D0)
      ANGS3(4) = -1.D0/SQRT(XX(9)**2 - YY2(9))
      ANGS3(5) = XX(9)
      ANGS3(6) = XX(6)
      ANGS3(7) = -XX(8)/(XX(8)**2 - YY2(8))**(1.5D0)
      ANGS3(8) = -1.D0/SQRT(XX(8)**2 - YY2(8))
      ANGS3(9) = -XX(6)/(XX(6)**2 - YY2(6))**(1.5D0)
      ANGS3(10) = -1.D0/SQRT(XX(6)**2 - YY2(6))

      MGQLLT = 0.D0
      MGQLLT = MGQLLT + N*CF**2 * (  - 4*S**(-1)*T1*MG2 + 4*S**(-1)*MS2
     +    *MG2 - 4*S**(-1)*MG2**2 + 4*MG2 )
     +
      MGQLLT = MGQLLT + N**2*CF * (  - 4*MG2 )
     +
      MGQLLT = MGQLLT + ANGS3(5)*N*CF**2 * ( 4*S**(-1)*MG2 )
     +
      MGQLLT = MGQLLT + ANGS3(7)*N**2*CF * (  - 4*S*MS2*MG2 - 4*S*
     +    MG2**2 - 8*T1*MG2**2 + 4*U1*MS2*MG2 - 4*U1*MG2**2 + 8*MS2*
     +    MG2**2 - 8*MG2**3 )
     +
      MGQLLT = MGQLLT + ANGS3(8)*N**2*CF * ( 8*S**(-1)*T1*MS2*MG2 + 8*
     +    S**(-1)*MS2*MG2**2 - 8*S**(-1)*MS2**2*MG2 - 4*S*MG2 - 4*T1*
     +    MG2 + 8*MS2*MG2 - 8*MG2**2 )
     +
      MGQLLT = MGQLLT + ANGS3(9)*N*CF**2 * (  - 8*T1*MS2*MG2 - 8*MS2*
     +    MG2**2 + 8*MS2**2*MG2 )
     +
      MGQLLT = MGQLLT + ANGS3(10)*N*CF**2 * (  - 8*S**(-1)*T1*MS2*MG2
     +     + 8*S**(-1)*T1*MG2**2 - 4*S**(-1)*U1*MS2*MG2 + 4*S**(-1)*U1*
     +    MG2**2 - 16*S**(-1)*MS2*MG2**2 + 8*S**(-1)*MS2**2*MG2 + 8*
     +    S**(-1)*MG2**3 - 4*T1*MG2 + 12*MS2*MG2 - 4*MG2**2 )
     +
      MGQLLT = MGQLLT + ANGS3(10)*N**2*CF * ( 8*S**(-1)*T1*MS2*MG2 + 8*
     +    S**(-1)*MS2*MG2**2 - 8*S**(-1)*MS2**2*MG2 )

      MGQLRT = 0.D0
      MGQLRT = MGQLRT + N*CF**2 * ( 4*S**(-1)*T1*MG2 + 4*S**(-1)*MS2*
     +    MG2 - 4*S**(-1)*MS2**2 - 4*MS2 )
     +
      MGQLRT = MGQLRT + N**2*CF * ( 4*MS2 )
     +
      MGQLRT = MGQLRT + ANGS3(5)*N*CF**2 * (  - 4*S**(-1)*MS2 )
     +
      MGQLRT = MGQLRT + ANGS3(7)*N**2*CF * ( 4*S*MS2*MG2 + 4*S*MG2**2
     +     + 8*T1*MG2**2 - 4*U1*MS2*MG2 + 4*U1*MG2**2 + 8*MS2*MG2**2 - 
     +    8*MS2**2*MG2 )
     +
      MGQLRT = MGQLRT + ANGS3(8)*N**2*CF * (  - 8*S**(-1)*T1*MS2*MG2 - 
     +    8*S**(-1)*MS2**2*MG2 + 8*S**(-1)*MS2**3 + 4*S*MS2 + 4*T1*MG2
     +     + 8*MS2*MG2 - 8*MS2**2 )
     +
      MGQLRT = MGQLRT + ANGS3(9)*N*CF**2 * ( 8*T1*MS2*MG2 + 8*MS2**2*
     +    MG2 - 8*MS2**3 )
     +
      MGQLRT = MGQLRT + ANGS3(10)*N*CF**2 * ( 8*S**(-1)*T1*MS2*MG2 - 8*
     +    S**(-1)*T1*MG2**2 + 4*S**(-1)*U1*MS2*MG2 - 4*S**(-1)*U1*
     +    MG2**2 - 8*S**(-1)*MS2*MG2**2 + 16*S**(-1)*MS2**2*MG2 - 8*
     +    S**(-1)*MS2**3 + 4*T1*MG2 + 4*MS2*MG2 - 8*MS2**2 - 4*MG2**2 )
     +
      MGQLRT = MGQLRT + ANGS3(10)*N**2*CF * (  - 8*S**(-1)*T1*MS2*MG2
     +     - 8*S**(-1)*MS2**2*MG2 + 8*S**(-1)*MS2**3 )

      MGPLLT = MGQLLT

      MGPLRT = MGQLRT

      M2GQT = 2.D0*0.5D0*(NS -2.D0)*(MGPLLT +MGPLRT) 
     +     + 2.D0*0.5D0*MGQLLT + 1*MGQLRT

      DSSGQT = ALPHAS**3 *AVG *M2GQT *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END



      REAL*8 FUNCTION DSSGQU(ALPHAS,S,T1,S4,MS,MG,EPS)
C***  DOUBLE-POLE PART OF CROSS SECTIONS FOR G +Q -> SQ + SQ +QB
C***  THE 1/S3G/S4G PART
C***  SUMMED OVER LL +RR +LR +RL
      IMPLICIT NONE
      REAL*8 ALPHAS,S,S4,T1,TG,U1,MS,MG,MS2,MG2,M2,CF,EPS
      REAL*8 ANGDEF(1:11), ANA(2:3,1:9), ANB(2:3,1:9), ANC(2:3,1:9)
      REAL*8 M2GQU, MGQLLU
      REAL*8 NS,CONV,N,AVG
      REAL*8 K4M1P1, K4P1P0, K4P1P1
      REAL*8 CANG4(1:4),KPI

      NS = 6.D0
      CONV = 389379660.D0
      N = 3.D0
      CF = (N**2-1.D0)/N/2.D0
      AVG = (1.D0/2.D0)**2 /N /(N**2 -1.D0)

      KPI = - 4.D0*ATAN(1.D0)
      MS2 = MS**2
      MG2 = MG**2
      M2 = MG2 -MS2
      TG = T1 -M2
      U1 = S4 -S -T1

      ANGDEF(1) = 2.D0*SQRT(S4 +MS2)
      ANGDEF(2) = (S +U1)/ANGDEF(1)
      ANGDEF(3) = (S +T1)/ANGDEF(1)
      ANGDEF(4) = S4/ANGDEF(1)
      ANGDEF(5) = (S4 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(6) = -(T1 +U1 +2.D0*MS2)/ANGDEF(1)
      ANGDEF(7) = SQRT((T1 +U1)**2 -4.D0*MS2*S)/ANGDEF(1)
      ANGDEF(8) = 
     +  (T1*S4 -S*(U1+2.D0*MS2))/(S+T1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(9) = SQRT(1 -ANGDEF(8)**2)
      ANGDEF(10) =
     +  (U1*S4 -S*(T1+2.D0*MS2))/(S+U1)/SQRT((T1+U1)**2-4.D0*MS2*S)
      ANGDEF(11) = SQRT(1 -ANGDEF(10)**2)


      ANA(2,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(2,1) = -2.D0*ANGDEF(4)*ANGDEF(7)
      ANC(2,1) =  0.D0
      ANA(2,2) = +2.D0*ANGDEF(5)*ANGDEF(6) +2.D0*MS2
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

      ANA(3,1) = +2.D0*ANGDEF(4)*ANGDEF(6)
      ANB(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(10)
      ANC(3,1) = -2.D0*ANGDEF(4)*ANGDEF(7)*ANGDEF(11)
      ANA(3,4) = -2.D0*ANGDEF(2)*ANGDEF(4)
      ANB(3,4) = -ANA(3,4)
      ANC(3,4) =  0.D0
      ANA(3,6) = -2.D0*ANGDEF(5)*ANGDEF(2)
      ANB(3,6) = -ANB(3,4)
      ANC(3,6) = -ANC(3,4)
      ANA(3,7) = +ANA(3,1) -M2
      ANB(3,7) = +ANB(3,1)
      ANC(3,7) = +ANC(3,1)

      CANG4(1) = K4M1P1(ANA(3,6),ANB(3,6),ANA(3,7),ANB(3,7),ANC(3,7))
      CANG4(2) = K4P1P0(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
      CANG4(3) = K4P1P1(ANA(2,7),ANB(2,7),ANA(2,8),ANB(2,8),ANC(2,8))
      CANG4(4) = K4P1P1(ANA(2,7),ANB(2,7),ANA(2,6),ANB(2,6),ANC(2,6))


      MGQLLU = 0.D0
      MGQLLU = MGQLLU + CANG4(1)*N*CF*KPI * (  - 2*TG**(-1)*S**(-1)*U1*
     +    MG2 - 2*TG**(-1)*S**(-1)*MS2*MG2 + 2*TG**(-1)*S**(-1)*MG2**2
     +     + 2*S**(-1)*T1*U1**(-1)*MG2 )
     +
      MGQLLU = MGQLLU + CANG4(1)*CF**2*KPI * (  - 4*S**(-1)*T1*U1**(-1)
     +    *MG2 + 4*S**(-1)*MG2 )
     +
      MGQLLU = MGQLLU + CANG4(2)*N*CF*KPI * (  - 2*TG**(-1)*S**(-1)*U1*
     +    MS2*MG2 + 2*TG**(-1)*S**(-1)*U1*MG2**2 + 16*TG**(-1)*S**(-1)*
     +    MS2*MG2**2 - 8*TG**(-1)*S**(-1)*MS2**2*MG2 - 8*TG**(-1)*
     +    S**(-1)*MG2**3 + 2*TG**(-1)*U1*MG2 + 2*TG**(-1)*MS2*MG2 + 6*
     +    TG**(-1)*MG2**2 + 4*S**(-1)*T1*U1**(-1)*MS2*MG2 - 4*S**(-1)*
     +    T1*U1**(-1)*MG2**2 - 4*S**(-1)*T1*MG2 - 4*S**(-1)*T1**2*
     +    U1**(-1)*MG2 + 6*S**(-1)*MS2*MG2 - 6*S**(-1)*MG2**2 + 4*T1*
     +    U1**(-1)*MG2 + 6*U1**(-1)*MS2*MG2 + 2*U1**(-1)*MG2**2 + 2*MG2
     +     )
     +
      MGQLLU = MGQLLU + CANG4(2)*CF**2*KPI * (  - 8*S**(-1)*T1*U1**(-1)
     +    *MS2*MG2 + 8*S**(-1)*T1*U1**(-1)*MG2**2 + 4*S**(-1)*T1*MG2 + 
     +    8*S**(-1)*T1**2*U1**(-1)*MG2 + 4*S**(-1)*U1*MG2 + 4*S**(-1)*
     +    MS2*MG2 - 4*S**(-1)*MG2**2 - 4*T1*U1**(-1)*MG2 - 12*U1**(-1)*
     +    MS2*MG2 - 4*U1**(-1)*MG2**2 )
     +
      MGQLLU = MGQLLU + CANG4(3)*N*CF*KPI * ( 4*S**(-1)*T1*MS2*MG2 - 4*
     +    S**(-1)*T1*MG2**2 + 2*S**(-1)*U1*MS2*MG2 - 2*S**(-1)*U1*
     +    MG2**2 + 8*S**(-1)*MS2*MG2**2 - 4*S**(-1)*MS2**2*MG2 - 4*
     +    S**(-1)*MG2**3 - 6*S*U1**(-1)*MS2*MG2 - 2*S*U1**(-1)*MG2**2
     +     - 2*S*MG2 - 4*T1*U1**(-1)*MS2*MG2 - 4*T1*U1**(-1)*MG2**2 + 4
     +    *U1**(-1)*MS2**2*MG2 - 4*U1**(-1)*MG2**3 + 2*U1*MG2 + 4*MS2*
     +    MG2 - 4*MG2**2 )
     +
      MGQLLU = MGQLLU + CANG4(4)*N*CF*KPI * (  - 4*TG**(-1)*S*MS2*MG2
     +     + 4*TG**(-1)*S*MG2**2 - 2*TG**(-1)*U1*MS2*MG2 + 2*TG**(-1)*
     +    U1*MG2**2 + 8*TG**(-1)*MS2**2*MG2 - 8*TG**(-1)*MG2**3 + 4*
     +    S**(-1)*T1*MS2*MG2 - 4*S**(-1)*T1*MG2**2 + 2*S**(-1)*U1*MS2*
     +    MG2 - 2*S**(-1)*U1*MG2**2 + 8*S**(-1)*MS2*MG2**2 - 4*S**(-1)*
     +    MS2**2*MG2 - 4*S**(-1)*MG2**3 + 2*S*MG2 - 16*T1*U1**(-1)*
     +    MG2**2 - 6*T1*MG2 - 8*T1**2*U1**(-1)*MG2 + 8*U1**(-1)*MS2**2*
     +    MG2 - 8*U1**(-1)*MG2**3 + 2*MS2*MG2 - 10*MG2**2 )
     +
      MGQLLU = MGQLLU + CANG4(4)*CF**2*KPI * (  - 8*S**(-1)*T1*MS2*MG2
     +     + 8*S**(-1)*T1*MG2**2 - 4*S**(-1)*U1*MS2*MG2 + 4*S**(-1)*U1*
     +    MG2**2 - 16*S**(-1)*MS2*MG2**2 + 8*S**(-1)*MS2**2*MG2 + 8*
     +    S**(-1)*MG2**3 + 16*T1*U1**(-1)*MG2**2 + 4*T1*MG2 + 8*T1**2*
     +    U1**(-1)*MG2 - 8*U1**(-1)*MS2**2*MG2 + 8*U1**(-1)*MG2**3 - 4*
     +    MS2*MG2 + 4*MG2**2 )



      M2GQU = 2.D0* 0.5D0*MGQLLU

      DSSGQU = ALPHAS**3 *AVG *M2GQU *S4/(S4 +MS2)/8.D0/S**2 *CONV
      RETURN
      END



