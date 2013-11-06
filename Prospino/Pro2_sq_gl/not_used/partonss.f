      PROGRAM PARTONSS
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      character*5 flav(3)

      data flav/'qq   ','qqp  ','gq   '/

      COMMON/IOUT/IPRINT
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/RESULT/S1,S2,S3,S4
      COMMON/IFLAVOR/IFLAVOR

      EXTERNAL BORN, SOVI, HARD, SOVISCA, HARDSCA

      CALL RSTART(12,34,56,78)

C --- PARAMETER TO CHANGE  -------------------------------------

C***  1:Q QBAR     2:QP QBAR    3:G Q 

      DO 100 IFLAVOR = 3,3
         
C***  THE MASSES:

      MS = 280.D0      
      MG = 200.D0
      MT = 175.D0

C***  START VALUE IN LOG(ETA), NUMBERS OF STEPS, STEP WIDTH IN LOG(ETA)
      
      ETALOG = -2D0
      INUM = 5
      STEP = 1D0


C***  THE NUMBER OF VEGAS CALLS ( DEFAULT = 500 )

      IBORN = 1000
      ISOVI = 500
      IHARD = 500

C***  PRINT VEGAS STATISTICS (10) OR NOT (0)

      IPRINT = 0


C --- CONSTANTS DON'T CHANGE PLEASE ------------------------------

      PI = 4.D0*DATAN(1.D0)
      CONV = 389379660.D0
      ALPHAS = 0.1D0
      SCA = 4.D0


C -------------------------------------------------------------

c$$$      IF (IFLAVOR.EQ.1) OPEN(10,FILE="partonss_qq.dat")
c$$$      IF (IFLAVOR.EQ.2) OPEN(10,FILE="partonss_qp.dat")
c$$$      IF (IFLAVOR.EQ.3) OPEN(10,FILE="partonss_gq.dat")

      PRINT 20, MS,MG,MT
      PRINT 30
c$$$      WRITE(10, 20) MS,MG,MT
c$$$      WRITE(10, 30)

 20   FORMAT(' MS = ',G10.4,' MG= ',G10.4,' MT = ',G10.4)
 30   FORMAT('  ETA        BORN          SOVI       REL     ',
     +   '   HARD       REL       SCALE')

C --- INTEGRATION BY VEGAS ------------------------------------


      DO 100  I = 1,INUM

        ETA = 10.D0**(ETALOG)
        ETALOG = ETALOG +STEP

        S = 4.D0*MS**2*(ETA + 1.D0)
        DEL = 1.D-4*MS**2 * (1 -4.D0*MS**2/S)

      IF (IFLAVOR.LE.2) THEN
         CALL INTEG(BORN,1,IBORN,3,5*IBORN,5,1.D-4,RES)
         RESBORN = MS**2/ALPHAS**2/CONV*S1
         DEVBORN = MS**2/ALPHAS**2/CONV*S2
         
         CALL INTEG(SOVI,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESSOVI = MS**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVSOVI = MS**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELSOVI = DABS(DEVSOVI/RESSOVI)
         
         CALL INTEG(SOVISCA,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESSOVISCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVSOVISCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2
         
         CALL INTEG(HARD,3,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARD = MS**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVHARD = MS**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELHARD = DABS(DEVHARD/RESHARD)
         
         CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARDSCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVHARDSCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2
      ELSE
         CALL INTEG(HARD,3,IHARD,5,5*IHARD,10,5.D-4,RES)
         RESHARD = MS**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVHARD = MS**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELHARD = DABS(DEVHARD/RESHARD)
         
         CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARDSCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVHARDSCA = MS**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

         RESBORN = 0D0
         RESSOVI = 0D0
         RELSOVI = 0D0
         RESSOVISCA = 0D0

      ENDIF

      RESSCALE = RESSOVISCA +RESHARDSCA

c$$$      PRINT 50,ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE
c$$$      WRITE(10,50)ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE

      PRINT 51,flav(iflavor),ETA,
     +     RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE

 50   FORMAT(G10.4,' ',G11.5,'  ',G11.5,' ',G8.2,'  ',G11.5,' ',
     +     G8.2,'  ',G11.5)

 51   FORMAT(A5,G9.3,' ',G11.5,'  ',G11.5,' ',F5.4,'  ',G11.5,' ',
     +     F5.4,'  ',G11.5)


 100  CONTINUE

      END

C ======================================================================

      DOUBLE PRECISION FUNCTION BORN(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(1)
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MS**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) BORN = SB*DSSQQB(ALPHAS,S,T1,MS,MG,1)
      IF (IFLAVOR.EQ.2) BORN = SB*DSSQQB(ALPHAS,S,T1,MS,MG,0)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION SOVI(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(1)
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MS**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) SOVI = SB*DSSQQV(ALPHAS,S,T1,MS,MG,MT,1)
      IF (IFLAVOR.EQ.2) SOVI = SB*DSSQQV(ALPHAS,S,T1,MS,MG,MT,0)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION HARD(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(3)
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      Y = VAR(2)
      Z = VAR(3)

      SB = S * DSQRT(1.D0 -4.D0*MS**2/S)
      T1 = DSQRT((DEL -S)**2 -4*MS**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4MAX = S +T1 +MS**2*S/T1
      S4 = ( S4MAX -DEL)*Y +DEL
      PREF = S4MAX*SB

C***  EPS = ( GAMMA(GLUINO)/MASS(GLUINO) )**2
      EPS = 1.D-8

      IF (IFLAVOR.EQ.3) THEN
         M2 = MG**2 -MS**2
         S4MAX0 = S - 2.D0*DSQRT(S*MS**2)
         GAM = DSQRT(EPS)*MG**2
         TYMIN = DATAN(-M2/GAM)
         TYMAX = DATAN((S4MAX0-M2)/GAM)
         TY = (TYMAX -TYMIN) * X + TYMIN 
         PREFTY = (TYMAX -TYMIN)/GAM 
         S40 = GAM*DTAN(TY) + M2
         S41 = M2

         PRET10 = DSQRT( (S-S40)**2 -4*MS**2*S )
         PRET11 = 1
         IF((S-S41)**2-4*MS**2*S.GT.0.D0)THEN
          PRET11 = DSQRT( (S-S41)**2 -4*MS**2*S )
         ENDIF
         T10 = -(S-S40)/2.D0 +(Y -0.5D0)*PRET10
         T11 = -(S-S41)/2.D0 +(Y -0.5D0)*PRET11

         PREF0 = PRET10 * PREFTY
         PREF1 = PRET11 * PREFTY
         S4G02 = (S40 -M2)**2 + EPS*MG**4

         S3MAX0 = S - 2.D0*DSQRT(S*MS**2)
         GAM = DSQRT(EPS)*MG**2
         TZMIN = DATAN(-M2/GAM)
         TZMAX = DATAN((S3MAX0-M2)/GAM)
         TZ = (TZMAX -TZMIN) * X + TZMIN 
         PREFTZ = (TZMAX -TZMIN)/GAM 
         S32 = GAM*DTAN(TZ) + M2
         S33 = M2

         PRES42 = S32/(S32+MS**2)*DSQRT((S-S32)**2 -4*MS**2*S)
         PRES43 = 1
         IF((S-S33)**2 -4*MS**2*S.GT.0.D0)THEN
          PRES43 = S33/(S33+MS**2)*DSQRT((S-S33)**2 -4*MS**2*S)
         ENDIF
         S42 = S32/2D0/(S32+MS**2)*(S -S32 -2*MS**2) +PRES42 *(Y -0.5D0)
         S43 = S33/2D0/(S33+MS**2)*(S -S33 -2*MS**2) +PRES43 *(Y -0.5D0)
         PRET12 = DSQRT( (S-S42)**2 -4*MS**2*S )
         PRET13 = 1
         IF((S-S43)**2 -4*MS**2*S.GT.0.D0)THEN
          PRET13 = DSQRT( (S-S43)**2 -4*MS**2*S )
         ENDIF
         T12 = -(S-S42)/2.D0 +(Z -0.5D0)*PRET12
         T13 = -(S-S43)/2.D0 +(Z -0.5D0)*PRET13
         PREX2 = S42/2D0/(S42+MS**2)*DSQRT((S-S42)**2-4*MS**2*S)
         PREX3 = 1
         IF((S-S43)**2 -4*MS**2*S.GT.0.D0)THEN
          PREX3 = S43/2D0/(S43+MS**2)*DSQRT((S-S43)**2-4*MS**2*S)
         ENDIF

         PREF2 = PRES42 *PRET12 /PREX2 * PREFTZ
         PREF3 = PRES43 *PRET13 /PREX3 * PREFTZ
         
         S3G2  = S32 -M2
         S3G3  = 0D0

         SB4 = S * DSQRT(1.D0 -4.D0*MS**2/S)
         T14 = DSQRT(S**2 -4*MS**2*S)*(X -1.D0/2.) -S/2.D0
         S4MAX4 = S +T14 +MS**2*S/T14
         PREF4 = SB4 

      END IF

      IF (IFLAVOR.EQ.1) HARD = 
     +     + PREF*DSSQQH(ALPHAS,S,T1,S4,MS,MG,1)
     +     + PREF*DSSQQD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,1)/(S4MAX-DEL)
      IF (IFLAVOR.EQ.2) HARD = 
     +     + PREF*DSSQQH(ALPHAS,S,T1,S4,MS,MG,0)
     +     + PREF*DSSQQD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,0)/(S4MAX-DEL)

      IF (IFLAVOR.EQ.3) THEN
         HARD1 = PREF0 * DSSGQH(ALPHAS,S,T10,S40,MS,MG,EPS)*S4G02
         IF ((MG.GT.MS).AND.(S.GT.(MS+MG)**2)) THEN
            HARD2=(PREF0 *DSSGQS(ALPHAS,S,T10,S40,MS,MG,EPS)
     +           - PREF1 *DSSGQS(ALPHAS,S,T11,M2 ,MS,MG,EPS) ) 
            HARD3 = 
     +           + PREF2 *DSSGQT(ALPHAS,S,T12,S42,S3G2,MS,MG,EPS)
     +           - PREF3 *DSSGQT(ALPHAS,S,T13,S43,0.D0,MS,MG,EPS)
         ELSE
            HARD2= PREF0 *DSSGQS(ALPHAS,S,T10,S40,MS,MG,EPS)
            HARD3= PREF2 *DSSGQT(ALPHAS,S,T12,S42,S3G2,MS,MG,EPS)
         END IF
         IF ((S4MAX4.GE.M2).AND.(MG.GT.MS).AND.(S.GT.(MS+MG)**2)) THEN
            HARD4 =PREF4*DSSGQU(ALPHAS,S,T14,+M2,MS,MG,EPS)
         ELSE
            HARD4 = 0D0
         END IF

         HARD = HARD1 +HARD2 +HARD3 +HARD4

      END IF

      IF (.NOT.(DABS(HARD).LT.1.D40)) THEN
         HARD = 0.D0
         PRINT *,'NAN'
      END IF

      RETURN
      END

C ======================================================================

    
      DOUBLE PRECISION FUNCTION SOVISCA(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(1)
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MS**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) SOVISCA = SB*DSSQQ1(ALPHAS,S,T1,MS,MG,SCA,1)
      IF (IFLAVOR.EQ.2) SOVISCA = SB*DSSQQ1(ALPHAS,S,T1,MS,MG,SCA,0)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION HARDSCA(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(2)
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      Y = VAR(2)
      SB = S * DSQRT(1.D0 -4.D0*MS**2/S)
      T1 = DSQRT((DEL -S)**2 -4*MS**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4MAX = S +T1 +MS**2*S/T1
      S4 = ( S4MAX -DEL)*Y +DEL
      PREF = (S4MAX -DEL)*(SB -S/SB*DEL)
      IF (IFLAVOR.EQ.1) HARDSCA = 
     +   + PREF*DSSQQ3(ALPHAS,S,T1,S4,MS,MG,SCA,1)
     +   + PREF*DSSQQ2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA,1)/(S4MAX-DEL)
      IF (IFLAVOR.EQ.2) HARDSCA = 
     +   + PREF*DSSQQ3(ALPHAS,S,T1,S4,MS,MG,SCA,0)
     +   + PREF*DSSQQ2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA,0)/(S4MAX-DEL)
      IF (IFLAVOR.EQ.3) HARDSCA = 
     +   + PREF*DSSGQ3(ALPHAS,S,T1,S4,MS,MG,SCA)

      RETURN
      END

C ======================================================================

