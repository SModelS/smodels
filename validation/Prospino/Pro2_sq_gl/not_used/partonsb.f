      PROGRAM PARTONSB
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      character*5 flav(5)

      data flav/'gg   ','qqb  ','qqbp ','gqb  ','qg   '/

      COMMON/IOUT/IPRINT
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/RESULT/S1,S2,S3,S4
      COMMON/IFLAVOR/IFLAVOR

      EXTERNAL BORN, SOVI, HARD, SOVISCA, HARDSCA

C --- PARAMETER TO CHANGE  -------------------------------------

C***  1:G G     2:Q QBAR     3:QP QBAR    4:G QBAR     5:Q G

      DO 100 IFLAVOR = 1,5
         

C***  THE MASSES:

      Ms = 280.D0      
      Mg = 200.D0
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

      CALL RSTART(12,34,56,78)

C -------------------------------------------------------------

c$$$      IF (IFLAVOR.EQ.1) OPEN(10,FILE="partonsb_gg.dat")
c$$$      IF (IFLAVOR.EQ.2) OPEN(10,FILE="partonsb_qb.dat")
c$$$      IF (IFLAVOR.EQ.3) OPEN(10,FILE="partonsb_qbp.dat")
c$$$      IF (IFLAVOR.EQ.4) OPEN(10,FILE="partonsb_gb.dat")
c$$$      IF (IFLAVOR.EQ.5) OPEN(10,FILE="partonsb_qg.dat")

      PRINT 20, MS,MG,MT
      PRINT 30
c$$$      WRITE(10, 20) MS,MG,MT
c$$$      WRITE(10, 30)

 20   FORMAT(' MS = ',G10.4,' MG= ',G10.4,' MT = ',G10.4)
 30   FORMAT('INIT  ETA        BORN          SOVI       REL     ',
     +   '   HARD       REL       SCALE')

C --- INTEGRATION BY VEGAS ------------------------------------


      DO 100  I = 1,INUM

        ETA = 10.D0**(ETALOG)
        ETALOG = ETALOG +STEP

        S = 4.D0*MS**2*(ETA + 1.D0)
        DEL = 1.D-4*MS**2 * (1 -4.D0*MS**2/S)

      IF (IFLAVOR.LE.3) THEN
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

      beta = dsqrt(1 - 4*ms**2/s)
      pi = 4D0*DATAN(1D0)

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
      U1 = -S -T1
      IF (IFLAVOR.EQ.1) BORN = SB*DSBGGB(ALPHAS,S,T1,MS,MG)
      IF (IFLAVOR.EQ.2) BORN = SB*DSBQBB(ALPHAS,S,T1,MS,MG,1)
      IF (IFLAVOR.EQ.3) BORN = SB*DSBQBB(ALPHAS,S,T1,MS,MG,0)
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
      IF (IFLAVOR.EQ.1) SOVI = SB*DSBGGV(ALPHAS,S,T1,MS,MG,MT)
      IF (IFLAVOR.EQ.2) SOVI = SB*DSBQBV(ALPHAS,S,T1,MS,MG,MT,1)
      IF (IFLAVOR.EQ.3) SOVI = SB*DSBQBV(ALPHAS,S,T1,MS,MG,MT,0)
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
      EPS = 1.D-5

      IF (IFLAVOR.EQ.4) THEN
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

      END IF

      IF (IFLAVOR.EQ.5) THEN
         M2 = MG**2 -MS**2
         S3MAX0 = S - 2.D0*DSQRT(S*MS**2)

         GAM = DSQRT(EPS)*MG**2
         TYMIN = DATAN(-M2/GAM)
         TYMAX = DATAN((S3MAX0-M2)/GAM)
         TY = (TYMAX -TYMIN) * X + TYMIN 
         PREFTY = (TYMAX -TYMIN)/GAM 
         S30 = GAM*DTAN(TY) + M2
         S31 = M2

         PRES40 = S30/(S30+MS**2)*DSQRT((S-S30)**2 -4*MS**2*S)
         PRES41 = S31/(S31+MS**2)*DSQRT((S-S31)**2 -4*MS**2*S)
         S40 = S30/2D0/(S30+MS**2)*(S -S30 -2*MS**2) +PRES40 *(Y -0.5D0)
         S41 = S31/2D0/(S31+MS**2)*(S -S31 -2*MS**2) +PRES41 *(Y -0.5D0)
         PRET10 = DSQRT( (S-S40)**2 -4*MS**2*S )
         PRET11 = 1
         IF((S-S41)**2 -4*MS**2*S.GT.0.D0)THEN
          PRET11 = DSQRT( (S-S41)**2 -4*MS**2*S )
         ENDIF
         T10 = -(S-S40)/2.D0 +(Z -0.5D0)*PRET10
         T11 = -(S-S41)/2.D0 +(Z -0.5D0)*PRET11
         PREX0 = S40/2D0/(S40+MS**2)*DSQRT((S-S40)**2-4*MS**2*S)
         PREX1 = 1
         IF((S-S41)**2 -4*MS**2*S.GT.0.D0)THEN
          PREX1 = S41/2D0/(S41+MS**2)*DSQRT((S-S41)**2-4*MS**2*S)
         ENDIF

         PREF0 = PRES40 *PRET10 /PREX0 * PREFTY
         PREF1 = PRES41 *PRET11 /PREX1 * PREFTY
         
         S3G0  = S30 -M2
         S3G1  = 0D0
      END IF


      IF (IFLAVOR.EQ.1)  HARD = PREF*
     +     ( DSBGGH(ALPHAS,S,T1,S4,MS,MG)
     +      + DSBGGD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX)/(S4MAX-DEL) )
      IF (IFLAVOR.EQ.2)  HARD = PREF* (
     +      + DSBQBH(ALPHAS,S,T1,S4,MS,MG,1,EPS)
     +      + DSBQBD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,1)/(S4MAX-DEL) )
      IF (IFLAVOR.EQ.3)  HARD = PREF* (
     +      + DSBQBH(ALPHAS,S,T1,S4,MS,MG,0,EPS)
     +      + DSBQBD(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,0)/(S4MAX-DEL) )

      IF (IFLAVOR.EQ.4) THEN
         HARD1 = PREF0 * DSBGBH(ALPHAS,S,T10,S40,MS,MG,EPS)*S4G02
         IF ((MG.GT.MS).AND.(S.GT.(MS+MG)**2)) THEN
            HARD2 = 
     +           +(PREF0 * DSBGBS(ALPHAS,S,T10,S40,MS,MG,EPS)
     +           - PREF1 * DSBGBS(ALPHAS,S,T11,S41,MS,MG,EPS) )
         ELSE
            HARD2 = 
     +           + PREF0 * DSBGBS(ALPHAS,S,T10,S40,MS,MG,EPS)
         END IF
         HARD = HARD1 + HARD2
      END IF

      IF (IFLAVOR.EQ.5)  THEN
         HARD1 = +PREF * DSBQGH(ALPHAS,S,T1,S4,MS,MG,EPS)
         IF ((MG.GT.MS).AND.(S.GT.(MS+MG)**2)) THEN
            HARD3 = 
     +           +PREF0 * DSBQGT(ALPHAS,S,T10,S40,S3G0,MS,MG,EPS)
     +           -PREF1 * DSBQGT(ALPHAS,S,T11,S41,S3G1,MS,MG,EPS)
         ELSE
            HARD3 = 
     +           +PREF0 * DSBQGT(ALPHAS,S,T10,S40,S3G0,MS,MG,EPS)
         END IF
         HARD = HARD1 +HARD3
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
      IF (IFLAVOR.EQ.1) SOVISCA = SB*DSBGG1(ALPHAS,S,T1,MS,MG,SCA)
      IF (IFLAVOR.EQ.2) SOVISCA = SB*DSBQB1(ALPHAS,S,T1,MS,MG,SCA,1)
      IF (IFLAVOR.EQ.3) SOVISCA = SB*DSBQB1(ALPHAS,S,T1,MS,MG,SCA,0)
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

      IF (IFLAVOR.EQ.1) HARDSCA = PREF* (
     +     + DSBGG2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA)/(S4MAX -DEL)
     +     + DSBGG3(ALPHAS,S,T1,S4,MS,MG,SCA) )
      IF (IFLAVOR.EQ.2) HARDSCA = PREF* (
     +     + DSBQB2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA,1)/(S4MAX -DEL)
     +     + DSBQB3(ALPHAS,S,T1,S4,MS,MG,SCA,1) )
      IF (IFLAVOR.EQ.3) HARDSCA = PREF* (
     +     + DSBQB2(ALPHAS,S,T1,S4,MS,MG,DEL,S4MAX,SCA,0)/(S4MAX -DEL)
     +     + DSBQB3(ALPHAS,S,T1,S4,MS,MG,SCA,0) )
      IF (IFLAVOR.EQ.4) HARDSCA = PREF* (
     +     + DSBGB3(ALPHAS,S,T1,S4,MS,MG,SCA) )
      IF (IFLAVOR.EQ.5) HARDSCA = PREF* (
     +     + DSBQG3(ALPHAS,S,T1,S4,MS,MG,SCA) )

      RETURN
      END

C ======================================================================

