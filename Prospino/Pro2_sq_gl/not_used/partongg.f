      PROGRAM PARTONGG
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      character*5 flav(3)

      data flav/'gg   ','qqb  ','gqb  '/

      COMMON/IOUT/IPRINT
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/RESULT/S1,S2,S3,S4
      COMMON/IFLAVOR/IFLAVOR

      EXTERNAL BORN, SOVI, HARD, SOVISCA, HARDSCA

      CALL RSTART(12,34,56,78)

C --- PARAMETER TO CHANGE  -------------------------------------

C***  1:G G     2:Q QBAR    3:G QBAR

      DO 100 IFLAVOR = 1,3

C***  THE MASSES:

      MG = 200.D0      
      MS = 280.D0 
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

c$$$      IF (IFLAVOR.EQ.1) OPEN(10,FILE="partongg_gg.dat")
c$$$      IF (IFLAVOR.EQ.2) OPEN(10,FILE="partongg_qb.dat")
c$$$      IF (IFLAVOR.EQ.3) OPEN(10,FILE="partongg_gb.dat")

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
         
         S = 4.D0*MG**2*(ETA + 1.D0)
         DEL = 1.D-6*MG**2*(1 -4.D0*MG**2/S)
         
         IF (IFLAVOR.LE.2) THEN
           CALL INTEG(BORN,1,IBORN,3,5*IBORN,5,1.D-4,RES)
           RESBORN = MG**2/ALPHAS**2/CONV*S1
           DEVBORN = MG**2/ALPHAS**2/CONV*S2
           
           CALL INTEG(SOVI,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESSOVI = MG**2/ALPHAS**3/4.D0/PI/CONV*S1
           DEVSOVI = MG**2/ALPHAS**3/4.D0/PI/CONV*S2
           RELSOVI = DABS(DEVSOVI/RESSOVI)
           
           CALL INTEG(SOVISCA,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESSOVISCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
           DEVSOVISCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2
           
           CALL INTEG(HARD,3,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESHARD = MG**2/ALPHAS**3/4.D0/PI/CONV*S1
           DEVHARD = MG**2/ALPHAS**3/4.D0/PI/CONV*S2
           RELHARD = DABS(DEVHARD/RESHARD)

           CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESHARDSCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
           DEVHARDSCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

        ELSE

           CALL INTEG(HARD,3,IHARD,5,5*IHARD,10,5.D-4,RES)
           RESHARD = MG**2/ALPHAS**3/4.D0/PI/CONV*S1
           DEVHARD = MG**2/ALPHAS**3/4.D0/PI/CONV*S2
           RELHARD = DABS(DEVHARD/RESHARD)

           CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESHARDSCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
           DEVHARDSCA = MG**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

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

 51   FORMAT(A5,G9.3,' ',G11.5,'  ',G11.5,' ',F5.4,'  ',G11.5,' ',
     +     F5.4,'  ',G11.5)


 50   FORMAT(G10.4,' ',G11.5,'  ',G11.5,' ',G8.2,'  ',G11.5,' ',
     +     G8.2,'  ',G11.5)

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
      SBG = S * DSQRT(1.D0 -4.D0*MG**2/S)
      TG = SBG*(X -1.D0/2.) -S/2.D0
      UG = -S -TG
      IF (IFLAVOR.EQ.1) BORN = SBG*DGGGGB(ALPHAS,S,TG,MS,MG,MT)
      IF (IFLAVOR.EQ.2) BORN = SBG*DGGQBB(ALPHAS,S,TG,MS,MG,MT)
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
      SBG = S * DSQRT(1.D0 -4.D0*MG**2/S)
      TG = SBG*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) SOVI = SBG*DGGGGV(ALPHAS,S,TG,MS,MG,MT)
      IF (IFLAVOR.EQ.2) SOVI = SBG*DGGQBV(ALPHAS,S,TG,MS,MG,MT)
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

      SBG = S * DSQRT(1.D0 -4.D0*MG**2/S)
      TG = DSQRT((DEL -S)**2 -4*MG**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4GMAX = S +TG +MG**2*S/TG
      S4G = ( S4GMAX -DEL)*Y +DEL
      PREF = S4GMAX*SBG

C***  EPS = ( GAMMA(SQUARK)/MASS(SQUARK) )**2
      EPS = 1.D-5

      IF (IFLAVOR.EQ.3) THEN      
         M2 = MG**2 -MS**2
         S4GMAX0 = S - 2.D0*DSQRT(S*MG**2)
         GAM = DSQRT(EPS)*MS**2
         TYMIN = DATAN(+M2/GAM)
         TYMAX = DATAN((S4GMAX0+M2)/GAM)
         TY = (TYMAX -TYMIN) * X + TYMIN 
         PREFTY = (TYMAX -TYMIN)/GAM 
         S4G0 = GAM*DTAN(TY) - M2
         S4G1 = - M2

         PRETG0 = DSQRT( (S-S4G0)**2 -4*MG**2*S )
         TG0 = -(S-S4G0)/2.D0 +(Y -0.5D0)*PRETG0
         PRETG1 = 1
         IF((S-S4G1)**2-4*MG**2*S.GT.0.D0)THEN
          PRETG1 = DSQRT( (S-S4G1)**2 -4*MG**2*S )
         ENDIF
         TG1 = -(S-S4G1)/2.D0 +(Y -0.5D0)*PRETG1

         PREF0 = PRETG0 * PREFTY
         PREF1 = PRETG1 * PREFTY
         S40  = S4G0 +M2
         S402 = (S4G0+m2)**2 +EPS*MS**4

         S3gMAX0 = S - 2.D0*DSQRT(S*Mg**2)
         GAM = DSQRT(EPS)*Ms**2
         TZMIN = DATAN(+M2/GAM)
         TZMAX = DATAN((S3gMAX0+M2)/GAM)
         TZ = (TZMAX -TZMIN) * X + TZMIN 
         PREFTZ = (TZMAX -TZMIN)/GAM 
         S3g2 = GAM*DTAN(TZ) - M2
         S3g3 = - M2

         PRES4G2 = S3G2/(S3G2+MG**2)*DSQRT((S-S3G2)**2 -4*MG**2*S)
         PRES4G3 = S3G3/(S3G3+MG**2)*DSQRT((S-S3G3)**2 -4*MG**2*S)
         S4G2 = S3G2/2D0/(S3G2+MG**2)*
     +        (S -S3G2 -2*MG**2) +PRES4G2 *(Y -0.5D0)
         S4G3 = S3G3/2D0/(S3G3+MG**2)*
     +        (S -S3G3 -2*MG**2) +PRES4G3 *(Y -0.5D0)
         PRETG2 = DSQRT( (S-S4G2)**2 -4*MG**2*S )
         PRETG3 = 1
         IF((S-S4G3)**2 -4*MG**2*S.GT.0.D0)THEN
          PRETG3 = DSQRT( (S-S4G3)**2 -4*MG**2*S )
         ENDIF
         TG2 = -(S-S4G2)/2.D0 +(Z -0.5D0)*PRETG2
         TG3 = -(S-S4G3)/2.D0 +(Z -0.5D0)*PRETG3
         PREX2 = S4G2/2D0/(S4G2+MG**2)*DSQRT((S-S4G2)**2-4*MG**2*S)
         PREX3 = 1
         IF((S-S4G3)**2 -4*MG**2*S.GT.0.D0)THEN
          PREX3 = S4G3/2D0/(S4G3+MG**2)*DSQRT((S-S4G3)**2-4*MG**2*S)
         ENDIF

         PREF2 = PRES4G2 *PRETG2 /PREX2 * PREFTZ
         PREF3 = PRES4G3 *PRETG3 /PREX3 * PREFTZ
         
         S32  = S3G2 +M2
         S33  = 0D0

         SBG4 = S * DSQRT(1.D0 -4.D0*MG**2/S)
         TG4 = DSQRT(S**2 -4*MG**2*S)*(X -1.D0/2.) -S/2.D0
         S4GMAX4 = S +TG4 +MG**2*S/TG4
         PREF4 = SBG4 

      END IF

      IF (IFLAVOR.EQ.1)  HARD = PREF* (
     +     DGGGGH(ALPHAS,S,TG,S4G,MS,MG)
     +     + DGGGGD(ALPHAS,S,TG,S4G,MS,MG,DEL,S4GMAX)/(S4GMAX-DEL) )
      IF (IFLAVOR.EQ.2)  HARD = PREF* (
     +     DGGQBH(ALPHAS,S,TG,S4G,MS,MG,EPS)
     +     + DGGQBD(ALPHAS,S,TG,S4G,MS,MG,DEL,S4GMAX)/(S4GMAX-DEL) )
      IF (IFLAVOR.EQ.3)  THEN
         HARD1 = PREF0*DGGGBH(ALPHAS,S,TG0,S4G0,MS,MG,EPS)*S402
         IF ((MS.GT.MG).AND.(S.GT.(MS+MG)**2)) THEN
            HARD2 = 
     +           + PREF0*DGGGBS(ALPHAS,S,TG0,S4G0,MS,MG,EPS)
     +           - PREF1*DGGGBS(ALPHAS,S,TG1,S4G1,MS,MG,EPS)
            HARD3 = 
     +           + PREF2*DGGGBT(ALPHAS,S,TG2,S4G2,S32,MS,MG,EPS)
     +           - PREF3*DGGGBT(ALPHAS,S,TG3,S4G3,S33,MS,MG,EPS)
         ELSE
            HARD2 = PREF0*DGGGBS(ALPHAS,S,TG0,S4G0,MS,MG,EPS)
            HARD3 = PREF2*DGGGBT(ALPHAS,S,TG2,S4G2,S32,MS,MG,EPS)
         END IF
         IF ((S4GMAX4.GE.-M2).AND.(MS.GT.MG).AND.(S.GT.(MS+MG)**2)) THEN
            HARD4 =PREF4*DGGGBU(ALPHAS,S,TG4,-M2,MS,MG,EPS)
         ELSE
            HARD4 = 0D0
         END IF
         HARD = HARD1 + HARD2 +HARD3 +HARD4
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
      SBG = S * DSQRT(1.D0 -4.D0*MG**2/S)
      TG = SBG*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) SOVISCA = SBG*DGGGG1(ALPHAS,S,TG,MS,MG,SCA)
      IF (IFLAVOR.EQ.2) SOVISCA = SBG*DGGQB1(ALPHAS,S,TG,MS,MG,SCA)
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
      SBG = S * DSQRT(1.D0 -4.D0*MG**2/S)
      TG = DSQRT((DEL -S)**2 -4*MG**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4GMAX = S +TG +MG**2*S/TG
      S4G = ( S4GMAX -DEL)*Y +DEL
      PREF = (S4GMAX -DEL)*(SBG -S/SBG*DEL)

      IF (IFLAVOR.EQ.1) HARDSCA = PREF* (
     +     + DGGGG2(ALPHAS,S,TG,S4G,MS,MG,DEL,S4GMAX,SCA)/(S4GMAX -DEL)
     +     + DGGGG3(ALPHAS,S,TG,S4G,MS,MG,SCA) )
      IF (IFLAVOR.EQ.2) HARDSCA = PREF* (
     +     + DGGQB2(ALPHAS,S,TG,S4G,MS,MG,DEL,S4GMAX,SCA)/(S4GMAX -DEL)
     +     + DGGQB3(ALPHAS,S,TG,S4G,MS,MG,SCA) )
      IF (IFLAVOR.EQ.3) HARDSCA = PREF* (
     +     + DGGGB3(ALPHAS,S,TG,S4G,MS,MG,SCA) )
      RETURN
      END

C ======================================================================

