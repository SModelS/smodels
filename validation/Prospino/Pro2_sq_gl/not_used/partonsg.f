      PROGRAM PARTONSG
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      character*5 flav(6)

      COMMON/IOUT/IPRINT
      COMMON/CONST/S,ALPHAS,MS,MG,MT,SCA
      COMMON/DELTA/DEL
      COMMON/RESULT/S1,S2,S3,S4
      COMMON/IFLAVOR/IFLAVOR

      data flav/'qg   ','gg   ','qqb  ','qqbp ','qq   ','qqp  '/

      EXTERNAL BORN, SOVI, HARD, SOVISCA, HARDSCA

      CALL RSTART(12,34,56,78)

C --- PARAMETER TO CHANGE  -------------------------------------

C***  1:Q G   2:G G   3:Q QBAR   4:QP QBAR   5:Q Q   6:Q QP

      DO 100 IFLAVOR = 2,2
         

C***  THE MASSES:

      MS = 280.D0 
      MG = 200.D0      
      MT = 175.D0

C***  START VALUE IN LOG(ETA), NUMBERS OF STEPS, STEP WIDTH IN LOG(ETA)
      
      ETALOG = -3D0
      INUM = 31
      STEP = 0.2D0

      ETALOG = -0.54D0
      INUM = 20
      STEP = 0.01D0


C***  THE NUMBER OF VEGAS CALLS ( DEFAULT = 500 )

      IBORN = 1000
      ISOVI = 500
      IHARD = 500

C***  PRINT VEGAS STATISTICS OR NOT

      IPRINT = 0


C --- CONSTANTS DON'T CHANGE PLEASE ------------------------------

      PI = 4.D0*DATAN(1.D0)
      CONV = 389379660.D0
      ALPHAS = 0.1D0
      SCA = 4.D0


C -------------------------------------------------------------

c$$$      IF (IFLAVOR.EQ.1) OPEN(10,FILE="partonsg_qg.dat")
c$$$      IF (IFLAVOR.EQ.2) OPEN(10,FILE="partonsg_gg.dat")
c$$$      IF (IFLAVOR.EQ.3) OPEN(10,FILE="partonsg_qb.dat")
c$$$      IF (IFLAVOR.EQ.4) OPEN(10,FILE="partonsg_qbp.dat")
c$$$      IF (IFLAVOR.EQ.5) OPEN(10,FILE="partonsg_qq.dat")
c$$$      IF (IFLAVOR.EQ.6) OPEN(10,FILE="partonsg_qp.dat")

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

        S = (MG +MS)**2*(ETA + 1.D0)
        DEL = 1.D-6*MG*MS*(1 -(MS+MG)**2/S)

        IF (IFLAVOR.EQ.1) THEN
           CALL INTEG(BORN,1,IBORN,3,5*IBORN,5,1.D-4,RES)
           RESBORN = (MS+MG)**2/4.D0/ALPHAS**2/CONV*S1
           DEVBORN = (MS+MG)**2/4.D0/ALPHAS**2/CONV*S2
           
           CALL INTEG(SOVI,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESSOVI = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S1
           DEVSOVI = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S2
           RELSOVI = DABS(DEVSOVI/RESSOVI)
           
           CALL INTEG(SOVISCA,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESSOVISCA=(MS+MG)**2/4D0/ALPHAS**3/4D0/PI/CONV/DLOG(SCA)*S1
           DEVSOVISCA=(MS+MG)**2/4D0/ALPHAS**3/4D0/PI/CONV/DLOG(SCA)*S2
           
           CALL INTEG(HARD,3,ISOVI,3,5*ISOVI,5,5.D-4,RES)
           RESHARD = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S1
           DEVHARD = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S2
           RELHARD = DABS(DEVHARD/RESHARD)

           CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
          RESHARDSCA=(MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
          DEVHARDSCA=(MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

        ELSE

           CALL INTEG(HARD,3,IHARD,5,5*IHARD,10,5.D-4,RES)
           RESHARD = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S1
           DEVHARD = (MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV*S2
           RELHARD = DABS(DEVHARD/RESHARD)
 
           CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
          RESHARDSCA=(MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
          DEVHARDSCA=(MS+MG)**2/4.D0/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

          RESBORN = 0D0
          RESSOVI = 0D0
          RELSOVI = 0D0
          RESSOVISCA = 0D0

      ENDIF

      RESSCALE = RESSOVISCA +RESHARDSCA

      PRINT 50,ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE
c$$$      WRITE(10,50)ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE
c$$$      PRINT 51,flav(iflavor),ETA,
c$$$     +     RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE

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
      XLAM = DSQRT((S -MG**2 -MS**2)**2 -4*MG**2*MS**2)
      TG = -S/2.D0*(1.D0 +(MG**2 -MS**2)/S -(2*X -1.D0)*XLAM/S)
      IF (IFLAVOR.EQ.1) BORN = XLAM*DSGQGB(ALPHAS,S,TG,MS,MG)
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
      XLAM = DSQRT((S -MG**2 -MS**2)**2 -4*MG**2*MS**2)
      TG = -S/2.D0*(1.D0 +(MG**2 -MS**2)/S -(2*X -1.D0)*XLAM/S)
      IF (IFLAVOR.EQ.1) SOVI = XLAM*DSGQGV(ALPHAS,S,TG,MS,MG,MT)
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
      IF (IFLAVOR.EQ.1) THEN
         XLAM = DSQRT((S -MG**2 -MS**2)**2 -4*MG**2*MS**2)
         TG = -(S +MG**2 -MS**2 -DEL)/2.D0 
     +        + (X -0.5D0) *DSQRT((S+MG**2-MS**2-DEL)**2 -4*MG**2*S)
         S4MAX = S +TG + MG**2 -MS**2 +MG**2*S/TG
         S4 = ( S4MAX -DEL)*Y +DEL
         PREF = S4MAX*XLAM
      ELSE
C***     EPSG = ( GAMMA(GLUINO)/MASS(GLUINO) )**2
C***     EPSS = ( GAMMA(SQUARK)/MASS(SQUARK) )**2
         EPSS = 1.D-10
         EPSG = 1.D-10
         M2 = MG**2 -MS**2
         S4MAX0 = S + M2 -2.D0 * DSQRT(S*MG**2)
         GAMG = DSQRT(EPSG)*MG**2
         TYMIN = DATAN(-M2/GAMG)
         TYMAX = DATAN((S4MAX0-M2)/GAMG)
         TY = (TYMAX -TYMIN) * X + TYMIN 
         PREFTY = (TYMAX -TYMIN)/GAMG 
         S40 = GAMG*DTAN(TY) + M2
         S4G0 = S40 -M2

         PRETG0 = DSQRT((S - S4G0)**2 - 4*S*MG**2)
         TG0 = - (S -S4G0)/2.D0 +(Y -0.5D0)*PRETG0
         PRETG1 = 1
         IF(S**2-4*MG**2*S.GT.0.D0)THEN
          PRETG1 = DSQRT(S**2 - 4*S*MG**2)
         ENDIF
         TG1 = - S/2.D0 +(Y -0.5D0)*PRETG1
         PREF0 = PREFTY * PRETG0
         PREF1 = PREFTY * PRETG1

         S3GMAX2 = S - M2 - 2*DSQRT(S*MS**2) 
         GAMS = DSQRT(EPSS)*MS**2
         TZMIN = DATAN(+M2/GAMS)
         TZMAX = DATAN((S3GMAX2+M2)/GAMS)
         TZ = (TZMAX -TZMIN) * X + TZMIN 
         PREFTZ = (TZMAX -TZMIN)/GAMS 
         S32 = GAMS*DTAN(TZ) 
         S33 = 0D0
         S3G2 = S32 -M2
         S3G3 = S33 -M2

         PRES42 = S3G2/(S32+MS**2)*DSQRT((S-S32)**2 -4*MS**2*S)
         PRES43 = 1
         IF((S-S33)**2 -4*MS**2*S.GT.0.D0)THEN
          PRES43 = S3G3/(S33+MS**2)*DSQRT((S-S33)**2 -4*MS**2*S)
         ENDIF
         S42 = S3G2/2D0/(S32+MS**2)*(S -S32 -2*MS**2) 
     +        +PRES42 *(Y -0.5D0)
         S43 = S3G3/2D0/(S33+MS**2)*(S -S33 -2*MS**2) 
     +        +PRES43 *(Y -0.5D0)
         S4G2 = S42 -M2
         S4G3 = S43 -M2

         PRETG2 = DSQRT( (S-S4G2)**2 -4*MG**2*S )
         PRETG3 = 1
         IF((S-S4G3)**2 -4*MG**2*S.GT.0.D0)THEN
          PRETG3 = DSQRT( (S-S4G3)**2 -4*MG**2*S )
         ENDIF
         TG2 = -(S-S4G2)/2.D0 +(Z -0.5D0)*PRETG2
         TG3 = -(S-S4G3)/2.D0 +(Z -0.5D0)*PRETG3

         PREX2 = S42/2D0/(S42+MS**2)*DSQRT((S-S4G2)**2-4*MG**2*S)
         PREX3 = 1
         IF((S-S4G3)**2-4*MG**2*S.GT.0.D0)THEN
          PREX3 = S43/2D0/(S43+MS**2)*DSQRT((S-S4G3)**2-4*MG**2*S)
         ENDIF

         PREF2 = PRES42 *PRETG2 /PREX2 * PREFTZ
         PREF3 = PRES43 *PRETG3 /PREX3 * PREFTZ
         
         S4G02 = S4G0**2 + EPSG*MG**4
      END IF

      IF (IFLAVOR.EQ.1) HARD = PREF* (
     +     + DSGQGH(ALPHAS,S,TG,S4,MS,MG)
     +     + DSGQGD(ALPHAS,S,TG,S4,MS,MG,DEL,S4MAX)/(S4MAX -DEL) )

      IF (IFLAVOR.EQ.2) THEN
         HARD1 = 0D0
     +        + PREF0 *DSGGGH(ALPHAS,S,TG0,S40,MS,MG,EPSS,EPSG) *S4G02
         IF ((S.GT.4*MG**2).AND.(MG.GT.MS)) THEN
            HARD2 = 0D0
     +           + PREF0 *DSGGGS(ALPHAS,S,TG0,S40,MS,MG)
     +           - PREF1 *DSGGGS(ALPHAS,S,TG1,M2 ,MS,MG) 
         ELSE
            HARD2 = 0D0
     +           + PREF0 *DSGGGS(ALPHAS,S,TG0,S40,MS,MG)
         END IF
         IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
            HARD3 = 0D0
     +           + PREF2*DSGGGT(ALPHAS,S,TG2,S42,S32,MS,MG)
     +           - PREF3*DSGGGT(ALPHAS,S,TG3,S43,0D0,MS,MG)
         ELSE
            HARD3 = 0D0
     +           + PREF2*DSGGGT(ALPHAS,S,TG2,S42,S32,MS,MG)
         END IF
         HARD = HARD1 + HARD2 +HARD3
      END IF

      IF (IFLAVOR.EQ.3) THEN
         HARD1 = 
     +        PREF0 *DSGQBH(ALPHAS,S,TG0,S40,MS,MG,1,EPSS,EPSG)*S4G02
         IF ((S.GT.4*MG**2).AND.(MG.GT.MS)) THEN
            HARD2 = 
     +           + PREF0 *DSGQBS(ALPHAS,S,TG0,S40,MS,MG,1)
     +           - PREF1 *DSGQBS(ALPHAS,S,TG1,M2 ,MS,MG,1)
         ELSE
            HARD2 =
     +           + PREF0 *DSGQBS(ALPHAS,S,TG0,S40,MS,MG,1)
         END IF
         IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
            HARD3 =
     +           + PREF2*DSGQBT(ALPHAS,S,TG2,S42,S32,MS,MG,1)
     +           - PREF3*DSGQBT(ALPHAS,S,TG3,S43,0D0,MS,MG,1)
         ELSE
            HARD3 =
     +           + PREF2*DSGQBT(ALPHAS,S,TG2,S42,S32,MS,MG,1)
         END IF
         HARD = HARD1 + HARD2 +HARD3
      END IF

      IF (IFLAVOR.EQ.4) THEN
         HARD1 = 
     +        PREF0* DSGQBH(ALPHAS,S,TG0,S40,MS,MG,0,EPSS,EPSG)*S4G02 
         HARD2 = 0.D0
         IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
            HARD3 =
     +           + PREF2*DSGQBT(ALPHAS,S,TG2,S42,S32,MS,MG,0)
     +           - PREF3*DSGQBT(ALPHAS,S,TG3,S43,0D0,MS,MG,0)
         ELSE
            HARD3 =
     +           + PREF2*DSGQBT(ALPHAS,S,TG2,S42,S32,MS,MG,0)
         END IF
         HARD = HARD1 + HARD2 +HARD3
      END IF

      IF (IFLAVOR.EQ.5) THEN
         HARD1 = 
     +        PREF0* DSGQQH(ALPHAS,S,TG0,S40,MS,MG,1,EPSS,EPSG) *S4G02
         HARD2 = 0.D0
         IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
            HARD3 =
     +           + PREF2*DSGQQT(ALPHAS,S,TG2,S42,S32,MS,MG,1)
     +           - PREF3*DSGQQT(ALPHAS,S,TG3,S43,0D0,MS,MG,1)
         ELSE
            HARD3 =
     +           + PREF2*DSGQQT(ALPHAS,S,TG2,S42,S32,MS,MG,1)
         END IF
         HARD = HARD1 + HARD2 +HARD3
      END IF

      IF (IFLAVOR.EQ.6) THEN
         HARD1 = 
     +        PREF0* DSGQQH(ALPHAS,S,TG0,S40,MS,MG,0,EPSS,EPSG) *S4G02
         HARD2 = 0.D0
         IF ((S.GT.4*MS**2).AND.(MS.GT.MG)) THEN
            HARD3 =
     +           + PREF2*DSGQQT(ALPHAS,S,TG2,S42,S32,MS,MG,0)
     +           - PREF3*DSGQQT(ALPHAS,S,TG3,S43,0D0,MS,MG,0)
         ELSE
            HARD3 =
     +           + PREF2*DSGQQT(ALPHAS,S,TG2,S42,S32,MS,MG,0)
         END IF
         HARD = HARD1 + HARD2 +HARD3
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
      XLAM = DSQRT((S -MG**2 -MS**2)**2 -4*MG**2*MS**2)
      TG = -S/2.D0*(1.D0 +(MG**2 -MS**2)/S -(2*X -1.D0)*XLAM/S)
      IF (IFLAVOR.EQ.1) SOVISCA = XLAM*DSGQG1(ALPHAS,S,TG,MS,MG,SCA)
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
      XLAM = DSQRT((S -MG**2 -MS**2)**2 -4*MG**2*MS**2)
      TG = -(S +MG**2 -MS**2 -DEL)/2.D0 
     +     + (X -0.5D0) *DSQRT((S+MG**2-MS**2-DEL)**2 -4*MG**2*S)
      S4MAX = S +TG + MG**2 -MS**2 +MG**2*S/TG
      S4 = ( S4MAX -DEL)*Y +DEL
      PREF = S4MAX*XLAM
      IF (IFLAVOR.EQ.1) HARDSCA = PREF* (
     +     + DSGQG2(ALPHAS,S,TG,S4,MS,MG,DEL,S4MAX,SCA)/(S4MAX -DEL)
     +     + DSGQG3(ALPHAS,S,TG,S4,MS,MG,SCA) )
      IF (IFLAVOR.EQ.2) HARDSCA = PREF* (
     +     + DSGGG3(ALPHAS,S,TG,S4,MS,MG,SCA) )
      IF (IFLAVOR.EQ.3) HARDSCA = PREF* (
     +     + DSGQB3(ALPHAS,S,TG,S4,MS,MG,SCA,1) )
      IF (IFLAVOR.EQ.4) HARDSCA = PREF* (
     +     + DSGQB3(ALPHAS,S,TG,S4,MS,MG,SCA,0) )
      IF (IFLAVOR.EQ.5) HARDSCA = PREF* (
     +     + DSGQQ3(ALPHAS,S,TG,S4,MS,MG,SCA,1) )
      IF (IFLAVOR.EQ.6) HARDSCA = PREF* (
     +     + DSGQQ3(ALPHAS,S,TG,S4,MS,MG,SCA,0) )

      RETURN
      END

C ======================================================================





 
