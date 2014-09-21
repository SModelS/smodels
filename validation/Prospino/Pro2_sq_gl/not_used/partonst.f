      PROGRAM PARTONSTOP
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)

      COMMON/IOUT/IPRINT
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/DELTA/DEL
      COMMON/RESULT/S1,S2,S3,S4
      COMMON/IFLAVOR/IFLAVOR

      EXTERNAL BORN, SOVI, HARD, SOVISCA, HARDSCA

C --- PARAMETER TO CHANGE  -------------------------------------

C***  1:G G     2:Q QBAR   3:G QBAR     4:Q G

      DO 100 IFLAVOR = 1,4
         

C***  THE MASSES:

      Mst1 = 90.D0      
      Mg = 200.D0
      MT = 175.D0
      MS = 280D0 +10D0
      MST2 = MS + 20D0
      SINT = -0.3D0

c$$$      Mst1 = 280.D0      
c$$$      Mg = 200.D0
c$$$      MT = 10.D0
c$$$      MS = MST1 
c$$$      MST2 = MST1 
c$$$      SINT = 0D0


C***  START VALUE IN LOG(ETA), NUMBERS OF STEPS, STEP WIDTH IN LOG(ETA)
      
      ETALOG = -2D0
      INUM = 3
      STEP = 2D0


C***  THE NUMBER OF VEGAS CALLS ( DEFAULT = 500 )

      IBORN = 500
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

      IF (IFLAVOR.EQ.1) OPEN(10,FILE="partonstop_gg.dat")
      IF (IFLAVOR.EQ.2) OPEN(10,FILE="partonstop_qb.dat")
      IF (IFLAVOR.EQ.3) OPEN(10,FILE="partonstop_gb.dat")
      IF (IFLAVOR.EQ.4) OPEN(10,FILE="partonstop_qg.dat")

      PRINT 20, MST1,MG,MT,ms,mst2,sint
      PRINT 30
      WRITE(10, 20) MST1,MG,MT,ms,mst2,sint
      WRITE(10, 30)

 20   FORMAT(' MST1= ',F6.1,' MG= ',F6.1,' MT= ',F6.1, ' MS= ',F6.1,
     +     ' MST2= ',F6.1,' sint= ',F7.3)
 30   FORMAT('  ETA        BORN          SOVI       REL     ',
     +   '   HARD       REL       SCALE')

C --- INTEGRATION BY VEGAS ------------------------------------


      DO 100  I = 1,INUM

        ETA = 10.D0**(ETALOG)
        ETALOG = ETALOG +STEP

        S = 4.D0*MST1**2*(ETA + 1.D0)
        DEL = 1.D-4*MST1**2 * (1 -4.D0*MST1**2/S)

      IF (IFLAVOR.LE.2) THEN
         CALL INTEG(BORN,1,IBORN,3,5*IBORN,5,1.D-4,RES)
         RESBORN = MST1**2/ALPHAS**2/CONV*S1
         DEVBORN = MST1**2/ALPHAS**2/CONV*S2
         
         CALL INTEG(SOVI,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESSOVI = MST1**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVSOVI = MST1**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELSOVI = DABS(DEVSOVI/RESSOVI)
         
         CALL INTEG(SOVISCA,1,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESSOVISCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVSOVISCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2
         
         CALL INTEG(HARD,3,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARD = MST1**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVHARD = MST1**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELHARD = DABS(DEVHARD/RESHARD)
         
         CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARDSCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVHARDSCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2
      ELSE
         CALL INTEG(HARD,3,IHARD,5,5*IHARD,10,5.D-4,RES)
         RESHARD = MST1**2/ALPHAS**3/4.D0/PI/CONV*S1
         DEVHARD = MST1**2/ALPHAS**3/4.D0/PI/CONV*S2
         RELHARD = DABS(DEVHARD/RESHARD)
         
         CALL INTEG(HARDSCA,2,ISOVI,3,5*ISOVI,5,5.D-4,RES)
         RESHARDSCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S1
         DEVHARDSCA = MST1**2/ALPHAS**3/4.D0/PI/CONV/DLOG(SCA)*S2

         RESBORN = 0D0
         RESSOVI = 0D0
         RELSOVI = 0D0
         RESSOVISCA = 0D0

      ENDIF

      RESSCALE = RESSOVISCA +RESHARDSCA

      PRINT 50,ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE
      WRITE(10,50)ETA,RESBORN,RESSOVI,RELSOVI,RESHARD,RELHARD, RESSCALE


      beta = dsqrt(1 - 4*MST1**2/s)
      pi = 4D0*DATAN(1D0)

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
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MST1**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      U1 = -S -T1
      IF (IFLAVOR.EQ.1) BORN = SB*DSTGGB(ALPHAS,S,T1,MST1,MG)
      IF (IFLAVOR.EQ.2) BORN = SB*DSTQBB(ALPHAS,S,T1,MST1,MG)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION SOVI(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(1)
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MST1**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) 
     +     SOVI = SB*DSTGGV(ALPHAS,S,T1,MST1,MG,MT,MS,MST2,SINT)
      IF (IFLAVOR.EQ.2) 
     +     SOVI = SB*DSTQBV(ALPHAS,S,T1,MST1,MG,MT,MS,MST2,SINT)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION HARD(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(3)
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/DELTA/DEL
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      Y = VAR(2)
      Z = VAR(3)

      SB = S * DSQRT(1.D0 -4.D0*MST1**2/S)
      T1 = DSQRT((DEL -S)**2 -4*MST1**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4MAX = S +T1 +MST1**2*S/T1
      S4 = ( S4MAX -DEL)*Y +DEL
      PREF = S4MAX*SB

      IF (IFLAVOR.EQ.1)  HARD = PREF*
     +     ( DSTGGH(ALPHAS,S,T1,S4,MST1,MG)
     +      + DSTGGD(ALPHAS,S,T1,S4,MST1,MG,DEL,S4MAX)/(S4MAX-DEL) )
      IF (IFLAVOR.EQ.2)  HARD = PREF* (
     +      + DSTQBH(ALPHAS,S,T1,S4,MST1,MG)
     +      + DSTQBD(ALPHAS,S,T1,S4,MST1,MG,DEL,S4MAX)/(S4MAX-DEL) )
      IF (IFLAVOR.EQ.3) HARD = PREF * 
     +     DSTGBH(ALPHAS,S,T1,S4,MST1,MG)
      IF (IFLAVOR.EQ.4) HARD = PREF * 
     +     DSTQGH(ALPHAS,S,T1,S4,MST1,MG)

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
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      SB = S * DSQRT(1.D0 -4.D0*MST1**2/S)
      T1 = SB*(X -1.D0/2.) -S/2.D0
      IF (IFLAVOR.EQ.1) SOVISCA = SB*DSTGG1(ALPHAS,S,T1,MST1,MG,SCA)
      IF (IFLAVOR.EQ.2) SOVISCA = SB*DSTQB1(ALPHAS,S,T1,MST1,MG,SCA)
      RETURN
      END

C ======================================================================
    
      DOUBLE PRECISION FUNCTION HARDSCA(VAR)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      IMPLICIT COMPLEX*16 (K)
      DIMENSION VAR(2)
      COMMON/CONST/S,ALPHAS,MST1,MG,MT,MS,MST2,SINT,SCA
      COMMON/DELTA/DEL
      COMMON/IFLAVOR/IFLAVOR
      X = VAR(1)
      Y = VAR(2)
      SB = S * DSQRT(1.D0 -4.D0*MST1**2/S)
      T1 = DSQRT((DEL -S)**2 -4*MST1**2*S)*(X -1.D0/2.) -(S -DEL)/2.D0
      S4MAX = S +T1 +MST1**2*S/T1
      S4 = ( S4MAX -DEL)*Y +DEL
      PREF = (S4MAX -DEL)*(SB -S/SB*DEL)

      IF (IFLAVOR.EQ.1) HARDSCA = PREF* (
     +     + DSTGG2(ALPHAS,S,T1,S4,MST1,MG,DEL,S4MAX,SCA)/(S4MAX -DEL)
     +     + DSTGG3(ALPHAS,S,T1,S4,MST1,MG,SCA) )
      IF (IFLAVOR.EQ.2) HARDSCA = PREF* (
     +     + DSTQB2(ALPHAS,S,T1,S4,MST1,MG,DEL,S4MAX,SCA)/(S4MAX -DEL)
     +     + DSTQB3(ALPHAS,S,T1,S4,MST1,MG,SCA) )
      IF (IFLAVOR.EQ.3) HARDSCA = PREF* (
     +     + DSTGB3(ALPHAS,S,T1,S4,MST1,MG,SCA) )
      IF (IFLAVOR.EQ.4) HARDSCA = PREF* (
     +     + DSTQG3(ALPHAS,S,T1,S4,MST1,MG,SCA) )

      RETURN
      END

C ======================================================================
