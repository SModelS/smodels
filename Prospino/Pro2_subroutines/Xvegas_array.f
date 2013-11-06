
C   25/02/93 404211828  MEMBER NAME  VEGASINT (MACROS)      FVS

************************************************************************
*                                                                      *
*     SUBROUTINE INTEG (CALLING VEGAS):                                *
*     ---------------------------------                                *
*                                                                      *
*     IDIM-DIMENSIONAL MONTE CARLO INTEGRATION OF FUNCTION "FUNCT"     *
*     1ST RUN (VEGAS) : IPOINT  CALLS, ITER  ITERATIONS                *
*     2ND RUN (VEGAS1): IPOINT1 CALLS, ITER1 ITERATIONS                *
*     RESULT "RES" IS RETURNED WITH DESIRED ACCURACY "ACC"             *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*     NEW VERSION WITH ARRAYS AS INPUT:                                *
*     1ST RUN (VEGAS) : IVEGAS(1) CALLS, IVEGAS(2) ITERATIONS          *
*     2ND RUN (VEGAS1) : IVEGAS(3) CALLS, IVEGAS(4) ITERATIONS         *
*     ACC FIXED TO 1.D-16                                              *
*                                                                      *
*                                                                      *
************************************************************************
      subroutine INTEG(fxn,idim,ivegas,ifast,res,rel)

      integer iprint, ivegas(1:4), ifast
      real*8  fxn,acc,res,rel

      integer INDO,IT,INXI
      real*8  SI,SI2,SWGT,SCHI,XI,SCALLS,D,DI
      COMMON/BVEG2/INDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     1,D(50,10),DI(50,10),INXI(50,10)

      real*8        S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4

      external fxn

      call RSTART(12,34,56,78)


      acc    = 1.d-4
      iprint = 10
      if (ifast .eq. 0) then 
         ipoint  = ivegas(1)
         iter    = ivegas(2)  
         ipoint1 = ivegas(3)
         iter1   = ivegas(4)
         call VEGAS(FXN,ACC,IDIM,IPOINT,ITER,IPRINT,0)
         call VEGAS1(FXN,ACC,IDIM,IPOINT1,ITER1,IPRINT,0)
      else if (ifast .eq. 1) then 
         ipoint1 = ivegas(3)
         iter1   = 3
         call VEGAS1(FXN,ACC,IDIM,IPOINT1,ITER1,IPRINT,0)
      end if 

      res = s1
      rel = abs(s2/s1)

      return
      end




      SUBROUTINE VEGAS(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 UNIV
      COMMON/BVEG2/INDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     1,D(50,10),DI(50,10),INXI(50,10)
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
c inserted by t. plehn 29.10.96 -> acc in VEGAS1      
      common/bveg3/acc
c end of change
      COMMON/RESULT/S1,S2,S3,S4
      common/vegas_2/wgt
      common/vegas_3/it_common
      EXTERNAL FXN
      DATA XL,XU/10*0.D0,10*1.D0/
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/1/
c replaces the hard coded cutoff 1.d-70 by joerg zunft
      parameter(tiny = 1.d-36)
c end of change 
      IPR=1
      IF(NPRN.GT.0)IPR=0
      INDO=1
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
      ENTRY VEGAS1(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      NOW=IGRAPH
C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
C---  IF(IGRAPH.GT.0)CALL INPLOT(NOW,F1,W)
      IT=0
c---  copy it to it_common 9/14/02 by tp
      it_common=it
      SI=0.D0
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
      ENTRY VEGAS2(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=(NCALL*0.5D0)**(1.D0/NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2)NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
      IF(ND.EQ.INDO) GO TO 8
      RC=INDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.D0
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6  I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      INDO=ND
      ACC=BCC
ctp remove output for Prospino, including moving tag `8' to next line
ctp8     IF(NPRN.NE.0.AND.NPRN.NE.10)PRINT 200,NDIM,CALLS,IT,ITMX
ctp     1,ACC,MDS,ND
ctp      IF(NPRN.EQ.10)PRINT 290,NDIM,CALLS,ITMX,ACC,MDS,ND
8     continue
      ENTRY VEGAS3(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
9     IT=IT+1
c---  copy it to it_common 9/14/02 by tp
      it_common=it
      TI=0.D0
      TSI=TI
C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
C---  IF(IGRAPH.GT.0)CALL REPLOT(NOW,F1,W)
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      INXI(I,J)=0
      D(I,J)=TI
10    DI(I,J)=TI
11    FB=0.D0
      F2B=FB
      K=0
12    K=K+1
      DO 121 J=1,NDIM
121   QRAN(J)=DBLE(UNIV())
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(KG(J)-QRAN(J))*DXG+ONE
      IA(J)=XN
      IAJ=IA(J)
      IAJ1=IAJ-1
      IF(IAJ.GT.1) GO TO 13
      XO=XI(IAJ,J)
      RC=(XN-IAJ)*XO
      GO TO 14
13    XO=XI(IAJ,J)-XI(IAJ1,J)
      RC=XI(IAJ1,J)+(XN-IAJ)*XO
14    X(J)=XL(J)+RC*DX(J)
15    WGT=WGT*XO*XND
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
C---  IF(IGRAPH.GT.0)CALL XPLOT(NOW,F1,W)
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      IAJ=IA(J)
      INXI(IAJ,J)=INXI(IAJ,J)+1
      DI(IAJ,J)=DI(IAJ,J)+F/CALLS
16    IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
      IF(K.LT.NPG) GO TO 12
      F2B=F2B*NPG
      F2B=SQRT(F2B)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
      IAJ=IA(J)
17    D(IAJ,J)=D(IAJ,J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
      IF (ABS(TSI).LT.tiny) TSI = SIGN(tiny,TSI)
C--------------------------------------------------------
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      SCALLS=SCALLS+CALLS
C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
      IF (ABS(SWGT).LT.tiny) SWGT = SIGN(tiny,SWGT)
      IF (ABS(SI2) .LT.tiny) SI2  = SIGN(tiny,SI2)
C--------------------------------------------------------
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=0.d0
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      SD=ONE/SD
C--------- CHANGE BY J. ZUNFT, 04/14/92-----------------
      IF (SD.LT.0.D0) THEN
      PRINT *,' SD  = ',SD,' < 0: SIGN CHANGED'
      SD = ABS(SD)
      END IF
      IF (TSI.LT.0.D0) THEN
      PRINT *,' TSI = ',TSI,' < 0: SIGN CHANGED'
ctp      TSI = ABS(TSI)
      tsi = tiny
      END IF
C-------------------------------------------------------
      SD=SQRT(SD)
      IF(NPRN.EQ.0) GO TO 21
      TSI=SQRT(TSI)
      IF(NPRN.NE.10)PRINT 201,IPR,IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.EQ.10)PRINT 203,IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
      PRINT 202,J
20    PRINT 204,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
C--------- CHANGE BY J. ZUNFT, 05/06/92------------------
21    IF (ABS(AVGI).LT.tiny) AVGI = SIGN(tiny,AVGI)
C--------------------------------------------------------
      IF(ABS(SD/AVGI).LE.ABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
C---  NEXT LINE TAKEN OUT ON 04/06/92 BY J. ZUNFT
C---  IF(IGRAPH.GT.0)CALL PLOTIT(NOW,F1,W)
C      DO 23 J=1,NDIM
C      XO=D(1,J)
C      XN=D(2,J)
C      D(1,J)=(XO+XN)*0.5D0
C      DT(J)=D(1,J)
C      DO 22 I=2,NDM
C      D(I,J)=XO+XN
C      XO=XN
C      XN=D(I+1,J)
C      D(I,J)=(D(I,J)+XN)/3.D0
C22    DT(J)=DT(J)+D(I,J)
C      D(ND,J)=(XN+XO)*0.5D0
C23    DT(J)=DT(J)+D(ND,J)
C-----THIS PART OF THE VEGAS-ALGORITHM IS UNSTABLE
C-----IT SHOULD BE REPLACED BY
      DO 23 J=1,NDIM
      DT(J)=0.D0
      DO 23 I=1,ND
      IF(INXI(I,J).GT.0)D(I,J)=D(I,J)/INXI(I,J)
23    DT(J)=DT(J)+D(I,J)
      DO 28 J=1,NDIM
      RC=0.D0
      DO 24 I=1,ND
      R(I)=0.D0
C--------- CHANGE BY J. ZUNFT, 04/15/92 ---------
C---  IF(D(I,J).LE.0.D0)GO TO 24
      IF(D(I,J).EQ.0.D0.OR.DT(J)/D(I,J).LE.0.D0) GO TO 24
C------------------------------------------------
      XO=DT(J)/D(I,J)
C--------- CHANGE BY J. ZUNFT, 05/06/92 ---------
ctp      IF (XO.EQ.1.D0) XO = 0.99999999999999999D0
c--------- change by t. plehn, 07/31/01 ---------
      IF (XO.EQ.1.D0) XO = 0.999999999D0
C------------------------------------------------
      R(I)=((XO-ONE)/XO/LOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
 285  XN=0.D0
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
C----------- CHANGE BY J. ZUNFT, 04/14/92 --------------
C---  XIN(I)=XN-(XN-XO)*DR/R(K)
      IF (DR.EQ.0.D0) THEN
      XIN(I)=XN
      ELSE
      XIN(I)=XN-(XN-XO)*DR/R(K)
      END IF
C-------------------------------------------------------
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.ABS(ACC).LT.ABS(SD/AVGI))GO TO 9
200   FORMAT(35H0INPUT PARAMETERS FOR VEGAS   NDIM=,I3
     1,8H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X
     2,6H  ACC=,G9.3/28X,6H  MDS=,I3,6H   ND=,I4//)
290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5
     1,6H  ACC=,G9.3,6H  MDS=,I3,6H   ND=,I4)
201   FORMAT(/I1,20HINTEGRATION BY VEGAS/13H0ITERATION NO,I3,
     114H.   INTEGRAL =,G14.8/20X,10HSTD DEV  =,G10.4/
     234H ACCUMULATED RESULTS.   INTEGRAL =,G14.8/
     324X,10HSTD DEV  =,G10.4 / 24X,18HCHI**2 PER ITN   =,G10.4)
202   FORMAT(14H0DATA FOR AXIS,I2 / 7X,1HX,7X,10H  DELT I  ,
     12X,11H CONVCE    ,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE
     2,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE     /)
204   FORMAT(1X,3G12.4,5X,3G12.4,5X,3G12.4)
203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
C---  NEXT 3 LINES TAKEN OUT ON 04/06/92 BY J. ZUNFT
C---  S1=AVGI
C---  S2=SD
C---  S3=CHI2A
      RETURN
      END





C-----------------------------------------------------------------------
C   INITIALIZING THE RANDOM NUMBER GENERATOR
C     IF OLD CONFIGURATION EXISTS THEN
C     CALL UREAD(11)
C     ELSE
C     CALL RSTART(12,34,56,78)
C     END IF
C-----------------------------------------------------------------------


C----------------------------------------------------------------------
C  A UNIVERSAL RANDOM NUMBER GENERATOR

        FUNCTION UNIV()
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,I,J
        UNIV=U(I)-U(J)
        IF(UNIV.LT.0.) UNIV=UNIV+1.
        U(I)=UNIV
        I=I-1
        IF(I.EQ.0) I=97
        J=J-1
        IF(J.EQ.0) J=97
        C=C-CD
        IF(C.LT.0.) C=C+CM
        UNIV=UNIV-C
        IF(UNIV.LT.0.) UNIV=UNIV+1
        RETURN
        END

C----------------------------------------------------------------------
C SAVING THE RANDOM NUMBER GENERATOR CONFIGURATION

        SUBROUTINE UWRITE(IUNIT)
        INTEGER IUNIT
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,I,J
ctp changed by tp: status
        OPEN(IUNIT,FILE='TESTRN',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IUNIT) U,C,CD,CM,I,J
        CLOSE(IUNIT)
        RETURN
        END


C----------------------------------------------------------------------
C READING IN THE RANDOM NUMBER GENERATOR CONFIGURATION

        SUBROUTINE UREAD(IUNIT)
        INTEGER IUNIT
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,I,J
ctp changed by tp: status
        OPEN(IUNIT,FILE='TESTRN',FORM='UNFORMATTED',STATUS='UNKNOWN')
        READ(IUNIT) U,C,CD,CM,I,J
        CLOSE(IUNIT)
        RETURN
        END

C----------------------------------------------------------------------
C INITIALIZING THE RANDOM NUMBER GENERATOR
C TO INITIALIZE CALL RSTART(12,34,56,78)

        SUBROUTINE RSTART(I,J,K,L)
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,ISTART,JSTART
        IF ((I.LT.0).OR.(I.GT.178 )) STOP 'FIRST SEED .LT.0 OR .GT.178'
        IF ((J.LT.0).OR.(J.GT.178 )) STOP 'SECOND SEED .LT.0 OR .GT.178'
        IF ((K.LT.0).OR.(K.GT.178 )) STOP 'THIRD SEED .LT.0 OR .GT.178'
        IF ((L.LT.0).OR.(L.GT.168 )) STOP 'FOURTH SEED .LT.0 OR .GT.168'
        IF ( (I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1) ) STOP
     &     'FIRST, SECOND AND THIRD SEEDS ARE ALL EQUAL TO 1'
        ISTART=97
        JSTART=33
        IDUM = I
        JDUM = J
        KDUM = K
        LDUM = L
        DO 2 II=1,97
        S=0.
        T=.5
        DO 3 JJ=1,24
          M=MOD(MOD(IDUM*JDUM,179)*KDUM,179)
          IDUM=JDUM
          JDUM=KDUM
          KDUM=M
          LDUM=MOD(53*LDUM+1,169)
          IF(MOD(LDUM*M,64).GE.32) S=S+T
3         T=.5*T
2         U(II)=S
        C=362436./16777216.
        CD=7654321./16777216.
        CM=16777213./16777216.
        RETURN
        END

