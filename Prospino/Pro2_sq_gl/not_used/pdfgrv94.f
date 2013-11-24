
********************************************************************
*                                                                  *
*   GRV PARTON DENSITIES, 1994 UPDATE  (GRID VERSION)              *
*                                                                  *
*   INPUT:   X  = Bjorken-x       (between  1.E-5  and  1)         *
*            Q2 = scale in GeV**2 (between  1  and   1.E6)         *
*             (for values outside this allowed range the program   * 
*              writes a warning and extrapolates to the x and      * 
*              Q2 values requested)                                * 
*                                                                  *
*            ISET = number of the parton set :                     *
*            ISET = 0    :   LO                                    *
*            ISET = ELSE :   NLO(MSbar)                            *
*                                                                  *
*   OUTPUT:  UV = u - u(bar), DV = d - d(bar), US = u(bar),        *
*            DS = d(bar), SS = s = s(bar), GL = gluon.             *
*            Always x times the distribution is returned in the    *
*            (modified) MS(bar) scheme (only 3 parton flavours,    *
*            see the reference below for a detailed discussion).   *
*                                                                  *
*   CHARM:   For applications where a massive treatment is not     *
*            possible since, e.g., the corresponding theoretical   *
*            expressions are not available, this programs returns  *
*            effective charm and bottom densities, CS = c =c(bar)  *
*            and BS = b = b(bar), as obtained from the previous    *
*            (1991) GRV parametrization.                           *
*                                                                  *   
*   REFERENCE:  M. Glueck, E. Reya and A. Vogt, DESY 94-206        *
*                 [published in Z. Phys. C67 (1995) 433]           *
*                                                                  *
*   COMMON:  The main program or the calling routine has to have   *
*            a common block  COMMON / INTINI / IINI , and the in-  *
*            teger variable  IINI  has always to be zero when      *
*            GRVPAR is called for the first time or when 'ISET'    *
*            has been changed.                                     *
*                                                                  *
*   COMMENTS, QUESTIONS ETC TO:  avogt@x4u2.desy.de                *
*                                                                  *
********************************************************************
*

      SUBROUTINE GRV94(X,SCALE,YPDF,ISET)
C***  ADDED BY RH TO HAVE THE SAME CONVENTIONS AS IN PDG
      IMPLICIT INTEGER (I,J)
      REAL*8 x, scale, YPDF(-6:6)
      COMMON / INTINI / IINI

      Q2=real(SCALE)**2
      IF(Q2.LT.5e0.OR.Q2.GT.2.e6)     PRINT 99
      IF(X.LT.1d-5.OR.X.GT.1.d0)            X = 1.D0           

      y = real(x)

      CALL GRVPAR(ISET,y,Q2,UP1,DOWN1,UP2,DOWN2,STRANGE,GLUON,CHARM,
     +     BOTTOM)
         
      YPDF(0) = dble( GLUON )
      YPDF(1) = dble( DOWN1 + DOWN2 )
      YPDF(2) = dble( UP1 + UP2 )
      YPDF(3) = dble( STRANGE )
      YPDF(4) = dble( CHARM )
      YPDF(5) = dble( BOTTOM )
      YPDF(6) = 0.D0
      YPDF(-1) = dble( DOWN2 )
      YPDF(-2) = dble( UP2 )
      YPDF(-3) = dble( STRANGE )
      YPDF(-4) = dble( CHARM )
      YPDF(-5) = dble( BOTTOM )
      YPDF(-6) = 0.D0

  99  FORMAT('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  FORMAT('  WARNING:   X  VALUE IS OUT OF RANGE   ')
      RETURN
      END


      SUBROUTINE GRVPAR (ISET, X, Q2, UV, DV, US, DS, SS, GL, CS, BS)
      PARAMETER (NPART=6, NX=49, NQ=25, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XUSF(NX,NQ), XDSF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          XCF(NX,NQ), XBF(NX,NQ), HEAVY (2,NQ,NX-1),
     3          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      COMMON / INTINI / IINI
      SAVE XUVF, XDVF, XUSF, XDSF, XSF, XGF, XCF, XBF, NA, ARRF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.0E0, 1.5E0, 2.5E0, 4.0E0, 6.4E0, 
     1           1.0E1, 1.5E1, 2.5E1, 4.0E1, 6.4E1, 
     2           1.0E2, 1.8E2, 3.2E2, 5.8E2, 
     3           1.0E3, 1.8E3, 3.2E3, 5.8E3, 
     4           1.0E4, 2.2E4, 4.6E4, 1.0E5, 2.2E5, 4.6E5, 1.E6 / 
       DATA XB / 1.E-5, 1.5E-5, 2.2E-5, 3.2E-5, 4.8E-5, 7.E-5,
     1           1.E-4, 1.5E-4, 2.2E-4, 3.2E-4, 4.8E-4, 7.E-4,
     2           1.E-3, 1.5E-3, 2.2E-3, 3.2E-3, 4.8E-3, 7.E-3,
     3           1.E-2, 1.5E-2, 2.2E-2, 3.2E-2, 5.0E-2, 7.5E-2,
     4           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,  
     6           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.1.0E-5) .OR. (X.GT.1.0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
       ENDIF
C***  UPPER VALUE FOR Q2 CHANGED BY RH
       IF ( (Q2.LT.1.0) .OR. (Q2.GT.2.E6) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    (COMMENT: FIRST NON-VANISHING NUMBER IN THE GRID)
      IF (IINI.NE.0) GOTO 16
      IF (ISET.EQ.0) THEN
        OPEN(11,FILE='grvl_94.grid')  !  3.128E-03
        OPEN(12,FILE='hql_91.grid')   !  7.936E-02 
      ELSE
        OPEN(11,FILE='grv_94.grid')   !  3.662E-03
        OPEN(12,FILE='hq_91.grid')    !  7.225E-02  
      END IF
      IINI = 1
      READ(11,89) LINE
      READ(12,89) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(11,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1            PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M)
      READ(12,94) HEAVY(1,N,M), HEAVY(2,N,M) 
  90  FORMAT (6(1PE10.3))
  94  FORMAT (2(1PE10.3))
  15  CONTINUE
      CLOSE(11)
      CLOSE(12)
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0**0.5)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.3)
        XUSF(IX,IQ) = PARTON(3,IQ,IX) /  XB1**7 * XB0**0.2 
        XDSF(IX,IQ) = PARTON(4,IQ,IX) /  XB1**7 * XB0**0.2 
        XSF(IX,IQ)  = PARTON(5,IQ,IX) /  XB1**7 * XB0**0.2 
        XGF(IX,IQ)  = PARTON(6,IQ,IX) /  XB1**5 * XB0**0.2
        XCF(IX,IQ)  = HEAVY(1,IQ,IX) /  XB1**7 * XB0**0.2 
        XBF(IX,IQ)  = HEAVY(2,IQ,IX) /  XB1**7 * XB0**0.2 
  20  CONTINUE
        XUVF(NX,IQ) = 0.E0
        XDVF(NX,IQ) = 0.E0
        XUSF(NX,IQ) = 0.E0
        XDSF(NX,IQ) = 0.E0
        XSF(NX,IQ)  = 0.E0
        XGF(NX,IQ)  = 0.E0
        XCF(NX,IQ)  = 0.E0
        XBF(NX,IQ)  = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = ALOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = ALOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = ALOG(X)
      XT(2) = ALOG(Q2)
      UV = FFINT(NARG,XT,NA,ARRF,XUVF) * (1.-X)**3 * X**0.5
      DV = FFINT(NARG,XT,NA,ARRF,XDVF) * (1.-X)**4 * X**0.3 
      US = FFINT(NARG,XT,NA,ARRF,XUSF) * (1.-X)**7 / X**0.2
      DS = FFINT(NARG,XT,NA,ARRF,XDSF) * (1.-X)**7 / X**0.2
      SS = FFINT(NARG,XT,NA,ARRF,XSF)  * (1.-X)**7 / X**0.2
      GL = FFINT(NARG,XT,NA,ARRF,XGF)  * (1.-X)**5 / X**0.2
      CS = FFINT(NARG,XT,NA,ARRF,XCF)  * (1.-X)**7 / X**0.2
      BS = FFINT(NARG,XT,NA,ARRF,XBF)  * (1.-X)**7 / X**0.2
 60   RETURN
      END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FFINT(NARG,ARG,NENT,ENT,TABLE)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FFINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FFINT=FFINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END




