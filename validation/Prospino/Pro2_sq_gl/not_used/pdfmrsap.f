C***************************************************************C
C								C
C     THIS IS A PACKAGE FOR THE NEW MRS(A') PARTON              C
C     DISTRIBUTION. THE MINIMUM Q^2  VALUE IS 5 GEV^2 AND THE   C
C     X RANGE IS, AS BEFORE 10^-5 < X < 1. MSBAR FACTORIZATION  C
C     IS USED. THE PACKAGE READS 1 GRID, WHICH IS IN A SEPARATE C
C     FILE (A'=mrsap.grid).                                     C  
C     NOTE THAT X TIMES THE PARTON DISTRIBUTION IS RETURNED,    C
C     Q IS THE SCALE IN GEV,                                    C
C     AND LAMBDA(MSBAR,NF=4) = 231 MEV FOR A'.                  C
C								C
C         THE REFERENCE IS :                                    C
C         A.D. MARTIN, R.G. ROBERTS AND W.J. STIRLING,          C
C         RAL PREPRINT  RAL-95-021 (1995)                       C
C                                                               C
C         COMMENTS TO : W.J.STIRLING@DURHAM.AC.UK               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE GRID IS 0.00341               C
C								C
C***************************************************************C

      SUBROUTINE MRSAP(X,SCALE,YPDF)
C***  ADDED BY RH TO HAVE THE SAME CONVENTIONS AS IN PDG
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I,J)
      DIMENSION YPDF(-6:6)
      Q2=SCALE**2
      IF(Q2.LT.5D0.OR.Q2.GT.1310720.D0)    PRINT 99
      IF(X.LT.1D-5.OR.X.GT.1D0)            X = 1.D0           
      CALL STRC30(X,SCALE,UP1,DOWN1,UP2,DOWN2,STRANGE,CHARM,BOTTOM,
     +     GLUON)
      YPDF(0) = GLUON
      YPDF(1) = DOWN1 + DOWN2
      YPDF(2) = UP1 + UP2
      YPDF(3) = STRANGE
      YPDF(4) = CHARM
      YPDF(5) = BOTTOM
      YPDF(6) = 0.D0
      YPDF(-1) = DOWN2
      YPDF(-2) = UP2
      YPDF(-3) = STRANGE
      YPDF(-4) = CHARM
      YPDF(-5) = BOTTOM
      YPDF(-6) = 0.D0
  99  FORMAT('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  FORMAT('  WARNING:   X  VALUE IS OUT OF RANGE   ')
      RETURN
      END



C
      SUBROUTINE STRC30(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)

C     THIS IS THE NEW  "APRIME" FIT -- FEB 1995 -- STANDARD Q^2 RANGE

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=47)
      PARAMETER(NTENTH=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
C***  ADDED BY RH (FOR SOME FORTRAN COMPILERS): 
      COMMON / GRID / F, INIT
C***  END OF ADDITION.
      DATA XX/1.D-5,2.D-5,4.D-5,6.D-5,8.D-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/
 
 
      XSAVE=X
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
C***  ADDED BY RH 
      OPEN(33,FILE='mrsap.grid')
C***  END OF ADDITION.
      DO 20 N=1,NX-1
      DO 20 M=1,19
      READ(33,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,NTENTH,K)/DLOG(F(I,NTENTH,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,NX,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.NTENTH) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,NTENTH,M)+B*F(I,NTENTH,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      X=XSAVE
 
      RETURN
      END
C

  
