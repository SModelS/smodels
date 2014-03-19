C============================================================================
C                CTEQ Parton Distribution Functions: Version 4.6
C                             June 21, 1996
C                   Modified: 10/17/96, 1/7/97, 1/15/97
C                             2/17/97, 2/21/97
C                   Last Modified on April 2, 1997
C
C   Ref[1]: "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF RECENT DEEP
C         INELASTIC SCATTERING AND INCLUSIVE JET DATA"
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C       Phys. Rev. D55, 1280 (1997)
C
C   Ref[2]: "CHARM PRODUCTION AND PARTON DISTRIBUTIONS"
C   By: H.L. Lai and W.K. Tung
C       MSU-HEP-61222, CTEQ-622, e-Print Archive: hep-ph/9701256
C       to appear in Z. Phys.
C
C   This package contains 13 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C Ref[1]
C   1    CTEQ4M   Standard MSbar scheme   0.116     298   202    cteq4m.tbl
C   2    CTEQ4D   Standard DIS scheme     0.116     298   202    cteq4d.tbl
C   3    CTEQ4L   Leading Order           0.132     236   181    cteq4l.tbl
C   4    CTEQ4A1  Alpha_s series          0.110     215   140    cteq4a1.tbl
C   5    CTEQ4A2  Alpha_s series          0.113     254   169    cteq4a2.tbl
C   6    CTEQ4A3            ( same as CTEQ4M )
C   7    CTEQ4A4  Alpha_s series          0.119     346   239    cteq4a4.tbl
C   8    CTEQ4A5  Alpha_s series          0.122     401   282    cteq4a5.tbl
C   9    CTEQ4HJ  High Jet                0.116     303   206    cteq4hj.tbl
C   10   CTEQ4LQ  Low Q0                  0.114     261   174    cteq4lq.tbl
C ---------------------------------------------------------------------------
C Ref[2]
C   11   CTEQ4HQ  Heavy Quark             0.116     298   202    cteq4hq.tbl
C   12   CTEQ4HQ1 Heavy Quark:Q0=1,Mc=1.3 0.116     298   202    cteq4hq1.tbl
C        (Improved version of CTEQ4HQ, recommended)
C   13   CTEQ4F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=385)   cteq4f3.tbl
C   14   CTEQ4F4  Nf=4 FixedFlavorNumber  0.111     292   XXX    cteq4f4.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ(4HQ1) for which Q starts at a lower value of 0.7(1.0) GeV.  
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=5.0 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq4(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq4Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ4F3 has, by definition, only 3 flavors and gluon;
C              CTEQ4F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq4Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq4Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq4Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq4Pdf: '
     >              , Iparton
         Endif
         Ctq4Pdf = 0D0
         Return
      Endif

      Ctq4Pdf = PartonX (Iparton, X, Q)
      if(Ctq4Pdf.lt.0.D0)  Ctq4Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X < Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q < Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq4 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=14)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq4m.tbl', 'cteq4d.tbl', 'cteq4l.tbl'
     > , 'cteq4a1.tbl', 'cteq4a2.tbl', 'cteq4m.tbl', 'cteq4a4.tbl'
     > , 'cteq4a5.tbl', 'cteq4hj.tbl', 'cteq4lq.tbl'
     > , 'cteq4hq.tbl', 'cteq4hq1.tbl', 'cteq4f3.tbl', 'cteq4f4.tbl' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq4 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq4!!'
      Stop
 101  Print*, Tablefile, ' cannot be opened '
      Print*, 'Please input the .tbl file:'
      Read (*,'(A)') Tablefile
      Goto 21
C                             ********************
      End

      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
