C============================================================================
C                CTEQ Parton Distribution Functions: version 6.0-6.6
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C                             August 6, 2003, v6.11
C                             December 12, 2004, v6.12
C                             December 4, 2006, v6.5 (CTEQ6.5M series added)
C                             March 23, 2007, v6.51 (CTEQ6.5S/C series added)
C                             April 24, 2007, v6.52 (minor improvement)
C                             March 30, 2008, v6.6 
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
C       JHEP 0310:046(2003), hep-ph/0303013
C
C   Ref[3]: "Neutrino dimuon Production and Strangeness Asymmetry of the Nucleon"
C       By: F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky, H.L. Lai, S. Kretzer, J.F. Owens, W.K. Tung
C       Eur. Phys. J. C40:145(2005), hep-ph/0312323
C
C   Ref[4]: "CTEQ6 Parton Distributions with Heavy Quark Mass Effects"
C       By: S. Kretzer, H.L. Lai, F. Olness, W.K. Tung
C       Phys. Rev. D69:114005(2004), hep-ph/0307022
C
C   Ref[5]: "Heavy Quark Mass Effects in Deep Inelastic Scattering and Global QCD Analysis"
C       By : W.K. Tung, H.L. Lai, A. Belyaev, J. Pumplin, D. Stump, C.-P. Yuan
C       JHEP 0702:053(2007), hep-ph/0611254
C
C   Ref[6]: "The Strange Parton Distribution of Nucleon: Global Analysis and Applications"
C       By : H.L. Lai, P. Nadolsky, J. Pumplin, D. Stump, W.K. Tung, C.-P. Yuan
C       JHEP 0704:089,2007, hep-ph/0702268
C
C   Ref[7]: "The Charm Content of the Nucleon"
C       By : J. Pumplin, H.L. Lai, W.K. Tung
C       Phys.Rev.D75:054029,2007, hep-ph/0701220

C   Ref[8]: "Implications of CTEQ global analysis for collider observables"
C       By : P. M. Nadolsky, H.-L. Lai, Q.-H. Cao, J. Huston, J. Pumplin, D. R. Stump, W.-K. Tung, C.-P. Yuan
C       arXiv:0802.0007 [hep-ph], submitted to Phys. Rev. D. 
C

C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
C   (4) 5 special sets for strangeness study from Ref[3].
C   (5) 1 special set for heavy quark study from Ref[4].
C   (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref[5].
C   (7) 8 sets of PDFs resulting from the strangeness study, Ref[6].
C   (8) 7 sets of PDFs resulting from the charm study, Ref[7].
C   (9) CTEQ6.6M and its 44 up/down eigenvector sets from Ref[8].
C  (10) Fits with nonperturbative charm from the study in  Ref[8].
C  (11) Fits with alternative values of the strong coupling strength from the study in Ref[8].


C  Details about the calling convention are:
C --------------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File   Ref
C ================================================================================
C Standard, "best-fit", sets:                 
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl    [1]
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl    [1]
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl    [1]
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl   [1]
C 200    CTEQ6.1M: updated CTEQ6M (see below, under "uncertainty" section)     [2]
C 400    CTEQ6.6M; the 2008 set (see below, under "uncertainty" section)       [8]
C
C --------------------------
C  Special sets with nonperturbative charm at Q_0=1.3 GeV from Ref [8]
C --------------------------
C 450    CTEQ6.6C1   BHPS model for IC    0.118     326   226    ctq66.c1.pds
C 451    CTEQ6.6C2   BHPS model for IC    0.118     326   226    ctq66.c2.pds
C 452    CTEQ6.6C3   Sea-like model       0.118     326   226    ctq66.c3.pds
C 453    CTEQ6.6C4   Sea-like model       0.118     326   226    ctq66.c4.pds
C     Momentum Fraction carried by c+cbar=2c at Q0=1.3 GeV:
C    Iset:     451  452   453   454 
C Mom. frac:  0.01 0.035  0.01  0.035


C --------------------------
C  Special CTEQ6.6 sets with alternative values of strong coupling strength [8]
C --------------------------
C 460    CTEQ6.6A1                        0.125           328    ctq66.a1.pds
C 461    CTEQ6.6A2                        0.122           281    ctq66.a2.pds
C 462    CTEQ6.6A3                        0.114           179    ctq66.a3.pds
C 463    CTEQ6.6A4                        0.112           159    ctq66.a4.pds

C --------------------------
C Special sets for strangeness study:  Ref.[3]
C --------------------------
C  11    CTEQ6A   Class A                 0.118     326   226    cteq6sa.pds
C  12    CTEQ6B   Class B                 0.118     326   226    cteq6sb.pds
C  13    CTEQ6C   Class C                 0.118     326   226    cteq6sc.pds
C  14    CTEQ6B+  Large [S-]              0.118     326   226    cteq6sb+.pds
C  15    CTEQ6B-  Negative [S-]           0.118     326   226    cteq6sb-.pds
C --------------------------
C Special set for Heavy Quark study:   Ref.[4]
C --------------------------
C  21    CTEQ6HQ                          0.118     326   226    cteq6hq.pds
C --------------------------
C Released sets for strangeness study:  Ref.[6]
C -------------------------- s=sbr
C  30    CTEQ6.5S0   Best-fit             0.118     326   226    ctq65.s+0.pds
C  31    CTEQ6.5S1   Low s+               0.118     326   226    ctq65.s+1.pds
C  32    CTEQ6.5S2   High s+              0.118     326   226    ctq65.s+2.pds
C  33    CTEQ6.5S3   Alt Low s+           0.118     326   226    ctq65.s+3.pds
C  34    CTEQ6.5S4   Alt High s+          0.118     326   226    ctq65.s+4.pds
C -------------------------- s!=sbr
C          strangeness asymmetry <x>_s-
C  35    CTEQ6.5S-0  Best-fit    0.0014    0.118     326   226    ctq65.s-0.pds
C  36    CTEQ6.5S-1  Low        -0.0010    0.118     326   226    ctq65.s-1.pds
C  37    CTEQ6.5S-2  High        0.0050    0.118     326   226    ctq65.s-2.pds
C --------------------------
C Released sets for charm study:  Ref.[7]
C --------------------------
C  40    CTEQ6.5C0   no intrinsic charm   0.118     326   226    ctq65.c0.pds
C  41    CTEQ6.5C1   BHPS model for IC    0.118     326   226    ctq65.c1.pds
C  42    CTEQ6.5C2   BHPS model for IC    0.118     326   226    ctq65.c2.pds
C  43    CTEQ6.5C3   Meson cloud model    0.118     326   226    ctq65.c3.pds
C  44    CTEQ6.5C4   Meson cloud model    0.118     326   226    ctq65.c4.pds
C  45    CTEQ6.5C5   Sea-like model       0.118     326   226    ctq65.c5.pds
C  46    CTEQ6.5C6   Sea-like model       0.118     326   226    ctq65.c6.pds
C
C     Momentum Fraction carried by c,cbar at Q0=1.3 GeV:
C    Iset:charm  ,cbar     | Iset:charm  ,cbar     | Iset:charm  ,cbar
C    41: 0.002857,0.002857 | 43: 0.003755,0.004817 | 45: 0.005714,0.005714
C    42: 0.010000,0.010000 | 44: 0.007259,0.009312 | 46: 0.012285,0.012285
C
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects, Ref[5]:  central fit: CTEQ6.5M (=CTEQ65.00)
C                              -----------------------
C  3xx  CTEQ65.xx  +/- sets               0.118     326   226    ctq65.xx.pds
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 300      is CTEQ65.00 (=CTEQ6.5M),
C             301/302 are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects and free strangeness, Ref[8]:  
C                central fit: CTEQ6.6M (=CTEQ66.00)
C                              -----------------------
C  4xx  CTEQ66.xx  +/- sets               0.118     326   226    ctq66.xx.pds
C        where xx = 01-44: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 400      is CTEQ66.00 (=CTEQ6.6M),
C             401/402 are CTEQ66.01/02, +/- sets of 1st eigenvector, ... etc.

C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for 
C    *  10^-8 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.6 series;
C    *  10^-7 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.5S/C series;
C    *  10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV) for CTEQ6, CTEQ6.1 series;
C
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   nadolsky@pa.msu.edu, pumplin@pa.msu.edu or tung@pa.msu.edu.
C
C===========================================================================

      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0d0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Ctq6Pdf = 0D0
        Return
      Endif

      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif

      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if (Ctq6Pdf.lt.0.D0) Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=8)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.', 'ctq66.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
      Data Isetmin4,Isetmax4 /400,444/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetmnC08,IsetmxC08 /450,453/
      Data IsetmnAS08,IsetmxAS08 /460,463/

      Data IsetHQ /21/
      Common /Setchange/ Isetch
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
        fmtpds=.true.

        If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                  Iset = 1,2,3 for 6m, 6d, 6l
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'.tbl'
        Elseif (Iset.eq.4) Then
C                                                             4  (2nd LO fit)
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'1.tbl'
        Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(1)//nn//'.tbl'
        Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(5)//nn(2:3)//'.tbl'
        Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15
          If(Iset.eq.11) then
            Tablefile=Flnm(6)//'a.pds'
          Elseif(Iset.eq.12) then
            Tablefile=Flnm(6)//'b.pds'
          Elseif(Iset.eq.13) then
            Tablefile=Flnm(6)//'c.pds'
          Elseif(Iset.eq.14) then
            Tablefile=Flnm(6)//'b+.pds'
          Elseif(Iset.eq.15) then
            Tablefile=Flnm(6)//'b-.pds'
          Endif
        Elseif (Iset.eq.IsetHQ) Then
C                                                               21
          TableFile='cteq6hq.pds'
        Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'s+'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37
          Is = Iset - 5
          write(nn,'(I2)') Is
          Tablefile=Flnm(7)//'s-'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'c'//nn(2:2)//'.pds'
        Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340
          write(nn,'(I3)') Iset
          Tablefile=Flnm(7)//nn(2:3)//'.pds'
        Elseif (Iset.ge.Isetmin4 .and. Iset.le.Isetmax4) Then
C                                                    (Cteq6.6)  400 - 444   
          write(nn,'(I3)') Iset
          Tablefile=Flnm(8)//nn(2:3)//'.pds'
        Elseif (Iset.ge.IsetmnC08 .and. Iset.le.IsetmxC08) Then
C                                                   (Cteq6.6C)  450 - 453
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.ge.IsetmnAS08 .and. Iset.le.IsetmxAS08) Then
C                                                   (Cteq6.6AS)  460 - 463
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'a'//nn(3:3)//'.pds'
        Else
          Print *, 'Invalid Iset number in SetCtq6 :', Iset
          Stop
        Endif
        IU= NextUn()
        Open(IU, File=Tablefile, Status='OLD', Err=100)
 21     Call Readpds (IU,fmtpds)
        Close (IU)
        Isetold=Iset
        Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >  //'in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
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
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END

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
