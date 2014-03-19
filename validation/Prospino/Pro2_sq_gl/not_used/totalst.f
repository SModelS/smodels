C**********************************************************************
C***                                                                ***
C***  THIS PROGRAM CALCULATES                                       ***
C***  THE TOTAL CROSS-SECTION AND DISTRIBUTIONS                     ***
C***  FOR STOP-ANTISTOP PRODUCTION AT HADRON COLLIDERS              ***
C***  INCLUDING FULL SUSY-QCD CORRECTIONS                           ***
C***                                                                ***
C***  WRITTEN BY: W. BEENAKKER, R. HOPKER AND T. PLEHN              ***
C***                                                                ***
C**********************************************************************

      PROGRAM TOTALST
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      COMMON/IOUT/IPRINT
      COMMON/CONST1/S,ENERGY,ALPHAS,MST1,MG,MT,MS,MST2,SIN2T
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT
      COMMON/CONST3/IPDFSET
      COMMON/CONST5/ILO,INLO,IONLYLO
      COMMON/CUT1/PTMIN,PTMAX
      COMMON/CUT2/YMIN,YMAX
      COMMON/FLAVOR/IFLAVOR,ITOTAL
      COMMON/CHARCONJ/ICHARCONJ

      COMMON/W50512/QCDL4,QCDL5

C ---------------------------------------------------------------
C --- INPUT PARAMETERS, CAN BE CHANGED  -------------------------
C ---------------------------------------------------------------

C***  FOR ST2 PAIR PRODUCTION REVERSE THE SIGN OF SIN2T AND INTERCHANGE
C***  THE MASSES MST1 AND MST2.

C***  THE MASSES (IN GEV)

      MST1 = 153.D0
      MST2 = 347.D0
      MG = 284.D0
      MS = 256.D0 
      SIN2T = -0.99D0

      MT = 175.D0

C***  THE COLLIDER TYPE ( P PBAR = 0, P P = 1 )
      
      ICOLL = 1

C***  THE CENTER OF MASS ENERGY (IN GEV)

      ENERGY = 14000D0

C***  THE SET OF PARTON DENSITIES (GRV = 0, MRSAP = 1, PDFLIB = 2)
C***  FOR PDFLIB PLEASE MAKE CHANGES IN SUBROUTINES INILO AND ININLO
C***  IN FILE INITPDF.F

      IPDFSET = 2

C***  THE INITIAL STATE 
C***  ALL = 0, G G = 1, Q Q = 2, G Q = 3, MAJOR=4, MINOR=5

      IFLAVOR = 0

C***  THE CROSS-SECTION WITH CUTS (0)
C***  THE TOTAL CROSS-SECTION WITHOUT CUTS IN A FASTER WAY (1)

      ITOTAL = 1

C***  THE SCALE FOR RENORMALIZATION AND FACTORIZATION
C***  ISCAPT = 0  --> SCALE = MST1 * SCAFAC 
C***  ISCAPT = 1  --> SCALE = SQRT(MST1**2 +PT**2) * SCAFAC
C***                  ONLY FOR ITOTAL = 0 
C***  DEFAULT FOR SCAFAC = 1.0

      ISCAPT = 0
      SCAFAC = 1.D0

C***  THE CUTS ON THE CROSS-SECTION IN PT (DEFAULT: 0,ENERGY )
C***  PT IS ONLY DEFINED FOR POSITIVE VALUES: PTMIN >= 0
C***  EQUAL LOWER AND UPPER CUT GIVES DSIGMA/DPT
C***  ITOTAL = 0 NECESSARY

      PTMIN = 0D0
      PTMAX = ENERGY

C***  THE CUTS ON THE CROSS-SECTION IN Y (DEFAULT: 0,+9.99)
C***  RAPIDITY IS ONLY DEFINED FOR POSITIVE VALUES: YMIN >=0
C***  EQUAL LOWER AND UPPER CUT GIVES DSIGMA/DY
C***  ITOTAL = 0 NECESSARY
  
      YMIN =  0.D0
      YMAX = +9.99D0

C***  ONLY FOR DISTRIBUTIONS ( ITOTAL = 0 )
C***  DISTINGUISH BETWEEN STOPS AND ANTISTOPS IN THE FINAL STATE
C***  DEFINES THE DIFFERENTIAL CROSS-SECTIONS WITH RESPECT TO
C***  ICHARCONJ =  1  <-- STOPS
C***  ICHARCONJ = -1  <-- ANTISTOPS
C***  ICHARCONJ =  0  <-- AVERAGE OF STOPS AND ANTISTOPS

      ICHARCONJ = 0

C***  IONLYLO = 0 CALCULATES BORN AND NLO CROSS-SECTIONS
C***  IONLYLO = 1 CALCULATES ONLY THE BORN CROSS-SECTION
      
      IONLYLO = 0

C***  THE NUMBER OF VEGAS CALLS ( DEFAULT = 1000 )

      ILO   = 1000
      INLO  = 500

C***  PRINT VEGAS STATISTICS (10) OR NOT (0)

      IPRINT = 0

C --- PRINT THE HEADER ----------------------------------------
      
      CALL PRIHEADST

C --- INITIALIZE VEGAS -----------------------------------------
     
      CALL RSTART(12,34,56,78)

C ---------------------------------------------------------------
C --- INTEGRATION BY VEGAS --------------------------------------
C --- CALCULATION OF CROSS-SECTIONS AND DISTRIBUTIONS -----------
C --- CAN BE CHANGED --------------------------------------------
C ---------------------------------------------------------------

      sin2t = 0.50D0
      do ims=160,1000,20
             
         mst1 = dble(ims)
         mst2 = mst1 + 200.D0
         ms   = mst1 + 100.D0
         mg   = mst1 + 250.D0
         
         CALL INTEGST(RESLO,ERRLO,RESNLO,ERRNLO)
         CALL PRIRESST(RESLO,ERRLO,RESNLO,ERRNLO)
         
      end do

      STOP
      END

C**********************************************************************
