C**********************************************************************
C***                                                                ***
C***  THIS PROGRAM CALCULATES                                       ***
C***  THE TOTAL CROSS-SECTION AND DISTRIBUTIONS                     ***
C***  FOR GLUINO-GLUINO PRODUCTION AT HADRON COLLIDERS              ***
C***  INCLUDING FULL SUSY-QCD CORRECTIONS                           ***
C***                                                                ***
C***  WRITTEN BY: W. BEENAKKER, R. HOPKER AND M. SPIRA              ***
C***                                                                ***
C**********************************************************************

      PROGRAM TOTALGG
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      COMMON/IOUT/IPRINT
      COMMON/CONST1/S,ENERGY,ALPHAS,MS,MG,MT
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT
      COMMON/CONST3/IPDFSET
      COMMON/CONST5/ILO,INLO,IONLYLO
      COMMON/CUT1/PTMIN,PTMAX
      COMMON/CUT2/YMIN,YMAX
      COMMON/FLAVOR/IFLAVOR,ITOTAL

C ---------------------------------------------------------------
C --- INPUT PARAMETERS, CAN BE CHANGED  -------------------------
C ---------------------------------------------------------------

C***  THE MASSES (IN GEV)

      MS = 280.D0      
      MG = 200.D0
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
C***  ISCAPT = 0  --> SCALE = MG * SCAFAC 
C***  ISCAPT = 1  --> SCALE = SQRT(MG**2 +PT**2) * SCAFAC
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
C***  RAPIDITY IS ONLY DEFINED FOR POSITIVE VALUES: YMIN >= 0
C***  EQUAL LOWER AND UPPER CUT GIVES DSIGMA/DY
C***  ITOTAL = 0 NECESSARY
  
      YMIN =  0.D0
      YMAX = +9.99D0

C***  IONLYLO = 0 CALCULATES BORN AND NLO CROSS-SECTIONS
C***  IONLYLO = 1 CALCULATES ONLY THE BORN CROSS-SECTION
      
      IONLYLO = 0

C***  THE NUMBER OF VEGAS CALLS ( DEFAULT = 1000 )

      ILO   = 1000
      INLO  = 500

C***  PRINT VEGAS STATISTICS (10) OR NOT (0)

      IPRINT = 0

C --- PRINT THE HEADER ----------------------------------------
      
      CALL PRIHEADGG

C --- INITIALIZE VEGAS -----------------------------------------

      CALL RSTART(12,34,56,78)

C ---------------------------------------------------------------
C --- INTEGRATION BY VEGAS --------------------------------------
C --- CALCULATION OF CROSS-SECTIONS AND DISTRIBUTIONS -----------
C --- CAN BE CHANGED --------------------------------------------
C ---------------------------------------------------------------

      do img=160,1000,20
         do ims=160,1000,20

            mg = dble(img)
            ms = dble(ims)

C***  CALCULATE THE CROSS-SECTIONS IN LO AND NLO
            CALL INTEGGG(RESLO,ERRLO,RESNLO,ERRNLO)

C***  PRINT THE CROSS-SECTIONS IN LO AND NLO AND THEIR RELATIVE ERRORS
            CALL PRIRESGG(RESLO,ERRLO,RESNLO,ERRNLO)

         end do
      end do

c$$$C***  CALCULATE THE DIFFERENTIAL CROSS-SECTION DSIGMA/DPT 
c$$$C***  WITH THE SCALE Q**2 = MG**2 + PT**2
c$$$C***  FOR PT = 50, 100, 150 GEV
c$$$      ITOTAL = 0
c$$$      ISCAPT = 1
c$$$      DO 100 I = 1,3
c$$$         PTMIN = 50.D0 * I
c$$$         PTMAX = PTMIN
c$$$         CALL INTEGGG(RESLO,ERRLO,RESNLO,ERRNLO)
c$$$         CALL PRIRESGG(RESLO,ERRLO,RESNLO,ERRNLO)
c$$$ 100  CONTINUE

      STOP
      END

C**********************************************************************


