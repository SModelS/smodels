C ======================================================================

      SUBROUTINE PRIHEADSS()

C***  PRINT THE ENERGY AND THE COLLIDER TYPE

      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      COMMON/CONST1/S,ENERGY,ALPHAS,MS,MG,MT
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT

      PRINT *
      PRINT 101
      PRINT *
      IF (ICOLL.EQ.0) PRINT 102,ENERGY
      IF (ICOLL.EQ.1) PRINT 103,ENERGY
      PRINT *
      PRINT 104
      PRINT 105
      PRINT 111

 101  FORMAT('CROSS-SECTIONS AND DISTRIBUTIONS FOR ',
     +     'SQUARK-SQUARK HADROPRODUCTION (IN PB)') 
 102  FORMAT('PROTON-ANTIPROTON COLLIDER WITH ENERGY = ',F6.0, ' GEV')
 103  FORMAT('PROTON-PROTON COLLIDER WITH ENERGY = ',F6.0,' GEV')
 104  FORMAT('IFLAVOR   SCAFAC ISCAPT                              ',
     +     'ERR/SIG           ERR/SIG')
 105  FORMAT('|   MS    MG   |  | PTMIN PTMAX YMIN YMAX   SIGMA_LO ',
     +     '   |     SIGMA_NLO   |')
 111  FORMAT('--------------------------------------------------------',
     +     '---------------------')

      RETURN
      END

C ======================================================================

      SUBROUTINE PRIRESSS(RESLO,ERRLO,RESNLO,ERRNLO)

C***  PRINT THE INITIAL STATE, SQUARK AND GLUINO MASSES, SCAFAC,ISCAPT,
C***  THE CUTS IN PT AND Y
C***  THE CROSS-SECTIONS IN LO AND NLO AND THE RELATIVE ERROR

      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      COMMON/CONST1/S,ENERGY,ALPHAS,MS,MG,MT
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT
      COMMON/CONST3/IPDFSET
      COMMON/CUT1/PTMIN,PTMAX
      COMMON/CUT2/YMIN,YMAX
      COMMON/FLAVOR/IFLAVOR,ITOTAL

      IF (ITOTAL.EQ.0)
     +     PRINT 111,IFLAVOR,MS,MG,SCAFAC,iscapt,PTMIN,PTMAX,
     +     YMIN,YMAX,RESLO,ERRLO,RESNLO,ERRNLO
      IF (ITOTAL.EQ.1)
     +     PRINT 111,IFLAVOR,MS,MG,SCAFAC,0,0D0,ENERGY,
     +     0D0,9.99D0,RESLO,ERRLO,RESNLO,ERRNLO

 111  FORMAT(I1,' ',2(F5.0,' '),F3.1,' ',I1,' ',F4.0,' ',
     +     F6.0,' ',2(F4.2,' '),'  ',
     +     G10.4,' ',F5.4,'  ',G10.4,' ',F5.4) 

      RETURN
      END

C ======================================================================

