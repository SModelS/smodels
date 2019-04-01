c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c        
c                         NNLL-fast,  version 1.1
c
c     Squark and gluino cross sections at NNLL+NNLO_Approx, including
c     NRQCD resummation of Coulomb corrections and bound state effects
c                               at LHC@13 TeV      
c
c A simple interpolation code, based on polynomial interpolation subroutines from
c Numerical Recipes, needs input of data grids. In the current version all data grids
c are for LHC production cross sections at 13 TeV collision energy. 
c
c All results in this version are obtained with PDF4LHC15 (mc) pdf sets.
c
c NNLL numbers obtained on the basis of work described in: 
c     *A. Kulesza, L. Motyka, Phys. Rev. Lett. 102 (2009) 111802, Phys. Rev. D80 (2009) 095004;
c     *W. Beenakker, S. Brensing, M. Krämer, A. Kulesza, E. Laenen and I. Niessen, JHEP 0912 (2009) 041;
c      JHEP 1008 (2010) 098; JHEP 1201 (2012) 076;
c     *W. Beenakker, T. Janssen, S. Lepoeter, M. Krämer, A. Kulesza, E. Laenen, I. Niessen, S. Thewes, T. Van Daal,
c      JHEP 1310 (2013) 120;
c     *W. Beenakker, Ch. Borschensky, M. Krämer, A. Kulesza, E. Laenen, V. Theeuwes, S. Thewes, JHEP 1412 (2014) 023;
c     *W. Beenakker, Ch. Borschensky, R. Heger, M. Krämer, A. Kulesza, E. Laenen, JHEP 1605 (2016) 153; 
c     *W. Beenakker, Ch. Borschensky, M. Krämer, A. Kulesza, E. Laenen, JHEP 1612 (2016) 133;
c     
c      
c NLO numbers (a part of NNLO_Approx) are the output of PROSPINO code:
c   *W. Beenakker, R. Hopker, M. Spira and P. M. Zerwas, Nucl. Phys. B492 (1997) 51.
c   *W. Beenakker, M. Krämer, T. Plehn, M. Spira, P.M. Zerwas, Nucl.Phys. B515 (1998) 3.
c
c After compiling the code, it should be called in the following way for all four squark and gluino
c production processes:
c
c ./name_of_the_executable process squark_mass gluino_mass
c
c or, for the case of gluino production in the limit of large squark masses (decoupling limit):
c
c ./name_of_the_executable gdcpl gluino_mass
c
c or, for the case of squark production in the limit of large gluino masses (decoupling limit):
c
c ./name_of_the_executable sdcpl squark_mass
c
c
c The command line attributes above can take the following values: 
c  - process: 
c      sb     squark-antisquark production,
c      ss     squark-squark production,
c      gg     gluino-gluino production,
c      sg     squark-gluino production,
c      st     stop-antistop (sbottom-antisbottom) production,
c  - masses for which the NNLL+NNLO_Approx cross section can be calculated:
c    for gluino-pair production:
c       500 GeV <= gluino mass <= 3000 GeV   and   500 GeV <= squark mass <= 3000 GeV,
c    for gluino-pair production in the decoupling limit (large squark mass):
c       500 GeV <= gluino mass <= 3000 GeV,
c    for squark-gluino  production:
c       500 GeV <= gluino mass <= 3000 GeV   and   500 GeV <= squark mass <= 3000 GeV,
c    for squark-antisquark and squark-pair production:
c       500 GeV <= squark mass <= 3000 GeV   and   500 GeV <= gluino mass <= 3000 GeV,
c    for squark production in the decoupling limit (large gluino mass):
c       500 GeV <= squark mass <= 3000 GeV,
c    and for stop-antistop production:
c       100 GeV <=  stop mass  <= 3000 GeV   and   500 GeV <= gluino mass <= 5000 GeV 
c       (other SUSY parameters: light flavour squarks and stop-2 considered as very heavy,
c                               stop mixing angle set according to CMSSM benchmark 
c                               point 40.2.5 as defined in arXiv:1109.3859.)  
c
c
c  After reading in the mass values, the code returns (in order of appearance) the values of:
c  - squark (stop) mass, gluino mass,
c  - NLO and NNLL+NNLO_Approx cross sections in pb, 
c  - upper and lower error on NNLL+NNLO_Approx due to scale mu variation between m/2 < mu < 2m, 
c    where m is the average mass of the pair-produced particles, also in pb
c  - combined upper and lower PDF+alpha_s uncertainty in percent (of the value for the NNLL+NNLO_Approx cross section),
c  - NNLL K-factor, with the NNLL K-factor defined as sigma_NNLL+NNLO_Approx/sigma_NLO,
c  - for stop production cross sections, additional upper and lower theoretical uncertainty
c    due to the residual dependence on the stop mixing angle around its central value.
c
c Comments:
c -- Results for production processes with a squark in the final state are obtained with 
c    5 squark flavours (mass degenerate), i.e. include contributions from sbottoms which
c    are treated as a light flavour
c -- Tests indicate interpolation accuracy indicate an error of up to ~2% for processes other than
c    gluino-pair production where the interpolation error can reach 4%. These maximal values are 
c    observed only at the edge of the grid (small masses of the final state), therefore not relevant
c    for most phenomenological applications. Away from these regions the interpolation accuracy is 
c    better than 1%.  
c -- Currently not a part of the standard output, but the code can also return interpolation errors 
c    on each quantity
c -- No leading order results are provided since PDF4LHC15 does not contain a LO set (the numbers
c    given in the gris are either calculated with NLO or NNLO pdfs)     
c
c 
c AK, 22/04/16
c DS, 18/07/16
c AK, 19/07/16
c CB, 21/07/16
c CB, 25/07/16
c CB, 05/04/17
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      PROGRAM NNLLFAST

      IMPLICIT NONE

      INTEGER NITP,NSQ,NGL,I,J,KSQ,KINSQ,KGL,KINGL
      INTEGER NSQMAX,NGLMAX
      REAL*8 MSQ,MGL,MSQGR,MGLGR,LOGR,LOGRPL,LOGRMI,
     &     NLOGR,NLOGRPL,NLOGRMI,PDFGRPL,PDFGRMI,ASGRPL,ASGRMI,
     &     LNLOGR,NLO,DNLO,PDFPL,DPDFPL,
     &     NNLLGR,NNLLGRPL,NNLLGRMI,LNNLLGR,LNNLLGRPL,LNNLLGRMI,NNLL,
     &     NNLLPL,NNLLMI,DNNLL,DNNLLPL,DNNLLMI,DMUPL,DMUMI,NLOPARPL,
     &     NLOPARMI,NNLLPARPL,NNLLPARMI,
     &     PARPL,DPARPL,PARMI,DPARMI
      REAL*8 MSQ_MIN,MSQ_MAX,MGL_MIN,MGL_MAX
 
      CHARACTER*5 CR,PR
    

      DIMENSION MSQGR(100),MGLGR(100),LOGR(100,100),LOGRPL(100,100),
     &     LOGRMI(100,100),NLOGR(100,100),NLOGRPL(100,100),
     &     NLOGRMI(100,100),PDFGRPL(100,100),PDFGRMI(100,100),
     &     ASGRPL(100,100),ASGRMI(100,100),
     &     LNLOGR(100,100),NNLLGR(100,100),
     &     NNLLGRPL(100,100), NNLLGRMI(100,100),LNNLLGR(100,100),
     &     LNNLLGRPL(100,100),LNNLLGRMI(100,100),NLOPARPL(100,100),
     &     NLOPARMI(100,100),NNLLPARPL(100,100),NNLLPARMI(100,100)


      COMMON/INPUT/PR
      COMMON/MASSES/MSQ,MGL
      COMMON/GRIDSIZE/NSQ,NGL
      COMMON/LIMITS/MSQ_MIN,MSQ_MAX,MGL_MIN,MGL_MAX 
      COMMON/INTERPLEVEL/NITP


c... initialization
      CALL INOUT(0)


c... initialize grids
       DO I=1,100
          MSQGR(I)=0.D0
          MGLGR(I)=0.D0
          DO J=1,100
              LOGR(I,J)=0.d0
              LOGRPL(I,J)=0.d0
              LOGRMI(I,J)=0.d0
              NLOGR(I,J)=0.d0
              NLOGRPL(I,J)=0.d0
              NLOGRMI(I,J)=0.d0
              NLOPARPL(I,J)=0.d0
              NLOPARMI(I,J)=0.d0
              LNLOGR(I,J)=0.d0
              PDFGRPL(I,J)=0.d0
              PDFGRMI(I,J)=0.d0
              ASGRPL(I,J)=0.d0
              ASGRMI(I,J)=0.d0
              NNLLGR(I,J)=0.d0
              NNLLGRPL(I,J)=0.d0
              NNLLGRMI(I,J)=0.d0
              NNLLPARPL(I,J)=0.d0
              NNLLPARMI(I,J)=0.d0
              LNNLLGR(I,J)=0.d0
              LNNLLGRPL(I,J)=0.d0
              LNNLLGRMI(I,J)=0.d0

           ENDDO
        ENDDO 

             
c... open grid files
      CALL INOUT(1)


c... read in data

        IF (PR.EQ.'st') THEN
            read(21,*)            
            DO I=1,NSQ
              DO J=1,NGL
                read(21,*) CR,MSQGR(I),MGLGR(J),LOGR(I,J),LOGRPL(I,J),
     &          LOGRMI(I,J),NLOGR(I,J),NLOGRPL(I,J),NLOGRMI(I,J),
     &          PDFGRPL(I,J),NLOPARPL(I,J),NLOPARMI(I,j),NNLLGR(I,J),
     &          NNLLGRPL(I,J),NNLLGRMI(I,J),ASGRPL(I,J),NNLLPARPL(I,J),
     &          NNLLPARMI(I,J)

                  PDFGRPL(I,J) = ASGRPL(I,J)
  
                  LNLOGR(I,J)=DLOG(NLOGR(I,J))
                  LNNLLGR(I,J)=DLOG(NNLLGR(I,J))
                  LNNLLGRPL(I,J)=DLOG(NNLLGR(I,J)+NNLLGRPL(I,J))
                  LNNLLGRMI(I,J)=DLOG(NNLLGR(I,J)+NNLLGRMI(I,J))
                ENDDO
            ENDDO   
        ELSE    
            read(21,*)            
            DO I=1,NSQ
                DO J=1,NGL
                  read(21,*) CR,MSQGR(I),MGLGR(J),LOGR(I,J),LOGRPL(I,J),
     &                 LOGRMI(I,J),NLOGR(I,J),NLOGRPL(I,J),NLOGRMI(I,J),
     &                 PDFGRPL(I,J),NNLLGR(I,J),NNLLGRPL(I,J),
     &                  NNLLGRMI(I,J),ASGRPL(I,J)

                  PDFGRPL(I,J) = ASGRPL(I,J)

                  LNLOGR(I,J)=DLOG(NLOGR(I,J))
                  LNNLLGR(I,J)=DLOG(NNLLGR(I,J))
                  LNNLLGRPL(I,J)=DLOG(NNLLGR(I,J)+NNLLGRPL(I,J))
                  LNNLLGRMI(I,J)=DLOG(NNLLGR(I,J)+NNLLGRMI(I,J))

                ENDDO
            ENDDO   
          ENDIF  



c...define level of interpolation
           CALL INOUT(2)

 
            
c... find the appropriate points on the grid

           NSQMAX=NSQ
           NGLMAX=NGL
   
     
           IF (PR.EQ.'gdcpl') THEN
              CALL LOCATE(MGLGR,NGLMAX,MGL,KGL)   
           ELSEIF (PR.EQ.'sdcpl') THEN
              CALL LOCATE(MSQGR,NSQMAX,MSQ,KSQ)
           ELSE
              CALL LOCATE(MGLGR,NGLMAX,MGL,KGL)
              CALL LOCATE(MSQGR,NSQMAX,MSQ,KSQ)
           ENDIF


c... calculate the offset where one has to start reading in the grid

           IF (PR.EQ.'gdcpl') THEN
              KINGL=MIN(MAX(KGL-(NITP-1)/2,1),NGLMAX+1-NITP)
              KINSQ=0
           ELSEIF (PR.EQ.'sdcpl') THEN
              KINGL=0
              KINSQ=MIN(MAX(KSQ-(NITP-1)/2,1),NSQMAX+1-NITP)
           ELSE
              KINGL=MIN(MAX(KGL-(NITP-1)/2,1),NGLMAX+1-NITP)
              KINSQ=MIN(MAX(KSQ-(NITP-1)/2,1),NSQMAX+1-NITP)
           ENDIF   

 
c... interpolate

            
c... for NLO:
            CALL GRIDINTP(MSQGR,MGLGR,LNLOGR,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NLO,DNLO)
            
c... for NNLL+NNLO_Approx (with scale variation):
            CALL GRIDINTP(MSQGR,MGLGR,LNNLLGR,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NNLL,DNNLL)
            CALL GRIDINTP(MSQGR,MGLGR,LNNLLGRPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NNLLPL,DNNLLPL)
            CALL GRIDINTP(MSQGR,MGLGR,LNNLLGRMI,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NNLLMI,DNNLLMI)

c... for pdf and as errors:
            CALL GRIDINTP(MSQGR,MGLGR,PDFGRPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,PDFPL,DPDFPL)

c... for stops for parameter variation:
            IF (PR.EQ.'st') THEN
             CALL GRIDINTP(MSQGR,MGLGR,NNLLPARPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,PARPL,DPARPL)
             CALL GRIDINTP(MSQGR,MGLGR,-NNLLPARMI,NITP,KINSQ,KINGL,
     &           MSQ,MGL,PARMI,DPARMI)
            ENDIF

            
c...  check the signs of the scale error:
            IF (NNLLPL.GE.NNLLMI) THEN
                     DMUPL=NNLLPL
                     DMUMI=NNLLMI
            ELSE
                     DMUPL=NNLLMI
                     DMUMI=NNLLPL
            ENDIF

               
c...  write out results:
        IF (PR.EQ.'st') THEN    
            WRITE(*,15) MSQ,MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           PARPL,-PARMI,DEXP(NNLL-NLO)
            WRITE(23,16) MSQ,MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           PARPL,-PARMI,DEXP(NNLL-NLO)

        ELSEIF (PR.EQ.'gdcpl') THEN    

            WRITE(*,17) MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)
            WRITE(23,18) MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)

        ELSEIF (PR.EQ.'sdcpl') THEN    

            WRITE(*,17) MSQ,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)
            WRITE(23,18) MSQ,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)

        ELSE

            WRITE(*,19) MSQ,MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)
            WRITE(23,20) MSQ,MGL,DEXP(NLO),DEXP(NNLL),
     &           (DEXP(DMUPL-NNLL)-1d0)*100.d0,
     &           (DEXP(DMUMI-NNLL)-1d0)*100.d0,PDFPL,-PDFPL,
     &           DEXP(NNLL-NLO)

        ENDIF    


 
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)


 15   FORMAT (' ',2(F7.0,2x),(E11.3),'   ',(E11.3),'     ',(f8.2),
     &     '  ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',(f8.2),
     &     '    ',(f8.2),' ',(f8.2))
 16   FORMAT (' ',2(F7.0,2x),(E11.3),'   ',(E11.3),'     ',(f8.2),
     &     '  ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',(f8.2),
     &     '    ',(f8.2),' ',(f8.2))
 17   FORMAT (' ',1(F7.0,2x),(E11.3),'     ',(E11.3),'      ',(f8.2),
     &     '    ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',1(f8.2))
 18   FORMAT (' ',1(F7.0,2x),(E11.3),'     ',(E11.3),'      ',(f8.2),
     &     '    ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',1(f8.2))
 19   FORMAT (' ',2(F7.0,2x),(E11.3),'     ',(E11.3),'      ',(f8.2),
     &     '    ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',1(f8.2))
 20   FORMAT (' ',2(F7.0,2x),(E11.3),'     ',(E11.3),'      ',(f8.2),
     &     '    ',(f8.2),'   ',(f8.2),'     ',(f8.2),'    ',1(f8.2))
      
      END

c.......................................................................

      SUBROUTINE INOUT(MODE)
      
      IMPLICIT NONE
      
      INTEGER MODE,NSQ,NGL
      
      REAL*8 MSQ,MGL,MSQLO,MSQUP,MGLLO,MGLUP
 
      CHARACTER*5 PR,MSQCHAR,MGLCHAR
      INTEGER NITP
      COMMON/INPUT/PR

      COMMON/MASSES/MSQ,MGL
      COMMON/GRIDSIZE/NSQ,NGL
      COMMON/LIMITS/MSQLO,MSQUP,MGLLO,MGLUP
      COMMON/INTERPLEVEL/NITP
      


      IF (MODE.EQ.0) THEN

      IF (IARGC().LT.2.OR.IARGC().GT.3) THEN
       WRITE(*,*) 'Usage:'
       WRITE(*,*) '  ./nnllfast process squark_mass gluino_mass'
       WRITE(*,*) 'or'
       WRITE(*,*) '  ./nnllfast gdcpl gluino_mass'
       WRITE(*,*) 'or'
       WRITE(*,*) '  ./nnllfast sdcpl squark_mass'
       WRITE(*,*)
       WRITE(*,*) 'where:'
       WRITE(*,*) ' - nnllfast is the name of the executable'
       WRITE(*,*) ' - process is one of the following:'
       WRITE(*,*) '    gg     gluino-gluino production'
       WRITE(*,*) '    sb     squark-antisquark production'
       WRITE(*,*) '    sg     squark-gluino production'
       WRITE(*,*) '    ss     squark-squark production'
       WRITE(*,*) '    st     stop-antistop (sbottom-antisbottom) produc
     &tion'
       WRITE(*,*) '    gdcpl  gluino production in the heavy squark limi
     &t'
       WRITE(*,*) '    sdcpl  squark production in the heavy gluino limi
     &t'
       WRITE(*,*) ' - squark_mass is the squark mass (for gg, sb, sg, ss
     &, sdcpl) or the light stop mass (for st)'
       WRITE(*,*) ' - gluino_mass is the gluino mass'
       WRITE(*,*)
       WRITE(*,*) 'The masses are in the following ranges:'
       WRITE(*,*) ' - for gg, sb, sg, ss:'
       WRITE(*,*) '    500 GeV <= squark mass <= 3000 GeV   and   500 Ge
     &V <= gluino mass <= 3000 GeV'
!       WRITE(*,*) '                   and'
!       WRITE(*,*) '    500 GeV <= gluino mass <= 3000 GeV'
       WRITE(*,*) ' - for st:'
       WRITE(*,*) '    100 GeV <=  stop mass  <= 3000 GeV   and   500 Ge
     &V <= gluino mass <= 5000 GeV'
!       WRITE(*,*) '                   and'
!       WRITE(*,*) '    500 GeV <= gluino mass <= 5000 GeV'
       WRITE(*,*) ' - for gdcpl:'
       WRITE(*,*) '    500 GeV <= gluino mass <= 3000 GeV'
       WRITE(*,*) ' - for sdcpl:'
       WRITE(*,*) '    500 GeV <= squark mass <= 3000 GeV'
       STOP
      END IF
         

c... read in command-line parameters
      
      CALL GETARG(1,PR)
      CALL GETARG(2,MSQCHAR)
      READ(MSQCHAR,*) MSQ
      IF ((PR.NE.'sdcpl').AND.(PR.NE.'gdcpl')) THEN
         CALL GETARG(3,MGLCHAR)
         READ(MGLCHAR,*) MGL
      ENDIF   
      IF (PR.EQ.'gdcpl') THEN
         MGL=MSQ
         MSQ=0d0
      ELSEIF (PR.EQ.'sdcpl') THEN
         MGL=0d0
      ENDIF
   

c... define lower and upper limits on squark and gluino masses

      IF (PR.EQ.'st') THEN
         MSQLO=1.d2
         MSQUP=3.d3
         MGLLO=5.d2
         MGLUP=5.d3
      ELSEIF (PR.EQ.'gdcpl') THEN
         MSQLO=0d0
         MSQUP=0d0
         MGLLO=5.d2
         MGLUP=3.d3
      ELSEIF (PR.EQ.'sdcpl') THEN
         MSQLO=5.d2
         MSQUP=3.d3
         MGLLO=0d0
         MGLUP=0d0
      ELSE
         MSQLO=5.d2
         MSQUP=3.d3
         MGLLO=5.d2
         MGLUP=3.d3
      ENDIF   


c... check the masses 

        IF (MSQ.LT.MSQLO) THEN
            write(*,*) 'TOO LOW SQUARK MASS'
            STOP
        ELSEIF (MSQ.GT.MSQUP) THEN
            write(*,*) 'TOO HIGH SQUARK MASS'
            STOP
        ELSEIF (MGL.LT.MGLLO) THEN
            write(*,*) 'TOO LOW GLUINO MASS'
            STOP
        ELSEIF (MGL.GT.MGLUP) THEN
            write(*,*) 'TOO HIGH GLUINO MASS'
            STOP
        ENDIF 
            
        
 
c... open output file

         IF (PR.EQ.'sb') THEN
            open(23,file="sb.out",status='unknown') 
         ELSEIF (PR.EQ.'ss') THEN
            open(23,file="ss.out",status='unknown') 
         ELSEIF (PR.EQ.'gg') THEN
            open(23,file="gg.out",status='unknown') 
         ELSEIF (PR.EQ.'sg') THEN
            open(23,file="sg.out",status='unknown') 
         ELSEIF (PR.EQ.'st') THEN
            open(23,file="st.out",status='unknown') 
         ELSEIF (PR.EQ.'gdcpl') THEN
            open(23,file="gdcpl.out",status='unknown') 
         ELSEIF (PR.EQ.'sdcpl') THEN
            open(23,file="sdcpl.out",status='unknown') 
         ELSE
            write(*,*) 'WRONG INPUT PROCESS'
            STOP
         ENDIF
         
c... write header for output 
c... stops:
       
      IF (PR.EQ.'st') THEN
         write(*,*)
         write(*,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for N
     &LO)'
         write(*,*) '# process: ', PR
         write(*,*) '# mst[GeV] mg[GeV]   NLO[pb]   NNLL+NNLO_app[pb]  d
     &_mu+[%]  d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]  d_par+[%]  d_par-[%]
     & K_NNLL'
         write(*,*) '---------------------------------------------------
     &------------------------------------------------------------------
     &-------'

         write(23,*)
         write(23,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for 
     &NLO)'
         write(23,*) '# process: ', PR
         write(23,*) '# mst[GeV] mg[GeV]   NLO[pb]   NNLL+NNLO_app[pb]  
     &d_mu+[%]  d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]  d_par+[%]  d_par-[%]
     &  K_NNLL'
         write(23,*) '--------------------------------------------------
     &------------------------------------------------------------------
     &-------'

c... gluinos in the squark decoupling limit:
       ELSEIF (PR.EQ.'gdcpl') THEN
         write(*,*)
         write(*,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for N
     &LO)'
         write(*,*) '# process: ', PR
         write(*,*) '# mg[GeV]   NLO[pb]     NNLL+NNLO_app[pb]   d_mu+[%
     &]    d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]   K_NNLL'

         write(*,*)'----------------------------------------------------
     &----------------------------------------------------------'

         write(23,*)
         write(23,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for 
     &NLO)'
         write(23,*) '# process: ', PR
         write(23,*) '# mg[GeV]  NLO[pb]     NNLL+NNLO_app[pb]   d_mu+[%
     &]    d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]   K_NNLL'

         write(23,*)'---------------------------------------------------
     &-----------------------------------------------------------'

c... squarks in the gluino decoupling limit:
       ELSEIF (PR.EQ.'sdcpl') THEN
         write(*,*)
         write(*,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for N
     &LO)'
         write(*,*) '# process: ', PR
         write(*,*) '# ms[GeV]   NLO[pb]     NNLL+NNLO_app[pb]   d_mu+[%
     &]    d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]   K_NNLL'

         write(*,*)'----------------------------------------------------
     &----------------------------------------------------------'

         write(23,*)
         write(23,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for 
     &NLO)'
         write(23,*) '# process: ', PR
         write(23,*) '# ms[GeV]  NLO[pb]     NNLL+NNLO_app[pb]   d_mu+[%
     &]    d_mu-[%]  d_pdfas+[%]  d_pdfas-[%]   K_NNLL'

         write(23,*)'---------------------------------------------------
     &-----------------------------------------------------------'

c... other processes:
       ELSE
         write(*,*)
         write(*,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for N
     &LO)'
         write(*,*) '# process: ', PR
         write(*,*) '# ms[GeV]  mg[GeV]   NLO[pb]     NNLL+NNLO_app[pb]
     &  d_mu+[%]    d_mu-[%]   d_pdfas+[%]  d_pdfas-[%]  K_NNLL'

         write(*,*)'----------------------------------------------------
     &----------------------------------------------------------'

         write(23,*)
         write(23,*) '# LHC @ 13 TeV, NNLO PDF4LHC15 (NLO PDF4LHC15 for 
     &NLO)'
         write(23,*) '# process: ', PR
         write(23,*) '# ms[GeV]  mg[GeV]   NLO[pb]     NNLL+NNLO_app[pb]
     &  d_mu+[%]    d_mu-[%]   d_pdfas+[%]  d_pdfas-[%]  K_NNLL'

         write(23,*)'---------------------------------------------------
     &-----------------------------------------------------------'

       ENDIF   



      ELSE IF (MODE.EQ.1) THEN 

c... open the grid file

         IF (PR.EQ.'sb') THEN
            open(21,file="sb_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'ss') THEN
            open(21,file="ss_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'gg') THEN
            open(21,file="gg_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'sg') THEN
            open(21,file="sg_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'st') THEN
            open(21,file="st_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'gdcpl') THEN
            open(21,file="gdcpl_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSEIF (PR.EQ.'sdcpl') THEN
            open(21,file="sdcpl_nnlonnll_pdf4lhc15_13TeV_wpresc.grid",
     &           status='unknown')
         ELSE   
            write(*,*) 'WRONG INPUT PROCESS'
            STOP
         ENDIF


c... define the size of the grid 

         IF (PR.EQ.'st') THEN
           NSQ=30
           NGL=46
         ELSEIF (PR.EQ.'gdcpl') THEN
           NSQ=1
           NGL=26
         ELSEIF (PR.EQ.'sdcpl') THEN
           NSQ=26
           NGL=1
         ELSE  
           NSQ=26
           NGL=26
         ENDIF 


      ELSE IF (MODE.EQ.2) THEN 
         

c... define level of interpolation   
      
        IF((PR.EQ.'ss').OR.(PR.EQ.'sg')
     &       .AND.(MSQ.GE.MSQLO).AND.(MSQ.LT.8.d2)) THEN
         NITP=6
        ELSEIF ((PR.EQ.'gg').AND.(MGL.GE.MGLLO).AND.(MGL.LT.8.d2)) THEN
         NITP=6 
        ELSE 
         NITP=4
        ENDIF

      ENDIF
      
      RETURN
      END

c.......................................................................

      SUBROUTINE GRIDINTP(MSQGR,MGLGR,DATAGR,NITP,KINSQ,KINGL,
     &                MSQ,MGL,INTPOUT,DINTPOUT)

      IMPLICIT NONE
      INTEGER NITP,KINSQ,KINGL,I,J
      REAL*8 MSQGR,MGLGR,DATAGR,MSQINGR,MGLINGR,DATAINGR,
     &       MSQ,MGL,INTPOUT,DINTPOUT,DATAINGRDCPL
      DIMENSION MSQGR(100),MGLGR(100),DATAGR(100,100)
      DIMENSION MSQINGR(NITP),MGLINGR(NITP),DATAINGR(NITP,NITP),
     &       DATAINGRDCPL(NITP)

      CHARACTER*5 PR
      COMMON/INPUT/PR

c... copy the grid

      IF (PR.EQ.'gdcpl') THEN
        DO I=1,NITP
          MGLINGR(I)=MGLGR(KINGL+I-1)
          DATAINGRDCPL(I)=DATAGR(1,KINGL+I-1)
        ENDDO  
c... interpolate (1-D grid)
        CALL POLINT(MGLINGR,DATAINGRDCPL,NITP,MGL,INTPOUT,DINTPOUT)
      ELSEIF (PR.EQ.'sdcpl') THEN
        DO I=1,NITP
          MSQINGR(I)=MSQGR(KINSQ+I-1) 
          DATAINGRDCPL(I)=DATAGR(KINSQ+I-1,1)
        ENDDO
c... interpolate (1-D grid)
        CALL POLINT(MSQINGR,DATAINGRDCPL,NITP,MSQ,INTPOUT,DINTPOUT)
      ELSE
        DO I=1,NITP
          MSQINGR(I)=MSQGR(KINSQ+I-1) 
          DO J=1,NITP
            MGLINGR(J)=MGLGR(KINGL+J-1)
            DATAINGR(I,J)=DATAGR(KINSQ+I-1,KINGL+J-1)
          ENDDO
        ENDDO
c... interpolate (2-D grid)
        CALL POLIN2(MSQINGR,MGLINGR,DATAINGR,
     &       NITP,NITP,MSQ,MGL,INTPOUT,DINTPOUT)
      ENDIF  

      RETURN
      END

c.......................................................................
c.... taken from Numerical Recipes:

      SUBROUTINE POLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY)
      
      IMPLICIT NONE
      INTEGER M,N,J,K,NMAX,MMAX
      REAL*8 X1A,X2A,YA,X1,X2,Y,DY,YNTMP,YMTMP
      PARAMETER (NMAX=20,MMAX=20)
      DIMENSION X1A(M),X2A(N),YA(M,N),YNTMP(NMAX),YMTMP(MMAX)

      DO 12 J=1,M
        DO 11 K=1,N
          YNTMP(K)=YA(J,K)
c          write(*,*) M,N,J,K,YNTMP(K),YA(J,K)
11      CONTINUE
        CALL POLINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
c         CALL RATINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
12    CONTINUE
      CALL POLINT(X1A,YMTMP,M,X1,Y,DY)
c      CALL RATINT(X1A,YMTMP,M,X1,Y,DY)
      RETURN
      END


      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

      IMPLICIT NONE
      INTEGER N,NS,M,I,NMAX
      REAL*8 XA,YA,X,Y,DY,DIF,DIFT,C,D,HO,HP,W,DEN

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
CMK          IF(DEN.EQ.0.)PAUSE
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

      SUBROUTINE LOCATE(XX,N,X,J)

      IMPLICIT NONE
      INTEGER N,J,JL,JU,JM
      REAL*8 XX,X
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END
