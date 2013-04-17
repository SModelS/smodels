c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c        
c                         NLL-fast,  version 2.1
c
c                 Squark and gluino (and stop) cross sections at NLL+NLO 
c
c A simple interpolation code, based on polynomial interpolation subroutines from
c Numerical Recipes, needs input of data grids. In the current version all data grids
c are for LHC production cross sections at 8 TeV collision energy.
c
c NLL numbers obtained on basis of work described in: 
c   *A. Kulesza, L. Motyka, Phys. Rev. Lett. 102(2009)111802, Phys. Rev. D80(2009)095004;
c   *W. Beenakker, S. Brensing, M. Kramer, A. Kulesza, E. Laenen and I. Niessen, JHEP 0912(2009)041, JHEP 1008(2010)098. 
c NLO numbers are the output of PROSPINO code:
c   *W. Beenakker, R. Hopker, M. Spira and P. M. Zerwas, Nucl. Phys. B492(1997)51.
c   *W. Beenakker, M. Krämer, T. Plehn, M. Spira, P.M. Zerwas, Nucl.Phys. B515 (1998) 3.
c
c After compiling the code, it should be called in the following way for all four squark and gluino
c production processes:
c
c ./name_of_the_executable process pdf squark_mass gluino_mass
c
c or, for the case of stop-antistop (sbottom-antisbottom) production:
c
c ./name_of_the_executable st pdf stop_mass
c
c or, for the case of gluino production in the limit of large squark masses (decoupling limit):
c
c ./name_of_the_executable gdcpl pdf gluino_mass
c
c or, for the case of squark production in the limit of large gluino masses (decoupling limit):
c
c ./name_of_the_executable sdcpl pdf squark_mass
c
c The command line attributes above can take the following values: 
c  - process: 
c      sb     squark-antisquark production,
c      ss     squark-squark production,
c      gg     gluino-gluino production,
c      sg     squark-gluino production,
c  - pdf:
c      cteq   CTEQ6.6M NLO (CTEQ6L1 for LO predictions) parton distribution functions,
c      mstw   MSTW2008NLO  (MSTW2008LO for LO predictions) parton distribution functions,
c  - masses for which the NLL+NLO cross section can be calculated:
c    for gluino-pair production:
c       200 GeV  <= gluino mass <= 2500 GeV    and   200 GeV <= squark mass <= 4500 GeV,
c    for gluino production in the decoupling limit (large squark mass)
c       200 GeV  <= gluino_mass <= 2500 GeV 
c    for squark-gluino  production:
c       200 GeV  <= gluino mass <= 2500 GeV    and   200 GeV <= squark mass <= 4500 GeV,
c    for squark-antisquark and squark-pair production:
c       200 GeV  <= squark mass <= 2500 GeV    and   200 GeV <= gluino mass <= 2500 GeV
c    for squark production in the decoupling limit (large gluino mass)
c       200 GeV  <= squark mass <= 2500 GeV
c    and for stop-antistop production:
c       100 GeV  <= stop mass   <= 2000 GeV  and stop_mass <= gluino mass 
c                                             (other SUSY parameters according to CMSSM benchmark 
c                                                    point 40.2.5 as defined in arXiv:1109.3859.)  
c
c  NOTE: option tot (active in versions 1.x) is disabled in this version.
c
c  After reading in the mass values, the code returns (in order of appearance) the values of:
c  - squark_mass, gluino_mass,
c  - LO, NLO and NLL+NLO cross sections in pb, 
c  - upper and lower error on NLL+NLO due to scale mu variation between m/2 < mu < 2m, 
c    where m is the average mass of the pair-produced particles, also in pb
c  - upper and lower 68% C.L. pdf error in percent (of the NLO value),
c  - upper and lower error due to alpha_s uncertainty, in percent (of the NLO value),
c  - NLO and NLL K-factors, with the NLL K-factor defined as sigma_NLL+NLO/sigma_NLO.
c
c Comments:
c -- Results for production processes with a squark in the final state are obtained with 
c    5 squark flavours (mass degenerate), i.e. include contributions from sbottoms which
c    are treated as a light flavour
c -- To obtain cross sections for sbottom-antisbottom production only, please use the 'st'   
c    option with sbottom_mass instead of stop_mass. (Contributions from initial-state bottom 
c    quarks are numerically negligible.)
c -- Test indicate interpolation accuracy indicate an error < 1%, apart from the case of very low
c    squark masses (squark mass close to 200 GeV) for squark-antisquark and squark-squark production,
c    and similarly gluino masses close to 200 GeV for gluino-pair production,
c    where the interpolation error can reach 2-3%.
c -- For stop-pair production the error due to variation of the mixing angle, gluino mass and the 
c    average mass of light-flavour squarks from the values corresponding to the point 40.2.5. can reach 
c    up to 4% when gluino mass is biger than stop_1 mass.
c -- Currently not a part of the standard output, but the code can also return:
c             - interpolation errors on each quantity
c             - upper and lower scale variation errors for LO and NLO (needed lines
c               in the code are commented out). 
c 
c AK, 05/08/11
c MK, 28/08/11
c AK, 28/10/11
c AK, 01/12/11
c AK, 30/12/11
c AK, 16/01/12
c AK, 02/02/12
c AK, 13/06/12
c AK, 30/06/12
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      PROGRAM NLLFAST
      IMPLICIT NONE

      INTEGER NITP,NSQ,NGL,NSQHM,NGLHM,I,J,IHM,KSQ,KINSQ,KGL,KINGL
      INTEGER I_SQ,J_GL,I_MSQ,J_MGL,NSQMAX,NGLMAX
      REAL*8 MSQ,MGL,MST,MSQGR,MGLGR,DUM,LOGR,LOGRPL,LOGRMI,
     &     NLOGR,NLOGRPL,NLOGRMI,PDFGRPL,PDFGRMI,ASGRPL,ASGRMI,
     &     LLOGR,LLOGRPL,LLOGRMI,LNLOGR,LNLOGRPL,LNLOGRMI,
     &     LO,DLO,LOPL,DLOPL,LOMI,DLOMI,NLO,DNLO,NLOPL,DNLOPL,NLOMI,
     &     DNLOMI,PDFPL,DPDFPL,PDFMI,DPDFMI,ASPL,DASPL,ASMI,DASMI,
     &     NLLGR,NLLGRPL,NLLGRMI,LNLLGR,LNLLGRPL,LNLLGRMI,NLL,NLLPL,
     &     NLLMI,DNLL,DNLLPL,DNLLMI
      REAL*8 MSQ_MIN,MSQ_MAX,MGL_MIN,MGL_MAX,MSQ_STEP,MGL_STEP
 
      REAL*8 MSQ_TEST(1:20), MGL_TEST(1:190)

      REAL*8 Y1(100,100),Y2(100,100),Y12(100,100),G1,G2,D1(4),D2(4),
     &     D12(4)
      CHARACTER*5 CR,PR

      DIMENSION MSQGR(100),MGLGR(100),LOGR(100,100),LOGRPL(100,100),
     &     LOGRMI(100,100),NLOGR(100,100),NLOGRPL(100,100),
     &     NLOGRMI(100,100),PDFGRPL(100,100),PDFGRMI(100,100),
     &     ASGRPL(100,100),ASGRMI(100,100),LLOGR(100,100),
     &     LLOGRPL(100,100),LLOGRMI(100,100),LNLOGR(100,100),
     &     LNLOGRPL(100,100),LNLOGRMI(100,100),NLLGR(100,100),
     &     NLLGRPL(100,100), NLLGRMI(100,100),LNLLGR(100,100),
     &     LNLLGRPL(100,100),LNLLGRMI(100,100)  
      
      COMMON/INPUT/PR
      COMMON/MASSES/MSQ,MGL
      COMMON/GRIDSIZE/NSQ,NGL,NSQHM,NGLHM
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
              LLOGR(I,J)=0.d0
              LLOGRPL(I,J)=0.d0
              LLOGRMI(I,J)=0.d0
              NLOGR(I,J)=0.d0
              NLOGRPL(I,J)=0.d0
              NLOGRMI(I,J)=0.d0
              LNLOGR(I,J)=0.d0
              LNLOGRPL(I,J)=0.d0
              LNLOGRMI(I,J)=0.d0
              PDFGRPL(I,J)=0.d0
              PDFGRMI(I,J)=0.d0
              ASGRPL(I,J)=0.d0
              ASGRMI(I,J)=0.d0
              NLLGR(I,J)=0.d0
              NLLGRPL(I,J)=0.d0
              NLLGRMI(I,J)=0.d0
              LNLLGR(I,J)=0.d0
              LNLLGRPL(I,J)=0.d0
              LNLLGRMI(I,J)=0.d0
           ENDDO
        ENDDO 

             
c... open grid files
      CALL INOUT(1)


c... read in data

            read(21,*)            
            DO I=1,NSQ
                read(21,*)
                read(21,*)
                DO J=1,NGL
                  read(21,*) CR,MSQGR(I),MGLGR(J),LOGR(I,J),LOGRPL(I,J),
     &                 LOGRMI(I,J),NLOGR(I,J),NLOGRPL(I,J),NLOGRMI(I,J),
     &                 PDFGRPL(I,J),PDFGRMI(I,J),ASGRPL(I,J),ASGRMI(I,J)
     &                 ,NLLGR(I,J),NLLGRPL(I,J),NLLGRMI(I,J)
            
  

                  LLOGR(I,J)=DLOG(LOGR(I,J))
                  LLOGRPL(I,J)=DLOG(LOGRPL(I,J))
                  LLOGRMI(I,J)=DLOG(DABS(LOGRMI(I,J)))
                  LNLOGR(I,J)=DLOG(NLOGR(I,J))
                  LNLOGRPL(I,J)=DLOG(NLOGRPL(I,J))
                  LNLOGRMI(I,J)=DLOG(DABS(NLOGRMI(I,J)))
                  LNLLGR(I,J)=DLOG(NLLGR(I,J))
                  LNLLGRPL(I,J)=DLOG(NLLGRPL(I,J))
                  LNLLGRMI(I,J)=DLOG(DABS(NLLGRMI(I,J)))
                ENDDO
            ENDDO   


            IF ((PR.EQ.'gg').OR.(PR.EQ.'sg')) THEN

              read(31,*)

              DO IHM=1,NSQHM 
               I=NSQ+IHM
 
               DO J=1,NGLHM
  
                  read(31,*) CR,MSQGR(I),MGLGR(J),LOGR(I,J),LOGRPL(I,J),
     &                 LOGRMI(I,J),NLOGR(I,J),NLOGRPL(I,J),NLOGRMI(I,J),
     &                 PDFGRPL(I,J),PDFGRMI(I,J),ASGRPL(I,J),ASGRMI(I,J)
     &                 ,NLLGR(I,J),NLLGRPL(I,J),NLLGRMI(I,J)


                  LLOGR(I,J)=DLOG(LOGR(I,J))
                  LLOGRPL(I,J)=DLOG(LOGRPL(I,J))
                  LLOGRMI(I,J)=DLOG(DABS(LOGRMI(I,J)))
                  LNLOGR(I,J)=DLOG(NLOGR(I,J))
                  LNLOGRPL(I,J)=DLOG(NLOGRPL(I,J))
                  LNLOGRMI(I,J)=DLOG(DABS(NLOGRMI(I,J)))
                  LNLLGR(I,J)=DLOG(NLLGR(I,J))
                  LNLLGRPL(I,J)=DLOG(NLLGRPL(I,J))
                  LNLLGRMI(I,J)=DLOG(DABS(NLLGRMI(I,J)))
               ENDDO
              ENDDO 
            ENDIF   



c...define level of interpolation
           CALL INOUT(2)

 
            
c... find the appropriate points on the grid

           IF ((PR.EQ.'gg').OR.(PR.EQ.'sg')) THEN
              IF (MSQ.LE.2.5d3) THEN
                NSQMAX=NSQ
                NGLMAX=NGL
              ELSE  
                NSQMAX=NSQ+NSQHM
                NGLMAX=NGLHM
              ENDIF
           ELSE
              NSQMAX=NSQ
              NGLMAX=NGL
           ENDIF   
   

            IF (PR.EQ.'gdcpl') THEN
               CALL LOCATE(MGLGR,NGLMAX,MGL,KGL)   
            ELSE   
               CALL LOCATE(MSQGR,NSQMAX,MSQ,KSQ)
               IF ((PR.NE.'st').AND.(PR.NE.'sdcpl')) 
     &           CALL LOCATE(MGLGR,NGLMAX,MGL,KGL)
            ENDIF
   

c... calculate the offset where one has to start reading in the grid
           IF (PR.EQ.'gdcpl') THEN
              KINGL=MIN(MAX(KGL-(NITP-1)/2,1),NGLMAX+1-NITP)
           ELSE   
              KINSQ=MIN(MAX(KSQ-(NITP-1)/2,1),NSQMAX+1-NITP)
              IF ((PR.NE.'st').AND.(PR.NE.'sdcpl'))
     &               KINGL=MIN(MAX(KGL-(NITP-1)/2,1),NGLMAX+1-NITP)
           ENDIF   
 
c... interpolate

c... for LO (with scale variation):
            CALL GRIDINTP(MSQGR,MGLGR,LLOGR,NITP,KINSQ,KINGL,
     &           MSQ,MGL,LO,DLO)
c$$$            CALL GRIDINTP(MSQGR,MGLGR,LLOGRPL,NITP,KINSQ,KINGL,
c$$$     &           MSQ,MGL,LOPL,DLOPL)
c$$$            CALL GRIDINTP(MSQGR,MGLGR,LLOGRMI,NITP,KINSQ,KINGL,
c$$$     &           MSQ,MGL,LOMI,DLOMI)
            
c... for NLO (with scale variation):
            CALL GRIDINTP(MSQGR,MGLGR,LNLOGR,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NLO,DNLO)
c$$$            CALL GRIDINTP(MSQGR,MGLGR,LNLOGRPL,NITP,KINSQ,KINGL,
c$$$     &           MSQ,MGL,NLOPL,DNLOPL)
c$$$            CALL GRIDINTP(MSQGR,MGLGR,LNLOGRMI,NITP,KINSQ,KINGL,
c$$$     &           MSQ,MGL,NLOMI,DNLOMI)
            
c... for NLL+NLO (with scale variation):
            CALL GRIDINTP(MSQGR,MGLGR,LNLLGR,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NLL,DNLL)
            CALL GRIDINTP(MSQGR,MGLGR,LNLLGRPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NLLPL,DNLLPL)
            CALL GRIDINTP(MSQGR,MGLGR,LNLLGRMI,NITP,KINSQ,KINGL,
     &           MSQ,MGL,NLLMI,DNLLMI)

c... for pdf and as errors:
            CALL GRIDINTP(MSQGR,MGLGR,PDFGRPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,PDFPL,DPDFPL)
            CALL GRIDINTP(MSQGR,MGLGR,PDFGRMI,NITP,KINSQ,KINGL,
     &           MSQ,MGL,PDFMI,DPDFMI)
            CALL GRIDINTP(MSQGR,MGLGR,ASGRPL,NITP,KINSQ,KINGL,
     &           MSQ,MGL,ASPL,DASPL)
            CALL GRIDINTP(MSQGR,MGLGR,ASGRMI,NITP,KINSQ,KINGL,
     &           MSQ,MGL,ASMI,DASMI)
            

c... write out results:

       IF (PR.EQ.'gdcpl') THEN
    
            WRITE(*,11) MGL,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)
            WRITE(23,12) MGL,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)


        ELSEIF ((PR.EQ.'st').OR.(PR.EQ.'sdcpl')) THEN
    
            WRITE(*,11) MSQ,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)
            WRITE(23,12) MSQ,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)

        ELSE

            WRITE(*,13) MSQ,MGL,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)
            WRITE(23,13) MSQ,MGL,DEXP(LO),DEXP(NLO),DEXP(NLL),
     &           DEXP(NLLPL),-DEXP(NLLMI),PDFPL,-PDFMI,ASPL,-ASMI,
     &           DEXP(NLO)/DEXP(LO),DEXP(NLL)/DEXP(NLO)


        ENDIF    

 
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)

 11   FORMAT (' ',1(F7.0,2x),3(E11.3),'  ',(E11.3),' ',(E11.3),(f8.1),
     &     '   ',(f8.1),'  ',(f8.1),'  ',(f8.1),'',2(f8.2))
 12   FORMAT (' ',1(F7.0,2x),3(E11.3),' ',(E11.3),' ',(E11.3),(f8.1),
     &     '   ',(f8.1),'  ',(f8.1),'  ',(f8.1),'',2(f8.2))
 13   FORMAT (' ',2(F7.0,2x),3(E11.3),'  ',(E11.3),' ',(E11.3),(f8.1),
     &     '   ',(f8.1),'  ',(f8.1),'  ',(f8.1),'',2(f8.2))
 14   FORMAT (' ',2(F7.0,2x),3(E11.3),' ',(E11.3),' ',(E11.3),(f8.1),
     &     '   ',(f8.1),'  ',(f8.1),'  ',(f8.1),'',2(f8.2))
      
      END

c.......................................................................

      SUBROUTINE INOUT(MODE)
      
      IMPLICIT NONE
      
      INTEGER MODE,NSQ,NGL,NSQHM,NGLHM
      
      REAL*8 MSQ,MGL,MSQLO,MSQUP,MGLLO,MGLUP
 
      CHARACTER*5 PR,MSQCHAR,MGLCHAR,PDF
      INTEGER NITP
      COMMON/INPUT/PR

      COMMON/MASSES/MSQ,MGL
      COMMON/GRIDSIZE/NSQ,NGL,NSQHM,NGLHM
      COMMON/LIMITS/MSQLO,MSQUP,MGLLO,MGLUP
      COMMON/INTERPLEVEL/NITP
      


      IF (MODE.EQ.0) THEN


         

c... read in command-line parameters
      
      CALL GETARG(1,PR)
      CALL GETARG(2,PDF)
      CALL GETARG(3,MSQCHAR)
      READ(MSQCHAR,*) MSQ      
      IF ((PR.NE.'st').AND.(PR.NE.'sdcpl').AND.(PR.NE.'gdcpl')) THEN
         CALL GETARG(4,MGLCHAR)
         READ(MGLCHAR,*) MGL
      ENDIF   
      IF (PR.EQ.'gdcpl') MGL=MSQ


c... for this version disable option tot
      IF (PR.EQ.'tot') THEN  
       WRITE(*,*) 'OPTION ',PR,' NOT AVAILABLE IN THE CURRENT VERSION'
       STOP 
      ENDIF   
         
c... define lower and upper limits on squark and gluino masses


      IF (PR.EQ.'st') THEN
         MSQLO=1.d2
         MSQUP=2.d3
         MGLLO=1.493d3
         MGLUP=1.493d3
      ELSEIF ((PR.EQ.'gg').OR.(PR.EQ.'sg')) THEN 
         MSQLO=2.d2
         MSQUP=4.5d3
         MGLLO=2.d2
         MGLUP=2.5d3  
       ELSE
         MSQLO=2.d2
         MSQUP=2.5d3
         MGLLO=2.d2
         MGLUP=2.5d3  
      ENDIF   



c... check the masses 


        IF ((PR.NE.'gdcpl').AND.((MSQ.LT.MSQLO).OR.(MSQ.GT.MSQUP))) THEN
            write(*,*) 'TOO LOW/HIGH SQUARK MASS'
            STOP
        ENDIF
         
        IF ((PR.NE.'st').AND.(PR.NE.'sdcpl')
     &            .AND.((MGL.LT.MGLLO).OR.(MGL.GT.MGLUP))) THEN
            write(*,*) 'TOO LOW/HIGH GLUINO MASS'
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
         ELSEIF (PR.EQ.'tot') THEN
            open(23,file="tot.out",status='unknown') 
         ELSEIF (PR.EQ.'sdcpl') THEN
            open(23,file="sdcpl.out",status='unknown') 
         ELSEIF (PR.EQ.'gdcpl') THEN
            open(23,file="gdcpl.out",status='unknown') 
         ELSE
            write(*,*) 'WRONG INPUT PROCESS'
            STOP
         ENDIF
         
c... write header for output 
c... stops and squarks in the decoupling limit:
       
      IF ((PR.EQ.'st').OR.(PR.EQ.'sdcpl')) THEN  
         write(*,*)
         IF (PDF.EQ.'mstw') THEN
          write(*,*) '# LHC @ 8 TeV, MSTW2008NLO (MSTW2008LO for LO)'
         ELSEIF (PDF.EQ.'cteq') THEN
          write(*,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO)' 
         ENDIF 
         write(*,*) '# process: ', PR
         write(*,*) '# ms[GeV]  LO[pb]     NLO[pb]    NLL+NLO[pb]  d_mu+[pb]
     & [pb] d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%]  K_NLO 
     & K_NLL'
         write(*,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------'

         write(23,*)
         IF (PDF.EQ.'mstw') THEN
          write(23,*) '# LHC @ 8 TeV, MSTW2008NLO (MSTW2008LO for LO)'
         ELSEIF (PDF.EQ.'cteq') THEN
          write(23,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO)' 
         ENDIF
         write(23,*) '# process: ', PR
         write(23,*) '# ms[GeV]  LO[pb]     NLO[pb]    NLL+NLO[pb]  d_mu+[pb]
     & [pb] d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%]  K_NLO 
     & K_NLL'
         write(23,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------'


c... gluinos in the decoupling limit:         

       ELSEIF  (PR.EQ.'gdcpl') THEN  
         write(*,*)
         IF (PDF.EQ.'mstw') THEN
          write(*,*) '# LHC @ 8 TeV, MSTW2008NLO (MSTW2008LO for LO)'
         ELSEIF (PDF.EQ.'cteq') THEN
          write(*,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO)' 
         ENDIF 
         write(*,*) '# process: ', PR
         write(*,*) '# mg[GeV]  LO[pb]     NLO[pb]    NLL+NLO[pb]  d_mu+[pb]
     & [pb] d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%]  K_NLO 
     & K_NLL'
         write(*,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------'

         write(23,*)
         IF (PDF.EQ.'mstw') THEN
          write(23,*) '# LHC @ 8 TeV, MSTW2008NLO(MSTW2008LO for LO) '
         ELSEIF (PDF.EQ.'cteq') THEN
          write(23,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO) ' 
         ENDIF
         write(23,*) '# process: ', PR
         write(23,*) '# mg[GeV]  LO[pb]     NLO[pb]    NLL+NLO[pb]  d_mu+[pb]
     & [pb] d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%]  K_NLO 
     & K_NLL'
         write(23,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------' 

c... other processes:
       ELSE
         write(*,*)
         IF (PDF.EQ.'mstw') THEN
          write(*,*) '# LHC @ 8 TeV, MSTW2008NLO (MSTW2008LO for LO)'
         ELSEIF (PDF.EQ.'cteq') THEN
          write(*,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO) ' 
         ENDIF
         write(*,*) '# process: ', PR
         write(*,*) '# ms[GeV]  mg[GeV]  LO[pb]     NLO[pb]    NLL+NLO[p
     &b]  d_mu+[pb]  d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%] 
     &  K_NLO   K_NLL'
         write(*,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------'

         write(23,*)
         IF (PDF.EQ.'mstw') THEN
          write(23,*) '# LHC @ 8 TeV, MSTW2008NLO (MSTW2008LO for LO)'
         ELSEIF (PDF.EQ.'cteq') THEN
          write(23,*) '# LHC @ 8 TeV, CTEQ6.6M (CTEQ6L1 for LO)' 
         ENDIF
         write(23,*) '# process: ', PR
         write(23,*) '# ms[GeV]  mg[GeV]  LO[pb]     NLO[pb]    NLL+NLO[p
     &b]  d_mu+[pb]  d_mu-[pb]   d_pdf+[%]  d_pdf-[%]  d_as+[%] d_as-[%] 
     &  K_NLO   K_NLL'
         write(23,*) '------------------------------------------------------
     &------------------------------------------------------------------
     &---------------'

       ENDIF   



      ELSE IF (MODE.EQ.1) THEN 

c... open the grid file
       IF (PDF.EQ.'mstw') THEN  
         IF (PR.EQ.'sb') THEN
            open(21,file="sb_nllnlo_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'ss') THEN
            open(21,file="ss_nllnlo_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'gg') THEN
            open(21,file="gg_nllnlo_mstw2008.grid",status='unknown')
            open(31,file="gg_nllnlo_hm_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'sg') THEN
            open(21,file="sg_nllnlo_mstw2008.grid",status='unknown') 
            open(31,file="sg_nllnlo_hm_mstw2008.grid",status='unknown') 
        ELSEIF (PR.EQ.'st') THEN
            open(21,file="st_nllnlo_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'tot') THEN
            open(21,file="tot_nllnlo_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'sdcpl') THEN
            open(21,file="sdcpl_nllnlo_mstw2008.grid",status='unknown') 
         ELSEIF (PR.EQ.'gdcpl') THEN
            open(21,file="gdcpl_nllnlo_mstw2008.grid",status='unknown') 
         ELSE   
            write(*,*) 'WRONG INPUT PROCESS'
            STOP
         ENDIF
        ELSEIF (PDF.EQ.'cteq') THEN
         IF (PR.EQ.'sb') THEN
            open(21,file="sb_nllnlo_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'ss') THEN
            open(21,file="ss_nllnlo_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'gg') THEN
            open(21,file="gg_nllnlo_cteq6.grid",status='unknown') 
            open(31,file="gg_nllnlo_hm_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'sg') THEN
            open(21,file="sg_nllnlo_cteq6.grid",status='unknown') 
            open(31,file="sg_nllnlo_hm_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'st') THEN
             open(21,file="st_nllnlo_cteq6.grid",status='unknown')
         ELSEIF (PR.EQ.'tot') THEN
            open(21,file="tot_nllnlo_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'sdcpl') THEN
            open(21,file="sdcpl_nllnlo_cteq6.grid",status='unknown') 
         ELSEIF (PR.EQ.'gdcpl') THEN
            open(21,file="gdcpl_nllnlo_cteq6.grid",status='unknown') 
         ELSE
            write(*,*) 'WRONG INPUT PROCESS'
            STOP
         ENDIF
        ENDIF 



c... define the size of the grid 

         IF (PR.EQ.'st') THEN
           NSQ=39
           NGL=1
         ELSEIF (PR.EQ.'sdcpl') THEN
           NSQ=24
           NGL=1
         ELSEIF (PR.EQ.'gdcpl') THEN             
           NSQ=1
           NGL=24
         ELSE  
           NSQ=24
           NGL=24
         ENDIF 

         IF ((PR.EQ.'gg').OR.(PR.EQ.'sg')) THEN
           NSQHM=20
           NGLHM=24
         ELSE
           NSQHM=0
           NGLHM=0
         ENDIF

  

      ELSE IF (MODE.EQ.2) THEN 
         

c... define level of interpolation
      
        IF((PR.EQ.'ss').OR.(PR.EQ.'sg')
     &       .AND.(MSQ.GE.MSQLO).AND.(MSQ.LT.5.d2)) THEN
         NITP=6
        ELSEIF ((PR.EQ.'gg').AND.(MGL.GE.MGLLO).AND.(MGL.LT.5.d2)) THEN
         NITP=6 
        ELSE 
         NITP=4
        ENDIF 

      end if
      
      RETURN
      END

c.......................................................................

      SUBROUTINE GRIDINTP(MSQGR,MGLGR,DATAGR,NITP,KINSQ,KINGL,
     &                MSQ,MGL,INTPOUT,DINTPOUT)

      IMPLICIT NONE
      INTEGER NITP,KINSQ,KINGL,I,J
      REAL*8 MSQGR,MGLGR,DATAGR,MSQINGR,MGLINGR,DATAINGR,
     &       MSQ,MGL,INTPOUT,DINTPOUT,DATAINGRST
      DIMENSION MSQGR(100),MGLGR(100),DATAGR(100,100),DATAINGRST(100)
      DIMENSION MSQINGR(NITP),MGLINGR(NITP),DATAINGR(NITP,NITP)

      CHARACTER*5 PR
      COMMON/INPUT/PR

c... copy the grid

      DO I=1,NITP
        IF (PR.EQ.'gdcpl') THEN
          MGLINGR(I)=MGLGR(KINGL+I-1)
          DATAINGRST(I)=DATAGR(1,KINGL+I-1)           
        ELSE 
          MSQINGR(I)=MSQGR(KINSQ+I-1) 
          IF ((PR.EQ.'st').OR.(PR.EQ.'sdcpl')) THEN
             DATAINGRST(I)=DATAGR(KINSQ+I-1,1)
          ELSE
             DO J=1,NITP
              MGLINGR(J)=MGLGR(KINGL+J-1)
              DATAINGR(I,J)=DATAGR(KINSQ+I-1,KINGL+J-1)
             ENDDO
          ENDIF
        ENDIF  
      ENDDO  

c... interpolate

      IF (PR.EQ.'gdcpl') THEN
         CALL POLINT(MGLINGR,DATAINGRST,NITP,MGL,INTPOUT,DINTPOUT)    
      ELSEIF ((PR.EQ.'st').OR.(PR.EQ.'sdcpl')) THEN
         CALL POLINT(MSQINGR,DATAINGRST,NITP,MSQ,INTPOUT,DINTPOUT)        
      ELSE
         CALL POLIN2(MSQINGR,MGLINGR,DATAINGR,
     &        NITP,NITP,MSQ,MGL,INTPOUT,DINTPOUT)
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
