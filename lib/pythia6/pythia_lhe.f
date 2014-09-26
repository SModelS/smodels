
C...Preamble: declarations.
      PROGRAM PYTHIA
C...All real arithmetic in double precision.
      IMPLICIT NONE

      CHARACTER*72 SIN
      COMMON/MYANLC/CHNAME, TOPDIR,FRZNAM
      CHARACTER*6 CHNAME
      CHARACTER*8 TOPDIR
      CHARACTER*10 CHNEV
      CHARACTER*60 FRZNAM
      CHARACTER*132 SLHANAME,LHENAME
      INTEGER I,NEV,ILIST,IMODE
      INTEGER ISLHA,ILHEU
      REAL*8 CMS
      REAL*8 WT4
      logical t,pythia2parton
      INTEGER I1,I2
      INTEGER IUNIT
      INTEGER FUNIT
      character*132 vfname
      logical slhaflag,lheflag
      real*8 pyp
      integer idel

      EXTERNAL PYDATA
      DO I=1,1000
       READ(*,'(A)',ERR=100) SIN
       IF(SIN.EQ.'end') goto 100
       CALL PYGIVE(SIN)
      ENDDO
100   READ(*,*) NEV,ILIST,IMODE,CMS
C-----------------------------------------------------------------
      IF(IMODE.EQ.1) THEN
       CALL PYINIT('CMS','p','pbar',CMS)
      ELSEIF(IMODE.EQ.0) THEN
       CALL PYINIT('CMS','p','p',CMS)
      ELSEIF(IMODE.EQ.2) THEN
       CALL PYINIT('CMS','e+','e-',CMS)
      ELSEIF(IMODE.EQ.3) THEN
       CALL PYINIT('CMS','gamma','gamma',CMS)
      ELSEIF(IMODE.EQ.4) THEN
       CALL PYINIT('CMS','e+','nu_e',CMS)
      ELSEIF(IMODE.LT.0) THEN
       CALL PYINIT('USER','','',0D0)
      ENDIF

C-----------------------------------------------------------------
      DO I=1,1000
       READ(*,'(A)',ERR=100) SIN
       IF(SIN.EQ.'end') goto 101
       CALL PYGIVE(SIN)
      ENDDO
101   CONTINUE          

      DO 200 I=1,NEV
       CALL PYUPEV
       IF(I.le.10) CALL PYLIST(7)
       IF(I.le.10) CALL PYLIST(1)
 200  CONTINUE
C/*
      IF(ILIST.NE.0) CALL PYLIST(12)
C-----------------------------------------------------------------
c      CALL USRDONE
C-----------------------------------------------------------------
      call pyupin
      CALL PYSTAT(1)
      CALL PYLHEF

      STOP
      END
