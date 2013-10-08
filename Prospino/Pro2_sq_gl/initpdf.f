C**********************************************************************
C***  THE SET OF PARTON DENSITIES AND THE CORRESPONDING VALUE       ***
C***  FOR ALPHAS ARE DEFINED IN LO (INILO) AND NLO (ININLO)         ***
ctp   changed: always call the same routine as in new processes     ***
C**********************************************************************

      REAL*8 FUNCTION ALPS(SCALE)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      COMMON/CONST4/ISET

ctp iset=0,1 for LO,NLO, set in inilo and ininlo
ctp alphas=alps(scale) set in hadron**.f
ctp alphas used in hadron**.f and matrix**.f (as argument)
ctp for new version make sure ALSINI is called first
      if (iset.eq.0) then 
         alps = ALPHAS(scale,1)
      else if (iset.eq.1) then 
         alps = ALPHAS(scale,2)
      end if

      RETURN
      END

C**********************************************************************      
      SUBROUTINE INILO(ivegas)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      integer   ivegas(1:4)
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT
      COMMON/CONST4/ISET
      COMMON / INTINI / IINI
ctp
      common/I_LOFAST/     i_lofast_common
      
      if (i_lofast_common .eq. 1) then
         ivegas(1) = 500
         ivegas(2) = 4
         ivegas(3) = 2000
         ivegas(4) = 4
      else
         ivegas(1) = 5000
         ivegas(2) = 8
         ivegas(3) = 20000
         ivegas(4) = 4
      end if
      
      ISET = 0
      IINI = 0
      call INIT_ALPHAS_ARG(iset)
      call INIT_PDF(iset)

      RETURN
      END

C**********************************************************************
      SUBROUTINE ININLO(ivegas)
      IMPLICIT REAL*8 (A-H,M-Z)
      IMPLICIT INTEGER (I,J)
      integer   ivegas(1:4)
      COMMON/CONST2/SCALE,SCAFAC,ICOLL,ISCAPT
      COMMON/CONST4/ISET
      COMMON / INTINI / IINI

      ivegas(1) = 2000
      ivegas(2) = 8
      ivegas(3) = 10000
      ivegas(4) = 8
      
      ISET = 1
      IINI = 0
      call INIT_ALPHAS_ARG(iset)
      call INIT_PDF(iset)

      RETURN
      END

C**********************************************************************
ctp independent subroutine, can be called from f77 part of the code
ctp copied from Xinitialize.f90
ctp calls go to Xget_pdf.f and to Xalphas_smooth.f
      subroutine INIT_ALPHAS_ARG(inlo) 
      
      implicit none

      integer  inlo
      integer  nq
      real*8   lam_dummy,acc,lam_qcd,mcq,mbq,mtq

      call GET_LAMBDA_QCD(inlo,lam_dummy)

      lam_qcd = lam_dummy
      acc = 1.e-8
      mcq = 1.5              
      mbq = 5.0
      mtq = 175000.0
      nq  = 5 

      call ALSINI(acc,lam_qcd,mcq,mbq,mtq,nq)

      return
      end

C**********************************************************************

      SUBROUTINE PARTONDF(X ,SCALE, XPDF )
      REAL*8 X,SCALE, XPDF(-6:6)
      INTEGER ISET,ii
      COMMON/CONST4/ISET

ctp in new version always call subroutine GET_PDF
ctp note that iset and inlo are the same thing in old/new prospino code
      call GET_PDF(iset,x,scale,xpdf)

ctp originally this was the PDFlib call, i.e. returning x*PDF
      do ii=5,6
         xpdf(+ii) = 0.D0
         xpdf(-ii) = 0.D0
      end do

      do ii=-4,4
         xpdf(ii) = x * xpdf(ii)
      end do

      RETURN
      END


