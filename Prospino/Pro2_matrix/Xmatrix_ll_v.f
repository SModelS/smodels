cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     LL_QBB(MASSIN,CS)                                                c
c                                                                      c
c     LL_QBV(MASSIN,CS,CV)                                             c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = t2                                                c
c       MASSIN(3)  = u2                                                c
c       MASSIN(4)  = t1                                                c
c       MASSIN(5)  = u1                                                c
c       MASSIN(6)  = mn1                                               c
c       MASSIN(7)  = mn2                                               c
c       MASSIN(8)  = mz                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(10) = mg                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = mu                                                c
c       MASSIN(13-20) not needed                                       c
c                                                                      c
c       CS(1:4)  TYPICAL S*S CHANNEL COUPLINGS                         c
c                                                                      c
c       N.B.: THE SYMMETRY FACTOR 1/2 NOT INCLUDED                     c
c                                                                      c
c                                                                      c
c    NEEDED MANDELSTAM VARIABLES :                                     c
c                                                                      c
c       Q(K1) + QB(K2) -> N1(P1) + N2(P2) [+G(K3)]                     c
c                                                                      c
c       S  = 2(K1.K2)                                                  c
c       S3 = 2(K3.P2)                                                  c
c       S4 = 2(K3.K1)                                                  c
c       S5 = 2(P1.P2) + M1^2 + M2^2                                    c
c       T1 = 2(K1.P1)                                                  c
c       U1 = 2(K2.P1)                                                  c
c       T2 = 2(K2.P2)                                                  c
c       U2 = 2(K1.P2)                                                  c
c       TP = 2(K2.K3)                                                  c
c       UP = 2(K1.K3)                                                  c
c                                                                      c
c                                                                      c
c  CHANGE 20.7.98: REMOVE ALL NON-S*S STRUCTURES FOR DRELL YAN         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function LL_QBB(massin,Csdum)

      implicit none 

      integer    n
      real*8     Csdum(4),massin(1:20),Cs(4),Nc
     &          ,s,m1,m2,t2,u2
      complex*16 QBB_s,QBB_tx,QBB_tm,QBB_ux,QBB_um 

      Nc = 3.D0

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      m1 = massin(6)
      m2 = massin(7)

      do n=1,4 
         Cs(n) = Csdum(n)
      end do

      QBB_s  = dcmplx(0.D0,0.D0)
      QBB_tm = dcmplx(0.D0,0.D0)
      QBB_tx = dcmplx(0.D0,0.D0)
      QBB_um = dcmplx(0.D0,0.D0)
      QBB_ux = dcmplx(0.D0,0.D0)

c               the form output 
      QBB_s =  32.D0 * Cs(1) * ( t2*u2 - s*m2**2 )/s**2

      LL_QBB = real( QBB_s + QBB_tx + QBB_tm + QBB_ux + QBB_um )/2.D0
c               the prefactors removed in the form program 
c                   except for alphas, which is cut 
      LL_QBB = LL_QBB * Nc 

c               the averaging factors
      LL_QBB = LL_QBB /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      LL_QBB = LL_QBB * (abs(m1)+abs(m2))**2/4.D0

      end


c --------------------------------------------------------------------
      subroutine BORN_BARTS_LL(massin,Cs,
     &                b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2)

      implicit none 

      real*8     massin(1:20),Cs(4),Nc
     &          ,s,m2,t2,u2
     &          ,b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2
      complex*16 QBB_ss,QBB_st,QBB_su
     &          ,QBB_txs,QBB_txt,QBB_tmu,QBB_tms
     &          ,QBB_uxs,QBB_uxu,QBB_umt,QBB_ums 

      Nc = 3.D0

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      m2 = massin(7)

c               the form output 
c               the form output 
      QBB_ss =  32.D0 * Cs(1) * ( t2*u2 - s*m2**2 )/s**2

      QBB_st = dcmplx(0.D0,0.D0)
      QBB_su = dcmplx(0.D0,0.D0)

      QBB_txs = dcmplx(0.D0,0.D0)
      QBB_txt = dcmplx(0.D0,0.D0)

      QBB_tmu = dcmplx(0.D0,0.D0)
      QBB_tms = dcmplx(0.D0,0.D0)

      QBB_uxs = dcmplx(0.D0,0.D0)
      QBB_uxu = dcmplx(0.D0,0.D0)

      QBB_umt = dcmplx(0.D0,0.D0)
      QBB_ums = dcmplx(0.D0,0.D0)

c               check if the entries are really not complex 
      if (abs(aimag(QBB_ss )).gt.1.D-12) print *,'QBB_ss : ',QBB_ss
      if (abs(aimag(QBB_st )).gt.1.D-12) print *,'QBB_st : ',QBB_st
      if (abs(aimag(QBB_su )).gt.1.D-12) print *,'QBB_su : ',QBB_su
      if (abs(aimag(QBB_txs)).gt.1.D-12) print *,'QBB_txs: ',QBB_txs
      if (abs(aimag(QBB_txt)).gt.1.D-12) print *,'QBB_txt: ',QBB_txt
      if (abs(aimag(QBB_tmu)).gt.1.D-12) print *,'QBB_tmu: ',QBB_tmu
      if (abs(aimag(QBB_tms)).gt.1.D-12) print *,'QBB_tms: ',QBB_tms
      if (abs(aimag(QBB_uxs)).gt.1.D-12) print *,'QBB_uxs: ',QBB_uxs
      if (abs(aimag(QBB_uxu)).gt.1.D-12) print *,'QBB_uxu: ',QBB_uxu
      if (abs(aimag(QBB_umt)).gt.1.D-12) print *,'QBB_umt: ',QBB_umt
      if (abs(aimag(QBB_ums)).gt.1.D-12) print *,'QBB_ums: ',QBB_ums
      
      b_s   = Nc*real( QBB_ss  + QBB_st   + QBB_su )
      b_tx  = Nc*real( QBB_txs + QBB_txt )
      b_tm  = Nc*real( QBB_tmu + QBB_tms )
      b_ux  = Nc*real( QBB_uxs + QBB_uxu )
      b_um  = Nc*real( QBB_umt + QBB_ums )
      b_tx4 = Nc*real( QBB_txs + QBB_txt )
      b_tm4 = Nc*real( QBB_tmu + QBB_tms )
      b_ux2 = Nc*real( QBB_uxs + QBB_uxu )
      b_um2 = Nc*real( QBB_umt + QBB_ums )
 
      end


c --------------------------------------------------------------------
      real*8 function LL_QBV(massin,Csdum,Cv)

      implicit none 

      integer    n
      real*8     Csdum(4),massin(1:20),Pi,Nc,Cf,finite
     &          ,s,t,u,m1,m2,mg,ms,mav2
     &          ,t1,u1,t2,u2,ts,us,tg,ug,kaellen,Cs(4)
     &          ,SCB(1:7,1:5),SCBP(1),SCC(1:6,1:4),SCD(1:2,1:2)
     &          ,b_s,b_t,b_u,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2
      complex*16 Cv(4),QBV
     &          ,Cupu,Cupt,Clou,Clot

      Pi = 4.D0*atan(1.D0)
      Nc = 3.D0
      Cf = 4.D0/3.D0 

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
      m1 = massin(6)
      m2 = massin(7)
      mg = massin(10)
      ms = massin(11)

      Clot = Cv(1) 
      Cupt = Cv(2) 
      Cupu = Cv(3) 
      Clou = Cv(4) 

      t   = t2 + m2**2 
      u   = u2 + m2**2 
      ts  = t2 + m2**2 - ms**2 
      us  = u2 + m2**2 - ms**2 
      tg  = t2 + m2**2 - mg**2 
      ug  = u2 + m2**2 - mg**2 

      kaellen = s**2 + m1**4 + m2**4
     &         - 2*( s*m1**2 + s*m2**2 + m1**2*m2**2 )

      mav2 = (abs(m1)+abs(m2))**2/4.D0

c               finite part from virtual subtraction matrix element 
      finite = - 3.D0/2.D0 * Pi**2/6.D0 
     &         - 3.D0/2.D0 * log(s/mav2) 
     &         + log(s/mav2)**2 /2.D0

c               the scalar functions 
      call SCALAR_ARRAY_NN(massin,SCB,SCBP,SCC,SCD)

      do n=1,4 
         Cs(n) = Csdum(n)
      end do

c               the born type structures 
      call BORN_BARTS_LL(massin,Cs,
     &                b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2)
      b_t  = b_tx  + b_tm
      b_u  = b_ux  + b_um

c               insert the form output 
      QBV =
     +  + finite*b_s * ( 4 )
      QBV = QBV + finite*b_t * ( 4 )
      QBV = QBV + finite*b_u * ( 4 )
      QBV = QBV + b_s * (  - 2 )
      QBV = QBV + b_t * (  - 2 )
      QBV = QBV + b_u * (  - 2 )
      QBV = QBV + SCB(1,1)*Cupt*b_t * (  - 4*m2*mg*t2**(-1) )
      QBV = QBV + SCB(1,1)*Clot*b_t * (  - 4*m1*mg*t1**(-1) )
      QBV = QBV + SCB(1,1)*b_t * (  - 4*m1**2*ts**(-1) - 4*m2**2*
     +    ts**(-1) + 4*mg**2*ts**(-1) + 4*s*ts**(-1) + 4*u*ts**(-1) )
      QBV = QBV + SCB(1,1)*b_um2 * (  - 4*s**(-1)*us )
      QBV = QBV + SCB(1,1)*b_ux2 * ( 4*m1**2*s*t1**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m1**2*s*u1**(-1)*t2**(-1)*u2**(-1)*us + 2*
     +    m1**2*u1**(-1)*u2**(-1)*us + 2*m1**2*t2**(-1)*u2**(-1)*us + 2
     +    *m2**2*s*u1**(-1)*t2**(-1)*u2**(-1)*us + 2*m2**2*u1**(-1)*
     +    u2**(-1)*us - 2*m2**2*t2**(-1)*u2**(-1)*us - 4*u2**(-1)*us )
      QBV = QBV + SCB(1,2)*Cupu*b_u * (  - 4*m1*mg*u1**(-1) )
      QBV = QBV + SCB(1,2)*Clou*b_u * (  - 4*m2*mg*u2**(-1) )
      QBV = QBV + SCB(1,2)*b_u * ( 4*mg**2*us**(-1) - 4*u*us**(-1)
     +     )
      QBV = QBV + SCB(1,2)*b_tm4 * (  - 4*s**(-1)*ts )
      QBV = QBV + SCB(1,2)*b_tx4 * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     + 4*m2**2*s*t1**(-1)*u1**(-1)*t2**(-1)*ts + 4*m2**2*s*
     +    t1**(-1)*t2**(-1)*u2**(-1)*ts + 2*m2**2*t1**(-1)*t2**(-1)*ts
     +     + 2*s*t1**(-1)*u1**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts - 2*s*
     +    u1**(-1)*t2**(-1)*ts - 2*s**2*t1**(-1)*u1**(-1)*t2**(-1)*ts
     +     - 2*t1**(-1)*ts + 2*u1**(-1)*ts )
      QBV = QBV + SCB(1,3)*b_t * (  - 4*mg**2*ts**(-1) + 4*ms**2*
     +    ts**(-1) )
      QBV = QBV + SCB(1,3)*b_u * (  - 4*mg**2*us**(-1) + 4*ms**2*
     +    us**(-1) )
      QBV = QBV + SCB(3,1)*b_t * (  - 4 - 4*m1**2*t1**(-1) + 4*
     +    m1**2*ts**(-1) - 4*m2**2*t2**(-1) + 4*m2**2*ts**(-1) + 4*
     +    ms**2*ts**(-1) - 4*s*ts**(-1) - 4*u*ts**(-1) )
      QBV = QBV + SCB(3,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCB(3,1)*b_tx * ( 3*m1**2*t1**(-1)*t2**(-1)*ts
     +     + m2**2*t1**(-1)*t2**(-1)*ts + t1**(-1)*ts + 3*t2**(-1)*ts )
      QBV = QBV + SCB(3,2)*b_u * (  - 4 - 4*m1**2*u1**(-1) - 4*
     +    m2**2*u2**(-1) + 4*ms**2*us**(-1) + 4*u*us**(-1) )
      QBV = QBV + SCB(3,2)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCB(3,2)*b_ux * ( 3*m1**2*u1**(-1)*u2**(-1)*us
     +     + m2**2*u1**(-1)*u2**(-1)*us + u1**(-1)*us + 3*u2**(-1)*us )
      QBV = QBV + SCB(3,3)*b_t * (  - 8*ms**2*ts**(-1) )
      QBV = QBV + SCB(3,3)*b_u * (  - 8*ms**2*us**(-1) )
      QBV = QBV + SCB(3,4)*Cupu*b_u * ( 4*m1*mg*u1**(-1) )
      QBV = QBV + SCB(3,4)*Clot*b_t * ( 4*m1*mg*t1**(-1) )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tm * ( 4*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*s**(-1)*u*ts + 4*m1**2*ts + 4*m2**2*
     +    s**(-1)*u*ts + 8*m2**2*ts - 4*m2**4*s**(-1)*ts - 4*s*ts - 4*u
     +    *ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tx * ( 6*m1**2*m2**2*s*
     +    t1**(-1)*t2**(-1)*ts + 4*m1**2*m2**2*t2**(-1)*ts - 4*m1**2*
     +    m2**4*t1**(-1)*t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts + 2*m1**2*
     +    s*t2**(-1)*ts + 8*m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 2*m1**4*
     +    s*t1**(-1)*t2**(-1)*ts - 4*m1**4*t2**(-1)*ts - 4*m1**6*
     +    t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_um * (  - 4*m1**2*m2**2
     +    *s**(-1)*us + 4*m1**2*s**(-1)*u*us - 4*m1**2*us - 4*m2**2*
     +    s**(-1)*u*us - 8*m2**2*us + 4*m2**4*s**(-1)*us + 4*s*us + 4*u
     +    *us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_ux * ( 10*m1**2*m2**2*s
     +    *u1**(-1)*u2**(-1)*us + 2*m1**2*m2**2*u1**(-1)*us + 2*m1**2*
     +    m2**2*u2**(-1)*us - 3*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 4*
     +    m1**2*s*u1**(-1)*us + 3*m1**4*m2**2*u1**(-1)*u2**(-1)*us - 3*
     +    m1**4*u1**(-1)*us - m1**4*u2**(-1)*us - m1**6*u1**(-1)*
     +    u2**(-1)*us - 2*m2**2*s*u1**(-1)*us + 2*m2**2*s*u2**(-1)*us
     +     - 2*m2**4*s*u1**(-1)*u2**(-1)*us + m2**4*u1**(-1)*us - m2**4
     +    *u2**(-1)*us + m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_um2 * ( 4*m1**2*m2**2*
     +    s**(-1)*us - 4*m1**2*s**(-1)*u*us + 4*m1**2*us + 4*m2**2*
     +    s**(-1)*u*us + 8*m2**2*us - 4*m2**4*s**(-1)*us - 4*s*us - 4*u
     +    *us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_ux2 * (  - 4*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 4*m1**2*m2**2*u2**(-1)*us + 6*
     +    m1**2*m2**4*u1**(-1)*u2**(-1)*us + 4*m1**2*s*u1**(-1)*us - 8*
     +    m1**2*s*u2**(-1)*us - 6*m1**4*m2**2*u1**(-1)*u2**(-1)*us - 8*
     +    m1**4*s*u1**(-1)*u2**(-1)*us + 2*m1**4*u1**(-1)*us + 2*m1**4*
     +    u2**(-1)*us + 2*m1**6*u1**(-1)*u2**(-1)*us + 4*m2**2*s*
     +    u1**(-1)*us - 4*m2**2*s*u2**(-1)*us + 4*m2**4*s*u1**(-1)*
     +    u2**(-1)*us - 2*m2**4*u1**(-1)*us + 2*m2**4*u2**(-1)*us - 2*
     +    m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tm4 * (  - 4*m1**2*
     +    m2**2*s**(-1)*ts + 4*m1**2*s**(-1)*u*ts - 4*m1**2*ts - 4*
     +    m2**2*s**(-1)*u*ts - 8*m2**2*ts + 4*m2**4*s**(-1)*ts + 4*s*ts
     +     + 4*u*ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tx4 * (  - 8*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**2*t2**(-1)*ts + 4*
     +    m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts - 6*
     +    m1**2*s*t2**(-1)*ts - 8*m1**4*m2**2*t1**(-1)*t2**(-1)*ts - 6*
     +    m1**4*s*t1**(-1)*t2**(-1)*ts + 4*m1**4*t2**(-1)*ts + 4*m1**6*
     +    t1**(-1)*t2**(-1)*ts + 6*m2**2*s*t1**(-1)*ts - 6*m2**2*s*
     +    t2**(-1)*ts + 6*m2**4*s*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_t * ( 4*m1**2*t1**(-1) )
      QBV = QBV + SCB(3,4)*b_u * ( 4*m1**2*u1**(-1) )
      QBV = QBV + SCB(3,4)*b_tm * ( 4*s**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_tx * ( m1**2*t1**(-1)*t2**(-1)*ts - 
     +    m2**2*t1**(-1)*t2**(-1)*ts - t1**(-1)*ts + t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_ux * ( m1**2*u1**(-1)*u2**(-1)*us - 
     +    m2**2*u1**(-1)*u2**(-1)*us - u1**(-1)*us + u2**(-1)*us )
      QBV = QBV + SCB(3,4)*b_um2 * ( 4*s**(-1)*us )
      QBV = QBV + SCB(3,4)*b_ux2 * (  - 4*m1**2*s*t1**(-1)*
     +    u1**(-1)*u2**(-1)*us - 4*m1**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*b_tx4 * ( 2*m1**2*t1**(-1)*u1**(-1)*ts
     +     - 4*m1**2*t1**(-1)*t2**(-1)*ts - 4*m2**2*s*t1**(-1)*u1**(-1)
     +    *t2**(-1)*ts - 2*m2**2*t1**(-1)*u1**(-1)*ts + 2*s*t1**(-1)*
     +    t2**(-1)*ts + 2*s*u1**(-1)*t2**(-1)*ts + 2*s**2*t1**(-1)*
     +    u1**(-1)*t2**(-1)*ts + 2*t1**(-1)*ts - 2*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*Cupt*b_t * ( 4*m2*mg*t2**(-1) )
      QBV = QBV + SCB(3,5)*Clou*b_u * ( 4*m2*mg*u2**(-1) )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tm * (  - 4*m1**2*m2**2
     +    *s**(-1)*ts + 4*m1**2*s**(-1)*u*ts - 4*m2**2*s**(-1)*u*ts - 4
     +    *m2**2*ts + 4*m2**4*s**(-1)*ts - 4*u*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tx * ( 6*m1**2*m2**2*s*
     +    t1**(-1)*t2**(-1)*ts + 6*m1**2*m2**2*t1**(-1)*ts - 2*m1**2*
     +    m2**2*t2**(-1)*ts + 9*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*s*t1**(-1)*ts - 2*m1**2*s*t2**(-1)*ts - 9*m1**4*m2**2*
     +    t1**(-1)*t2**(-1)*ts - 2*m1**4*s*t1**(-1)*t2**(-1)*ts - 3*
     +    m1**4*t1**(-1)*ts + 3*m1**4*t2**(-1)*ts + 3*m1**6*t1**(-1)*
     +    t2**(-1)*ts + 4*m2**2*s*t1**(-1)*ts + 4*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts - 3*m2**4*t1**(-1)*ts - m2**4*t2**(-1)*ts - 3*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_um * ( 4*m1**2*m2**2*
     +    s**(-1)*us - 4*m1**2*s**(-1)*u*us + 4*m2**2*s**(-1)*u*us + 4*
     +    m2**2*us - 4*m2**4*s**(-1)*us + 4*u*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_ux * ( 2*m1**2*m2**2*s*
     +    u1**(-1)*u2**(-1)*us + 4*m1**2*m2**2*u1**(-1)*us + 8*m1**2*
     +    m2**4*u1**(-1)*u2**(-1)*us - 4*m1**4*m2**2*u1**(-1)*u2**(-1)*
     +    us + 6*m2**2*s*u1**(-1)*us - 2*m2**2*s*u2**(-1)*us + 6*m2**4*
     +    s*u1**(-1)*u2**(-1)*us - 4*m2**4*u1**(-1)*us - 4*m2**6*
     +    u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_um2 * (  - 4*m1**2*
     +    m2**2*s**(-1)*us + 4*m1**2*s**(-1)*u*us - 4*m2**2*s**(-1)*u*
     +    us - 4*m2**2*us + 4*m2**4*s**(-1)*us - 4*u*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_ux2 * (  - 8*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 4*m1**2*m2**2*u1**(-1)*us - 8*
     +    m1**2*m2**4*u1**(-1)*u2**(-1)*us - 6*m1**2*s*u1**(-1)*us + 6*
     +    m1**2*s*u2**(-1)*us + 4*m1**4*m2**2*u1**(-1)*u2**(-1)*us + 6*
     +    m1**4*s*u1**(-1)*u2**(-1)*us - 6*m2**2*s*u1**(-1)*us + 2*
     +    m2**2*s*u2**(-1)*us - 6*m2**4*s*u1**(-1)*u2**(-1)*us + 4*
     +    m2**4*u1**(-1)*us + 4*m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tm4 * ( 4*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*s**(-1)*u*ts + 4*m2**2*s**(-1)*u*ts + 4*
     +    m2**2*ts - 4*m2**4*s**(-1)*ts + 4*u*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tx4 * (  - 4*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts - 6*
     +    m1**2*m2**4*t1**(-1)*t2**(-1)*ts - 4*m1**2*s*t1**(-1)*ts + 4*
     +    m1**2*s*t2**(-1)*ts + 6*m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 4*
     +    m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*t1**(-1)*ts - 2*m1**4*
     +    t2**(-1)*ts - 2*m1**6*t1**(-1)*t2**(-1)*ts - 8*m2**2*s*
     +    t1**(-1)*ts + 4*m2**2*s*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**4*t1**(-1)*ts + 2*m2**4*t2**(-1)*ts + 2*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_t * ( 4*m2**2*t2**(-1) )
      QBV = QBV + SCB(3,5)*b_u * ( 4*m2**2*u2**(-1) )
      QBV = QBV + SCB(3,5)*b_tx * (  - 3*m1**2*t1**(-1)*t2**(-1)*
     +    ts + 3*m2**2*t1**(-1)*t2**(-1)*ts + 3*t1**(-1)*ts - 3*
     +    t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_um * ( 4*s**(-1)*us )
      QBV = QBV + SCB(3,5)*b_ux * (  - 3*m1**2*u1**(-1)*u2**(-1)*
     +    us + 3*m2**2*u1**(-1)*u2**(-1)*us + 3*u1**(-1)*us - 3*
     +    u2**(-1)*us )
      QBV = QBV + SCB(3,5)*b_ux2 * (  - 2*m1**2*s*u1**(-1)*
     +    t2**(-1)*u2**(-1)*us - 2*m1**2*t2**(-1)*u2**(-1)*us - 2*m2**2
     +    *s*u1**(-1)*t2**(-1)*u2**(-1)*us - 4*m2**2*u1**(-1)*u2**(-1)*
     +    us + 2*m2**2*t2**(-1)*u2**(-1)*us - 2*u1**(-1)*us + 2*
     +    u2**(-1)*us )
      QBV = QBV + SCB(3,5)*b_tm4 * ( 4*s**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_tx4 * (  - 4*m2**2*s*t1**(-1)*
     +    t2**(-1)*u2**(-1)*ts - 4*m2**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_um2 * (  - 4*m1**2*us
     +     - 4*m2**2*us + 4*s*us + 8*u*us )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_ux2 * ( 12*m1**2*m2**2*
     +    s*u1**(-1)*u2**(-1)*us + 4*m1**2*m2**2*u1**(-1)*us + 4*m1**2*
     +    m2**2*u2**(-1)*us + 2*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 2*
     +    m1**2*s*u1**(-1)*us + 2*m1**2*s*u2**(-1)*us + 2*m1**4*m2**2*
     +    u1**(-1)*u2**(-1)*us + 2*m1**4*s*u1**(-1)*u2**(-1)*us - 2*
     +    m1**4*u1**(-1)*us - 2*m1**4*u2**(-1)*us - 2*m1**6*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*s*u1**(-1)*us + 2*m2**2*s*u2**(-1)*us
     +     + 2*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**4*u1**(-1)*us - 2*
     +    m2**4*u2**(-1)*us - 2*m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_tm4 * ( 4*m1**2*ts + 4*
     +    m2**2*ts - 4*s*ts - 8*u*ts )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_tx4 * ( 12*m1**2*m2**2*
     +    s*t1**(-1)*t2**(-1)*ts + 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*
     +    m2**2*t2**(-1)*ts + 2*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*s*t1**(-1)*ts + 2*m1**2*s*t2**(-1)*ts + 2*m1**4*m2**2*
     +    t1**(-1)*t2**(-1)*ts + 2*m1**4*s*t1**(-1)*t2**(-1)*ts - 2*
     +    m1**4*t1**(-1)*ts - 2*m1**4*t2**(-1)*ts - 2*m1**6*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*s*t1**(-1)*ts + 2*m2**2*s*t2**(-1)*ts
     +     + 2*m2**4*s*t1**(-1)*t2**(-1)*ts - 2*m2**4*t1**(-1)*ts - 2*
     +    m2**4*t2**(-1)*ts - 2*m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(4,1)*b_s * ( 2 + 4*mg**2*s**(-1) - 4*ms**2*
     +    s**(-1) )
      QBV = QBV + SCB(4,1)*b_ux2 * ( 2*m1**2*u1**(-1)*u2**(-1)*us
     +     + 2*m2**2*u1**(-1)*u2**(-1)*us + 2*u1**(-1)*us + 2*u2**(-1)*
     +    us )
      QBV = QBV + SCB(4,1)*b_tx4 * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     + 2*m2**2*t1**(-1)*t2**(-1)*ts + 2*t1**(-1)*ts + 2*t2**(-1)*
     +    ts )
      QBV = QBV + SCB(6,1)*b_s * (  - 2 - 4*mg**2*s**(-1) + 4*
     +    ms**2*s**(-1) )
      QBV = QBV + SCB(6,1)*b_t * (  - 2 )
      QBV = QBV + SCB(6,1)*b_u * (  - 2 )
      QBV = QBV + SCB(6,1)*b_tx4 * (  - 2*m1**2*t1**(-1)*u1**(-1)*
     +    ts + 2*m2**2*t1**(-1)*u1**(-1)*ts - 2*s*t1**(-1)*u1**(-1)*ts
     +     - 2*t1**(-1)*ts - 2*u1**(-1)*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tm * (  - 4*m1**2*ts - 
     +    4*m2**2*ts + 4*s*ts + 8*u*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tx * (  - 12*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 6*m1**2*m2**2*t1**(-1)*ts - 2*
     +    m1**2*m2**2*t2**(-1)*ts - 5*m1**2*m2**4*t1**(-1)*t2**(-1)*ts
     +     - 4*m1**2*s*t1**(-1)*ts + m1**4*m2**2*t1**(-1)*t2**(-1)*ts
     +     + 3*m1**4*t1**(-1)*ts + m1**4*t2**(-1)*ts + m1**6*t1**(-1)*
     +    t2**(-1)*ts - 4*m2**2*s*t1**(-1)*ts - 4*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts + 3*m2**4*t1**(-1)*ts + m2**4*t2**(-1)*ts + 3*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_um * ( 4*m1**2*us + 4*
     +    m2**2*us - 4*s*us - 8*u*us )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_ux * (  - 12*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 6*m1**2*m2**2*u1**(-1)*us - 2*
     +    m1**2*m2**2*u2**(-1)*us - 5*m1**2*m2**4*u1**(-1)*u2**(-1)*us
     +     - 4*m1**2*s*u1**(-1)*us + m1**4*m2**2*u1**(-1)*u2**(-1)*us
     +     + 3*m1**4*u1**(-1)*us + m1**4*u2**(-1)*us + m1**6*u1**(-1)*
     +    u2**(-1)*us - 4*m2**2*s*u1**(-1)*us - 4*m2**4*s*u1**(-1)*
     +    u2**(-1)*us + 3*m2**4*u1**(-1)*us + m2**4*u2**(-1)*us + 3*
     +    m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(7,1)*b_s * (  - 6 )
      QBV = QBV + SCB(7,1)*b_tx * (  - m1**2*t1**(-1)*t2**(-1)*ts
     +     - 3*m2**2*t1**(-1)*t2**(-1)*ts - 3*t1**(-1)*ts - t2**(-1)*ts
     +     )
      QBV = QBV + SCB(7,1)*b_ux * (  - m1**2*u1**(-1)*u2**(-1)*us
     +     - 3*m2**2*u1**(-1)*u2**(-1)*us - 3*u1**(-1)*us - u2**(-1)*us
     +     )
      QBV = QBV + SCBP(1)*b_s * (  - 2*mg**2 + 2*ms**2 )
      QBV = QBV + SCBP(1)*b_t * (  - 2*mg**2 + 2*ms**2 )
      QBV = QBV + SCBP(1)*b_u * (  - 2*mg**2 + 2*ms**2 )
      QBV = QBV + SCC(1,1)*Clot*b_t * (  - 4*m1*mg*ms**2*s**(-1)*
     +    t1**(-1) - 4*m1*mg*s**(-1) + 4*m1*mg**3*s**(-1)*t1**(-1) )
      QBV = QBV + SCC(1,1)*b_um2 * ( 2*m1**2*s**(-2)*us + 4*mg**2*
     +    s**(-2)*us - 4*ms**2*s**(-2)*us - 2*s**(-2)*t*us )
      QBV = QBV + SCC(1,1)*b_ux2 * (  - 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 4*m1**2*mg**2*t1**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m1**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 4*
     +    m1**2*ms**2*t1**(-1)*u1**(-1)*u2**(-1)*us - 2*m2**2*mg**2*
     +    s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*s**(-1)*u1**(-1)
     +    *u2**(-1)*us + 2*mg**2*s**(-1)*u2**(-1)*us + 2*mg**2*u1**(-1)
     +    *u2**(-1)*us - 2*ms**2*s**(-1)*u1**(-1)*us - 2*ms**2*s**(-1)*
     +    u2**(-1)*us - 4*ms**2*u1**(-1)*u2**(-1)*us + 2*s*u1**(-1)*
     +    u2**(-1)*us + 2*u1**(-1)*us )
      QBV = QBV + SCC(1,2)*Cupt*b_t * (  - 4*m2*mg*ms**2*s**(-1)*
     +    t2**(-1) - 4*m2*mg*s**(-1) + 4*m2*mg**3*s**(-1)*t2**(-1) )
      QBV = QBV + SCC(1,2)*b_um2 * ( 2*m2**2*s**(-2)*us + 4*mg**2*
     +    s**(-2)*us - 4*ms**2*s**(-2)*us - 2*s**(-2)*t*us )
      QBV = QBV + SCC(1,2)*b_ux2 * (  - 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 2*m1**2*mg**2*s**(-1)*t2**(-1)*
     +    u2**(-1)*us - 2*m1**2*mg**2*u1**(-1)*t2**(-1)*u2**(-1)*us + 2
     +    *m1**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m1**2*ms**2*
     +    s**(-1)*t2**(-1)*u2**(-1)*us + 2*m1**2*ms**2*u1**(-1)*
     +    t2**(-1)*u2**(-1)*us + 2*m1**2*s**(-1)*u2**(-1)*us - 2*m2**2*
     +    mg**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*mg**2*s**(-1)*
     +    t2**(-1)*u2**(-1)*us - 2*m2**2*mg**2*u1**(-1)*t2**(-1)*
     +    u2**(-1)*us + 2*m2**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*
     +    m2**2*ms**2*s**(-1)*t2**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     +    u1**(-1)*t2**(-1)*u2**(-1)*us + 2*mg**2*s**(-1)*u2**(-1)*us
     +     + 2*mg**2*u1**(-1)*u2**(-1)*us - 4*ms**2*s**(-1)*u2**(-1)*us
     +     - 4*ms**2*u1**(-1)*u2**(-1)*us - 2*s**(-1)*t*u2**(-1)*us - 2
     +    *s**(-1)*us + 2*s*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(1,3)*Cupu*b_u * (  - 4*m1*mg*ms**2*s**(-1)*
     +    u1**(-1) - 4*m1*mg*s**(-1) + 4*m1*mg**3*s**(-1)*u1**(-1) )
      QBV = QBV + SCC(1,3)*b_tm4 * ( 2*m1**2*s**(-2)*ts + 4*mg**2*
     +    s**(-2)*ts - 4*ms**2*s**(-2)*ts - 2*s**(-2)*u*ts )
      QBV = QBV + SCC(1,3)*b_tx4 * (  - 2*m1**2*m2**2*s**(-1)*
     +    t1**(-1)*u1**(-1)*ts + 2*m1**2*mg**2*s**(-1)*t1**(-1)*
     +    u1**(-1)*ts - 2*m1**2*mg**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 2*m1**2*s**(-1)*
     +    t1**(-1)*ts + 2*m1**2*s**(-1)*u1**(-1)*ts + 2*m1**2*t1**(-1)*
     +    u1**(-1)*ts + 2*m1**4*s**(-1)*t1**(-1)*u1**(-1)*ts - 2*m2**2*
     +    mg**2*s**(-1)*t1**(-1)*u1**(-1)*ts - 2*m2**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts - 4*m2**2*mg**2*t1**(-1)*u1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 4*
     +    m2**2*ms**2*t1**(-1)*u1**(-1)*t2**(-1)*ts + 2*m2**2*s**(-1)*
     +    t1**(-1)*ts + 2*mg**2*s**(-1)*t1**(-1)*ts + 2*mg**2*s*
     +    t1**(-1)*u1**(-1)*t2**(-1)*ts + 4*mg**2*t1**(-1)*t2**(-1)*ts
     +     + 2*mg**2*u1**(-1)*t2**(-1)*ts - 2*ms**2*s**(-1)*t1**(-1)*ts
     +     + 2*ms**2*s**(-1)*u1**(-1)*ts - 2*ms**2*s*t1**(-1)*u1**(-1)*
     +    t2**(-1)*ts + 2*ms**2*t1**(-1)*u1**(-1)*ts - 6*ms**2*t1**(-1)
     +    *t2**(-1)*ts )
      QBV = QBV + SCC(1,3)*b_tx4 * (  - 2*ms**2*u1**(-1)*t2**(-1)*
     +    ts - 2*s**(-1)*u*t1**(-1)*ts - 2*s**(-1)*ts + 2*s*t1**(-1)*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(1,4)*Clou*b_u * (  - 4*m2*mg*ms**2*s**(-1)*
     +    u2**(-1) - 4*m2*mg*s**(-1) + 4*m2*mg**3*s**(-1)*u2**(-1) )
      QBV = QBV + SCC(1,4)*b_tm4 * ( 2*m2**2*s**(-2)*ts + 4*mg**2*
     +    s**(-2)*ts - 4*ms**2*s**(-2)*ts - 2*s**(-2)*u*ts )
      QBV = QBV + SCC(1,4)*b_tx4 * (  - 2*m1**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*s**(-1)*t1**(-1)*
     +    t2**(-1)*ts - 2*m2**2*mg**2*s**(-1)*t1**(-1)*t2**(-1)*ts - 4*
     +    m2**2*mg**2*t1**(-1)*t2**(-1)*u2**(-1)*ts + 2*m2**2*ms**2*
     +    s**(-1)*t1**(-1)*t2**(-1)*ts + 4*m2**2*ms**2*t1**(-1)*
     +    t2**(-1)*u2**(-1)*ts + 2*mg**2*s**(-1)*t1**(-1)*ts + 2*mg**2*
     +    t1**(-1)*t2**(-1)*ts - 2*ms**2*s**(-1)*t1**(-1)*ts - 2*ms**2*
     +    s**(-1)*t2**(-1)*ts - 4*ms**2*t1**(-1)*t2**(-1)*ts + 2*s*
     +    t1**(-1)*t2**(-1)*ts + 2*t2**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_t * (  - 4*m1**2*t1**(-1) + 4*ms**2*
     +    t1**(-1) )
      QBV = QBV + SCC(2,1)*b_tm * ( 4*t1**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_tx * ( m1**2*t1**(-1)*t2**(-1)*ts + 
     +    m2**2*t1**(-1)*t2**(-1)*ts - 2*ms**2*t1**(-1)*t2**(-1)*ts + 3
     +    *t1**(-1)*ts - t2**(-1)*ts )
      QBV = QBV + SCC(2,2)*b_t * (  - 4*m2**2*t2**(-1) + 4*ms**2*
     +    t2**(-1) )
      QBV = QBV + SCC(2,2)*b_tm * ( 4*t2**(-1)*ts )
      QBV = QBV + SCC(2,2)*b_tx * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     - 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts + 4*t2**(-1)*
     +    ts )
      QBV = QBV + SCC(2,3)*b_u * (  - 4*m1**2*u1**(-1) + 4*ms**2*
     +    u1**(-1) )
      QBV = QBV + SCC(2,3)*b_um * ( 4*u1**(-1)*us )
      QBV = QBV + SCC(2,3)*b_ux * (  - m1**2*u1**(-1)*u2**(-1)*us
     +     + 3*m2**2*u1**(-1)*u2**(-1)*us - 2*ms**2*u1**(-1)*u2**(-1)*
     +    us + 5*u1**(-1)*us - 3*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_u * (  - 4*m2**2*u2**(-1) + 4*ms**2*
     +    u2**(-1) )
      QBV = QBV + SCC(2,4)*b_um * ( 4*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_ux * ( 2*m2**2*u1**(-1)*u2**(-1)*us
     +     - 2*ms**2*u1**(-1)*u2**(-1)*us + 2*u2**(-1)*us )
      QBV = QBV + SCC(3,1)*b_s * (  - 8*mg**2*ms**2*s**(-2) + 4*
     +    mg**2*s**(-1) + 4*mg**4*s**(-2) + 4*ms**4*s**(-2) )
      QBV = QBV + SCC(3,1)*b_um2 * (  - 2*s**(-1)*us )
      QBV = QBV + SCC(3,1)*b_ux2 * (  - 2*mg**2*u1**(-1)*u2**(-1)*
     +    us + 4*ms**2*u1**(-1)*u2**(-1)*us - 2*s*u1**(-1)*u2**(-1)*us
     +     )
      QBV = QBV + SCC(3,1)*b_tm4 * (  - 2*s**(-1)*ts )
      QBV = QBV + SCC(3,1)*b_tx4 * (  - 2*mg**2*t1**(-1)*t2**(-1)*
     +    ts + 4*ms**2*t1**(-1)*t2**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts
     +     )
      QBV = QBV + SCC(4,1)*b_s * (  - 4 )
      QBV = QBV + SCC(4,1)*b_tx * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     - 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts )
      QBV = QBV + SCC(4,1)*b_ux * ( 2*m2**2*u1**(-1)*u2**(-1)*us
     +     - 2*ms**2*u1**(-1)*u2**(-1)*us - 2*u2**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tm * ( 8*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*ms**2*s**(-1)*ts - 4*m1**2*s**(-1)*u*ts
     +     - 4*m2**2*ms**2*s**(-1)*ts - 4*m2**2*s**(-1)*u*ts + 8*ms**2*
     +    s**(-1)*u*ts + 4*ms**2*ts + 4*u*ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * (  - 4*m1**2*m2**2
     +    *ms**2*t1**(-1)*t2**(-1)*ts - 3*m1**2*m2**2*s*t1**(-1)*
     +    t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*m2**2*
     +    t2**(-1)*ts - 11*m1**2*m2**4*t1**(-1)*t2**(-1)*ts - m1**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*t1**(-1)*ts + 2*
     +    m1**2*ms**2*t2**(-1)*ts + m1**2*s*t1**(-1)*ts + 3*m1**2*s*
     +    t2**(-1)*ts - 2*m1**2*u*t1**(-1)*ts - 2*m1**2*ts + 13*m1**4*
     +    m2**2*t1**(-1)*t2**(-1)*ts + 2*m1**4*ms**2*t1**(-1)*t2**(-1)*
     +    ts + 3*m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*t1**(-1)*ts - 5
     +    *m1**4*t2**(-1)*ts - 5*m1**6*t1**(-1)*t2**(-1)*ts - 3*m2**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + m2**2*ms**2*t1**(-1)*ts + 2*
     +    m2**2*ms**2*t2**(-1)*ts - 3*m2**2*s*t1**(-1)*ts + 4*m2**2*s*
     +    t2**(-1)*ts - m2**2*u*t1**(-1)*ts - 2*m2**2*ts + 2*m2**4*
     +    ms**2*t1**(-1)*t2**(-1)*ts - 4*m2**4*s*t1**(-1)*t2**(-1)*ts
     +     + 2*m2**4*t1**(-1)*ts - 3*m2**4*t2**(-1)*ts + 3*m2**6*
     +    t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * (  - 2*ms**2*s*
     +    t1**(-1)*ts - ms**2*s*t2**(-1)*ts + ms**2*u*t1**(-1)*ts + 
     +    ms**2*ts + 2*s*u*t1**(-1)*ts + 3*s*ts - u*ts - u**2*t1**(-1)*
     +    ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_um * ( 4*m1**2*ms**2*
     +    s**(-1)*us + 4*m1**2*s**(-1)*u*us + 8*m1**2*us - 4*m1**4*
     +    s**(-1)*us + 4*m2**2*ms**2*s**(-1)*us + 4*m2**2*s**(-1)*u*us
     +     + 8*m2**2*us - 4*m2**4*s**(-1)*us - 8*ms**2*s**(-1)*u*us - 4
     +    *ms**2*us - 4*s*us - 4*u*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * (  - 4*m1**2*m2**2
     +    *ms**2*u1**(-1)*u2**(-1)*us - 9*m1**2*m2**2*s*u1**(-1)*
     +    u2**(-1)*us - 3*m1**2*m2**2*u1**(-1)*us - 7*m1**2*m2**2*
     +    u2**(-1)*us + 7*m1**2*m2**4*u1**(-1)*u2**(-1)*us - m1**2*
     +    ms**2*s*u1**(-1)*u2**(-1)*us + 3*m1**2*ms**2*u1**(-1)*us + 2*
     +    m1**2*ms**2*u2**(-1)*us - 3*m1**2*s*u1**(-1)*us + 2*m1**2*s*
     +    u2**(-1)*us + m1**2*u*u2**(-1)*us - m1**2*us - 5*m1**4*m2**2*
     +    u1**(-1)*u2**(-1)*us + 2*m1**4*ms**2*u1**(-1)*u2**(-1)*us + 
     +    m1**6*u1**(-1)*u2**(-1)*us - 3*m2**2*ms**2*s*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*ms**2*u1**(-1)*us + m2**2*ms**2*
     +    u2**(-1)*us + 3*m2**2*s*u1**(-1)*us - 3*m2**2*s*u2**(-1)*us
     +     + m2**2*u*u1**(-1)*us - 2*m2**2*us + 2*m2**4*ms**2*u1**(-1)*
     +    u2**(-1)*us + 5*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**4*
     +    u1**(-1)*us + m2**4*u2**(-1)*us - 3*m2**6*u1**(-1)*u2**(-1)*
     +    us - 3*ms**2*s*u1**(-1)*us - ms**2*s*u2**(-1)*us - ms**2*u*
     +    u1**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * ( ms**2*u*u2**(-1)
     +    *us + s*u*u1**(-1)*us - s*us + s**2*u1**(-1)*us - s**2*
     +    u2**(-1)*us - u*us + u**2*u2**(-1)*us )
      QBV = QBV + SCC(5,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_tx * (  - m1**2*s**(-1)*t2**(-1)*ts
     +     + 2*m1**2*t1**(-1)*t2**(-1)*ts - 4*m2**2*t1**(-1)*t2**(-1)*
     +    ts + 2*ms**2*s**(-1)*t1**(-1)*ts + 2*ms**2*s**(-1)*t2**(-1)*
     +    ts + 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*s**(-1)*u*t1**(-1)*ts
     +     - s**(-1)*u*t2**(-1)*ts - 3*s**(-1)*ts - t1**(-1)*ts + 3*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCC(5,1)*b_ux * (  - 2*m1**2*s**(-1)*u2**(-1)*us
     +     - 2*m2**2*s**(-1)*u1**(-1)*us - 3*m2**2*s**(-1)*u2**(-1)*us
     +     - 2*m2**2*u1**(-1)*u2**(-1)*us + 2*ms**2*s**(-1)*u1**(-1)*us
     +     + 2*ms**2*s**(-1)*u2**(-1)*us + 2*ms**2*u1**(-1)*u2**(-1)*us
     +     + 3*s**(-1)*u*u2**(-1)*us - 3*s**(-1)*us - u1**(-1)*us + 3*
     +    u2**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_um2 * ( 8*m1**2*m2**2*
     +    s**(-1)*us + 4*m1**2*ms**2*s**(-1)*us - 4*m1**2*s**(-1)*u*us
     +     + 4*m2**2*ms**2*s**(-1)*us - 4*m2**2*s**(-1)*u*us - 8*ms**2*
     +    s**(-1)*u*us - 4*ms**2*us + 4*u*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 4*m1**2*
     +    m2**2*ms**2*u1**(-1)*u2**(-1)*us + 6*m1**2*m2**2*s*u1**(-1)*
     +    u2**(-1)*us + 6*m1**2*m2**2*u1**(-1)*us - 2*m1**2*m2**2*
     +    u2**(-1)*us + 10*m1**2*m2**4*u1**(-1)*u2**(-1)*us - 2*m1**2*
     +    ms**2*s*u1**(-1)*u2**(-1)*us + 2*m1**2*ms**2*u1**(-1)*us + 2*
     +    m1**2*ms**2*u2**(-1)*us - 6*m1**2*s*u1**(-1)*us - 8*m1**2*s*
     +    u2**(-1)*us - 2*m1**2*us - 8*m1**4*m2**2*u1**(-1)*u2**(-1)*us
     +     + 2*m1**4*ms**2*u1**(-1)*u2**(-1)*us - 8*m1**4*s*u1**(-1)*
     +    u2**(-1)*us + 2*m1**4*u2**(-1)*us + 2*m1**6*u1**(-1)*u2**(-1)
     +    *us - 2*m2**2*ms**2*s*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     +    u1**(-1)*us + 2*m2**2*ms**2*u2**(-1)*us - 2*m2**2*s*u1**(-1)*
     +    us - 6*m2**2*s*u2**(-1)*us - 4*m2**2*u*u1**(-1)*us + 4*m2**2*
     +    us + 2*m2**4*ms**2*u1**(-1)*u2**(-1)*us + 6*m2**4*s*u1**(-1)*
     +    u2**(-1)*us + 4*m2**4*u2**(-1)*us - 4*m2**6*u1**(-1)*u2**(-1)
     +    *us - 2*ms**2*s*u1**(-1)*us - 2*ms**2*s*u2**(-1)*us + 6*s*u*
     +    u1**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 6*s*us + 4*
     +    s**2*u1**(-1)*us - 2*u*us + 2*u**2*u1**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tm4 * (  - 4*m1**2*
     +    ms**2*s**(-1)*ts + 4*m1**2*s**(-1)*u*ts + 8*m1**2*ts - 4*
     +    m1**4*s**(-1)*ts - 4*m2**2*ms**2*s**(-1)*ts + 4*m2**2*s**(-1)
     +    *u*ts + 8*m2**2*ts - 4*m2**4*s**(-1)*ts + 8*ms**2*s**(-1)*u*
     +    ts + 4*ms**2*ts - 4*s*ts - 4*u*ts )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tx4 * (  - 4*m1**2*
     +    m2**2*ms**2*t1**(-1)*t2**(-1)*ts + 6*m1**2*m2**2*s*t1**(-1)*
     +    t2**(-1)*ts - 2*m1**2*m2**2*t1**(-1)*ts + 6*m1**2*m2**2*
     +    t2**(-1)*ts - 8*m1**2*m2**4*t1**(-1)*t2**(-1)*ts - 2*m1**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*t1**(-1)*ts + 2*
     +    m1**2*ms**2*t2**(-1)*ts - 6*m1**2*s*t1**(-1)*ts + 4*m1**2*s*
     +    t2**(-1)*ts + 2*m1**2*ts + 10*m1**4*m2**2*t1**(-1)*t2**(-1)*
     +    ts + 2*m1**4*ms**2*t1**(-1)*t2**(-1)*ts + 6*m1**4*s*t1**(-1)*
     +    t2**(-1)*ts + 4*m1**4*t1**(-1)*ts - 2*m1**4*t2**(-1)*ts - 4*
     +    m1**6*t1**(-1)*t2**(-1)*ts - 2*m2**2*ms**2*s*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**2*t1**(-1)*ts + 2*m2**2*ms**2*
     +    t2**(-1)*ts - 8*m2**2*s*t1**(-1)*ts - 4*m2**2*s*t2**(-1)*ts
     +     - 4*m2**2*u*t2**(-1)*ts - 4*m2**2*ts + 2*m2**4*ms**2*
     +    t1**(-1)*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*t2**(-1)*ts + 2*
     +    m2**4*t1**(-1)*ts + 2*m2**4*t2**(-1)*ts + 2*m2**6*t1**(-1)*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tx4 * (  - 2*ms**2*s*
     +    t1**(-1)*ts - 2*ms**2*s*t2**(-1)*ts - 2*s*u*t2**(-1)*ts - 4*s
     +    *ts + 2*u*ts + 2*u**2*t2**(-1)*ts )
      QBV = QBV + SCC(6,1)*b_um2 * ( 2*m1**2*s**(-2)*us + 2*m2**2*
     +    s**(-2)*us - 4*mg**2*s**(-2)*us + 4*ms**2*s**(-2)*us - 4*
     +    s**(-2)*u*us - 2*s**(-1)*us )
      QBV = QBV + SCC(6,1)*b_ux2 * ( 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 2*m1**2*ms**2*s**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*mg**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*
     +    m2**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 4*m2**2*u1**(-1)*
     +    u2**(-1)*us - 2*mg**2*u1**(-1)*u2**(-1)*us + 2*ms**2*s**(-1)*
     +    u1**(-1)*us + 2*ms**2*s**(-1)*u2**(-1)*us + 4*ms**2*u1**(-1)*
     +    u2**(-1)*us - 2*s*u1**(-1)*u2**(-1)*us - 2*u1**(-1)*us - 2*
     +    u2**(-1)*us )
      QBV = QBV + SCC(6,1)*b_tm4 * (  - 2*m1**2*s**(-2)*ts - 2*
     +    m2**2*s**(-2)*ts - 4*mg**2*s**(-2)*ts + 4*ms**2*s**(-2)*ts + 
     +    4*s**(-2)*u*ts + 2*s**(-1)*ts )
      QBV = QBV + SCC(6,1)*b_tx4 * ( 2*m1**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts - 2*m1**2*ms**2*s**(-1)*t1**(-1)*
     +    t2**(-1)*ts + 4*m1**2*t1**(-1)*t2**(-1)*ts + 2*m2**2*mg**2*
     +    s**(-1)*t1**(-1)*t2**(-1)*ts - 2*m2**2*ms**2*s**(-1)*t1**(-1)
     +    *t2**(-1)*ts - 2*mg**2*t1**(-1)*t2**(-1)*ts + 2*ms**2*s**(-1)
     +    *t1**(-1)*ts + 2*ms**2*s**(-1)*t2**(-1)*ts + 4*ms**2*t1**(-1)
     +    *t2**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts - 2*
     +    t2**(-1)*ts )
      QBV = QBV + SCD(1,1)*b_tm * (  - 4 )
      QBV = QBV + SCD(1,1)*b_tx * (  - 2 + 4*m1**2*ms**2*t1**(-1)*
     +    t2**(-1) + 2*m1**2*t1**(-1) - 2*m1**2*t2**(-1) - 2*m1**4*
     +    t1**(-1)*t2**(-1) - 2*ms**2*t1**(-1) + 2*ms**2*t2**(-1) - 2*
     +    ms**4*t1**(-1)*t2**(-1) )
      QBV = QBV + SCD(1,2)*b_um * (  - 4 )
      QBV = QBV + SCD(1,2)*b_ux * (  - 3 + m1**2*u2**(-1) + 4*
     +    m2**2*ms**2*u1**(-1)*u2**(-1) - 2*m2**2*u1**(-1) + 2*m2**2*
     +    u2**(-1) - 2*m2**4*u1**(-1)*u2**(-1) + 2*ms**2*u1**(-1) - 2*
     +    ms**2*u2**(-1) - 2*ms**4*u1**(-1)*u2**(-1) - s*u2**(-1) - t*
     +    u2**(-1) )
      QBV = QBV + SCD(2,1)*b_um2 * (  - 2*m1**2*mg**2*s**(-2)*
     +    tg**(-1)*us + 2*m1**2*ms**2*s**(-2)*tg**(-1)*us - 2*m2**2*
     +    mg**2*s**(-2)*tg**(-1)*us + 2*m2**2*ms**2*s**(-2)*tg**(-1)*us
     +     + 4*mg**2*ms**2*s**(-2)*tg**(-1)*us + 4*mg**2*s**(-2)*us - 4
     +    *ms**2*s**(-2)*us - 4*ms**4*s**(-2)*tg**(-1)*us + 2*s**(-1)*
     +    us )
      QBV = QBV + SCD(2,1)*b_ux2 * (  - 4*m1**2*mg**2*ms**2*
     +    s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 4*m1**2*mg**2*
     +    u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m1**2*mg**4*s**(-1)*
     +    u1**(-1)*tg**(-1)*u2**(-1)*us - 4*m1**2*ms**2*u1**(-1)*
     +    tg**(-1)*u2**(-1)*us + 2*m1**2*ms**4*s**(-1)*u1**(-1)*
     +    tg**(-1)*u2**(-1)*us + 2*m1**2*s*u1**(-1)*tg**(-1)*u2**(-1)*
     +    us - 4*m2**2*mg**2*ms**2*s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*
     +    us + 2*m2**2*mg**2*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*
     +    mg**4*s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*ms**4*
     +    s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*mg**2*ms**2*s**(-1)
     +    *u1**(-1)*tg**(-1)*us + 2*mg**2*ms**2*s**(-1)*tg**(-1)*
     +    u2**(-1)*us + 8*mg**2*ms**2*u1**(-1)*tg**(-1)*u2**(-1)*us - 4
     +    *mg**2*s*u1**(-1)*tg**(-1)*u2**(-1)*us - 2*mg**2*u1**(-1)*
     +    tg**(-1)*us - 2*mg**4*u1**(-1)*tg**(-1)*u2**(-1)*us + 8*ms**2
     +    *s*u1**(-1)*tg**(-1)*u2**(-1)*us + 6*ms**2*u1**(-1)*tg**(-1)*
     +    us )
      QBV = QBV + SCD(2,1)*b_ux2 * ( 2*ms**2*tg**(-1)*u2**(-1)*us
     +     - 2*ms**4*s**(-1)*u1**(-1)*tg**(-1)*us - 2*ms**4*s**(-1)*
     +    tg**(-1)*u2**(-1)*us - 8*ms**4*u1**(-1)*tg**(-1)*u2**(-1)*us
     +     - 2*s*u1**(-1)*tg**(-1)*us - 2*s**2*u1**(-1)*tg**(-1)*
     +    u2**(-1)*us )
      QBV = QBV + SCD(2,2)*b_tm4 * (  - 2*m1**2*mg**2*s**(-2)*
     +    ug**(-1)*ts + 2*m1**2*ms**2*s**(-2)*ug**(-1)*ts - 2*m2**2*
     +    mg**2*s**(-2)*ug**(-1)*ts + 2*m2**2*ms**2*s**(-2)*ug**(-1)*ts
     +     + 4*mg**2*ms**2*s**(-2)*ug**(-1)*ts + 4*mg**2*s**(-2)*ts - 4
     +    *ms**2*s**(-2)*ts - 4*ms**4*s**(-2)*ug**(-1)*ts + 2*s**(-1)*
     +    ts )
      QBV = QBV + SCD(2,2)*b_tx4 * (  - 4*m1**2*mg**2*ms**2*
     +    s**(-1)*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*mg**2*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*mg**4*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*ms**4*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts - 4*m2**2*mg**2*ms**2*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 4*m2**2*mg**2*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts + 2*m2**2*mg**4*s**(-1)*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 4*m2**2*ms**2*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**4*s**(-1)*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts + 2*m2**2*s*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*
     +    mg**2*ms**2*s**(-1)*t1**(-1)*ug**(-1)*ts + 2*mg**2*ms**2*
     +    s**(-1)*ug**(-1)*t2**(-1)*ts + 8*mg**2*ms**2*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 4*mg**2*s*t1**(-1)*ug**(-1)*t2**(-1)*
     +    ts - 2*mg**2*ug**(-1)*t2**(-1)*ts - 2*mg**4*t1**(-1)*ug**(-1)
     +    *t2**(-1)*ts + 8*ms**2*s*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*
     +    ms**2*t1**(-1)*ug**(-1)*ts )
      QBV = QBV + SCD(2,2)*b_tx4 * ( 6*ms**2*ug**(-1)*t2**(-1)*ts
     +     - 2*ms**4*s**(-1)*t1**(-1)*ug**(-1)*ts - 2*ms**4*s**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 8*ms**4*t1**(-1)*ug**(-1)*t2**(-1)*ts
     +     - 2*s*ug**(-1)*t2**(-1)*ts - 2*s**2*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts )

c               function has be real at the end !!!
c               normalization factor !!!
      LL_QBV = real(QBV) / 8.D0

c               the prefactors removed in the form program 
c                   except for alphas, which is cut 
      LL_QBV = LL_QBV * Cf /Pi

c               the averaging factors
      LL_QBV = LL_QBV /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      LL_QBV = LL_QBV * (abs(m1)+abs(m2))**2/4.D0 /4.D0/Pi

      end


