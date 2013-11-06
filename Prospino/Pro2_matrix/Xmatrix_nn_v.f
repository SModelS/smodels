cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NN_QBB(IQ,IOUT,MASSIN,CS,CT,CL,CR)                               c
c                                                                      c
c     NN_QBV(IQ,IOUT,MASSIN,CS,CT,CL,CR,CV)                            c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       IQ   = +-1 FOR UP AND DOWN QUARKS                              c
c                                                                      c
c       IOUT = 1,2,3,4 FOR NN,C+C-,NC+,NC-                             c
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
c       CT(1:4)  TYPICAL S*T CHANNEL COUPLINGS [COMPLEX]               c
c       CL(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CR(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CV(1:4)  FRACTIONS OF COUPLINGS NEEDED FOR VIRTUAL [COMPLEX]   c
c                -> CLOT,CUPT,CUPU,CLOU                                c
c                                                                      c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ctp: note that only interference terms are not commented out!!

c --------------------------------------------------------------------
c     momentum assignment attached to incoming q:    k1-p1-C(1) with C(1) the chargino is t channel
c                                              qbar: k2-p2-C(2) with C(2) the neutralino is t channel
c               nqs/cqs coupings for NC+ (u-dbar):
c                t channel with C(1)-chargino-u-sd;    C(2)-neutralino-dbar-sd; s-down
c                u channel with C(3)-chargino-dbar-su; C(4)-neutralino-u-sd;    s-up  
c               nqs/cqs coupings for NC- (d-ubar):
c                t channel with C(1)-chargino-d-su;    C(2)-neutralino-ubar-su; s-up
c                u channel with C(3)-chargino-ubar-sd; C(4)-neutralino-d-sd;    s-down

      real*8 function NN_QBB_NG(iq,iout,massin,
     &                          Csdum,Ctdum,Cldum,Crdum,mst,msu)

      implicit none 

      integer    iq,iout,n
      real*8     Csdum(4),massin(1:20),Cs(4),Nc,only_u,only_d
     &          ,s,m1,m2,mx,mst(-1:1),msu(-1:1)
     &          ,t1,u1,t2,u2,tsl,usl,tsr,usr,sz
      complex*16 Ctdum(4),Cldum(4),Crdum(4)
     &          ,Ct(4),Ctc(4),Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_s,QBB_tx,QBB_tm,QBB_ux,QBB_um 

      Nc = 3.D0

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
      m1 = massin(6)
      m2 = massin(7)
      mx = massin(8)

      tsl = t2 + m2**2 - mst(-1)**2 
      tsr = t2 + m2**2 - mst(+1)**2 
      usl = u2 + m2**2 - msu(-1)**2 
      usr = u2 + m2**2 - msu( 1)**2 
      sz = s - mx**2 

      do n=1,4 
         Cs(n) = Csdum(n)
         Ct(n) = Ctdum(n)
         Cl(n) = Cldum(n)
         Cr(n) = Crdum(n)
      end do

c               the quark flavors in the t and u channel 
      only_u = ( 1.D0 + iq )/2.D0
      only_d = ( 1.D0 - iq )/2.D0

c               CC case : u,d type quarks in t,u channel 
      if (iout.eq.2) then 
         Cl(1) = only_u * Cl(1) 
         Cl(2) = only_u * Cl(2) 
         Cl(3) = only_d * Cl(3) 
         Cl(4) = only_d * Cl(4) 
c               CN+ case : only u,dbar incoming state 
      else if (iout.eq.3) then 
         do n=1,4 
            Cs(n) = only_u * Cs(n)
            Ct(n) = only_u * Ct(n)
            Cl(n) = only_u * Cl(n)
            Cr(n) = only_u * Cr(n)
         end do
c               CN- case : only ubar,d incoming state 
      else if (iout.eq.4) then 
         do n=1,4 
            Cs(n) = only_d * Cs(n)
            Ct(n) = only_d * Ct(n)
            Cl(n) = only_d * Cl(n)
            Cr(n) = only_d * Cr(n)
         end do
      end if 

      do n=1,4 
         Ctc(n) =  conjg( Ct(n) )
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      QBB_s  = dcmplx(0.D0,0.D0)
      QBB_tm = dcmplx(0.D0,0.D0)
      QBB_tx = dcmplx(0.D0,0.D0)
      QBB_um = dcmplx(0.D0,0.D0)
      QBB_ux = dcmplx(0.D0,0.D0)

c               the form output 
      QBB_s =         s**(-2) * (
     +     + 16*Cs(1)*t1*t2
     +     + 16*Cs(1)*u1*u2
     +     )
      QBB_s = QBB_s + s**(-1)* (                     ! s*t channel: C(1)*Cc(2)
     +     + 16*Cl(1)*Clc(2)*Ct(1)*t1*t2/tsl
     +     + 16*Cr(1)*Crc(2)*Ct(3)*t1*t2/tsr
     +     )
      QBB_s = QBB_s + s**(-1)* (                     ! s*u channel: C(4)*Cc(3)
     +     - 16*Cl(4)*Clc(3)*Ct(2)*u1*u2/usl
     +     - 16*Cr(4)*Crc(3)*Ct(4)*u1*u2/usr
     +     )
      QBB_s = QBB_s + s**(-1)*sz**(-1) * (
     +     - 32*Cs(3)*t1*t2
     +     + 32*Cs(3)*u1*u2
     +     )
      QBB_s = QBB_s + s**(-1) * (
     +     + 32*Cs(2)*m1*m2
     +     )
      QBB_s = QBB_s + (                              ! s*t channel: C(1)*Cc(2)
     +     + 16*Cl(1)*Clc(2)*Ct(2)*m1*m2/tsl
     +     + 16*Cr(1)*Crc(2)*Ct(4)*m1*m2/tsr
     +     )
      QBB_s = QBB_s + (                              ! s*u channel: C(4)*Cc(3)
     +     - 16*Cl(4)*Clc(3)*Ct(1)*m1*m2/usl
     +     - 16*Cr(4)*Crc(3)*Ct(3)*m1*m2/usr
     +     )

      QBB_tx =          s**(-1)* (                   ! t*s channel: C(2)*Cc(1)
     +      + 16*Cl(2)*Clc(1)*Ctc(1)*t1*t2/tsl
     +      + 16*Cr(2)*Crc(1)*Ctc(3)*t1*t2/tsr
     +     )
      QBB_tx = QBB_tx + (                            ! t*t channel: C(1)*Cc(2)*C(2)*Cc(1)
     +      + 32*Cl(1)*Cl(2)*Clc(1)*Clc(2)*t1*t2/tsl**2
     +      + 32*Cr(1)*Cr(2)*Crc(1)*Crc(2)*t1*t2/tsr**2
     +     )

      QBB_tm =          s* (                         ! t*u channel: C(2)*Cc(1)*C(4)*Cc(3)
     +      - 32*Cl(2)*Cl(4)*Clc(1)*Clc(3)*m1*m2/tsl/usl
     +      - 32*Cr(2)*Cr(4)*Crc(1)*Crc(3)*m1*m2/tsr/usr
     +     )
      QBB_tm = QBB_tm + (                            ! t*s channel: C(2)*Cc(1)
     +      + 16*Cl(2)*Clc(1)*Ctc(2)*m1*m2/tsl
     +      + 16*Cr(2)*Crc(1)*Ctc(4)*m1*m2/tsr
     +     )

      QBB_ux =          s**(-1)* (                   ! u*s channel: C(3)*Cc(4)
     +      - 16*Cl(3)*Clc(4)*Ctc(2)*u1*u2/usl
     +      - 16*Cr(3)*Crc(4)*Ctc(4)*u1*u2/usr
     +     )
      QBB_ux = QBB_ux + (                            ! u*u channel: C(3)*Cc(4)*C(4)*Cc(3)
     +      + 32*Cl(3)*Cl(4)*Clc(3)*Clc(4)*u1*u2/usl**2
     +      + 32*Cr(3)*Cr(4)*Crc(3)*Crc(4)*u1*u2/usr**2
     +     )

      QBB_um =          s* (                         ! t*u channel: C(1)*Cc(2)*C(3)*Cc(4)
     +      - 32*Cl(1)*Cl(3)*Clc(2)*Clc(4)*m1*m2/tsl/usl
     +      - 32*Cr(1)*Cr(3)*Crc(2)*Crc(4)*m1*m2/tsr/usr
     +     )
      QBB_um = QBB_um + (                            ! s*u channel: C(3)*Cc(4)
     +      - 16*Cl(3)*Clc(4)*Ctc(1)*m1*m2/usl
     +      - 16*Cr(3)*Crc(4)*Ctc(3)*m1*m2/usr
     +     )

      NN_QBB_NG = real( QBB_s + QBB_tx + QBB_tm + QBB_ux + QBB_um )/2.D0
c               the prefactors removed in the form program 
c                   except for alphas, which is cut 
      NN_QBB_NG = NN_QBB_NG * Nc 

c               the averaging factors
      NN_QBB_NG = NN_QBB_NG /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NN_QBB_NG = NN_QBB_NG * (abs(m1)+abs(m2))**2/4.D0

      end


c --------------------------------------------------------------------
c remember squark masses massin(11) in ts and us need not be the same
      real*8 function NN_QBB(iq,iout,massin,Csdum,Ctdum,Cldum,Crdum)

      implicit none 

      integer    iq,iout,n
      real*8     Csdum(4),massin(1:20),Cs(4),Nc,only_u,only_d
     &          ,s,m1,m2,mx,ms,t1,u1,t2,u2,ts,us,sz
      complex*16 Ctdum(4),Cldum(4),Crdum(4)
     &          ,Ct(4),Ctc(4),Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_s,QBB_tx,QBB_tm,QBB_ux,QBB_um 

      Nc = 3.D0

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
      m1 = massin(6)
      m2 = massin(7)
      mx = massin(8)
      ms = massin(11)

      ts = t2 + m2**2 - ms**2 
      us = u2 + m2**2 - ms**2 
      sz = s - mx**2 

      do n=1,4 
         Cs(n) = Csdum(n)
         Ct(n) = Ctdum(n)
         Cl(n) = Cldum(n)
         Cr(n) = Crdum(n)
      end do

c               the quark flavors in the t and u channel 
      only_u = ( 1.D0 + iq )/2.D0
      only_d = ( 1.D0 - iq )/2.D0

c               CC case : u,d type quarks in t,u channel 
      if (iout.eq.2) then 
         Cl(1) = only_u * Cl(1) 
         Cl(2) = only_u * Cl(2) 
         Cl(3) = only_d * Cl(3) 
         Cl(4) = only_d * Cl(4) 
c               CN+ case : only u,dbar incoming state 
      else if (iout.eq.3) then 
         do n=1,4 
            Cs(n) = only_u * Cs(n)
            Ct(n) = only_u * Ct(n)
            Cl(n) = only_u * Cl(n)
            Cr(n) = only_u * Cr(n)
         end do
c               CN- case : only ubar,d incoming state 
      else if (iout.eq.4) then 
         do n=1,4 
            Cs(n) = only_d * Cs(n)
            Ct(n) = only_d * Ct(n)
            Cl(n) = only_d * Cl(n)
            Cr(n) = only_d * Cr(n)
         end do
      end if 

      do n=1,4 
         Ctc(n) =  conjg( Ct(n) )
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      QBB_s  = dcmplx(0.D0,0.D0)
      QBB_tm = dcmplx(0.D0,0.D0)
      QBB_tx = dcmplx(0.D0,0.D0)
      QBB_um = dcmplx(0.D0,0.D0)
      QBB_ux = dcmplx(0.D0,0.D0)

c               the form output 
      QBB_s =         s**(-2) * (
     +     + 16*Cs(1)*t1*t2
     +     + 16*Cs(1)*u1*u2
     +     )
      QBB_s = QBB_s + s**(-1)*ts**(-1) * (           ! s*t channel: C(1)*Cc(2)
     +     + 16*Cl(1)*Clc(2)*Ct(1)*t1*t2
     +     + 16*Cr(1)*Crc(2)*Ct(3)*t1*t2
     +     )
      QBB_s = QBB_s + s**(-1)*us**(-1) * (           ! s*u channel: C(4)*Cc(3)
     +     - 16*Cl(4)*Clc(3)*Ct(2)*u1*u2
     +     - 16*Cr(4)*Crc(3)*Ct(4)*u1*u2
     +     )
      QBB_s = QBB_s + s**(-1)*sz**(-1) * (
     +     - 32*Cs(3)*t1*t2
     +     + 32*Cs(3)*u1*u2
     +     )
      QBB_s = QBB_s + s**(-1) * (
     +     + 32*Cs(2)*m1*m2
     +     )
      QBB_s = QBB_s + ts**(-1) * (                   ! s*t channel: C(1)*Cc(2)
     +     + 16*Cl(1)*Clc(2)*Ct(2)*m1*m2
     +     + 16*Cr(1)*Crc(2)*Ct(4)*m1*m2
     +     )
      QBB_s = QBB_s + us**(-1) * (                   ! s*u channel: C(4)*Cc(3)
     +     - 16*Cl(4)*Clc(3)*Ct(1)*m1*m2
     +     - 16*Cr(4)*Crc(3)*Ct(3)*m1*m2
     +     )

      QBB_tx =          s**(-1)*ts**(-1) * (         ! t*s channel: C(2)*Cc(1)
     +      + 16*Cl(2)*Clc(1)*Ctc(1)*t1*t2
     +      + 16*Cr(2)*Crc(1)*Ctc(3)*t1*t2
     +     )
      QBB_tx = QBB_tx + ts**(-2) * (                 ! t*t channel: C(1)*Cc(2)*C(2)*Cc(1)
     +      + 32*Cl(1)*Cl(2)*Clc(1)*Clc(2)*t1*t2
     +      + 32*Cr(1)*Cr(2)*Crc(1)*Crc(2)*t1*t2
     +     )

      QBB_tm =          s*ts**(-1)*us**(-1) * (      ! t*u channel: C(2)*Cc(1)*C(4)*Cc(3)
     +      - 32*Cl(2)*Cl(4)*Clc(1)*Clc(3)*m1*m2
     +      - 32*Cr(2)*Cr(4)*Crc(1)*Crc(3)*m1*m2
     +     )
      QBB_tm = QBB_tm + ts**(-1) * (                 ! t*s channel: C(2)*Cc(1)
     +      + 16*Cl(2)*Clc(1)*Ctc(2)*m1*m2
     +      + 16*Cr(2)*Crc(1)*Ctc(4)*m1*m2
     +     )

      QBB_ux =          s**(-1)*us**(-1) * (         ! u*s channel: C(3)*Cc(4)
     +      - 16*Cl(3)*Clc(4)*Ctc(2)*u1*u2
     +      - 16*Cr(3)*Crc(4)*Ctc(4)*u1*u2
     +     )
      QBB_ux = QBB_ux + us**(-2) * (                 ! u*u channel: C(3)*Cc(4)*C(4)*Cc(3)
     +      + 32*Cl(3)*Cl(4)*Clc(3)*Clc(4)*u1*u2
     +      + 32*Cr(3)*Cr(4)*Crc(3)*Crc(4)*u1*u2
     +     )

      QBB_um =          s*ts**(-1)*us**(-1) * (      ! t*u channel: C(1)*Cc(2)*C(3)*Cc(4)
     +      - 32*Cl(1)*Cl(3)*Clc(2)*Clc(4)*m1*m2
     +      - 32*Cr(1)*Cr(3)*Crc(2)*Crc(4)*m1*m2
     +     )
      QBB_um = QBB_um + us**(-1) * (                 ! s*u channel: C(3)*Cc(4)
     +      - 16*Cl(3)*Clc(4)*Ctc(1)*m1*m2
     +      - 16*Cr(3)*Crc(4)*Ctc(3)*m1*m2
     +     )

      NN_QBB = real( QBB_s + QBB_tx + QBB_tm + QBB_ux + QBB_um )/2.D0
c               the prefactors removed in the form program 
c                   except for alphas, which is cut 
      NN_QBB = NN_QBB * Nc 

c               the averaging factors
      NN_QBB = NN_QBB /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NN_QBB = NN_QBB * (abs(m1)+abs(m2))**2/4.D0

      end


c --------------------------------------------------------------------
      subroutine BORN_PARTS_NN(iq,iout,massin,Cs,Ct,Cl,Cr,
     &                b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2)

      implicit none 

      integer    iq,iout,n
      real*8     massin(1:20),Cs(4),Nc,only_u,only_d
     &          ,s,m1,m2,mx,ms,t1,u1,t2,u2,ts,us,sz
     &          ,b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2
      complex*16 Ct(4),Ctc(4),Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_ss,QBB_st,QBB_su
     &          ,QBB_txs,QBB_txt,QBB_tmu,QBB_tms
     &          ,QBB_uxs,QBB_uxu,QBB_umt,QBB_ums 

      Nc = 3.D0

      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
      m1 = massin(6)
      m2 = massin(7)
      mx = massin(8)
      ms = massin(11)

      ts = t2 + m2**2 - ms**2 
      us = u2 + m2**2 - ms**2 
      sz = s - mx**2 

      do n=1,4 
         Ctc(n) =  conjg( Ct(n) )
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

c               the form output 
      QBB_ss =          s**(-2) * (
     +     + 16*Cs(1)*t1*t2
     +     + 16*Cs(1)*u1*u2
     +     )
      QBB_ss = QBB_ss + s**(-1) * (
     +     + 32*Cs(2)*m1*m2
     +     )
      QBB_ss = QBB_ss + s**(-1)*sz**(-1) * (
     +     - 32*Cs(3)*t1*t2
     +     + 32*Cs(3)*u1*u2
     +     )

      QBB_st =          s**(-1)*ts**(-1) * (
     +     + 16*Cl(1)*Clc(2)*Ct(1)*t1*t2
     +     + 16*Cr(1)*Crc(2)*Ct(3)*t1*t2
     +     )
      QBB_st = QBB_st + ts**(-1) * (
     +     + 16*Cl(1)*Clc(2)*Ct(2)*m1*m2
     +     + 16*Cr(1)*Crc(2)*Ct(4)*m1*m2
     +     )

      QBB_su =          s**(-1)*us**(-1) * (
     +     - 16*Cl(4)*Clc(3)*Ct(2)*u1*u2
     +     - 16*Cr(4)*Crc(3)*Ct(4)*u1*u2
     +     )
      QBB_su = QBB_su + us**(-1) * (
     +     - 16*Cl(4)*Clc(3)*Ct(1)*m1*m2
     +     - 16*Cr(4)*Crc(3)*Ct(3)*m1*m2
     +     )

      QBB_txs =         s**(-1)*ts**(-1) * (
     +      + 16*Cl(2)*Clc(1)*Ctc(1)*t1*t2
     +      + 16*Cr(2)*Crc(1)*Ctc(3)*t1*t2
     +     )

      QBB_txt =         ts**(-2) * (
     +      + 32*Cl(1)*Cl(2)*Clc(1)*Clc(2)*t1*t2
     +      + 32*Cr(1)*Cr(2)*Crc(1)*Crc(2)*t1*t2
     +     )

      QBB_tmu =         s*ts**(-1)*us**(-1) * (
     +      - 32*Cl(2)*Cl(4)*Clc(1)*Clc(3)*m1*m2
     +      - 32*Cr(2)*Cr(4)*Crc(1)*Crc(3)*m1*m2
     +     )

      QBB_tms =         ts**(-1) * (
     +      + 16*Cl(2)*Clc(1)*Ctc(2)*m1*m2
     +      + 16*Cr(2)*Crc(1)*Ctc(4)*m1*m2
     +     )

      QBB_uxs =         s**(-1)*us**(-1) * (
     +      - 16*Cl(3)*Clc(4)*Ctc(2)*u1*u2
     +      - 16*Cr(3)*Crc(4)*Ctc(4)*u1*u2
     +     )

      QBB_uxu =         us**(-2) * (
     +      + 32*Cl(3)*Cl(4)*Clc(3)*Clc(4)*u1*u2
     +      + 32*Cr(3)*Cr(4)*Crc(3)*Crc(4)*u1*u2
     +     )

      QBB_umt =         s*ts**(-1)*us**(-1) * (
     +      - 32*Cl(1)*Cl(3)*Clc(2)*Clc(4)*m1*m2
     +      - 32*Cr(1)*Cr(3)*Crc(2)*Crc(4)*m1*m2
     +     )

      QBB_ums =         us**(-1) * (
     +      - 16*Cl(3)*Clc(4)*Ctc(1)*m1*m2
     +      - 16*Cr(3)*Crc(4)*Ctc(3)*m1*m2
     +     )

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
      
c               the quark flavors in the t and u channel 
      only_u = ( 1.D0 + iq )/2.D0
      only_d = ( 1.D0 - iq )/2.D0

c               NN case : anything goes
      if (iout.eq.1) then 
         b_s   = Nc*real( QBB_ss  + QBB_st   + QBB_su )
         b_tx  = Nc*real( QBB_txs + QBB_txt )
         b_tm  = Nc*real( QBB_tmu + QBB_tms )
         b_ux  = Nc*real( QBB_uxs + QBB_uxu )
         b_um  = Nc*real( QBB_umt + QBB_ums )
         b_tx4 = Nc*real( QBB_txs + QBB_txt )
         b_tm4 = Nc*real( QBB_tmu + QBB_tms )
         b_ux2 = Nc*real( QBB_uxs + QBB_uxu )
         b_um2 = Nc*real( QBB_umt + QBB_ums )
c               CC case : u,d type quarks in t,u channel 
      else if (iout.eq.2) then 
         b_s   = Nc*real( QBB_ss  + only_u*QBB_st   + only_d*QBB_su )
         b_tx  = only_u*Nc*real( QBB_txs + QBB_txt )
         b_tm  = only_u*Nc*real( only_d*QBB_tmu + QBB_tms )
         b_ux  = only_d*Nc*real( QBB_uxs + QBB_uxu )
         b_um  = only_d*Nc*real( only_u*QBB_umt + QBB_ums )
         b_tx4 = only_d*Nc*real( QBB_txs + only_u*QBB_txt )
         b_tm4 = only_d*Nc*real( QBB_tmu + QBB_tms )
         b_ux2 = only_u*Nc*real( QBB_uxs + only_d*QBB_uxu )
         b_um2 = only_u*Nc*real( QBB_umt + QBB_ums )
c               CN+ case : only u,dbar incoming state 
      else if (iout.eq.3) then 
         b_s   = only_u*Nc*real( QBB_ss  + QBB_st   + QBB_su )
         b_tx  = only_u*Nc*real( QBB_txs + QBB_txt )
         b_tm  = only_u*Nc*real( QBB_tmu + QBB_tms )
         b_ux  = only_u*Nc*real( QBB_uxs + QBB_uxu )
         b_um  = only_u*Nc*real( QBB_umt + QBB_ums )
         if ((Cl(2).eq.0.D0).or.(Clc(4).eq.0.D0)) then
          b_tx4 = 0.D0
          b_tm4 = 0.D0
          b_ux2 = 0.D0
          b_um2 = 0.D0
         else
          b_tx4 = only_u*Nc*real( Cl(4)/Cl(2)   * (QBB_txs + QBB_txt) )
          b_tm4 = only_u*Nc*real( Cl(4)/Cl(2)   * (QBB_tmu + QBB_tms) )
          b_ux2 = only_u*Nc*real( Clc(2)/Clc(4) * (QBB_uxs + QBB_uxu) )
          b_um2 = only_u*Nc*real( Clc(2)/Clc(4) * (QBB_umt + QBB_ums) )
         end if
c               CN- case : only ubar,d incoming state 
      else if (iout.eq.4) then 
         b_s   = only_d*Nc*real( QBB_ss  + QBB_st   + QBB_su )
         b_tx  = only_d*Nc*real( QBB_txs + QBB_txt )
         b_tm  = only_d*Nc*real( QBB_tmu + QBB_tms )
         b_ux  = only_d*Nc*real( QBB_uxs + QBB_uxu )
         b_um  = only_d*Nc*real( QBB_umt + QBB_ums )
         if ((Cl(2).eq.0.D0).or.(Clc(4).eq.0.D0)) then
          b_tx4 = 0.D0
          b_tm4 = 0.D0
          b_ux2 = 0.D0
          b_um2 = 0.D0
         else
          b_tx4 = only_d*Nc*real( Clc(4)/Clc(2) * (QBB_txs + QBB_txt) )
          b_tm4 = only_d*Nc*real( Clc(4)/Clc(2) * (QBB_tmu + QBB_tms) )
          b_ux2 = only_d*Nc*real( Cl(2)/Cl(4)   * (QBB_uxs + QBB_uxu) )
          b_um2 = only_d*Nc*real( Cl(2)/Cl(4)   * (QBB_umt + QBB_ums) )
         end if 
      end if
 
      end


c --------------------------------------------------------------------
      real*8 function NN_QBV(iq,iout,massin,Csdum,Ctdum,Cldum,Crdum,Cv)

      implicit none 

      integer    iq,iout,n
      real*8     Csdum(4),massin(1:20),Pi,Nc,Cf,finite,msbar
     &          ,s,t,u,m1,m2,mg,ms,mav2
     &          ,t1,u1,t2,u2,ts,us,tg,ug,kaellen,Cs(4)
     &          ,SCB(1:7,1:5),SCBP(1),SCC(1:6,1:4),SCD(1:2,1:2)
     &          ,b_s,b_t,b_u,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2
      complex*16 Ctdum(4),Cldum(4),Crdum(4)
     &          ,Ct(4),Cl(4),Cr(4),Cv(4),QBV,QBS
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

c               MSbar coupling shift to maintain SUSY in qsN couplings
      msbar = 1.d0

c               the scalar functions 
      call SCALAR_ARRAY_NN(massin,SCB,SCBP,SCC,SCD)

      do n=1,4 
         Cs(n) = Csdum(n)
         Ct(n) = Ctdum(n)
         Cl(n) = Cldum(n)
         Cr(n) = Crdum(n)
      end do

c               the born type structures 
      call BORN_PARTS_NN(iq,iout,massin,Cs,Ct,Cl,Cr,
     &                b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2)
      b_t  = b_tx  + b_tm
      b_u  = b_ux  + b_um

      QBV =
     &  + b_u * (  - 2*msbar )
      QBV = QBV + b_t * (  - 2*msbar )
      QBV = QBV + b_s * (  - 4 )
      QBV = QBV + finite*b_u * ( 4 )
      QBV = QBV + finite*b_t * ( 4 )
      QBV = QBV + finite*b_s * ( 4 )
      QBV = QBV + SCB(3,1)*b_tx * ( t1**(-1)*ts + 3*t2**(-1)*ts + 
     &    m2**2*t1**(-1)*t2**(-1)*ts + 3*m1**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCB(3,1)*b_t * (  - 4 - 4*u*ts**(-1) - 4*s*
     &    ts**(-1) + 4*ms**2*ts**(-1) - 4*m2**2*t2**(-1) + 4*m2**2*
     &    ts**(-1) - 4*m1**2*t1**(-1) + 4*m1**2*ts**(-1) )
      QBV = QBV + SCB(3,2)*b_ux * ( u1**(-1)*us + 3*u2**(-1)*us + 
     &    m2**2*u1**(-1)*u2**(-1)*us + 3*m1**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,2)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCB(3,2)*b_u * (  - 4 + 4*u*us**(-1) + 4*ms**2*
     &    us**(-1) - 4*m2**2*u2**(-1) - 4*m1**2*u1**(-1) )
      QBV = QBV + SCB(3,3)*b_u * (  - 8*ms**2*us**(-1) )
      QBV = QBV + SCB(3,3)*b_t * (  - 8*ms**2*ts**(-1) )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_ux * (  - 2*m2**2*s*
     &    u1**(-1)*us + 2*m2**2*s*u2**(-1)*us + m2**4*u1**(-1)*us - 
     &    m2**4*u2**(-1)*us - 2*m2**4*s*u1**(-1)*u2**(-1)*us + m2**6*
     &    u1**(-1)*u2**(-1)*us + 4*m1**2*s*u1**(-1)*us + 2*m1**2*m2**2*
     &    u1**(-1)*us + 2*m1**2*m2**2*u2**(-1)*us + 10*m1**2*m2**2*s*
     &    u1**(-1)*u2**(-1)*us - 3*m1**2*m2**4*u1**(-1)*u2**(-1)*us - 3
     &    *m1**4*u1**(-1)*us - m1**4*u2**(-1)*us + 3*m1**4*m2**2*
     &    u1**(-1)*u2**(-1)*us - m1**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_um * ( 4*u*us + 4*s*us
     &     - 4*m2**2*s**(-1)*u*us - 8*m2**2*us + 4*m2**4*s**(-1)*us + 4
     &    *m1**2*s**(-1)*u*us - 4*m1**2*us - 4*m1**2*m2**2*s**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tx * ( 2*m1**2*s*
     &    t1**(-1)*ts + 2*m1**2*s*t2**(-1)*ts + 4*m1**2*m2**2*t2**(-1)*
     &    ts + 6*m1**2*m2**2*s*t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**4*
     &    t1**(-1)*t2**(-1)*ts - 4*m1**4*t2**(-1)*ts + 2*m1**4*s*
     &    t1**(-1)*t2**(-1)*ts + 8*m1**4*m2**2*t1**(-1)*t2**(-1)*ts - 4
     &    *m1**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tm * (  - 4*u*ts - 4*s*
     &    ts + 4*m2**2*s**(-1)*u*ts + 8*m2**2*ts - 4*m2**4*s**(-1)*ts
     &     - 4*m1**2*s**(-1)*u*ts + 4*m1**2*ts + 4*m1**2*m2**2*s**(-1)*
     &    ts )
      QBV = QBV + SCB(3,4)*b_ux * (  - u1**(-1)*us + u2**(-1)*us
     &     - m2**2*u1**(-1)*u2**(-1)*us + m1**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*b_tx * (  - t1**(-1)*ts + t2**(-1)*ts
     &     - m2**2*t1**(-1)*t2**(-1)*ts + m1**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_tm * ( 4*s**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_u * ( 4*m1**2*u1**(-1) )
      QBV = QBV + SCB(3,4)*b_t * ( 4*m1**2*t1**(-1) )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_ux * ( 6*m2**2*s*
     &    u1**(-1)*us - 2*m2**2*s*u2**(-1)*us - 4*m2**4*u1**(-1)*us + 6
     &    *m2**4*s*u1**(-1)*u2**(-1)*us - 4*m2**6*u1**(-1)*u2**(-1)*us
     &     + 4*m1**2*m2**2*u1**(-1)*us + 2*m1**2*m2**2*s*u1**(-1)*
     &    u2**(-1)*us + 8*m1**2*m2**4*u1**(-1)*u2**(-1)*us - 4*m1**4*
     &    m2**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_um * ( 4*u*us + 4*m2**2
     &    *s**(-1)*u*us + 4*m2**2*us - 4*m2**4*s**(-1)*us - 4*m1**2*
     &    s**(-1)*u*us + 4*m1**2*m2**2*s**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tx * ( 4*m2**2*s*
     &    t1**(-1)*ts - 3*m2**4*t1**(-1)*ts - m2**4*t2**(-1)*ts + 4*
     &    m2**4*s*t1**(-1)*t2**(-1)*ts - 3*m2**6*t1**(-1)*t2**(-1)*ts
     &     + 2*m1**2*s*t1**(-1)*ts - 2*m1**2*s*t2**(-1)*ts + 6*m1**2*
     &    m2**2*t1**(-1)*ts - 2*m1**2*m2**2*t2**(-1)*ts + 6*m1**2*m2**2
     &    *s*t1**(-1)*t2**(-1)*ts + 9*m1**2*m2**4*t1**(-1)*t2**(-1)*ts
     &     - 3*m1**4*t1**(-1)*ts + 3*m1**4*t2**(-1)*ts - 2*m1**4*s*
     &    t1**(-1)*t2**(-1)*ts - 9*m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 3
     &    *m1**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tm * (  - 4*u*ts - 4*
     &    m2**2*s**(-1)*u*ts - 4*m2**2*ts + 4*m2**4*s**(-1)*ts + 4*
     &    m1**2*s**(-1)*u*ts - 4*m1**2*m2**2*s**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_ux * ( 3*u1**(-1)*us - 3*u2**(-1)*us
     &     + 3*m2**2*u1**(-1)*u2**(-1)*us - 3*m1**2*u1**(-1)*u2**(-1)*
     &    us )
      QBV = QBV + SCB(3,5)*b_um * ( 4*s**(-1)*us )
      QBV = QBV + SCB(3,5)*b_tx * ( 3*t1**(-1)*ts - 3*t2**(-1)*ts
     &     + 3*m2**2*t1**(-1)*t2**(-1)*ts - 3*m1**2*t1**(-1)*t2**(-1)*
     &    ts )
      QBV = QBV + SCB(3,5)*b_u * ( 4*m2**2*u2**(-1) )
      QBV = QBV + SCB(3,5)*b_t * ( 4*m2**2*t2**(-1) )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_ux * (  - 4*m2**2*s*
     &    u1**(-1)*us + 3*m2**4*u1**(-1)*us + m2**4*u2**(-1)*us - 4*
     &    m2**4*s*u1**(-1)*u2**(-1)*us + 3*m2**6*u1**(-1)*u2**(-1)*us
     &     - 4*m1**2*s*u1**(-1)*us - 6*m1**2*m2**2*u1**(-1)*us - 2*
     &    m1**2*m2**2*u2**(-1)*us - 12*m1**2*m2**2*s*u1**(-1)*u2**(-1)*
     &    us - 5*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 3*m1**4*u1**(-1)*us
     &     + m1**4*u2**(-1)*us + m1**4*m2**2*u1**(-1)*u2**(-1)*us + 
     &    m1**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_um * (  - 8*u*us - 4*s*
     &    us + 4*m2**2*us + 4*m1**2*us )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tx * (  - 4*m2**2*s*
     &    t1**(-1)*ts + 3*m2**4*t1**(-1)*ts + m2**4*t2**(-1)*ts - 4*
     &    m2**4*s*t1**(-1)*t2**(-1)*ts + 3*m2**6*t1**(-1)*t2**(-1)*ts
     &     - 4*m1**2*s*t1**(-1)*ts - 6*m1**2*m2**2*t1**(-1)*ts - 2*
     &    m1**2*m2**2*t2**(-1)*ts - 12*m1**2*m2**2*s*t1**(-1)*t2**(-1)*
     &    ts - 5*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 3*m1**4*t1**(-1)*ts
     &     + m1**4*t2**(-1)*ts + m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 
     &    m1**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tm * ( 8*u*ts + 4*s*ts
     &     - 4*m2**2*ts - 4*m1**2*ts )
      QBV = QBV + SCB(7,1)*b_ux * (  - 3*u1**(-1)*us - u2**(-1)*us
     &     - 3*m2**2*u1**(-1)*u2**(-1)*us - m1**2*u1**(-1)*u2**(-1)*us
     &     )
      QBV = QBV + SCB(7,1)*b_tx * (  - 3*t1**(-1)*ts - t2**(-1)*ts
     &     - 3*m2**2*t1**(-1)*t2**(-1)*ts - m1**2*t1**(-1)*t2**(-1)*ts
     &     )
      QBV = QBV + SCB(7,1)*b_s * (  - 6 )
      QBV = QBV + SCC(2,1)*b_tx * ( 3*t1**(-1)*ts - t2**(-1)*ts - 
     &    2*ms**2*t1**(-1)*t2**(-1)*ts + m2**2*t1**(-1)*t2**(-1)*ts + 
     &    m1**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_tm * ( 4*t1**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_t * ( 4*ms**2*t1**(-1) - 4*m1**2*
     &    t1**(-1) )
      QBV = QBV + SCC(2,2)*b_tx * (  - 2*t1**(-1)*ts + 4*t2**(-1)*
     &    ts - 2*ms**2*t1**(-1)*t2**(-1)*ts + 2*m1**2*t1**(-1)*t2**(-1)
     &    *ts )
      QBV = QBV + SCC(2,2)*b_tm * ( 4*t2**(-1)*ts )
      QBV = QBV + SCC(2,2)*b_t * ( 4*ms**2*t2**(-1) - 4*m2**2*
     &    t2**(-1) )
      QBV = QBV + SCC(2,3)*b_ux * ( 5*u1**(-1)*us - 3*u2**(-1)*us
     &     - 2*ms**2*u1**(-1)*u2**(-1)*us + 3*m2**2*u1**(-1)*u2**(-1)*
     &    us - m1**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(2,3)*b_um * ( 4*u1**(-1)*us )
      QBV = QBV + SCC(2,3)*b_u * ( 4*ms**2*u1**(-1) - 4*m1**2*
     &    u1**(-1) )
      QBV = QBV + SCC(2,4)*b_ux * ( 2*u2**(-1)*us - 2*ms**2*
     &    u1**(-1)*u2**(-1)*us + 2*m2**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_um * ( 4*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_u * ( 4*ms**2*u2**(-1) - 4*m2**2*
     &    u2**(-1) )
      QBV = QBV + SCC(4,1)*b_ux * (  - 2*u2**(-1)*us - 2*ms**2*
     &    u1**(-1)*u2**(-1)*us + 2*m2**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(4,1)*b_tx * (  - 2*t1**(-1)*ts - 2*ms**2*
     &    t1**(-1)*t2**(-1)*ts + 2*m1**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCC(4,1)*b_s * (  - 4 )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * (  - u*us + u**2*
     &    u2**(-1)*us - s*us + s*u*u1**(-1)*us + s**2*u1**(-1)*us - 
     &    s**2*u2**(-1)*us - ms**2*u*u1**(-1)*us + ms**2*u*u2**(-1)*us
     &     - 3*ms**2*s*u1**(-1)*us - ms**2*s*u2**(-1)*us - 2*m2**2*us
     &     + m2**2*u*u1**(-1)*us + 3*m2**2*s*u1**(-1)*us - 3*m2**2*s*
     &    u2**(-1)*us + 2*m2**2*ms**2*u1**(-1)*us + m2**2*ms**2*
     &    u2**(-1)*us - 3*m2**2*ms**2*s*u1**(-1)*u2**(-1)*us - 2*m2**4*
     &    u1**(-1)*us + m2**4*u2**(-1)*us + 5*m2**4*s*u1**(-1)*u2**(-1)
     &    *us + 2*m2**4*ms**2*u1**(-1)*u2**(-1)*us - 3*m2**6*u1**(-1)*
     &    u2**(-1)*us - m1**2*us + m1**2*u*u2**(-1)*us - 3*m1**2*s*
     &    u1**(-1)*us + 2*m1**2*s*u2**(-1)*us + 3*m1**2*ms**2*u1**(-1)*
     &    us + 2*m1**2*ms**2*u2**(-1)*us - m1**2*ms**2*s*u1**(-1)*
     &    u2**(-1)*us - 3*m1**2*m2**2*u1**(-1)*us - 7*m1**2*m2**2*
     &    u2**(-1)*us - 9*m1**2*m2**2*s*u1**(-1)*u2**(-1)*us - 4*m1**2*
     &    m2**2*ms**2*u1**(-1)*u2**(-1)*us + 7*m1**2*m2**4*u1**(-1)*
     &    u2**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * ( 2*m1**4*ms**2*
     &    u1**(-1)*u2**(-1)*us - 5*m1**4*m2**2*u1**(-1)*u2**(-1)*us + 
     &    m1**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_um * (  - 4*u*us - 4*s*
     &    us - 8*ms**2*s**(-1)*u*us - 4*ms**2*us + 4*m2**2*s**(-1)*u*us
     &     + 8*m2**2*us + 4*m2**2*ms**2*s**(-1)*us - 4*m2**4*s**(-1)*us
     &     + 4*m1**2*s**(-1)*u*us + 8*m1**2*us + 4*m1**2*ms**2*s**(-1)*
     &    us - 4*m1**4*s**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * (  - u*ts - u**2*
     &    t1**(-1)*ts + 3*s*ts + 2*s*u*t1**(-1)*ts + ms**2*ts + ms**2*u
     &    *t1**(-1)*ts - 2*ms**2*s*t1**(-1)*ts - ms**2*s*t2**(-1)*ts - 
     &    2*m2**2*ts - m2**2*u*t1**(-1)*ts - 3*m2**2*s*t1**(-1)*ts + 4*
     &    m2**2*s*t2**(-1)*ts + m2**2*ms**2*t1**(-1)*ts + 2*m2**2*ms**2
     &    *t2**(-1)*ts - 3*m2**2*ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m2**4
     &    *t1**(-1)*ts - 3*m2**4*t2**(-1)*ts - 4*m2**4*s*t1**(-1)*
     &    t2**(-1)*ts + 2*m2**4*ms**2*t1**(-1)*t2**(-1)*ts + 3*m2**6*
     &    t1**(-1)*t2**(-1)*ts - 2*m1**2*ts - 2*m1**2*u*t1**(-1)*ts + 
     &    m1**2*s*t1**(-1)*ts + 3*m1**2*s*t2**(-1)*ts + 2*m1**2*ms**2*
     &    t1**(-1)*ts + 2*m1**2*ms**2*t2**(-1)*ts - m1**2*ms**2*s*
     &    t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*
     &    m2**2*t2**(-1)*ts - 3*m1**2*m2**2*s*t1**(-1)*t2**(-1)*ts - 4*
     &    m1**2*m2**2*ms**2*t1**(-1)*t2**(-1)*ts - 11*m1**2*m2**4*
     &    t1**(-1)*t2**(-1)*ts + 2*m1**4*t1**(-1)*ts - 5*m1**4*t2**(-1)
     &    *ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * ( 3*m1**4*s*
     &    t1**(-1)*t2**(-1)*ts + 2*m1**4*ms**2*t1**(-1)*t2**(-1)*ts + 
     &    13*m1**4*m2**2*t1**(-1)*t2**(-1)*ts - 5*m1**6*t1**(-1)*
     &    t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tm * ( 4*u*ts + 8*ms**2
     &    *s**(-1)*u*ts + 4*ms**2*ts - 4*m2**2*s**(-1)*u*ts - 4*m2**2*
     &    ms**2*s**(-1)*ts - 4*m1**2*s**(-1)*u*ts - 4*m1**2*ms**2*
     &    s**(-1)*ts + 8*m1**2*m2**2*s**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_ux * (  - 3*s**(-1)*us + 3*s**(-1)*u*
     &    u2**(-1)*us - u1**(-1)*us + 3*u2**(-1)*us + 2*ms**2*s**(-1)*
     &    u1**(-1)*us + 2*ms**2*s**(-1)*u2**(-1)*us + 2*ms**2*u1**(-1)*
     &    u2**(-1)*us - 2*m2**2*s**(-1)*u1**(-1)*us - 3*m2**2*s**(-1)*
     &    u2**(-1)*us - 2*m2**2*u1**(-1)*u2**(-1)*us - 2*m1**2*s**(-1)*
     &    u2**(-1)*us )
      QBV = QBV + SCC(5,1)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCC(5,1)*b_tx * (  - 3*s**(-1)*ts - 2*s**(-1)*u*
     &    t1**(-1)*ts - s**(-1)*u*t2**(-1)*ts - t1**(-1)*ts + 3*
     &    t2**(-1)*ts + 2*ms**2*s**(-1)*t1**(-1)*ts + 2*ms**2*s**(-1)*
     &    t2**(-1)*ts + 2*ms**2*t1**(-1)*t2**(-1)*ts - 4*m2**2*t1**(-1)
     &    *t2**(-1)*ts - m1**2*s**(-1)*t2**(-1)*ts + 2*m1**2*t1**(-1)*
     &    t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCD(1,1)*b_tx * (  - 2 - 2*ms**2*t1**(-1) + 2*
     &    ms**2*t2**(-1) - 2*ms**4*t1**(-1)*t2**(-1) + 2*m1**2*t1**(-1)
     &     - 2*m1**2*t2**(-1) + 4*m1**2*ms**2*t1**(-1)*t2**(-1) - 2*
     &    m1**4*t1**(-1)*t2**(-1) )
      QBV = QBV + SCD(1,1)*b_tm * (  - 4 )
      QBV = QBV + SCD(1,2)*b_ux * (  - 3 - t*u2**(-1) - s*u2**(-1)
     &     + 2*ms**2*u1**(-1) - 2*ms**2*u2**(-1) - 2*ms**4*u1**(-1)*
     &    u2**(-1) - 2*m2**2*u1**(-1) + 2*m2**2*u2**(-1) + 4*m2**2*
     &    ms**2*u1**(-1)*u2**(-1) - 2*m2**4*u1**(-1)*u2**(-1) + m1**2*
     &    u2**(-1) )
      QBV = QBV + SCD(1,2)*b_um * (  - 4 )

      QBS =
     &  + b_s * ( 2 )
      QBS = QBS + SCB(1,1)*b_ux2 * (  - 4*u2**(-1)*us + 2*m2**2*
     &    u1**(-1)*u2**(-1)*us - 2*m2**2*t2**(-1)*u2**(-1)*us + 2*m2**2
     &    *s*u1**(-1)*t2**(-1)*u2**(-1)*us + 2*m1**2*u1**(-1)*u2**(-1)*
     &    us + 2*m1**2*t2**(-1)*u2**(-1)*us + 4*m1**2*s*t1**(-1)*
     &    u1**(-1)*u2**(-1)*us + 2*m1**2*s*u1**(-1)*t2**(-1)*u2**(-1)*
     &    us )
      QBS = QBS + SCB(1,1)*b_um2 * (  - 4*s**(-1)*us )
      QBS = QBS + SCB(1,1)*b_t * ( 4*u*ts**(-1) + 4*s*ts**(-1) + 4
     &    *mg**2*ts**(-1) - 4*m2**2*ts**(-1) - 4*m1**2*ts**(-1) )
      QBS = QBS + SCB(1,1)*Clot*b_t * (  - 4*m1*mg*t1**(-1) )
      QBS = QBS + SCB(1,1)*Cupt*b_t * (  - 4*m2*mg*t2**(-1) )
      QBS = QBS + SCB(1,2)*b_tx4 * (  - 2*t1**(-1)*ts + 2*u1**(-1)
     &    *ts + 2*s*t1**(-1)*u1**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts - 2
     &    *s*u1**(-1)*t2**(-1)*ts - 2*s**2*t1**(-1)*u1**(-1)*t2**(-1)*
     &    ts + 2*m2**2*t1**(-1)*t2**(-1)*ts + 4*m2**2*s*t1**(-1)*
     &    u1**(-1)*t2**(-1)*ts + 4*m2**2*s*t1**(-1)*t2**(-1)*u2**(-1)*
     &    ts + 2*m1**2*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCB(1,2)*b_tm4 * (  - 4*s**(-1)*ts )
      QBS = QBS + SCB(1,2)*b_u * (  - 4*u*us**(-1) + 4*mg**2*
     &    us**(-1) )
      QBS = QBS + SCB(1,2)*Clou*b_u * (  - 4*m2*mg*u2**(-1) )
      QBS = QBS + SCB(1,2)*Cupu*b_u * (  - 4*m1*mg*u1**(-1) )
      QBS = QBS + SCB(1,3)*b_u * ( 4*ms**2*us**(-1) - 4*mg**2*
     &    us**(-1) )
      QBS = QBS + SCB(1,3)*b_t * ( 4*ms**2*ts**(-1) - 4*mg**2*
     &    ts**(-1) )
      QBS = QBS + SCB(3,4)*kaellen**(-1)*b_tx4 * ( 6*m2**2*s*
     &    t1**(-1)*ts - 6*m2**2*s*t2**(-1)*ts + 6*m2**4*s*t1**(-1)*
     &    t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts - 6*m1**2*s*t2**(-1)*ts
     &     - 4*m1**2*m2**2*t2**(-1)*ts - 8*m1**2*m2**2*s*t1**(-1)*
     &    t2**(-1)*ts + 4*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 4*m1**4*
     &    t2**(-1)*ts - 6*m1**4*s*t1**(-1)*t2**(-1)*ts - 8*m1**4*m2**2*
     &    t1**(-1)*t2**(-1)*ts + 4*m1**6*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCB(3,4)*kaellen**(-1)*b_tm4 * ( 4*u*ts + 4*s*ts
     &     - 4*m2**2*s**(-1)*u*ts - 8*m2**2*ts + 4*m2**4*s**(-1)*ts + 4
     &    *m1**2*s**(-1)*u*ts - 4*m1**2*ts - 4*m1**2*m2**2*s**(-1)*ts )
      QBS = QBS + SCB(3,4)*kaellen**(-1)*b_ux2 * ( 4*m2**2*s*
     &    u1**(-1)*us - 4*m2**2*s*u2**(-1)*us - 2*m2**4*u1**(-1)*us + 2
     &    *m2**4*u2**(-1)*us + 4*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**6
     &    *u1**(-1)*u2**(-1)*us + 4*m1**2*s*u1**(-1)*us - 8*m1**2*s*
     &    u2**(-1)*us - 4*m1**2*m2**2*u2**(-1)*us - 4*m1**2*m2**2*s*
     &    u1**(-1)*u2**(-1)*us + 6*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 2
     &    *m1**4*u1**(-1)*us + 2*m1**4*u2**(-1)*us - 8*m1**4*s*u1**(-1)
     &    *u2**(-1)*us - 6*m1**4*m2**2*u1**(-1)*u2**(-1)*us + 2*m1**6*
     &    u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCB(3,4)*kaellen**(-1)*b_um2 * (  - 4*u*us - 4*s
     &    *us + 4*m2**2*s**(-1)*u*us + 8*m2**2*us - 4*m2**4*s**(-1)*us
     &     - 4*m1**2*s**(-1)*u*us + 4*m1**2*us + 4*m1**2*m2**2*s**(-1)*
     &    us )
      QBS = QBS + SCB(3,4)*b_tx4 * ( 2*t1**(-1)*ts - 2*t2**(-1)*ts
     &     + 2*s*t1**(-1)*t2**(-1)*ts + 2*s*u1**(-1)*t2**(-1)*ts + 2*
     &    s**2*t1**(-1)*u1**(-1)*t2**(-1)*ts - 2*m2**2*t1**(-1)*
     &    u1**(-1)*ts - 4*m2**2*s*t1**(-1)*u1**(-1)*t2**(-1)*ts + 2*
     &    m1**2*t1**(-1)*u1**(-1)*ts - 4*m1**2*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCB(3,4)*b_ux2 * (  - 4*m1**2*u1**(-1)*u2**(-1)*
     &    us - 4*m1**2*s*t1**(-1)*u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCB(3,4)*b_um2 * ( 4*s**(-1)*us )
      QBS = QBS + SCB(3,4)*Clot*b_t * ( 4*m1*mg*t1**(-1) )
      QBS = QBS + SCB(3,4)*Cupu*b_u * ( 4*m1*mg*u1**(-1) )
      QBS = QBS + SCB(3,5)*kaellen**(-1)*b_tx4 * (  - 8*m2**2*s*
     &    t1**(-1)*ts + 4*m2**2*s*t2**(-1)*ts + 2*m2**4*t1**(-1)*ts + 2
     &    *m2**4*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*t2**(-1)*ts + 2*m2**6
     &    *t1**(-1)*t2**(-1)*ts - 4*m1**2*s*t1**(-1)*ts + 4*m1**2*s*
     &    t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts - 4*m1**2*m2**2*s*
     &    t1**(-1)*t2**(-1)*ts - 6*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2
     &    *m1**4*t1**(-1)*ts - 2*m1**4*t2**(-1)*ts + 4*m1**4*s*t1**(-1)
     &    *t2**(-1)*ts + 6*m1**4*m2**2*t1**(-1)*t2**(-1)*ts - 2*m1**6*
     &    t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCB(3,5)*kaellen**(-1)*b_tm4 * ( 4*u*ts + 4*
     &    m2**2*s**(-1)*u*ts + 4*m2**2*ts - 4*m2**4*s**(-1)*ts - 4*
     &    m1**2*s**(-1)*u*ts + 4*m1**2*m2**2*s**(-1)*ts )
      QBS = QBS + SCB(3,5)*kaellen**(-1)*b_ux2 * (  - 6*m2**2*s*
     &    u1**(-1)*us + 2*m2**2*s*u2**(-1)*us + 4*m2**4*u1**(-1)*us - 6
     &    *m2**4*s*u1**(-1)*u2**(-1)*us + 4*m2**6*u1**(-1)*u2**(-1)*us
     &     - 6*m1**2*s*u1**(-1)*us + 6*m1**2*s*u2**(-1)*us - 4*m1**2*
     &    m2**2*u1**(-1)*us - 8*m1**2*m2**2*s*u1**(-1)*u2**(-1)*us - 8*
     &    m1**2*m2**4*u1**(-1)*u2**(-1)*us + 6*m1**4*s*u1**(-1)*
     &    u2**(-1)*us + 4*m1**4*m2**2*u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCB(3,5)*kaellen**(-1)*b_um2 * (  - 4*u*us - 4*
     &    m2**2*s**(-1)*u*us - 4*m2**2*us + 4*m2**4*s**(-1)*us + 4*
     &    m1**2*s**(-1)*u*us - 4*m1**2*m2**2*s**(-1)*us )
      QBS = QBS + SCB(3,5)*b_tx4 * (  - 4*m2**2*t1**(-1)*t2**(-1)*
     &    ts - 4*m2**2*s*t1**(-1)*t2**(-1)*u2**(-1)*ts )
      QBS = QBS + SCB(3,5)*b_tm4 * ( 4*s**(-1)*ts )
      QBS = QBS + SCB(3,5)*b_ux2 * (  - 2*u1**(-1)*us + 2*u2**(-1)
     &    *us - 4*m2**2*u1**(-1)*u2**(-1)*us + 2*m2**2*t2**(-1)*
     &    u2**(-1)*us - 2*m2**2*s*u1**(-1)*t2**(-1)*u2**(-1)*us - 2*
     &    m1**2*t2**(-1)*u2**(-1)*us - 2*m1**2*s*u1**(-1)*t2**(-1)*
     &    u2**(-1)*us )
      QBS = QBS + SCB(3,5)*Clou*b_u * ( 4*m2*mg*u2**(-1) )
      QBS = QBS + SCB(3,5)*Cupt*b_t * ( 4*m2*mg*t2**(-1) )
      QBS = QBS + SCB(4,1)*kaellen**(-1)*b_tx4 * ( 2*m2**2*s*
     &    t1**(-1)*ts + 2*m2**2*s*t2**(-1)*ts - 2*m2**4*t1**(-1)*ts - 2
     &    *m2**4*t2**(-1)*ts + 2*m2**4*s*t1**(-1)*t2**(-1)*ts - 2*m2**6
     &    *t1**(-1)*t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts + 2*m1**2*s*
     &    t2**(-1)*ts + 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*m2**2*
     &    t2**(-1)*ts + 12*m1**2*m2**2*s*t1**(-1)*t2**(-1)*ts + 2*m1**2
     &    *m2**4*t1**(-1)*t2**(-1)*ts - 2*m1**4*t1**(-1)*ts - 2*m1**4*
     &    t2**(-1)*ts + 2*m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*m2**2*
     &    t1**(-1)*t2**(-1)*ts - 2*m1**6*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCB(4,1)*kaellen**(-1)*b_tm4 * (  - 8*u*ts - 4*s
     &    *ts + 4*m2**2*ts + 4*m1**2*ts )
      QBS = QBS + SCB(4,1)*kaellen**(-1)*b_ux2 * ( 2*m2**2*s*
     &    u1**(-1)*us + 2*m2**2*s*u2**(-1)*us - 2*m2**4*u1**(-1)*us - 2
     &    *m2**4*u2**(-1)*us + 2*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**6
     &    *u1**(-1)*u2**(-1)*us + 2*m1**2*s*u1**(-1)*us + 2*m1**2*s*
     &    u2**(-1)*us + 4*m1**2*m2**2*u1**(-1)*us + 4*m1**2*m2**2*
     &    u2**(-1)*us + 12*m1**2*m2**2*s*u1**(-1)*u2**(-1)*us + 2*m1**2
     &    *m2**4*u1**(-1)*u2**(-1)*us - 2*m1**4*u1**(-1)*us - 2*m1**4*
     &    u2**(-1)*us + 2*m1**4*s*u1**(-1)*u2**(-1)*us + 2*m1**4*m2**2*
     &    u1**(-1)*u2**(-1)*us - 2*m1**6*u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCB(4,1)*kaellen**(-1)*b_um2 * ( 8*u*us + 4*s*us
     &     - 4*m2**2*us - 4*m1**2*us )
      QBS = QBS + SCB(4,1)*b_tx4 * ( 2*t1**(-1)*ts + 2*t2**(-1)*ts
     &     + 2*m2**2*t1**(-1)*t2**(-1)*ts + 2*m1**2*t1**(-1)*t2**(-1)*
     &    ts )
      QBS = QBS + SCB(4,1)*b_ux2 * ( 2*u1**(-1)*us + 2*u2**(-1)*us
     &     + 2*m2**2*u1**(-1)*u2**(-1)*us + 2*m1**2*u1**(-1)*u2**(-1)*
     &    us )
      QBS = QBS + SCB(4,1)*b_s * ( 2 - 4*ms**2*s**(-1) + 4*mg**2*
     &    s**(-1) )
      QBS = QBS + SCB(6,1)*b_tx4 * (  - 2*t1**(-1)*ts - 2*u1**(-1)
     &    *ts - 2*s*t1**(-1)*u1**(-1)*ts + 2*m2**2*t1**(-1)*u1**(-1)*ts
     &     - 2*m1**2*t1**(-1)*u1**(-1)*ts )
      QBS = QBS + SCB(6,1)*b_u * (  - 2 )
      QBS = QBS + SCB(6,1)*b_t * (  - 2 )
      QBS = QBS + SCB(6,1)*b_s * (  - 2 + 4*ms**2*s**(-1) - 4*
     &    mg**2*s**(-1) )
      QBS = QBS + SCBP(1)*b_u * ( 2*ms**2 - 2*mg**2 )
      QBS = QBS + SCBP(1)*b_t * ( 2*ms**2 - 2*mg**2 )
      QBS = QBS + SCBP(1)*b_s * ( 2*ms**2 - 2*mg**2 )
      QBS = QBS + SCC(1,1)*b_ux2 * ( 2*u1**(-1)*us + 2*s*u1**(-1)*
     &    u2**(-1)*us - 2*ms**2*s**(-1)*u1**(-1)*us - 2*ms**2*s**(-1)*
     &    u2**(-1)*us - 4*ms**2*u1**(-1)*u2**(-1)*us + 2*mg**2*s**(-1)*
     &    u2**(-1)*us + 2*mg**2*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     &    s**(-1)*u1**(-1)*u2**(-1)*us - 2*m2**2*mg**2*s**(-1)*u1**(-1)
     &    *u2**(-1)*us + 2*m1**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 4
     &    *m1**2*ms**2*t1**(-1)*u1**(-1)*u2**(-1)*us - 2*m1**2*mg**2*
     &    s**(-1)*u1**(-1)*u2**(-1)*us - 4*m1**2*mg**2*t1**(-1)*
     &    u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCC(1,1)*b_um2 * (  - 2*s**(-2)*t*us - 4*ms**2*
     &    s**(-2)*us + 4*mg**2*s**(-2)*us + 2*m1**2*s**(-2)*us )
      QBS = QBS + SCC(1,1)*Clot*b_t * (  - 4*m1*mg*s**(-1) - 4*m1*
     &    mg*ms**2*s**(-1)*t1**(-1) + 4*m1*mg**3*s**(-1)*t1**(-1) )
      QBS = QBS + SCC(1,2)*b_ux2 * (  - 2*s**(-1)*us - 2*s**(-1)*t
     &    *u2**(-1)*us + 2*s*u1**(-1)*u2**(-1)*us - 4*ms**2*s**(-1)*
     &    u2**(-1)*us - 4*ms**2*u1**(-1)*u2**(-1)*us + 2*mg**2*s**(-1)*
     &    u2**(-1)*us + 2*mg**2*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     &    s**(-1)*u1**(-1)*u2**(-1)*us - 2*m2**2*ms**2*s**(-1)*t2**(-1)
     &    *u2**(-1)*us + 2*m2**2*ms**2*u1**(-1)*t2**(-1)*u2**(-1)*us - 
     &    2*m2**2*mg**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*mg**2*
     &    s**(-1)*t2**(-1)*u2**(-1)*us - 2*m2**2*mg**2*u1**(-1)*
     &    t2**(-1)*u2**(-1)*us + 2*m1**2*s**(-1)*u2**(-1)*us + 2*m1**2*
     &    ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m1**2*ms**2*s**(-1)*
     &    t2**(-1)*u2**(-1)*us + 2*m1**2*ms**2*u1**(-1)*t2**(-1)*
     &    u2**(-1)*us - 2*m1**2*mg**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*
     &    m1**2*mg**2*s**(-1)*t2**(-1)*u2**(-1)*us - 2*m1**2*mg**2*
     &    u1**(-1)*t2**(-1)*u2**(-1)*us )
      QBS = QBS + SCC(1,2)*b_um2 * (  - 2*s**(-2)*t*us - 4*ms**2*
     &    s**(-2)*us + 4*mg**2*s**(-2)*us + 2*m2**2*s**(-2)*us )
      QBS = QBS + SCC(1,2)*Cupt*b_t * (  - 4*m2*mg*s**(-1) - 4*m2*
     &    mg*ms**2*s**(-1)*t2**(-1) + 4*m2*mg**3*s**(-1)*t2**(-1) )
      QBS = QBS + SCC(1,3)*b_tx4 * (  - 2*s**(-1)*ts - 2*s**(-1)*u
     &    *t1**(-1)*ts + 2*s*t1**(-1)*t2**(-1)*ts - 2*ms**2*s**(-1)*
     &    t1**(-1)*ts + 2*ms**2*s**(-1)*u1**(-1)*ts + 2*ms**2*t1**(-1)*
     &    u1**(-1)*ts - 6*ms**2*t1**(-1)*t2**(-1)*ts - 2*ms**2*u1**(-1)
     &    *t2**(-1)*ts - 2*ms**2*s*t1**(-1)*u1**(-1)*t2**(-1)*ts + 2*
     &    mg**2*s**(-1)*t1**(-1)*ts + 4*mg**2*t1**(-1)*t2**(-1)*ts + 2*
     &    mg**2*u1**(-1)*t2**(-1)*ts + 2*mg**2*s*t1**(-1)*u1**(-1)*
     &    t2**(-1)*ts + 2*m2**2*s**(-1)*t1**(-1)*ts + 2*m2**2*ms**2*
     &    s**(-1)*t1**(-1)*t2**(-1)*ts + 4*m2**2*ms**2*t1**(-1)*
     &    u1**(-1)*t2**(-1)*ts - 2*m2**2*mg**2*s**(-1)*t1**(-1)*
     &    u1**(-1)*ts - 2*m2**2*mg**2*s**(-1)*t1**(-1)*t2**(-1)*ts - 4*
     &    m2**2*mg**2*t1**(-1)*u1**(-1)*t2**(-1)*ts + 2*m1**2*s**(-1)*
     &    t1**(-1)*ts + 2*m1**2*s**(-1)*u1**(-1)*ts + 2*m1**2*t1**(-1)*
     &    u1**(-1)*ts + 2*m1**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 2*
     &    m1**2*mg**2*s**(-1)*t1**(-1)*u1**(-1)*ts - 2*m1**2*mg**2*
     &    s**(-1)*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCC(1,3)*b_tx4 * (  - 2*m1**2*m2**2*s**(-1)*
     &    t1**(-1)*u1**(-1)*ts + 2*m1**4*s**(-1)*t1**(-1)*u1**(-1)*ts )
      QBS = QBS + SCC(1,3)*b_tm4 * (  - 2*s**(-2)*u*ts - 4*ms**2*
     &    s**(-2)*ts + 4*mg**2*s**(-2)*ts + 2*m1**2*s**(-2)*ts )
      QBS = QBS + SCC(1,3)*Cupu*b_u * (  - 4*m1*mg*s**(-1) - 4*m1*
     &    mg*ms**2*s**(-1)*u1**(-1) + 4*m1*mg**3*s**(-1)*u1**(-1) )
      QBS = QBS + SCC(1,4)*b_tx4 * ( 2*t2**(-1)*ts + 2*s*t1**(-1)*
     &    t2**(-1)*ts - 2*ms**2*s**(-1)*t1**(-1)*ts - 2*ms**2*s**(-1)*
     &    t2**(-1)*ts - 4*ms**2*t1**(-1)*t2**(-1)*ts + 2*mg**2*s**(-1)*
     &    t1**(-1)*ts + 2*mg**2*t1**(-1)*t2**(-1)*ts + 2*m2**2*ms**2*
     &    s**(-1)*t1**(-1)*t2**(-1)*ts + 4*m2**2*ms**2*t1**(-1)*
     &    t2**(-1)*u2**(-1)*ts - 2*m2**2*mg**2*s**(-1)*t1**(-1)*
     &    t2**(-1)*ts - 4*m2**2*mg**2*t1**(-1)*t2**(-1)*u2**(-1)*ts + 2
     &    *m1**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts - 2*m1**2*mg**2*
     &    s**(-1)*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCC(1,4)*b_tm4 * (  - 2*s**(-2)*u*ts - 4*ms**2*
     &    s**(-2)*ts + 4*mg**2*s**(-2)*ts + 2*m2**2*s**(-2)*ts )
      QBS = QBS + SCC(1,4)*Clou*b_u * (  - 4*m2*mg*s**(-1) - 4*m2*
     &    mg*ms**2*s**(-1)*u2**(-1) + 4*m2*mg**3*s**(-1)*u2**(-1) )
      QBS = QBS + SCC(3,1)*b_tx4 * (  - 2*s*t1**(-1)*t2**(-1)*ts
     &     + 4*ms**2*t1**(-1)*t2**(-1)*ts - 2*mg**2*t1**(-1)*t2**(-1)*
     &    ts )
      QBS = QBS + SCC(3,1)*b_tm4 * (  - 2*s**(-1)*ts )
      QBS = QBS + SCC(3,1)*b_ux2 * (  - 2*s*u1**(-1)*u2**(-1)*us
     &     + 4*ms**2*u1**(-1)*u2**(-1)*us - 2*mg**2*u1**(-1)*u2**(-1)*
     &    us )
      QBS = QBS + SCC(3,1)*b_um2 * (  - 2*s**(-1)*us )
      QBS = QBS + SCC(3,1)*b_s * ( 4*ms**4*s**(-2) + 4*mg**2*
     &    s**(-1) - 8*mg**2*ms**2*s**(-2) + 4*mg**4*s**(-2) )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_tx4 * ( 2*u*ts + 2*u**2
     &    *t2**(-1)*ts - 4*s*ts - 2*s*u*t2**(-1)*ts - 2*ms**2*s*
     &    t1**(-1)*ts - 2*ms**2*s*t2**(-1)*ts - 4*m2**2*ts - 4*m2**2*u*
     &    t2**(-1)*ts - 8*m2**2*s*t1**(-1)*ts - 4*m2**2*s*t2**(-1)*ts
     &     + 2*m2**2*ms**2*t1**(-1)*ts + 2*m2**2*ms**2*t2**(-1)*ts - 2*
     &    m2**2*ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m2**4*t1**(-1)*ts + 2*
     &    m2**4*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*t2**(-1)*ts + 2*m2**4*
     &    ms**2*t1**(-1)*t2**(-1)*ts + 2*m2**6*t1**(-1)*t2**(-1)*ts + 2
     &    *m1**2*ts - 6*m1**2*s*t1**(-1)*ts + 4*m1**2*s*t2**(-1)*ts + 2
     &    *m1**2*ms**2*t1**(-1)*ts + 2*m1**2*ms**2*t2**(-1)*ts - 2*
     &    m1**2*ms**2*s*t1**(-1)*t2**(-1)*ts - 2*m1**2*m2**2*t1**(-1)*
     &    ts + 6*m1**2*m2**2*t2**(-1)*ts + 6*m1**2*m2**2*s*t1**(-1)*
     &    t2**(-1)*ts - 4*m1**2*m2**2*ms**2*t1**(-1)*t2**(-1)*ts - 8*
     &    m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 4*m1**4*t1**(-1)*ts - 2*
     &    m1**4*t2**(-1)*ts + 6*m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*
     &    ms**2*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_tx4 * ( 10*m1**4*m2**2*
     &    t1**(-1)*t2**(-1)*ts - 4*m1**6*t1**(-1)*t2**(-1)*ts )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_tm4 * (  - 4*u*ts - 4*s
     &    *ts + 8*ms**2*s**(-1)*u*ts + 4*ms**2*ts + 4*m2**2*s**(-1)*u*
     &    ts + 8*m2**2*ts - 4*m2**2*ms**2*s**(-1)*ts - 4*m2**4*s**(-1)*
     &    ts + 4*m1**2*s**(-1)*u*ts + 8*m1**2*ts - 4*m1**2*ms**2*
     &    s**(-1)*ts - 4*m1**4*s**(-1)*ts )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 2*u*us + 2*
     &    u**2*u1**(-1)*us - 6*s*us + 6*s*u*u1**(-1)*us + 4*s**2*
     &    u1**(-1)*us - 2*ms**2*s*u1**(-1)*us - 2*ms**2*s*u2**(-1)*us
     &     + 4*m2**2*us - 4*m2**2*u*u1**(-1)*us - 2*m2**2*s*u1**(-1)*us
     &     - 6*m2**2*s*u2**(-1)*us + 2*m2**2*ms**2*u1**(-1)*us + 2*
     &    m2**2*ms**2*u2**(-1)*us - 2*m2**2*ms**2*s*u1**(-1)*u2**(-1)*
     &    us + 4*m2**4*u2**(-1)*us + 6*m2**4*s*u1**(-1)*u2**(-1)*us + 2
     &    *m2**4*ms**2*u1**(-1)*u2**(-1)*us - 4*m2**6*u1**(-1)*u2**(-1)
     &    *us - 2*m1**2*us - 6*m1**2*s*u1**(-1)*us - 8*m1**2*s*u2**(-1)
     &    *us + 2*m1**2*ms**2*u1**(-1)*us + 2*m1**2*ms**2*u2**(-1)*us
     &     - 2*m1**2*ms**2*s*u1**(-1)*u2**(-1)*us + 6*m1**2*m2**2*
     &    u1**(-1)*us - 2*m1**2*m2**2*u2**(-1)*us + 6*m1**2*m2**2*s*
     &    u1**(-1)*u2**(-1)*us - 4*m1**2*m2**2*ms**2*u1**(-1)*u2**(-1)*
     &    us + 10*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 2*m1**4*u2**(-1)*
     &    us - 8*m1**4*s*u1**(-1)*u2**(-1)*us + 2*m1**4*ms**2*u1**(-1)*
     &    u2**(-1)*us )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 8*m1**4*
     &    m2**2*u1**(-1)*u2**(-1)*us + 2*m1**6*u1**(-1)*u2**(-1)*us )
      QBS = QBS + SCC(6,1)*kaellen**(-1)*b_um2 * ( 4*u*us - 8*
     &    ms**2*s**(-1)*u*us - 4*ms**2*us - 4*m2**2*s**(-1)*u*us + 4*
     &    m2**2*ms**2*s**(-1)*us - 4*m1**2*s**(-1)*u*us + 4*m1**2*ms**2
     &    *s**(-1)*us + 8*m1**2*m2**2*s**(-1)*us )
      QBS = QBS + SCC(6,1)*b_tx4 * (  - 2*t1**(-1)*ts - 2*t2**(-1)
     &    *ts - 2*s*t1**(-1)*t2**(-1)*ts + 2*ms**2*s**(-1)*t1**(-1)*ts
     &     + 2*ms**2*s**(-1)*t2**(-1)*ts + 4*ms**2*t1**(-1)*t2**(-1)*ts
     &     - 2*mg**2*t1**(-1)*t2**(-1)*ts - 2*m2**2*ms**2*s**(-1)*
     &    t1**(-1)*t2**(-1)*ts + 2*m2**2*mg**2*s**(-1)*t1**(-1)*
     &    t2**(-1)*ts + 4*m1**2*t1**(-1)*t2**(-1)*ts - 2*m1**2*ms**2*
     &    s**(-1)*t1**(-1)*t2**(-1)*ts + 2*m1**2*mg**2*s**(-1)*t1**(-1)
     &    *t2**(-1)*ts )
      QBS = QBS + SCC(6,1)*b_tm4 * ( 4*s**(-2)*u*ts + 2*s**(-1)*ts
     &     + 4*ms**2*s**(-2)*ts - 4*mg**2*s**(-2)*ts - 2*m2**2*s**(-2)*
     &    ts - 2*m1**2*s**(-2)*ts )
      QBS = QBS + SCC(6,1)*b_ux2 * (  - 2*u1**(-1)*us - 2*u2**(-1)
     &    *us - 2*s*u1**(-1)*u2**(-1)*us + 2*ms**2*s**(-1)*u1**(-1)*us
     &     + 2*ms**2*s**(-1)*u2**(-1)*us + 4*ms**2*u1**(-1)*u2**(-1)*us
     &     - 2*mg**2*u1**(-1)*u2**(-1)*us + 4*m2**2*u1**(-1)*u2**(-1)*
     &    us - 2*m2**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*
     &    mg**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*m1**2*ms**2*s**(-1)*
     &    u1**(-1)*u2**(-1)*us + 2*m1**2*mg**2*s**(-1)*u1**(-1)*
     &    u2**(-1)*us )
      QBS = QBS + SCC(6,1)*b_um2 * (  - 4*s**(-2)*u*us - 2*s**(-1)
     &    *us + 4*ms**2*s**(-2)*us - 4*mg**2*s**(-2)*us + 2*m2**2*
     &    s**(-2)*us + 2*m1**2*s**(-2)*us )
      QBS = QBS + SCD(2,1)*b_ux2 * (  - 2*s*u1**(-1)*tg**(-1)*us
     &     - 2*s**2*u1**(-1)*tg**(-1)*u2**(-1)*us + 6*ms**2*u1**(-1)*
     &    tg**(-1)*us + 2*ms**2*tg**(-1)*u2**(-1)*us + 8*ms**2*s*
     &    u1**(-1)*tg**(-1)*u2**(-1)*us - 2*ms**4*s**(-1)*u1**(-1)*
     &    tg**(-1)*us - 2*ms**4*s**(-1)*tg**(-1)*u2**(-1)*us - 8*ms**4*
     &    u1**(-1)*tg**(-1)*u2**(-1)*us - 2*mg**2*u1**(-1)*tg**(-1)*us
     &     - 4*mg**2*s*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*mg**2*ms**2*
     &    s**(-1)*u1**(-1)*tg**(-1)*us + 2*mg**2*ms**2*s**(-1)*tg**(-1)
     &    *u2**(-1)*us + 8*mg**2*ms**2*u1**(-1)*tg**(-1)*u2**(-1)*us - 
     &    2*mg**4*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*ms**4*s**(-1)
     &    *u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*mg**2*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us - 4*m2**2*mg**2*ms**2*s**(-1)*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us + 2*m2**2*mg**4*s**(-1)*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us + 2*m1**2*s*u1**(-1)*tg**(-1)*u2**(-1)*
     &    us - 4*m1**2*ms**2*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m1**2*
     &    ms**4*s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us )
      QBS = QBS + SCD(2,1)*b_ux2 * ( 4*m1**2*mg**2*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us - 4*m1**2*mg**2*ms**2*s**(-1)*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us + 2*m1**2*mg**4*s**(-1)*u1**(-1)*
     &    tg**(-1)*u2**(-1)*us )
      QBS = QBS + SCD(2,1)*b_um2 * ( 2*s**(-1)*us - 4*ms**2*
     &    s**(-2)*us - 4*ms**4*s**(-2)*tg**(-1)*us + 4*mg**2*s**(-2)*us
     &     + 4*mg**2*ms**2*s**(-2)*tg**(-1)*us + 2*m2**2*ms**2*s**(-2)*
     &    tg**(-1)*us - 2*m2**2*mg**2*s**(-2)*tg**(-1)*us + 2*m1**2*
     &    ms**2*s**(-2)*tg**(-1)*us - 2*m1**2*mg**2*s**(-2)*tg**(-1)*us
     &     )
      QBS = QBS + SCD(2,2)*b_tx4 * (  - 2*s*ug**(-1)*t2**(-1)*ts
     &     - 2*s**2*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*ms**2*t1**(-1)*
     &    ug**(-1)*ts + 6*ms**2*ug**(-1)*t2**(-1)*ts + 8*ms**2*s*
     &    t1**(-1)*ug**(-1)*t2**(-1)*ts - 2*ms**4*s**(-1)*t1**(-1)*
     &    ug**(-1)*ts - 2*ms**4*s**(-1)*ug**(-1)*t2**(-1)*ts - 8*ms**4*
     &    t1**(-1)*ug**(-1)*t2**(-1)*ts - 2*mg**2*ug**(-1)*t2**(-1)*ts
     &     - 4*mg**2*s*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*mg**2*ms**2*
     &    s**(-1)*t1**(-1)*ug**(-1)*ts + 2*mg**2*ms**2*s**(-1)*ug**(-1)
     &    *t2**(-1)*ts + 8*mg**2*ms**2*t1**(-1)*ug**(-1)*t2**(-1)*ts - 
     &    2*mg**4*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m2**2*s*t1**(-1)*
     &    ug**(-1)*t2**(-1)*ts - 4*m2**2*ms**2*t1**(-1)*ug**(-1)*
     &    t2**(-1)*ts + 2*m2**2*ms**4*s**(-1)*t1**(-1)*ug**(-1)*
     &    t2**(-1)*ts + 4*m2**2*mg**2*t1**(-1)*ug**(-1)*t2**(-1)*ts - 4
     &    *m2**2*mg**2*ms**2*s**(-1)*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*
     &    m2**2*mg**4*s**(-1)*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*
     &    ms**4*s**(-1)*t1**(-1)*ug**(-1)*t2**(-1)*ts )
      QBS = QBS + SCD(2,2)*b_tx4 * ( 2*m1**2*mg**2*t1**(-1)*
     &    ug**(-1)*t2**(-1)*ts - 4*m1**2*mg**2*ms**2*s**(-1)*t1**(-1)*
     &    ug**(-1)*t2**(-1)*ts + 2*m1**2*mg**4*s**(-1)*t1**(-1)*
     &    ug**(-1)*t2**(-1)*ts )
      QBS = QBS + SCD(2,2)*b_tm4 * ( 2*s**(-1)*ts - 4*ms**2*
     &    s**(-2)*ts - 4*ms**4*s**(-2)*ug**(-1)*ts + 4*mg**2*s**(-2)*ts
     &     + 4*mg**2*ms**2*s**(-2)*ug**(-1)*ts + 2*m2**2*ms**2*s**(-2)*
     &    ug**(-1)*ts - 2*m2**2*mg**2*s**(-2)*ug**(-1)*ts + 2*m1**2*
     &    ms**2*s**(-2)*ug**(-1)*ts - 2*m1**2*mg**2*s**(-2)*ug**(-1)*ts
     &     )

c               function has be real at the end !!!
c               normalization factor !!!
      NN_QBV = real(QBV+QBS) / 8.D0

c               the prefactors removed in the form program 
c                   except for alphas, which is cut 
      NN_QBV = NN_QBV * Cf /Pi

c               the averaging factors
      NN_QBV = NN_QBV /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NN_QBV = NN_QBV * (abs(m1)+abs(m2))**2/4.D0 /4.D0/Pi

      end


