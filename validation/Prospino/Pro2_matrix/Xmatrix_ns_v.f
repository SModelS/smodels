cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NS_QGB(MASSIN,C)                                                 c
c     NS_QGV(MASSIN,C)                                                 c
c     NS_QGD(MASSIN,C)                                                 c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = tg                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = m2                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(20) = del                                               c
c       MASSIN(3-5,8,10,13-19,21-30) not needed                        c
c                                                                      c
c       CL(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CR(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CV(1:4)  FRACTIONS OF COUPLINGS NEEDED FOR VIRTUAL [COMPLEX]   c
c                -> CLOT,CUPT,CUPU,CLOU                                c
c                                                                      c
c                                                                      c
c    ALL PHASE SPACE FACTORS INCOUDED TO GIVE                          c
c        s^2 d sig/(dtg ds4)                                           c
c                                                                      c
c    NEEDED MANDELSTAM VARIABLES :                                     c
c                                                                      c
c       Q(K1) + G(K2) -> SQ(P1) + NJ(P2) [+G(K3)]                      c
c                                                                      c
c       S  = 2(K1.K2)                                                  c
c       S3 = 2(K3.P2)                                                  c
c       S4 = 2(K3.K1)                                                  c
c       S5 = 2(P1.P2) + M1^2 + M2^2                                    c
c       T1 = 2(K1.P1)                                                  c
c       U1 = 2(K2.P1)                                                  c
c       TG = 2(K2.P2)                                                  c
c       UG = 2(K1.P2)                                                  c
c       TP = 2(K2.K3)                                                  c
c       UP = 2(K1.K3)                                                  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
c note that m1 (external squark) and ms (u-channel squark) should 
c be the same, because of the QCD ssg coupling for the t channel 
c 
      real*8 function NS_QGB(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),pi,Nc,Cf
     &          ,s,m1,m2,ms,t1,u1,t2,u2,us,gs
      complex*16 Ar(4),Arc(4),MMborn_s,MMborn_u

      pi= 4.D0 * atan(1.D0)

      Nc = 3.D0 
      Cf = 4.D0/3.D0

      s  = massin(1)
      t2 = massin(2)
      m1 = massin(6)
      m2 = massin(7)
      ms = massin(6)
ctp      ms = massin(11)

      t1 = t2 + m2**2 - m1**2

      u1 = - s - t2
      u2 = u1 + m1**2 - m2**2 
      us = u1 + m1**2 - ms**2 

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

      gs = 1.D0

c               my form output 
c  MMborn_x = 2*M_x*M_born => divide by factor 2 after adding
c    this expression from my FORM output
c    which in turn requires my and not Wim's coupling structures...
      MMborn_s = 0.D0
     +
      MMborn_s = MMborn_s + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     - 8*s**(-1)*t2
     +     + 4*s*us**(-1)
     +     )

      MMborn_u = 0.D0
     +
      MMborn_u = MMborn_u + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     + 16*m1**2*u2*us**(-2)
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     + 4*s*us**(-1)
     +     + 4*u1*u2*us**(-2)
     +     - 4*u2*us**(-1)
     +     )

      NS_QGB = real( MMborn_s + MMborn_u )/2.D0

c               the phase space except for 1/s**2 
      NS_QGB = NS_QGB / ( 16.D0 * pi )

c               the averaging factors
      NS_QGB = NS_QGB /4.D0 /Nc /(Nc**2-1.D0)

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NS_QGB = NS_QGB * (abs(m1)+abs(m2))**2/4.D0

      end


c --------------------------------------------------------------------
c extract the log(Delta) terms from the virtual+soft scaling function 
c everything else copied from NS_QGV 
      real*8 function NS_QGD(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,Nc,CF,CA,Sn,gs
     &          ,six,born,QF,Del,log_del,log2_del
     &          ,s,m1,m2,ms,t1,u1,t2,u2,tj,uj,us,gqFORT,gqFORT1,s4,s4p
      complex*16 Ar(4),Arc(4),MMborn_s,MMborn_u

      Pi    = 4.D0*atan(1.D0)
      Nc    = 3.D0
      CF    = 4.D0/3.D0 
      CA    = 3.D0
      Sn    = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      six = 6.D0

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      tj  = massin(2)
      s4  = massin(3)
      m1  = massin(6)
      m2  = massin(7)
      ms  = massin(11)
      Del = massin(20)
      s4p = massin(21)

c               the logaritms for linear s4 integration 
      log_del  =  log(s4p/m1**2)/(s4p-Del) - 1.D0/s4
      log2_del =  log(s4p/m1**2)**2/(s4p-Del) 
     &          - 2.D0*log(s4/m1**2)/s4

c               born kinematics built in
      t1 = tj + m2**2 - m1**2

      u1 = - s - tj
      uj = u1 + m1**2 - m2**2 

c               oly for my born term 
      us = u1 + m1**2 - ms**2 

c               all conventions twice ?!?
      t2 = tj
      u2 = uj

c               the factorization/renormalization scale 
      QF = massin(13) 

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               this is a copy from the born term 
      MMborn_s = 0.D0
     +
      MMborn_s = MMborn_s + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     - 8*s**(-1)*t2
     +     + 4*s*us**(-1)
     +     )

      MMborn_u = 0.D0
     +
      MMborn_u = MMborn_u + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     + 16*m1**2*u2*us**(-2)
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     + 4*s*us**(-1)
     +     + 4*u1*u2*us**(-2)
     +     - 4*u2*us**(-1)
     +     )
      
      born = real( MMborn_s + MMborn_u )/2.D0

c               'grep cdel' and replace log(Delta) terms 
      NS_QGD = 0.D0

      end




c --------------------------------------------------------------------
      real*8 function NS_QGV(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,Nc,CF,CA,Sn,nf,zeta2,gs
     &          ,six,born,rootlam,QF,QR
     &          ,s,m1,m2,mg,ms,mt,mgs2,mjs2,us,tj,uj,t1,t2,u1,u2
     &          ,gqFORT1,gqFORT2,gqFORT
     &          ,SCA(1:10),SCB(1:10,1:5),SCBP(10)
     &          ,SCC(1:20,1:5),SCD(1:10,1:2)
     &          ,Li2,set1fac,set2fac,boxfac
      complex*16 CSPEN,Ar(4),Arc(4),MMborn_s,MMborn_u

      external CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

      Pi    = 4.D0*atan(1.D0)
      zeta2 = Pi**2/6.D0
      nf    = 6.D0
      Nc    = 3.D0
      CF    = 4.D0/3.D0 
      CA    = 3.D0
      Sn    = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      six = 6.D0

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      tj  = massin(2)
      m1  = massin(6)
      m2  = massin(7)
      mt  = massin(9)
      mg  = massin(10)
      ms  = massin(11)

      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2

c               born kinematics built in
      t1 = tj + m2**2 - m1**2

      u1 = - s - tj
      uj = u1 + m1**2 - m2**2 

c               oly for my born term 
      us = u1 + m1**2 - ms**2 

c               all conventions twice ?!?
      t2 = tj
      u2 = uj

      rootlam = sqrt( s**2 +m1**4 +m2**4
     &           - 2*( s*m1**2 + s*m2**2 + m1**2*m2**2 ))

c               the factorization/renormalization scale 
      QR = massin(12)
      QF = massin(13) 

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

c               the scalar functions 
      call SCALAR_ARRAY_NS(massin,SCA,SCB,SCBP,SCC,SCD)

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               this is a copy from the born term 
      MMborn_s = 0.D0
     +
      MMborn_s = MMborn_s + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     - 8*s**(-1)*t2
     +     + 4*s*us**(-1)
     +     )

      MMborn_u = 0.D0
     +
      MMborn_u = MMborn_u + Ar(1)*Arc(1)*Nc*Cf*gs**2 * (
     +     + 16*m1**2*u2*us**(-2)
     +     - 4*m1**2*us**(-1)
     +     - 4*m2**2*us**(-1)
     +     - 4*s**(-1)*t1*t2*us**(-1)
     +     - 8*s**(-1)*t1*u2*us**(-1)
     +     - 4*s**(-1)*u1*u2*us**(-1)
     +     + 4*s*us**(-1)
     +     + 4*u1*u2*us**(-2)
     +     - 4*u2*us**(-1)
     +     )
      
      born = real( MMborn_s + MMborn_u )/2.D0

c               from wim's notes, but with CC = AA/2 
      set1fac = Nc*Cf * gs**2 * real( Ar(1)*Arc(1) )/2.D0
      set2fac = Nc*Cf * gs**2 * real( Ar(1)*Ar(2) )/2.D0 * m2*mg
      boxfac  = t1*u1-mgs2*tj

c               insert the form output 
c
c   commented using cdel: scale dependence
c   commented using csof: finite soft contributions 
c
c   change in the phase space logs etc: mi -> mx
c
c   change the scalar integrals : cc means not needed here, copied from N-gl
c
cc   A0fin(mt)                   -> SCA(1)
c    A0fin(ms)                   -> SCA(2)
c    A0fin(mg)                   -> SCA(3)
c
cc   B0fin(p1 - g1,ms,0)         -> SCB(3,1)
cc   B0fin(p1 - g1,mg,0)         -> SCB(1,1)
c    B0fin(p1 - g2,ms,0)         -> SCB(3,2)
c    B0fin(p1 - g2,mg,0)         -> SCB(1,2)
cc   B0fin(g1 + g2,ms,ms)        -> SCB(4,1)
c    B0fin(g1 + g2,mg,ms)        -> SCB(5,1)  new
c    B0fin(g1 + g2,0,0)          -> SCB(7,1)
c    B0fin(p1,ms,0)              -> SCB(3,4)  renamed for m1^2=ms^2
cc   B0fin(pj,mt,ms)             -> SCB(8,1)
cc   B0fin(pj,ms,0)              -> SCB(3,5)
cc   B0fin(pj,mg,0)              -> SCB(1,4)
c    B0fin(g1,ms,mg)             -> SCB(6,1)
c    B0fin(g1,mt,mt)             -> SCB(6,2)  new 
c    B0fin(g1,ms,ms)             -> SCB(6,3)  new 
c    B0fin(g1,mg,mg)             -> SCB(6,4)  new 
cc   B0fin(ps,ms,0)              -> SCB(3,3)
c    B0fin(ps,mg,0)              -> SCB(1,3)  renamed for m1^2=ms^2
c
cc   C0fin(p1,pj,ms,0,ms)        -> SCC(6,1)
c    C0fin(p1,pj,0,ms,0)         -> SCC(5,1)
c    C0fin(p1,-g1,ms,0,0)        -> SCC(2,1)
cc   C0fin(p1,-g1,0,ms,mg)       -> SCC(1,1)
c    C0fin(p1,-g2,ms,0,0)        -> SCC(2,3)
cc   C0fin(p1,-g2,0,ms,mg)       -> SCC(1,3)
c    C0fin(pj,-g1,ms,0,0)        -> SCC(2,4)
cc   C0fin(pj,-g1,mg,0,0)        -> SCC(7,4)
c    C0fin(pj,-g1,0,ms,mg)       -> SCC(1,4)
cc   C0fin(pj,-g1,0,mg,ms)       -> SCC(8,4)
c    C0fin(pj,-g2,ms,0,0)        -> SCC(2,2)
cc   C0fin(pj,-g2,mg,0,0)        -> SCC(7,2)
cc   C0fin(pj,-g2,0,ms,mg)       -> SCC(1,2)
cc   C0fin(pj,-g2,0,mg,ms)       -> SCC(8,2)
cc   C0fin(g1,g2,ms,mg,ms)       -> SCC(3,1)
c    C0fin(g1,g2,0,0,0)          -> SCC(4,1)
c    C0fin(p1,-g2,0,ms,ms)       -> SCC(9,3)  new
c    C0fin(pj,-g2,0,ms,ms)       -> SCC(9,2)  new
c    C0fin(p1,pj,mg,0,ms)        -> SCC(10,1) new
c    C0fin(p1,-g2,mg,0,0)        -> SCC(2,5)  new
c    C0fin(p1,-g1,0,mg,ms)       -> SCC(8,1)  new
c    C0fin(p1,-g2,0,mg,mg)       -> SCC(11,3) new 
c    C0fin(g1,g2,ms,mg,mg)       -> SCC(12,1) new 
c    C0fin(g1,g2,mg,ms,ms)       -> SCC(13,1) new

c
cc   D0fin(p1,pj,-g1,ms,0,ms,mg) -> SCD(2,2)
c    D0fin(p1,pj,-g1,0,ms,0,0)   -> SCD(1,2)
cc   D0fin(pj,p1,-g1,ms,0,ms,mg) -> SCD(2,1)
c    D0fin(pj,p1,-g1,0,ms,0,0)   -> SCD(1,1)
c    D0fin(pj,-g2,p1,mg,0,0,ms)  -> SCD(3,1)
c    D0fin(pj,-g2,p1,0,mg,ms,0)  -> SCD(3,2)
c    D0fin(pj,-g2,p1,0,ms,ms,0)  -> SCD(4,1)  new
c    D0fin(p1,pj,-g2,mg,0,ms,ms) -> SCD(5,1)  new
c    D0fin(pj,p1,-g2,ms,0,mg,mg) -> SCD(6,2)  new 
c    D0fin(pj,-g2,p1,ms,0,0,mg)  -> SCD(7,1)  new 
c
cc   B0pfin(pj,mt,ms)            -> SCBP(2)
cc   B0pfin(pj,mg,0)             -> SCBP(3)
cc   B0pfin(pj,0,ms)             -> SCBP(4)
c    B0pfin(g1,ms,mg)            -> SCBP(1)
c    B0pfin(p1,ms,0)             -> SCBP(5)  new
c    B0pfin(p1,0,mg)             -> SCBP(6)  new

      gqFORT1 = 0.D0
      gqFORT2 = 0.D0

c               the prefactor for the scaling functions 
      NS_QGV = 0.D0

      end


