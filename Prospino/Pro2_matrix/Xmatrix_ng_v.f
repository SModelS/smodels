cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NG_QBB(MASSIN,CL,CR)                                             c
c     NG_QBV(MASSIN,CL,CR,CV)                                          c
c     NG_QBD(MASSIN,CL,CR)                                             c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = tg                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = mg                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(20) = del                                               c
c       MASSIN(21) = s4p                                               c
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
c       Q(K1) + QB(K2) -> NI(P1) + GL(P2) [+G(K3)]                     c
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
      real*8 function NG_QBB_NG(massin,Cl,Cr,mst,msu)

      implicit none 

      integer    n
      real*8     massin(1:30),pi,Nc,Cf
     &          ,s,m1,mg,t1,u1,tg,ug,tsl,tsr,usl,usr
     &          ,mst(-1:1),msu(-1:1)
      complex*16 Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_tx,QBB_tm,QBB_ux,QBB_um 

      pi= 4.D0 * atan(1.D0)

      Nc = 3.D0
      Cf = 4.D0/3.D0

      s  = massin(1)
      tg = massin(2)
      m1 = massin(6)
      mg = massin(7)

      t1 = tg + mg**2 - m1**2
      tsl = tg + mg**2 - mst(-1)**2 
      tsr = tg + mg**2 - mst(+1)**2 

      u1 = - s - tg
      ug = u1 + m1**2 - mg**2 
      usl = u1 + m1**2 - msu(-1)**2 
      usr = u1 + m1**2 - msu(+1)**2 

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      QBB_tm = dcmplx(0.D0,0.D0)
      QBB_tx = dcmplx(0.D0,0.D0)
      QBB_um = dcmplx(0.D0,0.D0)
      QBB_ux = dcmplx(0.D0,0.D0)

c               wim's form output 
      QBB_tx =   16.D0 * Nc*Cf * t1*tg 
     &       * ( Cl(1)*Clc(1)/tsl**2 + Cr(1)*Crc(1)/tsr**2 ) 

      QBB_tm = - 16.D0 * Nc*Cf * s*m1*mg 
     &       * ( Clc(1)*Clc(3)/tsl/usl + Crc(1)*Crc(3)/tsr/usr )

      QBB_ux =   16.D0 * Nc*Cf * u1*ug 
     &       * ( Cl(3)*Clc(3)/usl**2 + Cr(3)*Crc(3)/usr**2 )

      QBB_um = - 16.D0 * Nc*Cf * s*m1*mg 
     &       * ( Cl(1)*Cl(3)/tsl/usl + Cr(1)*Cr(3)/tsr/usr )
      
      NG_QBB_NG = real( QBB_tx + QBB_tm + QBB_ux + QBB_um )/2.D0

c               the phase space except for 1/s**2 
      NG_QBB_NG = NG_QBB_NG / ( 16.D0 * pi )

c               the averaging factors
      NG_QBB_NG = NG_QBB_NG /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NG_QBB_NG = NG_QBB_NG * (abs(m1)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_QBB(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),pi,Nc,Cf
     &          ,s,m1,mg,ms,t1,u1,tg,ug,ts,us
      complex*16 Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_tx,QBB_tm,QBB_ux,QBB_um 

      pi= 4.D0 * atan(1.D0)

      Nc = 3.D0
      Cf = 4.D0/3.D0

      s  = massin(1)
      tg = massin(2)
      m1 = massin(6)
      mg = massin(7)
      ms = massin(11)

      t1 = tg + mg**2 - m1**2
      ts = tg + mg**2 - ms**2 

      u1 = - s - tg
      ug = u1 + m1**2 - mg**2 
      us = u1 + m1**2 - ms**2 

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      QBB_tm = dcmplx(0.D0,0.D0)
      QBB_tx = dcmplx(0.D0,0.D0)
      QBB_um = dcmplx(0.D0,0.D0)
      QBB_ux = dcmplx(0.D0,0.D0)

c               wim's form output 
      QBB_tx =   16.D0 * Nc*Cf * t1*tg * ts**(-2) 
     &       * ( Cl(1)*Clc(1) + Cr(1)*Crc(1) ) 

      QBB_tm = - 16.D0 * Nc*Cf * s*m1*mg * ts**(-1)*us**(-1) 
     &       * ( Clc(1)*Clc(3) + Crc(1)*Crc(3) )

      QBB_ux =   16.D0 * Nc*Cf * u1*ug * us**(-2)
     &       * ( Cl(3)*Clc(3) + Cr(3)*Crc(3) )

      QBB_um = - 16.D0 * Nc*Cf * s*m1*mg * ts**(-1)*us**(-1)
     &       * ( Cl(1)*Cl(3) + Cr(1)*Cr(3) )
      
      NG_QBB = real( QBB_tx + QBB_tm + QBB_ux + QBB_um )/2.D0

c               the phase space except for 1/s**2 
      NG_QBB = NG_QBB / ( 16.D0 * pi )

c               the averaging factors
      NG_QBB = NG_QBB /4.D0 /Nc**2

c               the prefactor for the scaling functions 
c                   alpha is cut out of the typical couplings 
      NG_QBB = NG_QBB * (abs(m1)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      subroutine BORN_PARTS_NG(massin,Cl,Cr,b_tx,b_tm,b_ux,b_um)

      implicit none 

      integer    n
      real*8     massin(1:30),Nc,Cf
     &          ,s,m1,mg,ms,t1,u1,tg,ug,ts,us
     &          ,b_tx,b_tm,b_ux,b_um
      complex*16 Cl(4),Clc(4),Cr(4),Crc(4)
     &          ,QBB_tx,QBB_tm,QBB_ux,QBB_um

      Nc = 3.D0
      Cf = 4.D0/3.D0

      s  = massin(1)
      tg = massin(2)
      m1 = massin(6)
      mg = massin(7)
      ms = massin(11)

c               born kinematics built in
      t1 = tg + mg**2 - m1**2
      ts = tg + mg**2 - ms**2 

      u1 = - s - tg
      ug = u1 + m1**2 - mg**2 
      us = u1 + m1**2 - ms**2 

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

c               wim's form output 
      QBB_tx =   16.D0 * t1*tg * ts**(-2) 
     &       * ( Cl(1)*Clc(1) + Cr(1)*Crc(1) ) 

      QBB_tm = - 16.D0 * s*m1*mg * ts**(-1)*us**(-1) 
     &       * ( Clc(1)*Clc(3) + Crc(1)*Crc(3) )

      QBB_ux =   16.D0 * u1*ug * us**(-2)
     &       * ( Cl(3)*Clc(3) + Cr(3)*Crc(3) )

      QBB_um = - 16.D0 * s*m1*mg * ts**(-1)*us**(-1)
     &       * ( Cl(1)*Cl(3) + Cr(1)*Cr(3) )
      
      b_tx  = Nc*Cf*real( QBB_tx )
      b_tm  = Nc*Cf*real( QBB_tm )
      b_ux  = Nc*Cf*real( QBB_ux )
      b_um  = Nc*Cf*real( QBB_um )

c               check if the entries are really not complex 
      if (abs(aimag(QBB_tx)).gt.1.D-12) print *,'QBB_tx: ',QBB_tx
      if (abs(aimag(QBB_tm)).gt.1.D-12) print *,'QBB_tm: ',QBB_tm
      if (abs(aimag(QBB_ux)).gt.1.D-12) print *,'QBB_ux: ',QBB_ux
      if (abs(aimag(QBB_um)).gt.1.D-12) print *,'QBB_um: ',QBB_um

      end


c --------------------------------------------------------------------
c extract the log(Delta) terms from the virtual+soft scaling function 
c everything else copied from NG_QBV 
      real*8 function NG_QBD(massin,Cl,Cr)

      implicit none 

      real*8     massin(1:30),Pi,Nc,CF,CA,Sn,gs
     &          ,six,born,mu,QF,Del,log_del,log2_del
     &          ,s,mi,mx,mg,ti,ui,tg,ug,qqFORT,s4,s4p
     &          ,b_tx,b_tm,b_ux,b_um
      complex*16 Cl(4),Cr(4)

      Pi    = 4.D0*atan(1.D0)
      Nc    = 3.D0
      CF    = 4.D0/3.D0 
      CA    = 3.D0
      Sn    = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      six = 6.D0

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      tg  = massin(2)
      s4  = massin(3)
      mi  = massin(6)
      mg  = massin(7)
      Del = massin(20)
      s4p = massin(21)

c               the logaritms for linear s4 integration 
      log_del  =  log(s4p/mi**2)/(s4p-Del) - 1.D0/s4
      log2_del =  log(s4p/mi**2)**2/(s4p-Del) 
     &          - 2.D0*log(s4/mi**2)/s4

c               born kinematics built in as for BORN_PARTS_NG 
      ti = tg + mg**2 - mi**2

      ui = - s - tg 
      ug = ui + mi**2 - mg**2 

c               the absolute value for everything from phase space
      mx = abs(mi)

c               to compare to klasen: mi used in the renormalization of alpha_s
      mu = mx
ctp      mu = (mg+mx)/2.D0
      
c               the factorization scale 
      QF = massin(13) 

c               the born type structures 
      call BORN_PARTS_NG(massin,Cl,Cr,b_tx,b_tm,b_ux,b_um)

c               wim's conventions 
      born     = ( b_tx + b_tm + b_ux + b_um )/2.D0 
      
c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               'grep cdel' and replace log(Delta) terms 
      qqFORT = 0.D0
      qqFORT = qqFORT + log(mi**(-2)*s)*log_del*gs**2*Sn*CA
     + *born*six**(-1) * (  - 24 )
      qqFORT = qqFORT + log(mi**(-2)*s)*log_del*gs**2*Sn*CF
     + *born*six**(-1) * ( 48 )
      qqFORT = qqFORT + log( - mi**(-2)*ti)*log_del*gs**2*
     + Sn*CF*born*six**(-1) * (  - 48 )
      qqFORT = qqFORT + log( - mi**(-2)*ui)*log_del*gs**2*
     + Sn*CF*born*six**(-1) * (  - 48 )
      qqFORT = qqFORT + log_del*gs**2*Sn*CA*born*six**(-1)
     +  * (  - 24 )
      qqFORT = qqFORT + log2_del*gs**2*Sn*
     + CF*born*six**(-1) * ( 48 )
      qqFORT = qqFORT + log_del*log(mu**(-2)*QF**2)*gs**2*
     + Sn*CF*born*six**(-1) * (  - 48 )
      qqFORT = qqFORT + log_del*log( - mx**(-1)*mg**(-1)*tg
     + )*gs**2*Sn*CA*born*six**(-1) * ( 24 )
      qqFORT = qqFORT + log_del*log( - mx**(-1)*mg**(-1)*ug
     + )*gs**2*Sn*CA*born*six**(-1) * ( 24 )

c               the phase space except for 1/s**2 
      qqFORT = qqFORT / ( 16.D0 * pi )

c               the averaging factors
      qqFORT = qqFORT /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      NG_QBD = qqFORT * (abs(mi)+mg)**2/4.D0 

      end




c --------------------------------------------------------------------
      real*8 function NG_QBV(massin,Cl,Cr,Cv)

      implicit none 

      real*8     massin(1:30),Pi,Nc,CF,CA,Sn,nf,zeta2,gs
     &          ,six,born,rootlam,xs,QF,QR,mu
     &          ,tbornint,tintmimg,tintnonmimg
     &          ,ubornint,uintmimg,uintnonmimg
     &          ,s,mi,mx,mg,ms,mt,ti,ui,ts,us,tg,ug,qqFORT,Cupu,Clot
     &          ,SCA(1:3),SCB(1:8,1:5),SCBP(4),SCC(1:8,1:4),SCD(1:3,1:2)
     &          ,Li2,b_tx,b_tm,b_ux,b_um
      complex*16 CSPEN,Cl(4),Cr(4),Cv(4)

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
      tg  = massin(2)
      mi  = massin(6)
      mg  = massin(7)
      mt  = massin(9)
      ms  = massin(11)

c               born kinematics built in
      ti = tg + mg**2 - mi**2
      ts = tg + mg**2 - ms**2 

      ui = - s - tg
      ug = ui + mi**2 - mg**2 
      us = ui + mi**2 - ms**2 

c               the absolute value for everything from phase space
      mx = abs(mi)

c               to compare to klasen: mi used in the renormalization of alpha_s
      mu = mx
ckla      mu = (mg+mx)/2.D0
      
      rootlam = sqrt( s**2 +mi**4 +mg**4
     &                - 2*( s*mi**2 + s*mg**2 + mi**2*mg**2 ))

      xs = (sqrt( 1.D0 - 4.D0*abs(mi)*mg/(s-(abs(mi)-mg)**2) ) - 1.D0)
     &    /(sqrt( 1.D0 - 4.D0*abs(mi)*mg/(s-(abs(mi)-mg)**2) ) + 1.D0)

      Clot = real( Cv(1) ) 
      Cupu = real( Cv(3) ) 

c               the factorization/renormalization scale 
      QR = massin(12)
      QF = massin(13) 

c               the scalar functions 
      call SCALAR_ARRAY_NG(massin,SCA,SCB,SCBP,SCC,SCD)

c               the born type structures 
      call BORN_PARTS_NG(massin,Cl,Cr,b_tx,b_tm,b_ux,b_um)

c               wim's conventions 
      born     = ( b_tx + b_tm + b_ux + b_um )/2.D0 
      tbornint    = b_tx + b_tm 
      tintmimg    = b_tm
      tintnonmimg = b_tx
      ubornint    = b_ux + b_um
      uintmimg    = b_um
      uintnonmimg = b_ux
      
c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               insert the form output 
c
c   commented using cdel: scale dependence
c   commented using csof: finite soft contributions 
c
c   change in the phase space logs etc: mi -> mx
c
c   change the scalar integrals : 
c
c    A0fin(mt)                   -> SCA(1)
c    A0fin(ms)                   -> SCA(2)
c    A0fin(mg)                   -> SCA(3)
c
c    B0fin(pi - g1,ms,0)         -> SCB(3,1)
c    B0fin(pi - g1,mg,0)         -> SCB(1,1)
c    B0fin(pi - g2,ms,0)         -> SCB(3,2)
c    B0fin(pi - g2,mg,0)         -> SCB(1,2)
c    B0fin(g1 + g2,ms,ms)        -> SCB(4,1)
c    B0fin(g1 + g2,0,0)          -> SCB(7,1)
c    B0fin(pi,ms,0)              -> SCB(3,4)
c    B0fin(p2,mt,ms)             -> SCB(8,1)  new 
c    B0fin(p2,ms,0)              -> SCB(3,5)
c    B0fin(p2,mg,0)              -> SCB(1,4)  new 
c    B0fin(g1,ms,mg)             -> SCB(6,1)
c    B0fin(ps,ms,0)              -> SCB(3,3)
c    B0fin(ps,mg,0)              -> SCB(1,3)
c
c    C0fin(pi,p2,ms,0,ms)        -> SCC(6,1)
c    C0fin(pi,p2,0,ms,0)         -> SCC(5,1)
c    C0fin(pi,-g1,ms,0,0)        -> SCC(2,1)
c    C0fin(pi,-g1,0,ms,mg)       -> SCC(1,1)
c    C0fin(pi,-g2,ms,0,0)        -> SCC(2,3)
c    C0fin(pi,-g2,0,ms,mg)       -> SCC(1,3)
c    C0fin(p2,-g1,ms,0,0)        -> SCC(2,4)
c    C0fin(p2,-g1,mg,0,0)        -> SCC(7,4)  new
c    C0fin(p2,-g1,0,ms,mg)       -> SCC(1,4)
c    C0fin(p2,-g1,0,mg,ms)       -> SCC(8,4)  new
c    C0fin(p2,-g2,ms,0,0)        -> SCC(2,2)
c    C0fin(p2,-g2,mg,0,0)        -> SCC(7,2)  new 
c    C0fin(p2,-g2,0,ms,mg)       -> SCC(1,2)
c    C0fin(p2,-g2,0,mg,ms)       -> SCC(8,2)  new 
c    C0fin(g1,g2,ms,mg,ms)       -> SCC(3,1)
c    C0fin(g1,g2,0,0,0)          -> SCC(4,1)
c
c    D0fin(pi,p2,-g1,ms,0,ms,mg) -> SCD(2,2)
c    D0fin(pi,p2,-g1,0,ms,0,0)   -> SCD(1,2)
c    D0fin(p2,pi,-g1,ms,0,ms,mg) -> SCD(2,1)
c    D0fin(p2,pi,-g1,0,ms,0,0)   -> SCD(1,1)
c    D0fin(p2,-g2,pi,mg,0,0,ms)  -> SCD(3,1)  new
c    D0fin(p2,-g2,pi,0,mg,ms,0)  -> SCD(3,2)  new
c
c    B0pfin(p2,mt,ms)            -> SCBP(2)  new
c    B0pfin(p2,mg,0)             -> SCBP(3)  new
c    B0pfin(p2,0,ms)             -> SCBP(4)  new
c    B0pfin(g1,ms,mg)            -> SCBP(1)

      qqFORT = 0.D0
      qqFORT = gs**2*Sn*CA*tbornint*six**(-1) * (  - 6 )
      qqFORT = qqFORT + gs**2*Sn*CA*ubornint*six**(-1) * (  - 6 )
      qqFORT = qqFORT + gs**2*Sn*CA*born*six**(-1) * ( 14 )
      qqFORT = qqFORT + gs**2*Sn*CF*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + gs**2*Sn*CF*born*six**(-1) * (  - 36*zeta2 )
      qqFORT = qqFORT + SCA(1)*gs**2*Sn*born*six**(-1) * (  - 6*
     +    mg**(-2) )
      qqFORT = qqFORT + SCA(2)*gs**2*Sn*nf*born*six**(-1) * ( 6*
     +    mg**(-2) )
      qqFORT = qqFORT + SCA(3)*gs**2*Sn*CA*born*six**(-1) * (  - 6*
     +    mg**(-2) )
      qqFORT = qqFORT + SCB(3,1)*gs**2*Sn*CF*tbornint*
     + six**(-1) * ( 24*ms**2*ts**(-1) - 12*mi**2*ti**(-1) - 12*mg**2*
     +    tg**(-1) )
      qqFORT = qqFORT + SCB(3,1)*gs**2*Sn*CF*tintmimg*
     + six**(-1) * (  - 12*s**(-1)*ts )
      qqFORT = qqFORT + SCB(3,1)*gs**2*Sn*CF*tintnonmimg*
     + six**(-1) * ( 12*mg**2*ti**(-1)*tg**(-1)*ts + 12*ti**(-1)*ts )
      qqFORT = qqFORT + SCB(1,1)*gs**2*Sn*CF*Clot*tbornint*
     + six**(-1) * (  - 12*mi*mg*ti**(-1) )
      qqFORT = qqFORT + SCB(1,1)*gs**2*Sn*CF*tbornint*
     + six**(-1) * (  - 12 - 12*ms**2*ts**(-1) - 12*mg**2*tg**(-1) + 12
     +    *mg**2*ts**(-1) )
      qqFORT = qqFORT + SCB(1,1)*gs**2*Sn*CF*uintmimg*
     + six**(-1) * (  - 12*s**(-1)*us )
      qqFORT = qqFORT + SCB(1,1)*gs**2*Sn*CF*uintnonmimg*
     + six**(-1) * (  - 12*mi**2*ti**(-1)*ui**(-1)*us - 12*mg**2*
     +    tg**(-1)*ug**(-1)*us - 12*mg**2*ui**(-1)*ug**(-1)*us - 12*
     +    ui**(-1)*us )
      qqFORT = qqFORT + SCB(3,2)*gs**2*Sn*CF*ubornint*
     + six**(-1) * ( 24*ms**2*us**(-1) - 12*mi**2*ui**(-1) - 12*mg**2*
     +    ug**(-1) )
      qqFORT = qqFORT + SCB(3,2)*gs**2*Sn*CF*uintmimg*
     + six**(-1) * (  - 12*s**(-1)*us )
      qqFORT = qqFORT + SCB(3,2)*gs**2*Sn*CF*uintnonmimg*
     + six**(-1) * ( 12*mg**2*ui**(-1)*ug**(-1)*us + 12*ui**(-1)*us )
      qqFORT = qqFORT + SCB(1,2)*gs**2*Sn*CF*Cupu*ubornint*
     + six**(-1) * (  - 12*mi*mg*ui**(-1) )
      qqFORT = qqFORT + SCB(1,2)*gs**2*Sn*CF*ubornint*
     + six**(-1) * (  - 12 - 12*ms**2*us**(-1) - 12*mg**2*ug**(-1) + 12
     +    *mg**2*us**(-1) )
      qqFORT = qqFORT + SCB(1,2)*gs**2*Sn*CF*tintmimg*
     + six**(-1) * (  - 12*s**(-1)*ts )
      qqFORT = qqFORT + SCB(1,2)*gs**2*Sn*CF*tintnonmimg*
     + six**(-1) * (  - 12*mi**2*ti**(-1)*ui**(-1)*ts - 12*mg**2*
     +    ti**(-1)*tg**(-1)*ts - 12*mg**2*tg**(-1)*ug**(-1)*ts - 12*
     +    ti**(-1)*ts )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*tintmimg
     + *rootlam**(-2)*six**(-1) * ( 6*mi**2*ts - 6*mg**2*ts - 6*s*ts - 
     +    12*tg*ts )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*rootlam**(-2)*six**(-1) * (  - 6*mi**2*mg**2*
     +    ti**(-1)*ts - 6*mi**2*mg**2*tg**(-1)*ts - 6*mi**2*s*ti**(-1)*
     +    ts + 6*mi**4*ti**(-1)*ts - 24*mg**2*s*ti**(-1)*ts + 18*mg**2*
     +    s*tg**(-1)*ts - 24*mg**4*s*ti**(-1)*tg**(-1)*ts + 6*mg**4*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * (  - 6*mg**2*ti**(-1)*tg**(-1)*ts - 6*
     +    ti**(-1)*ts )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*uintmimg
     + *rootlam**(-2)*six**(-1) * ( 6*mi**2*us - 6*mg**2*us - 6*s*us - 
     +    12*ug*us )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*rootlam**(-2)*six**(-1) * (  - 6*mi**2*mg**2*
     +    ui**(-1)*us - 6*mi**2*mg**2*ug**(-1)*us - 6*mi**2*s*ui**(-1)*
     +    us + 6*mi**4*ui**(-1)*us - 24*mg**2*s*ui**(-1)*us + 18*mg**2*
     +    s*ug**(-1)*us - 24*mg**4*s*ui**(-1)*ug**(-1)*us + 6*mg**4*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCB(4,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * (  - 6*mg**2*ui**(-1)*ug**(-1)*us - 6*
     +    ui**(-1)*us )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*tintmimg*
     + rootlam**(-2)*six**(-1) * (  - 6*mi**2*ts + 6*mg**2*ts + 6*s*ts
     +     + 12*tg*ts )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*rootlam**(-2)*six**(-1) * ( 6*mi**2*mg**2*ti**(-1)*
     +    ts + 6*mi**2*mg**2*tg**(-1)*ts + 6*mi**2*s*ti**(-1)*ts - 6*
     +    mi**4*ti**(-1)*ts + 24*mg**2*s*ti**(-1)*ts - 18*mg**2*s*
     +    tg**(-1)*ts + 24*mg**4*s*ti**(-1)*tg**(-1)*ts - 6*mg**4*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 6*mg**2*ti**(-1)*tg**(-1)*ts + 6*
     +    ti**(-1)*ts )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*uintmimg*
     + rootlam**(-2)*six**(-1) * (  - 6*mi**2*us + 6*mg**2*us + 6*s*us
     +     + 12*ug*us )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*rootlam**(-2)*six**(-1) * ( 6*mi**2*mg**2*ui**(-1)*
     +    us + 6*mi**2*mg**2*ug**(-1)*us + 6*mi**2*s*ui**(-1)*us - 6*
     +    mi**4*ui**(-1)*us + 24*mg**2*s*ui**(-1)*us - 18*mg**2*s*
     +    ug**(-1)*us + 24*mg**4*s*ui**(-1)*ug**(-1)*us - 6*mg**4*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCB(7,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 6*mg**2*ui**(-1)*ug**(-1)*us + 6*
     +    ui**(-1)*us )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*Clot*tbornint*
     + six**(-1) * ( 12*mi*mg*ti**(-1) )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*Cupu*ubornint*
     + six**(-1) * ( 12*mi*mg*ui**(-1) )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*tbornint*six**(-1)
     +  * ( 12*mi**2*ti**(-1) )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*ubornint*six**(-1)
     +  * ( 12*mi**2*ui**(-1) )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*tintmimg*six**(-1)
     +  * ( 12*s**(-1)*ts )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*tintnonmimg*
     + six**(-1) * ( 12*mi**2*ti**(-1)*ui**(-1)*ts )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*uintmimg*six**(-1)
     +  * ( 12*s**(-1)*us )
      qqFORT = qqFORT + SCB(3,4)*gs**2*Sn*CF*uintnonmimg*
     + six**(-1) * ( 12*mi**2*ti**(-1)*ui**(-1)*us )
      qqFORT = qqFORT + SCB(8,1)*gs**2*Sn*born*six**(-1) * (  - 
     +    6 + 6*mt**2*mg**(-2) - 6*ms**2*mg**(-2) )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*tbornint*
     + six**(-1) * (  - 12*mg**2*tg**(-1) )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*ubornint*
     + six**(-1) * (  - 12*mg**2*ug**(-1) )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*tintmimg*
     + six**(-1) * (  - 6*s**(-1)*ts )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*tintnonmimg*
     + six**(-1) * (  - 6*mg**2*tg**(-1)*ug**(-1)*ts )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*uintmimg*
     + six**(-1) * (  - 6*s**(-1)*us )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*CA**(-1)*uintnonmimg*
     + six**(-1) * (  - 6*mg**2*tg**(-1)*ug**(-1)*us )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*nf*born*six**(-1) * ( 
     +     - 6 - 6*ms**2*mg**(-2) )
      qqFORT = qqFORT + SCB(3,5)*gs**2*Sn*born*six**(-1) * ( 6 + 
     +    6*ms**2*mg**(-2) )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*tbornint*six**(-1)
     +  * ( 12 + 12*mg**2*tg**(-1) )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*ubornint*six**(-1)
     +  * ( 12 + 12*mg**2*ug**(-1) )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*tintmimg*six**(-1)
     +  * ( 6*s**(-1)*ts )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*tintnonmimg*
     + six**(-1) * ( 6*mg**2*tg**(-1)*ug**(-1)*ts )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*uintmimg*six**(-1)
     +  * ( 6*s**(-1)*us )
      qqFORT = qqFORT + SCB(1,4)*gs**2*Sn*CA*uintnonmimg*
     + six**(-1) * ( 6*mg**2*tg**(-1)*ug**(-1)*us )
      qqFORT = qqFORT + SCB(6,1)*gs**2*Sn*CF*born*six**(-1) * ( 
     +     - 12 )
      qqFORT = qqFORT + SCB(3,3)*gs**2*Sn*CF*tbornint*six**(-1)
     +  * (  - 24*ms**2*ts**(-1) )
      qqFORT = qqFORT + SCB(3,3)*gs**2*Sn*CF*ubornint*six**(-1)
     +  * (  - 24*ms**2*us**(-1) )
      qqFORT = qqFORT + SCB(1,3)*gs**2*Sn*CF*tbornint*six**(-1)
     +  * ( 12*ms**2*ts**(-1) - 12*mg**2*ts**(-1) )
      qqFORT = qqFORT + SCB(1,3)*gs**2*Sn*CF*ubornint*six**(-1)
     +  * ( 12*ms**2*us**(-1) - 12*mg**2*us**(-1) )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*tintmimg
     + *rootlam**(-2)*six**(-1) * (  - 6*ms**2*mi**2*ts + 6*ms**2*mg**2
     +    *ts + 6*ms**2*s*ts + 12*ms**2*tg*ts - 6*mi**2*mg**2*ts + 6*
     +    mi**2*tg*ts - 6*mg**2*s*ts + 6*mg**2*tg*ts + 6*mg**4*ts - 6*s
     +    *tg*ts )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*tintmimg
     + *six**(-1) * (  - 6*ms**2*s**(-1)*ts - 3*mi**2*s**(-1)*ts + 9*
     +    mg**2*s**(-1)*ts + 6*s**(-1)*tg*ts + 3*ts )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*rootlam**(-2)*six**(-1) * ( 6*ms**2*mi**2*mg**2*
     +    ti**(-1)*ts + 6*ms**2*mi**2*mg**2*tg**(-1)*ts + 6*ms**2*mi**2
     +    *s*ti**(-1)*ts - 6*ms**2*mi**4*ti**(-1)*ts + 24*ms**2*mg**2*s
     +    *ti**(-1)*ts - 18*ms**2*mg**2*s*tg**(-1)*ts + 24*ms**2*mg**4*
     +    s*ti**(-1)*tg**(-1)*ts - 6*ms**2*mg**4*tg**(-1)*ts - 18*mi**2
     +    *mg**2*s*ti**(-1)*ts + 6*mi**2*mg**2*s*tg**(-1)*ts - 6*mi**2*
     +    mg**4*ti**(-1)*ts + 6*mi**2*mg**4*tg**(-1)*ts + 6*mi**4*mg**2
     +    *ti**(-1)*ts - 6*mi**4*mg**2*tg**(-1)*ts - 24*mg**4*s*
     +    ti**(-1)*ts + 24*mg**4*s*tg**(-1)*ts - 24*mg**6*s*ti**(-1)*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 12*ms**2*mg**2*ti**(-1)*tg**(-1)*ts - 
     +    6*ms**2*s*ti**(-1)*tg**(-1)*ts + 6*ms**2*ti**(-1)*ts - 6*
     +    ms**2*tg**(-1)*ts - 3*mg**2*s*ti**(-1)*tg**(-1)*ts - 9*mg**2*
     +    ti**(-1)*ts + 9*mg**2*tg**(-1)*ts - 12*mg**4*ti**(-1)*
     +    tg**(-1)*ts - 3*s*ti**(-1)*ts + 3*s*tg**(-1)*ts + 3*s**2*
     +    ti**(-1)*tg**(-1)*ts )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*uintmimg
     + *rootlam**(-2)*six**(-1) * (  - 6*ms**2*mi**2*us + 6*ms**2*mg**2
     +    *us + 6*ms**2*s*us + 12*ms**2*ug*us - 6*mi**2*mg**2*us + 6*
     +    mi**2*ug*us - 6*mg**2*s*us + 6*mg**2*ug*us + 6*mg**4*us - 6*s
     +    *ug*us )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*uintmimg
     + *six**(-1) * (  - 6*ms**2*s**(-1)*us - 3*mi**2*s**(-1)*us + 9*
     +    mg**2*s**(-1)*us + 6*s**(-1)*ug*us + 3*us )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*rootlam**(-2)*six**(-1) * ( 6*ms**2*mi**2*mg**2*
     +    ui**(-1)*us + 6*ms**2*mi**2*mg**2*ug**(-1)*us + 6*ms**2*mi**2
     +    *s*ui**(-1)*us - 6*ms**2*mi**4*ui**(-1)*us + 24*ms**2*mg**2*s
     +    *ui**(-1)*us - 18*ms**2*mg**2*s*ug**(-1)*us + 24*ms**2*mg**4*
     +    s*ui**(-1)*ug**(-1)*us - 6*ms**2*mg**4*ug**(-1)*us - 18*mi**2
     +    *mg**2*s*ui**(-1)*us + 6*mi**2*mg**2*s*ug**(-1)*us - 6*mi**2*
     +    mg**4*ui**(-1)*us + 6*mi**2*mg**4*ug**(-1)*us + 6*mi**4*mg**2
     +    *ui**(-1)*us - 6*mi**4*mg**2*ug**(-1)*us - 24*mg**4*s*
     +    ui**(-1)*us + 24*mg**4*s*ug**(-1)*us - 24*mg**6*s*ui**(-1)*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCC(6,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 12*ms**2*mg**2*ui**(-1)*ug**(-1)*us - 
     +    6*ms**2*s*ui**(-1)*ug**(-1)*us + 6*ms**2*ui**(-1)*us - 6*
     +    ms**2*ug**(-1)*us - 3*mg**2*s*ui**(-1)*ug**(-1)*us - 9*mg**2*
     +    ui**(-1)*us + 9*mg**2*ug**(-1)*us - 12*mg**4*ui**(-1)*
     +    ug**(-1)*us - 3*s*ui**(-1)*us + 3*s*ug**(-1)*us + 3*s**2*
     +    ui**(-1)*ug**(-1)*us )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*tintmimg*
     + rootlam**(-2)*six**(-1) * (  - 6*ms**2*mi**2*ts + 6*ms**2*mg**2*
     +    ts + 6*ms**2*s*ts + 12*ms**2*tg*ts + 6*mi**2*mg**2*ts - 6*
     +    mi**2*tg*ts + 6*mg**2*s*ts - 6*mg**2*tg*ts - 6*mg**4*ts + 6*s
     +    *tg*ts )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*tintmimg*
     + six**(-1) * ( 12*ts )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*rootlam**(-2)*six**(-1) * ( 6*ms**2*mi**2*mg**2*
     +    ti**(-1)*ts + 6*ms**2*mi**2*mg**2*tg**(-1)*ts + 6*ms**2*mi**2
     +    *s*ti**(-1)*ts - 6*ms**2*mi**4*ti**(-1)*ts + 24*ms**2*mg**2*s
     +    *ti**(-1)*ts - 18*ms**2*mg**2*s*tg**(-1)*ts + 24*ms**2*mg**4*
     +    s*ti**(-1)*tg**(-1)*ts - 6*ms**2*mg**4*tg**(-1)*ts + 18*mi**2
     +    *mg**2*s*ti**(-1)*ts - 6*mi**2*mg**2*s*tg**(-1)*ts + 6*mi**2*
     +    mg**4*ti**(-1)*ts - 6*mi**2*mg**4*tg**(-1)*ts - 6*mi**4*mg**2
     +    *ti**(-1)*ts + 6*mi**4*mg**2*tg**(-1)*ts + 24*mg**4*s*
     +    ti**(-1)*ts - 24*mg**4*s*tg**(-1)*ts + 24*mg**6*s*ti**(-1)*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 6*ms**2*mg**2*ti**(-1)*tg**(-1)*ts - 3
     +    *ms**2*s*ti**(-1)*tg**(-1)*ts + 3*ms**2*ti**(-1)*ts - 3*ms**2
     +    *tg**(-1)*ts + 3*mi**2*tg**(-1)*ts + 3*mg**2*s*ti**(-1)*
     +    tg**(-1)*ts + 9*mg**2*ti**(-1)*ts - 6*mg**2*tg**(-1)*ts + 6*
     +    mg**4*ti**(-1)*tg**(-1)*ts - 3*s*tg**(-1)*ts )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*uintmimg*
     + rootlam**(-2)*six**(-1) * (  - 6*ms**2*mi**2*us + 6*ms**2*mg**2*
     +    us + 6*ms**2*s*us + 12*ms**2*ug*us + 6*mi**2*mg**2*us - 6*
     +    mi**2*ug*us + 6*mg**2*s*us - 6*mg**2*ug*us - 6*mg**4*us + 6*s
     +    *ug*us )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*uintmimg*
     + six**(-1) * ( 12*us )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*rootlam**(-2)*six**(-1) * ( 6*ms**2*mi**2*mg**2*
     +    ui**(-1)*us + 6*ms**2*mi**2*mg**2*ug**(-1)*us + 6*ms**2*mi**2
     +    *s*ui**(-1)*us - 6*ms**2*mi**4*ui**(-1)*us + 24*ms**2*mg**2*s
     +    *ui**(-1)*us - 18*ms**2*mg**2*s*ug**(-1)*us + 24*ms**2*mg**4*
     +    s*ui**(-1)*ug**(-1)*us - 6*ms**2*mg**4*ug**(-1)*us + 18*mi**2
     +    *mg**2*s*ui**(-1)*us - 6*mi**2*mg**2*s*ug**(-1)*us + 6*mi**2*
     +    mg**4*ui**(-1)*us - 6*mi**2*mg**4*ug**(-1)*us - 6*mi**4*mg**2
     +    *ui**(-1)*us + 6*mi**4*mg**2*ug**(-1)*us + 24*mg**4*s*
     +    ui**(-1)*us - 24*mg**4*s*ug**(-1)*us + 24*mg**6*s*ui**(-1)*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCC(5,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 6*ms**2*mg**2*ui**(-1)*ug**(-1)*us - 3
     +    *ms**2*s*ui**(-1)*ug**(-1)*us + 3*ms**2*ui**(-1)*us - 3*ms**2
     +    *ug**(-1)*us + 3*mi**2*ug**(-1)*us + 3*mg**2*s*ui**(-1)*
     +    ug**(-1)*us + 9*mg**2*ui**(-1)*us - 6*mg**2*ug**(-1)*us + 6*
     +    mg**4*ui**(-1)*ug**(-1)*us - 3*s*ug**(-1)*us )
      qqFORT = qqFORT + SCC(2,1)*gs**2*Sn*CA*tintmimg*
     + six**(-1) * (  - 3*mi**2*s**(-1)*ts + 3*mg**2*s**(-1)*ts + 3*
     +    s**(-1)*tg*ts )
      qqFORT = qqFORT + SCC(2,1)*gs**2*Sn*CF*tbornint*
     + six**(-1) * ( 12*ms**2 - 12*mi**2 )
      qqFORT = qqFORT + SCC(2,1)*gs**2*Sn*CF*tintmimg*
     + six**(-1) * ( 12*ts )
      qqFORT = qqFORT + SCC(2,1)*gs**2*Sn*CF*tintnonmimg*
     + six**(-1) * (  - 6*ms**2*tg**(-1)*ts + 6*mi**2*tg**(-1)*ts + 6*
     +    ts )
      qqFORT = qqFORT + SCC(1,1)*gs**2*Sn*CA*uintnonmimg*
     + six**(-1) * ( 3*ms**2*s*ui**(-1)*ug**(-1)*us + 3*ms**2*ui**(-1)*
     +    us - 3*s*ui**(-1)*us - 3*s*ug**(-1)*us - 3*s**2*ui**(-1)*
     +    ug**(-1)*us - 3*us )
      qqFORT = qqFORT + SCC(1,1)*gs**2*Sn*CF*Clot*tbornint
     + *six**(-1) * (  - 12*ms**2*mi*mg*ti**(-1) - 12*mi*mg + 12*mi*
     +    mg**3*ti**(-1) )
      qqFORT = qqFORT + SCC(1,1)*gs**2*Sn*CF*uintmimg*
     + six**(-1) * (  - 12*ms**2*s**(-1)*us + 12*mg**2*s**(-1)*us + 6*
     +    s**(-1)*ug*us + 6*us )
      qqFORT = qqFORT + SCC(1,1)*gs**2*Sn*CF*uintnonmimg*
     + six**(-1) * (  - 12*ms**2*mi**2*ti**(-1)*ui**(-1)*us - 12*ms**2*
     +    s*ui**(-1)*ug**(-1)*us - 12*ms**2*ui**(-1)*us + 12*mi**2*
     +    mg**2*ti**(-1)*ui**(-1)*us + 6*mg**2*s*ui**(-1)*ug**(-1)*us
     +     + 6*mg**2*ui**(-1)*us + 6*s*ui**(-1)*us + 6*s**2*ui**(-1)*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCC(2,3)*gs**2*Sn*CA*uintmimg*
     + six**(-1) * (  - 3*mi**2*s**(-1)*us + 3*mg**2*s**(-1)*us + 3*
     +    s**(-1)*ug*us )
      qqFORT = qqFORT + SCC(2,3)*gs**2*Sn*CF*ubornint*
     + six**(-1) * ( 12*ms**2 - 12*mi**2 )
      qqFORT = qqFORT + SCC(2,3)*gs**2*Sn*CF*uintmimg*
     + six**(-1) * ( 12*us )
      qqFORT = qqFORT + SCC(2,3)*gs**2*Sn*CF*uintnonmimg*
     + six**(-1) * (  - 6*ms**2*ug**(-1)*us + 6*mi**2*ug**(-1)*us + 6*
     +    us )
      qqFORT = qqFORT + SCC(1,3)*gs**2*Sn*CA*tintnonmimg*
     + six**(-1) * ( 3*ms**2*s*ti**(-1)*tg**(-1)*ts + 3*ms**2*ti**(-1)*
     +    ts - 3*s*ti**(-1)*ts - 3*s*tg**(-1)*ts - 3*s**2*ti**(-1)*
     +    tg**(-1)*ts - 3*ts )
      qqFORT = qqFORT + SCC(1,3)*gs**2*Sn*CF*Cupu*ubornint
     + *six**(-1) * (  - 12*ms**2*mi*mg*ui**(-1) - 12*mi*mg + 12*mi*
     +    mg**3*ui**(-1) )
      qqFORT = qqFORT + SCC(1,3)*gs**2*Sn*CF*tintmimg*
     + six**(-1) * (  - 12*ms**2*s**(-1)*ts + 12*mg**2*s**(-1)*ts + 6*
     +    s**(-1)*tg*ts + 6*ts )
      qqFORT = qqFORT + SCC(1,3)*gs**2*Sn*CF*tintnonmimg*
     + six**(-1) * (  - 12*ms**2*mi**2*ti**(-1)*ui**(-1)*ts - 12*ms**2*
     +    s*ti**(-1)*tg**(-1)*ts - 12*ms**2*ti**(-1)*ts + 12*mi**2*
     +    mg**2*ti**(-1)*ui**(-1)*ts + 6*mg**2*s*ti**(-1)*tg**(-1)*ts
     +     + 6*mg**2*ti**(-1)*ts + 6*s*ti**(-1)*ts + 6*s**2*ti**(-1)*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCC(2,4)*gs**2*Sn*CA**(-1)*ubornint
     + *six**(-1) * (  - 6*ms**2 + 6*mg**2 )
      qqFORT = qqFORT + SCC(2,4)*gs**2*Sn*CA**(-1)*uintmimg
     + *six**(-1) * (  - 6*us )
      qqFORT = qqFORT + SCC(2,4)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 3*ms**2*ui**(-1)*us - 3*mg**2*ui**(-1)
     +    *us - 3*us )
      qqFORT = qqFORT + SCC(7,4)*gs**2*Sn*CA*ubornint*
     + six**(-1) * (  - 6*ug )
      qqFORT = qqFORT + SCC(7,4)*gs**2*Sn*CA*tintmimg*
     + six**(-1) * ( 3*mi**2*s**(-1)*ts - 3*mg**2*s**(-1)*ts - 3*
     +    s**(-1)*tg*ts - 3*ts )
      qqFORT = qqFORT + SCC(7,4)*gs**2*Sn*CA*tintnonmimg*
     + six**(-1) * ( 3*ms**2*s*ti**(-1)*tg**(-1)*ts + 3*ms**2*tg**(-1)*
     +    ts - 3*mi**2*tg**(-1)*ts - 3*mg**2*s*ti**(-1)*tg**(-1)*ts + 3
     +    *s*tg**(-1)*ts + 3*ts )
      qqFORT = qqFORT + SCC(1,4)*gs**2*Sn*CA**(-1)*
     + ubornint*six**(-1) * ( 6*ms**2*mg**2*ug**(-1) + 6*mg**2 - 6*
     +    mg**4*ug**(-1) )
      qqFORT = qqFORT + SCC(1,4)*gs**2*Sn*CA**(-1)*
     + tintmimg*six**(-1) * ( 6*ms**2*s**(-1)*ts + 3*mi**2*s**(-1)*ts
     +     - 9*mg**2*s**(-1)*ts - 3*s**(-1)*tg*ts - 3*ts )
      qqFORT = qqFORT + SCC(1,4)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 6*ms**2*mg**2*tg**(-1)*ug**(-1)*ts + 6
     +    *ms**2*s*ti**(-1)*tg**(-1)*ts + 6*ms**2*tg**(-1)*ts - 3*mg**2
     +    *s*ti**(-1)*tg**(-1)*ts - 3*mg**2*tg**(-1)*ts - 6*mg**4*
     +    tg**(-1)*ug**(-1)*ts - 3*s*tg**(-1)*ts - 3*s**2*ti**(-1)*
     +    tg**(-1)*ts )
      qqFORT = qqFORT + SCC(8,4)*gs**2*Sn*CA*ubornint*
     + six**(-1) * ( 6*ms**2*mg**2*ug**(-1) + 6*ms**2 - 6*mg**4*
     +    ug**(-1) )
      qqFORT = qqFORT + SCC(8,4)*gs**2*Sn*CA*uintmimg*
     + six**(-1) * ( 6*ms**2*s**(-1)*us - 6*mg**2*s**(-1)*us - 3*
     +    s**(-1)*ug*us )
      qqFORT = qqFORT + SCC(8,4)*gs**2*Sn*CA*uintnonmimg*
     + six**(-1) * (  - 6*ms**2*mg**2*ui**(-1)*ug**(-1)*us - 3*ms**2*
     +    ui**(-1)*us + 3*mg**2*ui**(-1)*us + 6*mg**4*ui**(-1)*ug**(-1)
     +    *us + 3*us )
      qqFORT = qqFORT + SCC(2,2)*gs**2*Sn*CA**(-1)*tbornint
     + *six**(-1) * (  - 6*ms**2 + 6*mg**2 )
      qqFORT = qqFORT + SCC(2,2)*gs**2*Sn*CA**(-1)*tintmimg
     + *six**(-1) * (  - 6*ts )
      qqFORT = qqFORT + SCC(2,2)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 3*ms**2*ti**(-1)*ts - 3*mg**2*ti**(-1)
     +    *ts - 3*ts )
      qqFORT = qqFORT + SCC(7,2)*gs**2*Sn*CA*tbornint*
     + six**(-1) * (  - 6*tg )
      qqFORT = qqFORT + SCC(7,2)*gs**2*Sn*CA*uintmimg*
     + six**(-1) * ( 3*mi**2*s**(-1)*us - 3*mg**2*s**(-1)*us - 3*
     +    s**(-1)*ug*us - 3*us )
      qqFORT = qqFORT + SCC(7,2)*gs**2*Sn*CA*uintnonmimg*
     + six**(-1) * ( 3*ms**2*s*ui**(-1)*ug**(-1)*us + 3*ms**2*ug**(-1)*
     +    us - 3*mi**2*ug**(-1)*us - 3*mg**2*s*ui**(-1)*ug**(-1)*us + 3
     +    *s*ug**(-1)*us + 3*us )
      qqFORT = qqFORT + SCC(1,2)*gs**2*Sn*CA**(-1)*
     + tbornint*six**(-1) * ( 6*ms**2*mg**2*tg**(-1) + 6*mg**2 - 6*
     +    mg**4*tg**(-1) )
      qqFORT = qqFORT + SCC(1,2)*gs**2*Sn*CA**(-1)*
     + uintmimg*six**(-1) * ( 6*ms**2*s**(-1)*us + 3*mi**2*s**(-1)*us
     +     - 9*mg**2*s**(-1)*us - 3*s**(-1)*ug*us - 3*us )
      qqFORT = qqFORT + SCC(1,2)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 6*ms**2*mg**2*tg**(-1)*ug**(-1)*us + 6
     +    *ms**2*s*ui**(-1)*ug**(-1)*us + 6*ms**2*ug**(-1)*us - 3*mg**2
     +    *s*ui**(-1)*ug**(-1)*us - 3*mg**2*ug**(-1)*us - 6*mg**4*
     +    tg**(-1)*ug**(-1)*us - 3*s*ug**(-1)*us - 3*s**2*ui**(-1)*
     +    ug**(-1)*us )
      qqFORT = qqFORT + SCC(8,2)*gs**2*Sn*CA*tbornint*
     + six**(-1) * ( 6*ms**2*mg**2*tg**(-1) + 6*ms**2 - 6*mg**4*
     +    tg**(-1) )
      qqFORT = qqFORT + SCC(8,2)*gs**2*Sn*CA*tintmimg*
     + six**(-1) * ( 6*ms**2*s**(-1)*ts - 6*mg**2*s**(-1)*ts - 3*
     +    s**(-1)*tg*ts )
      qqFORT = qqFORT + SCC(8,2)*gs**2*Sn*CA*tintnonmimg*
     + six**(-1) * (  - 6*ms**2*mg**2*ti**(-1)*tg**(-1)*ts - 3*ms**2*
     +    ti**(-1)*ts + 3*mg**2*ti**(-1)*ts + 6*mg**4*ti**(-1)*tg**(-1)
     +    *ts + 3*ts )
      qqFORT = qqFORT + SCC(3,1)*gs**2*Sn*CA**(-1)*
     + tintmimg*six**(-1) * ( 3*ts )
      qqFORT = qqFORT + SCC(3,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * (  - 6*ms**2*s*ti**(-1)*tg**(-1)*ts + 3*
     +    mg**2*s*ti**(-1)*tg**(-1)*ts + 3*s**2*ti**(-1)*tg**(-1)*ts )
      qqFORT = qqFORT + SCC(3,1)*gs**2*Sn*CA**(-1)*
     + uintmimg*six**(-1) * ( 3*us )
      qqFORT = qqFORT + SCC(3,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * (  - 6*ms**2*s*ui**(-1)*ug**(-1)*us + 3*
     +    mg**2*s*ui**(-1)*ug**(-1)*us + 3*s**2*ui**(-1)*ug**(-1)*us )
      qqFORT = qqFORT + SCC(4,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 3*ms**2*s*ti**(-1)*tg**(-1)*ts - 3*
     +    mg**2*s*ti**(-1)*tg**(-1)*ts + 3*s*tg**(-1)*ts )
      qqFORT = qqFORT + SCC(4,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 3*ms**2*s*ui**(-1)*ug**(-1)*us - 3*
     +    mg**2*s*ui**(-1)*ug**(-1)*us + 3*s*ug**(-1)*us )
      qqFORT = qqFORT + SCD(2,2)*gs**2*Sn*CA**(-1)*
     + tintmimg*six**(-1) * ( 3*ms**2*mi**2*s**(-1)*ts - 15*ms**2*mg**2
     +    *s**(-1)*ts - 6*ms**2*s**(-1)*tg*ts - 6*ms**2*ts + 6*ms**4*
     +    s**(-1)*ts - 3*mi**2*mg**2*s**(-1)*ts - 3*mi**2*ts + 6*mg**2*
     +    s**(-1)*tg*ts + 9*mg**2*ts + 9*mg**4*s**(-1)*ts + 3*s*ts + 3*
     +    tg*ts )
      qqFORT = qqFORT + SCD(2,2)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * (  - 6*ms**2*mg**2*s*ti**(-1)*tg**(-1)*
     +    ts + 3*ms**2*mg**2*ti**(-1)*ts - 9*ms**2*mg**2*tg**(-1)*ts + 
     +    12*ms**2*mg**4*ti**(-1)*tg**(-1)*ts - 3*ms**2*s*ti**(-1)*ts
     +     - 9*ms**2*s*tg**(-1)*ts - 12*ms**2*s**2*ti**(-1)*tg**(-1)*ts
     +     - 6*ms**4*mg**2*ti**(-1)*tg**(-1)*ts + 12*ms**4*s*ti**(-1)*
     +    tg**(-1)*ts + 6*ms**4*tg**(-1)*ts - 3*mg**2*s*ti**(-1)*ts + 6
     +    *mg**2*s*tg**(-1)*ts + 3*mg**2*s**2*ti**(-1)*tg**(-1)*ts - 6*
     +    mg**4*s*ti**(-1)*tg**(-1)*ts - 3*mg**4*ti**(-1)*ts + 3*mg**4*
     +    tg**(-1)*ts - 6*mg**6*ti**(-1)*tg**(-1)*ts + 3*s**2*tg**(-1)*
     +    ts + 3*s**3*ti**(-1)*tg**(-1)*ts )
      qqFORT = qqFORT + SCD(1,2)*gs**2*Sn*CA**(-1)*
     + uintmimg*six**(-1) * ( 6*s*us )
      qqFORT = qqFORT + SCD(1,2)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * ( 3*s*ui**(-1)*ug**(-1)*us**3 - 3*s*
     +    ui**(-1)*us**2 - 3*s*ug**(-1)*us**2 + 6*s*us )
      qqFORT = qqFORT + SCD(2,1)*gs**2*Sn*CA**(-1)*
     + uintmimg*six**(-1) * ( 3*ms**2*mi**2*s**(-1)*us - 15*ms**2*mg**2
     +    *s**(-1)*us - 6*ms**2*s**(-1)*ug*us - 6*ms**2*us + 6*ms**4*
     +    s**(-1)*us - 3*mi**2*mg**2*s**(-1)*us - 3*mi**2*us + 6*mg**2*
     +    s**(-1)*ug*us + 9*mg**2*us + 9*mg**4*s**(-1)*us + 3*s*us + 3*
     +    ug*us )
      qqFORT = qqFORT + SCD(2,1)*gs**2*Sn*CA**(-1)*
     + uintnonmimg*six**(-1) * (  - 6*ms**2*mg**2*s*ui**(-1)*ug**(-1)*
     +    us + 3*ms**2*mg**2*ui**(-1)*us - 9*ms**2*mg**2*ug**(-1)*us + 
     +    12*ms**2*mg**4*ui**(-1)*ug**(-1)*us - 3*ms**2*s*ui**(-1)*us
     +     - 9*ms**2*s*ug**(-1)*us - 12*ms**2*s**2*ui**(-1)*ug**(-1)*us
     +     - 6*ms**4*mg**2*ui**(-1)*ug**(-1)*us + 12*ms**4*s*ui**(-1)*
     +    ug**(-1)*us + 6*ms**4*ug**(-1)*us - 3*mg**2*s*ui**(-1)*us + 6
     +    *mg**2*s*ug**(-1)*us + 3*mg**2*s**2*ui**(-1)*ug**(-1)*us - 6*
     +    mg**4*s*ui**(-1)*ug**(-1)*us - 3*mg**4*ui**(-1)*us + 3*mg**4*
     +    ug**(-1)*us - 6*mg**6*ui**(-1)*ug**(-1)*us + 3*s**2*ug**(-1)*
     +    us + 3*s**3*ui**(-1)*ug**(-1)*us )
      qqFORT = qqFORT + SCD(1,1)*gs**2*Sn*CA**(-1)*
     + tintmimg*six**(-1) * ( 6*s*ts )
      qqFORT = qqFORT + SCD(1,1)*gs**2*Sn*CA**(-1)*
     + tintnonmimg*six**(-1) * ( 3*s*ti**(-1)*tg**(-1)*ts**3 - 3*s*
     +    ti**(-1)*ts**2 - 3*s*tg**(-1)*ts**2 + 6*s*ts )
      qqFORT = qqFORT + SCD(3,1)*gs**2*Sn*CA*uintmimg
     + *six**(-1) * (  - 3*mi**2*s**(-1)*us**2 - 6*mi**2*us + 3*mg**2*
     +    s**(-1)*us**2 + 6*mg**2*us + 3*s**(-1)*ug*us**2 + 6*s*us + 6*
     +    ug*us + 3*us**2 )
      qqFORT = qqFORT + SCD(3,1)*gs**2*Sn*CA*
     + uintnonmimg*six**(-1) * ( 3*mi**2*ug**(-1)*us**2 - 6*mi**2*us - 
     +    3*mg**2*ug**(-1)*us**2 + 6*mg**2*us + 3*s*ui**(-1)*ug**(-1)*
     +    us**3 - 3*s*ui**(-1)*us**2 - 3*s*ug**(-1)*us**2 + 6*s*us + 3*
     +    ug**(-1)*us**3 + 6*ug*us - 6*us**2 )
      qqFORT = qqFORT + SCD(3,2)*gs**2*Sn*CA*tintmimg
     + *six**(-1) * (  - 3*mi**2*s**(-1)*ts**2 - 6*mi**2*ts + 3*mg**2*
     +    s**(-1)*ts**2 + 6*mg**2*ts + 3*s**(-1)*tg*ts**2 + 6*s*ts + 6*
     +    tg*ts + 3*ts**2 )
      qqFORT = qqFORT + SCD(3,2)*gs**2*Sn*CA*
     + tintnonmimg*six**(-1) * ( 3*mi**2*tg**(-1)*ts**2 - 6*mi**2*ts - 
     +    3*mg**2*tg**(-1)*ts**2 + 6*mg**2*ts + 3*s*ti**(-1)*tg**(-1)*
     +    ts**3 - 3*s*ti**(-1)*ts**2 - 3*s*tg**(-1)*ts**2 + 6*s*ts + 3*
     +    tg**(-1)*ts**3 + 6*tg*ts - 6*ts**2 )
      qqFORT = qqFORT + SCBP(2)*gs**2*Sn*born*six**(-1) * ( 
     +     - 12*mt**2 + 12*ms**2 - 12*mg**2 )
      qqFORT = qqFORT + SCBP(3)*gs**2*Sn*CA*born*six**(-1) * ( 
     +    24*mg**2 )
      qqFORT = qqFORT + SCBP(4)*gs**2*Sn*nf*born*six**(-1) * ( 
     +    12*ms**2 - 12*mg**2 )
      qqFORT = qqFORT + SCBP(4)*gs**2*Sn*born*six**(-1) * (  - 
     +    12*ms**2 + 12*mg**2 )
      qqFORT = qqFORT + SCBP(1)*gs**2*Sn*CF*born*six**(-1)
     +  * ( 12*ms**2 - 12*mg**2 )
csof      qqFORT = qqFORT + Li2(1 - mi**(-2)*s**(-1)*ti*ui)*gs**2*Sn*CA*
csof     + born*six**(-1) * ( 12 )
csof      qqFORT = qqFORT + Li2(1 - mi**(-2)*s**(-1)*ti*ui)*gs**2*Sn*CF*
csof     + born*six**(-1) * (  - 24 )
csof      qqFORT = qqFORT + Li2(1 + mx**(-1)*mg*ti*ug**(-1)*xs**(-1))*gs**2
csof     + *Sn*CA*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + Li2(1 + mx**(-1)*mg*tg**(-1)*ui*xs**(-1))*gs**2
csof     + *Sn*CA*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + Li2(1 + mx*mg**(-1)*ti**(-1)*ug*xs**(-1))*gs**2
csof     + *Sn*CA*born*six**(-1) * ( 12 )
csof      qqFORT = qqFORT + Li2(1 + mx*mg**(-1)*tg*ui**(-1)*xs**(-1))*gs**2
csof     + *Sn*CA*born*six**(-1) * ( 12 )
      qqFORT = qqFORT + log(mt**(-2)*QR**2)*gs**2*Sn*born*six**(-1)
     +  * ( 4 )
      qqFORT = qqFORT + log(ms**(-2)*QR**2)*gs**2*Sn*nf*born*six**(-1)
     +  * ( 2 )
cdel      qqFORT = qqFORT + log(mi**(-2)*s)*log(mi**(-2)*Del)*gs**2*Sn*CA
cdel     + *born*six**(-1) * (  - 24 )
cdel      qqFORT = qqFORT + log(mi**(-2)*s)*log(mi**(-2)*Del)*gs**2*Sn*CF
cdel     + *born*six**(-1) * ( 48 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ti)*log( - mi**(-2)*ti)*gs**2
csof     + *Sn*CA*born*six**(-1) * ( 6 )
cdel      qqFORT = qqFORT + log( - mi**(-2)*ti)*log(mi**(-2)*Del)*gs**2*
cdel     + Sn*CF*born*six**(-1) * (  - 48 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ti)*log(mu**(-2)*QF**2)*gs**2
csof     + *Sn*CF*born*six**(-1) * ( 24 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ti)*log( - mx**(-1)*mg**(-1)*
csof     + ug)*gs**2*Sn*CA*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ui)*log( - mi**(-2)*ui)*gs**2
csof     + *Sn*CA*born*six**(-1) * ( 6 )
cdel      qqFORT = qqFORT + log( - mi**(-2)*ui)*log(mi**(-2)*Del)*gs**2*
cdel     + Sn*CF*born*six**(-1) * (  - 48 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ui)*log(mu**(-2)*QF**2)*gs**2
csof     + *Sn*CF*born*six**(-1) * ( 24 )
csof      qqFORT = qqFORT + log( - mi**(-2)*ui)*log( - mx**(-1)*mg**(-1)*
csof     + tg)*gs**2*Sn*CA*born*six**(-1) * (  - 12 )
cdel      qqFORT = qqFORT + log(mi**(-2)*Del)*gs**2*Sn*CA*born*six**(-1)
cdel     +  * (  - 24 )
cdel      qqFORT = qqFORT + log(mi**(-2)*Del)*log(mi**(-2)*Del)*gs**2*Sn*
cdel     + CF*born*six**(-1) * ( 48 )
cdel      qqFORT = qqFORT + log(mi**(-2)*Del)*log(mu**(-2)*QF**2)*gs**2*
cdel     + Sn*CF*born*six**(-1) * (  - 48 )
cdel      qqFORT = qqFORT + log(mi**(-2)*Del)*log( - mx**(-1)*mg**(-1)*tg
cdel     + )*gs**2*Sn*CA*born*six**(-1) * ( 24 )
cdel      qqFORT = qqFORT + log(mi**(-2)*Del)*log( - mx**(-1)*mg**(-1)*ug
cdel     + )*gs**2*Sn*CA*born*six**(-1) * ( 24 )
csof      qqFORT = qqFORT + log(mu**(-2)*QF**2)*gs**2*Sn*CF*born*six**(-1)
csof     +  * (  - 36 )
      qqFORT = qqFORT + log(mu**(-2)*QR**2)*gs**2*Sn*CA*born*six**(-1)
     +  * ( 18 )
      qqFORT = qqFORT + log(mu**(-2)*QR**2)*gs**2*Sn*nf*born*six**(-1)
     +  * (  - 6 )
csof      qqFORT = qqFORT + log( - mx**(-1)*mg**(-1)*tg)*log( - mx**(-1)*
csof     + mg**(-1)*tg)*gs**2*Sn*CA*born*six**(-1) * ( 6 )
csof      qqFORT = qqFORT + log( - mx**(-1)*mg**(-1)*ug)*log( - mx**(-1)*
csof     + mg**(-1)*ug)*gs**2*Sn*CA*born*six**(-1) * ( 6 )
      qqFORT = qqFORT + log(mg**(-2)*QR**2)*gs**2*Sn*CA*born*six**(-1)
     +  * ( 4 )
csof      qqFORT = qqFORT + log( - xs)*gs**2*Sn*CA*born*rootlam**(-1)*
csof     + six**(-1) * ( 24*mi**2 + 24*mg**2 - 24*s )
csof      qqFORT = qqFORT + log( - xs)*log( - mi**(-2)*ti)*gs**2*Sn*CA*
csof     + born*six**(-1) * ( 12 )
csof      qqFORT = qqFORT + log( - xs)*log( - mi**(-2)*ui)*gs**2*Sn*CA*
csof     + born*six**(-1) * ( 12 )
csof      qqFORT = qqFORT + log( - xs)*log( - mx**(-1)*mg**(-1)*tg)*gs**2
csof     + *Sn*CA*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + log( - xs)*log( - mx**(-1)*mg**(-1)*ug)*gs**2
csof     + *Sn*CA*born*six**(-1) * (  - 12 )
csof      qqFORT = qqFORT + log( - xs)*log( - xs)*gs**2*Sn*CA*born*
csof     + six**(-1) * (  - 12 )

c  soft contribution 'grep csof'
      qqFORT = qqFORT + gs**2*Sn*CF*born*six**(-1) * (  - 36*zeta2 )
      qqFORT = qqFORT + Li2(1 - mi**(-2)*s**(-1)*ti*ui)*gs**2*Sn*CA*
     + born*six**(-1) * ( 12 )
      qqFORT = qqFORT + Li2(1 - mi**(-2)*s**(-1)*ti*ui)*gs**2*Sn*CF*
     + born*six**(-1) * (  - 24 )
      qqFORT = qqFORT + Li2(1 + mx**(-1)*mg*ti*ug**(-1)*xs**(-1))*gs**2
     + *Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + Li2(1 + mx**(-1)*mg*tg**(-1)*ui*xs**(-1))*gs**2
     + *Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + Li2(1 + mx*mg**(-1)*ti**(-1)*ug*xs**(-1))*gs**2
     + *Sn*CA*born*six**(-1) * ( 12 )
      qqFORT = qqFORT + Li2(1 + mx*mg**(-1)*tg*ui**(-1)*xs**(-1))*gs**2
     + *Sn*CA*born*six**(-1) * ( 12 )
      qqFORT = qqFORT + log( - mi**(-2)*ti)*log( - mi**(-2)*ti)*gs**2
     + *Sn*CA*born*six**(-1) * ( 6 )
      qqFORT = qqFORT + log( - mi**(-2)*ti)*log(mu**(-2)*QF**2)*gs**2
     + *Sn*CF*born*six**(-1) * ( 24 )
      qqFORT = qqFORT + log( - mi**(-2)*ti)*log( - mx**(-1)*mg**(-1)*
     + ug)*gs**2*Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + log( - mi**(-2)*ui)*log( - mi**(-2)*ui)*gs**2
     + *Sn*CA*born*six**(-1) * ( 6 )
      qqFORT = qqFORT + log( - mi**(-2)*ui)*log(mu**(-2)*QF**2)*gs**2
     + *Sn*CF*born*six**(-1) * ( 24 )
      qqFORT = qqFORT + log( - mi**(-2)*ui)*log( - mx**(-1)*mg**(-1)*
     + tg)*gs**2*Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + log(mu**(-2)*QF**2)*gs**2*Sn*CF*born*six**(-1)
     +  * (  - 36 )
      qqFORT = qqFORT + log( - mx**(-1)*mg**(-1)*tg)*log( - mx**(-1)*
     + mg**(-1)*tg)*gs**2*Sn*CA*born*six**(-1) * ( 6 )
      qqFORT = qqFORT + log( - mx**(-1)*mg**(-1)*ug)*log( - mx**(-1)*
     + mg**(-1)*ug)*gs**2*Sn*CA*born*six**(-1) * ( 6 )
      qqFORT = qqFORT + log( - xs)*gs**2*Sn*CA*born*rootlam**(-1)*
     + six**(-1) * ( 24*mi**2 + 24*mg**2 - 24*s )
      qqFORT = qqFORT + log( - xs)*log( - mi**(-2)*ti)*gs**2*Sn*CA*
     + born*six**(-1) * ( 12 )
      qqFORT = qqFORT + log( - xs)*log( - mi**(-2)*ui)*gs**2*Sn*CA*
     + born*six**(-1) * ( 12 )
      qqFORT = qqFORT + log( - xs)*log( - mx**(-1)*mg**(-1)*tg)*gs**2
     + *Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + log( - xs)*log( - mx**(-1)*mg**(-1)*ug)*gs**2
     + *Sn*CA*born*six**(-1) * (  - 12 )
      qqFORT = qqFORT + log( - xs)*log( - xs)*gs**2*Sn*CA*born*
     + six**(-1) * (  - 12 )

c               the phase space except for 1/s**2 
      qqFORT = qqFORT / ( 16.D0 * pi )

c               the averaging factors
      qqFORT = qqFORT /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      NG_QBV = qqFORT * (abs(mi)+mg)**2/4.D0

      end


