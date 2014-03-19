cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     HT_QGS(MASSIN,C)    SUSY CORRECTIONS                             c
c     HT_QGS1(MASSIN,C)   RESUMMATION EFFECT OF DELTA MB               c
c     HT_QGS2(MASSIN,C)   RESUMMATION OF DELTA MB ALONE                c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = t2                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = m2                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(20) = del                                               c
c       MASSIN(3-5,8,10,13-19,21-30) not needed                        c
c                                                                      c
c    ALL PHASE SPACE FACTORS INCOUDED TO GIVE                          c
c        s^2 d sig/(dtg ds4)                                           c
c                                                                      c
c    NEEDED MANDELSTAM VARIABLES :                                     c
c                                                                      c
c       Q(K1) + G(K2) -> T(P1) + H(P2) [+G(K3)]                        c
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
c    replacements in log file:                                         c
c                                                                      c
c       []     -> ()                                                   c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function HT_QGS1(massin,h_l,h_r,lumi,mu_tgb,alphas)

      implicit none 

      real*8     massin(1:30),lumi(1:3),h_l,h_r,Pi,Nc,Cf,alphas
     &          ,prefac,born,h_l2,h_r2,mb,xs2b,mu_tgb,alphas_b
     &          ,s,m1,m2,t1,t2,u1
     &          ,mg,msb1,msb2
     &          ,tag_delta_mb,tag_mb,B02
     &          ,susy_q
     &          ,MM_s
      complex*16 CSPEN

      external B02,CSPEN

      Pi     = 4.D0*atan(1.D0)
      Nc     = 3.D0
      Cf     = 4.D0/3.D0 
      prefac = 1.D0/(16.D0 * Pi**2)

      h_l2 = h_l**2
      h_r2 = h_r**2

      prefac = 1.D0/(16.D0*Pi**2)

      alphas_b = 1.D0

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      t2  = massin(2)
      m1  = massin(6)
      m2  = massin(7)
      mg   = massin(20)
      msb1 = massin(24)
      msb2 = massin(25)

c               born kinematics built in
      t1 = t2 + m2**2 - m1**2
      u1 = - s - t2
      
c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      susy_q   = 1.D0

c               avoid a possible mismatch with the running masses 
c               so make sure it cancels in the tag_mb fomula
      mb   = 2.5D0
      xs2b = 2.D0*mb*mu_tgb/(msb1**2-msb2**2)

c               replace B_fin with B02
      tag_mb = 0.D0

      tag_mb = tag_mb + Cf*xs2b*mg*mb**(-1)*Pi*alphas*prefac * (
     +     - 4*B02(0.D0,msb1,mg,m1**2)
     +     + 4*B02(0.D0,msb2,mg,m1**2)
     +     )

      tag_delta_mb = 1.D0/(1.D0+tag_mb)**2
ctp      tag_delta_mb = (1.D0-tag_mb)**2
ctp      tag_delta_mb = 1.D0-2.D0*tag_mb

c               want to compute K_SUSY as a correction to 2HDM
      tag_delta_mb = tag_delta_mb - 1.D0

      born =
     +  + (h_r2+h_l2)*Nc*Cf*Pi*alphas_b * (  - 8 - 8*s**(-1)*t1 - 8*
     +    s**(-1)*t1**2*u1**(-1) - 4*s**(-1)*u1 + 8*s*m1**2*u1**(-2) - 
     +    4*s*u1**(-1) + 8*m1**2*t1*u1**(-2) + 8*m1**2*u1**(-1) - 8*t1*
     +    u1**(-1) )

      MM_s =
     +  + susy_q * (  tag_delta_mb*born )

c               the phase space except for 1/s**2 
      HT_QGS1 = MM_s / ( 16.D0 * pi )

c               the averaging factors
      HT_QGS1 = HT_QGS1 /4.D0 /Nc /(Nc**2-1.D0)

c               the luminosity
      HT_QGS1 = HT_QGS1 *lumi(1)

      return
      end

c --------------------------------------------------------------------
c this is without the resummation, should be called either
c as a function of mu_tgb or as a function of (-Ab)
c dependent if resummable or remainder corrections are wanted
ctp sign of that one checked 11/8/04
      real*8 function HT_QGS2(massin,h_l,h_r,lumi,mu_tgb,alphas)

      implicit none 
      
      real*8     massin(1:30),lumi(1:3),h_l,h_r,Pi,Nc,Cf,alphas
     &          ,prefac,born,h_l2,h_r2,mb,xs2b,mu_tgb,alphas_b
     &          ,s,m1,m2,t1,t2,u1
     &          ,mg,msb1,msb2
     &          ,tag_delta_mb,tag_mb,B02
     &          ,susy_q
     &          ,MM_s
      complex*16 CSPEN
      logical lcheck_deltab
      parameter ( lcheck_deltab = .false. )

      external B02,CSPEN

      Pi     = 4.D0*atan(1.D0)
      Nc     = 3.D0
      Cf     = 4.D0/3.D0 
      prefac = 1.D0/(16.D0 * Pi**2)

      h_l2 = h_l**2
      h_r2 = h_r**2

      prefac = 1.D0/(16.D0*Pi**2)

      alphas_b = 1.D0

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      t2  = massin(2)
      m1  = massin(6)
      m2  = massin(7)
      mg   = massin(20)
      msb1 = massin(24)
      msb2 = massin(25)

c               born kinematics built in
      t1 = t2 + m2**2 - m1**2
      u1 = - s - t2
      
c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      susy_q   = 1.D0

c               avoid a possible mismatch with the running masses 
c               so make sure it cancels in the tag_mb fomula
      mb   = 2.5D0
      xs2b = 2.D0*mb*mu_tgb/(msb1**2-msb2**2)

      tag_mb = 0.D0

c               replace B_fin with B02
c               mass dimensions with 1/mb is m^0
      tag_mb = tag_mb + Cf*xs2b*mg*mb**(-1)*Pi*alphas*prefac * (
     +     - 4*B02(0.D0,msb1,mg,m1**2)
     +     + 4*B02(0.D0,msb2,mg,m1**2)
     +     )

      if ( lcheck_deltab ) then 
         print*, " mb      = ",mb
         print*, " mg      = ",mg
         print*, " msb1    = ",msb1
         print*, " msb2    = ",msb2
         print*, " xs2b    = ",xs2b
         print*, " alphas  = ",alphas
         print*, " B1 - B2 = ",B02(0.D0,msb1,mg,m1**2)
     +                       - B02(0.D0,msb2,mg,m1**2)
         print*, " factors = ",Cf*Pi*prefac 
         print*, " Delta_b = ",tag_mb
         stop
      end if 

c               assume cross section to be proportional to yb**2
      tag_delta_mb = (1.D0-tag_mb)**2

c               want to just add SUSY part to 2HDM
      tag_delta_mb = tag_delta_mb - 1

      born =
     +  + (h_r2+h_l2)*Nc*Cf*Pi*alphas_b * (  - 8 - 8*s**(-1)*t1 - 8*
     +    s**(-1)*t1**2*u1**(-1) - 4*s**(-1)*u1 + 8*s*m1**2*u1**(-2) - 
     +    4*s*u1**(-1) + 8*m1**2*t1*u1**(-2) + 8*m1**2*u1**(-1) - 8*t1*
     +    u1**(-1) )

      MM_s =
     +  + susy_q * (  tag_delta_mb*born )

c               the phase space except for 1/s**2 
      HT_QGS2 = MM_s / ( 16.D0 * pi )

c               the averaging factors
      HT_QGS2 = HT_QGS2 /4.D0 /Nc /(Nc**2-1.D0)

c               the luminosity
      HT_QGS2 = HT_QGS2 *lumi(1)

      return
      end

c --------------------------------------------------------------------
      real*8 function HT_QGS(massin,h_l,h_r,h_s,msb,mst,lumi)

      implicit none 
      
      real*8     massin(1:30),lumi(1:3),h_l,h_r,Pi,Nc,Cf,Ns,alphas
     &          ,prefac,kaellen,born,qr,log_qr,h_l2,h_r2
     &          ,s,m1,m2,t1,t2,u1,u2,u
     &          ,mg,msx,mst1,mst2,msb1,msb2,msusy
     &          ,sb2,cb2,st2,ct2,sb,st,cb,ct,s2b,c2b,s2t,c2t
     &          ,tag_delta_mb,h_s(1:2,1:2),B02,msb(2,2),mst(2,2)
     &          ,SCA(1:10,1:10),SCB(1:10,1:10),SCBP(1:10,1:10)
     &          ,SCC(1:10,1:10),SCD(1:10,1:10)
     &          ,susy_g,susy_q,susy_t,susy_qi,susy_ti
     &          ,susy_qqg,susy_ttg,susy_bth,susy_box
     &          ,MM_s
      complex*16 CSPEN

      external B02,CSPEN

      Pi     = 4.D0*atan(1.D0)
      Nc     = 3.D0
      Ns     = 6.D0
      Cf     = 4.D0/3.D0 
      prefac = 1.D0/(16.D0 * Pi**2)

      h_l2 = h_l**2
      h_r2 = h_r**2

      prefac = 1.D0/(16.D0*Pi**2)

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      t2  = massin(2)
      m1  = massin(6)
      m2  = massin(7)
      mg   = massin(20)
      msx  = massin(21)
      mst1 = massin(22)
      mst2 = massin(23)
      msb1 = massin(24)
      msb2 = massin(25)

c      print*, " h_l ",h_l
c      print*, " h_r ",h_r
c      print*, " h_s ",h_s
c      do i=6,7,1
c         print*, " HT_QGS ",i,massin(i)
c      end do
c      do i=20,29,1
c         print*, " HT_QGS ",i,massin(i)
c      end do

      st = mst(1,2)
      ct = mst(1,1)
      st2 = mst(1,2)**2 
      ct2 = mst(1,1)**2 
      s2t = 2.D0*mst(1,1)*mst(1,2)
      c2t = mst(1,1)**2 - mst(1,2)**2

      sb = msb(1,2)
      cb = msb(1,1)
      sb2 = msb(1,2)**2
      cb2 = msb(1,1)**2 
      s2b = 2.D0*msb(1,1)*mst(1,2)
      c2b = msb(1,1)**2 - msb(1,2)**2

c               born kinematics built in
      t1 = t2 + m2**2 - m1**2
      u1 = - s - t2
      u2 = u1 + m1**2 - m2**2 
      
      u = u2 + m2**2 

      kaellen = s**2 +m1**4 +m2**4 - 2*( s*m1**2+s*m2**2+m1**2*m2**2 )

c               the factorization/renormalization scale 
      qr = massin(12)
      log_qr = log(qr**2/m1**2)

c               the scalar functions 
      call SCALAR_ARRAY_HT_SUSY(massin,SCA,SCB,SCBP,SCC,SCD)

c               alpha_s is cut, re-appears as nlo later 
      alphas   = 1.D0
      susy_g   = 1.D0
      susy_q   = 1.D0
      susy_t   = 1.D0
      susy_qi  = 1.D0
      susy_ti  = 1.D0
      susy_qqg = 1.D0
      susy_ttg = 1.D0
      susy_bth = 1.D0
      susy_box = 1.D0

c               remove that one 
      tag_delta_mb = 0.D0

      born =
     +  + (h_r2+h_l2)*Nc*Cf*Pi*alphas * (  - 8 - 8*s**(-1)*t1 - 8*
     +    s**(-1)*t1**2*u1**(-1) - 4*s**(-1)*u1 + 8*s*m1**2*u1**(-2) - 
     +    4*s*u1**(-1) + 8*m1**2*t1*u1**(-2) + 8*m1**2*u1**(-1) - 8*t1*
     +    u1**(-1) )

      MM_s =
     +  + susy_q * (  - tag_delta_mb*born )
      MM_s = MM_s + (h_r2+h_l2)*Nc*Cf**2*Pi**2*alphas**2*prefac * (  - 
     +    64 - 32*s**(-1)*u1 - 32*s*u1**(-1) )
      MM_s = MM_s + (h_r2+h_l2)*Nc**2*Cf*Pi**2*alphas**2*prefac * ( 64
     +     - 16*m1**2*s*u1**(-2) - 16*m1**2*u1**(-1) + 16*s**(-1)*t1 + 
     +    32*s**(-1)*u1 + 32*s*u1**(-1) + 16*t1*u1**(-1) )
      MM_s = MM_s + (h_r2+h_l2)*Cf*Pi**2*alphas**2*prefac * (  - 16*
     +    m1**2*s*u1**(-2) - 16*m1**2*u1**(-1) + 16*s**(-1)*t1 + 16*t1*
     +    u1**(-1) )
      MM_s = MM_s + Nc*Pi*alphas*prefac * (  - 8./3.*born*log_qr )
      MM_s = MM_s + Cf*Pi*alphas*prefac * (  - 8*born*log_qr )
      MM_s = MM_s + Pi*alphas*prefac * (  - 4./3.*(Ns-2)*born*log_qr - 
     +    8./3.*born*log_qr )
      MM_s = MM_s + SCA(2,1)*susy_t*Cf*Pi*alphas*prefac * (  - 4*
     +    m1**(-2)*born )
      MM_s = MM_s + SCA(2,1)*susy_qi*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*s**(-1)*u1**(-1) - 32*s**(-2)*t1
     +     - 32*s**(-2)*t1**2*u1**(-1) - 32*s**(-2)*u1 - 32*s**(-1)*t1*
     +    u1**(-1) - 32*s**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 64*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) - 
     +    64*m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) - 128*m1**2*s*
     +    u1**(-3) + 64*m1**2*s*u1**(-2)*u**(-1) - 128*m1**2*t1*
     +    u1**(-3) + 32*m1**2*t1*u1**(-2)*u**(-1) - 64*m1**2*u1**(-2)
     +     + 128*m1**4*s*u1**(-3)*u**(-1) + 128*m1**4*t1*u1**(-3)*
     +    u**(-1) + 64*m1**4*u1**(-2)*u**(-1) + 64*s**(-1)*t1*u1**(-1)
     +     - 32*s**(-1)*t1*u**(-1) + 64*s**(-1)*t1**2*u1**(-2) - 32*
     +    s**(-1)*t1**2*u1**(-1)*u**(-1) - 32*s*u1**(-1)*u**(-1) + 32*
     +    t1*u1**(-2) - 32*t1*u1**(-1)*u**(-1) + 32*u1**(-1) - 32*
     +    u**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 16*(h_r2*cb2+h_l2*sb2)*s**(-2)*t1 - 24*
     +    (h_r2*cb2+h_l2*sb2)*s**(-2)*t1**2*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*s**(-2)*u1 - 8*(h_r2*cb2+h_l2*sb2)*
     +    s**(-1)*t1*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*s**(-2)*t1 - 24*
     +    (h_r2*sb2+h_l2*cb2)*s**(-2)*t1**2*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*s**(-2)*u1 - 8*(h_r2*sb2+h_l2*cb2)*
     +    s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_qqg*Cf*Pi**2*alphas**2*prefac * ( 16*
     +    (h_r2*cb2+h_l2*sb2)*s**(-2)*t1 + 24*(h_r2*cb2+h_l2*sb2)*
     +    s**(-2)*t1**2*u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*s**(-2)*u1 + 
     +    8*(h_r2*cb2+h_l2*sb2)*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*s**(-2)*t1 + 24*(h_r2*sb2+h_l2*cb2)*
     +    s**(-2)*t1**2*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*s**(-2)*u1 + 
     +    8*(h_r2*sb2+h_l2*cb2)*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1)
     +     + 16*m1**2*s*u1**(-2)*u**(-1) + 32*m1**2*t1*u1**(-2)*u**(-1)
     +     - 16*s**(-1)*t1*u1**(-1) - 16*s**(-1)*t1**2*u1**(-1)*u**(-1)
     +     + 32*s*u1**(-2) - 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-1)*
     +    u**(-1) + 16*u1**(-1) - 16*u**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) - 
     +    16*m1**2*s*u1**(-2)*u**(-1) - 32*m1**2*t1*u1**(-2)*u**(-1) + 
     +    16*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1**2*u1**(-1)*u**(-1) - 
     +    32*s*u1**(-2) + 16*s*u1**(-1)*u**(-1) + 16*t1*u1**(-1)*
     +    u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1)
     +     + 16*m1**2*s*u1**(-2)*u**(-1) + 32*m1**2*t1*u1**(-2)*u**(-1)
     +     - 16*s**(-1)*t1*u1**(-1) - 16*s**(-1)*t1**2*u1**(-1)*u**(-1)
     +     + 32*s*u1**(-2) - 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-1)*
     +    u**(-1) + 16*u1**(-1) - 16*u**(-1) )
      MM_s = MM_s + SCA(2,1)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) - 
     +    16*m1**2*s*u1**(-2)*u**(-1) - 32*m1**2*t1*u1**(-2)*u**(-1) + 
     +    16*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1**2*u1**(-1)*u**(-1) - 
     +    32*s*u1**(-2) + 16*s*u1**(-1)*u**(-1) + 16*t1*u1**(-1)*
     +    u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(3,1)*susy_t*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2t*m1**(-2)*s**(-1)*t1 - 16*c2t*
     +    m1**(-2)*s**(-1)*t1**2*u1**(-1) - 8*c2t*m1**(-2)*s**(-1)*u1
     +     - 8*c2t*m1**(-2)*s*u1**(-1) - 16*c2t*m1**(-2)*t1*u1**(-1) - 
     +    16*c2t*m1**(-2) - 16*c2t*s*u1**(-2) )
      MM_s = MM_s + SCA(3,1)*susy_t*Cf*Pi*alphas*prefac * ( 2*m1**(-2)*
     +    born )
      MM_s = MM_s + SCA(3,1)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*s*u1**(-3) - 
     +    32*m1**2*s*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3) - 16*m1**2
     +    *t1*u1**(-2)*u**(-1) + 32*m1**2*u1**(-2) - 64*m1**4*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1) - 32*m1**4*
     +    u1**(-2)*u**(-1) - 32*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1*
     +    u**(-1) - 32*s**(-1)*t1**2*u1**(-2) + 16*s**(-1)*t1**2*
     +    u1**(-1)*u**(-1) + 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-2) + 
     +    16*t1*u1**(-1)*u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(3,1)*susy_ti*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2t*m1**2*t1*u1**(-2)*u**(-1) - 16*c2t*
     +    s**(-1)*t1*u**(-1) - 16*c2t*s**(-1)*t1**2*u1**(-1)*u**(-1) - 
     +    16*c2t*s*u1**(-1)*u**(-1) - 16*c2t*t1*u1**(-1)*u**(-1) - 16*
     +    c2t*u**(-1) )
      MM_s = MM_s + SCA(3,1)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 16*m1**2*s*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3)
     +     - 48*m1**2*t1*u1**(-2)*u**(-1) + 16*m1**2*u1**(-1)*u**(-1)
     +     - 32*m1**4*s*u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1)
     +     - 16*s**(-1)*t1**2*u1**(-2) + 16*s**(-1)*t1**2*u1**(-1)*
     +    u**(-1) - 32*s*u1**(-2) + 16*s*u1**(-1)*u**(-1) + 16*t1*
     +    u1**(-1)*u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(3,1)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1)
     +     + 16*m1**2*s*u1**(-2)*u**(-1) - 64*m1**2*t1*u1**(-3) + 48*
     +    m1**2*t1*u1**(-2)*u**(-1) - 16*m1**2*u1**(-1)*u**(-1) + 32*
     +    m1**4*s*u1**(-3)*u**(-1) + 64*m1**4*t1*u1**(-3)*u**(-1) + 16*
     +    s**(-1)*t1**2*u1**(-2) - 16*s**(-1)*t1**2*u1**(-1)*u**(-1) + 
     +    32*s*u1**(-2) - 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-1)*
     +    u**(-1) + 16*u1**(-1) - 16*u**(-1) )
      MM_s = MM_s + SCA(3,1)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1*u1**(-1)*
     +    u**(-1) - 16*m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) - 64*m1**2*
     +    t1*u1**(-3) + 16*m1**2*t1*u1**(-2)*u**(-1) - 16*m1**2*
     +    u1**(-1)*u**(-1) + 32*m1**4*s*u1**(-3)*u**(-1) + 64*m1**4*t1*
     +    u1**(-3)*u**(-1) + 16*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1**2*
     +    u1**(-2) )
      MM_s = MM_s + SCA(3,1)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) + 16*
     +    m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3)
     +     - 16*m1**2*t1*u1**(-2)*u**(-1) + 16*m1**2*u1**(-1)*u**(-1)
     +     - 32*m1**4*s*u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1)
     +     - 16*s**(-1)*t1*u1**(-1) - 16*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCA(3,2)*susy_t*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2t*m1**(-2)*s**(-1)*t1 + 16*c2t*
     +    m1**(-2)*s**(-1)*t1**2*u1**(-1) + 8*c2t*m1**(-2)*s**(-1)*u1
     +     + 8*c2t*m1**(-2)*s*u1**(-1) + 16*c2t*m1**(-2)*t1*u1**(-1) + 
     +    16*c2t*m1**(-2) + 16*c2t*s*u1**(-2) )
      MM_s = MM_s + SCA(3,2)*susy_t*Cf*Pi*alphas*prefac * ( 2*m1**(-2)*
     +    born )
      MM_s = MM_s + SCA(3,2)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*s*u1**(-3) - 
     +    32*m1**2*s*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3) - 16*m1**2
     +    *t1*u1**(-2)*u**(-1) + 32*m1**2*u1**(-2) - 64*m1**4*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1) - 32*m1**4*
     +    u1**(-2)*u**(-1) - 32*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1*
     +    u**(-1) - 32*s**(-1)*t1**2*u1**(-2) + 16*s**(-1)*t1**2*
     +    u1**(-1)*u**(-1) + 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-2) + 
     +    16*t1*u1**(-1)*u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(3,2)*susy_ti*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2t*m1**2*t1*u1**(-2)*u**(-1) + 16*
     +    c2t*s**(-1)*t1*u**(-1) + 16*c2t*s**(-1)*t1**2*u1**(-1)*
     +    u**(-1) + 16*c2t*s*u1**(-1)*u**(-1) + 16*c2t*t1*u1**(-1)*
     +    u**(-1) + 16*c2t*u**(-1) )
      MM_s = MM_s + SCA(3,2)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1*u1**(-1)*
     +    u**(-1) - 16*m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) - 64*m1**2*
     +    t1*u1**(-3) + 16*m1**2*t1*u1**(-2)*u**(-1) - 16*m1**2*
     +    u1**(-1)*u**(-1) + 32*m1**4*s*u1**(-3)*u**(-1) + 64*m1**4*t1*
     +    u1**(-3)*u**(-1) + 16*s**(-1)*t1*u1**(-1) + 16*s**(-1)*t1**2*
     +    u1**(-2) )
      MM_s = MM_s + SCA(3,2)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1)*u**(-1) + 16*
     +    m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3)
     +     - 16*m1**2*t1*u1**(-2)*u**(-1) + 16*m1**2*u1**(-1)*u**(-1)
     +     - 32*m1**4*s*u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1)
     +     - 16*s**(-1)*t1*u1**(-1) - 16*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCA(3,2)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 16*m1**2*s*u1**(-2)*u**(-1) + 64*m1**2*t1*u1**(-3)
     +     - 48*m1**2*t1*u1**(-2)*u**(-1) + 16*m1**2*u1**(-1)*u**(-1)
     +     - 32*m1**4*s*u1**(-3)*u**(-1) - 64*m1**4*t1*u1**(-3)*u**(-1)
     +     - 16*s**(-1)*t1**2*u1**(-2) + 16*s**(-1)*t1**2*u1**(-1)*
     +    u**(-1) - 32*s*u1**(-2) + 16*s*u1**(-1)*u**(-1) + 16*t1*
     +    u1**(-1)*u**(-1) - 16*u1**(-1) + 16*u**(-1) )
      MM_s = MM_s + SCA(3,2)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*t1**2*u1**(-2)*u**(-1)
     +     + 16*m1**2*s*u1**(-2)*u**(-1) - 64*m1**2*t1*u1**(-3) + 48*
     +    m1**2*t1*u1**(-2)*u**(-1) - 16*m1**2*u1**(-1)*u**(-1) + 32*
     +    m1**4*s*u1**(-3)*u**(-1) + 64*m1**4*t1*u1**(-3)*u**(-1) + 16*
     +    s**(-1)*t1**2*u1**(-2) - 16*s**(-1)*t1**2*u1**(-1)*u**(-1) + 
     +    32*s*u1**(-2) - 16*s*u1**(-1)*u**(-1) - 16*t1*u1**(-1)*
     +    u**(-1) + 16*u1**(-1) - 16*u**(-1) )
      MM_s = MM_s + SCA(4,1)*susy_qi*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*u1**(-1) + 16*s**(-2)*
     +    t1 + 16*s**(-2)*t1**2*u1**(-1) + 16*s**(-2)*u1 + 16*s**(-1)*
     +    t1*u1**(-1) + 16*s**(-1) )
      MM_s = MM_s + SCA(4,1)*susy_qi*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2b*m1**2*s**(-1)*u1**(-1) + 16*c2b*
     +    s**(-2)*t1 + 16*c2b*s**(-2)*t1**2*u1**(-1) + 16*c2b*s**(-2)*
     +    u1 + 16*c2b*s**(-1)*t1*u1**(-1) + 16*c2b*s**(-1) )
      MM_s = MM_s + SCA(4,1)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * ( 16*(h_r2*sb2+h_l2*cb2)*s**(-2)*t1 + 24*(h_r2*sb2+h_l2*cb2)*
     +    s**(-2)*t1**2*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*s**(-2)*u1 + 
     +    8*(h_r2*sb2+h_l2*cb2)*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCA(4,1)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    16*(h_r2*sb2+h_l2*cb2)*s**(-2)*t1 - 24*(h_r2*sb2+h_l2*cb2)*
     +    s**(-2)*t1**2*u1**(-1) + 16*(h_r2*sb2+h_l2*cb2)*s**(-2)*u1 - 
     +    8*(h_r2*sb2+h_l2*cb2)*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCA(4,2)*susy_qi*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*s**(-1)*u1**(-1) + 16*s**(-2)*
     +    t1 + 16*s**(-2)*t1**2*u1**(-1) + 16*s**(-2)*u1 + 16*s**(-1)*
     +    t1*u1**(-1) + 16*s**(-1) )
      MM_s = MM_s + SCA(4,2)*susy_qi*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2b*m1**2*s**(-1)*u1**(-1) - 16*c2b*
     +    s**(-2)*t1 - 16*c2b*s**(-2)*t1**2*u1**(-1) - 16*c2b*s**(-2)*
     +    u1 - 16*c2b*s**(-1)*t1*u1**(-1) - 16*c2b*s**(-1) )
      MM_s = MM_s + SCA(4,2)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * ( 16*(h_r2*cb2+h_l2*sb2)*s**(-2)*t1 + 24*(h_r2*cb2+h_l2*sb2)*
     +    s**(-2)*t1**2*u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*s**(-2)*u1 + 
     +    8*(h_r2*cb2+h_l2*sb2)*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCA(4,2)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    16*(h_r2*cb2+h_l2*sb2)*s**(-2)*t1 - 24*(h_r2*cb2+h_l2*sb2)*
     +    s**(-2)*t1**2*u1**(-1) + 16*(h_r2*cb2+h_l2*sb2)*s**(-2)*u1 - 
     +    8*(h_r2*cb2+h_l2*sb2)*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(1,2)*susy_g*Nc*Pi*alphas*prefac * (  - 8./3.*
     +    born )
      MM_s = MM_s + SCB(1,3)*susy_g*Pi*alphas*prefac * (  - 4./3.*
     +    (Ns-2)*born )
      MM_s = MM_s + SCB(1,4)*susy_g*Pi*alphas*prefac * (  - 2./3.*born
     +     )
      MM_s = MM_s + SCB(1,5)*susy_g*Pi*alphas*prefac * (  - 2./3.*born
     +     )
      MM_s = MM_s + SCB(1,6)*susy_g*Pi*alphas*prefac * (  - 2./3.*born
     +     )
      MM_s = MM_s + SCB(1,6)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*m1**2*t1*u1**(-2) + 32*m1**2*u1**(-1) - 
     +    64*m1**4*s*u1**(-3) - 32*s**(-1)*t1 - 32*s**(-1)*t1**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(1,6)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*s*u1**(-2) - 32*m1**2*t1*u1**(-2)
     +     + 64*m1**4*s*u1**(-3) )
      MM_s = MM_s + SCB(1,6)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*cb*m1*s*t2**(-2) + 32*h_l*ct*cb*m1*s**2*
     +    u1**(-1)*t2**(-2) - 32*h_l*ct*cb*m1*u1**(-1) - 32*h_l*ct*cb*
     +    m1*t2**(-1) + 32*h_r*st*sb*m1*s*t2**(-2) + 32*h_r*st*sb*m1*
     +    s**2*u1**(-1)*t2**(-2) - 32*h_r*st*sb*m1*u1**(-1) - 32*h_r*st
     +    *sb*m1*t2**(-1) )
      MM_s = MM_s + SCB(1,6)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*sb*m1*s*t2**(-2) - 32*h_l*ct*sb*m1*s**2*
     +    u1**(-1)*t2**(-2) + 32*h_l*ct*sb*m1*u1**(-1) + 32*h_l*ct*sb*
     +    m1*t2**(-1) + 32*h_r*st*cb*m1*s*t2**(-2) + 32*h_r*st*cb*m1*
     +    s**2*u1**(-1)*t2**(-2) - 32*h_r*st*cb*m1*u1**(-1) - 32*h_r*st
     +    *cb*m1*t2**(-1) )
      MM_s = MM_s + SCB(1,7)*susy_g*Pi*alphas*prefac * (  - 2./3.*born
     +     )
      MM_s = MM_s + SCB(1,7)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*s*u1**(-2) - 32*m1**2*t1*u1**(-2)
     +     + 64*m1**4*s*u1**(-3) )
      MM_s = MM_s + SCB(1,7)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*m1**2*t1*u1**(-2) + 32*m1**2*u1**(-1) - 
     +    64*m1**4*s*u1**(-3) - 32*s**(-1)*t1 - 32*s**(-1)*t1**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(1,7)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*s*t2**(-2) - 32*h_l*st*cb*m1*s**2*
     +    u1**(-1)*t2**(-2) + 32*h_l*st*cb*m1*u1**(-1) + 32*h_l*st*cb*
     +    m1*t2**(-1) + 32*h_r*ct*sb*m1*s*t2**(-2) + 32*h_r*ct*sb*m1*
     +    s**2*u1**(-1)*t2**(-2) - 32*h_r*ct*sb*m1*u1**(-1) - 32*h_r*ct
     +    *sb*m1*t2**(-1) )
      MM_s = MM_s + SCB(1,7)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*s*t2**(-2) + 32*h_l*st*sb*m1*s**2*
     +    u1**(-1)*t2**(-2) - 32*h_l*st*sb*m1*u1**(-1) - 32*h_l*st*sb*
     +    m1*t2**(-1) + 32*h_r*ct*cb*m1*s*t2**(-2) + 32*h_r*ct*cb*m1*
     +    s**2*u1**(-1)*t2**(-2) - 32*h_r*ct*cb*m1*u1**(-1) - 32*h_r*ct
     +    *cb*m1*t2**(-1) )
      MM_s = MM_s + SCB(1,8)*susy_q*(h_r2+h_l2)**(-1)*(h_r2-h_l2)*Cf*Pi
     + *alphas*prefac * (  - 4*c2b*born )
      MM_s = MM_s + SCB(1,8)*susy_q*Cf*Pi*alphas*prefac * (  - 2*born )
      MM_s = MM_s + SCB(1,8)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * ( 16*(h_r2*sb2+h_l2*cb2)*m1**2*mg**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*msb1**2*s**(-1)*u1**(-1) + 40*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*t1**2*u1**(-1) - 32*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*u1 + 24*(h_r2*sb2+h_l2*cb2)
     +    *mg**2*s**(-1)*t1*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*mg**2*
     +    s**(-1) - 40*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*t1**2*
     +    u1**(-1) + 32*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*u1 - 24*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1) + 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(1,8)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    16*(h_r2*sb2+h_l2*cb2)*m1**2*mg**2*s**(-1)*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*msb1**2*s**(-1)*u1**(-1) - 40*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*t1**2*u1**(-1) + 32*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*u1 - 24*(h_r2*sb2+h_l2*cb2)
     +    *mg**2*s**(-1)*t1*u1**(-1) + 16*(h_r2*sb2+h_l2*cb2)*mg**2*
     +    s**(-1) + 40*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*t1**2*
     +    u1**(-1) - 32*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*u1 + 24*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1) - 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(1,9)*susy_q*(h_r2+h_l2)**(-1)*(h_r2-h_l2)*Cf*Pi
     + *alphas*prefac * ( 4*c2b*born )
      MM_s = MM_s + SCB(1,9)*susy_q*Cf*Pi*alphas*prefac * (  - 2*born )
      MM_s = MM_s + SCB(1,9)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * ( 16*(h_r2*cb2+h_l2*sb2)*m1**2*mg**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*msb2**2*s**(-1)*u1**(-1) + 40*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*t1**2*u1**(-1) - 32*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*u1 + 24*(h_r2*cb2+h_l2*sb2)
     +    *mg**2*s**(-1)*t1*u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*mg**2*
     +    s**(-1) - 40*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*t1**2*
     +    u1**(-1) + 32*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*u1 - 24*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1) - 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(1,9)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    16*(h_r2*cb2+h_l2*sb2)*m1**2*mg**2*s**(-1)*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*msb2**2*s**(-1)*u1**(-1) - 40*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*t1**2*u1**(-1) + 32*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*u1 - 24*(h_r2*cb2+h_l2*sb2)
     +    *mg**2*s**(-1)*t1*u1**(-1) + 16*(h_r2*cb2+h_l2*sb2)*mg**2*
     +    s**(-1) + 40*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*t1**2*
     +    u1**(-1) - 32*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*u1 + 24*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1) + 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(2,2)*susy_t*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2t*m1**(-2)*mg**2*s**(-1)*t1 - 16*
     +    c2t*m1**(-2)*mg**2*s**(-1)*t1**2*u1**(-1) - 8*c2t*m1**(-2)*
     +    mg**2*s**(-1)*u1 - 8*c2t*m1**(-2)*mg**2*s*u1**(-1) - 16*c2t*
     +    m1**(-2)*mg**2*t1*u1**(-1) - 16*c2t*m1**(-2)*mg**2 + 16*c2t*
     +    m1**(-2)*mst1**2*s**(-1)*t1 + 16*c2t*m1**(-2)*mst1**2*s**(-1)
     +    *t1**2*u1**(-1) + 8*c2t*m1**(-2)*mst1**2*s**(-1)*u1 + 8*c2t*
     +    m1**(-2)*mst1**2*s*u1**(-1) + 16*c2t*m1**(-2)*mst1**2*t1*
     +    u1**(-1) + 16*c2t*m1**(-2)*mst1**2 - 16*c2t*m1**2*s*u1**(-2)
     +     - 16*c2t*mg**2*s*u1**(-2) + 16*c2t*mst1**2*s*u1**(-2) - 16*
     +    c2t*s**(-1)*t1 - 16*c2t*s**(-1)*t1**2*u1**(-1) - 8*c2t*
     +    s**(-1)*u1 - 8*c2t*s*u1**(-1) - 16*c2t*t1*u1**(-1) - 16*c2t )
      MM_s = MM_s + SCB(2,2)*susy_t*Cf*Pi*alphas*prefac * ( 2*m1**(-2)*
     +    mg**2*born - 2*m1**(-2)*mst1**2*born - 2*born )
      MM_s = MM_s + SCB(2,2)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 64*m1**2*mg**2*s*u1**(-3) + 64*m1**2*mg**2*
     +    t1*u1**(-3) + 32*m1**2*mg**2*u1**(-2) - 64*m1**2*mst1**2*s*
     +    u1**(-3) - 64*m1**2*mst1**2*t1*u1**(-3) - 32*m1**2*mst1**2*
     +    u1**(-2) - 32*m1**2*s**(-1)*t1*u1**(-1) - 32*m1**2*s**(-1)*
     +    t1**2*u1**(-2) - 16*m1**2*t1*u1**(-2) - 16*m1**2*u1**(-1) + 
     +    64*m1**4*s*u1**(-3) + 64*m1**4*t1*u1**(-3) + 32*m1**4*
     +    u1**(-2) - 32*mg**2*s**(-1)*t1*u1**(-1) - 32*mg**2*s**(-1)*
     +    t1**2*u1**(-2) - 16*mg**2*t1*u1**(-2) - 16*mg**2*u1**(-1) + 
     +    32*mst1**2*s**(-1)*t1*u1**(-1) + 32*mst1**2*s**(-1)*t1**2*
     +    u1**(-2) + 16*mst1**2*t1*u1**(-2) + 16*mst1**2*u1**(-1) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) - 32*
     +    s2t*m1*mg*s**(-1)*t1**2*u1**(-2) + 16*s2t*m1*mg*s*u1**(-2) - 
     +    16*s2t*m1*mg*t1*u1**(-2) + 64*s2t*m1**3*mg*s*u1**(-3) + 64*
     +    s2t*m1**3*mg*t1*u1**(-3) + 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * ( 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) + 32*s2t*m1*mg*
     +    s**(-1)*t1**2*u1**(-2) - 16*s2t*m1*mg*s*u1**(-2) + 16*s2t*m1*
     +    mg*t1*u1**(-2) - 64*s2t*m1**3*mg*s*u1**(-3) - 64*s2t*m1**3*mg
     +    *t1*u1**(-3) - 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 48*m1**2*mg**2*s*u1**(-3) - 64*
     +    m1**2*mg**2*t1*u1**(-3) - 32*m1**2*mg**2*u1**(-2) + 48*m1**2*
     +    mst1**2*s*u1**(-3) + 64*m1**2*mst1**2*t1*u1**(-3) + 32*m1**2*
     +    mst1**2*u1**(-2) + 16*m1**2*s**(-1)*t1*u1**(-1) + 16*m1**2*
     +    s**(-1)*t1**2*u1**(-2) - 32*m1**2*s*u1**(-2) - 32*m1**2*t1*
     +    u1**(-2) - 16*m1**2*u1**(-1) + 32*m1**4*mg**2*s*u1**(-4) - 32
     +    *m1**4*mst1**2*s*u1**(-4) - 16*m1**4*s*u1**(-3) - 64*m1**4*t1
     +    *u1**(-3) - 32*m1**4*u1**(-2) + 32*m1**6*s*u1**(-4) + 16*
     +    mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*s**(-1)*t1**2*u1**(-2)
     +     - 32*mg**2*s*u1**(-2) - 16*mg**2*u1**(-1) - 16*mst1**2*
     +    s**(-1)*t1*u1**(-1) - 16*mst1**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    mst1**2*s*u1**(-2) + 16*mst1**2*u1**(-1) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 48*m1**2*mg**2*s*u1**(-3) + 64*m1**2*mg**2*
     +    t1*u1**(-3) + 32*m1**2*mg**2*u1**(-2) - 48*m1**2*mst1**2*s*
     +    u1**(-3) - 64*m1**2*mst1**2*t1*u1**(-3) - 32*m1**2*mst1**2*
     +    u1**(-2) - 48*m1**2*s**(-1)*t1*u1**(-1) - 48*m1**2*s**(-1)*
     +    t1**2*u1**(-2) - 32*m1**2*s*u1**(-2) - 16*m1**2*u1**(-1) - 32
     +    *m1**4*mg**2*s*u1**(-4) + 32*m1**4*mst1**2*s*u1**(-4) - 48*
     +    m1**4*s*u1**(-3) + 128*m1**4*t1*u1**(-3) + 32*m1**4*u1**(-2)
     +     - 96*m1**6*s*u1**(-4) - 16*mg**2*s**(-1)*t1*u1**(-1) - 16*
     +    mg**2*s**(-1)*t1**2*u1**(-2) + 32*mg**2*s*u1**(-2) + 16*mg**2
     +    *u1**(-1) + 16*mst1**2*s**(-1)*t1*u1**(-1) + 16*mst1**2*
     +    s**(-1)*t1**2*u1**(-2) - 32*mst1**2*s*u1**(-2) - 16*mst1**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1) + 16*
     +    m1**2*s**(-1)*t1**2*u1**(-2) + 32*m1**2*s*u1**(-2) + 32*m1**2
     +    *t1*u1**(-2) + 32*m1**2*u1**(-1) - 32*m1**4*mg**2*s*u1**(-4)
     +     + 32*m1**4*mst1**2*s*u1**(-4) - 32*m1**4*s*u1**(-3) - 32*
     +    m1**6*s*u1**(-4) + 16*mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*
     +    s**(-1)*t1**2*u1**(-2) - 16*mst1**2*s**(-1)*t1*u1**(-1) - 16*
     +    mst1**2*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1) + 16*m1**2*
     +    s**(-1)*t1**2*u1**(-2) + 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4
     +    *mst1**2*s*u1**(-4) + 64*m1**4*s*u1**(-3) - 64*m1**4*t1*
     +    u1**(-3) + 96*m1**6*s*u1**(-4) - 16*mg**2*s**(-1)*t1*u1**(-1)
     +     - 16*mg**2*s**(-1)*t1**2*u1**(-2) + 16*mst1**2*s**(-1)*t1*
     +    u1**(-1) + 16*mst1**2*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*ct*cb*m1*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 + 64*h_l*ct*cb*m1*
     +    m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*
     +    t1 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 96*h_l*ct
     +    *cb*m1*m2**2*t1*u1**(-1) - 16*h_l*ct*cb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*ct*cb*m1*s**(-1)*t1*u2 - 32*h_l*ct*cb*m1*
     +    s**(-1)*t1**2 - 80*h_l*ct*cb*m1*s*t1*u1**(-1) - 64*h_l*ct*cb*
     +    m1*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*t1 - 64*h_l*ct*cb*m1*
     +    t1**2*u1**(-1) - 32*h_l*ct*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 96*h_l*ct*cb*m1**3*s**(-2)*t1*u2 + 32*h_l*ct*cb*m1**3*
     +    s**(-2)*t1**2 + 64*h_l*ct*cb*m1**3*s**(-2)*u2**2 + 192*h_l*ct
     +    *cb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*s**(-1)
     +    *t1 + 64*h_l*ct*cb*m1**3*s**(-1)*t1**2*u1**(-1) + 128*h_l*ct*
     +    cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*h_l*ct*cb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*cb*m1**3*t1*u1**(-1) + 160*
     +    h_l*ct*cb*m1**3*u1**(-1)*u2 + 16*h_l*ct*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*ct*cb*m1**5*s**(-1)*u1**(-1)*u2 + 32*h_r*st
     +    *sb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*st*sb*m1*m2**2*s**(-2)*
     +    t1**2 + 64*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r
     +    *st*sb*m1*m2**2*s**(-1)*t1 + 64*h_r*st*sb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*st*sb*m1*m2**2*t1*u1**(-1) - 16*h_r*
     +    st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 32*h_r*st*sb*m1*s**(-1)*
     +    t1*u2 - 32*h_r*st*sb*m1*s**(-1)*t1**2 - 80*h_r*st*sb*m1*s*t1*
     +    u1**(-1) - 64*h_r*st*sb*m1*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1*
     +    t1 - 64*h_r*st*sb*m1*t1**2*u1**(-1) - 32*h_r*st*sb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*sb*m1**3*s**(-2)*t1*u2
     +     + 32*h_r*st*sb*m1**3*s**(-2)*t1**2 + 64*h_r*st*sb*m1**3*
     +    s**(-2)*u2**2 + 192*h_r*st*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_r*st*sb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 128*h_r*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*
     +    h_r*st*sb*m1**3*s**(-1)*u2 + 64*h_r*st*sb*m1**3*t1*u1**(-1)
     +     + 160*h_r*st*sb*m1**3*u1**(-1)*u2 + 16*h_r*st*sb*m1**5*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*st*sb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*ct*cb*m1*m2**2*s**(-2)*t1*u2 - 32
     +    *h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 - 64*h_l*ct*cb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*m2**2*s**(-1)*t1 - 
     +    64*h_l*ct*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 96*h_l*ct*cb*
     +    m1*m2**2*t1*u1**(-1) + 16*h_l*ct*cb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*ct*cb*m1*s**(-1)*t1*u2 + 32*h_l*ct*cb*m1*
     +    s**(-1)*t1**2 + 80*h_l*ct*cb*m1*s*t1*u1**(-1) + 64*h_l*ct*cb*
     +    m1*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*t1 + 64*h_l*ct*cb*m1*
     +    t1**2*u1**(-1) + 32*h_l*ct*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 96*h_l*ct*cb*m1**3*s**(-2)*t1*u2 - 32*h_l*ct*cb*m1**3*
     +    s**(-2)*t1**2 - 64*h_l*ct*cb*m1**3*s**(-2)*u2**2 - 192*h_l*ct
     +    *cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*s**(-1)
     +    *t1 - 64*h_l*ct*cb*m1**3*s**(-1)*t1**2*u1**(-1) - 128*h_l*ct*
     +    cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*h_l*ct*cb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*cb*m1**3*t1*u1**(-1) - 160*h_l
     +    *ct*cb*m1**3*u1**(-1)*u2 - 16*h_l*ct*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*ct*cb*m1**5*s**(-1)*u1**(-1)*u2 - 32*h_r*st
     +    *sb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*st*sb*m1*m2**2*s**(-2)*
     +    t1**2 - 64*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r
     +    *st*sb*m1*m2**2*s**(-1)*t1 - 64*h_r*st*sb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_r*st*sb*m1*m2**2*t1*u1**(-1) + 16*h_r*
     +    st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 32*h_r*st*sb*m1*s**(-1)*
     +    t1*u2 + 32*h_r*st*sb*m1*s**(-1)*t1**2 + 80*h_r*st*sb*m1*s*t1*
     +    u1**(-1) + 64*h_r*st*sb*m1*t1*u1**(-1)*u2 + 64*h_r*st*sb*m1*
     +    t1 + 64*h_r*st*sb*m1*t1**2*u1**(-1) + 32*h_r*st*sb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*sb*m1**3*s**(-2)*t1*u2
     +     - 32*h_r*st*sb*m1**3*s**(-2)*t1**2 - 64*h_r*st*sb*m1**3*
     +    s**(-2)*u2**2 - 192*h_r*st*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_r*st*sb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 128*h_r*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*
     +    h_r*st*sb*m1**3*s**(-1)*u2 - 64*h_r*st*sb*m1**3*t1*u1**(-1)
     +     - 160*h_r*st*sb*m1**3*u1**(-1)*u2 - 16*h_r*st*sb*m1**5*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*sb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*cb*m1*s**(-2)*t1 + 32*h_l*ct*cb*m1*s**(-2)
     +    *u2 + 16*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*m1*
     +    s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*m1*s**(-2)*t1 + 32*h_r*st*
     +    sb*m1*s**(-2)*u2 + 16*h_r*st*sb*m1*s**(-1)*t1*u1**(-1) + 32*
     +    h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1*s**(-2)*t1 - 32*h_l*ct*cb*m1*
     +    s**(-2)*u2 - 16*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*ct*
     +    cb*m1*s**(-1)*u1**(-1)*u2 - 32*h_r*st*sb*m1*s**(-2)*t1 - 32*
     +    h_r*st*sb*m1*s**(-2)*u2 - 16*h_r*st*sb*m1*s**(-1)*t1*u1**(-1)
     +     - 32*h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*cb*m1*s**(-1)*t1*u1 + 32*
     +    h_l*ct*cb*m1*s**(-1)*t1**2 + 32*h_l*ct*cb*m1*s**(-1)*u1**2 - 
     +    128*h_l*ct*cb*m1**3 + 64*h_r*st*sb*m1*s**(-1)*t1*u1 + 32*h_r*
     +    st*sb*m1*s**(-1)*t1**2 + 32*h_r*st*sb*m1*s**(-1)*u1**2 - 128*
     +    h_r*st*sb*m1**3 )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*ct*
     +    cb*m1*s**(-1) + 64*h_l*ct*cb*m1**3*u1**(-2) - 32*h_r*st*sb*m1
     +    *s**(-1)*t1*u1**(-1) - 32*h_r*st*sb*m1*s**(-1) + 64*h_r*st*sb
     +    *m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*ct*cb*
     +    m1**3*u1**(-2) + 32*h_r*st*sb*m1*s**(-1)*t1*u1**(-1) - 64*h_r
     +    *st*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*ct*sb*m1*m2**2*s**(-2)*t1*
     +    u2 - 32*h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 - 64*h_l*ct*sb*m1*
     +    m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*
     +    t1 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 96*h_l*ct
     +    *sb*m1*m2**2*t1*u1**(-1) + 16*h_l*ct*sb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*ct*sb*m1*s**(-1)*t1*u2 + 32*h_l*ct*sb*m1*
     +    s**(-1)*t1**2 + 80*h_l*ct*sb*m1*s*t1*u1**(-1) + 64*h_l*ct*sb*
     +    m1*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*t1 + 64*h_l*ct*sb*m1*
     +    t1**2*u1**(-1) + 32*h_l*ct*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 96*h_l*ct*sb*m1**3*s**(-2)*t1*u2 - 32*h_l*ct*sb*m1**3*
     +    s**(-2)*t1**2 - 64*h_l*ct*sb*m1**3*s**(-2)*u2**2 - 192*h_l*ct
     +    *sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*s**(-1)
     +    *t1 - 64*h_l*ct*sb*m1**3*s**(-1)*t1**2*u1**(-1) - 128*h_l*ct*
     +    sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*h_l*ct*sb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*sb*m1**3*t1*u1**(-1) - 
     +    160*h_l*ct*sb*m1**3*u1**(-1)*u2 - 16*h_l*ct*sb*m1**5*s**(-1)*
     +    t1*u1**(-1) - 32*h_l*ct*sb*m1**5*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *st*cb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*st*cb*m1*m2**2*s**(-2)
     +    *t1**2 + 64*h_r*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*st*cb*m1*m2**2*s**(-1)*t1 + 64*h_r*st*cb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) + 96*h_r*st*cb*m1*m2**2*t1*u1**(-1) - 16*h_r*
     +    st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 32*h_r*st*cb*m1*s**(-1)*
     +    t1*u2 - 32*h_r*st*cb*m1*s**(-1)*t1**2 - 80*h_r*st*cb*m1*s*t1*
     +    u1**(-1) - 64*h_r*st*cb*m1*t1*u1**(-1)*u2 - 64*h_r*st*cb*m1*
     +    t1 - 64*h_r*st*cb*m1*t1**2*u1**(-1) - 32*h_r*st*cb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*cb*m1**3*s**(-2)*t1*u2
     +     + 32*h_r*st*cb*m1**3*s**(-2)*t1**2 + 64*h_r*st*cb*m1**3*
     +    s**(-2)*u2**2 + 192*h_r*st*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_r*st*cb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 128*h_r*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*
     +    h_r*st*cb*m1**3*s**(-1)*u2 + 64*h_r*st*cb*m1**3*t1*u1**(-1)
     +     + 160*h_r*st*cb*m1**3*u1**(-1)*u2 + 16*h_r*st*cb*m1**5*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*st*cb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*ct*sb*m1*m2**2*s**(-2)*t1*u2 + 32*
     +    h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 + 64*h_l*ct*sb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*m2**2*s**(-1)*t1 + 
     +    64*h_l*ct*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 96*h_l*ct*sb*
     +    m1*m2**2*t1*u1**(-1) - 16*h_l*ct*sb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*ct*sb*m1*s**(-1)*t1*u2 - 32*h_l*ct*sb*m1*
     +    s**(-1)*t1**2 - 80*h_l*ct*sb*m1*s*t1*u1**(-1) - 64*h_l*ct*sb*
     +    m1*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*t1 - 64*h_l*ct*sb*m1*
     +    t1**2*u1**(-1) - 32*h_l*ct*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 96*h_l*ct*sb*m1**3*s**(-2)*t1*u2 + 32*h_l*ct*sb*m1**3*
     +    s**(-2)*t1**2 + 64*h_l*ct*sb*m1**3*s**(-2)*u2**2 + 192*h_l*ct
     +    *sb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*s**(-1)
     +    *t1 + 64*h_l*ct*sb*m1**3*s**(-1)*t1**2*u1**(-1) + 128*h_l*ct*
     +    sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*h_l*ct*sb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*sb*m1**3*t1*u1**(-1) + 160*h_l*ct
     +    *sb*m1**3*u1**(-1)*u2 + 16*h_l*ct*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*ct*sb*m1**5*s**(-1)*u1**(-1)*u2 - 32*h_r*st
     +    *cb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*st*cb*m1*m2**2*s**(-2)*
     +    t1**2 - 64*h_r*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r
     +    *st*cb*m1*m2**2*s**(-1)*t1 - 64*h_r*st*cb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_r*st*cb*m1*m2**2*t1*u1**(-1) + 16*h_r*
     +    st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 32*h_r*st*cb*m1*s**(-1)*
     +    t1*u2 + 32*h_r*st*cb*m1*s**(-1)*t1**2 + 80*h_r*st*cb*m1*s*t1*
     +    u1**(-1) + 64*h_r*st*cb*m1*t1*u1**(-1)*u2 + 64*h_r*st*cb*m1*
     +    t1 + 64*h_r*st*cb*m1*t1**2*u1**(-1) + 32*h_r*st*cb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*cb*m1**3*s**(-2)*t1*u2
     +     - 32*h_r*st*cb*m1**3*s**(-2)*t1**2 - 64*h_r*st*cb*m1**3*
     +    s**(-2)*u2**2 - 192*h_r*st*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_r*st*cb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 128*h_r*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*
     +    h_r*st*cb*m1**3*s**(-1)*u2 - 64*h_r*st*cb*m1**3*t1*u1**(-1)
     +     - 160*h_r*st*cb*m1**3*u1**(-1)*u2 - 16*h_r*st*cb*m1**5*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*cb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*sb*m1*s**(-2)*t1 - 32*h_l*ct*sb*m1*
     +    s**(-2)*u2 - 16*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*ct*
     +    sb*m1*s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*m1*s**(-2)*t1 + 32*
     +    h_r*st*cb*m1*s**(-2)*u2 + 16*h_r*st*cb*m1*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1*s**(-2)*t1 + 32*h_l*ct*sb*m1*s**(-2)*
     +    u2 + 16*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1*
     +    s**(-1)*u1**(-1)*u2 - 32*h_r*st*cb*m1*s**(-2)*t1 - 32*h_r*st*
     +    cb*m1*s**(-2)*u2 - 16*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*sb*m1*s**(-1)*t1*u1 - 32
     +    *h_l*ct*sb*m1*s**(-1)*t1**2 - 32*h_l*ct*sb*m1*s**(-1)*u1**2
     +     + 128*h_l*ct*sb*m1**3 + 64*h_r*st*cb*m1*s**(-1)*t1*u1 + 32*
     +    h_r*st*cb*m1*s**(-1)*t1**2 + 32*h_r*st*cb*m1*s**(-1)*u1**2 - 
     +    128*h_r*st*cb*m1**3 )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*
     +    m1*s**(-1) - 64*h_l*ct*sb*m1**3*u1**(-2) - 32*h_r*st*cb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*cb*m1*s**(-1) + 64*h_r*st*cb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*ct*sb
     +    *m1**3*u1**(-2) + 32*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) - 64*
     +    h_r*st*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*susy_t*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2t*m1**(-2)*mg**2*s**(-1)*t1 + 16*c2t*
     +    m1**(-2)*mg**2*s**(-1)*t1**2*u1**(-1) + 8*c2t*m1**(-2)*mg**2*
     +    s**(-1)*u1 + 8*c2t*m1**(-2)*mg**2*s*u1**(-1) + 16*c2t*
     +    m1**(-2)*mg**2*t1*u1**(-1) + 16*c2t*m1**(-2)*mg**2 - 16*c2t*
     +    m1**(-2)*mst2**2*s**(-1)*t1 - 16*c2t*m1**(-2)*mst2**2*s**(-1)
     +    *t1**2*u1**(-1) - 8*c2t*m1**(-2)*mst2**2*s**(-1)*u1 - 8*c2t*
     +    m1**(-2)*mst2**2*s*u1**(-1) - 16*c2t*m1**(-2)*mst2**2*t1*
     +    u1**(-1) - 16*c2t*m1**(-2)*mst2**2 + 16*c2t*m1**2*s*u1**(-2)
     +     + 16*c2t*mg**2*s*u1**(-2) - 16*c2t*mst2**2*s*u1**(-2) + 16*
     +    c2t*s**(-1)*t1 + 16*c2t*s**(-1)*t1**2*u1**(-1) + 8*c2t*
     +    s**(-1)*u1 + 8*c2t*s*u1**(-1) + 16*c2t*t1*u1**(-1) + 16*c2t )
      MM_s = MM_s + SCB(2,3)*susy_t*Cf*Pi*alphas*prefac * ( 2*m1**(-2)*
     +    mg**2*born - 2*m1**(-2)*mst2**2*born - 2*born )
      MM_s = MM_s + SCB(2,3)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 64*m1**2*mg**2*s*u1**(-3) + 64*m1**2*mg**2*
     +    t1*u1**(-3) + 32*m1**2*mg**2*u1**(-2) - 64*m1**2*mst2**2*s*
     +    u1**(-3) - 64*m1**2*mst2**2*t1*u1**(-3) - 32*m1**2*mst2**2*
     +    u1**(-2) - 32*m1**2*s**(-1)*t1*u1**(-1) - 32*m1**2*s**(-1)*
     +    t1**2*u1**(-2) - 16*m1**2*t1*u1**(-2) - 16*m1**2*u1**(-1) + 
     +    64*m1**4*s*u1**(-3) + 64*m1**4*t1*u1**(-3) + 32*m1**4*
     +    u1**(-2) - 32*mg**2*s**(-1)*t1*u1**(-1) - 32*mg**2*s**(-1)*
     +    t1**2*u1**(-2) - 16*mg**2*t1*u1**(-2) - 16*mg**2*u1**(-1) + 
     +    32*mst2**2*s**(-1)*t1*u1**(-1) + 32*mst2**2*s**(-1)*t1**2*
     +    u1**(-2) + 16*mst2**2*t1*u1**(-2) + 16*mst2**2*u1**(-1) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * ( 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) + 32*s2t*
     +    m1*mg*s**(-1)*t1**2*u1**(-2) - 16*s2t*m1*mg*s*u1**(-2) + 16*
     +    s2t*m1*mg*t1*u1**(-2) - 64*s2t*m1**3*mg*s*u1**(-3) - 64*s2t*
     +    m1**3*mg*t1*u1**(-3) - 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) - 32*s2t*m1*mg*
     +    s**(-1)*t1**2*u1**(-2) + 16*s2t*m1*mg*s*u1**(-2) - 16*s2t*m1*
     +    mg*t1*u1**(-2) + 64*s2t*m1**3*mg*s*u1**(-3) + 64*s2t*m1**3*mg
     +    *t1*u1**(-3) + 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1) + 16*
     +    m1**2*s**(-1)*t1**2*u1**(-2) + 32*m1**2*s*u1**(-2) + 32*m1**2
     +    *t1*u1**(-2) + 32*m1**2*u1**(-1) - 32*m1**4*mg**2*s*u1**(-4)
     +     + 32*m1**4*mst2**2*s*u1**(-4) - 32*m1**4*s*u1**(-3) - 32*
     +    m1**6*s*u1**(-4) + 16*mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*
     +    s**(-1)*t1**2*u1**(-2) - 16*mst2**2*s**(-1)*t1*u1**(-1) - 16*
     +    mst2**2*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*s**(-1)*t1*u1**(-1) + 16*m1**2*
     +    s**(-1)*t1**2*u1**(-2) + 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4
     +    *mst2**2*s*u1**(-4) + 64*m1**4*s*u1**(-3) - 64*m1**4*t1*
     +    u1**(-3) + 96*m1**6*s*u1**(-4) - 16*mg**2*s**(-1)*t1*u1**(-1)
     +     - 16*mg**2*s**(-1)*t1**2*u1**(-2) + 16*mst2**2*s**(-1)*t1*
     +    u1**(-1) + 16*mst2**2*s**(-1)*t1**2*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 48*m1**2*mg**2*s*u1**(-3) - 64*
     +    m1**2*mg**2*t1*u1**(-3) - 32*m1**2*mg**2*u1**(-2) + 48*m1**2*
     +    mst2**2*s*u1**(-3) + 64*m1**2*mst2**2*t1*u1**(-3) + 32*m1**2*
     +    mst2**2*u1**(-2) + 16*m1**2*s**(-1)*t1*u1**(-1) + 16*m1**2*
     +    s**(-1)*t1**2*u1**(-2) - 32*m1**2*s*u1**(-2) - 32*m1**2*t1*
     +    u1**(-2) - 16*m1**2*u1**(-1) + 32*m1**4*mg**2*s*u1**(-4) - 32
     +    *m1**4*mst2**2*s*u1**(-4) - 16*m1**4*s*u1**(-3) - 64*m1**4*t1
     +    *u1**(-3) - 32*m1**4*u1**(-2) + 32*m1**6*s*u1**(-4) + 16*
     +    mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*s**(-1)*t1**2*u1**(-2)
     +     - 32*mg**2*s*u1**(-2) - 16*mg**2*u1**(-1) - 16*mst2**2*
     +    s**(-1)*t1*u1**(-1) - 16*mst2**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    mst2**2*s*u1**(-2) + 16*mst2**2*u1**(-1) )
      MM_s = MM_s + SCB(2,3)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 48*m1**2*mg**2*s*u1**(-3) + 64*m1**2*mg**2*
     +    t1*u1**(-3) + 32*m1**2*mg**2*u1**(-2) - 48*m1**2*mst2**2*s*
     +    u1**(-3) - 64*m1**2*mst2**2*t1*u1**(-3) - 32*m1**2*mst2**2*
     +    u1**(-2) - 48*m1**2*s**(-1)*t1*u1**(-1) - 48*m1**2*s**(-1)*
     +    t1**2*u1**(-2) - 32*m1**2*s*u1**(-2) - 16*m1**2*u1**(-1) - 32
     +    *m1**4*mg**2*s*u1**(-4) + 32*m1**4*mst2**2*s*u1**(-4) - 48*
     +    m1**4*s*u1**(-3) + 128*m1**4*t1*u1**(-3) + 32*m1**4*u1**(-2)
     +     - 96*m1**6*s*u1**(-4) - 16*mg**2*s**(-1)*t1*u1**(-1) - 16*
     +    mg**2*s**(-1)*t1**2*u1**(-2) + 32*mg**2*s*u1**(-2) + 16*mg**2
     +    *u1**(-1) + 16*mst2**2*s**(-1)*t1*u1**(-1) + 16*mst2**2*
     +    s**(-1)*t1**2*u1**(-2) - 32*mst2**2*s*u1**(-2) - 16*mst2**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*st*cb*m1*m2**2*s**(-2)*t1*
     +    u2 - 32*h_l*st*cb*m1*m2**2*s**(-2)*t1**2 - 64*h_l*st*cb*m1*
     +    m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*m2**2*s**(-1)*
     +    t1 - 64*h_l*st*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 96*h_l*st
     +    *cb*m1*m2**2*t1*u1**(-1) + 16*h_l*st*cb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*st*cb*m1*s**(-1)*t1*u2 + 32*h_l*st*cb*m1*
     +    s**(-1)*t1**2 + 80*h_l*st*cb*m1*s*t1*u1**(-1) + 64*h_l*st*cb*
     +    m1*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*t1 + 64*h_l*st*cb*m1*
     +    t1**2*u1**(-1) + 32*h_l*st*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 96*h_l*st*cb*m1**3*s**(-2)*t1*u2 - 32*h_l*st*cb*m1**3*
     +    s**(-2)*t1**2 - 64*h_l*st*cb*m1**3*s**(-2)*u2**2 - 192*h_l*st
     +    *cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*s**(-1)
     +    *t1 - 64*h_l*st*cb*m1**3*s**(-1)*t1**2*u1**(-1) - 128*h_l*st*
     +    cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*h_l*st*cb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*cb*m1**3*t1*u1**(-1) - 
     +    160*h_l*st*cb*m1**3*u1**(-1)*u2 - 16*h_l*st*cb*m1**5*s**(-1)*
     +    t1*u1**(-1) - 32*h_l*st*cb*m1**5*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *ct*sb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*ct*sb*m1*m2**2*s**(-2)
     +    *t1**2 + 64*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*sb*m1*m2**2*s**(-1)*t1 + 64*h_r*ct*sb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) + 96*h_r*ct*sb*m1*m2**2*t1*u1**(-1) - 16*h_r*
     +    ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb*m1*s**(-1)*
     +    t1*u2 - 32*h_r*ct*sb*m1*s**(-1)*t1**2 - 80*h_r*ct*sb*m1*s*t1*
     +    u1**(-1) - 64*h_r*ct*sb*m1*t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1*
     +    t1 - 64*h_r*ct*sb*m1*t1**2*u1**(-1) - 32*h_r*ct*sb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*sb*m1**3*s**(-2)*t1*u2
     +     + 32*h_r*ct*sb*m1**3*s**(-2)*t1**2 + 64*h_r*ct*sb*m1**3*
     +    s**(-2)*u2**2 + 192*h_r*ct*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_r*ct*sb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 128*h_r*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*
     +    h_r*ct*sb*m1**3*s**(-1)*u2 + 64*h_r*ct*sb*m1**3*t1*u1**(-1)
     +     + 160*h_r*ct*sb*m1**3*u1**(-1)*u2 + 16*h_r*ct*sb*m1**5*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*ct*sb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*st*cb*m1*m2**2*s**(-2)*t1*u2 + 32*
     +    h_l*st*cb*m1*m2**2*s**(-2)*t1**2 + 64*h_l*st*cb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*m2**2*s**(-1)*t1 + 
     +    64*h_l*st*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 96*h_l*st*cb*
     +    m1*m2**2*t1*u1**(-1) - 16*h_l*st*cb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*st*cb*m1*s**(-1)*t1*u2 - 32*h_l*st*cb*m1*
     +    s**(-1)*t1**2 - 80*h_l*st*cb*m1*s*t1*u1**(-1) - 64*h_l*st*cb*
     +    m1*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*t1 - 64*h_l*st*cb*m1*
     +    t1**2*u1**(-1) - 32*h_l*st*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 96*h_l*st*cb*m1**3*s**(-2)*t1*u2 + 32*h_l*st*cb*m1**3*
     +    s**(-2)*t1**2 + 64*h_l*st*cb*m1**3*s**(-2)*u2**2 + 192*h_l*st
     +    *cb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1**3*s**(-1)
     +    *t1 + 64*h_l*st*cb*m1**3*s**(-1)*t1**2*u1**(-1) + 128*h_l*st*
     +    cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*h_l*st*cb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*cb*m1**3*t1*u1**(-1) + 160*h_l*st
     +    *cb*m1**3*u1**(-1)*u2 + 16*h_l*st*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*st*cb*m1**5*s**(-1)*u1**(-1)*u2 - 32*h_r*ct
     +    *sb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*ct*sb*m1*m2**2*s**(-2)*
     +    t1**2 - 64*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r
     +    *ct*sb*m1*m2**2*s**(-1)*t1 - 64*h_r*ct*sb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_r*ct*sb*m1*m2**2*t1*u1**(-1) + 16*h_r*
     +    ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 32*h_r*ct*sb*m1*s**(-1)*
     +    t1*u2 + 32*h_r*ct*sb*m1*s**(-1)*t1**2 + 80*h_r*ct*sb*m1*s*t1*
     +    u1**(-1) + 64*h_r*ct*sb*m1*t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1*
     +    t1 + 64*h_r*ct*sb*m1*t1**2*u1**(-1) + 32*h_r*ct*sb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*sb*m1**3*s**(-2)*t1*u2
     +     - 32*h_r*ct*sb*m1**3*s**(-2)*t1**2 - 64*h_r*ct*sb*m1**3*
     +    s**(-2)*u2**2 - 192*h_r*ct*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_r*ct*sb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 128*h_r*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*
     +    h_r*ct*sb*m1**3*s**(-1)*u2 - 64*h_r*ct*sb*m1**3*t1*u1**(-1)
     +     - 160*h_r*ct*sb*m1**3*u1**(-1)*u2 - 16*h_r*ct*sb*m1**5*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*s**(-2)*t1 - 32*h_l*st*cb*m1*
     +    s**(-2)*u2 - 16*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*st*
     +    cb*m1*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*m1*s**(-2)*t1 + 32*
     +    h_r*ct*sb*m1*s**(-2)*u2 + 16*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*s**(-2)*t1 + 32*h_l*st*cb*m1*s**(-2)*
     +    u2 + 16*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1*
     +    s**(-1)*u1**(-1)*u2 - 32*h_r*ct*sb*m1*s**(-2)*t1 - 32*h_r*ct*
     +    sb*m1*s**(-2)*u2 - 16*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*cb*m1*s**(-1)*t1*u1 - 32
     +    *h_l*st*cb*m1*s**(-1)*t1**2 - 32*h_l*st*cb*m1*s**(-1)*u1**2
     +     + 128*h_l*st*cb*m1**3 + 64*h_r*ct*sb*m1*s**(-1)*t1*u1 + 32*
     +    h_r*ct*sb*m1*s**(-1)*t1**2 + 32*h_r*ct*sb*m1*s**(-1)*u1**2 - 
     +    128*h_r*ct*sb*m1**3 )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*
     +    m1*s**(-1) - 64*h_l*st*cb*m1**3*u1**(-2) - 32*h_r*ct*sb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb*m1*s**(-1) + 64*h_r*ct*sb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*st*cb
     +    *m1**3*u1**(-2) + 32*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) - 64*
     +    h_r*ct*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*st*sb*m1*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*st*sb*m1*m2**2*s**(-2)*t1**2 + 64*h_l*st*sb*m1*
     +    m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*m2**2*s**(-1)*
     +    t1 + 64*h_l*st*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 96*h_l*st
     +    *sb*m1*m2**2*t1*u1**(-1) - 16*h_l*st*sb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*st*sb*m1*s**(-1)*t1*u2 - 32*h_l*st*sb*m1*
     +    s**(-1)*t1**2 - 80*h_l*st*sb*m1*s*t1*u1**(-1) - 64*h_l*st*sb*
     +    m1*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*t1 - 64*h_l*st*sb*m1*
     +    t1**2*u1**(-1) - 32*h_l*st*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 96*h_l*st*sb*m1**3*s**(-2)*t1*u2 + 32*h_l*st*sb*m1**3*
     +    s**(-2)*t1**2 + 64*h_l*st*sb*m1**3*s**(-2)*u2**2 + 192*h_l*st
     +    *sb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1**3*s**(-1)
     +    *t1 + 64*h_l*st*sb*m1**3*s**(-1)*t1**2*u1**(-1) + 128*h_l*st*
     +    sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*h_l*st*sb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*sb*m1**3*t1*u1**(-1) + 160*
     +    h_l*st*sb*m1**3*u1**(-1)*u2 + 16*h_l*st*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*st*sb*m1**5*s**(-1)*u1**(-1)*u2 + 32*h_r*ct
     +    *cb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*ct*cb*m1*m2**2*s**(-2)*
     +    t1**2 + 64*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r
     +    *ct*cb*m1*m2**2*s**(-1)*t1 + 64*h_r*ct*cb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*ct*cb*m1*m2**2*t1*u1**(-1) - 16*h_r*
     +    ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb*m1*s**(-1)*
     +    t1*u2 - 32*h_r*ct*cb*m1*s**(-1)*t1**2 - 80*h_r*ct*cb*m1*s*t1*
     +    u1**(-1) - 64*h_r*ct*cb*m1*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1*
     +    t1 - 64*h_r*ct*cb*m1*t1**2*u1**(-1) - 32*h_r*ct*cb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*cb*m1**3*s**(-2)*t1*u2
     +     + 32*h_r*ct*cb*m1**3*s**(-2)*t1**2 + 64*h_r*ct*cb*m1**3*
     +    s**(-2)*u2**2 + 192*h_r*ct*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_r*ct*cb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 128*h_r*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 128*
     +    h_r*ct*cb*m1**3*s**(-1)*u2 + 64*h_r*ct*cb*m1**3*t1*u1**(-1)
     +     + 160*h_r*ct*cb*m1**3*u1**(-1)*u2 + 16*h_r*ct*cb*m1**5*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*ct*cb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*st*sb*m1*m2**2*s**(-2)*t1*u2 - 32
     +    *h_l*st*sb*m1*m2**2*s**(-2)*t1**2 - 64*h_l*st*sb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*m2**2*s**(-1)*t1 - 
     +    64*h_l*st*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 96*h_l*st*sb*
     +    m1*m2**2*t1*u1**(-1) + 16*h_l*st*sb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*st*sb*m1*s**(-1)*t1*u2 + 32*h_l*st*sb*m1*
     +    s**(-1)*t1**2 + 80*h_l*st*sb*m1*s*t1*u1**(-1) + 64*h_l*st*sb*
     +    m1*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*t1 + 64*h_l*st*sb*m1*
     +    t1**2*u1**(-1) + 32*h_l*st*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 96*h_l*st*sb*m1**3*s**(-2)*t1*u2 - 32*h_l*st*sb*m1**3*
     +    s**(-2)*t1**2 - 64*h_l*st*sb*m1**3*s**(-2)*u2**2 - 192*h_l*st
     +    *sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1**3*s**(-1)
     +    *t1 - 64*h_l*st*sb*m1**3*s**(-1)*t1**2*u1**(-1) - 128*h_l*st*
     +    sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*h_l*st*sb*m1**3*s**(-1)
     +    *u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*sb*m1**3*t1*u1**(-1) - 160*h_l
     +    *st*sb*m1**3*u1**(-1)*u2 - 16*h_l*st*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*st*sb*m1**5*s**(-1)*u1**(-1)*u2 - 32*h_r*ct
     +    *cb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*ct*cb*m1*m2**2*s**(-2)*
     +    t1**2 - 64*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r
     +    *ct*cb*m1*m2**2*s**(-1)*t1 - 64*h_r*ct*cb*m1*m2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_r*ct*cb*m1*m2**2*t1*u1**(-1) + 16*h_r*
     +    ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 32*h_r*ct*cb*m1*s**(-1)*
     +    t1*u2 + 32*h_r*ct*cb*m1*s**(-1)*t1**2 + 80*h_r*ct*cb*m1*s*t1*
     +    u1**(-1) + 64*h_r*ct*cb*m1*t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1*
     +    t1 + 64*h_r*ct*cb*m1*t1**2*u1**(-1) + 32*h_r*ct*cb*m1**3*
     +    m2**2*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*cb*m1**3*s**(-2)*t1*u2
     +     - 32*h_r*ct*cb*m1**3*s**(-2)*t1**2 - 64*h_r*ct*cb*m1**3*
     +    s**(-2)*u2**2 - 192*h_r*ct*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_r*ct*cb*m1**3*s**(-1)*t1 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 128*h_r*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 128*
     +    h_r*ct*cb*m1**3*s**(-1)*u2 - 64*h_r*ct*cb*m1**3*t1*u1**(-1)
     +     - 160*h_r*ct*cb*m1**3*u1**(-1)*u2 - 16*h_r*ct*cb*m1**5*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb*m1**5*s**(-1)*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*s**(-2)*t1 + 32*h_l*st*sb*m1*s**(-2)
     +    *u2 + 16*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1*
     +    s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m1*s**(-2)*t1 + 32*h_r*ct*
     +    cb*m1*s**(-2)*u2 + 16*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1) + 32*
     +    h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*s**(-2)*t1 - 32*h_l*st*sb*m1*
     +    s**(-2)*u2 - 16*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*st*
     +    sb*m1*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*cb*m1*s**(-2)*t1 - 32*
     +    h_r*ct*cb*m1*s**(-2)*u2 - 16*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1)
     +     - 32*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*sb*m1*s**(-1)*t1*u1 + 32*
     +    h_l*st*sb*m1*s**(-1)*t1**2 + 32*h_l*st*sb*m1*s**(-1)*u1**2 - 
     +    128*h_l*st*sb*m1**3 + 64*h_r*ct*cb*m1*s**(-1)*t1*u1 + 32*h_r*
     +    ct*cb*m1*s**(-1)*t1**2 + 32*h_r*ct*cb*m1*s**(-1)*u1**2 - 128*
     +    h_r*ct*cb*m1**3 )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*st*
     +    sb*m1*s**(-1) + 64*h_l*st*sb*m1**3*u1**(-2) - 32*h_r*ct*cb*m1
     +    *s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb*m1*s**(-1) + 64*h_r*ct*cb
     +    *m1**3*u1**(-2) )
      MM_s = MM_s + SCB(2,3)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*st*sb*
     +    m1**3*u1**(-2) + 32*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1) - 64*h_r
     +    *ct*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 96*h_l*ct*cb*m1*m2**2*s**(-2)*t1*
     +    u2 - 64*h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*ct*cb*m1*
     +    m2**2*s**(-2)*u2**2 - 192*h_l*ct*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 128*h_l*ct*cb*m1*m2**2*s**(-1)*t1 - 128*h_l*ct*
     +    cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2 - 
     +    160*h_l*ct*cb*m1*m2**2*t1*u1**(-1) - 96*h_l*ct*cb*m1*m2**2*
     +    u1**(-1)*u2 + 32*h_l*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*m1*
     +    s**(-1)*t1*u2 + 32*h_l*ct*cb*m1*s**(-1)*u2**2 + 80*h_l*ct*cb*
     +    m1*s*u1**(-1)*u2 + 64*h_l*ct*cb*m1*t1*u1**(-1)*u2 + 64*h_l*ct
     +    *cb*m1*u1**(-1)*u2**2 + 64*h_l*ct*cb*m1*u2 - 32*h_l*ct*cb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*ct*cb*m1**3*s**(-2)*
     +    t1*u2 - 32*h_l*ct*cb*m1**3*s**(-2)*u2**2 - 64*h_l*ct*cb*m1**3
     +    *s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*ct*cb*m1**3*s**(-1)*u2 - 64*h_l*ct*cb
     +    *m1**3*u1**(-1)*u2 - 16*h_l*ct*cb*m1**5*s**(-1)*u1**(-1)*u2
     +     - 96*h_r*st*sb*m1*m2**2*s**(-2)*t1*u2 - 64*h_r*st*sb*m1*
     +    m2**2*s**(-2)*t1**2 - 32*h_r*st*sb*m1*m2**2*s**(-2)*u2**2 - 
     +    192*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_r*st*sb
     +    *m1*m2**2*s**(-1)*t1 - 128*h_r*st*sb*m1*m2**2*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_r*st*sb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*
     +    h_r*st*sb*m1*m2**2*s**(-1)*u2 - 160*h_r*st*sb*m1*m2**2*t1*
     +    u1**(-1) - 96*h_r*st*sb*m1*m2**2*u1**(-1)*u2 + 32*h_r*st*sb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2 + 32*h_r*st*sb*m1*s**(-1)*t1*u2 + 32*h_r*st*sb*m1
     +    *s**(-1)*u2**2 + 80*h_r*st*sb*m1*s*u1**(-1)*u2 + 64*h_r*st*sb
     +    *m1*t1*u1**(-1)*u2 + 64*h_r*st*sb*m1*u1**(-1)*u2**2 + 64*h_r*
     +    st*sb*m1*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_r*st*sb*m1**3*m2**2*s**(-1)*
     +    t1*u1**(-1) - 32*h_r*st*sb*m1**3*s**(-2)*t1*u2 - 32*h_r*st*sb
     +    *m1**3*s**(-2)*u2**2 - 64*h_r*st*sb*m1**3*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_r*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*
     +    sb*m1**3*s**(-1)*u2 - 64*h_r*st*sb*m1**3*u1**(-1)*u2 - 16*h_r
     +    *st*sb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*h_l*ct*cb*m1*m2**2*s**(-2)*t1*u2 + 64*
     +    h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*ct*cb*m1*m2**2*
     +    s**(-2)*u2**2 + 192*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     + 128*h_l*ct*cb*m1*m2**2*s**(-1)*t1 + 128*h_l*ct*cb*m1*m2**2
     +    *s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2 + 160*h_l*
     +    ct*cb*m1*m2**2*t1*u1**(-1) + 96*h_l*ct*cb*m1*m2**2*u1**(-1)*
     +    u2 - 32*h_l*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb
     +    *m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*m1*s**(-1)*t1*u2
     +     - 32*h_l*ct*cb*m1*s**(-1)*u2**2 - 80*h_l*ct*cb*m1*s*u1**(-1)
     +    *u2 - 64*h_l*ct*cb*m1*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*
     +    u1**(-1)*u2**2 - 64*h_l*ct*cb*m1*u2 + 32*h_l*ct*cb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*m1**3*s**(-2)*t1*u2
     +     + 32*h_l*ct*cb*m1**3*s**(-2)*u2**2 + 64*h_l*ct*cb*m1**3*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*ct*cb*m1**3*s**(-1)*u2 + 64*h_l*ct*cb*m1**3*
     +    u1**(-1)*u2 + 16*h_l*ct*cb*m1**5*s**(-1)*u1**(-1)*u2 + 96*h_r
     +    *st*sb*m1*m2**2*s**(-2)*t1*u2 + 64*h_r*st*sb*m1*m2**2*s**(-2)
     +    *t1**2 + 32*h_r*st*sb*m1*m2**2*s**(-2)*u2**2 + 192*h_r*st*sb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_r*st*sb*m1*m2**2*
     +    s**(-1)*t1 + 128*h_r*st*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*st*sb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb*
     +    m1*m2**2*s**(-1)*u2 + 160*h_r*st*sb*m1*m2**2*t1*u1**(-1) + 96
     +    *h_r*st*sb*m1*m2**2*u1**(-1)*u2 - 32*h_r*st*sb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 - 32*h_r*st*sb*m1*s**(-1)*t1*u2 - 32*h_r*st*sb*m1*s**(-1)*
     +    u2**2 - 80*h_r*st*sb*m1*s*u1**(-1)*u2 - 64*h_r*st*sb*m1*t1*
     +    u1**(-1)*u2 - 64*h_r*st*sb*m1*u1**(-1)*u2**2 - 64*h_r*st*sb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*st*sb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*st*sb*m1**3*s**(-2)*t1*u2 + 32*h_r*st*sb*
     +    m1**3*s**(-2)*u2**2 + 64*h_r*st*sb*m1**3*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_r*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb
     +    *m1**3*s**(-1)*u2 + 64*h_r*st*sb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    st*sb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_l*ct*cb*m1*s**(-2)*t1 - 32*h_l*ct*cb*m1*
     +    s**(-2)*u2 - 96*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 112*h_l*ct
     +    *cb*m1*s**(-1)*u1**(-1)*u2 - 96*h_l*ct*cb*m1*s*u1**(-2) - 96*
     +    h_l*ct*cb*m1*t1*u1**(-2) - 96*h_l*ct*cb*m1*u1**(-2)*u2 - 96*
     +    h_l*ct*cb*m1*u1**(-1) - 32*h_l*ct*cb*m1**3*s**(-1)*u1**(-1)
     +     - 64*h_l*ct*cb*m1**3*u1**(-2) - 32*h_r*st*sb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2**(-1) - 32*h_r*st*sb*m1*s**(-2)*t1 - 
     +    32*h_r*st*sb*m1*s**(-2)*u2 - 96*h_r*st*sb*m1*s**(-1)*t1*
     +    u1**(-1) - 112*h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*st*
     +    sb*m1*s*u1**(-2) - 96*h_r*st*sb*m1*t1*u1**(-2) - 96*h_r*st*sb
     +    *m1*u1**(-2)*u2 - 96*h_r*st*sb*m1*u1**(-1) - 32*h_r*st*sb*
     +    m1**3*s**(-1)*u1**(-1) - 64*h_r*st*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 
     +    32*h_l*ct*cb*m1*s**(-2)*t1 + 32*h_l*ct*cb*m1*s**(-2)*u2 + 96*
     +    h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) + 112*h_l*ct*cb*m1*s**(-1)*
     +    u1**(-1)*u2 + 96*h_l*ct*cb*m1*s*u1**(-2) + 96*h_l*ct*cb*m1*t1
     +    *u1**(-2) + 96*h_l*ct*cb*m1*u1**(-2)*u2 + 96*h_l*ct*cb*m1*
     +    u1**(-1) + 32*h_l*ct*cb*m1**3*s**(-1)*u1**(-1) + 64*h_l*ct*cb
     +    *m1**3*u1**(-2) + 32*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_r*st*sb*m1*s**(-2)*t1 + 32*h_r*st*sb*m1*
     +    s**(-2)*u2 + 96*h_r*st*sb*m1*s**(-1)*t1*u1**(-1) + 112*h_r*st
     +    *sb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*st*sb*m1*s*u1**(-2) + 96*
     +    h_r*st*sb*m1*t1*u1**(-2) + 96*h_r*st*sb*m1*u1**(-2)*u2 + 96*
     +    h_r*st*sb*m1*u1**(-1) + 32*h_r*st*sb*m1**3*s**(-1)*u1**(-1)
     +     + 64*h_r*st*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1*s**(-1)*t1*u1 - 32
     +    *h_l*ct*cb*m1*s**(-1)*t1**2 - 32*h_l*ct*cb*m1*s**(-1)*u1**2
     +     + 128*h_l*ct*cb*m1**3 - 64*h_r*st*sb*m1*s**(-1)*t1*u1 - 32*
     +    h_r*st*sb*m1*s**(-1)*t1**2 - 32*h_r*st*sb*m1*s**(-1)*u1**2 + 
     +    128*h_r*st*sb*m1**3 )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*cb*m1*u1**(-1) - 32*h_r*st*sb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(3,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*cb*m1*s**(-1) - 32*h_l*ct*cb*m1*s*t2**(-2)
     +     - 32*h_l*ct*cb*m1*s**2*u1**(-1)*t2**(-2) + 64*h_l*ct*cb*m1*
     +    u1**(-1) + 32*h_l*ct*cb*m1*t2**(-1) + 32*h_r*st*sb*m1*s**(-1)
     +     - 32*h_r*st*sb*m1*s*t2**(-2) - 32*h_r*st*sb*m1*s**2*u1**(-1)
     +    *t2**(-2) + 64*h_r*st*sb*m1*u1**(-1) + 32*h_r*st*sb*m1*
     +    t2**(-1) )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_l*ct*sb*m1*m2**2*s**(-2)*t1*u2
     +     + 64*h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*ct*sb*m1*
     +    m2**2*s**(-2)*u2**2 + 192*h_l*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 128*h_l*ct*sb*m1*m2**2*s**(-1)*t1 + 128*h_l*ct*
     +    sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2 + 
     +    160*h_l*ct*sb*m1*m2**2*t1*u1**(-1) + 96*h_l*ct*sb*m1*m2**2*
     +    u1**(-1)*u2 - 32*h_l*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*m1*
     +    s**(-1)*t1*u2 - 32*h_l*ct*sb*m1*s**(-1)*u2**2 - 80*h_l*ct*sb*
     +    m1*s*u1**(-1)*u2 - 64*h_l*ct*sb*m1*t1*u1**(-1)*u2 - 64*h_l*ct
     +    *sb*m1*u1**(-1)*u2**2 - 64*h_l*ct*sb*m1*u2 + 32*h_l*ct*sb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1**3*s**(-2)*
     +    t1*u2 + 32*h_l*ct*sb*m1**3*s**(-2)*u2**2 + 64*h_l*ct*sb*m1**3
     +    *s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1**3*s**(-1)*u1**(-1)*
     +    u2**2 + 64*h_l*ct*sb*m1**3*s**(-1)*u2 + 64*h_l*ct*sb*m1**3*
     +    u1**(-1)*u2 + 16*h_l*ct*sb*m1**5*s**(-1)*u1**(-1)*u2 - 96*h_r
     +    *st*cb*m1*m2**2*s**(-2)*t1*u2 - 64*h_r*st*cb*m1*m2**2*s**(-2)
     +    *t1**2 - 32*h_r*st*cb*m1*m2**2*s**(-2)*u2**2 - 192*h_r*st*cb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_r*st*cb*m1*m2**2*
     +    s**(-1)*t1 - 128*h_r*st*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*st*cb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*cb*
     +    m1*m2**2*s**(-1)*u2 - 160*h_r*st*cb*m1*m2**2*t1*u1**(-1) - 96
     +    *h_r*st*cb*m1*m2**2*u1**(-1)*u2 + 32*h_r*st*cb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 + 32*h_r*st*cb*m1*s**(-1)*t1*u2 + 32*h_r*st*cb*m1*s**(-1)*
     +    u2**2 + 80*h_r*st*cb*m1*s*u1**(-1)*u2 + 64*h_r*st*cb*m1*t1*
     +    u1**(-1)*u2 + 64*h_r*st*cb*m1*u1**(-1)*u2**2 + 64*h_r*st*cb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_r*st*cb*m1**3*m2**2*s**(-1)*
     +    t1*u1**(-1) - 32*h_r*st*cb*m1**3*s**(-2)*t1*u2 - 32*h_r*st*cb
     +    *m1**3*s**(-2)*u2**2 - 64*h_r*st*cb*m1**3*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_r*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*
     +    cb*m1**3*s**(-1)*u2 - 64*h_r*st*cb*m1**3*u1**(-1)*u2 - 16*h_r
     +    *st*cb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 96*h_l*ct*sb*m1*m2**2*s**(-2)*t1*u2 - 64
     +    *h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*ct*sb*m1*m2**2*
     +    s**(-2)*u2**2 - 192*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     - 128*h_l*ct*sb*m1*m2**2*s**(-1)*t1 - 128*h_l*ct*sb*m1*m2**2
     +    *s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2 - 160*h_l*
     +    ct*sb*m1*m2**2*t1*u1**(-1) - 96*h_l*ct*sb*m1*m2**2*u1**(-1)*
     +    u2 + 32*h_l*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb
     +    *m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*sb*m1*s**(-1)*t1*u2
     +     + 32*h_l*ct*sb*m1*s**(-1)*u2**2 + 80*h_l*ct*sb*m1*s*u1**(-1)
     +    *u2 + 64*h_l*ct*sb*m1*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*
     +    u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*u2 - 32*h_l*ct*sb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*ct*sb*m1**3*s**(-2)*t1*u2
     +     - 32*h_l*ct*sb*m1**3*s**(-2)*u2**2 - 64*h_l*ct*sb*m1**3*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*sb*m1**3*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_l*ct*sb*m1**3*s**(-1)*u2 - 64*h_l*ct*sb*m1**3*
     +    u1**(-1)*u2 - 16*h_l*ct*sb*m1**5*s**(-1)*u1**(-1)*u2 + 96*h_r
     +    *st*cb*m1*m2**2*s**(-2)*t1*u2 + 64*h_r*st*cb*m1*m2**2*s**(-2)
     +    *t1**2 + 32*h_r*st*cb*m1*m2**2*s**(-2)*u2**2 + 192*h_r*st*cb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_r*st*cb*m1*m2**2*
     +    s**(-1)*t1 + 128*h_r*st*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*st*cb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*
     +    m1*m2**2*s**(-1)*u2 + 160*h_r*st*cb*m1*m2**2*t1*u1**(-1) + 96
     +    *h_r*st*cb*m1*m2**2*u1**(-1)*u2 - 32*h_r*st*cb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 - 32*h_r*st*cb*m1*s**(-1)*t1*u2 - 32*h_r*st*cb*m1*s**(-1)*
     +    u2**2 - 80*h_r*st*cb*m1*s*u1**(-1)*u2 - 64*h_r*st*cb*m1*t1*
     +    u1**(-1)*u2 - 64*h_r*st*cb*m1*u1**(-1)*u2**2 - 64*h_r*st*cb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*st*cb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*st*cb*m1**3*s**(-2)*t1*u2 + 32*h_r*st*cb*
     +    m1**3*s**(-2)*u2**2 + 64*h_r*st*cb*m1**3*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_r*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb
     +    *m1**3*s**(-1)*u2 + 64*h_r*st*cb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    st*cb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 32*h_l*ct*sb*m1*s**(-2)*t1 + 32*h_l*ct*sb*m1*s**(-2)*u2 + 
     +    96*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 112*h_l*ct*sb*m1*
     +    s**(-1)*u1**(-1)*u2 + 96*h_l*ct*sb*m1*s*u1**(-2) + 96*h_l*ct*
     +    sb*m1*t1*u1**(-2) + 96*h_l*ct*sb*m1*u1**(-2)*u2 + 96*h_l*ct*
     +    sb*m1*u1**(-1) + 32*h_l*ct*sb*m1**3*s**(-1)*u1**(-1) + 64*h_l
     +    *ct*sb*m1**3*u1**(-2) - 32*h_r*st*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) - 32*h_r*st*cb*m1*s**(-2)*t1 - 32*h_r*st*cb
     +    *m1*s**(-2)*u2 - 96*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) - 112*
     +    h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*st*cb*m1*s*u1**(-2)
     +     - 96*h_r*st*cb*m1*t1*u1**(-2) - 96*h_r*st*cb*m1*u1**(-2)*u2
     +     - 96*h_r*st*cb*m1*u1**(-1) - 32*h_r*st*cb*m1**3*s**(-1)*
     +    u1**(-1) - 64*h_r*st*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 32*h_l*ct*sb*m1*s**(-2)*t1 - 32*h_l*ct*sb*m1*s**(-2)*u2 - 
     +    96*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) - 112*h_l*ct*sb*m1*
     +    s**(-1)*u1**(-1)*u2 - 96*h_l*ct*sb*m1*s*u1**(-2) - 96*h_l*ct*
     +    sb*m1*t1*u1**(-2) - 96*h_l*ct*sb*m1*u1**(-2)*u2 - 96*h_l*ct*
     +    sb*m1*u1**(-1) - 32*h_l*ct*sb*m1**3*s**(-1)*u1**(-1) - 64*h_l
     +    *ct*sb*m1**3*u1**(-2) + 32*h_r*st*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 32*h_r*st*cb*m1*s**(-2)*t1 + 32*h_r*st*cb
     +    *m1*s**(-2)*u2 + 96*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) + 112*
     +    h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*st*cb*m1*s*u1**(-2)
     +     + 96*h_r*st*cb*m1*t1*u1**(-2) + 96*h_r*st*cb*m1*u1**(-2)*u2
     +     + 96*h_r*st*cb*m1*u1**(-1) + 32*h_r*st*cb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*st*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1*s**(-1)*t1*u1 + 32*
     +    h_l*ct*sb*m1*s**(-1)*t1**2 + 32*h_l*ct*sb*m1*s**(-1)*u1**2 - 
     +    128*h_l*ct*sb*m1**3 - 64*h_r*st*cb*m1*s**(-1)*t1*u1 - 32*h_r*
     +    st*cb*m1*s**(-1)*t1**2 - 32*h_r*st*cb*m1*s**(-1)*u1**2 + 128*
     +    h_r*st*cb*m1**3 )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*sb*m1*u1**(-1) - 32*h_r*st*cb*m1*u1**(-1)
     +     )
      MM_s = MM_s + SCB(3,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*sb*m1*s**(-1) + 32*h_l*ct*sb*m1*s*
     +    t2**(-2) + 32*h_l*ct*sb*m1*s**2*u1**(-1)*t2**(-2) - 64*h_l*ct
     +    *sb*m1*u1**(-1) - 32*h_l*ct*sb*m1*t2**(-1) + 32*h_r*st*cb*m1*
     +    s**(-1) - 32*h_r*st*cb*m1*s*t2**(-2) - 32*h_r*st*cb*m1*s**2*
     +    u1**(-1)*t2**(-2) + 64*h_r*st*cb*m1*u1**(-1) + 32*h_r*st*cb*
     +    m1*t2**(-1) )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_l*st*cb*m1*m2**2*s**(-2)*t1*u2
     +     + 64*h_l*st*cb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*st*cb*m1*
     +    m2**2*s**(-2)*u2**2 + 192*h_l*st*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 128*h_l*st*cb*m1*m2**2*s**(-1)*t1 + 128*h_l*st*
     +    cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*cb*m1*m2**2*s**(-1)*u2 + 
     +    160*h_l*st*cb*m1*m2**2*t1*u1**(-1) + 96*h_l*st*cb*m1*m2**2*
     +    u1**(-1)*u2 - 32*h_l*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*st*cb*m1*
     +    s**(-1)*t1*u2 - 32*h_l*st*cb*m1*s**(-1)*u2**2 - 80*h_l*st*cb*
     +    m1*s*u1**(-1)*u2 - 64*h_l*st*cb*m1*t1*u1**(-1)*u2 - 64*h_l*st
     +    *cb*m1*u1**(-1)*u2**2 - 64*h_l*st*cb*m1*u2 + 32*h_l*st*cb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1**3*s**(-2)*
     +    t1*u2 + 32*h_l*st*cb*m1**3*s**(-2)*u2**2 + 64*h_l*st*cb*m1**3
     +    *s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1**3*s**(-1)*u1**(-1)*
     +    u2**2 + 64*h_l*st*cb*m1**3*s**(-1)*u2 + 64*h_l*st*cb*m1**3*
     +    u1**(-1)*u2 + 16*h_l*st*cb*m1**5*s**(-1)*u1**(-1)*u2 - 96*h_r
     +    *ct*sb*m1*m2**2*s**(-2)*t1*u2 - 64*h_r*ct*sb*m1*m2**2*s**(-2)
     +    *t1**2 - 32*h_r*ct*sb*m1*m2**2*s**(-2)*u2**2 - 192*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_r*ct*sb*m1*m2**2*
     +    s**(-1)*t1 - 128*h_r*ct*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*ct*sb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*u2 - 160*h_r*ct*sb*m1*m2**2*t1*u1**(-1) - 96
     +    *h_r*ct*sb*m1*m2**2*u1**(-1)*u2 + 32*h_r*ct*sb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 + 32*h_r*ct*sb*m1*s**(-1)*t1*u2 + 32*h_r*ct*sb*m1*s**(-1)*
     +    u2**2 + 80*h_r*ct*sb*m1*s*u1**(-1)*u2 + 64*h_r*ct*sb*m1*t1*
     +    u1**(-1)*u2 + 64*h_r*ct*sb*m1*u1**(-1)*u2**2 + 64*h_r*ct*sb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_r*ct*sb*m1**3*m2**2*s**(-1)*
     +    t1*u1**(-1) - 32*h_r*ct*sb*m1**3*s**(-2)*t1*u2 - 32*h_r*ct*sb
     +    *m1**3*s**(-2)*u2**2 - 64*h_r*ct*sb*m1**3*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_r*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*
     +    sb*m1**3*s**(-1)*u2 - 64*h_r*ct*sb*m1**3*u1**(-1)*u2 - 16*h_r
     +    *ct*sb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 96*h_l*st*cb*m1*m2**2*s**(-2)*t1*u2 - 64
     +    *h_l*st*cb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*st*cb*m1*m2**2*
     +    s**(-2)*u2**2 - 192*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     - 128*h_l*st*cb*m1*m2**2*s**(-1)*t1 - 128*h_l*st*cb*m1*m2**2
     +    *s**(-1)*t1**2*u1**(-1) - 64*h_l*st*cb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*st*cb*m1*m2**2*s**(-1)*u2 - 160*h_l*
     +    st*cb*m1*m2**2*t1*u1**(-1) - 96*h_l*st*cb*m1*m2**2*u1**(-1)*
     +    u2 + 32*h_l*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_l*st*cb
     +    *m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*st*cb*m1*s**(-1)*t1*u2
     +     + 32*h_l*st*cb*m1*s**(-1)*u2**2 + 80*h_l*st*cb*m1*s*u1**(-1)
     +    *u2 + 64*h_l*st*cb*m1*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*
     +    u1**(-1)*u2**2 + 64*h_l*st*cb*m1*u2 - 32*h_l*st*cb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*st*cb*m1**3*s**(-2)*t1*u2
     +     - 32*h_l*st*cb*m1**3*s**(-2)*u2**2 - 64*h_l*st*cb*m1**3*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*cb*m1**3*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_l*st*cb*m1**3*s**(-1)*u2 - 64*h_l*st*cb*m1**3*
     +    u1**(-1)*u2 - 16*h_l*st*cb*m1**5*s**(-1)*u1**(-1)*u2 + 96*h_r
     +    *ct*sb*m1*m2**2*s**(-2)*t1*u2 + 64*h_r*ct*sb*m1*m2**2*s**(-2)
     +    *t1**2 + 32*h_r*ct*sb*m1*m2**2*s**(-2)*u2**2 + 192*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_r*ct*sb*m1*m2**2*
     +    s**(-1)*t1 + 128*h_r*ct*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*sb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*u2 + 160*h_r*ct*sb*m1*m2**2*t1*u1**(-1) + 96
     +    *h_r*ct*sb*m1*m2**2*u1**(-1)*u2 - 32*h_r*ct*sb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 - 32*h_r*ct*sb*m1*s**(-1)*t1*u2 - 32*h_r*ct*sb*m1*s**(-1)*
     +    u2**2 - 80*h_r*ct*sb*m1*s*u1**(-1)*u2 - 64*h_r*ct*sb*m1*t1*
     +    u1**(-1)*u2 - 64*h_r*ct*sb*m1*u1**(-1)*u2**2 - 64*h_r*ct*sb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*ct*sb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*ct*sb*m1**3*s**(-2)*t1*u2 + 32*h_r*ct*sb*
     +    m1**3*s**(-2)*u2**2 + 64*h_r*ct*sb*m1**3*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_r*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb
     +    *m1**3*s**(-1)*u2 + 64*h_r*ct*sb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    ct*sb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 32*h_l*st*cb*m1*s**(-2)*t1 + 32*h_l*st*cb*m1*s**(-2)*u2 + 
     +    96*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 112*h_l*st*cb*m1*
     +    s**(-1)*u1**(-1)*u2 + 96*h_l*st*cb*m1*s*u1**(-2) + 96*h_l*st*
     +    cb*m1*t1*u1**(-2) + 96*h_l*st*cb*m1*u1**(-2)*u2 + 96*h_l*st*
     +    cb*m1*u1**(-1) + 32*h_l*st*cb*m1**3*s**(-1)*u1**(-1) + 64*h_l
     +    *st*cb*m1**3*u1**(-2) - 32*h_r*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) - 32*h_r*ct*sb*m1*s**(-2)*t1 - 32*h_r*ct*sb
     +    *m1*s**(-2)*u2 - 96*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) - 112*
     +    h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*sb*m1*s*u1**(-2)
     +     - 96*h_r*ct*sb*m1*t1*u1**(-2) - 96*h_r*ct*sb*m1*u1**(-2)*u2
     +     - 96*h_r*ct*sb*m1*u1**(-1) - 32*h_r*ct*sb*m1**3*s**(-1)*
     +    u1**(-1) - 64*h_r*ct*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 32*h_l*st*cb*m1*s**(-2)*t1 - 32*h_l*st*cb*m1*s**(-2)*u2 - 
     +    96*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) - 112*h_l*st*cb*m1*
     +    s**(-1)*u1**(-1)*u2 - 96*h_l*st*cb*m1*s*u1**(-2) - 96*h_l*st*
     +    cb*m1*t1*u1**(-2) - 96*h_l*st*cb*m1*u1**(-2)*u2 - 96*h_l*st*
     +    cb*m1*u1**(-1) - 32*h_l*st*cb*m1**3*s**(-1)*u1**(-1) - 64*h_l
     +    *st*cb*m1**3*u1**(-2) + 32*h_r*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 32*h_r*ct*sb*m1*s**(-2)*t1 + 32*h_r*ct*sb
     +    *m1*s**(-2)*u2 + 96*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) + 112*
     +    h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*sb*m1*s*u1**(-2)
     +     + 96*h_r*ct*sb*m1*t1*u1**(-2) + 96*h_r*ct*sb*m1*u1**(-2)*u2
     +     + 96*h_r*ct*sb*m1*u1**(-1) + 32*h_r*ct*sb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*ct*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1*s**(-1)*t1*u1 + 32*
     +    h_l*st*cb*m1*s**(-1)*t1**2 + 32*h_l*st*cb*m1*s**(-1)*u1**2 - 
     +    128*h_l*st*cb*m1**3 - 64*h_r*ct*sb*m1*s**(-1)*t1*u1 - 32*h_r*
     +    ct*sb*m1*s**(-1)*t1**2 - 32*h_r*ct*sb*m1*s**(-1)*u1**2 + 128*
     +    h_r*ct*sb*m1**3 )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*u1**(-1) - 32*h_r*ct*sb*m1*u1**(-1)
     +     )
      MM_s = MM_s + SCB(3,4)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*s**(-1) + 32*h_l*st*cb*m1*s*
     +    t2**(-2) + 32*h_l*st*cb*m1*s**2*u1**(-1)*t2**(-2) - 64*h_l*st
     +    *cb*m1*u1**(-1) - 32*h_l*st*cb*m1*t2**(-1) + 32*h_r*ct*sb*m1*
     +    s**(-1) - 32*h_r*ct*sb*m1*s*t2**(-2) - 32*h_r*ct*sb*m1*s**2*
     +    u1**(-1)*t2**(-2) + 64*h_r*ct*sb*m1*u1**(-1) + 32*h_r*ct*sb*
     +    m1*t2**(-1) )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 96*h_l*st*sb*m1*m2**2*s**(-2)*t1*
     +    u2 - 64*h_l*st*sb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*st*sb*m1*
     +    m2**2*s**(-2)*u2**2 - 192*h_l*st*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 128*h_l*st*sb*m1*m2**2*s**(-1)*t1 - 128*h_l*st*
     +    sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*sb*m1*m2**2*s**(-1)*u2 - 
     +    160*h_l*st*sb*m1*m2**2*t1*u1**(-1) - 96*h_l*st*sb*m1*m2**2*
     +    u1**(-1)*u2 + 32*h_l*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*st*sb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*st*sb*m1*
     +    s**(-1)*t1*u2 + 32*h_l*st*sb*m1*s**(-1)*u2**2 + 80*h_l*st*sb*
     +    m1*s*u1**(-1)*u2 + 64*h_l*st*sb*m1*t1*u1**(-1)*u2 + 64*h_l*st
     +    *sb*m1*u1**(-1)*u2**2 + 64*h_l*st*sb*m1*u2 - 32*h_l*st*sb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*st*sb*m1**3*s**(-2)*
     +    t1*u2 - 32*h_l*st*sb*m1**3*s**(-2)*u2**2 - 64*h_l*st*sb*m1**3
     +    *s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*st*sb*m1**3*s**(-1)*u2 - 64*h_l*st*sb
     +    *m1**3*u1**(-1)*u2 - 16*h_l*st*sb*m1**5*s**(-1)*u1**(-1)*u2
     +     - 96*h_r*ct*cb*m1*m2**2*s**(-2)*t1*u2 - 64*h_r*ct*cb*m1*
     +    m2**2*s**(-2)*t1**2 - 32*h_r*ct*cb*m1*m2**2*s**(-2)*u2**2 - 
     +    192*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_r*ct*cb
     +    *m1*m2**2*s**(-1)*t1 - 128*h_r*ct*cb*m1*m2**2*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_r*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*
     +    h_r*ct*cb*m1*m2**2*s**(-1)*u2 - 160*h_r*ct*cb*m1*m2**2*t1*
     +    u1**(-1) - 96*h_r*ct*cb*m1*m2**2*u1**(-1)*u2 + 32*h_r*ct*cb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2 + 32*h_r*ct*cb*m1*s**(-1)*t1*u2 + 32*h_r*ct*cb*m1
     +    *s**(-1)*u2**2 + 80*h_r*ct*cb*m1*s*u1**(-1)*u2 + 64*h_r*ct*cb
     +    *m1*t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1*u1**(-1)*u2**2 + 64*h_r*
     +    ct*cb*m1*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_r*ct*cb*m1**3*m2**2*s**(-1)*
     +    t1*u1**(-1) - 32*h_r*ct*cb*m1**3*s**(-2)*t1*u2 - 32*h_r*ct*cb
     +    *m1**3*s**(-2)*u2**2 - 64*h_r*ct*cb*m1**3*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_r*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*
     +    cb*m1**3*s**(-1)*u2 - 64*h_r*ct*cb*m1**3*u1**(-1)*u2 - 16*h_r
     +    *ct*cb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*h_l*st*sb*m1*m2**2*s**(-2)*t1*u2 + 64*
     +    h_l*st*sb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*st*sb*m1*m2**2*
     +    s**(-2)*u2**2 + 192*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     + 128*h_l*st*sb*m1*m2**2*s**(-1)*t1 + 128*h_l*st*sb*m1*m2**2
     +    *s**(-1)*t1**2*u1**(-1) + 64*h_l*st*sb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*st*sb*m1*m2**2*s**(-1)*u2 + 160*h_l*
     +    st*sb*m1*m2**2*t1*u1**(-1) + 96*h_l*st*sb*m1*m2**2*u1**(-1)*
     +    u2 - 32*h_l*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb
     +    *m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*st*sb*m1*s**(-1)*t1*u2
     +     - 32*h_l*st*sb*m1*s**(-1)*u2**2 - 80*h_l*st*sb*m1*s*u1**(-1)
     +    *u2 - 64*h_l*st*sb*m1*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*
     +    u1**(-1)*u2**2 - 64*h_l*st*sb*m1*u2 + 32*h_l*st*sb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1**3*s**(-2)*t1*u2
     +     + 32*h_l*st*sb*m1**3*s**(-2)*u2**2 + 64*h_l*st*sb*m1**3*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*st*sb*m1**3*s**(-1)*u2 + 64*h_l*st*sb*m1**3*
     +    u1**(-1)*u2 + 16*h_l*st*sb*m1**5*s**(-1)*u1**(-1)*u2 + 96*h_r
     +    *ct*cb*m1*m2**2*s**(-2)*t1*u2 + 64*h_r*ct*cb*m1*m2**2*s**(-2)
     +    *t1**2 + 32*h_r*ct*cb*m1*m2**2*s**(-2)*u2**2 + 192*h_r*ct*cb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_r*ct*cb*m1*m2**2*
     +    s**(-1)*t1 + 128*h_r*ct*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb*
     +    m1*m2**2*s**(-1)*u2 + 160*h_r*ct*cb*m1*m2**2*t1*u1**(-1) + 96
     +    *h_r*ct*cb*m1*m2**2*u1**(-1)*u2 - 32*h_r*ct*cb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*
     +    u2 - 32*h_r*ct*cb*m1*s**(-1)*t1*u2 - 32*h_r*ct*cb*m1*s**(-1)*
     +    u2**2 - 80*h_r*ct*cb*m1*s*u1**(-1)*u2 - 64*h_r*ct*cb*m1*t1*
     +    u1**(-1)*u2 - 64*h_r*ct*cb*m1*u1**(-1)*u2**2 - 64*h_r*ct*cb*
     +    m1*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*ct*cb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*ct*cb*m1**3*s**(-2)*t1*u2 + 32*h_r*ct*cb*
     +    m1**3*s**(-2)*u2**2 + 64*h_r*ct*cb*m1**3*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_r*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb
     +    *m1**3*s**(-1)*u2 + 64*h_r*ct*cb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    ct*cb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_l*st*sb*m1*s**(-2)*t1 - 32*h_l*st*sb*m1*
     +    s**(-2)*u2 - 96*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 112*h_l*st
     +    *sb*m1*s**(-1)*u1**(-1)*u2 - 96*h_l*st*sb*m1*s*u1**(-2) - 96*
     +    h_l*st*sb*m1*t1*u1**(-2) - 96*h_l*st*sb*m1*u1**(-2)*u2 - 96*
     +    h_l*st*sb*m1*u1**(-1) - 32*h_l*st*sb*m1**3*s**(-1)*u1**(-1)
     +     - 64*h_l*st*sb*m1**3*u1**(-2) - 32*h_r*ct*cb*m1*m2**2*
     +    s**(-1)*t1*u1**(-1)*u2**(-1) - 32*h_r*ct*cb*m1*s**(-2)*t1 - 
     +    32*h_r*ct*cb*m1*s**(-2)*u2 - 96*h_r*ct*cb*m1*s**(-1)*t1*
     +    u1**(-1) - 112*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*
     +    cb*m1*s*u1**(-2) - 96*h_r*ct*cb*m1*t1*u1**(-2) - 96*h_r*ct*cb
     +    *m1*u1**(-2)*u2 - 96*h_r*ct*cb*m1*u1**(-1) - 32*h_r*ct*cb*
     +    m1**3*s**(-1)*u1**(-1) - 64*h_r*ct*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 
     +    32*h_l*st*sb*m1*s**(-2)*t1 + 32*h_l*st*sb*m1*s**(-2)*u2 + 96*
     +    h_l*st*sb*m1*s**(-1)*t1*u1**(-1) + 112*h_l*st*sb*m1*s**(-1)*
     +    u1**(-1)*u2 + 96*h_l*st*sb*m1*s*u1**(-2) + 96*h_l*st*sb*m1*t1
     +    *u1**(-2) + 96*h_l*st*sb*m1*u1**(-2)*u2 + 96*h_l*st*sb*m1*
     +    u1**(-1) + 32*h_l*st*sb*m1**3*s**(-1)*u1**(-1) + 64*h_l*st*sb
     +    *m1**3*u1**(-2) + 32*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_r*ct*cb*m1*s**(-2)*t1 + 32*h_r*ct*cb*m1*
     +    s**(-2)*u2 + 96*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1) + 112*h_r*ct
     +    *cb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*cb*m1*s*u1**(-2) + 96*
     +    h_r*ct*cb*m1*t1*u1**(-2) + 96*h_r*ct*cb*m1*u1**(-2)*u2 + 96*
     +    h_r*ct*cb*m1*u1**(-1) + 32*h_r*ct*cb*m1**3*s**(-1)*u1**(-1)
     +     + 64*h_r*ct*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1*s**(-1)*t1*u1 - 32
     +    *h_l*st*sb*m1*s**(-1)*t1**2 - 32*h_l*st*sb*m1*s**(-1)*u1**2
     +     + 128*h_l*st*sb*m1**3 - 64*h_r*ct*cb*m1*s**(-1)*t1*u1 - 32*
     +    h_r*ct*cb*m1*s**(-1)*t1**2 - 32*h_r*ct*cb*m1*s**(-1)*u1**2 + 
     +    128*h_r*ct*cb*m1**3 )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*u1**(-1) - 32*h_r*ct*cb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(3,5)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*s**(-1) - 32*h_l*st*sb*m1*s*t2**(-2)
     +     - 32*h_l*st*sb*m1*s**2*u1**(-1)*t2**(-2) + 64*h_l*st*sb*m1*
     +    u1**(-1) + 32*h_l*st*sb*m1*t2**(-1) + 32*h_r*ct*cb*m1*s**(-1)
     +     - 32*h_r*ct*cb*m1*s*t2**(-2) - 32*h_r*ct*cb*m1*s**2*u1**(-1)
     +    *t2**(-2) + 64*h_r*ct*cb*m1*u1**(-1) + 32*h_r*ct*cb*m1*
     +    t2**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qi*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16 - 16*m1**2*mg**2*s**(-1)*u1**(-1) + 16*
     +    m1**2*msb1**2*s**(-1)*u1**(-1) - 16*m1**2*u1**(-1) + 16*mg**2
     +    *s**(-2)*t1 + 16*mg**2*s**(-2)*t1**2*u1**(-1) + 16*mg**2*
     +    s**(-2)*u1 + 16*mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*s**(-1)
     +     - 16*msb1**2*s**(-2)*t1 - 16*msb1**2*s**(-2)*t1**2*u1**(-1)
     +     - 16*msb1**2*s**(-2)*u1 - 16*msb1**2*s**(-1)*t1*u1**(-1) - 
     +    16*msb1**2*s**(-1) + 16*s**(-1)*t1 + 16*s**(-1)*t1**2*
     +    u1**(-1) + 16*s**(-1)*u1 + 16*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qi*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2b*m1**2*mg**2*s**(-1)*u1**(-1) + 16
     +    *c2b*m1**2*msb1**2*s**(-1)*u1**(-1) - 16*c2b*m1**2*u1**(-1)
     +     + 16*c2b*mg**2*s**(-2)*t1 + 16*c2b*mg**2*s**(-2)*t1**2*
     +    u1**(-1) + 16*c2b*mg**2*s**(-2)*u1 + 16*c2b*mg**2*s**(-1)*t1*
     +    u1**(-1) + 16*c2b*mg**2*s**(-1) - 16*c2b*msb1**2*s**(-2)*t1
     +     - 16*c2b*msb1**2*s**(-2)*t1**2*u1**(-1) - 16*c2b*msb1**2*
     +    s**(-2)*u1 - 16*c2b*msb1**2*s**(-1)*t1*u1**(-1) - 16*c2b*
     +    msb1**2*s**(-1) + 16*c2b*s**(-1)*t1 + 16*c2b*s**(-1)*t1**2*
     +    u1**(-1) + 16*c2b*s**(-1)*u1 + 16*c2b*t1*u1**(-1) + 16*c2b )
      MM_s = MM_s + SCB(4,2)*susy_qi*Nc*Cf**2*Pi**2*alphas**2*prefac
     +  * ( 64*h_l*h_r*s2b*m1*mg*s**(-1)*t1*u1**(-1) + 64*h_l*h_r*s2b*
     +    m1*mg*s**(-1) + 128*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 16*(h_r2*sb2+h_l2*cb2)*m1**2*mg**2*s**(-1)*u1**(-1) + 16
     +    *(h_r2*sb2+h_l2*cb2)*m1**2*msb1**2*s**(-1)*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*u1**(-1) + 16*(h_r2*sb2+h_l2*cb2)*
     +    mg**2*s**(-2)*t1 - 16*(h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*t1**2
     +    *u1**(-1) + 16*(h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*u1 - 16*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-1) - 16*(h_r2*sb2+h_l2*cb2)*
     +    msb1**2*s**(-2)*t1 + 16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*
     +    t1**2*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*u1 + 
     +    16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1) - 16*(h_r2*sb2+h_l2*cb2)*
     +    s**(-1)*t1 - 16*(h_r2*sb2+h_l2*cb2)*s**(-1)*t1**2*u1**(-1) - 
     +    16*(h_r2*sb2+h_l2*cb2)*s**(-1)*u1 - 16*(h_r2*sb2+h_l2*cb2)*t1
     +    *u1**(-1) - 16*(h_r2*sb2+h_l2*cb2) - 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qqg*Cf*Pi**2*alphas**2*prefac * ( 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*mg**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*msb1**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*m1**2*u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*
     +    mg**2*s**(-2)*t1 + 16*(h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*t1**2
     +    *u1**(-1) - 16*(h_r2*sb2+h_l2*cb2)*mg**2*s**(-2)*u1 + 16*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-1) + 16*(h_r2*sb2+h_l2*cb2)*
     +    msb1**2*s**(-2)*t1 - 16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*
     +    t1**2*u1**(-1) + 16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-2)*u1 - 
     +    16*(h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1) + 16*(h_r2*sb2+h_l2*cb2)*
     +    s**(-1)*t1 + 16*(h_r2*sb2+h_l2*cb2)*s**(-1)*t1**2*u1**(-1) + 
     +    16*(h_r2*sb2+h_l2*cb2)*s**(-1)*u1 + 16*(h_r2*sb2+h_l2*cb2)*t1
     +    *u1**(-1) + 16*(h_r2*sb2+h_l2*cb2) + 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*susy_qqg*Cf*Pi**2*alphas**2*prefac * ( 32*
     +    h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*cb*m1*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*ct*cb*m1*
     +    m2**2*s**(-2)*u2**2 + 128*h_l*ct*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*t1 + 64*h_l*ct*cb
     +    *m1*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2 + 
     +    64*h_l*ct*cb*m1*m2**2*t1*u1**(-1) + 96*h_l*ct*cb*m1*m2**2*
     +    u1**(-1)*u2 - 16*h_l*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*m1*
     +    s**(-1)*t1**2 - 32*h_l*ct*cb*m1*s**(-1)*u2**2 + 80*h_l*ct*cb*
     +    m1*s*t1*u1**(-1) - 80*h_l*ct*cb*m1*s*u1**(-1)*u2 + 64*h_l*ct*
     +    cb*m1*t1 + 64*h_l*ct*cb*m1*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*
     +    u1**(-1)*u2**2 - 64*h_l*ct*cb*m1*u2 + 32*h_l*ct*cb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1**3*s**(-2)*t1*u2
     +     - 32*h_l*ct*cb*m1**3*s**(-2)*t1**2 - 32*h_l*ct*cb*m1**3*
     +    s**(-2)*u2**2 - 128*h_l*ct*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_l*ct*cb*m1**3*s**(-1)*t1 - 64*h_l*ct*cb*m1**3*s**(-1)*
     +    t1**2*u1**(-1) - 64*h_l*ct*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 
     +    64*h_l*ct*cb*m1**3*s**(-1)*u2 - 64*h_l*ct*cb*m1**3*t1*
     +    u1**(-1) - 96*h_l*ct*cb*m1**3*u1**(-1)*u2 - 16*h_l*ct*cb*
     +    m1**5*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 + 64*h_r*st*sb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*st
     +    *sb*m1*m2**2*s**(-2)*t1**2 + 32*h_r*st*sb*m1*m2**2*s**(-2)*
     +    u2**2 + 128*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*st*sb*m1*m2**2*s**(-1)*t1 + 64*h_r*st*sb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) + 64*h_r*st*sb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 + 64*h_r*st*sb*m1*m2**2*s**(-1)*u2 + 64*h_r*st*sb*m1*
     +    m2**2*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_r*st*sb*m1*m2**2*u1**(-1)*u2 - 
     +    16*h_r*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*
     +    m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*m1*s**(-1)*t1**2 - 
     +    32*h_r*st*sb*m1*s**(-1)*u2**2 + 80*h_r*st*sb*m1*s*t1*u1**(-1)
     +     - 80*h_r*st*sb*m1*s*u1**(-1)*u2 + 64*h_r*st*sb*m1*t1 + 64*
     +    h_r*st*sb*m1*t1**2*u1**(-1) - 64*h_r*st*sb*m1*u1**(-1)*u2**2
     +     - 64*h_r*st*sb*m1*u2 + 32*h_r*st*sb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*st*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 - 64*
     +    h_r*st*sb*m1**3*s**(-2)*t1*u2 - 32*h_r*st*sb*m1**3*s**(-2)*
     +    t1**2 - 32*h_r*st*sb*m1**3*s**(-2)*u2**2 - 128*h_r*st*sb*
     +    m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1**3*s**(-1)*t1
     +     - 64*h_r*st*sb*m1**3*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb*
     +    m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*sb*m1**3*s**(-1)*u2
     +     - 64*h_r*st*sb*m1**3*t1*u1**(-1) - 96*h_r*st*sb*m1**3*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*st*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*st*sb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*cb*m1*m2**2*s**(-2)*t1*u2 - 32
     +    *h_l*ct*cb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*ct*cb*m1*m2**2*
     +    s**(-2)*u2**2 - 128*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_l*ct*cb*m1*m2**2*s**(-1)*t1 - 64*h_l*ct*cb*m1*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2 - 64*h_l*ct
     +    *cb*m1*m2**2*t1*u1**(-1) - 96*h_l*ct*cb*m1*m2**2*u1**(-1)*u2
     +     + 16*h_l*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*m1*s**(-1)*t1**2
     +     + 32*h_l*ct*cb*m1*s**(-1)*u2**2 - 80*h_l*ct*cb*m1*s*t1*
     +    u1**(-1) + 80*h_l*ct*cb*m1*s*u1**(-1)*u2 - 64*h_l*ct*cb*m1*t1
     +     - 64*h_l*ct*cb*m1*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*u1**(-1)*
     +    u2**2 + 64*h_l*ct*cb*m1*u2 - 32*h_l*ct*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1) - 32*h_l*ct*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 64*h_l*ct*cb*m1**3*s**(-2)*t1*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*ct*cb*m1**3*s**(-2)*t1**2 + 32*h_l*
     +    ct*cb*m1**3*s**(-2)*u2**2 + 128*h_l*ct*cb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*s**(-1)*t1 + 64*h_l*ct*cb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*ct*cb*m1**3*s**(-1)*u2 + 64*h_l*ct*cb
     +    *m1**3*t1*u1**(-1) + 96*h_l*ct*cb*m1**3*u1**(-1)*u2 + 16*h_l*
     +    ct*cb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 - 64*h_r*st*sb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*st
     +    *sb*m1*m2**2*s**(-2)*t1**2 - 32*h_r*st*sb*m1*m2**2*s**(-2)*
     +    u2**2 - 128*h_r*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_r*st*sb*m1*m2**2*s**(-1)*t1 - 64*h_r*st*sb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) - 64*h_r*st*sb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_r*st*sb*m1*m2**2*s**(-1)*u2 - 64*h_r*st*sb*m1*
     +    m2**2*t1*u1**(-1) - 96*h_r*st*sb*m1*m2**2*u1**(-1)*u2 + 16*
     +    h_r*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_r*st*sb*m1*m2**4*s**(-1)*u1**(-1)*u2
     +     - 32*h_r*st*sb*m1*s**(-1)*t1**2 + 32*h_r*st*sb*m1*s**(-1)*
     +    u2**2 - 80*h_r*st*sb*m1*s*t1*u1**(-1) + 80*h_r*st*sb*m1*s*
     +    u1**(-1)*u2 - 64*h_r*st*sb*m1*t1 - 64*h_r*st*sb*m1*t1**2*
     +    u1**(-1) + 64*h_r*st*sb*m1*u1**(-1)*u2**2 + 64*h_r*st*sb*m1*
     +    u2 - 32*h_r*st*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_r*st
     +    *sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*st*sb*m1**3*
     +    s**(-2)*t1*u2 + 32*h_r*st*sb*m1**3*s**(-2)*t1**2 + 32*h_r*st*
     +    sb*m1**3*s**(-2)*u2**2 + 128*h_r*st*sb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_r*st*sb*m1**3*s**(-1)*t1 + 64*h_r*st*sb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*st*sb*m1**3*s**(-1)*u2 + 64*h_r*st*sb
     +    *m1**3*t1*u1**(-1) + 96*h_r*st*sb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    st*sb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*
     +    m1*s**(-1)*u1**(-1)*u2 + 16*h_r*st*sb*m1*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb
     +    *m1*s**(-1)*u1**(-1)*u2 - 16*h_r*st*sb*m1*s**(-1)*t1*u1**(-1)
     +     - 16*h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*cb*m1*u1**(-1) + 32*h_r*st*sb*m1*u1**(-1)
     +     )
      MM_s = MM_s + SCB(4,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1*u1**(-1) - 32*h_r*st*sb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*cb*m1*m2**2*s**(-2)*t1*
     +    u2 - 32*h_l*st*cb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*st*cb*m1*
     +    m2**2*s**(-2)*u2**2 - 128*h_l*st*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_l*st*cb*m1*m2**2*s**(-1)*t1 - 64*h_l*st*cb
     +    *m1*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*cb*m1*m2**2*s**(-1)*u2 - 
     +    64*h_l*st*cb*m1*m2**2*t1*u1**(-1) - 96*h_l*st*cb*m1*m2**2*
     +    u1**(-1)*u2 + 16*h_l*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*st*cb*m1*
     +    s**(-1)*t1**2 + 32*h_l*st*cb*m1*s**(-1)*u2**2 - 80*h_l*st*cb*
     +    m1*s*t1*u1**(-1) + 80*h_l*st*cb*m1*s*u1**(-1)*u2 - 64*h_l*st*
     +    cb*m1*t1 - 64*h_l*st*cb*m1*t1**2*u1**(-1) + 64*h_l*st*cb*m1*
     +    u1**(-1)*u2**2 + 64*h_l*st*cb*m1*u2 - 32*h_l*st*cb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*st*cb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1**3*s**(-2)*t1*u2 + 32
     +    *h_l*st*cb*m1**3*s**(-2)*t1**2 + 32*h_l*st*cb*m1**3*s**(-2)*
     +    u2**2 + 128*h_l*st*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*
     +    st*cb*m1**3*s**(-1)*t1 + 64*h_l*st*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_l
     +    *st*cb*m1**3*s**(-1)*u2 + 64*h_l*st*cb*m1**3*t1*u1**(-1) + 96
     +    *h_l*st*cb*m1**3*u1**(-1)*u2 + 16*h_l*st*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*st*cb*m1**5*s**(-1)*u1**(-1)*u2 + 64*h_r*ct
     +    *sb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*ct*sb*m1*m2**2*s**(-2)*
     +    t1**2 + 32*h_r*ct*sb*m1*m2**2*s**(-2)*u2**2 + 128*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1*m2**2*
     +    s**(-1)*t1 + 64*h_r*ct*sb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*sb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*u2 + 64*h_r*ct*sb*m1*m2**2*t1*u1**(-1) + 96*
     +    h_r*ct*sb*m1*m2**2*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*ct*sb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *ct*sb*m1*s**(-1)*t1**2 - 32*h_r*ct*sb*m1*s**(-1)*u2**2 + 80*
     +    h_r*ct*sb*m1*s*t1*u1**(-1) - 80*h_r*ct*sb*m1*s*u1**(-1)*u2 + 
     +    64*h_r*ct*sb*m1*t1 + 64*h_r*ct*sb*m1*t1**2*u1**(-1) - 64*h_r*
     +    ct*sb*m1*u1**(-1)*u2**2 - 64*h_r*ct*sb*m1*u2 + 32*h_r*ct*sb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) + 32*h_r*ct*sb*m1**3*m2**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*s**(-2)*t1*u2 - 32*
     +    h_r*ct*sb*m1**3*s**(-2)*t1**2 - 32*h_r*ct*sb*m1**3*s**(-2)*
     +    u2**2 - 128*h_r*ct*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    ct*sb*m1**3*s**(-1)*t1 - 64*h_r*ct*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_r*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r
     +    *ct*sb*m1**3*s**(-1)*u2 - 64*h_r*ct*sb*m1**3*t1*u1**(-1) - 96
     +    *h_r*ct*sb*m1**3*u1**(-1)*u2 - 16*h_r*ct*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*ct*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*cb*m1*m2**2*s**(-2)*t1*u2 + 32*
     +    h_l*st*cb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*st*cb*m1*m2**2*
     +    s**(-2)*u2**2 + 128*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     + 64*h_l*st*cb*m1*m2**2*s**(-1)*t1 + 64*h_l*st*cb*m1*m2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_l*st*cb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*st*cb*m1*m2**2*s**(-1)*u2 + 64*h_l*st
     +    *cb*m1*m2**2*t1*u1**(-1) + 96*h_l*st*cb*m1*m2**2*u1**(-1)*u2
     +     - 16*h_l*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*st*cb*m1*s**(-1)*t1**2
     +     - 32*h_l*st*cb*m1*s**(-1)*u2**2 + 80*h_l*st*cb*m1*s*t1*
     +    u1**(-1) - 80*h_l*st*cb*m1*s*u1**(-1)*u2 + 64*h_l*st*cb*m1*t1
     +     + 64*h_l*st*cb*m1*t1**2*u1**(-1) - 64*h_l*st*cb*m1*u1**(-1)*
     +    u2**2 - 64*h_l*st*cb*m1*u2 + 32*h_l*st*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1) + 32*h_l*st*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 64*h_l*st*cb*m1**3*s**(-2)*t1*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*st*cb*m1**3*s**(-2)*t1**2 - 32*
     +    h_l*st*cb*m1**3*s**(-2)*u2**2 - 128*h_l*st*cb*m1**3*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*s**(-1)*t1 - 64*h_l*st*cb
     +    *m1**3*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*cb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*st*cb*m1**3*s**(-1)*u2 - 64*h_l*st*cb
     +    *m1**3*t1*u1**(-1) - 96*h_l*st*cb*m1**3*u1**(-1)*u2 - 16*h_l*
     +    st*cb*m1**5*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 - 64*h_r*ct*sb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*ct
     +    *sb*m1*m2**2*s**(-2)*t1**2 - 32*h_r*ct*sb*m1*m2**2*s**(-2)*
     +    u2**2 - 128*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_r*ct*sb*m1*m2**2*s**(-1)*t1 - 64*h_r*ct*sb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) - 64*h_r*ct*sb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_r*ct*sb*m1*m2**2*s**(-1)*u2 - 64*h_r*ct*sb*m1*
     +    m2**2*t1*u1**(-1) - 96*h_r*ct*sb*m1*m2**2*u1**(-1)*u2 + 16*
     +    h_r*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_r*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2
     +     - 32*h_r*ct*sb*m1*s**(-1)*t1**2 + 32*h_r*ct*sb*m1*s**(-1)*
     +    u2**2 - 80*h_r*ct*sb*m1*s*t1*u1**(-1) + 80*h_r*ct*sb*m1*s*
     +    u1**(-1)*u2 - 64*h_r*ct*sb*m1*t1 - 64*h_r*ct*sb*m1*t1**2*
     +    u1**(-1) + 64*h_r*ct*sb*m1*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1*
     +    u2 - 32*h_r*ct*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_r*ct
     +    *sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*ct*sb*m1**3*
     +    s**(-2)*t1*u2 + 32*h_r*ct*sb*m1**3*s**(-2)*t1**2 + 32*h_r*ct*
     +    sb*m1**3*s**(-2)*u2**2 + 128*h_r*ct*sb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_r*ct*sb*m1**3*s**(-1)*t1 + 64*h_r*ct*sb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*sb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*ct*sb*m1**3*s**(-1)*u2 + 64*h_r*ct*sb
     +    *m1**3*t1*u1**(-1) + 96*h_r*ct*sb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    ct*sb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) - 16*h_l*st*
     +    cb*m1*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*sb*m1*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 16*h_l*st*cb*m1
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) - 
     +    16*h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*u1**(-1) + 32*h_r*ct*sb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(4,2)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*u1**(-1) - 32*h_r*ct*sb*m1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qi*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16 - 16*m1**2*mg**2*s**(-1)*u1**(-1) + 16*
     +    m1**2*msb2**2*s**(-1)*u1**(-1) - 16*m1**2*u1**(-1) + 16*mg**2
     +    *s**(-2)*t1 + 16*mg**2*s**(-2)*t1**2*u1**(-1) + 16*mg**2*
     +    s**(-2)*u1 + 16*mg**2*s**(-1)*t1*u1**(-1) + 16*mg**2*s**(-1)
     +     - 16*msb2**2*s**(-2)*t1 - 16*msb2**2*s**(-2)*t1**2*u1**(-1)
     +     - 16*msb2**2*s**(-2)*u1 - 16*msb2**2*s**(-1)*t1*u1**(-1) - 
     +    16*msb2**2*s**(-1) + 16*s**(-1)*t1 + 16*s**(-1)*t1**2*
     +    u1**(-1) + 16*s**(-1)*u1 + 16*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qi*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2b*m1**2*mg**2*s**(-1)*u1**(-1) - 16*
     +    c2b*m1**2*msb2**2*s**(-1)*u1**(-1) + 16*c2b*m1**2*u1**(-1) - 
     +    16*c2b*mg**2*s**(-2)*t1 - 16*c2b*mg**2*s**(-2)*t1**2*u1**(-1)
     +     - 16*c2b*mg**2*s**(-2)*u1 - 16*c2b*mg**2*s**(-1)*t1*u1**(-1)
     +     - 16*c2b*mg**2*s**(-1) + 16*c2b*msb2**2*s**(-2)*t1 + 16*c2b*
     +    msb2**2*s**(-2)*t1**2*u1**(-1) + 16*c2b*msb2**2*s**(-2)*u1 + 
     +    16*c2b*msb2**2*s**(-1)*t1*u1**(-1) + 16*c2b*msb2**2*s**(-1)
     +     - 16*c2b*s**(-1)*t1 - 16*c2b*s**(-1)*t1**2*u1**(-1) - 16*c2b
     +    *s**(-1)*u1 - 16*c2b*t1*u1**(-1) - 16*c2b )
      MM_s = MM_s + SCB(4,3)*susy_qi*Nc*Cf**2*Pi**2*alphas**2*prefac
     +  * (  - 64*h_l*h_r*s2b*m1*mg*s**(-1)*t1*u1**(-1) - 64*h_l*h_r*
     +    s2b*m1*mg*s**(-1) - 128*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 16*(h_r2*cb2+h_l2*sb2)*m1**2*mg**2*s**(-1)*u1**(-1) + 16
     +    *(h_r2*cb2+h_l2*sb2)*m1**2*msb2**2*s**(-1)*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*u1**(-1) + 16*(h_r2*cb2+h_l2*sb2)*
     +    mg**2*s**(-2)*t1 - 16*(h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*t1**2
     +    *u1**(-1) + 16*(h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*u1 - 16*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-1) - 16*(h_r2*cb2+h_l2*sb2)*
     +    msb2**2*s**(-2)*t1 + 16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*
     +    t1**2*u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*u1 + 
     +    16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1) - 16*(h_r2*cb2+h_l2*sb2)*
     +    s**(-1)*t1 - 16*(h_r2*cb2+h_l2*sb2)*s**(-1)*t1**2*u1**(-1) - 
     +    16*(h_r2*cb2+h_l2*sb2)*s**(-1)*u1 - 16*(h_r2*cb2+h_l2*sb2)*t1
     +    *u1**(-1) - 16*(h_r2*cb2+h_l2*sb2) + 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * ( 32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qqg*Cf*Pi**2*alphas**2*prefac * ( 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*mg**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*msb2**2*s**(-1)*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*m1**2*u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*
     +    mg**2*s**(-2)*t1 + 16*(h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*t1**2
     +    *u1**(-1) - 16*(h_r2*cb2+h_l2*sb2)*mg**2*s**(-2)*u1 + 16*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-1)*t1*u1**(-1) - 16*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-1) + 16*(h_r2*cb2+h_l2*sb2)*
     +    msb2**2*s**(-2)*t1 - 16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*
     +    t1**2*u1**(-1) + 16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-2)*u1 - 
     +    16*(h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*t1*u1**(-1) + 16*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1) + 16*(h_r2*cb2+h_l2*sb2)*
     +    s**(-1)*t1 + 16*(h_r2*cb2+h_l2*sb2)*s**(-1)*t1**2*u1**(-1) + 
     +    16*(h_r2*cb2+h_l2*sb2)*s**(-1)*u1 + 16*(h_r2*cb2+h_l2*sb2)*t1
     +    *u1**(-1) + 16*(h_r2*cb2+h_l2*sb2) - 32*h_l*h_r*s2b*m1*mg*
     +    s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    32*h_l*h_r*s2b*m1*mg*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*sb*m1*m2**2*s**(-2)*t1*
     +    u2 - 32*h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*ct*sb*m1*
     +    m2**2*s**(-2)*u2**2 - 128*h_l*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*t1 - 64*h_l*ct*sb
     +    *m1*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2 - 
     +    64*h_l*ct*sb*m1*m2**2*t1*u1**(-1) - 96*h_l*ct*sb*m1*m2**2*
     +    u1**(-1)*u2 + 16*h_l*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*m1*
     +    s**(-1)*t1**2 + 32*h_l*ct*sb*m1*s**(-1)*u2**2 - 80*h_l*ct*sb*
     +    m1*s*t1*u1**(-1) + 80*h_l*ct*sb*m1*s*u1**(-1)*u2 - 64*h_l*ct*
     +    sb*m1*t1 - 64*h_l*ct*sb*m1*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*
     +    u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*u2 - 32*h_l*ct*sb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) - 32*h_l*ct*sb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1**3*s**(-2)*t1*u2 + 32
     +    *h_l*ct*sb*m1**3*s**(-2)*t1**2 + 32*h_l*ct*sb*m1**3*s**(-2)*
     +    u2**2 + 128*h_l*ct*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*
     +    ct*sb*m1**3*s**(-1)*t1 + 64*h_l*ct*sb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*ct*sb*m1**3*s**(-1)*u1**(-1)*u2**2 + 64*h_l
     +    *ct*sb*m1**3*s**(-1)*u2 + 64*h_l*ct*sb*m1**3*t1*u1**(-1) + 96
     +    *h_l*ct*sb*m1**3*u1**(-1)*u2 + 16*h_l*ct*sb*m1**5*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*ct*sb*m1**5*s**(-1)*u1**(-1)*u2 + 64*h_r*st
     +    *cb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*st*cb*m1*m2**2*s**(-2)*
     +    t1**2 + 32*h_r*st*cb*m1*m2**2*s**(-2)*u2**2 + 128*h_r*st*cb*
     +    m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*cb*m1*m2**2*
     +    s**(-1)*t1 + 64*h_r*st*cb*m1*m2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*st*cb*m1*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*
     +    m1*m2**2*s**(-1)*u2 + 64*h_r*st*cb*m1*m2**2*t1*u1**(-1) + 96*
     +    h_r*st*cb*m1*m2**2*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*st*cb*m1*m2**4*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *st*cb*m1*s**(-1)*t1**2 - 32*h_r*st*cb*m1*s**(-1)*u2**2 + 80*
     +    h_r*st*cb*m1*s*t1*u1**(-1) - 80*h_r*st*cb*m1*s*u1**(-1)*u2 + 
     +    64*h_r*st*cb*m1*t1 + 64*h_r*st*cb*m1*t1**2*u1**(-1) - 64*h_r*
     +    st*cb*m1*u1**(-1)*u2**2 - 64*h_r*st*cb*m1*u2 + 32*h_r*st*cb*
     +    m1**3*m2**2*s**(-1)*t1*u1**(-1) + 32*h_r*st*cb*m1**3*m2**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*st*cb*m1**3*s**(-2)*t1*u2 - 32*
     +    h_r*st*cb*m1**3*s**(-2)*t1**2 - 32*h_r*st*cb*m1**3*s**(-2)*
     +    u2**2 - 128*h_r*st*cb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    st*cb*m1**3*s**(-1)*t1 - 64*h_r*st*cb*m1**3*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_r*st*cb*m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r
     +    *st*cb*m1**3*s**(-1)*u2 - 64*h_r*st*cb*m1**3*t1*u1**(-1) - 96
     +    *h_r*st*cb*m1**3*u1**(-1)*u2 - 16*h_r*st*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*st*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*sb*m1*m2**2*s**(-2)*t1*u2 + 32*
     +    h_l*ct*sb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*ct*sb*m1*m2**2*
     +    s**(-2)*u2**2 + 128*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     + 64*h_l*ct*sb*m1*m2**2*s**(-1)*t1 + 64*h_l*ct*sb*m1*m2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2 + 64*h_l*ct
     +    *sb*m1*m2**2*t1*u1**(-1) + 96*h_l*ct*sb*m1*m2**2*u1**(-1)*u2
     +     - 16*h_l*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*sb*m1*s**(-1)*t1**2
     +     - 32*h_l*ct*sb*m1*s**(-1)*u2**2 + 80*h_l*ct*sb*m1*s*t1*
     +    u1**(-1) - 80*h_l*ct*sb*m1*s*u1**(-1)*u2 + 64*h_l*ct*sb*m1*t1
     +     + 64*h_l*ct*sb*m1*t1**2*u1**(-1) - 64*h_l*ct*sb*m1*u1**(-1)*
     +    u2**2 - 64*h_l*ct*sb*m1*u2 + 32*h_l*ct*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1) + 32*h_l*ct*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     - 64*h_l*ct*sb*m1**3*s**(-2)*t1*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*ct*sb*m1**3*s**(-2)*t1**2 - 32*
     +    h_l*ct*sb*m1**3*s**(-2)*u2**2 - 128*h_l*ct*sb*m1**3*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*s**(-1)*t1 - 64*h_l*ct*sb
     +    *m1**3*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*ct*sb*m1**3*s**(-1)*u2 - 64*h_l*ct*sb
     +    *m1**3*t1*u1**(-1) - 96*h_l*ct*sb*m1**3*u1**(-1)*u2 - 16*h_l*
     +    ct*sb*m1**5*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 - 64*h_r*st*cb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*st
     +    *cb*m1*m2**2*s**(-2)*t1**2 - 32*h_r*st*cb*m1*m2**2*s**(-2)*
     +    u2**2 - 128*h_r*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_r*st*cb*m1*m2**2*s**(-1)*t1 - 64*h_r*st*cb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) - 64*h_r*st*cb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_r*st*cb*m1*m2**2*s**(-1)*u2 - 64*h_r*st*cb*m1*
     +    m2**2*t1*u1**(-1) - 96*h_r*st*cb*m1*m2**2*u1**(-1)*u2 + 16*
     +    h_r*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_r*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2
     +     - 32*h_r*st*cb*m1*s**(-1)*t1**2 + 32*h_r*st*cb*m1*s**(-1)*
     +    u2**2 - 80*h_r*st*cb*m1*s*t1*u1**(-1) + 80*h_r*st*cb*m1*s*
     +    u1**(-1)*u2 - 64*h_r*st*cb*m1*t1 - 64*h_r*st*cb*m1*t1**2*
     +    u1**(-1) + 64*h_r*st*cb*m1*u1**(-1)*u2**2 + 64*h_r*st*cb*m1*
     +    u2 - 32*h_r*st*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_r*st
     +    *cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*st*cb*m1**3*
     +    s**(-2)*t1*u2 + 32*h_r*st*cb*m1**3*s**(-2)*t1**2 + 32*h_r*st*
     +    cb*m1**3*s**(-2)*u2**2 + 128*h_r*st*cb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_r*st*cb*m1**3*s**(-1)*t1 + 64*h_r*st*cb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*st*cb*m1**3*s**(-1)*u2 + 64*h_r*st*cb
     +    *m1**3*t1*u1**(-1) + 96*h_r*st*cb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    st*cb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) - 16*h_l*ct*
     +    sb*m1*s**(-1)*u1**(-1)*u2 + 16*h_r*st*cb*m1*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb*m1
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) - 
     +    16*h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*sb*m1*u1**(-1) + 32*h_r*st*cb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1*u1**(-1) - 32*h_r*st*cb*m1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*sb*m1*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*st*sb*m1*m2**2*s**(-2)*t1**2 + 32*h_l*st*sb*m1*
     +    m2**2*s**(-2)*u2**2 + 128*h_l*st*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*st*sb*m1*m2**2*s**(-1)*t1 + 64*h_l*st*sb
     +    *m1*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*sb*m1*m2**2*s**(-1)*u2 + 
     +    64*h_l*st*sb*m1*m2**2*t1*u1**(-1) + 96*h_l*st*sb*m1*m2**2*
     +    u1**(-1)*u2 - 16*h_l*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*st*sb*m1*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_l*st*sb*m1*
     +    s**(-1)*t1**2 - 32*h_l*st*sb*m1*s**(-1)*u2**2 + 80*h_l*st*sb*
     +    m1*s*t1*u1**(-1) - 80*h_l*st*sb*m1*s*u1**(-1)*u2 + 64*h_l*st*
     +    sb*m1*t1 + 64*h_l*st*sb*m1*t1**2*u1**(-1) - 64*h_l*st*sb*m1*
     +    u1**(-1)*u2**2 - 64*h_l*st*sb*m1*u2 + 32*h_l*st*sb*m1**3*
     +    m2**2*s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1**3*s**(-2)*t1*u2
     +     - 32*h_l*st*sb*m1**3*s**(-2)*t1**2 - 32*h_l*st*sb*m1**3*
     +    s**(-2)*u2**2 - 128*h_l*st*sb*m1**3*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_l*st*sb*m1**3*s**(-1)*t1 - 64*h_l*st*sb*m1**3*s**(-1)*
     +    t1**2*u1**(-1) - 64*h_l*st*sb*m1**3*s**(-1)*u1**(-1)*u2**2 - 
     +    64*h_l*st*sb*m1**3*s**(-1)*u2 - 64*h_l*st*sb*m1**3*t1*
     +    u1**(-1) - 96*h_l*st*sb*m1**3*u1**(-1)*u2 - 16*h_l*st*sb*
     +    m1**5*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 + 64*h_r*ct*cb*m1*m2**2*s**(-2)*t1*u2 + 32*h_r*ct
     +    *cb*m1*m2**2*s**(-2)*t1**2 + 32*h_r*ct*cb*m1*m2**2*s**(-2)*
     +    u2**2 + 128*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*cb*m1*m2**2*s**(-1)*t1 + 64*h_r*ct*cb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) + 64*h_r*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 + 64*h_r*ct*cb*m1*m2**2*s**(-1)*u2 + 64*h_r*ct*cb*m1*
     +    m2**2*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_r*ct*cb*m1*m2**2*u1**(-1)*u2 - 
     +    16*h_r*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1*
     +    m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m1*s**(-1)*t1**2 - 
     +    32*h_r*ct*cb*m1*s**(-1)*u2**2 + 80*h_r*ct*cb*m1*s*t1*u1**(-1)
     +     - 80*h_r*ct*cb*m1*s*u1**(-1)*u2 + 64*h_r*ct*cb*m1*t1 + 64*
     +    h_r*ct*cb*m1*t1**2*u1**(-1) - 64*h_r*ct*cb*m1*u1**(-1)*u2**2
     +     - 64*h_r*ct*cb*m1*u2 + 32*h_r*ct*cb*m1**3*m2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*ct*cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 - 64*
     +    h_r*ct*cb*m1**3*s**(-2)*t1*u2 - 32*h_r*ct*cb*m1**3*s**(-2)*
     +    t1**2 - 32*h_r*ct*cb*m1**3*s**(-2)*u2**2 - 128*h_r*ct*cb*
     +    m1**3*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*s**(-1)*t1
     +     - 64*h_r*ct*cb*m1**3*s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*cb*
     +    m1**3*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*cb*m1**3*s**(-1)*u2
     +     - 64*h_r*ct*cb*m1**3*t1*u1**(-1) - 96*h_r*ct*cb*m1**3*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*ct*cb*m1**5*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*ct*cb*m1**5*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*sb*m1*m2**2*s**(-2)*t1*u2 - 32
     +    *h_l*st*sb*m1*m2**2*s**(-2)*t1**2 - 32*h_l*st*sb*m1*m2**2*
     +    s**(-2)*u2**2 - 128*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_l*st*sb*m1*m2**2*s**(-1)*t1 - 64*h_l*st*sb*m1*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_l*st*sb*m1*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_l*st*sb*m1*m2**2*s**(-1)*u2 - 64*h_l*st
     +    *sb*m1*m2**2*t1*u1**(-1) - 96*h_l*st*sb*m1*m2**2*u1**(-1)*u2
     +     + 16*h_l*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*st*sb*m1*s**(-1)*t1**2
     +     + 32*h_l*st*sb*m1*s**(-1)*u2**2 - 80*h_l*st*sb*m1*s*t1*
     +    u1**(-1) + 80*h_l*st*sb*m1*s*u1**(-1)*u2 - 64*h_l*st*sb*m1*t1
     +     - 64*h_l*st*sb*m1*t1**2*u1**(-1) + 64*h_l*st*sb*m1*u1**(-1)*
     +    u2**2 + 64*h_l*st*sb*m1*u2 - 32*h_l*st*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1) - 32*h_l*st*sb*m1**3*m2**2*s**(-1)*u1**(-1)*u2
     +     + 64*h_l*st*sb*m1**3*s**(-2)*t1*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*st*sb*m1**3*s**(-2)*t1**2 + 32*h_l*
     +    st*sb*m1**3*s**(-2)*u2**2 + 128*h_l*st*sb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*st*sb*m1**3*s**(-1)*t1 + 64*h_l*st*sb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*sb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*st*sb*m1**3*s**(-1)*u2 + 64*h_l*st*sb
     +    *m1**3*t1*u1**(-1) + 96*h_l*st*sb*m1**3*u1**(-1)*u2 + 16*h_l*
     +    st*sb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*m1**5*s**(-1)*
     +    u1**(-1)*u2 - 64*h_r*ct*cb*m1*m2**2*s**(-2)*t1*u2 - 32*h_r*ct
     +    *cb*m1*m2**2*s**(-2)*t1**2 - 32*h_r*ct*cb*m1*m2**2*s**(-2)*
     +    u2**2 - 128*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_r*ct*cb*m1*m2**2*s**(-1)*t1 - 64*h_r*ct*cb*m1*m2**2*s**(-1)
     +    *t1**2*u1**(-1) - 64*h_r*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_r*ct*cb*m1*m2**2*s**(-1)*u2 - 64*h_r*ct*cb*m1*
     +    m2**2*t1*u1**(-1) - 96*h_r*ct*cb*m1*m2**2*u1**(-1)*u2 + 16*
     +    h_r*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_r*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*u2
     +     - 32*h_r*ct*cb*m1*s**(-1)*t1**2 + 32*h_r*ct*cb*m1*s**(-1)*
     +    u2**2 - 80*h_r*ct*cb*m1*s*t1*u1**(-1) + 80*h_r*ct*cb*m1*s*
     +    u1**(-1)*u2 - 64*h_r*ct*cb*m1*t1 - 64*h_r*ct*cb*m1*t1**2*
     +    u1**(-1) + 64*h_r*ct*cb*m1*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1*
     +    u2 - 32*h_r*ct*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1) - 32*h_r*ct
     +    *cb*m1**3*m2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*ct*cb*m1**3*
     +    s**(-2)*t1*u2 + 32*h_r*ct*cb*m1**3*s**(-2)*t1**2 + 32*h_r*ct*
     +    cb*m1**3*s**(-2)*u2**2 + 128*h_r*ct*cb*m1**3*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_r*ct*cb*m1**3*s**(-1)*t1 + 64*h_r*ct*cb*
     +    m1**3*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*cb*m1**3*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*ct*cb*m1**3*s**(-1)*u2 + 64*h_r*ct*cb
     +    *m1**3*t1*u1**(-1) + 96*h_r*ct*cb*m1**3*u1**(-1)*u2 + 16*h_r*
     +    ct*cb*m1**5*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*m1**5*s**(-1)*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*
     +    m1*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb
     +    *m1*s**(-1)*u1**(-1)*u2 - 16*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1)
     +     - 16*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*u1**(-1) + 32*h_r*ct*cb*m1*u1**(-1)
     +     )
      MM_s = MM_s + SCB(4,3)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*u1**(-1) - 32*h_r*ct*cb*m1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16 - 64*s2t*m1*mg*s**(-1)*t1*u1**(-1) - 64*
     +    s2t*m1*mg*s**(-1)*t1**2*u1**(-2) - 32*s2t*m1*mg*t1*u1**(-2)
     +     - 32*s2t*m1*mg*u1**(-1) + 128*s2t*m1**3*mg*s*u1**(-3) + 128*
     +    s2t*m1**3*mg*t1*u1**(-3) + 64*s2t*m1**3*mg*u1**(-2) + 32*
     +    m1**2*mg**2*s**(-1)*t1*u1**(-1)*u**(-1) + 32*m1**2*mg**2*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) - 32*m1**2*mg**2*s*u1**(-2)*
     +    u**(-1) - 16*m1**2*mg**2*t1*u1**(-2)*u**(-1) - 32*m1**2*
     +    mst1**2*s**(-1)*t1*u1**(-1)*u**(-1) - 32*m1**2*mst1**2*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) + 32*m1**2*mst1**2*s*u1**(-2)*
     +    u**(-1) + 16*m1**2*mst1**2*t1*u1**(-2)*u**(-1) + 32*m1**2*
     +    s**(-1)*t1*u1**(-1) + 32*m1**2*s**(-1)*t1**2*u1**(-2) - 32*
     +    m1**2*s*u1**(-2) - 16*m1**2*t1*u1**(-2) - 64*m1**4*mg**2*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1) - 32*
     +    m1**4*mg**2*u1**(-2)*u**(-1) + 64*m1**4*mst1**2*s*u1**(-3)*
     +    u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 64*m1**4*mst1**2*t1*u1**(-3)*u**(-1) + 32*
     +    m1**4*mst1**2*u1**(-2)*u**(-1) - 64*m1**4*s*u1**(-3) - 64*
     +    m1**4*t1*u1**(-3) - 32*m1**4*u1**(-2) + 16*mg**2*s**(-1)*t1*
     +    u**(-1) + 16*mg**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*mg**2*
     +    s*u1**(-1)*u**(-1) + 16*mg**2*t1*u1**(-1)*u**(-1) + 16*mg**2*
     +    u**(-1) - 16*mst1**2*s**(-1)*t1*u**(-1) - 16*mst1**2*s**(-1)*
     +    t1**2*u1**(-1)*u**(-1) - 16*mst1**2*s*u1**(-1)*u**(-1) - 16*
     +    mst1**2*t1*u1**(-1)*u**(-1) - 16*mst1**2*u**(-1) + 16*s**(-1)
     +    *t1 + 16*s**(-1)*t1**2*u1**(-1) + 16*s*u1**(-1) + 16*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ti*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16*c2t*m1**2*mg**2*t1*u1**(-2)*u**(-1) - 16
     +    *c2t*m1**2*mst1**2*t1*u1**(-2)*u**(-1) + 16*c2t*m1**2*t1*
     +    u1**(-2) - 16*c2t*mg**2*s**(-1)*t1*u**(-1) - 16*c2t*mg**2*
     +    s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*c2t*mg**2*s*u1**(-1)*
     +    u**(-1) - 16*c2t*mg**2*t1*u1**(-1)*u**(-1) - 16*c2t*mg**2*
     +    u**(-1) + 16*c2t*mst1**2*s**(-1)*t1*u**(-1) + 16*c2t*mst1**2*
     +    s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*c2t*mst1**2*s*u1**(-1)*
     +    u**(-1) + 16*c2t*mst1**2*t1*u1**(-1)*u**(-1) + 16*c2t*mst1**2
     +    *u**(-1) - 16*c2t*s**(-1)*t1 - 16*c2t*s**(-1)*t1**2*u1**(-1)
     +     - 16*c2t*s*u1**(-1) - 16*c2t*t1*u1**(-1) - 16*c2t )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * ( 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) + 32*s2t*
     +    m1*mg*s**(-1)*t1**2*u1**(-2) - 16*s2t*m1*mg*s*u1**(-2) + 16*
     +    s2t*m1*mg*t1*u1**(-2) - 64*s2t*m1**3*mg*s*u1**(-3) - 64*s2t*
     +    m1**3*mg*t1*u1**(-3) - 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) - 32*s2t*m1*mg*
     +    s**(-1)*t1**2*u1**(-2) + 16*s2t*m1*mg*s*u1**(-2) - 16*s2t*m1*
     +    mg*t1*u1**(-2) + 64*s2t*m1**3*mg*s*u1**(-3) + 64*s2t*m1**3*mg
     +    *t1*u1**(-3) + 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32 + 16*m1**2*mg**2*s**(-1)*t1**2*
     +    u1**(-2)*u**(-1) + 48*m1**2*mg**2*s*u1**(-3) - 16*m1**2*mg**2
     +    *s*u1**(-2)*u**(-1) + 128*m1**2*mg**2*t1*u1**(-3) - 48*m1**2*
     +    mg**2*t1*u1**(-2)*u**(-1) + 32*m1**2*mg**2*u1**(-2) + 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) - 16*m1**2*mst1**2*s**(-1)*t1**2
     +    *u1**(-2)*u**(-1) - 48*m1**2*mst1**2*s*u1**(-3) + 16*m1**2*
     +    mst1**2*s*u1**(-2)*u**(-1) - 128*m1**2*mst1**2*t1*u1**(-3) + 
     +    48*m1**2*mst1**2*t1*u1**(-2)*u**(-1) - 32*m1**2*mst1**2*
     +    u1**(-2) - 16*m1**2*mst1**2*u1**(-1)*u**(-1) - 16*m1**2*
     +    s**(-1)*t1*u1**(-1) - 32*m1**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 48*m1**2*s*u1**(-2) + 
     +    96*m1**2*t1*u1**(-2) - 32*m1**2*t1*u1**(-1)*u**(-1) + 32*
     +    m1**2*u**(-1) - 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4*mg**2*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1) + 32*
     +    m1**4*mst1**2*s*u1**(-4) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*m1**4*mst1**2*s*u1**(-3)*u**(-1)
     +     + 64*m1**4*mst1**2*t1*u1**(-3)*u**(-1) + 16*m1**4*s**(-1)*
     +    t1**2*u1**(-2)*u**(-1) + 48*m1**4*s*u1**(-3) - 48*m1**4*s*
     +    u1**(-2)*u**(-1) + 128*m1**4*t1*u1**(-3) - 112*m1**4*t1*
     +    u1**(-2)*u**(-1) + 32*m1**4*u1**(-2) + 16*m1**4*u1**(-1)*
     +    u**(-1) - 32*m1**6*s*u1**(-4) - 32*m1**6*s*u1**(-3)*u**(-1)
     +     - 64*m1**6*t1*u1**(-3)*u**(-1) - 16*mg**2*s**(-1)*t1*
     +    u1**(-1) - 32*mg**2*s**(-1)*t1**2*u1**(-2) + 16*mg**2*s**(-1)
     +    *t1**2*u1**(-1)*u**(-1) + 16*mg**2*s*u1**(-1)*u**(-1) + 16*
     +    mg**2*t1*u1**(-1)*u**(-1) + 16*mg**2*u**(-1) + 16*mst1**2*
     +    s**(-1)*t1*u1**(-1) + 32*mst1**2*s**(-1)*t1**2*u1**(-2) - 16*
     +    mst1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*mst1**2*s*
     +    u1**(-1)*u**(-1) - 16*mst1**2*t1*u1**(-1)*u**(-1) - 16*
     +    mst1**2*u**(-1) - 16*s**(-1)*t1 - 32*s**(-1)*t1**2*u1**(-1)
     +     + 16*s**(-1)*t1**2*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*s*u1**(-1) + 16*s*u**(-1) - 32*
     +    t1*u1**(-1) + 16*t1*u**(-1) + 16*u1*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 48*m1**2*mg**2*s*u1**(-3) + 16*m1**2*mg**2*s*
     +    u1**(-2)*u**(-1) - 128*m1**2*mg**2*t1*u1**(-3) + 48*m1**2*
     +    mg**2*t1*u1**(-2)*u**(-1) - 32*m1**2*mg**2*u1**(-2) - 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) + 16*m1**2*mst1**2*s**(-1)*t1**2
     +    *u1**(-2)*u**(-1) + 48*m1**2*mst1**2*s*u1**(-3) - 16*m1**2*
     +    mst1**2*s*u1**(-2)*u**(-1) + 128*m1**2*mst1**2*t1*u1**(-3) - 
     +    48*m1**2*mst1**2*t1*u1**(-2)*u**(-1) + 32*m1**2*mst1**2*
     +    u1**(-2) + 16*m1**2*mst1**2*u1**(-1)*u**(-1) + 48*m1**2*
     +    s**(-1)*t1*u1**(-1) + 32*m1**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 48*m1**2*s*u1**(-2) - 
     +    64*m1**2*t1*u1**(-2) - 32*m1**2*t1*u1**(-1)*u**(-1) - 32*
     +    m1**2*u1**(-1) + 32*m1**2*u**(-1) + 32*m1**4*mg**2*s*u1**(-4)
     +     + 32*m1**4*mg**2*s*u1**(-3)*u**(-1) + 64*m1**4*mg**2*t1*
     +    u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*m1**4*mst1**2*s*u1**(-4) - 32*m1**4*
     +    mst1**2*s*u1**(-3)*u**(-1) - 64*m1**4*mst1**2*t1*u1**(-3)*
     +    u**(-1) + 16*m1**4*s**(-1)*t1**2*u1**(-2)*u**(-1) + 144*m1**4
     +    *s*u1**(-3) - 48*m1**4*s*u1**(-2)*u**(-1) - 64*m1**4*t1*
     +    u1**(-3) - 112*m1**4*t1*u1**(-2)*u**(-1) - 32*m1**4*u1**(-2)
     +     + 16*m1**4*u1**(-1)*u**(-1) + 96*m1**6*s*u1**(-4) - 32*m1**6
     +    *s*u1**(-3)*u**(-1) - 64*m1**6*t1*u1**(-3)*u**(-1) + 16*mg**2
     +    *s**(-1)*t1*u1**(-1) + 32*mg**2*s**(-1)*t1**2*u1**(-2) - 16*
     +    mg**2*s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*mg**2*s*u1**(-1)*
     +    u**(-1) - 16*mg**2*t1*u1**(-1)*u**(-1) - 16*mg**2*u**(-1) - 
     +    16*mst1**2*s**(-1)*t1*u1**(-1) - 32*mst1**2*s**(-1)*t1**2*
     +    u1**(-2) + 16*mst1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*
     +    mst1**2*s*u1**(-1)*u**(-1) + 16*mst1**2*t1*u1**(-1)*u**(-1)
     +     + 16*mst1**2*u**(-1) + 48*s**(-1)*t1 + 32*s**(-1)*t1**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*s**(-1)*t1**2*u**(-1) + 16*s*u**(-1) + 
     +    16*t1*u**(-1) + 16*u1*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*m1**2*mg**2*s**(-1)*t1*u1**(-1)
     +    *u**(-1) - 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*u**(-1) - 64
     +    *m1**2*mg**2*t1*u1**(-3) + 16*m1**2*mg**2*t1*u1**(-2)*u**(-1)
     +     - 16*m1**2*mg**2*u1**(-1)*u**(-1) + 16*m1**2*mst1**2*s**(-1)
     +    *t1*u1**(-1)*u**(-1) + 16*m1**2*mst1**2*s**(-1)*t1**2*
     +    u1**(-2)*u**(-1) + 64*m1**2*mst1**2*t1*u1**(-3) - 16*m1**2*
     +    mst1**2*t1*u1**(-2)*u**(-1) + 16*m1**2*mst1**2*u1**(-1)*
     +    u**(-1) - 16*m1**2*s**(-1)*t1*u**(-1) - 16*m1**2*s**(-1)*
     +    t1**2*u1**(-1)*u**(-1) - 32*m1**2*t1*u1**(-2) + 16*m1**2*t1*
     +    u1**(-1)*u**(-1) - 16*m1**2*u**(-1) + 32*m1**4*mg**2*s*
     +    u1**(-4) + 32*m1**4*mg**2*s*u1**(-3)*u**(-1) + 64*m1**4*mg**2
     +    *t1*u1**(-3)*u**(-1) - 32*m1**4*mst1**2*s*u1**(-4) - 32*m1**4
     +    *mst1**2*s*u1**(-3)*u**(-1) - 64*m1**4*mst1**2*t1*u1**(-3)*
     +    u**(-1) - 16*m1**4*s**(-1)*t1*u1**(-1)*u**(-1) - 16*m1**4*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*m1**4*s*u1**(-2)*u**(-1) - 64*
     +    m1**4*t1*u1**(-3) + 80*m1**4*t1*u1**(-2)*u**(-1) - 16*m1**4*
     +    u1**(-1)*u**(-1) + 32*m1**6*s*u1**(-4) + 32*m1**6*s*u1**(-3)*
     +    u**(-1) + 64*m1**6*t1*u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*mg**2*s**(-1)*t1*u1**(-1)*u**(-1)
     +     + 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*
     +    mg**2*t1*u1**(-3) - 16*m1**2*mg**2*t1*u1**(-2)*u**(-1) + 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) - 16*m1**2*mst1**2*s**(-1)*t1*
     +    u1**(-1)*u**(-1) - 16*m1**2*mst1**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 64*m1**2*mst1**2*t1*u1**(-3) + 16*m1**2*mst1**2*t1*
     +    u1**(-2)*u**(-1) - 16*m1**2*mst1**2*u1**(-1)*u**(-1) - 16*
     +    m1**2*s**(-1)*t1*u**(-1) - 16*m1**2*s**(-1)*t1**2*u1**(-1)*
     +    u**(-1) - 64*m1**2*s*u1**(-2) + 16*m1**2*t1*u1**(-1)*u**(-1)
     +     - 16*m1**2*u**(-1) - 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4*
     +    mg**2*s*u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1)
     +     + 32*m1**4*mst1**2*s*u1**(-4) + 32*m1**4*mst1**2*s*u1**(-3)*
     +    u**(-1) + 64*m1**4*mst1**2*t1*u1**(-3)*u**(-1) - 16*m1**4*
     +    s**(-1)*t1*u1**(-1)*u**(-1) - 16*m1**4*s**(-1)*t1**2*u1**(-2)
     +    *u**(-1) )
      MM_s = MM_s + SCB(5,4)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 160*m1**4*s*u1**(-3) + 32*m1**4*s*
     +    u1**(-2)*u**(-1) + 80*m1**4*t1*u1**(-2)*u**(-1) - 16*m1**4*
     +    u1**(-1)*u**(-1) - 96*m1**6*s*u1**(-4) + 32*m1**6*s*u1**(-3)*
     +    u**(-1) + 64*m1**6*t1*u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,4)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 64*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*ct*cb*m1*
     +    s**(-1)*u1**(-1)*u2 + 96*h_l*ct*cb*m1*s*u1**(-2) + 96*h_l*ct*
     +    cb*m1*t1*u1**(-2) + 96*h_l*ct*cb*m1*u1**(-2)*u2 + 96*h_l*ct*
     +    cb*m1*u1**(-1) + 32*h_l*ct*cb*m1**3*s**(-1)*u1**(-1) + 64*h_l
     +    *ct*cb*m1**3*u1**(-2) + 32*h_r*st*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 64*h_r*st*sb*m1*s**(-1)*t1*u1**(-1) + 64*
     +    h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*st*sb*m1*s*u1**(-2)
     +     + 96*h_r*st*sb*m1*t1*u1**(-2) + 96*h_r*st*sb*m1*u1**(-2)*u2
     +     + 96*h_r*st*sb*m1*u1**(-1) + 32*h_r*st*sb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*st*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 64*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*ct*cb*m1*
     +    s**(-1)*u1**(-1)*u2 - 96*h_l*ct*cb*m1*s*u1**(-2) - 96*h_l*ct*
     +    cb*m1*t1*u1**(-2) - 96*h_l*ct*cb*m1*u1**(-2)*u2 - 96*h_l*ct*
     +    cb*m1*u1**(-1) - 32*h_l*ct*cb*m1**3*s**(-1)*u1**(-1) - 64*h_l
     +    *ct*cb*m1**3*u1**(-2) - 32*h_r*st*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) - 64*h_r*st*sb*m1*s**(-1)*t1*u1**(-1) - 64*
     +    h_r*st*sb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*st*sb*m1*s*u1**(-2)
     +     - 96*h_r*st*sb*m1*t1*u1**(-2) - 96*h_r*st*sb*m1*u1**(-2)*u2
     +     - 96*h_r*st*sb*m1*u1**(-1) - 32*h_r*st*sb*m1**3*s**(-1)*
     +    u1**(-1) - 64*h_r*st*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*
     +    m1*s**(-1) - 64*h_l*ct*cb*m1**3*u1**(-2) + 32*h_r*st*sb*m1*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*st*sb*m1*s**(-1) - 64*h_r*st*sb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*ct*cb
     +    *m1*s**(-1) + 64*h_l*ct*cb*m1**3*u1**(-2) - 32*h_r*st*sb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*sb*m1*s**(-1) + 64*h_r*st*sb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*ct*sb
     +    *m1*s**(-1)*u1**(-1)*u2 - 96*h_l*ct*sb*m1*s*u1**(-2) - 96*h_l
     +    *ct*sb*m1*t1*u1**(-2) - 96*h_l*ct*sb*m1*u1**(-2)*u2 - 96*h_l*
     +    ct*sb*m1*u1**(-1) - 32*h_l*ct*sb*m1**3*s**(-1)*u1**(-1) - 64*
     +    h_l*ct*sb*m1**3*u1**(-2) + 32*h_r*st*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 64*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) + 64*
     +    h_r*st*cb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*st*cb*m1*s*u1**(-2)
     +     + 96*h_r*st*cb*m1*t1*u1**(-2) + 96*h_r*st*cb*m1*u1**(-2)*u2
     +     + 96*h_r*st*cb*m1*u1**(-1) + 32*h_r*st*cb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*st*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 
     +    64*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*ct*sb*m1*s**(-1)
     +    *u1**(-1)*u2 + 96*h_l*ct*sb*m1*s*u1**(-2) + 96*h_l*ct*sb*m1*
     +    t1*u1**(-2) + 96*h_l*ct*sb*m1*u1**(-2)*u2 + 96*h_l*ct*sb*m1*
     +    u1**(-1) + 32*h_l*ct*sb*m1**3*s**(-1)*u1**(-1) + 64*h_l*ct*sb
     +    *m1**3*u1**(-2) - 32*h_r*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_r*st*cb*m1*s**(-1)*t1*u1**(-1) - 64*h_r*st*cb
     +    *m1*s**(-1)*u1**(-1)*u2 - 96*h_r*st*cb*m1*s*u1**(-2) - 96*h_r
     +    *st*cb*m1*t1*u1**(-2) - 96*h_r*st*cb*m1*u1**(-2)*u2 - 96*h_r*
     +    st*cb*m1*u1**(-1) - 32*h_r*st*cb*m1**3*s**(-1)*u1**(-1) - 64*
     +    h_r*st*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*ct*
     +    sb*m1*s**(-1) + 64*h_l*ct*sb*m1**3*u1**(-2) + 32*h_r*st*cb*m1
     +    *s**(-1)*t1*u1**(-1) + 32*h_r*st*cb*m1*s**(-1) - 64*h_r*st*cb
     +    *m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,4)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1
     +    *s**(-1) - 64*h_l*ct*sb*m1**3*u1**(-2) - 32*h_r*st*cb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*cb*m1*s**(-1) + 64*h_r*st*cb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 16 + 64*s2t*m1*mg*s**(-1)*t1*u1**(-1) + 64*
     +    s2t*m1*mg*s**(-1)*t1**2*u1**(-2) + 32*s2t*m1*mg*t1*u1**(-2)
     +     + 32*s2t*m1*mg*u1**(-1) - 128*s2t*m1**3*mg*s*u1**(-3) - 128*
     +    s2t*m1**3*mg*t1*u1**(-3) - 64*s2t*m1**3*mg*u1**(-2) + 32*
     +    m1**2*mg**2*s**(-1)*t1*u1**(-1)*u**(-1) + 32*m1**2*mg**2*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) - 32*m1**2*mg**2*s*u1**(-2)*
     +    u**(-1) - 16*m1**2*mg**2*t1*u1**(-2)*u**(-1) - 32*m1**2*
     +    mst2**2*s**(-1)*t1*u1**(-1)*u**(-1) - 32*m1**2*mst2**2*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) + 32*m1**2*mst2**2*s*u1**(-2)*
     +    u**(-1) + 16*m1**2*mst2**2*t1*u1**(-2)*u**(-1) + 32*m1**2*
     +    s**(-1)*t1*u1**(-1) + 32*m1**2*s**(-1)*t1**2*u1**(-2) - 32*
     +    m1**2*s*u1**(-2) - 16*m1**2*t1*u1**(-2) - 64*m1**4*mg**2*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1) - 32*
     +    m1**4*mg**2*u1**(-2)*u**(-1) + 64*m1**4*mst2**2*s*u1**(-3)*
     +    u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ti*(h_r2+h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * ( 64*m1**4*mst2**2*t1*u1**(-3)*u**(-1) + 32*
     +    m1**4*mst2**2*u1**(-2)*u**(-1) - 64*m1**4*s*u1**(-3) - 64*
     +    m1**4*t1*u1**(-3) - 32*m1**4*u1**(-2) + 16*mg**2*s**(-1)*t1*
     +    u**(-1) + 16*mg**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*mg**2*
     +    s*u1**(-1)*u**(-1) + 16*mg**2*t1*u1**(-1)*u**(-1) + 16*mg**2*
     +    u**(-1) - 16*mst2**2*s**(-1)*t1*u**(-1) - 16*mst2**2*s**(-1)*
     +    t1**2*u1**(-1)*u**(-1) - 16*mst2**2*s*u1**(-1)*u**(-1) - 16*
     +    mst2**2*t1*u1**(-1)*u**(-1) - 16*mst2**2*u**(-1) + 16*s**(-1)
     +    *t1 + 16*s**(-1)*t1**2*u1**(-1) + 16*s*u1**(-1) + 16*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ti*(h_r2-h_l2)*Nc*Cf**2*Pi**2*
     + alphas**2*prefac * (  - 16*c2t*m1**2*mg**2*t1*u1**(-2)*u**(-1)
     +     + 16*c2t*m1**2*mst2**2*t1*u1**(-2)*u**(-1) - 16*c2t*m1**2*t1
     +    *u1**(-2) + 16*c2t*mg**2*s**(-1)*t1*u**(-1) + 16*c2t*mg**2*
     +    s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*c2t*mg**2*s*u1**(-1)*
     +    u**(-1) + 16*c2t*mg**2*t1*u1**(-1)*u**(-1) + 16*c2t*mg**2*
     +    u**(-1) - 16*c2t*mst2**2*s**(-1)*t1*u**(-1) - 16*c2t*mst2**2*
     +    s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*c2t*mst2**2*s*u1**(-1)*
     +    u**(-1) - 16*c2t*mst2**2*t1*u1**(-1)*u**(-1) - 16*c2t*mst2**2
     +    *u**(-1) + 16*c2t*s**(-1)*t1 + 16*c2t*s**(-1)*t1**2*u1**(-1)
     +     + 16*c2t*s*u1**(-1) + 16*c2t*t1*u1**(-1) + 16*c2t )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) - 32*
     +    s2t*m1*mg*s**(-1)*t1**2*u1**(-2) + 16*s2t*m1*mg*s*u1**(-2) - 
     +    16*s2t*m1*mg*t1*u1**(-2) + 64*s2t*m1**3*mg*s*u1**(-3) + 64*
     +    s2t*m1**3*mg*t1*u1**(-3) + 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * ( 32*s2t*m1*mg*s**(-1)*t1*u1**(-1) + 32*s2t*m1*mg*
     +    s**(-1)*t1**2*u1**(-2) - 16*s2t*m1*mg*s*u1**(-2) + 16*s2t*m1*
     +    mg*t1*u1**(-2) - 64*s2t*m1**3*mg*s*u1**(-3) - 64*s2t*m1**3*mg
     +    *t1*u1**(-3) - 32*s2t*m1**3*mg*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*m1**2*mg**2*s**(-1)*t1*u1**(-1)
     +    *u**(-1) - 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*u**(-1) - 64
     +    *m1**2*mg**2*t1*u1**(-3) + 16*m1**2*mg**2*t1*u1**(-2)*u**(-1)
     +     - 16*m1**2*mg**2*u1**(-1)*u**(-1) + 16*m1**2*mst2**2*s**(-1)
     +    *t1*u1**(-1)*u**(-1) + 16*m1**2*mst2**2*s**(-1)*t1**2*
     +    u1**(-2)*u**(-1) + 64*m1**2*mst2**2*t1*u1**(-3) - 16*m1**2*
     +    mst2**2*t1*u1**(-2)*u**(-1) + 16*m1**2*mst2**2*u1**(-1)*
     +    u**(-1) - 16*m1**2*s**(-1)*t1*u**(-1) - 16*m1**2*s**(-1)*
     +    t1**2*u1**(-1)*u**(-1) - 32*m1**2*t1*u1**(-2) + 16*m1**2*t1*
     +    u1**(-1)*u**(-1) - 16*m1**2*u**(-1) + 32*m1**4*mg**2*s*
     +    u1**(-4) + 32*m1**4*mg**2*s*u1**(-3)*u**(-1) + 64*m1**4*mg**2
     +    *t1*u1**(-3)*u**(-1) - 32*m1**4*mst2**2*s*u1**(-4) - 32*m1**4
     +    *mst2**2*s*u1**(-3)*u**(-1) - 64*m1**4*mst2**2*t1*u1**(-3)*
     +    u**(-1) - 16*m1**4*s**(-1)*t1*u1**(-1)*u**(-1) - 16*m1**4*
     +    s**(-1)*t1**2*u1**(-2)*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*m1**4*s*u1**(-2)*u**(-1) - 64*
     +    m1**4*t1*u1**(-3) + 80*m1**4*t1*u1**(-2)*u**(-1) - 16*m1**4*
     +    u1**(-1)*u**(-1) + 32*m1**6*s*u1**(-4) + 32*m1**6*s*u1**(-3)*
     +    u**(-1) + 64*m1**6*t1*u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*m1**2*mg**2*s**(-1)*t1*u1**(-1)*u**(-1)
     +     + 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*u**(-1) + 64*m1**2*
     +    mg**2*t1*u1**(-3) - 16*m1**2*mg**2*t1*u1**(-2)*u**(-1) + 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) - 16*m1**2*mst2**2*s**(-1)*t1*
     +    u1**(-1)*u**(-1) - 16*m1**2*mst2**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 64*m1**2*mst2**2*t1*u1**(-3) + 16*m1**2*mst2**2*t1*
     +    u1**(-2)*u**(-1) - 16*m1**2*mst2**2*u1**(-1)*u**(-1) - 16*
     +    m1**2*s**(-1)*t1*u**(-1) - 16*m1**2*s**(-1)*t1**2*u1**(-1)*
     +    u**(-1) - 64*m1**2*s*u1**(-2) + 16*m1**2*t1*u1**(-1)*u**(-1)
     +     - 16*m1**2*u**(-1) - 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4*
     +    mg**2*s*u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1)
     +     + 32*m1**4*mst2**2*s*u1**(-4) + 32*m1**4*mst2**2*s*u1**(-3)*
     +    u**(-1) + 64*m1**4*mst2**2*t1*u1**(-3)*u**(-1) - 16*m1**4*
     +    s**(-1)*t1*u1**(-1)*u**(-1) - 16*m1**4*s**(-1)*t1**2*u1**(-2)
     +    *u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 160*m1**4*s*u1**(-3) + 32*m1**4*s*
     +    u1**(-2)*u**(-1) + 80*m1**4*t1*u1**(-2)*u**(-1) - 16*m1**4*
     +    u1**(-1)*u**(-1) - 96*m1**6*s*u1**(-4) + 32*m1**6*s*u1**(-3)*
     +    u**(-1) + 64*m1**6*t1*u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32 + 16*m1**2*mg**2*s**(-1)*t1**2*
     +    u1**(-2)*u**(-1) + 48*m1**2*mg**2*s*u1**(-3) - 16*m1**2*mg**2
     +    *s*u1**(-2)*u**(-1) + 128*m1**2*mg**2*t1*u1**(-3) - 48*m1**2*
     +    mg**2*t1*u1**(-2)*u**(-1) + 32*m1**2*mg**2*u1**(-2) + 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) - 16*m1**2*mst2**2*s**(-1)*t1**2
     +    *u1**(-2)*u**(-1) - 48*m1**2*mst2**2*s*u1**(-3) + 16*m1**2*
     +    mst2**2*s*u1**(-2)*u**(-1) - 128*m1**2*mst2**2*t1*u1**(-3) + 
     +    48*m1**2*mst2**2*t1*u1**(-2)*u**(-1) - 32*m1**2*mst2**2*
     +    u1**(-2) - 16*m1**2*mst2**2*u1**(-1)*u**(-1) - 16*m1**2*
     +    s**(-1)*t1*u1**(-1) - 32*m1**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 48*m1**2*s*u1**(-2) + 
     +    96*m1**2*t1*u1**(-2) - 32*m1**2*t1*u1**(-1)*u**(-1) + 32*
     +    m1**2*u**(-1) - 32*m1**4*mg**2*s*u1**(-4) - 32*m1**4*mg**2*s*
     +    u1**(-3)*u**(-1) - 64*m1**4*mg**2*t1*u1**(-3)*u**(-1) + 32*
     +    m1**4*mst2**2*s*u1**(-4) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*m1**4*mst2**2*s*u1**(-3)*u**(-1)
     +     + 64*m1**4*mst2**2*t1*u1**(-3)*u**(-1) + 16*m1**4*s**(-1)*
     +    t1**2*u1**(-2)*u**(-1) + 48*m1**4*s*u1**(-3) - 48*m1**4*s*
     +    u1**(-2)*u**(-1) + 128*m1**4*t1*u1**(-3) - 112*m1**4*t1*
     +    u1**(-2)*u**(-1) + 32*m1**4*u1**(-2) + 16*m1**4*u1**(-1)*
     +    u**(-1) - 32*m1**6*s*u1**(-4) - 32*m1**6*s*u1**(-3)*u**(-1)
     +     - 64*m1**6*t1*u1**(-3)*u**(-1) - 16*mg**2*s**(-1)*t1*
     +    u1**(-1) - 32*mg**2*s**(-1)*t1**2*u1**(-2) + 16*mg**2*s**(-1)
     +    *t1**2*u1**(-1)*u**(-1) + 16*mg**2*s*u1**(-1)*u**(-1) + 16*
     +    mg**2*t1*u1**(-1)*u**(-1) + 16*mg**2*u**(-1) + 16*mst2**2*
     +    s**(-1)*t1*u1**(-1) + 32*mst2**2*s**(-1)*t1**2*u1**(-2) - 16*
     +    mst2**2*s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*mst2**2*s*
     +    u1**(-1)*u**(-1) - 16*mst2**2*t1*u1**(-1)*u**(-1) - 16*
     +    mst2**2*u**(-1) - 16*s**(-1)*t1 - 32*s**(-1)*t1**2*u1**(-1)
     +     + 16*s**(-1)*t1**2*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*s*u1**(-1) + 16*s*u**(-1) - 32*
     +    t1*u1**(-1) + 16*t1*u**(-1) + 16*u1*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*m1**2*mg**2*s**(-1)*t1**2*u1**(-2)*
     +    u**(-1) - 48*m1**2*mg**2*s*u1**(-3) + 16*m1**2*mg**2*s*
     +    u1**(-2)*u**(-1) - 128*m1**2*mg**2*t1*u1**(-3) + 48*m1**2*
     +    mg**2*t1*u1**(-2)*u**(-1) - 32*m1**2*mg**2*u1**(-2) - 16*
     +    m1**2*mg**2*u1**(-1)*u**(-1) + 16*m1**2*mst2**2*s**(-1)*t1**2
     +    *u1**(-2)*u**(-1) + 48*m1**2*mst2**2*s*u1**(-3) - 16*m1**2*
     +    mst2**2*s*u1**(-2)*u**(-1) + 128*m1**2*mst2**2*t1*u1**(-3) - 
     +    48*m1**2*mst2**2*t1*u1**(-2)*u**(-1) + 32*m1**2*mst2**2*
     +    u1**(-2) + 16*m1**2*mst2**2*u1**(-1)*u**(-1) + 48*m1**2*
     +    s**(-1)*t1*u1**(-1) + 32*m1**2*s**(-1)*t1**2*u1**(-2) + 32*
     +    m1**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 48*m1**2*s*u1**(-2) - 
     +    64*m1**2*t1*u1**(-2) - 32*m1**2*t1*u1**(-1)*u**(-1) - 32*
     +    m1**2*u1**(-1) + 32*m1**2*u**(-1) + 32*m1**4*mg**2*s*u1**(-4)
     +     + 32*m1**4*mg**2*s*u1**(-3)*u**(-1) + 64*m1**4*mg**2*t1*
     +    u1**(-3)*u**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*m1**4*mst2**2*s*u1**(-4) - 32*m1**4*
     +    mst2**2*s*u1**(-3)*u**(-1) - 64*m1**4*mst2**2*t1*u1**(-3)*
     +    u**(-1) + 16*m1**4*s**(-1)*t1**2*u1**(-2)*u**(-1) + 144*m1**4
     +    *s*u1**(-3) - 48*m1**4*s*u1**(-2)*u**(-1) - 64*m1**4*t1*
     +    u1**(-3) - 112*m1**4*t1*u1**(-2)*u**(-1) - 32*m1**4*u1**(-2)
     +     + 16*m1**4*u1**(-1)*u**(-1) + 96*m1**6*s*u1**(-4) - 32*m1**6
     +    *s*u1**(-3)*u**(-1) - 64*m1**6*t1*u1**(-3)*u**(-1) + 16*mg**2
     +    *s**(-1)*t1*u1**(-1) + 32*mg**2*s**(-1)*t1**2*u1**(-2) - 16*
     +    mg**2*s**(-1)*t1**2*u1**(-1)*u**(-1) - 16*mg**2*s*u1**(-1)*
     +    u**(-1) - 16*mg**2*t1*u1**(-1)*u**(-1) - 16*mg**2*u**(-1) - 
     +    16*mst2**2*s**(-1)*t1*u1**(-1) - 32*mst2**2*s**(-1)*t1**2*
     +    u1**(-2) + 16*mst2**2*s**(-1)*t1**2*u1**(-1)*u**(-1) + 16*
     +    mst2**2*s*u1**(-1)*u**(-1) + 16*mst2**2*t1*u1**(-1)*u**(-1)
     +     + 16*mst2**2*u**(-1) + 48*s**(-1)*t1 + 32*s**(-1)*t1**2*
     +    u1**(-1) )
      MM_s = MM_s + SCB(5,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*s**(-1)*t1**2*u**(-1) + 16*s*u**(-1) + 
     +    16*t1*u**(-1) + 16*u1*u**(-1) )
      MM_s = MM_s + SCB(5,5)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*st*cb
     +    *m1*s**(-1)*u1**(-1)*u2 - 96*h_l*st*cb*m1*s*u1**(-2) - 96*h_l
     +    *st*cb*m1*t1*u1**(-2) - 96*h_l*st*cb*m1*u1**(-2)*u2 - 96*h_l*
     +    st*cb*m1*u1**(-1) - 32*h_l*st*cb*m1**3*s**(-1)*u1**(-1) - 64*
     +    h_l*st*cb*m1**3*u1**(-2) + 32*h_r*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 64*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) + 64*
     +    h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*sb*m1*s*u1**(-2)
     +     + 96*h_r*ct*sb*m1*t1*u1**(-2) + 96*h_r*ct*sb*m1*u1**(-2)*u2
     +     + 96*h_r*ct*sb*m1*u1**(-1) + 32*h_r*ct*sb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*ct*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 
     +    64*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*st*cb*m1*s**(-1)
     +    *u1**(-1)*u2 + 96*h_l*st*cb*m1*s*u1**(-2) + 96*h_l*st*cb*m1*
     +    t1*u1**(-2) + 96*h_l*st*cb*m1*u1**(-2)*u2 + 96*h_l*st*cb*m1*
     +    u1**(-1) + 32*h_l*st*cb*m1**3*s**(-1)*u1**(-1) + 64*h_l*st*cb
     +    *m1**3*u1**(-2) - 32*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1) - 64*h_r*ct*sb
     +    *m1*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*sb*m1*s*u1**(-2) - 96*h_r
     +    *ct*sb*m1*t1*u1**(-2) - 96*h_r*ct*sb*m1*u1**(-2)*u2 - 96*h_r*
     +    ct*sb*m1*u1**(-1) - 32*h_r*ct*sb*m1**3*s**(-1)*u1**(-1) - 64*
     +    h_r*ct*sb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*st*
     +    cb*m1*s**(-1) + 64*h_l*st*cb*m1**3*u1**(-2) + 32*h_r*ct*sb*m1
     +    *s**(-1)*t1*u1**(-1) + 32*h_r*ct*sb*m1*s**(-1) - 64*h_r*ct*sb
     +    *m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1
     +    *s**(-1) - 64*h_l*st*cb*m1**3*u1**(-2) - 32*h_r*ct*sb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb*m1*s**(-1) + 64*h_r*ct*sb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 64*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) + 64*h_l*st*sb*m1*
     +    s**(-1)*u1**(-1)*u2 + 96*h_l*st*sb*m1*s*u1**(-2) + 96*h_l*st*
     +    sb*m1*t1*u1**(-2) + 96*h_l*st*sb*m1*u1**(-2)*u2 + 96*h_l*st*
     +    sb*m1*u1**(-1) + 32*h_l*st*sb*m1**3*s**(-1)*u1**(-1) + 64*h_l
     +    *st*sb*m1**3*u1**(-2) + 32*h_r*ct*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) + 64*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1) + 64*
     +    h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 + 96*h_r*ct*cb*m1*s*u1**(-2)
     +     + 96*h_r*ct*cb*m1*t1*u1**(-2) + 96*h_r*ct*cb*m1*u1**(-2)*u2
     +     + 96*h_r*ct*cb*m1*u1**(-1) + 32*h_r*ct*cb*m1**3*s**(-1)*
     +    u1**(-1) + 64*h_r*ct*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 64*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 64*h_l*st*sb*m1*
     +    s**(-1)*u1**(-1)*u2 - 96*h_l*st*sb*m1*s*u1**(-2) - 96*h_l*st*
     +    sb*m1*t1*u1**(-2) - 96*h_l*st*sb*m1*u1**(-2)*u2 - 96*h_l*st*
     +    sb*m1*u1**(-1) - 32*h_l*st*sb*m1**3*s**(-1)*u1**(-1) - 64*h_l
     +    *st*sb*m1**3*u1**(-2) - 32*h_r*ct*cb*m1*m2**2*s**(-1)*t1*
     +    u1**(-1)*u2**(-1) - 64*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1) - 64*
     +    h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2 - 96*h_r*ct*cb*m1*s*u1**(-2)
     +     - 96*h_r*ct*cb*m1*t1*u1**(-2) - 96*h_r*ct*cb*m1*u1**(-2)*u2
     +     - 96*h_r*ct*cb*m1*u1**(-1) - 32*h_r*ct*cb*m1**3*s**(-1)*
     +    u1**(-1) - 64*h_r*ct*cb*m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*
     +    m1*s**(-1) - 64*h_l*st*sb*m1**3*u1**(-2) + 32*h_r*ct*cb*m1*
     +    s**(-1)*t1*u1**(-1) + 32*h_r*ct*cb*m1*s**(-1) - 64*h_r*ct*cb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCB(5,5)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*s**(-1)*t1*u1**(-1) - 32*h_l*st*sb
     +    *m1*s**(-1) + 64*h_l*st*sb*m1**3*u1**(-2) - 32*h_r*ct*cb*m1*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb*m1*s**(-1) + 64*h_r*ct*cb*
     +    m1**3*u1**(-2) )
      MM_s = MM_s + SCBP(1,1)*susy_q*(h_r2+h_l2)**(-1)*(h_r2-h_l2)*Cf*
     + Pi*alphas*prefac * (  - 4*c2b*mg**2*born + 4*c2b*msb1**2*born )
      MM_s = MM_s + SCBP(1,1)*susy_q*Cf*Pi*alphas*prefac * (  - 2*mg**2
     +    *born + 2*msb1**2*born )
      MM_s = MM_s + SCBP(1,2)*susy_q*(h_r2+h_l2)**(-1)*(h_r2-h_l2)*Cf*
     + Pi*alphas*prefac * ( 4*c2b*mg**2*born - 4*c2b*msb2**2*born )
      MM_s = MM_s + SCBP(1,2)*susy_q*Cf*Pi*alphas*prefac * (  - 2*mg**2
     +    *born + 2*msb2**2*born )
      MM_s = MM_s + SCBP(2,2)*susy_t*Cf*Pi*alphas*prefac * ( 8*s2t*m1*
     +    mg*born - 4*m1**2*born - 4*mg**2*born + 4*mst1**2*born )
      MM_s = MM_s + SCBP(2,3)*susy_t*Cf*Pi*alphas*prefac * (  - 8*s2t*
     +    m1*mg*born - 4*m1**2*born - 4*mg**2*born + 4*mst2**2*born )
      MM_s = MM_s + SCC(1,2)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    32*(h_r2*sb2+h_l2*cb2)*m1**2*msb1**2*u1**(-1) + 32*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*t1 + 32*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2*s**(-1)*u1 + 32*
     +    (h_r2*sb2+h_l2*cb2)*msb1**2 )
      MM_s = MM_s + SCC(1,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*cb*mg*s*u1**(-1) + 32*h_l*st*cb*mg*t1*
     +    u1**(-1) + 16*h_l*st*cb*mg + 32*h_l*ct*cb*m1*mg**2*u1**(-1)
     +     - 32*h_l*ct*cb*m1*msb1**2*u1**(-1) + 32*h_r*st*sb*m1*mg**2*
     +    u1**(-1) - 32*h_r*st*sb*m1*msb1**2*u1**(-1) + 16*h_r*ct*sb*mg
     +    *s*u1**(-1) + 32*h_r*ct*sb*mg*t1*u1**(-1) + 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(1,2)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*mg**2*u1**(-1) + 32*h_l*st*cb*m1*
     +    msb1**2*u1**(-1) + 16*h_l*ct*cb*mg*s*u1**(-1) + 32*h_l*ct*cb*
     +    mg*t1*u1**(-1) + 16*h_l*ct*cb*mg - 16*h_r*st*sb*mg*s*u1**(-1)
     +     - 32*h_r*st*sb*mg*t1*u1**(-1) - 16*h_r*st*sb*mg + 32*h_r*ct*
     +    sb*m1*mg**2*u1**(-1) - 32*h_r*ct*sb*m1*msb1**2*u1**(-1) )
      MM_s = MM_s + SCC(1,3)*susy_qqg*Cf*Pi**2*alphas**2*prefac * (  - 
     +    32*(h_r2*cb2+h_l2*sb2)*m1**2*msb2**2*u1**(-1) + 32*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*t1 + 32*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2*s**(-1)*u1 + 32*
     +    (h_r2*cb2+h_l2*sb2)*msb2**2 )
      MM_s = MM_s + SCC(1,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*sb*mg*s*u1**(-1) - 32*h_l*st*sb*mg*t1*
     +    u1**(-1) - 16*h_l*st*sb*mg - 32*h_l*ct*sb*m1*mg**2*u1**(-1)
     +     + 32*h_l*ct*sb*m1*msb2**2*u1**(-1) + 32*h_r*st*cb*m1*mg**2*
     +    u1**(-1) - 32*h_r*st*cb*m1*msb2**2*u1**(-1) + 16*h_r*ct*cb*mg
     +    *s*u1**(-1) + 32*h_r*ct*cb*mg*t1*u1**(-1) + 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(1,3)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*mg**2*u1**(-1) - 32*h_l*st*sb*m1*
     +    msb2**2*u1**(-1) - 16*h_l*ct*sb*mg*s*u1**(-1) - 32*h_l*ct*sb*
     +    mg*t1*u1**(-1) - 16*h_l*ct*sb*mg - 16*h_r*st*cb*mg*s*u1**(-1)
     +     - 32*h_r*st*cb*mg*t1*u1**(-1) - 16*h_r*st*cb*mg + 32*h_r*ct*
     +    cb*m1*mg**2*u1**(-1) - 32*h_r*ct*cb*m1*msb2**2*u1**(-1) )
      MM_s = MM_s + SCC(1,4)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 32*(h_r2*sb2+h_l2*cb2)*m1**2*mg**2*u1**(-1) + 32*
     +    (h_r2*sb2+h_l2*cb2)*mg**2*s**(-1)*t1 + 32*(h_r2*sb2+h_l2*cb2)
     +    *mg**2*s**(-1)*u1 + 32*(h_r2*sb2+h_l2*cb2)*mg**2 + 32*h_l*h_r
     +    *s2b*m1*mg*s*u1**(-1) + 32*h_l*h_r*s2b*m1*mg )
      MM_s = MM_s + SCC(1,4)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*st*cb*mg*s*u1**(-1) + 32*h_l*st*cb*mg*t1*
     +    u1**(-1) + 16*h_l*st*cb*mg + 32*h_l*ct*cb*m1*mg**2*u1**(-1)
     +     - 32*h_l*ct*cb*m1*msb1**2*u1**(-1) + 16*h_l*ct*cb*m1*s*
     +    u1**(-1) + 32*h_r*st*sb*m1*mg**2*u1**(-1) - 32*h_r*st*sb*m1*
     +    msb1**2*u1**(-1) + 16*h_r*st*sb*m1*s*u1**(-1) + 16*h_r*ct*sb*
     +    mg*s*u1**(-1) + 32*h_r*ct*sb*mg*t1*u1**(-1) + 16*h_r*ct*sb*mg
     +     )
      MM_s = MM_s + SCC(1,4)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*mg**2*u1**(-1) + 32*h_l*st*cb*m1*
     +    msb1**2*u1**(-1) - 16*h_l*st*cb*m1*s*u1**(-1) + 16*h_l*ct*cb*
     +    mg*s*u1**(-1) + 32*h_l*ct*cb*mg*t1*u1**(-1) + 16*h_l*ct*cb*mg
     +     - 16*h_r*st*sb*mg*s*u1**(-1) - 32*h_r*st*sb*mg*t1*u1**(-1)
     +     - 16*h_r*st*sb*mg + 32*h_r*ct*sb*m1*mg**2*u1**(-1) - 32*h_r*
     +    ct*sb*m1*msb1**2*u1**(-1) + 16*h_r*ct*sb*m1*s*u1**(-1) )
      MM_s = MM_s + SCC(1,5)*susy_qqg*Nc**2*Cf*Pi**2*alphas**2*prefac
     +  * (  - 32*(h_r2*cb2+h_l2*sb2)*m1**2*mg**2*u1**(-1) + 32*
     +    (h_r2*cb2+h_l2*sb2)*mg**2*s**(-1)*t1 + 32*(h_r2*cb2+h_l2*sb2)
     +    *mg**2*s**(-1)*u1 + 32*(h_r2*cb2+h_l2*sb2)*mg**2 - 32*h_l*h_r
     +    *s2b*m1*mg*s*u1**(-1) - 32*h_l*h_r*s2b*m1*mg )
      MM_s = MM_s + SCC(1,5)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*st*sb*mg*s*u1**(-1) - 32*h_l*st*sb*mg*t1*
     +    u1**(-1) - 16*h_l*st*sb*mg - 32*h_l*ct*sb*m1*mg**2*u1**(-1)
     +     + 32*h_l*ct*sb*m1*msb2**2*u1**(-1) - 16*h_l*ct*sb*m1*s*
     +    u1**(-1) + 32*h_r*st*cb*m1*mg**2*u1**(-1) - 32*h_r*st*cb*m1*
     +    msb2**2*u1**(-1) + 16*h_r*st*cb*m1*s*u1**(-1) + 16*h_r*ct*cb*
     +    mg*s*u1**(-1) + 32*h_r*ct*cb*mg*t1*u1**(-1) + 16*h_r*ct*cb*mg
     +     )
      MM_s = MM_s + SCC(1,5)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*mg**2*u1**(-1) - 32*h_l*st*sb*m1*
     +    msb2**2*u1**(-1) + 16*h_l*st*sb*m1*s*u1**(-1) - 16*h_l*ct*sb*
     +    mg*s*u1**(-1) - 32*h_l*ct*sb*mg*t1*u1**(-1) - 16*h_l*ct*sb*mg
     +     - 16*h_r*st*cb*mg*s*u1**(-1) - 32*h_r*st*cb*mg*t1*u1**(-1)
     +     - 16*h_r*st*cb*mg + 32*h_r*ct*cb*m1*mg**2*u1**(-1) - 32*h_r*
     +    ct*cb*m1*msb2**2*u1**(-1) + 16*h_r*ct*cb*m1*s*u1**(-1) )
      MM_s = MM_s + SCC(2,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*mg*s**(-1)*t1 - 32*h_l*st*cb*mg*t1*
     +    u1**(-1) - 16*h_l*ct*cb*m1*s**(-1)*t1 - 16*h_l*ct*cb*m1*t1*
     +    u1**(-1) - 16*h_r*st*sb*m1*s**(-1)*t1 - 16*h_r*st*sb*m1*t1*
     +    u1**(-1) + 32*h_r*ct*sb*mg*s**(-1)*t1 - 32*h_r*ct*sb*mg*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(2,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*mg*s**(-1)*t1 + 32*h_l*st*sb*mg*t1*
     +    u1**(-1) + 16*h_l*ct*sb*m1*s**(-1)*t1 + 16*h_l*ct*sb*m1*t1*
     +    u1**(-1) - 16*h_r*st*cb*m1*s**(-1)*t1 - 16*h_r*st*cb*m1*t1*
     +    u1**(-1) + 32*h_r*ct*cb*mg*s**(-1)*t1 - 32*h_r*ct*cb*mg*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(2,4)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*cb*m1*s**(-1)*t1 + 16*h_l*st*cb*m1*t1*
     +    u1**(-1) + 32*h_l*ct*cb*mg*s**(-1)*t1 - 32*h_l*ct*cb*mg*t1*
     +    u1**(-1) - 32*h_r*st*sb*mg*s**(-1)*t1 + 32*h_r*st*sb*mg*t1*
     +    u1**(-1) - 16*h_r*ct*sb*m1*s**(-1)*t1 - 16*h_r*ct*sb*m1*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(2,5)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*sb*m1*s**(-1)*t1 - 16*h_l*st*sb*m1*t1*
     +    u1**(-1) - 32*h_l*ct*sb*mg*s**(-1)*t1 + 32*h_l*ct*sb*mg*t1*
     +    u1**(-1) - 32*h_r*st*cb*mg*s**(-1)*t1 + 32*h_r*st*cb*mg*t1*
     +    u1**(-1) - 16*h_r*ct*cb*m1*s**(-1)*t1 - 16*h_r*ct*cb*m1*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(3,3)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * ( 32*s2t*m1*mg*s**(-1)*t1 + 32*s2t*m1*mg*s**(-1)*t1**2*
     +    u1**(-1) - 16*s2t*m1*mg*s*u1**(-1) + 16*s2t*m1*mg*t1*u1**(-1)
     +     - 64*s2t*m1**3*mg*s*u1**(-2) - 64*s2t*m1**3*mg*t1*u1**(-2)
     +     - 32*s2t*m1**3*mg*u1**(-1) )
      MM_s = MM_s + SCC(3,3)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*m1**2*mg**2*t1*u1**(-2) + 32*m1**2*mg**2
     +    *u1**(-1) - 32*m1**2*mst1**2*t1*u1**(-2) + 32*m1**2*s*
     +    u1**(-1) - 32*m1**2*t1*u1**(-1) - 64*m1**4*mg**2*s*u1**(-3)
     +     + 96*m1**4*s*u1**(-2) - 32*m1**4*t1*u1**(-2) + 64*m1**6*s*
     +    u1**(-3) - 32*mg**2*s**(-1)*t1 - 32*mg**2*s**(-1)*t1**2*
     +    u1**(-1) + 32*mst1**2*s**(-1)*t1 + 32*mst1**2*s**(-1)*t1**2*
     +    u1**(-1) + 32*mst1**2*s*u1**(-1) + 32*mst1**2*t1*u1**(-1) + 
     +    32*mst1**2 )
      MM_s = MM_s + SCC(3,3)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*mg**2*s*u1**(-2) - 32*m1**2*mg**2*
     +    t1*u1**(-2) - 64*m1**2*mst1**2*s*u1**(-2) - 32*m1**2*mst1**2*
     +    t1*u1**(-2) - 32*m1**2*mst1**2*u1**(-1) - 32*m1**2*s**(-1)*t1
     +     - 32*m1**2*s**(-1)*t1**2*u1**(-1) + 64*m1**4*mg**2*s*
     +    u1**(-3) + 96*m1**4*t1*u1**(-2) + 32*m1**4*u1**(-1) - 64*
     +    m1**6*s*u1**(-3) )
      MM_s = MM_s + SCC(3,3)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*mg*s**(-1)*t1 - 16*h_l*st*cb*mg*
     +    s**(-1)*u1 - 16*h_l*st*cb*mg + 32*h_l*ct*cb*m1*mg**2*s**(-1)
     +     - 32*h_l*ct*cb*m1*msb1**2*s**(-1) - 16*h_l*ct*cb*m1*s**(-1)*
     +    u1 - 16*h_l*ct*cb*m1 + 32*h_r*st*sb*m1*mg**2*s**(-1) - 32*h_r
     +    *st*sb*m1*msb1**2*s**(-1) - 16*h_r*st*sb*m1*s**(-1)*u1 - 16*
     +    h_r*st*sb*m1 - 32*h_r*ct*sb*mg*s**(-1)*t1 - 16*h_r*ct*sb*mg*
     +    s**(-1)*u1 - 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(3,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*mg*s**(-1)*t1 + 16*h_l*st*sb*mg*s**(-1)*
     +    u1 + 16*h_l*st*sb*mg - 32*h_l*ct*sb*m1*mg**2*s**(-1) + 32*h_l
     +    *ct*sb*m1*msb2**2*s**(-1) + 16*h_l*ct*sb*m1*s**(-1)*u1 + 16*
     +    h_l*ct*sb*m1 + 32*h_r*st*cb*m1*mg**2*s**(-1) - 32*h_r*st*cb*
     +    m1*msb2**2*s**(-1) - 16*h_r*st*cb*m1*s**(-1)*u1 - 16*h_r*st*
     +    cb*m1 - 32*h_r*ct*cb*mg*s**(-1)*t1 - 16*h_r*ct*cb*mg*s**(-1)*
     +    u1 - 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(3,4)*susy_ttg*(h_r2+h_l2)*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*s2t*m1*mg*s**(-1)*t1 - 32*s2t*m1*mg*s**(-1)*
     +    t1**2*u1**(-1) + 16*s2t*m1*mg*s*u1**(-1) - 16*s2t*m1*mg*t1*
     +    u1**(-1) + 64*s2t*m1**3*mg*s*u1**(-2) + 64*s2t*m1**3*mg*t1*
     +    u1**(-2) + 32*s2t*m1**3*mg*u1**(-1) )
      MM_s = MM_s + SCC(3,4)*susy_ttg*(h_r2*ct2+h_l2*st2)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*m1**2*mg**2*s*u1**(-2) - 32*m1**2*mg**2*
     +    t1*u1**(-2) - 64*m1**2*mst2**2*s*u1**(-2) - 32*m1**2*mst2**2*
     +    t1*u1**(-2) - 32*m1**2*mst2**2*u1**(-1) - 32*m1**2*s**(-1)*t1
     +     - 32*m1**2*s**(-1)*t1**2*u1**(-1) + 64*m1**4*mg**2*s*
     +    u1**(-3) + 96*m1**4*t1*u1**(-2) + 32*m1**4*u1**(-1) - 64*
     +    m1**6*s*u1**(-3) )
      MM_s = MM_s + SCC(3,4)*susy_ttg*(h_r2*st2+h_l2*ct2)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*m1**2*mg**2*t1*u1**(-2) + 32*m1**2*mg**2
     +    *u1**(-1) - 32*m1**2*mst2**2*t1*u1**(-2) + 32*m1**2*s*
     +    u1**(-1) - 32*m1**2*t1*u1**(-1) - 64*m1**4*mg**2*s*u1**(-3)
     +     + 96*m1**4*s*u1**(-2) - 32*m1**4*t1*u1**(-2) + 64*m1**6*s*
     +    u1**(-3) - 32*mg**2*s**(-1)*t1 - 32*mg**2*s**(-1)*t1**2*
     +    u1**(-1) + 32*mst2**2*s**(-1)*t1 + 32*mst2**2*s**(-1)*t1**2*
     +    u1**(-1) + 32*mst2**2*s*u1**(-1) + 32*mst2**2*t1*u1**(-1) + 
     +    32*mst2**2 )
      MM_s = MM_s + SCC(3,4)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*mg**2*s**(-1) + 32*h_l*st*cb*m1*
     +    msb1**2*s**(-1) + 16*h_l*st*cb*m1*s**(-1)*u1 + 16*h_l*st*cb*
     +    m1 - 32*h_l*ct*cb*mg*s**(-1)*t1 - 16*h_l*ct*cb*mg*s**(-1)*u1
     +     - 16*h_l*ct*cb*mg + 32*h_r*st*sb*mg*s**(-1)*t1 + 16*h_r*st*
     +    sb*mg*s**(-1)*u1 + 16*h_r*st*sb*mg + 32*h_r*ct*sb*m1*mg**2*
     +    s**(-1) - 32*h_r*ct*sb*m1*msb1**2*s**(-1) - 16*h_r*ct*sb*m1*
     +    s**(-1)*u1 - 16*h_r*ct*sb*m1 )
      MM_s = MM_s + SCC(3,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*mg**2*s**(-1) - 32*h_l*st*sb*m1*
     +    msb2**2*s**(-1) - 16*h_l*st*sb*m1*s**(-1)*u1 - 16*h_l*st*sb*
     +    m1 + 32*h_l*ct*sb*mg*s**(-1)*t1 + 16*h_l*ct*sb*mg*s**(-1)*u1
     +     + 16*h_l*ct*sb*mg + 32*h_r*st*cb*mg*s**(-1)*t1 + 16*h_r*st*
     +    cb*mg*s**(-1)*u1 + 16*h_r*st*cb*mg + 32*h_r*ct*cb*m1*mg**2*
     +    s**(-1) - 32*h_r*ct*cb*m1*msb2**2*s**(-1) - 16*h_r*ct*cb*m1*
     +    s**(-1)*u1 - 16*h_r*ct*cb*m1 )
      MM_s = MM_s + SCC(3,5)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*s2t*m1*mg*s*u1**(-1) - 16*s2t*m1*mg )
      MM_s = MM_s + SCC(3,5)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*m1**2*mg**2*t1*u1**(-2) + 32*m1**2
     +    *mg**2*u1**(-1) - 64*m1**4*mg**2*s*u1**(-3) + 32*mg**2*s*
     +    u1**(-1) + 32*mg**2*t1*u1**(-1) + 32*mg**2 )
      MM_s = MM_s + SCC(3,5)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*m1**2*mg**2*s*u1**(-2) - 64*
     +    m1**2*mg**2*t1*u1**(-2) - 32*m1**2*mg**2*u1**(-1) + 64*m1**4*
     +    mg**2*s*u1**(-3) )
      MM_s = MM_s + SCC(3,5)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*mg*s**(-1)*t1 + 16*h_l*st*cb*mg*s**(-1)
     +    *u1 + 16*h_l*st*cb*mg + 32*h_l*ct*cb*m1*mg**2*s**(-1) - 32*
     +    h_l*ct*cb*m1*msb1**2*s**(-1) + 16*h_l*ct*cb*m1 + 32*h_r*st*sb
     +    *m1*mg**2*s**(-1) - 32*h_r*st*sb*m1*msb1**2*s**(-1) + 16*h_r*
     +    st*sb*m1 + 32*h_r*ct*sb*mg*s**(-1)*t1 + 16*h_r*ct*sb*mg*
     +    s**(-1)*u1 + 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(3,5)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*mg*s**(-1)*t1 - 16*h_l*st*sb*mg*
     +    s**(-1)*u1 - 16*h_l*st*sb*mg - 32*h_l*ct*sb*m1*mg**2*s**(-1)
     +     + 32*h_l*ct*sb*m1*msb2**2*s**(-1) - 16*h_l*ct*sb*m1 + 32*h_r
     +    *st*cb*m1*mg**2*s**(-1) - 32*h_r*st*cb*m1*msb2**2*s**(-1) + 
     +    16*h_r*st*cb*m1 + 32*h_r*ct*cb*mg*s**(-1)*t1 + 16*h_r*ct*cb*
     +    mg*s**(-1)*u1 + 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(3,6)*susy_ttg*(h_r2+h_l2)*Nc**2*Cf*Pi**2*
     + alphas**2*prefac * ( 16*s2t*m1*mg*s*u1**(-1) + 16*s2t*m1*mg )
      MM_s = MM_s + SCC(3,6)*susy_ttg*(h_r2*ct2+h_l2*st2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*m1**2*mg**2*s*u1**(-2) - 64*
     +    m1**2*mg**2*t1*u1**(-2) - 32*m1**2*mg**2*u1**(-1) + 64*m1**4*
     +    mg**2*s*u1**(-3) )
      MM_s = MM_s + SCC(3,6)*susy_ttg*(h_r2*st2+h_l2*ct2)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*m1**2*mg**2*t1*u1**(-2) + 32*m1**2
     +    *mg**2*u1**(-1) - 64*m1**4*mg**2*s*u1**(-3) + 32*mg**2*s*
     +    u1**(-1) + 32*mg**2*t1*u1**(-1) + 32*mg**2 )
      MM_s = MM_s + SCC(3,6)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*cb*m1*mg**2*s**(-1) + 32*h_l*st*cb*m1*
     +    msb1**2*s**(-1) - 16*h_l*st*cb*m1 + 32*h_l*ct*cb*mg*s**(-1)*
     +    t1 + 16*h_l*ct*cb*mg*s**(-1)*u1 + 16*h_l*ct*cb*mg - 32*h_r*st
     +    *sb*mg*s**(-1)*t1 - 16*h_r*st*sb*mg*s**(-1)*u1 - 16*h_r*st*sb
     +    *mg + 32*h_r*ct*sb*m1*mg**2*s**(-1) - 32*h_r*ct*sb*m1*msb1**2
     +    *s**(-1) + 16*h_r*ct*sb*m1 )
      MM_s = MM_s + SCC(3,6)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*sb*m1*mg**2*s**(-1) - 32*h_l*st*sb*m1*
     +    msb2**2*s**(-1) + 16*h_l*st*sb*m1 - 32*h_l*ct*sb*mg*s**(-1)*
     +    t1 - 16*h_l*ct*sb*mg*s**(-1)*u1 - 16*h_l*ct*sb*mg - 32*h_r*st
     +    *cb*mg*s**(-1)*t1 - 16*h_r*st*cb*mg*s**(-1)*u1 - 16*h_r*st*cb
     +    *mg + 32*h_r*ct*cb*m1*mg**2*s**(-1) - 32*h_r*ct*cb*m1*msb2**2
     +    *s**(-1) + 16*h_r*ct*cb*m1 )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_l*
     +    st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 96*h_l*st*cb*m1**2*mg*
     +    u1**(-2)*u2 - 32*h_l*st*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*
     +    h_l*st*cb*m2**2*mg*u1**(-2)*u2 - 32*h_l*st*cb*mg*s**(-1)*t1*
     +    u1**(-1)*u2 - 32*h_l*st*cb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_l
     +    *st*cb*mg*s*t1*u1**(-2) + 32*h_l*st*cb*mg*s*u1**(-2)*u2 - 32*
     +    h_l*st*cb*mg*t1*u1**(-2)*u2 - 64*h_l*st*cb*mg*t1*u1**(-1) - 
     +    64*h_l*st*cb*mg*t1**2*u1**(-2) + 32*h_l*st*cb*mg*u1**(-1)*u2
     +     - 32*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 32*h_l*ct*cb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 64*h_l*ct*cb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_l*ct*cb*
     +    m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_l*ct*cb*m1*mg**2*s*
     +    u1**(-2) - 96*h_l*ct*cb*m1*mg**2*t1*u1**(-2) - 96*h_l*ct*cb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_l*ct*cb*m1*mg**2*u1**(-1) + 64*
     +    h_l*ct*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 64*h_l*ct*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 + 96*h_l
     +    *ct*cb*m1*msb1**2*s*u1**(-2) + 96*h_l*ct*cb*m1*msb1**2*t1*
     +    u1**(-2) + 96*h_l*ct*cb*m1*msb1**2*u1**(-2)*u2 + 96*h_l*ct*cb
     +    *m1*msb1**2*u1**(-1) - 32*h_l*ct*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) - 64*h_l*ct*cb*m1**3*mg**2*u1**(-2) + 32*h_l*ct*cb*
     +    m1**3*msb1**2*s**(-1)*u1**(-1) + 64*h_l*ct*cb*m1**3*msb1**2*
     +    u1**(-2) - 32*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_r*st*sb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_r
     +    *st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*sb*m1*mg**2*s
     +    *u1**(-2) - 96*h_r*st*sb*m1*mg**2*t1*u1**(-2) - 96*h_r*st*sb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_r*st*sb*m1*mg**2*u1**(-1) + 64*
     +    h_r*st*sb*m1*msb1**2*s**(-1)*t1*u1**(-1) + 64*h_r*st*sb*m1*
     +    msb1**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*sb*m1*msb1**2*s*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 96*h_r*st*sb*m1*msb1**2*t1*u1**(-2) + 96*h_r*st*sb*
     +    m1*msb1**2*u1**(-2)*u2 + 96*h_r*st*sb*m1*msb1**2*u1**(-1) - 
     +    32*h_r*st*sb*m1**3*mg**2*s**(-1)*u1**(-1) - 64*h_r*st*sb*
     +    m1**3*mg**2*u1**(-2) + 32*h_r*st*sb*m1**3*msb1**2*s**(-1)*
     +    u1**(-1) + 64*h_r*st*sb*m1**3*msb1**2*u1**(-2) + 32*h_r*ct*sb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_r*ct*sb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 + 96*h_r*ct*sb*m1**2*mg*u1**(-2)*u2 - 32*h_r*ct*
     +    sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*sb*m2**2*mg*
     +    u1**(-2)*u2 - 32*h_r*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r
     +    *ct*sb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*sb*mg*s*t1*
     +    u1**(-2) + 32*h_r*ct*sb*mg*s*u1**(-2)*u2 - 32*h_r*ct*sb*mg*t1
     +    *u1**(-2)*u2 - 64*h_r*ct*sb*mg*t1*u1**(-1) - 64*h_r*ct*sb*mg*
     +    t1**2*u1**(-2) + 32*h_r*ct*sb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_l
     +    *st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 96*h_l*st*cb*m1**2*mg*
     +    u1**(-2)*u2 + 32*h_l*st*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*
     +    h_l*st*cb*m2**2*mg*u1**(-2)*u2 + 32*h_l*st*cb*mg*s**(-1)*t1*
     +    u1**(-1)*u2 + 32*h_l*st*cb*mg*s**(-1)*t1**2*u1**(-1) + 64*h_l
     +    *st*cb*mg*s*t1*u1**(-2) - 32*h_l*st*cb*mg*s*u1**(-2)*u2 + 32*
     +    h_l*st*cb*mg*t1*u1**(-2)*u2 + 64*h_l*st*cb*mg*t1*u1**(-1) + 
     +    64*h_l*st*cb*mg*t1**2*u1**(-2) - 32*h_l*st*cb*mg*u1**(-1)*u2
     +     + 32*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 32*h_l*ct*cb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 64*h_l*ct*cb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_l*ct*cb*
     +    m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_l*ct*cb*m1*mg**2*s*
     +    u1**(-2) + 96*h_l*ct*cb*m1*mg**2*t1*u1**(-2) + 96*h_l*ct*cb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_l*ct*cb*m1*mg**2*u1**(-1) - 64*
     +    h_l*ct*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 64*h_l*ct*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 96*
     +    h_l*ct*cb*m1*msb1**2*s*u1**(-2) - 96*h_l*ct*cb*m1*msb1**2*t1*
     +    u1**(-2) - 96*h_l*ct*cb*m1*msb1**2*u1**(-2)*u2 - 96*h_l*ct*cb
     +    *m1*msb1**2*u1**(-1) + 32*h_l*ct*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) + 64*h_l*ct*cb*m1**3*mg**2*u1**(-2) - 32*h_l*ct*cb*
     +    m1**3*msb1**2*s**(-1)*u1**(-1) - 64*h_l*ct*cb*m1**3*msb1**2*
     +    u1**(-2) + 32*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 64*h_r*st*sb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_r
     +    *st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*sb*m1*mg**2*s
     +    *u1**(-2) + 96*h_r*st*sb*m1*mg**2*t1*u1**(-2) + 96*h_r*st*sb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_r*st*sb*m1*mg**2*u1**(-1) - 64*
     +    h_r*st*sb*m1*msb1**2*s**(-1)*t1*u1**(-1) - 64*h_r*st*sb*m1*
     +    msb1**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*sb*m1*msb1**2*s*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 96*h_r*st*sb*m1*msb1**2*t1*u1**(-2) - 96*h_r*st*sb
     +    *m1*msb1**2*u1**(-2)*u2 - 96*h_r*st*sb*m1*msb1**2*u1**(-1) + 
     +    32*h_r*st*sb*m1**3*mg**2*s**(-1)*u1**(-1) + 64*h_r*st*sb*
     +    m1**3*mg**2*u1**(-2) - 32*h_r*st*sb*m1**3*msb1**2*s**(-1)*
     +    u1**(-1) - 64*h_r*st*sb*m1**3*msb1**2*u1**(-2) - 32*h_r*ct*sb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_r*ct*sb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 - 96*h_r*ct*sb*m1**2*mg*u1**(-2)*u2 + 32*h_r*ct*
     +    sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*m2**2*mg*
     +    u1**(-2)*u2 + 32*h_r*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r
     +    *ct*sb*mg*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*sb*mg*s*t1*
     +    u1**(-2) - 32*h_r*ct*sb*mg*s*u1**(-2)*u2 + 32*h_r*ct*sb*mg*t1
     +    *u1**(-2)*u2 + 64*h_r*ct*sb*mg*t1*u1**(-1) + 64*h_r*ct*sb*mg*
     +    t1**2*u1**(-2) - 32*h_r*ct*sb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*st*cb*mg*s**(-1)*t1 - 16*h_l*st*cb*mg*s*
     +    u1**(-1) - 16*h_l*st*cb*mg*t1*u1**(-1) - 16*h_l*st*cb*mg - 32
     +    *h_l*ct*cb*m1*mg**2*s**(-1) - 32*h_l*ct*cb*m1*mg**2*u1**(-1)
     +     + 32*h_l*ct*cb*m1*msb1**2*s**(-1) + 32*h_l*ct*cb*m1*msb1**2*
     +    u1**(-1) - 16*h_l*ct*cb*m1*s*u1**(-1) - 16*h_l*ct*cb*m1*t1*
     +    u1**(-1) - 32*h_r*st*sb*m1*mg**2*s**(-1) - 32*h_r*st*sb*m1*
     +    mg**2*u1**(-1) + 32*h_r*st*sb*m1*msb1**2*s**(-1) + 32*h_r*st*
     +    sb*m1*msb1**2*u1**(-1) - 16*h_r*st*sb*m1*s*u1**(-1) - 16*h_r*
     +    st*sb*m1*t1*u1**(-1) - 16*h_r*ct*sb*mg*s**(-1)*t1 - 16*h_r*ct
     +    *sb*mg*s*u1**(-1) - 16*h_r*ct*sb*mg*t1*u1**(-1) - 16*h_r*ct*
     +    sb*mg )
      MM_s = MM_s + SCC(4,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*cb*mg*s**(-1)*t1 - 16*h_l*st*cb*mg*s*
     +    u1**(-1) - 16*h_l*st*cb*mg*t1*u1**(-1) - 16*h_l*st*cb*mg + 32
     +    *h_l*ct*cb*m1*mg**2*s**(-1) + 32*h_l*ct*cb*m1*mg**2*u1**(-1)
     +     - 32*h_l*ct*cb*m1*msb1**2*s**(-1) - 32*h_l*ct*cb*m1*msb1**2*
     +    u1**(-1) + 16*h_l*ct*cb*m1*s**(-1)*t1 + 16*h_l*ct*cb*m1*s*
     +    u1**(-1) + 16*h_l*ct*cb*m1*t1*u1**(-1) + 16*h_l*ct*cb*m1 + 32
     +    *h_r*st*sb*m1*mg**2*s**(-1) + 32*h_r*st*sb*m1*mg**2*u1**(-1)
     +     - 32*h_r*st*sb*m1*msb1**2*s**(-1) - 32*h_r*st*sb*m1*msb1**2*
     +    u1**(-1) + 16*h_r*st*sb*m1*s**(-1)*t1 + 16*h_r*st*sb*m1*s*
     +    u1**(-1) + 16*h_r*st*sb*m1*t1*u1**(-1) + 16*h_r*st*sb*m1 - 16
     +    *h_r*ct*sb*mg*s**(-1)*t1 - 16*h_r*ct*sb*mg*s*u1**(-1) - 16*
     +    h_r*ct*sb*mg*t1*u1**(-1) - 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1**2*mg*s**(-1)*t1*u1**(-1) - 64*
     +    h_l*st*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 96*h_l*st*sb*m1**2*
     +    mg*u1**(-2)*u2 + 32*h_l*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 
     +    32*h_l*st*sb*m2**2*mg*u1**(-2)*u2 + 32*h_l*st*sb*mg*s**(-1)*
     +    t1*u1**(-1)*u2 + 32*h_l*st*sb*mg*s**(-1)*t1**2*u1**(-1) + 64*
     +    h_l*st*sb*mg*s*t1*u1**(-2) - 32*h_l*st*sb*mg*s*u1**(-2)*u2 + 
     +    32*h_l*st*sb*mg*t1*u1**(-2)*u2 + 64*h_l*st*sb*mg*t1*u1**(-1)
     +     + 64*h_l*st*sb*mg*t1**2*u1**(-2) - 32*h_l*st*sb*mg*u1**(-1)*
     +    u2 + 32*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 32*h_l*ct*sb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 64*h_l*ct*sb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_l*ct*sb*
     +    m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_l*ct*sb*m1*mg**2*s*
     +    u1**(-2) + 96*h_l*ct*sb*m1*mg**2*t1*u1**(-2) + 96*h_l*ct*sb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_l*ct*sb*m1*mg**2*u1**(-1) - 64*
     +    h_l*ct*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 64*h_l*ct*sb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 96*
     +    h_l*ct*sb*m1*msb2**2*s*u1**(-2) - 96*h_l*ct*sb*m1*msb2**2*t1*
     +    u1**(-2) - 96*h_l*ct*sb*m1*msb2**2*u1**(-2)*u2 - 96*h_l*ct*sb
     +    *m1*msb2**2*u1**(-1) + 32*h_l*ct*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) + 64*h_l*ct*sb*m1**3*mg**2*u1**(-2) - 32*h_l*ct*sb*
     +    m1**3*msb2**2*s**(-1)*u1**(-1) - 64*h_l*ct*sb*m1**3*msb2**2*
     +    u1**(-2) - 32*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_r*st*cb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_r
     +    *st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*cb*m1*mg**2*s
     +    *u1**(-2) - 96*h_r*st*cb*m1*mg**2*t1*u1**(-2) - 96*h_r*st*cb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_r*st*cb*m1*mg**2*u1**(-1) + 64*
     +    h_r*st*cb*m1*msb2**2*s**(-1)*t1*u1**(-1) + 64*h_r*st*cb*m1*
     +    msb2**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*cb*m1*msb2**2*s*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 96*h_r*st*cb*m1*msb2**2*t1*u1**(-2) + 96*h_r*st*cb*
     +    m1*msb2**2*u1**(-2)*u2 + 96*h_r*st*cb*m1*msb2**2*u1**(-1) - 
     +    32*h_r*st*cb*m1**3*mg**2*s**(-1)*u1**(-1) - 64*h_r*st*cb*
     +    m1**3*mg**2*u1**(-2) + 32*h_r*st*cb*m1**3*msb2**2*s**(-1)*
     +    u1**(-1) + 64*h_r*st*cb*m1**3*msb2**2*u1**(-2) + 32*h_r*ct*cb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_r*ct*cb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 + 96*h_r*ct*cb*m1**2*mg*u1**(-2)*u2 - 32*h_r*ct*
     +    cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*cb*m2**2*mg*
     +    u1**(-2)*u2 - 32*h_r*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r
     +    *ct*cb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*cb*mg*s*t1*
     +    u1**(-2) + 32*h_r*ct*cb*mg*s*u1**(-2)*u2 - 32*h_r*ct*cb*mg*t1
     +    *u1**(-2)*u2 - 64*h_r*ct*cb*mg*t1*u1**(-1) - 64*h_r*ct*cb*mg*
     +    t1**2*u1**(-2) + 32*h_r*ct*cb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_l*st
     +    *sb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 96*h_l*st*sb*m1**2*mg*
     +    u1**(-2)*u2 - 32*h_l*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*
     +    h_l*st*sb*m2**2*mg*u1**(-2)*u2 - 32*h_l*st*sb*mg*s**(-1)*t1*
     +    u1**(-1)*u2 - 32*h_l*st*sb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_l
     +    *st*sb*mg*s*t1*u1**(-2) + 32*h_l*st*sb*mg*s*u1**(-2)*u2 - 32*
     +    h_l*st*sb*mg*t1*u1**(-2)*u2 - 64*h_l*st*sb*mg*t1*u1**(-1) - 
     +    64*h_l*st*sb*mg*t1**2*u1**(-2) + 32*h_l*st*sb*mg*u1**(-1)*u2
     +     - 32*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     + 32*h_l*ct*sb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*u2**(-1)
     +     - 64*h_l*ct*sb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_l*ct*sb*
     +    m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_l*ct*sb*m1*mg**2*s*
     +    u1**(-2) - 96*h_l*ct*sb*m1*mg**2*t1*u1**(-2) - 96*h_l*ct*sb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_l*ct*sb*m1*mg**2*u1**(-1) + 64*
     +    h_l*ct*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 64*h_l*ct*sb*m1*msb2**2*s**(-1)*u1**(-1)*u2 + 96*h_l*
     +    ct*sb*m1*msb2**2*s*u1**(-2) + 96*h_l*ct*sb*m1*msb2**2*t1*
     +    u1**(-2) + 96*h_l*ct*sb*m1*msb2**2*u1**(-2)*u2 + 96*h_l*ct*sb
     +    *m1*msb2**2*u1**(-1) - 32*h_l*ct*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) - 64*h_l*ct*sb*m1**3*mg**2*u1**(-2) + 32*h_l*ct*sb*
     +    m1**3*msb2**2*s**(-1)*u1**(-1) + 64*h_l*ct*sb*m1**3*msb2**2*
     +    u1**(-2) + 32*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 64*h_r*st*cb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_r
     +    *st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_r*st*cb*m1*mg**2*s
     +    *u1**(-2) + 96*h_r*st*cb*m1*mg**2*t1*u1**(-2) + 96*h_r*st*cb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_r*st*cb*m1*mg**2*u1**(-1) - 64*
     +    h_r*st*cb*m1*msb2**2*s**(-1)*t1*u1**(-1) - 64*h_r*st*cb*m1*
     +    msb2**2*s**(-1)*u1**(-1)*u2 - 96*h_r*st*cb*m1*msb2**2*s*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 96*h_r*st*cb*m1*msb2**2*t1*u1**(-2) - 96*h_r*st*cb
     +    *m1*msb2**2*u1**(-2)*u2 - 96*h_r*st*cb*m1*msb2**2*u1**(-1) + 
     +    32*h_r*st*cb*m1**3*mg**2*s**(-1)*u1**(-1) + 64*h_r*st*cb*
     +    m1**3*mg**2*u1**(-2) - 32*h_r*st*cb*m1**3*msb2**2*s**(-1)*
     +    u1**(-1) - 64*h_r*st*cb*m1**3*msb2**2*u1**(-2) - 32*h_r*ct*cb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_r*ct*cb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 - 96*h_r*ct*cb*m1**2*mg*u1**(-2)*u2 + 32*h_r*ct*
     +    cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m2**2*mg*
     +    u1**(-2)*u2 + 32*h_r*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r
     +    *ct*cb*mg*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*cb*mg*s*t1*
     +    u1**(-2) - 32*h_r*ct*cb*mg*s*u1**(-2)*u2 + 32*h_r*ct*cb*mg*t1
     +    *u1**(-2)*u2 + 64*h_r*ct*cb*mg*t1*u1**(-1) + 64*h_r*ct*cb*mg*
     +    t1**2*u1**(-2) - 32*h_r*ct*cb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*st*sb*mg*s**(-1)*t1 + 16*h_l*st*sb*mg*s*
     +    u1**(-1) + 16*h_l*st*sb*mg*t1*u1**(-1) + 16*h_l*st*sb*mg + 32
     +    *h_l*ct*sb*m1*mg**2*s**(-1) + 32*h_l*ct*sb*m1*mg**2*u1**(-1)
     +     - 32*h_l*ct*sb*m1*msb2**2*s**(-1) - 32*h_l*ct*sb*m1*msb2**2*
     +    u1**(-1) + 16*h_l*ct*sb*m1*s*u1**(-1) + 16*h_l*ct*sb*m1*t1*
     +    u1**(-1) - 32*h_r*st*cb*m1*mg**2*s**(-1) - 32*h_r*st*cb*m1*
     +    mg**2*u1**(-1) + 32*h_r*st*cb*m1*msb2**2*s**(-1) + 32*h_r*st*
     +    cb*m1*msb2**2*u1**(-1) - 16*h_r*st*cb*m1*s*u1**(-1) - 16*h_r*
     +    st*cb*m1*t1*u1**(-1) - 16*h_r*ct*cb*mg*s**(-1)*t1 - 16*h_r*ct
     +    *cb*mg*s*u1**(-1) - 16*h_r*ct*cb*mg*t1*u1**(-1) - 16*h_r*ct*
     +    cb*mg )
      MM_s = MM_s + SCC(4,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*sb*mg*s**(-1)*t1 + 16*h_l*st*sb*mg*s*
     +    u1**(-1) + 16*h_l*st*sb*mg*t1*u1**(-1) + 16*h_l*st*sb*mg - 32
     +    *h_l*ct*sb*m1*mg**2*s**(-1) - 32*h_l*ct*sb*m1*mg**2*u1**(-1)
     +     + 32*h_l*ct*sb*m1*msb2**2*s**(-1) + 32*h_l*ct*sb*m1*msb2**2*
     +    u1**(-1) - 16*h_l*ct*sb*m1*s**(-1)*t1 - 16*h_l*ct*sb*m1*s*
     +    u1**(-1) - 16*h_l*ct*sb*m1*t1*u1**(-1) - 16*h_l*ct*sb*m1 + 32
     +    *h_r*st*cb*m1*mg**2*s**(-1) + 32*h_r*st*cb*m1*mg**2*u1**(-1)
     +     - 32*h_r*st*cb*m1*msb2**2*s**(-1) - 32*h_r*st*cb*m1*msb2**2*
     +    u1**(-1) + 16*h_r*st*cb*m1*s**(-1)*t1 + 16*h_r*st*cb*m1*s*
     +    u1**(-1) + 16*h_r*st*cb*m1*t1*u1**(-1) + 16*h_r*st*cb*m1 - 16
     +    *h_r*ct*cb*mg*s**(-1)*t1 - 16*h_r*ct*cb*mg*s*u1**(-1) - 16*
     +    h_r*ct*cb*mg*t1*u1**(-1) - 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_l*st*cb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 64*h_l*st*cb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_l
     +    *st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_l*st*cb*m1*mg**2*s
     +    *u1**(-2) + 96*h_l*st*cb*m1*mg**2*t1*u1**(-2) + 96*h_l*st*cb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_l*st*cb*m1*mg**2*u1**(-1) - 64*
     +    h_l*st*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) - 64*h_l*st*cb*m1*
     +    msb1**2*s**(-1)*u1**(-1)*u2 - 96*h_l*st*cb*m1*msb1**2*s*
     +    u1**(-2) - 96*h_l*st*cb*m1*msb1**2*t1*u1**(-2) - 96*h_l*st*cb
     +    *m1*msb1**2*u1**(-2)*u2 - 96*h_l*st*cb*m1*msb1**2*u1**(-1) + 
     +    32*h_l*st*cb*m1**3*mg**2*s**(-1)*u1**(-1) + 64*h_l*st*cb*
     +    m1**3*mg**2*u1**(-2) - 32*h_l*st*cb*m1**3*msb1**2*s**(-1)*
     +    u1**(-1) - 64*h_l*st*cb*m1**3*msb1**2*u1**(-2) + 32*h_l*ct*cb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_l*ct*cb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 96*h_l*ct*cb*m1**2*mg*u1**(-2)*u2 - 32*h_l*ct*cb*
     +    m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*m2**2*mg*u1**(-2)
     +    *u2 - 32*h_l*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*cb*
     +    mg*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*mg*s*t1*u1**(-2) + 
     +    32*h_l*ct*cb*mg*s*u1**(-2)*u2 - 32*h_l*ct*cb*mg*t1*u1**(-2)*
     +    u2 - 64*h_l*ct*cb*mg*t1*u1**(-1) - 64*h_l*ct*cb*mg*t1**2*
     +    u1**(-2) + 32*h_l*ct*cb*mg*u1**(-1)*u2 - 32*h_r*st*sb*m1**2*
     +    mg*s**(-1)*t1*u1**(-1) - 64*h_r*st*sb*m1**2*mg*s**(-1)*
     +    u1**(-1)*u2 - 96*h_r*st*sb*m1**2*mg*u1**(-2)*u2 + 32*h_r*st*
     +    sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*m2**2*mg*
     +    u1**(-2)*u2 + 32*h_r*st*sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r
     +    *st*sb*mg*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*mg*s*t1*
     +    u1**(-2) - 32*h_r*st*sb*mg*s*u1**(-2)*u2 + 32*h_r*st*sb*mg*t1
     +    *u1**(-2)*u2 + 64*h_r*st*sb*mg*t1*u1**(-1) + 64*h_r*st*sb*mg*
     +    t1**2*u1**(-2) )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_r*st*sb*mg*u1**(-1)*u2 - 32*h_r*ct*sb*m1*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 32*h_r*ct*sb*m1*
     +    m2**2*msb1**2*s**(-1)*t1*u1**(-1)*u2**(-1) - 64*h_r*ct*sb*m1*
     +    mg**2*s**(-1)*t1*u1**(-1) - 64*h_r*ct*sb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 96*h_r*ct*sb*m1*mg**2*s*u1**(-2) - 96*h_r*ct*sb
     +    *m1*mg**2*t1*u1**(-2) - 96*h_r*ct*sb*m1*mg**2*u1**(-2)*u2 - 
     +    96*h_r*ct*sb*m1*mg**2*u1**(-1) + 64*h_r*ct*sb*m1*msb1**2*
     +    s**(-1)*t1*u1**(-1) + 64*h_r*ct*sb*m1*msb1**2*s**(-1)*
     +    u1**(-1)*u2 + 96*h_r*ct*sb*m1*msb1**2*s*u1**(-2) + 96*h_r*ct*
     +    sb*m1*msb1**2*t1*u1**(-2) + 96*h_r*ct*sb*m1*msb1**2*u1**(-2)*
     +    u2 + 96*h_r*ct*sb*m1*msb1**2*u1**(-1) - 32*h_r*ct*sb*m1**3*
     +    mg**2*s**(-1)*u1**(-1) - 64*h_r*ct*sb*m1**3*mg**2*u1**(-2) + 
     +    32*h_r*ct*sb*m1**3*msb1**2*s**(-1)*u1**(-1) + 64*h_r*ct*sb*
     +    m1**3*msb1**2*u1**(-2) )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_l*st*cb*m1*m2**2*msb1**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_l*st*cb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_l
     +    *st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_l*st*cb*m1*mg**2*s
     +    *u1**(-2) - 96*h_l*st*cb*m1*mg**2*t1*u1**(-2) - 96*h_l*st*cb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_l*st*cb*m1*mg**2*u1**(-1) + 64*
     +    h_l*st*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) + 64*h_l*st*cb*m1*
     +    msb1**2*s**(-1)*u1**(-1)*u2 + 96*h_l*st*cb*m1*msb1**2*s*
     +    u1**(-2) + 96*h_l*st*cb*m1*msb1**2*t1*u1**(-2) + 96*h_l*st*cb
     +    *m1*msb1**2*u1**(-2)*u2 + 96*h_l*st*cb*m1*msb1**2*u1**(-1) - 
     +    32*h_l*st*cb*m1**3*mg**2*s**(-1)*u1**(-1) - 64*h_l*st*cb*
     +    m1**3*mg**2*u1**(-2) + 32*h_l*st*cb*m1**3*msb1**2*s**(-1)*
     +    u1**(-1) + 64*h_l*st*cb*m1**3*msb1**2*u1**(-2) - 32*h_l*ct*cb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_l*ct*cb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 96*h_l*ct*cb*m1**2*mg*u1**(-2)*u2 + 32*h_l*ct*cb*
     +    m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*m2**2*mg*u1**(-2)
     +    *u2 + 32*h_l*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*cb*
     +    mg*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*mg*s*t1*u1**(-2) - 
     +    32*h_l*ct*cb*mg*s*u1**(-2)*u2 + 32*h_l*ct*cb*mg*t1*u1**(-2)*
     +    u2 + 64*h_l*ct*cb*mg*t1*u1**(-1) + 64*h_l*ct*cb*mg*t1**2*
     +    u1**(-2) - 32*h_l*ct*cb*mg*u1**(-1)*u2 + 32*h_r*st*sb*m1**2*
     +    mg*s**(-1)*t1*u1**(-1) + 64*h_r*st*sb*m1**2*mg*s**(-1)*
     +    u1**(-1)*u2 + 96*h_r*st*sb*m1**2*mg*u1**(-2)*u2 - 32*h_r*st*
     +    sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*st*sb*m2**2*mg*
     +    u1**(-2)*u2 - 32*h_r*st*sb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r
     +    *st*sb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb*mg*s*t1*
     +    u1**(-2) + 32*h_r*st*sb*mg*s*u1**(-2)*u2 - 32*h_r*st*sb*mg*t1
     +    *u1**(-2)*u2 - 64*h_r*st*sb*mg*t1*u1**(-1) - 64*h_r*st*sb*mg*
     +    t1**2*u1**(-2) )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*st*sb*mg*u1**(-1)*u2 + 32*h_r*ct*sb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2**(-1) - 32*h_r*ct*sb*m1*m2**2*
     +    msb1**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 64*h_r*ct*sb*m1*mg**2*
     +    s**(-1)*t1*u1**(-1) + 64*h_r*ct*sb*m1*mg**2*s**(-1)*u1**(-1)*
     +    u2 + 96*h_r*ct*sb*m1*mg**2*s*u1**(-2) + 96*h_r*ct*sb*m1*mg**2
     +    *t1*u1**(-2) + 96*h_r*ct*sb*m1*mg**2*u1**(-2)*u2 + 96*h_r*ct*
     +    sb*m1*mg**2*u1**(-1) - 64*h_r*ct*sb*m1*msb1**2*s**(-1)*t1*
     +    u1**(-1) - 64*h_r*ct*sb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 96*
     +    h_r*ct*sb*m1*msb1**2*s*u1**(-2) - 96*h_r*ct*sb*m1*msb1**2*t1*
     +    u1**(-2) - 96*h_r*ct*sb*m1*msb1**2*u1**(-2)*u2 - 96*h_r*ct*sb
     +    *m1*msb1**2*u1**(-1) + 32*h_r*ct*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) + 64*h_r*ct*sb*m1**3*mg**2*u1**(-2) - 32*h_r*ct*sb*
     +    m1**3*msb1**2*s**(-1)*u1**(-1) - 64*h_r*ct*sb*m1**3*msb1**2*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*mg**2*s**(-1) + 32*h_l*st*cb*m1*
     +    mg**2*u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s**(-1) - 32*h_l*st*
     +    cb*m1*msb1**2*u1**(-1) + 16*h_l*st*cb*m1*s*u1**(-1) + 16*h_l*
     +    st*cb*m1*t1*u1**(-1) - 16*h_l*ct*cb*mg*s**(-1)*t1 - 16*h_l*ct
     +    *cb*mg*s*u1**(-1) - 16*h_l*ct*cb*mg*t1*u1**(-1) - 16*h_l*ct*
     +    cb*mg + 16*h_r*st*sb*mg*s**(-1)*t1 + 16*h_r*st*sb*mg*s*
     +    u1**(-1) + 16*h_r*st*sb*mg*t1*u1**(-1) + 16*h_r*st*sb*mg - 32
     +    *h_r*ct*sb*m1*mg**2*s**(-1) - 32*h_r*ct*sb*m1*mg**2*u1**(-1)
     +     + 32*h_r*ct*sb*m1*msb1**2*s**(-1) + 32*h_r*ct*sb*m1*msb1**2*
     +    u1**(-1) - 16*h_r*ct*sb*m1*s*u1**(-1) - 16*h_r*ct*sb*m1*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(4,4)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*mg**2*s**(-1) - 32*h_l*st*cb*m1*
     +    mg**2*u1**(-1) + 32*h_l*st*cb*m1*msb1**2*s**(-1) + 32*h_l*st*
     +    cb*m1*msb1**2*u1**(-1) - 16*h_l*st*cb*m1*s**(-1)*t1 - 16*h_l*
     +    st*cb*m1*s*u1**(-1) - 16*h_l*st*cb*m1*t1*u1**(-1) - 16*h_l*st
     +    *cb*m1 - 16*h_l*ct*cb*mg*s**(-1)*t1 - 16*h_l*ct*cb*mg*s*
     +    u1**(-1) - 16*h_l*ct*cb*mg*t1*u1**(-1) - 16*h_l*ct*cb*mg + 16
     +    *h_r*st*sb*mg*s**(-1)*t1 + 16*h_r*st*sb*mg*s*u1**(-1) + 16*
     +    h_r*st*sb*mg*t1*u1**(-1) + 16*h_r*st*sb*mg + 32*h_r*ct*sb*m1*
     +    mg**2*s**(-1) + 32*h_r*ct*sb*m1*mg**2*u1**(-1) - 32*h_r*ct*sb
     +    *m1*msb1**2*s**(-1) - 32*h_r*ct*sb*m1*msb1**2*u1**(-1) + 16*
     +    h_r*ct*sb*m1*s**(-1)*t1 + 16*h_r*ct*sb*m1*s*u1**(-1) + 16*h_r
     +    *ct*sb*m1*t1*u1**(-1) + 16*h_r*ct*sb*m1 )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 32*h_l*st*sb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 64*h_l*st*sb*m1*mg**2*s**(-1)*t1*u1**(-1) - 64*h_l
     +    *st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 96*h_l*st*sb*m1*mg**2*s
     +    *u1**(-2) - 96*h_l*st*sb*m1*mg**2*t1*u1**(-2) - 96*h_l*st*sb*
     +    m1*mg**2*u1**(-2)*u2 - 96*h_l*st*sb*m1*mg**2*u1**(-1) + 64*
     +    h_l*st*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) + 64*h_l*st*sb*m1*
     +    msb2**2*s**(-1)*u1**(-1)*u2 + 96*h_l*st*sb*m1*msb2**2*s*
     +    u1**(-2) + 96*h_l*st*sb*m1*msb2**2*t1*u1**(-2) + 96*h_l*st*sb
     +    *m1*msb2**2*u1**(-2)*u2 + 96*h_l*st*sb*m1*msb2**2*u1**(-1) - 
     +    32*h_l*st*sb*m1**3*mg**2*s**(-1)*u1**(-1) - 64*h_l*st*sb*
     +    m1**3*mg**2*u1**(-2) + 32*h_l*st*sb*m1**3*msb2**2*s**(-1)*
     +    u1**(-1) + 64*h_l*st*sb*m1**3*msb2**2*u1**(-2) - 32*h_l*ct*sb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_l*ct*sb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 96*h_l*ct*sb*m1**2*mg*u1**(-2)*u2 + 32*h_l*ct*sb*
     +    m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*sb*m2**2*mg*u1**(-2)
     +    *u2 + 32*h_l*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*sb*
     +    mg*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*mg*s*t1*u1**(-2) - 
     +    32*h_l*ct*sb*mg*s*u1**(-2)*u2 + 32*h_l*ct*sb*mg*t1*u1**(-2)*
     +    u2 + 64*h_l*ct*sb*mg*t1*u1**(-1) + 64*h_l*ct*sb*mg*t1**2*
     +    u1**(-2) - 32*h_l*ct*sb*mg*u1**(-1)*u2 - 32*h_r*st*cb*m1**2*
     +    mg*s**(-1)*t1*u1**(-1) - 64*h_r*st*cb*m1**2*mg*s**(-1)*
     +    u1**(-1)*u2 - 96*h_r*st*cb*m1**2*mg*u1**(-2)*u2 + 32*h_r*st*
     +    cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*m2**2*mg*
     +    u1**(-2)*u2 + 32*h_r*st*cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r
     +    *st*cb*mg*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*mg*s*t1*
     +    u1**(-2) - 32*h_r*st*cb*mg*s*u1**(-2)*u2 + 32*h_r*st*cb*mg*t1
     +    *u1**(-2)*u2 + 64*h_r*st*cb*mg*t1*u1**(-1) + 64*h_r*st*cb*mg*
     +    t1**2*u1**(-2) )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_r*st*cb*mg*u1**(-1)*u2 - 32*h_r*ct*cb*m1*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 32*h_r*ct*cb*m1*
     +    m2**2*msb2**2*s**(-1)*t1*u1**(-1)*u2**(-1) - 64*h_r*ct*cb*m1*
     +    mg**2*s**(-1)*t1*u1**(-1) - 64*h_r*ct*cb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 96*h_r*ct*cb*m1*mg**2*s*u1**(-2) - 96*h_r*ct*cb
     +    *m1*mg**2*t1*u1**(-2) - 96*h_r*ct*cb*m1*mg**2*u1**(-2)*u2 - 
     +    96*h_r*ct*cb*m1*mg**2*u1**(-1) + 64*h_r*ct*cb*m1*msb2**2*
     +    s**(-1)*t1*u1**(-1) + 64*h_r*ct*cb*m1*msb2**2*s**(-1)*
     +    u1**(-1)*u2 + 96*h_r*ct*cb*m1*msb2**2*s*u1**(-2) + 96*h_r*ct*
     +    cb*m1*msb2**2*t1*u1**(-2) + 96*h_r*ct*cb*m1*msb2**2*u1**(-2)*
     +    u2 + 96*h_r*ct*cb*m1*msb2**2*u1**(-1) - 32*h_r*ct*cb*m1**3*
     +    mg**2*s**(-1)*u1**(-1) - 64*h_r*ct*cb*m1**3*mg**2*u1**(-2) + 
     +    32*h_r*ct*cb*m1**3*msb2**2*s**(-1)*u1**(-1) + 64*h_r*ct*cb*
     +    m1**3*msb2**2*u1**(-2) )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) - 32*h_l*st*sb*m1*m2**2*msb2**2*s**(-1)*t1*u1**(-1)*
     +    u2**(-1) + 64*h_l*st*sb*m1*mg**2*s**(-1)*t1*u1**(-1) + 64*h_l
     +    *st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 96*h_l*st*sb*m1*mg**2*s
     +    *u1**(-2) + 96*h_l*st*sb*m1*mg**2*t1*u1**(-2) + 96*h_l*st*sb*
     +    m1*mg**2*u1**(-2)*u2 + 96*h_l*st*sb*m1*mg**2*u1**(-1) - 64*
     +    h_l*st*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) - 64*h_l*st*sb*m1*
     +    msb2**2*s**(-1)*u1**(-1)*u2 - 96*h_l*st*sb*m1*msb2**2*s*
     +    u1**(-2) - 96*h_l*st*sb*m1*msb2**2*t1*u1**(-2) - 96*h_l*st*sb
     +    *m1*msb2**2*u1**(-2)*u2 - 96*h_l*st*sb*m1*msb2**2*u1**(-1) + 
     +    32*h_l*st*sb*m1**3*mg**2*s**(-1)*u1**(-1) + 64*h_l*st*sb*
     +    m1**3*mg**2*u1**(-2) - 32*h_l*st*sb*m1**3*msb2**2*s**(-1)*
     +    u1**(-1) - 64*h_l*st*sb*m1**3*msb2**2*u1**(-2) + 32*h_l*ct*sb
     +    *m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_l*ct*sb*m1**2*mg*s**(-1)
     +    *u1**(-1)*u2 )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 96*h_l*ct*sb*m1**2*mg*u1**(-2)*u2 - 32*h_l*ct*sb*
     +    m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*m2**2*mg*u1**(-2)
     +    *u2 - 32*h_l*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*sb*
     +    mg*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*mg*s*t1*u1**(-2) + 
     +    32*h_l*ct*sb*mg*s*u1**(-2)*u2 - 32*h_l*ct*sb*mg*t1*u1**(-2)*
     +    u2 - 64*h_l*ct*sb*mg*t1*u1**(-1) - 64*h_l*ct*sb*mg*t1**2*
     +    u1**(-2) + 32*h_l*ct*sb*mg*u1**(-1)*u2 + 32*h_r*st*cb*m1**2*
     +    mg*s**(-1)*t1*u1**(-1) + 64*h_r*st*cb*m1**2*mg*s**(-1)*
     +    u1**(-1)*u2 + 96*h_r*st*cb*m1**2*mg*u1**(-2)*u2 - 32*h_r*st*
     +    cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*st*cb*m2**2*mg*
     +    u1**(-2)*u2 - 32*h_r*st*cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r
     +    *st*cb*mg*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*cb*mg*s*t1*
     +    u1**(-2) + 32*h_r*st*cb*mg*s*u1**(-2)*u2 - 32*h_r*st*cb*mg*t1
     +    *u1**(-2)*u2 - 64*h_r*st*cb*mg*t1*u1**(-1) - 64*h_r*st*cb*mg*
     +    t1**2*u1**(-2) )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*st*cb*mg*u1**(-1)*u2 + 32*h_r*ct*cb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2**(-1) - 32*h_r*ct*cb*m1*m2**2*
     +    msb2**2*s**(-1)*t1*u1**(-1)*u2**(-1) + 64*h_r*ct*cb*m1*mg**2*
     +    s**(-1)*t1*u1**(-1) + 64*h_r*ct*cb*m1*mg**2*s**(-1)*u1**(-1)*
     +    u2 + 96*h_r*ct*cb*m1*mg**2*s*u1**(-2) + 96*h_r*ct*cb*m1*mg**2
     +    *t1*u1**(-2) + 96*h_r*ct*cb*m1*mg**2*u1**(-2)*u2 + 96*h_r*ct*
     +    cb*m1*mg**2*u1**(-1) - 64*h_r*ct*cb*m1*msb2**2*s**(-1)*t1*
     +    u1**(-1) - 64*h_r*ct*cb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 96*
     +    h_r*ct*cb*m1*msb2**2*s*u1**(-2) - 96*h_r*ct*cb*m1*msb2**2*t1*
     +    u1**(-2) - 96*h_r*ct*cb*m1*msb2**2*u1**(-2)*u2 - 96*h_r*ct*cb
     +    *m1*msb2**2*u1**(-1) + 32*h_r*ct*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1) + 64*h_r*ct*cb*m1**3*mg**2*u1**(-2) - 32*h_r*ct*cb*
     +    m1**3*msb2**2*s**(-1)*u1**(-1) - 64*h_r*ct*cb*m1**3*msb2**2*
     +    u1**(-2) )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*mg**2*s**(-1) - 32*h_l*st*sb*m1*
     +    mg**2*u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s**(-1) + 32*h_l*st*
     +    sb*m1*msb2**2*u1**(-1) - 16*h_l*st*sb*m1*s*u1**(-1) - 16*h_l*
     +    st*sb*m1*t1*u1**(-1) + 16*h_l*ct*sb*mg*s**(-1)*t1 + 16*h_l*ct
     +    *sb*mg*s*u1**(-1) + 16*h_l*ct*sb*mg*t1*u1**(-1) + 16*h_l*ct*
     +    sb*mg + 16*h_r*st*cb*mg*s**(-1)*t1 + 16*h_r*st*cb*mg*s*
     +    u1**(-1) + 16*h_r*st*cb*mg*t1*u1**(-1) + 16*h_r*st*cb*mg - 32
     +    *h_r*ct*cb*m1*mg**2*s**(-1) - 32*h_r*ct*cb*m1*mg**2*u1**(-1)
     +     + 32*h_r*ct*cb*m1*msb2**2*s**(-1) + 32*h_r*ct*cb*m1*msb2**2*
     +    u1**(-1) - 16*h_r*ct*cb*m1*s*u1**(-1) - 16*h_r*ct*cb*m1*t1*
     +    u1**(-1) )
      MM_s = MM_s + SCC(4,5)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*mg**2*s**(-1) + 32*h_l*st*sb*m1*mg**2
     +    *u1**(-1) - 32*h_l*st*sb*m1*msb2**2*s**(-1) - 32*h_l*st*sb*m1
     +    *msb2**2*u1**(-1) + 16*h_l*st*sb*m1*s**(-1)*t1 + 16*h_l*st*sb
     +    *m1*s*u1**(-1) + 16*h_l*st*sb*m1*t1*u1**(-1) + 16*h_l*st*sb*
     +    m1 + 16*h_l*ct*sb*mg*s**(-1)*t1 + 16*h_l*ct*sb*mg*s*u1**(-1)
     +     + 16*h_l*ct*sb*mg*t1*u1**(-1) + 16*h_l*ct*sb*mg + 16*h_r*st*
     +    cb*mg*s**(-1)*t1 + 16*h_r*st*cb*mg*s*u1**(-1) + 16*h_r*st*cb*
     +    mg*t1*u1**(-1) + 16*h_r*st*cb*mg + 32*h_r*ct*cb*m1*mg**2*
     +    s**(-1) + 32*h_r*ct*cb*m1*mg**2*u1**(-1) - 32*h_r*ct*cb*m1*
     +    msb2**2*s**(-1) - 32*h_r*ct*cb*m1*msb2**2*u1**(-1) + 16*h_r*
     +    ct*cb*m1*s**(-1)*t1 + 16*h_r*ct*cb*m1*s*u1**(-1) + 16*h_r*ct*
     +    cb*m1*t1*u1**(-1) + 16*h_r*ct*cb*m1 )
      MM_s = MM_s + SCC(5,3)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*mg*s**(-1)*t1 + 16*h_l*st*cb*mg*s**(-1)*
     +    u1 + 16*h_l*st*cb*mg*s*u1**(-1) + 32*h_l*st*cb*mg*t1*u1**(-1)
     +     + 32*h_l*st*cb*mg - 32*h_l*ct*cb*m1*mg**2*s**(-1) - 32*h_l*
     +    ct*cb*m1*mg**2*u1**(-1) + 32*h_l*ct*cb*m1*msb1**2*s**(-1) + 
     +    32*h_l*ct*cb*m1*msb1**2*s*t2**(-2) + 32*h_l*ct*cb*m1*msb1**2*
     +    s**2*u1**(-1)*t2**(-2) - 32*h_l*ct*cb*m1*msb1**2*t2**(-1) - 
     +    32*h_l*ct*cb*m1*mst1**2*s*t2**(-2) - 32*h_l*ct*cb*m1*mst1**2*
     +    s**2*u1**(-1)*t2**(-2) + 32*h_l*ct*cb*m1*mst1**2*u1**(-1) + 
     +    32*h_l*ct*cb*m1*mst1**2*t2**(-1) + 16*h_l*ct*cb*m1*s**(-1)*u1
     +     - 32*h_l*ct*cb*m1*s*t1*t2**(-2) + 48*h_l*ct*cb*m1*s*u1**(-1)
     +     + 32*h_l*ct*cb*m1*s*t2**(-1) - 32*h_l*ct*cb*m1*s**2*t1*
     +    u1**(-1)*t2**(-2) - 32*h_l*ct*cb*m1*s**2*t2**(-2) - 32*h_l*ct
     +    *cb*m1*s**3*u1**(-1)*t2**(-2) + 32*h_l*ct*cb*m1*t1*u1**(-1)
     +     + 32*h_l*ct*cb*m1*t1*t2**(-1) + 32*h_l*ct*cb*m1 - 32*h_l*ct*
     +    cb*m1**3*s*t2**(-2) )
      MM_s = MM_s + SCC(5,3)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1**3*s**2*u1**(-1)*t2**(-2) + 32*h_l
     +    *ct*cb*m1**3*u1**(-1) + 32*h_l*ct*cb*m1**3*t2**(-1) - 32*h_r*
     +    st*sb*m1*mg**2*s**(-1) - 32*h_r*st*sb*m1*mg**2*u1**(-1) + 32*
     +    h_r*st*sb*m1*msb1**2*s**(-1) + 32*h_r*st*sb*m1*msb1**2*s*
     +    t2**(-2) + 32*h_r*st*sb*m1*msb1**2*s**2*u1**(-1)*t2**(-2) - 
     +    32*h_r*st*sb*m1*msb1**2*t2**(-1) - 32*h_r*st*sb*m1*mst1**2*s*
     +    t2**(-2) - 32*h_r*st*sb*m1*mst1**2*s**2*u1**(-1)*t2**(-2) + 
     +    32*h_r*st*sb*m1*mst1**2*u1**(-1) + 32*h_r*st*sb*m1*mst1**2*
     +    t2**(-1) + 16*h_r*st*sb*m1*s**(-1)*u1 - 32*h_r*st*sb*m1*s*t1*
     +    t2**(-2) + 48*h_r*st*sb*m1*s*u1**(-1) + 32*h_r*st*sb*m1*s*
     +    t2**(-1) - 32*h_r*st*sb*m1*s**2*t1*u1**(-1)*t2**(-2) - 32*h_r
     +    *st*sb*m1*s**2*t2**(-2) - 32*h_r*st*sb*m1*s**3*u1**(-1)*
     +    t2**(-2) + 32*h_r*st*sb*m1*t1*u1**(-1) + 32*h_r*st*sb*m1*t1*
     +    t2**(-1) + 32*h_r*st*sb*m1 - 32*h_r*st*sb*m1**3*s*t2**(-2) - 
     +    32*h_r*st*sb*m1**3*s**2*u1**(-1)*t2**(-2) )
      MM_s = MM_s + SCC(5,3)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*st*sb*m1**3*u1**(-1) + 32*h_r*st*sb*m1**3*
     +    t2**(-1) + 32*h_r*ct*sb*mg*s**(-1)*t1 + 16*h_r*ct*sb*mg*
     +    s**(-1)*u1 + 16*h_r*ct*sb*mg*s*u1**(-1) + 32*h_r*ct*sb*mg*t1*
     +    u1**(-1) + 32*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(5,4)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*mg*s**(-1)*t1 - 16*h_l*st*sb*mg*
     +    s**(-1)*u1 - 16*h_l*st*sb*mg*s*u1**(-1) - 32*h_l*st*sb*mg*t1*
     +    u1**(-1) - 32*h_l*st*sb*mg + 32*h_l*ct*sb*m1*mg**2*s**(-1) + 
     +    32*h_l*ct*sb*m1*mg**2*u1**(-1) - 32*h_l*ct*sb*m1*msb2**2*
     +    s**(-1) - 32*h_l*ct*sb*m1*msb2**2*s*t2**(-2) - 32*h_l*ct*sb*
     +    m1*msb2**2*s**2*u1**(-1)*t2**(-2) + 32*h_l*ct*sb*m1*msb2**2*
     +    t2**(-1) + 32*h_l*ct*sb*m1*mst1**2*s*t2**(-2) + 32*h_l*ct*sb*
     +    m1*mst1**2*s**2*u1**(-1)*t2**(-2) - 32*h_l*ct*sb*m1*mst1**2*
     +    u1**(-1) - 32*h_l*ct*sb*m1*mst1**2*t2**(-1) - 16*h_l*ct*sb*m1
     +    *s**(-1)*u1 + 32*h_l*ct*sb*m1*s*t1*t2**(-2) - 48*h_l*ct*sb*m1
     +    *s*u1**(-1) - 32*h_l*ct*sb*m1*s*t2**(-1) + 32*h_l*ct*sb*m1*
     +    s**2*t1*u1**(-1)*t2**(-2) + 32*h_l*ct*sb*m1*s**2*t2**(-2) + 
     +    32*h_l*ct*sb*m1*s**3*u1**(-1)*t2**(-2) - 32*h_l*ct*sb*m1*t1*
     +    u1**(-1) - 32*h_l*ct*sb*m1*t1*t2**(-1) - 32*h_l*ct*sb*m1 + 32
     +    *h_l*ct*sb*m1**3*s*t2**(-2) )
      MM_s = MM_s + SCC(5,4)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1**3*s**2*u1**(-1)*t2**(-2) - 32*h_l*ct
     +    *sb*m1**3*u1**(-1) - 32*h_l*ct*sb*m1**3*t2**(-1) - 32*h_r*st*
     +    cb*m1*mg**2*s**(-1) - 32*h_r*st*cb*m1*mg**2*u1**(-1) + 32*h_r
     +    *st*cb*m1*msb2**2*s**(-1) + 32*h_r*st*cb*m1*msb2**2*s*
     +    t2**(-2) + 32*h_r*st*cb*m1*msb2**2*s**2*u1**(-1)*t2**(-2) - 
     +    32*h_r*st*cb*m1*msb2**2*t2**(-1) - 32*h_r*st*cb*m1*mst1**2*s*
     +    t2**(-2) - 32*h_r*st*cb*m1*mst1**2*s**2*u1**(-1)*t2**(-2) + 
     +    32*h_r*st*cb*m1*mst1**2*u1**(-1) + 32*h_r*st*cb*m1*mst1**2*
     +    t2**(-1) + 16*h_r*st*cb*m1*s**(-1)*u1 - 32*h_r*st*cb*m1*s*t1*
     +    t2**(-2) + 48*h_r*st*cb*m1*s*u1**(-1) + 32*h_r*st*cb*m1*s*
     +    t2**(-1) - 32*h_r*st*cb*m1*s**2*t1*u1**(-1)*t2**(-2) - 32*h_r
     +    *st*cb*m1*s**2*t2**(-2) - 32*h_r*st*cb*m1*s**3*u1**(-1)*
     +    t2**(-2) + 32*h_r*st*cb*m1*t1*u1**(-1) + 32*h_r*st*cb*m1*t1*
     +    t2**(-1) + 32*h_r*st*cb*m1 - 32*h_r*st*cb*m1**3*s*t2**(-2) - 
     +    32*h_r*st*cb*m1**3*s**2*u1**(-1)*t2**(-2) )
      MM_s = MM_s + SCC(5,4)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*st*cb*m1**3*u1**(-1) + 32*h_r*st*cb*m1**3*
     +    t2**(-1) + 32*h_r*ct*cb*mg*s**(-1)*t1 + 16*h_r*ct*cb*mg*
     +    s**(-1)*u1 + 16*h_r*ct*cb*mg*s*u1**(-1) + 32*h_r*ct*cb*mg*t1*
     +    u1**(-1) + 32*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(5,5)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*mg**2*s**(-1) + 32*h_l*st*cb*m1*mg**2
     +    *u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s**(-1) - 32*h_l*st*cb*m1
     +    *msb1**2*s*t2**(-2) - 32*h_l*st*cb*m1*msb1**2*s**2*u1**(-1)*
     +    t2**(-2) + 32*h_l*st*cb*m1*msb1**2*t2**(-1) + 32*h_l*st*cb*m1
     +    *mst2**2*s*t2**(-2) + 32*h_l*st*cb*m1*mst2**2*s**2*u1**(-1)*
     +    t2**(-2) - 32*h_l*st*cb*m1*mst2**2*u1**(-1) - 32*h_l*st*cb*m1
     +    *mst2**2*t2**(-1) - 16*h_l*st*cb*m1*s**(-1)*u1 + 32*h_l*st*cb
     +    *m1*s*t1*t2**(-2) - 48*h_l*st*cb*m1*s*u1**(-1) - 32*h_l*st*cb
     +    *m1*s*t2**(-1) + 32*h_l*st*cb*m1*s**2*t1*u1**(-1)*t2**(-2) + 
     +    32*h_l*st*cb*m1*s**2*t2**(-2) + 32*h_l*st*cb*m1*s**3*u1**(-1)
     +    *t2**(-2) - 32*h_l*st*cb*m1*t1*u1**(-1) - 32*h_l*st*cb*m1*t1*
     +    t2**(-1) - 32*h_l*st*cb*m1 + 32*h_l*st*cb*m1**3*s*t2**(-2) + 
     +    32*h_l*st*cb*m1**3*s**2*u1**(-1)*t2**(-2) - 32*h_l*st*cb*
     +    m1**3*u1**(-1) - 32*h_l*st*cb*m1**3*t2**(-1) + 32*h_l*ct*cb*
     +    mg*s**(-1)*t1 )
      MM_s = MM_s + SCC(5,5)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*ct*cb*mg*s**(-1)*u1 + 16*h_l*ct*cb*mg*s*
     +    u1**(-1) + 32*h_l*ct*cb*mg*t1*u1**(-1) + 32*h_l*ct*cb*mg - 32
     +    *h_r*st*sb*mg*s**(-1)*t1 - 16*h_r*st*sb*mg*s**(-1)*u1 - 16*
     +    h_r*st*sb*mg*s*u1**(-1) - 32*h_r*st*sb*mg*t1*u1**(-1) - 32*
     +    h_r*st*sb*mg - 32*h_r*ct*sb*m1*mg**2*s**(-1) - 32*h_r*ct*sb*
     +    m1*mg**2*u1**(-1) + 32*h_r*ct*sb*m1*msb1**2*s**(-1) + 32*h_r*
     +    ct*sb*m1*msb1**2*s*t2**(-2) + 32*h_r*ct*sb*m1*msb1**2*s**2*
     +    u1**(-1)*t2**(-2) - 32*h_r*ct*sb*m1*msb1**2*t2**(-1) - 32*h_r
     +    *ct*sb*m1*mst2**2*s*t2**(-2) - 32*h_r*ct*sb*m1*mst2**2*s**2*
     +    u1**(-1)*t2**(-2) + 32*h_r*ct*sb*m1*mst2**2*u1**(-1) + 32*h_r
     +    *ct*sb*m1*mst2**2*t2**(-1) + 16*h_r*ct*sb*m1*s**(-1)*u1 - 32*
     +    h_r*ct*sb*m1*s*t1*t2**(-2) + 48*h_r*ct*sb*m1*s*u1**(-1) + 32*
     +    h_r*ct*sb*m1*s*t2**(-1) - 32*h_r*ct*sb*m1*s**2*t1*u1**(-1)*
     +    t2**(-2) - 32*h_r*ct*sb*m1*s**2*t2**(-2) - 32*h_r*ct*sb*m1*
     +    s**3*u1**(-1)*t2**(-2) )
      MM_s = MM_s + SCC(5,5)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*ct*sb*m1*t1*u1**(-1) + 32*h_r*ct*sb*m1*t1*
     +    t2**(-1) + 32*h_r*ct*sb*m1 - 32*h_r*ct*sb*m1**3*s*t2**(-2) - 
     +    32*h_r*ct*sb*m1**3*s**2*u1**(-1)*t2**(-2) + 32*h_r*ct*sb*
     +    m1**3*u1**(-1) + 32*h_r*ct*sb*m1**3*t2**(-1) )
      MM_s = MM_s + SCC(5,6)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*mg**2*s**(-1) - 32*h_l*st*sb*m1*
     +    mg**2*u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s**(-1) + 32*h_l*st*
     +    sb*m1*msb2**2*s*t2**(-2) + 32*h_l*st*sb*m1*msb2**2*s**2*
     +    u1**(-1)*t2**(-2) - 32*h_l*st*sb*m1*msb2**2*t2**(-1) - 32*h_l
     +    *st*sb*m1*mst2**2*s*t2**(-2) - 32*h_l*st*sb*m1*mst2**2*s**2*
     +    u1**(-1)*t2**(-2) + 32*h_l*st*sb*m1*mst2**2*u1**(-1) + 32*h_l
     +    *st*sb*m1*mst2**2*t2**(-1) + 16*h_l*st*sb*m1*s**(-1)*u1 - 32*
     +    h_l*st*sb*m1*s*t1*t2**(-2) + 48*h_l*st*sb*m1*s*u1**(-1) + 32*
     +    h_l*st*sb*m1*s*t2**(-1) - 32*h_l*st*sb*m1*s**2*t1*u1**(-1)*
     +    t2**(-2) - 32*h_l*st*sb*m1*s**2*t2**(-2) - 32*h_l*st*sb*m1*
     +    s**3*u1**(-1)*t2**(-2) + 32*h_l*st*sb*m1*t1*u1**(-1) + 32*h_l
     +    *st*sb*m1*t1*t2**(-1) + 32*h_l*st*sb*m1 - 32*h_l*st*sb*m1**3*
     +    s*t2**(-2) - 32*h_l*st*sb*m1**3*s**2*u1**(-1)*t2**(-2) + 32*
     +    h_l*st*sb*m1**3*u1**(-1) + 32*h_l*st*sb*m1**3*t2**(-1) - 32*
     +    h_l*ct*sb*mg*s**(-1)*t1 )
      MM_s = MM_s + SCC(5,6)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*ct*sb*mg*s**(-1)*u1 - 16*h_l*ct*sb*mg*s*
     +    u1**(-1) - 32*h_l*ct*sb*mg*t1*u1**(-1) - 32*h_l*ct*sb*mg - 32
     +    *h_r*st*cb*mg*s**(-1)*t1 - 16*h_r*st*cb*mg*s**(-1)*u1 - 16*
     +    h_r*st*cb*mg*s*u1**(-1) - 32*h_r*st*cb*mg*t1*u1**(-1) - 32*
     +    h_r*st*cb*mg - 32*h_r*ct*cb*m1*mg**2*s**(-1) - 32*h_r*ct*cb*
     +    m1*mg**2*u1**(-1) + 32*h_r*ct*cb*m1*msb2**2*s**(-1) + 32*h_r*
     +    ct*cb*m1*msb2**2*s*t2**(-2) + 32*h_r*ct*cb*m1*msb2**2*s**2*
     +    u1**(-1)*t2**(-2) - 32*h_r*ct*cb*m1*msb2**2*t2**(-1) - 32*h_r
     +    *ct*cb*m1*mst2**2*s*t2**(-2) - 32*h_r*ct*cb*m1*mst2**2*s**2*
     +    u1**(-1)*t2**(-2) + 32*h_r*ct*cb*m1*mst2**2*u1**(-1) + 32*h_r
     +    *ct*cb*m1*mst2**2*t2**(-1) + 16*h_r*ct*cb*m1*s**(-1)*u1 - 32*
     +    h_r*ct*cb*m1*s*t1*t2**(-2) + 48*h_r*ct*cb*m1*s*u1**(-1) + 32*
     +    h_r*ct*cb*m1*s*t2**(-1) - 32*h_r*ct*cb*m1*s**2*t1*u1**(-1)*
     +    t2**(-2) - 32*h_r*ct*cb*m1*s**2*t2**(-2) - 32*h_r*ct*cb*m1*
     +    s**3*u1**(-1)*t2**(-2) )
      MM_s = MM_s + SCC(5,6)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_r*ct*cb*m1*t1*u1**(-1) + 32*h_r*ct*cb*m1*t1*
     +    t2**(-1) + 32*h_r*ct*cb*m1 - 32*h_r*ct*cb*m1**3*s*t2**(-2) - 
     +    32*h_r*ct*cb*m1**3*s**2*u1**(-1)*t2**(-2) + 32*h_r*ct*cb*
     +    m1**3*u1**(-1) + 32*h_r*ct*cb*m1**3*t2**(-1) )
      MM_s = MM_s + SCC(5,7)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*mg*s**(-1)*t1 - 16*h_l*st*cb*mg*
     +    s**(-1)*u1 - 16*h_l*st*cb*mg*s*u1**(-1) - 32*h_l*st*cb*mg*t1*
     +    u1**(-1) - 32*h_l*st*cb*mg - 32*h_l*ct*cb*m1*mg**2*s**(-1) - 
     +    32*h_l*ct*cb*m1*mg**2*u1**(-1) + 32*h_l*ct*cb*m1*msb1**2*
     +    s**(-1) + 32*h_l*ct*cb*m1*msb1**2*u1**(-1) - 32*h_r*st*sb*m1*
     +    mg**2*s**(-1) - 32*h_r*st*sb*m1*mg**2*u1**(-1) + 32*h_r*st*sb
     +    *m1*msb1**2*s**(-1) + 32*h_r*st*sb*m1*msb1**2*u1**(-1) - 32*
     +    h_r*ct*sb*mg*s**(-1)*t1 - 16*h_r*ct*sb*mg*s**(-1)*u1 - 16*h_r
     +    *ct*sb*mg*s*u1**(-1) - 32*h_r*ct*sb*mg*t1*u1**(-1) - 32*h_r*
     +    ct*sb*mg )
      MM_s = MM_s + SCC(5,8)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*cb*m1*mg**2*s**(-1) + 32*h_l*st*cb*m1*mg**2
     +    *u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s**(-1) - 32*h_l*st*cb*m1
     +    *msb1**2*u1**(-1) - 32*h_l*ct*cb*mg*s**(-1)*t1 - 16*h_l*ct*cb
     +    *mg*s**(-1)*u1 - 16*h_l*ct*cb*mg*s*u1**(-1) - 32*h_l*ct*cb*mg
     +    *t1*u1**(-1) - 32*h_l*ct*cb*mg + 32*h_r*st*sb*mg*s**(-1)*t1
     +     + 16*h_r*st*sb*mg*s**(-1)*u1 + 16*h_r*st*sb*mg*s*u1**(-1) + 
     +    32*h_r*st*sb*mg*t1*u1**(-1) + 32*h_r*st*sb*mg - 32*h_r*ct*sb*
     +    m1*mg**2*s**(-1) - 32*h_r*ct*sb*m1*mg**2*u1**(-1) + 32*h_r*ct
     +    *sb*m1*msb1**2*s**(-1) + 32*h_r*ct*sb*m1*msb1**2*u1**(-1) )
      MM_s = MM_s + SCC(5,9)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*mg*s**(-1)*t1 + 16*h_l*st*sb*mg*s**(-1)*
     +    u1 + 16*h_l*st*sb*mg*s*u1**(-1) + 32*h_l*st*sb*mg*t1*u1**(-1)
     +     + 32*h_l*st*sb*mg + 32*h_l*ct*sb*m1*mg**2*s**(-1) + 32*h_l*
     +    ct*sb*m1*mg**2*u1**(-1) - 32*h_l*ct*sb*m1*msb2**2*s**(-1) - 
     +    32*h_l*ct*sb*m1*msb2**2*u1**(-1) - 32*h_r*st*cb*m1*mg**2*
     +    s**(-1) - 32*h_r*st*cb*m1*mg**2*u1**(-1) + 32*h_r*st*cb*m1*
     +    msb2**2*s**(-1) + 32*h_r*st*cb*m1*msb2**2*u1**(-1) - 32*h_r*
     +    ct*cb*mg*s**(-1)*t1 - 16*h_r*ct*cb*mg*s**(-1)*u1 - 16*h_r*ct*
     +    cb*mg*s*u1**(-1) - 32*h_r*ct*cb*mg*t1*u1**(-1) - 32*h_r*ct*cb
     +    *mg )
      MM_s = MM_s + SCC(5,10)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*sb*m1*mg**2*s**(-1) - 32*h_l*st*sb*m1*
     +    mg**2*u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s**(-1) + 32*h_l*st*
     +    sb*m1*msb2**2*u1**(-1) + 32*h_l*ct*sb*mg*s**(-1)*t1 + 16*h_l*
     +    ct*sb*mg*s**(-1)*u1 + 16*h_l*ct*sb*mg*s*u1**(-1) + 32*h_l*ct*
     +    sb*mg*t1*u1**(-1) + 32*h_l*ct*sb*mg + 32*h_r*st*cb*mg*s**(-1)
     +    *t1 + 16*h_r*st*cb*mg*s**(-1)*u1 + 16*h_r*st*cb*mg*s*u1**(-1)
     +     + 32*h_r*st*cb*mg*t1*u1**(-1) + 32*h_r*st*cb*mg - 32*h_r*ct*
     +    cb*m1*mg**2*s**(-1) - 32*h_r*ct*cb*m1*mg**2*u1**(-1) + 32*h_r
     +    *ct*cb*m1*msb2**2*s**(-1) + 32*h_r*ct*cb*m1*msb2**2*u1**(-1)
     +     )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 96*h_l*ct*cb*m1*m2**2*mg**2*
     +    s**(-2)*t1*u2 - 64*h_l*ct*cb*m1*m2**2*mg**2*s**(-2)*t1**2 - 
     +    32*h_l*ct*cb*m1*m2**2*mg**2*s**(-2)*u2**2 - 192*h_l*ct*cb*m1*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_l*ct*cb*m1*m2**2*
     +    mg**2*s**(-1)*t1 - 128*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*t1**2
     +    *u1**(-1) - 64*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_l*ct*
     +    cb*m1*m2**2*mg**2*t1*u1**(-1) - 96*h_l*ct*cb*m1*m2**2*mg**2*
     +    u1**(-1)*u2 + 32*h_l*ct*cb*m1*m2**2*msb1**2*s**(-2)*t1*u2 + 
     +    32*h_l*ct*cb*m1*m2**2*msb1**2*s**(-2)*t1**2 + 64*h_l*ct*cb*m1
     +    *m2**2*msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*m2**2
     +    *msb1**2*s**(-1)*t1 + 64*h_l*ct*cb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_l*ct*cb*m1*m2**2*msb1**2*t1*u1**(-1) + 
     +    64*h_l*ct*cb*m1*m2**2*mst1**2*s**(-2)*t1*u2 + 32*h_l*ct*cb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*ct*cb*m1*m2**2*mst1**2*s**(-2)
     +    *u2**2 + 128*h_l*ct*cb*m1*m2**2*mst1**2*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_l*ct*cb*m1*m2**2*mst1**2*s**(-1)*t1 + 64*h_l*ct*cb*
     +    m1*m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*
     +    m2**2*mst1**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*cb*m1*m2**2*
     +    mst1**2*s**(-1)*u2 + 64*h_l*ct*cb*m1*m2**2*mst1**2*t1*
     +    u1**(-1) + 96*h_l*ct*cb*m1*m2**2*mst1**2*u1**(-1)*u2 + 96*h_l
     +    *ct*cb*m1*m2**2*s**(-1)*t1*u2 + 32*h_l*ct*cb*m1*m2**2*s**(-1)
     +    *t1**2 + 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2**2 + 80*h_l*ct*cb*
     +    m1*m2**2*s*t1*u1**(-1) + 176*h_l*ct*cb*m1*m2**2*s*u1**(-1)*u2
     +     + 192*h_l*ct*cb*m1*m2**2*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*
     +    m2**2*t1 + 64*h_l*ct*cb*m1*m2**2*t1**2*u1**(-1) + 128*h_l*ct*
     +    cb*m1*m2**2*u1**(-1)*u2**2 + 128*h_l*ct*cb*m1*m2**2*u2 + 32*
     +    h_l*ct*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_l*ct*cb*m1*m2**4*msb1**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*m1*m2**4*mst1**2*s**(-1)*
     +    t1*u1**(-1) - 16*h_l*ct*cb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_l*ct*cb*m1*m2**4*s**(-2)*t1*u2 - 32*h_l*ct*cb*m1*
     +    m2**4*s**(-2)*t1**2 - 32*h_l*ct*cb*m1*m2**4*s**(-2)*u2**2 - 
     +    128*h_l*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*
     +    m1*m2**4*s**(-1)*t1 - 64*h_l*ct*cb*m1*m2**4*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_l*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 - 64*
     +    h_l*ct*cb*m1*m2**4*s**(-1)*u2 - 96*h_l*ct*cb*m1*m2**4*t1*
     +    u1**(-1) - 112*h_l*ct*cb*m1*m2**4*u1**(-1)*u2 + 16*h_l*ct*cb*
     +    m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*m1*m2**6*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*ct*cb*m1*mg**2*s**(-1)*t1*u2 + 32*h_l*ct
     +    *cb*m1*mg**2*s**(-1)*u2**2 + 80*h_l*ct*cb*m1*mg**2*s*u1**(-1)
     +    *u2 + 64*h_l*ct*cb*m1*mg**2*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*
     +    mg**2*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*cb*m1*mg**2*u2 - 32*h_l*ct*
     +    cb*m1*msb1**2*s**(-1)*t1*u2 - 32*h_l*ct*cb*m1*msb1**2*s**(-1)
     +    *t1**2 - 80*h_l*ct*cb*m1*msb1**2*s*t1*u1**(-1) - 64*h_l*ct*cb
     +    *m1*msb1**2*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*msb1**2*t1 - 64*
     +    h_l*ct*cb*m1*msb1**2*t1**2*u1**(-1) + 32*h_l*ct*cb*m1*mst1**2
     +    *s**(-1)*t1**2 - 32*h_l*ct*cb*m1*mst1**2*s**(-1)*u2**2 + 80*
     +    h_l*ct*cb*m1*mst1**2*s*t1*u1**(-1) - 80*h_l*ct*cb*m1*mst1**2*
     +    s*u1**(-1)*u2 + 64*h_l*ct*cb*m1*mst1**2*t1 + 64*h_l*ct*cb*m1*
     +    mst1**2*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*mst1**2*u1**(-1)*
     +    u2**2 - 64*h_l*ct*cb*m1*mst1**2*u2 - 64*h_l*ct*cb*m1*s*t1*
     +    u1**(-1)*u2 - 64*h_l*ct*cb*m1*s*u1**(-1)*u2**2 - 64*h_l*ct*cb
     +    *m1*s*u2 - 80*h_l*ct*cb*m1*s**2*u1**(-1)*u2 - 32*h_l*ct*cb*m1
     +    *t1*u2 - 32*h_l*ct*cb*m1*u2**2 - 32*h_l*ct*cb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) - 32*h_l*ct*cb*m1**3*m2**2*msb1**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*ct*cb*m1**3*m2**2*mst1**2*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*m1**3*m2**2*mst1**2*
     +    s**(-1)*u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*ct*cb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_l*ct*cb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_l*ct*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_l*ct*cb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*cb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_l*ct*cb*m1**3*m2**2*t1*u1**(-1) + 96*h_l*
     +    ct*cb*m1**3*m2**2*u1**(-1)*u2 - 32*h_l*ct*cb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*ct*cb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*ct*cb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_l
     +    *ct*cb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_l*ct*cb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1**3*mg**2*s**(-1)*
     +    u2 - 64*h_l*ct*cb*m1**3*mg**2*u1**(-1)*u2 + 96*h_l*ct*cb*
     +    m1**3*msb1**2*s**(-2)*t1*u2 + 32*h_l*ct*cb*m1**3*msb1**2*
     +    s**(-2)*t1**2 + 64*h_l*ct*cb*m1**3*msb1**2*s**(-2)*u2**2 + 
     +    192*h_l*ct*cb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*
     +    ct*cb*m1**3*msb1**2*s**(-1)*t1 + 64*h_l*ct*cb*m1**3*msb1**2*
     +    s**(-1)*t1**2*u1**(-1) + 128*h_l*ct*cb*m1**3*msb1**2*s**(-1)*
     +    u1**(-1)*u2**2 + 128*h_l*ct*cb*m1**3*msb1**2*s**(-1)*u2 + 64*
     +    h_l*ct*cb*m1**3*msb1**2*t1*u1**(-1) + 160*h_l*ct*cb*m1**3*
     +    msb1**2*u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*mst1**2*s**(-2)*t1*
     +    u2 - 32*h_l*ct*cb*m1**3*mst1**2*s**(-2)*t1**2 - 32*h_l*ct*cb*
     +    m1**3*mst1**2*s**(-2)*u2**2 - 128*h_l*ct*cb*m1**3*mst1**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)*
     +    t1 - 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)*t1**2*u1**(-1) - 64*
     +    h_l*ct*cb*m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)
     +    *u2 - 64*h_l*ct*cb*m1**3*mst1**2*t1*u1**(-1) - 96*h_l*ct*cb*
     +    m1**3*mst1**2*u1**(-1)*u2 + 32*h_l*ct*cb*m1**3*s**(-1)*t1*u2
     +     + 32*h_l*ct*cb*m1**3*s**(-1)*u2**2 + 64*h_l*ct*cb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*t1*u1**(-1)*u2 + 64*h_l*ct*
     +    cb*m1**3*u1**(-1)*u2**2 + 64*h_l*ct*cb*m1**3*u2 + 16*h_l*ct*
     +    cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_l*ct*cb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_l*ct*cb*m1**5*msb1**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_l*ct*cb*m1**5*msb1**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct
     +    *cb*m1**5*mst1**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*m1**5*
     +    mst1**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct*cb*m1**5*u1**(-1)*u2
     +     - 96*h_r*st*sb*m1*m2**2*mg**2*s**(-2)*t1*u2 - 64*h_r*st*sb*
     +    m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_r*st*sb*m1*m2**2*mg**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 192*h_r*st*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_r*st*sb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*
     +    st*sb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_r*st*sb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_r*st*sb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_r*st*sb*m1*m2**2*msb1**2*s**(-2)*t1*u2 + 32*h_r*st*sb
     +    *m1*m2**2*msb1**2*s**(-2)*t1**2 + 64*h_r*st*sb*m1*m2**2*
     +    msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*sb*m1*m2**2*
     +    msb1**2*s**(-1)*t1 + 64*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*st*sb*m1*m2**2*msb1**2*t1*u1**(-1) + 
     +    64*h_r*st*sb*m1*m2**2*mst1**2*s**(-2)*t1*u2 + 32*h_r*st*sb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 + 32*h_r*st*sb*m1*m2**2*mst1**2*
     +    s**(-2)*u2**2 + 128*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*t1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)
     +    *t1 + 64*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r
     +    *st*sb*m1*m2**2*mst1**2*s**(-1)*u2 + 64*h_r*st*sb*m1*m2**2*
     +    mst1**2*t1*u1**(-1) + 96*h_r*st*sb*m1*m2**2*mst1**2*u1**(-1)*
     +    u2 + 96*h_r*st*sb*m1*m2**2*s**(-1)*t1*u2 + 32*h_r*st*sb*m1*
     +    m2**2*s**(-1)*t1**2 + 64*h_r*st*sb*m1*m2**2*s**(-1)*u2**2 + 
     +    80*h_r*st*sb*m1*m2**2*s*t1*u1**(-1) + 176*h_r*st*sb*m1*m2**2*
     +    s*u1**(-1)*u2 + 192*h_r*st*sb*m1*m2**2*t1*u1**(-1)*u2 + 64*
     +    h_r*st*sb*m1*m2**2*t1 + 64*h_r*st*sb*m1*m2**2*t1**2*u1**(-1)
     +     + 128*h_r*st*sb*m1*m2**2*u1**(-1)*u2**2 + 128*h_r*st*sb*m1*
     +    m2**2*u2 + 32*h_r*st*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_r*st*sb*m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st*
     +    sb*m1*m2**4*msb1**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*
     +    m2**4*mst1**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*st*sb*m1*m2**4*mst1**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*st*sb*m1*m2**4*s**(-2)*t1*u2 - 
     +    32*h_r*st*sb*m1*m2**4*s**(-2)*t1**2 - 32*h_r*st*sb*m1*m2**4*
     +    s**(-2)*u2**2 - 128*h_r*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_r*st*sb*m1*m2**4*s**(-1)*t1 - 64*h_r*st*sb*m1*m2**4*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*sb*m1*m2**4*s**(-1)*u2 - 96*h_r*st
     +    *sb*m1*m2**4*t1*u1**(-1) - 112*h_r*st*sb*m1*m2**4*u1**(-1)*u2
     +     + 16*h_r*st*sb*m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*
     +    m1*m2**6*s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*m1*mg**2*s**(-1)*
     +    t1*u2 + 32*h_r*st*sb*m1*mg**2*s**(-1)*u2**2 + 80*h_r*st*sb*m1
     +    *mg**2*s*u1**(-1)*u2 + 64*h_r*st*sb*m1*mg**2*t1*u1**(-1)*u2
     +     + 64*h_r*st*sb*m1*mg**2*u1**(-1)*u2**2 + 64*h_r*st*sb*m1*
     +    mg**2*u2 - 32*h_r*st*sb*m1*msb1**2*s**(-1)*t1*u2 - 32*h_r*st*
     +    sb*m1*msb1**2*s**(-1)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 80*h_r*st*sb*m1*msb1**2*s*t1*
     +    u1**(-1) - 64*h_r*st*sb*m1*msb1**2*t1*u1**(-1)*u2 - 64*h_r*st
     +    *sb*m1*msb1**2*t1 - 64*h_r*st*sb*m1*msb1**2*t1**2*u1**(-1) + 
     +    32*h_r*st*sb*m1*mst1**2*s**(-1)*t1**2 - 32*h_r*st*sb*m1*
     +    mst1**2*s**(-1)*u2**2 + 80*h_r*st*sb*m1*mst1**2*s*t1*u1**(-1)
     +     - 80*h_r*st*sb*m1*mst1**2*s*u1**(-1)*u2 + 64*h_r*st*sb*m1*
     +    mst1**2*t1 + 64*h_r*st*sb*m1*mst1**2*t1**2*u1**(-1) - 64*h_r*
     +    st*sb*m1*mst1**2*u1**(-1)*u2**2 - 64*h_r*st*sb*m1*mst1**2*u2
     +     - 64*h_r*st*sb*m1*s*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1*s*
     +    u1**(-1)*u2**2 - 64*h_r*st*sb*m1*s*u2 - 80*h_r*st*sb*m1*s**2*
     +    u1**(-1)*u2 - 32*h_r*st*sb*m1*t1*u2 - 32*h_r*st*sb*m1*u2**2
     +     - 32*h_r*st*sb*m1**3*m2**2*mg**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*st*sb*m1**3*m2**2*msb1**2*s**(-1)*u1**(-1)*u2 + 32*h_r*st
     +    *sb*m1**3*m2**2*mst1**2*s**(-1)*t1*u1**(-1) + 32*h_r*st*sb*
     +    m1**3*m2**2*mst1**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*sb*m1**3*m2**2*s**(-2)*t1*
     +    u2 + 32*h_r*st*sb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_r*st*sb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_r*st*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_r*st*sb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_r*st*sb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_r*st*sb*m1**3*m2**2*t1*u1**(-1) + 96*h_r*
     +    st*sb*m1**3*m2**2*u1**(-1)*u2 - 32*h_r*st*sb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*sb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*st*sb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_r
     +    *st*sb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_r*st*sb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*sb*m1**3*mg**2*s**(-1)*u2 - 64*h_r
     +    *st*sb*m1**3*mg**2*u1**(-1)*u2 + 96*h_r*st*sb*m1**3*msb1**2*
     +    s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*st*sb*m1**3*msb1**2*s**(-2)*
     +    t1**2 + 64*h_r*st*sb*m1**3*msb1**2*s**(-2)*u2**2 + 192*h_r*st
     +    *sb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*sb*m1**3
     +    *msb1**2*s**(-1)*t1 + 64*h_r*st*sb*m1**3*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 128*h_r*st*sb*m1**3*msb1**2*s**(-1)*u1**(-1)
     +    *u2**2 + 128*h_r*st*sb*m1**3*msb1**2*s**(-1)*u2 + 64*h_r*st*
     +    sb*m1**3*msb1**2*t1*u1**(-1) + 160*h_r*st*sb*m1**3*msb1**2*
     +    u1**(-1)*u2 - 64*h_r*st*sb*m1**3*mst1**2*s**(-2)*t1*u2 - 32*
     +    h_r*st*sb*m1**3*mst1**2*s**(-2)*t1**2 - 32*h_r*st*sb*m1**3*
     +    mst1**2*s**(-2)*u2**2 - 128*h_r*st*sb*m1**3*mst1**2*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_r*st*sb*m1**3*mst1**2*s**(-1)*t1 - 64*
     +    h_r*st*sb*m1**3*mst1**2*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb
     +    *m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*sb*m1**3*
     +    mst1**2*s**(-1)*u2 - 64*h_r*st*sb*m1**3*mst1**2*t1*u1**(-1)
     +     - 96*h_r*st*sb*m1**3*mst1**2*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*st*sb*m1**3*s**(-1)*t1*u2 + 32
     +    *h_r*st*sb*m1**3*s**(-1)*u2**2 + 64*h_r*st*sb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_r*st*sb*m1**3*t1*u1**(-1)*u2 + 64*h_r*st*
     +    sb*m1**3*u1**(-1)*u2**2 + 64*h_r*st*sb*m1**3*u2 + 16*h_r*st*
     +    sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*st*sb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*sb*m1**5*msb1**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*st*sb*m1**5*msb1**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st
     +    *sb*m1**5*mst1**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1**5*
     +    mst1**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*sb*m1**5*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*h_l*ct*cb*m1*m2**2*mg**2*s**(-2)*t1*u2
     +     + 64*h_l*ct*cb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_l*ct*cb*
     +    m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_l*ct*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 128*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)
     +    *t1 + 128*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_l*ct*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*
     +    ct*cb*m1*m2**2*mg**2*s**(-1)*u2 + 160*h_l*ct*cb*m1*m2**2*
     +    mg**2*t1*u1**(-1) + 96*h_l*ct*cb*m1*m2**2*mg**2*u1**(-1)*u2
     +     - 32*h_l*ct*cb*m1*m2**2*msb1**2*s**(-2)*t1*u2 - 32*h_l*ct*cb
     +    *m1*m2**2*msb1**2*s**(-2)*t1**2 - 64*h_l*ct*cb*m1*m2**2*
     +    msb1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*m2**2*
     +    msb1**2*s**(-1)*t1 - 64*h_l*ct*cb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_l*ct*cb*m1*m2**2*msb1**2*t1*u1**(-1) - 
     +    64*h_l*ct*cb*m1*m2**2*mst1**2*s**(-2)*t1*u2 - 32*h_l*ct*cb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*ct*cb*m1*m2**2*mst1**2*s**(-2)*
     +    u2**2 - 128*h_l*ct*cb*m1*m2**2*mst1**2*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_l*ct*cb*m1*m2**2*mst1**2*s**(-1)*t1 - 64*h_l*ct*cb*m1
     +    *m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*m2**2
     +    *mst1**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*cb*m1*m2**2*
     +    mst1**2*s**(-1)*u2 - 64*h_l*ct*cb*m1*m2**2*mst1**2*t1*
     +    u1**(-1) - 96*h_l*ct*cb*m1*m2**2*mst1**2*u1**(-1)*u2 - 96*h_l
     +    *ct*cb*m1*m2**2*s**(-1)*t1*u2 - 32*h_l*ct*cb*m1*m2**2*s**(-1)
     +    *t1**2 - 64*h_l*ct*cb*m1*m2**2*s**(-1)*u2**2 - 80*h_l*ct*cb*
     +    m1*m2**2*s*t1*u1**(-1) - 176*h_l*ct*cb*m1*m2**2*s*u1**(-1)*u2
     +     - 192*h_l*ct*cb*m1*m2**2*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*
     +    m2**2*t1 - 64*h_l*ct*cb*m1*m2**2*t1**2*u1**(-1) - 128*h_l*ct*
     +    cb*m1*m2**2*u1**(-1)*u2**2 - 128*h_l*ct*cb*m1*m2**2*u2 - 32*
     +    h_l*ct*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_l*ct*cb*m1*m2**4*msb1**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*ct*cb*m1*m2**4*mst1**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_l*ct*cb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*u2 + 64*h_l
     +    *ct*cb*m1*m2**4*s**(-2)*t1*u2 + 32*h_l*ct*cb*m1*m2**4*s**(-2)
     +    *t1**2 + 32*h_l*ct*cb*m1*m2**4*s**(-2)*u2**2 + 128*h_l*ct*cb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*m2**4*
     +    s**(-1)*t1 + 64*h_l*ct*cb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_l*ct*cb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*cb*
     +    m1*m2**4*s**(-1)*u2 + 96*h_l*ct*cb*m1*m2**4*t1*u1**(-1) + 112
     +    *h_l*ct*cb*m1*m2**4*u1**(-1)*u2 - 16*h_l*ct*cb*m1*m2**6*
     +    s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*m1*m2**6*s**(-1)*u1**(-1)*
     +    u2 - 32*h_l*ct*cb*m1*mg**2*s**(-1)*t1*u2 - 32*h_l*ct*cb*m1*
     +    mg**2*s**(-1)*u2**2 - 80*h_l*ct*cb*m1*mg**2*s*u1**(-1)*u2 - 
     +    64*h_l*ct*cb*m1*mg**2*t1*u1**(-1)*u2 - 64*h_l*ct*cb*m1*mg**2*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*cb*m1*mg**2*u2 + 32*h_l*ct*cb*
     +    m1*msb1**2*s**(-1)*t1*u2 + 32*h_l*ct*cb*m1*msb1**2*s**(-1)*
     +    t1**2 + 80*h_l*ct*cb*m1*msb1**2*s*t1*u1**(-1) + 64*h_l*ct*cb*
     +    m1*msb1**2*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1*msb1**2*t1 + 64*
     +    h_l*ct*cb*m1*msb1**2*t1**2*u1**(-1) - 32*h_l*ct*cb*m1*mst1**2
     +    *s**(-1)*t1**2 + 32*h_l*ct*cb*m1*mst1**2*s**(-1)*u2**2 - 80*
     +    h_l*ct*cb*m1*mst1**2*s*t1*u1**(-1) + 80*h_l*ct*cb*m1*mst1**2*
     +    s*u1**(-1)*u2 - 64*h_l*ct*cb*m1*mst1**2*t1 - 64*h_l*ct*cb*m1*
     +    mst1**2*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*mst1**2*u1**(-1)*
     +    u2**2 + 64*h_l*ct*cb*m1*mst1**2*u2 + 64*h_l*ct*cb*m1*s*t1*
     +    u1**(-1)*u2 + 64*h_l*ct*cb*m1*s*u1**(-1)*u2**2 + 64*h_l*ct*cb
     +    *m1*s*u2 + 80*h_l*ct*cb*m1*s**2*u1**(-1)*u2 + 32*h_l*ct*cb*m1
     +    *t1*u2 + 32*h_l*ct*cb*m1*u2**2 + 32*h_l*ct*cb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) + 32*h_l*ct*cb*m1**3*m2**2*msb1**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*ct*cb*m1**3*m2**2*mst1**2*s**(-1)
     +    *t1*u1**(-1) - 32*h_l*ct*cb*m1**3*m2**2*mst1**2*s**(-1)*
     +    u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_l
     +    *ct*cb*m1**3*m2**2*s**(-2)*t1**2 - 32*h_l*ct*cb*m1**3*m2**2*
     +    s**(-2)*u2**2 - 128*h_l*ct*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_l*ct*cb*m1**3*m2**2*s**(-1)*t1 - 64*h_l*ct*cb*
     +    m1**3*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*m1**3*m2**2
     +    *s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*cb*m1**3*m2**2*s**(-1)*u2
     +     - 96*h_l*ct*cb*m1**3*m2**2*t1*u1**(-1) - 96*h_l*ct*cb*m1**3*
     +    m2**2*u1**(-1)*u2 + 32*h_l*ct*cb*m1**3*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*ct*cb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*
     +    h_l*ct*cb*m1**3*mg**2*s**(-2)*t1*u2 + 32*h_l*ct*cb*m1**3*
     +    mg**2*s**(-2)*u2**2 + 64*h_l*ct*cb*m1**3*mg**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*ct*cb*m1**3*mg**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*cb*m1**3*mg**2*u1**(-1)*u2 - 96*
     +    h_l*ct*cb*m1**3*msb1**2*s**(-2)*t1*u2 - 32*h_l*ct*cb*m1**3*
     +    msb1**2*s**(-2)*t1**2 - 64*h_l*ct*cb*m1**3*msb1**2*s**(-2)*
     +    u2**2 - 192*h_l*ct*cb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_l*ct*cb*m1**3*msb1**2*s**(-1)*t1 - 64*h_l*ct*cb*m1**3*
     +    msb1**2*s**(-1)*t1**2*u1**(-1) - 128*h_l*ct*cb*m1**3*msb1**2*
     +    s**(-1)*u1**(-1)*u2**2 - 128*h_l*ct*cb*m1**3*msb1**2*s**(-1)*
     +    u2 - 64*h_l*ct*cb*m1**3*msb1**2*t1*u1**(-1) - 160*h_l*ct*cb*
     +    m1**3*msb1**2*u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*mst1**2*
     +    s**(-2)*t1*u2 + 32*h_l*ct*cb*m1**3*mst1**2*s**(-2)*t1**2 + 32
     +    *h_l*ct*cb*m1**3*mst1**2*s**(-2)*u2**2 + 128*h_l*ct*cb*m1**3*
     +    mst1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*cb*m1**3*mst1**2*
     +    s**(-1)*t1 + 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*ct*cb*m1**3*mst1**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*cb*m1**3*mst1**2*t1*u1**(-1) + 96
     +    *h_l*ct*cb*m1**3*mst1**2*u1**(-1)*u2 - 32*h_l*ct*cb*m1**3*
     +    s**(-1)*t1*u2 - 32*h_l*ct*cb*m1**3*s**(-1)*u2**2 - 64*h_l*ct*
     +    cb*m1**3*s*u1**(-1)*u2 - 64*h_l*ct*cb*m1**3*t1*u1**(-1)*u2 - 
     +    64*h_l*ct*cb*m1**3*u1**(-1)*u2**2 - 64*h_l*ct*cb*m1**3*u2 - 
     +    16*h_l*ct*cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*
     +    m1**5*m2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct*cb*m1**5*mg**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_l*ct*cb*m1**5*msb1**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*ct*cb*m1**5*msb1**2*s**(-1)*u1**(-1)*u2 + 
     +    16*h_l*ct*cb*m1**5*mst1**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb
     +    *m1**5*mst1**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct*cb*m1**5*
     +    u1**(-1)*u2 + 96*h_r*st*sb*m1*m2**2*mg**2*s**(-2)*t1*u2 + 64*
     +    h_r*st*sb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_r*st*sb*m1*
     +    m2**2*mg**2*s**(-2)*u2**2 + 192*h_r*st*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 128*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*t1 + 
     +    128*h_r*st*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*
     +    st*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb*m1
     +    *m2**2*mg**2*s**(-1)*u2 + 160*h_r*st*sb*m1*m2**2*mg**2*t1*
     +    u1**(-1) + 96*h_r*st*sb*m1*m2**2*mg**2*u1**(-1)*u2 - 32*h_r*
     +    st*sb*m1*m2**2*msb1**2*s**(-2)*t1*u2 - 32*h_r*st*sb*m1*m2**2*
     +    msb1**2*s**(-2)*t1**2 - 64*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)*t1 - 
     +    64*h_r*st*sb*m1*m2**2*msb1**2*s**(-1)*t1**2*u1**(-1) - 96*h_r
     +    *st*sb*m1*m2**2*msb1**2*t1*u1**(-1) - 64*h_r*st*sb*m1*m2**2*
     +    mst1**2*s**(-2)*t1*u2 - 32*h_r*st*sb*m1*m2**2*mst1**2*s**(-2)
     +    *t1**2 - 32*h_r*st*sb*m1*m2**2*mst1**2*s**(-2)*u2**2 - 128*
     +    h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*st
     +    *sb*m1*m2**2*mst1**2*s**(-1)*t1 - 64*h_r*st*sb*m1*m2**2*
     +    mst1**2*s**(-1)*t1**2*u1**(-1) )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*sb*m1*m2**2*mst1**2*s**(-1)*u2 - 
     +    64*h_r*st*sb*m1*m2**2*mst1**2*t1*u1**(-1) - 96*h_r*st*sb*m1*
     +    m2**2*mst1**2*u1**(-1)*u2 - 96*h_r*st*sb*m1*m2**2*s**(-1)*t1*
     +    u2 - 32*h_r*st*sb*m1*m2**2*s**(-1)*t1**2 - 64*h_r*st*sb*m1*
     +    m2**2*s**(-1)*u2**2 - 80*h_r*st*sb*m1*m2**2*s*t1*u1**(-1) - 
     +    176*h_r*st*sb*m1*m2**2*s*u1**(-1)*u2 - 192*h_r*st*sb*m1*m2**2
     +    *t1*u1**(-1)*u2 - 64*h_r*st*sb*m1*m2**2*t1 - 64*h_r*st*sb*m1*
     +    m2**2*t1**2*u1**(-1) - 128*h_r*st*sb*m1*m2**2*u1**(-1)*u2**2
     +     - 128*h_r*st*sb*m1*m2**2*u2 - 32*h_r*st*sb*m1*m2**4*mg**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*m2**4*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*sb*m1*m2**4*msb1**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*st*sb*m1*m2**4*mst1**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*st*sb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*u2 + 64*h_r
     +    *st*sb*m1*m2**4*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*st*sb*m1*m2**4*s**(-2)*t1**2 + 32*
     +    h_r*st*sb*m1*m2**4*s**(-2)*u2**2 + 128*h_r*st*sb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*sb*m1*m2**4*s**(-1)*t1 + 
     +    64*h_r*st*sb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb*m1*m2**4*
     +    s**(-1)*u2 + 96*h_r*st*sb*m1*m2**4*t1*u1**(-1) + 112*h_r*st*
     +    sb*m1*m2**4*u1**(-1)*u2 - 16*h_r*st*sb*m1*m2**6*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*st*sb*m1*m2**6*s**(-1)*u1**(-1)*u2 - 32*h_r
     +    *st*sb*m1*mg**2*s**(-1)*t1*u2 - 32*h_r*st*sb*m1*mg**2*s**(-1)
     +    *u2**2 - 80*h_r*st*sb*m1*mg**2*s*u1**(-1)*u2 - 64*h_r*st*sb*
     +    m1*mg**2*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1*mg**2*u1**(-1)*
     +    u2**2 - 64*h_r*st*sb*m1*mg**2*u2 + 32*h_r*st*sb*m1*msb1**2*
     +    s**(-1)*t1*u2 + 32*h_r*st*sb*m1*msb1**2*s**(-1)*t1**2 + 80*
     +    h_r*st*sb*m1*msb1**2*s*t1*u1**(-1) + 64*h_r*st*sb*m1*msb1**2*
     +    t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_r*st*sb*m1*msb1**2*t1 + 64*h_r*st*sb*
     +    m1*msb1**2*t1**2*u1**(-1) - 32*h_r*st*sb*m1*mst1**2*s**(-1)*
     +    t1**2 + 32*h_r*st*sb*m1*mst1**2*s**(-1)*u2**2 - 80*h_r*st*sb*
     +    m1*mst1**2*s*t1*u1**(-1) + 80*h_r*st*sb*m1*mst1**2*s*u1**(-1)
     +    *u2 - 64*h_r*st*sb*m1*mst1**2*t1 - 64*h_r*st*sb*m1*mst1**2*
     +    t1**2*u1**(-1) + 64*h_r*st*sb*m1*mst1**2*u1**(-1)*u2**2 + 64*
     +    h_r*st*sb*m1*mst1**2*u2 + 64*h_r*st*sb*m1*s*t1*u1**(-1)*u2 + 
     +    64*h_r*st*sb*m1*s*u1**(-1)*u2**2 + 64*h_r*st*sb*m1*s*u2 + 80*
     +    h_r*st*sb*m1*s**2*u1**(-1)*u2 + 32*h_r*st*sb*m1*t1*u2 + 32*
     +    h_r*st*sb*m1*u2**2 + 32*h_r*st*sb*m1**3*m2**2*mg**2*s**(-1)*
     +    t1*u1**(-1) + 32*h_r*st*sb*m1**3*m2**2*msb1**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*st*sb*m1**3*m2**2*mst1**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*st*sb*m1**3*m2**2*mst1**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_r*st*sb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_r*st*sb*
     +    m1**3*m2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_r*st*sb*m1**3*m2**2*s**(-2)*u2**2
     +     - 128*h_r*st*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    st*sb*m1**3*m2**2*s**(-1)*t1 - 64*h_r*st*sb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*sb*m1**3*m2**2*s**(-1)*u2 - 96*h_r
     +    *st*sb*m1**3*m2**2*t1*u1**(-1) - 96*h_r*st*sb*m1**3*m2**2*
     +    u1**(-1)*u2 + 32*h_r*st*sb*m1**3*m2**4*s**(-1)*t1*u1**(-1) + 
     +    32*h_r*st*sb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*
     +    m1**3*mg**2*s**(-2)*t1*u2 + 32*h_r*st*sb*m1**3*mg**2*s**(-2)*
     +    u2**2 + 64*h_r*st*sb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*st*sb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*sb*
     +    m1**3*mg**2*s**(-1)*u2 + 64*h_r*st*sb*m1**3*mg**2*u1**(-1)*u2
     +     - 96*h_r*st*sb*m1**3*msb1**2*s**(-2)*t1*u2 - 32*h_r*st*sb*
     +    m1**3*msb1**2*s**(-2)*t1**2 - 64*h_r*st*sb*m1**3*msb1**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 192*h_r*st*sb*m1**3*msb1**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_r*st*sb*m1**3*msb1**2*s**(-1)*t1 - 64*h_r*
     +    st*sb*m1**3*msb1**2*s**(-1)*t1**2*u1**(-1) - 128*h_r*st*sb*
     +    m1**3*msb1**2*s**(-1)*u1**(-1)*u2**2 - 128*h_r*st*sb*m1**3*
     +    msb1**2*s**(-1)*u2 - 64*h_r*st*sb*m1**3*msb1**2*t1*u1**(-1)
     +     - 160*h_r*st*sb*m1**3*msb1**2*u1**(-1)*u2 + 64*h_r*st*sb*
     +    m1**3*mst1**2*s**(-2)*t1*u2 + 32*h_r*st*sb*m1**3*mst1**2*
     +    s**(-2)*t1**2 + 32*h_r*st*sb*m1**3*mst1**2*s**(-2)*u2**2 + 
     +    128*h_r*st*sb*m1**3*mst1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*
     +    st*sb*m1**3*mst1**2*s**(-1)*t1 + 64*h_r*st*sb*m1**3*mst1**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*m1**3*mst1**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*st*sb*m1**3*mst1**2*s**(-1)*u2 + 64*
     +    h_r*st*sb*m1**3*mst1**2*t1*u1**(-1) + 96*h_r*st*sb*m1**3*
     +    mst1**2*u1**(-1)*u2 - 32*h_r*st*sb*m1**3*s**(-1)*t1*u2 - 32*
     +    h_r*st*sb*m1**3*s**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*sb*m1**3*s*u1**(-1)*u2 - 64*
     +    h_r*st*sb*m1**3*t1*u1**(-1)*u2 - 64*h_r*st*sb*m1**3*u1**(-1)*
     +    u2**2 - 64*h_r*st*sb*m1**3*u2 - 16*h_r*st*sb*m1**5*m2**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1**5*m2**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*sb*m1**5*mg**2*s**(-1)*u1**(-1)*u2 - 
     +    16*h_r*st*sb*m1**5*msb1**2*s**(-1)*t1*u1**(-1) - 32*h_r*st*sb
     +    *m1**5*msb1**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*sb*m1**5*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1**5*mst1**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_r*st*sb*m1**5*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1**2*mg*s**(-2)*t1 + 32*h_l*st*cb*
     +    m1**2*mg*s**(-2)*u2 + 32*h_l*st*cb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) + 64*h_l*st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l
     +    *st*cb*m2**2*mg*s**(-2)*t1 - 32*h_l*st*cb*m2**2*mg*s**(-2)*u2
     +     - 32*h_l*st*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*st*cb*
     +    mg*s**(-2)*t1*u2 - 32*h_l*st*cb*mg*s**(-2)*t1**2 - 64*h_l*st*
     +    cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*st*cb*mg*s**(-1)*t1 - 
     +    64*h_l*st*cb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_l*st*cb*mg*
     +    s**(-1)*u2 - 64*h_l*st*cb*mg*t1*u1**(-1) + 32*h_l*st*cb*mg*
     +    u1**(-1)*u2 - 16*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*m1*
     +    mg**2*s**(-2)*t1 - 32*h_l*ct*cb*m1*mg**2*s**(-2)*u2 - 32*h_l*
     +    ct*cb*m1*mg**2*s**(-1)*t1*u1**(-1) - 48*h_l*ct*cb*m1*mg**2*
     +    s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*m1*msb1**2*s**(-2)*t1 + 32
     +    *h_l*ct*cb*m1*msb1**2*s**(-2)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*ct*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) + 32*h_l
     +    *ct*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct*cb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*cb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 + 64*h_l*ct*cb*m1*s**(-2)*t1*u2 + 32*h_l*ct*cb*
     +    m1*s**(-2)*t1**2 + 32*h_l*ct*cb*m1*s**(-2)*u2**2 + 128*h_l*ct
     +    *cb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*cb*m1*s**(-1)*t1 + 
     +    64*h_l*ct*cb*m1*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*cb*m1*
     +    s**(-1)*u1**(-1)*u2**2 + 32*h_l*ct*cb*m1*s**(-1)*u2 + 64*h_l*
     +    ct*cb*m1*t1*u1**(-1) + 48*h_l*ct*cb*m1*u1**(-1)*u2 - 16*h_r*
     +    st*sb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2 - 32*h_r*st*sb*m1*mg**2*s**(-2)*t1 - 32*
     +    h_r*st*sb*m1*mg**2*s**(-2)*u2 - 32*h_r*st*sb*m1*mg**2*s**(-1)
     +    *t1*u1**(-1) - 48*h_r*st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 32
     +    *h_r*st*sb*m1*msb1**2*s**(-2)*t1 + 32*h_r*st*sb*m1*msb1**2*
     +    s**(-2)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_r*st*sb*m1*msb1**2*s**(-1)*t1*u1**(-1) + 32*h_r
     +    *st*sb*m1*msb1**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*sb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 + 64*h_r*st*sb*m1*s**(-2)*t1*u2 + 32*h_r*st*sb*
     +    m1*s**(-2)*t1**2 + 32*h_r*st*sb*m1*s**(-2)*u2**2 + 128*h_r*st
     +    *sb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*st*sb*m1*s**(-1)*t1 + 
     +    64*h_r*st*sb*m1*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*sb*m1*
     +    s**(-1)*u1**(-1)*u2**2 + 32*h_r*st*sb*m1*s**(-1)*u2 + 64*h_r*
     +    st*sb*m1*t1*u1**(-1) + 48*h_r*st*sb*m1*u1**(-1)*u2 + 32*h_r*
     +    ct*sb*m1**2*mg*s**(-2)*t1 + 32*h_r*ct*sb*m1**2*mg*s**(-2)*u2
     +     + 32*h_r*ct*sb*m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_r*ct*sb*
     +    m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*sb*m2**2*mg*s**(-2)*
     +    t1 - 32*h_r*ct*sb*m2**2*mg*s**(-2)*u2 - 32*h_r*ct*sb*m2**2*mg
     +    *s**(-1)*u1**(-1)*u2 - 32*h_r*ct*sb*mg*s**(-2)*t1*u2 - 32*h_r
     +    *ct*sb*mg*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 64*h_r*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*
     +    ct*sb*mg*s**(-1)*t1 - 64*h_r*ct*sb*mg*s**(-1)*t1**2*u1**(-1)
     +     + 32*h_r*ct*sb*mg*s**(-1)*u2 - 64*h_r*ct*sb*mg*t1*u1**(-1)
     +     + 32*h_r*ct*sb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1**2*mg*s**(-2)*t1 - 32*h_l*st*cb*
     +    m1**2*mg*s**(-2)*u2 - 32*h_l*st*cb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_l*st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l
     +    *st*cb*m2**2*mg*s**(-2)*t1 + 32*h_l*st*cb*m2**2*mg*s**(-2)*u2
     +     + 32*h_l*st*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*st*cb*
     +    mg*s**(-2)*t1*u2 + 32*h_l*st*cb*mg*s**(-2)*t1**2 + 64*h_l*st*
     +    cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*st*cb*mg*s**(-1)*t1 + 
     +    64*h_l*st*cb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_l*st*cb*mg*
     +    s**(-1)*u2 + 64*h_l*st*cb*mg*t1*u1**(-1) - 32*h_l*st*cb*mg*
     +    u1**(-1)*u2 + 16*h_l*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*ct*cb*m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*m1*
     +    mg**2*s**(-2)*t1 + 32*h_l*ct*cb*m1*mg**2*s**(-2)*u2 + 32*h_l*
     +    ct*cb*m1*mg**2*s**(-1)*t1*u1**(-1) + 48*h_l*ct*cb*m1*mg**2*
     +    s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*m1*msb1**2*s**(-2)*t1 - 32
     +    *h_l*ct*cb*m1*msb1**2*s**(-2)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*ct*cb*m1*msb1**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_l*ct*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct*cb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*cb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 - 64*h_l*ct*cb*m1*s**(-2)*t1*u2 - 32*h_l*ct*cb*
     +    m1*s**(-2)*t1**2 - 32*h_l*ct*cb*m1*s**(-2)*u2**2 - 128*h_l*ct
     +    *cb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*cb*m1*s**(-1)*t1 - 
     +    64*h_l*ct*cb*m1*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*cb*m1*
     +    s**(-1)*u1**(-1)*u2**2 - 32*h_l*ct*cb*m1*s**(-1)*u2 - 64*h_l*
     +    ct*cb*m1*t1*u1**(-1) - 48*h_l*ct*cb*m1*u1**(-1)*u2 + 16*h_r*
     +    st*sb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*sb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*m1*mg**2*s**(-2)*t1 + 32*
     +    h_r*st*sb*m1*mg**2*s**(-2)*u2 + 32*h_r*st*sb*m1*mg**2*s**(-1)
     +    *t1*u1**(-1) + 48*h_r*st*sb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 32
     +    *h_r*st*sb*m1*msb1**2*s**(-2)*t1 - 32*h_r*st*sb*m1*msb1**2*
     +    s**(-2)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*sb*m1*msb1**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*st*sb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st*sb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*sb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 - 64*h_r*st*sb*m1*s**(-2)*t1*u2 - 32*h_r*st*sb*
     +    m1*s**(-2)*t1**2 - 32*h_r*st*sb*m1*s**(-2)*u2**2 - 128*h_r*st
     +    *sb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*st*sb*m1*s**(-1)*t1 - 
     +    64*h_r*st*sb*m1*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*sb*m1*
     +    s**(-1)*u1**(-1)*u2**2 - 32*h_r*st*sb*m1*s**(-1)*u2 - 64*h_r*
     +    st*sb*m1*t1*u1**(-1) - 48*h_r*st*sb*m1*u1**(-1)*u2 - 32*h_r*
     +    ct*sb*m1**2*mg*s**(-2)*t1 - 32*h_r*ct*sb*m1**2*mg*s**(-2)*u2
     +     - 32*h_r*ct*sb*m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_r*ct*sb*
     +    m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*m2**2*mg*s**(-2)*
     +    t1 + 32*h_r*ct*sb*m2**2*mg*s**(-2)*u2 + 32*h_r*ct*sb*m2**2*mg
     +    *s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*mg*s**(-2)*t1*u2 + 32*h_r
     +    *ct*sb*mg*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 64*h_r*ct*sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*ct*sb
     +    *mg*s**(-1)*t1 + 64*h_r*ct*sb*mg*s**(-1)*t1**2*u1**(-1) - 32*
     +    h_r*ct*sb*mg*s**(-1)*u2 + 64*h_r*ct*sb*mg*t1*u1**(-1) - 32*
     +    h_r*ct*sb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*cb*m1*mg**2*s**(-1)*t1*
     +    u1 - 32*h_l*ct*cb*m1*mg**2*s**(-1)*t1**2 - 32*h_l*ct*cb*m1*
     +    mg**2*s**(-1)*u1**2 + 64*h_l*ct*cb*m1*msb1**2*s**(-1)*t1*u1
     +     + 32*h_l*ct*cb*m1*msb1**2*s**(-1)*t1**2 + 32*h_l*ct*cb*m1*
     +    msb1**2*s**(-1)*u1**2 + 64*h_l*ct*cb*m1*t1*u1 + 32*h_l*ct*cb*
     +    m1*t1**2 + 32*h_l*ct*cb*m1*u1**2 + 128*h_l*ct*cb*m1**3*mg**2
     +     - 128*h_l*ct*cb*m1**3*msb1**2 - 128*h_l*ct*cb*m1**3*s - 64*
     +    h_r*st*sb*m1*mg**2*s**(-1)*t1*u1 - 32*h_r*st*sb*m1*mg**2*
     +    s**(-1)*t1**2 - 32*h_r*st*sb*m1*mg**2*s**(-1)*u1**2 + 64*h_r*
     +    st*sb*m1*msb1**2*s**(-1)*t1*u1 + 32*h_r*st*sb*m1*msb1**2*
     +    s**(-1)*t1**2 + 32*h_r*st*sb*m1*msb1**2*s**(-1)*u1**2 + 64*
     +    h_r*st*sb*m1*t1*u1 + 32*h_r*st*sb*m1*t1**2 + 32*h_r*st*sb*m1*
     +    u1**2 + 128*h_r*st*sb*m1**3*mg**2 - 128*h_r*st*sb*m1**3*
     +    msb1**2 - 128*h_r*st*sb*m1**3*s )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*st*cb*mg*s**(-1)*t1 - 16*h_l*st*cb*mg*
     +    s**(-1)*u1 - 16*h_l*st*cb*mg*t1*u1**(-1) - 16*h_l*st*cb*mg - 
     +    32*h_l*ct*cb*m1*mg**2*u1**(-1) + 32*h_l*ct*cb*m1*mst1**2*
     +    u1**(-1) - 32*h_l*ct*cb*m1*s*u1**(-1) - 16*h_l*ct*cb*m1*t1*
     +    u1**(-1) - 80*h_l*ct*cb*m1 - 32*h_l*ct*cb*m1**3*u1**(-1) - 32
     +    *h_r*st*sb*m1*mg**2*u1**(-1) + 32*h_r*st*sb*m1*mst1**2*
     +    u1**(-1) - 32*h_r*st*sb*m1*s*u1**(-1) - 16*h_r*st*sb*m1*t1*
     +    u1**(-1) - 80*h_r*st*sb*m1 - 32*h_r*st*sb*m1**3*u1**(-1) - 16
     +    *h_r*ct*sb*mg*s**(-1)*t1 - 16*h_r*ct*sb*mg*s**(-1)*u1 - 16*
     +    h_r*ct*sb*mg*t1*u1**(-1) - 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(6,2)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*cb*mg*s**(-1)*t1 + 16*h_l*st*cb*mg*s**(-1)*
     +    u1 + 16*h_l*st*cb*mg*t1*u1**(-1) + 16*h_l*st*cb*mg + 32*h_l*
     +    ct*cb*m1*mg**2*s**(-1) + 32*h_l*ct*cb*m1*mg**2*u1**(-1) - 32*
     +    h_l*ct*cb*m1*msb1**2*s**(-1) - 32*h_l*ct*cb*m1*mst1**2*
     +    u1**(-1) + 32*h_l*ct*cb*m1*s*u1**(-1) + 32*h_l*ct*cb*m1*t1*
     +    u1**(-1) + 32*h_l*ct*cb*m1 + 32*h_l*ct*cb*m1**3*u1**(-1) + 32
     +    *h_r*st*sb*m1*mg**2*s**(-1) + 32*h_r*st*sb*m1*mg**2*u1**(-1)
     +     - 32*h_r*st*sb*m1*msb1**2*s**(-1) - 32*h_r*st*sb*m1*mst1**2*
     +    u1**(-1) + 32*h_r*st*sb*m1*s*u1**(-1) + 32*h_r*st*sb*m1*t1*
     +    u1**(-1) + 32*h_r*st*sb*m1 + 32*h_r*st*sb*m1**3*u1**(-1) + 16
     +    *h_r*ct*sb*mg*s**(-1)*t1 + 16*h_r*ct*sb*mg*s**(-1)*u1 + 16*
     +    h_r*ct*sb*mg*t1*u1**(-1) + 16*h_r*ct*sb*mg )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_l*ct*sb*m1*m2**2*mg**2*s**(-2)*
     +    t1*u2 + 64*h_l*ct*sb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_l*ct
     +    *sb*m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_l*ct*sb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_l*ct*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1 + 128*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*u2 + 160*h_l*ct*sb*m1*
     +    m2**2*mg**2*t1*u1**(-1) + 96*h_l*ct*sb*m1*m2**2*mg**2*
     +    u1**(-1)*u2 - 32*h_l*ct*sb*m1*m2**2*msb2**2*s**(-2)*t1*u2 - 
     +    32*h_l*ct*sb*m1*m2**2*msb2**2*s**(-2)*t1**2 - 64*h_l*ct*sb*m1
     +    *m2**2*msb2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*m2**2
     +    *msb2**2*s**(-1)*t1 - 64*h_l*ct*sb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_l*ct*sb*m1*m2**2*msb2**2*t1*u1**(-1) - 
     +    64*h_l*ct*sb*m1*m2**2*mst1**2*s**(-2)*t1*u2 - 32*h_l*ct*sb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*ct*sb*m1*m2**2*mst1**2*
     +    s**(-2)*u2**2 - 128*h_l*ct*sb*m1*m2**2*mst1**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_l*ct*sb*m1*m2**2*mst1**2*s**(-1)*t1 - 64*
     +    h_l*ct*sb*m1*m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct
     +    *sb*m1*m2**2*mst1**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*sb*m1
     +    *m2**2*mst1**2*s**(-1)*u2 - 64*h_l*ct*sb*m1*m2**2*mst1**2*t1*
     +    u1**(-1) - 96*h_l*ct*sb*m1*m2**2*mst1**2*u1**(-1)*u2 - 96*h_l
     +    *ct*sb*m1*m2**2*s**(-1)*t1*u2 - 32*h_l*ct*sb*m1*m2**2*s**(-1)
     +    *t1**2 - 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2**2 - 80*h_l*ct*sb*
     +    m1*m2**2*s*t1*u1**(-1) - 176*h_l*ct*sb*m1*m2**2*s*u1**(-1)*u2
     +     - 192*h_l*ct*sb*m1*m2**2*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*
     +    m2**2*t1 - 64*h_l*ct*sb*m1*m2**2*t1**2*u1**(-1) - 128*h_l*ct*
     +    sb*m1*m2**2*u1**(-1)*u2**2 - 128*h_l*ct*sb*m1*m2**2*u2 - 32*
     +    h_l*ct*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*h_l*ct*sb*m1*m2**4*msb2**2*s**(-1)
     +    *t1*u1**(-1) + 16*h_l*ct*sb*m1*m2**4*mst1**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*ct*sb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*u2
     +     + 64*h_l*ct*sb*m1*m2**4*s**(-2)*t1*u2 + 32*h_l*ct*sb*m1*
     +    m2**4*s**(-2)*t1**2 + 32*h_l*ct*sb*m1*m2**4*s**(-2)*u2**2 + 
     +    128*h_l*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*sb*
     +    m1*m2**4*s**(-1)*t1 + 64*h_l*ct*sb*m1*m2**4*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*
     +    h_l*ct*sb*m1*m2**4*s**(-1)*u2 + 96*h_l*ct*sb*m1*m2**4*t1*
     +    u1**(-1) + 112*h_l*ct*sb*m1*m2**4*u1**(-1)*u2 - 16*h_l*ct*sb*
     +    m1*m2**6*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*m1*m2**6*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*ct*sb*m1*mg**2*s**(-1)*t1*u2 - 32*h_l*ct
     +    *sb*m1*mg**2*s**(-1)*u2**2 - 80*h_l*ct*sb*m1*mg**2*s*u1**(-1)
     +    *u2 - 64*h_l*ct*sb*m1*mg**2*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*
     +    mg**2*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*ct*sb*m1*mg**2*u2 + 32*h_l*
     +    ct*sb*m1*msb2**2*s**(-1)*t1*u2 + 32*h_l*ct*sb*m1*msb2**2*
     +    s**(-1)*t1**2 + 80*h_l*ct*sb*m1*msb2**2*s*t1*u1**(-1) + 64*
     +    h_l*ct*sb*m1*msb2**2*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*msb2**2
     +    *t1 + 64*h_l*ct*sb*m1*msb2**2*t1**2*u1**(-1) - 32*h_l*ct*sb*
     +    m1*mst1**2*s**(-1)*t1**2 + 32*h_l*ct*sb*m1*mst1**2*s**(-1)*
     +    u2**2 - 80*h_l*ct*sb*m1*mst1**2*s*t1*u1**(-1) + 80*h_l*ct*sb*
     +    m1*mst1**2*s*u1**(-1)*u2 - 64*h_l*ct*sb*m1*mst1**2*t1 - 64*
     +    h_l*ct*sb*m1*mst1**2*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*mst1**2
     +    *u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*mst1**2*u2 + 64*h_l*ct*sb*
     +    m1*s*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*s*u1**(-1)*u2**2 + 64*
     +    h_l*ct*sb*m1*s*u2 + 80*h_l*ct*sb*m1*s**2*u1**(-1)*u2 + 32*h_l
     +    *ct*sb*m1*t1*u2 + 32*h_l*ct*sb*m1*u2**2 + 32*h_l*ct*sb*m1**3*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1**3*m2**2*
     +    msb2**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*ct*sb*m1**3*m2**2*mst1**2*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*ct*sb*m1**3*m2**2*mst1**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*m2**2*s**(-2)*t1*u2
     +     - 32*h_l*ct*sb*m1**3*m2**2*s**(-2)*t1**2 - 32*h_l*ct*sb*
     +    m1**3*m2**2*s**(-2)*u2**2 - 128*h_l*ct*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*m2**2*s**(-1)*t1 - 64*
     +    h_l*ct*sb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*sb*m1**3*m2**2
     +    *s**(-1)*u2 - 96*h_l*ct*sb*m1**3*m2**2*t1*u1**(-1) - 96*h_l*
     +    ct*sb*m1**3*m2**2*u1**(-1)*u2 + 32*h_l*ct*sb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*ct*sb*m1**3*mg**2*s**(-2)*t1*u2 + 32*h_l
     +    *ct*sb*m1**3*mg**2*s**(-2)*u2**2 + 64*h_l*ct*sb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1**3*mg**2*s**(-1)*u2
     +     + 64*h_l*ct*sb*m1**3*mg**2*u1**(-1)*u2 - 96*h_l*ct*sb*m1**3*
     +    msb2**2*s**(-2)*t1*u2 - 32*h_l*ct*sb*m1**3*msb2**2*s**(-2)*
     +    t1**2 - 64*h_l*ct*sb*m1**3*msb2**2*s**(-2)*u2**2 - 192*h_l*ct
     +    *sb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3
     +    *msb2**2*s**(-1)*t1 - 64*h_l*ct*sb*m1**3*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) - 128*h_l*ct*sb*m1**3*msb2**2*s**(-1)*u1**(-1)
     +    *u2**2 - 128*h_l*ct*sb*m1**3*msb2**2*s**(-1)*u2 - 64*h_l*ct*
     +    sb*m1**3*msb2**2*t1*u1**(-1) - 160*h_l*ct*sb*m1**3*msb2**2*
     +    u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*mst1**2*s**(-2)*t1*u2 + 32*
     +    h_l*ct*sb*m1**3*mst1**2*s**(-2)*t1**2 + 32*h_l*ct*sb*m1**3*
     +    mst1**2*s**(-2)*u2**2 + 128*h_l*ct*sb*m1**3*mst1**2*s**(-1)*
     +    t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*mst1**2*s**(-1)*t1 + 64*
     +    h_l*ct*sb*m1**3*mst1**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb
     +    *m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1**3*mst1**2*s**(-1)*u2
     +     + 64*h_l*ct*sb*m1**3*mst1**2*t1*u1**(-1) + 96*h_l*ct*sb*
     +    m1**3*mst1**2*u1**(-1)*u2 - 32*h_l*ct*sb*m1**3*s**(-1)*t1*u2
     +     - 32*h_l*ct*sb*m1**3*s**(-1)*u2**2 - 64*h_l*ct*sb*m1**3*s*
     +    u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*t1*u1**(-1)*u2 - 64*h_l*ct*
     +    sb*m1**3*u1**(-1)*u2**2 - 64*h_l*ct*sb*m1**3*u2 - 16*h_l*ct*
     +    sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 + 16*h_l*ct*sb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 16*h_l*ct*sb*m1**5*msb2**2*s**(-1)*t1*u1**(-1)
     +     - 32*h_l*ct*sb*m1**5*msb2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct
     +    *sb*m1**5*mst1**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb*m1**5*
     +    mst1**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct*sb*m1**5*u1**(-1)*u2
     +     - 96*h_r*st*cb*m1*m2**2*mg**2*s**(-2)*t1*u2 - 64*h_r*st*cb*
     +    m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_r*st*cb*m1*m2**2*mg**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 192*h_r*st*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_r*st*cb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*
     +    st*cb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_r*st*cb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_r*st*cb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_r*st*cb*m1*m2**2*msb2**2*s**(-2)*t1*u2 + 32*h_r*st*cb
     +    *m1*m2**2*msb2**2*s**(-2)*t1**2 + 64*h_r*st*cb*m1*m2**2*
     +    msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*cb*m1*m2**2*
     +    msb2**2*s**(-1)*t1 + 64*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*st*cb*m1*m2**2*msb2**2*t1*u1**(-1) + 
     +    64*h_r*st*cb*m1*m2**2*mst1**2*s**(-2)*t1*u2 + 32*h_r*st*cb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 + 32*h_r*st*cb*m1*m2**2*mst1**2*
     +    s**(-2)*u2**2 + 128*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*t1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)
     +    *t1 + 64*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r
     +    *st*cb*m1*m2**2*mst1**2*s**(-1)*u2 + 64*h_r*st*cb*m1*m2**2*
     +    mst1**2*t1*u1**(-1) + 96*h_r*st*cb*m1*m2**2*mst1**2*u1**(-1)*
     +    u2 + 96*h_r*st*cb*m1*m2**2*s**(-1)*t1*u2 + 32*h_r*st*cb*m1*
     +    m2**2*s**(-1)*t1**2 + 64*h_r*st*cb*m1*m2**2*s**(-1)*u2**2 + 
     +    80*h_r*st*cb*m1*m2**2*s*t1*u1**(-1) + 176*h_r*st*cb*m1*m2**2*
     +    s*u1**(-1)*u2 + 192*h_r*st*cb*m1*m2**2*t1*u1**(-1)*u2 + 64*
     +    h_r*st*cb*m1*m2**2*t1 + 64*h_r*st*cb*m1*m2**2*t1**2*u1**(-1)
     +     + 128*h_r*st*cb*m1*m2**2*u1**(-1)*u2**2 + 128*h_r*st*cb*m1*
     +    m2**2*u2 + 32*h_r*st*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_r*st*cb*m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st*
     +    cb*m1*m2**4*msb2**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1*
     +    m2**4*mst1**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*st*cb*m1*m2**4*mst1**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*st*cb*m1*m2**4*s**(-2)*t1*u2 - 
     +    32*h_r*st*cb*m1*m2**4*s**(-2)*t1**2 - 32*h_r*st*cb*m1*m2**4*
     +    s**(-2)*u2**2 - 128*h_r*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_r*st*cb*m1*m2**4*s**(-1)*t1 - 64*h_r*st*cb*m1*m2**4*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*st*cb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*cb*m1*m2**4*s**(-1)*u2 - 96*h_r*st
     +    *cb*m1*m2**4*t1*u1**(-1) - 112*h_r*st*cb*m1*m2**4*u1**(-1)*u2
     +     + 16*h_r*st*cb*m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*
     +    m1*m2**6*s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*m1*mg**2*s**(-1)*
     +    t1*u2 + 32*h_r*st*cb*m1*mg**2*s**(-1)*u2**2 + 80*h_r*st*cb*m1
     +    *mg**2*s*u1**(-1)*u2 + 64*h_r*st*cb*m1*mg**2*t1*u1**(-1)*u2
     +     + 64*h_r*st*cb*m1*mg**2*u1**(-1)*u2**2 + 64*h_r*st*cb*m1*
     +    mg**2*u2 - 32*h_r*st*cb*m1*msb2**2*s**(-1)*t1*u2 - 32*h_r*st*
     +    cb*m1*msb2**2*s**(-1)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 80*h_r*st*cb*m1*msb2**2*s*t1*
     +    u1**(-1) - 64*h_r*st*cb*m1*msb2**2*t1*u1**(-1)*u2 - 64*h_r*st
     +    *cb*m1*msb2**2*t1 - 64*h_r*st*cb*m1*msb2**2*t1**2*u1**(-1) + 
     +    32*h_r*st*cb*m1*mst1**2*s**(-1)*t1**2 - 32*h_r*st*cb*m1*
     +    mst1**2*s**(-1)*u2**2 + 80*h_r*st*cb*m1*mst1**2*s*t1*u1**(-1)
     +     - 80*h_r*st*cb*m1*mst1**2*s*u1**(-1)*u2 + 64*h_r*st*cb*m1*
     +    mst1**2*t1 + 64*h_r*st*cb*m1*mst1**2*t1**2*u1**(-1) - 64*h_r*
     +    st*cb*m1*mst1**2*u1**(-1)*u2**2 - 64*h_r*st*cb*m1*mst1**2*u2
     +     - 64*h_r*st*cb*m1*s*t1*u1**(-1)*u2 - 64*h_r*st*cb*m1*s*
     +    u1**(-1)*u2**2 - 64*h_r*st*cb*m1*s*u2 - 80*h_r*st*cb*m1*s**2*
     +    u1**(-1)*u2 - 32*h_r*st*cb*m1*t1*u2 - 32*h_r*st*cb*m1*u2**2
     +     - 32*h_r*st*cb*m1**3*m2**2*mg**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*st*cb*m1**3*m2**2*msb2**2*s**(-1)*u1**(-1)*u2 + 32*h_r*st
     +    *cb*m1**3*m2**2*mst1**2*s**(-1)*t1*u1**(-1) + 32*h_r*st*cb*
     +    m1**3*m2**2*mst1**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*st*cb*m1**3*m2**2*s**(-2)*t1*
     +    u2 + 32*h_r*st*cb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_r*st*cb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_r*st*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_r*st*cb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_r*st*cb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_r*st*cb*m1**3*m2**2*t1*u1**(-1) + 96*h_r*
     +    st*cb*m1**3*m2**2*u1**(-1)*u2 - 32*h_r*st*cb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*st*cb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*st*cb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_r
     +    *st*cb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_r*st*cb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_r*st*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*cb*m1**3*mg**2*s**(-1)*u2 - 64*h_r
     +    *st*cb*m1**3*mg**2*u1**(-1)*u2 + 96*h_r*st*cb*m1**3*msb2**2*
     +    s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*st*cb*m1**3*msb2**2*s**(-2)*
     +    t1**2 + 64*h_r*st*cb*m1**3*msb2**2*s**(-2)*u2**2 + 192*h_r*st
     +    *cb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*cb*m1**3
     +    *msb2**2*s**(-1)*t1 + 64*h_r*st*cb*m1**3*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 128*h_r*st*cb*m1**3*msb2**2*s**(-1)*u1**(-1)
     +    *u2**2 + 128*h_r*st*cb*m1**3*msb2**2*s**(-1)*u2 + 64*h_r*st*
     +    cb*m1**3*msb2**2*t1*u1**(-1) + 160*h_r*st*cb*m1**3*msb2**2*
     +    u1**(-1)*u2 - 64*h_r*st*cb*m1**3*mst1**2*s**(-2)*t1*u2 - 32*
     +    h_r*st*cb*m1**3*mst1**2*s**(-2)*t1**2 - 32*h_r*st*cb*m1**3*
     +    mst1**2*s**(-2)*u2**2 - 128*h_r*st*cb*m1**3*mst1**2*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_r*st*cb*m1**3*mst1**2*s**(-1)*t1 - 64*
     +    h_r*st*cb*m1**3*mst1**2*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*cb
     +    *m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*st*cb*m1**3*
     +    mst1**2*s**(-1)*u2 - 64*h_r*st*cb*m1**3*mst1**2*t1*u1**(-1)
     +     - 96*h_r*st*cb*m1**3*mst1**2*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*st*cb*m1**3*s**(-1)*t1*u2 + 32
     +    *h_r*st*cb*m1**3*s**(-1)*u2**2 + 64*h_r*st*cb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_r*st*cb*m1**3*t1*u1**(-1)*u2 + 64*h_r*st*
     +    cb*m1**3*u1**(-1)*u2**2 + 64*h_r*st*cb*m1**3*u2 + 16*h_r*st*
     +    cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*st*cb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*cb*m1**5*msb2**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*st*cb*m1**5*msb2**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st
     +    *cb*m1**5*mst1**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1**5*
     +    mst1**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*cb*m1**5*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 96*h_l*ct*sb*m1*m2**2*mg**2*s**(-2)*t1*
     +    u2 - 64*h_l*ct*sb*m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_l*ct*sb
     +    *m1*m2**2*mg**2*s**(-2)*u2**2 - 192*h_l*ct*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*ct*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*
     +    ct*sb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_l*ct*sb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_l*ct*sb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_l*ct*sb*m1*m2**2*msb2**2*s**(-2)*t1*u2 + 32*h_l*ct*sb
     +    *m1*m2**2*msb2**2*s**(-2)*t1**2 + 64*h_l*ct*sb*m1*m2**2*
     +    msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*m2**2*
     +    msb2**2*s**(-1)*t1 + 64*h_l*ct*sb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_l*ct*sb*m1*m2**2*msb2**2*t1*u1**(-1) + 
     +    64*h_l*ct*sb*m1*m2**2*mst1**2*s**(-2)*t1*u2 + 32*h_l*ct*sb*m1
     +    *m2**2*mst1**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*ct*sb*m1*m2**2*mst1**2*s**(-2)*u2**2
     +     + 128*h_l*ct*sb*m1*m2**2*mst1**2*s**(-1)*t1*u1**(-1)*u2 + 64
     +    *h_l*ct*sb*m1*m2**2*mst1**2*s**(-1)*t1 + 64*h_l*ct*sb*m1*
     +    m2**2*mst1**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*m2**2*
     +    mst1**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*ct*sb*m1*m2**2*
     +    mst1**2*s**(-1)*u2 + 64*h_l*ct*sb*m1*m2**2*mst1**2*t1*
     +    u1**(-1) + 96*h_l*ct*sb*m1*m2**2*mst1**2*u1**(-1)*u2 + 96*h_l
     +    *ct*sb*m1*m2**2*s**(-1)*t1*u2 + 32*h_l*ct*sb*m1*m2**2*s**(-1)
     +    *t1**2 + 64*h_l*ct*sb*m1*m2**2*s**(-1)*u2**2 + 80*h_l*ct*sb*
     +    m1*m2**2*s*t1*u1**(-1) + 176*h_l*ct*sb*m1*m2**2*s*u1**(-1)*u2
     +     + 192*h_l*ct*sb*m1*m2**2*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*
     +    m2**2*t1 + 64*h_l*ct*sb*m1*m2**2*t1**2*u1**(-1) + 128*h_l*ct*
     +    sb*m1*m2**2*u1**(-1)*u2**2 + 128*h_l*ct*sb*m1*m2**2*u2 + 32*
     +    h_l*ct*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*h_l*ct*sb*m1*m2**4*msb2**2*s**(-1)*t1
     +    *u1**(-1) - 16*h_l*ct*sb*m1*m2**4*mst1**2*s**(-1)*t1*u1**(-1)
     +     - 16*h_l*ct*sb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*u2 - 64*h_l
     +    *ct*sb*m1*m2**4*s**(-2)*t1*u2 - 32*h_l*ct*sb*m1*m2**4*s**(-2)
     +    *t1**2 - 32*h_l*ct*sb*m1*m2**4*s**(-2)*u2**2 - 128*h_l*ct*sb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*m2**4*
     +    s**(-1)*t1 - 64*h_l*ct*sb*m1*m2**4*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*ct*sb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*sb*
     +    m1*m2**4*s**(-1)*u2 - 96*h_l*ct*sb*m1*m2**4*t1*u1**(-1) - 112
     +    *h_l*ct*sb*m1*m2**4*u1**(-1)*u2 + 16*h_l*ct*sb*m1*m2**6*
     +    s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb*m1*m2**6*s**(-1)*u1**(-1)*
     +    u2 + 32*h_l*ct*sb*m1*mg**2*s**(-1)*t1*u2 + 32*h_l*ct*sb*m1*
     +    mg**2*s**(-1)*u2**2 + 80*h_l*ct*sb*m1*mg**2*s*u1**(-1)*u2 + 
     +    64*h_l*ct*sb*m1*mg**2*t1*u1**(-1)*u2 + 64*h_l*ct*sb*m1*mg**2*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*ct*sb*m1*mg**2*u2 - 32*h_l*ct*sb*m1*
     +    msb2**2*s**(-1)*t1*u2 - 32*h_l*ct*sb*m1*msb2**2*s**(-1)*t1**2
     +     - 80*h_l*ct*sb*m1*msb2**2*s*t1*u1**(-1) - 64*h_l*ct*sb*m1*
     +    msb2**2*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1*msb2**2*t1 - 64*h_l*
     +    ct*sb*m1*msb2**2*t1**2*u1**(-1) + 32*h_l*ct*sb*m1*mst1**2*
     +    s**(-1)*t1**2 - 32*h_l*ct*sb*m1*mst1**2*s**(-1)*u2**2 + 80*
     +    h_l*ct*sb*m1*mst1**2*s*t1*u1**(-1) - 80*h_l*ct*sb*m1*mst1**2*
     +    s*u1**(-1)*u2 + 64*h_l*ct*sb*m1*mst1**2*t1 + 64*h_l*ct*sb*m1*
     +    mst1**2*t1**2*u1**(-1) - 64*h_l*ct*sb*m1*mst1**2*u1**(-1)*
     +    u2**2 - 64*h_l*ct*sb*m1*mst1**2*u2 - 64*h_l*ct*sb*m1*s*t1*
     +    u1**(-1)*u2 - 64*h_l*ct*sb*m1*s*u1**(-1)*u2**2 - 64*h_l*ct*sb
     +    *m1*s*u2 - 80*h_l*ct*sb*m1*s**2*u1**(-1)*u2 - 32*h_l*ct*sb*m1
     +    *t1*u2 - 32*h_l*ct*sb*m1*u2**2 - 32*h_l*ct*sb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) - 32*h_l*ct*sb*m1**3*m2**2*msb2**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*ct*sb*m1**3*m2**2*mst1**2*s**(-1)*t1
     +    *u1**(-1) + 32*h_l*ct*sb*m1**3*m2**2*mst1**2*s**(-1)*u1**(-1)
     +    *u2 + 64*h_l*ct*sb*m1**3*m2**2*s**(-2)*t1*u2 + 32*h_l*ct*sb*
     +    m1**3*m2**2*s**(-2)*t1**2 + 32*h_l*ct*sb*m1**3*m2**2*s**(-2)*
     +    u2**2 + 128*h_l*ct*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64
     +    *h_l*ct*sb*m1**3*m2**2*s**(-1)*t1 + 64*h_l*ct*sb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*ct*sb*m1**3*m2**2*s**(-1)*u2 + 96*h_l
     +    *ct*sb*m1**3*m2**2*t1*u1**(-1) + 96*h_l*ct*sb*m1**3*m2**2*
     +    u1**(-1)*u2 - 32*h_l*ct*sb*m1**3*m2**4*s**(-1)*t1*u1**(-1) - 
     +    32*h_l*ct*sb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*
     +    m1**3*mg**2*s**(-2)*t1*u2 - 32*h_l*ct*sb*m1**3*mg**2*s**(-2)*
     +    u2**2 - 64*h_l*ct*sb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_l*ct*sb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*ct*sb*
     +    m1**3*mg**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*sb*m1**3*mg**2*u1**(-1)*u2 + 
     +    96*h_l*ct*sb*m1**3*msb2**2*s**(-2)*t1*u2 + 32*h_l*ct*sb*m1**3
     +    *msb2**2*s**(-2)*t1**2 + 64*h_l*ct*sb*m1**3*msb2**2*s**(-2)*
     +    u2**2 + 192*h_l*ct*sb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_l*ct*sb*m1**3*msb2**2*s**(-1)*t1 + 64*h_l*ct*sb*m1**3*
     +    msb2**2*s**(-1)*t1**2*u1**(-1) + 128*h_l*ct*sb*m1**3*msb2**2*
     +    s**(-1)*u1**(-1)*u2**2 + 128*h_l*ct*sb*m1**3*msb2**2*s**(-1)*
     +    u2 + 64*h_l*ct*sb*m1**3*msb2**2*t1*u1**(-1) + 160*h_l*ct*sb*
     +    m1**3*msb2**2*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*mst1**2*
     +    s**(-2)*t1*u2 - 32*h_l*ct*sb*m1**3*mst1**2*s**(-2)*t1**2 - 32
     +    *h_l*ct*sb*m1**3*mst1**2*s**(-2)*u2**2 - 128*h_l*ct*sb*m1**3*
     +    mst1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*ct*sb*m1**3*mst1**2*
     +    s**(-1)*t1 - 64*h_l*ct*sb*m1**3*mst1**2*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_l*ct*sb*m1**3*mst1**2*s**(-1)*u1**(-1)*u2**2
     +     - 64*h_l*ct*sb*m1**3*mst1**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*ct*sb*m1**3*mst1**2*t1*u1**(-1)
     +     - 96*h_l*ct*sb*m1**3*mst1**2*u1**(-1)*u2 + 32*h_l*ct*sb*
     +    m1**3*s**(-1)*t1*u2 + 32*h_l*ct*sb*m1**3*s**(-1)*u2**2 + 64*
     +    h_l*ct*sb*m1**3*s*u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*t1*
     +    u1**(-1)*u2 + 64*h_l*ct*sb*m1**3*u1**(-1)*u2**2 + 64*h_l*ct*
     +    sb*m1**3*u2 + 16*h_l*ct*sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_l*ct*sb*m1**5*m2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct*sb*
     +    m1**5*mg**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct*sb*m1**5*msb2**2*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*ct*sb*m1**5*msb2**2*s**(-1)*
     +    u1**(-1)*u2 - 16*h_l*ct*sb*m1**5*mst1**2*s**(-1)*t1*u1**(-1)
     +     - 16*h_l*ct*sb*m1**5*mst1**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct
     +    *sb*m1**5*u1**(-1)*u2 + 96*h_r*st*cb*m1*m2**2*mg**2*s**(-2)*
     +    t1*u2 + 64*h_r*st*cb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_r*st
     +    *cb*m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_r*st*cb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 128*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*t1 + 
     +    128*h_r*st*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*
     +    st*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*m1
     +    *m2**2*mg**2*s**(-1)*u2 + 160*h_r*st*cb*m1*m2**2*mg**2*t1*
     +    u1**(-1) + 96*h_r*st*cb*m1*m2**2*mg**2*u1**(-1)*u2 - 32*h_r*
     +    st*cb*m1*m2**2*msb2**2*s**(-2)*t1*u2 - 32*h_r*st*cb*m1*m2**2*
     +    msb2**2*s**(-2)*t1**2 - 64*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)*t1 - 
     +    64*h_r*st*cb*m1*m2**2*msb2**2*s**(-1)*t1**2*u1**(-1) - 96*h_r
     +    *st*cb*m1*m2**2*msb2**2*t1*u1**(-1) - 64*h_r*st*cb*m1*m2**2*
     +    mst1**2*s**(-2)*t1*u2 - 32*h_r*st*cb*m1*m2**2*mst1**2*s**(-2)
     +    *t1**2 - 32*h_r*st*cb*m1*m2**2*mst1**2*s**(-2)*u2**2 - 128*
     +    h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*st
     +    *cb*m1*m2**2*mst1**2*s**(-1)*t1 - 64*h_r*st*cb*m1*m2**2*
     +    mst1**2*s**(-1)*t1**2*u1**(-1) )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*cb*m1*m2**2*mst1**2*s**(-1)*u2 - 
     +    64*h_r*st*cb*m1*m2**2*mst1**2*t1*u1**(-1) - 96*h_r*st*cb*m1*
     +    m2**2*mst1**2*u1**(-1)*u2 - 96*h_r*st*cb*m1*m2**2*s**(-1)*t1*
     +    u2 - 32*h_r*st*cb*m1*m2**2*s**(-1)*t1**2 - 64*h_r*st*cb*m1*
     +    m2**2*s**(-1)*u2**2 - 80*h_r*st*cb*m1*m2**2*s*t1*u1**(-1) - 
     +    176*h_r*st*cb*m1*m2**2*s*u1**(-1)*u2 - 192*h_r*st*cb*m1*m2**2
     +    *t1*u1**(-1)*u2 - 64*h_r*st*cb*m1*m2**2*t1 - 64*h_r*st*cb*m1*
     +    m2**2*t1**2*u1**(-1) - 128*h_r*st*cb*m1*m2**2*u1**(-1)*u2**2
     +     - 128*h_r*st*cb*m1*m2**2*u2 - 32*h_r*st*cb*m1*m2**4*mg**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1*m2**4*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*cb*m1*m2**4*msb2**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*st*cb*m1*m2**4*mst1**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*st*cb*m1*m2**4*mst1**2*s**(-1)*u1**(-1)*u2 + 64*h_r
     +    *st*cb*m1*m2**4*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*st*cb*m1*m2**4*s**(-2)*t1**2 + 32*
     +    h_r*st*cb*m1*m2**4*s**(-2)*u2**2 + 128*h_r*st*cb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_r*st*cb*m1*m2**4*s**(-1)*t1 + 
     +    64*h_r*st*cb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*m1*m2**4*
     +    s**(-1)*u2 + 96*h_r*st*cb*m1*m2**4*t1*u1**(-1) + 112*h_r*st*
     +    cb*m1*m2**4*u1**(-1)*u2 - 16*h_r*st*cb*m1*m2**6*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*st*cb*m1*m2**6*s**(-1)*u1**(-1)*u2 - 32*h_r
     +    *st*cb*m1*mg**2*s**(-1)*t1*u2 - 32*h_r*st*cb*m1*mg**2*s**(-1)
     +    *u2**2 - 80*h_r*st*cb*m1*mg**2*s*u1**(-1)*u2 - 64*h_r*st*cb*
     +    m1*mg**2*t1*u1**(-1)*u2 - 64*h_r*st*cb*m1*mg**2*u1**(-1)*
     +    u2**2 - 64*h_r*st*cb*m1*mg**2*u2 + 32*h_r*st*cb*m1*msb2**2*
     +    s**(-1)*t1*u2 + 32*h_r*st*cb*m1*msb2**2*s**(-1)*t1**2 + 80*
     +    h_r*st*cb*m1*msb2**2*s*t1*u1**(-1) + 64*h_r*st*cb*m1*msb2**2*
     +    t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_r*st*cb*m1*msb2**2*t1 + 64*h_r*st*cb*
     +    m1*msb2**2*t1**2*u1**(-1) - 32*h_r*st*cb*m1*mst1**2*s**(-1)*
     +    t1**2 + 32*h_r*st*cb*m1*mst1**2*s**(-1)*u2**2 - 80*h_r*st*cb*
     +    m1*mst1**2*s*t1*u1**(-1) + 80*h_r*st*cb*m1*mst1**2*s*u1**(-1)
     +    *u2 - 64*h_r*st*cb*m1*mst1**2*t1 - 64*h_r*st*cb*m1*mst1**2*
     +    t1**2*u1**(-1) + 64*h_r*st*cb*m1*mst1**2*u1**(-1)*u2**2 + 64*
     +    h_r*st*cb*m1*mst1**2*u2 + 64*h_r*st*cb*m1*s*t1*u1**(-1)*u2 + 
     +    64*h_r*st*cb*m1*s*u1**(-1)*u2**2 + 64*h_r*st*cb*m1*s*u2 + 80*
     +    h_r*st*cb*m1*s**2*u1**(-1)*u2 + 32*h_r*st*cb*m1*t1*u2 + 32*
     +    h_r*st*cb*m1*u2**2 + 32*h_r*st*cb*m1**3*m2**2*mg**2*s**(-1)*
     +    t1*u1**(-1) + 32*h_r*st*cb*m1**3*m2**2*msb2**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*st*cb*m1**3*m2**2*mst1**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*st*cb*m1**3*m2**2*mst1**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_r*st*cb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_r*st*cb*
     +    m1**3*m2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_r*st*cb*m1**3*m2**2*s**(-2)*u2**2
     +     - 128*h_r*st*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    st*cb*m1**3*m2**2*s**(-1)*t1 - 64*h_r*st*cb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*st*cb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*st*cb*m1**3*m2**2*s**(-1)*u2 - 96*h_r
     +    *st*cb*m1**3*m2**2*t1*u1**(-1) - 96*h_r*st*cb*m1**3*m2**2*
     +    u1**(-1)*u2 + 32*h_r*st*cb*m1**3*m2**4*s**(-1)*t1*u1**(-1) + 
     +    32*h_r*st*cb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*
     +    m1**3*mg**2*s**(-2)*t1*u2 + 32*h_r*st*cb*m1**3*mg**2*s**(-2)*
     +    u2**2 + 64*h_r*st*cb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*st*cb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*st*cb*
     +    m1**3*mg**2*s**(-1)*u2 + 64*h_r*st*cb*m1**3*mg**2*u1**(-1)*u2
     +     - 96*h_r*st*cb*m1**3*msb2**2*s**(-2)*t1*u2 - 32*h_r*st*cb*
     +    m1**3*msb2**2*s**(-2)*t1**2 - 64*h_r*st*cb*m1**3*msb2**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 192*h_r*st*cb*m1**3*msb2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_r*st*cb*m1**3*msb2**2*s**(-1)*t1 - 64*h_r*
     +    st*cb*m1**3*msb2**2*s**(-1)*t1**2*u1**(-1) - 128*h_r*st*cb*
     +    m1**3*msb2**2*s**(-1)*u1**(-1)*u2**2 - 128*h_r*st*cb*m1**3*
     +    msb2**2*s**(-1)*u2 - 64*h_r*st*cb*m1**3*msb2**2*t1*u1**(-1)
     +     - 160*h_r*st*cb*m1**3*msb2**2*u1**(-1)*u2 + 64*h_r*st*cb*
     +    m1**3*mst1**2*s**(-2)*t1*u2 + 32*h_r*st*cb*m1**3*mst1**2*
     +    s**(-2)*t1**2 + 32*h_r*st*cb*m1**3*mst1**2*s**(-2)*u2**2 + 
     +    128*h_r*st*cb*m1**3*mst1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*
     +    st*cb*m1**3*mst1**2*s**(-1)*t1 + 64*h_r*st*cb*m1**3*mst1**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*m1**3*mst1**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*st*cb*m1**3*mst1**2*s**(-1)*u2 + 64*
     +    h_r*st*cb*m1**3*mst1**2*t1*u1**(-1) + 96*h_r*st*cb*m1**3*
     +    mst1**2*u1**(-1)*u2 - 32*h_r*st*cb*m1**3*s**(-1)*t1*u2 - 32*
     +    h_r*st*cb*m1**3*s**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*st*cb*m1**3*s*u1**(-1)*u2 - 64*
     +    h_r*st*cb*m1**3*t1*u1**(-1)*u2 - 64*h_r*st*cb*m1**3*u1**(-1)*
     +    u2**2 - 64*h_r*st*cb*m1**3*u2 - 16*h_r*st*cb*m1**5*m2**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1**5*m2**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*st*cb*m1**5*mg**2*s**(-1)*u1**(-1)*u2 - 
     +    16*h_r*st*cb*m1**5*msb2**2*s**(-1)*t1*u1**(-1) - 32*h_r*st*cb
     +    *m1**5*msb2**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*cb*m1**5*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1**5*mst1**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_r*st*cb*m1**5*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1**2*mg*s**(-2)*t1 - 32*h_l*st*sb*
     +    m1**2*mg*s**(-2)*u2 - 32*h_l*st*sb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_l*st*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l
     +    *st*sb*m2**2*mg*s**(-2)*t1 + 32*h_l*st*sb*m2**2*mg*s**(-2)*u2
     +     + 32*h_l*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*st*sb*
     +    mg*s**(-2)*t1*u2 + 32*h_l*st*sb*mg*s**(-2)*t1**2 + 64*h_l*st*
     +    sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*st*sb*mg*s**(-1)*t1 + 
     +    64*h_l*st*sb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_l*st*sb*mg*
     +    s**(-1)*u2 + 64*h_l*st*sb*mg*t1*u1**(-1) - 32*h_l*st*sb*mg*
     +    u1**(-1)*u2 + 16*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*
     +    h_l*ct*sb*m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*sb*m1*
     +    mg**2*s**(-2)*t1 + 32*h_l*ct*sb*m1*mg**2*s**(-2)*u2 + 32*h_l*
     +    ct*sb*m1*mg**2*s**(-1)*t1*u1**(-1) + 48*h_l*ct*sb*m1*mg**2*
     +    s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*m1*msb2**2*s**(-2)*t1 - 32
     +    *h_l*ct*sb*m1*msb2**2*s**(-2)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*ct*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_l*ct*sb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*ct*sb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 - 64*h_l*ct*sb*m1*s**(-2)*t1*u2 - 32*h_l*ct*sb*
     +    m1*s**(-2)*t1**2 - 32*h_l*ct*sb*m1*s**(-2)*u2**2 - 128*h_l*ct
     +    *sb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*sb*m1*s**(-1)*t1 - 
     +    64*h_l*ct*sb*m1*s**(-1)*t1**2*u1**(-1) - 64*h_l*ct*sb*m1*
     +    s**(-1)*u1**(-1)*u2**2 - 32*h_l*ct*sb*m1*s**(-1)*u2 - 64*h_l*
     +    ct*sb*m1*t1*u1**(-1) - 48*h_l*ct*sb*m1*u1**(-1)*u2 - 16*h_r*
     +    st*cb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2 - 32*h_r*st*cb*m1*mg**2*s**(-2)*t1 - 32*
     +    h_r*st*cb*m1*mg**2*s**(-2)*u2 - 32*h_r*st*cb*m1*mg**2*s**(-1)
     +    *t1*u1**(-1) - 48*h_r*st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 + 32
     +    *h_r*st*cb*m1*msb2**2*s**(-2)*t1 + 32*h_r*st*cb*m1*msb2**2*
     +    s**(-2)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_r*st*cb*m1*msb2**2*s**(-1)*t1*u1**(-1) + 32*h_r
     +    *st*cb*m1*msb2**2*s**(-1)*u1**(-1)*u2 + 16*h_r*st*cb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 + 64*h_r*st*cb*m1*s**(-2)*t1*u2 + 32*h_r*st*cb*
     +    m1*s**(-2)*t1**2 + 32*h_r*st*cb*m1*s**(-2)*u2**2 + 128*h_r*st
     +    *cb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*st*cb*m1*s**(-1)*t1 + 
     +    64*h_r*st*cb*m1*s**(-1)*t1**2*u1**(-1) + 64*h_r*st*cb*m1*
     +    s**(-1)*u1**(-1)*u2**2 + 32*h_r*st*cb*m1*s**(-1)*u2 + 64*h_r*
     +    st*cb*m1*t1*u1**(-1) + 48*h_r*st*cb*m1*u1**(-1)*u2 + 32*h_r*
     +    ct*cb*m1**2*mg*s**(-2)*t1 + 32*h_r*ct*cb*m1**2*mg*s**(-2)*u2
     +     + 32*h_r*ct*cb*m1**2*mg*s**(-1)*t1*u1**(-1) + 64*h_r*ct*cb*
     +    m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*cb*m2**2*mg*s**(-2)*
     +    t1 - 32*h_r*ct*cb*m2**2*mg*s**(-2)*u2 - 32*h_r*ct*cb*m2**2*mg
     +    *s**(-1)*u1**(-1)*u2 - 32*h_r*ct*cb*mg*s**(-2)*t1*u2 - 32*h_r
     +    *ct*cb*mg*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 64*h_r*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*
     +    ct*cb*mg*s**(-1)*t1 - 64*h_r*ct*cb*mg*s**(-1)*t1**2*u1**(-1)
     +     + 32*h_r*ct*cb*mg*s**(-1)*u2 - 64*h_r*ct*cb*mg*t1*u1**(-1)
     +     + 32*h_r*ct*cb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1**2*mg*s**(-2)*t1 + 32*h_l*st*sb*m1**2
     +    *mg*s**(-2)*u2 + 32*h_l*st*sb*m1**2*mg*s**(-1)*t1*u1**(-1) + 
     +    64*h_l*st*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*st*sb*
     +    m2**2*mg*s**(-2)*t1 - 32*h_l*st*sb*m2**2*mg*s**(-2)*u2 - 32*
     +    h_l*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*st*sb*mg*
     +    s**(-2)*t1*u2 - 32*h_l*st*sb*mg*s**(-2)*t1**2 - 64*h_l*st*sb*
     +    mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*st*sb*mg*s**(-1)*t1 - 64*
     +    h_l*st*sb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_l*st*sb*mg*s**(-1)
     +    *u2 - 64*h_l*st*sb*mg*t1*u1**(-1) + 32*h_l*st*sb*mg*u1**(-1)*
     +    u2 - 16*h_l*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l*ct*sb
     +    *m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*m1*mg**2*s**(-2)
     +    *t1 - 32*h_l*ct*sb*m1*mg**2*s**(-2)*u2 - 32*h_l*ct*sb*m1*
     +    mg**2*s**(-1)*t1*u1**(-1) - 48*h_l*ct*sb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*ct*sb*m1*msb2**2*s**(-2)*t1 + 32*h_l*ct*
     +    sb*m1*msb2**2*s**(-2)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*ct*sb*m1*msb2**2*s**(-1)*t1*u1**(-1) + 32*h_l*
     +    ct*sb*m1*msb2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*ct*sb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) + 16*h_l*ct*sb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 + 64*h_l*ct*sb*m1*s**(-2)*t1*u2 + 32*h_l*ct*sb*
     +    m1*s**(-2)*t1**2 + 32*h_l*ct*sb*m1*s**(-2)*u2**2 + 128*h_l*ct
     +    *sb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*sb*m1*s**(-1)*t1 + 
     +    64*h_l*ct*sb*m1*s**(-1)*t1**2*u1**(-1) + 64*h_l*ct*sb*m1*
     +    s**(-1)*u1**(-1)*u2**2 + 32*h_l*ct*sb*m1*s**(-1)*u2 + 64*h_l*
     +    ct*sb*m1*t1*u1**(-1) + 48*h_l*ct*sb*m1*u1**(-1)*u2 + 16*h_r*
     +    st*cb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*st*cb*m1*m2**2*
     +    s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*m1*mg**2*s**(-2)*t1 + 32*
     +    h_r*st*cb*m1*mg**2*s**(-2)*u2 + 32*h_r*st*cb*m1*mg**2*s**(-1)
     +    *t1*u1**(-1) + 48*h_r*st*cb*m1*mg**2*s**(-1)*u1**(-1)*u2 - 32
     +    *h_r*st*cb*m1*msb2**2*s**(-2)*t1 - 32*h_r*st*cb*m1*msb2**2*
     +    s**(-2)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*cb*m1*msb2**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*st*cb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 16*h_r*st*cb*m1*
     +    mst1**2*s**(-1)*t1*u1**(-1) - 16*h_r*st*cb*m1*mst1**2*s**(-1)
     +    *u1**(-1)*u2 - 64*h_r*st*cb*m1*s**(-2)*t1*u2 - 32*h_r*st*cb*
     +    m1*s**(-2)*t1**2 - 32*h_r*st*cb*m1*s**(-2)*u2**2 - 128*h_r*st
     +    *cb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*st*cb*m1*s**(-1)*t1 - 
     +    64*h_r*st*cb*m1*s**(-1)*t1**2*u1**(-1) - 64*h_r*st*cb*m1*
     +    s**(-1)*u1**(-1)*u2**2 - 32*h_r*st*cb*m1*s**(-1)*u2 - 64*h_r*
     +    st*cb*m1*t1*u1**(-1) - 48*h_r*st*cb*m1*u1**(-1)*u2 - 32*h_r*
     +    ct*cb*m1**2*mg*s**(-2)*t1 - 32*h_r*ct*cb*m1**2*mg*s**(-2)*u2
     +     - 32*h_r*ct*cb*m1**2*mg*s**(-1)*t1*u1**(-1) - 64*h_r*ct*cb*
     +    m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m2**2*mg*s**(-2)*
     +    t1 + 32*h_r*ct*cb*m2**2*mg*s**(-2)*u2 + 32*h_r*ct*cb*m2**2*mg
     +    *s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*mg*s**(-2)*t1*u2 + 32*h_r
     +    *ct*cb*mg*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 64*h_r*ct*cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*ct*cb
     +    *mg*s**(-1)*t1 + 64*h_r*ct*cb*mg*s**(-1)*t1**2*u1**(-1) - 32*
     +    h_r*ct*cb*mg*s**(-1)*u2 + 64*h_r*ct*cb*mg*t1*u1**(-1) - 32*
     +    h_r*ct*cb*mg*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*ct*sb*m1*mg**2*s**(-1)*t1*u1
     +     + 32*h_l*ct*sb*m1*mg**2*s**(-1)*t1**2 + 32*h_l*ct*sb*m1*
     +    mg**2*s**(-1)*u1**2 - 64*h_l*ct*sb*m1*msb2**2*s**(-1)*t1*u1
     +     - 32*h_l*ct*sb*m1*msb2**2*s**(-1)*t1**2 - 32*h_l*ct*sb*m1*
     +    msb2**2*s**(-1)*u1**2 - 64*h_l*ct*sb*m1*t1*u1 - 32*h_l*ct*sb*
     +    m1*t1**2 - 32*h_l*ct*sb*m1*u1**2 - 128*h_l*ct*sb*m1**3*mg**2
     +     + 128*h_l*ct*sb*m1**3*msb2**2 + 128*h_l*ct*sb*m1**3*s - 64*
     +    h_r*st*cb*m1*mg**2*s**(-1)*t1*u1 - 32*h_r*st*cb*m1*mg**2*
     +    s**(-1)*t1**2 - 32*h_r*st*cb*m1*mg**2*s**(-1)*u1**2 + 64*h_r*
     +    st*cb*m1*msb2**2*s**(-1)*t1*u1 + 32*h_r*st*cb*m1*msb2**2*
     +    s**(-1)*t1**2 + 32*h_r*st*cb*m1*msb2**2*s**(-1)*u1**2 + 64*
     +    h_r*st*cb*m1*t1*u1 + 32*h_r*st*cb*m1*t1**2 + 32*h_r*st*cb*m1*
     +    u1**2 + 128*h_r*st*cb*m1**3*mg**2 - 128*h_r*st*cb*m1**3*
     +    msb2**2 - 128*h_r*st*cb*m1**3*s )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*st*sb*mg*s**(-1)*t1 + 16*h_l*st*sb*mg*s**(-1)
     +    *u1 + 16*h_l*st*sb*mg*t1*u1**(-1) + 16*h_l*st*sb*mg + 32*h_l*
     +    ct*sb*m1*mg**2*u1**(-1) - 32*h_l*ct*sb*m1*mst1**2*u1**(-1) + 
     +    32*h_l*ct*sb*m1*s*u1**(-1) + 16*h_l*ct*sb*m1*t1*u1**(-1) + 80
     +    *h_l*ct*sb*m1 + 32*h_l*ct*sb*m1**3*u1**(-1) - 32*h_r*st*cb*m1
     +    *mg**2*u1**(-1) + 32*h_r*st*cb*m1*mst1**2*u1**(-1) - 32*h_r*
     +    st*cb*m1*s*u1**(-1) - 16*h_r*st*cb*m1*t1*u1**(-1) - 80*h_r*st
     +    *cb*m1 - 32*h_r*st*cb*m1**3*u1**(-1) - 16*h_r*ct*cb*mg*
     +    s**(-1)*t1 - 16*h_r*ct*cb*mg*s**(-1)*u1 - 16*h_r*ct*cb*mg*t1*
     +    u1**(-1) - 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(6,3)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*sb*mg*s**(-1)*t1 - 16*h_l*st*sb*mg*
     +    s**(-1)*u1 - 16*h_l*st*sb*mg*t1*u1**(-1) - 16*h_l*st*sb*mg - 
     +    32*h_l*ct*sb*m1*mg**2*s**(-1) - 32*h_l*ct*sb*m1*mg**2*
     +    u1**(-1) + 32*h_l*ct*sb*m1*msb2**2*s**(-1) + 32*h_l*ct*sb*m1*
     +    mst1**2*u1**(-1) - 32*h_l*ct*sb*m1*s*u1**(-1) - 32*h_l*ct*sb*
     +    m1*t1*u1**(-1) - 32*h_l*ct*sb*m1 - 32*h_l*ct*sb*m1**3*
     +    u1**(-1) + 32*h_r*st*cb*m1*mg**2*s**(-1) + 32*h_r*st*cb*m1*
     +    mg**2*u1**(-1) - 32*h_r*st*cb*m1*msb2**2*s**(-1) - 32*h_r*st*
     +    cb*m1*mst1**2*u1**(-1) + 32*h_r*st*cb*m1*s*u1**(-1) + 32*h_r*
     +    st*cb*m1*t1*u1**(-1) + 32*h_r*st*cb*m1 + 32*h_r*st*cb*m1**3*
     +    u1**(-1) + 16*h_r*ct*cb*mg*s**(-1)*t1 + 16*h_r*ct*cb*mg*
     +    s**(-1)*u1 + 16*h_r*ct*cb*mg*t1*u1**(-1) + 16*h_r*ct*cb*mg )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 96*h_l*st*cb*m1*m2**2*mg**2*s**(-2)*
     +    t1*u2 + 64*h_l*st*cb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_l*st
     +    *cb*m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_l*st*cb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2 + 128*h_l*st*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1 + 128*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*u2 + 160*h_l*st*cb*m1*
     +    m2**2*mg**2*t1*u1**(-1) + 96*h_l*st*cb*m1*m2**2*mg**2*
     +    u1**(-1)*u2 - 32*h_l*st*cb*m1*m2**2*msb1**2*s**(-2)*t1*u2 - 
     +    32*h_l*st*cb*m1*m2**2*msb1**2*s**(-2)*t1**2 - 64*h_l*st*cb*m1
     +    *m2**2*msb1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*m2**2
     +    *msb1**2*s**(-1)*t1 - 64*h_l*st*cb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_l*st*cb*m1*m2**2*msb1**2*t1*u1**(-1) - 
     +    64*h_l*st*cb*m1*m2**2*mst2**2*s**(-2)*t1*u2 - 32*h_l*st*cb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*st*cb*m1*m2**2*mst2**2*
     +    s**(-2)*u2**2 - 128*h_l*st*cb*m1*m2**2*mst2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_l*st*cb*m1*m2**2*mst2**2*s**(-1)*t1 - 64*
     +    h_l*st*cb*m1*m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st
     +    *cb*m1*m2**2*mst2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*cb*m1
     +    *m2**2*mst2**2*s**(-1)*u2 - 64*h_l*st*cb*m1*m2**2*mst2**2*t1*
     +    u1**(-1) - 96*h_l*st*cb*m1*m2**2*mst2**2*u1**(-1)*u2 - 96*h_l
     +    *st*cb*m1*m2**2*s**(-1)*t1*u2 - 32*h_l*st*cb*m1*m2**2*s**(-1)
     +    *t1**2 - 64*h_l*st*cb*m1*m2**2*s**(-1)*u2**2 - 80*h_l*st*cb*
     +    m1*m2**2*s*t1*u1**(-1) - 176*h_l*st*cb*m1*m2**2*s*u1**(-1)*u2
     +     - 192*h_l*st*cb*m1*m2**2*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*
     +    m2**2*t1 - 64*h_l*st*cb*m1*m2**2*t1**2*u1**(-1) - 128*h_l*st*
     +    cb*m1*m2**2*u1**(-1)*u2**2 - 128*h_l*st*cb*m1*m2**2*u2 - 32*
     +    h_l*st*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 16*h_l*st*cb*m1*m2**4*msb1**2*s**(-1)
     +    *t1*u1**(-1) + 16*h_l*st*cb*m1*m2**4*mst2**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*st*cb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*u2
     +     + 64*h_l*st*cb*m1*m2**4*s**(-2)*t1*u2 + 32*h_l*st*cb*m1*
     +    m2**4*s**(-2)*t1**2 + 32*h_l*st*cb*m1*m2**4*s**(-2)*u2**2 + 
     +    128*h_l*st*cb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*cb*
     +    m1*m2**4*s**(-1)*t1 + 64*h_l*st*cb*m1*m2**4*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*
     +    h_l*st*cb*m1*m2**4*s**(-1)*u2 + 96*h_l*st*cb*m1*m2**4*t1*
     +    u1**(-1) + 112*h_l*st*cb*m1*m2**4*u1**(-1)*u2 - 16*h_l*st*cb*
     +    m1*m2**6*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*m1*m2**6*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*st*cb*m1*mg**2*s**(-1)*t1*u2 - 32*h_l*st
     +    *cb*m1*mg**2*s**(-1)*u2**2 - 80*h_l*st*cb*m1*mg**2*s*u1**(-1)
     +    *u2 - 64*h_l*st*cb*m1*mg**2*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*
     +    mg**2*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*cb*m1*mg**2*u2 + 32*h_l*
     +    st*cb*m1*msb1**2*s**(-1)*t1*u2 + 32*h_l*st*cb*m1*msb1**2*
     +    s**(-1)*t1**2 + 80*h_l*st*cb*m1*msb1**2*s*t1*u1**(-1) + 64*
     +    h_l*st*cb*m1*msb1**2*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*msb1**2
     +    *t1 + 64*h_l*st*cb*m1*msb1**2*t1**2*u1**(-1) - 32*h_l*st*cb*
     +    m1*mst2**2*s**(-1)*t1**2 + 32*h_l*st*cb*m1*mst2**2*s**(-1)*
     +    u2**2 - 80*h_l*st*cb*m1*mst2**2*s*t1*u1**(-1) + 80*h_l*st*cb*
     +    m1*mst2**2*s*u1**(-1)*u2 - 64*h_l*st*cb*m1*mst2**2*t1 - 64*
     +    h_l*st*cb*m1*mst2**2*t1**2*u1**(-1) + 64*h_l*st*cb*m1*mst2**2
     +    *u1**(-1)*u2**2 + 64*h_l*st*cb*m1*mst2**2*u2 + 64*h_l*st*cb*
     +    m1*s*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*s*u1**(-1)*u2**2 + 64*
     +    h_l*st*cb*m1*s*u2 + 80*h_l*st*cb*m1*s**2*u1**(-1)*u2 + 32*h_l
     +    *st*cb*m1*t1*u2 + 32*h_l*st*cb*m1*u2**2 + 32*h_l*st*cb*m1**3*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1**3*m2**2*
     +    msb1**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 32*h_l*st*cb*m1**3*m2**2*mst2**2*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*st*cb*m1**3*m2**2*mst2**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*m2**2*s**(-2)*t1*u2
     +     - 32*h_l*st*cb*m1**3*m2**2*s**(-2)*t1**2 - 32*h_l*st*cb*
     +    m1**3*m2**2*s**(-2)*u2**2 - 128*h_l*st*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*m2**2*s**(-1)*t1 - 64*
     +    h_l*st*cb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*cb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*cb*m1**3*m2**2
     +    *s**(-1)*u2 - 96*h_l*st*cb*m1**3*m2**2*t1*u1**(-1) - 96*h_l*
     +    st*cb*m1**3*m2**2*u1**(-1)*u2 + 32*h_l*st*cb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*st*cb*m1**3*mg**2*s**(-2)*t1*u2 + 32*h_l
     +    *st*cb*m1**3*mg**2*s**(-2)*u2**2 + 64*h_l*st*cb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1**3*mg**2*s**(-1)*u2
     +     + 64*h_l*st*cb*m1**3*mg**2*u1**(-1)*u2 - 96*h_l*st*cb*m1**3*
     +    msb1**2*s**(-2)*t1*u2 - 32*h_l*st*cb*m1**3*msb1**2*s**(-2)*
     +    t1**2 - 64*h_l*st*cb*m1**3*msb1**2*s**(-2)*u2**2 - 192*h_l*st
     +    *cb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1**3
     +    *msb1**2*s**(-1)*t1 - 64*h_l*st*cb*m1**3*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) - 128*h_l*st*cb*m1**3*msb1**2*s**(-1)*u1**(-1)
     +    *u2**2 - 128*h_l*st*cb*m1**3*msb1**2*s**(-1)*u2 - 64*h_l*st*
     +    cb*m1**3*msb1**2*t1*u1**(-1) - 160*h_l*st*cb*m1**3*msb1**2*
     +    u1**(-1)*u2 + 64*h_l*st*cb*m1**3*mst2**2*s**(-2)*t1*u2 + 32*
     +    h_l*st*cb*m1**3*mst2**2*s**(-2)*t1**2 + 32*h_l*st*cb*m1**3*
     +    mst2**2*s**(-2)*u2**2 + 128*h_l*st*cb*m1**3*mst2**2*s**(-1)*
     +    t1*u1**(-1)*u2 + 64*h_l*st*cb*m1**3*mst2**2*s**(-1)*t1 + 64*
     +    h_l*st*cb*m1**3*mst2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*cb
     +    *m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1**3*mst2**2*s**(-1)*u2
     +     + 64*h_l*st*cb*m1**3*mst2**2*t1*u1**(-1) + 96*h_l*st*cb*
     +    m1**3*mst2**2*u1**(-1)*u2 - 32*h_l*st*cb*m1**3*s**(-1)*t1*u2
     +     - 32*h_l*st*cb*m1**3*s**(-1)*u2**2 - 64*h_l*st*cb*m1**3*s*
     +    u1**(-1)*u2 - 64*h_l*st*cb*m1**3*t1*u1**(-1)*u2 - 64*h_l*st*
     +    cb*m1**3*u1**(-1)*u2**2 - 64*h_l*st*cb*m1**3*u2 - 16*h_l*st*
     +    cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 + 16*h_l*st*cb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 16*h_l*st*cb*m1**5*msb1**2*s**(-1)*t1*u1**(-1)
     +     - 32*h_l*st*cb*m1**5*msb1**2*s**(-1)*u1**(-1)*u2 + 16*h_l*st
     +    *cb*m1**5*mst2**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*cb*m1**5*
     +    mst2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*st*cb*m1**5*u1**(-1)*u2
     +     - 96*h_r*ct*sb*m1*m2**2*mg**2*s**(-2)*t1*u2 - 64*h_r*ct*sb*
     +    m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_r*ct*sb*m1*m2**2*mg**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 192*h_r*ct*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_r*ct*sb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_r*ct*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*ct*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*
     +    ct*sb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_r*ct*sb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_r*ct*sb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_r*ct*sb*m1*m2**2*msb1**2*s**(-2)*t1*u2 + 32*h_r*ct*sb
     +    *m1*m2**2*msb1**2*s**(-2)*t1**2 + 64*h_r*ct*sb*m1*m2**2*
     +    msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1*m2**2*
     +    msb1**2*s**(-1)*t1 + 64*h_r*ct*sb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*ct*sb*m1*m2**2*msb1**2*t1*u1**(-1) + 
     +    64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-2)*t1*u2 + 32*h_r*ct*sb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 + 32*h_r*ct*sb*m1*m2**2*mst2**2*
     +    s**(-2)*u2**2 + 128*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*t1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)
     +    *t1 + 64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r
     +    *ct*sb*m1*m2**2*mst2**2*s**(-1)*u2 + 64*h_r*ct*sb*m1*m2**2*
     +    mst2**2*t1*u1**(-1) + 96*h_r*ct*sb*m1*m2**2*mst2**2*u1**(-1)*
     +    u2 + 96*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u2 + 32*h_r*ct*sb*m1*
     +    m2**2*s**(-1)*t1**2 + 64*h_r*ct*sb*m1*m2**2*s**(-1)*u2**2 + 
     +    80*h_r*ct*sb*m1*m2**2*s*t1*u1**(-1) + 176*h_r*ct*sb*m1*m2**2*
     +    s*u1**(-1)*u2 + 192*h_r*ct*sb*m1*m2**2*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*sb*m1*m2**2*t1 + 64*h_r*ct*sb*m1*m2**2*t1**2*u1**(-1)
     +     + 128*h_r*ct*sb*m1*m2**2*u1**(-1)*u2**2 + 128*h_r*ct*sb*m1*
     +    m2**2*u2 + 32*h_r*ct*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_r*ct*sb*m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 - 16*h_r*ct*
     +    sb*m1*m2**4*msb1**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1*
     +    m2**4*mst2**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*ct*sb*m1*m2**4*mst2**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*ct*sb*m1*m2**4*s**(-2)*t1*u2 - 
     +    32*h_r*ct*sb*m1*m2**4*s**(-2)*t1**2 - 32*h_r*ct*sb*m1*m2**4*
     +    s**(-2)*u2**2 - 128*h_r*ct*sb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_r*ct*sb*m1*m2**4*s**(-1)*t1 - 64*h_r*ct*sb*m1*m2**4*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*sb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*sb*m1*m2**4*s**(-1)*u2 - 96*h_r*ct
     +    *sb*m1*m2**4*t1*u1**(-1) - 112*h_r*ct*sb*m1*m2**4*u1**(-1)*u2
     +     + 16*h_r*ct*sb*m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*
     +    m1*m2**6*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*m1*mg**2*s**(-1)*
     +    t1*u2 + 32*h_r*ct*sb*m1*mg**2*s**(-1)*u2**2 + 80*h_r*ct*sb*m1
     +    *mg**2*s*u1**(-1)*u2 + 64*h_r*ct*sb*m1*mg**2*t1*u1**(-1)*u2
     +     + 64*h_r*ct*sb*m1*mg**2*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1*
     +    mg**2*u2 - 32*h_r*ct*sb*m1*msb1**2*s**(-1)*t1*u2 - 32*h_r*ct*
     +    sb*m1*msb1**2*s**(-1)*t1**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 80*h_r*ct*sb*m1*msb1**2*s*t1*
     +    u1**(-1) - 64*h_r*ct*sb*m1*msb1**2*t1*u1**(-1)*u2 - 64*h_r*ct
     +    *sb*m1*msb1**2*t1 - 64*h_r*ct*sb*m1*msb1**2*t1**2*u1**(-1) + 
     +    32*h_r*ct*sb*m1*mst2**2*s**(-1)*t1**2 - 32*h_r*ct*sb*m1*
     +    mst2**2*s**(-1)*u2**2 + 80*h_r*ct*sb*m1*mst2**2*s*t1*u1**(-1)
     +     - 80*h_r*ct*sb*m1*mst2**2*s*u1**(-1)*u2 + 64*h_r*ct*sb*m1*
     +    mst2**2*t1 + 64*h_r*ct*sb*m1*mst2**2*t1**2*u1**(-1) - 64*h_r*
     +    ct*sb*m1*mst2**2*u1**(-1)*u2**2 - 64*h_r*ct*sb*m1*mst2**2*u2
     +     - 64*h_r*ct*sb*m1*s*t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1*s*
     +    u1**(-1)*u2**2 - 64*h_r*ct*sb*m1*s*u2 - 80*h_r*ct*sb*m1*s**2*
     +    u1**(-1)*u2 - 32*h_r*ct*sb*m1*t1*u2 - 32*h_r*ct*sb*m1*u2**2
     +     - 32*h_r*ct*sb*m1**3*m2**2*mg**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*ct*sb*m1**3*m2**2*msb1**2*s**(-1)*u1**(-1)*u2 + 32*h_r*ct
     +    *sb*m1**3*m2**2*mst2**2*s**(-1)*t1*u1**(-1) + 32*h_r*ct*sb*
     +    m1**3*m2**2*mst2**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*sb*m1**3*m2**2*s**(-2)*t1*
     +    u2 + 32*h_r*ct*sb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_r*ct*sb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_r*ct*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_r*ct*sb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*sb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_r*ct*sb*m1**3*m2**2*t1*u1**(-1) + 96*h_r*
     +    ct*sb*m1**3*m2**2*u1**(-1)*u2 - 32*h_r*ct*sb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*ct*sb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_r
     +    *ct*sb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_r*ct*sb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*sb*m1**3*mg**2*s**(-1)*u2 - 64*h_r
     +    *ct*sb*m1**3*mg**2*u1**(-1)*u2 + 96*h_r*ct*sb*m1**3*msb1**2*
     +    s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*ct*sb*m1**3*msb1**2*s**(-2)*
     +    t1**2 + 64*h_r*ct*sb*m1**3*msb1**2*s**(-2)*u2**2 + 192*h_r*ct
     +    *sb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1**3
     +    *msb1**2*s**(-1)*t1 + 64*h_r*ct*sb*m1**3*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 128*h_r*ct*sb*m1**3*msb1**2*s**(-1)*u1**(-1)
     +    *u2**2 + 128*h_r*ct*sb*m1**3*msb1**2*s**(-1)*u2 + 64*h_r*ct*
     +    sb*m1**3*msb1**2*t1*u1**(-1) + 160*h_r*ct*sb*m1**3*msb1**2*
     +    u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*mst2**2*s**(-2)*t1*u2 - 32*
     +    h_r*ct*sb*m1**3*mst2**2*s**(-2)*t1**2 - 32*h_r*ct*sb*m1**3*
     +    mst2**2*s**(-2)*u2**2 - 128*h_r*ct*sb*m1**3*mst2**2*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*mst2**2*s**(-1)*t1 - 64*
     +    h_r*ct*sb*m1**3*mst2**2*s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*sb
     +    *m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*sb*m1**3*
     +    mst2**2*s**(-1)*u2 - 64*h_r*ct*sb*m1**3*mst2**2*t1*u1**(-1)
     +     - 96*h_r*ct*sb*m1**3*mst2**2*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*ct*sb*m1**3*s**(-1)*t1*u2 + 32
     +    *h_r*ct*sb*m1**3*s**(-1)*u2**2 + 64*h_r*ct*sb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_r*ct*sb*m1**3*t1*u1**(-1)*u2 + 64*h_r*ct*
     +    sb*m1**3*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1**3*u2 + 16*h_r*ct*
     +    sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*ct*sb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*sb*m1**5*msb1**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*ct*sb*m1**5*msb1**2*s**(-1)*u1**(-1)*u2 - 16*h_r*ct
     +    *sb*m1**5*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1**5*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*sb*m1**5*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 96*h_l*st*cb*m1*m2**2*mg**2*s**(-2)*t1*
     +    u2 - 64*h_l*st*cb*m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_l*st*cb
     +    *m1*m2**2*mg**2*s**(-2)*u2**2 - 192*h_l*st*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_l*st*cb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*st*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*
     +    st*cb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_l*st*cb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_l*st*cb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_l*st*cb*m1*m2**2*msb1**2*s**(-2)*t1*u2 + 32*h_l*st*cb
     +    *m1*m2**2*msb1**2*s**(-2)*t1**2 + 64*h_l*st*cb*m1*m2**2*
     +    msb1**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*m2**2*
     +    msb1**2*s**(-1)*t1 + 64*h_l*st*cb*m1*m2**2*msb1**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_l*st*cb*m1*m2**2*msb1**2*t1*u1**(-1) + 
     +    64*h_l*st*cb*m1*m2**2*mst2**2*s**(-2)*t1*u2 + 32*h_l*st*cb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*st*cb*m1*m2**2*mst2**2*s**(-2)*u2**2
     +     + 128*h_l*st*cb*m1*m2**2*mst2**2*s**(-1)*t1*u1**(-1)*u2 + 64
     +    *h_l*st*cb*m1*m2**2*mst2**2*s**(-1)*t1 + 64*h_l*st*cb*m1*
     +    m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*cb*m1*m2**2*
     +    mst2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*cb*m1*m2**2*
     +    mst2**2*s**(-1)*u2 + 64*h_l*st*cb*m1*m2**2*mst2**2*t1*
     +    u1**(-1) + 96*h_l*st*cb*m1*m2**2*mst2**2*u1**(-1)*u2 + 96*h_l
     +    *st*cb*m1*m2**2*s**(-1)*t1*u2 + 32*h_l*st*cb*m1*m2**2*s**(-1)
     +    *t1**2 + 64*h_l*st*cb*m1*m2**2*s**(-1)*u2**2 + 80*h_l*st*cb*
     +    m1*m2**2*s*t1*u1**(-1) + 176*h_l*st*cb*m1*m2**2*s*u1**(-1)*u2
     +     + 192*h_l*st*cb*m1*m2**2*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*
     +    m2**2*t1 + 64*h_l*st*cb*m1*m2**2*t1**2*u1**(-1) + 128*h_l*st*
     +    cb*m1*m2**2*u1**(-1)*u2**2 + 128*h_l*st*cb*m1*m2**2*u2 + 32*
     +    h_l*st*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*cb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 16*h_l*st*cb*m1*m2**4*msb1**2*s**(-1)*t1
     +    *u1**(-1) - 16*h_l*st*cb*m1*m2**4*mst2**2*s**(-1)*t1*u1**(-1)
     +     - 16*h_l*st*cb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*u2 - 64*h_l
     +    *st*cb*m1*m2**4*s**(-2)*t1*u2 - 32*h_l*st*cb*m1*m2**4*s**(-2)
     +    *t1**2 - 32*h_l*st*cb*m1*m2**4*s**(-2)*u2**2 - 128*h_l*st*cb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*m2**4*
     +    s**(-1)*t1 - 64*h_l*st*cb*m1*m2**4*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*st*cb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*cb*
     +    m1*m2**4*s**(-1)*u2 - 96*h_l*st*cb*m1*m2**4*t1*u1**(-1) - 112
     +    *h_l*st*cb*m1*m2**4*u1**(-1)*u2 + 16*h_l*st*cb*m1*m2**6*
     +    s**(-1)*t1*u1**(-1) + 16*h_l*st*cb*m1*m2**6*s**(-1)*u1**(-1)*
     +    u2 + 32*h_l*st*cb*m1*mg**2*s**(-1)*t1*u2 + 32*h_l*st*cb*m1*
     +    mg**2*s**(-1)*u2**2 + 80*h_l*st*cb*m1*mg**2*s*u1**(-1)*u2 + 
     +    64*h_l*st*cb*m1*mg**2*t1*u1**(-1)*u2 + 64*h_l*st*cb*m1*mg**2*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*cb*m1*mg**2*u2 - 32*h_l*st*cb*m1*
     +    msb1**2*s**(-1)*t1*u2 - 32*h_l*st*cb*m1*msb1**2*s**(-1)*t1**2
     +     - 80*h_l*st*cb*m1*msb1**2*s*t1*u1**(-1) - 64*h_l*st*cb*m1*
     +    msb1**2*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1*msb1**2*t1 - 64*h_l*
     +    st*cb*m1*msb1**2*t1**2*u1**(-1) + 32*h_l*st*cb*m1*mst2**2*
     +    s**(-1)*t1**2 - 32*h_l*st*cb*m1*mst2**2*s**(-1)*u2**2 + 80*
     +    h_l*st*cb*m1*mst2**2*s*t1*u1**(-1) - 80*h_l*st*cb*m1*mst2**2*
     +    s*u1**(-1)*u2 + 64*h_l*st*cb*m1*mst2**2*t1 + 64*h_l*st*cb*m1*
     +    mst2**2*t1**2*u1**(-1) - 64*h_l*st*cb*m1*mst2**2*u1**(-1)*
     +    u2**2 - 64*h_l*st*cb*m1*mst2**2*u2 - 64*h_l*st*cb*m1*s*t1*
     +    u1**(-1)*u2 - 64*h_l*st*cb*m1*s*u1**(-1)*u2**2 - 64*h_l*st*cb
     +    *m1*s*u2 - 80*h_l*st*cb*m1*s**2*u1**(-1)*u2 - 32*h_l*st*cb*m1
     +    *t1*u2 - 32*h_l*st*cb*m1*u2**2 - 32*h_l*st*cb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) - 32*h_l*st*cb*m1**3*m2**2*msb1**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_l*st*cb*m1**3*m2**2*mst2**2*s**(-1)*t1
     +    *u1**(-1) + 32*h_l*st*cb*m1**3*m2**2*mst2**2*s**(-1)*u1**(-1)
     +    *u2 + 64*h_l*st*cb*m1**3*m2**2*s**(-2)*t1*u2 + 32*h_l*st*cb*
     +    m1**3*m2**2*s**(-2)*t1**2 + 32*h_l*st*cb*m1**3*m2**2*s**(-2)*
     +    u2**2 + 128*h_l*st*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 + 64
     +    *h_l*st*cb*m1**3*m2**2*s**(-1)*t1 + 64*h_l*st*cb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_l*st*cb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_l*st*cb*m1**3*m2**2*s**(-1)*u2 + 96*h_l
     +    *st*cb*m1**3*m2**2*t1*u1**(-1) + 96*h_l*st*cb*m1**3*m2**2*
     +    u1**(-1)*u2 - 32*h_l*st*cb*m1**3*m2**4*s**(-1)*t1*u1**(-1) - 
     +    32*h_l*st*cb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 - 32*h_l*st*cb*
     +    m1**3*mg**2*s**(-2)*t1*u2 - 32*h_l*st*cb*m1**3*mg**2*s**(-2)*
     +    u2**2 - 64*h_l*st*cb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 - 64*
     +    h_l*st*cb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*cb*
     +    m1**3*mg**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*cb*m1**3*mg**2*u1**(-1)*u2 + 
     +    96*h_l*st*cb*m1**3*msb1**2*s**(-2)*t1*u2 + 32*h_l*st*cb*m1**3
     +    *msb1**2*s**(-2)*t1**2 + 64*h_l*st*cb*m1**3*msb1**2*s**(-2)*
     +    u2**2 + 192*h_l*st*cb*m1**3*msb1**2*s**(-1)*t1*u1**(-1)*u2 + 
     +    64*h_l*st*cb*m1**3*msb1**2*s**(-1)*t1 + 64*h_l*st*cb*m1**3*
     +    msb1**2*s**(-1)*t1**2*u1**(-1) + 128*h_l*st*cb*m1**3*msb1**2*
     +    s**(-1)*u1**(-1)*u2**2 + 128*h_l*st*cb*m1**3*msb1**2*s**(-1)*
     +    u2 + 64*h_l*st*cb*m1**3*msb1**2*t1*u1**(-1) + 160*h_l*st*cb*
     +    m1**3*msb1**2*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*mst2**2*
     +    s**(-2)*t1*u2 - 32*h_l*st*cb*m1**3*mst2**2*s**(-2)*t1**2 - 32
     +    *h_l*st*cb*m1**3*mst2**2*s**(-2)*u2**2 - 128*h_l*st*cb*m1**3*
     +    mst2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*cb*m1**3*mst2**2*
     +    s**(-1)*t1 - 64*h_l*st*cb*m1**3*mst2**2*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_l*st*cb*m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2
     +     - 64*h_l*st*cb*m1**3*mst2**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*cb*m1**3*mst2**2*t1*u1**(-1)
     +     - 96*h_l*st*cb*m1**3*mst2**2*u1**(-1)*u2 + 32*h_l*st*cb*
     +    m1**3*s**(-1)*t1*u2 + 32*h_l*st*cb*m1**3*s**(-1)*u2**2 + 64*
     +    h_l*st*cb*m1**3*s*u1**(-1)*u2 + 64*h_l*st*cb*m1**3*t1*
     +    u1**(-1)*u2 + 64*h_l*st*cb*m1**3*u1**(-1)*u2**2 + 64*h_l*st*
     +    cb*m1**3*u2 + 16*h_l*st*cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_l*st*cb*m1**5*m2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*st*cb*
     +    m1**5*mg**2*s**(-1)*u1**(-1)*u2 + 16*h_l*st*cb*m1**5*msb1**2*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*st*cb*m1**5*msb1**2*s**(-1)*
     +    u1**(-1)*u2 - 16*h_l*st*cb*m1**5*mst2**2*s**(-1)*t1*u1**(-1)
     +     - 16*h_l*st*cb*m1**5*mst2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*st
     +    *cb*m1**5*u1**(-1)*u2 + 96*h_r*ct*sb*m1*m2**2*mg**2*s**(-2)*
     +    t1*u2 + 64*h_r*ct*sb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_r*ct
     +    *sb*m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_r*ct*sb*m1*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 128*h_r*ct*sb*m1*m2**2*mg**2*s**(-1)*t1 + 
     +    128*h_r*ct*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*
     +    ct*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1
     +    *m2**2*mg**2*s**(-1)*u2 + 160*h_r*ct*sb*m1*m2**2*mg**2*t1*
     +    u1**(-1) + 96*h_r*ct*sb*m1*m2**2*mg**2*u1**(-1)*u2 - 32*h_r*
     +    ct*sb*m1*m2**2*msb1**2*s**(-2)*t1*u2 - 32*h_r*ct*sb*m1*m2**2*
     +    msb1**2*s**(-2)*t1**2 - 64*h_r*ct*sb*m1*m2**2*msb1**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1*m2**2*msb1**2*s**(-1)*t1 - 
     +    64*h_r*ct*sb*m1*m2**2*msb1**2*s**(-1)*t1**2*u1**(-1) - 96*h_r
     +    *ct*sb*m1*m2**2*msb1**2*t1*u1**(-1) - 64*h_r*ct*sb*m1*m2**2*
     +    mst2**2*s**(-2)*t1*u2 - 32*h_r*ct*sb*m1*m2**2*mst2**2*s**(-2)
     +    *t1**2 - 32*h_r*ct*sb*m1*m2**2*mst2**2*s**(-2)*u2**2 - 128*
     +    h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*ct
     +    *sb*m1*m2**2*mst2**2*s**(-1)*t1 - 64*h_r*ct*sb*m1*m2**2*
     +    mst2**2*s**(-1)*t1**2*u1**(-1) )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*sb*m1*m2**2*mst2**2*s**(-1)*u2 - 
     +    64*h_r*ct*sb*m1*m2**2*mst2**2*t1*u1**(-1) - 96*h_r*ct*sb*m1*
     +    m2**2*mst2**2*u1**(-1)*u2 - 96*h_r*ct*sb*m1*m2**2*s**(-1)*t1*
     +    u2 - 32*h_r*ct*sb*m1*m2**2*s**(-1)*t1**2 - 64*h_r*ct*sb*m1*
     +    m2**2*s**(-1)*u2**2 - 80*h_r*ct*sb*m1*m2**2*s*t1*u1**(-1) - 
     +    176*h_r*ct*sb*m1*m2**2*s*u1**(-1)*u2 - 192*h_r*ct*sb*m1*m2**2
     +    *t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1*m2**2*t1 - 64*h_r*ct*sb*m1*
     +    m2**2*t1**2*u1**(-1) - 128*h_r*ct*sb*m1*m2**2*u1**(-1)*u2**2
     +     - 128*h_r*ct*sb*m1*m2**2*u2 - 32*h_r*ct*sb*m1*m2**4*mg**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1*m2**4*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*sb*m1*m2**4*msb1**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*ct*sb*m1*m2**4*mst2**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*ct*sb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_r
     +    *ct*sb*m1*m2**4*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*ct*sb*m1*m2**4*s**(-2)*t1**2 + 32*
     +    h_r*ct*sb*m1*m2**4*s**(-2)*u2**2 + 128*h_r*ct*sb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*sb*m1*m2**4*s**(-1)*t1 + 
     +    64*h_r*ct*sb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*sb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1*m2**4*
     +    s**(-1)*u2 + 96*h_r*ct*sb*m1*m2**4*t1*u1**(-1) + 112*h_r*ct*
     +    sb*m1*m2**4*u1**(-1)*u2 - 16*h_r*ct*sb*m1*m2**6*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*ct*sb*m1*m2**6*s**(-1)*u1**(-1)*u2 - 32*h_r
     +    *ct*sb*m1*mg**2*s**(-1)*t1*u2 - 32*h_r*ct*sb*m1*mg**2*s**(-1)
     +    *u2**2 - 80*h_r*ct*sb*m1*mg**2*s*u1**(-1)*u2 - 64*h_r*ct*sb*
     +    m1*mg**2*t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1*mg**2*u1**(-1)*
     +    u2**2 - 64*h_r*ct*sb*m1*mg**2*u2 + 32*h_r*ct*sb*m1*msb1**2*
     +    s**(-1)*t1*u2 + 32*h_r*ct*sb*m1*msb1**2*s**(-1)*t1**2 + 80*
     +    h_r*ct*sb*m1*msb1**2*s*t1*u1**(-1) + 64*h_r*ct*sb*m1*msb1**2*
     +    t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_r*ct*sb*m1*msb1**2*t1 + 64*h_r*ct*sb*
     +    m1*msb1**2*t1**2*u1**(-1) - 32*h_r*ct*sb*m1*mst2**2*s**(-1)*
     +    t1**2 + 32*h_r*ct*sb*m1*mst2**2*s**(-1)*u2**2 - 80*h_r*ct*sb*
     +    m1*mst2**2*s*t1*u1**(-1) + 80*h_r*ct*sb*m1*mst2**2*s*u1**(-1)
     +    *u2 - 64*h_r*ct*sb*m1*mst2**2*t1 - 64*h_r*ct*sb*m1*mst2**2*
     +    t1**2*u1**(-1) + 64*h_r*ct*sb*m1*mst2**2*u1**(-1)*u2**2 + 64*
     +    h_r*ct*sb*m1*mst2**2*u2 + 64*h_r*ct*sb*m1*s*t1*u1**(-1)*u2 + 
     +    64*h_r*ct*sb*m1*s*u1**(-1)*u2**2 + 64*h_r*ct*sb*m1*s*u2 + 80*
     +    h_r*ct*sb*m1*s**2*u1**(-1)*u2 + 32*h_r*ct*sb*m1*t1*u2 + 32*
     +    h_r*ct*sb*m1*u2**2 + 32*h_r*ct*sb*m1**3*m2**2*mg**2*s**(-1)*
     +    t1*u1**(-1) + 32*h_r*ct*sb*m1**3*m2**2*msb1**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*ct*sb*m1**3*m2**2*mst2**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*ct*sb*m1**3*m2**2*mst2**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_r*ct*sb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_r*ct*sb*
     +    m1**3*m2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_r*ct*sb*m1**3*m2**2*s**(-2)*u2**2
     +     - 128*h_r*ct*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    ct*sb*m1**3*m2**2*s**(-1)*t1 - 64*h_r*ct*sb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*sb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*sb*m1**3*m2**2*s**(-1)*u2 - 96*h_r
     +    *ct*sb*m1**3*m2**2*t1*u1**(-1) - 96*h_r*ct*sb*m1**3*m2**2*
     +    u1**(-1)*u2 + 32*h_r*ct*sb*m1**3*m2**4*s**(-1)*t1*u1**(-1) + 
     +    32*h_r*ct*sb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*
     +    m1**3*mg**2*s**(-2)*t1*u2 + 32*h_r*ct*sb*m1**3*mg**2*s**(-2)*
     +    u2**2 + 64*h_r*ct*sb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*sb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*sb*
     +    m1**3*mg**2*s**(-1)*u2 + 64*h_r*ct*sb*m1**3*mg**2*u1**(-1)*u2
     +     - 96*h_r*ct*sb*m1**3*msb1**2*s**(-2)*t1*u2 - 32*h_r*ct*sb*
     +    m1**3*msb1**2*s**(-2)*t1**2 - 64*h_r*ct*sb*m1**3*msb1**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 192*h_r*ct*sb*m1**3*msb1**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*msb1**2*s**(-1)*t1 - 64*h_r*
     +    ct*sb*m1**3*msb1**2*s**(-1)*t1**2*u1**(-1) - 128*h_r*ct*sb*
     +    m1**3*msb1**2*s**(-1)*u1**(-1)*u2**2 - 128*h_r*ct*sb*m1**3*
     +    msb1**2*s**(-1)*u2 - 64*h_r*ct*sb*m1**3*msb1**2*t1*u1**(-1)
     +     - 160*h_r*ct*sb*m1**3*msb1**2*u1**(-1)*u2 + 64*h_r*ct*sb*
     +    m1**3*mst2**2*s**(-2)*t1*u2 + 32*h_r*ct*sb*m1**3*mst2**2*
     +    s**(-2)*t1**2 + 32*h_r*ct*sb*m1**3*mst2**2*s**(-2)*u2**2 + 
     +    128*h_r*ct*sb*m1**3*mst2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*
     +    ct*sb*m1**3*mst2**2*s**(-1)*t1 + 64*h_r*ct*sb*m1**3*mst2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*sb*m1**3*mst2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*ct*sb*m1**3*mst2**2*s**(-1)*u2 + 64*
     +    h_r*ct*sb*m1**3*mst2**2*t1*u1**(-1) + 96*h_r*ct*sb*m1**3*
     +    mst2**2*u1**(-1)*u2 - 32*h_r*ct*sb*m1**3*s**(-1)*t1*u2 - 32*
     +    h_r*ct*sb*m1**3*s**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*sb*m1**3*s*u1**(-1)*u2 - 64*
     +    h_r*ct*sb*m1**3*t1*u1**(-1)*u2 - 64*h_r*ct*sb*m1**3*u1**(-1)*
     +    u2**2 - 64*h_r*ct*sb*m1**3*u2 - 16*h_r*ct*sb*m1**5*m2**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1**5*m2**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*sb*m1**5*mg**2*s**(-1)*u1**(-1)*u2 - 
     +    16*h_r*ct*sb*m1**5*msb1**2*s**(-1)*t1*u1**(-1) - 32*h_r*ct*sb
     +    *m1**5*msb1**2*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*sb*m1**5*
     +    mst2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*m1**5*mst2**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_r*ct*sb*m1**5*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_l*
     +    st*cb*m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_l*st*cb*m1*mg**2*
     +    s**(-2)*t1 + 32*h_l*st*cb*m1*mg**2*s**(-2)*u2 + 32*h_l*st*cb*
     +    m1*mg**2*s**(-1)*t1*u1**(-1) + 48*h_l*st*cb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*st*cb*m1*msb1**2*s**(-2)*t1 - 32*h_l*st*
     +    cb*m1*msb1**2*s**(-2)*u2 - 16*h_l*st*cb*m1*msb1**2*s**(-1)*t1
     +    *u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 16*
     +    h_l*st*cb*m1*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*cb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 - 64*h_l*st*cb*m1*s**(-2)*t1*u2
     +     - 32*h_l*st*cb*m1*s**(-2)*t1**2 - 32*h_l*st*cb*m1*s**(-2)*
     +    u2**2 - 128*h_l*st*cb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*st*
     +    cb*m1*s**(-1)*t1 - 64*h_l*st*cb*m1*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*st*cb*m1*s**(-1)*u1**(-1)*u2**2 - 32*h_l*st*cb*m1*
     +    s**(-1)*u2 - 64*h_l*st*cb*m1*t1*u1**(-1) - 48*h_l*st*cb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*cb*m1**2*mg*s**(-2)*t1 + 32*h_l*ct*cb*
     +    m1**2*mg*s**(-2)*u2 + 32*h_l*ct*cb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) + 64*h_l*ct*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l
     +    *ct*cb*m2**2*mg*s**(-2)*t1 - 32*h_l*ct*cb*m2**2*mg*s**(-2)*u2
     +     - 32*h_l*ct*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*cb*
     +    mg*s**(-2)*t1*u2 - 32*h_l*ct*cb*mg*s**(-2)*t1**2 - 64*h_l*ct*
     +    cb*mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*cb*mg*s**(-1)*t1 - 
     +    64*h_l*ct*cb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_l*ct*cb*mg*
     +    s**(-1)*u2 - 64*h_l*ct*cb*mg*t1*u1**(-1) + 32*h_l*ct*cb*mg*
     +    u1**(-1)*u2 - 32*h_r*st*sb*m1**2*mg*s**(-2)*t1 - 32*h_r*st*sb
     +    *m1**2*mg*s**(-2)*u2 - 32*h_r*st*sb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_r*st*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *st*sb*m2**2*mg*s**(-2)*t1 + 32*h_r*st*sb*m2**2*mg*s**(-2)*u2
     +     + 32*h_r*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*st*sb*
     +    mg*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_r*st*sb*mg*s**(-2)*t1**2 + 64*h_r*st*sb*mg*
     +    s**(-1)*t1*u1**(-1)*u2 + 32*h_r*st*sb*mg*s**(-1)*t1 + 64*h_r*
     +    st*sb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_r*st*sb*mg*s**(-1)*u2
     +     + 64*h_r*st*sb*mg*t1*u1**(-1) - 32*h_r*st*sb*mg*u1**(-1)*u2
     +     - 16*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*sb*m1*mg**2*s**(-2)*
     +    t1 - 32*h_r*ct*sb*m1*mg**2*s**(-2)*u2 - 32*h_r*ct*sb*m1*mg**2
     +    *s**(-1)*t1*u1**(-1) - 48*h_r*ct*sb*m1*mg**2*s**(-1)*u1**(-1)
     +    *u2 + 32*h_r*ct*sb*m1*msb1**2*s**(-2)*t1 + 32*h_r*ct*sb*m1*
     +    msb1**2*s**(-2)*u2 + 16*h_r*ct*sb*m1*msb1**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*ct*sb*m1*msb1**2*s**(-1)*u1**(-1)*u2 + 16*
     +    h_r*ct*sb*m1*mst2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*ct*sb*m1*s**(-2)*t1*u2
     +     + 32*h_r*ct*sb*m1*s**(-2)*t1**2 + 32*h_r*ct*sb*m1*s**(-2)*
     +    u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 128*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*ct*
     +    sb*m1*s**(-1)*t1 + 64*h_r*ct*sb*m1*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2**2 + 32*h_r*ct*sb*m1*
     +    s**(-1)*u2 + 64*h_r*ct*sb*m1*t1*u1**(-1) + 48*h_r*ct*sb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*cb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l
     +    *st*cb*m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_l*st*cb*m1*mg**2*
     +    s**(-2)*t1 - 32*h_l*st*cb*m1*mg**2*s**(-2)*u2 - 32*h_l*st*cb*
     +    m1*mg**2*s**(-1)*t1*u1**(-1) - 48*h_l*st*cb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*st*cb*m1*msb1**2*s**(-2)*t1 + 32*h_l*st*
     +    cb*m1*msb1**2*s**(-2)*u2 + 16*h_l*st*cb*m1*msb1**2*s**(-1)*t1
     +    *u1**(-1) + 32*h_l*st*cb*m1*msb1**2*s**(-1)*u1**(-1)*u2 + 16*
     +    h_l*st*cb*m1*mst2**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*cb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_l*st*cb*m1*s**(-2)*t1*u2
     +     + 32*h_l*st*cb*m1*s**(-2)*t1**2 + 32*h_l*st*cb*m1*s**(-2)*
     +    u2**2 + 128*h_l*st*cb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*st*
     +    cb*m1*s**(-1)*t1 + 64*h_l*st*cb*m1*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_l*st*cb*m1*s**(-1)*u1**(-1)*u2**2 + 32*h_l*st*cb*m1*
     +    s**(-1)*u2 + 64*h_l*st*cb*m1*t1*u1**(-1) + 48*h_l*st*cb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1**2*mg*s**(-2)*t1 - 32*h_l*ct*cb*
     +    m1**2*mg*s**(-2)*u2 - 32*h_l*ct*cb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_l*ct*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l
     +    *ct*cb*m2**2*mg*s**(-2)*t1 + 32*h_l*ct*cb*m2**2*mg*s**(-2)*u2
     +     + 32*h_l*ct*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*cb*
     +    mg*s**(-2)*t1*u2 + 32*h_l*ct*cb*mg*s**(-2)*t1**2 + 64*h_l*ct*
     +    cb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*cb*mg*s**(-1)*t1 + 
     +    64*h_l*ct*cb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_l*ct*cb*mg*
     +    s**(-1)*u2 + 64*h_l*ct*cb*mg*t1*u1**(-1) - 32*h_l*ct*cb*mg*
     +    u1**(-1)*u2 + 32*h_r*st*sb*m1**2*mg*s**(-2)*t1 + 32*h_r*st*sb
     +    *m1**2*mg*s**(-2)*u2 + 32*h_r*st*sb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) + 64*h_r*st*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r
     +    *st*sb*m2**2*mg*s**(-2)*t1 - 32*h_r*st*sb*m2**2*mg*s**(-2)*u2
     +     - 32*h_r*st*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*st*sb*
     +    mg*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_r*st*sb*mg*s**(-2)*t1**2 - 64*h_r*st*sb*mg*
     +    s**(-1)*t1*u1**(-1)*u2 - 32*h_r*st*sb*mg*s**(-1)*t1 - 64*h_r*
     +    st*sb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_r*st*sb*mg*s**(-1)*u2
     +     - 64*h_r*st*sb*mg*t1*u1**(-1) + 32*h_r*st*sb*mg*u1**(-1)*u2
     +     + 16*h_r*ct*sb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*sb*
     +    m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*sb*m1*mg**2*s**(-2)*
     +    t1 + 32*h_r*ct*sb*m1*mg**2*s**(-2)*u2 + 32*h_r*ct*sb*m1*mg**2
     +    *s**(-1)*t1*u1**(-1) + 48*h_r*ct*sb*m1*mg**2*s**(-1)*u1**(-1)
     +    *u2 - 32*h_r*ct*sb*m1*msb1**2*s**(-2)*t1 - 32*h_r*ct*sb*m1*
     +    msb1**2*s**(-2)*u2 - 16*h_r*ct*sb*m1*msb1**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*ct*sb*m1*msb1**2*s**(-1)*u1**(-1)*u2 - 16*
     +    h_r*ct*sb*m1*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*sb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 - 64*h_r*ct*sb*m1*s**(-2)*t1*u2
     +     - 32*h_r*ct*sb*m1*s**(-2)*t1**2 - 32*h_r*ct*sb*m1*s**(-2)*
     +    u2**2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 128*h_r*ct*sb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*
     +    ct*sb*m1*s**(-1)*t1 - 64*h_r*ct*sb*m1*s**(-1)*t1**2*u1**(-1)
     +     - 64*h_r*ct*sb*m1*s**(-1)*u1**(-1)*u2**2 - 32*h_r*ct*sb*m1*
     +    s**(-1)*u2 - 64*h_r*ct*sb*m1*t1*u1**(-1) - 48*h_r*ct*sb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*cb*m1*mg**2*s**(-1)*t1*u1
     +     + 32*h_l*st*cb*m1*mg**2*s**(-1)*t1**2 + 32*h_l*st*cb*m1*
     +    mg**2*s**(-1)*u1**2 - 64*h_l*st*cb*m1*msb1**2*s**(-1)*t1*u1
     +     - 32*h_l*st*cb*m1*msb1**2*s**(-1)*t1**2 - 32*h_l*st*cb*m1*
     +    msb1**2*s**(-1)*u1**2 - 64*h_l*st*cb*m1*t1*u1 - 32*h_l*st*cb*
     +    m1*t1**2 - 32*h_l*st*cb*m1*u1**2 - 128*h_l*st*cb*m1**3*mg**2
     +     + 128*h_l*st*cb*m1**3*msb1**2 + 128*h_l*st*cb*m1**3*s - 64*
     +    h_r*ct*sb*m1*mg**2*s**(-1)*t1*u1 - 32*h_r*ct*sb*m1*mg**2*
     +    s**(-1)*t1**2 - 32*h_r*ct*sb*m1*mg**2*s**(-1)*u1**2 + 64*h_r*
     +    ct*sb*m1*msb1**2*s**(-1)*t1*u1 + 32*h_r*ct*sb*m1*msb1**2*
     +    s**(-1)*t1**2 + 32*h_r*ct*sb*m1*msb1**2*s**(-1)*u1**2 + 64*
     +    h_r*ct*sb*m1*t1*u1 + 32*h_r*ct*sb*m1*t1**2 + 32*h_r*ct*sb*m1*
     +    u1**2 + 128*h_r*ct*sb*m1**3*mg**2 - 128*h_r*ct*sb*m1**3*
     +    msb1**2 - 128*h_r*ct*sb*m1**3*s )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*st*cb*m1*mg**2*u1**(-1) - 32*h_l*st*cb*m1*
     +    mst2**2*u1**(-1) + 32*h_l*st*cb*m1*s*u1**(-1) + 16*h_l*st*cb*
     +    m1*t1*u1**(-1) + 80*h_l*st*cb*m1 + 32*h_l*st*cb*m1**3*
     +    u1**(-1) - 16*h_l*ct*cb*mg*s**(-1)*t1 - 16*h_l*ct*cb*mg*
     +    s**(-1)*u1 - 16*h_l*ct*cb*mg*t1*u1**(-1) - 16*h_l*ct*cb*mg + 
     +    16*h_r*st*sb*mg*s**(-1)*t1 + 16*h_r*st*sb*mg*s**(-1)*u1 + 16*
     +    h_r*st*sb*mg*t1*u1**(-1) + 16*h_r*st*sb*mg - 32*h_r*ct*sb*m1*
     +    mg**2*u1**(-1) + 32*h_r*ct*sb*m1*mst2**2*u1**(-1) - 32*h_r*ct
     +    *sb*m1*s*u1**(-1) - 16*h_r*ct*sb*m1*t1*u1**(-1) - 80*h_r*ct*
     +    sb*m1 - 32*h_r*ct*sb*m1**3*u1**(-1) )
      MM_s = MM_s + SCC(6,4)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*st*cb*m1*mg**2*s**(-1) - 32*h_l*st*cb*m1*
     +    mg**2*u1**(-1) + 32*h_l*st*cb*m1*msb1**2*s**(-1) + 32*h_l*st*
     +    cb*m1*mst2**2*u1**(-1) - 32*h_l*st*cb*m1*s*u1**(-1) - 32*h_l*
     +    st*cb*m1*t1*u1**(-1) - 32*h_l*st*cb*m1 - 32*h_l*st*cb*m1**3*
     +    u1**(-1) + 16*h_l*ct*cb*mg*s**(-1)*t1 + 16*h_l*ct*cb*mg*
     +    s**(-1)*u1 + 16*h_l*ct*cb*mg*t1*u1**(-1) + 16*h_l*ct*cb*mg - 
     +    16*h_r*st*sb*mg*s**(-1)*t1 - 16*h_r*st*sb*mg*s**(-1)*u1 - 16*
     +    h_r*st*sb*mg*t1*u1**(-1) - 16*h_r*st*sb*mg + 32*h_r*ct*sb*m1*
     +    mg**2*s**(-1) + 32*h_r*ct*sb*m1*mg**2*u1**(-1) - 32*h_r*ct*sb
     +    *m1*msb1**2*s**(-1) - 32*h_r*ct*sb*m1*mst2**2*u1**(-1) + 32*
     +    h_r*ct*sb*m1*s*u1**(-1) + 32*h_r*ct*sb*m1*t1*u1**(-1) + 32*
     +    h_r*ct*sb*m1 + 32*h_r*ct*sb*m1**3*u1**(-1) )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 96*h_l*st*sb*m1*m2**2*mg**2*
     +    s**(-2)*t1*u2 - 64*h_l*st*sb*m1*m2**2*mg**2*s**(-2)*t1**2 - 
     +    32*h_l*st*sb*m1*m2**2*mg**2*s**(-2)*u2**2 - 192*h_l*st*sb*m1*
     +    m2**2*mg**2*s**(-1)*t1*u1**(-1)*u2 - 128*h_l*st*sb*m1*m2**2*
     +    mg**2*s**(-1)*t1 - 128*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*t1**2
     +    *u1**(-1) - 64*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*
     +    u2**2 - 64*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_l*st*
     +    sb*m1*m2**2*mg**2*t1*u1**(-1) - 96*h_l*st*sb*m1*m2**2*mg**2*
     +    u1**(-1)*u2 + 32*h_l*st*sb*m1*m2**2*msb2**2*s**(-2)*t1*u2 + 
     +    32*h_l*st*sb*m1*m2**2*msb2**2*s**(-2)*t1**2 + 64*h_l*st*sb*m1
     +    *m2**2*msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*m2**2
     +    *msb2**2*s**(-1)*t1 + 64*h_l*st*sb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_l*st*sb*m1*m2**2*msb2**2*t1*u1**(-1) + 
     +    64*h_l*st*sb*m1*m2**2*mst2**2*s**(-2)*t1*u2 + 32*h_l*st*sb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*st*sb*m1*m2**2*mst2**2*s**(-2)
     +    *u2**2 + 128*h_l*st*sb*m1*m2**2*mst2**2*s**(-1)*t1*u1**(-1)*
     +    u2 + 64*h_l*st*sb*m1*m2**2*mst2**2*s**(-1)*t1 + 64*h_l*st*sb*
     +    m1*m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*sb*m1*
     +    m2**2*mst2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*sb*m1*m2**2*
     +    mst2**2*s**(-1)*u2 + 64*h_l*st*sb*m1*m2**2*mst2**2*t1*
     +    u1**(-1) + 96*h_l*st*sb*m1*m2**2*mst2**2*u1**(-1)*u2 + 96*h_l
     +    *st*sb*m1*m2**2*s**(-1)*t1*u2 + 32*h_l*st*sb*m1*m2**2*s**(-1)
     +    *t1**2 + 64*h_l*st*sb*m1*m2**2*s**(-1)*u2**2 + 80*h_l*st*sb*
     +    m1*m2**2*s*t1*u1**(-1) + 176*h_l*st*sb*m1*m2**2*s*u1**(-1)*u2
     +     + 192*h_l*st*sb*m1*m2**2*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*
     +    m2**2*t1 + 64*h_l*st*sb*m1*m2**2*t1**2*u1**(-1) + 128*h_l*st*
     +    sb*m1*m2**2*u1**(-1)*u2**2 + 128*h_l*st*sb*m1*m2**2*u2 + 32*
     +    h_l*st*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_l*st*sb*m1*m2**4*msb2**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*m1*m2**4*mst2**2*s**(-1)*
     +    t1*u1**(-1) - 16*h_l*st*sb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_l*st*sb*m1*m2**4*s**(-2)*t1*u2 - 32*h_l*st*sb*m1*
     +    m2**4*s**(-2)*t1**2 - 32*h_l*st*sb*m1*m2**4*s**(-2)*u2**2 - 
     +    128*h_l*st*sb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*
     +    m1*m2**4*s**(-1)*t1 - 64*h_l*st*sb*m1*m2**4*s**(-1)*t1**2*
     +    u1**(-1) - 64*h_l*st*sb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 - 64*
     +    h_l*st*sb*m1*m2**4*s**(-1)*u2 - 96*h_l*st*sb*m1*m2**4*t1*
     +    u1**(-1) - 112*h_l*st*sb*m1*m2**4*u1**(-1)*u2 + 16*h_l*st*sb*
     +    m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*m1*m2**6*s**(-1)*
     +    u1**(-1)*u2 + 32*h_l*st*sb*m1*mg**2*s**(-1)*t1*u2 + 32*h_l*st
     +    *sb*m1*mg**2*s**(-1)*u2**2 + 80*h_l*st*sb*m1*mg**2*s*u1**(-1)
     +    *u2 + 64*h_l*st*sb*m1*mg**2*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*
     +    mg**2*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_l*st*sb*m1*mg**2*u2 - 32*h_l*st*
     +    sb*m1*msb2**2*s**(-1)*t1*u2 - 32*h_l*st*sb*m1*msb2**2*s**(-1)
     +    *t1**2 - 80*h_l*st*sb*m1*msb2**2*s*t1*u1**(-1) - 64*h_l*st*sb
     +    *m1*msb2**2*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*msb2**2*t1 - 64*
     +    h_l*st*sb*m1*msb2**2*t1**2*u1**(-1) + 32*h_l*st*sb*m1*mst2**2
     +    *s**(-1)*t1**2 - 32*h_l*st*sb*m1*mst2**2*s**(-1)*u2**2 + 80*
     +    h_l*st*sb*m1*mst2**2*s*t1*u1**(-1) - 80*h_l*st*sb*m1*mst2**2*
     +    s*u1**(-1)*u2 + 64*h_l*st*sb*m1*mst2**2*t1 + 64*h_l*st*sb*m1*
     +    mst2**2*t1**2*u1**(-1) - 64*h_l*st*sb*m1*mst2**2*u1**(-1)*
     +    u2**2 - 64*h_l*st*sb*m1*mst2**2*u2 - 64*h_l*st*sb*m1*s*t1*
     +    u1**(-1)*u2 - 64*h_l*st*sb*m1*s*u1**(-1)*u2**2 - 64*h_l*st*sb
     +    *m1*s*u2 - 80*h_l*st*sb*m1*s**2*u1**(-1)*u2 - 32*h_l*st*sb*m1
     +    *t1*u2 - 32*h_l*st*sb*m1*u2**2 - 32*h_l*st*sb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) - 32*h_l*st*sb*m1**3*m2**2*msb2**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_l*st*sb*m1**3*m2**2*mst2**2*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1**3*m2**2*mst2**2*
     +    s**(-1)*u1**(-1)*u2 + 64*h_l*st*sb*m1**3*m2**2*s**(-2)*t1*u2
     +     + 32*h_l*st*sb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_l*st*sb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_l*st*sb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_l*st*sb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_l*st*sb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_l*st*sb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*sb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_l*st*sb*m1**3*m2**2*t1*u1**(-1) + 96*h_l*
     +    st*sb*m1**3*m2**2*u1**(-1)*u2 - 32*h_l*st*sb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_l*st*sb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*st*sb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_l
     +    *st*sb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_l*st*sb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1**3*mg**2*s**(-1)*
     +    u2 - 64*h_l*st*sb*m1**3*mg**2*u1**(-1)*u2 + 96*h_l*st*sb*
     +    m1**3*msb2**2*s**(-2)*t1*u2 + 32*h_l*st*sb*m1**3*msb2**2*
     +    s**(-2)*t1**2 + 64*h_l*st*sb*m1**3*msb2**2*s**(-2)*u2**2 + 
     +    192*h_l*st*sb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*
     +    st*sb*m1**3*msb2**2*s**(-1)*t1 + 64*h_l*st*sb*m1**3*msb2**2*
     +    s**(-1)*t1**2*u1**(-1) + 128*h_l*st*sb*m1**3*msb2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 128*h_l*st*sb*m1**3*msb2**2*s**(-1)*u2 + 64*
     +    h_l*st*sb*m1**3*msb2**2*t1*u1**(-1) + 160*h_l*st*sb*m1**3*
     +    msb2**2*u1**(-1)*u2 - 64*h_l*st*sb*m1**3*mst2**2*s**(-2)*t1*
     +    u2 - 32*h_l*st*sb*m1**3*mst2**2*s**(-2)*t1**2 - 32*h_l*st*sb*
     +    m1**3*mst2**2*s**(-2)*u2**2 - 128*h_l*st*sb*m1**3*mst2**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1**3*mst2**2*s**(-1)*
     +    t1 - 64*h_l*st*sb*m1**3*mst2**2*s**(-1)*t1**2*u1**(-1) - 64*
     +    h_l*st*sb*m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1**3*mst2**2*s**(-1)
     +    *u2 - 64*h_l*st*sb*m1**3*mst2**2*t1*u1**(-1) - 96*h_l*st*sb*
     +    m1**3*mst2**2*u1**(-1)*u2 + 32*h_l*st*sb*m1**3*s**(-1)*t1*u2
     +     + 32*h_l*st*sb*m1**3*s**(-1)*u2**2 + 64*h_l*st*sb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_l*st*sb*m1**3*t1*u1**(-1)*u2 + 64*h_l*st*
     +    sb*m1**3*u1**(-1)*u2**2 + 64*h_l*st*sb*m1**3*u2 + 16*h_l*st*
     +    sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_l*st*sb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_l*st*sb*m1**5*msb2**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_l*st*sb*m1**5*msb2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*st
     +    *sb*m1**5*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*m1**5*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*st*sb*m1**5*u1**(-1)*u2
     +     - 96*h_r*ct*cb*m1*m2**2*mg**2*s**(-2)*t1*u2 - 64*h_r*ct*cb*
     +    m1*m2**2*mg**2*s**(-2)*t1**2 - 32*h_r*ct*cb*m1*m2**2*mg**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 192*h_r*ct*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 128*h_r*ct*cb*m1*m2**2*mg**2*s**(-1)
     +    *t1 - 128*h_r*ct*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_r*ct*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*
     +    ct*cb*m1*m2**2*mg**2*s**(-1)*u2 - 160*h_r*ct*cb*m1*m2**2*
     +    mg**2*t1*u1**(-1) - 96*h_r*ct*cb*m1*m2**2*mg**2*u1**(-1)*u2
     +     + 32*h_r*ct*cb*m1*m2**2*msb2**2*s**(-2)*t1*u2 + 32*h_r*ct*cb
     +    *m1*m2**2*msb2**2*s**(-2)*t1**2 + 64*h_r*ct*cb*m1*m2**2*
     +    msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1*m2**2*
     +    msb2**2*s**(-1)*t1 + 64*h_r*ct*cb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 96*h_r*ct*cb*m1*m2**2*msb2**2*t1*u1**(-1) + 
     +    64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-2)*t1*u2 + 32*h_r*ct*cb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 + 32*h_r*ct*cb*m1*m2**2*mst2**2*
     +    s**(-2)*u2**2 + 128*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*t1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)
     +    *t1 + 64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r
     +    *ct*cb*m1*m2**2*mst2**2*s**(-1)*u2 + 64*h_r*ct*cb*m1*m2**2*
     +    mst2**2*t1*u1**(-1) + 96*h_r*ct*cb*m1*m2**2*mst2**2*u1**(-1)*
     +    u2 + 96*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u2 + 32*h_r*ct*cb*m1*
     +    m2**2*s**(-1)*t1**2 + 64*h_r*ct*cb*m1*m2**2*s**(-1)*u2**2 + 
     +    80*h_r*ct*cb*m1*m2**2*s*t1*u1**(-1) + 176*h_r*ct*cb*m1*m2**2*
     +    s*u1**(-1)*u2 + 192*h_r*ct*cb*m1*m2**2*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*cb*m1*m2**2*t1 + 64*h_r*ct*cb*m1*m2**2*t1**2*u1**(-1)
     +     + 128*h_r*ct*cb*m1*m2**2*u1**(-1)*u2**2 + 128*h_r*ct*cb*m1*
     +    m2**2*u2 + 32*h_r*ct*cb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_r*ct*cb*m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 - 16*h_r*ct*
     +    cb*m1*m2**4*msb2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1*
     +    m2**4*mst2**2*s**(-1)*t1*u1**(-1) )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 16*h_r*ct*cb*m1*m2**4*mst2**2*
     +    s**(-1)*u1**(-1)*u2 - 64*h_r*ct*cb*m1*m2**4*s**(-2)*t1*u2 - 
     +    32*h_r*ct*cb*m1*m2**4*s**(-2)*t1**2 - 32*h_r*ct*cb*m1*m2**4*
     +    s**(-2)*u2**2 - 128*h_r*ct*cb*m1*m2**4*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_r*ct*cb*m1*m2**4*s**(-1)*t1 - 64*h_r*ct*cb*m1*m2**4*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*cb*m1*m2**4*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*cb*m1*m2**4*s**(-1)*u2 - 96*h_r*ct
     +    *cb*m1*m2**4*t1*u1**(-1) - 112*h_r*ct*cb*m1*m2**4*u1**(-1)*u2
     +     + 16*h_r*ct*cb*m1*m2**6*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*
     +    m1*m2**6*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m1*mg**2*s**(-1)*
     +    t1*u2 + 32*h_r*ct*cb*m1*mg**2*s**(-1)*u2**2 + 80*h_r*ct*cb*m1
     +    *mg**2*s*u1**(-1)*u2 + 64*h_r*ct*cb*m1*mg**2*t1*u1**(-1)*u2
     +     + 64*h_r*ct*cb*m1*mg**2*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1*
     +    mg**2*u2 - 32*h_r*ct*cb*m1*msb2**2*s**(-1)*t1*u2 - 32*h_r*ct*
     +    cb*m1*msb2**2*s**(-1)*t1**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 80*h_r*ct*cb*m1*msb2**2*s*t1*
     +    u1**(-1) - 64*h_r*ct*cb*m1*msb2**2*t1*u1**(-1)*u2 - 64*h_r*ct
     +    *cb*m1*msb2**2*t1 - 64*h_r*ct*cb*m1*msb2**2*t1**2*u1**(-1) + 
     +    32*h_r*ct*cb*m1*mst2**2*s**(-1)*t1**2 - 32*h_r*ct*cb*m1*
     +    mst2**2*s**(-1)*u2**2 + 80*h_r*ct*cb*m1*mst2**2*s*t1*u1**(-1)
     +     - 80*h_r*ct*cb*m1*mst2**2*s*u1**(-1)*u2 + 64*h_r*ct*cb*m1*
     +    mst2**2*t1 + 64*h_r*ct*cb*m1*mst2**2*t1**2*u1**(-1) - 64*h_r*
     +    ct*cb*m1*mst2**2*u1**(-1)*u2**2 - 64*h_r*ct*cb*m1*mst2**2*u2
     +     - 64*h_r*ct*cb*m1*s*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1*s*
     +    u1**(-1)*u2**2 - 64*h_r*ct*cb*m1*s*u2 - 80*h_r*ct*cb*m1*s**2*
     +    u1**(-1)*u2 - 32*h_r*ct*cb*m1*t1*u2 - 32*h_r*ct*cb*m1*u2**2
     +     - 32*h_r*ct*cb*m1**3*m2**2*mg**2*s**(-1)*t1*u1**(-1) - 32*
     +    h_r*ct*cb*m1**3*m2**2*msb2**2*s**(-1)*u1**(-1)*u2 + 32*h_r*ct
     +    *cb*m1**3*m2**2*mst2**2*s**(-1)*t1*u1**(-1) + 32*h_r*ct*cb*
     +    m1**3*m2**2*mst2**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 64*h_r*ct*cb*m1**3*m2**2*s**(-2)*t1*
     +    u2 + 32*h_r*ct*cb*m1**3*m2**2*s**(-2)*t1**2 + 32*h_r*ct*cb*
     +    m1**3*m2**2*s**(-2)*u2**2 + 128*h_r*ct*cb*m1**3*m2**2*s**(-1)
     +    *t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1**3*m2**2*s**(-1)*t1 + 64*
     +    h_r*ct*cb*m1**3*m2**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*cb*
     +    m1**3*m2**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1**3*m2**2
     +    *s**(-1)*u2 + 96*h_r*ct*cb*m1**3*m2**2*t1*u1**(-1) + 96*h_r*
     +    ct*cb*m1**3*m2**2*u1**(-1)*u2 - 32*h_r*ct*cb*m1**3*m2**4*
     +    s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb*m1**3*m2**4*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*ct*cb*m1**3*mg**2*s**(-2)*t1*u2 - 32*h_r
     +    *ct*cb*m1**3*mg**2*s**(-2)*u2**2 - 64*h_r*ct*cb*m1**3*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*mg**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*cb*m1**3*mg**2*s**(-1)*u2 - 64*h_r
     +    *ct*cb*m1**3*mg**2*u1**(-1)*u2 + 96*h_r*ct*cb*m1**3*msb2**2*
     +    s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*ct*cb*m1**3*msb2**2*s**(-2)*
     +    t1**2 + 64*h_r*ct*cb*m1**3*msb2**2*s**(-2)*u2**2 + 192*h_r*ct
     +    *cb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1**3
     +    *msb2**2*s**(-1)*t1 + 64*h_r*ct*cb*m1**3*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) + 128*h_r*ct*cb*m1**3*msb2**2*s**(-1)*u1**(-1)
     +    *u2**2 + 128*h_r*ct*cb*m1**3*msb2**2*s**(-1)*u2 + 64*h_r*ct*
     +    cb*m1**3*msb2**2*t1*u1**(-1) + 160*h_r*ct*cb*m1**3*msb2**2*
     +    u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*mst2**2*s**(-2)*t1*u2 - 32*
     +    h_r*ct*cb*m1**3*mst2**2*s**(-2)*t1**2 - 32*h_r*ct*cb*m1**3*
     +    mst2**2*s**(-2)*u2**2 - 128*h_r*ct*cb*m1**3*mst2**2*s**(-1)*
     +    t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*mst2**2*s**(-1)*t1 - 64*
     +    h_r*ct*cb*m1**3*mst2**2*s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*cb
     +    *m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_r*ct*cb*m1**3*
     +    mst2**2*s**(-1)*u2 - 64*h_r*ct*cb*m1**3*mst2**2*t1*u1**(-1)
     +     - 96*h_r*ct*cb*m1**3*mst2**2*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * ( 32*h_r*ct*cb*m1**3*s**(-1)*t1*u2 + 32
     +    *h_r*ct*cb*m1**3*s**(-1)*u2**2 + 64*h_r*ct*cb*m1**3*s*
     +    u1**(-1)*u2 + 64*h_r*ct*cb*m1**3*t1*u1**(-1)*u2 + 64*h_r*ct*
     +    cb*m1**3*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1**3*u2 + 16*h_r*ct*
     +    cb*m1**5*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*m1**5*m2**2
     +    *s**(-1)*u1**(-1)*u2 - 16*h_r*ct*cb*m1**5*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*cb*m1**5*msb2**2*s**(-1)*t1*u1**(-1)
     +     + 32*h_r*ct*cb*m1**5*msb2**2*s**(-1)*u1**(-1)*u2 - 16*h_r*ct
     +    *cb*m1**5*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1**5*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*cb*m1**5*u1**(-1)*u2
     +     )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 96*h_l*st*sb*m1*m2**2*mg**2*s**(-2)*t1*u2
     +     + 64*h_l*st*sb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_l*st*sb*
     +    m1*m2**2*mg**2*s**(-2)*u2**2 + 192*h_l*st*sb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 + 128*h_l*st*sb*m1*m2**2*mg**2*s**(-1)
     +    *t1 + 128*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_l*st*sb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_l*
     +    st*sb*m1*m2**2*mg**2*s**(-1)*u2 + 160*h_l*st*sb*m1*m2**2*
     +    mg**2*t1*u1**(-1) + 96*h_l*st*sb*m1*m2**2*mg**2*u1**(-1)*u2
     +     - 32*h_l*st*sb*m1*m2**2*msb2**2*s**(-2)*t1*u2 - 32*h_l*st*sb
     +    *m1*m2**2*msb2**2*s**(-2)*t1**2 - 64*h_l*st*sb*m1*m2**2*
     +    msb2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*m2**2*
     +    msb2**2*s**(-1)*t1 - 64*h_l*st*sb*m1*m2**2*msb2**2*s**(-1)*
     +    t1**2*u1**(-1) - 96*h_l*st*sb*m1*m2**2*msb2**2*t1*u1**(-1) - 
     +    64*h_l*st*sb*m1*m2**2*mst2**2*s**(-2)*t1*u2 - 32*h_l*st*sb*m1
     +    *m2**2*mst2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*st*sb*m1*m2**2*mst2**2*s**(-2)*
     +    u2**2 - 128*h_l*st*sb*m1*m2**2*mst2**2*s**(-1)*t1*u1**(-1)*u2
     +     - 64*h_l*st*sb*m1*m2**2*mst2**2*s**(-1)*t1 - 64*h_l*st*sb*m1
     +    *m2**2*mst2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*sb*m1*m2**2
     +    *mst2**2*s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*sb*m1*m2**2*
     +    mst2**2*s**(-1)*u2 - 64*h_l*st*sb*m1*m2**2*mst2**2*t1*
     +    u1**(-1) - 96*h_l*st*sb*m1*m2**2*mst2**2*u1**(-1)*u2 - 96*h_l
     +    *st*sb*m1*m2**2*s**(-1)*t1*u2 - 32*h_l*st*sb*m1*m2**2*s**(-1)
     +    *t1**2 - 64*h_l*st*sb*m1*m2**2*s**(-1)*u2**2 - 80*h_l*st*sb*
     +    m1*m2**2*s*t1*u1**(-1) - 176*h_l*st*sb*m1*m2**2*s*u1**(-1)*u2
     +     - 192*h_l*st*sb*m1*m2**2*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*
     +    m2**2*t1 - 64*h_l*st*sb*m1*m2**2*t1**2*u1**(-1) - 128*h_l*st*
     +    sb*m1*m2**2*u1**(-1)*u2**2 - 128*h_l*st*sb*m1*m2**2*u2 - 32*
     +    h_l*st*sb*m1*m2**4*mg**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*
     +    m1*m2**4*mg**2*s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 16*h_l*st*sb*m1*m2**4*msb2**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_l*st*sb*m1*m2**4*mst2**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_l*st*sb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_l
     +    *st*sb*m1*m2**4*s**(-2)*t1*u2 + 32*h_l*st*sb*m1*m2**4*s**(-2)
     +    *t1**2 + 32*h_l*st*sb*m1*m2**4*s**(-2)*u2**2 + 128*h_l*st*sb*
     +    m1*m2**4*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*m2**4*
     +    s**(-1)*t1 + 64*h_l*st*sb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_l*st*sb*m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_l*st*sb*
     +    m1*m2**4*s**(-1)*u2 + 96*h_l*st*sb*m1*m2**4*t1*u1**(-1) + 112
     +    *h_l*st*sb*m1*m2**4*u1**(-1)*u2 - 16*h_l*st*sb*m1*m2**6*
     +    s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*m1*m2**6*s**(-1)*u1**(-1)*
     +    u2 - 32*h_l*st*sb*m1*mg**2*s**(-1)*t1*u2 - 32*h_l*st*sb*m1*
     +    mg**2*s**(-1)*u2**2 - 80*h_l*st*sb*m1*mg**2*s*u1**(-1)*u2 - 
     +    64*h_l*st*sb*m1*mg**2*t1*u1**(-1)*u2 - 64*h_l*st*sb*m1*mg**2*
     +    u1**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_l*st*sb*m1*mg**2*u2 + 32*h_l*st*sb*
     +    m1*msb2**2*s**(-1)*t1*u2 + 32*h_l*st*sb*m1*msb2**2*s**(-1)*
     +    t1**2 + 80*h_l*st*sb*m1*msb2**2*s*t1*u1**(-1) + 64*h_l*st*sb*
     +    m1*msb2**2*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1*msb2**2*t1 + 64*
     +    h_l*st*sb*m1*msb2**2*t1**2*u1**(-1) - 32*h_l*st*sb*m1*mst2**2
     +    *s**(-1)*t1**2 + 32*h_l*st*sb*m1*mst2**2*s**(-1)*u2**2 - 80*
     +    h_l*st*sb*m1*mst2**2*s*t1*u1**(-1) + 80*h_l*st*sb*m1*mst2**2*
     +    s*u1**(-1)*u2 - 64*h_l*st*sb*m1*mst2**2*t1 - 64*h_l*st*sb*m1*
     +    mst2**2*t1**2*u1**(-1) + 64*h_l*st*sb*m1*mst2**2*u1**(-1)*
     +    u2**2 + 64*h_l*st*sb*m1*mst2**2*u2 + 64*h_l*st*sb*m1*s*t1*
     +    u1**(-1)*u2 + 64*h_l*st*sb*m1*s*u1**(-1)*u2**2 + 64*h_l*st*sb
     +    *m1*s*u2 + 80*h_l*st*sb*m1*s**2*u1**(-1)*u2 + 32*h_l*st*sb*m1
     +    *t1*u2 + 32*h_l*st*sb*m1*u2**2 + 32*h_l*st*sb*m1**3*m2**2*
     +    mg**2*s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1**3*m2**2*msb2**2*
     +    s**(-1)*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_l*st*sb*m1**3*m2**2*mst2**2*s**(-1)
     +    *t1*u1**(-1) - 32*h_l*st*sb*m1**3*m2**2*mst2**2*s**(-1)*
     +    u1**(-1)*u2 - 64*h_l*st*sb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_l
     +    *st*sb*m1**3*m2**2*s**(-2)*t1**2 - 32*h_l*st*sb*m1**3*m2**2*
     +    s**(-2)*u2**2 - 128*h_l*st*sb*m1**3*m2**2*s**(-1)*t1*u1**(-1)
     +    *u2 - 64*h_l*st*sb*m1**3*m2**2*s**(-1)*t1 - 64*h_l*st*sb*
     +    m1**3*m2**2*s**(-1)*t1**2*u1**(-1) - 64*h_l*st*sb*m1**3*m2**2
     +    *s**(-1)*u1**(-1)*u2**2 - 64*h_l*st*sb*m1**3*m2**2*s**(-1)*u2
     +     - 96*h_l*st*sb*m1**3*m2**2*t1*u1**(-1) - 96*h_l*st*sb*m1**3*
     +    m2**2*u1**(-1)*u2 + 32*h_l*st*sb*m1**3*m2**4*s**(-1)*t1*
     +    u1**(-1) + 32*h_l*st*sb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*
     +    h_l*st*sb*m1**3*mg**2*s**(-2)*t1*u2 + 32*h_l*st*sb*m1**3*
     +    mg**2*s**(-2)*u2**2 + 64*h_l*st*sb*m1**3*mg**2*s**(-1)*t1*
     +    u1**(-1)*u2 + 64*h_l*st*sb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*st*sb*m1**3*mg**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*sb*m1**3*mg**2*u1**(-1)*u2 - 96*
     +    h_l*st*sb*m1**3*msb2**2*s**(-2)*t1*u2 - 32*h_l*st*sb*m1**3*
     +    msb2**2*s**(-2)*t1**2 - 64*h_l*st*sb*m1**3*msb2**2*s**(-2)*
     +    u2**2 - 192*h_l*st*sb*m1**3*msb2**2*s**(-1)*t1*u1**(-1)*u2 - 
     +    64*h_l*st*sb*m1**3*msb2**2*s**(-1)*t1 - 64*h_l*st*sb*m1**3*
     +    msb2**2*s**(-1)*t1**2*u1**(-1) - 128*h_l*st*sb*m1**3*msb2**2*
     +    s**(-1)*u1**(-1)*u2**2 - 128*h_l*st*sb*m1**3*msb2**2*s**(-1)*
     +    u2 - 64*h_l*st*sb*m1**3*msb2**2*t1*u1**(-1) - 160*h_l*st*sb*
     +    m1**3*msb2**2*u1**(-1)*u2 + 64*h_l*st*sb*m1**3*mst2**2*
     +    s**(-2)*t1*u2 + 32*h_l*st*sb*m1**3*mst2**2*s**(-2)*t1**2 + 32
     +    *h_l*st*sb*m1**3*mst2**2*s**(-2)*u2**2 + 128*h_l*st*sb*m1**3*
     +    mst2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_l*st*sb*m1**3*mst2**2*
     +    s**(-1)*t1 + 64*h_l*st*sb*m1**3*mst2**2*s**(-1)*t1**2*
     +    u1**(-1) + 64*h_l*st*sb*m1**3*mst2**2*s**(-1)*u1**(-1)*u2**2
     +     + 64*h_l*st*sb*m1**3*mst2**2*s**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_l*st*sb*m1**3*mst2**2*t1*u1**(-1) + 96
     +    *h_l*st*sb*m1**3*mst2**2*u1**(-1)*u2 - 32*h_l*st*sb*m1**3*
     +    s**(-1)*t1*u2 - 32*h_l*st*sb*m1**3*s**(-1)*u2**2 - 64*h_l*st*
     +    sb*m1**3*s*u1**(-1)*u2 - 64*h_l*st*sb*m1**3*t1*u1**(-1)*u2 - 
     +    64*h_l*st*sb*m1**3*u1**(-1)*u2**2 - 64*h_l*st*sb*m1**3*u2 - 
     +    16*h_l*st*sb*m1**5*m2**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*
     +    m1**5*m2**2*s**(-1)*u1**(-1)*u2 + 16*h_l*st*sb*m1**5*mg**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_l*st*sb*m1**5*msb2**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_l*st*sb*m1**5*msb2**2*s**(-1)*u1**(-1)*u2 + 
     +    16*h_l*st*sb*m1**5*mst2**2*s**(-1)*t1*u1**(-1) + 16*h_l*st*sb
     +    *m1**5*mst2**2*s**(-1)*u1**(-1)*u2 - 16*h_l*st*sb*m1**5*
     +    u1**(-1)*u2 + 96*h_r*ct*cb*m1*m2**2*mg**2*s**(-2)*t1*u2 + 64*
     +    h_r*ct*cb*m1*m2**2*mg**2*s**(-2)*t1**2 + 32*h_r*ct*cb*m1*
     +    m2**2*mg**2*s**(-2)*u2**2 + 192*h_r*ct*cb*m1*m2**2*mg**2*
     +    s**(-1)*t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 128*h_r*ct*cb*m1*m2**2*mg**2*s**(-1)*t1 + 
     +    128*h_r*ct*cb*m1*m2**2*mg**2*s**(-1)*t1**2*u1**(-1) + 64*h_r*
     +    ct*cb*m1*m2**2*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1
     +    *m2**2*mg**2*s**(-1)*u2 + 160*h_r*ct*cb*m1*m2**2*mg**2*t1*
     +    u1**(-1) + 96*h_r*ct*cb*m1*m2**2*mg**2*u1**(-1)*u2 - 32*h_r*
     +    ct*cb*m1*m2**2*msb2**2*s**(-2)*t1*u2 - 32*h_r*ct*cb*m1*m2**2*
     +    msb2**2*s**(-2)*t1**2 - 64*h_r*ct*cb*m1*m2**2*msb2**2*s**(-1)
     +    *t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1*m2**2*msb2**2*s**(-1)*t1 - 
     +    64*h_r*ct*cb*m1*m2**2*msb2**2*s**(-1)*t1**2*u1**(-1) - 96*h_r
     +    *ct*cb*m1*m2**2*msb2**2*t1*u1**(-1) - 64*h_r*ct*cb*m1*m2**2*
     +    mst2**2*s**(-2)*t1*u2 - 32*h_r*ct*cb*m1*m2**2*mst2**2*s**(-2)
     +    *t1**2 - 32*h_r*ct*cb*m1*m2**2*mst2**2*s**(-2)*u2**2 - 128*
     +    h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*ct
     +    *cb*m1*m2**2*mst2**2*s**(-1)*t1 - 64*h_r*ct*cb*m1*m2**2*
     +    mst2**2*s**(-1)*t1**2*u1**(-1) )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*cb*m1*m2**2*mst2**2*s**(-1)*u2 - 
     +    64*h_r*ct*cb*m1*m2**2*mst2**2*t1*u1**(-1) - 96*h_r*ct*cb*m1*
     +    m2**2*mst2**2*u1**(-1)*u2 - 96*h_r*ct*cb*m1*m2**2*s**(-1)*t1*
     +    u2 - 32*h_r*ct*cb*m1*m2**2*s**(-1)*t1**2 - 64*h_r*ct*cb*m1*
     +    m2**2*s**(-1)*u2**2 - 80*h_r*ct*cb*m1*m2**2*s*t1*u1**(-1) - 
     +    176*h_r*ct*cb*m1*m2**2*s*u1**(-1)*u2 - 192*h_r*ct*cb*m1*m2**2
     +    *t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1*m2**2*t1 - 64*h_r*ct*cb*m1*
     +    m2**2*t1**2*u1**(-1) - 128*h_r*ct*cb*m1*m2**2*u1**(-1)*u2**2
     +     - 128*h_r*ct*cb*m1*m2**2*u2 - 32*h_r*ct*cb*m1*m2**4*mg**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1*m2**4*mg**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*cb*m1*m2**4*msb2**2*s**(-1)*t1*
     +    u1**(-1) + 16*h_r*ct*cb*m1*m2**4*mst2**2*s**(-1)*t1*u1**(-1)
     +     + 16*h_r*ct*cb*m1*m2**4*mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_r
     +    *ct*cb*m1*m2**4*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 32*h_r*ct*cb*m1*m2**4*s**(-2)*t1**2 + 32*
     +    h_r*ct*cb*m1*m2**4*s**(-2)*u2**2 + 128*h_r*ct*cb*m1*m2**4*
     +    s**(-1)*t1*u1**(-1)*u2 + 64*h_r*ct*cb*m1*m2**4*s**(-1)*t1 + 
     +    64*h_r*ct*cb*m1*m2**4*s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*cb*
     +    m1*m2**4*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1*m2**4*
     +    s**(-1)*u2 + 96*h_r*ct*cb*m1*m2**4*t1*u1**(-1) + 112*h_r*ct*
     +    cb*m1*m2**4*u1**(-1)*u2 - 16*h_r*ct*cb*m1*m2**6*s**(-1)*t1*
     +    u1**(-1) - 16*h_r*ct*cb*m1*m2**6*s**(-1)*u1**(-1)*u2 - 32*h_r
     +    *ct*cb*m1*mg**2*s**(-1)*t1*u2 - 32*h_r*ct*cb*m1*mg**2*s**(-1)
     +    *u2**2 - 80*h_r*ct*cb*m1*mg**2*s*u1**(-1)*u2 - 64*h_r*ct*cb*
     +    m1*mg**2*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1*mg**2*u1**(-1)*
     +    u2**2 - 64*h_r*ct*cb*m1*mg**2*u2 + 32*h_r*ct*cb*m1*msb2**2*
     +    s**(-1)*t1*u2 + 32*h_r*ct*cb*m1*msb2**2*s**(-1)*t1**2 + 80*
     +    h_r*ct*cb*m1*msb2**2*s*t1*u1**(-1) + 64*h_r*ct*cb*m1*msb2**2*
     +    t1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * ( 64*h_r*ct*cb*m1*msb2**2*t1 + 64*h_r*ct*cb*
     +    m1*msb2**2*t1**2*u1**(-1) - 32*h_r*ct*cb*m1*mst2**2*s**(-1)*
     +    t1**2 + 32*h_r*ct*cb*m1*mst2**2*s**(-1)*u2**2 - 80*h_r*ct*cb*
     +    m1*mst2**2*s*t1*u1**(-1) + 80*h_r*ct*cb*m1*mst2**2*s*u1**(-1)
     +    *u2 - 64*h_r*ct*cb*m1*mst2**2*t1 - 64*h_r*ct*cb*m1*mst2**2*
     +    t1**2*u1**(-1) + 64*h_r*ct*cb*m1*mst2**2*u1**(-1)*u2**2 + 64*
     +    h_r*ct*cb*m1*mst2**2*u2 + 64*h_r*ct*cb*m1*s*t1*u1**(-1)*u2 + 
     +    64*h_r*ct*cb*m1*s*u1**(-1)*u2**2 + 64*h_r*ct*cb*m1*s*u2 + 80*
     +    h_r*ct*cb*m1*s**2*u1**(-1)*u2 + 32*h_r*ct*cb*m1*t1*u2 + 32*
     +    h_r*ct*cb*m1*u2**2 + 32*h_r*ct*cb*m1**3*m2**2*mg**2*s**(-1)*
     +    t1*u1**(-1) + 32*h_r*ct*cb*m1**3*m2**2*msb2**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_r*ct*cb*m1**3*m2**2*mst2**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*ct*cb*m1**3*m2**2*mst2**2*s**(-1)*u1**(-1)*
     +    u2 - 64*h_r*ct*cb*m1**3*m2**2*s**(-2)*t1*u2 - 32*h_r*ct*cb*
     +    m1**3*m2**2*s**(-2)*t1**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 32*h_r*ct*cb*m1**3*m2**2*s**(-2)*u2**2
     +     - 128*h_r*ct*cb*m1**3*m2**2*s**(-1)*t1*u1**(-1)*u2 - 64*h_r*
     +    ct*cb*m1**3*m2**2*s**(-1)*t1 - 64*h_r*ct*cb*m1**3*m2**2*
     +    s**(-1)*t1**2*u1**(-1) - 64*h_r*ct*cb*m1**3*m2**2*s**(-1)*
     +    u1**(-1)*u2**2 - 64*h_r*ct*cb*m1**3*m2**2*s**(-1)*u2 - 96*h_r
     +    *ct*cb*m1**3*m2**2*t1*u1**(-1) - 96*h_r*ct*cb*m1**3*m2**2*
     +    u1**(-1)*u2 + 32*h_r*ct*cb*m1**3*m2**4*s**(-1)*t1*u1**(-1) + 
     +    32*h_r*ct*cb*m1**3*m2**4*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*
     +    m1**3*mg**2*s**(-2)*t1*u2 + 32*h_r*ct*cb*m1**3*mg**2*s**(-2)*
     +    u2**2 + 64*h_r*ct*cb*m1**3*mg**2*s**(-1)*t1*u1**(-1)*u2 + 64*
     +    h_r*ct*cb*m1**3*mg**2*s**(-1)*u1**(-1)*u2**2 + 64*h_r*ct*cb*
     +    m1**3*mg**2*s**(-1)*u2 + 64*h_r*ct*cb*m1**3*mg**2*u1**(-1)*u2
     +     - 96*h_r*ct*cb*m1**3*msb2**2*s**(-2)*t1*u2 - 32*h_r*ct*cb*
     +    m1**3*msb2**2*s**(-2)*t1**2 - 64*h_r*ct*cb*m1**3*msb2**2*
     +    s**(-2)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 192*h_r*ct*cb*m1**3*msb2**2*s**(-1)*t1*
     +    u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*msb2**2*s**(-1)*t1 - 64*h_r*
     +    ct*cb*m1**3*msb2**2*s**(-1)*t1**2*u1**(-1) - 128*h_r*ct*cb*
     +    m1**3*msb2**2*s**(-1)*u1**(-1)*u2**2 - 128*h_r*ct*cb*m1**3*
     +    msb2**2*s**(-1)*u2 - 64*h_r*ct*cb*m1**3*msb2**2*t1*u1**(-1)
     +     - 160*h_r*ct*cb*m1**3*msb2**2*u1**(-1)*u2 + 64*h_r*ct*cb*
     +    m1**3*mst2**2*s**(-2)*t1*u2 + 32*h_r*ct*cb*m1**3*mst2**2*
     +    s**(-2)*t1**2 + 32*h_r*ct*cb*m1**3*mst2**2*s**(-2)*u2**2 + 
     +    128*h_r*ct*cb*m1**3*mst2**2*s**(-1)*t1*u1**(-1)*u2 + 64*h_r*
     +    ct*cb*m1**3*mst2**2*s**(-1)*t1 + 64*h_r*ct*cb*m1**3*mst2**2*
     +    s**(-1)*t1**2*u1**(-1) + 64*h_r*ct*cb*m1**3*mst2**2*s**(-1)*
     +    u1**(-1)*u2**2 + 64*h_r*ct*cb*m1**3*mst2**2*s**(-1)*u2 + 64*
     +    h_r*ct*cb*m1**3*mst2**2*t1*u1**(-1) + 96*h_r*ct*cb*m1**3*
     +    mst2**2*u1**(-1)*u2 - 32*h_r*ct*cb*m1**3*s**(-1)*t1*u2 - 32*
     +    h_r*ct*cb*m1**3*s**(-1)*u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*kaellen**(-1)*Cf*Pi**2*
     + alphas**2*prefac * (  - 64*h_r*ct*cb*m1**3*s*u1**(-1)*u2 - 64*
     +    h_r*ct*cb*m1**3*t1*u1**(-1)*u2 - 64*h_r*ct*cb*m1**3*u1**(-1)*
     +    u2**2 - 64*h_r*ct*cb*m1**3*u2 - 16*h_r*ct*cb*m1**5*m2**2*
     +    s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1**5*m2**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_r*ct*cb*m1**5*mg**2*s**(-1)*u1**(-1)*u2 - 
     +    16*h_r*ct*cb*m1**5*msb2**2*s**(-1)*t1*u1**(-1) - 32*h_r*ct*cb
     +    *m1**5*msb2**2*s**(-1)*u1**(-1)*u2 + 16*h_r*ct*cb*m1**5*
     +    mst2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*m1**5*mst2**2*
     +    s**(-1)*u1**(-1)*u2 - 16*h_r*ct*cb*m1**5*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*
     +    h_l*st*sb*m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_l*st*sb*m1*
     +    mg**2*s**(-2)*t1 - 32*h_l*st*sb*m1*mg**2*s**(-2)*u2 - 32*h_l*
     +    st*sb*m1*mg**2*s**(-1)*t1*u1**(-1) - 48*h_l*st*sb*m1*mg**2*
     +    s**(-1)*u1**(-1)*u2 + 32*h_l*st*sb*m1*msb2**2*s**(-2)*t1 + 32
     +    *h_l*st*sb*m1*msb2**2*s**(-2)*u2 + 16*h_l*st*sb*m1*msb2**2*
     +    s**(-1)*t1*u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s**(-1)*
     +    u1**(-1)*u2 + 16*h_l*st*sb*m1*mst2**2*s**(-1)*t1*u1**(-1) + 
     +    16*h_l*st*sb*m1*mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_l*st*sb*m1
     +    *s**(-2)*t1*u2 + 32*h_l*st*sb*m1*s**(-2)*t1**2 + 32*h_l*st*sb
     +    *m1*s**(-2)*u2**2 + 128*h_l*st*sb*m1*s**(-1)*t1*u1**(-1)*u2
     +     + 32*h_l*st*sb*m1*s**(-1)*t1 + 64*h_l*st*sb*m1*s**(-1)*t1**2
     +    *u1**(-1) + 64*h_l*st*sb*m1*s**(-1)*u1**(-1)*u2**2 + 32*h_l*
     +    st*sb*m1*s**(-1)*u2 + 64*h_l*st*sb*m1*t1*u1**(-1) + 48*h_l*st
     +    *sb*m1*u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*sb*m1**2*mg*s**(-2)*t1 - 32*h_l*ct*sb*
     +    m1**2*mg*s**(-2)*u2 - 32*h_l*ct*sb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_l*ct*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l
     +    *ct*sb*m2**2*mg*s**(-2)*t1 + 32*h_l*ct*sb*m2**2*mg*s**(-2)*u2
     +     + 32*h_l*ct*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_l*ct*sb*
     +    mg*s**(-2)*t1*u2 + 32*h_l*ct*sb*mg*s**(-2)*t1**2 + 64*h_l*ct*
     +    sb*mg*s**(-1)*t1*u1**(-1)*u2 + 32*h_l*ct*sb*mg*s**(-1)*t1 + 
     +    64*h_l*ct*sb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_l*ct*sb*mg*
     +    s**(-1)*u2 + 64*h_l*ct*sb*mg*t1*u1**(-1) - 32*h_l*ct*sb*mg*
     +    u1**(-1)*u2 - 32*h_r*st*cb*m1**2*mg*s**(-2)*t1 - 32*h_r*st*cb
     +    *m1**2*mg*s**(-2)*u2 - 32*h_r*st*cb*m1**2*mg*s**(-1)*t1*
     +    u1**(-1) - 64*h_r*st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r
     +    *st*cb*m2**2*mg*s**(-2)*t1 + 32*h_r*st*cb*m2**2*mg*s**(-2)*u2
     +     + 32*h_r*st*cb*m2**2*mg*s**(-1)*u1**(-1)*u2 + 32*h_r*st*cb*
     +    mg*s**(-2)*t1*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_r*st*cb*mg*s**(-2)*t1**2 + 64*h_r*st*cb*mg*
     +    s**(-1)*t1*u1**(-1)*u2 + 32*h_r*st*cb*mg*s**(-1)*t1 + 64*h_r*
     +    st*cb*mg*s**(-1)*t1**2*u1**(-1) - 32*h_r*st*cb*mg*s**(-1)*u2
     +     + 64*h_r*st*cb*mg*t1*u1**(-1) - 32*h_r*st*cb*mg*u1**(-1)*u2
     +     - 16*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*
     +    m1*m2**2*s**(-1)*u1**(-1)*u2 - 32*h_r*ct*cb*m1*mg**2*s**(-2)*
     +    t1 - 32*h_r*ct*cb*m1*mg**2*s**(-2)*u2 - 32*h_r*ct*cb*m1*mg**2
     +    *s**(-1)*t1*u1**(-1) - 48*h_r*ct*cb*m1*mg**2*s**(-1)*u1**(-1)
     +    *u2 + 32*h_r*ct*cb*m1*msb2**2*s**(-2)*t1 + 32*h_r*ct*cb*m1*
     +    msb2**2*s**(-2)*u2 + 16*h_r*ct*cb*m1*msb2**2*s**(-1)*t1*
     +    u1**(-1) + 32*h_r*ct*cb*m1*msb2**2*s**(-1)*u1**(-1)*u2 + 16*
     +    h_r*ct*cb*m1*mst2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 + 64*h_r*ct*cb*m1*s**(-2)*t1*u2
     +     + 32*h_r*ct*cb*m1*s**(-2)*t1**2 + 32*h_r*ct*cb*m1*s**(-2)*
     +    u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 128*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1)*u2 + 32*h_r*ct*
     +    cb*m1*s**(-1)*t1 + 64*h_r*ct*cb*m1*s**(-1)*t1**2*u1**(-1) + 
     +    64*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2**2 + 32*h_r*ct*cb*m1*
     +    s**(-1)*u2 + 64*h_r*ct*cb*m1*t1*u1**(-1) + 48*h_r*ct*cb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*sb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_l*st
     +    *sb*m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_l*st*sb*m1*mg**2*
     +    s**(-2)*t1 + 32*h_l*st*sb*m1*mg**2*s**(-2)*u2 + 32*h_l*st*sb*
     +    m1*mg**2*s**(-1)*t1*u1**(-1) + 48*h_l*st*sb*m1*mg**2*s**(-1)*
     +    u1**(-1)*u2 - 32*h_l*st*sb*m1*msb2**2*s**(-2)*t1 - 32*h_l*st*
     +    sb*m1*msb2**2*s**(-2)*u2 - 16*h_l*st*sb*m1*msb2**2*s**(-1)*t1
     +    *u1**(-1) - 32*h_l*st*sb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 16*
     +    h_l*st*sb*m1*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_l*st*sb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 - 64*h_l*st*sb*m1*s**(-2)*t1*u2
     +     - 32*h_l*st*sb*m1*s**(-2)*t1**2 - 32*h_l*st*sb*m1*s**(-2)*
     +    u2**2 - 128*h_l*st*sb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*st*
     +    sb*m1*s**(-1)*t1 - 64*h_l*st*sb*m1*s**(-1)*t1**2*u1**(-1) - 
     +    64*h_l*st*sb*m1*s**(-1)*u1**(-1)*u2**2 - 32*h_l*st*sb*m1*
     +    s**(-1)*u2 - 64*h_l*st*sb*m1*t1*u1**(-1) - 48*h_l*st*sb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1**2*mg*s**(-2)*t1 + 32*h_l*ct*sb*m1**2
     +    *mg*s**(-2)*u2 + 32*h_l*ct*sb*m1**2*mg*s**(-1)*t1*u1**(-1) + 
     +    64*h_l*ct*sb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*
     +    m2**2*mg*s**(-2)*t1 - 32*h_l*ct*sb*m2**2*mg*s**(-2)*u2 - 32*
     +    h_l*ct*sb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_l*ct*sb*mg*
     +    s**(-2)*t1*u2 - 32*h_l*ct*sb*mg*s**(-2)*t1**2 - 64*h_l*ct*sb*
     +    mg*s**(-1)*t1*u1**(-1)*u2 - 32*h_l*ct*sb*mg*s**(-1)*t1 - 64*
     +    h_l*ct*sb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_l*ct*sb*mg*s**(-1)
     +    *u2 - 64*h_l*ct*sb*mg*t1*u1**(-1) + 32*h_l*ct*sb*mg*u1**(-1)*
     +    u2 + 32*h_r*st*cb*m1**2*mg*s**(-2)*t1 + 32*h_r*st*cb*m1**2*mg
     +    *s**(-2)*u2 + 32*h_r*st*cb*m1**2*mg*s**(-1)*t1*u1**(-1) + 64*
     +    h_r*st*cb*m1**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*st*cb*m2**2*
     +    mg*s**(-2)*t1 - 32*h_r*st*cb*m2**2*mg*s**(-2)*u2 - 32*h_r*st*
     +    cb*m2**2*mg*s**(-1)*u1**(-1)*u2 - 32*h_r*st*cb*mg*s**(-2)*t1*
     +    u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_r*st*cb*mg*s**(-2)*t1**2 - 64*h_r*st*cb*mg*
     +    s**(-1)*t1*u1**(-1)*u2 - 32*h_r*st*cb*mg*s**(-1)*t1 - 64*h_r*
     +    st*cb*mg*s**(-1)*t1**2*u1**(-1) + 32*h_r*st*cb*mg*s**(-1)*u2
     +     - 64*h_r*st*cb*mg*t1*u1**(-1) + 32*h_r*st*cb*mg*u1**(-1)*u2
     +     + 16*h_r*ct*cb*m1*m2**2*s**(-1)*t1*u1**(-1) + 16*h_r*ct*cb*
     +    m1*m2**2*s**(-1)*u1**(-1)*u2 + 32*h_r*ct*cb*m1*mg**2*s**(-2)*
     +    t1 + 32*h_r*ct*cb*m1*mg**2*s**(-2)*u2 + 32*h_r*ct*cb*m1*mg**2
     +    *s**(-1)*t1*u1**(-1) + 48*h_r*ct*cb*m1*mg**2*s**(-1)*u1**(-1)
     +    *u2 - 32*h_r*ct*cb*m1*msb2**2*s**(-2)*t1 - 32*h_r*ct*cb*m1*
     +    msb2**2*s**(-2)*u2 - 16*h_r*ct*cb*m1*msb2**2*s**(-1)*t1*
     +    u1**(-1) - 32*h_r*ct*cb*m1*msb2**2*s**(-1)*u1**(-1)*u2 - 16*
     +    h_r*ct*cb*m1*mst2**2*s**(-1)*t1*u1**(-1) - 16*h_r*ct*cb*m1*
     +    mst2**2*s**(-1)*u1**(-1)*u2 - 64*h_r*ct*cb*m1*s**(-2)*t1*u2
     +     - 32*h_r*ct*cb*m1*s**(-2)*t1**2 - 32*h_r*ct*cb*m1*s**(-2)*
     +    u2**2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_bth*Cf*Pi**2*alphas**2*
     + prefac * (  - 128*h_r*ct*cb*m1*s**(-1)*t1*u1**(-1)*u2 - 32*h_r*
     +    ct*cb*m1*s**(-1)*t1 - 64*h_r*ct*cb*m1*s**(-1)*t1**2*u1**(-1)
     +     - 64*h_r*ct*cb*m1*s**(-1)*u1**(-1)*u2**2 - 32*h_r*ct*cb*m1*
     +    s**(-1)*u2 - 64*h_r*ct*cb*m1*t1*u1**(-1) - 48*h_r*ct*cb*m1*
     +    u1**(-1)*u2 )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_box*kaellen**(-1)*Nc**2*Cf*
     + Pi**2*alphas**2*prefac * (  - 64*h_l*st*sb*m1*mg**2*s**(-1)*t1*
     +    u1 - 32*h_l*st*sb*m1*mg**2*s**(-1)*t1**2 - 32*h_l*st*sb*m1*
     +    mg**2*s**(-1)*u1**2 + 64*h_l*st*sb*m1*msb2**2*s**(-1)*t1*u1
     +     + 32*h_l*st*sb*m1*msb2**2*s**(-1)*t1**2 + 32*h_l*st*sb*m1*
     +    msb2**2*s**(-1)*u1**2 + 64*h_l*st*sb*m1*t1*u1 + 32*h_l*st*sb*
     +    m1*t1**2 + 32*h_l*st*sb*m1*u1**2 + 128*h_l*st*sb*m1**3*mg**2
     +     - 128*h_l*st*sb*m1**3*msb2**2 - 128*h_l*st*sb*m1**3*s - 64*
     +    h_r*ct*cb*m1*mg**2*s**(-1)*t1*u1 - 32*h_r*ct*cb*m1*mg**2*
     +    s**(-1)*t1**2 - 32*h_r*ct*cb*m1*mg**2*s**(-1)*u1**2 + 64*h_r*
     +    ct*cb*m1*msb2**2*s**(-1)*t1*u1 + 32*h_r*ct*cb*m1*msb2**2*
     +    s**(-1)*t1**2 + 32*h_r*ct*cb*m1*msb2**2*s**(-1)*u1**2 + 64*
     +    h_r*ct*cb*m1*t1*u1 + 32*h_r*ct*cb*m1*t1**2 + 32*h_r*ct*cb*m1*
     +    u1**2 + 128*h_r*ct*cb*m1**3*mg**2 - 128*h_r*ct*cb*m1**3*
     +    msb2**2 - 128*h_r*ct*cb*m1**3*s )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*st*sb*m1*mg**2*u1**(-1) + 32*h_l*st*sb*m1*
     +    mst2**2*u1**(-1) - 32*h_l*st*sb*m1*s*u1**(-1) - 16*h_l*st*sb*
     +    m1*t1*u1**(-1) - 80*h_l*st*sb*m1 - 32*h_l*st*sb*m1**3*
     +    u1**(-1) + 16*h_l*ct*sb*mg*s**(-1)*t1 + 16*h_l*ct*sb*mg*
     +    s**(-1)*u1 + 16*h_l*ct*sb*mg*t1*u1**(-1) + 16*h_l*ct*sb*mg + 
     +    16*h_r*st*cb*mg*s**(-1)*t1 + 16*h_r*st*cb*mg*s**(-1)*u1 + 16*
     +    h_r*st*cb*mg*t1*u1**(-1) + 16*h_r*st*cb*mg - 32*h_r*ct*cb*m1*
     +    mg**2*u1**(-1) + 32*h_r*ct*cb*m1*mst2**2*u1**(-1) - 32*h_r*ct
     +    *cb*m1*s*u1**(-1) - 16*h_r*ct*cb*m1*t1*u1**(-1) - 80*h_r*ct*
     +    cb*m1 - 32*h_r*ct*cb*m1**3*u1**(-1) )
      MM_s = MM_s + SCC(6,5)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*st*sb*m1*mg**2*s**(-1) + 32*h_l*st*sb*m1*mg**2
     +    *u1**(-1) - 32*h_l*st*sb*m1*msb2**2*s**(-1) - 32*h_l*st*sb*m1
     +    *mst2**2*u1**(-1) + 32*h_l*st*sb*m1*s*u1**(-1) + 32*h_l*st*sb
     +    *m1*t1*u1**(-1) + 32*h_l*st*sb*m1 + 32*h_l*st*sb*m1**3*
     +    u1**(-1) - 16*h_l*ct*sb*mg*s**(-1)*t1 - 16*h_l*ct*sb*mg*
     +    s**(-1)*u1 - 16*h_l*ct*sb*mg*t1*u1**(-1) - 16*h_l*ct*sb*mg - 
     +    16*h_r*st*cb*mg*s**(-1)*t1 - 16*h_r*st*cb*mg*s**(-1)*u1 - 16*
     +    h_r*st*cb*mg*t1*u1**(-1) - 16*h_r*st*cb*mg + 32*h_r*ct*cb*m1*
     +    mg**2*s**(-1) + 32*h_r*ct*cb*m1*mg**2*u1**(-1) - 32*h_r*ct*cb
     +    *m1*msb2**2*s**(-1) - 32*h_r*ct*cb*m1*mst2**2*u1**(-1) + 32*
     +    h_r*ct*cb*m1*s*u1**(-1) + 32*h_r*ct*cb*m1*t1*u1**(-1) + 32*
     +    h_r*ct*cb*m1 + 32*h_r*ct*cb*m1**3*u1**(-1) )
      MM_s = MM_s + SCD(2,1)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 48*h_l*st*cb*m1**2*mg*s*u1**(-1) - 32*h_l*st*cb*
     +    m1**2*mg*t1*u1**(-1) - 16*h_l*st*cb*m1**2*mg + 32*h_l*st*cb*
     +    mg*msb1**2*s**(-1)*t1 + 16*h_l*st*cb*mg*msb1**2*s**(-1)*u1 + 
     +    16*h_l*st*cb*mg*msb1**2 + 16*h_l*st*cb*mg*mst1**2*s*u1**(-1)
     +     + 32*h_l*st*cb*mg*mst1**2*t1*u1**(-1) + 16*h_l*st*cb*mg*
     +    mst1**2 + 16*h_l*st*cb*mg*s + 16*h_l*st*cb*mg*u1 - 32*h_l*st*
     +    cb*mg**3*s**(-1)*t1 - 16*h_l*st*cb*mg**3*s**(-1)*u1 - 16*h_l*
     +    st*cb*mg**3*s*u1**(-1) - 32*h_l*st*cb*mg**3*t1*u1**(-1) - 32*
     +    h_l*st*cb*mg**3 + 64*h_l*ct*cb*m1*mg**2*msb1**2*s**(-1) + 32*
     +    h_l*ct*cb*m1*mg**2*msb1**2*u1**(-1) + 32*h_l*ct*cb*m1*mg**2*
     +    mst1**2*u1**(-1) + 16*h_l*ct*cb*m1*mg**2*s*u1**(-1) - 16*h_l*
     +    ct*cb*m1*mg**2 - 32*h_l*ct*cb*m1*mg**4*s**(-1) - 32*h_l*ct*cb
     +    *m1*mg**4*u1**(-1) - 32*h_l*ct*cb*m1*msb1**2*mst1**2*u1**(-1)
     +     + 48*h_l*ct*cb*m1*msb1**2 - 32*h_l*ct*cb*m1*msb1**4*s**(-1)
     +     + 16*h_l*ct*cb*m1*mst1**2*s*u1**(-1) )
      MM_s = MM_s + SCD(2,1)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 16*h_l*ct*cb*m1*s - 32*h_l*ct*cb*m1**3*mg**2*
     +    u1**(-1) + 32*h_l*ct*cb*m1**3*msb1**2*u1**(-1) - 16*h_l*ct*cb
     +    *m1**3*s*u1**(-1) + 64*h_r*st*sb*m1*mg**2*msb1**2*s**(-1) + 
     +    32*h_r*st*sb*m1*mg**2*msb1**2*u1**(-1) + 32*h_r*st*sb*m1*
     +    mg**2*mst1**2*u1**(-1) + 16*h_r*st*sb*m1*mg**2*s*u1**(-1) - 
     +    16*h_r*st*sb*m1*mg**2 - 32*h_r*st*sb*m1*mg**4*s**(-1) - 32*
     +    h_r*st*sb*m1*mg**4*u1**(-1) - 32*h_r*st*sb*m1*msb1**2*mst1**2
     +    *u1**(-1) + 48*h_r*st*sb*m1*msb1**2 - 32*h_r*st*sb*m1*msb1**4
     +    *s**(-1) + 16*h_r*st*sb*m1*mst1**2*s*u1**(-1) - 16*h_r*st*sb*
     +    m1*s - 32*h_r*st*sb*m1**3*mg**2*u1**(-1) + 32*h_r*st*sb*m1**3
     +    *msb1**2*u1**(-1) - 16*h_r*st*sb*m1**3*s*u1**(-1) - 48*h_r*ct
     +    *sb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*sb*m1**2*mg*t1*u1**(-1)
     +     - 16*h_r*ct*sb*m1**2*mg + 32*h_r*ct*sb*mg*msb1**2*s**(-1)*t1
     +     + 16*h_r*ct*sb*mg*msb1**2*s**(-1)*u1 + 16*h_r*ct*sb*mg*
     +    msb1**2 )
      MM_s = MM_s + SCD(2,1)*h_s(1,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_r*ct*sb*mg*mst1**2*s*u1**(-1) + 32*h_r*ct*sb*mg
     +    *mst1**2*t1*u1**(-1) + 16*h_r*ct*sb*mg*mst1**2 + 16*h_r*ct*sb
     +    *mg*s + 16*h_r*ct*sb*mg*u1 - 32*h_r*ct*sb*mg**3*s**(-1)*t1 - 
     +    16*h_r*ct*sb*mg**3*s**(-1)*u1 - 16*h_r*ct*sb*mg**3*s*u1**(-1)
     +     - 32*h_r*ct*sb*mg**3*t1*u1**(-1) - 32*h_r*ct*sb*mg**3 )
      MM_s = MM_s + SCD(2,2)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 48*h_l*st*sb*m1**2*mg*s*u1**(-1) + 32*h_l*st*sb*
     +    m1**2*mg*t1*u1**(-1) + 16*h_l*st*sb*m1**2*mg - 32*h_l*st*sb*
     +    mg*msb2**2*s**(-1)*t1 - 16*h_l*st*sb*mg*msb2**2*s**(-1)*u1 - 
     +    16*h_l*st*sb*mg*msb2**2 - 16*h_l*st*sb*mg*mst1**2*s*u1**(-1)
     +     - 32*h_l*st*sb*mg*mst1**2*t1*u1**(-1) - 16*h_l*st*sb*mg*
     +    mst1**2 - 16*h_l*st*sb*mg*s - 16*h_l*st*sb*mg*u1 + 32*h_l*st*
     +    sb*mg**3*s**(-1)*t1 + 16*h_l*st*sb*mg**3*s**(-1)*u1 + 16*h_l*
     +    st*sb*mg**3*s*u1**(-1) + 32*h_l*st*sb*mg**3*t1*u1**(-1) + 32*
     +    h_l*st*sb*mg**3 - 64*h_l*ct*sb*m1*mg**2*msb2**2*s**(-1) - 32*
     +    h_l*ct*sb*m1*mg**2*msb2**2*u1**(-1) - 32*h_l*ct*sb*m1*mg**2*
     +    mst1**2*u1**(-1) - 16*h_l*ct*sb*m1*mg**2*s*u1**(-1) + 16*h_l*
     +    ct*sb*m1*mg**2 + 32*h_l*ct*sb*m1*mg**4*s**(-1) + 32*h_l*ct*sb
     +    *m1*mg**4*u1**(-1) + 32*h_l*ct*sb*m1*msb2**2*mst1**2*u1**(-1)
     +     - 48*h_l*ct*sb*m1*msb2**2 + 32*h_l*ct*sb*m1*msb2**4*s**(-1)
     +     - 16*h_l*ct*sb*m1*mst1**2*s*u1**(-1) )
      MM_s = MM_s + SCD(2,2)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_l*ct*sb*m1*s + 32*h_l*ct*sb*m1**3*mg**2*
     +    u1**(-1) - 32*h_l*ct*sb*m1**3*msb2**2*u1**(-1) + 16*h_l*ct*sb
     +    *m1**3*s*u1**(-1) + 64*h_r*st*cb*m1*mg**2*msb2**2*s**(-1) + 
     +    32*h_r*st*cb*m1*mg**2*msb2**2*u1**(-1) + 32*h_r*st*cb*m1*
     +    mg**2*mst1**2*u1**(-1) + 16*h_r*st*cb*m1*mg**2*s*u1**(-1) - 
     +    16*h_r*st*cb*m1*mg**2 - 32*h_r*st*cb*m1*mg**4*s**(-1) - 32*
     +    h_r*st*cb*m1*mg**4*u1**(-1) - 32*h_r*st*cb*m1*msb2**2*mst1**2
     +    *u1**(-1) + 48*h_r*st*cb*m1*msb2**2 - 32*h_r*st*cb*m1*msb2**4
     +    *s**(-1) + 16*h_r*st*cb*m1*mst1**2*s*u1**(-1) - 16*h_r*st*cb*
     +    m1*s - 32*h_r*st*cb*m1**3*mg**2*u1**(-1) + 32*h_r*st*cb*m1**3
     +    *msb2**2*u1**(-1) - 16*h_r*st*cb*m1**3*s*u1**(-1) - 48*h_r*ct
     +    *cb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*cb*m1**2*mg*t1*u1**(-1)
     +     - 16*h_r*ct*cb*m1**2*mg + 32*h_r*ct*cb*mg*msb2**2*s**(-1)*t1
     +     + 16*h_r*ct*cb*mg*msb2**2*s**(-1)*u1 + 16*h_r*ct*cb*mg*
     +    msb2**2 )
      MM_s = MM_s + SCD(2,2)*h_s(1,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 16*h_r*ct*cb*mg*mst1**2*s*u1**(-1) + 32*h_r*ct*cb*mg
     +    *mst1**2*t1*u1**(-1) + 16*h_r*ct*cb*mg*mst1**2 + 16*h_r*ct*cb
     +    *mg*s + 16*h_r*ct*cb*mg*u1 - 32*h_r*ct*cb*mg**3*s**(-1)*t1 - 
     +    16*h_r*ct*cb*mg**3*s**(-1)*u1 - 16*h_r*ct*cb*mg**3*s*u1**(-1)
     +     - 32*h_r*ct*cb*mg**3*t1*u1**(-1) - 32*h_r*ct*cb*mg**3 )
      MM_s = MM_s + SCD(2,3)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 64*h_l*st*cb*m1*mg**2*msb1**2*s**(-1) - 32*h_l*st
     +    *cb*m1*mg**2*msb1**2*u1**(-1) - 32*h_l*st*cb*m1*mg**2*mst2**2
     +    *u1**(-1) - 16*h_l*st*cb*m1*mg**2*s*u1**(-1) + 16*h_l*st*cb*
     +    m1*mg**2 + 32*h_l*st*cb*m1*mg**4*s**(-1) + 32*h_l*st*cb*m1*
     +    mg**4*u1**(-1) + 32*h_l*st*cb*m1*msb1**2*mst2**2*u1**(-1) - 
     +    48*h_l*st*cb*m1*msb1**2 + 32*h_l*st*cb*m1*msb1**4*s**(-1) - 
     +    16*h_l*st*cb*m1*mst2**2*s*u1**(-1) + 16*h_l*st*cb*m1*s + 32*
     +    h_l*st*cb*m1**3*mg**2*u1**(-1) - 32*h_l*st*cb*m1**3*msb1**2*
     +    u1**(-1) + 16*h_l*st*cb*m1**3*s*u1**(-1) - 48*h_l*ct*cb*m1**2
     +    *mg*s*u1**(-1) - 32*h_l*ct*cb*m1**2*mg*t1*u1**(-1) - 16*h_l*
     +    ct*cb*m1**2*mg + 32*h_l*ct*cb*mg*msb1**2*s**(-1)*t1 + 16*h_l*
     +    ct*cb*mg*msb1**2*s**(-1)*u1 + 16*h_l*ct*cb*mg*msb1**2 + 16*
     +    h_l*ct*cb*mg*mst2**2*s*u1**(-1) + 32*h_l*ct*cb*mg*mst2**2*t1*
     +    u1**(-1) + 16*h_l*ct*cb*mg*mst2**2 + 16*h_l*ct*cb*mg*s + 16*
     +    h_l*ct*cb*mg*u1 )
      MM_s = MM_s + SCD(2,3)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_l*ct*cb*mg**3*s**(-1)*t1 - 16*h_l*ct*cb*
     +    mg**3*s**(-1)*u1 - 16*h_l*ct*cb*mg**3*s*u1**(-1) - 32*h_l*ct*
     +    cb*mg**3*t1*u1**(-1) - 32*h_l*ct*cb*mg**3 + 48*h_r*st*sb*
     +    m1**2*mg*s*u1**(-1) + 32*h_r*st*sb*m1**2*mg*t1*u1**(-1) + 16*
     +    h_r*st*sb*m1**2*mg - 32*h_r*st*sb*mg*msb1**2*s**(-1)*t1 - 16*
     +    h_r*st*sb*mg*msb1**2*s**(-1)*u1 - 16*h_r*st*sb*mg*msb1**2 - 
     +    16*h_r*st*sb*mg*mst2**2*s*u1**(-1) - 32*h_r*st*sb*mg*mst2**2*
     +    t1*u1**(-1) - 16*h_r*st*sb*mg*mst2**2 - 16*h_r*st*sb*mg*s - 
     +    16*h_r*st*sb*mg*u1 + 32*h_r*st*sb*mg**3*s**(-1)*t1 + 16*h_r*
     +    st*sb*mg**3*s**(-1)*u1 + 16*h_r*st*sb*mg**3*s*u1**(-1) + 32*
     +    h_r*st*sb*mg**3*t1*u1**(-1) + 32*h_r*st*sb*mg**3 + 64*h_r*ct*
     +    sb*m1*mg**2*msb1**2*s**(-1) + 32*h_r*ct*sb*m1*mg**2*msb1**2*
     +    u1**(-1) + 32*h_r*ct*sb*m1*mg**2*mst2**2*u1**(-1) + 16*h_r*ct
     +    *sb*m1*mg**2*s*u1**(-1) - 16*h_r*ct*sb*m1*mg**2 - 32*h_r*ct*
     +    sb*m1*mg**4*s**(-1) )
      MM_s = MM_s + SCD(2,3)*h_s(2,1)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_r*ct*sb*m1*mg**4*u1**(-1) - 32*h_r*ct*sb*m1*
     +    msb1**2*mst2**2*u1**(-1) + 48*h_r*ct*sb*m1*msb1**2 - 32*h_r*
     +    ct*sb*m1*msb1**4*s**(-1) + 16*h_r*ct*sb*m1*mst2**2*s*u1**(-1)
     +     - 16*h_r*ct*sb*m1*s - 32*h_r*ct*sb*m1**3*mg**2*u1**(-1) + 32
     +    *h_r*ct*sb*m1**3*msb1**2*u1**(-1) - 16*h_r*ct*sb*m1**3*s*
     +    u1**(-1) )
      MM_s = MM_s + SCD(2,4)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 64*h_l*st*sb*m1*mg**2*msb2**2*s**(-1) + 32*h_l*st*sb
     +    *m1*mg**2*msb2**2*u1**(-1) + 32*h_l*st*sb*m1*mg**2*mst2**2*
     +    u1**(-1) + 16*h_l*st*sb*m1*mg**2*s*u1**(-1) - 16*h_l*st*sb*m1
     +    *mg**2 - 32*h_l*st*sb*m1*mg**4*s**(-1) - 32*h_l*st*sb*m1*
     +    mg**4*u1**(-1) - 32*h_l*st*sb*m1*msb2**2*mst2**2*u1**(-1) + 
     +    48*h_l*st*sb*m1*msb2**2 - 32*h_l*st*sb*m1*msb2**4*s**(-1) + 
     +    16*h_l*st*sb*m1*mst2**2*s*u1**(-1) - 16*h_l*st*sb*m1*s - 32*
     +    h_l*st*sb*m1**3*mg**2*u1**(-1) + 32*h_l*st*sb*m1**3*msb2**2*
     +    u1**(-1) - 16*h_l*st*sb*m1**3*s*u1**(-1) + 48*h_l*ct*sb*m1**2
     +    *mg*s*u1**(-1) + 32*h_l*ct*sb*m1**2*mg*t1*u1**(-1) + 16*h_l*
     +    ct*sb*m1**2*mg - 32*h_l*ct*sb*mg*msb2**2*s**(-1)*t1 - 16*h_l*
     +    ct*sb*mg*msb2**2*s**(-1)*u1 - 16*h_l*ct*sb*mg*msb2**2 - 16*
     +    h_l*ct*sb*mg*mst2**2*s*u1**(-1) - 32*h_l*ct*sb*mg*mst2**2*t1*
     +    u1**(-1) - 16*h_l*ct*sb*mg*mst2**2 - 16*h_l*ct*sb*mg*s - 16*
     +    h_l*ct*sb*mg*u1 )
      MM_s = MM_s + SCD(2,4)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * ( 32*h_l*ct*sb*mg**3*s**(-1)*t1 + 16*h_l*ct*sb*mg**3*
     +    s**(-1)*u1 + 16*h_l*ct*sb*mg**3*s*u1**(-1) + 32*h_l*ct*sb*
     +    mg**3*t1*u1**(-1) + 32*h_l*ct*sb*mg**3 + 48*h_r*st*cb*m1**2*
     +    mg*s*u1**(-1) + 32*h_r*st*cb*m1**2*mg*t1*u1**(-1) + 16*h_r*st
     +    *cb*m1**2*mg - 32*h_r*st*cb*mg*msb2**2*s**(-1)*t1 - 16*h_r*st
     +    *cb*mg*msb2**2*s**(-1)*u1 - 16*h_r*st*cb*mg*msb2**2 - 16*h_r*
     +    st*cb*mg*mst2**2*s*u1**(-1) - 32*h_r*st*cb*mg*mst2**2*t1*
     +    u1**(-1) - 16*h_r*st*cb*mg*mst2**2 - 16*h_r*st*cb*mg*s - 16*
     +    h_r*st*cb*mg*u1 + 32*h_r*st*cb*mg**3*s**(-1)*t1 + 16*h_r*st*
     +    cb*mg**3*s**(-1)*u1 + 16*h_r*st*cb*mg**3*s*u1**(-1) + 32*h_r*
     +    st*cb*mg**3*t1*u1**(-1) + 32*h_r*st*cb*mg**3 + 64*h_r*ct*cb*
     +    m1*mg**2*msb2**2*s**(-1) + 32*h_r*ct*cb*m1*mg**2*msb2**2*
     +    u1**(-1) + 32*h_r*ct*cb*m1*mg**2*mst2**2*u1**(-1) + 16*h_r*ct
     +    *cb*m1*mg**2*s*u1**(-1) - 16*h_r*ct*cb*m1*mg**2 - 32*h_r*ct*
     +    cb*m1*mg**4*s**(-1) )
      MM_s = MM_s + SCD(2,4)*h_s(2,2)*susy_box*Nc**2*Cf*Pi**2*alphas**2
     + *prefac * (  - 32*h_r*ct*cb*m1*mg**4*u1**(-1) - 32*h_r*ct*cb*m1*
     +    msb2**2*mst2**2*u1**(-1) + 48*h_r*ct*cb*m1*msb2**2 - 32*h_r*
     +    ct*cb*m1*msb2**4*s**(-1) + 16*h_r*ct*cb*m1*mst2**2*s*u1**(-1)
     +     - 16*h_r*ct*cb*m1*s - 32*h_r*ct*cb*m1**3*mg**2*u1**(-1) + 32
     +    *h_r*ct*cb*m1**3*msb2**2*u1**(-1) - 16*h_r*ct*cb*m1**3*s*
     +    u1**(-1) )
      MM_s = MM_s + SCD(3,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*cb*m1**2*mg*s*u1**(-1) - 32*h_l*st*cb*
     +    m1**2*mg*t1*u1**(-1) - 16*h_l*st*cb*m1**2*mg + 32*h_l*st*cb*
     +    mg*msb1**2*s**(-1)*t1 + 16*h_l*st*cb*mg*msb1**2*s**(-1)*u1 + 
     +    16*h_l*st*cb*mg*msb1**2 + 16*h_l*st*cb*mg*mst1**2*s*u1**(-1)
     +     + 32*h_l*st*cb*mg*mst1**2*t1*u1**(-1) + 16*h_l*st*cb*mg*
     +    mst1**2 - 16*h_l*st*cb*mg*s*t1*u1**(-1) - 16*h_l*st*cb*mg*t1
     +     - 32*h_l*st*cb*mg*t1**2*u1**(-1) - 32*h_l*st*cb*mg**3*
     +    s**(-1)*t1 - 16*h_l*st*cb*mg**3*s**(-1)*u1 - 16*h_l*st*cb*
     +    mg**3*s*u1**(-1) - 32*h_l*st*cb*mg**3*t1*u1**(-1) - 32*h_l*st
     +    *cb*mg**3 + 64*h_l*ct*cb*m1*mg**2*msb1**2*s**(-1) + 32*h_l*ct
     +    *cb*m1*mg**2*msb1**2*u1**(-1) + 32*h_l*ct*cb*m1*mg**2*mst1**2
     +    *u1**(-1) - 32*h_l*ct*cb*m1*mg**2*t1*u1**(-1) - 32*h_l*ct*cb*
     +    m1*mg**4*s**(-1) - 32*h_l*ct*cb*m1*mg**4*u1**(-1) - 32*h_l*ct
     +    *cb*m1*msb1**2*mst1**2*u1**(-1) + 32*h_l*ct*cb*m1*msb1**2*s*
     +    u1**(-1) )
      MM_s = MM_s + SCD(3,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*cb*m1*msb1**2*t1*u1**(-1) + 32*h_l*ct*cb*m1
     +    *msb1**2 - 32*h_l*ct*cb*m1*msb1**4*s**(-1) - 32*h_l*ct*cb*
     +    m1**3*mg**2*u1**(-1) + 32*h_l*ct*cb*m1**3*msb1**2*u1**(-1) + 
     +    64*h_r*st*sb*m1*mg**2*msb1**2*s**(-1) + 32*h_r*st*sb*m1*mg**2
     +    *msb1**2*u1**(-1) + 32*h_r*st*sb*m1*mg**2*mst1**2*u1**(-1) - 
     +    32*h_r*st*sb*m1*mg**2*t1*u1**(-1) - 32*h_r*st*sb*m1*mg**4*
     +    s**(-1) - 32*h_r*st*sb*m1*mg**4*u1**(-1) - 32*h_r*st*sb*m1*
     +    msb1**2*mst1**2*u1**(-1) + 32*h_r*st*sb*m1*msb1**2*s*u1**(-1)
     +     + 32*h_r*st*sb*m1*msb1**2*t1*u1**(-1) + 32*h_r*st*sb*m1*
     +    msb1**2 - 32*h_r*st*sb*m1*msb1**4*s**(-1) - 32*h_r*st*sb*
     +    m1**3*mg**2*u1**(-1) + 32*h_r*st*sb*m1**3*msb1**2*u1**(-1) - 
     +    16*h_r*ct*sb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*sb*m1**2*mg*t1*
     +    u1**(-1) - 16*h_r*ct*sb*m1**2*mg + 32*h_r*ct*sb*mg*msb1**2*
     +    s**(-1)*t1 + 16*h_r*ct*sb*mg*msb1**2*s**(-1)*u1 + 16*h_r*ct*
     +    sb*mg*msb1**2 )
      MM_s = MM_s + SCD(3,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_r*ct*sb*mg*mst1**2*s*u1**(-1) + 32*h_r*ct*sb*mg*
     +    mst1**2*t1*u1**(-1) + 16*h_r*ct*sb*mg*mst1**2 - 16*h_r*ct*sb*
     +    mg*s*t1*u1**(-1) - 16*h_r*ct*sb*mg*t1 - 32*h_r*ct*sb*mg*t1**2
     +    *u1**(-1) - 32*h_r*ct*sb*mg**3*s**(-1)*t1 - 16*h_r*ct*sb*
     +    mg**3*s**(-1)*u1 - 16*h_r*ct*sb*mg**3*s*u1**(-1) - 32*h_r*ct*
     +    sb*mg**3*t1*u1**(-1) - 32*h_r*ct*sb*mg**3 )
      MM_s = MM_s + SCD(3,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*sb*m1**2*mg*s*u1**(-1) + 32*h_l*st*sb*m1**2
     +    *mg*t1*u1**(-1) + 16*h_l*st*sb*m1**2*mg - 32*h_l*st*sb*mg*
     +    msb2**2*s**(-1)*t1 - 16*h_l*st*sb*mg*msb2**2*s**(-1)*u1 - 16*
     +    h_l*st*sb*mg*msb2**2 - 16*h_l*st*sb*mg*mst1**2*s*u1**(-1) - 
     +    32*h_l*st*sb*mg*mst1**2*t1*u1**(-1) - 16*h_l*st*sb*mg*mst1**2
     +     + 16*h_l*st*sb*mg*s*t1*u1**(-1) + 16*h_l*st*sb*mg*t1 + 32*
     +    h_l*st*sb*mg*t1**2*u1**(-1) + 32*h_l*st*sb*mg**3*s**(-1)*t1
     +     + 16*h_l*st*sb*mg**3*s**(-1)*u1 + 16*h_l*st*sb*mg**3*s*
     +    u1**(-1) + 32*h_l*st*sb*mg**3*t1*u1**(-1) + 32*h_l*st*sb*
     +    mg**3 - 64*h_l*ct*sb*m1*mg**2*msb2**2*s**(-1) - 32*h_l*ct*sb*
     +    m1*mg**2*msb2**2*u1**(-1) - 32*h_l*ct*sb*m1*mg**2*mst1**2*
     +    u1**(-1) + 32*h_l*ct*sb*m1*mg**2*t1*u1**(-1) + 32*h_l*ct*sb*
     +    m1*mg**4*s**(-1) + 32*h_l*ct*sb*m1*mg**4*u1**(-1) + 32*h_l*ct
     +    *sb*m1*msb2**2*mst1**2*u1**(-1) - 32*h_l*ct*sb*m1*msb2**2*s*
     +    u1**(-1) )
      MM_s = MM_s + SCD(3,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*sb*m1*msb2**2*t1*u1**(-1) - 32*h_l*ct*sb
     +    *m1*msb2**2 + 32*h_l*ct*sb*m1*msb2**4*s**(-1) + 32*h_l*ct*sb*
     +    m1**3*mg**2*u1**(-1) - 32*h_l*ct*sb*m1**3*msb2**2*u1**(-1) + 
     +    64*h_r*st*cb*m1*mg**2*msb2**2*s**(-1) + 32*h_r*st*cb*m1*mg**2
     +    *msb2**2*u1**(-1) + 32*h_r*st*cb*m1*mg**2*mst1**2*u1**(-1) - 
     +    32*h_r*st*cb*m1*mg**2*t1*u1**(-1) - 32*h_r*st*cb*m1*mg**4*
     +    s**(-1) - 32*h_r*st*cb*m1*mg**4*u1**(-1) - 32*h_r*st*cb*m1*
     +    msb2**2*mst1**2*u1**(-1) + 32*h_r*st*cb*m1*msb2**2*s*u1**(-1)
     +     + 32*h_r*st*cb*m1*msb2**2*t1*u1**(-1) + 32*h_r*st*cb*m1*
     +    msb2**2 - 32*h_r*st*cb*m1*msb2**4*s**(-1) - 32*h_r*st*cb*
     +    m1**3*mg**2*u1**(-1) + 32*h_r*st*cb*m1**3*msb2**2*u1**(-1) - 
     +    16*h_r*ct*cb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*cb*m1**2*mg*t1*
     +    u1**(-1) - 16*h_r*ct*cb*m1**2*mg + 32*h_r*ct*cb*mg*msb2**2*
     +    s**(-1)*t1 + 16*h_r*ct*cb*mg*msb2**2*s**(-1)*u1 + 16*h_r*ct*
     +    cb*mg*msb2**2 )
      MM_s = MM_s + SCD(3,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_r*ct*cb*mg*mst1**2*s*u1**(-1) + 32*h_r*ct*cb*mg*
     +    mst1**2*t1*u1**(-1) + 16*h_r*ct*cb*mg*mst1**2 - 16*h_r*ct*cb*
     +    mg*s*t1*u1**(-1) - 16*h_r*ct*cb*mg*t1 - 32*h_r*ct*cb*mg*t1**2
     +    *u1**(-1) - 32*h_r*ct*cb*mg**3*s**(-1)*t1 - 16*h_r*ct*cb*
     +    mg**3*s**(-1)*u1 - 16*h_r*ct*cb*mg**3*s*u1**(-1) - 32*h_r*ct*
     +    cb*mg**3*t1*u1**(-1) - 32*h_r*ct*cb*mg**3 )
      MM_s = MM_s + SCD(3,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 64*h_l*st*cb*m1*mg**2*msb1**2*s**(-1) - 32*h_l*st*
     +    cb*m1*mg**2*msb1**2*u1**(-1) - 32*h_l*st*cb*m1*mg**2*mst2**2*
     +    u1**(-1) + 32*h_l*st*cb*m1*mg**2*t1*u1**(-1) + 32*h_l*st*cb*
     +    m1*mg**4*s**(-1) + 32*h_l*st*cb*m1*mg**4*u1**(-1) + 32*h_l*st
     +    *cb*m1*msb1**2*mst2**2*u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s*
     +    u1**(-1) - 32*h_l*st*cb*m1*msb1**2*t1*u1**(-1) - 32*h_l*st*cb
     +    *m1*msb1**2 + 32*h_l*st*cb*m1*msb1**4*s**(-1) + 32*h_l*st*cb*
     +    m1**3*mg**2*u1**(-1) - 32*h_l*st*cb*m1**3*msb1**2*u1**(-1) - 
     +    16*h_l*ct*cb*m1**2*mg*s*u1**(-1) - 32*h_l*ct*cb*m1**2*mg*t1*
     +    u1**(-1) - 16*h_l*ct*cb*m1**2*mg + 32*h_l*ct*cb*mg*msb1**2*
     +    s**(-1)*t1 + 16*h_l*ct*cb*mg*msb1**2*s**(-1)*u1 + 16*h_l*ct*
     +    cb*mg*msb1**2 + 16*h_l*ct*cb*mg*mst2**2*s*u1**(-1) + 32*h_l*
     +    ct*cb*mg*mst2**2*t1*u1**(-1) + 16*h_l*ct*cb*mg*mst2**2 - 16*
     +    h_l*ct*cb*mg*s*t1*u1**(-1) - 16*h_l*ct*cb*mg*t1 - 32*h_l*ct*
     +    cb*mg*t1**2*u1**(-1) )
      MM_s = MM_s + SCD(3,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*mg**3*s**(-1)*t1 - 16*h_l*ct*cb*mg**3
     +    *s**(-1)*u1 - 16*h_l*ct*cb*mg**3*s*u1**(-1) - 32*h_l*ct*cb*
     +    mg**3*t1*u1**(-1) - 32*h_l*ct*cb*mg**3 + 16*h_r*st*sb*m1**2*
     +    mg*s*u1**(-1) + 32*h_r*st*sb*m1**2*mg*t1*u1**(-1) + 16*h_r*st
     +    *sb*m1**2*mg - 32*h_r*st*sb*mg*msb1**2*s**(-1)*t1 - 16*h_r*st
     +    *sb*mg*msb1**2*s**(-1)*u1 - 16*h_r*st*sb*mg*msb1**2 - 16*h_r*
     +    st*sb*mg*mst2**2*s*u1**(-1) - 32*h_r*st*sb*mg*mst2**2*t1*
     +    u1**(-1) - 16*h_r*st*sb*mg*mst2**2 + 16*h_r*st*sb*mg*s*t1*
     +    u1**(-1) + 16*h_r*st*sb*mg*t1 + 32*h_r*st*sb*mg*t1**2*
     +    u1**(-1) + 32*h_r*st*sb*mg**3*s**(-1)*t1 + 16*h_r*st*sb*mg**3
     +    *s**(-1)*u1 + 16*h_r*st*sb*mg**3*s*u1**(-1) + 32*h_r*st*sb*
     +    mg**3*t1*u1**(-1) + 32*h_r*st*sb*mg**3 + 64*h_r*ct*sb*m1*
     +    mg**2*msb1**2*s**(-1) + 32*h_r*ct*sb*m1*mg**2*msb1**2*
     +    u1**(-1) + 32*h_r*ct*sb*m1*mg**2*mst2**2*u1**(-1) - 32*h_r*ct
     +    *sb*m1*mg**2*t1*u1**(-1) )
      MM_s = MM_s + SCD(3,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_r*ct*sb*m1*mg**4*s**(-1) - 32*h_r*ct*sb*m1*
     +    mg**4*u1**(-1) - 32*h_r*ct*sb*m1*msb1**2*mst2**2*u1**(-1) + 
     +    32*h_r*ct*sb*m1*msb1**2*s*u1**(-1) + 32*h_r*ct*sb*m1*msb1**2*
     +    t1*u1**(-1) + 32*h_r*ct*sb*m1*msb1**2 - 32*h_r*ct*sb*m1*
     +    msb1**4*s**(-1) - 32*h_r*ct*sb*m1**3*mg**2*u1**(-1) + 32*h_r*
     +    ct*sb*m1**3*msb1**2*u1**(-1) )
      MM_s = MM_s + SCD(3,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 64*h_l*st*sb*m1*mg**2*msb2**2*s**(-1) + 32*h_l*st*sb*
     +    m1*mg**2*msb2**2*u1**(-1) + 32*h_l*st*sb*m1*mg**2*mst2**2*
     +    u1**(-1) - 32*h_l*st*sb*m1*mg**2*t1*u1**(-1) - 32*h_l*st*sb*
     +    m1*mg**4*s**(-1) - 32*h_l*st*sb*m1*mg**4*u1**(-1) - 32*h_l*st
     +    *sb*m1*msb2**2*mst2**2*u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s*
     +    u1**(-1) + 32*h_l*st*sb*m1*msb2**2*t1*u1**(-1) + 32*h_l*st*sb
     +    *m1*msb2**2 - 32*h_l*st*sb*m1*msb2**4*s**(-1) - 32*h_l*st*sb*
     +    m1**3*mg**2*u1**(-1) + 32*h_l*st*sb*m1**3*msb2**2*u1**(-1) + 
     +    16*h_l*ct*sb*m1**2*mg*s*u1**(-1) + 32*h_l*ct*sb*m1**2*mg*t1*
     +    u1**(-1) + 16*h_l*ct*sb*m1**2*mg - 32*h_l*ct*sb*mg*msb2**2*
     +    s**(-1)*t1 - 16*h_l*ct*sb*mg*msb2**2*s**(-1)*u1 - 16*h_l*ct*
     +    sb*mg*msb2**2 - 16*h_l*ct*sb*mg*mst2**2*s*u1**(-1) - 32*h_l*
     +    ct*sb*mg*mst2**2*t1*u1**(-1) - 16*h_l*ct*sb*mg*mst2**2 + 16*
     +    h_l*ct*sb*mg*s*t1*u1**(-1) + 16*h_l*ct*sb*mg*t1 + 32*h_l*ct*
     +    sb*mg*t1**2*u1**(-1) )
      MM_s = MM_s + SCD(3,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*mg**3*s**(-1)*t1 + 16*h_l*ct*sb*mg**3*
     +    s**(-1)*u1 + 16*h_l*ct*sb*mg**3*s*u1**(-1) + 32*h_l*ct*sb*
     +    mg**3*t1*u1**(-1) + 32*h_l*ct*sb*mg**3 + 16*h_r*st*cb*m1**2*
     +    mg*s*u1**(-1) + 32*h_r*st*cb*m1**2*mg*t1*u1**(-1) + 16*h_r*st
     +    *cb*m1**2*mg - 32*h_r*st*cb*mg*msb2**2*s**(-1)*t1 - 16*h_r*st
     +    *cb*mg*msb2**2*s**(-1)*u1 - 16*h_r*st*cb*mg*msb2**2 - 16*h_r*
     +    st*cb*mg*mst2**2*s*u1**(-1) - 32*h_r*st*cb*mg*mst2**2*t1*
     +    u1**(-1) - 16*h_r*st*cb*mg*mst2**2 + 16*h_r*st*cb*mg*s*t1*
     +    u1**(-1) + 16*h_r*st*cb*mg*t1 + 32*h_r*st*cb*mg*t1**2*
     +    u1**(-1) + 32*h_r*st*cb*mg**3*s**(-1)*t1 + 16*h_r*st*cb*mg**3
     +    *s**(-1)*u1 + 16*h_r*st*cb*mg**3*s*u1**(-1) + 32*h_r*st*cb*
     +    mg**3*t1*u1**(-1) + 32*h_r*st*cb*mg**3 + 64*h_r*ct*cb*m1*
     +    mg**2*msb2**2*s**(-1) + 32*h_r*ct*cb*m1*mg**2*msb2**2*
     +    u1**(-1) + 32*h_r*ct*cb*m1*mg**2*mst2**2*u1**(-1) - 32*h_r*ct
     +    *cb*m1*mg**2*t1*u1**(-1) )
      MM_s = MM_s + SCD(3,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_r*ct*cb*m1*mg**4*s**(-1) - 32*h_r*ct*cb*m1*
     +    mg**4*u1**(-1) - 32*h_r*ct*cb*m1*msb2**2*mst2**2*u1**(-1) + 
     +    32*h_r*ct*cb*m1*msb2**2*s*u1**(-1) + 32*h_r*ct*cb*m1*msb2**2*
     +    t1*u1**(-1) + 32*h_r*ct*cb*m1*msb2**2 - 32*h_r*ct*cb*m1*
     +    msb2**4*s**(-1) - 32*h_r*ct*cb*m1**3*mg**2*u1**(-1) + 32*h_r*
     +    ct*cb*m1**3*msb2**2*u1**(-1) )
      MM_s = MM_s + SCD(4,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*st*cb*m1**2*mg*s*u1**(-1) - 32*h_l*st*cb*
     +    m1**2*mg*t1*u1**(-1) - 16*h_l*st*cb*m1**2*mg - 32*h_l*st*cb*
     +    mg*msb1**2*s**(-1)*t1 - 16*h_l*st*cb*mg*msb1**2*s**(-1)*u1 - 
     +    16*h_l*st*cb*mg*msb1**2 - 16*h_l*st*cb*mg*mst1**2*s*u1**(-1)
     +     - 32*h_l*st*cb*mg*mst1**2*t1*u1**(-1) - 16*h_l*st*cb*mg*
     +    mst1**2 + 16*h_l*st*cb*mg*s**(-1)*t1*u1 + 32*h_l*st*cb*mg*
     +    s**(-1)*t1**2 + 16*h_l*st*cb*mg*t1 + 32*h_l*st*cb*mg**3*
     +    s**(-1)*t1 + 16*h_l*st*cb*mg**3*s**(-1)*u1 + 16*h_l*st*cb*
     +    mg**3*s*u1**(-1) + 32*h_l*st*cb*mg**3*t1*u1**(-1) + 32*h_l*st
     +    *cb*mg**3 + 64*h_l*ct*cb*m1*mg**2*msb1**2*s**(-1) + 32*h_l*ct
     +    *cb*m1*mg**2*msb1**2*u1**(-1) + 32*h_l*ct*cb*m1*mg**2*mst1**2
     +    *u1**(-1) - 32*h_l*ct*cb*m1*mg**2*s**(-1)*t1 + 16*h_l*ct*cb*
     +    m1*mg**2*s**(-1)*u1 + 16*h_l*ct*cb*m1*mg**2*s*u1**(-1) + 32*
     +    h_l*ct*cb*m1*mg**2 - 32*h_l*ct*cb*m1*mg**4*s**(-1) - 32*h_l*
     +    ct*cb*m1*mg**4*u1**(-1) )
      MM_s = MM_s + SCD(4,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 32*h_l*ct*cb*m1*msb1**2*mst1**2*u1**(-1) + 32*h_l*
     +    ct*cb*m1*msb1**2*s**(-1)*t1 - 16*h_l*ct*cb*m1*msb1**2*s**(-1)
     +    *u1 - 16*h_l*ct*cb*m1*msb1**2 - 32*h_l*ct*cb*m1*msb1**4*
     +    s**(-1) + 16*h_l*ct*cb*m1*mst1**2*s*u1**(-1) + 16*h_l*ct*cb*
     +    m1*mst1**2 + 16*h_l*ct*cb*m1*s**(-1)*t1*u1 + 16*h_l*ct*cb*m1*
     +    t1 + 32*h_l*ct*cb*m1**3*mg**2*u1**(-1) - 32*h_l*ct*cb*m1**3*
     +    msb1**2*u1**(-1) - 16*h_l*ct*cb*m1**3*s*u1**(-1) - 16*h_l*ct*
     +    cb*m1**3 + 64*h_r*st*sb*m1*mg**2*msb1**2*s**(-1) + 32*h_r*st*
     +    sb*m1*mg**2*msb1**2*u1**(-1) + 32*h_r*st*sb*m1*mg**2*mst1**2*
     +    u1**(-1) - 32*h_r*st*sb*m1*mg**2*s**(-1)*t1 + 16*h_r*st*sb*m1
     +    *mg**2*s**(-1)*u1 + 16*h_r*st*sb*m1*mg**2*s*u1**(-1) + 32*h_r
     +    *st*sb*m1*mg**2 - 32*h_r*st*sb*m1*mg**4*s**(-1) - 32*h_r*st*
     +    sb*m1*mg**4*u1**(-1) - 32*h_r*st*sb*m1*msb1**2*mst1**2*
     +    u1**(-1) + 32*h_r*st*sb*m1*msb1**2*s**(-1)*t1 - 16*h_r*st*sb*
     +    m1*msb1**2*s**(-1)*u1 )
      MM_s = MM_s + SCD(4,1)*h_s(1,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*sb*m1*msb1**2 - 32*h_r*st*sb*m1*msb1**4*
     +    s**(-1) + 16*h_r*st*sb*m1*mst1**2*s*u1**(-1) + 16*h_r*st*sb*
     +    m1*mst1**2 + 16*h_r*st*sb*m1*s**(-1)*t1*u1 + 16*h_r*st*sb*m1*
     +    t1 + 32*h_r*st*sb*m1**3*mg**2*u1**(-1) - 32*h_r*st*sb*m1**3*
     +    msb1**2*u1**(-1) - 16*h_r*st*sb*m1**3*s*u1**(-1) - 16*h_r*st*
     +    sb*m1**3 - 16*h_r*ct*sb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*sb*
     +    m1**2*mg*t1*u1**(-1) - 16*h_r*ct*sb*m1**2*mg - 32*h_r*ct*sb*
     +    mg*msb1**2*s**(-1)*t1 - 16*h_r*ct*sb*mg*msb1**2*s**(-1)*u1 - 
     +    16*h_r*ct*sb*mg*msb1**2 - 16*h_r*ct*sb*mg*mst1**2*s*u1**(-1)
     +     - 32*h_r*ct*sb*mg*mst1**2*t1*u1**(-1) - 16*h_r*ct*sb*mg*
     +    mst1**2 + 16*h_r*ct*sb*mg*s**(-1)*t1*u1 + 32*h_r*ct*sb*mg*
     +    s**(-1)*t1**2 + 16*h_r*ct*sb*mg*t1 + 32*h_r*ct*sb*mg**3*
     +    s**(-1)*t1 + 16*h_r*ct*sb*mg**3*s**(-1)*u1 + 16*h_r*ct*sb*
     +    mg**3*s*u1**(-1) + 32*h_r*ct*sb*mg**3*t1*u1**(-1) + 32*h_r*ct
     +    *sb*mg**3 )
      MM_s = MM_s + SCD(4,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*st*sb*m1**2*mg*s*u1**(-1) + 32*h_l*st*sb*m1**2
     +    *mg*t1*u1**(-1) + 16*h_l*st*sb*m1**2*mg + 32*h_l*st*sb*mg*
     +    msb2**2*s**(-1)*t1 + 16*h_l*st*sb*mg*msb2**2*s**(-1)*u1 + 16*
     +    h_l*st*sb*mg*msb2**2 + 16*h_l*st*sb*mg*mst1**2*s*u1**(-1) + 
     +    32*h_l*st*sb*mg*mst1**2*t1*u1**(-1) + 16*h_l*st*sb*mg*mst1**2
     +     - 16*h_l*st*sb*mg*s**(-1)*t1*u1 - 32*h_l*st*sb*mg*s**(-1)*
     +    t1**2 - 16*h_l*st*sb*mg*t1 - 32*h_l*st*sb*mg**3*s**(-1)*t1 - 
     +    16*h_l*st*sb*mg**3*s**(-1)*u1 - 16*h_l*st*sb*mg**3*s*u1**(-1)
     +     - 32*h_l*st*sb*mg**3*t1*u1**(-1) - 32*h_l*st*sb*mg**3 - 64*
     +    h_l*ct*sb*m1*mg**2*msb2**2*s**(-1) - 32*h_l*ct*sb*m1*mg**2*
     +    msb2**2*u1**(-1) - 32*h_l*ct*sb*m1*mg**2*mst1**2*u1**(-1) + 
     +    32*h_l*ct*sb*m1*mg**2*s**(-1)*t1 - 16*h_l*ct*sb*m1*mg**2*
     +    s**(-1)*u1 - 16*h_l*ct*sb*m1*mg**2*s*u1**(-1) - 32*h_l*ct*sb*
     +    m1*mg**2 + 32*h_l*ct*sb*m1*mg**4*s**(-1) + 32*h_l*ct*sb*m1*
     +    mg**4*u1**(-1) )
      MM_s = MM_s + SCD(4,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 32*h_l*ct*sb*m1*msb2**2*mst1**2*u1**(-1) - 32*h_l*ct*
     +    sb*m1*msb2**2*s**(-1)*t1 + 16*h_l*ct*sb*m1*msb2**2*s**(-1)*u1
     +     + 16*h_l*ct*sb*m1*msb2**2 + 32*h_l*ct*sb*m1*msb2**4*s**(-1)
     +     - 16*h_l*ct*sb*m1*mst1**2*s*u1**(-1) - 16*h_l*ct*sb*m1*
     +    mst1**2 - 16*h_l*ct*sb*m1*s**(-1)*t1*u1 - 16*h_l*ct*sb*m1*t1
     +     - 32*h_l*ct*sb*m1**3*mg**2*u1**(-1) + 32*h_l*ct*sb*m1**3*
     +    msb2**2*u1**(-1) + 16*h_l*ct*sb*m1**3*s*u1**(-1) + 16*h_l*ct*
     +    sb*m1**3 + 64*h_r*st*cb*m1*mg**2*msb2**2*s**(-1) + 32*h_r*st*
     +    cb*m1*mg**2*msb2**2*u1**(-1) + 32*h_r*st*cb*m1*mg**2*mst1**2*
     +    u1**(-1) - 32*h_r*st*cb*m1*mg**2*s**(-1)*t1 + 16*h_r*st*cb*m1
     +    *mg**2*s**(-1)*u1 + 16*h_r*st*cb*m1*mg**2*s*u1**(-1) + 32*h_r
     +    *st*cb*m1*mg**2 - 32*h_r*st*cb*m1*mg**4*s**(-1) - 32*h_r*st*
     +    cb*m1*mg**4*u1**(-1) - 32*h_r*st*cb*m1*msb2**2*mst1**2*
     +    u1**(-1) + 32*h_r*st*cb*m1*msb2**2*s**(-1)*t1 - 16*h_r*st*cb*
     +    m1*msb2**2*s**(-1)*u1 )
      MM_s = MM_s + SCD(4,2)*h_s(1,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*cb*m1*msb2**2 - 32*h_r*st*cb*m1*msb2**4*
     +    s**(-1) + 16*h_r*st*cb*m1*mst1**2*s*u1**(-1) + 16*h_r*st*cb*
     +    m1*mst1**2 + 16*h_r*st*cb*m1*s**(-1)*t1*u1 + 16*h_r*st*cb*m1*
     +    t1 + 32*h_r*st*cb*m1**3*mg**2*u1**(-1) - 32*h_r*st*cb*m1**3*
     +    msb2**2*u1**(-1) - 16*h_r*st*cb*m1**3*s*u1**(-1) - 16*h_r*st*
     +    cb*m1**3 - 16*h_r*ct*cb*m1**2*mg*s*u1**(-1) - 32*h_r*ct*cb*
     +    m1**2*mg*t1*u1**(-1) - 16*h_r*ct*cb*m1**2*mg - 32*h_r*ct*cb*
     +    mg*msb2**2*s**(-1)*t1 - 16*h_r*ct*cb*mg*msb2**2*s**(-1)*u1 - 
     +    16*h_r*ct*cb*mg*msb2**2 - 16*h_r*ct*cb*mg*mst1**2*s*u1**(-1)
     +     - 32*h_r*ct*cb*mg*mst1**2*t1*u1**(-1) - 16*h_r*ct*cb*mg*
     +    mst1**2 + 16*h_r*ct*cb*mg*s**(-1)*t1*u1 + 32*h_r*ct*cb*mg*
     +    s**(-1)*t1**2 + 16*h_r*ct*cb*mg*t1 + 32*h_r*ct*cb*mg**3*
     +    s**(-1)*t1 + 16*h_r*ct*cb*mg**3*s**(-1)*u1 + 16*h_r*ct*cb*
     +    mg**3*s*u1**(-1) + 32*h_r*ct*cb*mg**3*t1*u1**(-1) + 32*h_r*ct
     +    *cb*mg**3 )
      MM_s = MM_s + SCD(4,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 64*h_l*st*cb*m1*mg**2*msb1**2*s**(-1) - 32*h_l*st*
     +    cb*m1*mg**2*msb1**2*u1**(-1) - 32*h_l*st*cb*m1*mg**2*mst2**2*
     +    u1**(-1) + 32*h_l*st*cb*m1*mg**2*s**(-1)*t1 - 16*h_l*st*cb*m1
     +    *mg**2*s**(-1)*u1 - 16*h_l*st*cb*m1*mg**2*s*u1**(-1) - 32*h_l
     +    *st*cb*m1*mg**2 + 32*h_l*st*cb*m1*mg**4*s**(-1) + 32*h_l*st*
     +    cb*m1*mg**4*u1**(-1) + 32*h_l*st*cb*m1*msb1**2*mst2**2*
     +    u1**(-1) - 32*h_l*st*cb*m1*msb1**2*s**(-1)*t1 + 16*h_l*st*cb*
     +    m1*msb1**2*s**(-1)*u1 + 16*h_l*st*cb*m1*msb1**2 + 32*h_l*st*
     +    cb*m1*msb1**4*s**(-1) - 16*h_l*st*cb*m1*mst2**2*s*u1**(-1) - 
     +    16*h_l*st*cb*m1*mst2**2 - 16*h_l*st*cb*m1*s**(-1)*t1*u1 - 16*
     +    h_l*st*cb*m1*t1 - 32*h_l*st*cb*m1**3*mg**2*u1**(-1) + 32*h_l*
     +    st*cb*m1**3*msb1**2*u1**(-1) + 16*h_l*st*cb*m1**3*s*u1**(-1)
     +     + 16*h_l*st*cb*m1**3 - 16*h_l*ct*cb*m1**2*mg*s*u1**(-1) - 32
     +    *h_l*ct*cb*m1**2*mg*t1*u1**(-1) - 16*h_l*ct*cb*m1**2*mg - 32*
     +    h_l*ct*cb*mg*msb1**2*s**(-1)*t1 )
      MM_s = MM_s + SCD(4,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_l*ct*cb*mg*msb1**2*s**(-1)*u1 - 16*h_l*ct*cb*
     +    mg*msb1**2 - 16*h_l*ct*cb*mg*mst2**2*s*u1**(-1) - 32*h_l*ct*
     +    cb*mg*mst2**2*t1*u1**(-1) - 16*h_l*ct*cb*mg*mst2**2 + 16*h_l*
     +    ct*cb*mg*s**(-1)*t1*u1 + 32*h_l*ct*cb*mg*s**(-1)*t1**2 + 16*
     +    h_l*ct*cb*mg*t1 + 32*h_l*ct*cb*mg**3*s**(-1)*t1 + 16*h_l*ct*
     +    cb*mg**3*s**(-1)*u1 + 16*h_l*ct*cb*mg**3*s*u1**(-1) + 32*h_l*
     +    ct*cb*mg**3*t1*u1**(-1) + 32*h_l*ct*cb*mg**3 + 16*h_r*st*sb*
     +    m1**2*mg*s*u1**(-1) + 32*h_r*st*sb*m1**2*mg*t1*u1**(-1) + 16*
     +    h_r*st*sb*m1**2*mg + 32*h_r*st*sb*mg*msb1**2*s**(-1)*t1 + 16*
     +    h_r*st*sb*mg*msb1**2*s**(-1)*u1 + 16*h_r*st*sb*mg*msb1**2 + 
     +    16*h_r*st*sb*mg*mst2**2*s*u1**(-1) + 32*h_r*st*sb*mg*mst2**2*
     +    t1*u1**(-1) + 16*h_r*st*sb*mg*mst2**2 - 16*h_r*st*sb*mg*
     +    s**(-1)*t1*u1 - 32*h_r*st*sb*mg*s**(-1)*t1**2 - 16*h_r*st*sb*
     +    mg*t1 - 32*h_r*st*sb*mg**3*s**(-1)*t1 - 16*h_r*st*sb*mg**3*
     +    s**(-1)*u1 )
      MM_s = MM_s + SCD(4,3)*h_s(2,1)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*sb*mg**3*s*u1**(-1) - 32*h_r*st*sb*mg**3
     +    *t1*u1**(-1) - 32*h_r*st*sb*mg**3 + 64*h_r*ct*sb*m1*mg**2*
     +    msb1**2*s**(-1) + 32*h_r*ct*sb*m1*mg**2*msb1**2*u1**(-1) + 32
     +    *h_r*ct*sb*m1*mg**2*mst2**2*u1**(-1) - 32*h_r*ct*sb*m1*mg**2*
     +    s**(-1)*t1 + 16*h_r*ct*sb*m1*mg**2*s**(-1)*u1 + 16*h_r*ct*sb*
     +    m1*mg**2*s*u1**(-1) + 32*h_r*ct*sb*m1*mg**2 - 32*h_r*ct*sb*m1
     +    *mg**4*s**(-1) - 32*h_r*ct*sb*m1*mg**4*u1**(-1) - 32*h_r*ct*
     +    sb*m1*msb1**2*mst2**2*u1**(-1) + 32*h_r*ct*sb*m1*msb1**2*
     +    s**(-1)*t1 - 16*h_r*ct*sb*m1*msb1**2*s**(-1)*u1 - 16*h_r*ct*
     +    sb*m1*msb1**2 - 32*h_r*ct*sb*m1*msb1**4*s**(-1) + 16*h_r*ct*
     +    sb*m1*mst2**2*s*u1**(-1) + 16*h_r*ct*sb*m1*mst2**2 + 16*h_r*
     +    ct*sb*m1*s**(-1)*t1*u1 + 16*h_r*ct*sb*m1*t1 + 32*h_r*ct*sb*
     +    m1**3*mg**2*u1**(-1) - 32*h_r*ct*sb*m1**3*msb1**2*u1**(-1) - 
     +    16*h_r*ct*sb*m1**3*s*u1**(-1) - 16*h_r*ct*sb*m1**3 )
      MM_s = MM_s + SCD(4,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 64*h_l*st*sb*m1*mg**2*msb2**2*s**(-1) + 32*h_l*st*sb*
     +    m1*mg**2*msb2**2*u1**(-1) + 32*h_l*st*sb*m1*mg**2*mst2**2*
     +    u1**(-1) - 32*h_l*st*sb*m1*mg**2*s**(-1)*t1 + 16*h_l*st*sb*m1
     +    *mg**2*s**(-1)*u1 + 16*h_l*st*sb*m1*mg**2*s*u1**(-1) + 32*h_l
     +    *st*sb*m1*mg**2 - 32*h_l*st*sb*m1*mg**4*s**(-1) - 32*h_l*st*
     +    sb*m1*mg**4*u1**(-1) - 32*h_l*st*sb*m1*msb2**2*mst2**2*
     +    u1**(-1) + 32*h_l*st*sb*m1*msb2**2*s**(-1)*t1 - 16*h_l*st*sb*
     +    m1*msb2**2*s**(-1)*u1 - 16*h_l*st*sb*m1*msb2**2 - 32*h_l*st*
     +    sb*m1*msb2**4*s**(-1) + 16*h_l*st*sb*m1*mst2**2*s*u1**(-1) + 
     +    16*h_l*st*sb*m1*mst2**2 + 16*h_l*st*sb*m1*s**(-1)*t1*u1 + 16*
     +    h_l*st*sb*m1*t1 + 32*h_l*st*sb*m1**3*mg**2*u1**(-1) - 32*h_l*
     +    st*sb*m1**3*msb2**2*u1**(-1) - 16*h_l*st*sb*m1**3*s*u1**(-1)
     +     - 16*h_l*st*sb*m1**3 + 16*h_l*ct*sb*m1**2*mg*s*u1**(-1) + 32
     +    *h_l*ct*sb*m1**2*mg*t1*u1**(-1) + 16*h_l*ct*sb*m1**2*mg + 32*
     +    h_l*ct*sb*mg*msb2**2*s**(-1)*t1 )
      MM_s = MM_s + SCD(4,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * ( 16*h_l*ct*sb*mg*msb2**2*s**(-1)*u1 + 16*h_l*ct*sb*mg*
     +    msb2**2 + 16*h_l*ct*sb*mg*mst2**2*s*u1**(-1) + 32*h_l*ct*sb*
     +    mg*mst2**2*t1*u1**(-1) + 16*h_l*ct*sb*mg*mst2**2 - 16*h_l*ct*
     +    sb*mg*s**(-1)*t1*u1 - 32*h_l*ct*sb*mg*s**(-1)*t1**2 - 16*h_l*
     +    ct*sb*mg*t1 - 32*h_l*ct*sb*mg**3*s**(-1)*t1 - 16*h_l*ct*sb*
     +    mg**3*s**(-1)*u1 - 16*h_l*ct*sb*mg**3*s*u1**(-1) - 32*h_l*ct*
     +    sb*mg**3*t1*u1**(-1) - 32*h_l*ct*sb*mg**3 + 16*h_r*st*cb*
     +    m1**2*mg*s*u1**(-1) + 32*h_r*st*cb*m1**2*mg*t1*u1**(-1) + 16*
     +    h_r*st*cb*m1**2*mg + 32*h_r*st*cb*mg*msb2**2*s**(-1)*t1 + 16*
     +    h_r*st*cb*mg*msb2**2*s**(-1)*u1 + 16*h_r*st*cb*mg*msb2**2 + 
     +    16*h_r*st*cb*mg*mst2**2*s*u1**(-1) + 32*h_r*st*cb*mg*mst2**2*
     +    t1*u1**(-1) + 16*h_r*st*cb*mg*mst2**2 - 16*h_r*st*cb*mg*
     +    s**(-1)*t1*u1 - 32*h_r*st*cb*mg*s**(-1)*t1**2 - 16*h_r*st*cb*
     +    mg*t1 - 32*h_r*st*cb*mg**3*s**(-1)*t1 - 16*h_r*st*cb*mg**3*
     +    s**(-1)*u1 )
      MM_s = MM_s + SCD(4,4)*h_s(2,2)*susy_box*Cf*Pi**2*alphas**2*
     + prefac * (  - 16*h_r*st*cb*mg**3*s*u1**(-1) - 32*h_r*st*cb*mg**3
     +    *t1*u1**(-1) - 32*h_r*st*cb*mg**3 + 64*h_r*ct*cb*m1*mg**2*
     +    msb2**2*s**(-1) + 32*h_r*ct*cb*m1*mg**2*msb2**2*u1**(-1) + 32
     +    *h_r*ct*cb*m1*mg**2*mst2**2*u1**(-1) - 32*h_r*ct*cb*m1*mg**2*
     +    s**(-1)*t1 + 16*h_r*ct*cb*m1*mg**2*s**(-1)*u1 + 16*h_r*ct*cb*
     +    m1*mg**2*s*u1**(-1) + 32*h_r*ct*cb*m1*mg**2 - 32*h_r*ct*cb*m1
     +    *mg**4*s**(-1) - 32*h_r*ct*cb*m1*mg**4*u1**(-1) - 32*h_r*ct*
     +    cb*m1*msb2**2*mst2**2*u1**(-1) + 32*h_r*ct*cb*m1*msb2**2*
     +    s**(-1)*t1 - 16*h_r*ct*cb*m1*msb2**2*s**(-1)*u1 - 16*h_r*ct*
     +    cb*m1*msb2**2 - 32*h_r*ct*cb*m1*msb2**4*s**(-1) + 16*h_r*ct*
     +    cb*m1*mst2**2*s*u1**(-1) + 16*h_r*ct*cb*m1*mst2**2 + 16*h_r*
     +    ct*cb*m1*s**(-1)*t1*u1 + 16*h_r*ct*cb*m1*t1 + 32*h_r*ct*cb*
     +    m1**3*mg**2*u1**(-1) - 32*h_r*ct*cb*m1**3*msb2**2*u1**(-1) - 
     +    16*h_r*ct*cb*m1**3*s*u1**(-1) - 16*h_r*ct*cb*m1**3 )

c         decouple the heavy flavors from alpha_s
      MM_s = MM_s + born * prefac * 4.D0*Pi*alphas *
     &      (  (Ns-2.D0)/3.D0  * log(qr**2/msx**2)
     &       + 1.D0/6.D0       * log(qr**2/mst1**2)
     &       + 1.D0/6.D0       * log(qr**2/mst2**2)
     &       + 1.D0/6.D0       * log(qr**2/msb1**2)
     &       + 1.D0/6.D0       * log(qr**2/msb2**2)
     &       + 2.D0/3.D0 * Nc  * log(qr**2/mg**2)   )
ctp     &       + 2.D0/3.D0       * log(qr**2/mt**2) 

c         decouple the heavy flavors from running mass
      msusy = ( mg+msb1+msb2+mst1+mst2 )/5.D0

      MM_s = MM_s + born * prefac * 4.D0*Pi*alphas *
     &              2.D0 * Cf  * log(qr**2/msusy**2)   

c               the phase space except for 1/s**2 
      HT_QGS = MM_s / ( 16.D0 * pi )

c               the averaging factors
      HT_QGS = HT_QGS /4.D0 /Nc /(Nc**2-1.D0)

c               the luminosity
      HT_QGS = HT_QGS *lumi(1)

      return
      end



