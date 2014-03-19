cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     HH_QGH(MASSIN,C)                                                 c
c                                                                      c
c     DYFACT and DYFACU NEDDED FROM Xmatrix_hh_qb.f                    c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = t2                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(6)  = m2                                                c
c       MASSIN(7)  = mt                                                c
c       MASSIN(8)  = mz                                                c
c       MASSIN(9)  = mh1                                               c
c       MASSIN(10) = mh2                                               c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c                                                                      c
c       C(1:20)  ALL COUPLINGS                                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function HH_QGOS(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Nc,Cf,alphas
     &          ,m1,m2,mt
     &          ,hl,hr 
     &          ,s,s4,t2,u2
     &          ,hardfac,theta_s3
     &          ,MMqgos,ANGfin(0:12,0:12,-2:2,-2:2)

      Pi    = 4.D0*atan(1.D0)
      Nc    = 3.D0
      Cf    = 4.D0/3.D0

      s   = massin(1)
      t2  = massin(2)
      s4  = massin(3)
      m1  = massin(6)
      m2  = massin(6)
      mt  = massin(7)

c               real kinematics built in
      u2  = s4 - s - t2 - m1**2 + m1**2 

      hl      = C(4)
      hr      = C(5)

c               set gs=1 
      alphas = 1.D0/(4.D0*Pi) 

      hardfac = 1.D0

c               the s3 regularization
      theta_s3 = 0.D0
      if ((m2.lt.mt).and.(s.gt.(m1+mt)**2)) theta_s3 = 1.D0 

c               the angular functions 
      call ANGULAR_ARRAY_HH_QGOS(massin,theta_s3,ANGfin)

c               form output
      MMqgos =
     &  + ANGfin(8,0,-2,0)*Nc*Cf*Pi*alphas*hardfac * ( 8*hl**2*hr**2*
     &    m1**2*mt**2*s**(-1) - 4*hl**2*hr**2*mt**2*s**(-1)*u2 - 8*
     &    hl**2*hr**2*mt**2 - 8*hl**2*hr**2*mt**4*s**(-1) + 4*hl**4*
     &    m1**2*mt**2*s**(-1) + 4*hl**4*m1**2 - 4*hl**4*m1**4*s**(-1)
     &     + 2*hl**4*mt**2*s**(-1)*u2 + 4*hr**4*m1**2*mt**2*s**(-1) + 4
     &    *hr**4*m1**2 - 4*hr**4*m1**4*s**(-1) + 2*hr**4*mt**2*s**(-1)*
     &    u2 )
      MMqgos = MMqgos + ANGfin(8,4,-2,-2)*Nc*Cf*Pi*alphas*hardfac * ( 
     &     - 4*hl**2*hr**2*m1**2*mt**2*s + 4*hl**2*hr**2*m1**2*mt**2*t2
     &     + 8*hl**2*hr**2*m1**2*mt**4 - 4*hl**2*hr**2*mt**4*s - 4*
     &    hl**2*hr**2*mt**4*t2 - 8*hl**2*hr**2*mt**4*u2 - 8*hl**2*hr**2
     &    *mt**6 + 2*hl**4*m1**2*mt**2*s - 2*hl**4*m1**2*mt**2*t2 + 4*
     &    hl**4*m1**2*mt**4 - 4*hl**4*m1**4*mt**2 + 2*hl**4*mt**4*s + 2
     &    *hl**4*mt**4*t2 + 4*hl**4*mt**4*u2 + 2*hr**4*m1**2*mt**2*s - 
     &    2*hr**4*m1**2*mt**2*t2 + 4*hr**4*m1**2*mt**4 - 4*hr**4*m1**4*
     &    mt**2 + 2*hr**4*mt**4*s + 2*hr**4*mt**4*t2 + 4*hr**4*mt**4*u2
     &     )
      MMqgos = MMqgos + ANGfin(8,4,-2,-1)*Nc*Cf*Pi*alphas*hardfac * ( 4
     &    *hl**2*hr**2*m1**2*mt**2*s**(-1)*t2 + 8*hl**2*hr**2*m1**2*
     &    mt**2*s**(-1)*u2 + 4*hl**2*hr**2*m1**2*mt**2 + 16*hl**2*hr**2
     &    *m1**2*mt**4*s**(-1) - 8*hl**2*hr**2*m1**4*mt**2*s**(-1) - 4*
     &    hl**2*hr**2*mt**2*s - 4*hl**2*hr**2*mt**2*u2 - 4*hl**2*hr**2*
     &    mt**4*s**(-1)*t2 - 8*hl**2*hr**2*mt**4*s**(-1)*u2 - 12*hl**2*
     &    hr**2*mt**4 - 8*hl**2*hr**2*mt**6*s**(-1) - 2*hl**4*m1**2*
     &    mt**2*s**(-1)*t2 - 4*hl**4*m1**2*mt**2*s**(-1)*u2 + 6*hl**4*
     &    m1**2*mt**2 + 4*hl**4*m1**2*mt**4*s**(-1) + 2*hl**4*m1**2*s
     &     - 8*hl**4*m1**4*mt**2*s**(-1) - 4*hl**4*m1**4 + 4*hl**4*
     &    m1**6*s**(-1) + 2*hl**4*mt**2*u2 + 2*hl**4*mt**4*s**(-1)*t2
     &     + 4*hl**4*mt**4*s**(-1)*u2 + 2*hl**4*mt**4 - 2*hr**4*m1**2*
     &    mt**2*s**(-1)*t2 - 4*hr**4*m1**2*mt**2*s**(-1)*u2 + 6*hr**4*
     &    m1**2*mt**2 + 4*hr**4*m1**2*mt**4*s**(-1) + 2*hr**4*m1**2*s
     &     - 8*hr**4*m1**4*mt**2*s**(-1) - 4*hr**4*m1**4 + 4*hr**4*
     &    m1**6*s**(-1) )
      MMqgos = MMqgos + ANGfin(8,4,-2,-1)*Nc*Cf*Pi*alphas*hardfac * ( 2
     &    *hr**4*mt**2*u2 + 2*hr**4*mt**4*s**(-1)*t2 + 4*hr**4*mt**4*
     &    s**(-1)*u2 + 2*hr**4*mt**4 )
      MMqgos = MMqgos + ANGfin(8,4,-2,1)*Nc*Cf*Pi*alphas*hardfac * ( 
     &     - 4*hl**2*hr**2*mt**2*s**(-1) + 2*hl**4*m1**2*s**(-1) + 2*
     &    hr**4*m1**2*s**(-1) )

c               the phase space except for 1/s**2 
      HH_QGOS = MMqgos / ( 16.D0 * pi**2 )**2 / 2.D0*s4/(s4+m1**2)

c               the averaging factors
      HH_QGOS = HH_QGOS /4.D0 /Nc/(Nc**2-1.D0)

c               the prefactor for the scaling functions 
      HH_QGOS = HH_QGOS * (m1+m2)**2/4.D0 

      end

c --------------------------------------------------------------------
      real*8 function HH_QGH(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Nc,Cf,sqrt2,alphas
     &          ,m1,m2,mt,mz,mh1,mh2,m12,mt2,mz2,mh12,mh22
     &          ,ssp,ssz,hl,hr 
     &          ,h1,h2,lambda1,lambda2
     &          ,lq,rq,pq
     &          ,lq2,rq2,pq2
     &          ,s,s4,t2,topfac,s3fac,u2
     &          ,hardfac,logall,theta_s3
     &          ,dyfact,dyfacu
     &          ,MMcrossed5,ANGfin(0:12,0:12,-2:2,-2:2)
     &          ,pt2

      external dyfact, dyfacu

c          common block to talk to DYFACT and DYFACU
      common/HH_QBB_INT/s,t2,u2

      Pi    = 4.D0*atan(1.D0)
      sqrt2 = sqrt(2.D0)
      Nc    = 3.D0
      Cf    = 4.D0/3.D0

      s   = massin(1)
      t2  = massin(2)
      s4  = massin(3)
      m1  = massin(6)
      m2  = massin(6)
      mt  = massin(7)
      mz  = massin(8)
      mh1 = massin(9)
      mh2 = massin(10)

c               real kinematics built in
      u2  = s4 - s - t2 - m1**2 + m1**2 

      lq      = C(1)
      rq      = C(2)
      pq      = C(3)

      hl      = C(4)
      hr      = C(5)
      h1      = C(6) 
      h2      = C(7)

      ssp     = C(8)
      ssz     = C(9) 
      lambda1 = C(10)
      lambda2 = C(11)

      lq2 = lq**2 
      rq2 = rq**2 
      pq2 = pq**2 

c               mass squares more convenient
      m12  = m1**2
      mt2  = mt**2
      mz2  = mz**2
      mh12 = mh1**2
      mh22 = mh2**2

c               logall = logqf - log(s4/m1^2) + log(1+m1^2/s4)
ctp      logall = log(massin(13)**2) 
      logall = log(massin(13)**2/m1**2) 
     &        - log(s4/m1**2) + log(1.D0+m1**2/s4)

c               set gs=1 
      alphas = 1.D0/(4.D0*Pi) 

      hardfac = 1.D0

c               remaining special function without argument 
      topfac = 1.D0 + (mt2-m12)*(s+t2)/t2/u2
      topfac = 1.D0/topfac

      s3fac  = 1.D0 + (mt2-m12)*(s+t2)/t2/s4
      s3fac  = 1.D0/s3fac

      pt2 = (t2*u2-s*m12)/s 

c               the s3 regularization
      theta_s3 = 0.D0
      if ((m2.lt.mt).and.(s.gt.(m1+mt)**2)) theta_s3 = 1.D0 

c               the angular functions 
      call ANGULAR_ARRAY_HH_QG(massin,theta_s3,ANGfin)

c               form output
      MMcrossed5 =
     &  + (1+m12/s4)*pq*ssp*Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac * (  - 
     &    16*s3fac*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2 - 16*
     &    s3fac*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1) - 32*
     &    s3fac*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1) + 32*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*t2**(-1)*u2 + 48*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**2*t2**(-1) + 48*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2 + 32*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*u2 - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*s4**(-1)*Pi**2*
     & alphas*hardfac * (  - 32*s3fac*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2 + 16*s3fac*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*
     &    t2 + 16*s3fac*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1)
     &     + 32*s3fac*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2 + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s*t2**(-1)*u2 + 48*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s + 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**2*t2**(-1) + 48*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2 + 32*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*u2 )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*s4*Pi**2*alphas
     & *hardfac * (  - 16*(u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*
     &    u2**(-1) - 16*(u2-m12+mt2)**(-1)*hl**2*u2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*u2**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*u2**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*u2**(-1) - 128*lq*ssz*mz**(-2)*
     &    s**(-1)*t2*u2**(-1) - 128*lq*ssz*mz**(-2)*u2**(-1) - 128*rq*
     &    ssz*mz**(-2)*s**(-1)*t2*u2**(-1) - 128*rq*ssz*mz**(-2)*
     &    u2**(-1) + 128*dyfact(mz)*(s+t2)*lq*ssz*mz**(-2)*s**(-1)*
     &    u2**(-1) + 128*dyfact(mz)*(s+t2)*rq*ssz*mz**(-2)*s**(-1)*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2)**(-1)*hl**2*s*t2**(-1)*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2)**(-1)*hl**2*u2**(-1) - 16
     &    *topfac*(u2-m12+mt2)**(-1)*hr**2*s*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*s4*Pi**2*alphas
     & *hardfac * (  - 16*topfac*(u2-m12+mt2)**(-1)*hr**2*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac*logall * (  - 256*dyfact(mz)*(s+t2)**(-1)*lq*ssz*pt2*
     &    s**(-1) - 256*dyfact(mz)*(s+t2)**(-1)*rq*ssz*pt2*s**(-1) - 
     &    128*dyfact(mz)*(s+t2)*lq*ssz*pt2*s**(-1)*u2**(-2) - 128*
     &    dyfact(mz)*(s+t2)*rq*ssz*pt2*s**(-1)*u2**(-2) - 256*dyfact(mz
     &    )*lq*ssz*pt2*s**(-1)*u2**(-1) - 256*dyfact(mz)*rq*ssz*pt2*
     &    s**(-1)*u2**(-1) - 64*topfac*(s+t2)**(-1)*hl**2*pt2*t2**(-1)
     &     - 64*topfac*(s+t2)**(-1)*hr**2*pt2*t2**(-1) - 32*topfac*
     &    (s+t2)*hl**2*pt2*t2**(-1)*u2**(-2) - 32*topfac*(s+t2)*hr**2*
     &    pt2*t2**(-1)*u2**(-2) - 64*topfac*hl**2*pt2*t2**(-1)*u2**(-1)
     &     - 64*topfac*hr**2*pt2*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac * ( 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2*
     &    u2**(-1) - 32*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2*
     &    u2**(-1) - 32*(u2-m12+mt2)**(-1)*hl**2*mt**2*u2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2*u2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*s*u2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*t2*u2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2*u2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2*u2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*u2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2*u2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*s*u2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*t2*u2**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2*u2**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*u2**(-1) + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac * ( 16*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2*
     &    u2**(-1) + 16*(t2+u2-m12+mt2)**(-1)*hr**2*s*u2**(-1) + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2*u2**(-1) - 128*lq*ssz*m1**2*
     &    mz**(-2)*u2**(-1) + 128*lq*ssz*mz**(-2)*s**(-1)*t2 + 128*lq*
     &    ssz*mz**(-2)*s**(-1)*t2**2*u2**(-1) + 128*lq*ssz*mz**(-2)*s*
     &    u2**(-1) + 256*lq*ssz*mz**(-2)*t2*u2**(-1) + 128*lq*ssz*
     &    mz**(-2) - 128*rq*ssz*m1**2*mz**(-2)*u2**(-1) + 128*rq*ssz*
     &    mz**(-2)*s**(-1)*t2 + 128*rq*ssz*mz**(-2)*s**(-1)*t2**2*
     &    u2**(-1) + 128*rq*ssz*mz**(-2)*s*u2**(-1) + 256*rq*ssz*
     &    mz**(-2)*t2*u2**(-1) + 128*rq*ssz*mz**(-2) + 32*hl**2*s**(-1)
     &    *t2*u2**(-1) + 32*hl**2*u2**(-1) + 32*hr**2*s**(-1)*t2*
     &    u2**(-1) + 32*hr**2*u2**(-1) - 256*dyfact(mz)*(s+t2)**(-1)*lq
     &    *ssz*pt2*s**(-1) - 256*dyfact(mz)*(s+t2)**(-1)*rq*ssz*pt2*
     &    s**(-1) - 128*dyfact(mz)*(s+t2)*lq*ssz*pt2*s**(-1)*u2**(-2)
     &     - 128*dyfact(mz)*(s+t2)*lq*ssz*mz**(-2)*s**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac * ( 128*dyfact(mz)*(s+t2)*lq*ssz*mz**(-2)*u2**(-1) - 128
     &    *dyfact(mz)*(s+t2)*lq*ssz*s**(-1)*u2**(-1) - 128*dyfact(mz)*
     &    (s+t2)*rq*ssz*pt2*s**(-1)*u2**(-2) - 128*dyfact(mz)*(s+t2)*rq
     &    *ssz*mz**(-2)*s**(-1) + 128*dyfact(mz)*(s+t2)*rq*ssz*mz**(-2)
     &    *u2**(-1) - 128*dyfact(mz)*(s+t2)*rq*ssz*s**(-1)*u2**(-1) - 
     &    128*dyfact(mz)*(s+t2)**2*lq*ssz*mz**(-2)*s**(-1)*u2**(-1) + 
     &    128*dyfact(mz)*(s+t2)**2*lq*ssz*s**(-2)*u2**(-1) - 128*
     &    dyfact(mz)*(s+t2)**2*rq*ssz*mz**(-2)*s**(-1)*u2**(-1) + 128*
     &    dyfact(mz)*(s+t2)**2*rq*ssz*s**(-2)*u2**(-1) - 256*dyfact(mz)
     &    *lq*ssz*pt2*s**(-1)*u2**(-1) + 128*dyfact(mz)*lq*ssz*m1**2*
     &    mz**(-2)*u2**(-1) - 128*dyfact(mz)*lq*ssz*mz**(-2)*s*u2**(-1)
     &     - 128*dyfact(mz)*lq*ssz*mz**(-2)*t2*u2**(-1) - 256*dyfact(mz
     &    )*rq*ssz*pt2*s**(-1)*u2**(-1) + 128*dyfact(mz)*rq*ssz*m1**2*
     &    mz**(-2)*u2**(-1) - 128*dyfact(mz)*rq*ssz*mz**(-2)*s*u2**(-1)
     &     - 128*dyfact(mz)*rq*ssz*mz**(-2)*t2*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac * (  - 64*topfac*(s+t2)**(-1)*hl**2*pt2*t2**(-1) - 64*
     &    topfac*(s+t2)**(-1)*hr**2*pt2*t2**(-1) - 32*topfac*(s+t2)*
     &    hl**2*pt2*t2**(-1)*u2**(-2) - 32*topfac*(s+t2)*hr**2*pt2*
     &    t2**(-1)*u2**(-2) - 16*topfac*(u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*s**(-1)*t2*u2**(-1) - 16*topfac*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2*u2**(-1) - 16*
     &    topfac*(u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1)*u2**(-1)
     &     + 16*topfac*(u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 + 16*
     &    topfac*(u2-m12+mt2)**(-1)*hl**2*s*t2**(-1) + 32*topfac*
     &    (u2-m12+mt2)**(-1)*hl**2*s*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2)**(-1)*hl**2*s**2*t2**(-1)*u2**(-1) + 16*topfac
     &    *(u2-m12+mt2)**(-1)*hl**2*t2*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2)**(-1)*hl**2 - 16*topfac*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**(-1)*t2*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*pq*ssp*Nc*Cf*Pi**2*alphas*
     & hardfac * (  - 16*topfac*(u2-m12+mt2)**(-1)*hr**2*m1**2*s*
     &    t2**(-1)*u2**(-1) + 16*topfac*(u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*s**(-1)*t2*u2**(-1) - 16*topfac*(u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2 + 16*topfac*
     &    (u2-m12+mt2)**(-1)*hr**2*s*t2**(-1) + 32*topfac*
     &    (u2-m12+mt2)**(-1)*hr**2*s*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2)**(-1)*hr**2*s**2*t2**(-1)*u2**(-1) + 16*topfac
     &    *(u2-m12+mt2)**(-1)*hr**2*t2*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2)**(-1)*hr**2 - 64*topfac*hl**2*pt2*t2**(-1)*
     &    u2**(-1) - 64*topfac*hr**2*pt2*t2**(-1)*u2**(-1) - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*t2**(-1) - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s*t2**(-1) - 16*s3fac*
     &    (t2+u2-m12+mt2)**(-1)*hr**2 )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*ssp**2*pq2*Nc*Cf*s4*Pi**2*
     & alphas*hardfac * (  - 128*s**(-2)*t2**2*u2**(-2) - 256*s**(-1)*
     &    t2*u2**(-2) - 128*u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*ssp**2*pq2*Nc*Cf*Pi**2*
     & alphas*hardfac*logall * (  - 256*(s+t2)**(-1)*pt2*s**(-1) - 128*
     &    (s+t2)*pt2*s**(-1)*u2**(-2) - 256*pt2*s**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*ssp**2*pq2*Nc*Cf*Pi**2*
     & alphas*hardfac * (  - 256*(s+t2)**(-1)*pt2*s**(-1) - 128*(s+t2)*
     &    pt2*s**(-1)*u2**(-2) - 256*pt2*s**(-1)*u2**(-1) - 128*m1**2*
     &    s**(-1)*t2*u2**(-2) - 128*m1**2*u2**(-2) + 256*s**(-2)*t2**2*
     &    u2**(-1) + 128*s**(-2)*t2**3*u2**(-2) + 384*s**(-1)*t2*
     &    u2**(-1) + 384*s**(-1)*t2**2*u2**(-2) + 128*s*u2**(-2) + 384*
     &    t2*u2**(-2) + 128*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*s4**(-1)*Pi**2*alphas*
     & hardfac * (  - 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s**(-1)*t2 - 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*s*t2**(-1) - 32*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2 + 
     &    16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*
     &    t2**(-1) + 32*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**2 + 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2*u2 + 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*s**(-1)*t2**2 + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2**(-1)*u2 + 48*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2*t2**(-1) + 48*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 32*s3fac
     &    *(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*s4**(-1)*Pi**2*alphas*
     & hardfac * (  - 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s**(-1)*t2 - 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*m1**2*s*t2**(-1) - 32*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2 + 
     &    16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*
     &    t2**(-1) + 32*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**2 + 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*t2*u2 + 16*s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*s**(-1)*t2**2 + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2**(-1)*u2 + 48*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s + 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2*t2**(-1) + 48*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 32*s3fac
     &    *(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*s4**(-1)*Pi**2*alphas*
     & hardfac * ( 4*s3fac*hl**4*m1**2*s**(-1) + 4*s3fac*hl**4*m1**2*s*
     &    t2**(-2) + 8*s3fac*hl**4*m1**2*t2**(-1) - 4*s3fac*hl**4*mt**2
     &    *s**(-1) - 4*s3fac*hl**4*mt**2*s*t2**(-2) - 8*s3fac*hl**4*
     &    mt**2*t2**(-1) - 4*s3fac*hl**4*s**(-1)*t2 - 4*s3fac*hl**4*
     &    s**(-1)*u2 - 4*s3fac*hl**4*s*t2**(-2)*u2 - 12*s3fac*hl**4*s*
     &    t2**(-1) - 4*s3fac*hl**4*s**2*t2**(-2) - 8*s3fac*hl**4*
     &    t2**(-1)*u2 - 12*s3fac*hl**4 + 4*s3fac*hr**4*m1**2*s**(-1) + 
     &    4*s3fac*hr**4*m1**2*s*t2**(-2) + 8*s3fac*hr**4*m1**2*t2**(-1)
     &     - 4*s3fac*hr**4*mt**2*s**(-1) - 4*s3fac*hr**4*mt**2*s*
     &    t2**(-2) - 8*s3fac*hr**4*mt**2*t2**(-1) - 4*s3fac*hr**4*
     &    s**(-1)*t2 - 4*s3fac*hr**4*s**(-1)*u2 - 4*s3fac*hr**4*s*
     &    t2**(-2)*u2 - 12*s3fac*hr**4*s*t2**(-1) - 4*s3fac*hr**4*s**2*
     &    t2**(-2) - 8*s3fac*hr**4*t2**(-1)*u2 - 12*s3fac*hr**4 )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*s4*Pi**2*alphas*
     & hardfac * (  - 16*dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)*lq
     &    *hl**2*ssz*s**(-1)*u2**(-1) - 16*dyfact(mz)*(s+t2)*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2**(-1) - 16*
     &    dyfact(mz)*(s+t2)*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*u2**(-1) - 16*dyfact(mz)*(s+t2)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2**(-1) - 
     &    64*dyfact(mz)**2*(s+t2)**2*ssz**2*lq2*s**(-2)*u2**(-2) - 64*
     &    dyfact(mz)**2*(s+t2)**2*ssz**2*rq2*s**(-2)*u2**(-2) - 16*
     &    topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2**(-1)*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s
     &    *t2**(-1)*u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*u2**(-1) - 4*topfac*hl**4*s*t2**(-2)*u2**(-1) - 4*
     &    topfac*hl**4*t2**(-1)*u2**(-1) - 4*topfac*hr**4*s*t2**(-2)*
     &    u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*s4*Pi**2*alphas*
     & hardfac * (  - 4*topfac*hr**4*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac*
     & logall * (  - 128*dyfact(mz)**2*(s+t2)**(-1)*ssz**2*lq2*pt2*
     &    s**(-1) - 128*dyfact(mz)**2*(s+t2)**(-1)*ssz**2*rq2*pt2*
     &    s**(-1) - 64*dyfact(mz)**2*(s+t2)*ssz**2*lq2*pt2*s**(-1)*
     &    u2**(-2) - 64*dyfact(mz)**2*(s+t2)*ssz**2*rq2*pt2*s**(-1)*
     &    u2**(-2) - 128*dyfact(mz)**2*ssz**2*lq2*pt2*s**(-1)*u2**(-1)
     &     - 128*dyfact(mz)**2*ssz**2*rq2*pt2*s**(-1)*u2**(-1) - 64*
     &    dyfact(mz)*topfac*(s+t2)**(-1)*lq*hl**2*ssz*pt2*t2**(-1) - 64
     &    *dyfact(mz)*topfac*(s+t2)**(-1)*rq*hr**2*ssz*pt2*t2**(-1) - 
     &    32*dyfact(mz)*topfac*(s+t2)*lq*hl**2*ssz*pt2*t2**(-1)*
     &    u2**(-2) - 32*dyfact(mz)*topfac*(s+t2)*rq*hr**2*ssz*pt2*
     &    t2**(-1)*u2**(-2) - 64*dyfact(mz)*topfac*lq*hl**2*ssz*pt2*
     &    t2**(-1)*u2**(-1) - 64*dyfact(mz)*topfac*rq*hr**2*ssz*pt2*
     &    t2**(-1)*u2**(-1) - 32*dyfact(mh1)**2*(s+t2)**(-1)*h1**2*
     &    lambda1**2*s**(-1) - 16*dyfact(mh1)**2*(s+t2)*h1**2*
     &    lambda1**2*s**(-1)*u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac*
     & logall * (  - 32*dyfact(mh1)**2*h1**2*lambda1**2*s**(-1)*
     &    u2**(-1) - 64*dyfact(mh1)*dyfact(mh2)*(s+t2)**(-1)*h1*h2*
     &    lambda1*lambda2*s**(-1) - 32*dyfact(mh1)*dyfact(mh2)*(s+t2)*
     &    h1*h2*lambda1*lambda2*s**(-1)*u2**(-2) - 64*dyfact(mh1)*
     &    dyfact(mh2)*h1*h2*lambda1*lambda2*s**(-1)*u2**(-1) - 64*
     &    dyfact(mh1)*topfac*(s+t2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt*t2**(-1) - 32*dyfact(mh1)*topfac*(s+t2)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*t2**(-1)*u2**(-2) - 64*dyfact(mh1)*topfac*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*u2**(-1) - 32*dyfact(
     &    mh2)**2*(s+t2)**(-1)*h2**2*lambda2**2*s**(-1) - 16*dyfact(mh2
     &    )**2*(s+t2)*h2**2*lambda2**2*s**(-1)*u2**(-2) - 32*dyfact(mh2
     &    )**2*h2**2*lambda2**2*s**(-1)*u2**(-1) - 64*dyfact(mh2)*
     &    topfac*(s+t2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2**(-1)
     &     - 32*dyfact(mh2)*topfac*(s+t2)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*t2**(-1)*u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac*
     & logall * (  - 64*dyfact(mh2)*topfac*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt*t2**(-1)*u2**(-1) - 16*topfac**2*(s+t2)**(-1)*hl**2*hr**2
     &    *mt**2*s*t2**(-2) - 8*topfac**2*(s+t2)**(-1)*hl**4*pt2*s*
     &    t2**(-2) - 8*topfac**2*(s+t2)**(-1)*hr**4*pt2*s*t2**(-2) - 8*
     &    topfac**2*(s+t2)*hl**2*hr**2*mt**2*s*t2**(-2)*u2**(-2) - 4*
     &    topfac**2*(s+t2)*hl**4*pt2*s*t2**(-2)*u2**(-2) - 4*topfac**2*
     &    (s+t2)*hr**4*pt2*s*t2**(-2)*u2**(-2) - 16*topfac**2*hl**2*
     &    hr**2*mt**2*s*t2**(-2)*u2**(-1) - 8*topfac**2*hl**4*pt2*s*
     &    t2**(-2)*u2**(-1) - 8*topfac**2*hr**4*pt2*s*t2**(-2)*u2**(-1)
     &     )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * ( 32*dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**2*s**(-1)*u2**(-1) - 32*dyfact(mz)*(s+t2)*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*u2**(-1)
     &     - 16*dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*u2**(-1) + 32*dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)
     &    *rq*hr**2*ssz*m1**2*s**(-1)*u2**(-1) - 32*dyfact(mz)*(s+t2)*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*u2**(-1)
     &     - 16*dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*u2**(-1) - 16*dyfact(mz)*(s+t2)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2**(-1) - 16*
     &    dyfact(mz)*(s+t2)*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    u2**(-1) + 32*dyfact(mz)*(s+t2)*lq*hl**2*ssz*s**(-1)*u2**(-1)
     &     + 32*dyfact(mz)*(s+t2)*rq*hr**2*ssz*s**(-1)*u2**(-1) - 64*
     &    dyfact(mz)*(s+t2)*ssz**2*lq2*s**(-1)*u2**(-1) - 64*dyfact(mz)
     &    *(s+t2)*ssz**2*rq2*s**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 16*dyfact(mz)*(s+t2)**2*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mz**2*s**(-2)*u2**(-1) + 16*dyfact(mz)*(s+t2)**2*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*u2**(-1) - 16*
     &    dyfact(mz)*(s+t2)**2*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mz**2*s**(-2)*u2**(-1) + 16*dyfact(mz)*(s+t2)**2*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2**(-1) - 16*
     &    dyfact(mz)*(s+t2)**2*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mz**2*s**(-2)*u2**(-1) + 16*dyfact(mz)*(s+t2)**2*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*u2**(-1) - 
     &    16*dyfact(mz)*(s+t2)**2*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*mz**2*s**(-2)*u2**(-1) + 16*dyfact(mz)*(s+t2)**2*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2**(-1) + 
     &    64*dyfact(mz)*(s+t2)**2*ssz**2*lq2*s**(-2)*u2**(-1) + 64*
     &    dyfact(mz)*(s+t2)**2*ssz**2*rq2*s**(-2)*u2**(-1) - 32*dyfact(
     &    mz)*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * ( 16*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*
     &    u2**(-1) + 16*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*t2*u2**(-1) - 32*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2**(-1) + 16*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2**(-1) + 16*
     &    dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2*u2**(-1)
     &     + 16*dyfact(mz)*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*
     &    u2**(-1) + 16*dyfact(mz)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*t2*u2**(-1) + 16*dyfact(mz)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2**(-1) + 16*
     &    dyfact(mz)*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2*
     &    u2**(-1) - 128*dyfact(mz)**2*(s+t2)**(-1)*ssz**2*lq2*pt2*
     &    s**(-1) - 128*dyfact(mz)**2*(s+t2)**(-1)*ssz**2*rq2*pt2*
     &    s**(-1) - 64*dyfact(mz)**2*(s+t2)*ssz**2*lq2*pt2*s**(-1)*
     &    u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 64*dyfact(mz)**2*(s+t2)*ssz**2*lq2*m1**2*s**(-1)*
     &    u2**(-2) - 64*dyfact(mz)**2*(s+t2)*ssz**2*lq2*u2**(-2) - 64*
     &    dyfact(mz)**2*(s+t2)*ssz**2*rq2*pt2*s**(-1)*u2**(-2) - 64*
     &    dyfact(mz)**2*(s+t2)*ssz**2*rq2*m1**2*s**(-1)*u2**(-2) - 64*
     &    dyfact(mz)**2*(s+t2)*ssz**2*rq2*u2**(-2) + 64*dyfact(mz)**2*
     &    (s+t2)**2*ssz**2*lq2*mz**2*s**(-2)*u2**(-2) + 64*dyfact(mz)**
     &    2*(s+t2)**2*ssz**2*lq2*s**(-2)*u2**(-1) + 64*dyfact(mz)**2*
     &    (s+t2)**2*ssz**2*rq2*mz**2*s**(-2)*u2**(-2) + 64*dyfact(mz)**
     &    2*(s+t2)**2*ssz**2*rq2*s**(-2)*u2**(-1) - 64*dyfact(mz)**2*
     &    (s+t2)**3*ssz**2*lq2*mz**2*s**(-3)*u2**(-2) + 64*dyfact(mz)**
     &    2*(s+t2)**3*ssz**2*lq2*s**(-2)*u2**(-2) - 64*dyfact(mz)**2*
     &    (s+t2)**3*ssz**2*rq2*mz**2*s**(-3)*u2**(-2) + 64*dyfact(mz)**
     &    2*(s+t2)**3*ssz**2*rq2*s**(-2)*u2**(-2) - 128*dyfact(mz)**2*
     &    ssz**2*lq2*pt2*s**(-1)*u2**(-1) + 64*dyfact(mz)**2*ssz**2*lq2
     &    *s*u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * ( 64*dyfact(mz)**2*ssz**2*lq2*t2*u2**(-2) - 128*dyfact(mz)**2
     &    *ssz**2*rq2*pt2*s**(-1)*u2**(-1) + 64*dyfact(mz)**2*ssz**2*
     &    rq2*s*u2**(-2) + 64*dyfact(mz)**2*ssz**2*rq2*t2*u2**(-2) - 64
     &    *dyfact(mz)*topfac*(s+t2)**(-1)*lq*hl**2*ssz*pt2*t2**(-1) - 
     &    64*dyfact(mz)*topfac*(s+t2)**(-1)*rq*hr**2*ssz*pt2*t2**(-1)
     &     - 32*dyfact(mz)*topfac*(s+t2)*lq*hl**2*ssz*pt2*t2**(-1)*
     &    u2**(-2) - 32*dyfact(mz)*topfac*(s+t2)*rq*hr**2*ssz*pt2*
     &    t2**(-1)*u2**(-2) - 64*dyfact(mz)*topfac*lq*hl**2*ssz*pt2*
     &    t2**(-1)*u2**(-1) - 64*dyfact(mz)*topfac*rq*hr**2*ssz*pt2*
     &    t2**(-1)*u2**(-1) - 32*dyfact(mh1)*(mh12-mh22)**(-1)*h1*h2*
     &    lambda1*lambda2*u2**(-1) + 32*dyfact(mh1)*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    u2**(-1) - 32*dyfact(mh1)**2*(s+t2)**(-1)*h1**2*lambda1**2*
     &    s**(-1) - 32*dyfact(mh1)**2*h1**2*lambda1**2*s**(-1)*u2**(-1)
     &     - 64*dyfact(mh1)*dyfact(mh2)*(s+t2)**(-1)*h1*h2*lambda1*
     &    lambda2*s**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 32*dyfact(mh1)*dyfact(mh2)*(s+t2)*h1*h2*lambda1*lambda2*
     &    s**(-1)*u2**(-2) - 64*dyfact(mh1)*dyfact(mh2)*h1*h2*lambda1*
     &    lambda2*s**(-1)*u2**(-1) - 64*dyfact(mh1)*topfac*(s+t2)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1) - 32*dyfact(mh1)*
     &    topfac*(s+t2)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*
     &    u2**(-2) - 64*dyfact(mh1)*topfac*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt*t2**(-1)*u2**(-1) + 32*dyfact(mh2)*(mh12-mh22)**(-1)*h1
     &    *h2*lambda1*lambda2*u2**(-1) + 32*dyfact(mh2)*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    u2**(-1) - 32*dyfact(mh2)**2*(s+t2)**(-1)*h2**2*lambda2**2*
     &    s**(-1) - 32*dyfact(mh2)**2*h2**2*lambda2**2*s**(-1)*u2**(-1)
     &     - 64*dyfact(mh2)*topfac*(s+t2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2**(-1) - 32*dyfact(mh2)*topfac*(s+t2)*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt*t2**(-1)*u2**(-2) - 64*dyfact(mh2)*
     &    topfac*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 16*topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*t2*u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*t2 + 16*topfac
     &    *(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2**(-1) + 32*
     &    topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*u2**(-1) + 16
     &    *topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2*t2**(-1)*
     &    u2**(-1) + 16*topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    t2*u2**(-1) + 32*topfac*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz - 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*t2*u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * ( 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**(-1)*t2*u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2 + 16*topfac
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2**(-1) + 32*
     &    topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2**(-1) + 16
     &    *topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2*t2**(-1)*
     &    u2**(-1) + 16*topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    t2*u2**(-1) + 32*topfac*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz + 32*topfac*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2**(-1)*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*
     &    t2**(-1)*u2**(-1) - 4*topfac*hl**4*m1**2*s**(-1)*u2**(-1) - 4
     &    *topfac*hl**4*m1**2*s*t2**(-2)*u2**(-1) - 8*topfac*hl**4*
     &    m1**2*t2**(-1)*u2**(-1) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * ( 4*topfac*hl**4*mt**2*s**(-1)*u2**(-1) + 4*topfac*hl**4*
     &    mt**2*s*t2**(-2)*u2**(-1) + 8*topfac*hl**4*mt**2*t2**(-1)*
     &    u2**(-1) + 4*topfac*hl**4*s**(-1) + 4*topfac*hl**4*s*t2**(-2)
     &     + 12*topfac*hl**4*s*t2**(-1)*u2**(-1) + 4*topfac*hl**4*s**2*
     &    t2**(-2)*u2**(-1) + 8*topfac*hl**4*t2**(-1) + 8*topfac*hl**4*
     &    u2**(-1) - 4*topfac*hr**4*m1**2*s**(-1)*u2**(-1) - 4*topfac*
     &    hr**4*m1**2*s*t2**(-2)*u2**(-1) - 8*topfac*hr**4*m1**2*
     &    t2**(-1)*u2**(-1) + 4*topfac*hr**4*mt**2*s**(-1)*u2**(-1) + 4
     &    *topfac*hr**4*mt**2*s*t2**(-2)*u2**(-1) + 8*topfac*hr**4*
     &    mt**2*t2**(-1)*u2**(-1) + 4*topfac*hr**4*s**(-1) + 4*topfac*
     &    hr**4*s*t2**(-2) + 12*topfac*hr**4*s*t2**(-1)*u2**(-1) + 4*
     &    topfac*hr**4*s**2*t2**(-2)*u2**(-1) + 8*topfac*hr**4*t2**(-1)
     &     + 8*topfac*hr**4*u2**(-1) - 16*topfac**2*(s+t2)**(-1)*hl**2*
     &    hr**2*mt**2*s*t2**(-2) - 8*topfac**2*(s+t2)**(-1)*hl**4*pt2*s
     &    *t2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 8*topfac**2*(s+t2)**(-1)*hr**4*pt2*s*t2**(-2) - 8*topfac
     &    **2*(s+t2)*hl**2*hr**2*mt**2*s*t2**(-2)*u2**(-2) - 4*topfac**
     &    2*(s+t2)*hl**4*pt2*s*t2**(-2)*u2**(-2) - 4*topfac**2*(s+t2)*
     &    hr**4*pt2*s*t2**(-2)*u2**(-2) - 16*topfac**2*hl**2*hr**2*
     &    mt**2*s*t2**(-2)*u2**(-1) + 8*topfac**2*hl**2*hr**2*mt**2*s*
     &    t2**(-1)*u2**(-2) + 8*topfac**2*hl**2*hr**2*mt**2*s**2*
     &    t2**(-2)*u2**(-2) - 8*topfac**2*hl**4*pt2*s*t2**(-2)*u2**(-1)
     &     + 4*topfac**2*hl**4*m1**2*s*t2**(-1)*u2**(-2) + 4*topfac**2*
     &    hl**4*m1**2*u2**(-2) - 8*topfac**2*hl**4*mt**2*s*t2**(-1)*
     &    u2**(-2) - 4*topfac**2*hl**4*mt**2*s**2*t2**(-2)*u2**(-2) - 4
     &    *topfac**2*hl**4*mt**2*u2**(-2) - 8*topfac**2*hr**4*pt2*s*
     &    t2**(-2)*u2**(-1) + 4*topfac**2*hr**4*m1**2*s*t2**(-1)*
     &    u2**(-2) + 4*topfac**2*hr**4*m1**2*u2**(-2) - 8*topfac**2*
     &    hr**4*mt**2*s*t2**(-1)*u2**(-2) - 4*topfac**2*hr**4*mt**2*
     &    s**2*t2**(-2)*u2**(-2) )
      MMcrossed5 = MMcrossed5 + (1+m12/s4)*Nc*Cf*Pi**2*alphas*hardfac
     &  * (  - 4*topfac**2*hr**4*mt**2*u2**(-2) - 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2**(-1) - 16*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz - 16*s3fac*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2**(-1) - 16*
     &    s3fac*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 4*s3fac*
     &    hl**4*s*t2**(-2) + 4*s3fac*hl**4*t2**(-1) + 4*s3fac*hr**4*s*
     &    t2**(-2) + 4*s3fac*hr**4*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1) - 16
     &    *(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2 + 128*lq*ssz*s**(-1)
     &     + 128*rq*ssz*s**(-1) - 16*hl**2*s**(-1) - 16*hr**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * ( 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*mt**2 + 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*u2
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**4 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *m1**2*mt**2 + 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**4 - 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1) - 
     &    16*(u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2 + 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**4*s**(-1) - 16*(u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 16*(u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2
     &     - 16*(u2-m12+mt2)**(-1)*hr**2*m1**2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1) - 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2 + 8*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*u2 - 16*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 24*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*
     &    t2 - 16*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2 + 8*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*u2 + 128*lq*ssz*m1**2*s**(-1) - 64*lq*ssz*s**(-1)*t2 - 
     &    64*lq*ssz + 128*rq*ssz*m1**2*s**(-1) - 64*rq*ssz*s**(-1)*t2
     &     - 64*rq*ssz + 8*hl**2*s**(-1)*t2 + 8*hl**2 + 8*hr**2*s**(-1)
     &    *t2 + 8*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*ssp**2*pq2*Nc*Cf*s4*Pi
     & *alphas*hardfac * ( 128*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * (  - 64 + 128*m1**2*s**(-1) - 64*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mz**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*u2 - 16*lq*hl**2*ssz*s**(-1) - 16*rq*hr**2*ssz*
     &    s**(-1) + 4*hl**4*s**(-1) + 4*hr**4*s**(-1) + 64*ssz**2*lq2*
     &    s**(-1) + 64*ssz**2*rq2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh1**2*s**(-1) + 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2
     &    *mh2**2*s**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*u2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mz**2*s**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s**(-1) + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3 - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s**(-1) + 16*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1) + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**(-1)*u2 + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *t2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s**(-1) + 16*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s**(-1) + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*s**(-1) - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*
     &    s**(-1) - 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2
     &    *s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mz**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    t2**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*
     &    s**(-1) - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2
     &    *s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mz**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    t2**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1) + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*t2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s**(-1) - 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2*s**(-1) + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**(-1)*u2 + 16*(t2+u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*s**(-1) + 8*lq*hl**2*ssz*
     &    s**(-1)*t2 + 8*lq*hl**2*ssz + 8*rq*hr**2*ssz*s**(-1)*t2 + 8*
     &    rq*hr**2*ssz - 8*hl**2*hr**2*mt**2*s**(-1) + 4*hl**4*m1**2*
     &    s**(-1) - 2*hl**4*s**(-1)*t2 - 2*hl**4 + 4*hr**4*m1**2*
     &    s**(-1) - 2*hr**4*s**(-1)*t2 - 2*hr**4 - 16*h1**2*lambda1**2*
     &    s**(-1) - 16*h2**2*lambda2**2*s**(-1) + 64*ssz**2*lq2*m1**2*
     &    s**(-1) - 32*ssz**2*lq2*s**(-1)*t2 - 32*ssz**2*lq2 + 64*
     &    ssz**2*rq2*m1**2*s**(-1) - 32*ssz**2*rq2*s**(-1)*t2 - 32*
     &    ssz**2*rq2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*hl**2*hr**2*m1**2*mt**2*t2**(-1) - 4*hl**2*hr**2*
     &    mt**2*s*t2**(-1) - 8*hl**2*hr**2*mt**2*t2**(-1)*u2 - 8*hl**2*
     &    hr**2*mt**4*t2**(-1) + 12*hl**4*m1**2*mt**2*t2**(-1) + 8*
     &    hl**4*m1**2*s*t2**(-1) + 8*hl**4*m1**2*t2**(-1)*u2 + 4*hl**4*
     &    m1**2 - 8*hl**4*m1**4*t2**(-1) - 6*hl**4*mt**2*s*t2**(-1) - 4
     &    *hl**4*mt**2*t2**(-1)*u2 - 4*hl**4*mt**2 - 4*hl**4*mt**4*
     &    t2**(-1) + 12*hr**4*m1**2*mt**2*t2**(-1) + 8*hr**4*m1**2*s*
     &    t2**(-1) + 8*hr**4*m1**2*t2**(-1)*u2 + 4*hr**4*m1**2 - 8*
     &    hr**4*m1**4*t2**(-1) - 6*hr**4*mt**2*s*t2**(-1) - 4*hr**4*
     &    mt**2*t2**(-1)*u2 - 4*hr**4*mt**2 - 4*hr**4*mt**4*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-2)*Nc*Cf*s4**2*Pi*alphas
     & *hardfac * (  - 4*hl**4*m1**2*t2**(-1) + 4*hl**4*mt**2*t2**(-1)
     &     - 4*hr**4*m1**2*t2**(-1) + 4*hr**4*mt**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 12*hl**2*hr**2*m1**2*mt**2*s*t2**(-1) - 16*hl**2*
     &    hr**2*m1**2*mt**2*t2**(-1)*u2 - 4*hl**2*hr**2*m1**2*mt**2 - 
     &    16*hl**2*hr**2*m1**2*mt**4*t2**(-1) + 8*hl**2*hr**2*m1**4*
     &    mt**2*t2**(-1) + 12*hl**2*hr**2*mt**2*s*t2**(-1)*u2 + 4*hl**2
     &    *hr**2*mt**2*s**2*t2**(-1) + 8*hl**2*hr**2*mt**2*t2**(-1)*
     &    u2**2 + 4*hl**2*hr**2*mt**2*u2 + 12*hl**2*hr**2*mt**4*s*
     &    t2**(-1) + 16*hl**2*hr**2*mt**4*t2**(-1)*u2 + 4*hl**2*hr**2*
     &    mt**4 + 8*hl**2*hr**2*mt**6*t2**(-1) - 10*hl**4*m1**2*mt**2*s
     &    *t2**(-1) - 8*hl**4*m1**2*mt**2*t2**(-1)*u2 - 6*hl**4*m1**2*
     &    mt**2 - 4*hl**4*m1**2*mt**4*t2**(-1) - 8*hl**4*m1**2*s*
     &    t2**(-1)*u2 - 4*hl**4*m1**2*s - 4*hl**4*m1**2*s**2*t2**(-1)
     &     - 4*hl**4*m1**2*t2**(-1)*u2**2 - 2*hl**4*m1**2*t2 - 4*hl**4*
     &    m1**2*u2 + 8*hl**4*m1**4*mt**2*t2**(-1) + 8*hl**4*m1**4*s*
     &    t2**(-1) + 8*hl**4*m1**4*t2**(-1)*u2 + 4*hl**4*m1**4 - 4*
     &    hl**4*m1**6*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 2*hl**4*mt**2*s*t2**(-1)*u2 + 4*hl**4*mt**2*s + 2*
     &    hl**4*mt**2*s**2*t2**(-1) + 2*hl**4*mt**2*t2 + 2*hl**4*mt**2*
     &    u2 + 2*hl**4*mt**4*s*t2**(-1) + 2*hl**4*mt**4 - 10*hr**4*
     &    m1**2*mt**2*s*t2**(-1) - 8*hr**4*m1**2*mt**2*t2**(-1)*u2 - 6*
     &    hr**4*m1**2*mt**2 - 4*hr**4*m1**2*mt**4*t2**(-1) - 8*hr**4*
     &    m1**2*s*t2**(-1)*u2 - 4*hr**4*m1**2*s - 4*hr**4*m1**2*s**2*
     &    t2**(-1) - 4*hr**4*m1**2*t2**(-1)*u2**2 - 2*hr**4*m1**2*t2 - 
     &    4*hr**4*m1**2*u2 + 8*hr**4*m1**4*mt**2*t2**(-1) + 8*hr**4*
     &    m1**4*s*t2**(-1) + 8*hr**4*m1**4*t2**(-1)*u2 + 4*hr**4*m1**4
     &     - 4*hr**4*m1**6*t2**(-1) + 2*hr**4*mt**2*s*t2**(-1)*u2 + 4*
     &    hr**4*mt**2*s + 2*hr**4*mt**2*s**2*t2**(-1) + 2*hr**4*mt**2*
     &    t2 + 2*hr**4*mt**2*u2 + 2*hr**4*mt**4*s*t2**(-1) + 2*hr**4*
     &    mt**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1) - 48*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*t2**(-1)
     &     - 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2 - 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**4*s**(-1) + 32*(u2-m12+mt2)**(-1)*hl**2*m1**4*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 - 16
     &    *(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2 + 24*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2 - 16*(u2-m12+mt2)**(-1)*
     &    hl**2*mt**4*s**(-1) + 16*(u2-m12+mt2)**(-1)*hl**2*mt**4*
     &    t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 + 
     &    8*(u2-m12+mt2)**(-1)*hl**2*s + 8*(u2-m12+mt2)**(-1)*hl**2
     &    *t2 + 8*(u2-m12+mt2)**(-1)*hl**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) - 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2 - 16*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**4*s**(-1) + 32*(u2-m12+mt2)**(-1)*hr**2*m1**4*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 - 16
     &    *(u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2 + 24*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*hr**2*mt**2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**4*t2**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s + 8*(u2-m12+mt2)**(-1)*hr**2*
     &    t2 + 8*(u2-m12+mt2)**(-1)*hr**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*s4**2*Pi
     & *alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1)
     &     - 16*(u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hl**2 + 16*(u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*t2**(-1) - 16*(u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1)
     &     - 8*(u2-m12+mt2)**(-1)*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*t2 - 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*u2 + 40*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s*
     &    t2**(-1) + 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*t2**(-1)
     &    *u2 - 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*s**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*t2**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2**2 + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1)*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s + 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*s**2*t2**(-1) + 16*(u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*t2**(-1)*u2**2 + 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*t2
     &     - 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 32*(u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*t2 + 32
     &    *(u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*u2 - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*t2**(-1)*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4 - 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**6*s**(-1) + 16*(u2-m12+mt2)**(-1)*hl**2*m1**6*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2*u2
     &     - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s - 8*(u2-m12+mt2)**(-1)*
     &    hl**2*mt**2*s**2*t2**(-1) - 8*(u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*t2 + 8*(u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1)*t2 - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**4*s*t2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*t2 - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 40*(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s*
     &    t2**(-1) + 32*(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2**(-1)
     &    *u2 - 8*(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*s**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*t2**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2**2 + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1)*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s + 16*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**2*t2**(-1) + 16*(u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*t2**(-1)*u2**2 + 8*(u2-m12+mt2)**(-1)*hr**2*m1**2*t2
     &     - 8*(u2-m12+mt2)**(-1)*hr**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*s**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*t2**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 32*(u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*
     &    u2 - 32*(u2-m12+mt2)**(-1)*hr**2*m1**4*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*t2**(-1)*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4 - 16*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**6*s**(-1) + 16*(u2-m12+mt2)**(-1)*hr**2*m1**6*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2*u2
     &     - 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s - 8*(u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*s**2*t2**(-1) - 8*(u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*t2 + 8*(u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*t2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**4*s*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*s**(-1) - 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s**(-1)*t2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s**(-1)*u2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*t2**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    t2**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**(-1)*t2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**(-1)*u2 + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    t2**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    t2*u2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1) - 
     &    48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*
     &    t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*t2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*u2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s**(-1) + 
     &    32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**(-1)*u2 + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1) + 
     &    16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2**(-1) + 16*(u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t2**(-1) + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*t2**(-1) + 8*hl**2*hr**2*m1**2*mt**2*
     &    t2**(-2) - 4*hl**2*hr**2*mt**2*s*t2**(-2) - 8*hl**2*hr**2*
     &    mt**2*t2**(-2)*u2 - 8*hl**2*hr**2*mt**2*t2**(-1) - 8*hl**2*
     &    hr**2*mt**4*t2**(-2) + 8*hl**4*m1**2*mt**2*s**(-1)*t2**(-1)
     &     + 12*hl**4*m1**2*mt**2*t2**(-2) + 4*hl**4*m1**2*s**(-1)*
     &    t2**(-1)*u2 - 2*hl**4*m1**2*s**(-1) + 8*hl**4*m1**2*s*
     &    t2**(-2) + 8*hl**4*m1**2*t2**(-2)*u2 + 16*hl**4*m1**2*
     &    t2**(-1) - 4*hl**4*m1**4*s**(-1)*t2**(-1) - 8*hl**4*m1**4*
     &    t2**(-2) - 4*hl**4*mt**2*s**(-1)*t2**(-1)*u2 + 2*hl**4*mt**2*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 6*hl**4*mt**2*s*t2**(-2) - 4*hl**4*mt**2*t2**(-2)
     &    *u2 - 12*hl**4*mt**2*t2**(-1) - 4*hl**4*mt**4*s**(-1)*
     &    t2**(-1) - 4*hl**4*mt**4*t2**(-2) + 2*hl**4*s**(-1)*u2 - 4*
     &    hl**4*s*t2**(-1) - 2*hl**4*t2**(-1)*u2 - 2*hl**4 + 8*hr**4*
     &    m1**2*mt**2*s**(-1)*t2**(-1) + 12*hr**4*m1**2*mt**2*t2**(-2)
     &     + 4*hr**4*m1**2*s**(-1)*t2**(-1)*u2 - 2*hr**4*m1**2*s**(-1)
     &     + 8*hr**4*m1**2*s*t2**(-2) + 8*hr**4*m1**2*t2**(-2)*u2 + 16*
     &    hr**4*m1**2*t2**(-1) - 4*hr**4*m1**4*s**(-1)*t2**(-1) - 8*
     &    hr**4*m1**4*t2**(-2) - 4*hr**4*mt**2*s**(-1)*t2**(-1)*u2 + 2*
     &    hr**4*mt**2*s**(-1) - 6*hr**4*mt**2*s*t2**(-2) - 4*hr**4*
     &    mt**2*t2**(-2)*u2 - 12*hr**4*mt**2*t2**(-1) - 4*hr**4*mt**4*
     &    s**(-1)*t2**(-1) - 4*hr**4*mt**4*t2**(-2) + 2*hr**4*s**(-1)*
     &    u2 - 4*hr**4*s*t2**(-1) - 2*hr**4*t2**(-1)*u2 - 2*hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*s4**2*Pi*alphas
     & *hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz + 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2**(-1) - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2**(-1) - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz - 4*hl**4*m1**2*
     &    t2**(-2) + 4*hl**4*mt**2*t2**(-2) + 2*hl**4*t2**(-1) - 4*
     &    hr**4*m1**2*t2**(-2) + 4*hr**4*mt**2*t2**(-2) + 2*hr**4*
     &    t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*s**(-1)*t2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1)*u2 + 40*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*s*t2**(-1) + 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*t2**(-1)*
     &    u2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 - 
     &    16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*s**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*
     &    t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*t2*u2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s**(-1)*u2**2 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s*t2**(-1)*u2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**2*s**2*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*t2**(-1)*u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*s**(-1) - 
     &    32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*
     &    t2**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s**(-1)*t2 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s**(-1)*u2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    t2**(-1)*u2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*s**(-1) + 
     &    16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2*u2 - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*t2**(-1)*u2
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1)*t2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s*t2**(-1)
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*
     &    s**(-1)*t2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*s**(-1)*u2 + 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*s*t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*t2**(-1)*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**4*s**(-1) + 
     &    16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**4*
     &    t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*t2*u2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*u2**2 + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s*t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s + 
     &    16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**2*t2**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2**(-1)*
     &    u2**2 + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 - 8
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*s**(-1) - 
     &    32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*
     &    t2**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1)*t2 + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1)*u2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    t2**(-1)*u2 + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6*s**(-1) + 
     &    16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*
     &    t2**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &    *s - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2*
     &    t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*t2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s*t2**(-1)
     &     - 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *m1**2*mt*s**(-1)*u2 + 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*m1**2*mt*s*t2**(-1) + 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1)*u2 - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt - 32*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt**3*s**(-1) + 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt**3*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**4*mt*s**(-1) - 16*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*m1**4*mt*t2**(-1) + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**(-1)*u2**2 - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*t2**(-1)*u2 + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s
     &     - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt*s**2*t2**(-1) - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*t2**(-1)*u2**2 + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2
     &     + 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt**3*s**(-1)*u2 - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s*t2**(-1) - 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3 + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**5*s**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**5
     &    *t2**(-1) - 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s**(-1)*u2 + 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*s*t2**(-1) + 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*t2**(-1)*u2 - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt - 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) + 32*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt**3*t2**(-1) + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**4*mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**4*mt*t2**(-1) + 16*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)*u2**2 - 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*
     &    t2**(-1)*u2 + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**2*t2**(-1) - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    t2**(-1)*u2**2 + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*u2 + 32*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3*s**(-1)*u2 - 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *s*t2**(-1) - 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*t2**(-1)*u2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**5*s**(-1) - 16*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**5*t2**(-1) - 16*hl**2*hr**2*
     &    m1**2*mt**2*s**(-1)*t2**(-1)*u2 - 12*hl**2*hr**2*m1**2*mt**2*
     &    s*t2**(-2) - 16*hl**2*hr**2*m1**2*mt**2*t2**(-2)*u2 - 28*
     &    hl**2*hr**2*m1**2*mt**2*t2**(-1) - 16*hl**2*hr**2*m1**2*mt**4
     &    *s**(-1)*t2**(-1) - 16*hl**2*hr**2*m1**2*mt**4*t2**(-2) + 8*
     &    hl**2*hr**2*m1**4*mt**2*s**(-1)*t2**(-1) + 8*hl**2*hr**2*
     &    m1**4*mt**2*t2**(-2) + 8*hl**2*hr**2*mt**2*s**(-1)*t2**(-1)*
     &    u2**2 + 12*hl**2*hr**2*mt**2*s*t2**(-2)*u2 + 16*hl**2*hr**2*
     &    mt**2*s*t2**(-1) + 4*hl**2*hr**2*mt**2*s**2*t2**(-2) + 8*
     &    hl**2*hr**2*mt**2*t2**(-2)*u2**2 + 28*hl**2*hr**2*mt**2*
     &    t2**(-1)*u2 + 4*hl**2*hr**2*mt**2 + 16*hl**2*hr**2*mt**4*
     &    s**(-1)*t2**(-1)*u2 + 12*hl**2*hr**2*mt**4*s*t2**(-2) + 16*
     &    hl**2*hr**2*mt**4*t2**(-2)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 28*hl**2*hr**2*mt**4*t2**(-1) + 8*hl**2*hr**2*mt**6*
     &    s**(-1)*t2**(-1) + 8*hl**2*hr**2*mt**6*t2**(-2) - 8*hl**4*
     &    m1**2*mt**2*s**(-1)*t2**(-1)*u2 - 4*hl**4*m1**2*mt**2*s**(-1)
     &     - 10*hl**4*m1**2*mt**2*s*t2**(-2) - 8*hl**4*m1**2*mt**2*
     &    t2**(-2)*u2 - 22*hl**4*m1**2*mt**2*t2**(-1) - 4*hl**4*m1**2*
     &    mt**4*s**(-1)*t2**(-1) - 4*hl**4*m1**2*mt**4*t2**(-2) - 4*
     &    hl**4*m1**2*s**(-1)*t2**(-1)*u2**2 - 2*hl**4*m1**2*s**(-1)*u2
     &     - 8*hl**4*m1**2*s*t2**(-2)*u2 - 16*hl**4*m1**2*s*t2**(-1) - 
     &    4*hl**4*m1**2*s**2*t2**(-2) - 4*hl**4*m1**2*t2**(-2)*u2**2 - 
     &    18*hl**4*m1**2*t2**(-1)*u2 - 8*hl**4*m1**2 + 8*hl**4*m1**4*
     &    mt**2*s**(-1)*t2**(-1) + 8*hl**4*m1**4*mt**2*t2**(-2) + 8*
     &    hl**4*m1**4*s**(-1)*t2**(-1)*u2 + 2*hl**4*m1**4*s**(-1) + 8*
     &    hl**4*m1**4*s*t2**(-2) + 8*hl**4*m1**4*t2**(-2)*u2 + 18*hl**4
     &    *m1**4*t2**(-1) - 4*hl**4*m1**6*s**(-1)*t2**(-1) - 4*hl**4*
     &    m1**6*t2**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 2*hl**4*mt**2*s**(-1)*u2 + 2*hl**4*mt**2*s*t2**(-2)*
     &    u2 + 8*hl**4*mt**2*s*t2**(-1) + 2*hl**4*mt**2*s**2*t2**(-2)
     &     + 4*hl**4*mt**2*t2**(-1)*u2 + 6*hl**4*mt**2 + 2*hl**4*mt**4*
     &    s**(-1) + 2*hl**4*mt**4*s*t2**(-2) + 4*hl**4*mt**4*t2**(-1)
     &     + 2*hl**4*s*t2**(-1)*u2 + 4*hl**4*s + 2*hl**4*s**2*t2**(-1)
     &     + 2*hl**4*t2 + 2*hl**4*u2 - 8*hr**4*m1**2*mt**2*s**(-1)*
     &    t2**(-1)*u2 - 4*hr**4*m1**2*mt**2*s**(-1) - 10*hr**4*m1**2*
     &    mt**2*s*t2**(-2) - 8*hr**4*m1**2*mt**2*t2**(-2)*u2 - 22*hr**4
     &    *m1**2*mt**2*t2**(-1) - 4*hr**4*m1**2*mt**4*s**(-1)*t2**(-1)
     &     - 4*hr**4*m1**2*mt**4*t2**(-2) - 4*hr**4*m1**2*s**(-1)*
     &    t2**(-1)*u2**2 - 2*hr**4*m1**2*s**(-1)*u2 - 8*hr**4*m1**2*s*
     &    t2**(-2)*u2 - 16*hr**4*m1**2*s*t2**(-1) - 4*hr**4*m1**2*s**2*
     &    t2**(-2) - 4*hr**4*m1**2*t2**(-2)*u2**2 - 18*hr**4*m1**2*
     &    t2**(-1)*u2 - 8*hr**4*m1**2 + 8*hr**4*m1**4*mt**2*s**(-1)*
     &    t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*hr**4*m1**4*mt**2*t2**(-2) + 8*hr**4*m1**4*s**(-1)
     &    *t2**(-1)*u2 + 2*hr**4*m1**4*s**(-1) + 8*hr**4*m1**4*s*
     &    t2**(-2) + 8*hr**4*m1**4*t2**(-2)*u2 + 18*hr**4*m1**4*
     &    t2**(-1) - 4*hr**4*m1**6*s**(-1)*t2**(-1) - 4*hr**4*m1**6*
     &    t2**(-2) + 2*hr**4*mt**2*s**(-1)*u2 + 2*hr**4*mt**2*s*
     &    t2**(-2)*u2 + 8*hr**4*mt**2*s*t2**(-1) + 2*hr**4*mt**2*s**2*
     &    t2**(-2) + 4*hr**4*mt**2*t2**(-1)*u2 + 6*hr**4*mt**2 + 2*
     &    hr**4*mt**4*s**(-1) + 2*hr**4*mt**4*s*t2**(-2) + 4*hr**4*
     &    mt**4*t2**(-1) + 2*hr**4*s*t2**(-1)*u2 + 4*hr**4*s + 2*hr**4*
     &    s**2*t2**(-1) + 2*hr**4*t2 + 2*hr**4*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-2)*ssp**2*pq2*Nc*Cf*s4*
     & Pi*alphas*hardfac * (  - 64*s - 64*t2 - 64*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-2)*ssp**2*pq2*Nc*Cf*
     & s4**2*Pi*alphas*hardfac * ( 64 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-2)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * ( 64*m1**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2
     &     + 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*s + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*mt**2*t2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2
     &     + 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*s + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*t2 + 8*(u2-m12+mt2)**(-1)*hl**2*m1**2 - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2 + 8*(u2-m12+mt2)**(-1)*
     &    hl**2*s + 8*(u2-m12+mt2)**(-1)*hl**2*t2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2 - 8*(u2-m12+mt2)**(-1)*
     &    hr**2*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*hr**2*s + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*t2 - 16*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*s - 8*(t2+u2-m12+mt2)**(-1)*hl**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*u2 - 16*(t2+u2-m12+mt2)**(-1)
     &    *hr**2*s - 8*(t2+u2-m12+mt2)**(-1)*hr**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*u2 + 64*lq*ssz*mz**(-2)*s + 64*
     &    lq*ssz*mz**(-2)*t2 + 64*lq*ssz*mz**(-2)*u2 + 64*rq*ssz*
     &    mz**(-2)*s + 64*rq*ssz*mz**(-2)*t2 + 64*rq*ssz*mz**(-2)*u2 + 
     &    8*hl**2 + 8*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*pq*ssp*Nc*Cf*s4**2*Pi
     & *alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 - 16*(u2-m12+mt2)**(-1)
     &    *(t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2 - 16
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2 - 8
     &    *(u2-m12+mt2)**(-1)*hl**2 - 8*(u2-m12+mt2)**(-1)*hr**2 + 
     &    8*(t2+u2-m12+mt2)**(-1)*hl**2 + 8*(t2+u2-m12+mt2)**(-1)*
     &    hr**2 - 64*lq*ssz*mz**(-2) - 64*rq*ssz*mz**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**2
     &     - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2
     &    *s*t2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*s**2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s*t2 + 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**2 - 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2
     &     - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2
     &    *s**2 - 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*s + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*t2 - 8*(u2-m12+mt2)**(-1)*
     &    hl**2*mt**2*s - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*t2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s + 8*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*s - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2*u2 - 64*lq*ssz*m1**2*
     &    mz**(-2)*s - 64*rq*ssz*m1**2*mz**(-2)*s + 8*hl**2*s + 8*hl**2
     &    *t2 + 8*hr**2*s + 8*hr**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*ssp**2*pq2*Nc*Cf*s4*
     & Pi*alphas*hardfac * ( 64*s**(-1)*t2 + 128*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,5,-1,-1)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * (  - 128*m1**2 - 64*s**(-1)*t2*u2 - 64*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s - 
     &    48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    t2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*u2 + 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**4 + 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*mt**2*s + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*u2
     &     + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**4 - 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *m1**2*mt**2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s - 48*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4 + 24
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s
     &     + 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*t2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*u2 + 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4 + 48*(u2-m12+mt2)**(-1)
     &    *hl**2*m1**2*mt**2*t2**(-1) + 32*(u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*s*t2**(-1) + 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    t2**(-1)*u2 + 48*(u2-m12+mt2)**(-1)*hl**2*m1**2 - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*t2**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1)*u2 - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*hl**2*mt**4*
     &    t2**(-1) + 48*(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2**(-1)
     &     + 32*(u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1)*u2 + 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2 - 32*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**4*t2**(-1) - 24*(u2-m12+mt2)**(-1)*hr**2*mt**2*s*
     &    t2**(-1) - 16*(u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1)*u2 - 
     &    32*(u2-m12+mt2)**(-1)*hr**2*mt**2 - 16*(u2-m12+mt2)**(-1)
     &    *hr**2*mt**4*t2**(-1) + 32*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*mt**2*s**(-1) + 24*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    s**(-1)*t2 + 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*
     &    u2 + 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) - 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2 - 
     &    16*(t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1) - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s - 24*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*t2 - 8*(t2+u2-m12+mt2)**(-1)*hl**2*u2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) - 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 8*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2
     &    *u2 - 8*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s - 24*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*t2 - 8*(t2+u2-m12+mt2)**(-1)*hr**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*s4**2*Pi
     & *alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 - 16*(u2-m12+mt2)**(-1)
     &    *(t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2 - 16
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2 - 
     &    16*(u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1) + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 + 8*(t2+u2-m12+mt2)**(-1)*
     &    hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 40*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s + 56*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    mt**2*t2 + 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*mt**2*u2 + 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4 + 48*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2
     &     + 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*s*u2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*s**2 + 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    t2**2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2
     &    *m1**2*u2**2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s - 48*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**4*t2
     &     - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**4*u2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**6 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*mt**2*s*t2 - 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s*u2 - 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**2
     &     - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2
     &    *t2*u2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2
     &    *mt**2*t2**2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*mt**4*s - 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*t2 + 40*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 56*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**2*u2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*mt**4 + 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s*u2
     &     + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*s**2 + 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*t2*u2 + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2**2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    u2**2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *m1**4*mt**2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*t2 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4*u2
     &     + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**6 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *mt**2*s*t2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*s*u2 - 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**2 - 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*t2*
     &    u2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*t2**2 - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**4*s - 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*t2 - 40*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 56*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2
     &     - 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1)*u2 - 48*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s - 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*s**2*t2**(-1) - 16*(u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*t2**(-1)*u2**2 - 32*(u2-m12+mt2)**(-1)*hl**2*m1**2*t2
     &     - 48*(u2-m12+mt2)**(-1)*hl**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*t2**(-1)*u2 + 48*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4 - 16*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**6*t2**(-1) + 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2)**(-1)*hl**2*mt**2*s + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s**2*t2**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**4*s*t2**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**4 - 40*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*mt**2*s*t2**(-1) - 32*(u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*mt**2*t2**(-1)*u2 - 56*(u2-m12+mt2)**(-1)*hr**2*m1**2
     &    *mt**2 - 16*(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*t2**(-1)
     &     - 32*(u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1)*u2 - 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s - 16*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**2*t2**(-1) - 16*(u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*t2**(-1)*u2**2 - 32*(u2-m12+mt2)**(-1)*hr**2*m1**2*t2
     &     - 48*(u2-m12+mt2)**(-1)*hr**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*t2**(-1)*u2 + 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*hr**2*m1**6*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1)*u2
     &     + 16*(u2-m12+mt2)**(-1)*hr**2*mt**2*s + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**2*t2**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*t2 + 8*(u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*u2 + 8*(u2-m12+mt2)**(-1)*hr**2*mt**4*s*
     &    t2**(-1) + 8*(u2-m12+mt2)**(-1)*hr**2*mt**4 - 48*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1)*t2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1)*u2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*s**(-1) - 40*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2*u2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 48*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2
     &     - 24*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*u2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*s**(-1) + 40*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*t2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**6*s**(-1) + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**3
     &     + 24*(t2+u2-m12+mt2)**(-1)*hl**2*s*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*t2**2 - 48*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*t2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*u2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*s**(-1) - 40*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2*u2 - 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s - 48*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 24*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2
     &     + 32*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*s**(-1) + 40*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*t2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**6*s**(-1) + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 24*(t2+u2-m12+mt2)**(-1)*hr**2*s*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*t2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*u2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4 + 24*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*u2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 - 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s - 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*t2 - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4 + 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4 + 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2**(-1)*u2 + 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2
     &     - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*t2**(-1)
     &     - 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*t2**(-1)
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2**(-1)*
     &    u2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*t2**(-1) + 48*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2**(-1)
     &     + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*t2**(-1)
     &     + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2**(-1)*
     &    u2 + 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*t2**(-1) - 24*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*t2**(-1) - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2**(-1)*u2 - 32
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*
     &    t2**(-1) - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2 + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2 + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3 + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t2**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*t2**(-1) - 16*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*m1**2*mt + 8*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s + 16*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*t2 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*u2 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2**(-1) - 16*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*t2**(-1)*u2 - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt - 
     &    16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*t2**(-1) + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**2*mt**2*s**(-1) + 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*s**(-1) - 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2 - 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*u2
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s - 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 + 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1)
     &     + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)
     &    *t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1) - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s**(-1)*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*t2*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*t2**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 8*hl**2*hr**2*
     &    m1**2*mt**2*t2**(-2) + 4*hl**2*hr**2*mt**2*s*t2**(-2) + 8*
     &    hl**2*hr**2*mt**2*t2**(-2)*u2 + 8*hl**2*hr**2*mt**2*t2**(-1)
     &     + 8*hl**2*hr**2*mt**4*t2**(-2) - 8*hl**4*m1**2*mt**2*s**(-1)
     &    *t2**(-1) - 12*hl**4*m1**2*mt**2*t2**(-2) - 4*hl**4*m1**2*
     &    s**(-1)*t2**(-1)*u2 - 6*hl**4*m1**2*s**(-1) - 8*hl**4*m1**2*s
     &    *t2**(-2) - 8*hl**4*m1**2*t2**(-2)*u2 - 16*hl**4*m1**2*
     &    t2**(-1) + 4*hl**4*m1**4*s**(-1)*t2**(-1) + 8*hl**4*m1**4*
     &    t2**(-2) + 4*hl**4*mt**2*s**(-1)*t2**(-1)*u2 + 6*hl**4*mt**2*
     &    s**(-1) + 6*hl**4*mt**2*s*t2**(-2) + 4*hl**4*mt**2*t2**(-2)*
     &    u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 12*hl**4*mt**2*t2**(-1) + 4*hl**4*mt**4*s**(-1)*
     &    t2**(-1) + 4*hl**4*mt**4*t2**(-2) + 2*hl**4*s**(-1)*t2 + 2*
     &    hl**4*s**(-1)*u2 + 4*hl**4*s*t2**(-1) + 2*hl**4*t2**(-1)*u2
     &     + 6*hl**4 - 8*hr**4*m1**2*mt**2*s**(-1)*t2**(-1) - 12*hr**4*
     &    m1**2*mt**2*t2**(-2) - 4*hr**4*m1**2*s**(-1)*t2**(-1)*u2 - 6*
     &    hr**4*m1**2*s**(-1) - 8*hr**4*m1**2*s*t2**(-2) - 8*hr**4*
     &    m1**2*t2**(-2)*u2 - 16*hr**4*m1**2*t2**(-1) + 4*hr**4*m1**4*
     &    s**(-1)*t2**(-1) + 8*hr**4*m1**4*t2**(-2) + 4*hr**4*mt**2*
     &    s**(-1)*t2**(-1)*u2 + 6*hr**4*mt**2*s**(-1) + 6*hr**4*mt**2*s
     &    *t2**(-2) + 4*hr**4*mt**2*t2**(-2)*u2 + 12*hr**4*mt**2*
     &    t2**(-1) + 4*hr**4*mt**4*s**(-1)*t2**(-1) + 4*hr**4*mt**4*
     &    t2**(-2) + 2*hr**4*s**(-1)*t2 + 2*hr**4*s**(-1)*u2 + 4*hr**4*
     &    s*t2**(-1) + 2*hr**4*t2**(-1)*u2 + 6*hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*s4**2*Pi*alphas
     & *hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*t2**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*t2**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &     + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 4*hl**4*
     &    m1**2*t2**(-2) - 4*hl**4*mt**2*t2**(-2) - 2*hl**4*t2**(-1) + 
     &    4*hr**4*m1**2*t2**(-2) - 4*hr**4*mt**2*t2**(-2) - 2*hr**4*
     &    t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 40*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*s + 56*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*t2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*u2 + 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**4 + 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s*t2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s*u2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**2 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*t2*u2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*s - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*t2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*u2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s*t2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*t2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**4*t2 + 40*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s + 56*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*t2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*u2 + 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**4 + 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*t2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s*u2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**2 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2**2 - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*s - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*t2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*u2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s*t2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2**2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4*t2 - 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2*s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*t2**(-1)*u2 - 56*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*t2**(-1)
     &     - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s*t2**(-1)
     &    *u2 - 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**2*t2**(-1) - 
     &    16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2**(-1)*
     &    u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2 - 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*t2**(-1)
     &     + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*s*t2**(-1)
     &     + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*t2**(-1)*
     &    u2 + 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*t2**(-1)*u2 + 
     &    16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2*t2**(-1) + 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 - 40*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s*t2**(-1)
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*t2**(-1)*u2 - 56*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *m1**2*mt**2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**4*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s*t2**(-1)*u2 - 48*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *m1**2*s**2*t2**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*t2**(-1)*u2**2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*t2 - 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*u2 + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*mt**2*t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*s*t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*t2**(-1)*u2 + 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**4 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6*
     &    t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &    *s + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2*
     &    t2**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s*t2**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4 + 24*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s + 40*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*t2 + 32*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*u2 + 32*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*m1**2*mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**4*mt - 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*t2 - 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*u2 - 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**2 - 40*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2*u2 - 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2**2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s - 40*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*t2 - 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*u2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**5 - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s*t2**(-1) - 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1)*u2 - 40*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt - 32*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt**3*t2**(-1) + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4
     &    *mt*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2**(-1)*u2 + 32*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**2*t2**(-1) + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*t2**(-1)*u2**2 + 24*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2
     &     + 40*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt*u2 + 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s*t2**(-1) + 32*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*t2**(-1)*u2 + 40*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &     + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt**5*t2**(-1) + 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 40*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*t2 + 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*u2 + 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt**3 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**4*mt - 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s*t2 - 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s*u2 - 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 40*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *t2*u2 - 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *t2**2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2**2 - 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s - 40*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*t2 - 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*u2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**5 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*t2**(-1) - 32*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1)*u2 - 40*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt - 32*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**3*t2**(-1) + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4
     &    *mt*t2**(-1) + 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s*t2**(-1)*u2 + 32*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    s**2*t2**(-1) + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*t2**(-1)*u2**2 + 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2
     &     + 40*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s*t2**(-1) + 32*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*t2**(-1)*u2 + 40*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &     + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt**5*t2**(-1) - 48*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2*s**(-1)*t2 - 32*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*mt**2*s**(-1)*u2 - 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*s**(-1)
     &     - 40*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**(-1)
     &    *t2*u2 - 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*t2**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s**(-1)*u2**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 48*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*t2 - 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2
     &    *u2 + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    mt**2*s**(-1) + 40*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s**(-1)*t2 + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4*s**(-1)*u2 + 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**6*s**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s**(-1)*t2*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s + 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1)*t2 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    t2**2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2**3 + 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2*u2 + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2**2 - 48*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1)
     &    *t2 - 32*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*s**(-1)*u2 - 32*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mt**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mt**4*s**(-1) - 40*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s**(-1)*t2*u2 - 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)*t2**2
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)
     &    *u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s - 48*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2 - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 + 
     &    32*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*
     &    s**(-1) + 40*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1)*t2 + 32*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s**(-1)*u2 + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**6*s**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s**(-1)*t2*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s + 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2**2*u2
     &     + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2**3
     &     + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2*u2 + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2**2 - 32*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*t2 - 32*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**2*mt*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt - 32*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**4*mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1)*t2*u2 + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*t2**2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s**(-1)*u2**2 + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*u2 + 32*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mt**3*s**(-1)*t2 + 32*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s**(-1)*u2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3 + 16*(t2+u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt**5*s**(-1) - 32*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*m1**2*mt*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt - 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**4*mt*s**(-1) + 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*t2*u2 + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**(-1)*t2**2 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*u2**2 + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *t2 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*u2 + 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2
     &    *lambda2*sqrt2**(-1)*mt**3*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**5*s**(-1) + 16*hl**2*hr**2*m1**2*mt**2*
     &    s**(-1)*t2**(-1)*u2 + 16*hl**2*hr**2*m1**2*mt**2*s**(-1) + 12
     &    *hl**2*hr**2*m1**2*mt**2*s*t2**(-2) + 16*hl**2*hr**2*m1**2*
     &    mt**2*t2**(-2)*u2 + 28*hl**2*hr**2*m1**2*mt**2*t2**(-1) + 16*
     &    hl**2*hr**2*m1**2*mt**4*s**(-1)*t2**(-1) + 16*hl**2*hr**2*
     &    m1**2*mt**4*t2**(-2) - 8*hl**2*hr**2*m1**4*mt**2*s**(-1)*
     &    t2**(-1) - 8*hl**2*hr**2*m1**4*mt**2*t2**(-2) - 8*hl**2*hr**2
     &    *mt**2*s**(-1)*t2**(-1)*u2**2 - 8*hl**2*hr**2*mt**2*s**(-1)*
     &    t2 - 16*hl**2*hr**2*mt**2*s**(-1)*u2 - 12*hl**2*hr**2*mt**2*s
     &    *t2**(-2)*u2 - 16*hl**2*hr**2*mt**2*s*t2**(-1) - 4*hl**2*
     &    hr**2*mt**2*s**2*t2**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hl**2*hr**2*mt**2*t2**(-2)*u2**2 - 28*hl**2*
     &    hr**2*mt**2*t2**(-1)*u2 - 20*hl**2*hr**2*mt**2 - 16*hl**2*
     &    hr**2*mt**4*s**(-1)*t2**(-1)*u2 - 16*hl**2*hr**2*mt**4*
     &    s**(-1) - 12*hl**2*hr**2*mt**4*s*t2**(-2) - 16*hl**2*hr**2*
     &    mt**4*t2**(-2)*u2 - 28*hl**2*hr**2*mt**4*t2**(-1) - 8*hl**2*
     &    hr**2*mt**6*s**(-1)*t2**(-1) - 8*hl**2*hr**2*mt**6*t2**(-2)
     &     + 8*hl**4*m1**2*mt**2*s**(-1)*t2**(-1)*u2 + 12*hl**4*m1**2*
     &    mt**2*s**(-1) + 10*hl**4*m1**2*mt**2*s*t2**(-2) + 8*hl**4*
     &    m1**2*mt**2*t2**(-2)*u2 + 22*hl**4*m1**2*mt**2*t2**(-1) + 4*
     &    hl**4*m1**2*mt**4*s**(-1)*t2**(-1) + 4*hl**4*m1**2*mt**4*
     &    t2**(-2) + 4*hl**4*m1**2*s**(-1)*t2**(-1)*u2**2 + 8*hl**4*
     &    m1**2*s**(-1)*t2 + 10*hl**4*m1**2*s**(-1)*u2 + 8*hl**4*m1**2*
     &    s*t2**(-2)*u2 + 16*hl**4*m1**2*s*t2**(-1) + 4*hl**4*m1**2*
     &    s**2*t2**(-2) + 4*hl**4*m1**2*t2**(-2)*u2**2 + 18*hl**4*m1**2
     &    *t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 20*hl**4*m1**2 - 8*hl**4*m1**4*mt**2*s**(-1)*
     &    t2**(-1) - 8*hl**4*m1**4*mt**2*t2**(-2) - 8*hl**4*m1**4*
     &    s**(-1)*t2**(-1)*u2 - 10*hl**4*m1**4*s**(-1) - 8*hl**4*m1**4*
     &    s*t2**(-2) - 8*hl**4*m1**4*t2**(-2)*u2 - 18*hl**4*m1**4*
     &    t2**(-1) + 4*hl**4*m1**6*s**(-1)*t2**(-1) + 4*hl**4*m1**6*
     &    t2**(-2) - 4*hl**4*mt**2*s**(-1)*t2 - 2*hl**4*mt**2*s**(-1)*
     &    u2 - 2*hl**4*mt**2*s*t2**(-2)*u2 - 8*hl**4*mt**2*s*t2**(-1)
     &     - 2*hl**4*mt**2*s**2*t2**(-2) - 4*hl**4*mt**2*t2**(-1)*u2 - 
     &    10*hl**4*mt**2 - 2*hl**4*mt**4*s**(-1) - 2*hl**4*mt**4*s*
     &    t2**(-2) - 4*hl**4*mt**4*t2**(-1) - 2*hl**4*s**(-1)*t2*u2 - 2
     &    *hl**4*s**(-1)*t2**2 - 2*hl**4*s*t2**(-1)*u2 - 6*hl**4*s - 2*
     &    hl**4*s**2*t2**(-1) - 6*hl**4*t2 - 4*hl**4*u2 + 8*hr**4*m1**2
     &    *mt**2*s**(-1)*t2**(-1)*u2 + 12*hr**4*m1**2*mt**2*s**(-1) + 
     &    10*hr**4*m1**2*mt**2*s*t2**(-2) + 8*hr**4*m1**2*mt**2*
     &    t2**(-2)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 22*hr**4*m1**2*mt**2*t2**(-1) + 4*hr**4*m1**2*mt**4*
     &    s**(-1)*t2**(-1) + 4*hr**4*m1**2*mt**4*t2**(-2) + 4*hr**4*
     &    m1**2*s**(-1)*t2**(-1)*u2**2 + 8*hr**4*m1**2*s**(-1)*t2 + 10*
     &    hr**4*m1**2*s**(-1)*u2 + 8*hr**4*m1**2*s*t2**(-2)*u2 + 16*
     &    hr**4*m1**2*s*t2**(-1) + 4*hr**4*m1**2*s**2*t2**(-2) + 4*
     &    hr**4*m1**2*t2**(-2)*u2**2 + 18*hr**4*m1**2*t2**(-1)*u2 + 20*
     &    hr**4*m1**2 - 8*hr**4*m1**4*mt**2*s**(-1)*t2**(-1) - 8*hr**4*
     &    m1**4*mt**2*t2**(-2) - 8*hr**4*m1**4*s**(-1)*t2**(-1)*u2 - 10
     &    *hr**4*m1**4*s**(-1) - 8*hr**4*m1**4*s*t2**(-2) - 8*hr**4*
     &    m1**4*t2**(-2)*u2 - 18*hr**4*m1**4*t2**(-1) + 4*hr**4*m1**6*
     &    s**(-1)*t2**(-1) + 4*hr**4*m1**6*t2**(-2) - 4*hr**4*mt**2*
     &    s**(-1)*t2 - 2*hr**4*mt**2*s**(-1)*u2 - 2*hr**4*mt**2*s*
     &    t2**(-2)*u2 - 8*hr**4*mt**2*s*t2**(-1) - 2*hr**4*mt**2*s**2*
     &    t2**(-2) - 4*hr**4*mt**2*t2**(-1)*u2 - 10*hr**4*mt**2 - 2*
     &    hr**4*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 2*hr**4*mt**4*s*t2**(-2) - 4*hr**4*mt**4*t2**(-1)
     &     - 2*hr**4*s**(-1)*t2*u2 - 2*hr**4*s**(-1)*t2**2 - 2*hr**4*s*
     &    t2**(-1)*u2 - 6*hr**4*s - 2*hr**4*s**2*t2**(-1) - 6*hr**4*t2
     &     - 4*hr**4*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*ssz**2*lq2*mz**2*s**(-1)*t2 + 64*ssz**2*lq2*mz**2
     &    *s**(-1)*u2 + 64*ssz**2*lq2*mz**4*s**(-1) - 32*ssz**2*lq2*s
     &     - 32*ssz**2*lq2*t2 - 32*ssz**2*lq2*u2 + 32*ssz**2*rq2*mz**2*
     &    s**(-1)*t2 + 64*ssz**2*rq2*mz**2*s**(-1)*u2 + 64*ssz**2*rq2*
     &    mz**4*s**(-1) - 32*ssz**2*rq2*s - 32*ssz**2*rq2*t2 - 32*
     &    ssz**2*rq2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-2)*Nc*Cf*s4**2*Pi*
     & alphas*hardfac * ( 32*ssz**2*lq2 + 32*ssz**2*rq2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssz**2*lq2*m1**2*mz**2 + 64*ssz**2*lq2*m1**2*
     &    mz**4*s**(-1) + 32*ssz**2*lq2*m1**2*s - 32*ssz**2*lq2*mz**2*
     &    s**(-1)*t2*u2 - 32*ssz**2*lq2*mz**2*u2 - 32*ssz**2*lq2*mz**4*
     &    s**(-1)*t2 - 32*ssz**2*lq2*mz**4 - 64*ssz**2*rq2*m1**2*mz**2
     &     + 64*ssz**2*rq2*m1**2*mz**4*s**(-1) + 32*ssz**2*rq2*m1**2*s
     &     - 32*ssz**2*rq2*mz**2*s**(-1)*t2*u2 - 32*ssz**2*rq2*mz**2*u2
     &     - 32*ssz**2*rq2*mz**4*s**(-1)*t2 - 32*ssz**2*rq2*mz**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 64*lq*ssz*mz**(-2)*s - 64*lq*ssz*mz**(-2)*
     &    t2 - 64*lq*ssz*mz**(-2)*u2 + 128*lq*ssz*mz**2*s**(-1) + 64*lq
     &    *ssz*s**(-1)*t2 + 128*lq*ssz*s**(-1)*u2 - 64*rq*ssz*mz**(-2)*
     &    s - 64*rq*ssz*mz**(-2)*t2 - 64*rq*ssz*mz**(-2)*u2 + 128*rq*
     &    ssz*mz**2*s**(-1) + 64*rq*ssz*s**(-1)*t2 + 128*rq*ssz*s**(-1)
     &    *u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*pq*ssp*Nc*Cf*s4**2*
     & Pi*alphas*hardfac * ( 64*lq*ssz*mz**(-2) + 64*rq*ssz*mz**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 64*lq*ssz*m1**2*mz**(-2)*s + 128*lq*ssz*m1**2
     &    *mz**2*s**(-1) - 128*lq*ssz*m1**2 - 64*lq*ssz*mz**2*s**(-1)*
     &    t2 - 64*lq*ssz*mz**2 - 64*lq*ssz*s**(-1)*t2*u2 - 64*lq*ssz*u2
     &     + 64*rq*ssz*m1**2*mz**(-2)*s + 128*rq*ssz*m1**2*mz**2*
     &    s**(-1) - 128*rq*ssz*m1**2 - 64*rq*ssz*mz**2*s**(-1)*t2 - 64*
     &    rq*ssz*mz**2 - 64*rq*ssz*s**(-1)*t2*u2 - 64*rq*ssz*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2 + 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*t2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*mz**2*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    s**(-1)*t2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*s**(-1) + 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*
     &    s**(-1)*t2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*s**(-1)*t2 - 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*s**(-1)*u2
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**4*s**(-1) - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1)*t2 - 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1)*u2
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s - 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 16*lq*hl**2*
     &    ssz*mz**2*s**(-1) + 8*lq*hl**2*ssz - 16*rq*hr**2*ssz*mz**2*
     &    s**(-1) + 8*rq*hr**2*ssz + 128*ssz**2*lq2*mz**2*s**(-1) + 32*
     &    ssz**2*lq2*s**(-1)*t2 + 64*ssz**2*lq2*s**(-1)*u2 + 128*ssz**2
     &    *rq2*mz**2*s**(-1) + 32*ssz**2*rq2*s**(-1)*t2 + 64*ssz**2*rq2
     &    *s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*s4**2*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*s - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2*t2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**4 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s*t2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2*s + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s*t2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*t2 + 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**4 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s**2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2*t2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1)*t2 + 24*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mz**4*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*
     &    s**(-1)*t2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    mz**2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*s**(-1)*t2
     &     + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 - 16
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**4*s**(-1)
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*t2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**4*
     &    s**(-1) + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    s**(-1)*t2*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mz**2*s - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*u2 + 8
     &    *(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**4*s**(-1)*t2
     &     + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**4 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**4*s**(-1)
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*
     &    s**(-1)*t2*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mz**2*s - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*
     &    t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*u2 + 8
     &    *(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**4*s**(-1)*t2
     &     + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**4 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2*u2 + 8*lq*hl**2*
     &    ssz*mz**2*s**(-1)*t2 + 8*lq*hl**2*ssz*mz**2 + 8*lq*hl**2*ssz*
     &    s + 8*lq*hl**2*ssz*t2 + 8*rq*hr**2*ssz*mz**2*s**(-1)*t2 + 8*
     &    rq*hr**2*ssz*mz**2 + 8*rq*hr**2*ssz*s + 8*rq*hr**2*ssz*t2 + 
     &    128*ssz**2*lq2*m1**2*mz**2*s**(-1) - 64*ssz**2*lq2*m1**2 - 64
     &    *ssz**2*lq2*mz**2*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssz**2*lq2*mz**2 - 32*ssz**2*lq2*s**(-1)*t2*u2
     &     - 32*ssz**2*lq2*u2 + 128*ssz**2*rq2*m1**2*mz**2*s**(-1) - 64
     &    *ssz**2*rq2*m1**2 - 64*ssz**2*rq2*mz**2*s**(-1)*t2 - 64*
     &    ssz**2*rq2*mz**2 - 32*ssz**2*rq2*s**(-1)*t2*u2 - 32*ssz**2*
     &    rq2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,11,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h1**2*lambda1**2*mh1**2 - 16*h1**2*lambda1**2*
     &    mh1**4*s**(-1) - 8*h1**2*lambda1**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,11,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2
     &     - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**4*
     &    s**(-1) - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 24
     &    *(u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s + 8*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*t2 - 16*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*mh1**4 - 8*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*s*t2 - 8*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*s**2 - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    mh1**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**4*s**(-1) + 16*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*s - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**4*s**(-1) - 32*h1**2*lambda1**2*mh1**2*
     &    s**(-1) + 16*h1**2*lambda1**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,12,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h2**2*lambda2**2*mh2**2 - 16*h2**2*lambda2**2*
     &    mh2**4*s**(-1) - 8*h2**2*lambda2**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,12,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s )
      MMcrossed5 = MMcrossed5 + ANGfin(1,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2 + 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh2**4*
     &    s**(-1) + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 24
     &    *(u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s + 8*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*t2 - 16*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mh2**4 - 8*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s*t2 - 8*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s**2 - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    mh2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(1,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**4*s**(-1) + 16*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s - 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**4*s**(-1) - 32*h2**2*lambda2**2*mh2**2*
     &    s**(-1) + 16*h2**2*lambda2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 80*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s + 
     &    16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    t2 + 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*u2 - 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**4 - 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*mt**2*s - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*u2
     &     - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**4 + 80*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *m1**2*mt**2 + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2 - 48*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4 - 24
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*t2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*u2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4 - 16*(u2-m12+mt2)**(-1)
     &    *hl**2*m1**2*mt**2*s**(-1) + 8*(u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*s**(-1)*t2 + 8*(u2-m12+mt2)**(-1)*hl**2*m1**2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2 + 8*(u2-m12+mt2)**(-1)*
     &    hl**2*mt**4*s**(-1) - 8*(u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2
     &    *u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 8*(u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2**2
     &     - 8*(u2-m12+mt2)**(-1)*hl**2*s - 8*(u2-m12+mt2)**(-1)*
     &    hl**2*t2 - 16*(u2-m12+mt2)**(-1)*hl**2*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2 + 8*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**4*s**(-1) - 8*(u2-m12+mt2)**(-1)*hr**2*mt**2*
     &    s**(-1)*t2 - 8*(u2-m12+mt2)**(-1)*hr**2*mt**2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2**2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s - 8*(u2-m12+mt2)**(-1)*hr**2*
     &    t2 - 16*(u2-m12+mt2)**(-1)*hr**2*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 
     &    16*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 8*hl**2*m1**2*s**(-1)
     &     - 8*hl**2*mt**2*s**(-1) + 8*hr**2*m1**2*s**(-1) - 8*hr**2*
     &    mt**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*s4**2*Pi
     & *alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 16*(u2-m12+mt2)**(-1)
     &    *(t2+u2-m12+mt2)**(-1)*hl**2*mt**2 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2 + 16
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 8
     &    *(u2-m12+mt2)**(-1)*hl**2 + 8*(u2-m12+mt2)**(-1)*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 64*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s - 24*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    mt**2*t2 - 80*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*mt**2*u2 - 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2
     &     - 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*s*u2 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*s**2 - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2*u2 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    u2**2 + 72*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2
     &    *m1**4*mt**2 + 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*t2 + 64*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**4*u2
     &     - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**6 + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*s*t2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*mt**2*s*u2 + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2*
     &    u2 + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*u2**2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*mt**4*s + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**4*u2
     &     + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**6
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 64*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s - 24*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**2*t2 - 80*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*mt**2*u2 - 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2
     &     - 48*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*s*u2 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**2 - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2*u2 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    u2**2 + 72*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *m1**4*mt**2 + 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*t2 + 64*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4*u2
     &     - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**6 + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*s*t2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*s*u2 + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*t2*
     &    u2 + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    mt**2*u2**2 + 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hr**2*mt**4*s + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*t2 + 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**4*u2
     &     + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**6
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*t2 + 64*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*u2 + 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 + 40*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*s**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2*u2 + 24*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2**2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s - 8*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*t2 + 16*(u2-m12+mt2)**(-1)*hl**2*m1**2*u2 - 56*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*t2 - 48*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*u2 - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4 + 24*(u2-m12+mt2)**(-1)*
     &    hl**2*m1**6*s**(-1) - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*
     &    s**(-1)*t2*u2 - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*
     &    u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2 - 8*(u2-m12+mt2)**(-1)*
     &    hl**2*mt**4*s**(-1)*t2 - 16*(u2-m12+mt2)**(-1)*hl**2*mt**4*
     &    s**(-1)*u2 - 8*(u2-m12+mt2)**(-1)*hl**2*mt**6*s**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*t2 + 64*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2 + 40*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*s**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2*u2 + 24*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2**2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s - 8*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*t2 + 16*(u2-m12+mt2)**(-1)*hr**2*m1**2*u2 - 56*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*t2 - 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*hr**2*m1**4 + 24*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**6*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2*u2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2**2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s + 8*(u2-m12+mt2)**(-1)*
     &    hr**2*mt**2*t2 - 8*(u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*
     &    t2 - 16*(u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*u2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**6*s**(-1) - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*t2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)
     &    *t2*u2 + 8*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2**2
     &     - 8*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,5,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)
     &    *u2**2 - 8*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4 - 32*hl**2*m1**2*mt**2*
     &    s**(-1) - 16*hl**2*m1**2*s**(-1)*u2 - 16*hl**2*m1**2 + 24*
     &    hl**2*m1**4*s**(-1) + 8*hl**2*mt**4*s**(-1) - 32*hr**2*m1**2*
     &    mt**2*s**(-1) - 16*hr**2*m1**2*s**(-1)*u2 - 16*hr**2*m1**2 + 
     &    24*hr**2*m1**4*s**(-1) + 8*hr**2*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*mt**2*s - 4*hl**4*m1**2*mt**2 + 2*
     &    hl**4*mt**2*s + 4*hl**4*mt**4 - 4*hr**4*m1**2*mt**2 + 2*hr**4
     &    *mt**2*s + 4*hr**4*mt**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*m1**2*mt**2*s + 4*hl**2*hr**2*m1**2
     &    *mt**2*t2 + 8*hl**2*hr**2*m1**2*mt**4 + 4*hl**2*hr**2*mt**2*s
     &    *t2 + 4*hl**2*hr**2*mt**2*s*u2 + 4*hl**2*hr**2*mt**2*s**2 - 4
     &    *hl**2*hr**2*mt**4*s - 4*hl**2*hr**2*mt**4*t2 - 8*hl**2*hr**2
     &    *mt**4*u2 - 8*hl**2*hr**2*mt**6 + 6*hl**4*m1**2*mt**2*s + 2*
     &    hl**4*m1**2*mt**2*t2 + 4*hl**4*m1**2*mt**2*u2 + 4*hl**4*m1**2
     &    *mt**4 - 4*hl**4*m1**4*mt**2 - 2*hl**4*mt**2*s*t2 - 2*hl**4*
     &    mt**2*s*u2 - 2*hl**4*mt**2*s**2 - 2*hl**4*mt**4*s - 2*hl**4*
     &    mt**4*t2 + 6*hr**4*m1**2*mt**2*s + 2*hr**4*m1**2*mt**2*t2 + 4
     &    *hr**4*m1**2*mt**2*u2 + 4*hr**4*m1**2*mt**4 - 4*hr**4*m1**4*
     &    mt**2 - 2*hr**4*mt**2*s*t2 - 2*hr**4*mt**2*s*u2 - 2*hr**4*
     &    mt**2*s**2 - 2*hr**4*mt**4*s - 2*hr**4*mt**4*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*hl**2*hr**2*m1**2*mt**2*t2**(-1) - 4*hl**2*hr**2*
     &    mt**2*s*t2**(-1) - 8*hl**2*hr**2*mt**2*t2**(-1)*u2 - 8*hl**2*
     &    hr**2*mt**4*t2**(-1) + 12*hl**4*m1**2*mt**2*t2**(-1) + 8*
     &    hl**4*m1**2*s*t2**(-1) + 8*hl**4*m1**2*t2**(-1)*u2 + 4*hl**4*
     &    m1**2 - 8*hl**4*m1**4*t2**(-1) - 6*hl**4*mt**2*s*t2**(-1) - 4
     &    *hl**4*mt**2*t2**(-1)*u2 - 4*hl**4*mt**2 - 4*hl**4*mt**4*
     &    t2**(-1) + 12*hr**4*m1**2*mt**2*t2**(-1) + 8*hr**4*m1**2*s*
     &    t2**(-1) + 8*hr**4*m1**2*t2**(-1)*u2 + 4*hr**4*m1**2 - 8*
     &    hr**4*m1**4*t2**(-1) - 6*hr**4*mt**2*s*t2**(-1) - 4*hr**4*
     &    mt**2*t2**(-1)*u2 - 4*hr**4*mt**2 - 4*hr**4*mt**4*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-1)*Nc*Cf*s4**2*Pi*alphas
     & *hardfac * (  - 4*hl**4*m1**2*t2**(-1) + 4*hl**4*mt**2*t2**(-1)
     &     - 4*hr**4*m1**2*t2**(-1) + 4*hr**4*mt**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 12*hl**2*hr**2*m1**2*mt**2*s*t2**(-1) - 16*hl**2*
     &    hr**2*m1**2*mt**2*t2**(-1)*u2 - 4*hl**2*hr**2*m1**2*mt**2 - 
     &    16*hl**2*hr**2*m1**2*mt**4*t2**(-1) + 8*hl**2*hr**2*m1**4*
     &    mt**2*t2**(-1) + 12*hl**2*hr**2*mt**2*s*t2**(-1)*u2 + 4*hl**2
     &    *hr**2*mt**2*s**2*t2**(-1) + 8*hl**2*hr**2*mt**2*t2**(-1)*
     &    u2**2 + 4*hl**2*hr**2*mt**2*u2 + 12*hl**2*hr**2*mt**4*s*
     &    t2**(-1) + 16*hl**2*hr**2*mt**4*t2**(-1)*u2 + 4*hl**2*hr**2*
     &    mt**4 + 8*hl**2*hr**2*mt**6*t2**(-1) - 10*hl**4*m1**2*mt**2*s
     &    *t2**(-1) - 8*hl**4*m1**2*mt**2*t2**(-1)*u2 - 6*hl**4*m1**2*
     &    mt**2 - 4*hl**4*m1**2*mt**4*t2**(-1) - 8*hl**4*m1**2*s*
     &    t2**(-1)*u2 - 4*hl**4*m1**2*s - 4*hl**4*m1**2*s**2*t2**(-1)
     &     - 4*hl**4*m1**2*t2**(-1)*u2**2 - 2*hl**4*m1**2*t2 - 4*hl**4*
     &    m1**2*u2 + 8*hl**4*m1**4*mt**2*t2**(-1) + 8*hl**4*m1**4*s*
     &    t2**(-1) + 8*hl**4*m1**4*t2**(-1)*u2 + 4*hl**4*m1**4 - 4*
     &    hl**4*m1**6*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-2,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 2*hl**4*mt**2*s*t2**(-1)*u2 + 4*hl**4*mt**2*s + 2*
     &    hl**4*mt**2*s**2*t2**(-1) + 2*hl**4*mt**2*t2 + 2*hl**4*mt**2*
     &    u2 + 2*hl**4*mt**4*s*t2**(-1) + 2*hl**4*mt**4 - 10*hr**4*
     &    m1**2*mt**2*s*t2**(-1) - 8*hr**4*m1**2*mt**2*t2**(-1)*u2 - 6*
     &    hr**4*m1**2*mt**2 - 4*hr**4*m1**2*mt**4*t2**(-1) - 8*hr**4*
     &    m1**2*s*t2**(-1)*u2 - 4*hr**4*m1**2*s - 4*hr**4*m1**2*s**2*
     &    t2**(-1) - 4*hr**4*m1**2*t2**(-1)*u2**2 - 2*hr**4*m1**2*t2 - 
     &    4*hr**4*m1**2*u2 + 8*hr**4*m1**4*mt**2*t2**(-1) + 8*hr**4*
     &    m1**4*s*t2**(-1) + 8*hr**4*m1**4*t2**(-1)*u2 + 4*hr**4*m1**4
     &     - 4*hr**4*m1**6*t2**(-1) + 2*hr**4*mt**2*s*t2**(-1)*u2 + 4*
     &    hr**4*mt**2*s + 2*hr**4*mt**2*s**2*t2**(-1) + 2*hr**4*mt**2*
     &    t2 + 2*hr**4*mt**2*u2 + 2*hr**4*mt**4*s*t2**(-1) + 2*hr**4*
     &    mt**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*mt**2 - 8*hl**4*m1**2*mt**2*s**(-1)
     &     - 4*hl**4*m1**2 + 4*hl**4*m1**4*s**(-1) + 6*hl**4*mt**2 + 4*
     &    hl**4*mt**4*s**(-1) + 2*hl**4*s - 8*hr**4*m1**2*mt**2*s**(-1)
     &     - 4*hr**4*m1**2 + 4*hr**4*m1**4*s**(-1) + 6*hr**4*mt**2 + 4*
     &    hr**4*mt**4*s**(-1) + 2*hr**4*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*hl**2*hr**2*m1**2*mt**2*s**(-1)*t2 + 8*hl**2*hr**2
     &    *m1**2*mt**2*s**(-1)*u2 + 4*hl**2*hr**2*m1**2*mt**2 + 16*
     &    hl**2*hr**2*m1**2*mt**4*s**(-1) - 8*hl**2*hr**2*m1**4*mt**2*
     &    s**(-1) + 4*hl**2*hr**2*mt**2*t2 - 4*hl**2*hr**2*mt**4*
     &    s**(-1)*t2 - 8*hl**2*hr**2*mt**4*s**(-1)*u2 - 12*hl**2*hr**2*
     &    mt**4 - 8*hl**2*hr**2*mt**6*s**(-1) + 6*hl**4*m1**2*mt**2*
     &    s**(-1)*t2 + 4*hl**4*m1**2*mt**2*s**(-1)*u2 + 14*hl**4*m1**2*
     &    mt**2 + 4*hl**4*m1**2*mt**4*s**(-1) + 6*hl**4*m1**2*s + 4*
     &    hl**4*m1**2*t2 + 4*hl**4*m1**2*u2 - 8*hl**4*m1**4*mt**2*
     &    s**(-1) - 4*hl**4*m1**4*s**(-1)*t2 - 4*hl**4*m1**4*s**(-1)*u2
     &     - 8*hl**4*m1**4 + 4*hl**4*m1**6*s**(-1) - 6*hl**4*mt**2*s - 
     &    6*hl**4*mt**2*t2 - 4*hl**4*mt**2*u2 - 2*hl**4*mt**4*s**(-1)*
     &    t2 - 2*hl**4*mt**4 - 2*hl**4*s*t2 - 2*hl**4*s*u2 - 2*hl**4*
     &    s**2 + 6*hr**4*m1**2*mt**2*s**(-1)*t2 + 4*hr**4*m1**2*mt**2*
     &    s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 14*hr**4*m1**2*mt**2 + 4*hr**4*m1**2*mt**4*s**(-1)
     &     + 6*hr**4*m1**2*s + 4*hr**4*m1**2*t2 + 4*hr**4*m1**2*u2 - 8*
     &    hr**4*m1**4*mt**2*s**(-1) - 4*hr**4*m1**4*s**(-1)*t2 - 4*
     &    hr**4*m1**4*s**(-1)*u2 - 8*hr**4*m1**4 + 4*hr**4*m1**6*
     &    s**(-1) - 6*hr**4*mt**2*s - 6*hr**4*mt**2*t2 - 4*hr**4*mt**2*
     &    u2 - 2*hr**4*mt**4*s**(-1)*t2 - 2*hr**4*mt**4 - 2*hr**4*s*t2
     &     - 2*hr**4*s*u2 - 2*hr**4*s**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 48*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s + 
     &    32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    u2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**4 - 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2
     &    *mt**2*s - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*mt**2*u2 - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4 + 48*(u2-m12+mt2)**(-1)
     &    *(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s + 
     &    32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    u2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**4 - 24*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2
     &    *mt**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*u2 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**4 - 48
     &    *(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1)*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*t2**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*s*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**4*t2**(-1) - 48*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1)*u2 + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*t2**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 16*(u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1)*
     &    u2 + 16*(u2-m12+mt2)**(-1)*hr**2*mt**4*t2**(-1) + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) + 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1) + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) + 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*s4**2*Pi
     & *alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 16*(u2-m12+mt2)**(-1)
     &    *(t2+u2-m12+mt2)**(-1)*hl**2*mt**2 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2 + 16
     &    *(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 
     &    16*(u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 40*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    mt**2*u2 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*mt**4 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s*u2 - 16*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**2
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**2*u2**2 + 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hl**2*m1**4*mt**2 + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**4*u2
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    m1**6 + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*
     &    mt**2*s*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s - 
     &    40*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**2*s - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*mt**2*u2 - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4 - 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s*u2
     &     - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*s**2 - 16*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*u2**2 + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s + 
     &    32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**4*
     &    u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**6 + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s*u2 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**2
     &     + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**4
     &    *s + 40*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s*t2**(-1) + 
     &    32*(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**2*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*t2**(-1)*u2**2 - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**6*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s*
     &    t2**(-1)*u2 - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**2*
     &    t2**(-1) - 8*(u2-m12+mt2)**(-1)*hl**2*mt**4*s*t2**(-1) + 40
     &    *(u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**4*t2**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**2*t2**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*t2**(-1)*u2**2 - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s*t2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**6*t2**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s*t2**(-1)*u2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 8*(u2-m12+mt2)**(-1)*hr**2*mt**4*s*
     &    t2**(-1) - 24*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*t2 - 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*
     &    s**(-1)*u2 - 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**4*
     &    s**(-1) + 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*u2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*mt**2*s**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**6*s**(-1) - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1)*t2 - 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*
     &    mt**4*s**(-1) + 16*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s + 16
     &    *(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*u2 + 32*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*mt**2*s**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**6*s**(-1) - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4 - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*u2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4 - 24*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4 - 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*t2**(-1)*u2 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**4*t2**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*mt**4*t2**(-1) - 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mt**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s
     &    *t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2**(-1)*u2 + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4
     &    *t2**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s
     &    *t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4
     &    *t2**(-1) + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt - 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t2**(-1) + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    t2**(-1)*u2 + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*t2**(-1) + 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt - 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2**(-1) + 16*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *t2**(-1) + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4*s**(-1) + 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**4*s**(-1) + 32*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*s**(-1) + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt + 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt + 8*hl**2*hr**2*m1**2*mt**2*t2**(-2)
     &     - 4*hl**2*hr**2*mt**2*s*t2**(-2) - 8*hl**2*hr**2*mt**2*
     &    t2**(-2)*u2 - 8*hl**2*hr**2*mt**2*t2**(-1) - 8*hl**2*hr**2*
     &    mt**4*t2**(-2) + 8*hl**4*m1**2*mt**2*s**(-1)*t2**(-1) + 12*
     &    hl**4*m1**2*mt**2*t2**(-2) + 4*hl**4*m1**2*s**(-1)*t2**(-1)*
     &    u2 - 2*hl**4*m1**2*s**(-1) + 8*hl**4*m1**2*s*t2**(-2) + 8*
     &    hl**4*m1**2*t2**(-2)*u2 + 16*hl**4*m1**2*t2**(-1) - 4*hl**4*
     &    m1**4*s**(-1)*t2**(-1) - 8*hl**4*m1**4*t2**(-2) - 4*hl**4*
     &    mt**2*s**(-1)*t2**(-1)*u2 + 2*hl**4*mt**2*s**(-1) - 6*hl**4*
     &    mt**2*s*t2**(-2) - 4*hl**4*mt**2*t2**(-2)*u2 - 12*hl**4*mt**2
     &    *t2**(-1) - 4*hl**4*mt**4*s**(-1)*t2**(-1) - 4*hl**4*mt**4*
     &    t2**(-2) - 4*hl**4*s*t2**(-1) - 2*hl**4*t2**(-1)*u2 - 4*hl**4
     &     + 8*hr**4*m1**2*mt**2*s**(-1)*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 12*hr**4*m1**2*mt**2*t2**(-2) + 4*hr**4*m1**2*
     &    s**(-1)*t2**(-1)*u2 - 2*hr**4*m1**2*s**(-1) + 8*hr**4*m1**2*s
     &    *t2**(-2) + 8*hr**4*m1**2*t2**(-2)*u2 + 16*hr**4*m1**2*
     &    t2**(-1) - 4*hr**4*m1**4*s**(-1)*t2**(-1) - 8*hr**4*m1**4*
     &    t2**(-2) - 4*hr**4*mt**2*s**(-1)*t2**(-1)*u2 + 2*hr**4*mt**2*
     &    s**(-1) - 6*hr**4*mt**2*s*t2**(-2) - 4*hr**4*mt**2*t2**(-2)*
     &    u2 - 12*hr**4*mt**2*t2**(-1) - 4*hr**4*mt**4*s**(-1)*t2**(-1)
     &     - 4*hr**4*mt**4*t2**(-2) - 4*hr**4*s*t2**(-1) - 2*hr**4*
     &    t2**(-1)*u2 - 4*hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*s4**2*Pi*alphas
     & *hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*t2**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*t2**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*t2**(-1) - 4*hl**4*m1**2*t2**(-2) + 4*hl**4*mt**2*
     &    t2**(-2) + 2*hl**4*t2**(-1) - 4*hr**4*m1**2*t2**(-2) + 4*
     &    hr**4*mt**2*t2**(-2) + 2*hr**4*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 40*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*s - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*u2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s*u2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*u2**2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*s + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**6 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s - 40*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*s - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*u2 - 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**4 - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s**2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2**2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*mt**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*u2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4*s + 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**2*mt**2*s*t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*t2**(-1)*u2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*t2**(-1)
     &     + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s*t2**(-1)
     &    *u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**2
     &    *t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2**(-1)*u2**2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*mt**2*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4*s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4*t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**6*t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s*t2**(-1)*u2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**2*t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*mt**4*s*t2**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*m1**2*mt**2*s*t2**(-1) + 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2**(-1)*
     &    u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**4*
     &    t2**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*
     &    t2**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**2
     &    *t2**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2**(-1)*u2**2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*mt**2*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*s*t2**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*t2**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**6*t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*s*t2**(-1)*u2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**2*t2**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*mt**4*s*t2**(-1) - 24*(u2-m12+mt2+mh12)**(-1)
     &    *(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s - 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*u2 - 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**4*mt + 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*u2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**2 + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2**2 + 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s + 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*u2 + 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**5 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s*t2**(-1) + 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1)*u2 + 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt**3*t2**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4
     &    *mt*t2**(-1) - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s*t2**(-1)*u2 - 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*t2**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    t2**(-1)*u2**2 - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s*t2**(-1) - 32*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *t2**(-1)*u2 - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt**5*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s - 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*u2 - 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt**3 + 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**4*mt + 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s*u2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**2 + 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s + 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*u2 + 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**5 + 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*t2**(-1) + 32*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*t2**(-1)*u2 + 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt**3*t2**(-1) - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4
     &    *mt*t2**(-1) - 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s*t2**(-1)*u2 - 8*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*t2**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2**(-1)*u2**2 - 24*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*s*t2**(-1) - 32*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *t2**(-1)*u2 - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt**5*t2**(-1) - 24*(t2+u2-m12+mt2+mz2)**(-1)
     &    *lq*hl**2*ssz*m1**2*mt**2*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*s**(-1)
     &    *u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**4*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    u2 + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2
     &    *s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4
     &    *s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**6*s**(-1) - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    s**(-1)*t2 - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**4*s**(-1)
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 + 32*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*s**(-1)
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s**(-1)
     &    *t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1)*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**6*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &    *s + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 + 8
     &    *(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*t2
     &     - 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s**(-1)*t2 - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*u2 + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**2*mt - 32*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt**3*s**(-1) + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**4*mt*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s - 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*t2 - 24*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mt*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s**(-1)*u2 - 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3 + 16*(t2+u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt**5*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*t2 - 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*m1**2*mt*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt - 32*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**4*mt*s**(-1) - 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*t2 - 24*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt*u2 + 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s**(-1)*t2 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**5*s**(-1) - 16*hl**2*hr**2*m1**2*mt**2*
     &    s**(-1)*t2**(-1)*u2 - 12*hl**2*hr**2*m1**2*mt**2*s*t2**(-2)
     &     - 16*hl**2*hr**2*m1**2*mt**2*t2**(-2)*u2 - 28*hl**2*hr**2*
     &    m1**2*mt**2*t2**(-1) - 16*hl**2*hr**2*m1**2*mt**4*s**(-1)*
     &    t2**(-1) - 16*hl**2*hr**2*m1**2*mt**4*t2**(-2) + 8*hl**2*
     &    hr**2*m1**4*mt**2*s**(-1)*t2**(-1) + 8*hl**2*hr**2*m1**4*
     &    mt**2*t2**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*hl**2*hr**2*mt**2*s**(-1)*t2**(-1)*u2**2 + 4*hl**2
     &    *hr**2*mt**2*s**(-1)*u2 + 12*hl**2*hr**2*mt**2*s*t2**(-2)*u2
     &     + 16*hl**2*hr**2*mt**2*s*t2**(-1) + 4*hl**2*hr**2*mt**2*s**2
     &    *t2**(-2) + 8*hl**2*hr**2*mt**2*t2**(-2)*u2**2 + 28*hl**2*
     &    hr**2*mt**2*t2**(-1)*u2 + 8*hl**2*hr**2*mt**2 + 16*hl**2*
     &    hr**2*mt**4*s**(-1)*t2**(-1)*u2 + 12*hl**2*hr**2*mt**4*s*
     &    t2**(-2) + 16*hl**2*hr**2*mt**4*t2**(-2)*u2 + 28*hl**2*hr**2*
     &    mt**4*t2**(-1) + 8*hl**2*hr**2*mt**6*s**(-1)*t2**(-1) + 8*
     &    hl**2*hr**2*mt**6*t2**(-2) - 8*hl**4*m1**2*mt**2*s**(-1)*
     &    t2**(-1)*u2 - 4*hl**4*m1**2*mt**2*s**(-1) - 10*hl**4*m1**2*
     &    mt**2*s*t2**(-2) - 8*hl**4*m1**2*mt**2*t2**(-2)*u2 - 22*hl**4
     &    *m1**2*mt**2*t2**(-1) - 4*hl**4*m1**2*mt**4*s**(-1)*t2**(-1)
     &     - 4*hl**4*m1**2*mt**4*t2**(-2) - 4*hl**4*m1**2*s**(-1)*
     &    t2**(-1)*u2**2 - 4*hl**4*m1**2*s**(-1)*u2 - 8*hl**4*m1**2*s*
     &    t2**(-2)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*hl**4*m1**2*s*t2**(-1) - 4*hl**4*m1**2*s**2*
     &    t2**(-2) - 4*hl**4*m1**2*t2**(-2)*u2**2 - 18*hl**4*m1**2*
     &    t2**(-1)*u2 - 10*hl**4*m1**2 + 8*hl**4*m1**4*mt**2*s**(-1)*
     &    t2**(-1) + 8*hl**4*m1**4*mt**2*t2**(-2) + 8*hl**4*m1**4*
     &    s**(-1)*t2**(-1)*u2 + 2*hl**4*m1**4*s**(-1) + 8*hl**4*m1**4*s
     &    *t2**(-2) + 8*hl**4*m1**4*t2**(-2)*u2 + 18*hl**4*m1**4*
     &    t2**(-1) - 4*hl**4*m1**6*s**(-1)*t2**(-1) - 4*hl**4*m1**6*
     &    t2**(-2) + 2*hl**4*mt**2*s**(-1)*u2 + 2*hl**4*mt**2*s*
     &    t2**(-2)*u2 + 8*hl**4*mt**2*s*t2**(-1) + 2*hl**4*mt**2*s**2*
     &    t2**(-2) + 4*hl**4*mt**2*t2**(-1)*u2 + 6*hl**4*mt**2 + 2*
     &    hl**4*mt**4*s**(-1) + 2*hl**4*mt**4*s*t2**(-2) + 4*hl**4*
     &    mt**4*t2**(-1) + 2*hl**4*s*t2**(-1)*u2 + 4*hl**4*s + 2*hl**4*
     &    s**2*t2**(-1) + 2*hl**4*t2 + 2*hl**4*u2 - 8*hr**4*m1**2*mt**2
     &    *s**(-1)*t2**(-1)*u2 - 4*hr**4*m1**2*mt**2*s**(-1) - 10*hr**4
     &    *m1**2*mt**2*s*t2**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hr**4*m1**2*mt**2*t2**(-2)*u2 - 22*hr**4*m1**2*
     &    mt**2*t2**(-1) - 4*hr**4*m1**2*mt**4*s**(-1)*t2**(-1) - 4*
     &    hr**4*m1**2*mt**4*t2**(-2) - 4*hr**4*m1**2*s**(-1)*t2**(-1)*
     &    u2**2 - 4*hr**4*m1**2*s**(-1)*u2 - 8*hr**4*m1**2*s*t2**(-2)*
     &    u2 - 16*hr**4*m1**2*s*t2**(-1) - 4*hr**4*m1**2*s**2*t2**(-2)
     &     - 4*hr**4*m1**2*t2**(-2)*u2**2 - 18*hr**4*m1**2*t2**(-1)*u2
     &     - 10*hr**4*m1**2 + 8*hr**4*m1**4*mt**2*s**(-1)*t2**(-1) + 8*
     &    hr**4*m1**4*mt**2*t2**(-2) + 8*hr**4*m1**4*s**(-1)*t2**(-1)*
     &    u2 + 2*hr**4*m1**4*s**(-1) + 8*hr**4*m1**4*s*t2**(-2) + 8*
     &    hr**4*m1**4*t2**(-2)*u2 + 18*hr**4*m1**4*t2**(-1) - 4*hr**4*
     &    m1**6*s**(-1)*t2**(-1) - 4*hr**4*m1**6*t2**(-2) + 2*hr**4*
     &    mt**2*s**(-1)*u2 + 2*hr**4*mt**2*s*t2**(-2)*u2 + 8*hr**4*
     &    mt**2*s*t2**(-1) + 2*hr**4*mt**2*s**2*t2**(-2) + 4*hr**4*
     &    mt**2*t2**(-1)*u2 + 6*hr**4*mt**2 + 2*hr**4*mt**4*s**(-1) + 2
     &    *hr**4*mt**4*s*t2**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,8,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*hr**4*mt**4*t2**(-1) + 2*hr**4*s*t2**(-1)*u2 + 4*
     &    hr**4*s + 2*hr**4*s**2*t2**(-1) + 2*hr**4*t2 + 2*hr**4*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 80*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*t2 + 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 - 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2 - 24*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 + 80*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2 - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2 - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*t2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*mz**2*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s**(-1)*t2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    s**(-1)*u2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    u2**2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1) + 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*s**(-1)
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)*t2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s**(-1) - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2 - 8
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1) - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1)*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    t2*u2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    u2**2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**(-1)
     &    *t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    mz**2*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s**(-1)*u2 + 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mz**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*t2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s**(-1)*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*mz**2*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*mt**2*s**(-1)*u2 + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 + 8*lq*hl**2*
     &    ssz*m1**2*s**(-1) - 8*lq*hl**2*ssz*mt**2*s**(-1) + 8*rq*hr**2
     &    *ssz*m1**2*s**(-1) - 8*rq*hr**2*ssz*mt**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*s4**2*Pi*
     & alphas*hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz + 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*mz**2
     &     - 64*(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*mt**2*s - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*t2 - 80
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*u2 - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*u2 - 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s*t2 - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*u2**2 + 72*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*mz**2 + 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*s + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*t2 + 64*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**6 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*t2*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**4*mz**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**4*t2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**6 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*mz**2
     &     - 64*(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*m1**2*mt**2*s - 24*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2 - 80
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*u2 - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**2*s - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*u2 - 16
     &    *(u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s*t2 - 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s*u2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2**2 + 72*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*mz**2 + 48*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4*t2 + 64*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**6 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4*mz**2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4*t2 + 16*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**6 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*mz**2*s**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2*s**(-1)*t2 + 64*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*s**(-1)*u2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 + 40*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*s**(-1) + 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*s**(-1)*
     &    u2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2 + 8
     &    *(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**(-1)*t2*u2
     &     + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s**(-1)*
     &    u2**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 - 56*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mz**2*s**(-1)
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s**(-1)*t2 - 48*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    s**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4
     &     + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*s**(-1)*
     &    u2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2 - 8
     &    *(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2*u2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*
     &    u2**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*mz**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1)*t2 - 
     &    16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s**(-1)*u2 - 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**6*s**(-1) + 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*mz**2*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*s**(-1)*t2 + 64*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1)*u2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 + 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**2*mt**4*s**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**2*s**(-1)*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)*t2*u2 + 
     &    24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)*u2**2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 + 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 - 56*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mz**2*s**(-1)
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*s**(-1)*t2
     &     )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    s**(-1)*u2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4
     &     + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*
     &    u2 + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2 - 8
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2*u2
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*
     &    u2**2 + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*mz**2*s**(-1) - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*t2 - 
     &    16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*u2 - 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**6*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*mz**2*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2
     &    *mt**2*s**(-1)*t2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2*s**(-1)*u2 + 24*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*mt**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq
     &    *hl**2*ssz*m1**2*mz**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*t2 + 32*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*u2 - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4
     &    *mz**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**4 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    mz**2*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*mz**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**2*s**(-1)*t2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**(-1)*u2**2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2
     &    *t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**4*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**4 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*mz**2*s**(-1) - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq
     &    *hr**2*ssz*m1**2*mt**2*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1)
     &    *u2 + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mz**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s
     &     + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 + 32
     &    *(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mz**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    mz**2*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*mz**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**2*s**(-1)*t2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s**(-1)*u2**2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*t2 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s**(-1)*u2
     &     - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4 - 32*lq*
     &    hl**2*ssz*m1**2*mt**2*s**(-1) - 16*lq*hl**2*ssz*m1**2*s**(-1)
     &    *u2 - 16*lq*hl**2*ssz*m1**2 + 24*lq*hl**2*ssz*m1**4*s**(-1)
     &     + 8*lq*hl**2*ssz*mt**4*s**(-1) - 32*rq*hr**2*ssz*m1**2*mt**2
     &    *s**(-1) - 16*rq*hr**2*ssz*m1**2*s**(-1)*u2 - 16*rq*hr**2*ssz
     &    *m1**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*rq*hr**2*ssz*m1**4*s**(-1) + 8*rq*hr**2*ssz*mt**4
     &    *s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt - 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3 + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*mh1**2 - 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s - 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*t2 - 48*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*u2 - 48*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt**3 + 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**4*mt + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2*u2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*t2 + 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*u2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2*u2 + 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2**2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*mh1**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s + 8*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*t2 + 48*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*u2 + 24*(u2-m12+mt2+mh12)**(-1)*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**5 + 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s**(-1)*u2 + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt + 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4
     &    *mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1)*u2**2 - 16*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2
     &     - 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt**3*s**(-1)*u2 - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3 - 16*(u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**5*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*mh1**2*s**(-1) + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*u2 + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**2*mt - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *mh1**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2 - 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr
     &    *h1*lambda1*sqrt2**(-1)*mt*s**(-1)*t2*u2 - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*u2**2 - 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s - 8*(t2+u2-m12+mt2+mh12)**(-1)*hl
     &    *hr*h1*lambda1*sqrt2**(-1)*mt*t2 - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2 - 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*mh1**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*s**(-1)*t2 - 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s**(-1)
     &     + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1)*u2 + 16*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt**3*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt - 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 + 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*mh2**2 - 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s - 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*t2 - 48*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*u2 - 48*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt**3 + 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**4*mt + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2*u2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s*t2 + 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s*u2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *t2*u2 + 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2**2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*mh2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s + 8*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*t2 + 48*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*u2 + 24*(u2-m12+mt2+mh22)**(-1)*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**5 + 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s**(-1)*u2 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt + 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt**3*s**(-1) - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4
     &    *mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**(-1)*u2**2 - 16*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s - 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*u2
     &     - 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt**3*s**(-1)*u2 - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3 - 16*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**5*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*mh2**2*s**(-1) + 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*u2 + 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*m1**2*mt - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *mh2**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2 - 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr
     &    *h2*lambda2*sqrt2**(-1)*mt*s**(-1)*t2*u2 - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**(-1)*u2**2 - 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s - 8*(t2+u2-m12+mt2+mh22)**(-1)*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*mt*t2 - 16*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2 - 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*mh2**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s**(-1)*t2 - 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3 )
      MMcrossed5 = MMcrossed5 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s**(-1)
     &     + 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)*u2 + 16*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt**3*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-2,0)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * ( 64*m1**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 8*(u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1) - 
     &    8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 + 8*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 + 8*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2
     &     - 16*(t2+u2-m12+mt2)**(-1)*hr**2 + 8*hl**2*s**(-1) + 8*
     &    hr**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 40*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*mt**2 - 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*u2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*m1**4 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*u2
     &     + 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hl**2*mt**4
     &     - 40*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2*mt**2 - 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)
     &    *hr**2*m1**2*u2 + 32*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4 + 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*u2 + 8*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**4 + 32
     &    *(u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 - 24*
     &    (u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 8*(u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2 - 
     &    8*(u2-m12+mt2)**(-1)*hl**2*mt**4*s**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 - 24*
     &    (u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2 - 8*
     &    (u2-m12+mt2)**(-1)*hr**2*mt**4*s**(-1) - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * ( 8*(t2+u2-m12+mt2)**(-1)*hl**2*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*u2 - 16*(t2+u2-m12+mt2)**(-1)
     &    *hr**2*m1**2*s**(-1)*u2 + 16*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    m1**2 + 8*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*u2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s + 8*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*t2 + 8*(t2+u2-m12+mt2)**(-1)*hr**2*u2 - 64*lq*ssz*
     &    m1**2*mz**(-2)*s**(-1)*u2 - 64*rq*ssz*m1**2*mz**(-2)*s**(-1)*
     &    u2 - 24*hl**2*m1**2*s**(-1) + 8*hl**2*mt**2*s**(-1) - 8*hl**2
     &    *s**(-1)*t2 - 8*hl**2*s**(-1)*u2 - 24*hr**2*m1**2*s**(-1) + 8
     &    *hr**2*mt**2*s**(-1) - 8*hr**2*s**(-1)*t2 - 8*hr**2*s**(-1)*
     &    u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*ssp**2*pq2*Nc*Cf*s4*Pi
     & *alphas*hardfac * (  - 64*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,0,-1,0)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * ( 64 - 64*m1**2*s**(-1) + 64*s**(-1)*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(5,7,-2,1)*ssp**2*pq2*Nc*Cf*s4*Pi
     & *alphas*hardfac * ( 64*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,7,-2,1)*ssp**2*pq2*Nc*Cf*Pi*
     & alphas*hardfac * (  - 64 + 64*m1**2*s**(-1) - 64*s**(-1)*t2 - 64
     &    *s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(5,7,-1,1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1) - 64*lq*ssz*
     &    mz**(-2)*s**(-1) - 64*rq*ssz*mz**(-2)*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(5,7,-1,1)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)
     &     + 16*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 - 16*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**(-1) + 16*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    s**(-1)*t2 + 16*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2 + 16
     &    *(t2+u2-m12+mt2)**(-1)*hr**2 - 64*lq*ssz*m1**2*mz**(-2)*
     &    s**(-1) + 64*lq*ssz*mz**(-2)*s**(-1)*t2 + 64*lq*ssz*mz**(-2)*
     &    s**(-1)*u2 + 64*lq*ssz*mz**(-2) - 64*rq*ssz*m1**2*mz**(-2)*
     &    s**(-1) + 64*rq*ssz*mz**(-2)*s**(-1)*t2 + 64*rq*ssz*mz**(-2)*
     &    s**(-1)*u2 + 64*rq*ssz*mz**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 2*hl**4*s**(-1) + 2*hr**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*mt**2*s**(-1) + 2*hl**4*m1**2*
     &    s**(-1) - 2*hl**4*s**(-1)*t2 - 2*hl**4*s**(-1)*u2 - 2*hl**4
     &     + 2*hr**4*m1**2*s**(-1) - 2*hr**4*s**(-1)*t2 - 2*hr**4*
     &    s**(-1)*u2 - 2*hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-1)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * (  - 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)
     &     + 16*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 - 16*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**(-1) + 16*(t2+u2-m12+mt2)**(-1)*hr**2*
     &    s**(-1)*t2 + 16*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2 + 16
     &    *(t2+u2-m12+mt2)**(-1)*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(7,8,1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz - 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1) + 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1) + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-2,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 2*hl**4*m1**2*s**(-1) + 2*hl**4*mt**2*s**(-1) + 4
     &    *hl**4 - 2*hr**4*m1**2*s**(-1) + 2*hr**4*mt**2*s**(-1) + 4*
     &    hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*m1**2*mt**2*s**(-1) + 4*hl**2*hr**2
     &    *mt**4*s**(-1) - 2*hl**4*m1**2*mt**2*s**(-1) + 2*hl**4*m1**2*
     &    s**(-1)*t2 + 2*hl**4*m1**2*s**(-1)*u2 + 2*hl**4*m1**2 + 2*
     &    hl**4*m1**4*s**(-1) - 2*hl**4*mt**2*s**(-1)*t2 - 2*hl**4*
     &    mt**2*s**(-1)*u2 - 2*hl**4*mt**2 - 4*hl**4*s - 4*hl**4*t2 - 4
     &    *hl**4*u2 - 2*hr**4*m1**2*mt**2*s**(-1) + 2*hr**4*m1**2*
     &    s**(-1)*t2 + 2*hr**4*m1**2*s**(-1)*u2 + 2*hr**4*m1**2 + 2*
     &    hr**4*m1**4*s**(-1) - 2*hr**4*mt**2*s**(-1)*t2 - 2*hr**4*
     &    mt**2*s**(-1)*u2 - 2*hr**4*mt**2 - 4*hr**4*s - 4*hr**4*t2 - 4
     &    *hr**4*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 16*(t2+u2-m12+mt2)**(-1)*hl**2*m1**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2 + 16*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*m1**2*s**(-1) - 16*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2*
     &    s**(-1) - 8*(t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * ( 32*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*m1**2*t2 - 8*(u2-m12+mt2)**(-1)*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*t2 + 32*
     &    (u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*m1**2*t2
     &     - 8*(u2-m12+mt2)**(-1)*(t2+u2-m12+mt2)**(-1)*hr**2*mt**2
     &    *t2 - 32*(u2-m12+mt2)**(-1)*hl**2*m1**2 + 8*
     &    (u2-m12+mt2)**(-1)*hl**2*mt**2 - 32*(u2-m12+mt2)**(-1)*
     &    hr**2*m1**2 + 8*(u2-m12+mt2)**(-1)*hr**2*mt**2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*mt**2*s**(-1) - 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*t2 - 24*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*m1**4*s**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*pq*ssp*Nc*Cf*Pi*alphas
     & *hardfac * ( 8*(t2+u2-m12+mt2)**(-1)*hl**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hl**2*s + 24*(t2+u2-m12+mt2)**(-1)*
     &    hl**2*t2 + 24*(t2+u2-m12+mt2)**(-1)*hl**2*u2 - 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*mt**2*s**(-1) - 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*t2 - 24*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*m1**4*s**(-1) + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*mt**2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2)**(-1)*hr**2*s + 24*(t2+u2-m12+mt2)**(-1)*
     &    hr**2*t2 + 24*(t2+u2-m12+mt2)**(-1)*hr**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**(-1) - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)
     &    *t2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*u2
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1) - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 2*hl**4*s**(-1)
     &     + 2*hr**4*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*t2 + 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*t2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*t2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2 - 
     &    32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 - 24*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*t2 + 24*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt - 
     &    24*(u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**2*mt**2*s**(-1) - 24*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*t2 - 24*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*s**(-1) + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**(-1)*t2 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*t2*u2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*t2**2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s**(-1)
     &     - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)
     &    *t2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *s**(-1)*t2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1)*t2**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s + 24*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 24*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 16*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1) + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s**(-1)*t2 + 24*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*u2 + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt )
      MMcrossed5 = MMcrossed5 + ANGfin(8,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s**(-1) - 16*(t2+u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s**(-1) + 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**(-1)*t2 + 24*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*u2 + 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &     + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s**(-1) - 4*hl**2*hr**2*mt**2*s**(-1) + 2*
     &    hl**4*m1**2*s**(-1) - 2*hl**4*s**(-1)*t2 - 2*hl**4 + 2*hr**4*
     &    m1**2*s**(-1) - 2*hr**4*s**(-1)*t2 - 2*hr**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,4,-2,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*hl**2*hr**2*m1**2*mt**2*s - 4*hl**2*hr**2*m1**2*
     &    mt**2*t2 - 8*hl**2*hr**2*m1**2*mt**4 + 4*hl**2*hr**2*mt**4*s
     &     + 4*hl**2*hr**2*mt**4*t2 + 8*hl**2*hr**2*mt**4*u2 + 8*hl**2*
     &    hr**2*mt**6 - 2*hl**4*m1**2*mt**2*s + 2*hl**4*m1**2*mt**2*t2
     &     - 4*hl**4*m1**2*mt**4 + 4*hl**4*m1**4*mt**2 - 2*hl**4*mt**4*
     &    s - 2*hl**4*mt**4*t2 - 4*hl**4*mt**4*u2 - 2*hr**4*m1**2*mt**2
     &    *s + 2*hr**4*m1**2*mt**2*t2 - 4*hr**4*m1**2*mt**4 + 4*hr**4*
     &    m1**4*mt**2 - 2*hr**4*mt**4*s - 2*hr**4*mt**4*t2 - 4*hr**4*
     &    mt**4*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,4,-2,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*m1**2*mt**2*s**(-1)*t2 - 8*hl**2*
     &    hr**2*m1**2*mt**2*s**(-1)*u2 - 4*hl**2*hr**2*m1**2*mt**2 - 16
     &    *hl**2*hr**2*m1**2*mt**4*s**(-1) + 8*hl**2*hr**2*m1**4*mt**2*
     &    s**(-1) + 4*hl**2*hr**2*mt**2*s + 4*hl**2*hr**2*mt**2*u2 + 4*
     &    hl**2*hr**2*mt**4*s**(-1)*t2 + 8*hl**2*hr**2*mt**4*s**(-1)*u2
     &     + 12*hl**2*hr**2*mt**4 + 8*hl**2*hr**2*mt**6*s**(-1) + 2*
     &    hl**4*m1**2*mt**2*s**(-1)*t2 + 4*hl**4*m1**2*mt**2*s**(-1)*u2
     &     - 6*hl**4*m1**2*mt**2 - 4*hl**4*m1**2*mt**4*s**(-1) - 2*
     &    hl**4*m1**2*s + 8*hl**4*m1**4*mt**2*s**(-1) + 4*hl**4*m1**4
     &     - 4*hl**4*m1**6*s**(-1) - 2*hl**4*mt**2*u2 - 2*hl**4*mt**4*
     &    s**(-1)*t2 - 4*hl**4*mt**4*s**(-1)*u2 - 2*hl**4*mt**4 + 2*
     &    hr**4*m1**2*mt**2*s**(-1)*t2 + 4*hr**4*m1**2*mt**2*s**(-1)*u2
     &     - 6*hr**4*m1**2*mt**2 - 4*hr**4*m1**2*mt**4*s**(-1) - 2*
     &    hr**4*m1**2*s + 8*hr**4*m1**4*mt**2*s**(-1) + 4*hr**4*m1**4
     &     - 4*hr**4*m1**6*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(8,4,-2,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 2*hr**4*mt**2*u2 - 2*hr**4*mt**4*s**(-1)*t2 - 4*
     &    hr**4*mt**4*s**(-1)*u2 - 2*hr**4*mt**4 )
      MMcrossed5 = MMcrossed5 + ANGfin(8,4,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*hl**2*hr**2*mt**2*s**(-1) - 2*hl**4*m1**2*s**(-1)
     &     - 2*hr**4*m1**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-2,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*ssz**2*lq2*mz**2*s**(-1) - 32*ssz**2*rq2*mz**2
     &    *s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*ssz**2*lq2*m1**2*mz**2*s**(-1) + 32*ssz**2*lq2
     &    *m1**2*s**(-1)*u2 + 32*ssz**2*lq2*mz**2*s**(-1)*t2 + 32*
     &    ssz**2*lq2*mz**2 - 32*ssz**2*rq2*m1**2*mz**2*s**(-1) + 32*
     &    ssz**2*rq2*m1**2*s**(-1)*u2 + 32*ssz**2*rq2*mz**2*s**(-1)*t2
     &     + 32*ssz**2*rq2*mz**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * (  - 64*lq*ssz*s**(-1) - 64*rq*ssz*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 64*lq*ssz*m1**2*mz**(-2)*s**(-1)*u2 - 64*lq*
     &    ssz*m1**2*s**(-1) + 64*lq*ssz*s**(-1)*t2 + 64*lq*ssz + 64*rq*
     &    ssz*m1**2*mz**(-2)*s**(-1)*u2 - 64*rq*ssz*m1**2*s**(-1) + 64*
     &    rq*ssz*s**(-1)*t2 + 64*rq*ssz )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*
     &    s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*u2
     &     + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1) - 8
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1) - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1) - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*s**(-1) + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*u2 - 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1) + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz + 8*
     &    lq*hl**2*ssz*s**(-1) + 8*rq*hr**2*ssz*s**(-1) - 32*ssz**2*lq2
     &    *s**(-1) - 32*ssz**2*rq2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 40*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2 - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*u2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4 - 40*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2 - 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2 + 8*(u2-m12+mt2+mz2)**(-1)*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*(t2+u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**4 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*u2 - 24*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*mt**2*s**(-1)*u2 - 8*(u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*mt**4*s**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mt**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mz**2*s**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*u2 - 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*mz**2*s**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*s**(-1)*u2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**4*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s**(-1)*u2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*mt**2 - 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mz**2*s**(-1)*t2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mz**2*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mz**2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    t2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**(-1)*
     &    u2**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*s**(-1)
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1)
     &    *u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*t2 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**(-1)*u2
     &     - 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 - 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1)*t2 + 
     &    8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s**(-1)*u2
     &     - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    t2*u2 + 8*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*
     &    u2**2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2 + 8*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 - 24*lq*hl**2*
     &    ssz*m1**2*s**(-1) + 8*lq*hl**2*ssz*mt**2*s**(-1) - 8*lq*hl**2
     &    *ssz*s**(-1)*t2 - 8*lq*hl**2*ssz*s**(-1)*u2 - 24*rq*hr**2*ssz
     &    *m1**2*s**(-1) + 8*rq*hr**2*ssz*mt**2*s**(-1) - 8*rq*hr**2*
     &    ssz*s**(-1)*t2 - 8*rq*hr**2*ssz*s**(-1)*u2 - 32*ssz**2*lq2*
     &    m1**2*s**(-1) + 32*ssz**2*lq2*s**(-1)*t2 + 32*ssz**2*lq2 - 32
     &    *ssz**2*rq2*m1**2*s**(-1) + 32*ssz**2*rq2*s**(-1)*t2 + 32*
     &    ssz**2*rq2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-2,1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*ssz**2*lq2*s**(-1) + 32*ssz**2*rq2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*ssz**2*lq2*m1**2*s**(-1) - 32*ssz**2*lq2*s**(-1)*
     &    t2 - 32*ssz**2*lq2*s**(-1)*u2 - 32*ssz**2*lq2 + 32*ssz**2*rq2
     &    *m1**2*s**(-1) - 32*ssz**2*rq2*s**(-1)*t2 - 32*ssz**2*rq2*
     &    s**(-1)*u2 - 32*ssz**2*rq2 )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-1,1)*pq*ssp*Nc*Cf*s4*Pi*
     & alphas*hardfac * ( 64*lq*ssz*mz**(-2)*s**(-1) + 64*rq*ssz*
     &    mz**(-2)*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-1,1)*pq*ssp*Nc*Cf*Pi*
     & alphas*hardfac * ( 64*lq*ssz*m1**2*mz**(-2)*s**(-1) - 64*lq*ssz*
     &    mz**(-2)*s**(-1)*t2 - 64*lq*ssz*mz**(-2)*s**(-1)*u2 - 64*lq*
     &    ssz*mz**(-2) + 64*rq*ssz*m1**2*mz**(-2)*s**(-1) - 64*rq*ssz*
     &    mz**(-2)*s**(-1)*t2 - 64*rq*ssz*mz**(-2)*s**(-1)*u2 - 64*rq*
     &    ssz*mz**(-2) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-1,1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1) - 16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(10,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s**(-1) + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*t2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    s**(-1)*u2 + 16*(t2+u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz - 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s**(-1) + 
     &    16*(t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*t2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**(-1)*u2 + 16*
     &    (t2+u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz )
      MMcrossed5 = MMcrossed5 + ANGfin(11,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h1**2*lambda1**2*mh1**2*s**(-1) - 8*h1**2*
     &    lambda1**2*s**(-1)*u2 - 16*h1**2*lambda1**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(11,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2*
     &    s**(-1) - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    s**(-1)*u2 - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2 - 
     &    24*(u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt + 24*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2 + 24*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*u2 + 24*
     &    (u2-m12+mt2+mh12)**(-1)*(t2+u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**3 + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*s**(-1) - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2*s**(-1) - 16*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(11,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s**(-1) - 16*(t2+u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*t2 + 8*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s**(-1)*u2 + 8*
     &    (t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &     + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1) + 8*h1**2*
     &    lambda1**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(11,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*h1**2*lambda1**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(11,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    s**(-1) + 16*(t2+u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(12,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h2**2*lambda2**2*mh2**2*s**(-1) - 8*h2**2*
     &    lambda2**2*s**(-1)*u2 - 16*h2**2*lambda2**2 )
      MMcrossed5 = MMcrossed5 + ANGfin(12,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2*s**(-1) + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2
     &    *s**(-1)*u2 + 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2 - 
     &    24*(u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt + 24*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2 + 24*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*u2 + 24*
     &    (u2-m12+mt2+mh22)**(-1)*(t2+u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*s**(-1) - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*s**(-1) - 16*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)*u2 )
      MMcrossed5 = MMcrossed5 + ANGfin(12,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s**(-1) - 16*(t2+u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s**(-1) - 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**(-1)*t2 + 8*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*u2 + 8*
     &    (t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &     + 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1) + 8*h2**2*
     &    lambda2**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(12,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*h2**2*lambda2**2*s**(-1) )
      MMcrossed5 = MMcrossed5 + ANGfin(12,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**(-1)
     &     + 16*(t2+u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**(-1) )

c               the phase space except for 1/s**2 
      HH_QGH = MMcrossed5 / ( 16.D0 * pi**2 )**2 / 2.D0*s4/(s4+m1**2)

c               the averaging factors
      HH_QGH = HH_QGH /4.D0 /Nc/(Nc**2-1.D0)

c               the prefactor for the scaling functions 
      HH_QGH = HH_QGH * (m1+m2)**2/4.D0 

      end

      































