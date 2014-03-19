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

c --------------------------------------------------------------------
      real*8 function HH_GBOS(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Nc,Cf,alphas
     &          ,m1,m2,mt,gamt
     &          ,hl,hr 
     &          ,s,s4,t2,t2t,s4t
     &          ,hardfac,theta_s3,theta_s4
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
      gamt = massin(25)

c               real kinematics built in
      t2t = t2 + m2**2 - mt**2
      s4t = s4 + m1**2 - mt**2

c               os subtraction 
      theta_s4 = 0.D0
      if ((m1.lt.mt).and.(s.gt.(m2+mt)**2)) theta_s4 = 1.D0 
      if (theta_s4.eq.1.D0) then 
         if ( s4t .gt. 0.D0 ) then
            s4t =  sqrt( s4t**2 + mt**2*gamt**2 )
         else 
            s4t = -sqrt( s4t**2 + mt**2*gamt**2 )
         end if 
      end if 

      hl      = C(4)
      hr      = C(5)

c               set gs=1 
      alphas = 1.D0/(4.D0*Pi) 

      hardfac = 1.D0

c               the s3 regularization
      theta_s3 = 0.D0

c               the angular functions 
      call ANGULAR_ARRAY_HH_GB(massin,theta_s3,ANGfin)

c               form output
      MMqgos =
     &  + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*alphas*hardfac * ( 4*hl**2
     &    *hr**2*m1**2*mt**2*s**(-1) - 4*hl**2*hr**2*m1**2*mt**2*s*
     &    t2t**(-2) - 4*hl**2*hr**2*mt**2*s**(-1)*t2t - 4*hl**2*hr**2*
     &    mt**2*s*t2t**(-1) - 8*hl**2*hr**2*mt**2 - 4*hl**2*hr**2*mt**4
     &    *s**(-1) - 4*hl**2*hr**2*mt**4*s*t2t**(-2) - 8*hl**2*hr**2*
     &    mt**4*t2t**(-1) + 6*hl**4*m1**2*mt**2*s**(-1) + 2*hl**4*m1**2
     &    *mt**2*s*t2t**(-2) + 8*hl**4*m1**2*mt**2*t2t**(-1) + 12*hl**4
     &    *m1**2*mt**4*s**(-1)*t2t**(-1) + 8*hl**4*m1**2*mt**4*
     &    t2t**(-2) + 2*hl**4*m1**2*s**(-1)*t2t + 2*hl**4*m1**2*s*
     &    t2t**(-1) + 4*hl**4*m1**2 - 12*hl**4*m1**4*mt**2*s**(-1)*
     &    t2t**(-1) - 4*hl**4*m1**4*mt**2*t2t**(-2) - 4*hl**4*m1**4*
     &    s**(-1) - 4*hl**4*m1**4*t2t**(-1) + 4*hl**4*m1**6*s**(-1)*
     &    t2t**(-1) - 2*hl**4*mt**4*s**(-1) + 2*hl**4*mt**4*s*t2t**(-2)
     &     - 4*hl**4*mt**6*s**(-1)*t2t**(-1) - 4*hl**4*mt**6*t2t**(-2)
     &     + 6*hr**4*m1**2*mt**2*s**(-1) + 2*hr**4*m1**2*mt**2*s*
     &    t2t**(-2) )
      MMqgos = MMqgos + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*alphas*
     & hardfac * ( 8*hr**4*m1**2*mt**2*t2t**(-1) + 12*hr**4*m1**2*mt**4
     &    *s**(-1)*t2t**(-1) + 8*hr**4*m1**2*mt**4*t2t**(-2) + 2*hr**4*
     &    m1**2*s**(-1)*t2t + 2*hr**4*m1**2*s*t2t**(-1) + 4*hr**4*m1**2
     &     - 12*hr**4*m1**4*mt**2*s**(-1)*t2t**(-1) - 4*hr**4*m1**4*
     &    mt**2*t2t**(-2) - 4*hr**4*m1**4*s**(-1) - 4*hr**4*m1**4*
     &    t2t**(-1) + 4*hr**4*m1**6*s**(-1)*t2t**(-1) - 2*hr**4*mt**4*
     &    s**(-1) + 2*hr**4*mt**4*s*t2t**(-2) - 4*hr**4*mt**6*s**(-1)*
     &    t2t**(-1) - 4*hr**4*mt**6*t2t**(-2) )
      MMqgos = MMqgos + ANGfin(5,0,1,0)*Nc*Cf*s4t**(-2)*Pi*alphas*
     & hardfac * (  - 8*hl**2*hr**2*m1**2*mt**2*s**(-1)*t2t**(-1) + 4*
     &    hl**2*hr**2*mt**2*s**(-1) + 4*hl**2*hr**2*mt**2*t2t**(-1) + 8
     &    *hl**2*hr**2*mt**4*s**(-1)*t2t**(-1) + 8*hl**2*hr**2*mt**4*
     &    t2t**(-2) + 4*hl**4*m1**2*mt**2*s**(-1)*t2t**(-1) - 2*hl**4*
     &    mt**2*s**(-1) - 2*hl**4*mt**2*t2t**(-1) - 4*hl**4*mt**4*
     &    s**(-1)*t2t**(-1) - 4*hl**4*mt**4*t2t**(-2) + 4*hr**4*m1**2*
     &    mt**2*s**(-1)*t2t**(-1) - 2*hr**4*mt**2*s**(-1) - 2*hr**4*
     &    mt**2*t2t**(-1) - 4*hr**4*mt**4*s**(-1)*t2t**(-1) - 4*hr**4*
     &    mt**4*t2t**(-2) )
      MMqgos = MMqgos + ANGfin(7,0,1,0)*Nc*Cf*s4t**(-2)*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*m1**2*mt**2*s**(-1)*t2t**(-1) + 4*
     &    hl**2*hr**2*m1**2*mt**2*t2t**(-2) + 4*hl**2*hr**2*mt**2*
     &    s**(-1) + 4*hl**2*hr**2*mt**2*t2t**(-1) + 4*hl**2*hr**2*mt**4
     &    *s**(-1)*t2t**(-1) + 4*hl**2*hr**2*mt**4*t2t**(-2) + 2*hl**4*
     &    m1**2*mt**2*s**(-1)*t2t**(-1) - 2*hl**4*m1**2*mt**2*t2t**(-2)
     &     - 2*hl**4*mt**2*s**(-1) - 2*hl**4*mt**2*t2t**(-1) - 2*hl**4*
     &    mt**4*s**(-1)*t2t**(-1) - 2*hl**4*mt**4*t2t**(-2) + 2*hr**4*
     &    m1**2*mt**2*s**(-1)*t2t**(-1) - 2*hr**4*m1**2*mt**2*t2t**(-2)
     &     - 2*hr**4*mt**2*s**(-1) - 2*hr**4*mt**2*t2t**(-1) - 2*hr**4*
     &    mt**4*s**(-1)*t2t**(-1) - 2*hr**4*mt**4*t2t**(-2) )

c               the phase space except for 1/s**2 
      HH_GBOS = MMqgos / ( 16.D0 * pi**2 )**2 / 2.D0*s4/(s4+m1**2)

c               the averaging factors
      HH_GBOS = HH_GBOS /4.D0 /Nc/(Nc**2-1.D0)

c               the prefactor for the scaling functions 
      HH_GBOS = HH_GBOS * (m1+m2)**2/4.D0 

      end

c --------------------------------------------------------------------
      real*8 function HH_GBH(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Nc,Cf,sqrt2,alphas
     &          ,m1,m2,mt,mz,mh1,mh2,m12,mh12,mh22,gamt
     &          ,ssp,ssz,hl,hr 
     &          ,h1,h2,lambda1,lambda2
     &          ,lq,rq,pq
     &          ,lq2,rq2,pq2
     &          ,s,s4,t2,u2,t2t,s4t
     &          ,hardfac,logall,theta_s3,theta_s4
     &          ,dyfact,dyfacu
     &          ,MMcrossed4,ANGfin(0:12,0:12,-2:2,-2:2)

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
      gamt = massin(25)

c               real kinematics built in
      u2  = s4 - s - t2 - m1**2 + m1**2 
      t2t = t2 + m2**2 - mt**2
      s4t = s4 + m1**2 - mt**2

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

c               os subtraction 
      theta_s4 = 0.D0
      if ((m1.lt.mt).and.(s.gt.(m2+mt)**2)) theta_s4 = 1.D0 
      if (theta_s4.eq.1.D0) then 
         if ( s4t .gt. 0.D0 ) then
            s4t =  sqrt( s4t**2 + mt**2*gamt**2 )
         else 
            s4t = -sqrt( s4t**2 + mt**2*gamt**2 )
         end if
      end if 

c               mass squares more convenient
      m12  = m1**2
      mh12 = mh1**2
      mh22 = mh2**2

c               logall = logqf - log(s4/m1^2) + log(1+m1^2/s4)
      logall = log(massin(13)**2/m1**2) 
     &        - log(s4/m1**2) + log(1.D0+m1**2/s4)

c               set gs=1 
      alphas = 1.D0/(4.D0*Pi) 

      hardfac = 1.D0

c               the s3 regularization
      theta_s3 = 0.D0

c               the angular functions 
      call ANGULAR_ARRAY_HH_GB(massin,theta_s3,ANGfin)

c               form output
      MMcrossed4 =
     &  + Nc*Cf*s4t**(-1)*Pi**2*alphas*hardfac * ( 4*(s+u2)**(-2)*
     &    (1+m12/s4)*hl**4*s*t2*u2*t2t**(-1) + 4*(s+u2)**(-2)*
     &    (1+m12/s4)*hl**4*t2*u2**2*t2t**(-1) + 4*(s+u2)**(-2)*
     &    (1+m12/s4)*hr**4*s*t2*u2*t2t**(-1) + 4*(s+u2)**(-2)*
     &    (1+m12/s4)*hr**4*t2*u2**2*t2t**(-1) - 16*(s+u2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**2 - 16*(s+u2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*u2 - 16*(s+u2)**(-1)*(1+m12/s4)*pq*
     &    hr**2*ssp*s**(-1)*u2**2 - 16*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2
     &    *ssp*u2 - 16*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*
     &    u2**2 - 16*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*u2 - 16*
     &    (s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*u2**2 - 16*
     &    (s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*u2 + 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hl**4*m1**2*u2*t2t**(-1) - 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hl**4*mt**2*u2*t2t**(-1) - 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hl**4*u2 + 4*(s+u2)**(-1)*(1+m12/s4)*hr**4*m1**2*
     &    u2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*s4t**(-1)*Pi**2*alphas*hardfac
     &  * (  - 4*(s+u2)**(-1)*(1+m12/s4)*hr**4*mt**2*u2*t2t**(-1) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*u2 + 16*(1+m12/s4)*pq*hl**2*ssp
     &    *s**(-1)*u2 + 16*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2 - 16*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*ssz*mz**2*s**(-1)*
     &    t2**(-1) + 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*ssz*
     &    s**(-1) - 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*mz**2*
     &    s**(-1)*t2**(-1) + 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*
     &    ssz*s**(-1) + 16*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*lq*hl**2*ssz
     &    *mz**2*s**(-2)*t2**(-1) + 16*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*
     &    rq*hr**2*ssz*mz**2*s**(-2)*t2**(-1) - 16*dyfacu(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz - 16*dyfacu(mz)*(1+m12/s4)*rq*hr**2*
     &    ssz )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 16*(s+u2)**(-3)*(1+m12/s4)*hl**2*hr**2*mt**2*s*t2**2*
     &    t2t**(-2) + 8*(s+u2)**(-3)*(1+m12/s4)*hl**4*m1**2*s*t2**2*
     &    t2t**(-2) - 8*(s+u2)**(-3)*(1+m12/s4)*hl**4*t2**3*u2*
     &    t2t**(-2) + 8*(s+u2)**(-3)*(1+m12/s4)*hr**4*m1**2*s*t2**2*
     &    t2t**(-2) - 8*(s+u2)**(-3)*(1+m12/s4)*hr**4*t2**3*u2*
     &    t2t**(-2) - 64*(s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2*
     &    t2t**(-1) + 64*(s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*
     &    t2**2*u2*t2t**(-1) - 64*(s+u2)**(-2)*(1+m12/s4)*pq*hr**2*ssp*
     &    m1**2*t2*t2t**(-1) + 64*(s+u2)**(-2)*(1+m12/s4)*pq*hr**2*ssp*
     &    s**(-1)*t2**2*u2*t2t**(-1) - 16*(s+u2)**(-2)*(1+m12/s4)*hl**2
     &    *hr**2*mt**2*s*t2*t2t**(-2) + 8*(s+u2)**(-2)*(1+m12/s4)*hl**4
     &    *m1**2*s*t2*t2t**(-2) - 8*(s+u2)**(-2)*(1+m12/s4)*hl**4*t2**2
     &    *u2*t2t**(-2) + 8*(s+u2)**(-2)*(1+m12/s4)*hr**4*m1**2*s*t2*
     &    t2t**(-2) - 8*(s+u2)**(-2)*(1+m12/s4)*hr**4*t2**2*u2*
     &    t2t**(-2) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 64*(s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2t**(-1) + 
     &    64*(s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2*u2*
     &    t2t**(-1) - 64*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) + 64*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*
     &    t2*u2*t2t**(-1) - 8*(s+u2)**(-1)*(1+m12/s4)*hl**2*hr**2*mt**2
     &    *s*t2t**(-2) + 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*m1**2*s*
     &    t2t**(-2) - 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*t2*u2*t2t**(-2)
     &     + 4*(s+u2)**(-1)*(1+m12/s4)*hr**4*m1**2*s*t2t**(-2) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*t2*u2*t2t**(-2) + 256*
     &    (s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*m1**2*s**(-1) - 256*
     &    (s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*s**(-2)*t2*u2 + 128*(s+u2)
     &    *(1+m12/s4)*ssp**2*pq2*m1**2*s**(-1)*t2**(-2) - 128*(s+u2)*
     &    (1+m12/s4)*ssp**2*pq2*s**(-2)*t2**(-1)*u2 - 32*(1+m12/s4)*pq*
     &    hl**2*ssp*m1**2*t2**(-1)*t2t**(-1) + 32*(1+m12/s4)*pq*hl**2*
     &    ssp*s**(-1)*u2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 32*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2**(-1)*t2t**(-1) + 32*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2*t2t**(-1) + 256*(1+m12/s4)
     &    *ssp**2*pq2*m1**2*s**(-1)*t2**(-1) - 256*(1+m12/s4)*ssp**2*
     &    pq2*s**(-2)*u2 - 64*dyfacu(mz)*(s+u2)**(-2)*(1+m12/s4)*lq*
     &    hl**2*ssz*m1**2*t2*t2t**(-1) - 64*dyfacu(mz)*(s+u2)**(-2)*
     &    (1+m12/s4)*lq*hl**2*ssz*t2**2*t2t**(-1) - 64*dyfacu(mz)*
     &    (s+u2)**(-2)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2*t2t**(-1) - 64*
     &    dyfacu(mz)*(s+u2)**(-2)*(1+m12/s4)*rq*hr**2*ssz*t2**2*
     &    t2t**(-1) + 256*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*pq*lq*ssz*
     &    ssp*m1**2*s**(-1) + 256*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*pq
     &    *lq*ssz*ssp*s**(-1)*t2 + 256*dyfacu(mz)*(s+u2)**(-1)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1) + 256*dyfacu(mz)*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2 - 64*dyfacu(
     &    mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2t**(-1) + 64
     &    *dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*
     &    t2**2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 64*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*t2*
     &    t2t**(-1) - 64*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*rq*hr**2*
     &    ssz*m1**2*t2t**(-1) + 64*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*
     &    rq*hr**2*ssz*s**(-1)*t2**2*t2t**(-1) - 64*dyfacu(mz)*
     &    (s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*t2*t2t**(-1) + 128*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-2) - 256*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-2) + 128*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-1)*t2**(-1) + 128*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*rq*ssz
     &    *ssp*m1**2*s**(-1)*t2**(-2) - 256*dyfacu(mz)*(s+u2)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-2) + 128*dyfacu(mz)*(s+u2)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2**(-1) + 32*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2t**(-1) + 32*dyfacu(
     &    mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2t**(-1) - 128*
     &    dyfacu(mz)*(s+u2)**2*(1+m12/s4)*pq*lq*ssz*ssp*s**(-2)*
     &    t2**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 128*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*
     &    t2**(-1) + 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*
     &    s**(-1)*t2**(-1) - 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-2)*t2 + 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*s**(-1)
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-1) - 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*t2
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1) - 32*
     &    dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2**(-1)*t2t**(-1)
     &     + 64*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2*t2t**(-1)
     &     - 32*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*t2t**(-1) - 32*
     &    dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2**(-1)*t2t**(-1)
     &     + 64*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2*t2t**(-1)
     &     - 32*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*t2t**(-1) + 128*
     &    dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*lq2*m1**2*
     &    s**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &    128*dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*lq2*s**(-1)*
     &    t2 + 128*dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*rq2*
     &    m1**2*s**(-1) + 128*dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*
     &    ssz**2*rq2*s**(-1)*t2 + 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*
     &    ssz**2*lq2*m1**2*s**(-1)*t2**(-2) - 128*dyfacu(mz)**2*(s+u2)*
     &    (1+m12/s4)*ssz**2*lq2*s**(-2) + 64*dyfacu(mz)**2*(s+u2)*
     &    (1+m12/s4)*ssz**2*lq2*s**(-1)*t2**(-1) + 64*dyfacu(mz)**2*
     &    (s+u2)*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*t2**(-2) - 128*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-2) + 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-1)*t2**(-1)
     &     - 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*lq2*s**(-2)*
     &    t2**(-1) - 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*rq2*
     &    s**(-2)*t2**(-1) + 128*dyfacu(mz)**2*(1+m12/s4)*ssz**2*lq2*
     &    m1**2*s**(-1)*t2**(-1) - 128*dyfacu(mz)**2*(1+m12/s4)*ssz**2*
     &    lq2*s**(-2)*t2 )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &    128*dyfacu(mz)**2*(1+m12/s4)*ssz**2*lq2*s**(-1) + 128*dyfacu(
     &    mz)**2*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*t2**(-1) - 128*
     &    dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-2)*t2 + 128*dyfacu(
     &    mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-1) + 64*dyfacu(mh1)*
     &    (s+u2)**(-2)*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2*
     &    t2t**(-1) + 64*dyfacu(mh1)*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*t2t**(-1) + 32*dyfacu(mh1)*(1+m12/s4)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*t2t**(-1) - 32*
     &    dyfacu(mh1)**2*(s+u2)**(-1)*(1+m12/s4)*h1**2*lambda1**2*
     &    s**(-1) - 16*dyfacu(mh1)**2*(s+u2)*(1+m12/s4)*h1**2*
     &    lambda1**2*s**(-1)*t2**(-2) - 32*dyfacu(mh1)**2*(1+m12/s4)*
     &    h1**2*lambda1**2*s**(-1)*t2**(-1) - 64*dyfacu(mh1)*dyfacu(mh2
     &    )*(s+u2)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*s**(-1) - 32*
     &    dyfacu(mh1)*dyfacu(mh2)*(s+u2)*(1+m12/s4)*h1*h2*lambda1*
     &    lambda2*s**(-1)*t2**(-2) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac*logall * ( 
     &     - 64*dyfacu(mh1)*dyfacu(mh2)*(1+m12/s4)*h1*h2*lambda1*
     &    lambda2*s**(-1)*t2**(-1) + 64*dyfacu(mh2)*(s+u2)**(-2)*
     &    (1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2*t2t**(-1) + 64*
     &    dyfacu(mh2)*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2t**(-1) + 32*dyfacu(mh2)*(1+m12/s4)*hl*hr*h2
     &    *lambda2*sqrt2**(-1)*mt*t2**(-1)*t2t**(-1) - 32*dyfacu(mh2)**
     &    2*(s+u2)**(-1)*(1+m12/s4)*h2**2*lambda2**2*s**(-1) - 16*
     &    dyfacu(mh2)**2*(s+u2)*(1+m12/s4)*h2**2*lambda2**2*s**(-1)*
     &    t2**(-2) - 32*dyfacu(mh2)**2*(1+m12/s4)*h2**2*lambda2**2*
     &    s**(-1)*t2**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 16*
     &    (s+u2)**(-3)*(1+m12/s4)*hl**2*hr**2*mt**2*s*t2**2*t2t**(-2)
     &     + 8*(s+u2)**(-3)*(1+m12/s4)*hl**4*m1**2*s*t2**2*t2t**(-2) - 
     &    8*(s+u2)**(-3)*(1+m12/s4)*hl**4*t2**3*u2*t2t**(-2) + 8*
     &    (s+u2)**(-3)*(1+m12/s4)*hr**4*m1**2*s*t2**2*t2t**(-2) - 8*
     &    (s+u2)**(-3)*(1+m12/s4)*hr**4*t2**3*u2*t2t**(-2) - 64*
     &    (s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2*t2t**(-1) + 64*
     &    (s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2**2*u2*
     &    t2t**(-1) - 64*(s+u2)**(-2)*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2*
     &    t2t**(-1) + 64*(s+u2)**(-2)*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*
     &    t2**2*u2*t2t**(-1) - 16*(s+u2)**(-2)*(1+m12/s4)*hl**2*hr**2*
     &    mt**2*s*t2*t2t**(-2) + 8*(s+u2)**(-2)*(1+m12/s4)*hl**4*m1**2*
     &    s*t2*t2t**(-2) - 8*(s+u2)**(-2)*(1+m12/s4)*hl**4*t2**2*u2*
     &    t2t**(-2) + 8*(s+u2)**(-2)*(1+m12/s4)*hr**4*m1**2*s*t2*
     &    t2t**(-2) - 8*(s+u2)**(-2)*(1+m12/s4)*hr**4*t2**2*u2*
     &    t2t**(-2) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 64*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2t**(-1) + 64*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2*u2*t2t**(-1)
     &     - 16*(s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**2*
     &    t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*u2*
     &    t2t**(-1) - 64*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) + 64*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*
     &    t2*u2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*
     &    s**(-1)*u2**2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*pq*hr**2
     &    *ssp*u2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    s**(-1)*u2**2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2
     &    *ssz*u2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*
     &    s**(-1)*u2**2*t2t**(-1) - 16*(s+u2)**(-1)*(1+m12/s4)*rq*hr**2
     &    *ssz*u2*t2t**(-1) - 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*m1**2*u2*
     &    t2t**(-2) + 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*mt**2*u2*
     &    t2t**(-2) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*t2*u2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*u2*t2t**(-1) - 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hr**4*m1**2*u2*t2t**(-2) + 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hr**4*mt**2*u2*t2t**(-2) - 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hr**4*t2*u2*t2t**(-2) + 4*(s+u2)**(-1)*(1+m12/s4)*
     &    hr**4*u2*t2t**(-1) + 256*(s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*
     &    m1**2*s**(-1) - 256*(s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*
     &    s**(-2)*t2*u2 + 128*(s+u2)*(1+m12/s4)*ssp**2*pq2*m1**2*
     &    s**(-1)*t2**(-2) - 128*(s+u2)*(1+m12/s4)*ssp**2*pq2*s**(-2)*
     &    t2**(-1)*u2 - 128*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*mz**(-2)*
     &    t2**(-1) - 128*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*mz**(-2)*
     &    t2**(-1) + 16*(1+m12/s4)*pq*hl**2*ssp*m1**2*s**(-1)*t2**(-1)*
     &    u2*t2t**(-1) - 16*(1+m12/s4)*pq*hl**2*ssp*mt**2*s**(-1)*
     &    t2**(-1)*u2*t2t**(-1) - 16*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*
     &    t2**(-1)*u2 )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * ( 32*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2*t2t**(-1) + 16*(1+m12/s4)*
     &    pq*hr**2*ssp*m1**2*s**(-1)*t2**(-1)*u2*t2t**(-1) - 16*
     &    (1+m12/s4)*pq*hr**2*ssp*mt**2*s**(-1)*t2**(-1)*u2*t2t**(-1)
     &     - 16*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*t2**(-1)*u2 + 32*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2*t2t**(-1) - 128*(1+m12/s4)
     &    *ssp**2*pq2*m1**2*s**(-1)*t2**(-2)*u2 + 256*(1+m12/s4)*ssp**2
     &    *pq2*m1**2*s**(-1)*t2**(-1) - 128*(1+m12/s4)*ssp**2*pq2*m1**2
     &    *t2**(-2) + 128*(1+m12/s4)*ssp**2*pq2*s**(-2)*t2**(-1)*u2**2
     &     - 256*(1+m12/s4)*ssp**2*pq2*s**(-2)*u2 + 128*(1+m12/s4)*
     &    ssp**2*pq2*s**(-1)*t2**(-1)*u2 - 64*dyfacu(mz)*(s+u2)**(-2)*
     &    (1+m12/s4)*lq*hl**2*ssz*m1**2*t2*t2t**(-1) - 64*dyfacu(mz)*
     &    (s+u2)**(-2)*(1+m12/s4)*lq*hl**2*ssz*t2**2*t2t**(-1) - 64*
     &    dyfacu(mz)*(s+u2)**(-2)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2*
     &    t2t**(-1) - 64*dyfacu(mz)*(s+u2)**(-2)*(1+m12/s4)*rq*hr**2*
     &    ssz*t2**2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * ( 256*
     &    dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*
     &    s**(-1) + 256*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*pq*lq*ssz*
     &    ssp*s**(-1)*t2 + 256*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*pq*rq
     &    *ssz*ssp*m1**2*s**(-1) + 256*dyfacu(mz)*(s+u2)**(-1)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2 - 64*dyfacu(mz)*
     &    (s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2t**(-1) + 64*
     &    dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2**2
     &    *t2t**(-1) - 64*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*
     &    ssz*t2*t2t**(-1) - 64*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*rq*
     &    hr**2*ssz*m1**2*t2t**(-1) + 64*dyfacu(mz)*(s+u2)**(-1)*
     &    (1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2**2*t2t**(-1) - 64*dyfacu(
     &    mz)*(s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*t2*t2t**(-1) + 128*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-2) - 256*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-2) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * ( 128*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-2) - 256*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    s**(-2) + 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*
     &    lq*hl**2*ssz*mt**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)
     &    *(s+u2)*(1+m12/s4)*lq*hl**2*ssz*mz**2*s**(-1)*t2**(-1)*
     &    t2t**(-1) - 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*ssz*
     &    s**(-1)*t2**(-1) + 32*dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*
     &    ssz*s**(-1)*t2t**(-1) + 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*
     &    hr**2*ssz*m1**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*rq*hr**2*ssz*mt**2*s**(-1)*t2**(-1)*
     &    t2t**(-1) - 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*
     &    mz**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*(s+u2)*
     &    (1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2**(-1) + 32*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 64*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*ssz**2*lq2*s**(-1)*t2**(-1) - 64
     &    *dyfacu(mz)*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-1)*t2**(-1) + 
     &    16*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*lq*hl**2*ssz*mz**2*s**(-2)
     &    *t2**(-1)*t2t**(-1) + 16*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*rq*
     &    hr**2*ssz*mz**2*s**(-2)*t2**(-1)*t2t**(-1) + 64*dyfacu(mz)*
     &    (s+u2)**2*(1+m12/s4)*ssz**2*lq2*s**(-2)*t2**(-1) + 64*dyfacu(
     &    mz)*(s+u2)**2*(1+m12/s4)*ssz**2*rq2*s**(-2)*t2**(-1) + 128*
     &    dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*mz**(-2)*t2**(-1)
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-1) - 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*s**(-2)*t2
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*s**(-1) + 128*
     &    dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*mz**(-2)*t2**(-1)
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-1) - 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*t2
     &     + 256*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 16*
     &    dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2**(-1)*t2t**(-1)
     &     + 16*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*mt**2*t2**(-1)*
     &    t2t**(-1) + 64*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2*
     &    t2t**(-1) + 16*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*t2**(-1) - 
     &    32*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*t2t**(-1) - 16*dyfacu(
     &    mz)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2**(-1)*t2t**(-1) + 16*
     &    dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*mt**2*t2**(-1)*t2t**(-1)
     &     + 64*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2*t2t**(-1)
     &     + 16*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*t2**(-1) - 32*
     &    dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*t2t**(-1) + 128*dyfacu(mz)
     &    **2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*lq2*m1**2*s**(-1) + 128*
     &    dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*lq2*s**(-1)*t2
     &     + 128*dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*rq2*m1**2
     &    *s**(-1) + 128*dyfacu(mz)**2*(s+u2)**(-1)*(1+m12/s4)*ssz**2*
     &    rq2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 128*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*lq2*s**(-2) + 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*lq2*s**(-1)*t2**(-1)
     &     - 128*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-2) + 
     &    64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-1)*
     &    t2**(-1) + 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*lq2*
     &    mz**2*s**(-2)*t2**(-2) - 64*dyfacu(mz)**2*(s+u2)**2*
     &    (1+m12/s4)*ssz**2*lq2*s**(-2)*t2**(-1) + 64*dyfacu(mz)**2*
     &    (s+u2)**2*(1+m12/s4)*ssz**2*rq2*mz**2*s**(-2)*t2**(-2) - 64*
     &    dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*rq2*s**(-2)*
     &    t2**(-1) - 64*dyfacu(mz)**2*(s+u2)**3*(1+m12/s4)*ssz**2*lq2*
     &    mz**2*s**(-3)*t2**(-2) - 64*dyfacu(mz)**2*(s+u2)**3*
     &    (1+m12/s4)*ssz**2*rq2*mz**2*s**(-3)*t2**(-2) + 128*dyfacu(mz)
     &    **2*(1+m12/s4)*ssz**2*lq2*m1**2*s**(-1)*t2**(-1) - 128*
     &    dyfacu(mz)**2*(1+m12/s4)*ssz**2*lq2*s**(-2)*t2 + 128*dyfacu(
     &    mz)**2*(1+m12/s4)*ssz**2*lq2*s**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * ( 128*
     &    dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*t2**(-1) - 
     &    128*dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-2)*t2 + 128*
     &    dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-1) + 64*dyfacu(mh1)*
     &    (s+u2)**(-2)*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2*
     &    t2t**(-1) + 64*dyfacu(mh1)*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*t2t**(-1) - 32*dyfacu(mh1)*
     &    (mh12-mh22)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*t2**(-1)
     &     - 32*dyfacu(mh1)**2*(s+u2)**(-1)*(1+m12/s4)*h1**2*lambda1**2
     &    *s**(-1) - 32*dyfacu(mh1)**2*(1+m12/s4)*h1**2*lambda1**2*
     &    s**(-1)*t2**(-1) - 64*dyfacu(mh1)*dyfacu(mh2)*(s+u2)**(-1)*
     &    (1+m12/s4)*h1*h2*lambda1*lambda2*s**(-1) - 32*dyfacu(mh1)*
     &    dyfacu(mh2)*(s+u2)*(1+m12/s4)*h1*h2*lambda1*lambda2*s**(-1)*
     &    t2**(-2) - 64*dyfacu(mh1)*dyfacu(mh2)*(1+m12/s4)*h1*h2*
     &    lambda1*lambda2*s**(-1)*t2**(-1) + 64*dyfacu(mh2)*
     &    (s+u2)**(-2)*(1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2*
     &    t2t**(-1) )
      MMcrossed4 = MMcrossed4 + Nc*Cf*Pi**2*alphas*hardfac * ( 64*
     &    dyfacu(mh2)*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2t**(-1) + 32*dyfacu(mh2)*(mh12-mh22)**(-1)
     &    *(1+m12/s4)*h1*h2*lambda1*lambda2*t2**(-1) - 32*dyfacu(mh2)**
     &    2*(s+u2)**(-1)*(1+m12/s4)*h2**2*lambda2**2*s**(-1) - 32*
     &    dyfacu(mh2)**2*(1+m12/s4)*h2**2*lambda2**2*s**(-1)*t2**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * (  - 4*hl**2*hr**2*m1**2*mt**2*s**(-1)*u2*
     &    t2t**(-1) + 4*hl**2*hr**2*m1**2*mt**2*s*t2t**(-2) + 4*hl**2*
     &    hr**2*m1**2*mt**2*u2*t2t**(-2) + 4*hl**2*hr**2*mt**2*s**(-1)*
     &    t2 + 8*hl**2*hr**2*mt**2*s**(-1)*u2 + 4*hl**2*hr**2*mt**2*
     &    s**(-1)*u2**2*t2t**(-1) + 8*hl**2*hr**2*mt**2*s*u2*t2t**(-2)
     &     + 12*hl**2*hr**2*mt**2*s*t2t**(-1) + 4*hl**2*hr**2*mt**2*
     &    s**2*t2t**(-2) + 16*hl**2*hr**2*mt**2*u2*t2t**(-1) + 4*hl**2*
     &    hr**2*mt**2*u2**2*t2t**(-2) + 12*hl**2*hr**2*mt**2 + 4*hl**2*
     &    hr**2*mt**4*s**(-1)*u2*t2t**(-1) + 4*hl**2*hr**2*mt**4*s*
     &    t2t**(-2) + 4*hl**2*hr**2*mt**4*u2*t2t**(-2) + 8*hl**2*hr**2*
     &    mt**4*t2t**(-1) + 18*hl**4*m1**2*mt**2*s**(-1)*u2*t2t**(-1)
     &     + 4*hl**4*m1**2*mt**2*s**(-1) + 2*hl**4*m1**2*mt**2*s*
     &    t2t**(-2) + 2*hl**4*m1**2*mt**2*u2*t2t**(-2) + 20*hl**4*m1**2
     &    *mt**2*t2t**(-1) + 2*hl**4*m1**2*s**(-1)*t2 + 6*hl**4*m1**2*
     &    s**(-1)*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * ( 4*hl**4*m1**2*s**(-1)*u2**2*t2t**(-1) - 4*
     &    hl**4*m1**2*s*u2*t2t**(-2) + 2*hl**4*m1**2*s*t2t**(-1) - 2*
     &    hl**4*m1**2*s**2*t2t**(-2) + 6*hl**4*m1**2*u2*t2t**(-1) - 2*
     &    hl**4*m1**2*u2**2*t2t**(-2) + 6*hl**4*m1**2 - 8*hl**4*m1**4*
     &    s**(-1)*u2*t2t**(-1) - 2*hl**4*m1**4*s**(-1) - 8*hl**4*m1**4*
     &    t2t**(-1) - 4*hl**4*mt**2*s**(-1)*t2 - 10*hl**4*mt**2*s**(-1)
     &    *u2 - 6*hl**4*mt**2*s**(-1)*u2**2*t2t**(-1) - 8*hl**4*mt**2*s
     &    *t2t**(-1) - 14*hl**4*mt**2*u2*t2t**(-1) - 12*hl**4*mt**2 - 
     &    10*hl**4*mt**4*s**(-1)*u2*t2t**(-1) - 2*hl**4*mt**4*s**(-1)
     &     - 6*hl**4*mt**4*s*t2t**(-2) - 6*hl**4*mt**4*u2*t2t**(-2) - 
     &    16*hl**4*mt**4*t2t**(-1) - 8*hl**4*s**(-1)*t2*u2 - 4*hl**4*
     &    s**(-1)*t2**2 - 4*hl**4*s**(-1)*u2**2 - 4*hl**4*s - 8*hl**4*
     &    t2 - 8*hl**4*u2 + 18*hr**4*m1**2*mt**2*s**(-1)*u2*t2t**(-1)
     &     + 4*hr**4*m1**2*mt**2*s**(-1) + 2*hr**4*m1**2*mt**2*s*
     &    t2t**(-2) )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * ( 2*hr**4*m1**2*mt**2*u2*t2t**(-2) + 20*hr**4*
     &    m1**2*mt**2*t2t**(-1) + 2*hr**4*m1**2*s**(-1)*t2 + 6*hr**4*
     &    m1**2*s**(-1)*u2 + 4*hr**4*m1**2*s**(-1)*u2**2*t2t**(-1) - 4*
     &    hr**4*m1**2*s*u2*t2t**(-2) + 2*hr**4*m1**2*s*t2t**(-1) - 2*
     &    hr**4*m1**2*s**2*t2t**(-2) + 6*hr**4*m1**2*u2*t2t**(-1) - 2*
     &    hr**4*m1**2*u2**2*t2t**(-2) + 6*hr**4*m1**2 - 8*hr**4*m1**4*
     &    s**(-1)*u2*t2t**(-1) - 2*hr**4*m1**4*s**(-1) - 8*hr**4*m1**4*
     &    t2t**(-1) - 4*hr**4*mt**2*s**(-1)*t2 - 10*hr**4*mt**2*s**(-1)
     &    *u2 - 6*hr**4*mt**2*s**(-1)*u2**2*t2t**(-1) - 8*hr**4*mt**2*s
     &    *t2t**(-1) - 14*hr**4*mt**2*u2*t2t**(-1) - 12*hr**4*mt**2 - 
     &    10*hr**4*mt**4*s**(-1)*u2*t2t**(-1) - 2*hr**4*mt**4*s**(-1)
     &     - 6*hr**4*mt**4*s*t2t**(-2) - 6*hr**4*mt**4*u2*t2t**(-2) - 
     &    16*hr**4*mt**4*t2t**(-1) - 8*hr**4*s**(-1)*t2*u2 - 4*hr**4*
     &    s**(-1)*t2**2 - 4*hr**4*s**(-1)*u2**2 - 4*hr**4*s - 8*hr**4*
     &    t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * (  - 8*hr**4*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*pq*hl**2*ssp*m1**2*mt**2*s**(-1)*t2t**(-1)
     &     + 16*pq*hl**2*ssp*m1**2*s**(-1)*u2*t2t**(-1) - 16*pq*hl**2*
     &    ssp*m1**2*s**(-1) + 16*pq*hl**2*ssp*m1**2*t2t**(-1) - 24*pq*
     &    hl**2*ssp*mt**2*s**(-1)*u2*t2t**(-1) - 16*pq*hl**2*ssp*mt**2*
     &    s**(-1) - 8*pq*hl**2*ssp*mt**2*t2t**(-1) - 16*pq*hl**2*ssp*
     &    mt**4*s**(-1)*t2t**(-1) - 24*pq*hl**2*ssp*s**(-1)*t2 - 16*pq*
     &    hl**2*ssp*s**(-1)*u2 - 8*pq*hl**2*ssp + 16*pq*hr**2*ssp*m1**2
     &    *mt**2*s**(-1)*t2t**(-1) + 16*pq*hr**2*ssp*m1**2*s**(-1)*u2*
     &    t2t**(-1) - 16*pq*hr**2*ssp*m1**2*s**(-1) + 16*pq*hr**2*ssp*
     &    m1**2*t2t**(-1) - 24*pq*hr**2*ssp*mt**2*s**(-1)*u2*t2t**(-1)
     &     - 16*pq*hr**2*ssp*mt**2*s**(-1) - 8*pq*hr**2*ssp*mt**2*
     &    t2t**(-1) - 16*pq*hr**2*ssp*mt**4*s**(-1)*t2t**(-1) - 24*pq*
     &    hr**2*ssp*s**(-1)*t2 - 16*pq*hr**2*ssp*s**(-1)*u2 - 8*pq*
     &    hr**2*ssp + 16*lq*hl**2*ssz*m1**2*mt**2*s**(-1)*t2t**(-1) + 
     &    16*lq*hl**2*ssz*m1**2*s**(-1)*u2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*lq*hl**2*ssz*m1**2*s**(-1) + 16*lq*
     &    hl**2*ssz*m1**2*t2t**(-1) - 24*lq*hl**2*ssz*mt**2*s**(-1)*u2*
     &    t2t**(-1) - 16*lq*hl**2*ssz*mt**2*s**(-1) - 8*lq*hl**2*ssz*
     &    mt**2*t2t**(-1) - 16*lq*hl**2*ssz*mt**4*s**(-1)*t2t**(-1) - 
     &    24*lq*hl**2*ssz*s**(-1)*t2 - 16*lq*hl**2*ssz*s**(-1)*u2 - 8*
     &    lq*hl**2*ssz + 16*rq*hr**2*ssz*m1**2*mt**2*s**(-1)*t2t**(-1)
     &     + 16*rq*hr**2*ssz*m1**2*s**(-1)*u2*t2t**(-1) - 16*rq*hr**2*
     &    ssz*m1**2*s**(-1) + 16*rq*hr**2*ssz*m1**2*t2t**(-1) - 24*rq*
     &    hr**2*ssz*mt**2*s**(-1)*u2*t2t**(-1) - 16*rq*hr**2*ssz*mt**2*
     &    s**(-1) - 8*rq*hr**2*ssz*mt**2*t2t**(-1) - 16*rq*hr**2*ssz*
     &    mt**4*s**(-1)*t2t**(-1) - 24*rq*hr**2*ssz*s**(-1)*t2 - 16*rq*
     &    hr**2*ssz*s**(-1)*u2 - 8*rq*hr**2*ssz - 16*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s**(-1)*t2t**(-1) + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1)*u2*t2t**(-1) + 32*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2t**(-1)
     &     + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*s**(-1)*t2t**(-1) - 
     &    16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s**(-1)*t2t**(-1) + 
     &    8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)*u2*t2t**(-1) + 32*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2t**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*s**(-1)*t2t**(-1) + 4*hl**2*hr**2*m1**2*mt**2*t2t**(-2)
     &     - 8*hl**2*hr**2*mt**2*s*t2t**(-2) - 4*hl**2*hr**2*mt**2*u2*
     &    t2t**(-2) - 12*hl**2*hr**2*mt**2*t2t**(-1) - 4*hl**2*hr**2*
     &    mt**4*t2t**(-2) + 2*hl**4*m1**2*mt**2*t2t**(-2) + 2*hl**4*
     &    m1**2*s**(-1)*u2*t2t**(-1) + 4*hl**4*m1**2*s*t2t**(-2) + 4*
     &    hl**4*m1**2*u2*t2t**(-2) + 2*hl**4*m1**2*t2t**(-1) - 4*hl**4*
     &    m1**4*t2t**(-2) - 2*hl**4*mt**2*s**(-1)*u2*t2t**(-1) - 2*
     &    hl**4*mt**2*u2*t2t**(-2) + 4*hl**4*mt**2*t2t**(-1) + 2*hl**4*
     &    mt**4*t2t**(-2) )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 2*hl**4*s**(-1)*u2 - 2*hl**4*s**(-1)*u2**2
     &    *t2t**(-1) - 2*hl**4*u2*t2t**(-1) + 2*hl**4 + 2*hr**4*m1**2*
     &    mt**2*t2t**(-2) + 2*hr**4*m1**2*s**(-1)*u2*t2t**(-1) + 4*
     &    hr**4*m1**2*s*t2t**(-2) + 4*hr**4*m1**2*u2*t2t**(-2) + 2*
     &    hr**4*m1**2*t2t**(-1) - 4*hr**4*m1**4*t2t**(-2) - 2*hr**4*
     &    mt**2*s**(-1)*u2*t2t**(-1) - 2*hr**4*mt**2*u2*t2t**(-2) + 4*
     &    hr**4*mt**2*t2t**(-1) + 2*hr**4*mt**4*t2t**(-2) - 2*hr**4*
     &    s**(-1)*u2 - 2*hr**4*s**(-1)*u2**2*t2t**(-1) - 2*hr**4*u2*
     &    t2t**(-1) + 2*hr**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(0,0,0,0)*Nc*Cf*Pi*alphas*hardfac
     &  * (  - 8*pq*hl**2*ssp*m1**2*s**(-1)*t2t**(-1) - 8*pq*hl**2*ssp*
     &    mt**2*s**(-1)*t2t**(-1) - 8*pq*hl**2*ssp*s**(-1)*u2*t2t**(-1)
     &     - 8*pq*hl**2*ssp*s**(-1) - 8*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t2t**(-1) - 8*pq*hr**2*ssp*mt**2*s**(-1)*t2t**(-1) - 8*pq*
     &    hr**2*ssp*s**(-1)*u2*t2t**(-1) - 8*pq*hr**2*ssp*s**(-1) - 8*
     &    lq*hl**2*ssz*m1**2*s**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*
     &    s**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*s**(-1)*u2*t2t**(-1) - 8*
     &    lq*hl**2*ssz*s**(-1) - 8*rq*hr**2*ssz*m1**2*s**(-1)*t2t**(-1)
     &     - 8*rq*hr**2*ssz*mt**2*s**(-1)*t2t**(-1) - 8*rq*hr**2*ssz*
     &    s**(-1)*u2*t2t**(-1) - 8*rq*hr**2*ssz*s**(-1) + 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s**(-1)*t2t**(-1) + 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**(-1)*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*pq*hl**2*ssp*m1**2*mt**2*t2t**(-1) - 16*pq
     &    *hl**2*ssp*m1**2 - 8*pq*hl**2*ssp*mt**2*s*t2t**(-1) - 8*pq*
     &    hl**2*ssp*mt**2*u2*t2t**(-1) - 16*pq*hl**2*ssp*mt**2 - 16*pq*
     &    hl**2*ssp*mt**4*t2t**(-1) + 16*pq*hl**2*ssp*s**(-1)*t2*u2 + 
     &    16*pq*hl**2*ssp*s**(-1)*t2**2 + 8*pq*hl**2*ssp*s**(-1)*u2**2
     &     - 8*pq*hl**2*ssp*s - 8*pq*hl**2*ssp*t2 + 16*pq*hr**2*ssp*
     &    m1**2*mt**2*t2t**(-1) - 16*pq*hr**2*ssp*m1**2 - 8*pq*hr**2*
     &    ssp*mt**2*s*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*u2*t2t**(-1) - 
     &    16*pq*hr**2*ssp*mt**2 - 16*pq*hr**2*ssp*mt**4*t2t**(-1) + 16*
     &    pq*hr**2*ssp*s**(-1)*t2*u2 + 16*pq*hr**2*ssp*s**(-1)*t2**2 + 
     &    8*pq*hr**2*ssp*s**(-1)*u2**2 - 8*pq*hr**2*ssp*s - 8*pq*hr**2*
     &    ssp*t2 + 16*lq*hl**2*ssz*m1**2*mt**2*t2t**(-1) + 16*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1) - 16*lq*hl**2*ssz*m1**2 - 8*lq*hl**2*
     &    ssz*mt**2*s*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*u2*t2t**(-1) - 
     &    16*lq*hl**2*ssz*mt**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*lq*hl**2*ssz*mt**4*t2t**(-1) + 16*lq*
     &    hl**2*ssz*mz**2*s**(-1)*t2 + 8*lq*hl**2*ssz*mz**2*s**(-1)*u2
     &     + 8*lq*hl**2*ssz*mz**2 + 16*lq*hl**2*ssz*s**(-1)*t2*u2 + 16*
     &    lq*hl**2*ssz*s**(-1)*t2**2 + 8*lq*hl**2*ssz*s**(-1)*u2**2 - 8
     &    *lq*hl**2*ssz*s - 8*lq*hl**2*ssz*t2 + 16*rq*hr**2*ssz*m1**2*
     &    mt**2*t2t**(-1) + 16*rq*hr**2*ssz*m1**2*mz**2*s**(-1) - 16*rq
     &    *hr**2*ssz*m1**2 - 8*rq*hr**2*ssz*mt**2*s*t2t**(-1) - 8*rq*
     &    hr**2*ssz*mt**2*u2*t2t**(-1) - 16*rq*hr**2*ssz*mt**2 - 16*rq*
     &    hr**2*ssz*mt**4*t2t**(-1) + 16*rq*hr**2*ssz*mz**2*s**(-1)*t2
     &     + 8*rq*hr**2*ssz*mz**2*s**(-1)*u2 + 8*rq*hr**2*ssz*mz**2 + 
     &    16*rq*hr**2*ssz*s**(-1)*t2*u2 + 16*rq*hr**2*ssz*s**(-1)*t2**2
     &     + 8*rq*hr**2*ssz*s**(-1)*u2**2 - 8*rq*hr**2*ssz*s - 8*rq*
     &    hr**2*ssz*t2 - 16*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    t2t**(-1) - 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1)
     &     + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2*
     &    t2t**(-1) + 32*hl*hr*h1*lambda1*sqrt2**(-1)*mt + 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*t2t**(-1) - 16*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) - 16*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*s**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*u2*t2t**(-1) + 32*hl*hr*h2*lambda2*sqrt2**(-1)*mt + 16*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3*t2t**(-1) + 4*hl**2*hr**2*
     &    m1**2*mt**2*s*t2t**(-2) - 4*hl**2*hr**2*mt**2*s*t2t**(-1) - 4
     &    *hl**2*hr**2*mt**4*s*t2t**(-2) + 6*hl**4*m1**2*mt**2*s*
     &    t2t**(-2) + 8*hl**4*m1**2*mt**2*u2*t2t**(-2) + 24*hl**4*m1**2
     &    *mt**2*t2t**(-1) + 12*hl**4*m1**2*mt**4*t2t**(-2) + 6*hl**4*
     &    m1**2*s*t2t**(-1) + 6*hl**4*m1**2*u2*t2t**(-1) + 8*hl**4*
     &    m1**2 - 12*hl**4*m1**4*mt**2*t2t**(-2) - 4*hl**4*m1**4*s*
     &    t2t**(-2) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 4*hl**4*m1**4*u2*t2t**(-2) - 12*hl**4*
     &    m1**4*t2t**(-1) + 4*hl**4*m1**6*t2t**(-2) - 4*hl**4*mt**2*s*
     &    t2t**(-1) - 6*hl**4*mt**2*u2*t2t**(-1) - 8*hl**4*mt**2 - 2*
     &    hl**4*mt**4*s*t2t**(-2) - 4*hl**4*mt**4*u2*t2t**(-2) - 12*
     &    hl**4*mt**4*t2t**(-1) - 4*hl**4*mt**6*t2t**(-2) - 2*hl**4*s
     &     - 4*hl**4*t2 - 2*hl**4*u2 + 6*hr**4*m1**2*mt**2*s*t2t**(-2)
     &     + 8*hr**4*m1**2*mt**2*u2*t2t**(-2) + 24*hr**4*m1**2*mt**2*
     &    t2t**(-1) + 12*hr**4*m1**2*mt**4*t2t**(-2) + 6*hr**4*m1**2*s*
     &    t2t**(-1) + 6*hr**4*m1**2*u2*t2t**(-1) + 8*hr**4*m1**2 - 12*
     &    hr**4*m1**4*mt**2*t2t**(-2) - 4*hr**4*m1**4*s*t2t**(-2) - 4*
     &    hr**4*m1**4*u2*t2t**(-2) - 12*hr**4*m1**4*t2t**(-1) + 4*hr**4
     &    *m1**6*t2t**(-2) - 4*hr**4*mt**2*s*t2t**(-1) - 6*hr**4*mt**2*
     &    u2*t2t**(-1) - 8*hr**4*mt**2 - 2*hr**4*mt**4*s*t2t**(-2) - 4*
     &    hr**4*mt**4*u2*t2t**(-2) - 12*hr**4*mt**4*t2t**(-1) - 4*hr**4
     &    *mt**6*t2t**(-2) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 2*hr**4*s - 4*hr**4*t2 - 2*hr**4*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh1**2*s**(-1) + 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2
     &    *mh2**2*s**(-1) + 128*pq*lq*ssz*ssp*m1**2*s**(-1) + 128*pq*lq
     &    *ssz*ssp*s**(-1)*t2 + 64*pq*lq*ssz*ssp*s**(-1)*u2 + 64*pq*lq*
     &    ssz*ssp + 128*pq*rq*ssz*ssp*m1**2*s**(-1) + 128*pq*rq*ssz*ssp
     &    *s**(-1)*t2 + 64*pq*rq*ssz*ssp*s**(-1)*u2 + 64*pq*rq*ssz*ssp
     &     - 32*pq*hl**2*ssp*m1**2*mt**2*s**(-1)*t2t**(-1) - 16*pq*
     &    hl**2*ssp*m1**2*s**(-1)*u2*t2t**(-1) - 16*pq*hl**2*ssp*m1**2*
     &    s**(-1) - 24*pq*hl**2*ssp*m1**2*t2t**(-1) + 16*pq*hl**2*ssp*
     &    m1**4*s**(-1)*t2t**(-1) + 16*pq*hl**2*ssp*mt**2*s**(-1)*u2*
     &    t2t**(-1) + 16*pq*hl**2*ssp*mt**2*s**(-1) + 8*pq*hl**2*ssp*
     &    mt**2*t2t**(-1) + 16*pq*hl**2*ssp*mt**4*s**(-1)*t2t**(-1) + 
     &    16*pq*hl**2*ssp*s**(-1)*t2 + 16*pq*hl**2*ssp*s**(-1)*u2 + 8*
     &    pq*hl**2*ssp*s**(-1)*u2**2*t2t**(-1) + 8*pq*hl**2*ssp*u2*
     &    t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*pq*hl**2*ssp - 32*pq*hr**2*ssp*m1**2*mt**2*s**(-1)
     &    *t2t**(-1) - 16*pq*hr**2*ssp*m1**2*s**(-1)*u2*t2t**(-1) - 16*
     &    pq*hr**2*ssp*m1**2*s**(-1) - 24*pq*hr**2*ssp*m1**2*t2t**(-1)
     &     + 16*pq*hr**2*ssp*m1**4*s**(-1)*t2t**(-1) + 16*pq*hr**2*ssp*
     &    mt**2*s**(-1)*u2*t2t**(-1) + 16*pq*hr**2*ssp*mt**2*s**(-1) + 
     &    8*pq*hr**2*ssp*mt**2*t2t**(-1) + 16*pq*hr**2*ssp*mt**4*
     &    s**(-1)*t2t**(-1) + 16*pq*hr**2*ssp*s**(-1)*t2 + 16*pq*hr**2*
     &    ssp*s**(-1)*u2 + 8*pq*hr**2*ssp*s**(-1)*u2**2*t2t**(-1) + 8*
     &    pq*hr**2*ssp*u2*t2t**(-1) + 8*pq*hr**2*ssp - 32*lq*hl**2*ssz*
     &    m1**2*mt**2*s**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*m1**2*s**(-1)
     &    *u2*t2t**(-1) - 16*lq*hl**2*ssz*m1**2*s**(-1) - 24*lq*hl**2*
     &    ssz*m1**2*t2t**(-1) + 16*lq*hl**2*ssz*m1**4*s**(-1)*t2t**(-1)
     &     + 16*lq*hl**2*ssz*mt**2*mz**2*s**(-1)*t2t**(-1) + 16*lq*
     &    hl**2*ssz*mt**2*s**(-1)*u2*t2t**(-1) + 16*lq*hl**2*ssz*mt**2*
     &    s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*lq*hl**2*ssz*mt**2*t2t**(-1) + 16*lq*hl**2*ssz*
     &    mt**4*s**(-1)*t2t**(-1) + 8*lq*hl**2*ssz*mz**2*s**(-1)*u2*
     &    t2t**(-1) + 16*lq*hl**2*ssz*mz**2*s**(-1) + 8*lq*hl**2*ssz*
     &    mz**2*t2t**(-1) + 16*lq*hl**2*ssz*s**(-1)*t2 + 16*lq*hl**2*
     &    ssz*s**(-1)*u2 + 8*lq*hl**2*ssz*s**(-1)*u2**2*t2t**(-1) + 8*
     &    lq*hl**2*ssz*u2*t2t**(-1) + 8*lq*hl**2*ssz - 32*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*s**(-1)
     &    *u2*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*s**(-1) - 24*rq*hr**2*
     &    ssz*m1**2*t2t**(-1) + 16*rq*hr**2*ssz*m1**4*s**(-1)*t2t**(-1)
     &     + 16*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*t2t**(-1) + 16*rq*
     &    hr**2*ssz*mt**2*s**(-1)*u2*t2t**(-1) + 16*rq*hr**2*ssz*mt**2*
     &    s**(-1) + 8*rq*hr**2*ssz*mt**2*t2t**(-1) + 16*rq*hr**2*ssz*
     &    mt**4*s**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*mz**2*s**(-1)*u2*
     &    t2t**(-1) + 16*rq*hr**2*ssz*mz**2*s**(-1) + 8*rq*hr**2*ssz*
     &    mz**2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*rq*hr**2*ssz*s**(-1)*t2 + 16*rq*hr**2*ssz*s**(-1)
     &    *u2 + 8*rq*hr**2*ssz*s**(-1)*u2**2*t2t**(-1) + 8*rq*hr**2*ssz
     &    *u2*t2t**(-1) + 8*rq*hr**2*ssz - 16*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2*s**(-1)*t2t**(-1) + 16*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*t2t**(-1) - 16*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*mh2**2*s**(-1)*t2t**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt*t2t**(-1) - 4*hl**2*hr**2*mt**2*s*t2t**(-2) + 2*hl**4*
     &    m1**2*s*t2t**(-2) + 2*hl**4*m1**2*u2*t2t**(-2) - 2*hl**4*
     &    mt**2*u2*t2t**(-2) - 2*hl**4*u2*t2t**(-1) + 2*hr**4*m1**2*s*
     &    t2t**(-2) + 2*hr**4*m1**2*u2*t2t**(-2) - 2*hr**4*mt**2*u2*
     &    t2t**(-2) - 2*hr**4*u2*t2t**(-1) - 16*h1**2*lambda1**2*
     &    s**(-1) - 16*h2**2*lambda2**2*s**(-1) + 64*ssz**2*lq2*m1**2*
     &    s**(-1) + 64*ssz**2*lq2*s**(-1)*t2 + 32*ssz**2*lq2*s**(-1)*u2
     &     + 32*ssz**2*lq2 + 64*ssz**2*rq2*m1**2*s**(-1) + 64*ssz**2*
     &    rq2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*ssz**2*rq2*s**(-1)*u2 + 32*ssz**2*rq2 + 128*
     &    ssp**2*pq2*m1**2*s**(-1) + 128*ssp**2*pq2*s**(-1)*t2 + 64*
     &    ssp**2*pq2*s**(-1)*u2 + 64*ssp**2*pq2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*ssp**2*pq2*m1**2*s )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 24*pq*hl**2*ssp*m1**2*mt**2*s*t2t**(-1) + 32*
     &    pq*hl**2*ssp*m1**2*mt**2*u2*t2t**(-1) + 32*pq*hl**2*ssp*m1**2
     &    *mt**2 + 48*pq*hl**2*ssp*m1**2*mt**4*t2t**(-1) + 16*pq*hl**2*
     &    ssp*m1**2*s + 16*pq*hl**2*ssp*m1**2*t2 + 16*pq*hl**2*ssp*
     &    m1**2*u2 - 48*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1) - 16*pq*
     &    hl**2*ssp*m1**4*s*t2t**(-1) - 16*pq*hl**2*ssp*m1**4*u2*
     &    t2t**(-1) - 16*pq*hl**2*ssp*m1**4 + 16*pq*hl**2*ssp*m1**6*
     &    t2t**(-1) - 8*pq*hl**2*ssp*mt**2*s - 16*pq*hl**2*ssp*mt**2*t2
     &     - 16*pq*hl**2*ssp*mt**2*u2 - 8*pq*hl**2*ssp*mt**4*s*
     &    t2t**(-1) - 16*pq*hl**2*ssp*mt**4*u2*t2t**(-1) - 16*pq*hl**2*
     &    ssp*mt**4 - 16*pq*hl**2*ssp*mt**6*t2t**(-1) - 8*pq*hl**2*ssp*
     &    s*t2 - 8*pq*hl**2*ssp*t2*u2 - 16*pq*hl**2*ssp*t2**2 + 24*pq*
     &    hr**2*ssp*m1**2*mt**2*s*t2t**(-1) + 32*pq*hr**2*ssp*m1**2*
     &    mt**2*u2*t2t**(-1) + 32*pq*hr**2*ssp*m1**2*mt**2 + 48*pq*
     &    hr**2*ssp*m1**2*mt**4*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*pq*hr**2*ssp*m1**2*s + 16*pq*hr**2*ssp*
     &    m1**2*t2 + 16*pq*hr**2*ssp*m1**2*u2 - 48*pq*hr**2*ssp*m1**4*
     &    mt**2*t2t**(-1) - 16*pq*hr**2*ssp*m1**4*s*t2t**(-1) - 16*pq*
     &    hr**2*ssp*m1**4*u2*t2t**(-1) - 16*pq*hr**2*ssp*m1**4 + 16*pq*
     &    hr**2*ssp*m1**6*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*s - 16*pq*
     &    hr**2*ssp*mt**2*t2 - 16*pq*hr**2*ssp*mt**2*u2 - 8*pq*hr**2*
     &    ssp*mt**4*s*t2t**(-1) - 16*pq*hr**2*ssp*mt**4*u2*t2t**(-1) - 
     &    16*pq*hr**2*ssp*mt**4 - 16*pq*hr**2*ssp*mt**6*t2t**(-1) - 8*
     &    pq*hr**2*ssp*s*t2 - 8*pq*hr**2*ssp*t2*u2 - 16*pq*hr**2*ssp*
     &    t2**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s - 64*pq*rq*ssz*
     &    ssp*m1**2*mz**(-2)*s + 16*pq*hl**2*ssp*m1**2*s*t2t**(-1) + 8*
     &    pq*hl**2*ssp*m1**2*u2*t2t**(-1) - 8*pq*hl**2*ssp*mt**2*u2*
     &    t2t**(-1) - 8*pq*hl**2*ssp*u2 + 16*pq*hr**2*ssp*m1**2*s*
     &    t2t**(-1) + 8*pq*hr**2*ssp*m1**2*u2*t2t**(-1) - 8*pq*hr**2*
     &    ssp*mt**2*u2*t2t**(-1) - 8*pq*hr**2*ssp*u2 - 128*ssp**2*pq2*
     &    m1**2 + 128*ssp**2*pq2*s**(-1)*t2*u2 + 128*ssp**2*pq2*s**(-1)
     &    *t2**2 + 64*ssp**2*pq2*s**(-1)*u2**2 + 64*ssp**2*pq2*t2 + 64*
     &    ssp**2*pq2*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*pq*hl**2*ssp*m1**2*s**(-1) + 16*pq*hl**2*
     &    ssp*s**(-1)*t2 + 8*pq*hl**2*ssp*s**(-1)*u2 + 8*pq*hl**2*ssp
     &     + 16*pq*hr**2*ssp*m1**2*s**(-1) + 16*pq*hr**2*ssp*s**(-1)*t2
     &     + 8*pq*hr**2*ssp*s**(-1)*u2 + 8*pq*hr**2*ssp + 16*lq*hl**2*
     &    ssz*m1**2*s**(-1) + 16*lq*hl**2*ssz*s**(-1)*t2 + 8*lq*hl**2*
     &    ssz*s**(-1)*u2 + 8*lq*hl**2*ssz + 16*rq*hr**2*ssz*m1**2*
     &    s**(-1) + 16*rq*hr**2*ssz*s**(-1)*t2 + 8*rq*hr**2*ssz*s**(-1)
     &    *u2 + 8*rq*hr**2*ssz - 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**(-1) - 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1) - 8*
     &    hl**2*hr**2*m1**2*mt**2*t2t**(-2) + 4*hl**2*hr**2*mt**2*s*
     &    t2t**(-2) + 4*hl**2*hr**2*mt**2*u2*t2t**(-2) + 16*hl**2*hr**2
     &    *mt**2*t2t**(-1) + 8*hl**2*hr**2*mt**4*t2t**(-2) - 8*hl**4*
     &    m1**2*mt**2*s**(-1)*t2t**(-1) + 4*hl**4*m1**2*mt**2*t2t**(-2)
     &     - 4*hl**4*m1**2*s**(-1)*u2*t2t**(-1) - 4*hl**4*m1**2*s**(-1)
     &     - 2*hl**4*m1**2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 4*hl**4*m1**4*s**(-1)*t2t**(-1) + 4*hl**4*
     &    mt**2*s**(-1)*u2*t2t**(-1) + 4*hl**4*mt**2*s**(-1) - 2*hl**4*
     &    mt**2*s*t2t**(-2) - 2*hl**4*mt**2*u2*t2t**(-2) - 6*hl**4*
     &    mt**2*t2t**(-1) + 4*hl**4*mt**4*s**(-1)*t2t**(-1) - 4*hl**4*
     &    mt**4*t2t**(-2) + 4*hl**4*s**(-1)*t2 + 4*hl**4*s**(-1)*u2 + 2
     &    *hl**4*s**(-1)*u2**2*t2t**(-1) - 2*hl**4*s*t2t**(-1) - 2*
     &    hl**4 - 8*hr**4*m1**2*mt**2*s**(-1)*t2t**(-1) + 4*hr**4*m1**2
     &    *mt**2*t2t**(-2) - 4*hr**4*m1**2*s**(-1)*u2*t2t**(-1) - 4*
     &    hr**4*m1**2*s**(-1) - 2*hr**4*m1**2*t2t**(-1) + 4*hr**4*m1**4
     &    *s**(-1)*t2t**(-1) + 4*hr**4*mt**2*s**(-1)*u2*t2t**(-1) + 4*
     &    hr**4*mt**2*s**(-1) - 2*hr**4*mt**2*s*t2t**(-2) - 2*hr**4*
     &    mt**2*u2*t2t**(-2) - 6*hr**4*mt**2*t2t**(-1) + 4*hr**4*mt**4*
     &    s**(-1)*t2t**(-1) - 4*hr**4*mt**4*t2t**(-2) + 4*hr**4*s**(-1)
     &    *t2 + 4*hr**4*s**(-1)*u2 + 2*hr**4*s**(-1)*u2**2*t2t**(-1) - 
     &    2*hr**4*s*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 2*hr**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*pq*hl**2*ssp*mt**2*s**(-1)*t2t**(-1) + 8*pq*hl**2
     &    *ssp*s**(-1)*u2*t2t**(-1) + 16*pq*hl**2*ssp*s**(-1) + 8*pq*
     &    hl**2*ssp*t2t**(-1) + 16*pq*hr**2*ssp*mt**2*s**(-1)*t2t**(-1)
     &     + 8*pq*hr**2*ssp*s**(-1)*u2*t2t**(-1) + 16*pq*hr**2*ssp*
     &    s**(-1) + 8*pq*hr**2*ssp*t2t**(-1) + 16*lq*hl**2*ssz*mt**2*
     &    s**(-1)*t2t**(-1) + 8*lq*hl**2*ssz*s**(-1)*u2*t2t**(-1) + 16*
     &    lq*hl**2*ssz*s**(-1) + 8*lq*hl**2*ssz*t2t**(-1) + 16*rq*hr**2
     &    *ssz*mt**2*s**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*s**(-1)*u2*
     &    t2t**(-1) + 16*rq*hr**2*ssz*s**(-1) + 8*rq*hr**2*ssz*
     &    t2t**(-1) - 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1)*
     &    t2t**(-1) - 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)*
     &    t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,5,-1,2)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 8*hl**2*hr**2*mt**2*s**(-1)*t2t**(-1) + 4*
     &    hl**4*mt**2*s**(-1)*t2t**(-1) + 2*hl**4*s**(-1)*u2*t2t**(-1)
     &     + 4*hl**4*s**(-1) + 2*hl**4*t2t**(-1) + 4*hr**4*mt**2*
     &    s**(-1)*t2t**(-1) + 2*hr**4*s**(-1)*u2*t2t**(-1) + 4*hr**4*
     &    s**(-1) + 2*hr**4*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssz**2*lq2*m1**2*mz**2 + 64*ssz**2*lq2*m1**2*
     &    mz**4*s**(-1) + 32*ssz**2*lq2*m1**2*s + 64*ssz**2*lq2*mz**2*
     &    s**(-1)*t2*u2 + 64*ssz**2*lq2*mz**2*s**(-1)*t2**2 + 32*ssz**2
     &    *lq2*mz**2*s**(-1)*u2**2 + 32*ssz**2*lq2*mz**2*t2 + 32*ssz**2
     &    *lq2*mz**2*u2 + 64*ssz**2*lq2*mz**4*s**(-1)*t2 + 32*ssz**2*
     &    lq2*mz**4*s**(-1)*u2 + 32*ssz**2*lq2*mz**4 - 64*ssz**2*rq2*
     &    m1**2*mz**2 + 64*ssz**2*rq2*m1**2*mz**4*s**(-1) + 32*ssz**2*
     &    rq2*m1**2*s + 64*ssz**2*rq2*mz**2*s**(-1)*t2*u2 + 64*ssz**2*
     &    rq2*mz**2*s**(-1)*t2**2 + 32*ssz**2*rq2*mz**2*s**(-1)*u2**2
     &     + 32*ssz**2*rq2*mz**2*t2 + 32*ssz**2*rq2*mz**2*u2 + 64*
     &    ssz**2*rq2*mz**4*s**(-1)*t2 + 32*ssz**2*rq2*mz**4*s**(-1)*u2
     &     + 32*ssz**2*rq2*mz**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*lq*hl**2*ssz*m1**2*mt**2*mz**2*t2t**(-1)
     &     + 24*lq*hl**2*ssz*m1**2*mt**2*s*t2t**(-1) + 32*lq*hl**2*ssz*
     &    m1**2*mt**2*u2*t2t**(-1) + 32*lq*hl**2*ssz*m1**2*mt**2 + 48*
     &    lq*hl**2*ssz*m1**2*mt**4*t2t**(-1) - 16*lq*hl**2*ssz*m1**2*
     &    mz**2 + 16*lq*hl**2*ssz*m1**2*mz**4*s**(-1) + 16*lq*hl**2*ssz
     &    *m1**2*s + 16*lq*hl**2*ssz*m1**2*t2 + 16*lq*hl**2*ssz*m1**2*
     &    u2 - 48*lq*hl**2*ssz*m1**4*mt**2*t2t**(-1) - 16*lq*hl**2*ssz*
     &    m1**4*s*t2t**(-1) - 16*lq*hl**2*ssz*m1**4*u2*t2t**(-1) - 16*
     &    lq*hl**2*ssz*m1**4 + 16*lq*hl**2*ssz*m1**6*t2t**(-1) - 8*lq*
     &    hl**2*ssz*mt**2*mz**2*s*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*
     &    mz**2*u2*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*mz**2 - 8*lq*hl**2
     &    *ssz*mt**2*s - 16*lq*hl**2*ssz*mt**2*t2 - 16*lq*hl**2*ssz*
     &    mt**2*u2 - 16*lq*hl**2*ssz*mt**4*mz**2*t2t**(-1) - 8*lq*hl**2
     &    *ssz*mt**4*s*t2t**(-1) - 16*lq*hl**2*ssz*mt**4*u2*t2t**(-1)
     &     - 16*lq*hl**2*ssz*mt**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*lq*hl**2*ssz*mt**6*t2t**(-1) + 16*lq*
     &    hl**2*ssz*mz**2*s**(-1)*t2*u2 + 16*lq*hl**2*ssz*mz**2*s**(-1)
     &    *t2**2 + 8*lq*hl**2*ssz*mz**2*s**(-1)*u2**2 - 8*lq*hl**2*ssz*
     &    mz**2*s - 8*lq*hl**2*ssz*mz**2*t2 + 16*lq*hl**2*ssz*mz**4*
     &    s**(-1)*t2 + 8*lq*hl**2*ssz*mz**4*s**(-1)*u2 + 8*lq*hl**2*ssz
     &    *mz**4 - 8*lq*hl**2*ssz*s*t2 - 8*lq*hl**2*ssz*t2*u2 - 16*lq*
     &    hl**2*ssz*t2**2 + 16*rq*hr**2*ssz*m1**2*mt**2*mz**2*t2t**(-1)
     &     + 24*rq*hr**2*ssz*m1**2*mt**2*s*t2t**(-1) + 32*rq*hr**2*ssz*
     &    m1**2*mt**2*u2*t2t**(-1) + 32*rq*hr**2*ssz*m1**2*mt**2 + 48*
     &    rq*hr**2*ssz*m1**2*mt**4*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*
     &    mz**2 + 16*rq*hr**2*ssz*m1**2*mz**4*s**(-1) + 16*rq*hr**2*ssz
     &    *m1**2*s + 16*rq*hr**2*ssz*m1**2*t2 + 16*rq*hr**2*ssz*m1**2*
     &    u2 - 48*rq*hr**2*ssz*m1**4*mt**2*t2t**(-1) - 16*rq*hr**2*ssz*
     &    m1**4*s*t2t**(-1) - 16*rq*hr**2*ssz*m1**4*u2*t2t**(-1) - 16*
     &    rq*hr**2*ssz*m1**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*rq*hr**2*ssz*m1**6*t2t**(-1) - 8*rq*hr**2*
     &    ssz*mt**2*mz**2*s*t2t**(-1) - 8*rq*hr**2*ssz*mt**2*mz**2*u2*
     &    t2t**(-1) - 16*rq*hr**2*ssz*mt**2*mz**2 - 8*rq*hr**2*ssz*
     &    mt**2*s - 16*rq*hr**2*ssz*mt**2*t2 - 16*rq*hr**2*ssz*mt**2*u2
     &     - 16*rq*hr**2*ssz*mt**4*mz**2*t2t**(-1) - 8*rq*hr**2*ssz*
     &    mt**4*s*t2t**(-1) - 16*rq*hr**2*ssz*mt**4*u2*t2t**(-1) - 16*
     &    rq*hr**2*ssz*mt**4 - 16*rq*hr**2*ssz*mt**6*t2t**(-1) + 16*rq*
     &    hr**2*ssz*mz**2*s**(-1)*t2*u2 + 16*rq*hr**2*ssz*mz**2*s**(-1)
     &    *t2**2 + 8*rq*hr**2*ssz*mz**2*s**(-1)*u2**2 - 8*rq*hr**2*ssz*
     &    mz**2*s - 8*rq*hr**2*ssz*mz**2*t2 + 16*rq*hr**2*ssz*mz**4*
     &    s**(-1)*t2 + 8*rq*hr**2*ssz*mz**4*s**(-1)*u2 + 8*rq*hr**2*ssz
     &    *mz**4 - 8*rq*hr**2*ssz*s*t2 - 8*rq*hr**2*ssz*t2*u2 - 16*rq*
     &    hr**2*ssz*t2**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s + 128*pq*lq*ssz*
     &    ssp*m1**2*mz**2*s**(-1) - 128*pq*lq*ssz*ssp*m1**2 + 128*pq*lq
     &    *ssz*ssp*mz**2*s**(-1)*t2 + 64*pq*lq*ssz*ssp*mz**2*s**(-1)*u2
     &     + 64*pq*lq*ssz*ssp*mz**2 + 128*pq*lq*ssz*ssp*s**(-1)*t2*u2
     &     + 128*pq*lq*ssz*ssp*s**(-1)*t2**2 + 64*pq*lq*ssz*ssp*s**(-1)
     &    *u2**2 + 64*pq*lq*ssz*ssp*t2 + 64*pq*lq*ssz*ssp*u2 + 64*pq*rq
     &    *ssz*ssp*m1**2*mz**(-2)*s + 128*pq*rq*ssz*ssp*m1**2*mz**2*
     &    s**(-1) - 128*pq*rq*ssz*ssp*m1**2 + 128*pq*rq*ssz*ssp*mz**2*
     &    s**(-1)*t2 + 64*pq*rq*ssz*ssp*mz**2*s**(-1)*u2 + 64*pq*rq*ssz
     &    *ssp*mz**2 + 128*pq*rq*ssz*ssp*s**(-1)*t2*u2 + 128*pq*rq*ssz*
     &    ssp*s**(-1)*t2**2 + 64*pq*rq*ssz*ssp*s**(-1)*u2**2 + 64*pq*rq
     &    *ssz*ssp*t2 + 64*pq*rq*ssz*ssp*u2 - 32*lq*hl**2*ssz*m1**2*
     &    mt**2*mz**2*s**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*m1**2*mz**2*
     &    s**(-1)*u2*t2t**(-1) - 16*lq*hl**2*ssz*m1**2*mz**2*s**(-1) - 
     &    24*lq*hl**2*ssz*m1**2*mz**2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*lq*hl**2*ssz*m1**2*s*t2t**(-1) + 8*lq*hl**2*ssz*
     &    m1**2*u2*t2t**(-1) + 16*lq*hl**2*ssz*m1**4*mz**2*s**(-1)*
     &    t2t**(-1) + 16*lq*hl**2*ssz*mt**2*mz**2*s**(-1)*u2*t2t**(-1)
     &     + 16*lq*hl**2*ssz*mt**2*mz**2*s**(-1) + 8*lq*hl**2*ssz*mt**2
     &    *mz**2*t2t**(-1) + 16*lq*hl**2*ssz*mt**2*mz**4*s**(-1)*
     &    t2t**(-1) - 8*lq*hl**2*ssz*mt**2*u2*t2t**(-1) + 16*lq*hl**2*
     &    ssz*mt**4*mz**2*s**(-1)*t2t**(-1) + 16*lq*hl**2*ssz*mz**2*
     &    s**(-1)*t2 + 16*lq*hl**2*ssz*mz**2*s**(-1)*u2 + 8*lq*hl**2*
     &    ssz*mz**2*s**(-1)*u2**2*t2t**(-1) + 8*lq*hl**2*ssz*mz**2*u2*
     &    t2t**(-1) + 8*lq*hl**2*ssz*mz**2 + 8*lq*hl**2*ssz*mz**4*
     &    s**(-1)*u2*t2t**(-1) + 16*lq*hl**2*ssz*mz**4*s**(-1) + 8*lq*
     &    hl**2*ssz*mz**4*t2t**(-1) - 8*lq*hl**2*ssz*u2 - 32*rq*hr**2*
     &    ssz*m1**2*mt**2*mz**2*s**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*
     &    m1**2*mz**2*s**(-1)*u2*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*
     &    mz**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*rq*hr**2*ssz*m1**2*mz**2*t2t**(-1) + 16*rq*
     &    hr**2*ssz*m1**2*s*t2t**(-1) + 8*rq*hr**2*ssz*m1**2*u2*
     &    t2t**(-1) + 16*rq*hr**2*ssz*m1**4*mz**2*s**(-1)*t2t**(-1) + 
     &    16*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*u2*t2t**(-1) + 16*rq*
     &    hr**2*ssz*mt**2*mz**2*s**(-1) + 8*rq*hr**2*ssz*mt**2*mz**2*
     &    t2t**(-1) + 16*rq*hr**2*ssz*mt**2*mz**4*s**(-1)*t2t**(-1) - 8
     &    *rq*hr**2*ssz*mt**2*u2*t2t**(-1) + 16*rq*hr**2*ssz*mt**4*
     &    mz**2*s**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*mz**2*s**(-1)*t2 + 
     &    16*rq*hr**2*ssz*mz**2*s**(-1)*u2 + 8*rq*hr**2*ssz*mz**2*
     &    s**(-1)*u2**2*t2t**(-1) + 8*rq*hr**2*ssz*mz**2*u2*t2t**(-1)
     &     + 8*rq*hr**2*ssz*mz**2 + 8*rq*hr**2*ssz*mz**4*s**(-1)*u2*
     &    t2t**(-1) + 16*rq*hr**2*ssz*mz**4*s**(-1) + 8*rq*hr**2*ssz*
     &    mz**4*t2t**(-1) - 8*rq*hr**2*ssz*u2 + 128*ssz**2*lq2*m1**2*
     &    mz**2*s**(-1) - 64*ssz**2*lq2*m1**2 + 128*ssz**2*lq2*mz**2*
     &    s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*ssz**2*lq2*mz**2*s**(-1)*u2 + 64*ssz**2*lq2*mz**2
     &     + 64*ssz**2*lq2*s**(-1)*t2*u2 + 64*ssz**2*lq2*s**(-1)*t2**2
     &     + 32*ssz**2*lq2*s**(-1)*u2**2 + 32*ssz**2*lq2*t2 + 32*ssz**2
     &    *lq2*u2 + 128*ssz**2*rq2*m1**2*mz**2*s**(-1) - 64*ssz**2*rq2*
     &    m1**2 + 128*ssz**2*rq2*mz**2*s**(-1)*t2 + 64*ssz**2*rq2*mz**2
     &    *s**(-1)*u2 + 64*ssz**2*rq2*mz**2 + 64*ssz**2*rq2*s**(-1)*t2*
     &    u2 + 64*ssz**2*rq2*s**(-1)*t2**2 + 32*ssz**2*rq2*s**(-1)*
     &    u2**2 + 32*ssz**2*rq2*t2 + 32*ssz**2*rq2*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,11,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h1**2*lambda1**2*mh1**2 - 16*h1**2*lambda1**2*
     &    mh1**4*s**(-1) - 8*h1**2*lambda1**2*s )
      MMcrossed4 = MMcrossed4 + ANGfin(2,11,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    mh1**2*t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s*
     &    t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s*
     &    t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*u2*
     &    t2t**(-1) + 32*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2 - 16*hl
     &    *hr*h1*lambda1*sqrt2**(-1)*mt*mh1**4*s**(-1) - 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s + 16*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*mh1**2*t2t**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2
     &     - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**4*
     &    s**(-1) - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 16
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*t2t**(-1) - 16*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*mh1**4*s**(-1)*t2t**(-1) - 16*hl*hr
     &    *h1*lambda1*sqrt2**(-1)*mt*s*t2t**(-1) - 32*h1**2*lambda1**2*
     &    mh1**2*s**(-1) + 16*h1**2*lambda1**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(2,12,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h2**2*lambda2**2*mh2**2 - 16*h2**2*lambda2**2*
     &    mh2**4*s**(-1) - 8*h2**2*lambda2**2*s )
      MMcrossed4 = MMcrossed4 + ANGfin(2,12,-1,-1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    mh2**2*t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s*
     &    t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s*
     &    t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*u2*
     &    t2t**(-1) + 32*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2 - 16*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*mt*mh2**4*s**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s + 16*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*mh2**2*t2t**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(2,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2 + 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh2**4*
     &    s**(-1) + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 16
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*t2t**(-1) - 16*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt*mh2**4*s**(-1)*t2t**(-1) - 16*hl*hr
     &    *h2*lambda2*sqrt2**(-1)*mt*s*t2t**(-1) - 32*h2**2*lambda2**2*
     &    mh2**2*s**(-1) + 16*h2**2*lambda2**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssp**2*pq2*m1**2*s**(-1)*u2 - 64*ssp**2*pq2*
     &    m1**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 72*pq*hl**2*ssp*m1**2*mt**2*s**(-1)*u2*
     &    t2t**(-1) + 32*pq*hl**2*ssp*m1**2*mt**2*s**(-1) + 56*pq*hl**2
     &    *ssp*m1**2*mt**2*t2t**(-1) + 48*pq*hl**2*ssp*m1**2*mt**4*
     &    s**(-1)*t2t**(-1) + 16*pq*hl**2*ssp*m1**2*s**(-1)*t2 + 16*pq*
     &    hl**2*ssp*m1**2*s**(-1)*u2 + 16*pq*hl**2*ssp*m1**2*s**(-1)*
     &    u2**2*t2t**(-1) + 16*pq*hl**2*ssp*m1**2*s*t2t**(-1) + 32*pq*
     &    hl**2*ssp*m1**2*u2*t2t**(-1) + 32*pq*hl**2*ssp*m1**2 - 48*pq*
     &    hl**2*ssp*m1**4*mt**2*s**(-1)*t2t**(-1) - 32*pq*hl**2*ssp*
     &    m1**4*s**(-1)*u2*t2t**(-1) - 16*pq*hl**2*ssp*m1**4*s**(-1) - 
     &    48*pq*hl**2*ssp*m1**4*t2t**(-1) + 16*pq*hl**2*ssp*m1**6*
     &    s**(-1)*t2t**(-1) - 16*pq*hl**2*ssp*mt**2*s**(-1)*t2 - 40*pq*
     &    hl**2*ssp*mt**2*s**(-1)*u2 - 24*pq*hl**2*ssp*mt**2*s**(-1)*
     &    u2**2*t2t**(-1) - 16*pq*hl**2*ssp*mt**2*u2*t2t**(-1) - 8*pq*
     &    hl**2*ssp*mt**2 - 40*pq*hl**2*ssp*mt**4*s**(-1)*u2*t2t**(-1)
     &     - 16*pq*hl**2*ssp*mt**4*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 8*pq*hl**2*ssp*mt**4*t2t**(-1) - 16*pq*
     &    hl**2*ssp*mt**6*s**(-1)*t2t**(-1) - 32*pq*hl**2*ssp*s**(-1)*
     &    t2*u2 - 16*pq*hl**2*ssp*s**(-1)*t2**2 - 16*pq*hl**2*ssp*
     &    s**(-1)*u2**2 - 8*pq*hl**2*ssp*t2 - 16*pq*hl**2*ssp*u2 + 72*
     &    pq*hr**2*ssp*m1**2*mt**2*s**(-1)*u2*t2t**(-1) + 32*pq*hr**2*
     &    ssp*m1**2*mt**2*s**(-1) + 56*pq*hr**2*ssp*m1**2*mt**2*
     &    t2t**(-1) + 48*pq*hr**2*ssp*m1**2*mt**4*s**(-1)*t2t**(-1) + 
     &    16*pq*hr**2*ssp*m1**2*s**(-1)*t2 + 16*pq*hr**2*ssp*m1**2*
     &    s**(-1)*u2 + 16*pq*hr**2*ssp*m1**2*s**(-1)*u2**2*t2t**(-1) + 
     &    16*pq*hr**2*ssp*m1**2*s*t2t**(-1) + 32*pq*hr**2*ssp*m1**2*u2*
     &    t2t**(-1) + 32*pq*hr**2*ssp*m1**2 - 48*pq*hr**2*ssp*m1**4*
     &    mt**2*s**(-1)*t2t**(-1) - 32*pq*hr**2*ssp*m1**4*s**(-1)*u2*
     &    t2t**(-1) - 16*pq*hr**2*ssp*m1**4*s**(-1) - 48*pq*hr**2*ssp*
     &    m1**4*t2t**(-1) + 16*pq*hr**2*ssp*m1**6*s**(-1)*t2t**(-1) - 
     &    16*pq*hr**2*ssp*mt**2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 40*pq*hr**2*ssp*mt**2*s**(-1)*u2 - 24*pq*
     &    hr**2*ssp*mt**2*s**(-1)*u2**2*t2t**(-1) - 16*pq*hr**2*ssp*
     &    mt**2*u2*t2t**(-1) - 8*pq*hr**2*ssp*mt**2 - 40*pq*hr**2*ssp*
     &    mt**4*s**(-1)*u2*t2t**(-1) - 16*pq*hr**2*ssp*mt**4*s**(-1) - 
     &    8*pq*hr**2*ssp*mt**4*t2t**(-1) - 16*pq*hr**2*ssp*mt**6*
     &    s**(-1)*t2t**(-1) - 32*pq*hr**2*ssp*s**(-1)*t2*u2 - 16*pq*
     &    hr**2*ssp*s**(-1)*t2**2 - 16*pq*hr**2*ssp*s**(-1)*u2**2 - 8*
     &    pq*hr**2*ssp*t2 - 16*pq*hr**2*ssp*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s**(-1)*u2 + 64*pq*
     &    lq*ssz*ssp*m1**2*mz**(-2) + 64*pq*rq*ssz*ssp*m1**2*mz**(-2)*
     &    s**(-1)*u2 + 64*pq*rq*ssz*ssp*m1**2*mz**(-2) + 8*pq*hl**2*ssp
     &    *m1**2*s**(-1)*u2*t2t**(-1) - 8*pq*hl**2*ssp*mt**2*s**(-1)*u2
     &    *t2t**(-1) - 8*pq*hl**2*ssp*s**(-1)*u2 - 8*pq*hl**2*ssp*
     &    s**(-1)*u2**2*t2t**(-1) + 8*pq*hr**2*ssp*m1**2*s**(-1)*u2*
     &    t2t**(-1) - 8*pq*hr**2*ssp*mt**2*s**(-1)*u2*t2t**(-1) - 8*pq*
     &    hr**2*ssp*s**(-1)*u2 - 8*pq*hr**2*ssp*s**(-1)*u2**2*t2t**(-1)
     &     - 128*ssp**2*pq2*m1**2*s**(-1) - 64*ssp**2*pq2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,1,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * ( 4*hl**2*hr**2*mt**2*s**(-1)*u2*t2t**(-1) + 4*
     &    hl**2*hr**2*mt**2*s**(-1) + 4*hl**2*hr**2*mt**2*s*t2t**(-2)
     &     + 4*hl**2*hr**2*mt**2*u2*t2t**(-2) + 8*hl**2*hr**2*mt**2*
     &    t2t**(-1) + 4*hl**4*m1**2*s**(-1)*u2*t2t**(-1) + 4*hl**4*
     &    m1**2*t2t**(-1) - 6*hl**4*mt**2*s**(-1)*u2*t2t**(-1) - 2*
     &    hl**4*mt**2*s**(-1) - 2*hl**4*mt**2*s*t2t**(-2) - 2*hl**4*
     &    mt**2*u2*t2t**(-2) - 8*hl**4*mt**2*t2t**(-1) - 4*hl**4*
     &    s**(-1)*t2 - 4*hl**4*s**(-1)*u2 - 2*hl**4*s*t2t**(-1) - 2*
     &    hl**4*u2*t2t**(-1) - 6*hl**4 + 4*hr**4*m1**2*s**(-1)*u2*
     &    t2t**(-1) + 4*hr**4*m1**2*t2t**(-1) - 6*hr**4*mt**2*s**(-1)*
     &    u2*t2t**(-1) - 2*hr**4*mt**2*s**(-1) - 2*hr**4*mt**2*s*
     &    t2t**(-2) - 2*hr**4*mt**2*u2*t2t**(-2) - 8*hr**4*mt**2*
     &    t2t**(-1) - 4*hr**4*s**(-1)*t2 - 4*hr**4*s**(-1)*u2 - 2*hr**4
     &    *s*t2t**(-1) - 2*hr**4*u2*t2t**(-1) - 6*hr**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,0,1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 8*hl**2*hr**2*mt**2*s**(-1)*t2t**(-1) - 2*
     &    hl**4*m1**2*s**(-1)*t2t**(-1) - 2*hl**4*mt**2*s**(-1)*
     &    t2t**(-1) - 2*hl**4*s**(-1)*u2*t2t**(-1) - 2*hl**4*s**(-1) - 
     &    2*hr**4*m1**2*s**(-1)*t2t**(-1) - 2*hr**4*mt**2*s**(-1)*
     &    t2t**(-1) - 2*hr**4*s**(-1)*u2*t2t**(-1) - 2*hr**4*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(5,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssp**2*pq2*m1**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(5,7,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 8*pq*hl**2*ssp*m1**2*mt**2*s**(-1)*t2t**(-1)
     &     - 16*pq*hl**2*ssp*m1**2*s**(-1) - 8*pq*hl**2*ssp*mt**2*
     &    s**(-1)*u2*t2t**(-1) - 8*pq*hl**2*ssp*mt**2*s**(-1) - 8*pq*
     &    hl**2*ssp*mt**4*s**(-1)*t2t**(-1) - 8*pq*hl**2*ssp*s**(-1)*t2
     &     - 8*pq*hl**2*ssp*s**(-1)*u2 + 8*pq*hr**2*ssp*m1**2*mt**2*
     &    s**(-1)*t2t**(-1) - 16*pq*hr**2*ssp*m1**2*s**(-1) - 8*pq*
     &    hr**2*ssp*mt**2*s**(-1)*u2*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*
     &    s**(-1) - 8*pq*hr**2*ssp*mt**4*s**(-1)*t2t**(-1) - 8*pq*hr**2
     &    *ssp*s**(-1)*t2 - 8*pq*hr**2*ssp*s**(-1)*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(5,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s**(-1) + 64*pq*rq*
     &    ssz*ssp*m1**2*mz**(-2)*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(7,0,1,0)*Nc*Cf*s4t**(-2)*Pi*
     & alphas*hardfac * ( 4*hl**2*hr**2*mt**2*s**(-1)*u2*t2t**(-1) + 4*
     &    hl**2*hr**2*mt**2*s**(-1) + 4*hl**2*hr**2*mt**2*s*t2t**(-2)
     &     + 4*hl**2*hr**2*mt**2*u2*t2t**(-2) + 8*hl**2*hr**2*mt**2*
     &    t2t**(-1) - 2*hl**4*m1**2*s**(-1) - 2*hl**4*m1**2*s*t2t**(-2)
     &     - 2*hl**4*m1**2*u2*t2t**(-2) - 2*hl**4*m1**2*t2t**(-1) - 2*
     &    hl**4*mt**2*s**(-1)*u2*t2t**(-1) - 2*hl**4*mt**2*t2t**(-1) - 
     &    2*hl**4*s**(-1)*t2 - 2*hl**4*s**(-1)*u2 - 2*hl**4 - 2*hr**4*
     &    m1**2*s**(-1) - 2*hr**4*m1**2*s*t2t**(-2) - 2*hr**4*m1**2*u2*
     &    t2t**(-2) - 2*hr**4*m1**2*t2t**(-1) - 2*hr**4*mt**2*s**(-1)*
     &    u2*t2t**(-1) - 2*hr**4*mt**2*t2t**(-1) - 2*hr**4*s**(-1)*t2
     &     - 2*hr**4*s**(-1)*u2 - 2*hr**4 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssz**2*lq2*m1**2*mz**2*s**(-1) - 32*ssz**2*lq2
     &    *m1**2*s**(-1)*u2 - 32*ssz**2*lq2*m1**2 - 32*ssz**2*lq2*mz**2
     &    *s**(-1)*t2 - 64*ssz**2*rq2*m1**2*mz**2*s**(-1) - 32*ssz**2*
     &    rq2*m1**2*s**(-1)*u2 - 32*ssz**2*rq2*m1**2 - 32*ssz**2*rq2*
     &    mz**2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 16*lq*hl**2*ssz*m1**2*mt**2*mz**2*s**(-1)*
     &    t2t**(-1) + 72*lq*hl**2*ssz*m1**2*mt**2*s**(-1)*u2*t2t**(-1)
     &     + 32*lq*hl**2*ssz*m1**2*mt**2*s**(-1) + 56*lq*hl**2*ssz*
     &    m1**2*mt**2*t2t**(-1) + 48*lq*hl**2*ssz*m1**2*mt**4*s**(-1)*
     &    t2t**(-1) + 16*lq*hl**2*ssz*m1**2*mz**2*s**(-1)*u2*t2t**(-1)
     &     - 16*lq*hl**2*ssz*m1**2*mz**2*s**(-1) + 16*lq*hl**2*ssz*
     &    m1**2*mz**2*t2t**(-1) + 16*lq*hl**2*ssz*m1**2*s**(-1)*t2 + 16
     &    *lq*hl**2*ssz*m1**2*s**(-1)*u2 + 16*lq*hl**2*ssz*m1**2*
     &    s**(-1)*u2**2*t2t**(-1) + 16*lq*hl**2*ssz*m1**2*s*t2t**(-1)
     &     + 32*lq*hl**2*ssz*m1**2*u2*t2t**(-1) + 32*lq*hl**2*ssz*m1**2
     &     - 48*lq*hl**2*ssz*m1**4*mt**2*s**(-1)*t2t**(-1) - 32*lq*
     &    hl**2*ssz*m1**4*s**(-1)*u2*t2t**(-1) - 16*lq*hl**2*ssz*m1**4*
     &    s**(-1) - 48*lq*hl**2*ssz*m1**4*t2t**(-1) + 16*lq*hl**2*ssz*
     &    m1**6*s**(-1)*t2t**(-1) - 24*lq*hl**2*ssz*mt**2*mz**2*s**(-1)
     &    *u2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*lq*hl**2*ssz*mt**2*mz**2*s**(-1) - 8*lq
     &    *hl**2*ssz*mt**2*mz**2*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*
     &    s**(-1)*t2 - 40*lq*hl**2*ssz*mt**2*s**(-1)*u2 - 24*lq*hl**2*
     &    ssz*mt**2*s**(-1)*u2**2*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*u2*
     &    t2t**(-1) - 8*lq*hl**2*ssz*mt**2 - 16*lq*hl**2*ssz*mt**4*
     &    mz**2*s**(-1)*t2t**(-1) - 40*lq*hl**2*ssz*mt**4*s**(-1)*u2*
     &    t2t**(-1) - 16*lq*hl**2*ssz*mt**4*s**(-1) - 8*lq*hl**2*ssz*
     &    mt**4*t2t**(-1) - 16*lq*hl**2*ssz*mt**6*s**(-1)*t2t**(-1) - 
     &    24*lq*hl**2*ssz*mz**2*s**(-1)*t2 - 16*lq*hl**2*ssz*mz**2*
     &    s**(-1)*u2 - 8*lq*hl**2*ssz*mz**2 - 32*lq*hl**2*ssz*s**(-1)*
     &    t2*u2 - 16*lq*hl**2*ssz*s**(-1)*t2**2 - 16*lq*hl**2*ssz*
     &    s**(-1)*u2**2 - 8*lq*hl**2*ssz*t2 - 16*lq*hl**2*ssz*u2 + 16*
     &    rq*hr**2*ssz*m1**2*mt**2*mz**2*s**(-1)*t2t**(-1) + 72*rq*
     &    hr**2*ssz*m1**2*mt**2*s**(-1)*u2*t2t**(-1) + 32*rq*hr**2*ssz*
     &    m1**2*mt**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 56*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) + 48*rq
     &    *hr**2*ssz*m1**2*mt**4*s**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*
     &    m1**2*mz**2*s**(-1)*u2*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*
     &    mz**2*s**(-1) + 16*rq*hr**2*ssz*m1**2*mz**2*t2t**(-1) + 16*rq
     &    *hr**2*ssz*m1**2*s**(-1)*t2 + 16*rq*hr**2*ssz*m1**2*s**(-1)*
     &    u2 + 16*rq*hr**2*ssz*m1**2*s**(-1)*u2**2*t2t**(-1) + 16*rq*
     &    hr**2*ssz*m1**2*s*t2t**(-1) + 32*rq*hr**2*ssz*m1**2*u2*
     &    t2t**(-1) + 32*rq*hr**2*ssz*m1**2 - 48*rq*hr**2*ssz*m1**4*
     &    mt**2*s**(-1)*t2t**(-1) - 32*rq*hr**2*ssz*m1**4*s**(-1)*u2*
     &    t2t**(-1) - 16*rq*hr**2*ssz*m1**4*s**(-1) - 48*rq*hr**2*ssz*
     &    m1**4*t2t**(-1) + 16*rq*hr**2*ssz*m1**6*s**(-1)*t2t**(-1) - 
     &    24*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*u2*t2t**(-1) - 16*rq*
     &    hr**2*ssz*mt**2*mz**2*s**(-1) - 8*rq*hr**2*ssz*mt**2*mz**2*
     &    t2t**(-1) - 16*rq*hr**2*ssz*mt**2*s**(-1)*t2 - 40*rq*hr**2*
     &    ssz*mt**2*s**(-1)*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 24*rq*hr**2*ssz*mt**2*s**(-1)*u2**2*
     &    t2t**(-1) - 16*rq*hr**2*ssz*mt**2*u2*t2t**(-1) - 8*rq*hr**2*
     &    ssz*mt**2 - 16*rq*hr**2*ssz*mt**4*mz**2*s**(-1)*t2t**(-1) - 
     &    40*rq*hr**2*ssz*mt**4*s**(-1)*u2*t2t**(-1) - 16*rq*hr**2*ssz*
     &    mt**4*s**(-1) - 8*rq*hr**2*ssz*mt**4*t2t**(-1) - 16*rq*hr**2*
     &    ssz*mt**6*s**(-1)*t2t**(-1) - 24*rq*hr**2*ssz*mz**2*s**(-1)*
     &    t2 - 16*rq*hr**2*ssz*mz**2*s**(-1)*u2 - 8*rq*hr**2*ssz*mz**2
     &     - 32*rq*hr**2*ssz*s**(-1)*t2*u2 - 16*rq*hr**2*ssz*s**(-1)*
     &    t2**2 - 16*rq*hr**2*ssz*s**(-1)*u2**2 - 8*rq*hr**2*ssz*t2 - 
     &    16*rq*hr**2*ssz*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s**(-1)*u2 - 64*
     &    pq*lq*ssz*ssp*m1**2*mz**(-2) - 128*pq*lq*ssz*ssp*m1**2*
     &    s**(-1) - 64*pq*lq*ssz*ssp*s**(-1)*t2 - 64*pq*rq*ssz*ssp*
     &    m1**2*mz**(-2)*s**(-1)*u2 - 64*pq*rq*ssz*ssp*m1**2*mz**(-2)
     &     - 128*pq*rq*ssz*ssp*m1**2*s**(-1) - 64*pq*rq*ssz*ssp*s**(-1)
     &    *t2 - 8*lq*hl**2*ssz*m1**2*mz**2*s**(-1)*t2t**(-1) + 8*lq*
     &    hl**2*ssz*m1**2*s**(-1)*u2*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*
     &    mz**2*s**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*s**(-1)*u2*
     &    t2t**(-1) - 8*lq*hl**2*ssz*mz**2*s**(-1)*u2*t2t**(-1) - 8*lq*
     &    hl**2*ssz*mz**2*s**(-1) - 8*lq*hl**2*ssz*s**(-1)*u2 - 8*lq*
     &    hl**2*ssz*s**(-1)*u2**2*t2t**(-1) - 8*rq*hr**2*ssz*m1**2*
     &    mz**2*s**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*m1**2*s**(-1)*u2*
     &    t2t**(-1) - 8*rq*hr**2*ssz*mt**2*mz**2*s**(-1)*t2t**(-1) - 8*
     &    rq*hr**2*ssz*mt**2*s**(-1)*u2*t2t**(-1) - 8*rq*hr**2*ssz*
     &    mz**2*s**(-1)*u2*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*rq*hr**2*ssz*mz**2*s**(-1) - 8*rq*hr**2*ssz*
     &    s**(-1)*u2 - 8*rq*hr**2*ssz*s**(-1)*u2**2*t2t**(-1) - 64*
     &    ssz**2*lq2*m1**2*s**(-1) - 32*ssz**2*lq2*s**(-1)*t2 - 64*
     &    ssz**2*rq2*m1**2*s**(-1) - 32*ssz**2*rq2*s**(-1)*t2 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*ssz**2*lq2*m1**2*s**(-1) - 32*ssz**2*rq2*m1**2
     &    *s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(10,7,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * ( 8*lq*hl**2*ssz*m1**2*mt**2*s**(-1)*t2t**(-1)
     &     - 16*lq*hl**2*ssz*m1**2*s**(-1) - 8*lq*hl**2*ssz*mt**2*
     &    s**(-1)*u2*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*s**(-1) - 8*lq*
     &    hl**2*ssz*mt**4*s**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*s**(-1)*t2
     &     - 8*lq*hl**2*ssz*s**(-1)*u2 + 8*rq*hr**2*ssz*m1**2*mt**2*
     &    s**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*s**(-1) - 8*rq*
     &    hr**2*ssz*mt**2*s**(-1)*u2*t2t**(-1) - 8*rq*hr**2*ssz*mt**2*
     &    s**(-1) - 8*rq*hr**2*ssz*mt**4*s**(-1)*t2t**(-1) - 8*rq*hr**2
     &    *ssz*s**(-1)*t2 - 8*rq*hr**2*ssz*s**(-1)*u2 )
      MMcrossed4 = MMcrossed4 + ANGfin(10,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*m1**2*mz**(-2)*s**(-1) - 64*pq*
     &    rq*ssz*ssp*m1**2*mz**(-2)*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(11,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h1**2*lambda1**2*mh1**2*s**(-1) + 8*h1**2*
     &    lambda1**2*s**(-1)*u2 - 8*h1**2*lambda1**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(11,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    mh1**2*s**(-1)*t2t**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*u2*t2t**(-1) + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2*s**(-1)*u2*t2t**(-1) + 32*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1) + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2*t2t**(-1) + 24*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**(-1)*u2 + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s**(-1)*u2**2*t2t**(-1) + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *u2*t2t**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt + 16*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt**3*mh1**2*s**(-1)*t2t**(-1) + 8*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**3*s**(-1)*u2*t2t**(-1) - 8*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**3*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(11,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2*
     &    s**(-1) + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    s**(-1)*u2 - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2 + 
     &    16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1)*t2t**(-1)
     &     - 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2t**(-1) + 16*h1**2*
     &    lambda1**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(11,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h1**2*lambda1**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(11,7,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 8*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    s**(-1)*t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1)
     &    *u2*t2t**(-1) + 24*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**(-1) + 
     &    16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2t**(-1) + 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s**(-1)*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(11,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**(-1)
     &     )
      MMcrossed4 = MMcrossed4 + ANGfin(12,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*h2**2*lambda2**2*mh2**2*s**(-1) + 8*h2**2*
     &    lambda2**2*s**(-1)*u2 - 8*h2**2*lambda2**2 )
      MMcrossed4 = MMcrossed4 + ANGfin(12,0,-1,0)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    mh2**2*s**(-1)*t2t**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    m1**2*mt*s**(-1)*u2*t2t**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*s**(-1)*u2*t2t**(-1) + 32*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*mh2**2*s**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*t2t**(-1) + 24*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**(-1)*u2 + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *s**(-1)*u2**2*t2t**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt
     &    *u2*t2t**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt + 16*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt**3*mh2**2*s**(-1)*t2t**(-1) + 8*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3*s**(-1)*u2*t2t**(-1) - 8*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(12,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2*s**(-1) - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2
     &    *s**(-1)*u2 + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2 + 
     &    16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s**(-1)*t2t**(-1)
     &     - 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2t**(-1) + 16*h2**2*
     &    lambda2**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(12,7,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h2**2*lambda2**2*s**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(12,7,-1,1)*Nc*Cf*s4t**(-1)*Pi*
     & alphas*hardfac * (  - 8*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    s**(-1)*t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1)
     &    *u2*t2t**(-1) + 24*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**(-1) + 
     &    16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2t**(-1) + 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*s**(-1)*t2t**(-1) )
      MMcrossed4 = MMcrossed4 + ANGfin(12,7,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    s**(-1) )

c               the phase space except for 1/s**2 
      HH_GBH = MMcrossed4 / ( 16.D0 * pi**2 )**2 / 2.D0*s4/(s4+m1**2)

c               the averaging factors
      HH_GBH = HH_GBH /4.D0 /Nc/(Nc**2-1.D0)

c               the prefactor for the scaling functions 
      HH_GBH = HH_GBH * (m1+m2)**2/4.D0 

      end

      































