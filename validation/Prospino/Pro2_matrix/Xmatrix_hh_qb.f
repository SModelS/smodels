cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     HH_QBH(MASSIN,C)                                                 c
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
      real*8 function HH_QBH(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Nc,Cf,sqrt2,alphas
     &          ,m1,m2,mt,mz,mh1,mh2,m12,mt2,mz2,mh12,mh22
     &          ,ssp,ssz,hl,hr 
     &          ,h1,h2,lambda1,lambda2
     &          ,lq,rq,pq
     &          ,lq2,rq2,pq2
     &          ,s,s1,s2,sz,s4,t2,t2t,topfac,u2,epsbw
     &          ,hardfac,logall
     &          ,dyfact,dyfacu
     &          ,MMpartial4,ANGfin(0:12,0:12,-2:2,-2:2)

      external dyfact, dyfacu

c          common block to talk to DYFACT and DYFACU
      common/HH_QBB_INT/s,t2,u2

      Pi    = 4.D0*atan(1.D0)
      sqrt2 = sqrt(2.D0)
      Nc    = 3.D0
      Cf    = 4.D0/3.D0

      s     = massin(1)
      t2    = massin(2)
      s4    = massin(3)
      m1    = massin(6)
      m2    = massin(6)
      mt    = massin(7)
      mz    = massin(8)
      mh1   = massin(9)
      mh2   = massin(10)
      epsbw = massin(26)

c               real kinematics built in
      u2  = s4 - s - t2 - m1**2 + m1**2 
      t2t = t2+m1**2-mt**2

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
      logall = log(massin(13)**2/m1**2) 
     &        - log(s4/m1**2) + log(1.D0+m1**2/s4)

c               denominators 
      sz = s - mz2
      s1 = s - mh12
      s2 = s - mh22
      
      if ( sz .gt. 0.D0 ) then
         sz =  sqrt( sz**2 + mz**4 *epsbw**2 )
      else 
         sz = -sqrt( sz**2 + mz**4 *epsbw**2 )
      end if

      if ( s1 .gt. 0.D0 ) then
         s1 =  sqrt( s1**2 + mh1**4*epsbw**2 )
      else 
         s1 = -sqrt( s1**2 + mh1**4*epsbw**2 )
      end if

      if ( s2 .gt. 0.D0 ) then
         s2 =  sqrt( s2**2 + mh2**4*epsbw**2 )
      else 
         s2 = -sqrt( s2**2 + mh2**4*epsbw**2 )
      end if

c               set gs=1 
      alphas = 1.D0/(4.D0*Pi) 

      hardfac = 1.D0

c               remaining special function without argument 
      topfac = 1.D0+(mt2-m12)*(s+t2)/t2/u2
      topfac = 1.D0/topfac

c               the angular functions 
      call ANGULAR_ARRAY_HH_QB(massin,ANGfin)

c               form output
      MMpartial4 =
     &  + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*logall * ( 128*
     &    (s+t2)**(-1)*(1+m12/s4)*ssp**2*pq2*m1**2*s**(-1)*t2 + 128*
     &    (s+t2)**(-1)*(1+m12/s4)*ssp**2*pq2*m1**2 - 128*(s+t2)**(-1)*
     &    (1+m12/s4)*ssp**2*pq2*s**(-2)*t2**2*u2 - 128*(s+t2)**(-1)*
     &    (1+m12/s4)*ssp**2*pq2*s**(-1)*t2*u2 + 128*(s+t2)*(1+m12/s4)*
     &    ssp**2*pq2*m1**2*s**(-1)*t2*u2**(-2) + 128*(s+t2)*(1+m12/s4)*
     &    ssp**2*pq2*m1**2*u2**(-2) - 128*(s+t2)*(1+m12/s4)*ssp**2*pq2*
     &    s**(-2)*t2**2*u2**(-1) - 128*(s+t2)*(1+m12/s4)*ssp**2*pq2*
     &    s**(-1)*t2*u2**(-1) - 8*(s+u2)**(-3)*(1+m12/s4)*hl**2*hr**2*
     &    mt**2*s*t2**2*u2*t2t**(-2) - 8*(s+u2)**(-3)*(1+m12/s4)*hl**2*
     &    hr**2*mt**2*s**2*t2**2*t2t**(-2) + 4*(s+u2)**(-3)*(1+m12/s4)*
     &    hl**4*m1**2*s*t2**2*u2*t2t**(-2) + 4*(s+u2)**(-3)*(1+m12/s4)*
     &    hl**4*m1**2*s**2*t2**2*t2t**(-2) - 4*(s+u2)**(-3)*(1+m12/s4)*
     &    hl**4*s*t2**3*u2*t2t**(-2) - 4*(s+u2)**(-3)*(1+m12/s4)*hl**4*
     &    t2**3*u2**2*t2t**(-2) + 4*(s+u2)**(-3)*(1+m12/s4)*hr**4*m1**2
     &    *s*t2**2*u2*t2t**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * ( 4*(s+u2)**(-3)*(1+m12/s4)*hr**4*m1**2*s**2*t2**2*
     &    t2t**(-2) - 4*(s+u2)**(-3)*(1+m12/s4)*hr**4*s*t2**3*u2*
     &    t2t**(-2) - 4*(s+u2)**(-3)*(1+m12/s4)*hr**4*t2**3*u2**2*
     &    t2t**(-2) - 32*(s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*s*
     &    t2*t2t**(-1) - 32*(s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*
     &    t2*u2*t2t**(-1) + 32*(s+u2)**(-2)*(1+m12/s4)*pq*hl**2*ssp*
     &    s**(-1)*t2**2*u2**2*t2t**(-1) + 32*(s+u2)**(-2)*(1+m12/s4)*pq
     &    *hl**2*ssp*t2**2*u2*t2t**(-1) - 32*(s+u2)**(-2)*(1+m12/s4)*pq
     &    *hr**2*ssp*m1**2*s*t2*t2t**(-1) - 32*(s+u2)**(-2)*(1+m12/s4)*
     &    pq*hr**2*ssp*m1**2*t2*u2*t2t**(-1) + 32*(s+u2)**(-2)*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*t2**2*u2**2*t2t**(-1) + 32*
     &    (s+u2)**(-2)*(1+m12/s4)*pq*hr**2*ssp*t2**2*u2*t2t**(-1) - 8*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**2*hr**2*mt**2*s*u2*t2t**(-2) - 8*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**2*hr**2*mt**2*s**2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*m1**2*s*u2*t2t**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * ( 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*m1**2*s**2*t2t**(-2)
     &     - 4*(s+u2)**(-1)*(1+m12/s4)*hl**4*s*t2*u2*t2t**(-2) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*t2*u2**2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*m1**2*s*u2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*m1**2*s**2*t2t**(-2) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*s*t2*u2*t2t**(-2) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hr**4*t2*u2**2*t2t**(-2) + 128*
     &    (s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*m1**2*s**(-1)*u2 + 128*
     &    (s+u2)**(-1)*(1+m12/s4)*ssp**2*pq2*m1**2 - 128*(s+u2)**(-1)*
     &    (1+m12/s4)*ssp**2*pq2*s**(-2)*t2*u2**2 - 128*(s+u2)**(-1)*
     &    (1+m12/s4)*ssp**2*pq2*s**(-1)*t2*u2 + 128*(s+u2)*(1+m12/s4)*
     &    ssp**2*pq2*m1**2*s**(-1)*t2**(-2)*u2 + 128*(s+u2)*(1+m12/s4)*
     &    ssp**2*pq2*m1**2*t2**(-2) - 128*(s+u2)*(1+m12/s4)*ssp**2*pq2*
     &    s**(-2)*t2**(-1)*u2**2 - 128*(s+u2)*(1+m12/s4)*ssp**2*pq2*
     &    s**(-1)*t2**(-1)*u2 )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 32*(1+m12/s4)*pq*hl**2*ssp*m1**2*s*t2**(-1)*
     &    t2t**(-1) - 32*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2**(-1)*u2*
     &    t2t**(-1) + 32*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**2*
     &    t2t**(-1) + 32*(1+m12/s4)*pq*hl**2*ssp*u2*t2t**(-1) - 32*
     &    (1+m12/s4)*pq*hr**2*ssp*m1**2*s*t2**(-1)*t2t**(-1) - 32*
     &    (1+m12/s4)*pq*hr**2*ssp*m1**2*t2**(-1)*u2*t2t**(-1) + 32*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2**2*t2t**(-1) + 32*
     &    (1+m12/s4)*pq*hr**2*ssp*u2*t2t**(-1) - 128*dyfact(mz)*(s+t2)*
     &    (1+m12/s4)*pq*lq*ssz*ssp*s**(-2)*u2 - 128*dyfact(mz)*(s+t2)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*u2 + 128*dyfact(mz)*
     &    (s+t2)**2*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1)*u2**(-2) + 
     &    128*dyfact(mz)*(s+t2)**2*(1+m12/s4)*pq*lq*ssz*ssp*s**(-1)*
     &    u2**(-1) + 128*dyfact(mz)*(s+t2)**2*(1+m12/s4)*pq*rq*ssz*ssp*
     &    m1**2*s**(-1)*u2**(-2) + 128*dyfact(mz)*(s+t2)**2*(1+m12/s4)*
     &    pq*rq*ssz*ssp*s**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 128*dyfact(mz)*(s+t2)**3*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-2)*u2**(-1) - 128*dyfact(mz)*(s+t2)**3*(1+m12/s4)*pq*rq*
     &    ssz*ssp*s**(-2)*u2**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*lq*
     &    ssz*ssp*m1**2*s**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*lq*ssz*
     &    ssp*s**(-1)*u2 + 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    m1**2*s**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    s**(-1)*u2 - 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*lq2*
     &    s**(-2)*u2 - 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*rq2*
     &    s**(-2)*u2 + 64*dyfact(mz)**2*(s+t2)**2*(1+m12/s4)*ssz**2*lq2
     &    *m1**2*s**(-1)*u2**(-2) + 64*dyfact(mz)**2*(s+t2)**2*
     &    (1+m12/s4)*ssz**2*lq2*s**(-1)*u2**(-1) + 64*dyfact(mz)**2*
     &    (s+t2)**2*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*u2**(-2) + 64*
     &    dyfact(mz)**2*(s+t2)**2*(1+m12/s4)*ssz**2*rq2*s**(-1)*
     &    u2**(-1) - 64*dyfact(mz)**2*(s+t2)**3*(1+m12/s4)*ssz**2*lq2*
     &    s**(-2)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 64*dyfact(mz)**2*(s+t2)**3*(1+m12/s4)*ssz**2*rq2*
     &    s**(-2)*u2**(-1) + 64*dyfact(mz)**2*(1+m12/s4)*ssz**2*lq2*
     &    m1**2*s**(-1) + 64*dyfact(mz)**2*(1+m12/s4)*ssz**2*lq2*
     &    s**(-1)*u2 + 64*dyfact(mz)**2*(1+m12/s4)*ssz**2*rq2*m1**2*
     &    s**(-1) + 64*dyfact(mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-1)*u2
     &     + 32*dyfact(mz)*topfac*(s+t2)**2*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*t2**(-1)*u2**(-2) - 32*dyfact(mz)*topfac*(s+t2)**2*
     &    (1+m12/s4)*lq*hl**2*ssz*s**(-1)*u2**(-1) + 32*dyfact(mz)*
     &    topfac*(s+t2)**2*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2**(-1)*
     &    u2**(-2) - 32*dyfact(mz)*topfac*(s+t2)**2*(1+m12/s4)*rq*hr**2
     &    *ssz*s**(-1)*u2**(-1) + 32*dyfact(mz)*topfac*(1+m12/s4)*lq*
     &    hl**2*ssz*m1**2*t2**(-1) - 32*dyfact(mz)*topfac*(1+m12/s4)*lq
     &    *hl**2*ssz*s**(-1)*u2 + 32*dyfact(mz)*topfac*(1+m12/s4)*rq*
     &    hr**2*ssz*m1**2*t2**(-1) - 32*dyfact(mz)*topfac*(1+m12/s4)*rq
     &    *hr**2*ssz*s**(-1)*u2 )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 16*dyfact(mh1)**2*(s+t2)**2*(1+m12/s4)*h1**2*
     &    lambda1**2*s**(-1)*u2**(-2) - 16*dyfact(mh1)**2*(1+m12/s4)*
     &    h1**2*lambda1**2*s**(-1) - 32*dyfact(mh1)*dyfact(mh2)*
     &    (s+t2)**2*(1+m12/s4)*h1*h2*lambda1*lambda2*s**(-1)*u2**(-2)
     &     - 32*dyfact(mh1)*dyfact(mh2)*(1+m12/s4)*h1*h2*lambda1*
     &    lambda2*s**(-1) - 32*dyfact(mh1)*topfac*(s+t2)**2*(1+m12/s4)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*u2**(-2) - 32*
     &    dyfact(mh1)*topfac*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *t2**(-1) - 16*dyfact(mh2)**2*(s+t2)**2*(1+m12/s4)*h2**2*
     &    lambda2**2*s**(-1)*u2**(-2) - 16*dyfact(mh2)**2*(1+m12/s4)*
     &    h2**2*lambda2**2*s**(-1) - 32*dyfact(mh2)*topfac*(s+t2)**2*
     &    (1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2**(-1)*u2**(-2)
     &     - 32*dyfact(mh2)*topfac*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2**(-1) - 32*dyfacu(mz)*(s+u2)**(-1)*
     &    (1+m12/s4)*lq*hl**2*ssz*m1**2*t2*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 32*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz
     &    *t2**2*t2t**(-1) - 32*dyfacu(mz)*(s+u2)**(-1)*(1+m12/s4)*rq*
     &    hr**2*ssz*m1**2*t2*t2t**(-1) - 32*dyfacu(mz)*(s+u2)**(-1)*
     &    (1+m12/s4)*rq*hr**2*ssz*t2**2*t2t**(-1) - 128*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*s**(-2)*t2 - 128*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*t2 - 32*dyfacu(mz)*
     &    (s+u2)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2**(-1)*t2t**(-1) - 32*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*lq*hl**2*ssz*t2t**(-1) - 32*
     &    dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2**(-1)*
     &    t2t**(-1) - 32*dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*
     &    t2t**(-1) + 128*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*pq*lq*ssz*ssp
     &    *m1**2*s**(-1)*t2**(-2) + 128*dyfacu(mz)*(s+u2)**2*(1+m12/s4)
     &    *pq*lq*ssz*ssp*s**(-1)*t2**(-1) + 128*dyfacu(mz)*(s+u2)**2*
     &    (1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1)*t2**(-2) + 128*dyfacu(
     &    mz)*(s+u2)**2*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * ( 32*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*lq*hl**2*ssz*
     &    s**(-1)*t2t**(-1) + 32*dyfacu(mz)*(s+u2)**2*(1+m12/s4)*rq*
     &    hr**2*ssz*s**(-1)*t2t**(-1) - 128*dyfacu(mz)*(s+u2)**3*
     &    (1+m12/s4)*pq*lq*ssz*ssp*s**(-2)*t2**(-1) - 128*dyfacu(mz)*
     &    (s+u2)**3*(1+m12/s4)*pq*rq*ssz*ssp*s**(-2)*t2**(-1) + 128*
     &    dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1) + 128*
     &    dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*s**(-1)*t2 + 128*dyfacu(
     &    mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1) + 128*dyfacu(mz)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2 + 32*dyfacu(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2**2*t2t**(-1) + 32*dyfacu(
     &    mz)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2**2*t2t**(-1) - 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*lq2*s**(-2)*t2 - 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-2)*t2 + 64*
     &    dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*lq2*m1**2*s**(-1)*
     &    t2**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * ( 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*lq2*
     &    s**(-1)*t2**(-1) + 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*
     &    ssz**2*rq2*m1**2*s**(-1)*t2**(-2) + 64*dyfacu(mz)**2*
     &    (s+u2)**2*(1+m12/s4)*ssz**2*rq2*s**(-1)*t2**(-1) - 64*dyfacu(
     &    mz)**2*(s+u2)**3*(1+m12/s4)*ssz**2*lq2*s**(-2)*t2**(-1) - 64*
     &    dyfacu(mz)**2*(s+u2)**3*(1+m12/s4)*ssz**2*rq2*s**(-2)*
     &    t2**(-1) + 64*dyfacu(mz)**2*(1+m12/s4)*ssz**2*lq2*m1**2*
     &    s**(-1) + 64*dyfacu(mz)**2*(1+m12/s4)*ssz**2*lq2*s**(-1)*t2
     &     + 64*dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1) + 64*
     &    dyfacu(mz)**2*(1+m12/s4)*ssz**2*rq2*s**(-1)*t2 + 32*dyfacu(
     &    mh1)*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    t2*t2t**(-1) + 32*dyfacu(mh1)*(s+u2)*(1+m12/s4)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*t2**(-1)*t2t**(-1) - 16*dyfacu(mh1)**2
     &    *(s+u2)**2*(1+m12/s4)*h1**2*lambda1**2*s**(-1)*t2**(-2) - 16*
     &    dyfacu(mh1)**2*(1+m12/s4)*h1**2*lambda1**2*s**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 32*dyfacu(mh1)*dyfacu(mh2)*(s+u2)**2*(1+m12/s4)*h1
     &    *h2*lambda1*lambda2*s**(-1)*t2**(-2) - 32*dyfacu(mh1)*dyfacu(
     &    mh2)*(1+m12/s4)*h1*h2*lambda1*lambda2*s**(-1) + 32*dyfacu(mh2
     &    )*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2*
     &    t2t**(-1) + 32*dyfacu(mh2)*(s+u2)*(1+m12/s4)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*t2**(-1)*t2t**(-1) - 16*dyfacu(mh2)**2*
     &    (s+u2)**2*(1+m12/s4)*h2**2*lambda2**2*s**(-1)*t2**(-2) - 16*
     &    dyfacu(mh2)**2*(1+m12/s4)*h2**2*lambda2**2*s**(-1) + 32*
     &    topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*s*t2**(-1)
     &     + 32*topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2 - 32*
     &    topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2*u2 - 
     &    32*topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*u2 + 32*topfac
     &    *(s+t2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*s*t2**(-1) + 32*
     &    topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2 - 32*topfac
     &    *(s+t2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*t2*u2 )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 32*topfac*(s+t2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*u2
     &     + 32*topfac*(s+t2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*s*t2**(-1)*
     &    u2**(-2) + 32*topfac*(s+t2)*(1+m12/s4)*pq*hl**2*ssp*m1**2*
     &    u2**(-2) - 32*topfac*(s+t2)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*
     &    t2*u2**(-1) - 32*topfac*(s+t2)*(1+m12/s4)*pq*hl**2*ssp*
     &    u2**(-1) + 32*topfac*(s+t2)*(1+m12/s4)*pq*hr**2*ssp*m1**2*s*
     &    t2**(-1)*u2**(-2) + 32*topfac*(s+t2)*(1+m12/s4)*pq*hr**2*ssp*
     &    m1**2*u2**(-2) - 32*topfac*(s+t2)*(1+m12/s4)*pq*hr**2*ssp*
     &    s**(-1)*t2*u2**(-1) - 32*topfac*(s+t2)*(1+m12/s4)*pq*hr**2*
     &    ssp*u2**(-1) - 8*topfac**2*(s+t2)**(-1)*(1+m12/s4)*hl**2*
     &    hr**2*mt**2*s*t2**(-1) - 8*topfac**2*(s+t2)**(-1)*(1+m12/s4)*
     &    hl**2*hr**2*mt**2*s**2*t2**(-2) + 4*topfac**2*(s+t2)**(-1)*
     &    (1+m12/s4)*hl**4*m1**2*s*t2**(-1) + 4*topfac**2*(s+t2)**(-1)*
     &    (1+m12/s4)*hl**4*m1**2*s**2*t2**(-2) - 4*topfac**2*
     &    (s+t2)**(-1)*(1+m12/s4)*hl**4*s*t2**(-1)*u2 )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 4*topfac**2*(s+t2)**(-1)*(1+m12/s4)*hl**4*u2 + 4*
     &    topfac**2*(s+t2)**(-1)*(1+m12/s4)*hr**4*m1**2*s*t2**(-1) + 4*
     &    topfac**2*(s+t2)**(-1)*(1+m12/s4)*hr**4*m1**2*s**2*t2**(-2)
     &     - 4*topfac**2*(s+t2)**(-1)*(1+m12/s4)*hr**4*s*t2**(-1)*u2 - 
     &    4*topfac**2*(s+t2)**(-1)*(1+m12/s4)*hr**4*u2 - 8*topfac**2*
     &    (s+t2)*(1+m12/s4)*hl**2*hr**2*mt**2*s*t2**(-1)*u2**(-2) - 8*
     &    topfac**2*(s+t2)*(1+m12/s4)*hl**2*hr**2*mt**2*s**2*t2**(-2)*
     &    u2**(-2) + 4*topfac**2*(s+t2)*(1+m12/s4)*hl**4*m1**2*s*
     &    t2**(-1)*u2**(-2) + 4*topfac**2*(s+t2)*(1+m12/s4)*hl**4*m1**2
     &    *s**2*t2**(-2)*u2**(-2) - 4*topfac**2*(s+t2)*(1+m12/s4)*hl**4
     &    *s*t2**(-1)*u2**(-1) - 4*topfac**2*(s+t2)*(1+m12/s4)*hl**4*
     &    u2**(-1) + 4*topfac**2*(s+t2)*(1+m12/s4)*hr**4*m1**2*s*
     &    t2**(-1)*u2**(-2) + 4*topfac**2*(s+t2)*(1+m12/s4)*hr**4*m1**2
     &    *s**2*t2**(-2)*u2**(-2) - 4*topfac**2*(s+t2)*(1+m12/s4)*hr**4
     &    *s*t2**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**(-1)*Pi**2*alphas*hardfac*
     & logall * (  - 4*topfac**2*(s+t2)*(1+m12/s4)*hr**4*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 8*
     &    (s+u2)**(-2)*(1+m12/s4)*hl**2*hr**2*mt**2*s*t2t**(-2) - 4*
     &    (s+u2)**(-2)*(1+m12/s4)*hl**4*m1**2*s*t2t**(-2) - 4*
     &    (s+u2)**(-2)*(1+m12/s4)*hr**4*m1**2*s*t2t**(-2) - 16*
     &    (s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*
     &    m1**2*t2t**(-1) + 16*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*mt**2*t2t**(-1) - 16*(s+u2)**(-1)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) + 16*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *pq*hr**2*ssp*mt**2*t2t**(-1) - 16*(s+u2)**(-1)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*
     &    t2t**(-1) + 16*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *lq*hl**2*ssz*mt**2*t2t**(-1) - 16*(s+u2)**(-1)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*
     &    t2t**(-1) + 16*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *rq*hr**2*ssz*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 4*
     &    (s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hl**4*m1**2*
     &    t2t**(-1) + 4*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*
     &    hl**4*mt**2*t2t**(-1) - 4*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)
     &    *(1+m12/s4)*hr**4*m1**2*t2t**(-1) + 4*(s+u2)**(-1)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*hr**4*mt**2*t2t**(-1) + 16*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*t2t**(-1) + 16*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*t2t**(-1) + 16*
     &    (s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*t2t**(-1) + 16*
     &    (s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*t2t**(-1) - 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*m1**2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*mt**2*t2t**(-2) + 4*
     &    (s+u2)**(-1)*(1+m12/s4)*hl**4*s*t2t**(-2) - 4*(s+u2)**(-1)*
     &    (1+m12/s4)*hr**4*m1**2*t2t**(-2) + 4*(s+u2)**(-1)*(1+m12/s4)*
     &    hr**4*mt**2*t2t**(-2) + 4*(s+u2)**(-1)*(1+m12/s4)*hr**4*s*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 32*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*s**(-1)*
     &    u2**(-1) - 32*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*
     &    mt**2*s**(-1)*u2**(-1) + 16*(u2-m12+mt2)**(-1)*(1+m12/s4)*
     &    pq*hl**2*ssp*s**(-1)*t2*u2**(-1) + 16*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*u2**(-1) + 32*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hr**2*ssp*m1**2*s**(-1)*u2**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*mt**2*s**(-1)*
     &    u2**(-1) + 16*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*
     &    s**(-1)*t2*u2**(-1) + 16*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*
     &    hr**2*ssp*u2**(-1) + 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*
     &    hl**2*ssp*m1**2*s**(-1)*t2**(-1) - 16*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*m1**2*t2**(-1)*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*mt**2*s**(-1)*
     &    t2**(-1) + 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*
     &    mt**2*t2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 32*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t2**(-1) - 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*
     &    m1**2*t2**(-1)*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hr**2*ssp*mt**2*s**(-1)*t2**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*mt**2*t2**(-1)
     &    *t2t**(-1) + 32*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*
     &    ssp*s**(-1) + 32*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*
     &    ssz*sz**(-1) + 32*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*
     &    ssz*sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hl**4*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hr**4*
     &    t2t**(-1) + 128*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)*
     &    t2**(-1)*u2 + 128*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)*
     &    t2*u2**(-1) + 256*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)
     &     + 128*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*t2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 128*
     &    (1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*u2**(-1) - 256*(1+m12/s4)*
     &    pq*lq*ssz*ssp*s**(-1)*sz**(-1) + 128*(1+m12/s4)*pq*rq*ssz*ssp
     &    *mz**(-2)*s**(-1)*t2**(-1)*u2 + 128*(1+m12/s4)*pq*rq*ssz*ssp*
     &    mz**(-2)*s**(-1)*t2*u2**(-1) + 256*(1+m12/s4)*pq*rq*ssz*ssp*
     &    mz**(-2)*s**(-1) + 128*(1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*
     &    t2**(-1) + 128*(1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*u2**(-1) - 
     &    256*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*sz**(-1) + 32*(1+m12/s4)
     &    *pq*hl**2*ssp*m1**2*s**(-1)*t2**(-1)*t2t**(-1) + 16*
     &    (1+m12/s4)*pq*hl**2*ssp*m1**2*s**(-1)*u2**(-1)*t2t**(-1) - 32
     &    *(1+m12/s4)*pq*hl**2*ssp*mt**2*s**(-1)*t2**(-1)*t2t**(-1) - 
     &    16*(1+m12/s4)*pq*hl**2*ssp*mt**2*s**(-1)*u2**(-1)*t2t**(-1)
     &     - 16*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2**(-1)*u2*t2t**(-1)
     &     - 16*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2**(-1) + 16*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**(-1) - 32*(1+m12/s4)*pq*
     &    hl**2*ssp*s**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 32
     &    *(1+m12/s4)*pq*hl**2*ssp*t2**(-1)*t2t**(-1) - 16*(1+m12/s4)*
     &    pq*hl**2*ssp*u2**(-1)*t2t**(-1) + 32*(1+m12/s4)*pq*hr**2*ssp*
     &    m1**2*s**(-1)*t2**(-1)*t2t**(-1) + 16*(1+m12/s4)*pq*hr**2*ssp
     &    *m1**2*s**(-1)*u2**(-1)*t2t**(-1) - 32*(1+m12/s4)*pq*hr**2*
     &    ssp*mt**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*(1+m12/s4)*pq*hr**2
     &    *ssp*mt**2*s**(-1)*u2**(-1)*t2t**(-1) - 16*(1+m12/s4)*pq*
     &    hr**2*ssp*s**(-1)*t2**(-1)*u2*t2t**(-1) - 16*(1+m12/s4)*pq*
     &    hr**2*ssp*s**(-1)*t2**(-1) + 16*(1+m12/s4)*pq*hr**2*ssp*
     &    s**(-1)*u2**(-1) - 32*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*
     &    t2t**(-1) - 32*(1+m12/s4)*pq*hr**2*ssp*t2**(-1)*t2t**(-1) - 
     &    16*(1+m12/s4)*pq*hr**2*ssp*u2**(-1)*t2t**(-1) - 32*(1+m12/s4)
     &    *lq*hl**2*ssz*sz**(-1)*t2t**(-1) - 32*(1+m12/s4)*rq*hr**2*ssz
     &    *sz**(-1)*t2t**(-1) - 128*(1+m12/s4)*ssz**2*lq2*sz**(-2) - 
     &    128*(1+m12/s4)*ssz**2*rq2*sz**(-2) + 128*(1+m12/s4)*ssp**2*
     &    pq2*s**(-2)*t2**(-2)*u2**2 )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 128*
     &    (1+m12/s4)*ssp**2*pq2*s**(-2)*t2**(-1)*u2 + 128*(1+m12/s4)*
     &    ssp**2*pq2*s**(-2)*t2*u2**(-1) + 128*(1+m12/s4)*ssp**2*pq2*
     &    s**(-2)*t2**2*u2**(-2) - 256*(1+m12/s4)*ssp**2*pq2*s**(-2) + 
     &    256*(1+m12/s4)*ssp**2*pq2*s**(-1)*t2**(-2)*u2 + 256*
     &    (1+m12/s4)*ssp**2*pq2*s**(-1)*t2*u2**(-2) + 128*(1+m12/s4)*
     &    ssp**2*pq2*t2**(-2) + 128*(1+m12/s4)*ssp**2*pq2*u2**(-2) + 16
     &    *dyfact(mz)*(s+t2)*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*
     &    hl**2*ssz*s**(-1)*u2**(-1) + 16*dyfact(mz)*(s+t2)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*
     &    u2**(-1) - 128*dyfact(mz)*(s+t2)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    mz**(-2)*s**(-1)*u2**(-1) - 128*dyfact(mz)*(s+t2)*(1+m12/s4)*
     &    pq*rq*ssz*ssp*mz**(-2)*s**(-1)*u2**(-1) + 32*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*
     &    s**(-1)*u2**(-1) - 32*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*
     &    (1+m12/s4)*lq*hl**2*ssz*mt**2*s**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 32*
     &    dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*
     &    m1**2*s**(-1)*u2**(-1) - 32*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*mt**2*
     &    s**(-1)*u2**(-1) - 128*dyfact(mz)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    mz**(-2)*s**(-1) - 128*dyfact(mz)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-1)*u2**(-1) - 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    mz**(-2)*s**(-1) - 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    s**(-1)*u2**(-1) + 16*dyfact(mz)*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*s**(-1)*u2**(-1)*t2t**(-1) - 16*dyfact(mz)*(1+m12/s4)*
     &    lq*hl**2*ssz*mt**2*s**(-1)*u2**(-1)*t2t**(-1) + 16*dyfact(mz)
     &    *(1+m12/s4)*lq*hl**2*ssz*s**(-1)*u2**(-1) - 16*dyfact(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz*u2**(-1)*t2t**(-1) + 16*dyfact(mz)*
     &    (1+m12/s4)*rq*hr**2*ssz*m1**2*s**(-1)*u2**(-1)*t2t**(-1) - 16
     &    *dyfact(mz)*(1+m12/s4)*rq*hr**2*ssz*mt**2*s**(-1)*u2**(-1)*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 16*
     &    dyfact(mz)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*u2**(-1) - 16*
     &    dyfact(mz)*(1+m12/s4)*rq*hr**2*ssz*u2**(-1)*t2t**(-1) - 64*
     &    dyfact(mz)*(1+m12/s4)*ssz**2*lq2*s**(-1)*u2**(-1) - 64*
     &    dyfact(mz)*(1+m12/s4)*ssz**2*rq2*s**(-1)*u2**(-1) + 64*
     &    dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*lq2*mz**2*s**(-2)*
     &    u2**(-2) + 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*lq2*
     &    s**(-2)*u2**(-1) + 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*
     &    rq2*mz**2*s**(-2)*u2**(-2) + 64*dyfact(mz)**2*(s+t2)*
     &    (1+m12/s4)*ssz**2*rq2*s**(-2)*u2**(-1) + 64*dyfact(mz)**2*
     &    (s+t2)**2*(1+m12/s4)*ssz**2*lq2*s**(-2)*u2**(-2) + 64*dyfact(
     &    mz)**2*(s+t2)**2*(1+m12/s4)*ssz**2*rq2*s**(-2)*u2**(-2) - 128
     &    *dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)*
     &    t2**(-1) - 128*dyfacu(mz)*(s+u2)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    mz**(-2)*s**(-1)*t2**(-1) - 16*dyfacu(mz)*(s+u2)*(1+m12/s4)*
     &    lq*hl**2*ssz*s**(-1)*t2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 16
     &    *dyfacu(mz)*(s+u2)*(1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2**(-1)*
     &    t2t**(-1) + 16*dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*
     &    lq*hl**2*ssz*m1**2*mz**2*s**(-1)*t2**(-1)*t2t**(-1) + 16*
     &    dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*s**(-1)*t2**(-1) - 16*dyfacu(mz)*(s+u2-m12+mt2)**(-1)
     &    *(1+m12/s4)*lq*hl**2*ssz*m1**2*t2**(-1)*t2t**(-1) - 16*
     &    dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    mt**2*mz**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*mt**2*s**(-1)*
     &    t2**(-1) + 16*dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq
     &    *hl**2*ssz*mt**2*t2**(-1)*t2t**(-1) + 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*mz**2*
     &    s**(-1)*t2**(-1)*t2t**(-1) + 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*s**(-1)*
     &    t2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 16
     &    *dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*
     &    m1**2*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*mt**2*mz**2*
     &    s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*mt**2*s**(-1)*
     &    t2**(-1) + 16*dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq
     &    *hr**2*ssz*mt**2*t2**(-1)*t2t**(-1) - 128*dyfacu(mz)*
     &    (1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1) - 128*dyfacu(mz)*
     &    (1+m12/s4)*pq*lq*ssz*ssp*s**(-1)*t2**(-1) - 128*dyfacu(mz)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*s**(-1) - 128*dyfacu(mz)*
     &    (1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*t2**(-1) + 32*dyfacu(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz*m1**2*s**(-1)*t2**(-1)*t2t**(-1) - 32
     &    *dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*mt**2*s**(-1)*t2**(-1)*
     &    t2t**(-1) - 16*dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*mz**2*
     &    s**(-1)*t2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 16
     &    *dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*s**(-1)*t2**(-1) - 16*
     &    dyfacu(mz)*(1+m12/s4)*lq*hl**2*ssz*t2**(-1)*t2t**(-1) + 32*
     &    dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*m1**2*s**(-1)*t2**(-1)*
     &    t2t**(-1) - 32*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*mt**2*
     &    s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*(1+m12/s4)*rq*
     &    hr**2*ssz*mz**2*s**(-1)*t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*
     &    (1+m12/s4)*rq*hr**2*ssz*s**(-1)*t2**(-1) - 16*dyfacu(mz)*
     &    (1+m12/s4)*rq*hr**2*ssz*t2**(-1)*t2t**(-1) - 64*dyfacu(mz)*
     &    (1+m12/s4)*ssz**2*lq2*s**(-1)*t2**(-1) - 64*dyfacu(mz)*
     &    (1+m12/s4)*ssz**2*rq2*s**(-1)*t2**(-1) + 64*dyfacu(mz)**2*
     &    (s+u2)*(1+m12/s4)*ssz**2*lq2*mz**2*s**(-2)*t2**(-2) + 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*lq2*s**(-2)*t2**(-1)
     &     + 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*mz**2*
     &    s**(-2)*t2**(-2) + 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*
     &    rq2*s**(-2)*t2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * ( 64*
     &    dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*lq2*s**(-2)*
     &    t2**(-2) + 64*dyfacu(mz)**2*(s+u2)**2*(1+m12/s4)*ssz**2*rq2*
     &    s**(-2)*t2**(-2) + 32*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*
     &    pq*hl**2*ssp*m1**2*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*mt**2*t2**(-1)*
     &    u2**(-1) + 16*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2
     &    *ssp*s*t2**(-1)*u2**(-1) + 16*topfac*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2**(-1)*
     &    u2**(-1) - 32*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2
     &    *ssp*mt**2*t2**(-1)*u2**(-1) + 16*topfac*(u2-m12+mt2)**(-1)
     &    *(1+m12/s4)*pq*hr**2*ssp*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*u2**(-1) + 32*
     &    topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*t2**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 32
     &    *topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    mt**2*t2**(-1)*u2**(-1) + 16*topfac*(u2-m12+mt2+mz2)**(-1)
     &    *(1+m12/s4)*lq*hl**2*ssz*s*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*u2**(-1) + 
     &    32*topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*
     &    m1**2*t2**(-1)*u2**(-1) - 32*topfac*(u2-m12+mt2+mz2)**(-1)
     &    *(1+m12/s4)*rq*hr**2*ssz*mt**2*t2**(-1)*u2**(-1) + 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*s*t2**(-1)*
     &    u2**(-1) + 16*topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*
     &    hr**2*ssz*u2**(-1) + 4*topfac*(1+m12/s4)*hl**4*m1**2*t2**(-1)
     &    *u2**(-1)*t2t**(-1) - 4*topfac*(1+m12/s4)*hl**4*mt**2*
     &    t2**(-1)*u2**(-1)*t2t**(-1) - 4*topfac*(1+m12/s4)*hl**4*s*
     &    t2**(-1)*u2**(-1)*t2t**(-1) + 4*topfac*(1+m12/s4)*hr**4*m1**2
     &    *t2**(-1)*u2**(-1)*t2t**(-1) - 4*topfac*(1+m12/s4)*hr**4*
     &    mt**2*t2**(-1)*u2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4*Pi**2*alphas*hardfac * (  - 4*
     &    topfac*(1+m12/s4)*hr**4*s*t2**(-1)*u2**(-1)*t2t**(-1) + 4*
     &    topfac**2*(1+m12/s4)*hl**4*m1**2*s*t2**(-2)*u2**(-2) + 4*
     &    topfac**2*(1+m12/s4)*hl**4*m1**2*t2**(-1)*u2**(-2) - 4*topfac
     &    **2*(1+m12/s4)*hl**4*mt**2*s*t2**(-2)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hl**4*mt**2*t2**(-1)*u2**(-2) + 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*s*t2**(-2)*u2**(-2) + 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*t2**(-1)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*mt**2*s*t2**(-2)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*mt**2*t2**(-1)*u2**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 4*
     &    (s+u2)**(-2)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hl**4*m1**2*s*
     &    t2t**(-2) - 4*(s+u2)**(-2)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*
     &    hl**4*mt**2*s*t2t**(-2) + 4*(s+u2)**(-2)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*hr**4*m1**2*s*t2t**(-2) - 4
     &    *(s+u2)**(-2)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hr**4*mt**2*s
     &    *t2t**(-2) - 4*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *hl**4*s*t2t**(-2) - 4*(s+u2)**(-1)*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*hr**4*s*t2t**(-2) - 16*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2**(-1)
     &     - 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*s**(-1)*
     &    t2**(-1) + 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*
     &    t2**(-1)*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*
     &    hr**2*ssp*s**(-1)*t2**(-1) + 16*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hr**2*ssp*t2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 
     &     - 128*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)*t2**(-1) - 
     &    128*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*s**(-1)*u2**(-1) - 128*
     &    (1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*s**(-1)*t2**(-1) - 128*
     &    (1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*s**(-1)*u2**(-1) + 16*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*t2**(-1)*t2t**(-1) + 16*
     &    (1+m12/s4)*pq*hl**2*ssp*s**(-1)*u2**(-1)*t2t**(-1) + 16*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*t2**(-1)*t2t**(-1) + 16*
     &    (1+m12/s4)*pq*hr**2*ssp*s**(-1)*u2**(-1)*t2t**(-1) - 128*
     &    (1+m12/s4)*ssp**2*pq2*s**(-2)*t2**(-2)*u2 + 128*(1+m12/s4)*
     &    ssp**2*pq2*s**(-2)*t2**(-1) - 128*(1+m12/s4)*ssp**2*pq2*
     &    s**(-2)*t2*u2**(-2) + 128*(1+m12/s4)*ssp**2*pq2*s**(-2)*
     &    u2**(-1) - 128*(1+m12/s4)*ssp**2*pq2*s**(-1)*t2**(-2) - 128*
     &    (1+m12/s4)*ssp**2*pq2*s**(-1)*u2**(-2) - 16*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*u2**(-1)*
     &    sz**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 
     &     - 16*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*
     &    hr**2*ssz*u2**(-1)*sz**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*lq
     &    *ssz*ssp*mz**(-2)*s**(-1)*u2**(-1) + 128*dyfact(mz)*
     &    (1+m12/s4)*pq*lq*ssz*ssp*s**(-1)*u2**(-1)*sz**(-1) + 128*
     &    dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*mz**(-2)*s**(-1)*u2**(-1)
     &     + 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*s**(-1)*u2**(-1)*
     &    sz**(-1) + 16*dyfact(mz)*(1+m12/s4)*lq*hl**2*ssz*u2**(-1)*
     &    sz**(-1)*t2t**(-1) + 16*dyfact(mz)*(1+m12/s4)*rq*hr**2*ssz*
     &    u2**(-1)*sz**(-1)*t2t**(-1) + 64*dyfact(mz)*(1+m12/s4)*ssz**2
     &    *lq2*u2**(-1)*sz**(-2) + 64*dyfact(mz)*(1+m12/s4)*ssz**2*rq2*
     &    u2**(-1)*sz**(-2) - 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2
     &    *lq2*s**(-1)*u2**(-2)*sz**(-1) - 64*dyfact(mz)**2*(s+t2)*
     &    (1+m12/s4)*ssz**2*rq2*s**(-1)*u2**(-2)*sz**(-1) - 16*dyfacu(
     &    mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*t2**(-1)*
     &    sz**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 16
     &    *dyfacu(mz)*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    t2**(-1)*t2t**(-1) - 16*dyfacu(mz)*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*rq*hr**2*ssz*t2**(-1)*sz**(-1) + 16*dyfacu(mz)*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*t2**(-1)*
     &    t2t**(-1) + 128*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*mz**(-2)*
     &    s**(-1)*t2**(-1) + 128*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    s**(-1)*t2**(-1)*sz**(-1) + 128*dyfacu(mz)*(1+m12/s4)*pq*rq*
     &    ssz*ssp*mz**(-2)*s**(-1)*t2**(-1) + 128*dyfacu(mz)*(1+m12/s4)
     &    *pq*rq*ssz*ssp*s**(-1)*t2**(-1)*sz**(-1) + 16*dyfacu(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz*t2**(-1)*sz**(-1)*t2t**(-1) + 16*
     &    dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*t2**(-1)*sz**(-1)*
     &    t2t**(-1) + 64*dyfacu(mz)*(1+m12/s4)*ssz**2*lq2*t2**(-1)*
     &    sz**(-2) + 64*dyfacu(mz)*(1+m12/s4)*ssz**2*rq2*t2**(-1)*
     &    sz**(-2) - 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*lq2*
     &    s**(-1)*t2**(-2)*sz**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 
     &     - 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*s**(-1)*
     &    t2**(-2)*sz**(-1) - 16*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *pq*hl**2*ssp*t2**(-1)*u2**(-1) - 16*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*t2**(-1)*
     &    u2**(-1) - 16*topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*
     &    hl**2*ssz*s*t2**(-1)*u2**(-1)*sz**(-1) - 16*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*s*t2**(-1)*
     &    u2**(-1)*sz**(-1) + 16*topfac*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*t2**(-1)*u2**(-1) + 16*topfac*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*t2**(-1)*
     &    u2**(-1) + 16*topfac*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*lq*
     &    hl**2*ssz*s*t2**(-1)*u2**(-1)*sz**(-1) + 16*topfac*
     &    (s+u2-m12+mt2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*s*t2**(-1)*
     &    u2**(-1)*sz**(-1) + 4*topfac*(s+u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*hl**4*s*t2**(-1)*u2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*s4**2*Pi**2*alphas*hardfac * ( 4*
     &    topfac*(s+u2-m12+mt2)**(-1)*(1+m12/s4)*hr**4*s*t2**(-1)*
     &    u2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 32*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2t**(-1) + 32*
     &    (s+u2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2t**(-1) + 32*
     &    (s+u2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*t2t**(-1) + 32*
     &    (s+u2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*t2t**(-1) - 32*
     &    (s+u2)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    t2t**(-1) - 32*(s+u2)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2t**(-1) - 32*(u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *pq*hl**2*ssp*m1**2*u2**(-1) - 32*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hr**2*ssp*m1**2*u2**(-1) - 128*(1+m12/s4)*pq*lq
     &    *ssz*ssp*m1**2*mz**(-2)*t2**(-1) - 128*(1+m12/s4)*pq*lq*ssz*
     &    ssp*m1**2*mz**(-2)*u2**(-1) - 128*(1+m12/s4)*pq*rq*ssz*ssp*
     &    m1**2*mz**(-2)*t2**(-1) - 128*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*
     &    mz**(-2)*u2**(-1) + 32*(1+m12/s4)*pq*hl**2*ssp*m1**2*t2**(-1)
     &    *t2t**(-1) + 32*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2**(-1)*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 128*
     &    (1+m12/s4)*ssp**2*pq2*m1**2*s**(-1)*t2**(-2)*u2 - 128*
     &    (1+m12/s4)*ssp**2*pq2*m1**2*s**(-1)*t2**(-1) - 128*(1+m12/s4)
     &    *ssp**2*pq2*m1**2*s**(-1)*t2*u2**(-2) - 128*(1+m12/s4)*ssp**2
     &    *pq2*m1**2*s**(-1)*u2**(-1) - 128*(1+m12/s4)*ssp**2*pq2*m1**2
     &    *t2**(-2) - 128*(1+m12/s4)*ssp**2*pq2*m1**2*u2**(-2) + 32*
     &    dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*mz**2*s**(-1)*u2**(-1) - 32*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*
     &    u2**(-1) + 32*dyfact(mz)*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)
     &    *rq*hr**2*ssz*m1**2*mz**2*s**(-1)*u2**(-1) - 32*dyfact(mz)*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*
     &    u2**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*
     &    mz**(-2)*u2**(-1) - 128*dyfact(mz)*(1+m12/s4)*pq*lq*ssz*ssp*
     &    m1**2*s**(-1)*u2**(-1) + 128*dyfact(mz)*(1+m12/s4)*pq*rq*ssz*
     &    ssp*m1**2*mz**(-2)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 128*
     &    dyfact(mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*s**(-1)*u2**(-1) - 
     &    64*dyfact(mz)*(1+m12/s4)*ssz**2*lq2*m1**2*s**(-1)*u2**(-1) - 
     &    64*dyfact(mz)*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*u2**(-1) + 
     &    64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*lq2*m1**2*mz**2*
     &    s**(-2)*u2**(-2) - 64*dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*
     &    lq2*m1**2*s**(-1)*u2**(-2) + 64*dyfact(mz)**2*(s+t2)*
     &    (1+m12/s4)*ssz**2*rq2*m1**2*mz**2*s**(-2)*u2**(-2) - 64*
     &    dyfact(mz)**2*(s+t2)*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*
     &    u2**(-2) + 32*dyfact(mh1)*(mh12-mh22)**(-1)*(1+m12/s4)*h1*
     &    h2*lambda1*lambda2*mh1**2*s**(-1)*u2**(-1) - 32*dyfact(mh1)*
     &    (mh12-mh22)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*u2**(-1)
     &     - 32*dyfact(mh1)*(u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr
     &    *h1*lambda1*sqrt2**(-1)*mt*mh1**2*s**(-1)*u2**(-1) + 32*
     &    dyfact(mh1)*(u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 16*
     &    dyfact(mh1)*(1+m12/s4)*h1**2*lambda1**2*s**(-1)*u2**(-1) - 16
     &    *dyfact(mh1)**2*(s+t2)*(1+m12/s4)*h1**2*lambda1**2*mh1**2*
     &    s**(-2)*u2**(-2) + 16*dyfact(mh1)**2*(s+t2)*(1+m12/s4)*h1**2*
     &    lambda1**2*s**(-1)*u2**(-2) - 32*dyfact(mh2)*
     &    (mh12-mh22)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*mh2**2*
     &    s**(-1)*u2**(-1) + 32*dyfact(mh2)*(mh12-mh22)**(-1)*
     &    (1+m12/s4)*h1*h2*lambda1*lambda2*u2**(-1) - 32*dyfact(mh2)*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*s**(-1)*u2**(-1) + 32*dyfact(mh2)*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*u2**(-1) + 16*dyfact(mh2)*(1+m12/s4)*h2**2*
     &    lambda2**2*s**(-1)*u2**(-1) - 16*dyfact(mh2)**2*(s+t2)*
     &    (1+m12/s4)*h2**2*lambda2**2*mh2**2*s**(-2)*u2**(-2) + 16*
     &    dyfact(mh2)**2*(s+t2)*(1+m12/s4)*h2**2*lambda2**2*s**(-1)*
     &    u2**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 128*
     &    dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*mz**(-2)*t2**(-1)
     &     - 128*dyfacu(mz)*(1+m12/s4)*pq*lq*ssz*ssp*m1**2*s**(-1)*
     &    t2**(-1) + 128*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*m1**2*
     &    mz**(-2)*t2**(-1) - 128*dyfacu(mz)*(1+m12/s4)*pq*rq*ssz*ssp*
     &    m1**2*s**(-1)*t2**(-1) - 32*dyfacu(mz)*(1+m12/s4)*lq*hl**2*
     &    ssz*m1**2*mz**2*s**(-1)*t2**(-1)*t2t**(-1) + 32*dyfacu(mz)*
     &    (1+m12/s4)*lq*hl**2*ssz*m1**2*t2**(-1)*t2t**(-1) - 32*dyfacu(
     &    mz)*(1+m12/s4)*rq*hr**2*ssz*m1**2*mz**2*s**(-1)*t2**(-1)*
     &    t2t**(-1) + 32*dyfacu(mz)*(1+m12/s4)*rq*hr**2*ssz*m1**2*
     &    t2**(-1)*t2t**(-1) - 64*dyfacu(mz)*(1+m12/s4)*ssz**2*lq2*
     &    m1**2*s**(-1)*t2**(-1) - 64*dyfacu(mz)*(1+m12/s4)*ssz**2*rq2*
     &    m1**2*s**(-1)*t2**(-1) + 64*dyfacu(mz)**2*(s+u2)*(1+m12/s4)*
     &    ssz**2*lq2*m1**2*mz**2*s**(-2)*t2**(-2) - 64*dyfacu(mz)**2*
     &    (s+u2)*(1+m12/s4)*ssz**2*lq2*m1**2*s**(-1)*t2**(-2) + 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*m1**2*mz**2*
     &    s**(-2)*t2**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 64*
     &    dyfacu(mz)**2*(s+u2)*(1+m12/s4)*ssz**2*rq2*m1**2*s**(-1)*
     &    t2**(-2) + 32*dyfacu(mh1)*(mh12-mh22)**(-1)*(1+m12/s4)*h1*
     &    h2*lambda1*lambda2*mh1**2*s**(-1)*t2**(-1) - 32*dyfacu(mh1)*
     &    (mh12-mh22)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*t2**(-1)
     &     + 32*dyfacu(mh1)*(1+m12/s4)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    mh1**2*s**(-1)*t2**(-1)*t2t**(-1) - 32*dyfacu(mh1)*(1+m12/s4)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2**(-1)*t2t**(-1) + 16*
     &    dyfacu(mh1)*(1+m12/s4)*h1**2*lambda1**2*s**(-1)*t2**(-1) - 16
     &    *dyfacu(mh1)**2*(s+u2)*(1+m12/s4)*h1**2*lambda1**2*mh1**2*
     &    s**(-2)*t2**(-2) + 16*dyfacu(mh1)**2*(s+u2)*(1+m12/s4)*h1**2*
     &    lambda1**2*s**(-1)*t2**(-2) - 32*dyfacu(mh2)*
     &    (mh12-mh22)**(-1)*(1+m12/s4)*h1*h2*lambda1*lambda2*mh2**2*
     &    s**(-1)*t2**(-1) + 32*dyfacu(mh2)*(mh12-mh22)**(-1)*
     &    (1+m12/s4)*h1*h2*lambda1*lambda2*t2**(-1) + 32*dyfacu(mh2)*
     &    (1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*s**(-1)*
     &    t2**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 32*
     &    dyfacu(mh2)*(1+m12/s4)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    t2**(-1)*t2t**(-1) + 16*dyfacu(mh2)*(1+m12/s4)*h2**2*
     &    lambda2**2*s**(-1)*t2**(-1) - 16*dyfacu(mh2)**2*(s+u2)*
     &    (1+m12/s4)*h2**2*lambda2**2*mh2**2*s**(-2)*t2**(-2) + 16*
     &    dyfacu(mh2)**2*(s+u2)*(1+m12/s4)*h2**2*lambda2**2*s**(-1)*
     &    t2**(-2) - 32*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2
     &    *ssp*m1**2*mt**2*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hl**2*ssp*m1**2*s*t2**(-1)
     &    *u2**(-1) - 32*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*
     &    hl**2*ssp*m1**2*t2**(-1) + 32*topfac*(u2-m12+mt2)**(-1)*
     &    (1+m12/s4)*pq*hl**2*ssp*m1**4*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*mt**2*
     &    t2**(-1)*u2**(-1) - 32*topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)
     &    *pq*hr**2*ssp*m1**2*s*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**2*t2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 32*
     &    topfac*(u2-m12+mt2)**(-1)*(1+m12/s4)*pq*hr**2*ssp*m1**4*
     &    t2**(-1)*u2**(-1) - 32*topfac*(u2-m12+mt2+mz2)**(-1)*
     &    (1+m12/s4)*lq*hl**2*ssz*m1**2*mt**2*t2**(-1)*u2**(-1) - 32*
     &    topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*
     &    m1**2*s*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*hl**2*ssz*m1**2*
     &    t2**(-1) + 32*topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*lq*
     &    hl**2*ssz*m1**4*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*mt**2
     &    *t2**(-1)*u2**(-1) - 32*topfac*(u2-m12+mt2+mz2)**(-1)*
     &    (1+m12/s4)*rq*hr**2*ssz*m1**2*s*t2**(-1)*u2**(-1) - 32*topfac
     &    *(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*hr**2*ssz*m1**2*
     &    t2**(-1) + 32*topfac*(u2-m12+mt2+mz2)**(-1)*(1+m12/s4)*rq*
     &    hr**2*ssz*m1**4*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 32*
     &    topfac*(u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s*t2**(-1)*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*t2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh12)**(-1)*(1+m12/s4)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*t2**(-1)*u2**(-1) - 32*topfac*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2**(-1)*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2**(-1)*u2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*t2**(-1) + 32*topfac*
     &    (u2-m12+mt2+mh22)**(-1)*(1+m12/s4)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*t2**(-1)*u2**(-1) - 8*topfac*(1+m12/s4)*
     &    hl**2*hr**2*mt**2*t2**(-1)*u2**(-1) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * ( 4*topfac
     &    *(1+m12/s4)*hl**4*m1**2*t2**(-1)*u2**(-1) + 4*topfac*
     &    (1+m12/s4)*hr**4*m1**2*t2**(-1)*u2**(-1) - 8*topfac**2*
     &    (1+m12/s4)*hl**2*hr**2*m1**2*mt**2*s*t2**(-2)*u2**(-2) - 8*
     &    topfac**2*(1+m12/s4)*hl**2*hr**2*m1**2*mt**2*t2**(-1)*
     &    u2**(-2) + 8*topfac**2*(1+m12/s4)*hl**2*hr**2*mt**2*s*
     &    t2**(-2)*u2**(-1) + 8*topfac**2*(1+m12/s4)*hl**2*hr**2*mt**2*
     &    s*t2**(-1)*u2**(-2) + 8*topfac**2*(1+m12/s4)*hl**2*hr**2*
     &    mt**2*s**2*t2**(-2)*u2**(-2) + 8*topfac**2*(1+m12/s4)*hl**2*
     &    hr**2*mt**2*t2**(-1)*u2**(-1) + 8*topfac**2*(1+m12/s4)*hl**2*
     &    hr**2*mt**4*s*t2**(-2)*u2**(-2) + 8*topfac**2*(1+m12/s4)*
     &    hl**2*hr**2*mt**4*t2**(-1)*u2**(-2) - 4*topfac**2*(1+m12/s4)*
     &    hl**4*m1**2*mt**2*s*t2**(-2)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hl**4*m1**2*mt**2*t2**(-1)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hl**4*m1**2*s*t2**(-2)*u2**(-1) - 4*topfac**2*
     &    (1+m12/s4)*hl**4*m1**2*s*t2**(-1)*u2**(-2) )
      MMpartial4 = MMpartial4 + Nc*Cf*Pi**2*alphas*hardfac * (  - 4*
     &    topfac**2*(1+m12/s4)*hl**4*m1**2*s**2*t2**(-2)*u2**(-2) - 4*
     &    topfac**2*(1+m12/s4)*hl**4*m1**2*t2**(-1)*u2**(-1) + 4*topfac
     &    **2*(1+m12/s4)*hl**4*m1**4*s*t2**(-2)*u2**(-2) + 4*topfac**2*
     &    (1+m12/s4)*hl**4*m1**4*t2**(-1)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*mt**2*s*t2**(-2)*u2**(-2) - 4*topfac**
     &    2*(1+m12/s4)*hr**4*m1**2*mt**2*t2**(-1)*u2**(-2) - 4*topfac**
     &    2*(1+m12/s4)*hr**4*m1**2*s*t2**(-2)*u2**(-1) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*s*t2**(-1)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*s**2*t2**(-2)*u2**(-2) - 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**2*t2**(-1)*u2**(-1) + 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**4*s*t2**(-2)*u2**(-2) + 4*topfac**2*
     &    (1+m12/s4)*hr**4*m1**4*t2**(-1)*u2**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(0,0,0,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-2) + 2*hl**4*
     &    t2t**(-2) + 2*hr**4*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(0,0,0,0)*Nc*Cf*Pi*alphas*hardfac
     &  * (  - 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*m1**2*mt**2*
     &    t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**4*
     &    t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-1) + 
     &    2*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-2) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*t2t**(-2) - 2*hl**4*m1**2*
     &    t2t**(-2) + 2*hl**4*mt**2*t2t**(-2) - 2*hl**4*s*t2t**(-2) - 2
     &    *hl**4*u2*t2t**(-2) - 2*hl**4*t2t**(-1) - 2*hr**4*m1**2*
     &    t2t**(-2) + 2*hr**4*mt**2*t2t**(-2) - 2*hr**4*s*t2t**(-2) - 2
     &    *hr**4*u2*t2t**(-2) - 2*hr**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*s + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*s + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s*t2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s**2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s*t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s**2*sz**(-1)
     &     + 2*(s+u2-m12+mt2)**(-1)*hl**4*s + 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*s**2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*s + 2*(s+u2-m12+mt2)**(-1)*
     &    hr**4*s**2*t2t**(-1) - 128*pq*lq*ssz*ssp*s*sz**(-1) - 64*pq*
     &    lq*ssz*ssp*t2*sz**(-1) - 64*pq*lq*ssz*ssp*u2*sz**(-1) - 128*
     &    pq*rq*ssz*ssp*s*sz**(-1) - 64*pq*rq*ssz*ssp*t2*sz**(-1) - 64*
     &    pq*rq*ssz*ssp*u2*sz**(-1) + 8*pq*hl**2*ssp*m1**2*t2t**(-1) - 
     &    8*pq*hl**2*ssp*mt**2*t2t**(-1) - 16*pq*hl**2*ssp*s*t2t**(-1)
     &     - 8*pq*hl**2*ssp*u2*t2t**(-1) + 8*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) - 8*pq*hr**2*ssp*mt**2*t2t**(-1) - 16*pq*hr**2*ssp*
     &    s*t2t**(-1) - 8*pq*hr**2*ssp*u2*t2t**(-1) + 8*lq*hl**2*ssz*
     &    m1**2*s*sz**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*mt**2*s*sz**(-1)*
     &    t2t**(-1) - 8*lq*hl**2*ssz*s*u2*sz**(-1)*t2t**(-1) - 16*lq*
     &    hl**2*ssz*s**2*sz**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*rq*hr**2*ssz*m1**2*s*sz**(-1)*t2t**(-1) - 8*rq*
     &    hr**2*ssz*mt**2*s*sz**(-1)*t2t**(-1) - 8*rq*hr**2*ssz*s*u2*
     &    sz**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*s**2*sz**(-1)*t2t**(-1)
     &     + 2*hl**4*s*t2t**(-1) + 2*hr**4*s*t2t**(-1) - 32*ssz**2*lq2*
     &    s*t2*sz**(-2) - 32*ssz**2*lq2*s*u2*sz**(-2) - 64*ssz**2*lq2*
     &    s**2*sz**(-2) - 32*ssz**2*rq2*s*t2*sz**(-2) - 32*ssz**2*rq2*s
     &    *u2*sz**(-2) - 64*ssz**2*rq2*s**2*sz**(-2) - 64*ssp**2*pq2*
     &    s**(-1)*t2 - 64*ssp**2*pq2*s**(-1)*u2 - 128*ssp**2*pq2 )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2
     &     + 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s**2*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s**2*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*
     &    s1**(-1) - 16*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**2*s2**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**2*s**2*
     &    t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*s*
     &    t2t**(-1) - 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**4*s*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*s**2*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*s*t2t**(-1) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*s*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s**2*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*s*t2t**(-1) - 128*pq*lq*
     &    ssz*ssp*m1**2*s*sz**(-1) + 128*pq*lq*ssz*ssp*s*t2*sz**(-1) + 
     &    128*pq*lq*ssz*ssp*s*u2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 128*pq*lq*ssz*ssp*s**2*sz**(-1) + 128*pq*lq*ssz*ssp*
     &    t2*u2*sz**(-1) - 128*pq*rq*ssz*ssp*m1**2*s*sz**(-1) + 128*pq*
     &    rq*ssz*ssp*s*t2*sz**(-1) + 128*pq*rq*ssz*ssp*s*u2*sz**(-1) + 
     &    128*pq*rq*ssz*ssp*s**2*sz**(-1) + 128*pq*rq*ssz*ssp*t2*u2*
     &    sz**(-1) - 32*pq*hl**2*ssp*m1**2*s*t2t**(-1) - 16*pq*hl**2*
     &    ssp*m1**2*u2*t2t**(-1) + 16*pq*hl**2*ssp*mt**2*s*t2t**(-1) + 
     &    16*pq*hl**2*ssp*mt**2*u2*t2t**(-1) + 16*pq*hl**2*ssp*s*u2*
     &    t2t**(-1) + 16*pq*hl**2*ssp*s**2*t2t**(-1) - 16*pq*hl**2*ssp*
     &    t2 + 16*pq*hl**2*ssp*u2 - 32*pq*hr**2*ssp*m1**2*s*t2t**(-1)
     &     - 16*pq*hr**2*ssp*m1**2*u2*t2t**(-1) + 16*pq*hr**2*ssp*mt**2
     &    *s*t2t**(-1) + 16*pq*hr**2*ssp*mt**2*u2*t2t**(-1) + 16*pq*
     &    hr**2*ssp*s*u2*t2t**(-1) + 16*pq*hr**2*ssp*s**2*t2t**(-1) - 
     &    16*pq*hr**2*ssp*t2 + 16*pq*hr**2*ssp*u2 - 16*lq*hl**2*ssz*
     &    m1**2*s*u2*sz**(-1)*t2t**(-1) - 32*lq*hl**2*ssz*m1**2*s**2*
     &    sz**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*lq*hl**2*ssz*mt**2*s*u2*sz**(-1)*t2t**(-1) + 16*
     &    lq*hl**2*ssz*mt**2*s**2*sz**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*
     &    s*t2*sz**(-1) + 16*lq*hl**2*ssz*s*u2*sz**(-1) + 16*lq*hl**2*
     &    ssz*s**2*u2*sz**(-1)*t2t**(-1) + 16*lq*hl**2*ssz*s**3*
     &    sz**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*m1**2*s*u2*sz**(-1)*
     &    t2t**(-1) - 32*rq*hr**2*ssz*m1**2*s**2*sz**(-1)*t2t**(-1) + 
     &    16*rq*hr**2*ssz*mt**2*s*u2*sz**(-1)*t2t**(-1) + 16*rq*hr**2*
     &    ssz*mt**2*s**2*sz**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*s*t2*
     &    sz**(-1) + 16*rq*hr**2*ssz*s*u2*sz**(-1) + 16*rq*hr**2*ssz*
     &    s**2*u2*sz**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*s**3*sz**(-1)*
     &    t2t**(-1) + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*s1**(-1)*
     &    t2t**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*s2**(-1)*
     &    t2t**(-1) + 4*hl**4*m1**2*s*t2t**(-1) - 4*hl**4*mt**2*s*
     &    t2t**(-1) - 4*hl**4*s - 4*hl**4*s**2*t2t**(-1) + 4*hr**4*
     &    m1**2*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hr**4*mt**2*s*t2t**(-1) - 4*hr**4*s - 4*hr**4*
     &    s**2*t2t**(-1) + 32*h1*h2*lambda1*lambda2*s**2*s1**(-1)*
     &    s2**(-1) + 16*h1**2*lambda1**2*s**2*s1**(-2) + 16*h2**2*
     &    lambda2**2*s**2*s2**(-2) - 64*ssz**2*lq2*m1**2*s**2*sz**(-2)
     &     + 64*ssz**2*lq2*s*t2*u2*sz**(-2) + 64*ssz**2*lq2*s**2*t2*
     &    sz**(-2) + 64*ssz**2*lq2*s**2*u2*sz**(-2) + 64*ssz**2*lq2*
     &    s**3*sz**(-2) - 64*ssz**2*rq2*m1**2*s**2*sz**(-2) + 64*ssz**2
     &    *rq2*s*t2*u2*sz**(-2) + 64*ssz**2*rq2*s**2*t2*sz**(-2) + 64*
     &    ssz**2*rq2*s**2*u2*sz**(-2) + 64*ssz**2*rq2*s**3*sz**(-2) - 
     &    128*ssp**2*pq2*m1**2 + 128*ssp**2*pq2*s**(-1)*t2*u2 + 128*
     &    ssp**2*pq2*s + 128*ssp**2*pq2*t2 + 128*ssp**2*pq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(1,2,1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-2) + 2*hl**4*
     &    t2t**(-2) + 2*hr**4*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(1,2,1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*m1**2*mt**2*
     &    t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**4*
     &    t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-1) + 
     &    2*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-2) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*t2t**(-2) - 2*hl**4*m1**2*
     &    t2t**(-2) + 2*hl**4*mt**2*t2t**(-2) - 2*hl**4*s*t2t**(-2) - 2
     &    *hl**4*u2*t2t**(-2) - 2*hl**4*t2t**(-1) - 2*hr**4*m1**2*
     &    t2t**(-2) + 2*hr**4*mt**2*t2t**(-2) - 2*hr**4*s*t2t**(-2) - 2
     &    *hr**4*u2*t2t**(-2) - 2*hr**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 4*hl**2*hr**2*mt**2*s*t2t**(-1) + 8*hl**4*m1**2*
     &    mt**2*t2t**(-1) + 2*hl**4*m1**2 - 4*hl**4*m1**4*t2t**(-1) - 2
     &    *hl**4*mt**2*s*t2t**(-1) - 2*hl**4*mt**2 - 4*hl**4*mt**4*
     &    t2t**(-1) + 8*hr**4*m1**2*mt**2*t2t**(-1) + 2*hr**4*m1**2 - 4
     &    *hr**4*m1**4*t2t**(-1) - 2*hr**4*mt**2*s*t2t**(-1) - 2*hr**4*
     &    mt**2 - 4*hr**4*mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*hl**2*hr**2*m1**2*mt**2*s*t2t**(-1) - 4*hl**2*
     &    hr**2*m1**2*mt**2*u2*t2t**(-1) - 8*hl**2*hr**2*m1**2*mt**4*
     &    t2t**(-1) + 4*hl**2*hr**2*m1**4*mt**2*t2t**(-1) - 4*hl**2*
     &    hr**2*mt**2*s*u2*t2t**(-1) - 8*hl**2*hr**2*mt**2*s - 4*hl**2*
     &    hr**2*mt**2*s**2*t2t**(-1) - 8*hl**2*hr**2*mt**4*s*t2t**(-1)
     &     + 4*hl**2*hr**2*mt**4*u2*t2t**(-1) + 4*hl**2*hr**2*mt**6*
     &    t2t**(-1) - 12*hl**4*m1**2*mt**2*s*t2t**(-1) - 6*hl**4*m1**2*
     &    mt**2*u2*t2t**(-1) - 12*hl**4*m1**2*mt**2 - 20*hl**4*m1**2*
     &    mt**4*t2t**(-1) - 2*hl**4*m1**2*s - 4*hl**4*m1**2*t2 - 2*
     &    hl**4*m1**2*u2 + 22*hl**4*m1**4*mt**2*t2t**(-1) + 4*hl**4*
     &    m1**4*s*t2t**(-1) + 4*hl**4*m1**4*u2*t2t**(-1) + 6*hl**4*
     &    m1**4 - 8*hl**4*m1**6*t2t**(-1) + 2*hl**4*mt**2*s*u2*
     &    t2t**(-1) + 6*hl**4*mt**2*s + 2*hl**4*mt**2*s**2*t2t**(-1) + 
     &    4*hl**4*mt**2*t2 + 2*hl**4*mt**2*u2 + 8*hl**4*mt**4*s*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 2*hl**4*mt**4*u2*t2t**(-1) + 6*hl**4*mt**4 + 6*hl**4
     &    *mt**6*t2t**(-1) - 12*hr**4*m1**2*mt**2*s*t2t**(-1) - 6*hr**4
     &    *m1**2*mt**2*u2*t2t**(-1) - 12*hr**4*m1**2*mt**2 - 20*hr**4*
     &    m1**2*mt**4*t2t**(-1) - 2*hr**4*m1**2*s - 4*hr**4*m1**2*t2 - 
     &    2*hr**4*m1**2*u2 + 22*hr**4*m1**4*mt**2*t2t**(-1) + 4*hr**4*
     &    m1**4*s*t2t**(-1) + 4*hr**4*m1**4*u2*t2t**(-1) + 6*hr**4*
     &    m1**4 - 8*hr**4*m1**6*t2t**(-1) + 2*hr**4*mt**2*s*u2*
     &    t2t**(-1) + 6*hr**4*mt**2*s + 2*hr**4*mt**2*s**2*t2t**(-1) + 
     &    4*hr**4*mt**2*t2 + 2*hr**4*mt**2*u2 + 8*hr**4*mt**4*s*
     &    t2t**(-1) + 2*hr**4*mt**4*u2*t2t**(-1) + 6*hr**4*mt**4 + 6*
     &    hr**4*mt**6*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s - 8*(u2-m12+mt2)**(-1)*
     &    pq*hl**2*ssp*t2 - 32*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2
     &     + 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s - 8*(u2-m12+mt2)**(-1)*
     &    pq*hr**2*ssp*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**4*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s
     &    *t2*sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2*
     &    sz**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2*
     &    sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2*
     &    sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*s + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*s + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s*t2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s**2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*
     &    sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s*t2*
     &    sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s**2*
     &    sz**(-1) + 2*(s+u2-m12+mt2)**(-1)*hl**4*s + 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*s**2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*s + 2*(s+u2-m12+mt2)**(-1)*
     &    hr**4*s**2*t2t**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*
     &    t2t**(-1) + 4*hl**4*m1**2*t2t**(-1) - 4*hl**4*mt**2*t2t**(-1)
     &     - 2*hl**4*s*t2t**(-1) - 2*hl**4 + 4*hr**4*m1**2*t2t**(-1) - 
     &    4*hr**4*mt**2*t2t**(-1) - 2*hr**4*s*t2t**(-1) - 2*hr**4 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 48*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*s*
     &    t2t**(-1) + 24*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    u2*t2t**(-1) + 40*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**2 + 80*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**4*
     &    t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*s + 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 - 88*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*s*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*u2*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**6*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*u2*t2t**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s**2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2 - 32
     &    *(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*s*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*u2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4 - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**6*t2t**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s*t2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2**2 + 48*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*s*t2t**(-1) + 
     &    24*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*u2*t2t**(-1)
     &     + 40*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2 + 80*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*s + 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 - 88*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*s*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*u2*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**6*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*u2*t2t**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s**2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2 - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*s*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*u2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4 - 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**6*t2t**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s*t2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s**2 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2**2 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*s*
     &    t2t**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*u2*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2 + 80*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**4*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2*s*t2*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*t2 - 88*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    mt**2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *m1**4*u2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s
     &    *t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*
     &    u2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**2*sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*mt**2*t2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2
     &    *u2 - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s*
     &    t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*u2
     &    *t2t**(-1) - 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4
     &     - 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**6*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2**2 + 48*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*u2*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*mt**2 + 80*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**4*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**2*s*t2*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*t2 - 88*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    mt**2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *m1**4*u2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4 + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**6*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s
     &    *t2*sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &    *s*u2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**2*sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**2*t2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &    *u2 - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s*
     &    t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*u2
     &    *t2t**(-1) - 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4
     &     - 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**6*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s**2*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s**2*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*
     &    s1**(-1) - 16*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**2*s2**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl**2
     &    *hr**2*mt**2*s**2*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl**4*
     &    m1**2*mt**2*s*t2t**(-1) - 4*(s+u2-m12+mt2)**(-1)*hl**4*
     &    m1**2*s + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**4*s*t2t**(-1) + 
     &    4*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*s )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*s**2*t2t**(-1)
     &     + 4*(s+u2-m12+mt2)**(-1)*hl**4*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*s*t2t**(-1) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*s*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s**2*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*s*t2t**(-1) - 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*s*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt**3*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4
     &    *mt*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*u2*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**2*s1**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s**2*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2 + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*u2*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt**3 - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**5
     &    *t2t**(-1) - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2
     &    *lambda2*sqrt2**(-1)*m1**2*mt**3*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4
     &    *mt*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s*u2*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    s**2*s2**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s**2*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*u2 + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*u2*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt**3 )
      MMpartial4 = MMpartial4 + ANGfin(1,4,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**5*t2t**(-1) - 4*hl**2*hr**2*m1**2*mt**2*
     &    t2t**(-1) + 4*hl**2*hr**2*mt**2*s*t2t**(-1) + 4*hl**2*hr**2*
     &    mt**4*t2t**(-1) - 10*hl**4*m1**2*mt**2*t2t**(-1) - 6*hl**4*
     &    m1**2*s*t2t**(-1) - 2*hl**4*m1**2*u2*t2t**(-1) - 4*hl**4*
     &    m1**2 + 6*hl**4*m1**4*t2t**(-1) + 4*hl**4*mt**2*s*t2t**(-1)
     &     + 2*hl**4*mt**2*u2*t2t**(-1) + 4*hl**4*mt**2 + 4*hl**4*mt**4
     &    *t2t**(-1) + 2*hl**4*s*u2*t2t**(-1) + 2*hl**4*s + 2*hl**4*t2
     &     + 2*hl**4*u2 - 10*hr**4*m1**2*mt**2*t2t**(-1) - 6*hr**4*
     &    m1**2*s*t2t**(-1) - 2*hr**4*m1**2*u2*t2t**(-1) - 4*hr**4*
     &    m1**2 + 6*hr**4*m1**4*t2t**(-1) + 4*hr**4*mt**2*s*t2t**(-1)
     &     + 2*hr**4*mt**2*u2*t2t**(-1) + 4*hr**4*mt**2 + 4*hr**4*mt**4
     &    *t2t**(-1) + 2*hr**4*s*u2*t2t**(-1) + 2*hr**4*s + 2*hr**4*t2
     &     + 2*hr**4*u2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*ssp**2*pq2*s + 64*ssp**2*pq2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*ssp**2*pq2*m1**2*s - 128*ssp**2*pq2*s*t2 - 64*
     &    ssp**2*pq2*s*u2 - 64*ssp**2*pq2*s**2 - 64*ssp**2*pq2*t2*u2 - 
     &    64*ssp**2*pq2*t2**2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s - 8*(u2-m12+mt2)**(-1)*
     &    pq*hl**2*ssp*t2 - 32*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2
     &     + 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s - 8*(u2-m12+mt2)**(-1)*
     &    pq*hr**2*ssp*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*mz**(-2)*s - 64*pq*lq*ssz*ssp*
     &    mz**(-2)*t2 - 64*pq*rq*ssz*ssp*mz**(-2)*s - 64*pq*rq*ssz*ssp*
     &    mz**(-2)*t2 + 8*pq*hl**2*ssp*m1**2*t2t**(-1) - 8*pq*hl**2*ssp
     &    *mt**2*t2t**(-1) + 8*pq*hl**2*ssp*u2*t2t**(-1) + 8*pq*hr**2*
     &    ssp*m1**2*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*t2t**(-1) + 8*pq*
     &    hr**2*ssp*u2*t2t**(-1) + 64*ssp**2*pq2*s**(-1)*t2 + 64*ssp**2
     &    *pq2*s**(-1)*u2 + 64*ssp**2*pq2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 40*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*s*
     &    t2t**(-1) + 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2
     &     + 48*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**4*t2t**(-1)
     &     + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*s + 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 - 48*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*s*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**6*t2t**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s**2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*s*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**6*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*s*t2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2**2 + 40*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*s*t2t**(-1) + 
     &    32*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2 + 48*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**4*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*s + 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 - 48*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*mt**2*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*s*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**6*t2t**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s**2*t2t**(-1) - 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*s*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**6*t2t**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s*t2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s**2 + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2**2 - 64*pq*lq*ssz*ssp*
     &    m1**2*mz**(-2)*s + 128*pq*lq*ssz*ssp*mz**(-2)*s*t2 + 64*pq*lq
     &    *ssz*ssp*mz**(-2)*s*u2 + 64*pq*lq*ssz*ssp*mz**(-2)*s**2 + 64*
     &    pq*lq*ssz*ssp*mz**(-2)*t2*u2 + 64*pq*lq*ssz*ssp*mz**(-2)*
     &    t2**2 - 64*pq*rq*ssz*ssp*m1**2*mz**(-2)*s + 128*pq*rq*ssz*ssp
     &    *mz**(-2)*s*t2 + 64*pq*rq*ssz*ssp*mz**(-2)*s*u2 + 64*pq*rq*
     &    ssz*ssp*mz**(-2)*s**2 + 64*pq*rq*ssz*ssp*mz**(-2)*t2*u2 + 64*
     &    pq*rq*ssz*ssp*mz**(-2)*t2**2 - 16*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 8*pq*hl**2*ssp*m1**2*s*t2t**(-1) + 8*pq*hl**2*ssp
     &    *m1**2*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*pq*hl**2*ssp*m1**2 + 8*pq*hl**2*ssp*m1**4*
     &    t2t**(-1) + 8*pq*hl**2*ssp*mt**2*s*t2t**(-1) - 8*pq*hl**2*ssp
     &    *mt**2*u2*t2t**(-1) + 8*pq*hl**2*ssp*mt**2 + 8*pq*hl**2*ssp*
     &    mt**4*t2t**(-1) - 8*pq*hl**2*ssp*s*u2*t2t**(-1) + 16*pq*hl**2
     &    *ssp*s + 16*pq*hl**2*ssp*t2 - 8*pq*hl**2*ssp*u2 - 16*pq*hr**2
     &    *ssp*m1**2*mt**2*t2t**(-1) - 8*pq*hr**2*ssp*m1**2*s*t2t**(-1)
     &     + 8*pq*hr**2*ssp*m1**2*u2*t2t**(-1) - 8*pq*hr**2*ssp*m1**2
     &     + 8*pq*hr**2*ssp*m1**4*t2t**(-1) + 8*pq*hr**2*ssp*mt**2*s*
     &    t2t**(-1) - 8*pq*hr**2*ssp*mt**2*u2*t2t**(-1) + 8*pq*hr**2*
     &    ssp*mt**2 + 8*pq*hr**2*ssp*mt**4*t2t**(-1) - 8*pq*hr**2*ssp*s
     &    *u2*t2t**(-1) + 16*pq*hr**2*ssp*s + 16*pq*hr**2*ssp*t2 - 8*pq
     &    *hr**2*ssp*u2 + 64*ssp**2*pq2*m1**2 - 128*ssp**2*pq2*s**(-1)*
     &    t2*u2 - 64*ssp**2*pq2*s - 64*ssp**2*pq2*t2 - 128*ssp**2*pq2*
     &    u2 )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*ssz**2*lq2*mz**2*s*sz**(-1) - 32*ssz**2*lq2*mz**2
     &     + 32*ssz**2*lq2*s*t2*sz**(-1) + 32*ssz**2*lq2*s*u2*sz**(-1)
     &     + 32*ssz**2*lq2*s**2*sz**(-1) - 32*ssz**2*lq2*u2 + 32*ssz**2
     &    *rq2*mz**2*s*sz**(-1) - 32*ssz**2*rq2*mz**2 + 32*ssz**2*rq2*s
     &    *t2*sz**(-1) + 32*ssz**2*rq2*s*u2*sz**(-1) + 32*ssz**2*rq2*
     &    s**2*sz**(-1) - 32*ssz**2*rq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*ssz**2*lq2*m1**2*mz**2 - 32*ssz**2*lq2*m1**2*s
     &     + 64*ssz**2*lq2*m1**2*s**2*sz**(-1) - 64*ssz**2*lq2*mz**2*s*
     &    t2*sz**(-1) + 32*ssz**2*lq2*mz**2*s - 64*ssz**2*lq2*mz**2*
     &    s**2*sz**(-1) + 32*ssz**2*lq2*mz**2*t2 - 64*ssz**2*lq2*s*t2*
     &    u2*sz**(-1) - 64*ssz**2*lq2*s*t2 + 32*ssz**2*lq2*s*u2 - 64*
     &    ssz**2*lq2*s**2*u2*sz**(-1) - 32*ssz**2*lq2*s**2 + 32*ssz**2*
     &    lq2*t2*u2 - 32*ssz**2*lq2*t2**2 - 32*ssz**2*rq2*m1**2*mz**2
     &     - 32*ssz**2*rq2*m1**2*s + 64*ssz**2*rq2*m1**2*s**2*sz**(-1)
     &     - 64*ssz**2*rq2*mz**2*s*t2*sz**(-1) + 32*ssz**2*rq2*mz**2*s
     &     - 64*ssz**2*rq2*mz**2*s**2*sz**(-1) + 32*ssz**2*rq2*mz**2*t2
     &     - 64*ssz**2*rq2*s*t2*u2*sz**(-1) - 64*ssz**2*rq2*s*t2 + 32*
     &    ssz**2*rq2*s*u2 - 64*ssz**2*rq2*s**2*u2*sz**(-1) - 32*ssz**2*
     &    rq2*s**2 + 32*ssz**2*rq2*t2*u2 - 32*ssz**2*rq2*t2**2 )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*s*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**4*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s
     &    *t2*sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2*
     &    sz**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2*
     &    sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2*
     &    sz**(-1) + 64*pq*lq*ssz*ssp*mz**(-2)*s + 64*pq*lq*ssz*ssp*
     &    mz**(-2)*t2 + 64*pq*lq*ssz*ssp*mz**2*sz**(-1) + 64*pq*lq*ssz*
     &    ssp*s*sz**(-1) + 64*pq*lq*ssz*ssp*t2*sz**(-1) + 64*pq*lq*ssz*
     &    ssp*u2*sz**(-1) + 64*pq*rq*ssz*ssp*mz**(-2)*s + 64*pq*rq*ssz*
     &    ssp*mz**(-2)*t2 + 64*pq*rq*ssz*ssp*mz**2*sz**(-1) + 64*pq*rq*
     &    ssz*ssp*s*sz**(-1) + 64*pq*rq*ssz*ssp*t2*sz**(-1) + 64*pq*rq*
     &    ssz*ssp*u2*sz**(-1) - 8*lq*hl**2*ssz*m1**2*s*sz**(-1)*
     &    t2t**(-1) + 16*lq*hl**2*ssz*m1**2*t2t**(-1) + 8*lq*hl**2*ssz*
     &    mt**2*s*sz**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*t2t**(-1)
     &     + 8*lq*hl**2*ssz*mz**2*s*sz**(-1)*t2t**(-1) + 8*lq*hl**2*ssz
     &    *s*u2*sz**(-1)*t2t**(-1) - 8*lq*hl**2*ssz*s*t2t**(-1) + 8*lq*
     &    hl**2*ssz*s**2*sz**(-1)*t2t**(-1) - 8*rq*hr**2*ssz*m1**2*s*
     &    sz**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*rq*hr**2*ssz*m1**2*t2t**(-1) + 8*rq*hr**2*ssz*
     &    mt**2*s*sz**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*mt**2*t2t**(-1)
     &     + 8*rq*hr**2*ssz*mz**2*s*sz**(-1)*t2t**(-1) + 8*rq*hr**2*ssz
     &    *s*u2*sz**(-1)*t2t**(-1) - 8*rq*hr**2*ssz*s*t2t**(-1) + 8*rq*
     &    hr**2*ssz*s**2*sz**(-1)*t2t**(-1) + 32*ssz**2*lq2*mz**2*s*
     &    sz**(-2) + 32*ssz**2*lq2*s*t2*sz**(-2) + 32*ssz**2*lq2*s*u2*
     &    sz**(-2) + 32*ssz**2*lq2*s*sz**(-1) + 32*ssz**2*lq2*s**2*
     &    sz**(-2) - 32*ssz**2*lq2 + 32*ssz**2*rq2*mz**2*s*sz**(-2) + 
     &    32*ssz**2*rq2*s*t2*sz**(-2) + 32*ssz**2*rq2*s*u2*sz**(-2) + 
     &    32*ssz**2*rq2*s*sz**(-1) + 32*ssz**2*rq2*s**2*sz**(-2) - 32*
     &    ssz**2*rq2 )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*mz**2*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2*mt**2*s*t2t**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*m1**2*mt**2 + 48*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**4*t2t**(-1) + 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s*t2*sz**(-1) + 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 48*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mz**2*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*s
     &    *t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**6*t2t**(-1)
     &     + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*s*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*t2*
     &    sz**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2*
     &    sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*t2 + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    mz**2*t2t**(-1) - 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**4*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**4 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**6*
     &    t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2 + 
     &    8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2**2 - 24*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*mz**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*s*t2t**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *m1**2*mt**2 + 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**4*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**2*s*t2*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*t2 - 48*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    mt**2*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*mz**2*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**4*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**6*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*mz**2*s*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*mt**2*s*t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2*
     &    sz**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*t2 + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*
     &    mz**2*t2t**(-1) - 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**4*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**4 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**6*
     &    t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2 + 
     &    8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2 + 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2**2 + 64*pq*lq*ssz*
     &    ssp*m1**2*mz**(-2)*s + 128*pq*lq*ssz*ssp*m1**2*s*sz**(-1) - 
     &    64*pq*lq*ssz*ssp*m1**2 - 128*pq*lq*ssz*ssp*mz**(-2)*s*t2 - 64
     &    *pq*lq*ssz*ssp*mz**(-2)*s*u2 - 64*pq*lq*ssz*ssp*mz**(-2)*s**2
     &     - 64*pq*lq*ssz*ssp*mz**(-2)*t2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*mz**(-2)*t2**2 - 128*pq*lq*ssz*
     &    ssp*mz**2*s*sz**(-1) - 128*pq*lq*ssz*ssp*mz**2*t2*sz**(-1) - 
     &    128*pq*lq*ssz*ssp*s*u2*sz**(-1) - 64*pq*lq*ssz*ssp*s - 128*pq
     &    *lq*ssz*ssp*t2*u2*sz**(-1) - 64*pq*lq*ssz*ssp*t2 + 64*pq*rq*
     &    ssz*ssp*m1**2*mz**(-2)*s + 128*pq*rq*ssz*ssp*m1**2*s*sz**(-1)
     &     - 64*pq*rq*ssz*ssp*m1**2 - 128*pq*rq*ssz*ssp*mz**(-2)*s*t2
     &     - 64*pq*rq*ssz*ssp*mz**(-2)*s*u2 - 64*pq*rq*ssz*ssp*mz**(-2)
     &    *s**2 - 64*pq*rq*ssz*ssp*mz**(-2)*t2*u2 - 64*pq*rq*ssz*ssp*
     &    mz**(-2)*t2**2 - 128*pq*rq*ssz*ssp*mz**2*s*sz**(-1) - 128*pq*
     &    rq*ssz*ssp*mz**2*t2*sz**(-1) - 128*pq*rq*ssz*ssp*s*u2*
     &    sz**(-1) - 64*pq*rq*ssz*ssp*s - 128*pq*rq*ssz*ssp*t2*u2*
     &    sz**(-1) - 64*pq*rq*ssz*ssp*t2 - 16*lq*hl**2*ssz*m1**2*mt**2*
     &    t2t**(-1) + 16*lq*hl**2*ssz*m1**2*mz**2*s*sz**(-1)*t2t**(-1)
     &     + 16*lq*hl**2*ssz*m1**2*s*u2*sz**(-1)*t2t**(-1) - 24*lq*
     &    hl**2*ssz*m1**2*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*lq*hl**2*ssz*m1**2*s**2*sz**(-1)*t2t**(-1) - 8*lq
     &    *hl**2*ssz*m1**2*u2*t2t**(-1) - 8*lq*hl**2*ssz*m1**2 + 8*lq*
     &    hl**2*ssz*m1**4*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*mz**2*s*
     &    sz**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*s*u2*sz**(-1)*
     &    t2t**(-1) + 8*lq*hl**2*ssz*mt**2*s*t2t**(-1) + 8*lq*hl**2*ssz
     &    *mt**2*u2*t2t**(-1) + 8*lq*hl**2*ssz*mt**2 + 8*lq*hl**2*ssz*
     &    mt**4*t2t**(-1) - 16*lq*hl**2*ssz*mz**2*s*sz**(-1) - 16*lq*
     &    hl**2*ssz*mz**2*s**2*sz**(-1)*t2t**(-1) + 16*lq*hl**2*ssz*s*
     &    t2*sz**(-1) - 16*lq*hl**2*ssz*s*u2*sz**(-1) + 8*lq*hl**2*ssz*
     &    s*u2*t2t**(-1) - 16*lq*hl**2*ssz*s**2*u2*sz**(-1)*t2t**(-1)
     &     + 16*lq*hl**2*ssz*s**2*sz**(-1) + 8*lq*hl**2*ssz*u2 - 16*rq*
     &    hr**2*ssz*m1**2*mt**2*t2t**(-1) + 16*rq*hr**2*ssz*m1**2*mz**2
     &    *s*sz**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*m1**2*s*u2*sz**(-1)*
     &    t2t**(-1) - 24*rq*hr**2*ssz*m1**2*s*t2t**(-1) + 16*rq*hr**2*
     &    ssz*m1**2*s**2*sz**(-1)*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*rq*hr**2*ssz*m1**2*u2*t2t**(-1) - 8*rq*hr**2*
     &    ssz*m1**2 + 8*rq*hr**2*ssz*m1**4*t2t**(-1) - 16*rq*hr**2*ssz*
     &    mt**2*mz**2*s*sz**(-1)*t2t**(-1) - 16*rq*hr**2*ssz*mt**2*s*u2
     &    *sz**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*mt**2*s*t2t**(-1) + 8*rq
     &    *hr**2*ssz*mt**2*u2*t2t**(-1) + 8*rq*hr**2*ssz*mt**2 + 8*rq*
     &    hr**2*ssz*mt**4*t2t**(-1) - 16*rq*hr**2*ssz*mz**2*s*sz**(-1)
     &     - 16*rq*hr**2*ssz*mz**2*s**2*sz**(-1)*t2t**(-1) + 16*rq*
     &    hr**2*ssz*s*t2*sz**(-1) - 16*rq*hr**2*ssz*s*u2*sz**(-1) + 8*
     &    rq*hr**2*ssz*s*u2*t2t**(-1) - 16*rq*hr**2*ssz*s**2*u2*
     &    sz**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*s**2*sz**(-1) + 8*rq*
     &    hr**2*ssz*u2 + 64*ssz**2*lq2*m1**2*s**2*sz**(-2) - 32*ssz**2*
     &    lq2*m1**2 - 64*ssz**2*lq2*mz**2*s*t2*sz**(-2) - 64*ssz**2*lq2
     &    *mz**2*s**2*sz**(-2) - 64*ssz**2*lq2*s*t2*u2*sz**(-2) - 64*
     &    ssz**2*lq2*s*t2*sz**(-1) + 32*ssz**2*lq2*s - 64*ssz**2*lq2*
     &    s**2*u2*sz**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(1,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 64*ssz**2*lq2*s**2*sz**(-1) + 32*ssz**2*lq2*t2 + 
     &    64*ssz**2*rq2*m1**2*s**2*sz**(-2) - 32*ssz**2*rq2*m1**2 - 64*
     &    ssz**2*rq2*mz**2*s*t2*sz**(-2) - 64*ssz**2*rq2*mz**2*s**2*
     &    sz**(-2) - 64*ssz**2*rq2*s*t2*u2*sz**(-2) - 64*ssz**2*rq2*s*
     &    t2*sz**(-1) + 32*ssz**2*rq2*s - 64*ssz**2*rq2*s**2*u2*
     &    sz**(-2) - 64*ssz**2*rq2*s**2*sz**(-1) + 32*ssz**2*rq2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(1,11,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h1**2*lambda1**2*mh1**2 + 8*h1**2*lambda1**2*s - 
     &    16*h1**2*lambda1**2*s**2*s1**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,11,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2
     &     + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s - 32*
     &    (mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**2*s1**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*mh1**2*t2t**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    mh1**2*s*t2t**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*mh1**2 + 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    s**2*s1**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s**2*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*mh1**2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t2t**(-1) - 16*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*s1**(-1)*t2t**(-1) - 16*
     &    h1**2*lambda1**2*s**2*s1**(-2) + 8*h1**2*lambda1**2 )
      MMpartial4 = MMpartial4 + ANGfin(1,12,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h2**2*lambda2**2*mh2**2 + 8*h2**2*lambda2**2*s - 
     &    16*h2**2*lambda2**2*s**2*s2**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,12,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2 - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 32*
     &    (mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**2*s1**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*mh2**2*t2t**(-1) - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    mh2**2*s*t2t**(-1) - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*mh2**2 + 8*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    s**2*s2**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s**2*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*mh2**2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(1,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t2t**(-1) - 16*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*s2**(-1)*t2t**(-1) - 32*
     &    h1*h2*lambda1*lambda2*s**2*s1**(-1)*s2**(-1) - 16*h2**2*
     &    lambda2**2*s**2*s2**(-2) + 8*h2**2*lambda2**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*
     &    t2t**(-1) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**2*s*
     &    t2t**(-2) + 8*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*s*t2t**(-2)
     &     + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-1) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**4*t2t**(-2) - 6*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*s*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 4*(s+u2-m12+mt2)**(-1)*hl**4*mt**4*t2t**(-2) + 
     &    2*(s+u2-m12+mt2)**(-1)*hl**4*s*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-1) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*t2t**(-2) - 6*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-1) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*s*t2t**(-1) - 2*hl**4*m1**2*
     &    t2t**(-2) + 2*hl**4*mt**2*t2t**(-2) + 4*hl**4*s*t2t**(-2) + 2
     &    *hl**4*t2t**(-1) - 2*hr**4*m1**2*t2t**(-2) + 2*hr**4*mt**2*
     &    t2t**(-2) + 4*hr**4*s*t2t**(-2) + 2*hr**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*
     &    t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**3*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl**2*
     &    hr**2*m1**2*mt**2*s*t2t**(-2) - 8*(s+u2-m12+mt2)**(-1)*
     &    hl**2*hr**2*mt**2*s*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*
     &    hl**2*hr**2*mt**4*s*t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*
     &    hl**4*m1**2*mt**2*s*t2t**(-2) + 12*(s+u2-m12+mt2)**(-1)*
     &    hl**4*m1**2*mt**4*t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*
     &    m1**2*s*t2t**(-1) - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2 - 12
     &    *(s+u2-m12+mt2)**(-1)*hl**4*m1**4*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**6*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 6*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*s*t2t**(-1) + 2
     &    *(s+u2-m12+mt2)**(-1)*hl**4*mt**2 - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*s*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**6*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*s*t2t**(-2) + 12*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**4*t2t**(-2) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2 - 12*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**6*t2t**(-2) + 6*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2 - 4*(s+u2-m12+mt2)**(-1)
     &    *hr**4*mt**4*s*t2t**(-2) - 4*(s+u2-m12+mt2)**(-1)*hr**4*
     &    mt**6*t2t**(-2) - 16*pq*hl**2*ssp*m1**2*t2t**(-1) + 8*pq*
     &    hl**2*ssp*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*pq*hl**2*ssp - 16*pq*hr**2*ssp*m1**2*t2t**(-1)
     &     + 8*pq*hr**2*ssp*mt**2*t2t**(-1) - 8*pq*hr**2*ssp - 16*lq*
     &    hl**2*ssz*m1**2*t2t**(-1) + 8*lq*hl**2*ssz*mt**2*t2t**(-1) - 
     &    8*lq*hl**2*ssz - 16*rq*hr**2*ssz*m1**2*t2t**(-1) + 8*rq*hr**2
     &    *ssz*mt**2*t2t**(-1) - 8*rq*hr**2*ssz + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    t2t**(-1) - 8*hl**2*hr**2*mt**2*s*t2t**(-2) - 8*hl**4*m1**2*
     &    mt**2*t2t**(-2) + 4*hl**4*m1**2*s*t2t**(-2) + 4*hl**4*m1**2*
     &    u2*t2t**(-2) + 4*hl**4*m1**4*t2t**(-2) - 4*hl**4*mt**2*u2*
     &    t2t**(-2) + 4*hl**4*mt**4*t2t**(-2) - 4*hl**4*s*u2*t2t**(-2)
     &     - 8*hl**4*s*t2t**(-1) - 4*hl**4*s**2*t2t**(-2) - 2*hl**4*u2*
     &    t2t**(-1) - 2*hl**4 - 8*hr**4*m1**2*mt**2*t2t**(-2) + 4*hr**4
     &    *m1**2*s*t2t**(-2) + 4*hr**4*m1**2*u2*t2t**(-2) + 4*hr**4*
     &    m1**4*t2t**(-2) - 4*hr**4*mt**2*u2*t2t**(-2) + 4*hr**4*mt**4*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(2,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*hr**4*s*u2*t2t**(-2) - 8*hr**4*s*t2t**(-1) - 4*
     &    hr**4*s**2*t2t**(-2) - 2*hr**4*u2*t2t**(-1) - 2*hr**4 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*ssp**2*pq2*s + 64*ssp**2*pq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*ssp**2*pq2*m1**2*s - 64*ssp**2*pq2*s*t2 - 128*
     &    ssp**2*pq2*s*u2 - 64*ssp**2*pq2*s**2 - 64*ssp**2*pq2*t2*u2 - 
     &    64*ssp**2*pq2*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2 + 32*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) + 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2 - 64*pq*lq*ssz*ssp*
     &    mz**(-2)*s )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 64*pq*lq*ssz*ssp*mz**(-2)*u2 - 64*pq*rq*ssz*ssp*
     &    mz**(-2)*s - 64*pq*rq*ssz*ssp*mz**(-2)*u2 - 8*pq*hl**2*ssp*
     &    m1**2*t2t**(-1) + 8*pq*hl**2*ssp*mt**2*t2t**(-1) + 16*pq*
     &    hl**2*ssp*s*t2t**(-1) + 8*pq*hl**2*ssp*u2*t2t**(-1) - 8*pq*
     &    hr**2*ssp*m1**2*t2t**(-1) + 8*pq*hr**2*ssp*mt**2*t2t**(-1) + 
     &    16*pq*hr**2*ssp*s*t2t**(-1) + 8*pq*hr**2*ssp*u2*t2t**(-1) + 
     &    64*ssp**2*pq2*s**(-1)*t2 + 64*ssp**2*pq2*s**(-1)*u2 + 64*
     &    ssp**2*pq2 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*s*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**2 + 48*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**4*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 - 
     &    48*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4 + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**6*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**6*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*s*t2t**(-1)
     &     + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2 + 48*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**4*t2t**(-1) + 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 48*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*
     &    mt**2*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4
     &     + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**6*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**6*t2t**(-1) - 64*pq*
     &    lq*ssz*ssp*m1**2*mz**(-2)*s + 64*pq*lq*ssz*ssp*mz**(-2)*s*t2
     &     + 128*pq*lq*ssz*ssp*mz**(-2)*s*u2 + 64*pq*lq*ssz*ssp*
     &    mz**(-2)*s**2 + 64*pq*lq*ssz*ssp*mz**(-2)*t2*u2 + 64*pq*lq*
     &    ssz*ssp*mz**(-2)*u2**2 - 64*pq*rq*ssz*ssp*m1**2*mz**(-2)*s + 
     &    64*pq*rq*ssz*ssp*mz**(-2)*s*t2 + 128*pq*rq*ssz*ssp*mz**(-2)*s
     &    *u2 + 64*pq*rq*ssz*ssp*mz**(-2)*s**2 + 64*pq*rq*ssz*ssp*
     &    mz**(-2)*t2*u2 + 64*pq*rq*ssz*ssp*mz**(-2)*u2**2 - 32*pq*
     &    hl**2*ssp*m1**2*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 40*pq*hl**2*ssp*m1**2*s*t2t**(-1) + 24*pq*hl**2*ssp*
     &    m1**2*u2*t2t**(-1) - 8*pq*hl**2*ssp*m1**2 + 16*pq*hl**2*ssp*
     &    m1**4*t2t**(-1) - 16*pq*hl**2*ssp*mt**2*s*t2t**(-1) - 24*pq*
     &    hl**2*ssp*mt**2*u2*t2t**(-1) + 8*pq*hl**2*ssp*mt**2 + 16*pq*
     &    hl**2*ssp*mt**4*t2t**(-1) - 24*pq*hl**2*ssp*s*u2*t2t**(-1) - 
     &    16*pq*hl**2*ssp*s - 16*pq*hl**2*ssp*s**2*t2t**(-1) + 8*pq*
     &    hl**2*ssp*t2 - 16*pq*hl**2*ssp*u2 - 8*pq*hl**2*ssp*u2**2*
     &    t2t**(-1) - 32*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) + 40*pq*
     &    hr**2*ssp*m1**2*s*t2t**(-1) + 24*pq*hr**2*ssp*m1**2*u2*
     &    t2t**(-1) - 8*pq*hr**2*ssp*m1**2 + 16*pq*hr**2*ssp*m1**4*
     &    t2t**(-1) - 16*pq*hr**2*ssp*mt**2*s*t2t**(-1) - 24*pq*hr**2*
     &    ssp*mt**2*u2*t2t**(-1) + 8*pq*hr**2*ssp*mt**2 + 16*pq*hr**2*
     &    ssp*mt**4*t2t**(-1) - 24*pq*hr**2*ssp*s*u2*t2t**(-1) - 16*pq*
     &    hr**2*ssp*s - 16*pq*hr**2*ssp*s**2*t2t**(-1) + 8*pq*hr**2*ssp
     &    *t2 )
      MMpartial4 = MMpartial4 + ANGfin(2,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*pq*hr**2*ssp*u2 - 8*pq*hr**2*ssp*u2**2*
     &    t2t**(-1) + 64*ssp**2*pq2*m1**2 - 128*ssp**2*pq2*s**(-1)*t2*
     &    u2 - 64*ssp**2*pq2*s - 128*ssp**2*pq2*t2 - 64*ssp**2*pq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-2)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*ssz**2*lq2*mz**2*s*sz**(-1) - 32*ssz**2*lq2*mz**2
     &     + 32*ssz**2*lq2*s*t2*sz**(-1) + 32*ssz**2*lq2*s*u2*sz**(-1)
     &     + 32*ssz**2*lq2*s**2*sz**(-1) - 32*ssz**2*lq2*t2 + 32*ssz**2
     &    *rq2*mz**2*s*sz**(-1) - 32*ssz**2*rq2*mz**2 + 32*ssz**2*rq2*s
     &    *t2*sz**(-1) + 32*ssz**2*rq2*s*u2*sz**(-1) + 32*ssz**2*rq2*
     &    s**2*sz**(-1) - 32*ssz**2*rq2*t2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*ssz**2*lq2*m1**2*mz**2 - 32*ssz**2*lq2*m1**2*s
     &     + 64*ssz**2*lq2*m1**2*s**2*sz**(-1) - 64*ssz**2*lq2*mz**2*s*
     &    u2*sz**(-1) + 32*ssz**2*lq2*mz**2*s - 64*ssz**2*lq2*mz**2*
     &    s**2*sz**(-1) + 32*ssz**2*lq2*mz**2*u2 - 64*ssz**2*lq2*s*t2*
     &    u2*sz**(-1) + 32*ssz**2*lq2*s*t2 - 64*ssz**2*lq2*s*u2 - 64*
     &    ssz**2*lq2*s**2*t2*sz**(-1) - 32*ssz**2*lq2*s**2 + 32*ssz**2*
     &    lq2*t2*u2 - 32*ssz**2*lq2*u2**2 - 32*ssz**2*rq2*m1**2*mz**2
     &     - 32*ssz**2*rq2*m1**2*s + 64*ssz**2*rq2*m1**2*s**2*sz**(-1)
     &     - 64*ssz**2*rq2*mz**2*s*u2*sz**(-1) + 32*ssz**2*rq2*mz**2*s
     &     - 64*ssz**2*rq2*mz**2*s**2*sz**(-1) + 32*ssz**2*rq2*mz**2*u2
     &     - 64*ssz**2*rq2*s*t2*u2*sz**(-1) + 32*ssz**2*rq2*s*t2 - 64*
     &    ssz**2*rq2*s*u2 - 64*ssz**2*rq2*s**2*t2*sz**(-1) - 32*ssz**2*
     &    rq2*s**2 + 32*ssz**2*rq2*t2*u2 - 32*ssz**2*rq2*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 32*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mz**2*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2
     &    *s*sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2 - 
     &    16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*t2t**(-1) + 8
     &    *(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mz**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s*t2*sz**(-1) + 32*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) + 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*t2t**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*t2t**(-1) + 8
     &    *(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mz**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s*t2*sz**(-1) + 64*pq*lq*
     &    ssz*ssp*mz**(-2)*s + 64*pq*lq*ssz*ssp*mz**(-2)*u2 + 64*pq*lq*
     &    ssz*ssp*mz**2*sz**(-1) + 64*pq*lq*ssz*ssp*s*sz**(-1) + 64*pq*
     &    lq*ssz*ssp*t2*sz**(-1) + 64*pq*lq*ssz*ssp*u2*sz**(-1) + 64*pq
     &    *rq*ssz*ssp*mz**(-2)*s + 64*pq*rq*ssz*ssp*mz**(-2)*u2 + 64*pq
     &    *rq*ssz*ssp*mz**2*sz**(-1) + 64*pq*rq*ssz*ssp*s*sz**(-1) + 64
     &    *pq*rq*ssz*ssp*t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*pq*rq*ssz*ssp*u2*sz**(-1) - 8*lq*hl**2*ssz*m1**2*
     &    s*sz**(-1)*t2t**(-1) + 8*lq*hl**2*ssz*mt**2*s*sz**(-1)*
     &    t2t**(-1) + 8*lq*hl**2*ssz*mz**2*s*sz**(-1)*t2t**(-1) + 8*lq*
     &    hl**2*ssz*s*u2*sz**(-1)*t2t**(-1) + 8*lq*hl**2*ssz*s*
     &    t2t**(-1) + 8*lq*hl**2*ssz*s**2*sz**(-1)*t2t**(-1) - 8*rq*
     &    hr**2*ssz*m1**2*s*sz**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*mt**2*s
     &    *sz**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*mz**2*s*sz**(-1)*
     &    t2t**(-1) + 8*rq*hr**2*ssz*s*u2*sz**(-1)*t2t**(-1) + 8*rq*
     &    hr**2*ssz*s*t2t**(-1) + 8*rq*hr**2*ssz*s**2*sz**(-1)*
     &    t2t**(-1) + 32*ssz**2*lq2*mz**2*s*sz**(-2) + 32*ssz**2*lq2*s*
     &    t2*sz**(-2) + 32*ssz**2*lq2*s*u2*sz**(-2) + 32*ssz**2*lq2*s*
     &    sz**(-1) + 32*ssz**2*lq2*s**2*sz**(-2) - 32*ssz**2*lq2 + 32*
     &    ssz**2*rq2*mz**2*s*sz**(-2) + 32*ssz**2*rq2*s*t2*sz**(-2) + 
     &    32*ssz**2*rq2*s*u2*sz**(-2) + 32*ssz**2*rq2*s*sz**(-1) + 32*
     &    ssz**2*rq2*s**2*sz**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 32*ssz**2*rq2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*
     &    mz**2*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2
     &    *mt**2*s*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2 + 48*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**4*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mz**2*s*sz**(-1) - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mz**2 + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s*t2*sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s**2*
     &    sz**(-1) - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 
     &    48*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*t2t**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4 + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**6*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*s*sz**(-1) + 
     &    8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*mz**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*t2
     &    *sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*t2 - 
     &    8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*mz**2*t2t**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*s*t2t**(-1) - 
     &    8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4 - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**6*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*mz**2*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2
     &    *s*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2 + 48*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**4*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mz**2*s*sz**(-1) - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mz**2 + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*
     &    t2*sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s
     &     - 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s**2*sz**(-1)
     &     )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*t2 - 
     &    48*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*t2t**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4 + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**6*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s*sz**(-1) + 
     &    8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2 - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*mz**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4 - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**6*t2t**(-1) + 64*pq*
     &    lq*ssz*ssp*m1**2*mz**(-2)*s + 128*pq*lq*ssz*ssp*m1**2*s*
     &    sz**(-1) - 64*pq*lq*ssz*ssp*m1**2 - 64*pq*lq*ssz*ssp*mz**(-2)
     &    *s*t2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 128*pq*lq*ssz*ssp*mz**(-2)*s*u2 - 64*pq*lq*ssz*
     &    ssp*mz**(-2)*s**2 - 64*pq*lq*ssz*ssp*mz**(-2)*t2*u2 - 64*pq*
     &    lq*ssz*ssp*mz**(-2)*u2**2 - 128*pq*lq*ssz*ssp*mz**2*s*
     &    sz**(-1) - 128*pq*lq*ssz*ssp*mz**2*u2*sz**(-1) - 128*pq*lq*
     &    ssz*ssp*s*t2*sz**(-1) - 64*pq*lq*ssz*ssp*s - 128*pq*lq*ssz*
     &    ssp*t2*u2*sz**(-1) - 64*pq*lq*ssz*ssp*u2 + 64*pq*rq*ssz*ssp*
     &    m1**2*mz**(-2)*s + 128*pq*rq*ssz*ssp*m1**2*s*sz**(-1) - 64*pq
     &    *rq*ssz*ssp*m1**2 - 64*pq*rq*ssz*ssp*mz**(-2)*s*t2 - 128*pq*
     &    rq*ssz*ssp*mz**(-2)*s*u2 - 64*pq*rq*ssz*ssp*mz**(-2)*s**2 - 
     &    64*pq*rq*ssz*ssp*mz**(-2)*t2*u2 - 64*pq*rq*ssz*ssp*mz**(-2)*
     &    u2**2 - 128*pq*rq*ssz*ssp*mz**2*s*sz**(-1) - 128*pq*rq*ssz*
     &    ssp*mz**2*u2*sz**(-1) - 128*pq*rq*ssz*ssp*s*t2*sz**(-1) - 64*
     &    pq*rq*ssz*ssp*s - 128*pq*rq*ssz*ssp*t2*u2*sz**(-1) - 64*pq*rq
     &    *ssz*ssp*u2 - 32*lq*hl**2*ssz*m1**2*mt**2*t2t**(-1) - 16*lq*
     &    hl**2*ssz*m1**2*mz**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*lq*hl**2*ssz*m1**2*s*u2*sz**(-1)*t2t**(-1) + 8*lq
     &    *hl**2*ssz*m1**2*s*t2t**(-1) + 32*lq*hl**2*ssz*m1**2*s**2*
     &    sz**(-1)*t2t**(-1) + 8*lq*hl**2*ssz*m1**2*u2*t2t**(-1) - 8*lq
     &    *hl**2*ssz*m1**2 + 16*lq*hl**2*ssz*m1**4*t2t**(-1) + 8*lq*
     &    hl**2*ssz*mt**2*mz**2*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*s*u2*
     &    sz**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*mt**2*s**2*sz**(-1)*
     &    t2t**(-1) - 8*lq*hl**2*ssz*mt**2*u2*t2t**(-1) + 8*lq*hl**2*
     &    ssz*mt**2 + 16*lq*hl**2*ssz*mt**4*t2t**(-1) - 16*lq*hl**2*ssz
     &    *mz**2*s*u2*sz**(-1)*t2t**(-1) + 16*lq*hl**2*ssz*mz**2*s*
     &    sz**(-1) - 16*lq*hl**2*ssz*mz**2*s**2*sz**(-1)*t2t**(-1) - 8*
     &    lq*hl**2*ssz*mz**2 + 16*lq*hl**2*ssz*s*t2*sz**(-1) - 16*lq*
     &    hl**2*ssz*s*u2*sz**(-1) - 24*lq*hl**2*ssz*s*u2*t2t**(-1) - 16
     &    *lq*hl**2*ssz*s**2*sz**(-1) - 16*lq*hl**2*ssz*s**2*t2t**(-1)
     &     - 8*lq*hl**2*ssz*t2 - 8*lq*hl**2*ssz*u2**2*t2t**(-1) - 32*rq
     &    *hr**2*ssz*m1**2*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*rq*hr**2*ssz*m1**2*mz**2*t2t**(-1) + 16*rq*
     &    hr**2*ssz*m1**2*s*u2*sz**(-1)*t2t**(-1) + 8*rq*hr**2*ssz*
     &    m1**2*s*t2t**(-1) + 32*rq*hr**2*ssz*m1**2*s**2*sz**(-1)*
     &    t2t**(-1) + 8*rq*hr**2*ssz*m1**2*u2*t2t**(-1) - 8*rq*hr**2*
     &    ssz*m1**2 + 16*rq*hr**2*ssz*m1**4*t2t**(-1) + 8*rq*hr**2*ssz*
     &    mt**2*mz**2*t2t**(-1) - 16*rq*hr**2*ssz*mt**2*s*u2*sz**(-1)*
     &    t2t**(-1) - 16*rq*hr**2*ssz*mt**2*s**2*sz**(-1)*t2t**(-1) - 8
     &    *rq*hr**2*ssz*mt**2*u2*t2t**(-1) + 8*rq*hr**2*ssz*mt**2 + 16*
     &    rq*hr**2*ssz*mt**4*t2t**(-1) - 16*rq*hr**2*ssz*mz**2*s*u2*
     &    sz**(-1)*t2t**(-1) + 16*rq*hr**2*ssz*mz**2*s*sz**(-1) - 16*rq
     &    *hr**2*ssz*mz**2*s**2*sz**(-1)*t2t**(-1) - 8*rq*hr**2*ssz*
     &    mz**2 + 16*rq*hr**2*ssz*s*t2*sz**(-1) - 16*rq*hr**2*ssz*s*u2*
     &    sz**(-1) - 24*rq*hr**2*ssz*s*u2*t2t**(-1) - 16*rq*hr**2*ssz*
     &    s**2*sz**(-1) - 16*rq*hr**2*ssz*s**2*t2t**(-1) - 8*rq*hr**2*
     &    ssz*t2 )
      MMpartial4 = MMpartial4 + ANGfin(2,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*rq*hr**2*ssz*u2**2*t2t**(-1) + 64*ssz**2*lq2*
     &    m1**2*s**2*sz**(-2) - 32*ssz**2*lq2*m1**2 - 64*ssz**2*lq2*
     &    mz**2*s*u2*sz**(-2) - 64*ssz**2*lq2*mz**2*s**2*sz**(-2) - 64*
     &    ssz**2*lq2*s*t2*u2*sz**(-2) - 64*ssz**2*lq2*s*u2*sz**(-1) + 
     &    32*ssz**2*lq2*s - 64*ssz**2*lq2*s**2*t2*sz**(-2) - 64*ssz**2*
     &    lq2*s**2*sz**(-1) + 32*ssz**2*lq2*u2 + 64*ssz**2*rq2*m1**2*
     &    s**2*sz**(-2) - 32*ssz**2*rq2*m1**2 - 64*ssz**2*rq2*mz**2*s*
     &    u2*sz**(-2) - 64*ssz**2*rq2*mz**2*s**2*sz**(-2) - 64*ssz**2*
     &    rq2*s*t2*u2*sz**(-2) - 64*ssz**2*rq2*s*u2*sz**(-1) + 32*
     &    ssz**2*rq2*s - 64*ssz**2*rq2*s**2*t2*sz**(-2) - 64*ssz**2*rq2
     &    *s**2*sz**(-1) + 32*ssz**2*rq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(2,11,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h1**2*lambda1**2*mh1**2 + 8*h1**2*lambda1**2*s - 
     &    16*h1**2*lambda1**2*s**2*s1**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,11,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*mh1**2
     &     + 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s - 32*
     &    (mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**2*s1**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    mh1**2*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s + 16*(s+u2-m12+mt2)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*s1**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*
     &    mh1**2*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s*t2t**(-1) + 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*mh1**2*t2t**(-1) - 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2t**(-1) - 16*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt*s**2*s1**(-1)*t2t**(-1) - 16*h1**2*lambda1**2*s**2*
     &    s1**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(2,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h1**2*lambda1**2 )
      MMpartial4 = MMpartial4 + ANGfin(2,12,-1,-2)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*h2**2*lambda2**2*mh2**2 + 8*h2**2*lambda2**2*s - 
     &    16*h2**2*lambda2**2*s**2*s2**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,12,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*
     &    mh2**2 - 16*(mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s + 32*
     &    (mh12-mh22)**(-1)*h1*h2*lambda1*lambda2*s**2*s1**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    mh2**2*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s + 16*(s+u2-m12+mt2)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*s2**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*
     &    mh2**2*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s*t2t**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mh2**2*t2t**(-1) - 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2t**(-1) - 16*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*s**2*s2**(-1)*t2t**(-1) - 32*h1*h2*lambda1*lambda2*s**2*
     &    s1**(-1)*s2**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(2,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*h2**2*lambda2**2*s**2*s2**(-2) + 8*h2**2*
     &    lambda2**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-2,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 4*hl**2*hr**2*mt**2*s*t2t**(-2) - 4*hl**4*m1**2*
     &    mt**2*t2t**(-2) - 2*hl**4*m1**2*u2*t2t**(-2) + 2*hl**4*m1**4*
     &    t2t**(-2) + 2*hl**4*mt**2*s*t2t**(-2) + 2*hl**4*mt**2*u2*
     &    t2t**(-2) + 2*hl**4*mt**4*t2t**(-2) - 4*hr**4*m1**2*mt**2*
     &    t2t**(-2) - 2*hr**4*m1**2*u2*t2t**(-2) + 2*hr**4*m1**4*
     &    t2t**(-2) + 2*hr**4*mt**2*s*t2t**(-2) + 2*hr**4*mt**2*u2*
     &    t2t**(-2) + 2*hr**4*mt**4*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hl**2*hr**2*m1**2*mt**2*s*t2t**(-2) + 4*hl**2*
     &    hr**2*m1**2*mt**2*t2t**(-1) + 4*hl**2*hr**2*mt**2*s*u2*
     &    t2t**(-2) + 4*hl**2*hr**2*mt**2*s**2*t2t**(-2) - 4*hl**2*
     &    hr**2*mt**4*t2t**(-1) + 8*hl**4*m1**2*mt**2*s*t2t**(-2) + 8*
     &    hl**4*m1**2*mt**2*u2*t2t**(-2) - 10*hl**4*m1**2*mt**2*
     &    t2t**(-1) - 2*hl**4*m1**2*mt**4*t2t**(-2) + 2*hl**4*m1**2*s*
     &    u2*t2t**(-2) + 2*hl**4*m1**2*u2**2*t2t**(-2) - 4*hl**4*m1**2
     &     - 2*hl**4*m1**4*mt**2*t2t**(-2) - 2*hl**4*m1**4*s*t2t**(-2)
     &     - 4*hl**4*m1**4*u2*t2t**(-2) + 4*hl**4*m1**4*t2t**(-1) + 2*
     &    hl**4*m1**6*t2t**(-2) - 4*hl**4*mt**2*s*u2*t2t**(-2) - 2*
     &    hl**4*mt**2*s**2*t2t**(-2) - 2*hl**4*mt**2*u2**2*t2t**(-2) + 
     &    4*hl**4*mt**2 - 2*hl**4*mt**4*s*t2t**(-2) - 4*hl**4*mt**4*u2*
     &    t2t**(-2) + 6*hl**4*mt**4*t2t**(-1) + 2*hl**4*mt**6*t2t**(-2)
     &     + 8*hr**4*m1**2*mt**2*s*t2t**(-2) + 8*hr**4*m1**2*mt**2*u2*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 10*hr**4*m1**2*mt**2*t2t**(-1) - 2*hr**4*m1**2*
     &    mt**4*t2t**(-2) + 2*hr**4*m1**2*s*u2*t2t**(-2) + 2*hr**4*
     &    m1**2*u2**2*t2t**(-2) - 4*hr**4*m1**2 - 2*hr**4*m1**4*mt**2*
     &    t2t**(-2) - 2*hr**4*m1**4*s*t2t**(-2) - 4*hr**4*m1**4*u2*
     &    t2t**(-2) + 4*hr**4*m1**4*t2t**(-1) + 2*hr**4*m1**6*t2t**(-2)
     &     - 4*hr**4*mt**2*s*u2*t2t**(-2) - 2*hr**4*mt**2*s**2*
     &    t2t**(-2) - 2*hr**4*mt**2*u2**2*t2t**(-2) + 4*hr**4*mt**2 - 2
     &    *hr**4*mt**4*s*t2t**(-2) - 4*hr**4*mt**4*u2*t2t**(-2) + 6*
     &    hr**4*mt**4*t2t**(-1) + 2*hr**4*mt**6*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*
     &    t2t**(-1) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**2*s*
     &    t2t**(-2) + 16*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*s*t2t**(-2)
     &     + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*u2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**4*t2t**(-2) - 6*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*s*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 4*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*u2*t2t**(-2)
     &     - 4*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*s*t2t**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*u2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*t2t**(-2) - 6*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*u2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*s*t2t**(-1) - 4*hl**4*m1**2*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 4*hl**4*mt**2*t2t**(-2) + 2*hl**4*s*t2t**(-2) + 2*
     &    hl**4*u2*t2t**(-2) - 4*hr**4*m1**2*t2t**(-2) + 4*hr**4*mt**2*
     &    t2t**(-2) + 2*hr**4*s*t2t**(-2) + 2*hr**4*u2*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2*
     &    t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*t2t**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*u2 - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*t2t**(-1)
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*
     &    t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    u2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1)
     &     - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1)
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*
     &    t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    u2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*t2t**(-1)
     &     - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) - 8
     &    *(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2 + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2 - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*m1**2*mt**2*
     &    s*t2t**(-2) - 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*m1**2*
     &    mt**2*u2*t2t**(-2) - 8*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    m1**2*mt**4*t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    m1**4*mt**2*t2t**(-2) - 8*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    mt**2*s*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    mt**4*s*t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    mt**4*u2*t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*
     &    mt**6*t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    s*t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*u2*
     &    t2t**(-2) - 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**4*
     &    t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*s*t2t**(-1)
     &     - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*u2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*m1**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 14*(s+u2-m12+mt2)**(-1)*hl**4*m1**4*mt**2*
     &    t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**4*t2t**(-1) + 
     &    4*(s+u2-m12+mt2)**(-1)*hl**4*m1**6*t2t**(-2) + 6*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*s*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*u2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2 - 4*(s+u2-m12+mt2)**(-1)
     &    *hl**4*mt**4*s*t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*
     &    mt**4*u2*t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*mt**4*
     &    t2t**(-1) - 6*(s+u2-m12+mt2)**(-1)*hl**4*mt**6*t2t**(-2) + 
     &    4*(s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*s*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*u2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**4*t2t**(-2) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*s*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 2*(s+u2-m12+mt2)**(-1)*hr**4*m1**2 - 14*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*mt**2*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**4*t2t**(-1) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**6*t2t**(-2) + 6*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*s*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*u2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2 - 4*(s+u2-m12+mt2)**(-1)
     &    *hr**4*mt**4*s*t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hr**4*
     &    mt**4*u2*t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hr**4*mt**4*
     &    t2t**(-1) - 6*(s+u2-m12+mt2)**(-1)*hr**4*mt**6*t2t**(-2) - 
     &    8*pq*hl**2*ssp*mt**2*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*
     &    t2t**(-1) - 8*lq*hl**2*ssz*mt**2*t2t**(-1) - 8*rq*hr**2*ssz*
     &    mt**2*t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*t2t**(-1)
     &     + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*t2t**(-1) + 4*hl**2*
     &    hr**2*m1**2*mt**2*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hl**2*hr**2*mt**2*s*t2t**(-2) - 4*hl**2*hr**2*
     &    mt**4*t2t**(-2) - 14*hl**4*m1**2*mt**2*t2t**(-2) + 4*hl**4*
     &    m1**2*s*t2t**(-2) + 4*hl**4*m1**2*u2*t2t**(-2) - 2*hl**4*
     &    m1**2*t2t**(-1) + 6*hl**4*m1**4*t2t**(-2) - 4*hl**4*mt**2*u2*
     &    t2t**(-2) + 2*hl**4*mt**2*t2t**(-1) + 8*hl**4*mt**4*t2t**(-2)
     &     - 4*hl**4*s*u2*t2t**(-2) - 4*hl**4*s*t2t**(-1) - 2*hl**4*
     &    s**2*t2t**(-2) - 2*hl**4*u2*t2t**(-1) - 2*hl**4*u2**2*
     &    t2t**(-2) - 14*hr**4*m1**2*mt**2*t2t**(-2) + 4*hr**4*m1**2*s*
     &    t2t**(-2) + 4*hr**4*m1**2*u2*t2t**(-2) - 2*hr**4*m1**2*
     &    t2t**(-1) + 6*hr**4*m1**4*t2t**(-2) - 4*hr**4*mt**2*u2*
     &    t2t**(-2) + 2*hr**4*mt**2*t2t**(-1) + 8*hr**4*mt**4*t2t**(-2)
     &     - 4*hr**4*s*u2*t2t**(-2) - 4*hr**4*s*t2t**(-1) - 2*hr**4*
     &    s**2*t2t**(-2) - 2*hr**4*u2*t2t**(-1) - 2*hr**4*u2**2*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-2,1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 2*hl**4*m1**2*t2t**(-2) + 2*hl**4*mt**2*t2t**(-2)
     &     - 2*hr**4*m1**2*t2t**(-2) + 2*hr**4*mt**2*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-2,1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 4*hl**2*hr**2*m1**2*mt**2*t2t**(-2) - 4*hl**2*hr**2*
     &    mt**4*t2t**(-2) + 2*hl**4*m1**2*mt**2*t2t**(-2) + 2*hl**4*
     &    m1**2*s*t2t**(-2) + 2*hl**4*m1**2*u2*t2t**(-2) - 2*hl**4*
     &    m1**4*t2t**(-2) - 2*hl**4*mt**2*s*t2t**(-2) - 2*hl**4*mt**2*
     &    u2*t2t**(-2) + 2*hr**4*m1**2*mt**2*t2t**(-2) + 2*hr**4*m1**2*
     &    s*t2t**(-2) + 2*hr**4*m1**2*u2*t2t**(-2) - 2*hr**4*m1**4*
     &    t2t**(-2) - 2*hr**4*mt**2*s*t2t**(-2) - 2*hr**4*mt**2*u2*
     &    t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) + 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2*
     &    t2t**(-1) - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*
     &    t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2*
     &    t2t**(-1) - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*
     &    t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*s + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2 + 8*(u2-m12+mt2)**(-1)
     &    *pq*hl**2*ssp*u2 + 64*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2
     &    *mt**2*t2t**(-1) + 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2
     &    *u2*t2t**(-1) - 32*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*
     &    t2t**(-1) - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*
     &    t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2*
     &    t2t**(-1) - 32*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*s + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2 + 8*(u2-m12+mt2)**(-1)
     &    *pq*hr**2*ssp*u2 + 64*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*
     &    m1**2*mt**2*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*
     &    ssp*m1**2*s*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*
     &    ssp*m1**2*u2*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*
     &    ssp*m1**2 - 32*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*
     &    t2t**(-1) - 24*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 - 32
     &    *(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*s - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*u2 + 64*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*s*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2*
     &    t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 - 32
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) - 24*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 - 32*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*s - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*u2 - 16*pq*hl**2*ssp*
     &    m1**2*t2t**(-1) + 16*pq*hl**2*ssp*mt**2*t2t**(-1) - 16*pq*
     &    hr**2*ssp*m1**2*t2t**(-1) + 16*pq*hr**2*ssp*mt**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 88*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    s*t2t**(-1) - 200*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**2*u2*t2t**(-1) - 104*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*
     &    m1**2*mt**2 - 256*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**4*t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    s*u2*t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*s
     &     - 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 - 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2 - 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2**2*t2t**(-1) + 272
     &    *(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1) + 32
     &    *(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*s*t2t**(-1) + 112*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*u2*t2t**(-1) + 48*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4 - 96*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**6*t2t**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s**2*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2 + 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2**2*t2t**(-1) + 56*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*s*t2t**(-1) + 88*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*u2*t2t**(-1) + 56*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4 + 80*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**6*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s*t2 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s*u2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*s**2 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2*u2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*t2**2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 88*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*
     &    s*t2t**(-1) - 200*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    mt**2*u2*t2t**(-1) - 104*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*
     &    m1**2*mt**2 - 256*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    mt**4*t2t**(-1) - 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    s*u2*t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*s
     &     - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 - 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2 - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2**2*t2t**(-1) + 272
     &    *(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*mt**2*t2t**(-1) + 32
     &    *(u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*s*t2t**(-1) + 112*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*u2*t2t**(-1) + 48*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4 - 96*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**6*t2t**(-1) + 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s**2*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2 + 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2**2*t2t**(-1) + 56*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*s*t2t**(-1) + 88*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*u2*t2t**(-1) + 56*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4 + 80*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**6*t2t**(-1) - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s*t2 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s*u2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*s**2 - 16*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2*u2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*t2**2 - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*s
     &    *t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**2*u2*t2t**(-1) + 32*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*
     &    m1**2*mt**2 + 64*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*
     &    mt**4*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2
     &    *s + 8*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2 - 56*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*mt**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4 + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**6*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*s - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2 - 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4 - 24
     &    *(s+u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**6*t2t**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*s*t2t**(-1)
     &     + 8*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*u2*
     &    t2t**(-1) + 32*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*
     &    mt**2 + 64*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**4*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*s + 8
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2 - 56*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*mt**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4 + 16*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**6*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*s - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*s*
     &    t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*u2*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4 - 24
     &    *(s+u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**6*t2t**(-1) + 40*pq*
     &    hl**2*ssp*m1**2*mt**2*t2t**(-1) + 16*pq*hl**2*ssp*m1**2*u2*
     &    t2t**(-1) - 16*pq*hl**2*ssp*m1**4*t2t**(-1) - 8*pq*hl**2*ssp*
     &    mt**2*s*t2t**(-1) - 8*pq*hl**2*ssp*mt**2*u2*t2t**(-1) - 8*pq*
     &    hl**2*ssp*mt**2 - 24*pq*hl**2*ssp*mt**4*t2t**(-1) + 40*pq*
     &    hr**2*ssp*m1**2*mt**2*t2t**(-1) + 16*pq*hr**2*ssp*m1**2*u2*
     &    t2t**(-1) - 16*pq*hr**2*ssp*m1**4*t2t**(-1) - 8*pq*hr**2*ssp*
     &    mt**2*s*t2t**(-1) - 8*pq*hr**2*ssp*mt**2*u2*t2t**(-1) - 8*pq*
     &    hr**2*ssp*mt**2 - 24*pq*hr**2*ssp*mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 4*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-2) + 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-2) - 4*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-2) + 2*hl**4*
     &    t2t**(-2) + 2*hr**4*t2t**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(4,5,-1,1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*m1**2*mt**2*
     &    t2t**(-2) + 4*(s+u2-m12+mt2)**(-1)*hl**2*hr**2*mt**4*
     &    t2t**(-2) + 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*mt**2*
     &    t2t**(-2) - 2*(s+u2-m12+mt2)**(-1)*hl**4*m1**2*t2t**(-1) + 
     &    2*(s+u2-m12+mt2)**(-1)*hl**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hl**4*mt**4*t2t**(-2) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*mt**2*t2t**(-2) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*m1**2*t2t**(-1) + 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**2*t2t**(-1) - 2*
     &    (s+u2-m12+mt2)**(-1)*hr**4*mt**4*t2t**(-2) - 2*hl**4*m1**2*
     &    t2t**(-2) + 2*hl**4*mt**2*t2t**(-2) - 2*hl**4*s*t2t**(-2) - 2
     &    *hl**4*u2*t2t**(-2) - 2*hl**4*t2t**(-1) - 2*hr**4*m1**2*
     &    t2t**(-2) + 2*hr**4*mt**2*t2t**(-2) - 2*hr**4*s*t2t**(-2) - 2
     &    *hr**4*u2*t2t**(-2) - 2*hr**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mz**2*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2*u2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)
     &    *lq*hl**2*ssz*mt**2*s*t2t**(-1) - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2*t2t**(-1) - 
     &    32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1) + 8
     &    *(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*s*sz**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2*sz**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*u2*sz**(-1) + 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2*sz**(-1) + 64*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2*
     &    t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    mz**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*s*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**2*u2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mt**4*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*mz**2*s*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*s*t2*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    s*u2*sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2
     &    *sz**(-1) + 64*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mz**2*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*
     &    ssz*m1**2*u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2 - 32*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*mz**2*t2t**(-1) - 
     &    24*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*t2t**(-1) - 16
     &    *(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*u2*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2 - 32*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mz**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s*t2*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s*u2*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*s**2*sz**(-1) + 64*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) + 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*t2t**(-1)
     &     + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*s*t2t**(-1)
     &     + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*u2*t2t**(-1)
     &     )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2 - 32*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*t2t**(-1) - 
     &    24*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2t**(-1) - 16
     &    *(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*u2*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2 - 32*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mz**2*s*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s*t2*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s*u2*sz**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*s**2*sz**(-1) - 16*lq*
     &    hl**2*ssz*m1**2*t2t**(-1) + 16*lq*hl**2*ssz*mt**2*t2t**(-1)
     &     - 16*rq*hr**2*ssz*m1**2*t2t**(-1) + 16*rq*hr**2*ssz*mt**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 152*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*mz**2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2*mt**2*s*sz**(-1) - 88*(u2-m12+mt2+mz2)**(-1)*lq
     &    *hl**2*ssz*m1**2*mt**2*s*t2t**(-1) - 200*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*u2*
     &    t2t**(-1) - 72*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2 - 256*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**4*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mz**2*s*sz**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mz**2*s*t2t**(-1) - 48*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*u2*
     &    t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mz**2 - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*mz**4
     &    *t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    s*t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*s
     &    *u2*sz**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*s*u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    t2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*u2**2*t2t**(-1)
     &     + 272*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*
     &    t2t**(-1) + 80*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**4*
     &    mz**2*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s*sz**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*s*t2t**(-1) + 112*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*
     &    ssz*m1**4*u2*t2t**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**4 - 96*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**6*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*mz**2*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    mz**2*s*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz
     &    *mt**2*mz**2*u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*mz**2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*mt**2*mz**4*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2*s*t2*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*
     &    lq*hl**2*ssz*mt**2*s*u2*sz**(-1) + 32*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s*u2*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*s**2*
     &    sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*
     &    s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*t2 + 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2
     &     + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**2*u2**2*
     &    t2t**(-1) + 72*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    mz**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s*
     &    sz**(-1) + 56*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*s*
     &    t2t**(-1) + 88*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**4*
     &    u2*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**4 + 80*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mt**6*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*s
     &     - 16*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*mz**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*s**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*t2**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 152*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*mz**2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*m1**2*mt**2*s*sz**(-1) - 88*(u2-m12+mt2+mz2)**(-1)*rq
     &    *hr**2*ssz*m1**2*mt**2*s*t2t**(-1) - 200*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*u2*
     &    t2t**(-1) - 72*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2 - 256*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**4*t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mz**2*s*sz**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**2*s*t2t**(-1) - 48*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*u2*
     &    t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mz**2 - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**4
     &    *t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    s*t2*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*s
     &    *u2*sz**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**2*s*u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    t2 - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2 - 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*u2**2*t2t**(-1)
     &     + 272*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*
     &    t2t**(-1) + 80*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*
     &    mz**2*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s*sz**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**4*s*t2t**(-1) + 112*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*
     &    ssz*m1**4*u2*t2t**(-1) + 32*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*m1**4 - 96*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    m1**6*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*mz**2*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    mz**2*s*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz
     &    *mt**2*mz**2*u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*mz**2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2
     &    *ssz*mt**2*mz**4*t2t**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*rq*
     &    hr**2*ssz*mt**2*s*t2*sz**(-1) + 16*(u2-m12+mt2+mz2)**(-1)*
     &    rq*hr**2*ssz*mt**2*s*u2*sz**(-1) + 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s*u2*t2t**(-1)
     &     + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*s**2*
     &    sz**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*
     &    s**2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**2*t2 + 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2
     &     + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2**2*
     &    t2t**(-1) + 72*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*
     &    mz**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s*
     &    sz**(-1) + 56*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*s*
     &    t2t**(-1) + 88*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*
     &    u2*t2t**(-1) + 40*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*
     &    mt**4 + 80*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**6*
     &    t2t**(-1) - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*s
     &     - 16*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**4 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*t2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*s**2 - 16*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2*u2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*t2**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2**2 )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**2*
     &    mz**2*t2t**(-1) + 32*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mt**2*s*sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*
     &    ssz*m1**2*mt**2*s*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*lq*
     &    hl**2*ssz*m1**2*mt**2*u2*t2t**(-1) + 64*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mt**4*t2t**(-1) + 
     &    16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mz**2*s*sz**(-1)
     &     - 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*mz**2 + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*t2*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s*u2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*s - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*t2 - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**2*u2 - 56*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4*mt**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**4*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*m1**6*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*
     &    mz**2*s*sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*
     &    mt**2*mz**2 - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*
     &    t2*sz**(-1) - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s*
     &    u2*sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s
     &     - 16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*s**2*sz**(-1)
     &     + 8*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**2*u2 - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*mz**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*s*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**4*u2*t2t**(-1) - 24*
     &    (s+u2-m12+mt2)**(-1)*lq*hl**2*ssz*mt**6*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*mz**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 32*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*s
     &    *sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*
     &    mt**2*s*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*
     &    m1**2*mt**2*u2*t2t**(-1) + 64*(s+u2-m12+mt2)**(-1)*rq*hr**2
     &    *ssz*m1**2*mt**4*t2t**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*
     &    hr**2*ssz*m1**2*mz**2*s*sz**(-1) - 8*(s+u2-m12+mt2)**(-1)*
     &    rq*hr**2*ssz*m1**2*mz**2 + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2
     &    *ssz*m1**2*s*t2*sz**(-1) + 16*(s+u2-m12+mt2)**(-1)*rq*hr**2
     &    *ssz*m1**2*s*u2*sz**(-1) + 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*
     &    ssz*m1**2*s - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*t2
     &     - 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**2*u2 - 56*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4*mt**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**4*s*sz**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*m1**6*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*mz**2 - 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*t2*sz**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s*u2*sz**(-1) + 
     &    8*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*s**2*sz**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*t2 + 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**2*u2 - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*mz**2*t2t**(-1) - 
     &    16*(s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*s*sz**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*s*t2t**(-1) - 8*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**4*u2*t2t**(-1) - 24*
     &    (s+u2-m12+mt2)**(-1)*rq*hr**2*ssz*mt**6*t2t**(-1) + 40*lq*
     &    hl**2*ssz*m1**2*mt**2*t2t**(-1) + 16*lq*hl**2*ssz*m1**2*u2*
     &    t2t**(-1) - 16*lq*hl**2*ssz*m1**4*t2t**(-1) - 8*lq*hl**2*ssz*
     &    mt**2*mz**2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,10,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*lq*hl**2*ssz*mt**2*s*t2t**(-1) - 8*lq*hl**2*ssz
     &    *mt**2*u2*t2t**(-1) - 8*lq*hl**2*ssz*mt**2 - 24*lq*hl**2*ssz*
     &    mt**4*t2t**(-1) + 40*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1) + 16*
     &    rq*hr**2*ssz*m1**2*u2*t2t**(-1) - 16*rq*hr**2*ssz*m1**4*
     &    t2t**(-1) - 8*rq*hr**2*ssz*mt**2*mz**2*t2t**(-1) - 8*rq*hr**2
     &    *ssz*mt**2*s*t2t**(-1) - 8*rq*hr**2*ssz*mt**2*u2*t2t**(-1) - 
     &    8*rq*hr**2*ssz*mt**2 - 24*rq*hr**2*ssz*mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,11,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*mh1**2*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    s*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*u2*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt**3*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4*mt*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s + 16*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s**2*s1**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**3*mh1**2*t2t**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*s*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*u2*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt**5*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*mh1**2*t2t**(-1) + 24*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt*s*t2t**(-1) - 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*u2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2
     &    *mt - 32*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt**3*t2t**(-1) + 16*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*m1**4
     &    *mt*t2t**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*mh1**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    mh1**2*u2*t2t**(-1) - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s
     &     )
      MMpartial4 = MMpartial4 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s**2*s1**(-1) - 8*(u2-m12+mt2+mh12)**(-1)*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2
     &     + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt*u2**2*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*mh1**2*t2t**(-1) - 24*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) + 24*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt**3*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3
     &     + 16*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *mt**5*t2t**(-1) + 24*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    t2t**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,11,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t2t**(-1) - 8
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2*t2t**(-1) + 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt - 24*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt**3*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,12,-1,-1)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*mh2**2*t2t**(-1) - 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    s*t2t**(-1) - 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*u2*t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**3*t2t**(-1) + 8*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4*mt*
     &    t2t**(-1) - 16*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s + 16*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*s**2*s2**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**3*mh2**2*t2t**(-1) + 16*
     &    (s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*s*
     &    t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*u2*t2t**(-1) + 8*(s+u2-m12+mt2)**(-1)*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt**5*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*mh2**2*t2t**(-1) + 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt*s*t2t**(-1) - 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*u2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2
     &    *mt - 32*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt**3*t2t**(-1) + 16*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*m1**4
     &    *mt*t2t**(-1) - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt*mh2**2*s*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*
     &    mh2**2*u2*t2t**(-1) - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s
     &     )
      MMpartial4 = MMpartial4 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s**2*s2**(-1) - 8*(u2-m12+mt2+mh22)**(-1)*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*u2
     &     + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*u2**2*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*mh2**2*t2t**(-1) - 24*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &    *s*t2t**(-1) + 24*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt**3*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3
     &     + 16*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt**5*t2t**(-1) + 24*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    t2t**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*
     &    t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(4,12,-1,-1)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t2t**(-1) - 8
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*u2*t2t**(-1) + 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt - 24*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**3*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(5,0,-2,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 128*ssp**2*pq2 )
      MMpartial4 = MMpartial4 + ANGfin(5,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 128*ssp**2*pq2*m1**2 - 128*ssp**2*pq2*s - 128*ssp**2
     &    *pq2*t2 - 128*ssp**2*pq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(5,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * (  - 128*pq*lq*ssz*ssp*mz**(-2) - 128*pq*rq*ssz*ssp*
     &    mz**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(5,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 56*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*mt**2*
     &    t2t**(-1) - 32*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2*u2*
     &    t2t**(-1) + 8*(u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**2 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*m1**4*t2t**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**2 + 24*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*mt**4*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hl**2*ssp*u2 - 56*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*mt**2*t2t**(-1) - 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2*u2*t2t**(-1) + 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**2 + 32*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*m1**4*t2t**(-1) + 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2*u2*t2t**(-1) - 8*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**2 + 24*
     &    (u2-m12+mt2)**(-1)*pq*hr**2*ssp*mt**4*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(5,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2)**(-1)*pq*hr**2*ssp*u2 - 128*pq*
     &    lq*ssz*ssp*m1**2*mz**(-2) + 128*pq*lq*ssz*ssp*mz**(-2)*s + 
     &    128*pq*lq*ssz*ssp*mz**(-2)*t2 + 128*pq*lq*ssz*ssp*mz**(-2)*u2
     &     - 128*pq*rq*ssz*ssp*m1**2*mz**(-2) + 128*pq*rq*ssz*ssp*
     &    mz**(-2)*s + 128*pq*rq*ssz*ssp*mz**(-2)*t2 + 128*pq*rq*ssz*
     &    ssp*mz**(-2)*u2 + 32*pq*hl**2*ssp*m1**2*t2t**(-1) - 24*pq*
     &    hl**2*ssp*mt**2*t2t**(-1) + 8*pq*hl**2*ssp + 32*pq*hr**2*ssp*
     &    m1**2*t2t**(-1) - 24*pq*hr**2*ssp*mt**2*t2t**(-1) + 8*pq*
     &    hr**2*ssp )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-2,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 64*ssz**2*lq2 + 64*ssz**2*rq2 )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-2,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 64*ssz**2*lq2*m1**2 - 64*ssz**2*lq2*mz**2*s*sz**(-1)
     &     - 128*ssz**2*lq2*s + 64*ssz**2*lq2*s**2*sz**(-1) - 64*ssz**2
     &    *lq2*t2 - 64*ssz**2*lq2*u2 + 64*ssz**2*rq2*m1**2 - 64*ssz**2*
     &    rq2*mz**2*s*sz**(-1) - 128*ssz**2*rq2*s + 64*ssz**2*rq2*s**2*
     &    sz**(-1) - 64*ssz**2*rq2*t2 - 64*ssz**2*rq2*u2 )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-1,0)*Nc*Cf*s4*Pi*alphas*
     & hardfac * ( 128*pq*lq*ssz*ssp*mz**(-2) + 128*pq*rq*ssz*ssp*
     &    mz**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 56*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*m1**2*
     &    mt**2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**2*mz**2*t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*m1**2*u2*t2t**(-1) + 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*m1**2 + 32*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    m1**4*t2t**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**2*mz**2*t2t**(-1) + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2
     &    *ssz*mt**2*u2*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*
     &    hl**2*ssz*mt**2 + 24*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mt**4*t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*
     &    mz**2 - 8*(u2-m12+mt2+mz2)**(-1)*lq*hl**2*ssz*u2 - 56*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mt**2*t2t**(-1)
     &     - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*mz**2*
     &    t2t**(-1) - 32*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2*
     &    u2*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * ( 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**2 + 32*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*m1**4*t2t**(-1) + 24*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*mz**2*t2t**(-1)
     &     + 24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2*u2*
     &    t2t**(-1) - 8*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**2 + 
     &    24*(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mt**4*t2t**(-1) - 8
     &    *(u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*mz**2 - 8*
     &    (u2-m12+mt2+mz2)**(-1)*rq*hr**2*ssz*u2 + 128*pq*lq*ssz*ssp
     &    *m1**2*mz**(-2) - 128*pq*lq*ssz*ssp*mz**(-2)*s - 128*pq*lq*
     &    ssz*ssp*mz**(-2)*t2 - 128*pq*lq*ssz*ssp*mz**(-2)*u2 - 128*pq*
     &    lq*ssz*ssp*mz**2*sz**(-1) + 128*pq*lq*ssz*ssp*s*sz**(-1) - 
     &    128*pq*lq*ssz*ssp + 128*pq*rq*ssz*ssp*m1**2*mz**(-2) - 128*pq
     &    *rq*ssz*ssp*mz**(-2)*s - 128*pq*rq*ssz*ssp*mz**(-2)*t2 - 128*
     &    pq*rq*ssz*ssp*mz**(-2)*u2 - 128*pq*rq*ssz*ssp*mz**2*sz**(-1)
     &     + 128*pq*rq*ssz*ssp*s*sz**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(10,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 128*pq*rq*ssz*ssp + 32*lq*hl**2*ssz*m1**2*
     &    t2t**(-1) - 24*lq*hl**2*ssz*mt**2*t2t**(-1) - 16*lq*hl**2*ssz
     &    *mz**2*s*sz**(-1)*t2t**(-1) - 16*lq*hl**2*ssz*s*t2t**(-1) + 
     &    16*lq*hl**2*ssz*s**2*sz**(-1)*t2t**(-1) + 8*lq*hl**2*ssz + 32
     &    *rq*hr**2*ssz*m1**2*t2t**(-1) - 24*rq*hr**2*ssz*mt**2*
     &    t2t**(-1) - 16*rq*hr**2*ssz*mz**2*s*sz**(-1)*t2t**(-1) - 16*
     &    rq*hr**2*ssz*s*t2t**(-1) + 16*rq*hr**2*ssz*s**2*sz**(-1)*
     &    t2t**(-1) + 8*rq*hr**2*ssz - 64*ssz**2*lq2*mz**2*s*sz**(-2)
     &     - 64*ssz**2*lq2*s*sz**(-1) + 64*ssz**2*lq2*s**2*sz**(-2) - 
     &    64*ssz**2*rq2*mz**2*s*sz**(-2) - 64*ssz**2*rq2*s*sz**(-1) + 
     &    64*ssz**2*rq2*s**2*sz**(-2) )
      MMpartial4 = MMpartial4 + ANGfin(11,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*mh1**2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*sqrt2**(-1)*mt*u2
     &    *t2t**(-1) + 8*(u2-m12+mt2+mh12)**(-1)*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*t2t**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    mt*t2t**(-1) )
      MMpartial4 = MMpartial4 + ANGfin(12,0,-1,0)*Nc*Cf*Pi*alphas*
     & hardfac * (  - 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)
     &    *hl*hr*h2*lambda2*sqrt2**(-1)*mt*mh2**2*t2t**(-1) + 8*
     &    (u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*sqrt2**(-1)*mt*u2
     &    *t2t**(-1) + 8*(u2-m12+mt2+mh22)**(-1)*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*t2t**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt*t2t**(-1) )

c               the phase space except for 1/s**2 
      HH_QBH = MMpartial4 / ( 16.D0 * pi**2 )**2 / 2.D0*s4/(s4+m1**2)

c               the averaging factors
      HH_QBH = HH_QBH /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      HH_QBH = HH_QBH * (m1+m2)**2/4.D0 

      end

      

c ---------------------------------------
      real*8 function DYFACT(mz)
      
      implicit none

      real*8 mz 
      real*8            s,t2,u2 
      common/HH_QBB_INT/s,t2,u2
      
      DYFACT = 1.D0+mz**2*(s+t2)/s/u2 
      DYFACT = 1.D0/DYFACT

      end 

c ---------------------------------------
      real*8 function DYFACU(mz)
      
      implicit none

      real*8 mz 
      real*8            s,t2,u2 
      common/HH_QBB_INT/s,t2,u2
      
      DYFACU = 1.D0+mz**2*(s+u2)/s/t2 
      DYFACU = 1.D0/DYFACU

      end 
































