c --------------------------------------------------------------------
      real*8 function HH_QBB(massin,C)

      implicit none 

      real*8     massin(1:30),Pi,Nc,C(1:20)
     &          ,mz,mh1,mh2,mt
     &          ,ssp,ssz,hl,hr 
     &          ,h1,h2,lambda1,lambda2
     &          ,lq,rq,pq
     &          ,lq2,rq2,pq2
     &          ,s,sz,s1,s2,m1,m2,t2,u2,tx,pt2,sqrt2,epsbw
     &          ,MMb

      Pi= 4.D0 * atan(1.D0)
      sqrt2 = sqrt(2.D0) 

      Nc = 3.D0 

      s     = massin(1)
      t2    = massin(2)
      m1    = massin(6)
      m2    = massin(6)
      mt    = massin(7)
      mz    = massin(8)
      mh1   = massin(9)
      mh2   = massin(10)
      epsbw = massin(26)

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

      sz = s-mz**2
      s1 = s-mh1**2
      s2 = s-mh2**2

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

      tx = t2+m2**2-mt**2
      u2 = -s-t2

      pt2 = ( t2*u2 - s*m2**2 )/s 

c this a copy from the form output 
      MMb =
     &  + pt2*Nc*s**(-1) * ( 8*ssp**2*pq2 )
      MMb = MMb + pt2*Nc*s*tx**(-2) * ( 1./4.*hl**4 + 1./4.*hr**4 )
      MMb = MMb + pt2*Nc*s*tx**(-1)*sz**(-1) * ( 2*lq*hl**2*ssz + 2*rq*
     &    hr**2*ssz )
      MMb = MMb + pt2*Nc*s*sz**(-2) * ( 4*ssz**2*lq2 + 4*ssz**2*rq2 )
      MMb = MMb + pt2*Nc*tx**(-1) * ( 2*pq*hl**2*ssp + 2*pq*hr**2*ssp )
      MMb = MMb + pt2*Nc*sz**(-1) * ( 8*pq*lq*ssz*ssp + 8*pq*rq*ssz*ssp
     &     )
      MMb = MMb + Nc*mt*s*tx**(-1)*s1**(-1) * ( 2*hl*hr*h1*lambda1*
     &    sqrt2**(-1) )
      MMb = MMb + Nc*mt*s*tx**(-1)*s2**(-1) * ( 2*hl*hr*h2*lambda2*
     &    sqrt2**(-1) )
      MMb = MMb + Nc*mt**2*s*tx**(-2) * ( 1./2.*hl**2*hr**2 )
      MMb = MMb + Nc*s*s1**(-2) * ( h1**2*lambda1**2 )
      MMb = MMb + Nc*s*s1**(-1)*s2**(-1) * ( 2*h1*h2*lambda1*lambda2 )
      MMb = MMb + Nc*s*s2**(-2) * ( h2**2*lambda2**2 )
      
c               the phase space except for 1/s**2 
      HH_QBB = MMb / ( 16.D0 * Pi )

c               the averaging factors
      HH_QBB = HH_QBB /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      HH_QBB = HH_QBB * (m1+m2)**2/4.D0

      end

c --------------------------------------------------------------------
c extract the log(Delta) terms from the virtual+soft scaling function 
c everything else copied from HH_QGV 
      real*8 function HH_QBD(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,Cf,alphas,prefac
     &          ,logqf,delta,logdel,log2del,bornfunc
     &          ,s,m1
     &          ,t1,u1,t2,u2,MMd,s4,s4p,HH_QBB

      Pi     = 4.D0 * atan(1.D0)
      prefac = 1.D0/(16.D0*Pi**2)

      Cf = 4.D0/3.D0

      s     = massin(1)
      t2    = massin(2)
      s4    = massin(3)
      m1    = massin(6)
      delta = massin(20)
      s4p   = massin(21)

      u2 = -s-t2

      t1 = t2 
      u1 = u2 

c               the logaritms for linear s4 integration 
      logdel  =  log(s4p/m1**2)/(s4p-delta) - 1.D0/s4
      log2del =  log(s4p/m1**2)**2/(s4p-delta) 
     &          - 2.D0*log(s4/m1**2)/s4

c               the factorization/renormalization scale 
      logqf = log( massin(13)**2/m1**2 )

      bornfunc = HH_QBB(massin,C)

c               set gs=1 
      alphas = 1/(4.D0*Pi)

c               insert the form output
      MMd = 0.D0

      MMd = MMd + Cf*Pi*alphas*prefac*bornfunc*logdel * (
     &     - 32*logqf
     &     + 32*log(m1**(-2)*s)
     &     - 32*log( - m1**(-2)*t1)
     &     - 32*log( - m1**(-2)*u1)
     &     )
      MMd = MMd + Cf*Pi*alphas*prefac*bornfunc*log2del * (
     &     + 32
     &     )

      HH_QBD = MMd 

      end

c --------------------------------------------------------------------
      real*8 function HH_QBV(massin,C)

      implicit none 

      real*8     massin(1:30),C(1:20),Pi,sqrt2,Nc,Cf,zeta2,alphas
     &          ,logqr,logqf,kaellen,prefac,epsbw
     &          ,h1,h2,hl,hr,lambda1,lambda2
     &          ,lq,lq2,rq,rq2,pq,pq2,ssp,ssz
     &          ,m1,m2,mt,mz,mh1,mh2
     &          ,s,t,tx,t1,t2,u1,u2,s1,s2,sz,pt2
     &          ,SCA(1:10),SCB(1:10,1:6),SCBP(10)
     &          ,SCC(1:20,1:9),SCD(1:10,1:8)
     &          ,MMv,Li2
      complex*16 CSPEN

      external CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

      Pi     = 4.D0 * atan(1.D0)
      sqrt2  = sqrt(2.D0) 
      prefac = 1.D0/(16.D0*Pi**2)
      zeta2  = Pi**2/6.D0

      Nc = 3.D0 
      Cf = 4.D0/3.D0

      s     = massin(1)
      t2    = massin(2)
      m1    = massin(6)
      m2    = massin(6)
      mt    = massin(7)
      mz    = massin(8)
      mh1   = massin(9)
      mh2   = massin(10)
      epsbw = massin(26)

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

      sz = s-mz**2
      s1 = s-mh1**2
      s2 = s-mh2**2

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

      u2 = -s-t2
      tx = t2+m2**2-mt**2

      t1 = t2 
      u1 = u2 
      t  = t2+m2**2

      pt2 = ( t2*u2 - s*m2**2 )/s 
      kaellen = s**2 +m1**4 +m2**4 - 2*( s*m1**2+s*m2**2+m1**2*m2**2 )

c               the factorization/renormalization scale 
      logqr = log( massin(12)**2/m1**2 )
      logqf = log( massin(13)**2/m1**2 )

c               the scalar functions 
      call SCALAR_ARRAY_HH(massin,SCA,SCB,SCBP,SCC,SCD)

c               set gs=1 
      alphas = 1/(4.D0*Pi)

      MMv = 0.d0
c               insert the form output 
c               insert abs into log(-1+mt^2/m1^2)
      MMv =
     &  + zeta2*Nc*Cf*Pi*alphas*prefac * (  - 192*pq*lq*ssz*ssp*pt2*
     &    sz**(-1) - 192*pq*rq*ssz*ssp*pt2*sz**(-1) - 48*pq*hl**2*ssp*
     &    pt2*tx**(-1) - 48*pq*hr**2*ssp*pt2*tx**(-1) - 48*lq*hl**2*ssz
     &    *pt2*s*tx**(-1)*sz**(-1) - 48*rq*hr**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) - 48*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s1**(-1) - 48*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s2**(-1) - 12*hl**2*hr**2*mt**2*s*tx**(-2) - 6*hl**4*pt2*s*
     &    tx**(-2) - 6*hr**4*pt2*s*tx**(-2) - 48*h1*h2*lambda1*lambda2*
     &    s*s1**(-1)*s2**(-1) - 24*h1**2*lambda1**2*s*s1**(-2) - 24*
     &    h2**2*lambda2**2*s*s2**(-2) - 96*ssz**2*lq2*pt2*s*sz**(-2) - 
     &    96*ssz**2*rq2*pt2*s*sz**(-2) - 192*ssp**2*pq2*pt2*s**(-1) )
      MMv = MMv + Nc*Cf*Pi*alphas*prefac*logqf * (  - 192*pq*lq*ssz*ssp
     &    *pt2*sz**(-1) - 192*pq*rq*ssz*ssp*pt2*sz**(-1) - 48*pq*hl**2*
     &    ssp*pt2*tx**(-1) - 48*pq*hr**2*ssp*pt2*tx**(-1) - 48*lq*hl**2
     &    *ssz*pt2*s*tx**(-1)*sz**(-1) - 48*rq*hr**2*ssz*pt2*s*tx**(-1)
     &    *sz**(-1) - 48*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s1**(-1) - 48*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s2**(-1) - 12*hl**2*hr**2*mt**2*s*tx**(-2) - 6*hl**4*pt2*s*
     &    tx**(-2) - 6*hr**4*pt2*s*tx**(-2) - 48*h1*h2*lambda1*lambda2*
     &    s*s1**(-1)*s2**(-1) - 24*h1**2*lambda1**2*s*s1**(-2) - 24*
     &    h2**2*lambda2**2*s*s2**(-2) - 96*ssz**2*lq2*pt2*s*sz**(-2) - 
     &    96*ssz**2*rq2*pt2*s*sz**(-2) - 192*ssp**2*pq2*pt2*s**(-1) )
      MMv = MMv + Nc*Cf*Pi*alphas*prefac * ( 128*pq*lq*ssz*ssp*m1**2*
     &    sz**(-1) - 128*pq*lq*ssz*ssp*s**(-1)*t1*u1*sz**(-1) + 128*pq*
     &    rq*ssz*ssp*m1**2*sz**(-1) - 128*pq*rq*ssz*ssp*s**(-1)*t1*u1*
     &    sz**(-1) + 8*pq*hl**2*ssp*pt2*m1**2*tx**(-2) - 8*pq*hl**2*ssp
     &    *pt2*mt**2*tx**(-2) + 8*pq*hl**2*ssp*pt2*t1*tx**(-2) + 48*pq*
     &    hl**2*ssp*pt2*tx**(-1)*logqr + 48*pq*hl**2*ssp*m1**2*tx**(-1)
     &     - 8*pq*hl**2*ssp*s**(-1)*t1**(-2)*u1**3 - 48*pq*hl**2*ssp*
     &    s**(-1)*t1*u1*tx**(-1) + 8*pq*hl**2*ssp*s**(-1)*u1 - 8*pq*
     &    hl**2*ssp*s*t1**(-2)*u1 - 16*pq*hl**2*ssp*t1**(-2)*u1**2 + 8*
     &    pq*hr**2*ssp*pt2*m1**2*tx**(-2) - 8*pq*hr**2*ssp*pt2*mt**2*
     &    tx**(-2) + 8*pq*hr**2*ssp*pt2*t1*tx**(-2) + 48*pq*hr**2*ssp*
     &    pt2*tx**(-1)*logqr + 48*pq*hr**2*ssp*m1**2*tx**(-1) - 8*pq*
     &    hr**2*ssp*s**(-1)*t1**(-2)*u1**3 - 48*pq*hr**2*ssp*s**(-1)*t1
     &    *u1*tx**(-1) + 8*pq*hr**2*ssp*s**(-1)*u1 - 8*pq*hr**2*ssp*s*
     &    t1**(-2)*u1 - 16*pq*hr**2*ssp*t1**(-2)*u1**2 + 8*lq*hl**2*ssz
     &    *pt2*m1**2*s*tx**(-2)*sz**(-1) )
      MMv = MMv + Nc*Cf*Pi*alphas*prefac * (  - 8*lq*hl**2*ssz*pt2*
     &    mt**2*s*tx**(-2)*sz**(-1) + 8*lq*hl**2*ssz*pt2*s*t1*tx**(-2)*
     &    sz**(-1) + 48*lq*hl**2*ssz*pt2*s*tx**(-1)*sz**(-1)*logqr + 48
     &    *lq*hl**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*s*
     &    t1**(-2)*u1**2*sz**(-1) - 8*lq*hl**2*ssz*s**2*t1**(-2)*u1*
     &    sz**(-1) - 8*lq*hl**2*ssz*t1**(-2)*u1**3*sz**(-1) - 48*lq*
     &    hl**2*ssz*t1*u1*tx**(-1)*sz**(-1) + 8*lq*hl**2*ssz*u1*
     &    sz**(-1) + 8*rq*hr**2*ssz*pt2*m1**2*s*tx**(-2)*sz**(-1) - 8*
     &    rq*hr**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz*
     &    pt2*s*t1*tx**(-2)*sz**(-1) + 48*rq*hr**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1)*logqr + 48*rq*hr**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 
     &    16*rq*hr**2*ssz*s*t1**(-2)*u1**2*sz**(-1) - 8*rq*hr**2*ssz*
     &    s**2*t1**(-2)*u1*sz**(-1) - 8*rq*hr**2*ssz*t1**(-2)*u1**3*
     &    sz**(-1) - 48*rq*hr**2*ssz*t1*u1*tx**(-1)*sz**(-1) + 8*rq*
     &    hr**2*ssz*u1*sz**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*
     &    mt*s*tx**(-2)*s1**(-1) )
      MMv = MMv + Nc*Cf*Pi*alphas*prefac * ( 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*t1*tx**(-2)*s1**(-1) + 72*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1)*logqr - 48*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) - 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s*tx**(-2)*s1**(-1) + 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s2**(-1) + 8*hl*hr*h2
     &    *lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) + 72*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1)*logqr - 48*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt**3*s*tx**(-2)*s2**(-1) + 4*hl**2*hr**2
     &    *m1**2*mt**2*s*tx**(-3) + 4*hl**2*hr**2*mt**2*s*t1*tx**(-3)
     &     + 24*hl**2*hr**2*mt**2*s*tx**(-2)*logqr - 16*hl**2*hr**2*
     &    mt**2*s*tx**(-2) - 4*hl**2*hr**2*mt**4*s*tx**(-3) + 2*hl**4*
     &    pt2*m1**2*s*tx**(-3) - 2*hl**4*pt2*mt**2*s*tx**(-3) + 2*hl**4
     &    *pt2*s*t1*tx**(-3) + 12*hl**4*pt2*s*tx**(-2)*logqr + 8*hl**4*
     &    m1**2*s*tx**(-2) )
      MMv = MMv + Nc*Cf*Pi*alphas*prefac * (  - 4*hl**4*s*t1**(-2)*
     &    u1**2*tx**(-1) - 2*hl**4*s**2*t1**(-2)*u1*tx**(-1) - 2*hl**4*
     &    t1**(-2)*u1**3*tx**(-1) - 8*hl**4*t1*u1*tx**(-2) + 2*hl**4*u1
     &    *tx**(-1) + 2*hr**4*pt2*m1**2*s*tx**(-3) - 2*hr**4*pt2*mt**2*
     &    s*tx**(-3) + 2*hr**4*pt2*s*t1*tx**(-3) + 12*hr**4*pt2*s*
     &    tx**(-2)*logqr + 8*hr**4*m1**2*s*tx**(-2) - 4*hr**4*s*
     &    t1**(-2)*u1**2*tx**(-1) - 2*hr**4*s**2*t1**(-2)*u1*tx**(-1)
     &     - 2*hr**4*t1**(-2)*u1**3*tx**(-1) - 8*hr**4*t1*u1*tx**(-2)
     &     + 2*hr**4*u1*tx**(-1) + 48*h1*h2*lambda1*lambda2*s*s1**(-1)*
     &    s2**(-1)*logqr - 32*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1)
     &     + 24*h1**2*lambda1**2*s*s1**(-2)*logqr - 16*h1**2*lambda1**2
     &    *s*s1**(-2) + 24*h2**2*lambda2**2*s*s2**(-2)*logqr - 16*h2**2
     &    *lambda2**2*s*s2**(-2) + 64*ssz**2*lq2*m1**2*s*sz**(-2) - 64*
     &    ssz**2*lq2*t1*u1*sz**(-2) + 64*ssz**2*rq2*m1**2*s*sz**(-2) - 
     &    64*ssz**2*rq2*t1*u1*sz**(-2) + 128*ssp**2*pq2*m1**2*s**(-1)
     &     - 128*ssp**2*pq2*s**(-2)*t1*u1 )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * ( 32*pq*hl**2*ssp*m1**2*mt**2*s**(-1)*tx**(-1) - 8*pq*hl**2*
     &    ssp*m1**2*s**(-1)*t1**(-3)*u1**3 - 32*pq*hl**2*ssp*m1**2*
     &    s**(-1)*t1*tx**(-1) + 8*pq*hl**2*ssp*m1**2*s**(-1) - 8*pq*
     &    hl**2*ssp*m1**2*s*t1**(-3)*u1 - 16*pq*hl**2*ssp*m1**2*
     &    t1**(-3)*u1**2 - 8*pq*hl**2*ssp*m1**2*t1**(-1) - 16*pq*hl**2*
     &    ssp*m1**4*s**(-1)*tx**(-1) + 8*pq*hl**2*ssp*mt**2*s**(-1)*
     &    t1**(-3)*u1**3 + 32*pq*hl**2*ssp*mt**2*s**(-1)*t1*tx**(-1) - 
     &    8*pq*hl**2*ssp*mt**2*s**(-1) + 8*pq*hl**2*ssp*mt**2*s*
     &    t1**(-3)*u1 + 16*pq*hl**2*ssp*mt**2*t1**(-3)*u1**2 + 8*pq*
     &    hl**2*ssp*mt**2*t1**(-1) - 16*pq*hl**2*ssp*mt**4*s**(-1)*
     &    tx**(-1) - 8*pq*hl**2*ssp*s**(-1)*t1**(-3)*u1**4 - 8*pq*hl**2
     &    *ssp*s**(-1)*t1**(-2)*u1**3 + 8*pq*hl**2*ssp*s**(-1)*t1**(-1)
     &    *u1**2 + 8*pq*hl**2*ssp*s**(-1)*t1 - 16*pq*hl**2*ssp*s**(-1)*
     &    t1**2*tx**(-1) - 24*pq*hl**2*ssp*s*t1**(-3)*u1**2 - 8*pq*
     &    hl**2*ssp*s*t1**(-2)*u1 )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * (  - 8*pq*hl**2*ssp*s**2*t1**(-3)*u1 - 24*pq*hl**2*ssp*
     &    t1**(-3)*u1**3 - 16*pq*hl**2*ssp*t1**(-2)*u1**2 + 8*pq*hl**2*
     &    ssp*t1**(-1)*u1 - 8*pq*hl**2*ssp + 32*pq*hr**2*ssp*m1**2*
     &    mt**2*s**(-1)*tx**(-1) - 8*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t1**(-3)*u1**3 - 32*pq*hr**2*ssp*m1**2*s**(-1)*t1*tx**(-1) + 
     &    8*pq*hr**2*ssp*m1**2*s**(-1) - 8*pq*hr**2*ssp*m1**2*s*
     &    t1**(-3)*u1 - 16*pq*hr**2*ssp*m1**2*t1**(-3)*u1**2 - 8*pq*
     &    hr**2*ssp*m1**2*t1**(-1) - 16*pq*hr**2*ssp*m1**4*s**(-1)*
     &    tx**(-1) + 8*pq*hr**2*ssp*mt**2*s**(-1)*t1**(-3)*u1**3 + 32*
     &    pq*hr**2*ssp*mt**2*s**(-1)*t1*tx**(-1) - 8*pq*hr**2*ssp*mt**2
     &    *s**(-1) + 8*pq*hr**2*ssp*mt**2*s*t1**(-3)*u1 + 16*pq*hr**2*
     &    ssp*mt**2*t1**(-3)*u1**2 + 8*pq*hr**2*ssp*mt**2*t1**(-1) - 16
     &    *pq*hr**2*ssp*mt**4*s**(-1)*tx**(-1) - 8*pq*hr**2*ssp*s**(-1)
     &    *t1**(-3)*u1**4 - 8*pq*hr**2*ssp*s**(-1)*t1**(-2)*u1**3 + 8*
     &    pq*hr**2*ssp*s**(-1)*t1**(-1)*u1**2 )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * ( 8*pq*hr**2*ssp*s**(-1)*t1 - 16*pq*hr**2*ssp*s**(-1)*t1**2*
     &    tx**(-1) - 24*pq*hr**2*ssp*s*t1**(-3)*u1**2 - 8*pq*hr**2*ssp*
     &    s*t1**(-2)*u1 - 8*pq*hr**2*ssp*s**2*t1**(-3)*u1 - 24*pq*hr**2
     &    *ssp*t1**(-3)*u1**3 - 16*pq*hr**2*ssp*t1**(-2)*u1**2 + 8*pq*
     &    hr**2*ssp*t1**(-1)*u1 - 8*pq*hr**2*ssp + 32*lq*hl**2*ssz*
     &    m1**2*mt**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**2*s*
     &    t1**(-3)*u1**2*sz**(-1) - 8*lq*hl**2*ssz*m1**2*s*t1**(-1)*
     &    sz**(-1) - 8*lq*hl**2*ssz*m1**2*s**2*t1**(-3)*u1*sz**(-1) - 8
     &    *lq*hl**2*ssz*m1**2*t1**(-3)*u1**3*sz**(-1) - 32*lq*hl**2*ssz
     &    *m1**2*t1*tx**(-1)*sz**(-1) + 8*lq*hl**2*ssz*m1**2*sz**(-1)
     &     - 16*lq*hl**2*ssz*m1**4*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*
     &    mt**2*s*t1**(-3)*u1**2*sz**(-1) + 8*lq*hl**2*ssz*mt**2*s*
     &    t1**(-1)*sz**(-1) + 8*lq*hl**2*ssz*mt**2*s**2*t1**(-3)*u1*
     &    sz**(-1) + 8*lq*hl**2*ssz*mt**2*t1**(-3)*u1**3*sz**(-1) + 32*
     &    lq*hl**2*ssz*mt**2*t1*tx**(-1)*sz**(-1) )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * (  - 8*lq*hl**2*ssz*mt**2*sz**(-1) - 16*lq*hl**2*ssz*mt**4*
     &    tx**(-1)*sz**(-1) - 24*lq*hl**2*ssz*s*t1**(-3)*u1**3*sz**(-1)
     &     - 16*lq*hl**2*ssz*s*t1**(-2)*u1**2*sz**(-1) + 8*lq*hl**2*ssz
     &    *s*t1**(-1)*u1*sz**(-1) - 8*lq*hl**2*ssz*s*sz**(-1) - 24*lq*
     &    hl**2*ssz*s**2*t1**(-3)*u1**2*sz**(-1) - 8*lq*hl**2*ssz*s**2*
     &    t1**(-2)*u1*sz**(-1) - 8*lq*hl**2*ssz*s**3*t1**(-3)*u1*
     &    sz**(-1) - 8*lq*hl**2*ssz*t1**(-3)*u1**4*sz**(-1) - 8*lq*
     &    hl**2*ssz*t1**(-2)*u1**3*sz**(-1) + 8*lq*hl**2*ssz*t1**(-1)*
     &    u1**2*sz**(-1) + 8*lq*hl**2*ssz*t1*sz**(-1) - 16*lq*hl**2*ssz
     &    *t1**2*tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*m1**2*mt**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*m1**2*s*t1**(-3)*u1**2*
     &    sz**(-1) - 8*rq*hr**2*ssz*m1**2*s*t1**(-1)*sz**(-1) - 8*rq*
     &    hr**2*ssz*m1**2*s**2*t1**(-3)*u1*sz**(-1) - 8*rq*hr**2*ssz*
     &    m1**2*t1**(-3)*u1**3*sz**(-1) - 32*rq*hr**2*ssz*m1**2*t1*
     &    tx**(-1)*sz**(-1) )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * ( 8*rq*hr**2*ssz*m1**2*sz**(-1) - 16*rq*hr**2*ssz*m1**4*
     &    tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*mt**2*s*t1**(-3)*u1**2*
     &    sz**(-1) + 8*rq*hr**2*ssz*mt**2*s*t1**(-1)*sz**(-1) + 8*rq*
     &    hr**2*ssz*mt**2*s**2*t1**(-3)*u1*sz**(-1) + 8*rq*hr**2*ssz*
     &    mt**2*t1**(-3)*u1**3*sz**(-1) + 32*rq*hr**2*ssz*mt**2*t1*
     &    tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*mt**2*sz**(-1) - 16*rq*
     &    hr**2*ssz*mt**4*tx**(-1)*sz**(-1) - 24*rq*hr**2*ssz*s*
     &    t1**(-3)*u1**3*sz**(-1) - 16*rq*hr**2*ssz*s*t1**(-2)*u1**2*
     &    sz**(-1) + 8*rq*hr**2*ssz*s*t1**(-1)*u1*sz**(-1) - 8*rq*hr**2
     &    *ssz*s*sz**(-1) - 24*rq*hr**2*ssz*s**2*t1**(-3)*u1**2*
     &    sz**(-1) - 8*rq*hr**2*ssz*s**2*t1**(-2)*u1*sz**(-1) - 8*rq*
     &    hr**2*ssz*s**3*t1**(-3)*u1*sz**(-1) - 8*rq*hr**2*ssz*t1**(-3)
     &    *u1**4*sz**(-1) - 8*rq*hr**2*ssz*t1**(-2)*u1**3*sz**(-1) + 8*
     &    rq*hr**2*ssz*t1**(-1)*u1**2*sz**(-1) + 8*rq*hr**2*ssz*t1*
     &    sz**(-1) )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * (  - 16*rq*hr**2*ssz*t1**2*tx**(-1)*sz**(-1) + 8*hl**4*m1**2*
     &    mt**2*tx**(-2) - 8*hl**4*m1**2*t1*tx**(-2) + 2*hl**4*m1**2*
     &    tx**(-1) - 4*hl**4*m1**4*tx**(-2) + 8*hl**4*mt**2*t1*tx**(-2)
     &     - 2*hl**4*mt**2*tx**(-1) - 4*hl**4*mt**4*tx**(-2) - 4*hl**4*
     &    s*t1**(-3)*u1**2 - 6*hl**4*s*t1**(-3)*u1**3*tx**(-1) + 2*
     &    hl**4*s*t1**(-1)*u1*tx**(-1) - 2*hl**4*s*t1**(-1) - 2*hl**4*
     &    s**2*t1**(-3)*u1 - 6*hl**4*s**2*t1**(-3)*u1**2*tx**(-1) - 2*
     &    hl**4*s**3*t1**(-3)*u1*tx**(-1) - 2*hl**4*t1**(-3)*u1**3 - 2*
     &    hl**4*t1**(-3)*u1**4*tx**(-1) + 2*hl**4*t1**(-1)*u1**2*
     &    tx**(-1) + 2*hl**4*t1*tx**(-1) - 4*hl**4*t1**2*tx**(-2) + 8*
     &    hr**4*m1**2*mt**2*tx**(-2) - 8*hr**4*m1**2*t1*tx**(-2) + 2*
     &    hr**4*m1**2*tx**(-1) - 4*hr**4*m1**4*tx**(-2) + 8*hr**4*mt**2
     &    *t1*tx**(-2) - 2*hr**4*mt**2*tx**(-1) - 4*hr**4*mt**4*
     &    tx**(-2) - 4*hr**4*s*t1**(-3)*u1**2 - 6*hr**4*s*t1**(-3)*
     &    u1**3*tx**(-1) )
      MMv = MMv + log(abs(-1+m1**(-2)*mt**2))*Nc*Cf*Pi*alphas*prefac
     &  * ( 2*hr**4*s*t1**(-1)*u1*tx**(-1) - 2*hr**4*s*t1**(-1) - 2*
     &    hr**4*s**2*t1**(-3)*u1 - 6*hr**4*s**2*t1**(-3)*u1**2*tx**(-1)
     &     - 2*hr**4*s**3*t1**(-3)*u1*tx**(-1) - 2*hr**4*t1**(-3)*u1**3
     &     - 2*hr**4*t1**(-3)*u1**4*tx**(-1) + 2*hr**4*t1**(-1)*u1**2*
     &    tx**(-1) + 2*hr**4*t1*tx**(-1) - 4*hr**4*t1**2*tx**(-2) )
      MMv = MMv + log(m1**(-2)*s)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*
     &    hl**2*ssp*m1**2*mt**2*s**(-1)*tx**(-1) + 16*pq*hl**2*ssp*
     &    m1**2*s**(-1)*t1*tx**(-1) - 8*pq*hl**2*ssp*m1**2*s**(-1) + 8*
     &    pq*hl**2*ssp*m1**4*s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*mt**2*
     &    s**(-1)*t1*tx**(-1) + 8*pq*hl**2*ssp*mt**2*s**(-1) + 8*pq*
     &    hl**2*ssp*mt**4*s**(-1)*tx**(-1) - 8*pq*hl**2*ssp*s**(-1)*t1
     &     + 8*pq*hl**2*ssp*s**(-1)*t1**2*tx**(-1) - 16*pq*hr**2*ssp*
     &    m1**2*mt**2*s**(-1)*tx**(-1) + 16*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t1*tx**(-1) - 8*pq*hr**2*ssp*m1**2*s**(-1) + 8*pq*hr**2*ssp*
     &    m1**4*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp*mt**2*s**(-1)*t1*
     &    tx**(-1) + 8*pq*hr**2*ssp*mt**2*s**(-1) + 8*pq*hr**2*ssp*
     &    mt**4*s**(-1)*tx**(-1) - 8*pq*hr**2*ssp*s**(-1)*t1 + 8*pq*
     &    hr**2*ssp*s**(-1)*t1**2*tx**(-1) - 16*lq*hl**2*ssz*m1**2*
     &    mt**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*m1**2*t1*tx**(-1)*
     &    sz**(-1) - 8*lq*hl**2*ssz*m1**2*sz**(-1) + 8*lq*hl**2*ssz*
     &    m1**4*tx**(-1)*sz**(-1) )
      MMv = MMv + log(m1**(-2)*s)*Nc*Cf*Pi*alphas*prefac * (  - 16*lq*
     &    hl**2*ssz*mt**2*t1*tx**(-1)*sz**(-1) + 8*lq*hl**2*ssz*mt**2*
     &    sz**(-1) + 8*lq*hl**2*ssz*mt**4*tx**(-1)*sz**(-1) - 8*lq*
     &    hl**2*ssz*t1*sz**(-1) + 8*lq*hl**2*ssz*t1**2*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*m1**2*mt**2*tx**(-1)*sz**(-1) + 16
     &    *rq*hr**2*ssz*m1**2*t1*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*
     &    m1**2*sz**(-1) + 8*rq*hr**2*ssz*m1**4*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*mt**2*t1*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*
     &    mt**2*sz**(-1) + 8*rq*hr**2*ssz*mt**4*tx**(-1)*sz**(-1) - 8*
     &    rq*hr**2*ssz*t1*sz**(-1) + 8*rq*hr**2*ssz*t1**2*tx**(-1)*
     &    sz**(-1) - 4*hl**4*m1**2*mt**2*tx**(-2) + 4*hl**4*m1**2*t1*
     &    tx**(-2) - 2*hl**4*m1**2*tx**(-1) + 2*hl**4*m1**4*tx**(-2) - 
     &    4*hl**4*mt**2*t1*tx**(-2) + 2*hl**4*mt**2*tx**(-1) + 2*hl**4*
     &    mt**4*tx**(-2) - 2*hl**4*t1*tx**(-1) + 2*hl**4*t1**2*tx**(-2)
     &     - 4*hr**4*m1**2*mt**2*tx**(-2) + 4*hr**4*m1**2*t1*tx**(-2)
     &     - 2*hr**4*m1**2*tx**(-1) )
      MMv = MMv + log(m1**(-2)*s)*Nc*Cf*Pi*alphas*prefac * ( 2*hr**4*
     &    m1**4*tx**(-2) - 4*hr**4*mt**2*t1*tx**(-2) + 2*hr**4*mt**2*
     &    tx**(-1) + 2*hr**4*mt**4*tx**(-2) - 2*hr**4*t1*tx**(-1) + 2*
     &    hr**4*t1**2*tx**(-2) )
      MMv = MMv + log( - m1**(-2)*t1)*Nc*Cf*Pi*alphas*prefac*logqf
     &  * ( 128*pq*lq*ssz*ssp*pt2*sz**(-1) + 128*pq*rq*ssz*ssp*pt2*
     &    sz**(-1) + 32*pq*hl**2*ssp*pt2*tx**(-1) + 32*pq*hr**2*ssp*pt2
     &    *tx**(-1) + 32*lq*hl**2*ssz*pt2*s*tx**(-1)*sz**(-1) + 32*rq*
     &    hr**2*ssz*pt2*s*tx**(-1)*sz**(-1) + 32*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) + 32*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) + 8*hl**2*hr**2*mt**2*s*
     &    tx**(-2) + 4*hl**4*pt2*s*tx**(-2) + 4*hr**4*pt2*s*tx**(-2) + 
     &    32*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1) + 16*h1**2*
     &    lambda1**2*s*s1**(-2) + 16*h2**2*lambda2**2*s*s2**(-2) + 64*
     &    ssz**2*lq2*pt2*s*sz**(-2) + 64*ssz**2*rq2*pt2*s*sz**(-2) + 
     &    128*ssp**2*pq2*pt2*s**(-1) )
      MMv = MMv + log( - m1**(-2)*u1)*Nc*Cf*Pi*alphas*prefac*logqf
     &  * ( 128*pq*lq*ssz*ssp*pt2*sz**(-1) + 128*pq*rq*ssz*ssp*pt2*
     &    sz**(-1) + 32*pq*hl**2*ssp*pt2*tx**(-1) + 32*pq*hr**2*ssp*pt2
     &    *tx**(-1) + 32*lq*hl**2*ssz*pt2*s*tx**(-1)*sz**(-1) + 32*rq*
     &    hr**2*ssz*pt2*s*tx**(-1)*sz**(-1) + 32*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) + 32*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) + 8*hl**2*hr**2*mt**2*s*
     &    tx**(-2) + 4*hl**4*pt2*s*tx**(-2) + 4*hr**4*pt2*s*tx**(-2) + 
     &    32*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1) + 16*h1**2*
     &    lambda1**2*s*s1**(-2) + 16*h2**2*lambda2**2*s*s2**(-2) + 64*
     &    ssz**2*lq2*pt2*s*sz**(-2) + 64*ssz**2*rq2*pt2*s*sz**(-2) + 
     &    128*ssp**2*pq2*pt2*s**(-1) )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * (  - 32
     &    *pq*hl**2*ssp*m1**2*mt**2*s**(-1)*tx**(-1) + 8*pq*hl**2*ssp*
     &    m1**2*s**(-1)*t1**(-3)*u1**3 + 32*pq*hl**2*ssp*m1**2*s**(-1)*
     &    t1*tx**(-1) - 8*pq*hl**2*ssp*m1**2*s**(-1) + 8*pq*hl**2*ssp*
     &    m1**2*s*t1**(-3)*u1 + 16*pq*hl**2*ssp*m1**2*t1**(-3)*u1**2 + 
     &    8*pq*hl**2*ssp*m1**2*t1**(-1) + 16*pq*hl**2*ssp*m1**4*s**(-1)
     &    *tx**(-1) - 8*pq*hl**2*ssp*mt**2*s**(-1)*t1**(-3)*u1**3 - 32*
     &    pq*hl**2*ssp*mt**2*s**(-1)*t1*tx**(-1) + 8*pq*hl**2*ssp*mt**2
     &    *s**(-1) - 8*pq*hl**2*ssp*mt**2*s*t1**(-3)*u1 - 16*pq*hl**2*
     &    ssp*mt**2*t1**(-3)*u1**2 - 8*pq*hl**2*ssp*mt**2*t1**(-1) + 16
     &    *pq*hl**2*ssp*mt**4*s**(-1)*tx**(-1) + 8*pq*hl**2*ssp*s**(-1)
     &    *t1**(-3)*u1**4 + 8*pq*hl**2*ssp*s**(-1)*t1**(-2)*u1**3 - 8*
     &    pq*hl**2*ssp*s**(-1)*t1**(-1)*u1**2 - 8*pq*hl**2*ssp*s**(-1)*
     &    t1 + 16*pq*hl**2*ssp*s**(-1)*t1**2*tx**(-1) + 24*pq*hl**2*ssp
     &    *s*t1**(-3)*u1**2 + 8*pq*hl**2*ssp*s*t1**(-2)*u1 + 8*pq*hl**2
     &    *ssp*s**2*t1**(-3)*u1 )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * ( 24*pq
     &    *hl**2*ssp*t1**(-3)*u1**3 + 16*pq*hl**2*ssp*t1**(-2)*u1**2 - 
     &    8*pq*hl**2*ssp*t1**(-1)*u1 + 8*pq*hl**2*ssp - 32*pq*hr**2*ssp
     &    *m1**2*mt**2*s**(-1)*tx**(-1) + 8*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t1**(-3)*u1**3 + 32*pq*hr**2*ssp*m1**2*s**(-1)*t1*tx**(-1) - 
     &    8*pq*hr**2*ssp*m1**2*s**(-1) + 8*pq*hr**2*ssp*m1**2*s*
     &    t1**(-3)*u1 + 16*pq*hr**2*ssp*m1**2*t1**(-3)*u1**2 + 8*pq*
     &    hr**2*ssp*m1**2*t1**(-1) + 16*pq*hr**2*ssp*m1**4*s**(-1)*
     &    tx**(-1) - 8*pq*hr**2*ssp*mt**2*s**(-1)*t1**(-3)*u1**3 - 32*
     &    pq*hr**2*ssp*mt**2*s**(-1)*t1*tx**(-1) + 8*pq*hr**2*ssp*mt**2
     &    *s**(-1) - 8*pq*hr**2*ssp*mt**2*s*t1**(-3)*u1 - 16*pq*hr**2*
     &    ssp*mt**2*t1**(-3)*u1**2 - 8*pq*hr**2*ssp*mt**2*t1**(-1) + 16
     &    *pq*hr**2*ssp*mt**4*s**(-1)*tx**(-1) + 8*pq*hr**2*ssp*s**(-1)
     &    *t1**(-3)*u1**4 + 8*pq*hr**2*ssp*s**(-1)*t1**(-2)*u1**3 - 8*
     &    pq*hr**2*ssp*s**(-1)*t1**(-1)*u1**2 - 8*pq*hr**2*ssp*s**(-1)*
     &    t1 )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * ( 16*pq
     &    *hr**2*ssp*s**(-1)*t1**2*tx**(-1) + 24*pq*hr**2*ssp*s*
     &    t1**(-3)*u1**2 + 8*pq*hr**2*ssp*s*t1**(-2)*u1 + 8*pq*hr**2*
     &    ssp*s**2*t1**(-3)*u1 + 24*pq*hr**2*ssp*t1**(-3)*u1**3 + 16*pq
     &    *hr**2*ssp*t1**(-2)*u1**2 - 8*pq*hr**2*ssp*t1**(-1)*u1 + 8*pq
     &    *hr**2*ssp - 32*lq*hl**2*ssz*m1**2*mt**2*tx**(-1)*sz**(-1) + 
     &    16*lq*hl**2*ssz*m1**2*s*t1**(-3)*u1**2*sz**(-1) + 8*lq*hl**2*
     &    ssz*m1**2*s*t1**(-1)*sz**(-1) + 8*lq*hl**2*ssz*m1**2*s**2*
     &    t1**(-3)*u1*sz**(-1) + 8*lq*hl**2*ssz*m1**2*t1**(-3)*u1**3*
     &    sz**(-1) + 32*lq*hl**2*ssz*m1**2*t1*tx**(-1)*sz**(-1) - 8*lq*
     &    hl**2*ssz*m1**2*sz**(-1) + 16*lq*hl**2*ssz*m1**4*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*mt**2*s*t1**(-3)*u1**2*sz**(-1) - 
     &    8*lq*hl**2*ssz*mt**2*s*t1**(-1)*sz**(-1) - 8*lq*hl**2*ssz*
     &    mt**2*s**2*t1**(-3)*u1*sz**(-1) - 8*lq*hl**2*ssz*mt**2*
     &    t1**(-3)*u1**3*sz**(-1) - 32*lq*hl**2*ssz*mt**2*t1*tx**(-1)*
     &    sz**(-1) )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * ( 8*lq*
     &    hl**2*ssz*mt**2*sz**(-1) + 16*lq*hl**2*ssz*mt**4*tx**(-1)*
     &    sz**(-1) + 24*lq*hl**2*ssz*s*t1**(-3)*u1**3*sz**(-1) + 16*lq*
     &    hl**2*ssz*s*t1**(-2)*u1**2*sz**(-1) - 8*lq*hl**2*ssz*s*
     &    t1**(-1)*u1*sz**(-1) + 8*lq*hl**2*ssz*s*sz**(-1) + 24*lq*
     &    hl**2*ssz*s**2*t1**(-3)*u1**2*sz**(-1) + 8*lq*hl**2*ssz*s**2*
     &    t1**(-2)*u1*sz**(-1) + 8*lq*hl**2*ssz*s**3*t1**(-3)*u1*
     &    sz**(-1) + 8*lq*hl**2*ssz*t1**(-3)*u1**4*sz**(-1) + 8*lq*
     &    hl**2*ssz*t1**(-2)*u1**3*sz**(-1) - 8*lq*hl**2*ssz*t1**(-1)*
     &    u1**2*sz**(-1) - 8*lq*hl**2*ssz*t1*sz**(-1) + 16*lq*hl**2*ssz
     &    *t1**2*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*m1**2*mt**2*
     &    tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s*t1**(-3)*u1**2*
     &    sz**(-1) + 8*rq*hr**2*ssz*m1**2*s*t1**(-1)*sz**(-1) + 8*rq*
     &    hr**2*ssz*m1**2*s**2*t1**(-3)*u1*sz**(-1) + 8*rq*hr**2*ssz*
     &    m1**2*t1**(-3)*u1**3*sz**(-1) + 32*rq*hr**2*ssz*m1**2*t1*
     &    tx**(-1)*sz**(-1) )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * (  - 8*
     &    rq*hr**2*ssz*m1**2*sz**(-1) + 16*rq*hr**2*ssz*m1**4*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*mt**2*s*t1**(-3)*u1**2*sz**(-1) - 
     &    8*rq*hr**2*ssz*mt**2*s*t1**(-1)*sz**(-1) - 8*rq*hr**2*ssz*
     &    mt**2*s**2*t1**(-3)*u1*sz**(-1) - 8*rq*hr**2*ssz*mt**2*
     &    t1**(-3)*u1**3*sz**(-1) - 32*rq*hr**2*ssz*mt**2*t1*tx**(-1)*
     &    sz**(-1) + 8*rq*hr**2*ssz*mt**2*sz**(-1) + 16*rq*hr**2*ssz*
     &    mt**4*tx**(-1)*sz**(-1) + 24*rq*hr**2*ssz*s*t1**(-3)*u1**3*
     &    sz**(-1) + 16*rq*hr**2*ssz*s*t1**(-2)*u1**2*sz**(-1) - 8*rq*
     &    hr**2*ssz*s*t1**(-1)*u1*sz**(-1) + 8*rq*hr**2*ssz*s*sz**(-1)
     &     + 24*rq*hr**2*ssz*s**2*t1**(-3)*u1**2*sz**(-1) + 8*rq*hr**2*
     &    ssz*s**2*t1**(-2)*u1*sz**(-1) + 8*rq*hr**2*ssz*s**3*t1**(-3)*
     &    u1*sz**(-1) + 8*rq*hr**2*ssz*t1**(-3)*u1**4*sz**(-1) + 8*rq*
     &    hr**2*ssz*t1**(-2)*u1**3*sz**(-1) - 8*rq*hr**2*ssz*t1**(-1)*
     &    u1**2*sz**(-1) - 8*rq*hr**2*ssz*t1*sz**(-1) + 16*rq*hr**2*ssz
     &    *t1**2*tx**(-1)*sz**(-1) )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * (  - 8*
     &    hl**4*m1**2*mt**2*tx**(-2) + 8*hl**4*m1**2*t1*tx**(-2) - 2*
     &    hl**4*m1**2*tx**(-1) + 4*hl**4*m1**4*tx**(-2) - 8*hl**4*mt**2
     &    *t1*tx**(-2) + 2*hl**4*mt**2*tx**(-1) + 4*hl**4*mt**4*
     &    tx**(-2) + 4*hl**4*s*t1**(-3)*u1**2 + 6*hl**4*s*t1**(-3)*
     &    u1**3*tx**(-1) - 2*hl**4*s*t1**(-1)*u1*tx**(-1) + 2*hl**4*s*
     &    t1**(-1) + 2*hl**4*s**2*t1**(-3)*u1 + 6*hl**4*s**2*t1**(-3)*
     &    u1**2*tx**(-1) + 2*hl**4*s**3*t1**(-3)*u1*tx**(-1) + 2*hl**4*
     &    t1**(-3)*u1**3 + 2*hl**4*t1**(-3)*u1**4*tx**(-1) - 2*hl**4*
     &    t1**(-1)*u1**2*tx**(-1) - 2*hl**4*t1*tx**(-1) + 4*hl**4*t1**2
     &    *tx**(-2) - 8*hr**4*m1**2*mt**2*tx**(-2) + 8*hr**4*m1**2*t1*
     &    tx**(-2) - 2*hr**4*m1**2*tx**(-1) + 4*hr**4*m1**4*tx**(-2) - 
     &    8*hr**4*mt**2*t1*tx**(-2) + 2*hr**4*mt**2*tx**(-1) + 4*hr**4*
     &    mt**4*tx**(-2) + 4*hr**4*s*t1**(-3)*u1**2 + 6*hr**4*s*
     &    t1**(-3)*u1**3*tx**(-1) - 2*hr**4*s*t1**(-1)*u1*tx**(-1) + 2*
     &    hr**4*s*t1**(-1) )
      MMv = MMv + log( - m1**(-2)*tx)*Nc*Cf*Pi*alphas*prefac * ( 2*
     &    hr**4*s**2*t1**(-3)*u1 + 6*hr**4*s**2*t1**(-3)*u1**2*tx**(-1)
     &     + 2*hr**4*s**3*t1**(-3)*u1*tx**(-1) + 2*hr**4*t1**(-3)*u1**3
     &     + 2*hr**4*t1**(-3)*u1**4*tx**(-1) - 2*hr**4*t1**(-1)*u1**2*
     &    tx**(-1) - 2*hr**4*t1*tx**(-1) + 4*hr**4*t1**2*tx**(-2) )
      MMv = MMv + Li2(1 - m1**(-2)*s**(-1)*t1*u1)*Nc*Cf*Pi*alphas*
     & prefac * (  - 128*pq*lq*ssz*ssp*pt2*sz**(-1) - 128*pq*rq*ssz*ssp
     &    *pt2*sz**(-1) - 32*pq*hl**2*ssp*pt2*tx**(-1) - 32*pq*hr**2*
     &    ssp*pt2*tx**(-1) - 32*lq*hl**2*ssz*pt2*s*tx**(-1)*sz**(-1) - 
     &    32*rq*hr**2*ssz*pt2*s*tx**(-1)*sz**(-1) - 32*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) - 32*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 8*hl**2*hr**2*mt**2*s*
     &    tx**(-2) - 4*hl**4*pt2*s*tx**(-2) - 4*hr**4*pt2*s*tx**(-2) - 
     &    32*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1) - 16*h1**2*
     &    lambda1**2*s*s1**(-2) - 16*h2**2*lambda2**2*s*s2**(-2) - 64*
     &    ssz**2*lq2*pt2*s*sz**(-2) - 64*ssz**2*rq2*pt2*s*sz**(-2) - 
     &    128*ssp**2*pq2*pt2*s**(-1) )
      MMv = MMv + SCA(1)*Nc*Cf*Pi*alphas*prefac * ( 8*pq*hl**2*ssp*pt2*
     &    m1**2*tx**(-2)*t**(-1) + 8*pq*hl**2*ssp*pt2*mt**2*tx**(-2)*
     &    t**(-1) + 8*pq*hl**2*ssp*pt2*t1*tx**(-2)*t**(-1) - 16*pq*
     &    hl**2*ssp*pt2*tx**(-2) + 8*pq*hr**2*ssp*pt2*m1**2*tx**(-2)*
     &    t**(-1) + 8*pq*hr**2*ssp*pt2*mt**2*tx**(-2)*t**(-1) + 8*pq*
     &    hr**2*ssp*pt2*t1*tx**(-2)*t**(-1) - 16*pq*hr**2*ssp*pt2*
     &    tx**(-2) + 8*lq*hl**2*ssz*pt2*m1**2*s*tx**(-2)*t**(-1)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*mt**2*s*tx**(-2)*t**(-1)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*s*t1*tx**(-2)*t**(-1)*sz**(-1)
     &     - 16*lq*hl**2*ssz*pt2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz*
     &    pt2*m1**2*s*tx**(-2)*t**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*s*t1*
     &    tx**(-2)*t**(-1)*sz**(-1) - 16*rq*hr**2*ssz*pt2*s*tx**(-2)*
     &    sz**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt**(-1)*s*
     &    tx**(-2)*s1**(-1) + 16*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*
     &    s*tx**(-2)*t**(-1)*s1**(-1) )
      MMv = MMv + SCA(1)*Nc*Cf*Pi*alphas*prefac * (  - 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**(-1)*s*t1*tx**(-2)*s1**(-1) + 16*hl*
     &    hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*s1**(-1)
     &     - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*tx**(-2)*s1**(-1) - 8*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**(-1)*s*tx**(-2)*
     &    s2**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s*
     &    tx**(-2)*t**(-1)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**(-1)*s*t1*tx**(-2)*s2**(-1) + 16*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*s2**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*tx**(-2)*s2**(-1) + 8*hl**2*hr**2*
     &    m1**2*mt**2*s*tx**(-3)*t**(-1) - 4*hl**2*hr**2*m1**2*s*
     &    tx**(-3) + 8*hl**2*hr**2*mt**2*s*t1*tx**(-3)*t**(-1) - 4*
     &    hl**2*hr**2*mt**2*s*tx**(-3) - 4*hl**2*hr**2*s*t1*tx**(-3) + 
     &    2*hl**4*pt2*m1**2*s*tx**(-3)*t**(-1) + 2*hl**4*pt2*mt**2*s*
     &    tx**(-3)*t**(-1) + 2*hl**4*pt2*s*t1*tx**(-3)*t**(-1) - 4*
     &    hl**4*pt2*s*tx**(-3) )
      MMv = MMv + SCA(1)*Nc*Cf*Pi*alphas*prefac * ( 2*hr**4*pt2*m1**2*s
     &    *tx**(-3)*t**(-1) + 2*hr**4*pt2*mt**2*s*tx**(-3)*t**(-1) + 2*
     &    hr**4*pt2*s*t1*tx**(-3)*t**(-1) - 4*hr**4*pt2*s*tx**(-3) )
      MMv = MMv + SCB(2,1)*kaellen**(-1)*Nc*Cf*Pi*alphas*prefac * (  - 
     &    32*pq*hl**2*ssp*m1**2*t1 + 8*pq*hl**2*ssp*s**(-1)*t1*u1**2 + 
     &    16*pq*hl**2*ssp*s**(-1)*t1**2*u1 + 8*pq*hl**2*ssp*s**(-1)*
     &    t1**3 - 32*pq*hr**2*ssp*m1**2*t1 + 8*pq*hr**2*ssp*s**(-1)*t1*
     &    u1**2 + 16*pq*hr**2*ssp*s**(-1)*t1**2*u1 + 8*pq*hr**2*ssp*
     &    s**(-1)*t1**3 - 32*lq*hl**2*ssz*m1**2*s*t1*sz**(-1) + 8*lq*
     &    hl**2*ssz*t1*u1**2*sz**(-1) + 16*lq*hl**2*ssz*t1**2*u1*
     &    sz**(-1) + 8*lq*hl**2*ssz*t1**3*sz**(-1) - 32*rq*hr**2*ssz*
     &    m1**2*s*t1*sz**(-1) + 8*rq*hr**2*ssz*t1*u1**2*sz**(-1) + 16*
     &    rq*hr**2*ssz*t1**2*u1*sz**(-1) + 8*rq*hr**2*ssz*t1**3*
     &    sz**(-1) - 8*hl**4*m1**2*s*t1*tx**(-1) + 2*hl**4*t1*u1**2*
     &    tx**(-1) + 4*hl**4*t1**2*u1*tx**(-1) + 2*hl**4*t1**3*tx**(-1)
     &     - 8*hr**4*m1**2*s*t1*tx**(-1) + 2*hr**4*t1*u1**2*tx**(-1) + 
     &    4*hr**4*t1**2*u1*tx**(-1) + 2*hr**4*t1**3*tx**(-1) )
      MMv = MMv + SCB(2,1)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hl**2*ssp*
     &    m1**2*t1**(-1) - 16*pq*hl**2*ssp*m1**2*tx**(-1) - 16*pq*hl**2
     &    *ssp*m1**4*t1**(-1)*tx**(-1) + 16*pq*hl**2*ssp*mt**2*s**(-1)*
     &    u1*tx**(-1) - 8*pq*hl**2*ssp*s**(-1)*t1 + 16*pq*hr**2*ssp*
     &    m1**2*t1**(-1) - 16*pq*hr**2*ssp*m1**2*tx**(-1) - 16*pq*hr**2
     &    *ssp*m1**4*t1**(-1)*tx**(-1) + 16*pq*hr**2*ssp*mt**2*s**(-1)*
     &    u1*tx**(-1) - 8*pq*hr**2*ssp*s**(-1)*t1 + 16*lq*hl**2*ssz*
     &    m1**2*s*t1**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**2*s*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*m1**4*s*t1**(-1)*tx**(-1)*sz**(-1)
     &     + 16*lq*hl**2*ssz*mt**2*u1*tx**(-1)*sz**(-1) - 8*lq*hl**2*
     &    ssz*t1*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s*t1**(-1)*sz**(-1)
     &     - 16*rq*hr**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 16*rq*hr**2*
     &    ssz*m1**4*s*t1**(-1)*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    mt**2*u1*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*t1*sz**(-1) + 16*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) )
      MMv = MMv + SCB(2,1)*Nc*Cf*Pi*alphas*prefac * ( 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) + 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*t1**(-1)*tx**(-1)*s2**(-1) + 
     &    16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 8*
     &    hl**2*hr**2*m1**2*s*t1**(-1)*tx**(-1) + 8*hl**2*hr**2*m1**2*s
     &    *tx**(-2) + 8*hl**2*hr**2*m1**4*s*t1**(-1)*tx**(-2) + 8*hl**2
     &    *hr**2*mt**2*s*tx**(-2) + 4*hl**4*m1**2*s*t1**(-1)*tx**(-1)
     &     - 4*hl**4*m1**2*s*tx**(-2) - 4*hl**4*m1**4*s*t1**(-1)*
     &    tx**(-2) + 4*hl**4*mt**2*u1*tx**(-2) - 2*hl**4*t1*tx**(-1) + 
     &    4*hr**4*m1**2*s*t1**(-1)*tx**(-1) - 4*hr**4*m1**2*s*tx**(-2)
     &     - 4*hr**4*m1**4*s*t1**(-1)*tx**(-2) + 4*hr**4*mt**2*u1*
     &    tx**(-2) - 2*hr**4*t1*tx**(-1) )
      MMv = MMv + SCB(3,1)*kaellen**(-1)*Nc*Cf*Pi*alphas*prefac * ( 32*
     &    pq*hl**2*ssp*m1**2*t1 - 8*pq*hl**2*ssp*s**(-1)*t1*u1**2 - 16*
     &    pq*hl**2*ssp*s**(-1)*t1**2*u1 - 8*pq*hl**2*ssp*s**(-1)*t1**3
     &     + 32*pq*hr**2*ssp*m1**2*t1 - 8*pq*hr**2*ssp*s**(-1)*t1*u1**2
     &     - 16*pq*hr**2*ssp*s**(-1)*t1**2*u1 - 8*pq*hr**2*ssp*s**(-1)*
     &    t1**3 + 32*lq*hl**2*ssz*m1**2*s*t1*sz**(-1) - 8*lq*hl**2*ssz*
     &    t1*u1**2*sz**(-1) - 16*lq*hl**2*ssz*t1**2*u1*sz**(-1) - 8*lq*
     &    hl**2*ssz*t1**3*sz**(-1) + 32*rq*hr**2*ssz*m1**2*s*t1*
     &    sz**(-1) - 8*rq*hr**2*ssz*t1*u1**2*sz**(-1) - 16*rq*hr**2*ssz
     &    *t1**2*u1*sz**(-1) - 8*rq*hr**2*ssz*t1**3*sz**(-1) + 8*hl**4*
     &    m1**2*s*t1*tx**(-1) - 2*hl**4*t1*u1**2*tx**(-1) - 4*hl**4*
     &    t1**2*u1*tx**(-1) - 2*hl**4*t1**3*tx**(-1) + 8*hr**4*m1**2*s*
     &    t1*tx**(-1) - 2*hr**4*t1*u1**2*tx**(-1) - 4*hr**4*t1**2*u1*
     &    tx**(-1) - 2*hr**4*t1**3*tx**(-1) )
      MMv = MMv + SCB(3,1)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hl**2*ssp*
     &    m1**2*t1**(-1) - 16*pq*hl**2*ssp*m1**2*tx**(-1) - 16*pq*hl**2
     &    *ssp*m1**4*t1**(-1)*tx**(-1) + 16*pq*hl**2*ssp*mt**2*s**(-1)*
     &    u1*tx**(-1) + 8*pq*hl**2*ssp*s**(-1)*t1**(-2)*u1**3 + 8*pq*
     &    hl**2*ssp*s**(-1)*t1 - 8*pq*hl**2*ssp*s**(-1)*u1 + 8*pq*hl**2
     &    *ssp*s*t1**(-2)*u1 + 16*pq*hl**2*ssp*t1**(-2)*u1**2 + 16*pq*
     &    hr**2*ssp*m1**2*t1**(-1) - 16*pq*hr**2*ssp*m1**2*tx**(-1) - 
     &    16*pq*hr**2*ssp*m1**4*t1**(-1)*tx**(-1) + 16*pq*hr**2*ssp*
     &    mt**2*s**(-1)*u1*tx**(-1) + 8*pq*hr**2*ssp*s**(-1)*t1**(-2)*
     &    u1**3 + 8*pq*hr**2*ssp*s**(-1)*t1 - 8*pq*hr**2*ssp*s**(-1)*u1
     &     + 8*pq*hr**2*ssp*s*t1**(-2)*u1 + 16*pq*hr**2*ssp*t1**(-2)*
     &    u1**2 + 16*lq*hl**2*ssz*m1**2*s*t1**(-1)*sz**(-1) - 16*lq*
     &    hl**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**4*s
     &    *t1**(-1)*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*mt**2*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*s*t1**(-2)*u1**2*sz**(-1)
     &     + 8*lq*hl**2*ssz*s**2*t1**(-2)*u1*sz**(-1) )
      MMv = MMv + SCB(3,1)*Nc*Cf*Pi*alphas*prefac * ( 8*lq*hl**2*ssz*
     &    t1**(-2)*u1**3*sz**(-1) + 8*lq*hl**2*ssz*t1*sz**(-1) - 8*lq*
     &    hl**2*ssz*u1*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s*t1**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*m1**4*s*t1**(-1)*tx**(-1)*sz**(-1) + 16*rq*hr**2*
     &    ssz*mt**2*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*s*t1**(-2)*
     &    u1**2*sz**(-1) + 8*rq*hr**2*ssz*s**2*t1**(-2)*u1*sz**(-1) + 8
     &    *rq*hr**2*ssz*t1**(-2)*u1**3*sz**(-1) + 8*rq*hr**2*ssz*t1*
     &    sz**(-1) - 8*rq*hr**2*ssz*u1*sz**(-1) + 16*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) + 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*t1**(-1)*tx**(-1)*s2**(-1) + 
     &    16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 8*
     &    hl**2*hr**2*m1**2*s*t1**(-1)*tx**(-1) + 8*hl**2*hr**2*m1**2*s
     &    *tx**(-2) + 8*hl**2*hr**2*m1**4*s*t1**(-1)*tx**(-2) + 8*hl**2
     &    *hr**2*mt**2*s*tx**(-2) )
      MMv = MMv + SCB(3,1)*Nc*Cf*Pi*alphas*prefac * ( 4*hl**4*m1**2*s*
     &    t1**(-1)*tx**(-1) - 4*hl**4*m1**2*s*tx**(-2) - 4*hl**4*m1**4*
     &    s*t1**(-1)*tx**(-2) + 4*hl**4*mt**2*u1*tx**(-2) + 4*hl**4*s*
     &    t1**(-2)*u1**2*tx**(-1) + 2*hl**4*s**2*t1**(-2)*u1*tx**(-1)
     &     + 2*hl**4*t1**(-2)*u1**3*tx**(-1) + 2*hl**4*t1*tx**(-1) - 2*
     &    hl**4*u1*tx**(-1) + 4*hr**4*m1**2*s*t1**(-1)*tx**(-1) - 4*
     &    hr**4*m1**2*s*tx**(-2) - 4*hr**4*m1**4*s*t1**(-1)*tx**(-2) + 
     &    4*hr**4*mt**2*u1*tx**(-2) + 4*hr**4*s*t1**(-2)*u1**2*tx**(-1)
     &     + 2*hr**4*s**2*t1**(-2)*u1*tx**(-1) + 2*hr**4*t1**(-2)*u1**3
     &    *tx**(-1) + 2*hr**4*t1*tx**(-1) - 2*hr**4*u1*tx**(-1) )
      MMv = MMv + SCB(4,1)*Nc*Cf*Pi*alphas*prefac * ( 192*pq*lq*ssz*ssp
     &    *m1**2*sz**(-1) - 192*pq*lq*ssz*ssp*s**(-1)*t1*u1*sz**(-1) + 
     &    192*pq*rq*ssz*ssp*m1**2*sz**(-1) - 192*pq*rq*ssz*ssp*s**(-1)*
     &    t1*u1*sz**(-1) + 24*pq*hl**2*ssp*m1**2*tx**(-1) - 24*pq*hl**2
     &    *ssp*s**(-1)*t1*u1*tx**(-1) + 24*pq*hr**2*ssp*m1**2*tx**(-1)
     &     - 24*pq*hr**2*ssp*s**(-1)*t1*u1*tx**(-1) + 24*lq*hl**2*ssz*
     &    m1**2*s*tx**(-1)*sz**(-1) - 24*lq*hl**2*ssz*t1*u1*tx**(-1)*
     &    sz**(-1) + 24*rq*hr**2*ssz*m1**2*s*tx**(-1)*sz**(-1) - 24*rq*
     &    hr**2*ssz*t1*u1*tx**(-1)*sz**(-1) + 96*ssz**2*lq2*m1**2*s*
     &    sz**(-2) - 96*ssz**2*lq2*t1*u1*sz**(-2) + 96*ssz**2*rq2*m1**2
     &    *s*sz**(-2) - 96*ssz**2*rq2*t1*u1*sz**(-2) + 192*ssp**2*pq2*
     &    m1**2*s**(-1) - 192*ssp**2*pq2*s**(-2)*t1*u1 )
      MMv = MMv + SCB(6,1)*Nc*Cf*Pi*alphas*prefac * (  - 8*pq*hl**2*ssp
     &    *pt2*m1**2*mt**2*tx**(-2)*t**(-1) - 8*pq*hl**2*ssp*pt2*m1**2*
     &    tx**(-2) - 8*pq*hl**2*ssp*pt2*mt**2*t1*tx**(-2)*t**(-1) + 56*
     &    pq*hl**2*ssp*pt2*mt**2*tx**(-2) - 8*pq*hl**2*ssp*pt2*mt**4*
     &    tx**(-2)*t**(-1) - 8*pq*hl**2*ssp*pt2*t1*tx**(-2) - 32*pq*
     &    hl**2*ssp*m1**2*t1**(-1) + 32*pq*hl**2*ssp*m1**4*t1**(-1)*
     &    tx**(-1) - 32*pq*hl**2*ssp*mt**2*s**(-1)*u1*tx**(-1) + 32*pq*
     &    hl**2*ssp*s**(-1)*t1*u1*tx**(-1) - 8*pq*hr**2*ssp*pt2*m1**2*
     &    mt**2*tx**(-2)*t**(-1) - 8*pq*hr**2*ssp*pt2*m1**2*tx**(-2) - 
     &    8*pq*hr**2*ssp*pt2*mt**2*t1*tx**(-2)*t**(-1) + 56*pq*hr**2*
     &    ssp*pt2*mt**2*tx**(-2) - 8*pq*hr**2*ssp*pt2*mt**4*tx**(-2)*
     &    t**(-1) - 8*pq*hr**2*ssp*pt2*t1*tx**(-2) - 32*pq*hr**2*ssp*
     &    m1**2*t1**(-1) + 32*pq*hr**2*ssp*m1**4*t1**(-1)*tx**(-1) - 32
     &    *pq*hr**2*ssp*mt**2*s**(-1)*u1*tx**(-1) + 32*pq*hr**2*ssp*
     &    s**(-1)*t1*u1*tx**(-1) - 8*lq*hl**2*ssz*pt2*m1**2*mt**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) )
      MMv = MMv + SCB(6,1)*Nc*Cf*Pi*alphas*prefac * (  - 8*lq*hl**2*ssz
     &    *pt2*m1**2*s*tx**(-2)*sz**(-1) - 8*lq*hl**2*ssz*pt2*mt**2*s*
     &    t1*tx**(-2)*t**(-1)*sz**(-1) + 56*lq*hl**2*ssz*pt2*mt**2*s*
     &    tx**(-2)*sz**(-1) - 8*lq*hl**2*ssz*pt2*mt**4*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 8*lq*hl**2*ssz*pt2*s*t1*tx**(-2)*sz**(-1)
     &     - 32*lq*hl**2*ssz*m1**2*s*t1**(-1)*sz**(-1) + 32*lq*hl**2*
     &    ssz*m1**4*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*
     &    mt**2*u1*tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*t1*u1*tx**(-1)*
     &    sz**(-1) - 8*rq*hr**2*ssz*pt2*m1**2*mt**2*s*tx**(-2)*t**(-1)*
     &    sz**(-1) - 8*rq*hr**2*ssz*pt2*m1**2*s*tx**(-2)*sz**(-1) - 8*
     &    rq*hr**2*ssz*pt2*mt**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 56*rq
     &    *hr**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) - 8*rq*hr**2*ssz*pt2
     &    *mt**4*s*tx**(-2)*t**(-1)*sz**(-1) - 8*rq*hr**2*ssz*pt2*s*t1*
     &    tx**(-2)*sz**(-1) - 32*rq*hr**2*ssz*m1**2*s*t1**(-1)*sz**(-1)
     &     + 32*rq*hr**2*ssz*m1**4*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*rq
     &    *hr**2*ssz*mt**2*u1*tx**(-1)*sz**(-1) )
      MMv = MMv + SCB(6,1)*Nc*Cf*Pi*alphas*prefac * ( 32*rq*hr**2*ssz*
     &    t1*u1*tx**(-1)*sz**(-1) - 32*hl*hr*h1*lambda1*sqrt2**(-1)*
     &    m1**2*mt*s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) - 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt**3*s*tx**(-2)*t**(-1)*s1**(-1)
     &     + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*s1**(-1)
     &     - 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*s*t1*tx**(-2)*
     &    t**(-1)*s1**(-1) + 32*hl*hr*h1*lambda1*sqrt2**(-1)*mt**3*s*
     &    tx**(-2)*s1**(-1) - 32*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    s*t1**(-1)*tx**(-1)*s2**(-1) + 16*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s2**(-1) - 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt**3*s*tx**(-2)*t**(-1)*s2**(-1)
     &     + 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1)
     &     - 16*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*s*t1*tx**(-2)*
     &    t**(-1)*s2**(-1) + 32*hl*hr*h2*lambda2*sqrt2**(-1)*mt**3*s*
     &    tx**(-2)*s2**(-1) )
      MMv = MMv + SCB(6,1)*Nc*Cf*Pi*alphas*prefac * ( 8*hl**2*hr**2*
     &    m1**2*mt**2*s*tx**(-3) - 8*hl**2*hr**2*m1**2*mt**4*s*tx**(-3)
     &    *t**(-1) + 16*hl**2*hr**2*m1**2*s*t1**(-1)*tx**(-1) - 16*
     &    hl**2*hr**2*m1**2*s*tx**(-2) - 16*hl**2*hr**2*m1**4*s*
     &    t1**(-1)*tx**(-2) + 8*hl**2*hr**2*mt**2*s*t1*tx**(-3) - 8*
     &    hl**2*hr**2*mt**4*s*t1*tx**(-3)*t**(-1) + 16*hl**2*hr**2*
     &    mt**4*s*tx**(-3) - 2*hl**4*pt2*m1**2*mt**2*s*tx**(-3)*t**(-1)
     &     - 2*hl**4*pt2*m1**2*s*tx**(-3) - 2*hl**4*pt2*mt**2*s*t1*
     &    tx**(-3)*t**(-1) + 14*hl**4*pt2*mt**2*s*tx**(-3) - 2*hl**4*
     &    pt2*mt**4*s*tx**(-3)*t**(-1) - 2*hl**4*pt2*s*t1*tx**(-3) - 8*
     &    hl**4*m1**2*s*t1**(-1)*tx**(-1) + 8*hl**4*m1**4*s*t1**(-1)*
     &    tx**(-2) - 8*hl**4*mt**2*u1*tx**(-2) + 8*hl**4*t1*u1*tx**(-2)
     &     - 2*hr**4*pt2*m1**2*mt**2*s*tx**(-3)*t**(-1) - 2*hr**4*pt2*
     &    m1**2*s*tx**(-3) - 2*hr**4*pt2*mt**2*s*t1*tx**(-3)*t**(-1) + 
     &    14*hr**4*pt2*mt**2*s*tx**(-3) - 2*hr**4*pt2*mt**4*s*tx**(-3)*
     &    t**(-1) )
      MMv = MMv + SCB(6,1)*Nc*Cf*Pi*alphas*prefac * (  - 2*hr**4*pt2*s*
     &    t1*tx**(-3) - 8*hr**4*m1**2*s*t1**(-1)*tx**(-1) + 8*hr**4*
     &    m1**4*s*t1**(-1)*tx**(-2) - 8*hr**4*mt**2*u1*tx**(-2) + 8*
     &    hr**4*t1*u1*tx**(-2) )
      MMv = MMv + SCB(7,1)*Nc*Cf*Pi*alphas*prefac * (  - 32*pq*hl**2*
     &    ssp*pt2*mt**2*tx**(-2) - 32*pq*hr**2*ssp*pt2*mt**2*tx**(-2)
     &     - 32*lq*hl**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) - 32*rq*
     &    hr**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) - 16*hl*hr*h1*lambda1
     &    *sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) - 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*s1**(-1) - 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt**3*s*tx**(-2)*s1**(-1) - 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s2**(-1) - 16*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) - 16*hl*hr*
     &    h2*lambda2*sqrt2**(-1)*mt**3*s*tx**(-2)*s2**(-1) - 8*hl**2*
     &    hr**2*m1**2*mt**2*s*tx**(-3) - 8*hl**2*hr**2*mt**2*s*t1*
     &    tx**(-3) - 8*hl**2*hr**2*mt**4*s*tx**(-3) - 8*hl**4*pt2*mt**2
     &    *s*tx**(-3) - 8*hr**4*pt2*mt**2*s*tx**(-3) )
      MMv = MMv + SCC(1,1)*Nc*Cf*Pi*alphas*prefac * ( 128*pq*lq*ssz*ssp
     &    *m1**2*s*sz**(-1) - 128*pq*lq*ssz*ssp*t1*u1*sz**(-1) + 128*pq
     &    *rq*ssz*ssp*m1**2*s*sz**(-1) - 128*pq*rq*ssz*ssp*t1*u1*
     &    sz**(-1) + 16*pq*hl**2*ssp*m1**2*s*tx**(-1) - 8*pq*hl**2*ssp*
     &    m1**2 + 8*pq*hl**2*ssp*mt**2 - 16*pq*hl**2*ssp*t1*u1*tx**(-1)
     &     - 8*pq*hl**2*ssp*u1 + 16*pq*hr**2*ssp*m1**2*s*tx**(-1) - 8*
     &    pq*hr**2*ssp*m1**2 + 8*pq*hr**2*ssp*mt**2 - 16*pq*hr**2*ssp*
     &    t1*u1*tx**(-1) - 8*pq*hr**2*ssp*u1 - 8*lq*hl**2*ssz*m1**2*s*
     &    sz**(-1) + 16*lq*hl**2*ssz*m1**2*s**2*tx**(-1)*sz**(-1) + 8*
     &    lq*hl**2*ssz*mt**2*s*sz**(-1) - 16*lq*hl**2*ssz*s*t1*u1*
     &    tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*s*u1*sz**(-1) - 8*rq*hr**2
     &    *ssz*m1**2*s*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s**2*tx**(-1)*
     &    sz**(-1) + 8*rq*hr**2*ssz*mt**2*s*sz**(-1) - 16*rq*hr**2*ssz*
     &    s*t1*u1*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*s*u1*sz**(-1) - 16
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*tx**(-1)*s1**(-1) - 16*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s**2*tx**(-1)*s2**(-1) )
      MMv = MMv + SCC(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 2*hl**4*m1**2*
     &    s*tx**(-1) + 2*hl**4*mt**2*s*tx**(-1) - 2*hl**4*s*u1*tx**(-1)
     &     - 2*hr**4*m1**2*s*tx**(-1) + 2*hr**4*mt**2*s*tx**(-1) - 2*
     &    hr**4*s*u1*tx**(-1) - 32*h1*h2*lambda1*lambda2*s**2*s1**(-1)*
     &    s2**(-1) - 16*h1**2*lambda1**2*s**2*s1**(-2) - 16*h2**2*
     &    lambda2**2*s**2*s2**(-2) + 64*ssz**2*lq2*m1**2*s**2*sz**(-2)
     &     - 64*ssz**2*lq2*s*t1*u1*sz**(-2) + 64*ssz**2*rq2*m1**2*s**2*
     &    sz**(-2) - 64*ssz**2*rq2*s*t1*u1*sz**(-2) + 128*ssp**2*pq2*
     &    m1**2 - 128*ssp**2*pq2*s**(-1)*t1*u1 )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * (  - 32*pq*hl**2*
     &    ssp*m1**2*mt**2*tx**(-1) + 8*pq*hl**2*ssp*m1**2*s**(-1)*
     &    t1**(-2)*u1**3 - 32*pq*hl**2*ssp*m1**2*s**(-1)*t1*u1*tx**(-1)
     &     - 8*pq*hl**2*ssp*m1**2*s**(-1)*t1 + 8*pq*hl**2*ssp*m1**2*s*
     &    t1**(-2)*u1 + 16*pq*hl**2*ssp*m1**2*t1**(-2)*u1**2 - 24*pq*
     &    hl**2*ssp*m1**2 + 32*pq*hl**2*ssp*m1**4*tx**(-1) - 8*pq*hl**2
     &    *ssp*mt**2*s**(-1)*t1**(-2)*u1**3 + 32*pq*hl**2*ssp*mt**2*
     &    s**(-1)*t1*u1*tx**(-1) + 8*pq*hl**2*ssp*mt**2*s**(-1)*t1 - 8*
     &    pq*hl**2*ssp*mt**2*s*t1**(-2)*u1 - 16*pq*hl**2*ssp*mt**2*
     &    t1**(-2)*u1**2 - 8*pq*hl**2*ssp*mt**2 + 8*pq*hl**2*ssp*
     &    s**(-1)*t1**(-2)*u1**4 + 8*pq*hl**2*ssp*s**(-1)*t1**(-1)*
     &    u1**3 + 16*pq*hl**2*ssp*s**(-1)*t1*u1 + 24*pq*hl**2*ssp*s*
     &    t1**(-2)*u1**2 + 8*pq*hl**2*ssp*s*t1**(-1)*u1 + 8*pq*hl**2*
     &    ssp*s**2*t1**(-2)*u1 + 24*pq*hl**2*ssp*t1**(-2)*u1**3 + 16*pq
     &    *hl**2*ssp*t1**(-1)*u1**2 - 32*pq*hr**2*ssp*m1**2*mt**2*
     &    tx**(-1) )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * ( 8*pq*hr**2*ssp*
     &    m1**2*s**(-1)*t1**(-2)*u1**3 - 32*pq*hr**2*ssp*m1**2*s**(-1)*
     &    t1*u1*tx**(-1) - 8*pq*hr**2*ssp*m1**2*s**(-1)*t1 + 8*pq*hr**2
     &    *ssp*m1**2*s*t1**(-2)*u1 + 16*pq*hr**2*ssp*m1**2*t1**(-2)*
     &    u1**2 - 24*pq*hr**2*ssp*m1**2 + 32*pq*hr**2*ssp*m1**4*
     &    tx**(-1) - 8*pq*hr**2*ssp*mt**2*s**(-1)*t1**(-2)*u1**3 + 32*
     &    pq*hr**2*ssp*mt**2*s**(-1)*t1*u1*tx**(-1) + 8*pq*hr**2*ssp*
     &    mt**2*s**(-1)*t1 - 8*pq*hr**2*ssp*mt**2*s*t1**(-2)*u1 - 16*pq
     &    *hr**2*ssp*mt**2*t1**(-2)*u1**2 - 8*pq*hr**2*ssp*mt**2 + 8*pq
     &    *hr**2*ssp*s**(-1)*t1**(-2)*u1**4 + 8*pq*hr**2*ssp*s**(-1)*
     &    t1**(-1)*u1**3 + 16*pq*hr**2*ssp*s**(-1)*t1*u1 + 24*pq*hr**2*
     &    ssp*s*t1**(-2)*u1**2 + 8*pq*hr**2*ssp*s*t1**(-1)*u1 + 8*pq*
     &    hr**2*ssp*s**2*t1**(-2)*u1 + 24*pq*hr**2*ssp*t1**(-2)*u1**3
     &     + 16*pq*hr**2*ssp*t1**(-1)*u1**2 - 32*lq*hl**2*ssz*m1**2*
     &    mt**2*s*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*m1**2*s*t1**(-2)*
     &    u1**2*sz**(-1) )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * (  - 24*lq*hl**2*
     &    ssz*m1**2*s*sz**(-1) + 8*lq*hl**2*ssz*m1**2*s**2*t1**(-2)*u1*
     &    sz**(-1) + 8*lq*hl**2*ssz*m1**2*t1**(-2)*u1**3*sz**(-1) - 32*
     &    lq*hl**2*ssz*m1**2*t1*u1*tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*
     &    m1**2*t1*sz**(-1) + 32*lq*hl**2*ssz*m1**4*s*tx**(-1)*sz**(-1)
     &     - 16*lq*hl**2*ssz*mt**2*s*t1**(-2)*u1**2*sz**(-1) - 8*lq*
     &    hl**2*ssz*mt**2*s*sz**(-1) - 8*lq*hl**2*ssz*mt**2*s**2*
     &    t1**(-2)*u1*sz**(-1) - 8*lq*hl**2*ssz*mt**2*t1**(-2)*u1**3*
     &    sz**(-1) + 32*lq*hl**2*ssz*mt**2*t1*u1*tx**(-1)*sz**(-1) + 8*
     &    lq*hl**2*ssz*mt**2*t1*sz**(-1) + 24*lq*hl**2*ssz*s*t1**(-2)*
     &    u1**3*sz**(-1) + 16*lq*hl**2*ssz*s*t1**(-1)*u1**2*sz**(-1) + 
     &    24*lq*hl**2*ssz*s**2*t1**(-2)*u1**2*sz**(-1) + 8*lq*hl**2*ssz
     &    *s**2*t1**(-1)*u1*sz**(-1) + 8*lq*hl**2*ssz*s**3*t1**(-2)*u1*
     &    sz**(-1) + 8*lq*hl**2*ssz*t1**(-2)*u1**4*sz**(-1) + 8*lq*
     &    hl**2*ssz*t1**(-1)*u1**3*sz**(-1) + 16*lq*hl**2*ssz*t1*u1*
     &    sz**(-1) )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*hr**2*
     &    ssz*m1**2*mt**2*s*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s
     &    *t1**(-2)*u1**2*sz**(-1) - 24*rq*hr**2*ssz*m1**2*s*sz**(-1)
     &     + 8*rq*hr**2*ssz*m1**2*s**2*t1**(-2)*u1*sz**(-1) + 8*rq*
     &    hr**2*ssz*m1**2*t1**(-2)*u1**3*sz**(-1) - 32*rq*hr**2*ssz*
     &    m1**2*t1*u1*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*m1**2*t1*
     &    sz**(-1) + 32*rq*hr**2*ssz*m1**4*s*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*mt**2*s*t1**(-2)*u1**2*sz**(-1) - 8*rq*hr**2*ssz*
     &    mt**2*s*sz**(-1) - 8*rq*hr**2*ssz*mt**2*s**2*t1**(-2)*u1*
     &    sz**(-1) - 8*rq*hr**2*ssz*mt**2*t1**(-2)*u1**3*sz**(-1) + 32*
     &    rq*hr**2*ssz*mt**2*t1*u1*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*
     &    mt**2*t1*sz**(-1) + 24*rq*hr**2*ssz*s*t1**(-2)*u1**3*sz**(-1)
     &     + 16*rq*hr**2*ssz*s*t1**(-1)*u1**2*sz**(-1) + 24*rq*hr**2*
     &    ssz*s**2*t1**(-2)*u1**2*sz**(-1) + 8*rq*hr**2*ssz*s**2*
     &    t1**(-1)*u1*sz**(-1) + 8*rq*hr**2*ssz*s**3*t1**(-2)*u1*
     &    sz**(-1) )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * ( 8*rq*hr**2*ssz*
     &    t1**(-2)*u1**4*sz**(-1) + 8*rq*hr**2*ssz*t1**(-1)*u1**3*
     &    sz**(-1) + 16*rq*hr**2*ssz*t1*u1*sz**(-1) - 32*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-1)*s1**(-1) + 32*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*s1**(-1) + 32*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*mt**3*s*tx**(-1)*s1**(-1) - 32*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mt*s*tx**(-1)*s2**(-1) + 32*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*s2**(-1) + 32*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt**3*s*tx**(-1)*s2**(-1) - 16*hl**2*hr**2*m1**2*
     &    mt**2*s*tx**(-2) + 16*hl**2*hr**2*mt**2*s*tx**(-1) + 16*hl**2
     &    *hr**2*mt**4*s*tx**(-2) - 8*hl**4*m1**2*mt**2*s*tx**(-2) - 6*
     &    hl**4*m1**2*s*tx**(-1) - 8*hl**4*m1**2*t1*u1*tx**(-2) - 2*
     &    hl**4*m1**2*t1*tx**(-1) + 8*hl**4*m1**4*s*tx**(-2) - 2*hl**4*
     &    mt**2*s*tx**(-1) + 8*hl**4*mt**2*t1*u1*tx**(-2) + 2*hl**4*
     &    mt**2*t1*tx**(-1) + 4*hl**4*s*t1**(-2)*u1**2 + 6*hl**4*s*
     &    t1**(-2)*u1**3*tx**(-1) )
      MMv = MMv + SCC(2,1)*Nc*Cf*Pi*alphas*prefac * ( 2*hl**4*s**2*
     &    t1**(-2)*u1 + 6*hl**4*s**2*t1**(-2)*u1**2*tx**(-1) + 2*hl**4*
     &    s**3*t1**(-2)*u1*tx**(-1) + 2*hl**4*t1**(-2)*u1**3 + 2*hl**4*
     &    t1**(-2)*u1**4*tx**(-1) + 4*hl**4*t1*u1*tx**(-1) - 8*hr**4*
     &    m1**2*mt**2*s*tx**(-2) - 6*hr**4*m1**2*s*tx**(-1) - 8*hr**4*
     &    m1**2*t1*u1*tx**(-2) - 2*hr**4*m1**2*t1*tx**(-1) + 8*hr**4*
     &    m1**4*s*tx**(-2) - 2*hr**4*mt**2*s*tx**(-1) + 8*hr**4*mt**2*
     &    t1*u1*tx**(-2) + 2*hr**4*mt**2*t1*tx**(-1) + 4*hr**4*s*
     &    t1**(-2)*u1**2 + 6*hr**4*s*t1**(-2)*u1**3*tx**(-1) + 2*hr**4*
     &    s**2*t1**(-2)*u1 + 6*hr**4*s**2*t1**(-2)*u1**2*tx**(-1) + 2*
     &    hr**4*s**3*t1**(-2)*u1*tx**(-1) + 2*hr**4*t1**(-2)*u1**3 + 2*
     &    hr**4*t1**(-2)*u1**4*tx**(-1) + 4*hr**4*t1*u1*tx**(-1) )
      MMv = MMv + SCC(6,1)*kaellen**(-1)*Nc*Cf*Pi*alphas*prefac * (  - 
     &    32*pq*hl**2*ssp*m1**2*s*t1 + 8*pq*hl**2*ssp*t1*u1**2 + 16*pq*
     &    hl**2*ssp*t1**2*u1 + 8*pq*hl**2*ssp*t1**3 - 32*pq*hr**2*ssp*
     &    m1**2*s*t1 + 8*pq*hr**2*ssp*t1*u1**2 + 16*pq*hr**2*ssp*t1**2*
     &    u1 + 8*pq*hr**2*ssp*t1**3 - 32*lq*hl**2*ssz*m1**2*s**2*t1*
     &    sz**(-1) + 8*lq*hl**2*ssz*s*t1*u1**2*sz**(-1) + 16*lq*hl**2*
     &    ssz*s*t1**2*u1*sz**(-1) + 8*lq*hl**2*ssz*s*t1**3*sz**(-1) - 
     &    32*rq*hr**2*ssz*m1**2*s**2*t1*sz**(-1) + 8*rq*hr**2*ssz*s*t1*
     &    u1**2*sz**(-1) + 16*rq*hr**2*ssz*s*t1**2*u1*sz**(-1) + 8*rq*
     &    hr**2*ssz*s*t1**3*sz**(-1) - 8*hl**4*m1**2*s**2*t1*tx**(-1)
     &     + 2*hl**4*s*t1*u1**2*tx**(-1) + 4*hl**4*s*t1**2*u1*tx**(-1)
     &     + 2*hl**4*s*t1**3*tx**(-1) - 8*hr**4*m1**2*s**2*t1*tx**(-1)
     &     + 2*hr**4*s*t1*u1**2*tx**(-1) + 4*hr**4*s*t1**2*u1*tx**(-1)
     &     + 2*hr**4*s*t1**3*tx**(-1) )
      MMv = MMv + SCC(6,1)*Nc*Cf*Pi*alphas*prefac * ( 8*pq*hl**2*ssp*
     &    m1**2*s**(-1)*t1 - 8*pq*hl**2*ssp*m1**2*s**(-1)*u1 + 16*pq*
     &    hl**2*ssp*m1**2 - 8*pq*hl**2*ssp*mt**2*s**(-1)*t1 + 8*pq*
     &    hl**2*ssp*mt**2*s**(-1)*u1 - 8*pq*hl**2*ssp*s**(-1)*t1*u1 - 8
     &    *pq*hl**2*ssp*s**(-1)*u1**2 - 8*pq*hl**2*ssp*t1 + 8*pq*hr**2*
     &    ssp*m1**2*s**(-1)*t1 - 8*pq*hr**2*ssp*m1**2*s**(-1)*u1 + 16*
     &    pq*hr**2*ssp*m1**2 - 8*pq*hr**2*ssp*mt**2*s**(-1)*t1 + 8*pq*
     &    hr**2*ssp*mt**2*s**(-1)*u1 - 8*pq*hr**2*ssp*s**(-1)*t1*u1 - 8
     &    *pq*hr**2*ssp*s**(-1)*u1**2 - 8*pq*hr**2*ssp*t1 + 16*lq*hl**2
     &    *ssz*m1**2*s*sz**(-1) + 8*lq*hl**2*ssz*m1**2*t1*sz**(-1) - 8*
     &    lq*hl**2*ssz*m1**2*u1*sz**(-1) - 8*lq*hl**2*ssz*mt**2*t1*
     &    sz**(-1) + 8*lq*hl**2*ssz*mt**2*u1*sz**(-1) - 8*lq*hl**2*ssz*
     &    s*t1*sz**(-1) - 8*lq*hl**2*ssz*t1*u1*sz**(-1) - 8*lq*hl**2*
     &    ssz*u1**2*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s*sz**(-1) + 8*rq*
     &    hr**2*ssz*m1**2*t1*sz**(-1) - 8*rq*hr**2*ssz*m1**2*u1*
     &    sz**(-1) )
      MMv = MMv + SCC(6,1)*Nc*Cf*Pi*alphas*prefac * (  - 8*rq*hr**2*ssz
     &    *mt**2*t1*sz**(-1) + 8*rq*hr**2*ssz*mt**2*u1*sz**(-1) - 8*rq*
     &    hr**2*ssz*s*t1*sz**(-1) - 8*rq*hr**2*ssz*t1*u1*sz**(-1) - 8*
     &    rq*hr**2*ssz*u1**2*sz**(-1) + 4*hl**4*m1**2*s*tx**(-1) + 2*
     &    hl**4*m1**2*t1*tx**(-1) - 2*hl**4*m1**2*u1*tx**(-1) - 2*hl**4
     &    *mt**2*t1*tx**(-1) + 2*hl**4*mt**2*u1*tx**(-1) - 2*hl**4*s*t1
     &    *tx**(-1) - 2*hl**4*t1*u1*tx**(-1) - 2*hl**4*u1**2*tx**(-1)
     &     + 4*hr**4*m1**2*s*tx**(-1) + 2*hr**4*m1**2*t1*tx**(-1) - 2*
     &    hr**4*m1**2*u1*tx**(-1) - 2*hr**4*mt**2*t1*tx**(-1) + 2*hr**4
     &    *mt**2*u1*tx**(-1) - 2*hr**4*s*t1*tx**(-1) - 2*hr**4*t1*u1*
     &    tx**(-1) - 2*hr**4*u1**2*tx**(-1) )
      MMv = MMv + SCD(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hl**2*
     &    ssp*m1**2*mt**2 + 16*pq*hl**2*ssp*m1**2*s + 8*pq*hl**2*ssp*
     &    m1**2*t1 + 8*pq*hl**2*ssp*m1**2*u1 + 8*pq*hl**2*ssp*m1**4 - 8
     &    *pq*hl**2*ssp*mt**2*t1 - 8*pq*hl**2*ssp*mt**2*u1 + 8*pq*hl**2
     &    *ssp*mt**4 - 8*pq*hl**2*ssp*t1*u1 - 16*pq*hr**2*ssp*m1**2*
     &    mt**2 + 16*pq*hr**2*ssp*m1**2*s + 8*pq*hr**2*ssp*m1**2*t1 + 8
     &    *pq*hr**2*ssp*m1**2*u1 + 8*pq*hr**2*ssp*m1**4 - 8*pq*hr**2*
     &    ssp*mt**2*t1 - 8*pq*hr**2*ssp*mt**2*u1 + 8*pq*hr**2*ssp*mt**4
     &     - 8*pq*hr**2*ssp*t1*u1 - 16*lq*hl**2*ssz*m1**2*mt**2*s*
     &    sz**(-1) + 8*lq*hl**2*ssz*m1**2*s*t1*sz**(-1) + 8*lq*hl**2*
     &    ssz*m1**2*s*u1*sz**(-1) + 16*lq*hl**2*ssz*m1**2*s**2*sz**(-1)
     &     + 8*lq*hl**2*ssz*m1**4*s*sz**(-1) - 8*lq*hl**2*ssz*mt**2*s*
     &    t1*sz**(-1) - 8*lq*hl**2*ssz*mt**2*s*u1*sz**(-1) + 8*lq*hl**2
     &    *ssz*mt**4*s*sz**(-1) - 8*lq*hl**2*ssz*s*t1*u1*sz**(-1) - 16*
     &    rq*hr**2*ssz*m1**2*mt**2*s*sz**(-1) + 8*rq*hr**2*ssz*m1**2*s*
     &    t1*sz**(-1) )
      MMv = MMv + SCD(1,4)*Nc*Cf*Pi*alphas*prefac * ( 8*rq*hr**2*ssz*
     &    m1**2*s*u1*sz**(-1) + 16*rq*hr**2*ssz*m1**2*s**2*sz**(-1) + 8
     &    *rq*hr**2*ssz*m1**4*s*sz**(-1) - 8*rq*hr**2*ssz*mt**2*s*t1*
     &    sz**(-1) - 8*rq*hr**2*ssz*mt**2*s*u1*sz**(-1) + 8*rq*hr**2*
     &    ssz*mt**4*s*sz**(-1) - 8*rq*hr**2*ssz*s*t1*u1*sz**(-1) - 16*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt*s**2*s1**(-1) - 16*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s**2*s2**(-1) - 8*hl**2*hr**2*mt**2*
     &    s**2*tx**(-1) - 4*hl**4*m1**2*mt**2*s*tx**(-1) + 2*hl**4*
     &    m1**2*s*t1*tx**(-1) + 2*hl**4*m1**2*s*u1*tx**(-1) + 4*hl**4*
     &    m1**2*s**2*tx**(-1) + 2*hl**4*m1**4*s*tx**(-1) - 2*hl**4*
     &    mt**2*s*t1*tx**(-1) - 2*hl**4*mt**2*s*u1*tx**(-1) + 2*hl**4*
     &    mt**4*s*tx**(-1) - 2*hl**4*s*t1*u1*tx**(-1) - 4*hr**4*m1**2*
     &    mt**2*s*tx**(-1) + 2*hr**4*m1**2*s*t1*tx**(-1) + 2*hr**4*
     &    m1**2*s*u1*tx**(-1) + 4*hr**4*m1**2*s**2*tx**(-1) + 2*hr**4*
     &    m1**4*s*tx**(-1) - 2*hr**4*mt**2*s*t1*tx**(-1) - 2*hr**4*
     &    mt**2*s*u1*tx**(-1) )
      MMv = MMv + SCD(1,4)*Nc*Cf*Pi*alphas*prefac * ( 2*hr**4*mt**4*s*
     &    tx**(-1) - 2*hr**4*s*t1*u1*tx**(-1) )

c               the phase space except for 1/s**2 
      HH_QBV = MMv / ( 16.D0 * Pi )

c               the averaging factors
      HH_QBV = HH_QBV /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      HH_QBV = HH_QBV * (m1+m2)**2/4.D0

      end
