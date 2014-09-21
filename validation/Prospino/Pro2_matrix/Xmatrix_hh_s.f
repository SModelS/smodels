c --------------------------------------------------------------------
      real*8 function HH_QBSY(massin,C,mu_tgb,alphas)

      implicit none 

      real*8     massin(1:30),Pi,prefac,Nc,Cf,C(1:20),mu_tgb,alphas
     &          ,mz,mh1,mh2,mt
     &          ,ssp,ssz,hl,hr 
     &          ,h1,h2,lambda1,lambda2
     &          ,lq,rq,pq
     &          ,lq2,rq2,pq2
     &          ,s,sz,s1,s2,m1,m2,t2,u2,tx,pt2,sqrt2,epsbw
     &          ,MMb,tag_mb,xs2b,mg,msb1,msb2,mb,B02
     &          ,tag_delta_mb
      logical lcheck_deltab
      parameter ( lcheck_deltab = .false. )

      Pi= 4.D0 * atan(1.D0)
      prefac = 1.D0/(16.D0*Pi**2)
      sqrt2 = sqrt(2.D0) 

      Nc = 3.D0 
      Cf = 4.D0/3.D0

      s   = massin(1)
      t2  = massin(2)
      m1  = massin(6)
      m2  = massin(6)
      mt  = massin(7)
      mz  = massin(8)
      mh1 = massin(9)
      mh2 = massin(10)
      mg   = massin(15)
      msb1 = massin(16)
      msb2 = massin(17)
      epsbw = massin(26)

      lq      = C(1)
      rq      = C(2)
      pq      = C(3)

c               hl=C(4) is the bottom Yukawa mb*tgb/mw
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

c               resum the bottom Yukawa contribution
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

c               note that MM proportional to C(4)^2 in large tan(beta) limit 
      tag_delta_mb = 1.D0/(1.D0+tag_mb)**4

c               want to compute K_SUSY as a correction to 2HDM
      tag_delta_mb = tag_delta_mb - 1.D0

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
      HH_QBSY = MMb / ( 16.D0 * Pi ) * tag_delta_mb

c               the averaging factors
      HH_QBSY = HH_QBSY /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      HH_QBSY = HH_QBSY * (m1+m2)**2/4.D0

      end

c --------------------------------------------------------------------
      real*8 function HH_QBS(massin,C)

      implicit none 

      integer    i1,i2
      real*8     massin(1:30),C(1:20),Pi,sqrt2,Nc,Cf,alphas
     &          ,logqr,prefac
     &          ,h1,h2,hl,hr,lambda1,lambda2
     &          ,lq,lq2,rq,rq2,pq,pq2,ssp,ssz
     &          ,m1,m2,mt,mz,mh1,mh2
     &          ,s,t,tx,t1,t2,u1,u2,s1,s2,sz,pt2
     &          ,mg,msb1,msb2,mst1,mst2
     &          ,sb,cb,s2b,c2b,st,ct,s2t,c2t
     &          ,hss(2,2),hhss(2,2),yuk1,yuk2
     &          ,SCA(1:10),SCB(1:10,1:6),SCBP(10)
     &          ,SCC(1:20,1:9),SCD(1:10,1:8)
     &          ,MMs

      Pi     = 4.D0 * atan(1.D0)
      sqrt2  = sqrt(2.D0) 
      prefac = 1.D0/(16.D0*Pi**2)

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
      mg    = massin(15)
      msb1  = massin(16)
      msb2  = massin(17)
      mst1  = massin(18)
      mst2  = massin(19)
      s2b   = massin(27)
      s2t   = massin(28)

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
      
ctp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ctp for now...
      do i1 = 1,2
         do i2  = 1,2
            hhss(i1,i2) = 0.D0
            Hss(i1,i2)  = 0.D0
         end do
      end do
      yuk1 = 0.D0
      yuk2 = 0.D0

      lq2 = lq**2 
      rq2 = rq**2 
      pq2 = pq**2 

      sz = s-mz**2
      s1 = s-mh1**2
      s2 = s-mh2**2

      u2 = -s-t2
      tx = t2+m2**2-mt**2

      t1 = t2 
      u1 = u2 
      t  = t2+m2**2

      pt2 = ( t2*u2 - s*m2**2 )/s 

c               squark mixing angles
      c2t = sqrt( 1.D0-s2t**2 )
      st  = sqrt((1.D0+c2t)/2.D0)
      if (s2t<0.D0) st = -st
      ct  = sqrt((1.D0-c2t)/2.D0)
 
      c2b = sqrt( 1.D0-s2b**2 )
      sb  = sqrt((1.D0+c2b)/2.D0)
      if (s2b<0.D0) sb = -sb
      cb  = sqrt((1.D0-c2b)/2.D0)

c               the factorization/renormalization scale 
      logqr = log( massin(12)**2/m1**2 )

c               the scalar functions 
      call SCALAR_ARRAY_HH(massin,SCA,SCB,SCBP,SCC,SCD)
c      print*, " test 1 ",scb(1,1:2)
c      print*, " test 2 ",scb(2,3:5)
c      print*, " test 4 ",scb(4,4:6)
c      print*, " test 6 ",scb(6,2:3)
c      print*, " test 7 ",scb(7,2:3)

c               set gs=1 
      alphas = 1/(4.D0*Pi)

      MMs = 0.d0
c               insert the form output 

      MMs =
     &  + Nc*Cf*Pi*alphas*prefac*logqr * (  - 16*pq*hl**2*ssp*pt2*
     &    tx**(-1) - 16*pq*hr**2*ssp*pt2*tx**(-1) - 16*lq*hl**2*ssz*pt2
     &    *s*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) - 24*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s1**(-1) - 24*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-1)*
     &    s2**(-1) - 8*hl**2*hr**2*mt**2*s*tx**(-2) - 4*hl**4*pt2*s*
     &    tx**(-2) - 4*hr**4*pt2*s*tx**(-2) - 16*h1*h2*lambda1*lambda2*
     &    s*s1**(-1)*s2**(-1) - 8*h1**2*lambda1**2*s*s1**(-2) - 8*h2**2
     &    *lambda2**2*s*s2**(-2) )
      MMs = MMs + Nc*Cf*Pi*alphas*prefac * (  - 64*pq*lq*ssz*ssp*m1**2*
     &    sz**(-1) + 64*pq*lq*ssz*ssp*s**(-1)*t1*u1*sz**(-1) - 64*pq*rq
     &    *ssz*ssp*m1**2*sz**(-1) + 64*pq*rq*ssz*ssp*s**(-1)*t1*u1*
     &    sz**(-1) - 8*pq*hl**2*ssp*m1**2*tx**(-1) + 8*pq*hl**2*ssp*
     &    s**(-1)*t1*u1*tx**(-1) - 8*pq*hr**2*ssp*m1**2*tx**(-1) + 8*pq
     &    *hr**2*ssp*s**(-1)*t1*u1*tx**(-1) - 8*lq*hl**2*ssz*m1**2*s*
     &    tx**(-1)*sz**(-1) + 8*lq*hl**2*ssz*t1*u1*tx**(-1)*sz**(-1) - 
     &    8*rq*hr**2*ssz*m1**2*s*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*t1*
     &    u1*tx**(-1)*sz**(-1) - 32*ssz**2*lq2*m1**2*s*sz**(-2) + 32*
     &    ssz**2*lq2*t1*u1*sz**(-2) - 32*ssz**2*rq2*m1**2*s*sz**(-2) + 
     &    32*ssz**2*rq2*t1*u1*sz**(-2) - 64*ssp**2*pq2*m1**2*s**(-1) + 
     &    64*ssp**2*pq2*s**(-2)*t1*u1 )
      MMs = MMs + SCA(2)*Nc*Cf*Pi*alphas*prefac * (  - 128*pq*lq*ssz*
     &    ssp*s**(-2)*t1*u1*sz**(-1) - 32*pq*lq*ssz*ssp*s**(-2)*t1**2*
     &    sz**(-1) - 32*pq*lq*ssz*ssp*s**(-2)*u1**2*sz**(-1) - 16*pq*lq
     &    *ssz*ssp*s**(-1)*t1*sz**(-1) - 16*pq*lq*ssz*ssp*s**(-1)*u1*
     &    sz**(-1) - 128*pq*rq*ssz*ssp*s**(-2)*t1*u1*sz**(-1) - 32*pq*
     &    rq*ssz*ssp*s**(-2)*t1**2*sz**(-1) - 32*pq*rq*ssz*ssp*s**(-2)*
     &    u1**2*sz**(-1) - 16*pq*rq*ssz*ssp*s**(-1)*t1*sz**(-1) - 16*pq
     &    *rq*ssz*ssp*s**(-1)*u1*sz**(-1) + 8*pq*hl**2*ssp*pt2*m1**2*
     &    tx**(-2)*t**(-1) + 8*pq*hl**2*ssp*pt2*mt**2*tx**(-2)*t**(-1)
     &     + 8*pq*hl**2*ssp*pt2*t1*tx**(-2)*t**(-1) - 16*pq*hl**2*ssp*
     &    pt2*tx**(-2) - 16*pq*hl**2*ssp*s**(-2)*t1*u1*tx**(-1) - 4*pq*
     &    hl**2*ssp*s**(-2)*t1**2*tx**(-1) - 4*pq*hl**2*ssp*s**(-2)*
     &    u1**2*tx**(-1) - 2*pq*hl**2*ssp*s**(-1)*t1*tx**(-1) - 2*pq*
     &    hl**2*ssp*s**(-1)*u1*tx**(-1) + 8*pq*hr**2*ssp*pt2*m1**2*
     &    tx**(-2)*t**(-1) + 8*pq*hr**2*ssp*pt2*mt**2*tx**(-2)*t**(-1)
     &     + 8*pq*hr**2*ssp*pt2*t1*tx**(-2)*t**(-1) )
      MMs = MMs + SCA(2)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hr**2*ssp*
     &    pt2*tx**(-2) - 16*pq*hr**2*ssp*s**(-2)*t1*u1*tx**(-1) - 4*pq*
     &    hr**2*ssp*s**(-2)*t1**2*tx**(-1) - 4*pq*hr**2*ssp*s**(-2)*
     &    u1**2*tx**(-1) - 2*pq*hr**2*ssp*s**(-1)*t1*tx**(-1) - 2*pq*
     &    hr**2*ssp*s**(-1)*u1*tx**(-1) + 8*lq*hl**2*ssz*pt2*m1**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) + 8*lq*hl**2*ssz*pt2*mt**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) + 8*lq*hl**2*ssz*pt2*s*t1*tx**(-2)*
     &    t**(-1)*sz**(-1) - 16*lq*hl**2*ssz*pt2*s*tx**(-2)*sz**(-1) - 
     &    16*lq*hl**2*ssz*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 4*lq*hl**2*
     &    ssz*s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*lq*hl**2*ssz*s**(-1)*
     &    u1**2*tx**(-1)*sz**(-1) - 2*lq*hl**2*ssz*t1*tx**(-1)*sz**(-1)
     &     - 2*lq*hl**2*ssz*u1*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*
     &    m1**2*s*tx**(-2)*t**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*mt**2*
     &    s*tx**(-2)*t**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*s*t1*
     &    tx**(-2)*t**(-1)*sz**(-1) - 16*rq*hr**2*ssz*pt2*s*tx**(-2)*
     &    sz**(-1) )
      MMs = MMs + SCA(2)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*ssz*
     &    s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*s**(-1)*
     &    t1**2*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*s**(-1)*u1**2*
     &    tx**(-1)*sz**(-1) - 2*rq*hr**2*ssz*t1*tx**(-1)*sz**(-1) - 2*
     &    rq*hr**2*ssz*u1*tx**(-1)*sz**(-1) - 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mt**(-1)*s*tx**(-2)*s1**(-1) + 16*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*t**(-1)*s1**(-1) - 8*
     &    hl*hr*h1*lambda1*sqrt2**(-1)*mt**(-1)*s*t1*tx**(-2)*s1**(-1)
     &     + 16*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*
     &    s1**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*tx**(-2)*
     &    s1**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**(-1)*s*
     &    tx**(-2)*s2**(-1) + 16*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*
     &    s*tx**(-2)*t**(-1)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mt**(-1)*s*t1*tx**(-2)*s2**(-1) + 16*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*s2**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*tx**(-2)*s2**(-1) )
      MMs = MMs + SCA(2)*Nc*Cf*Pi*alphas*prefac * ( 8*hl**2*hr**2*m1**2
     &    *mt**2*s*tx**(-3)*t**(-1) - 4*hl**2*hr**2*m1**2*s*tx**(-3) + 
     &    8*hl**2*hr**2*mt**2*s*t1*tx**(-3)*t**(-1) - 4*hl**2*hr**2*
     &    mt**2*s*tx**(-3) - 4*hl**2*hr**2*s*t1*tx**(-3) + 2*hl**4*pt2*
     &    m1**2*s*tx**(-3)*t**(-1) + 2*hl**4*pt2*mt**2*s*tx**(-3)*
     &    t**(-1) + 2*hl**4*pt2*s*t1*tx**(-3)*t**(-1) - 4*hl**4*pt2*s*
     &    tx**(-3) + 2*hr**4*pt2*m1**2*s*tx**(-3)*t**(-1) + 2*hr**4*pt2
     &    *mt**2*s*tx**(-3)*t**(-1) + 2*hr**4*pt2*s*t1*tx**(-3)*t**(-1)
     &     - 4*hr**4*pt2*s*tx**(-3) - 64*ssz**2*lq2*s**(-1)*t1*u1*
     &    sz**(-2) - 16*ssz**2*lq2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*
     &    lq2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*lq2*t1*sz**(-2) - 8*
     &    ssz**2*lq2*u1*sz**(-2) - 64*ssz**2*rq2*s**(-1)*t1*u1*sz**(-2)
     &     - 16*ssz**2*rq2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*rq2*
     &    s**(-1)*u1**2*sz**(-2) - 8*ssz**2*rq2*t1*sz**(-2) - 8*ssz**2*
     &    rq2*u1*sz**(-2) - 128*ssp**2*pq2*s**(-3)*t1*u1 - 32*ssp**2*
     &    pq2*s**(-3)*t1**2 )
      MMs = MMs + SCA(2)*Nc*Cf*Pi*alphas*prefac * (  - 32*ssp**2*pq2*
     &    s**(-3)*u1**2 - 16*ssp**2*pq2*s**(-2)*t1 - 16*ssp**2*pq2*
     &    s**(-2)*u1 )
      MMs = MMs + SCA(3)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hl**2*ssp*
     &    pt2*c2t*m1**2*tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*c2t*mt**2
     &    *tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*c2t*t1*tx**(-2)*
     &    t**(-1) - 4*pq*hl**2*ssp*pt2*m1**2*tx**(-2)*t**(-1) - 4*pq*
     &    hl**2*ssp*pt2*mt**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*t1*
     &    tx**(-2)*t**(-1) + 8*pq*hl**2*ssp*pt2*tx**(-2) + 4*pq*hr**2*
     &    ssp*pt2*c2t*m1**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*c2t*
     &    mt**2*tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*c2t*t1*tx**(-2)*
     &    t**(-1) - 4*pq*hr**2*ssp*pt2*m1**2*tx**(-2)*t**(-1) - 4*pq*
     &    hr**2*ssp*pt2*mt**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*t1*
     &    tx**(-2)*t**(-1) + 8*pq*hr**2*ssp*pt2*tx**(-2) - 4*lq*hl**2*
     &    ssz*pt2*c2t*m1**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*
     &    ssz*pt2*c2t*mt**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*
     &    ssz*pt2*c2t*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*
     &    pt2*m1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) )
      MMs = MMs + SCA(3)*Nc*Cf*Pi*alphas*prefac * (  - 4*lq*hl**2*ssz*
     &    pt2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 8*lq*hl**2*ssz*pt2*s*
     &    tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*m1**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*mt**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*s*t1*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*m1**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mt**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*s*t1*tx**(-2)*t**(-1)*
     &    sz**(-1) + 8*rq*hr**2*ssz*pt2*s*tx**(-2)*sz**(-1) + 4*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**2*mt**(-1)*s*tx**(-2)*s1**(-1) - 8
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*t**(-1)*
     &    s1**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt**(-1)*s*t1*
     &    tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*
     &    tx**(-2)*t**(-1)*s1**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*tx**(-2)*s1**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*
     &    mt**(-1)*s*tx**(-2)*s2**(-1) )
      MMs = MMs + SCA(3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*t**(-1)*s2**(-1) + 4*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt**(-1)*s*t1*tx**(-2)*s2**(-1)
     &     - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*
     &    s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-2)*
     &    s2**(-1) - 4*hl**2*hr**2*m1**2*mt**2*s*tx**(-3)*t**(-1) + 2*
     &    hl**2*hr**2*m1**2*s*tx**(-3) - 4*hl**2*hr**2*mt**2*s*t1*
     &    tx**(-3)*t**(-1) + 2*hl**2*hr**2*mt**2*s*tx**(-3) + 2*hl**2*
     &    hr**2*s*t1*tx**(-3) - hl**4*pt2*c2t*m1**2*s*tx**(-3)*t**(-1)
     &     + hl**4*pt2*c2t*mt**2*s*tx**(-3)*t**(-1) - hl**4*pt2*c2t*s*
     &    t1*tx**(-3)*t**(-1) - hl**4*pt2*m1**2*s*tx**(-3)*t**(-1) - 
     &    hl**4*pt2*mt**2*s*tx**(-3)*t**(-1) - hl**4*pt2*s*t1*tx**(-3)*
     &    t**(-1) + 2*hl**4*pt2*s*tx**(-3) + hr**4*pt2*c2t*m1**2*s*
     &    tx**(-3)*t**(-1) - hr**4*pt2*c2t*mt**2*s*tx**(-3)*t**(-1) + 
     &    hr**4*pt2*c2t*s*t1*tx**(-3)*t**(-1) - hr**4*pt2*m1**2*s*
     &    tx**(-3)*t**(-1) )
      MMs = MMs + SCA(3)*Nc*Cf*Pi*alphas*prefac * (  - hr**4*pt2*mt**2*
     &    s*tx**(-3)*t**(-1) - hr**4*pt2*s*t1*tx**(-3)*t**(-1) + 2*
     &    hr**4*pt2*s*tx**(-3) )
      MMs = MMs + SCA(4)*Nc*Cf*Pi*alphas*prefac * ( 4*pq*hl**2*ssp*pt2*
     &    c2t*m1**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*c2t*mt**2*
     &    tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*c2t*t1*tx**(-2)*t**(-1)
     &     - 4*pq*hl**2*ssp*pt2*m1**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp
     &    *pt2*mt**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*t1*tx**(-2)*
     &    t**(-1) + 8*pq*hl**2*ssp*pt2*tx**(-2) - 4*pq*hr**2*ssp*pt2*
     &    c2t*m1**2*tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*c2t*mt**2*
     &    tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*c2t*t1*tx**(-2)*t**(-1)
     &     - 4*pq*hr**2*ssp*pt2*m1**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp
     &    *pt2*mt**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*t1*tx**(-2)*
     &    t**(-1) + 8*pq*hr**2*ssp*pt2*tx**(-2) + 4*lq*hl**2*ssz*pt2*
     &    c2t*m1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*
     &    c2t*mt**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*
     &    c2t*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*m1**2
     &    *s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*mt**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) )
      MMs = MMs + SCA(4)*Nc*Cf*Pi*alphas*prefac * (  - 4*lq*hl**2*ssz*
     &    pt2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 8*lq*hl**2*ssz*pt2*s*
     &    tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*m1**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*mt**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*s*t1*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*m1**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mt**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*s*t1*tx**(-2)*t**(-1)*
     &    sz**(-1) + 8*rq*hr**2*ssz*pt2*s*tx**(-2)*sz**(-1) + 4*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*m1**2*mt**(-1)*s*tx**(-2)*s1**(-1) - 8
     &    *hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*t**(-1)*
     &    s1**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt**(-1)*s*t1*
     &    tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*
     &    tx**(-2)*t**(-1)*s1**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt
     &    *s*tx**(-2)*s1**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*
     &    mt**(-1)*s*tx**(-2)*s2**(-1) )
      MMs = MMs + SCA(4)*Nc*Cf*Pi*alphas*prefac * (  - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*t**(-1)*s2**(-1) + 4*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt**(-1)*s*t1*tx**(-2)*s2**(-1)
     &     - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*t**(-1)*
     &    s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*tx**(-2)*
     &    s2**(-1) - 4*hl**2*hr**2*m1**2*mt**2*s*tx**(-3)*t**(-1) + 2*
     &    hl**2*hr**2*m1**2*s*tx**(-3) - 4*hl**2*hr**2*mt**2*s*t1*
     &    tx**(-3)*t**(-1) + 2*hl**2*hr**2*mt**2*s*tx**(-3) + 2*hl**2*
     &    hr**2*s*t1*tx**(-3) + hl**4*pt2*c2t*m1**2*s*tx**(-3)*t**(-1)
     &     - hl**4*pt2*c2t*mt**2*s*tx**(-3)*t**(-1) + hl**4*pt2*c2t*s*
     &    t1*tx**(-3)*t**(-1) - hl**4*pt2*m1**2*s*tx**(-3)*t**(-1) - 
     &    hl**4*pt2*mt**2*s*tx**(-3)*t**(-1) - hl**4*pt2*s*t1*tx**(-3)*
     &    t**(-1) + 2*hl**4*pt2*s*tx**(-3) - hr**4*pt2*c2t*m1**2*s*
     &    tx**(-3)*t**(-1) + hr**4*pt2*c2t*mt**2*s*tx**(-3)*t**(-1) - 
     &    hr**4*pt2*c2t*s*t1*tx**(-3)*t**(-1) - hr**4*pt2*m1**2*s*
     &    tx**(-3)*t**(-1) )
      MMs = MMs + SCA(4)*Nc*Cf*Pi*alphas*prefac * (  - hr**4*pt2*mt**2*
     &    s*tx**(-3)*t**(-1) - hr**4*pt2*s*t1*tx**(-3)*t**(-1) + 2*
     &    hr**4*pt2*s*tx**(-3) )
      MMs = MMs + SCA(5)*Nc*Cf*Pi*alphas*prefac * ( 128*pq*lq*ssz*ssp*
     &    cb**2*s**(-2)*t1*u1*sz**(-1) + 32*pq*lq*ssz*ssp*cb**2*s**(-2)
     &    *t1**2*sz**(-1) + 32*pq*lq*ssz*ssp*cb**2*s**(-2)*u1**2*
     &    sz**(-1) + 16*pq*lq*ssz*ssp*cb**2*s**(-1)*t1*sz**(-1) + 16*pq
     &    *lq*ssz*ssp*cb**2*s**(-1)*u1*sz**(-1) - 128*pq*rq*ssz*ssp*
     &    cb**2*s**(-2)*t1*u1*sz**(-1) - 32*pq*rq*ssz*ssp*cb**2*s**(-2)
     &    *t1**2*sz**(-1) - 32*pq*rq*ssz*ssp*cb**2*s**(-2)*u1**2*
     &    sz**(-1) - 16*pq*rq*ssz*ssp*cb**2*s**(-1)*t1*sz**(-1) - 16*pq
     &    *rq*ssz*ssp*cb**2*s**(-1)*u1*sz**(-1) + 128*pq*rq*ssz*ssp*
     &    s**(-2)*t1*u1*sz**(-1) + 32*pq*rq*ssz*ssp*s**(-2)*t1**2*
     &    sz**(-1) + 32*pq*rq*ssz*ssp*s**(-2)*u1**2*sz**(-1) + 16*pq*rq
     &    *ssz*ssp*s**(-1)*t1*sz**(-1) + 16*pq*rq*ssz*ssp*s**(-1)*u1*
     &    sz**(-1) + 16*pq*hl**2*ssp*cb**2*s**(-2)*t1*u1*tx**(-1) + 4*
     &    pq*hl**2*ssp*cb**2*s**(-2)*t1**2*tx**(-1) + 4*pq*hl**2*ssp*
     &    cb**2*s**(-2)*u1**2*tx**(-1) + 2*pq*hl**2*ssp*cb**2*s**(-1)*
     &    t1*tx**(-1) )
      MMs = MMs + SCA(5)*Nc*Cf*Pi*alphas*prefac * ( 2*pq*hl**2*ssp*
     &    cb**2*s**(-1)*u1*tx**(-1) - 16*pq*hr**2*ssp*cb**2*s**(-2)*t1*
     &    u1*tx**(-1) - 4*pq*hr**2*ssp*cb**2*s**(-2)*t1**2*tx**(-1) - 4
     &    *pq*hr**2*ssp*cb**2*s**(-2)*u1**2*tx**(-1) - 2*pq*hr**2*ssp*
     &    cb**2*s**(-1)*t1*tx**(-1) - 2*pq*hr**2*ssp*cb**2*s**(-1)*u1*
     &    tx**(-1) + 16*pq*hr**2*ssp*s**(-2)*t1*u1*tx**(-1) + 4*pq*
     &    hr**2*ssp*s**(-2)*t1**2*tx**(-1) + 4*pq*hr**2*ssp*s**(-2)*
     &    u1**2*tx**(-1) + 2*pq*hr**2*ssp*s**(-1)*t1*tx**(-1) + 2*pq*
     &    hr**2*ssp*s**(-1)*u1*tx**(-1) + 16*lq*hl**2*ssz*cb**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) + 4*lq*hl**2*ssz*cb**2*s**(-1)*t1**2
     &    *tx**(-1)*sz**(-1) + 4*lq*hl**2*ssz*cb**2*s**(-1)*u1**2*
     &    tx**(-1)*sz**(-1) + 2*lq*hl**2*ssz*cb**2*t1*tx**(-1)*sz**(-1)
     &     + 2*lq*hl**2*ssz*cb**2*u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*
     &    ssz*cb**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*
     &    cb**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*cb**2*
     &    s**(-1)*u1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCA(5)*Nc*Cf*Pi*alphas*prefac * (  - 2*rq*hr**2*ssz*
     &    cb**2*t1*tx**(-1)*sz**(-1) - 2*rq*hr**2*ssz*cb**2*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hr**2*ssz*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     + 4*rq*hr**2*ssz*s**(-1)*t1**2*tx**(-1)*sz**(-1) + 4*rq*
     &    hr**2*ssz*s**(-1)*u1**2*tx**(-1)*sz**(-1) + 2*rq*hr**2*ssz*t1
     &    *tx**(-1)*sz**(-1) + 2*rq*hr**2*ssz*u1*tx**(-1)*sz**(-1) + 64
     &    *ssz**2*lq2*cb**2*s**(-1)*t1*u1*sz**(-2) + 16*ssz**2*lq2*
     &    cb**2*s**(-1)*t1**2*sz**(-2) + 16*ssz**2*lq2*cb**2*s**(-1)*
     &    u1**2*sz**(-2) + 8*ssz**2*lq2*cb**2*t1*sz**(-2) + 8*ssz**2*
     &    lq2*cb**2*u1*sz**(-2) - 64*ssz**2*rq2*cb**2*s**(-1)*t1*u1*
     &    sz**(-2) - 16*ssz**2*rq2*cb**2*s**(-1)*t1**2*sz**(-2) - 16*
     &    ssz**2*rq2*cb**2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*rq2*cb**2*
     &    t1*sz**(-2) - 8*ssz**2*rq2*cb**2*u1*sz**(-2) + 64*ssz**2*rq2*
     &    s**(-1)*t1*u1*sz**(-2) + 16*ssz**2*rq2*s**(-1)*t1**2*sz**(-2)
     &     + 16*ssz**2*rq2*s**(-1)*u1**2*sz**(-2) + 8*ssz**2*rq2*t1*
     &    sz**(-2) )
      MMs = MMs + SCA(5)*Nc*Cf*Pi*alphas*prefac * ( 8*ssz**2*rq2*u1*
     &    sz**(-2) + 64*ssp**2*pq2*s**(-3)*t1*u1 + 16*ssp**2*pq2*
     &    s**(-3)*t1**2 + 16*ssp**2*pq2*s**(-3)*u1**2 + 8*ssp**2*pq2*
     &    s**(-2)*t1 + 8*ssp**2*pq2*s**(-2)*u1 )
      MMs = MMs + SCA(6)*Nc*Cf*Pi*alphas*prefac * (  - 128*pq*lq*ssz*
     &    ssp*cb**2*s**(-2)*t1*u1*sz**(-1) - 32*pq*lq*ssz*ssp*cb**2*
     &    s**(-2)*t1**2*sz**(-1) - 32*pq*lq*ssz*ssp*cb**2*s**(-2)*u1**2
     &    *sz**(-1) - 16*pq*lq*ssz*ssp*cb**2*s**(-1)*t1*sz**(-1) - 16*
     &    pq*lq*ssz*ssp*cb**2*s**(-1)*u1*sz**(-1) + 128*pq*lq*ssz*ssp*
     &    s**(-2)*t1*u1*sz**(-1) + 32*pq*lq*ssz*ssp*s**(-2)*t1**2*
     &    sz**(-1) + 32*pq*lq*ssz*ssp*s**(-2)*u1**2*sz**(-1) + 16*pq*lq
     &    *ssz*ssp*s**(-1)*t1*sz**(-1) + 16*pq*lq*ssz*ssp*s**(-1)*u1*
     &    sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*s**(-2)*t1*u1*sz**(-1) + 
     &    32*pq*rq*ssz*ssp*cb**2*s**(-2)*t1**2*sz**(-1) + 32*pq*rq*ssz*
     &    ssp*cb**2*s**(-2)*u1**2*sz**(-1) + 16*pq*rq*ssz*ssp*cb**2*
     &    s**(-1)*t1*sz**(-1) + 16*pq*rq*ssz*ssp*cb**2*s**(-1)*u1*
     &    sz**(-1) - 16*pq*hl**2*ssp*cb**2*s**(-2)*t1*u1*tx**(-1) - 4*
     &    pq*hl**2*ssp*cb**2*s**(-2)*t1**2*tx**(-1) - 4*pq*hl**2*ssp*
     &    cb**2*s**(-2)*u1**2*tx**(-1) - 2*pq*hl**2*ssp*cb**2*s**(-1)*
     &    t1*tx**(-1) )
      MMs = MMs + SCA(6)*Nc*Cf*Pi*alphas*prefac * (  - 2*pq*hl**2*ssp*
     &    cb**2*s**(-1)*u1*tx**(-1) + 16*pq*hl**2*ssp*s**(-2)*t1*u1*
     &    tx**(-1) + 4*pq*hl**2*ssp*s**(-2)*t1**2*tx**(-1) + 4*pq*hl**2
     &    *ssp*s**(-2)*u1**2*tx**(-1) + 2*pq*hl**2*ssp*s**(-1)*t1*
     &    tx**(-1) + 2*pq*hl**2*ssp*s**(-1)*u1*tx**(-1) + 16*pq*hr**2*
     &    ssp*cb**2*s**(-2)*t1*u1*tx**(-1) + 4*pq*hr**2*ssp*cb**2*
     &    s**(-2)*t1**2*tx**(-1) + 4*pq*hr**2*ssp*cb**2*s**(-2)*u1**2*
     &    tx**(-1) + 2*pq*hr**2*ssp*cb**2*s**(-1)*t1*tx**(-1) + 2*pq*
     &    hr**2*ssp*cb**2*s**(-1)*u1*tx**(-1) - 16*lq*hl**2*ssz*cb**2*
     &    s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 4*lq*hl**2*ssz*cb**2*
     &    s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*lq*hl**2*ssz*cb**2*
     &    s**(-1)*u1**2*tx**(-1)*sz**(-1) - 2*lq*hl**2*ssz*cb**2*t1*
     &    tx**(-1)*sz**(-1) - 2*lq*hl**2*ssz*cb**2*u1*tx**(-1)*sz**(-1)
     &     + 16*lq*hl**2*ssz*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 4*lq*
     &    hl**2*ssz*s**(-1)*t1**2*tx**(-1)*sz**(-1) + 4*lq*hl**2*ssz*
     &    s**(-1)*u1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCA(6)*Nc*Cf*Pi*alphas*prefac * ( 2*lq*hl**2*ssz*t1*
     &    tx**(-1)*sz**(-1) + 2*lq*hl**2*ssz*u1*tx**(-1)*sz**(-1) + 16*
     &    rq*hr**2*ssz*cb**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 4*rq*
     &    hr**2*ssz*cb**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) + 4*rq*hr**2*
     &    ssz*cb**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) + 2*rq*hr**2*ssz*
     &    cb**2*t1*tx**(-1)*sz**(-1) + 2*rq*hr**2*ssz*cb**2*u1*tx**(-1)
     &    *sz**(-1) - 64*ssz**2*lq2*cb**2*s**(-1)*t1*u1*sz**(-2) - 16*
     &    ssz**2*lq2*cb**2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*lq2*cb**2
     &    *s**(-1)*u1**2*sz**(-2) - 8*ssz**2*lq2*cb**2*t1*sz**(-2) - 8*
     &    ssz**2*lq2*cb**2*u1*sz**(-2) + 64*ssz**2*lq2*s**(-1)*t1*u1*
     &    sz**(-2) + 16*ssz**2*lq2*s**(-1)*t1**2*sz**(-2) + 16*ssz**2*
     &    lq2*s**(-1)*u1**2*sz**(-2) + 8*ssz**2*lq2*t1*sz**(-2) + 8*
     &    ssz**2*lq2*u1*sz**(-2) + 64*ssz**2*rq2*cb**2*s**(-1)*t1*u1*
     &    sz**(-2) + 16*ssz**2*rq2*cb**2*s**(-1)*t1**2*sz**(-2) + 16*
     &    ssz**2*rq2*cb**2*s**(-1)*u1**2*sz**(-2) + 8*ssz**2*rq2*cb**2*
     &    t1*sz**(-2) )
      MMs = MMs + SCA(6)*Nc*Cf*Pi*alphas*prefac * ( 8*ssz**2*rq2*cb**2*
     &    u1*sz**(-2) + 64*ssp**2*pq2*s**(-3)*t1*u1 + 16*ssp**2*pq2*
     &    s**(-3)*t1**2 + 16*ssp**2*pq2*s**(-3)*u1**2 + 8*ssp**2*pq2*
     &    s**(-2)*t1 + 8*ssp**2*pq2*s**(-2)*u1 )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * ( 64*pq*lq*ssz*ssp*
     &    pt2*c2b*sz**(-1) - 32*pq*lq*ssz*ssp*pt2*sz**(-1) + 128*pq*lq*
     &    ssz*ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) - 128*pq*lq*ssz*
     &    ssp*cb**2*m1**2*msb1**2*s**(-1)*sz**(-1) + 32*pq*lq*ssz*ssp*
     &    cb**2*mg**2*s**(-2)*t1**2*sz**(-1) + 32*pq*lq*ssz*ssp*cb**2*
     &    mg**2*s**(-2)*u1**2*sz**(-1) + 16*pq*lq*ssz*ssp*cb**2*mg**2*
     &    s**(-1)*t1*sz**(-1) + 16*pq*lq*ssz*ssp*cb**2*mg**2*s**(-1)*u1
     &    *sz**(-1) - 32*pq*lq*ssz*ssp*cb**2*msb1**2*s**(-2)*t1**2*
     &    sz**(-1) - 32*pq*lq*ssz*ssp*cb**2*msb1**2*s**(-2)*u1**2*
     &    sz**(-1) - 16*pq*lq*ssz*ssp*cb**2*msb1**2*s**(-1)*t1*sz**(-1)
     &     - 16*pq*lq*ssz*ssp*cb**2*msb1**2*s**(-1)*u1*sz**(-1) - 64*pq
     &    *rq*ssz*ssp*pt2*c2b*sz**(-1) - 32*pq*rq*ssz*ssp*pt2*sz**(-1)
     &     - 128*pq*rq*ssz*ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) + 128
     &    *pq*rq*ssz*ssp*cb**2*m1**2*msb1**2*s**(-1)*sz**(-1) - 32*pq*
     &    rq*ssz*ssp*cb**2*mg**2*s**(-2)*t1**2*sz**(-1) - 32*pq*rq*ssz*
     &    ssp*cb**2*mg**2*s**(-2)*u1**2*sz**(-1) )
ctp      print*, " dummy 1 "
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*rq*ssz*
     &    ssp*cb**2*mg**2*s**(-1)*t1*sz**(-1) - 16*pq*rq*ssz*ssp*cb**2*
     &    mg**2*s**(-1)*u1*sz**(-1) + 32*pq*rq*ssz*ssp*cb**2*msb1**2*
     &    s**(-2)*t1**2*sz**(-1) + 32*pq*rq*ssz*ssp*cb**2*msb1**2*
     &    s**(-2)*u1**2*sz**(-1) + 16*pq*rq*ssz*ssp*cb**2*msb1**2*
     &    s**(-1)*t1*sz**(-1) + 16*pq*rq*ssz*ssp*cb**2*msb1**2*s**(-1)*
     &    u1*sz**(-1) + 128*pq*rq*ssz*ssp*m1**2*mg**2*s**(-1)*sz**(-1)
     &     - 128*pq*rq*ssz*ssp*m1**2*msb1**2*s**(-1)*sz**(-1) + 32*pq*
     &    rq*ssz*ssp*mg**2*s**(-2)*t1**2*sz**(-1) + 32*pq*rq*ssz*ssp*
     &    mg**2*s**(-2)*u1**2*sz**(-1) + 16*pq*rq*ssz*ssp*mg**2*s**(-1)
     &    *t1*sz**(-1) + 16*pq*rq*ssz*ssp*mg**2*s**(-1)*u1*sz**(-1) - 
     &    32*pq*rq*ssz*ssp*msb1**2*s**(-2)*t1**2*sz**(-1) - 32*pq*rq*
     &    ssz*ssp*msb1**2*s**(-2)*u1**2*sz**(-1) - 16*pq*rq*ssz*ssp*
     &    msb1**2*s**(-1)*t1*sz**(-1) - 16*pq*rq*ssz*ssp*msb1**2*
     &    s**(-1)*u1*sz**(-1) - 32*pq*hl*hr*ssp*sb*cb*mg*mt*s**(-1)*t1*
     &    tx**(-1) )
ctp      print*, " dummy 1 "
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * ( 32*pq*hl*hr*ssp*
     &    sb*cb*mg*mt*s**(-1)*u1*tx**(-1) + 16*pq*hl**2*ssp*pt2*c2b*
     &    tx**(-1) - 8*pq*hl**2*ssp*pt2*tx**(-1) + 16*pq*hl**2*ssp*
     &    cb**2*m1**2*mg**2*s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*cb**2*
     &    m1**2*msb1**2*s**(-1)*tx**(-1) + 4*pq*hl**2*ssp*cb**2*mg**2*
     &    s**(-2)*t1**2*tx**(-1) + 4*pq*hl**2*ssp*cb**2*mg**2*s**(-2)*
     &    u1**2*tx**(-1) + 2*pq*hl**2*ssp*cb**2*mg**2*s**(-1)*t1*
     &    tx**(-1) + 2*pq*hl**2*ssp*cb**2*mg**2*s**(-1)*u1*tx**(-1) - 4
     &    *pq*hl**2*ssp*cb**2*msb1**2*s**(-2)*t1**2*tx**(-1) - 4*pq*
     &    hl**2*ssp*cb**2*msb1**2*s**(-2)*u1**2*tx**(-1) - 2*pq*hl**2*
     &    ssp*cb**2*msb1**2*s**(-1)*t1*tx**(-1) - 2*pq*hl**2*ssp*cb**2*
     &    msb1**2*s**(-1)*u1*tx**(-1) - 16*pq*hr**2*ssp*pt2*c2b*
     &    tx**(-1) - 8*pq*hr**2*ssp*pt2*tx**(-1) - 16*pq*hr**2*ssp*
     &    cb**2*m1**2*mg**2*s**(-1)*tx**(-1) + 16*pq*hr**2*ssp*cb**2*
     &    m1**2*msb1**2*s**(-1)*tx**(-1) - 4*pq*hr**2*ssp*cb**2*mg**2*
     &    s**(-2)*t1**2*tx**(-1) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hr**2*ssp
     &    *cb**2*mg**2*s**(-2)*u1**2*tx**(-1) - 2*pq*hr**2*ssp*cb**2*
     &    mg**2*s**(-1)*t1*tx**(-1) - 2*pq*hr**2*ssp*cb**2*mg**2*
     &    s**(-1)*u1*tx**(-1) + 4*pq*hr**2*ssp*cb**2*msb1**2*s**(-2)*
     &    t1**2*tx**(-1) + 4*pq*hr**2*ssp*cb**2*msb1**2*s**(-2)*u1**2*
     &    tx**(-1) + 2*pq*hr**2*ssp*cb**2*msb1**2*s**(-1)*t1*tx**(-1)
     &     + 2*pq*hr**2*ssp*cb**2*msb1**2*s**(-1)*u1*tx**(-1) + 16*pq*
     &    hr**2*ssp*m1**2*mg**2*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp*
     &    m1**2*msb1**2*s**(-1)*tx**(-1) + 4*pq*hr**2*ssp*mg**2*s**(-2)
     &    *t1**2*tx**(-1) + 4*pq*hr**2*ssp*mg**2*s**(-2)*u1**2*tx**(-1)
     &     + 2*pq*hr**2*ssp*mg**2*s**(-1)*t1*tx**(-1) + 2*pq*hr**2*ssp*
     &    mg**2*s**(-1)*u1*tx**(-1) - 4*pq*hr**2*ssp*msb1**2*s**(-2)*
     &    t1**2*tx**(-1) - 4*pq*hr**2*ssp*msb1**2*s**(-2)*u1**2*
     &    tx**(-1) - 2*pq*hr**2*ssp*msb1**2*s**(-1)*t1*tx**(-1) - 2*pq*
     &    hr**2*ssp*msb1**2*s**(-1)*u1*tx**(-1) - 64*pq*h1*ssp*lambda1*
     &    sb*cb*sqrt2**(-1)*mg*s**(-1)*t1*s1**(-1) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * ( 64*pq*h1*ssp*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s**(-1)*u1*s1**(-1) - 64*pq*h2*
     &    ssp*lambda2*sb*cb*sqrt2**(-1)*mg*s**(-1)*t1*s2**(-1) + 64*pq*
     &    h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg*s**(-1)*u1*s2**(-1) - 16*
     &    lq*hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) + 16*lq*hl*hr*
     &    ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*pt2*
     &    c2b*s*tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2*msb1**2*tx**(-1)*
     &    sz**(-1) + 4*lq*hl**2*ssz*cb**2*mg**2*s**(-1)*t1**2*tx**(-1)*
     &    sz**(-1) + 4*lq*hl**2*ssz*cb**2*mg**2*s**(-1)*u1**2*tx**(-1)*
     &    sz**(-1) + 2*lq*hl**2*ssz*cb**2*mg**2*t1*tx**(-1)*sz**(-1) + 
     &    2*lq*hl**2*ssz*cb**2*mg**2*u1*tx**(-1)*sz**(-1) - 4*lq*hl**2*
     &    ssz*cb**2*msb1**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*lq*
     &    hl**2*ssz*cb**2*msb1**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) - 2*
     &    lq*hl**2*ssz*cb**2*msb1**2*t1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 2*lq*hl**2*ssz
     &    *cb**2*msb1**2*u1*tx**(-1)*sz**(-1) - 32*lq*h1*ssz*lambda1*sb
     &    *cb*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) + 32*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) - 32*lq*h2*
     &    ssz*lambda2*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) + 32*lq
     &    *h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 
     &    16*rq*hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) + 16*rq*hl*
     &    hr*ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*pt2
     &    *c2b*s*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**2*m1**2*msb1**2*tx**(-1)*
     &    sz**(-1) - 4*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1**2*tx**(-1)*
     &    sz**(-1) - 4*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*u1**2*tx**(-1)*
     &    sz**(-1) - 2*rq*hr**2*ssz*cb**2*mg**2*t1*tx**(-1)*sz**(-1) - 
     &    2*rq*hr**2*ssz*cb**2*mg**2*u1*tx**(-1)*sz**(-1) + 4*rq*hr**2*
     &    ssz*cb**2*msb1**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * ( 4*rq*hr**2*ssz*
     &    cb**2*msb1**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) + 2*rq*hr**2*
     &    ssz*cb**2*msb1**2*t1*tx**(-1)*sz**(-1) + 2*rq*hr**2*ssz*cb**2
     &    *msb1**2*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*m1**2*msb1**2*tx**(-1)*
     &    sz**(-1) + 4*rq*hr**2*ssz*mg**2*s**(-1)*t1**2*tx**(-1)*
     &    sz**(-1) + 4*rq*hr**2*ssz*mg**2*s**(-1)*u1**2*tx**(-1)*
     &    sz**(-1) + 2*rq*hr**2*ssz*mg**2*t1*tx**(-1)*sz**(-1) + 2*rq*
     &    hr**2*ssz*mg**2*u1*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*msb1**2
     &    *s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*rq*hr**2*ssz*msb1**2*
     &    s**(-1)*u1**2*tx**(-1)*sz**(-1) - 2*rq*hr**2*ssz*msb1**2*t1*
     &    tx**(-1)*sz**(-1) - 2*rq*hr**2*ssz*msb1**2*u1*tx**(-1)*
     &    sz**(-1) - 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*t1*
     &    sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg
     &    *u1*sz**(-1)*s1**(-1) - 32*rq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * ( 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 8*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 2*hl**2*hr**2*
     &    mt**2*s*tx**(-2) + 2*hl**4*pt2*c2b*s*tx**(-2) - hl**4*pt2*s*
     &    tx**(-2) - 2*hr**4*pt2*c2b*s*tx**(-2) - hr**4*pt2*s*tx**(-2)
     &     - 8*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1) - 4*h1**2*
     &    lambda1**2*s*s1**(-2) - 4*h2**2*lambda2**2*s*s2**(-2) + 32*
     &    ssz**2*lq2*pt2*c2b*s*sz**(-2) - 16*ssz**2*lq2*pt2*s*sz**(-2)
     &     + 64*ssz**2*lq2*cb**2*m1**2*mg**2*sz**(-2) - 64*ssz**2*lq2*
     &    cb**2*m1**2*msb1**2*sz**(-2) + 16*ssz**2*lq2*cb**2*mg**2*
     &    s**(-1)*t1**2*sz**(-2) + 16*ssz**2*lq2*cb**2*mg**2*s**(-1)*
     &    u1**2*sz**(-2) + 8*ssz**2*lq2*cb**2*mg**2*t1*sz**(-2) + 8*
     &    ssz**2*lq2*cb**2*mg**2*u1*sz**(-2) - 16*ssz**2*lq2*cb**2*
     &    msb1**2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*lq2*cb**2*msb1**2*
     &    s**(-1)*u1**2*sz**(-2) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 8*ssz**2*lq2*
     &    cb**2*msb1**2*t1*sz**(-2) - 8*ssz**2*lq2*cb**2*msb1**2*u1*
     &    sz**(-2) - 32*ssz**2*rq2*pt2*c2b*s*sz**(-2) - 16*ssz**2*rq2*
     &    pt2*s*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*mg**2*sz**(-2) + 
     &    64*ssz**2*rq2*cb**2*m1**2*msb1**2*sz**(-2) - 16*ssz**2*rq2*
     &    cb**2*mg**2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*rq2*cb**2*
     &    mg**2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*rq2*cb**2*mg**2*t1*
     &    sz**(-2) - 8*ssz**2*rq2*cb**2*mg**2*u1*sz**(-2) + 16*ssz**2*
     &    rq2*cb**2*msb1**2*s**(-1)*t1**2*sz**(-2) + 16*ssz**2*rq2*
     &    cb**2*msb1**2*s**(-1)*u1**2*sz**(-2) + 8*ssz**2*rq2*cb**2*
     &    msb1**2*t1*sz**(-2) + 8*ssz**2*rq2*cb**2*msb1**2*u1*sz**(-2)
     &     + 64*ssz**2*rq2*m1**2*mg**2*sz**(-2) - 64*ssz**2*rq2*m1**2*
     &    msb1**2*sz**(-2) + 16*ssz**2*rq2*mg**2*s**(-1)*t1**2*sz**(-2)
     &     + 16*ssz**2*rq2*mg**2*s**(-1)*u1**2*sz**(-2) + 8*ssz**2*rq2*
     &    mg**2*t1*sz**(-2) + 8*ssz**2*rq2*mg**2*u1*sz**(-2) - 16*
     &    ssz**2*rq2*msb1**2*s**(-1)*t1**2*sz**(-2) )
      MMs = MMs + SCB(1,1)*Nc*Cf*Pi*alphas*prefac * (  - 16*ssz**2*rq2*
     &    msb1**2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*rq2*msb1**2*t1*
     &    sz**(-2) - 8*ssz**2*rq2*msb1**2*u1*sz**(-2) - 32*ssp**2*pq2*
     &    pt2*s**(-1) + 64*ssp**2*pq2*m1**2*mg**2*s**(-2) - 64*ssp**2*
     &    pq2*m1**2*msb1**2*s**(-2) + 16*ssp**2*pq2*mg**2*s**(-3)*t1**2
     &     + 16*ssp**2*pq2*mg**2*s**(-3)*u1**2 + 8*ssp**2*pq2*mg**2*
     &    s**(-2)*t1 + 8*ssp**2*pq2*mg**2*s**(-2)*u1 - 16*ssp**2*pq2*
     &    msb1**2*s**(-3)*t1**2 - 16*ssp**2*pq2*msb1**2*s**(-3)*u1**2
     &     - 8*ssp**2*pq2*msb1**2*s**(-2)*t1 - 8*ssp**2*pq2*msb1**2*
     &    s**(-2)*u1 )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 64*pq*lq*ssz*
     &    ssp*pt2*c2b*sz**(-1) - 32*pq*lq*ssz*ssp*pt2*sz**(-1) - 128*pq
     &    *lq*ssz*ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) + 128*pq*lq*
     &    ssz*ssp*cb**2*m1**2*msb2**2*s**(-1)*sz**(-1) - 32*pq*lq*ssz*
     &    ssp*cb**2*mg**2*s**(-2)*t1**2*sz**(-1) - 32*pq*lq*ssz*ssp*
     &    cb**2*mg**2*s**(-2)*u1**2*sz**(-1) - 16*pq*lq*ssz*ssp*cb**2*
     &    mg**2*s**(-1)*t1*sz**(-1) - 16*pq*lq*ssz*ssp*cb**2*mg**2*
     &    s**(-1)*u1*sz**(-1) + 32*pq*lq*ssz*ssp*cb**2*msb2**2*s**(-2)*
     &    t1**2*sz**(-1) + 32*pq*lq*ssz*ssp*cb**2*msb2**2*s**(-2)*u1**2
     &    *sz**(-1) + 16*pq*lq*ssz*ssp*cb**2*msb2**2*s**(-1)*t1*
     &    sz**(-1) + 16*pq*lq*ssz*ssp*cb**2*msb2**2*s**(-1)*u1*sz**(-1)
     &     + 128*pq*lq*ssz*ssp*m1**2*mg**2*s**(-1)*sz**(-1) - 128*pq*lq
     &    *ssz*ssp*m1**2*msb2**2*s**(-1)*sz**(-1) + 32*pq*lq*ssz*ssp*
     &    mg**2*s**(-2)*t1**2*sz**(-1) + 32*pq*lq*ssz*ssp*mg**2*s**(-2)
     &    *u1**2*sz**(-1) + 16*pq*lq*ssz*ssp*mg**2*s**(-1)*t1*sz**(-1)
     &     + 16*pq*lq*ssz*ssp*mg**2*s**(-1)*u1*sz**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*pq*lq*ssz*
     &    ssp*msb2**2*s**(-2)*t1**2*sz**(-1) - 32*pq*lq*ssz*ssp*msb2**2
     &    *s**(-2)*u1**2*sz**(-1) - 16*pq*lq*ssz*ssp*msb2**2*s**(-1)*t1
     &    *sz**(-1) - 16*pq*lq*ssz*ssp*msb2**2*s**(-1)*u1*sz**(-1) + 64
     &    *pq*rq*ssz*ssp*pt2*c2b*sz**(-1) - 32*pq*rq*ssz*ssp*pt2*
     &    sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*m1**2*mg**2*s**(-1)*
     &    sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*m1**2*msb2**2*s**(-1)*
     &    sz**(-1) + 32*pq*rq*ssz*ssp*cb**2*mg**2*s**(-2)*t1**2*
     &    sz**(-1) + 32*pq*rq*ssz*ssp*cb**2*mg**2*s**(-2)*u1**2*
     &    sz**(-1) + 16*pq*rq*ssz*ssp*cb**2*mg**2*s**(-1)*t1*sz**(-1)
     &     + 16*pq*rq*ssz*ssp*cb**2*mg**2*s**(-1)*u1*sz**(-1) - 32*pq*
     &    rq*ssz*ssp*cb**2*msb2**2*s**(-2)*t1**2*sz**(-1) - 32*pq*rq*
     &    ssz*ssp*cb**2*msb2**2*s**(-2)*u1**2*sz**(-1) - 16*pq*rq*ssz*
     &    ssp*cb**2*msb2**2*s**(-1)*t1*sz**(-1) - 16*pq*rq*ssz*ssp*
     &    cb**2*msb2**2*s**(-1)*u1*sz**(-1) + 32*pq*hl*hr*ssp*sb*cb*mg*
     &    mt*s**(-1)*t1*tx**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*pq*hl*hr*
     &    ssp*sb*cb*mg*mt*s**(-1)*u1*tx**(-1) - 16*pq*hl**2*ssp*pt2*c2b
     &    *tx**(-1) - 8*pq*hl**2*ssp*pt2*tx**(-1) - 16*pq*hl**2*ssp*
     &    cb**2*m1**2*mg**2*s**(-1)*tx**(-1) + 16*pq*hl**2*ssp*cb**2*
     &    m1**2*msb2**2*s**(-1)*tx**(-1) - 4*pq*hl**2*ssp*cb**2*mg**2*
     &    s**(-2)*t1**2*tx**(-1) - 4*pq*hl**2*ssp*cb**2*mg**2*s**(-2)*
     &    u1**2*tx**(-1) - 2*pq*hl**2*ssp*cb**2*mg**2*s**(-1)*t1*
     &    tx**(-1) - 2*pq*hl**2*ssp*cb**2*mg**2*s**(-1)*u1*tx**(-1) + 4
     &    *pq*hl**2*ssp*cb**2*msb2**2*s**(-2)*t1**2*tx**(-1) + 4*pq*
     &    hl**2*ssp*cb**2*msb2**2*s**(-2)*u1**2*tx**(-1) + 2*pq*hl**2*
     &    ssp*cb**2*msb2**2*s**(-1)*t1*tx**(-1) + 2*pq*hl**2*ssp*cb**2*
     &    msb2**2*s**(-1)*u1*tx**(-1) + 16*pq*hl**2*ssp*m1**2*mg**2*
     &    s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*m1**2*msb2**2*s**(-1)*
     &    tx**(-1) + 4*pq*hl**2*ssp*mg**2*s**(-2)*t1**2*tx**(-1) + 4*pq
     &    *hl**2*ssp*mg**2*s**(-2)*u1**2*tx**(-1) + 2*pq*hl**2*ssp*
     &    mg**2*s**(-1)*t1*tx**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * ( 2*pq*hl**2*ssp*
     &    mg**2*s**(-1)*u1*tx**(-1) - 4*pq*hl**2*ssp*msb2**2*s**(-2)*
     &    t1**2*tx**(-1) - 4*pq*hl**2*ssp*msb2**2*s**(-2)*u1**2*
     &    tx**(-1) - 2*pq*hl**2*ssp*msb2**2*s**(-1)*t1*tx**(-1) - 2*pq*
     &    hl**2*ssp*msb2**2*s**(-1)*u1*tx**(-1) + 16*pq*hr**2*ssp*pt2*
     &    c2b*tx**(-1) - 8*pq*hr**2*ssp*pt2*tx**(-1) + 16*pq*hr**2*ssp*
     &    cb**2*m1**2*mg**2*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp*cb**2*
     &    m1**2*msb2**2*s**(-1)*tx**(-1) + 4*pq*hr**2*ssp*cb**2*mg**2*
     &    s**(-2)*t1**2*tx**(-1) + 4*pq*hr**2*ssp*cb**2*mg**2*s**(-2)*
     &    u1**2*tx**(-1) + 2*pq*hr**2*ssp*cb**2*mg**2*s**(-1)*t1*
     &    tx**(-1) + 2*pq*hr**2*ssp*cb**2*mg**2*s**(-1)*u1*tx**(-1) - 4
     &    *pq*hr**2*ssp*cb**2*msb2**2*s**(-2)*t1**2*tx**(-1) - 4*pq*
     &    hr**2*ssp*cb**2*msb2**2*s**(-2)*u1**2*tx**(-1) - 2*pq*hr**2*
     &    ssp*cb**2*msb2**2*s**(-1)*t1*tx**(-1) - 2*pq*hr**2*ssp*cb**2*
     &    msb2**2*s**(-1)*u1*tx**(-1) + 64*pq*h1*ssp*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*s**(-1)*t1*s1**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 64*pq*h1*ssp*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s**(-1)*u1*s1**(-1) + 64*pq*h2*
     &    ssp*lambda2*sb*cb*sqrt2**(-1)*mg*s**(-1)*t1*s2**(-1) - 64*pq*
     &    h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg*s**(-1)*u1*s2**(-1) + 16*
     &    lq*hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) - 16*lq*hl*hr*
     &    ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*pt2*
     &    c2b*s*tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2*msb2**2*tx**(-1)*
     &    sz**(-1) - 4*lq*hl**2*ssz*cb**2*mg**2*s**(-1)*t1**2*tx**(-1)*
     &    sz**(-1) - 4*lq*hl**2*ssz*cb**2*mg**2*s**(-1)*u1**2*tx**(-1)*
     &    sz**(-1) - 2*lq*hl**2*ssz*cb**2*mg**2*t1*tx**(-1)*sz**(-1) - 
     &    2*lq*hl**2*ssz*cb**2*mg**2*u1*tx**(-1)*sz**(-1) + 4*lq*hl**2*
     &    ssz*cb**2*msb2**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) + 4*lq*
     &    hl**2*ssz*cb**2*msb2**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) + 2*
     &    lq*hl**2*ssz*cb**2*msb2**2*t1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * ( 2*lq*hl**2*ssz*
     &    cb**2*msb2**2*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*m1**2*
     &    mg**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**2*msb2**2*
     &    tx**(-1)*sz**(-1) + 4*lq*hl**2*ssz*mg**2*s**(-1)*t1**2*
     &    tx**(-1)*sz**(-1) + 4*lq*hl**2*ssz*mg**2*s**(-1)*u1**2*
     &    tx**(-1)*sz**(-1) + 2*lq*hl**2*ssz*mg**2*t1*tx**(-1)*sz**(-1)
     &     + 2*lq*hl**2*ssz*mg**2*u1*tx**(-1)*sz**(-1) - 4*lq*hl**2*ssz
     &    *msb2**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*lq*hl**2*ssz*
     &    msb2**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) - 2*lq*hl**2*ssz*
     &    msb2**2*t1*tx**(-1)*sz**(-1) - 2*lq*hl**2*ssz*msb2**2*u1*
     &    tx**(-1)*sz**(-1) + 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg
     &    *t1*sz**(-1)*s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) + 32*lq*h2*ssz*lambda2*sb
     &    *cb*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) + 16*rq*hl*
     &    hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hl*hr*
     &    ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*pt2*
     &    c2b*s*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*pt2*s*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*cb**2*m1**2*msb2**2*tx**(-1)*
     &    sz**(-1) + 4*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1**2*tx**(-1)*
     &    sz**(-1) + 4*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*u1**2*tx**(-1)*
     &    sz**(-1) + 2*rq*hr**2*ssz*cb**2*mg**2*t1*tx**(-1)*sz**(-1) + 
     &    2*rq*hr**2*ssz*cb**2*mg**2*u1*tx**(-1)*sz**(-1) - 4*rq*hr**2*
     &    ssz*cb**2*msb2**2*s**(-1)*t1**2*tx**(-1)*sz**(-1) - 4*rq*
     &    hr**2*ssz*cb**2*msb2**2*s**(-1)*u1**2*tx**(-1)*sz**(-1) - 2*
     &    rq*hr**2*ssz*cb**2*msb2**2*t1*tx**(-1)*sz**(-1) - 2*rq*hr**2*
     &    ssz*cb**2*msb2**2*u1*tx**(-1)*sz**(-1) + 32*rq*h1*ssz*lambda1
     &    *sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) + 32*rq*h2*
     &    ssz*lambda2*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 8*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*tx**(-1)*s1**(-1) - 8*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*mt*s*tx**(-1)*s2**(-1) - 2*hl**2*hr**2*
     &    mt**2*s*tx**(-2) - 2*hl**4*pt2*c2b*s*tx**(-2) - hl**4*pt2*s*
     &    tx**(-2) + 2*hr**4*pt2*c2b*s*tx**(-2) - hr**4*pt2*s*tx**(-2)
     &     - 8*h1*h2*lambda1*lambda2*s*s1**(-1)*s2**(-1) - 4*h1**2*
     &    lambda1**2*s*s1**(-2) - 4*h2**2*lambda2**2*s*s2**(-2) - 32*
     &    ssz**2*lq2*pt2*c2b*s*sz**(-2) - 16*ssz**2*lq2*pt2*s*sz**(-2)
     &     - 64*ssz**2*lq2*cb**2*m1**2*mg**2*sz**(-2) + 64*ssz**2*lq2*
     &    cb**2*m1**2*msb2**2*sz**(-2) - 16*ssz**2*lq2*cb**2*mg**2*
     &    s**(-1)*t1**2*sz**(-2) - 16*ssz**2*lq2*cb**2*mg**2*s**(-1)*
     &    u1**2*sz**(-2) - 8*ssz**2*lq2*cb**2*mg**2*t1*sz**(-2) - 8*
     &    ssz**2*lq2*cb**2*mg**2*u1*sz**(-2) + 16*ssz**2*lq2*cb**2*
     &    msb2**2*s**(-1)*t1**2*sz**(-2) + 16*ssz**2*lq2*cb**2*msb2**2*
     &    s**(-1)*u1**2*sz**(-2) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * ( 8*ssz**2*lq2*
     &    cb**2*msb2**2*t1*sz**(-2) + 8*ssz**2*lq2*cb**2*msb2**2*u1*
     &    sz**(-2) + 64*ssz**2*lq2*m1**2*mg**2*sz**(-2) - 64*ssz**2*lq2
     &    *m1**2*msb2**2*sz**(-2) + 16*ssz**2*lq2*mg**2*s**(-1)*t1**2*
     &    sz**(-2) + 16*ssz**2*lq2*mg**2*s**(-1)*u1**2*sz**(-2) + 8*
     &    ssz**2*lq2*mg**2*t1*sz**(-2) + 8*ssz**2*lq2*mg**2*u1*sz**(-2)
     &     - 16*ssz**2*lq2*msb2**2*s**(-1)*t1**2*sz**(-2) - 16*ssz**2*
     &    lq2*msb2**2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*lq2*msb2**2*t1*
     &    sz**(-2) - 8*ssz**2*lq2*msb2**2*u1*sz**(-2) + 32*ssz**2*rq2*
     &    pt2*c2b*s*sz**(-2) - 16*ssz**2*rq2*pt2*s*sz**(-2) + 64*ssz**2
     &    *rq2*cb**2*m1**2*mg**2*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*
     &    msb2**2*sz**(-2) + 16*ssz**2*rq2*cb**2*mg**2*s**(-1)*t1**2*
     &    sz**(-2) + 16*ssz**2*rq2*cb**2*mg**2*s**(-1)*u1**2*sz**(-2)
     &     + 8*ssz**2*rq2*cb**2*mg**2*t1*sz**(-2) + 8*ssz**2*rq2*cb**2*
     &    mg**2*u1*sz**(-2) - 16*ssz**2*rq2*cb**2*msb2**2*s**(-1)*t1**2
     &    *sz**(-2) )
      MMs = MMs + SCB(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*ssz**2*rq2*
     &    cb**2*msb2**2*s**(-1)*u1**2*sz**(-2) - 8*ssz**2*rq2*cb**2*
     &    msb2**2*t1*sz**(-2) - 8*ssz**2*rq2*cb**2*msb2**2*u1*sz**(-2)
     &     - 32*ssp**2*pq2*pt2*s**(-1) + 64*ssp**2*pq2*m1**2*mg**2*
     &    s**(-2) - 64*ssp**2*pq2*m1**2*msb2**2*s**(-2) + 16*ssp**2*pq2
     &    *mg**2*s**(-3)*t1**2 + 16*ssp**2*pq2*mg**2*s**(-3)*u1**2 + 8*
     &    ssp**2*pq2*mg**2*s**(-2)*t1 + 8*ssp**2*pq2*mg**2*s**(-2)*u1
     &     - 16*ssp**2*pq2*msb2**2*s**(-3)*t1**2 - 16*ssp**2*pq2*
     &    msb2**2*s**(-3)*u1**2 - 8*ssp**2*pq2*msb2**2*s**(-2)*t1 - 8*
     &    ssp**2*pq2*msb2**2*s**(-2)*u1 )
      MMs = MMs + SCB(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*hss(1,1)*pq
     &    *hl*ssp*ct*cb*m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hl*
     &    ssp*ct*cb*mt*s**(-1)*u1*tx**(-1) - 32*hss(1,1)*pq*hr*ssp*st*
     &    sb*m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hr*ssp*st*sb*
     &    mt*s**(-1)*u1*tx**(-1) - 32*hss(1,1)*lq*hl*ssz*ct*cb*m1**2*mt
     &    *s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,1)*lq*hl*ssz*ct*cb*
     &    mt*u1*tx**(-1)*sz**(-1) - 32*hss(1,1)*rq*hr*ssz*st*sb*m1**2*
     &    mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,1)*rq*hr*ssz*st*sb
     &    *mt*u1*tx**(-1)*sz**(-1) + 8*hss(1,1)*hl*hr**2*ct*cb*m1**2*mt
     &    *s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hl*hr**2*ct*cb*mt*s*
     &    tx**(-2) + 16*hss(1,1)*hl*h1*lambda1*st*sb*sqrt2**(-1)*m1**2*
     &    s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(1,1)*hl*h1*lambda1*st*
     &    sb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(1,1)*hl*h2*
     &    lambda2*st*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1)
     &     + 16*hss(1,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*s*tx**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCB(2,2)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,1)*hl**2*
     &    hr*st*sb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hl**2*hr*
     &    st*sb*mt*s*tx**(-2) - 8*hss(1,1)*hl**3*ct*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) + 8*hss(1,1)*hl**3*ct*cb*mt*u1*tx**(-2) + 
     &    16*hss(1,1)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) + 16*hss(1,1)*hr*h1*lambda1*ct*cb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(1,1)*hr*h2*lambda2*
     &    ct*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(1,1)*hr*h2*lambda2*ct*cb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     - 8*hss(1,1)*hr**3*st*sb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(1,1)*hr**3*st*sb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(2,3)*Nc*Cf*Pi*alphas*prefac * ( 32*hss(2,1)*pq*hl
     &    *ssp*st*cb*m1**2*mt*t1**(-1)*tx**(-1) - 32*hss(2,1)*pq*hl*ssp
     &    *st*cb*mt*s**(-1)*u1*tx**(-1) - 32*hss(2,1)*pq*hr*ssp*ct*sb*
     &    m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(2,1)*pq*hr*ssp*ct*sb*mt*
     &    s**(-1)*u1*tx**(-1) + 32*hss(2,1)*lq*hl*ssz*st*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(2,1)*lq*hl*ssz*st*cb*mt*
     &    u1*tx**(-1)*sz**(-1) - 32*hss(2,1)*rq*hr*ssz*ct*sb*m1**2*mt*s
     &    *t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,1)*rq*hr*ssz*ct*sb*mt*
     &    u1*tx**(-1)*sz**(-1) - 8*hss(2,1)*hl*hr**2*st*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(2,1)*hl*hr**2*st*cb*mt*s*tx**(-2)
     &     + 16*hss(2,1)*hl*h1*lambda1*ct*sb*sqrt2**(-1)*m1**2*s*
     &    t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(2,1)*hl*h1*lambda1*ct*sb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(2,1)*hl*h2*lambda2*
     &    ct*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(2,1)*hl*h2*lambda2*ct*sb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     + 8*hss(2,1)*hl**2*hr*ct*sb*m1**2*mt*s*t1**(-1)*tx**(-2) )
      MMs = MMs + SCB(2,3)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,1)*hl**2*
     &    hr*ct*sb*mt*s*tx**(-2) + 8*hss(2,1)*hl**3*st*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(2,1)*hl**3*st*cb*mt*u1*tx**(-2) - 
     &    16*hss(2,1)*hr*h1*lambda1*st*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) - 16*hss(2,1)*hr*h1*lambda1*st*cb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) - 16*hss(2,1)*hr*h2*lambda2*
     &    st*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*
     &    hss(2,1)*hr*h2*lambda2*st*cb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     - 8*hss(2,1)*hr**3*ct*sb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(2,1)*hr**3*ct*sb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(2,4)*Nc*Cf*Pi*alphas*prefac * ( 32*hss(1,2)*pq*hl
     &    *ssp*ct*sb*m1**2*mt*t1**(-1)*tx**(-1) - 32*hss(1,2)*pq*hl*ssp
     &    *ct*sb*mt*s**(-1)*u1*tx**(-1) - 32*hss(1,2)*pq*hr*ssp*st*cb*
     &    m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(1,2)*pq*hr*ssp*st*cb*mt*
     &    s**(-1)*u1*tx**(-1) + 32*hss(1,2)*lq*hl*ssz*ct*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,2)*lq*hl*ssz*ct*sb*mt*
     &    u1*tx**(-1)*sz**(-1) - 32*hss(1,2)*rq*hr*ssz*st*cb*m1**2*mt*s
     &    *t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,2)*rq*hr*ssz*st*cb*mt*
     &    u1*tx**(-1)*sz**(-1) - 8*hss(1,2)*hl*hr**2*ct*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,2)*hl*hr**2*ct*sb*mt*s*tx**(-2)
     &     + 16*hss(1,2)*hl*h1*lambda1*st*cb*sqrt2**(-1)*m1**2*s*
     &    t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(1,2)*hl*h1*lambda1*st*cb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(1,2)*hl*h2*lambda2*
     &    st*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(1,2)*hl*h2*lambda2*st*cb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     + 8*hss(1,2)*hl**2*hr*st*cb*m1**2*mt*s*t1**(-1)*tx**(-2) )
      MMs = MMs + SCB(2,4)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hl**2*
     &    hr*st*cb*mt*s*tx**(-2) + 8*hss(1,2)*hl**3*ct*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,2)*hl**3*ct*sb*mt*u1*tx**(-2) - 
     &    16*hss(1,2)*hr*h1*lambda1*ct*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) - 16*hss(1,2)*hr*h1*lambda1*ct*sb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) - 16*hss(1,2)*hr*h2*lambda2*
     &    ct*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*
     &    hss(1,2)*hr*h2*lambda2*ct*sb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     - 8*hss(1,2)*hr**3*st*cb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(1,2)*hr**3*st*cb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(2,5)*Nc*Cf*Pi*alphas*prefac * (  - 32*hss(2,2)*pq
     &    *hl*ssp*st*sb*m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(2,2)*pq*hl*
     &    ssp*st*sb*mt*s**(-1)*u1*tx**(-1) - 32*hss(2,2)*pq*hr*ssp*ct*
     &    cb*m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(2,2)*pq*hr*ssp*ct*cb*
     &    mt*s**(-1)*u1*tx**(-1) - 32*hss(2,2)*lq*hl*ssz*st*sb*m1**2*mt
     &    *s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,2)*lq*hl*ssz*st*sb*
     &    mt*u1*tx**(-1)*sz**(-1) - 32*hss(2,2)*rq*hr*ssz*ct*cb*m1**2*
     &    mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,2)*rq*hr*ssz*ct*cb
     &    *mt*u1*tx**(-1)*sz**(-1) + 8*hss(2,2)*hl*hr**2*st*sb*m1**2*mt
     &    *s*t1**(-1)*tx**(-2) + 8*hss(2,2)*hl*hr**2*st*sb*mt*s*
     &    tx**(-2) + 16*hss(2,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*m1**2*
     &    s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(2,2)*hl*h1*lambda1*ct*
     &    cb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(2,2)*hl*h2*
     &    lambda2*ct*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1)
     &     + 16*hss(2,2)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*s*tx**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCB(2,5)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)*hl**2*
     &    hr*ct*cb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(2,2)*hl**2*hr*
     &    ct*cb*mt*s*tx**(-2) - 8*hss(2,2)*hl**3*st*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) + 8*hss(2,2)*hl**3*st*sb*mt*u1*tx**(-2) + 
     &    16*hss(2,2)*hr*h1*lambda1*st*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) + 16*hss(2,2)*hr*h1*lambda1*st*sb*
     &    sqrt2**(-1)*s*tx**(-1)*s1**(-1) + 16*hss(2,2)*hr*h2*lambda2*
     &    st*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(2,2)*hr*h2*lambda2*st*sb*sqrt2**(-1)*s*tx**(-1)*s2**(-1)
     &     - 8*hss(2,2)*hr**3*ct*cb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(2,2)*hr**3*ct*cb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * (  - 128*pq*lq*ssz*
     &    ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) + 128*pq*lq*ssz*ssp*
     &    cb**2*m1**2*msb1**2*s**(-1)*sz**(-1) - 64*pq*lq*ssz*ssp*cb**2
     &    *m1**2*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*mg**2*s**(-2)*t1*u1
     &    *sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*msb1**2*s**(-2)*t1*u1*
     &    sz**(-1) + 64*pq*lq*ssz*ssp*cb**2*s**(-1)*t1*u1*sz**(-1) + 
     &    128*pq*rq*ssz*ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) - 128*pq
     &    *rq*ssz*ssp*cb**2*m1**2*msb1**2*s**(-1)*sz**(-1) + 64*pq*rq*
     &    ssz*ssp*cb**2*m1**2*sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*mg**2*
     &    s**(-2)*t1*u1*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*msb1**2*
     &    s**(-2)*t1*u1*sz**(-1) - 64*pq*rq*ssz*ssp*cb**2*s**(-1)*t1*u1
     &    *sz**(-1) - 128*pq*rq*ssz*ssp*m1**2*mg**2*s**(-1)*sz**(-1) + 
     &    128*pq*rq*ssz*ssp*m1**2*msb1**2*s**(-1)*sz**(-1) - 64*pq*rq*
     &    ssz*ssp*m1**2*sz**(-1) + 128*pq*rq*ssz*ssp*mg**2*s**(-2)*t1*
     &    u1*sz**(-1) - 128*pq*rq*ssz*ssp*msb1**2*s**(-2)*t1*u1*
     &    sz**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 64*pq*rq*ssz*ssp*
     &    s**(-1)*t1*u1*sz**(-1) + 32*pq*hl*hr*ssp*sb*cb*mg*mt*s**(-1)*
     &    t1*tx**(-1) - 32*pq*hl*hr*ssp*sb*cb*mg*mt*s**(-1)*u1*tx**(-1)
     &     - 16*pq*hl**2*ssp*cb**2*m1**2*mg**2*s**(-1)*tx**(-1) + 16*pq
     &    *hl**2*ssp*cb**2*m1**2*msb1**2*s**(-1)*tx**(-1) - 8*pq*hl**2*
     &    ssp*cb**2*m1**2*tx**(-1) + 16*pq*hl**2*ssp*cb**2*mg**2*
     &    s**(-2)*t1*u1*tx**(-1) - 16*pq*hl**2*ssp*cb**2*msb1**2*
     &    s**(-2)*t1*u1*tx**(-1) + 8*pq*hl**2*ssp*cb**2*s**(-1)*t1*u1*
     &    tx**(-1) + 16*pq*hr**2*ssp*cb**2*m1**2*mg**2*s**(-1)*tx**(-1)
     &     - 16*pq*hr**2*ssp*cb**2*m1**2*msb1**2*s**(-1)*tx**(-1) + 8*
     &    pq*hr**2*ssp*cb**2*m1**2*tx**(-1) - 16*pq*hr**2*ssp*cb**2*
     &    mg**2*s**(-2)*t1*u1*tx**(-1) + 16*pq*hr**2*ssp*cb**2*msb1**2*
     &    s**(-2)*t1*u1*tx**(-1) - 8*pq*hr**2*ssp*cb**2*s**(-1)*t1*u1*
     &    tx**(-1) - 16*pq*hr**2*ssp*m1**2*mg**2*s**(-1)*tx**(-1) + 16*
     &    pq*hr**2*ssp*m1**2*msb1**2*s**(-1)*tx**(-1) - 8*pq*hr**2*ssp*
     &    m1**2*tx**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hr**2*ssp*
     &    mg**2*s**(-2)*t1*u1*tx**(-1) - 16*pq*hr**2*ssp*msb1**2*
     &    s**(-2)*t1*u1*tx**(-1) + 8*pq*hr**2*ssp*s**(-1)*t1*u1*
     &    tx**(-1) + 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*s**(-1)*
     &    t1*s1**(-1) - 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*
     &    s**(-1)*u1*s1**(-1) + 64*pq*h2*ssp*lambda2*sb*cb*sqrt2**(-1)*
     &    mg*s**(-1)*t1*s2**(-1) - 64*pq*h2*ssp*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s**(-1)*u1*s2**(-1) - 128*lq*rq*ssz**2*cb**2*
     &    m1**2*mg**2*sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*msb1**2*
     &    sz**(-2) - 64*lq*rq*ssz**2*cb**2*m1**2*s*sz**(-2) + 128*lq*rq
     &    *ssz**2*cb**2*mg**2*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*ssz**2
     &    *cb**2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*lq*rq*ssz**2*cb**2
     &    *t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*mg**2*sz**(-2)
     &     - 128*lq*rq*ssz**2*cb**4*m1**2*msb1**2*sz**(-2) + 64*lq*rq*
     &    ssz**2*cb**4*m1**2*s*sz**(-2) - 128*lq*rq*ssz**2*cb**4*mg**2*
     &    s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 128*lq*rq*ssz**2*
     &    cb**4*msb1**2*s**(-1)*t1*u1*sz**(-2) - 64*lq*rq*ssz**2*cb**4*
     &    t1*u1*sz**(-2) + 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*t1*tx**(-1)*
     &    sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*cb**4*m1**2*msb1**2*tx**(-1)*
     &    sz**(-1) - 8*lq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*sz**(-1) + 
     &    16*lq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*lq*hl**2*ssz*cb**4*msb1**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 8*lq*hl**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    lq*hr**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*sz**(-1) + 16*lq*
     &    hr**2*ssz*cb**2*m1**2*msb1**2*tx**(-1)*sz**(-1) - 8*lq*hr**2*
     &    ssz*cb**2*m1**2*s*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*
     &    mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*lq*hr**2*ssz*
     &    cb**2*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**4*m1**2*mg**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*
     &    m1**2*msb1**2*tx**(-1)*sz**(-1) + 8*lq*hr**2*ssz*cb**4*m1**2*
     &    s*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*mg**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*msb1**2*s**(-1)*
     &    t1*u1*tx**(-1)*sz**(-1) - 8*lq*hr**2*ssz*cb**4*t1*u1*tx**(-1)
     &    *sz**(-1) + 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*t1*
     &    sz**(-1)*s1**(-1) - 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*u1*sz**(-1)*s1**(-1) + 64*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb
     &    *cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) + 32*rq*hl*hr*ssz*
     &    sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb*mg*
     &    mt*u1*tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*t1*
     &    tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*u1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hl**2*ssz*
     &    cb**2*m1**2*msb1**2*tx**(-1)*sz**(-1) - 8*rq*hl**2*ssz*cb**2*
     &    m1**2*s*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*mg**2*
     &    s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*
     &    msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*rq*hl**2*ssz*
     &    cb**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*
     &    mg**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*msb1**2
     &    *tx**(-1)*sz**(-1) + 8*rq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) - 16*rq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hl**2*ssz*cb**4*msb1**2*s**(-1)*t1*u1*
     &    tx**(-1)*sz**(-1) - 8*rq*hl**2*ssz*cb**4*t1*u1*tx**(-1)*
     &    sz**(-1) + 32*rq*hr**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) - 32*rq*hr**2*ssz*cb**2*m1**2*msb1**2*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**2*m1**2*s*tx**(-1)*sz**(-1) - 
     &    32*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     + 32*rq*hr**2*ssz*cb**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*
     &    ssz*cb**2*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*
     &    m1**2*mg**2*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*m1**2*
     &    msb1**2*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*cb**4*m1**2*s*
     &    tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*cb**4*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*m1**2*mg**2*tx**(-1)*sz**(-1) + 16
     &    *rq*hr**2*ssz*m1**2*msb1**2*tx**(-1)*sz**(-1) - 8*rq*hr**2*
     &    ssz*m1**2*s*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*mg**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*t1*u1*tx**(-1)*
     &    sz**(-1) + 64*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*t1*
     &    sz**(-1)*s1**(-1) - 64*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg
     &    *u1*sz**(-1)*s1**(-1) - 64*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * ( 64*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) + 64*rq*
     &    h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) - 64
     &    *rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1)
     &     - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*
     &    s2**(-1) + 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*u1*
     &    sz**(-1)*s2**(-1) - 64*ssz**2*lq2*cb**4*m1**2*mg**2*sz**(-2)
     &     + 64*ssz**2*lq2*cb**4*m1**2*msb1**2*sz**(-2) - 32*ssz**2*lq2
     &    *cb**4*m1**2*s*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*s**(-1)*
     &    t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**4*msb1**2*s**(-1)*t1*u1*
     &    sz**(-2) + 32*ssz**2*lq2*cb**4*t1*u1*sz**(-2) + 128*ssz**2*
     &    rq2*cb**2*m1**2*mg**2*sz**(-2) - 128*ssz**2*rq2*cb**2*m1**2*
     &    msb1**2*sz**(-2) + 64*ssz**2*rq2*cb**2*m1**2*s*sz**(-2) - 128
     &    *ssz**2*rq2*cb**2*mg**2*s**(-1)*t1*u1*sz**(-2) + 128*ssz**2*
     &    rq2*cb**2*msb1**2*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*rq2*
     &    cb**2*t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,4)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*rq2*
     &    cb**4*m1**2*mg**2*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*
     &    msb1**2*sz**(-2) - 32*ssz**2*rq2*cb**4*m1**2*s*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*mg**2*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*rq2
     &    *cb**4*msb1**2*s**(-1)*t1*u1*sz**(-2) + 32*ssz**2*rq2*cb**4*
     &    t1*u1*sz**(-2) - 64*ssz**2*rq2*m1**2*mg**2*sz**(-2) + 64*
     &    ssz**2*rq2*m1**2*msb1**2*sz**(-2) - 32*ssz**2*rq2*m1**2*s*
     &    sz**(-2) + 64*ssz**2*rq2*mg**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*rq2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 32*ssz**2*rq2*t1*
     &    u1*sz**(-2) - 64*ssp**2*pq2*m1**2*mg**2*s**(-2) + 64*ssp**2*
     &    pq2*m1**2*msb1**2*s**(-2) - 32*ssp**2*pq2*m1**2*s**(-1) + 64*
     &    ssp**2*pq2*mg**2*s**(-3)*t1*u1 - 64*ssp**2*pq2*msb1**2*
     &    s**(-3)*t1*u1 + 32*ssp**2*pq2*s**(-2)*t1*u1 )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * ( 256*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**2*sz**(-2) - 128*lq*rq*ssz**2*cb**2*m1**2*
     &    msb1**2*sz**(-2) - 128*lq*rq*ssz**2*cb**2*m1**2*msb2**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*s*sz**(-2) - 256*lq*
     &    rq*ssz**2*cb**2*mg**2*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**2*msb2**2*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*
     &    ssz**2*cb**2*t1*u1*sz**(-2) - 256*lq*rq*ssz**2*cb**4*m1**2*
     &    mg**2*sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*msb1**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*msb2**2*sz**(-2) - 
     &    128*lq*rq*ssz**2*cb**4*m1**2*s*sz**(-2) + 256*lq*rq*ssz**2*
     &    cb**4*mg**2*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*
     &    msb1**2*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*
     &    msb2**2*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*t1*u1
     &    *sz**(-2) + 32*lq*hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1)
     &     - 32*lq*hl*hr*ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * (  - 64*lq*hl*hr*
     &    ssz*sb*cb**3*mg*mt*t1*tx**(-1)*sz**(-1) + 64*lq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*u1*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**2*
     &    m1**2*mg**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2*
     &    msb1**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2*s*
     &    tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*cb**2*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*msb2**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*t1*u1*
     &    tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*cb**4*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*msb1**2*
     &    tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*msb2**2*
     &    tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) - 32*lq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hl**2*ssz*
     &    cb**4*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*
     &    ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*
     &    hl**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) + 32*lq*hr**2*ssz*
     &    cb**2*m1**2*mg**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*
     &    m1**2*msb1**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2
     &    *msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*s*
     &    tx**(-1)*sz**(-1) - 32*lq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*msb2**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*t1*u1*
     &    tx**(-1)*sz**(-1) - 32*lq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*msb1**2*
     &    tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*msb2**2*
     &    tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * ( 32*lq*hr**2*ssz*
     &    cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz
     &    *cb**4*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*
     &    ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*
     &    hr**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) + 64*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) - 64*lq*h1*
     &    ssz*lambda1*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) - 128*
     &    lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*
     &    s1**(-1) + 128*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*u1*
     &    sz**(-1)*s1**(-1) + 64*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg
     &    *t1*sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 128*lq*h2*ssz*lambda2*
     &    sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) + 128*lq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 32*rq*
     &    hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*
     &    sb*cb*mg*mt*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * ( 64*rq*hl*hr*ssz*
     &    sb*cb**3*mg*mt*t1*tx**(-1)*sz**(-1) - 64*rq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*u1*tx**(-1)*sz**(-1) + 32*rq*hl**2*ssz*cb**2*
     &    m1**2*mg**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*
     &    msb1**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*
     &    msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*m1**2*s*
     &    tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*cb**2*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*msb2**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*t1*u1*
     &    tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*cb**4*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*msb1**2*
     &    tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*msb2**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) + 32*rq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hl**2*
     &    ssz*cb**4*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*
     &    hl**2*ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    rq*hl**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*
     &    cb**2*m1**2*mg**2*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**2*
     &    m1**2*msb1**2*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**2*m1**2
     &    *msb2**2*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**2*m1**2*s*
     &    tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**2*msb1**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**2*msb2**2*s**(-1)
     &    *t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**2*t1*u1*
     &    tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*m1**2*msb1**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*m1**2*msb2**2*
     &    tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*hr**2*
     &    ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2
     &    *ssz*cb**4*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) - 64*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) + 64*rq*h1*
     &    ssz*lambda1*sb*cb*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) + 128*
     &    rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*
     &    s1**(-1) - 128*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*u1*
     &    sz**(-1)*s1**(-1) - 64*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg
     &    *t1*sz**(-1)*s2**(-1) + 64*rq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) + 128*rq*h2*ssz*lambda2*
     &    sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*s2**(-1) - 128*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 128*
     &    ssz**2*lq2*cb**2*m1**2*mg**2*sz**(-2) + 64*ssz**2*lq2*cb**2*
     &    m1**2*msb1**2*sz**(-2) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * ( 64*ssz**2*lq2*
     &    cb**2*m1**2*msb2**2*sz**(-2) - 64*ssz**2*lq2*cb**2*m1**2*s*
     &    sz**(-2) + 128*ssz**2*lq2*cb**2*mg**2*s**(-1)*t1*u1*sz**(-2)
     &     - 64*ssz**2*lq2*cb**2*msb1**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*lq2*cb**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*
     &    lq2*cb**2*t1*u1*sz**(-2) + 128*ssz**2*lq2*cb**4*m1**2*mg**2*
     &    sz**(-2) - 64*ssz**2*lq2*cb**4*m1**2*msb1**2*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*m1**2*msb2**2*sz**(-2) + 64*ssz**2*lq2*cb**4
     &    *m1**2*s*sz**(-2) - 128*ssz**2*lq2*cb**4*mg**2*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*lq2*cb**4*msb1**2*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*lq2*cb**4*msb2**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*t1*u1*sz**(-2) - 128*ssz**2*rq2*cb**2*m1**2*
     &    mg**2*sz**(-2) + 64*ssz**2*rq2*cb**2*m1**2*msb1**2*sz**(-2)
     &     + 64*ssz**2*rq2*cb**2*m1**2*msb2**2*sz**(-2) - 64*ssz**2*rq2
     &    *cb**2*m1**2*s*sz**(-2) + 128*ssz**2*rq2*cb**2*mg**2*s**(-1)*
     &    t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,5)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*rq2*
     &    cb**2*msb1**2*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*rq2*cb**2*
     &    msb2**2*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*rq2*cb**2*t1*u1*
     &    sz**(-2) + 128*ssz**2*rq2*cb**4*m1**2*mg**2*sz**(-2) - 64*
     &    ssz**2*rq2*cb**4*m1**2*msb1**2*sz**(-2) - 64*ssz**2*rq2*cb**4
     &    *m1**2*msb2**2*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*s*
     &    sz**(-2) - 128*ssz**2*rq2*cb**4*mg**2*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*rq2*cb**4*msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*msb2**2*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*
     &    rq2*cb**4*t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 128*pq*lq*ssz*ssp
     &    *cb**2*m1**2*mg**2*s**(-1)*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2
     &    *m1**2*msb2**2*s**(-1)*sz**(-1) + 64*pq*lq*ssz*ssp*cb**2*
     &    m1**2*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*mg**2*s**(-2)*t1*u1*
     &    sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*msb2**2*s**(-2)*t1*u1*
     &    sz**(-1) - 64*pq*lq*ssz*ssp*cb**2*s**(-1)*t1*u1*sz**(-1) - 
     &    128*pq*lq*ssz*ssp*m1**2*mg**2*s**(-1)*sz**(-1) + 128*pq*lq*
     &    ssz*ssp*m1**2*msb2**2*s**(-1)*sz**(-1) - 64*pq*lq*ssz*ssp*
     &    m1**2*sz**(-1) + 128*pq*lq*ssz*ssp*mg**2*s**(-2)*t1*u1*
     &    sz**(-1) - 128*pq*lq*ssz*ssp*msb2**2*s**(-2)*t1*u1*sz**(-1)
     &     + 64*pq*lq*ssz*ssp*s**(-1)*t1*u1*sz**(-1) - 128*pq*rq*ssz*
     &    ssp*cb**2*m1**2*mg**2*s**(-1)*sz**(-1) + 128*pq*rq*ssz*ssp*
     &    cb**2*m1**2*msb2**2*s**(-1)*sz**(-1) - 64*pq*rq*ssz*ssp*cb**2
     &    *m1**2*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*mg**2*s**(-2)*t1*u1
     &    *sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*msb2**2*s**(-2)*t1*u1*
     &    sz**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 64*pq*rq*ssz*ssp*
     &    cb**2*s**(-1)*t1*u1*sz**(-1) - 32*pq*hl*hr*ssp*sb*cb*mg*mt*
     &    s**(-1)*t1*tx**(-1) + 32*pq*hl*hr*ssp*sb*cb*mg*mt*s**(-1)*u1*
     &    tx**(-1) + 16*pq*hl**2*ssp*cb**2*m1**2*mg**2*s**(-1)*tx**(-1)
     &     - 16*pq*hl**2*ssp*cb**2*m1**2*msb2**2*s**(-1)*tx**(-1) + 8*
     &    pq*hl**2*ssp*cb**2*m1**2*tx**(-1) - 16*pq*hl**2*ssp*cb**2*
     &    mg**2*s**(-2)*t1*u1*tx**(-1) + 16*pq*hl**2*ssp*cb**2*msb2**2*
     &    s**(-2)*t1*u1*tx**(-1) - 8*pq*hl**2*ssp*cb**2*s**(-1)*t1*u1*
     &    tx**(-1) - 16*pq*hl**2*ssp*m1**2*mg**2*s**(-1)*tx**(-1) + 16*
     &    pq*hl**2*ssp*m1**2*msb2**2*s**(-1)*tx**(-1) - 8*pq*hl**2*ssp*
     &    m1**2*tx**(-1) + 16*pq*hl**2*ssp*mg**2*s**(-2)*t1*u1*tx**(-1)
     &     - 16*pq*hl**2*ssp*msb2**2*s**(-2)*t1*u1*tx**(-1) + 8*pq*
     &    hl**2*ssp*s**(-1)*t1*u1*tx**(-1) - 16*pq*hr**2*ssp*cb**2*
     &    m1**2*mg**2*s**(-1)*tx**(-1) + 16*pq*hr**2*ssp*cb**2*m1**2*
     &    msb2**2*s**(-1)*tx**(-1) - 8*pq*hr**2*ssp*cb**2*m1**2*
     &    tx**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hr**2*ssp*
     &    cb**2*mg**2*s**(-2)*t1*u1*tx**(-1) - 16*pq*hr**2*ssp*cb**2*
     &    msb2**2*s**(-2)*t1*u1*tx**(-1) + 8*pq*hr**2*ssp*cb**2*s**(-1)
     &    *t1*u1*tx**(-1) - 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*
     &    s**(-1)*t1*s1**(-1) + 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*
     &    mg*s**(-1)*u1*s1**(-1) - 64*pq*h2*ssp*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s**(-1)*t1*s2**(-1) + 64*pq*h2*ssp*lambda2*sb*
     &    cb*sqrt2**(-1)*mg*s**(-1)*u1*s2**(-1) - 128*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**2*sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*
     &    msb2**2*sz**(-2) - 64*lq*rq*ssz**2*cb**2*m1**2*s*sz**(-2) + 
     &    128*lq*rq*ssz**2*cb**2*mg**2*s**(-1)*t1*u1*sz**(-2) - 128*lq*
     &    rq*ssz**2*cb**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 64*lq*rq*
     &    ssz**2*cb**2*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*
     &    mg**2*sz**(-2) - 128*lq*rq*ssz**2*cb**4*m1**2*msb2**2*
     &    sz**(-2) + 64*lq*rq*ssz**2*cb**4*m1**2*s*sz**(-2) - 128*lq*rq
     &    *ssz**2*cb**4*mg**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 128*lq*rq*ssz**2*
     &    cb**4*msb2**2*s**(-1)*t1*u1*sz**(-2) - 64*lq*rq*ssz**2*cb**4*
     &    t1*u1*sz**(-2) - 32*lq*hl*hr*ssz*sb*cb*mg*mt*t1*tx**(-1)*
     &    sz**(-1) + 32*lq*hl*hr*ssz*sb*cb*mg*mt*u1*tx**(-1)*sz**(-1)
     &     + 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*t1*tx**(-1)*sz**(-1) - 32*
     &    lq*hl*hr*ssz*sb*cb**3*mg*mt*u1*tx**(-1)*sz**(-1) + 32*lq*
     &    hl**2*ssz*cb**2*m1**2*mg**2*tx**(-1)*sz**(-1) - 32*lq*hl**2*
     &    ssz*cb**2*m1**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*
     &    cb**2*m1**2*s*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**2*mg**2
     &    *s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*cb**2*
     &    msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*
     &    cb**2*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*
     &    mg**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**4*m1**2*msb2**2
     &    *tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * (  - 16*lq*hl**2*
     &    ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*lq*
     &    hl**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*
     &    m1**2*mg**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*m1**2*msb2**2
     &    *tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*m1**2*s*tx**(-1)*sz**(-1)
     &     + 16*lq*hl**2*ssz*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16
     &    *lq*hl**2*ssz*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*lq*
     &    hl**2*ssz*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*
     &    m1**2*mg**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 8*lq*hr**2*ssz*cb**2*m1**2*s*
     &    tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*mg**2*s**(-1)*t1*u1
     &    *tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*msb2**2*s**(-1)*t1
     &    *u1*tx**(-1)*sz**(-1) + 8*lq*hr**2*ssz*cb**2*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*msb2**2*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 8*lq*hr**2*ssz*
     &    cb**4*m1**2*s*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*mg**2
     &    *s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*
     &    msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 8*lq*hr**2*ssz*
     &    cb**4*t1*u1*tx**(-1)*sz**(-1) - 64*lq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1*sb
     &    *cb*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) - 64*lq*
     &    h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1)
     &     - 64*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*t1*sz**(-1)*
     &    s2**(-1) + 64*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*u1*
     &    sz**(-1)*s2**(-1) + 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)
     &    *mg*t1*sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*u1*sz**(-1)*s2**(-1) - 32*rq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*t1*tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*
     &    mg*mt*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hl**2*
     &    ssz*cb**2*m1**2*mg**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**2*m1**2*msb2**2*tx**(-1)*sz**(-1) - 8*rq*hl**2*ssz*cb**2*
     &    m1**2*s*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*mg**2*
     &    s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*
     &    msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*rq*hl**2*ssz*
     &    cb**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*
     &    mg**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*msb2**2
     &    *tx**(-1)*sz**(-1) + 8*rq*hl**2*ssz*cb**4*m1**2*s*tx**(-1)*
     &    sz**(-1) - 16*rq*hl**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hl**2*ssz*cb**4*msb2**2*s**(-1)*t1*u1*
     &    tx**(-1)*sz**(-1) - 8*rq*hl**2*ssz*cb**4*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*rq*hr**2*ssz*cb**4*m1**2*mg**2*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**4*m1**2*msb2**2*tx**(-1)*
     &    sz**(-1) - 8*rq*hr**2*ssz*cb**4*m1**2*s*tx**(-1)*sz**(-1) + 
     &    16*rq*hr**2*ssz*cb**4*mg**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*
     &    ssz*cb**4*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 8*rq*
     &    hr**2*ssz*cb**4*t1*u1*tx**(-1)*sz**(-1) - 64*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*s1**(-1) + 64*rq*
     &    h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*u1*sz**(-1)*s1**(-1)
     &     - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*t1*sz**(-1)*
     &    s2**(-1) + 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*u1*
     &    sz**(-1)*s2**(-1) + 128*ssz**2*lq2*cb**2*m1**2*mg**2*sz**(-2)
     &     - 128*ssz**2*lq2*cb**2*m1**2*msb2**2*sz**(-2) + 64*ssz**2*
     &    lq2*cb**2*m1**2*s*sz**(-2) - 128*ssz**2*lq2*cb**2*mg**2*
     &    s**(-1)*t1*u1*sz**(-2) + 128*ssz**2*lq2*cb**2*msb2**2*s**(-1)
     &    *t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**2*t1*u1*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*m1**2*mg**2*sz**(-2) + 64*ssz**2*lq2*cb**4*
     &    m1**2*msb2**2*sz**(-2) - 32*ssz**2*lq2*cb**4*m1**2*s*sz**(-2)
     &     + 64*ssz**2*lq2*cb**4*mg**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*msb2**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCB(4,6)*Nc*Cf*Pi*alphas*prefac * ( 32*ssz**2*lq2*
     &    cb**4*t1*u1*sz**(-2) - 64*ssz**2*lq2*m1**2*mg**2*sz**(-2) + 
     &    64*ssz**2*lq2*m1**2*msb2**2*sz**(-2) - 32*ssz**2*lq2*m1**2*s*
     &    sz**(-2) + 64*ssz**2*lq2*mg**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*lq2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 32*ssz**2*lq2*t1*
     &    u1*sz**(-2) - 64*ssz**2*rq2*cb**4*m1**2*mg**2*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*m1**2*msb2**2*sz**(-2) - 32*ssz**2*rq2*cb**4
     &    *m1**2*s*sz**(-2) + 64*ssz**2*rq2*cb**4*mg**2*s**(-1)*t1*u1*
     &    sz**(-2) - 64*ssz**2*rq2*cb**4*msb2**2*s**(-1)*t1*u1*sz**(-2)
     &     + 32*ssz**2*rq2*cb**4*t1*u1*sz**(-2) - 64*ssp**2*pq2*m1**2*
     &    mg**2*s**(-2) + 64*ssp**2*pq2*m1**2*msb2**2*s**(-2) - 32*
     &    ssp**2*pq2*m1**2*s**(-1) + 64*ssp**2*pq2*mg**2*s**(-3)*t1*u1
     &     - 64*ssp**2*pq2*msb2**2*s**(-3)*t1*u1 + 32*ssp**2*pq2*
     &    s**(-2)*t1*u1 )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hl**2*ssp*
     &    pt2*s2t*mg*mt*tx**(-2) - 4*pq*hl**2*ssp*pt2*c2t*m1**2*mg**2*
     &    tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*c2t*m1**2*mst1**2*
     &    tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*c2t*m1**2*tx**(-2) + 4*
     &    pq*hl**2*ssp*pt2*c2t*mg**2*mt**2*tx**(-2)*t**(-1) - 4*pq*
     &    hl**2*ssp*pt2*c2t*mg**2*t1*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*
     &    pt2*c2t*mt**2*mst1**2*tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*
     &    c2t*mt**2*tx**(-2) + 4*pq*hl**2*ssp*pt2*c2t*mst1**2*t1*
     &    tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*c2t*t1*tx**(-2) - 4*pq*
     &    hl**2*ssp*pt2*m1**2*mg**2*tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*
     &    pt2*m1**2*mst1**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*m1**2
     &    *tx**(-2) - 4*pq*hl**2*ssp*pt2*mg**2*mt**2*tx**(-2)*t**(-1)
     &     - 4*pq*hl**2*ssp*pt2*mg**2*t1*tx**(-2)*t**(-1) + 4*pq*hl**2*
     &    ssp*pt2*mt**2*mst1**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*
     &    mt**2*tx**(-2) + 4*pq*hl**2*ssp*pt2*mst1**2*t1*tx**(-2)*
     &    t**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hl**2*ssp
     &    *pt2*t1*tx**(-2) + 16*pq*hr**2*ssp*pt2*s2t*mg*mt*tx**(-2) + 4
     &    *pq*hr**2*ssp*pt2*c2t*m1**2*mg**2*tx**(-2)*t**(-1) - 4*pq*
     &    hr**2*ssp*pt2*c2t*m1**2*mst1**2*tx**(-2)*t**(-1) + 4*pq*hr**2
     &    *ssp*pt2*c2t*m1**2*tx**(-2) - 4*pq*hr**2*ssp*pt2*c2t*mg**2*
     &    mt**2*tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*c2t*mg**2*t1*
     &    tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*c2t*mt**2*mst1**2*
     &    tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*c2t*mt**2*tx**(-2) - 4*
     &    pq*hr**2*ssp*pt2*c2t*mst1**2*t1*tx**(-2)*t**(-1) + 4*pq*hr**2
     &    *ssp*pt2*c2t*t1*tx**(-2) - 4*pq*hr**2*ssp*pt2*m1**2*mg**2*
     &    tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*m1**2*mst1**2*tx**(-2)*
     &    t**(-1) - 4*pq*hr**2*ssp*pt2*m1**2*tx**(-2) - 4*pq*hr**2*ssp*
     &    pt2*mg**2*mt**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*mg**2*
     &    t1*tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*mt**2*mst1**2*
     &    tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*mt**2*tx**(-2) + 4*pq*
     &    hr**2*ssp*pt2*mst1**2*t1*tx**(-2)*t**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hr**2*ssp
     &    *pt2*t1*tx**(-2) + 16*lq*hl**2*ssz*pt2*s2t*mg*mt*s*tx**(-2)*
     &    sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*m1**2*mg**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*m1**2*mst1**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*m1**2*s*
     &    tx**(-2)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*mg**2*mt**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*mg**2*s*t1
     &    *tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*mt**2*
     &    mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*
     &    mt**2*s*tx**(-2)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*mst1**2*s*
     &    t1*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*s*t1*
     &    tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*m1**2*mg**2*s*tx**(-2)
     &    *t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*m1**2*mst1**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*m1**2*s*
     &    tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*mg**2*mt**2*s*tx**(-2)
     &    *t**(-1)*sz**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 4*lq*hl**2*ssz
     &    *pt2*mg**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*
     &    pt2*mt**2*mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*
     &    ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) + 4*lq*hl**2*ssz*pt2*
     &    mst1**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*s
     &    *t1*tx**(-2)*sz**(-1) + 16*rq*hr**2*ssz*pt2*s2t*mg*mt*s*
     &    tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*m1**2*mg**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*m1**2*
     &    mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*
     &    m1**2*s*tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*mg**2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*
     &    mg**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t
     &    *mt**2*mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*
     &    pt2*c2t*mt**2*s*tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*
     &    mst1**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*
     &    c2t*s*t1*tx**(-2)*sz**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 4*rq*hr**2*ssz
     &    *pt2*m1**2*mg**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz
     &    *pt2*m1**2*mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*
     &    ssz*pt2*m1**2*s*tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mg**2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mg**2*
     &    s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*mt**2*
     &    mst1**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*
     &    mt**2*s*tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*mst1**2*s*t1*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*s*t1*tx**(-2)*
     &    sz**(-1) + 8*hl*hr*h1*lambda1*s2t*sqrt2**(-1)*m1**2*mg*s*
     &    tx**(-2)*s1**(-1) + 8*hl*hr*h1*lambda1*s2t*sqrt2**(-1)*mg*
     &    mt**2*s*tx**(-2)*s1**(-1) + 8*hl*hr*h1*lambda1*s2t*
     &    sqrt2**(-1)*mg*s*t1*tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mg**2*mt*s*tx**(-2)*t**(-1)*s1**(-1) + 8*hl
     &    *hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*mst1**2*s*tx**(-2)*
     &    t**(-1)*s1**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) - 8*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mg**2*mt*s*t1*tx**(-2)*t**(-1)*s1**(-1)
     &     + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mst1**2*s*t1*tx**(-2)*
     &    t**(-1)*s1**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*
     &    tx**(-2)*s1**(-1) + 8*hl*hr*h2*lambda2*s2t*sqrt2**(-1)*m1**2*
     &    mg*s*tx**(-2)*s2**(-1) + 8*hl*hr*h2*lambda2*s2t*sqrt2**(-1)*
     &    mg*mt**2*s*tx**(-2)*s2**(-1) + 8*hl*hr*h2*lambda2*s2t*
     &    sqrt2**(-1)*mg*s*t1*tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mg**2*mt*s*tx**(-2)*t**(-1)*s2**(-1) + 8*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*mst1**2*s*tx**(-2)*
     &    t**(-1)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s*
     &    tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mg**2*mt*s
     &    *t1*tx**(-2)*t**(-1)*s2**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mst1**2*s*t1*tx**(-2)*t**(-1)*s2**(-1) - 8*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * ( 4*hl**2*hr**2*s2t
     &    *m1**2*mg*mt*s*tx**(-3) + 4*hl**2*hr**2*s2t*mg*mt*s*t1*
     &    tx**(-3) + 4*hl**2*hr**2*s2t*mg*mt**3*s*tx**(-3) - 4*hl**2*
     &    hr**2*m1**2*mg**2*mt**2*s*tx**(-3)*t**(-1) + 4*hl**2*hr**2*
     &    m1**2*mt**2*mst1**2*s*tx**(-3)*t**(-1) - 4*hl**2*hr**2*m1**2*
     &    mt**2*s*tx**(-3) - 4*hl**2*hr**2*mg**2*mt**2*s*t1*tx**(-3)*
     &    t**(-1) + 4*hl**2*hr**2*mt**2*mst1**2*s*t1*tx**(-3)*t**(-1)
     &     - 4*hl**2*hr**2*mt**2*s*t1*tx**(-3) + 4*hl**4*pt2*s2t*mg*mt*
     &    s*tx**(-3) - hl**4*pt2*c2t*m1**2*mg**2*s*tx**(-3)*t**(-1) + 
     &    hl**4*pt2*c2t*m1**2*mst1**2*s*tx**(-3)*t**(-1) - hl**4*pt2*
     &    c2t*m1**2*s*tx**(-3) + hl**4*pt2*c2t*mg**2*mt**2*s*tx**(-3)*
     &    t**(-1) - hl**4*pt2*c2t*mg**2*s*t1*tx**(-3)*t**(-1) - hl**4*
     &    pt2*c2t*mt**2*mst1**2*s*tx**(-3)*t**(-1) + hl**4*pt2*c2t*
     &    mt**2*s*tx**(-3) + hl**4*pt2*c2t*mst1**2*s*t1*tx**(-3)*
     &    t**(-1) - hl**4*pt2*c2t*s*t1*tx**(-3) - hl**4*pt2*m1**2*mg**2
     &    *s*tx**(-3)*t**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * ( hl**4*pt2*m1**2*
     &    mst1**2*s*tx**(-3)*t**(-1) - hl**4*pt2*m1**2*s*tx**(-3) - 
     &    hl**4*pt2*mg**2*mt**2*s*tx**(-3)*t**(-1) - hl**4*pt2*mg**2*s*
     &    t1*tx**(-3)*t**(-1) + hl**4*pt2*mt**2*mst1**2*s*tx**(-3)*
     &    t**(-1) - hl**4*pt2*mt**2*s*tx**(-3) + hl**4*pt2*mst1**2*s*t1
     &    *tx**(-3)*t**(-1) - hl**4*pt2*s*t1*tx**(-3) + 4*hr**4*pt2*s2t
     &    *mg*mt*s*tx**(-3) + hr**4*pt2*c2t*m1**2*mg**2*s*tx**(-3)*
     &    t**(-1) - hr**4*pt2*c2t*m1**2*mst1**2*s*tx**(-3)*t**(-1) + 
     &    hr**4*pt2*c2t*m1**2*s*tx**(-3) - hr**4*pt2*c2t*mg**2*mt**2*s*
     &    tx**(-3)*t**(-1) + hr**4*pt2*c2t*mg**2*s*t1*tx**(-3)*t**(-1)
     &     + hr**4*pt2*c2t*mt**2*mst1**2*s*tx**(-3)*t**(-1) - hr**4*pt2
     &    *c2t*mt**2*s*tx**(-3) - hr**4*pt2*c2t*mst1**2*s*t1*tx**(-3)*
     &    t**(-1) + hr**4*pt2*c2t*s*t1*tx**(-3) - hr**4*pt2*m1**2*mg**2
     &    *s*tx**(-3)*t**(-1) + hr**4*pt2*m1**2*mst1**2*s*tx**(-3)*
     &    t**(-1) - hr**4*pt2*m1**2*s*tx**(-3) - hr**4*pt2*mg**2*mt**2*
     &    s*tx**(-3)*t**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - hr**4*pt2*
     &    mg**2*s*t1*tx**(-3)*t**(-1) + hr**4*pt2*mt**2*mst1**2*s*
     &    tx**(-3)*t**(-1) - hr**4*pt2*mt**2*s*tx**(-3) + hr**4*pt2*
     &    mst1**2*s*t1*tx**(-3)*t**(-1) - hr**4*pt2*s*t1*tx**(-3) + 32*
     &    hss(1,1)*pq*hl*ssp*ct*cb*m1**2*mt*t1**(-1)*tx**(-1) - 32*hss(
     &    1,1)*pq*hl*ssp*ct*cb*mt*s**(-1)*u1*tx**(-1) + 32*hss(1,1)*pq*
     &    hr*ssp*st*sb*m1**2*mt*t1**(-1)*tx**(-1) - 32*hss(1,1)*pq*hr*
     &    ssp*st*sb*mt*s**(-1)*u1*tx**(-1) + 32*hss(1,1)*lq*hl*ssz*ct*
     &    cb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,1)*lq*hl*
     &    ssz*ct*cb*mt*u1*tx**(-1)*sz**(-1) + 32*hss(1,1)*rq*hr*ssz*st*
     &    sb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,1)*rq*hr*
     &    ssz*st*sb*mt*u1*tx**(-1)*sz**(-1) - 8*hss(1,1)*hl*hr**2*ct*cb
     &    *m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(1,1)*hl*hr**2*ct*cb*mt*
     &    s*tx**(-2) - 16*hss(1,1)*hl*h1*lambda1*st*sb*sqrt2**(-1)*
     &    m1**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(1,1)*hl*h1*
     &    lambda1*st*sb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(1,1)*hl
     &    *h2*lambda2*st*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) - 16*hss(1,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*s*
     &    tx**(-1)*s2**(-1) - 8*hss(1,1)*hl**2*hr*st*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,1)*hl**2*hr*st*sb*mt*s*tx**(-2)
     &     + 8*hss(1,1)*hl**3*ct*cb*m1**2*mt*s*t1**(-1)*tx**(-2) - 8*
     &    hss(1,1)*hl**3*ct*cb*mt*u1*tx**(-2) - 16*hss(1,1)*hr*h1*
     &    lambda1*ct*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s1**(-1)
     &     - 16*hss(1,1)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*s*tx**(-1)*
     &    s1**(-1) - 16*hss(1,1)*hr*h2*lambda2*ct*cb*sqrt2**(-1)*m1**2*
     &    s*t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(1,1)*hr*h2*lambda2*ct*
     &    cb*sqrt2**(-1)*s*tx**(-1)*s2**(-1) + 8*hss(1,1)*hr**3*st*sb*
     &    m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(1,1)*hr**3*st*sb*mt*u1*
     &    tx**(-2) - 32*hss(1,2)*pq*hl*ssp*ct*sb*m1**2*mt*t1**(-1)*
     &    tx**(-1) + 32*hss(1,2)*pq*hl*ssp*ct*sb*mt*s**(-1)*u1*tx**(-1)
     &     + 32*hss(1,2)*pq*hr*ssp*st*cb*m1**2*mt*t1**(-1)*tx**(-1) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*hss(1,2)*pq
     &    *hr*ssp*st*cb*mt*s**(-1)*u1*tx**(-1) - 32*hss(1,2)*lq*hl*ssz*
     &    ct*sb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,2)*lq*
     &    hl*ssz*ct*sb*mt*u1*tx**(-1)*sz**(-1) + 32*hss(1,2)*rq*hr*ssz*
     &    st*cb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,2)*rq*
     &    hr*ssz*st*cb*mt*u1*tx**(-1)*sz**(-1) + 8*hss(1,2)*hl*hr**2*ct
     &    *sb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(1,2)*hl*hr**2*ct*sb*
     &    mt*s*tx**(-2) - 16*hss(1,2)*hl*h1*lambda1*st*cb*sqrt2**(-1)*
     &    m1**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(1,2)*hl*h1*
     &    lambda1*st*cb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) - 16*hss(1,2)*
     &    hl*h2*lambda2*st*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) - 16*hss(1,2)*hl*h2*lambda2*st*cb*sqrt2**(-1)*s*
     &    tx**(-1)*s2**(-1) - 8*hss(1,2)*hl**2*hr*st*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,2)*hl**2*hr*st*cb*mt*s*tx**(-2)
     &     - 8*hss(1,2)*hl**3*ct*sb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(1,2)*hl**3*ct*sb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(6,2)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(1,2)*hr*h1
     &    *lambda1*ct*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s1**(-1)
     &     + 16*hss(1,2)*hr*h1*lambda1*ct*sb*sqrt2**(-1)*s*tx**(-1)*
     &    s1**(-1) + 16*hss(1,2)*hr*h2*lambda2*ct*sb*sqrt2**(-1)*m1**2*
     &    s*t1**(-1)*tx**(-1)*s2**(-1) + 16*hss(1,2)*hr*h2*lambda2*ct*
     &    sb*sqrt2**(-1)*s*tx**(-1)*s2**(-1) + 8*hss(1,2)*hr**3*st*cb*
     &    m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(1,2)*hr**3*st*cb*mt*u1*
     &    tx**(-2) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hl**2*
     &    ssp*pt2*s2t*mg*mt*tx**(-2) + 4*pq*hl**2*ssp*pt2*c2t*m1**2*
     &    mg**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*c2t*m1**2*mst2**2
     &    *tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*c2t*m1**2*tx**(-2) - 4
     &    *pq*hl**2*ssp*pt2*c2t*mg**2*mt**2*tx**(-2)*t**(-1) + 4*pq*
     &    hl**2*ssp*pt2*c2t*mg**2*t1*tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*
     &    pt2*c2t*mt**2*mst2**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*
     &    c2t*mt**2*tx**(-2) - 4*pq*hl**2*ssp*pt2*c2t*mst2**2*t1*
     &    tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*pt2*c2t*t1*tx**(-2) - 4*pq*
     &    hl**2*ssp*pt2*m1**2*mg**2*tx**(-2)*t**(-1) + 4*pq*hl**2*ssp*
     &    pt2*m1**2*mst2**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*m1**2
     &    *tx**(-2) - 4*pq*hl**2*ssp*pt2*mg**2*mt**2*tx**(-2)*t**(-1)
     &     - 4*pq*hl**2*ssp*pt2*mg**2*t1*tx**(-2)*t**(-1) + 4*pq*hl**2*
     &    ssp*pt2*mt**2*mst2**2*tx**(-2)*t**(-1) - 4*pq*hl**2*ssp*pt2*
     &    mt**2*tx**(-2) + 4*pq*hl**2*ssp*pt2*mst2**2*t1*tx**(-2)*
     &    t**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hl**2*ssp
     &    *pt2*t1*tx**(-2) - 16*pq*hr**2*ssp*pt2*s2t*mg*mt*tx**(-2) - 4
     &    *pq*hr**2*ssp*pt2*c2t*m1**2*mg**2*tx**(-2)*t**(-1) + 4*pq*
     &    hr**2*ssp*pt2*c2t*m1**2*mst2**2*tx**(-2)*t**(-1) - 4*pq*hr**2
     &    *ssp*pt2*c2t*m1**2*tx**(-2) + 4*pq*hr**2*ssp*pt2*c2t*mg**2*
     &    mt**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*c2t*mg**2*t1*
     &    tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*c2t*mt**2*mst2**2*
     &    tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*c2t*mt**2*tx**(-2) + 4*
     &    pq*hr**2*ssp*pt2*c2t*mst2**2*t1*tx**(-2)*t**(-1) - 4*pq*hr**2
     &    *ssp*pt2*c2t*t1*tx**(-2) - 4*pq*hr**2*ssp*pt2*m1**2*mg**2*
     &    tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*m1**2*mst2**2*tx**(-2)*
     &    t**(-1) - 4*pq*hr**2*ssp*pt2*m1**2*tx**(-2) - 4*pq*hr**2*ssp*
     &    pt2*mg**2*mt**2*tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*mg**2*
     &    t1*tx**(-2)*t**(-1) + 4*pq*hr**2*ssp*pt2*mt**2*mst2**2*
     &    tx**(-2)*t**(-1) - 4*pq*hr**2*ssp*pt2*mt**2*tx**(-2) + 4*pq*
     &    hr**2*ssp*pt2*mst2**2*t1*tx**(-2)*t**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*pq*hr**2*ssp
     &    *pt2*t1*tx**(-2) - 16*lq*hl**2*ssz*pt2*s2t*mg*mt*s*tx**(-2)*
     &    sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*m1**2*mg**2*s*tx**(-2)*
     &    t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*m1**2*mst2**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*m1**2*s*
     &    tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*mg**2*mt**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*mg**2*s*t1
     &    *tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*mt**2*
     &    mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*
     &    mt**2*s*tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*c2t*mst2**2*s*
     &    t1*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*c2t*s*t1*
     &    tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*m1**2*mg**2*s*tx**(-2)
     &    *t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*pt2*m1**2*mst2**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*m1**2*s*
     &    tx**(-2)*sz**(-1) - 4*lq*hl**2*ssz*pt2*mg**2*mt**2*s*tx**(-2)
     &    *t**(-1)*sz**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*lq*hl**2*ssz
     &    *pt2*mg**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*lq*hl**2*ssz*
     &    pt2*mt**2*mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*
     &    ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) + 4*lq*hl**2*ssz*pt2*
     &    mst2**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*lq*hl**2*ssz*pt2*s
     &    *t1*tx**(-2)*sz**(-1) - 16*rq*hr**2*ssz*pt2*s2t*mg*mt*s*
     &    tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*m1**2*mg**2*s*
     &    tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*m1**2*
     &    mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*
     &    m1**2*s*tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*mg**2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t*
     &    mg**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*c2t
     &    *mt**2*mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*
     &    pt2*c2t*mt**2*s*tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*c2t*
     &    mst2**2*s*t1*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*
     &    c2t*s*t1*tx**(-2)*sz**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*rq*hr**2*ssz
     &    *pt2*m1**2*mg**2*s*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz
     &    *pt2*m1**2*mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*
     &    ssz*pt2*m1**2*s*tx**(-2)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mg**2*
     &    mt**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*mg**2*
     &    s*t1*tx**(-2)*t**(-1)*sz**(-1) + 4*rq*hr**2*ssz*pt2*mt**2*
     &    mst2**2*s*tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*
     &    mt**2*s*tx**(-2)*sz**(-1) + 4*rq*hr**2*ssz*pt2*mst2**2*s*t1*
     &    tx**(-2)*t**(-1)*sz**(-1) - 4*rq*hr**2*ssz*pt2*s*t1*tx**(-2)*
     &    sz**(-1) - 8*hl*hr*h1*lambda1*s2t*sqrt2**(-1)*m1**2*mg*s*
     &    tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*s2t*sqrt2**(-1)*mg*
     &    mt**2*s*tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*s2t*
     &    sqrt2**(-1)*mg*s*t1*tx**(-2)*s1**(-1) - 8*hl*hr*h1*lambda1*
     &    sqrt2**(-1)*m1**2*mg**2*mt*s*tx**(-2)*t**(-1)*s1**(-1) + 8*hl
     &    *hr*h1*lambda1*sqrt2**(-1)*m1**2*mt*mst2**2*s*tx**(-2)*
     &    t**(-1)*s1**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) - 8*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mg**2*mt*s*t1*tx**(-2)*t**(-1)*s1**(-1)
     &     + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*mst2**2*s*t1*tx**(-2)*
     &    t**(-1)*s1**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*s*t1*
     &    tx**(-2)*s1**(-1) - 8*hl*hr*h2*lambda2*s2t*sqrt2**(-1)*m1**2*
     &    mg*s*tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*s2t*sqrt2**(-1)*
     &    mg*mt**2*s*tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*s2t*
     &    sqrt2**(-1)*mg*s*t1*tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*m1**2*mg**2*mt*s*tx**(-2)*t**(-1)*s2**(-1) + 8*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*mst2**2*s*tx**(-2)*
     &    t**(-1)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s*
     &    tx**(-2)*s2**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*mg**2*mt*s
     &    *t1*tx**(-2)*t**(-1)*s2**(-1) + 8*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mt*mst2**2*s*t1*tx**(-2)*t**(-1)*s2**(-1) - 8*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*hl**2*hr**2*
     &    s2t*m1**2*mg*mt*s*tx**(-3) - 4*hl**2*hr**2*s2t*mg*mt*s*t1*
     &    tx**(-3) - 4*hl**2*hr**2*s2t*mg*mt**3*s*tx**(-3) - 4*hl**2*
     &    hr**2*m1**2*mg**2*mt**2*s*tx**(-3)*t**(-1) + 4*hl**2*hr**2*
     &    m1**2*mt**2*mst2**2*s*tx**(-3)*t**(-1) - 4*hl**2*hr**2*m1**2*
     &    mt**2*s*tx**(-3) - 4*hl**2*hr**2*mg**2*mt**2*s*t1*tx**(-3)*
     &    t**(-1) + 4*hl**2*hr**2*mt**2*mst2**2*s*t1*tx**(-3)*t**(-1)
     &     - 4*hl**2*hr**2*mt**2*s*t1*tx**(-3) - 4*hl**4*pt2*s2t*mg*mt*
     &    s*tx**(-3) + hl**4*pt2*c2t*m1**2*mg**2*s*tx**(-3)*t**(-1) - 
     &    hl**4*pt2*c2t*m1**2*mst2**2*s*tx**(-3)*t**(-1) + hl**4*pt2*
     &    c2t*m1**2*s*tx**(-3) - hl**4*pt2*c2t*mg**2*mt**2*s*tx**(-3)*
     &    t**(-1) + hl**4*pt2*c2t*mg**2*s*t1*tx**(-3)*t**(-1) + hl**4*
     &    pt2*c2t*mt**2*mst2**2*s*tx**(-3)*t**(-1) - hl**4*pt2*c2t*
     &    mt**2*s*tx**(-3) - hl**4*pt2*c2t*mst2**2*s*t1*tx**(-3)*
     &    t**(-1) + hl**4*pt2*c2t*s*t1*tx**(-3) - hl**4*pt2*m1**2*mg**2
     &    *s*tx**(-3)*t**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * ( hl**4*pt2*m1**2*
     &    mst2**2*s*tx**(-3)*t**(-1) - hl**4*pt2*m1**2*s*tx**(-3) - 
     &    hl**4*pt2*mg**2*mt**2*s*tx**(-3)*t**(-1) - hl**4*pt2*mg**2*s*
     &    t1*tx**(-3)*t**(-1) + hl**4*pt2*mt**2*mst2**2*s*tx**(-3)*
     &    t**(-1) - hl**4*pt2*mt**2*s*tx**(-3) + hl**4*pt2*mst2**2*s*t1
     &    *tx**(-3)*t**(-1) - hl**4*pt2*s*t1*tx**(-3) - 4*hr**4*pt2*s2t
     &    *mg*mt*s*tx**(-3) - hr**4*pt2*c2t*m1**2*mg**2*s*tx**(-3)*
     &    t**(-1) + hr**4*pt2*c2t*m1**2*mst2**2*s*tx**(-3)*t**(-1) - 
     &    hr**4*pt2*c2t*m1**2*s*tx**(-3) + hr**4*pt2*c2t*mg**2*mt**2*s*
     &    tx**(-3)*t**(-1) - hr**4*pt2*c2t*mg**2*s*t1*tx**(-3)*t**(-1)
     &     - hr**4*pt2*c2t*mt**2*mst2**2*s*tx**(-3)*t**(-1) + hr**4*pt2
     &    *c2t*mt**2*s*tx**(-3) + hr**4*pt2*c2t*mst2**2*s*t1*tx**(-3)*
     &    t**(-1) - hr**4*pt2*c2t*s*t1*tx**(-3) - hr**4*pt2*m1**2*mg**2
     &    *s*tx**(-3)*t**(-1) + hr**4*pt2*m1**2*mst2**2*s*tx**(-3)*
     &    t**(-1) - hr**4*pt2*m1**2*s*tx**(-3) - hr**4*pt2*mg**2*mt**2*
     &    s*tx**(-3)*t**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - hr**4*pt2*
     &    mg**2*s*t1*tx**(-3)*t**(-1) + hr**4*pt2*mt**2*mst2**2*s*
     &    tx**(-3)*t**(-1) - hr**4*pt2*mt**2*s*tx**(-3) + hr**4*pt2*
     &    mst2**2*s*t1*tx**(-3)*t**(-1) - hr**4*pt2*s*t1*tx**(-3) - 32*
     &    hss(2,1)*pq*hl*ssp*st*cb*m1**2*mt*t1**(-1)*tx**(-1) + 32*hss(
     &    2,1)*pq*hl*ssp*st*cb*mt*s**(-1)*u1*tx**(-1) + 32*hss(2,1)*pq*
     &    hr*ssp*ct*sb*m1**2*mt*t1**(-1)*tx**(-1) - 32*hss(2,1)*pq*hr*
     &    ssp*ct*sb*mt*s**(-1)*u1*tx**(-1) - 32*hss(2,1)*lq*hl*ssz*st*
     &    cb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,1)*lq*hl*
     &    ssz*st*cb*mt*u1*tx**(-1)*sz**(-1) + 32*hss(2,1)*rq*hr*ssz*ct*
     &    sb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(2,1)*rq*hr*
     &    ssz*ct*sb*mt*u1*tx**(-1)*sz**(-1) + 8*hss(2,1)*hl*hr**2*st*cb
     &    *m1**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(2,1)*hl*hr**2*st*cb*mt*
     &    s*tx**(-2) - 16*hss(2,1)*hl*h1*lambda1*ct*sb*sqrt2**(-1)*
     &    m1**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(2,1)*hl*h1*
     &    lambda1*ct*sb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(2,1)*hl
     &    *h2*lambda2*ct*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) - 16*hss(2,1)*hl*h2*lambda2*ct*sb*sqrt2**(-1)*s*
     &    tx**(-1)*s2**(-1) - 8*hss(2,1)*hl**2*hr*ct*sb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(2,1)*hl**2*hr*ct*sb*mt*s*tx**(-2)
     &     - 8*hss(2,1)*hl**3*st*cb*m1**2*mt*s*t1**(-1)*tx**(-2) + 8*
     &    hss(2,1)*hl**3*st*cb*mt*u1*tx**(-2) + 16*hss(2,1)*hr*h1*
     &    lambda1*st*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s1**(-1)
     &     + 16*hss(2,1)*hr*h1*lambda1*st*cb*sqrt2**(-1)*s*tx**(-1)*
     &    s1**(-1) + 16*hss(2,1)*hr*h2*lambda2*st*cb*sqrt2**(-1)*m1**2*
     &    s*t1**(-1)*tx**(-1)*s2**(-1) + 16*hss(2,1)*hr*h2*lambda2*st*
     &    cb*sqrt2**(-1)*s*tx**(-1)*s2**(-1) + 8*hss(2,1)*hr**3*ct*sb*
     &    m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(2,1)*hr**3*ct*sb*mt*u1*
     &    tx**(-2) + 32*hss(2,2)*pq*hl*ssp*st*sb*m1**2*mt*t1**(-1)*
     &    tx**(-1) - 32*hss(2,2)*pq*hl*ssp*st*sb*mt*s**(-1)*u1*tx**(-1)
     &     + 32*hss(2,2)*pq*hr*ssp*ct*cb*m1**2*mt*t1**(-1)*tx**(-1) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 32*hss(2,2)*pq
     &    *hr*ssp*ct*cb*mt*s**(-1)*u1*tx**(-1) + 32*hss(2,2)*lq*hl*ssz*
     &    st*sb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(2,2)*lq*
     &    hl*ssz*st*sb*mt*u1*tx**(-1)*sz**(-1) + 32*hss(2,2)*rq*hr*ssz*
     &    ct*cb*m1**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(2,2)*rq*
     &    hr*ssz*ct*cb*mt*u1*tx**(-1)*sz**(-1) - 8*hss(2,2)*hl*hr**2*st
     &    *sb*m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(2,2)*hl*hr**2*st*sb*
     &    mt*s*tx**(-2) - 16*hss(2,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*
     &    m1**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(2,2)*hl*h1*
     &    lambda1*ct*cb*sqrt2**(-1)*s*tx**(-1)*s1**(-1) - 16*hss(2,2)*
     &    hl*h2*lambda2*ct*cb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) - 16*hss(2,2)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*s*
     &    tx**(-1)*s2**(-1) - 8*hss(2,2)*hl**2*hr*ct*cb*m1**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(2,2)*hl**2*hr*ct*cb*mt*s*tx**(-2)
     &     + 8*hss(2,2)*hl**3*st*sb*m1**2*mt*s*t1**(-1)*tx**(-2) - 8*
     &    hss(2,2)*hl**3*st*sb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(6,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(2,2)*hr
     &    *h1*lambda1*st*sb*sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) - 16*hss(2,2)*hr*h1*lambda1*st*sb*sqrt2**(-1)*s*
     &    tx**(-1)*s1**(-1) - 16*hss(2,2)*hr*h2*lambda2*st*sb*
     &    sqrt2**(-1)*m1**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(2,2)*
     &    hr*h2*lambda2*st*sb*sqrt2**(-1)*s*tx**(-1)*s2**(-1) + 8*hss(2
     &    ,2)*hr**3*ct*cb*m1**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(2,2)*
     &    hr**3*ct*cb*mt*u1*tx**(-2) )
      MMs = MMs + SCB(7,2)*Nc*Cf*Pi*alphas*prefac * ( 8*pq*hl**2*ssp*
     &    pt2*mg**2*tx**(-2) + 8*pq*hl**2*ssp*pt2*mt**2*tx**(-2) - 8*pq
     &    *hl**2*ssp*pt2*mst1**2*tx**(-2) + 8*pq*hr**2*ssp*pt2*mg**2*
     &    tx**(-2) + 8*pq*hr**2*ssp*pt2*mt**2*tx**(-2) - 8*pq*hr**2*ssp
     &    *pt2*mst1**2*tx**(-2) + 8*lq*hl**2*ssz*pt2*mg**2*s*tx**(-2)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) - 8*
     &    lq*hl**2*ssz*pt2*mst1**2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz
     &    *pt2*mg**2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz*pt2*mt**2*s*
     &    tx**(-2)*sz**(-1) - 8*rq*hr**2*ssz*pt2*mst1**2*s*tx**(-2)*
     &    sz**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mg**2*
     &    mt**(-1)*s*tx**(-2)*s1**(-1) - 4*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *m1**2*mt**(-1)*mst1**2*s*tx**(-2)*s1**(-1) + 4*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) + 4*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mg**2*mt**(-1)*s*t1*tx**(-2)*s1**(-1) + 
     &    4*hl*hr*h1*lambda1*sqrt2**(-1)*mg**2*mt*s*tx**(-2)*s1**(-1)
     &     - 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt**(-1)*mst1**2*s*t1*
     &    tx**(-2)*s1**(-1) )
      MMs = MMs + SCB(7,2)*Nc*Cf*Pi*alphas*prefac * (  - 4*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*mst1**2*s*tx**(-2)*s1**(-1) + 4*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*s1**(-1) + 4*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mt**3*s*tx**(-2)*s1**(-1) + 4*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mg**2*mt**(-1)*s*tx**(-2)*s2**(-1)
     &     - 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**(-1)*mst1**2*s*
     &    tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s
     &    *tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*mg**2*
     &    mt**(-1)*s*t1*tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mg**2*mt*s*tx**(-2)*s2**(-1) - 4*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt**(-1)*mst1**2*s*t1*tx**(-2)*s2**(-1) - 4*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mst1**2*s*tx**(-2)*s2**(-1) + 4*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) + 4*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*mt**3*s*tx**(-2)*s2**(-1) + 2*
     &    hl**2*hr**2*m1**2*mg**2*s*tx**(-3) + 2*hl**2*hr**2*m1**2*
     &    mt**2*s*tx**(-3) )
      MMs = MMs + SCB(7,2)*Nc*Cf*Pi*alphas*prefac * (  - 2*hl**2*hr**2*
     &    m1**2*mst1**2*s*tx**(-3) + 2*hl**2*hr**2*mg**2*mt**2*s*
     &    tx**(-3) + 2*hl**2*hr**2*mg**2*s*t1*tx**(-3) - 2*hl**2*hr**2*
     &    mt**2*mst1**2*s*tx**(-3) + 2*hl**2*hr**2*mt**2*s*t1*tx**(-3)
     &     + 2*hl**2*hr**2*mt**4*s*tx**(-3) - 2*hl**2*hr**2*mst1**2*s*
     &    t1*tx**(-3) + 2*hl**4*pt2*mg**2*s*tx**(-3) + 2*hl**4*pt2*
     &    mt**2*s*tx**(-3) - 2*hl**4*pt2*mst1**2*s*tx**(-3) + 2*hr**4*
     &    pt2*mg**2*s*tx**(-3) + 2*hr**4*pt2*mt**2*s*tx**(-3) - 2*hr**4
     &    *pt2*mst1**2*s*tx**(-3) )
      MMs = MMs + SCB(7,3)*Nc*Cf*Pi*alphas*prefac * ( 8*pq*hl**2*ssp*
     &    pt2*mg**2*tx**(-2) + 8*pq*hl**2*ssp*pt2*mt**2*tx**(-2) - 8*pq
     &    *hl**2*ssp*pt2*mst2**2*tx**(-2) + 8*pq*hr**2*ssp*pt2*mg**2*
     &    tx**(-2) + 8*pq*hr**2*ssp*pt2*mt**2*tx**(-2) - 8*pq*hr**2*ssp
     &    *pt2*mst2**2*tx**(-2) + 8*lq*hl**2*ssz*pt2*mg**2*s*tx**(-2)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*mt**2*s*tx**(-2)*sz**(-1) - 8*
     &    lq*hl**2*ssz*pt2*mst2**2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz
     &    *pt2*mg**2*s*tx**(-2)*sz**(-1) + 8*rq*hr**2*ssz*pt2*mt**2*s*
     &    tx**(-2)*sz**(-1) - 8*rq*hr**2*ssz*pt2*mst2**2*s*tx**(-2)*
     &    sz**(-1) + 4*hl*hr*h1*lambda1*sqrt2**(-1)*m1**2*mg**2*
     &    mt**(-1)*s*tx**(-2)*s1**(-1) - 4*hl*hr*h1*lambda1*sqrt2**(-1)
     &    *m1**2*mt**(-1)*mst2**2*s*tx**(-2)*s1**(-1) + 4*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*m1**2*mt*s*tx**(-2)*s1**(-1) + 4*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mg**2*mt**(-1)*s*t1*tx**(-2)*s1**(-1) + 
     &    4*hl*hr*h1*lambda1*sqrt2**(-1)*mg**2*mt*s*tx**(-2)*s1**(-1)
     &     - 4*hl*hr*h1*lambda1*sqrt2**(-1)*mt**(-1)*mst2**2*s*t1*
     &    tx**(-2)*s1**(-1) )
      MMs = MMs + SCB(7,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*hl*hr*h1*
     &    lambda1*sqrt2**(-1)*mt*mst2**2*s*tx**(-2)*s1**(-1) + 4*hl*hr*
     &    h1*lambda1*sqrt2**(-1)*mt*s*t1*tx**(-2)*s1**(-1) + 4*hl*hr*h1
     &    *lambda1*sqrt2**(-1)*mt**3*s*tx**(-2)*s1**(-1) + 4*hl*hr*h2*
     &    lambda2*sqrt2**(-1)*m1**2*mg**2*mt**(-1)*s*tx**(-2)*s2**(-1)
     &     - 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt**(-1)*mst2**2*s*
     &    tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*m1**2*mt*s
     &    *tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*sqrt2**(-1)*mg**2*
     &    mt**(-1)*s*t1*tx**(-2)*s2**(-1) + 4*hl*hr*h2*lambda2*
     &    sqrt2**(-1)*mg**2*mt*s*tx**(-2)*s2**(-1) - 4*hl*hr*h2*lambda2
     &    *sqrt2**(-1)*mt**(-1)*mst2**2*s*t1*tx**(-2)*s2**(-1) - 4*hl*
     &    hr*h2*lambda2*sqrt2**(-1)*mt*mst2**2*s*tx**(-2)*s2**(-1) + 4*
     &    hl*hr*h2*lambda2*sqrt2**(-1)*mt*s*t1*tx**(-2)*s2**(-1) + 4*hl
     &    *hr*h2*lambda2*sqrt2**(-1)*mt**3*s*tx**(-2)*s2**(-1) + 2*
     &    hl**2*hr**2*m1**2*mg**2*s*tx**(-3) + 2*hl**2*hr**2*m1**2*
     &    mt**2*s*tx**(-3) )
      MMs = MMs + SCB(7,3)*Nc*Cf*Pi*alphas*prefac * (  - 2*hl**2*hr**2*
     &    m1**2*mst2**2*s*tx**(-3) + 2*hl**2*hr**2*mg**2*mt**2*s*
     &    tx**(-3) + 2*hl**2*hr**2*mg**2*s*t1*tx**(-3) - 2*hl**2*hr**2*
     &    mt**2*mst2**2*s*tx**(-3) + 2*hl**2*hr**2*mt**2*s*t1*tx**(-3)
     &     + 2*hl**2*hr**2*mt**4*s*tx**(-3) - 2*hl**2*hr**2*mst2**2*s*
     &    t1*tx**(-3) + 2*hl**4*pt2*mg**2*s*tx**(-3) + 2*hl**4*pt2*
     &    mt**2*s*tx**(-3) - 2*hl**4*pt2*mst2**2*s*tx**(-3) + 2*hr**4*
     &    pt2*mg**2*s*tx**(-3) + 2*hr**4*pt2*mt**2*s*tx**(-3) - 2*hr**4
     &    *pt2*mst2**2*s*tx**(-3) )
ctp      print*, " SCB ",MMs
      MMs = MMs + SCBP(1)*Nc*Cf*Pi*alphas*prefac * ( 64*pq*lq*ssz*ssp
     &    *pt2*c2b*mg**2*sz**(-1) - 64*pq*lq*ssz*ssp*pt2*c2b*msb1**2*
     &    sz**(-1) - 32*pq*lq*ssz*ssp*pt2*mg**2*sz**(-1) + 32*pq*lq*ssz
     &    *ssp*pt2*msb1**2*sz**(-1) - 64*pq*rq*ssz*ssp*pt2*c2b*mg**2*
     &    sz**(-1) + 64*pq*rq*ssz*ssp*pt2*c2b*msb1**2*sz**(-1) - 32*pq*
     &    rq*ssz*ssp*pt2*mg**2*sz**(-1) + 32*pq*rq*ssz*ssp*pt2*msb1**2*
     &    sz**(-1) + 16*pq*hl**2*ssp*pt2*c2b*mg**2*tx**(-1) - 16*pq*
     &    hl**2*ssp*pt2*c2b*msb1**2*tx**(-1) - 8*pq*hl**2*ssp*pt2*mg**2
     &    *tx**(-1) + 8*pq*hl**2*ssp*pt2*msb1**2*tx**(-1) - 16*pq*hr**2
     &    *ssp*pt2*c2b*mg**2*tx**(-1) + 16*pq*hr**2*ssp*pt2*c2b*msb1**2
     &    *tx**(-1) - 8*pq*hr**2*ssp*pt2*mg**2*tx**(-1) + 8*pq*hr**2*
     &    ssp*pt2*msb1**2*tx**(-1) + 16*lq*hl**2*ssz*pt2*c2b*mg**2*s*
     &    tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*pt2*c2b*msb1**2*s*
     &    tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*pt2*mg**2*s*tx**(-1)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*msb1**2*s*tx**(-1)*sz**(-1) - 
     &    16*rq*hr**2*ssz*pt2*c2b*mg**2*s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCBP(1)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hr**2*ssz*
     &    pt2*c2b*msb1**2*s*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*pt2*
     &    mg**2*s*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*msb1**2*s*
     &    tx**(-1)*sz**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mg**2*mt*s
     &    *tx**(-1)*s1**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    msb1**2*s*tx**(-1)*s1**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mg**2*mt*s*tx**(-1)*s2**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt*msb1**2*s*tx**(-1)*s2**(-1) - 2*hl**2*hr**2*mg**2*mt**2*s
     &    *tx**(-2) + 2*hl**2*hr**2*mt**2*msb1**2*s*tx**(-2) + 2*hl**4*
     &    pt2*c2b*mg**2*s*tx**(-2) - 2*hl**4*pt2*c2b*msb1**2*s*tx**(-2)
     &     - hl**4*pt2*mg**2*s*tx**(-2) + hl**4*pt2*msb1**2*s*tx**(-2)
     &     - 2*hr**4*pt2*c2b*mg**2*s*tx**(-2) + 2*hr**4*pt2*c2b*msb1**2
     &    *s*tx**(-2) - hr**4*pt2*mg**2*s*tx**(-2) + hr**4*pt2*msb1**2*
     &    s*tx**(-2) - 8*h1*h2*lambda1*lambda2*mg**2*s*s1**(-1)*
     &    s2**(-1) + 8*h1*h2*lambda1*lambda2*msb1**2*s*s1**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCBP(1)*Nc*Cf*Pi*alphas*prefac * (  - 4*h1**2*
     &    lambda1**2*mg**2*s*s1**(-2) + 4*h1**2*lambda1**2*msb1**2*s*
     &    s1**(-2) - 4*h2**2*lambda2**2*mg**2*s*s2**(-2) + 4*h2**2*
     &    lambda2**2*msb1**2*s*s2**(-2) + 32*ssz**2*lq2*pt2*c2b*mg**2*s
     &    *sz**(-2) - 32*ssz**2*lq2*pt2*c2b*msb1**2*s*sz**(-2) - 16*
     &    ssz**2*lq2*pt2*mg**2*s*sz**(-2) + 16*ssz**2*lq2*pt2*msb1**2*s
     &    *sz**(-2) - 32*ssz**2*rq2*pt2*c2b*mg**2*s*sz**(-2) + 32*
     &    ssz**2*rq2*pt2*c2b*msb1**2*s*sz**(-2) - 16*ssz**2*rq2*pt2*
     &    mg**2*s*sz**(-2) + 16*ssz**2*rq2*pt2*msb1**2*s*sz**(-2) - 32*
     &    ssp**2*pq2*pt2*mg**2*s**(-1) + 32*ssp**2*pq2*pt2*msb1**2*
     &    s**(-1) )
      MMs = MMs + SCBP(2)*Nc*Cf*Pi*alphas*prefac * (  - 64*pq*lq*ssz*
     &    ssp*pt2*c2b*mg**2*sz**(-1) + 64*pq*lq*ssz*ssp*pt2*c2b*msb2**2
     &    *sz**(-1) - 32*pq*lq*ssz*ssp*pt2*mg**2*sz**(-1) + 32*pq*lq*
     &    ssz*ssp*pt2*msb2**2*sz**(-1) + 64*pq*rq*ssz*ssp*pt2*c2b*mg**2
     &    *sz**(-1) - 64*pq*rq*ssz*ssp*pt2*c2b*msb2**2*sz**(-1) - 32*pq
     &    *rq*ssz*ssp*pt2*mg**2*sz**(-1) + 32*pq*rq*ssz*ssp*pt2*msb2**2
     &    *sz**(-1) - 16*pq*hl**2*ssp*pt2*c2b*mg**2*tx**(-1) + 16*pq*
     &    hl**2*ssp*pt2*c2b*msb2**2*tx**(-1) - 8*pq*hl**2*ssp*pt2*mg**2
     &    *tx**(-1) + 8*pq*hl**2*ssp*pt2*msb2**2*tx**(-1) + 16*pq*hr**2
     &    *ssp*pt2*c2b*mg**2*tx**(-1) - 16*pq*hr**2*ssp*pt2*c2b*msb2**2
     &    *tx**(-1) - 8*pq*hr**2*ssp*pt2*mg**2*tx**(-1) + 8*pq*hr**2*
     &    ssp*pt2*msb2**2*tx**(-1) - 16*lq*hl**2*ssz*pt2*c2b*mg**2*s*
     &    tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*pt2*c2b*msb2**2*s*
     &    tx**(-1)*sz**(-1) - 8*lq*hl**2*ssz*pt2*mg**2*s*tx**(-1)*
     &    sz**(-1) + 8*lq*hl**2*ssz*pt2*msb2**2*s*tx**(-1)*sz**(-1) + 
     &    16*rq*hr**2*ssz*pt2*c2b*mg**2*s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCBP(2)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*
     &    ssz*pt2*c2b*msb2**2*s*tx**(-1)*sz**(-1) - 8*rq*hr**2*ssz*pt2*
     &    mg**2*s*tx**(-1)*sz**(-1) + 8*rq*hr**2*ssz*pt2*msb2**2*s*
     &    tx**(-1)*sz**(-1) - 8*hl*hr*h1*lambda1*sqrt2**(-1)*mg**2*mt*s
     &    *tx**(-1)*s1**(-1) + 8*hl*hr*h1*lambda1*sqrt2**(-1)*mt*
     &    msb2**2*s*tx**(-1)*s1**(-1) - 8*hl*hr*h2*lambda2*sqrt2**(-1)*
     &    mg**2*mt*s*tx**(-1)*s2**(-1) + 8*hl*hr*h2*lambda2*sqrt2**(-1)
     &    *mt*msb2**2*s*tx**(-1)*s2**(-1) - 2*hl**2*hr**2*mg**2*mt**2*s
     &    *tx**(-2) + 2*hl**2*hr**2*mt**2*msb2**2*s*tx**(-2) - 2*hl**4*
     &    pt2*c2b*mg**2*s*tx**(-2) + 2*hl**4*pt2*c2b*msb2**2*s*tx**(-2)
     &     - hl**4*pt2*mg**2*s*tx**(-2) + hl**4*pt2*msb2**2*s*tx**(-2)
     &     + 2*hr**4*pt2*c2b*mg**2*s*tx**(-2) - 2*hr**4*pt2*c2b*msb2**2
     &    *s*tx**(-2) - hr**4*pt2*mg**2*s*tx**(-2) + hr**4*pt2*msb2**2*
     &    s*tx**(-2) - 8*h1*h2*lambda1*lambda2*mg**2*s*s1**(-1)*
     &    s2**(-1) + 8*h1*h2*lambda1*lambda2*msb2**2*s*s1**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCBP(2)*Nc*Cf*Pi*alphas*prefac * (  - 4*h1**2*
     &    lambda1**2*mg**2*s*s1**(-2) + 4*h1**2*lambda1**2*msb2**2*s*
     &    s1**(-2) - 4*h2**2*lambda2**2*mg**2*s*s2**(-2) + 4*h2**2*
     &    lambda2**2*msb2**2*s*s2**(-2) - 32*ssz**2*lq2*pt2*c2b*mg**2*s
     &    *sz**(-2) + 32*ssz**2*lq2*pt2*c2b*msb2**2*s*sz**(-2) - 16*
     &    ssz**2*lq2*pt2*mg**2*s*sz**(-2) + 16*ssz**2*lq2*pt2*msb2**2*s
     &    *sz**(-2) + 32*ssz**2*rq2*pt2*c2b*mg**2*s*sz**(-2) - 32*
     &    ssz**2*rq2*pt2*c2b*msb2**2*s*sz**(-2) - 16*ssz**2*rq2*pt2*
     &    mg**2*s*sz**(-2) + 16*ssz**2*rq2*pt2*msb2**2*s*sz**(-2) - 32*
     &    ssp**2*pq2*pt2*mg**2*s**(-1) + 32*ssp**2*pq2*pt2*msb2**2*
     &    s**(-1) )
ctp      print*, " SCBP ",MMs
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 256*pq*lq*ssz*ssp
     &    *cb**2*m1**2*mg**2*msb1**2*s**(-1)*sz**(-1) - 128*pq*lq*ssz*
     &    ssp*cb**2*m1**2*mg**2*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*
     &    m1**2*mg**4*s**(-1)*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*m1**2*
     &    msb1**4*s**(-1)*sz**(-1) - 256*pq*lq*ssz*ssp*cb**2*mg**2*
     &    msb1**2*s**(-2)*t1*u1*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*
     &    mg**2*s**(-1)*t1*u1*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*mg**4*
     &    s**(-2)*t1*u1*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*msb1**4*
     &    s**(-2)*t1*u1*sz**(-1) - 256*pq*rq*ssz*ssp*cb**2*m1**2*mg**2*
     &    msb1**2*s**(-1)*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*m1**2*
     &    mg**2*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*m1**2*mg**4*s**(-1)*
     &    sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*m1**2*msb1**4*s**(-1)*
     &    sz**(-1) + 256*pq*rq*ssz*ssp*cb**2*mg**2*msb1**2*s**(-2)*t1*
     &    u1*sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*mg**2*s**(-1)*t1*u1*
     &    sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*mg**4*s**(-2)*t1*u1*
     &    sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 128*pq*rq*ssz*
     &    ssp*cb**2*msb1**4*s**(-2)*t1*u1*sz**(-1) + 256*pq*rq*ssz*ssp*
     &    m1**2*mg**2*msb1**2*s**(-1)*sz**(-1) - 128*pq*rq*ssz*ssp*
     &    m1**2*mg**2*sz**(-1) - 128*pq*rq*ssz*ssp*m1**2*mg**4*s**(-1)*
     &    sz**(-1) - 128*pq*rq*ssz*ssp*m1**2*msb1**4*s**(-1)*sz**(-1)
     &     - 256*pq*rq*ssz*ssp*mg**2*msb1**2*s**(-2)*t1*u1*sz**(-1) + 
     &    128*pq*rq*ssz*ssp*mg**2*s**(-1)*t1*u1*sz**(-1) + 128*pq*rq*
     &    ssz*ssp*mg**4*s**(-2)*t1*u1*sz**(-1) + 128*pq*rq*ssz*ssp*
     &    msb1**4*s**(-2)*t1*u1*sz**(-1) - 32*pq*hl*hr*ssp*sb*cb*mg*mt*
     &    msb1**2*s**(-1)*t1*tx**(-1) + 32*pq*hl*hr*ssp*sb*cb*mg*mt*
     &    msb1**2*s**(-1)*u1*tx**(-1) + 16*pq*hl*hr*ssp*sb*cb*mg*mt*t1*
     &    tx**(-1) - 16*pq*hl*hr*ssp*sb*cb*mg*mt*u1*tx**(-1) + 32*pq*hl
     &    *hr*ssp*sb*cb*mg**3*mt*s**(-1)*t1*tx**(-1) - 32*pq*hl*hr*ssp*
     &    sb*cb*mg**3*mt*s**(-1)*u1*tx**(-1) + 32*pq*hl**2*ssp*cb**2*
     &    m1**2*mg**2*msb1**2*s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*cb**2*
     &    m1**2*mg**2*tx**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hl**2*
     &    ssp*cb**2*m1**2*mg**4*s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*
     &    cb**2*m1**2*msb1**4*s**(-1)*tx**(-1) - 32*pq*hl**2*ssp*cb**2*
     &    mg**2*msb1**2*s**(-2)*t1*u1*tx**(-1) + 16*pq*hl**2*ssp*cb**2*
     &    mg**2*s**(-1)*t1*u1*tx**(-1) + 16*pq*hl**2*ssp*cb**2*mg**4*
     &    s**(-2)*t1*u1*tx**(-1) + 16*pq*hl**2*ssp*cb**2*msb1**4*
     &    s**(-2)*t1*u1*tx**(-1) - 32*pq*hr**2*ssp*cb**2*m1**2*mg**2*
     &    msb1**2*s**(-1)*tx**(-1) + 16*pq*hr**2*ssp*cb**2*m1**2*mg**2*
     &    tx**(-1) + 16*pq*hr**2*ssp*cb**2*m1**2*mg**4*s**(-1)*tx**(-1)
     &     + 16*pq*hr**2*ssp*cb**2*m1**2*msb1**4*s**(-1)*tx**(-1) + 32*
     &    pq*hr**2*ssp*cb**2*mg**2*msb1**2*s**(-2)*t1*u1*tx**(-1) - 16*
     &    pq*hr**2*ssp*cb**2*mg**2*s**(-1)*t1*u1*tx**(-1) - 16*pq*hr**2
     &    *ssp*cb**2*mg**4*s**(-2)*t1*u1*tx**(-1) - 16*pq*hr**2*ssp*
     &    cb**2*msb1**4*s**(-2)*t1*u1*tx**(-1) + 32*pq*hr**2*ssp*m1**2*
     &    mg**2*msb1**2*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp*m1**2*mg**2*
     &    tx**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hr**2*
     &    ssp*m1**2*mg**4*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp*m1**2*
     &    msb1**4*s**(-1)*tx**(-1) - 32*pq*hr**2*ssp*mg**2*msb1**2*
     &    s**(-2)*t1*u1*tx**(-1) + 16*pq*hr**2*ssp*mg**2*s**(-1)*t1*u1*
     &    tx**(-1) + 16*pq*hr**2*ssp*mg**4*s**(-2)*t1*u1*tx**(-1) + 16*
     &    pq*hr**2*ssp*msb1**4*s**(-2)*t1*u1*tx**(-1) - 64*pq*h1*ssp*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*s**(-1)*t1*s1**(-1) + 64
     &    *pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*s**(-1)*u1*
     &    s1**(-1) + 32*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*t1*
     &    s1**(-1) - 32*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*u1*
     &    s1**(-1) + 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg**3*
     &    s**(-1)*t1*s1**(-1) - 64*pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*
     &    mg**3*s**(-1)*u1*s1**(-1) - 64*pq*h2*ssp*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*msb1**2*s**(-1)*t1*s2**(-1) + 64*pq*h2*ssp*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*s**(-1)*u1*s2**(-1) + 32
     &    *pq*h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg*t1*s2**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 32*pq*h2*ssp*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*u1*s2**(-1) + 64*pq*h2*ssp*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*t1*s2**(-1) - 64*pq*
     &    h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*u1*s2**(-1) + 
     &    256*lq*rq*ssz**2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) - 128*lq*
     &    rq*ssz**2*cb**2*m1**2*mg**2*s*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**4*sz**(-2) - 128*lq*rq*ssz**2*cb**2*m1**2*
     &    msb1**4*sz**(-2) - 256*lq*rq*ssz**2*cb**2*mg**2*msb1**2*
     &    s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*t1*u1*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**4*s**(-1)*t1*u1*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*msb1**4*s**(-1)*t1*u1*
     &    sz**(-2) - 256*lq*rq*ssz**2*cb**4*m1**2*mg**2*msb1**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*mg**2*s*sz**(-2) + 
     &    128*lq*rq*ssz**2*cb**4*m1**2*mg**4*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**4*m1**2*msb1**4*sz**(-2) + 256*lq*rq*ssz**2*cb**4*
     &    mg**2*msb1**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**4*mg**2*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*msb1**4
     &    *s**(-1)*t1*u1*sz**(-2) - 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*
     &    msb1**2*t1*tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*mg*mt
     &    *msb1**2*u1*tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb**3*mg*
     &    mt*s*t1*tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*s*
     &    u1*tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*cb**4*m1**2*mg**2*msb1**2
     &    *tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*mg**2*s*
     &    tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*mg**4*
     &    tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*m1**2*msb1**4*
     &    tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**4*mg**2*msb1**2*
     &    s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**4*mg**2
     &    *t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hl**2*ssz*
     &    cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz
     &    *cb**4*msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*lq*hr**2*
     &    ssz*cb**2*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) - 16*lq*hr**2
     &    *ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) - 32*lq*hr**2*ssz*cb**2*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*
     &    cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2
     &    *msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*lq*hr**2*ssz*
     &    cb**4*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz
     &    *cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) + 32*lq*hr**2*ssz*cb**4*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*lq*hr**2*
     &    ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*
     &    cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz
     &    *cb**4*msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 64*lq*h1*ssz
     &    *lambda1*sb*cb**3*sqrt2**(-1)*mg*msb1**2*t1*sz**(-1)*s1**(-1)
     &     + 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb1**2*u1*
     &    sz**(-1)*s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1*
     &    sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) - 64*lq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1)
     &     - 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*t1*
     &    sz**(-1)*s2**(-1) + 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)
     &    *mg*msb1**2*u1*sz**(-1)*s2**(-1) + 32*lq*h2*ssz*lambda2*sb*
     &    cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 64*lq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) - 64*
     &    lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*
     &    s2**(-1) - 32*rq*hl*hr*ssz*sb*cb*mg*mt*msb1**2*t1*tx**(-1)*
     &    sz**(-1) + 32*rq*hl*hr*ssz*sb*cb*mg*mt*msb1**2*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hl*hr*ssz*sb*cb*mg*mt*s*t1*tx**(-1)*sz**(-1)
     &     - 16*rq*hl*hr*ssz*sb*cb*mg*mt*s*u1*tx**(-1)*sz**(-1) + 32*rq
     &    *hl*hr*ssz*sb*cb*mg**3*mt*t1*tx**(-1)*sz**(-1) - 32*rq*hl*hr*
     &    ssz*sb*cb*mg**3*mt*u1*tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*msb1**2*t1*tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb
     &    *cb**3*mg*mt*msb1**2*u1*tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*
     &    sb*cb**3*mg*mt*s*t1*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*s*u1*tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3
     &    *mg**3*mt*t1*tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*
     &    mg**3*mt*u1*tx**(-1)*sz**(-1) + 32*rq*hl**2*ssz*cb**2*m1**2*
     &    mg**2*msb1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hl**2*
     &    ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*cb**2*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2
     &    *msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*
     &    cb**4*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz
     &    *cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) + 32*rq*hl**2*ssz*cb**4*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4
     &    *msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 64*rq*hr**2*
     &    ssz*cb**2*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 32*rq*hr**2
     &    *ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*cb**2*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) + 64*rq*hr**2*ssz*cb**2*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*
     &    cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*cb**2*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*cb**2
     &    *msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*rq*hr**2*ssz*
     &    cb**4*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz
     &    *cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*
     &    m1**2*msb1**4*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*cb**4*mg**2
     &    *msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hr**2*ssz*
     &    cb**4*msb1**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*rq*hr**2*
     &    ssz*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*
     &    m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*m1**2*msb1**4*tx**(-1)*
     &    sz**(-1) - 32*rq*hr**2*ssz*mg**2*msb1**2*s**(-1)*t1*u1*
     &    tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*mg**4*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*msb1**4*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 64*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*
     &    t1*sz**(-1)*s1**(-1) + 64*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)
     &    *mg*msb1**2*u1*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb
     &    *sqrt2**(-1)*mg*s*t1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1
     &    *sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 64*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) - 64*rq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 64*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*msb1**2*t1*sz**(-1)*s1**(-1)
     &     - 64*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb1**2*u1*
     &    sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) - 64*rq*h1*ssz*lambda1*
     &    sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) + 64*rq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1)
     &     - 64*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*t1*
     &    sz**(-1)*s2**(-1) + 64*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg
     &    *msb1**2*u1*sz**(-1)*s2**(-1) + 32*rq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) - 32*rq*h2*ssz*lambda2*
     &    sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 64*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) - 64*rq*
     &    h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1)
     &     + 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*t1*
     &    sz**(-1)*s2**(-1) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 64*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*u1*sz**(-1)*s2**(-1)
     &     - 32*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s2**(-1) + 32*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1
     &    *sz**(-1)*s2**(-1) - 64*rq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 64*rq*h2*ssz*lambda2
     &    *sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) + 16*hl*hr*
     &    yuk1*lambda1*sb*cb*mg*mt*s*tx**(-1)*s1**(-1) + 16*hl*hr*yuk2*
     &    lambda2*sb*cb*mg*mt*s*tx**(-1)*s2**(-1) + 32*h1*yuk1*
     &    lambda1**2*sb*cb*sqrt2**(-1)*mg*s*s1**(-2) + 32*h1*yuk2*
     &    lambda1*lambda2*sb*cb*sqrt2**(-1)*mg*s*s1**(-1)*s2**(-1) + 32
     &    *h2*yuk1*lambda1*lambda2*sb*cb*sqrt2**(-1)*mg*s*s1**(-1)*
     &    s2**(-1) + 32*h2*yuk2*lambda2**2*sb*cb*sqrt2**(-1)*mg*s*
     &    s2**(-2) + 128*ssz**2*lq2*cb**4*m1**2*mg**2*msb1**2*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*lq2
     &    *cb**4*m1**2*mg**4*sz**(-2) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*lq2*
     &    cb**4*m1**2*msb1**4*sz**(-2) - 128*ssz**2*lq2*cb**4*mg**2*
     &    msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*t1
     &    *u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**4*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*lq2*cb**4*msb1**4*s**(-1)*t1*u1*sz**(-2)
     &     - 256*ssz**2*rq2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) + 128*
     &    ssz**2*rq2*cb**2*m1**2*mg**2*s*sz**(-2) + 128*ssz**2*rq2*
     &    cb**2*m1**2*mg**4*sz**(-2) + 128*ssz**2*rq2*cb**2*m1**2*
     &    msb1**4*sz**(-2) + 256*ssz**2*rq2*cb**2*mg**2*msb1**2*s**(-1)
     &    *t1*u1*sz**(-2) - 128*ssz**2*rq2*cb**2*mg**2*t1*u1*sz**(-2)
     &     - 128*ssz**2*rq2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2) - 128*
     &    ssz**2*rq2*cb**2*msb1**4*s**(-1)*t1*u1*sz**(-2) + 128*ssz**2*
     &    rq2*cb**4*m1**2*mg**2*msb1**2*sz**(-2) - 64*ssz**2*rq2*cb**4*
     &    m1**2*mg**2*s*sz**(-2) - 64*ssz**2*rq2*cb**4*m1**2*mg**4*
     &    sz**(-2) - 64*ssz**2*rq2*cb**4*m1**2*msb1**4*sz**(-2) - 128*
     &    ssz**2*rq2*cb**4*mg**2*msb1**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * ( 64*ssz**2*rq2*
     &    cb**4*mg**2*t1*u1*sz**(-2) + 64*ssz**2*rq2*cb**4*mg**4*
     &    s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*rq2*cb**4*msb1**4*s**(-1)*
     &    t1*u1*sz**(-2) + 128*ssz**2*rq2*m1**2*mg**2*msb1**2*sz**(-2)
     &     - 64*ssz**2*rq2*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*rq2*m1**2
     &    *mg**4*sz**(-2) - 64*ssz**2*rq2*m1**2*msb1**4*sz**(-2) - 128*
     &    ssz**2*rq2*mg**2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*
     &    rq2*mg**2*t1*u1*sz**(-2) + 64*ssz**2*rq2*mg**4*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*rq2*msb1**4*s**(-1)*t1*u1*sz**(-2) + 128
     &    *ssp**2*pq2*m1**2*mg**2*msb1**2*s**(-2) - 64*ssp**2*pq2*m1**2
     &    *mg**2*s**(-1) - 64*ssp**2*pq2*m1**2*mg**4*s**(-2) - 64*
     &    ssp**2*pq2*m1**2*msb1**4*s**(-2) - 128*ssp**2*pq2*mg**2*
     &    msb1**2*s**(-3)*t1*u1 + 64*ssp**2*pq2*mg**2*s**(-2)*t1*u1 + 
     &    64*ssp**2*pq2*mg**4*s**(-3)*t1*u1 + 64*ssp**2*pq2*msb1**4*
     &    s**(-3)*t1*u1 - 8*hss(1,1)**2*pq*ssp*sb**2 - 8*hss(1,1)**2*pq
     &    *ssp*cb**2 )
      MMs = MMs + SCC(1,2)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)**2*
     &    lq*ssz*cb**2*s*sz**(-1) - 8*hss(1,1)**2*rq*ssz*sb**2*s*
     &    sz**(-1) - 2*hss(1,1)**2*hl**2*cb**2*s*tx**(-1) - 2*hss(1,1)
     &    **2*hr**2*sb**2*s*tx**(-1) - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2
     &     - 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2 - 8*hss(1,2)*hss(2,1)*lq*
     &    ssz*cb**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*s*
     &    sz**(-1) - 2*hss(1,2)*hss(2,1)*hl**2*cb**2*s*tx**(-1) - 2*
     &    hss(1,2)*hss(2,1)*hr**2*sb**2*s*tx**(-1) - 8*hhss(1,1)*hl*hr*
     &    sb*cb*mg*mt*s*tx**(-1) - 16*hhss(1,1)*h1*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*s*s1**(-1) - 16*hhss(1,1)*h2*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**2*msb2**2*sz**(-2) + 128*lq*rq*ssz**2*cb**2*
     &    m1**2*mg**2*s*sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*mg**4*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*msb1**2*msb2**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*msb1**2*s**(-1)*t1*u1
     &    *sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*msb2**2*s**(-1)*t1*
     &    u1*sz**(-2) - 128*lq*rq*ssz**2*cb**2*mg**2*t1*u1*sz**(-2) - 
     &    128*lq*rq*ssz**2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2) - 128*lq*
     &    rq*ssz**2*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 128*
     &    lq*rq*ssz**2*cb**4*m1**2*mg**2*msb1**2*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**4*m1**2*mg**2*msb2**2*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**4*m1**2*mg**2*s*sz**(-2) - 128*lq*rq*ssz**2*cb**4*m1**2*
     &    mg**4*sz**(-2) - 128*lq*rq*ssz**2*cb**4*m1**2*msb1**2*msb2**2
     &    *sz**(-2) - 128*lq*rq*ssz**2*cb**4*mg**2*msb1**2*s**(-1)*t1*
     &    u1*sz**(-2) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq
     &    *ssz**2*cb**4*mg**2*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*msb1**2
     &    *msb2**2*s**(-1)*t1*u1*sz**(-2) + 8*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb1**2*s*tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb1**2*u1*tx**(-1)*sz**(-1) - 8*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb2**2*s*tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb2**2*t1*tx**(-1)*sz**(-1) + 8*lq*hl*hr*ssz*sb*cb*mg*mt*s*
     &    t1*tx**(-1)*sz**(-1) - 8*lq*hl*hr*ssz*sb*cb*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*s*
     &    tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*u1
     &    *tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*s
     &    *tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 32*lq*hl*hr*ssz*
     &    sb*cb**3*mg*mt*msb2**2*t1*tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz
     &    *sb*cb**3*mg*mt*s*t1*tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*
     &    cb**3*mg*mt*s*u1*tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3
     &    *mg**3*mt*t1*tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*
     &    mg**3*mt*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2*
     &    mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*m1**2
     &    *mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*
     &    m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2
     &    *mg**4*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2*
     &    msb1**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*
     &    mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*
     &    ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hl**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*lq*
     &    hl**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq
     &    *hl**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*lq*hl**2*
     &    ssz*cb**4*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) - 16*lq*hl**2
     &    *ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*
     &    hl**2*ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hl**2
     &    *ssz*cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*
     &    cb**4*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*
     &    ssz*cb**4*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hl**2*ssz*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) - 16*lq*hl**2*ssz*cb**4*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**2
     &    *msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*
     &    mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**2*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*
     &    ssz*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hr**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) - 16*lq*hr**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*mg**2
     &    *msb2**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*
     &    mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*msb1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*mg**2*
     &    msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*
     &    cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*
     &    hr**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz
     &    *cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*
     &    u1*sz**(-1)*s1**(-1) - 16*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)
     &    *mg*msb2**2*s*sz**(-1)*s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s1**(-1) + 16*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*s1**(-1) - 16*lq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 
     &    32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*
     &    sz**(-1)*s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*msb1**2*s*sz**(-1)*s1**(-1) - 64*lq*h1*ssz*lambda1*sb*
     &    cb**3*sqrt2**(-1)*mg*msb1**2*u1*sz**(-1)*s1**(-1) + 32*lq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 64*lq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s1**(-1)
     &     - 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*s*u1
     &    *sz**(-1)*s1**(-1) - 64*lq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1
     &    *sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1) + 16*lq*h2*
     &    ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*s2**(-1)
     &     + 32*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*u1*
     &    sz**(-1)*s2**(-1) - 16*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg
     &    *msb2**2*s*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s2**(-1) + 16*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) - 16*lq*
     &    h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 
     &    32*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 32*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) - 32*lq*
     &    h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*u1*sz**(-1)*s2**(-1) + 32*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s2**(-1) + 64*lq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s2**(-1)
     &     - 32*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s2**(-1) + 32*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1
     &    *sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 64*lq*h2*ssz*lambda2
     &    *sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) - 8*rq*hl*hr
     &    *ssz*sb*cb*mg*mt*msb1**2*s*tx**(-1)*sz**(-1) - 16*rq*hl*hr*
     &    ssz*sb*cb*mg*mt*msb1**2*u1*tx**(-1)*sz**(-1) + 8*rq*hl*hr*ssz
     &    *sb*cb*mg*mt*msb2**2*s*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb
     &    *cb*mg*mt*msb2**2*t1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*rq*hl*hr*ssz
     &    *sb*cb*mg*mt*s*t1*tx**(-1)*sz**(-1) + 8*rq*hl*hr*ssz*sb*cb*mg
     &    *mt*s*u1*tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb*mg**3*mt*
     &    t1*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*s*
     &    tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*u1
     &    *tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*s
     &    *tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*
     &    t1*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*t1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*mg**2*msb1**2
     &    *tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*mg**2*
     &    msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*m1**2*mg**2
     &    *s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hl**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*
     &    m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hl**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hl**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    rq*hl**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16
     &    *rq*hl**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*msb1**2*tx**(-1)
     &    *sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*msb2**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*s*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*mg**4*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*msb1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*mg**2*
     &    msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hl**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4
     &    *msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**2*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*rq
     &    *hr**2*ssz*cb**2*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hr**2*
     &    ssz*cb**2*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hr**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hr**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hr**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*
     &    ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hr**2
     &    *ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*
     &    m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    cb**4*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hr**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16
     &    *rq*hr**2*ssz*cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*
     &    s*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*
     &    mg*msb1**2*u1*sz**(-1)*s1**(-1) + 16*rq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s1**(-1) - 
     &    16*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) - 32*rq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1)
     &     + 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*
     &    s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*s*sz**(-1)*s1**(-1) + 64*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*msb1**2*u1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s1**(-1)
     &     - 64*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*
     &    sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 64*rq*h1*ssz*lambda1*
     &    sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) - 64*rq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1)
     &     - 16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*s*
     &    sz**(-1)*s2**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*u1*sz**(-1)*s2**(-1) + 
     &    16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*
     &    s2**(-1) + 32*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb2**2*
     &    t1*sz**(-1)*s2**(-1) - 16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s2**(-1) + 16*rq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) - 32*rq*h2*ssz*lambda2*
     &    sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) + 32*rq*
     &    h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s2**(-1) + 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*u1*sz**(-1)*s2**(-1) - 32*rq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s2**(-1) - 64*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s2**(-1)
     &     + 32*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s2**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 64*
     &    rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s2**(-1) - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*u1
     &    *sz**(-1)*s2**(-1) + 64*ssz**2*lq2*cb**2*m1**2*mg**2*msb1**2*
     &    sz**(-2) + 64*ssz**2*lq2*cb**2*m1**2*mg**2*msb2**2*sz**(-2)
     &     - 64*ssz**2*lq2*cb**2*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*lq2
     &    *cb**2*m1**2*mg**4*sz**(-2) - 64*ssz**2*lq2*cb**2*m1**2*
     &    msb1**2*msb2**2*sz**(-2) - 64*ssz**2*lq2*cb**2*mg**2*msb1**2*
     &    s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**2*mg**2*msb2**2*
     &    s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**2*mg**2*t1*u1*
     &    sz**(-2) + 64*ssz**2*lq2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*lq2*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*m1**2*mg**2*msb1**2*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*m1**2*mg**2*msb2**2*sz**(-2) + 64*ssz**2*lq2
     &    *cb**4*m1**2*mg**2*s*sz**(-2) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 64*ssz**2*lq2*
     &    cb**4*m1**2*mg**4*sz**(-2) + 64*ssz**2*lq2*cb**4*m1**2*
     &    msb1**2*msb2**2*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*msb1**2*
     &    s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*msb2**2*
     &    s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**4*mg**2*t1*u1*
     &    sz**(-2) - 64*ssz**2*lq2*cb**4*mg**4*s**(-1)*t1*u1*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*rq2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) + 64*
     &    ssz**2*rq2*cb**2*m1**2*mg**2*msb2**2*sz**(-2) - 64*ssz**2*rq2
     &    *cb**2*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*
     &    mg**4*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*msb1**2*msb2**2*
     &    sz**(-2) - 64*ssz**2*rq2*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*
     &    sz**(-2) - 64*ssz**2*rq2*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*rq2*cb**2*mg**2*t1*u1*sz**(-2) + 64*
     &    ssz**2*rq2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*rq2
     &    *cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*rq2*
     &    cb**4*m1**2*mg**2*msb1**2*sz**(-2) - 64*ssz**2*rq2*cb**4*
     &    m1**2*mg**2*msb2**2*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*
     &    mg**2*s*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*mg**4*sz**(-2)
     &     + 64*ssz**2*rq2*cb**4*m1**2*msb1**2*msb2**2*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*mg**2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*rq2*cb**4*mg**2*t1*u1*sz**(-2) - 64*ssz**2*rq2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*rq2*cb**4*msb1**2*
     &    msb2**2*s**(-1)*t1*u1*sz**(-2) + 8*hss(1,1)*hss(2,1)*lq*ssz*
     &    sb*cb*s*sz**(-1) - 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*s*
     &    sz**(-1) + 2*hss(1,1)*hss(2,1)*hl**2*sb*cb*s*tx**(-1) - 2*
     &    hss(1,1)*hss(2,1)*hr**2*sb*cb*s*tx**(-1) + 8*hss(2,1)*hss(2,2
     &    )*lq*ssz*sb*cb*s*sz**(-1) - 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*
     &    s*sz**(-1) + 2*hss(2,1)*hss(2,2)*hl**2*sb*cb*s*tx**(-1) - 2*
     &    hss(2,1)*hss(2,2)*hr**2*sb*cb*s*tx**(-1) )
      MMs = MMs + SCC(1,3)*Nc*Cf*Pi*alphas*prefac * ( 4*hhss(1,1)*hl*hr
     &    *sb**2*mg*mt*s*tx**(-1) - 4*hhss(1,1)*hl*hr*cb**2*mg*mt*s*
     &    tx**(-1) + 8*hhss(1,1)*h1*lambda1*sb**2*sqrt2**(-1)*mg*s*
     &    s1**(-1) - 8*hhss(1,1)*h1*lambda1*cb**2*sqrt2**(-1)*mg*s*
     &    s1**(-1) + 8*hhss(1,1)*h2*lambda2*sb**2*sqrt2**(-1)*mg*s*
     &    s2**(-1) - 8*hhss(1,1)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*
     &    s2**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**2*msb2**2*sz**(-2) + 128*lq*rq*ssz**2*cb**2*
     &    m1**2*mg**2*s*sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*mg**4*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*m1**2*msb1**2*msb2**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*msb1**2*s**(-1)*t1*u1
     &    *sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*msb2**2*s**(-1)*t1*
     &    u1*sz**(-2) - 128*lq*rq*ssz**2*cb**2*mg**2*t1*u1*sz**(-2) - 
     &    128*lq*rq*ssz**2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2) - 128*lq*
     &    rq*ssz**2*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 128*
     &    lq*rq*ssz**2*cb**4*m1**2*mg**2*msb1**2*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**4*m1**2*mg**2*msb2**2*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**4*m1**2*mg**2*s*sz**(-2) - 128*lq*rq*ssz**2*cb**4*m1**2*
     &    mg**4*sz**(-2) - 128*lq*rq*ssz**2*cb**4*m1**2*msb1**2*msb2**2
     &    *sz**(-2) - 128*lq*rq*ssz**2*cb**4*mg**2*msb1**2*s**(-1)*t1*
     &    u1*sz**(-2) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq
     &    *ssz**2*cb**4*mg**2*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**4*msb1**2
     &    *msb2**2*s**(-1)*t1*u1*sz**(-2) - 8*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb1**2*s*tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb1**2*t1*tx**(-1)*sz**(-1) + 8*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb2**2*s*tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb*mg*mt*
     &    msb2**2*u1*tx**(-1)*sz**(-1) + 8*lq*hl*hr*ssz*sb*cb*mg*mt*s*
     &    t1*tx**(-1)*sz**(-1) - 8*lq*hl*hr*ssz*sb*cb*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*s*
     &    tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*t1
     &    *tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*s
     &    *tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 32*lq*hl*hr*
     &    ssz*sb*cb**3*mg*mt*msb2**2*u1*tx**(-1)*sz**(-1) - 16*lq*hl*hr
     &    *ssz*sb*cb**3*mg*mt*s*t1*tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*
     &    sb*cb**3*mg*mt*s*u1*tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*
     &    cb**3*mg**3*mt*t1*tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*
     &    cb**3*mg**3*mt*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2*
     &    m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**2
     &    *m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*
     &    cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2
     &    *m1**2*mg**4*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*m1**2*
     &    msb1**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**2*
     &    mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hl**2*
     &    ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hl**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*lq*
     &    hl**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq
     &    *hl**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 16*lq*hl**2*
     &    ssz*cb**4*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) - 16*lq*hl**2
     &    *ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*
     &    hl**2*ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hl**2
     &    *ssz*cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*
     &    cb**4*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hl**2*
     &    ssz*cb**4*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hl**2*ssz*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hl**2*ssz*cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) - 16*lq*hl**2*ssz*cb**4*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**2
     &    *msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*
     &    mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**2*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*
     &    ssz*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*
     &    lq*hr**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) - 16*lq*hr**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*m1**2*mg**2
     &    *msb2**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*
     &    mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*m1**2*msb1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*mg**2*
     &    msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*
     &    cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*
     &    hr**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz
     &    *cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*
     &    t1*sz**(-1)*s1**(-1) + 16*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)
     &    *mg*msb2**2*s*sz**(-1)*s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s1**(-1) + 16*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*s1**(-1) - 16*lq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 
     &    32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*
     &    sz**(-1)*s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*msb1**2*s*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1*sb*
     &    cb**3*sqrt2**(-1)*mg*msb1**2*t1*sz**(-1)*s1**(-1) - 32*lq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 64*lq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s1**(-1)
     &     - 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*s*u1
     &    *sz**(-1)*s1**(-1) - 64*lq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1
     &    *sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1) - 16*lq*h2*
     &    ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*s2**(-1)
     &     - 32*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*t1*
     &    sz**(-1)*s2**(-1) + 16*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg
     &    *msb2**2*s*sz**(-1)*s2**(-1) + 32*lq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s2**(-1) + 16*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) - 16*lq*
     &    h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 
     &    32*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 32*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) + 32*lq*
     &    h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s2**(-1) + 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*t1*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s2**(-1) - 64*lq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s2**(-1)
     &     - 32*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s2**(-1) + 32*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1
     &    *sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 64*lq*h2*ssz*lambda2
     &    *sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) + 8*rq*hl*hr
     &    *ssz*sb*cb*mg*mt*msb1**2*s*tx**(-1)*sz**(-1) + 16*rq*hl*hr*
     &    ssz*sb*cb*mg*mt*msb1**2*t1*tx**(-1)*sz**(-1) - 8*rq*hl*hr*ssz
     &    *sb*cb*mg*mt*msb2**2*s*tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb
     &    *cb*mg*mt*msb2**2*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 8*rq*hl*hr*ssz
     &    *sb*cb*mg*mt*s*t1*tx**(-1)*sz**(-1) + 8*rq*hl*hr*ssz*sb*cb*mg
     &    *mt*s*u1*tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb*mg**3*mt*
     &    t1*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*s*
     &    tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb1**2*t1
     &    *tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*s
     &    *tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*
     &    u1*tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*t1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*mg**2*msb1**2
     &    *tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*m1**2*mg**2*
     &    msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*m1**2*mg**2
     &    *s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hl**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*
     &    m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hl**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hl**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    rq*hl**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16
     &    *rq*hl**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*msb1**2*tx**(-1)
     &    *sz**(-1) + 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*msb2**2*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*mg**2*s*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*mg**4*
     &    tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*m1**2*msb1**2*
     &    msb2**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*mg**2*
     &    msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*hl**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4
     &    *msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**2*m1**2*mg**2*msb1**2*tx**(-1)*sz**(-1) + 16*rq
     &    *hr**2*ssz*cb**2*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hr**2*
     &    ssz*cb**2*m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hr**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hr**2*ssz*cb**2*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*hr**2*ssz*cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)
     &    *sz**(-1) + 16*rq*hr**2*ssz*cb**2*msb1**2*msb2**2*s**(-1)*t1*
     &    u1*tx**(-1)*sz**(-1) - 16*rq*hr**2*ssz*cb**4*m1**2*mg**2*
     &    msb1**2*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hr**2*
     &    ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hr**2
     &    *ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*cb**4*
     &    m1**2*msb1**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hr**2*ssz*
     &    cb**4*mg**2*msb1**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1)
     &     - 16*rq*hr**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*
     &    rq*hr**2*ssz*cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16
     &    *rq*hr**2*ssz*cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb1**2*
     &    s*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*
     &    mg*msb1**2*t1*sz**(-1)*s1**(-1) - 16*rq*h1*ssz*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s1**(-1) - 
     &    16*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 16*rq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) - 32*rq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1)
     &     + 32*rq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*
     &    s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*s*sz**(-1)*s1**(-1) - 64*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*msb1**2*t1*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s1**(-1)
     &     + 64*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*u1*
     &    sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 64*rq*h1*ssz*lambda1*
     &    sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) - 64*rq*h1*
     &    ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1)
     &     + 16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*s*
     &    sz**(-1)*s2**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*msb1**2*t1*sz**(-1)*s2**(-1) - 
     &    16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*
     &    s2**(-1) - 32*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*msb2**2*
     &    u1*sz**(-1)*s2**(-1) - 16*rq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)
     &    *mg*s*t1*sz**(-1)*s2**(-1) + 16*rq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) - 32*rq*h2*ssz*lambda2*
     &    sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 32*rq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) - 32*rq*
     &    h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb1**2*s*sz**(-1)*
     &    s2**(-1) - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb1**2*t1*sz**(-1)*s2**(-1) + 32*rq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*msb2**2*s*sz**(-1)*s2**(-1) + 64*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s2**(-1)
     &     + 32*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)
     &    *s2**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 32*rq*h2*ssz*
     &    lambda2*sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 64*
     &    rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s2**(-1) - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*u1
     &    *sz**(-1)*s2**(-1) + 64*ssz**2*lq2*cb**2*m1**2*mg**2*msb1**2*
     &    sz**(-2) + 64*ssz**2*lq2*cb**2*m1**2*mg**2*msb2**2*sz**(-2)
     &     - 64*ssz**2*lq2*cb**2*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*lq2
     &    *cb**2*m1**2*mg**4*sz**(-2) - 64*ssz**2*lq2*cb**2*m1**2*
     &    msb1**2*msb2**2*sz**(-2) - 64*ssz**2*lq2*cb**2*mg**2*msb1**2*
     &    s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**2*mg**2*msb2**2*
     &    s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**2*mg**2*t1*u1*
     &    sz**(-2) + 64*ssz**2*lq2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*lq2*cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*m1**2*mg**2*msb1**2*sz**(-2) - 64*
     &    ssz**2*lq2*cb**4*m1**2*mg**2*msb2**2*sz**(-2) + 64*ssz**2*lq2
     &    *cb**4*m1**2*mg**2*s*sz**(-2) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 64*ssz**2*lq2*
     &    cb**4*m1**2*mg**4*sz**(-2) + 64*ssz**2*lq2*cb**4*m1**2*
     &    msb1**2*msb2**2*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*msb1**2*
     &    s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*msb2**2*
     &    s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*lq2*cb**4*mg**2*t1*u1*
     &    sz**(-2) - 64*ssz**2*lq2*cb**4*mg**4*s**(-1)*t1*u1*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2)
     &     + 64*ssz**2*rq2*cb**2*m1**2*mg**2*msb1**2*sz**(-2) + 64*
     &    ssz**2*rq2*cb**2*m1**2*mg**2*msb2**2*sz**(-2) - 64*ssz**2*rq2
     &    *cb**2*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*
     &    mg**4*sz**(-2) - 64*ssz**2*rq2*cb**2*m1**2*msb1**2*msb2**2*
     &    sz**(-2) - 64*ssz**2*rq2*cb**2*mg**2*msb1**2*s**(-1)*t1*u1*
     &    sz**(-2) - 64*ssz**2*rq2*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*rq2*cb**2*mg**2*t1*u1*sz**(-2) + 64*
     &    ssz**2*rq2*cb**2*mg**4*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*rq2
     &    *cb**2*msb1**2*msb2**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*rq2*
     &    cb**4*m1**2*mg**2*msb1**2*sz**(-2) - 64*ssz**2*rq2*cb**4*
     &    m1**2*mg**2*msb2**2*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*
     &    mg**2*s*sz**(-2) + 64*ssz**2*rq2*cb**4*m1**2*mg**4*sz**(-2)
     &     + 64*ssz**2*rq2*cb**4*m1**2*msb1**2*msb2**2*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*mg**2*msb1**2*s**(-1)*t1*u1*sz**(-2) + 64*
     &    ssz**2*rq2*cb**4*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) - 64*
     &    ssz**2*rq2*cb**4*mg**2*t1*u1*sz**(-2) - 64*ssz**2*rq2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) - 64*ssz**2*rq2*cb**4*msb1**2*
     &    msb2**2*s**(-1)*t1*u1*sz**(-2) + 8*hss(1,1)*hss(1,2)*lq*ssz*
     &    sb*cb*s*sz**(-1) - 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*s*
     &    sz**(-1) + 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*s*tx**(-1) - 2*
     &    hss(1,1)*hss(1,2)*hr**2*sb*cb*s*tx**(-1) + 8*hss(1,2)*hss(2,2
     &    )*lq*ssz*sb*cb*s*sz**(-1) - 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*
     &    s*sz**(-1) + 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*s*tx**(-1) - 2*
     &    hss(1,2)*hss(2,2)*hr**2*sb*cb*s*tx**(-1) )
      MMs = MMs + SCC(1,4)*Nc*Cf*Pi*alphas*prefac * ( 4*hhss(1,2)*hl*hr
     &    *sb**2*mg*mt*s*tx**(-1) - 4*hhss(1,2)*hl*hr*cb**2*mg*mt*s*
     &    tx**(-1) + 8*hhss(1,2)*h1*lambda1*sb**2*sqrt2**(-1)*mg*s*
     &    s1**(-1) - 8*hhss(1,2)*h1*lambda1*cb**2*sqrt2**(-1)*mg*s*
     &    s1**(-1) + 8*hhss(1,2)*h2*lambda2*sb**2*sqrt2**(-1)*mg*s*
     &    s2**(-1) - 8*hhss(1,2)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*
     &    s2**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 256*pq*lq*ssz*
     &    ssp*cb**2*m1**2*mg**2*msb2**2*s**(-1)*sz**(-1) + 128*pq*lq*
     &    ssz*ssp*cb**2*m1**2*mg**2*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*
     &    m1**2*mg**4*s**(-1)*sz**(-1) + 128*pq*lq*ssz*ssp*cb**2*m1**2*
     &    msb2**4*s**(-1)*sz**(-1) + 256*pq*lq*ssz*ssp*cb**2*mg**2*
     &    msb2**2*s**(-2)*t1*u1*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*
     &    mg**2*s**(-1)*t1*u1*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*mg**4*
     &    s**(-2)*t1*u1*sz**(-1) - 128*pq*lq*ssz*ssp*cb**2*msb2**4*
     &    s**(-2)*t1*u1*sz**(-1) + 256*pq*lq*ssz*ssp*m1**2*mg**2*
     &    msb2**2*s**(-1)*sz**(-1) - 128*pq*lq*ssz*ssp*m1**2*mg**2*
     &    sz**(-1) - 128*pq*lq*ssz*ssp*m1**2*mg**4*s**(-1)*sz**(-1) - 
     &    128*pq*lq*ssz*ssp*m1**2*msb2**4*s**(-1)*sz**(-1) - 256*pq*lq*
     &    ssz*ssp*mg**2*msb2**2*s**(-2)*t1*u1*sz**(-1) + 128*pq*lq*ssz*
     &    ssp*mg**2*s**(-1)*t1*u1*sz**(-1) + 128*pq*lq*ssz*ssp*mg**4*
     &    s**(-2)*t1*u1*sz**(-1) + 128*pq*lq*ssz*ssp*msb2**4*s**(-2)*t1
     &    *u1*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 256*pq*rq*ssz*ssp
     &    *cb**2*m1**2*mg**2*msb2**2*s**(-1)*sz**(-1) - 128*pq*rq*ssz*
     &    ssp*cb**2*m1**2*mg**2*sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*
     &    m1**2*mg**4*s**(-1)*sz**(-1) - 128*pq*rq*ssz*ssp*cb**2*m1**2*
     &    msb2**4*s**(-1)*sz**(-1) - 256*pq*rq*ssz*ssp*cb**2*mg**2*
     &    msb2**2*s**(-2)*t1*u1*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*
     &    mg**2*s**(-1)*t1*u1*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*mg**4*
     &    s**(-2)*t1*u1*sz**(-1) + 128*pq*rq*ssz*ssp*cb**2*msb2**4*
     &    s**(-2)*t1*u1*sz**(-1) + 32*pq*hl*hr*ssp*sb*cb*mg*mt*msb2**2*
     &    s**(-1)*t1*tx**(-1) - 32*pq*hl*hr*ssp*sb*cb*mg*mt*msb2**2*
     &    s**(-1)*u1*tx**(-1) - 16*pq*hl*hr*ssp*sb*cb*mg*mt*t1*tx**(-1)
     &     + 16*pq*hl*hr*ssp*sb*cb*mg*mt*u1*tx**(-1) - 32*pq*hl*hr*ssp*
     &    sb*cb*mg**3*mt*s**(-1)*t1*tx**(-1) + 32*pq*hl*hr*ssp*sb*cb*
     &    mg**3*mt*s**(-1)*u1*tx**(-1) - 32*pq*hl**2*ssp*cb**2*m1**2*
     &    mg**2*msb2**2*s**(-1)*tx**(-1) + 16*pq*hl**2*ssp*cb**2*m1**2*
     &    mg**2*tx**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 16*pq*hl**2*ssp*
     &    cb**2*m1**2*mg**4*s**(-1)*tx**(-1) + 16*pq*hl**2*ssp*cb**2*
     &    m1**2*msb2**4*s**(-1)*tx**(-1) + 32*pq*hl**2*ssp*cb**2*mg**2*
     &    msb2**2*s**(-2)*t1*u1*tx**(-1) - 16*pq*hl**2*ssp*cb**2*mg**2*
     &    s**(-1)*t1*u1*tx**(-1) - 16*pq*hl**2*ssp*cb**2*mg**4*s**(-2)*
     &    t1*u1*tx**(-1) - 16*pq*hl**2*ssp*cb**2*msb2**4*s**(-2)*t1*u1*
     &    tx**(-1) + 32*pq*hl**2*ssp*m1**2*mg**2*msb2**2*s**(-1)*
     &    tx**(-1) - 16*pq*hl**2*ssp*m1**2*mg**2*tx**(-1) - 16*pq*hl**2
     &    *ssp*m1**2*mg**4*s**(-1)*tx**(-1) - 16*pq*hl**2*ssp*m1**2*
     &    msb2**4*s**(-1)*tx**(-1) - 32*pq*hl**2*ssp*mg**2*msb2**2*
     &    s**(-2)*t1*u1*tx**(-1) + 16*pq*hl**2*ssp*mg**2*s**(-1)*t1*u1*
     &    tx**(-1) + 16*pq*hl**2*ssp*mg**4*s**(-2)*t1*u1*tx**(-1) + 16*
     &    pq*hl**2*ssp*msb2**4*s**(-2)*t1*u1*tx**(-1) + 32*pq*hr**2*ssp
     &    *cb**2*m1**2*mg**2*msb2**2*s**(-1)*tx**(-1) - 16*pq*hr**2*ssp
     &    *cb**2*m1**2*mg**2*tx**(-1) - 16*pq*hr**2*ssp*cb**2*m1**2*
     &    mg**4*s**(-1)*tx**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 16*pq*hr**2*
     &    ssp*cb**2*m1**2*msb2**4*s**(-1)*tx**(-1) - 32*pq*hr**2*ssp*
     &    cb**2*mg**2*msb2**2*s**(-2)*t1*u1*tx**(-1) + 16*pq*hr**2*ssp*
     &    cb**2*mg**2*s**(-1)*t1*u1*tx**(-1) + 16*pq*hr**2*ssp*cb**2*
     &    mg**4*s**(-2)*t1*u1*tx**(-1) + 16*pq*hr**2*ssp*cb**2*msb2**4*
     &    s**(-2)*t1*u1*tx**(-1) + 64*pq*h1*ssp*lambda1*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*s**(-1)*t1*s1**(-1) - 64*pq*h1*ssp*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*msb2**2*s**(-1)*u1*s1**(-1) - 32
     &    *pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*t1*s1**(-1) + 32*pq*
     &    h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg*u1*s1**(-1) - 64*pq*h1*
     &    ssp*lambda1*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*t1*s1**(-1) + 64*
     &    pq*h1*ssp*lambda1*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*u1*s1**(-1)
     &     + 64*pq*h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg*msb2**2*s**(-1)*
     &    t1*s2**(-1) - 64*pq*h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg*
     &    msb2**2*s**(-1)*u1*s2**(-1) - 32*pq*h2*ssp*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*t1*s2**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 32*pq*h2*ssp*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*u1*s2**(-1) - 64*pq*h2*ssp*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*t1*s2**(-1) + 64*pq*
     &    h2*ssp*lambda2*sb*cb*sqrt2**(-1)*mg**3*s**(-1)*u1*s2**(-1) + 
     &    256*lq*rq*ssz**2*cb**2*m1**2*mg**2*msb2**2*sz**(-2) - 128*lq*
     &    rq*ssz**2*cb**2*m1**2*mg**2*s*sz**(-2) - 128*lq*rq*ssz**2*
     &    cb**2*m1**2*mg**4*sz**(-2) - 128*lq*rq*ssz**2*cb**2*m1**2*
     &    msb2**4*sz**(-2) - 256*lq*rq*ssz**2*cb**2*mg**2*msb2**2*
     &    s**(-1)*t1*u1*sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**2*t1*u1*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*mg**4*s**(-1)*t1*u1*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**2*msb2**4*s**(-1)*t1*u1*
     &    sz**(-2) - 256*lq*rq*ssz**2*cb**4*m1**2*mg**2*msb2**2*
     &    sz**(-2) + 128*lq*rq*ssz**2*cb**4*m1**2*mg**2*s*sz**(-2) + 
     &    128*lq*rq*ssz**2*cb**4*m1**2*mg**4*sz**(-2) + 128*lq*rq*
     &    ssz**2*cb**4*m1**2*msb2**4*sz**(-2) + 256*lq*rq*ssz**2*cb**4*
     &    mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 128*lq*rq*
     &    ssz**2*cb**4*mg**2*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) - 128*lq*rq*ssz**2*cb**4*msb2**4
     &    *s**(-1)*t1*u1*sz**(-2) + 32*lq*hl*hr*ssz*sb*cb*mg*mt*msb2**2
     &    *t1*tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*cb*mg*mt*msb2**2*
     &    u1*tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb*mg*mt*s*t1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb*mg*mt*s*u1*tx**(-1)
     &    *sz**(-1) - 32*lq*hl*hr*ssz*sb*cb*mg**3*mt*t1*tx**(-1)*
     &    sz**(-1) + 32*lq*hl*hr*ssz*sb*cb*mg**3*mt*u1*tx**(-1)*
     &    sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*t1*tx**(-1)
     &    *sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*s*t1*
     &    tx**(-1)*sz**(-1) - 16*lq*hl*hr*ssz*sb*cb**3*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) + 32*lq*hl*hr*ssz*sb*cb**3*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) - 32*lq*hl*hr*ssz*sb*cb**3*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 64*lq*hl**2*
     &    ssz*cb**2*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 32*lq*hl**2
     &    *ssz*cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*cb**2*
     &    m1**2*msb2**4*tx**(-1)*sz**(-1) + 64*lq*hl**2*ssz*cb**2*mg**2
     &    *msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*
     &    cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**2*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**2
     &    *msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*lq*hl**2*ssz*
     &    cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz
     &    *cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*cb**4*
     &    m1**2*msb2**4*tx**(-1)*sz**(-1) - 32*lq*hl**2*ssz*cb**4*mg**2
     &    *msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hl**2*ssz*
     &    cb**4*msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*lq*hl**2*
     &    ssz*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*
     &    m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**2*mg**4
     &    *tx**(-1)*sz**(-1) - 16*lq*hl**2*ssz*m1**2*msb2**4*tx**(-1)*
     &    sz**(-1) - 32*lq*hl**2*ssz*mg**2*msb2**2*s**(-1)*t1*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hl**2*ssz*mg**2*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*mg**4*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 16*lq*hl**2*ssz*msb2**4*s**(-1)*t1*u1*tx**(-1)*
     &    sz**(-1) + 32*lq*hr**2*ssz*cb**2*m1**2*mg**2*msb2**2*tx**(-1)
     &    *sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**2*s*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*mg**4*tx**(-1)*
     &    sz**(-1) - 16*lq*hr**2*ssz*cb**2*m1**2*msb2**4*tx**(-1)*
     &    sz**(-1) - 32*lq*hr**2*ssz*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*
     &    tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**2*mg**2*t1*u1*
     &    tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 16*lq*hr**2*ssz*
     &    cb**2*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz
     &    *cb**2*msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*lq*hr**2*
     &    ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*lq*hr**2
     &    *ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*lq*hr**2*ssz*cb**4*
     &    m1**2*msb2**4*tx**(-1)*sz**(-1) + 32*lq*hr**2*ssz*cb**4*mg**2
     &    *msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*lq*hr**2*ssz*cb**4
     &    *msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 64*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*s1**(-1) - 
     &    64*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)
     &    *s1**(-1) - 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg*s*t1*
     &    sz**(-1)*s1**(-1) + 32*lq*h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg
     &    *s*u1*sz**(-1)*s1**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 64*lq*h1*ssz*
     &    lambda1*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) + 64*lq*
     &    h1*ssz*lambda1*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1)
     &     - 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*
     &    sz**(-1)*s1**(-1) + 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)
     &    *mg*msb2**2*u1*sz**(-1)*s1**(-1) + 32*lq*h1*ssz*lambda1*sb*
     &    cb**3*sqrt2**(-1)*mg*s*t1*sz**(-1)*s1**(-1) - 32*lq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) + 64*
     &    lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s1**(-1) - 64*lq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1
     &    *sz**(-1)*s1**(-1) + 64*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*
     &    mg*msb2**2*t1*sz**(-1)*s2**(-1) - 64*lq*h2*ssz*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*msb2**2*u1*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) + 32*lq*
     &    h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) - 
     &    64*lq*h2*ssz*lambda2*sb*cb*sqrt2**(-1)*mg**3*t1*sz**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 64*lq*h2*ssz*
     &    lambda2*sb*cb*sqrt2**(-1)*mg**3*u1*sz**(-1)*s2**(-1) - 64*lq*
     &    h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*
     &    s2**(-1) + 64*lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb2**2*u1*sz**(-1)*s2**(-1) + 32*lq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) - 32*lq*h2*ssz*lambda2*
     &    sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) + 64*lq*h2*ssz
     &    *lambda2*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) - 64
     &    *lq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*
     &    s2**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*t1*tx**(-1)
     &    *sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg*mt*msb2**2*u1*
     &    tx**(-1)*sz**(-1) - 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*t1*
     &    tx**(-1)*sz**(-1) + 16*rq*hl*hr*ssz*sb*cb**3*mg*mt*s*u1*
     &    tx**(-1)*sz**(-1) - 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*t1*
     &    tx**(-1)*sz**(-1) + 32*rq*hl*hr*ssz*sb*cb**3*mg**3*mt*u1*
     &    tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 32*rq*hl**2*ssz*
     &    cb**2*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz
     &    *cb**2*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**2*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**2*
     &    m1**2*msb2**4*tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*cb**2*mg**2
     &    *msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**2*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**2
     &    *msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 32*rq*hl**2*ssz*
     &    cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz
     &    *cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*
     &    cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) + 16*rq*hl**2*ssz*cb**4*
     &    m1**2*msb2**4*tx**(-1)*sz**(-1) + 32*rq*hl**2*ssz*cb**4*mg**2
     &    *msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*
     &    cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) - 16*rq*hl**2*ssz*cb**4*
     &    mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 16*rq*hl**2*
     &    ssz*cb**4*msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 32*rq*
     &    hr**2*ssz*cb**4*m1**2*mg**2*msb2**2*tx**(-1)*sz**(-1) - 16*rq
     &    *hr**2*ssz*cb**4*m1**2*mg**2*s*tx**(-1)*sz**(-1) - 16*rq*
     &    hr**2*ssz*cb**4*m1**2*mg**4*tx**(-1)*sz**(-1) - 16*rq*hr**2*
     &    ssz*cb**4*m1**2*msb2**4*tx**(-1)*sz**(-1) - 32*rq*hr**2*ssz*
     &    cb**4*mg**2*msb2**2*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*
     &    hr**2*ssz*cb**4*mg**2*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2*
     &    ssz*cb**4*mg**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 16*rq*hr**2
     &    *ssz*cb**4*msb2**4*s**(-1)*t1*u1*tx**(-1)*sz**(-1) + 64*rq*h1
     &    *ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)*
     &    s1**(-1) - 64*rq*h1*ssz*lambda1*sb*cb**3*sqrt2**(-1)*mg*
     &    msb2**2*u1*sz**(-1)*s1**(-1) - 32*rq*h1*ssz*lambda1*sb*cb**3*
     &    sqrt2**(-1)*mg*s*t1*sz**(-1)*s1**(-1) + 32*rq*h1*ssz*lambda1*
     &    sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s1**(-1) - 64*rq*h1*ssz
     &    *lambda1*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s1**(-1) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 64*rq*h1*ssz*
     &    lambda1*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*s1**(-1) + 64*
     &    rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*msb2**2*t1*sz**(-1)
     &    *s2**(-1) - 64*rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg*
     &    msb2**2*u1*sz**(-1)*s2**(-1) - 32*rq*h2*ssz*lambda2*sb*cb**3*
     &    sqrt2**(-1)*mg*s*t1*sz**(-1)*s2**(-1) + 32*rq*h2*ssz*lambda2*
     &    sb*cb**3*sqrt2**(-1)*mg*s*u1*sz**(-1)*s2**(-1) - 64*rq*h2*ssz
     &    *lambda2*sb*cb**3*sqrt2**(-1)*mg**3*t1*sz**(-1)*s2**(-1) + 64
     &    *rq*h2*ssz*lambda2*sb*cb**3*sqrt2**(-1)*mg**3*u1*sz**(-1)*
     &    s2**(-1) - 16*hl*hr*yuk1*lambda1*sb*cb*mg*mt*s*tx**(-1)*
     &    s1**(-1) - 16*hl*hr*yuk2*lambda2*sb*cb*mg*mt*s*tx**(-1)*
     &    s2**(-1) - 32*h1*yuk1*lambda1**2*sb*cb*sqrt2**(-1)*mg*s*
     &    s1**(-2) - 32*h1*yuk2*lambda1*lambda2*sb*cb*sqrt2**(-1)*mg*s*
     &    s1**(-1)*s2**(-1) - 32*h2*yuk1*lambda1*lambda2*sb*cb*
     &    sqrt2**(-1)*mg*s*s1**(-1)*s2**(-1) - 32*h2*yuk2*lambda2**2*sb
     &    *cb*sqrt2**(-1)*mg*s*s2**(-2) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 256*ssz**2*lq2
     &    *cb**2*m1**2*mg**2*msb2**2*sz**(-2) + 128*ssz**2*lq2*cb**2*
     &    m1**2*mg**2*s*sz**(-2) + 128*ssz**2*lq2*cb**2*m1**2*mg**4*
     &    sz**(-2) + 128*ssz**2*lq2*cb**2*m1**2*msb2**4*sz**(-2) + 256*
     &    ssz**2*lq2*cb**2*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) - 128*
     &    ssz**2*lq2*cb**2*mg**2*t1*u1*sz**(-2) - 128*ssz**2*lq2*cb**2*
     &    mg**4*s**(-1)*t1*u1*sz**(-2) - 128*ssz**2*lq2*cb**2*msb2**4*
     &    s**(-1)*t1*u1*sz**(-2) + 128*ssz**2*lq2*cb**4*m1**2*mg**2*
     &    msb2**2*sz**(-2) - 64*ssz**2*lq2*cb**4*m1**2*mg**2*s*sz**(-2)
     &     - 64*ssz**2*lq2*cb**4*m1**2*mg**4*sz**(-2) - 64*ssz**2*lq2*
     &    cb**4*m1**2*msb2**4*sz**(-2) - 128*ssz**2*lq2*cb**4*mg**2*
     &    msb2**2*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**2*t1
     &    *u1*sz**(-2) + 64*ssz**2*lq2*cb**4*mg**4*s**(-1)*t1*u1*
     &    sz**(-2) + 64*ssz**2*lq2*cb**4*msb2**4*s**(-1)*t1*u1*sz**(-2)
     &     + 128*ssz**2*lq2*m1**2*mg**2*msb2**2*sz**(-2) - 64*ssz**2*
     &    lq2*m1**2*mg**2*s*sz**(-2) )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * (  - 64*ssz**2*lq2*
     &    m1**2*mg**4*sz**(-2) - 64*ssz**2*lq2*m1**2*msb2**4*sz**(-2)
     &     - 128*ssz**2*lq2*mg**2*msb2**2*s**(-1)*t1*u1*sz**(-2) + 64*
     &    ssz**2*lq2*mg**2*t1*u1*sz**(-2) + 64*ssz**2*lq2*mg**4*s**(-1)
     &    *t1*u1*sz**(-2) + 64*ssz**2*lq2*msb2**4*s**(-1)*t1*u1*
     &    sz**(-2) + 128*ssz**2*rq2*cb**4*m1**2*mg**2*msb2**2*sz**(-2)
     &     - 64*ssz**2*rq2*cb**4*m1**2*mg**2*s*sz**(-2) - 64*ssz**2*rq2
     &    *cb**4*m1**2*mg**4*sz**(-2) - 64*ssz**2*rq2*cb**4*m1**2*
     &    msb2**4*sz**(-2) - 128*ssz**2*rq2*cb**4*mg**2*msb2**2*s**(-1)
     &    *t1*u1*sz**(-2) + 64*ssz**2*rq2*cb**4*mg**2*t1*u1*sz**(-2) + 
     &    64*ssz**2*rq2*cb**4*mg**4*s**(-1)*t1*u1*sz**(-2) + 64*ssz**2*
     &    rq2*cb**4*msb2**4*s**(-1)*t1*u1*sz**(-2) + 128*ssp**2*pq2*
     &    m1**2*mg**2*msb2**2*s**(-2) - 64*ssp**2*pq2*m1**2*mg**2*
     &    s**(-1) - 64*ssp**2*pq2*m1**2*mg**4*s**(-2) - 64*ssp**2*pq2*
     &    m1**2*msb2**4*s**(-2) - 128*ssp**2*pq2*mg**2*msb2**2*s**(-3)*
     &    t1*u1 )
      MMs = MMs + SCC(1,5)*Nc*Cf*Pi*alphas*prefac * ( 64*ssp**2*pq2*
     &    mg**2*s**(-2)*t1*u1 + 64*ssp**2*pq2*mg**4*s**(-3)*t1*u1 + 64*
     &    ssp**2*pq2*msb2**4*s**(-3)*t1*u1 - 8*hss(1,2)*hss(2,1)*pq*ssp
     &    *sb**2 - 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2 - 8*hss(1,2)*hss(2,
     &    1)*lq*ssz*sb**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2
     &    *s*sz**(-1) - 2*hss(1,2)*hss(2,1)*hl**2*sb**2*s*tx**(-1) - 2*
     &    hss(1,2)*hss(2,1)*hr**2*cb**2*s*tx**(-1) - 8*hss(2,2)**2*pq*
     &    ssp*sb**2 - 8*hss(2,2)**2*pq*ssp*cb**2 - 8*hss(2,2)**2*lq*ssz
     &    *sb**2*s*sz**(-1) - 8*hss(2,2)**2*rq*ssz*cb**2*s*sz**(-1) - 2
     &    *hss(2,2)**2*hl**2*sb**2*s*tx**(-1) - 2*hss(2,2)**2*hr**2*
     &    cb**2*s*tx**(-1) + 8*hhss(1,2)*hl*hr*sb*cb*mg*mt*s*tx**(-1)
     &     + 16*hhss(1,2)*h1*lambda1*sb*cb*sqrt2**(-1)*mg*s*s1**(-1) + 
     &    16*hhss(1,2)*h2*lambda2*sb*cb*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(1,1)*pq*hl
     &    *ssp*st*cb*m1**2*mg*tx**(-1) - 16*hss(1,1)*pq*hl*ssp*st*cb*mg
     &    *s**(-1)*t1*u1*tx**(-1) + 16*hss(1,1)*pq*hl*ssp*ct*sb*m1**2*
     &    mg*tx**(-1) - 16*hss(1,1)*pq*hl*ssp*ct*sb*mg*s**(-1)*t1*u1*
     &    tx**(-1) - 32*hss(1,1)*pq*hl*ssp*ct*cb*m1**2*mg**2*mt*
     &    t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hl*ssp*ct*cb*m1**2*mt*
     &    msb1**2*t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hl*ssp*ct*cb*mg**2
     &    *mt*s**(-1)*u1*tx**(-1) - 32*hss(1,1)*pq*hl*ssp*ct*cb*mt*
     &    msb1**2*s**(-1)*u1*tx**(-1) - 32*hss(1,1)*pq*hr*ssp*st*sb*
     &    m1**2*mg**2*mt*t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hr*ssp*st*
     &    sb*m1**2*mt*msb1**2*t1**(-1)*tx**(-1) + 32*hss(1,1)*pq*hr*ssp
     &    *st*sb*mg**2*mt*s**(-1)*u1*tx**(-1) - 32*hss(1,1)*pq*hr*ssp*
     &    st*sb*mt*msb1**2*s**(-1)*u1*tx**(-1) + 16*hss(1,1)*pq*hr*ssp*
     &    st*cb*m1**2*mg*tx**(-1) - 16*hss(1,1)*pq*hr*ssp*st*cb*mg*
     &    s**(-1)*t1*u1*tx**(-1) + 16*hss(1,1)*pq*hr*ssp*ct*sb*m1**2*mg
     &    *tx**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(1,1)*pq
     &    *hr*ssp*ct*sb*mg*s**(-1)*t1*u1*tx**(-1) + 16*hss(1,1)*lq*hl*
     &    ssz*st*cb*m1**2*mg*s*tx**(-1)*sz**(-1) - 16*hss(1,1)*lq*hl*
     &    ssz*st*cb*mg*t1*u1*tx**(-1)*sz**(-1) + 16*hss(1,1)*lq*hl*ssz*
     &    ct*sb*m1**2*mg*s*tx**(-1)*sz**(-1) - 16*hss(1,1)*lq*hl*ssz*ct
     &    *sb*mg*t1*u1*tx**(-1)*sz**(-1) - 32*hss(1,1)*lq*hl*ssz*ct*cb*
     &    m1**2*mg**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,1)*lq*
     &    hl*ssz*ct*cb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-1)*sz**(-1) + 
     &    32*hss(1,1)*lq*hl*ssz*ct*cb*mg**2*mt*u1*tx**(-1)*sz**(-1) - 
     &    32*hss(1,1)*lq*hl*ssz*ct*cb*mt*msb1**2*u1*tx**(-1)*sz**(-1)
     &     - 32*hss(1,1)*rq*hr*ssz*st*sb*m1**2*mg**2*mt*s*t1**(-1)*
     &    tx**(-1)*sz**(-1) + 32*hss(1,1)*rq*hr*ssz*st*sb*m1**2*mt*
     &    msb1**2*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,1)*rq*hr*ssz*
     &    st*sb*mg**2*mt*u1*tx**(-1)*sz**(-1) - 32*hss(1,1)*rq*hr*ssz*
     &    st*sb*mt*msb1**2*u1*tx**(-1)*sz**(-1) + 16*hss(1,1)*rq*hr*ssz
     &    *st*cb*m1**2*mg*s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(1,1)*rq
     &    *hr*ssz*st*cb*mg*t1*u1*tx**(-1)*sz**(-1) + 16*hss(1,1)*rq*hr*
     &    ssz*ct*sb*m1**2*mg*s*tx**(-1)*sz**(-1) - 16*hss(1,1)*rq*hr*
     &    ssz*ct*sb*mg*t1*u1*tx**(-1)*sz**(-1) - 4*hss(1,1)*hl*hr**2*st
     &    *cb*mg*mt**2*s*tx**(-2) - 4*hss(1,1)*hl*hr**2*ct*sb*mg*mt**2*
     &    s*tx**(-2) + 8*hss(1,1)*hl*hr**2*ct*cb*m1**2*mg**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,1)*hl*hr**2*ct*cb*m1**2*mt*
     &    msb1**2*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hl*hr**2*ct*cb*mg**2
     &    *mt*s*tx**(-2) - 8*hss(1,1)*hl*hr**2*ct*cb*mt*msb1**2*s*
     &    tx**(-2) + 16*hss(1,1)*hl*h1*lambda1*st*sb*sqrt2**(-1)*m1**2*
     &    mg**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(1,1)*hl*h1*
     &    lambda1*st*sb*sqrt2**(-1)*m1**2*msb1**2*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) + 16*hss(1,1)*hl*h1*lambda1*st*sb*sqrt2**(-1)*mg**2*
     &    s*tx**(-1)*s1**(-1) - 16*hss(1,1)*hl*h1*lambda1*st*sb*
     &    sqrt2**(-1)*msb1**2*s*tx**(-1)*s1**(-1) - 8*hss(1,1)*hl*h1*
     &    lambda1*st*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*hl*
     &    h1*lambda1*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) + 16*
     &    hss(1,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*m1**2*mg**2*s*
     &    t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(1,1)*hl*h2*lambda2*st*sb*
     &    sqrt2**(-1)*m1**2*msb1**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(1,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*mg**2*s*tx**(-1)*
     &    s2**(-1) - 16*hss(1,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*
     &    msb1**2*s*tx**(-1)*s2**(-1) - 8*hss(1,1)*hl*h2*lambda2*st*cb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) - 8*hss(1,1)*hl*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) + 8*hss(1
     &    ,1)*hl**2*hr*st*sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) - 8*
     &    hss(1,1)*hl**2*hr*st*sb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2)
     &     + 8*hss(1,1)*hl**2*hr*st*sb*mg**2*mt*s*tx**(-2) - 8*hss(1,1)
     &    *hl**2*hr*st*sb*mt*msb1**2*s*tx**(-2) - 4*hss(1,1)*hl**2*hr*
     &    st*cb*mg*mt**2*s*tx**(-2) - 4*hss(1,1)*hl**2*hr*ct*sb*mg*
     &    mt**2*s*tx**(-2) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * ( 4*hss(1,1)*hl**3*
     &    st*cb*m1**2*mg*s*tx**(-2) - 4*hss(1,1)*hl**3*st*cb*mg*t1*u1*
     &    tx**(-2) + 4*hss(1,1)*hl**3*ct*sb*m1**2*mg*s*tx**(-2) - 4*
     &    hss(1,1)*hl**3*ct*sb*mg*t1*u1*tx**(-2) - 8*hss(1,1)*hl**3*ct*
     &    cb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hl**3*ct*
     &    cb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hl**3*ct
     &    *cb*mg**2*mt*u1*tx**(-2) - 8*hss(1,1)*hl**3*ct*cb*mt*msb1**2*
     &    u1*tx**(-2) - 8*hss(1,1)*hr*h1*lambda1*st*cb*sqrt2**(-1)*mg*
     &    mt*s*tx**(-1)*s1**(-1) - 8*hss(1,1)*hr*h1*lambda1*ct*sb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) + 16*hss(1,1)*hr*h1*
     &    lambda1*ct*cb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) - 16*hss(1,1)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*m1**2*
     &    msb1**2*s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(1,1)*hr*h1*
     &    lambda1*ct*cb*sqrt2**(-1)*mg**2*s*tx**(-1)*s1**(-1) - 16*hss(
     &    1,1)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*msb1**2*s*tx**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*hr*
     &    h2*lambda2*st*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) - 8*
     &    hss(1,1)*hr*h2*lambda2*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s2**(-1) + 16*hss(1,1)*hr*h2*lambda2*ct*cb*sqrt2**(-1)*m1**2*
     &    mg**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(1,1)*hr*h2*
     &    lambda2*ct*cb*sqrt2**(-1)*m1**2*msb1**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) + 16*hss(1,1)*hr*h2*lambda2*ct*cb*sqrt2**(-1)*mg**2*
     &    s*tx**(-1)*s2**(-1) - 16*hss(1,1)*hr*h2*lambda2*ct*cb*
     &    sqrt2**(-1)*msb1**2*s*tx**(-1)*s2**(-1) - 8*hss(1,1)*hr**3*st
     &    *sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hr**3*st*
     &    sb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2) + 8*hss(1,1)*hr**3*st
     &    *sb*mg**2*mt*u1*tx**(-2) - 8*hss(1,1)*hr**3*st*sb*mt*msb1**2*
     &    u1*tx**(-2) + 4*hss(1,1)*hr**3*st*cb*m1**2*mg*s*tx**(-2) - 4*
     &    hss(1,1)*hr**3*st*cb*mg*t1*u1*tx**(-2) + 4*hss(1,1)*hr**3*ct*
     &    sb*m1**2*mg*s*tx**(-2) - 4*hss(1,1)*hr**3*ct*sb*mg*t1*u1*
     &    tx**(-2) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)**2*
     &    pq*ssp*sb**2*s**(-1)*t1 + 8*hss(1,1)**2*pq*ssp*sb**2*s**(-1)*
     &    u1 + 8*hss(1,1)**2*pq*ssp*sb**2 - 8*hss(1,1)**2*pq*ssp*cb**2*
     &    s**(-1)*t1 + 8*hss(1,1)**2*pq*ssp*cb**2*s**(-1)*u1 + 8*hss(1,
     &    1)**2*pq*ssp*cb**2 + 8*hss(1,1)**2*lq*ssz*cb**2*s*sz**(-1) - 
     &    8*hss(1,1)**2*lq*ssz*cb**2*t1*sz**(-1) + 8*hss(1,1)**2*lq*ssz
     &    *cb**2*u1*sz**(-1) + 8*hss(1,1)**2*rq*ssz*sb**2*s*sz**(-1) - 
     &    8*hss(1,1)**2*rq*ssz*sb**2*t1*sz**(-1) + 8*hss(1,1)**2*rq*ssz
     &    *sb**2*u1*sz**(-1) + 2*hss(1,1)**2*hl**2*cb**2*s*tx**(-1) - 2
     &    *hss(1,1)**2*hl**2*cb**2*t1*tx**(-1) + 2*hss(1,1)**2*hl**2*
     &    cb**2*u1*tx**(-1) + 2*hss(1,1)**2*hr**2*sb**2*s*tx**(-1) - 2*
     &    hss(1,1)**2*hr**2*sb**2*t1*tx**(-1) + 2*hss(1,1)**2*hr**2*
     &    sb**2*u1*tx**(-1) - 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*s*
     &    sz**(-1) - 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*u1*sz**(-1) + 8*
     &    hss(1,1)*hss(1,2)*rq*ssz*sb*cb*s*sz**(-1) + 8*hss(1,1)*hss(1,
     &    2)*rq*ssz*sb*cb*u1*sz**(-1) )
      MMs = MMs + SCC(2,2)*Nc*Cf*Pi*alphas*prefac * (  - 2*hss(1,1)*
     &    hss(1,2)*hl**2*sb*cb*s*tx**(-1) - 2*hss(1,1)*hss(1,2)*hl**2*
     &    sb*cb*u1*tx**(-1) + 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*s*
     &    tx**(-1) + 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*u1*tx**(-1) + 8*
     &    hss(1,1)*hss(2,1)*lq*ssz*sb*cb*t1*sz**(-1) - 8*hss(1,1)*hss(2
     &    ,1)*rq*ssz*sb*cb*t1*sz**(-1) + 2*hss(1,1)*hss(2,1)*hl**2*sb*
     &    cb*t1*tx**(-1) - 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*t1*tx**(-1)
     &     )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,2)*
     &    hss(2,1)*pq*ssp*sb**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*pq*ssp
     &    *sb**2*s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2 - 8*hss(
     &    1,2)*hss(2,1)*pq*ssp*cb**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*
     &    pq*ssp*cb**2*s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2 + 
     &    8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*s*sz**(-1) - 8*hss(1,2)*hss(
     &    2,1)*lq*ssz*cb**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*
     &    cb**2*u1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*s*
     &    sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*t1*sz**(-1) + 8*
     &    hss(1,2)*hss(2,1)*rq*ssz*sb**2*u1*sz**(-1) + 2*hss(1,2)*hss(2
     &    ,1)*hl**2*cb**2*s*tx**(-1) - 2*hss(1,2)*hss(2,1)*hl**2*cb**2*
     &    t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*cb**2*u1*tx**(-1) + 2
     &    *hss(1,2)*hss(2,1)*hr**2*sb**2*s*tx**(-1) - 2*hss(1,2)*hss(2,
     &    1)*hr**2*sb**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*sb**2*
     &    u1*tx**(-1) - 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*s*sz**(-1) - 8
     &    *hss(1,2)*hss(2,2)*lq*ssz*sb*cb*u1*sz**(-1) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    2)*rq*ssz*sb*cb*s*sz**(-1) + 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb
     &    *u1*sz**(-1) - 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*s*tx**(-1) - 2
     &    *hss(1,2)*hss(2,2)*hl**2*sb*cb*u1*tx**(-1) + 2*hss(1,2)*hss(2
     &    ,2)*hr**2*sb*cb*s*tx**(-1) + 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*
     &    u1*tx**(-1) - 16*hss(2,1)*pq*hl*ssp*st*sb*m1**2*mg*tx**(-1)
     &     + 16*hss(2,1)*pq*hl*ssp*st*sb*mg*s**(-1)*t1*u1*tx**(-1) + 32
     &    *hss(2,1)*pq*hl*ssp*st*cb*m1**2*mg**2*mt*t1**(-1)*tx**(-1) - 
     &    32*hss(2,1)*pq*hl*ssp*st*cb*m1**2*mt*msb1**2*t1**(-1)*
     &    tx**(-1) - 32*hss(2,1)*pq*hl*ssp*st*cb*mg**2*mt*s**(-1)*u1*
     &    tx**(-1) + 32*hss(2,1)*pq*hl*ssp*st*cb*mt*msb1**2*s**(-1)*u1*
     &    tx**(-1) + 16*hss(2,1)*pq*hl*ssp*ct*cb*m1**2*mg*tx**(-1) - 16
     &    *hss(2,1)*pq*hl*ssp*ct*cb*mg*s**(-1)*t1*u1*tx**(-1) - 16*hss(
     &    2,1)*pq*hr*ssp*st*sb*m1**2*mg*tx**(-1) + 16*hss(2,1)*pq*hr*
     &    ssp*st*sb*mg*s**(-1)*t1*u1*tx**(-1) - 32*hss(2,1)*pq*hr*ssp*
     &    ct*sb*m1**2*mg**2*mt*t1**(-1)*tx**(-1) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * ( 32*hss(2,1)*pq*hr
     &    *ssp*ct*sb*m1**2*mt*msb1**2*t1**(-1)*tx**(-1) + 32*hss(2,1)*
     &    pq*hr*ssp*ct*sb*mg**2*mt*s**(-1)*u1*tx**(-1) - 32*hss(2,1)*pq
     &    *hr*ssp*ct*sb*mt*msb1**2*s**(-1)*u1*tx**(-1) + 16*hss(2,1)*pq
     &    *hr*ssp*ct*cb*m1**2*mg*tx**(-1) - 16*hss(2,1)*pq*hr*ssp*ct*cb
     &    *mg*s**(-1)*t1*u1*tx**(-1) - 16*hss(2,1)*lq*hl*ssz*st*sb*
     &    m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(2,1)*lq*hl*ssz*st*sb*mg
     &    *t1*u1*tx**(-1)*sz**(-1) + 32*hss(2,1)*lq*hl*ssz*st*cb*m1**2*
     &    mg**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(2,1)*lq*hl*ssz
     &    *st*cb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*
     &    hss(2,1)*lq*hl*ssz*st*cb*mg**2*mt*u1*tx**(-1)*sz**(-1) + 32*
     &    hss(2,1)*lq*hl*ssz*st*cb*mt*msb1**2*u1*tx**(-1)*sz**(-1) + 16
     &    *hss(2,1)*lq*hl*ssz*ct*cb*m1**2*mg*s*tx**(-1)*sz**(-1) - 16*
     &    hss(2,1)*lq*hl*ssz*ct*cb*mg*t1*u1*tx**(-1)*sz**(-1) - 16*hss(
     &    2,1)*rq*hr*ssz*st*sb*m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(2,
     &    1)*rq*hr*ssz*st*sb*mg*t1*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * (  - 32*hss(2,1)*rq
     &    *hr*ssz*ct*sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 
     &    32*hss(2,1)*rq*hr*ssz*ct*sb*m1**2*mt*msb1**2*s*t1**(-1)*
     &    tx**(-1)*sz**(-1) + 32*hss(2,1)*rq*hr*ssz*ct*sb*mg**2*mt*u1*
     &    tx**(-1)*sz**(-1) - 32*hss(2,1)*rq*hr*ssz*ct*sb*mt*msb1**2*u1
     &    *tx**(-1)*sz**(-1) + 16*hss(2,1)*rq*hr*ssz*ct*cb*m1**2*mg*s*
     &    tx**(-1)*sz**(-1) - 16*hss(2,1)*rq*hr*ssz*ct*cb*mg*t1*u1*
     &    tx**(-1)*sz**(-1) + 4*hss(2,1)*hl*hr**2*st*sb*mg*mt**2*s*
     &    tx**(-2) - 8*hss(2,1)*hl*hr**2*st*cb*m1**2*mg**2*mt*s*
     &    t1**(-1)*tx**(-2) + 8*hss(2,1)*hl*hr**2*st*cb*m1**2*mt*
     &    msb1**2*s*t1**(-1)*tx**(-2) - 8*hss(2,1)*hl*hr**2*st*cb*mg**2
     &    *mt*s*tx**(-2) + 8*hss(2,1)*hl*hr**2*st*cb*mt*msb1**2*s*
     &    tx**(-2) - 4*hss(2,1)*hl*hr**2*ct*cb*mg*mt**2*s*tx**(-2) + 8*
     &    hss(2,1)*hl*h1*lambda1*st*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s1**(-1) + 16*hss(2,1)*hl*h1*lambda1*ct*sb*sqrt2**(-1)*m1**2*
     &    mg**2*s*t1**(-1)*tx**(-1)*s1**(-1) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(2,1)*hl
     &    *h1*lambda1*ct*sb*sqrt2**(-1)*m1**2*msb1**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) + 16*hss(2,1)*hl*h1*lambda1*ct*sb*
     &    sqrt2**(-1)*mg**2*s*tx**(-1)*s1**(-1) - 16*hss(2,1)*hl*h1*
     &    lambda1*ct*sb*sqrt2**(-1)*msb1**2*s*tx**(-1)*s1**(-1) - 8*
     &    hss(2,1)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s1**(-1) + 8*hss(2,1)*hl*h2*lambda2*st*sb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s2**(-1) + 16*hss(2,1)*hl*h2*lambda2*ct*sb*
     &    sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*
     &    hss(2,1)*hl*h2*lambda2*ct*sb*sqrt2**(-1)*m1**2*msb1**2*s*
     &    t1**(-1)*tx**(-1)*s2**(-1) + 16*hss(2,1)*hl*h2*lambda2*ct*sb*
     &    sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) - 16*hss(2,1)*hl*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*msb1**2*s*tx**(-1)*s2**(-1) - 8*
     &    hss(2,1)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s2**(-1) + 4*hss(2,1)*hl**2*hr*st*sb*mg*mt**2*s*tx**(-2) + 8*
     &    hss(2,1)*hl**2*hr*ct*sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(2,1)*
     &    hl**2*hr*ct*sb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2) + 8*hss(2
     &    ,1)*hl**2*hr*ct*sb*mg**2*mt*s*tx**(-2) - 8*hss(2,1)*hl**2*hr*
     &    ct*sb*mt*msb1**2*s*tx**(-2) - 4*hss(2,1)*hl**2*hr*ct*cb*mg*
     &    mt**2*s*tx**(-2) - 4*hss(2,1)*hl**3*st*sb*m1**2*mg*s*tx**(-2)
     &     + 4*hss(2,1)*hl**3*st*sb*mg*t1*u1*tx**(-2) + 8*hss(2,1)*
     &    hl**3*st*cb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(2,1)*
     &    hl**3*st*cb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2) - 8*hss(2,1)
     &    *hl**3*st*cb*mg**2*mt*u1*tx**(-2) + 8*hss(2,1)*hl**3*st*cb*mt
     &    *msb1**2*u1*tx**(-2) + 4*hss(2,1)*hl**3*ct*cb*m1**2*mg*s*
     &    tx**(-2) - 4*hss(2,1)*hl**3*ct*cb*mg*t1*u1*tx**(-2) + 8*hss(2
     &    ,1)*hr*h1*lambda1*st*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1)
     &     - 16*hss(2,1)*hr*h1*lambda1*st*cb*sqrt2**(-1)*m1**2*mg**2*s*
     &    t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(2,1)*hr*h1*lambda1*st*cb*
     &    sqrt2**(-1)*m1**2*msb1**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*
     &    hss(2,1)*hr*h1*lambda1*st*cb*sqrt2**(-1)*mg**2*s*tx**(-1)*
     &    s1**(-1) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(2,1)*hr*h1
     &    *lambda1*st*cb*sqrt2**(-1)*msb1**2*s*tx**(-1)*s1**(-1) - 8*
     &    hss(2,1)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s1**(-1) + 8*hss(2,1)*hr*h2*lambda2*st*sb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s2**(-1) - 16*hss(2,1)*hr*h2*lambda2*st*cb*
     &    sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*
     &    hss(2,1)*hr*h2*lambda2*st*cb*sqrt2**(-1)*m1**2*msb1**2*s*
     &    t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(2,1)*hr*h2*lambda2*st*cb*
     &    sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) + 16*hss(2,1)*hr*h2*
     &    lambda2*st*cb*sqrt2**(-1)*msb1**2*s*tx**(-1)*s2**(-1) - 8*
     &    hss(2,1)*hr*h2*lambda2*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s2**(-1) - 4*hss(2,1)*hr**3*st*sb*m1**2*mg*s*tx**(-2) + 4*
     &    hss(2,1)*hr**3*st*sb*mg*t1*u1*tx**(-2) - 8*hss(2,1)*hr**3*ct*
     &    sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(2,1)*hr**3*ct*
     &    sb*m1**2*mt*msb1**2*s*t1**(-1)*tx**(-2) + 8*hss(2,1)*hr**3*ct
     &    *sb*mg**2*mt*u1*tx**(-2) )
      MMs = MMs + SCC(2,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(2,1)*
     &    hr**3*ct*sb*mt*msb1**2*u1*tx**(-2) + 4*hss(2,1)*hr**3*ct*cb*
     &    m1**2*mg*s*tx**(-2) - 4*hss(2,1)*hr**3*ct*cb*mg*t1*u1*
     &    tx**(-2) + 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*t1*sz**(-1) - 8*
     &    hss(2,1)*hss(2,2)*rq*ssz*sb*cb*t1*sz**(-1) + 2*hss(2,1)*hss(2
     &    ,2)*hl**2*sb*cb*t1*tx**(-1) - 2*hss(2,1)*hss(2,2)*hr**2*sb*cb
     &    *t1*tx**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,1)*hss(1,
     &    2)*lq*ssz*sb*cb*t1*sz**(-1) - 8*hss(1,1)*hss(1,2)*rq*ssz*sb*
     &    cb*t1*sz**(-1) + 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*t1*tx**(-1)
     &     - 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*t1*tx**(-1) - 8*hss(1,1)*
     &    hss(2,1)*lq*ssz*sb*cb*s*sz**(-1) - 8*hss(1,1)*hss(2,1)*lq*ssz
     &    *sb*cb*u1*sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*s*
     &    sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(1,1)*hss(2,1)*hl**2*sb*cb*s*tx**(-1) - 2*hss(1,1)*hss(2,1
     &    )*hl**2*sb*cb*u1*tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*s
     &    *tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*u1*tx**(-1) - 16*
     &    hss(1,2)*pq*hl*ssp*st*sb*m1**2*mg*tx**(-1) + 16*hss(1,2)*pq*
     &    hl*ssp*st*sb*mg*s**(-1)*t1*u1*tx**(-1) + 32*hss(1,2)*pq*hl*
     &    ssp*ct*sb*m1**2*mg**2*mt*t1**(-1)*tx**(-1) - 32*hss(1,2)*pq*
     &    hl*ssp*ct*sb*m1**2*mt*msb2**2*t1**(-1)*tx**(-1) - 32*hss(1,2)
     &    *pq*hl*ssp*ct*sb*mg**2*mt*s**(-1)*u1*tx**(-1) + 32*hss(1,2)*
     &    pq*hl*ssp*ct*sb*mt*msb2**2*s**(-1)*u1*tx**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(1,2)*pq*hl
     &    *ssp*ct*cb*m1**2*mg*tx**(-1) - 16*hss(1,2)*pq*hl*ssp*ct*cb*mg
     &    *s**(-1)*t1*u1*tx**(-1) - 16*hss(1,2)*pq*hr*ssp*st*sb*m1**2*
     &    mg*tx**(-1) + 16*hss(1,2)*pq*hr*ssp*st*sb*mg*s**(-1)*t1*u1*
     &    tx**(-1) - 32*hss(1,2)*pq*hr*ssp*st*cb*m1**2*mg**2*mt*
     &    t1**(-1)*tx**(-1) + 32*hss(1,2)*pq*hr*ssp*st*cb*m1**2*mt*
     &    msb2**2*t1**(-1)*tx**(-1) + 32*hss(1,2)*pq*hr*ssp*st*cb*mg**2
     &    *mt*s**(-1)*u1*tx**(-1) - 32*hss(1,2)*pq*hr*ssp*st*cb*mt*
     &    msb2**2*s**(-1)*u1*tx**(-1) + 16*hss(1,2)*pq*hr*ssp*ct*cb*
     &    m1**2*mg*tx**(-1) - 16*hss(1,2)*pq*hr*ssp*ct*cb*mg*s**(-1)*t1
     &    *u1*tx**(-1) - 16*hss(1,2)*lq*hl*ssz*st*sb*m1**2*mg*s*
     &    tx**(-1)*sz**(-1) + 16*hss(1,2)*lq*hl*ssz*st*sb*mg*t1*u1*
     &    tx**(-1)*sz**(-1) + 32*hss(1,2)*lq*hl*ssz*ct*sb*m1**2*mg**2*
     &    mt*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,2)*lq*hl*ssz*ct*sb
     &    *m1**2*mt*msb2**2*s*t1**(-1)*tx**(-1)*sz**(-1) - 32*hss(1,2)*
     &    lq*hl*ssz*ct*sb*mg**2*mt*u1*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * ( 32*hss(1,2)*lq*hl
     &    *ssz*ct*sb*mt*msb2**2*u1*tx**(-1)*sz**(-1) + 16*hss(1,2)*lq*
     &    hl*ssz*ct*cb*m1**2*mg*s*tx**(-1)*sz**(-1) - 16*hss(1,2)*lq*hl
     &    *ssz*ct*cb*mg*t1*u1*tx**(-1)*sz**(-1) - 16*hss(1,2)*rq*hr*ssz
     &    *st*sb*m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(1,2)*rq*hr*ssz*
     &    st*sb*mg*t1*u1*tx**(-1)*sz**(-1) - 32*hss(1,2)*rq*hr*ssz*st*
     &    cb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(1,2)*
     &    rq*hr*ssz*st*cb*m1**2*mt*msb2**2*s*t1**(-1)*tx**(-1)*sz**(-1)
     &     + 32*hss(1,2)*rq*hr*ssz*st*cb*mg**2*mt*u1*tx**(-1)*sz**(-1)
     &     - 32*hss(1,2)*rq*hr*ssz*st*cb*mt*msb2**2*u1*tx**(-1)*
     &    sz**(-1) + 16*hss(1,2)*rq*hr*ssz*ct*cb*m1**2*mg*s*tx**(-1)*
     &    sz**(-1) - 16*hss(1,2)*rq*hr*ssz*ct*cb*mg*t1*u1*tx**(-1)*
     &    sz**(-1) + 4*hss(1,2)*hl*hr**2*st*sb*mg*mt**2*s*tx**(-2) - 8*
     &    hss(1,2)*hl*hr**2*ct*sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 
     &    8*hss(1,2)*hl*hr**2*ct*sb*m1**2*mt*msb2**2*s*t1**(-1)*
     &    tx**(-2) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,2)*hl*
     &    hr**2*ct*sb*mg**2*mt*s*tx**(-2) + 8*hss(1,2)*hl*hr**2*ct*sb*
     &    mt*msb2**2*s*tx**(-2) - 4*hss(1,2)*hl*hr**2*ct*cb*mg*mt**2*s*
     &    tx**(-2) + 8*hss(1,2)*hl*h1*lambda1*st*sb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s1**(-1) + 16*hss(1,2)*hl*h1*lambda1*st*cb*
     &    sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*
     &    hss(1,2)*hl*h1*lambda1*st*cb*sqrt2**(-1)*m1**2*msb2**2*s*
     &    t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(1,2)*hl*h1*lambda1*st*cb*
     &    sqrt2**(-1)*mg**2*s*tx**(-1)*s1**(-1) - 16*hss(1,2)*hl*h1*
     &    lambda1*st*cb*sqrt2**(-1)*msb2**2*s*tx**(-1)*s1**(-1) - 8*
     &    hss(1,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*
     &    s1**(-1) + 8*hss(1,2)*hl*h2*lambda2*st*sb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s2**(-1) + 16*hss(1,2)*hl*h2*lambda2*st*cb*
     &    sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*
     &    hss(1,2)*hl*h2*lambda2*st*cb*sqrt2**(-1)*m1**2*msb2**2*s*
     &    t1**(-1)*tx**(-1)*s2**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(1,2)*hl*h2
     &    *lambda2*st*cb*sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) - 16*
     &    hss(1,2)*hl*h2*lambda2*st*cb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s2**(-1) - 8*hss(1,2)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s2**(-1) + 4*hss(1,2)*hl**2*hr*st*sb*mg*mt**2*s*
     &    tx**(-2) + 8*hss(1,2)*hl**2*hr*st*cb*m1**2*mg**2*mt*s*
     &    t1**(-1)*tx**(-2) - 8*hss(1,2)*hl**2*hr*st*cb*m1**2*mt*
     &    msb2**2*s*t1**(-1)*tx**(-2) + 8*hss(1,2)*hl**2*hr*st*cb*mg**2
     &    *mt*s*tx**(-2) - 8*hss(1,2)*hl**2*hr*st*cb*mt*msb2**2*s*
     &    tx**(-2) - 4*hss(1,2)*hl**2*hr*ct*cb*mg*mt**2*s*tx**(-2) - 4*
     &    hss(1,2)*hl**3*st*sb*m1**2*mg*s*tx**(-2) + 4*hss(1,2)*hl**3*
     &    st*sb*mg*t1*u1*tx**(-2) + 8*hss(1,2)*hl**3*ct*sb*m1**2*mg**2*
     &    mt*s*t1**(-1)*tx**(-2) - 8*hss(1,2)*hl**3*ct*sb*m1**2*mt*
     &    msb2**2*s*t1**(-1)*tx**(-2) - 8*hss(1,2)*hl**3*ct*sb*mg**2*mt
     &    *u1*tx**(-2) + 8*hss(1,2)*hl**3*ct*sb*mt*msb2**2*u1*tx**(-2)
     &     + 4*hss(1,2)*hl**3*ct*cb*m1**2*mg*s*tx**(-2) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * (  - 4*hss(1,2)*
     &    hl**3*ct*cb*mg*t1*u1*tx**(-2) + 8*hss(1,2)*hr*h1*lambda1*st*
     &    sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) - 16*hss(1,2)*hr*h1*
     &    lambda1*ct*sb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) + 16*hss(1,2)*hr*h1*lambda1*ct*sb*sqrt2**(-1)*m1**2*
     &    msb2**2*s*t1**(-1)*tx**(-1)*s1**(-1) - 16*hss(1,2)*hr*h1*
     &    lambda1*ct*sb*sqrt2**(-1)*mg**2*s*tx**(-1)*s1**(-1) + 16*hss(
     &    1,2)*hr*h1*lambda1*ct*sb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s1**(-1) - 8*hss(1,2)*hr*h1*lambda1*ct*cb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s1**(-1) + 8*hss(1,2)*hr*h2*lambda2*st*sb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) - 16*hss(1,2)*hr*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) + 16*hss(1,2)*hr*h2*lambda2*ct*sb*sqrt2**(-1)*m1**2*
     &    msb2**2*s*t1**(-1)*tx**(-1)*s2**(-1) - 16*hss(1,2)*hr*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) + 16*hss(
     &    1,2)*hr*h2*lambda2*ct*sb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s2**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,2)*hr*
     &    h2*lambda2*ct*cb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) - 4*
     &    hss(1,2)*hr**3*st*sb*m1**2*mg*s*tx**(-2) + 4*hss(1,2)*hr**3*
     &    st*sb*mg*t1*u1*tx**(-2) - 8*hss(1,2)*hr**3*st*cb*m1**2*mg**2*
     &    mt*s*t1**(-1)*tx**(-2) + 8*hss(1,2)*hr**3*st*cb*m1**2*mt*
     &    msb2**2*s*t1**(-1)*tx**(-2) + 8*hss(1,2)*hr**3*st*cb*mg**2*mt
     &    *u1*tx**(-2) - 8*hss(1,2)*hr**3*st*cb*mt*msb2**2*u1*tx**(-2)
     &     + 4*hss(1,2)*hr**3*ct*cb*m1**2*mg*s*tx**(-2) - 4*hss(1,2)*
     &    hr**3*ct*cb*mg*t1*u1*tx**(-2) - 8*hss(1,2)*hss(2,1)*pq*ssp*
     &    sb**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*s**(-1)*
     &    u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2 - 8*hss(1,2)*hss(2,1)*
     &    pq*ssp*cb**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*
     &    s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2 + 8*hss(1,2)*
     &    hss(2,1)*lq*ssz*sb**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*ssz
     &    *sb**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*u1*
     &    sz**(-1) )
      MMs = MMs + SCC(2,4)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*rq*ssz*cb**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2
     &    *t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*u1*sz**(-1)
     &     + 2*hss(1,2)*hss(2,1)*hl**2*sb**2*s*tx**(-1) - 2*hss(1,2)*
     &    hss(2,1)*hl**2*sb**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*
     &    sb**2*u1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*cb**2*s*
     &    tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*cb**2*t1*tx**(-1) + 2*
     &    hss(1,2)*hss(2,1)*hr**2*cb**2*u1*tx**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    2)*lq*ssz*sb*cb*t1*sz**(-1) - 8*hss(1,2)*hss(2,2)*rq*ssz*sb*
     &    cb*t1*sz**(-1) + 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*t1*tx**(-1)
     &     - 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*t1*tx**(-1) - 8*hss(2,1)*
     &    hss(2,2)*lq*ssz*sb*cb*s*sz**(-1) - 8*hss(2,1)*hss(2,2)*lq*ssz
     &    *sb*cb*u1*sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*s*
     &    sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(2,1)*hss(2,2)*hl**2*sb*cb*s*tx**(-1) - 2*hss(2,1)*hss(2,2
     &    )*hl**2*sb*cb*u1*tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*s
     &    *tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*u1*tx**(-1) - 32*
     &    hss(2,2)*pq*hl*ssp*st*sb*m1**2*mg**2*mt*t1**(-1)*tx**(-1) + 
     &    32*hss(2,2)*pq*hl*ssp*st*sb*m1**2*mt*msb2**2*t1**(-1)*
     &    tx**(-1) + 32*hss(2,2)*pq*hl*ssp*st*sb*mg**2*mt*s**(-1)*u1*
     &    tx**(-1) - 32*hss(2,2)*pq*hl*ssp*st*sb*mt*msb2**2*s**(-1)*u1*
     &    tx**(-1) - 16*hss(2,2)*pq*hl*ssp*st*cb*m1**2*mg*tx**(-1) + 16
     &    *hss(2,2)*pq*hl*ssp*st*cb*mg*s**(-1)*t1*u1*tx**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * (  - 16*hss(2,2)*pq
     &    *hl*ssp*ct*sb*m1**2*mg*tx**(-1) + 16*hss(2,2)*pq*hl*ssp*ct*sb
     &    *mg*s**(-1)*t1*u1*tx**(-1) - 16*hss(2,2)*pq*hr*ssp*st*cb*
     &    m1**2*mg*tx**(-1) + 16*hss(2,2)*pq*hr*ssp*st*cb*mg*s**(-1)*t1
     &    *u1*tx**(-1) - 16*hss(2,2)*pq*hr*ssp*ct*sb*m1**2*mg*tx**(-1)
     &     + 16*hss(2,2)*pq*hr*ssp*ct*sb*mg*s**(-1)*t1*u1*tx**(-1) - 32
     &    *hss(2,2)*pq*hr*ssp*ct*cb*m1**2*mg**2*mt*t1**(-1)*tx**(-1) + 
     &    32*hss(2,2)*pq*hr*ssp*ct*cb*m1**2*mt*msb2**2*t1**(-1)*
     &    tx**(-1) + 32*hss(2,2)*pq*hr*ssp*ct*cb*mg**2*mt*s**(-1)*u1*
     &    tx**(-1) - 32*hss(2,2)*pq*hr*ssp*ct*cb*mt*msb2**2*s**(-1)*u1*
     &    tx**(-1) - 32*hss(2,2)*lq*hl*ssz*st*sb*m1**2*mg**2*mt*s*
     &    t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,2)*lq*hl*ssz*st*sb*
     &    m1**2*mt*msb2**2*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,2)*
     &    lq*hl*ssz*st*sb*mg**2*mt*u1*tx**(-1)*sz**(-1) - 32*hss(2,2)*
     &    lq*hl*ssz*st*sb*mt*msb2**2*u1*tx**(-1)*sz**(-1) - 16*hss(2,2)
     &    *lq*hl*ssz*st*cb*m1**2*mg*s*tx**(-1)*sz**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(2,2)*lq*hl
     &    *ssz*st*cb*mg*t1*u1*tx**(-1)*sz**(-1) - 16*hss(2,2)*lq*hl*ssz
     &    *ct*sb*m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(2,2)*lq*hl*ssz*
     &    ct*sb*mg*t1*u1*tx**(-1)*sz**(-1) - 16*hss(2,2)*rq*hr*ssz*st*
     &    cb*m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(2,2)*rq*hr*ssz*st*cb
     &    *mg*t1*u1*tx**(-1)*sz**(-1) - 16*hss(2,2)*rq*hr*ssz*ct*sb*
     &    m1**2*mg*s*tx**(-1)*sz**(-1) + 16*hss(2,2)*rq*hr*ssz*ct*sb*mg
     &    *t1*u1*tx**(-1)*sz**(-1) - 32*hss(2,2)*rq*hr*ssz*ct*cb*m1**2*
     &    mg**2*mt*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*hss(2,2)*rq*hr*ssz
     &    *ct*cb*m1**2*mt*msb2**2*s*t1**(-1)*tx**(-1)*sz**(-1) + 32*
     &    hss(2,2)*rq*hr*ssz*ct*cb*mg**2*mt*u1*tx**(-1)*sz**(-1) - 32*
     &    hss(2,2)*rq*hr*ssz*ct*cb*mt*msb2**2*u1*tx**(-1)*sz**(-1) + 8*
     &    hss(2,2)*hl*hr**2*st*sb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) - 
     &    8*hss(2,2)*hl*hr**2*st*sb*m1**2*mt*msb2**2*s*t1**(-1)*
     &    tx**(-2) + 8*hss(2,2)*hl*hr**2*st*sb*mg**2*mt*s*tx**(-2) - 8*
     &    hss(2,2)*hl*hr**2*st*sb*mt*msb2**2*s*tx**(-2) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 4*hss(2,2)*hl*
     &    hr**2*st*cb*mg*mt**2*s*tx**(-2) + 4*hss(2,2)*hl*hr**2*ct*sb*
     &    mg*mt**2*s*tx**(-2) + 8*hss(2,2)*hl*h1*lambda1*st*cb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) + 8*hss(2,2)*hl*h1*
     &    lambda1*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) + 16*hss(
     &    2,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*
     &    tx**(-1)*s1**(-1) - 16*hss(2,2)*hl*h1*lambda1*ct*cb*
     &    sqrt2**(-1)*m1**2*msb2**2*s*t1**(-1)*tx**(-1)*s1**(-1) + 16*
     &    hss(2,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*mg**2*s*tx**(-1)*
     &    s1**(-1) - 16*hss(2,2)*hl*h1*lambda1*ct*cb*sqrt2**(-1)*
     &    msb2**2*s*tx**(-1)*s1**(-1) + 8*hss(2,2)*hl*h2*lambda2*st*cb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) + 8*hss(2,2)*hl*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) + 16*hss(
     &    2,2)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*
     &    tx**(-1)*s2**(-1) - 16*hss(2,2)*hl*h2*lambda2*ct*cb*
     &    sqrt2**(-1)*m1**2*msb2**2*s*t1**(-1)*tx**(-1)*s2**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(2,2)*hl*h2
     &    *lambda2*ct*cb*sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) - 16*
     &    hss(2,2)*hl*h2*lambda2*ct*cb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s2**(-1) + 4*hss(2,2)*hl**2*hr*st*cb*mg*mt**2*s*tx**(-2) + 4*
     &    hss(2,2)*hl**2*hr*ct*sb*mg*mt**2*s*tx**(-2) + 8*hss(2,2)*
     &    hl**2*hr*ct*cb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) - 8*hss(2,2
     &    )*hl**2*hr*ct*cb*m1**2*mt*msb2**2*s*t1**(-1)*tx**(-2) + 8*
     &    hss(2,2)*hl**2*hr*ct*cb*mg**2*mt*s*tx**(-2) - 8*hss(2,2)*
     &    hl**2*hr*ct*cb*mt*msb2**2*s*tx**(-2) - 8*hss(2,2)*hl**3*st*sb
     &    *m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(2,2)*hl**3*st*sb*
     &    m1**2*mt*msb2**2*s*t1**(-1)*tx**(-2) + 8*hss(2,2)*hl**3*st*sb
     &    *mg**2*mt*u1*tx**(-2) - 8*hss(2,2)*hl**3*st*sb*mt*msb2**2*u1*
     &    tx**(-2) - 4*hss(2,2)*hl**3*st*cb*m1**2*mg*s*tx**(-2) + 4*
     &    hss(2,2)*hl**3*st*cb*mg*t1*u1*tx**(-2) - 4*hss(2,2)*hl**3*ct*
     &    sb*m1**2*mg*s*tx**(-2) + 4*hss(2,2)*hl**3*ct*sb*mg*t1*u1*
     &    tx**(-2) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 16*hss(2,2)*hr*h1
     &    *lambda1*st*sb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*
     &    s1**(-1) - 16*hss(2,2)*hr*h1*lambda1*st*sb*sqrt2**(-1)*m1**2*
     &    msb2**2*s*t1**(-1)*tx**(-1)*s1**(-1) + 16*hss(2,2)*hr*h1*
     &    lambda1*st*sb*sqrt2**(-1)*mg**2*s*tx**(-1)*s1**(-1) - 16*hss(
     &    2,2)*hr*h1*lambda1*st*sb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s1**(-1) + 8*hss(2,2)*hr*h1*lambda1*st*cb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s1**(-1) + 8*hss(2,2)*hr*h1*lambda1*ct*sb*
     &    sqrt2**(-1)*mg*mt*s*tx**(-1)*s1**(-1) + 16*hss(2,2)*hr*h2*
     &    lambda2*st*sb*sqrt2**(-1)*m1**2*mg**2*s*t1**(-1)*tx**(-1)*
     &    s2**(-1) - 16*hss(2,2)*hr*h2*lambda2*st*sb*sqrt2**(-1)*m1**2*
     &    msb2**2*s*t1**(-1)*tx**(-1)*s2**(-1) + 16*hss(2,2)*hr*h2*
     &    lambda2*st*sb*sqrt2**(-1)*mg**2*s*tx**(-1)*s2**(-1) - 16*hss(
     &    2,2)*hr*h2*lambda2*st*sb*sqrt2**(-1)*msb2**2*s*tx**(-1)*
     &    s2**(-1) + 8*hss(2,2)*hr*h2*lambda2*st*cb*sqrt2**(-1)*mg*mt*s
     &    *tx**(-1)*s2**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)*hr*h2*
     &    lambda2*ct*sb*sqrt2**(-1)*mg*mt*s*tx**(-1)*s2**(-1) - 4*hss(2
     &    ,2)*hr**3*st*cb*m1**2*mg*s*tx**(-2) + 4*hss(2,2)*hr**3*st*cb*
     &    mg*t1*u1*tx**(-2) - 4*hss(2,2)*hr**3*ct*sb*m1**2*mg*s*
     &    tx**(-2) + 4*hss(2,2)*hr**3*ct*sb*mg*t1*u1*tx**(-2) - 8*hss(2
     &    ,2)*hr**3*ct*cb*m1**2*mg**2*mt*s*t1**(-1)*tx**(-2) + 8*hss(2,
     &    2)*hr**3*ct*cb*m1**2*mt*msb2**2*s*t1**(-1)*tx**(-2) + 8*hss(2
     &    ,2)*hr**3*ct*cb*mg**2*mt*u1*tx**(-2) - 8*hss(2,2)*hr**3*ct*cb
     &    *mt*msb2**2*u1*tx**(-2) - 8*hss(2,2)**2*pq*ssp*sb**2*s**(-1)*
     &    t1 + 8*hss(2,2)**2*pq*ssp*sb**2*s**(-1)*u1 + 8*hss(2,2)**2*pq
     &    *ssp*sb**2 - 8*hss(2,2)**2*pq*ssp*cb**2*s**(-1)*t1 + 8*hss(2,
     &    2)**2*pq*ssp*cb**2*s**(-1)*u1 + 8*hss(2,2)**2*pq*ssp*cb**2 + 
     &    8*hss(2,2)**2*lq*ssz*sb**2*s*sz**(-1) - 8*hss(2,2)**2*lq*ssz*
     &    sb**2*t1*sz**(-1) + 8*hss(2,2)**2*lq*ssz*sb**2*u1*sz**(-1) + 
     &    8*hss(2,2)**2*rq*ssz*cb**2*s*sz**(-1) - 8*hss(2,2)**2*rq*ssz*
     &    cb**2*t1*sz**(-1) )
      MMs = MMs + SCC(2,5)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)**2*rq*
     &    ssz*cb**2*u1*sz**(-1) + 2*hss(2,2)**2*hl**2*sb**2*s*tx**(-1)
     &     - 2*hss(2,2)**2*hl**2*sb**2*t1*tx**(-1) + 2*hss(2,2)**2*
     &    hl**2*sb**2*u1*tx**(-1) + 2*hss(2,2)**2*hr**2*cb**2*s*
     &    tx**(-1) - 2*hss(2,2)**2*hr**2*cb**2*t1*tx**(-1) + 2*hss(2,2)
     &    **2*hr**2*cb**2*u1*tx**(-1) )
      MMs = MMs + SCC(6,2)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,1)**2*pq*
     &    ssp*sb**2*s**(-1)*t1 - 8*hss(1,1)**2*pq*ssp*sb**2*s**(-1)*u1
     &     + 8*hss(1,1)**2*pq*ssp*cb**2*s**(-1)*t1 - 8*hss(1,1)**2*pq*
     &    ssp*cb**2*s**(-1)*u1 + 8*hss(1,1)**2*lq*ssz*cb**2*t1*sz**(-1)
     &     - 8*hss(1,1)**2*lq*ssz*cb**2*u1*sz**(-1) + 8*hss(1,1)**2*rq*
     &    ssz*sb**2*t1*sz**(-1) - 8*hss(1,1)**2*rq*ssz*sb**2*u1*
     &    sz**(-1) + 2*hss(1,1)**2*hl**2*cb**2*t1*tx**(-1) - 2*hss(1,1)
     &    **2*hl**2*cb**2*u1*tx**(-1) + 2*hss(1,1)**2*hr**2*sb**2*t1*
     &    tx**(-1) - 2*hss(1,1)**2*hr**2*sb**2*u1*tx**(-1) )
      MMs = MMs + SCC(6,3)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*sb**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2
     &    *s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*s**(-1)*t1 - 8
     &    *hss(1,2)*hss(2,1)*pq*ssp*cb**2*s**(-1)*u1 + 8*hss(1,2)*hss(2
     &    ,1)*lq*ssz*cb**2*t1*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*ssz*
     &    cb**2*u1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*t1*
     &    sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*u1*sz**(-1) + 2*
     &    hss(1,2)*hss(2,1)*hl**2*cb**2*t1*tx**(-1) - 2*hss(1,2)*hss(2,
     &    1)*hl**2*cb**2*u1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*sb**2*
     &    t1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*sb**2*u1*tx**(-1) )
      MMs = MMs + SCC(6,4)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*
     &    hss(2,1)*lq*ssz*sb*cb*t1*sz**(-1) + 8*hss(1,1)*hss(2,1)*lq*
     &    ssz*sb*cb*u1*sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*t1*
     &    sz**(-1) - 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(1,1)*hss(2,1)*hl**2*sb*cb*t1*tx**(-1) + 2*hss(1,1)*hss(2,
     &    1)*hl**2*sb*cb*u1*tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*
     &    t1*tx**(-1) - 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*u1*tx**(-1) )
      MMs = MMs + SCC(6,5)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(2,1)*
     &    hss(2,2)*lq*ssz*sb*cb*t1*sz**(-1) + 8*hss(2,1)*hss(2,2)*lq*
     &    ssz*sb*cb*u1*sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*t1*
     &    sz**(-1) - 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(2,1)*hss(2,2)*hl**2*sb*cb*t1*tx**(-1) + 2*hss(2,1)*hss(2,
     &    2)*hl**2*sb*cb*u1*tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*
     &    t1*tx**(-1) - 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*u1*tx**(-1) )
      MMs = MMs + SCC(6,6)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*
     &    hss(1,2)*lq*ssz*sb*cb*t1*sz**(-1) + 8*hss(1,1)*hss(1,2)*lq*
     &    ssz*sb*cb*u1*sz**(-1) + 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*t1*
     &    sz**(-1) - 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(1,1)*hss(1,2)*hl**2*sb*cb*t1*tx**(-1) + 2*hss(1,1)*hss(1,
     &    2)*hl**2*sb*cb*u1*tx**(-1) + 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*
     &    t1*tx**(-1) - 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*u1*tx**(-1) )
      MMs = MMs + SCC(6,7)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,2)*
     &    hss(2,2)*lq*ssz*sb*cb*t1*sz**(-1) + 8*hss(1,2)*hss(2,2)*lq*
     &    ssz*sb*cb*u1*sz**(-1) + 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*t1*
     &    sz**(-1) - 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*u1*sz**(-1) - 2*
     &    hss(1,2)*hss(2,2)*hl**2*sb*cb*t1*tx**(-1) + 2*hss(1,2)*hss(2,
     &    2)*hl**2*sb*cb*u1*tx**(-1) + 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*
     &    t1*tx**(-1) - 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*u1*tx**(-1) )
      MMs = MMs + SCC(6,8)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*sb**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2
     &    *s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*s**(-1)*t1 - 8
     &    *hss(1,2)*hss(2,1)*pq*ssp*cb**2*s**(-1)*u1 + 8*hss(1,2)*hss(2
     &    ,1)*lq*ssz*sb**2*t1*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*ssz*
     &    sb**2*u1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*t1*
     &    sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*u1*sz**(-1) + 2*
     &    hss(1,2)*hss(2,1)*hl**2*sb**2*t1*tx**(-1) - 2*hss(1,2)*hss(2,
     &    1)*hl**2*sb**2*u1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*cb**2*
     &    t1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*cb**2*u1*tx**(-1) )
      MMs = MMs + SCC(6,9)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)**2*pq*
     &    ssp*sb**2*s**(-1)*t1 - 8*hss(2,2)**2*pq*ssp*sb**2*s**(-1)*u1
     &     + 8*hss(2,2)**2*pq*ssp*cb**2*s**(-1)*t1 - 8*hss(2,2)**2*pq*
     &    ssp*cb**2*s**(-1)*u1 + 8*hss(2,2)**2*lq*ssz*sb**2*t1*sz**(-1)
     &     - 8*hss(2,2)**2*lq*ssz*sb**2*u1*sz**(-1) + 8*hss(2,2)**2*rq*
     &    ssz*cb**2*t1*sz**(-1) - 8*hss(2,2)**2*rq*ssz*cb**2*u1*
     &    sz**(-1) + 2*hss(2,2)**2*hl**2*sb**2*t1*tx**(-1) - 2*hss(2,2)
     &    **2*hl**2*sb**2*u1*tx**(-1) + 2*hss(2,2)**2*hr**2*cb**2*t1*
     &    tx**(-1) - 2*hss(2,2)**2*hr**2*cb**2*u1*tx**(-1) )
ctp      print*, " SCC ",MMs
      MMs = MMs + SCD(3,1)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,1)**2*pq*
     &    ssp*sb**2*m1**2 + 8*hss(1,1)**2*pq*ssp*sb**2*mg**2*s**(-1)*t1
     &     - 8*hss(1,1)**2*pq*ssp*sb**2*mg**2*s**(-1)*u1 - 8*hss(1,1)**
     &    2*pq*ssp*sb**2*msb1**2*s**(-1)*t1 + 8*hss(1,1)**2*pq*ssp*
     &    sb**2*msb1**2*s**(-1)*u1 + 8*hss(1,1)**2*pq*ssp*sb**2*msb1**2
     &     - 8*hss(1,1)**2*pq*ssp*sb**2*mst1**2 + 8*hss(1,1)**2*pq*ssp*
     &    sb**2*t1 + 8*hss(1,1)**2*pq*ssp*cb**2*m1**2 + 8*hss(1,1)**2*
     &    pq*ssp*cb**2*mg**2*s**(-1)*t1 - 8*hss(1,1)**2*pq*ssp*cb**2*
     &    mg**2*s**(-1)*u1 - 8*hss(1,1)**2*pq*ssp*cb**2*msb1**2*s**(-1)
     &    *t1 + 8*hss(1,1)**2*pq*ssp*cb**2*msb1**2*s**(-1)*u1 + 8*hss(1
     &    ,1)**2*pq*ssp*cb**2*msb1**2 - 8*hss(1,1)**2*pq*ssp*cb**2*
     &    mst1**2 + 8*hss(1,1)**2*pq*ssp*cb**2*t1 + 8*hss(1,1)**2*lq*
     &    ssz*cb**2*m1**2*s*sz**(-1) + 8*hss(1,1)**2*lq*ssz*cb**2*mg**2
     &    *t1*sz**(-1) - 8*hss(1,1)**2*lq*ssz*cb**2*mg**2*u1*sz**(-1)
     &     + 8*hss(1,1)**2*lq*ssz*cb**2*msb1**2*s*sz**(-1) - 8*hss(1,1)
     &    **2*lq*ssz*cb**2*msb1**2*t1*sz**(-1) )
      MMs = MMs + SCD(3,1)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,1)**2*lq*
     &    ssz*cb**2*msb1**2*u1*sz**(-1) - 8*hss(1,1)**2*lq*ssz*cb**2*
     &    mst1**2*s*sz**(-1) + 8*hss(1,1)**2*lq*ssz*cb**2*s*t1*sz**(-1)
     &     + 8*hss(1,1)**2*rq*ssz*sb**2*m1**2*s*sz**(-1) + 8*hss(1,1)**
     &    2*rq*ssz*sb**2*mg**2*t1*sz**(-1) - 8*hss(1,1)**2*rq*ssz*sb**2
     &    *mg**2*u1*sz**(-1) + 8*hss(1,1)**2*rq*ssz*sb**2*msb1**2*s*
     &    sz**(-1) - 8*hss(1,1)**2*rq*ssz*sb**2*msb1**2*t1*sz**(-1) + 8
     &    *hss(1,1)**2*rq*ssz*sb**2*msb1**2*u1*sz**(-1) - 8*hss(1,1)**2
     &    *rq*ssz*sb**2*mst1**2*s*sz**(-1) + 8*hss(1,1)**2*rq*ssz*sb**2
     &    *s*t1*sz**(-1) + 8*hss(1,1)**2*hl*hr*sb*cb*mg*mt*s*tx**(-1)
     &     + 2*hss(1,1)**2*hl**2*cb**2*m1**2*s*tx**(-1) + 2*hss(1,1)**2
     &    *hl**2*cb**2*mg**2*t1*tx**(-1) - 2*hss(1,1)**2*hl**2*cb**2*
     &    mg**2*u1*tx**(-1) + 2*hss(1,1)**2*hl**2*cb**2*msb1**2*s*
     &    tx**(-1) - 2*hss(1,1)**2*hl**2*cb**2*msb1**2*t1*tx**(-1) + 2*
     &    hss(1,1)**2*hl**2*cb**2*msb1**2*u1*tx**(-1) - 2*hss(1,1)**2*
     &    hl**2*cb**2*mst1**2*s*tx**(-1) )
      MMs = MMs + SCD(3,1)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,1)**2*
     &    hl**2*cb**2*s*t1*tx**(-1) + 2*hss(1,1)**2*hr**2*sb**2*m1**2*s
     &    *tx**(-1) + 2*hss(1,1)**2*hr**2*sb**2*mg**2*t1*tx**(-1) - 2*
     &    hss(1,1)**2*hr**2*sb**2*mg**2*u1*tx**(-1) + 2*hss(1,1)**2*
     &    hr**2*sb**2*msb1**2*s*tx**(-1) - 2*hss(1,1)**2*hr**2*sb**2*
     &    msb1**2*t1*tx**(-1) + 2*hss(1,1)**2*hr**2*sb**2*msb1**2*u1*
     &    tx**(-1) - 2*hss(1,1)**2*hr**2*sb**2*mst1**2*s*tx**(-1) + 2*
     &    hss(1,1)**2*hr**2*sb**2*s*t1*tx**(-1) + 16*hss(1,1)**2*h1*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*s1**(-1) + 16*hss(1,1)**2*h2*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCD(3,2)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*sb**2*m1**2 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*
     &    mg**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*mg**2*
     &    s**(-1)*u1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb1**2*s**(-1)
     &    *t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb1**2*s**(-1)*u1 + 8
     &    *hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb1**2 - 8*hss(1,2)*hss(2,1)
     &    *pq*ssp*sb**2*mst2**2 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*t1
     &     + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*m1**2 + 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*cb**2*mg**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp
     &    *cb**2*mg**2*s**(-1)*u1 - 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*
     &    msb1**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*msb1**2
     &    *s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*msb1**2 - 8*
     &    hss(1,2)*hss(2,1)*pq*ssp*cb**2*mst2**2 + 8*hss(1,2)*hss(2,1)*
     &    pq*ssp*cb**2*t1 + 8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*m1**2*s*
     &    sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*mg**2*t1*sz**(-1)
     &     - 8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*mg**2*u1*sz**(-1) )
      MMs = MMs + SCD(3,2)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*lq*ssz*cb**2*msb1**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*
     &    ssz*cb**2*msb1**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*
     &    cb**2*msb1**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*
     &    mst2**2*s*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*cb**2*s*t1*
     &    sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*m1**2*s*sz**(-1)
     &     + 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*mg**2*t1*sz**(-1) - 8*
     &    hss(1,2)*hss(2,1)*rq*ssz*sb**2*mg**2*u1*sz**(-1) + 8*hss(1,2)
     &    *hss(2,1)*rq*ssz*sb**2*msb1**2*s*sz**(-1) - 8*hss(1,2)*hss(2,
     &    1)*rq*ssz*sb**2*msb1**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*
     &    ssz*sb**2*msb1**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*
     &    sb**2*mst2**2*s*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*sb**2*s
     &    *t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*hl*hr*sb*cb*mg*mt*s*
     &    tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*cb**2*m1**2*s*tx**(-1)
     &     + 2*hss(1,2)*hss(2,1)*hl**2*cb**2*mg**2*t1*tx**(-1) - 2*hss(
     &    1,2)*hss(2,1)*hl**2*cb**2*mg**2*u1*tx**(-1) )
      MMs = MMs + SCD(3,2)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,2)*hss(2,
     &    1)*hl**2*cb**2*msb1**2*s*tx**(-1) - 2*hss(1,2)*hss(2,1)*hl**2
     &    *cb**2*msb1**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*cb**2*
     &    msb1**2*u1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hl**2*cb**2*mst2**2
     &    *s*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*cb**2*s*t1*tx**(-1)
     &     + 2*hss(1,2)*hss(2,1)*hr**2*sb**2*m1**2*s*tx**(-1) + 2*hss(1
     &    ,2)*hss(2,1)*hr**2*sb**2*mg**2*t1*tx**(-1) - 2*hss(1,2)*hss(2
     &    ,1)*hr**2*sb**2*mg**2*u1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2
     &    *sb**2*msb1**2*s*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*sb**2*
     &    msb1**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*sb**2*msb1**2
     &    *u1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*sb**2*mst2**2*s*
     &    tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*sb**2*s*t1*tx**(-1) + 16
     &    *hss(1,2)*hss(2,1)*h1*lambda1*sb*cb*sqrt2**(-1)*mg*s*s1**(-1)
     &     + 16*hss(1,2)*hss(2,1)*h2*lambda2*sb*cb*sqrt2**(-1)*mg*s*
     &    s2**(-1) )
      MMs = MMs + SCD(3,3)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*
     &    hss(1,2)*lq*ssz*sb*cb*m1**2*s*sz**(-1) - 8*hss(1,1)*hss(1,2)*
     &    lq*ssz*sb*cb*mg**2*t1*sz**(-1) + 8*hss(1,1)*hss(1,2)*lq*ssz*
     &    sb*cb*mg**2*u1*sz**(-1) + 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*
     &    msb1**2*t1*sz**(-1) - 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*
     &    msb2**2*s*sz**(-1) - 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*msb2**2
     &    *u1*sz**(-1) + 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*mst1**2*s*
     &    sz**(-1) - 8*hss(1,1)*hss(1,2)*lq*ssz*sb*cb*s*t1*sz**(-1) + 8
     &    *hss(1,1)*hss(1,2)*rq*ssz*sb*cb*m1**2*s*sz**(-1) + 8*hss(1,1)
     &    *hss(1,2)*rq*ssz*sb*cb*mg**2*t1*sz**(-1) - 8*hss(1,1)*hss(1,2
     &    )*rq*ssz*sb*cb*mg**2*u1*sz**(-1) - 8*hss(1,1)*hss(1,2)*rq*ssz
     &    *sb*cb*msb1**2*t1*sz**(-1) + 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb
     &    *msb2**2*s*sz**(-1) + 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*
     &    msb2**2*u1*sz**(-1) - 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*
     &    mst1**2*s*sz**(-1) + 8*hss(1,1)*hss(1,2)*rq*ssz*sb*cb*s*t1*
     &    sz**(-1) )
      MMs = MMs + SCD(3,3)*Nc*Cf*Pi*alphas*prefac * (  - 4*hss(1,1)*
     &    hss(1,2)*hl*hr*sb**2*mg*mt*s*tx**(-1) + 4*hss(1,1)*hss(1,2)*
     &    hl*hr*cb**2*mg*mt*s*tx**(-1) - 2*hss(1,1)*hss(1,2)*hl**2*sb*
     &    cb*m1**2*s*tx**(-1) - 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*mg**2*
     &    t1*tx**(-1) + 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*mg**2*u1*
     &    tx**(-1) + 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*msb1**2*t1*
     &    tx**(-1) - 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*msb2**2*s*tx**(-1)
     &     - 2*hss(1,1)*hss(1,2)*hl**2*sb*cb*msb2**2*u1*tx**(-1) + 2*
     &    hss(1,1)*hss(1,2)*hl**2*sb*cb*mst1**2*s*tx**(-1) - 2*hss(1,1)
     &    *hss(1,2)*hl**2*sb*cb*s*t1*tx**(-1) + 2*hss(1,1)*hss(1,2)*
     &    hr**2*sb*cb*m1**2*s*tx**(-1) + 2*hss(1,1)*hss(1,2)*hr**2*sb*
     &    cb*mg**2*t1*tx**(-1) - 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*mg**2*
     &    u1*tx**(-1) - 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*msb1**2*t1*
     &    tx**(-1) + 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*msb2**2*s*tx**(-1)
     &     + 2*hss(1,1)*hss(1,2)*hr**2*sb*cb*msb2**2*u1*tx**(-1) - 2*
     &    hss(1,1)*hss(1,2)*hr**2*sb*cb*mst1**2*s*tx**(-1) )
      MMs = MMs + SCD(3,3)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,1)*hss(1,
     &    2)*hr**2*sb*cb*s*t1*tx**(-1) - 8*hss(1,1)*hss(1,2)*h1*lambda1
     &    *sb**2*sqrt2**(-1)*mg*s*s1**(-1) + 8*hss(1,1)*hss(1,2)*h1*
     &    lambda1*cb**2*sqrt2**(-1)*mg*s*s1**(-1) - 8*hss(1,1)*hss(1,2)
     &    *h2*lambda2*sb**2*sqrt2**(-1)*mg*s*s2**(-1) + 8*hss(1,1)*hss(
     &    1,2)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCD(3,4)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,2)*
     &    hss(2,2)*lq*ssz*sb*cb*m1**2*s*sz**(-1) - 8*hss(1,2)*hss(2,2)*
     &    lq*ssz*sb*cb*mg**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,2)*lq*ssz*
     &    sb*cb*mg**2*u1*sz**(-1) + 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*
     &    msb1**2*t1*sz**(-1) - 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*
     &    msb2**2*s*sz**(-1) - 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*msb2**2
     &    *u1*sz**(-1) + 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*mst2**2*s*
     &    sz**(-1) - 8*hss(1,2)*hss(2,2)*lq*ssz*sb*cb*s*t1*sz**(-1) + 8
     &    *hss(1,2)*hss(2,2)*rq*ssz*sb*cb*m1**2*s*sz**(-1) + 8*hss(1,2)
     &    *hss(2,2)*rq*ssz*sb*cb*mg**2*t1*sz**(-1) - 8*hss(1,2)*hss(2,2
     &    )*rq*ssz*sb*cb*mg**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,2)*rq*ssz
     &    *sb*cb*msb1**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb
     &    *msb2**2*s*sz**(-1) + 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*
     &    msb2**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*
     &    mst2**2*s*sz**(-1) + 8*hss(1,2)*hss(2,2)*rq*ssz*sb*cb*s*t1*
     &    sz**(-1) )
      MMs = MMs + SCD(3,4)*Nc*Cf*Pi*alphas*prefac * (  - 4*hss(1,2)*
     &    hss(2,2)*hl*hr*sb**2*mg*mt*s*tx**(-1) + 4*hss(1,2)*hss(2,2)*
     &    hl*hr*cb**2*mg*mt*s*tx**(-1) - 2*hss(1,2)*hss(2,2)*hl**2*sb*
     &    cb*m1**2*s*tx**(-1) - 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*mg**2*
     &    t1*tx**(-1) + 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*mg**2*u1*
     &    tx**(-1) + 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*msb1**2*t1*
     &    tx**(-1) - 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*msb2**2*s*tx**(-1)
     &     - 2*hss(1,2)*hss(2,2)*hl**2*sb*cb*msb2**2*u1*tx**(-1) + 2*
     &    hss(1,2)*hss(2,2)*hl**2*sb*cb*mst2**2*s*tx**(-1) - 2*hss(1,2)
     &    *hss(2,2)*hl**2*sb*cb*s*t1*tx**(-1) + 2*hss(1,2)*hss(2,2)*
     &    hr**2*sb*cb*m1**2*s*tx**(-1) + 2*hss(1,2)*hss(2,2)*hr**2*sb*
     &    cb*mg**2*t1*tx**(-1) - 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*mg**2*
     &    u1*tx**(-1) - 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*msb1**2*t1*
     &    tx**(-1) + 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*msb2**2*s*tx**(-1)
     &     + 2*hss(1,2)*hss(2,2)*hr**2*sb*cb*msb2**2*u1*tx**(-1) - 2*
     &    hss(1,2)*hss(2,2)*hr**2*sb*cb*mst2**2*s*tx**(-1) )
      MMs = MMs + SCD(3,4)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,2)*hss(2,
     &    2)*hr**2*sb*cb*s*t1*tx**(-1) - 8*hss(1,2)*hss(2,2)*h1*lambda1
     &    *sb**2*sqrt2**(-1)*mg*s*s1**(-1) + 8*hss(1,2)*hss(2,2)*h1*
     &    lambda1*cb**2*sqrt2**(-1)*mg*s*s1**(-1) - 8*hss(1,2)*hss(2,2)
     &    *h2*lambda2*sb**2*sqrt2**(-1)*mg*s*s2**(-1) + 8*hss(1,2)*hss(
     &    2,2)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCD(3,5)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(1,1)*
     &    hss(2,1)*lq*ssz*sb*cb*m1**2*s*sz**(-1) - 8*hss(1,1)*hss(2,1)*
     &    lq*ssz*sb*cb*mg**2*t1*sz**(-1) + 8*hss(1,1)*hss(2,1)*lq*ssz*
     &    sb*cb*mg**2*u1*sz**(-1) - 8*hss(1,1)*hss(2,1)*lq*ssz*sb*cb*
     &    msb1**2*s*sz**(-1) - 8*hss(1,1)*hss(2,1)*lq*ssz*sb*cb*msb1**2
     &    *u1*sz**(-1) + 8*hss(1,1)*hss(2,1)*lq*ssz*sb*cb*msb2**2*t1*
     &    sz**(-1) + 8*hss(1,1)*hss(2,1)*lq*ssz*sb*cb*mst1**2*s*
     &    sz**(-1) - 8*hss(1,1)*hss(2,1)*lq*ssz*sb*cb*s*t1*sz**(-1) + 8
     &    *hss(1,1)*hss(2,1)*rq*ssz*sb*cb*m1**2*s*sz**(-1) + 8*hss(1,1)
     &    *hss(2,1)*rq*ssz*sb*cb*mg**2*t1*sz**(-1) - 8*hss(1,1)*hss(2,1
     &    )*rq*ssz*sb*cb*mg**2*u1*sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz
     &    *sb*cb*msb1**2*s*sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*
     &    msb1**2*u1*sz**(-1) - 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*
     &    msb2**2*t1*sz**(-1) - 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*
     &    mst1**2*s*sz**(-1) + 8*hss(1,1)*hss(2,1)*rq*ssz*sb*cb*s*t1*
     &    sz**(-1) )
      MMs = MMs + SCD(3,5)*Nc*Cf*Pi*alphas*prefac * (  - 4*hss(1,1)*
     &    hss(2,1)*hl*hr*sb**2*mg*mt*s*tx**(-1) + 4*hss(1,1)*hss(2,1)*
     &    hl*hr*cb**2*mg*mt*s*tx**(-1) - 2*hss(1,1)*hss(2,1)*hl**2*sb*
     &    cb*m1**2*s*tx**(-1) - 2*hss(1,1)*hss(2,1)*hl**2*sb*cb*mg**2*
     &    t1*tx**(-1) + 2*hss(1,1)*hss(2,1)*hl**2*sb*cb*mg**2*u1*
     &    tx**(-1) - 2*hss(1,1)*hss(2,1)*hl**2*sb*cb*msb1**2*s*tx**(-1)
     &     - 2*hss(1,1)*hss(2,1)*hl**2*sb*cb*msb1**2*u1*tx**(-1) + 2*
     &    hss(1,1)*hss(2,1)*hl**2*sb*cb*msb2**2*t1*tx**(-1) + 2*hss(1,1
     &    )*hss(2,1)*hl**2*sb*cb*mst1**2*s*tx**(-1) - 2*hss(1,1)*hss(2,
     &    1)*hl**2*sb*cb*s*t1*tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*
     &    cb*m1**2*s*tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*mg**2*
     &    t1*tx**(-1) - 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*mg**2*u1*
     &    tx**(-1) + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*msb1**2*s*tx**(-1)
     &     + 2*hss(1,1)*hss(2,1)*hr**2*sb*cb*msb1**2*u1*tx**(-1) - 2*
     &    hss(1,1)*hss(2,1)*hr**2*sb*cb*msb2**2*t1*tx**(-1) - 2*hss(1,1
     &    )*hss(2,1)*hr**2*sb*cb*mst1**2*s*tx**(-1) )
      MMs = MMs + SCD(3,5)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,1)*hss(2,
     &    1)*hr**2*sb*cb*s*t1*tx**(-1) - 8*hss(1,1)*hss(2,1)*h1*lambda1
     &    *sb**2*sqrt2**(-1)*mg*s*s1**(-1) + 8*hss(1,1)*hss(2,1)*h1*
     &    lambda1*cb**2*sqrt2**(-1)*mg*s*s1**(-1) - 8*hss(1,1)*hss(2,1)
     &    *h2*lambda2*sb**2*sqrt2**(-1)*mg*s*s2**(-1) + 8*hss(1,1)*hss(
     &    2,1)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCD(3,6)*Nc*Cf*Pi*alphas*prefac * (  - 8*hss(2,1)*
     &    hss(2,2)*lq*ssz*sb*cb*m1**2*s*sz**(-1) - 8*hss(2,1)*hss(2,2)*
     &    lq*ssz*sb*cb*mg**2*t1*sz**(-1) + 8*hss(2,1)*hss(2,2)*lq*ssz*
     &    sb*cb*mg**2*u1*sz**(-1) - 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*
     &    msb1**2*s*sz**(-1) - 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*msb1**2
     &    *u1*sz**(-1) + 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*msb2**2*t1*
     &    sz**(-1) + 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*mst2**2*s*
     &    sz**(-1) - 8*hss(2,1)*hss(2,2)*lq*ssz*sb*cb*s*t1*sz**(-1) + 8
     &    *hss(2,1)*hss(2,2)*rq*ssz*sb*cb*m1**2*s*sz**(-1) + 8*hss(2,1)
     &    *hss(2,2)*rq*ssz*sb*cb*mg**2*t1*sz**(-1) - 8*hss(2,1)*hss(2,2
     &    )*rq*ssz*sb*cb*mg**2*u1*sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz
     &    *sb*cb*msb1**2*s*sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*
     &    msb1**2*u1*sz**(-1) - 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*
     &    msb2**2*t1*sz**(-1) - 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*
     &    mst2**2*s*sz**(-1) + 8*hss(2,1)*hss(2,2)*rq*ssz*sb*cb*s*t1*
     &    sz**(-1) )
      MMs = MMs + SCD(3,6)*Nc*Cf*Pi*alphas*prefac * (  - 4*hss(2,1)*
     &    hss(2,2)*hl*hr*sb**2*mg*mt*s*tx**(-1) + 4*hss(2,1)*hss(2,2)*
     &    hl*hr*cb**2*mg*mt*s*tx**(-1) - 2*hss(2,1)*hss(2,2)*hl**2*sb*
     &    cb*m1**2*s*tx**(-1) - 2*hss(2,1)*hss(2,2)*hl**2*sb*cb*mg**2*
     &    t1*tx**(-1) + 2*hss(2,1)*hss(2,2)*hl**2*sb*cb*mg**2*u1*
     &    tx**(-1) - 2*hss(2,1)*hss(2,2)*hl**2*sb*cb*msb1**2*s*tx**(-1)
     &     - 2*hss(2,1)*hss(2,2)*hl**2*sb*cb*msb1**2*u1*tx**(-1) + 2*
     &    hss(2,1)*hss(2,2)*hl**2*sb*cb*msb2**2*t1*tx**(-1) + 2*hss(2,1
     &    )*hss(2,2)*hl**2*sb*cb*mst2**2*s*tx**(-1) - 2*hss(2,1)*hss(2,
     &    2)*hl**2*sb*cb*s*t1*tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*
     &    cb*m1**2*s*tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*mg**2*
     &    t1*tx**(-1) - 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*mg**2*u1*
     &    tx**(-1) + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*msb1**2*s*tx**(-1)
     &     + 2*hss(2,1)*hss(2,2)*hr**2*sb*cb*msb1**2*u1*tx**(-1) - 2*
     &    hss(2,1)*hss(2,2)*hr**2*sb*cb*msb2**2*t1*tx**(-1) - 2*hss(2,1
     &    )*hss(2,2)*hr**2*sb*cb*mst2**2*s*tx**(-1) )
      MMs = MMs + SCD(3,6)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(2,1)*hss(2,
     &    2)*hr**2*sb*cb*s*t1*tx**(-1) - 8*hss(2,1)*hss(2,2)*h1*lambda1
     &    *sb**2*sqrt2**(-1)*mg*s*s1**(-1) + 8*hss(2,1)*hss(2,2)*h1*
     &    lambda1*cb**2*sqrt2**(-1)*mg*s*s1**(-1) - 8*hss(2,1)*hss(2,2)
     &    *h2*lambda2*sb**2*sqrt2**(-1)*mg*s*s2**(-1) + 8*hss(2,1)*hss(
     &    2,2)*h2*lambda2*cb**2*sqrt2**(-1)*mg*s*s2**(-1) )
      MMs = MMs + SCD(3,7)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*sb**2*m1**2 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*
     &    mg**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*mg**2*
     &    s**(-1)*u1 - 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb2**2*s**(-1)
     &    *t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb2**2*s**(-1)*u1 + 8
     &    *hss(1,2)*hss(2,1)*pq*ssp*sb**2*msb2**2 - 8*hss(1,2)*hss(2,1)
     &    *pq*ssp*sb**2*mst1**2 + 8*hss(1,2)*hss(2,1)*pq*ssp*sb**2*t1
     &     + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*m1**2 + 8*hss(1,2)*hss(2,
     &    1)*pq*ssp*cb**2*mg**2*s**(-1)*t1 - 8*hss(1,2)*hss(2,1)*pq*ssp
     &    *cb**2*mg**2*s**(-1)*u1 - 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*
     &    msb2**2*s**(-1)*t1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*msb2**2
     &    *s**(-1)*u1 + 8*hss(1,2)*hss(2,1)*pq*ssp*cb**2*msb2**2 - 8*
     &    hss(1,2)*hss(2,1)*pq*ssp*cb**2*mst1**2 + 8*hss(1,2)*hss(2,1)*
     &    pq*ssp*cb**2*t1 + 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*m1**2*s*
     &    sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*mg**2*t1*sz**(-1)
     &     - 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*mg**2*u1*sz**(-1) )
      MMs = MMs + SCD(3,7)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(1,2)*hss(2,
     &    1)*lq*ssz*sb**2*msb2**2*s*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*
     &    ssz*sb**2*msb2**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*
     &    sb**2*msb2**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*
     &    mst1**2*s*sz**(-1) + 8*hss(1,2)*hss(2,1)*lq*ssz*sb**2*s*t1*
     &    sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*m1**2*s*sz**(-1)
     &     + 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*mg**2*t1*sz**(-1) - 8*
     &    hss(1,2)*hss(2,1)*rq*ssz*cb**2*mg**2*u1*sz**(-1) + 8*hss(1,2)
     &    *hss(2,1)*rq*ssz*cb**2*msb2**2*s*sz**(-1) - 8*hss(1,2)*hss(2,
     &    1)*rq*ssz*cb**2*msb2**2*t1*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*
     &    ssz*cb**2*msb2**2*u1*sz**(-1) - 8*hss(1,2)*hss(2,1)*rq*ssz*
     &    cb**2*mst1**2*s*sz**(-1) + 8*hss(1,2)*hss(2,1)*rq*ssz*cb**2*s
     &    *t1*sz**(-1) - 8*hss(1,2)*hss(2,1)*hl*hr*sb*cb*mg*mt*s*
     &    tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*sb**2*m1**2*s*tx**(-1)
     &     + 2*hss(1,2)*hss(2,1)*hl**2*sb**2*mg**2*t1*tx**(-1) - 2*hss(
     &    1,2)*hss(2,1)*hl**2*sb**2*mg**2*u1*tx**(-1) )
      MMs = MMs + SCD(3,7)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(1,2)*hss(2,
     &    1)*hl**2*sb**2*msb2**2*s*tx**(-1) - 2*hss(1,2)*hss(2,1)*hl**2
     &    *sb**2*msb2**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*sb**2*
     &    msb2**2*u1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hl**2*sb**2*mst1**2
     &    *s*tx**(-1) + 2*hss(1,2)*hss(2,1)*hl**2*sb**2*s*t1*tx**(-1)
     &     + 2*hss(1,2)*hss(2,1)*hr**2*cb**2*m1**2*s*tx**(-1) + 2*hss(1
     &    ,2)*hss(2,1)*hr**2*cb**2*mg**2*t1*tx**(-1) - 2*hss(1,2)*hss(2
     &    ,1)*hr**2*cb**2*mg**2*u1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2
     &    *cb**2*msb2**2*s*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*cb**2*
     &    msb2**2*t1*tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*cb**2*msb2**2
     &    *u1*tx**(-1) - 2*hss(1,2)*hss(2,1)*hr**2*cb**2*mst1**2*s*
     &    tx**(-1) + 2*hss(1,2)*hss(2,1)*hr**2*cb**2*s*t1*tx**(-1) - 16
     &    *hss(1,2)*hss(2,1)*h1*lambda1*sb*cb*sqrt2**(-1)*mg*s*s1**(-1)
     &     - 16*hss(1,2)*hss(2,1)*h2*lambda2*sb*cb*sqrt2**(-1)*mg*s*
     &    s2**(-1) )
      MMs = MMs + SCD(3,8)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)**2*pq*
     &    ssp*sb**2*m1**2 + 8*hss(2,2)**2*pq*ssp*sb**2*mg**2*s**(-1)*t1
     &     - 8*hss(2,2)**2*pq*ssp*sb**2*mg**2*s**(-1)*u1 - 8*hss(2,2)**
     &    2*pq*ssp*sb**2*msb2**2*s**(-1)*t1 + 8*hss(2,2)**2*pq*ssp*
     &    sb**2*msb2**2*s**(-1)*u1 + 8*hss(2,2)**2*pq*ssp*sb**2*msb2**2
     &     - 8*hss(2,2)**2*pq*ssp*sb**2*mst2**2 + 8*hss(2,2)**2*pq*ssp*
     &    sb**2*t1 + 8*hss(2,2)**2*pq*ssp*cb**2*m1**2 + 8*hss(2,2)**2*
     &    pq*ssp*cb**2*mg**2*s**(-1)*t1 - 8*hss(2,2)**2*pq*ssp*cb**2*
     &    mg**2*s**(-1)*u1 - 8*hss(2,2)**2*pq*ssp*cb**2*msb2**2*s**(-1)
     &    *t1 + 8*hss(2,2)**2*pq*ssp*cb**2*msb2**2*s**(-1)*u1 + 8*hss(2
     &    ,2)**2*pq*ssp*cb**2*msb2**2 - 8*hss(2,2)**2*pq*ssp*cb**2*
     &    mst2**2 + 8*hss(2,2)**2*pq*ssp*cb**2*t1 + 8*hss(2,2)**2*lq*
     &    ssz*sb**2*m1**2*s*sz**(-1) + 8*hss(2,2)**2*lq*ssz*sb**2*mg**2
     &    *t1*sz**(-1) - 8*hss(2,2)**2*lq*ssz*sb**2*mg**2*u1*sz**(-1)
     &     + 8*hss(2,2)**2*lq*ssz*sb**2*msb2**2*s*sz**(-1) - 8*hss(2,2)
     &    **2*lq*ssz*sb**2*msb2**2*t1*sz**(-1) )
      MMs = MMs + SCD(3,8)*Nc*Cf*Pi*alphas*prefac * ( 8*hss(2,2)**2*lq*
     &    ssz*sb**2*msb2**2*u1*sz**(-1) - 8*hss(2,2)**2*lq*ssz*sb**2*
     &    mst2**2*s*sz**(-1) + 8*hss(2,2)**2*lq*ssz*sb**2*s*t1*sz**(-1)
     &     + 8*hss(2,2)**2*rq*ssz*cb**2*m1**2*s*sz**(-1) + 8*hss(2,2)**
     &    2*rq*ssz*cb**2*mg**2*t1*sz**(-1) - 8*hss(2,2)**2*rq*ssz*cb**2
     &    *mg**2*u1*sz**(-1) + 8*hss(2,2)**2*rq*ssz*cb**2*msb2**2*s*
     &    sz**(-1) - 8*hss(2,2)**2*rq*ssz*cb**2*msb2**2*t1*sz**(-1) + 8
     &    *hss(2,2)**2*rq*ssz*cb**2*msb2**2*u1*sz**(-1) - 8*hss(2,2)**2
     &    *rq*ssz*cb**2*mst2**2*s*sz**(-1) + 8*hss(2,2)**2*rq*ssz*cb**2
     &    *s*t1*sz**(-1) - 8*hss(2,2)**2*hl*hr*sb*cb*mg*mt*s*tx**(-1)
     &     + 2*hss(2,2)**2*hl**2*sb**2*m1**2*s*tx**(-1) + 2*hss(2,2)**2
     &    *hl**2*sb**2*mg**2*t1*tx**(-1) - 2*hss(2,2)**2*hl**2*sb**2*
     &    mg**2*u1*tx**(-1) + 2*hss(2,2)**2*hl**2*sb**2*msb2**2*s*
     &    tx**(-1) - 2*hss(2,2)**2*hl**2*sb**2*msb2**2*t1*tx**(-1) + 2*
     &    hss(2,2)**2*hl**2*sb**2*msb2**2*u1*tx**(-1) - 2*hss(2,2)**2*
     &    hl**2*sb**2*mst2**2*s*tx**(-1) )
      MMs = MMs + SCD(3,8)*Nc*Cf*Pi*alphas*prefac * ( 2*hss(2,2)**2*
     &    hl**2*sb**2*s*t1*tx**(-1) + 2*hss(2,2)**2*hr**2*cb**2*m1**2*s
     &    *tx**(-1) + 2*hss(2,2)**2*hr**2*cb**2*mg**2*t1*tx**(-1) - 2*
     &    hss(2,2)**2*hr**2*cb**2*mg**2*u1*tx**(-1) + 2*hss(2,2)**2*
     &    hr**2*cb**2*msb2**2*s*tx**(-1) - 2*hss(2,2)**2*hr**2*cb**2*
     &    msb2**2*t1*tx**(-1) + 2*hss(2,2)**2*hr**2*cb**2*msb2**2*u1*
     &    tx**(-1) - 2*hss(2,2)**2*hr**2*cb**2*mst2**2*s*tx**(-1) + 2*
     &    hss(2,2)**2*hr**2*cb**2*s*t1*tx**(-1) - 16*hss(2,2)**2*h1*
     &    lambda1*sb*cb*sqrt2**(-1)*mg*s*s1**(-1) - 16*hss(2,2)**2*h2*
     &    lambda2*sb*cb*sqrt2**(-1)*mg*s*s2**(-1) )
ctp      print*, " SCD ",MMs

c               the phase space except for 1/s**2 
      HH_QBS = MMs / ( 16.D0 * Pi )

c               the averaging factors
      HH_QBS = HH_QBS /4.D0 /Nc**2

c               the prefactor for the scaling functions 
      HH_QBS = HH_QBS * (m1+m2)**2/4.D0

      end
