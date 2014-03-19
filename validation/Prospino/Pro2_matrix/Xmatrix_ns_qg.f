cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NS_QGH(MASSIN,C)                                                  c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = t2                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = m2                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(4-5,8,10,13-30) not needed                              c
c                                                                      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function NS_QGH(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,gs
     &          ,QF,s,s4,m1,m2,ms,t1,tj,u1,uj,u1hat
     &          ,sptj,spuj,hardfac,mjs2,coeff
     &          ,gqH,NSANG(0:9,0:9,-2:2,-2:2)
      complex*16 Ar(4),Arc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      tj  = massin(2)
      s4  = massin(3)
      m1  = massin(6)
      m2  = massin(7)
      ms  = massin(11)

c               real kinematics built in
      t1 = tj + m2**2 - m1**2

      u1 = s4 - s - tj 
      uj = u1 + m1**2 - m2**2 

      hardfac = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2

c               some denominators following wim's notes 
      sptj  = s + tj 
      spuj  = s + uj
      u1hat = -( tj*uj - mjs2*s )/(s+uj)

c               the angular functions 
      call ANGULAR_ARRAY_NS_QG(massin,NSANG)

c               the factorization/renormalization scale 
      QF = massin(13) 

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2 
      coeff = gs**4 * real( Ar(1)*Arc(1) )/2.D0 * Sn**2/2.D0

c               wim's form output, replace mjs**2 -> mjs2 
      gqH = s4**(-1)*pi*CA**(-2)*coeff * (  - 4 - 4*s**(-1)*u1 )
      gqH = gqH + s4**(-1)*pi*CA**(-1)*CF*coeff * (  - 8 - 8*s**(-1)*u1
     +     )
      gqH = gqH + s4**(-1)*pi*coeff * ( 4 + 4*s**(-1)*u1 )
      gqH = gqH + pi*u1hat**(-2)*CA**(-1)*CF*spuj**(-2)*coeff * ( 8*
     +    ms**2*mjs2**(-1)*tj*u1**2 - 8*ms**2*mjs2*u1 - 8*ms**2*tj*u1
     +     + 8*ms**2*u1**2 )
      gqH = gqH + pi*u1hat**(-1)*CA**(-1)*CF*spuj**(-1)*coeff * ( 16*
     +    ms**2*mjs2**(-1)*u1 - 8*ms**2 + 8*mjs2*tj**(-1)*u1 + 8*u1 )
      gqH = gqH + pi*CA**(-2)*sptj**(-1)*coeff * ( 2*s**(-1)*u1 )
      gqH = gqH + pi*CA**(-2)*coeff * (  - 2*s**(-1)*tj*u1**(-1) + 2*
     +    s**(-1) - 2*u1**(-1) )
      gqH = gqH + pi*CA**(-1)*CF*sptj**(-1)*coeff * ( 4*s**(-1)*u1 )
      gqH = gqH + pi*CA**(-1)*CF*coeff * ( 8*ms**2*mjs2**(-1)*tj**(-1)
     +     - 8*ms**2*tj**(-1)*t1**(-1) + 4*mjs2*s**(-1)*t1**(-2)*u1
     +     + 4*mjs2*t1**(-2) - 4*mjs2**2*s**(-1)*t1**(-2) - 8*s**(-1)*
     +    tj**(-1)*u1 - 4*s**(-1)*tj*u1**(-1) + 4*s**(-1)*t1**(-1)*u1
     +     + 4*t1**(-1) - 4*u1**(-1) )
      gqH = gqH + pi*sptj**(-1)*coeff * (  - 2*s**(-1)*u1 )
      gqH = gqH + pi*coeff * (  - 16*ms**2*mjs2*s*u1**(-4) - 16*ms**2
     +    *mjs2*tj*u1**(-4) - 16*ms**2*mjs2*u1**(-3) + 16*mjs2*
     +    s**(-1)*tj*u1**(-2) + 16*mjs2*s**(-1)*tj**2*u1**(-3) + 16*
     +    mjs2*tj*u1**(-3) - 16*mjs2**2*s**(-1)*tj*u1**(-3) - 16*
     +    mjs2**2*s**(-1)*tj**2*u1**(-4) - 16*mjs2**2*s*u1**(-4) - 32*
     +    mjs2**2*tj*u1**(-4) - 16*mjs2**2*u1**(-3) + 2*s**(-1)*tj*
     +    u1**(-1) - 2*s**(-1) + 2*u1**(-1) )
      gqH = gqH + NSANG(0,0,0,0)*s4**(-2)*CA**(-1)*CF*hardfac**(-1)*
     + coeff * ( 8*ms**2*mjs2*s**(-1) + 8*ms**2*mjs2*u1**(-1) - 8*
     +    ms**2*mjs2**2*s**(-1)*u1**(-1) - 4*ms**2*s**(-1)*u1 - 4*ms**2
     +     + 8*ms**4*mjs2*u1**(-2) )
      gqH = gqH + NSANG(0,0,0,0)*s4**(-1)*CA**(-1)*CF*hardfac**(-1)*
     + coeff * ( 2*ms**2*mjs2*s**(-1)*u1**(-1) + 8*ms**2*mjs2*
     +    u1**(-2) - 2*ms**2*s**(-1) - 2*ms**2*u1**(-1) - 2*mjs2*
     +    s**(-1) - 2*mjs2**2*s**(-1)*u1**(-1) - 2*tj*u1**(-1) )
      gqH = gqH + NSANG(0,0,0,0)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 2*
     +    mjs2*s**(-1)*u1**(-1) + 4*mjs2*u1**(-2) )
      gqH = gqH + NSANG(1,0,-2,0)*hardfac**(-1)*coeff * ( 8*mjs2*s*tj
     +    *u1**(-2) + 4*mjs2*s**2*u1**(-2) + 4*mjs2*tj**2*u1**(-2)
     +     )
      gqH = gqH + NSANG(1,0,-1,0)*s4**(-1)*hardfac**(-1)*coeff * ( 4*
     +    ms**2*mjs2*s**(-1) + 4*ms**2*mjs2*u1**(-1) - 4*ms**2*
     +    mjs2**2*s**(-1)*u1**(-1) + mjs2*s**(-1)*u1 )
      gqH = gqH + NSANG(1,0,-1,0)*sptj**(-1)*hardfac**(-1)*coeff * ( 4*
     +    ms**2*mjs2*u1**(-1) + 4*mjs2 - 2*u1 )
      gqH = gqH + NSANG(1,0,-1,0)*hardfac**(-1)*coeff * (  - 3 + 8*
     +    ms**2*mjs2*u1**(-2) - mjs2*s**(-1)*tj*u1**(-1) + mjs2*
     +    s**(-1) + 4*mjs2*s*u1**(-2) + 4*mjs2*tj*u1**(-2) + 9*
     +    mjs2*u1**(-1) + s**(-1)*tj + s**(-1)*tj**2*u1**(-1) - 2*s*
     +    u1**(-1) - tj*u1**(-1) )
      gqH = gqH + NSANG(1,2,-1,-1)*sptj**(-1)*hardfac**(-1)*coeff * ( 
     +     - 4*ms**2*mjs2*tj*u1**(-1) - 2*mjs2*tj + 2*mjs2*u1 - 2
     +    *mjs2**2 + tj*u1 - u1**2 )
      gqH = gqH + NSANG(1,2,-1,-1)*hardfac**(-1)*coeff * ( 4*ms**2*
     +    mjs2*u1**(-1) + 2*mjs2*tj*u1**(-1) + 2*mjs2 + 2*mjs2**2*
     +    u1**(-1) + s*tj*u1**(-1) + tj**2*u1**(-1) - u1 )
      gqH = gqH + NSANG(1,5,-2,-2)*hardfac**(-1)*coeff * ( 4*mjs2*
     +    s**2 )
      gqH = gqH + NSANG(1,5,-2,-1)*hardfac**(-1)*coeff * ( 8*mjs2*s*
     +    tj*u1**(-1) + 8*mjs2*s**2*u1**(-1) )
      gqH = gqH + NSANG(1,5,-1,-2)*hardfac**(-1)*coeff * ( 4*mjs2*s )
      gqH = gqH + NSANG(1,5,-1,-1)*s4**(-1)*hardfac**(-1)*coeff * ( 4*
     +    ms**2*mjs2*tj*u1**(-1) + 4*ms**2*mjs2 + 2*mjs2*tj + 2*
     +    mjs2**2 - tj*u1 )
      gqH = gqH + NSANG(1,5,-1,-1)*sptj**(-1)*hardfac**(-1)*coeff * ( 4
     +    *ms**2*mjs2*tj*u1**(-1) + 2*mjs2*tj + 2*mjs2*u1 - 2*
     +    mjs2**2 - tj*u1 - u1**2 )
      gqH = gqH + NSANG(1,5,-1,-1)*hardfac**(-1)*coeff * (  - 8*ms**2*
     +    mjs2*u1**(-1) + 8*mjs2*s*u1**(-1) + 4*mjs2*tj*u1**(-1)
     +     + 8*mjs2 - 16*mjs2**2*u1**(-1) - 4*s*tj*u1**(-1) - 2*s - 4*
     +    tj - 4*tj**2*u1**(-1) - 3*u1 )
      gqH = gqH + NSANG(1,7,-1,1)*s4**(-1)*hardfac**(-1)*coeff * ( 4*
     +    ms**2*mjs2*s**(-1)*u1**(-1) + 2*mjs2*s**(-1) - s**(-1)*u1
     +     )
      gqH = gqH + NSANG(1,7,-1,1)*hardfac**(-1)*coeff * ( 2*mjs2*
     +    s**(-1)*u1**(-1) - s**(-1)*tj*u1**(-1) - s**(-1) - u1**(-1) )
      gqH = gqH + NSANG(2,0,-1,0)*s4**(-1)*CA**(-2)*hardfac**(-1)*coeff
     +  * (  - ms**2*mjs2*s**(-1) - ms**2*mjs2*u1**(-1) + ms**2*
     +    s**(-1)*u1 - ms**2*tj*u1**(-1) + ms**2 + mjs2*s**(-1)*u1 - 
     +    2*mjs2*tj*u1**(-1) + mjs2 - mjs2**2*s**(-1) - mjs2**2*
     +    u1**(-1) - tj**2*u1**(-1) )
      gqH = gqH + NSANG(2,0,-1,0)*CA**(-2)*hardfac**(-1)*coeff * ( 
     +    ms**2*u1**(-1) + mjs2*u1**(-1) + tj*u1**(-1) )
      gqH = gqH + NSANG(2,5,-1,-2)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +     - 2*mjs2*s )
      gqH = gqH + NSANG(2,5,-1,-1)*s4**(-1)*CA**(-2)*hardfac**(-1)*
     + coeff * (  - ms**2*mjs2*tj*u1**(-1) - ms**2*mjs2 - ms**2*tj
     +     - ms**2*tj**2*u1**(-1) - mjs2*tj - 4*mjs2*tj**2*u1**(-1)
     +   - 5*mjs2**2*tj*u1**(-1) - mjs2**2 - 2*mjs2**3*u1**(-1) - tj**3*
     +    u1**(-1) )
      gqH = gqH + NSANG(2,5,-1,-1)*CA**(-2)*hardfac**(-1)*coeff * ( 
     +    ms**2*mjs2*u1**(-1) + ms**2*tj*u1**(-1) + mjs2*tj*
     +    u1**(-1) + mjs2**2*u1**(-1) )
      gqH = gqH + NSANG(2,5,-1,-1)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +    2*ms**2*mjs2*s*tj**(-1)*u1**(-1) - 4*ms**2*s*tj**(-1) + 2*
     +    ms**2*s*u1**(-1) - 2*mjs2*s*tj**(-1) + 2*mjs2*s*u1**(-1)
     +     - 2*mjs2*tj**(-1)*u1 - 6*mjs2*tj*u1**(-1) + 2*mjs2**2*s*
     +  tj**(-1)*u1**(-1) + 4*mjs2**2*tj**(-1) - 8*mjs2**2*u1**(-1) - 4
     +    *mjs2**3*tj**(-1)*u1**(-1) - 2*tj**2*u1**(-1) )
      gqH = gqH + NSANG(2,6,-1,-2)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +    8*ms**2*mjs2 - 8*ms**2*mjs2**2*u1**(-1) - 4*ms**2*u1 )
      gqH = gqH + NSANG(2,6,-1,-1)*s4**(-1)*CA**(-2)*hardfac**(-1)*
     + coeff * ( ms**2*mjs2*s**(-1)*u1 + 4*ms**2*mjs2*tj*u1**(-1)
     +     + 4*ms**2*mjs2**2*u1**(-1) - ms**2*s**(-1)*u1**2 - ms**2*u1
     +     - 4*mjs2*s**(-1)*u1**2 + 2*mjs2*tj - 4*mjs2*u1 + 5*
     +    mjs2**2*s**(-1)*u1 + 4*mjs2**2 - 2*mjs2**3*s**(-1) + s**(-1)*
     +    u1**3 - tj*u1 + u1**2 )
      gqH = gqH + NSANG(2,6,-1,-1)*CA**(-2)*sptj**(-1)*hardfac**(-1)*
     + coeff * (  - 4*ms**2*mjs2*tj*u1**(-1) - 2*mjs2*tj + 2*mjs2
     +    *u1 - 2*mjs2**2 + tj*u1 - u1**2 )
      gqH = gqH + NSANG(2,6,-1,-1)*CA**(-2)*hardfac**(-1)*coeff * (  - 
     +    ms**2*mjs2*s**(-1) + ms**2*s**(-1)*u1 + mjs2*s**(-1)*u1
     +     - mjs2**2*s**(-1) )
      gqH = gqH + NSANG(2,6,-1,-1)*CA**(-1)*CF*sptj**(-1)*hardfac**(-1)
     + *coeff * (  - 8*ms**2*mjs2*tj*u1**(-1) - 4*mjs2*tj + 4*
     +    mjs2*u1 - 4*mjs2**2 + 2*tj*u1 - 2*u1**2 )
      gqH = gqH + NSANG(2,6,-1,-1)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +    2*ms**2*mjs2*s*tj**(-1)*u1**(-1) + 10*ms**2*mjs2*u1**(-1)
     +     - 4*ms**2*s*tj**(-1) - 2*ms**2 - 2*mjs2*s*tj**(-1) - 2*
     +    mjs2*tj**(-1)*u1 + 2*mjs2 + 2*mjs2**2*s*tj**(-1)*u1**(-1)
     +   + 4*mjs2**2*tj**(-1) - 2*mjs2**2*u1**(-1) - 4*mjs2**3*tj**(-1)*
     +    u1**(-1) - 2*u1 )
      gqH = gqH + NSANG(5,0,-2,0)*CA**(-2)*hardfac**(-1)*coeff * ( 2*
     +    mjs2 )
      gqH = gqH + NSANG(5,0,-2,0)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 2
     +    *mjs2*s**(-1)*u1 + 6*mjs2 - 2*mjs2**2*s**(-1) )
      gqH = gqH + NSANG(5,0,-1,0)*s4**(-1)*CA**(-2)*hardfac**(-1)*coeff
     +  * ( ms**2*mjs2*u1**(-1) + 3*ms**2*tj*u1**(-1) + ms**2 + 2*
     +    mjs2*tj*u1**(-1) + mjs2 + mjs2**2*u1**(-1) + tj + tj**2*
     +    u1**(-1) )
      gqH = gqH + NSANG(5,0,-1,0)*s4**(-1)*CA**(-1)*CF*hardfac**(-1)*
     + coeff * ( 6*ms**2*mjs2*s**(-1) - 2*ms**2*mjs2*u1**(-1) - 2*
     +    ms**2*mjs2**2*s**(-1)*u1**(-1) - 4*ms**2*s**(-1)*u1 + 4*ms**2*
     +    tj*u1**(-1) - 4*ms**2 - 2*mjs2*tj*u1**(-1) - 4*mjs2**2*
     +    u1**(-1) + 2*mjs2**3*s**(-1)*u1**(-1) + 2*tj )
      gqH = gqH + NSANG(5,0,-1,0)*CA**(-2)*hardfac**(-1)*coeff * ( 
     +    ms**2*mjs2*s**(-1)*u1**(-1) - ms**2*s**(-1)*tj*u1**(-1) - 
     +    ms**2*s**(-1) - 8*ms**2*u1**(-1) - mjs2*s**(-1)*tj*u1**(-1)
     +     - mjs2*s**(-1) - 6*mjs2*u1**(-1) + mjs2**2*s**(-1)*
     +    u1**(-1) - 2*tj*u1**(-1) )
      gqH = gqH + NSANG(5,0,-1,0)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +     - 2 + 2*ms**2*mjs2*s**(-1)*u1**(-1) - 2*ms**2*s**(-1) - 8*
     +    ms**2*u1**(-1) + 6*mjs2*s**(-1)*tj*u1**(-1) + 4*mjs2*
     +    s**(-1) + 2*mjs2*u1**(-1) - 6*mjs2**2*s**(-1)*u1**(-1) - 4*
     +    s**(-1)*tj - 2*s**(-1)*tj**2*u1**(-1) - 2*s**(-1)*u1 - 6*tj*
     +    u1**(-1) )
      gqH = gqH + NSANG(5,6,-1,-1)*CA**(-2)*sptj**(-1)*hardfac**(-1)*
     + coeff * ( 4*ms**2*mjs2*tj*u1**(-1) + 2*mjs2*tj + 2*mjs2*u1
     +     - 2*mjs2**2 - tj*u1 - u1**2 )
      gqH = gqH + NSANG(5,6,-1,-1)*CA**(-2)*hardfac**(-1)*coeff * ( 
     +    ms**2*mjs2*s**(-1)*tj*u1**(-1) + ms**2*mjs2*u1**(-1) - 
     +    ms**2*s**(-1)*tj - ms**2*s**(-1)*tj**2*u1**(-1) - 3*ms**2*tj*
     +    u1**(-1) - 2*ms**2 - 7*mjs2*s**(-1)*tj - 4*mjs2*s**(-1)*
     +    tj**2*u1**(-1) - 3*mjs2*s**(-1)*u1 - 2*mjs2*tj*u1**(-1)
     +     - mjs2 + 5*mjs2**2*s**(-1)*tj*u1**(-1) + 4*mjs2**2*s**(-1)
     +    + mjs2**2*u1**(-1) - 2*mjs2**3*s**(-1)*u1**(-1) + 3*s**(-1)*tj
     +    *u1 + 3*s**(-1)*tj**2 + s**(-1)*tj**3*u1**(-1) + s**(-1)*
     +    u1**2 + tj + tj**2*u1**(-1) )
      gqH = gqH + NSANG(5,6,-1,-1)*CA**(-1)*CF*sptj**(-1)*hardfac**(-1)
     + *coeff * ( 8*ms**2*mjs2*tj*u1**(-1) + 4*mjs2*tj + 4*mjs2*
     +    u1 - 4*mjs2**2 - 2*tj*u1 - 2*u1**2 )
      gqH = gqH + NSANG(5,6,-1,-1)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +     - 2*ms**2*mjs2*s*tj**(-1)*u1**(-1) - 10*ms**2*mjs2*
     +    u1**(-1) + 4*ms**2*s*tj**(-1) + 6*ms**2*s*u1**(-1) + 2*ms**2*
     +    tj*u1**(-1) + 2*ms**2 + 2*mjs2*s*tj**(-1) + 4*mjs2*s*
     +    u1**(-1) + 2*mjs2*tj**(-1)*u1 + 8*mjs2*tj*u1**(-1) + 10*
     +    mjs2 - 2*mjs2**2*s*tj**(-1)*u1**(-1) - 4*mjs2**2*tj**(-1) - 
     +    10*mjs2**2*u1**(-1) + 4*mjs2**3*tj**(-1)*u1**(-1) - 2*s*tj*
     +    u1**(-1) - 2*s - 6*tj - 2*tj**2*u1**(-1) - 4*u1 )
      gqH = gqH + NSANG(5,7,-2,1)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 2
     +    *mjs2*s**(-1) )
      gqH = gqH + NSANG(5,7,-1,1)*s4**(-1)*CA**(-1)*CF*hardfac**(-1)*
     + coeff * ( 2*ms**2*mjs2*s**(-1)*u1**(-1) - 4*ms**2*s**(-1) - 6*
     +    ms**2*u1**(-1) - 2*mjs2*s**(-1) - 4*mjs2*u1**(-1) + 2*
     +    mjs2**2*s**(-1)*u1**(-1) - 2*tj*u1**(-1) )
      gqH = gqH + NSANG(5,7,-1,1)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +     - 2*ms**2*s**(-1)*u1**(-1) - 2*mjs2*s**(-1)*u1**(-1) )
      gqH = gqH + NSANG(6,0,-2,0)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 8
     +    *ms**2*mjs2*u1**(-1) - 4*ms**2 + 8*ms**4*mjs2*u1**(-2) )
      gqH = gqH + NSANG(6,0,-1,0)*s4**(-1)*CA**(-2)*hardfac**(-1)*coeff
     +  * (  - 4*ms**2*mjs2*u1**(-1) - ms**2*s**(-1)*u1 + 2*ms**2 - 8
     +    *ms**4*mjs2*u1**(-2) - 3*mjs2*s**(-1)*u1 - 2*mjs2 + 2*
     +    mjs2**2*s**(-1) + s**(-1)*u1**2 + u1 )
      gqH = gqH + NSANG(6,0,-1,0)*CA**(-2)*sptj**(-1)*hardfac**(-1)*
     + coeff * ( 4*ms**2*mjs2*u1**(-1) + 4*mjs2 - 2*u1 )
      gqH = gqH + NSANG(6,0,-1,0)*CA**(-2)*hardfac**(-1)*coeff * (  - 1
     +     - ms**2*s**(-1)*tj*u1**(-1) + ms**2*s**(-1) - 3*ms**2*
     +    u1**(-1) - 3*mjs2*s**(-1)*tj*u1**(-1) - mjs2*s**(-1) + 
     +    mjs2*u1**(-1) + 2*mjs2**2*s**(-1)*u1**(-1) + 2*s**(-1)*tj + 
     +    s**(-1)*tj**2*u1**(-1) + s**(-1)*u1 - s*u1**(-1) )
      gqH = gqH + NSANG(6,0,-1,0)*CA**(-1)*CF*sptj**(-1)*hardfac**(-1)*
     + coeff * ( 8*ms**2*mjs2*u1**(-1) + 8*mjs2 - 4*u1 )
      gqH = gqH + NSANG(6,0,-1,0)*CA**(-1)*CF*hardfac**(-1)*coeff * ( 
     +     - 4 + 8*ms**2*mjs2*u1**(-2) + 8*mjs2*u1**(-1) - 2*s*
     +    u1**(-1) - 2*tj*u1**(-1) )
      gqH = gqH + NSANG(6,7,-1,1)*s4**(-1)*CA**(-2)*hardfac**(-1)*coeff
     +  * ( 4*ms**2*mjs2*s**(-1)*u1**(-1) + 2*mjs2*s**(-1) - 
     +    s**(-1)*u1 )
      gqH = gqH + NSANG(6,7,-1,1)*CA**(-2)*hardfac**(-1)*coeff * ( 2*
     +    mjs2*s**(-1)*u1**(-1) - s**(-1)*tj*u1**(-1) - s**(-1) - 
     +    u1**(-1) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4**(-1)*pi*u1hat**(-2)*
     + CA**(-1)*CF*coeff * ( 16*ms**2*mjs2 )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4**(-1)*pi*u1hat**(-1)*
     + CA**(-1)*CF*coeff * ( 16*mjs2*s**(-1)*tj*t1**(-1)*uj + 16*
     +    mjs2*tj*t1**(-1) + 16*mjs2**2*s**(-1)*t1**(-1)*uj + 16*
     +    mjs2**2*t1**(-1) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4**(-1)*pi*CA**(-1)*CF*
     + coeff * (  - 8*s**(-1)*tj*t1**(-1)*uj - 8*tj*t1**(-1) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4**(-1)*pi*coeff * ( 16*
     +    ms**2*mjs2*u1**(-2) - 16*mjs2*s**(-1)*tj*u1**(-1) + 16*
     +    mjs2**2*s**(-1)*tj*u1**(-2) + 16*mjs2**2*u1**(-2) + 8*s**(-1)*
     +    tj )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4*pi*sptj**(-2)*coeff * ( 
     +    16*ms**2*mjs2*u1**(-2) - 16*mjs2*s**(-1)*tj*u1**(-1) + 8*
     +    s**(-1)*tj )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4*pi*sptj**(-1)*coeff * ( 
     +    16*mjs2**2*s**(-1)*u1**(-2) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*s4*pi*coeff * ( 16*ms**2*
     +    mjs2*u1**(-4) - 16*mjs2*s**(-1)*tj*u1**(-3) + 16*mjs2**2*
     +    s**(-1)*tj*u1**(-4) + 16*mjs2**2*u1**(-4) + 8*s**(-1)*tj*
     +    u1**(-2) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*pi*u1hat**(-2)*CA**(-1)*CF*
     + spuj**(-1)*coeff * (  - 8*ms**2*mjs2 )
      gqH = gqH + log(1 + ms**2*s4**(-1))*pi*u1hat**(-2)*CA**(-1)*CF*
     + coeff * (  - 8*ms**2*mjs2*t1**(-1) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*pi*u1hat**(-1)*CA**(-1)*CF*
     + coeff * (  - 8*mjs2*s**(-1)*tj*t1**(-2)*uj - 8*mjs2*s**(-1)*
     +    tj*t1**(-1) - 8*mjs2*tj*t1**(-2) - 8*mjs2**2*s**(-1)*
     +   t1**(-2)*uj - 8*mjs2**2*s**(-1)*t1**(-1) - 8*mjs2**2*t1**(-2) )
      gqH = gqH + log(1 + ms**2*s4**(-1))*pi*CA**(-1)*CF*coeff * ( 4*
     +    s**(-1)*tj*t1**(-2)*uj + 4*s**(-1)*tj*t1**(-1) + 4*tj*
     +    t1**(-2) )
      gqH = gqH + log(ms**(-2)*s4)*s4**(-1)*pi*u1hat**(-2)*CA**(-1)*CF
     + *coeff * (  - 16*ms**2*mjs2 )
      gqH = gqH + log(ms**(-2)*s4)*s4**(-1)*pi*u1hat**(-1)*CA**(-1)*CF
     + *coeff * (  - 16*mjs2*s**(-1)*tj*t1**(-1)*uj - 16*mjs2*tj*
     +   t1**(-1) - 16*mjs2**2*s**(-1)*t1**(-1)*uj - 16*mjs2**2*t1**(-1)
     +     )
      gqH = gqH + log(ms**(-2)*s4)*s4**(-1)*pi*CA**(-1)*CF*coeff * ( 8
     +    *s**(-1)*tj*t1**(-1)*uj + 8*tj*t1**(-1) )
      gqH = gqH + log(ms**(-2)*s4)*s4**(-1)*pi*coeff * (  - 16*ms**2*
     +    mjs2*u1**(-2) + 16*mjs2*s**(-1)*tj*u1**(-1) - 16*mjs2**2*
     +    s**(-1)*tj*u1**(-2) - 16*mjs2**2*u1**(-2) - 8*s**(-1)*tj )
      gqH = gqH + log(ms**(-2)*s4)*s4*pi*sptj**(-2)*coeff * (  - 16*
     +    ms**2*mjs2*u1**(-2) + 16*mjs2*s**(-1)*tj*u1**(-1) - 8*
     +    s**(-1)*tj )
      gqH = gqH + log(ms**(-2)*s4)*s4*pi*sptj**(-1)*coeff * (  - 16*
     +    mjs2**2*s**(-1)*u1**(-2) )
      gqH = gqH + log(ms**(-2)*s4)*s4*pi*coeff * (  - 16*ms**2*mjs2*
     +    u1**(-4) + 16*mjs2*s**(-1)*tj*u1**(-3) - 16*mjs2**2*s**(-1)*
     +    tj*u1**(-4) - 16*mjs2**2*u1**(-4) - 8*s**(-1)*tj*u1**(-2) )
      gqH = gqH + log(ms**(-2)*s4)*pi*u1hat**(-2)*CA**(-1)*CF*
     + spuj**(-1)*coeff * ( 8*ms**2*mjs2 )
      gqH = gqH + log(ms**(-2)*s4)*pi*u1hat**(-2)*CA**(-1)*CF*coeff
     +  * ( 8*ms**2*mjs2*t1**(-1) )
      gqH = gqH + log(ms**(-2)*s4)*pi*u1hat**(-1)*CA**(-1)*CF*coeff
     +  * ( 8*mjs2*s**(-1)*tj*t1**(-2)*uj + 8*mjs2*s**(-1)*tj*
     +    t1**(-1) + 8*mjs2*tj*t1**(-2) + 8*mjs2**2*s**(-1)*t1**(-2)*
     +    uj + 8*mjs2**2*s**(-1)*t1**(-1) + 8*mjs2**2*t1**(-2) )
      gqH = gqH + log(ms**(-2)*s4)*pi*CA**(-1)*CF*coeff * (  - 4*
     +    s**(-1)*tj*t1**(-2)*uj - 4*s**(-1)*tj*t1**(-1) - 4*tj*
     +    t1**(-2) )
      gqH = gqH + log(ms**(-2)*QF**2)*s4**(-1)*pi*u1hat**(-2)*CA**(-1)
     + *CF*coeff * ( 16*ms**2*mjs2 )
      gqH = gqH + log(ms**(-2)*QF**2)*s4**(-1)*pi*u1hat**(-1)*CA**(-1)
     + *CF*coeff * ( 16*mjs2*s**(-1)*tj*t1**(-1)*uj + 16*mjs2*tj*
     +   t1**(-1) + 16*mjs2**2*s**(-1)*t1**(-1)*uj + 16*mjs2**2*t1**(-1)
     +     )
      gqH = gqH + log(ms**(-2)*QF**2)*s4**(-1)*pi*CA**(-1)*CF*coeff
     +  * (  - 8*s**(-1)*tj*t1**(-1)*uj - 8*tj*t1**(-1) )
      gqH = gqH + log(ms**(-2)*QF**2)*s4**(-1)*pi*coeff * ( 16*ms**2*
     +    mjs2*u1**(-2) - 16*mjs2*s**(-1)*tj*u1**(-1) + 16*mjs2**2*
     +    s**(-1)*tj*u1**(-2) + 16*mjs2**2*u1**(-2) + 8*s**(-1)*tj )
      gqH = gqH + log(ms**(-2)*QF**2)*s4*pi*sptj**(-2)*coeff * ( 16*
     +    ms**2*mjs2*u1**(-2) - 16*mjs2*s**(-1)*tj*u1**(-1) + 8*
     +    s**(-1)*tj )
      gqH = gqH + log(ms**(-2)*QF**2)*s4*pi*sptj**(-1)*coeff * ( 16*
     +    mjs2**2*s**(-1)*u1**(-2) )
      gqH = gqH + log(ms**(-2)*QF**2)*s4*pi*coeff * ( 16*ms**2*mjs2*
     +    u1**(-4) - 16*mjs2*s**(-1)*tj*u1**(-3) + 16*mjs2**2*s**(-1)*
     +    tj*u1**(-4) + 16*mjs2**2*u1**(-4) + 8*s**(-1)*tj*u1**(-2) )
      gqH = gqH + log(ms**(-2)*QF**2)*pi*u1hat**(-2)*CA**(-1)*CF*
     + spuj**(-1)*coeff * (  - 8*ms**2*mjs2 )
      gqH = gqH + log(ms**(-2)*QF**2)*pi*u1hat**(-2)*CA**(-1)*CF*coeff
     +  * (  - 8*ms**2*mjs2*t1**(-1) )
      gqH = gqH + log(ms**(-2)*QF**2)*pi*u1hat**(-1)*CA**(-1)*CF*coeff
     +  * (  - 8*mjs2*s**(-1)*tj*t1**(-2)*uj - 8*mjs2*s**(-1)*tj*
     +    t1**(-1) - 8*mjs2*tj*t1**(-2) - 8*mjs2**2*s**(-1)*t1**(-2)*
     +    uj - 8*mjs2**2*s**(-1)*t1**(-1) - 8*mjs2**2*t1**(-2) )
      gqH = gqH + log(ms**(-2)*QF**2)*pi*CA**(-1)*CF*coeff * ( 4*
     +    s**(-1)*tj*t1**(-2)*uj + 4*s**(-1)*tj*t1**(-1) + 4*tj*
     +    t1**(-2) )


c               the prefactor for the scaling functions 
      NS_QGH = gqH * (abs(m1)+abs(m2))**2/4.D0

      end
