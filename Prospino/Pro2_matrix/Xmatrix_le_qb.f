cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     LE_QB(MASSIN,LUMI)                                                c
c     LE_QBOS(MASSIN,LUMI)                                              c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = t2                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(25) = gamma_s                                           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function LE_QBOS(massin,lumi)

      implicit none 

      real*8     massin(1:30),lumi(1:3),pi,Cf,Nc,alphas
     &          ,Ar2,st,MM_s3s
     &          ,s,s4,m1,theta_s3,t2,u2
     &          ,ANG_fin(0:9,0:9,-2:2,-2:2)

      pi = 4.D0*atan(1.D0)
      Cf = 4.D0/3.D0
      Nc = 3.D0

      s    = massin(1)
      t2   = massin(2)
      s4   = massin(3)
      m1   = massin(6)

c               real kinematics built in
      u2  = s4 - s - t2 + m1**2

c               the s3 regularization
      theta_s3 = 0.D0
      if (s.gt.4.D0*m1**2) theta_s3 = 1.D0 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      Ar2    = 1.D0
      st     = 1.D0
      alphas = 1.D0/4.D0/pi

c               the angular functions 
      call ANGULAR_ARRAY_LE_QBOS(massin,theta_s3,ANG_fin)

c               the form output
c               no replacements
      MM_s3s =
     +  + lumi(3)*ANG_fin(7,8,1,-2)*Ar2*Nc*Cf*st**2*Pi**2*alphas**2
     +  * ( 256*s**(-2)*m1**2*t2 + 256*s**(-2)*m1**2*u2 - 256*s**(-2)*
     +    m1**2*s4 - 256*s**(-2)*m1**4 + 128*s**(-1)*m1**2 )
      MM_s3s = MM_s3s + lumi(3)*ANG_fin(7,8,2,-2)*Ar2*Nc*Cf*st**2*Pi**2
     + *alphas**2 * (  - 128*s**(-2)*m1**2 )
      MM_s3s = MM_s3s + lumi(3)*ANG_fin(8,0,-2,0)*Ar2*Nc*Cf*st**2*Pi**2
     + *alphas**2 * (  - 256*s**(-2)*m1**2*t2*u2 + 256*s**(-2)*m1**2*t2
     +    *s4 - 128*s**(-2)*m1**2*t2**2 + 256*s**(-2)*m1**2*u2*s4 - 128
     +    *s**(-2)*m1**2*u2**2 - 128*s**(-2)*m1**2*s4**2 + 256*s**(-2)*
     +    m1**4*t2 + 256*s**(-2)*m1**4*u2 - 256*s**(-2)*m1**4*s4 - 128*
     +    s**(-2)*m1**6 - 128*s**(-1)*m1**2*t2 - 128*s**(-1)*m1**2*u2
     +     + 128*s**(-1)*m1**2*s4 )

c               the prefactor for the scaling functions 
      LE_QBOS = MM_s3s * m1**2/4.D0

c               the phase space except for 1/s**2 
      LE_QBOS = LE_QBOS / ( 16.D0 * pi**2 )**2 / 2.D0 * s4/(s4+m1**2)

c               the averaging factors
      LE_QBOS = LE_QBOS /4.D0 /Nc**2

      return
      end

c-----------------------------------------------------------------------
      real*8 function LE_QB(massin,lumi)

      implicit none 

      real*8     massin(1:30),lumi(1:3),pi,Cf,Nc,alphas
     &          ,Ar2,st,log_qf,MM_qb,log_all
     &          ,s,t2,s4,m1,m12,u2,u2s,s4j
     &          ,hardfac
     &          ,ANG_fin(0:9,0:9,-2:2,-2:2)

      pi = 4.D0*atan(1.D0)
      Cf = 4.D0/3.D0
      Nc = 3.D0

      s    = massin(1)
      t2   = massin(2)
      s4   = massin(3)
      m1   = massin(6)
      m12  = m1**2 

c               the factorization scale 
      log_qf = log( massin(13)**2/m1**2 ) 
      log_all = log_qf - log(s4/m1**2) + log(1.D0+m1**2/s4)

c               real kinematics built in
      u2  = s4 - s - t2 + m1**2 
      u2s = u2 - m1**2 

      s4j = s4 + m1**2

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      Ar2    = 1.D0
      st     = 1.D0
      alphas = 1.D0/4.D0/pi
      hardfac = 1.D0

c               the angular functions 
      call ANGULAR_ARRAY_LE_QB(massin,ANG_fin)

c               the form output
c               replacements: [] -> ()
      MM_qb =
     +  + lumi(1)*Ar2*Nc*Cf*st**2*Pi**3*alphas**2*hardfac*log_all * ( 
     +    256*s**(-1)*m1**2*t2*s4**(-1)*u2s**(-1) + 512*s**(-1)*m1**2*
     +    t2*u2s**(-2) + 256*s**(-1)*m1**2*t2**2*s4**(-1)*u2s**(-2) + 
     +    512*s**(-1)*m1**2*t2**2*u2s**(-3) + 128*s**(-1)*m1**2*
     +    s4**(-1) + 256*s**(-1)*m1**2*u2s**(-1) + 512*s**(-1)*m1**4*t2
     +    *s4**(-1)*u2s**(-2) + 512*s**(-1)*m1**4*t2*u2s**(-3) + 512*
     +    s**(-1)*m1**4*t2**2*s4**(-1)*u2s**(-3) + 512*s**(-1)*m1**4*
     +    t2**2*u2s**(-4) + 256*s**(-1)*m1**4*s4**(-1)*u2s**(-1) + 256*
     +    s**(-1)*m1**4*u2s**(-2) + 512*s**(-1)*m1**6*t2*s4**(-1)*
     +    u2s**(-3) + 512*s**(-1)*m1**6*t2**2*s4**(-1)*u2s**(-4) + 256*
     +    s**(-1)*m1**6*s4**(-1)*u2s**(-2) + 256*s**(-1)*t2*u2s**(-1)
     +     + 256*s**(-1)*t2**2*u2s**(-2) + 128*s**(-1) - 128*
     +    (s+t2)**(-1)*m1**2*s4**(-1) - 256*(s+t2)**(-1)*m1**2*
     +    u2s**(-1) - 256*(s+t2)**(-1)*m1**4*s4**(-1)*u2s**(-1) - 256*
     +    (s+t2)**(-1)*m1**4*u2s**(-2) - 256*(s+t2)**(-1)*m1**6*
     +    s4**(-1)*u2s**(-2) )
      MM_qb = MM_qb + lumi(1)*Ar2*Nc*Cf*st**2*Pi**3*alphas**2*hardfac*
     + log_all * (  - 128*(s+t2)**(-1) + 256*m1**2*t2*s4**(-1)*
     +    u2s**(-2) + 512*m1**2*t2*u2s**(-3) + 512*m1**4*t2*s4**(-1)*
     +    u2s**(-3) + 512*m1**4*t2*u2s**(-4) + 512*m1**6*t2*s4**(-1)*
     +    u2s**(-4) + 256*t2*u2s**(-2) )
      MM_qb = MM_qb + lumi(1)*Ar2*Nc*Cf*st**2*Pi**3*alphas**2*hardfac
     +  * (  - 512*s**(-1)*m1**2*t2*u2s**(-2) - 512*s**(-1)*m1**2*t2**2
     +    *u2s**(-3) - 128*s**(-1)*m1**2*s4**(-1) - 256*s**(-1)*m1**2*
     +    u2s**(-1) - 512*s**(-1)*m1**4*t2*s4**(-1)*u2s**(-2) - 512*
     +    s**(-1)*m1**4*t2*u2s**(-3) - 512*s**(-1)*m1**4*t2**2*s4**(-1)
     +    *u2s**(-3) - 512*s**(-1)*m1**4*t2**2*u2s**(-4) - 256*s**(-1)*
     +    m1**4*s4**(-1)*u2s**(-1) - 256*s**(-1)*m1**4*u2s**(-2) - 512*
     +    s**(-1)*m1**6*t2*s4**(-1)*u2s**(-3) - 512*s**(-1)*m1**6*t2**2
     +    *s4**(-1)*u2s**(-4) - 256*s**(-1)*m1**6*s4**(-1)*u2s**(-2) - 
     +    128*s**(-1) + 128*(s+t2)**(-1)*m1**2*s4**(-1) + 256*
     +    (s+t2)**(-1)*m1**2*u2s**(-1) + 256*(s+t2)**(-1)*m1**4*
     +    s4**(-1)*u2s**(-1) + 256*(s+t2)**(-1)*m1**4*u2s**(-2) + 256*
     +    (s+t2)**(-1)*m1**6*s4**(-1)*u2s**(-2) + 128*(s+t2)**(-1) - 
     +    512*m1**2*t2*u2s**(-3) - 512*m1**4*t2*s4**(-1)*u2s**(-3) - 
     +    512*m1**4*t2*u2s**(-4) - 512*m1**6*t2*s4**(-1)*u2s**(-4) )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,0,-2,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 256*s*m1**2*t2*u2s**(-2) - 128*s**2*
     +    m1**2*u2s**(-2) - 128*m1**2*t2**2*u2s**(-2) )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,0,-1,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64 - 128*s*m1**2*u2s**(-2) - 64*s*
     +    u2s**(-1) - 128*m1**2*t2*u2s**(-2) - 128*m1**2*u2s**(-1) - 
     +    128*m1**4*u2s**(-2) - 64*t2*u2s**(-1) )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,5,-2,-2)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**2*m1**2 )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,5,-2,-1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 256*s*m1**2*t2*u2s**(-1) - 256*s**2*
     +    m1**2*u2s**(-1) )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,5,-1,-2)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s*m1**2 )
      MM_qb = MM_qb + lumi(1)*ANG_fin(1,5,-1,-1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 256*s*m1**2*u2s**(-1) - 128*s*t2*
     +    u2s**(-1) - 64*s - 128*m1**2*t2*u2s**(-1) - 64*m1**2 - 128*
     +    m1**4*u2s**(-1) - 128*t2 - 128*t2**2*u2s**(-1) - 64*u2 )
      MM_qb = MM_qb + lumi(1)*ANG_fin(5,0,-2,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64*m1**2 )
      MM_qb = MM_qb + lumi(1)*ANG_fin(5,0,-1,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64*m1**2*u2s**(-1) - 64*t2*u2s**(-1) )
      MM_qb = MM_qb + lumi(2)*Ar2*Cf*st**2*Pi**3*alphas**2*hardfac * ( 
     +     - 128*s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**2*t2*u2
     +     - 128*s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**2*t2**2
     +    *u2*s4**(-1) - 128*s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)
     +    *m1**2*t2**2 - 128*s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)
     +    *m1**2*t2**3*s4**(-1) - 128*s**(-1)*(t2+u2)**(-1)*
     +    (t2*s4j+s*m12)**(-1)*m1**4*t2*u2*s4**(-1) - 128*s**(-1)*
     +    (t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**4*t2**2*s4**(-1) - 128
     +    *s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*t2**2*u2 - 128*
     +    s**(-1)*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*t2**3 - 256*
     +    s**(-1)*(t2+u2)**(-1)*m1**2*t2*u2s**(-1) - 128*s**(-1)*
     +    (t2+u2)**(-1)*m1**2*t2**2*s4**(-1)*u2s**(-1) - 256*s**(-1)*
     +    (t2+u2)**(-1)*m1**4*t2*s4**(-1)*u2s**(-1) - 128*s**(-1)*
     +    (t2+u2)**(-1)*m1**4*u2s**(-1) - 128*s**(-1)*(t2+u2)**(-1)*
     +    m1**6*s4**(-1)*u2s**(-1) - 128*s**(-1)*(t2+u2)**(-1)*t2**2*
     +    u2s**(-1) )
      MM_qb = MM_qb + lumi(2)*Ar2*Cf*st**2*Pi**3*alphas**2*hardfac * ( 
     +    128*s**(-1)*m1**2*t2*s4**(-1)*u2s**(-1) + 128*s**(-1)*m1**2*
     +    u2s**(-1) + 128*s**(-1)*m1**4*s4**(-1)*u2s**(-1) + 128*
     +    s**(-1)*t2*u2s**(-1) - 256*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)
     +    *m1**2*t2 - 128*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**2*
     +    t2**2*s4**(-1) - 256*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**4
     +    *t2*s4**(-1) - 128*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**4
     +     - 128*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*m1**6*s4**(-1) - 
     +    128*(t2+u2)**(-1)*(t2*s4j+s*m12)**(-1)*t2**2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(1,0,-1,0)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64*s**(-1)*(t2+u2)**(-1)*t2*u2 - 64*
     +    s**(-1)*(t2+u2)**(-1)*t2**2 + 64*s**(-1)*t2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(1,5,-1,-2)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s*m1**2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(1,5,-1,-1)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64*s*(t2+u2)**(-1)*m1**2 - 64*s*
     +    (t2+u2)**(-1)*t2 + 128*s*m1**2*u2s**(-1) - 64*(t2+u2)**(-1)*
     +    m1**2*t2 - 64*(t2+u2)**(-1)*m1**2*u2 - 64*(t2+u2)**(-1)*m1**4
     +     + 64*(t2+u2)**(-1)*t2*u2 + 128*m1**2*t2*u2s**(-1) + 64*m1**2
     +     - 64*t2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(1,8,-1,-1)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 64*s**(-1)*(t2+u2)**(-1)*m1**2*t2*u2 - 
     +    64*s**(-1)*(t2+u2)**(-1)*m1**2*t2**2 - 64*s**(-1)*
     +    (t2+u2)**(-1)*t2**2*u2 - 64*s**(-1)*(t2+u2)**(-1)*t2**3 + 64*
     +    s**(-1)*m1**2*t2 + 64*s**(-1)*m1**4*t2*u2s**(-1) + 64*s**(-1)
     +    *t2**2 + 64*s**(-1)*t2**3*u2s**(-1) - 64*s*(t2+u2)**(-1)*
     +    m1**2 - 64*s*(t2+u2)**(-1)*t2 + 64*s*m1**2*u2s**(-1) + 64*s*
     +    t2*u2s**(-1) - 128*(t2+u2)**(-1)*m1**2*t2 - 128*(t2+u2)**(-1)
     +    *m1**2*u2 - 64*(t2+u2)**(-1)*m1**4 - 64*(t2+u2)**(-1)*t2*u2
     +     - 128*(t2+u2)**(-1)*t2**2 + 64*m1**2*t2*u2s**(-1) + 128*
     +    m1**2 + 64*m1**4*u2s**(-1) + 64*t2 + 128*t2**2*u2s**(-1) )
      MM_qb = MM_qb + lumi(2)*ANG_fin(5,0,-2,0)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s**(-1)*m1**2*u2 + 256*m1**2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(5,0,-1,0)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 64 - 64*s**(-1)*(t2+u2)**(-1)*m1**2*t2 - 
     +    64*s**(-1)*(t2+u2)**(-1)*m1**2*u2 + 64*s**(-1)*m1**2*t2*
     +    u2s**(-1) + 128*s**(-1)*m1**2 + 128*s**(-1)*m1**4*u2s**(-1)
     +     + 64*s**(-1)*t2 + 64*s**(-1)*t2**2*u2s**(-1) - 192*
     +    (t2+u2)**(-1)*m1**2 - 128*(t2+u2)**(-1)*t2 - 64*(t2+u2)**(-1)
     +    *u2 + 192*m1**2*u2s**(-1) + 64*t2*u2s**(-1) )
      MM_qb = MM_qb + lumi(2)*ANG_fin(5,7,-2,1)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s**(-1)*m1**2 )
      MM_qb = MM_qb + lumi(2)*ANG_fin(5,7,-1,1)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**(-1)*(t2+u2)**(-1)*m1**2 - 64*
     +    s**(-1)*(t2+u2)**(-1)*t2 - 64*s**(-1)*(t2+u2)**(-1)*u2 + 64*
     +    s**(-1)*m1**2*u2s**(-1) + 64*s**(-1)*t2*u2s**(-1) + 64*
     +    s**(-1) )
      MM_qb = MM_qb + lumi(2)*ANG_fin(7,8,1,-1)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**(-1)*(t2+u2)**(-1)*m1**2 - 64*
     +    s**(-1)*(t2+u2)**(-1)*t2 - 64*s**(-1)*(t2+u2)**(-1)*u2 + 64*
     +    s**(-1)*m1**2*u2s**(-1) + 64*s**(-1)*t2*u2s**(-1) + 64*
     +    s**(-1) )
      MM_qb = MM_qb + lumi(2)*ANG_fin(8,0,-1,0)*Ar2*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 64 - 64*s**(-1)*(t2+u2)**(-1)*m1**2*t2 - 
     +    64*s**(-1)*(t2+u2)**(-1)*m1**2*u2 - 64*s**(-1)*(t2+u2)**(-1)*
     +    t2*u2 - 64*s**(-1)*(t2+u2)**(-1)*t2**2 + 64*s**(-1)*m1**2 + 
     +    64*s**(-1)*m1**4*u2s**(-1) + 64*s**(-1)*t2 + 64*s**(-1)*t2**2
     +    *u2s**(-1) - 192*(t2+u2)**(-1)*m1**2 - 128*(t2+u2)**(-1)*t2
     +     - 64*(t2+u2)**(-1)*u2 + 128*m1**2*u2s**(-1) + 128*t2*
     +    u2s**(-1) )
      MM_qb = MM_qb + lumi(3)*ANG_fin(0,0,0,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 64*s**(-2)*(t2+u2)**(-1)*t2*u2 + 64*
     +    s**(-2)*(t2+u2)**(-1)*u2**2 - 64*s**(-2)*u2 )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,0,-2,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**(-2)*m1**2*u2**2 - 128*s**(-1)*
     +    m1**2*u2 - 64*m1**2 )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,0,-1,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s**(-2)*(t2+u2)**(-1)*m1**2*t2*u2 + 
     +    128*s**(-2)*(t2+u2)**(-1)*m1**2*u2**2 + 128*s**(-2)*
     +    (t2+u2)**(-1)*t2*u2**2 + 64*s**(-2)*(t2+u2)**(-1)*t2**2*u2 + 
     +    64*s**(-2)*(t2+u2)**(-1)*u2**3 - 128*s**(-2)*m1**2*u2 - 64*
     +    s**(-2)*t2*u2 - 64*s**(-2)*u2**2 + 192*s**(-1)*(t2+u2)**(-1)*
     +    m1**2*t2 + 320*s**(-1)*(t2+u2)**(-1)*m1**2*u2 + 128*s**(-1)*
     +    (t2+u2)**(-1)*m1**4 + 128*s**(-1)*(t2+u2)**(-1)*t2*u2 + 64*
     +    s**(-1)*(t2+u2)**(-1)*t2**2 + 64*s**(-1)*(t2+u2)**(-1)*u2**2
     +     - 320*s**(-1)*m1**2 - 64*s**(-1)*u2 + 64*(t2+u2)**(-1)*m1**2
     +     + 64*(t2+u2)**(-1)*t2 )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,7,-2,1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 256*s**(-2)*m1**2*u2 - 128*s**(-1)*
     +    m1**2 )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,7,-2,2)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**(-2)*m1**2 )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,7,-1,1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s**(-2)*(t2+u2)**(-1)*m1**2*t2 + 384*
     +    s**(-2)*(t2+u2)**(-1)*m1**2*u2 + 256*s**(-2)*(t2+u2)**(-1)*t2
     +    *u2 + 64*s**(-2)*(t2+u2)**(-1)*t2**2 + 192*s**(-2)*
     +    (t2+u2)**(-1)*u2**2 - 128*s**(-2)*m1**2 - 64*s**(-2)*t2 - 192
     +    *s**(-2)*u2 + 256*s**(-1)*(t2+u2)**(-1)*m1**2 + 192*s**(-1)*
     +    (t2+u2)**(-1)*t2 + 64*s**(-1)*(t2+u2)**(-1)*u2 - 128*s**(-1)
     +     )
      MM_qb = MM_qb + lumi(3)*ANG_fin(5,7,-1,2)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 256*s**(-2)*(t2+u2)**(-1)*m1**2 + 128*
     +    s**(-2)*(t2+u2)**(-1)*t2 + 128*s**(-2)*(t2+u2)**(-1)*u2 - 128
     +    *s**(-2) )
      MM_qb = MM_qb + lumi(3)*ANG_fin(7,8,1,-1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * (  - 128*s**(-2)*(t2+u2)**(-1)*m1**2*t2 + 
     +    128*s**(-2)*(t2+u2)**(-1)*m1**2*u2 + 128*s**(-2)*m1**2 + 256*
     +    s**(-1)*(t2+u2)**(-1)*m1**2 + 192*s**(-1)*(t2+u2)**(-1)*t2 + 
     +    64*s**(-1)*(t2+u2)**(-1)*u2 - 128*s**(-1) )
      MM_qb = MM_qb + lumi(3)*ANG_fin(7,8,2,-1)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 256*s**(-2)*(t2+u2)**(-1)*m1**2 + 128*
     +    s**(-2)*(t2+u2)**(-1)*t2 + 128*s**(-2)*(t2+u2)**(-1)*u2 - 128
     +    *s**(-2) )
      MM_qb = MM_qb + lumi(3)*ANG_fin(8,0,-1,0)*Ar2*Nc*Cf*st**2*Pi**2*
     + alphas**2*hardfac * ( 128*s**(-1)*(t2+u2)**(-1)*m1**2*u2 + 128*
     +    s**(-1)*(t2+u2)**(-1)*m1**4 - 128*s**(-1)*m1**2 + 64*
     +    (t2+u2)**(-1)*m1**2 + 64*(t2+u2)**(-1)*t2 )

c               the prefactor for the scaling functions 
      LE_QB = MM_qb * m1**2/4.D0

c               the phase space except for 1/s**2 
      LE_QB = LE_QB / ( 16.D0 * pi**2 )**2 / 2.D0 * s4/(s4+m1**2)
      
c               the averaging factors
      LE_QB = LE_QB /4.D0 /Nc**2

      return
      end


