cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NS_GGH(MASSIN,C)                                                 c
c     NS_GGOS(MASSIN,C)                                                c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = tg                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = mg                                                c
c       MASSIN(9)  = mt                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(12) = qr                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(25) = gamma_s                                           c
c       MASSIN(26) = yesno_s4                                          c
c       MASSIN(27) = yesno_s3                                          c
c       MASSIN(30) = kinematic factor s4/(s4-m1**2)                    c
c       MASSIN(4-5,8,10,13-30) not needed                              c
c                                                                      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function NS_GGOS(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,gs,theta_s3
     &          ,s,s4,m1,m2,ms,s3j,hardfac,mjs2,coegg
     &          ,subs3fac,ggHsub,NSANG(0:9,0:9,-2:2,-2:2)
      complex*16 Ar(4),Arc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               here the neutralino mass has to be positive
      s   = massin(1)
      s4  = massin(3)
      s3j = massin(4)
      m1  = massin(6)
      m2  = massin(7)
      ms  = massin(11)

c               the s3 regularization
      theta_s3 = 0.D0
      if ((s.gt.(ms+abs(m1))**2).and.(ms.gt.abs(m2))) theta_s3 = 1.D0 

c               real kinematics built in
      subs3fac = 1.D0 + m2**2/s3j
      hardfac  = 1.D0 + m1**2/s4

      mjs2 = m2**2 - m1**2

c               the angular functions 
      call ANGULAR_ARRAY_NS_GGOS(massin,theta_s3,NSANG)

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2
      coegg = gs**4 * real( Ar(1)*Arc(1) )/2.D0 * Sn**2/4.D0

c               wim's form output
c               the prefactor for the scaling functions 
      NS_GGOS = 0.D0

      end




c --------------------------------------------------------------------
      real*8 function NS_GGH(massin,Ar)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,gs,coegg
     &          ,QF,s,s4,m1,m2,mjs2,ms,tj,t1,uj,u1,t1hat,u1hat
     &          ,sptj,spuj,tjpuj,hardfac
     &          ,ggH,NSANG(0:9,0:9,-2:2,-2:2)
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
      t1hat = -( tj*uj - mjs2*s )/(s+tj)
      u1hat = -( tj*uj - mjs2*s )/(s+uj)
      sptj  = s + tj
      spuj  = s + uj
      tjpuj = tj + uj 

c               the angular functions 
      call ANGULAR_ARRAY_NS_GG(massin,NSANG)

c               the factorization/renormalization scale 
      QF = massin(13) 

      do n=1,4 
         Arc(n) =  conjg( Ar(n) )
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2
      coegg = gs**4 * real( Ar(1)*Arc(1) )/2.D0 * Sn**2/4.D0

c               wim's form output

c               the prefactor for the scaling functions 
      NS_GGH = 0.D0

      end






