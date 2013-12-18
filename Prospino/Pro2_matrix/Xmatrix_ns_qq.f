cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NNS_QH(MASSIN,C)                                                 c
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
      real*8 function NS_QQOS(massin,Ar,sq_in,q1in_in,q2in_in,qout_in)

      implicit none 

      integer    n1,n2,sq_in,q1in_in,q2in_in,qout_in
     &          ,sq,q1in,q2in,qout
      real*8     massin(1:30),Pi,CF,N,Sn,gs,theta_s3
     &          ,s,s4,m1,m2,ms,s3j,hardfac,mjs2,mgs2,coeqq
     &          ,mg,mgs,sfac1 
     &          ,subs3fac,qqHsub,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,df,abssq,Re
      complex*16 Ar(1:5,1:5,4)
     &          ,Cjp(1:5,1:5)
     &          ,Cjpf(1:5,1:5)

      save NSANG

      external df,abssq,Re

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      N  = 3.D0 
      Sn = 1.D0/(16.D0 * Pi**2)

c               only flavor matters, not (anti-) particle
      sq   = abs( sq_in )
      q1in = abs( q1in_in ) 
      q2in = abs( q2in_in ) 
      qout = abs( qout_in ) 

c               here the neutralino mass has to be positive?!?!
      s   = massin(1)
      s4  = massin(3)
      s3j = massin(4)
      m1  = massin(6)
      m2  = massin(7)
      mg  = massin(10)
      ms  = massin(11)

c               the s3 regularization
      theta_s3 = 0.D0
      if ((s.gt.(ms+abs(m1))**2).and.(ms.gt.abs(m2))) theta_s3 = 1.D0 

c               real kinematics built in
      subs3fac = 1.D0 + m2**2/s3j
      hardfac  = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

      sfac1 = s + 2.D0*mgs2 

c               these appear only with even powers higher than 2
      mgs  = sqrt( abs(mgs2) )

c               the angular functions, called twice, no icount
      call ANGULAR_ARRAY_NS_QQOS(massin,theta_s3,NSANG)

      do n1=1,5,1 
         do n2=1,5,1
            Cjp(n1,n2)   = Ar(n1,n2,2) 
            Cjpf(n1,n2)  = Ar(n1,n2,4) 
         end do
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2 
      coeqq = gs**4 * Sn**2/2.D0


c               the prefactor for the scaling functions 
      NS_QQOS = 0.D0

      end


c --------------------------------------------------------------------
      real*8 function NS_QQH(massin,Ar
     &                      ,sq_in,q1in_in,q2in_in,qout_in,icall,icount)

      implicit none 

      integer    n1,n2,n3,n4,icall,icount
     &          ,sq_in,q1in_in,q2in_in,qout_in
     &          ,sq,q1in,q2in,qout
      real*8     massin(1:30),Pi,CF,Sn,N,gs
     &          ,QF,s,s4,m1,m2,mj,ms,mg,t1,tj,u1,uj,ug,tg
     &          ,sptj,hardfac,mjs2,mgs2,coeqq,sfac1,sfac2 
     &          ,ujpmgs2,tjmmgs2,ujmmgs2,tjpmgs2,tjpuj,s4mug,s4mtg
     &          ,spuj
     &          ,qqH,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,NSANG1(0:9,0:9,-2:2,-2:2)
     &          ,NSANG2(0:9,0:9,-2:2,-2:2)
     &          ,df,abssq,Re
      complex*16 Ar(1:5,1:5,4)
     &          ,Cj(1:5,1:5),Cjcc(1:5,1:5),Cjp(1:5,1:5)
     &          ,Cjf(1:5,1:5),Cjfcc(1:5,1:5),Cjpf(1:5,1:5)
      
      save NSANG
      save NSANG1
      save NSANG2

      external df,abssq,Re
      
      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      Sn = 1.D0/(16.D0 * Pi**2)
      N  = 3.D0 

c               only flavor matters, not (anti-) particle
      sq   = abs( sq_in )
      q1in = abs( q1in_in ) 
      q2in = abs( q2in_in ) 
      qout = abs( qout_in ) 

c               here the neutralino mass has to be positive?!?!
      s    = massin(1)
      tj   = massin(2)
      s4   = massin(3)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(10)
      ms   = massin(11)
      
      mj  = m2 

c               real kinematics built in
      t1 = tj + m2**2 - m1**2
      tg = tj + m2**2 - mg**2 

      u1  = s4 - s - tj 
      uj  = u1 + m1**2 - m2**2 
      ug  = u1 + m1**2 - mg**2

      hardfac = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

c               some denominators following wim's notes 
      s4mug   = s4 - ug 
      s4mtg   = s4 - tg 
      sptj    = s + tj 
      spuj    = s + uj 
      tjpuj   = tj + uj 
      tjmmgs2 = tj - mgs2 
      tjpmgs2 = tj + mgs2 
      ujmmgs2 = uj - mgs2
      ujpmgs2 = uj + mgs2

      sfac1 = s + 2.D0*mgs2 
      sfac2 = s + 2.D0*mgs2 - mjs2 

c               the angular functions 
      if (icall.eq.0) then 
         if (icount.eq.0) call ANGULAR_ARRAY_NS_QQ(massin,NSANG)
      else if (icall.eq.1) then 
         if (icount.eq.0) call ANGULAR_ARRAY_NS_QQ(massin,NSANG1)
         do n1=0,9
            do n2=0,9
               do n3=-2,2
                  do n4=-2,2
                     NSANG(n1,n2,n3,n4) = NSANG1(n1,n2,n3,n4)
                  end do
               end do
            end do
         end do
      else if (icall.eq.2) then 
         if (icount.eq.0) call ANGULAR_ARRAY_NS_QQ(massin,NSANG2)
         do n1=0,9
            do n2=0,9
               do n3=-2,2
                  do n4=-2,2
                     NSANG(n1,n2,n3,n4) = NSANG2(n1,n2,n3,n4)
                  end do
               end do
            end do
         end do
      end if 

c               the factorization/renormalization scale 
      QF = massin(13) 

      do n1=1,5,1 
         do n2=1,5,1
            Cj(n1,n2)    = Ar(n1,n2,1) 
            Cjp(n1,n2)   = Ar(n1,n2,2) 
            Cjf(n1,n2)   = Ar(n1,n2,3) 
            Cjpf(n1,n2)  = Ar(n1,n2,4) 
            Cjcc(n1,n2)  = conjg( Cj(n1,n2) )
            Cjfcc(n1,n2) = conjg( Cjf(n1,n2) )
         end do
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2 
      coeqq = gs**4 * Sn**2/2.D0

c               wim's form output
c               df() with arguments sq,q1in,q2in,qout

c               the prefactor for the scaling functions 
      NS_QQH = 0.D0

      end


c-----------------------------------------------------------------
c some couplings a la wim numerically 
      real*8 function abssq(Cj,i1,i2)

      implicit none 

      integer    i1,i2
      complex*16 Cj(1:5,1:5)

      abssq = abs( Cj(i1,i2) )**2 

      end 


      real*8 function Re(Cj1,i1,i2,Cj2,i3,i4) 
      
      implicit none 

      integer    i1,i2,i3,i4
      complex*16 Cj1(1:5,1:5),Cj2(1:5,1:5)

      Re = real( Cj1(i1,i2) * Cj2(i3,i4) )

      end 


c-----------------------------------------------------------------
c numerical kronecker delta, only affects the absolute values
c   i.e. checks for flavor, not for charge
      real*8 function df(i1,i2) 
      
      implicit none 

      integer i1,i2

      if ( abs(i1).eq.abs(i2) ) then 
         df = 1.D0
      else 
         df = 0.D0 
      end if

      end 








