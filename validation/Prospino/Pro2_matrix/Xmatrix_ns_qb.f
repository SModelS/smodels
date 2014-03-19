cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NS_QBH(MASSIN,C)                                                 c
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
      real*8 function NS_QBPI2(massin,Ar,sq_in,qin_in,qbin_in,qbout_in)

      implicit none 

      integer    n1,n2
     &          ,sq_in,qin_in,qbin_in,qbout_in
     &          ,sq,qin,qbin,qbout
      real*8     massin(1:30),Pi,CF,Sn,N,gs,ms
     &          ,s,s4,m1,m2,mg,t1,u1,mj,tj,gamg
     &          ,mjs2,mgs2,coeqq,hardfac,theta_s3,theta_s4
     &          ,df,abssq,Re
      complex*16 Ar(1:5,1:5,4),s4g,qqbHs3s4g,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,Cj(1:5,1:5),Cjcc(1:5,1:5),Cjp(1:5,1:5)
     &          ,Cjf(1:5,1:5),Cjfcc(1:5,1:5),Cjpf(1:5,1:5)

      external df,abssq,Re

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      Sn = 1.D0/(16.D0 * Pi**2)
      N  = 3.D0 

c               only flavor matters, not (anti-) particle
      sq    = abs( sq_in )
      qin   = abs( qin_in ) 
      qbin  = abs( qbin_in ) 
      qbout = abs( qbout_in ) 

c               here the neutralino mass has to be positive?!?!
      s    = massin(1)
      tj   = massin(2)
      s4   = massin(3)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(10)
      ms   = massin(11)
      gamg = massin(26)
      mj   = m2

c               real kinematics built in
      t1 = tj + m2**2 - m1**2
      u1  = s4 - s - tj 
      s4g = dcmplx( s4 + m1**2 - mg**2, 0.D0 ) 

c               the s3 regularization
      theta_s3 = 0.D0
      theta_s4 = 0.D0
      if ((s.gt.(ms+abs(m1))**2).and.(ms.gt.abs(m2))) theta_s3 = 1.D0 
      if ((s.ge.(mg+abs(m2))**2).and.(mg.ge.m1))      theta_s4 = 1.D0 

      if      ( (theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0) ) then 
         s4g  = dcmplx(real(s4g),mg*gamg)
      else if ( (theta_s3.eq.0.D0).and.(theta_s4.eq.1.D0) ) then 
         s4g = ( s4g**2 + mg**2*gamg**2 )/s4g
      end if

      hardfac  = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

c               the angular functions 
      call ANGULAR_ARRAY_NS_QBPI(massin,theta_s3,NSANG)

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

c               wim's form output, replace mjs**2 -> mjs2 
c               df(sq,qin)   
c               df(sq,qbout) 
c               df(sq,qbin)  
c               df(qbin,qin)  
c               df(qbin,qbout)

c               the prefactor for the scaling functions 
      NS_QBPI2 = 0.D0

      end


c --------------------------------------------------------------------
      real*8 function NS_QBOS3(massin,Ar,sq_in,qin_in,qbin_in,qbout_in)

      implicit none 

      integer    n1,n2
     &          ,sq_in,qin_in,qbin_in,qbout_in
     &          ,sq,qin,qbin,qbout
      real*8     massin(1:30),Pi,CF,Sn,N,gs,ms
     &          ,s,s4,m1,m2,mg,s3j
     &          ,mjs2,mgs2,coeqq,hardfac,subs3fac,theta_s3
     &          ,qqbHsub,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,df,abssq,Re
      complex*16 Ar(1:5,1:5,4)
     &          ,Cj(1:5,1:5)
     &          ,Cjf(1:5,1:5)

      external df,abssq,Re

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      Sn = 1.D0/(16.D0 * Pi**2)
      N  = 3.D0 

c               only flavor matters, not (anti-) particle
      sq    = abs( sq_in )
      qin   = abs( qin_in ) 
      qbin  = abs( qbin_in ) 
      qbout = abs( qbout_in ) 

c               here the neutralino mass has to be positive?!?!
      s    = massin(1)
      s4   = massin(3)
      s3j  = massin(4)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(10)
      ms   = massin(11)
      
c               the s3 regularization
      theta_s3 = 0.D0
      if ((s.gt.(ms+abs(m1))**2).and.(ms.gt.abs(m2))) theta_s3 = 1.D0 

      subs3fac = 1.D0 + m2**2/s3j
      hardfac  = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

c               the angular functions 
      call ANGULAR_ARRAY_NS_QBOS(massin,theta_s3,NSANG)

      do n1=1,5,1 
         do n2=1,5,1
            Cj(n1,n2)    = Ar(n1,n2,1) 
            Cjf(n1,n2)   = Ar(n1,n2,3) 
         end do
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2 
      coeqq = gs**4 * Sn**2/2.D0

c               wim's form output, replace mjs**2 -> mjs2 
c               df(sq,qin)   
c               df(sq,qbout) 
c               df(sq,qbin)  
c               df(qbin,qin)  
c               df(qbin,qbout)

c               the prefactor for the scaling functions 
      NS_QBOS3 = 0.D0

      end


c --------------------------------------------------------------------
      real*8 function NS_QBOS4(massin,Ar,sq_in,qin_in,qbin_in,qbout_in)

      implicit none 

      integer    n1,n2
     &          ,sq_in,qin_in,qbin_in,qbout_in
     &          ,sq,qin,qbin,qbout
      real*8     massin(1:30),Pi,CF,Sn,N,gs,gamg
     &          ,s,s4,m1,m2,mj,mg,t1,tj,u1,s4g
     &          ,mjs2,mgs2,coeqq
     &          ,qqbHsub,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,df,abssq,Re
      complex*16 Ar(1:5,1:5,4)
     &          ,Cj(1:5,1:5),Cjp(1:5,1:5)
     &          ,Cjf(1:5,1:5),Cjpf(1:5,1:5)

      external df,abssq,Re

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      Sn = 1.D0/(16.D0 * Pi**2)
      N  = 3.D0 

c               only flavor matters, not (anti-) particle
      sq    = abs( sq_in )
      qin   = abs( qin_in ) 
      qbin  = abs( qbin_in ) 
      qbout = abs( qbout_in ) 

c               here the neutralino mass has to be positive?!?!
      s    = massin(1)
      tj   = massin(2)
      s4   = massin(3)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(10)
      gamg = massin(26)
      
      mj  = m2 

c               real kinematics built in
      t1 = tj + m2**2 - m1**2

      u1  = s4 - s - tj 
      s4g = s4 + m1**2 - mg**2 

      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

c               the s4 regularization
      if ((s.ge.(mg+abs(m2))**2).and.(mg.ge.m1)) then 
         s4g  = sqrt( s4g**2 + mg**2*gamg**2 )
      end if

c               the angular functions 
      NSANG(0,0,0,0) = 2.D0*pi

      do n1=1,5,1 
         do n2=1,5,1
            Cj(n1,n2)    = Ar(n1,n2,1) 
            Cjp(n1,n2)   = Ar(n1,n2,2) 
            Cjf(n1,n2)   = Ar(n1,n2,3) 
            Cjpf(n1,n2)  = Ar(n1,n2,4) 
         end do
      end do

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
      gs = 1.D0

c               coupling matrices, taking into account CC = AA/2 
      coeqq = gs**4 * Sn**2/2.D0

c               wim's form output, replace mjs**2 -> mjs2 
c               df(sq,qin)   
c               df(sq,qbout) 
c               df(sq,qbin)  
c               df(qbin,qin)  
c               df(qbin,qbout)

c               the prefactor for the scaling functions 
      NS_QBOS4 = 0.D0

      end


c --------------------------------------------------------------------
      real*8 function NS_QBH(massin,Ar
     &                      ,sq_in,qin_in,qbin_in,qbout_in,icall,icount)

      implicit none 

      integer    n1,n2,n3,n4,icall,icount
     &          ,sq_in,qin_in,qbin_in,qbout_in
     &          ,sq,qin,qbin,qbout
      real*8     massin(1:30),Pi,CF,Sn,N,gs,gamg
     &          ,QF,s,s4,m1,m2,mj,ms,mg,t1,tj,u1,uj,s4g
     &          ,sptj,hardfac,mjs2,mgs2,coeqq
     &          ,ujpmgs2,tjmmgs2,tjpuj
     &          ,qqbH,NSANG(0:9,0:9,-2:2,-2:2)
     &          ,NSANG1(0:9,0:9,-2:2,-2:2)
     &          ,NSANG2(0:9,0:9,-2:2,-2:2)
     &          ,df,abssq,Re,theta_s3
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
      sq    = abs( sq_in )
      qin   = abs( qin_in ) 
      qbin  = abs( qbin_in ) 
      qbout = abs( qbout_in ) 

c               here the neutralino mass has to be positive?!?!
      s    = massin(1)
      tj   = massin(2)
      s4   = massin(3)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(10)
      ms   = massin(11)
      gamg = massin(26)
      
      mj  = m2 

c               real kinematics built in
      t1 = tj + m2**2 - m1**2

      u1  = s4 - s - tj 
      uj  = u1 + m1**2 - m2**2 
      s4g = s4 + m1**2 - mg**2 

      hardfac = 1.D0 + m1**2/s4
      mjs2 = m2**2 - m1**2
      mgs2 = mg**2 - m1**2 

c               some denominators following wim's notes 
      sptj    = s + tj 
      tjpuj   = tj + uj 
      tjmmgs2 = tj - mgs2 
      ujpmgs2 = uj + mgs2

c               the s3 regularization
      theta_s3 = 0.D0
      if ((s.gt.(ms+abs(m1))**2).and.(ms.gt.abs(m2))) theta_s3 = 1.D0 

c               the explicit s4g denominator [hauptwert]
      if ((s.ge.(mg+abs(m2))**2).and.(mg.ge.m1)) then 
         s4g = ( s4g**2 + mg**2*gamg**2 )/s4g
      end if

c               the angular functions 
      if (icall.eq.0) then 
         if (icount.eq.0) 
     &      call ANGULAR_ARRAY_NS_QB(massin,theta_s3,NSANG)
      else if (icall.eq.1) then 
         if (icount.eq.0) 
     &      call ANGULAR_ARRAY_NS_QB(massin,theta_s3,NSANG1)
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
         if (icount.eq.0) 
     &      call ANGULAR_ARRAY_NS_QB(massin,theta_s3,NSANG2)
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

c               wim's form output, replace mjs**2 -> mjs2 
c               df(sq,qin)   
c               df(sq,qbout) 
c               df(sq,qbin)  
c               df(qbin,qin)  
c               df(qbin,qbout)

c               the prefactor for the scaling functions 
      NS_QBH = 0.D0

      end






