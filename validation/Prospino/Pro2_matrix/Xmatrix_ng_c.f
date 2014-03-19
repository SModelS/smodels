cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NG_QGH(MASSIN,CL,CR)                                             c
c     NG_QGOS1(MASSIN,CL,CR)                                           c
c     NG_QGOS2(MASSIN,CL,CR)                                           c
c     NG_GBH(MASSIN,CL,CR)                                             c
c     NG_GBOS1(MASSIN,CL,CR)                                           c
c     NG_GBOS2(MASSIN,CL,CR)                                           c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = tg                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(4)  = s3                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = mg                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(13) = qf                                                c
c       MASSIN(25) = gamma_s                                           c
c       MASSIN(26) = yesno_s4                                          c
c       MASSIN(27) = yesno_s3                                          c
c                                                                      c
c       CL(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CR(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function NG_QGOS1(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc1,struc2,gams
     &          ,eight,s,tg,s4,s3,mi,mg,ms,ui,us,s4i,hardfac,s4s
     &          ,mgs2,mis2,IMANG(0:9,0:9,-2:2,-2:2)
     &          ,subs3fac,qgHsub,NGANG(0:9,0:9,-2:2,-2:2)
     &          ,theta_s3,theta_s4,yesno_s4,yesno_s3,dum
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      s3   = massin(4)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)
      yesno_s4 = massin(26)
      yesno_s3 = massin(27)
      
c               the s4 regularization
      theta_s4 = 0.D0
      s4s     = s4 + mi**2 - ms**2 
      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         theta_s4 = 1.D0
         dum = s4s**2 + ms**2*gams**2
         if (dum.ge.0.D0) then 
            s4s = sign(1.D0,s4s) * sqrt( dum )
         else 
            print*, " NG_QGOS1: serious problem with s4s ",s4s
            call HARD_STOP 
         end if 
      end if 

c               the s3 subtraction flag
      theta_s3 = 0.D0
      if ((s.ge.(ms+abs(mi))**2).and.(ms.ge.abs(mg))) then 
         theta_s3 = 1.D0
      end if 

c               real kinematics built in
      ui = s4 - s - tg 
      us = ui + mi**2 - ms**2 

      s4i     = s4
      if (s4i.lt.0.D0) s4i = 0.D0

c               some factors from the form file
      mis2 = mi**2 - ms**2 
      mgs2 = mg**2 - ms**2 

c               wim's coupling structures 
      struc1 = real( Cl(1)*Clc(1) + Cr(1)*Crc(1) )
      struc2 = real( Cl(3)*Clc(3) + Cr(3)*Crc(3) )

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 re-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the angular functions 
      call ANGULAR_ARRAY_NG_S(1,theta_s3,theta_s4,massin,NGANG,IMANG)

c               wim's form output
c      (ii) replaced: mgs**2 -> mgs2; mis**2 -> mis2
      qgHsub = 0.D0
      if (yesno_s4.eq.1.D0) then 
      qgHsub = qgHsub +  
     + NGANG(0,0,0,0)*s4s**(-2)*e**2*gs**4*struc2*CA**(-1)*CF*Sn**2*
     + eight**(-1) * (  - 8*s**(-1)*tg*ms**(-2)*mis2**2 + 16*s**(-1)*
     +    us**(-1)*ms**(-2)*mis2**2*mgs2**2 
     +    - 16*s**(-1)*ms**(-2)*mis2**2*
     +    mgs2 - 16*us**(-2)*mis2**2*mgs2 - 16*us**(-1)*ms**(-2)*
     +    mis2**2*mgs2 )
      qgHsub = qgHsub + NGANG(0,0,0,0)*s4s**(-2)*e**2*gs**4*struc2*
     + Sn**2*eight**(-1) * (  - 8*s**(-1)*tg**(-1)*ms**(-2)*mis2**2*
     +    mgs2**2 - 8*s**(-1)*us**(-1)*ms**(-2)*mis2**2*mgs2**2 - 8*s*
     +    tg**(-1)*ms**(-2)*mis2**2 
     +    - 16*tg**(-2)*ms**(-2)*mis2**2*mgs2**2
     +     - 16*tg**(-2)*mis2**2*mgs2 - 8*tg**(-1)*us**(-1)*ms**(-2)*
     +    mis2**2*mgs2**2 - 16*tg**(-1)*us**(-1)*mis2**2*mgs2 - 16*
     +    tg**(-1)*ms**(-2)*mis2**2*mgs2 - 8*ms**(-2)*mis2**2 )
      end if 
      if (yesno_s3.eq.1.D0) then 
      hardfac  = 1.D0 + mi**2/s4i
      subs3fac = 1.D0 + mg**2/s3
      qgHsub = qgHsub + NGANG(4,0,-2,0)*e**2*gs**4*struc1*CA**(-1)*CF*
     + Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 16*s**(-1)*
     +    ms**(-2)*mis2*mgs2**2 )
      qgHsub = qgHsub + NGANG(4,6,-2,1)*e**2*gs**4*struc1*CA**(-1)*CF*
     + Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 8*s**(-1)*
     +    ms**(-2)*mgs2**2 )
      qgHsub = qgHsub + NGANG(4,9,-2,-2)*e**2*gs**4*struc1*CA**(-1)*CF*
     + Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 16*mis2*mgs2**2
     +     )
      qgHsub = qgHsub + NGANG(4,9,-2,-1)*e**2*gs**4*struc1*CA**(-1)*CF*
     + Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * ( 16*s**(-1)*ms**(-2)
     +    *mis2**2*mgs2**2 - 16*ms**(-2)*mis2*mgs2**2 )
      end if

c               the usual kinematical factor 
c      qgHsub = qgHsub * kin

c               the prefactor for the scaling functions 
      NG_QGOS1 = qgHsub * (abs(mi)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_GBOS1(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc1,struc2,gams,ts
     &          ,eight,s,tg,s3,s4,mi,mg,ms,ug,ui,s4i,hardfac,s4s
     &          ,mgs2,mis2,subs3fac,gqbHsub,NGANG(0:9,0:9,-2:2,-2:2)
     &          ,yesno_s4,yesno_s3,theta_s3,theta_s4
     &          ,IMANG(0:9,0:9,-2:2,-2:2),dum
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      s3   = massin(4)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)
      yesno_s4 = massin(26)
      yesno_s3 = massin(27)
      
c               the s4 regularization
      theta_s4 = 0.D0
      s4s     = s4 + mi**2 - ms**2 
      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         theta_s4 = 1.D0
         dum = s4s**2 + ms**2*gams**2
         if (dum.ge.0.D0) then 
            s4s = sign(1.D0,s4s) * sqrt( dum )
         else 
            print*, " NG_GBOS1: serious problem with s4s ",s4s
            call HARD_STOP 
         end if 
      end if 

c               the s3 subtraction flag
      theta_s3 = 0.D0
      if ((s.ge.(ms+abs(mi))**2).and.(ms.ge.abs(mg))) then 
         theta_s3 = 1.D0
      end if 

c               real kinematics built in
      ts = tg + mg**2 - ms**2 

      ui = s4 - s - tg 
      ug = ui + mi**2 - mg**2 

      s4i     = s4
      if (s4i.lt.0.D0) s4i = 0.D0

c               some factors from the form file
      mis2 = mi**2 - ms**2 
      mgs2 = mg**2 - ms**2 

c               wim's coupling structures 
      struc1 = real( Cl(1)*Clc(1) + Cr(1)*Crc(1) )
      struc2 = real( Cl(3)*Clc(3) + Cr(3)*Crc(3) )

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 re-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the angular functions 
c               for theta_s3=1 choose the one dimensional integrals  
      call ANGULAR_ARRAY_NG_S(1,theta_s3,theta_s4,massin,NGANG,IMANG)

c               wim's form output
c      (ii) replaced: mgs**2 -> mgs2; mis**2 -> mis2
      gqbHsub = 0.D0
      if (yesno_s4.eq.1.D0) then 
      gqbHsub = gqbHsub +  
     + NGANG(0,0,0,0)*s4s**(-2)*e**2*gs**4*struc1*CA**(-1)*CF*Sn**2*
     + eight**(-1) * (  - 8*s**(-1)*ug*ms**(-2)*mis2**2 + 16*s**(-1)*
     +    ts**(-1)*ms**(-2)*mis2**2*mgs2**2 
     +    - 16*s**(-1)*ms**(-2)*mis2**2*
     +    mgs2 - 16*ts**(-2)*mis2**2*mgs2 - 16*ts**(-1)*ms**(-2)*
     +    mis2**2*mgs2 )
      gqbHsub = gqbHsub + NGANG(0,0,0,0)*s4s**(-2)*e**2*gs**4*struc1*
     + Sn**2*eight**(-1) * (  - 8*s**(-1)*ug**(-1)*ms**(-2)*mis2**2*
     +    mgs2**2 - 8*s**(-1)*ts**(-1)*ms**(-2)*mis2**2*mgs2**2 - 8*s*
     +    ug**(-1)*ms**(-2)*mis2**2 
     +    - 16*ug**(-2)*ms**(-2)*mis2**2*mgs2**2
     +     - 16*ug**(-2)*mis2**2*mgs2 - 8*ug**(-1)*ts**(-1)*ms**(-2)*
     +    mis2**2*mgs2**2 - 16*ug**(-1)*ts**(-1)*mis2**2*mgs2 - 16*
     +    ug**(-1)*ms**(-2)*mis2**2*mgs2 - 8*ms**(-2)*mis2**2 )
      end if 
      if (yesno_s3.eq.1.D0) then 
      hardfac  = 1.D0 + mi**2/s4i
      subs3fac = 1.D0 + mg**2/s3
      gqbHsub = gqbHsub + NGANG(4,0,-2,0)*e**2*gs**4*struc2*CA**(-1)*CF
     + *Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 16*s**(-1)*
     +    ms**(-2)*mis2*mgs2**2 )
      gqbHsub = gqbHsub + NGANG(4,7,-2,1)*e**2*gs**4*struc2*CA**(-1)*CF
     + *Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 8*s**(-1)*
     +    ms**(-2)*mgs2**2 )
      gqbHsub = gqbHsub + NGANG(4,8,-2,-2)*e**2*gs**4*struc2*CA**(-1)*
     + CF*Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * (  - 16*mis2*
     +    mgs2**2 )
      gqbHsub = gqbHsub + NGANG(4,8,-2,-1)*e**2*gs**4*struc2*CA**(-1)*
     + CF*Sn**2*hardfac**(-1)*subs3fac*eight**(-1) * ( 16*s**(-1)*
     +    ms**(-2)*mis2**2*mgs2**2 - 16*ms**(-2)*mis2*mgs2**2 )
      end if

c               the usual kinematical factor 
c      gqbHsub = gqbHsub * kin

c               the prefactor for the scaling functions 
      NG_GBOS1 = gqbHsub * (abs(mi)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_QGOS2(massin,Cl,Cr)

      implicit none 

      integer    n,n1,n2,n3,n4
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc3,gams
     &          ,eight,s,tg,s4,mi,mg,ms,ui,us,ug,s4i,hardfac
     &          ,RGANG(0:9,0:9,-2:2,-2:2),IMANG(0:9,0:9,-2:2,-2:2)
     &          ,theta_s3,theta_s4
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4),qgHs3s4,s4s,s4s_d
     &          ,NGANG(0:9,0:9,-2:2,-2:2)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)

c               real kinematics built in
      ui = s4 - s - tg 
      us = ui + mi**2 - ms**2 
      ug = ui + mi**2 - mg**2 

c               the denominator regularized 
      s4i     = s4
      if (s4i.lt.0.D0) s4i = 0.D0

      hardfac = 1.D0 + mi**2/s4i

c               wim's coupling structures 
      struc3 = mi*mg*real(  Cl(1)*Cl(3) + Clc(1)*Clc(3) 
     &                     + Cr(1)*Cr(3) + Crc(1)*Crc(3) ) 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 re-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the s4 regularization
      theta_s4 = 0.D0
      s4s  = dcmplx( s4+mi**2-ms**2 )
      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         theta_s4 = 1.D0
         s4s_d  = s4s/( s4s**2 + ms**2*gams**2 )
      else 
         s4s_d = 1.D0/s4s
      end if

c               the s3 subtraction flag
      theta_s3 = 0.D0
      if ((s.ge.(ms+abs(mi))**2).and.(ms.ge.abs(mg))) then 
         theta_s3 = 1.D0
      end if 

      if ((theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0)) then
         s4s     = dcmplx(s4+mi**2-ms**2,ms*gams)  
      end if

c               the angular functions 
      call ANGULAR_ARRAY_NG_S(2,theta_s3,theta_s4,massin,RGANG,IMANG)
      
      do n1=0,9,1
         do n2=0,9,1
            do n3=-2,2,1
               do n4=-2,2,1
                  NGANG(n1,n2,n3,n4) = RGANG(n1,n2,n3,n4)
                  if ((theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0)) then                     
                     NGANG(n1,n2,n3,n4) = dcmplx(RGANG(n1,n2,n3,n4)
     &                                          ,IMANG(n1,n2,n3,n4))
                  end if
               end do
            end do
         end do
      end do

c               wim's form output
      qgHs3s4 = 
     +  NGANG(4,0,-1,0)*s4s_d*e**2*gs**4*struc3*CA**(-1)*CF*Sn**2
     + *hardfac**(-1)*eight**(-1) * ( 4 + 8*s**(-1)*tg*ug*us**(-1) - 4*
     +    s**(-1)*tg + 4*s**(-1)*ug*us**(-1)*ms**2 - 4*s**(-1)*ug*
     +    us**(-1)*mi**2 + 4*s**(-1)*ug + 8*s**(-1)*ug**2*us**(-1) + 12
     +    *s*us**(-1) + 22*ug*us**(-1) - 6*ui*us**(-1) + 4*us**(-1)*
     +    ms**2 - 14*us**(-1)*mi**2 - 6*us**(-1)*mg**2 )
      qgHs3s4 = qgHs3s4 + NGANG(4,0,-1,0)*s4s_d*e**2*gs**4*struc3*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 2 - 2*s**(-1)*tg**(-1)*ug
     +    *ms**2 + 2*s**(-1)*tg**(-1)*ug*mi**2 - 4*s**(-1)*tg**(-1)*
     +    ug**2 - 4*s**(-1)*tg*ug*us**(-1) - 2*s**(-1)*ug*us**(-1)*
     +    ms**2 + 2*s**(-1)*ug*us**(-1)*mi**2 - 4*s**(-1)*ug - 4*
     +    s**(-1)*ug**2*us**(-1) - 4*s*tg**(-1) - 6*s*us**(-1) - 8*
     +    tg**(-1)*ug + 2*tg**(-1)*ui - 2*tg**(-1)*ms**2 + 6*tg**(-1)*
     +    mi**2 + 4*tg**(-1)*mg**2 - 11*ug*us**(-1) + 3*ui*us**(-1) - 2
     +    *us**(-1)*ms**2 + 7*us**(-1)*mi**2 + 3*us**(-1)*mg**2 )
      qgHs3s4 = qgHs3s4 + NGANG(4,7,-1,1)*s4s_d*e**2*gs**4*struc3*
     + CA**(-1)*CF*Sn**2*hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*tg*
     +    us**(-1) + 4*s**(-1)*ug*us**(-1) + 8*s**(-1) + 8*us**(-1) )
      qgHs3s4 = qgHs3s4 + NGANG(4,7,-1,1)*s4s_d*e**2*gs**4*struc3*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 2*s**(-1)*tg**(-1)*ug - 2
     +    *s**(-1)*tg*us**(-1) - 2*s**(-1)*ug*us**(-1) - 2*s**(-1) - 4*
     +    us**(-1) )
      qgHs3s4 = qgHs3s4 + NGANG(4,9,-1,-1)*s4s_d*e**2*gs**4*struc3*
     + CA**(-1)*CF*Sn**2*hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*tg*
     +    ms**2 - 4*s**(-1)*tg*mi**2 + 8*s**(-1)*ug*ms**2 - 8*s**(-1)*
     +    ug*mi**2 - 16*s**(-1)*ms**2*mi**2 + 8*s**(-1)*ms**4 + 8*
     +    s**(-1)*mi**4 - 4*s*tg*us**(-1) + 8*s*ug*us**(-1) - 4*s*ui*
     +    us**(-1) + 4*s*us**(-1)*ms**2 - 16*s*us**(-1)*mi**2 - 4*s*
     +    us**(-1)*mg**2 + 4*s + 4*s**2*us**(-1) - 4*tg*ug*us**(-1) - 4
     +    *tg*us**(-1)*ms**2 + 4*tg*us**(-1)*mi**2 + 4*tg - 4*ug*ui*
     +    us**(-1) + 4*ug*us**(-1)*ms**2 - 16*ug*us**(-1)*mi**2 - 4*ug*
     +    us**(-1)*mg**2 + 2*ug + 4*ug**2*us**(-1) - 4*ui*us**(-1)*
     +    ms**2 + 4*ui*us**(-1)*mi**2 + 2*ui - 12*us**(-1)*ms**2*mi**2
     +     - 4*us**(-1)*ms**2*mg**2 + 4*us**(-1)*mi**2*mg**2 + 12*
     +    us**(-1)*mi**4 + 8*ms**2 - 10*mi**2 + 2*mg**2 )
      qgHs3s4 = qgHs3s4 + NGANG(4,9,-1,-1)*s4s_d*e**2*gs**4*struc3*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 3*s*tg**(-1)*ug + 3*s*
     +    tg**(-1)*ui + 2*s*tg**(-1)*ms**2 + 3*s*tg**(-1)*mi**2 + 3*s*
     +    tg**(-1)*mg**2 + 2*s*tg*us**(-1) - 4*s*ug*us**(-1) + 2*s*ui*
     +    us**(-1) - 2*s*us**(-1)*ms**2 + 8*s*us**(-1)*mi**2 + 2*s*
     +    us**(-1)*mg**2 + 4*s - 2*s**2*us**(-1) + 2*tg**(-1)*ug*ui + 2
     +    *tg**(-1)*ug*ms**2 + 4*tg**(-1)*ug*mi**2 + 2*tg**(-1)*ug*
     +    mg**2 - 2*tg**(-1)*ug**2 + 2*tg**(-1)*ui*ms**2 - 2*tg**(-1)*
     +    ui*mi**2 - 2*tg**(-1)*ms**2*mi**2 + 2*tg**(-1)*ms**2*mg**2 + 
     +    4*tg**(-1)*ms**4 - 2*tg**(-1)*mi**2*mg**2 - 2*tg**(-1)*mi**4
     +     + 2*tg*ug*us**(-1) + 2*tg*us**(-1)*ms**2 - 2*tg*us**(-1)*
     +    mi**2 + 2*ug*ui*us**(-1) - 2*ug*us**(-1)*ms**2 + 8*ug*
     +    us**(-1)*mi**2 + 2*ug*us**(-1)*mg**2 + 2*ug - 2*ug**2*
     +    us**(-1) + 2*ui*us**(-1)*ms**2 - 2*ui*us**(-1)*mi**2 + 6*
     +    us**(-1)*ms**2*mi**2 + 2*us**(-1)*ms**2*mg**2 - 2*us**(-1)*
     +    mi**2*mg**2 )
      qgHs3s4 = qgHs3s4 + NGANG(4,9,-1,-1)*s4s_d*e**2*gs**4*struc3*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 6*us**(-1)*mi**4 + 4*
     +    ms**2 - 4*mi**2 )

c               the prefactor for the scaling functions 
      NG_QGOS2 = real(qgHs3s4) * (abs(mi)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_GBOS2(massin,Cl,Cr)

      implicit none 

      integer    n,n1,n2,n3,n4
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc3,gams
     &          ,eight,s,tg,s4,mi,mg,ms,ui,ts,ug,s4i,hardfac
     &          ,RGANG(0:9,0:9,-2:2,-2:2),IMANG(0:9,0:9,-2:2,-2:2)
     &          ,theta_s3,theta_s4
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4),gqbHs3s4,s4s,s4s_d
     &          ,NGANG(0:9,0:9,-2:2,-2:2)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)

c               real kinematics built in
      ts = tg + mg**2 - ms**2 

      ui = s4 - s - tg 
      ug = ui + mi**2 - mg**2 

c               the denominator regularized 
      s4i     = s4
      if (s4i.lt.0.D0) s4i = 0.D0

      hardfac = 1.D0 + mi**2/s4i

c               wim's coupling structures 
      struc3 = mi*mg*real(  Cl(1)*Cl(3) + Clc(1)*Clc(3) 
     &                     + Cr(1)*Cr(3) + Crc(1)*Crc(3) ) 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 re-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the s4 regularization
      theta_s4 = 0.D0
      s4s  = dcmplx( s4+mi**2-ms**2 )

      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         theta_s4 = 1.D0
         s4s_d  = s4s/( s4s**2 + ms**2*gams**2 )
      else 
         s4s_d = 1.D0/s4s
      end if

c               the s3 subtraction flag
      theta_s3 = 0.D0
      if ((s.ge.(ms+abs(mi))**2).and.(ms.ge.abs(mg))) then 
         theta_s3 = 1.D0
      end if 

      if ((theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0)) then
         s4s     = dcmplx(s4+mi**2-ms**2,ms*gams)  
      end if

c               the angular functions 
      call ANGULAR_ARRAY_NG_S(2,theta_s3,theta_s4,massin,RGANG,IMANG)

      do n1=0,9,1
         do n2=0,9,1
            do n3=-2,2,1
               do n4=-2,2,1
                  NGANG(n1,n2,n3,n4) = RGANG(n1,n2,n3,n4)
                  if ((theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0)) then
                     NGANG(n1,n2,n3,n4) = dcmplx(RGANG(n1,n2,n3,n4)
     &                                          ,IMANG(n1,n2,n3,n4))
                  end if
               end do
            end do
         end do
      end do

c               wim's form output
      gqbHs3s4 = 
     +  NGANG(4,0,-1,0)*s4s_d*e**2*gs**4*struc3*CA**(-1)*CF*Sn**2
     + *hardfac**(-1)*eight**(-1) * (  - 4 + 8*s**(-1)*tg*ug*ts**(-1)
     +     + 8*s**(-1)*tg*ts**(-1)*ms**2 - 8*s**(-1)*tg*ts**(-1)*mi**2
     +     + 4*s**(-1)*tg + 8*s**(-1)*tg**2*ts**(-1) + 4*s**(-1)*ug*
     +    ts**(-1)*ms**2 - 4*s**(-1)*ug*ts**(-1)*mi**2 - 4*s**(-1)*ug
     +     + 8*s**(-1)*ms**2 - 8*s**(-1)*mi**2 + 4*s*ts**(-1) + 12*tg*
     +    ts**(-1) + 2*ug*ts**(-1) - 6*ui*ts**(-1) + 12*ts**(-1)*ms**2
     +     - 22*ts**(-1)*mi**2 - 6*ts**(-1)*mg**2 )
      gqbHs3s4 = gqbHs3s4 + NGANG(4,0,-1,0)*s4s_d*e**2*gs**4*struc3
     + *Sn**2*hardfac**(-1)*eight**(-1) * (  - 2 - 4*s**(-1)*tg*
     +    ug**(-1)*ms**2 + 4*s**(-1)*tg*ug**(-1)*mi**2 - 4*s**(-1)*tg*
     +    ug*ts**(-1) - 4*s**(-1)*tg*ts**(-1)*ms**2 + 4*s**(-1)*tg*
     +    ts**(-1)*mi**2 - 4*s**(-1)*tg - 4*s**(-1)*tg**2*ug**(-1) - 4*
     +    s**(-1)*tg**2*ts**(-1) - 2*s**(-1)*ug*ts**(-1)*ms**2 + 2*
     +    s**(-1)*ug*ts**(-1)*mi**2 - 2*s**(-1)*ms**2 + 2*s**(-1)*mi**2
     +     - 4*s*ug**(-1) - 2*s*ts**(-1) - 4*tg*ug**(-1) - 6*tg*
     +    ts**(-1) + 2*ug**(-1)*ui - 2*ug**(-1)*ms**2 + 6*ug**(-1)*
     +    mi**2 + 4*ug**(-1)*mg**2 - ug*ts**(-1) + 3*ui*ts**(-1) - 6*
     +    ts**(-1)*ms**2 + 11*ts**(-1)*mi**2 + 3*ts**(-1)*mg**2 )
      gqbHs3s4 = gqbHs3s4 + NGANG(4,7,-1,1)*s4s_d*e**2*gs**4*struc3
     + *CA**(-1)*CF*Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg
     +    *ts**(-1) - 4*s**(-1)*ug*ts**(-1) - 8*s**(-1) - 8*ts**(-1) )
      gqbHs3s4 = gqbHs3s4 + NGANG(4,7,-1,1)*s4s_d*e**2*gs**4*struc3
     + *Sn**2*hardfac**(-1)*eight**(-1) * ( 2*s**(-1)*tg*ug**(-1) + 2*
     +    s**(-1)*tg*ts**(-1) + 2*s**(-1)*ug*ts**(-1) + 2*s**(-1) + 4*
     +    ts**(-1) )
      gqbHs3s4 = gqbHs3s4 + NGANG(4,8,-1,-1)*s4s_d*e**2*gs**4*
     + struc3*CA**(-1)*CF*Sn**2*hardfac**(-1)*eight**(-1) * ( 8*s**(-1)
     +    *tg*ms**2 - 8*s**(-1)*tg*mi**2 + 4*s**(-1)*ug*ms**2 - 4*
     +    s**(-1)*ug*mi**2 - 16*s**(-1)*ms**2*mi**2 + 8*s**(-1)*ms**4
     +     + 8*s**(-1)*mi**4 + 4*s*tg*ts**(-1) - 4*s*ui*ts**(-1) + 4*s*
     +    ts**(-1)*ms**2 - 16*s*ts**(-1)*mi**2 - 4*s*ts**(-1)*mg**2 + 4
     +    *s + 4*s**2*ts**(-1) - 4*tg*ui*ts**(-1) - 12*tg*ts**(-1)*
     +    mi**2 - 4*tg*ts**(-1)*mg**2 + 4*tg + 2*ug - 4*ui*ts**(-1)*
     +    ms**2 + 4*ui*ts**(-1)*mi**2 + 2*ui - 12*ts**(-1)*ms**2*mi**2
     +     - 4*ts**(-1)*ms**2*mg**2 + 4*ts**(-1)*mi**2*mg**2 + 12*
     +    ts**(-1)*mi**4 + 8*ms**2 - 10*mi**2 + 2*mg**2 )
      gqbHs3s4 = gqbHs3s4 + NGANG(4,8,-1,-1)*s4s_d*e**2*gs**4*
     + struc3*Sn**2*hardfac**(-1)*eight**(-1) * (  - 2*s*tg*ts**(-1) + 
     +    3*s*ug**(-1)*ui + 2*s*ug**(-1)*ms**2 + 3*s*ug**(-1)*mi**2 + 3
     +    *s*ug**(-1)*mg**2 + 2*s*ui*ts**(-1) - 2*s*ts**(-1)*ms**2 + 8*
     +    s*ts**(-1)*mi**2 + 2*s*ts**(-1)*mg**2 + s - 2*s**2*ts**(-1)
     +     + 2*tg*ug**(-1)*ui + 4*tg*ug**(-1)*ms**2 + 2*tg*ug**(-1)*
     +    mi**2 + 2*tg*ug**(-1)*mg**2 + 2*tg*ui*ts**(-1) + 6*tg*
     +    ts**(-1)*mi**2 + 2*tg*ts**(-1)*mg**2 + 2*ug**(-1)*ui*ms**2 - 
     +    2*ug**(-1)*ui*mi**2 - 2*ug**(-1)*ms**2*mi**2 + 2*ug**(-1)*
     +    ms**2*mg**2 + 4*ug**(-1)*ms**4 - 2*ug**(-1)*mi**2*mg**2 - 2*
     +    ug**(-1)*mi**4 + 2*ui*ts**(-1)*ms**2 - 2*ui*ts**(-1)*mi**2 + 
     +    6*ts**(-1)*ms**2*mi**2 + 2*ts**(-1)*ms**2*mg**2 - 2*ts**(-1)*
     +    mi**2*mg**2 - 6*ts**(-1)*mi**4 + 2*ms**2 - 2*mi**2 )


c               the prefactor for the scaling functions 
      NG_GBOS2 = real(gqbHs3s4) * (abs(mi)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_QGH(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc1,struc2,struc3
     &          ,eight,QF,s,tg,s4,mi,mg,ms,ti,ts,ui,us,ug,tshat
     &          ,sptg,s4i,hardfac,s4s,s4s_d,gams
     &          ,qgH,NGANG(0:9,0:9,-2:2,-2:2)
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)

c               the factorization scale 
      QF = massin(13) 

c               real kinematics built in
      ti = tg + mg**2 - mi**2
      ts = tg + mg**2 - ms**2 

      ui = s4 - s - tg 
      us = ui + mi**2 - ms**2 
      ug = ui + mi**2 - mg**2 

      s4i     = s4 
      if (s4i.lt.0.D0) s4i = 0.D0

      s4s     = s4 + mi**2 - ms**2 

      hardfac = 1.D0 + mi**2/s4i

c               some denominators following wim's notes 
      sptg   = s + tg 

c               following wim's notes  
      tshat = mg**2 - ms**2 - ui*tg/(s+tg) 

c               wim's coupling structures 
      struc1 = real( Cl(1)*Clc(1) + Cr(1)*Crc(1) )
      struc2 = real( Cl(3)*Clc(3) + Cr(3)*Crc(3) )
      struc3 = mi*mg*real(  Cl(1)*Cl(3) + Clc(1)*Clc(3) 
     &                     + Cr(1)*Cr(3) + Crc(1)*Crc(3) ) 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 ra-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the s4s denominator [hauptwert]
      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         s4s_d  = s4s/( s4s**2 + ms**2*gams**2 )
      else 
         s4s_d = 1.D0/s4s
      end if

c               the angular functions 
      call ANGULAR_ARRAY_NG_C(massin,NGANG)

c               wim's form output
c               color factors are Cf/Nc^2, Cf
c               replace s4s**(-1) -> s4s_d 
      qgH = tshat**(-2)*e**2*gs**4*struc1*CA**(-2)*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * (  - 8*tg*ms**2 + 8*tg*mi**2 )
      qgH = qgH + tshat**(-2)*e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-4)*eight**(-1) * ( 32*tg**2*ui**3 )
      qgH = qgH + tshat**(-2)*e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-3)*eight**(-1) * ( 32*tg*ui**2*mi**2 - 32*tg*ui**2*mg**2
     +     + 32*tg**2*ui**2 )
      qgH = qgH + tshat**(-2)*e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-2)*eight**(-1) * ( 32*tg*ui*mi**2 - 32*tg*ui*mg**2 + 16*
     +    tg**2*ui )
      qgH = qgH + tshat**(-2)*e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-1)*eight**(-1) * ( 16*tg*mi**2 - 16*tg*mg**2 )
      qgH = qgH + tshat**(-2)*e**2*gs**4*struc1*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * ( 8*tg*ms**2 - 8*tg*mi**2 )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc1*CA**(-2)*pi*Sn**2*
     + sptg**(-1)*eight**(-1) * (  - 8*tg )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc1*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * ( 8*tg )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc3*CA**(-2)*pi*Sn**2*
     + sptg**(-1)*eight**(-1) * ( 8*s*us**(-1) )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-3)*eight**(-1) * ( 32*s*ui**2*us**(-1) )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-2)*eight**(-1) * ( 32*s*ui*us**(-1) )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + sptg**(-1)*eight**(-1) * ( 16*s*us**(-1) )
      qgH = qgH + tshat**(-1)*e**2*gs**4*struc3*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * (  - 8*s*us**(-1) )
      qgH = qgH + e**2*gs**4*struc2*CA**(-2)*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * ( 8*s*ug*us**(-2) + 8*tg*ug*us**(-2) )
      qgH = qgH + e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*sptg**(-2)*
     + eight**(-1) * (  - 32*ug*ui*us**(-2)*mi**2 + 32*ug*ui*us**(-2)*
     +    mg**2 + 32*ug**2*ui*us**(-2) )
      qgH = qgH + e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*sptg**(-1)*
     + eight**(-1) * (  - 32*ug*us**(-2)*mi**2 + 32*ug*us**(-2)*mg**2
     +     + 32*ug**2*us**(-2) )
      qgH = qgH + e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*eight**(-1)
     +  * (  - 16*ug*ui**(-1)*us**(-2)*mi**2 + 16*ug*ui**(-1)*us**(-2)*
     +    mg**2 + 16*ug**2*ui**(-1)*us**(-2) )
      qgH = qgH + e**2*gs**4*struc2*pi*Sn**2*sptg**(-1)*eight**(-1)
     +  * (  - 8*s*ug*us**(-2) - 8*tg*ug*us**(-2) )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc2*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 4 + 4*s**(-1)*tg*ug*
     +    us**(-1) + 4*s**(-1)*ug*ui*us**(-1) + 4*s**(-1)*ug*s4i**(-1)*
     +    ms**(-2)*mi**4 - 4*s**(-1)*ug*s4i**(-1)*mi**2 - 4*s**(-1)*ug
     +     + 4*s**(-1)*ug**2*us**(-1) + 8*s**(-1)*s4i**(-1)*us**(-1)*
     +    ms**(-2)*mi**4*mg**4 + 16*s**(-1)*s4i**(-1)*us**(-1)*ms**2*
     +    mi**2*mg**2 + 8*s**(-1)*s4i**(-1)*us**(-1)*ms**2*mi**4 - 8*
     +    s**(-1)*s4i**(-1)*us**(-1)*ms**4*mi**2 - 8*s**(-1)*s4i**(-1)*
     +    us**(-1)*mi**2*mg**4 - 16*s**(-1)*s4i**(-1)*us**(-1)*mi**4*
     +    mg**2 - 4*s**(-1)*s4i**(-1)*ms**(-2)*mi**4*mg**2 - 4*s**(-1)*
     +    s4i**(-1)*ms**2*mi**2 + 4*s**(-1)*s4i**(-1)*mi**2*mg**2 + 4*
     +    s**(-1)*s4i**(-1)*mi**4 + 8*s**(-1)*us**(-1)*ms**2*mi**2 + 24
     +    *s**(-1)*us**(-1)*ms**2*mg**2 - 16*s**(-1)*us**(-1)*ms**4 - 8
     +    *s**(-1)*us**(-1)*mi**2*mg**2 - 8*s**(-1)*us**(-1)*mg**4 - 12
     +    *s**(-1)*ms**2 + 8*s**(-1)*mi**2 + 4*s**(-1)*mg**2 + 4*s*ug*
     +    us**(-2) )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc2*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * ( 4*s*us**(-1) + 4*tg*ug*
     +    us**(-2) + 4*tg*us**(-1) + 2*ug*ui*us**(-2) + 6*ug*us**(-2)*
     +    mi**2 + 2*ug*us**(-2)*mg**2 + 8*ug*us**(-1) + 2*ug**2*
     +    us**(-2) + 4*ui*us**(-1) + 8*s4i**(-1)*us**(-2)*ms**2*mi**2*
     +    mg**2 + 8*s4i**(-1)*us**(-2)*ms**2*mi**4 - 8*s4i**(-1)*
     +    us**(-2)*ms**4*mi**2 - 8*s4i**(-1)*us**(-2)*mi**4*mg**2 - 8*
     +    s4i**(-1)*us**(-1)*ms**(-2)*mi**4*mg**2 - 8*s4i**(-1)*
     +    us**(-1)*ms**2*mi**2 + 8*s4i**(-1)*us**(-1)*mi**2*mg**2 + 8*
     +    s4i**(-1)*us**(-1)*mi**4 + 4*s4i**(-1)*ms**(-2)*mi**4 - 4*
     +    s4i**(-1)*mi**2 + 4*us**(-2)*ms**2*mi**2 + 12*us**(-2)*ms**2*
     +    mg**2 - 12*us**(-2)*ms**4 - 4*us**(-2)*mi**2*mg**2 - 12*
     +    us**(-1)*ms**2 + 4*us**(-1)*mi**2 + 8*us**(-1)*mg**2 )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4 + 8*s**(-1)*tg**(-1)*
     +    s4i**(-1)*ms**(-2)*mi**4*mg**4 + 16*s**(-1)*tg**(-1)*
     +    s4i**(-1)*ms**2*mi**2*mg**2 + 8*s**(-1)*tg**(-1)*s4i**(-1)*
     +    ms**2*mi**4 - 8*s**(-1)*tg**(-1)*s4i**(-1)*ms**4*mi**2 - 8*
     +    s**(-1)*tg**(-1)*s4i**(-1)*mi**2*mg**4 - 16*s**(-1)*tg**(-1)*
     +    s4i**(-1)*mi**4*mg**2 + 16*s**(-1)*tg**(-1)*ms**2*mi**2 + 32*
     +    s**(-1)*tg**(-1)*ms**2*mg**2 - 24*s**(-1)*tg**(-1)*ms**4 - 16
     +    *s**(-1)*tg**(-1)*mi**2*mg**2 - 8*s**(-1)*tg**(-1)*mg**4 - 4*
     +    s**(-1)*ug*s4i**(-1)*ms**(-2)*mi**4 + 4*s**(-1)*ug*s4i**(-1)*
     +    mi**2 + 4*s**(-1)*ug + 4*s**(-1)*s4i**(-1)*ms**(-2)*mi**4*
     +    mg**2 + 4*s**(-1)*s4i**(-1)*ms**2*mi**2 - 4*s**(-1)*s4i**(-1)
     +    *mi**2*mg**2 - 4*s**(-1)*s4i**(-1)*mi**4 + 12*s**(-1)*ms**2
     +     - 8*s**(-1)*mi**2 - 4*s**(-1)*mg**2 + 4*s*tg**(-1)*ui*
     +    us**(-1) + 8*s*tg**(-1)*s4i**(-1)*ms**(-2)*mi**4 - 8*s*
     +    tg**(-1)*s4i**(-1)*mi**2 )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s*tg**(-1) - 4*s*ug*us**(-2)
     +     + 4*s*us**(-1) + 4*s**2*tg**(-1)*us**(-1) + 16*tg**(-2)*
     +    s4i**(-1)*ms**(-2)*mi**4*mg**4 + 16*tg**(-2)*s4i**(-1)*ms**2*
     +    mi**2*mg**2 - 16*tg**(-2)*s4i**(-1)*mi**2*mg**4 - 16*tg**(-2)
     +    *s4i**(-1)*mi**4*mg**2 + 32*tg**(-2)*ms**2*mg**2 - 16*
     +    tg**(-2)*mi**2*mg**2 - 16*tg**(-2)*mg**4 - 4*tg**(-1)*ug*
     +    us**(-1)*mi**2 - 4*tg**(-1)*ug*us**(-1)*mg**2 + 8*tg**(-1)*
     +    s4i**(-1)*us**(-1)*ms**(-2)*mi**4*mg**4 - 8*tg**(-1)*
     +    s4i**(-1)*us**(-1)*ms**2*mi**4 + 8*tg**(-1)*s4i**(-1)*
     +    us**(-1)*ms**4*mi**2 - 8*tg**(-1)*s4i**(-1)*us**(-1)*mi**2*
     +    mg**4 + 16*tg**(-1)*s4i**(-1)*ms**(-2)*mi**4*mg**2 + 16*
     +    tg**(-1)*s4i**(-1)*ms**2*mi**2 - 16*tg**(-1)*s4i**(-1)*mi**2*
     +    mg**2 - 16*tg**(-1)*s4i**(-1)*mi**4 - 8*tg**(-1)*us**(-1)*
     +    ms**2*mi**2 + 16*tg**(-1)*us**(-1)*ms**4 - 8*tg**(-1)*
     +    us**(-1)*mg**4 )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 32*tg**(-1)*ms**2 - 16*tg**(-1)*
     +    mi**2 - 16*tg**(-1)*mg**2 - 4*tg*ug*us**(-2) - 2*ug*ui*
     +    us**(-2) - 6*ug*us**(-2)*mi**2 - 2*ug*us**(-2)*mg**2 - 4*ug*
     +    us**(-1) - 2*ug**2*us**(-2) - 8*s4i**(-1)*us**(-2)*ms**2*
     +    mi**2*mg**2 - 8*s4i**(-1)*us**(-2)*ms**2*mi**4 + 8*s4i**(-1)*
     +    us**(-2)*ms**4*mi**2 + 8*s4i**(-1)*us**(-2)*mi**4*mg**2 + 8*
     +    s4i**(-1)*us**(-1)*ms**(-2)*mi**4*mg**2 + 8*s4i**(-1)*
     +    us**(-1)*ms**2*mi**2 - 8*s4i**(-1)*us**(-1)*mi**2*mg**2 - 8*
     +    s4i**(-1)*us**(-1)*mi**4 + 4*s4i**(-1)*ms**(-2)*mi**4 - 4*
     +    s4i**(-1)*mi**2 - 4*us**(-2)*ms**2*mi**2 - 12*us**(-2)*ms**2*
     +    mg**2 + 12*us**(-2)*ms**4 + 4*us**(-2)*mi**2*mg**2 + 12*
     +    us**(-1)*ms**2 - 4*us**(-1)*mi**2 - 8*us**(-1)*mg**2 )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc3*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 2*s**(-1)*ug*us**(-1) + 4
     +    *s**(-1) - 2*us**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 6*s**(-1)*tg**(-1)*ug - 2*
     +    s**(-1) + 2*s*tg**(-1)*us**(-1) + 16*tg**(-2)*mg**2 + 2*
     +    tg**(-1)*ug*us**(-1) + 2*tg**(-1)*ui*us**(-1) + 2*tg**(-1)*
     +    us**(-1)*mi**2 + 6*tg**(-1)*us**(-1)*mg**2 + 8*tg**(-1) + 6*
     +    us**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*s4s*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1) + 4*tg**(-2)*ti + 4*
     +    tg**(-2)*mi**2 - 4*tg**(-2)*mg**2 - 4*tg**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)*
     + eight**(-1) * (  - 4*s**(-1)*tg**(-1)*ug + 16*tg**(-2)*mg**2 + 8
     +    *tg**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*s4i**(-1)*ms**(-2)*
     +    mi**4 + 4*s**(-1)*s4i**(-1)*mi**2 - 8*s**(-1)*us**(-1)*ms**2
     +     + 8*s**(-1)*us**(-1)*mg**2 - 4*s**(-1) - 4*ug*us**(-2) - 4*
     +    us**(-2)*ms**2 + 4*us**(-2)*mg**2 - 4*us**(-1) )
      qgH = qgH + NGANG(0,0,0,0)*e**2*gs**4*struc2*Sn**2*hardfac**(-1)*
     + eight**(-1) * (  - 24*s**(-1)*tg**(-1)*ms**2 + 8*s**(-1)*
     +    tg**(-1)*mi**2 + 16*s**(-1)*tg**(-1)*mg**2 + 4*s**(-1)*
     +    s4i**(-1)*ms**(-2)*mi**4 - 4*s**(-1)*s4i**(-1)*mi**2 + 4*
     +    s**(-1) + 16*tg**(-2)*mg**2 + 8*tg**(-1)*us**(-1)*ms**2 + 16*
     +    tg**(-1) + 4*ug*us**(-2) + 4*us**(-2)*ms**2 - 4*us**(-2)*
     +    mg**2 + 4*us**(-1) )
      qgH = qgH + NGANG(1,0,-1,0)*s4s_d*e**2*gs**4*struc2*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ug*ui*
     +    us**(-1) - 4*s**(-1)*tg*ug**2*us**(-1) - 4*s**(-1)*tg**2*ug*
     +    us**(-1) - 8*s**(-1)*ug**2*ui*us**(-1) + 2*s*tg*ug*us**(-2)
     +     - 4*s*tg*us**(-1) - 2*s*ug*ti*us**(-2) - 6*s*ug*ui*us**(-2)
     +     - 2*s*ug**2*us**(-2) + 4*s*ti*us**(-1) - 2*tg*ug*ti*us**(-2)
     +     - 2*tg*ug*ui*us**(-2) - 4*tg*ug*us**(-1) - 2*tg*ug**2*
     +    us**(-2) + 4*tg*ti*us**(-1) - 4*tg*ui*us**(-1) + 2*tg**2*ug*
     +    us**(-2) - 4*tg**2*us**(-1) - 4*ug*ti*ui*us**(-2) - 4*ug*ui*
     +    us**(-1) - 4*ug*ui**2*us**(-2) - 4*ug**2*ui*us**(-2) - 4*
     +    ug**2*us**(-1) + 4*ti*ui*us**(-1) )
      qgH = qgH + NGANG(1,0,-1,0)*s4s_d*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s*tg**(-1)*ug*ti*us**(-1) + 4*s*
     +    tg**(-1)*ti*ui*us**(-1) - 2*s*tg*ug*us**(-2) - 4*s*tg*
     +    us**(-1) + 2*s*ug*ti*us**(-2) + 6*s*ug*ui*us**(-2) - 8*s*ug*
     +    us**(-1) + 2*s*ug**2*us**(-2) + 4*s*ti*us**(-1) - 4*s*ui*
     +    us**(-1) + 4*s**2*tg**(-1)*ti*us**(-1) - 4*s**2*us**(-1) + 8*
     +    tg**(-1)*ug*ti*ui*us**(-1) + 2*tg*ug*ti*us**(-2) + 2*tg*ug*ui
     +    *us**(-2) - 8*tg*ug*us**(-1) + 2*tg*ug**2*us**(-2) - 2*tg**2*
     +    ug*us**(-2) + 4*ug*ti*ui*us**(-2) + 4*ug*ti*us**(-1) - 8*ug*
     +    ui*us**(-1) + 4*ug*ui**2*us**(-2) + 4*ug**2*ui*us**(-2) )
      qgH = qgH + NGANG(1,0,-1,0)*s4s_d*e**2*gs**4*struc3*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 4 - 8*s**(-1)*ug - 4*
     +    s**(-1)*ms**2 + 4*s**(-1)*mi**2 - 2*s*us**(-1) - 2*ug*
     +    us**(-1) - 2*ti*us**(-1) - 2*ui*us**(-1) )
      qgH = qgH + NGANG(1,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2 - 2*s*tg**(-1) + 2*s*us**(-1)
     +     - 4*tg**(-1)*ug + 4*tg**(-1)*ti - 4*tg**(-1)*ms**2 + 4*
     +    tg**(-1)*mi**2 + 2*ug*us**(-1) + 2*ti*us**(-1) + 2*ui*
     +    us**(-1) )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4 - 8*s**(-1)*tg**(-1)*ug*ti - 8*
     +    s**(-1)*tg**(-1)*ug*mi**2 + 8*s**(-1)*tg**(-1)*ug*mg**2 - 16*
     +    s**(-1)*tg**(-1)*ti*ms**2 + 8*s**(-1)*tg**(-1)*ti*mi**2 + 8*
     +    s**(-1)*tg**(-1)*ti*mg**2 - 16*s**(-1)*tg**(-1)*ms**2*mi**2
     +     + 16*s**(-1)*tg**(-1)*ms**2*mg**2 + 8*s**(-1)*tg**(-1)*mi**4
     +     - 8*s**(-1)*tg**(-1)*mg**4 - 4*tg**(-2)*ug*ti - 4*tg**(-2)*
     +    ug*mi**2 + 4*tg**(-2)*ug*mg**2 + 4*tg**(-2)*ti*ui - 12*
     +    tg**(-2)*ti*ms**2 + 8*tg**(-2)*ti*mi**2 + 4*tg**(-2)*ti*mg**2
     +     + 4*tg**(-2)*ti**2 + 4*tg**(-2)*ui*mi**2 - 4*tg**(-2)*ui*
     +    mg**2 - 4*tg**(-2)*ts*ms**2 + 4*tg**(-2)*ts*mi**2 - 8*
     +    tg**(-2)*ms**2*mi**2 + 16*tg**(-2)*ms**2*mg**2 - 4*tg**(-2)*
     +    ms**4 + 4*tg**(-2)*mi**4 - 8*tg**(-2)*mg**4 + 4*tg**(-1)*ug
     +     - 12*tg**(-1)*ti - 4*tg**(-1)*ui - 4*tg**(-1)*ts + 12*
     +    tg**(-1)*ms**2 - 16*tg**(-1)*mi**2 + 4*tg**(-1)*mg**2 )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 8 + 8*tg**(-1)*ti + 8*tg**(-1)*ts - 8*
     +    tg**(-1)*ms**2 + 8*tg**(-1)*mi**2 )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s*ug*us**(-2) - 4*tg*ug*
     +    us**(-2) )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 4*s*ug*us**(-2) + 4*tg*ug*us**(-2) )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ug*us**(-1) - 4*
     +    s**(-1)*ti*us**(-1) - 4*s**(-1)*us**(-1)*ms**2 + 4*s**(-1)*
     +    us**(-1)*mg**2 - 4*tg**(-1)*ti*us**(-1) - 4*tg**(-1)*us**(-1)
     +    *mi**2 + 4*tg**(-1)*us**(-1)*mg**2 - 4*us**(-1) )
      qgH = qgH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 2*s*tg**(-1)*us**(-1) - 4*tg**(-1)*ug*
     +    us**(-1) + 4*tg**(-1)*ti*us**(-1) - 4*tg**(-1)*us**(-1)*ms**2
     +     + 4*tg**(-1)*us**(-1)*mi**2 - 2*us**(-1) )
      qgH = qgH + NGANG(1,4,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg**(-1)*ug*ti*ms**2
     +     + 8*s**(-1)*tg**(-1)*ug*ti*mg**2 + 16*s**(-1)*tg**(-1)*ug*
     +    ms**2*mg**2 - 8*s**(-1)*tg**(-1)*ug*ms**4 - 8*s**(-1)*
     +    tg**(-1)*ug*mg**4 + 32*s**(-1)*tg**(-1)*ti*ms**2*mg**2 - 16*
     +    s**(-1)*tg**(-1)*ti*ms**4 - 16*s**(-1)*tg**(-1)*ti*mg**4 - 8*
     +    s**(-1)*tg**(-1)*ti**2*ms**2 + 8*s**(-1)*tg**(-1)*ti**2*mg**2
     +     - 24*s**(-1)*tg**(-1)*ms**2*mg**4 + 24*s**(-1)*tg**(-1)*
     +    ms**4*mg**2 - 8*s**(-1)*tg**(-1)*ms**6 + 8*s**(-1)*tg**(-1)*
     +    mg**6 - 4*s**(-1)*ug*ti - 4*s**(-1)*ti*ms**2 + 4*s**(-1)*ti*
     +    mg**2 - 4*s**(-1)*ti**2 - 4*tg**(-2)*ug*ti*ms**2 + 4*tg**(-2)
     +    *ug*ti*mg**2 + 8*tg**(-2)*ug*ms**2*mg**2 - 4*tg**(-2)*ug*
     +    ms**4 - 4*tg**(-2)*ug*mg**4 + 4*tg**(-2)*ti*ui*ms**2 - 4*
     +    tg**(-2)*ti*ui*mg**2 + 24*tg**(-2)*ti*ms**2*mg**2 - 12*
     +    tg**(-2)*ti*ms**4 - 12*tg**(-2)*ti*mg**4 - 4*tg**(-2)*ti**2*
     +    ms**2 )
      qgH = qgH + NGANG(1,4,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*tg**(-2)*ti**2*mg**2 - 8*
     +    tg**(-2)*ui*ms**2*mg**2 + 4*tg**(-2)*ui*ms**4 + 4*tg**(-2)*ui
     +    *mg**4 - 24*tg**(-2)*ms**2*mg**4 + 24*tg**(-2)*ms**4*mg**2 - 
     +    8*tg**(-2)*ms**6 + 8*tg**(-2)*mg**6 - 2*tg**(-1)*ug*ti + 2*
     +    tg**(-1)*ug*ms**2 - 2*tg**(-1)*ug*mg**2 + 2*tg**(-1)*ti*ui - 
     +    6*tg**(-1)*ti*ms**2 + 6*tg**(-1)*ti*mg**2 - 2*tg**(-1)*ti**2
     +     - 2*tg**(-1)*ui*ms**2 + 2*tg**(-1)*ui*mg**2 + 8*tg**(-1)*
     +    ms**2*mg**2 - 4*tg**(-1)*ms**4 - 4*tg**(-1)*mg**4 - 2*ti - 2*
     +    ms**2 + 2*mg**2 )
      qgH = qgH + NGANG(1,4,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ug*ti*us**(-1) - 8*
     +    s**(-1)*ug*us**(-1)*ms**2 + 8*s**(-1)*ug*us**(-1)*mg**2 - 4*
     +    s**(-1)*ug**2*us**(-1) - 8*s**(-1)*ti*us**(-1)*ms**2 + 8*
     +    s**(-1)*ti*us**(-1)*mg**2 - 4*s**(-1)*ti**2*us**(-1) + 8*
     +    s**(-1)*us**(-1)*ms**2*mg**2 - 4*s**(-1)*us**(-1)*ms**4 - 4*
     +    s**(-1)*us**(-1)*mg**4 - s*tg**(-1)*ug*us**(-1) - s*tg**(-1)*
     +    ti*us**(-1) + s*tg**(-1)*ui*us**(-1) - 2*s*tg**(-1)*us**(-1)*
     +    ms**2 + 2*s*tg**(-1)*us**(-1)*mg**2 - s*us**(-1) - 4*tg**(-1)
     +    *ug*ti*us**(-1) + 2*tg**(-1)*ug*ui*us**(-1) - 6*tg**(-1)*ug*
     +    us**(-1)*ms**2 + 6*tg**(-1)*ug*us**(-1)*mg**2 - 2*tg**(-1)*
     +    ug**2*us**(-1) + 2*tg**(-1)*ti*ui*us**(-1) - 6*tg**(-1)*ti*
     +    us**(-1)*ms**2 + 6*tg**(-1)*ti*us**(-1)*mg**2 - 2*tg**(-1)*
     +    ti**2*us**(-1) + 2*tg**(-1)*ui*us**(-1)*ms**2 - 2*tg**(-1)*ui
     +    *us**(-1)*mg**2 + 8*tg**(-1)*us**(-1)*ms**2*mg**2 - 4*
     +    tg**(-1)*us**(-1)*ms**4 )
      qgH = qgH + NGANG(1,4,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*tg**(-1)*us**(-1)*mg**4 - 4*
     +    ug*us**(-1) - 4*ti*us**(-1) - 4*us**(-1)*ms**2 + 4*us**(-1)*
     +    mg**2 )
      qgH = qgH + NGANG(1,7,-1,1)*s4s_d*e**2*gs**4*struc2*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ug*us**(-1)
     +     - 4*s**(-1)*tg*ui*us**(-1) - 4*s**(-1)*tg**2*us**(-1) - 8*
     +    s**(-1)*ug*ui*us**(-1) - 4*s*us**(-1) - 8*tg*us**(-1) - 4*ug*
     +    us**(-1) - 4*ui*us**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*s4s_d*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s*tg**(-1)*ug*us**(-1) - 4*s*
     +    tg**(-1)*ui*us**(-1) - 8*s*us**(-1) - 4*s**2*tg**(-1)*
     +    us**(-1) - 8*tg**(-1)*ug*ui*us**(-1) - 4*tg*us**(-1) - 4*ug*
     +    us**(-1) - 4*ui*us**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*s4s_d*e**2*gs**4*struc3*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*tg**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg**(-1)*ti - 8*
     +    s**(-1)*tg**(-1)*mi**2 + 8*s**(-1)*tg**(-1)*mg**2 - 8*
     +    tg**(-2)*ti - 8*tg**(-2)*mi**2 + 8*tg**(-2)*mg**2 + 8*
     +    tg**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 8*tg**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*us**(-1) )
      qgH = qgH + NGANG(1,7,-1,1)*e**2*gs**4*struc3*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 4*tg**(-1)*us**(-1) )
      qgH = qgH + NGANG(1,9,-1,-2)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*tg**(-1)*ug*ts*ms**2 + 4*
     +    tg**(-1)*ug*ts*mi**2 + 4*tg**(-1)*ui*ts*ms**2 - 4*tg**(-1)*ui
     +    *ts*mi**2 + 8*tg**(-1)*ts*ms**2*mi**2 - 4*tg**(-1)*ts*ms**4
     +     - 4*tg**(-1)*ts*mi**4 + 4*tg**(-1)*ts**2*ms**2 - 4*tg**(-1)*
     +    ts**2*mi**2 + 2*tg*ms**2 - 2*tg*mi**2 + 2*ug*ms**2 - 2*ug*
     +    mi**2 - 2*ui*ms**2 + 2*ui*mi**2 - 2*ts*ms**2 + 2*ts*mi**2 - 4
     +    *ms**2*mi**2 + 2*ms**4 + 2*mi**4 )
      qgH = qgH + NGANG(1,9,-1,-2)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*tg**(-1)*ts**2*ms**2 + 8*
     +    tg**(-1)*ts**2*mi**2 - 4*tg*ms**2 + 4*tg*mi**2 + 8*ts*ms**2
     +     - 8*ts*mi**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*s4s_d*e**2*gs**4*struc3*CA**(-2)
     + *Sn**2*hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ug*ms**2 + 8*
     +    s**(-1)*ug*mi**2 - 4*s**(-1)*ug**2 + 8*s**(-1)*ms**2*mi**2 - 
     +    4*s**(-1)*ms**4 - 4*s**(-1)*mi**4 + s*tg*us**(-1) - 3*s*ug*
     +    us**(-1) - s*ui*us**(-1) - s*ts*us**(-1) - 3*s*us**(-1)*ms**2
     +     + 3*s*us**(-1)*mi**2 - 2*ug*ui*us**(-1) - 2*ug*ts*us**(-1)
     +     - 4*ug*us**(-1)*ms**2 + 4*ug*us**(-1)*mi**2 - 4*ug - 2*ug**2
     +    *us**(-1) - 2*ui*us**(-1)*ms**2 + 2*ui*us**(-1)*mi**2 - 2*ts*
     +    us**(-1)*ms**2 + 2*ts*us**(-1)*mi**2 + 4*us**(-1)*ms**2*mi**2
     +     - 2*us**(-1)*ms**4 - 2*us**(-1)*mi**4 - 4*ms**2 + 4*mi**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*s*tg**(-1)*ts - s*tg*us**(-1) + 
     +    3*s*ug*us**(-1) + s*ui*us**(-1) + s*ts*us**(-1) + 3*s*
     +    us**(-1)*ms**2 - 3*s*us**(-1)*mi**2 - 2*s + 4*tg**(-1)*ug*ts
     +     + 4*tg**(-1)*ts*ms**2 - 4*tg**(-1)*ts*mi**2 + 2*ug*ui*
     +    us**(-1) + 2*ug*ts*us**(-1) + 4*ug*us**(-1)*ms**2 - 4*ug*
     +    us**(-1)*mi**2 - 2*ug + 2*ug**2*us**(-1) + 2*ui*us**(-1)*
     +    ms**2 - 2*ui*us**(-1)*mi**2 + 2*ts*us**(-1)*ms**2 - 2*ts*
     +    us**(-1)*mi**2 - 4*us**(-1)*ms**2*mi**2 + 2*us**(-1)*ms**4 + 
     +    2*us**(-1)*mi**4 - 2*ms**2 + 2*mi**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg**(-1)*ug*ts*ms**2
     +     + 8*s**(-1)*tg**(-1)*ug*ts*mi**2 + 16*s**(-1)*tg**(-1)*ts*
     +    ms**2*mi**2 - 8*s**(-1)*tg**(-1)*ts*ms**4 - 8*s**(-1)*
     +    tg**(-1)*ts*mi**4 + 4*s**(-1)*ug*ts + 4*s**(-1)*ug*ms**2 - 4*
     +    s**(-1)*ug*mi**2 + 4*s**(-1)*ts*ms**2 - 4*s**(-1)*ts*mi**2 - 
     +    8*s**(-1)*ms**2*mi**2 + 4*s**(-1)*ms**4 + 4*s**(-1)*mi**4 - 4
     +    *tg**(-2)*ug*ts*ms**2 + 4*tg**(-2)*ug*ts*mi**2 + 4*tg**(-2)*
     +    ui*ts*ms**2 - 4*tg**(-2)*ui*ts*mi**2 + 8*tg**(-2)*ts*ms**2*
     +    mi**2 - 4*tg**(-2)*ts*ms**4 - 4*tg**(-2)*ts*mi**4 + 4*
     +    tg**(-2)*ts**2*ms**2 - 4*tg**(-2)*ts**2*mi**2 + 6*tg**(-1)*ug
     +    *ms**2 - 6*tg**(-1)*ug*mi**2 - 6*tg**(-1)*ui*ms**2 + 6*
     +    tg**(-1)*ui*mi**2 - 18*tg**(-1)*ts*ms**2 + 18*tg**(-1)*ts*
     +    mi**2 + 4*tg**(-1)*ts**2 - 12*tg**(-1)*ms**2*mi**2 + 6*
     +    tg**(-1)*ms**4 + 6*tg**(-1)*mi**4 + 2*tg + 2*ug - 2*ui - 2*ts
     +     + 8*ms**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*mi**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 16*tg**(-1)*ts*ms**2 - 16*tg**(-1)
     +    *ts*mi**2 - 8*tg**(-1)*ts**2 - 4*tg + 8*ts - 8*ms**2 + 8*
     +    mi**2 )
      qgH = qgH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( s*tg**(-1)*ug*us**(-1) - s*
     +    tg**(-1)*ui*us**(-1) - s*tg**(-1)*ts*us**(-1) + s*tg**(-1)*
     +    us**(-1)*ms**2 - s*tg**(-1)*us**(-1)*mi**2 - 3*s*us**(-1) - 2
     +    *tg**(-1)*ug*ui*us**(-1) - 2*tg**(-1)*ug*ts*us**(-1) + 4*
     +    tg**(-1)*ug*us**(-1)*ms**2 - 4*tg**(-1)*ug*us**(-1)*mi**2 + 2
     +    *tg**(-1)*ug**2*us**(-1) - 2*tg**(-1)*ui*us**(-1)*ms**2 + 2*
     +    tg**(-1)*ui*us**(-1)*mi**2 - 2*tg**(-1)*ts*us**(-1)*ms**2 + 2
     +    *tg**(-1)*ts*us**(-1)*mi**2 - 4*tg**(-1)*us**(-1)*ms**2*mi**2
     +     + 2*tg**(-1)*us**(-1)*ms**4 + 2*tg**(-1)*us**(-1)*mi**4 )
      qgH = qgH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*s*tg**(-1)*ts*us**(-1) + 2*s*
     +    us**(-1) + 4*tg**(-1)*ug*ts*us**(-1) + 4*tg**(-1)*ts*us**(-1)
     +    *ms**2 - 4*tg**(-1)*ts*us**(-1)*mi**2 - 2*ug*us**(-1) - 2*
     +    us**(-1)*ms**2 + 2*us**(-1)*mi**2 )
      qgH = qgH + NGANG(3,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*ms**(-2)*mi**2*mg**2
     +     - 4*s**(-1)*ms**(-2)*mg**4 + 8*s**(-1)*mg**2 + 4*ms**(-2)*
     +    mg**2 )
      qgH = qgH + NGANG(3,0,-1,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 4*s**(-1)*ms**(-2)*mi**2*mg**2 + 4*s**(-1)*
     +    ms**(-2)*mg**4 - 8*s**(-1)*mg**2 - 4*ms**(-2)*mg**2 )
      qgH = qgH + NGANG(3,7,-1,1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*ms**(-2)*mg**2 )
      qgH = qgH + NGANG(3,7,-1,1)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 4*s**(-1)*ms**(-2)*mg**2 )
      qgH = qgH + NGANG(3,9,-1,-2)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*ms**2*mg**2 - 8*mi**2*mg**2 )
      qgH = qgH + NGANG(3,9,-1,-2)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*ms**2*mg**2 + 8*mi**2*mg**2 )
      qgH = qgH + NGANG(3,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s**(-1)*ms**(-2)*mi**4*mg**2 + 8
     +    *s**(-1)*ms**2*mg**2 - 16*s**(-1)*mi**2*mg**2 - 8*ms**(-2)*
     +    mi**2*mg**2 + 8*mg**2 )
      qgH = qgH + NGANG(3,9,-1,-1)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ms**(-2)*mi**4*mg**2
     +     - 8*s**(-1)*ms**2*mg**2 + 16*s**(-1)*mi**2*mg**2 + 8*
     +    ms**(-2)*mi**2*mg**2 - 8*mg**2 )
      qgH = qgH + NGANG(4,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2 - 8*s**(-1)*tg**(-1)*ti*ms**2
     +     + 8*s**(-1)*tg**(-1)*ti*mg**2 - 8*s**(-1)*tg**(-1)*ms**2*
     +    mi**2 + 8*s**(-1)*tg**(-1)*ms**2*mg**2 + 8*s**(-1)*tg**(-1)*
     +    mi**2*mg**2 - 8*s**(-1)*tg**(-1)*mg**4 - 4*s**(-1)*ug - 4*
     +    s**(-1)*ti + 4*s**(-1)*ms**(-2)*mi**2*mg**2 - 4*s**(-1)*ms**2
     +     - 4*tg**(-2)*ti*ms**2 + 4*tg**(-2)*ti*mg**2 - 4*tg**(-2)*
     +    ms**2*mi**2 + 4*tg**(-2)*ms**2*mg**2 + 4*tg**(-2)*mi**2*mg**2
     +     - 4*tg**(-2)*mg**4 - 2*tg**(-1)*ti + 4*tg**(-1)*ms**2 - 2*
     +    tg**(-1)*mi**2 - 2*tg**(-1)*mg**2 - 4*ms**(-2)*mg**2 )
      qgH = qgH + NGANG(4,0,-1,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 4*s**(-1)*tg**(-1)*ug*ms**2 + 4*s**(-1)*
     +    tg**(-1)*ug*mg**2 - 4*s**(-1)*ms**(-2)*mi**2*mg**2 + 4*
     +    s**(-1)*ms**2 - 2*tg**(-1)*ug + 2*tg**(-1)*ui + 4*tg**(-1)*
     +    ms**2 + 2*tg**(-1)*mi**2 + 2*tg**(-1)*mg**2 + 4*ms**(-2)*
     +    mg**2 )
      qgH = qgH + NGANG(4,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*ug*us**(-1) - 4*
     +    s**(-1)*ti*us**(-1) - 4*s**(-1)*us**(-1)*ms**2 + 4*s**(-1)*
     +    us**(-1)*mg**2 - 2*tg**(-1)*ti*us**(-1) - 2*tg**(-1)*us**(-1)
     +    *mi**2 + 2*tg**(-1)*us**(-1)*mg**2 )
      qgH = qgH + NGANG(4,7,-1,1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*ms**(-2)*mg**2 - 4*
     +    s**(-1) )
      qgH = qgH + NGANG(4,7,-1,1)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 4*s**(-1)*tg**(-1)*ug + 4*s**(-1)*ms**(-2)*
     +    mg**2 )
      qgH = qgH + NGANG(4,9,-1,-2)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*tg**(-1)*ug*ms**2*mi**2 - 4*
     +    tg**(-1)*ug*ms**2*mg**2 + 4*tg**(-1)*ug*ms**4 + 4*tg**(-1)*ug
     +    *mi**2*mg**2 + 4*tg**(-1)*ui*ms**2*mi**2 + 4*tg**(-1)*ui*
     +    ms**2*mg**2 - 4*tg**(-1)*ui*ms**4 - 4*tg**(-1)*ui*mi**2*mg**2
     +     + 16*tg**(-1)*ms**2*mi**2*mg**2 + 4*tg**(-1)*ms**2*mi**4 + 4
     +    *tg**(-1)*ms**2*mg**4 - 12*tg**(-1)*ms**4*mi**2 - 12*tg**(-1)
     +    *ms**4*mg**2 + 8*tg**(-1)*ms**6 - 4*tg**(-1)*mi**2*mg**4 - 4*
     +    tg**(-1)*mi**4*mg**2 + 8*ms**2*mi**2 - 8*ms**4 )
      qgH = qgH + NGANG(4,9,-1,-2)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*tg**(-1)*ug*ms**2*mi**2 + 4*
     +    tg**(-1)*ug*ms**2*mg**2 - 4*tg**(-1)*ug*ms**4 - 4*tg**(-1)*ug
     +    *mi**2*mg**2 - 4*tg**(-1)*ui*ms**2*mi**2 - 4*tg**(-1)*ui*
     +    ms**2*mg**2 + 4*tg**(-1)*ui*ms**4 + 4*tg**(-1)*ui*mi**2*mg**2
     +     - 4*tg**(-1)*ms**2*mi**4 - 4*tg**(-1)*ms**2*mg**4 - 4*
     +    tg**(-1)*ms**4*mi**2 - 4*tg**(-1)*ms**4*mg**2 + 8*tg**(-1)*
     +    ms**6 + 4*tg**(-1)*mi**2*mg**4 + 4*tg**(-1)*mi**4*mg**2 - 8*
     +    ms**2*mi**2 + 8*ms**4 )
      qgH = qgH + NGANG(4,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg**(-1)*ug*ms**2*
     +    mi**2 - 8*s**(-1)*tg**(-1)*ug*ms**2*mg**2 + 8*s**(-1)*
     +    tg**(-1)*ug*ms**4 + 8*s**(-1)*tg**(-1)*ug*mi**2*mg**2 + 16*
     +    s**(-1)*tg**(-1)*ms**2*mi**2*mg**2 + 8*s**(-1)*tg**(-1)*ms**2
     +    *mi**4 - 16*s**(-1)*tg**(-1)*ms**4*mi**2 - 8*s**(-1)*tg**(-1)
     +    *ms**4*mg**2 + 8*s**(-1)*tg**(-1)*ms**6 - 8*s**(-1)*tg**(-1)*
     +    mi**4*mg**2 - 8*s**(-1)*ms**(-2)*mi**4*mg**2 + 12*s**(-1)*
     +    ms**2*mi**2 - 4*s**(-1)*ms**2*mg**2 - 8*s**(-1)*ms**4 + 12*
     +    s**(-1)*mi**2*mg**2 - 4*s**(-1)*mi**4 - 4*tg**(-2)*ug*ms**2*
     +    mi**2 - 4*tg**(-2)*ug*ms**2*mg**2 + 4*tg**(-2)*ug*ms**4 + 4*
     +    tg**(-2)*ug*mi**2*mg**2 + 4*tg**(-2)*ui*ms**2*mi**2 + 4*
     +    tg**(-2)*ui*ms**2*mg**2 - 4*tg**(-2)*ui*ms**4 - 4*tg**(-2)*ui
     +    *mi**2*mg**2 + 16*tg**(-2)*ms**2*mi**2*mg**2 + 4*tg**(-2)*
     +    ms**2*mi**4 + 4*tg**(-2)*ms**2*mg**4 - 12*tg**(-2)*ms**4*
     +    mi**2 )
      qgH = qgH + NGANG(4,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 12*tg**(-2)*ms**4*mg**2 + 8*
     +    tg**(-2)*ms**6 - 4*tg**(-2)*mi**2*mg**4 - 4*tg**(-2)*mi**4*
     +    mg**2 + 2*tg**(-1)*ug*ms**2 - 2*tg**(-1)*ug*mi**2 - 2*
     +    tg**(-1)*ui*ms**2 + 2*tg**(-1)*ui*mi**2 - 14*tg**(-1)*ms**2*
     +    mi**2 - 18*tg**(-1)*ms**2*mg**2 + 16*tg**(-1)*ms**4 + 10*
     +    tg**(-1)*mi**2*mg**2 + 2*tg**(-1)*mi**4 + 4*tg**(-1)*mg**4 + 
     +    8*ms**(-2)*mi**2*mg**2 + 4*mi**2 - 4*mg**2 )
      qgH = qgH + NGANG(4,9,-1,-1)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s**(-1)*tg**(-1)*ug*ms**2*mi**2
     +     + 8*s**(-1)*tg**(-1)*ug*ms**2*mg**2 - 8*s**(-1)*tg**(-1)*ug*
     +    ms**4 - 8*s**(-1)*tg**(-1)*ug*mi**2*mg**2 + 8*s**(-1)*
     +    ms**(-2)*mi**4*mg**2 - 12*s**(-1)*ms**2*mi**2 + 4*s**(-1)*
     +    ms**2*mg**2 + 8*s**(-1)*ms**4 - 12*s**(-1)*mi**2*mg**2 + 4*
     +    s**(-1)*mi**4 - 2*tg**(-1)*ug*ms**2 + 2*tg**(-1)*ug*mi**2 + 2
     +    *tg**(-1)*ui*ms**2 - 2*tg**(-1)*ui*mi**2 - 2*tg**(-1)*ms**2*
     +    mi**2 + 2*tg**(-1)*ms**2*mg**2 + 8*tg**(-1)*ms**4 - 2*
     +    tg**(-1)*mi**2*mg**2 - 2*tg**(-1)*mi**4 - 4*tg**(-1)*mg**4 - 
     +    8*ms**(-2)*mi**2*mg**2 - 4*mi**2 + 4*mg**2 )
      qgH = qgH + NGANG(4,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( s*tg**(-1)*ug*us**(-1) - s*
     +    tg**(-1)*ui*us**(-1) + 2*s*tg**(-1)*us**(-1)*ms**2 - s*
     +    tg**(-1)*us**(-1)*mi**2 - s*tg**(-1)*us**(-1)*mg**2 - 2*
     +    tg**(-1)*ug*ui*us**(-1) + 6*tg**(-1)*ug*us**(-1)*ms**2 - 4*
     +    tg**(-1)*ug*us**(-1)*mi**2 - 2*tg**(-1)*ug*us**(-1)*mg**2 + 2
     +    *tg**(-1)*ug**2*us**(-1) - 2*tg**(-1)*ui*us**(-1)*ms**2 + 2*
     +    tg**(-1)*ui*us**(-1)*mi**2 - 6*tg**(-1)*us**(-1)*ms**2*mi**2
     +     - 2*tg**(-1)*us**(-1)*ms**2*mg**2 + 4*tg**(-1)*us**(-1)*
     +    ms**4 + 2*tg**(-1)*us**(-1)*mi**2*mg**2 + 2*tg**(-1)*us**(-1)
     +    *mi**4 - 2*ug*us**(-1) )
      qgH = qgH + NGANG(6,9,1,-2)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 8*tg**(-1)*ms**2 - 8*tg**(-1)*mi**2 )
      qgH = qgH + NGANG(6,9,1,-1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2*s**(-1)*tg**(-1)*ug + 2*s*
     +    tg**(-1)*us**(-1) + 2*tg**(-1)*ug*us**(-1) - 2*tg**(-1) )
      qgH = qgH + NGANG(6,9,1,-1)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * (  - 4*s**(-1)*tg**(-1)*ug + 4*tg**(-1) )
      qgH = qgH + NGANG(9,0,-2,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*tg**(-1)*ts*ms**2 - 4*tg**(-1)*
     +    ts*mi**2 - 4*tg**(-1)*ms**2*mi**2 - 4*tg**(-1)*ms**2*mg**2 + 
     +    4*tg**(-1)*ms**4 + 4*tg**(-1)*mi**2*mg**2 - 4*ms**2 + 4*mi**2
     +     )
      qgH = qgH + NGANG(9,0,-2,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 8*s*tg**(-1)*ms**2 - 8*s*tg**(-1)*mi**2 - 16*
     +    tg**(-2)*ms**2*mi**2*mg**2 - 16*tg**(-2)*ms**2*mg**4 + 16*
     +    tg**(-2)*ms**4*mg**2 + 16*tg**(-2)*mi**2*mg**4 - 2*tg**(-1)*
     +    ug*ms**2 + 2*tg**(-1)*ug*mi**2 + 2*tg**(-1)*ui*ms**2 - 2*
     +    tg**(-1)*ui*mi**2 - 8*tg**(-1)*ts*ms**2 + 8*tg**(-1)*ts*mi**2
     +     - 6*tg**(-1)*ms**2*mi**2 - 10*tg**(-1)*ms**2*mg**2 + 8*
     +    tg**(-1)*ms**4 + 10*tg**(-1)*mi**2*mg**2 - 2*tg**(-1)*mi**4
     +     + 8*ms**2 - 8*mi**2 )
      qgH = qgH + NGANG(9,0,-1,0)*s4s_d*e**2*gs**4*struc3*CA**(-2)*
     + Sn**2*hardfac**(-1)*eight**(-1) * ( 2 + 2*s**(-1)*ug + 2*s*
     +    us**(-1) + 2*ug*us**(-1) )
      qgH = qgH + NGANG(9,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg**(-1)*ug*ms**2 + 4
     +    *s**(-1)*tg**(-1)*ug*mi**2 - 4*s**(-1)*tg**(-1)*ug**2 - 2*
     +    s**(-1)*ug + 2*s**(-1)*ms**2 - 2*s**(-1)*mi**2 + 16*s*
     +    tg**(-2)*mg**2 + 5*s*tg**(-1)*ug*us**(-1) + s*tg**(-1)*ui*
     +    us**(-1) + 2*s*tg**(-1)*us**(-1)*ms**2 + s*tg**(-1)*us**(-1)*
     +    mi**2 + 5*s*tg**(-1)*us**(-1)*mg**2 + 10*s*tg**(-1) + 4*s*
     +    us**(-1) + 2*s**2*tg**(-1)*us**(-1) + 16*tg**(-2)*ug*mg**2 + 
     +    16*tg**(-2)*ms**2*mg**2 - 16*tg**(-2)*mi**2*mg**2 + 2*
     +    tg**(-1)*ug*ui*us**(-1) + 2*tg**(-1)*ug*us**(-1)*ms**2 + 6*
     +    tg**(-1)*ug*us**(-1)*mg**2 + 6*tg**(-1)*ug + 2*tg**(-1)*ug**2
     +    *us**(-1) + 2*tg**(-1)*ui*us**(-1)*ms**2 - 2*tg**(-1)*ui*
     +    us**(-1)*mi**2 + 2*tg**(-1)*us**(-1)*ms**2*mi**2 + 6*tg**(-1)
     +    *us**(-1)*ms**2*mg**2 - 6*tg**(-1)*us**(-1)*mi**2*mg**2 - 2*
     +    tg**(-1)*us**(-1)*mi**4 + 14*tg**(-1)*ms**2 - 16*tg**(-1)*
     +    mi**2 )
      qgH = qgH + NGANG(9,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*tg**(-1)*mg**2 + 6*ug*us**(-1)
     +     + 6*us**(-1)*ms**2 - 6*us**(-1)*mi**2 )
      qgH = qgH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*ug + 4*tg**(-2)*ts*ms**2
     +     - 4*tg**(-2)*ts*mi**2 - 4*tg**(-2)*ms**2*mi**2 - 4*tg**(-2)*
     +    ms**2*mg**2 + 4*tg**(-2)*ms**4 + 4*tg**(-2)*mi**2*mg**2 + 4*
     +    tg**(-1)*ts + 4*tg**(-1)*mi**2 - 4*tg**(-1)*mg**2 )
      qgH = qgH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 4 - 8*s**(-1)*tg**(-1)*ug*ms**2 + 4*s**(-1)*
     +    tg**(-1)*ug*mi**2 + 4*s**(-1)*tg**(-1)*ug*mg**2 + 4*s**(-1)*
     +    ms**2 - 4*s**(-1)*mi**2 + 4*s*tg**(-1) + 32*tg**(-2)*ms**2*
     +    mg**2 - 16*tg**(-2)*mi**2*mg**2 - 16*tg**(-2)*mg**4 - 2*
     +    tg**(-1)*ug - 2*tg**(-1)*ui - 8*tg**(-1)*ts + 12*tg**(-1)*
     +    ms**2 - 14*tg**(-1)*mi**2 - 6*tg**(-1)*mg**2 )
      qgH = qgH + NGANG(9,0,-1,0)*e**2*gs**4*struc3*Sn**2*hardfac**(-1)
     + *eight**(-1) * ( 2*s*tg**(-1)*us**(-1) - 2*tg**(-1)*ug*us**(-1)
     +     )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF*pi*Sn**2*sptg**(-4)*eight**(-1) * ( 32*tg**2*
     +    ui**3 )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * ( 32*tg*
     +    ui**2*mi**2 - 32*tg*ui**2*mg**2 + 32*tg**2*ui**2 )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * ( 32*tg*ui*
     +    mi**2 - 32*tg*ui*mg**2 + 16*tg**2*ui )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * ( 16*tg*
     +    mi**2 - 16*tg*mg**2 )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * ( 32*s*
     +    ui**2*us**(-1) )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * ( 32*s*ui*
     +    us**(-1) )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * ( 16*s*
     +    us**(-1) )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc2*CA**(-1)*
     + CF*pi*Sn**2*sptg**(-2)*eight**(-1) * (  - 32*ug*ui*us**(-2)*
     +    mi**2 + 32*ug*ui*us**(-2)*mg**2 + 32*ug**2*ui*us**(-2) )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc2*CA**(-1)*
     + CF*pi*Sn**2*sptg**(-1)*eight**(-1) * (  - 32*ug*us**(-2)*mi**2
     +     + 32*ug*us**(-2)*mg**2 + 32*ug**2*us**(-2) )
      qgH = qgH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc2*CA**(-1)*
     + CF*pi*Sn**2*eight**(-1) * (  - 16*ug*ui**(-1)*us**(-2)*mi**2 + 
     +    16*ug*ui**(-1)*us**(-2)*mg**2 + 16*ug**2*ui**(-1)*us**(-2) )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-4)*eight**(-1) * (  - 32*tg**2*
     +    ui**3 )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * (  - 32*tg*ui**2*
     +    mi**2 + 32*tg*ui**2*mg**2 - 32*tg**2*ui**2 )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * (  - 32*tg*ui*
     +    mi**2 + 32*tg*ui*mg**2 - 16*tg**2*ui )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * (  - 16*tg*mi**2
     +     + 16*tg*mg**2 )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * (  - 32*s*ui**2*
     +    us**(-1) )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * (  - 32*s*ui*
     +    us**(-1) )
      qgH = qgH + log(s4i*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * (  - 16*s*us**(-1)
     +     )
      qgH = qgH + log(s4i*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi*
     + Sn**2*sptg**(-2)*eight**(-1) * ( 32*ug*ui*us**(-2)*mi**2 - 32*ug
     +    *ui*us**(-2)*mg**2 - 32*ug**2*ui*us**(-2) )
      qgH = qgH + log(s4i*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi*
     + Sn**2*sptg**(-1)*eight**(-1) * ( 32*ug*us**(-2)*mi**2 - 32*ug*
     +    us**(-2)*mg**2 - 32*ug**2*us**(-2) )
      qgH = qgH + log(s4i*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi*
     + Sn**2*eight**(-1) * ( 16*ug*ui**(-1)*us**(-2)*mi**2 - 16*ug*
     +    ui**(-1)*us**(-2)*mg**2 - 16*ug**2*ui**(-1)*us**(-2) )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-4)*eight**(-1) * ( 32*tg**2*ui**3 )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * ( 32*tg*ui**2*
     +    mi**2 - 32*tg*ui**2*mg**2 + 32*tg**2*ui**2 )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * ( 32*tg*ui*mi**2
     +     - 32*tg*ui*mg**2 + 16*tg**2*ui )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * ( 16*tg*mi**2 - 16
     +    *tg*mg**2 )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-3)*eight**(-1) * ( 32*s*ui**2*
     +    us**(-1) )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-2)*eight**(-1) * ( 32*s*ui*us**(-1)
     +     )
      qgH = qgH + log(QF**2*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*sptg**(-1)*eight**(-1) * ( 16*s*us**(-1) )
      qgH = qgH + log(QF**2*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi
     + *Sn**2*sptg**(-2)*eight**(-1) * (  - 32*ug*ui*us**(-2)*mi**2 + 
     +    32*ug*ui*us**(-2)*mg**2 + 32*ug**2*ui*us**(-2) )
      qgH = qgH + log(QF**2*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi
     + *Sn**2*sptg**(-1)*eight**(-1) * (  - 32*ug*us**(-2)*mi**2 + 32*
     +    ug*us**(-2)*mg**2 + 32*ug**2*us**(-2) )
      qgH = qgH + log(QF**2*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF*pi
     + *Sn**2*eight**(-1) * (  - 16*ug*ui**(-1)*us**(-2)*mi**2 + 16*ug*
     +    ui**(-1)*us**(-2)*mg**2 + 16*ug**2*ui**(-1)*us**(-2) )

ctp checking the integrals
c      qgH = Sn**2*hardfac**(-1) * NGANG(0,0,0,0)
c      qgH = 2.D0*pi

c               the prefactor for the scaling functions 
      NG_QGH = qgH * (abs(mi)+mg)**2/4.D0

      end


c --------------------------------------------------------------------
      real*8 function NG_GBH(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc1,struc2,struc3
     &          ,eight,QF,s,tg,s4,mi,mg,ms,ti,ts,ui,us,ug,ushat
     &          ,spug,s4i,hardfac,s4s,s4s_d,gams
     &          ,gqbH,NGANG(0:9,0:9,-2:2,-2:2)
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      eight = 8.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s    = massin(1)
      tg   = massin(2)
      s4   = massin(3)
      mi   = massin(6)
      mg   = massin(7)
      ms   = massin(11)
      gams = massin(25)

c               the factorization scale 
      QF = massin(13) 

c               real kinematics built in
      ti = tg + mg**2 - mi**2
      ts = tg + mg**2 - ms**2 

      ui = s4 - s - tg 
      us = ui + mi**2 - ms**2 
      ug = ui + mi**2 - mg**2 

      s4i     = s4 
      if (s4i.lt.0.D0) s4i = 0.D0

      s4s     = s4 + mi**2 - ms**2 

      hardfac = 1.D0 + mi**2/s4i

c               some denominators following wim's notes 
      spug   = s + ug 

c               following wim's notes  
      ushat = mg**2 - ms**2 - ti*ug/(s+ug)

c               wim's coupling structures 
      struc1 = real( Cl(1)*Clc(1) + Cr(1)*Crc(1) )
      struc2 = real( Cl(3)*Clc(3) + Cr(3)*Crc(3) )
      struc3 = mi*mg*real(  Cl(1)*Cl(3) + Clc(1)*Clc(3) 
     &                    + Cr(1)*Cr(3) + Crc(1)*Crc(3) ) 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 ra-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the s4s denominator [hauptwert]
      if ((s.ge.(ms+abs(mg))**2).and.(ms.ge.abs(mi))) then 
         s4s_d  = s4s/( s4s**2 + ms**2*gams**2 )
      else 
         s4s_d = 1.D0/s4s
      end if

c               the angular functions 
      call ANGULAR_ARRAY_NG_C(massin,NGANG)

c               wim's form output
c               color factors are Cf/Nc^2, Cf
c               replace s4s**(-1) -> s4s_d 
      gqbH = ushat**(-2)*e**2*gs**4*struc2*CA**(-2)*pi*Sn**2*spug**(-1)*
     + eight**(-1) * (  - 8*ug*ms**2 + 8*ug*mi**2 )
      gqbH = gqbH + ushat**(-2)*e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*
     + spug**(-4)*eight**(-1) * ( 32*ug**2*ti**3 )
      gqbH = gqbH + ushat**(-2)*e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*
     + spug**(-3)*eight**(-1) * ( 32*ug*ti**2*mi**2 - 32*ug*ti**2*mg**2
     +     + 32*ug**2*ti**2 )
      gqbH = gqbH + ushat**(-2)*e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*
     + spug**(-2)*eight**(-1) * ( 32*ug*ti*mi**2 - 32*ug*ti*mg**2 + 16*
     +    ug**2*ti )
      gqbH = gqbH + ushat**(-2)*e**2*gs**4*struc2*CA**(-1)*CF*pi*Sn**2*
     + spug**(-1)*eight**(-1) * ( 16*ug*mi**2 - 16*ug*mg**2 )
      gqbH = gqbH + ushat**(-2)*e**2*gs**4*struc2*pi*Sn**2*spug**(-1)*
     + eight**(-1) * ( 8*ug*ms**2 - 8*ug*mi**2 )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc2*CA**(-2)*pi*Sn**2*
     + spug**(-1)*eight**(-1) * (  - 8*ug )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc2*pi*Sn**2*spug**(-1)*
     + eight**(-1) * ( 8*ug )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc3*CA**(-2)*pi*Sn**2*
     + spug**(-1)*eight**(-1) * ( 8*s*ts**(-1) )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + spug**(-3)*eight**(-1) * ( 32*s*ti**2*ts**(-1) )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + spug**(-2)*eight**(-1) * ( 32*s*ti*ts**(-1) )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc3*CA**(-1)*CF*pi*Sn**2*
     + spug**(-1)*eight**(-1) * ( 16*s*ts**(-1) )
      gqbH = gqbH + ushat**(-1)*e**2*gs**4*struc3*pi*Sn**2*spug**(-1)*
     + eight**(-1) * (  - 8*s*ts**(-1) )
      gqbH = gqbH + e**2*gs**4*struc1*CA**(-2)*pi*Sn**2*spug**(-1)*
     + eight**(-1) * ( 8*s*tg*ts**(-2) + 8*tg*ug*ts**(-2) )
      gqbH = gqbH + e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*spug**(-2)*
     + eight**(-1) * (  - 32*tg*ti*ts**(-2)*mi**2 + 32*tg*ti*ts**(-2)*
     +    mg**2 + 32*tg**2*ti*ts**(-2) )
      gqbH = gqbH + e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*spug**(-1)*
     + eight**(-1) * (  - 32*tg*ts**(-2)*mi**2 + 32*tg*ts**(-2)*mg**2
     +     + 32*tg**2*ts**(-2) )
      gqbH = gqbH + e**2*gs**4*struc1*CA**(-1)*CF*pi*Sn**2*eight**(-1)
     +  * (  - 16*tg*ti**(-1)*ts**(-2)*mi**2 + 16*tg*ti**(-1)*ts**(-2)*
     +    mg**2 + 16*tg**2*ti**(-1)*ts**(-2) )
      gqbH = gqbH + e**2*gs**4*struc1*pi*Sn**2*spug**(-1)*eight**(-1)
     +  * (  - 8*s*tg*ts**(-2) - 8*tg*ug*ts**(-2) )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc1*CA**(-2)
     + *Sn**2*hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*tg*ui*ts**(-1) + 
     +    8*s**(-1)*tg**2*ts**(-1) - 4*s**(-1)*ug*s4i**(-1)*ms**(-2)*
     +    mi**4 + 4*s**(-1)*ug*s4i**(-1)*mi**2 + 4*s**(-1)*ug + 8*
     +    s**(-1)*s4i**(-1)*ts**(-1)*ms**(-2)*mi**4*mg**4 + 16*s**(-1)*
     +    s4i**(-1)*ts**(-1)*ms**2*mi**2*mg**2 + 8*s**(-1)*s4i**(-1)*
     +    ts**(-1)*ms**2*mi**4 - 8*s**(-1)*s4i**(-1)*ts**(-1)*ms**4*
     +    mi**2 - 8*s**(-1)*s4i**(-1)*ts**(-1)*mi**2*mg**4 - 16*s**(-1)
     +    *s4i**(-1)*ts**(-1)*mi**4*mg**2 - 8*s**(-1)*s4i**(-1)*
     +    ms**(-2)*mi**4*mg**2 - 8*s**(-1)*s4i**(-1)*ms**2*mi**2 + 8*
     +    s**(-1)*s4i**(-1)*mi**2*mg**2 + 8*s**(-1)*s4i**(-1)*mi**4 + 8
     +    *s**(-1)*ts**(-1)*ms**2*mi**2 + 24*s**(-1)*ts**(-1)*ms**2*
     +    mg**2 - 16*s**(-1)*ts**(-1)*ms**4 - 8*s**(-1)*ts**(-1)*mi**2*
     +    mg**2 - 8*s**(-1)*ts**(-1)*mg**4 - 16*s**(-1)*ms**2 + 8*
     +    s**(-1)*mi**2 + 8*s**(-1)*mg**2 + 4*s*tg*ts**(-2) + 4*s*
     +    ts**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc1*CA**(-2)
     + *Sn**2*hardfac**(-1)*eight**(-1) * ( 2*tg*ug*ts**(-2) + 2*tg*ui*
     +    ts**(-2) + 6*tg*ts**(-2)*mi**2 + 2*tg*ts**(-2)*mg**2 + 12*tg*
     +    ts**(-1) + 4*tg**2*ts**(-2) + 4*ui*ts**(-1) + 8*s4i**(-1)*
     +    ts**(-2)*ms**2*mi**2*mg**2 + 8*s4i**(-1)*ts**(-2)*ms**2*mi**4
     +     - 8*s4i**(-1)*ts**(-2)*ms**4*mi**2 - 8*s4i**(-1)*ts**(-2)*
     +    mi**4*mg**2 - 8*s4i**(-1)*ts**(-1)*ms**(-2)*mi**4*mg**2 - 8*
     +    s4i**(-1)*ts**(-1)*ms**2*mi**2 + 8*s4i**(-1)*ts**(-1)*mi**2*
     +    mg**2 + 8*s4i**(-1)*ts**(-1)*mi**4 + 4*ts**(-2)*ms**2*mi**2
     +     + 12*ts**(-2)*ms**2*mg**2 - 12*ts**(-2)*ms**4 - 4*ts**(-2)*
     +    mi**2*mg**2 - 12*ts**(-1)*ms**2 + 4*ts**(-1)*mi**2 + 8*
     +    ts**(-1)*mg**2 )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8 + 8*s**(-1)*ug**(-1)*
     +    s4i**(-1)*ms**(-2)*mi**4*mg**4 + 16*s**(-1)*ug**(-1)*
     +    s4i**(-1)*ms**2*mi**2*mg**2 + 8*s**(-1)*ug**(-1)*s4i**(-1)*
     +    ms**2*mi**4 - 8*s**(-1)*ug**(-1)*s4i**(-1)*ms**4*mi**2 - 8*
     +    s**(-1)*ug**(-1)*s4i**(-1)*mi**2*mg**4 - 16*s**(-1)*ug**(-1)*
     +    s4i**(-1)*mi**4*mg**2 + 16*s**(-1)*ug**(-1)*ms**2*mi**2 + 32*
     +    s**(-1)*ug**(-1)*ms**2*mg**2 - 24*s**(-1)*ug**(-1)*ms**4 - 16
     +    *s**(-1)*ug**(-1)*mi**2*mg**2 - 8*s**(-1)*ug**(-1)*mg**4 + 4*
     +    s**(-1)*ug*s4i**(-1)*ms**(-2)*mi**4 - 4*s**(-1)*ug*s4i**(-1)*
     +    mi**2 - 4*s**(-1)*ug + 8*s**(-1)*s4i**(-1)*ms**(-2)*mi**4*
     +    mg**2 + 8*s**(-1)*s4i**(-1)*ms**2*mi**2 - 8*s**(-1)*s4i**(-1)
     +    *mi**2*mg**2 - 8*s**(-1)*s4i**(-1)*mi**4 + 16*s**(-1)*ms**2
     +     - 8*s**(-1)*mi**2 - 8*s**(-1)*mg**2 + 4*s*tg*ug**(-1)*
     +    ts**(-1) - 4*s*tg*ts**(-2) + 4*s*ug**(-1)*ui*ts**(-1) + 8*s*
     +    ug**(-1)*s4i**(-1)*ms**(-2)*mi**4 )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s*ug**(-1)*s4i**(-1)*mi**2 - 
     +    8*s*ug**(-1) + 4*s**2*ug**(-1)*ts**(-1) - 4*tg*ug**(-1)*
     +    ts**(-1)*mi**2 - 4*tg*ug**(-1)*ts**(-1)*mg**2 - 2*tg*ug*
     +    ts**(-2) - 2*tg*ui*ts**(-2) - 6*tg*ts**(-2)*mi**2 - 2*tg*
     +    ts**(-2)*mg**2 - 4*tg*ts**(-1) - 4*tg**2*ts**(-2) + 16*
     +    ug**(-2)*s4i**(-1)*ms**(-2)*mi**4*mg**4 + 16*ug**(-2)*
     +    s4i**(-1)*ms**2*mi**2*mg**2 - 16*ug**(-2)*s4i**(-1)*mi**2*
     +    mg**4 - 16*ug**(-2)*s4i**(-1)*mi**4*mg**2 + 32*ug**(-2)*ms**2
     +    *mg**2 - 16*ug**(-2)*mi**2*mg**2 - 16*ug**(-2)*mg**4 + 8*
     +    ug**(-1)*s4i**(-1)*ts**(-1)*ms**(-2)*mi**4*mg**4 - 8*ug**(-1)
     +    *s4i**(-1)*ts**(-1)*ms**2*mi**4 + 8*ug**(-1)*s4i**(-1)*
     +    ts**(-1)*ms**4*mi**2 - 8*ug**(-1)*s4i**(-1)*ts**(-1)*mi**2*
     +    mg**4 + 16*ug**(-1)*s4i**(-1)*ms**(-2)*mi**4*mg**2 + 16*
     +    ug**(-1)*s4i**(-1)*ms**2*mi**2 - 16*ug**(-1)*s4i**(-1)*mi**2*
     +    mg**2 )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 16*ug**(-1)*s4i**(-1)*mi**4 - 8
     +    *ug**(-1)*ts**(-1)*ms**2*mi**2 + 16*ug**(-1)*ts**(-1)*ms**4
     +     - 8*ug**(-1)*ts**(-1)*mg**4 + 32*ug**(-1)*ms**2 - 16*
     +    ug**(-1)*mi**2 - 16*ug**(-1)*mg**2 - 8*s4i**(-1)*ts**(-2)*
     +    ms**2*mi**2*mg**2 - 8*s4i**(-1)*ts**(-2)*ms**2*mi**4 + 8*
     +    s4i**(-1)*ts**(-2)*ms**4*mi**2 + 8*s4i**(-1)*ts**(-2)*mi**4*
     +    mg**2 + 8*s4i**(-1)*ts**(-1)*ms**(-2)*mi**4*mg**2 + 8*
     +    s4i**(-1)*ts**(-1)*ms**2*mi**2 - 8*s4i**(-1)*ts**(-1)*mi**2*
     +    mg**2 - 8*s4i**(-1)*ts**(-1)*mi**4 + 8*s4i**(-1)*ms**(-2)*
     +    mi**4 - 8*s4i**(-1)*mi**2 - 4*ts**(-2)*ms**2*mi**2 - 12*
     +    ts**(-2)*ms**2*mg**2 + 12*ts**(-2)*ms**4 + 4*ts**(-2)*mi**2*
     +    mg**2 + 12*ts**(-1)*ms**2 - 4*ts**(-1)*mi**2 - 8*ts**(-1)*
     +    mg**2 )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc3*CA**(-2)
     + *Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ts**(-1) - 
     +    2*s**(-1)*ug*ts**(-1) - 6*ts**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ug**(-1) + 2*s*
     +    ug**(-1)*ts**(-1) + 4*tg*ug**(-1)*ts**(-1) + 16*ug**(-2)*
     +    mg**2 + 2*ug**(-1)*ui*ts**(-1) + 2*ug**(-1)*ts**(-1)*mi**2 + 
     +    6*ug**(-1)*ts**(-1)*mg**2 + 8*ug**(-1) + 4*ts**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*s4s*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ug**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ts**(-1)*ms**2 + 8*
     +    s**(-1)*ts**(-1)*mg**2 - 8*s**(-1) - 4*tg*ts**(-2) - 4*
     +    ts**(-2)*ms**2 + 4*ts**(-2)*mg**2 - 4*ts**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 24*s**(-1)*ug**(-1)*ms**2 + 8*
     +    s**(-1)*ug**(-1)*mi**2 + 16*s**(-1)*ug**(-1)*mg**2 + 8*
     +    s**(-1) + 4*tg*ts**(-2) + 16*ug**(-2)*mg**2 + 8*ug**(-1)*
     +    ts**(-1)*ms**2 + 16*ug**(-1) + 4*ts**(-2)*ms**2 - 4*ts**(-2)*
     +    mg**2 + 4*ts**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*ug**(-2)*ui + 4*ug**(-2)*mi**2
     +     - 4*ug**(-2)*mg**2 - 4*ug**(-1) )
      gqbH = gqbH + NGANG(0,0,0,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ug**(-1) + 16*
     +    ug**(-2)*mg**2 + 8*ug**(-1) )
      gqbH = gqbH + NGANG(2,0,-1,0)*s4s_d*e**2*gs**4*struc1*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ug*
     +    ui*ts**(-1) - 8*s**(-1)*tg**2*ui*ts**(-1) - 8*s**(-1)*tg**3*
     +    ts**(-1) + 8*s*tg*ug*ts**(-2) - 8*s*tg*ui*ts**(-2) - 8*s*
     +    tg**2*ts**(-2) - 4*s*ug*ts**(-1) + 4*s*ui*ts**(-1) + 8*tg*ug*
     +    ui*ts**(-2) - 4*tg*ug*ts**(-1) - 8*tg*ui**2*ts**(-2) + 8*
     +    tg**2*ug*ts**(-2) - 16*tg**2*ui*ts**(-2) - 8*tg**2*ts**(-1)
     +     - 8*tg**3*ts**(-2) - 4*ug*ui*ts**(-1) + 4*ui**2*ts**(-1) )
      gqbH = gqbH + NGANG(2,0,-1,0)*s4s_d*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s*tg*ug**(-1)*ui*ts**(-1) - 8*s*
     +    tg*ug*ts**(-2) + 8*s*tg*ui*ts**(-2) - 12*s*tg*ts**(-1) + 8*s*
     +    tg**2*ts**(-2) + 4*s*ug**(-1)*ui**2*ts**(-1) - 4*s*ui*
     +    ts**(-1) + 4*s**2*ug**(-1)*ui*ts**(-1) - 4*s**2*ts**(-1) + 8*
     +    tg*ug**(-1)*ui**2*ts**(-1) - 8*tg*ug*ui*ts**(-2) - 12*tg*ui*
     +    ts**(-1) + 8*tg*ui**2*ts**(-2) + 8*tg**2*ug**(-1)*ui*ts**(-1)
     +     - 8*tg**2*ug*ts**(-2) + 16*tg**2*ui*ts**(-2) - 8*tg**2*
     +    ts**(-1) + 8*tg**3*ts**(-2) )
      gqbH = gqbH + NGANG(2,0,-1,0)*s4s_d*e**2*gs**4*struc3*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * (  - 4 - 8*s**(-1)*tg
     +     - 4*s**(-1)*ms**2 + 4*s**(-1)*mi**2 - 2*s*ts**(-1) - 4*tg*
     +    ts**(-1) + 2*ug*ts**(-1) - 4*ui*ts**(-1) )
      gqbH = gqbH + NGANG(2,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2 - 2*s*ug**(-1) + 2*s*ts**(-1)
     +     - 4*tg*ug**(-1) + 4*tg*ts**(-1) + 4*ug**(-1)*ui - 4*ug**(-1)
     +    *ms**2 + 4*ug**(-1)*mi**2 - 2*ug*ts**(-1) + 4*ui*ts**(-1) )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s*tg*ts**(-2) - 4*tg*ug*
     +    ts**(-2) )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s*tg*ts**(-2) + 4*tg*ug*ts**(-2)
     +     )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8 - 8*s**(-1)*tg*ug**(-1)*ui - 8*
     +    s**(-1)*tg*ug**(-1)*mi**2 + 8*s**(-1)*tg*ug**(-1)*mg**2 - 16*
     +    s**(-1)*ug**(-1)*ui*ms**2 + 8*s**(-1)*ug**(-1)*ui*mi**2 + 8*
     +    s**(-1)*ug**(-1)*ui*mg**2 - 16*s**(-1)*ug**(-1)*ms**2*mi**2
     +     + 16*s**(-1)*ug**(-1)*ms**2*mg**2 + 8*s**(-1)*ug**(-1)*mi**4
     +     - 8*s**(-1)*ug**(-1)*mg**4 - 12*ug**(-2)*ui*ms**2 + 12*
     +    ug**(-2)*ui*mi**2 + 8*ug**(-2)*ui**2 - 4*ug**(-2)*us*ms**2 + 
     +    4*ug**(-2)*us*mi**2 - 8*ug**(-2)*ms**2*mi**2 + 16*ug**(-2)*
     +    ms**2*mg**2 - 4*ug**(-2)*ms**4 + 4*ug**(-2)*mi**4 - 8*
     +    ug**(-2)*mg**4 - 20*ug**(-1)*ui - 4*ug**(-1)*us + 12*ug**(-1)
     +    *ms**2 - 20*ug**(-1)*mi**2 + 8*ug**(-1)*mg**2 )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8 + 8*ug**(-1)*ui + 8*ug**(-1)*
     +    us - 8*ug**(-1)*ms**2 + 8*ug**(-1)*mi**2 )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ts**(-1) - 4*
     +    s**(-1)*ui*ts**(-1) - 4*s**(-1)*ts**(-1)*ms**2 + 4*s**(-1)*
     +    ts**(-1)*mg**2 - 4*ug**(-1)*ui*ts**(-1) - 4*ug**(-1)*ts**(-1)
     +    *mi**2 + 4*ug**(-1)*ts**(-1)*mg**2 - 4*ts**(-1) )
      gqbH = gqbH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2*s*ug**(-1)*ts**(-1) - 4*tg*
     +    ug**(-1)*ts**(-1) + 4*ug**(-1)*ui*ts**(-1) - 4*ug**(-1)*
     +    ts**(-1)*ms**2 + 4*ug**(-1)*ts**(-1)*mi**2 - 2*ts**(-1) )
      gqbH = gqbH + NGANG(2,4,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ug**(-1)*ui*ms**2
     +     + 8*s**(-1)*tg*ug**(-1)*ui*mg**2 + 16*s**(-1)*tg*ug**(-1)*
     +    ms**2*mg**2 - 8*s**(-1)*tg*ug**(-1)*ms**4 - 8*s**(-1)*tg*
     +    ug**(-1)*mg**4 - 4*s**(-1)*tg*ui + 32*s**(-1)*ug**(-1)*ui*
     +    ms**2*mg**2 - 16*s**(-1)*ug**(-1)*ui*ms**4 - 16*s**(-1)*
     +    ug**(-1)*ui*mg**4 - 8*s**(-1)*ug**(-1)*ui**2*ms**2 + 8*
     +    s**(-1)*ug**(-1)*ui**2*mg**2 - 24*s**(-1)*ug**(-1)*ms**2*
     +    mg**4 + 24*s**(-1)*ug**(-1)*ms**4*mg**2 - 8*s**(-1)*ug**(-1)*
     +    ms**6 + 8*s**(-1)*ug**(-1)*mg**6 - 4*s**(-1)*ui*ms**2 + 4*
     +    s**(-1)*ui*mg**2 - 4*s**(-1)*ui**2 + 16*ug**(-2)*ui*ms**2*
     +    mg**2 - 8*ug**(-2)*ui*ms**4 - 8*ug**(-2)*ui*mg**4 - 24*
     +    ug**(-2)*ms**2*mg**4 + 24*ug**(-2)*ms**4*mg**2 - 8*ug**(-2)*
     +    ms**6 + 8*ug**(-2)*mg**6 - 12*ug**(-1)*ui*ms**2 + 12*ug**(-1)
     +    *ui*mg**2 + 16*ug**(-1)*ms**2*mg**2 - 8*ug**(-1)*ms**4 - 8*
     +    ug**(-1)*mg**4 )
      gqbH = gqbH + NGANG(2,4,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*ui )
      gqbH = gqbH + NGANG(2,4,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ui*ts**(-1) - 8*
     +    s**(-1)*tg*ts**(-1)*ms**2 + 8*s**(-1)*tg*ts**(-1)*mg**2 - 4*
     +    s**(-1)*tg**2*ts**(-1) - 8*s**(-1)*ui*ts**(-1)*ms**2 + 8*
     +    s**(-1)*ui*ts**(-1)*mg**2 - 4*s**(-1)*ui**2*ts**(-1) + 8*
     +    s**(-1)*ts**(-1)*ms**2*mg**2 - 4*s**(-1)*ts**(-1)*ms**4 - 4*
     +    s**(-1)*ts**(-1)*mg**4 - 2*s*ug**(-1)*ts**(-1)*ms**2 + 2*s*
     +    ug**(-1)*ts**(-1)*mg**2 - 2*s*ts**(-1) - 4*tg*ug**(-1)*
     +    ts**(-1)*ms**2 + 4*tg*ug**(-1)*ts**(-1)*mg**2 - 6*tg*ts**(-1)
     +     - 4*ug**(-1)*ui*ts**(-1)*ms**2 + 4*ug**(-1)*ui*ts**(-1)*
     +    mg**2 + 8*ug**(-1)*ts**(-1)*ms**2*mg**2 - 4*ug**(-1)*ts**(-1)
     +    *ms**4 - 4*ug**(-1)*ts**(-1)*mg**4 - 6*ui*ts**(-1) - 6*
     +    ts**(-1)*ms**2 + 6*ts**(-1)*mg**2 )
      gqbH = gqbH + NGANG(2,6,-1,1)*s4s_d*e**2*gs**4*struc1*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ui*
     +    ts**(-1) - 8*s**(-1)*tg**2*ts**(-1) - 4*s**(-1)*ug*ui*
     +    ts**(-1) - 4*s*ts**(-1) - 8*tg*ts**(-1) - 4*ug*ts**(-1) - 4*
     +    ui*ts**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*s4s_d*e**2*gs**4*struc1*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s*tg*ug**(-1)*ts**(-1) - 4*s*
     +    ug**(-1)*ui*ts**(-1) - 4*s*ts**(-1) - 4*s**2*ug**(-1)*
     +    ts**(-1) - 8*tg*ug**(-1)*ui*ts**(-1) - 8*tg**2*ug**(-1)*
     +    ts**(-1) - 4*ui*ts**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*s4s_d*e**2*gs**4*struc3*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * (  - 4*s**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*ug**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ug**(-1)*ui - 8*
     +    s**(-1)*ug**(-1)*mi**2 + 8*s**(-1)*ug**(-1)*mg**2 - 8*
     +    ug**(-2)*ui - 8*ug**(-2)*mi**2 + 8*ug**(-2)*mg**2 + 8*
     +    ug**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*ug**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*ts**(-1) )
      gqbH = gqbH + NGANG(2,6,-1,1)*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*ug**(-1)*ts**(-1) )
      gqbH = gqbH + NGANG(2,8,-1,-2)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*ug**(-1)*ui*us*ms**2 - 4*
     +    ug**(-1)*ui*us*mi**2 + 8*ug**(-1)*us*ms**2*mi**2 - 4*ug**(-1)
     +    *us*ms**4 - 4*ug**(-1)*us*mi**4 + 4*ug**(-1)*us**2*ms**2 - 4*
     +    ug**(-1)*us**2*mi**2 + 4*ug*ms**2 - 4*ug*mi**2 - 2*ui*ms**2
     +     + 2*ui*mi**2 - 6*us*ms**2 + 6*us*mi**2 - 4*ms**2*mi**2 + 2*
     +    ms**4 + 2*mi**4 )
      gqbH = gqbH + NGANG(2,8,-1,-2)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*ug**(-1)*us**2*ms**2 + 8*
     +    ug**(-1)*us**2*mi**2 - 4*ug*ms**2 + 4*ug*mi**2 + 8*us*ms**2
     +     - 8*us*mi**2 )
      gqbH = gqbH + NGANG(2,8,-1,-1)*s4s_d*e**2*gs**4*struc3*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*
     +    ms**2 + 8*s**(-1)*tg*mi**2 - 4*s**(-1)*tg**2 + 8*s**(-1)*
     +    ms**2*mi**2 - 4*s**(-1)*ms**4 - 4*s**(-1)*mi**4 - 4*s*tg*
     +    ts**(-1) + 2*s*ug*ts**(-1) - s*ui*ts**(-1) - s*ts**(-1)*us - 
     +    3*s*ts**(-1)*ms**2 + 3*s*ts**(-1)*mi**2 + 2*tg*ug*ts**(-1) - 
     +    2*tg*ui*ts**(-1) - 2*tg*ts**(-1)*us - 6*tg*ts**(-1)*ms**2 + 6
     +    *tg*ts**(-1)*mi**2 - 4*tg - 4*tg**2*ts**(-1) + 2*ug*ts**(-1)*
     +    ms**2 - 2*ug*ts**(-1)*mi**2 - 2*ui*ts**(-1)*ms**2 + 2*ui*
     +    ts**(-1)*mi**2 - 2*ts**(-1)*us*ms**2 + 2*ts**(-1)*us*mi**2 + 
     +    4*ts**(-1)*ms**2*mi**2 - 2*ts**(-1)*ms**4 - 2*ts**(-1)*mi**4
     +     - 4*ms**2 + 4*mi**2 )
      gqbH = gqbH + NGANG(2,8,-1,-1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s*tg*ts**(-1) + 2*s*ug**(-1)*us
     +     - 2*s*ug*ts**(-1) + s*ui*ts**(-1) + s*ts**(-1)*us + 3*s*
     +    ts**(-1)*ms**2 - 3*s*ts**(-1)*mi**2 - 2*s + 4*tg*ug**(-1)*us
     +     - 2*tg*ug*ts**(-1) + 2*tg*ui*ts**(-1) + 2*tg*ts**(-1)*us + 6
     +    *tg*ts**(-1)*ms**2 - 6*tg*ts**(-1)*mi**2 - 2*tg + 4*tg**2*
     +    ts**(-1) + 4*ug**(-1)*us*ms**2 - 4*ug**(-1)*us*mi**2 - 2*ug*
     +    ts**(-1)*ms**2 + 2*ug*ts**(-1)*mi**2 + 2*ui*ts**(-1)*ms**2 - 
     +    2*ui*ts**(-1)*mi**2 + 2*ts**(-1)*us*ms**2 - 2*ts**(-1)*us*
     +    mi**2 - 4*ts**(-1)*ms**2*mi**2 + 2*ts**(-1)*ms**4 + 2*
     +    ts**(-1)*mi**4 - 2*ms**2 + 2*mi**2 )
      gqbH = gqbH + NGANG(2,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ug**(-1)*us*ms**2
     +     + 8*s**(-1)*tg*ug**(-1)*us*mi**2 + 4*s**(-1)*tg*us + 4*
     +    s**(-1)*tg*ms**2 - 4*s**(-1)*tg*mi**2 + 16*s**(-1)*ug**(-1)*
     +    us*ms**2*mi**2 - 8*s**(-1)*ug**(-1)*us*ms**4 - 8*s**(-1)*
     +    ug**(-1)*us*mi**4 + 4*s**(-1)*us*ms**2 - 4*s**(-1)*us*mi**2
     +     - 8*s**(-1)*ms**2*mi**2 + 4*s**(-1)*ms**4 + 4*s**(-1)*mi**4
     +     + 4*ug**(-2)*ui*us*ms**2 - 4*ug**(-2)*ui*us*mi**2 + 8*
     +    ug**(-2)*us*ms**2*mi**2 - 4*ug**(-2)*us*ms**4 - 4*ug**(-2)*us
     +    *mi**4 + 4*ug**(-2)*us**2*ms**2 - 4*ug**(-2)*us**2*mi**2 - 6*
     +    ug**(-1)*ui*ms**2 + 6*ug**(-1)*ui*mi**2 - 22*ug**(-1)*us*
     +    ms**2 + 22*ug**(-1)*us*mi**2 + 4*ug**(-1)*us**2 - 12*ug**(-1)
     +    *ms**2*mi**2 + 6*ug**(-1)*ms**4 + 6*ug**(-1)*mi**4 + 4*ug - 2
     +    *ui - 2*us + 14*ms**2 - 14*mi**2 )
      gqbH = gqbH + NGANG(2,8,-1,-1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 16*ug**(-1)*us*ms**2 - 16*ug**(-1)
     +    *us*mi**2 - 8*ug**(-1)*us**2 - 4*ug + 8*us - 8*ms**2 + 8*
     +    mi**2 )
      gqbH = gqbH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - s*ug**(-1)*ui*ts**(-1) - s*
     +    ug**(-1)*ts**(-1)*us + s*ug**(-1)*ts**(-1)*ms**2 - s*ug**(-1)
     +    *ts**(-1)*mi**2 - 2*s*ts**(-1) - 2*tg*ug**(-1)*ui*ts**(-1) - 
     +    2*tg*ug**(-1)*ts**(-1)*us + 2*tg*ug**(-1)*ts**(-1)*ms**2 - 2*
     +    tg*ug**(-1)*ts**(-1)*mi**2 + 2*tg*ts**(-1) - 2*ug**(-1)*ui*
     +    ts**(-1)*ms**2 + 2*ug**(-1)*ui*ts**(-1)*mi**2 - 2*ug**(-1)*
     +    ts**(-1)*us*ms**2 + 2*ug**(-1)*ts**(-1)*us*mi**2 - 4*ug**(-1)
     +    *ts**(-1)*ms**2*mi**2 + 2*ug**(-1)*ts**(-1)*ms**4 + 2*
     +    ug**(-1)*ts**(-1)*mi**4 + 2*ts**(-1)*ms**2 - 2*ts**(-1)*mi**2
     +     )
      gqbH = gqbH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*s*ug**(-1)*ts**(-1)*us + 2*s*
     +    ts**(-1) + 4*tg*ug**(-1)*ts**(-1)*us - 2*tg*ts**(-1) + 4*
     +    ug**(-1)*ts**(-1)*us*ms**2 - 4*ug**(-1)*ts**(-1)*us*mi**2 - 2
     +    *ts**(-1)*ms**2 + 2*ts**(-1)*mi**2 )
      gqbH = gqbH + NGANG(3,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ms**(-2)*mi**2*mg**2
     +     + 8*s**(-1)*mg**2 )
      gqbH = gqbH + NGANG(3,0,-1,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s**(-1)*ms**(-2)*mi**2*mg**2 - 8
     +    *s**(-1)*mg**2 )
      gqbH = gqbH + NGANG(3,7,-1,1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*ms**(-2)*mg**2 )
      gqbH = gqbH + NGANG(3,7,-1,1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*ms**(-2)*mg**2 )
      gqbH = gqbH + NGANG(3,8,-1,-2)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*ms**2*mg**2 - 8*mi**2*mg**2 )
      gqbH = gqbH + NGANG(3,8,-1,-2)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*ms**2*mg**2 + 8*mi**2*mg**2 )
      gqbH = gqbH + NGANG(3,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s**(-1)*ms**(-2)*mi**4*mg**2 + 8
     +    *s**(-1)*ms**2*mg**2 - 16*s**(-1)*mi**2*mg**2 - 8*ms**(-2)*
     +    mi**2*mg**2 + 8*mg**2 )
      gqbH = gqbH + NGANG(3,8,-1,-1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*ms**(-2)*mi**4*mg**2
     +     - 8*s**(-1)*ms**2*mg**2 + 16*s**(-1)*mi**2*mg**2 + 8*
     +    ms**(-2)*mi**2*mg**2 - 8*mg**2 )
      gqbH = gqbH + NGANG(4,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2 - 4*s**(-1)*tg - 8*s**(-1)*
     +    ug**(-1)*ui*ms**2 + 8*s**(-1)*ug**(-1)*ui*mg**2 - 8*s**(-1)*
     +    ug**(-1)*ms**2*mi**2 + 8*s**(-1)*ug**(-1)*ms**2*mg**2 + 8*
     +    s**(-1)*ug**(-1)*mi**2*mg**2 - 8*s**(-1)*ug**(-1)*mg**4 - 4*
     +    s**(-1)*ui + 8*s**(-1)*ms**(-2)*mi**2*mg**2 - 8*s**(-1)*ms**2
     +     + 4*s**(-1)*mi**2 - 4*s**(-1)*mg**2 - 4*ug**(-2)*ui*ms**2 + 
     +    4*ug**(-2)*ui*mg**2 - 4*ug**(-2)*ms**2*mi**2 + 4*ug**(-2)*
     +    ms**2*mg**2 + 4*ug**(-2)*mi**2*mg**2 - 4*ug**(-2)*mg**4 - 2*
     +    ug**(-1)*ui + 4*ug**(-1)*ms**2 - 2*ug**(-1)*mi**2 - 2*
     +    ug**(-1)*mg**2 )
      gqbH = gqbH + NGANG(4,0,-1,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2 - 8*s**(-1)*tg*ug**(-1)*ms**2
     +     + 4*s**(-1)*tg*ug**(-1)*mi**2 + 4*s**(-1)*tg*ug**(-1)*mg**2
     +     - 8*s**(-1)*ms**(-2)*mi**2*mg**2 + 4*s**(-1)*ms**2 + 4*
     +    s**(-1)*mg**2 + 4*tg*ug**(-1) + 2*ug**(-1)*ui + 4*ug**(-1)*
     +    ms**2 + 2*ug**(-1)*mi**2 + 2*ug**(-1)*mg**2 )
      gqbH = gqbH + NGANG(4,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ts**(-1) - 4*
     +    s**(-1)*ui*ts**(-1) - 4*s**(-1)*ts**(-1)*ms**2 + 4*s**(-1)*
     +    ts**(-1)*mg**2 - 2*ug**(-1)*ui*ts**(-1) - 2*ug**(-1)*ts**(-1)
     +    *mi**2 + 2*ug**(-1)*ts**(-1)*mg**2 )
      gqbH = gqbH + NGANG(4,7,-1,1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*ms**(-2)*mg**2 + 4*
     +    s**(-1) )
      gqbH = gqbH + NGANG(4,7,-1,1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*tg*ug**(-1) - 4*s**(-1)*
     +    ms**(-2)*mg**2 )
      gqbH = gqbH + NGANG(4,8,-1,-2)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*ug**(-1)*ui*ms**2*mi**2 + 4*
     +    ug**(-1)*ui*ms**2*mg**2 - 4*ug**(-1)*ui*ms**4 - 4*ug**(-1)*ui
     +    *mi**2*mg**2 + 16*ug**(-1)*ms**2*mi**2*mg**2 + 4*ug**(-1)*
     +    ms**2*mi**4 + 4*ug**(-1)*ms**2*mg**4 - 12*ug**(-1)*ms**4*
     +    mi**2 - 12*ug**(-1)*ms**4*mg**2 + 8*ug**(-1)*ms**6 - 4*
     +    ug**(-1)*mi**2*mg**4 - 4*ug**(-1)*mi**4*mg**2 + 4*ms**2*mi**2
     +     - 4*ms**2*mg**2 - 4*ms**4 + 4*mi**2*mg**2 )
      gqbH = gqbH + NGANG(4,8,-1,-2)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*ug**(-1)*ui*ms**2*mi**2 - 4*
     +    ug**(-1)*ui*ms**2*mg**2 + 4*ug**(-1)*ui*ms**4 + 4*ug**(-1)*ui
     +    *mi**2*mg**2 - 4*ug**(-1)*ms**2*mi**4 - 4*ug**(-1)*ms**2*
     +    mg**4 - 4*ug**(-1)*ms**4*mi**2 - 4*ug**(-1)*ms**4*mg**2 + 8*
     +    ug**(-1)*ms**6 + 4*ug**(-1)*mi**2*mg**4 + 4*ug**(-1)*mi**4*
     +    mg**2 - 4*ms**2*mi**2 + 4*ms**2*mg**2 + 4*ms**4 - 4*mi**2*
     +    mg**2 )
      gqbH = gqbH + NGANG(4,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 8*s**(-1)*tg*ug**(-1)*ms**2*
     +    mi**2 - 8*s**(-1)*tg*ug**(-1)*ms**2*mg**2 + 8*s**(-1)*tg*
     +    ug**(-1)*ms**4 + 8*s**(-1)*tg*ug**(-1)*mi**2*mg**2 + 16*
     +    s**(-1)*ug**(-1)*ms**2*mi**2*mg**2 + 8*s**(-1)*ug**(-1)*ms**2
     +    *mi**4 - 16*s**(-1)*ug**(-1)*ms**4*mi**2 - 8*s**(-1)*ug**(-1)
     +    *ms**4*mg**2 + 8*s**(-1)*ug**(-1)*ms**6 - 8*s**(-1)*ug**(-1)*
     +    mi**4*mg**2 - 8*s**(-1)*ms**(-2)*mi**4*mg**2 + 12*s**(-1)*
     +    ms**2*mi**2 - 4*s**(-1)*ms**2*mg**2 - 8*s**(-1)*ms**4 + 12*
     +    s**(-1)*mi**2*mg**2 - 4*s**(-1)*mi**4 + 4*ug**(-2)*ui*ms**2*
     +    mi**2 + 4*ug**(-2)*ui*ms**2*mg**2 - 4*ug**(-2)*ui*ms**4 - 4*
     +    ug**(-2)*ui*mi**2*mg**2 + 16*ug**(-2)*ms**2*mi**2*mg**2 + 4*
     +    ug**(-2)*ms**2*mi**4 + 4*ug**(-2)*ms**2*mg**4 - 12*ug**(-2)*
     +    ms**4*mi**2 - 12*ug**(-2)*ms**4*mg**2 + 8*ug**(-2)*ms**6 - 4*
     +    ug**(-2)*mi**2*mg**4 - 4*ug**(-2)*mi**4*mg**2 - 2*ug**(-1)*ui
     +    *ms**2 )
      gqbH = gqbH + NGANG(4,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*ug**(-1)*ui*mi**2 - 18*ug**(-1)*
     +    ms**2*mi**2 - 22*ug**(-1)*ms**2*mg**2 + 20*ug**(-1)*ms**4 + 
     +    14*ug**(-1)*mi**2*mg**2 + 2*ug**(-1)*mi**4 + 4*ug**(-1)*mg**4
     +     + 8*ms**(-2)*mi**2*mg**2 + 2*ms**2 + 2*mi**2 - 4*mg**2 )
      gqbH = gqbH + NGANG(4,8,-1,-1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s**(-1)*tg*ug**(-1)*ms**2*mi**2
     +     + 8*s**(-1)*tg*ug**(-1)*ms**2*mg**2 - 8*s**(-1)*tg*ug**(-1)*
     +    ms**4 - 8*s**(-1)*tg*ug**(-1)*mi**2*mg**2 + 8*s**(-1)*
     +    ms**(-2)*mi**4*mg**2 - 12*s**(-1)*ms**2*mi**2 + 4*s**(-1)*
     +    ms**2*mg**2 + 8*s**(-1)*ms**4 - 12*s**(-1)*mi**2*mg**2 + 4*
     +    s**(-1)*mi**4 + 2*ug**(-1)*ui*ms**2 - 2*ug**(-1)*ui*mi**2 - 2
     +    *ug**(-1)*ms**2*mi**2 + 2*ug**(-1)*ms**2*mg**2 + 8*ug**(-1)*
     +    ms**4 - 2*ug**(-1)*mi**2*mg**2 - 2*ug**(-1)*mi**4 - 4*
     +    ug**(-1)*mg**4 - 8*ms**(-2)*mi**2*mg**2 - 2*ms**2 - 2*mi**2
     +     + 4*mg**2 )
      gqbH = gqbH + NGANG(4,8,-1,-1)*e**2*gs**4*struc3*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - s*ug**(-1)*ui*ts**(-1) + 2*s*
     +    ug**(-1)*ts**(-1)*ms**2 - s*ug**(-1)*ts**(-1)*mi**2 - s*
     +    ug**(-1)*ts**(-1)*mg**2 + s*ts**(-1) - 2*tg*ug**(-1)*ui*
     +    ts**(-1) + 4*tg*ug**(-1)*ts**(-1)*ms**2 - 2*tg*ug**(-1)*
     +    ts**(-1)*mi**2 - 2*tg*ug**(-1)*ts**(-1)*mg**2 - 2*ug**(-1)*ui
     +    *ts**(-1)*ms**2 + 2*ug**(-1)*ui*ts**(-1)*mi**2 - 6*ug**(-1)*
     +    ts**(-1)*ms**2*mi**2 - 2*ug**(-1)*ts**(-1)*ms**2*mg**2 + 4*
     +    ug**(-1)*ts**(-1)*ms**4 + 2*ug**(-1)*ts**(-1)*mi**2*mg**2 + 2
     +    *ug**(-1)*ts**(-1)*mi**4 + 2*ts**(-1)*ms**2 - 2*ts**(-1)*
     +    mi**2 )
      gqbH = gqbH + NGANG(7,8,1,-2)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*ug**(-1)*ms**2 - 8*ug**(-1)*
     +    mi**2 )
      gqbH = gqbH + NGANG(7,8,1,-1)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 2*s**(-1)*tg*ug**(-1) + 2*s*
     +    ug**(-1)*ts**(-1) + 2*tg*ug**(-1)*ts**(-1) - 2*ug**(-1) )
      gqbH = gqbH + NGANG(7,8,1,-1)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ug**(-1) + 4*
     +    ug**(-1) )
      gqbH = gqbH + NGANG(8,0,-2,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*ug**(-1)*us*ms**2 - 4*ug**(-1)*
     +    us*mi**2 - 4*ug**(-1)*ms**2*mi**2 - 4*ug**(-1)*ms**2*mg**2 + 
     +    4*ug**(-1)*ms**4 + 4*ug**(-1)*mi**2*mg**2 - 4*ms**2 + 4*mi**2
     +     )
      gqbH = gqbH + NGANG(8,0,-2,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 8*s*ug**(-1)*ms**2 - 8*s*ug**(-1)*
     +    mi**2 - 16*ug**(-2)*ms**2*mi**2*mg**2 - 16*ug**(-2)*ms**2*
     +    mg**4 + 16*ug**(-2)*ms**4*mg**2 + 16*ug**(-2)*mi**2*mg**4 + 2
     +    *ug**(-1)*ui*ms**2 - 2*ug**(-1)*ui*mi**2 - 8*ug**(-1)*us*
     +    ms**2 + 8*ug**(-1)*us*mi**2 - 6*ug**(-1)*ms**2*mi**2 - 10*
     +    ug**(-1)*ms**2*mg**2 + 8*ug**(-1)*ms**4 + 10*ug**(-1)*mi**2*
     +    mg**2 - 2*ug**(-1)*mi**4 + 6*ms**2 - 6*mi**2 )
      gqbH = gqbH + NGANG(8,0,-1,0)*s4s_d*e**2*gs**4*struc3*
     + CA**(-2)*Sn**2*hardfac**(-1)*eight**(-1) * ( 2 + 2*s**(-1)*tg + 
     +    2*s*ts**(-1) + 2*tg*ts**(-1) )
      gqbH = gqbH + NGANG(8,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 4*s**(-1)*tg*ug**(-1)*ms**2 + 4
     +    *s**(-1)*tg*ug**(-1)*mi**2 - 2*s**(-1)*tg - 4*s**(-1)*tg**2*
     +    ug**(-1) + 2*s**(-1)*ms**2 - 2*s**(-1)*mi**2 + 6*s*tg*
     +    ug**(-1)*ts**(-1) + 16*s*ug**(-2)*mg**2 + s*ug**(-1)*ui*
     +    ts**(-1) + 2*s*ug**(-1)*ts**(-1)*ms**2 + s*ug**(-1)*ts**(-1)*
     +    mi**2 + 5*s*ug**(-1)*ts**(-1)*mg**2 + 10*s*ug**(-1) + 3*s*
     +    ts**(-1) + 2*s**2*ug**(-1)*ts**(-1) + 16*tg*ug**(-2)*mg**2 + 
     +    2*tg*ug**(-1)*ui*ts**(-1) + 4*tg*ug**(-1)*ts**(-1)*ms**2 - 2*
     +    tg*ug**(-1)*ts**(-1)*mi**2 + 6*tg*ug**(-1)*ts**(-1)*mg**2 + 6
     +    *tg*ug**(-1) + 4*tg*ts**(-1) + 4*tg**2*ug**(-1)*ts**(-1) + 16
     +    *ug**(-2)*ms**2*mg**2 - 16*ug**(-2)*mi**2*mg**2 + 2*ug**(-1)*
     +    ui*ts**(-1)*ms**2 - 2*ug**(-1)*ui*ts**(-1)*mi**2 + 2*ug**(-1)
     +    *ts**(-1)*ms**2*mi**2 + 6*ug**(-1)*ts**(-1)*ms**2*mg**2 - 6*
     +    ug**(-1)*ts**(-1)*mi**2*mg**2 - 2*ug**(-1)*ts**(-1)*mi**4 + 
     +    14*ug**(-1)*ms**2 )
      gqbH = gqbH + NGANG(8,0,-1,0)*s4s_d*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * (  - 16*ug**(-1)*mi**2 + 2*ug**(-1)*
     +    mg**2 + 4*ts**(-1)*ms**2 - 4*ts**(-1)*mi**2 )
      gqbH = gqbH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 4*s**(-1)*tg + 4*ug**(-2)*us*ms**2
     +     - 4*ug**(-2)*us*mi**2 - 4*ug**(-2)*ms**2*mi**2 - 4*ug**(-2)*
     +    ms**2*mg**2 + 4*ug**(-2)*ms**4 + 4*ug**(-2)*mi**2*mg**2 + 4*
     +    ug**(-1)*us + 4*ug**(-1)*mi**2 - 4*ug**(-1)*mg**2 )
      gqbH = gqbH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 6 - 8*s**(-1)*tg*ug**(-1)*ms**2 + 
     +    4*s**(-1)*tg*ug**(-1)*mi**2 + 4*s**(-1)*tg*ug**(-1)*mg**2 + 4
     +    *s**(-1)*ms**2 - 4*s**(-1)*mi**2 + 4*s*ug**(-1) - 4*tg*
     +    ug**(-1) + 32*ug**(-2)*ms**2*mg**2 - 16*ug**(-2)*mi**2*mg**2
     +     - 16*ug**(-2)*mg**4 - 2*ug**(-1)*ui - 8*ug**(-1)*us + 12*
     +    ug**(-1)*ms**2 - 14*ug**(-1)*mi**2 - 6*ug**(-1)*mg**2 )
      gqbH = gqbH + NGANG(8,0,-1,0)*e**2*gs**4*struc3*Sn**2*
     + hardfac**(-1)*eight**(-1) * ( 2*s*ug**(-1)*ts**(-1) - 2*tg*
     +    ug**(-1)*ts**(-1) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF*pi*Sn**2*spug**(-4)*eight**(-1) * ( 32*ug**2*
     +    ti**3 )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * ( 32*ug*
     +    ti**2*mi**2 - 32*ug*ti**2*mg**2 + 32*ug**2*ti**2 )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * ( 32*ug*ti*
     +    mi**2 - 32*ug*ti*mg**2 + 16*ug**2*ti )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * ( 16*ug*
     +    mi**2 - 16*ug*mg**2 )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * ( 32*s*
     +    ti**2*ts**(-1) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * ( 32*s*ti*
     +    ts**(-1) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * ( 16*s*
     +    ts**(-1) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * (  - 32*tg*ti*
     +    ts**(-2)*mi**2 + 32*tg*ti*ts**(-2)*mg**2 + 32*tg**2*ti*
     +    ts**(-2) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * (  - 32*tg*
     +    ts**(-2)*mi**2 + 32*tg*ts**(-2)*mg**2 + 32*tg**2*ts**(-2) )
      gqbH = gqbH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc1*
     + CA**(-1)*CF*pi*Sn**2*eight**(-1) * (  - 16*tg*ti**(-1)*ts**(-2)*
     +    mi**2 + 16*tg*ti**(-1)*ts**(-2)*mg**2 + 16*tg**2*ti**(-1)*
     +    ts**(-2) )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-4)*eight**(-1) * (  - 32*ug**2*
     +    ti**3 )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * (  - 32*ug*ti**2*
     +    mi**2 + 32*ug*ti**2*mg**2 - 32*ug**2*ti**2 )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * (  - 32*ug*ti*
     +    mi**2 + 32*ug*ti*mg**2 - 16*ug**2*ti )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * (  - 16*ug*mi**2
     +     + 16*ug*mg**2 )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * (  - 32*s*ti**2*
     +    ts**(-1) )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * (  - 32*s*ti*
     +    ts**(-1) )
      gqbH = gqbH + log(s4i*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * (  - 16*s*ts**(-1)
     +     )
      gqbH = gqbH + log(s4i*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*pi
     + *Sn**2*spug**(-2)*eight**(-1) * ( 32*tg*ti*ts**(-2)*mi**2 - 32*
     +    tg*ti*ts**(-2)*mg**2 - 32*tg**2*ti*ts**(-2) )
      gqbH = gqbH + log(s4i*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*pi
     + *Sn**2*spug**(-1)*eight**(-1) * ( 32*tg*ts**(-2)*mi**2 - 32*tg*
     +    ts**(-2)*mg**2 - 32*tg**2*ts**(-2) )
      gqbH = gqbH + log(s4i*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*pi
     + *Sn**2*eight**(-1) * ( 16*tg*ti**(-1)*ts**(-2)*mi**2 - 16*tg*
     +    ti**(-1)*ts**(-2)*mg**2 - 16*tg**2*ti**(-1)*ts**(-2) )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-4)*eight**(-1) * ( 32*ug**2*ti**3 )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * ( 32*ug*ti**2*
     +    mi**2 - 32*ug*ti**2*mg**2 + 32*ug**2*ti**2 )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * ( 32*ug*ti*mi**2
     +     - 32*ug*ti*mg**2 + 16*ug**2*ti )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * ( 16*ug*mi**2 - 16
     +    *ug*mg**2 )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-3)*eight**(-1) * ( 32*s*ti**2*
     +    ts**(-1) )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-2)*eight**(-1) * ( 32*s*ti*ts**(-1)
     +     )
      gqbH = gqbH + log(QF**2*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF*pi*Sn**2*spug**(-1)*eight**(-1) * ( 16*s*ts**(-1) )
      gqbH = gqbH + log(QF**2*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*
     + pi*Sn**2*spug**(-2)*eight**(-1) * (  - 32*tg*ti*ts**(-2)*mi**2
     +     + 32*tg*ti*ts**(-2)*mg**2 + 32*tg**2*ti*ts**(-2) )
      gqbH = gqbH + log(QF**2*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*
     + pi*Sn**2*spug**(-1)*eight**(-1) * (  - 32*tg*ts**(-2)*mi**2 + 32
     +    *tg*ts**(-2)*mg**2 + 32*tg**2*ts**(-2) )
      gqbH = gqbH + log(QF**2*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF*
     + pi*Sn**2*eight**(-1) * (  - 16*tg*ti**(-1)*ts**(-2)*mi**2 + 16*
     +    tg*ti**(-1)*ts**(-2)*mg**2 + 16*tg**2*ti**(-1)*ts**(-2) )

c               the prefactor for the scaling functions 
      NG_GBH = gqbH * (abs(mi)+mg)**2/4.D0

      end
