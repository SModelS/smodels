cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     THE SCALING FUNCTIONS                                            c
c                                                                      c
c     NG_QBH(MASSIN,CL,CR)                                             c
c                                                                      c
c     INPUT :                                                          c
c                                                                      c
c       MASSIN(1)  = s                                                 c
c       MASSIN(2)  = tg                                                c
c       MASSIN(3)  = s4                                                c
c       MASSIN(6)  = m1                                                c
c       MASSIN(7)  = mg                                                c
c       MASSIN(11) = ms                                                c
c       MASSIN(13) = qf                                                c
c                                                                      c
c       CL(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c       CR(1:4)  HIGGSINO/GAUGINO-QUARK-SQUARK COUPLING [COMPLEX]      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
      real*8 function NG_QBH(massin,Cl,Cr)

      implicit none 

      integer    n
      real*8     massin(1:30),Pi,CF,CA,Sn,e,gs,struc1,struc2,struc3
     &          ,four,QF,s,tg,s4,mi,mg,ms,ti,ts,ui,us,ug,tshat,ushat
     &          ,sptg,spug,s4i,s4imts,s4imus,sfac,hardfac
     &          ,qqH,NGANG(0:9,0:9,-2:2,-2:2)
      complex*16 Cl(4),Cr(4),Clc(4),Crc(4)

      pi = 4.D0*atan(1.D0)
      CF = 4.D0/3.D0
      CA = 3.D0
      Sn = 1.D0/(16.D0 * Pi**2)

c               the denominators 
      four = 4.D0

      do n=1,4 
         Clc(n) =  conjg( Cl(n) )
         Crc(n) =  conjg( Cr(n) )
      end do

      s  = massin(1)
      tg = massin(2)
      s4 = massin(3)
      mi = massin(6)
      mg = massin(7)
      ms = massin(11)

c               the factorization scale 
      QF = massin(13) 

c               real kinematics built in
      ti = tg + mg**2 - mi**2
      ts = tg + mg**2 - ms**2 

      ui = s4 - s - tg 
      us = ui + mi**2 - ms**2 
      ug = ui + mi**2 - mg**2 

      s4i     = s4 
      hardfac = 1.D0 + mi**2/s4i

c               some denominators following wim's notes 
      sptg   = s + tg 
      spug   = s + ug 
      s4imts = s + ug + ms**2 - mi**2 
      s4imus = s + tg + ms**2 - mi**2 
      sfac   = s + 2.D0*ms**2 - mi**2 - mg**2 

c               following wim's notes  
      tshat = mg**2 - ms**2 - ui*tg/(s+tg) 
      ushat = mg**2 - ms**2 - ti*ug/(s+ug)

c               wim's coupling structures 
      struc1 = real( Cl(1)*Clc(1) + Cr(1)*Crc(1) )
      struc2 = real( Cl(3)*Clc(3) + Cr(3)*Crc(3) )
      struc3 = mi*mg*real(  Cl(1)*Cl(3) + Clc(1)*Clc(3) 
     &                     + Cr(1)*Cr(3) + Crc(1)*Crc(3) ) 

c               gs**2=4*pi*alpha_s is cut, re-appears as nlo later 
c               gs**2*e**2 ra-appears in the over-all factor
      gs = 1.D0
      e  = 1.D0

c               the angular functions 
      call ANGULAR_ARRAY_NG_R(massin,NGANG)

c               wim's form output
c               color factors are Cf/Nc^2, Cf
      qqH =tshat**(-2)*e**2*gs**4*struc1*CA**(-2)*CF*pi*Sn**2*sptg**(-1)
     + *four**(-1) * (  - 8*ts*ms**2 + 8*ts*mi**2 )
      qqH = qqH + tshat**(-2)*e**2*gs**4*struc1*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * ( 8*ts*ms**2 - 8*ts*mi**2 )
      qqH = qqH + tshat**(-1)*e**2*gs**4*struc1*CA**(-2)*CF*pi*Sn**2*
     + sptg**(-1)*four**(-1) * (  - 8*ts + 8*ms**2 - 8*mi**2 )
      qqH = qqH + tshat**(-1)*e**2*gs**4*struc1*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * ( 8*ts - 8*ms**2 + 8*mi**2 )
      qqH = qqH + tshat**(-1)*e**2*gs**4*struc3*CA**(-2)*CF*pi*Sn**2*
     + sptg**(-1)*four**(-1) * ( 8*s*us**(-1) + 8*ug*us**(-1) + 8*
     +    us**(-1)*ms**2 - 8*us**(-1)*mi**2 )
      qqH = qqH + tshat**(-1)*e**2*gs**4*struc3*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * (  - 8*s*us**(-1) - 8*ug*us**(-1) - 8*us**(-1)*
     +    ms**2 + 8*us**(-1)*mi**2 )
      qqH = qqH + ushat**(-2)*e**2*gs**4*struc2*CA**(-2)*CF*pi*Sn**2*
     + spug**(-1)*four**(-1) * (  - 8*us*ms**2 + 8*us*mi**2 )
      qqH = qqH + ushat**(-2)*e**2*gs**4*struc2*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * ( 8*us*ms**2 - 8*us*mi**2 )
      qqH = qqH + ushat**(-1)*e**2*gs**4*struc2*CA**(-2)*CF*pi*Sn**2*
     + spug**(-1)*four**(-1) * (  - 8*us + 8*ms**2 - 8*mi**2 )
      qqH = qqH + ushat**(-1)*e**2*gs**4*struc2*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * ( 8*us - 8*ms**2 + 8*mi**2 )
      qqH = qqH + ushat**(-1)*e**2*gs**4*struc3*CA**(-2)*CF*pi*Sn**2*
     + spug**(-1)*four**(-1) * ( 8*s*ts**(-1) + 8*tg*ts**(-1) + 8*
     +    ts**(-1)*ms**2 - 8*ts**(-1)*mi**2 )
      qqH = qqH + ushat**(-1)*e**2*gs**4*struc3*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * (  - 8*s*ts**(-1) - 8*tg*ts**(-1) - 8*ts**(-1)*
     +    ms**2 + 8*ts**(-1)*mi**2 )
      qqH = qqH + e**2*gs**4*struc1*CA**(-2)*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * ( 8 )
      qqH = qqH + e**2*gs**4*struc1*CA**(-2)*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * ( 8*s*tg*ts**(-2) + 8*tg*ui*ts**(-2) + 8*tg**2*
     +    ts**(-2) )
      qqH = qqH + e**2*gs**4*struc1*CF*pi*Sn**2*sptg**(-1)*four**(-1)
     +  * (  - 8 )
      qqH = qqH + e**2*gs**4*struc1*CF*pi*Sn**2*spug**(-1)*four**(-1)
     +  * (  - 8*s*tg*ts**(-2) - 8*tg*ui*ts**(-2) - 8*tg**2*ts**(-2) )
      qqH = qqH + e**2*gs**4*struc2*CA**(-2)*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * ( 8*s*ug*us**(-2) + 8*tg*ug*us**(-2) + 8*ug*ui*
     +    us**(-2) )
      qqH = qqH + e**2*gs**4*struc2*CA**(-2)*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * ( 8 )
      qqH = qqH + e**2*gs**4*struc2*CF*pi*Sn**2*sptg**(-1)*four**(-1)
     +  * (  - 8*s*ug*us**(-2) - 8*tg*ug*us**(-2) - 8*ug*ui*us**(-2) )
      qqH = qqH + e**2*gs**4*struc2*CF*pi*Sn**2*spug**(-1)*four**(-1)
     +  * (  - 8 )
      qqH = qqH + e**2*gs**4*struc3*CA**(-2)*CF*pi*Sn**2*sptg**(-1)*
     + four**(-1) * ( 8*us**(-1) )
      qqH = qqH + e**2*gs**4*struc3*CA**(-2)*CF*pi*Sn**2*spug**(-1)*
     + four**(-1) * ( 8*ts**(-1) )
      qqH = qqH + e**2*gs**4*struc3*CF*pi*Sn**2*sptg**(-1)*four**(-1)
     +  * (  - 8*us**(-1) )
      qqH = qqH + e**2*gs**4*struc3*CF*pi*Sn**2*spug**(-1)*four**(-1)
     +  * (  - 8*ts**(-1) )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 2*s*tg*ts**(-2) - 2
     +    *tg*ug*ts**(-2) - 2*tg*ts**(-2)*ms**2 + 2*tg*ts**(-2)*mi**2 )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*tg*ts**(-2) )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 2*s*tg*ts**(-2) + 2*tg
     +    *ug*ts**(-2) + 2*tg*ts**(-2)*ms**2 - 2*tg*ts**(-2)*mi**2 )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*tg*ts**(-2) )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 2*s*ug*us**(-2) - 2
     +    *tg*ug*us**(-2) - 2*ug*us**(-2)*ms**2 + 2*ug*us**(-2)*mi**2 )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*ug*us**(-2) )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 2*s*ug*us**(-2) + 2*tg
     +    *ug*us**(-2) + 2*ug*us**(-2)*ms**2 - 2*ug*us**(-2)*mi**2 )
      qqH = qqH + NGANG(0,0,0,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*ug*us**(-2) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 4*s*ug*ts**(-1) + 4*s*
     +    ts**(-1)*ms**2 - 4*s*ts**(-1)*mi**2 + 4*s**2*ts**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4 + 4*tg*ts**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4 )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 2*s*tg*ug*us**(-2)
     +     + 8*s*tg*us**(-1) + 2*s*ug*ti*us**(-2) + 6*s*ug*ui*us**(-2)
     +     - 4*s*ug*us**(-1) + 2*s*ug**2*us**(-2) + 4*s**2*us**(-1) + 2
     +    *tg*ug*ti*us**(-2) + 6*tg*ug*ui*us**(-2) - 4*tg*ug*us**(-1)
     +     + 2*tg*ug**2*us**(-2) - 2*tg**2*ug*us**(-2) + 4*tg**2*
     +    us**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4*s*ug*us**(-2) - 4*tg*ug*
     +    us**(-2) - 4*ug*ui*us**(-2) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 2*s*tg*ug*us**(-2) - 8
     +    *s*tg*us**(-1) - 2*s*ug*ti*us**(-2) - 6*s*ug*ui*us**(-2) + 4*
     +    s*ug*us**(-1) - 2*s*ug**2*us**(-2) - 4*s**2*us**(-1) - 2*tg*
     +    ug*ti*us**(-2) - 6*tg*ug*ui*us**(-2) + 4*tg*ug*us**(-1) - 2*
     +    tg*ug**2*us**(-2) + 2*tg**2*ug*us**(-2) - 4*tg**2*us**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*s*ug*us**(-2) + 4*tg*ug*us**(-2)
     +     + 4*ug*ui*us**(-2) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 2*s*us**(-1) + 2*tg*
     +    us**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ts**(-1)*us**(-1) + 2*tg*
     +    ts**(-1)*us**(-1) - 4*us**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 2*s*us**(-1) - 2*tg
     +    *us**(-1) )
      qqH = qqH + NGANG(1,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*us**(-1) )
      qqH = qqH + NGANG(1,2,-1,-1)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 4*s*tg*ug*ts**(-1)
     +     + 4*s*tg*ui*ts**(-1) + 4*s*tg**2*ts**(-1) + 4*s*ug*ui*
     +    ts**(-1) + 4*s**2*ug*ts**(-1) + 4*s**2*ui*ts**(-1) + 4*s**3*
     +    ts**(-1) )
      qqH = qqH + NGANG(1,2,-1,-1)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 4*s*tg*ug*us**(-1)
     +     + 4*s*tg*ui*us**(-1) + 4*s*tg**2*us**(-1) + 4*s*ug*ui*
     +    us**(-1) + 8*s**2*tg*us**(-1) - 4*s**2*ug*us**(-1) + 4*s**2*
     +    ui*us**(-1) + 4*s**3*us**(-1) )
      qqH = qqH + NGANG(1,2,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*s4imus**(-1)*four**(-1) * ( 4*s**2 )
      qqH = qqH + NGANG(1,2,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*s**2*ts**(-1)*us**(-1) )
      qqH = qqH + NGANG(1,3,-1,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*tg**2*ti*ts**(-2) )
      qqH = qqH + NGANG(1,3,-1,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 4*s*tg*ti*us**(-1)
     +     - 4*s*tg**2*us**(-1) - 4*s**2*tg*us**(-1) - 4*tg*ug*ui*
     +    us**(-1) - 4*tg**2*ti*us**(-1) )
      qqH = qqH + NGANG(1,3,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 2*s*tg*ts**(-1) - 2*tg
     +    *ug*ts**(-1) - 2*tg*ti*ts**(-1) )
      qqH = qqH + NGANG(1,3,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*s*tg*ts**(-1)*us**(-1) + 2*tg*
     +    ug*ts**(-1)*us**(-1) + 2*tg*ti*ts**(-1)*us**(-1) )
      qqH = qqH + NGANG(1,9,-1,-2)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*tg*ug*ts**(-1)*ms**2 - 2*tg*ug*
     +    ts**(-1)*mi**2 - 2*tg*ui*ts**(-1)*ms**2 + 2*tg*ui*ts**(-1)*
     +    mi**2 - 4*tg*ts**(-1)*ms**2*mi**2 + 2*tg*ts**(-1)*ms**4 + 2*
     +    tg*ts**(-1)*mi**4 - 6*tg*ms**2 + 6*tg*mi**2 + 6*tg**2*
     +    ts**(-1)*ms**2 - 6*tg**2*ts**(-1)*mi**2 + 4*ts*ms**2 - 4*ts*
     +    mi**2 )
      qqH = qqH + NGANG(1,9,-1,-2)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*tg*ms**2 - 8*tg*mi**2 - 8*tg**2*
     +    ts**(-1)*ms**2 + 8*tg**2*ts**(-1)*mi**2 - 4*ts*ms**2 + 4*ts*
     +    mi**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 4*s*tg*ui*ts**(-1) + 8
     +    *s*tg*ts**(-1)*ms**2 - 8*s*tg*ts**(-1)*mi**2 + 4*s*tg**2*
     +    ts**(-1) - 4*s*ms**2 + 4*s*mi**2 + 4*s**2*tg*ts**(-1) )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4*s*tg*ts**(-1) + 2*tg*ug*
     +    ts**(-1) - 2*tg*ui*ts**(-1) + 6*tg*ts**(-1)*ms**2 - 6*tg*
     +    ts**(-1)*mi**2 - 6*tg + 2*tg**2*ts**(-1) + 4*ts - 4*ms**2 + 4
     +    *mi**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*tg - 8*tg**2*ts**(-2)*ms**2 + 8*
     +    tg**2*ts**(-2)*mi**2 - 8*tg**2*ts**(-1) - 4*ts + 4*ms**2 - 4*
     +    mi**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*s4imus**(-1)*four**(-1) * ( 4*s**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - s*tg*us**(-1) + 3*s
     +    *ug*us**(-1) + s*ui*us**(-1) + s*ts*us**(-1) + 3*s*us**(-1)*
     +    ms**2 - 3*s*us**(-1)*mi**2 - 4*s + 2*tg*ug*us**(-1) + 2*tg*
     +    us**(-1)*ms**2 - 2*tg*us**(-1)*mi**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - s*tg*ts**(-1)*us**(-1) - s*ug*
     +    ts**(-1)*us**(-1) + s*ui*ts**(-1)*us**(-1) - s*ts**(-1)*
     +    us**(-1)*ms**2 + s*ts**(-1)*us**(-1)*mi**2 - 3*s*us**(-1) + 2
     +    *tg*ug*ts**(-1)*us**(-1) + 2*tg*ts**(-1)*us**(-1)*ms**2 - 2*
     +    tg*ts**(-1)*us**(-1)*mi**2 - 4*ug*us**(-1) - 4*us**(-1)*ms**2
     +     + 4*us**(-1)*mi**2 )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 2*s*tg*ts**(-1) + s
     +    *tg*us**(-1) - 3*s*ug*us**(-1) - s*ui*us**(-1) - s*ts*
     +    us**(-1) - 3*s*us**(-1)*ms**2 + 3*s*us**(-1)*mi**2 + 2*s + 2*
     +    tg*ug*ts**(-1) - 2*tg*ug*us**(-1) + 2*tg*ts**(-1)*ms**2 - 2*
     +    tg*ts**(-1)*mi**2 - 2*tg*us**(-1)*ms**2 + 2*tg*us**(-1)*mi**2
     +     )
      qqH = qqH + NGANG(1,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*tg*ts**(-1)*us**(-1) + 2*s*
     +    us**(-1) - 2*tg*ug*ts**(-1)*us**(-1) - 2*tg*ts**(-1)*us**(-1)
     +    *ms**2 + 2*tg*ts**(-1)*us**(-1)*mi**2 + 4*ug*us**(-1) + 4*
     +    us**(-1)*ms**2 - 4*us**(-1)*mi**2 )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 8*s*tg*ug*ts**(-2)
     +     + 8*s*tg*ui*ts**(-2) - 4*s*tg*ts**(-1) + 8*s*tg**2*ts**(-2)
     +     + 8*s*ug*ts**(-1) + 4*s**2*ts**(-1) + 8*tg*ug*ui*ts**(-2) - 
     +    4*tg*ug*ts**(-1) - 8*tg*ug**2*ts**(-2) + 8*tg**2*ug*ts**(-2)
     +     + 4*ug**2*ts**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4*s*tg*ts**(-2) - 4*tg*ui*
     +    ts**(-2) - 4*tg**2*ts**(-2) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 8*s*tg*ug*ts**(-2) - 8
     +    *s*tg*ui*ts**(-2) + 4*s*tg*ts**(-1) - 8*s*tg**2*ts**(-2) - 8*
     +    s*ug*ts**(-1) - 4*s**2*ts**(-1) - 8*tg*ug*ui*ts**(-2) + 4*tg*
     +    ug*ts**(-1) + 8*tg*ug**2*ts**(-2) - 8*tg**2*ug*ts**(-2) - 4*
     +    ug**2*ts**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*s*tg*ts**(-2) + 4*tg*ui*ts**(-2)
     +     + 4*tg**2*ts**(-2) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 4*s*tg*us**(-1) + 4*s*
     +    us**(-1)*ms**2 - 4*s*us**(-1)*mi**2 + 4*s**2*us**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4 + 4*ug*us**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4 )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 2*s*ts**(-1) + 2*ug*
     +    ts**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ts**(-1)*us**(-1) + 2*ug*
     +    ts**(-1)*us**(-1) - 4*ts**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 2*s*ts**(-1) - 2*ug
     +    *ts**(-1) )
      qqH = qqH + NGANG(2,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*ts**(-1) )
      qqH = qqH + NGANG(2,3,-1,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 4*s*ug*ui*ts**(-1)
     +     - 4*s*ug**2*ts**(-1) - 4*s**2*ug*ts**(-1) - 4*tg*ug*ui*
     +    ts**(-1) + 4*tg*ug**2*ts**(-1) - 4*tg**2*ug*ts**(-1) - 4*
     +    ug**2*ui*ts**(-1) )
      qqH = qqH + NGANG(2,3,-1,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*ug**2*ui*us**(-2) )
      qqH = qqH + NGANG(2,3,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 2*s*ug*us**(-1) - 2*tg
     +    *ug*us**(-1) - 2*ug*ui*us**(-1) )
      qqH = qqH + NGANG(2,3,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*s*ug*ts**(-1)*us**(-1) + 2*tg*
     +    ug*ts**(-1)*us**(-1) + 2*ug*ui*ts**(-1)*us**(-1) )
      qqH = qqH + NGANG(2,8,-1,-2)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*ug*ui*us**(-1)*ms**2 + 2*ug*ui
     +    *us**(-1)*mi**2 - 4*ug*us**(-1)*ms**2*mi**2 + 2*ug*us**(-1)*
     +    ms**4 + 2*ug*us**(-1)*mi**4 - 6*ug*ms**2 + 6*ug*mi**2 + 8*
     +    ug**2*us**(-1)*ms**2 - 8*ug**2*us**(-1)*mi**2 + 4*us*ms**2 - 
     +    4*us*mi**2 )
      qqH = qqH + NGANG(2,8,-1,-2)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*ug*ms**2 - 8*ug*mi**2 - 8*ug**2*
     +    us**(-1)*ms**2 + 8*ug**2*us**(-1)*mi**2 - 4*us*ms**2 + 4*us*
     +    mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 4*s*tg*ug*us**(-1) + 4
     +    *s*ug*ui*us**(-1) + 8*s*ug*us**(-1)*ms**2 - 8*s*ug*us**(-1)*
     +    mi**2 - 4*s*ms**2 + 4*s*mi**2 + 4*s**2*ug*us**(-1) )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4*s*ug*us**(-1) - 2*ug*ui*
     +    us**(-1) + 6*ug*us**(-1)*ms**2 - 6*ug*us**(-1)*mi**2 - 6*ug
     +     + 4*ug**2*us**(-1) + 4*us - 4*ms**2 + 4*mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*ug - 8*ug**2*us**(-2)*ms**2 + 8*
     +    ug**2*us**(-2)*mi**2 - 8*ug**2*us**(-1) - 4*us + 4*ms**2 - 4*
     +    mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*s4imus**(-1)*four**(-1) * ( 4*s**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 4*s*tg*ts**(-1) - 2*s*
     +    ug*ts**(-1) + s*ui*ts**(-1) + s*ts**(-1)*us + 3*s*ts**(-1)*
     +    ms**2 - 3*s*ts**(-1)*mi**2 - 4*s + 2*tg*ug*ts**(-1) + 2*ug*
     +    ts**(-1)*ms**2 - 2*ug*ts**(-1)*mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*s*ug*ts**(-1)*us**(-1) + s*ui*
     +    ts**(-1)*us**(-1) - s*ts**(-1)*us**(-1)*ms**2 + s*ts**(-1)*
     +    us**(-1)*mi**2 - 3*s*ts**(-1) + 2*tg*ug*ts**(-1)*us**(-1) - 4
     +    *tg*ts**(-1) + 2*ug*ts**(-1)*us**(-1)*ms**2 - 2*ug*ts**(-1)*
     +    us**(-1)*mi**2 - 4*ts**(-1)*ms**2 + 4*ts**(-1)*mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 4*s*tg*ts**(-1) + 2
     +    *s*ug*ts**(-1) - 2*s*ug*us**(-1) - s*ui*ts**(-1) - s*ts**(-1)
     +    *us - 3*s*ts**(-1)*ms**2 + 3*s*ts**(-1)*mi**2 + 2*s - 2*tg*ug
     +    *ts**(-1) + 2*tg*ug*us**(-1) - 2*ug*ts**(-1)*ms**2 + 2*ug*
     +    ts**(-1)*mi**2 + 2*ug*us**(-1)*ms**2 - 2*ug*us**(-1)*mi**2 )
      qqH = qqH + NGANG(2,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ug*ts**(-1)*us**(-1) + 2*s*
     +    ts**(-1) - 2*tg*ug*ts**(-1)*us**(-1) + 4*tg*ts**(-1) - 2*ug*
     +    ts**(-1)*us**(-1)*ms**2 + 2*ug*ts**(-1)*us**(-1)*mi**2 + 4*
     +    ts**(-1)*ms**2 - 4*ts**(-1)*mi**2 )
      qqH = qqH + NGANG(3,0,-2,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 16*mg**2 )
      qqH = qqH + NGANG(3,0,-2,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 16*mg**2 )
      qqH = qqH + NGANG(3,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 4*s*ug*ts**(-1) - 4
     +    *ug*ts**(-1)*ms**2 + 4*ug*ts**(-1)*mi**2 - 4*ug**2*ts**(-1) )
      qqH = qqH + NGANG(3,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8 - 4*tg*ts**(-1) + 2*ug*ts**(-1)
     +     - 2*ui*ts**(-1) - 2*ts**(-1)*mi**2 - 6*ts**(-1)*mg**2 )
      qqH = qqH + NGANG(3,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 4*s*tg*us**(-1) - 4
     +    *tg*us**(-1)*ms**2 + 4*tg*us**(-1)*mi**2 - 4*tg**2*us**(-1) )
      qqH = qqH + NGANG(3,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8 - 2*ug*us**(-1) - 2*ui*us**(-1)
     +     - 2*us**(-1)*mi**2 - 6*us**(-1)*mg**2 )
      qqH = qqH + NGANG(3,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*ts**(-1) + 4*us**(-1) )
      qqH = qqH + NGANG(3,8,-2,-2)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 16*ms**2*mi**2*mg**2 + 16*ms**2*
     +    mg**4 - 16*ms**4*mg**2 - 16*mi**2*mg**4 )
      qqH = qqH + NGANG(3,8,-2,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 32*ms**2*mg**2 + 16*mi**2*mg**2
     +     + 16*mg**4 )
      qqH = qqH + NGANG(3,8,-2,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-1)*four**(-1) * (  - 16*s*mg**2 )
      qqH = qqH + NGANG(3,8,-1,-2)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 6*ug*us**(-1)*ms**2*mi**2 - 2*ug*
     +    us**(-1)*ms**2*mg**2 - 6*ug*us**(-1)*ms**4 + 2*ug*us**(-1)*
     +    mi**2*mg**2 - 8*ug**2*us**(-1)*ms**2 + 8*ug**2*us**(-1)*mi**2
     +     + 2*ui*us**(-1)*ms**2*mi**2 + 2*ui*us**(-1)*ms**2*mg**2 - 2*
     +    ui*us**(-1)*ms**4 - 2*ui*us**(-1)*mi**2*mg**2 + 8*us**(-1)*
     +    ms**2*mi**2*mg**2 + 2*us**(-1)*ms**2*mi**4 + 6*us**(-1)*ms**2
     +    *mg**4 - 2*us**(-1)*ms**4*mi**2 - 6*us**(-1)*ms**4*mg**2 - 6*
     +    us**(-1)*mi**2*mg**4 - 2*us**(-1)*mi**4*mg**2 - 8*ms**2*mi**2
     +     + 8*ms**4 )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 4*s*tg*ug*us**(-1)
     +     - 4*tg*ug*ui*us**(-1) - 4*tg*ug*us**(-1)*ms**2 + 4*tg*ug*
     +    us**(-1)*mi**2 + 4*tg*us**(-1)*ms**2*mi**2 + 4*tg*us**(-1)*
     +    ms**2*mg**2 - 4*tg*us**(-1)*ms**4 - 4*tg*us**(-1)*mi**2*mg**2
     +     - 4*tg**2*ug*us**(-1) )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*tg*ug*us**(-1) + 4*ug*ui*us**(-1)
     +     - 8*ug*us**(-1)*ms**2 + 6*ug*us**(-1)*mi**2 + 2*ug*us**(-1)*
     +    mg**2 - 8*ug**2*us**(-2)*ms**2 + 8*ug**2*us**(-2)*mi**2 - 8*
     +    ug**2*us**(-1) - 4*ui*us**(-1)*ms**2 + 2*ui*us**(-1)*mi**2 + 
     +    2*ui*us**(-1)*mg**2 - 4*us**(-1)*ms**2*mi**2 - 12*us**(-1)*
     +    ms**2*mg**2 + 8*us**(-1)*mi**2*mg**2 + 2*us**(-1)*mi**4 + 6*
     +    us**(-1)*mg**4 + 16*ms**2 - 8*mi**2 )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 4*s*tg*ts**(-1) - 2*tg
     +    *ug*ts**(-1) + 2*tg*ts**(-1)*ms**2 - 2*tg*ts**(-1)*mg**2 )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-2)*four**(-1) * (  - 16*s*mg**2 )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-1)*four**(-1) * (  - 8*s*tg*ts**(-1) - s*
     +    ug*ts**(-1) - 3*s*ug*us**(-1) - s*ui*ts**(-1) - s*ui*us**(-1)
     +     - s*ts**(-1)*mi**2 - 7*s*ts**(-1)*mg**2 - s*us**(-1)*mi**2
     +     - 7*s*us**(-1)*mg**2 - 8*s + 2*tg*ug*ts**(-1) + 2*tg*ug*
     +    us**(-1) - 4*tg*ts**(-1)*ms**2 + 4*tg*ts**(-1)*mg**2 + 2*tg*
     +    us**(-1)*ms**2 - 2*tg*us**(-1)*mg**2 - 2*tg - 2*ug*ts**(-1)*
     +    ms**2 + 2*ug*ts**(-1)*mi**2 + 4*ug*us**(-1)*ms**2 - 4*ug*
     +    us**(-1)*mi**2 - 2*ug )
      qqH = qqH + NGANG(3,8,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ug*ts**(-1)*us**(-1) - 2*tg*ug*
     +    ts**(-1)*us**(-1) + 6*tg*ts**(-1) - 2*ug*ts**(-1)*us**(-1)*
     +    ms**2 + 2*ug*ts**(-1)*us**(-1)*mi**2 + 4*ts**(-1)*ms**2 - 4*
     +    ts**(-1)*mi**2 )
      qqH = qqH + NGANG(3,9,-2,-2)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 16*ms**2*mi**2*mg**2 + 16*ms**2*
     +    mg**4 - 16*ms**4*mg**2 - 16*mi**2*mg**4 )
      qqH = qqH + NGANG(3,9,-2,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 32*ms**2*mg**2 + 16*mi**2*mg**2
     +     + 16*mg**4 )
      qqH = qqH + NGANG(3,9,-2,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-1)*four**(-1) * (  - 16*s*mg**2 )
      qqH = qqH + NGANG(3,9,-1,-2)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 8*tg*ts**(-1)*ms**2*mi**2 - 8*tg*
     +    ts**(-1)*ms**4 - 8*tg**2*ts**(-1)*ms**2 + 8*tg**2*ts**(-1)*
     +    mi**2 - 2*ug*ts**(-1)*ms**2*mi**2 - 2*ug*ts**(-1)*ms**2*mg**2
     +     + 2*ug*ts**(-1)*ms**4 + 2*ug*ts**(-1)*mi**2*mg**2 + 2*ui*
     +    ts**(-1)*ms**2*mi**2 + 2*ui*ts**(-1)*ms**2*mg**2 - 2*ui*
     +    ts**(-1)*ms**4 - 2*ui*ts**(-1)*mi**2*mg**2 + 8*ts**(-1)*ms**2
     +    *mi**2*mg**2 + 2*ts**(-1)*ms**2*mi**4 + 6*ts**(-1)*ms**2*
     +    mg**4 - 2*ts**(-1)*ms**4*mi**2 - 6*ts**(-1)*ms**4*mg**2 - 6*
     +    ts**(-1)*mi**2*mg**4 - 2*ts**(-1)*mi**4*mg**2 - 8*ms**2*mi**2
     +     + 8*ms**4 )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 4*s*tg*ug*ts**(-1)
     +     - 4*tg*ug*ui*ts**(-1) - 4*tg*ug*ts**(-1)*ms**2 + 4*tg*ug*
     +    ts**(-1)*mi**2 - 4*tg**2*ug*ts**(-1) + 4*ug*ts**(-1)*ms**2*
     +    mi**2 + 4*ug*ts**(-1)*ms**2*mg**2 - 4*ug*ts**(-1)*ms**4 - 4*
     +    ug*ts**(-1)*mi**2*mg**2 )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4*tg*ui*ts**(-1) - 12*tg*ts**(-1)*
     +    ms**2 + 8*tg*ts**(-1)*mi**2 + 4*tg*ts**(-1)*mg**2 - 8*tg**2*
     +    ts**(-2)*ms**2 + 8*tg**2*ts**(-2)*mi**2 - 4*tg**2*ts**(-1) + 
     +    4*ug*ts**(-1)*ms**2 - 2*ug*ts**(-1)*mi**2 - 2*ug*ts**(-1)*
     +    mg**2 - 4*ui*ts**(-1)*ms**2 + 2*ui*ts**(-1)*mi**2 + 2*ui*
     +    ts**(-1)*mg**2 - 4*ts**(-1)*ms**2*mi**2 - 12*ts**(-1)*ms**2*
     +    mg**2 + 8*ts**(-1)*mi**2*mg**2 + 2*ts**(-1)*mi**4 + 6*
     +    ts**(-1)*mg**4 + 16*ms**2 - 8*mi**2 )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 4*s*ug*us**(-1) - 2*tg
     +    *ug*us**(-1) + 2*ug*us**(-1)*ms**2 - 2*ug*us**(-1)*mg**2 )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-2)*four**(-1) * (  - 16*s*mg**2 )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-1)*four**(-1) * (  - 4*s*tg*ts**(-1) - 2*s
     +    *tg*us**(-1) + s*ug*ts**(-1) - 7*s*ug*us**(-1) - s*ui*
     +    ts**(-1) - s*ui*us**(-1) - s*ts**(-1)*mi**2 - 7*s*ts**(-1)*
     +    mg**2 - s*us**(-1)*mi**2 - 7*s*us**(-1)*mg**2 - 8*s + 2*tg*ug
     +    *ts**(-1) + 2*tg*ug*us**(-1) + 4*tg*ts**(-1)*ms**2 - 4*tg*
     +    ts**(-1)*mi**2 - 2*tg*us**(-1)*ms**2 + 2*tg*us**(-1)*mi**2 - 
     +    2*tg + 2*ug*ts**(-1)*ms**2 - 2*ug*ts**(-1)*mg**2 - 4*ug*
     +    us**(-1)*ms**2 + 4*ug*us**(-1)*mg**2 - 2*ug )
      qqH = qqH + NGANG(3,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*tg*ts**(-1)*us**(-1) - 2*tg*ug*
     +    ts**(-1)*us**(-1) - 2*tg*ts**(-1)*us**(-1)*ms**2 + 2*tg*
     +    ts**(-1)*us**(-1)*mi**2 + 6*ug*us**(-1) + 4*us**(-1)*ms**2 - 
     +    4*us**(-1)*mi**2 )
      qqH = qqH + NGANG(8,0,-2,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*ug*ui*us**(-2)*ms**2 - 2*ug*ui*
     +    us**(-2)*mi**2 - 2*ug*us**(-2)*ms**2*mi**2 + 2*ug*us**(-2)*
     +    ms**2*mg**2 + 4*ug*us**(-2)*ms**4 - 2*ug*us**(-2)*mi**2*mg**2
     +     - 2*ug*us**(-2)*mi**4 - 2*ug*us**(-1)*ms**2 + 2*ug*us**(-1)*
     +    mi**2 + 2*ug**2*us**(-2)*ms**2 - 2*ug**2*us**(-2)*mi**2 - 2*
     +    ui*us**(-1)*ms**2 + 2*ui*us**(-1)*mi**2 - 2*us**(-1)*ms**2*
     +    mi**2 - 6*us**(-1)*ms**2*mg**2 + 6*us**(-1)*mi**2*mg**2 + 2*
     +    us**(-1)*mi**4 + 4*ms**2 - 4*mi**2 )
      qqH = qqH + NGANG(8,0,-2,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*ug*ui*us**(-2)*ms**2 + 2*ug*ui
     +    *us**(-2)*mi**2 + 2*ug*us**(-2)*ms**2*mi**2 - 2*ug*us**(-2)*
     +    ms**2*mg**2 - 4*ug*us**(-2)*ms**4 + 2*ug*us**(-2)*mi**2*mg**2
     +     + 2*ug*us**(-2)*mi**4 + 2*ug*us**(-1)*ms**2 - 2*ug*us**(-1)*
     +    mi**2 - 2*ug**2*us**(-2)*ms**2 + 2*ug**2*us**(-2)*mi**2 + 2*
     +    ui*us**(-1)*ms**2 - 2*ui*us**(-1)*mi**2 + 2*us**(-1)*ms**2*
     +    mi**2 + 6*us**(-1)*ms**2*mg**2 - 6*us**(-1)*mi**2*mg**2 - 2*
     +    us**(-1)*mi**4 - 4*ms**2 + 4*mi**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 4*s*tg*ug*us**(-2) + 2
     +    *s*ug*ui*us**(-2) - 4*s*us**(-1)*ms**2 + 4*s*us**(-1)*mi**2
     +     + 2*s**2*ug*us**(-2) + 2*tg*ug*ui*us**(-2) - 4*tg*us**(-1)*
     +    ms**2 + 4*tg*us**(-1)*mi**2 + 2*tg**2*ug*us**(-2) - 4*ug*ui*
     +    us**(-2)*ms**2 + 4*ug*ui*us**(-2)*mi**2 + 6*ug*us**(-2)*ms**2
     +    *mi**2 - 2*ug*us**(-2)*ms**2*mg**2 - 2*ug*us**(-2)*ms**4 + 2*
     +    ug*us**(-2)*mi**2*mg**2 - 4*ug*us**(-2)*mi**4 + 4*ug*us**(-1)
     +    *ms**2 - 4*ug*us**(-1)*mi**2 - 2*ug**2*us**(-2)*ms**2 + 2*
     +    ug**2*us**(-2)*mi**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4 - 2*s*ug*us**(-2) - 2*tg*ug*
     +    us**(-2) + 6*ug*us**(-2)*ms**2 - 8*ug*us**(-2)*mi**2 + 2*ug*
     +    us**(-2)*mg**2 - 6*ug*us**(-1) + 2*ug**2*us**(-2) - 2*ui*
     +    us**(-1) - 2*us**(-1)*mi**2 - 6*us**(-1)*mg**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 4*s*tg*ug*us**(-2)
     +     - 2*s*ug*ui*us**(-2) + 4*s*us**(-1)*ms**2 - 4*s*us**(-1)*
     +    mi**2 - 2*s**2*ug*us**(-2) - 2*tg*ug*ui*us**(-2) + 4*tg*
     +    us**(-1)*ms**2 - 4*tg*us**(-1)*mi**2 - 2*tg**2*ug*us**(-2) + 
     +    4*ug*ui*us**(-2)*ms**2 - 4*ug*ui*us**(-2)*mi**2 - 6*ug*
     +    us**(-2)*ms**2*mi**2 + 2*ug*us**(-2)*ms**2*mg**2 + 2*ug*
     +    us**(-2)*ms**4 - 2*ug*us**(-2)*mi**2*mg**2 + 4*ug*us**(-2)*
     +    mi**4 - 4*ug*us**(-1)*ms**2 + 4*ug*us**(-1)*mi**2 + 2*ug**2*
     +    us**(-2)*ms**2 - 2*ug**2*us**(-2)*mi**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc2*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4 + 2*s*ug*us**(-2) + 2*tg*ug*
     +    us**(-2) - 6*ug*us**(-2)*ms**2 + 8*ug*us**(-2)*mi**2 - 2*ug*
     +    us**(-2)*mg**2 + 6*ug*us**(-1) - 2*ug**2*us**(-2) + 2*ui*
     +    us**(-1) + 2*us**(-1)*mi**2 + 6*us**(-1)*mg**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( s*us**(-1) + tg*
     +    us**(-1) - us**(-1)*ms**2 + us**(-1)*mi**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ts**(-1)*us**(-1) + 2*tg*
     +    ts**(-1)*us**(-1) + us**(-1) )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - s*us**(-1) - tg*
     +    us**(-1) + us**(-1)*ms**2 - us**(-1)*mi**2 )
      qqH = qqH + NGANG(8,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - us**(-1) )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*s4imus**(-1)*four**(-1) * ( 4*s**2 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 3*s*tg*ts**(-1) - s*ug
     +    *ts**(-1) + s*ui*ts**(-1) + s*ts**(-1)*ms**2 - 2*s*ts**(-1)*
     +    mi**2 + s*ts**(-1)*mg**2 - 4*s + tg*ug*ts**(-1) - tg*ts**(-1)
     +    *ms**2 + tg*ts**(-1)*mi**2 + ug*ts**(-1)*ms**2 - ug*ts**(-1)*
     +    mi**2 + 2*ts**(-1)*ms**2*mi**2 - ts**(-1)*ms**4 - ts**(-1)*
     +    mi**4 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * ( 2*s*ug*us**(-1) + s*ui
     +    *us**(-1) + s*us**(-1)*ms**2 - 2*s*us**(-1)*mi**2 + s*
     +    us**(-1)*mg**2 - 4*s + tg*ug*us**(-1) + tg*us**(-1)*ms**2 - 
     +    tg*us**(-1)*mi**2 - ug*us**(-1)*ms**2 + ug*us**(-1)*mi**2 + 2
     +    *us**(-1)*ms**2*mi**2 - us**(-1)*ms**4 - us**(-1)*mi**4 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*tg*ts**(-1)*us**(-1) + 2*s*ui*
     +    ts**(-1)*us**(-1) + 4*s*ts**(-1)*us**(-1)*ms**2 + 2*s*
     +    ts**(-1)*us**(-1)*mi**2 + 2*s*ts**(-1)*us**(-1)*mg**2 + 4*
     +    s**2*ts**(-1)*us**(-1) - tg*ts**(-1) - ug*us**(-1) + ts**(-1)
     +    *ms**2 - ts**(-1)*mi**2 + us**(-1)*ms**2 - us**(-1)*mi**2 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 3*s*tg*ts**(-1) + s
     +    *ug*ts**(-1) - 2*s*ug*us**(-1) - s*ui*ts**(-1) - s*ts**(-1)*
     +    ms**2 + 2*s*ts**(-1)*mi**2 - s*ts**(-1)*mg**2 + 2*s - tg*ug*
     +    ts**(-1) + 2*tg*ug*us**(-1) + tg*ts**(-1)*ms**2 - tg*ts**(-1)
     +    *mi**2 - ug*ts**(-1)*ms**2 + ug*ts**(-1)*mi**2 + 2*ug*
     +    us**(-1)*ms**2 - 2*ug*us**(-1)*mi**2 - 2*ts**(-1)*ms**2*mi**2
     +     + ts**(-1)*ms**4 + ts**(-1)*mi**4 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imus**(-1)*four**(-1) * (  - 2*s*tg*ts**(-1) - 2
     +    *s*ug*us**(-1) - s*ui*us**(-1) - s*us**(-1)*ms**2 + 2*s*
     +    us**(-1)*mi**2 - s*us**(-1)*mg**2 + 2*s + 2*tg*ug*ts**(-1) - 
     +    tg*ug*us**(-1) + 2*tg*ts**(-1)*ms**2 - 2*tg*ts**(-1)*mi**2 - 
     +    tg*us**(-1)*ms**2 + tg*us**(-1)*mi**2 + ug*us**(-1)*ms**2 - 
     +    ug*us**(-1)*mi**2 - 2*us**(-1)*ms**2*mi**2 + us**(-1)*ms**4
     +     + us**(-1)*mi**4 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-2)*four**(-1) * ( 16*s*mg**2 )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*sfac**(-1)*four**(-1) * ( 4*s*tg*ts**(-1) + 2*s*tg
     +    *us**(-1) + s*ug*ts**(-1) + 3*s*ug*us**(-1) + s*ui*ts**(-1)
     +     + s*ui*us**(-1) + 6*s*ts**(-1)*ms**2 - 2*s*ts**(-1)*mi**2 + 
     +    4*s*ts**(-1)*mg**2 + 6*s*us**(-1)*ms**2 - 2*s*us**(-1)*mi**2
     +     + 4*s*us**(-1)*mg**2 + 8*s + 3*s**2*ts**(-1) + 3*s**2*
     +    us**(-1) - 2*tg*ug*ts**(-1) - 2*tg*ug*us**(-1) - 4*tg*
     +    ts**(-1)*ms**2 + 4*tg*ts**(-1)*mi**2 + 2*tg*us**(-1)*ms**2 - 
     +    2*tg*us**(-1)*mi**2 + 2*tg + 2*ug*ts**(-1)*ms**2 - 2*ug*
     +    ts**(-1)*mi**2 - 4*ug*us**(-1)*ms**2 + 4*ug*us**(-1)*mi**2 + 
     +    2*ug )
      qqH = qqH + NGANG(8,9,-1,-1)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( s*ts**(-1) + s*us**(-1) + tg*
     +    ts**(-1) + ug*us**(-1) - ts**(-1)*ms**2 + ts**(-1)*mi**2 - 
     +    us**(-1)*ms**2 + us**(-1)*mi**2 )
      qqH = qqH + NGANG(9,0,-2,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 2*tg*ug*ts**(-2)*ms**2 + 2*tg*ug
     +    *ts**(-2)*mi**2 + 2*tg*ui*ts**(-2)*ms**2 - 2*tg*ui*ts**(-2)*
     +    mi**2 - 2*tg*ts**(-2)*ms**2*mi**2 + 2*tg*ts**(-2)*ms**2*mg**2
     +     + 4*tg*ts**(-2)*ms**4 - 2*tg*ts**(-2)*mi**2*mg**2 - 2*tg*
     +    ts**(-2)*mi**4 - 4*tg*ts**(-1)*ms**2 + 4*tg*ts**(-1)*mi**2 + 
     +    4*tg**2*ts**(-2)*ms**2 - 4*tg**2*ts**(-2)*mi**2 + 2*ug*
     +    ts**(-1)*ms**2 - 2*ug*ts**(-1)*mi**2 - 2*ui*ts**(-1)*ms**2 + 
     +    2*ui*ts**(-1)*mi**2 - 2*ts**(-1)*ms**2*mi**2 - 6*ts**(-1)*
     +    ms**2*mg**2 + 6*ts**(-1)*mi**2*mg**2 + 2*ts**(-1)*mi**4 + 4*
     +    ms**2 - 4*mi**2 )
      qqH = qqH + NGANG(9,0,-2,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*tg*ug*ts**(-2)*ms**2 - 2*tg*ug*
     +    ts**(-2)*mi**2 - 2*tg*ui*ts**(-2)*ms**2 + 2*tg*ui*ts**(-2)*
     +    mi**2 + 2*tg*ts**(-2)*ms**2*mi**2 - 2*tg*ts**(-2)*ms**2*mg**2
     +     - 4*tg*ts**(-2)*ms**4 + 2*tg*ts**(-2)*mi**2*mg**2 + 2*tg*
     +    ts**(-2)*mi**4 + 4*tg*ts**(-1)*ms**2 - 4*tg*ts**(-1)*mi**2 - 
     +    4*tg**2*ts**(-2)*ms**2 + 4*tg**2*ts**(-2)*mi**2 - 2*ug*
     +    ts**(-1)*ms**2 + 2*ug*ts**(-1)*mi**2 + 2*ui*ts**(-1)*ms**2 - 
     +    2*ui*ts**(-1)*mi**2 + 2*ts**(-1)*ms**2*mi**2 + 6*ts**(-1)*
     +    ms**2*mg**2 - 6*ts**(-1)*mi**2*mg**2 - 2*ts**(-1)*mi**4 - 4*
     +    ms**2 + 4*mi**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( 2*s*tg*ug*ts**(-2) + 2
     +    *s*tg*ui*ts**(-2) + 2*s*tg**2*ts**(-2) - 4*s*ts**(-1)*ms**2
     +     + 4*s*ts**(-1)*mi**2 + 2*s**2*tg*ts**(-2) + 2*tg*ug*ui*
     +    ts**(-2) + 4*tg*ug*ts**(-2)*ms**2 - 4*tg*ug*ts**(-2)*mi**2 - 
     +    4*tg*ui*ts**(-2)*ms**2 + 4*tg*ui*ts**(-2)*mi**2 + 6*tg*
     +    ts**(-2)*ms**2*mi**2 - 2*tg*ts**(-2)*ms**2*mg**2 - 2*tg*
     +    ts**(-2)*ms**4 + 2*tg*ts**(-2)*mi**2*mg**2 - 4*tg*ts**(-2)*
     +    mi**4 + 4*tg*ts**(-1)*ms**2 - 4*tg*ts**(-1)*mi**2 + 2*tg**2*
     +    ug*ts**(-2) - 6*tg**2*ts**(-2)*ms**2 + 6*tg**2*ts**(-2)*mi**2
     +     - 4*ug*ts**(-1)*ms**2 + 4*ug*ts**(-1)*mi**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 4 - 2*s*tg*ts**(-2) - 2*tg*ug*
     +    ts**(-2) + 6*tg*ts**(-2)*ms**2 - 8*tg*ts**(-2)*mi**2 + 2*tg*
     +    ts**(-2)*mg**2 - 8*tg*ts**(-1) + 2*tg**2*ts**(-2) + 2*ug*
     +    ts**(-1) - 2*ui*ts**(-1) - 2*ts**(-1)*mi**2 - 6*ts**(-1)*
     +    mg**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - 2*s*tg*ug*ts**(-2)
     +     - 2*s*tg*ui*ts**(-2) - 2*s*tg**2*ts**(-2) + 4*s*ts**(-1)*
     +    ms**2 - 4*s*ts**(-1)*mi**2 - 2*s**2*tg*ts**(-2) - 2*tg*ug*ui*
     +    ts**(-2) - 4*tg*ug*ts**(-2)*ms**2 + 4*tg*ug*ts**(-2)*mi**2 + 
     +    4*tg*ui*ts**(-2)*ms**2 - 4*tg*ui*ts**(-2)*mi**2 - 6*tg*
     +    ts**(-2)*ms**2*mi**2 + 2*tg*ts**(-2)*ms**2*mg**2 + 2*tg*
     +    ts**(-2)*ms**4 - 2*tg*ts**(-2)*mi**2*mg**2 + 4*tg*ts**(-2)*
     +    mi**4 - 4*tg*ts**(-1)*ms**2 + 4*tg*ts**(-1)*mi**2 - 2*tg**2*
     +    ug*ts**(-2) + 6*tg**2*ts**(-2)*ms**2 - 6*tg**2*ts**(-2)*mi**2
     +     + 4*ug*ts**(-1)*ms**2 - 4*ug*ts**(-1)*mi**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc1*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - 4 + 2*s*tg*ts**(-2) + 2*tg*ug*
     +    ts**(-2) - 6*tg*ts**(-2)*ms**2 + 8*tg*ts**(-2)*mi**2 - 2*tg*
     +    ts**(-2)*mg**2 + 8*tg*ts**(-1) - 2*tg**2*ts**(-2) - 2*ug*
     +    ts**(-1) + 2*ui*ts**(-1) + 2*ts**(-1)*mi**2 + 6*ts**(-1)*
     +    mg**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * ( s*ts**(-1) + ug*
     +    ts**(-1) - ts**(-1)*ms**2 + ts**(-1)*mi**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc3*CA**(-2)*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * ( 2*s*ts**(-1)*us**(-1) + 2*ug*
     +    ts**(-1)*us**(-1) + ts**(-1) )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*s4imts**(-1)*four**(-1) * (  - s*ts**(-1) - ug*
     +    ts**(-1) + ts**(-1)*ms**2 - ts**(-1)*mi**2 )
      qqH = qqH + NGANG(9,0,-1,0)*e**2*gs**4*struc3*CF*Sn**2*
     + hardfac**(-1)*four**(-1) * (  - ts**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*tshat**(-2)*e**2*
     + gs**4*struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 
     +     - 32*tg**2*ui**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*tshat**(-2)*e**2*
     + gs**4*struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 
     +     - 32*tg*ui*mi**2 + 32*tg*ui*mg**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*tshat**(-1)*e**2*
     + gs**4*struc3*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 
     +     - 32*s*ui*us**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*ushat**(-2)*e**2*
     + gs**4*struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 
     +     - 32*ug**2*ti**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*ushat**(-2)*e**2*
     + gs**4*struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 
     +     - 32*ug*ti*mi**2 + 32*ug*ti*mg**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*ushat**(-1)*e**2*
     + gs**4*struc3*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 
     +     - 32*s*ti*ts**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*e**2*gs**4*struc1
     + *CA**(-1)*CF**2*pi*Sn**2*four**(-1) * ( 32*tg*ts**(-2)*mi**2 - 
     +    32*tg*ts**(-2)*mg**2 - 32*tg**2*ts**(-2) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*s4i**(-1)*e**2*gs**4*struc2
     + *CA**(-1)*CF**2*pi*Sn**2*four**(-1) * ( 32*ug*us**(-2)*mi**2 - 
     +    32*ug*us**(-2)*mg**2 - 32*ug**2*us**(-2) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-3)*four**(-1) * ( 16*
     +    tg**2*ui**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 16*tg*
     +    ui*mi**2 - 16*tg*ui*mg**2 + 16*tg**2*ui )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 16*tg*
     +    mi**2 - 16*tg*mg**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 16*s*ui
     +    *us**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 16*s*
     +    us**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-3)*four**(-1) * ( 16*
     +    ug**2*ti**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 16*ug*
     +    ti*mi**2 - 16*ug*ti*mg**2 + 16*ug**2*ti )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 16*ug*
     +    mi**2 - 16*ug*mg**2 )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 16*s*ti
     +    *ts**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 16*s*
     +    ts**(-1) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc1*CA**(-1)*
     + CF**2*pi*Sn**2*spug**(-1)*four**(-1) * (  - 16*tg*ts**(-2)*mi**2
     +     + 16*tg*ts**(-2)*mg**2 + 16*tg**2*ts**(-2) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc1*CA**(-1)*
     + CF**2*pi*Sn**2*four**(-1) * (  - 16*tg*ti**(-1)*ts**(-2)*mi**2
     +     + 16*tg*ti**(-1)*ts**(-2)*mg**2 + 16*tg**2*ti**(-1)*ts**(-2)
     +     )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc2*CA**(-1)*
     + CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * (  - 16*ug*us**(-2)*mi**2
     +     + 16*ug*us**(-2)*mg**2 + 16*ug**2*us**(-2) )
      qqH = qqH + log(1 + s4i**(-1)*mi**2)*e**2*gs**4*struc2*CA**(-1)*
     + CF**2*pi*Sn**2*four**(-1) * (  - 16*ug*ui**(-1)*us**(-2)*mi**2
     +     + 16*ug*ui**(-1)*us**(-2)*mg**2 + 16*ug**2*ui**(-1)*us**(-2)
     +     )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 32*
     +    tg**2*ui**2 )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*tshat**(-2)*e**2*gs**4*
     + struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 32*tg*
     +    ui*mi**2 - 32*tg*ui*mg**2 )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*tshat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 32*s*ui
     +    *us**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 32*
     +    ug**2*ti**2 )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*ushat**(-2)*e**2*gs**4*
     + struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 32*ug*
     +    ti*mi**2 - 32*ug*ti*mg**2 )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*ushat**(-1)*e**2*gs**4*
     + struc3*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 32*s*ti
     +    *ts**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*four**(-1) * (  - 32*tg*ts**(-2)*mi**2
     +     + 32*tg*ts**(-2)*mg**2 + 32*tg**2*ts**(-2) )
      qqH = qqH + log(s4i*mi**(-2))*s4i**(-1)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*four**(-1) * (  - 32*ug*us**(-2)*mi**2
     +     + 32*ug*us**(-2)*mg**2 + 32*ug**2*us**(-2) )
      qqH = qqH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-3)*four**(-1) * (  - 16*tg**2*
     +    ui**2 )
      qqH = qqH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * (  - 16*tg*ui*
     +    mi**2 + 16*tg*ui*mg**2 - 16*tg**2*ui )
      qqH = qqH + log(s4i*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * (  - 16*tg*mi**2
     +     + 16*tg*mg**2 )
      qqH = qqH + log(s4i*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * (  - 16*s*ui*
     +    us**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * (  - 16*s*
     +    us**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-3)*four**(-1) * (  - 16*ug**2*
     +    ti**2 )
      qqH = qqH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * (  - 16*ug*ti*
     +    mi**2 + 16*ug*ti*mg**2 - 16*ug**2*ti )
      qqH = qqH + log(s4i*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * (  - 16*ug*mi**2
     +     + 16*ug*mg**2 )
      qqH = qqH + log(s4i*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * (  - 16*s*ti*
     +    ts**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * (  - 16*s*
     +    ts**(-1) )
      qqH = qqH + log(s4i*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF**2*
     + pi*Sn**2*spug**(-1)*four**(-1) * ( 16*tg*ts**(-2)*mi**2 - 16*tg*
     +    ts**(-2)*mg**2 - 16*tg**2*ts**(-2) )
      qqH = qqH + log(s4i*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF**2*
     + pi*Sn**2*four**(-1) * ( 16*tg*ti**(-1)*ts**(-2)*mi**2 - 16*tg*
     +    ti**(-1)*ts**(-2)*mg**2 - 16*tg**2*ti**(-1)*ts**(-2) )
      qqH = qqH + log(s4i*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF**2*
     + pi*Sn**2*sptg**(-1)*four**(-1) * ( 16*ug*us**(-2)*mi**2 - 16*ug*
     +    us**(-2)*mg**2 - 16*ug**2*us**(-2) )
      qqH = qqH + log(s4i*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF**2*
     + pi*Sn**2*four**(-1) * ( 16*ug*ui**(-1)*us**(-2)*mi**2 - 16*ug*
     +    ui**(-1)*us**(-2)*mg**2 - 16*ug**2*ui**(-1)*us**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*tshat**(-2)*e**2*gs**4
     + *struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * (  - 32*
     +    tg**2*ui**2 )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*tshat**(-2)*e**2*gs**4
     + *struc1*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * (  - 32*
     +    tg*ui*mi**2 + 32*tg*ui*mg**2 )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*tshat**(-1)*e**2*gs**4
     + *struc3*CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * (  - 32*
     +    s*ui*us**(-1) )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*ushat**(-2)*e**2*gs**4
     + *struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * (  - 32*
     +    ug**2*ti**2 )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*ushat**(-2)*e**2*gs**4
     + *struc2*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * (  - 32*
     +    ug*ti*mi**2 + 32*ug*ti*mg**2 )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*ushat**(-1)*e**2*gs**4
     + *struc3*CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * (  - 32*
     +    s*ti*ts**(-1) )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*four**(-1) * ( 32*tg*ts**(-2)*mi**2 - 32
     +    *tg*ts**(-2)*mg**2 - 32*tg**2*ts**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*s4i**(-1)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*four**(-1) * ( 32*ug*us**(-2)*mi**2 - 32
     +    *ug*us**(-2)*mg**2 - 32*ug**2*us**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-3)*four**(-1) * ( 16*tg**2*ui**2
     +     )
      qqH = qqH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 16*tg*ui*mi**2
     +     - 16*tg*ui*mg**2 + 16*tg**2*ui )
      qqH = qqH + log(QF**2*mi**(-2))*tshat**(-2)*e**2*gs**4*struc1*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 16*tg*mi**2 - 
     +    16*tg*mg**2 )
      qqH = qqH + log(QF**2*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-2)*four**(-1) * ( 16*s*ui*
     +    us**(-1) )
      qqH = qqH + log(QF**2*mi**(-2))*tshat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*sptg**(-1)*four**(-1) * ( 16*s*us**(-1)
     +     )
      qqH = qqH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-3)*four**(-1) * ( 16*ug**2*ti**2
     +     )
      qqH = qqH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 16*ug*ti*mi**2
     +     - 16*ug*ti*mg**2 + 16*ug**2*ti )
      qqH = qqH + log(QF**2*mi**(-2))*ushat**(-2)*e**2*gs**4*struc2*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 16*ug*mi**2 - 
     +    16*ug*mg**2 )
      qqH = qqH + log(QF**2*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-2)*four**(-1) * ( 16*s*ti*
     +    ts**(-1) )
      qqH = qqH + log(QF**2*mi**(-2))*ushat**(-1)*e**2*gs**4*struc3*
     + CA**(-1)*CF**2*pi*Sn**2*spug**(-1)*four**(-1) * ( 16*s*ts**(-1)
     +     )
      qqH = qqH + log(QF**2*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF**2
     + *pi*Sn**2*spug**(-1)*four**(-1) * (  - 16*tg*ts**(-2)*mi**2 + 16
     +    *tg*ts**(-2)*mg**2 + 16*tg**2*ts**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*e**2*gs**4*struc1*CA**(-1)*CF**2
     + *pi*Sn**2*four**(-1) * (  - 16*tg*ti**(-1)*ts**(-2)*mi**2 + 16*
     +    tg*ti**(-1)*ts**(-2)*mg**2 + 16*tg**2*ti**(-1)*ts**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF**2
     + *pi*Sn**2*sptg**(-1)*four**(-1) * (  - 16*ug*us**(-2)*mi**2 + 16
     +    *ug*us**(-2)*mg**2 + 16*ug**2*us**(-2) )
      qqH = qqH + log(QF**2*mi**(-2))*e**2*gs**4*struc2*CA**(-1)*CF**2
     + *pi*Sn**2*four**(-1) * (  - 16*ug*ui**(-1)*us**(-2)*mi**2 + 16*
     +    ug*ui**(-1)*us**(-2)*mg**2 + 16*ug**2*ui**(-1)*us**(-2) )

c               the prefactor for the scaling functions 
      NG_QBH = qqH * (abs(mi)+mg)**2/4.D0

      end






































