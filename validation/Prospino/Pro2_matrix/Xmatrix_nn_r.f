cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  all real qb scaling functions : 
c
c    NN_QBF1(iq,iout,massin,Cs,Ct,Cl,Cr)
c    NN_QBF2(iq,iout,massin,Cs,Ct,Cl,Cr)
c    NN_QBH(massin,mkraemer)
c     using the subroutines : qqd132_nn(massin,mkraemer)
c                             qqd321_nn(massin,mkraemer)
c    
c  additional subroutine Born(massin,mkraemer) as a check 
c
c  common blocks rewritten as an array :
c   mkraemer(99):
c
c    mkraemer(1)  = mg1 
c    mkraemer(2)  = mg2 
c    mkraemer(3)  = mg 
c    mkraemer(4)  = ms 
c    mkraemer(11) = mx 
c    mkraemer(12) = V1  
c    mkraemer(13) = V2  
c    mkraemer(14) = A2  
c    mkraemer(15) = V1w  
c    mkraemer(16) = V2w 
c    mkraemer(17) = A2w 
c    mkraemer(21) = Cl1   
c    mkraemer(22) = Cl2   
c    mkraemer(23) = Cr1   
c    mkraemer(24) = Cr2   
c    mkraemer(25) = Ctl1   
c    mkraemer(26) = Ctl2   
c    mkraemer(27) = Ctr1   
c    mkraemer(28) = Ctr2   
c
c instead of :
c      common/susymasses/mg1,mg2,mg,ms
c      common/couplings/mx,V1,V2,A2,V1w,V2w,A2w,
c     .                 Cl1,Cl2,Cr1,Cr2,Ctl1,Ctl2,Ctr1,Ctr2
c
c i.e. no common blocks left
c      function qqreal re-written as NN_QBH including the dipole terms 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Hallo Tilman, 
c
c  hier meine subroutinen fuer born, qqbar und die beiden dipol-Terme.
c  Der Beitrag von der reellen Abstrahlung ist qqreal-qqd132_nn-qqd231_nn.
c  Uebergeben werden die Kopplungen wie von Wim definiert mit 
c  Axial -> - Axial. CL1,Cl2,... sind die t-Kanal, Ctl1,... die u-Kanal 
c  Kopplungen. Ich habe folgende Kombinationen definiert: 
c
c id Cr1^2*Cr2^2 = C4_t - Cl1^2*Cl2^2;
c id Ctr1^2*Ctr2^2 = C4_u - Ctl1^2*Ctl2^2;
c id Cr1*Cr2*Ctr1*Ctr2 = C4_tu - Cl1*Cl2*Ctl1*Ctl2;
c id A2w^2*A2^2 = C1_22 - A2w^2*V2^2;
c id V1w^2*V1^2 = C2_11;
c id V2w^2*A2^2 = C2_22 - V2w^2*V2^2;
c id V1w*V2w*V1*V2 = C2_12;
c id V2*A2*V2w*A2w = C3_22/4;
c id V1*A2*V1w*A2w = C3_12;
c id A2w*Cr1*Cr2*V2 = C5_2 - A2w*(Cr1*Cr2*A2-Cl1*Cl2*(V2-A2));
c id V1*V1w*Cr1*Cr2 = C6_1 - V1*V1w*Cl1*Cl2;
c id V2w*Cr1*Cr2*V2 = C6_2 - V2w*(Cr1*Cr2*A2+Cl1*Cl2*(V2-A2));
c id A2w*Ctr1*Ctr2*V2 = Ct5_2 - A2w*(Ctr1*Ctr2*A2-Ctl1*Ctl2*(V2-A2));
c id V1*V1w*Ctr1*Ctr2 = Ct6_1 - V1*V1w*Ctl1*Ctl2;
c id V2w*Ctr1*Ctr2*V2 = Ct6_2 - V2w*(Ctr1*Ctr2*A2+Ctl1*Ctl2*(V2-A2));
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c --------------------------------------------------------------------
c               the finite term for the total cross section 
      real*8 function NN_QBF1(iq,iout,massin,Cs,Ct,Cl,Cr)

      implicit none 

      integer    n,iout,iq
      real*8     massin(30),pi,NN_QBB
     &          ,x,sx,mu2,Cf,shift(30),D1,D2,Cs(4)
      complex*16 Ct(4),Cl(4),Cr(4) 

c               the + functions 
      D1(x) = log(1.D0-x)/(1.D0-x)
      D2(x) = (1.D0+x**2)/(1.D0-x)

      pi = 4.D0 * atan(1.D0) 
      Cf = 4.D0/3.D0      

      mu2 = massin(12)**2
      x   = massin(18)
      sx  = massin(13)

c               the born term for s -> s*x 
      do n=1,30 
         shift(n) = massin(n) 
      end do 

c               rescaling of all born varibles 
      shift(1) = massin(13)
      shift(2) = massin(14)
      shift(3) = massin(15) 
      shift(4) = massin(16)
      shift(5) = massin(17)

c               the complete kernel without +presciptions
      NN_QBF1 = - D2(x) * log(mu2/sx*x) 
     &         + 4.D0 * D1(x)
     &         - 2.D0 * (1.D0+x) * log(1.D0-x) 
     &         + ( 1.D0 - x )

c               shifted born term 
      NN_QBF1 = NN_QBF1 * NN_QBB(iq,iout,shift,Cs,Ct,Cl,Cr)

      NN_QBF1 = NN_QBF1 * Cf / (4.D0*pi**2)

      return
      end


c --------------------------------------------------------------------
      real*8 function NN_QBF2(iq,iout,massin,Cs,Ct,Cl,Cr)

      implicit none 

      integer    iout,iq
      real*8     massin(30),pi,NN_QBB,Cs(4)
     &          ,x,dum,m1,m2,s,mu2,Cf,D1,D2,D1t,D2t,meas,tau
      complex*16 Ct(4),Cl(4),Cr(4)

c               the + functions 
      D1(x) = log(1.D0-x)/(1.D0-x)
      D2(x) = (1.D0+x**2)/(1.D0-x)

      pi = 4.D0 * atan(1.D0) 
      Cf = 4.D0/3.D0

      s    = massin(1)
      m1   = massin(6)
      m2   = massin(7)
      mu2  = massin(12)**2
      x    = massin(18)
      meas = massin(19)
      tau  = massin(20)

c               the integrated + functions 
      dum = 1.D0 - tau
      D1t = - log(dum)**2/2.D0 
      D2t = - ( 3.D0/2.D0 - 2.D0*dum + dum**2/2.D0 + 2.D0*log(dum) )

      NN_QBF2 =  ( - D2(x) * log(mu2/s) 
     &            + 4.D0 * D1(x)        )
     &        + ( - D2t * log(mu2/s) 
     &            + 4.D0 * D1t          )/meas 

c               non-shifted born term 
      NN_QBF2 = NN_QBF2 * NN_QBB(iq,iout,massin,Cs,Ct,Cl,Cr)

      NN_QBF2 = NN_QBF2 * Cf / (4.D0*pi**2)
      
      return
      end


***********************************************************************
c$      real*8 function Born(massin,mkraemer)
c$      implicit real*8 (a-h,j-z)
c$      dimension massin(5),mkraemer(99)
c$ctp      common/susymasses/mg1,mg2,mg,ms
c$ctp      common/couplings/mx,V1,V2,A2,V1w,V2w,A2w,
c$ctp     .                 Cl1,Cl2,Cr1,Cr2,Ctl1,Ctl2,Ctr1,Ctr2
c$*
c$      mg1  = mkraemer(1)
c$      mg2  = mkraemer(2)
c$      mg   = mkraemer(3)
c$      ms   = mkraemer(4)
c$      mx   = mkraemer(11)
c$      V1   = mkraemer(12)
c$      V2   = mkraemer(13)
c$      A2   = mkraemer(14)
c$      V1w  = mkraemer(15)
c$      V2w  = mkraemer(16)
c$      A2w  = mkraemer(17)
c$      Cl1  = mkraemer(21)
c$      Cl2  = mkraemer(22)
c$      Cr1  = mkraemer(23)
c$      Cr2  = mkraemer(24)
c$      Ctl1 = mkraemer(25)
c$      Ctl2 = mkraemer(26)
c$      Ctr1 = mkraemer(27)
c$      Ctr2 = mkraemer(28)
c$*      
c$      s  = massin(1)
c$      t2 = massin(2)
c$      u2 = massin(3)
c$      t1 = massin(4)
c$      u1 = massin(5)
c$*
c$      mx2=mx**2
c$      t = t2 + mg2**2
c$      u = u2 + mg2**2
c$      ts  = t - ms**2 
c$      us  = u - ms**2 
c$*
c$      C4_t = Cl1**2*Cl2**2 + Cr1**2*Cr2**2
c$      C4_u = Ctl1**2*Ctl2**2 + Ctr1**2*Ctr2**2
c$      C4_tu = Cl1*Cl2*Ctl1*Ctl2 + Cr1*Cr2*Ctr1*Ctr2
c$      C1_22 = A2w**2*V2**2 + A2w**2*A2**2
c$      C2_11 = V1w**2*V1**2
c$      C2_12 = V1w*V2w*V1*V2
c$      C2_22 = V2w**2*A2**2 + V2w**2*V2**2 
c$      C3_22 = 4d0*V2*A2*V2w*A2w  
c$      C3_12 = V1*A2*V1w*A2w
c$      C5_2 = A2w*(Cr1*Cr2*A2-Cl1*Cl2*(V2-A2)) + A2w*Cr1*Cr2*V2
c$      C6_1 = V1*V1w*Cl1*Cl2 + V1*V1w*Cr1*Cr2
c$      C6_2 = V2w*(Cr1*Cr2*A2+Cl1*Cl2*(V2-A2)) + V2w*Cr1*Cr2*V2
c$      Ct5_2 = A2w*(Ctr1*Ctr2*A2-Ctl1*Ctl2*(V2-A2)) + A2w*Ctr1*Ctr2*V2 
c$      Ct6_1 = V1*V1w*Ctl1*Ctl2 + V1*V1w*Ctr1*Ctr2
c$      Ct6_2 = V2w*(Ctr1*Ctr2*A2+Ctl1*Ctl2*(V2-A2)) + V2w*Ctr1*Ctr2*V2
c$*
c$      msq_s =
c$     +  + s**(-2)*C2_11 * ( 8*t1*t2 + 8*u1*u2 )
c$      msq_s = msq_s + s**(-1)*(s-mx2)**(-1)*C2_12 * ( 16*t1*t2 + 16*u1*
c$     +    u2 )
c$      msq_s = msq_s + s**(-1)*(s-mx2)**(-1)*C3_12 * (  - 16*t1*t2 + 16*
c$     +    u1*u2 )
c$      msq_s = msq_s + s**(-1)*C2_11 * ( 16*mg1*mg2 )
c$      msq_s = msq_s + s*(s-mx2)**(-2)*C1_22 * (  - 16*mg1*mg2 )
c$      msq_s = msq_s + s*(s-mx2)**(-2)*C2_22 * ( 16*mg1*mg2 )
c$      msq_s = msq_s + (s-mx2)**(-2)*C1_22 * ( 8*t1*t2 + 8*u1*u2 )
c$      msq_s = msq_s + (s-mx2)**(-2)*C2_22 * ( 8*t1*t2 + 8*u1*u2 )
c$      msq_s = msq_s + (s-mx2)**(-2)*C3_22 * (  - 8*t1*t2 + 8*u1*u2 )
c$      msq_s = msq_s + (s-mx2)**(-1)*C2_12 * ( 32*mg1*mg2 )
c$      msq_t =
c$     +  + ts**(-2)*C4_t * ( 16*t1*t2 )
c$*
c$      msq_u =
c$     +  + us**(-2)*C4_u * ( 16*u1*u2 )
c$*
c$      msq_st =
c$     +  + s**(-1)*ts**(-1)*C6_1 * ( 16*t1*t2 )
c$      msq_st = msq_st + s*ts**(-1)*(s-mx2)**(-1)*C5_2 * ( 16*mg1*mg2 )
c$      msq_st = msq_st + s*ts**(-1)*(s-mx2)**(-1)*C6_2 * ( 16*mg1*mg2 )
c$      msq_st = msq_st + ts**(-1)*(s-mx2)**(-1)*C5_2 * (  - 16*t1*t2 )
c$      msq_st = msq_st + ts**(-1)*(s-mx2)**(-1)*C6_2 * ( 16*t1*t2 )
c$      msq_st = msq_st + ts**(-1)*C6_1 * ( 16*mg1*mg2 )
c$*
c$      msq_su =
c$     +  + s**(-1)*us**(-1)*Ct6_1 * (  - 16*u1*u2 )
c$      msq_su = msq_su + s*us**(-1)*(s-mx2)**(-1)*Ct5_2 * ( 16*mg1*mg2 )
c$      msq_su = msq_su + s*us**(-1)*(s-mx2)**(-1)*Ct6_2 * (  - 16*mg1*
c$     +    mg2 )
c$      msq_su = msq_su + us**(-1)*(s-mx2)**(-1)*Ct5_2 * (  - 16*u1*u2 )
c$      msq_su = msq_su + us**(-1)*(s-mx2)**(-1)*Ct6_2 * (  - 16*u1*u2 )
c$      msq_su = msq_su + us**(-1)*Ct6_1 * (  - 16*mg1*mg2 )
c$*
c$      msq_tu =
c$     +  + s*ts**(-1)*us**(-1)*C4_tu * (  - 32*mg1*mg2 )
c$*
c$      msq=3d0*(msq_s+msq_t+msq_u+msq_st+msq_su+msq_tu) ! 3 = colour
c$C --- colour and spin average
c$      Born=msq/4d0/3d0**2
c$ctp
c$      Born = Born * (abs(mg1)+abs(mg2))**2/4.D0
c$*
c$      return
c$      end
************************************************************************
      real*8 function NN_QBH(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/couplings/mx,V1,V2,A2,V1w,V2w,A2w,
ctp     .                 Cl1,Cl2,Cr1,Cr2,Ctl1,Ctl2,Ctr1,Ctr2
*
      mg1  = mkraemer(1)
      mg2  = mkraemer(2)
      mg   = mkraemer(3)
      ms   = mkraemer(4)
      mx   = mkraemer(11)
      V1   = mkraemer(12)
      V2   = mkraemer(13)
      A2   = mkraemer(14)
      V1w  = mkraemer(15)
      V2w  = mkraemer(16)
      A2w  = mkraemer(17)
      Cl1  = mkraemer(21)
      Cl2  = mkraemer(22)
      Cr1  = mkraemer(23)
      Cr2  = mkraemer(24)
      Ctl1 = mkraemer(25)
      Ctl2 = mkraemer(26)
      Ctr1 = mkraemer(27)
      Ctr2 = mkraemer(28)
*
      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
*
      mx2=mx**2
      t = t2 + mg2**2
      u = u2 + mg2**2
      ts  = t - ms**2 
      us  = u - ms**2 
      tpr = - s - t - u1 + mg2**2
      upr = - s - u - t1 + mg2**2
      s5 = s + tpr + upr
      u6 = - s - t - tpr + mg2**2
      u7 = - s - u - upr + mg2**2
*
      C4_t = Cl1**2*Cl2**2 + Cr1**2*Cr2**2
      C4_u = Ctl1**2*Ctl2**2 + Ctr1**2*Ctr2**2
      C4_tu = Cl1*Cl2*Ctl1*Ctl2 + Cr1*Cr2*Ctr1*Ctr2
      C1_22 = A2w**2*V2**2 + A2w**2*A2**2
      C2_11 = V1w**2*V1**2
      C2_12 = V1w*V2w*V1*V2
      C2_22 = V2w**2*A2**2 + V2w**2*V2**2 
      C3_22 = 4d0*V2*A2*V2w*A2w  
      C3_12 = V1*A2*V1w*A2w
      C5_2 = A2w*(Cr1*Cr2*A2-Cl1*Cl2*(V2-A2)) + A2w*Cr1*Cr2*V2
      C6_1 = V1*V1w*Cl1*Cl2 + V1*V1w*Cr1*Cr2
      C6_2 = V2w*(Cr1*Cr2*A2+Cl1*Cl2*(V2-A2)) + V2w*Cr1*Cr2*V2
      Ct5_2 = A2w*(Ctr1*Ctr2*A2-Ctl1*Ctl2*(V2-A2)) + A2w*Ctr1*Ctr2*V2 
      Ct6_1 = V1*V1w*Ctl1*Ctl2 + V1*V1w*Ctr1*Ctr2
      Ct6_2 = V2w*(Ctr1*Ctr2*A2+Ctl1*Ctl2*(V2-A2)) + V2w*Ctr1*Ctr2*V2
C --- the form output for |M(q+qb->g1+g2+g)|**2
C --- spin and colour average, colour factors included
C --- factor gs2 = 4*pi*alphas NOT included
      msq_s1 =
     +  + (s5-mx2)**(-2)*(tpr)**(-1)*(upr)**(-1)*C1_22 * (  - 64./
     +    9.*mg1*mg2*s**2 + 64./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*s*t
     +     - 32./9.*mg1**2*s*u - 32./9.*mg1**2*s**2 - 32./9.*mg2**2*s*t
     +     - 32./9.*mg2**2*s*u - 32./9.*mg2**2*s**2 + 32./9.*s*t**2 + 
     +    32./9.*s*u**2 + 32./9.*s**2*t + 32./9.*s**2*u + 32./9.*s**3 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(tpr)**(-1)*(upr)**(-1)*
     + C2_22 * ( 64./9.*mg1*mg2*s**2 + 64./9.*mg1**2*mg2**2*s - 32./9.*
     +    mg1**2*s*t - 32./9.*mg1**2*s*u - 32./9.*mg1**2*s**2 - 32./9.*
     +    mg2**2*s*t - 32./9.*mg2**2*s*u - 32./9.*mg2**2*s**2 + 32./9.*
     +    s*t**2 + 32./9.*s*u**2 + 32./9.*s**2*t + 32./9.*s**2*u + 32./
     +    9.*s**3 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(tpr)**(-1)*(upr)**(-1)*
     + C3_22 * ( 32./9.*s**2*t - 32./9.*s**2*u )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(tpr)**(-1)*C1_22 * (  - 
     +    64./9.*mg1*mg2*s - 32./9.*mg1*mg2*upr + 32./9.*mg1**2*mg2**2
     +     - 32./9.*mg1**2*s - 32./9.*mg1**2*u - 16./9.*mg1**2*upr - 64.
     +    /9.*mg2**2*s - 64./9.*mg2**2*t - 32./9.*mg2**2*u - 16./9.*
     +    mg2**2*upr + 32./9.*mg2**4 + 32./9.*s*t + 64./9.*s*u + 16./3.
     +    *s*upr + 64./9.*s**2 + 32./9.*t**2 + 32./9.*u*upr + 32./9.*
     +    u**2 + 16./9.*upr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(tpr)**(-1)*C2_22 * ( 64./
     +    9.*mg1*mg2*s + 32./9.*mg1*mg2*upr + 32./9.*mg1**2*mg2**2 - 32.
     +    /9.*mg1**2*s - 32./9.*mg1**2*u - 16./9.*mg1**2*upr - 64./9.*
     +    mg2**2*s - 64./9.*mg2**2*t - 32./9.*mg2**2*u - 16./9.*mg2**2*
     +    upr + 32./9.*mg2**4 + 32./9.*s*t + 64./9.*s*u + 16./3.*s*upr
     +     + 64./9.*s**2 + 32./9.*t**2 + 32./9.*u*upr + 32./9.*u**2 + 
     +    16./9.*upr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(tpr)**(-1)*C3_22 * ( 32./
     +    9.*mg1**2*s + 16./9.*mg1**2*upr + 16./9.*mg2**2*upr + 32./9.*
     +    s*t - 64./9.*s*u - 16./3.*s*upr - 32./9.*s**2 - 32./9.*u*upr
     +     - 16./9.*upr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(upr)**(-1)*C1_22 * (  - 
     +    64./9.*mg1*mg2*s - 32./9.*mg1*mg2*tpr + 32./9.*mg1**2*mg2**2
     +     - 32./9.*mg1**2*s - 32./9.*mg1**2*t - 16./9.*mg1**2*tpr - 64.
     +    /9.*mg2**2*s - 32./9.*mg2**2*t - 64./9.*mg2**2*u - 16./9.*
     +    mg2**2*tpr + 32./9.*mg2**4 + 64./9.*s*t + 32./9.*s*u + 16./3.
     +    *s*tpr + 64./9.*s**2 + 32./9.*t*tpr + 32./9.*t**2 + 32./9.*
     +    u**2 + 16./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(upr)**(-1)*C2_22 * ( 64./
     +    9.*mg1*mg2*s + 32./9.*mg1*mg2*tpr + 32./9.*mg1**2*mg2**2 - 32.
     +    /9.*mg1**2*s - 32./9.*mg1**2*t - 16./9.*mg1**2*tpr - 64./9.*
     +    mg2**2*s - 32./9.*mg2**2*t - 64./9.*mg2**2*u - 16./9.*mg2**2*
     +    tpr + 32./9.*mg2**4 + 64./9.*s*t + 32./9.*s*u + 16./3.*s*tpr
     +     + 64./9.*s**2 + 32./9.*t*tpr + 32./9.*t**2 + 32./9.*u**2 + 
     +    16./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*(upr)**(-1)*C3_22 * (  - 
     +    32./9.*mg1**2*s - 16./9.*mg1**2*tpr - 16./9.*mg2**2*tpr + 64./
     +    9.*s*t - 32./9.*s*u + 16./3.*s*tpr + 32./9.*s**2 + 32./9.*t*
     +    tpr + 16./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*C1_22 * (  - 64./9.*mg2**2
     +     + 64./9.*s + 32./9.*t + 32./9.*u + 16./9.*tpr + 16./9.*upr )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*C2_22 * (  - 64./9.*mg2**2
     +     + 64./9.*s + 32./9.*t + 32./9.*u + 16./9.*tpr + 16./9.*upr )
      msq_s1 = msq_s1 + (s5-mx2)**(-2)*C3_22 * ( 32./9.*t - 32./9.
     +    *u + 16./9.*tpr - 16./9.*upr )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (s5)**(-1)*C2_12 * ( 128./9.*mg1*mg2*s**2 + 128./9.*mg1**2
     +    *mg2**2*s - 64./9.*mg1**2*s*t - 64./9.*mg1**2*s*u - 64./9.*
     +    mg1**2*s**2 - 64./9.*mg2**2*s*t - 64./9.*mg2**2*s*u - 64./9.*
     +    mg2**2*s**2 + 64./9.*s*t**2 + 64./9.*s*u**2 + 64./9.*s**2*t
     +     + 64./9.*s**2*u + 64./9.*s**3 )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (s5)**(-1)*C3_12 * ( 64./9.*s**2*t - 64./9.*s**2*u )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(tpr)**(-1)*(s5)**(-1)*
     + C2_12 * ( 128./9.*mg1*mg2*s + 64./9.*mg1*mg2*upr + 64./9.*mg1**2
     +    *mg2**2 - 64./9.*mg1**2*s - 64./9.*mg1**2*u - 32./9.*mg1**2*
     +    upr - 128./9.*mg2**2*s - 128./9.*mg2**2*t - 64./9.*mg2**2*u
     +     - 32./9.*mg2**2*upr + 64./9.*mg2**4 + 64./9.*s*t + 128./9.*s
     +    *u + 32./3.*s*upr + 128./9.*s**2 + 64./9.*t**2 + 64./9.*u*upr
     +     + 64./9.*u**2 + 32./9.*upr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(tpr)**(-1)*(s5)**(-1)*
     + C3_12 * ( 64./9.*mg1**2*s + 32./9.*mg1**2*upr + 32./9.*mg2**2*
     +    upr + 64./9.*s*t - 128./9.*s*u - 32./3.*s*upr - 64./9.*s**2
     +     - 64./9.*u*upr - 32./9.*upr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(upr)**(-1)*(s5)**(-1)*
     + C2_12 * ( 128./9.*mg1*mg2*s + 64./9.*mg1*mg2*tpr + 64./9.*mg1**2
     +    *mg2**2 - 64./9.*mg1**2*s - 64./9.*mg1**2*t - 32./9.*mg1**2*
     +    tpr - 128./9.*mg2**2*s - 64./9.*mg2**2*t - 128./9.*mg2**2*u
     +     - 32./9.*mg2**2*tpr + 64./9.*mg2**4 + 128./9.*s*t + 64./9.*s
     +    *u + 32./3.*s*tpr + 128./9.*s**2 + 64./9.*t*tpr + 64./9.*t**2
     +     + 64./9.*u**2 + 32./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(upr)**(-1)*(s5)**(-1)*
     + C3_12 * (  - 64./9.*mg1**2*s - 32./9.*mg1**2*tpr - 32./9.*mg2**2
     +    *tpr + 128./9.*s*t - 64./9.*s*u + 32./3.*s*tpr + 64./9.*s**2
     +     + 64./9.*t*tpr + 32./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(s5)**(-1)*C2_12 * (  - 
     +    128./9.*mg2**2 + 128./9.*s + 64./9.*t + 64./9.*u + 32./9.*tpr
     +     + 32./9.*upr )
      msq_s1 = msq_s1 + (s5-mx2)**(-1)*(s5)**(-1)*C3_12 * ( 64./9.
     +    *t - 64./9.*u + 32./9.*tpr - 32./9.*upr )
      msq_s1 = msq_s1 + (tpr)**(-1)*(upr)**(-1)*(s5)**(-2)*C2_11
     +  * ( 64./9.*mg1*mg2*s**2 + 64./9.*mg1**2*mg2**2*s - 32./9.*
     +    mg1**2*s*t - 32./9.*mg1**2*s*u - 32./9.*mg1**2*s**2 - 32./9.*
     +    mg2**2*s*t - 32./9.*mg2**2*s*u - 32./9.*mg2**2*s**2 + 32./9.*
     +    s*t**2 + 32./9.*s*u**2 + 32./9.*s**2*t + 32./9.*s**2*u + 32./
     +    9.*s**3 )
      msq_s1 = msq_s1 + (tpr)**(-1)*(s5)**(-2)*C2_11 * ( 64./9.*
     +    mg1*mg2*s + 32./9.*mg1*mg2*upr + 32./9.*mg1**2*mg2**2 - 32./9.
     +    *mg1**2*s - 32./9.*mg1**2*u - 16./9.*mg1**2*upr - 64./9.*
     +    mg2**2*s - 64./9.*mg2**2*t - 32./9.*mg2**2*u - 16./9.*mg2**2*
     +    upr + 32./9.*mg2**4 + 32./9.*s*t + 64./9.*s*u + 16./3.*s*upr
     +     + 64./9.*s**2 + 32./9.*t**2 + 32./9.*u*upr + 32./9.*u**2 + 
     +    16./9.*upr**2 )
      msq_s1 = msq_s1 + (upr)**(-1)*(s5)**(-2)*C2_11 * ( 64./9.*
     +    mg1*mg2*s + 32./9.*mg1*mg2*tpr + 32./9.*mg1**2*mg2**2 - 32./9.
     +    *mg1**2*s - 32./9.*mg1**2*t - 16./9.*mg1**2*tpr - 64./9.*
     +    mg2**2*s - 32./9.*mg2**2*t - 64./9.*mg2**2*u - 16./9.*mg2**2*
     +    tpr + 32./9.*mg2**4 + 64./9.*s*t + 32./9.*s*u + 16./3.*s*tpr
     +     + 64./9.*s**2 + 32./9.*t*tpr + 32./9.*t**2 + 32./9.*u**2 + 
     +    16./9.*tpr**2 )
      msq_s1 = msq_s1 + (s5)**(-2)*C2_11 * (  - 64./9.*mg2**2 + 
     +    64./9.*s + 32./9.*t + 32./9.*u + 16./9.*tpr + 16./9.*upr )


      msq_t1 =
     +  + (ts)**(-2)*(upr)**(-1)*(u7+mg1**2-ms**2)**(-1)*C4_t
     +  * ( 64./9.*mg1**2*mg2**2*s + 64./9.*mg1**2*mg2**2*t + 64./9.*
     +    mg1**2*mg2**2*u - 64./9.*mg1**2*mg2**4 - 64./9.*mg1**2*s*t - 
     +    64./9.*mg1**2*t*u - 64./9.*mg2**2*s*t - 64./9.*mg2**2*t*u - 
     +    64./9.*mg2**2*t**2 + 64./9.*mg2**4*t + 64./9.*s*t**2 + 64./9.
     +    *t**2*u )
      msq_t1 = msq_t1 + (ts)**(-2)*(upr)**(-1)*C4_t * ( 32./9.
     +    *mg1**2*mg2**2 - 32./9.*mg1**2*t - 32./9.*mg2**2*s - 64./9.*
     +    mg2**2*t - 32./9.*mg2**2*u + 32./9.*mg2**4 + 32./9.*s*t + 32./
     +    9.*t*u + 32./9.*t**2 )
      msq_t1 = msq_t1 + (ts)**(-2)*(u7+mg1**2-ms**2)**(-2)*C4_t
     +  * (  - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t - 32./9.
     +    *mg1**2*mg2**2*u - 32./9.*mg1**2*mg2**2*upr + 32./9.*mg1**2*
     +    mg2**4 + 32./9.*mg1**2*s*t + 32./9.*mg1**2*t*u + 32./9.*
     +    mg1**2*t*upr + 32./9.*mg2**2*s*t + 64./9.*mg2**2*s*u + 64./9.
     +    *mg2**2*s*upr + 32./9.*mg2**2*s**2 + 32./9.*mg2**2*t*u + 32./
     +    9.*mg2**2*t*upr - 32./9.*mg2**2*t**2 + 64./9.*mg2**2*u*upr + 
     +    32./9.*mg2**2*u**2 + 32./9.*mg2**2*upr**2 - 64./9.*mg2**4*s
     +     - 64./9.*mg2**4*u - 64./9.*mg2**4*upr + 32./9.*mg2**6 - 64./
     +    9.*s*t*u - 64./9.*s*t*upr + 32./9.*s*t**2 - 32./9.*s**2*t - 
     +    64./9.*t*u*upr - 32./9.*t*u**2 - 32./9.*t*upr**2 + 32./9.*
     +    t**2*u + 32./9.*t**2*upr )
      msq_t1 = msq_t1 + (ts)**(-2)*(u7+mg1**2-ms**2)**(-1)*C4_t
     +  * ( 32./9.*mg2**2*s - 32./9.*mg2**2*t + 32./9.*mg2**2*u + 32./9.
     +    *mg2**2*upr - 32./9.*mg2**4 - 32./9.*s*t - 32./9.*t*u - 32./9.
     +    *t*upr + 64./9.*t**2 )
      msq_t1 = msq_t1 + (ts)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_t * ( 64./9.*mg1**2*mg2**2*s - 32.
     +    /9.*mg1**2*s*t - 32./9.*mg1**2*s*u - 32./9.*mg1**2*s**2 - 32./
     +    9.*mg2**2*s*t - 32./9.*mg2**2*s*u - 32./9.*mg2**2*s**2 + 32./
     +    9.*s*t**2 + 32./9.*s*u**2 + 64./9.*s**2*u + 32./9.*s**3 )
      msq_t1 = msq_t1 + (ts)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-2)*C4_t * (  - 64./9.*mg1**2*mg2**2*s
     +     - 64./9.*mg1**2*mg2**2*t - 64./9.*mg1**2*mg2**2*u - 64./9.*
     +    mg1**2*mg2**2*upr + 64./9.*mg1**2*mg2**4 + 64./9.*mg1**2*s*t
     +     + 64./9.*mg1**2*t*u + 64./9.*mg1**2*t*upr + 64./9.*mg2**2*s*
     +    t + 128./9.*mg2**2*s*u + 128./9.*mg2**2*s*upr + 64./9.*mg2**2
     +    *s**2 + 64./9.*mg2**2*t*u + 64./9.*mg2**2*t*upr + 128./9.*
     +    mg2**2*u*upr + 64./9.*mg2**2*u**2 + 64./9.*mg2**2*upr**2 - 64.
     +    /9.*mg2**4*s - 64./9.*mg2**4*u - 64./9.*mg2**4*upr - 128./9.*
     +    s*t*u - 128./9.*s*t*upr - 64./9.*s**2*t - 128./9.*t*u*upr - 
     +    64./9.*t*u**2 - 64./9.*t*upr**2 )
      msq_t1 = msq_t1 + (ts)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_t * (  - 32./9.*mg1**2*s - 32./9.
     +    *mg2**2*t + 32./9.*mg2**2*u + 32./9.*mg2**2*upr - 32./9.*s*t
     +     + 64./9.*s*u + 32./9.*s*upr + 64./9.*s**2 - 32./9.*t*u - 32./
     +    9.*t*upr + 32./9.*t**2 )
      msq_t1 = msq_t1 + (ts)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_t * (  - 32./9.*mg2**2*s + 32./9.
     +    *mg2**2*t - 32./9.*mg2**2*u - 32./9.*s*t + 64./9.*s*u + 32./9.
     +    *s**2 - 32./9.*t*u + 32./9.*u**2 )
      msq_t1 = msq_t1 + (ts)**(-1)*(u7+mg1**2-ms**2)**(-2)*C4_t
     +  * (  - 32./9.*mg2**2*s + 32./9.*mg2**2*t - 32./9.*mg2**2*u - 32.
     +    /9.*mg2**2*upr + 32./9.*mg2**4 - 32./9.*s*t - 32./9.*t*u - 32.
     +    /9.*t*upr )
      msq_t1 = msq_t1 + (ts)**(-1)*(u7+mg1**2-ms**2)**(-1)*C4_t
     +  * (  - 32./9.*mg2**2 + 32./9.*s + 32./9.*u )
      msq_t1 = msq_t1 + (tpr)**(-1)*(u7+mg1**2-ms**2)**(-2)*C4_t
     +  * ( 32./9.*mg1**2*mg2**2 - 32./9.*mg1**2*s - 32./9.*mg1**2*u - 
     +    32./9.*mg1**2*upr - 64./9.*mg2**2*s - 32./9.*mg2**2*t - 64./9.
     +    *mg2**2*u - 64./9.*mg2**2*upr + 32./9.*mg2**4 + 32./9.*s*t + 
     +    64./9.*s*u + 64./9.*s*upr + 32./9.*s**2 + 32./9.*t*u + 32./9.
     +    *t*upr + 64./9.*u*upr + 32./9.*u**2 + 32./9.*upr**2 )
      msq_t1 = msq_t1 + (u7+mg1**2-ms**2)**(-2)*C4_t * (  - 32./9.*
     +    mg2**2 + 32./9.*s + 32./9.*u + 32./9.*upr )


      msq_u1 =
     +  + (us)**(-2)*(tpr)**(-1)*(u6+mg1**2-ms**2)**(-1)*C4_u
     +  * ( 64./9.*mg1**2*mg2**2*s + 64./9.*mg1**2*mg2**2*t + 64./9.*
     +    mg1**2*mg2**2*u - 64./9.*mg1**2*mg2**4 - 64./9.*mg1**2*s*u - 
     +    64./9.*mg1**2*t*u - 64./9.*mg2**2*s*u - 64./9.*mg2**2*t*u - 
     +    64./9.*mg2**2*u**2 + 64./9.*mg2**4*u + 64./9.*s*u**2 + 64./9.
     +    *t*u**2 )
      msq_u1 = msq_u1 + (us)**(-2)*(tpr)**(-1)*C4_u * ( 32./9.
     +    *mg1**2*mg2**2 - 32./9.*mg1**2*u - 32./9.*mg2**2*s - 32./9.*
     +    mg2**2*t - 64./9.*mg2**2*u + 32./9.*mg2**4 + 32./9.*s*u + 32./
     +    9.*t*u + 32./9.*u**2 )
      msq_u1 = msq_u1 + (us)**(-2)*(u6+mg1**2-ms**2)**(-2)*C4_u
     +  * (  - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t - 32./9.
     +    *mg1**2*mg2**2*u - 32./9.*mg1**2*mg2**2*tpr + 32./9.*mg1**2*
     +    mg2**4 + 32./9.*mg1**2*s*u + 32./9.*mg1**2*t*u + 32./9.*
     +    mg1**2*u*tpr + 64./9.*mg2**2*s*t + 32./9.*mg2**2*s*u + 64./9.
     +    *mg2**2*s*tpr + 32./9.*mg2**2*s**2 + 32./9.*mg2**2*t*u + 64./
     +    9.*mg2**2*t*tpr + 32./9.*mg2**2*t**2 + 32./9.*mg2**2*u*tpr - 
     +    32./9.*mg2**2*u**2 + 32./9.*mg2**2*tpr**2 - 64./9.*mg2**4*s
     +     - 64./9.*mg2**4*t - 64./9.*mg2**4*tpr + 32./9.*mg2**6 - 64./
     +    9.*s*t*u - 64./9.*s*u*tpr + 32./9.*s*u**2 - 32./9.*s**2*u - 
     +    64./9.*t*u*tpr + 32./9.*t*u**2 - 32./9.*t**2*u - 32./9.*u*
     +    tpr**2 + 32./9.*u**2*tpr )
      msq_u1 = msq_u1 + (us)**(-2)*(u6+mg1**2-ms**2)**(-1)*C4_u
     +  * ( 32./9.*mg2**2*s + 32./9.*mg2**2*t - 32./9.*mg2**2*u + 32./9.
     +    *mg2**2*tpr - 32./9.*mg2**4 - 32./9.*s*u - 32./9.*t*u - 32./9.
     +    *u*tpr + 64./9.*u**2 )
      msq_u1 = msq_u1 + (us)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_u * ( 64./9.*mg1**2*mg2**2*s - 32.
     +    /9.*mg1**2*s*t - 32./9.*mg1**2*s*u - 32./9.*mg1**2*s**2 - 32./
     +    9.*mg2**2*s*t - 32./9.*mg2**2*s*u - 32./9.*mg2**2*s**2 + 32./
     +    9.*s*t**2 + 32./9.*s*u**2 + 64./9.*s**2*t + 32./9.*s**3 )
      msq_u1 = msq_u1 + (us)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_u * (  - 32./9.*mg2**2*s - 32./9.
     +    *mg2**2*t + 32./9.*mg2**2*u + 64./9.*s*t - 32./9.*s*u + 32./9.
     +    *s**2 - 32./9.*t*u + 32./9.*t**2 )
      msq_u1 = msq_u1 + (us)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-2)*C4_u * (  - 64./9.*mg1**2*mg2**2*s
     +     - 64./9.*mg1**2*mg2**2*t - 64./9.*mg1**2*mg2**2*u - 64./9.*
     +    mg1**2*mg2**2*tpr + 64./9.*mg1**2*mg2**4 + 64./9.*mg1**2*s*u
     +     + 64./9.*mg1**2*t*u + 64./9.*mg1**2*u*tpr + 128./9.*mg2**2*s
     +    *t + 64./9.*mg2**2*s*u + 128./9.*mg2**2*s*tpr + 64./9.*mg2**2
     +    *s**2 + 64./9.*mg2**2*t*u + 128./9.*mg2**2*t*tpr + 64./9.*
     +    mg2**2*t**2 + 64./9.*mg2**2*u*tpr + 64./9.*mg2**2*tpr**2 - 64.
     +    /9.*mg2**4*s - 64./9.*mg2**4*t - 64./9.*mg2**4*tpr - 128./9.*
     +    s*t*u - 128./9.*s*u*tpr - 64./9.*s**2*u - 128./9.*t*u*tpr - 
     +    64./9.*t**2*u - 64./9.*u*tpr**2 )
      msq_u1 = msq_u1 + (us)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_u * (  - 32./9.*mg1**2*s + 32./9.
     +    *mg2**2*t - 32./9.*mg2**2*u + 32./9.*mg2**2*tpr + 64./9.*s*t
     +     - 32./9.*s*u + 32./9.*s*tpr + 64./9.*s**2 - 32./9.*t*u - 32./
     +    9.*u*tpr + 32./9.*u**2 )
      msq_u1 = msq_u1 + (us)**(-1)*(u6+mg1**2-ms**2)**(-2)*C4_u
     +  * (  - 32./9.*mg2**2*s - 32./9.*mg2**2*t + 32./9.*mg2**2*u - 32.
     +    /9.*mg2**2*tpr + 32./9.*mg2**4 - 32./9.*s*u - 32./9.*t*u - 32.
     +    /9.*u*tpr )
      msq_u1 = msq_u1 + (us)**(-1)*(u6+mg1**2-ms**2)**(-1)*C4_u
     +  * (  - 32./9.*mg2**2 + 32./9.*s + 32./9.*t )
      msq_u1 = msq_u1 + (upr)**(-1)*(u6+mg1**2-ms**2)**(-2)*C4_u
     +  * ( 32./9.*mg1**2*mg2**2 - 32./9.*mg1**2*s - 32./9.*mg1**2*t - 
     +    32./9.*mg1**2*tpr - 64./9.*mg2**2*s - 64./9.*mg2**2*t - 32./9.
     +    *mg2**2*u - 64./9.*mg2**2*tpr + 32./9.*mg2**4 + 64./9.*s*t + 
     +    32./9.*s*u + 64./9.*s*tpr + 32./9.*s**2 + 32./9.*t*u + 64./9.
     +    *t*tpr + 32./9.*t**2 + 32./9.*u*tpr + 32./9.*tpr**2 )
      msq_u1 = msq_u1 + (u6+mg1**2-ms**2)**(-2)*C4_u * (  - 32./9.*
     +    mg2**2 + 32./9.*s + 32./9.*t + 32./9.*tpr )


      msq_st =
     +  + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + C5_2 * ( 32./9.*mg1*mg2*s**2 - 32./9.*mg1**2*mg2**2*s + 16./9.*
     +    mg1**2*s*t + 16./9.*mg1**2*s*u + 16./9.*mg1**2*s**2 + 16./9.*
     +    mg2**2*s*t + 16./9.*mg2**2*s*u + 16./9.*mg2**2*s**2 - 16./9.*
     +    s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*u - 16./9.*s**3 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (upr)**(-1)*C6_2 * ( 32./9.*mg1*mg2*s**2 + 32./9.*mg1**2*
     +    mg2**2*s - 16./9.*mg1**2*s*t - 16./9.*mg1**2*s*u - 16./9.*
     +    mg1**2*s**2 - 16./9.*mg2**2*s*t - 16./9.*mg2**2*s*u - 16./9.*
     +    mg2**2*s**2 + 16./9.*s*t**2 + 16./9.*s*u**2 + 32./9.*s**2*u
     +     + 16./9.*s**3 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * (  - 16./9.*mg1*mg2*s*t + 16./
     +    9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*upr + 16./9.*mg1*mg2*s**2
     +     - 16./9.*mg1*mg2*t*upr + 16./9.*mg1*mg2**3*s + 16./9.*mg1*
     +    mg2**3*upr + 32./9.*mg1**2*mg2**2*s + 32./9.*mg1**2*mg2**2*t
     +     + 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**2*upr - 32./9.
     +    *mg1**2*mg2**4 - 32./9.*mg1**2*s*t - 32./9.*mg1**2*t*u - 32./
     +    9.*mg1**2*t*upr - 16./9.*mg1**3*mg2*s - 32./9.*mg2**2*s*t - 
     +    64./9.*mg2**2*s*u - 64./9.*mg2**2*s*upr - 32./9.*mg2**2*s**2
     +     - 32./9.*mg2**2*t*u - 32./9.*mg2**2*t*upr - 64./9.*mg2**2*u*
     +    upr - 32./9.*mg2**2*u**2 - 32./9.*mg2**2*upr**2 + 32./9.*
     +    mg2**4*s + 32./9.*mg2**4*u + 32./9.*mg2**4*upr + 64./9.*s*t*u
     +     + 64./9.*s*t*upr + 32./9.*s**2*t + 64./9.*t*u*upr + 32./9.*t
     +    *u**2 + 32./9.*t*upr**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * (  - 16./9.*mg1*mg2*s*t + 16./
     +    9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*upr + 16./9.*mg1*mg2*s**2
     +     - 16./9.*mg1*mg2*t*upr + 16./9.*mg1*mg2**3*s + 16./9.*mg1*
     +    mg2**3*upr - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t
     +     - 32./9.*mg1**2*mg2**2*u - 32./9.*mg1**2*mg2**2*upr + 32./9.
     +    *mg1**2*mg2**4 + 32./9.*mg1**2*s*t + 32./9.*mg1**2*t*u + 32./
     +    9.*mg1**2*t*upr - 16./9.*mg1**3*mg2*s + 32./9.*mg2**2*s*t + 
     +    64./9.*mg2**2*s*u + 64./9.*mg2**2*s*upr + 32./9.*mg2**2*s**2
     +     + 32./9.*mg2**2*t*u + 32./9.*mg2**2*t*upr + 64./9.*mg2**2*u*
     +    upr + 32./9.*mg2**2*u**2 + 32./9.*mg2**2*upr**2 - 32./9.*
     +    mg2**4*s - 32./9.*mg2**4*u - 32./9.*mg2**4*upr - 64./9.*s*t*u
     +     - 64./9.*s*t*upr - 32./9.*s**2*t - 64./9.*t*u*upr - 32./9.*t
     +    *u**2 - 32./9.*t*upr**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     +  C5_2 * ( 32./9.*mg1*mg2*s + 16./9.*mg1**2*s + 16./9.*mg2**2*t
     +     - 16./9.*mg2**2*u - 16./9.*mg2**2*upr + 16./9.*s*t - 32./9.*
     +    s*u - 16./9.*s*upr - 32./9.*s**2 + 16./9.*t*u + 16./9.*t*upr
     +     - 16./9.*t**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     +  C6_2 * ( 32./9.*mg1*mg2*s - 16./9.*mg1**2*s - 16./9.*mg2**2*t
     +     + 16./9.*mg2**2*u + 16./9.*mg2**2*upr - 16./9.*s*t + 32./9.*
     +    s*u + 16./9.*s*upr + 32./9.*s**2 - 16./9.*t*u - 16./9.*t*upr
     +     + 16./9.*t**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * (  - 16./9.*mg1*mg2*s*t + 16./
     +    9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*tpr + 16./9.*mg1*mg2*s**2
     +     + 16./9.*mg1*mg2*u*tpr - 16./9.*mg1*mg2**3*s - 16./9.*mg1*
     +    mg2**3*tpr - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t
     +     - 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**4 + 32./9.*
     +    mg1**2*s*t + 32./9.*mg1**2*t*u + 16./9.*mg1**3*mg2*s + 32./9.
     +    *mg2**2*s*t + 32./9.*mg2**2*t*u + 32./9.*mg2**2*t**2 - 32./9.
     +    *mg2**4*t - 32./9.*s*t**2 - 32./9.*t**2*u )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * (  - 16./9.*mg1*mg2*s*t + 16./
     +    9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*tpr + 16./9.*mg1*mg2*s**2
     +     + 16./9.*mg1*mg2*u*tpr - 16./9.*mg1*mg2**3*s - 16./9.*mg1*
     +    mg2**3*tpr + 32./9.*mg1**2*mg2**2*s + 32./9.*mg1**2*mg2**2*t
     +     + 32./9.*mg1**2*mg2**2*u - 32./9.*mg1**2*mg2**4 - 32./9.*
     +    mg1**2*s*t - 32./9.*mg1**2*t*u + 16./9.*mg1**3*mg2*s - 32./9.
     +    *mg2**2*s*t - 32./9.*mg2**2*t*u - 32./9.*mg2**2*t**2 + 32./9.
     +    *mg2**4*t + 32./9.*s*t**2 + 32./9.*t**2*u )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     +  C5_2 * ( 32./9.*mg1*mg2*s + 32./9.*mg1*mg2*tpr - 32./9.*mg1**2*
     +    mg2**2 + 32./9.*mg1**2*t + 16./3.*mg2**2*s + 16./3.*mg2**2*t
     +     + 16./3.*mg2**2*u - 32./9.*mg2**4 - 16./9.*s*t - 32./9.*s*u
     +     - 16./9.*s**2 - 16./9.*t*u - 32./9.*t**2 - 16./9.*u**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     +  C6_2 * ( 32./9.*mg1*mg2*s + 32./9.*mg1*mg2*tpr + 32./9.*mg1**2*
     +    mg2**2 - 32./9.*mg1**2*t - 16./3.*mg2**2*s - 16./3.*mg2**2*t
     +     - 16./3.*mg2**2*u + 32./9.*mg2**4 + 16./9.*s*t + 32./9.*s*u
     +     + 16./9.*s**2 + 16./9.*t*u + 32./9.*t**2 + 16./9.*u**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * (  - 16./9.*mg1*mg2*t + 16./9.
     +    *mg1*mg2*u + 32./9.*s*t + 32./9.*t*u + 32./9.*t*upr - 32./9.*
     +    t**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * (  - 16./9.*mg1*mg2*t + 16./9.
     +    *mg1*mg2*u - 32./9.*s*t - 32./9.*t*u - 32./9.*t*upr + 32./9.*
     +    t**2 )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*C5_2 * ( 16.
     +    /9.*mg2**2 - 16./9.*s - 16./9.*u )
      msq_st = msq_st + (ts)**(-1)*(s5-mx2)**(-1)*C6_2 * ( 
     +     - 16./9.*mg2**2 + 16./9.*s + 16./9.*u )
      msq_st = msq_st + (ts)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (s5)**(-1)*C6_1 * ( 32./9.*mg1*mg2*s**2 + 32./9.*mg1**2*
     +    mg2**2*s - 16./9.*mg1**2*s*t - 16./9.*mg1**2*s*u - 16./9.*
     +    mg1**2*s**2 - 16./9.*mg2**2*s*t - 16./9.*mg2**2*s*u - 16./9.*
     +    mg2**2*s**2 + 16./9.*s*t**2 + 16./9.*s*u**2 + 32./9.*s**2*u
     +     + 16./9.*s**3 )
      msq_st = msq_st + (ts)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*(s5)**(-1)*C6_1 * (  - 16./9.*mg1*
     +    mg2*s*t + 16./9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*upr + 16./9.*
     +    mg1*mg2*s**2 - 16./9.*mg1*mg2*t*upr + 16./9.*mg1*mg2**3*s + 
     +    16./9.*mg1*mg2**3*upr - 32./9.*mg1**2*mg2**2*s - 32./9.*
     +    mg1**2*mg2**2*t - 32./9.*mg1**2*mg2**2*u - 32./9.*mg1**2*
     +    mg2**2*upr + 32./9.*mg1**2*mg2**4 + 32./9.*mg1**2*s*t + 32./9.
     +    *mg1**2*t*u + 32./9.*mg1**2*t*upr - 16./9.*mg1**3*mg2*s + 32./
     +    9.*mg2**2*s*t + 64./9.*mg2**2*s*u + 64./9.*mg2**2*s*upr + 32./
     +    9.*mg2**2*s**2 + 32./9.*mg2**2*t*u + 32./9.*mg2**2*t*upr + 64.
     +    /9.*mg2**2*u*upr + 32./9.*mg2**2*u**2 + 32./9.*mg2**2*upr**2
     +     - 32./9.*mg2**4*s - 32./9.*mg2**4*u - 32./9.*mg2**4*upr - 64.
     +    /9.*s*t*u - 64./9.*s*t*upr - 32./9.*s**2*t - 64./9.*t*u*upr
     +     - 32./9.*t*u**2 - 32./9.*t*upr**2 )
      msq_st = msq_st + (ts)**(-1)*(tpr)**(-1)*(s5)**(-1)*
     + C6_1 * ( 32./9.*mg1*mg2*s - 16./9.*mg1**2*s - 16./9.*mg2**2*t + 
     +    16./9.*mg2**2*u + 16./9.*mg2**2*upr - 16./9.*s*t + 32./9.*s*u
     +     + 16./9.*s*upr + 32./9.*s**2 - 16./9.*t*u - 16./9.*t*upr + 
     +    16./9.*t**2 )
      msq_st = msq_st + (ts)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*(s5)**(-1)*C6_1 * (  - 16./9.*mg1*
     +    mg2*s*t + 16./9.*mg1*mg2*s*u + 16./9.*mg1*mg2*s*tpr + 16./9.*
     +    mg1*mg2*s**2 + 16./9.*mg1*mg2*u*tpr - 16./9.*mg1*mg2**3*s - 
     +    16./9.*mg1*mg2**3*tpr + 32./9.*mg1**2*mg2**2*s + 32./9.*
     +    mg1**2*mg2**2*t + 32./9.*mg1**2*mg2**2*u - 32./9.*mg1**2*
     +    mg2**4 - 32./9.*mg1**2*s*t - 32./9.*mg1**2*t*u + 16./9.*
     +    mg1**3*mg2*s - 32./9.*mg2**2*s*t - 32./9.*mg2**2*t*u - 32./9.
     +    *mg2**2*t**2 + 32./9.*mg2**4*t + 32./9.*s*t**2 + 32./9.*t**2*
     +    u )
      msq_st = msq_st + (ts)**(-1)*(upr)**(-1)*(s5)**(-1)*
     + C6_1 * ( 32./9.*mg1*mg2*s + 32./9.*mg1*mg2*tpr + 32./9.*mg1**2*
     +    mg2**2 - 32./9.*mg1**2*t - 16./3.*mg2**2*s - 16./3.*mg2**2*t
     +     - 16./3.*mg2**2*u + 32./9.*mg2**4 + 16./9.*s*t + 32./9.*s*u
     +     + 16./9.*s**2 + 16./9.*t*u + 32./9.*t**2 + 16./9.*u**2 )
      msq_st = msq_st + (ts)**(-1)*(u7+mg1**2-ms**2)**(-1)*(s5)**(-1)
     + *C6_1 * (  - 16./9.*mg1*mg2*t + 16./9.*mg1*mg2*u - 32./9.*
     +    s*t - 32./9.*t*u - 32./9.*t*upr + 32./9.*t**2 )
      msq_st = msq_st + (ts)**(-1)*(s5)**(-1)*C6_1 * (  - 16./
     +    9.*mg2**2 + 16./9.*s + 16./9.*u )
      msq_st = msq_st + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * ( 32./9.*mg1*mg2*s**2 - 32./9.
     +    *mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*u + 16.
     +    /9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*u + 16./
     +    9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*
     +    u - 16./9.*s**3 )
      msq_st = msq_st + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * ( 32./9.*mg1*mg2*s**2 + 32./9.
     +    *mg1**2*mg2**2*s - 16./9.*mg1**2*s*t - 16./9.*mg1**2*s*u - 16.
     +    /9.*mg1**2*s**2 - 16./9.*mg2**2*s*t - 16./9.*mg2**2*s*u - 16./
     +    9.*mg2**2*s**2 + 16./9.*s*t**2 + 16./9.*s*u**2 + 32./9.*s**2*
     +    u + 16./9.*s**3 )
      msq_st = msq_st + (s5-mx2)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * ( 32./9.*mg1*mg2*s + 32./9.*
     +    mg1*mg2*upr - 32./9.*mg1**2*mg2**2 + 16./3.*mg1**2*s + 32./9.
     +    *mg1**2*u + 32./9.*mg1**2*upr + 64./9.*mg2**2*s + 16./3.*
     +    mg2**2*t + 16./3.*mg2**2*u + 16./3.*mg2**2*upr - 32./9.*
     +    mg2**4 - 16./9.*s*t - 32./3.*s*u - 80./9.*s*upr - 64./9.*s**2
     +     - 16./9.*t*u - 16./9.*t*upr - 16./9.*t**2 - 64./9.*u*upr - 
     +    32./9.*u**2 - 32./9.*upr**2 )
      msq_st = msq_st + (s5-mx2)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * ( 32./9.*mg1*mg2*s + 32./9.*
     +    mg1*mg2*upr + 32./9.*mg1**2*mg2**2 - 16./3.*mg1**2*s - 32./9.
     +    *mg1**2*u - 32./9.*mg1**2*upr - 64./9.*mg2**2*s - 16./3.*
     +    mg2**2*t - 16./3.*mg2**2*u - 16./3.*mg2**2*upr + 32./9.*
     +    mg2**4 + 16./9.*s*t + 32./3.*s*u + 80./9.*s*upr + 64./9.*s**2
     +     + 16./9.*t*u + 16./9.*t*upr + 16./9.*t**2 + 64./9.*u*upr + 
     +    32./9.*u**2 + 32./9.*upr**2 )
      msq_st = msq_st + (s5-mx2)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C5_2 * ( 32./9.*mg1*mg2*s + 16./9.*
     +    mg2**2*s - 16./9.*mg2**2*t + 16./9.*mg2**2*u + 16./9.*s*t - 
     +    32./9.*s*u - 16./9.*s**2 + 16./9.*t*u - 16./9.*u**2 )
      msq_st = msq_st + (s5-mx2)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C6_2 * ( 32./9.*mg1*mg2*s - 16./9.*
     +    mg2**2*s + 16./9.*mg2**2*t - 16./9.*mg2**2*u - 16./9.*s*t + 
     +    32./9.*s*u + 16./9.*s**2 - 16./9.*t*u + 16./9.*u**2 )
      msq_st = msq_st + (s5-mx2)**(-1)*(u7+mg1**2-ms**2)**(-1)*C5_2
     +  * ( 16./3.*mg2**2 - 16./3.*s - 16./3.*u - 32./9.*upr )
      msq_st = msq_st + (s5-mx2)**(-1)*(u7+mg1**2-ms**2)**(-1)*C6_2
     +  * (  - 16./3.*mg2**2 + 16./3.*s + 16./3.*u + 32./9.*upr )
      msq_st = msq_st + (tpr)**(-1)*(upr)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (s5)**(-1)*C6_1 * ( 32./9.*mg1*mg2*s**2 + 32./9.*mg1**2*
     +    mg2**2*s - 16./9.*mg1**2*s*t - 16./9.*mg1**2*s*u - 16./9.*
     +    mg1**2*s**2 - 16./9.*mg2**2*s*t - 16./9.*mg2**2*s*u - 16./9.*
     +    mg2**2*s**2 + 16./9.*s*t**2 + 16./9.*s*u**2 + 32./9.*s**2*u
     +     + 16./9.*s**3 )
      msq_st = msq_st + (tpr)**(-1)*(u7+mg1**2-ms**2)**(-1)*(s5)**(-1)*
     + C6_1 * ( 32./9.*mg1*mg2*s + 32./9.*mg1*mg2*upr + 32./9.*
     +    mg1**2*mg2**2 - 16./3.*mg1**2*s - 32./9.*mg1**2*u - 32./9.*
     +    mg1**2*upr - 64./9.*mg2**2*s - 16./3.*mg2**2*t - 16./3.*
     +    mg2**2*u - 16./3.*mg2**2*upr + 32./9.*mg2**4 + 16./9.*s*t + 
     +    32./3.*s*u + 80./9.*s*upr + 64./9.*s**2 + 16./9.*t*u + 16./9.
     +    *t*upr + 16./9.*t**2 + 64./9.*u*upr + 32./9.*u**2 + 32./9.*
     +    upr**2 )
      msq_st = msq_st + (upr)**(-1)*(u7+mg1**2-ms**2)**(-1)*(s5)**(-1)*
     + C6_1 * ( 32./9.*mg1*mg2*s - 16./9.*mg2**2*s + 16./9.*
     +    mg2**2*t - 16./9.*mg2**2*u - 16./9.*s*t + 32./9.*s*u + 16./9.
     +    *s**2 - 16./9.*t*u + 16./9.*u**2 )
      msq_st = msq_st + (u7+mg1**2-ms**2)**(-1)*(s5)**(-1)*C6_1
     +  * (  - 16./3.*mg2**2 + 16./3.*s + 16./3.*u + 32./9.*upr )

      msq_su =
     +  + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + Ct5_2 * ( 32./9.*mg1*mg2*s**2 - 32./9.*mg1**2*mg2**2*s + 16./9.*
     +    mg1**2*s*t + 16./9.*mg1**2*s*u + 16./9.*mg1**2*s**2 + 16./9.*
     +    mg2**2*s*t + 16./9.*mg2**2*s*u + 16./9.*mg2**2*s**2 - 16./9.*
     +    s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*t - 16./9.*s**3 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (upr)**(-1)*Ct6_2 * (  - 32./9.*mg1*mg2*s**2 - 32./9.*
     +    mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*u + 16./
     +    9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*u + 16./
     +    9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*
     +    t - 16./9.*s**3 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 16./9.*mg1*mg2*s*t - 16./9.
     +    *mg1*mg2*s*u + 16./9.*mg1*mg2*s*upr + 16./9.*mg1*mg2*s**2 + 
     +    16./9.*mg1*mg2*t*upr - 16./9.*mg1*mg2**3*s - 16./9.*mg1*
     +    mg2**3*upr - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t
     +     - 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**4 + 32./9.*
     +    mg1**2*s*u + 32./9.*mg1**2*t*u + 16./9.*mg1**3*mg2*s + 32./9.
     +    *mg2**2*s*u + 32./9.*mg2**2*t*u + 32./9.*mg2**2*u**2 - 32./9.
     +    *mg2**4*u - 32./9.*s*u**2 - 32./9.*t*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 16./9.*mg1*mg2*s*t + 16.
     +    /9.*mg1*mg2*s*u - 16./9.*mg1*mg2*s*upr - 16./9.*mg1*mg2*s**2
     +     - 16./9.*mg1*mg2*t*upr + 16./9.*mg1*mg2**3*s + 16./9.*mg1*
     +    mg2**3*upr - 32./9.*mg1**2*mg2**2*s - 32./9.*mg1**2*mg2**2*t
     +     - 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**4 + 32./9.*
     +    mg1**2*s*u + 32./9.*mg1**2*t*u - 16./9.*mg1**3*mg2*s + 32./9.
     +    *mg2**2*s*u + 32./9.*mg2**2*t*u + 32./9.*mg2**2*u**2 - 32./9.
     +    *mg2**4*u - 32./9.*s*u**2 - 32./9.*t*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)*
     +  Ct5_2 * ( 32./9.*mg1*mg2*s + 32./9.*mg1*mg2*upr - 32./9.*mg1**2
     +    *mg2**2 + 32./9.*mg1**2*u + 16./3.*mg2**2*s + 16./3.*mg2**2*t
     +     + 16./3.*mg2**2*u - 32./9.*mg2**4 - 32./9.*s*t - 16./9.*s*u
     +     - 16./9.*s**2 - 16./9.*t*u - 16./9.*t**2 - 32./9.*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(tpr)**(-1)
     + *Ct6_2 * (  - 32./9.*mg1*mg2*s - 32./9.*mg1*mg2*upr - 32./9.*
     +    mg1**2*mg2**2 + 32./9.*mg1**2*u + 16./3.*mg2**2*s + 16./3.*
     +    mg2**2*t + 16./3.*mg2**2*u - 32./9.*mg2**4 - 32./9.*s*t - 16./
     +    9.*s*u - 16./9.*s**2 - 16./9.*t*u - 16./9.*t**2 - 32./9.*u**2
     +     )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 16./9.*mg1*mg2*s*t - 16./9.
     +    *mg1*mg2*s*u + 16./9.*mg1*mg2*s*tpr + 16./9.*mg1*mg2*s**2 - 
     +    16./9.*mg1*mg2*u*tpr + 16./9.*mg1*mg2**3*s + 16./9.*mg1*
     +    mg2**3*tpr + 32./9.*mg1**2*mg2**2*s + 32./9.*mg1**2*mg2**2*t
     +     + 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**2*tpr - 32./9.
     +    *mg1**2*mg2**4 - 32./9.*mg1**2*s*u - 32./9.*mg1**2*t*u - 32./
     +    9.*mg1**2*u*tpr - 16./9.*mg1**3*mg2*s - 64./9.*mg2**2*s*t - 
     +    32./9.*mg2**2*s*u - 64./9.*mg2**2*s*tpr - 32./9.*mg2**2*s**2
     +     - 32./9.*mg2**2*t*u - 64./9.*mg2**2*t*tpr - 32./9.*mg2**2*
     +    t**2 - 32./9.*mg2**2*u*tpr - 32./9.*mg2**2*tpr**2 + 32./9.*
     +    mg2**4*s + 32./9.*mg2**4*t + 32./9.*mg2**4*tpr + 64./9.*s*t*u
     +     + 64./9.*s*u*tpr + 32./9.*s**2*u + 64./9.*t*u*tpr + 32./9.*
     +    t**2*u + 32./9.*u*tpr**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 16./9.*mg1*mg2*s*t + 16.
     +    /9.*mg1*mg2*s*u - 16./9.*mg1*mg2*s*tpr - 16./9.*mg1*mg2*s**2
     +     + 16./9.*mg1*mg2*u*tpr - 16./9.*mg1*mg2**3*s - 16./9.*mg1*
     +    mg2**3*tpr + 32./9.*mg1**2*mg2**2*s + 32./9.*mg1**2*mg2**2*t
     +     + 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*mg2**2*tpr - 32./9.
     +    *mg1**2*mg2**4 - 32./9.*mg1**2*s*u - 32./9.*mg1**2*t*u - 32./
     +    9.*mg1**2*u*tpr + 16./9.*mg1**3*mg2*s - 64./9.*mg2**2*s*t - 
     +    32./9.*mg2**2*s*u - 64./9.*mg2**2*s*tpr - 32./9.*mg2**2*s**2
     +     - 32./9.*mg2**2*t*u - 64./9.*mg2**2*t*tpr - 32./9.*mg2**2*
     +    t**2 - 32./9.*mg2**2*u*tpr - 32./9.*mg2**2*tpr**2 + 32./9.*
     +    mg2**4*s + 32./9.*mg2**4*t + 32./9.*mg2**4*tpr + 64./9.*s*t*u
     +     + 64./9.*s*u*tpr + 32./9.*s**2*u + 64./9.*t*u*tpr + 32./9.*
     +    t**2*u + 32./9.*u*tpr**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)
     + *Ct5_2 * ( 32./9.*mg1*mg2*s + 16./9.*mg1**2*s - 16./9.*mg2**2*t
     +     + 16./9.*mg2**2*u - 16./9.*mg2**2*tpr - 32./9.*s*t + 16./9.*
     +    s*u - 16./9.*s*tpr - 32./9.*s**2 + 16./9.*t*u + 16./9.*u*tpr
     +     - 16./9.*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*(upr)**(-1)
     + *Ct6_2 * (  - 32./9.*mg1*mg2*s + 16./9.*mg1**2*s - 16./9.*mg2**2
     +    *t + 16./9.*mg2**2*u - 16./9.*mg2**2*tpr - 32./9.*s*t + 16./9.
     +    *s*u - 16./9.*s*tpr - 32./9.*s**2 + 16./9.*t*u + 16./9.*u*tpr
     +     - 16./9.*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 16./9.*mg1*mg2*t - 16./9.*
     +    mg1*mg2*u + 32./9.*s*u + 32./9.*t*u + 32./9.*u*tpr - 32./9.*
     +    u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 16./9.*mg1*mg2*t + 16./
     +    9.*mg1*mg2*u + 32./9.*s*u + 32./9.*t*u + 32./9.*u*tpr - 32./9.
     +    *u**2 )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*Ct5_2 * ( 
     +    16./9.*mg2**2 - 16./9.*s - 16./9.*t )
      msq_su = msq_su + (us)**(-1)*(s5-mx2)**(-1)*Ct6_2 * ( 
     +    16./9.*mg2**2 - 16./9.*s - 16./9.*t )
      msq_su = msq_su + (us)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (s5)**(-1)*Ct6_1 * (  - 32./9.*mg1*mg2*s**2 - 32./9.*
     +    mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*u + 16./
     +    9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*u + 16./
     +    9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*
     +    t - 16./9.*s**3 )
      msq_su = msq_su + (us)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*(s5)**(-1)*Ct6_1 * (  - 16./9.*mg1*
     +    mg2*s*t + 16./9.*mg1*mg2*s*u - 16./9.*mg1*mg2*s*upr - 16./9.*
     +    mg1*mg2*s**2 - 16./9.*mg1*mg2*t*upr + 16./9.*mg1*mg2**3*s + 
     +    16./9.*mg1*mg2**3*upr - 32./9.*mg1**2*mg2**2*s - 32./9.*
     +    mg1**2*mg2**2*t - 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*
     +    mg2**4 + 32./9.*mg1**2*s*u + 32./9.*mg1**2*t*u - 16./9.*
     +    mg1**3*mg2*s + 32./9.*mg2**2*s*u + 32./9.*mg2**2*t*u + 32./9.
     +    *mg2**2*u**2 - 32./9.*mg2**4*u - 32./9.*s*u**2 - 32./9.*t*
     +    u**2 )
      msq_su = msq_su + (us)**(-1)*(tpr)**(-1)*(s5)**(-1)*
     + Ct6_1 * (  - 32./9.*mg1*mg2*s - 32./9.*mg1*mg2*upr - 32./9.*
     +    mg1**2*mg2**2 + 32./9.*mg1**2*u + 16./3.*mg2**2*s + 16./3.*
     +    mg2**2*t + 16./3.*mg2**2*u - 32./9.*mg2**4 - 32./9.*s*t - 16./
     +    9.*s*u - 16./9.*s**2 - 16./9.*t*u - 16./9.*t**2 - 32./9.*u**2
     +     )
      msq_su = msq_su + (us)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*(s5)**(-1)*Ct6_1 * (  - 16./9.*mg1*
     +    mg2*s*t + 16./9.*mg1*mg2*s*u - 16./9.*mg1*mg2*s*tpr - 16./9.*
     +    mg1*mg2*s**2 + 16./9.*mg1*mg2*u*tpr - 16./9.*mg1*mg2**3*s - 
     +    16./9.*mg1*mg2**3*tpr + 32./9.*mg1**2*mg2**2*s + 32./9.*
     +    mg1**2*mg2**2*t + 32./9.*mg1**2*mg2**2*u + 32./9.*mg1**2*
     +    mg2**2*tpr - 32./9.*mg1**2*mg2**4 - 32./9.*mg1**2*s*u - 32./9.
     +    *mg1**2*t*u - 32./9.*mg1**2*u*tpr + 16./9.*mg1**3*mg2*s - 64./
     +    9.*mg2**2*s*t - 32./9.*mg2**2*s*u - 64./9.*mg2**2*s*tpr - 32./
     +    9.*mg2**2*s**2 - 32./9.*mg2**2*t*u - 64./9.*mg2**2*t*tpr - 32.
     +    /9.*mg2**2*t**2 - 32./9.*mg2**2*u*tpr - 32./9.*mg2**2*tpr**2
     +     + 32./9.*mg2**4*s + 32./9.*mg2**4*t + 32./9.*mg2**4*tpr + 64.
     +    /9.*s*t*u + 64./9.*s*u*tpr + 32./9.*s**2*u + 64./9.*t*u*tpr
     +     + 32./9.*t**2*u + 32./9.*u*tpr**2 )
      msq_su = msq_su + (us)**(-1)*(upr)**(-1)*(s5)**(-1)*
     + Ct6_1 * (  - 32./9.*mg1*mg2*s + 16./9.*mg1**2*s - 16./9.*mg2**2*
     +    t + 16./9.*mg2**2*u - 16./9.*mg2**2*tpr - 32./9.*s*t + 16./9.
     +    *s*u - 16./9.*s*tpr - 32./9.*s**2 + 16./9.*t*u + 16./9.*u*tpr
     +     - 16./9.*u**2 )
      msq_su = msq_su + (us)**(-1)*(u6+mg1**2-ms**2)**(-1)*(s5)**(-1)
     + *Ct6_1 * (  - 16./9.*mg1*mg2*t + 16./9.*mg1*mg2*u + 32./9.
     +    *s*u + 32./9.*t*u + 32./9.*u*tpr - 32./9.*u**2 )
      msq_su = msq_su + (us)**(-1)*(s5)**(-1)*Ct6_1 * ( 16./9.
     +    *mg2**2 - 16./9.*s - 16./9.*t )
      msq_su = msq_su + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 32./9.*mg1*mg2*s**2 - 32./
     +    9.*mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*u + 
     +    16./9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*u + 
     +    16./9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./9.*
     +    s**2*t - 16./9.*s**3 )
      msq_su = msq_su + (s5-mx2)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 32./9.*mg1*mg2*s**2 - 
     +    32./9.*mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*
     +    u + 16./9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*
     +    u + 16./9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./
     +    9.*s**2*t - 16./9.*s**3 )
      msq_su = msq_su + (s5-mx2)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 32./9.*mg1*mg2*s + 16./9.*
     +    mg2**2*s + 16./9.*mg2**2*t - 16./9.*mg2**2*u - 32./9.*s*t + 
     +    16./9.*s*u - 16./9.*s**2 + 16./9.*t*u - 16./9.*t**2 )
      msq_su = msq_su + (s5-mx2)**(-1)*(tpr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 32./9.*mg1*mg2*s + 16./
     +    9.*mg2**2*s + 16./9.*mg2**2*t - 16./9.*mg2**2*u - 32./9.*s*t
     +     + 16./9.*s*u - 16./9.*s**2 + 16./9.*t*u - 16./9.*t**2 )
      msq_su = msq_su + (s5-mx2)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct5_2 * ( 32./9.*mg1*mg2*s + 32./9.*
     +    mg1*mg2*tpr - 32./9.*mg1**2*mg2**2 + 16./3.*mg1**2*s + 32./9.
     +    *mg1**2*t + 32./9.*mg1**2*tpr + 64./9.*mg2**2*s + 16./3.*
     +    mg2**2*t + 16./3.*mg2**2*u + 16./3.*mg2**2*tpr - 32./9.*
     +    mg2**4 - 32./3.*s*t - 16./9.*s*u - 80./9.*s*tpr - 64./9.*s**2
     +     - 16./9.*t*u - 64./9.*t*tpr - 32./9.*t**2 - 16./9.*u*tpr - 
     +    16./9.*u**2 - 32./9.*tpr**2 )
      msq_su = msq_su + (s5-mx2)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*Ct6_2 * (  - 32./9.*mg1*mg2*s - 32./
     +    9.*mg1*mg2*tpr - 32./9.*mg1**2*mg2**2 + 16./3.*mg1**2*s + 32./
     +    9.*mg1**2*t + 32./9.*mg1**2*tpr + 64./9.*mg2**2*s + 16./3.*
     +    mg2**2*t + 16./3.*mg2**2*u + 16./3.*mg2**2*tpr - 32./9.*
     +    mg2**4 - 32./3.*s*t - 16./9.*s*u - 80./9.*s*tpr - 64./9.*s**2
     +     - 16./9.*t*u - 64./9.*t*tpr - 32./9.*t**2 - 16./9.*u*tpr - 
     +    16./9.*u**2 - 32./9.*tpr**2 )
      msq_su = msq_su + (s5-mx2)**(-1)*(u6+mg1**2-ms**2)**(-1)*
     + Ct5_2 * ( 16./3.*mg2**2 - 16./3.*s - 16./3.*t - 32./9.*tpr )
      msq_su = msq_su + (s5-mx2)**(-1)*(u6+mg1**2-ms**2)**(-1)*
     + Ct6_2 * ( 16./3.*mg2**2 - 16./3.*s - 16./3.*t - 32./9.*tpr )
      msq_su = msq_su + (tpr)**(-1)*(upr)**(-1)*(u6+mg1**2-ms**2)**(-1)*
     + (s5)**(-1)*Ct6_1 * (  - 32./9.*mg1*mg2*s**2 - 32./9.*
     +    mg1**2*mg2**2*s + 16./9.*mg1**2*s*t + 16./9.*mg1**2*s*u + 16./
     +    9.*mg1**2*s**2 + 16./9.*mg2**2*s*t + 16./9.*mg2**2*s*u + 16./
     +    9.*mg2**2*s**2 - 16./9.*s*t**2 - 16./9.*s*u**2 - 32./9.*s**2*
     +    t - 16./9.*s**3 )
      msq_su = msq_su + (tpr)**(-1)*(u6+mg1**2-ms**2)**(-1)*(s5)**(-1)*
     + Ct6_1 * (  - 32./9.*mg1*mg2*s + 16./9.*mg2**2*s + 16./9.*
     +    mg2**2*t - 16./9.*mg2**2*u - 32./9.*s*t + 16./9.*s*u - 16./9.
     +    *s**2 + 16./9.*t*u - 16./9.*t**2 )
      msq_su = msq_su + (upr)**(-1)*(u6+mg1**2-ms**2)**(-1)*(s5)**(-1)*
     + Ct6_1 * (  - 32./9.*mg1*mg2*s - 32./9.*mg1*mg2*tpr - 32./9.
     +    *mg1**2*mg2**2 + 16./3.*mg1**2*s + 32./9.*mg1**2*t + 32./9.*
     +    mg1**2*tpr + 64./9.*mg2**2*s + 16./3.*mg2**2*t + 16./3.*
     +    mg2**2*u + 16./3.*mg2**2*tpr - 32./9.*mg2**4 - 32./3.*s*t - 
     +    16./9.*s*u - 80./9.*s*tpr - 64./9.*s**2 - 16./9.*t*u - 64./9.
     +    *t*tpr - 32./9.*t**2 - 16./9.*u*tpr - 16./9.*u**2 - 32./9.*
     +    tpr**2 )
      msq_su = msq_su + (u6+mg1**2-ms**2)**(-1)*(s5)**(-1)*Ct6_1
     +  * ( 16./3.*mg2**2 - 16./3.*s - 16./3.*t - 32./9.*tpr )

      msq_tu =
     +  + (ts)**(-1)*(us)**(-1)*(tpr)**(-1)*(upr)**(-1)*
     + C4_tu * (  - 64./9.*mg1*mg2*s**2 )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_tu * ( 32./9.*mg1*mg2*s*t - 32./9.
     +    *mg1*mg2*s*u - 32./9.*mg1*mg2*s*upr - 32./9.*mg1*mg2*s**2 + 
     +    32./9.*mg1*mg2*t*upr - 32./9.*mg1*mg2**3*s - 32./9.*mg1*
     +    mg2**3*upr + 32./9.*mg1**3*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*(tpr)**(-1)
     + *C4_tu * (  - 64./9.*mg1*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 32./9.*mg1*mg2*s*t + 32.
     +    /9.*mg1*mg2*s*u - 32./9.*mg1*mg2*s*tpr - 32./9.*mg1*mg2*s**2
     +     + 32./9.*mg1*mg2*u*tpr - 32./9.*mg1*mg2**3*s - 32./9.*mg1*
     +    mg2**3*tpr + 32./9.*mg1**3*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*(upr)**(-1)
     + *C4_tu * (  - 64./9.*mg1*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*(u6+mg1**2-ms**2)**(-1)*C4_tu * ( 32./
     +    9.*mg1*mg2*s*tpr + 32./9.*mg1*mg2*s*upr - 64./9.*mg1*mg2**3*s
     +     - 64./9.*mg1**3*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_tu * (  - 32./9.*mg1*mg2*u + 32./
     +    9.*mg1*mg2**3 )
      msq_tu = msq_tu + (ts)**(-1)*(us)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 32./9.*mg1*mg2*t + 32./
     +    9.*mg1*mg2**3 )
      msq_tu = msq_tu + (ts)**(-1)*(upr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*(u6+mg1**2-ms**2)**(-1)*C4_tu * ( 32./
     +    9.*mg1*mg2*s*t - 32./9.*mg1*mg2*s*u - 32./9.*mg1*mg2*s*tpr - 
     +    32./9.*mg1*mg2*s**2 - 32./9.*mg1*mg2*u*tpr + 32./9.*mg1*
     +    mg2**3*s + 32./9.*mg1*mg2**3*tpr - 32./9.*mg1**3*mg2*s )
      msq_tu = msq_tu + (ts)**(-1)*(upr)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 64./9.*mg1*mg2*tpr )
      msq_tu = msq_tu + (ts)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * ( 32./9.*mg1*mg2*t - 32./9.*
     +    mg1*mg2**3 )
      msq_tu = msq_tu + (us)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*(u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 
     +    32./9.*mg1*mg2*s*t + 32./9.*mg1*mg2*s*u - 32./9.*mg1*mg2*s*
     +    upr - 32./9.*mg1*mg2*s**2 - 32./9.*mg1*mg2*t*upr + 32./9.*mg1
     +    *mg2**3*s + 32./9.*mg1*mg2**3*upr - 32./9.*mg1**3*mg2*s )
      msq_tu = msq_tu + (us)**(-1)*(tpr)**(-1)*
     + (u7+mg1**2-ms**2)**(-1)*C4_tu * (  - 64./9.*mg1*mg2*upr )
      msq_tu = msq_tu + (us)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * ( 32./9.*mg1*mg2*u - 32./9.*
     +    mg1*mg2**3 )
      msq_tu = msq_tu + (tpr)**(-1)*(upr)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 64./9.*mg1*mg2*s**2 )
      msq_tu = msq_tu + (tpr)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 64./9.*mg1*mg2*s )
      msq_tu = msq_tu + (upr)**(-1)*(u7+mg1**2-ms**2)**(-1)*
     + (u6+mg1**2-ms**2)**(-1)*C4_tu * (  - 64./9.*mg1*mg2*s )
*
      qqreal=msq_s1+msq_t1+msq_u1+msq_st+msq_su+msq_tu
ctp
      NN_QBH = ( qqreal
     &        -qqd132_nn(massin,mkraemer)
     &        -qqd231_nn(massin,mkraemer)) * (abs(mg1)+abs(mg2))**2/4.D0
*
      return
      end


************************************************************************
      real*8 function qqd132_nn(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/couplings/mx,V1,V2,A2,V1w,V2w,A2w,
ctp     .                 Cl1,Cl2,Cr1,Cr2,Ctl1,Ctl2,Ctr1,Ctr2
*
      mg1  = mkraemer(1)
      mg2  = mkraemer(2)
      mg   = mkraemer(3)
      ms   = mkraemer(4)
      mx   = mkraemer(11)
      V1   = mkraemer(12)
      V2   = mkraemer(13)
      A2   = mkraemer(14)
      V1w  = mkraemer(15)
      V2w  = mkraemer(16)
      A2w  = mkraemer(17)
      Cl1  = mkraemer(21)
      Cl2  = mkraemer(22)
      Cr1  = mkraemer(23)
      Cr2  = mkraemer(24)
      Ctl1 = mkraemer(25)
      Ctl2 = mkraemer(26)
      Ctr1 = mkraemer(27)
      Ctr2 = mkraemer(28)
*
      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
*
      t = t2 + mg2**2
      u = u2 + mg2**2
      tpr = - s - t - u1 + mg2**2
      upr = - s - u - t1 + mg2**2
*
      C4_t = Cl1**2*Cl2**2 + Cr1**2*Cr2**2
      C4_u = Ctl1**2*Ctl2**2 + Ctr1**2*Ctr2**2
      C4_tu = Cl1*Cl2*Ctl1*Ctl2 + Cr1*Cr2*Ctr1*Ctr2
      C1_22 = A2w**2*V2**2 + A2w**2*A2**2
      C2_11 = V1w**2*V1**2
      C2_12 = V1w*V2w*V1*V2
      C2_22 = V2w**2*A2**2 + V2w**2*V2**2 
      C3_22 = 4d0*V2*A2*V2w*A2w  
      C3_12 = V1*A2*V1w*A2w
      C5_2 = A2w*(Cr1*Cr2*A2-Cl1*Cl2*(V2-A2)) + A2w*Cr1*Cr2*V2
      C6_1 = V1*V1w*Cl1*Cl2 + V1*V1w*Cr1*Cr2
      C6_2 = V2w*(Cr1*Cr2*A2+Cl1*Cl2*(V2-A2)) + V2w*Cr1*Cr2*V2
      Ct5_2 = A2w*(Ctr1*Ctr2*A2-Ctl1*Ctl2*(V2-A2)) + A2w*Ctr1*Ctr2*V2 
      Ct6_1 = V1*V1w*Ctl1*Ctl2 + V1*V1w*Ctr1*Ctr2
      Ct6_2 = V2w*(Ctr1*Ctr2*A2+Ctl1*Ctl2*(V2-A2)) + V2w*Ctr1*Ctr2*V2
*
      gs=1d0
      cf=4d0/3d0
C --- the form output for D132_nn
C --- spin and colour average, colour factors included
C --- factor gs2 = 4*pi*alphas NOT included
      sg132=tpr+upr+s
      sx132=tpr+upr+s-mx**2
      ts132=((2.0*mg1**2*s**2+2.0*mg1**2*s*tpr+3.0*mg1**2*s*upr+mg1**
     . 2*tpr*upr+mg1**2*upr**2+2.0*mg2**2*s**2+4.0*mg2**2*s*tpr+mg2**
     . 2*s*upr+2.0*mg2**2*tpr**2+2.0*mg2**2*tpr*upr-4.0*ms**2*s**2-
     . 4.0*ms**2*s*tpr-4.0*ms**2*s*upr-ms**2*tpr*upr-ms**2*upr**2-2.0
     . *s**3+2.0*s**2*t-4.0*s**2*tpr-2.0*s**2*u-5.0*s**2*upr+2.0*s*t*
     . tpr+3.0*s*t*upr-2.0*s*tpr**2-4.0*s*tpr*u-6.0*s*tpr*upr-3.0*s*u
     . *upr-4.0*s*upr**2+t*tpr*upr+t*upr**2-2.0*tpr**2*u-tpr**2*upr-
     . 3.0*tpr*u*upr-2.0*tpr*upr**2-u*upr**2-upr**3)*(4.0*s**2+4.0*s*
     . tpr+4.0*s*upr+tpr*upr+upr**2)*(s+tpr+upr))/(((tpr+upr+s)*(3.0*
     . s+upr)+(s+tpr)*s)**2*(tpr+upr+s))
      us132=((2.0*mg1**2*s**2+2.0*mg1**2*s*tpr+mg1**2*s*upr+2.0*mg2**
     . 2*s**2+3.0*mg2**2*s*upr-2.0*mg2**2*tpr**2-mg2**2*tpr*upr+mg2**
     . 2*upr**2-4.0*ms**2*s**2-4.0*ms**2*s*tpr-4.0*ms**2*s*upr-ms**2*
     . tpr*upr-ms**2*upr**2-2.0*s**3-2.0*s**2*t-4.0*s**2*tpr+2.0*s**2
     . *u-3.0*s**2*upr-2.0*s*t*tpr-3.0*s*t*upr-2.0*s*tpr**2+4.0*s*tpr
     . *u-3.0*s*tpr*upr+3.0*s*u*upr-s*upr**2-t*tpr*upr-t*upr**2+2.0*
     . tpr**2*u+3.0*tpr*u*upr+u*upr**2)*(4.0*s**2+4.0*s*tpr+4.0*s*upr
     . +tpr*upr+upr**2))/((tpr+upr+s)*(3.0*s+upr)+(s+tpr)*s)**2
      v132=(-2.0*((tpr+upr+2.0*s)*(tpr+upr)+2.0*s**2)*cf*gs**2)/((tpr
     . +upr)*s)
      ans10=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+
     . upr)*s*upr)*mg2**2+(2.0*s**2+2.0*s*t+2.0*s*tpr-2.0*s*u+s*upr+t
     . *upr-2.0*tpr*u-u*upr)*(2.0*s+tpr+upr)*upr)*mg1**2-(((4.0*s**2+
     . 4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2+(2.0*s+2.0*tpr+upr)
     . *(2.0*s+upr)*mg1**3*s)*mg1+(2.0*mg2**2*s**2+4.0*mg2**2*s*tpr+
     . mg2**2*s*upr+2.0*mg2**2*tpr**2+2.0*mg2**2*tpr*upr+2.0*s**3+2.0
     . *s**2*t+4.0*s**2*tpr-2.0*s**2*u+3.0*s**2*upr+2.0*s*t*tpr+3.0*s
     . *t*upr+2.0*s*tpr**2-4.0*s*tpr*u+3.0*s*tpr*upr-3.0*s*u*upr+s*
     . upr**2+t*tpr*upr+t*upr**2-2.0*tpr**2*u-3.0*tpr*u*upr-u*upr**2)
     . *(2.0*mg2**2*s-2.0*mg2**2*tpr+mg2**2*upr-2.0*s**2-2.0*s*t-2.0*
     . s*tpr+2.0*s*u-s*upr-t*upr+2.0*tpr*u+u*upr)))*ct5_2*sg132*ts132
      ans9=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+upr
     . )*s*upr)*mg2**2+(2.0*s**2+2.0*s*t+2.0*s*tpr-2.0*s*u+s*upr+t*
     . upr-2.0*tpr*u-u*upr)*(2.0*s+tpr+upr)*upr)*mg1**2+((4.0*s**2+
     . 4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2-(2.0*s+2.0*tpr+upr)
     . *(2.0*s+upr)*mg1**3*s)*mg1-(2.0*mg2**2*s**2+4.0*mg2**2*s*tpr+
     . mg2**2*s*upr+2.0*mg2**2*tpr**2+2.0*mg2**2*tpr*upr+2.0*s**3+2.0
     . *s**2*t+4.0*s**2*tpr-2.0*s**2*u+3.0*s**2*upr+2.0*s*t*tpr+3.0*s
     . *t*upr+2.0*s*tpr**2-4.0*s*tpr*u+3.0*s*tpr*upr-3.0*s*u*upr+s*
     . upr**2+t*tpr*upr+t*upr**2-2.0*tpr**2*u-3.0*tpr*u*upr-u*upr**2)
     . *(2.0*mg2**2*s-2.0*mg2**2*tpr+mg2**2*upr-2.0*s**2-2.0*s*t-2.0*
     . s*tpr+2.0*s*u-s*upr-t*upr+2.0*tpr*u+u*upr))*(ct6_1*sx132+ct6_2
     . *sg132)*ts132+ans10
      ans8=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+upr
     . )*s*upr)*mg2**2-(2.0*s**2-2.0*s*t+2.0*s*tpr+2.0*s*u+3.0*s*upr-
     . t*upr+2.0*tpr*u+tpr*upr+u*upr+upr**2)*(2.0*s+tpr+upr)*upr)*mg1
     . **2-(((4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2+(
     . 2.0*s+2.0*tpr+upr)*(2.0*s+upr)*mg1**3*s)*mg1+(2.0*mg2**2*s**2+
     . 4.0*mg2**2*s*tpr+mg2**2*s*upr+2.0*mg2**2*tpr**2+2.0*mg2**2*tpr
     . *upr-2.0*s**3+2.0*s**2*t-4.0*s**2*tpr-2.0*s**2*u-5.0*s**2*upr+
     . 2.0*s*t*tpr+3.0*s*t*upr-2.0*s*tpr**2-4.0*s*tpr*u-6.0*s*tpr*upr
     . -3.0*s*u*upr-4.0*s*upr**2+t*tpr*upr+t*upr**2-2.0*tpr**2*u-tpr
     . **2*upr-3.0*tpr*u*upr-2.0*tpr*upr**2-u*upr**2-upr**3)*(2.0*mg2
     . **2*s-2.0*mg2**2*tpr+mg2**2*upr+2.0*s**2-2.0*s*t+2.0*s*tpr+2.0
     . *s*u+3.0*s*upr-t*upr+2.0*tpr*u+tpr*upr+u*upr+upr**2)))*c5_2*
     . sg132*us132+ans9
      ans7=((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+upr)
     . *s*upr)*mg2**2-(2.0*s**2-2.0*s*t+2.0*s*tpr+2.0*s*u+3.0*s*upr-t
     . *upr+2.0*tpr*u+tpr*upr+u*upr+upr**2)*(2.0*s+tpr+upr)*upr)*mg1
     . **2+((4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2-(2.0
     . *s+2.0*tpr+upr)*(2.0*s+upr)*mg1**3*s)*mg1-(2.0*mg2**2*s**2+4.0
     . *mg2**2*s*tpr+mg2**2*s*upr+2.0*mg2**2*tpr**2+2.0*mg2**2*tpr*
     . upr-2.0*s**3+2.0*s**2*t-4.0*s**2*tpr-2.0*s**2*u-5.0*s**2*upr+
     . 2.0*s*t*tpr+3.0*s*t*upr-2.0*s*tpr**2-4.0*s*tpr*u-6.0*s*tpr*upr
     . -3.0*s*u*upr-4.0*s*upr**2+t*tpr*upr+t*upr**2-2.0*tpr**2*u-tpr
     . **2*upr-3.0*tpr*u*upr-2.0*tpr*upr**2-u*upr**2-upr**3)*(2.0*mg2
     . **2*s-2.0*mg2**2*tpr+mg2**2*upr+2.0*s**2-2.0*s*t+2.0*s*tpr+2.0
     . *s*u+3.0*s*upr-t*upr+2.0*tpr*u+tpr*upr+u*upr+upr**2))*(c6_1*
     . sx132+c6_2*sg132)*us132-2.0*((tpr+upr+s)*(3.0*s+upr)+(s+tpr)*s
     . )**2*c4_tu*mg1*mg2*sg132*sx132+ans8
      ans6=ans7*us132
      ans11=(((4.0*u+upr)*tpr+(3.0*u+upr)*upr-t*upr+(4.0*u+upr)*s-(
     . 4.0*tpr+upr+4.0*s)*mg2**2-mg1**2*upr)*(tpr+upr+s)+(2.0*(2.0*(
     . tpr+upr)+s)*s+(2.0*tpr+upr)*(tpr+upr))*mg2**2+((2.0*s+upr)*(
     . tpr+upr)+2.0*s**2)*mg1**2-(2.0*(u+upr+tpr+t)*s+(tpr+upr)*(2.0*
     . u+upr)+2.0*s**2)*(tpr+upr+s))*((2.0*tpr+upr)*u-t*upr-2.0*s**2+
     . (2.0*u-upr-2.0*tpr-2.0*t)*s-(2.0*tpr-upr-2.0*s)*mg2**2-(2.0*s+
     . upr)*mg1**2)*c4_u*sg132*sx132*ts132
      ans5=ans6+ans11
      ans4=2.0*ans5*sg132*sx132*ts132
      ans18=-((2.0*(2.0*s**2+4.0*s*tpr+s*upr+2.0*tpr**2+2.0*tpr*upr)*
     . (2.0*s-2.0*tpr+upr)*mg2**2-(4.0*s*t-4.0*s*u-2.0*s*upr+2.0*t*
     . upr-4.0*tpr*u-tpr*upr-2.0*u*upr-upr**2)*(4.0*s*tpr-2.0*s*upr+
     . 4.0*tpr**2+3.0*tpr*upr-upr**2))*mg2**2-(8.0*s**4+16.0*s**3*tpr
     . +16.0*s**3*upr+8.0*s**2*t**2-16.0*s**2*t*u-8.0*s**2*t*upr+8.0*
     . s**2*tpr**2+20.0*s**2*tpr*upr+8.0*s**2*u**2+8.0*s**2*u*upr+
     . 14.0*s**2*upr**2+8.0*s*t**2*upr-16.0*s*t*tpr*u-4.0*s*t*tpr*upr
     . -16.0*s*t*u*upr-8.0*s*t*upr**2+4.0*s*tpr**2*upr+16.0*s*tpr*u**
     . 2+12.0*s*tpr*u*upr+10.0*s*tpr*upr**2+8.0*s*u**2*upr+8.0*s*u*
     . upr**2+6.0*s*upr**3+2.0*t**2*upr**2-8.0*t*tpr*u*upr-2.0*t*tpr*
     . upr**2-4.0*t*u*upr**2-2.0*t*upr**3+8.0*tpr**2*u**2+4.0*tpr**2*
     . u*upr+tpr**2*upr**2+8.0*tpr*u**2*upr+6.0*tpr*u*upr**2+2.0*tpr*
     . upr**3+2.0*u**2*upr**2+2.0*u*upr**3+upr**4)*(s+tpr+upr)-2.0*((
     . 4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2-(2.0*s+2.0
     . *tpr+upr)*(2.0*s+upr)*mg1**3*s)*mg1)
      ans17=(4.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+upr)
     . *s*upr)*mg2**2+(4.0*s*t-4.0*s*u-2.0*s*upr+2.0*t*upr-4.0*tpr*u-
     . tpr*upr-2.0*u*upr-upr**2)*(2.0*s+tpr+upr)*upr)*mg1**2+ans18
      ans19=((2.0*c2_12*sx132+c2_22*sg132)*sg132+c2_11*sx132**2)
      ans16=ans17*ans19
      ans23=-((2.0*(2.0*s**2+4.0*s*tpr+s*upr+2.0*tpr**2+2.0*tpr*upr)*
     . (2.0*s-2.0*tpr+upr)*mg2**2-(4.0*s*t-4.0*s*u-2.0*s*upr+2.0*t*
     . upr-4.0*tpr*u-tpr*upr-2.0*u*upr-upr**2)*(4.0*s*tpr-2.0*s*upr+
     . 4.0*tpr**2+3.0*tpr*upr-upr**2))*mg2**2-(8.0*s**4+16.0*s**3*tpr
     . +16.0*s**3*upr+8.0*s**2*t**2-16.0*s**2*t*u-8.0*s**2*t*upr+8.0*
     . s**2*tpr**2+20.0*s**2*tpr*upr+8.0*s**2*u**2+8.0*s**2*u*upr+
     . 14.0*s**2*upr**2+8.0*s*t**2*upr-16.0*s*t*tpr*u-4.0*s*t*tpr*upr
     . -16.0*s*t*u*upr-8.0*s*t*upr**2+4.0*s*tpr**2*upr+16.0*s*tpr*u**
     . 2+12.0*s*tpr*u*upr+10.0*s*tpr*upr**2+8.0*s*u**2*upr+8.0*s*u*
     . upr**2+6.0*s*upr**3+2.0*t**2*upr**2-8.0*t*tpr*u*upr-2.0*t*tpr*
     . upr**2-4.0*t*u*upr**2-2.0*t*upr**3+8.0*tpr**2*u**2+4.0*tpr**2*
     . u*upr+tpr**2*upr**2+8.0*tpr*u**2*upr+6.0*tpr*u*upr**2+2.0*tpr*
     . upr**3+2.0*u**2*upr**2+2.0*u*upr**3+upr**4)*(s+tpr+upr)+2.0*((
     . 4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr*upr+upr**2)**2*mg2+(2.0*s+2.0
     . *tpr+upr)*(2.0*s+upr)*mg1**3*s)*mg1)
      ans22=(4.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(4.0*tpr+upr)
     . *s*upr)*mg2**2+(4.0*s*t-4.0*s*u-2.0*s*upr+2.0*t*upr-4.0*tpr*u-
     . tpr*upr-2.0*u*upr-upr**2)*(2.0*s+tpr+upr)*upr)*mg1**2+ans23
      ans24=sg132**2
      ans21=ans22*c1_22*ans24
      ans20=(((4.0*tpr-upr)*(tpr+upr)+2.0*(2.0*tpr-upr)*s)*mg2**2+(
     . 4.0*s*t-4.0*s*u-2.0*s*upr+2.0*t*upr-4.0*tpr*u-tpr*upr-2.0*u*
     . upr-upr**2)*(s+tpr+upr)+(tpr+upr+2.0*s)*mg1**2*upr)*((4.0*s+
     . upr)*(tpr+upr)+4.0*s**2)*(2.0*c3_12*sx132+c3_22*sg132)*sg132+
     . ans21
      ans15=ans16+ans20
      ans25=ts132**2
      ans14=ans15*ans25
      ans26=2.0*(((4.0*u+5.0*upr+4.0*tpr+4.0*s)*s+(3.0*u+2.0*upr)*upr
     . +2.0*(2.0*u+upr)*tpr-t*upr-(4.0*tpr+upr+4.0*s)*mg2**2-mg1**2*
     . upr)*(tpr+upr+s)+(2.0*(2.0*(tpr+upr)+s)*s+(2.0*tpr+upr)*(tpr+
     . upr))*mg2**2+((2.0*s+upr)*(tpr+upr)+2.0*s**2)*mg1**2-(2.0*(u+
     . upr+tpr+t)*s+(tpr+upr)*(2.0*u+upr)+2.0*s**2)*(tpr+upr+s))*((
     . 2.0*u+upr)*tpr+(u+upr)*upr-t*upr+2.0*s**2+(2.0*u+3.0*upr+2.0*
     . tpr-2.0*t)*s-(2.0*tpr-upr-2.0*s)*mg2**2-(2.0*s+upr)*mg1**2)*
     . c4_t*sg132**2*sx132**2
      ans13=ans14+ans26
      ans27=us132**2
      ans12=ans13*ans27
      ans3=ans4+ans12
      ans2=2.0*ans3*s*v132
      ans1=-ans2
      d132=ans1/(3.0*((tpr+upr+s)*(3.0*s+upr)+(s+tpr)*s)**2*sg132**2*
     . sx132**2*ts132**2*upr*us132**2)
*
      qqd132_nn=d132
*
      return
      end
************************************************************************
      real*8 function qqd231_nn(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/couplings/mx,V1,V2,A2,V1w,V2w,A2w,
ctp     .                 Cl1,Cl2,Cr1,Cr2,Ctl1,Ctl2,Ctr1,Ctr2
*
      mg1  = mkraemer(1)
      mg2  = mkraemer(2)
      mg   = mkraemer(3)
      ms   = mkraemer(4)
      mx   = mkraemer(11)
      V1   = mkraemer(12)
      V2   = mkraemer(13)
      A2   = mkraemer(14)
      V1w  = mkraemer(15)
      V2w  = mkraemer(16)
      A2w  = mkraemer(17)
      Cl1  = mkraemer(21)
      Cl2  = mkraemer(22)
      Cr1  = mkraemer(23)
      Cr2  = mkraemer(24)
      Ctl1 = mkraemer(25)
      Ctl2 = mkraemer(26)
      Ctr1 = mkraemer(27)
      Ctr2 = mkraemer(28)
*
      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
*
      t = t2 + mg2**2
      u = u2 + mg2**2
      tpr = - s - t - u1 + mg2**2
      upr = - s - u - t1 + mg2**2
*
      C4_t = Cl1**2*Cl2**2 + Cr1**2*Cr2**2
      C4_u = Ctl1**2*Ctl2**2 + Ctr1**2*Ctr2**2
      C4_tu = Cl1*Cl2*Ctl1*Ctl2 + Cr1*Cr2*Ctr1*Ctr2
      C1_22 = A2w**2*V2**2 + A2w**2*A2**2
      C2_11 = V1w**2*V1**2
      C2_12 = V1w*V2w*V1*V2
      C2_22 = V2w**2*A2**2 + V2w**2*V2**2 
      C3_22 = 4d0*V2*A2*V2w*A2w  
      C3_12 = V1*A2*V1w*A2w
      C5_2 = A2w*(Cr1*Cr2*A2-Cl1*Cl2*(V2-A2)) + A2w*Cr1*Cr2*V2
      C6_1 = V1*V1w*Cl1*Cl2 + V1*V1w*Cr1*Cr2
      C6_2 = V2w*(Cr1*Cr2*A2+Cl1*Cl2*(V2-A2)) + V2w*Cr1*Cr2*V2
      Ct5_2 = A2w*(Ctr1*Ctr2*A2-Ctl1*Ctl2*(V2-A2)) + A2w*Ctr1*Ctr2*V2 
      Ct6_1 = V1*V1w*Ctl1*Ctl2 + V1*V1w*Ctr1*Ctr2
      Ct6_2 = V2w*(Ctr1*Ctr2*A2+Ctl1*Ctl2*(V2-A2)) + V2w*Ctr1*Ctr2*V2
*
      gs=1d0
      cf=4d0/3d0
C --- the form output for D231
C --- spin and colour average, colour factors included
C --- factor gs2 = 4*pi*alphas NOT included
      sg231=tpr+upr+s
      sx231=tpr+upr+s-mx**2
      ts231=((2.0*mg1**2*s**2+mg1**2*s*tpr+2.0*mg1**2*s*upr+2.0*mg2**
     . 2*s**2+3.0*mg2**2*s*tpr+mg2**2*tpr**2-mg2**2*tpr*upr-2.0*mg2**
     . 2*upr**2-4.0*ms**2*s**2-4.0*ms**2*s*tpr-4.0*ms**2*s*upr-ms**2*
     . tpr**2-ms**2*tpr*upr-2.0*s**3+2.0*s**2*t-3.0*s**2*tpr-2.0*s**2
     . *u-4.0*s**2*upr+3.0*s*t*tpr+4.0*s*t*upr-s*tpr**2-3.0*s*tpr*u-
     . 3.0*s*tpr*upr-2.0*s*u*upr-2.0*s*upr**2+t*tpr**2+3.0*t*tpr*upr+
     . 2.0*t*upr**2-tpr**2*u-tpr*u*upr)*(4.0*s**2+4.0*s*tpr+4.0*s*upr
     . +tpr**2+tpr*upr)*(s+tpr+upr))/(((tpr+2.0*upr+2.0*s)*s+(tpr+upr
     . +s)*(2.0*s+tpr))**2*(tpr+upr+s))
      us231=((2.0*mg1**2*s**2+3.0*mg1**2*s*tpr+2.0*mg1**2*s*upr+mg1**
     . 2*tpr**2+mg1**2*tpr*upr+2.0*mg2**2*s**2+mg2**2*s*tpr+4.0*mg2**
     . 2*s*upr+2.0*mg2**2*tpr*upr+2.0*mg2**2*upr**2-4.0*ms**2*s**2-
     . 4.0*ms**2*s*tpr-4.0*ms**2*s*upr-ms**2*tpr**2-ms**2*tpr*upr-2.0
     . *s**3-2.0*s**2*t-5.0*s**2*tpr+2.0*s**2*u-4.0*s**2*upr-3.0*s*t*
     . tpr-4.0*s*t*upr-4.0*s*tpr**2+3.0*s*tpr*u-6.0*s*tpr*upr+2.0*s*u
     . *upr-2.0*s*upr**2-t*tpr**2-3.0*t*tpr*upr-2.0*t*upr**2-tpr**3+
     . tpr**2*u-2.0*tpr**2*upr+tpr*u*upr-tpr*upr**2)*(4.0*s**2+4.0*s*
     . tpr+4.0*s*upr+tpr**2+tpr*upr))/((tpr+2.0*upr+2.0*s)*s+(tpr+upr
     . +s)*(2.0*s+tpr))**2
      v231=(-2.0*((tpr+upr+2.0*s)*(tpr+upr)+2.0*s**2)*cf*gs**2)/((tpr
     . +upr)*s)
      ans10=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*
     . upr)*s*tpr)*mg2**2-(2.0*s**2+2.0*s*t+3.0*s*tpr-2.0*s*u+2.0*s*
     . upr+t*tpr+2.0*t*upr+tpr**2-tpr*u+tpr*upr)*(2.0*s+tpr+upr)*tpr)
     . *mg1**2-(((4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2
     . +(2.0*s+tpr+2.0*upr)*(2.0*s+tpr)*mg1**3*s)*mg1+(2.0*mg2**2*s**
     . 2+mg2**2*s*tpr+4.0*mg2**2*s*upr+2.0*mg2**2*tpr*upr+2.0*mg2**2*
     . upr**2-2.0*s**3-2.0*s**2*t-5.0*s**2*tpr+2.0*s**2*u-4.0*s**2*
     . upr-3.0*s*t*tpr-4.0*s*t*upr-4.0*s*tpr**2+3.0*s*tpr*u-6.0*s*tpr
     . *upr+2.0*s*u*upr-2.0*s*upr**2-t*tpr**2-3.0*t*tpr*upr-2.0*t*upr
     . **2-tpr**3+tpr**2*u-2.0*tpr**2*upr+tpr*u*upr-tpr*upr**2)*(2.0*
     . mg2**2*s+mg2**2*tpr-2.0*mg2**2*upr+2.0*s**2+2.0*s*t+3.0*s*tpr-
     . 2.0*s*u+2.0*s*upr+t*tpr+2.0*t*upr+tpr**2-tpr*u+tpr*upr)))*
     . ct5_2*sg231*ts231
      ans9=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*upr
     . )*s*tpr)*mg2**2-(2.0*s**2+2.0*s*t+3.0*s*tpr-2.0*s*u+2.0*s*upr+
     . t*tpr+2.0*t*upr+tpr**2-tpr*u+tpr*upr)*(2.0*s+tpr+upr)*tpr)*mg1
     . **2+((4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2-(2.0
     . *s+tpr+2.0*upr)*(2.0*s+tpr)*mg1**3*s)*mg1-(2.0*mg2**2*s**2+mg2
     . **2*s*tpr+4.0*mg2**2*s*upr+2.0*mg2**2*tpr*upr+2.0*mg2**2*upr**
     . 2-2.0*s**3-2.0*s**2*t-5.0*s**2*tpr+2.0*s**2*u-4.0*s**2*upr-3.0
     . *s*t*tpr-4.0*s*t*upr-4.0*s*tpr**2+3.0*s*tpr*u-6.0*s*tpr*upr+
     . 2.0*s*u*upr-2.0*s*upr**2-t*tpr**2-3.0*t*tpr*upr-2.0*t*upr**2-
     . tpr**3+tpr**2*u-2.0*tpr**2*upr+tpr*u*upr-tpr*upr**2)*(2.0*mg2
     . **2*s+mg2**2*tpr-2.0*mg2**2*upr+2.0*s**2+2.0*s*t+3.0*s*tpr-2.0
     . *s*u+2.0*s*upr+t*tpr+2.0*t*upr+tpr**2-tpr*u+tpr*upr))*(ct6_1*
     . sx231+ct6_2*sg231)*ts231+ans10
      ans8=-((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*upr
     . )*s*tpr)*mg2**2+(2.0*s**2-2.0*s*t+s*tpr+2.0*s*u+2.0*s*upr-t*
     . tpr-2.0*t*upr+tpr*u)*(2.0*s+tpr+upr)*tpr)*mg1**2-(((4.0*s**2+
     . 4.0*s*tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2+(2.0*s+tpr+2.0*upr)
     . *(2.0*s+tpr)*mg1**3*s)*mg1+(2.0*mg2**2*s**2+mg2**2*s*tpr+4.0*
     . mg2**2*s*upr+2.0*mg2**2*tpr*upr+2.0*mg2**2*upr**2+2.0*s**3-2.0
     . *s**2*t+3.0*s**2*tpr+2.0*s**2*u+4.0*s**2*upr-3.0*s*t*tpr-4.0*s
     . *t*upr+s*tpr**2+3.0*s*tpr*u+3.0*s*tpr*upr+2.0*s*u*upr+2.0*s*
     . upr**2-t*tpr**2-3.0*t*tpr*upr-2.0*t*upr**2+tpr**2*u+tpr*u*upr)
     . *(2.0*mg2**2*s+mg2**2*tpr-2.0*mg2**2*upr-2.0*s**2+2.0*s*t-s*
     . tpr-2.0*s*u-2.0*s*upr+t*tpr+2.0*t*upr-tpr*u)))*c5_2*sg231*
     . us231+ans9
      ans7=((2.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*upr)
     . *s*tpr)*mg2**2+(2.0*s**2-2.0*s*t+s*tpr+2.0*s*u+2.0*s*upr-t*tpr
     . -2.0*t*upr+tpr*u)*(2.0*s+tpr+upr)*tpr)*mg1**2+((4.0*s**2+4.0*s
     . *tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2-(2.0*s+tpr+2.0*upr)*(2.0
     . *s+tpr)*mg1**3*s)*mg1-(2.0*mg2**2*s**2+mg2**2*s*tpr+4.0*mg2**2
     . *s*upr+2.0*mg2**2*tpr*upr+2.0*mg2**2*upr**2+2.0*s**3-2.0*s**2*
     . t+3.0*s**2*tpr+2.0*s**2*u+4.0*s**2*upr-3.0*s*t*tpr-4.0*s*t*upr
     . +s*tpr**2+3.0*s*tpr*u+3.0*s*tpr*upr+2.0*s*u*upr+2.0*s*upr**2-t
     . *tpr**2-3.0*t*tpr*upr-2.0*t*upr**2+tpr**2*u+tpr*u*upr)*(2.0*
     . mg2**2*s+mg2**2*tpr-2.0*mg2**2*upr-2.0*s**2+2.0*s*t-s*tpr-2.0*
     . s*u-2.0*s*upr+t*tpr+2.0*t*upr-tpr*u))*(c6_1*sx231+c6_2*sg231)*
     . us231-2.0*((tpr+2.0*upr+2.0*s)*s+(tpr+upr+s)*(2.0*s+tpr))**2*
     . c4_tu*mg1*mg2*sg231*sx231+ans8
      ans6=ans7*us231
      ans11=-(((5.0*tpr+4.0*upr+4.0*t+4.0*s)*s-((u-2.0*upr-2.0*tpr)*
     . tpr-(3.0*tpr+4.0*upr)*t)-(tpr+4.0*upr+4.0*s)*mg2**2-mg1**2*tpr
     . )*(tpr+upr+s)+(2.0*(2.0*(tpr+upr)+s)*s+(tpr+2.0*upr)*(tpr+upr)
     . )*mg2**2+((2.0*s+tpr)*(tpr+upr)+2.0*s**2)*mg1**2-(2.0*(u+upr+
     . tpr+t)*s+(2.0*t+tpr)*(tpr+upr)+2.0*s**2)*(tpr+upr+s))*((u-upr-
     . tpr)*tpr-(tpr+2.0*upr)*t-2.0*s**2+(2.0*(u-upr)-3.0*tpr-2.0*t)*
     . s-(tpr-2.0*upr+2.0*s)*mg2**2+(2.0*s+tpr)*mg1**2)*c4_u*sg231*
     . sx231*ts231
      ans5=ans6+ans11
      ans4=2.0*ans5*sg231*sx231*ts231
      ans18=-((2.0*(2.0*s**2+s*tpr+4.0*s*upr+2.0*tpr*upr+2.0*upr**2)*
     . (2.0*s+tpr-2.0*upr)*mg2**2-(4.0*s*t+2.0*s*tpr-4.0*s*u+2.0*t*
     . tpr+4.0*t*upr+tpr**2-2.0*tpr*u+tpr*upr)*(2.0*s*tpr-4.0*s*upr+
     . tpr**2-3.0*tpr*upr-4.0*upr**2))*mg2**2-(8.0*s**4+16.0*s**3*tpr
     . +16.0*s**3*upr+8.0*s**2*t**2+8.0*s**2*t*tpr-16.0*s**2*t*u+14.0
     . *s**2*tpr**2-8.0*s**2*tpr*u+20.0*s**2*tpr*upr+8.0*s**2*u**2+
     . 8.0*s**2*upr**2+8.0*s*t**2*tpr+16.0*s*t**2*upr+8.0*s*t*tpr**2-
     . 16.0*s*t*tpr*u+12.0*s*t*tpr*upr-16.0*s*t*u*upr+6.0*s*tpr**3-
     . 8.0*s*tpr**2*u+10.0*s*tpr**2*upr+8.0*s*tpr*u**2-4.0*s*tpr*u*
     . upr+4.0*s*tpr*upr**2+2.0*t**2*tpr**2+8.0*t**2*tpr*upr+8.0*t**2
     . *upr**2+2.0*t*tpr**3-4.0*t*tpr**2*u+6.0*t*tpr**2*upr-8.0*t*tpr
     . *u*upr+4.0*t*tpr*upr**2+tpr**4-2.0*tpr**3*u+2.0*tpr**3*upr+2.0
     . *tpr**2*u**2-2.0*tpr**2*u*upr+tpr**2*upr**2)*(s+tpr+upr)-2.0*(
     . (4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2-(2.0*s+
     . tpr+2.0*upr)*(2.0*s+tpr)*mg1**3*s)*mg1)
      ans17=(4.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*upr)
     . *s*tpr)*mg2**2-(4.0*s*t+2.0*s*tpr-4.0*s*u+2.0*t*tpr+4.0*t*upr+
     . tpr**2-2.0*tpr*u+tpr*upr)*(2.0*s+tpr+upr)*tpr)*mg1**2+ans18
      ans19=((2.0*c2_12*sx231+c2_22*sg231)*sg231+c2_11*sx231**2)
      ans16=ans17*ans19
      ans23=-((2.0*(2.0*s**2+s*tpr+4.0*s*upr+2.0*tpr*upr+2.0*upr**2)*
     . (2.0*s+tpr-2.0*upr)*mg2**2-(4.0*s*t+2.0*s*tpr-4.0*s*u+2.0*t*
     . tpr+4.0*t*upr+tpr**2-2.0*tpr*u+tpr*upr)*(2.0*s*tpr-4.0*s*upr+
     . tpr**2-3.0*tpr*upr-4.0*upr**2))*mg2**2-(8.0*s**4+16.0*s**3*tpr
     . +16.0*s**3*upr+8.0*s**2*t**2+8.0*s**2*t*tpr-16.0*s**2*t*u+14.0
     . *s**2*tpr**2-8.0*s**2*tpr*u+20.0*s**2*tpr*upr+8.0*s**2*u**2+
     . 8.0*s**2*upr**2+8.0*s*t**2*tpr+16.0*s*t**2*upr+8.0*s*t*tpr**2-
     . 16.0*s*t*tpr*u+12.0*s*t*tpr*upr-16.0*s*t*u*upr+6.0*s*tpr**3-
     . 8.0*s*tpr**2*u+10.0*s*tpr**2*upr+8.0*s*tpr*u**2-4.0*s*tpr*u*
     . upr+4.0*s*tpr*upr**2+2.0*t**2*tpr**2+8.0*t**2*tpr*upr+8.0*t**2
     . *upr**2+2.0*t*tpr**3-4.0*t*tpr**2*u+6.0*t*tpr**2*upr-8.0*t*tpr
     . *u*upr+4.0*t*tpr*upr**2+tpr**4-2.0*tpr**3*u+2.0*tpr**3*upr+2.0
     . *tpr**2*u**2-2.0*tpr**2*u*upr+tpr**2*upr**2)*(s+tpr+upr)+2.0*(
     . (4.0*s**2+4.0*s*tpr+4.0*s*upr+tpr**2+tpr*upr)**2*mg2+(2.0*s+
     . tpr+2.0*upr)*(2.0*s+tpr)*mg1**3*s)*mg1)
      ans22=(4.0*((4.0*s**2+tpr*upr)*(tpr+upr)+4.0*s**3+(tpr+4.0*upr)
     . *s*tpr)*mg2**2-(4.0*s*t+2.0*s*tpr-4.0*s*u+2.0*t*tpr+4.0*t*upr+
     . tpr**2-2.0*tpr*u+tpr*upr)*(2.0*s+tpr+upr)*tpr)*mg1**2+ans23
      ans24=sg231**2
      ans21=ans22*c1_22*ans24
      ans20=(((tpr+upr)*(tpr-4.0*upr)+2.0*(tpr-2.0*upr)*s)*mg2**2+(
     . 4.0*s*t+2.0*s*tpr-4.0*s*u+2.0*t*tpr+4.0*t*upr+tpr**2-2.0*tpr*u
     . +tpr*upr)*(s+tpr+upr)-(tpr+upr+2.0*s)*mg1**2*tpr)*((4.0*s+tpr)
     . *(tpr+upr)+4.0*s**2)*(2.0*c3_12*sx231+c3_22*sg231)*sg231+ans21
      ans15=ans16+ans20
      ans25=ts231**2
      ans14=ans15*ans25
      ans26=-2.0*(((u-upr-tpr)*tpr-(3.0*tpr+4.0*upr)*t-(4.0*t+tpr)*s+
     . (tpr+4.0*upr+4.0*s)*mg2**2+mg1**2*tpr)*(tpr+upr+s)-((2.0*(2.0*
     . (tpr+upr)+s)*s+(tpr+2.0*upr)*(tpr+upr))*mg2**2+((2.0*s+tpr)*(
     . tpr+upr)+2.0*s**2)*mg1**2-(2.0*(u+upr+tpr+t)*s+(2.0*t+tpr)*(
     . tpr+upr)+2.0*s**2)*(tpr+upr+s)))*((tpr+2.0*upr)*t-tpr*u-2.0*s
     . **2-(2.0*(u+upr)+tpr-2.0*t)*s+(tpr-2.0*upr+2.0*s)*mg2**2-(2.0*
     . s+tpr)*mg1**2)*c4_t*sg231**2*sx231**2
      ans13=ans14+ans26
      ans27=us231**2
      ans12=ans13*ans27
      ans3=ans4+ans12
      ans2=2.0*ans3*s*v231
      ans1=-ans2
      d231=ans1/(3.0*((tpr+2.0*upr+2.0*s)*s+(tpr+upr+s)*(2.0*s+tpr))
     . **2*sg231**2*sx231**2*tpr*ts231**2*us231**2)
*
      qqd231_nn=d231
*
      return
      end



