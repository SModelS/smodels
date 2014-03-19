cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   all qg scaling functions :
c
c      LL_QGF(massin,Cs,Ct,Cl,Cr)
c      LL_QGH(massin,mkraemer)
c      LL_QGHSUB(massin,mkraemer)
c     using the subroutines : qgd231_ll(massin,mkraemer)
c                             qbgd132_ll(massin,mkraemer)
c
c   common block mkraemer including the coupings defined 
c    as in Xmatrix_gen_r.f 
c
c   n.b. all functions concerning qb commented out for inclusive 
c    cross sections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --------------------------------------------------------------------
c               the finite term for the crossed channel (qg)
c               born term symmetric under (k1<->k2) => also (gb)
      real*8 function LL_QGF(massin,Cs)

      implicit none 

      integer    n 
      real*8     massin(30),pi,LL_QBB
     &          ,x,sx,mu2,shift(30),Cs(4)

      pi = 4.D0 * atan(1.D0) 

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

c      a la spirix:
c mf  :  wqg = x*(1-x)
c        pqg = ( (1-x)**2 + x**2 )/2 
c mat :  wqg = ( 1 + 2*x - 3*x**2 )/4 
c        pqg = 0
c DY  :  wqg = ( 1 + 6*x - 7*x**2 )/4
c        pqg = ( (1-x)**2 + x**2 )/2   

c               the rescaled part of the finite terms       
      LL_QGF = ( - (x**2+(1.D0-x)**2)/4.D0 * log(mu2/(1.D0-x)**2/sx*x)
     &          + x * (1.D0-x)/2.D0   )
     &        * LL_QBB(shift,Cs)

c               the drell-yan type matrix element
      LL_QGF = LL_QGF 
     &      + ( (1.D0+2.D0*x-3.D0*x**2)/8.D0 )
     &        * LL_QBB(shift,Cs)


ctp    check the whole drell-yan part
ctp     &          + ( 1 + 6*x - 7*x**2 )/8.D0 ) 
ctp     &          + x * (1.D0-x)/2.D0   )

      LL_QGF = LL_QGF / (4.D0*pi**2)

      return
      end


************************************************************************
      real*8 function LL_QGH(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
      complex*16 s4sc,s3sc,qg_s1,qg_t1,qg_u1,qg_st,qg_su,qg_tu 
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/squarkw/lambda
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
      lambda = massin(30)
*
      mg12=mg1**2
      ms2=ms**2
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
      s3 = - s - t - u - tpr - upr + mg1**2 + mg2**2        
      s4 = s + t + u - mg1**2 - mg2**2
      s3s = s3 + mg2**2 - ms**2
      s4s = s4 + mg1**2 - ms**2
*
      s5 = s + tpr + upr
      u6 = - s - t - tpr + mg2**2
      u7 = - s - u - upr + mg2**2
*
      u7den=u7+mg1**2-ms**2
      u6den=u6+mg1**2-ms**2
      s5x=s5-mx2
C --- introduce complex s4s, s3s to perform PV integration 
C --- s4s --> s4s + i*ms*Gamma
C --- s3s --> s3s + i*ms*Gamma
      Gamma = lambda*ms
      s4sc = s4s + ms*Gamma*dcmplx(0d0,1d0)
      s3sc = s3s - ms*Gamma*dcmplx(0d0,1d0)
      s4s2 = s4s**2 + ms**2*Gamma**2
      s3s2 = s3s**2 + ms**2*Gamma**2
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
C --- the form output for |M(q+g->g1+g2+q)|**2
C --- spin and colour average, colour factors included
      gs=1d0
C --- factor gs2 = 4*pi*alphas NOT included
      ans5=((((tpr+2.0*upr)*tpr+2.0*(u+upr)*upr+2.0*(tpr+upr)*s3s+(
     . 2.0*(u+upr)+s)*s+2.0*(tpr+upr)*ms**2-2.0*(tpr+2.0*upr+s)*mg2**
     . 2)*mg1+2.0*((s+2.0*upr)*s+tpr**2+2.0*tpr*upr+2.0*upr**2)*mg2)*
     . mg1+(2.0*(2.0*tpr+upr)*(u+upr)+tpr**2+2.0*(tpr+upr)*s3s+(2.0*(
     . u+2.0*upr+2.0*s3s)+s)*s+2.0*(tpr+upr+2.0*s)*ms**2)*mg2**2-((
     . 2.0*ms**4+4.0*ms**2*s3s+2.0*ms**2*tpr+2.0*ms**2*upr+s**2+2.0*s
     . *u+2.0*s*upr+2.0*s3s**2+2.0*s3s*tpr+2.0*s3s*upr+tpr**2+2.0*tpr
     . *upr+2.0*u**2+2.0*u*upr+2.0*upr**2)*(s+tpr+upr)+2.0*mg2**4*s)-
     . 2.0*(mg2**2-2.0*s)*mg2**2*tpr)*c1_22
      ans4=((((tpr+2.0*upr)*tpr+2.0*(u+upr)*upr+2.0*(tpr+upr)*s3s+(
     . 2.0*(u+upr)+s)*s+2.0*(tpr+upr)*ms**2-2.0*(tpr+2.0*upr+s)*mg2**
     . 2)*mg1-2.0*((s+2.0*upr)*s+tpr**2+2.0*tpr*upr+2.0*upr**2)*mg2)*
     . mg1+(2.0*(2.0*tpr+upr)*(u+upr)+tpr**2+2.0*(tpr+upr)*s3s+(2.0*(
     . u+2.0*upr+2.0*s3s)+s)*s+2.0*(tpr+upr+2.0*s)*ms**2)*mg2**2-((
     . 2.0*ms**4+4.0*ms**2*s3s+2.0*ms**2*tpr+2.0*ms**2*upr+s**2+2.0*s
     . *u+2.0*s*upr+2.0*s3s**2+2.0*s3s*tpr+2.0*s3s*upr+tpr**2+2.0*tpr
     . *upr+2.0*u**2+2.0*u*upr+2.0*upr**2)*(s+tpr+upr)+2.0*mg2**4*s)-
     . 2.0*(mg2**2-2.0*s)*mg2**2*tpr)*c2_22+(((2.0*u+3.0*upr+tpr)*s+s
     . **2+4.0*u*upr)*s-(tpr**2+2.0*tpr*upr-2.0*u*upr+2.0*(tpr+upr)*
     . s3s)*(tpr+upr)-2.0*(tpr+upr)**2*ms**2-(s+tpr)*(s-tpr)*mg2**2-(
     . s+tpr+2.0*upr)*(s-tpr)*mg1**2-2.0*(s3s*upr-tpr*u+ms**2*upr)*s-
     . (tpr**2-2.0*upr**2+2.0*s3s*tpr+2.0*ms**2*tpr)*s)*c3_22+ans5
      ans6=s5**2
      ans3=ans4*ans6
      ans7=(((((tpr+2.0*upr)*tpr+2.0*(u+upr)*upr+2.0*(tpr+upr)*s3s+(
     . 2.0*(u+upr)+s)*s+2.0*(tpr+upr)*ms**2-2.0*(tpr+2.0*upr+s)*mg2**
     . 2)*mg1-2.0*((s+2.0*upr)*s+tpr**2+2.0*tpr*upr+2.0*upr**2)*mg2)*
     . mg1+(2.0*(2.0*tpr+upr)*(u+upr)+tpr**2+2.0*(tpr+upr)*s3s+(2.0*(
     . u+2.0*upr+2.0*s3s)+s)*s+2.0*(tpr+upr+2.0*s)*ms**2)*mg2**2-((
     . 2.0*ms**4+4.0*ms**2*s3s+2.0*ms**2*tpr+2.0*ms**2*upr+s**2+2.0*s
     . *u+2.0*s*upr+2.0*s3s**2+2.0*s3s*tpr+2.0*s3s*upr+tpr**2+2.0*tpr
     . *upr+2.0*u**2+2.0*u*upr+2.0*upr**2)*(s+tpr+upr)+2.0*mg2**4*s)-
     . 2.0*(mg2**2-2.0*s)*mg2**2*tpr)*(c2_11*s5x+2.0*c2_12*s5)+2.0*((
     . (2.0*u+3.0*upr+tpr)*s+s**2+4.0*u*upr)*s-(tpr**2+2.0*tpr*upr-
     . 2.0*u*upr+2.0*(tpr+upr)*s3s)*(tpr+upr)-2.0*(tpr+upr)**2*ms**2-
     . (s+tpr)*(s-tpr)*mg2**2-(s+tpr+2.0*upr)*(s-tpr)*mg1**2-2.0*(s3s
     . *upr-tpr*u+ms**2*upr)*s-(tpr**2-2.0*upr**2+2.0*s3s*tpr+2.0*ms
     . **2*tpr)*s)*c3_12*s5)*s5x
      ans2=ans3+ans7
      ans8=gs**2
      ans1=2.0*ans2*ans8
      qg_s1=ans1/(3.0*s*s5**2*s5x**2*tpr)

      ans3=-((2.0*mg1**2*mg2**2*s3s+mg1**2*mg2**2*tpr-2.0*mg1**2*ms**
     . 2*s3s-mg1**2*ms**2*tpr-mg1**2*s3s**2-mg1**2*s3s*tpr+mg2**4*tpr
     . -2.0*mg2**2*s*s3s-mg2**2*s*tpr+mg2**2*s3s**2+mg2**2*s3s*tpr-
     . 2.0*mg2**2*s3s*u-2.0*mg2**2*s3s*upr-mg2**2*tpr*u-mg2**2*tpr*
     . upr-ms**4*tpr+2.0*ms**2*s*s3s+ms**2*s*tpr-ms**2*s3s**2-ms**2*
     . s3s*tpr+2.0*ms**2*s3s*u+2.0*ms**2*s3s*upr+ms**2*tpr*u+ms**2*
     . tpr*upr+s*s3s**2+s*s3s*tpr-s3s**3-s3s**2*tpr+s3s**2*u+s3s**2*
     . upr+s3s*tpr*u+s3s*tpr*upr)*(mg2**2-s-u-upr)*s+(mg1**2+mg2**2-
     . ms**2-s3s-u-upr)*(mg2**2-ms**2-s3s)*tpr*u7den**2)
      ans2=((((3.0*tpr+upr)*s3s+2.0*(u+upr)*tpr+(s3s+tpr)*s+2.0*ms**2
     . *tpr)*ms**2-(((s3s+tpr)*u+s*tpr)*s-((3.0*tpr+upr)*(u+upr)+s3s*
     . upr)*s3s)+mg2**2*s*tpr)*mg2**2-((((u+upr)*tpr+s3s*upr)*s3s+(
     . tpr+upr)*(u+upr)**2-((s3s-upr)*s+s3s*u)*s)*s3s+(2.0*(u+upr)*
     . tpr+s3s*upr)*ms**4+((3.0*(u+upr)*tpr+2.0*s3s*upr)*s3s-(s3s+tpr
     . )*s**2-((u+upr)*s3s+tpr*u)*s)*ms**2)+(((2.0*tpr+upr)*(u+upr)+
     . s3s*upr+s*upr)*s3s+(2.0*(u+upr)*tpr+s3s*upr)*ms**2-2.0*((tpr+
     . upr)*s3s+tpr*upr+ms**2*tpr-mg2**2*tpr)*mg2**2)*mg1**2+(s3s**2-
     . tpr*upr)*mg2**2*s+(s*tpr-2.0*u*upr)*s*s3s-(s**2-s3s*tpr)*mg2**
     . 2*s3s-2.0*(mg1**2*mg2**2*u+s*s3s**2)*tpr-(4.0*s3s-upr+2.0*ms**
     . 2)*ms**2*s*tpr-(s3s**2-s3s*upr+2.0*upr**2+ms**4+2.0*ms**2*s3s)
     . *s*s3s-2.0*((ms**2+s3s)*mg2**2-s*s3s)*mg2**2*tpr)*u7den+ans3
      ans4=gs**2
      ans1=4.0*c4_t*ans2*ans4
      qg_t1=ans1/(3.0*s*s3s2*tpr*u7den**2)

      ans4=-((2.0*mg1**2*mg2**2*s4s-mg1**2*mg2**2*tpr-2.0*mg1**2*ms**
     . 2*s4s+mg1**2*ms**2*tpr-2.0*mg1**2*s3s*s4s+mg1**2*s3s*tpr-mg1**
     . 2*s4s**2-2.0*mg1**2*s4s*upr+mg1**2*tpr**2+mg1**2*tpr*upr-mg2**
     . 4*tpr+2.0*mg2**2*ms**2*tpr+2.0*mg2**2*s3s*tpr-mg2**2*s4s**2+
     . mg2**2*s4s*tpr-2.0*mg2**2*s4s*u+2.0*mg2**2*tpr**2-mg2**2*tpr*u
     . +2.0*mg2**2*tpr*upr-ms**4*tpr-2.0*ms**2*s3s*tpr+ms**2*s4s**2-
     . ms**2*s4s*tpr+2.0*ms**2*s4s*u-2.0*ms**2*tpr**2+ms**2*tpr*u-2.0
     . *ms**2*tpr*upr-s3s**2*tpr+s3s*s4s**2-s3s*s4s*tpr+2.0*s3s*s4s*u
     . -2.0*s3s*tpr**2+s3s*tpr*u-2.0*s3s*tpr*upr+s4s**2*u+s4s**2*upr-
     . s4s*tpr**2+2.0*s4s*tpr*u-s4s*tpr*upr+2.0*s4s*u*upr-tpr**3+tpr
     . **2*u-2.0*tpr**2*upr+tpr*u*upr-tpr*upr**2)*(mg2**2-u)*s-(mg1**
     . 2+mg2**2-ms**2-s-s3s-tpr-u-upr)*(mg2**2-ms**2-s3s-tpr-upr)*tpr
     . *us**2)
      ans3=(((2.0*(u+2.0*upr+2.0*tpr)*tpr+(tpr-upr)*s4s+4.0*s3s*tpr-(
     . s4s+tpr)*s+2.0*ms**2*tpr)*ms**2+(2.0*(tpr+u+upr)*tpr+(tpr-u-
     . upr)*s4s)*(tpr+upr)+2.0*s3s**2*tpr+(2.0*(u+2.0*upr+2.0*tpr)*
     . tpr+(tpr-upr)*s4s)*s3s+(u-upr-tpr-s3s)*(s4s+tpr)*s)*mg2**2-(((
     . (tpr+upr)*tpr+s4s*upr+(s4s+tpr)*s3s)*s+(2.0*s3s+s4s+2.0*tpr+
     . 2.0*upr)*(s3s+tpr+upr)*tpr+(4.0*(tpr+upr+s3s)*tpr+s*s4s+2.0*ms
     . **2*tpr)*ms**2)*u+2.0*(tpr+upr+s3s)*mg2**4*tpr)-((u+upr+tpr)*
     . s4s*upr-2.0*(tpr+upr)*tpr*u+(s4s*upr-2.0*tpr*u)*s3s+(s4s*upr-
     . 2.0*tpr*u)*ms**2+2.0*((u+upr+tpr)*tpr-s4s*upr+s3s*tpr)*mg2**2)
     . *mg1**2-(2.0*ms**2-s)*mg2**4*tpr+2.0*(mg2+ms)*(mg2-ms)*mg1**2*
     . mg2**2*tpr+((2.0*(tpr+upr)+s3s)*s3s*upr+(tpr*upr+u**2+upr**2)*
     . (tpr+upr)-((u-2.0*upr)*tpr-2.0*upr**2-2.0*s3s*upr)*ms**2+ms**4
     . *upr)*s4s+(((tpr+2.0*upr)*s4s-tpr*u+2.0*s3s*s4s)*ms**2+(s3s+
     . tpr+upr)*(s3s+upr)*s4s+ms**4*s4s)*s)*us+ans4
      ans5=gs**2
      ans2=4.0*c4_u*ans3*ans5
      ans1=-ans2
      qg_u1=ans1/(3.0*s*s4s2*tpr*us**2)

      ans8=((2.0*s3s+upr)*s3s+2.0*(2.0*s3s-upr)*ms**2+2.0*ms**4)*s*
     . tpr+((s3s*upr-2.0*u**2)*ms**2+s3s**2*upr)*s+((s3s**2+2.0*upr**
     . 2)*s3s+2.0*(s3s+upr)*(s3s-upr)*ms**2+ms**4*s3s)*s+(2.0*s**2-s*
     . s3s+4.0*u*upr+2.0*ms**2*u)*mg2**2*s+(2.0*(ms**2+s3s)*mg2**2-
     . 3.0*s*s3s)*mg2**2*tpr+((tpr+upr)*u+s*s3s-ms**2*s)*mg1*mg2*s
      ans7=((2.0*((tpr+upr)*s3s+(u+upr)*tpr-(u+upr+s)*s-(s-tpr)*ms**2
     . +(s-tpr)*mg2**2)*mg2**2-(((2.0*tpr+upr)*(u+upr)+s3s*upr+s*upr)
     . *s3s+(2.0*(u+upr)*tpr+s3s*upr-2.0*(u+upr+s)*s)*ms**2))*mg1+((s
     . *upr-s3s*tpr)*s+(2.0*s3s*upr+tpr**2)*upr-ms**2*s*tpr+(s+tpr+
     . upr)*(s-tpr)*mg2**2-(s-tpr)*mg1**2*upr)*mg2)*mg1+(((u+upr)*tpr
     . +s3s*upr)*s3s+(tpr+upr)*(u+upr)**2)*s3s+(2.0*(u+upr)*tpr+s3s*
     . upr)*ms**4-(2.0*(tpr+2.0*u+s)*s**2-(3.0*(u+upr)*tpr+2.0*s3s*
     . upr)*s3s)*ms**2-2.0*(s+u)*mg2**4*s-(((3.0*tpr+upr)*s3s+2.0*(u+
     . upr)*tpr-(2.0*s-3.0*s3s)*s+2.0*ms**2*tpr)*ms**2+((3.0*tpr+upr)
     . *(u+upr)+(tpr+upr)*s3s)*s3s-(4.0*(u+upr)*s-3.0*s3s**2)*s)*mg2
     . **2+(tpr+upr)*mg1*mg2*tpr*u+2.0*(mg2**2-upr)*(s3s-upr)*mg2**2*
     . s-(2.0*ms**2-s3s)*s*tpr*u+(s3s+upr-ms**2)*(s+tpr)*mg1*mg2*upr+
     . ((s3s+upr)*s3s+(s3s-4.0*upr)*ms**2)*s**2+((s3s+2.0*upr)*s3s+(
     . s3s-4.0*upr)*ms**2)*s*u-((s3s-2.0*u)*u-2.0*ms**2*upr)*mg2**2*s
     . +ans8
      ans6=ans7*c6_2
      ans11=-(ms**2*s3s-2.0*tpr*upr)*ms**2*s+(s3s+upr-ms**2)*(s+tpr)*
     . mg1*mg2*upr+(3.0*s3s**2+3.0*s3s*tpr-2.0*upr**2+3.0*ms**2*s3s)*
     . mg2**2*s+(2.0*(s**2+2.0*u*upr)-(s3s-2.0*tpr)*s)*ms**2*s+((tpr+
     . upr)*u+s*s3s-ms**2*s)*mg1*mg2*s
      ans10=(((2.0*(u+upr)*tpr+s3s*upr-2.0*(s+u)*s)*ms**2+((2.0*tpr+
     . upr)*(u+upr)+s3s*upr)*s3s-2.0*(s-tpr)*mg2**4-2.0*((tpr+upr)*
     . s3s+(u+upr)*tpr-(u+upr+s)*s-(s-tpr)*ms**2)*mg2**2)*mg1+((s*upr
     . -s3s*tpr)*s+(2.0*s3s*upr+tpr**2)*upr-ms**2*s*tpr+(s+tpr+upr)*(
     . s-tpr)*mg2**2-(s-tpr)*mg1**2*upr)*mg2)*mg1-((((u+upr)*tpr+s3s*
     . upr)*s3s+(tpr+upr)*(u+upr)**2+(s3s+upr)*s**2+((u+upr+2.0*tpr+
     . s3s)*s3s+(2.0*(u+upr)+tpr)*upr)*s)*s3s+(2.0*((s3s+2.0*tpr)*s3s
     . -2.0*(u+upr)*s)*s+(3.0*(u+upr)*tpr+2.0*s3s*upr)*s3s+(2.0*(u+
     . upr)*tpr+s3s*upr+2.0*s*tpr)*ms**2)*ms**2-(((3.0*tpr+upr)*(u+
     . upr)+(tpr+upr)*s3s)*s3s-2.0*(s+u+2.0*upr)*(s+u)*s+((3.0*tpr+
     . upr)*s3s+2.0*(u+upr)*tpr+2.0*ms**2*tpr)*ms**2+2.0*((u+upr-s3s)
     . *s+s**2-s3s*tpr-ms**2*tpr)*mg2**2)*mg2**2)+(tpr+upr)*mg1*mg2*
     . tpr*u-(s3s-2.0*upr)*ms**2*s*upr-(s3s-2.0*tpr)*ms**2*s*u-2.0*(
     . ms**2-s3s)*mg2**2*s*upr-((s+u)*mg2**2+mg1**2*upr)*(2.0*ms**2-
     . s3s)*s+(2.0*ms**2*u-s3s*tpr)*s*u+ans11
      ans9=ans10*c5_2
      ans14=-((2.0*(2.0*tpr+upr)*s3s+(u+upr)*tpr+2.0*s*s3s+(2.0*tpr+
     . upr+s)*ms**2)*ms**2+((tpr+upr)*(u+upr)+s3s*tpr)*(u+upr)+(2.0*
     . tpr+upr)*s3s**2+((tpr+2.0*upr)*upr+s3s**2)*s-((u+upr+s3s)*(3.0
     . *tpr+upr)-s**2+(3.0*tpr+upr)*ms**2-2.0*mg2**2*tpr)*mg2**2-(((
     . 2.0*tpr+upr)*s3s+(u+upr)*upr+s*upr+(2.0*tpr+upr)*ms**2-2.0*(
     . tpr+upr)*mg2**2)*mg1+2.0*(tpr**2+tpr*upr+upr**2+s*upr)*mg2)*
     . mg1-(s3s-2.0*upr+ms**2)*s*u-(s3s-upr+ms**2)*s**2-(s3s*upr-tpr*
     . u+ms**2*upr)*s-(tpr-u+s3s+ms**2)*mg2**2*s)*c5_2
      ans13=((2.0*(2.0*tpr+upr)*s3s+(u+upr)*tpr-(u+upr+s)*s+(2.0*tpr+
     . upr)*ms**2)*ms**2+((tpr+upr)*(u+upr)+s3s*tpr)*(u+upr)+(2.0*tpr
     . +upr)*s3s**2+s**2*upr-((s3s-u-s)*s+(u+upr+s3s)*(3.0*tpr+upr)+(
     . 3.0*tpr+upr+s)*ms**2)*mg2**2-(((2.0*tpr+upr)*s3s+(u+upr)*upr+s
     . *upr+(2.0*tpr+upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2
     . +tpr*upr+upr**2+s*upr)*mg2)*mg1-(s3s-tpr)*s*u-(s*s3s-2.0*u*upr
     . )*s+(2.0*mg2**2-s)*mg2**2*tpr+((s3s-upr)*s3s+(tpr+2.0*upr)*upr
     . +ms**4+2.0*ms**2*s3s)*s)*c6_2+ans14
      ans12=ans13*u7den
      ans5=ans6+ans9+ans12
      ans4=ans5*s5
      ans17=((2.0*s3s+upr)*s3s+2.0*(2.0*s3s-upr)*ms**2+2.0*ms**4)*s*
     . tpr+((s3s*upr-2.0*u**2)*ms**2+s3s**2*upr)*s+((s3s**2+2.0*upr**
     . 2)*s3s+2.0*(s3s+upr)*(s3s-upr)*ms**2+ms**4*s3s)*s+(2.0*s**2-s*
     . s3s+4.0*u*upr+2.0*ms**2*u)*mg2**2*s+(2.0*(ms**2+s3s)*mg2**2-
     . 3.0*s*s3s)*mg2**2*tpr+((tpr+upr)*u+s*s3s-ms**2*s)*mg1*mg2*s+((
     . 2.0*(2.0*tpr+upr)*s3s+(u+upr)*tpr-(u+upr+s)*s+(2.0*tpr+upr)*ms
     . **2)*ms**2+((tpr+upr)*(u+upr)+s3s*tpr)*(u+upr)+(2.0*tpr+upr)*
     . s3s**2+s**2*upr-((s3s-u-s)*s+(u+upr+s3s)*(3.0*tpr+upr)+(3.0*
     . tpr+upr+s)*ms**2)*mg2**2-(((2.0*tpr+upr)*s3s+(u+upr)*upr+s*upr
     . +(2.0*tpr+upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+tpr
     . *upr+upr**2+s*upr)*mg2)*mg1-(s3s-tpr)*s*u-(s*s3s-2.0*u*upr)*s+
     . (2.0*mg2**2-s)*mg2**2*tpr+((s3s-upr)*s3s+(tpr+2.0*upr)*upr+ms
     . **4+2.0*ms**2*s3s)*s)*u7den
      ans16=((2.0*((tpr+upr)*s3s+(u+upr)*tpr-(u+upr+s)*s-(s-tpr)*ms**
     . 2+(s-tpr)*mg2**2)*mg2**2-(((2.0*tpr+upr)*(u+upr)+s3s*upr+s*upr
     . )*s3s+(2.0*(u+upr)*tpr+s3s*upr-2.0*(u+upr+s)*s)*ms**2))*mg1+((
     . s*upr-s3s*tpr)*s+(2.0*s3s*upr+tpr**2)*upr-ms**2*s*tpr+(s+tpr+
     . upr)*(s-tpr)*mg2**2-(s-tpr)*mg1**2*upr)*mg2)*mg1+(((u+upr)*tpr
     . +s3s*upr)*s3s+(tpr+upr)*(u+upr)**2)*s3s+(2.0*(u+upr)*tpr+s3s*
     . upr)*ms**4-(2.0*(tpr+2.0*u+s)*s**2-(3.0*(u+upr)*tpr+2.0*s3s*
     . upr)*s3s)*ms**2-2.0*(s+u)*mg2**4*s-(((3.0*tpr+upr)*s3s+2.0*(u+
     . upr)*tpr-(2.0*s-3.0*s3s)*s+2.0*ms**2*tpr)*ms**2+((3.0*tpr+upr)
     . *(u+upr)+(tpr+upr)*s3s)*s3s-(4.0*(u+upr)*s-3.0*s3s**2)*s)*mg2
     . **2+(tpr+upr)*mg1*mg2*tpr*u+2.0*(mg2**2-upr)*(s3s-upr)*mg2**2*
     . s-(2.0*ms**2-s3s)*s*tpr*u+(s3s+upr-ms**2)*(s+tpr)*mg1*mg2*upr+
     . ((s3s+upr)*s3s+(s3s-4.0*upr)*ms**2)*s**2+((s3s+2.0*upr)*s3s+(
     . s3s-4.0*upr)*ms**2)*s*u-((s3s-2.0*u)*u-2.0*ms**2*upr)*mg2**2*s
     . +ans17
      ans15=ans16*c6_1*s5x
      ans3=ans4+ans15
      ans18=gs**2
      ans2=2.0*ans3*ans18
      ans1=-ans2
      qg_st=ans1/(3.0*s*s3sc*s5*s5x*tpr*u7den)

      ans8=((2.0*((u+upr+tpr)*tpr-s4s*upr+s3s*tpr-(u+upr+s4s+s3s)*s-(
     . s-tpr)*ms**2+(s-tpr)*mg2**2)*mg2**2+(u+upr+tpr)*s4s*upr-2.0*(
     . tpr+upr)*tpr*u+(s4s*upr-2.0*tpr*u)*s3s+2.0*(s4s+upr+s3s)*s*u+(
     . s4s*upr-2.0*tpr*u+2.0*s*u)*ms**2)*mg1-(((2.0*s4s+upr)*upr+(tpr
     . +upr)*s3s+(s3s+upr)*s)*s+((2.0*s4s+tpr)*(tpr+upr)+s3s*tpr)*upr
     . +(s+tpr)*(s+upr)*ms**2-(s+tpr+upr)*(s-tpr)*mg2**2+(s-tpr)*mg1
     . **2*upr)*mg2)*mg1+(tpr+upr)*mg1*mg2*tpr*u-(2.0*ms**2-s4s)*mg2
     . **2*tpr*u+((tpr+upr)*u-2.0*s*s4s)*mg1*mg2*s+(2.0*(tpr+upr+s3s)
     . *tpr-(s4s-2.0*tpr)*ms**2)*s*u
      ans7=(((u+upr)*upr-tpr**2)*s4s-2.0*(tpr+u+upr)*(tpr+upr)*tpr-
     . 2.0*s3s**2*tpr-(2.0*(u+2.0*upr+2.0*tpr)*tpr+(tpr-upr)*s4s)*s3s
     . +((3.0*(u+upr)+tpr)*s4s+2.0*(u+upr)*u+(3.0*s4s+2.0*u)*s3s)*s-(
     . 4.0*(tpr+upr)*tpr+(tpr-upr)*s4s+4.0*s3s*tpr-(3.0*s4s+2.0*u)*s+
     . 2.0*ms**2*tpr)*ms**2+2.0*((tpr+upr+s3s)*tpr-(s4s+u)*s+ms**2*
     . tpr)*mg2**2)*mg2**2-((((u+2.0*upr+tpr)*s4s+2.0*u**2+s3s*s4s)*
     . s3s+(2.0*u**2+u*upr+upr**2+tpr*upr)*s4s+2.0*(tpr+upr)*u**2)*s-
     . (((tpr*u-tpr*upr-u**2-upr**2)*s4s+2.0*(tpr+upr)*tpr*u)*(tpr+
     . upr)-(s4s*upr-2.0*tpr*u)*s3s**2+(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u)*s3s)-(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u-2.0*(s4s*upr-2.0*tpr*u)*s3s-((tpr+2.0*
     . upr)*s4s+2.0*u**2+2.0*s3s*s4s)*s-(s4s*upr-2.0*tpr*u+s*s4s)*ms
     . **2)*ms**2)+ans8
      ans6=ans7*ct6_2
      ans11=((2.0*((u+upr+tpr)*tpr-s4s*upr+s3s*tpr-(u+upr+s4s+s3s)*s-
     . (s-tpr)*ms**2+(s-tpr)*mg2**2)*mg2**2+(u+upr+tpr)*s4s*upr-2.0*(
     . tpr+upr)*tpr*u+(s4s*upr-2.0*tpr*u)*s3s+2.0*(s4s+upr+s3s)*s*u+(
     . s4s*upr-2.0*tpr*u+2.0*s*u)*ms**2)*mg1+(((2.0*s4s+upr+s3s)*s-u*
     . upr)*s-(((u-upr)*tpr-2.0*s4s*upr)*(tpr+upr)-s3s*tpr*upr)+(s**2
     . +tpr*upr)*ms**2-(s+tpr+upr)*(s-tpr)*mg2**2+(s-tpr)*mg1**2*upr)
     . *mg2)*mg1-(2.0*ms**2-s4s)*mg2**2*tpr*u+(2.0*(tpr+upr+s3s)*tpr-
     . (s4s-2.0*tpr)*ms**2)*s*u-(tpr*u-upr**2-2.0*s4s*upr-(tpr+upr)*
     . s3s-(tpr+upr)*ms**2)*mg1*mg2*s
      ans10=(((u+upr)*upr-tpr**2)*s4s-2.0*(tpr+u+upr)*(tpr+upr)*tpr-
     . 2.0*s3s**2*tpr-(2.0*(u+2.0*upr+2.0*tpr)*tpr+(tpr-upr)*s4s)*s3s
     . +((3.0*(u+upr)+tpr)*s4s+2.0*(u+upr)*u+(3.0*s4s+2.0*u)*s3s)*s-(
     . 4.0*(tpr+upr)*tpr+(tpr-upr)*s4s+4.0*s3s*tpr-(3.0*s4s+2.0*u)*s+
     . 2.0*ms**2*tpr)*ms**2+2.0*((tpr+upr+s3s)*tpr-(s4s+u)*s+ms**2*
     . tpr)*mg2**2)*mg2**2-((((u+2.0*upr+tpr)*s4s+2.0*u**2+s3s*s4s)*
     . s3s+(2.0*u**2+u*upr+upr**2+tpr*upr)*s4s+2.0*(tpr+upr)*u**2)*s-
     . (((tpr*u-tpr*upr-u**2-upr**2)*s4s+2.0*(tpr+upr)*tpr*u)*(tpr+
     . upr)-(s4s*upr-2.0*tpr*u)*s3s**2+(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u)*s3s)-(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u-2.0*(s4s*upr-2.0*tpr*u)*s3s-((tpr+2.0*
     . upr)*s4s+2.0*u**2+2.0*s3s*s4s)*s-(s4s*upr-2.0*tpr*u+s*s4s)*ms
     . **2)*ms**2)+ans11
      ans9=ans10*ct5_2
      ans15=((2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)*tpr+2.0*(2.0*tpr+
     . upr)*s3s-(u-2.0*upr-3.0*tpr-2.0*s3s)*s+(2.0*tpr+upr+s)*ms**2)*
     . ms**2-(((u-2.0*upr-3.0*tpr-s3s)*s3s-((2.0*tpr+3.0*upr)*tpr-(u-
     . upr)*upr))*s-((2.0*tpr**2+tpr*u+3.0*tpr*upr+u**2+upr**2)*(tpr+
     . upr)+(2.0*tpr+upr)*s3s**2+(2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)
     . *tpr)*s3s))-((tpr+u+upr+s3s+ms**2)*(3.0*tpr+upr)-2.0*mg2**2*
     . tpr)*mg2**2-(((2.0*tpr+3.0*upr)*tpr+(u+upr)*upr+(2.0*tpr+upr)*
     . s3s+(2.0*tpr+upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+
     . tpr*upr+upr**2+s*upr)*mg2)*mg1+(u-upr-3.0*tpr-s3s-ms**2)*mg2**
     . 2*s)*ct6_2
      ans14=((2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)*tpr+2.0*(2.0*tpr+
     . upr)*s3s-(u-2.0*upr-3.0*tpr-2.0*s3s)*s+(2.0*tpr+upr+s)*ms**2)*
     . ms**2-(((u-2.0*upr-3.0*tpr-s3s)*s3s-((2.0*tpr+3.0*upr)*tpr-(u-
     . upr)*upr))*s-((2.0*tpr**2+tpr*u+3.0*tpr*upr+u**2+upr**2)*(tpr+
     . upr)+(2.0*tpr+upr)*s3s**2+(2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)
     . *tpr)*s3s))-((tpr+u+upr+s3s+ms**2)*(3.0*tpr+upr)-2.0*mg2**2*
     . tpr)*mg2**2-(((2.0*tpr+3.0*upr)*tpr+(u+upr)*upr+(2.0*tpr+upr)*
     . s3s+(2.0*tpr+upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1+2.0*(tpr**2+
     . tpr*upr+upr**2+s*upr)*mg2)*mg1+(u-upr-3.0*tpr-s3s-ms**2)*mg2**
     . 2*s)*ct5_2+ans15
      ans13=ans14*us
      ans12=-ans13
      ans5=ans6+ans9+ans12
      ans4=ans5*s5
      ans19=-((2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)*tpr+2.0*(2.0*tpr+
     . upr)*s3s-(u-2.0*upr-3.0*tpr-2.0*s3s)*s+(2.0*tpr+upr+s)*ms**2)*
     . ms**2-(((u-2.0*upr-3.0*tpr-s3s)*s3s-((2.0*tpr+3.0*upr)*tpr-(u-
     . upr)*upr))*s-((2.0*tpr**2+tpr*u+3.0*tpr*upr+u**2+upr**2)*(tpr+
     . upr)+(2.0*tpr+upr)*s3s**2+(2.0*(2.0*tpr**2+upr**2)+(u+6.0*upr)
     . *tpr)*s3s))-((tpr+u+upr+s3s+ms**2)*(3.0*tpr+upr)-2.0*mg2**2*
     . tpr)*mg2**2-(((2.0*tpr+3.0*upr)*tpr+(u+upr)*upr+(2.0*tpr+upr)*
     . s3s+(2.0*tpr+upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+
     . tpr*upr+upr**2+s*upr)*mg2)*mg1+(u-upr-3.0*tpr-s3s-ms**2)*mg2**
     . 2*s)*us
      ans18=((2.0*((u+upr+tpr)*tpr-s4s*upr+s3s*tpr-(u+upr+s4s+s3s)*s-
     . (s-tpr)*ms**2+(s-tpr)*mg2**2)*mg2**2+(u+upr+tpr)*s4s*upr-2.0*(
     . tpr+upr)*tpr*u+(s4s*upr-2.0*tpr*u)*s3s+2.0*(s4s+upr+s3s)*s*u+(
     . s4s*upr-2.0*tpr*u+2.0*s*u)*ms**2)*mg1-(((2.0*s4s+upr)*upr+(tpr
     . +upr)*s3s+(s3s+upr)*s)*s+((2.0*s4s+tpr)*(tpr+upr)+s3s*tpr)*upr
     . +(s+tpr)*(s+upr)*ms**2-(s+tpr+upr)*(s-tpr)*mg2**2+(s-tpr)*mg1
     . **2*upr)*mg2)*mg1+(tpr+upr)*mg1*mg2*tpr*u-(2.0*ms**2-s4s)*mg2
     . **2*tpr*u+((tpr+upr)*u-2.0*s*s4s)*mg1*mg2*s+(2.0*(tpr+upr+s3s)
     . *tpr-(s4s-2.0*tpr)*ms**2)*s*u+ans19
      ans17=(((u+upr)*upr-tpr**2)*s4s-2.0*(tpr+u+upr)*(tpr+upr)*tpr-
     . 2.0*s3s**2*tpr-(2.0*(u+2.0*upr+2.0*tpr)*tpr+(tpr-upr)*s4s)*s3s
     . +((3.0*(u+upr)+tpr)*s4s+2.0*(u+upr)*u+(3.0*s4s+2.0*u)*s3s)*s-(
     . 4.0*(tpr+upr)*tpr+(tpr-upr)*s4s+4.0*s3s*tpr-(3.0*s4s+2.0*u)*s+
     . 2.0*ms**2*tpr)*ms**2+2.0*((tpr+upr+s3s)*tpr-(s4s+u)*s+ms**2*
     . tpr)*mg2**2)*mg2**2-((((u+2.0*upr+tpr)*s4s+2.0*u**2+s3s*s4s)*
     . s3s+(2.0*u**2+u*upr+upr**2+tpr*upr)*s4s+2.0*(tpr+upr)*u**2)*s-
     . (((tpr*u-tpr*upr-u**2-upr**2)*s4s+2.0*(tpr+upr)*tpr*u)*(tpr+
     . upr)-(s4s*upr-2.0*tpr*u)*s3s**2+(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u)*s3s)-(((u-2.0*upr)*tpr-2.0*upr**2)*
     . s4s+4.0*(tpr+upr)*tpr*u-2.0*(s4s*upr-2.0*tpr*u)*s3s-((tpr+2.0*
     . upr)*s4s+2.0*u**2+2.0*s3s*s4s)*s-(s4s*upr-2.0*tpr*u+s*s4s)*ms
     . **2)*ms**2)+ans18
      ans16=ans17*ct6_1*s5x
      ans3=ans4+ans16
      ans20=gs**2
      ans2=2.0*ans3*ans20
      ans1=-ans2
      qg_su=ans1/(3.0*s*s4sc*s5*s5x*tpr*us)

      qg_tu=(-4.0*c4_tu*((((u-upr)*upr+tpr*u+s4s*upr-s3s*upr)*s3s-(((
     . u+upr)*upr+tpr*u)*s4s-tpr**2*upr)-((s4s+upr+s3s)*s3s+(s4s-tpr)
     . *upr)*s-(s+upr)*(s3s-s4s)*ms**2-((s4s+2.0*tpr-s3s)*upr-(s3s-
     . s4s)*s)*mg2**2+(s4s-2.0*tpr-s3s)*mg1**2*upr-(s3s-s4s)*mg2**2*
     . tpr)*s-((((tpr+2.0*upr)*s3s+tpr*upr)*upr-(tpr+upr-s)*mg2**2*
     . tpr+mg1**2*tpr*upr+(tpr+upr)*tpr*u-(ms**2-tpr)*tpr*upr-((tpr-
     . 2.0*upr)*s3s+ms**2*tpr)*s+2.0*tpr**2*u7den)*us+(((2.0*s4s+tpr)
     . *(tpr+upr)+s3s*tpr)*upr+(s3s*tpr+2.0*s4s*upr)*s+(s+upr)*ms**2*
     . tpr+mg2**2*tpr**2-(tpr+upr)*tpr*u-((s-upr)*mg2**2+mg1**2*upr)*
     . tpr)*u7den))*gs**2*mg1*mg2)/(3.0*s*s3sc*s4sc*tpr*u7den*us)
*
      qgreal=real(qg_s1+qg_t1+qg_u1+qg_st+qg_su+qg_tu)
ctp 
      LL_QGH = ( qgreal
     &         -qgd231_ll(massin,mkraemer))* (abs(mg1)+abs(mg2))**2/4.D0
*
      return
      end

************************************************************************
      real*8 function LL_QGHSUB(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/squarkw/lambda
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
      lambda = massin(30)
*
      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
*
      t = t2 + mg2**2
      u = u2 + mg2**2
      ts  = t - ms**2 
      us  = u - ms**2 
      tpr = - s - t - u1 + mg2**2
      upr = - s - u - t1 + mg2**2
      s3 = - s - t - u - tpr - upr + mg1**2 + mg2**2        
      s4 = s + t + u - mg1**2 - mg2**2
      s3s = s3 + mg2**2 - ms**2
      s4s = s4 + mg1**2 - ms**2
*
C --- introduce complex s4s, s3s to perform PV integration 
C --- s4s --> s4s + i*ms*Gamma
C --- s3s --> s3s + i*ms*Gamma
      Gamma = lambda*ms
      s4s2 = s4s**2 + ms**2*Gamma**2
      s3s2 = s3s**2 + ms**2*Gamma**2
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
C --- the form output for |M(q+g->g1+g2+q)|**2
C --- spin and colour average, colour factors included
      gs=1d0
C --- factor gs2 = 4*pi*alphas NOT included
      qg_s4=(-4.0*(mg1+ms)*(mg1-ms)*(mg2**2*ms**2*s-2.0*mg2**2*ms**2*
     . us-mg2**2*s**2-mg2**2*s*t+2.0*mg2**2*t*us+mg2**2*us**2+2.0*ms
     . **4*s+2.0*ms**4*us-3.0*ms**2*s**2-3.0*ms**2*s*t-3.0*ms**2*s*us
     . -2.0*ms**2*t*us+s**3+2.0*s**2*t+s**2*us+s*t**2+s*t*us-t*us**2)
     . *c4_u*gs**2)/(3.0*s*s4s2*us**2)
      qg_s3=(-4.0*(((t+tpr+2.0*ms**2)*(t+tpr)+mg2**4+2.0*ms**4-2.0*(t
     . +tpr+ms**2)*mg2**2)*(t+tpr+s-mg2**2)+2.0*(t+tpr-mg2**2)*mg1**4
     . -2.0*((t+tpr+s)*(t+tpr)+mg2**4+(2.0*(t+tpr)+s)*ms**2-(2.0*(t+
     . tpr)+s+2.0*ms**2)*mg2**2)*mg1**2)*(mg2+ms)*(mg2-ms)*c4_t*gs**2
     . )/(3.0*(t+tpr-mg2**2)**2*s*s3s2)
*
      qgres=qg_s4+qg_s3
ctp 
      LL_QGHSUB = qgres * (abs(mg1)+abs(mg2))**2/4.D0
*
      return
      end
************************************************************************
      real*8 function qgd231_ll(massin,mkraemer)
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
C --- the form output for D132
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
      v231=((2.0*(tpr+upr+s)*(tpr+upr)+s**2)*gs**2)/s**2
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
      qgd231_ll=d231
*
      return
      end

************************************************************************
      real*8 function LL_GBH(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
      complex*16 s4sc,s3sc,qbg_s1,qbg_t1,qbg_u1,qbg_st,qbg_su,qbg_tu 
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/squarkw/lambda
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
      lambda = massin(30)
*
      mg12=mg1**2
      ms2=ms**2
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
      s3 = - s - t - u - tpr - upr + mg1**2 + mg2**2        
      s4 = s + t + u - mg1**2 - mg2**2
      s3s = s3 + mg2**2 - ms**2
      s4s = s4 + mg1**2 - ms**2
*
      s5 = s + tpr + upr
      u6 = - s - t - tpr + mg2**2
      u7 = - s - u - upr + mg2**2
*
      u7den=u7+mg1**2-ms**2
      u6den=u6+mg1**2-ms**2
      s5x=s5-mx2
C --- introduce complex s4s, s3s to perform PV integration 
C --- s4s --> s4s + i*ms*Gamma
C --- s3s --> s3s + i*ms*Gamma
      Gamma = lambda*ms
      s4sc = s4s + ms*Gamma*dcmplx(0d0,1d0)
      s3sc = s3s - ms*Gamma*dcmplx(0d0,1d0)
      s4s2 = s4s**2 + ms**2*Gamma**2
      s3s2 = s3s**2 + ms**2*Gamma**2
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
C --- the form output for |M(q+g->g1+g2+q)|**2
C --- spin and colour average, colour factors included
      gs=1d0
C --- factor gs2 = 4*pi*alphas NOT included
      ans5=(((2.0*tpr**2+2.0*tpr*upr+upr**2+2.0*t*tpr+2.0*(tpr+upr)*
     . s3s+(2.0*(t+tpr)+s)*s+2.0*(tpr+upr)*ms**2-2.0*(2.0*tpr+upr+s)*
     . mg2**2)*mg1-2.0*((s+2.0*tpr)*s+2.0*tpr**2+2.0*tpr*upr+upr**2)*
     . mg2)*mg1+((2.0*(t+2.0*tpr+2.0*s3s)+s)*s+2.0*tpr**2+4.0*tpr*upr
     . +upr**2+2.0*(tpr+2.0*upr)*t+2.0*(tpr+upr)*s3s+2.0*(tpr+upr+2.0
     . *s)*ms**2)*mg2**2-((2.0*ms**4+4.0*ms**2*s3s+2.0*ms**2*tpr+2.0*
     . ms**2*upr+s**2+2.0*s*t+2.0*s*tpr+2.0*s3s**2+2.0*s3s*tpr+2.0*
     . s3s*upr+2.0*t**2+2.0*t*tpr+2.0*tpr**2+2.0*tpr*upr+upr**2)*(s+
     . tpr+upr)+2.0*mg2**4*s)-2.0*(mg2**2-2.0*s)*mg2**2*upr)*c2_22-((
     . 3.0*tpr+upr+2.0*t+s)*s**2-((2.0*tpr+upr)*upr-2.0*t*tpr+2.0*(
     . tpr+upr)*s3s)*(tpr+upr)+(2.0*(2.0*tpr+upr)*t-upr**2-2.0*(tpr+
     . upr)*s3s)*s-2.0*((tpr+upr)**2+s*upr)*ms**2-(s+upr)*(s-upr)*mg2
     . **2-(s+2.0*tpr)*mg1**2*s+(2.0*tpr+upr)*mg1**2*upr-2.0*(ms**2-
     . tpr)*s*tpr)*c3_22
      ans4=(((2.0*tpr**2+2.0*tpr*upr+upr**2+2.0*t*tpr+2.0*(tpr+upr)*
     . s3s+(2.0*(t+tpr)+s)*s+2.0*(tpr+upr)*ms**2-2.0*(2.0*tpr+upr+s)*
     . mg2**2)*mg1+2.0*((s+2.0*tpr)*s+2.0*tpr**2+2.0*tpr*upr+upr**2)*
     . mg2)*mg1+((2.0*(t+2.0*tpr+2.0*s3s)+s)*s+2.0*tpr**2+4.0*tpr*upr
     . +upr**2+2.0*(tpr+2.0*upr)*t+2.0*(tpr+upr)*s3s+2.0*(tpr+upr+2.0
     . *s)*ms**2)*mg2**2-((2.0*ms**4+4.0*ms**2*s3s+2.0*ms**2*tpr+2.0*
     . ms**2*upr+s**2+2.0*s*t+2.0*s*tpr+2.0*s3s**2+2.0*s3s*tpr+2.0*
     . s3s*upr+2.0*t**2+2.0*t*tpr+2.0*tpr**2+2.0*tpr*upr+upr**2)*(s+
     . tpr+upr)+2.0*mg2**4*s)-2.0*(mg2**2-2.0*s)*mg2**2*upr)*c1_22+
     . ans5
      ans6=s5**2
      ans3=ans4*ans6
      ans7=-(2.0*((3.0*tpr+upr+2.0*t+s)*s**2-((2.0*tpr+upr)*upr-2.0*t
     . *tpr+2.0*(tpr+upr)*s3s)*(tpr+upr)+(2.0*(2.0*tpr+upr)*t-upr**2-
     . 2.0*(tpr+upr)*s3s)*s-2.0*((tpr+upr)**2+s*upr)*ms**2-(s+upr)*(s
     . -upr)*mg2**2-(s+2.0*tpr)*mg1**2*s+(2.0*tpr+upr)*mg1**2*upr-2.0
     . *(ms**2-tpr)*s*tpr)*c3_12*s5-(((2.0*tpr**2+2.0*tpr*upr+upr**2+
     . 2.0*t*tpr+2.0*(tpr+upr)*s3s+(2.0*(t+tpr)+s)*s+2.0*(tpr+upr)*ms
     . **2-2.0*(2.0*tpr+upr+s)*mg2**2)*mg1-2.0*((s+2.0*tpr)*s+2.0*tpr
     . **2+2.0*tpr*upr+upr**2)*mg2)*mg1+((2.0*(t+2.0*tpr+2.0*s3s)+s)*
     . s+2.0*tpr**2+4.0*tpr*upr+upr**2+2.0*(tpr+2.0*upr)*t+2.0*(tpr+
     . upr)*s3s+2.0*(tpr+upr+2.0*s)*ms**2)*mg2**2-((2.0*ms**4+4.0*ms
     . **2*s3s+2.0*ms**2*tpr+2.0*ms**2*upr+s**2+2.0*s*t+2.0*s*tpr+2.0
     . *s3s**2+2.0*s3s*tpr+2.0*s3s*upr+2.0*t**2+2.0*t*tpr+2.0*tpr**2+
     . 2.0*tpr*upr+upr**2)*(s+tpr+upr)+2.0*mg2**4*s)-2.0*(mg2**2-2.0*
     . s)*mg2**2*upr)*(c2_11*s5x+2.0*c2_12*s5))*s5x
      ans2=ans3+ans7
      ans8=gs**2
      ans1=2.0*ans2*ans8
      qbg_s1=ans1/(3.0*s*s5**2*s5x**2*upr)

      ans4=-((2.0*mg1**2*mg2**2*s4s-mg1**2*mg2**2*upr-2.0*mg1**2*ms**
     . 2*s4s+mg1**2*ms**2*upr-2.0*mg1**2*s3s*s4s+mg1**2*s3s*upr-mg1**
     . 2*s4s**2-2.0*mg1**2*s4s*tpr+mg1**2*tpr*upr+mg1**2*upr**2-mg2**
     . 4*upr+2.0*mg2**2*ms**2*upr+2.0*mg2**2*s3s*upr-mg2**2*s4s**2-
     . 2.0*mg2**2*s4s*t+mg2**2*s4s*upr-mg2**2*t*upr+2.0*mg2**2*tpr*
     . upr+2.0*mg2**2*upr**2-ms**4*upr-2.0*ms**2*s3s*upr+ms**2*s4s**2
     . +2.0*ms**2*s4s*t-ms**2*s4s*upr+ms**2*t*upr-2.0*ms**2*tpr*upr-
     . 2.0*ms**2*upr**2-s3s**2*upr+s3s*s4s**2+2.0*s3s*s4s*t-s3s*s4s*
     . upr+s3s*t*upr-2.0*s3s*tpr*upr-2.0*s3s*upr**2+s4s**2*t+s4s**2*
     . tpr+2.0*s4s*t*tpr+2.0*s4s*t*upr-s4s*tpr*upr-s4s*upr**2+t*tpr*
     . upr+t*upr**2-tpr**2*upr-2.0*tpr*upr**2-upr**3)*(mg2**2-t)*s-(
     . mg1**2+mg2**2-ms**2-s-s3s-t-tpr-upr)*(mg2**2-ms**2-s3s-tpr-upr
     . )*ts**2*upr)
      ans3=(((2.0*(2.0*(tpr+upr)+t)*upr-(tpr-upr)*s4s+4.0*s3s*upr-(
     . s4s+upr)*s+2.0*ms**2*upr)*ms**2+(2.0*(tpr+upr+t)*upr-(tpr-upr+
     . t)*s4s)*(tpr+upr)+2.0*s3s**2*upr+(2.0*(2.0*(tpr+upr)+t)*upr-(
     . tpr-upr)*s4s)*s3s-(tpr+upr-t+s3s)*(s4s+upr)*s)*mg2**2-((((tpr+
     . upr)*upr+s4s*tpr+(s4s+upr)*s3s)*s+(2.0*s3s+s4s+2.0*tpr+2.0*upr
     . )*(s3s+tpr+upr)*upr+(4.0*(tpr+upr+s3s)*upr+s*s4s+2.0*ms**2*upr
     . )*ms**2)*t+2.0*(tpr+upr+s3s)*mg2**4*upr)-((tpr+upr+t)*s4s*tpr-
     . 2.0*(tpr+upr)*t*upr+(s4s*tpr-2.0*t*upr)*s3s+(s4s*tpr-2.0*t*upr
     . )*ms**2+2.0*((tpr+upr+t)*upr-s4s*tpr+s3s*upr)*mg2**2)*mg1**2-(
     . 2.0*ms**2-s)*mg2**4*upr+2.0*(mg2+ms)*(mg2-ms)*mg1**2*mg2**2*
     . upr+((2.0*(tpr+upr)*tpr-t*upr+2.0*s3s*tpr)*ms**2+((tpr+upr)*
     . tpr+t**2)*(tpr+upr)+(2.0*(tpr+upr)+s3s)*s3s*tpr+ms**4*tpr)*s4s
     . +(((2.0*tpr+upr)*s4s-t*upr+2.0*s3s*s4s)*ms**2+(s3s+tpr+upr)*(
     . s3s+tpr)*s4s+ms**4*s4s)*s)*ts+ans4
      ans5=gs**2
      ans2=4.0*c4_t*ans3*ans5
      ans1=-ans2
      qbg_t1=ans1/(3.0*s*s4s2*ts**2*upr)


      ans4=(2.0*mg1**2*mg2**2*s3s+mg1**2*mg2**2*upr-2.0*mg1**2*ms**2*
     . s3s-mg1**2*ms**2*upr-mg1**2*s3s**2-mg1**2*s3s*upr+mg2**4*upr-
     . 2.0*mg2**2*s*s3s-mg2**2*s*upr+mg2**2*s3s**2-2.0*mg2**2*s3s*t-
     . 2.0*mg2**2*s3s*tpr+mg2**2*s3s*upr-mg2**2*t*upr-mg2**2*tpr*upr-
     . ms**4*upr+2.0*ms**2*s*s3s+ms**2*s*upr-ms**2*s3s**2+2.0*ms**2*
     . s3s*t+2.0*ms**2*s3s*tpr-ms**2*s3s*upr+ms**2*t*upr+ms**2*tpr*
     . upr+s*s3s**2+s*s3s*upr-s3s**3+s3s**2*t+s3s**2*tpr-s3s**2*upr+
     . s3s*t*upr+s3s*tpr*upr)*(mg2**2-s-t-tpr)*s+(mg1**2+mg2**2-ms**2
     . -s3s-t-tpr)*(mg2**2-ms**2-s3s)*u6den**2*upr
      ans3=(((3.0*(t+tpr)*upr+2.0*s3s*tpr)*s3s-(s3s+upr)*s**2-((t+tpr
     . )*s3s+t*upr)*s+(2.0*(t+tpr)*upr+s3s*tpr)*ms**2)*ms**2+((t**2+
     . tpr**2)*(tpr+upr)+2.0*t*tpr**2+((t+tpr)*upr+s3s*tpr)*s3s-((s3s
     . -tpr)*s+(t+tpr)*s3s)*s)*s3s-(((t+tpr)*(tpr+3.0*upr)+s3s*tpr)*
     . s3s-s**2*upr+((s3s-t)*s3s-(t+tpr)*upr)*s+(2.0*(t+tpr)*upr+(tpr
     . +3.0*upr)*s3s+(s3s+upr)*s+2.0*ms**2*upr)*ms**2)*mg2**2-(((t+
     . tpr)*(tpr+2.0*upr)+s3s*tpr+s*tpr)*s3s+(2.0*(t+tpr)*upr+s3s*tpr
     . )*ms**2-2.0*((tpr+upr)*s3s+tpr*upr+ms**2*upr-mg2**2*upr)*mg2**
     . 2)*mg1**2-(s*upr-2.0*t*tpr)*s*s3s+(s**2-s3s*upr)*mg2**2*s3s-(
     . mg2**4*s-2.0*s3s*t*tpr)*upr+2.0*(mg1**2*mg2**2*t+s*s3s**2)*upr
     . +((2.0*s3s**2+4.0*s3s*upr-tpr*upr)*ms**2+(s3s**2+2.0*tpr**2)*
     . s3s+(s3s+2.0*upr)*ms**4)*s+2.0*((ms**2+s3s)*mg2**2-s*s3s)*mg2
     . **2*upr)*u6den+ans4
      ans5=gs**2
      ans2=4.0*c4_u*ans3*ans5
      ans1=-ans2
      qbg_u1=ans1/(3.0*s*s3s2*u6den**2*upr)

      ans7=-(((tpr+upr)*tpr+t**2)*(tpr+upr)+(2.0*(tpr+upr)+s3s)*s3s*
     . tpr+s*t*tpr+(2.0*(tpr+upr)*tpr-t*upr+2.0*s3s*tpr)*ms**2+ms**4*
     . tpr)*s4s-(((2.0*tpr+upr)*s4s+2.0*t**2+s3s*s4s)*s3s+((tpr+upr)*
     . tpr+2.0*t**2)*s4s+2.0*(tpr+upr)*t**2+(2.0*(t-upr)*t+(2.0*tpr+
     . upr)*s4s+2.0*s3s*s4s)*ms**2+ms**4*s4s)*s
      ans6=((2.0*((tpr+upr+t)*upr-s4s*tpr+s3s*upr-(t+tpr+s4s+s3s)*s-(
     . s-upr)*ms**2+(s-upr)*mg2**2)*mg2**2+(tpr+upr+t)*s4s*tpr-2.0*(
     . tpr+upr)*t*upr+(s4s*tpr-2.0*t*upr)*s3s+2.0*(s4s+tpr+s3s)*s*t+(
     . s4s*tpr-2.0*t*upr+2.0*s*t)*ms**2)*mg1+(((t-tpr)*upr-2.0*s4s*
     . tpr+s*t)*tpr-(2.0*s4s+tpr+s3s)*s**2-ms**2*s**2+(s+tpr+upr)*(s-
     . upr)*mg2**2-(s-upr)*mg1**2*tpr)*mg2)*mg1-(((2.0*(2.0*(tpr+upr)
     . +t)*upr-(tpr-upr)*s4s+4.0*s3s*upr-(3.0*s4s+2.0*t)*s+2.0*ms**2*
     . upr)*ms**2+(2.0*(tpr+upr+t)*upr-(tpr-upr+t)*s4s)*(tpr+upr)+2.0
     . *s3s**2*upr+(2.0*(2.0*(tpr+upr)+t)*upr-(tpr-upr)*s4s)*s3s-((
     . 3.0*tpr+upr+3.0*t)*s4s+2.0*(t+tpr)*t+(3.0*s4s+2.0*t)*s3s)*s)*
     . mg2**2-(2.0*((tpr+upr+s3s)*upr-(s4s+t)*s+ms**2*upr)*mg2**4+(
     . 2.0*(2.0*(tpr+upr+s3s)+ms**2)*ms**2+(2.0*s3s+s4s+2.0*tpr+2.0*
     . upr)*(s3s+tpr+upr)+2.0*(s3s+upr)*s)*t*upr))-(s3s*s4s-2.0*tpr*
     . upr+ms**2*s4s)*s*t+((t-tpr)*upr-2.0*s4s*tpr-s3s*tpr-ms**2*tpr)
     . *mg1*mg2*upr+(t*upr-tpr**2-2.0*s4s*tpr-(tpr+upr)*s3s-(tpr+upr)
     . *ms**2)*mg1*mg2*s+ans7
      ans5=ans6*c6_2
      ans10=-(((tpr+2.0*upr)*(tpr+upr)+t**2)*(tpr+upr)+(2.0*(tpr+upr)
     . +s3s)*(tpr+2.0*upr)*s3s-(s3s+tpr)*s*t+(2.0*(tpr+upr+s3s)+ms**2
     . )*(tpr+2.0*upr)*ms**2-((tpr+upr+t+s3s)*(tpr+3.0*upr)-(t-tpr-
     . s3s)*s+(tpr+3.0*upr+s)*ms**2)*mg2**2-(((tpr+2.0*upr)*(tpr+upr)
     . +t*tpr+(tpr+2.0*upr)*s3s+(tpr+2.0*upr)*ms**2-2.0*(tpr+upr)*mg2
     . **2)*mg1-2.0*(tpr**2+tpr*upr+upr**2+s*tpr)*mg2)*mg1+((tpr+upr+
     . s3s+ms**2)*t+2.0*mg2**4-3.0*mg2**2*s)*upr+((2.0*tpr+3.0*upr-t+
     . 2.0*s3s)*ms**2+(s3s+tpr+2.0*upr)*(s3s+tpr+upr)+ms**4)*s)*c6_2
      ans9=(((2.0*tpr+3.0*upr-t+2.0*s3s)*s+2.0*(tpr+2.0*upr)*(tpr+upr
     . )+t*upr+2.0*(tpr+2.0*upr)*s3s+(tpr+2.0*upr+s)*ms**2)*ms**2+(
     . 2.0*(tpr+2.0*upr)*(tpr+upr)+t*upr+(tpr+2.0*upr)*s3s)*s3s+((t+
     . upr)*t+(tpr+2.0*upr)*(tpr+upr))*(tpr+upr)+((tpr+2.0*upr)*(tpr+
     . upr)-t*tpr+s3s**2+(2.0*tpr+3.0*upr-t)*s3s)*s-((tpr+upr+t+s3s+
     . ms**2)*(tpr+3.0*upr)-2.0*mg2**2*upr)*mg2**2-(((tpr+2.0*upr)*(
     . tpr+upr)+t*tpr+(tpr+2.0*upr)*s3s+(tpr+2.0*upr)*ms**2-2.0*(tpr+
     . upr)*mg2**2)*mg1+2.0*(tpr**2+tpr*upr+upr**2+s*tpr)*mg2)*mg1-(
     . tpr+3.0*upr-t+s3s+ms**2)*mg2**2*s)*c5_2+ans10
      ans8=ans9*ts
      ans13=-((2.0*((tpr+upr+t)*upr-s4s*tpr+s3s*upr-(t+tpr+s4s+s3s)*s
     . -(s-upr)*ms**2+(s-upr)*mg2**2)*mg2**2+(tpr+upr+t)*s4s*tpr-2.0*
     . (tpr+upr)*t*upr+(s4s*tpr-2.0*t*upr)*s3s+2.0*(s4s+tpr+s3s)*s*t+
     . (s4s*tpr-2.0*t*upr+2.0*s*t)*ms**2)*mg1-(((t-tpr)*upr-2.0*s4s*
     . tpr+s*t)*tpr-(2.0*s4s+tpr+s3s)*s**2-ms**2*s**2+(s+tpr+upr)*(s-
     . upr)*mg2**2-(s-upr)*mg1**2*tpr)*mg2)*mg1+(s4s-2.0*upr)*s*t*tpr
     . +(2.0*ms**2-s4s)*mg2**2*t*upr-(2.0*(s3s+upr)*upr-(s4s-2.0*upr)
     . *ms**2)*s*t+((t-tpr)*upr-2.0*s4s*tpr-s3s*tpr-ms**2*tpr)*mg1*
     . mg2*upr+(t*upr-tpr**2-2.0*s4s*tpr-(tpr+upr)*s3s-(tpr+upr)*ms**
     . 2)*mg1*mg2*s
      ans12=((2.0*(2.0*(tpr+upr)+t)*upr-(tpr-upr)*s4s+2.0*s3s*upr)*
     . s3s-(((tpr+upr)*(tpr-upr)+t*tpr)*s4s-2.0*(tpr+upr+t)*(tpr+upr)
     . *upr)-((3.0*tpr+upr+3.0*t)*s4s+2.0*(t+tpr)*t+(3.0*s4s+2.0*t)*
     . s3s)*s+(4.0*(tpr+upr)*upr-(tpr-upr)*s4s+4.0*s3s*upr-(3.0*s4s+
     . 2.0*t)*s+2.0*ms**2*upr)*ms**2-2.0*((tpr+upr+s3s)*upr-(s4s+t)*s
     . +ms**2*upr)*mg2**2)*mg2**2+((2.0*(tpr+upr)*tpr-t*upr)*s4s-4.0*
     . (tpr+upr)*t*upr+2.0*(s4s*tpr-2.0*t*upr)*s3s+((2.0*tpr+upr)*s4s
     . +2.0*t**2+2.0*s3s*s4s)*s+(s4s*tpr-2.0*t*upr+s*s4s)*ms**2)*ms**
     . 2+((2.0*(tpr+upr)*tpr-t*upr)*s4s-4.0*(tpr+upr)*t*upr+(s4s*tpr-
     . 2.0*t*upr)*s3s)*s3s+(((t-upr)*t+(tpr+upr)*tpr)*s4s-2.0*(tpr+
     . upr)*t*upr)*(tpr+upr)+(((tpr+upr)*tpr+2.0*t**2)*s4s+2.0*(tpr+
     . upr)*t**2+s3s**2*s4s+((2.0*tpr+upr+t)*s4s+2.0*t**2)*s3s)*s+
     . ans13
      ans11=ans12*c5_2
      ans4=ans5+ans8+ans11
      ans3=ans4*s5
      ans16=-(((tpr+upr)*tpr+t**2)*(tpr+upr)+(2.0*(tpr+upr)+s3s)*s3s*
     . tpr+s*t*tpr+(2.0*(tpr+upr)*tpr-t*upr+2.0*s3s*tpr)*ms**2+ms**4*
     . tpr)*s4s-(((2.0*tpr+upr)*s4s+2.0*t**2+s3s*s4s)*s3s+((tpr+upr)*
     . tpr+2.0*t**2)*s4s+2.0*(tpr+upr)*t**2+(2.0*(t-upr)*t+(2.0*tpr+
     . upr)*s4s+2.0*s3s*s4s)*ms**2+ms**4*s4s)*s-(((tpr+2.0*upr)*(tpr+
     . upr)+t**2)*(tpr+upr)+(2.0*(tpr+upr)+s3s)*(tpr+2.0*upr)*s3s-(
     . s3s+tpr)*s*t+(2.0*(tpr+upr+s3s)+ms**2)*(tpr+2.0*upr)*ms**2-((
     . tpr+upr+t+s3s)*(tpr+3.0*upr)-(t-tpr-s3s)*s+(tpr+3.0*upr+s)*ms
     . **2)*mg2**2-(((tpr+2.0*upr)*(tpr+upr)+t*tpr+(tpr+2.0*upr)*s3s+
     . (tpr+2.0*upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+tpr*
     . upr+upr**2+s*tpr)*mg2)*mg1+((tpr+upr+s3s+ms**2)*t+2.0*mg2**4-
     . 3.0*mg2**2*s)*upr+((2.0*tpr+3.0*upr-t+2.0*s3s)*ms**2+(s3s+tpr+
     . 2.0*upr)*(s3s+tpr+upr)+ms**4)*s)*ts
      ans15=((2.0*((tpr+upr+t)*upr-s4s*tpr+s3s*upr-(t+tpr+s4s+s3s)*s-
     . (s-upr)*ms**2+(s-upr)*mg2**2)*mg2**2+(tpr+upr+t)*s4s*tpr-2.0*(
     . tpr+upr)*t*upr+(s4s*tpr-2.0*t*upr)*s3s+2.0*(s4s+tpr+s3s)*s*t+(
     . s4s*tpr-2.0*t*upr+2.0*s*t)*ms**2)*mg1+(((t-tpr)*upr-2.0*s4s*
     . tpr+s*t)*tpr-(2.0*s4s+tpr+s3s)*s**2-ms**2*s**2+(s+tpr+upr)*(s-
     . upr)*mg2**2-(s-upr)*mg1**2*tpr)*mg2)*mg1-(((2.0*(2.0*(tpr+upr)
     . +t)*upr-(tpr-upr)*s4s+4.0*s3s*upr-(3.0*s4s+2.0*t)*s+2.0*ms**2*
     . upr)*ms**2+(2.0*(tpr+upr+t)*upr-(tpr-upr+t)*s4s)*(tpr+upr)+2.0
     . *s3s**2*upr+(2.0*(2.0*(tpr+upr)+t)*upr-(tpr-upr)*s4s)*s3s-((
     . 3.0*tpr+upr+3.0*t)*s4s+2.0*(t+tpr)*t+(3.0*s4s+2.0*t)*s3s)*s)*
     . mg2**2-(2.0*((tpr+upr+s3s)*upr-(s4s+t)*s+ms**2*upr)*mg2**4+(
     . 2.0*(2.0*(tpr+upr+s3s)+ms**2)*ms**2+(2.0*s3s+s4s+2.0*tpr+2.0*
     . upr)*(s3s+tpr+upr)+2.0*(s3s+upr)*s)*t*upr))-(s3s*s4s-2.0*tpr*
     . upr+ms**2*s4s)*s*t+((t-tpr)*upr-2.0*s4s*tpr-s3s*tpr-ms**2*tpr)
     . *mg1*mg2*upr+(t*upr-tpr**2-2.0*s4s*tpr-(tpr+upr)*s3s-(tpr+upr)
     . *ms**2)*mg1*mg2*s+ans16
      ans14=ans15*c6_1*s5x
      ans2=ans3+ans14
      ans17=gs**2
      ans1=2.0*ans2*ans17
      qbg_st=ans1/(3.0*s*s4sc*s5*s5x*ts*upr)


      ans7=-(s3s+tpr-ms**2)*mg1*mg2*s*tpr-(3.0*s3s**2+3.0*s3s*upr-2.0
     . *tpr**2+3.0*ms**2*s3s)*mg2**2*s-(2.0*(s**2+2.0*t*tpr)-(s3s-2.0
     . *upr)*s)*ms**2*s-((tpr+upr)*t+s*s3s-ms**2*s)*mg1*mg2*s
      ans6=((2.0*((t+tpr)*upr+(tpr+upr)*s3s-(t+tpr+s)*s-(s-upr)*ms**2
     . +(s-upr)*mg2**2)*mg2**2-((2.0*(t+tpr)*upr+s3s*tpr-2.0*(s+t)*s)
     . *ms**2+((t+tpr)*(tpr+2.0*upr)+s3s*tpr)*s3s))*mg1-((s*tpr-s3s*
     . upr)*s+(2.0*s3s+upr)*tpr**2-ms**2*s*upr+(s+tpr+upr)*(s-upr)*
     . mg2**2-(s-upr)*mg1**2*tpr)*mg2)*mg1+(2.0*((s3s+2.0*upr)*s3s-
     . 2.0*(t+tpr)*s)*s+(3.0*(t+tpr)*upr+2.0*s3s*tpr)*s3s+(2.0*(t+tpr
     . )*upr+s3s*tpr+(s3s+2.0*upr)*s)*ms**2)*ms**2+((t+tpr)**2*(tpr+
     . upr)+(t+tpr)*s3s*upr+s3s**2*tpr+(s3s+tpr)*s**2+(2.0*(t+tpr)*
     . tpr+s3s**2+(tpr+2.0*upr+t)*s3s)*s)*s3s-(((t+tpr)*(tpr+3.0*upr)
     . +(tpr+upr)*s3s)*s3s-2.0*(s+t+2.0*tpr)*(s+t)*s+(2.0*(t+tpr)*upr
     . +(tpr+3.0*upr)*s3s+2.0*ms**2*upr)*ms**2+2.0*((t+tpr-s3s)*s+s**
     . 2-s3s*upr-ms**2*upr)*mg2**2)*mg2**2+(s3s-2.0*upr)*(t+tpr)*ms**
     . 2*s+(mg1*upr+2.0*mg2*s)*(ms**2-s3s)*mg2*tpr+((s+t)*mg2**2+mg1
     . **2*tpr)*(2.0*ms**2-s3s)*s-(2.0*ms**2*tpr-s3s*upr)*s*tpr-(2.0*
     . ms**2*t-s3s*upr)*s*t-((tpr+upr)*t+tpr*upr)*mg1*mg2*upr+ans7
      ans5=ans6*ct5_2
      ans11=((((t+tpr)*tpr+(tpr+2.0*upr)*s3s+s*tpr+(tpr+2.0*upr)*ms**
     . 2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+tpr*upr+upr**2+s*tpr)*
     . mg2)*mg1-(((t+tpr)*upr+2.0*(tpr+2.0*upr)*s3s+2.0*s*s3s+(tpr+
     . 2.0*upr+s)*ms**2)*ms**2+(t+tpr)**2*(tpr+upr)+(t+tpr)*s3s*upr+(
     . tpr+2.0*upr)*s3s**2+((2.0*tpr+upr)*tpr+s3s**2)*s-((t+tpr+s3s)*
     . (tpr+3.0*upr)-s**2+(tpr+3.0*upr)*ms**2-2.0*mg2**2*upr)*mg2**2)
     . +(s3s-2.0*tpr+ms**2)*s*t+(s3s-tpr+ms**2)*s**2+(s3s*tpr-t*upr+
     . ms**2*tpr)*s-(t-upr-s3s-ms**2)*mg2**2*s)*ct6_2
      ans10=((((t+tpr)*tpr+(tpr+2.0*upr)*s3s+s*tpr+(tpr+2.0*upr)*ms**
     . 2-2.0*(tpr+upr)*mg2**2)*mg1+2.0*(tpr**2+tpr*upr+upr**2+s*tpr)*
     . mg2)*mg1-(((t+tpr)*upr+2.0*(tpr+2.0*upr)*s3s+2.0*s*s3s+(tpr+
     . 2.0*upr+s)*ms**2)*ms**2+(t+tpr)**2*(tpr+upr)+(t+tpr)*s3s*upr+(
     . tpr+2.0*upr)*s3s**2+((2.0*tpr+upr)*tpr+s3s**2)*s-((t+tpr+s3s)*
     . (tpr+3.0*upr)-s**2+(tpr+3.0*upr)*ms**2-2.0*mg2**2*upr)*mg2**2)
     . +(s3s-2.0*tpr+ms**2)*s*t+(s3s-tpr+ms**2)*s**2+(s3s*tpr-t*upr+
     . ms**2*tpr)*s-(t-upr-s3s-ms**2)*mg2**2*s)*ct5_2+ans11
      ans9=ans10*u6den
      ans8=-ans9
      ans14=-(s3s*upr-tpr**2+ms**2*upr)*mg1*mg2*s-(3.0*s3s**2+3.0*s3s
     . *upr-2.0*tpr**2+3.0*ms**2*s3s)*mg2**2*s-(2.0*(s**2+2.0*t*tpr)-
     . (s3s-2.0*upr)*s)*ms**2*s+((tpr+upr)*t+tpr*upr+s3s*tpr-ms**2*
     . tpr)*mg1*mg2*upr
      ans13=((2.0*((t+tpr)*upr+(tpr+upr)*s3s-(t+tpr+s)*s-(s-upr)*ms**
     . 2+(s-upr)*mg2**2)*mg2**2-((2.0*(t+tpr)*upr+s3s*tpr-2.0*(s+t)*s
     . )*ms**2+((t+tpr)*(tpr+2.0*upr)+s3s*tpr)*s3s))*mg1+(((s3s+t)*
     . tpr+s*s3s)*s+(2.0*s3s+upr)*tpr**2+(s+tpr+upr)*(s-upr)*mg2**2-(
     . s-upr)*mg1**2*tpr)*mg2)*mg1+(2.0*((s3s+2.0*upr)*s3s-2.0*(t+tpr
     . )*s)*s+(3.0*(t+tpr)*upr+2.0*s3s*tpr)*s3s+(2.0*(t+tpr)*upr+s3s*
     . tpr+(s3s+2.0*upr)*s)*ms**2)*ms**2+((t+tpr)**2*(tpr+upr)+(t+tpr
     . )*s3s*upr+s3s**2*tpr+(s3s+tpr)*s**2+(2.0*(t+tpr)*tpr+s3s**2+(
     . tpr+2.0*upr+t)*s3s)*s)*s3s-(((t+tpr)*(tpr+3.0*upr)+(tpr+upr)*
     . s3s)*s3s-2.0*(s+t+2.0*tpr)*(s+t)*s+(2.0*(t+tpr)*upr+(tpr+3.0*
     . upr)*s3s+2.0*ms**2*upr)*ms**2+2.0*((t+tpr-s3s)*s+s**2-s3s*upr-
     . ms**2*upr)*mg2**2)*mg2**2+(s3s-2.0*upr)*(t+tpr)*ms**2*s-(ms**2
     . -tpr)*mg1*mg2*s**2+2.0*(ms**2-s3s)*mg2**2*s*tpr+((s+t)*mg2**2+
     . mg1**2*tpr)*(2.0*ms**2-s3s)*s-(ms**2*tpr-t*upr)*mg1*mg2*s-(2.0
     . *ms**2*tpr-s3s*upr)*s*tpr-(2.0*ms**2*t-s3s*upr)*s*t+ans14
      ans12=ans13*ct6_2
      ans4=ans5+ans8+ans12
      ans3=ans4*s5
      ans17=-(s3s*upr-tpr**2+ms**2*upr)*mg1*mg2*s-(3.0*s3s**2+3.0*s3s
     . *upr-2.0*tpr**2+3.0*ms**2*s3s)*mg2**2*s-(2.0*(s**2+2.0*t*tpr)-
     . (s3s-2.0*upr)*s)*ms**2*s+((tpr+upr)*t+tpr*upr+s3s*tpr-ms**2*
     . tpr)*mg1*mg2*upr-((((t+tpr)*tpr+(tpr+2.0*upr)*s3s+s*tpr+(tpr+
     . 2.0*upr)*ms**2-2.0*(tpr+upr)*mg2**2)*mg1-2.0*(tpr**2+tpr*upr+
     . upr**2+s*tpr)*mg2)*mg1-(((t+tpr)*upr+2.0*(tpr+2.0*upr)*s3s+2.0
     . *s*s3s+(tpr+2.0*upr+s)*ms**2)*ms**2+(t+tpr)**2*(tpr+upr)+(t+
     . tpr)*s3s*upr+(tpr+2.0*upr)*s3s**2+((2.0*tpr+upr)*tpr+s3s**2)*s
     . -((t+tpr+s3s)*(tpr+3.0*upr)-s**2+(tpr+3.0*upr)*ms**2-2.0*mg2**
     . 2*upr)*mg2**2)+(s3s-2.0*tpr+ms**2)*s*t+(s3s-tpr+ms**2)*s**2+(
     . s3s*tpr-t*upr+ms**2*tpr)*s-(t-upr-s3s-ms**2)*mg2**2*s)*u6den
      ans16=((2.0*((t+tpr)*upr+(tpr+upr)*s3s-(t+tpr+s)*s-(s-upr)*ms**
     . 2+(s-upr)*mg2**2)*mg2**2-((2.0*(t+tpr)*upr+s3s*tpr-2.0*(s+t)*s
     . )*ms**2+((t+tpr)*(tpr+2.0*upr)+s3s*tpr)*s3s))*mg1+(((s3s+t)*
     . tpr+s*s3s)*s+(2.0*s3s+upr)*tpr**2+(s+tpr+upr)*(s-upr)*mg2**2-(
     . s-upr)*mg1**2*tpr)*mg2)*mg1+(2.0*((s3s+2.0*upr)*s3s-2.0*(t+tpr
     . )*s)*s+(3.0*(t+tpr)*upr+2.0*s3s*tpr)*s3s+(2.0*(t+tpr)*upr+s3s*
     . tpr+(s3s+2.0*upr)*s)*ms**2)*ms**2+((t+tpr)**2*(tpr+upr)+(t+tpr
     . )*s3s*upr+s3s**2*tpr+(s3s+tpr)*s**2+(2.0*(t+tpr)*tpr+s3s**2+(
     . tpr+2.0*upr+t)*s3s)*s)*s3s-(((t+tpr)*(tpr+3.0*upr)+(tpr+upr)*
     . s3s)*s3s-2.0*(s+t+2.0*tpr)*(s+t)*s+(2.0*(t+tpr)*upr+(tpr+3.0*
     . upr)*s3s+2.0*ms**2*upr)*ms**2+2.0*((t+tpr-s3s)*s+s**2-s3s*upr-
     . ms**2*upr)*mg2**2)*mg2**2+(s3s-2.0*upr)*(t+tpr)*ms**2*s-(ms**2
     . -tpr)*mg1*mg2*s**2+2.0*(ms**2-s3s)*mg2**2*s*tpr+((s+t)*mg2**2+
     . mg1**2*tpr)*(2.0*ms**2-s3s)*s-(ms**2*tpr-t*upr)*mg1*mg2*s-(2.0
     . *ms**2*tpr-s3s*upr)*s*tpr-(2.0*ms**2*t-s3s*upr)*s*t+ans17
      ans15=ans16*ct6_1*s5x
      ans2=ans3+ans15
      ans18=gs**2
      ans1=2.0*ans2*ans18
      qbg_su=ans1/(3.0*s*s3sc*s5*s5x*u6den*upr)

      qbg_tu=(-4.0*c4_tu*(((t*upr-tpr**2+s4s*tpr-s3s*tpr)*s3s-(((tpr+
     . upr)*t+tpr**2)*s4s-tpr*upr**2)-((s3s+tpr)*s3s+(s4s-upr)*tpr)*s
     . -(s+tpr)*(s3s-s4s)*ms**2-((s4s+2.0*upr-s3s)*tpr-(s3s-s4s)*s)*
     . mg2**2+(s4s-2.0*upr-s3s)*mg1**2*tpr-(s3s-s4s)*mg2**2*upr-(s*
     . s4s-t*tpr)*s3s)*s+((tpr+upr-s)*mg2**2*upr-(2.0*s3s+upr)*tpr**2
     . -mg1**2*tpr*upr-((2.0*tpr-upr)*s3s-ms**2*upr)*s-((tpr+upr)*t+
     . tpr*upr+s3s*tpr-ms**2*tpr)*upr-2.0*u6den*upr**2)*ts-((s3s*upr+
     . 2.0*s4s*tpr)*s+(2.0*s4s+upr)*tpr**2+ms**2*s*upr+mg2**2*tpr*upr
     . -((s-upr)*mg2**2+mg1**2*tpr)*upr-((tpr+upr)*t-tpr*upr-2.0*s4s*
     . tpr-s3s*tpr-ms**2*tpr)*upr)*u6den)*gs**2*mg1*mg2)/(3.0*s*s3sc*
     . s4sc*ts*u6den*upr)
*
      qbgreal=real(qbg_s1+qbg_t1+qbg_u1+qbg_st+qbg_su+qbg_tu)
ctp 
      LL_GBH = ( qbgreal
     &     -qbgd132_ll(massin,mkraemer))* (abs(mg1)+abs(mg2))**2/4.D0
*
      return
      end
************************************************************************
      real*8 function qbgd132_ll(massin,mkraemer)
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
C --- the form output for D132
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
      v132=((2.0*(tpr+upr+s)*(tpr+upr)+s**2)*gs**2)/s**2
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
      qbgd132_ll=d132
*
      return
      end
************************************************************************
      real*8 function LL_GBHSUB(massin,mkraemer)
      implicit real*8 (a-h,j-z)
      dimension massin(30),mkraemer(99)
ctp      common/susymasses/mg1,mg2,mg,ms
ctp      common/squarkw/lambda
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
      lambda = massin(30)
*
      s  = massin(1)
      t2 = massin(2)
      u2 = massin(3)
      t1 = massin(4)
      u1 = massin(5)
*
      t = t2 + mg2**2
      u = u2 + mg2**2
      ts  = t - ms**2 
      us  = u - ms**2 
      tpr = - s - t - u1 + mg2**2
      upr = - s - u - t1 + mg2**2
      s3 = - s - t - u - tpr - upr + mg1**2 + mg2**2        
      s4 = s + t + u - mg1**2 - mg2**2
      s3s = s3 + mg2**2 - ms**2
      s4s = s4 + mg1**2 - ms**2
*
C --- introduce complex s4s, s3s to perform PV integration 
C --- s4s --> s4s + i*ms*Gamma
C --- s3s --> s3s + i*ms*Gamma
      Gamma = lambda*ms
      s4s2 = s4s**2 + ms**2*Gamma**2
      s3s2 = s3s**2 + ms**2*Gamma**2
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
C --- the form output for |M(q+g->g1+g2+q)|**2
C --- spin and colour average, colour factors included
      gs=1d0
C --- factor gs2 = 4*pi*alphas NOT included
      qbg_s4=(4.0*(((s+2.0*ts)*ms**2+(t+ts)*s-2.0*(mg2**2-t)*ts)*mg2
     . **2-((t**2-t*ts+ts**2)*s+t*ts**2+((2.0*t-ts)*ts+s*t)*ms**2))*(
     . mg1+ms)*(mg1-ms)*c4_t*gs**2)/(3.0*s*s4s2*ts**2)
      qbg_s3=(4.0*((t+tpr+2.0*s+3.0*ms**2-mg2**2-mg1**2)*mg1**4-(((t+
     . tpr+s)**2+2.0*(t+tpr+s)*ms**2+mg2**4+3.0*ms**4-2.0*(t+tpr+s+ms
     . **2)*mg2**2)*mg1**2-(t+tpr+ms**2-mg2**2)*((t+tpr+s)**2-2.0*(t+
     . tpr+s)*mg2**2+mg2**4+ms**4)))*(mg2+ms)*(mg2-ms)*c4_u*gs**2)/(
     . 3.0*(t+tpr+s+ms**2-mg2**2-mg1**2)**2*s*s3s2)
*
      qbgres=qbg_s4+qbg_s3
ctp 
      LL_GBHSUB = qbgres * (abs(mg1)+abs(mg2))**2/4.D0
*
      return
      end
************************************************************************


