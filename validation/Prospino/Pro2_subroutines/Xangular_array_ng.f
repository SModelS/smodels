ctp 7/27/02 reverse the arguments in arg
ctp debug mode mostly to check completeness of integrals 
ctp by initalization to large values 
c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_NG_R(massin,NGANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      REAL*8 NGANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
     
      n        = 0
      theta_s3 = 0.D0
      call MAND_TO_ANG_NG(n,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then 
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NGANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if

c              normalization 
      NGANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]        
      NGANG(1,0,-1, 0) = 0.D0
c               tp[1], up[2]   -> (A)
      NGANG(1,2,-1,-1) = ANG2(1,2,-1,-1,1)
c               tp[1] ,s3[3]   -> (A)
      NGANG(1,3,-1,-1) = ANG2(1,3,-1,-1,1)
c               tp[1] ,u7s[9]  -> (A)     
      NGANG(1,9,-1,-1) = ANG2(1,9,-1,-1,1)
      NGANG(1,9,-1,-2) = ANG2(1,9,-1,-2,1)
c               up[2]
      NGANG(2,0,-1, 0) = 0.D0
c               up[2], s3[3]   -> (B)   
      NGANG(2,3,-1,-1) = ANG2(2,3,-1,-1,2)
c               up[2], u6s[8]  -> (B)
      NGANG(2,8,-1,-1) = ANG2(2,8,-1,-1,2)
      NGANG(2,8,-1,-2) = ANG2(2,8,-1,-2,2)
c               s3[3]          -> (C)
      NGANG(3,0,-1,0)  = ANG2(3,0,-1, 0,3)
      NGANG(3,0,-2,0)  = ANG2(3,0,-2, 0,3)
c               s3[3], u6s[8]  -> (C)
      NGANG(3,8,-1,-1) = ANG2(3,8,-1,-1,3)
      NGANG(3,8,-1,-2) = ANG2(3,8,-1,-2,3)
      NGANG(3,8,-2,-1) = ANG2(3,8,-2,-1,3)
      NGANG(3,8,-2,-2) = ANG2(3,8,-2,-2,3)
c               s3[3], u7s[9]  -> (C)
      NGANG(3,9,-1,-1) = ANG2(3,9,-1,-1,3)
      NGANG(3,9,-1,-2) = ANG2(3,9,-1,-2,3)
      NGANG(3,9,-2,-1) = ANG2(3,9,-2,-1,3)
      NGANG(3,9,-2,-2) = ANG2(3,9,-2,-2,3)
c               u6s[8]         -> (A)
      NGANG(8,0,-1,0)  = ANG2(8,0,-1, 0,1)
      NGANG(8,0,-2,0)  = ANG2(8,0,-2, 0,1)
c               u6s[8], u7s[9] -> (A)
      NGANG(8,9,-1,-1) = ANG2(8,9,-1,-1,1)
c               u7s[9]         -> (B)
      NGANG(9,0,-1,0)  = ANG2(9,0,-1, 0,2)
      NGANG(9,0,-2,0)  = ANG2(9,0,-2, 0,2)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_NG_C(massin,NGANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5
      REAL*8 pi
      REAL*8 gams,theta_s3
      REAL*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      REAL*8 CBP1P1_NEW
      REAL*8 NGANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))

      pi = 4.D0 * atan(1.D0)
      gams = massin(25)

      n        = 0
      theta_s3 = 0.D0
      call MAND_TO_ANG_NG(n,massin,arg,theta_s3,arg_x,arg_y)

      if (ldebug) then
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NGANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c            needed for NNQGH
c              normalization 
      NGANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]        
      NGANG(1,0,-1, 0) = 0.D0
c               tp[1] ,s3s[4]  -> (A)
      NGANG(1,4,-1,-1) = pi*CBP1P1_NEW(gams,arg(1,1,1),arg(1,4,1))
      NGANG(1,4,1,-1)  = ANG2(1,4, 1,-1,1)
c               tp[1], u7[7]   -> (A)
      NGANG(1,7,-1, 1) = ANG2(1,7,-1, 1,1)
      NGANG(1,7, 1,-1) = ANG2(1,7, 1,-1,1)
c               tp[1] ,u7s[9]  -> (A)     
      NGANG(1,9,-1,-1) = ANG2(1,9,-1,-1,1)
      NGANG(1,9,-1,-2) = ANG2(1,9,-1,-2,1)
      NGANG(1,9,1,-1)  = ANG2(1,9, 1,-1,1)
c               s3[3]          -> (C)
      NGANG(3,0,-1,0)  = ANG2(3,0,-1, 0,3)
      NGANG(3,0, 1,0)  = ANG2(3,0, 1, 0,3)
c               s3[3], u7[7]   -> u7[7], s3[3] -> (B)
      NGANG(3,7,-1, 1) = ANG2(7,3, 1,-1,2)
c               s3[3], u7s[9]  -> (C)
      NGANG(3,9,-1,-1) = ANG2(3,9,-1,-1,3)
      NGANG(3,9,-1,-2) = ANG2(3,9,-1,-2,3)
c               s3s[4]         -> (C)
      NGANG(4,0,-1, 0) = ANG2(4,0,-1, 0,3)
      NGANG(4,0,-2, 0) = ANG2(4,0,-2, 0,3)
c               s3s[4], u6[6]  -> (C)
      NGANG(4,6,-2, 1) = ANG2(4,6,-2, 1,3)
c               s3s[4], u7[7]  -> u7[7], s3s[4] -> (B)
      NGANG(4,7,-1, 1) = ANG2(7,4, 1,-1,2)
c               s3s[4], u7s[9] -> (C)
      NGANG(4,9,-1,-1) = ANG2(4,9,-1,-1,3)
      NGANG(4,9,-1,-2) = ANG2(4,9,-1,-2,3)
      NGANG(4,9,-2,-1) = ANG2(4,9,-2,-1,3)
      NGANG(4,9,-2,-2) = ANG2(4,9,-2,-2,3)
c               u6[6] ,u7s[9]  -> (A)     
      NGANG(6,9, 1,-1) = ANG2(6,9, 1,-1,1)
      NGANG(6,9, 1,-2) = ANG2(6,9, 1,-2,1)
c               u7s[9]         -> (B)
      NGANG(9,0,-1,0)  = ANG2(9,0,-1, 0,2)
      NGANG(9,0,-2,0)  = ANG2(9,0,-2, 0,2)

c            in addition needed for NNGBH
c               up[2]
      NGANG(2,0,-1, 0) = 0.D0
c               up[2] ,s3s[4]  -> (B)
      NGANG(2,4,-1,-1) = pi*CBP1P1_NEW(gams,arg(1,2,2),arg(1,4,2))
c               up[2], u6[6]   -> (B)
      NGANG(2,6,-1, 1) = ANG2(2,6,-1, 1,2)
c               tp[2] ,u6s[8]  -> (B)     
      NGANG(2,8,-1,-1) = ANG2(2,8,-1,-1,2)
      NGANG(2,8,-1,-2) = ANG2(2,8,-1,-2,2)
c               s3[3], u6s[8]  -> (C)
      NGANG(3,8,-1,-1) = ANG2(3,8,-1,-1,3)
      NGANG(3,8,-1,-2) = ANG2(3,8,-1,-2,3)
c               s3s[4], u6s[8] -> (C)
      NGANG(4,8,-1,-1) = ANG2(4,8,-1,-1,3)
      NGANG(4,8,-1,-2) = ANG2(4,8,-1,-2,3)
c               u7[7] ,u6s[8]  -> (B)     
      NGANG(7,8, 1,-1) = ANG2(7,8, 1,-1,2)
      NGANG(7,8, 1,-2) = ANG2(7,8, 1,-2,2)
c               u6s[8]         -> (A)
      NGANG(8,0,-1,0)  = ANG2(8,0,-1, 0,1)
      NGANG(8,0,-2,0)  = ANG2(8,0,-2, 0,1)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_NG_S
     &           (n,theta_s3,theta_s4,massin,NGANG,IMANG)

      IMPLICIT NONE

      integer n1,n2,n3,n4,n5,n
      REAL*8 s3,s3s,theta_s3,theta_s4,yesno_s3,yesno_s4,m2,ms,gams
      REAL*8 massin(1:30)
      REAL*8 arg(1:3,0:9,1:3), arg_x(0:9), arg_y(0:9)
      REAL*8 NGANG(0:9,0:9,-2:2,-2:2)
      real*8 IMANG(0:9,0:9,-2:2,-2:2)
      real*8 ANG2_EXT,ANG2,ANG1_EXT,ANG1,AIM2_EXT,AIM2
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
      AIM2(n4,n5,n1,n2,n3)=AIM2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
      ANG1(n1,n2)         =ANG1_EXT(n2,arg_x(n1),arg_y(n1))

      s3   = massin(4)
      m2   = massin(7)
      ms   = massin(11)
      gams = massin(25)
      yesno_s4 = massin(26)
      yesno_s3 = massin(27)

      call MAND_TO_ANG_NG(n,massin,arg,theta_s3,arg_x,arg_y)

      if (ldebug) then 
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NGANG(n1,n2,n3,n4) = 9.99999D+99
                     IMANG(n1,n2,n3,n4) = 0.D0
                  end do
               end do
            end do
         end do
      end if

      if (n.eq.1) then 
c            needed for NNQGOS1 and NNGBOS1
         if (yesno_s3.eq.1.D0) then 
c            variable needed here 
            s3s = s3 + m2**2 - ms**2
            s3s = s3s**2 + ms**2*gams**2
            if (s3s.ge.0.D0) then 
               s3s  = sqrt( s3s )
            else 
               print*, " ANGULAR_ARRAY_NG_S: s3s sqrt fixed ",s3s
               s3s = 0.D0
            end if 
ctp            s3s = sign(1.D0,s3s) * sqrt( s3s**2 + ms**2*gams**2 )
c               s3s[4]         -> (C)
            NGANG(4,0,-2, 0) = ANG1(0, 0)/s3s**2
c               s3s[4], u6[6]  -> (C)
            NGANG(4,6,-2, 1) = ANG1(6, 1)/s3s**2
c               s3s[4], u7[7]  -> (C)
            NGANG(4,7,-2, 1) = ANG1(7, 1)/s3s**2
c               s3s[4], u6s[8] -> (C)
            NGANG(4,8,-2,-2) = ANG1(8,-2)/s3s**2
            NGANG(4,8,-2,-1) = ANG1(8,-1)/s3s**2
c               s3s[4], u7s[9] -> (C)
            NGANG(4,9,-2,-2) = ANG1(9,-2)/s3s**2
            NGANG(4,9,-2,-1) = ANG1(9,-1)/s3s**2
c              normalization 
         end if 
         if (yesno_s4.eq.1.D0) then 
            NGANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
         end if 
      else if (n.eq.2) then 
c            needed for NNQGOS2 and NNGBOS2
c               s3s[4]         -> (C)
         NGANG(4,0,-1, 0) = ANG2(4,0,-1, 0,3)
c               s3s[4], u7[7]  -> u7[7], s3s[4] -> (B)
         NGANG(4,7,-1, 1) = ANG2(7,4, 1,-1,2)
c               s3s[4], u6s[8] -> (C)
         NGANG(4,8,-1,-1) = ANG2(4,8,-1,-1,3)
c               s3s[4], u7s[9] -> (C)
         NGANG(4,9,-1,-1) = ANG2(4,9,-1,-1,3)
         if ((theta_s3.eq.1.D0).and.(theta_s4.eq.1.D0)) then 
c               s3s[4]         -> (C)
            IMANG(4,0,-1, 0) = AIM2(4,0,-1, 0,3)
c               s3s[4], u7[7]  -> (C)
            IMANG(4,7,-1, 1) = AIM2(4,7,-1, 1,3)
c               s3s[4], u7s[9] -> (C)
            IMANG(4,9,-1,-1) = AIM2(4,9,-1,-1,3)
c               s3s[4], u6s[8] -> (C)
            IMANG(4,8,-1,-1) = AIM2(4,8,-1,-1,3)
         else 
            IMANG(4,0,-1, 0) = 0.D0
            IMANG(4,7,-1, 1) = 0.D0
            IMANG(4,9,-1,-1) = 0.D0
            IMANG(4,8,-1,-1) = 0.D0
         end if 
      end if

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all the soft variables in the denominator :                      c
c         first argument 1,2,3 the reference frame A,B,C               c
c         second argument 1,...,9 the mandelstam variables:            c
c                                                                      c
c              1    tp  = (k2-k3)^2                                    c
c              2    up  = (k1-k3)^2                                    c
c              3    s3  = (k3-p2)^2 - m2^2                             c
c              4   s3s  = (k3-p2)^2 - ms^2                             c
c              5    s5  = (p1+p2)^2        is cm connected to s3[3]    c
c              6    u6  = (k2-p1)^2 - m1^2 is cm connected to tp[1]    c
c              7    u7  = (k1-p1)^2 - m1^2 is cm connected to up[2]    c
c              8    u6s = (k2-p1)^2 - ms^2                             c
c              9    u7s = (k1-p1)^2 - ms^2                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MAND_TO_ANG_NG(n,massin,arg,theta_s3,arg_x,arg_y)

      integer n,n1,n2,n3
      real*8  massin(1:30),arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
     &       ,theta_s3,s,t2,s4,m1,m2,ms,s3,s3s
     &       ,gams,yesno_s3
     &       ,m12,m22,u1,u2,s42
     &       ,norm,w1,w2,w3,e1,e2,p,cpa,cpb,spa,spb,c1
     &       ,dum1,dum2 

      logical ldebug 
      parameter( ldebug=.false. )
           
c               initialize the output array
      do n1=1,3
         do n2=0,9
            do n3=1,3
               arg(n1,n2,n3) = 0.D0
            end do
         end do
      end do

      do n3=0,9
         arg_x(n3)     = 0.D0
         arg_y(n3)     = 0.D0
      end do

c               define the different variables
      s    = massin(1)
      t2   = massin(2) 
      s4   = massin(3)
      m1   = massin(6) 
      m2   = massin(7)
      ms   = massin(11)
      gams = massin(25)
      yesno_s3 = massin(27)

c              the two outgoing masses set to m1,m2 [see my thesis]
      m12 = m1**2
      m22 = m2**2
      u1  = s4 - s - t2
      u2  = u1 + m1**2 - m2**2 
      s42 = s4 + m1**2 - m2**2  

c              define in the three reference frames (A,B,C):
      norm = s4+m12
      if (norm.ge.0.D0) then 
         norm = 2.D0*sqrt(norm)
      else 
         print*, " MAND_TO_ANG_NG: serious problem with s4 ",s4
         call HARD_STOP
      end if
      w1  = (s + u2)                        /norm
      w2  = (s + t2)                        /norm
      w3  =  s4                             /norm
      e1  = (s4 + 2.D0*m12)                 /norm
      e2  = -(t2 + u2 + 2.D0*m22)           /norm
      p   = (t2 + u2)**2 - 4.D0*m22*s
      if (p.ge.0.D0) then 
         p = sqrt(p)/norm
      else 
         print*, " MAND_TO_ANG_NG: serious problem with p ",p
         call HARD_STOP
      end if 

c              reference frames (A,B)
      cpa =(t2*s42 -s*(u2+2.D0*m22))/(s+t2)/sqrt((t2+u2)**2-4.D0*m22*s)
      cpb =(u2*s42 -s*(t2+2.D0*m22))/(s+u2)/sqrt((u2+t2)**2-4.D0*m22*s)

      if ((abs(cpa).gt.1.D0).and.(abs(cpa).lt.1.1D0)) then 
         cpa = sign(1.D0,cpa)
      end if 

      if ((abs(cpb).gt.1.D0).and.(abs(cpb).lt.1.1D0)) then 
         cpb = sign(1.D0,cpb)
      end if

      spa = sqrt(1.D0 - cpa**2)
      spb = sqrt(1.D0 - cpb**2)

c              reference frame (A): k2||z, cm(k3,p1)
c              two terms   :    tp,u6,u6s
c              a=-b        :    tp
c              A^2=B^2+C^2 :    up

      arg(1,1,1) = -2.D0*w2*w3
      arg(2,1,1) = -arg(1,1,1)
      arg(3,1,1) =  0.D0

      arg(1,2,1) = -2.D0*w1*w3
      arg(2,2,1) = +2.D0*w3*p*cpa-2.D0*w2*w3
      arg(3,2,1) = +2.D0*w3*p*spa

      arg(1,3,1) = +2.D0*w3*e2
      arg(2,3,1) = -2.D0*w3*p*cpa
      arg(3,3,1) = -2.D0*w3*p*spa 

c              reference frame (B): k1||z, cm(k3,p1)
c              two terms   :    up,u7,u7s
c              a=-b        :    up
c              A^2=B^2+C^2 :    tp

      arg(1,1,2) = -2.D0*w2*w3
      arg(2,1,2) = +2.D0*w3*p*cpb-2.D0*w1*w3
      arg(3,1,2) = +2.D0*w3*p*spb

      arg(1,2,2) = -2.D0*w1*w3
      arg(2,2,2) = -arg(1,2,2)
      arg(3,2,2) =  0.D0

      arg(1,3,2) = +2.D0*w3*e2
      arg(2,3,2) = -2.D0*w3*p*cpb
      arg(3,3,2) = -2.D0*w3*p*spb

c              reference frame (C): p2||z, cm(k3,p1)
c              two terms   :    s3,s5,s3s
c              a=-b        :    -
c              A^2=B^2+C^2 :    tp,up

      arg(1,1,3) = -2.D0*w2*w3
      arg(2,1,3) = +2.D0*w2*w3*cpa
      arg(3,1,3) = +2.D0*w2*w3*spa

      arg(1,2,3) = -2.D0*w1*w3
      arg(2,2,3) = -2.D0*w2*w3*cpa+2.D0*w3*p
      arg(3,2,3) = -2.D0*w2*w3*spa

      arg(1,3,3) = +2.D0*w3*e2
      arg(2,3,3) = -2.D0*w3*p
      arg(3,3,3) =  0.D0

      do n1=1,3 
c            cm(k3,p1) system s3 -> s5 
         arg(1,5,n1) = +2.D0*e1*e2 + m12 + m22
         arg(2,5,n1) = -arg(2,3,n1)
         arg(3,5,n1) = -arg(3,3,n1)
c            cm(k3,p1) system tp -> u6 
         arg(1,6,n1) = -2.D0*w2*e1
         arg(2,6,n1) = -arg(2,1,n1)
         arg(3,6,n1) = -arg(3,1,n1)
c            cm(k3,p1) system up -> u7 
         arg(1,7,n1) = -2.D0*w1*e1
         arg(2,7,n1) = -arg(2,2,n1)
         arg(3,7,n1) = -arg(3,2,n1)
c            derive u6 -> u6s
         arg(1,8,n1) = arg(1,6,n1) + m1**2 - ms**2 
         arg(2,8,n1) = arg(2,6,n1)
         arg(3,8,n1) = arg(3,6,n1)
c            derive u7 -> u7s
         arg(1,9,n1) = arg(1,7,n1) + m1**2 - ms**2 
         arg(2,9,n1) = arg(2,7,n1)
         arg(3,9,n1) = arg(3,7,n1)
c            derive s3 -> s3s including regularization 
         if ((n.eq.1).and.(theta_s3.eq.1.D0)) then
            arg(1,4,n1) = arg(1,3,n1) + m2**2 - ms**2 - gams**2 
         else 
            arg(1,4,n1) = arg(1,3,n1) + m2**2 - ms**2
         end if
         arg(2,4,n1) = arg(2,3,n1)
         arg(3,4,n1) = arg(3,3,n1)
      end do

c              special case of reference frame (C)
c              only for theta_s3=1 and asked for s3 subtraction 
c              only for the intgerals including u6s,u7s
      if ((yesno_s3.eq.1.D0).and.(n.eq.1)) then 
         s3  = massin(4)
         s3s = s3 + m2**2 - ms**2
         s3s = s3s**2 + ms**2*gams**2 
         if (s3s.ge.0.D0) then 
            s3s  = sqrt( s3s )
         else 
            print*, " MAD_TO_ANG: s3s sqrt fixed ",s3s
            s3s = 0.D0
         end if 

         c1 = ( 2.D0*w3*e2 - s3 )/(2.D0*p*w3)
         if (abs(c1).gt.1.D0) then 
            print*, " MAND_TO_ANG_NG: c1 fixed ",c1
            c1 = 0.D0
         end if 

         do n1=0,9
            arg_x(n1) = arg(1,n1,3) + arg(2,n1,3)*c1
            arg_y(n1) = arg(3,n1,3)*sqrt(1.D0-c1**2)
         end do
      end if


c         check of the structure in debug mode 
      if (ldebug) then 
         do n1=1,3,1
            do n2=1,9,1
               if (abs(arg(3,n2,n1)).le.1.D-16) then
                  print*, " only a,b   : ",n1,n2,arg(3,n2,n1)
               end if
            end do
         end do

c         special case a=-b 
         do n1=1,3,1
            do n2=1,9,1
               dum1 = arg(1,n2,n1) + arg(2,n2,n1)
               if ((abs(arg(3,n2,n1)).le.1.D-16).and.
     &             (abs(dum1).le.1.D-16))         then
                  print*, " a = -b     : ",n1,n2,dum1
               end if
            end do
         end do

c         special case A^2=B^2+C^2 
c           n.b. the dum2 cutoff has to be bigger than 1.e-12
         do n1=1,3,1
            do n2=1,9,1
               dum1 = arg(1,n2,n1) + arg(2,n2,n1)
               dum2 = (arg(1,n2,n1)**2-arg(2,n2,n1)**2-arg(3,n2,n1)**2)
               dum2 = dum2/arg(1,n2,n1)**2
               if ((abs(dum1).gt.1.D-16).and.
     &             (abs(dum2).le.1.D-12))         then
                  print*, "A^2=B^2+C^2 : ",n1,n2,dum2
               end if
            end do
         end do
      end if 

      return 
      end 





































