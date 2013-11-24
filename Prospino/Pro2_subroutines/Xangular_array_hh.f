c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HH_QB(massin,HHANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5,ndim
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
      REAL*8 HHANG(0:12,0:12,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
     
      n        = 0
      ndim     = 2 
      theta_s3 = 0.D0
      call MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then 
         do n1=0,12,1
            do n2=0,12,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     HHANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if

c              normalization 
      HHANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1], up[2]  -> (A) or (B)
      HHANG(1,2,-1,-1) = ANG2(1,2,-1,-1,1)
      HHANG(1,2, 1,-1) = ANG2(2,1,-1, 1,2)
c               tp[1] ,t1t[4]  -> (A)
      HHANG(1,4,-1,-2) = ANG2(1,4,-1,-2,1)
      HHANG(1,4,-1,-1) = ANG2(1,4,-1,-1,1)
c               tp[1], s5[5]  -> (A)
      HHANG(1, 5,-1,-2) = ANG2(1, 5,-1,-2,1)
      HHANG(1, 5,-1,-1) = ANG2(1, 5,-1,-1,1)
      HHANG(1,10,-1,-2) = ANG2(1,10,-1,-2,1)
      HHANG(1,10,-1,-1) = ANG2(1,10,-1,-1,1)
      HHANG(1,11,-1,-2) = ANG2(1,11,-1,-2,1)
      HHANG(1,11,-1,-1) = ANG2(1,11,-1,-1,1)
      HHANG(1,12,-1,-2) = ANG2(1,12,-1,-2,1)
      HHANG(1,12,-1,-1) = ANG2(1,12,-1,-1,1)
c               up[2]
      HHANG(2,0,-1, 0) = 0.D0
c               up[1], s5[5]  -> (B)
      HHANG(2, 5,-1,-2) = ANG2(2, 5,-1,-2,2)
      HHANG(2, 5,-1,-1) = ANG2(2, 5,-1,-1,2)
      HHANG(2,10,-1,-2) = ANG2(2,10,-1,-2,2)
      HHANG(2,10,-1,-1) = ANG2(2,10,-1,-1,2)
      HHANG(2,11,-1,-2) = ANG2(2,11,-1,-2,2)
      HHANG(2,11,-1,-1) = ANG2(2,11,-1,-1,2)
      HHANG(2,12,-1,-2) = ANG2(2,12,-1,-2,2)
      HHANG(2,12,-1,-1) = ANG2(2,12,-1,-1,2)
c               t1t[4]        -> (B)
c               (check those)
      HHANG(4,0,-2, 0) = ANG2(4,0,-2, 0,2)
      HHANG(4,0,-1, 0) = ANG2(4,0,-1, 0,2)
c               t1t[4], s5[5]  -> (B)
      HHANG(4, 5,-2, 1) = ANG2(4, 5,-2, 1,2) 
      HHANG(4, 5,-1, 1) = ANG2(4, 5,-1, 1,2) 
      HHANG(4, 5,-1,-1) = ANG2(4, 5,-1,-1,2) 
      HHANG(4,10,-1,-1) = ANG2(4,10,-1,-1,2) 
      HHANG(4,11,-1,-1) = ANG2(4,11,-1,-1,2) 
      HHANG(4,12,-1,-1) = ANG2(4,12,-1,-1,2) 
c               s5[5]          -> (C)
      HHANG( 5,0,-2, 0) = ANG2( 5,0,-2, 0,3)
      HHANG( 5,0,-1, 0) = ANG2( 5,0,-1, 0,3)
      HHANG(10,0,-2, 0) = ANG2(10,0,-2, 0,3)
      HHANG(10,0,-1, 0) = ANG2(10,0,-1, 0,3)
      HHANG(11,0,-1, 0) = ANG2(11,0,-1, 0,3)
      HHANG(12,0,-1, 0) = ANG2(12,0,-1, 0,3)

      return
      end

c ----------------------------------------------------------------------
c remember one complication: replace divergent integral 
c like 1/tp/s3s with C1P1P1_NEW routine
      subroutine ANGULAR_ARRAY_HH_QG(massin,theta_s3,HHANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5,ndim
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
      REAL*8 HHANG(0:12,0:12,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
     
      n        = 0
      ndim     = 2 
      call MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)

      if (ldebug) then 
         do n1=0,12,1
            do n2=0,12,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     HHANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if

c              normalization 
      HHANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]        
      HHANG(1,0,-1, 0) = 0.D0
c               tp[1] ,t1t[4]  -> (A)
      HHANG(1,4,-1,-2) = ANG2(1,4,-1,-2,1)
      HHANG(1,4,-1,-1) = ANG2(1,4,-1,-1,1)
c               tp[1], s5[5]  -> (A)
      HHANG(1, 5,-1,-2) = ANG2(1, 5,-1,-2,1)
      HHANG(1, 5,-1,-1) = ANG2(1, 5,-1,-1,1)
c               tp[1], s3t[8]  -> (A)
      HHANG(1, 8,-1,-1) = ANG2(1, 8,-1,-1,1)
c               tp[1], s5z[10]  -> (A)
      HHANG(1,10,-1,-2) = ANG2(1,10,-1,-2,1)
      HHANG(1,10,-1,-1) = ANG2(1,10,-1,-1,1)
      HHANG(1,11,-1,-2) = ANG2(1,11,-1,-2,1)
      HHANG(1,11,-1,-1) = ANG2(1,11,-1,-1,1)
      HHANG(1,12,-1,-2) = ANG2(1,12,-1,-2,1)
      HHANG(1,12,-1,-1) = ANG2(1,12,-1,-1,1)
c               t1t[4], s5[5]  -> (B)
      HHANG(4,5,-1,-1) = ANG2(4,5,-1,-1,2) 
c               t1t[4], s3t[8]  -> (B)
      HHANG(4,8,-2,-2) = ANG2(4,8,-2,-2,2) 
      HHANG(4,8,-2,-1) = ANG2(4,8,-2,-1,2) 
      HHANG(4,8,-1,-2) = ANG2(4,8,-1,-2,2) 
      HHANG(4,8,-1,-1) = ANG2(4,8,-1,-1,2) 
      HHANG(4,8, 1,-2) = ANG2(4,8, 1,-2,2) 
      HHANG(8,4,-2,-2) = HHANG(4,8,-2,-2) 
      HHANG(8,4,-2,-1) = HHANG(4,8,-1,-2) 
      HHANG(8,4,-2, 1) = HHANG(4,8, 1,-2) 
c               t1t[4], s5z[10]  -> (B)
      HHANG(4,10,-1,-1) = ANG2(4,10,-1,-1,2) 
      HHANG(4,11,-1,-1) = ANG2(4,11,-1,-1,2) 
      HHANG(4,12,-1,-1) = ANG2(4,12,-1,-1,2) 
c               s5[5]          -> (C)
      HHANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)
      HHANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
c               s5[5], t1[7]   -> (C) for example
      HHANG(5,7,-2, 1) = ANG2(5,7,-2, 1,3)
ctp      HHANG(5,7,-2, 1) = ANG2(7,5, 1,-2,2)
      HHANG(5,7,-1, 1) = ANG2(5,7,-1, 1,3)
c               t1[7], s3t[8]  -> (C) for example 
      HHANG(7,8, 1,-2) = ANG2(7,8, 1,-2,2)
      HHANG(7,8, 1,-1) = ANG2(7,8, 1,-1,2)
      HHANG(7,0,-1, 0) = ANG2(7,0,-1, 0,2)
c               s3t[8]         -> (C)
      HHANG(8,0,-2, 0) = ANG2(8,0,-2, 0,3)
      HHANG(8,0,-1, 0) = ANG2(8,0,-1, 0,3)
c               s5z[10]          -> (C)
      HHANG(10,0,-2, 0) = ANG2(10,0,-2, 0,3)
      HHANG(10,0,-1, 0) = ANG2(10,0,-1, 0,3)
      HHANG(11,0,-2, 0) = ANG2(11,0,-2, 0,3)
      HHANG(11,0,-1, 0) = ANG2(11,0,-1, 0,3)
      HHANG(12,0,-2, 0) = ANG2(12,0,-2, 0,3)
      HHANG(12,0,-1, 0) = ANG2(12,0,-1, 0,3)
c               s5z[10], t1[7]   -> (C) for example
      HHANG(10,7,-2, 1) = ANG2(10,7,-2, 1,3)
      HHANG(10,7,-1, 1) = ANG2(10,7,-1, 1,3)
      HHANG(11,7,-2, 1) = ANG2(11,7,-2, 1,3)
      HHANG(11,7,-1, 1) = ANG2(11,7,-1, 1,3)
      HHANG(12,7,-2, 1) = ANG2(12,7,-2, 1,3)
      HHANG(12,7,-1, 1) = ANG2(12,7,-1, 1,3)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HH_QGOS(massin,theta_s3,HHANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,ndim
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
      REAL*8 HHANG(0:12,0:12,-2:2,-2:2),massin(1:30)
      real*8 m2,mt,gamt,s3,s3t
      REAL*8 ANG1_EXT,ANG1

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG1(n1,n2)         =ANG1_EXT(n2,arg_x(n1),arg_y(n1))
     
c            n should not matter
      n        = 0
      ndim     = 1 
      call MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      s3   = massin(4)
      m2   = massin(6)
      mt   = massin(7)
      gamt = massin(25)
      
      if (ldebug) then 
         do n1=0,12,1
            do n2=0,12,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     HHANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c            the correct denominator for a numerical s3 integration   
      s3t = s3 + m2**2 - mt**2 
      s3t = sign(1.D0,s3t) * sqrt( s3t**2 + mt**2*gamt**2 )
      
c               n.b. no fixed cms for ANG1
c               s3t[8], t1t[4]
      HHANG(8,4,-2,-2) = ANG1(4,-2)/s3t**2
      HHANG(8,4,-2,-1) = ANG1(4,-1)/s3t**2
      HHANG(8,4,-2, 1) = ANG1(4, 1)/s3t**2
c               s3t[8]          -> (C)
      HHANG(8,0,-1, 0) = ANG1(0, 0)/s3t
      HHANG(8,0,-2, 0) = ANG1(0, 0)/s3t**2

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HH_GB(massin,theta_s3,HHANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5,ndim
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
      REAL*8 HHANG(0:12,0:12,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
     
      n        = 0
      ndim     = 2 
      call MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then 
         do n1=0,12,1
            do n2=0,12,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     HHANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if

c              normalization 
      HHANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               up[2]
      HHANG(2,0,-1, 0) = 0.D0
c               up[2], s5[5]  -> (B)
      HHANG(2, 5,-1,-2) = ANG2(2, 5,-1,-2,2)
      HHANG(2, 5,-1,-1) = ANG2(2, 5,-1,-1,2)
      HHANG(2, 5,-1, 1) = ANG2(2, 5,-1, 1,2)
      HHANG(2, 5,-1, 2) = ANG2(2, 5,-1, 2,2)
c               up[2], s5z[10]  -> (B)
      HHANG(2,10,-1,-2) = ANG2(2,10,-1,-2,2)
      HHANG(2,10,-1,-1) = ANG2(2,10,-1,-1,2)
      HHANG(2,11,-1,-2) = ANG2(2,11,-1,-2,2)
      HHANG(2,11,-1,-1) = ANG2(2,11,-1,-1,2)
      HHANG(2,12,-1,-2) = ANG2(2,12,-1,-2,2)
      HHANG(2,12,-1,-1) = ANG2(2,12,-1,-1,2)
c               s5[5]          -> (C)
      HHANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)
      HHANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
      HHANG(5,0, 1, 0) = ANG2(5,0, 1, 0,3)
c               s5[5], t1[7]   -> (C) for example
      HHANG(5,7,-2, 1) = ANG2(5,7,-2, 1,3)
      HHANG(5,7,-1, 1) = ANG2(5,7,-1, 1,3)
c               t1t[7]          -> (B)
      HHANG(7,0, 1, 0) = ANG2(7,0, 1, 0,2)
c               s5z[10]         -> (C)
      HHANG(10,0,-2, 0) = ANG2(10,0,-2, 0,3)
      HHANG(10,0,-1, 0) = ANG2(10,0,-1, 0,3)
      HHANG(11,0,-2, 0) = ANG2(11,0,-2, 0,3)
      HHANG(11,0,-1, 0) = ANG2(11,0,-1, 0,3)
      HHANG(12,0,-2, 0) = ANG2(12,0,-2, 0,3)
      HHANG(12,0,-1, 0) = ANG2(12,0,-1, 0,3)
c               s5z[10], t1[7]   -> (C) for example
      HHANG(10,7,-2, 1) = ANG2(10,7,-2, 1,3)
      HHANG(10,7,-1, 1) = ANG2(10,7,-1, 1,3)
      HHANG(11,7,-2, 1) = ANG2(11,7,-2, 1,3)
      HHANG(11,7,-1, 1) = ANG2(11,7,-1, 1,3)
      HHANG(12,7,-2, 1) = ANG2(12,7,-2, 1,3)
      HHANG(12,7,-1, 1) = ANG2(12,7,-1, 1,3)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HH_GBOS(massin,theta_s3,HHANG)

      IMPLICIT NONE

      integer n,n1,n2,n3,n4,n5,ndim
      REAL*8 theta_s3
      REAL*8 arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
      REAL*8 HHANG(0:12,0:12,-2:2,-2:2),massin(1:30)
      REAL*8 ANG2_EXT,ANG2

      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
     
      n        = 0
      ndim     = 2 
      call MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then 
         do n1=0,12,1
            do n2=0,12,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     HHANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if

c              normalization 
      HHANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               s5[5]          -> (C)
      HHANG(5,0, 1, 0) = ANG2(5,0, 1, 0,3)
c               t1t[7]          -> (B)
      HHANG(7,0, 1, 0) = ANG2(7,0, 1, 0,2)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all the soft variables in the denominator :                      c
c         first argument 1,2,3 the reference frame A,B,C               c
c         second argument 1,...,12 the mandelstam variables:           c
c                                                                      c
c              1    tp  = (k2-k3)^2                                    c
c              2    up  = (k1-k3)^2                                    c
c              3    s3  = (k3-p2)^2 - m2^2                             c
c              4    t1t = t1+m1^2-mt^2                                 c
c              5    s5  = (p1+p2)^2        is cm connected to s3[3]    c
c              6    u1  = (k2-p1)^2 - m1^2 is cm connected to tp[1]    c
c              7    t1  = (k1-p1)^2 - m1^2 is cm connected to up[2]    c
c              8    s3t = s3+m2^2-mt^2                                 c
c              9    u1t = u1+m1^2-mt^2                                 c
c             10    s5z = s5-mz^2                                      c
c             11    s51 = s5-mh1^2                                     c
c             12    s52 = s5-mh2^2                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MAND_TO_ANG_HH(n,ndim,massin,arg,theta_s3,arg_x,arg_y)

      integer n,n1,n2,n3,ndim
      real*8  massin(1:30),arg(1:3,0:12,1:3),arg_x(0:12),arg_y(0:12)
     &       ,theta_s3,s,t2,s4,m1,m2,mt,mz,mh1,mh2,s3,s3t
     &       ,gamt
     &       ,m12,m22,u1,u2,s42
     &       ,norm,w1,w2,w3,e1,e2,p,cpa,cpb,spa,spb,c1
     &       ,dum1,dum2 

      logical ldebug 
      parameter( ldebug=.false. )
           
c               initialize the output array
      do n1=1,3
         do n2=0,12
            do n3=1,3
               arg(n1,n2,n3) = 0.D0
            end do
         end do
      end do

      do n3=0,12
         arg_x(n3)     = 0.D0
         arg_y(n3)     = 0.D0
      end do

c               define the different variables
      s     = massin(1)
      t2    = massin(2) 
      s4    = massin(3)
      m1    = massin(6) 
      m2    = massin(6)
      mt    = massin(7)
      mz    = massin(8)
      mh1   = massin(9)
      mh2   = massin(10)
      gamt  = massin(25)

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
         print*, " MAND_TO_ANG_HH: serious problem with s4 ",s4
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
         print*, " MAND_TO_ANG_HH: serious problem with p ",p
         call HARD_STOP
      end if 

c              reference frames (A,B)
      cpa =(t2*s42 -s*(u2+2.D0*m22))/(s+t2)/sqrt((t2+u2)**2-4.D0*m22*s)
      cpb =(u2*s42 -s*(t2+2.D0*m22))/(s+u2)/sqrt((u2+t2)**2-4.D0*m22*s)

      if ((abs(cpa).gt.1.D0).and.(abs(cpa).lt.1.1D0)) then 
         if ( cpa .gt. 0.D0 ) then
            cpa =  1.D0
         else 
            cpa = -1.D0
         end if
      end if 

      if ((abs(cpb).gt.1.D0).and.(abs(cpb).lt.1.1D0)) then 
         if ( cpb .gt. 0.D0 ) then
            cpb =  1.D0
         else 
            cpb = -1.D0
         end if
      end if

      spa = sqrt(1.D0 - cpa**2)
      spb = sqrt(1.D0 - cpb**2)

c              reference frame (A): k2||z, cm(k3,p1)
c              two terms   :    tp,u1,u1t
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
c              two terms   :    up,t1,t1t
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
c              two terms   :    s3,s5,s3t
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
c            cm(k3,p1) system tp -> u1 
         arg(1,6,n1) = -2.D0*w2*e1
         arg(2,6,n1) = -arg(2,1,n1)
         arg(3,6,n1) = -arg(3,1,n1)
c            cm(k3,p1) system up -> t1 
         arg(1,7,n1) = -2.D0*w1*e1
         arg(2,7,n1) = -arg(2,2,n1)
         arg(3,7,n1) = -arg(3,2,n1)
c            derive t1 -> t1t
         arg(1,4,n1) = arg(1,7,n1) + m1**2 - mt**2 
         arg(2,4,n1) = arg(2,7,n1)
         arg(3,4,n1) = arg(3,7,n1)
c            derive s3 -> s3t
         arg(1,8,n1) = arg(1,3,n1) + m2**2 - mt**2 
         arg(2,8,n1) = arg(2,3,n1)
         arg(3,8,n1) = arg(3,3,n1)
c            derive u1 -> u1t
         arg(1,9,n1) = arg(1,6,n1) + m1**2 - mt**2 
         arg(2,9,n1) = arg(2,6,n1)
         arg(3,9,n1) = arg(3,6,n1)
c            derive s5 -> s5z
         arg(1,10,n1) = arg(1,5,n1) - mz**2 
         arg(2,10,n1) = arg(2,5,n1)
         arg(3,10,n1) = arg(3,5,n1)
c            derive s5 -> s51
         arg(1,11,n1) = arg(1,5,n1) - mh1**2 
         arg(2,11,n1) = arg(2,5,n1)
         arg(3,11,n1) = arg(3,5,n1)
c            derive s5 -> s52
         arg(1,12,n1) = arg(1,5,n1) - mh2**2 
         arg(2,12,n1) = arg(2,5,n1)
         arg(3,12,n1) = arg(3,5,n1)
c            move s3t, works like cutoff around s3t=0, s3t>0 ensured
c         if (theta_s3.eq.1.D0) then 
c            dummy = arg(1,8,n1) 
c            arg(1,8,n1)=sign(1.D0,dummy)*sqrt(dummy**2+gamt**2*mt**2) 
c         end if 
      end do

c              special case of reference frame (C)
c              only for theta_s3=1 and asked for s3 subtraction 
c              only for the integrals including t1t,u1t
      if (ndim.eq.1) then 
         s3  = massin(4)
         s3t = s3 + m2**2 - mt**2
         if (theta_s3.eq.1.D0)  then
            s3t = sign(1.D0,s3t) * sqrt( s3t**2 + mt**2*gamt**2 )
         end if
         c1 = ( 2.D0*w3*e2 - s3 )/(2.D0*p*w3)
         if ( (c1.lt.-1.D0).or.(c1.gt.1.D0) ) then 
            print*," MAND_TO_ANG_HH: problem with c1 ",c1 
            call HARD_STOP
         end if 
         do n1=0,12
            arg_x(n1) = arg(1,n1,3) + arg(2,n1,3)*c1
            arg_y(n1) = arg(3,n1,3)*sqrt(1.D0-c1**2)
         end do
      end if

c         check of the structure in debug mode 
      if ( (ndim.eq.2).and.(ldebug) ) then 
         do n1=1,3,1
            do n2=1,12,1
               if (abs(arg(3,n2,n1)).le.1.D-16) then
                  print*, " only a,b   : ",n1,n2,arg(3,n2,n1)
               end if
            end do
         end do

c         special case a=-b 
         do n1=1,3,1
            do n2=1,12,1
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
            do n2=1,12,1
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





































