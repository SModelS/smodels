ctp 7/27 reverse the arguments of array arg
ctp n.b. no changes in Xangular_basic.f since there it is
ctp      three dimensional only 
c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_QG(massin,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,n5,ndim
      real*8 theta_s3
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG2_EXT,ANG2
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
c      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1:3,n4,n3),arg(1:3,n5,n3))
     
c            n=0 fixes the mandelstam variables to u7j/t1j,u6j/u1j
      n        = 0 
      ndim     = 2 
      theta_s3 = 0.D0
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then 
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c              normalization 
      NSANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]          -> (A)
      NSANG(1,0,-1, 0) = 0.D0
      NSANG(1,0,-2, 0) = ANG2(1,0,-2, 0,1)
c               tp[1], up[2]   -> (A)
      NSANG(1,2,-1,-1) = ANG2(1,2,-1,-1,1)
c               tp[1], s5[5]   -> (A)
      NSANG(1,5,-1,-1) = ANG2(1,5,-1,-1,1)
      NSANG(1,5,-1,-2) = ANG2(1,5,-1,-2,1)
      NSANG(1,5,-2,-1) = ANG2(1,5,-2,-1,1)
      NSANG(1,5,-2,-2) = ANG2(1,5,-2,-2,1)
c               tp[1], u7[7]   -> (A)
      NSANG(1,7,-1, 1) = ANG2(1,7,-1, 1,1)
      NSANG(1,7,-1, 2) = ANG2(1,7,-1, 2,1)
c               up[2]
      NSANG(2,0,-1, 0) = 0.D0
c               up[2], s5[5]   -> (B)      
      NSANG(2,5,-1,-1) = ANG2(2,5,-1,-1,2)
      NSANG(2,5,-1,-2) = ANG2(2,5,-1,-2,2)
c               up[2], u6[6]   -> (B)      
      NSANG(2,6,-1,-1) = ANG2(2,6,-1,-1,2)
      NSANG(2,6,-1,-2) = ANG2(2,6,-1,-2,2)
c               s5[5]          -> (C)
      NSANG(5,0, 1, 0) = ANG2(5,0, 1, 0,3)
      NSANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
      NSANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)
c               s5[5], u6[6]   -> (C)
      NSANG(5,6,-2, 1) = ANG2(5,6,-2, 1,3)
      NSANG(5,6,-1,-1) = ANG2(5,6,-1,-1,3)
      NSANG(5,6,-1, 1) = ANG2(5,6,-1, 1,3)
ctp      NSANG(5,6,-1, 2) = ANG2(5,6,-1, 2,3)
      NSANG(5,6,-1, 2) = ANG2(6,5, 2,-1,1)
c               u6[6]          -> (A)
      NSANG(6,0,-1, 0) = ANG2(6,0,-1, 0,1)
      NSANG(6,0,-2, 0) = ANG2(6,0,-2, 0,1)
c               u7[7]          -> (B)
      NSANG(7,0, 1, 0) = ANG2(7,0, 1, 0,2)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_GG(massin,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,n5,ndim
      real*8 theta_s3
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG2_EXT,ANG2
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
c      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1:3,n4,n3),arg(1:3,n5,n3))
     
c            n=0 fixes the mandelstam variables to u7j/t1j,u6j/u1j
      n        = 0 
      ndim     = 2 
      theta_s3 = 0.D0
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c              normalization 
      NSANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]
      NSANG(1,0,-1, 0) = 0.D0
c               tp[1], up[2]   -> (A)
      NSANG(1,2,-1,-1) = ANG2(1,2,-1,-1,1)
c               tp[1], s5[5]   -> (A)
      NSANG(1,5,-1,-1) = ANG2(1,5,-1,-1,1)
      NSANG(1,5,-1,-2) = ANG2(1,5,-1,-2,1)
c               tp[1], u6[6]   -> (A)
      NSANG(1,6,-1, 1) = ANG2(1,6,-1, 1,1)
c               tp[1], u7[7]   -> (A)
      NSANG(1,7,-1, 1) = ANG2(1,7,-1, 1,1)
      NSANG(1,7,-1,-1) = ANG2(1,7,-1,-1,1)
      NSANG(1,7,-1,-2) = ANG2(1,7,-1,-2,1)
c               tp[1], s3[8]   -> (A)
      NSANG(1,8,-1,-1) = ANG2(1,8,-1,-1,1)
c               up[2]
      NSANG(2,0,-1, 0) = 0.D0
c               up[2], s5[5]   -> (B)      
      NSANG(2,5,-1,-1) = ANG2(2,5,-1,-1,2)
      NSANG(2,5,-1,-2) = ANG2(2,5,-1,-2,2)
c               up[2], u6[6]   -> (B)      
      NSANG(2,6,-1, 1) = ANG2(2,6,-1, 1,2)
      NSANG(2,6,-1,-1) = ANG2(2,6,-1,-1,2)
      NSANG(2,6,-1,-2) = ANG2(2,6,-1,-2,2)
c               up[2], s3[8]   -> (B)      
      NSANG(2,8,-1,-1) = ANG2(2,8,-1,-1,2)
c               s5[5]          -> (C)
      NSANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
      NSANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)
c               s5[5], u6[6]   -> u6[6], s5[5] -> (A)
ctp      NSANG(5,6,-1, 2) = ANG2(5,6,-1, 2,3)
      NSANG(5,6,-1, 2) = ANG2(6,5, 2,-1,1)
c               s5[5], u6[6]   -> (C) 
      NSANG(5,6,-1, 1) = ANG2(5,6,-1, 1,3)
      NSANG(5,6,-1,-1) = ANG2(5,6,-1,-1,3)
      NSANG(5,6,-2, 1) = ANG2(5,6,-2, 1,3)
      NSANG(5,6,-2, 2) = ANG2(5,6,-2, 2,3)
c               s5[5], u7[7]   -> (C)
      NSANG(5,7,-1,-1) = ANG2(5,7,-1,-1,3)
      NSANG(5,7, 1,-1) = ANG2(5,7, 1,-1,3)
c               u6[6]          -> (A)
      NSANG(6,0,-1, 0) = ANG2(6,0,-1, 0,1)
      NSANG(6,0, 1, 0) = ANG2(6,0, 1, 0,1)
c               u6[6], u[7]    -> (A)
      NSANG(6,7,-1,-1) = ANG2(6,7,-1,-1,1)
c               u6[6], s3[8]   -> s3[8], u6[6] -> (C)
ctp      NSANG(6,8, 2,-1) = ANG2(8,6,-1, 2,3)
      NSANG(6,8, 2,-1) = ANG2(6,8, 2,-1,1)
      NSANG(6,8, 2,-2) = ANG2(8,6,-2, 2,3)
      NSANG(6,8, 1,-1) = ANG2(8,6,-1, 1,3)
      NSANG(6,8, 1,-2) = ANG2(8,6,-2, 1,3)
      NSANG(6,8,-1,-2) = ANG2(8,6,-2,-1,3)
      NSANG(6,8,-1,-1) = ANG2(8,6,-1,-1,3)
      NSANG(6,8,-2,-1) = ANG2(8,6,-1,-2,3)
      NSANG(6,8,-2,-2) = ANG2(8,6,-2,-2,3)
c               u7[7]          -> (B)
      NSANG(7,0,-1, 0) = ANG2(7,0,-1, 0,2)
c               u7[7], s3[8]   -> (B)
      NSANG(7,8,-1,-1) = ANG2(7,8,-1,-1,2)
      NSANG(7,8,-1,-2) = ANG2(7,8,-1,-2,2)
      NSANG(7,8,-2,-2) = ANG2(7,8,-2,-2,2)
      NSANG(7,8,-2,-1) = ANG2(7,8,-2,-1,2)
c               s3[8]          -> (C)
      NSANG(8,0,-1, 0) = ANG2(8,0,-1, 0,3)
      NSANG(8,0,-2, 0) = ANG2(8,0,-2, 0,3)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_GGOS(massin,theta_s3,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,ndim
      real*8 theta_s3,s3,s3t,mt,m2,gamt
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG1_EXT,ANG1
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG1(n1,n2)         =ANG1_EXT(n2,arg_x(n1),arg_y(n1))
     
c            n should not matter
      n        = 0
      ndim     = 1 
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      s3   = massin(4)
      mt   = massin(6)
      m2   = massin(7)
      gamt = massin(8)
      
      if (ldebug) then 
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c            the correct denominator for a numerical s3 integration   
      s3t = s3 + m2**2 - mt**2 
      s3t = sign(1.D0,s3t) * sqrt( s3t**2 + mt**2*gamt**2 )

c               n.b. no fixed cms for ANG1
      NSANG(6,8, 2,-2) = ANG1(6, 2)/s3t**2
      NSANG(6,8, 1,-2) = ANG1(6, 1)/s3t**2
      NSANG(6,8,-2,-2) = ANG1(6,-2)/s3t**2
      NSANG(6,8,-1,-2) = ANG1(6,-1)/s3t**2
      NSANG(7,8,-2,-2) = ANG1(7,-2)/s3t**2
      NSANG(7,8,-1,-2) = ANG1(7,-1)/s3t**2
c               s3[8]          -> (C)
      NSANG(8,0,-2, 0) = ANG1(0, 0)/s3t**2

      return
      end 
      
c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_QB(massin,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,n5,ndim
      real*8 theta_s3
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG2_EXT,ANG2
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
c      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1:3,n4,n3),arg(1:3,n5,n3))

c            n=0 fixes the mandelstam variables to u7j/t1j,u6j/u1j
      n        = 0 
      ndim     = 2 
      theta_s3 = 0.D0
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c              normalization 
      NSANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]
      NSANG(1,0,-1, 0) = 0.D0
      NSANG(1,0,-2, 0) = ANG2(1,0,-2, 0,1)
c               tp[1], s5[5]   -> (A)
      NSANG(1,5,-1,-1) = ANG2(1,5,-1,-1,1)
      NSANG(1,5,-1,-2) = ANG2(1,5,-1,-2,1)
      NSANG(1,5,-2,-1) = ANG2(1,5,-2,-1,1)
      NSANG(1,5,-2,-2) = ANG2(1,5,-2,-2,1)
c               tp[1], u6[6]   -> (A)
      NSANG(1,6,-1, 1) = ANG2(1,6,-1, 1,1)
c               tp[1], u7[7]   -> (A)
      NSANG(1,7,-1, 1) = ANG2(1,7,-1, 1,1)
c               tp[1], s3[8]   -> (A)
      NSANG(1,8,-1,-1) = ANG2(1,8,-1,-1,1)
c               s5[5]          -> (C)
      NSANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
      NSANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)
c               s5[5], u6[6]   -> u6[6], s5[5] -> (A)
      NSANG(5,6,-1, 2) = ANG2(6,5, 2,-1,1)
c               s5[5], u6[6]   -> (C) 
      NSANG(5,6,-1, 1) = ANG2(5,6,-1, 1,3)
      NSANG(5,6,-2, 1) = ANG2(5,6,-2, 1,3)
      NSANG(5,6,-2, 2) = ANG2(5,6,-2, 2,3)
c               u6[6]          -> (A)
      NSANG(6,0, 1, 0) = ANG2(6,0, 1, 0,1)
c               u6[6], u[7]    -> (A)
      NSANG(6,7,-1,-1) = ANG2(6,7,-1,-1,1)
c               u6[6], s3[8]   -> s3[8], u6[6] -> (C)
      NSANG(6,8, 1,-2) = ANG2(8,6,-2, 1,3)
      NSANG(6,8, 1,-1) = ANG2(8,6,-1, 1,3)
      NSANG(6,8, 2,-2) = ANG2(8,6,-2, 2,3)
      NSANG(6,8, 2,-1) = ANG2(6,8, 2,-1,1)
c               s3[8]          -> (C)
      NSANG(8,0,-1, 0) = ANG2(8,0,-1, 0,3)
      NSANG(8,0,-2, 0) = ANG2(8,0,-2, 0,3)

      return
      end

c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_QBOS(massin,theta_s3,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,ndim
      real*8 theta_s3,s3,s3t,mt,m2,gamt
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG1_EXT,ANG1
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG1(n1,n2)         =ANG1_EXT(n2,arg_x(n1),arg_y(n1))
     
c            n should not matter
      n        = 0
      ndim     = 1 
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      s3   = massin(4)
      mt   = massin(6)
      m2   = massin(7)
      gamt = massin(8)
      
      if (ldebug) then 
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c            the correct denominator for a numerical s3 integration   
      s3t = s3 + m2**2 - mt**2 
      s3t = sign(1.D0,s3t) * sqrt( s3t**2 + mt**2*gamt**2 )

c               n.b. no fixed cms for ANG1
      NSANG(0,0, 0, 0) = ANG1(0, 0)
      NSANG(6,8, 1,-2) = ANG1(6, 1)/s3t**2
      NSANG(6,8, 2,-2) = ANG1(6, 2)/s3t**2
c               s3[8]          -> (C)
      NSANG(8,0,-2, 0) = ANG1(0, 0)/s3t**2

      return
      end 
      
c ----------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_HT_QQ(massin,NSANG)

      implicit none

      integer n,n1,n2,n3,n4,n5,ndim
      real*8 theta_s3
      real*8 arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
      real*8 NSANG(0:9,0:9,-2:2,-2:2),massin(1:30)
      real*8 ANG2_EXT,ANG2
           
      logical ldebug 
      parameter( ldebug=.false. )
           
c            to make the conversion simple  
      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1,n4,n3),arg(1,n5,n3))
c      ANG2(n4,n5,n1,n2,n3)=ANG2_EXT(n1,n2,arg(1:3,n4,n3),arg(1:3,n5,n3))
     
c            n=0 fixes the mandelstam variables to u7j/t1j,u6j/u1j
      n        = 0 
      ndim     = 2 
      theta_s3 = 0.D0
      call MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)
      
      if (ldebug) then
         do n1=0,9,1
            do n2=0,9,1
               do n3=-2,2,1
                  do n4=-2,2,1
                     NSANG(n1,n2,n3,n4) = 9.99999D+99
                  end do
               end do
            end do
         end do
      end if 

c              normalization 
      NSANG(0,0, 0, 0) = ANG2(0,0, 0, 0,1)
c               tp[1]
      NSANG(1,0,-1, 0) = 0.D0
      NSANG(1,0,-2, 0) = ANG2(1,0,-2, 0,1)
c               tp[1], up[2]   -> (A)
      NSANG(1,2,-1,-1) = ANG2(1,2,-1,-1,1)
c               tp[1], s5[5]   -> (A)
      NSANG(1,5,-1,-1) = ANG2(1,5,-1,-1,1)
      NSANG(1,5,-1,-2) = ANG2(1,5,-1,-2,1)
      NSANG(1,5,-2,-1) = ANG2(1,5,-2,-1,1)
      NSANG(1,5,-2,-2) = ANG2(1,5,-2,-2,1)
c               tp[1], u6[6]   -> (A)
      NSANG(1,6,-1, 1) = ANG2(1,6,-1, 1,1)
c               tp[1], u7[7]   -> (A)
      NSANG(1,7,-1, 1) = ANG2(1,7,-1, 1,1)
c               up[2]          0> (B)
      NSANG(2,0,-1, 0) = 0.D0
      NSANG(2,0,-2, 0) = ANG2(2,0,-2, 0,2)
c               up[2], s5[5]   -> (B)      
      NSANG(2,5,-1,-1) = ANG2(2,5,-1,-1,2)
      NSANG(2,5,-1,-2) = ANG2(2,5,-1,-2,2)
      NSANG(2,5,-2,-1) = ANG2(2,5,-2,-1,2)
      NSANG(2,5,-2,-2) = ANG2(2,5,-2,-2,2)
c               up[2], s3[8]   -> (B)      
      NSANG(2,8,-1, 1) = ANG2(2,8,-1, 1,2)
c               s5[5]          -> (C)
      NSANG(5,0,-1, 0) = ANG2(5,0,-1, 0,3)
      NSANG(5,0,-2, 0) = ANG2(5,0,-2, 0,3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c      new?    me    wim     definition                                c
c                                                                      c
c      for qg and gg channel [n=0]                                     c
c           1  tp  =  tp  = (k2-k3)^2                                  c
c           2  up  =  up  = (k1-k3)^2                                  c
c           3  s3  = s3j  = (k3-p2)^2 - m2^2                           c
c        *  4  t1j = u7j  = (k1-p1)^2 - m1^2                           c
c           5  s5  =  s5  = (p1+p2)^2        is cm connected to s3[3]  c
c           6  u1  =  u6  = (k2-p1)^2 - m1^2 is cm connected to tp[1]  c
c           7  t1  =  u7  = (k1-p1)^2 - m1^2 is cm connected to up[2]  c
c        *  8  s3s =  s3  = (k3-p2)^2 - m1^2                           c
c        *  9  u1j =  u6j = (k2-p1)^2 - m1^2                           c
c                                                                      c
c      different for crossed channels [n=1]                            c
c        *  4  t1g =  u7g = (k1-p1)^2 - mg^2 is never needed           c
c        *  9  u1g =  u6g = (k1-p1)^2 - mg^2 is never needed           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MAND_TO_ANG_HT(n,ndim,massin,arg,theta_s3,arg_x,arg_y)

      integer n,n1,n2,n3,ndim
      real*8  massin(1:30),arg(1:3,0:9,1:3),arg_x(0:9),arg_y(0:9)
     &       ,theta_s3,s,t2,s4,m1,m2,mg,s3,s3s
     &       ,gams
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
      gams = massin(8)

c              the two outgoing masses set to m1,m2 [see my thesis]
      m12 = m1**2
      m22 = m2**2
      u1  = s4 - s - t2
      u2  = u1 + m1**2 - m2**2 
      s42 = s4 + m1**2 - m2**2  

c              define in the three reference frames (A,B,C):
c               (1) and calculate the angular polynomial for tp,up,s3
      norm = 2.D0*sqrt(s4+m12)
      w1  = (s + u2)                        /norm
      w2  = (s + t2)                        /norm
      w3  =  s4                             /norm
      e1  = (s4 + 2.D0*m12)                 /norm
      e2  = -(t2 + u2 + 2.D0*m22)           /norm
      p   = sqrt((t2 + u2)**2 - 4.D0*m22*s)/norm

c              reference frames (A,B)
      cpa =(t2*s42 -s*(u2+2.D0*m22))/(s+t2)/sqrt((t2+u2)**2-4.D0*m22*s)
      cpb =(u2*s42 -s*(t2+2.D0*m22))/(s+u2)/sqrt((u2+t2)**2-4.D0*m22*s)

      spa = 1.D0 - cpa**2
      if (spa.ge.0.D0) then 
         spa = sqrt(spa)
      else 
         print*," MAND_TO_ANG_HT: fixing problem, spa**2= ",spa
         spa = 1.D-8
      end if 

      spb = 1.D0 - cpb**2
      if (spb.ge.0.D0) then 
         spb = sqrt(spb)
      else 
         print*," MAND_TO_ANG_HT: fixing problem, spb**2= ",spb
         spb = 1.D-8
      end if 
        
c              reference frame (A): k2||z, cm(k3,p1)
c              two terms   :    tp,u6,u6s,u6j,u6g,...
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
c              two terms   :    up,u7,u7s,u7j,u7g,...
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

c            calculate all the other angular polynomials 
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
c            derive u6 -> u6j,u6g
c$$$         if (n.eq.0) then 
            arg(1,9,n1) = arg(1,6,n1) + m1**2 - m2**2 
c$$$         else if (n.eq.1) then 
c$$$            arg(1,9,n1) = arg(1,6,n1) + m1**2 - mg**2 
c$$$         end if 
         arg(2,9,n1) = arg(2,6,n1)
         arg(3,9,n1) = arg(3,6,n1)
c            derive u7 -> u7j,u7g
c$$$         if (n.eq.0) then 
            arg(1,4,n1) = arg(1,7,n1) + m1**2 - m2**2 
c$$$         else if (n.eq.1) then 
c$$$            arg(1,4,n1) = arg(1,7,n1) + m1**2 - mg**2 
c$$$         end if 
         arg(2,4,n1) = arg(2,7,n1)
         arg(3,4,n1) = arg(3,7,n1)
c            derive s3 -> s3s, includes regularization 
         arg(1,8,n1) = arg(1,3,n1) + m2**2 - m1**2
         arg(2,8,n1) = arg(2,3,n1)
         arg(3,8,n1) = arg(3,3,n1)
      end do

c              special case of reference frame (C)
c              only for theta_s3=1 and asked for s3 subtraction 
c              only for the integrals including u6s,u7s
      if (ndim.eq.1) then 
         s3  = massin(4)
         s3s = s3 + m2**2 - m1**2
         s3s  = sqrt( s3s**2 + m1**2*gams**2 )
         c1 = ( 2.D0*w3*e2 - s3 )/(2.D0*p*w3)
         do n1=0,9
            arg_x(n1) = arg(1,n1,3) + arg(2,n1,3)*c1
            arg_y(n1) = arg(3,n1,3)*sqrt(1.D0-c1**2)
         end do
      end if


c         check of the structure in debug mode 
      if ( (ndim.eq.2).and.(ldebug) ) then 
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

