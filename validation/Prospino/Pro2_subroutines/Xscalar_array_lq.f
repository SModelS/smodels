CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  TWO ROUTINES TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS C
C     FOR THE LEPTOQUARK PAIR HADROPRODUCTION                         C
C                                                                     C
C     SCALAR_ARRAY_B(MASSIN,SK3B0A,SK3B0B,SK3B0C,SK3B0D,SK3B0P)       C
C     SCALAR_ARRAY_C(MASSIN,SK3C0A,SK3C0B,SK3C0C,SK3D0)               C
C                                                                     C
C                                                                     C
C    INPUT  : MASSIN(1:2,1:6)                                         C
C                                                                     C
C             CONTAINING IN THE FIRST ROW  : S,T,U,T1,U1              C
C                               SECOND ROW : STOP-1 MASS SQUARED      C
C                                            [STOP-2 MASS SQUARED]    C
C                                            [SQUARK MASS SQUARED]    C
C                                            [GLUINO MASS SQUARED]    C
C                                            TOP MASS SQUARED         C
C                                            MASS SCALE SQUARED       C
C                                                                     C
C                                                                     C
C    OUTPUT : SK3B0A(1:5)         B_0(K1,M,M)                         C
C             SK3B0B(1:6)         B_0(K1+K2,M,M)                      C
C             SK3B0C(1:2)         B_0(P1,M1,M2)                       C
C             SK3B0D(1:2,1:2)     B_0(K1+P1,M1,M2)                    C
C             SK3B0P(1:3)         B_0'(K1/P1,M1,M2)                   C
C                                                                     C
C             SK3C0A(1:8)         C_0(K1,K2,M1,M2,M3)                 C
C             SK3C0B(1:4)         C_0(P1,P2,M1,M2,M3)                 C
C             SK3C0C(1:5,1:2)     C_0(K1,P1,M1,M2,M3)                 C
C             SK3D0(1:7,1:2)      D_0(........)                       C
C                                                                     C
C   N.B. : THE INTEGRAL MEASURE IS D^NQ/(I*Pi^2)                      C
C          ALL VARIABLES REAL*8, ONLY THE REAL PARTS OF THE INTEGRALS C
C                                                                     C
C                                                                     C
C    NECESSARY SUBROUTINES :  XSPENCE.F                               C
C                             XWOPOINT_2.F                            C
C                             XTHREEPOINT_2.F                         C
C                             XDENNER_2.F                             C
C                                                                     C
C    CALLED FUNCTIONS : B_FIN                                         C
C                       BP_FIN                                        C
C                       C_FIN                                         C
C                       D_FIN                                         C
C                       LIR2,LIC2                                     C
C                                                                     C
C                                                                     C
C    LAST MODIFIED : 26.03.97                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_LQ_B(m2,sba,sbb,sbc,sbd,sbp)
      
      implicit none 

      integer n1,n2
      real*8 m2(1:2,1:6)
      real*8 sba(1:5),sbb(1:6),sbc(1:2),sbd(1:2,1:2)
      real*8 sbp(1:3)

      real*8 m1(1:2,1:6)
      real*8 B02,BP02

      m2(1,6) = 1.D-16

      do n1=1,5
         sba(n1) = 1.D40 
      end do

      do n1=1,6 
         sbb(n1) = 1.D40
      end do

      do n1=1,2
         sbc(n1) = 1.D40
      end do
         
      do n1=1,2 
         do n2=1,2 
            sbd(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,3
         sbp(n1) = 1.D40
      end do

      do n1=1,6
         m1(2,n1) = sqrt(abs(m2(2,n1)))
      end do

c      sba(1) = B02(0.D0,m1(2,4),m1(2,4),m2(2,6)) 
      sba(2) = B02(0.D0,m1(2,5),m1(2,5),m2(2,6)) 
c      sba(3) = B02(0.D0,m1(2,3),m1(2,4),m2(2,6)) 
c      sba(4) = B02(0.D0,m1(2,3),m1(2,3),m2(2,6)) 
c      sba(5) = B02(0.D0,m1(2,2),m1(2,2),m2(2,6)) 

      sbb(1) = B02(m2(1,1),0.D0,0.D0,m2(2,6))
      sbb(2) = B02(m2(1,1),m1(2,1),m1(2,1),m2(2,6))
c      sbb(3) = B02(m2(1,1),m1(2,4),m1(2,4),m2(2,6))
      sbb(4) = B02(m2(1,1),m1(2,5),m1(2,5),m2(2,6))
c      sbb(5) = B02(m2(1,1),m1(2,3),m1(2,3),m2(2,6))
c      sbb(6) = B02(m2(1,1),m1(2,2),m1(2,2),m2(2,6))

      sbc(1) = B02(m2(2,1),m1(2,1),0.D0,m2(2,6))
c      sbc(2) = B02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))

      sbd(1,1) = B02(m2(1,2),m1(2,1),0.D0,m2(2,6))
c      sbd(2,1) = B02(m2(1,2),m1(2,4),m1(2,5),m2(2,6))

      sbd(1,2) = B02(m2(1,3),m1(2,1),0.D0,m2(2,6))
c      sbd(2,2) = B02(m2(1,3),m1(2,4),m1(2,5),m2(2,6))

      sbp(1) = BP02(m2(2,1),m1(2,1),0.D0,m2(2,6))
c      sbp(2) = BP02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))
c      sbp(3) = BP02(0.D0,m1(2,3),m1(2,4),m2(2,6))

      return 
      end 

c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_LQ_C(m2,sca,scb,scc,sd)
      
      implicit none 

      integer n1,n2 
      real*8 m2(1:2,1:6), s
      real*8 sca(1:8),scb(1:4),scc(1:5,1:2)
      real*8 sd(1:7,1:2)

      real*8 Lir2
c      real*8 C_fin_w,C_fin_d,D_fin,eps
      real*8 zeta2,xw1,xx1
      complex*16 cw1,cw5,cx1,cx5
      complex*16 CSPEN
c      complex*16 Lic2,cw2,cw3,cw4,cx2,cx3,cx4
      
c               real part of the spence function included in D04
      Lir2(s) = real( CSPEN(dcmplx(s)) )
c      Lic2(s) =       CSPEN(dcmplx(s))

c      eps = 1.D-4      
      
      do n1=1,8
         sca(n1) = 1.D40
      end do

      do n1=1,4 
         scb(n1) = 1.D40 
      end do

      do n1=1,5 
         do n2=1,2 
            scc(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,7
         do n2=1,2
            sd(n1,n2) = 1.D40 
         end do
      end do

      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0
      xw1 = sqrt(1.D0-4.D0*m2(2,1)/m2(1,1))
      xx1 = ( 1.D0 - xw1 )/( 1.D0 + xw1 )

      cw1 = sqrt( dcmplx( 1.D0-4.D0*m2(2,1)/m2(1,1)) )
      cx1 = ( 1.D0 - cw1 )/( 1.D0 + cw1 )
c      cw2 = sqrt( dcmplx( 1.D0-4.D0*m2(2,2)/m2(1,1)) )
c      cx2 = ( 1.D0 - cw2 )/( 1.D0 + cw2 )
c      cw3 = sqrt( dcmplx( 1.D0-4.D0*m2(2,3)/m2(1,1)) )
c      cx3 = ( 1.D0 - cw3 )/( 1.D0 + cw3 )
c      cw4 = sqrt( dcmplx( 1.D0-4.D0*m2(2,4)/m2(1,1)) )
c      cx4 = ( 1.D0 - cw4 )/( 1.D0 + cw4 )
      cw5 = sqrt( dcmplx( 1.D0-4.D0*m2(2,5)/m2(1,1)) )
      cx5 = ( 1.D0 - cw5 )/( 1.D0 + cw5 )

      scc(1,1) = (  log(-m2(1,4)/m2(2,1))**2 
     &            + Lir2(m2(1,2)/m2(2,1))
     &            + zeta2/4.D0              )/m2(1,4)
      scc(1,2) = (  log(-m2(1,5)/m2(2,1))**2 
     &            + Lir2(m2(1,3)/m2(2,1))
     &            + zeta2/4.D0              )/m2(1,5)

c   the three point functions with two vanishing external masses 
c   according to roland resp. wim
      sca(1) = (  log(m2(1,1)/m2(2,1))**2 /2.D0
     &          - 7.D0/2.D0*zeta2                )/m2(1,1)
      sca(2) = real( log(-cx1)**2 ) /(2.D0*m2(1,1))
c      sca(3) = real( log(-cx4)**2 ) /(2.D0*m2(1,1))
      sca(4) = real( log(-cx5)**2 ) /(2.D0*m2(1,1))
c      sca(5) = real(  Lic2(1.D0 + m2(2,4)*cx3/m2(2,3))
c     &               + Lic2(1.D0 + m2(2,4)/cx3/m2(2,3))
c     &               - 2.D0*Lir2(1.D0 - m2(2,4)/m2(2,3))
c     &               + log(-cx3)**2                   )/m2(1,1)
c      sca(6) = real(  Lic2(1.D0 + m2(2,3)*cx4/m2(2,4))
c     &               + Lic2(1.D0 + m2(2,3)/cx4/m2(2,4))
c     &               - 2.D0*Lir2(1.D0 - m2(2,3)/m2(2,4))
c     &               + log(-cx4)**2                   )/m2(1,1)
c      sca(7) = real( log(-cx3)**2 ) /(2.D0*m2(1,1))
c      sca(8) = real( log(-cx2)**2 ) /(2.D0*m2(1,1))

      scb(1) = (  - 2.D0*log(xx1)*log(1.D0-xx1)
     &            - 2.D0*Lir2(xx1) 
     &            + log(xx1)**2/2.D0
     &            - 4.D0*zeta2                  )/m2(1,1)/xw1
      scb(2) = (  2.D0 * Lir2(-xx1) 
     &          + log(xx1)**2 /2.D0 
     &          + zeta2              )/(m2(1,1)*xw1)
c      scb(3) = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,4),m2(2,5),m2(2,4))
c      scb(4) = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,5),m2(2,4),m2(2,5))

c      scc(2,1) = C_fin_d(m2(1,2),m2(2,1),eps,m2(2,5),m2(2,4),m2(2,3))
      scc(3,1) = ( - Lir2(m2(1,2)/m2(2,1)) + zeta2 )/m2(1,4)
c      scc(4,1) = C_fin_d(m2(1,2),m2(2,1),eps,m2(2,5),m2(2,4),m2(2,4))
c      scc(5,1) = C_fin_d(m2(1,2),m2(2,1),eps,m2(2,4),m2(2,5),m2(2,5))

c      scc(2,2) = C_fin_d(m2(1,3),m2(2,1),eps,m2(2,5),m2(2,4),m2(2,3))
      scc(3,2) = ( - Lir2(m2(1,3)/m2(2,1)) + zeta2 )/m2(1,5)
c      scc(4,2) = C_fin_d(m2(1,3),m2(2,1),eps,m2(2,5),m2(2,4),m2(2,4))
c      scc(5,2) = C_fin_d(m2(1,3),m2(2,1),eps,m2(2,4),m2(2,5),m2(2,5))

c   the special four point function according to prd40 
      sd(1,1) = (  2.D0*log(m2(1,1)/m2(2,1))*log(-m2(1,4)/m2(2,1)) 
     &           - 4.D0*zeta2        )/(m2(1,1)*m2(1,4))
      sd(2,1) = (- 2.D0*log(xx1)*log(1.D0-xx1) 
     &           + 2.D0*log(xx1)*log(1.D0+xx1) 
     &           - 2.D0*log(xx1)*log(-m2(1,4)/m2(2,1))
     &           - 2.D0*Lir2(xx1) 
     &           + 2.D0*Lir2(-xx1)
     &           - 3.D0*zeta2        )/(m2(1,1)*m2(1,4)*xw1)
      sd(3,1) = (  2.D0*log(-m2(1,4)/m2(2,1))*log(-m2(1,5)/m2(2,1))
     &           - 7.D0*zeta2/2.D0     )/(m2(1,4)*m2(1,5))

      sd(1,2) = (  2.D0*log(m2(1,1)/m2(2,1))*log(-m2(1,5)/m2(2,1)) 
     &           - 4.D0*zeta2        )/(m2(1,1)*m2(1,5))
      sd(2,2) = (- 2.D0*log(xx1)*log(1.D0-xx1) 
     &           + 2.D0*log(xx1)*log(1.D0+xx1) 
     &           - 2.D0*log(xx1)*log(-m2(1,5)/m2(2,1))
     &           - 2.D0*Lir2(xx1) 
     &           + 2.D0*Lir2(-xx1)
     &           - 3.D0*zeta2        )/(m2(1,1)*m2(1,5)*xw1)
      sd(3,2) = (  2.D0*log(-m2(1,5)/m2(2,1))*log(-m2(1,4)/m2(2,1))
     &           - 7.D0*zeta2/2.D0     )/(m2(1,5)*m2(1,4))

c   the general four point functions
c      sd(4,1) = D_fin(m2(2,1),eps,eps,m2(2,1),m2(1,2),m2(1,1),
c     &                m2(2,4),m2(2,5),m2(2,5),m2(2,5)) 
c      sd(5,1) = D_fin(eps,eps,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
c     &                m2(2,4),m2(2,4),m2(2,4),m2(2,5)) 
c      sd(6,1) = D_fin(eps,m2(2,1),eps,m2(2,1),m2(1,3),m2(1,2),
c     &                m2(2,5),m2(2,5),m2(2,4),m2(2,4)) 
c      sd(7,1) = D_fin(eps,eps,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
c     &                m2(2,4),m2(2,3),m2(2,4),m2(2,5)) 

c      sd(4,2) = D_fin(m2(2,1),eps,eps,m2(2,1),m2(1,3),m2(1,1),
c     &                m2(2,4),m2(2,5),m2(2,5),m2(2,5)) 
c      sd(5,2) = D_fin(eps,eps,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
c     &                m2(2,4),m2(2,4),m2(2,4),m2(2,5)) 
c      sd(6,2) = D_fin(eps,m2(2,1),eps,m2(2,1),m2(1,2),m2(1,3),
c     &                m2(2,5),m2(2,5),m2(2,4),m2(2,4)) 
c      sd(7,2) = D_fin(eps,eps,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
c     &                m2(2,4),m2(2,3),m2(2,4),m2(2,5)) 

      return 
      end 


c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c                                                              c
c$$$c all the other routines :                                     c
c$$$c                                                              c
c$$$c      B_fin(s,xm1**2,xm2**2,xmu**2)                           c
c$$$c      Bp_fin(s,xm1**2,xm2**2,xmu**2)                          c
c$$$c      C_fin(p0**2,p1**2,p2**2,m0**2,m1**2,m2**2)              c
c$$$c      D_fin(p1**2,p2**2,p3**2,p4,p12,p23,m12,m22,m32,m42)     c
c$$$c      Lic2(x),Lir2(x)                                         c
c$$$c                                                              c
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      function B_fin(s,m12,m22,mu2)
c$$$
c$$$      implicit none
c$$$      
c$$$      real*8 B_fin,s,m12,m22,mu2,m1,m2,B02
c$$$
c$$$      m1 = sqrt(m12) 
c$$$      m2 = sqrt(m22) 
c$$$
c$$$      B_fin = B02(s,m1,m2,mu2)
c$$$
c$$$      return
c$$$      end
c$$$
c$$$c---------------------------------------------------------------
c$$$      function Bp_fin(s,m12,m22,mu2)
c$$$
c$$$      implicit none
c$$$      
c$$$      real*8 Bp_fin,s,m12,m22,mu2,m1,m2,Bp02
c$$$
c$$$      m1 = sqrt(m12) 
c$$$      m2 = sqrt(m22) 
c$$$
c$$$      Bp_fin = BP02(s,m1,m2,mu2)
c$$$
c$$$      return
c$$$      end
