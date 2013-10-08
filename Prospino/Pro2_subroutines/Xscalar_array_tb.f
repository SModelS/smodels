CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C ROUTINES TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS      C
C     FOR THE STOP PAIR HADROPRODUCTION WITH QUARK INCOMING STATE     C
C                                                                     C
C     SCALAR_ARRAY_TB_B_GS(MASSIN,SK3B0A,SK3B0B,SK3B0C,SK3B0D,SK3B0P) C
C     SCALAR_ARRAY_TB_C_GS(MASSIN,SK3C0A,SK3C0B,SK3C0C,SK3D0)         C
C     SCALAR_ARRAY_TB_B_QS(MASSIN,SK3B0A,SK3B0B,SK3B0C,SK3B0D,SK3B0P) C
C     SCALAR_ARRAY_TB_C_QS(MASSIN,SK3C0A,SK3C0B,SK3C0C,SK3D0)         C
C                                                                     C
C COMMENTED OUT IF FUNCTION IS NOT NEEDED FOR GIVEN INITIAL STATE     C
C                                                                     C
C                                                                     C
C    INPUT  : MASSIN(1:2,1:6)                                         C
C                                                                     C
C             CONTAINING IN THE FIRST ROW  : S,T,U,T1,U1              C
C                               SECOND ROW : STOP-1 MASS SQUARED      C
C                                            STOP-2 MASS SQUARED      C
C                                            SQUARK MASS SQUARED      C
C                                            GLUINO MASS SQUARED      C
C                                            TOP MASS SQUARED         C
C                                            MASS SCALE SQUARED       C
C                                            SBOTTOM-1 MASS SQUARED   C
C                                            SBOTTOM-2 MASS SQUARED   C
C                                                                     C
C                                                                     C
C    OUTPUT : SK3B0A(1:7)         B_0(K1,M,M)                         C
C             SK3B0B(1:8)         B_0(K1+K2,M,M)                      C
C             SK3B0C(1:3)         B_0(P1,M1,M2)                       C
C             SK3B0D(1:3,1:2)     B_0(K1+P1,M1,M2)                    C
C             SK3B0P(1:4)         B_0'(K1/P1,M1,M2)                   C
C                                                                     C
C             SK3C0A(1:10)        C_0(K1,K2,M1,M2,M3)                 C
C             SK3C0B(1:6)         C_0(P1,P2,M1,M2,M3)                 C
C             SK3C0C(1:8,1:2)     C_0(K1,P1,M1,M2,M3)                 C
C             SK3D0(1:11,1:2)     D_0(........)                       C
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
C    LAST MODIFIED : 28.01.08 (SBOTTOM INTEGRALS)                     C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_TB_B_QS(m2,sba,sbb,sbc,sbd,sbp)
      
      implicit none 

      integer n1,n2
      real*8 m2(1:2,1:8)
      real*8 sba(1:7),sbb(1:8),sbc(1:3),sbd(1:3,1:2)
      real*8 sbp(1:4)

      real*8 m1(1:2,1:8)
      real*8 B02,BP02

      m2(1,6) = 1.D-16

      do n1=1,7
         sba(n1) = 1.D40 
      end do

      do n1=1,8 
         sbb(n1) = 1.D40
      end do

      do n1=1,3
         sbc(n1) = 1.D40
      end do
         
      do n1=1,3 
         do n2=1,2 
            sbd(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,4
         sbp(n1) = 1.D40
      end do

      do n1=1,8
         m1(2,n1) = sqrt(abs(m2(2,n1)))
      end do

      sba(1) = B02(0.D0,m1(2,4),m1(2,4),m2(2,6)) 
      sba(2) = B02(0.D0,m1(2,5),m1(2,5),m2(2,6)) 
      sba(3) = B02(0.D0,m1(2,3),m1(2,4),m2(2,6)) 
      sba(4) = B02(0.D0,m1(2,3),m1(2,3),m2(2,6)) 
      sba(5) = B02(0.D0,m1(2,2),m1(2,2),m2(2,6)) 
      sba(6) = B02(0.D0,m1(2,7),m1(2,7),m2(2,6)) 
      sba(7) = B02(0.D0,m1(2,8),m1(2,8),m2(2,6)) 

      sbb(1) = B02(m2(1,1),0.D0,0.D0,m2(2,6))
      sbb(2) = B02(m2(1,1),m1(2,1),m1(2,1),m2(2,6))
      sbb(3) = B02(m2(1,1),m1(2,4),m1(2,4),m2(2,6))
      sbb(4) = B02(m2(1,1),m1(2,5),m1(2,5),m2(2,6))
      sbb(5) = B02(m2(1,1),m1(2,3),m1(2,3),m2(2,6))
      sbb(6) = B02(m2(1,1),m1(2,2),m1(2,2),m2(2,6))
      sbb(7) = B02(m2(1,1),m1(2,7),m1(2,7),m2(2,6))
      sbb(8) = B02(m2(1,1),m1(2,8),m1(2,8),m2(2,6))

      sbc(1) = B02(m2(2,1),m1(2,1),0.D0,m2(2,6))
      sbc(2) = B02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))
      sbc(3) = B02(m2(2,1),m1(2,4),0.D0   ,m2(2,6))

c      sbd(1,1) = B02(m2(1,2),m1(2,1),0.D0,   m2(2,6))
c      sbd(2,1) = B02(m2(1,2),m1(2,4),m1(2,5),m2(2,6))
c      sbd(3,1) = B02(m2(1,2),m1(2,4),0.D0,   m2(2,6))

c      sbd(1,2) = B02(m2(1,3),m1(2,1),0.D0,   m2(2,6))
c      sbd(2,2) = B02(m2(1,3),m1(2,4),m1(2,5),m2(2,6))
c      sbd(3,2) = B02(m2(1,3),m1(2,4),0.D0,   m2(2,6))

      sbp(1) = BP02(m2(2,1),m1(2,1),0.D0,m2(2,6))
      sbp(2) = BP02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))
      sbp(3) = BP02(0.D0,m1(2,3),m1(2,4),m2(2,6))
      sbp(4) = BP02(m2(2,1),m1(2,4),0.D0   ,m2(2,6))

      return 
      end 

c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_TB_C_QS(m2,sca,scb,scc,sd)
      
      implicit none 

      integer n1,n2 
      real*8 m2(1:2,1:8)
      real*8 sca(1:10),scb(1:6),scc(1:8,1:2)
      real*8 sd(1:11,1:2)

      real*8 m1(1:2,1:8)
      real*8 Lir2,eps1
      real*8 zeta2,xw1,xx1,s
ctp      real*8 dummy1, dummy2
      complex*16 cw1,cw2,cw3,cw4,cw5,cx1,cx2,cx3,cx4,cx5
      complex*16 C03,D04,CSPEN

      Lir2(s) = real( CSPEN(dcmplx(s)) )

      eps1 = 1.D-8

      do n1=1,10
         sca(n1) = 1.D40
      end do

      do n1=1,6 
         scb(n1) = 1.D40 
      end do

      do n1=1,8 
         do n2=1,2 
            scc(n1,n2) = 1.D40 
         end do
      end do
            
      do n1=1,11
         do n2=1,2
            sd(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,8
         m1(2,n1) = sqrt(abs(m2(2,n1)))
      end do

      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0
      xw1 = sqrt(1.D0-4.D0*m2(2,1)/m2(1,1))
      xx1 = ( 1.D0 - xw1 )/( 1.D0 + xw1 )

      cw1 = sqrt( dcmplx( 1.D0-4.D0*m2(2,1)/m2(1,1)) )
      cx1 = ( 1.D0 - cw1 )/( 1.D0 + cw1 )
      cw2 = sqrt( dcmplx( 1.D0-4.D0*m2(2,2)/m2(1,1)) )
      cx2 = ( 1.D0 - cw2 )/( 1.D0 + cw2 )
      cw3 = sqrt( dcmplx( 1.D0-4.D0*m2(2,3)/m2(1,1)) )
      cx3 = ( 1.D0 - cw3 )/( 1.D0 + cw3 )
      cw4 = sqrt( dcmplx( 1.D0-4.D0*m2(2,4)/m2(1,1)) )
      cx4 = ( 1.D0 - cw4 )/( 1.D0 + cw4 )
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
c      sca(2) = real( log(-cx1)**2 ) /(2.D0*m2(1,1))
c      sca(3) = real( log(-cx4)**2 ) /(2.D0*m2(1,1))
c      sca(4) = real( log(-cx5)**2 ) /(2.D0*m2(1,1))
      sca(5) = real(  CSPEN(1.D0 + m2(2,4)*cx3/m2(2,3))
     &               + CSPEN(1.D0 + m2(2,4)/cx3/m2(2,3))
     &               - 2.D0*Lir2(1.D0 - m2(2,4)/m2(2,3))
     &               + log(-cx3)**2                   )/m2(1,1)
      sca(6) = real(  CSPEN(1.D0 + m2(2,3)*cx4/m2(2,4))
     &               + CSPEN(1.D0 + m2(2,3)/cx4/m2(2,4))
     &               - 2.D0*Lir2(1.D0 - m2(2,3)/m2(2,4))
     &               + log(-cx4)**2                   )/m2(1,1)
c      sca(7) = real( log(-cx3)**2 ) /(2.D0*m2(1,1))
c      sca(8) = real( log(-cx2)**2 ) /(2.D0*m2(1,1))

      scb(1) = (  - 2.D0*log(xx1)*log(1.D0-xx1)
     &            - 2.D0*Lir2(xx1) 
     &            + log(xx1)**2/2.D0
     &            - 4.D0*zeta2                  )/m2(1,1)/xw1
      scb(2) = (  2.D0 * Lir2(-xx1) 
     &          + log(xx1)**2 /2.D0 
     &          + zeta2              )/(m2(1,1)*xw1)
      scb(3)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,4),m1(2,5),m1(2,4)))
      scb(4)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,5),m1(2,4),m1(2,5)))
      scb(5)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,4),0.D0   ,m1(2,4)))
      scb(6)=real(C03(m2(2,1),m2(2,1),m2(1,1),0.D0   ,m1(2,4),0.D0   ))

ctp      dummy1 = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,4),m2(2,5),m2(2,4))
ctp      dummy2 = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,5),m2(2,4),m2(2,5))
ctp      if ( abs(scb(3)-dummy1)/abs(scb(3)+dummy1) .gt. 1.e-6 ) then
ctp         print*, " C_QS: problem 1 ",dummy1,scb(3)
ctp      end if
ctp      if ( abs(scb(4)-dummy2)/abs(scb(4)+dummy2) .gt. 1.e-6 ) then
ctp         print*, " C_QS: problem 2 ",dummy2,scb(4)
ctp      end if

ctp      scc(2,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,3)))
      scc(2,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,3),m1(2,5),m1(2,4)))
      scc(6,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,3),0.D0   ,m1(2,4)))
c      scc(3,1) = ( - Lir2(m2(1,2)/m2(2,1)) + zeta2 )/m2(1,4)
ctp      scc(4,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,4)))
ctp      scc(5,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,5)))
c      scc(4,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,4),m1(2,3),m1(2,4)))
c      scc(5,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,5)))

ctp      scc(2,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,3)))
      scc(2,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,3),m1(2,5),m1(2,4)))
      scc(6,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,3),0.D0   ,m1(2,4)))
c      scc(3,2) = ( - Lir2(m2(1,3)/m2(2,1)) + zeta2 )/m2(1,5)
ctp      scc(4,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,4)))
ctp      scc(5,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,5)))
c      scc(4,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,4)))
c      scc(5,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,5)))

c   the special four point function according to prd40 
      sd(1,1) = (  2.D0*log(m2(1,1)/m2(2,1))*log(-m2(1,4)/m2(2,1)) 
     &           - 4.D0*zeta2        )/(m2(1,1)*m2(1,4))
c      sd(2,1) = (- 2.D0*log(xx1)*log(1.D0-xx1) 
c     &           + 2.D0*log(xx1)*log(1.D0+xx1) 
c     &           - 2.D0*log(xx1)*log(-m2(1,4)/m2(2,1))
c     &           - 2.D0*Lir2(xx1) 
c     &           + 2.D0*Lir2(-xx1)
c     &           - 3.D0*zeta2        )/(m2(1,1)*m2(1,4)*xw1)
c      sd(3,1) = (  2.D0*log(-m2(1,4)/m2(2,1))*log(-m2(1,5)/m2(2,1))
c     &           - 7.D0*zeta2/2.D0     )/(m2(1,4)*m2(1,5))

      sd(1,2) = (  2.D0*log(m2(1,1)/m2(2,1))*log(-m2(1,5)/m2(2,1)) 
     &           - 4.D0*zeta2        )/(m2(1,1)*m2(1,5))
c      sd(2,2) = (- 2.D0*log(xx1)*log(1.D0-xx1) 
c     &           + 2.D0*log(xx1)*log(1.D0+xx1) 
c     &           - 2.D0*log(xx1)*log(-m2(1,5)/m2(2,1))
c     &           - 2.D0*Lir2(xx1) 
c     &           + 2.D0*Lir2(-xx1)
c     &           - 3.D0*zeta2        )/(m2(1,1)*m2(1,5)*xw1)
c      sd(3,2) = (  2.D0*log(-m2(1,5)/m2(2,1))*log(-m2(1,4)/m2(2,1))
c     &           - 7.D0*zeta2/2.D0     )/(m2(1,5)*m2(1,4))

c   the general four point functions
c      sd(4,1) = real(D04(m2(2,1),eps1,eps1,m2(2,1),m2(1,2),m2(1,1),
c     &                   m1(2,4),m1(2,5),m1(2,5),m1(2,5)))
c      sd(5,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
c     &                   m1(2,4),m1(2,4),m1(2,4),m1(2,5)))
c      sd(6,1) = real(D04(eps1,m2(2,1),eps1,m2(2,1),m2(1,3),m2(1,2),
c     &                   m1(2,5),m1(2,5),m1(2,4),m1(2,4)))
      sd(7,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
     &                   m1(2,4),m1(2,3),m1(2,4),m1(2,5)))
      sd(8,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
     &                   m1(2,4),m1(2,3),m1(2,4),0.D0   ))

c      sd(4,2) = real(D04(m2(2,1),eps1,eps1,m2(2,1),m2(1,3),m2(1,1),
c     &                   m1(2,4),m1(2,5),m1(2,5),m1(2,5))) 
c      sd(5,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
c     &                   m1(2,4),m1(2,4),m1(2,4),m1(2,5)))
c      sd(6,2) = real(D04(eps1,m2(2,1),eps1,m2(2,1),m2(1,2),m2(1,3),
c     &                   m1(2,5),m1(2,5),m1(2,4),m1(2,4)))
      sd(7,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
     &                   m1(2,4),m1(2,3),m1(2,4),m1(2,5)))
      sd(8,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
     &                   m1(2,4),m1(2,3),m1(2,4),0.D0   ))

      return 
      end 

c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_TB_B_GS(m2,sba,sbb,sbc,sbd,sbp)
      
      implicit none 

      integer n1,n2
      real*8 m2(1:2,1:8)
      real*8 sba(1:7),sbb(1:8),sbc(1:3),sbd(1:3,1:2)
      real*8 sbp(1:4)

      real*8 m1(1:2,1:8)
      real*8 B02,BP02

      m2(1,6) = 1.D-16

      do n1=1,7
         sba(n1) = 1.D40       
      end do

      do n1=1,8 
         sbb(n1) = 1.D40
      end do

      do n1=1,3
         sbc(n1) = 1.D40
      end do

      do n1=1,3 
         do n2=1,2 
            sbd(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,4
         sbp(n1) = 1.D40
      end do

      do n1=1,8
         m1(2,n1) = sqrt(abs(m2(2,n1)))
      end do

      sba(1) = B02(0.D0,m1(2,4),m1(2,4),m2(2,6)) 
      sba(2) = B02(0.D0,m1(2,5),m1(2,5),m2(2,6)) 
c      sba(3) = B02(0.D0,m1(2,3),m1(2,4),m2(2,6)) 
      sba(4) = B02(0.D0,m1(2,3),m1(2,3),m2(2,6)) 
      sba(5) = B02(0.D0,m1(2,2),m1(2,2),m2(2,6)) 
      sba(6) = B02(0.D0,m1(2,7),m1(2,7),m2(2,6)) 
      sba(7) = B02(0.D0,m1(2,8),m1(2,8),m2(2,6)) 

      sbb(1) = B02(m2(1,1),0.D0,0.D0,m2(2,6))
      sbb(2) = B02(m2(1,1),m1(2,1),m1(2,1),m2(2,6))
      sbb(3) = B02(m2(1,1),m1(2,4),m1(2,4),m2(2,6))
      sbb(4) = B02(m2(1,1),m1(2,5),m1(2,5),m2(2,6))
      sbb(5) = B02(m2(1,1),m1(2,3),m1(2,3),m2(2,6))
      sbb(6) = B02(m2(1,1),m1(2,2),m1(2,2),m2(2,6))
      sbb(7) = B02(m2(1,1),m1(2,7),m1(2,7),m2(2,6))
      sbb(8) = B02(m2(1,1),m1(2,8),m1(2,8),m2(2,6))

      sbc(1) = B02(m2(2,1),m1(2,1),0.D0,m2(2,6))
      sbc(2) = B02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))
      sbc(3) = B02(m2(2,1),m1(2,4),0.D0,m2(2,6))

      sbd(1,1) = B02(m2(1,2),m1(2,1),0.D0,   m2(2,6))
      sbd(2,1) = B02(m2(1,2),m1(2,4),m1(2,5),m2(2,6))
      sbd(3,1) = B02(m2(1,2),m1(2,4),0.D0,   m2(2,6))

      sbd(1,2) = B02(m2(1,3),m1(2,1),0.D0,   m2(2,6))
      sbd(2,2) = B02(m2(1,3),m1(2,4),m1(2,5),m2(2,6))
      sbd(3,2) = B02(m2(1,3),m1(2,4),0.D0,   m2(2,6))

      sbp(1) = BP02(m2(2,1),m1(2,1),0.D0,m2(2,6))
      sbp(2) = BP02(m2(2,1),m1(2,4),m1(2,5),m2(2,6))
c      sbp(3) = BP02(0.D0,m1(2,3),m1(2,4),m2(2,6))
      sbp(4) = BP02(m2(2,1),m1(2,4),0.D0,m2(2,6))

      return 
      end 

c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_TB_C_GS(m2,sca,scb,scc,sd)
      
      implicit none 

      integer n1,n2 
      real*8 m2(1:2,1:8)
      real*8 sca(1:10),scb(1:6),scc(1:8,1:2)
      real*8 sd(1:11,1:2)

      real*8 m1(1:2,1:6)
      real*8 Lir2,t,m12,ms2,mu2,mg2,dumc_1,dumd_1,dumd_2,eps1,epsi
      real*8 tg,ug,mdiff
      real*8 zeta2,xw1,xx1,s
ctp      real*8 dummy1, dummy2
      complex*16 cw1,cw2,cw3,cw4,cw5,cx1,cx2,cx3,cx4,cx5
      complex*16 cw7,cw8,cx7,cx8,msc2,mgc2,sc,mdiffc
      complex*16 C03,D04,CSPEN
ctp      complex*16 Cc
      
      Lir2(s) = real( CSPEN(dcmplx(s)) )
c               divergent three point functions, roland's (0.34) 
c               copied from neutralino-gluino case:
c               SCC(2,1) = C_FIN(P1,K1,MS,0,0) 
c                        = dumc_1(t,m12,ms2,msc2,mu2)
      dumc_1(t,m12,ms2,msc2,mu2) =  
     &  (   Lir2(t/ms2) - Lir2(m12/ms2)
     &    + real( log( 1.D0 - t/msc2     )**2 )  
     &    - real( log( 1.D0 - m12/msc2 )**2 )
     &    + log(ms2/mu2) * log( abs(-(t-ms2)/(ms2-m12) ) ) )
     &  /(t-m12)

c               divergent four point function, roland's (0.41) gen. 
c               copied from neutralino-gluino case:
c               SCD(1,2) = D_FIN(K2,K1,P2,0,0,0,MS) 
c                        = dumd_1(u,mg2,m12,ms2,msc2,s,sc,mu2,zeta2)
      dumd_1(t,m12,mg2,ms2,msc2,s,sc,mu2,zeta2) =
     &  ( - 2.D0 * Lir2( 1.D0 + (ms2-m12)/(t-ms2) )
     &    - 2.D0 * Lir2( 1.D0 + (ms2-mg2)/(t-ms2) )
     &    - Lir2( 1.D0 + (ms2-m12)*(ms2-mg2)/s/ms2 )
     &    - 3.D0/2.D0 * zeta2
     &    - real( log( 1.D0 + (msc2-m12)*(msc2-mg2)/s/msc2 )
     &             * (   log( -(msc2-m12)*(msc2-mg2)/s/msc2 )
     &                 - log( (msc2-m12)/mu2)
     &                 - log( (msc2-mg2)/mu2)
     &                 + log( -s*msc2/mu2**2)         )              )
     &    + log( s/mu2  )**2 /2.D0
     &    - log( s/ms2  )**2 /2.D0
     &    + 2.D0*real( log( -sc/mu2 )
     &                 *log( -(t-msc2)/msc2) )
     &    - real( log( (msc2-m12)/mu2 )
     &           *log( (msc2-m12)/msc2) )
     &    - real( log( (msc2-mg2)/mu2 )
     &           *log( (msc2-mg2)/msc2) )                        )
     &  /s/(t-ms2)
     
      dumd_2(tg,ug,mg2,ms2,mdiff,mdiffc,mu2) =
     &  -1.D0/(tg*ug -mdiff**2) *(
     &             log(mdiff**2/mg2/mu2)* log( abs(mdiff**2/tg/ug) )
     &           + real( log( -mdiffc/tg )**2 )
     &           + real( log( -mdiffc/ug )**2 )
     &           + 4.d0*Lir2(1.D0 +mdiff/tg) 
     &           + 4.d0*Lir2(1.D0 +mdiff/ug)
     &           + 2.d0*Lir2(1.D0 -tg*ug/mdiff**2) )

      eps1 = 1.D-8
      epsi = 1.D-8

      mgc2 = m2(2,4) * dcmplx(1.D0,-epsi) 
      sc   = m2(1,1) * dcmplx(1.D0,epsi) 

      mdiff  = m2(2,4) - m2(2,1)
      mdiffc = mdiff * dcmplx(1.D0,epsi)
      tg    = m2(1,2) - m2(2,4)
      ug    = m2(1,3) - m2(2,4)
      
      do n1=1,10
         sca(n1) = 1.D40
      end do

      do n1=1,6 
         scb(n1) = 1.D40 
      end do

      do n1=1,8 
         do n2=1,2 
            scc(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,11
         do n2=1,2
            sd(n1,n2) = 1.D40 
         end do
      end do

      do n1=1,6
         m1(2,n1) = sqrt(abs(m2(2,n1)))
      end do

      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0
      xw1 = sqrt(1.D0-4.D0*m2(2,1)/m2(1,1))
      xx1 = ( 1.D0 - xw1 )/( 1.D0 + xw1 )

      cw1 = sqrt( dcmplx( 1.D0-4.D0*m2(2,1)/m2(1,1)) )
      cx1 = ( 1.D0 - cw1 )/( 1.D0 + cw1 )
      cw2 = sqrt( dcmplx( 1.D0-4.D0*m2(2,2)/m2(1,1)) )
      cx2 = ( 1.D0 - cw2 )/( 1.D0 + cw2 )
      cw3 = sqrt( dcmplx( 1.D0-4.D0*m2(2,3)/m2(1,1)) )
      cx3 = ( 1.D0 - cw3 )/( 1.D0 + cw3 )
      cw4 = sqrt( dcmplx( 1.D0-4.D0*m2(2,4)/m2(1,1)) )
      cx4 = ( 1.D0 - cw4 )/( 1.D0 + cw4 )
      cw5 = sqrt( dcmplx( 1.D0-4.D0*m2(2,5)/m2(1,1)) )
      cx5 = ( 1.D0 - cw5 )/( 1.D0 + cw5 )
      cw7 = sqrt( dcmplx( 1.D0-4.D0*m2(2,7)/m2(1,1)) )
      cx7 = ( 1.D0 - cw7 )/( 1.D0 + cw7 )
      cw8 = sqrt( dcmplx( 1.D0-4.D0*m2(2,8)/m2(1,1)) )
      cx8 = ( 1.D0 - cw8 )/( 1.D0 + cw8 )

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
      sca(3) = real( log(-cx4)**2 ) /(2.D0*m2(1,1))
      sca(4) = real( log(-cx5)**2 ) /(2.D0*m2(1,1))
c      sca(5) = real(  CSPEN(1.D0 + m2(2,4)*cx3/m2(2,3))
c     &               + CSPEN(1.D0 + m2(2,4)/cx3/m2(2,3))
c     &               - 2.D0*Lir2(1.D0 - m2(2,4)/m2(2,3))
c     &               + log(-cx3)**2                   )/m2(1,1)
c      sca(6) = real(  CSPEN(1.D0 + m2(2,3)*cx4/m2(2,4))
c     &               + CSPEN(1.D0 + m2(2,3)/cx4/m2(2,4))
c     &               - 2.D0*Lir2(1.D0 - m2(2,3)/m2(2,4))
c     &               + log(-cx4)**2                   )/m2(1,1)
      sca(7) = real( log(-cx3)**2 ) /(2.D0*m2(1,1))
      sca(8) = real( log(-cx2)**2 ) /(2.D0*m2(1,1))
      sca(9)  = real( log(-cx7)**2 ) /(2.D0*m2(1,1))
      sca(10) = real( log(-cx8)**2 ) /(2.D0*m2(1,1))

      scb(1) = (  - 2.D0*log(xx1)*log(1.D0-xx1)
     &            - 2.D0*Lir2(xx1) 
     &            + log(xx1)**2/2.D0
     &            - 4.D0*zeta2                  )/m2(1,1)/xw1
      scb(2) = (  2.D0 * Lir2(-xx1) 
     &          + log(xx1)**2 /2.D0 
     &          + zeta2              )/(m2(1,1)*xw1)
      scb(3)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,4),m1(2,5),m1(2,4)))
      scb(4)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,5),m1(2,4),m1(2,5)))
      scb(5)=real(C03(m2(2,1),m2(2,1),m2(1,1),m1(2,4),0.D0   ,m1(2,4)))
      scb(6)=real(C03(m2(2,1),m2(2,1),m2(1,1),0.D0   ,m1(2,4),0.D0   ))

ctp      dummy1 = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,4),m2(2,5),m2(2,4))
ctp      dummy2 = C_fin_w(m2(1,1),m2(2,1),m2(2,1),m2(2,5),m2(2,4),m2(2,5))
ctp      if ( abs(scb(3)-dummy1)/abs(scb(3)+dummy1) .gt. 1.e-6 ) then
ctp         print*, " C_QS: problem 1 ",dummy1,scb(3)
ctp      end if
ctp      if ( abs(scb(4)-dummy2)/abs(scb(4)+dummy2) .gt. 1.e-6 ) then
ctp         print*, " C_QS: problem 2 ",dummy2,scb(4)
ctp      end if

ctp      scc(2,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,3)))
c      scc(2,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,3),m1(2,5),m1(2,4)))
      scc(3,1) = ( - Lir2(m2(1,2)/m2(2,1)) + zeta2 )/m2(1,4)
ctp      scc(4,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,4)))
ctp      scc(5,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,5)))
      scc(4,1)=real(C03(m2(1,2),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,4)))
      scc(5,1)=real(C03(m2(1,2),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,5)))
ctp      scc(6,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,3),0.D0   ,m1(2,4)))
      scc(7,1) = real(C03(m2(1,2),m2(2,1),eps1,m1(2,4),0.D0   ,m1(2,4)))
c     Cfin(p1,k1,mg,0,0)
      scc(8,1) = dumc_1(m2(1,2),m2(2,1),m2(2,4),mgc2,m2(2,6))

ctp      scc(2,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,3)))
c      scc(2,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,3),m1(2,5),m1(2,4)))
      scc(3,2) = ( - Lir2(m2(1,3)/m2(2,1)) + zeta2 )/m2(1,5)
ctp      scc(4,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,4)))
ctp      scc(5,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,5)))
      scc(4,2)=real(C03(m2(1,3),m2(2,1),eps1,m1(2,4),m1(2,5),m1(2,4)))
      scc(5,2)=real(C03(m2(1,3),m2(2,1),eps1,m1(2,5),m1(2,4),m1(2,5)))
ctp      scc(6,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,3),0.D0   ,m1(2,4)))
      scc(7,2) = real(C03(m2(1,3),m2(2,1),eps1,m1(2,4),0.D0   ,m1(2,4)))
c     Cfin(p1,k1,mg,0,0)
      scc(8,2) = dumc_1(m2(1,3),m2(2,1),m2(2,4),mgc2,m2(2,6))

C  from old stop code: C(p02,p12,p22;m0,m1,m2)                        C
C                                                                     C
C             p1  -----\_  m1                                         C
C                      | \_                                           C
C                      |   \_                                         C
C                  m2  |    _----  p0 = p1+p2                         C
C                      |  _/                                          C
C                      |_/                                            C
C             p2  -----/   m0                                         C
C                                                                     C
C called as   C03(p02,p12,p22;m0,m1,m2)                               C
C             Cc (p02,p12,p22;m1,m2,m0)                               C
C                                                                     C

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
      sd(4,1) = real(D04(m2(2,1),eps1,eps1,m2(2,1),m2(1,2),m2(1,1),
     &                m1(2,4),m1(2,5),m1(2,5),m1(2,5)))
      sd(5,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
     &                m1(2,4),m1(2,4),m1(2,4),m1(2,5)))
      sd(6,1) = real(D04(eps1,m2(2,1),eps1,m2(2,1),m2(1,3),m2(1,2),
     &                m1(2,5),m1(2,5),m1(2,4),m1(2,4)))
c      sd(7,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
c     &                m1(2,4),m1(2,3),m1(2,4),m1(2,5))) 
      sd(8,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
     &                   m1(2,4),m1(2,4),m1(2,4),0.D0   ))
      sd(9,1) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,2),
     &                   m1(2,4),m1(2,3),m1(2,4),0.D0   ))
      sd(10,1) = dumd_1(m2(1,2),m2(2,1),m2(2,1),
     &                  m2(2,4),mgc2,m2(1,1),sc,m2(2,6),zeta2)
      sd(11,1) = dumd_2(tg,ug,m2(2,4),m2(2,1),mdiff,mdiffc,m2(2,6)) 
      
      sd(4,2) = real(D04(m2(2,1),eps1,eps1,m2(2,1),m2(1,3),m2(1,1),
     &                m1(2,4),m1(2,5),m1(2,5),m1(2,5)))
      sd(5,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
     &                m1(2,4),m1(2,4),m1(2,4),m1(2,5)))
      sd(6,2) = real(D04(eps1,m2(2,1),eps1,m2(2,1),m2(1,2),m2(1,3),
     &                m1(2,5),m1(2,5),m1(2,4),m1(2,4)))
c      sd(7,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
c     &                m1(2,4),m1(2,3),m1(2,4),m1(2,5))) 
      sd(8,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
     &                   m1(2,4),m1(2,3),m1(2,4),0.D0   ))
      sd(9,2) = real(D04(eps1,eps1,m2(2,1),m2(2,1),m2(1,1),m2(1,3),
     &                   m1(2,4),m1(2,4),m1(2,4),0.D0   ))
      sd(10,2) = dumd_1(m2(1,3),m2(2,1),m2(2,1),
     &                  m2(2,4),mgc2,m2(1,1),sc,m2(2,6),zeta2)
      sd(11,2) = sd(11,1)

      return 
      end 


