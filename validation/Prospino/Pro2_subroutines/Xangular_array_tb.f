CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  ROUTINE FILLING THE ARRAY INCLUDING ALL ANGULAR INTEGRALS          C
C   FOR THE STOP PAIR HADROPRODUCTION WITH ALL INCOMING               C
C   STATES, DERIVED FROM XANGULAR_ARRAY.F  [ HARD ]                   C
C      REFERENCE MOSTLY PRD 40 (1989) 54                              C
C                                                                     C
C                                                                     C
C     ANGULAR_ARRAY_TB_QS(HARDIN,ANG4)                                C
C     LOGAS_TB_QS(HARDIN,COLO1)                                       C
C     ANGULAR_ARRAY_TB_GS(HARDIN,ANG4)                                C
C     LOGAS_TB_GS(HARDIN,COLO1)                                       C
C     ANGULAR_ARRAY_TB_CS(HARDIN,ANG4)                                C
C     LOGAS_TB_CS(HARDIN,COLO1)                                       C
C                                                                     C
C                                                                     C
C  INPUT  :   HARDIN(1:5)    CONTAINING   S,T1,U1,S4,MS1^2            C
C                                                                     C
C  OUTPUT :   ANG4(1:71)     ALL THE FINITE PARTS OF THE ANGULAR      C
C                            INTEGRALS ACCORDING TO ROLAND            C
C             COLO1(1:9)     ALL THE TYPICAL LOGARITHMS               C
C                                                                     C
C  N.B. : THE INTEGRALS ARE OF THE FORM 1/S3 = PI*ANG4(37)            C
C         I.E. THE INTEGRATION INCLUDES THE MEASURE                   C
C                     ANG4(37) = INT D OMEGA_N /PI                    C
C                                                                     C
C                                                                     C
C  NEEDED FUNCTIONS :  IP2P0,IP1P0,                                   C
C                      IM1P2,IM1P1,IM2P2,IM2P1,                       C
C                      IBP1P2,IBP1P1,                                 C
C                      IHP1P1                                         C
C  NEEDED FUNCTIONS :  IP2P2,IP1P2,IP2P1,IP2P0,IP1P1,IP1P0,           C
C                      IM1P2,IM1P1,IM2P2,IM2P1,                       C
C                      IBP2P2,IBP2P1,IBP1P2,IBP1P1,IBP2P0,IBP1M1,     C
C                      IHP1P1                                         C
C                                                                     C
C                                                                     C
C                                                                     C
C  LAST MODIFIED : 15.11.96                                           C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_TB_QS(m2,ang4)

      implicit none 

      integer n

      real*8 m2(1:5),ang4(1:71)

      real*8 xd,w1,w2,w3,e1,e2,p,cp,sp
      real*8 aa(1:2,1:3)

      real*8 Ip2p0,Ip1p0
      real*8 Im1p2,Im1p1,Im2p2,Im2p1
      real*8 Ibp1p2,Ibp1p1
      real*8 Ihp1p1

      aa(1,3) = 9.99D-15

c -----------------------------------------
c the unnecessary parts of the ang4 array set to 1.D40
      do n =1,71 
         ang4(n) = 1.D40
      end do

c -----------------------------------------
c define coordinate system A : { k2||z }
      xd = 2.D0 * sqrt( m2(4) + m2(5) ) 

      w1 = ( m2(1) + m2(3) )/xd 
      w2 = ( m2(1) + m2(2) )/xd 
      w3 = m2(4)/xd
      e1 = ( m2(4) + 2.D0*m2(5) )/xd
      e2 = - ( m2(2) + m2(3) + 2.D0*m2(5) )/xd
      p  = sqrt( (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )/xd
      cp = (  m2(2)*m2(4) - m2(1)*( m2(3) + 2.D0*m2(5) ) )
     &        /(m2(1)+m2(2))
     &        /sqrt(  (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c$$$c the integrals containing {u6,s3}
c$$$      aa(1,1) = -2.D0*w2*e1
c$$$      aa(1,2) = -2.D0*w2*w3
c$$$      aa(2,1) = 2.D0*w3*e2
c$$$      aa(2,2) = -2.D0*w3*p*cp
c$$$      aa(2,3) = -2.D0*w3*p*sp
c$$$
c$$$      ang4(1)  = Ip2p2(aa)
c$$$      ang4(2)  = Ip2p1(aa)
c$$$      ang4(3)  = Ip1p2(aa)
c$$$      ang4(4)  = Ip1p1(aa)
c$$$  
c$$$c  the integrals containing {u6,s5}
c$$$      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
c$$$      aa(2,2) = 2.D0*w3*p*cp
c$$$      aa(2,3) = 2.D0*w3*p*sp 
c$$$      
c$$$      ang4(18) = Ip1p1(aa) 
c$$$
c$$$c  the integrals containing {u6,u7}
c$$$      aa(2,1) = -2.D0*e1*w1
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w2 )
c$$$      aa(2,3) = -2.D0*w3*p*sp
c$$$
c$$$      ang4(19) = Ip1p1(aa) 
c$$$      ang4(21) = Im1p1(aa) 
c$$$
c$$$c  the integrals containing {u6} 
c$$$      ang4(40) = Ip2p0(aa) 
c$$$      ang4(41) = Ip1p0(aa) 


c  the integrals containing {tp,up} 
      aa(1,1) = -2.D0*w2*w3
      aa(1,2) = 2.D0*w2*w3
      aa(2,1) = -2.D0*w1*w3
      aa(2,2) = 2.D0*w3*( p*cp - w2 )
      aa(2,3) = 2.D0*w3*p*sp
      
      ang4(48) = Ihp1p1(aa)

c  the integrals containing {tp,s3} 
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(49) = Ibp1p1(aa)
      
c  the integrals containing {tp,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

c$$$      ang4(51) = Ibp2p2(aa)
      ang4(52) = Ibp1p2(aa) 
c$$$      ang4(53) = Ibp2p1(aa) 
      ang4(54) = Ibp1p1(aa) 

c$$$c  the integrals containing {tp,u7}  
c$$$      aa(2,1) = -2.D0*e1*w1
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w2 )
c$$$      aa(2,3) = -2.D0*w3*p*sp
c$$$
c$$$      ang4(59) = Ibp2p2(aa) 
c$$$      ang4(60) = Ibp1p2(aa)
c$$$      ang4(61) = Ibp2p1(aa) 
c$$$      ang4(62) = Ibp1p1(aa) 
c$$$      ang4(63) = Ibp1m1(aa) 
c$$$
c$$$c  the integrals containing {tp}
c$$$
c$$$      ang4(70) = Ibp2p0(aa) 



c ---------------------------
c finite integrals in the coordinate system C : { p2||z }
c      n.b. the definitions of the array cm stay the same  

c  the integrals containing {s3}
      aa(1,1) = 2.D0*w3*e2
      aa(1,2) = -2.D0*w3*p

      ang4(36) = Ip2p0(aa) 
      ang4(37) = Ip1p0(aa)

c  the integrals containing {s5}
      aa(1,1) = 2.D0*( e1*e2 + m2(5) )
      aa(1,2) = 2.D0*w3*p

      ang4(38) = Ip2p0(aa) 
      ang4(39) = Ip1p0(aa)



c ---------------------------
c finite integrals in the coordinate system B : { k1||z }
c slight change in the definition of cm :
      cp  = ( m2(3)*m2(4) - m2(1)*( m2(2) + 2.D0*m2(5) ))
     &        /(m2(1)+m2(3))
     &        /sqrt( (m2(3)+m2(2))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c  the integrals containing {u7,s3}
      aa(1,1) = -2.D0*w1*e1
      aa(1,2) = -2.D0*w1*w3
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(5)  = Ip2p2(aa)
c$$$      ang4(6)  = Ip2p1(aa)
c$$$      ang4(7)  = Ip1p2(aa)
c$$$      ang4(8)  = Ip1p1(aa)
      ang4(9)  = Im2p2(aa)
      ang4(10) = Im1p2(aa)
      ang4(11) = Im2p1(aa)
      ang4(12) = Im1p1(aa)

c  the integrals containing {u7,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 
      
      ang4(13) = Im2p2(aa)
      ang4(14) = Im1p2(aa) 
      ang4(15) = Im2p1(aa) 
      ang4(16) = Im1p1(aa) 
c$$$      ang4(17) = Ip1p1(aa) 
      
c$$$c  the integrals containing {u7,u6}
c$$$      aa(2,1) = -2.D0*e1*w2
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
c$$$      aa(2,3) = -2.D0*w3*p*sp
c$$$
c$$$      ang4(20) = Im1p1(aa) 
c$$$
c$$$c  the integrals containing {u7} 
c$$$      ang4(43) = Ip2p0(aa) 
c$$$      ang4(44) = Ip1p0(aa) 

c  the integrals containing {up,u6} 
      aa(1,1) = -2.D0*w1*w3
      aa(1,2) = 2.D0*w1*w3
c$$$      aa(2,1) = -2.D0*e1*w2
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
c$$$      aa(2,3) = -2.D0*w3*p*sp
c$$$
c$$$      ang4(64) = Ibp2p2(aa) 
c$$$      ang4(65) = Ibp1p2(aa) 
c$$$      ang4(66) = Ibp2p1(aa) 
c$$$      ang4(67) = Ibp1p1(aa) 
c$$$      ang4(68) = Ibp1m1(aa) 
c$$$
c$$$c  the integrals containing {up}
c$$$      ang4(71) = Ibp2p0(aa) 


c  the integrals containing {up,s3} 
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(50) = Ibp1p1(aa) 

c  the integrals containing {up,s5} 
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

c$$$      ang4(55) = Ibp2p2(aa)
      ang4(56) = Ibp1p2(aa) 
c$$$      ang4(57) = Ibp2p1(aa) 
      ang4(58) = Ibp1p1(aa) 

      end 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  NEEDED FUNCTIONS :  IP2P2,IP1P2,IP2P1,IP2P0,IP1P1,IP1P0,           C
C                      IM1P2,IM1P1,IM2P2,IM2P1,                       C
C                      IBP2P2,IBP2P1,IBP1P2,IBP1P1,IBP2P0,IBP1M1,     C
C                      IHP1P1                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_TB_GS(m2,ang4)

      implicit none 

      integer n

      real*8 m2(1:5),ang4(1:71)

      real*8 xd,w1,w2,w3,e1,e2,p,cp,sp
      real*8 aa(1:2,1:3)

      real*8 Ip2p2,Ip1p2,Ip2p1,Ip2p0,Ip1p1,Ip1p0
      real*8 Im1p2,Im1p1,Im2p2,Im2p1
      real*8 Ibp2p2,Ibp2p1,Ibp1p2,Ibp1p1,Ibp2p0,Ibp1m1
      real*8 Ihp1p1

      aa(1,3) = 9.99D-15

c -----------------------------------------
c the unnecessary parts of the ang4 array set to 1.D40
      do n =1,71 
         ang4(n) = 1.D40
      end do

c -----------------------------------------
c define coordinate system A : { k2||z }
      xd = 2.D0 * sqrt( m2(4) + m2(5) ) 

      w1 = ( m2(1) + m2(3) )/xd 
      w2 = ( m2(1) + m2(2) )/xd 
      w3 = m2(4)/xd
      e1 = ( m2(4) + 2.D0*m2(5) )/xd
      e2 = - ( m2(2) + m2(3) + 2.D0*m2(5) )/xd
      p  = sqrt( (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )/xd
      cp = (  m2(2)*m2(4) - m2(1)*( m2(3) + 2.D0*m2(5) ) )
     &        /(m2(1)+m2(2))
     &        /sqrt(  (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c the integrals containing {u6,s3}
      aa(1,1) = -2.D0*w2*e1
      aa(1,2) = -2.D0*w2*w3
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(1)  = Ip2p2(aa)
      ang4(2)  = Ip2p1(aa)
      ang4(3)  = Ip1p2(aa)
      ang4(4)  = Ip1p1(aa)
  
c  the integrals containing {u6,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 
      
      ang4(18) = Ip1p1(aa) 

c  the integrals containing {u6,u7}
      aa(2,1) = -2.D0*e1*w1
      aa(2,2) = -2.D0*w3*( p*cp - w2 )
      aa(2,3) = -2.D0*w3*p*sp

      ang4(19) = Ip1p1(aa) 
      ang4(21) = Im1p1(aa) 

c  the integrals containing {u6} 
      ang4(40) = Ip2p0(aa) 
      ang4(41) = Ip1p0(aa) 


c  the integrals containing {tp,up} 
      aa(1,1) = -2.D0*w2*w3
      aa(1,2) = 2.D0*w2*w3
      aa(2,1) = -2.D0*w1*w3
      aa(2,2) = 2.D0*w3*( p*cp - w2 )
      aa(2,3) = 2.D0*w3*p*sp
      
      ang4(48) = Ihp1p1(aa)

c  the integrals containing {tp,s3} 
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(49) = Ibp1p1(aa)
      
c  the integrals containing {tp,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

      ang4(51) = Ibp2p2(aa)
      ang4(52) = Ibp1p2(aa) 
      ang4(53) = Ibp2p1(aa) 
      ang4(54) = Ibp1p1(aa) 

c  the integrals containing {tp,u7}  
      aa(2,1) = -2.D0*e1*w1
      aa(2,2) = -2.D0*w3*( p*cp - w2 )
      aa(2,3) = -2.D0*w3*p*sp

      ang4(59) = Ibp2p2(aa) 
      ang4(60) = Ibp1p2(aa)
      ang4(61) = Ibp2p1(aa) 
      ang4(62) = Ibp1p1(aa) 
      ang4(63) = Ibp1m1(aa) 

c  the integrals containing {tp}

      ang4(70) = Ibp2p0(aa) 



c ---------------------------
c finite integrals in the coordinate system C : { p2||z }
c      n.b. the definitions of the array cm stay the same  

c  the integrals containing {s3}
      aa(1,1) = 2.D0*w3*e2
      aa(1,2) = -2.D0*w3*p

      ang4(36) = Ip2p0(aa) 
      ang4(37) = Ip1p0(aa)

c  the integrals containing {s5}
      aa(1,1) = 2.D0*( e1*e2 + m2(5) )
      aa(1,2) = 2.D0*w3*p

      ang4(38) = Ip2p0(aa) 
      ang4(39) = Ip1p0(aa)



c ---------------------------
c finite integrals in the coordinate system B : { k1||z }
c slight change in the definition of cm :
      cp  = ( m2(3)*m2(4) - m2(1)*( m2(2) + 2.D0*m2(5) ))
     &        /(m2(1)+m2(3))
     &        /sqrt( (m2(3)+m2(2))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c  the integrals containing {u7,s3}
      aa(1,1) = -2.D0*w1*e1
      aa(1,2) = -2.D0*w1*w3
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(5)  = Ip2p2(aa)
      ang4(6)  = Ip2p1(aa)
      ang4(7)  = Ip1p2(aa)
      ang4(8)  = Ip1p1(aa)
      ang4(9)  = Im2p2(aa)
      ang4(10) = Im1p2(aa)
      ang4(11) = Im2p1(aa)
      ang4(12) = Im1p1(aa)

c  the integrals containing {u7,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 
      
      ang4(13) = Im2p2(aa)
      ang4(14) = Im1p2(aa) 
      ang4(15) = Im2p1(aa) 
      ang4(16) = Im1p1(aa) 
      ang4(17) = Ip1p1(aa) 
      
c  the integrals containing {u7,u6}
      aa(2,1) = -2.D0*e1*w2
      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
      aa(2,3) = -2.D0*w3*p*sp

      ang4(20) = Im1p1(aa) 

c  the integrals containing {u7} 
      ang4(43) = Ip2p0(aa) 
      ang4(44) = Ip1p0(aa) 

c  the integrals containing {up,u6} 
      aa(1,1) = -2.D0*w1*w3
      aa(1,2) = 2.D0*w1*w3
      aa(2,1) = -2.D0*e1*w2
      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
      aa(2,3) = -2.D0*w3*p*sp

      ang4(64) = Ibp2p2(aa) 
      ang4(65) = Ibp1p2(aa) 
      ang4(66) = Ibp2p1(aa) 
      ang4(67) = Ibp1p1(aa) 
      ang4(68) = Ibp1m1(aa) 

c  the integrals containing {up}
      ang4(71) = Ibp2p0(aa) 


c  the integrals containing {up,s3} 
      aa(2,1) = 2.D0*w3*e2
      aa(2,2) = -2.D0*w3*p*cp
      aa(2,3) = -2.D0*w3*p*sp

      ang4(50) = Ibp1p1(aa) 

c  the integrals containing {up,s5} 
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

      ang4(55) = Ibp2p2(aa)
      ang4(56) = Ibp1p2(aa) 
      ang4(57) = Ibp2p1(aa) 
      ang4(58) = Ibp1p1(aa) 

      end 



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C                                                                     C
C  NEEDED FUNCTIONS :  IP2P0,IP1P1,IP1P0,                             C
C                      IM1P2,                                         C
C                      IBP2P2,IBP2P1,IBP1P2,IBP1P1,IBP2P0,IBP1M1,     C
C                      IHP1P1                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine ANGULAR_ARRAY_TB_CS(m2,ang4)

      implicit none 

      integer n

      real*8 m2(1:5),ang4(1:71)

      real*8 xd,w1,w2,w3,e1,e2,p,cp,sp
      real*8 aa(1:2,1:3)

      real*8 Ip2p0,Ip1p1,Ip1p0
      real*8 Im1p2
      real*8 Ibp2p2,Ibp2p1,Ibp1p2,Ibp1p1,Ibp2p0,Ibp1m1
      real*8 Ihp1p1
c      real*8 Ip2p0,Ip1p1,Ip1p0
c      real*8 Im1p2
c      real*8 Ibp2p2,Ibp2p1,Ibp1p2,Ibp1p1,Ibp2p0,Ibp1m1
c      real*8 Ihp1p1

      aa(1,3) = 9.99D-15

c -----------------------------------------
c the unnecessary parts of the ang4 array set to 1.D40
      do n =1,71 
         ang4(n) = 1.D40
      end do

c -----------------------------------------
c define coordinate system A : { k2||z }
      xd = 2.D0 * sqrt( m2(4) + m2(5) ) 

      w1 = ( m2(1) + m2(3) )/xd 
      w2 = ( m2(1) + m2(2) )/xd 
      w3 = m2(4)/xd
      e1 = ( m2(4) + 2.D0*m2(5) )/xd
      e2 = - ( m2(2) + m2(3) + 2.D0*m2(5) )/xd
      p  = sqrt( (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )/xd
      cp = (  m2(2)*m2(4) - m2(1)*( m2(3) + 2.D0*m2(5) ) )
     &        /(m2(1)+m2(2))
     &        /sqrt(  (m2(2)+m2(3))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c the integrals containing {u6,s3}
      aa(1,1) = -2.D0*w2*e1
      aa(1,2) = -2.D0*w2*w3
c$$$      aa(2,1) = 2.D0*w3*e2
c$$$      aa(2,2) = -2.D0*w3*p*cp
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(1)  = Ip2p2(aa)
c$$$      ang4(2)  = Ip2p1(aa)
c$$$      ang4(3)  = Ip1p2(aa)
c$$$      ang4(4)  = Ip1p1(aa)
  
c  the integrals containing {u6,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 
      
      ang4(18) = Ip1p1(aa) 

c  the integrals containing {u6,u7}
c$$$      aa(2,1) = -2.D0*e1*w1
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w2 )
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(19) = Ip1p1(aa) 
c$$$      ang4(21) = Im1p1(aa) 

c  the integrals containing {u6} 
c$$$      ang4(40) = Ip2p0(aa) 
      ang4(41) = Ip1p0(aa) 


c  the integrals containing {tp,up} 
      aa(1,1) = -2.D0*w2*w3
      aa(1,2) = 2.D0*w2*w3
      aa(2,1) = -2.D0*w1*w3
      aa(2,2) = 2.D0*w3*( p*cp - w2 )
      aa(2,3) = 2.D0*w3*p*sp
      
      ang4(48) = Ihp1p1(aa)

c  the integrals containing {tp,s3} 
c$$$      aa(2,1) = 2.D0*w3*e2
c$$$      aa(2,2) = -2.D0*w3*p*cp
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(49) = Ibp1p1(aa)
      
c  the integrals containing {tp,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

      ang4(51) = Ibp2p2(aa)
      ang4(52) = Ibp1p2(aa) 
      ang4(53) = Ibp2p1(aa) 
      ang4(54) = Ibp1p1(aa) 

c  the integrals containing {tp,u7}  
      aa(2,1) = -2.D0*e1*w1
      aa(2,2) = -2.D0*w3*( p*cp - w2 )
      aa(2,3) = -2.D0*w3*p*sp

      ang4(59) = Ibp2p2(aa) 
      ang4(60) = Ibp1p2(aa)
      ang4(61) = Ibp2p1(aa) 
      ang4(62) = Ibp1p1(aa) 
      ang4(63) = Ibp1m1(aa) 

c  the integrals containing {tp}

      ang4(70) = Ibp2p0(aa) 



c ---------------------------
c finite integrals in the coordinate system C : { p2||z }
c      n.b. the definitions of the array cm stay the same  

c  the integrals containing {s3}
c$$$      aa(1,1) = 2.D0*w3*e2
c$$$      aa(1,2) = -2.D0*w3*p

c$$$      ang4(36) = Ip2p0(aa) 
c$$$      ang4(37) = Ip1p0(aa)

c  the integrals containing {s5}
      aa(1,1) = 2.D0*( e1*e2 + m2(5) )
      aa(1,2) = 2.D0*w3*p

      ang4(38) = Ip2p0(aa) 
      ang4(39) = Ip1p0(aa)



c ---------------------------
c finite integrals in the coordinate system B : { k1||z }
c slight change in the definition of cm :
      cp  = ( m2(3)*m2(4) - m2(1)*( m2(2) + 2.D0*m2(5) ))
     &        /(m2(1)+m2(3))
     &        /sqrt( (m2(3)+m2(2))**2 - 4.D0*m2(1)*m2(5) )
      sp = sqrt(1.D0 - cp**2)

c  the integrals containing {u7,s3}
      aa(1,1) = -2.D0*w1*e1
      aa(1,2) = -2.D0*w1*w3
c$$$      aa(2,1) = 2.D0*w3*e2
c$$$      aa(2,2) = -2.D0*w3*p*cp
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(5)  = Ip2p2(aa)
c$$$      ang4(6)  = Ip2p1(aa)
c$$$      ang4(7)  = Ip1p2(aa)
c$$$      ang4(8)  = Ip1p1(aa)
c$$$      ang4(9)  = Im2p2(aa)
c$$$      ang4(10) = Im1p2(aa)
c$$$      ang4(11) = Im2p1(aa)
c$$$      ang4(12) = Im1p1(aa)

c  the integrals containing {u7,s5}
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 
      
c$$$      ang4(13) = Im2p2(aa)
      ang4(14) = Im1p2(aa) 
c$$$      ang4(15) = Im2p1(aa) 
c$$$      ang4(16) = Im1p1(aa) 
      ang4(17) = Ip1p1(aa) 
      
c  the integrals containing {u7,u6}
c$$$      aa(2,1) = -2.D0*e1*w2
c$$$      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(20) = Im1p1(aa) 

c  the integrals containing {u7} 
c$$$      ang4(43) = Ip2p0(aa) 
      ang4(44) = Ip1p0(aa) 

c  the integrals containing {up,u6} 
      aa(1,1) = -2.D0*w1*w3
      aa(1,2) = 2.D0*w1*w3
      aa(2,1) = -2.D0*e1*w2
      aa(2,2) = -2.D0*w3*( p*cp - w1 ) 
      aa(2,3) = -2.D0*w3*p*sp

      ang4(64) = Ibp2p2(aa) 
      ang4(65) = Ibp1p2(aa) 
      ang4(66) = Ibp2p1(aa) 
      ang4(67) = Ibp1p1(aa) 
c$$$      ang4(68) = Ibp1m1(aa) 

c  the integrals containing {up}
      ang4(71) = Ibp2p0(aa) 


c  the integrals containing {up,s3} 
c$$$      aa(2,1) = 2.D0*w3*e2
c$$$      aa(2,2) = -2.D0*w3*p*cp
c$$$      aa(2,3) = -2.D0*w3*p*sp

c$$$      ang4(50) = Ibp1p1(aa) 

c  the integrals containing {up,s5} 
      aa(2,1) = 2.D0*( e1*e2 + m2(5) )
      aa(2,2) = 2.D0*w3*p*cp
      aa(2,3) = 2.D0*w3*p*sp 

      ang4(55) = Ibp2p2(aa)
      ang4(56) = Ibp1p2(aa) 
      ang4(57) = Ibp2p1(aa) 
      ang4(58) = Ibp1p1(aa) 

      end 

















































