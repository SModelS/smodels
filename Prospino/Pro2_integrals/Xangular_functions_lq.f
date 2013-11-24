CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  ROUTINE INCLUDING THE HARD INTEGRALS FROM THE BIBLE,               C
C   NECESSARY FOR THE ANGULAR INTEGRAL ROUTINES : XANGULAR_ARRAY*     C 
C                                                                     C
C                                                                     C
C  SUBROUTINE :  LOGAS(HARDIN,COLO1)                                  C
C                                                                     C
C  INPUT  :   HARDIN(1:5)    CONTAINING   S,T1,U1,S4,MS1^2            C
C                                                                     C
C  OUTPUT :   COLO1(1:9)     ALL THE TYPICAL LOGARITHMS               C
C                                                                     C
C  FUNCTIONS :  IP2P2,IP1P2,IP2P1,IP2P0,IP1P1,IP1P0,                  C
C               IM1P2,IM1P1,IM2P2,IM2P1,                              C
C               IBP2P2,IBP2P1,IBP1P2,IBP1P1,IBP2P0,IBP1M1,            C
C               IHP1P1                                                C
C                                                                     C
C  ALL THE INTEGRALS ACCORDING TO THE BIBLE,                          C 
C   FINITE PART DIVIDED BY A FACTOR PI                                C
C                                                                     C
C  LAST MODIFIED : 15.11.96                                           C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine LOGAS_LQ(m2,colo1)

      implicit none 

      real*8 m2(1:5),colo1(1:9)
      real*8 xw,xs 
      
      xw = sqrt(1.D0-4.D0*m2(5)/m2(1))
      xs = ( 1.D0-xw )/( 1.D0+xw )

      colo1(1) = log(xs)
      colo1(2) = 1.D40
      colo1(3) = log(m2(1)/m2(5))
      colo1(4) = log(-m2(2)/m2(5))
      colo1(5) = 1.D40
      colo1(6) = 1.D40
      colo1(7) = 1.D40
      colo1(8) = log(-m2(3)/m2(5))
      colo1(9) = log(m2(4)**2/(m2(5)+m2(4))/m2(5))

      end


c ---------------------------------
      real*8 function Ip2p2(aa)

      implicit none 
      real*8 aa(1:2,1:3),xx 

      xx = -2.D0*aa(1,1)*aa(1,2)*aa(2,1)*aa(2,2)
     &   + aa(1,1)**2 * ( aa(2,2)**2 + aa(2,3)**2 )
     &   + aa(1,2)**2 * ( aa(2,1)**2 - aa(2,3)**2 ) 
c$$$      xx =  ( aa(1,1)*aa(2,1) - aa(1,2)*aa(2,2) )**2 
c$$$     &    - ( aa(2,1)**2 - aa(2,2)**2 - aa(2,3)**2 )
c$$$     &      * ( aa(1,1)**2 - aa(1,2)**2 )

      Ip2p2 = 2.D0*aa(1,2)**2 / ( xx * (aa(1,1)**2-aa(1,2)**2) )
     &      + 2.D0*( aa(2,2)**2 + aa(2,3)**2 ) /
     &             (( aa(2,1)**2 - aa(2,2)**2 - aa(2,3)**2 ) * xx )
     &      - 6.D0*aa(1,2)**2*aa(2,3)**2 / xx**2 
     &      + (   aa(1,2)*aa(2,2) / xx**(3.D0/2.D0)
     &          + 3.D0*aa(1,2)*(aa(1,2)*aa(2,1)-aa(1,1)*aa(2,2))
     &                *(  aa(1,1)*( aa(2,2)**2 + aa(2,3)**2 )
     &                  - aa(1,2)*aa(2,1)*aa(2,2)            ) 
     &            / xx**(5.D0/2.D0)
     &         )* log( (aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)+sqrt(xx))
     &                 /(aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)-sqrt(xx)) )

      end


c ---------------------------------
      real*8 function Ip1p2(aa) 

      implicit none 
      real*8 aa(1:2,1:3),xx

      xx =  ( aa(1,1)*aa(2,1) - aa(1,2)*aa(2,2) )**2 
     &    - ( aa(2,1)**2 - aa(2,2)**2 - aa(2,3)**2 )
     &      * ( aa(1,1)**2 - aa(1,2)**2 )

      Ip1p2 =  (  2.D0*aa(1,1)*( aa(2,2)**2+aa(2,3)**2 ) 
     &          - 2.D0*aa(1,2)*aa(2,1)*aa(2,2)           )
     &          / ((aa(2,1)**2-aa(2,2)**2-aa(2,3)**2) * xx)
     &       + aa(1,2) * ( aa(1,2)*aa(2,1)-aa(1,1)*aa(2,2) )
     &         / xx**(3.D0/2.D0)
     &         * log( (aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)+sqrt(xx))
     &                /(aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)-sqrt(xx)) ) 

      end


c ---------------------------------
      real*8 function Ip2p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3),xx 

      xx =  ( aa(1,1)*aa(2,1) - aa(1,2)*aa(2,2) )**2 
     &    - ( aa(2,1)**2 - aa(2,2)**2 - aa(2,3)**2 )
     &      * ( aa(1,1)**2 - aa(1,2)**2 )

      Ip2p1 =  2.D0*aa(1,2)*( aa(1,2)*aa(2,1)-aa(1,1)*aa(2,2) )
     &          /( (aa(1,1)**2-aa(1,2)**2) * xx )
     &       + (  aa(1,1)*( aa(2,2)**2 + aa(2,3)**2 )
     &          - aa(1,2)*aa(2,1)*aa(2,2)           ) / xx**(3.D0/2.D0)
     &         * log( (aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)+sqrt(xx))
     &                /(aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)-sqrt(xx)) )

      end


c ---------------------------------
      real*8 function Ip2p0(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Ip2p0 = 2.D0 / ( aa(1,1)**2 - aa(1,2)**2 )  

      end


c ---------------------------------
      real*8 function Ip1p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3),xx 
      
      xx = -2.D0*aa(1,1)*aa(1,2)*aa(2,1)*aa(2,2)
     &   + aa(1,1)**2 * ( aa(2,2)**2 + aa(2,3)**2 )
     &   + aa(1,2)**2 * ( aa(2,1)**2 - aa(2,3)**2 ) 
c$$$      xx =  ( aa(1,1)*aa(2,1) - aa(1,2)*aa(2,2) )**2 
c$$$     &    - ( aa(2,1)**2 - aa(2,2)**2 - aa(2,3)**2 )
c$$$     &      * ( aa(1,1)**2 - aa(1,2)**2 )

      Ip1p1 = log( (aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)+sqrt(xx))
     &             /(aa(1,1)*aa(2,1)-aa(1,2)*aa(2,2)-sqrt(xx)) )
     &        / sqrt(xx) 

      end


c ---------------------------------
      real*8 function Ip1p0(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Ip1p0 = log((aa(1,1)+aa(1,2))/(aa(1,1)-aa(1,2)))/aa(1,2) 

      end


c ---------------------------------
      real*8 function Im1p2(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Im1p2 = 2.D0* (  aa(1,1)*(aa(2,2)**2+aa(2,3)**2)
     &               - aa(1,2)*aa(2,1)*aa(2,2)        )
     &         / ( ( aa(2,2)**2+aa(2,3)**2 )
     &            *( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 ) )
     &      + aa(1,2)*aa(2,2) / ( aa(2,2)**2+aa(2,3)**2 )**(3.D0/2.D0)
     &       * log( (aa(2,1) + sqrt(aa(2,2)**2+aa(2,3)**2))
     &              /(aa(2,1) - sqrt(aa(2,2)**2+aa(2,3)**2)) )

      end


c ---------------------------------
      real*8 function Im1p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Im1p1 =  2.D0*aa(1,2)*aa(2,2) / (aa(2,2)**2+aa(2,3)**2 )
     &       + (  aa(1,1) * (aa(2,2)**2+aa(2,3)**2)
     &          - aa(1,2)*aa(2,1)*aa(2,2)          )
     &         / (aa(2,2)**2+aa(2,3)**2)**(3.D0/2.D0)
     &         * log( (aa(2,1) + sqrt(aa(2,2)**2+aa(2,3)**2))
     &                /(aa(2,1) - sqrt(aa(2,2)**2+aa(2,3)**2)) )

      end


c ---------------------------------
      real*8 function Im2p2(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Im2p2 = 2.D0*aa(1,2)**2*( aa(2,2)**2-aa(2,3)**2 )
     &            /( aa(2,2)**2 + aa(2,3)**2 )**2
     &      + 2.D0*( aa(1,1)*( aa(2,2)**2+aa(2,3)**2 )
     &                - aa(1,2)*aa(2,1)*aa(2,2)          )**2 
     &               / ( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 )
     &               / ( aa(2,2)**2 + aa(2,3)**2 )**2 
     &      + ( 2.D0*aa(1,2)*aa(2,2)*( aa(1,1)*( aa(2,2)**2+aa(2,3)**2)
     &                                 -aa(1,2)*aa(2,1)*aa(2,2)  )
     &                      /( aa(2,2)**2+aa(2,3)**2 )**(5.D0/2.D0)
     &        + aa(1,2)**2*aa(2,1)*aa(2,3)**2
     &                 /( aa(2,2)**2 + aa(2,3)**2 )**(5.D0/2.D0)    )
     &       * log( (aa(2,1) + sqrt(aa(2,2)**2+aa(2,3)**2))
     &              /(aa(2,1) - sqrt(aa(2,2)**2+aa(2,3)**2)) )

      end


c ---------------------------------
      real*8 function Im2p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 
      
      Im2p1 =  4.D0*aa(1,1)*aa(1,2)*aa(2,2)/( aa(2,2)**2+aa(2,3)**2 )
     &       + aa(1,2)**2*aa(2,1)*( aa(2,3)**2 - 2.D0*aa(2,2)**2 )
     &                           /( aa(2,2)**2+aa(2,3)**2 )**2 
     &       + (  (  aa(1,1)*(  aa(2,2)**2+aa(2,3)**2 )
     &             - aa(1,2)*aa(2,1)*aa(2,2)            )**2
     &             / ( aa(2,2)**2+aa(2,3)**2 )**(5.D0/2.D0)
     &          - aa(1,2)**2*aa(2,3)**2
     &            * ( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 )
     &            / ( 2.D0*( aa(2,2)**2+aa(2,3)**2 )**(5.D0/2.D0) ) )
     &           * log( (aa(2,1) + sqrt(aa(2,2)**2+aa(2,3)**2))
     &                  /(aa(2,1) - sqrt(aa(2,2)**2+aa(2,3)**2)) ) 

      end


c ---------------------------------
      real*8 function Ibp2p2(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp2p2 = ( (  3.D0*aa(2,3)**2 / ( aa(2,1)+aa(2,2) )**2
     &            + 2.D0*aa(2,2) / ( aa(2,1)+aa(2,2) )      )
     &           * log( ( aa(2,1)+aa(2,2)) **2
     &                /( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 ) )
     &          - 8.D0*aa(2,3)**2/( aa(2,1)+aa(2,2) )**2 
     &          + 2.D0*( aa(2,2)**2+aa(2,3)**2 )
     &                /( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 )
     &          - 1.D0                                         )
     &         /( aa(1,1)**2*( aa(2,1)+aa(2,2) )**2 )

      end


c ---------------------------------
      real*8 function Ibp2p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp2p1 = (( aa(2,2)**2+aa(2,1)*aa(2,2)+aa(2,3)**2 )
     &            /( aa(2,1)+aa(2,2) )**2 
     &            * log( (aa(2,1)+aa(2,2))**2 
     &                   /(aa(2,1)**2-aa(2,2)**2-aa(2,3)**2) )
     &          - 2.D0*aa(2,3)**2/(aa(2,1)+aa(2,2))**2
     &          - 1.D0                                         )
     &         /( aa(1,1)**2*( aa(2,1)+aa(2,2) ) )

      end


c ---------------------------------
      real*8 function Ibp1p2(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp1p2 = (  log( (aa(2,1)+aa(2,2))**2 
     &                 /(aa(2,1)**2-aa(2,2)**2-aa(2,3)**2) )
     &          + 2.D0*( aa(2,2)**2+aa(2,3)**2+aa(2,1)*aa(2,2) )
     &                /( aa(2,1)**2-aa(2,2)**2-aa(2,3)**2 )     )
     &         /aa(1,1) / ( aa(2,1)+aa(2,2) )**2 

      end


c ---------------------------------
      real*8 function Ibp1p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp1p1 = log( (aa(2,1)+aa(2,2) )**2 
     &              /(aa(2,1)**2-aa(2,2)**2-aa(2,3)**2) )
     &         /( aa(1,1)*( aa(2,1)+aa(2,2) ) )

      end


c ---------------------------------
      real*8 function Ibp2p0(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp2p0 = - 1.D0/aa(1,1)**2 

      end


c ---------------------------------
      real*8 function Ibp1m1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ibp1m1 = -2.D0*aa(2,2)/aa(1,1)

      end


c ---------------------------------
      real*8 function Ihp1p1(aa) 
      
      implicit none 
      real*8 aa(1:2,1:3) 

      Ihp1p1 = 2.D0/(aa(1,1)*aa(2,1))
     &             *aa(2,1)/( aa(2,1)+aa(2,2) )
     &             *log( (aa(2,1)+aa(2,2))/(2.D0*aa(2,1)) )

      end





















