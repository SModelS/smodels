CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  ROUTINE FILLING THE ARRAY INCLUDING ALL ANGULAR INTEGRALS          C
C   FOR THE STOP PAIR HADROPRODUCTION   [ SOFT ]                      C
C                                                                     C
C                                                                     C
C     SOFT_ARRAY(SOFTIN,SOF1)                                         C
C                                                                     C
C                                                                     C
C  INPUT  :   SOFTIN(1:4)    CONTAINING   S,T1,U1,MS1^2               C
C                                                                     C
C  OUTPUT :   SOF1(1:8)     ALL THE FINITE PARTS OF THE SOFT          C
C                            INTEGRALS ACCORDING TO ROLAND            C
C                                                                     C
C  N.B. : THE INTEGRALS ARE OF THE FORM                               C
C                                                                     C
C         PI * DELTA^EE/MS1^EE * SOF1(1)                              C
C          = INT_0^DELTA     S4^(1-2*EE)/(S4+MS1^2)^(1-EE)            C
C            INT D(OMEGA_N)  1/S4^2                                   C
C                                                                     C
C  CALLED FUNCTIONS : LIR2(X)                                         C
C                                                                     C
C                                                                     C
C  NECESSARY SUBROUTINES : XSPENCE                                    C
C                                                                     C
C                                                                     C
C  LAST MODIFIED : 21.10.96                                           C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c --------------------------------------------------------------------
c --------------------------------------------------------------------
      subroutine SOFT_ARRAY_LQ(m2,sof1)

      implicit none 

      integer n
      real*8 m2(4),sof1(8)
      real*8 zeta2,xw,xx,Lir2,s
      complex*16 CSPEN

      Lir2(s) = real( CSPEN(dcmplx(s)) )     

      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0
      xw = sqrt( 1.D0 - 4.D0*m2(4)/m2(1) )
      xx = ( 1.D0 - xw )/( 1.D0 + xw )
      
c --------------------------
c  1/(s4^2)
      sof1(3) = - 2.D0/m2(4)  

c  1/(s3*s4)      
      sof1(1) = 2.D0/(m2(1)*xw)*(  2.D0*Lir2(xx) + 2.D0*Lir2(-xx) 
     &                           - log(xx)**2 
     &                           + 2.D0*log(xx)*log(1.D0-xx**2)
     &                           - zeta2                        )

c  1/(tp*s4)
      sof1(4) = - 3.D0*zeta2 / (2.D0*m2(3)) 

c  1/(up*s4) 
      sof1(5) = - 3.D0*zeta2 / (2.D0*m2(2)) 

c  1/(s3^2) 
      sof1(2) = 2.D0/m2(4)* 
     &          (m2(1)-2.D0*m2(4))*log(xx)/(m2(1)*xw) 

c  1/(tp*up)
      sof1(8) = 2.D0/m2(1)*( log(m2(2)*m2(3)/m2(1)/m2(4))**2 /2.D0
     &                     + Lir2(1.D0 - m2(1)*m2(4)/m2(2)/m2(3) )
     &                     - 3.D0/2.D0 * zeta2 ) 

c  1/(up*s3)
      sof1(7) = ( - 2.D0*log(xx)*log(m2(3)/m2(2))
     &            - log(xx)**2 
     &            + log(m2(3)/m2(2))**2
     &            + 2.D0*Lir2(1.D0-m2(3)/(m2(2)*xx))
     &            - 2.D0*Lir2(1.D0-m2(2)/(m2(3)*xx))
     &            - 3.D0/2.D0*zeta2                  )/m2(3)

c  1/(tp*s3)
      sof1(6) = ( - 2.D0*log(xx)*log(m2(2)/m2(3))
     &            - log(xx)**2 
     &            + log(m2(2)/m2(3))**2
     &            + 2.D0*Lir2(1.D0-m2(2)/(m2(3)*xx))
     &            - 2.D0*Lir2(1.D0-m2(3)/(m2(2)*xx))
     &            - 3.D0/2.D0*zeta2                  )/m2(2)

c adjust the right prefactor       
      
      do 996 n = 1,8
 996     sof1(n) = sof1(n)/8.D0


      end




















