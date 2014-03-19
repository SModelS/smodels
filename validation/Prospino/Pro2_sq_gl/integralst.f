      REAL*8 FUNCTION SL3B0A(MST12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MST12,M12,M22
      INTEGER INUM

      IF (INUM.EQ.1) 
     +     SL3B0A = -LOG(M12/MST12)
C***  FINITE PART OF B0(K1,M1,M1)

      IF (INUM.EQ.2) 
     +     SL3B0A = 1.D0 -M12/(M12 -M22)*LOG(M12/M22) -Log(m22/mst12) 
C***  FINITE PART OF B0(K1,M1,M2)

      RETURN
      END


      REAL*8 FUNCTION SL3B0B(S,MST12,M12,INUM)
      IMPLICIT NONE
      REAL*8 S,MST12,M12
      COMPLEX*16 KBETA,KX
      INTEGER INUM

      IF (INUM.EQ.1) 
     +     SL3B0B = 2.D0 - LOG(S/MST12)
C***  FINITE PART OF B0(K1+K2,0,0)

      IF (INUM.EQ.2) THEN
         KBETA  = SQRT( DCMPLX(1.D0 - 4.D0*M12/S)) 
         KX     = (DCMPLX(1.D0) -KBETA)*(DCMPLX(1.D0)+KBETA)**(-1)
         SL3B0B = 2.D0 -LOG(M12/MST12) +REAL(KBETA*LOG(-KX))
C***  FINITE PART OF B0(K1+K2,M1,M1)
      END IF

      RETURN
      END


      REAL*8 FUNCTION SL3B0C(MST12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MSt12,M12,M22
      COMPLEX*16 KDEL,KA,KB,KX1,KX2
      INTEGER INUM
      
      IF (INUM.EQ.1) 
     +     SL3B0C = 2.D0
C***  FINITE PART OF B0(-P2,0,MST1)

      IF (INUM.EQ.2) THEN
         KDEL = (0.D0,1.D-15)
         KA = dcmplx(1.D0 +(M12 -M22)/Mst12)
         KB = dcmplx(M12/Mst12) -KDEL
         KX1 = KA/2.D0 + SQRT(KA**2/4.D0 - KB)
         KX2 = KA/2.D0 - SQRT(KA**2/4.D0 - KB)
         SL3B0C = 2.D0 + real( 
     +        - (1.D0 -KX1)*LOG(1.D0-KX1) -KX1*LOG(-KX1) 
     +        - (1.D0 -KX2)*LOG(1.D0-KX2) -KX2*LOG(-KX2) )
C***  FINITE PART OF B0(-P2,M1,m2)
      END IF

      RETURN
      END


      REAL*8 FUNCTION SL3B0D(T,MSt12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MSt12,M12,M22,T,A,B,X1,X2
      INTEGER INUM

      IF (INUM.EQ.1) 
     +     SL3B0D = 2.D0 -(T-MSt12)/T*LOG((MSt12-T)/MSt12)
C***  FINITE PART OF B0(K2-P2,0,MST1)

      IF (INUM.EQ.2) THEN
         A = 1.D0 +(M12 -M22)/T
         B = M12/T
         X1 = A/2.D0 + SQRT(A**2/4.D0 - B)
         X2 = A/2.D0 - SQRT(A**2/4.D0 - B)
         SL3B0D = 2.D0 -LOG(-T/MST12) 
     +        - (1.D0 -X1)*LOG(ABS(1.D0-X1)) -X1*LOG(ABS(-X1)) 
     +        - (1.D0 -X2)*LOG(ABS(1.D0-X2)) -X2*LOG(ABS(-X2))
C***  FINITE PART OF B0(K2-P2,M1,M2)
      END IF

      RETURN
      END



      REAL*8 FUNCTION SL3BP(MSt12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MST12,M12,M22,M2
      COMPLEX*16 KA,KB,KX1,KX2,KDEL
      INTEGER INUM

      IF (INUM.EQ.1) 
     +     SL3BP = -1.D0/MST12
C***  FINITE PART OF B0P(-P2,0,MST1)

      IF (INUM.EQ.2) THEN
         KDEL = (0.D0,1.D-15)
         KA = dcmplx(1.D0 +(M12 -M22)/MST12)
         KB = dcmplx(M12/MST12) -KDEL
         KX1 = KA/2.D0 + SQRT(KA**2/4.D0 - KB)
         KX2 = KA/2.D0 - SQRT(KA**2/4.D0 - KB)
         SL3BP    = (-1D0 + real(
     +        + (kx1 -kx1**2)/(kx1-kx2)*log((kx1-1.D0)/kx1)
     +        - (kx2 -kx2**2)/(kx1-kx2)*log((kx2-1.D0)/kx2) ))/MST12
C***  FINITE PART OF B0P(-P2,M1,M2) 
      END IF

      IF (INUM.EQ.3) THEN
         M2 = M22 -M12         
         SL3BP  = -1.D0/2.D0/M2 +M22/M2**2 - M22*M12/M2**3*LOG(M22/M12)
C***  FINITE PART OF B0P(K1,M1,M2) 
      END IF

      RETURN
      END




      REAL*8 FUNCTION SL3C0A(S,MST12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MST12,M12,M22,S,ZETA2,pi,SPENCE
      COMPLEX*16 KX,KBEta,KSPENC
      INTEGER INUM

      IF (INUM.EQ.1) THEN
         pi = 4D0 * atan(1D0)
         ZETA2 = pi**2/6D0
         SL3C0A = 1.D0/S*( 0.5D0*LOG(S/MST12)**2 - 3.5D0*ZETA2)
      END IF
C***  FINITE PART OF C0(K1,K2,0,0,0)

      IF (INUM.EQ.2) THEN
         KBETA  = SQRT( DCMPLX(1.D0 - 4.D0*M12/S)) 
         KX     = (DCMPLX(1.D0) -KBETA)*(DCMPLX(1.D0)+KBETA)**(-1)
         SL3C0A = 1.D0/S * REAL( LOG(-KX)**2 )/2.D0
C***  FINITE PART OF C0(K1,K2,M1,M1,M1)
      END IF

      IF (INUM.EQ.3) THEN
         KBETA  = SQRT( DCMPLX(1.D0 - 4.D0*M12/S)) 
         KX     = (DCMPLX(1.D0) -KBETA)*(DCMPLX(1.D0)+KBETA)**(-1)
         SL3C0A = 1.D0/S * REAL( 
     +        + LOG(-KX)**2 +KSPENC(1.D0 +M22/M12*KX) 
     +        + KSPENC(1.D0 +M22/M12/KX) 
     +        -2.D0*SPENCE(1.D0 -M22/M12) )
C***  FINITE PART OF C0(K1,K2,M1,M2,M1)
      end if

      RETURN
      END


      REAL*8 FUNCTION SL3C0B(S,MSt12,M12,M22,INUM)
      IMPLICIT NONE
      REAL*8 MST12,M12,M22,S,SB,pi,ZETA2,SPENCE,C0GENS,beta,xs
      INTEGER INUM
      pi = 4D0 * atan(1D0)
      ZETA2 = pi**2/6D0
      beta = sqrt(1d0 - 4*mst12/s)
      sb = s *beta
      Xs = (1d0 -beta)/(1d0 +beta)

      IF (INUM.EQ.1) 
     +     SL3C0B = 1.D0/SB *(
     +     - 2.D0*SPENCE(XS) - 2.D0*LOG(XS)*LOG(1-XS)
     +     + LOG(XS)**2/2.D0 - 4.D0*ZETA2 )
C***  FINITE PART OF C0(P1,-P2,MST1,0,MST1)

      IF (INUM.EQ.2) 
     +     SL3C0B= 1.D0/SB*(2.D0*SPENCE(-XS) + 0.5D0*LOG(XS)**2 +ZETA2)
C***  FINITE PART OF C0(-P1,-P2,0,MST1,0)

      IF (INUM.EQ.3) 
     +     SL3C0B = C0GENS(S,MST12,M12,M22)
C***  FINITE PART OF C0(-P1,-P2,M2,M1,M2)

      RETURN
      END


      REAL*8 FUNCTION SL3C0C(T,MST12,M12,M22,M32,INUM)
      IMPLICIT NONE
      REAL*8 MST12,M12,M22,M32,T,T1,SPENCE,Pi,ZETA2,C0FIN,DEL
      INTEGER INUM
      pi = 4D0 * atan(1D0)
      ZETA2 = pi**2/6D0
      T1 = T -MST12
      DEL = 1.D-4

      IF (INUM.EQ.1) 
     +     SL3C0C = 1.D0/T1*
     +     (+ 0.25D0*ZETA2 + LOG(-T1/MST12)**2 +SPENCE(T/MST12))
C***  FINITE PART OF C0(K2,-P2,0,0,MST1)

      IF (INUM.EQ.2)
     +     SL3C0C = C0FIN(DEL,MST12,T,M12,M22,M32) 
C***  FINITE PART OF C0(K2,-P2,M1,M2,m3)

      IF (INUM.EQ.3) 
     +     SL3C0C = 1.D0/T1*(ZETA2 - SPENCE(T/MST12))
C***  FINITE PART OF C0(K2,-P2,MST1,MST1,0)

      RETURN
      END



      REAL*8 FUNCTION SL3D0(S,T,U,MSt12,M12,M22,M32,INUM)
      IMPLICIT NONE
      REAL*8 MST12,S,T,U,SB,beta,XS,SPENCE,ZETA2,pi,T1,U1
      REAL*8 M12,M22,M32,M1,M2,M3,Rd0REG
      INTEGER INUM
      pi = 4D0 * atan(1D0)
      ZETA2 = pi**2/6D0
      T1 = T -MST12
      U1 = U -MST12
      beta = sqrt(1d0 - 4*mst12/s)
      sb = s *beta
      Xs = (1d0 -beta)/(1d0 +beta)

      IF (INUM.EQ.1) 
     +     SL3D0 = 
     +     1.D0/S/T1*(2.D0*LOG(S/MST12)*LOG(-T1/MST12) -4.D0*ZETA2)
C***  FINITE PART OF D0(K1,K2,-P2,0,0,0,MST1)

      IF (INUM.EQ.2) 
     +     SL3D0 = 1.D0/SB/T1 *(
     +     -2.D0*LOG(XS)*LOG(1.D0-XS) +2.D0*LOG(XS)*LOG(1.D0+XS)
     +     -2.D0*LOG(XS)*LOG(-T1/MST12) -2.D0*SPENCE(XS)
     +     +2.D0*SPENCE(-XS) -3.D0*ZETA2 )
C***  FINITE PART OF D0(K1,K2,-P2,MSt1,MSt1,MSt1,0)

      IF (INUM.EQ.3) 
     +     SL3D0 = 1.D0/T1/U1 *(
     +     2.D0*LOG(-T1/MST12)*LOG(-U1/MST12) -7.D0/2.D0*ZETA2 )
C***  FINITE PART OF D0(K1,-P2,K2,MST1,MST1,0,0)

      IF (INUM.eq.4) THEN
         M1 = SQRT(M12)
         M2 = SQRT(M22)
         SL3D0 = RD0REG(MST12,0.D0,0.D0,MST12,S,T,M1,M2,M2,M2)
C***  FINITE PART OF D0(K1,K2,-P2,M2,M2,M2,M1)
      END IF

      IF (INUM.eq.5) THEN
         M1 = SQRT(M12)
         M2 = SQRT(M22)
         SL3D0 =  RD0REG(0.D0,MST12,0.D0,MST12,U,T,M1,M1,M2,M2)
C***  FINITE PART OF D0(K1,-P2,K2,M1,M1,M2,M2)
      END IF

      IF (INUM.eq.6) THEN
         M1 = SQRT(M12)
         M2 = SQRT(M22)
         M3 = SQRT(M32)
         SL3D0 = RD0REG(MST12,0.D0,0.D0,MST12,S,T,M1,M2,M3,M2)
C***  FINITE PART OF D0(K1,K2,-P2,M2,M3,M2,M1)
      END IF


      return
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     SUBROUTINE CALCULATING THE GENERAL MASSIVE THREE POINT FUNCTION C
C                                                                     C
C           Cc(p0**2,p1**2,p2**2,m0**2,m1**2,m2**2)                   C
C                                                                     C
C      ONLY THE REAL PART :                                           C
C                                                                     C
C           C(p0**2,p1**2,p2**2,m0**2,m1**2,m2**2)                    C

c ---------------------------------------------------------------------


      real*8 function C0fin(p02,p12,p22,m02,m12,m22)
      implicit real*8 (a-h,m,o-y), complex*16 (z) 
      complex*16 Cc 
c$$$  attention to different notation 

      C0fin = real(Cc(p02,p12,p22,m12,m22,m02))
      return
      end



c ---------------------------------------------------------------------
      complex*16 function Cc(p02,p12,p22,m02,m12,m22)
      implicit real*8 (a-h,m,o-y), complex*16 (z) 
      complex*16 Kspenc,detafct
      dimension m2(0:2),p2(0:2,0:2),zy0(0:2),zalpha(0:2)
      dimension zx(0:2,0:1),zy(0:2,0:1)

      zi = dcmplx(0.D0,1.D0)
      pi = 4D0*atan(1.D0)
      eps=1.D-14

* remember the masses in denner's article :
*     m(0) between p1 and p2 -> m1
*     m(1) between p1 and p0 -> m0
*     m(2) between p2 and p0 -> m2

      m2(0) = m12 
      m2(1) = m02 
      m2(2) = m22 

      p2(1,0) = p12 
      p2(2,0) = p22 
      p2(2,1) = p02 

      do 77, i=0,2
         p2(i,i) = 0 
         do 77, j=i+1,2 
             p2(i,j) = p2(j,i) 
 77   continue
       
      zkappa = sqrt( dcmplx( p2(1,0)**2
     &                        +p2(2,1)**2
     &                        +p2(2,0)**2
     &                        -2.D0*( p2(1,0)*p2(2,1)
     &                             +p2(1,0)*p2(2,0)
     &                             +p2(2,0)*p2(2,1) ) ))

       if (real(zi*zkappa) .ne. 0.D0) then
         write(6,*)'warning !!!'
       endif 

      do 78 i=0,2
         j=mod(i+1,3)
         k=mod(j+1,3)

         zy0(i) = ( p2(j,k) * ( p2(j,k)
     &                         -p2(k,i)
     &                         -p2(i,j) 
     &                         +2.D0*m2(i)
     &                         -m2(j)
     &                         -m2(k)   )
     &          - ( p2(k,i) 
     &             -p2(i,j) ) * ( m2(j) 
     &                           -m2(k) )
     &          + zkappa*( p2(j,k) 
     &                    -m2(j)
     &                    +m2(k)  )                  
     &                                    )/( 2.D0 * zkappa * p2(j,k) )

         if (zy0(i) .eq. 1.D0) then
           zy0(i) = zy0(i)*(1.D0+eps)
         endif

         zalpha(i) = sqrt( dcmplx( p2(j,k)**2
     &                              +m2(j)**2
     &                              +m2(k)**2 
     &                              -2.D0*( p2(j,k)*m2(j)
     &                                   +p2(j,k)*m2(k)
     &                                   +m2(j)*m2(k)  ) ))
     &             *( 1.D0 
     &               +zi * eps * p2(j,k) )

         do 78 ip=0,1
             zx(i,ip) = ( p2(j,k)
     &                   -m2(j)
     &                   +m2(k)
     &                   +(-1)**ip * zalpha(i) ) /( 2.D0 * p2(j,k) )

             zy(i,ip) = zy0(i) 
     &                 -zx(i,ip)
 78      continue   

      zc = dcmplx(0.D0,0.D0)

      do 87, i=0,2
         j=mod(i+1,3)
         k=mod(j+1,3)

         zc1 = dcmplx(0.D0,0.D0)
         do 88, ip=0,1
            zc2 = Kspenc( (zy0(i)-1.D0)/zy(i,ip) ) 
     &           -Kspenc( zy0(i)/zy(i,ip) ) 
     &           +detafct( 1D0-zx(i,ip) , 1.D0/zy(i,ip) ) 
     &                             * log( (zy0(i)-1.D0)/zy(i,ip) )
     &           -detafct( -zx(i,ip) , 1.D0/zy(i,ip) )
     &                             * log( zy0(i)/zy(i,ip) ) 

 88      zc1 = zc1 + zc2 

         zc1 = zc1 
     &       - ( detafct( -zx(i,0) , -zx(i,1) ) 
     &          -detafct( zy(i,0) , zy(i,1) )
     &          -2.D0*pi*zi*Dthetafct( -p2(j,k) )
     &                   *Dthetafct( -real( -zi*zy(i,0)*zy(i,1) ) )  )
     &       * log( (1.D0-zy0(i))/(-zy0(i)) ) 
 87   zc = zc + zc1 

      Cc = zc/zkappa

      return
      end


c -----------------------------------
      complex*16 function detafct(z1,z2)
      implicit real*8 (a-h,o-y), complex*16 (z)

      pi=4.D0*atan(1.D0)

      zi=dcmplx(0.D0,1.D0)
      x1=aimag(z1)
      x2=aimag(z2)
      x12=aimag(z1*z2)

      detafct=2.D0*pi*zi*( 
     &     Dthetafct(-x1) * Dthetafct(-x2) * Dthetafct(x12)
     &     -Dthetafct(x1)  * Dthetafct(x2)  * Dthetafct(-x12) )

      return
      end


c ----------------------------------------
      real*8 function Dthetafct(x)
      implicit real*8 (a-h,o-y), complex*16 (z)
  
      if (x.gt.1.d-99) then
        Dthetafct=1.D0
      else
        Dthetafct=0.D0
      endif

      return
      end            





















