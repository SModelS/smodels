CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C        SUBROUTINE CALCULATING THE FINITE REAL PART OF THE           C
C          GENERAL MASSIVE TWO POINT FUNCTION                         C
C                                                                     C
C           B02(P.P,M1,M2,MU**2)                                      C
C           BP02(P.P,M1,M2,MU**2)                                     C
C                                                                     C
C        LAST CHANGE:                                                 C
C         24.02.99 [TP] : B02P CASES                                  C
C         28.05.01 [TP] : B(0,M,0)                                    C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c ---------------------------------------------------------------------
      real*8 function B02(s,m1,m2,mu2)

      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2 

      m12 = m1**2 
      m22 = m2**2 

      zkappa=sqrt(dcmplx(s**2+m12**2+m22**2
     &                     -2.D0*(s*m12+s*m22+m12*m22)))

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            B02=-log(m12/mu2) 
         else
            if ((m12.ne.0.D0).and.(m22.ne.0.D0)) then 
               B02=1.D0 - m12/(m12-m22)*log(m12/mu2)
     &                  + m22/(m12-m22)*log(m22/mu2) 
            elseif (m12.eq.0.D0) then 
               B02=1.D0 - log(m22/mu2)
            elseif (m22.eq.0.D0) then
               B02=1.D0 - log(m12/mu2)
            end if 
         endif
      else 
         if ((m12.eq.0.D0).and.(m22.eq.0.D0)) then 
            B02=2.D0 - log(s/mu2)
         elseif ((m12.eq.s).and.(m22.eq.0.D0)) then 
            B02=2.D0 - log(m12/mu2)
         elseif ((m22.eq.s).and.(m12.eq.0.D0)) then 
            B02=2.D0 - log(m22/mu2)
         elseif (m12.eq.0.D0) then
            B02=2.D0 - (s-m22)/s*log( abs(m22-s)/m22 )
     &                 - log(m22/mu2)
         elseif (m22.eq.0.D0) then
            B02=2.D0 - (s-m12)/s*log( abs(m12-s)/m12 ) 
     &                 - log(m12/mu2)
         else
            x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
            x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
            B02=real( 2.D0 + log(mu2/m22) 
     &                       + x1*log(1.D0-1.D0/x1) 
     &                       + x2*log(1.D0-1.D0/x2))
         endif
      endif 

      return
      end



c ---------------------------------------------------------------------
      real*8 function BP02(s,m1,m2,mu2)
      
      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2
      
      m12 = m1**2
      m22 = m2**2 

      zkappa=sqrt(dcmplx(s**2+m12**2+m22**2
     &                    -2.D0*(s*m12+s*m22+m12*m22)))

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            BP02=1.D0/(6.D0*m12)
         else
            BP02=( (m12+m22)/2.D0 
     &        - m12*m22/(m12-m22)*log(m12/m22) )/(m12-m22)**2 
         endif
      elseif ((s.eq.m12).and.(m22.eq.0.D0)) then 
         BP02=( -1.D0 + log(m12/mu2)/2.D0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0.D0)) then 
         BP02=( -1.D0 + log(m22/mu2)/2.D0 )/m22
      elseif ((m12.eq.0.D0).and.(m22.ne.0.D0)) then 
         BP02=( -1.D0 - m22/s*log(abs(m22-s)/m22) )/s  
      elseif ((m22.eq.0.D0).and.(m12.ne.0.D0)) then 
         BP02=( -1.D0 - m12/s*log(abs(m12-s)/m12) )/s  
      else 
         x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
         x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
         BP02=real( -1.D0 + ( x1*(1.D0-x1)*log(1.D0-1.D0/x1)
     &                       - x2*(1.D0-x2)*log(1.D0-1.D0/x2) )  
     &                                                  /(x1-x2) )/s
      endif 

      return
      end




