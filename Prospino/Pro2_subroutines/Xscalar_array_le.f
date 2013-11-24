CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE ASSOCIATED STOP-LEPTON PRODUCTION                       C
C                                                                     C
C       SCALAR_ARRAY(MASSIN,SCA,SCB,SCBP,SCC,SCD)                     C
C                                                                     C
C    INPUT  : MASSIN(1:30)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = T2                                         C
C             MASSIN(6)  = M_S1                                       C
C             MASSIN(7)  = 0                                          C
C             MASSIN(8)  = M_GL                                       C
C             MASSIN(9)  = M_SQ                                       C
C             MASSIN(10) = M_S2                                       C
C             MASSIN(11) = M_T                                        C
C                                                                     C
C    OUTPUT:                                                          C
C                                                                     C
C       SCB(1:10,1:6)                                                 C
C                                                                     C
C        B_fin(p2,m1,0)         = SCB(3,1)                            C
C        B_fin(p1 + k2,m1,0)    = SCB(3,2)                            C
C        B_fin(p1,m1,0)         = SCB(3,4)                            C
C        B_fin(k1 + k2,mg,msx)  = SCB(5,1)                            C
C        B_fin(p1 + k2,mg,mt)   = SCB(5,2)                            C
C        B_fin(p1,mg,mt)        = SCB(5,3)                            C
C        B_fin(k1,mg,msx)       = SCB(6,1)                            C
C        B_fin(k1,mt,mt)        = SCB(6,2)                            C
C        B_fin(k1,m1,m1)        = SCB(6,3)                            C
C        B_fin(k1,mg,mg)        = SCB(6,4)                            C
C        B_fin(k1,msx,msx)      = SCB(6,5)                            C
C        B_fin(k1,ms2,ms2)      = SCB(6,6)                            C
C        B_fin(k1 + k2,0,0)     = SCB(7,1)                            C
C                                                                     C
C       SCBP(1:10)                                                    C
C                                                                     C
C        Bp_fin(k1,mg,msx)      = SCBP(1)                             C
C        Bp_fin(p1,m1,0)        = SCBP(5)                             C
C        Bp_fin(p1,mg,mt)       = SCBP(6)                             C
C                                                                     C
C       SCC(1:20,1:4)                                                 C
C                                                                     C
C        C_fin(p1,k1,m1,0,0)    = SCC(2,1)                            C
C        C_fin(p2,k2,m1,0,0)    = SCC(2,2)                            C
C        C_fin(p1,k2,m1,0,0)    = SCC(2,3)                            C
C        C_fin(p2,k1,m1,0,0)    = SCC(2,4)                            C
C        C_fin(k1,k2,0,0,0)     = SCC(4,1)                            C
C        C_fin(p1,p2,0,m1,0)    = SCC(5,1)                            C
C        C_fin(p2,k2,0,m1,m1)   = SCC(9,2)                            C
C        C_fin(p1,k2,0,m1,m1)   = SCC(9,3)                            C
C        C_fin(k1,k2,mg,msx,msx)= SCC(13,1)                           C 
C        C_fin(k1,k2,msx,mg,mg) = SCC(12,1)                           C 
C                                                                     C
C       SCD(1:10,1:2)                                                 C
C                                                                     C
C        D_fin(k1,k2,p2,0,0,0,m1)  = SCD(1,1)                         C
C        D_fin(k2,k1,p2,0,0,0,m1)  = SCD(1,2)                         C
C        D_fin(k1,p2,k2,0,0,m1,m1) = SCD(4,1)                         C
C                                                                     C
C        SOF1(s4^2)  = SOF1(1)                                        C
C        SOF1(s4*tp) = SOF1(2)                                        C
C        SOF1(s4*up) = SOF1(3)                                        C
C        SOF1(tp*up) = SOF1(4)                                        C
C
C                                                                     C
C   N.B. : THE INTEGRAL MEASURE IS D^NQ/(I*Pi^2)                      C
C          ALL VARIABLES REAL*8, ONLY THE REAL PARTS OF THE INTEGRALS C
C                                                                     C
C                                                                     C
C    NECESSARY SUBROUTINES :  XSPENCE.F                               C
C                             XWOPOINT_SPIRA.F                        C
C                             XTHREEPOINT_SPIRA.F                     C
C                             XDENNER_3.F                             C
C                                                                     C
C    CALLED FUNCTIONS : B02                                           C
C                       BP02                                          C
C                       C03                                           C
C                       D04                                           C
C                       CSPEN                                         C
C                                                                     C
C    LAST MODIFIED : 24.05.01                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine SCALAR_ARRAY_LE(massin,scb,scbp,scc,scd)
      
      implicit none 

      integer    n,n1,n2

      real*8     massin(1:30)
     &          ,scb(1:10,1:6),scbp(1:10)
     &          ,scc(1:20,1:4),scd(1:10,1:2)
     &          ,s,t,m1,m2,mg,msx,ms2,mt,u,mu2
     &          ,zeta2,eps1,eps2,epsi,Li2,B02,BP02
     &          ,dumc_1,dumc_3,dumd_1,dumd_3

      complex*16 C03,D04,CSPEN,m1c2,m2c2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               divergent three point functions, first roland's (0.34)
c                check: called for m1=0, always m2>0
      dumc_1(t,m1,m2,m2c2,mu2) =  
     &  (   Li2(t/m2**2) - Li2(m1**2/m2**2)
     &    + real( log( 1.D0 - t/m2c2     )**2 )  
     &    - real( log( 1.D0 - m1**2/m2c2 )**2 )
     &    + log(m2**2/mu2) * log( abs(-(t-m2**2)/(m2**2-m1**2) ) ) )
     &  /(t-m1**2)

c               same as above for m1=ms, roland's (0.33)
      dumc_3(t,m1,m1c2,mu2,zeta2) =  
     &  (   Li2(t/m1**2)
     &    + real( log( 1.D0 - t/m1c2 )**2 )  
     &    + 1.D0/4.D0 * log(mu2/m1**2)**2 
     &    - log(mu2/m1**2)*log(abs(-(t-m1**2)/m1**2)) 
     &    + zeta2/4.D0                            )
     &  /(t-m1**2)

c               divergent four point, roland's (0.46)   
c                check: called for m2=0, always m1>0
      dumd_1(t,m1,m2,m1c2,s,mu2,zeta2) =
     &  ( - 2.D0 * Li2( (m2**2-t)/(m1**2-t) )
     &    - real( log( (m1c2-m2**2)/m1/sqrt(mu2) )**2 )
     &    + 2.D0*log( (m1**2-t)/m1/sqrt(mu2) )*log(s/mu2) 
     &    - 13.D0 * zeta2 /4.D0 )
     &   /s/(t-m1**2)

c               roland's (0.49)      
c                check: called for m2=0, always m1>0
      dumd_3(t,u,m1,m2,mu2,zeta2) = 
     &  (  2.D0 * Li2((t-m2**2)/(m1**2-m2**2))
     &   - 2.D0 * Li2((m2**2-u)/(m1**2-u))
     &   + log((m1**2-t)/m1/sqrt(mu2))**2 
     &   + 2.D0 * log((m1**2-t)/m1/sqrt(mu2))
     &          * log(abs((m1**2-u)/(m1**2-m2**2)))
     &   - 6.D0*zeta2/8.D0                           )
     &   /(t-m1**2)/(u-m1**2)

c               constants 
      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0

c               define the different variables and fix the scale 
      s   = massin(1)
      m1  = massin(6)
      m2  = 0.D0
      mg  = massin(8)
      msx = massin(9)
      ms2 = massin(10)
      mt  = massin(11)

c               fix the scale to the value in the overall prefactor 
      mu2 = m1**2

      t   = massin(2) + m2**2 
      u   = m1**2 + m2**2 - s - t 

c               cut-off for external momenta in C03 and D04 [non-zero 1.D-7]
      eps1 = 1.D-7

c               cut-off for masses in C03 and D04 [1.D-8]
      eps2 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               set everything to 9.99999D+9
      do n1=1,10,1
         do n2=1,6,1
            scb(n1,n2) = 9.99999D+9
         end do
      end do
      
      do n=1,10,1
         scbp(n) = 9.99999D+9 
      end do

      do n1=1,20,1
         do n2=1,4,1
            scc(n1,n2) = 9.99999D+9
         end do
      end do

      do n1=1,10,1
         do n2=1,2,1
            scd(n1,n2) = 9.99999D+9
         end do
      end do

c               the two point functions
      scb(3,1) = B02(m2**2,m1  ,0.D0,mu2)
      scb(3,2) = B02(u    ,m1  ,0.D0,mu2)
      scb(3,4) = B02(m1**2,m1  ,0.D0,mu2)
      scb(5,1) = B02(s    ,mg  ,msx ,mu2)
      scb(5,2) = B02(u    ,mg  ,mt  ,mu2)
      scb(5,3) = B02(m1**2,mg  ,mt  ,mu2)
      scb(6,1) = B02(0.D0 ,mg  ,msx ,mu2)
      scb(6,2) = B02(0.D0 ,mt  ,mt  ,mu2)
      scb(6,3) = B02(0.D0 ,m1  ,m1  ,mu2)
      scb(6,4) = B02(0.D0 ,mg  ,mg  ,mu2)
      scb(6,5) = B02(0.D0 ,msx ,msx ,mu2)
      scb(6,6) = B02(0.D0 ,ms2 ,ms2 ,mu2)
      scb(7,1) = B02(s    ,0.D0,0.D0,mu2)

c               derived two point function
      scbp(1) = BP02(0.D0 ,mg  ,msx ,mu2)
      scbp(5) = BP02(m1**2,m1  ,0.D0,mu2)
      scbp(6) = BP02(m1**2,mg  ,mt  ,mu2)

c               finite three point functions 
      scc(5,1)  = real( C03(m1**2,m2**2,s,eps2,m1 ,eps2) )
      scc(9,2)  = real( C03(m2**2,eps1 ,t,eps2,m1 ,m1  ) )
      scc(9,3)  = real( C03(m1**2,eps1 ,u,eps2,m1 ,m1  ) )
      scc(12,1) = real( C03(eps1 ,eps1 ,s,msx ,mg ,mg  ) )
      scc(13,1) = real( C03(eps1 ,eps1 ,s,mg  ,msx,msx ) )

c               divergent three point functions 
      m1c2 = m1**2 * dcmplx(1.D0,-epsi) 

      scc(4,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
      scc(2,1) = dumc_3(t,m1,m1c2,mu2,zeta2)
      scc(2,3) = dumc_3(u,m1,m1c2,mu2,zeta2)
      scc(2,2) = dumc_1(t,m2,m1,m1c2,mu2)
      scc(2,4) = dumc_1(u,m2,m1,m1c2,mu2)

c               divergent four point functions 
      scd(1,1) = dumd_1(t,m1,m2,m1c2,s,mu2,zeta2)
      scd(1,2) = dumd_1(u,m1,m2,m1c2,s,mu2,zeta2)
      scd(4,1) = dumd_3(t,u,m1,m2,mu2,zeta2)

      return 
      end 




