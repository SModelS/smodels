CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE NEUTRALINO GLUINO ASSOCIATED HADROPRODUCTION            C
C                                                                     C
C       SCALAR_ARRAY_NS(MASSIN,SCA,SCB,SCBP,SCC,SCD)                  C
C                                                                     C
C    INPUT  : MASSIN(1:30)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = T2                                         C
C             MASSIN(6)  = MS                                         C
C             MASSIN(7)  = MN                                         C
C             MASSIN(9)  = MT                                         C
C             MASSIN(10) = MG                                         C
C             MASSIN(11) = MS                                         C
C                                                                     C
C    OUTPUT : SCA(1:10)                                               C
C                                                                     C
C             SCA(2) = A_FIN(MS)                                      C
C             SCA(3) = A_FIN(MG)                                      C
C                                                                     C
C             SCB(1:10,1:5)                                           C
C                                                                     C
C             SCB(1,2) = B_FIN(U,MG,0)                                C
C             SCB(1,3) = B_FIN(MS**2,MG,0)                            C
C             SCB(3,2) = B_FIN(U,MS,0)                                C
C             SCB(3,4) = B_FIN(MS**2,MS,0)                            C
C             SCB(3,5) = B_FIN(MN**2,MS,0)                            C
C             SCB(5,1) = B_FIN(S,MG,MS)                               C
C             SCB(6,1) = B_FIN(0,MG,MS)                               C
C             SCB(6,2) = B_FIN(0,MT,MT)                               C
C             SCB(6,3) = B_FIN(0,MS,MS)                               C
C             SCB(6,4) = B_FIN(0,MG,MG)                               C
C             SCB(7,1) = B_FIN(S,0,0)                                 C
C                                                                     C
C             SCBP(1:10)                                              C
C                                                                     C
C             SCBP(1) = BP_FIN(0,MG,MS)                               C
C             SCBP(5) = BP_FIN(MS**2,MS,0)                            C
C             SCBP(6) = BP_FIN(MS**2,MG,0)                            C
C                                                                     C
C             SCC(1:20,1:5)                                           C
C                                                                     C
C             SCC(1,4) = C_FIN(P2,K1,0,MS,MG)                         C
C             SCC(2,1) = C_FIN(P1,K1,MS,0,0)                          C
C             SCC(2,2) = C_FIN(P2,K2,MS,0,0)                          C
C             SCC(2,3) = C_FIN(P1,K2,MS,0,0)                          C
C             SCC(2,4) = C_FIN(P2,K1,MS,0,0)                          C
C             SCC(2,5) = C_FIN(P1,K2,MG,0,0)                          C
C             SCC(4,1) = C_FIN(K1,K2,0,0,0)                           C
C             SCC(5,1) = C_FIN(P1,P2,0,MS,0)                          C
C             SCC(8,1) = C_FIN(P1,K1,0,MG,MS)                         C
C             SCC(9,2) = C_FIN(P2,K2,0,MS,MS)                         C
C             SCC(9,3) = C_FIN(P1,K2,0,MS,MS)                         C
C             SCC(10,1)= C_FIN(P1,P2,MG,0,MS)                         C
C             SCC(11,3)= C_FIN(P1,K2,0,MG,MG)                         C
C             SCC(12,1)= C_FIN(K1,K2,MS,MG,MG)                        C
C             SCC(13,1)= C_FIN(K1,K2,MG,MS,MS)                        C
C                                                                     C
C             SCD(1:10,1:2)                                           C
C                                                                     C
C             SCD(1,1) = D_FIN(K2,K1,P1,0,0,0,MS)                     C
C             SCD(1,2) = D_FIN(K2,K1,P2,0,0,0,MS)                     C
C             SCD(4,1) = D_FIN(K2,P1,K1,MS,MS,0,0)                    C
C             SCD(5,1) = D_FIN(K2,K1,P1,MS,MS,MG,0)                   C
C             SCD(6,2) = D_FIN(K2,K1,P2,MG,MG,MS,0)                   C
C             SCD(7,1) = D_FIN(K2,P1,K1,0,0,MG,MS)                    C
C                                                                     C
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
C    LAST MODIFIED : 23.02.99                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_NS(massin,sca,scb,scbp,scc,scd)
      
      implicit none 

      integer    n,n1,n2

      real*8     massin(1:30)
     &          ,sca(1:10),scb(1:10,1:5),scbp(1:10)
     &          ,scc(1:20,1:5),scd(1:10,1:2)
     &          ,s,t,m1,m2,mg,ms,mt,u,mu2
     &          ,zeta2,eps1,eps2,epsi,Li2,B02,BP02
     &          ,duma,dumc_1,dumc_3,dumd_1,dumd_3
     &          ,mjs2,msg2,mu,ts,ug,x1,x2

      complex*16 C03,D04,CSPEN,m1c2,m2c2,msc2,mgc2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               finite part of one point function 
      duma(mt,mu2) = mt**2*( 1.D0 - log(mt**2/mu2) )
      
c               divergent three point functions, first roland's (0.34)
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
      dumd_1(t,m1,m2,m1c2,s,mu2,zeta2) =
     &  ( - 2.D0 * Li2( (m2**2-t)/(m1**2-t) )
     &    - real( log( (m1c2-m2**2)/m1/sqrt(mu2) )**2 )
     &    + 2.D0*log( (m1**2-t)/m1/sqrt(mu2) )*log(s/mu2) 
     &    - 13.D0 * zeta2 /4.D0 )
     &   /s/(t-m1**2)

c               roland's (0.49)      
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
      s  = massin(1)
      m1 = abs(massin(6)) 
      m2 = abs(massin(7))
      mt = massin(9)
      mg = massin(10)
      ms = massin(11)

c               fix the scale to the value in the overall prefactor 
      mu2 = m1**2
      mu  = sqrt( mu2 )

      t    = massin(2) + m2**2 
      u    = m1**2 + m2**2 - s - t 
      ts   = t - ms**2 
      ug   = u - mg**2 
      mjs2 = m2**2 - ms**2 
      msg2 = ms**2 - mg**2 

c               cut-off for external momenta in C03 and D04 [non-zero 1.D-7]
      eps1 = 1.D-7

c               cut-off for masses in C03 and D04 [1.D-8]
      eps2 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               set everything to 9.99999D+9
      do n=1,3,1
         sca(n) = 9.99999D+9 
      end do

      do n1=1,10,1
         do n2=1,5,1
            scb(n1,n2) = 9.99999D+9
         end do
      end do
      
      do n=1,10,1
         scbp(n) = 9.99999D+9 
      end do

      do n1=1,20,1
         do n2=1,5,1
            scc(n1,n2) = 9.99999D+9
         end do
      end do

      do n1=1,10,1
         do n2=1,2,1
            scd(n1,n2) = 9.99999D+9
         end do
      end do

c               the one point functions (finite parts)
      sca(2) = duma(ms,mu2)
      sca(3) = duma(mg,mu2)

c               the two point functions (finite parts)
      scb(1,2) = B02(u,mg,0.D0,mu2)
      scb(1,3) = B02(m1**2,mg,0.D0,mu2)
      scb(3,2) = B02(u,ms,0.D0,mu2)
      scb(3,4) = B02(m1**2,ms,0.D0,mu2)
      scb(3,5) = B02(m2**2,ms,0.D0,mu2)
      scb(5,1) = B02(s,mg,ms,mu2)
      scb(6,1) = B02(0.D0,mg,ms,mu2)
      scb(6,2) = B02(0.D0,mt,mt,mu2)
      scb(6,3) = B02(0.D0,ms,ms,mu2)
      scb(6,4) = B02(0.D0,mg,mg,mu2)
      scb(7,1) = B02(s,0.D0,0.D0,mu2)

c               derived two point function
      scbp(1) = BP02(0.D0,mg,ms,mu2)
      scbp(5) = BP02(ms**2,ms,0.D0,mu2)
      scbp(6) = BP02(ms**2,mg,0.D0,mu2)

c               finite three point functions 
      scc(1,4)  = real( C03(m2**2,eps1,u,eps2,ms,mg) )
      scc(3,1)  = real( C03(eps1,eps1,s,ms,m2,ms) )
      scc(5,1)  = real( C03(m1**2,m2**2,s,eps2,ms,eps2) )
      scc(8,1)  = real( C03(m1**2,eps1,t,eps2,mg,ms) )
      scc(9,2)  = real( C03(m2**2,eps1,t,eps2,ms,ms) )
      scc(9,3)  = real( C03(m1**2,eps1,u,eps2,ms,ms) )
      scc(10,1) = real( C03(m1**2,m2**2,s,mg,eps2,ms) )
      scc(11,3) = real( C03(m1**2,eps1,u,eps2,mg,mg) )
      scc(12,1) = real( C03(eps1,eps1,s,ms,mg,mg) )
      scc(13,1) = real( C03(eps1,eps1,s,mg,ms,ms) )

c               divergent three point functions 
      msc2 = ms**2 * dcmplx(1.D0,-epsi) 
      mgc2 = mg**2 * dcmplx(1.D0,-epsi) 

      scc(4,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
      scc(2,1) = dumc_3(t,ms,msc2,mu2,zeta2)
      scc(2,2) = dumc_1(t,m2,ms,msc2,mu2)
      scc(2,3) = dumc_3(u,ms,msc2,mu2,zeta2)
      scc(2,4) = dumc_1(u,m2,ms,msc2,mu2)
      scc(2,5) = dumc_1(u,m1,mg,mgc2,mu2)

c               finite four point functions
      scd(5,1) = real( D04(eps1,eps1,m1**2,m2**2,s,t,ms,ms,mg,eps2) )
      scd(6,2) = real( D04(eps1,eps1,m2**2,m1**2,s,u,mg,mg,ms,eps2) )

c               divergent four point functions 
      scd(1,1) = dumd_1(t,ms,m2,msc2,s,mu2,zeta2)
      scd(1,2) = dumd_1(u,ms,m2,msc2,s,mu2,zeta2)
      scd(4,1) = dumd_3(t,u,m1,m2,mu2,zeta2)

c               the new integral
      x1 = mjs2/ug
      x2 = msg2/ts

c               note: only real part wanted!
      scd(7,1) =  2.D0 * log( -ts/ms/mu ) * log( abs(-ug/mg/mu) )
     &          - ( abs( ( log( dcmplx(-mjs2/ms/mu) ) )**2 )
     &             -abs( ( log( dcmplx(-msg2/mg/mu) ) )**2 ) )
     &          - ( log(mg/ms) )**2 
     &          - 2.D0 * Li2(1.D0-mjs2/ts) - 2.D0 * Li2(1.D0-msg2/ug)
     &          - 2.D0 * Li2(1.D0-msg2/ts) - 2.D0 * Li2(1.D0-mjs2/ug)
     &          - Li2(1.D0-ms**2/mg**2*msg2/ts)
     &          - Li2(1.D0-mg**2/ms**2*mjs2/ug)
     &          + 2.D0 *(   Li2(1.D0-x1*x2)
     &                   + log(1.D0-x1*x2) * (  log( abs(x1*x2) )
     &                                         - log( abs(x1) )
     &                                         - log( abs(x2) )   ) )
      scd(7,1) = scd(7,1)/( ts*ug - mjs2*msg2 ) 

      return 
      end 

