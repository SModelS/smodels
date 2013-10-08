CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE NEUTRALINO PAIR HADROPRODUCTION                         C
C                                                                     C
C       SCALAR_ARRAY_NN(MASSIN,SCB,SCBP,SCC,SCD)                      C
C                                                                     C
C    INPUT  : MASSIN(1:20)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = T2                                         C
C             MASSIN(3)  = U2                                         C
C             MASSIN(4)  = T1                                         C
C             MASSIN(5)  = U1                                         C
C             MASSIN(6)  = MN1                                        C
C             MASSIN(7)  = MN2                                        C
C             MASSIN(8)  = MZ                                         C
C             MASSIN(9)  = MT                                         C
C             MASSIN(10) = MG                                         C
C             MASSIN(11) = MS                                         C
C             MASSIN(12) = MU                                         C
C                                                                     C
C    OUTPUT : SCB(1:7,1:5)                                            C
C                                                                     C
C             SCB(1,1) = B_FIN(T,MG,0)                                C
C             SCB(1,2) = B_FIN(U,MG,0)                                C
C             SCB(1,3) = B_FIN(MS**2,MG,0)                            C
C             SCB(3,1) = B_FIN(T,MS,0)                                C
C             SCB(3,2) = B_FIN(U,MS,0)                                C
C             SCB(3,3) = B_FIN(MS,MS,0)                               C
C             SCB(3,4) = B_FIN(MN1**2,MS,0)                           C
C             SCB(3,5) = B_FIN(MN2**2,MS,0)                           C
C             SCB(4,1) = B_FIN(S,MS,MS)                               C
C             SCB(6,1) = B_FIN(0,MG,MS)                               C
C             SCB(7,1) = B_FIN(S,0,0)                                 C
C                                                                     C
C             SCBP(1:1)                                               C
C                                                                     C
C             SCBP(1) = BP_FIN(0,MG,MS)                               C
C                                                                     C
C             SCC(1:6,1:4)                                            C
C                                                                     C
C             SCC(1,1) = S  * C_FIN(P1,K1,0,MS,MG)                    C
C             SCC(1,2) = S  * C_FIN(P2,K2,0,MS,MG)                    C
C             SCC(1,3) = S  * C_FIN(P1,K2,0,MS,MG)                    C
C             SCC(1,4) = S  * C_FIN(P2,K1,0,MS,MG)                    C
C             SCC(2,1) = T1 * C_FIN(P1,K1,MS,0,0)                     C
C             SCC(2,2) = T2 * C_FIN(P2,K2,MS,0,0)                     C
C             SCC(2,3) = U1 * C_FIN(P1,K2,MS,0,0)                     C
C             SCC(2,4) = U2 * C_FIN(P2,K1,MS,0,0)                     C
C             SCC(3,1) = S  * C_FIN(K1,K2,MS,MG,MS)                   C
C             SCC(4,1) = S  * C_FIN(K1,K2,0,0,0)                      C
C             SCC(5,1) = S  * C_FIN(P1,P2,0,MS,0)                     C
C             SCC(6,1) = S  * C_FIN(P1,P2,MS,0,MS)                    C
C                                                                     C
C             SCD(1:2,1:2)                                            C
C                                                                     C
C             SCD(1,1) = S*TS * D_FIN(K2,K1,P1,0,0,0,MS)              C
C             SCD(1,2) = S*US * D_FIN(K2,K1,P2,0,0,0,MS)              C
C             SCD(2,1) = S*TG * D_FIN(K2,K1,P1,MS,MG,MS,0)            C
C             SCD(2,2) = S*UG * D_FIN(K2,K1,P2,MS,MG,MS,0)            C
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
C    LAST MODIFIED : 01.12.97                                         C
C                    08.05.09 INTRODUCED EPSI, NOTHING CHANGED        C
C                    08.05.09 INCREASED EPS1 FOR C                    C
C                             SOLVED PROBLEM IN SOMMERFELD REGION     C
C                             FOR 2 TEV SQUARK MASSES SHARP           C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_NN(massin,scb,scbp,scc,scd)
      
      implicit none 

      integer    n1,n2

      real*8     massin(1:20)
     &          ,scb(1:7,1:5),scbp(1:1),scc(1:6,1:4),scd(1:2,1:2)
     &          ,s,t,mn1,mn2,mg,ms,u,tg,ug,mu2 
     &          ,zeta2,eps1,eps2,epsi,Li2,B02,BP02,dumc,dumd
      complex*16 C03,D04,CSPEN,sc,msc2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               divergent three point function 
      dumc(t,mn1,mn2,ms,msc2,mu2) =  
     &    Li2(t/ms**2) - Li2(mn1**2/ms**2)
     &  + real( log( 1.D0 - t/msc2      )**2 )  
     &  - real( log( 1.D0 - mn1**2/msc2 )**2 )
     &  + log(ms**2/mu2) * log( abs( -(t-ms**2)/(ms**2-mn1**2) ) )

c               divergent four point function 
      dumd(t,mn1,mn2,ms,msc2,s,sc,mu2,zeta2) =
     &  - 2.D0 * Li2( 1.D0 + (ms**2-mn1**2)/(t-ms**2) )
     &  - 2.D0 * Li2( 1.D0 + (ms**2-mn2**2)/(t-ms**2) )
     &  - Li2( 1.D0 + (ms**2-mn1**2)*(ms**2-mn2**2)/s/ms**2 )
     &  - 3.D0/2.D0 * zeta2
     &  - real( log( 1.D0 + (msc2-mn1**2)*(msc2-mn2**2)/s/msc2 )
     &           * (   log( -(msc2-mn1**2)*(msc2-mn2**2)/s/msc2 )
     &               - log( (msc2-mn1**2)/mu2)
     &               - log( (msc2-mn2**2)/mu2)
     &               + log( -s*msc2/mu2**2)         )              )
     &  + log( s/mu2  )**2 /2.D0
     &  - log( s/ms**2)**2 /2.D0
     &  + 2.D0*real( log( -sc/mu2 )
     &               *log( -(t-msc2)/msc2) )
     &  - real( log( (msc2-mn1**2)/mu2 )
     &          *log( (msc2-mn1**2)/msc2) )
     &  - real( log( (msc2-mn2**2)/mu2 )
     &          *log( (msc2-mn2**2)/msc2) )
     
c               constants 
      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0

c               define the different variables and fix the scale 
      s   = massin(1)
      mn1 = massin(6) 
      mn2 = massin(7)
      mg  = massin(10)
      ms  = massin(11)
c      mu2 = massin(12)**2
      mu2 = ( abs(mn1) + abs(mn2) )**2 / 4.D0

      t   = massin(2) + mn2**2 

      u   = mn1**2 + mn2**2 - s - t 
      tg  = t - mg**2 
      ug  = u - mg**2 

c               cut-off for external momenta in C03 and D04 [non-zero 1.D-7]
      eps1 = 1.D-7

c               cut-off for masses in C03 and D04 [1.D-8]
      eps2 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               set everything to 9.D-9
      do 9999, n1=1,7,1
      do 9999, n2=1,5,1
 9999    scb(n1,n2) = 9.D-9

      scbp(1) = 9.D-9 

      do 9998 n1=1,6,1
      do 9998 n2=1,4,1
 9998    scc(n1,n2) = 9.D-9

      do 9997 n1=1,2,1
      do 9997 n2=1,2,1
 9997    scd(n1,n2) = 9.D-9

c               the two point functions (finite parts)
      scb(1,1) = B02(t,mg,0.D0,mu2)
      scb(1,2) = B02(u,mg,0.D0,mu2)
      scb(1,3) = B02(ms**2,mg,0.D0,mu2)
      scb(3,1) = B02(t,ms,0.D0,mu2)
      scb(3,2) = B02(u,ms,0.D0,mu2)
      scb(3,3) = B02(ms**2,ms,0.D0,mu2)
      scb(3,4) = B02(mn1**2,ms,0.D0,mu2)
      scb(3,5) = B02(mn2**2,ms,0.D0,mu2)
      scb(4,1) = B02(s,ms,ms,mu2)
      scb(6,1) = B02(0.D0,mg,ms,mu2)
      scb(7,1) = B02(s,0.D0,0.D0,mu2)

c               derived two point function
      scbp(1) = BP02(0.D0,mg,ms,mu2)

c               finite three point functions 
      scc(1,1) = s  * real( C03(mn1**2,eps1,t,eps2,ms,mg) )
      scc(1,2) = s  * real( C03(mn2**2,eps1,t,eps2,ms,mg) )
      scc(1,3) = s  * real( C03(mn1**2,eps1,u,eps2,ms,mg) )
      scc(1,4) = s  * real( C03(mn2**2,eps1,u,eps2,ms,mg) )
      scc(3,1) = s  * real( C03(eps1,eps1,s,ms,mg,ms) )
      scc(5,1) = s  * real( C03(mn1**2,mn2**2,s,eps2,ms,eps2) )
      scc(6,1) = s  * real( C03(mn1**2,mn2**2,s,ms,eps2,ms) )

c               divergent three point functions 
      msc2 = ms**2 * dcmplx(1.D0,-epsi) 

      scc(4,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )
      scc(2,1) = dumc(t,abs(mn1),abs(mn2),ms,msc2,mu2)
      scc(2,2) = dumc(t,abs(mn2),abs(mn1),ms,msc2,mu2)
      scc(2,3) = dumc(u,abs(mn1),abs(mn2),ms,msc2,mu2)
      scc(2,4) = dumc(u,abs(mn2),abs(mn1),ms,msc2,mu2)

c               finite four point functions
      scd(2,1) = s*tg
     &          *real( D04(eps1,eps1,mn1**2,mn2**2,s,t,ms,mg,ms,eps2) )
      scd(2,2) = s*ug
     &          *real( D04(eps1,eps1,mn2**2,mn1**2,s,u,ms,mg,ms,eps2) )


c               divergent four point functions 
      sc = s * dcmplx(1.D0,epsi) 

      scd(1,1) = dumd(t,abs(mn1),abs(mn2),ms,msc2,s,sc,mu2,zeta2)
      scd(1,2) = dumd(u,abs(mn2),abs(mn1),ms,msc2,s,sc,mu2,zeta2)

      return 
      end 

