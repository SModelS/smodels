CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE NEUTRALINO GLUINO ASSOCIATED HADROPRODUCTION            C
C                                                                     C
C       SCALAR_ARRAY_NG(MASSIN,SCA,SCB,SCBP,SCC,SCD)                  C
C                                                                     C
C    INPUT  : MASSIN(1:30)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = TG                                         C
C             MASSIN(6)  = MN                                         C
C             MASSIN(7)  = MG                                         C
C             MASSIN(9)  = MT                                         C
C             MASSIN(11) = MS                                         C
C                                                                     C
C    OUTPUT : SCA(1:3)                                                C
C                                                                     C
C             SCA(1) = A_FIN(MT)                                      C
C             SCA(2) = A_FIN(MS)                                      C
C             SCA(3) = A_FIN(MG)                                      C
C                                                                     C
C             SCB(1:8,1:5)                                            C
C                                                                     C
C             SCB(1,1) = B_FIN(T,MG,0)                                C
C             SCB(1,2) = B_FIN(U,MG,0)                                C
C             SCB(1,3) = B_FIN(MS**2,MG,0)                            C
C             SCB(1,4) = B_FIN(MG**2,MG,0)                            C
C             SCB(3,1) = B_FIN(T,MS,0)                                C
C             SCB(3,2) = B_FIN(U,MS,0)                                C
C             SCB(3,3) = B_FIN(MS**2,MS,0)                            C
C             SCB(3,4) = B_FIN(MN**2,MS,0)                            C
C             SCB(3,5) = B_FIN(MG**2,MS,0)                            C
C             SCB(4,1) = B_FIN(S,MS,MS)                               C
C             SCB(6,1) = B_FIN(0,MG,MS)                               C
C             SCB(7,1) = B_FIN(S,0,0)                                 C
C             SCB(8,1) = B_FIN(MG**2,MT,MS)                           C
C                                                                     C
C             SCBP(1:4)                                               C
C                                                                     C
C             SCBP(1) = BP_FIN(0,MG,MS)                               C
C             SCBP(2) = BP_FIN(MG**2,MT,MS)                           C
C             SCBP(3) = BP_FIN(MG**2,MG,0)                            C
C             SCBP(4) = BP_FIN(MG**2,MS,0)                            C
C                                                                     C
C             SCC(1:8,1:4)                                            C
C                                                                     C
C             SCC(1,1) = C_FIN(P1,K1,0,MS,MG)                         C
C             SCC(1,2) = C_FIN(P2,K2,0,MS,MG)                         C
C             SCC(1,3) = C_FIN(P1,K2,0,MS,MG)                         C
C             SCC(1,4) = C_FIN(P2,K1,0,MS,MG)                         C
C             SCC(2,1) = C_FIN(P1,K1,MS,0,0)                          C
C             SCC(2,2) = C_FIN(P2,K2,MS,0,0)                          C
C             SCC(2,3) = C_FIN(P1,K2,MS,0,0)                          C
C             SCC(2,4) = C_FIN(P2,K1,MS,0,0)                          C
C             SCC(3,1) = C_FIN(K1,K2,MS,MG,MS)                        C
C             SCC(4,1) = C_FIN(K1,K2,0,0,0)                           C
C             SCC(5,1) = C_FIN(P1,P2,0,MS,0)                          C
C             SCC(6,1) = C_FIN(P1,P2,MS,0,MS)                         C
C             SCC(7,2) = C_FIN(P2,K2,MG,0,0)                          C
C             SCC(7,4) = C_FIN(P2,K1,MG,0,0)                          C
C             SCC(8,2) = C_FIN(P2,K2,0,MG,MS)                         C
C             SCC(8,4) = C_FIN(P2,K1,0,MG,MS)                         C
C                                                                     C
C             SCD(1:3,1:2)                                            C
C                                                                     C
C             SCD(1,1) = D_FIN(K2,K1,P1,0,0,0,MS)                     C
C             SCD(1,2) = D_FIN(K2,K1,P2,0,0,0,MS)                     C
C             SCD(2,1) = D_FIN(K2,K1,P1,MS,MG,MS,0)                   C
C             SCD(2,2) = D_FIN(K2,K1,P2,MS,MG,MS,0)                   C
C             SCD(3,1) = D_FIN(K2,P1,K1,0,0,MS,MG)                    C
C             SCD(3,2) = D_FIN(K1,P1,K2,0,0,MS,MG)                    C
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
C    LAST MODIFIED : 08.06.99                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_NG(massin,sca,scb,scbp,scc,scd)
      
      implicit none 

      integer    n,n1,n2

      real*8     massin(1:30)
     &          ,sca(1:3),scb(1:8,1:5),scbp(1:4)
     &          ,scc(1:8,1:4),scd(1:3,1:2)
     &          ,s,t,m1,mg,ms,mt,u,mu,mu2
     &          ,zeta2,eps1,eps2,epsi,Li2,B02,BP02
     &          ,duma,dumc_1,dumc_2,dumd_1,dumd_2

      complex*16 C03,D04,CSPEN,sc,msc2,mgc2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               finite part of one point function 
      duma(mt,mu2) = mt**2*( 1.D0 - log(mt**2/mu2) )
      
c               divergent three point functions, roland's (0.34) 
      dumc_1(t,m1,ms,msc2,mu2) =  
     &  (   Li2(t/ms**2) - Li2(m1**2/ms**2)
     &    + real( log( 1.D0 - t/msc2     )**2 )  
     &    - real( log( 1.D0 - m1**2/msc2 )**2 )
     &    + log(ms**2/mu2) * log( abs(-(t-ms**2)/(ms**2-m1**2) ) ) )
     &  /(t-m1**2)

c               watch the scales, roland's (0.33) gen. 
      dumc_2(t,mg,mgc2,mu2,zeta2) = 
     &  (   real( log(-(t-mg**2)/mgc2)**2 )
     &    + 1.D0/4.D0 * log(mu2/mg**2)**2 
     &    - log(mu2/mg**2)*log(abs(-(t-mg**2)/mg**2)) 
     &    + Li2(t/mg**2) + zeta2/4.D0          )
     &  /(t-mg**2)

c               divergent four point function, roland's (0.41) gen. 
      dumd_1(t,m1,mg,ms,msc2,s,sc,mu2,zeta2) =
     &  ( - 2.D0 * Li2( 1.D0 + (ms**2-m1**2)/(t-ms**2) )
     &    - 2.D0 * Li2( 1.D0 + (ms**2-mg**2)/(t-ms**2) )
     &    - Li2( 1.D0 + (ms**2-m1**2)*(ms**2-mg**2)/s/ms**2 )
     &    - 3.D0/2.D0 * zeta2
     &    - real( log( 1.D0 + (msc2-m1**2)*(msc2-mg**2)/s/msc2 )
     &             * (   log( -(msc2-m1**2)*(msc2-mg**2)/s/msc2 )
     &                 - log( (msc2-m1**2)/mu2)
     &                 - log( (msc2-mg**2)/mu2)
     &                 + log( -s*msc2/mu2**2)         )              )
     &    + log( s/mu2  )**2 /2.D0
     &    - log( s/ms**2)**2 /2.D0
     &    + 2.D0*real( log( -sc/mu2 )
     &                 *log( -(t-msc2)/msc2) )
     &    - real( log( (msc2-m1**2)/mu2 )
     &            *log( (msc2-m1**2)/msc2) )
     &    - real( log( (msc2-mg**2)/mu2 )
     &            *log( (msc2-mg**2)/msc2) )                        )
     &  /s/(t-ms**2)
     
c                wim's integral 
      dumd_2(t,u,mg,mgc2,m1,ms,msc2,mu,zeta2) = 
     &  (   Li2( 1.D0 - (1.D0-t/mg**2)/(1.D0-m1**2/ms**2) )
     &    + Li2( 1.D0 - (mg**2-t)/(ms**2-m1**2) ) 
     &    - 2.D0 * Li2( 1.D0 - (ms**2-m1**2)/(ms**2-u) ) 
     &    + real( log( (mgc2-t)/mg/mu )**2 )
     &    + 2.D0 * real( log( (mgc2-t)/mg/mu )
     &                   *log( (msc2-u)/(ms**2-m1**2) ) )
     &    - 3.D0/4.D0 * zeta2                               )
     &  /(t-mg**2)/(u-ms**2)

c               constants 
      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0

c               define the different variables and fix the scale 
      s  = massin(1)
      m1 = massin(6) 
      mg = massin(7)
      mt = massin(9)
      ms = massin(11)
c               fix the scale to the value in the over all prefactor 
      mu2 = m1**2
      
      if (mu2.gt.0.D0) then 
         mu  = sqrt( mu2 )
      else 
         print*, " SCALAR_ARRAY_NG: problem with mu^2 ",mu2
         call HARD_STOP
      end if 
         
      t   = massin(2) + mg**2 
      u   = m1**2 + mg**2 - s - t 

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

      do n1=1,7,1
         do n2=1,5,1
            scb(n1,n2) = 9.99999D+9
         end do
      end do
      
      do n=1,4,1
         scbp(n) = 9.99999D+9 
      end do

      do n1=1,6,1
         do n2=1,4,1
            scc(n1,n2) = 9.99999D+9
         end do
      end do

      do n1=1,2,1
         do n2=1,2,1
            scd(n1,n2) = 9.99999D+9
         end do
      end do

c               the one point functions (finite parts)
      sca(1) = duma(mt,mu2)
      sca(2) = duma(ms,mu2)
      sca(3) = duma(mg,mu2)

c               the two point functions (finite parts)
      scb(1,1) = B02(t,mg,0.D0,mu2)
      scb(1,2) = B02(u,mg,0.D0,mu2)
      scb(1,3) = B02(ms**2,mg,0.D0,mu2)
      scb(1,4) = B02(mg**2,mg,0.D0,mu2)
      scb(3,1) = B02(t,ms,0.D0,mu2)
      scb(3,2) = B02(u,ms,0.D0,mu2)
      scb(3,3) = B02(ms**2,ms,0.D0,mu2)
      scb(3,4) = B02(m1**2,ms,0.D0,mu2)
      scb(3,5) = B02(mg**2,ms,0.D0,mu2)
      scb(4,1) = B02(s,ms,ms,mu2)
      scb(6,1) = B02(0.D0,mg,ms,mu2)
      scb(7,1) = B02(s,0.D0,0.D0,mu2)
      scb(8,1) = B02(mg**2,mt,ms,mu2)

c               derived two point function
      scbp(1) = BP02(0.D0,mg,ms,mu2)
      scbp(2) = BP02(mg**2,mt,ms,mu2)
      scbp(3) = BP02(mg**2,mg,0.D0,mu2)
      scbp(4) = BP02(mg**2,ms,0.D0,mu2)

c               finite three point functions 
      scc(1,1) = real( C03(m1**2,eps1,t,eps2,ms,mg) )
      scc(1,2) = real( C03(mg**2,eps1,t,eps2,ms,mg) )
      scc(1,3) = real( C03(m1**2,eps1,u,eps2,ms,mg) )
      scc(1,4) = real( C03(mg**2,eps1,u,eps2,ms,mg) )
      scc(3,1) = real( C03(eps1,eps1,s,ms,mg,ms) )
      scc(5,1) = real( C03(m1**2,mg**2,s,eps2,ms,eps2) )
      scc(6,1) = real( C03(m1**2,mg**2,s,ms,eps2,ms) )
      scc(8,2) = real( C03(mg**2,eps1,t,eps2,mg,ms) )
      scc(8,4) = real( C03(mg**2,eps1,u,eps2,mg,ms) )

c               divergent three point functions 
      msc2 = ms**2 * dcmplx(1.D0,-epsi) 
      mgc2 = mg**2 * dcmplx(1.D0,-epsi)

      scc(4,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
      scc(2,1) = dumc_1(t,abs(m1),ms,msc2,mu2)
      scc(2,2) = dumc_1(t,abs(mg),ms,msc2,mu2)
      scc(2,3) = dumc_1(u,abs(m1),ms,msc2,mu2)
      scc(2,4) = dumc_1(u,abs(mg),ms,msc2,mu2)
      scc(7,2) = dumc_2(t,mg,mgc2,mu2,zeta2)
      scc(7,4) = dumc_2(u,mg,mgc2,mu2,zeta2)

c               finite four point functions
      scd(2,1) = real( D04(eps1,eps1,m1**2,mg**2,s,t,ms,mg,ms,eps2) )
      scd(2,2) = real( D04(eps1,eps1,mg**2,m1**2,s,u,ms,mg,ms,eps2) )

c               divergent four point functions 
      sc = s * dcmplx(1.D0,epsi) 

      scd(1,1) = dumd_1(t,abs(m1),abs(mg),ms,msc2,s,sc,mu2,zeta2)
      scd(1,2) = dumd_1(u,abs(mg),abs(m1),ms,msc2,s,sc,mu2,zeta2)
      scd(3,1) = dumd_2(t,u,mg,mgc2,m1,ms,msc2,mu,zeta2) 
      scd(3,2) = dumd_2(u,t,mg,mgc2,m1,ms,msc2,mu,zeta2) 
      
      return 
      end 

