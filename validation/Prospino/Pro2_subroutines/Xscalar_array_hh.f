CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE CHARGED HIGGS PAIR PRODUCTION                           C
C                                                                     C
C       SCALAR_ARRAY_HH(MASSIN,SCA,SCB,SCBP,SCC,SCD)                  C
C                                                                     C
C    INPUT  : MASSIN(1:30)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = T2                                         C
C             MASSIN(6)  = MH                                         C
C             MASSIN(7)  = MT                                         C
C             MASSIN(8)  = MZ                                         C
C             MASSIN(9)  = MH1                                        C
C             MASSIN(10) = MH2                                        C
C             MASSIN(11) = MS                                         C
C                                                                     C
C    OUTPUT : SCA(1:10)                                               C
C             SCB(1:10,1:6)                                           C
C             SCBP(1:10)                                              C
C             SCC(1:20,1:9)                                           C
C             SCD(1:10,1:8)                                           C
C                                                                     C
C        Afin(mt)                      = SCA(1,1)                     C 
C        Bfin(p1,mt,0)                 = SCB(2,1)                     C 
C        Bfin(p2,mt,0)                 = SCB(3,1)                     C 
C        Bfin(k1+k2,0,0)               = SCB(4,1)                     C 
C        Bfin(p1+k2,mt,0)              = SCB(5,1)                     C 
C        Bfin(p1+k1,mt,0)              = SCB(6,1)                     C 
C        Bfin(mt,mt,0)                 = SCB(7,1)                     C 
C        Cfin(k1,k2,0,0,0)             = SCC(1,1)                     C 
C        Cfin(p1,k1,mt,0,0)            = SCC(2,1)                     C 
C        Cfin(p2,k2,mt,0,0)            = SCC(5,1)                     C 
C        Cfin(p1,p2,0,mt,0)            = SCC(6,1)                     C 
C        Dfin(k2,k1,p1,0,0,0,mt)       = SCD(1,4)                     C 
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
C    LAST MODIFIED : 02.09.04 (inset the forgotten B(mt)              C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_HH(massin,sca,scb,scbp,scc,scd)
      
      implicit none 

      integer    n,n1,n2

      real*8     massin(1:30)
     &          ,sca(1:10),scb(1:10,1:6),scbp(1:10)
     &          ,scc(1:20,1:9),scd(1:10,1:8)
     &          ,s,t,m1,m2,mt,u,mu2,ms
     &          ,zeta2,eps1,eps2,epsi,Li2,B02,BP02
     &          ,duma,dumc_1,dumd_2
     &          ,mg,msb1,msb2,mst1,mst2

      complex*16 C03,CSPEN,D04,m2c2,msc2,mtc2,sc

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
c$$$      dumc_3(t,m1,m1c2,mu2,zeta2) =  
c$$$     &  (   Li2(t/m1**2)
c$$$     &    + real( log( 1.D0 - t/m1c2**2 )**2 )  
c$$$     &    + 1.D0/4.D0 * log(mu2/m1**2)**2 
c$$$     &    - log(mu2/m1**2)*log(abs(-(t-m1**2)/m1**2)) 
c$$$     &    + zeta2/4.D0                            )
c$$$     &  /(t-m1**2)

c               divergent four point, roland's (0.46)   
c$$$      dumd_1(t,m1,m2,m1c2,s,mu2,zeta2) =
c$$$     &  ( - 2.D0 * Li2( (m2**2-t)/(m1**2-t) )
c$$$     &    - real( log( (m1c2-m2**2)/m1/sqrt(mu2) )**2 )
c$$$     &    + 2.D0*log( (m1**2-t)/m1/sqrt(mu2) )*log(s/mu2) 
c$$$     &    - 13.D0 * zeta2 /4.D0 )
c$$$     &   /s/(t-m1**2)

c               divergent four point function, roland's (0.41)
c               generalized to m1.ne.m2, copy from ng, mg->m2 
      dumd_2(t,m1,m2,ms,msc2,s,sc,mu2,zeta2) =
     &  ( - 2.D0 * Li2( 1.D0 + (ms**2-m1**2)/(t-ms**2) )
     &    - 2.D0 * Li2( 1.D0 + (ms**2-m2**2)/(t-ms**2) )
     &    - Li2( 1.D0 + (ms**2-m1**2)*(ms**2-m2**2)/s/ms**2 )
     &    - 3.D0/2.D0 * zeta2
     &    - real( log( 1.D0 + (msc2-m1**2)*(msc2-m2**2)/s/msc2 )
     &             * (   log( -(msc2-m1**2)*(msc2-m2**2)/s/msc2 )
     &                 - log( (msc2-m1**2)/mu2)
     &                 - log( (msc2-m2**2)/mu2)
     &                 + log( -s*msc2/mu2**2)         )              )
     &    + log( s/mu2  )**2 /2.D0
     &    - log( s/ms**2)**2 /2.D0
     &    + 2.D0*real( log( -sc/mu2 )
     &                 *log( -(t-msc2)/msc2) )
     &    - real( log( (msc2-m1**2)/mu2 )
     &            *log( (msc2-m1**2)/msc2) )
     &    - real( log( (msc2-m2**2)/mu2 )
     &            *log( (msc2-m2**2)/msc2) )                        )
     &  /s/(t-ms**2)
     
c               roland's (0.49)      
c$$$      dumd_3(t,u,m1,m2,mu2,zeta2) = 
c$$$     &  (  2.D0 * Li2((t-m2**2)/(m1**2-m2**2))
c$$$     &   - 2.D0 * Li2((m2**2-u)/(m1**2-u))
c$$$     &   + log((m1**2-t)/m1/sqrt(mu2))**2 
c$$$     &   + 2.D0 * log((m1**2-t)/m1/sqrt(mu2))
c$$$     &          * log(abs((m1**2-u)/(m1**2-m2**2)))
c$$$     &   - 6.D0*zeta2/8.D0                           )
c$$$     &   /(t-m1**2)/(u-m1**2)

c               constants 
      zeta2 = ( 4.D0*atan(1.D0) )**2 /6.D0

c               define the different variables and fix the scale 
      s  = massin(1)
      m1 = massin(6)
      m2 = massin(6)
      mt = massin(7)
      mg    = massin(15)
      msb1  = massin(16)
      msb2  = massin(17)
      mst1  = massin(18)
      mst2  = massin(19)

c               cut-off for external momenta in C03 and D04 [non-zero 1.D-7]
      eps1 = 1.D-7

c               cut-off for masses in C03 and D04 [1.D-8]
      eps2 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               needed for the loop particle masses
      mtc2 = mt**2 * dcmplx(1.D0,-epsi) 
      sc = s * dcmplx(1.D0,epsi) 

c               fix the scale to the value in the oval all prefactor 
      mu2 = m1**2

      t    = massin(2) + m2**2 
      u    = m1**2 + m2**2 - s - t 

c               set everything to 9.99999D+9
      do n=1,10,1
         sca(n) = 9.99999D+9 
      end do

      do n1=1,10,1
         do n2=1,6,1
            scb(n1,n2) = 9.99999D+9
         end do
      end do
      
      do n=1,10,1
         scbp(n) = 9.99999D+9 
      end do

      do n1=1,20,1
         do n2=1,9,1
            scc(n1,n2) = 9.99999D+9
         end do
      end do

      do n1=1,10,1
         do n2=1,8,1
            scd(n1,n2) = 9.99999D+9
         end do
      end do

c               the one point functions (finite parts)
      sca(1) = duma(mt,mu2)
      sca(2) = duma(mg,mu2)
      sca(3) = duma(mst1,mu2)
      sca(4) = duma(mst2,mu2)
      sca(5) = duma(msb1,mu2)
      sca(6) = duma(msb2,mu2)
      
c               the two point functions (finite parts)
      scb(1,1) = B02(0.D0 ,mg  ,msb1,mu2)
      scb(1,2) = B02(0.D0 ,mg  ,msb2,mu2)
      scb(2,1) = B02(m1**2,mt  ,0.D0,mu2)
      scb(2,2) = B02(m1**2,msb1,mst1,mu2)
      scb(2,3) = B02(m1**2,msb1,mst2,mu2)
      scb(2,4) = B02(m1**2,msb2,mst1,mu2)
      scb(2,5) = B02(m1**2,msb2,mst2,mu2)
      scb(3,1) = B02(m2**2,mt,  0.D0,mu2)
      scb(4,1) = B02(s,    0.D0,0.D0,mu2)
      scb(4,4) = B02(s,    msb1,msb1,mu2)
      scb(4,5) = B02(s,    msb1,msb2,mu2)
      scb(4,6) = B02(s,    msb2,msb2,mu2)
      scb(5,1) = B02(u,    mt,  0.D0,mu2)
      scb(6,1) = B02(t,    mt,  0.D0,mu2)
      scb(6,2) = B02(t,    mg,  mst1,mu2)
      scb(6,3) = B02(t,    mg,  mst2,mu2)
      scb(7,1) = B02(mt**2,mt,  0.D0,mu2)
      scb(7,2) = B02(mt**2,mg,  mst1,mu2)
      scb(7,3) = B02(mt**2,mg,  mst2,mu2)

      scbp(1) = BP02(0.D0,mg,msb1,mu2)
      scbp(2) = BP02(0.D0,mg,msb2,mu2)

c               three point functions 
      scc(6,1)  = real( C03(m1**2,m2**2,s ,eps2,mt  ,eps2) )

      scc(1,2)  = real( C03(eps1 ,eps1 ,s ,msb1,mg  ,msb1) )
      scc(1,3)  = real( C03(eps1 ,eps1 ,s ,msb1,mg  ,msb2) )
      scc(1,4)  = real( C03(eps1 ,eps1 ,s ,msb2,mg  ,msb1) )
      scc(1,5)  = real( C03(eps1 ,eps1 ,s ,msb2,mg  ,msb2) )

      scc(2,2)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb1) )
      scc(2,3)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb1) )
      scc(2,4)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb2) )
      scc(2,5)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb2) )

      scc(6,2)  = real( C03(m1**2,m1**2,s ,msb1,mst1,msb1) )
      scc(6,3)  = real( C03(m1**2,m1**2,s ,msb1,mst2,msb1) )
      scc(6,4)  = real( C03(m1**2,m1**2,s ,msb1,mst1,msb2) )
      scc(6,5)  = real( C03(m1**2,m1**2,s ,msb1,mst2,msb2) )
      scc(6,6)  = real( C03(m1**2,m1**2,s ,msb2,mst1,msb1) )
      scc(6,7)  = real( C03(m1**2,m1**2,s ,msb2,mst2,msb1) )
      scc(6,8)  = real( C03(m1**2,m1**2,s ,msb2,mst1,msb2) )
      scc(6,9)  = real( C03(m1**2,m1**2,s ,msb2,mst2,msb2) )

c               roland's (0.20)
      scc(1,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
c                        (0.34)
      scc(2,1) = dumc_1(t,m1,mt,mtc2,mu2)
      scc(5,1) = dumc_1(t,m1,mt,mtc2,mu2)

c               finite four point functions 
      scd(3,1) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb1,mg,msb1,mst1))
      scd(3,2) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb1,mg,msb1,mst2))
      scd(3,3) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb1,mg,msb2,mst1))
      scd(3,4) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb1,mg,msb2,mst2))
      scd(3,5) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb2,mg,msb1,mst1))
      scd(3,6) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb2,mg,msb1,mst2))
      scd(3,7) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb2,mg,msb2,mst1))
      scd(3,8) = real(D04(eps1,eps1,m1**2,m2**2,s,t,msb2,mg,msb2,mst2))

c               divergent four point functions 
      scd(1,4) = dumd_2(t,m1,m2,mt,mtc2,s,sc,mu2,zeta2)

      return 
      end 

