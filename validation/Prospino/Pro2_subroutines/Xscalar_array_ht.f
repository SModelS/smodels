CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C    ROUTINE TO FILL THE ARRAYS INCLUDING ALL THE SCALAR INTEGRALS    C
C     FOR THE TOP CHARGED HIGGS ASSOCIATED HADROPRODUCTION            C
C                                                                     C
C       SCALAR_ARRAY_HT_2HDM(MASSIN,SCA,SCB,SCBP,SCC,SCD)             C
C       SCALAR_ARRAY_HT_SUSY(MASSIN,SCA,SCB,SCBP,SCC,SCD)             C
C                                                                     C
C    INPUT  : MASSIN(1:50)                                            C
C                                                                     C
C             MASSIN(1)  = S                                          C
C             MASSIN(2)  = T2                                         C
C             MASSIN(6)  = M1                                         C
C             MASSIN(7)  = M2                                         C
C             MASSIN(20) = MG                                         C
C             MASSIN(21) = MS                                         C
C             MASSIN(22) = MST1                                       C
C             MASSIN(23) = MST2                                       C
C             MASSIN(24) = MSB2                                       C
C             MASSIN(25) = MSB2                                       C
C                                                                     C
C    OUTPUT : (copied from soft_virt.frm)                             C
C                                                                     C
C     A_fin(mt)   = SCA(1,1) ;
C     A_fin(mg)   = SCA(2,1) ;
C     A_fin(mst1) = SCA(3,1) ;
C     A_fin(mst2) = SCA(3,2) ;
C     A_fin(msb1) = SCA(4,1) ;
C     A_fin(msb2) = SCA(4,2) ;
C     B_fin(k1,mt,mt)      = SCB(1,1) ;
C     B_fin(k1,mg,mg)      = SCB(1,2) ;
C     B_fin(k1,ms,ms)      = SCB(1,3) ;
C     B_fin(k1,msb1,msb1)  = SCB(1,4) ;
C     B_fin(k1,msb2,msb2)  = SCB(1,5) ;
C     B_fin(k1,mst1,mst1)  = SCB(1,6) ;
C     B_fin(k1,mst2,mst2)  = SCB(1,7) ;
C     B_fin(k1,mg,msb1)    = SCB(1,8) ;
C     B_fin(k1,mg,msb2)    = SCB(1,9) ;
C     B_fin(p1,mt,0)       = SCB(2,1) ;
C     B_fin(p1,mg,mst1)    = SCB(2,2) ;
C     B_fin(p1,mg,mst2)    = SCB(2,3) ;
C     B_fin(p2,mt,0)       = SCB(3,1) ;
C     B_fin(p2,msb1,mst1)  = SCB(3,2) ;
C     B_fin(p2,msb2,mst1)  = SCB(3,3) ;
C     B_fin(p2,msb1,mst2)  = SCB(3,4) ;
C     B_fin(p2,msb2,mst2)  = SCB(3,5) ;
C     B_fin(k1+k2,0,0)     = SCB(4,1) ;
C     B_fin(k1+k2,mg,msb1) = SCB(4,2) ;
C     B_fin(k1+k2,mg,msb2) = SCB(4,3) ;
C     B_fin(p1+k2,mt,0)    = SCB(5,1) ;
C     B_fin(p1+k2,mg,msb1) = SCB(5,2) ;
C     B_fin(p1+k2,mg,msb2) = SCB(5,3) ;
C     B_fin(p1+k2,mg,mst1) = SCB(5,4) ;
C     B_fin(p1+k2,mg,mst2) = SCB(5,5) ;
C     Bp_fin(k1,mg,msb1)   = SCBP(1,1) ;
C     Bp_fin(k1,mg,msb2)   = SCBP(1,2) ;
C     Bp_fin(p1,mt,0)      = SCBP(2,1) ;
C     Bp_fin(p1,mg,mst1)   = SCBP(2,2) ;
C     Bp_fin(p1,mg,mst2)   = SCBP(2,3) ;
C     C_fin(k1,k2,0,0,0)          = SCC(1,1) ;
C     C_fin(k1,k2,mg,msb1,msb1)   = SCC(1,2) ; 
C     C_fin(k1,k2,mg,msb2,msb2)   = SCC(1,3) ; 
C     C_fin(k1,k2,msb1,mg,mg)     = SCC(1,4) ; 
C     C_fin(k1,k2,msb2,mg,mg)     = SCC(1,5) ; 
C     C_fin(p1,k1,mt,0,0)         = SCC(2,1) ;
C     C_fin(p1,k1,mst1,mg,msb1)   = SCC(2,2) ;
C     C_fin(p1,k1,mst1,mg,msb2)   = SCC(2,3) ;
C     C_fin(p1,k1,mst2,mg,msb1)   = SCC(2,4) ;
C     C_fin(p1,k1,mst2,mg,msb2)   = SCC(2,5) ;
C     C_fin(p1,k2,mt,0,0)         = SCC(3,1) ;
C     C_fin(p1,k2,0,mt,mt)        = SCC(3,2) ;
C     C_fin(p1,k2,mg,mst1,mst1)   = SCC(3,3) ;
C     C_fin(p1,k2,mg,mst2,mst2)   = SCC(3,4) ;
C     C_fin(p1,k2,mst1,mg,mg)     = SCC(3,5) ;
C     C_fin(p1,k2,mst2,mg,mg)     = SCC(3,6) ;
C     C_fin(p2,k1,mt,0,0)         = SCC(4,1) ;
C     C_fin(p2,k1,mst1,msb1,mg)   = SCC(4,2) ;
C     C_fin(p2,k1,mst1,msb2,mg)   = SCC(4,3) ;
C     C_fin(p2,k1,mst2,msb1,mg)   = SCC(4,4) ;
C     C_fin(p2,k1,mst2,msb2,mg)   = SCC(4,5) ;
C     C_fin(p2,k2,mt,0,0)         = SCC(5,1) ;
C     C_fin(p2,k2,0,mt,mt)        = SCC(5,2) ;
C     C_fin(p2,k2,msb1,mst1,mst1) = SCC(5,3) ;
C     C_fin(p2,k2,msb2,mst1,mst1) = SCC(5,4) ;
C     C_fin(p2,k2,msb1,mst2,mst2) = SCC(5,5) ;
C     C_fin(p2,k2,msb2,mst2,mst2) = SCC(5,6) ;
C     C_fin(p2,k2,mst1,msb1,msb1) = SCC(5,7) ;
C     C_fin(p2,k2,mst2,msb1,msb1) = SCC(5,8) ;
C     C_fin(p2,k2,mst1,msb2,msb2) = SCC(5,9) ;
C     C_fin(p2,k2,mst2,msb2,msb2) = SCC(5,10) ;
C     C_fin(p1,p2,0,mt,0)         = SCC(6,1) ;
C     C_fin(p1,p2,mg,mst1,msb1)   = SCC(6,2) ;
C     C_fin(p1,p2,mg,mst1,msb2)   = SCC(6,3) ;
C     C_fin(p1,p2,mg,mst2,msb1)   = SCC(6,4) ;
C     C_fin(p1,p2,mg,mst2,msb2)   = SCC(6,5) ;
C     D_fin(k1,k2,p2,0,0,0,m1)           = SCD(1,1) ;
C     D_fin(k2,k1,p2,0,0,0,m1)           = SCD(1,2) ;
C     D_fin(k2,p1,k1,m1,m1,0,0)          = SCD(1,3) ;
C     D_fin(k1,k2,p1,msb1,mg,mg,mst1)    = SCD(2,1) ;
C     D_fin(k1,k2,p1,msb2,mg,mg,mst1)    = SCD(2,2) ;
C     D_fin(k1,k2,p1,msb1,mg,mg,mst2)    = SCD(2,3) ;
C     D_fin(k1,k2,p1,msb2,mg,mg,mst2)    = SCD(2,4) ;
C     D_fin(k2,k1,p1,msb1,msb1,mg,mst1)  = SCD(3,1) ;
C     D_fin(k2,k1,p1,msb2,msb2,mg,mst1)  = SCD(3,2) ;
C     D_fin(k2,k1,p1,msb1,msb1,mg,mst2)  = SCD(3,3) ;
C     D_fin(k2,k1,p1,msb2,msb2,mg,mst2)  = SCD(3,4) ;
C     D_fin(k2,p1,k1,mst1,mst1,mg,msb1)  = SCD(4,1) ;
C     D_fin(k2,p1,k1,mst1,mst1,mg,msb2)  = SCD(4,2) ;
C     D_fin(k2,p1,k1,mst2,mst2,mg,msb1)  = SCD(4,3) ;
C     D_fin(k2,p1,k1,mst2,mst2,mg,msb2)  = SCD(4,4) ;
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
C    LAST MODIFIED : 29.03.02                                         C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_HT_2HDM(massin,sca,scb,scbp,scc,scd)
      
      implicit none 

      integer    n1,n2

      real*8     massin(1:30)
     &          ,sca(1:10,1:10),scb(1:10,1:10),scbp(1:10,1:10)
     &          ,scc(1:10,1:10),scd(1:10,1:10)
     &          ,s,t,m1,m2,mt,u,mu2
     &          ,zeta2,eps1,epsi,Li2,B02,BP02
     &          ,duma,dumc_1,dumc_3,dumd_1,dumd_3
     &          ,mg,ms,mst1,mst2,msb1,msb2

      complex*16 C03,D04,CSPEN,m1c2,m2c2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               finite part of one point function 
      duma(m1,mu2) = m1**2*( 1.D0 - log(m1**2/mu2) )
      
c               divergent three point functions, first roland's (0.34)
c               with internal mass m2 and external mass m1
      dumc_1(t,m1,m2,m2c2,mu2) =  
     &  (   Li2(t/m2**2) - Li2(m1**2/m2**2)
     &    + real( log( 1.D0 - t/m2c2     )**2 )  
     &    - real( log( 1.D0 - m1**2/m2c2 )**2 )
     &    + log(m2**2/mu2) * log( abs(-(t-m2**2)/(m2**2-m1**2) ) ) )
     &  /(t-m1**2)

c               second roland's (0.33)
      dumc_3(t,m1,m1c2,mu2,zeta2) =  
     &  (   Li2(t/m1**2)
     &    + real( log( 1.D0 - t/m1c2 )**2 )  
     &    + 1.D0/4.D0 * log(mu2/m1**2)**2 
     &    - log(mu2/m1**2)*log(abs(-(t-m1**2)/m1**2)) 
     &    + zeta2/4.D0                            )
     &  /(t-m1**2)

c               divergent four point functions, first roland's (0.46)   
c               with internal mass m1 and external masses m1,m2 
      dumd_1(t,m1,m2,m1c2,s,mu2,zeta2) =
     &  ( - 2.D0 * Li2( (m2**2-t)/(m1**2-t) )
     &    - real( log( (m1c2-m2**2)/m1/sqrt(mu2) )**2 )
     &    + 2.D0*log( (m1**2-t)/m1/sqrt(mu2) )*log(s/mu2) 
     &    - 13.D0 * zeta2 /4.D0 )
     &   /s/(t-m1**2)

c               second roland's (0.49)      
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
      s    = massin(1)
      m1   = massin(6)
      m2   = massin(7)
ctp      mg   = massin(20)
ctp      ms   = massin(21)
ctp      mst1 = massin(22)
ctp      mst2 = massin(23)
ctp      msb1 = massin(24)
ctp      msb2 = massin(25)

      mt   = m1

c               fix the scale to the value in the oval all prefactor 
      mu2 = m1**2

c               born kinematics built in 
      t    = massin(2) + m2**2 
      u    = m1**2 + m2**2 - s - t 

c               cut-off for C03 and D04, instead of 0.D0 [1.D-8]
      eps1 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               needed for the divergent three and four point functions 
      m1c2 = m1**2 * dcmplx(1.D0,-epsi) 

c               set everything to 9.99999D+9
      do n1=1,10,1
         do n2=1,10,1
            sca(n1,n2)  = 9.99999D+9 
            scb(n1,n2)  = 9.99999D+9
            scbp(n1,n2) = 9.99999D+9 
            scc(n1,n2)  = 9.99999D+9
            scd(n1,n2)  = 9.99999D+9
         end do
      end do

c               the one point functions (finite parts)
      sca(1,1) = duma(mt  ,mu2)
ctp      sca(2,1) = duma(mg  ,mu2)
ctp      sca(3,1) = duma(mst1,mu2)
ctp      sca(3,2) = duma(mst2,mu2)
ctp      sca(4,1) = duma(msb1,mu2)
ctp      sca(4,2) = duma(msb2,mu2)

c               the two point functions (finite parts)
      scb(1,1) = B02(0.D0 ,mt  ,mt   ,mu2)
ctp      scb(1,2) = B02(0.D0 ,mg  ,mg   ,mu2)
ctp      scb(1,3) = B02(0.D0 ,ms  ,ms   ,mu2)
ctp      scb(1,4) = B02(0.D0 ,msb1,msb1 ,mu2)
ctp      scb(1,5) = B02(0.D0 ,msb2,msb2 ,mu2)
ctp      scb(1,6) = B02(0.D0 ,mst1,mst1 ,mu2)
ctp      scb(1,7) = B02(0.D0 ,mst2,mst2 ,mu2)
ctp      scb(1,8) = B02(0.D0 ,mg  ,msb1 ,mu2)
ctp      scb(1,9) = B02(0.D0 ,mg  ,msb2 ,mu2)

      scb(2,1) = B02(m1**2,mt  ,0.D0 ,mu2)
ctp      scb(2,2) = B02(m1**2,mg  ,mst1 ,mu2)
ctp      scb(2,3) = B02(m1**2,mg  ,mst2 ,mu2)

      scb(3,1) = B02(m2**2,mt  ,0.D0 ,mu2)
ctp      scb(3,2) = B02(m2**2,msb1,mst1 ,mu2)
ctp      scb(3,3) = B02(m2**2,msb2,mst1 ,mu2)
ctp      scb(3,4) = B02(m2**2,msb1,mst2 ,mu2)
ctp      scb(3,5) = B02(m2**2,msb2,mst2 ,mu2)

      scb(4,1) = B02(s    ,0.D0,0.D0 ,mu2)
ctp      scb(4,2) = B02(s    ,mg  ,msb1 ,mu2)
ctp      scb(4,3) = B02(s    ,mg  ,msb2 ,mu2)

      scb(5,1) = B02(u    ,mt  ,0.D0 ,mu2)
ctp      scb(5,2) = B02(u    ,mg  ,msb1 ,mu2)
ctp      scb(5,3) = B02(u    ,mg  ,msb2 ,mu2)
ctp      scb(5,4) = B02(u    ,mg  ,mst1 ,mu2)
ctp      scb(5,5) = B02(u    ,mg  ,mst2 ,mu2)

c               derived two point function
ctp      scbp(1,1) = BP02(0.D0 ,mg  ,msb1 ,mu2)
ctp      scbp(1,2) = BP02(0.D0 ,mg  ,msb2 ,mu2)
      scbp(2,1) = BP02(m1**2,mt  ,0.D0 ,mu2)
ctp      scbp(2,2) = BP02(m1**2,mg  ,mst1 ,mu2)
ctp      scbp(2,3) = BP02(m1**2,mg  ,mst2 ,mu2)

c               finite three point functions 
ctp      scc(1,2)  = real( C03(eps1   ,eps1, s ,mg  ,msb1,msb1) )
ctp      scc(1,3)  = real( C03(eps1   ,eps1 ,s ,mg  ,msb2,msb2) )
ctp      scc(1,4)  = real( C03(eps1   ,eps1 ,s ,msb1,mg  ,mg  ) )
ctp      scc(1,5)  = real( C03(eps1   ,eps1 ,s ,msb2,mg  ,mg  ) )
      
ctp      scc(2,2)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb1) )
ctp      scc(2,3)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb2) )
ctp      scc(2,4)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb1) )
ctp      scc(2,5)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb2) )

      scc(3,2)  = real( C03(m1**2,eps1 ,u ,0.D0,mt  ,mt  ) )
ctp      scc(3,3)  = real( C03(m1**2,eps1 ,u ,mg  ,mst1,mst1) )
ctp      scc(3,4)  = real( C03(m1**2,eps1 ,u ,mg  ,mst2,mst2) )
ctp      scc(3,5)  = real( C03(m1**2,eps1 ,u ,mst1,mg  ,mg  ) )
ctp      scc(3,6)  = real( C03(m1**2,eps1 ,u ,mst2,mg  ,mg  ) )

ctp      scc(4,2)  = real( C03(m2**2,eps1 ,u ,mst1,msb1,mg  ) )
ctp      scc(4,3)  = real( C03(m2**2,eps1 ,u ,mst1,msb2,mg  ) )
ctp      scc(4,4)  = real( C03(m2**2,eps1 ,u ,mst2,msb1,mg  ) )
ctp      scc(4,5)  = real( C03(m2**2,eps1 ,u ,mst2,msb2,mg  ) )

      scc(5,2)  = real( C03(m2**2,eps1 ,t ,0.D0,mt  ,mt  ) )
ctp      scc(5,3)  = real( C03(m2**2,eps1 ,t ,msb1,mst1,mst1) )
ctp      scc(5,4)  = real( C03(m2**2,eps1 ,t ,msb2,mst1,mst1) )
ctp      scc(5,5)  = real( C03(m2**2,eps1 ,t ,msb1,mst2,mst2) )
ctp      scc(5,6)  = real( C03(m2**2,eps1 ,t ,msb2,mst2,mst2) )
ctp      scc(5,7)  = real( C03(m2**2,eps1 ,t ,mst1,msb1,msb1) )
ctp      scc(5,8)  = real( C03(m2**2,eps1 ,t ,mst2,msb1,msb1) )
ctp      scc(5,9)  = real( C03(m2**2,eps1 ,t ,mst1,msb2,msb2) )
ctp      scc(5,10) = real( C03(m2**2,eps1 ,t ,mst2,msb2,mst2) )

      scc(6,1)  = real( C03(m1**2,m2**2,s ,0.D0,mt  ,0.D0) )
ctp      scc(6,2)  = real( C03(m1**2,m2**2,s ,mg  ,mst1,msb1) )
ctp      scc(6,3)  = real( C03(m1**2,m2**2,s ,mg  ,mst1,msb2) )
ctp      scc(6,4)  = real( C03(m1**2,m2**2,s ,mg  ,mst2,msb1) )
ctp      scc(6,5)  = real( C03(m1**2,m2**2,s ,mg  ,mst2,msb2) )

c               roland's (0.20)
      scc(1,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
c                        (0.33)
      scc(2,1) = dumc_3(t,m1,m1c2,mu2,zeta2)
      scc(3,1) = dumc_3(u,m1,m1c2,mu2,zeta2)
c                        (0.34)
      scc(4,1) = dumc_1(u,m2,m1,m1c2,mu2)
      scc(5,1) = dumc_1(t,m2,m1,m1c2,mu2)

c               finite four point functions
ctp      scd(2,1) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb1,mg,mg,mst1) )
ctp      scd(2,2) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb2,mg,mg,mst1) )
ctp      scd(2,3) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb1,mg,mg,mst2) )
ctp      scd(2,4) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb2,mg,mg,mst2) )

ctp      scd(3,1) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb1,msb1,mg,mst1) )
ctp      scd(3,2) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb2,msb2,mg,mst1) )
ctp      scd(3,3) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb1,msb1,mg,mst2) )
ctp      scd(3,4) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb2,msb2,mg,mst2) )

ctp      scd(4,1) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst1,mst1,mg,msb1) )
ctp      scd(4,2) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst1,mst1,mg,msb2) )
ctp      scd(4,3) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst2,mst2,mg,msb1) )
ctp      scd(4,4) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst2,mst2,mg,msb2) )

c               divergent four point functions 
c               roland's (0.46)
      scd(1,1) = dumd_1(t,m1,m2,m1c2,s,mu2,zeta2)
      scd(1,2) = dumd_1(u,m1,m2,m1c2,s,mu2,zeta2)
c               roland's (0.49)
      scd(1,3) = dumd_3(u,t,m1,m2,mu2,zeta2)

      return 
      end 

c ---------------------------------------------------------------------
      subroutine SCALAR_ARRAY_HT_SUSY(massin,sca,scb,scbp,scc,scd)
      
      implicit none 

      integer    n1,n2

      real*8     massin(1:30)
     &          ,sca(1:10,1:10),scb(1:10,1:10),scbp(1:10,1:10)
     &          ,scc(1:10,1:10),scd(1:10,1:10)
     &          ,s,t,m1,m2,mt,u,mu2
     &          ,zeta2,eps1,epsi,Li2,B02,BP02
     &          ,duma,dumc_1,dumc_3,dumd_1,dumd_3
     &          ,mg,ms,mst1,mst2,msb1,msb2

      complex*16 C03,D04,CSPEN,m1c2,m2c2

      external B02,BP02,C03,D04,CSPEN

c               real part of the spence function included in D04
      Li2(s) = real( CSPEN(dcmplx(s)) )

c               finite part of one point function 
      duma(m1,mu2) = m1**2*( 1.D0 - log(m1**2/mu2) )
      
c               divergent three point functions, first roland's (0.34)
c               with internal mass m2 and external mass m1
      dumc_1(t,m1,m2,m2c2,mu2) =  
     &  (   Li2(t/m2**2) - Li2(m1**2/m2**2)
     &    + real( log( 1.D0 - t/m2c2     )**2 )  
     &    - real( log( 1.D0 - m1**2/m2c2 )**2 )
     &    + log(m2**2/mu2) * log( abs(-(t-m2**2)/(m2**2-m1**2) ) ) )
     &  /(t-m1**2)

c               second roland's (0.33)
      dumc_3(t,m1,m1c2,mu2,zeta2) =  
     &  (   Li2(t/m1**2)
     &    + real( log( 1.D0 - t/m1c2 )**2 )  
     &    + 1.D0/4.D0 * log(mu2/m1**2)**2 
     &    - log(mu2/m1**2)*log(abs(-(t-m1**2)/m1**2)) 
     &    + zeta2/4.D0                            )
     &  /(t-m1**2)

c               divergent four point functions, first roland's (0.46)   
c               with internal mass m1 and external masses m1,m2 
      dumd_1(t,m1,m2,m1c2,s,mu2,zeta2) =
     &  ( - 2.D0 * Li2( (m2**2-t)/(m1**2-t) )
     &    - real( log( (m1c2-m2**2)/m1/sqrt(mu2) )**2 )
     &    + 2.D0*log( (m1**2-t)/m1/sqrt(mu2) )*log(s/mu2) 
     &    - 13.D0 * zeta2 /4.D0 )
     &   /s/(t-m1**2)

c               second roland's (0.49)      
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
      s    = massin(1)
      m1   = massin(6)
      m2   = massin(7)
      mg   = massin(20)
      ms   = massin(21)
      mst1 = massin(22)
      mst2 = massin(23)
      msb1 = massin(24)
      msb2 = massin(25)

      mt   = m1

c               fix the scale to the value in the oval all prefactor 
      mu2 = m1**2

c               born kinematics built in 
      t    = massin(2) + m2**2 
      u    = m1**2 + m2**2 - s - t 

c               cut-off for C03 and D04, instead of 0.D0 [1.D-8]
      eps1 = 1.D-8

c               imaginary part for masses [1.D-8]
      epsi = 1.D-8

c               needed for the divergent three and four point functions 
      m1c2 = m1**2 * dcmplx(1.D0,-epsi) 

c               set everything to 9.99999D+9
      do n1=1,10,1
         do n2=1,10,1
            sca(n1,n2)  = 9.99999D+9 
            scb(n1,n2)  = 9.99999D+9
            scbp(n1,n2) = 9.99999D+9 
            scc(n1,n2)  = 9.99999D+9
            scd(n1,n2)  = 9.99999D+9
         end do
      end do

c               the one point functions (finite parts)
      sca(1,1) = duma(mt  ,mu2)
      sca(2,1) = duma(mg  ,mu2)
      sca(3,1) = duma(mst1,mu2)
      sca(3,2) = duma(mst2,mu2)
      sca(4,1) = duma(msb1,mu2)
      sca(4,2) = duma(msb2,mu2)

c               the two point functions (finite parts)
      scb(1,1) = B02(0.D0 ,mt  ,mt   ,mu2)
      scb(1,2) = B02(0.D0 ,mg  ,mg   ,mu2)
      scb(1,3) = B02(0.D0 ,ms  ,ms   ,mu2)
      scb(1,4) = B02(0.D0 ,msb1,msb1 ,mu2)
      scb(1,5) = B02(0.D0 ,msb2,msb2 ,mu2)
      scb(1,6) = B02(0.D0 ,mst1,mst1 ,mu2)
      scb(1,7) = B02(0.D0 ,mst2,mst2 ,mu2)
      scb(1,8) = B02(0.D0 ,mg  ,msb1 ,mu2)
      scb(1,9) = B02(0.D0 ,mg  ,msb2 ,mu2)

      scb(2,1) = B02(m1**2,mt  ,0.D0 ,mu2)
      scb(2,2) = B02(m1**2,mg  ,mst1 ,mu2)
      scb(2,3) = B02(m1**2,mg  ,mst2 ,mu2)

      scb(3,1) = B02(m2**2,mt  ,0.D0 ,mu2)
      scb(3,2) = B02(m2**2,msb1,mst1 ,mu2)
      scb(3,3) = B02(m2**2,msb2,mst1 ,mu2)
      scb(3,4) = B02(m2**2,msb1,mst2 ,mu2)
      scb(3,5) = B02(m2**2,msb2,mst2 ,mu2)

      scb(4,1) = B02(s    ,0.D0,0.D0 ,mu2)
      scb(4,2) = B02(s    ,mg  ,msb1 ,mu2)
      scb(4,3) = B02(s    ,mg  ,msb2 ,mu2)

      scb(5,1) = B02(u    ,mt  ,0.D0 ,mu2)
      scb(5,2) = B02(u    ,mg  ,msb1 ,mu2)
      scb(5,3) = B02(u    ,mg  ,msb2 ,mu2)
      scb(5,4) = B02(u    ,mg  ,mst1 ,mu2)
      scb(5,5) = B02(u    ,mg  ,mst2 ,mu2)

c               derived two point function
      scbp(1,1) = BP02(0.D0 ,mg  ,msb1 ,mu2)
      scbp(1,2) = BP02(0.D0 ,mg  ,msb2 ,mu2)
      scbp(2,1) = BP02(m1**2,mt  ,0.D0 ,mu2)
      scbp(2,2) = BP02(m1**2,mg  ,mst1 ,mu2)
      scbp(2,3) = BP02(m1**2,mg  ,mst2 ,mu2)

c               finite three point functions 
      scc(1,2)  = real( C03(eps1   ,eps1, s ,mg  ,msb1,msb1) )
      scc(1,3)  = real( C03(eps1   ,eps1 ,s ,mg  ,msb2,msb2) )
      scc(1,4)  = real( C03(eps1   ,eps1 ,s ,msb1,mg  ,mg  ) )
      scc(1,5)  = real( C03(eps1   ,eps1 ,s ,msb2,mg  ,mg  ) )
      
      scc(2,2)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb1) )
      scc(2,3)  = real( C03(m1**2,eps1 ,t ,mst1,mg  ,msb2) )
      scc(2,4)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb1) )
      scc(2,5)  = real( C03(m1**2,eps1 ,t ,mst2,mg  ,msb2) )

      scc(3,2)  = real( C03(m1**2,eps1 ,u ,0.D0,mt  ,mt  ) )
      scc(3,3)  = real( C03(m1**2,eps1 ,u ,mg  ,mst1,mst1) )
      scc(3,4)  = real( C03(m1**2,eps1 ,u ,mg  ,mst2,mst2) )
      scc(3,5)  = real( C03(m1**2,eps1 ,u ,mst1,mg  ,mg  ) )
      scc(3,6)  = real( C03(m1**2,eps1 ,u ,mst2,mg  ,mg  ) )

      scc(4,2)  = real( C03(m2**2,eps1 ,u ,mst1,msb1,mg  ) )
      scc(4,3)  = real( C03(m2**2,eps1 ,u ,mst1,msb2,mg  ) )
      scc(4,4)  = real( C03(m2**2,eps1 ,u ,mst2,msb1,mg  ) )
      scc(4,5)  = real( C03(m2**2,eps1 ,u ,mst2,msb2,mg  ) )

      scc(5,2)  = real( C03(m2**2,eps1 ,t ,0.D0,mt  ,mt  ) )
      scc(5,3)  = real( C03(m2**2,eps1 ,t ,msb1,mst1,mst1) )
      scc(5,4)  = real( C03(m2**2,eps1 ,t ,msb2,mst1,mst1) )
      scc(5,5)  = real( C03(m2**2,eps1 ,t ,msb1,mst2,mst2) )
      scc(5,6)  = real( C03(m2**2,eps1 ,t ,msb2,mst2,mst2) )
      scc(5,7)  = real( C03(m2**2,eps1 ,t ,mst1,msb1,msb1) )
      scc(5,8)  = real( C03(m2**2,eps1 ,t ,mst2,msb1,msb1) )
      scc(5,9)  = real( C03(m2**2,eps1 ,t ,mst1,msb2,msb2) )
      scc(5,10) = real( C03(m2**2,eps1 ,t ,mst2,msb2,mst2) )

      scc(6,1)  = real( C03(m1**2,m2**2,s ,0.D0,mt  ,0.D0) )
      scc(6,2)  = real( C03(m1**2,m2**2,s ,mg  ,mst1,msb1) )
      scc(6,3)  = real( C03(m1**2,m2**2,s ,mg  ,mst1,msb2) )
      scc(6,4)  = real( C03(m1**2,m2**2,s ,mg  ,mst2,msb1) )
      scc(6,5)  = real( C03(m1**2,m2**2,s ,mg  ,mst2,msb2) )

c               roland's (0.20)
      scc(1,1) = ( 1.D0/2.D0 * log(s/mu2)**2 - 7.D0/2.D0*zeta2 )/s
c                        (0.33)
      scc(2,1) = dumc_3(t,m1,m1c2,mu2,zeta2)
      scc(3,1) = dumc_3(u,m1,m1c2,mu2,zeta2)
c                        (0.34)
      scc(4,1) = dumc_1(u,m2,m1,m1c2,mu2)
      scc(5,1) = dumc_1(t,m2,m1,m1c2,mu2)

c               finite four point functions
      scd(2,1) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb1,mg,mg,mst1) )
      scd(2,2) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb2,mg,mg,mst1) )
      scd(2,3) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb1,mg,mg,mst2) )
      scd(2,4) = real( D04(eps1,eps1,m1**2,m2**2,s,u,msb2,mg,mg,mst2) )

      scd(3,1) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb1,msb1,mg,mst1))
      scd(3,2) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb2,msb2,mg,mst1))
      scd(3,3) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb1,msb1,mg,mst2))
      scd(3,4) = real( D04(eps1,eps1,m1**2,m2**2,s,t,msb2,msb2,mg,mst2))

      scd(4,1) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst1,mst1,mg,msb1))
      scd(4,2) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst1,mst1,mg,msb2))
      scd(4,3) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst2,mst2,mg,msb1))
      scd(4,4) = real( D04(eps1,m1**2,eps1,m2**2,u,t,mst2,mst2,mg,msb2))

c               divergent four point functions 
c               roland's (0.46)
      scd(1,1) = dumd_1(t,m1,m2,m1c2,s,mu2,zeta2)
      scd(1,2) = dumd_1(u,m1,m2,m1c2,s,mu2,zeta2)
c               roland's (0.49)
      scd(1,3) = dumd_3(u,t,m1,m2,mu2,zeta2)

      return 
      end 

