! ===========================================================================================================
module xx_integral_ht
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: COUPLING_HT
  public :: IFCT_HT_X12
contains
! ------------------------------
  function IFCT_HT_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 

    integer                               :: inlo
    real(kind=double), dimension(1:30)    :: massin,massin_1,massin_2
    real(kind=double), dimension(1:3)     :: lumi_ht
    real(kind=double), dimension(1:2,1:2) :: h_s
    real(kind=double)  :: m1,m2,s,beta,del_s4,h_l,h_r,mu_tgb
    real(kind=double)  :: alpha_s,nlo,ALPHAS,qr,qf,gamt,theta_s3,LUMI
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac,t2m,t2p,t2,t2_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: s4_1m,s4_1p,s4_1,s4_1_jac
    real(kind=double)  :: s3_1m,s3_1p,s3_1,s3_1_jac
    real(kind=double)  :: z3_1m,z3_1p,z3_1,z3_1_jac
    real(kind=double)  :: t2_1m,t2_1p,t2_1,t2_1_jac,beta1_1,beta2_1
    real(kind=double)  :: s4_2m,s4_2p,s4_2,s4_2_jac
    real(kind=double)  :: t2_2m,t2_2p,t2_2,t2_2_jac,s3_2,beta1_2,beta2_2,prop_s3
    real(kind=double)  :: HT_QGB,HT_QGV,HT_QGH,HT_QGD,HT_GGH,HT_QBH,HT_QQH,HT_QGS,HT_QGS1,HT_QGS2,HT_QBOS,HT_GGOS

    if (ii>9) then                                                             ! finish early 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut                  ! cut off the integration in general

    massin(1:30)   = 0.0                                                       ! initialize the massin arrays
    massin_1(1:30) = 0.0
    massin_2(1:30) = 0.0
 
    m1 = mt                                                                    ! assign the final state masses 
    m2 = mch
    mu_tgb = mu_susy * tan_b                                                   ! for the delta mb corrections
    
    x1m    = (m1+m2)**2 /sc                                                    ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = (m1+m2)**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                                           ! partonic cm energy

    if (isca==0) then                                                          ! renormalization/factorization scale
       qf = scafac  * (m1+m2)/4.0
       qr = scafac  * (m1+m2)/2.0
    else if (isca==1) then
       qf = scafac * sqrt(s)/4.0
       qr = scafac * sqrt(s)
    end if

    if (qf < 30.0) qf = 30.0                                                   ! too small scales are not good
    if (qr < 30.0) qr = 30.0

    nlo = ALPHAS(qr,2)                                                         ! different for example from hh case

    if (ii<=0) then                                                            ! coupling factor alpha_s not always nlo 
       inlo    = 0
       alpha_s = ALPHAS(qr,1)
    else if (ii>0) then 
       inlo    = 1 
       alpha_s = ALPHAS(qr,2)
    end if

    theta_s3 = 0.0                          
    gamt = ewi * m1                                                            ! needed for on-shell subtraction

    if ((s>4.0*m1**2).and.(m2<m1)) theta_s3 = 1.0 

    del_s4 = eps_sli * (m1+m2)**2                                                ! unit is m^2
!tp    s4p    = 0.0                                                               ! only for real corrections

    select case (ii)                                                           ! all the phase spaces 
    case(-1,0,1,2,3)                                                         ! born_lo, born_nlo, virt

       beta   = sqrt(1.0-(m1+m2)**2/s)                                      &  ! t2 integration
               *sqrt(1.0-(m1-m2)**2/s)      
       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )                
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

    case(4,7)                                                                  ! real(qg), crossed(qq)
       beta   = ((1.0-del_s4/s)**2-(m1+m2)**2/s)                            &  ! t2 integration shifted 
               *((1.0-del_s4/s)**2-(m1-m2)**2/s)
       if (beta<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta = sqrt( beta )
       end if

       t2m    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = del_s4                                                         ! s4 mapped to log
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2                            !  -> only a slight improvement 
       if (s4p<s4m) then
          n_faulty = n_faulty+1
          dsig = 0.0
          return
       end if
       s4     = s4m * (s4p/s4m)**var(4)
       s4_jac = s4 * log(s4p/s4m)                                              ! t2 first, since s4p explicitly needed

    case(5,6)                                                                  ! crossed(gg), crossed(qb)
       
       beta   = sqrt(1.0-(m1+m2)**2/s)                                      &  ! t2 integration shifted 
               *sqrt(1.0-(m1-m2)**2/s)
       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = 0.0
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2 
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m

       s3_1m    = 0.0                                                          ! for the s3 integration, complete phase space
       s3_1p    = s + m1**2 - m2**2 - 2.0 * sqrt(s*m1**2)
       s3_1     = var(5) * (s3_1p-s3_1m) + s3_1m
       s3_1_jac = s3_1p-s3_1m

       beta1_1    = (s-s3_1-m2**2+m1**2)**2 - 4.0*m1**2*s 
       if (beta1_1<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta1_1 = sqrt( beta1_1 )
       end if
       s4_1m    = s3_1/2.0/(s3_1+m2**2) * (s - s3_1 - m1**2 - m2**2 - beta1_1)
       s4_1p    = s3_1/2.0/(s3_1+m2**2) * (s - s3_1 - m1**2 - m2**2 + beta1_1)
       s4_1     = var(6) * (s4_1p-s4_1m) + s4_1m
       s4_1_jac = s4_1p-s4_1m
       
       beta2_1    = (s-s4_1-m1**2+m2**2)**2 - 4.0*s*m2**2  
       if (beta2_1<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta2_1 = sqrt( beta2_1 )
       end if
       t2_1m    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 + beta2_1 )
       t2_1p    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 - beta2_1 )
       t2_1     = var(7) * (t2_1p-t2_1m) + t2_1m 
       t2_1_jac = t2_1p-t2_1m 
       
       t2_1_jac = t2_1_jac * 2.0/s4_1*(s4_1+m1**2) / beta2_1                   ! including the over-all jacobian
       
    case(8,9)                                                                  ! on-shell(qb)

       if (theta_s3==0.0) then 
          dsig = 0.0
          return
       end if

       s3_1m    = 0.0                                                          ! for the s3 integration, complete phase space
       s3_1p    = s + m1**2 - m2**2 - 2.0 * sqrt(s*m1**2)
       z3_1m    = atan( (s3_1m + m2**2 - mt**2)/mt/gamt )
       z3_1p    = atan( (s3_1p + m2**2 - mt**2)/mt/gamt ) 
       z3_1     = var(3) * (z3_1p-z3_1m) + z3_1m 

       s3_1     = mt*gamt*tan(z3_1) + mt**2 - m2**2   
       s3_1_jac = ((s3_1+m2**2-mt**2)**2/mt/gamt+mt*gamt)*(z3_1p-z3_1m)

       beta1_1    = (s-s3_1-m2**2+m1**2)**2 - 4.0*m1**2*s 
       if (beta1_1<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta1_1 = sqrt( beta1_1 )
       end if
       s4_1m    = s3_1/2.0/(s3_1+m2**2) * (s - s3_1 - m1**2 - m2**2 - beta1_1)
       s4_1p    = s3_1/2.0/(s3_1+m2**2) * (s - s3_1 - m1**2 - m2**2 + beta1_1)
       s4_1     = var(4) * (s4_1p-s4_1m) + s4_1m
       s4_1_jac = s4_1p-s4_1m
       
       beta2_1    = (s-s4_1-m1**2+m2**2)**2 - 4.0*s*m2**2  
       if (beta2_1<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta2_1 = sqrt( beta2_1 )
       end if
       t2_1m    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 + beta2_1 )
       t2_1p    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 - beta2_1 )
       t2_1     = var(5) * (t2_1p-t2_1m) + t2_1m 
       t2_1_jac = t2_1p-t2_1m 
       
       t2_1_jac = t2_1_jac * 2.0/s4_1*(s4_1+m1**2) / beta2_1                   ! including the over-all jacobian
       
       s3_2 = m1**2 - m2**2                                                    ! go to constrained phase space 

       beta1_2    = (s-s3_2-m2**2+m1**2)**2 - 4.0*m1**2*s 
       if (beta1_2<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta1_2 = sqrt( beta1_2 )
       end if
       s4_2m    = s3_2/2.0/(s3_2+m2**2) * (s - s3_2 - m1**2 - m2**2 - beta1_2)
       s4_2p    = s3_2/2.0/(s3_2+m2**2) * (s - s3_2 - m1**2 - m2**2 + beta1_2)
       s4_2     = var(4) * (s4_2p-s4_2m) + s4_2m
       s4_2_jac = s4_2p-s4_2m
       
       beta2_2    = (s-s4_2-m1**2+m2**2)**2 - 4.0*s*m2**2  
       if (beta2_2<0.0) then 
          n_faulty = n_faulty+1             
          dsig = 0.0                        
          return
       else 
          beta2_2 = sqrt( beta2_2 )
       end if
       t2_2m    = -1.0/2.0 * ( s - s4_2 - m1**2 + m2**2 + beta2_2 )
       t2_2p    = -1.0/2.0 * ( s - s4_2 - m1**2 + m2**2 - beta2_2 )
       t2_2     = var(5) * (t2_2p-t2_2m) + t2_2m 
       t2_2_jac = t2_2p-t2_2m 
       
       t2_2_jac = t2_2_jac * 2.0/s4_2*(s4_2+m1**2) / beta2_2                   ! including the over-all jacobian
       
       prop_s3 = mt**2*gamt**2 / ((s3_1+m2**2-mt**2)**2+mt**2*gamt**2)         ! compensate for the wrong breit-wigner 

   end select

    massin(1)  = s                                                             ! assign the mass arrays
    massin(2)  = t2
    massin(3)  = s4                                                            ! s4 only for ii>2
    massin(6)  = m1
    massin(7)  = m2
    massin(12) = qr                                                            ! renormalization scale 
    massin(13) = qf                                                            ! factorization scale 
    massin(14) = del_s4                                                        ! slicing parameter only for ii=3
    massin(15) = s4p                                                           ! s4^max for the log(delta) terms 
    massin(20) = mg
    massin(21) = ms
    massin(22) = msq(-6)
    massin(23) = msq( 6)
    massin(24) = msq(-5)
    massin(25) = msq( 5)

    if ( (ii==5).or.(ii==6).or.(ii==8).or.(ii==9) ) then
       massin_1(1:30) = massin(1:30)
       massin_1(2)    = t2_1
       massin_1(3)    = s4_1
       massin_1(4)    = s3_1
       massin_1(8)    = gamt
    end if

    if ( (ii==8).or.(ii==9) ) then
       massin_2(1:30) = massin_1(1:30)
       massin_2(2)    = t2_2
       massin_2(3)    = s4_2
       massin_2(4)    = s3_2
    end if

    call COUPLING_HT(qr,h_l,h_r,h_s)

    lumi_ht(1:3) = 0.0                                                                              ! historic construct 
    select case (ii)
    case(-1,0,1,2,3,4)                           
       lumi_ht(1) = LUMI(inlo, 70,icoll,0, 0,x1,x2,qf)
    case(5,8)                                                                                       ! gg initial state (imode/iq=0 if no effect)
       lumi_ht(1) = LUMI(inlo, 50,icoll,0, 0,x1,x2,qf)
    case(6,9)                                                                                       ! qqbar initial state: q-qbar, q-bbar, b-bbar
       lumi_ht(1) = LUMI(inlo, 20,icoll,0,+1,x1,x2,qf)    &
                  + LUMI(inlo, 20,icoll,0,-1,x1,x2,qf)
       lumi_ht(2) = LUMI(inlo,101,icoll,0,+1,x1,x2,qf)    & 
                  + LUMI(inlo,101,icoll,0,-1,x1,x2,qf) 
       lumi_ht(3) = LUMI(inlo, 60,icoll,0, 0,x1,x2,qf) 
    case(7)                                                                                         ! qq initial state: q-b, b-b
       lumi_ht(1) = LUMI(inlo,112,icoll,0,+1,x1,x2,qf)    & 
                  + LUMI(inlo,112,icoll,0,-1,x1,x2,qf) 
       lumi_ht(2) = LUMI(inlo,111,icoll,0,+1,x1,x2,qf)    & 
                  + LUMI(inlo,111,icoll,0,-1,x1,x2,qf) 
       lumi_ht(3) = LUMI(inlo, 91,icoll,0, 0,x1,x2,qf) 
    end select

    dsig = 0.0 

    select case (ii)
    case(-1,0,1)                                                                                    ! born
       dsig =               HT_QGB(massin,h_l,h_r,lumi_ht)             * t2_jac
    case(2)                                                                                         ! virt+soft 2hdm
       dsig =       nlo *   HT_QGV(massin,h_l,h_r,lumi_ht)             * t2_jac
    case(3)                                                                                         ! susy
       dsig =               HT_QGS2(massin,h_l,h_r,lumi_ht,mu_tgb,nlo)  * t2_jac                    ! delta_b naive (ME for sign of a_b)
    case(4)                                                                                         ! real(qg) 
       dsig =       nlo * ( HT_QGH(massin,h_l,h_r,lumi_ht)                                         &
                           +HT_QGD(massin,h_l,h_r,lumi_ht) )           * t2_jac   * s4_jac 
    case(5)                                                                                         ! real(gg)
       dsig =       nlo *   HT_GGH(massin,h_l,h_r,lumi_ht)             * t2_jac   * s4_jac
       if (theta_s3==0.0)                                                                          &  
       dsig = dsig+ nlo *   HT_GGOS(massin_1,h_l,h_r,lumi_ht)          * t2_1_jac * s4_1_jac * s3_1_jac          
    case(6)                                                                                         ! real(qb)
       dsig =       nlo *   HT_QBH(massin,h_l,h_r,lumi_ht)             * t2_jac   * s4_jac
       if (theta_s3==0.0)                                                                          &  
       dsig = dsig+ nlo *   HT_QBOS(massin_1,h_l,h_r,lumi_ht)          * t2_1_jac * s4_1_jac * s3_1_jac          
    case(7)                                                                                         ! real(qq)
       dsig =       nlo *   HT_QQH(massin,h_l,h_r,lumi_ht)             * t2_jac   * s4_jac
    case(8)                                                                                         ! on-shell(gg)
       dsig =       nlo *   HT_GGOS(massin_1,h_l,h_r,lumi_ht)          * t2_1_jac * s4_1_jac * s3_1_jac &          
                  - nlo *   HT_GGOS(massin_2,h_l,h_r,lumi_ht) *prop_s3 * t2_2_jac * s4_2_jac * s3_1_jac          
    case(9)                                                                                        ! on-shell(qb)
       dsig =       nlo *   HT_QBOS(massin_1,h_l,h_r,lumi_ht)          * t2_1_jac * s4_1_jac * s3_1_jac &          
                  - nlo *   HT_QBOS(massin_2,h_l,h_r,lumi_ht) *prop_s3 * t2_2_jac * s4_2_jac * s3_1_jac          
    case default
       dsig = 0.0 
    end select
    
    dsig = dsig/s**2                                           ! remaining phase space factor
    dsig = dsig * x1_jac * x2_jac                              ! over-all jacobians 
    dsig = dsig * alpha_s                                      ! strong coupling in leading order 
    dsig = dsig * gevpb                                        ! result in pb 

    ii_done(ii) = 1

  end function IFCT_HT_X12
! ------------------------------
  subroutine COUPLING_HT(qr,h_l,h_r,h_s)
    real(kind=double),                     intent(in)  :: qr
    real(kind=double),                     intent(out) :: h_l,h_r
    real(kind=double), dimension(1:2,1:2), intent(out) :: h_s

    real(kind=double)                     :: g,sin_beta,cos_beta,sin_2beta,RUNM_EXT,yb,yt
    real(kind=double)                     :: xs2b,tag_mb,ALPHAS,B02,prefac,Cf,yukawa,zero
    real(kind=double), dimension(1:2,1:2) :: x_s

    sin_beta  = tan_b/sqrt(1.0 + tan_b**2)
    cos_beta  = sin_beta/tan_b
    sin_2beta = 2.0*sin_beta*cos_beta

    g = sqrt( 8.0*mw**2*gf /sqrt(2.0) )                                  ! on-shell scheme?

    if (ii<=0) then                                                      ! use running masses
       yt = RUNM_EXT(qr,6,1)
       yb = RUNM_EXT(qr,5,1)
    else 
       yt = RUNM_EXT(qr,6,2)
       yb = RUNM_EXT(qr,5,2)
    end if
       
    h_l = g*yt/mw /tan_b
    h_r = g*yb/mw *tan_b

    x_s(1,1) = g*mw*(sin_2beta - ( yb**2*tan_b                        &
                                  +yt**2/tan_b )/mw**2 )                 ! as defined the first argument the stop index L,R
    x_s(1,2) = g*yb/mw * (mu_susy - a_b*tan_b )
    x_s(2,1) = g*yt/mw * (mu_susy - a_t/tan_b )
    x_s(2,2) =-g*yt*yb/mw * ( tan_b + 1.0/tan_b )

    h_s = matmul( mst, x_s)
    h_s = matmul( h_s, transpose(msb) )
    
    prefac = 1.0/(16*pi**2)
    Cf     = 4.0/3.0
    zero   = 0.0
    
    xs2b   = 2.0*mu_susy*tan_b/(msq(-5)**2-msq(+5)**2)
    tag_mb = Cf*xs2b*mg*pi*ALPHAS(qr,2)*prefac                          &
             *  ( - 4*B02(zero,msq(-5:-5),mg,mt**2)                     &
                  + 4*B02(zero,msq(+5:+5),mg,mt**2) )                    ! from FORM file
    yukawa = g/sqrt(2.0)/mw * yb * tan_b/(1.0+tag_mb)
    
    if (abs(yukawa)>2.0*pi) then 
       print*, " COUPLING: Yukawa coupling exploding "
       h_l = 0.0
       h_r = 0.0
       h_s(1:2,1:2) = 0.0
    end if
    
  end subroutine COUPLING_HT

end module xx_integral_ht



