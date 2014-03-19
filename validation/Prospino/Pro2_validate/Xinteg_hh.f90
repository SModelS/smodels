! ===========================================================================================================
module xx_integral_hh
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: COUPLING_HH
  public :: IFCT_HH_X12
contains
! ------------------------------
  function IFCT_HH_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 
    integer                            :: inlo, iq
    real(kind=double), dimension(1:30) :: massin,massix,massin_1,massin_2
    real(kind=double), dimension(1:20) :: C
    real(kind=double)  :: m1,m2,s,qf,qr,beta,del_s4,gamt,beta1,beta2
    real(kind=double)  :: sw,alpha,alpha_s,mu_tgb,nlo,ALPHAS
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac,t2m,t2p,t2,t2_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: z4m,z4p,z4,s4s4,theta_s4,prop_s4,betas4,t2s4m,t2s4p,t2s4,t2s4_jac
    real(kind=double)  :: z3m,z3p,z3,s3s3,theta_s3,prop_s3,beta1s3,beta2s3,t2s3m,t2s3p,t2s3,t2s3_jac
    real(kind=double)  :: s3m,s3p,s3,s3_jac,betax,t2xm,t2xp,t2x,t2x_jac,s4s3m,s4s3p,s4s3,s4s3_jac
    real(kind=double)  :: s4xm,s4xp,s4x,s4x_jac
    real(kind=double)  :: s4_1m,s4_1p,s4_1_jac,s4_1
    real(kind=double)  :: t2_1m,t2_1p,t2_1_jac,t2_1
    real(kind=double)  :: s3_1m,s3_1p,s3_1_jac,s3_1
    real(kind=double)  :: z3_1m,z3_1p,z3_1
    real(kind=double)  :: beta1_1,beta2_1
    real(kind=double)  :: s4_2m,s4_2p,s4_2_jac,s4_2
    real(kind=double)  :: t2_2m,t2_2p,t2_2_jac,t2_2
    real(kind=double)  :: s3_2
    real(kind=double)  :: beta_1,beta1_2,beta2_2
    real(kind=double)  :: LUMI
    real(kind=double)  :: HH_QBB,HH_QBV,HH_QBSY,HH_QBS,HH_QBD,HH_QBH,HH_QGH,HH_QGOS,HH_GBH,HH_GBOS

    if (ii>8) then                                             ! finish early
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    iq = -1                                                    ! no impact for the bottom density

    massin(1:30)   = 0.0                                       ! initialize the massin arrays
    massix(1:30)   = 0.0
    massin_1(1:30) = 0.0
    massin_2(1:30) = 0.0

    m1 = mch                                                   ! assign the final state masses 
    m2 = mch
    mu_tgb = mu_susy * tan_b                                   ! for the delta mb corrections

    sw    = sqrt( 1.0 - mw**2/mz**2 )                          ! weak parameters in the on-shell scheme 
    alpha = sqrt(2.0) * mw**2 * sw**2 /pi * gf  

    x1m    = (m1+m2)**2 /sc                                    ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = (m1+m2)**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                           ! partonic cm energy

    if (iscaling==1) s = (m1+m2)**2 * (eta+1)                  ! overwrite integration for scaling fct

    if (isca==0) then                                          ! renormalization/factorization scale, factorization scale low 
       qf = scafac  * (m1+m2)/4.0
       qr = scafac  * (m1+m2)/2.0
    else if (isca==1) then
       qf = scafac * sqrt(s)/4.0
       qr = scafac * sqrt(s)
    end if

    if (qf < 30.0) qf = 30.0                                   ! too small scales are not good
    if (qr < 30.0) qr = 30.0

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       alpha_s = ALPHAS(qr,2)
       nlo     = 4.0 * pi * alpha_s
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! coupling factor alpha_s not always nlo 
       inlo = 0
    else if (ii>0) then 
       inlo = 1 
    end if

    theta_s4 = 0.0                                             ! the os subtraction theta function 
    theta_s3 = 0.0 
    gamt     = ewi * mt                                        ! always a small width needed for arcustan 

    if ((s>(mt+m2)**2).and.(mt>m1)) theta_s4 = 1.0             ! intermediate top
    if ((s>(mt+m1)**2).and.(mt>m2)) theta_s3 = 1.0             ! intermediate top

    del_s4 = eps_sli * min(m1**2,m2**2)                      & ! unit is m^2
                   * (1.0-(m1+m2)**2/s)*(1.0-(m1-m2)**2/s)     ! rescale the s4 cutoff 
    if (ii/=4) del_s4 = 0.0

    s4p    = 0.0                                                               ! only for real corrections 

    select case (ii)                                                           ! all the phase spaces 
    case(-1,0,1,2,3)                                                           ! born_lo, born_nlo, virt

       beta   = sqrt(1.0-(m1+m2)**2/s)                       &                 ! t2 integration
               *sqrt(1.0-(m1-m2)**2/s)      
       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )                
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

    case(4)                                                                    ! real(qb)
       
       beta   = ((1.0-del_s4/s)**2-(m1+m2)**2/s)              &                ! t2 integration shifted 
               *((1.0-del_s4/s)**2-(m1-m2)**2/s)                    
       if (beta>=0.0) then 
          beta = sqrt(beta)
       else 
          print*, " IFCT_HH_X12: sqrt(beta) not defined "
          beta = 0.0
       end if
       t2m    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = del_s4                                                         ! s4 mapped to log
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2                            !  -> only a slight improvement 
!       s4     = s4m * (s4p/s4m)**var(4)
!       s4_jac = s4 * log(s4p/s4m)
       if (s4p<s4m) then
          n_faulty = n_faulty+1
          dsig = 0.0
          return
       end if
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m 

    case(5)                                                                    ! crossed(qg), no os subtraction 
       
       beta   = ((1.0-del_s4/s)**2-(m1+m2)**2/s)              &                ! t2 integration shifted 
               *((1.0-del_s4/s)**2-(m1-m2)**2/s)                    
       if (beta>=0.0) then 
          beta = sqrt(beta)
       else 
          print*, " IFCT_HH_X12: sqrt(beta) not defined "
          beta = 0.0
       end if
       t2m    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = del_s4                                                         ! s4 mapped to log
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2                            !  -> only a slight improvement 
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
       
    case(6)                                                                    ! on-shell(qg)

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
       
       s3_2 = mt**2 - m2**2                                                    ! go to constrained phase space 

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

    case(7,8)                                                                  ! crossed(gb) including on-shell
                                                                       
       if ( (ii==8).and.(theta_s4==0.0) ) then 
          dsig = 0.0
          return
       end if

       s4m    = 0.0
       s4p    = s + m2**2 - m1**2 - 2.0 * sqrt(s*m2**2)

       z4m    = atan( (s4m + m1**2 - mt**2)/mt/gamt )                          ! s4 mapped to atan
       z4p    = atan( (s4p + m1**2 - mt**2)/mt/gamt )                          ! since s4t=s4+m1^2-mt^2 
       z4     = var(3) * (z4p-z4m) + z4m 

       s4     = mt*gamt*tan(z4) + mt**2 - m1**2   
       s4_jac = ((s4+m1**2-mt**2)**2/mt/gamt+mt*gamt)*(z4p-z4m) 

       beta   = sqrt( dabs( (s-s4-m1**2+m2**2)**2 - 4.0*s*m2**2 ) )            ! absolute value because of numerics
       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 + beta )                   ! t2 integration as usual
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 - beta )
       t2     = var(4) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 

       if (theta_s4==1.0) then                                                 ! restricted phase space for theta_s4=1 
          s4_1 = mt**2 - m1**2

          beta_1   = sqrt( (s-s4_1-m1**2+m2**2)**2 - 4.0*s*m2**2 )
          t2_1m    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 + beta_1 )
          t2_1p    = -1.0/2.0 * ( s - s4_1 - m1**2 + m2**2 - beta_1 )
          t2_1     = var(4) * (t2_1p-t2_1m) + t2_1m 
          t2_1_jac = t2_1p-t2_1m 

          prop_s4 = mt**2*gamt**2 / ((s4+m1**2-mt**2)**2+mt**2*gamt**2)        ! compensate for the wrong breit-wigner
       end if

    end select

    massin(1)  = s                                             ! assign the mass arrays
    massin(2)  = t2
    massin(3)  = s4                                            ! s4 only for 2->3 processes
    massin(6)  = m1
    massin(7)  = mt
    massin(8)  = mz
    massin(9)  = mh1
    massin(10) = mh2
    massin(12) = qr                                            ! renormalization scale 
    massin(13) = qf                                            ! factorization scale 
    massin(15) = mg
    massin(16) = msq(-5)
    massin(17) = msq( 5)
    massin(18) = msq(-6)
    massin(19) = msq( 6)
    massin(20) = del_s4                                        ! slicing parameter only fort ii=4
    massin(21) = s4p                                           ! s4^max for the log(delta) terms 
    massin(25) = gamt                                          ! width of the top
    massin(26) = 0.0                                           ! gamma/m of s channel particles in breit wigner

    if ( (ii==5).or.(ii==6).or.(ii==8) ) then                  ! s3 subtraction phase space
       massin_1(1:30) = massin(1:30) 
       massin_1(2)    = t2_1
       massin_1(3)    = s4_1
       massin_1(4)    = s3_1
    end if

    if (ii==6) then                                            ! restricted phace space for s3 subtraction
       massin_2(1:30) = massin(1:30)
       massin_2(2)    = t2_2
       massin_2(3)    = s4_2
       massin_2(4)    = s3_2
    end if

    call COUPLING_HH(qr,C)

    if ( (m1+m2) < mz )  C(9)  = 0.0                           ! switch off s-channel couplings  
    if ( (m1+m2) < mh1 ) C(10) = 0.0
    if ( (m1+m2) < mh2 ) C(11) = 0.0
        
    dsig = 0.0 

    select case (ii)
    case(-1,0,1)                                                                             ! born
       dsig = dsig +       LUMI(inlo,60,icoll,idub,iq,x1,x2,qf)/s**2                        &
                       *   HH_QBB(massin,C)     * t2_jac 
    case(2)                                                                                  ! virtual
       dsig = dsig + nlo * LUMI(inlo,60,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *   HH_QBV(massin,C)     * t2_jac 
    case(3)                                                                                  ! virtual-susy
       dsig = dsig +       LUMI(inlo,60,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *   HH_QBSY(massin,C,mu_tgb,alpha_s) * t2_jac 
    case(4)                                                                                  ! real(qb) 
       dsig = dsig + nlo * LUMI(inlo,60,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *(  HH_QBH(massin,C)                                                 &
                          +HH_QBD(massin,C) )   * t2_jac   * s4_jac 
    case(5)                                                                                  ! crossed(qg) 
       dsig = dsig + nlo * LUMI(inlo,70,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *   HH_QGH(massin,C)     * t2_jac   * s4_jac 
       if (theta_s3==0.0)                                                                   &  
       dsig = dsig + nlo * LUMI(inlo,70,icoll,idub,iq,x1,x2,qf)/s**2                        &  
                       *   HH_QGOS(massin_1,C)  * t2_1_jac * s4_1_jac * s3_1_jac 
    case(6)                                                                                  ! on-shell(qg)
       dsig = dsig + nlo * LUMI(inlo,70,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *(  HH_QGOS(massin_1,C)  * t2_1_jac * s4_1_jac * s3_1_jac            &
                         - HH_QGOS(massin_2,C)  * t2_2_jac * s4_2_jac * s3_1_jac * prop_s3 )
    case(7)                                                                                  ! crossed(gb) 
       dsig = dsig + nlo * LUMI(inlo,80,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                       *   HH_GBH(massin,C)     * t2_jac   * s4_jac 
       if (theta_s4==0.0)                                                                   &  
       dsig = dsig + nlo * LUMI(inlo,80,icoll,idub,iq,x1,x2,qf)/s**2                        &  
                       *   HH_GBOS(massin,C)    * t2_jac   * s4_jac 
    case(8)                                                                                  ! on-shell(qg)
       dsig = dsig + nlo * LUMI(inlo,80,icoll,idub,iq,x1,x2,qf)/s**2                        &  
                       *(  HH_GBOS(massin,C)    * t2_jac   * s4_jac                         &
                         - HH_GBOS(massin_1,C)  * t2_1_jac * s4_jac * prop_s4 ) 
    case default
       dsig = 0.0 
    end select

    if (iscaling==0)                                         & ! some jacobians
         dsig = dsig * x1_jac * x2_jac
    
    if (iscaling==0)                                         &  ! what is missing from the matrix elements 
        dsig = dsig * 4.0/(m1+m2)**2 * (4.D0*pi*alpha)**2 

    if (iscaling==0)                                         &  ! result in pb 
         dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_HH_X12
! ------------------------------
  subroutine COUPLING_HH(qr,C)
    real(kind=double),                  intent(in)   :: qr
    real(kind=double), dimension(1:20), intent(out)  :: C

    integer              :: ic 
    real(kind=double)    :: e,g 
    real(kind=double)    :: t3,qq,sw2,cw2,sw,cw,c2w,sqrt2
    real(kind=double)    :: yb,yt
    real(kind=double)    :: sa,ca,sb,cb,sbma,cbma,sbpa,cbpa,c2b
    real(kind=double)    :: RUNM_EXT

    sqrt2 = sqrt(2.0)

    t3 = -0.5                                                  ! bottom quark quantum numbers 
    qq = -1.0/3.0

    sw = sqrt( 1.0 - mw**2/mz**2 )                             ! weak mixing angle 
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)
    c2w = cw**2-sw**2

    sa = sin_a                                                 ! scalar higgs mixing angle 
    ca = cos_a

    cb   = 1.D0/dsqrt(1.D0+tan_b**2)                           ! beta as in tan(beta)
    sb   = tan_b*cb

    sbma = sb*ca - cb*sa
    cbma = cb*ca + sb*sa
    sbpa = sb*ca + cb*sa
    cbpa = cb*ca - sb*sa
    c2b  = cb**2-sb**2

    e = 1.0                                                    ! over-all factor
    g = e/sw

    if (ii<=0) then 
       yb = RUNM_EXT(qr,5,1)
       yt = RUNM_EXT(qr,6,1)
    else 
       yb = RUNM_EXT(qr,5,2)
       yt = RUNM_EXT(qr,6,2)
    end if

    C(1:20) = 0.0                                              ! initialize arrays
 
    C(1)  = e/sw/cw * ( t3 - qq*sw2 )                          ! lq 
    C(2)  = e/sw/cw * (    - qq*sw2 )                          ! rq
    C(3)  = e*qq                                               ! pq

    C(4)  =  g * yb*tan_b/mw                                   ! hl
    C(5)  =  g * yt/tan_b/mw                                   ! hr
    C(6)  =  g/sqrt2 * yb/mw * sa/cb                           ! h1 
    C(7)  = -g/sqrt2 * yb/mw * ca/cb                           ! h2

    C(8)  = -e                                                 ! ssp 
    C(9)  = -g*c2w/2.0/cw                                      ! ssz 
    C(10) = -g*( mw*sbma + mz/2.0/cw*c2b*sbpa )                ! lambda1 (checked with ME)
    C(11) = -g*( mw*cbpa - mz/2.0/cw*c2b*cbma )                ! lambda2 (checked with ME)

  end subroutine COUPLING_HH

end module xx_integral_hh



