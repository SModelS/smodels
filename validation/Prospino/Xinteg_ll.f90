! ===========================================================================================================
module xx_integral_ll
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: KIN_LL,COUPLING_LL
  public :: IFCT_LL_X12
contains
! ------------------------------
  function IFCT_LL_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 

    real(kind=double)                  :: dsig 
    real(kind=double), dimension(1:30) :: massin
    real(kind=double), dimension(1:99) :: mkraemer
    real(kind=double), dimension(2,2)  :: taumix
    real(kind=double), dimension(4)    :: Cs,Csx
    complex(kind=double), dimension(4) :: Cv
    integer            :: iq,inlo  
    real(kind=double)  :: m1,m2,delta_soft,s,beta,mu
    real(kind=double)  :: x1m,x1p,x1,x1_jac
    real(kind=double)  :: x2m,x2p,x2,x2_jac
    real(kind=double)  :: t2m,t2p,t2,t2_jac
    real(kind=double)  :: s4m,s4p,s4,s4_jac
    real(kind=double)  :: c1m,c1p,c1,c1_jac
    real(kind=double)  :: axm,axp,ax,ax_jac
    real(kind=double)  :: xm,xp,x,x_jac
    real(kind=double)  :: t2xm,t2xp,t2x,t2x_jac
    real(kind=double)  :: u1,u2,t1,s3,s5,tp,up
    real(kind=double)  :: sx,betax,u1x,u2x,t1x,mx
    real(kind=double)  :: sw,alpha,nlo,ALPHAS,LUMI
    real(kind=double)  :: LL_QBB,LL_QBV,LL_QBH,LL_QBF1,LL_QBF2,LL_QGH,LL_QGF
!tp needed for the scales check
!tp    real(kind=double)  :: shat,that,uhat

    if (ii>6) then                                             ! finish early for inclusive case 
       dsig = 0.0 
       return
    end if

    if ( (ii==3).or.(ii==5) ) then                             ! DY hard emission m-e included in massfac convolution
       dsig = 0.0
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    massin(1:30) = 0.0                                         ! initialize the massin array

    m1  = mass_s(1)                                            ! assign the final state masses 
    m2  = mass_s(2)
    taumix(1:2,1:2) = msl(1:2,1:2)

    if ( (abs(m1)+abs(m2))**2 > 0.98*sc ) then
!tp       print*, " collider energy not large enough ",m1,m2,sqrt(sc)
       dsig = 0.0
       return
    end if

    delta_soft = eps_sub * (m1+m2)**2/4.0            ! the soft cut-off rescaled 
    
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
    if (iscaling==1) s = (m1+mg)**2 * (eta+1)                  ! overwrite integration for scaling fct

    if (isca==0) then                                          ! renormalization/factorization scale
       mu = scafac * (m1+m2)/2.0
    else if (isca==1) then
       mu = scafac * sqrt(s)
    end if

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       nlo = 4.0 * pi * ALPHAS(mu,2)
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! coupling factor alpha_s not always nlo 
       inlo = 0
    else if (ii>0) then 
       inlo = 1 
    end if

    select case (ii)                                           ! all the phase spaces 
    case(-1,0,1,2)                                             ! t2 integration for born, virtual 

       beta = sqrt(1.0-(m1+m2)**2/s) * sqrt(1.0-(m1-m2)**2/s)

       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       u1  = - s - t2                                          ! born kinematic variables
       t1  = t2 + m2**2 - m1**2 
       u2  = u1 + m1**2 - m2**2 

!tp to check the isajet scales 2.*SHAT*THAT*UHAT/(SHAT**2+THAT**2+UHAT**2)
!tp       shat = s
!tp       that = t2 + m2**2
!tp       uhat = u2 + m2**2
!tp       mu = sqrt( 2.0*shat*that*uhat/(shat**2+that**2+uhat**2) )
!tp       mu = max( mu, m1+m2 )

    case(3,5)                                                  ! t2-s4-omega integration for real(qb)

       print*, " IFCT_LL_X12: this is only for testing! ",ii
       call HARD_STOP

       beta = sqrt(1.0-(m1+m2)**2/s) * sqrt(1.0-(m1-m2)**2/s)

       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = 0.0
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m

       c1m    = -1.0
       c1p    = 1.0
       c1     = var(5) * (c1p-c1m) + c1m 
       c1_jac = c1p-c1m 

       axm    = 0.0 
       axp    = pi 
       ax     = var(6) * (axp-axm) + axm 
       ax_jac = axp-axm

       call KIN_LL(m1,m2,s,t2,s4,c1,ax,s3,s5,u2,t1,u1,tp,up)   ! the real phase space 

       if ((abs(tp)<delta_soft).or.(abs(up)<delta_soft).or.(abs(s3)<delta_soft).or.(abs(s4)<delta_soft)) then 
          dsig = 0.0                                           ! cut the soft phase space 
          return 
       end if

    case(4,6)                                                               ! t2-x integration for mass factorization

       beta = sqrt(1.0-(m1+m2)**2/s) * sqrt(1.0-(m1-m2)**2/s)

       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       u1 = - s - t2  
       t1 = t2 + m2**2 - m1**2 
       u2 = u1 + m1**2 - m2**2 

       xm    = (m1+m2)**2 /s
       xp    = 1.0
       x     = var(4) * (xp-xm) + xm  
       x_jac = xp-xm 

       sx = s * x                                              ! shift in the variables s

       betax = sqrt(1.0-(m1+m2)**2/sx) * sqrt(1.0-(m1-m2)**2/sx)

       t2xm    = -1.0/2.0 * ( sx + m2**2 - m1**2 + sx*betax )  ! redo t2 integration 
       t2xp    = -1.0/2.0 * ( sx + m2**2 - m1**2 - sx*betax )
       t2x     = var(3) * (t2xp-t2xm) + t2xm
       t2x_jac = t2xp-t2xm 

       u1x = - sx - t2x                                        ! redo born type kinematics 
       t1x = t2x + m2**2 - m1**2 
       u2x = u1x + m1**2 - m2**2

    end select
      
    massin(1)  = s                                             ! assign the mass array
    massin(2)  = t2
    massin(3)  = u2 
    massin(4)  = t1 
    massin(5)  = u1 
    massin(6)  = m1
    massin(7)  = m2
    massin(9)  = mt
    massin(10) = mg
    massin(11) = ms
    massin(12) = mu
    if ((ii==4).or.(ii==6)) then                               ! additional entries for mass factorization  
       massin(13) = sx
       massin(14) = t2x 
       massin(15) = u2x 
       massin(16) = t1x 
       massin(17) = u1x 
       massin(18) = x
       massin(19) = xp-xm
       massin(20) = xm 
    end if
    massin(26) = 1.0                                           ! yes/no for s4 o-s subtraction
    massin(27) = 1.0                                           ! yes/no for s3 o-s subtraction 

    mkraemer(1:99) = 0.0                                       ! the array mkraemer for real correction 
    mkraemer(1)  = m1 
    mkraemer(2)  = m2 
    mkraemer(3)  = mg 
    mkraemer(4)  = ms 

    dsig = 0.0 
    do iq =-1,1,2                                              ! the loop for up-type and down-type quarks 

       call COUPLING_LL(s,m1,m2,iq,Cs,Cv,mx,taumix,mkraemer)
       if ((ii.eq.4).or.(ii.eq.6)) call COUPLING_LL(sx,m1,m2,iq,Csx,Cv,mx,taumix,mkraemer)

       massin(8) = mx
       
       select case (ii)
       case(-1,0,1)                                                                            ! born
          dsig = dsig +          LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * LL_QBB(massin,Cs)/s**2       * t2_jac
       case(2)                                                                                 ! virt
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * LL_QBV(massin,Cs,Cv)/s**2    * t2_jac
       case(3)                                                                                 ! matrix qb
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * LL_QBH(massin,mkraemer) /s**2                                   &
                            * s4 / (s4+m1**2) * t2_jac * s4_jac * c1_jac * ax_jac
       case(4)                                                                                 ! massfac qb
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                          * ( LL_QBF1(massin,Csx)/sx**2  * x_jac * t2x_jac                    &
                            - LL_QBF2(massin,Cs) /s**2   * x_jac * t2_jac  )
       case(5)                                                                                 ! matrix qg+gb
          dsig = dsig + nlo * (  LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                         &
                               + LUMI(inlo,40,icoll,idub,iq,x1,x2,mu) )                       &
                            * LL_QGH(massin,mkraemer) /s**2                                   &
                            * s4 / (s4+m1**2) * t2_jac * s4_jac * c1_jac * ax_jac
       case(6)                                                                                 ! massfac qg+gb, note that x fixes lumi
          dsig = dsig + nlo * (  LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                         &
                               + LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)  )                      &
                            * LL_QGF(massin,Csx)/sx**2   * x_jac * t2x_jac
       case default
          dsig = 0.0 
       end select
    end do
       
    if ( (ii==3).or.(ii==5) ) then                                    ! phase space without 1/s**2 
       dsig = dsig / ( 2.0 * 256.0 * pi**4 )
    else
       dsig = dsig / ( 16.0 * pi )
    end if
    
    if (iscaling==0) dsig = dsig * x1_jac * x2_jac

    if (iscaling==0) dsig = dsig * 4.0/(m1+m2)**2 * alpha**2 

    if (iscaling==0) dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_LL_X12
    
! ------------------------------
!
  subroutine KIN_LL(m1,m2,s,t2,s4,c1_in,ax,s3,s5,u2,t1,u1,tp,up)
    real(kind=double), intent(in)     :: m1,m2,s,t2,s4,c1_in,ax
    real(kind=double), intent(out)    :: s5,u2,t1,u1,tp,up
    real(kind=double), intent(inout)  :: s3 
    integer                           :: type
    real(kind=double)                 :: c1,s1,c2,norm,w1,w2,w3,e1,e2,p,cx,sx
    real(kind=double), dimension(0:5) :: test

    type = 0 
    
    norm = 2.0 * sqrt(s4 + m1**2) 
 
    u2 = s4 - s - t2 - m2**2 + m1**2 
      
    w1 = ( s + u2 )                       /norm        ! remember: same for coordinate systems A and C
    w2 = ( s + t2 )                       /norm
    w3 =  s4                              /norm
    e1 = ( s4 + 2.0*m1**2 )               /norm
    e2 = - ( t2 + u2 + 2.0*m2**2 )        /norm
    p  = sqrt( (t2+u2)**2 - 4.0*m2**2*s ) /norm
    cx = ( t2 * ( s4+m1**2-m2**2 ) - s*( u2+2.0*m2**2 ) )/ ( s+t2 ) / sqrt( (t2+u2)**2 - 4.0*m2**2*s )

    sx = sqrt( 1.0 - cx**2 )
      
    if (type==0) then 
       c1 = c1_in 
    else if (type==1) then 
       c1 = ( 2.0*w3*e2 - s3 )/(2.0*p*w3)
    else 
       print*, " KIN_LL: wrong input on type:",type
    end if

    s1 = sqrt( 1.0 - c1**2 )
    c2 = cos(ax)

    if (type==0) then 
       t1 = 2.0 * ( - w1*e1 - w3*s1*c2*p*sx - w3*p*cx*c1 + w3*w2*c1 )
       u1 = 2.0 * ( - w2*e1 - w2*w3*c1 ) 
       tp = 2.0 * ( - w2*w3 + w2*w3*c1 )
       up = 2.0 * ( - w1*w3 + w3*p*sx*s1*c2 + w3*p*cx*c1 - w2*w3*c1 )
       s3 = 2.0 * ( e2*w3  - w3*p*sx*s1*c2 - w3*p*cx*c1 )
       s5 = 2.0 * ( e1*e2 + w3*p*sx*s1*c2 + w3*p*c1*cx )+m1**2+m2**2 
    else if (type==1) then 
       t1 = 2.0 * ( - w1*e1 + w2*w3*s1*c2*sx + w2*w3*cx*c1 - w3*p*c1 )
       u1 = 2.0 * ( - w2*e1 - w2*w3*s1*c2*sx - w2*w3*c1*cx ) 
       tp = 2.0 * ( - w2*w3 + w2*w3*s1*c2*sx + w2*w3*c1*cx )
       up = 2.0 * ( - w1*w3 - w2*w3*sx*s1*c2 - w2*w3*cx*c1 + w3*p*c1 )
       s3 = 2.0 * (   w3*e2 - p*w3*c1 )
       s5 = 2.0 * ( e1*e2 + w3*p*c1 )+m1**2+m2**2 
    end if

    test(0) = s + t2 + u2 - m1**2 + m2**2 - s4     ! some checks of kinematic relations 
    test(1) = s5 + t2 + u2 + m2**2 - m1**2 + s3 
    test(2) = s5 + t1 + u1 + m1**2 - m2**2 + s4 
    test(3) = s + t1 + u2 + up 
    test(4) = s + t2 + u1 + tp
    test(5) = s + tp + up - s5

    if (test(0)>1.e-6) print *," KIN_LL: test0 better be zero "
    if (test(1)>1.e-6) print *," KIN_LL: test1 better be zero "
    if (test(2)>1.e-6) print *," KIN_LL: test2 better be zero "
    if (test(3)>1.e-6) print *," KIN_LL: test3 better be zero "
    if (test(4)>1.e-6) print *," KIN_LL: test4 better be zero "
    if (test(5)>1.e-6) print *," KIN_LL: test5 better be zero "

  end subroutine KIN_LL

! ------------------------------
! all general couplings : C(1) located at ( k1(q_in) , ipart1 )
!                         C(2) located at ( k2(q_out), ipart1 )
!                         C(3) located at ( k2(q_in) , ipart2 )
!                         C(4) located at ( k2(q_out), ipart2 )  for all Cl,Cr
!  mx   = mw,mz s channel gauge boson mass 
! ------------------------------
  subroutine COUPLING_LL(s,m1,m2,iq,Cs,Cv,mx,taumix,mkraemer)
    integer,                              intent(in)   :: iq
    real(kind=double),                    intent(in)   :: s,m1,m2
    real(kind=double), dimension(2,2),    intent(in)   :: taumix
    real(kind=double),                    intent(out)  :: mx
    real(kind=double), dimension(4),      intent(out)  :: Cs
    real(kind=double), dimension(99),     intent(out)  :: mkraemer
    complex(kind=double), dimension(4),   intent(out)  :: Cv

    real(kind=double)    :: t3,qq,t3w,qqw,v1,v1w,v2,a2
    real(kind=double)    :: sw,sw2,cw2,cw
    real(kind=double)    :: only_u,only_d,gamfac
    complex(kind=double) :: v2w,a2w

    t3(iq)  = dble(iq) /2.0                                                  ! quark quantum numbers 
    qq(iq)  = 2.0/3.0 + ( dble(iq) - 1.0 ) /2.0

    sw = sqrt( 1.0 - mw**2/mz**2 )
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)

    if ((ipart1<=3).or.((ipart1>=6).and.(ipart1<=9)).or.(ipart1==14) ) then ! neutral final state

       if ( (ipart1==3).or.(ipart1==9) ) then                               ! sneutrinos 
          t3w =  1.0/2.0 
          qqw =  0.0 
       else if (ipart1==8) then                                             ! mixed pair 
          t3w = -1.0/2.0
          qqw =  0.0 
       else if (ipart1==14) then                                            ! charged Higgs
          t3w =  0.0                                                        ! only dummy
          qqw =  1.0
       else                                                                 ! sleptons (right handed included)
          t3w = -1.0/2.0
          qqw = -1.0 
       end if

       mx = mz

       v1  = - qq(iq)                                                        ! pqq coupling 
       v1w = - qqw                                                           ! pll coupling
       v2  = - ( t3(iq) - 2.0 * sw2 * qq(iq)   )/(2.0*sw*cw)                 ! zqq couplings 
       a2  = -   t3(iq)                         /(2.0*sw*cw) 

       gamfac = 1.0
       select case(ipart1)
       case(0)
          v2w  = - ( t3w          - 2.0*sw2 * qqw )/(2.0*sw*cw) 
          a2w  = -   t3w                           /(2.0*sw*cw) 
          gamfac = 2.0 
       case(1,3,9)
          v2w  = - ( t3w         -      sw2 * qqw )/(2.0*sw*cw) 
          a2w  = v2w
       case(2)
          v2w  = - ( 0.0         -      sw2 * qqw )/(2.0*sw*cw) 
          a2w  = - v2w
       case(6)                                                               ! stau1-stau1*
          v2w  = - ( t3w*taumix(1,1)**2 - sw2 * qqw )/(2.0*sw*cw) 
          a2w  = v2w
       case(7)                                                               ! stau2-stau2*
          v2w  = - ( t3w*taumix(2,1)**2 - sw2 * qqw )/(2.0*sw*cw) 
          a2w  = v2w
       case(8)                                                               ! stau1-stau2* plus stau2-stau1*
          v2w = - t3w*taumix(1,1)*taumix(2,1) /(2.0*sw*cw)
          a2w  = v2w
       case(14)
          v2w  = -(cw2-sw2)/2.0                    /(2.0*sw*cw)
          a2w  = 0.0
       end select

    else if ( (ipart1==4).or.(ipart1==5).or.(ipart1>=10) ) then              ! charged final state

       mx = mw

       v1  = 0.0                                                                ! pqq coupling
       v1w = 0.0                                                                ! pcn coupling 
       v2  = - 1.0 /(2.0*sqrt(2.0)*sw)                                          ! wqq coupling
       a2  = - 1.0 /(2.0*sqrt(2.0)*sw) 

       gamfac = 1.0
       select case(ipart1)
       case(4,5)
          v2w  = - 1.0 /(2.0*sqrt(2.0)*sw)
       case(10,11)                                                               ! sntau-stau1
          v2w  = - taumix(1,1) * 1.0 /(2.0*sqrt(2.0)*sw)
       case(12,13)                                                               ! sntau-stau2
          v2w  = - taumix(2,1) * 1.0 /(2.0*sqrt(2.0)*sw)
       end select
       a2w = v2w

    else
       print*, " COUPLING_LL: ipart1 not set correctly ",ipart1
       call HARD_STOP
    end if

    if ( (abs(m1)+abs(m2)) < mz ) then
       print*, " COUPLINGS_LL: masses low, W/Z decays might be more suitable, decouple W/Z "
       v2   = 0.0
       a2   = 0.0
       v2w  = 0.0
       a2w  = 0.0
    end if

    Cs(1:4) = 0.0 
    Cv(1:4) = 0.0

    Cs(1) =  gamfac * (v1*v1w)**2 /2.0 + s * 2.0 * real( v1*v1w*v2*v2w )/(s-mx**2)  &  ! wim's conventions 
           + s**2*( v2**2+a2**2 )*( abs(v2w)**2+abs(a2w)**2 ) /(s-mx**2)**2             ! [v1,v2,v1w,a2 real]

    mkraemer(11) =  mx
    mkraemer(12) =  v1
    mkraemer(13) =  v2
    mkraemer(14) = -a2
    mkraemer(15) =  v1w
    mkraemer(16) =  real(v2w)
    mkraemer(17) = -real(a2w)
    mkraemer(21:28) = 0.0

    only_u = ( 1.0 + iq )/2.0                                                           ! quarks in the t/u channel 
    only_d = ( 1.0 - iq )/2.0
    
    select case (ipart1)
    case(4,10,12) 
       mkraemer(12:28) = only_u * mkraemer(12:28) 
       Cs(1)           = only_u * Cs(1)
    case(5,11,13)
       mkraemer(12:28) = only_d * mkraemer(12:28) 
       Cs(1)           = only_d * Cs(1)
    end select
    
    mkraemer(12:28) = 2.0 * sqrt(pi) * mkraemer(12:28)                                  ! rescaling with 2 and pi 
    
    Cs(1) = 16.0 * pi**2    * Cs(1)                                                     ! rescaling g^2 -> alpha

  end subroutine COUPLING_LL

end module xx_integral_ll

