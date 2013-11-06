! ===========================================================================================================
module xx_integral_tb
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  public :: IFCT_TB_X12
contains
! ------------------------------
  function IFCT_TB_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 

    real(kind=double)                  :: dsig 
    integer            :: iq,inlo  
    real(kind=double)  :: m1,m2,s2l,c2l,del_s4,s,beta,mu,als_dum,strong_coup,sca_frac,s2t,c2t
    real(kind=double)  :: mx1,mx2,c2b,ms2,susy
    real(kind=double)  :: x1m,x1p,x1,x1_jac
    real(kind=double)  :: x2m,x2p,x2,x2_jac
    real(kind=double)  :: t2m,t2p,t2,t2_jac
    real(kind=double)  :: s4m,s4p,s4,s4_jac
    real(kind=double)  :: nlo,ALPHAS,LUMI
    real(kind=double)  :: TB_GGB,TB_GGV,TB_GGH,TB_GGD,TB_GG1,TB_GG2,TB_GG3
    real(kind=double)  :: TB_QBB,TB_QBV,TB_QBH,TB_QBD,TB_QB1,TB_QB2,TB_QB3
    real(kind=double)  :: TB_QGH,TB_QG3,TB_GBH,TB_GB3,BB_QBV,BB_GGV
    real(kind=double)  :: TB_GGV_REF

    if (ii>4) then                                             ! finish early for inclusive case 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    if (ipart1==1) then                                        ! note this is not the stau syntax...
       m1  =  mass_s(1)
       m2  =  mass_s(1)
       s2t =  mass_s(3)
       c2t =  mass_s(4)
       ms2 =  mass_s(2)
    else if (ipart1==2) then
       m1  =  mass_s(2)
       m2  =  mass_s(2)
       s2t = -mass_s(3)
       c2t =  mass_s(4)
       ms2 =  mass_s(1)
    end if

    mx1 = mass_x(1)                                            ! also set the sbottom parameters
    mx2 = mass_x(2)
    c2b = mass_x(4) 

    susy = 1.0                                                 ! switch for SUSY corrections

    if ( (abs(m1)+abs(m2))**2 > 0.98*sc ) then
!tp       print*, " collider energy not large enough ",m1,m2,sqrt(sc)
       dsig = 0.0
       return
    end if

    als_dum = 1.0

    x1m    = (m1+m2)**2 /sc                                    ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = (m1+m2)**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                           ! partonic cm energy

    if (isca==0) then                                          ! renormalization/factorization scale
       mu = scafac * (m1+m2)/2.0
    else if (isca==1) then
       mu = scafac * sqrt(s)
    end if

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       nlo = ALPHAS(mu,2)
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! coupling factor alpha_s not always nlo 
       inlo = 0
       strong_coup = ALPHAS(mu,1)
    else if (ii>0) then 
       inlo = 1 
       strong_coup = ALPHAS(mu,2)
    end if

    if (ii==3) then 
       del_s4 = eps_sli * m1**2                               &! unit is m^2
                   * (1.0-(m1+m2)**2/s)*(1.0-(m1-m2)**2/s)     ! rescale the s4 cutoff to be small everywhere
    else 
       del_s4 = 0.0
    end if

    select case (ii)                                           ! all the phase spaces 
    case(-1,0,1,2)                                             ! t2 integration for born, virtual 

       beta = sqrt(1.0-(m1+m2)**2/s) * sqrt(1.0-(m1-m2)**2/s)

       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

    case(3,4)                                                  ! real emission

       beta   = sqrt((1.0-del_s4/s)**2-(m1+m2)**2/s)          &
               *sqrt((1.0-del_s4/s)**2-(m1-m2)**2/s)
       t2m    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s - del_s4 + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m

       s4m    = del_s4                                         ! linear mapping just fine
       s4p    = s + t2 + m2**2 - m1**2 + s*m2**2/t2 
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m

    end select
      
    dsig = 0.0 

    select case (ii)
    case(-1,0,1)                                                                            ! born
       dsig = dsig +          LUMI(inlo,50,icoll,idub,iq,x1,x2,mu)                         &
                            * TB_GGB(s,t2,m1)
       dsig = dsig +       (  LUMI(inlo,20,icoll,idub,-1,x1,x2,mu)                         &
                            + LUMI(inlo,20,icoll,idub,+1,x1,x2,mu) )                       &
                            * TB_QBB(s,t2,m1)
    case(2)                                                                                 ! virt
       if (final_state=='tb') then
       dsig = dsig + nlo *    LUMI(inlo,50,icoll,idub,iq,x1,x2,mu)                         &
                         * (  TB_GGV(s,t2,m1,mg,mt,ms,ms2,mx1,mx2,s2t,c2t,c2b,susy)        &
                            + TB_GG1(s,t2,m1,mu) )
       dsig = dsig + nlo * (  LUMI(inlo,20,icoll,idub,-1,x1,x2,mu)                         &
                             +LUMI(inlo,20,icoll,idub,+1,x1,x2,mu) )                       &
                         * (  TB_QBV(s,t2,m1,mg,mt,ms,ms2,mx1,mx2,s2t,susy)                &
                            + TB_QB1(s,t2,m1,mu) )
       else if (final_state=='bb') then
       dsig = dsig + nlo *    LUMI(inlo,50,icoll,idub,iq,x1,x2,mu)                         &
                         * (  BB_GGV(s,t2,m1,mg,mt,ms,ms2,mx1,mx2,s2t,c2t,c2b,susy)        &
                            + TB_GG1(s,t2,m1,mu) )
       dsig = dsig + nlo * (  LUMI(inlo,20,icoll,idub,-1,x1,x2,mu)                         &
                             +LUMI(inlo,20,icoll,idub,+1,x1,x2,mu) )                       &
                         * (  BB_QBV(s,t2,m1,mg,mt,ms,ms2,mx1,mx2,susy)                    &
                            + TB_QB1(s,t2,m1,mu) )
       end if 
    case(3) 
       dsig = dsig + nlo *    LUMI(inlo,50,icoll,idub,iq,x1,x2,mu)                         &! real emission
                         * (  TB_GGH(s,t2,s4,m1)    * s4_jac                               &
                            + TB_GG3(s,t2,s4,m1,mu) * s4_jac                               &
                            + TB_GGD(s,t2,s4,m1,del_s4,s4p)                                &
                            + TB_GG2(s,t2,s4,m1,del_s4,s4p,mu)  )
       dsig = dsig + nlo * (  LUMI(inlo,20,icoll,idub,-1,x1,x2,mu)                         &
                            + LUMI(inlo,20,icoll,idub,+1,x1,x2,mu) )                       &
                         * (  TB_QBH(s,t2,s4,m1)    * s4_jac                               &
                            + TB_QB3(s,t2,s4,m1,mu) * s4_jac                               &
                            + TB_QBD(s,t2,s4,m1,del_s4,s4p)                                &
                            + TB_QB2(s,t2,s4,m1,del_s4,s4p,mu)  )
    case(4)
       dsig = dsig + nlo * (  LUMI(inlo,30,icoll,idub,-1,x1,x2,mu)                         &! qg initial state
                            + LUMI(inlo,30,icoll,idub,+1,x1,x2,mu) )                       &
                         * (  TB_QGH(s,t2,s4,m1)    * s4_jac                               &
                            + TB_QG3(s,t2,s4,m1,mu) * s4_jac )
       dsig = dsig + nlo * (  LUMI(inlo,40,icoll,idub,-1,x1,x2,mu)                         &! gb initial state
                            + LUMI(inlo,40,icoll,idub,+1,x1,x2,mu) )                       &
                         * (  TB_GBH(s,t2,s4,m1)    * s4_jac                               &
                            + TB_GB3(s,t2,s4,m1,mu) * s4_jac )
    case default
       dsig = 0.0 
    end select

    if (iscaling==0) dsig = dsig * x1_jac * x2_jac * t2_jac

    if (iscaling==0) dsig = dsig * strong_coup**2

    if (iscaling==0) dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_TB_X12
    
end module xx_integral_tb

