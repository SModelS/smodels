! ===========================================================================================================
module xx_integral_ns
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: COUPLING_NS_SIM,COUPLING_NS_EXT
  public :: IFCT_NS_X12
contains
! ------------------------------
  function IFCT_NS_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 
    integer                            :: iq,ih,iq1,iq2,iqf,isq,jq1,jq2,icount1,icount2
    integer                            :: inlo,iswitch1
    real(kind=double), dimension(1:30) :: massin,massix,massix_s3,massin_s4
    real(kind=double), dimension(-6:6) :: charge
    real(kind=double)  :: m1,m2,s,qf,qr,beta,del_s4,gams,gamg,beta1,beta2
    real(kind=double)  :: sw,alpha,alpha_s,nlo,ALPHAS
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac,t2m,t2p,t2,t2_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: z4m,z4p,z4,s4s4,theta_s4,prop_s4,betas4,t2s4m,t2s4p,t2s4,t2s4_jac
    real(kind=double)  :: z3m,z3p,z3,s3s3,theta_s3,prop_s3,beta1s3,beta2s3,t2s3m,t2s3p,t2s3,t2s3_jac
    real(kind=double)  :: s3m,s3p,s3,s3_jac,t2xm,t2xp,t2x,t2x_jac,s4s3m,s4s3p,s4s3,s4s3_jac
    real(kind=double)  :: s4xm,s4xp,s4x,s4x_jac,mult,sum_charge,fs_charge,symmfac
    real(kind=double)  :: LUMI
    real(kind=double)  :: NS_QGB,NS_QGV,NS_QGH,NS_QGD,NS_GGH,NS_GGOS,NS_QQH,NS_QQOS,NS_QBH,NS_QBOS4,NS_QBOS3,NS_QBPI2
    real(kind=double), dimension(-6:6)           :: pdf1,pdf2
    real(kind=double), dimension(1:2)            :: lumi_ex
    complex(kind=double), dimension(1:4)         :: C1,C2
    complex(kind=double), dimension(1:5,1:5,1:4) :: Cx

    if (ii>1) then                                             ! finish early for inclusive case 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    massin(1:30) = 0.0                                         ! initialize the massin arrays
    massix(1:30) = 0.0

    charge(-6:6) = 99.0                                        ! charge of the quark flavors (according to cteq)
    charge( 1)   =  2.0/3.0                                    ! u quark
    charge(-1)   = -2.0/3.0                                    ! u antiquark
    charge( 2)   = -1.0/3.0                                    ! d quark
    charge(-2)   =  1.0/3.0                                    ! d antiquark
    charge( 3)   = -1.0/3.0                                    ! s quark
    charge(-3)   =  1.0/3.0                                    ! s antiquark
    charge( 4)   =  2.0/3.0                                    ! c quark
    charge(-4)   = -2.0/3.0                                    ! c antiquark
    charge( 5)   = -1.0/3.0                                    ! b quark
    charge(-5)   =  1.0/3.0                                    ! b antiquark

    select case (ipart1)
    case(1,2,3,4)
       fs_charge =   0.0
    case(5,6) 
       fs_charge = + 1.0
    case(7,8)
       fs_charge = - 1.0
    end select

    if (ii==-1) then                                           ! overwrite local squark-mass parameter
       if (i_ngtest==0) then 
          ms = msq(isquark1)
       else if (i_ngtest==1) then 
          ms = ms
       else 
          print*, " IFCT_NS_X12: i_ngtest not set "
          call HARD_STOP
       end if
    end if

    m1 = ms                                                    ! assign the final state masses 
    m2 = mass_n(ipart1)                                        !  n.b. m1=ms identified later in the code

    if ( (abs(m1)+abs(m2))**2 > 0.98*sc ) then
!tp       print*, " collider energy not large enough ",m1,m2,sqrt(sc)
       dsig = 0.0
       return
    end if

    sw    = sqrt( 1.0 - mw**2/mz**2 )                          ! weak parameters in the on-shell scheme 
    alpha = sqrt(2.0) * mw**2 * sw**2 /pi * gf  

    x1m    = (abs(m1)+abs(m2))**2 /sc                          ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = (abs(m1)+abs(m2))**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                           ! partonic cm energy

    if (iscaling==1) s = (abs(m1)+abs(m2))**2 * (eta+1)        ! overwrite integration for scaling fct

    if (isca==0) then                                          ! renormalization/factorization scale
       qf = scafac * (abs(m1)+abs(m2))/2.0
       qr = scafac * (abs(m1)+abs(m2))/2.0
    else if (isca==1) then
       qf = scafac * sqrt(s)
       qr = scafac * sqrt(s)
    end if

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       nlo = 4.0 * pi * ALPHAS(qr,2)
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! coupling factor alpha_s not always nlo 
       inlo = 0
       alpha_s = ALPHAS(qr,1)
    else if (ii>0) then 
       inlo = 1 
       alpha_s = ALPHAS(qr,2)
    end if

    theta_s4 = 0.0                                             ! the os subtraction theta function 
    theta_s3 = 0.0 
    gams     = ewi * ms                                        ! always a small width needed for arcustan 
    gamg     = ewi * mg

    if ((s>(mg+abs(m2))**2).and.(mg>    m1 )) theta_s4 = 1.0   ! intermediate gluino (gluino-neutralino associated) 
    if ((s>(ms+    m1 )**2).and.(ms>abs(m2))) theta_s3 = 1.0   ! intermediate squark (squark pair)

    del_s4 = eps_sli * min(m1**2,m2**2)                                    &    ! unit is m^2
                   * (1.0-(abs(m1)+abs(m2))**2/s)*(1.0-(abs(m1)-abs(m2))**2/s)  ! rescale the s4 cutoff 
    if (ii/=3) del_s4 = 0.0

    s4p    = 0.0                                                               ! only for real corrections needed 

    select case (ii)                                                           ! all the phase spaces 
    case(-1,0,1,2)                                                             ! born_lo, born_nlo, virt

       beta   = sqrt(1.0-(abs(m1)+abs(m2))**2/s)                            &  ! t2 integration
               *sqrt(1.0-(abs(m1)-abs(m2))**2/s)      
       t2m    = -1.0/2.0 * ( s + m2**2 - m1**2 + s*beta )                
       t2p    = -1.0/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

    case(3,4)                                                                ! real(qg), crossed(qbq), crossed(qb), no os subtraction 
       
       beta   = ((1.0-del_s4/s)**2-(abs(m1)+abs(m2))**2/s)                  &  ! t2 integration shifted 
               *((1.0-del_s4/s)**2-(abs(m1)-abs(m2))**2/s)                    
       if (beta>=0.0) then 
          beta = sqrt(beta)
       else 
          print*, " IFCT_NS_X12: sqrt(beta) not defined "
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
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m 

    case(6)                                                                    ! real(qb), real(qq) with os s4g subtraction
                                                                       
       s4m    = 0.0
       s4p    = s + m2**2 - m1**2 - 2.0 * sqrt(s*m2**2)

       z4m    = atan( (s4m + m1**2 - mg**2)/mg/gamg )                          ! s4 mapped to atan
       z4p    = atan( (s4p + m1**2 - mg**2)/mg/gamg )                          ! since s4g=s4+m1^2-mg^2 
       z4     = var(3) * (z4p-z4m) + z4m 

       s4     = mg*gamg*tan(z4) + mg**2 - m1**2   
       s4_jac = ((s4+m1**2-mg**2)**2/mg/gamg+mg*gamg)*(z4p-z4m) 

       beta   = sqrt( dabs( (s-s4-m1**2+m2**2)**2 - 4.0*s*m2**2 ) )            ! absolute value because of numerics in subtraction
       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 + beta )                   ! t2 integration as usual
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 - beta )
       t2     = var(4) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 

       if (theta_s4==1.0) then                                                 ! restricted phase space for theta_s4=1 
          s4s4 = mg**2 - m1**2

          betas4   = sqrt( (s-s4s4-m1**2+m2**2)**2 - 4.0*s*m2**2 )
          t2s4m    = -1.0/2.0 * ( s - s4s4 - m1**2 + m2**2 + betas4 )
          t2s4p    = -1.0/2.0 * ( s - s4s4 - m1**2 + m2**2 - betas4 )
          t2s4     = var(4) * (t2s4p-t2s4m) + t2s4m 
          t2s4_jac = t2s4p-t2s4m 

          prop_s4 = mg**2*gamg**2 / ((s4+m1**2-mg**2)**2+mg**2*gamg**2)        ! compensate for the wrong breit-wigner
       end if

    case(5,7)                                                                  ! real(gg), crossed(qq), crossed(qb) os subtraction s3s 
       
       s3m    = 0.0                                                            ! s3 mapped to atan
       s3p    = s + m1**2 - m2**2 - 2.0 * sqrt(s*m1**2)

       z3m    = atan( (s3m + m2**2 - ms**2)/ms/gams )
       z3p    = atan( (s3p + m2**2 - ms**2)/ms/gams ) 
       z3     = var(3) * (z3p-z3m) + z3m 

       s3     = ms*gams*tan(z3) + ms**2 - m2**2   
       s3_jac = ((s3+m2**2-ms**2)**2/ms/gams+ms*gams)*(z3p-z3m)
       
       beta1   = sqrt( (s-s3-m2**2+m1**2)**2 - 4.0*m1**2*s )
       s4xm    = s3/2.0/(s3+m2**2) * (s - s3 - m1**2 - m2**2 - beta1)
       s4xp    = s3/2.0/(s3+m2**2) * (s - s3 - m1**2 - m2**2 + beta1)
       s4x     = var(4) * (s4xp-s4xm) + s4xm
       s4x_jac = s4xp-s4xm

       beta2   = sqrt( (s-s4x-m1**2+m2**2)**2 - 4.0*m2**2*s ) 
       t2xm    = -1.0/2.0 * ( s - s4x - m1**2 + m2**2 + beta2 )                ! angular integration re-written
       t2xp    = -1.0/2.0 * ( s - s4x - m1**2 + m2**2 - beta2 )
       t2x     = var(5) * (t2xp-t2xm) + t2xm 
       t2x_jac = t2xp-t2xm 

       t2x_jac = t2x_jac * 2.0*(s4x+m1**2)/s4x/beta2                           ! the additional jacobian

       if (theta_s3==1.D0) then 
          s3s3 = ms**2 - m2**2                                                 ! restricted phase space for theta_s3=1 
         
          beta1s3   = sqrt( (s-s3s3-m2**2+m1**2)**2 - 4.0*m1**2*s )
          s4s3m    = s3s3/2.0/(s3s3+m2**2) * (s - s3s3 - m1**2 - m2**2 - beta1s3)
          s4s3p    = s3s3/2.0/(s3s3+m2**2) * (s - s3s3 - m1**2 - m2**2 + beta1s3)
          s4s3     = var(4) * (s4s3p-s4s3m) + s4s3m
          s4s3_jac = s4s3p-s4s3m
          
          beta2s3  = sqrt( (s-s4s3-m1**2+m2**2)**2 - 4.0*s*m2**2 ) 
          t2s3m    = -1.0/2.0 * ( s - s4s3 - m1**2 + m2**2 + beta2s3 )         ! angular integration re-written
          t2s3p    = -1.0/2.0 * ( s - s4s3 - m1**2 + m2**2 - beta2s3 )
          t2s3     = var(5) * (t2s3p-t2s3m) + t2s3m 
          t2s3_jac = t2s3p-t2s3m 

          t2s3_jac = t2s3_jac * 2.0*(s4s3+m1**2)/s4s3/beta2s3
          
          prop_s3 = ms**2*gams**2 / ((s3+m2**2-ms**2)**2+ms**2*gams**2)        ! compensate for the wrong breit-wigner 
       end if

    end select

    massin(1)  = s                                             ! assign the mass arrays
    massin(2)  = t2
    massin(3)  = s4                                            ! s4 only for ii=3,7,9
    massin(6)  = m1
    massin(7)  = m2
    massin(9)  = mt
    massin(10) = mg
    massin(11) = ms
    massin(12) = qr                                            ! renormalization scale 
    massin(13) = qf                                            ! factorization scale 
    massin(20) = del_s4                                        ! slicing parameter only fort ii=3
    massin(21) = s4p                                           ! s4^max for the log(delta) terms 
    massin(25) = gams                                          ! width of the squark
    massin(26) = gamg                                          ! width of the gluino

    if ( (ii==5).or.(ii==7) ) then                             ! for the s3s on-shell subtraction
       massix(1:30)    = massin(1:30)
       massix(2)       = t2x
       massix(3)       = s4x
       massix(4)       = s3
       if (theta_s3==1.D0) then 
          massix_s3(1:30) = massin(1:30)
          massix_s3(2)    = t2s3
          massix_s3(3)    = s4s3
          massix_s3(4)    = s3s3
       end if
    end if

    if (ii==6) then                                            ! for the s4g subtraction
       massin_s4(1:30) = massin(1:30)
       massin_s4(2)    = t2s4
       massin_s4(3)    = s4s4 
    end if

    if ( (ii==-1).or.(ii>5) ) then
       call GET_PDF(inlo,x1,qf,pdf1)
       call GET_PDF(inlo,x2,qf,pdf2)
    end if

    dsig = 0.0 

    if (ii<=5) then                                                                   ! as simple only for simple lumis

       do iq  = -1,+1,2                                                                ! loop for up-type and down-type quarks 
       do ih  = -1,+1,2                                                                ! loop for the squark helicities (-1=left, +1=right)

          if (iq==1) then                                                             ! multiplicity of squark flavors in the gg channel 
             mult = 2.0                              
          else if (iq==-1) then 
             mult = 3.0
          end if

          call COUPLING_NS_SIM(iq,ih,+1,C1)                                           ! charge conservation built into coupling routine (squark couplings)
          call COUPLING_NS_SIM(iq,ih,-1,C2)                                           ! antisquark couplings 

          select case (ii)
          case(-1)                                                                    ! born

             if (isquark1==0) cycle
             if ( (ih==+1).and.(isquark1<0) ) cycle                                   ! easiest: integrate isquark1 by skipping loop
             if ( (ih==-1).and.(isquark1>0) ) cycle
             if ( ipart1 < 5 ) then                                                   ! neutralinos
                iswitch1 = isquark1 
             else                                                                     ! charginos
                select case (isquark1)
                case(+1) 
                   iswitch1 = +2
                case(-1) 
                   iswitch1 = -2
                case(+2) 
                   iswitch1 = +1
                case(-2) 
                   iswitch1 = -1
                case(+3) 
                   iswitch1 = +4
                case(-3) 
                   iswitch1 = -4
                case(+4) 
                   iswitch1 = +3
                case(-4) 
                   iswitch1 = -3
                case(-5,5)
                   print*, " IFCT_NS_X12: problem with bottom case ",isquark1
                   call HARD_STOP
                end select
             end if

             if ( (iq==-1).and.((abs(iswitch1)==1).or.(abs(iswitch1)==4)) ) cycle
             if ( (iq==+1).and.((abs(iswitch1)==2).or.(abs(iswitch1)==3)) ) cycle
             
             if (icoll==0) then                                                       ! Tevatron   
                lumi_ex(1) = pdf1( abs(iswitch1))*pdf2(0) + pdf1(0)*pdf2(-abs(iswitch1))
                lumi_ex(2) = pdf1(-abs(iswitch1))*pdf2(0) + pdf1(0)*pdf2( abs(iswitch1))
             else                                                                     ! LHC
                lumi_ex(1) = pdf1( abs(iswitch1))*pdf2(0) + pdf1(0)*pdf2( abs(iswitch1))
                lumi_ex(2) = pdf1(-abs(iswitch1))*pdf2(0) + pdf1(0)*pdf2(-abs(iswitch1))
             end if

             if ( (C1(1)*C1(2))/=0.0 )                                               &! make sure this is not zero anyway
             dsig = dsig + lumi_ex(1) /s**2 * NS_QGB(massin,C1) * t2_jac              ! quark-gluon

             if ( (C2(1)*C2(2))/=0.0 )                                               &! antiquark-gluon
             dsig = dsig + lumi_ex(2) /s**2 * NS_QGB(massin,C2) * t2_jac
          case(0,1)                                                                   ! born
             if ( (C1(1)*C1(2))/=0.0 )                                               &! make sure this is not zero anyway
             dsig = dsig +       LUMI(inlo,30,icoll,idub,iq,x1,x2,qf) /s**2          &! quark-gluon
                         *    NS_QGB(massin,C1) * t2_jac 
             if ( (C2(1)*C2(2))/=0.0 )                                               &! antiquark-gluon
             dsig = dsig +       LUMI(inlo,40,icoll,idub,iq,x1,x2,qf) /s**2          &
                         *    NS_QGB(massin,C2) * t2_jac 
          case(2)                                                                     ! virt+soft (incl. scale) 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf) /s**2          & 
                         *    NS_QGV(massin,C1)  *t2_jac 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf) /s**2          & 
                         *    NS_QGV(massin,C2)  *t2_jac 
          case(3)                                                                     ! real(qg) 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf) /s**2          &
                         * (  NS_QGH(massin,C1)   *t2_jac  *s4_jac                   &
                            + NS_QGD(massin,C1)   *t2_jac  *s4_jac )
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf) /s**2          &
                         * (  NS_QGH(massin,C2)   *t2_jac  *s4_jac                   &
                            + NS_QGD(massin,C2)   *t2_jac  *s4_jac )
          case(4)                                                                     ! real(gg)
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGH(massin,C1)   *t2_jac  *s4_jac 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGH(massin,C2)   *t2_jac  *s4_jac 
          case(5)                                                                     ! on-shell(gg) 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGOS(massix,C1)   *s4x_jac *t2x_jac *s3_jac
             if (theta_s3==1.0) then 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig - nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGOS(massix_s3,C1)*s4s3_jac*t2s3_jac*s3_jac *prop_s3   
             end if  
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig + nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGOS(massix,C2)   *s4x_jac *t2x_jac *s3_jac
             if (theta_s3==1.0) then 
             if ( (C1(1)*C1(2))/=0.0 )                                               &
             dsig = dsig - nlo * LUMI(inlo,50,icoll,idub,iq,x1,x2,qf) /s**2 * mult   &
                         *    NS_GGOS(massix_s3,C2)*s4s3_jac*t2s3_jac*s3_jac *prop_s3   
             end if  
          case default
             dsig = 0.0 
          end select
       end do
       end do

    else if (ii>5) then                                             ! complicated lumis for qq and qb

       symmfac = 1.0
       icount1 = 0                                                  ! only call the angular integrals once 
       icount2 = 0
       dsig = 0.0

       do iq1 =-1,1,2                                               ! the flavors as defined by cteq
       do iq2 =-1,1,2                                               !  n.b. matrix element only needs flavors
       do iqf =-4,4,1                                               !       -> sign does not matter
       do isq = 1,4,1                                               ! always produce squark, not anti-squark
          
          if (iq1==iq2) symmfac = 0.5                               ! identical quarks in the qq initial state

          if ((iq1==0).or.(iq2==0).or.(iqf==0).or.(isq==0) ) cycle  ! that would be a gluon

          if (icoll==0) then                                        ! ppbar collider
             jq1 =  iq1 
             jq2 = -iq2 
          else if (icoll>0) then                                   ! pp collider
             jq1 =  iq1 
             jq2 =  iq2 
          end if

          sum_charge = charge(jq1)+charge(jq2)-charge(iqf)-charge(isq)-fs_charge 
          if ( abs(sum_charge)>0.01 ) cycle 
             
          do ih  = -1,1,2                                           ! squark helicity numbers  

             call COUPLING_NS_EXT(ih,Cx)

             select case (ii)
             case(6)                                                                                          ! real(qq) 
                if (iq1*iq2>0) then 
                dsig = dsig + nlo * (   pdf1(jq1)*pdf2( jq2)                                             &
                                      + pdf2(jq1)*pdf1( jq2) )/s**2 * symmfac                            &
                            * NS_QQH(massin,Cx,isq,iq1,iq2,iqf,0,icount1)   *t2_jac  *s4_jac     

                else if (iq1*iq2<0) then                  
                dsig = dsig + nlo * (   pdf1(jq1)*pdf2(-jq2)                                             &    ! real(qb)
                                      + pdf2(jq1)*pdf1(-jq2) )/s**2                                      &
                            * ( NS_QBH(massin,Cx,isq,iq1,iq2,iqf,0,icount2)                              &
                               +NS_QBPI2(massin,Cx,isq,iq1,iq2,iqf) ) *t2_jac  *s4_jac     
                dsig = dsig + nlo * (   pdf1(jq1)*pdf2(-jq2)                                             &    ! os-s4(qb)
                                      + pdf2(jq1)*pdf1(-jq2) )/s**2                                      &
                            * NS_QBOS4(massin,Cx,isq,iq1,iq2,iqf)        *t2_jac  *s4_jac     
                if (theta_s4==1.0) then                  
                dsig = dsig - nlo * (   pdf1(jq1)*pdf2(-jq2)                                             &
                                      + pdf2(jq1)*pdf1(-jq2) )/s**2                                      &
                            * NS_QBOS4(massin_s4,Cx,isq,iq1,iq2,iqf)     *t2s4_jac*s4_jac* prop_s4     
                end if 
                end if
             case(7)                                                                                          ! os-s3(qq) 
                if (iq1*iq2>0) then 
                dsig = dsig + nlo * (   pdf1(jq1)*pdf2( jq2)                                             &
                                      + pdf2(jq1)*pdf1( jq2) )/s**2 * symmfac                            &
                            * NS_QQOS(massix,Cx,isq,iq1,iq2,iqf)   *s4x_jac *t2x_jac *s3_jac
                if (theta_s3==1.0) then 
                dsig = dsig - nlo * (   pdf1(jq1)*pdf2( jq2)                                             &
                                      + pdf2(jq1)*pdf1( jq2) )/s**2 * symmfac                            &
                            * NS_QQOS(massix_s3,Cx,isq,iq1,iq2,iqf)*s4s3_jac*t2s3_jac*s3_jac * prop_s3
                end if

                else if (iq1*iq2<0) then                                                                      ! os-s3(qb)
                dsig = dsig + nlo * (   pdf1(jq1)*pdf2(-jq2)                                             &
                                      + pdf2(jq1)*pdf1(-jq2) )/s**2                                      &
                            * NS_QBOS3(massix,Cx,isq,iq1,iq2,iqf)   *s4x_jac *t2x_jac *s3_jac
                if (theta_s3==1.0) then 
                dsig = dsig - nlo * (   pdf1(jq1)*pdf2(-jq2)                                             &
                                      + pdf2(jq1)*pdf1(-jq2) )/s**2                                      &
                            * NS_QBOS3(massix_s3,Cx,isq,iq1,iq2,iqf)  *s4s3_jac*t2s3_jac*s3_jac * prop_s3
                end if
                end if
             case default
                dsig = 0.0 
             end select

          end do
       end do
       end do
       end do
       end do

       icount1 = icount1+1                                                                ! routines are only called once in iq... loops
       icount2 = icount2+1

    end if

    if (iscaling==0)                                                                   &  ! some jacobians
         dsig = dsig * x1_jac * x2_jac
    
    if (iscaling==0)                                                                   &  ! scaling function flag
         dsig = dsig * 4.0/(abs(m1)+abs(m2))**2 * 4.D0*pi*alpha * 4.D0*pi*alpha_s

    if (iscaling==0)                                                                   &  ! result in pb 
         dsig = dsig * gevpb

    if ( (ii==-1).and.(dsig==0.0) ) ii_done(-2) = 0                                       ! only needed for this channel

  end function IFCT_NS_X12


! ------------------------------
  subroutine COUPLING_NS_SIM(iq,ih,isq,C)
    integer,                              intent(in)   :: iq,ih,isq
    complex(kind=double), dimension(4),   intent(out)  :: C

    integer              :: ic 
    real(kind=double)    :: t3,qq,sw2,cw2,sw,cw
    complex(kind=double) :: zzc(4,4)

    t3(iq)  = dble(iq) /2.0                                                  ! quark quantum numbers 
    qq(iq)  = 2.0/3.0 + ( dble(iq) - 1.0 ) /2.0

    if ( (isq/=1).and.(isq/=-1) ) then                                       ! prevent garbage
       print*, " COUPLING_NS_SIM: problem with isq, return zero "
       C(1:4) = 0.0
       return
    end if

    if ( (iq/=1).and.(iq/=-1) ) then
       print*, " COUPLING_NS_SIM: problem with iq, return zero "
       C(1:4) = 0.0
       return
    end if

    sw = sqrt( 1.0 - mw**2/mz**2 )
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)

    zzc(1:4,1:4) =  conjg( zz(1:4,1:4) )                                     ! complex conjugation of zz=bq

    C(1:4) = 0.0                                                             ! initialize arrays

    ic = ipart1 - 4                                                          ! set chargino number 1 or 2  
    if (ic>2) ic = ic - 2

    select case (ipart1)                                                     ! only non-zero couplings 
    case(1,2,3,4)                                                            ! neutralino

       if (ih==-1) then 
          C(1) = ( zz(ipart1,1)*sw*(qq(iq)-t3(iq))+zz(ipart1,2)*cw*t3(iq))/(sw*cw)    ! nqs coupings [gauginos]
       else if (ih==+1) then  
          C(1) =  zzc(ipart1,1)*sw* qq(iq)                                /(sw*cw)    ! in my normalization 
       end if
       
       C(2) = C(1)                                                           ! what wim calls Cp

    case(5,6)                                                                ! chargino+

       if (ih==-1) then 
          if (isq==1) then                                                   ! squark production, LUMI(30)
             if (iq==1) then                                                 ! u-g incoming 
                C(1) = uu(ic,1)/(sqrt(2.0)*sw)
                C(2) = vv(ic,1)/(sqrt(2.0)*sw)
             else if (iq==-1) then                                           ! d-g incoming 
                C(1:2) = 0.0
             end if
          else if (isq==-1) then                                             ! antisquark production, LUMI(40)
             if (iq==-1) then                                                ! dbar-g incoming 
                C(1) = vv(ic,1)/(sqrt(2.0)*sw)
                C(2) = uu(ic,1)/(sqrt(2.0)*sw)
             else if (iq==1) then                                            ! ubar-g incoming 
                C(1:2) = 0.0
             end if
          end if
       else if (ih==+1) then
          C(1) = 0.0
          C(2) = 0.0
       end if

    case(7,8)                                                                ! chargino-

       if (ih==-1) then 
          if (isq==1) then                                                   ! squark production, LUMI(30)
             if (iq==-1) then                                                ! d-g incoming 
                C(1) = vv(ic,1)/(sqrt(2.0)*sw)
                C(2) = uu(ic,1)/(sqrt(2.0)*sw)
             else if (iq==1) then                                            ! u-g incoming 
                C(1:2) = 0.0
             end if
          else if (isq==-1) then                                             ! antisquark production, LUMI(40)
             if (iq==1) then                                                 ! ubar-g incoming 
                C(1) = uu(ic,1)/(sqrt(2.0)*sw)
                C(2) = vv(ic,1)/(sqrt(2.0)*sw)
             else if (iq==-1) then                                           ! dbar-g incoming 
                C(1:2) = 0.0
             end if
         end if
       else if (ih==+1) then
          C(1) = 0.0
          C(2) = 0.0
       end if

    end select
    
  end subroutine COUPLING_NS_SIM

! ------------------------------
!  Cx includes Cj,Cjp,Cjf,Cjpf
  subroutine COUPLING_NS_EXT(ih,Cx)
    integer,                                      intent(in)   :: ih
    complex(kind=double), dimension(1:5,1:5,1:4), intent(out)  :: Cx

    integer                            :: ic,n1,n2
    real(kind=double)                  :: sw2,cw2,sw,cw
    real(kind=double), dimension(1:5)  :: t3,qq
    complex(kind=double)               :: zzc(4,4)

    t3(1:5) = (/ +1.0,-1.0,-1.0,+1.0,-1.0 /)                                 ! quarks: u,d,s,c,b

    qq(1:5) = (/ +2.0,-1.0,-1.0,+2.0,-1.0 /)
    qq(1:5) = qq(1:5)/3.0

    sw = sqrt( 1.0 - mw**2/mz**2 )
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)

    zzc(1:4,1:4) =  conjg( zz(1:4,1:4) )                                     ! complex conjugation of zz=bw

    Cx(1:5,1:5,1:4) = 0.0                                                    ! initialize arrays

    ic = ipart1 - 4                                                          ! set chargino number 1 or 2  
    if (ic>2) ic = ic - 2

    select case (ipart1)                                                      ! only non-zero couplings 
    case(1,2,3,4)                                                            ! neutralino

       do n1=1,5,1
          do n2=1,5,1
             if (n1/=n2) Cx(n1,n2,1:4) = 0.0                                 ! only flavor diagonal couplings 
             if (ih==-1) then 
                Cx(n1,n2,1) = ( zz(ipart1,1)*sw*(qq(n1)-t3(n1))+zz(ipart1,2)*cw*t3(n1))/(sw*cw)    ! nqs coupings [gauginos]
                Cx(n1,n2,3) =  zzc(ipart1,1)*sw* qq(n1)                                /(sw*cw)    ! in my normalization 
             else if (ih==+1) then  
                Cx(n1,n2,1) =  zzc(ipart1,1)*sw* qq(n1)                                /(sw*cw)
                Cx(n1,n2,3) = ( zz(ipart1,1)*sw*(qq(n1)-t3(n1))+zz(ipart1,2)*cw*t3(n1))/(sw*cw)
             end if
           end do
       end do
       
       Cx(1:5,1:5,2) = Cx(1:5,1:5,1)                                         ! neutralinos majorana 
       Cx(1:5,1:5,4) = Cx(1:5,1:5,3)

    case(5,6)                                                                ! chargino+

       if (ih==-1) then                                                      ! only left coupling 
          Cx(1:5,1:5,1) = uu(ic,1)/(sqrt(2.0)*sw)
          Cx(1:5,1:5,2) = vv(ic,1)/(sqrt(2.0)*sw)
       end if

       if (ih==+1) then
          Cx(1:5,1:5,3) = uu(ic,1)/(sqrt(2.0)*sw)
          Cx(1:5,1:5,4) = vv(ic,1)/(sqrt(2.0)*sw)
       end if

    case(7,8)                                                                ! chargino-

       if (ih==-1) then                                                      ! only left coupling  
          Cx(1:5,1:5,1) = vv(ic,1)/(sqrt(2.0)*sw)
          Cx(1:5,1:5,2) = uu(ic,1)/(sqrt(2.0)*sw)
       end if

       if (ih==+1) then                                                      ! only left coupling  
          Cx(1:5,1:5,3) = vv(ic,1)/(sqrt(2.0)*sw)
          Cx(1:5,1:5,4) = uu(ic,1)/(sqrt(2.0)*sw)
       end if

    end select

  end subroutine COUPLING_NS_EXT

end module xx_integral_ns



