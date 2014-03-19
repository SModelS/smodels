! ===========================================================================================================
module xx_integral_le
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  public :: IFCT_LE_X12
  private :: LUMINOSITY_LE
contains
! ------------------------------
  function IFCT_LE_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 
    integer                            :: n
    real(kind=double), dimension(1:30) :: massin,massin_s3
    real(kind=double), dimension(3)    :: lumi
    real(kind=double)  :: m1,s,qf,qr,beta,del_s4,gams,beta1,beta2,sin_stop,mq
    real(kind=double)  :: sw,alpha_s,nlo,ALPHAS,zero1,zero2,u1,u2,t1
    real(kind=double)  :: s3s3m,s3s3p,s3s3_jac,s3s3_x
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac,t2m,t2p,t2,t2_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: z3m,z3p,z3,s3s3,theta_s3,prop_s3,beta1s3,beta2s3,t2s3m,t2s3p,t2s3,t2s3_jac
    real(kind=double)  :: s3m,s3p,s3,s3_jac,betax,t2xm,t2xp,t2x,t2x_jac,s4s3m,s4s3p,s4s3,s4s3_jac
    real(kind=double)  :: s4xm,s4xp,s4x,s4x_jac
    real(kind=double)  :: LE_QGB,LE_QGV,LE_QG,LE_QGD,LE_GG,LE_GGOS,LE_QQ,LE_QB,LE_QBOS

    if (ii>10) then                                            ! finish early 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    massin(1:30) = 0.0                                         ! initialize the massin arrays

    zero1  = 1.0                                               ! needed later to remove ps points in os subtraction 
    zero2  = 1.0

    m1 = mass_s(1)

    x1m    = m1**2 /sc                                         ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = m1**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                           ! partonic cm energy

    if (iscaling==1) s = m1**2 * (eta+1)                       ! overwrite integration for scaling fct

    theta_s3 = 0.0                                             ! the os subtraction theta function, only local
    if (s>4.0*m1**2) theta_s3 = 1.0                            ! intermediate stop

    gams   = ewi * m1                                          ! always a small width needed for arcustan 

    del_s4 = eps_sli * m1**2 * (1.0-m1**2/s)                   ! rescale the s4 cutoff 
    s4p    = 0.0                                               ! only for real corrections needed 
    s4     = 0.0                                               ! since it is used later on in massin
    s3     = 0.0                                               ! dto

    select case (ii)                                           ! all the phase spaces 
    case(-1,0,1,2)                                             ! born_lo, born_nlo, virt

       beta   = 1.0 - m1**2/s                                  ! t2 integration
       t2m    = -1.0/2.0 * ( s - m1**2 + s*beta )                
       t2p    = -1.0/2.0 * ( s - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

    case(3,4,7,8)                                              ! real(qg,gg,qq,qb), no os subtraction 
       
       beta   = (1.0-del_s4/s)**2 - m1**2/s                    ! t2 integration shifted 
       t2m    = -1.0/2.0 * ( s - del_s4 - m1**2 + s*beta )
       t2p    = -1.0/2.0 * ( s - del_s4 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       s4m    = del_s4                                         ! s4 mapped to log
       s4p    = s + t2 - m1**2                                 !  -> only a slight improvement 
       s4     = s4m * (s4p/s4m)**var(4)
       s4_jac = s4 * log(s4p/s4m)

    case(5,9)                                                  ! s3s on-shell subtraction
       
       if (theta_s3==1.0) then 
          dsig = 0.0
          return
       end if

       s3m    = 0.0                                            ! s3 mapped to atan
       s3p    = s + m1**2 - 2.0 * sqrt(s*m1**2)

       s3     = var(3) * (s3p-s3m) + s3m 
       s3_jac = s3p-s3m

       beta1   = (s-s3+m1**2)**2 - 4.0*m1**2*s 
       if ( beta1<0.0 ) then 
          n_faulty = n_faulty+1
          dsig = 0.0 
          return
       end if
          
       beta1   = sqrt( beta1 )
       s4m    = 1.0/2.0 * (s - s3 - m1**2 - beta1)
       s4p    = 1.0/2.0 * (s - s3 - m1**2 + beta1)
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m

       beta2   = abs( s-s4-m1**2 )
       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + beta2 )          ! angular integration re-written
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 - beta2 )
       t2     = var(5) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 

       t2_jac = t2_jac * 2.0*(s4+m1**2)/s4/beta2               ! the additional jacobian

    case(6,10)                                                 ! s3s on-shell subtraction
       
       if (theta_s3==0.0) then 
          dsig = 0.0
          return
       end if

       s3m    = 0.0                                            ! s3 mapped to atan
       s3p    = s + m1**2 - 2.0 * sqrt(s*m1**2)

       z3m    = atan( (s3m - m1**2)/m1/gams )
       z3p    = atan( (s3p - m1**2)/m1/gams ) 
       z3     = var(3) * (z3p-z3m) + z3m 

       s3     = m1*gams*tan(z3) + m1**2
       s3_jac = ((s3-m1**2)**2/m1/gams+m1*gams)*(z3p-z3m)

       beta1   = (s-s3+m1**2)**2 - 4.0*m1**2*s 
       if ( beta1<0.0 ) then 
          n_faulty = n_faulty+1
          dsig = 0.0 
          return
       end if
          
       beta1   = sqrt( beta1 )
       s4m    = 1.0/2.0 * (s - s3 - m1**2 - beta1)
       s4p    = 1.0/2.0 * (s - s3 - m1**2 + beta1)
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m

       beta2   = abs( s-s4-m1**2 )
       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + beta2 )          ! angular integration re-written
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 - beta2 )
       t2     = var(5) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 

       t2_jac = t2_jac * 2.0*(s4+m1**2)/s4/beta2               ! the additional jacobian

       s3s3     = m1**2                                        ! restricted phase space for theta_s3=1 

       beta1s3   = sqrt( (s-s3s3+m1**2)**2 - 4.0*m1**2*s )
       s4s3m    = 1.0/2.0 * (s - s3s3 - m1**2 - beta1s3)
       s4s3p    = 1.0/2.0 * (s - s3s3 - m1**2 + beta1s3)
       s4s3     = var(4) * (s4s3p-s4s3m) + s4s3m
       s4s3_jac = s4s3p-s4s3m
          
       beta2s3  = abs( s-s4s3-m1**2 ) 
       t2s3m    = -1.0/2.0 * ( s - s4s3 - m1**2 + beta2s3 )    ! angular integration re-written
       t2s3p    = -1.0/2.0 * ( s - s4s3 - m1**2 - beta2s3 )
       t2s3     = var(5) * (t2s3p-t2s3m) + t2s3m 
       t2s3_jac = t2s3p-t2s3m 

       t2s3_jac = t2s3_jac * 2.0*(s4s3+m1**2)/s4s3/beta2s3
       
       prop_s3  = m1**2*gams**2 / ((s3-m1**2)**2+m1**2*gams**2)! compensate for the on-shell breit-wigner

    end select

    if (isca==0) then                                          ! renormalization/factorization scale
       qr = scafac * m1
       qf = scafac * m1
    else
       print*, " IFCT1: isca not set correctly ",isca
       stop
    end if

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       nlo = 4.0 * pi * ALPHAS(qr,2)
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! coupling factor alpha_s not always nlo 
       alpha_s = ALPHAS(qr,1)
    else if (ii>0) then 
       alpha_s = ALPHAS(qr,2)
    end if

    massin(1)  = s                                             ! assign the mass arrays
    massin(2)  = t2
    massin(3)  = s4                                            ! s4 only for real phase space
    massin(4)  = s3
    massin(6)  = m1
    massin(7)  = 0.0
    massin(8)  = mg
    massin(9)  = ms
    massin(10) = msq(+6) 
    massin(11) = mt                                            ! top mass in loops (problem for sbottoms with SUSY and SM?)
    massin(12) = qr                                            ! renormalization scale 
    massin(13) = qf                                            ! factorization scale 
    massin(14) = mt                                            ! top mass in decoupled alpha_s
    massin(20) = del_s4                                        ! slicing parameter only for real(qg)
    massin(21) = s4p                                           ! s4^max for the log(delta) terms 
    massin(25) = gams                                          ! width of the stop only for ii=???

    if ( (ii==6).or.(ii==10) ) then                            ! for the s3s on-shell subtraction
       massin_s3(1:30) = massin(1:30)
       massin_s3(2)    = t2s3
       massin_s3(3)    = s4s3
       massin_s3(4)    = s3s3
    end if

    call LUMINOSITY_LE(ii,x1,x2,qf,lumi)

    sin_stop = 1.0                                                                        ! no stop2

    select case (ii)
    case(-1,0,1)                                                                          ! born
       dsig =           LE_QGB(massin,lumi)           *t2_jac                   /s**2 
    case(2)                                                                               ! virtual
       dsig = nlo*      LE_QGV(sin_stop,massin,lumi)  *t2_jac                   /s**2
    case(3)                                                                               ! real qg
       dsig = nlo*   (  LE_QG(massin,lumi)                                                  &       
                      + LE_QGD(massin,lumi)         ) *t2_jac  *s4_jac          /s**2    
    case(4)                                                                               ! real gg
       dsig = nlo*      LE_GG(massin,lumi)            *t2_jac  *s4_jac          /s**2    
    case(5)                                                                               ! on-shell(gg) 
       dsig = nlo*      LE_GGOS(massin,lumi)          *s4_jac  *t2_jac  *s3_jac /s**2
    case(6)                                                                               ! on-shell(gg) 
       dsig = nlo*zero1*LE_GGOS(massin,lumi)          *s4_jac  *t2_jac  *s3_jac /s**2          &
            - nlo*zero2*LE_GGOS(massin_s3,lumi)       *s4s3_jac*t2s3_jac*s3_jac /s**2 * prop_s3
    case(7)                                                                               ! real qq
       dsig = nlo*      LE_QQ(massin,lumi)            *t2_jac  *s4_jac          /s**2    
    case(8)                                                                               ! real qb
       dsig = nlo*      LE_QB(massin,lumi)            *t2_jac  *s4_jac          /s**2    
    case(9)                                                                               ! on-shell(qb) 
       dsig = nlo*      LE_QBOS(massin,lumi)          *s4_jac  *t2_jac  *s3_jac /s**2
    case(10)                                                                              ! on-shell(qb) 
       dsig = nlo*zero1*LE_QBOS(massin,lumi)          *s4_jac  *t2_jac  *s3_jac /s**2          &
            - nlo*zero2*LE_QBOS(massin_s3,lumi)       *s4s3_jac*t2s3_jac*s3_jac /s**2 * prop_s3
    case default
       dsig = 0.0 
    end select

    if (iscaling==0)                                                                   &  ! some jacobians
         dsig = dsig * x1_jac * x2_jac
    
    if (iscaling==0)                                                                   &  ! scaling function
         dsig = dsig * 4.0/m1**2 * 4.0*pi*alpha_s

    if (iscaling==0)                                                                   &  ! result in pb 
         dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_LE_X12

! ------------------------------
  subroutine LUMINOSITY_LE(ii,x1,x2,qf,lumi)

    integer,                         intent(in)  :: ii 
    real(kind=double),               intent(in)  :: x1,x2,qf
    real(kind=double), dimension(3), intent(out) :: lumi

    integer                                      :: i1
    real(kind=double), dimension(-6:6)           :: pdf1,pdf2,pdf3                       ! strictly quark densities for '+'

    lumi(1:3) = 0.0

    if (ii<=0) then                                                                 ! call structure functions for LQ+lepton
       call GET_PDF(0,x1,qf,pdf1)
       call GET_PDF(0,x2,qf,pdf2)
    else 
       call GET_PDF(1,x1,qf,pdf1)
       call GET_PDF(1,x2,qf,pdf2)
    end if

    if (icoll==0) then                                                              ! swap q and qbar for Tevatron
       do i1=-6,6,1
          pdf3(i1) = pdf2(-i1)
       end do
       pdf2(-6:6) = pdf3(-6:6)
    end if
    
    select case (ii)                                                       
    case(-1,0,1,2,3)                                                                ! born type (qg incoming)
          lumi(1)   =  pdf1(ifla_le) * pdf2(0)                   &
                     + pdf2(ifla_le) * pdf1(0)
          lumi(2:3) =  0.0
    case(4,5,6)                                                                     ! crossed (gg incoming)
       lumi(1)   =  pdf1(0) * pdf2(0)                                               ! note the factor 1/2 built in
       lumi(2:3) =  0.0
    case(7)                                                                         ! crossed (qq incoming) 
       lumi(1)   =  pdf1(ifla_le) * sum( pdf2(1:5) )          &
                  + pdf2(ifla_le) * sum( pdf1(1:5) )
       lumi(2)   =  pdf1(ifla_le) * sum( pdf2(1:5) )          &
                  + pdf2(ifla_le) * sum( pdf1(1:5) )
       lumi(3)   =  pdf1(ifla_le) * pdf2(ifla_le)    
    case(8,9,10)                                                                    ! crossed (qb incoming)  
          lumi(1)   =  pdf1(ifla_le) * sum( pdf2(-5:-1) )        &
                     + pdf2(ifla_le) * sum( pdf1(-5:-1) ) 
          lumi(2)   =  pdf1(ifla_le) * pdf2(-ifla_le)               &                     ! no factor 1/2 for q-qbar
                     + pdf2(ifla_le) * pdf1(-ifla_le)
          lumi(3)   =   pdf1( 1) * pdf2(-1) + pdf1(-1) * pdf2( 1) &
                     +  pdf1( 2) * pdf2(-2) + pdf1(-2) * pdf2( 2) &
                     +  pdf1( 3) * pdf2(-3) + pdf1(-3) * pdf2( 3) &
                     +  pdf1( 4) * pdf2(-4) + pdf1(-4) * pdf2( 4) &
                     +  pdf1( 5) * pdf2(-5) + pdf1(-5) * pdf2( 5) 
    end select

  end subroutine LUMINOSITY_LE

end module xx_integral_le






