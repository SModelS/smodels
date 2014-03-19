! ===========================================================================================================
module xx_integral_ng
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: COUPLING_NG
  public :: IFCT_NG_X12
contains
! ------------------------------
  function IFCT_NG_X12(dum) result(dsig)
    
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 
    integer                            :: inlo,iq,is,il1,il2
    real(kind=double), dimension(1:30) :: massin,massix,massix_s3,massin_s4
    real(kind=double), dimension(-6:6) :: pdf1,pdf2
    real(kind=double), dimension(-1:1) :: mst,msu
    real(kind=double)  :: m1,m1_sign,s,qf,qr,beta,del_s4,gams,beta1,beta2
    real(kind=double)  :: sw,alpha,alpha_s,nlo,ALPHAS,LUMI,lumi_ex
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac,tgm,tgp,tg,tg_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: z4m,z4p,z4,s4s4,theta_s4,prop_s4,betas4,tgs4m,tgs4p,tgs4,tgs4_jac
    real(kind=double)  :: z3m,z3p,z3,s3s3,theta_s3,prop_s3,beta1s3,beta2s3,tgs3m,tgs3p,tgs3,tgs3_jac
    real(kind=double)  :: s3m,s3p,s3,s3_jac,tgxm,tgxp,tgx,tgx_jac,s4s3m,s4s3p,s4s3,s4s3_jac
    real(kind=double)  :: s4xm,s4xp,s4x,s4x_jac
    real(kind=double)  :: NG_QBB,NG_QBB_NG,NG_QBV,NG_QBH,NG_QBD,NG_QGH,NG_GBH,NG_QGOS1,NG_QGOS2,NG_GBOS1,NG_GBOS2
    complex(kind=double), dimension(1:4) :: Cl,Cr,Cv

    if (ii>5) then                                             ! finish early for inclusive case 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    massin(1:30) = 0.0                                         ! initialize the massin arrays
    massix(1:30) = 0.0

    m1      = abs( mass_n(ipart1) )                            ! assign the final state masses 
    m1_sign = mass_n(ipart1)                                   ! this one has the sign included

    if ( (abs(m1)+abs(mg))**2 > 0.98*sc ) then
!tp       print*, " collider energy not large enough ",m1,mg,sqrt(sc)
       dsig = 0.0
       return
    end if

    sw    = sqrt( 1.0 - mw**2/mz**2 )                          ! weak parameters in the on-shell scheme 
    alpha = sqrt(2.0) * mw**2 * sw**2 /pi * gf  

    x1m    = (m1+mg)**2 /sc                                    ! x1-x2 integration, map x->log(x)
    x1p    = 1.0
    x1     = x1m * (x1p/x1m)**var(1)
    x1_jac = x1 * log(x1p/x1m)

    x2m    = (m1+mg)**2 /sc /x1
    x2p    = 1.0
    x2     = x2m * (x2p/x2m)**var(2)
    x2_jac = x2 * log(x2p/x2m)

    s = x1 * x2 * sc                                           ! partonic cm energy
    if (iscaling==1) s = (m1+mg)**2 * (eta+1)                  ! overwrite integration for scaling fct

    if (isca==0) then                                          ! renormalization/factorization scale
       qf = scafac * (m1+mg)/2.0
       qr = scafac * (m1+mg)/2.0
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

    if ((s>(ms+mg)**2).and.(ms>m1)) theta_s4 = 1.0 
    if ((s>(ms+m1)**2).and.(ms>mg)) theta_s3 = 1.0 

    del_s4 = eps_sli * min(m1**2,mg**2)                       &! unit is m^2
                   * (1.0-(m1+mg)**2/s)*(1.0-(m1-mg)**2/s)     ! rescale the s4 cutoff to be small everywhere 
    s4p    = 0.0                                               ! only for real corrections needed 

    select case (ii)                                           ! all the phase spaces 
    case(-1,0,1,2)                                             ! born_lo, born_nlo, virt

       beta   = sqrt(1.0-(m1+mg)**2/s) * sqrt(1.0-(m1-mg)**2/s)               ! tg integration
       tgm    = -1.0/2.0 * ( s + mg**2 - m1**2 + s*beta )                
       tgp    = -1.0/2.0 * ( s + mg**2 - m1**2 - s*beta )
       tg     = var(3) * (tgp-tgm) + tgm
       tg_jac = tgp-tgm 

    case(3)                                                                    ! real(qb), no os subtraction 
       
       beta   = sqrt((1.0-del_s4/s)**2-(m1+mg)**2/s)                   &       ! tg integration shifted 
               *sqrt((1.0-del_s4/s)**2-(m1-mg)**2/s)                    
       tgm    = -1.0/2.0 * ( s - del_s4 + mg**2 - m1**2 + s*beta )
       tgp    = -1.0/2.0 * ( s - del_s4 + mg**2 - m1**2 - s*beta )
       tg     = var(3) * (tgp-tgm) + tgm
       tg_jac = tgp-tgm 

       s4m    = del_s4                                                         ! s4 mapped to log
       s4p    = s + tg + mg**2 - m1**2 + s*mg**2/tg                            !  -> only a slight improvement 
       s4     = s4m * (s4p/s4m)**var(4)
       s4_jac = s4 * log(s4p/s4m)

    case(4)                                                                    ! real(qg+gb) w/o os subtraction

       s4m    = 0.0                                                            ! s4 mapped to atan
       s4p    = s + mg**2 - m1**2 - 2.0 * sqrt(s*mg**2)
       z4m    = atan( (s4m + m1**2 - ms**2)/ms/gams )
       z4p    = atan( (s4p + m1**2 - ms**2)/ms/gams ) 
       z4     = var(3) * (z4p-z4m) + z4m 

       s4     = ms*gams*tan(z4) + ms**2 - m1**2   
       s4_jac = ((s4+m1**2-ms**2)**2/ms/gams+ms*gams)*(z4p-z4m) 

       beta   = sqrt( (s-s4-m1**2+mg**2)**2 - 4.0*s*mg**2 )
       tgm    = -1.0/2.0 * ( s - s4 - m1**2 + mg**2 + beta )                   ! tg integration as usual
       tgp    = -1.0/2.0 * ( s - s4 - m1**2 + mg**2 - beta )
       tg     = var(4) * (tgp-tgm) + tgm 
       tg_jac = tgp-tgm 

    case(5)                                                                    ! real(qg+gb) os subtraction

       s4m    = 0.0                                                            ! s4 mapped to atan
       s4p    = s + mg**2 - m1**2 - 2.0 * sqrt(s*mg**2)
       z4m    = atan( (s4m + m1**2 - ms**2)/ms/gams )
       z4p    = atan( (s4p + m1**2 - ms**2)/ms/gams ) 
       z4     = var(3) * (z4p-z4m) + z4m 

       s4     = ms*gams*tan(z4) + ms**2 - m1**2   
       s4_jac = ((s4+m1**2-ms**2)**2/ms/gams+ms*gams)*(z4p-z4m) 
          
       beta   = sqrt( (s-s4-m1**2+mg**2)**2 - 4.0*s*mg**2 )
       tgm    = -1.0/2.0 * ( s - s4 - m1**2 + mg**2 + beta )                   ! tg integration as usual
       tgp    = -1.0/2.0 * ( s - s4 - m1**2 + mg**2 - beta )
       tg     = var(4) * (tgp-tgm) + tgm 
       tg_jac = tgp-tgm 

       if (theta_s4==1.0) then                                                 ! restricted phase space for theta_s4=1 
          s4s4 = ms**2 - m1**2

          betas4   = sqrt( (s-s4s4-m1**2+mg**2)**2 - 4.0*s*mg**2 )
          tgs4m    = -1.0/2.0 * ( s - s4s4 - m1**2 + mg**2 + betas4 )
          tgs4p    = -1.0/2.0 * ( s - s4s4 - m1**2 + mg**2 - betas4 )
          tgs4     = var(4) * (tgs4p-tgs4m) + tgs4m 
          tgs4_jac = tgs4p-tgs4m 

          prop_s4 = ms**2*gams**2 / ((s4+m1**2-ms**2)**2+ms**2*gams**2)        ! compensate for the wrong breit-wigner
       end if

       s3m    = 0.0                                                            ! s3 mapped to atan
       s3p    = s + m1**2 - mg**2 - 2.0 * sqrt(s*m1**2) 
       z3m    = atan( (s3m + mg**2 - ms**2)/ms/gams )
       z3p    = atan( (s3p + mg**2 - ms**2)/ms/gams ) 
       z3     = var(3) * (z3p-z3m) + z3m 
       
       s3     = ms*gams*tan(z3) + ms**2 - mg**2   
       s3_jac = ((s3+mg**2-ms**2)**2/ms/gams+ms*gams)*(z3p-z3m)

       beta1   = sqrt( (s-s3-mg**2+m1**2)**2 - 4.0*m1**2*s )
       s4xm    = s3/2.0/(s3+mg**2) * (s - s3 - m1**2 - mg**2 - beta1)
       s4xp    = s3/2.0/(s3+mg**2) * (s - s3 - m1**2 - mg**2 + beta1)
       s4x     = var(4) * (s4xp-s4xm) + s4xm
       s4x_jac = s4xp-s4xm
       
       beta2  = sqrt( (s-s4x-m1**2+mg**2)**2 - 4.0*mg**2*s ) 
       tgxm    = -1.0/2.0 * ( s - s4x - m1**2 + mg**2 + beta2 )                ! angular integration re-written
       tgxp    = -1.0/2.0 * ( s - s4x - m1**2 + mg**2 - beta2 )
       tgx     = var(5) * (tgxp-tgxm) + tgxm 
       tgx_jac = tgxp-tgxm 

       tgx_jac = tgx_jac * 2.0*(s4x+m1**2)/s4x/beta2                           ! the additional jacobian

       if (theta_s3==1.0) then                                                 ! the complete s3 type phase space
          s3s3 = ms**2 - mg**2                                                 ! restricted phase space for theta_s3=1 
         
          beta1s3   = sqrt( (s-s3s3-mg**2+m1**2)**2 - 4.0*m1**2*s )
          s4s3m    = s3s3/2.0/(s3s3+mg**2) * (s - s3s3 - m1**2 - mg**2 - beta1s3)
          s4s3p    = s3s3/2.0/(s3s3+mg**2) * (s - s3s3 - m1**2 - mg**2 + beta1s3)
          s4s3     = var(4) * (s4s3p-s4s3m) + s4s3m
          s4s3_jac = s4s3p-s4s3m

          beta2s3  = sqrt( (s-s4s3-m1**2+mg**2)**2 - 4.0*s*mg**2 ) 
          tgs3m    = -1.0/2.0 * ( s - s4s3 - m1**2 + mg**2 + beta2s3 )         ! angular integration re-written
          tgs3p    = -1.0/2.0 * ( s - s4s3 - m1**2 + mg**2 - beta2s3 )
          tgs3     = var(5) * (tgs3p-tgs3m) + tgs3m 
          tgs3_jac = tgs3p-tgs3m 

          tgs3_jac = tgs3_jac * 2.0*(s4s3+m1**2)/s4s3/beta2s3

          prop_s3 = ms**2*gams**2 / ((s3+mg**2-ms**2)**2+ms**2*gams**2)        ! compensate for the wrong breit-wigner 
       end if

    end select

    massin(1)  = s                                             ! assign the mass arrays
    massin(2)  = tg
    massin(3)  = s4                                            ! s4 only for ii=3,7,9
    massin(6)  = m1_sign
    massin(7)  = mg
    massin(9)  = mt
    massin(11) = ms
    massin(12) = qr                                            ! renormalization scale 
    massin(13) = qf                                            ! factorization scale 
    massin(20) = del_s4                                        ! slicing parameter only fort ii=3
    massin(21) = s4p                                           ! s4^max for the log(delta) terms 
    massin(25) = gams                                          ! width of the squark only for ii=7,9
    massin(26) = 1.0                                           ! to tag 1/s4^2 term only for ii=8,10 
    massin(27) = 0.0                                           ! to tag 1/s3^2 term only for ii=8,10

    if (ii==5) then                                            ! define restricted phase space 
       massix(1:30)    = massin(1:30)
       massix(2)       = tgx
       massix(3)       = s4x
       massix(4)       = s3
       massix(26)      = 0.0                                   ! to tag 1/s4^2 term 
       massix(27)      = 1.0                                   ! to tag 1/s3^2 term 
       if (theta_s3==1.0) then 
          massix_s3(1:30) = massin(1:30)
          massix_s3(2)    = tgs3
          massix_s3(3)    = s4s3
          massix_s3(4)    = s3s3
          massix_s3(26)   = 0.0
          massix_s3(27)   = 1.0
       end if
       if (theta_s4==1.0) then 
          massin_s4(1:30) = massin(1:30)
          massin_s4(2)    = tgs4
          massin_s4(3)    = s4s4 
          massin_s4(26)   = 1.0
          massin_s4(27)   = 0.0
       end if
    end if

    if (ii==-1) then
       call GET_PDF(inlo,x1,qf,pdf1)                           ! call structure functions once forever
       call GET_PDF(inlo,x2,qf,pdf2)
    end if

    dsig = 0.0 
    do iq = -1,1,2                                             ! the loop for up-type and down-type quarks 

       if ( ((ipart1==5).or.(ipart1==6)).and.(iq==-1) ) then   ! charge of final state has to be positive
          dsig = dsig + 0.0
          cycle
       end if
 
       if ( ((ipart1==7).or.(ipart1==8)).and.(iq==+1) ) then   ! charge of final state has to be negative
          dsig = dsig + 0.0
          cycle
       end if
 
       call COUPLING_NG(iq,Cl,Cr,Cv)

       select case (ii)
       case(-1)

          do is = 1,4,1                                                                         ! t-channel squarks
             
             if (ipart1<5) then                                                                 ! N case

                if ((is==2).or.(is==3)) then
                   if (iq==+1) cycle
                else if ((is==1).or.(is==4)) then
                   if (iq==-1) cycle
                end if
                il1 = is
                il2 = il1

                msu(-1:1:2) = msq(-is:is:2*is) 
                
             else if ( (ipart1==5).or.(ipart1==6) ) then                                        ! C+ case

                if ( (is==1).or.(is==4) ) cycle                                                 ! only sdown-type in t channel

                if (is==2) then                                                                 ! u channel attached to il1
                   il1 = 1                                                                      ! quark of il1: u
                   il2 = 2                                                                      ! antiquark of il2; dbar
                   msu(-1:1:2) = msq(-1:1:2)
                else if (is==3) then
                   il1 = 4
                   il2 = 3
                   msu(-1:1:2) = msq(-4:4:8)
                end if
                
             else if ( (ipart1==7).or.(ipart1==8) ) then                                        ! C- case
            
                if ( (is==2).or.(is==3) ) cycle                                                 ! only sup-type in t channel
                                
                if (is==1) then                                                                 ! u channel attached to il1
                   il1 = 2                                                                      ! quark of il1: d
                   il2 = 1                                                                      ! antiquark of il2: ubar
                   msu(-1:1:2) = msq(-2:2:4)
                else if (is==4) then
                   il1 = 3
                   il2 = 4
                   msu(-1:1:2) = msq(-3:3:6)
                end if
                
             end if
                
             if (i_ngtest==0) then 
                mst(-1:1:2) = msq(-is:is:2*is)                                                  ! this is the definiton of is
             else if (i_ngtest==1) then 
                mst(-1:1:2) = ms
                msu(-1:1:2) = ms
             else 
                print*, " IFCT_NG_X12: i_ngtest not set "
                call HARD_STOP
             end if

             if (icoll==0) then                                                                 ! Tevatron   
                lumi_ex =   pdf1( il1)*pdf2( il2) + pdf1(-il2)*pdf2(-il1)                       ! il1=q; il2=qbar
             else                                                                               ! LHC
                lumi_ex =   pdf1( il1)*pdf2(-il2) + pdf1(-il2)*pdf2( il1)                       ! il1=q; il2=qbar 
             end if

             dsig = dsig +  lumi_ex                                                           &
                            * NG_QBB_NG(massin,Cl,Cr,mst,msu)/s**2* tg_jac
          end do
       case(0,1)                                                                                ! born
          dsig = dsig +       LUMI(inlo,20,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * NG_QBB(massin,Cl,Cr)     *tg_jac 
       case(2)                                                                                  ! virt+soft (incl. scale) 
          dsig = dsig + nlo * LUMI(inlo,20,icoll,idub,iq,x1,x2,qf)/s**2                        & 
                            * NG_QBV(massin,Cl,Cr,Cv)  *tg_jac 
       case(3)                                                                                  ! real(qb) 
          dsig = dsig + nlo * LUMI(inlo,20,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_QBH(massin,Cl,Cr)     *tg_jac  *s4_jac                     &
                               + NG_QBD(massin,Cl,Cr)     *tg_jac  *s4_jac )
       case(4)                                                                                  ! real(qg)
          dsig = dsig + nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_QGH(massin,Cl,Cr)                                          &
                               + NG_QGOS2(massin,Cl,Cr) ) *tg_jac  *s4_jac

          dsig = dsig + nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_GBH(massin,Cl,Cr)                                          &
                               + NG_GBOS2(massin,Cl,Cr) ) *tg_jac  *s4_jac  
       case(5)                                                                                  ! on-shell(qg) 
          dsig = dsig + nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_QGOS1(massin,Cl,Cr)   *tg_jac  *s4_jac                     &  
                               + NG_QGOS1(massix,Cl,Cr)   *s4x_jac *tgx_jac *s3_jac )
          if (theta_s4==1.0)                                                                   &
          dsig = dsig - nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_QGOS1(massin_s4,Cl,Cr)*tgs4_jac*s4_jac          * prop_s4 )
          if (theta_s3==1.0)                                                                   & 
          dsig = dsig - nlo * LUMI(inlo,30,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_QGOS1(massix_s3,Cl,Cr)*s4s3_jac*tgs3_jac*s3_jac * prop_s3 )   

          dsig = dsig + nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_GBOS1(massin,Cl,Cr)   *tg_jac  *s4_jac                     &  
                               + NG_GBOS1(massix,Cl,Cr)   *s4x_jac *tgx_jac *s3_jac )
          if (theta_s4==1.0)                                                                   &
          dsig = dsig - nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_GBOS1(massin_s4,Cl,Cr)*tgs4_jac*s4_jac          * prop_s4 )   
          if (theta_s3==1.0)                                                                   &
          dsig = dsig - nlo * LUMI(inlo,40,icoll,idub,iq,x1,x2,qf)/s**2                        &
                            * (  NG_GBOS1(massix_s3,Cl,Cr)*s4s3_jac*tgs3_jac*s3_jac * prop_s3 )   
       case default
          dsig = 0.0 
       end select
    end do

    if (iscaling==0) dsig = dsig * x1_jac * x2_jac

    if (iscaling==0) dsig = dsig * 4.0/(m1+mg)**2 * 4.D0*pi*alpha * 4.D0*pi*alpha_s
    
    if (iscaling==0) dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_NG_X12


! ------------------------------
! all general couplings : C(1) located at ( k1(q_in) , N/C )
!                         C(2) located at ( k2(q_out), gl  )
!                         C(3) located at ( k2(q_in) , gl  )
!                         C(4) located at ( k2(q_out), N/C )  for all Cl,Cr
! ------------------------------
  subroutine COUPLING_NG(iq,Cl,Cr,Cv)
    integer,                              intent(in)   :: iq
    complex(kind=double), dimension(4),   intent(out)  :: Cl,Cr,Cv

    integer              :: ic 
    real(kind=double)    :: t3,qq,sw2,cw2,sw,cw
    complex(kind=double) :: zzc(4,4)

    t3(iq)  = dble(iq) /2.0                                                  ! quark quantum numbers 
    qq(iq)  = 2.0/3.0 + ( dble(iq) - 1.0 ) /2.0

    sw = sqrt( 1.0 - mw**2/mz**2 )
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)

    zzc(1:4,1:4) =  conjg( zz(1:4,1:4) )                                     ! complex conjugation of bw=zz

    Cl(1:4) = 0.0                                                            ! initialize arrays
    Cr(1:4) = 0.0
    Cv(1:4) = 0.0

    ic = ipart1 - 4                                                          ! set chargino number 1 or 2  
    if (ic>2) ic = ic - 2

    select case (ipart1)                                                     ! only non-zero couplings 
    case(1,2,3,4)                                                            ! neutralino and gluino 

       Cl(1) = ( zz(ipart1,1)*sw*(qq(iq)-t3(iq))+zz(ipart1,2)*cw*t3(iq))/(sqrt(2.0)*sw*cw) ! nqs coupings [gauginos]
       Cr(1) =  zzc(ipart1,1)*sw* qq(iq)                                /(sqrt(2.0)*sw*cw) ! in wim's normalization 
       Cl(3) = Cl(1)
       Cr(3) = Cr(1) 

    case(5,6)                                                                ! chargino+ and gluino

       Cl(1) = uu(ic,1)/(2.0*sw)                                             ! t channel with s-down for C+
       Cl(3) = vv(ic,1)/(2.0*sw)                                             ! u channel with s-up for C+

    case(7,8)                                                                ! chargino- and gluino

       Cl(3) = uu(ic,1)/(2.0*sw)                                             ! t channel with s-up for C-
       Cl(1) = vv(ic,1)/(2.0*sw)                                             ! u channel with s-down for C-

    end select
    
    Cv(1) = Cl(3)/conjg(Cl(1))                                               ! set {Clot,Cupt,Cupu,Clou}
    Cv(3) = conjg(Cl(1))/Cl(3)                                               !  [identical for L and R]

  end subroutine COUPLING_NG

end module xx_integral_ng



