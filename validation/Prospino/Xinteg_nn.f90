! ===========================================================================================================
module xx_integral_nn
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  implicit none 
  private :: KIN_NN,COUPLING_NN
  public :: IFCT_NN_X12
contains
! ------------------------------
  function IFCT_NN_X12(dum) result(dsig)
    real(kind=double), dimension(dim(ii)), intent(in) :: dum ! vegas integration variable
    real(kind=double), dimension(dim(ii))             :: var ! internal integration variable 
    real(kind=double)                  :: dsig 
    integer                            :: inlo,iq,iout,is,il1,il2
    real(kind=double), dimension(1:30) :: massin,massin_s3,massin_s4,massin_x3,massin_x4
    real(kind=double), dimension(1:99) :: mkraemer
    real(kind=double), dimension(1:4)  :: Cs,Csx
    real(kind=double), dimension(-6:6) :: pdf1,pdf2
    real(kind=double), dimension(-1:1) :: mst,msu
    real(kind=double)  :: m1,m2,delta_soft,symmfac,s,mu,beta,beta1,beta2,mx,gams
    real(kind=double)  :: sw,alpha,nlo,ALPHAS,LUMI,lumi_ex
    real(kind=double)  :: x1m,x1p,x1,x1_jac,x2m,x2p,x2,x2_jac
    real(kind=double)  :: t2m,t2p,t2,t2_jac,s4m,s4p,s4,s4_jac
    real(kind=double)  :: c1m,c1p,c1,c1_jac,axm,axp,ax,ax_jac
    real(kind=double)  :: u1,t1,u2,tp,up,s5
    real(kind=double)  :: z4m,z4p,z4,z3m,z3p,z3
    real(kind=double)  :: s4s4,betas4,t2s4m,t2s4p,t2s4,t2s4_jac,prop_s4,theta_s4,u2s3,t1s3,u1s3
    real(kind=double)  :: s4s3,beta1s3,beta2s3,t2s3m,t2s3p,t2s3,t2s3_jac,prop_s3,theta_s3,u2s4,t1s4,u1s4
    real(kind=double)  :: s3s4,s5s4,tps4,ups4
    real(kind=double)  :: s3s3,s5s3,tps3,ups3,s4s3m,s4s3p,s4s3_jac
    real(kind=double)  :: s3m,s3p,s3,s3_jac,betax,t2xm,t2xp,t2x,t2x_jac
    real(kind=double)  :: xm,xp,x,x_jac,sx,u1x,t1x,u2x
    real(kind=double)  :: NN_QBB,NN_QBV,NN_QBH,NN_QBF1,NN_QBF2,NN_QGH,NN_QGHSUB,NN_QGF,NN_GBH,NN_GBHSUB
    real(kind=double)  :: NN_QBB_NG
    complex(kind=double), dimension(1:4) :: Cl,Cr,Ct,Cv,Ctx
!tp needed for the scales check
!tp    real(kind=double)  :: shat,that,uhat

    if (ii>11) then                                            ! finish early for inclusive case 
       dsig = 0.0 
       return
    end if

    var(1:dim(ii)) = dum(1:dim(ii)) * ( 1.0 - 2.0*cut ) + cut  ! cut off the integration in general

    massin(1:30) = 0.0                                         ! initialize the massin array

    m1 = mass_n(ipart1)                                        ! assign the final state masses 
    m2 = mass_n(ipart2)
    if ( (abs(m1)+abs(m2)) < mz ) then
       print*, " IFCT_NN_X12: masses low, Z decays might be more suitable "
       call HARD_STOP
    end if

    if ( (abs(m1)+abs(m2))**2 > 0.98*sc ) then
!tp       print*, " collider energy not large enough ",m1,m2,sqrt(sc)
       dsig = 0.0
       return
    end if

    delta_soft = eps_sub * (abs(m1)+abs(m2))**2/4.0            ! the soft cut-off rescaled 
    
    if (ipart1==ipart2) then                                   ! phase scace symmetry factor 
       symmfac = 1.0/2.0
    else 
       symmfac = 1.0 
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
       mu = scafac*(abs(m1)+abs(m2))/2.0
    else if (isca==1) then
       mu = scafac*sqrt(s)
    end if

    if (iscaling==0) then                                      ! nlo factor [always nlo alpha_s]
       nlo = 4.0 * pi * ALPHAS(mu,2)
    else if (iscaling==1) then 
       nlo = 1.0
    end if

    if (ii<=0) then                                            ! note that there is no consistent alpha_s in leading oder
       inlo = 0
    else if (ii>0) then 
       inlo = 1 
    end if

    theta_s4 = 0.0                                             ! the os subtraction theta function 
    theta_s3 = 0.0 
    gams = ewi * ms 
    if ((s>(ms+abs(m2))**2).and.(ms>abs(m1))) theta_s4 = 1.0 
    if ((s>(ms+abs(m1))**2).and.(ms>abs(m2))) theta_s3 = 1.0 

    select case (ii)                                           ! different phase space integrations 
    case(-1,0,1,2)                                             ! t2 integration for born, virtual 

       beta = sqrt(1.0-(abs(m1)+abs(m2))**2/s) * sqrt(1.0-(abs(m1)-abs(m2))**2/s)

       t2m    = -1/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1/2.0 * ( s + m2**2 - m1**2 - s*beta )
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

    case(3)                                                    ! t2-s4-omega integration for real(qb)

       beta = sqrt(1.0-(abs(m1)+abs(m2))**2/s) * sqrt(1.0-(abs(m1)-abs(m2))**2/s)

       t2m    = -1/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1/2.0 * ( s + m2**2 - m1**2 - s*beta )
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

       call KIN_NN(0,m1,m2,s,t2,s4,c1,ax,s3,s5,u2,t1,u1,tp,up) ! the real phase space 

       if ((abs(tp)<delta_soft).or.(abs(up)<delta_soft).or.(abs(s3)<delta_soft).or.(abs(s4)<delta_soft)) then 
          dsig = 0.0                                           ! cut the soft phase space 
          return 
       end if

    case(4,5)                                                               ! t2-x integration for mass factorization

       beta = sqrt(1.0-(abs(m1)+abs(m2))**2/s) * sqrt(1.0-(abs(m1)-abs(m2))**2/s)

       t2m    = -1/2.0 * ( s + m2**2 - m1**2 + s*beta )
       t2p    = -1/2.0 * ( s + m2**2 - m1**2 - s*beta )
       t2     = var(3) * (t2p-t2m) + t2m
       t2_jac = t2p-t2m 

       u1  = - s - t2  
       t1  = t2 + m2**2 - m1**2 
       u2  = u1 + m1**2 - m2**2 

       xm    = (abs(m1)+abs(m2))**2 /s
       xp    = 1.0
       x     = var(4) * (xp-xm) + xm  
       x_jac = xp-xm 

       sx = s * x                                              ! shift in the variables s

       betax = sqrt(1.0-(abs(m1)+abs(m2))**2/sx) * sqrt(1.0-(abs(m1)-abs(m2))**2/sx)

       t2xm    = -1/2.0 * ( sx + m2**2 - m1**2 + sx*betax )    ! redo t2 integration 
       t2xp    = -1/2.0 * ( sx + m2**2 - m1**2 - sx*betax )
       t2x     = var(3) * (t2xp-t2xm) + t2xm
       t2x_jac = t2xp-t2xm 

       u1x = - sx - t2x                                        ! redo born type kinematics 
       t1x = t2x + m2**2 - m1**2 
       u2x = u1x + m1**2 - m2**2
       
    case(6,7,9,10)                                             ! z-t2-omega for real subtracted (qg,gb)

       s4m    = 0.0                                            ! s4 integration mapped to atan
       s4p    = s + m2**2 - m1**2 - 2.0 * sqrt(s*m2**2)
       z4m    = atan( (s4m + m1**2 - ms**2)/ms/gams )
       z4p    = atan( (s4p + m1**2 - ms**2)/ms/gams ) 
       z4     = var(3) * (z4p-z4m) + z4m 

       s4     = gams*ms*tan(z4) + ms**2 - m1**2   
       s4_jac = ((s4+m1**2-ms**2)**2/gams/ms+gams*ms)*(z4p-z4m) 

       beta = sqrt( (s-s4-m1**2+m2**2)**2 - 4.0*s*m2**2 )

       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 + beta )
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 - beta )
       t2     = var(4) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 

       c1m    = -1.0
       c1p    = 1.0
       c1     = var(5) * (c1p-c1m) + c1m 
       c1_jac = c1p-c1m 

       axm    = 0.0 
       axp    = pi 
       ax     = var(6) * (axp-axm) + axm 
       ax_jac = axp-axm 

       call KIN_NN(0,m1,m2,s,t2,s4,c1,ax,s3,s5,u2,t1,u1,tp,up)

       if (theta_s4==1.0) then                                    ! restricted phase space for theta_s4=1

          s4s4 = ms**2 - m1**2

          betas4 = sqrt( (s-s4s4-m1**2+m2**2)**2 - 4.0*s*m2**2 )

          t2s4m    = -1.0/2.0 * ( s - s4s4 - m1**2 + m2**2 + betas4 )
          t2s4p    = -1.0/2.0 * ( s - s4s4 - m1**2 + m2**2 - betas4 )
          t2s4     = var(4) * (t2s4p-t2s4m) + t2s4m 
          t2s4_jac = t2s4p-t2s4m 
          
          call KIN_NN(0,m1,m2,s,t2s4,s4s4,c1,ax,s3s4,s5s4,u2s4,t1s4,u1s4,tps4,ups4)

          prop_s4 = gams**2*ms**2 / ((s4+m1**2-ms**2)**2+gams**2*ms**2) ! compensate for the wrong breit-wigner

       end if
       
    case(8,11)

       if (theta_s3==0.0) then                                     ! what it's all about 
          dsig = 0.0 
          return
       end if

       s3m    = 0.0                                                ! re-written integration in general
       s3p    = s + m1**2 - m2**2 - 2.0 * sqrt(s*m1**2)
       z3m    = atan( (s3m + m2**2 - ms**2)/ms/gams )
       z3p    = atan( (s3p + m2**2 - ms**2)/ms/gams ) 
       z3     = var(3) * (z3p-z3m) + z3m 
       
       s3     = ms*gams*tan(z3) + ms**2 - m2**2   
       s3_jac = ((s3+m2**2-ms**2)**2/ms/gams+ms*gams)*(z3p-z3m)
       
       beta1   = sqrt( (s-s3-m2**2+m1**2)**2 - 4.0*m1**2*s )
       s4m    = s3/2.0/(s3+m2**2) * (s - s3 - m1**2 - m2**2 - beta1)
       s4p    = s3/2.0/(s3+m2**2) * (s - s3 - m1**2 - m2**2 + beta1)
       s4     = var(4) * (s4p-s4m) + s4m
       s4_jac = s4p-s4m
       
       beta2  = sqrt( (s-s4-m1**2+m2**2)**2 - 4.0*s*m2**2 ) 
       t2m    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 + beta2 )
       t2p    = -1.0/2.0 * ( s - s4 - m1**2 + m2**2 - beta2 )
       t2     = var(5) * (t2p-t2m) + t2m 
       t2_jac = t2p-t2m 
       
       t2_jac = t2_jac * 2.0/s4*(s4+m1**2) / beta2                          ! including the over-all jacobian
       
       axm    = 0.0 
       axp    = pi 
       ax     = var(6) * (axp-axm) + axm 
       ax_jac = axp-axm 

       call KIN_NN(1,m1,m2,s,t2,s4,c1,ax,s3,s5,u2,t1,u1,tp,up)
          
       s3s3 = ms**2 - m2**2

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
       
       t2s3_jac = t2s3_jac * 2.0/s4s3*(s4s3+m1**2) / beta2s3                ! including the over-all jacobian
       
       call KIN_NN(1,m1,m2,s,t2s3,s4s3,c1,ax,s3s3,s5s3,u2s3,t1s3,u1s3,tps3,ups3)
          
       prop_s3 = ms**2*gams**2 / ((s3+m2**2-ms**2)**2+ms**2*gams**2)        ! compensate for the wrong breit-wigner 

       if ((abs(tp)<delta_soft).or.(abs(up)<delta_soft).or.(abs(s3)<delta_soft).or.(abs(s4)<delta_soft)) then 
          dsig = 0.0                                                        ! cut the soft phase space  
          return
       end if

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
    if ((ii==4).or.(ii==5)) then                               ! additional entries for mass factorization  
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
    massin(30) = ewi                                           ! additional entry for crossed channel 

    mkraemer(1:99) = 0.0                                       ! the array mkraemer for real correction 
    mkraemer(1)  = m1 
    mkraemer(2)  = m2 
    mkraemer(3)  = mg 
    mkraemer(4)  = ms 

    if (ii>=6) then                                            ! define restricted phase space 
       if (theta_s4==1.0) then 
          massin_x4(1:30) = massin(1:30) 
          massin_x4(26)   = 1.0                                ! only the s4 os subtraction 
          massin_x4(27)   = 0.0
          massin_s4(1:30) = massin(1:30) 
          massin_s4(2)    = t2s4
          massin_s4(3)    = u2s4 
          massin_s4(4)    = t1s4 
          massin_s4(5)    = u1s4 
          massin_s4(26)   = 1.0
          massin_s4(27)   = 0.0
       end if
       if (theta_s3==1.0) then 
          massin_x3(1:30) = massin(1:30) 
          massin_x3(26)   = 0.0                                ! only the s4 os subtraction 
          massin_x3(27)   = 1.0
          massin_s3(1:30) = massin(1:30)
          massin_s3(2)    = t2s3
          massin_s3(3)    = u2s3
          massin_s3(4)    = t1s3
          massin_s3(5)    = u1s3
          massin_s3(26)   = 0.0
          massin_s3(27)   = 1.0
       end if
    end if

    if (ii==-1) then
       call GET_PDF(inlo,x1,mu,pdf1)                                      ! call structure functions once forever
       call GET_PDF(inlo,x2,mu,pdf2)
    end if

    dsig = 0.0 
    do iq =-1,1,2                                                         ! the loop for up-type and down-type quarks 

       if ( (((ipart1<=4).and.((ipart2==5).or.(ipart2==6))).or.    &      ! charge of final state positive
             ((ipart2<=4).and.((ipart1==5).or.(ipart2==6)))    )   &
           .and.(iq==-1) ) then
          dsig = dsig + 0.0
          cycle
       end if
 
       if ( (((ipart1<=4).and.((ipart2==7).or.(ipart2==8))).or.    &      ! charge of final state negative
             ((ipart2<=4).and.((ipart1==7).or.(ipart2==8)))    )   &
           .and.(iq==+1) ) then
          dsig = dsig + 0.0
          cycle
       end if
 
       call COUPLING_NN(s ,m1,m2,iq,Cs ,Ct ,Cl,Cr,Cv,mx,iout,mkraemer)
       if ((ii==4).or.(ii==5)) call COUPLING_NN(sx,m1,m2,iq,Csx,Ctx,Cl,Cr,Cv,mx,iout,mkraemer)

       if ( (ii>=0) .and. (2.0*ms/(abs(m1)+abs(m2))>20.0) ) then          ! no t channel couplings
          Cl(1:4) = 0.0                           
          Cr(1:4) = 0.0 
          Cv(1:4) = 0.0
       end if

       if ( (ii>=0) .and. ((2.0*ms/(abs(m1)+abs(m2))>100.0).and.(ms>1.e4)) ) then ! no subtraction of intermediate squarks
          massin(26) = 0.0                        
          massin(27) = 0.0
       end if

       massin(8) = mx                                                     ! the s channel particle mass 
       if (ii>=6) then 
          if (theta_s3==1.0) then 
             massin_s3(8) = mx
             massin_x3(8) = mx
          end if
          if (theta_s4==1.0) then 
             massin_s4(8) = mx
             massin_x4(8) = mx
          end if
       end if

       select case (ii)
       case(-1)

          if (imx==1) then 
             dsig = 0.0
             cycle
          end if

          do is = 1,4,1                                                                     ! t-channel squarks
             
             if ((ipart1<5).and.(ipart2<5)) then                                            ! NN case

                if ((is==2).or.(is==3)) then
                   if (iq==+1) cycle
                else if ((is==1).or.(is==4)) then
                   if (iq==-1) cycle
                end if
                il1 = is 
                il2 = il1 
                
                msu(-1:1:2) = msq(-is:is:2*is) 
                
             else if ((((ipart1>4).and.(ipart2>6)).or.((ipart1>6).and.(ipart2>4))).and.   &
                       (.not.((ipart1>6).and.(ipart2>6)))) then                             ! CC case
                
                if (is==2) then
                   if (iq==-1) cycle                                                        ! only t channel sd for u-ubar
                   il1 =  1                                                                 ! incoming quark tag
                else if (is==3) then
                   if (iq==-1) cycle                                                        ! only t channel ss for c-cbar
                   il1 =  4
                else if (is==1) then
                   if (iq==+1) cycle                                                        ! only u channel su for d-dbar
                   il1 =  2
                else if (is==4) then
                   if (iq==+1) cycle                                                        ! only u channel sc for s-sbar
                   il1 =  3
                end if
                il2 = il1 
                
                msu(-1:1:2) = msq(-is:is:2*is) 
                
             else if (((ipart1>=5).and.(ipart1<=6).and.(ipart2<=4)).or.                    &
                     ((ipart2>=5).and.(ipart2<=6).and.(ipart1<=4))     ) then               ! CN+ case
            
                if ( (is==1).or.(is==4) ) cycle                                             ! only sdown-type in t channel

                if (is==2) then                                                             ! u channel attached to il1
                   il1 = 1                                                                  ! quark of il1: u
                   il2 = 2                                                                  ! antiquark of il2; dbar
                   msu(-1:1:2) = msq(-1:1:2)
                else if (is==3) then
                   il1 = 4
                   il2 = 3
                   msu(-1:1:2) = msq(-4:4:8)
                end if
                
             else if (((ipart1>=7).and.(ipart1<=8).and.(ipart2<=4)).or.                    &
                      ((ipart2>=7).and.(ipart2<=8).and.(ipart1<=4))     ) then              ! CN- case
            
                if ( (is==2).or.(is==3) ) cycle                                             ! only sup-type in t channel
                                
                if (is==1) then                                                             ! u channel attached to il1
                   il1 = 2                                                                  ! quark of il1: d
                   il2 = 1                                                                  ! antiquark of il2: ubar
                   msu(-1:1:2) = msq(-2:2:4)
                else if (is==4) then
                   il1 = 3
                   il2 = 4
                   msu(-1:1:2) = msq(-3:3:6)
                end if
                
             end if

             if (i_ngtest==0) then 
                mst(-1:1:2) = msq(-is:is:2*is)                                                 ! this is the definiton of is
             else if (i_ngtest==1) then 
                mst(-1:1:2) = ms
                msu(-1:1:2) = ms
             else 
                print*, " IFCT_NN_X12: i_ngtest not set "
                call HARD_STOP
             end if

             if (icoll==0) then                                                             ! Tevatron   
                lumi_ex =   pdf1( il1)*pdf2( il2) + pdf1(-il2)*pdf2(-il1)                   ! il1=q; il2=qbar
             else                                                                           ! LHC
                lumi_ex =   pdf1( il1)*pdf2(-il2) + pdf1(-il2)*pdf2( il1)                   ! il1=q; il2=qbar 
             end if

             dsig = dsig +  lumi_ex                                                           &
                          * NN_QBB_NG(iq,iout,massin,Cs,Ct,Cl,Cr,mst,msu)/s**2       * t2_jac

          end do
       case(0,1)                                                                               ! born
          dsig = dsig +          LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * NN_QBB(iq,iout,massin,Cs,Ct,Cl,Cr)/s**2       * t2_jac
       case(2)                                                                                 ! virt
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * NN_QBV(iq,iout,massin,Cs,Ct,Cl,Cr,Cv)/s**2    * t2_jac
       case(3)                                                                                 ! matrix qb
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                            * NN_QBH(massin,mkraemer) /s**2                                   &
                            * s4 / (s4+m1**2) * t2_jac * s4_jac * c1_jac * ax_jac
       case(4)                                                                                 ! massfac qb
          dsig = dsig + nlo *    LUMI(inlo,20,icoll,idub,iq,x1,x2,mu)                         &
                          * ( NN_QBF1(iq,iout,massin,Csx,Ctx,Cl,Cr)/s**2/x * x_jac * t2x_jac  &
                            - NN_QBF2(iq,iout,massin,Cs,Ct,Cl,Cr)  /s**2   * x_jac * t2_jac  )
       case(5)                                                                                 ! massfac qg+gb
          dsig = dsig + nlo * (  LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                         &
                               + LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)  )                      &
                            * NN_QGF(iq,iout,massin,Csx,Ctx,Cl,Cr)/s**2/x  * x_jac * t2x_jac
       case(6)                                                                                 ! s4-matrix qg
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_QGH(massin,mkraemer)                                         &
                            - NN_QGHSUB(massin,mkraemer)) /s**2 * s4 / (s4+m1**2)             &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac
!          if ((theta_s4==0.0).and.(theta_s3==0.0)) then                                       ! moved to ii=7
!          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
!                            * NN_QGHSUB(massin,mkraemer)  /s**2 * s4 / (s4+m1**2)             &
!                                                        * t2_jac * s4_jac * c1_jac * ax_jac
!          end if 
       case(7)                                                                                 ! s4-matrix qg
          if ((theta_s4==0.0).and.(theta_s3==0.0)) then
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                            * NN_QGHSUB(massin,mkraemer)  /s**2 * s4 / (s4+m1**2)             &
                            * t2_jac * s4_jac * c1_jac * ax_jac
          else if ((theta_s4==1.0).and.(theta_s3==0.0)) then 
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_QGHSUB(massin,mkraemer)     /s**2 * s4   / (s4+m1**2)        &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac   &
                            - NN_QGHSUB(massin_s4,mkraemer)  /s**2 * s4s4 / (s4s4+m1**2)      &
                                            * prop_s4 * t2s4_jac * s4_jac * c1_jac * ax_jac )             
          else if ((theta_s4==1.0).and.(theta_s3==1.0)) then 
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_QGHSUB(massin_x4,mkraemer)  /s**2 * s4   / (s4+m1**2)        &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac   &
                            - NN_QGHSUB(massin_s4,mkraemer)  /s**2 * s4s4 / (s4s4+m1**2)      &
                                            * prop_s4 * t2s4_jac * s4_jac * c1_jac * ax_jac )             
          end if 
       case(8)                                                                                 ! s3-matrix qg
          if ((theta_s4==0.0).and.(theta_s3==1.0)) then
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_QGHSUB(massin,mkraemer)  /s**2 * s4   / (s4+m1**2)           &
                                                    * t2_jac   * s4_jac   * s3_jac * ax_jac   &
                            - NN_QGHSUB(massin_s3,mkraemer)  /s**2 * s4s3 / (s4s3+m1**2)      &
                                          * prop_s3 * t2s3_jac * s4s3_jac * s3_jac * ax_jac )             
          else if ((theta_s4==1.0).and.(theta_s3==1.0)) then
          dsig = dsig + nlo *   LUMI(inlo,30,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_QGHSUB(massin_x3,mkraemer)  /s**2 * s4   / (s4+m1**2)        &
                                                    * t2_jac   * s4_jac   * s3_jac * ax_jac   &
                            - NN_QGHSUB(massin_s3,mkraemer)  /s**2 * s4s3 / (s4s3+m1**2)      &
                                          * prop_s3 * t2s3_jac * s4s3_jac * s3_jac * ax_jac )             
          end if
       case(9)                                                                                 ! s4-matrix gb
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_GBH(massin,mkraemer)                                         &
                            - NN_GBHSUB(massin,mkraemer)) /s**2 * s4 / (s4+m1**2)             &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac
!          if ((theta_s4==0.0).and.(theta_s3==0.0)) then                                       ! moved to ii=7
!          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
!                          * NN_GBHSUB(massin,mkraemer)  /s**2 * s4 / (s4+m1**2)               &
!                                                        * t2_jac * s4_jac * c1_jac * ax_jac
!          end if
       case(10)                                                                                ! s4-matrix gb
          if ((theta_s4==0.0).and.(theta_s3==0.0)) then
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
                          * NN_GBHSUB(massin,mkraemer)  /s**2 * s4 / (s4+m1**2)               &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac
          else if ((theta_s4==1.0).and.(theta_s3==0.0)) then
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                           &
                          * ( NN_GBHSUB(massin,mkraemer)     /s**2 * s4   / (s4+m1**2)         &
                                                           * t2_jac * s4_jac * c1_jac * ax_jac &
                            - NN_GBHSUB(massin_s4,mkraemer)  /s**2 * s4s4 / (s4s4+m1**2)       &
                                                * prop_s4 * t2s4_jac * s4_jac * c1_jac * ax_jac )             
          else if ((theta_s4==1.0).and.(theta_s3==1.0)) then 
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_GBHSUB(massin_x4,mkraemer)  /s**2 * s4   / (s4+m1**2)        &
                                                        * t2_jac * s4_jac * c1_jac * ax_jac   &
                            - NN_GBHSUB(massin_s4,mkraemer)  /s**2 * s4s4 / (s4s4+m1**2)      &
                                            * prop_s4 * t2s4_jac * s4_jac * c1_jac * ax_jac )             
          end if
       case(11)                                                                                 ! matrix gb
          if ((theta_s4==0.0).and.(theta_s3==1.0)) then
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
                          * ( NN_GBHSUB(massin,mkraemer)  /s**2 * s4   / (s4+m1**2)           &
                                                   * t2_jac   * s4_jac   * s3_jac * ax_jac    &
                            - NN_GBHSUB(massin_s3,mkraemer)  /s**2 * s4s3 / (s4s3+m1**2)      &
                                          * prop_s3 * t2s3_jac * s4s3_jac * s3_jac * ax_jac )             
          else if ((theta_s4==1.0).and.(theta_s3==1.0)) then
          dsig = dsig + nlo *   LUMI(inlo,40,icoll,idub,iq,x1,x2,mu)                          &
                         * ( NN_GBHSUB(massin_x3,mkraemer)  /s**2 * s4   / (s4+m1**2)         &
                                                    * t2_jac   * s4_jac   * s3_jac * ax_jac   &
                           - NN_GBHSUB(massin_s3,mkraemer)  /s**2 * s4s3 / (s4s3+m1**2)       &
                                          * prop_s3 * t2s3_jac * s4s3_jac * s3_jac * ax_jac )             
          end if
       case default
          dsig = 0.0 
       end select
    end do

    dsig = dsig * symmfac                                   ! phase space symmetry factor 
    
    if ( (ii==3).or.(ii>=6) ) then                          ! phase space without 1/s**2 
       dsig = dsig / ( 2.0 * 256.0 * pi**4 )
    else
       dsig = dsig / ( 16.0 * pi )
    end if
    
    if (iscaling==0) dsig = dsig * x1_jac * x2_jac

    if (iscaling==0) dsig = dsig * 4.0/(abs(m1)+abs(m2))**2 * alpha**2 

    if (iscaling==0) dsig = dsig * gevpb

    ii_done(ii) = 1

  end function IFCT_NN_X12

! ------------------------------
! for type 0 : input m1,m2,s,t2,s4,c1,a2, cm system (A)
! for type 1 : input m1,m2,s,t2,s4,a2,s3, cm system (C)
!
  subroutine KIN_NN(type,m1,m2,s,t2,s4,c1_in,ax,s3,s5,u2,t1,u1,tp,up)
    integer,           intent(in)    :: type
    real(kind=double), intent(in)    :: m1,m2,s,t2,s4,c1_in,ax
    real(kind=double), intent(out)   :: s5,u2,t1,u1,tp,up
    real(kind=double), intent(inout) :: s3 
    real(kind=double)                 :: c1,s1,c2,norm,w1,w2,w3,e1,e2,p,cx,sx
    real(kind=double), dimension(0:5) :: test

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
       c1 = ( 2.D0*w3*e2 - s3 )/(2.D0*p*w3)
    else 
       print*, " KIN_NN: wrong input on type:",type
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

    if (test(0)>1.e-6) print *," KIN_NN: test0 better be zero "
    if (test(1)>1.e-6) print *," KIN_NN: test1 better be zero "
    if (test(2)>1.e-6) print *," KIN_NN: test2 better be zero "
    if (test(3)>1.e-6) print *," KIN_NN: test3 better be zero "
    if (test(4)>1.e-6) print *," KIN_NN: test4 better be zero "
    if (test(5)>1.e-6) print *," KIN_NN: test5 better be zero "

  end subroutine KIN_NN

! ------------------------------
! all general couplings : C(1) located at ( k1(q_in) , ipart1 )
!                         C(2) located at ( k2(q_out), ipart1 )
!                         C(3) located at ( k2(q_in) , ipart2 )
!                         C(4) located at ( k2(q_out), ipart2 )  for all Cl,Cr
!  iout = 1,2,3,4 for NN,CC,NC+,NC-
!  mx   = mw,mz s channel gauge boson mass 
! ------------------------------
  subroutine COUPLING_NN(s,m1,m2,iq,Cs,Ct,Cl,Cr,Cv,mx,iout,mkraemer)
    integer,                              intent(in)   :: iq
    real(kind=double),                    intent(in)   :: s,m1,m2
    integer,                              intent(out)  :: iout
    real(kind=double),                    intent(out)  :: mx
    real(kind=double), dimension(4),      intent(out)  :: Cs
    real(kind=double), dimension(99),     intent(out)  :: mkraemer
    complex(kind=double), dimension(4),   intent(out)  :: Ct,Cl,Cr,Cv

    integer  ::  ic1,ic2,in,ic,i1
    real(kind=double) :: t3,qq,sw2,cw2,sw,cw,v1,v2,a2,v1w,delta_ij,only_d,only_u 
    complex(kind=double) :: i,zzc(4,4),v2w,a2w,Clx(4),Crx(4),Cvx(4),Ctx(4)

    t3(iq)  = dble(iq) /2.0                                                  ! quark quantum numbers 
    qq(iq)  = 2.0/3.0 + ( dble(iq) - 1.0 ) /2.0

    i = (0.0,1.0) 

    sw = sqrt( 1.0 - mw**2/mz**2 )
    sw2 = sw**2
    cw2 = 1.0 - sw2
    cw  = sqrt(cw2)

    zzc(1:4,1:4) =  conjg( zz(1:4,1:4) )                                     ! complex conjugation of bw=zz

    if ((ipart1<5).and.(ipart2<5)) then                                            ! NN case 

       iout = 1 
       mx = mz

       v1  = - qq(iq)                                                        ! pqq coupling 
       v2  = - ( t3(iq) - 2.0 * sw2 * qq(iq)   )/(2.0*sw*cw)                 ! zqq couplings 
       a2  = -   t3(iq)                         /(2.0*sw*cw) 

       v1w = 0.0                                                             ! nnp coupling
       v2w = -i*aimag(zz(ipart1,3)*zzc(ipart2,3)-zz(ipart1,4)*zzc(ipart2,4))/(2.0*sw*cw) ! nnz couplings [higgsinos]
       a2w = -   real(zz(ipart1,3)*zzc(ipart2,3)-zz(ipart1,4)*zzc(ipart2,4))/(2.0*sw*cw)  

       Clx(1) = (zz(ipart1,1)*sw*(qq(iq)-t3(iq))+zz(ipart1,2)*cw*t3(iq))/(sqrt(2.0)*sw*cw) ! nqs coupings [gauginos]
       Clx(2) = (zz(ipart2,1)*sw*(qq(iq)-t3(iq))+zz(ipart2,2)*cw*t3(iq))/(sqrt(2.0)*sw*cw) 
       Crx(1) =  zzc(ipart1,1)*sw*qq(iq)                                /(sqrt(2.0)*sw*cw) 
       Crx(2) =  zzc(ipart2,1)*sw*qq(iq)                                /(sqrt(2.0)*sw*cw) 
       Clx(3:4) = Clx(1:2)
       Crx(3:4) = Crx(1:2) 
         
       Cvx(1:4) = conjg(Clx(1:4))/Clx(1:4)                                    ! set {Clot,Cupt,Cupu,Clou}=+-1 [identical for L and R]
         
    else if ( ((ipart1>4).and.(ipart2>6)).and.(.not.(ipart1>6)) ) then! CC case

       ic1 = ipart1 - 4                                                         ! ic1,ic2 = {1,2} in any case 
       ic2 = ipart2 - 4 
       if (ic1>2) ic1 = ic1 - 2 
       if (ic2>2) ic2 = ic2 - 2 

       delta_ij = 0.0 
       if (ic1.eq.ic2) delta_ij = 1.0

       iout = 2 
       mx = mz

       v1  = - qq(iq)                                                        ! pqq coupling
       v2  = - ( t3(iq) - 2.0 * sw2 * qq(iq)   )/(2.0*sw*cw)                 ! zqq couplings 
       a2  = -   t3(iq)                         /(2.0*sw*cw)

       v1w = - delta_ij                                                                    ! ccp coupling
       v2w = (2.0*delta_ij*(sw2-cw2)-uu(ic1,1)*uu(ic2,1)-vv(ic1,1)*vv(ic2,1))/(4.0*sw*cw)  ! ccz couplings 
       a2w = (                    uu(ic1,1)*uu(ic2,1)-vv(ic1,1)*vv(ic2,1))/(4.0*sw*cw)

       Clx(1) = uu(ic1,1)/(2.0*sw)                      ! cqs coupings [gauginos]
       Clx(2) = uu(ic2,1)/(2.0*sw)                      ! only u quarks in fermion number conserving vertices
       Clx(3) = vv(ic1,1)/(2.0*sw)
       Clx(4) = vv(ic2,1)/(2.0*sw)                      ! only d quarks in fermion number violating vertices 
       Crx(1:4) = (0.0,0.0)                             ! only left handed su(2) doublets 

       Cvx(1) = vv(ic1,1) / uu(ic1,1)                   ! set {Clot,Cupt,Cupu,Clou}
       Cvx(2) = vv(ic2,1) / uu(ic2,1)
       Cvx(3) = uu(ic1,1) / vv(ic1,1)
       Cvx(4) = uu(ic2,1) / vv(ic2,1)

    else if ((ipart2>4).and.(ipart1<=4)) then
       in = min(ipart1,ipart2)
       ic = max(ipart1,ipart2) - 4 
       if (ic>2) ic = ic - 2 

       mx = mw 

       v1  = 0.0                                                                ! pqq coupling
       v2  = - 1.0 /(2.0*sqrt(2.0)*sw)                                          ! wqq coupling
       a2  = - 1.0 /(2.0*sqrt(2.0)*sw) 

       v1w = 0.0                                                                ! pcn coupling 
       v2w = (  sqrt(2.0)*( zz(in,2)*uu(ic,1)+zzc(in,2)*vv(ic,1) )           & 
              + zz(in,3)*uu(ic,2) - zzc(in,4)*vv(ic,2) )/(2.0*sqrt(2.0)*sw)     ! wcn couplings
       a2w = (  sqrt(2.0)*( zzc(in,2)*vv(ic,1)-zz(in,2)*uu(ic,1) )           &
              - zz(in,3)*uu(ic,2) - zzc(in,4)*vv(ic,2) )/(2.0*sqrt(2.0)*sw)
         
!     momentum assignment attached to incoming q:    k1-p1-C(1) with C(1) the chargino(p1) is t channel
!                                              q:    k2-p1-C(3) with C(3) the chargino(p1) is u channel
!                                              qbar: k2-p2-C(2) with C(2) the neutralino(p2) is t channel
!                                              qbar: k1-p2-C(4) with C(4) the neutralino(p2) is u channel
!               nqs/cqs coupings for NC+ (u-dbar, iq=+1):
!                t channel with C(1)-chargino-u-sd;    C(2)-neutralino-dbar-sd; s-down
!                u channel with C(3)-chargino-dbar-su; C(4)-neutralino-u-sd;    s-up  
!               nqs/cqs coupings for NC- (d-ubar, iq=-1):
!                t channel with C(1)-chargino-d-su;    C(2)-neutralino-ubar-su; s-up
!                u channel with C(3)-chargino-ubar-sd; C(4)-neutralino-d-sd;    s-down
!               quantum numbers here adjusted for NC+ production: neutralino in t(u) channel coupled to d(u) quark
       Clx(1) = uu(ic,1)/(2.0*sw)                                                          ! t-channel chargino coupling
       Clx(2) = ( zz(in,1)*sw*( qq(-1)-t3(-1) ) + zz(in,2)*cw*t3(-1) )/(sqrt(2.0)*sw*cw)   ! t-channel neutralino coupling
       Clx(3) = vv(ic,1)/(2.0*sw)                                                          ! u-channel chargino coupling
       Clx(4) = ( zz(in,1)*sw*( qq(+1)-t3(+1) ) + zz(in,2)*cw*t3(+1) )/(sqrt(2.0)*sw*cw)   ! u-channel neutralino coupling

       Crx(1) = (0.0,0.0)
       Crx(2) =  zzc(in,1)*sw* qq(-1)/(sqrt(2.0)*sw*cw) 
       Crx(3) = (0.0,0.0)
       Crx(4) =  zzc(in,1)*sw* qq(+1)/(sqrt(2.0)*sw*cw) 

       Cvx(1) = vv(ic,1) / uu(ic,1)                                            ! set {Clot,Cupt,Cupu,Clou}
       Cvx(2) = conjg(Clx(2)) / Clx(2)                                         ! this is only required for degenerate squarks
       Cvx(3) = uu(ic,1) / vv(ic,1)
       Cvx(4) = conjg(Clx(4)) / Clx(4)        

       if ((ipart1<7).and.(ipart2<7)) then                                     ! NC+
          iout = 3 
       else if ((ipart1>6).or.(ipart2>6)) then                                 ! NC-
          iout = 4 
       end if

    else
       print*, " COUPLING_NN: something wrong ",ipart1,ipart2
    end if                                                                     ! end of NN,CC,Cn if construct

    if ( (abs(m1)+abs(m2)) < mz ) then
       print*, " COUPLINGS_NN: masses low, W/Z decays might be more suitable, decouple W/Z "
       v2   = 0.0
       a2   = 0.0
       v2w  = 0.0
       a2w  = 0.0
    end if

    Cs(1) = (v1*v1w)**2 + s * 2.0 * real( v1*v1w*v2*v2w )/(s-mx**2)         &  ! wim's conventions 
           + s**2*( v2**2+a2**2 )*( abs(v2w)**2+abs(a2w)**2 ) /(s-mx**2)**2    ! [v1,v2,v1w,a2 real]
    Cs(2) = (v1*v1w)**2 + s * 2.0 * real( v1*v1w*v2*v2w )/(s-mx**2)         &
           + s**2*( v2**2+a2**2 )*( abs(v2w)**2-abs(a2w)**2 ) /(s-mx**2)**2
    Cs(3) = ( 2.0 * real( v1*v1w*a2*a2w ) + s * 2.0*a2*v2                   &
           * 2.0*real( a2w*conjg(v2w) ) /(s-mx**2)         )/2.0
    Cs(4) = 0.0 

    Ctx(1) = v1*v1w + s*( v2+a2 )*( v2w-a2w ) /(s-mx**2) 
    Ctx(2) = v1*v1w + s*( v2+a2 )*( v2w+a2w ) /(s-mx**2) 
    Ctx(3) = v1*v1w + s*( v2-a2 )*( v2w+a2w ) /(s-mx**2) 
    Ctx(4) = v1*v1w + s*( v2-a2 )*( v2w-a2w ) /(s-mx**2) 

    if (iout==3) then                                                          ! NC+
       Cl(1:4) = Clx(1:4)
       Cr(1:4) = Crx(1:4)
       Ct(1:4) = Ctx(1:4)
       Cv(1:2) = Cvx(2:1:-1)
       Cv(3:4) = Cvx(4:3:-1) 
    else if (iout==4) then                                                     ! NC-
       Cl(1:3:2) = Clx(3:1:-2)
       Cl(2:4:2) = Clx(4:2:-2)
       Cr(1:3:2) = Crx(3:1:-2)
       Cr(2:4:2) = Crx(4:2:-2)
       Ct(1:2:1) =-Ctx(2:1:-1)
       Ct(3:4:1) =-Ctx(4:3:-1)
       Cv(1:2) = Cvx(2:1:-1)
       Cv(3:4) = Cvx(4:3:-1) 
    else 
       Cl(1:4) = Clx(1:4)
       Cr(1:4) = Crx(1:4)
       Ct(1:4) = Ctx(1:4)
       Cv(1:4) = Cvx(1:4)
    end if

    mkraemer(11) =  mx                                                       ! without mkraemer(1:4)
    mkraemer(12) =  v1
    mkraemer(13) =  v2
    mkraemer(14) = -a2
    mkraemer(15) =  v1w
    mkraemer(16) =  real(v2w)
    mkraemer(17) = -real(a2w)
    mkraemer(21) = -real(Clx(1))   
    mkraemer(22) = -real(Clx(2))  
    mkraemer(23) =  real(Crx(1))   
    mkraemer(24) =  real(Crx(2))   
    mkraemer(25) = -real(Clx(3))   
    mkraemer(26) = -real(Clx(4))  
    mkraemer(27) =  real(Crx(3))
    mkraemer(28) =  real(Crx(4))
    
    do i1=1,4
       if ( abs(aimag(Cl(i1))) > 1.e-8 ) then
          print*, " COUPLING_NN: problem with complex couplings ",i1,Cl(i1)
          call HARD_STOP
       end if
       if ( abs(aimag(Cr(i1))) > 1.e-8 ) then
          print*, " COUPLING_NN: problem with complex couplings ",i1,Cr(i1)
          call HARD_STOP
       end if
    end do

    only_u = ( 1.0 + iq )/2.0                                                ! incoming quarks
    only_d = ( 1.0 - iq )/2.0
    
    select case (iout)
    case(2)                                                        ! CC : u,d type quarks in t,u channel 
       mkraemer(21) = only_u * mkraemer(21)                        !      Ctl1
       mkraemer(22) = only_u * mkraemer(22)                        !      Ctl2
       mkraemer(25) = only_d * mkraemer(25)                        !      Ctr1
       mkraemer(26) = only_d * mkraemer(26)                        !      Ctr2
    case(3)                                                        ! CN+ : only u,dbar incoming state 
       mkraemer(12:28) = only_u * mkraemer(12:28)
    case(4)                                                        ! CN- : only ubar,d incoming state 
       mkraemer(12:28) = only_d * mkraemer(12:28)
    end select
    
    mkraemer(12:28) = 2.0 * sqrt(pi) * mkraemer(12:28)             ! rescaling with 2 and pi 
    
    Cs(1:3) = 16.0 * pi**2    * Cs(1:3)                            ! rescaling g^2 -> alpha
    Ct(1:4) =  4.0 * pi       * Ct(1:4)
    Cl(1:4) =  2.0 * sqrt(pi) * Cl(1:4)
    Cr(1:4) =  2.0 * sqrt(pi) * Cr(1:4)
    
  end subroutine COUPLING_NN

end module xx_integral_nn

