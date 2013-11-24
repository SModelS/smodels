! ===========================================================================================================
module xx_prospino_subroutine
  use xx_kinds
  use xx_public_variables
  use xx_pass_integ
  use xx_integral_ng
  use xx_integral_ns
  use xx_integral_nn
  use xx_integral_ll
  use xx_integral_tb
  use xx_integral_lq
  use xx_integral_le
  use xx_integral_hh
  use xx_integral_ht
  use xx_initialize
  use xx_in_out
  implicit none 
  private :: GET_UNITS_READ_WRITE
  public :: PROSPINO, PROSPINO_OPEN_CLOSE
contains
! ------------------------------
  subroutine PROSPINO(inlo,isq_ng_in,icoll_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in)
    
    implicit none

    integer,                  intent(in) :: inlo,isq_ng_in,icoll_in,i_error_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in
    character(len=2),         intent(in) :: final_state_in
  
    integer                              :: n0,n1,ns,ii_min,ii_max,ns_min,ns_max
    integer                              :: nin1,nin2,ndat1,ndat2,ndat3
    integer                              :: isquark1_min,isquark1_max
    integer                              :: isquark2_min,isquark2_max
    integer                              :: i_lofast,ifast,inext_sq
    integer,dimension(4)                 :: ivegas
    real(kind=double)                    :: run0,run1,kfac,ms1_print,ms2_print,kng
    real(kind=double), dimension(-1:22)  :: fct,rel
    real(kind=double), dimension(1:20)   :: unimass
    real(kind=double), dimension(0:99)   :: lowmass
    character(len=2)                     :: final_state_new
    logical                              :: lvalid_global,lvalid_iteration

    integer           :: ilo_common, inlo_common, ionlylo_common, idg_common, ing_common
    common/CONST5/       ilo_common, inlo_common, ionlylo_common, idg_common, ing_common        ! to fill idg_common and ing_common

    integer           :: i_lofast_common                                                        ! only communicating with initpdf.f
    common/I_LOFAST/     i_lofast_common

    fct(-1:22) = 0.0                                                                            ! initialize cross section and rel. error
    rel(-1:22) = 0.0 
    kfac       = 0.0
    kng        = 0.0
          
    isquark1_min = 0                                                                            ! initialize squark loop
    isquark1_max = 0
    isquark2_min = 0
    isquark2_max = 0
    inext_sq     = 0                                                                            ! initialize squark-sum indicator in loop

    icoll           = icoll_in                                                                  ! move input into public variables
    isq_ng          = isq_ng_in
    final_state     = final_state_in
    final_state_new = final_state_in                                                            ! copy for second select statement
!----------------------------------------------------------------------------------------
    i_ngtest = 0     ! test degenerate squark masses : default[0] , test[1]             !
!----------------------------------------------------------------------------------------
    if (i_ngtest == 1 ) print*, " PROSPINO: i_ngtest in testing mode ",i_ngtest
    
!----------------------------------------------------------------------------------------
    i_lofast = 0     ! run fast sq-gl Born terms : default[0] , test[1]                 !
!----------------------------------------------------------------------------------------
    if (i_lofast == 1 ) then 
       print*, " PROSPINO: i_lofast in testing mode, sorry ",i_lofast
       call HARD_STOP                                                                           ! mode only for testing
    end if
    i_lofast_common = i_lofast
    
    if ( (abs(isquark1_in) > 4).or.(abs(isquark2_in) > 4) ) then                                ! protect against    
       print*, " PROSPINO: for heavy-flavor squarks please use the final_state flag "
       call HARD_STOP
    end if
            
    if (inlo == 0) then 
       ii_max = 0
    else if (inlo == 1 ) then 
       ii_max = 11
    else
       print*, " PROSPINO: inlo not set correctly?, continue in LO mode "
       ii_max = 0
    end if

    if (i_error_in == 0) then
       ns_min = 0
       ns_max = 0 
    else if (i_error_in == 1) then
       ns_min = -10
       ns_max = +10
    else 
       print*, " PROSPINO: i_error_in not set correctly?, continue with central scale "
       ns_min = 0
       ns_max = 0 
    end if
    
    if (icoll == 0) then                                                                          ! set the collider energy Tevatron/LHC
       sc = 1960.0**2 
    else if (icoll == 1) then 
       sc = 14000.0**2 
    else if (icoll == 2) then 
       sc =  7000.0**2 
    else if (icoll == 3) then 
       sc =  8000.0**2 
    else 
       print*, " PROSPINO: icoll not set correctly?, continue LHC "
       sc = 14000.0**2 
    end if

!----------------------------------------------------------------------------------------
    ifla_le = 1       ! specify the incoming quark flavor for LQ+lepton                 !
                      !  d[2], s[3], b[5], u[1], c[4]                                   !
!----------------------------------------------------------------------------------------
    if (final_state == 'le') then
       print*, " LQ+lepton final state "
       print*, "  R violating coupling and mixing normalized to unity "  
       print*, "  incoming quark set in Xprospino_subroutine.f90 = ",ifla_le
    end if

    call GET_UNITS_READ_WRITE(nin1,nin2,ndat1,ndat2,ndat3)                                      ! obtain units from where the OPEN_CLOSE routine got them

    call INIT_GLOBAL(ipart1_in,ipart2_in,lvalid_global)                                         ! output lvalid_global, possibly skip s
    if (.not. lvalid_global) then 
       print*, " MAIN: problem inside INIT_GLOBAL, cycle "
       return
    end if

    do n1   = 0,0,200                                                                           ! can be used to shift MSSM masses
    do n0   = 0,0,150                                                                           ! in routine INIT_SUSY (Xinitialize.f90)
    do ns   = ns_min,ns_max,10                                                                  ! change both scales simultaneously

       run0 = n0
       run1 = n1
       scafac = 2.0**(ns/10.0)
       
       call INIT_SUSY(nin1,nin2,ipart1_in,ipart2_in,isq_ng_in,run0,run1,unimass,lowmass)        ! get spectrum and initialize mass arrays.
       call DAT3(ndat3,lowmass,unimass)                                                         ! print out variables and model parameters
    
       ms1_print = ms                                                                           ! default mode
       ms2_print = ms
       idg_common  = 1                                                                          ! make this the sq-gl default mode 
       if (isq_ng == 0) then                                                                    ! only mass degenerate case
          ing_common = 0                                                                        ! for sq-gl only
          ii_min = 0                                                                            ! for new channels only
          isquark1_min =  9
          isquark1_max =  9
          isquark2_min =  9
          isquark2_max =  9
       else if (isq_ng == 1) then                                                               ! also different squark masses
          ing_common =  1                                                                       ! for sq-gl only
          ii_min = -1                                                                           ! for new channels only
          if ( (isquark1_in == 0).and.(isquark2_in == 0).and.          &                        ! summing mode for two external squarks
                (final_state == 'ss')                    ) then
             isquark1_min = -4                                                                  ! including sbottoms in the final state
             isquark1_max =  4
             isquark2_min = -4
             isquark2_max =  4
             ms1_print = ms
             ms2_print = ms
          else if ( (isquark1_in == 0).and.(isquark2_in == 0).and.     &                        ! summing mode for two external squarks
                (final_state == 'sb')                    ) then
             isquark1_min = -5                                                                  ! including sbottoms in the final state
             isquark1_max =  5
             isquark2_min = -5
             isquark2_max =  5
             ms1_print = ms
             ms2_print = ms
          else if ( (isquark1_in == 0).and.(isquark2_in == 0).and.     &                        ! one external squark to sum over 
                    ((final_state == 'sg').or.                         &
                     (final_state == 'ns')                  ) ) then
             isquark1_min = -4                                                                  ! no sbottoms in initial state
             isquark1_max =  4
             isquark2_min =  9
             isquark2_max =  9
             ms1_print = ms
          else if ((final_state == 'sb').or.                           &                        ! run individual channels if they exist
                   (final_state == 'ss').or.                           &
                   (final_state == 'sg').or.                           &
                   (final_state == 'ns')                 ) then
             isquark1_min = isquark1_in
             isquark1_max = isquark1_in
             isquark2_min = isquark2_in
             isquark2_max = isquark2_in
             ms1_print = msq(isquark1_min)
             ms2_print = msq(isquark2_min)
          else                                                                                  ! default, nothing to sum over
             isquark1_min =  9
             isquark1_max =  9
             isquark2_min =  9
             isquark2_max =  9
          end if
       else 
          print*, " PROSPINO: isq_ng not set correctly? ",isq_ng
          call HARD_STOP
       end if

       fct(20:22) = 0.0                                                                         ! initialize in case of loops 
       do isquark1 = isquark1_min,isquark1_max,1                                                ! note this is always a loop
       do isquark2 = isquark2_min,isquark2_max,1                                                !    so `cycle' is the way out
          if ( (isquark1 == 0).or.(isquark2 == 0) ) cycle                                       ! just a tag for summing in isquark_in
          if ( (isquark1 /= isquark1_min) .or. (isquark2 /= isquark2_min) )  inext_sq = 1 

          print*, " running with  final_state = ",final_state
          print*, "               ipart*_in   = ",ipart1_in,ipart2_in
          print*, "               isquark*    = ",isquark1,isquark2
             
          select case (final_state)                                                             ! setup for old prospino
          case ('gg')
             call INIT_IFAST(1,inext_sq,ifast)
             call FILL_COMMONS_SQ_GL(inlo,scafac) 
             call INTEGGG( ifast, fct(-1:-1), rel(-1:-1), fct(0:0), rel(0:0), fct(20:20), rel(20:20) )
          case ('sb')
             idg_common = 0
             if ((isquark1 == isquark1_max).and.(isquark2 == isquark2_max)) idg_common = 1
             call INIT_IFAST(1,inext_sq,ifast)
             call FILL_COMMONS_SQ_GL(inlo,scafac) 
             call INTEGSB( ifast, fct(-1:-1), rel(-1:-1), fct(0:0), rel(0:0), fct(20:20), rel(20:20) )
          case ('sg')
             idg_common = 0
             if (isquark1 == isquark1_max) idg_common = 1
             call INIT_IFAST(1,inext_sq,ifast)
             call FILL_COMMONS_SQ_GL(inlo,scafac) 
             call INTEGSG( ifast, fct(-1:-1), rel(-1:-1), fct(0:0), rel(0:0), fct(20:20), rel(20:20) )
          case ('ss')
             if (isquark2 < isquark1) cycle
             idg_common = 0
             if ((isquark1 == isquark1_max).and.(isquark2 == isquark2_max)) idg_common = 1
             call INIT_IFAST(1,inext_sq,ifast)
             call FILL_COMMONS_SQ_GL(inlo,scafac)
             call INTEGSS( ifast, fct(-1:-1), rel(-1:-1), fct(0:0), rel(0:0), fct(20:20), rel(20:20) )
          case ('xx')
             call INIT_IFAST(1,inext_sq,ifast)
             call FILL_COMMONS_SQ_GL(inlo,scafac)
             call INTEGST( ifast, fct(0:0) , rel(0:0) , fct(20:20) , rel(20:20) )
          case default                                                                          ! below the new channels

!----------------------------------------------------------------------------------------------
!                  ii = -1    born consistent with non-degenerate squarks                     !
!                       0,1   born consistent[0] , nlo-type[1]                                !
!                       2     virtual correction                                              !
!                       3     real gluon emission                                             !
!                       4     subtraction term, crossed channels...                           !
!                       20    complete nlo (only sum of 1...9)                                !
!                       21    summed born term for non-degenerate squarks                     !
!                       22    summed born term multipled with squark-degenerate K factor      !
!----------------------------------------------------------------------------------------------
             do ii=ii_min,ii_max,1
                if ( (isquark1/=isquark1_max).and.(ii/=-1) ) cycle                              ! only run complete mode once
                
                print*, "               ii          = ",ii
                call INIT_ITERATION(lvalid_iteration)                                           ! always needed (pdf initialization...)
                if (.not. lvalid_iteration) cycle
                
                call INIT_ALPHAS                                                                ! initialize alpha_s, check Xget_pdf.f
                call INIT_VEGAS(i_lofast,ivegas)
                call INIT_IFAST(0,inext_sq,ifast)

                select case (final_state_new)                                                   ! the actual vegas integration
                case('ng')
                   call INTEG(IFCT_NG_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('ns')
                   call INTEG(IFCT_NS_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('nn')
                   call INTEG(IFCT_NN_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('ll')
                   call INTEG(IFCT_LL_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('tb','bb')
                   call INTEG(IFCT_TB_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('lq')
                   call INTEG(IFCT_LQ_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('le')
                   call INTEG(IFCT_LE_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('hh')
                   call INTEG(IFCT_HH_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                case('ht')
                   call INTEG(IFCT_HT_X12,dim(ii),ivegas,ifast,fct(ii),rel(ii))
                end select                                                                      ! close final_state_new
                print*, " MAIN: check number attempts stopped early ",n_faulty                  ! non-fatal errors kept track of 
                
                if (abs(fct(ii))<1.e-20) then                                                   ! reject obvious garbage
                   fct(ii) = 0.0 
                   rel(ii) = 0.0
                else if (abs(rel(ii))<1.e-20) then 
                   rel(ii) = 0.0
                end if
             end do                                                                             ! close the ii loop

             if ((final_state=='hh').or.(final_state=='ht')) then                               ! compute K factor for 2HDM, #22 not needed
                fct(22) = fct(3)
                fct(3)  = 0.0
             end if

             fct(20) = sum( fct(1:11) )                                                         ! sum all nlo contributions
             if ( fct(20)>1.e-20 )                                                            &
                  rel(20) = abs( sqrt( sum( rel(1:9)**2*fct(1:9)**2 ) )/fct(20) )               ! relative error added in quadrature

          end select                                                                            ! close final_state

          fct(21) = fct(21) + fct(-1)
          rel(21) = rel(21) + rel(-1)**2*fct(-1)**2

       end do                                                                                   ! close the isquark1 loop
       end do                                                                                   ! close the isquark2 loop

       if ((abs(fct(21)))>1.D-12) rel(21) = abs( sqrt( rel(21) )/fct(21) )                      ! just like formula for rel(20)
       if ((abs(fct(0)))>1.D-12)     kfac = fct(20)/fct(0)                                      ! compute k factor, avoid core dump
       
       if ((final_state/='hh').and.(final_state/='ht')) then                                    ! means for all SUSY channels
          fct(22) = fct(21) * kfac                                                              ! scale new LO result by K factor
          if ((abs(fct(0)))>1.D-12)     kng  = fct(21)/fct(0)                                   ! compute factor between two born results
       else                                                                                     ! fct(22) set before for Higgs channels
          fct(22) = fct(22) + fct(20)
          if ((abs(fct(20)))>1.D-12)    kng  = fct(22)/fct(20)                                  ! compute factor between two NLO results
       end if

       call DAT1(ndat1,run0,run1,ms1_print,ms2_print,fct,rel,kfac,kng)                          ! short results output file 
       call DAT2(ndat2,run0,run1,ms1_print,ms2_print,fct,rel,kfac,kng)                          ! long results output file 
       
    end do                                                                                      ! close the ns loop
    end do                                                                                      ! close the n0 loop
    end do                                                                                      ! close the n1 loop

    return
    
  end subroutine PROSPINO

! ------------------------------
! open (status=0) or close (status=1) all input and output files 
  subroutine PROSPINO_OPEN_CLOSE(status)
    
    implicit none

    integer, intent(in) :: status
    integer             :: nin1,nin2,ndat1,ndat2,ndat3

    character(len=20), parameter :: form1="(/,a121)"

    call GET_UNITS_READ_WRITE(nin1,nin2,ndat1,ndat2,ndat3)

    if (status == 0 ) then 
       open(unit=nin2,  file="prospino.in.les_houches", action="read" )                         ! open all files here
       open(unit=ndat1, file="prospino.dat",            action="write")
       open(unit=ndat2, file="prospino.dat2",           action="write")
       open(unit=ndat3, file="prospino.dat3",           action="write")
    else if (status == 1) then 
       if ((final_state=='hh').or.(final_state=='ht')) then
       write(ndat1, fmt=form1)                                                                 &
       "    i1 i2  dummy0 dummy1 scafac  m1    m2      angle LO[pb]   rel_error NLO[pb]   rel_error   K    NLO_SUSY[pb]           "
       else 
       write(ndat1, fmt=form1)                                                                 &
       "    i1 i2  dummy0 dummy1 scafac  m1    m2      angle LO[pb]   rel_error NLO[pb]   rel_error   K    LO_ms[pb] NLO_ms[pb]   "
       end if

       close(nin1)                                                                              ! close all input/output files
       close(nin2)
       close(ndat1)
       close(ndat2)
       close(ndat3)
    else 
       print*, " something wrong with input/output status "
       call HARD_STOP
    end if

  end subroutine PROSPINO_OPEN_CLOSE

! ------------------------------
! set these at one point in the program
  subroutine GET_UNITS_READ_WRITE(nin1,nin2,ndat1,ndat2,ndat3)
    
    implicit none

    integer, intent(out) :: nin1,nin2,ndat1,ndat2,ndat3

    nin1  = 95                                                                                  ! set all the input/output unit numbers
    nin2  = 96
    ndat1 = 97
    ndat2 = 98
    ndat3 = 99
    
  end subroutine GET_UNITS_READ_WRITE

end module xx_prospino_subroutine


